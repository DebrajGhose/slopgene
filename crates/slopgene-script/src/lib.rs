//! SlopGene scripting: a Rhai engine in which sequence documents are ordinary values.
//!
//! Scripts open files, inspect and edit documents, run analyses and save results. The same
//! engine runs from the command line (`slopgene run script.rhai`) and inside the desktop app.
//! Positions are 1-based and inclusive, as everywhere else in SlopGene.
//!
//! ```rhai
//! for path in glob("plasmids/**/*.dna") {
//!     let doc = open(path);
//!     let cuts = doc.cut_counts("EcoRI,BamHI");
//!     print(`${doc.name}\t${doc.length} bp\tEcoRI x${cuts.EcoRI}`);
//! }
//! ```

mod api;
pub mod files;

pub use api::{API, ApiEntry, api_json, api_text};

use rhai::{Array, Dynamic, Engine, EvalAltResult, Map, Scope};
use slopgene_core::digest::{circularize, digest, ligate};
use slopgene_core::enzymes::{Enzyme, EnzymeDb, find_sites};
use slopgene_core::io::{self, Format};
use slopgene_core::primers::{self, BindingRules, TmConditions};
use slopgene_core::protocol::{self, Change, NewFile, Orientation, Protocol, Run};
use slopgene_core::render::{self, ApproxMeasure, MapOptions, Theme};
use slopgene_core::search::{find_dna, find_protein};
use slopgene_core::translate::{GeneticCode, OrfOptions, find_orfs, translate, translate_feature};
use slopgene_core::*;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

type R<T> = std::result::Result<T, Box<EvalAltResult>>;

/// A sequence document inside a script, with the file it came from (if any).
#[derive(Clone, Debug)]
pub struct Doc {
    pub doc: Document,
    pub path: Option<PathBuf>,
}

fn err<T>(message: impl std::fmt::Display) -> R<T> {
    Err(message.to_string().into())
}

fn core<T>(result: slopgene_core::Result<T>) -> R<T> {
    result.map_err(|e| e.to_string().into())
}

fn to_dyn(value: serde_json::Value) -> Dynamic {
    rhai::serde::to_dynamic(value).unwrap_or(Dynamic::UNIT)
}

fn usize_of(n: i64, what: &str) -> R<usize> {
    usize::try_from(n).map_err(|_| format!("{what} must not be negative (got {n})").into())
}

/// A 1-based inclusive range on the document.
fn span_of(doc: &Document, start: i64, end: i64) -> R<Span> {
    core(Span::from_1based(usize_of(start, "start")?, usize_of(end, "end")?, doc.len(), doc.is_circular()))
}

fn db() -> &'static EnzymeDb {
    EnzymeDb::builtin()
}

/// Enzymes from a set name or a comma-separated list.
fn enzymes_from(spec: &str) -> R<Vec<&'static Enzyme>> {
    core(db().set(spec))
}

fn feature_dyn(f: &Feature, seq_len: usize, index: usize) -> Dynamic {
    let mut v = json::feature_to_value(f, seq_len);
    v["index"] = serde_json::json!(index + 1);
    to_dyn(v)
}

/// Find a feature by exact name (case-insensitive) or `#N`.
fn find_feature<'a>(doc: &'a Document, name: &str) -> R<(usize, &'a Feature)> {
    if let Some(n) = name.strip_prefix('#').and_then(|n| n.parse::<usize>().ok()) {
        return match doc.features.get(n.wrapping_sub(1)) {
            Some(f) => Ok((n - 1, f)),
            None => err(format!("there is no feature #{n}")),
        };
    }
    let mut found = doc.features.iter().enumerate().filter(|(_, f)| f.name.eq_ignore_ascii_case(name));
    match (found.next(), found.next()) {
        (Some(hit), None) => Ok(hit),
        (None, _) => err(format!("no feature named '{name}'")),
        (Some(_), Some(_)) => err(format!("several features are named '{name}'; use \"#N\"")),
    }
}

fn clean(bases: &str) -> R<Vec<u8>> {
    alphabet::clean_nucleotides(bases).map(|b| b.to_ascii_uppercase()).map_err(|e| e.into())
}

/// Print one line to standard output. When the reader has gone away (for example `| head`),
/// exit quietly instead of failing, as command-line tools are expected to.
pub fn stdout_line(text: &str) {
    use std::io::Write;
    let mut out = std::io::stdout().lock();
    if writeln!(out, "{text}").and_then(|_| out.flush()).is_err() {
        std::process::exit(0);
    }
}

/// Options for running a script.
#[derive(Clone, Default)]
pub struct RunOptions {
    /// Values passed to the script as the `ARGS` map.
    pub args: Vec<(String, String)>,
    /// Where `print` output goes. `None` prints to standard output.
    pub output: Option<Arc<Mutex<String>>>,
    /// A document made available to the script as `DOC` (the desktop app passes the open file).
    pub document: Option<Doc>,
}

/// Build an engine with the whole SlopGene API registered.
pub fn engine(output: Option<Arc<Mutex<String>>>) -> Engine {
    let mut engine = Engine::new();
    engine.set_max_expr_depths(256, 256);
    if output.is_none() {
        engine.on_print(stdout_line);
    }
    if let Some(buffer) = output.clone() {
        engine.on_print(move |text| {
            if let Ok(mut b) = buffer.lock() {
                b.push_str(text);
                b.push('\n');
            }
        });
        let debug_buffer = output;
        engine.on_debug(move |text, _, _| {
            if let Some(b) = &debug_buffer
                && let Ok(mut b) = b.lock()
            {
                b.push_str(text);
                b.push('\n');
            }
        });
    }
    register_documents(&mut engine);
    register_analysis(&mut engine);
    register_utilities(&mut engine);
    register_protocols(&mut engine);
    files::register(&mut engine);
    engine
}

/// Run a script. Returns the value of its last expression.
pub fn run(source: &str, options: RunOptions) -> std::result::Result<Dynamic, String> {
    let engine = engine(options.output.clone());
    let mut scope = Scope::new();
    let mut args = Map::new();
    for (k, v) in &options.args {
        args.insert(k.as_str().into(), v.clone().into());
    }
    scope.push_constant("ARGS", args);
    if let Some(doc) = options.document {
        scope.push("DOC", doc);
    }
    engine.eval_with_scope::<Dynamic>(&mut scope, source).map_err(|e| e.to_string())
}

/// Run a script with `DOC` bound to a document, and return the script's result together with
/// `DOC` as the script left it (the desktop app applies those changes to the open file).
pub fn run_with_document(source: &str, options: RunOptions) -> (std::result::Result<Dynamic, String>, Option<Doc>) {
    let engine = engine(options.output.clone());
    let mut scope = Scope::new();
    let mut args = Map::new();
    for (k, v) in &options.args {
        args.insert(k.as_str().into(), v.clone().into());
    }
    scope.push_constant("ARGS", args);
    if let Some(doc) = options.document {
        scope.push("DOC", doc);
    }
    let result = engine.eval_with_scope::<Dynamic>(&mut scope, source).map_err(|e| e.to_string());
    let doc = scope.get_value::<Doc>("DOC");
    (result, doc)
}

/// A script value as JSON (documents become SlopGene document JSON).
pub fn to_json_value(value: &Dynamic) -> serde_json::Value {
    if let Some(doc) = value.clone().try_cast::<Doc>() {
        return json::document_to_value(&doc.doc, false);
    }
    rhai::serde::from_dynamic::<serde_json::Value>(value).unwrap_or_else(|_| serde_json::Value::String(value.to_string()))
}

fn register_documents(engine: &mut Engine) {
    engine.register_type_with_name::<Doc>("Doc");

    // Opening and creating.
    engine.register_fn("open", |path: &str| -> R<Doc> {
        let doc = core(io::read_one(path))?;
        Ok(Doc { doc, path: Some(PathBuf::from(path)) })
    });
    engine.register_fn("open_all", |path: &str| -> R<Array> {
        let docs = core(io::read_path(path))?;
        Ok(docs.into_iter().map(|doc| Dynamic::from(Doc { doc, path: Some(PathBuf::from(path)) })).collect())
    });
    engine.register_fn("new_dna", |name: &str, seq: &str| -> R<Doc> { Ok(Doc { doc: Document::new_dna(name, clean(seq)?, Topology::Linear), path: None }) });
    engine.register_fn("new_dna", |name: &str, seq: &str, topology: &str| -> R<Doc> {
        let t = Topology::parse(topology).ok_or("topology must be linear or circular")?;
        Ok(Doc { doc: Document::new_dna(name, clean(seq)?, t), path: None })
    });
    engine.register_fn("from_json", |text: &str| -> R<Doc> {
        let value: serde_json::Value = serde_json::from_str(text).map_err(|e| e.to_string())?;
        Ok(Doc { doc: core(json::document_from_value(&value))?, path: None })
    });

    // Properties.
    engine.register_get("name", |d: &mut Doc| d.doc.name.clone());
    engine.register_set("name", |d: &mut Doc, v: &str| d.doc.name = v.to_string());
    engine.register_get("description", |d: &mut Doc| d.doc.meta.description.clone());
    engine.register_set("description", |d: &mut Doc, v: &str| d.doc.meta.description = v.to_string());
    engine.register_get("length", |d: &mut Doc| d.doc.len() as i64);
    engine.register_fn("len", |d: &mut Doc| d.doc.len() as i64);
    engine.register_get("seq", |d: &mut Doc| d.doc.seq_str());
    engine.register_get("topology", |d: &mut Doc| d.doc.topology.name().to_string());
    engine.register_get("molecule", |d: &mut Doc| d.doc.molecule.name().to_string());
    engine.register_get("is_circular", |d: &mut Doc| d.doc.is_circular());
    engine.register_get("gc", |d: &mut Doc| alphabet::gc_fraction(&d.doc.seq) * 100.0);
    engine.register_get("path", |d: &mut Doc| d.path.as_ref().map(|p| p.to_string_lossy().replace('\\', "/")).unwrap_or_default());
    engine.register_get("features", |d: &mut Doc| -> Array {
        let n = d.doc.len();
        d.doc.features.iter().enumerate().map(|(ii, f)| feature_dyn(f, n, ii)).collect()
    });
    engine.register_get("primers", |d: &mut Doc| -> Array {
        let n = d.doc.len();
        d.doc.primers.iter().map(|p| to_dyn(json::primer_to_value(p, n))).collect()
    });
    engine.register_fn("to_string", |d: &mut Doc| format!("Doc({}, {} bp, {})", d.doc.name, d.doc.len(), d.doc.topology.name()));
    engine.register_fn("to_debug", |d: &mut Doc| format!("Doc({}, {} bp, {})", d.doc.name, d.doc.len(), d.doc.topology.name()));
    engine.register_fn("to_map", |d: &mut Doc| to_dyn(json::document_to_value(&d.doc, true)));
    engine.register_fn("to_json", |d: &mut Doc| serde_json::to_string_pretty(&json::document_to_value(&d.doc, true)).unwrap_or_default());
    engine.register_fn("copy", |d: &mut Doc| d.clone());

    // Reading parts.
    engine.register_fn("subseq", |d: &mut Doc, start: i64, end: i64| -> R<String> {
        Ok(String::from_utf8_lossy(&d.doc.subseq(span_of(&d.doc, start, end)?)).into_owned())
    });
    engine.register_fn("feature", |d: &mut Doc, name: &str| -> R<Dynamic> {
        let (index, f) = find_feature(&d.doc, name)?;
        Ok(feature_dyn(f, d.doc.len(), index))
    });
    engine.register_fn("has_feature", |d: &mut Doc, name: &str| d.doc.features.iter().any(|f| f.name.eq_ignore_ascii_case(name)));
    engine.register_fn("features_of_type", |d: &mut Doc, kind: &str| -> Array {
        let n = d.doc.len();
        d.doc.features.iter().enumerate().filter(|(_, f)| f.kind.eq_ignore_ascii_case(kind)).map(|(ii, f)| feature_dyn(f, n, ii)).collect()
    });
    engine.register_fn("feature_seq", |d: &mut Doc, name: &str| -> R<String> {
        let (_, f) = find_feature(&d.doc, name)?;
        Ok(String::from_utf8_lossy(&f.extract(&d.doc.seq)).into_owned())
    });
    engine.register_fn("translate_feature", |d: &mut Doc, name: &str| -> R<String> {
        let (_, f) = find_feature(&d.doc, name)?;
        Ok(String::from_utf8_lossy(&translate_feature(f, &d.doc)).into_owned())
    });
    engine.register_fn("extract", |d: &mut Doc, start: i64, end: i64| -> R<Doc> {
        let span = span_of(&d.doc, start, end)?;
        Ok(Doc { doc: core(edit::extract(&d.doc, span))?, path: None })
    });

    // Editing in place. Primer binding sites are recomputed after sequence changes.
    let rules = BindingRules::default();
    engine.register_fn("insert", move |d: &mut Doc, pos: i64, bases: &str| -> R<()> {
        let pos = usize_of(pos, "position")?;
        if pos == 0 || pos > d.doc.len() + 1 {
            return err(format!("insert position must be 1..={} (the new bases start there)", d.doc.len() + 1));
        }
        core(edit::insert(&mut d.doc, pos - 1, &clean(bases)?))?;
        primers::refresh_binding_sites(&mut d.doc, &rules);
        Ok(())
    });
    engine.register_fn("delete", move |d: &mut Doc, start: i64, end: i64| -> R<()> {
        let span = span_of(&d.doc, start, end)?;
        core(edit::delete(&mut d.doc, span))?;
        primers::refresh_binding_sites(&mut d.doc, &rules);
        Ok(())
    });
    engine.register_fn("replace", move |d: &mut Doc, start: i64, end: i64, bases: &str| -> R<()> {
        let span = span_of(&d.doc, start, end)?;
        core(edit::replace(&mut d.doc, span, &clean(bases)?))?;
        primers::refresh_binding_sites(&mut d.doc, &rules);
        Ok(())
    });
    engine.register_fn("reverse_complement", |d: &mut Doc| edit::reverse_complement(&mut d.doc));
    engine.register_fn("invert", move |d: &mut Doc, start: i64, end: i64| -> R<()> {
        let span = span_of(&d.doc, start, end)?;
        core(edit::invert_region(&mut d.doc, span))?;
        primers::refresh_binding_sites(&mut d.doc, &rules);
        Ok(())
    });
    engine.register_fn("set_origin", |d: &mut Doc, pos: i64| -> R<()> {
        let pos = usize_of(pos, "position")?;
        if pos == 0 {
            return err("positions start at 1");
        }
        core(edit::set_origin(&mut d.doc, pos - 1))
    });
    engine.register_set("topology", |d: &mut Doc, t: &str| -> R<()> {
        let t = Topology::parse(t).ok_or("topology must be linear or circular")?;
        edit::set_topology(&mut d.doc, t);
        Ok(())
    });
    engine.register_fn("add_feature", |d: &mut Doc, name: &str, kind: &str, start: i64, end: i64, strand: &str| -> R<Dynamic> {
        let span = span_of(&d.doc, start, end)?;
        let strand = Strand::parse(strand).ok_or("strand must be forward, reverse, both or none")?;
        let feature = Feature::new(name, kind, span, strand);
        d.doc.features.push(feature);
        let index = d.doc.features.len() - 1;
        Ok(feature_dyn(&d.doc.features[index], d.doc.len(), index))
    });
    engine.register_fn("add_feature", |d: &mut Doc, spec: Map| -> R<Dynamic> {
        let get = |k: &str| spec.get(k).map(|v| v.to_string());
        let name = get("name").ok_or("add_feature needs a name")?;
        let kind = get("type").unwrap_or_else(|| "misc_feature".into());
        let start = spec.get("start").and_then(|v| v.as_int().ok()).ok_or("add_feature needs start")?;
        let end = spec.get("end").and_then(|v| v.as_int().ok()).ok_or("add_feature needs end")?;
        let strand = Strand::parse(&get("strand").unwrap_or_else(|| "forward".into())).ok_or("bad strand")?;
        let mut feature = Feature::new(name, kind, span_of(&d.doc, start, end)?, strand);
        if let Some(color) = get("color").and_then(|c| Color::parse_hex(&c)) {
            feature.set_color(Some(color));
        }
        if let Some(quals) = spec.get("qualifiers").and_then(|q| q.read_lock::<Map>().map(|m| m.clone())) {
            for (k, v) in quals {
                feature.set_qualifier(k.as_str(), QualValue::Text(v.to_string()));
            }
        }
        d.doc.features.push(feature);
        let index = d.doc.features.len() - 1;
        Ok(feature_dyn(&d.doc.features[index], d.doc.len(), index))
    });
    engine.register_fn("remove_feature", |d: &mut Doc, name: &str| -> i64 {
        let before = d.doc.features.len();
        d.doc.features.retain(|f| !f.name.eq_ignore_ascii_case(name));
        (before - d.doc.features.len()) as i64
    });
    engine.register_fn("update_feature", |d: &mut Doc, name: &str, changes: Map| -> R<()> {
        let (index, _) = find_feature(&d.doc, name)?;
        let seq_len = d.doc.len();
        let circular = d.doc.is_circular();
        let feature = &mut d.doc.features[index];
        for (key, value) in changes {
            match key.as_str() {
                "name" => feature.name = value.to_string(),
                "type" => feature.kind = value.to_string(),
                "strand" => feature.strand = Strand::parse(&value.to_string()).ok_or("bad strand")?,
                "color" => feature.set_color(Some(Color::parse_hex(&value.to_string()).ok_or("colour must be #rrggbb")?)),
                "range" => {
                    let span = core(slopgene_core::span::parse_range(&value.to_string(), seq_len, circular))?;
                    let color = feature.color();
                    feature.segments = vec![Segment::new(span)];
                    feature.set_color(color);
                }
                other => feature.set_qualifier(other, QualValue::Text(value.to_string())),
            }
        }
        Ok(())
    });
    engine.register_fn("add_primer", move |d: &mut Doc, name: &str, seq: &str| -> R<Dynamic> {
        let bases = clean(seq)?;
        d.doc.primers.push(Primer::new(name, String::from_utf8_lossy(&bases).into_owned()));
        primers::refresh_binding_sites(&mut d.doc, &rules);
        let p = d.doc.primers.last().unwrap();
        Ok(to_dyn(json::primer_to_value(p, d.doc.len())))
    });
    engine.register_fn("remove_primer", |d: &mut Doc, name: &str| -> i64 {
        let before = d.doc.primers.len();
        d.doc.primers.retain(|p| !p.name.eq_ignore_ascii_case(name));
        (before - d.doc.primers.len()) as i64
    });
    engine.register_fn("refresh_primers", move |d: &mut Doc| primers::refresh_binding_sites(&mut d.doc, &rules));

    // Saving. The format follows the extension; `save()` writes back to the file it came from.
    engine.register_fn("save", |d: &mut Doc| -> R<String> {
        let path = d.path.clone().ok_or("this document has no file yet; use save(path)")?;
        // Writing back must not drop what the file's format cannot hold, such as features added
        // to a FASTA file. Saving under a new name in another format is an explicit choice.
        let format = io::format_of(&path);
        let lost = io::losses(&d.doc, format);
        if !lost.is_empty() {
            let keep = path.with_extension(Format::SnapGene.extension(d.doc.molecule));
            return err(format!(
                "{} is a {} file, which cannot keep {}; nothing was written. Use save(\"{}\") to keep everything",
                path.display(),
                format.name(),
                lost.join(", "),
                keep.to_string_lossy().replace('\\', "/")
            ));
        }
        d.doc.meta.modified = Some(clock::today_iso());
        core(io::write_path(&d.doc, &path, None))?;
        Ok(path.to_string_lossy().into_owned())
    });
    engine.register_fn("save", |d: &mut Doc, path: &str| -> R<String> {
        d.doc.meta.modified = Some(clock::today_iso());
        core(io::write_path(&d.doc, path, None))?;
        d.path = Some(PathBuf::from(path));
        Ok(path.to_string())
    });
    engine.register_fn("save", |d: &mut Doc, path: &str, format: &str| -> R<String> {
        let format = Format::from_name(format).ok_or("format must be dna, genbank, fasta or json")?;
        d.doc.meta.modified = Some(clock::today_iso());
        core(io::write_path(&d.doc, path, Some(format)))?;
        Ok(path.to_string())
    });
    engine.register_fn("render", |d: &mut Doc, path: &str| -> R<String> { render_to(&d.doc, path, Map::new()) });
    engine.register_fn("render", |d: &mut Doc, path: &str, options: Map| -> R<String> { render_to(&d.doc, path, options) });
}

/// Draw a map to SVG or PNG. Options: width, height, enzymes (set or list, "none"), dark,
/// colors ("palette" or "file").
fn render_to(doc: &Document, path: &str, options: Map) -> R<String> {
    let num = |k: &str, default: f32| {
        options.get(k).and_then(|v| v.as_float().ok().map(|f| f as f32).or_else(|| v.as_int().ok().map(|i| i as f32))).unwrap_or(default)
    };
    let enzymes = options.get("enzymes").map(|v| v.to_string()).unwrap_or_else(|| "common".into());
    let marks = if enzymes.eq_ignore_ascii_case("none") { Vec::new() } else { render::unique_cutter_marks(doc, &enzymes_from(&enzymes)?) };
    let dark = options.get("dark").and_then(|v| v.as_bool().ok()).unwrap_or(false);
    let colors = options.get("colors").and_then(|v| render::ColorMode::parse(&v.to_string())).unwrap_or_default();
    let map_options = MapOptions {
        width: num("width", 900.0),
        height: num("height", 700.0),
        theme: if dark { Theme::dark() } else { Theme::light() },
        enzymes: marks,
        colors,
        ..MapOptions::default()
    };
    let scene = render::render_map(doc, &map_options, &ApproxMeasure);
    files::save_scene(&scene, std::path::Path::new(path), 2.0).map_err(|e| e.to_string())?;
    Ok(path.to_string())
}

fn register_analysis(engine: &mut Engine) {
    engine.register_fn("sites", |d: &mut Doc, spec: &str| -> R<Array> {
        let enzymes = enzymes_from(spec)?;
        let n = d.doc.len();
        Ok(find_sites(&d.doc.seq, d.doc.is_circular(), &enzymes)
            .into_iter()
            .filter(|s| !s.cuts.is_empty())
            .map(|s| {
                let (a, b) = Span::new(s.start, s.start + s.len).to_1based(n);
                to_dyn(serde_json::json!({ "enzyme": s.enzyme, "start": a, "end": b, "strand": s.strand.name(), "cut": s.cuts[0].pos.max(1) }))
            })
            .collect())
    });
    engine.register_fn("cut_counts", |d: &mut Doc, spec: &str| -> R<Map> {
        let enzymes = enzymes_from(spec)?;
        let mut out = Map::new();
        for (e, n) in slopgene_core::enzymes::cut_counts(&d.doc.seq, d.doc.is_circular(), &enzymes) {
            out.insert(e.name.as_str().into(), (n as i64).into());
        }
        Ok(out)
    });
    engine.register_fn("unique_cutters", |d: &mut Doc, spec: &str| -> R<Array> {
        let marks = render::unique_cutter_marks(&d.doc, &enzymes_from(spec)?);
        Ok(marks.into_iter().map(|m| to_dyn(serde_json::json!({ "enzyme": m.name, "cut": m.pos.max(1) }))).collect())
    });
    engine.register_fn("digest", |d: &mut Doc, spec: &str| -> R<Array> {
        let enzymes = enzymes_from(spec)?;
        let result = core(digest(&d.doc, &enzymes))?;
        let n = d.doc.len();
        Ok(result
            .fragments
            .into_iter()
            .map(|f| {
                let (a, b) = f.span.to_1based(n);
                let mut m = Map::new();
                m.insert("length".into(), (f.doc.len() as i64).into());
                m.insert("start".into(), (a as i64).into());
                m.insert("end".into(), (b as i64).into());
                m.insert("left_enzyme".into(), f.left_enzyme.clone().map(Dynamic::from).unwrap_or(Dynamic::UNIT));
                m.insert("right_enzyme".into(), f.right_enzyme.clone().map(Dynamic::from).unwrap_or(Dynamic::UNIT));
                m.insert("left_end".into(), f.doc.ends.left.describe().into());
                m.insert("right_end".into(), f.doc.ends.right.describe().into());
                m.insert("doc".into(), Dynamic::from(Doc { doc: f.doc, path: None }));
                Dynamic::from_map(m)
            })
            .collect())
    });
    engine.register_fn("nearest_cuts", |d: &mut Doc, start: i64, end: i64, spec: &str| -> R<Dynamic> { nearest_cuts_of(d, start, end, spec, Map::new()) });
    engine.register_fn("nearest_cuts", |d: &mut Doc, start: i64, end: i64, spec: &str, options: Map| -> R<Dynamic> {
        nearest_cuts_of(d, start, end, spec, options)
    });
    engine.register_fn("excision", |d: &mut Doc, start: i64, end: i64, left: &str, right: &str| -> R<Dynamic> {
        excision_of(d, start, end, left, right, Map::new())
    });
    engine.register_fn("excision", |d: &mut Doc, start: i64, end: i64, left: &str, right: &str, options: Map| -> R<Dynamic> {
        excision_of(d, start, end, left, right, options)
    });
    engine.register_fn("orfs", |d: &mut Doc| -> Array { orfs_of(&d.doc, 75) });
    engine.register_fn("orfs", |d: &mut Doc, min_codons: i64| -> R<Array> { Ok(orfs_of(&d.doc, usize_of(min_codons, "min_codons")?)) });
    engine.register_fn("find", |d: &mut Doc, pattern: &str| -> R<Array> { find_in(&d.doc, pattern, 0) });
    engine.register_fn("find", |d: &mut Doc, pattern: &str, mismatches: i64| -> R<Array> { find_in(&d.doc, pattern, usize_of(mismatches, "mismatches")?) });
    engine.register_fn("find_protein", |d: &mut Doc, motif: &str| -> R<Array> {
        let motif = alphabet::clean_protein(motif).map_err(|e| e.to_string())?;
        let n = d.doc.len();
        Ok(find_protein(&d.doc.seq, d.doc.is_circular(), &motif, GeneticCode::standard())
            .into_iter()
            .map(|h| {
                let (a, b) = h.span.to_1based(n);
                to_dyn(serde_json::json!({ "start": a, "end": b, "strand": h.strand.name() }))
            })
            .collect())
    });
    engine.register_fn("pcr", |d: &mut Doc, fwd: &str, rev: &str| -> R<Doc> {
        let pick = |text: &str, name: &str| -> R<Primer> {
            if let Some(p) = d.doc.primers.iter().find(|p| p.name.eq_ignore_ascii_case(text)) {
                return Ok(p.clone());
            }
            Ok(Primer::new(name, String::from_utf8_lossy(&clean(text)?).into_owned()))
        };
        let product = core(slopgene_core::pcr::pcr(&d.doc, &pick(fwd, "forward")?, &pick(rev, "reverse")?, &BindingRules::default()))?;
        Ok(Doc { doc: product.doc, path: None })
    });
    engine.register_fn("ligate", |a: Doc, b: Doc| -> R<Doc> { Ok(Doc { doc: core(ligate(&a.doc, &b.doc))?, path: None }) });
    engine.register_fn("restriction_clone", |vector: Doc, insert: Doc, spec: &str| -> R<Array> {
        let products = core(slopgene_core::cloning::restriction_clone(&vector.doc, &insert.doc, &enzymes_from(spec)?))?;
        Ok(products.into_iter().map(|p| Dynamic::from(Doc { doc: p.doc, path: None })).collect())
    });
    engine.register_fn("assemble", |fragments: Array| -> R<Doc> { assemble_docs(fragments, 15, true) });
    engine.register_fn("assemble", |fragments: Array, min_overlap: i64| -> R<Doc> { assemble_docs(fragments, usize_of(min_overlap, "min_overlap")?, true) });
    engine.register_fn("golden_gate", |parts: Array, enzyme_name: &str| -> R<Doc> {
        let enzyme = db().get(enzyme_name).ok_or_else(|| format!("unknown enzyme '{enzyme_name}'"))?;
        let docs = docs_of(parts)?;
        Ok(Doc { doc: core(slopgene_core::cloning::golden_gate(&docs, enzyme))?.doc, path: None })
    });
    engine.register_fn("annotate", |d: &mut Doc, library: Array| -> R<Array> { annotate_with(d, library, 95.0) });
    engine.register_fn("annotate", |d: &mut Doc, library: Array, min_identity: f64| -> R<Array> { annotate_with(d, library, min_identity) });
    engine.register_fn("align", |reference: Doc, read: Doc| -> R<Dynamic> {
        let aln =
            slopgene_core::align::align_read(&reference.doc, &read.doc.seq, read.doc.is_circular()).ok_or("the read does not align (no shared 15-mers)")?;
        Ok(to_dyn(aln.to_json()))
    });
    engine.register_fn("circularize", |a: Doc| -> R<Doc> { Ok(Doc { doc: core(circularize(&a.doc))?, path: None }) });
}

/// Add features from a library (an array of file paths or documents) to a document; returns
/// what was added.
fn annotate_with(d: &mut Doc, library: Array, min_identity: f64) -> R<Array> {
    let mut lib = slopgene_core::annotate::FeatureLibrary::new();
    for item in library {
        if let Some(doc) = item.clone().try_cast::<Doc>() {
            lib.add_document(&doc.doc, 30);
        } else {
            for doc in core(io::read_path(item.to_string()))? {
                lib.add_document(&doc, 30);
            }
        }
    }
    let hits = slopgene_core::annotate::annotate(&d.doc, &lib, min_identity / 100.0);
    let len = d.doc.len();
    let mut added = Array::new();
    for hit in hits {
        let index = d.doc.features.len();
        let mut value = json::feature_to_value(&hit.feature, len);
        value["identity"] = serde_json::json!(hit.identity * 100.0);
        value["index"] = serde_json::json!(index + 1);
        added.push(to_dyn(value));
        d.doc.features.push(hit.feature);
    }
    Ok(added)
}

/// Documents from a script array.
fn docs_of(items: Array) -> R<Vec<Document>> {
    items.into_iter().map(|v| v.try_cast::<Doc>().map(|d| d.doc).ok_or_else(|| "expected an array of documents".into())).collect()
}

/// The `vector` option of the cloning functions: the document an insert goes into.
fn vector_option(options: &Map) -> R<Option<Doc>> {
    match options.get("vector") {
        Some(value) => Ok(Some(value.clone().try_cast::<Doc>().ok_or("vector must be a document")?)),
        None => Ok(None),
    }
}

/// The enzymes that cut nearest a region of a document on each side, nearest first (see
/// `slopgene_core::excise`). Options: `vector` (the document the insert goes into, so each
/// enzyme says how often it cuts it) and `count` (how many enzymes on each side).
fn nearest_cuts_of(d: &Doc, start: i64, end: i64, spec: &str, options: Map) -> R<Dynamic> {
    let region = span_of(&d.doc, start, end)?;
    let enzymes = enzymes_from(spec)?;
    let count = match options.get("count") {
        Some(value) => usize_of(value.as_int().map_err(|_| "count must be a number")?, "count")?,
        None => slopgene_core::excise::DEFAULT_PER_SIDE,
    };
    let vector = vector_option(&options)?;
    let found = slopgene_core::excise::nearest_cuts(&d.doc, region, &enzymes, vector.as_ref().map(|v| &v.doc), count);
    Ok(to_dyn(found.to_json()))
}

/// What cutting a document with `left` before a region and `right` after it gives: the piece
/// that holds the region, with its ends. Option: `vector`, on which the pair is checked too.
fn excision_of(d: &Doc, start: i64, end: i64, left: &str, right: &str, options: Map) -> R<Dynamic> {
    let region = span_of(&d.doc, start, end)?;
    let enzyme = |name: &str| db().get(name).ok_or_else(|| -> Box<EvalAltResult> { format!("unknown enzyme '{name}'").into() });
    let vector = vector_option(&options)?;
    let cut = core(slopgene_core::excise::excision(&d.doc, region, enzyme(left)?, enzyme(right)?, vector.as_ref().map(|v| &v.doc)))?;
    Ok(to_dyn(cut.to_json(d.doc.len())))
}

fn assemble_docs(fragments: Array, min_overlap: usize, circular: bool) -> R<Doc> {
    let docs = docs_of(fragments)?;
    Ok(Doc { doc: core(slopgene_core::cloning::homology_assembly(&docs, min_overlap, circular))?.doc, path: None })
}

fn orfs_of(doc: &Document, min_codons: usize) -> Array {
    let n = doc.len();
    let options = OrfOptions { min_codons, ..OrfOptions::default() };
    find_orfs(&doc.seq, doc.is_circular(), &options)
        .into_iter()
        .map(|o| {
            let (a, b) = o.span.to_1based(n);
            to_dyn(serde_json::json!({ "start": a, "end": b, "strand": o.strand.name(), "frame": o.frame, "codons": o.codons }))
        })
        .collect()
}

fn find_in(doc: &Document, pattern: &str, mismatches: usize) -> R<Array> {
    let probe = clean(pattern)?;
    let n = doc.len();
    Ok(find_dna(&doc.seq, doc.is_circular(), &probe, mismatches, true)
        .into_iter()
        .map(|h| {
            let (a, b) = h.span.to_1based(n);
            to_dyn(serde_json::json!({ "start": a, "end": b, "strand": h.strand.name(), "mismatches": h.mismatches }))
        })
        .collect())
}

/// Cloning protocol files: run one, write a new one, and make the changes the app's bench makes.
/// Each returns the run as `slopgene protocol run --json` reports it.
fn register_protocols(engine: &mut Engine) {
    engine.register_fn("protocol_run", |path: &str| -> R<Dynamic> {
        let path = PathBuf::from(path);
        let (plan, ran) = core(protocol::run_file(&path, db()))?;
        Ok(run_value(&path, &plan, &ran, None))
    });
    engine.register_fn("protocol_new", |path: &str, vector: &str, insert: &str| -> R<Dynamic> { protocol_new(path, vector, insert, Map::new()) });
    engine.register_fn("protocol_new", protocol_new);
    engine.register_fn("protocol_add", |path: &str, map: &str| -> R<Dynamic> { protocol_add(path, map, Map::new()) });
    engine.register_fn("protocol_add", protocol_add);
    engine.register_fn("protocol_cut", |path: &str, step: &str, options: Map| -> R<Dynamic> {
        // A part of "" forgets the part, and enzymes of "" take the enzymes off.
        let part = text_option(&options, "part").map(|p| Some(p).filter(|p| !p.trim().is_empty()));
        let enzymes = text_option(&options, "enzymes").map(|e| enzyme_names(&e)).transpose()?;
        change_protocol(path, Change::Cut { step: step.into(), part, enzymes })
    });
    engine.register_fn("protocol_keep", |path: &str, pick: &str, piece: i64| -> R<Dynamic> {
        change_protocol(path, Change::Keep { pick: pick.into(), piece: usize_of(piece, "piece")? })
    });
    engine.register_fn("protocol_put", |path: &str, piece: &str, into: &str| -> R<Dynamic> {
        change_protocol(path, Change::Put { piece: piece.into(), into: into.into(), dephosphorylate: true })
    });
    engine.register_fn("protocol_take_out", |path: &str, piece: &str, join: &str| -> R<Dynamic> {
        change_protocol(path, Change::TakeOut { piece: piece.into(), join: join.into() })
    });
    engine.register_fn("protocol_remove_map", |path: &str, step: &str| -> R<Dynamic> { change_protocol(path, Change::RemoveMap { step: step.into() }) });
    engine.register_fn("protocol_leave_out", |path: &str, step: &str| -> R<Dynamic> { change_protocol(path, Change::LeaveOut { step: step.into() }) });
    engine.register_fn("protocol_dephosphorylate", |path: &str, piece: &str, join: &str| -> R<Dynamic> {
        change_protocol(path, Change::Dephosphorylate { piece: piece.into(), join: join.into() })
    });
    engine.register_fn("protocol_turn", |path: &str, piece: &str, join: &str, way: &str| -> R<Dynamic> {
        let way = Orientation::parse(way).ok_or_else(|| format!("the way is forward, reverse or either (got '{way}')"))?;
        change_protocol(path, Change::Turn { piece: piece.into(), join: join.into(), way })
    });
    engine.register_fn("protocol_copy_out", |path: &str, step: &str| -> R<Dynamic> { change_protocol(path, Change::CopyOut { step: step.into() }) });
    engine.register_fn("protocol_cut_out", |path: &str, step: &str| -> R<Dynamic> { change_protocol(path, Change::CutOut { step: step.into() }) });
    engine.register_fn("protocol_primers", |path: &str, step: &str, options: Map| -> R<Dynamic> {
        // A pair number keeps one of the pairs the run designs, forward and reverse give two
        // primers of one's own, tm aims the design at a melting temperature, and clear designs
        // them anew.
        let pair = match options.get("pair").filter(|v| !v.is_unit()) {
            Some(value) => Some(usize_of(value.as_int().map_err(|_| "pair must be a number")?, "pair")?),
            None => None,
        };
        let (forward, reverse) = (text_option(&options, "forward"), text_option(&options, "reverse"));
        let tm = number_option(&options, "tm")?;
        if pair.is_none() && forward.is_none() && reverse.is_none() && tm.is_none() && !flag_option(&options, "clear") {
            return Err("say what to do with the primers: pair, forward and reverse, tm, or clear: true".into());
        }
        change_protocol(path, Change::Primers { step: step.into(), pair, forward, reverse, tm })
    });
}

/// A protocol's run as a script value: what `slopgene protocol run --json` reports, with the
/// file's path and, after a change, what changed.
fn run_value(path: &Path, plan: &Protocol, ran: &Run, changed: Option<String>) -> Dynamic {
    let mut value = ran.to_json(plan);
    value["path"] = serde_json::json!(path.to_string_lossy().replace('\\', "/"));
    if let Some(changed) = changed {
        value["changed"] = serde_json::json!(changed);
    }
    to_dyn(value)
}

/// Make `change` to the protocol file at `path`, write it back and run it. A change that fails
/// leaves the file alone.
fn change_protocol(path: &str, change: Change) -> R<Dynamic> {
    let path = PathBuf::from(path);
    let (plan, ran, said) = core(protocol::change_file(&path, db(), &change))?;
    Ok(run_value(&path, &plan, &ran, Some(said)))
}

/// Enzyme names from a list such as "EcoRI,BamHI", each checked. "" stands for none.
fn enzyme_names(spec: &str) -> R<Vec<String>> {
    if spec.trim().is_empty() {
        return Ok(Vec::new());
    }
    Ok(enzymes_from(spec)?.iter().map(|e| e.name.clone()).collect())
}

/// Option `key` of an options map as text, when the map has it.
fn text_option(options: &Map, key: &str) -> Option<String> {
    options.get(key).filter(|v| !v.is_unit()).map(|v| v.to_string())
}

/// Option `key` of an options map as a switch, off when the map lacks it.
fn flag_option(options: &Map, key: &str) -> bool {
    options.get(key).and_then(|v| v.as_bool().ok()).unwrap_or(false)
}

/// Option `key` of an options map as a number, or none when the map lacks it.
fn number_option(options: &Map, key: &str) -> R<Option<f64>> {
    match options.get(key).filter(|v| !v.is_unit()) {
        None => Ok(None),
        Some(value) => value.as_float().or_else(|_| value.as_int().map(|n| n as f64)).map(Some).map_err(|_| format!("{key} must be a number").into()),
    }
}

/// Write a new restriction-cloning protocol file and run it. Options: `enzymes`,
/// `vector_enzymes`, `cut_out`, `replace`, `name`, `keep_phosphates`, `pcr`, `tm` and `force`, as
/// for `slopgene protocol new`.
fn protocol_new(path: &str, vector: &str, insert: &str, options: Map) -> R<Dynamic> {
    let new = NewFile {
        vector: vector.into(),
        insert: insert.into(),
        enzymes: enzyme_names(&text_option(&options, "enzymes").unwrap_or_default())?,
        vector_enzymes: enzyme_names(&text_option(&options, "vector_enzymes").unwrap_or_default())?,
        cut_out: text_option(&options, "cut_out"),
        replace: text_option(&options, "replace"),
        name: text_option(&options, "name"),
        dephosphorylate: !flag_option(&options, "keep_phosphates"),
        pcr: flag_option(&options, "pcr"),
        tm: number_option(&options, "tm")?,
    };
    let path = PathBuf::from(path);
    let (plan, ran) = core(protocol::write_new(&path, &new, flag_option(&options, "force"), db()))?;
    Ok(run_value(&path, &plan, &ran, None))
}

/// Add a map to cut a piece from, or with `pcr` to copy a piece from. Options: `cut_out` (the
/// part), `enzymes` and `pcr`.
fn protocol_add(path: &str, map: &str, options: Map) -> R<Dynamic> {
    let enzymes = enzyme_names(&text_option(&options, "enzymes").unwrap_or_default())?;
    let pcr = flag_option(&options, "pcr");
    change_protocol(path, Change::AddMap { map: map.into(), cut_out: text_option(&options, "cut_out"), enzymes, pcr })
}

fn register_utilities(engine: &mut Engine) {
    engine.register_fn("revcomp", |seq: &str| String::from_utf8_lossy(&alphabet::revcomp(seq.as_bytes())).into_owned());
    engine.register_fn("complement", |seq: &str| String::from_utf8_lossy(&alphabet::complement(seq.as_bytes())).into_owned());
    engine.register_fn("gc", |seq: &str| alphabet::gc_fraction(seq.as_bytes()) * 100.0);
    engine.register_fn("translate", |seq: &str| String::from_utf8_lossy(&translate(seq.as_bytes(), GeneticCode::standard(), false)).into_owned());
    engine.register_fn("translate", |seq: &str, table: i64| -> R<String> {
        let code = GeneticCode::get(table as u32).ok_or("unknown genetic code")?;
        Ok(String::from_utf8_lossy(&translate(seq.as_bytes(), code, false)).into_owned())
    });
    engine.register_fn("tm", |seq: &str| -> R<f64> { Ok(primers::melting_temperature(&clean(seq)?, &TmConditions::default())) });
    engine.register_fn("primer_report", |seq: &str| -> R<Dynamic> {
        let r = primers::analyze(&clean(seq)?, &TmConditions::default());
        Ok(to_dyn(serde_json::json!({
            "length": r.length, "tm": r.tm, "gc_percent": r.gc_percent, "self_dimer": r.self_dimer,
            "three_prime_dimer": r.three_prime_dimer, "hairpin_stem": r.hairpin_stem, "warnings": r.warnings,
        })))
    });
    engine.register_fn("enzyme", |name: &str| -> R<Dynamic> {
        let e = db().get(name).ok_or_else(|| format!("unknown enzyme '{name}'"))?;
        Ok(to_dyn(serde_json::json!({ "name": e.name, "site": e.site_str(), "site_with_cuts": e.site_with_cuts(), "suppliers": e.supplier_names() })))
    });
    engine.register_fn("enzymes", |spec: &str| -> R<Array> { Ok(enzymes_from(spec)?.into_iter().map(|e| Dynamic::from(e.name.clone())).collect()) });
    engine.register_fn("today", clock::today_iso);
    engine.register_fn("join", |items: Array, separator: &str| -> String { items.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(separator) });
    engine.register_fn("to_json", |value: Dynamic| -> R<String> {
        let v: serde_json::Value = rhai::serde::from_dynamic(&value)?;
        Ok(serde_json::to_string_pretty(&v).unwrap_or_default())
    });
    engine.register_fn("parse_json", |text: &str| -> R<Dynamic> {
        let v: serde_json::Value = serde_json::from_str(text).map_err(|e| e.to_string())?;
        Ok(to_dyn(v))
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn eval(src: &str) -> std::result::Result<Dynamic, String> {
        run(src, RunOptions::default())
    }

    #[test]
    fn build_edit_and_analyse() {
        let result = eval(
            r#"
            let d = new_dna("demo", "AAAAGAATTCAAAAAAAAAAGGATCCAAAA", "circular");
            d.add_feature("site1", "misc_feature", 5, 10, "forward");
            d.insert(1, "TT");
            let f = d.feature("site1");
            [d.length, f.start, f.end, d.cut_counts("EcoRI,BamHI").EcoRI, d.digest("EcoRI,BamHI").len()]
            "#,
        )
        .unwrap();
        let values: Vec<i64> = result.into_array().unwrap().into_iter().map(|v| v.as_int().unwrap()).collect();
        assert_eq!(values, vec![32, 7, 12, 1, 2]);
    }

    #[test]
    fn the_nearest_cuts_and_a_chosen_pair() {
        // A region of A and C between an EcoRI site and a BamHI site, and a vector with both.
        let result = eval(
            r#"
            let donor = new_dna("donor", "CACACACACACACACACGAATTCACACACACACACACACACACACAGGATCCACACACACACACACACAC", "circular");
            let vector = new_dna("vector", "ACACACACACACACGAATTCACACAGGATCCACACACACACACACACACACACACACACACACACACAC", "circular");
            let near = donor.nearest_cuts(22, 47, "common", #{ vector: vector, count: 1 });
            let piece = donor.excision(22, 47, "EcoRI", "BamHI", #{ vector: vector });
            [near.before[0].enzyme, near.before[0].distance, near.after[0].enzyme, near.after[0].vector_cuts, piece.directional, piece.left.cut, piece.vector.dropped]
            "#,
        )
        .unwrap();
        let values: Vec<String> = result.into_array().unwrap().into_iter().map(|v| v.to_string()).collect();
        assert_eq!(values, vec!["EcoRI", "3", "BamHI", "1", "true", "18", "11"]);
        let e = eval(r#"new_dna("x", "ACGTACGT").excision(2, 5, "EcoRI", "Nonsense")"#).unwrap_err();
        assert!(e.contains("unknown enzyme"), "{e}");
    }

    #[test]
    fn a_protocol_file_is_written_changed_and_run() {
        // A donor with geneX between an EcoRI site and a BamHI site, and a vector with both sites,
        // saved by the script into a folder of their own. A and C alone make no other sites.
        let folder = std::env::temp_dir().join(format!("slopgene-script-protocol-{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let filler = |n: usize| "CA".repeat(n);
        let donor = format!("{}GAATTC{}GGATCC{}", filler(50), filler(100), filler(50));
        let vector = format!("{}GAATTC{}GGATCC{}", filler(100), filler(15), filler(100));
        let script = format!(
            r#"
            let dir = "{dir}";
            let donor = new_dna("donor", "{donor}", "circular");
            donor.add_feature("geneX", "CDS", 120, 300, "forward");
            donor.save(`${{dir}}/donor.gb`);
            new_dna("vector", "{vector}", "circular").save(`${{dir}}/vector.gb`);
            let file = `${{dir}}/geneX.json`;
            let made = protocol_new(file, `${{dir}}/vector.gb`, `${{dir}}/donor.gb`, #{{ cut_out: "geneX", enzymes: "EcoRI,BamHI" }});
            let other = protocol_keep(file, "insert", 2);
            let back = protocol_keep(file, "insert", 1);
            let added = protocol_add(file, `${{dir}}/donor.gb`, #{{ cut_out: "geneX", enzymes: "EcoRI,BamHI" }});
            let put = protocol_put(file, "insert_2", "product");
            let out = protocol_take_out(file, "insert_2", "product");
            let removed = protocol_remove_map(file, "donor_2");
            let left = protocol_leave_out(file, "cip");
            let cip = protocol_dephosphorylate(file, "backbone", "product");
            let turned = protocol_turn(file, "insert", "product", "forward");
            let cut = protocol_cut(file, "cut_donor", #{{ enzymes: "" }});
            let ran = protocol_run(file);
            [made.steps[made.steps.len() - 1].status, other.changed.sub_string(0, 22), back.changed.sub_string(0, 22), added.changed, put.changed, out.changed, removed.changed, left.changed, cip.changed, turned.changed, cut.changed, ran.steps[1].status]
            "#,
            dir = folder.to_string_lossy().replace('\\', "/"),
        );
        let result = eval(&script).unwrap();
        let values: Vec<String> = result.into_array().unwrap().into_iter().map(|v| v.to_string()).collect();
        assert_eq!(
            values,
            vec![
                "done",
                "insert keeps piece 2 (",
                "insert keeps piece 1 (",
                "added donor as donor_2; insert_2 keeps its piece",
                "product takes insert_2",
                "took insert_2 out of product",
                "took donor_2 out, 3 step(s) in all",
                "left cip out",
                "cip takes the 5' phosphates off backbone before product",
                "insert goes into product as given",
                "cut_donor has no enzymes; cut_vector follows it",
                "needs_choice",
            ]
        );
        // A change that cannot be made says why.
        let e = eval(&format!(r#"protocol_put("{}/geneX.json", "insert", "cut_donor")"#, folder.to_string_lossy().replace('\\', "/"))).unwrap_err();
        assert!(e.contains("neither a map the protocol opens nor a join"), "{e}");
        let e = eval(&format!(r#"protocol_turn("{}/geneX.json", "insert", "product", "sideways")"#, folder.to_string_lossy().replace('\\', "/"))).unwrap_err();
        assert!(e.contains("forward, reverse or either"), "{e}");
        let _ = std::fs::remove_dir_all(&folder);
    }

    #[test]
    fn a_pcr_lane_is_made_and_changed_by_scripts() {
        // Random bases, so primers anneal in one place, without the sites of EcoRI, BamHI or
        // BglII, which the vector then carries on purpose.
        fn random_dna(n: usize, mut x: u64) -> String {
            (0..n)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    b"ACGT"[(x >> 33) as usize % 4] as char
                })
                .collect()
        }
        let clean = |seed: u64| (seed..).map(|s| random_dna(1000, s)).find(|s| !["GAATTC", "GGATCC", "AGATCT"].iter().any(|site| s.contains(site))).unwrap();
        let donor = clean(1);
        let vector = format!("{}GAATTC{}GGATCC{}", clean(2), &clean(3)[..20], clean(4));
        let folder = std::env::temp_dir().join(format!("slopgene-script-pcr-{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let script = format!(
            r#"
            let dir = "{dir}";
            let donor = new_dna("donor", "{donor}", "circular");
            donor.add_feature("geneX", "CDS", 301, 600, "forward");
            donor.save(`${{dir}}/donor.gb`);
            new_dna("vector", "{vector}", "circular").save(`${{dir}}/vector.gb`);
            let file = `${{dir}}/geneX.json`;
            let made = protocol_new(file, `${{dir}}/vector.gb`, `${{dir}}/donor.gb`, #{{ cut_out: "geneX", enzymes: "EcoRI,BamHI", pcr: true }});
            let copy = made.steps[5];
            let kept = protocol_primers(file, "pcr_donor", #{{ pair: 2 }});
            let own = protocol_primers(file, "pcr_donor", #{{ forward: copy.primers[0].sequence, reverse: copy.primers[1].sequence }});
            let anew = protocol_primers(file, "pcr_donor", #{{ clear: true }});
            let aimed = protocol_primers(file, "pcr_donor", #{{ tm: 56 }});
            let back = protocol_cut_out(file, "donor");
            let again = protocol_copy_out(file, "donor");
            let added = protocol_add(file, `${{dir}}/donor.gb`, #{{ cut_out: "geneX", enzymes: "EcoRI,BamHI", pcr: true }});
            [copy.type, copy.primers.len(), copy.primer_pairs[0].pair, kept.changed.sub_string(0, 29), own.changed, anew.changed, aimed.changed, back.changed, again.changed, added.changed, made.steps[made.steps.len() - 1].status]
            "#,
            dir = folder.to_string_lossy().replace('\\', "/"),
        );
        let result = eval(&script).unwrap();
        let values: Vec<String> = result.into_array().unwrap().into_iter().map(|v| v.to_string()).collect();
        assert_eq!(
            values,
            vec![
                "amplify",
                "2",
                "1",
                "pcr_donor keeps primer pair 2",
                "pcr_donor uses the primers given",
                "pcr_donor designs its primers anew",
                "pcr_donor designs its primers anew, aiming at a Tm of 56.0 °C",
                "cut_pcr cuts the part out of donor again, without a PCR",
                "pcr_donor copies the part out of donor by PCR, and cut_pcr trims its tails",
                "added donor as donor_2, copied out by PCR; insert_2 keeps its piece",
                "done",
            ]
        );
        let e = eval(&format!(r#"protocol_primers("{}/geneX.json", "pcr_donor", #{{}})"#, folder.to_string_lossy().replace('\\', "/"))).unwrap_err();
        assert!(e.contains("say what to do with the primers"), "{e}");
        let _ = std::fs::remove_dir_all(&folder);
    }

    #[test]
    fn errors_are_readable() {
        let e = eval(r#"new_dna("x", "ACGT").subseq(3, 9)"#).unwrap_err();
        assert!(e.contains("outside the sequence"), "{e}");
        let e = eval(r#"new_dna("x", "ACGT").feature("nope")"#).unwrap_err();
        assert!(e.contains("no feature named"), "{e}");
    }

    #[test]
    fn utilities() {
        assert_eq!(eval(r#"revcomp("AACG")"#).unwrap().into_string().unwrap(), "CGTT");
        let tm = eval(r#"tm("GTAAAACGACGGCCAGT")"#).unwrap().as_float().unwrap();
        assert!((50.0..58.0).contains(&tm));
        assert!(!api_text().is_empty());
    }

    #[test]
    fn print_goes_to_the_buffer() {
        let buffer = Arc::new(Mutex::new(String::new()));
        let _ = run(r#"print("hello"); print(new_dna("n", "ACGT").length)"#, RunOptions { output: Some(buffer.clone()), ..RunOptions::default() }).unwrap();
        assert_eq!(buffer.lock().unwrap().as_str(), "hello\n4\n");
    }
}
