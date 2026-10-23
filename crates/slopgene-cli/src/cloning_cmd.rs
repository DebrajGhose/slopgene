//! Cloning commands: restriction cloning, homology assembly (Gibson, HiFi, In-Fusion), Golden
//! Gate, and launching the desktop app.

use crate::enzyme_cmd::choose_enzymes;
use crate::util::*;
use anyhow::{Result, anyhow, bail};
use serde_json::{Value, json};
use slopgene_core::cloning::{Product, golden_gate, homology_assembly, restriction_clone};
use slopgene_core::excise::{self, NearCut};
use slopgene_core::render::thousands;
use slopgene_core::*;
use std::path::{Path, PathBuf};

/// What `slopgene excise` looks at: the region to cut out of a map, the vector it goes into,
/// and the pair of enzymes chosen, if any.
pub struct ExcisePlan<'a> {
    /// The map that carries the insert.
    pub path: &'a str,
    /// The region, as a range or a feature name.
    pub range: Option<&'a str>,
    pub feature: Option<&'a str>,
    /// The vector the insert goes into.
    pub vector: Option<&'a str>,
    /// The enzymes to choose from: a set name or a list.
    pub enzymes: &'a str,
    /// How many enzymes to list on each side.
    pub count: usize,
    /// The enzymes chosen to cut before and after the region.
    pub left: Option<&'a str>,
    pub right: Option<&'a str>,
}

/// The enzymes that cut nearest a region on each side, nearest first, or with `left` and
/// `right` the piece that pair cuts out (see [`slopgene_core::excise`]).
pub fn excise(plan: &ExcisePlan) -> Result<Output> {
    let db = enzyme_db();
    let donor = load(Path::new(plan.path))?;
    if plan.range.is_none() && plan.feature.is_none() {
        bail!("say which region to cut out: --range 1201..1950, or --feature geneX");
    }
    let region = resolve_region(&donor, plan.range, plan.feature)?;
    let vector = plan.vector.map(|p| load(Path::new(p))).transpose()?;
    let len = donor.len();

    let (first, last) = region.span.to_1based(len);
    let mut value = json!({
        "path": path_text(Path::new(plan.path)),
        "name": donor.name,
        "region": { "start": first, "end": last, "length": region.span.len(), "label": region.label },
    });
    if let Some(v) = &vector {
        value["vector"] = json!({ "path": path_text(Path::new(plan.vector.unwrap_or(""))), "name": v.name });
    }

    // A chosen pair: the piece it cuts out, with its ends, and how it opens the vector.
    if let (Some(left), Some(right)) = (plan.left, plan.right) {
        let enzyme = |name: &str| db.get(name).ok_or_else(|| anyhow!("unknown enzyme '{name}'"));
        let cut = excise::excision(&donor, region.span, enzyme(left)?, enzyme(right)?, vector.as_ref())?;
        value["excision"] = cut.to_json(len);
        let mut text = format!("{} of {}, cut with {} before it and {} after it\n", region.label, donor.name, left, right);
        let mut rows = vec![
            vec!["piece".to_string(), format!("{} ({} bp)", cut.fragment.format(len), thousands(cut.length))],
            vec!["ends".to_string(), cut.ends_text()],
            vec!["goes into a vector".to_string(), cut.direction_text().to_string()],
            vec!["bases beyond the region".to_string(), cut.extra_text()],
            vec!["other cuts in the map".to_string(), cut.other_cuts.to_string()],
        ];
        if let Some(fit) = &cut.vector {
            rows.push(vec!["vector".to_string(), fit.describe()]);
        }
        text.push_str(&table(&rows));
        return Ok(Output::new(value, text));
    }

    // Otherwise the nearest cuts on each side, for a person to choose one on each.
    let enzymes = choose_enzymes(db, plan.enzymes)?;
    let found = excise::nearest_cuts(&donor, region.span, &enzymes, vector.as_ref(), plan.count);
    let listed = found.to_json();
    for key in ["before", "after", "inside"] {
        value[key] = listed[key].clone();
    }
    value["enzymes"] = json!(plan.enzymes);
    let mut text = format!("{} of {}: the {} enzymes that cut nearest it on each side\n", region.label, donor.name, plan.enzymes);
    for (title, side) in [("before it", &found.before), ("after it", &found.after)] {
        text.push_str(&format!("\n{title}\n"));
        if side.is_empty() {
            text.push_str("none\n");
            continue;
        }
        let mut rows = vec![vec!["enzyme".to_string(), "distance".into(), "cuts after".into(), "end".into(), "cuts in map".into()]];
        if vector.is_some() {
            rows[0].push("cuts in vector".into());
        }
        for cut in side.iter() {
            rows.push(near_row(cut));
        }
        text.push_str(&table(&rows));
    }
    if !found.inside.is_empty() {
        text.push_str(&format!("\ncutting inside the region, so not usable: {}\n", found.inside.join(", ")));
    }
    text.push_str("\npick one enzyme on each side, then see what they cut out with --left A --right B\n");
    Ok(Output::new(value, text))
}

/// One cut as a row of the listing: the enzyme, its distance, where it cuts, its end, and how
/// often it cuts the map and the vector.
fn near_row(cut: &NearCut) -> Vec<String> {
    let mut row = vec![cut.enzyme.clone(), format!("{} bp", thousands(cut.distance)), thousands(cut.pos), cut.end_text(), cut.cuts.to_string()];
    if let Some(n) = cut.vector_cuts {
        row.push(n.to_string());
    }
    row
}

fn product_json(p: &Product) -> Value {
    json!({
        "name": p.doc.name,
        "length": p.doc.len(),
        "topology": p.doc.topology.name(),
        "features": p.doc.features.iter().map(|f| f.name.clone()).collect::<Vec<_>>(),
        "notes": p.notes,
    })
}

/// Save one or several products. Several products from one call get `_1`, `_2` suffixes.
fn save_products(products: &mut [Product], output: Option<&str>) -> Result<Vec<String>> {
    let Some(out) = output else { return Ok(Vec::new()) };
    let mut written = Vec::new();
    let count = products.len();
    for (ii, p) in products.iter_mut().enumerate() {
        let path = if count == 1 {
            PathBuf::from(out)
        } else {
            let base = Path::new(out);
            let stem = base.file_stem().and_then(|s| s.to_str()).unwrap_or("product");
            let ext = base.extension().and_then(|s| s.to_str()).unwrap_or("dna");
            base.with_file_name(format!("{stem}_{}.{ext}", ii + 1))
        };
        p.doc.meta.created = Some(today_iso());
        save(&p.doc, &path, None)?;
        written.push(path_text(&path));
    }
    Ok(written)
}

fn report(mut products: Vec<Product>, output: Option<&str>) -> Result<Output> {
    let written = save_products(&mut products, output)?;
    let mut text = String::new();
    for (ii, p) in products.iter().enumerate() {
        text.push_str(&format!("product {}: {} ({} bp, {})\n", ii + 1, p.doc.name, thousands(p.doc.len()), p.doc.topology.name()));
        for note in &p.notes {
            text.push_str(&format!("  {note}\n"));
        }
    }
    for w in &written {
        text.push_str(&format!("wrote {w}\n"));
    }
    let mut values: Vec<Value> = products.iter().map(product_json).collect();
    for (v, w) in values.iter_mut().zip(&written) {
        v["written"] = json!(w);
    }
    if output.is_none() {
        for (v, p) in values.iter_mut().zip(&products) {
            v["sequence"] = json!(p.doc.seq_str());
        }
    }
    Ok(Output::new(Value::Array(values), text))
}

pub fn clone_cmd(vector: &str, insert: &str, enzymes: &str, output: Option<&str>) -> Result<Output> {
    let db = enzyme_db();
    let chosen = choose_enzymes(db, enzymes)?;
    let vector = load(Path::new(vector))?;
    let insert = load(Path::new(insert))?;
    let products = restriction_clone(&vector, &insert, &chosen)?;
    report(products, output)
}

pub fn assemble(fragments: &[String], min_overlap: usize, linear: bool, output: Option<&str>) -> Result<Output> {
    if fragments.is_empty() {
        bail!("give the fragments in order");
    }
    let docs: Vec<Document> = fragments.iter().map(|f| load(Path::new(f))).collect::<Result<_>>()?;
    let product = homology_assembly(&docs, min_overlap, !linear)?;
    report(vec![product], output)
}

pub fn golden_gate_cmd(parts: &[String], enzyme: &str, output: Option<&str>) -> Result<Output> {
    let db = enzyme_db();
    let enzyme = db.get(enzyme).ok_or_else(|| anyhow::anyhow!("unknown enzyme '{enzyme}'"))?;
    let docs: Vec<Document> = parts.iter().map(|f| load(Path::new(f))).collect::<Result<_>>()?;
    let product = golden_gate(&docs, enzyme)?;
    report(vec![product], output)
}

/// Open files in the desktop app (`slopgene-gui`, next to this program or on the PATH).
pub fn view(files: &[String]) -> Result<Output> {
    let exe = if cfg!(windows) { "slopgene-gui.exe" } else { "slopgene-gui" };
    let sibling = std::env::current_exe().ok().and_then(|p| p.parent().map(|d| d.join(exe)));
    let program = match sibling {
        Some(p) if p.exists() => p,
        _ => PathBuf::from(exe),
    };
    std::process::Command::new(&program).args(files).spawn().map_err(|e| anyhow::anyhow!("could not start {}: {e}", program.display()))?;
    Ok(Output::new(json!({ "ok": true, "opened": files }), format!("opened {} file(s) in SlopGene", files.len())))
}

/// Annotate targets with features from a library of trusted, annotated files.
#[allow(clippy::too_many_arguments)]
pub fn annotate(
    targets: &[String],
    recursive: bool,
    library_files: &[String],
    min_identity: f64,
    min_length: usize,
    output: Option<&str>,
    in_place: bool,
    save_library: Option<&str>,
) -> Result<Output> {
    use slopgene_core::annotate::{FeatureLibrary, annotate};
    let mut library = FeatureLibrary::new();
    for path in expand_paths(library_files, recursive)? {
        for doc in load_all(&path)? {
            library.add_document(&doc, min_length);
        }
    }
    if library.is_empty() {
        bail!("the library files contain no features of at least {min_length} bp");
    }
    if let Some(out) = save_library {
        save(&library.to_document("feature library"), Path::new(out), None)?;
    }
    let targets = expand_paths(targets, recursive)?;
    if output.is_some() && targets.len() > 1 {
        bail!("-o names one output file; with several targets use --in-place");
    }
    let mut values = Vec::new();
    let mut text = format!("library: {} features\n", library.len());
    for path in &targets {
        let mut doc = load(path)?;
        let hits = annotate(&doc, &library, min_identity / 100.0);
        let len = doc.len();
        let mut rows = vec![vec!["feature".into(), "type".into(), "range".into(), "strand".into(), "identity".into(), "from".into()]];
        let mut found = Vec::new();
        for hit in &hits {
            let span = hit.feature.span(len);
            rows.push(vec![
                hit.feature.name.clone(),
                hit.feature.kind.clone(),
                span.format(len),
                hit.feature.strand.name().into(),
                format!("{:.1}%", hit.identity * 100.0),
                hit.source.clone(),
            ]);
            let (start, end) = span.to_1based(len);
            found.push(json!({ "name": hit.feature.name, "type": hit.feature.kind, "start": start, "end": end, "strand": hit.feature.strand.name(), "identity": (hit.identity * 1000.0).round() / 10.0, "source": hit.source }));
        }
        text.push_str(&format!("== {} ({} new features)\n", doc.name, hits.len()));
        if !hits.is_empty() {
            text.push_str(&table(&rows));
        }
        let mut value = json!({ "path": path_text(path), "name": doc.name, "added": found });
        let target = match (output, in_place) {
            (Some(o), _) => Some(PathBuf::from(o)),
            (None, true) => Some(path.clone()),
            _ => None,
        };
        if let Some(target) = target {
            doc.features.extend(hits.into_iter().map(|h| h.feature));
            if in_place && output.is_none() {
                refuse_lossy_in_place(&doc, &target)?;
            }
            doc.meta.modified = Some(today_iso());
            save(&doc, &target, None)?;
            value["written"] = json!(path_text(&target));
            text.push_str(&format!("wrote {}\n", path_text(&target)));
        }
        values.push(value);
    }
    Ok(Output::new(Value::Array(values), text))
}
