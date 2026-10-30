//! Restriction enzymes: listing the database, finding sites, digesting, virtual gels.

use crate::render_cmd::save_scene;
use crate::util::*;
use anyhow::{Result, anyhow, bail};
use rayon::prelude::*;
use serde_json::{Value, json};
use slopgene_core::digest::digest;
use slopgene_core::enzymes::{ENZYME_SETS, Effect, EndType, Enzyme, EnzymeDb, Methylase, blocked_groups, find_sites_in};
use slopgene_core::gel::{GelOptions, Ladder, Lane, render_gel};
use slopgene_core::render::thousands;
use slopgene_core::*;
use std::path::{Path, PathBuf};

fn enzyme_json(e: &Enzyme, db: &EnzymeDb) -> Value {
    let end = match e.end_type() {
        EndType::Blunt => json!({ "type": "blunt", "length": 0 }),
        EndType::FivePrime(n) => json!({ "type": "5'", "length": n }),
        EndType::ThreePrime(n) => json!({ "type": "3'", "length": n }),
    };
    let isoschizomers: Vec<&str> = db.enzymes.iter().filter(|o| o.site == e.site && o.name != e.name).map(|o| o.name.as_str()).collect();
    json!({
        "name": e.name,
        "site": e.site_str(),
        "site_with_cuts": e.site_with_cuts(),
        "length": e.site_len(),
        "ends": end,
        "type_iis": e.cuts_outside_site(),
        "suppliers": e.supplier_names(),
        "isoschizomers": isoschizomers,
        "methylation": methylation_json(e),
    })
}

/// The enzyme's response to each methylase, keyed by the methylase's name.
fn methylation_json(e: &Enzyme) -> Value {
    let mut map = serde_json::Map::new();
    for m in Methylase::ALL {
        map.insert(m.name().to_string(), json!(e.sensitivity.get(m).name()));
    }
    Value::Object(map)
}

/// The methylases that block, impair or are needed by the enzyme, in a few words, or nothing.
fn methylation_notes(e: &Enzyme) -> String {
    Methylase::ALL
        .iter()
        .filter_map(|m| {
            let effect = e.sensitivity.get(*m);
            matches!(effect, Effect::Blocked | Effect::Impaired | Effect::Variable | Effect::Required).then(|| format!("{} {}", m.name(), effect.describe()))
        })
        .collect::<Vec<_>>()
        .join(", ")
}

/// Pick enzymes from `--enzymes`: a set name (common, 6+, neb, golden-gate, ...) or a list.
pub fn choose_enzymes<'a>(db: &'a EnzymeDb, spec: &str) -> Result<Vec<&'a Enzyme>> {
    let enzymes = db.set(spec)?;
    if enzymes.is_empty() {
        bail!("no enzymes selected by '{spec}'");
    }
    Ok(enzymes)
}

pub fn list(set: Option<&str>, search: Option<&str>) -> Result<Output> {
    let db = enzyme_db();
    let mut enzymes = match set {
        Some(s) => choose_enzymes(db, s)?,
        None => db.enzymes.iter().collect(),
    };
    if let Some(q) = search {
        let q = q.to_ascii_uppercase();
        enzymes.retain(|e| e.name.to_ascii_uppercase().contains(&q) || e.site_str().contains(&q));
    }
    let mut rows = vec![vec!["enzyme".into(), "site".into(), "ends".into(), "suppliers".into(), "methylation".into()]];
    let mut values = Vec::new();
    for e in &enzymes {
        let ends = match e.end_type() {
            EndType::Blunt => "blunt".to_string(),
            EndType::FivePrime(n) => format!("5' +{n}"),
            EndType::ThreePrime(n) => format!("3' +{n}"),
        };
        rows.push(vec![e.name.clone(), e.site_with_cuts(), ends, e.suppliers.clone(), methylation_notes(e)]);
        values.push(enzyme_json(e, db));
    }
    let mut text = table(&rows);
    text.push_str(&format!("{} enzymes ({})\n", enzymes.len(), db.source));
    Ok(Output::new(Value::Array(values), text))
}

pub fn info(names: &[String]) -> Result<Output> {
    let db = enzyme_db();
    let enzymes = db.resolve(names)?;
    let mut text = String::new();
    let mut values = Vec::new();
    for e in enzymes {
        let v = enzyme_json(e, db);
        let methylation: Vec<String> = Methylase::ALL.iter().map(|m| format!("{} {}", m.name(), e.sensitivity.get(*m).describe())).collect();
        text.push_str(&format!(
            "{}  {}\n  ends: {}\n  suppliers: {}\n  isoschizomers: {}\n  methylation: {}\n",
            e.name,
            e.site_with_cuts(),
            v["ends"]["type"].as_str().unwrap_or(""),
            e.supplier_names().join(", "),
            v["isoschizomers"].as_array().map(|a| a.iter().filter_map(|x| x.as_str()).collect::<Vec<_>>().join(", ")).unwrap_or_default(),
            methylation.join(", ")
        ));
        values.push(v);
    }
    Ok(Output::new(Value::Array(values), text))
}

pub fn sets() -> Result<Output> {
    let db = enzyme_db();
    let mut rows = vec![vec!["set".into(), "enzymes".into(), "description".into()]];
    let mut values = Vec::new();
    for set in ENZYME_SETS {
        let count = db.enzymes.iter().filter(|e| (set.filter)(e)).count();
        rows.push(vec![set.name.into(), count.to_string(), set.description.into()]);
        values.push(json!({ "name": set.name, "enzymes": count, "description": set.description }));
    }
    Ok(Output::new(Value::Array(values), table(&rows)))
}

pub fn import_rebase(emboss_e: &str, emboss_r: &str, output: &str) -> Result<Output> {
    let e = std::fs::read_to_string(emboss_e)?;
    let r = std::fs::read_to_string(emboss_r)?;
    let mut db = EnzymeDb::from_rebase_emboss(&e, &r)?;
    // The EMBOSS files say nothing about methylation, so the known sensitivities carry over.
    db.keep_sensitivity_of(EnzymeDb::builtin());
    std::fs::write(output, db.to_tsv())?;
    let text = format!("wrote {} enzymes ({}) to {output}\nuse it with: set SLOPGENE_ENZYMES={output}", db.len(), db.source);
    Ok(Output::new(json!({ "enzymes": db.len(), "source": db.source, "output": output }), text))
}

#[allow(clippy::too_many_arguments)]
pub fn sites(
    paths: &[String],
    recursive: bool,
    spec: &str,
    unique: bool,
    min_cuts: Option<usize>,
    max_cuts: Option<usize>,
    noncutters: bool,
) -> Result<Output> {
    let db = enzyme_db();
    let enzymes = choose_enzymes(db, spec)?;
    let paths = expand_paths(paths, recursive)?;
    let many = paths.len() > 1;
    // Load and analyse every record in parallel; results keep the input order.
    let loaded: Vec<(PathBuf, Result<Vec<Document>>)> = paths.par_iter().map(|p| (p.clone(), load_all(p))).collect();
    let mut records = Vec::new();
    for (path, docs) in loaded {
        for doc in docs? {
            records.push((path.clone(), doc));
        }
    }
    let reports: Vec<(Value, String)> =
        records.par_iter().map(|(path, doc)| site_report(path, doc, &enzymes, unique, min_cuts, max_cuts, noncutters, many)).collect();
    let mut values = Vec::with_capacity(reports.len());
    let mut text = String::new();
    for (v, t) in reports {
        values.push(v);
        text.push_str(&t);
    }
    Ok(Output::new(Value::Array(values), text))
}

/// Sites of one document, as JSON and as text.
#[allow(clippy::too_many_arguments)]
fn site_report(
    path: &Path,
    doc: &Document,
    enzymes: &[&Enzyme],
    unique: bool,
    min_cuts: Option<usize>,
    max_cuts: Option<usize>,
    noncutters: bool,
    many: bool,
) -> (Value, String) {
    let seq_len = doc.len();
    let found = find_sites_in(doc, enzymes);
    let by_enzyme = slopgene_core::enzymes::cutting_sites_by_enzyme(&found);
    // Sites that the molecule's methylation keeps from cutting, with the reason for each.
    let mut blocked_sites: Vec<Value> = Vec::new();
    let blocked_text: Vec<String> = blocked_groups(&found).iter().map(|g| g.describe()).collect();
    for s in found.iter().filter(|s| s.blocked.is_some()) {
        let Some(block) = s.blocked else { continue };
        let at = s.blocked_cuts.first().map_or(s.start + 1, |c| c.pos.max(1));
        let (a, b) = Span::new(s.start, s.start + s.len).to_1based(seq_len);
        blocked_sites.push(json!({
            "enzyme": s.enzyme, "start": a, "end": b, "strand": s.strand.name(), "cut": at,
            "blocked_by": block.methylase.name(), "reason": block.describe(),
        }));
    }
    let mut rows = vec![vec!["enzyme".into(), "site".into(), "cuts".into(), "positions".into()]];
    let mut per_enzyme = Vec::new();
    let mut none = Vec::new();
    for e in enzymes {
        let hits: &[&slopgene_core::enzymes::Site] = by_enzyme.get(e.name.as_str()).map(|v| v.as_slice()).unwrap_or(&[]);
        let count = hits.len();
        if count == 0 {
            none.push(e.name.clone());
            continue;
        }
        if unique && count != 1 {
            continue;
        }
        if min_cuts.is_some_and(|m| count < m) || max_cuts.is_some_and(|m| count > m) {
            continue;
        }
        // Positions are top-strand cut positions, as SnapGene labels them.
        let positions: Vec<usize> = hits.iter().map(|s| s.cuts[0].pos.max(1)).collect();
        let site_list: Vec<Value> = hits
            .iter()
            .map(|s| {
                let (a, b) = Span::new(s.start, s.start + s.len).to_1based(seq_len);
                let mut site = json!({ "start": a, "end": b, "strand": s.strand.name(), "cut": s.cuts[0].pos.max(1) });
                if let Some(impaired) = s.impaired {
                    site["impaired"] = json!(impaired.describe());
                }
                site
            })
            .collect();
        rows.push(vec![e.name.clone(), e.site_with_cuts(), count.to_string(), positions.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(", ")]);
        per_enzyme.push(json!({ "name": e.name, "site": e.site_str(), "cuts": count, "positions": positions, "sites": site_list }));
    }
    let mut v = json!({
        "path": path_text(path),
        "name": doc.name,
        "length": seq_len,
        "topology": doc.topology.name(),
        "methylation": doc.methylation.describe(),
        "enzymes": per_enzyme,
        "blocked": blocked_sites,
    });
    if noncutters {
        v["noncutters"] = json!(none);
    }
    let mut text = String::new();
    if many {
        text.push_str(&format!("== {} ({} bp)\n", doc.name, thousands(seq_len)));
    }
    text.push_str(&table(&rows));
    if !blocked_text.is_empty() {
        text.push_str(&format!("blocked by methylation ({}): {}\n", doc.methylation.describe(), blocked_text.join("; ")));
    }
    if noncutters && !none.is_empty() {
        text.push_str(&format!("does not cut: {}\n", none.join(", ")));
    }
    if many {
        text.push('\n');
    }
    (v, text)
}

#[allow(clippy::too_many_arguments)]
pub fn digest_cmd(path: &str, spec: &str, out_dir: Option<&str>, gel: Option<&str>, ladder: &str, agarose: f64, uncut: bool) -> Result<Output> {
    let db = enzyme_db();
    let enzymes = choose_enzymes(db, spec)?;
    let path = PathBuf::from(path);
    let doc = load(&path)?;
    let result = digest(&doc, &enzymes)?;
    let seq_len = doc.len();

    let end_json = |o: Overhang, enzyme: &Option<String>| {
        let (kind, n) = match o {
            Overhang::Blunt => ("blunt", 0),
            Overhang::FivePrime(n) => ("5'", n),
            Overhang::ThreePrime(n) => ("3'", n),
        };
        json!({ "type": kind, "length": n, "enzyme": enzyme })
    };
    let mut rows = vec![vec!["#".into(), "length".into(), "from".into(), "to".into(), "left end".into(), "right end".into()]];
    let mut fragments = Vec::new();
    let mut written = Vec::new();
    for (ii, frag) in result.fragments.iter().enumerate() {
        let (a, b) = frag.span.to_1based(seq_len);
        let features: Vec<&str> = frag.doc.features.iter().map(|f| f.name.as_str()).collect();
        let mut v = json!({
            "index": ii + 1,
            "length": frag.len(),
            "start": a,
            "end": b,
            "left_end": end_json(frag.doc.ends.left, &frag.left_enzyme),
            "right_end": end_json(frag.doc.ends.right, &frag.right_enzyme),
            "features": features,
            "circular": frag.doc.is_circular(),
        });
        if let Some(dir) = out_dir {
            std::fs::create_dir_all(dir)?;
            let name = format!("{}_fragment{}_{}bp.dna", doc.name.replace([' ', '/', '\\'], "_"), ii + 1, frag.len());
            let target = Path::new(dir).join(name);
            save(&frag.doc, &target, None)?;
            v["written"] = json!(path_text(&target));
            written.push(path_text(&target));
        }
        rows.push(vec![
            (ii + 1).to_string(),
            thousands(frag.len()),
            frag.left_enzyme.clone().unwrap_or_else(|| "end".into()),
            frag.right_enzyme.clone().unwrap_or_else(|| "end".into()),
            frag.doc.ends.left.describe(),
            frag.doc.ends.right.describe(),
        ]);
        fragments.push(v);
    }

    let mut value = json!({
        "path": path_text(&path),
        "name": doc.name,
        "enzymes": enzymes.iter().map(|e| e.name.clone()).collect::<Vec<_>>(),
        "fragments": fragments,
        "noncutters": result.non_cutters,
        "methylation": doc.methylation.describe(),
        "blocked": result.blocked.iter().filter_map(|s| {
            let block = s.blocked?;
            let at = s.blocked_cuts.first().map_or(s.start + 1, |c| c.pos.max(1));
            Some(json!({ "enzyme": s.enzyme, "cut": at, "blocked_by": block.methylase.name(), "reason": block.describe() }))
        }).collect::<Vec<_>>(),
    });
    let mut text = table(&rows);
    if !result.non_cutters.is_empty() {
        text.push_str(&format!("does not cut: {}\n", result.non_cutters.join(", ")));
    }
    if !result.blocked.is_empty() {
        let blocked: Vec<String> = blocked_groups(&result.blocked).iter().map(|g| g.describe()).collect();
        text.push_str(&format!("blocked by methylation ({}): {}\n", doc.methylation.describe(), blocked.join("; ")));
    }
    if let Some(gel_path) = gel {
        let ladder = Ladder::get(ladder).ok_or_else(|| anyhow!("unknown ladder '{ladder}' (1kb, 1kb-plus, 100bp, generuler-1kb, generuler-1kb-plus)"))?;
        let label = enzymes.iter().map(|e| e.name.as_str()).collect::<Vec<_>>().join("+");
        let mut lanes = vec![Lane::ladder(ladder)];
        if uncut && doc.is_circular() {
            lanes.push(Lane::uncut_plasmid("uncut", seq_len));
        }
        let sizes: Vec<usize> = result.fragments.iter().map(|f| f.len()).collect();
        lanes.push(if result.fragments.len() == 1 && result.fragments[0].doc.is_circular() {
            Lane::uncut_plasmid(label, seq_len)
        } else {
            Lane::digest(label, &sizes)
        });
        let scene = render_gel(&lanes, &GelOptions { agarose_percent: agarose, ..GelOptions::default() });
        save_scene(&scene, Path::new(gel_path), 2.0)?;
        value["gel"] = json!(gel_path);
        text.push_str(&format!("gel image: {gel_path}\n"));
    }
    if !written.is_empty() {
        text.push_str(&format!("wrote {} fragment files to {}\n", written.len(), out_dir.unwrap_or("")));
    }
    Ok(Output::new(value, text))
}
