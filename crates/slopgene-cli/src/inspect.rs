//! Read-only commands: summaries, sequences, features, primers, searches, ORFs, translation.

use crate::util::*;
use anyhow::{Result, bail};
use rayon::prelude::*;
use serde_json::{Value, json};
use slopgene_core::render::thousands;
use slopgene_core::search::{find_dna, find_protein};
use slopgene_core::translate::{GeneticCode, OrfOptions, StartCodons, find_orfs, translate, translate_feature};
use slopgene_core::*;
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Load every record of every file in parallel, keeping input order.
fn load_many(paths: &[PathBuf]) -> Vec<(PathBuf, Result<Vec<Document>>)> {
    paths.par_iter().map(|p| (p.clone(), load_all(p))).collect()
}

/// One-line summary of a document, as JSON.
pub fn summary(path: &std::path::Path, doc: &Document) -> Value {
    let mut types: BTreeMap<String, usize> = BTreeMap::new();
    for f in &doc.features {
        *types.entry(f.kind.clone()).or_default() += 1;
    }
    let mut obj = json!({
        "path": path_text(path),
        "name": doc.name,
        "molecule": doc.molecule.name(),
        "topology": doc.topology.name(),
        "length": doc.len(),
        "features": doc.features.len(),
        "primers": doc.primers.len(),
        "feature_types": types,
        "description": doc.meta.description,
        "mass_da": (stats::document_mass(doc) * 10.0).round() / 10.0,
    });
    if doc.molecule.is_nucleic() {
        obj["gc_percent"] = json!((alphabet::gc_fraction(&doc.seq) * 1000.0).round() / 10.0);
    }
    if doc.molecule == Molecule::Dna && !doc.is_circular() && (doc.ends.left != Overhang::Blunt || doc.ends.right != Overhang::Blunt) {
        obj["ends"] = json!({ "left": doc.ends.left.describe(), "right": doc.ends.right.describe() });
    }
    obj
}

pub fn info(paths: &[String], recursive: bool) -> Result<Output> {
    let paths = expand_paths(paths, recursive)?;
    let mut rows = vec![vec!["file".into(), "name".into(), "length".into(), "topology".into(), "GC".into(), "features".into(), "primers".into()]];
    let mut values = Vec::new();
    let mut errors = Vec::new();
    for (path, result) in load_many(&paths) {
        match result {
            Ok(docs) => {
                for doc in &docs {
                    let s = summary(&path, doc);
                    rows.push(vec![
                        path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default(),
                        doc.name.clone(),
                        format!("{} {}", thousands(doc.len()), if doc.molecule == Molecule::Protein { "aa" } else { "bp" }),
                        doc.topology.name().into(),
                        s.get("gc_percent").map(|g| format!("{g}%")).unwrap_or_default(),
                        doc.features.len().to_string(),
                        doc.primers.len().to_string(),
                    ]);
                    values.push(s);
                }
            }
            Err(e) => {
                errors.push(format!("{e:#}"));
                values.push(json!({ "path": path_text(&path), "error": format!("{e:#}") }));
            }
        }
    }
    let mut text = table(&rows);
    for e in &errors {
        text.push_str(&format!("error: {e}\n"));
    }
    Ok(Output::new(Value::Array(values), text))
}

pub fn show(path: &str, include_sequence: bool) -> Result<Output> {
    let path = PathBuf::from(path);
    let doc = load(&path)?;
    let mut value = json::document_to_value(&doc, include_sequence);
    value["path"] = json!(path_text(&path));

    let seq_len = doc.len();
    let mut text = format!(
        "{}  {} {}  {}  {}\n",
        doc.name,
        thousands(seq_len),
        if doc.molecule == Molecule::Protein { "aa" } else { "bp" },
        doc.topology.name(),
        doc.molecule.name().to_uppercase()
    );
    if !doc.meta.description.is_empty() {
        text.push_str(&format!("{}\n", doc.meta.description.lines().next().unwrap_or("")));
    }
    if !doc.features.is_empty() {
        text.push('\n');
        text.push_str(&features_table(&doc));
    }
    if !doc.primers.is_empty() {
        text.push('\n');
        text.push_str(&primers_table(&doc));
    }
    Ok(Output::new(value, text))
}

fn features_table(doc: &Document) -> String {
    let seq_len = doc.len();
    let mut rows = vec![vec!["#".into(), "name".into(), "type".into(), "range".into(), "strand".into(), "length".into()]];
    for (ii, f) in doc.features.iter().enumerate() {
        rows.push(vec![(ii + 1).to_string(), f.name.clone(), f.kind.clone(), f.span(seq_len).format(seq_len), f.strand.name().into(), thousands(f.len())]);
    }
    table(&rows)
}

fn primers_table(doc: &Document) -> String {
    let seq_len = doc.len();
    let mut rows = vec![vec!["primer".into(), "length".into(), "binds".into(), "sequence".into()]];
    for p in &doc.primers {
        let sites: Vec<String> = p
            .sites
            .iter()
            .filter(|s| s.extra.attr("simplified") != Some("1"))
            .map(|s| format!("{} {}", s.span.format(seq_len), if s.strand.is_reverse() { "(-)" } else { "(+)" }))
            .collect();
        rows.push(vec![p.name.clone(), p.seq.len().to_string(), if sites.is_empty() { "-".into() } else { sites.join(", ") }, p.seq.clone()]);
    }
    table(&rows)
}

#[allow(clippy::too_many_arguments)]
pub fn seq(
    path: &str,
    range: Option<&str>,
    feature: Option<&str>,
    revcomp: bool,
    fasta: bool,
    translate_it: bool,
    table_id: u32,
    frame: usize,
) -> Result<Output> {
    let doc = load(&PathBuf::from(path))?;
    let region = resolve_region(&doc, range, feature)?;
    let mut bases = match &region.feature {
        Some(f) => f.extract(&doc.seq),
        None => doc.subseq(region.span),
    };
    if revcomp {
        bases = alphabet::revcomp(&bases);
    }
    let mut value = json!({
        "name": doc.name,
        "region": region.label,
        "length": bases.len(),
        "sequence": String::from_utf8_lossy(&bases),
    });
    let mut text = if fasta {
        let mut out = format!(">{} {}\n", doc.name.replace(' ', "_"), region.label);
        for chunk in bases.chunks(70) {
            out.push_str(&String::from_utf8_lossy(chunk));
            out.push('\n');
        }
        out
    } else {
        String::from_utf8_lossy(&bases).into_owned()
    };
    if translate_it {
        let code = GeneticCode::get(table_id).ok_or_else(|| anyhow::anyhow!("unknown genetic code {table_id}"))?;
        let start = frame.clamp(1, 3) - 1;
        let protein = translate(&bases[start.min(bases.len())..], code, false);
        let protein = String::from_utf8_lossy(&protein).into_owned();
        value["protein"] = json!(protein);
        text = protein;
    }
    Ok(Output::new(value, text))
}

pub fn features(paths: &[String], recursive: bool, kind: Option<&str>, name: Option<&str>) -> Result<Output> {
    let paths = expand_paths(paths, recursive)?;
    let mut values = Vec::new();
    let mut text = String::new();
    let many = paths.len() > 1;
    for (path, result) in load_many(&paths) {
        let docs = result?;
        for doc in &docs {
            let seq_len = doc.len();
            let mut rows = vec![vec!["#".into(), "name".into(), "type".into(), "range".into(), "strand".into(), "length".into()]];
            for (ii, f) in doc.features.iter().enumerate() {
                if kind.is_some_and(|k| !f.kind.eq_ignore_ascii_case(k)) {
                    continue;
                }
                if name.is_some_and(|n| !f.name.to_lowercase().contains(&n.to_lowercase())) {
                    continue;
                }
                let mut v = json::feature_to_value(f, seq_len);
                v["index"] = json!(ii + 1);
                v["path"] = json!(path_text(&path));
                v["document"] = json!(doc.name);
                values.push(v);
                rows.push(vec![
                    (ii + 1).to_string(),
                    f.name.clone(),
                    f.kind.clone(),
                    f.span(seq_len).format(seq_len),
                    f.strand.name().into(),
                    thousands(f.len()),
                ]);
            }
            if many {
                text.push_str(&format!("== {} ({})\n", doc.name, path_text(&path)));
            }
            text.push_str(&table(&rows));
            if many {
                text.push('\n');
            }
        }
    }
    Ok(Output::new(Value::Array(values), text))
}

pub fn primers(paths: &[String], recursive: bool, refresh: bool) -> Result<Output> {
    let paths = expand_paths(paths, recursive)?;
    let conditions = primers::TmConditions::default();
    let mut values = Vec::new();
    let mut text = String::new();
    for (path, result) in load_many(&paths) {
        for mut doc in result? {
            if refresh {
                primers::refresh_binding_sites(&mut doc, &primers::BindingRules::default());
            }
            let seq_len = doc.len();
            for p in &doc.primers {
                let report = primers::analyze(p.seq.as_bytes(), &conditions);
                let mut v = json::primer_to_value(p, seq_len);
                v["path"] = json!(path_text(&path));
                v["document"] = json!(doc.name);
                v["tm"] = json!((report.tm * 10.0).round() / 10.0);
                v["gc_percent"] = json!((report.gc_percent * 10.0).round() / 10.0);
                v["warnings"] = json!(report.warnings);
                values.push(v);
            }
            if !doc.primers.is_empty() {
                if paths.len() > 1 {
                    text.push_str(&format!("== {}\n", doc.name));
                }
                text.push_str(&primers_table(&doc));
            }
        }
    }
    Ok(Output::new(Value::Array(values), text))
}

pub fn find(paths: &[String], recursive: bool, pattern: &str, mismatches: usize, protein: bool, one_strand: bool, table_id: u32) -> Result<Output> {
    let paths = expand_paths(paths, recursive)?;
    let code = GeneticCode::get(table_id).ok_or_else(|| anyhow::anyhow!("unknown genetic code {table_id}"))?;
    let probe: Vec<u8> = if protein { alphabet::clean_protein(pattern).map_err(|e| anyhow::anyhow!(e))? } else { clean_dna(pattern)? };
    if probe.is_empty() {
        bail!("empty search pattern");
    }
    let mut values = Vec::new();
    let mut rows = vec![vec!["file".into(), "match".into(), "strand".into(), "mismatches".into(), "bases".into()]];
    for (path, result) in load_many(&paths) {
        for doc in result? {
            let hits = if protein {
                find_protein(&doc.seq, doc.is_circular(), &probe, code)
            } else {
                find_dna(&doc.seq, doc.is_circular(), &probe, mismatches, !one_strand)
            };
            let seq_len = doc.len();
            for hit in hits {
                let mut matched = hit.span.extract(&doc.seq);
                if hit.strand == Strand::Reverse {
                    matched = alphabet::revcomp(&matched);
                }
                let (start, end) = hit.span.to_1based(seq_len);
                values.push(json!({
                    "path": path_text(&path),
                    "document": doc.name,
                    "start": start,
                    "end": end,
                    "strand": hit.strand.name(),
                    "mismatches": hit.mismatches,
                    "matched": String::from_utf8_lossy(&matched),
                }));
                rows.push(vec![
                    doc.name.clone(),
                    hit.span.format(seq_len),
                    hit.strand.name().into(),
                    hit.mismatches.to_string(),
                    String::from_utf8_lossy(&matched).into_owned(),
                ]);
            }
        }
    }
    let text = if rows.len() == 1 { "no matches".to_string() } else { table(&rows) };
    Ok(Output::new(Value::Array(values), text))
}

pub fn orfs(path: &str, min_codons: usize, starts: &str, table_id: u32, one_strand: bool, with_protein: bool) -> Result<Output> {
    let doc = load(&PathBuf::from(path))?;
    let starts = match starts.to_ascii_lowercase().as_str() {
        "atg" => StartCodons::Atg,
        "alt" | "alternative" => StartCodons::Alternative,
        "any" | "stop-to-stop" => StartCodons::Any,
        other => bail!("--starts must be atg, alternative or any (got '{other}')"),
    };
    let options = OrfOptions { min_codons, starts, genetic_code: table_id, both_strands: !one_strand };
    let code = GeneticCode::get(table_id).ok_or_else(|| anyhow::anyhow!("unknown genetic code {table_id}"))?;
    let found = find_orfs(&doc.seq, doc.is_circular(), &options);
    let seq_len = doc.len();
    let mut values = Vec::new();
    let mut rows = vec![vec!["#".into(), "range".into(), "strand".into(), "frame".into(), "codons".into(), "kDa".into()]];
    for (ii, orf) in found.iter().enumerate() {
        let mut coding = orf.span.extract(&doc.seq);
        if orf.strand == Strand::Reverse {
            coding = alphabet::revcomp(&coding);
        }
        let protein = translate(&coding, code, true);
        let kda = stats::protein_mass(&protein) / 1000.0;
        let (start, end) = orf.span.to_1based(seq_len);
        let mut v = json!({
            "index": ii + 1, "start": start, "end": end, "strand": orf.strand.name(),
            "frame": orf.frame, "codons": orf.codons, "mass_kda": (kda * 100.0).round() / 100.0,
        });
        if with_protein {
            v["protein"] = json!(String::from_utf8_lossy(&protein));
        }
        values.push(v);
        rows.push(vec![
            (ii + 1).to_string(),
            orf.span.format(seq_len),
            orf.strand.name().into(),
            orf.frame.to_string(),
            orf.codons.to_string(),
            format!("{kda:.1}"),
        ]);
    }
    let text = if found.is_empty() { format!("no ORFs of at least {min_codons} codons") } else { table(&rows) };
    Ok(Output::new(Value::Array(values), text))
}

pub fn translate_cmd(
    path: Option<&str>,
    seq_text: Option<&str>,
    range: Option<&str>,
    feature: Option<&str>,
    table_id: u32,
    frame: usize,
    to_stop: bool,
) -> Result<Output> {
    let code = GeneticCode::get(table_id).ok_or_else(|| anyhow::anyhow!("unknown genetic code {table_id}"))?;
    let (label, protein) = match (path, seq_text) {
        (_, Some(text)) => {
            let bases = clean_dna(text)?;
            let start = frame.clamp(1, 3) - 1;
            ("sequence".to_string(), translate(&bases[start.min(bases.len())..], code, to_stop))
        }
        (Some(path), None) => {
            let doc = load(&PathBuf::from(path))?;
            let region = resolve_region(&doc, range, feature)?;
            match &region.feature {
                Some(f) if range.is_none() && f.is_coding() => {
                    let mut p = translate_feature(f, &doc);
                    if to_stop && let Some(stop) = p.iter().position(|&a| a == b'*') {
                        p.truncate(stop);
                    }
                    (region.label, p)
                }
                _ => {
                    let mut bases = match &region.feature {
                        Some(f) => f.extract(&doc.seq),
                        None => doc.subseq(region.span),
                    };
                    if region.strand == Strand::Reverse && region.feature.is_none() {
                        bases = alphabet::revcomp(&bases);
                    }
                    let start = frame.clamp(1, 3) - 1;
                    (region.label, translate(&bases[start.min(bases.len())..], code, to_stop))
                }
            }
        }
        (None, None) => bail!("give a file or --seq"),
    };
    let mass = stats::protein_mass(&protein);
    let value = json!({
        "region": label,
        "protein": String::from_utf8_lossy(&protein),
        "length": protein.len(),
        "mass_da": mass.round(),
        "isoelectric_point": (stats::isoelectric_point(&protein) * 100.0).round() / 100.0,
        "extinction_coefficient": stats::extinction_coefficient(&protein),
    });
    Ok(Output::new(value, String::from_utf8_lossy(&protein).into_owned()))
}
