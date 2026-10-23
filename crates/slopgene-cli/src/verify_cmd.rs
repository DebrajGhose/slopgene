//! `slopgene align` and `slopgene verify`: compare sequencing reads with reference maps.
//!
//! Alignment ignores where a circular read starts and which strand it is on, so a correct clone
//! comes back at 100% whatever its rotation. Changes of 50 bp or more (deletions, insertions,
//! duplications, inversions, replacements) are reported as structural changes, and identity is
//! measured outside them, so a clone that lost one element still reads as that plasmid.

use crate::util::*;
use anyhow::{Result, bail};
use rayon::prelude::*;
use serde_json::{Value, json};
use slopgene_core::align::{Alignment, KmerIndex, align_read_with, rank_references};
use slopgene_core::render::thousands;
use slopgene_core::*;
use std::path::{Path, PathBuf};

/// The structural changes as indented lines.
fn events_text(aln: &Alignment) -> String {
    aln.events.iter().map(|e| format!("structural change: {}\n", e.description)).collect()
}

/// The differences as a table, at most `limit` rows (the JSON output always has them all).
fn differences_text(aln: &Alignment, limit: usize) -> String {
    let mut rows = vec![vec!["position".into(), "change".into(), "features".into(), "effect".into(), "note".into()]];
    for d in aln.differences.iter().take(limit) {
        let change = match d.kind {
            align::DiffKind::Substitution => format!("{}>{}", d.reference, d.read),
            align::DiffKind::Insertion => format!("+{}", d.read),
            align::DiffKind::Deletion => format!("-{}", d.reference),
        };
        rows.push(vec![
            thousands(d.position + 1),
            change,
            if d.features.is_empty() { "-".into() } else { d.features.join(", ") },
            d.effects.join("; "),
            if d.homopolymer { "homopolymer".into() } else { String::new() },
        ]);
    }
    if rows.len() == 1 {
        return String::new();
    }
    let mut text = table(&rows);
    if aln.differences.len() > limit {
        text.push_str(&format!("... and {} more (see --json)\n", aln.differences.len() - limit));
    }
    text
}

/// Align one read to one reference. With `show` (columns per line), the output also holds the
/// aligned bases: written out in the text, and as two gapped rows in the JSON. With `map`, the
/// reference map is drawn to that file with the read as a ring and its differences marked.
pub fn align(reference: &str, read: &str, show: Option<usize>, map: Option<&str>) -> Result<Output> {
    let ref_doc = load(Path::new(reference))?;
    let read_doc = load(Path::new(read))?;
    let index = KmerIndex::new(&ref_doc.seq, ref_doc.is_circular());
    let Some(aln) = align_read_with(&ref_doc, &index, &read_doc.seq, read_doc.is_circular()) else {
        bail!("{} does not align to {} (no shared 15-mers)", read_doc.name, ref_doc.name);
    };
    let mut value = aln.to_json();
    value["reference"] = json!(ref_doc.name);
    value["read"] = json!(read_doc.name);
    let mut text = format!("{} vs {}: {}\n", read_doc.name, ref_doc.name, aln.summary());
    text.push_str(&format!("method: {}\n", align::method_line()));
    text.push_str(&events_text(&aln));
    text.push_str(&differences_text(&aln, 60));
    if let Some(path) = map {
        let options = render::MapOptions { alignments: vec![aln.track(&read_doc.name, true)], ..render::MapOptions::default() };
        let scene = render::render_map(&ref_doc, &options, &render::ApproxMeasure);
        crate::render_cmd::save_scene(&scene, Path::new(path), 2.0)?;
        value["map"] = json!(path_text(Path::new(path)));
        text.push_str(&format!("wrote {}\n", path_text(Path::new(path))));
    }
    if let Some(width) = show {
        let (upper, lower) = aln.rows();
        value["aligned"] = json!({ "reference": upper, "read": lower });
        text.push('\n');
        text.push_str(&aln.pairwise_text(width, &ref_doc.name, &read_doc.name));
    }
    Ok(Output::new(value, text))
}

/// For every read, find the best reference, align, and report.
pub fn verify(reads: &[String], refs: &[String], recursive: bool, top: usize) -> Result<Output> {
    let read_paths = expand_paths(reads, recursive)?;
    let ref_paths = expand_paths(refs, recursive)?;
    let references: Vec<(PathBuf, Document)> =
        ref_paths.par_iter().filter_map(|p| load(p).ok().map(|d| (p.clone(), d))).filter(|(_, d)| d.molecule == Molecule::Dna && !d.is_empty()).collect();
    if references.is_empty() {
        bail!("no usable reference sequences");
    }
    let indexes: Vec<(String, KmerIndex)> = references.par_iter().map(|(_, d)| (d.name.clone(), KmerIndex::new(&d.seq, d.is_circular()))).collect();

    let results: Vec<(PathBuf, Result<Value>, String)> = read_paths
        .par_iter()
        .map(|path| {
            let outcome = (|| -> Result<(Value, String)> {
                let read = load(path)?;
                let ranked = rank_references(&read.seq, &indexes);
                let (best, score) = ranked[0];
                let (ref_path, reference) = &references[best];
                // A read that shares no 15-mer with even the best reference is no match at all.
                let Some(aln) = align_read_with(reference, &indexes[best].1, &read.seq, read.is_circular()) else {
                    let value = json!({ "read": read.name, "read_path": path_text(path), "reference": reference.name,
                        "reference_path": path_text(ref_path), "verdict": "not this plasmid", "identity": 0.0, "events": [], "differences": [] });
                    return Ok((value, format!("{}\n  best: {}: not this plasmid (no shared sequence)\n", read.name, reference.name)));
                };
                let others: Vec<Value> = ranked
                    .iter()
                    .skip(1)
                    .take(top.saturating_sub(1))
                    .map(|&(ii, s)| json!({ "reference": references[ii].1.name, "kmer_containment": (s * 10000.0).round() / 100.0 }))
                    .collect();
                let mut value = aln.to_json();
                value["read"] = json!(read.name);
                value["read_path"] = json!(path_text(path));
                value["reference"] = json!(reference.name);
                value["reference_path"] = json!(path_text(ref_path));
                value["kmer_containment"] = json!((score * 10000.0).round() / 100.0);
                value["runners_up"] = Value::Array(others);
                let mut text = format!("{}\n  best: {}: {}\n", read.name, reference.name, aln.summary());
                // A read that is not this plasmid differs everywhere; listing that is noise.
                if aln.verdict() != "not this plasmid" {
                    let detail = events_text(&aln) + differences_text(&aln, 40).as_str();
                    for line in detail.lines() {
                        text.push_str("    ");
                        text.push_str(line);
                        text.push('\n');
                    }
                }
                Ok((value, text))
            })();
            match outcome {
                Ok((v, t)) => (path.clone(), Ok(v), t),
                Err(e) => (path.clone(), Err(e), String::new()),
            }
        })
        .collect();

    let mut values = Vec::new();
    let mut text = String::new();
    for (path, result, t) in results {
        match result {
            Ok(v) => {
                values.push(v);
                text.push_str(&t);
            }
            Err(e) => {
                values.push(json!({ "read_path": path_text(&path), "error": format!("{e:#}") }));
                text.push_str(&format!("{}: error: {e:#}\n", path_text(&path)));
            }
        }
    }
    Ok(Output::new(Value::Array(values), text))
}
