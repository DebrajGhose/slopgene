//! Primer and PCR commands: melting temperatures, PCR simulation, primer design, conversion.

use crate::util::*;
use anyhow::{Result, anyhow, bail};
use serde_json::{Value, json};
use slopgene_core::io::Format;
use slopgene_core::pcr::pcr;
use slopgene_core::primer_design::{self, DesignOptions};
use slopgene_core::primers::{self, BindingRules, TmConditions};
use slopgene_core::render::thousands;
use slopgene_core::span::parse_range;
use slopgene_core::*;
use std::path::{Path, PathBuf};

pub fn tm(seqs: &[String], conditions: TmConditions) -> Result<Output> {
    let mut values = Vec::new();
    let mut rows =
        vec![vec!["sequence".into(), "length".into(), "Tm".into(), "GC".into(), "self any".into(), "self end".into(), "3' ΔG".into(), "warnings".into()]];
    for s in seqs {
        let bases = clean_dna(s)?;
        let report = primers::analyze(&bases, &conditions);
        values.push(json!({
            "sequence": String::from_utf8_lossy(&bases),
            "length": report.length,
            "tm": (report.tm * 10.0).round() / 10.0,
            "gc_percent": (report.gc_percent * 10.0).round() / 10.0,
            "self_any": report.self_any,
            "self_end": report.self_end,
            "end_stability": (report.end_stability * 100.0).round() / 100.0,
            "self_dimer": report.self_dimer,
            "three_prime_dimer": report.three_prime_dimer,
            "hairpin_stem": report.hairpin_stem,
            "longest_run": report.longest_run,
            "warnings": report.warnings,
        }));
        rows.push(vec![
            String::from_utf8_lossy(&bases).into_owned(),
            report.length.to_string(),
            format!("{:.1} °C", report.tm),
            format!("{:.0}%", report.gc_percent),
            format!("{:.1}", report.self_any),
            format!("{:.1}", report.self_end),
            format!("{:.2}", report.end_stability),
            report.warnings.join("; "),
        ]);
    }
    Ok(Output::new(Value::Array(values), table(&rows)))
}

/// A primer given as a name (looked up in the template file) or as a sequence.
fn primer_arg(doc: &Document, text: &str, fallback_name: &str) -> Result<Primer> {
    if let Some(p) = doc.primers.iter().find(|p| p.name.eq_ignore_ascii_case(text)) {
        return Ok(p.clone());
    }
    match clean_dna(text) {
        Ok(bases) => Ok(Primer::new(fallback_name, String::from_utf8_lossy(&bases).into_owned())),
        Err(_) => {
            let names: Vec<&str> = doc.primers.iter().map(|p| p.name.as_str()).collect();
            bail!(
                "'{text}' is neither a primer in this file ({}) nor a DNA sequence",
                if names.is_empty() { "it has none".to_string() } else { names.join(", ") }
            )
        }
    }
}

pub fn pcr_cmd(template: &str, fwd: &str, rev: &str, output: Option<&str>, min_tm: f64, strict: bool) -> Result<Output> {
    let path = PathBuf::from(template);
    let doc = load(&path)?;
    let forward = primer_arg(&doc, fwd, "forward")?;
    let reverse = primer_arg(&doc, rev, "reverse")?;
    let rules = BindingRules { min_tm, allow_mismatch: !strict, ..BindingRules::default() };
    let product = pcr(&doc, &forward, &reverse, &rules)?;
    let seq_len = doc.len();
    let site_json = |b: &primers::Binding| {
        let (s, e) = b.span.to_1based(seq_len);
        json!({ "start": s, "end": e, "strand": b.strand.name(), "annealed": b.annealed, "mismatches": b.mismatches, "tm": (b.tm * 10.0).round() / 10.0 })
    };
    let mut value = json!({
        "template": path_text(&path),
        "product_length": product.doc.len(),
        "template_region": product.template_span.format(seq_len),
        "forward_site": site_json(&product.forward_site),
        "reverse_site": site_json(&product.reverse_site),
        "other_products": product.alternatives,
        "features": product.doc.features.iter().map(|f| f.name.clone()).collect::<Vec<_>>(),
    });
    let mut text = format!(
        "product {} bp from template {} (forward anneals {} nt, reverse {} nt)",
        thousands(product.doc.len()),
        product.template_span.format(seq_len),
        product.forward_site.annealed,
        product.reverse_site.annealed
    );
    if product.alternatives > 0 {
        text.push_str(&format!("\nwarning: {} other primer pairings would also give products", product.alternatives));
    }
    if let Some(out) = output {
        let mut out_doc = product.doc.clone();
        out_doc.meta.created = Some(today_iso());
        save(&out_doc, Path::new(out), None)?;
        value["written"] = json!(out);
        text.push_str(&format!("\nwrote {out}"));
    } else {
        value["sequence"] = json!(product.doc.seq_str());
    }
    Ok(Output::new(value, text))
}

#[allow(clippy::too_many_arguments)]
pub fn design_primers(path: &str, range: &str, target_tm: f64, fwd_tail: &str, rev_tail: &str, add: bool, output: Option<&str>) -> Result<Output> {
    let path = PathBuf::from(path);
    let mut doc = load(&path)?;
    let span = parse_range(range, doc.len(), doc.is_circular())?;
    let tail = |t: &str| -> Result<Vec<u8>> { if t.is_empty() { Ok(Vec::new()) } else { clean_dna(t) } };
    let options = DesignOptions { target_tm, forward_tail: tail(fwd_tail)?, reverse_tail: tail(rev_tail)?, count: 1, ..DesignOptions::default() };
    let pairs = primer_design::design_pairs(&doc, span, &options)?;
    let best = &pairs[0];
    let (fwd, rev) = best.primers(&doc.name, doc.len());
    // Each primer with Primer3's numbers for it.
    let describe = |p: &Primer, c: &primer_design::Candidate| {
        json!({
            "name": p.name, "sequence": p.seq, "length": p.seq.len(), "anneals": c.annealing().len(), "annealing": p.description,
            "tm": (c.tm * 10.0).round() / 10.0, "gc_percent": (c.gc_percent * 10.0).round() / 10.0,
            "self_any": c.self_any, "self_end": c.self_end, "end_stability": (c.end_stability * 100.0).round() / 100.0,
            "penalty": (c.penalty * 1000.0).round() / 1000.0, "warnings": c.warnings,
        })
    };
    let product_len = best.product_len;
    let mut value = json!({
        "forward": describe(&fwd, &best.forward), "reverse": describe(&rev, &best.reverse), "product_length": product_len,
        "pair": { "compl_any": best.compl_any, "compl_end": best.compl_end, "penalty": (best.penalty * 1000.0).round() / 1000.0, "warnings": best.warnings },
        "method": "Primer3",
    });
    let mut text =
        format!("{}  {}  ({})\n{}  {}  ({})\nproduct: {} bp", fwd.name, fwd.seq, fwd.description, rev.name, rev.seq, rev.description, thousands(product_len));
    let warnings: Vec<&String> = best.forward.warnings.iter().chain(&best.reverse.warnings).chain(&best.warnings).collect();
    if !warnings.is_empty() {
        text.push_str(&format!("\nwarning: {}", warnings.iter().map(|w| w.as_str()).collect::<Vec<_>>().join("; ")));
    }
    if add {
        doc.primers.push(fwd);
        doc.primers.push(rev);
        primers::refresh_binding_sites(&mut doc, &BindingRules::default());
        doc.meta.modified = Some(today_iso());
        let target = output.map(PathBuf::from).unwrap_or_else(|| path.clone());
        save(&doc, &target, None)?;
        value["written"] = json!(path_text(&target));
        text.push_str(&format!("\nadded both primers to {}", path_text(&target)));
    }
    Ok(Output::new(value, text))
}

pub fn convert(inputs: &[String], recursive: bool, to: Option<&str>, output: Option<&str>, out_dir: Option<&str>, force: bool) -> Result<Output> {
    let paths = expand_paths(inputs, recursive)?;
    let format = match (to, output) {
        (Some(t), _) => parse_format(t)?,
        (None, Some(o)) => Format::from_path(Path::new(o)).ok_or_else(|| anyhow!("cannot tell the format from '{o}'; add --to"))?,
        (None, None) => bail!("say which format to convert to with --to (dna, genbank, fasta, json)"),
    };
    if output.is_some() && paths.len() > 1 {
        bail!("-o names one output file; with several inputs use --out-dir");
    }
    let mut values = Vec::new();
    let mut lines = Vec::new();
    for path in &paths {
        let docs = load_all(path)?;
        for (ii, doc) in docs.iter().enumerate() {
            let ext = format.extension(doc.molecule);
            let target = match (output, out_dir) {
                (Some(o), _) => PathBuf::from(o),
                (None, dir) => {
                    let stem = if docs.len() > 1 {
                        format!("{}_{}", path.file_stem().and_then(|s| s.to_str()).unwrap_or("record"), ii + 1)
                    } else {
                        path.file_stem().and_then(|s| s.to_str()).unwrap_or("record").to_string()
                    };
                    let base = dir.map(PathBuf::from).unwrap_or_else(|| path.parent().map(Path::to_path_buf).unwrap_or_default());
                    base.join(format!("{stem}.{ext}"))
                }
            };
            if target.exists() && !force {
                if target == *path {
                    bail!("{} would overwrite its own input; choose --out-dir or -o", path_text(&target));
                }
                bail!("{} exists; pass --force to overwrite", path_text(&target));
            }
            if let Some(dir) = target.parent().filter(|d| !d.as_os_str().is_empty()) {
                std::fs::create_dir_all(dir)?;
            }
            save(doc, &target, Some(format))?;
            values.push(json!({ "input": path_text(path), "output": path_text(&target), "format": format.name(), "name": doc.name }));
            lines.push(format!("{} -> {}", path_text(path), path_text(&target)));
        }
    }
    Ok(Output::new(Value::Array(values), lines.join("\n")))
}
