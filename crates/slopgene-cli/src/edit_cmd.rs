//! `slopgene edit`: change a file in place (or write the result elsewhere).
//!
//! Every edit keeps annotations in step with the sequence (see `slopgene_core::edit`), stamps
//! the document as modified today, and is written atomically. `--dry-run` reports what would
//! change without writing.

use crate::util::*;
use anyhow::{Result, anyhow, bail};
use serde_json::json;
use slopgene_core::span::parse_range;
use slopgene_core::*;
use std::path::PathBuf;

/// One edit operation, parsed from the command line.
pub enum EditOp {
    Insert { at: usize, seq: String },
    Delete { range: String },
    Replace { range: String, seq: String },
    Revcomp,
    Invert { range: String },
    Origin { at: usize },
    Topology { topology: String },
    Methylation { state: String },
    Describe { text: String },
    AddFeature { name: String, kind: String, range: String, strand: String, color: Option<String>, qualifiers: Vec<String> },
    RemoveFeature { name: String, all: bool },
    UpdateFeature { name: String, new_name: Option<String>, kind: Option<String>, range: Option<String>, strand: Option<String>, color: Option<String> },
    AddPrimer { name: String, seq: String, description: Option<String> },
    RemovePrimer { name: String },
    RefreshPrimers,
}

impl EditOp {
    fn label(&self) -> &'static str {
        match self {
            EditOp::Insert { .. } => "insert",
            EditOp::Delete { .. } => "delete",
            EditOp::Replace { .. } => "replace",
            EditOp::Revcomp => "revcomp",
            EditOp::Invert { .. } => "invert",
            EditOp::Origin { .. } => "origin",
            EditOp::Topology { .. } => "topology",
            EditOp::Methylation { .. } => "methylation",
            EditOp::Describe { .. } => "describe",
            EditOp::AddFeature { .. } => "add-feature",
            EditOp::RemoveFeature { .. } => "remove-feature",
            EditOp::UpdateFeature { .. } => "update-feature",
            EditOp::AddPrimer { .. } => "add-primer",
            EditOp::RemovePrimer { .. } => "remove-primer",
            EditOp::RefreshPrimers => "refresh-primers",
        }
    }
}

fn parse_strand(text: &str) -> Result<Strand> {
    Strand::parse(text).ok_or_else(|| anyhow!("strand must be forward, reverse, both or none (got '{text}')"))
}

fn parse_color(text: &str) -> Result<Color> {
    Color::parse_hex(text).ok_or_else(|| anyhow!("colour must be #rrggbb (got '{text}')"))
}

/// Indices of features matching `name` (exact, case-insensitive), or `#N`.
fn feature_indices(doc: &Document, name: &str) -> Result<Vec<usize>> {
    if let Some(n) = name.strip_prefix('#').and_then(|n| n.parse::<usize>().ok()) {
        if n == 0 || n > doc.features.len() {
            bail!("there is no feature #{n} ({} features)", doc.features.len());
        }
        return Ok(vec![n - 1]);
    }
    let found: Vec<usize> = doc.features.iter().enumerate().filter(|(_, f)| f.name.eq_ignore_ascii_case(name)).map(|(ii, _)| ii).collect();
    if found.is_empty() {
        bail!("no feature named '{name}'");
    }
    Ok(found)
}

/// Apply one operation. Returns a short description of what happened.
pub fn apply(doc: &mut Document, op: &EditOp) -> Result<String> {
    let seq_len = doc.len();
    let circular = doc.is_circular();
    let rules = primers::BindingRules::default();
    let message = match op {
        EditOp::Insert { at, seq } => {
            let bases = clean_dna(seq)?;
            if *at == 0 || *at > seq_len + 1 {
                bail!("--at must be between 1 and {} (the new bases start at that position)", seq_len + 1);
            }
            edit::insert(doc, at - 1, &bases)?;
            primers::refresh_binding_sites(doc, &rules);
            format!("inserted {} bases at {at}", bases.len())
        }
        EditOp::Delete { range } => {
            let span = parse_range(range, seq_len, circular)?;
            edit::delete(doc, span)?;
            primers::refresh_binding_sites(doc, &rules);
            format!("deleted {} bases ({})", span.len(), span.format(seq_len))
        }
        EditOp::Replace { range, seq } => {
            let span = parse_range(range, seq_len, circular)?;
            let bases = clean_dna(seq)?;
            edit::replace(doc, span, &bases)?;
            primers::refresh_binding_sites(doc, &rules);
            format!("replaced {} ({} bases) with {} bases", span.format(seq_len), span.len(), bases.len())
        }
        EditOp::Revcomp => {
            edit::reverse_complement(doc);
            "reverse-complemented the whole sequence".to_string()
        }
        EditOp::Invert { range } => {
            let span = parse_range(range, seq_len, circular)?;
            edit::invert_region(doc, span)?;
            primers::refresh_binding_sites(doc, &rules);
            format!("inverted {}", span.format(seq_len))
        }
        EditOp::Origin { at } => {
            if *at == 0 || *at > seq_len {
                bail!("--at must be between 1 and {seq_len}");
            }
            edit::set_origin(doc, at - 1)?;
            format!("base {at} is now base 1")
        }
        EditOp::Topology { topology } => {
            let t = Topology::parse(topology).ok_or_else(|| anyhow!("topology must be linear or circular"))?;
            edit::set_topology(doc, t);
            format!("topology set to {}", t.name())
        }
        EditOp::Methylation { state } => {
            let methylation = doc.methylation.changed(state).map_err(|e| anyhow!(e))?;
            doc.methylation = methylation;
            format!("methylation set to {}", methylation.describe())
        }
        EditOp::Describe { text } => {
            doc.meta.description = text.clone();
            "description updated".to_string()
        }
        EditOp::AddFeature { name, kind, range, strand, color, qualifiers } => {
            let span = parse_range(range, seq_len, circular)?;
            let mut feature = Feature::new(name.clone(), kind.clone(), span, parse_strand(strand)?);
            if let Some(c) = color {
                feature.set_color(Some(parse_color(c)?));
            }
            if kind.eq_ignore_ascii_case("CDS") {
                for s in &mut feature.segments {
                    s.translated = Some(true);
                }
            }
            for q in qualifiers {
                let (key, value) = q.split_once('=').ok_or_else(|| anyhow!("qualifier must be key=value (got '{q}')"))?;
                let value = match value.trim().parse::<i64>() {
                    Ok(n) if ["codon_start", "transl_table"].contains(&key.trim()) => QualValue::Int(n),
                    _ => QualValue::Text(value.to_string()),
                };
                feature.qualifiers.push(Qualifier { key: key.trim().to_string(), values: vec![value] });
            }
            doc.features.push(feature);
            format!("added feature '{name}' ({kind}, {})", span.format(seq_len))
        }
        EditOp::RemoveFeature { name, all } => {
            let indices = feature_indices(doc, name)?;
            if indices.len() > 1 && !all {
                bail!("{} features are named '{name}'; pass --all to remove them all, or use #N", indices.len());
            }
            for &index in indices.iter().rev() {
                doc.features.remove(index);
            }
            format!("removed {} feature(s) named '{name}'", indices.len())
        }
        EditOp::UpdateFeature { name, new_name, kind, range, strand, color } => {
            let indices = feature_indices(doc, name)?;
            if indices.len() > 1 {
                bail!("{} features are named '{name}'; use #N to pick one", indices.len());
            }
            let feature = &mut doc.features[indices[0]];
            if let Some(n) = new_name {
                feature.name = n.clone();
            }
            if let Some(k) = kind {
                feature.kind = k.clone();
            }
            if let Some(r) = range {
                let span = parse_range(r, seq_len, circular)?;
                let color = feature.color();
                feature.segments = vec![Segment::new(span)];
                feature.set_color(color);
            }
            if let Some(s) = strand {
                feature.strand = parse_strand(s)?;
            }
            if let Some(c) = color {
                feature.set_color(Some(parse_color(c)?));
            }
            format!("updated feature '{name}'")
        }
        EditOp::AddPrimer { name, seq, description } => {
            let bases = clean_dna(seq)?;
            let mut primer = Primer::new(name.clone(), String::from_utf8_lossy(&bases).into_owned());
            primer.description = description.clone().unwrap_or_default();
            doc.primers.push(primer);
            primers::refresh_binding_sites(doc, &rules);
            let sites = doc.primers.last().map(|p| p.sites.len()).unwrap_or(0);
            format!("added primer '{name}' ({sites} binding site(s))")
        }
        EditOp::RemovePrimer { name } => {
            let before = doc.primers.len();
            doc.primers.retain(|p| !p.name.eq_ignore_ascii_case(name));
            if doc.primers.len() == before {
                bail!("no primer named '{name}'");
            }
            format!("removed primer '{name}'")
        }
        EditOp::RefreshPrimers => {
            primers::refresh_binding_sites(doc, &rules);
            let sites: usize = doc.primers.iter().map(|p| p.sites.len()).sum();
            format!("recomputed binding sites ({sites} in total)")
        }
    };
    Ok(message)
}

pub fn run(path: &str, op: EditOp, output: Option<&str>, dry_run: bool) -> Result<Output> {
    let path = PathBuf::from(path);
    let mut doc = load(&path)?;
    let before_len = doc.len();
    let before_features = doc.features.len();
    let message = apply(&mut doc, &op)?;
    doc.meta.modified = Some(today_iso());
    let target = output.map(PathBuf::from).unwrap_or_else(|| path.clone());
    if output.is_none() {
        refuse_lossy_in_place(&doc, &path)?;
    }
    if !dry_run {
        save(&doc, &target, None)?;
    }
    let value = json!({
        "ok": true,
        "operation": op.label(),
        "message": message,
        "file": path_text(&path),
        "written": if dry_run { serde_json::Value::Null } else { json!(path_text(&target)) },
        "length_before": before_len,
        "length_after": doc.len(),
        "features_before": before_features,
        "features_after": doc.features.len(),
        "primers": doc.primers.len(),
    });
    let text = if dry_run { format!("{message} (dry run, nothing written)") } else { format!("{message}; saved {}", path_text(&target)) };
    Ok(Output::new(value, text))
}
