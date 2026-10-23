//! Shared helpers: expanding file arguments, loading and saving documents, resolving regions,
//! and printing results as text or JSON.

use anyhow::{Context, Result, anyhow, bail};
use serde_json::Value;
use slopgene_core::enzymes::EnzymeDb;
use slopgene_core::io::{self, Format};
use slopgene_core::span::parse_range;
use slopgene_core::*;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// Expand file arguments into paths.
///
/// Each argument may be a file, a folder (its sequence files; with `recursive`, also those in
/// subfolders) or a glob pattern such as `plasmids/**/*.dna`. Patterns are expanded here rather
/// than by the shell, so they work the same in PowerShell, cmd and bash.
pub fn expand_paths(args: &[String], recursive: bool) -> Result<Vec<PathBuf>> {
    use slopgene_script::files::{glob, has_glob_chars, sequence_files};
    let mut out = Vec::new();
    for arg in args {
        let path = PathBuf::from(arg);
        if path.is_file() {
            out.push(path);
        } else if path.is_dir() {
            out.extend(sequence_files(arg, recursive).into_iter().map(PathBuf::from));
        } else if has_glob_chars(arg) {
            let found = glob(arg).map_err(|e| anyhow!(e))?;
            if found.is_empty() {
                bail!("no files match '{arg}'");
            }
            out.extend(found.into_iter().map(PathBuf::from));
        } else {
            bail!("'{arg}' does not exist");
        }
    }
    if out.is_empty() {
        bail!("no sequence files found");
    }
    Ok(out)
}

/// Read the first record of a file.
pub fn load(path: &Path) -> Result<Document> {
    io::read_one(path).with_context(|| format!("could not read {}", path.display()))
}

/// Read every record of a file.
pub fn load_all(path: &Path) -> Result<Vec<Document>> {
    io::read_path(path).with_context(|| format!("could not read {}", path.display()))
}

/// Save a document; the format comes from `format` or else the file extension.
pub fn save(doc: &Document, path: &Path, format: Option<Format>) -> Result<Format> {
    io::write_path(doc, path, format).with_context(|| format!("could not write {}", path.display()))
}

/// Refuse to write a document back into its own file when that file's format would drop part of
/// it (a feature added to a FASTA file, a primer added to a GenBank file). Nothing is written,
/// and the message says how to keep everything.
pub fn refuse_lossy_in_place(doc: &Document, path: &Path) -> Result<()> {
    let format = io::format_of(path);
    let lost = io::losses(doc, format);
    if lost.is_empty() {
        return Ok(());
    }
    let keep = path.with_extension(Format::SnapGene.extension(doc.molecule));
    bail!(
        "{} is a {} file, which cannot keep {}; nothing was written. Add -o {} to write a SnapGene file that keeps everything",
        path_text(path),
        format.name(),
        lost.join(", "),
        path_text(&keep)
    )
}

/// Parse a format name given on the command line.
pub fn parse_format(name: &str) -> Result<Format> {
    Format::from_name(name).ok_or_else(|| anyhow!("unknown format '{name}' (use dna, genbank, fasta or json)"))
}

/// The enzyme database: `SLOPGENE_ENZYMES` (a TSV file) when set, otherwise the built-in one.
pub fn enzyme_db() -> &'static EnzymeDb {
    static DB: OnceLock<EnzymeDb> = OnceLock::new();
    DB.get_or_init(|| {
        if let Ok(path) = std::env::var("SLOPGENE_ENZYMES") {
            match std::fs::read_to_string(&path).map_err(|e| e.to_string()).and_then(|t| EnzymeDb::from_tsv(&t).map_err(|e| e.to_string())) {
                Ok(db) => return db,
                Err(e) => eprintln!("warning: SLOPGENE_ENZYMES={path} could not be loaded ({e}); using the built-in enzymes"),
            }
        }
        EnzymeDb::builtin().clone()
    })
}

/// A region of a document chosen by `--range` or `--feature`, or the whole sequence.
pub struct Region {
    pub span: Span,
    /// Reverse for reverse-strand features (their sequence reads along the bottom strand).
    pub strand: Strand,
    pub label: String,
    /// The feature, when the region was chosen by feature name.
    pub feature: Option<Feature>,
}

pub fn resolve_region(doc: &Document, range: Option<&str>, feature: Option<&str>) -> Result<Region> {
    let seq_len = doc.len();
    if let Some(name) = feature {
        let found = find_feature(doc, name)?;
        let span = found.span(seq_len);
        return Ok(Region { span, strand: found.strand, label: format!("{} ({})", found.name, span.format(seq_len)), feature: Some(found.clone()) });
    }
    if let Some(text) = range {
        let span = parse_range(text, seq_len, doc.is_circular())?;
        return Ok(Region { span, strand: Strand::Forward, label: span.format(seq_len), feature: None });
    }
    Ok(Region { span: Span::new(0, seq_len), strand: Strand::Forward, label: format!("1..{seq_len}"), feature: None })
}

/// Find a feature by exact name (case-insensitive), or by `#N` index (1-based).
pub fn find_feature<'a>(doc: &'a Document, name: &str) -> Result<&'a Feature> {
    if let Some(index) = name.strip_prefix('#').and_then(|n| n.parse::<usize>().ok()) {
        return doc.features.get(index.wrapping_sub(1)).ok_or_else(|| anyhow!("there is no feature #{index} ({} features)", doc.features.len()));
    }
    let matches: Vec<&Feature> = doc.features.iter().filter(|f| f.name.eq_ignore_ascii_case(name)).collect();
    match matches.len() {
        1 => Ok(matches[0]),
        0 => {
            let names: Vec<&str> = doc.features.iter().map(|f| f.name.as_str()).take(30).collect();
            bail!("no feature named '{name}'. Features: {}", names.join(", "))
        }
        n => bail!("{n} features are named '{name}'; pick one with #N (see `slopgene features`)"),
    }
}

/// Today's date as ISO `YYYY-MM-DD` (UTC), for "last modified" stamps.
pub fn today_iso() -> String {
    slopgene_core::clock::today_iso()
}

/// How results are printed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Text,
    Json,
    JsonLines,
}

/// A command's result: a JSON value for machines and text for people.
pub struct Output {
    pub json: Value,
    pub text: String,
}

impl Output {
    pub fn new(json: Value, text: impl Into<String>) -> Output {
        Output { json, text: text.into() }
    }

    pub fn print(&self, mode: Mode) {
        let line = slopgene_script::stdout_line;
        match mode {
            Mode::Text => {
                if !self.text.is_empty() {
                    line(self.text.trim_end());
                }
            }
            Mode::Json => line(&serde_json::to_string_pretty(&self.json).unwrap_or_default()),
            Mode::JsonLines => match &self.json {
                Value::Array(items) => {
                    for item in items {
                        line(&serde_json::to_string(item).unwrap_or_default());
                    }
                }
                other => line(&serde_json::to_string(other).unwrap_or_default()),
            },
        }
    }
}

/// A fixed-width text table: the first row is the header.
pub fn table(rows: &[Vec<String>]) -> String {
    if rows.is_empty() {
        return String::new();
    }
    let columns = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    let mut widths = vec![0usize; columns];
    for row in rows {
        for (ii, cell) in row.iter().enumerate() {
            widths[ii] = widths[ii].max(cell.chars().count());
        }
    }
    let mut out = String::new();
    for (ii, row) in rows.iter().enumerate() {
        let mut line = String::new();
        for (iii, cell) in row.iter().enumerate() {
            if iii + 1 == row.len() {
                line.push_str(cell);
            } else {
                line.push_str(&format!("{cell:<width$}  ", width = widths[iii]));
            }
        }
        out.push_str(line.trim_end());
        out.push('\n');
        if ii == 0 && rows.len() > 1 {
            let rule: Vec<String> = widths.iter().map(|w| "-".repeat(*w)).collect();
            out.push_str(rule.join("  ").trim_end());
            out.push('\n');
        }
    }
    out
}

/// Clean a nucleotide sequence typed on the command line.
pub fn clean_dna(text: &str) -> Result<Vec<u8>> {
    let cleaned = slopgene_core::alphabet::clean_nucleotides(text).map_err(|e| anyhow!(e))?;
    if cleaned.is_empty() {
        bail!("the sequence is empty");
    }
    Ok(cleaned.to_ascii_uppercase())
}

/// Path as text with forward slashes, for JSON output that reads the same on every platform.
pub fn path_text(path: &Path) -> String {
    path.to_string_lossy().replace('\\', "/")
}
