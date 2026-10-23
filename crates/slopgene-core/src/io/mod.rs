//! Reading and writing sequence files.
//!
//! Formats are detected from content first and file extension second, so a GenBank file named
//! `.txt` still opens. Writes are atomic: the new content goes to a temporary file in the same
//! folder, which then replaces the target, so a crash or a concurrent reader never sees a
//! half-written file.

pub mod fasta;
pub mod genbank;
pub mod snapgene;

use crate::error::{Error, Result};
use crate::model::{Document, Ends, Molecule};
use std::path::Path;

/// A file format SlopGene can read and write.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Format {
    /// SnapGene `.dna` / `.rna` / `.prot` (and `.ftrs` feature libraries).
    SnapGene,
    GenBank,
    Fasta,
    /// SlopGene's JSON document schema (see [`crate::json`]).
    Json,
}

impl Format {
    pub const ALL: [Format; 4] = [Format::SnapGene, Format::GenBank, Format::Fasta, Format::Json];

    pub fn name(self) -> &'static str {
        match self {
            Format::SnapGene => "snapgene",
            Format::GenBank => "genbank",
            Format::Fasta => "fasta",
            Format::Json => "json",
        }
    }

    /// Parse a format name or extension: `dna`, `snapgene`, `gb`, `gbk`, `genbank`, `fasta`,
    /// `fa`, `json`, ...
    pub fn from_name(name: &str) -> Option<Format> {
        match name.trim().trim_start_matches('.').to_ascii_lowercase().as_str() {
            "snapgene" | "dna" | "rna" | "prot" | "ftrs" | "sg" => Some(Format::SnapGene),
            "genbank" | "gb" | "gbk" | "gbff" | "genbank-file" | "ape" | "str" => Some(Format::GenBank),
            "fasta" | "fa" | "fas" | "fna" | "faa" | "ffn" | "seq" | "fsa" => Some(Format::Fasta),
            "json" => Some(Format::Json),
            _ => None,
        }
    }

    pub fn from_path(path: &Path) -> Option<Format> {
        path.extension().and_then(|e| e.to_str()).and_then(Format::from_name)
    }

    /// Usual extension for this format and molecule, without the dot.
    pub fn extension(self, molecule: Molecule) -> &'static str {
        match (self, molecule) {
            (Format::SnapGene, Molecule::Dna) => "dna",
            (Format::SnapGene, Molecule::Rna) => "rna",
            (Format::SnapGene, Molecule::Protein) => "prot",
            (Format::GenBank, _) => "gb",
            (Format::Fasta, _) => "fasta",
            (Format::Json, _) => "json",
        }
    }
}

/// Work out the format from the bytes, falling back to the file extension.
pub fn detect(bytes: &[u8], path: Option<&Path>) -> Option<Format> {
    if snapgene::is_snapgene(bytes) {
        return Some(Format::SnapGene);
    }
    let head = String::from_utf8_lossy(&bytes[..bytes.len().min(4096)]);
    let start = head.trim_start_matches('\u{feff}').trim_start();
    if start.starts_with("LOCUS") {
        return Some(Format::GenBank);
    }
    if start.starts_with('>') || start.starts_with(';') {
        return Some(Format::Fasta);
    }
    if start.starts_with('{') {
        return Some(Format::Json);
    }
    path.and_then(Format::from_path)
}

/// Read every record in `bytes`.
pub fn read_bytes(bytes: &[u8], format: Option<Format>) -> Result<Vec<Document>> {
    let format = format.or_else(|| detect(bytes, None)).ok_or_else(|| Error::Unsupported("could not recognise the file format".into()))?;
    let text = || std::str::from_utf8(bytes).map(|s| s.to_string()).or_else(|_| Ok::<String, Error>(String::from_utf8_lossy(bytes).into_owned()));
    match format {
        Format::SnapGene => Ok(vec![snapgene::read(bytes)?]),
        Format::GenBank => genbank::read(&text()?),
        Format::Fasta => fasta::read(&text()?),
        Format::Json => {
            let value: serde_json::Value = serde_json::from_slice(bytes).map_err(|e| Error::Parse(format!("JSON: {e}")))?;
            crate::json::documents_from_value(&value)
        }
    }
}

/// Read every record in a file. Records from single-record files are named after the file,
/// as SnapGene does; records from multi-record files keep their own names.
pub fn read_path(path: impl AsRef<Path>) -> Result<Vec<Document>> {
    let path = path.as_ref();
    let bytes = std::fs::read(path).map_err(|e| Error::Io(std::io::Error::new(e.kind(), format!("{}: {e}", path.display()))))?;
    let format = detect(&bytes, Some(path)).ok_or_else(|| Error::Unsupported(format!("{}: could not recognise the file format", path.display())))?;
    let mut docs = read_bytes(&bytes, Some(format)).map_err(|e| match e {
        Error::Parse(msg) => Error::Parse(format!("{}: {msg}", path.display())),
        other => other,
    })?;
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("").to_string();
    if docs.len() == 1 {
        if format != Format::Json || docs[0].name.is_empty() {
            docs[0].name = stem;
        }
    } else {
        for (ii, doc) in docs.iter_mut().enumerate() {
            if doc.name.is_empty() {
                doc.name = format!("{stem}_{}", ii + 1);
            }
        }
    }
    Ok(docs)
}

/// Read the first record of a file.
pub fn read_one(path: impl AsRef<Path>) -> Result<Document> {
    let path = path.as_ref();
    read_path(path)?.into_iter().next().ok_or_else(|| Error::Parse(format!("{}: file contains no records", path.display())))
}

/// Encode a document in the given format.
pub fn to_bytes(doc: &Document, format: Format) -> Result<Vec<u8>> {
    Ok(match format {
        Format::SnapGene => snapgene::write(doc)?,
        Format::GenBank => genbank::write(doc).into_bytes(),
        Format::Fasta => fasta::write(doc, 70).into_bytes(),
        Format::Json => {
            let value = crate::json::document_to_value(doc, true);
            let mut text = serde_json::to_string_pretty(&value).map_err(|e| Error::Parse(e.to_string()))?;
            text.push('\n');
            text.into_bytes()
        }
    })
}

/// Write a document to `path`. The format is `format`, or else guessed from the extension
/// (SnapGene when there is none). Returns the format used.
/// What `doc` holds that `format` cannot store, in words for a message ("features", "primers",
/// ...). Empty when nothing would be lost. SnapGene files and SlopGene JSON keep everything
/// SlopGene models, so they are the formats to save annotated work in.
pub fn losses(doc: &Document, format: Format) -> Vec<&'static str> {
    let mut lost = Vec::new();
    match format {
        Format::SnapGene | Format::Json => {}
        Format::Fasta => {
            if !doc.features.is_empty() {
                lost.push("features");
            }
            if !doc.primers.is_empty() {
                lost.push("primers");
            }
            if doc.is_circular() {
                lost.push("the circular topology");
            }
        }
        Format::GenBank => {
            if !doc.primers.is_empty() {
                lost.push("primers (GenBank keeps only their binding sites, as features)");
            }
            if doc.ends != Ends::default() {
                lost.push("sticky ends");
            }
        }
    }
    lost
}

/// The format a path's extension names, SnapGene when it names none.
pub fn format_of(path: &Path) -> Format {
    Format::from_path(path).unwrap_or(Format::SnapGene)
}

pub fn write_path(doc: &Document, path: impl AsRef<Path>, format: Option<Format>) -> Result<Format> {
    let path = path.as_ref();
    let format = format.or_else(|| Format::from_path(path)).unwrap_or(Format::SnapGene);
    let bytes = to_bytes(doc, format)?;
    write_atomic(path, &bytes)?;
    Ok(format)
}

/// Replace `path` with `bytes` atomically (temporary file in the same folder, then rename).
pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    let dir = path.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
    let file_name = path.file_name().and_then(|n| n.to_str()).unwrap_or("out");
    let tmp = dir.join(format!(".{file_name}.slopgene-tmp-{}", std::process::id()));
    std::fs::write(&tmp, bytes)?;
    if let Err(e) = std::fs::rename(&tmp, path) {
        let _ = std::fs::remove_file(&tmp);
        return Err(Error::Io(e));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Feature, Primer, Strand, Topology};
    use crate::span::Span;

    #[test]
    fn losses_name_what_a_format_cannot_keep() {
        let mut doc = Document::new_dna("p", b"ACGTACGTACGTACGTACGT", Topology::Circular);
        assert!(losses(&doc, Format::Fasta).contains(&"the circular topology"));
        doc.features.push(Feature::new("f", "misc_feature", Span::new(2, 8), Strand::Forward));
        doc.primers.push(Primer::new("p1", "ACGTACGTAC"));
        assert_eq!(losses(&doc, Format::Fasta), vec!["features", "primers", "the circular topology"]);
        assert_eq!(losses(&doc, Format::GenBank).len(), 1, "GenBank keeps features, not primers");
        assert!(losses(&doc, Format::SnapGene).is_empty());
        assert!(losses(&doc, Format::Json).is_empty());
        assert_eq!(format_of(Path::new("a/b.gbk")), Format::GenBank);
        assert_eq!(format_of(Path::new("a/b")), Format::SnapGene);
    }
}
