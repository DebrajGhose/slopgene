//! FASTA files: one or more `>name description` headers, each followed by sequence lines.

use crate::error::{Result, parse_err};
use crate::model::*;

/// Read every record. Topology is not stored in FASTA, so records come back linear unless the
/// header says "circular" (a convention several plasmid tools use).
pub fn read(text: &str) -> Result<Vec<Document>> {
    let mut docs: Vec<Document> = Vec::new();
    for line in text.lines() {
        let line = line.trim_end();
        if let Some(header) = line.strip_prefix('>') {
            let header = header.trim();
            let (name, description) = header.split_once(char::is_whitespace).unwrap_or((header, ""));
            let mut doc = Document::new_dna(name, b"", Topology::Linear);
            doc.meta.description = description.trim().to_string();
            if header.to_ascii_lowercase().contains("circular") {
                doc.topology = Topology::Circular;
            }
            docs.push(doc);
        } else if line.starts_with(';') {
            // Old-style comment line.
            continue;
        } else if let Some(doc) = docs.last_mut() {
            doc.seq.extend(line.bytes().filter(|b| b.is_ascii_alphabetic() || *b == b'*' || *b == b'-'));
        } else if !line.trim().is_empty() {
            // Sequence before any header: treat the whole file as one unnamed record.
            let mut doc = Document::new_dna("", b"", Topology::Linear);
            doc.seq.extend(line.bytes().filter(|b| b.is_ascii_alphabetic()));
            docs.push(doc);
        }
    }
    if docs.is_empty() {
        return Err(parse_err("no FASTA records found"));
    }
    for doc in &mut docs {
        if looks_like_protein(&doc.seq) {
            doc.molecule = Molecule::Protein;
            doc.strandedness = Strandedness::Single;
        }
    }
    Ok(docs)
}

/// Guess protein from the letters used. Nucleic-acid sequences are almost entirely A, C, G,
/// T/U and N; ambiguity codes such as M or K are rare there but common amino acids, so they do
/// not count as nucleotides here.
pub fn looks_like_protein(seq: &[u8]) -> bool {
    let letters = seq.iter().filter(|b| b.is_ascii_alphabetic()).count();
    if letters == 0 {
        return false;
    }
    let nucleotide = seq.iter().filter(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U' | b'N')).count();
    (nucleotide as f64) < 0.9 * letters as f64
}

/// Write one record, wrapping the sequence at `width` characters.
pub fn write(doc: &Document, width: usize) -> String {
    let mut out = String::new();
    out.push('>');
    out.push_str(if doc.name.is_empty() { "unnamed" } else { &doc.name });
    if !doc.meta.description.is_empty() {
        out.push(' ');
        out.push_str(&doc.meta.description.replace('\n', " "));
    }
    out.push('\n');
    let width = width.max(1);
    for chunk in doc.seq.chunks(width) {
        out.push_str(std::str::from_utf8(chunk).unwrap_or(""));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multi_record_round_trip() {
        let text = ">a first one\nACGT\nAC\n>b\nMKV*\n";
        let docs = read(text).unwrap();
        assert_eq!(docs.len(), 2);
        assert_eq!(docs[0].name, "a");
        assert_eq!(docs[0].meta.description, "first one");
        assert_eq!(docs[0].seq, b"ACGTAC".to_vec());
        assert_eq!(docs[1].molecule, Molecule::Protein);
        let written = write(&docs[0], 4);
        assert_eq!(written, ">a first one\nACGT\nAC\n");
    }
}
