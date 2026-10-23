//! SnapGene opens what SlopGene writes.
//!
//! SlopGene writes three kinds of `.dna` file: a SnapGene file it edited, a document it made
//! from nothing, and a GenBank record it converted. SnapGene's own command-line converter
//! (`--convert "GenBank - Standard"`) then reads each one and writes it out as GenBank, and the
//! test checks that SnapGene saw the same sequence, topology and features SlopGene wrote.
//!
//! The test launches SnapGene, so it runs only when `SLOPGENE_SNAPGENE_EXE` names the
//! program (for example `C:\Program Files\SnapGene\SnapGene.exe`). The edited file starts from
//! SnapGene's sample project, found as in `snapgene_parity.rs`. Elsewhere the test passes
//! without checking anything.

use slopgene_core::*;
use std::path::{Path, PathBuf};
use std::process::Command;

fn snapgene() -> Option<PathBuf> {
    std::env::var("SLOPGENE_SNAPGENE_EXE").ok().map(PathBuf::from).filter(|p| p.exists())
}

fn samples() -> Option<PathBuf> {
    let dir = std::env::var("SLOPGENE_SNAPGENE_SAMPLES")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("C:\\Program Files\\SnapGene\\sampleData\\Sample project"));
    dir.exists().then_some(dir)
}

/// Have SnapGene convert a `.dna` file to GenBank, and read what it wrote.
fn read_back_through_snapgene(exe: &Path, dna: &Path) -> Document {
    let gb = dna.with_extension("snapgene.gb");
    let _ = std::fs::remove_file(&gb);
    let status = Command::new(exe).args(["--convert", "GenBank - Standard", "--input"]).arg(dna).arg("--output").arg(&gb).status().expect("SnapGene starts");
    assert!(status.success(), "SnapGene failed on {}", dna.display());
    assert!(gb.exists(), "SnapGene wrote nothing for {}: it could not open the file", dna.display());
    io::read_one(&gb).unwrap()
}

/// A feature as compared across the round trip: name, type, range, and direction when it has
/// one (GenBank cannot say "no direction").
fn key(f: &Feature, len: usize) -> (String, String, Span, Option<Strand>) {
    let directed = matches!(f.strand, Strand::Forward | Strand::Reverse).then_some(f.strand);
    (f.name.clone(), f.kind.clone(), f.span(len), directed)
}

fn assert_same(written: &Document, read: &Document, what: &str) {
    assert_eq!(written.seq.to_ascii_uppercase(), read.seq.to_ascii_uppercase(), "{what}: sequence");
    assert_eq!(written.topology, read.topology, "{what}: topology");
    let len = written.len();
    for feature in written.features.iter().filter(|f| f.kind != "primer_bind") {
        let ours = key(feature, len);
        // SnapGene's GenBank export spells out letters such as α and Ψ ("lacZ-alpha"), so a name
        // with such letters is matched by type, range and direction alone.
        let found = read.features.iter().any(|f| {
            let theirs = key(f, len);
            let same_name = ours.0 == theirs.0 || !ours.0.is_ascii();
            same_name && ours.1 == theirs.1 && ours.2 == theirs.2 && (ours.3.is_none() || theirs.3.is_none() || ours.3 == theirs.3)
        });
        assert!(found, "{what}: SnapGene did not read back feature {ours:?}");
    }
    for primer in &written.primers {
        assert!(read.features.iter().any(|f| f.kind == "primer_bind" && f.name == primer.name), "{what}: primer {} missing", primer.name);
    }
}

#[test]
fn snapgene_reads_files_slopgene_wrote() {
    let Some(exe) = snapgene() else { return };
    let dir = std::env::temp_dir().join(format!("slopgene-snapgene-opens-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let rules = primers::BindingRules::default();

    // 1. A SnapGene file edited by SlopGene: bases inserted and deleted, a feature and a primer
    //    added.
    if let Some(samples) = samples() {
        let mut doc = io::read_one(samples.join("2. Experimental Design").join("2. Restriction Cloning").join("1. pUC57.dna")).unwrap();
        edit::insert(&mut doc, 100, b"GAATTCGAATTC").unwrap();
        edit::delete(&mut doc, Span::new(1500, 1520)).unwrap();
        doc.features.push(Feature::new("test feature", "misc_feature", Span::new(200, 300), Strand::Reverse));
        doc.primers.push(Primer::new("testF", "GTAAAACGACGGCCAGT"));
        primers::refresh_binding_sites(&mut doc, &rules);
        let path = dir.join("edited.dna");
        io::write_path(&doc, &path, None).unwrap();
        assert_same(&doc, &read_back_through_snapgene(&exe, &path), "edited");
    }

    // 2. A document made from nothing, with a feature through the origin.
    let seq: Vec<u8> = (0..3000u32).map(|ii| b"ACGT"[((ii * 7 + ii / 13) % 4) as usize]).collect();
    let mut doc = Document::new_dna("created", seq, Topology::Circular);
    doc.features.push(Feature::new("orf A", "CDS", Span::new(0, 300), Strand::Forward));
    doc.features.push(Feature::new("site B", "misc_feature", Span::new(900, 1000), Strand::Reverse));
    doc.features.push(Feature::new("wrap C", "promoter", Span::new(2900, 3060), Strand::Forward));
    let path = dir.join("created.dna");
    io::write_path(&doc, &path, None).unwrap();
    assert_same(&doc, &read_back_through_snapgene(&exe, &path), "created");

    // 3. A GenBank record converted to SnapGene.
    let record = io::genbank::write(&doc);
    let converted = io::read_bytes(record.as_bytes(), Some(io::Format::GenBank)).unwrap().remove(0);
    let path = dir.join("converted.dna");
    io::write_path(&converted, &path, None).unwrap();
    assert_same(&converted, &read_back_through_snapgene(&exe, &path), "converted");

    let _ = std::fs::remove_dir_all(&dir);
}
