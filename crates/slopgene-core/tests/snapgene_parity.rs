//! Black-box parity with SnapGene: run the workflows of SnapGene's own sample project and check
//! that SlopGene produces the same molecules SnapGene saved as results.
//!
//! The sample files ship with SnapGene and are not part of this repository, so these tests run
//! only where SnapGene is installed (or `SLOPGENE_SNAPGENE_SAMPLES` points at the sample
//! project), and pass without checking anything elsewhere.

use slopgene_core::cloning::{golden_gate, restriction_clone};
use slopgene_core::enzymes::EnzymeDb;
use slopgene_core::pcr::pcr;
use slopgene_core::primers::BindingRules;
use slopgene_core::protocol;
use slopgene_core::*;
use std::path::PathBuf;

fn samples() -> Option<PathBuf> {
    let dir = std::env::var("SLOPGENE_SNAPGENE_SAMPLES")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("C:\\Program Files\\SnapGene\\sampleData\\Sample project"));
    dir.exists().then_some(dir)
}

use slopgene_core::alphabet::same_circle;

#[test]
fn pcr_matches_snapgene() {
    let Some(dir) = samples() else { return };
    let design = dir.join("2. Experimental Design").join("1. PCR");
    let template = io::read_one(design.join("1. GFP.dna")).unwrap();
    let expected = io::read_one(design.join("2. GFP Amplified.dna")).unwrap();
    let fwd = template.primers.iter().find(|p| p.name == "GFP_FOR").unwrap();
    let rev = template.primers.iter().find(|p| p.name == "GFP_REV").unwrap();
    let product = pcr(&template, fwd, rev, &BindingRules::default()).unwrap();
    assert_eq!(product.doc.seq.to_ascii_uppercase(), expected.seq.to_ascii_uppercase());
    let ours: Vec<(String, Span)> = product.doc.features.iter().map(|f| (f.name.clone(), f.span(product.doc.len()))).collect();
    let theirs: Vec<(String, Span)> = expected.features.iter().map(|f| (f.name.clone(), f.span(expected.len()))).collect();
    assert_eq!(ours, theirs);
}

#[test]
fn restriction_cloning_matches_snapgene() {
    let Some(dir) = samples() else { return };
    let design = dir.join("2. Experimental Design").join("2. Restriction Cloning");
    let vector = io::read_one(design.join("1. pUC57.dna")).unwrap();
    let insert = io::read_one(design.join("2. mTAGBFP SalI, KpnI.dna")).unwrap();
    let expected = io::read_one(design.join("3. pUC57-mTAGBFP2.dna")).unwrap();
    let enzymes = EnzymeDb::builtin().resolve(&["SalI", "KpnI"]).unwrap();
    let products = restriction_clone(&vector, &insert, &enzymes).unwrap();
    assert_eq!(products.len(), 1, "SalI + KpnI cloning is directional");
    assert!(same_circle(&products[0].doc.seq, &expected.seq));
}

/// Reads a protocol's maps from one folder.
struct Folder(PathBuf);

impl protocol::Loader for Folder {
    fn load(&self, path: &str) -> Result<Document> {
        io::read_one(self.0.join(path))
    }
}

#[test]
fn a_cloning_protocol_matches_snapgene() {
    let Some(dir) = samples() else { return };
    let design = dir.join("2. Experimental Design").join("2. Restriction Cloning");
    let expected = io::read_one(design.join("3. pUC57-mTAGBFP2.dna")).unwrap();
    // The cloning above as a protocol: cut both maps with SalI and KpnI, keep the backbone and
    // the insert, dephosphorylate the backbone and ligate.
    let plan = protocol::Protocol::from_json(
        r#"{"name": "pUC57-mTAGBFP2", "steps": [
        {"id": "vector", "type": "open", "params": {"path": "1. pUC57.dna"}},
        {"id": "cut_vector", "in": ["vector"], "type": "digest", "params": {"enzymes": ["SalI", "KpnI"]}},
        {"id": "backbone", "in": ["cut_vector"], "type": "pick", "params": {"largest": true}},
        {"id": "cip", "in": ["backbone"], "type": "dephosphorylate"},
        {"id": "donor", "type": "open", "params": {"path": "2. mTAGBFP SalI, KpnI.dna"}},
        {"id": "cut_donor", "in": ["donor"], "type": "digest", "params": {"enzymes": ["SalI", "KpnI"]}},
        {"id": "insert", "in": ["cut_donor"], "type": "pick", "params": {"ends": "SalI/KpnI"}},
        {"id": "product", "in": ["cip", "insert"], "type": "ligate"}]}"#,
    )
    .unwrap();
    let run = protocol::run(&plan, &Folder(design), EnzymeDb::builtin());
    let product = run.step("product").unwrap();
    assert_eq!(product.status, protocol::Status::Done, "{:?}", product.error);
    assert!(same_circle(&product.outputs[0].doc.seq, &expected.seq));
}

#[test]
fn golden_gate_matches_snapgene() {
    let Some(dir) = samples() else { return };
    let design = dir.join("2. Experimental Design").join("3. Golden Gate Assembly");
    let vector = io::read_one(design.join("1. pGGA.dna")).unwrap();
    let part = io::read_one(design.join("2. mCherry BsaI.dna")).unwrap();
    let expected = io::read_one(design.join("3. pGGA-mCherry.dna")).unwrap();
    let bsai = EnzymeDb::builtin().get("BsaI").unwrap();
    let product = golden_gate(&[vector, part], bsai).unwrap();
    assert!(same_circle(&product.doc.seq, &expected.seq));
}
