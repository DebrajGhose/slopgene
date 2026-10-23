//! Physical properties of sequences: molecular weight, isoelectric point, extinction.

use crate::alphabet::base_counts;
use crate::model::{Document, Molecule, Strandedness};

/// Average mass (Da) of a single-stranded DNA oligo with a 5' hydroxyl, from base counts:
/// sum of nucleotide residue masses minus the terminal phosphate (61.96 Da), the formula used
/// by oligo suppliers.
pub fn ssdna_mass(seq: &[u8]) -> f64 {
    let [a, c, g, t, _] = base_counts(seq);
    let mass = a as f64 * 313.21 + c as f64 * 289.18 + g as f64 * 329.21 + t as f64 * 304.20;
    if mass > 0.0 { mass - 61.96 } else { 0.0 }
}

/// Average mass (Da) of double-stranded DNA: both strands as oligos.
pub fn dsdna_mass(seq: &[u8], circular: bool) -> f64 {
    let rc = crate::alphabet::revcomp(seq);
    let linear = ssdna_mass(seq) + ssdna_mass(&rc);
    // Closing a circle removes one water per strand.
    if circular { linear - 2.0 * 18.02 } else { linear }
}

/// Average mass (Da) of single-stranded RNA with a 5' hydroxyl (synthetic RNA; a transcript
/// with a 5' triphosphate weighs about 159 Da more).
pub fn rna_mass(seq: &[u8]) -> f64 {
    let [a, c, g, u, _] = base_counts(seq);
    let mass = a as f64 * 329.21 + c as f64 * 305.18 + g as f64 * 345.21 + u as f64 * 306.17;
    if mass > 0.0 { mass - 61.96 } else { 0.0 }
}

/// Average residue mass of an amino acid (Da), the mass it adds inside a chain.
fn residue_mass(aa: u8) -> Option<f64> {
    Some(match aa.to_ascii_uppercase() {
        b'A' => 71.0788,
        b'R' => 156.1875,
        b'N' => 114.1038,
        b'D' => 115.0886,
        b'C' => 103.1388,
        b'E' => 129.1155,
        b'Q' => 128.1307,
        b'G' => 57.0519,
        b'H' => 137.1411,
        b'I' => 113.1594,
        b'L' => 113.1594,
        b'K' => 128.1741,
        b'M' => 131.1926,
        b'F' => 147.1766,
        b'P' => 97.1167,
        b'S' => 87.0782,
        b'T' => 101.1051,
        b'W' => 186.2132,
        b'Y' => 163.1760,
        b'V' => 99.1326,
        b'U' => 150.0388,
        b'O' => 237.3018,
        _ => return None,
    })
}

/// Average mass of a protein (Da): residue masses plus one water. Stops and unknown codes are
/// skipped.
pub fn protein_mass(protein: &[u8]) -> f64 {
    let residues: f64 = protein.iter().filter_map(|&aa| residue_mass(aa)).sum();
    if residues > 0.0 { residues + 18.015 } else { 0.0 }
}

/// Net charge of a protein at `ph`, from EMBOSS pKa values.
fn net_charge(protein: &[u8], ph: f64) -> f64 {
    let count = |aa: u8| protein.iter().filter(|&&b| b.to_ascii_uppercase() == aa).count() as f64;
    let positive = |pka: f64| 1.0 / (1.0 + 10f64.powf(ph - pka));
    let negative = |pka: f64| -1.0 / (1.0 + 10f64.powf(pka - ph));
    positive(8.6)
        + negative(3.6)
        + count(b'K') * positive(10.8)
        + count(b'R') * positive(12.5)
        + count(b'H') * positive(6.5)
        + count(b'D') * negative(3.9)
        + count(b'E') * negative(4.1)
        + count(b'C') * negative(8.5)
        + count(b'Y') * negative(10.1)
}

/// Isoelectric point: the pH where the net charge is zero, found by bisection.
pub fn isoelectric_point(protein: &[u8]) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 14.0f64);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if net_charge(protein, mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Molar extinction coefficient at 280 nm (M^-1 cm^-1), assuming every cysteine pair forms a
/// disulfide (Pace et al. 1995).
pub fn extinction_coefficient(protein: &[u8]) -> f64 {
    let count = |aa: u8| protein.iter().filter(|&&b| b.to_ascii_uppercase() == aa).count() as f64;
    count(b'W') * 5500.0 + count(b'Y') * 1490.0 + (count(b'C') / 2.0).floor() * 125.0
}

/// Mass of a document's molecule (Da).
pub fn document_mass(doc: &Document) -> f64 {
    match doc.molecule {
        Molecule::Protein => protein_mass(&doc.seq),
        Molecule::Rna => rna_mass(&doc.seq),
        Molecule::Dna if doc.strandedness == Strandedness::Single => ssdna_mass(&doc.seq),
        Molecule::Dna => dsdna_mass(&doc.seq, doc.is_circular()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn masses_are_plausible() {
        // A 20-mer oligo weighs about 6 kDa.
        let m = ssdna_mass(b"ACGTACGTACGTACGTACGT");
        assert!((6000.0..6300.0).contains(&m), "{m}");
        // pUC19-sized dsDNA: ~650 Da per bp.
        let ds = dsdna_mass(&vec![b'A'; 2686], true);
        assert!((1.6e6..1.7e6).contains(&ds), "{ds}");
    }

    #[test]
    fn protein_properties() {
        // GFP-like stretch: mass and pI in sensible ranges.
        let p = b"MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTF";
        let mass = protein_mass(p);
        assert!((6500.0..7500.0).contains(&mass), "{mass}");
        let pi = isoelectric_point(p);
        assert!((4.0..7.0).contains(&pi), "{pi}");
        assert_eq!(extinction_coefficient(b"WYCC"), 5500.0 + 1490.0 + 125.0);
    }
}
