//! Cloning simulations: restriction cloning, homology-based assembly (Gibson, NEBuilder HiFi,
//! In-Fusion) and Golden Gate assembly.
//!
//! Each simulation returns the product as a [`Document`] carrying the features of its parts,
//! plus notes describing what was joined, so results can be checked or reported.

use crate::alphabet;
use crate::digest::{circularize, digest, ends_compatible, ligate, merge_split_features};
use crate::edit;
use crate::enzymes::{Enzyme, find_sites_in};
use crate::error::{Error, Result};
use crate::model::*;
use crate::primers::{TmConditions, melting_temperature};
use crate::span::Span;

/// A simulated cloning product.
#[derive(Clone, Debug)]
pub struct Product {
    pub doc: Document,
    /// What was joined, in order, for the record.
    pub notes: Vec<String>,
}

/// A linear copy of `doc` flipped end to end (reverse complement), with its ends swapped.
fn flipped(doc: &Document) -> Document {
    let mut copy = doc.clone();
    edit::reverse_complement(&mut copy);
    copy.name = format!("{} (reversed)", doc.name);
    copy
}

/// Restriction cloning: cut `vector` and `insert` with the same enzymes, keep the vector
/// backbone (its largest fragment) and every insert fragment whose two ends fit it, then ligate
/// and close the circle. Both orientations are tried, so a non-directional cloning gives two
/// products.
pub fn restriction_clone(vector: &Document, insert: &Document, enzymes: &[&Enzyme]) -> Result<Vec<Product>> {
    let vector_cut = digest(vector, enzymes)?;
    if vector_cut.sites.is_empty() {
        return Err(Error::Invalid(format!("the enzymes do not cut the vector {}", vector.name)));
    }
    let backbone = vector_cut.fragments.iter().max_by_key(|f| f.len()).ok_or_else(|| Error::Invalid("the vector gave no fragments".into()))?;
    if backbone.doc.is_circular() {
        return Err(Error::Invalid("the vector is cut nowhere".into()));
    }
    let insert_cut = digest(insert, enzymes)?;
    // Insert pieces with enzyme-made ends on both sides (not the leftover ends of a PCR product).
    let candidates: Vec<&crate::digest::Fragment> =
        insert_cut.fragments.iter().filter(|f| f.left_enzyme.is_some() && f.right_enzyme.is_some() && !f.doc.is_circular()).collect();
    if candidates.is_empty() {
        return Err(Error::Invalid(format!("the enzymes do not release a fragment from {}", insert.name)));
    }

    let mut products = Vec::new();
    for candidate in candidates {
        for reversed in [false, true] {
            let piece = if reversed { flipped(&candidate.doc) } else { candidate.doc.clone() };
            if !(ends_compatible(&backbone.doc, &piece) && ends_compatible(&piece, &backbone.doc)) {
                continue;
            }
            let joined = ligate(&backbone.doc, &piece)?;
            let mut circle = circularize(&joined)?;
            // The product is a plasmid to transform, so it is methylated as E. coli leaves it.
            circle.methylation = Methylation::E_COLI;
            circle.name = format!("{} + {}", vector.name, insert.name);
            circle.meta.description =
                format!("{} cloned into {} with {}", insert.name, vector.name, enzymes.iter().map(|e| e.name.as_str()).collect::<Vec<_>>().join(" + "));
            circle.meta.synthetic = Some(true);
            let notes = vec![
                format!(
                    "vector backbone {} bp ({} to {})",
                    backbone.len(),
                    backbone.left_enzyme.clone().unwrap_or_default(),
                    backbone.right_enzyme.clone().unwrap_or_default()
                ),
                format!("insert {} bp{}", piece.len(), if reversed { ", reversed" } else { "" }),
                format!("product {} bp", circle.len()),
            ];
            products.push(Product { doc: circle, notes });
        }
    }
    if products.is_empty() {
        return Err(Error::Invalid(format!(
            "the insert's ends do not fit the vector backbone (backbone ends: {} / {})",
            backbone.doc.ends.left.describe(),
            backbone.doc.ends.right.describe()
        )));
    }
    if products.len() > 1 {
        for p in &mut products {
            p.notes.push("non-directional: the insert fits both ways round".into());
        }
    }
    Ok(products)
}

/// Longest `k` in `min..=max` such that the last `k` bases of `a` equal the first `k` of `b`.
fn terminal_overlap(a: &[u8], b: &[u8], min: usize, max: usize) -> Option<usize> {
    let max = max.min(a.len()).min(b.len());
    (min..=max).rev().find(|&k| a[a.len() - k..].eq_ignore_ascii_case(&b[..k]))
}

/// Homology-based assembly (Gibson, NEBuilder HiFi, In-Fusion): neighbouring fragments must
/// share an end overlap of at least `min_overlap` bases. A fragment that only matches in
/// reverse is flipped. With `circular`, the last fragment's end must also overlap the first
/// fragment's start, and the product is a circle.
pub fn homology_assembly(fragments: &[Document], min_overlap: usize, circular: bool) -> Result<Product> {
    if fragments.is_empty() {
        return Err(Error::Invalid("give at least one fragment".into()));
    }
    if fragments.len() == 1 && !circular {
        return Err(Error::Invalid("a single fragment can only be assembled into a circle".into()));
    }
    let max_overlap = 300;
    let mut notes = Vec::new();
    let mut product = fragments[0].clone();
    product.raw = RawData::default();
    if product.is_circular() {
        return Err(Error::Invalid(format!("{} is circular; linearise it first (digest or PCR)", product.name)));
    }
    product.ends = Ends::default();
    let conditions = TmConditions::default();
    // The fragment each junction follows, for the notes.
    let mut previous = fragments[0].name.as_str();
    for next in &fragments[1..] {
        if next.is_circular() {
            return Err(Error::Invalid(format!("{} is circular; linearise it first (digest or PCR)", next.name)));
        }
        let mut piece = next.clone();
        piece.ends = Ends::default();
        let mut k = terminal_overlap(&product.seq, &piece.seq, min_overlap, max_overlap);
        if k.is_none() {
            let reversed = flipped(&piece);
            if let Some(kr) = terminal_overlap(&product.seq, &reversed.seq, min_overlap, max_overlap) {
                piece = reversed;
                k = Some(kr);
                notes.push(format!("{} joined reversed", next.name));
            }
        }
        let k = k.ok_or_else(|| Error::Invalid(format!("no overlap of at least {min_overlap} bases between {} and {}", product.name, next.name)))?;
        let overlap = &piece.seq[..k];
        notes.push(format!("{previous} -> {}: {k} bp overlap, Tm {:.0} C", next.name, melting_temperature(overlap, &conditions)));
        previous = next.name.as_str();
        // The overlap is the end of `product`; append only what follows it.
        let offset = product.len() - k;
        product.seq.truncate(offset);
        product.seq.extend_from_slice(&piece.seq);
        for f in &piece.features {
            let mut f = f.clone();
            f.uid = new_uid();
            for s in &mut f.segments {
                s.span = Span::new(s.span.start + offset, s.span.end + offset);
            }
            product.features.push(f);
        }
        for p in &piece.primers {
            let mut p = p.clone();
            p.uid = new_uid();
            for site in &mut p.sites {
                site.span = Span::new(site.span.start + offset, site.span.end + offset);
            }
            product.primers.push(p);
        }
    }
    if circular {
        let first_len = fragments[0].len().min(product.len());
        let k = terminal_overlap(&product.seq, &product.seq[..first_len], min_overlap, max_overlap.min(product.len() / 2)).ok_or_else(|| {
            Error::Invalid(format!("the last fragment does not overlap the first by at least {min_overlap} bases, so the assembly cannot close"))
        })?;
        notes.push(format!("closing overlap {k} bp"));
        let new_len = product.len() - k;
        product.seq.truncate(new_len);
        product.topology = Topology::Circular;
        product.methylation = Methylation::E_COLI;
        // Features that ran into the duplicated overlap now wrap through the origin.
        for f in &mut product.features {
            for s in &mut f.segments {
                if s.span.start >= new_len {
                    s.span = Span::new(s.span.start - new_len, s.span.end - new_len);
                }
                s.span = Span::new(s.span.start, s.span.end.min(s.span.start + new_len));
            }
            f.normalize(new_len, true);
        }
        product.primers.iter_mut().for_each(|p| p.sites.retain(|s| s.span.end <= new_len + k));
    }
    merge_split_features(&mut product);
    product.name = fragments.iter().map(|f| f.name.as_str()).collect::<Vec<_>>().join(" + ");
    product.meta.synthetic = Some(true);
    product.meta.description = format!("Homology assembly of {}", product.name);
    notes.push(format!("product {} bp, {}", product.len(), product.topology.name()));
    Ok(Product { doc: product, notes })
}

/// Golden Gate assembly with one Type IIS enzyme. Each part is cut; from each, the piece with
/// enzyme-made ends and no site left in it is kept. Pieces are chained so every overhang meets
/// its partner (flipping pieces when needed) and closed into a circle.
pub fn golden_gate(parts: &[Document], enzyme: &Enzyme) -> Result<Product> {
    if parts.len() < 2 {
        return Err(Error::Invalid("Golden Gate needs a vector and at least one part".into()));
    }
    // Usable pieces per part, in both orientations.
    let mut options: Vec<Vec<Document>> = Vec::new();
    for part in parts {
        let cut = digest(part, &[enzyme])?;
        let pieces: Vec<Document> = cut
            .fragments
            .into_iter()
            .filter(|f| f.left_enzyme.is_some() && f.right_enzyme.is_some() && !f.doc.is_circular())
            .filter(|f| find_sites_in(&f.doc, &[enzyme]).is_empty())
            .map(|f| f.doc)
            .collect();
        if pieces.is_empty() {
            return Err(Error::Invalid(format!(
                "{} gives no {}-flanked piece free of {} sites (check the sites and for internal sites)",
                part.name, enzyme.name, enzyme.name
            )));
        }
        // Try each piece as given before reversed, so the product keeps the vector's orientation.
        let mut both = Vec::new();
        for p in pieces {
            let reversed = flipped(&p);
            both.push(p);
            both.push(reversed);
        }
        options.push(both);
    }

    // Depth-first search for an order: start from the first part (the vector), then repeatedly
    // add an unused part whose piece's left overhang matches the current right overhang.
    fn search(chain: &mut Vec<Document>, used: &mut Vec<bool>, options: &[Vec<Document>]) -> bool {
        if used.iter().all(|&u| u) {
            return ends_compatible(chain.last().unwrap(), &chain[0]);
        }
        for ii in 0..options.len() {
            if used[ii] {
                continue;
            }
            for piece in &options[ii] {
                if ends_compatible(chain.last().unwrap(), piece) {
                    used[ii] = true;
                    chain.push(piece.clone());
                    if search(chain, used, options) {
                        return true;
                    }
                    chain.pop();
                    used[ii] = false;
                }
            }
        }
        false
    }
    for start in &options[0] {
        let mut chain = vec![start.clone()];
        let mut used = vec![false; parts.len()];
        used[0] = true;
        if search(&mut chain, &mut used, &options) {
            let mut notes = Vec::new();
            for piece in &chain {
                notes.push(format!(
                    "{} ({} bp): {} | {}",
                    piece.name,
                    piece.len(),
                    String::from_utf8_lossy(&crate::digest::overhang_bases(piece, true)),
                    String::from_utf8_lossy(&crate::digest::overhang_bases(piece, false))
                ));
            }
            let mut joined = chain[0].clone();
            for piece in &chain[1..] {
                joined = ligate(&joined, piece)?;
            }
            let mut circle = circularize(&joined)?;
            circle.methylation = Methylation::E_COLI;
            circle.name = parts.iter().map(|p| p.name.as_str()).collect::<Vec<_>>().join(" + ");
            circle.meta.synthetic = Some(true);
            circle.meta.description = format!("Golden Gate assembly with {} of {}", enzyme.name, circle.name);
            notes.push(format!("product {} bp", circle.len()));
            return Ok(Product { doc: circle, notes });
        }
    }
    Err(Error::Invalid(format!("the {} overhangs of the parts do not form a closed chain", enzyme.name)))
}

/// Reverse-complement helper for building test parts and primers.
pub fn revcomp_text(text: &str) -> String {
    String::from_utf8_lossy(&alphabet::revcomp(text.as_bytes())).into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzymes::EnzymeDb;
    use crate::enzymes::find_sites;

    fn random_dna(n: usize, seed: u64) -> Vec<u8> {
        // A small linear congruential generator keeps tests deterministic.
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    /// Random sequence with no site of the given enzymes.
    fn clean_dna(n: usize, seed: u64, avoid: &[&Enzyme]) -> Vec<u8> {
        let mut s = seed;
        loop {
            let seq = random_dna(n, s);
            if find_sites(&seq, false, avoid).is_empty() {
                return seq;
            }
            s += 1;
        }
    }

    #[test]
    fn directional_restriction_cloning() {
        let db = EnzymeDb::builtin();
        let enzymes = db.resolve(&["EcoRI", "BamHI"]).unwrap();
        let mut vector_seq = clean_dna(300, 1, &enzymes);
        vector_seq.splice(100..100, b"GAATTCAAAAAAGGATCC".iter().copied());
        let vector = Document::new_dna("pVec", vector_seq, Topology::Circular);
        let core = clean_dna(120, 7, &enzymes);
        let mut insert_seq = b"CCGAATTC".to_vec();
        insert_seq.extend_from_slice(&core);
        insert_seq.extend_from_slice(b"GGATCCGG");
        let mut insert = Document::new_dna("ins", insert_seq, Topology::Linear);
        insert.features.push(Feature::new("gene", "CDS", Span::new(8, 128), Strand::Forward));

        let products = restriction_clone(&vector, &insert, &enzymes).unwrap();
        assert_eq!(products.len(), 1, "EcoRI + BamHI is directional");
        let product = &products[0].doc;
        assert!(product.is_circular());
        // The product contains EcoRI, the insert, then BamHI, once each.
        let mut expected = b"GAATTC".to_vec();
        expected.extend_from_slice(&core);
        expected.extend_from_slice(b"GGATCC");
        let doubled = [product.seq.clone(), product.seq.clone()].concat();
        assert!(doubled.windows(expected.len()).any(|w| w == expected.as_slice()));
        // 318 bp vector, minus the 12 bases between the EcoRI and BamHI cuts, plus 126 insert bases.
        assert_eq!(product.len(), 432);
        assert!(product.features.iter().any(|f| f.name == "gene" && f.len() == 120));
    }

    #[test]
    fn gibson_closes_a_circle() {
        let a = random_dna(200, 11);
        let b = random_dna(150, 12);
        let c = random_dna(180, 13);
        // Each fragment repeats the last 25 bases of the one before; the last overlaps the first.
        let mut f1 = c[155..].to_vec();
        f1.extend_from_slice(&a);
        let mut f2 = a[175..].to_vec();
        f2.extend_from_slice(&b);
        let mut f3 = b[125..].to_vec();
        f3.extend_from_slice(&c);
        let frags = vec![
            Document::new_dna("f1", f1, Topology::Linear),
            Document::new_dna("f2", crate::alphabet::revcomp(&f2), Topology::Linear),
            Document::new_dna("f3", f3, Topology::Linear),
        ];
        let product = homology_assembly(&frags, 20, true).unwrap();
        assert_eq!(product.doc.len(), 200 + 150 + 180);
        let expected = [a.clone(), b.clone(), c.clone()].concat();
        let doubled = [product.doc.seq.clone(), product.doc.seq.clone()].concat();
        assert!(doubled.windows(expected.len()).any(|w| w == expected.as_slice()));
        assert!(product.notes.iter().any(|n| n.contains("reversed")), "f2 was given reversed and must be flipped");
        // Each junction note names the two fragments it joins.
        assert!(product.notes.iter().any(|n| n.starts_with("f1 -> f2:")), "{:?}", product.notes);
        assert!(product.notes.iter().any(|n| n.starts_with("f2 -> f3:")), "{:?}", product.notes);
    }

    #[test]
    fn golden_gate_two_parts() {
        let db = EnzymeDb::builtin();
        let bsai = db.get("BsaI").unwrap();
        // Overhangs: vector right AATG -> part1 -> GCTT -> part2 -> vector left TTGA.
        let backbone = clean_dna(400, 21, &[bsai]);
        let mut vector = backbone.clone();
        // Dropout cassette with outward... BsaI sites so the backbone keeps TTGA and AATG ends.
        vector.extend_from_slice(b"AATGTGAGACCAAAAAAAAAAGGTCTCATTGA");
        let vector = Document::new_dna("pGG", vector, Topology::Circular);
        let core1 = clean_dna(90, 31, &[bsai]);
        let core2 = clean_dna(110, 41, &[bsai]);
        let part = |left: &[u8], core: &[u8], right: &[u8], name: &str| {
            let mut s = b"CCGGTCTCA".to_vec();
            s.extend_from_slice(left);
            s.extend_from_slice(core);
            s.extend_from_slice(right);
            s.extend_from_slice(b"TGAGACCGG");
            Document::new_dna(name, s, Topology::Linear)
        };
        let p1 = part(b"AATG", &core1, b"GCTT", "p1");
        let p2 = part(b"GCTT", &core2, b"TTGA", "p2");
        let product = golden_gate(&[vector, p2, p1], bsai).unwrap();
        assert!(product.doc.is_circular());
        let mut expected = b"AATG".to_vec();
        expected.extend_from_slice(&core1);
        expected.extend_from_slice(b"GCTT");
        expected.extend_from_slice(&core2);
        expected.extend_from_slice(b"TTGA");
        let doubled = [product.doc.seq.clone(), product.doc.seq.clone()].concat();
        assert!(doubled.windows(expected.len()).any(|w| w == expected.as_slice()), "{:?}", product.notes);
        assert!(find_sites(&product.doc.seq, true, &[bsai]).is_empty(), "no BsaI site survives");
    }
}
