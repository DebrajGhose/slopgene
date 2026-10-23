//! Restriction digests: cutting a molecule into fragments with exact sticky ends.
//!
//! A fragment is a [`Document`] whose sequence spans both strands, so overhanging bases are
//! part of the sequence and [`Ends`] says which strand they belong to. Two ends can be ligated
//! when they are the same kind (blunt, 5' or 3') and length and their overhang bases are the
//! same when read on the top strand, which is what complementary single strands look like there.

use crate::edit;
use crate::enzymes::{Cut, Enzyme, Site, find_sites};
use crate::error::{Result, invalid};
use crate::model::*;
use crate::span::Span;

/// One fragment of a digest.
#[derive(Clone, Debug)]
pub struct Fragment {
    /// The fragment as a linear document (sequence covers both strands, ends set).
    pub doc: Document,
    /// Where the fragment came from in the parent, as a span over the parent's sequence
    /// (union of both strands).
    pub span: Span,
    /// Enzymes that made the left and right ends (`None` for an original end of a linear parent).
    pub left_enzyme: Option<String>,
    pub right_enzyme: Option<String>,
}

impl Fragment {
    pub fn len(&self) -> usize {
        self.doc.len()
    }

    pub fn is_empty(&self) -> bool {
        self.doc.is_empty()
    }
}

/// Result of a digest.
#[derive(Clone, Debug)]
pub struct Digest {
    pub fragments: Vec<Fragment>,
    /// Every site that cut, sorted by position.
    pub sites: Vec<Site>,
    /// Enzymes that did not cut at all.
    pub non_cutters: Vec<String>,
}

/// Cut `doc` with `enzymes`. A circular molecule with no cuts comes back as one uncut circular
/// fragment; with `k` cuts it gives `k` fragments. Fragments are listed in sequence order,
/// starting from the first cut.
pub fn digest(doc: &Document, enzymes: &[&Enzyme]) -> Result<Digest> {
    if !doc.molecule.is_nucleic() {
        return Err(invalid("only DNA can be digested"));
    }
    let seq_len = doc.len();
    let circular = doc.is_circular();
    let sites: Vec<Site> = find_sites(&doc.seq, circular, enzymes).into_iter().filter(|s| !s.cuts.is_empty()).collect();
    let cutters: std::collections::HashSet<&str> = sites.iter().map(|s| s.enzyme.as_str()).collect();
    let non_cutters = enzymes.iter().filter(|e| !cutters.contains(e.name.as_str())).map(|e| e.name.clone()).collect();

    // All cuts, labelled with their enzyme. Identical cuts (isoschizomers chosen together)
    // count once.
    let mut cuts: Vec<(Cut, String)> = Vec::new();
    for site in &sites {
        for cut in &site.cuts {
            if !cuts.iter().any(|(c, _)| c == cut) {
                cuts.push((*cut, site.enzyme.clone()));
            }
        }
    }
    cuts.sort_by_key(|(c, _)| (c.pos, c.overhang));

    let mut fragments = Vec::new();
    if cuts.is_empty() {
        let mut whole = doc.clone();
        whole.raw = RawData::default();
        fragments.push(Fragment { doc: whole, span: Span::new(0, seq_len), left_enzyme: None, right_enzyme: None });
        return Ok(Digest { fragments, sites, non_cutters });
    }

    // Unwrapped boundaries of a fragment: its left end starts where the earlier strand of the
    // left cut begins, its right end stops where the later strand of the right cut ends.
    let make = |left: Option<(Cut, &str)>, right: Option<(Cut, &str)>, wrap_right: bool| -> Result<Fragment> {
        let (start, left_end) = match left {
            Some((cut, _)) => {
                let p = cut.pos as i64;
                let o = cut.overhang as i64;
                let end = match o {
                    0 => Overhang::Blunt,
                    o if o > 0 => Overhang::FivePrime(o as usize),
                    o => Overhang::ThreePrime((-o) as usize),
                };
                (p.min(p + o), end)
            }
            None => (0, doc.ends.left),
        };
        let (end, right_end) = match right {
            Some((cut, _)) => {
                let q = cut.pos as i64 + if wrap_right { seq_len as i64 } else { 0 };
                let o = cut.overhang as i64;
                let e = match o {
                    0 => Overhang::Blunt,
                    o if o > 0 => Overhang::FivePrime(o as usize),
                    o => Overhang::ThreePrime((-o) as usize),
                };
                (q.max(q + o), e)
            }
            None => (seq_len as i64, doc.ends.right),
        };
        if end <= start {
            return Err(invalid("cuts overlap so closely that a fragment has no bases"));
        }
        // Bring the start into 0..L (it can be negative when an overhang crosses the origin).
        let len = (end - start) as usize;
        let start = start.rem_euclid(seq_len as i64) as usize;
        let span = Span::new(start, start + len);
        let mut piece = edit::extract(doc, span)?;
        piece.ends = Ends {
            left: left_end,
            right: right_end,
            left_modification: Some("FivePrimePhosphorylated".into()),
            right_modification: Some("FivePrimePhosphorylated".into()),
        };
        let label = |side: Option<(Cut, &str)>| side.map(|(_, e)| e.to_string());
        piece.name = format!("{} {}-{} fragment", doc.name, label(left).unwrap_or_else(|| "end".into()), label(right).unwrap_or_else(|| "end".into()))
            .trim()
            .to_string();
        Ok(Fragment { doc: piece, span, left_enzyme: label(left), right_enzyme: label(right) })
    };

    if circular {
        for ii in 0..cuts.len() {
            let (left_cut, left_name) = &cuts[ii];
            let (right_cut, right_name) = &cuts[(ii + 1) % cuts.len()];
            // The fragment after the last cut runs through the origin to the first cut; with a
            // single cut the molecule is simply opened up.
            let wraps = ii + 1 == cuts.len();
            fragments.push(make(Some((*left_cut, left_name)), Some((*right_cut, right_name)), wraps)?);
        }
    } else {
        fragments.push(make(None, Some((cuts[0].0, &cuts[0].1)), false)?);
        for ii in 0..cuts.len() - 1 {
            fragments.push(make(Some((cuts[ii].0, &cuts[ii].1)), Some((cuts[ii + 1].0, &cuts[ii + 1].1)), false)?);
        }
        let last = cuts.len() - 1;
        fragments.push(make(Some((cuts[last].0, &cuts[last].1)), None, false)?);
    }
    Ok(Digest { fragments, sites, non_cutters })
}

/// Bases of an end's overhang, as read on the top strand.
pub fn overhang_bases(doc: &Document, left: bool) -> Vec<u8> {
    let end = if left { doc.ends.left } else { doc.ends.right };
    let n = end.len().min(doc.len());
    if left { doc.seq[..n].to_ascii_uppercase() } else { doc.seq[doc.len() - n..].to_ascii_uppercase() }
}

/// Whether the right end of `a` can be ligated to the left end of `b`.
pub fn ends_compatible(a: &Document, b: &Document) -> bool {
    let (right, left) = (a.ends.right, b.ends.left);
    let same_kind = matches!(
        (right, left),
        (Overhang::Blunt, Overhang::Blunt) | (Overhang::FivePrime(_), Overhang::FivePrime(_)) | (Overhang::ThreePrime(_), Overhang::ThreePrime(_))
    );
    same_kind && right.len() == left.len() && overhang_bases(a, false) == overhang_bases(b, true)
}

/// Ligate linear fragments in order (`a` then `b`), merging matching overhangs. Use
/// [`circularize`] to close a molecule on itself.
pub fn ligate(a: &Document, b: &Document) -> Result<Document> {
    if !ends_compatible(a, b) {
        return Err(invalid(format!("ends are not compatible: {} ({}) then {} ({})", a.name, a.ends.right.describe(), b.name, b.ends.left.describe())));
    }
    let overlap = a.ends.right.len();
    let mut out = a.clone();
    out.raw = RawData::default();
    let offset = a.len() - overlap;
    out.seq.truncate(offset);
    out.seq.extend_from_slice(&b.seq);
    for feature in &b.features {
        let mut f = feature.clone();
        f.uid = new_uid();
        for s in &mut f.segments {
            s.span = Span::new(s.span.start + offset, s.span.end + offset);
        }
        out.features.push(f);
    }
    for primer in &b.primers {
        let mut p = primer.clone();
        p.uid = new_uid();
        for site in &mut p.sites {
            site.span = Span::new(site.span.start + offset, site.span.end + offset);
        }
        out.primers.push(p);
    }
    out.ends.right = b.ends.right;
    out.ends.right_modification = b.ends.right_modification.clone();
    out.name = format!("{} + {}", a.name, b.name);
    merge_split_features(&mut out);
    Ok(out)
}

/// Close a linear molecule into a circle when its two ends are compatible.
pub fn circularize(doc: &Document) -> Result<Document> {
    if !ends_compatible(doc, doc) {
        return Err(invalid(format!("the ends of {} are not compatible ({} and {})", doc.name, doc.ends.left.describe(), doc.ends.right.describe())));
    }
    let overlap = doc.ends.right.len();
    let mut out = doc.clone();
    out.raw = RawData::default();
    let new_len = doc.len() - overlap;
    out.seq.truncate(new_len);
    out.topology = Topology::Circular;
    out.ends = Ends::default();
    // Annotations inside the dropped overhang copy now sit at the start of the circle.
    for feature in &mut out.features {
        for s in &mut feature.segments {
            if s.span.start >= new_len {
                s.span = Span::new(s.span.start - new_len, s.span.end - new_len);
            }
            s.span = Span::new(s.span.start, s.span.end.min(s.span.start + new_len));
        }
        feature.normalize(new_len, true);
    }
    for primer in &mut out.primers {
        primer.sites.retain(|s| s.span.end <= new_len + overlap);
        for site in &mut primer.sites {
            site.span = site.span.normalize(new_len);
        }
    }
    merge_split_features(&mut out);
    Ok(out)
}

/// After a ligation the two halves of a feature that the cut split can sit side by side. When
/// two features have the same name, type and strand and touch or overlap, join them back.
pub fn merge_split_features(doc: &mut Document) {
    let seq_len = doc.len();
    let circular = doc.is_circular();
    let mut ii = 0;
    while ii < doc.features.len() {
        let mut iii = ii + 1;
        let mut merged = false;
        while iii < doc.features.len() {
            let (a, b) = (&doc.features[ii], &doc.features[iii]);
            if a.name == b.name
                && a.kind == b.kind
                && a.strand == b.strand
                && a.segments.len() == 1
                && b.segments.len() == 1
                && let Some(span) = joined(a.segments[0].span, b.segments[0].span, seq_len, circular)
            {
                doc.features[ii].segments[0].span = span;
                doc.features.remove(iii);
                merged = true;
                continue;
            }
            iii += 1;
        }
        if !merged {
            ii += 1;
        }
    }
}

/// The span that covers both `a` and `b` when they touch or overlap, or `None` when a gap
/// separates them. On a circle either span may run through the origin, and the two can meet
/// across it, as the halves of a feature cut once do when the circle closes again.
fn joined(a: Span, b: Span, seq_len: usize, circular: bool) -> Option<Span> {
    let turn = seq_len.max(1) as i64;
    // Place b as written, one turn back and one turn on, and look for where it meets a.
    let shifts: &[i64] = if circular { &[0, -1, 1] } else { &[0] };
    let (a0, a1) = (a.start as i64, a.end as i64);
    for shift in shifts {
        let (b0, b1) = (b.start as i64 + shift * turn, b.end as i64 + shift * turn);
        if b0 <= a1 && a0 <= b1 {
            let start = a0.min(b0);
            let len = ((a1.max(b1) - start) as usize).min(seq_len);
            let start = start.rem_euclid(turn) as usize;
            return Some(Span::new(start, start + len));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzymes::EnzymeDb;

    fn plasmid_with_two_sites() -> Document {
        // EcoRI at 10 and BamHI at 40 in a 60 bp circle.
        let mut seq = vec![b'A'; 60];
        seq[10..16].copy_from_slice(b"GAATTC");
        seq[40..46].copy_from_slice(b"GGATCC");
        let mut doc = Document::new_dna("p", seq, Topology::Circular);
        doc.features.push(Feature::new("insert", "misc_feature", Span::new(20, 30), Strand::Forward));
        doc
    }

    #[test]
    fn circular_double_digest() {
        let db = EnzymeDb::builtin();
        let doc = plasmid_with_two_sites();
        let enzymes = db.resolve(&["EcoRI", "BamHI", "NotI"]).unwrap();
        let result = digest(&doc, &enzymes).unwrap();
        assert_eq!(result.non_cutters, vec!["NotI".to_string()]);
        assert_eq!(result.fragments.len(), 2);
        let small = &result.fragments[0];
        // EcoRI cuts the top strand at 11 (G^AATTC); BamHI at 41 (G^GATCC). The small fragment
        // spans 11..45: AATTC...G plus the GATC overhang on the bottom strand.
        assert_eq!(small.doc.ends.left, Overhang::FivePrime(4));
        assert_eq!(small.doc.ends.right, Overhang::FivePrime(4));
        assert_eq!(small.span, Span::new(11, 45));
        assert_eq!(&small.doc.seq[..4], b"AATT");
        assert_eq!(&small.doc.seq[small.len() - 4..], b"GATC");
        assert_eq!(small.doc.features.len(), 1);
        assert_eq!(small.doc.features[0].segments[0].span, Span::new(9, 19));
        // The large fragment runs through the origin and has the matching opposite ends.
        let large = &result.fragments[1];
        assert_eq!(large.span, Span::new(41, 75));
        assert!(ends_compatible(&small.doc, &large.doc));
        assert!(ends_compatible(&large.doc, &small.doc));
    }

    #[test]
    fn religation_restores_the_plasmid() {
        let db = EnzymeDb::builtin();
        let doc = plasmid_with_two_sites();
        let enzymes = db.resolve(&["EcoRI", "BamHI"]).unwrap();
        let result = digest(&doc, &enzymes).unwrap();
        let (small, large) = (&result.fragments[0].doc, &result.fragments[1].doc);
        let joined = ligate(small, large).unwrap();
        let circle = circularize(&joined).unwrap();
        assert_eq!(circle.len(), 60);
        // Same molecule, possibly with a different origin.
        let doubled = [doc.seq.clone(), doc.seq.clone()].concat();
        assert!(doubled.windows(60).any(|w| w == circle.seq.as_slice()));
    }

    #[test]
    fn one_cut_opens_the_circle_with_the_overhang_at_both_ends() {
        let db = EnzymeDb::builtin();
        // EcoRI cuts the 60 bp circle once, inside a feature that covers bases 5..20.
        let mut doc = plasmid_with_two_sites();
        doc.features.push(Feature::new("lacZ", "CDS", Span::new(5, 20), Strand::Forward));
        let result = digest(&doc, &db.resolve(&["EcoRI"]).unwrap()).unwrap();
        assert_eq!(result.fragments.len(), 1);
        // Both strands reach past the cut, so AATT starts the fragment and ends it.
        let opened = &result.fragments[0].doc;
        assert_eq!(opened.len(), 64);
        assert_eq!(&opened.seq[..4], b"AATT");
        assert_eq!(&opened.seq[60..], b"AATT");
        assert!(ends_compatible(opened, opened));
        // The cut leaves a part of lacZ at each end.
        let parts: Vec<Span> = opened.features.iter().filter(|f| f.name == "lacZ").map(|f| f.segments[0].span).collect();
        assert_eq!(parts, vec![Span::new(0, 9), Span::new(54, 64)]);
        // Closing the circle gives the plasmid back, with lacZ whole again.
        let closed = circularize(opened).unwrap();
        assert!(crate::alphabet::same_circle(&closed.seq, &doc.seq));
        let lacz: Vec<&Feature> = closed.features.iter().filter(|f| f.name == "lacZ").collect();
        assert_eq!(lacz.len(), 1);
        assert_eq!(lacz[0].extract(&closed.seq), doc.features[1].extract(&doc.seq));
    }

    #[test]
    fn linear_digest_keeps_outer_ends() {
        let db = EnzymeDb::builtin();
        let mut seq = vec![b'C'; 30];
        seq[10..16].copy_from_slice(b"GGTACC"); // KpnI: GGTAC^C, 3' overhang
        let doc = Document::new_dna("lin", seq, Topology::Linear);
        let result = digest(&doc, &db.resolve(&["KpnI"]).unwrap()).unwrap();
        assert_eq!(result.fragments.len(), 2);
        assert_eq!(result.fragments[0].doc.ends.left, Overhang::Blunt);
        assert_eq!(result.fragments[0].doc.ends.right, Overhang::ThreePrime(4));
        assert_eq!(result.fragments[1].doc.ends.left, Overhang::ThreePrime(4));
        assert_eq!(result.fragments[0].len(), 15);
        assert_eq!(result.fragments[1].len(), 19);
        let back = ligate(&result.fragments[0].doc, &result.fragments[1].doc).unwrap();
        assert_eq!(back.seq, doc.seq);
    }
}
