//! Sequence edits that keep annotations in step.
//!
//! Each edit changes the bases and then moves every feature segment and primer binding site so
//! it still covers the same bases. The rules follow from where each surviving base ends up in
//! the new sequence:
//!
//! * Inserting between two bases of a span grows the span; inserting at its first base pushes
//!   it right; inserting after its last base leaves it alone.
//! * Deleting removes the deleted bases from every span; a span with nothing left is dropped.
//! * Replacing bases with the same number of bases (a mutation) moves nothing.
//!
//! Positions are 0-based. An insertion position `pos` means "before base `pos`", so `pos` can be
//! anything from 0 to the sequence length.

use crate::alphabet;
use crate::error::{Result, invalid};
use crate::model::*;
use crate::span::Span;

/// Where a span lands after inserting `n` bases at `pos` into a sequence of length `seq_len`.
fn map_span_insert(span: Span, pos: usize, n: usize, seq_len: usize, circular: bool) -> Span {
    if span.is_empty() {
        let start = if span.start >= pos { span.start + n } else { span.start };
        return Span::new(start, start);
    }
    // The insertion lands inside the span when the insertion point falls strictly between two
    // of its bases. On a circular molecule that includes the origin for spans that wrap.
    let inside = (span.start < pos && pos < span.end) || (circular && span.start < pos + seq_len && pos + seq_len < span.end);
    let start = if span.start >= pos { span.start + n } else { span.start };
    let len = span.len() + if inside { n } else { 0 };
    Span::new(start, start + len).normalize(seq_len + n)
}

/// Where a span lands after deleting `del` (which may wrap). `None` when nothing survives.
fn map_span_delete(span: Span, del: Span, seq_len: usize) -> Option<Span> {
    let overlap: usize = overlap_len(span, del, seq_len);
    let remaining = span.len().saturating_sub(overlap);
    if remaining == 0 {
        return None;
    }
    // The first surviving base: the span's own start, or the base just after the deletion
    // when the deletion swallowed the start.
    let first = if del.contains(span.start % seq_len, seq_len) { del.end % seq_len } else { span.start % seq_len };
    let new_len = seq_len - del.len();
    let new_first = map_index_delete(first, del, seq_len);
    Some(Span::new(new_first, new_first + remaining).normalize(new_len))
}

/// New index of a base that survives deletion of `del`.
fn map_index_delete(k: usize, del: Span, seq_len: usize) -> usize {
    if del.end <= seq_len {
        if k < del.start { k } else { k - del.len() }
    } else {
        // A wrapping deletion removes the tail and the head; survivors start at del.end - L.
        k - (del.end - seq_len)
    }
}

/// Number of bases two spans share.
fn overlap_len(a: Span, b: Span, seq_len: usize) -> usize {
    let mut total = 0;
    for (a0, a1) in a.pieces(seq_len) {
        for (b0, b1) in b.pieces(seq_len) {
            let lo = a0.max(b0);
            let hi = a1.min(b1);
            if hi > lo {
                total += hi - lo;
            }
        }
    }
    total
}

/// Apply a span mapping to every feature segment and primer site, dropping what disappears.
fn remap(doc: &mut Document, mut map: impl FnMut(Span) -> Option<Span>) {
    for feature in &mut doc.features {
        feature.segments.retain_mut(|segment| match map(segment.span) {
            Some(span) => {
                segment.span = span;
                true
            }
            None => false,
        });
    }
    doc.features.retain(|f| !f.segments.is_empty());
    for primer in &mut doc.primers {
        primer.sites.retain_mut(|site| match map(site.span) {
            Some(span) => {
                site.span = span;
                true
            }
            None => false,
        });
    }
}

fn check_bases(doc: &Document, bases: &[u8]) -> Result<()> {
    let ok = match doc.molecule {
        Molecule::Protein => bases.iter().all(|b| b.is_ascii_alphabetic() || *b == b'*'),
        _ => bases.iter().all(|&b| alphabet::is_nucleotide(b)),
    };
    if ok { Ok(()) } else { Err(invalid("bases contain characters that are not sequence codes")) }
}

/// Insert `bases` before base `pos` (0-based; `pos == len` appends).
pub fn insert(doc: &mut Document, pos: usize, bases: &[u8]) -> Result<()> {
    let seq_len = doc.len();
    if pos > seq_len {
        return Err(invalid(format!("insert position {} is past the end ({seq_len} bases)", pos + 1)));
    }
    if bases.is_empty() {
        return Ok(());
    }
    check_bases(doc, bases)?;
    let n = bases.len();
    let circular = doc.is_circular();
    remap(doc, |span| Some(map_span_insert(span, pos, n, seq_len, circular)));
    doc.seq.splice(pos..pos, bases.iter().copied());
    Ok(())
}

/// Delete the bases of `span` (which may run through the origin of a circular sequence).
pub fn delete(doc: &mut Document, span: Span) -> Result<()> {
    let seq_len = doc.len();
    let span = span.normalize(seq_len);
    if span.is_empty() {
        return Ok(());
    }
    if span.end > seq_len && !doc.is_circular() {
        return Err(invalid("a deletion on a linear sequence cannot run past the end"));
    }
    if span.len() >= seq_len {
        return Err(invalid("cannot delete the whole sequence"));
    }
    remap(doc, |s| map_span_delete(s, span, seq_len));
    if span.end <= seq_len {
        doc.seq.drain(span.start..span.end);
    } else {
        // Remove the tail first so the head indices stay valid.
        doc.seq.truncate(span.start);
        doc.seq.drain(..span.end - seq_len);
    }
    Ok(())
}

/// Replace the bases of `span` with `bases`. Equal lengths overwrite in place and move no
/// annotation; otherwise this is a delete followed by an insert at the same place.
pub fn replace(doc: &mut Document, span: Span, bases: &[u8]) -> Result<()> {
    let seq_len = doc.len();
    let span = span.normalize(seq_len);
    check_bases(doc, bases)?;
    if bases.len() == span.len() {
        for (ii, &b) in bases.iter().enumerate() {
            let pos = (span.start + ii) % seq_len;
            doc.seq[pos] = b;
        }
        return Ok(());
    }
    if span.is_empty() {
        return insert(doc, span.start, bases);
    }
    // Replacing every base (select all, then paste) leaves nothing for annotations to sit on, so
    // they go, as they would with the deleted bases. Primers keep their sequence but lose sites.
    if span.len() >= seq_len {
        if bases.is_empty() {
            return Err(invalid("cannot delete the whole sequence"));
        }
        remap(doc, |_| None);
        doc.seq = bases.to_vec();
        return Ok(());
    }
    // Where the replacement goes once the old bases are gone.
    let at = if span.end > seq_len { 0 } else { span.start };
    delete(doc, span)?;
    let at = at.min(doc.len());
    insert(doc, at, bases)
}

/// Reverse-complement the whole molecule. Annotations flip with it.
pub fn reverse_complement(doc: &mut Document) {
    let seq_len = doc.len();
    if seq_len == 0 {
        return;
    }
    // A span s..e maps to (L - e)..(L - s); written with modular arithmetic so spans through
    // the origin stay correct.
    let flip = |span: Span| -> Span {
        let start = (2 * seq_len - span.end) % seq_len;
        Span::new(start, start + span.len())
    };
    for feature in &mut doc.features {
        for segment in &mut feature.segments {
            segment.span = flip(segment.span);
        }
        feature.segments.reverse();
        feature.strand = feature.strand.flipped();
    }
    for primer in &mut doc.primers {
        for site in &mut primer.sites {
            site.span = flip(site.span);
            site.strand = site.strand.flipped();
        }
    }
    doc.seq = if doc.molecule == Molecule::Rna { alphabet::revcomp_rna(&doc.seq) } else { alphabet::revcomp(&doc.seq) };
    std::mem::swap(&mut doc.ends.left, &mut doc.ends.right);
    std::mem::swap(&mut doc.ends.left_modification, &mut doc.ends.right_modification);
}

/// Reverse-complement only the bases of `span`, in place (for example to flip an insert).
/// Annotations wholly inside the span flip with it; others stay where they are.
pub fn invert_region(doc: &mut Document, span: Span) -> Result<()> {
    let seq_len = doc.len();
    let span = span.normalize(seq_len);
    if span.wraps(seq_len) {
        return Err(invalid("invert a region that does not run through the origin (rotate the origin first)"));
    }
    let (rs, re) = (span.start, span.end);
    let flipped = alphabet::revcomp(&doc.seq[rs..re]);
    doc.seq[rs..re].copy_from_slice(&flipped);
    let mirror = |s: Span| Span::new(rs + re - s.end, rs + re - s.start);
    for feature in &mut doc.features {
        let extent = feature.span(seq_len);
        if extent.wraps(seq_len) || extent.start < rs || extent.end > re {
            continue;
        }
        for segment in &mut feature.segments {
            segment.span = mirror(segment.span);
        }
        feature.segments.reverse();
        feature.strand = feature.strand.flipped();
    }
    for primer in &mut doc.primers {
        for site in &mut primer.sites {
            if !site.span.wraps(seq_len) && site.span.start >= rs && site.span.end <= re {
                site.span = mirror(site.span);
                site.strand = site.strand.flipped();
            }
        }
    }
    Ok(())
}

/// Make base `new_origin` (0-based) the first base of a circular sequence.
pub fn set_origin(doc: &mut Document, new_origin: usize) -> Result<()> {
    if !doc.is_circular() {
        return Err(invalid("only circular sequences have a movable origin"));
    }
    let seq_len = doc.len();
    if new_origin >= seq_len {
        return Err(invalid(format!("position {} is past the end", new_origin + 1)));
    }
    if new_origin == 0 {
        return Ok(());
    }
    let shift = |span: Span| -> Span {
        let start = (span.start + seq_len - new_origin) % seq_len;
        Span::new(start, start + span.len())
    };
    for feature in &mut doc.features {
        for segment in &mut feature.segments {
            segment.span = shift(segment.span);
        }
        feature.normalize(seq_len, true);
    }
    for primer in &mut doc.primers {
        for site in &mut primer.sites {
            site.span = shift(site.span);
        }
    }
    doc.seq.rotate_left(new_origin);
    Ok(())
}

/// Change topology. Going from circular to linear splits annotations that run through the
/// origin into two pieces (a feature becomes a join of its two ends).
pub fn set_topology(doc: &mut Document, topology: Topology) {
    if doc.topology == topology {
        return;
    }
    let seq_len = doc.len();
    if topology == Topology::Linear {
        for feature in &mut doc.features {
            let mut pieces = Vec::new();
            for segment in &feature.segments {
                for (a, b) in segment.span.pieces(seq_len) {
                    let mut piece = segment.clone();
                    piece.span = Span::new(a, b);
                    pieces.push(piece);
                }
            }
            pieces.sort_by_key(|s| s.span.start);
            feature.segments = pieces;
        }
        for primer in &mut doc.primers {
            primer.sites.retain(|s| !s.span.wraps(seq_len));
        }
    }
    doc.topology = topology;
    if topology == Topology::Circular {
        doc.ends = Ends::default();
    }
}

/// A new linear document holding the bases of `span`, with the annotations that overlap it
/// (clipped to the region).
///
/// On a circular molecule the region may run through the origin, and it may be longer than the
/// molecule (up to two turns). A circle opened by a single cut that leaves a sticky end needs
/// that: both strands reach past the cut, so the overhang bases sit at both ends.
///
/// A feature that meets the region in two separate places becomes two features, one per part.
/// That is what a cut through a feature leaves when the region runs from the cut round to it.
pub fn extract(doc: &Document, span: Span) -> Result<Document> {
    let seq_len = doc.len();
    let circular = doc.is_circular();
    if span.is_empty() || seq_len == 0 {
        return Err(invalid("cannot extract an empty region"));
    }
    // Bring the start into the sequence. A circular region keeps its length, up to two turns.
    let span = if circular {
        let start = span.start % seq_len;
        Span::new(start, start + span.len().min(2 * seq_len))
    } else {
        span.normalize(seq_len)
    };
    if span.end > seq_len && !circular {
        return Err(invalid("region runs past the end of a linear sequence"));
    }
    let mut out = Document {
        name: format!("{} {}", doc.name, span.format(seq_len)).trim().to_string(),
        molecule: doc.molecule,
        topology: Topology::Linear,
        strandedness: doc.strandedness,
        // The bases, read round the circle as often as the region needs.
        seq: (span.start..span.end).map(|ii| doc.seq[ii % seq_len]).collect(),
        methylation: doc.methylation,
        ..Document::default()
    };
    out.meta.synthetic = doc.meta.synthetic;
    out.meta.organism = doc.meta.organism.clone();

    // Map a span of the parent into region coordinates: every part of it inside the region,
    // clipped to the region. On a circle the span is tried one turn back, as written, and one
    // and two turns on, so the parts through the origin and the repeats of a region longer than
    // the molecule are all found. Parts that touch are one part.
    let turns: &[i64] = if circular { &[-1, 0, 1, 2] } else { &[0] };
    let clip = |s: Span| -> Vec<Span> {
        let (region_start, region_end) = (span.start as i64, span.end as i64);
        let mut parts: Vec<Span> = Vec::new();
        for turn in turns {
            let lo = (s.start as i64 + turn * seq_len as i64).max(region_start);
            let hi = (s.end as i64 + turn * seq_len as i64).min(region_end);
            if hi > lo {
                parts.push(Span::new((lo - region_start) as usize, (hi - region_start) as usize));
            }
        }
        parts.sort_by_key(|p| p.start);
        let mut joined: Vec<Span> = Vec::new();
        for part in parts {
            match joined.last_mut() {
                Some(last) if part.start <= last.end => last.end = last.end.max(part.end),
                _ => joined.push(part),
            }
        }
        joined
    };
    for feature in &doc.features {
        let parts: Vec<Vec<Span>> = feature.segments.iter().map(|segment| clip(segment.span)).collect();
        // A one-piece feature that meets the region twice was cut through: keep each part as a
        // feature of its own, as a digest does when the cut leaves the parts in two fragments.
        if feature.segments.len() == 1 && parts[0].len() > 1 {
            for part in &parts[0] {
                let mut copy = feature.clone();
                copy.uid = new_uid();
                copy.segments[0].span = *part;
                out.features.push(copy);
            }
            continue;
        }
        // Otherwise the feature keeps every part of every segment, in order along the region.
        let mut copy = feature.clone();
        copy.uid = new_uid();
        copy.segments = feature
            .segments
            .iter()
            .zip(&parts)
            .flat_map(|(segment, spans)| {
                spans.iter().map(move |part| {
                    let mut s = segment.clone();
                    s.span = *part;
                    s
                })
            })
            .collect();
        if copy.segments.is_empty() {
            continue;
        }
        copy.segments.sort_by_key(|s| s.span.start);
        out.features.push(copy);
    }
    for primer in &doc.primers {
        // A binding site comes along only where the region holds all of it.
        let sites: Vec<BindingSite> = primer
            .sites
            .iter()
            .flat_map(|site| clip(site.span).into_iter().filter(|part| part.len() == site.span.len()).map(|part| BindingSite { span: part, ..site.clone() }))
            .collect();
        if !sites.is_empty() {
            let mut copy = primer.clone();
            copy.uid = new_uid();
            copy.sites = sites;
            out.primers.push(copy);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc_with(seq: &[u8], circular: bool, spans: &[Span]) -> Document {
        let topology = if circular { Topology::Circular } else { Topology::Linear };
        let mut doc = Document::new_dna("t", seq, topology);
        for (ii, &span) in spans.iter().enumerate() {
            doc.features.push(Feature::new(format!("f{ii}"), "misc_feature", span, Strand::Forward));
        }
        doc
    }

    /// Check a feature still covers exactly the expected bases.
    fn feature_bases(doc: &Document, index: usize) -> Vec<u8> {
        doc.features[index].extract(&doc.seq)
    }

    #[test]
    fn insert_inside_grows_and_before_shifts() {
        let mut doc = doc_with(b"AAAACCCCGGGG", false, &[Span::new(4, 8), Span::new(8, 12)]);
        insert(&mut doc, 6, b"TT").unwrap();
        assert_eq!(doc.seq, b"AAAACCTTCCGGGG".to_vec());
        assert_eq!(feature_bases(&doc, 0), b"CCTTCC".to_vec());
        assert_eq!(feature_bases(&doc, 1), b"GGGG".to_vec());
        // Insert exactly at a feature start: the feature moves right, unchanged.
        insert(&mut doc, 10, b"A").unwrap();
        assert_eq!(feature_bases(&doc, 1), b"GGGG".to_vec());
        assert_eq!(feature_bases(&doc, 0), b"CCTTCC".to_vec());
    }

    #[test]
    fn insert_at_origin_inside_wrapping_feature() {
        // Feature covers GG + AA across the origin of a circular molecule.
        let mut doc = doc_with(b"AACCCCGG", true, &[Span::new(6, 10)]);
        insert(&mut doc, 0, b"T").unwrap();
        assert_eq!(doc.seq, b"TAACCCCGG".to_vec());
        assert_eq!(feature_bases(&doc, 0), b"GGTAA".to_vec());
        let end = doc.len();
        insert(&mut doc, end, b"C").unwrap();
        assert_eq!(feature_bases(&doc, 0), b"GGCTAA".to_vec());
    }

    #[test]
    fn delete_trims_and_drops() {
        let mut doc = doc_with(b"AAAACCCCGGGGTTTT", false, &[Span::new(2, 10), Span::new(4, 8), Span::new(12, 16)]);
        delete(&mut doc, Span::new(4, 8)).unwrap();
        assert_eq!(doc.seq, b"AAAAGGGGTTTT".to_vec());
        assert_eq!(doc.features.len(), 2, "the feature made only of deleted bases is gone");
        assert_eq!(feature_bases(&doc, 0), b"AAGG".to_vec());
        assert_eq!(feature_bases(&doc, 1), b"TTTT".to_vec());
    }

    #[test]
    fn delete_through_origin() {
        let mut doc = doc_with(b"ABCDEFGHIJ", true, &[Span::new(8, 13), Span::new(4, 6)]);
        delete(&mut doc, Span::new(9, 11)).unwrap(); // removes J and A
        assert_eq!(doc.seq, b"BCDEFGHI".to_vec());
        assert_eq!(feature_bases(&doc, 0), b"IBC".to_vec());
        assert_eq!(feature_bases(&doc, 1), b"EF".to_vec());
    }

    #[test]
    fn mutation_moves_nothing() {
        let mut doc = doc_with(b"AAAACCCC", false, &[Span::new(2, 6)]);
        replace(&mut doc, Span::new(3, 5), b"GG").unwrap();
        assert_eq!(doc.seq, b"AAAGGCCC".to_vec());
        assert_eq!(doc.features[0].segments[0].span, Span::new(2, 6));
    }

    #[test]
    fn revcomp_flips_features() {
        let mut doc = doc_with(b"AAAACCCGGT", true, &[Span::new(8, 12)]);
        let before = feature_bases(&doc, 0);
        reverse_complement(&mut doc);
        assert_eq!(doc.features[0].strand, Strand::Reverse);
        assert_eq!(feature_bases(&doc, 0), before, "a feature reads the same along its own strand");
    }

    #[test]
    fn origin_rotation_keeps_feature_bases() {
        let mut doc = doc_with(b"ABCDEFGHIJ", true, &[Span::new(1, 4), Span::new(8, 12)]);
        let before: Vec<Vec<u8>> = (0..2).map(|ii| feature_bases(&doc, ii)).collect();
        set_origin(&mut doc, 3).unwrap();
        assert_eq!(doc.seq, b"DEFGHIJABC".to_vec());
        for (ii, expected) in before.iter().enumerate() {
            assert_eq!(&feature_bases(&doc, ii), expected);
        }
    }

    #[test]
    fn extract_clips_features() {
        let doc = doc_with(b"AAAACCCCGGGGTTTT", false, &[Span::new(2, 6), Span::new(10, 14)]);
        let part = extract(&doc, Span::new(4, 12)).unwrap();
        assert_eq!(part.seq, b"CCCCGGGG".to_vec());
        assert_eq!(part.features.len(), 2);
        assert_eq!(part.features[0].segments[0].span, Span::new(0, 2));
        assert_eq!(part.features[1].segments[0].span, Span::new(6, 8));
    }

    #[test]
    fn replacing_every_base_swaps_the_sequence() {
        // Select all, then paste: the new bases replace the old ones and the old features go.
        let mut doc = doc_with(b"AAAACCCCGGGG", true, &[Span::new(2, 6)]);
        doc.primers.push(Primer::new("p", "AAAACC"));
        replace(&mut doc, Span::new(0, 12), b"TTTGGG").unwrap();
        assert_eq!(doc.seq, b"TTTGGG".to_vec());
        assert!(doc.features.is_empty());
        assert_eq!(doc.primers.len(), 1, "primers stay; only their sites go");
        // Replacing everything with nothing would leave an empty molecule, as a deletion would.
        assert!(replace(&mut doc, Span::new(0, 6), b"").is_err());
        assert_eq!(doc.seq, b"TTTGGG".to_vec());
    }

    #[test]
    fn extract_through_origin() {
        let doc = doc_with(b"ABCDEFGHIJ", true, &[Span::new(8, 12)]);
        let part = extract(&doc, Span::new(7, 13)).unwrap();
        assert_eq!(part.seq, b"HIJABC".to_vec());
        assert_eq!(part.features[0].segments[0].span, Span::new(1, 5));
    }

    #[test]
    fn extract_keeps_both_parts_of_a_feature_it_meets_twice() {
        // The region runs from base 6 round the origin to base 4, and the feature covers 2..8,
        // so the region holds both ends of the feature but not its middle.
        let doc = doc_with(b"ABCDEFGHIJ", true, &[Span::new(2, 8)]);
        let part = extract(&doc, Span::new(6, 14)).unwrap();
        assert_eq!(part.seq, b"GHIJABCD".to_vec());
        let spans: Vec<Span> = part.features.iter().map(|f| f.segments[0].span).collect();
        assert_eq!(spans, vec![Span::new(0, 2), Span::new(6, 8)]);
        // A feature through the origin is found by a region that does not run through it.
        let doc = doc_with(b"ABCDEFGHIJ", true, &[Span::new(8, 13)]);
        let part = extract(&doc, Span::new(1, 5)).unwrap();
        assert_eq!(part.features[0].segments[0].span, Span::new(0, 2));
        // A circular region may be longer than the molecule, and then repeats its first bases.
        let part = extract(&doc, Span::new(3, 15)).unwrap();
        assert_eq!(part.seq, b"DEFGHIJABCDE".to_vec());
    }
}
