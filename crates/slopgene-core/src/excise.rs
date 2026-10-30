//! Which enzymes cut a chosen region out of a molecule.
//!
//! A person marks the region on the map, as a feature or a stretch of bases. [`nearest_cuts`]
//! lists the enzymes that cut closest to it on each side, nearest first. Each enzyme appears once
//! on each side, at its cut nearest the region, with the end that cut leaves and how often the
//! enzyme cuts the molecule. An enzyme that cuts inside the region is listed apart, because it
//! would cut the region too. Nothing is paired or scored: the person picks one enzyme on each side.
//!
//! [`excision`] works out what a pick gives. It digests the molecule with the two enzymes and
//! finds the piece that holds the region, so the piece's ends, its length and the bases it carries
//! beyond the region are exact. With a vector it also says whether the two enzymes open the vector
//! with ends that take the piece.

use crate::digest::{digest, ends_compatible, overhang_bases};
use crate::enzymes::{Cut, Enzyme, cutting_sites_by_enzyme, find_sites_in};
use crate::error::{Result, invalid};
use crate::model::*;
use crate::render::thousands;
use crate::span::Span;
use serde_json::{Value, json};
use std::collections::HashMap;

/// How many enzymes [`nearest_cuts`] lists on each side unless asked otherwise.
pub const DEFAULT_PER_SIDE: usize = 8;

/// One enzyme's cut nearest the region on one side of it.
#[derive(Clone, Debug, PartialEq)]
pub struct NearCut {
    pub enzyme: String,
    /// The top-strand cut position: the number of bases before the cut, as maps label it.
    pub pos: usize,
    /// Top-strand bases between the cut and the region, 0 for a cut at the region's edge.
    pub distance: usize,
    /// The end this cut leaves, and its overhang bases read on the top strand.
    pub end: Overhang,
    pub bases: String,
    /// How many times the enzyme cuts the molecule in all.
    pub cuts: usize,
    /// How many times it cuts the vector, when a vector was given.
    pub vector_cuts: Option<usize>,
}

impl NearCut {
    /// The end in words: "5' CATG", "3' ACGT" or "blunt".
    pub fn end_text(&self) -> String {
        Excision::end_text(self.end, &self.bases)
    }

    /// The cut as JSON, with its 1-based position label (bases before the cut).
    pub fn to_json(&self) -> Value {
        json!({
            "enzyme": self.enzyme,
            "cut": self.pos,
            "distance": self.distance,
            "end": self.end_text(),
            "cuts": self.cuts,
            "vector_cuts": self.vector_cuts,
        })
    }
}

/// What [`nearest_cuts`] found.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NearestCuts {
    /// The enzymes that cut before the region, nearest first.
    pub before: Vec<NearCut>,
    /// The enzymes that cut after the region, nearest first.
    pub after: Vec<NearCut>,
    /// The enzymes that cut inside the region, which would cut it too.
    pub inside: Vec<String>,
}

impl NearestCuts {
    pub fn to_json(&self) -> Value {
        json!({
            "before": self.before.iter().map(NearCut::to_json).collect::<Vec<_>>(),
            "after": self.after.iter().map(NearCut::to_json).collect::<Vec<_>>(),
            "inside": self.inside,
        })
    }
}

/// A cut that bounds the fragment on one side.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FlankCut {
    pub enzyme: String,
    /// The top-strand cut position: the number of bases before the cut, as maps label it.
    pub pos: usize,
    /// Top-strand bases between the cut and the region, which the fragment carries along.
    pub extra: usize,
}

/// How a pair of enzymes opens the vector.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VectorFit {
    /// Each enzyme and its cut positions in the vector, for the enzymes that cut it few enough times.
    pub cuts: Vec<(String, usize)>,
    /// The top-strand bases between the two cuts, which the vector loses (0 for one cut).
    pub dropped: usize,
    /// Why the pair cannot open the vector for this insert, if it cannot.
    pub problem: Option<String>,
}

impl VectorFit {
    /// The fit in words, such as "cuts once each, at 396 and 417, dropping 21 bp".
    pub fn describe(&self) -> String {
        if let Some(problem) = &self.problem {
            return problem.clone();
        }
        let positions: Vec<String> = self.cuts.iter().map(|(_, p)| thousands(*p)).collect();
        match self.cuts.as_slice() {
            [_] => format!("cuts once, at {}", positions[0]),
            [(a, _), (b, _)] if a == b => format!("cuts twice, at {}, dropping {} bp", positions.join(" and "), thousands(self.dropped)),
            _ => format!("cuts once each, at {}, dropping {} bp", positions.join(" and "), thousands(self.dropped)),
        }
    }
}

/// What cutting with a chosen pair of enzymes gives: the piece that holds the region.
#[derive(Clone, Debug, PartialEq)]
pub struct Excision {
    pub left: FlankCut,
    pub right: FlankCut,
    /// The fragment over the donor's sequence, both strands together.
    pub fragment: Span,
    /// The fragment's length, both strands together.
    pub length: usize,
    pub left_end: Overhang,
    pub right_end: Overhang,
    /// The overhang bases of each end, read on the top strand.
    pub left_bases: String,
    pub right_bases: String,
    /// Whether the two ends differ, so the insert ligates one way round only.
    pub directional: bool,
    /// Cuts of the two enzymes elsewhere in the donor. The fragment survives them; they add
    /// bands to the gel.
    pub other_cuts: usize,
    /// How the enzymes open the vector, when one was given.
    pub vector: Option<VectorFit>,
}

impl Excision {
    /// The enzyme names, one when the same enzyme makes both ends.
    pub fn enzymes(&self) -> Vec<String> {
        if self.left.enzyme.eq_ignore_ascii_case(&self.right.enzyme) {
            vec![self.left.enzyme.clone()]
        } else {
            vec![self.left.enzyme.clone(), self.right.enzyme.clone()]
        }
    }

    /// The enzymes as a label: "EcoRI / BamHI", or "EcoRI" when it makes both ends.
    pub fn label(&self) -> String {
        self.enzymes().join(" / ")
    }

    /// One end in words: "5' AATT", "3' GTAC" or "blunt".
    pub fn end_text(end: Overhang, bases: &str) -> String {
        match end {
            Overhang::Blunt => "blunt".to_string(),
            Overhang::FivePrime(_) => format!("5' {bases}"),
            Overhang::ThreePrime(_) => format!("3' {bases}"),
        }
    }

    /// Both ends in words: "5' AATT / 5' GATC".
    pub fn ends_text(&self) -> String {
        format!("{} / {}", Excision::end_text(self.left_end, &self.left_bases), Excision::end_text(self.right_end, &self.right_bases))
    }

    /// Whether the insert goes in one way round or either way.
    pub fn direction_text(&self) -> &'static str {
        if self.directional { "one way round" } else { "either way round" }
    }

    /// The bases carried along in words: "+4 / +8 bp", or "exact" for none.
    pub fn extra_text(&self) -> String {
        if self.left.extra == 0 && self.right.extra == 0 { "exact".to_string() } else { format!("+{} / +{} bp", self.left.extra, self.right.extra) }
    }

    /// The excision on one line, for listings and tooltips.
    pub fn summary(&self, donor_len: usize) -> String {
        let mut text = format!(
            "{}: {} ({} bp), ends {}, {}, {}",
            self.label(),
            self.fragment.format(donor_len),
            thousands(self.length),
            self.ends_text(),
            self.direction_text(),
            self.extra_text()
        );
        if self.other_cuts > 0 {
            text.push_str(&format!(", {} other cut(s) in the donor", self.other_cuts));
        }
        if let Some(fit) = &self.vector {
            text.push_str(&format!("; vector: {}", fit.describe()));
        }
        text
    }

    /// The excision as JSON, with 1-based positions.
    pub fn to_json(&self, donor_len: usize) -> Value {
        let (start, end) = self.fragment.to_1based(donor_len);
        let cut =
            |c: &FlankCut, end: Overhang, bases: &str| json!({ "enzyme": c.enzyme, "cut": c.pos, "extra": c.extra, "end": Excision::end_text(end, bases) });
        let mut value = json!({
            "enzymes": self.enzymes(),
            "left": cut(&self.left, self.left_end, &self.left_bases),
            "right": cut(&self.right, self.right_end, &self.right_bases),
            "start": start,
            "end": end,
            "length": self.length,
            "directional": self.directional,
            "other_cuts": self.other_cuts,
        });
        if let Some(fit) = &self.vector {
            value["vector"] = json!({
                "cuts": fit.cuts.iter().map(|(e, p)| json!({ "enzyme": e, "cut": p })).collect::<Vec<_>>(),
                "dropped": fit.dropped,
                "problem": fit.problem,
            });
        }
        value
    }
}

/// A cut's distance to a span: `None` when the cut falls inside the span, and otherwise how
/// far it lies before the span's start and after its end. On a circle every outside cut lies
/// both before and after the span; on a linear molecule only one side applies. A cut at the
/// span's boundary is 0 bases before (or after) it.
fn distances(span: Span, pos: usize, len: usize, circular: bool) -> Option<(Option<usize>, Option<usize>)> {
    let span_len = span.len();
    if circular {
        let rel = (pos % len + len - span.start % len) % len;
        if rel > 0 && rel < span_len {
            return None;
        }
        let before = (len - rel) % len;
        let after = (rel + len - span_len) % len;
        return Some((Some(before), Some(after)));
    }
    if pos > span.start && pos < span.end {
        return None;
    }
    let before = (pos <= span.start).then(|| span.start - pos);
    let after = (pos >= span.end).then(|| pos - span.end);
    Some((before, after))
}

/// Where a cut's left and right boundaries fall in the sequence, both strands together: the
/// earlier strand's cut starts a fragment and the later strand's cut ends one.
fn cut_bounds(cut: &Cut, len: usize) -> (usize, usize) {
    let start = (cut.pos as i64 + i64::from(cut.overhang.min(0))).rem_euclid(len as i64) as usize;
    let end = (cut.pos as i64 + i64::from(cut.overhang.max(0))).rem_euclid(len as i64) as usize;
    (start, end)
}

/// The end one cut leaves, and its overhang bases read on the top strand. The overhang runs
/// between the two strands' cuts, which a circle may place across its origin.
fn cut_end(doc: &Document, cut: &Cut) -> (Overhang, String) {
    let end = match cut.overhang {
        0 => Overhang::Blunt,
        o if o > 0 => Overhang::FivePrime(o as usize),
        o => Overhang::ThreePrime(o.unsigned_abs() as usize),
    };
    let len = doc.len() as i64;
    let from = cut.pos as i64 + i64::from(cut.overhang.min(0));
    let bases: String = (0..end.len() as i64)
        .filter_map(|ii| {
            let at = from + ii;
            // A linear molecule has no bases past its ends; a circle continues through its origin.
            let at = if doc.is_circular() { at.rem_euclid(len) } else { at };
            (0..len).contains(&at).then(|| doc.seq[at as usize].to_ascii_uppercase() as char)
        })
        .collect();
    (end, bases)
}

/// Every distinct cut of each enzyme in a molecule, by enzyme name, in position order.
fn cuts_by_enzyme(doc: &Document, enzymes: &[&Enzyme]) -> HashMap<String, Vec<Cut>> {
    let sites = find_sites_in(doc, enzymes);
    let mut out: HashMap<String, Vec<Cut>> = HashMap::new();
    for (name, group) in cutting_sites_by_enzyme(&sites) {
        let mut cuts: Vec<Cut> = group.iter().flat_map(|s| s.cuts.iter().copied()).collect();
        cuts.sort();
        cuts.dedup();
        out.insert(name.to_string(), cuts);
    }
    out
}

/// Whether an end with this overhang can join one of `piece`'s ends, with the piece either way
/// round: the same kind and length of overhang, with bases that pair. The bases are read on
/// the top strand, as [`NearCut::bases`] and [`overhang_bases`] give them.
pub fn can_join(end: Overhang, bases: &str, piece: &Document) -> bool {
    if piece.is_circular() {
        return false;
    }
    let mine = bases.to_ascii_uppercase();
    let kind = |o: Overhang| std::mem::discriminant(&o);
    [(piece.ends.left, true), (piece.ends.right, false)].into_iter().any(|(other, left)| {
        let theirs = overhang_bases(piece, left);
        let turned = String::from_utf8_lossy(&crate::alphabet::revcomp(&theirs)).into_owned();
        kind(other) == kind(end) && other.len() == end.len() && (mine == String::from_utf8_lossy(&theirs) || mine == turned)
    })
}

/// The enzymes among `enzymes` that cut closest to `region` on each side, at most `per_side` on
/// each, nearest first. With a `vector`, each also says how often it cuts the vector.
pub fn nearest_cuts(doc: &Document, region: Span, enzymes: &[&Enzyme], vector: Option<&Document>, per_side: usize) -> NearestCuts {
    let mut found = NearestCuts::default();
    let len = doc.len();
    let circular = doc.is_circular();
    if len == 0 || !doc.molecule.is_nucleic() || region.is_empty() || region.len() >= len {
        return found;
    }
    let region = region.normalize(len);
    let in_donor = cuts_by_enzyme(doc, enzymes);
    let in_vector = vector.filter(|v| v.molecule.is_nucleic()).map(|v| cuts_by_enzyme(v, enzymes));

    // Each enzyme is placed by its cuts: set apart when one falls inside the region, and
    // otherwise listed on each side at its cut nearest the region there.
    for enzyme in enzymes {
        let Some(cuts) = in_donor.get(&enzyme.name) else { continue };
        let placed: Vec<Option<(Option<usize>, Option<usize>)>> = cuts.iter().map(|c| distances(region, c.pos, len, circular)).collect();
        if placed.iter().any(Option::is_none) {
            if !found.inside.contains(&enzyme.name) {
                found.inside.push(enzyme.name.clone());
            }
            continue;
        }
        // (distance, cut) of the nearest cut before the region and after it.
        let mut nearest: [Option<(usize, Cut)>; 2] = [None, None];
        for (cut, place) in cuts.iter().zip(&placed) {
            let Some((before, after)) = place else { continue };
            for (slot, distance) in nearest.iter_mut().zip([before, after]) {
                if let Some(d) = distance
                    && slot.is_none_or(|(best, _)| *d < best)
                {
                    *slot = Some((*d, *cut));
                }
            }
        }
        let vector_cuts = in_vector.as_ref().map(|counts| counts.get(&enzyme.name).map_or(0, Vec::len));
        for (side, slot) in [&mut found.before, &mut found.after].into_iter().zip(nearest) {
            let Some((distance, cut)) = slot else { continue };
            if side.iter().any(|c| c.enzyme == enzyme.name) {
                continue;
            }
            let (end, bases) = cut_end(doc, &cut);
            side.push(NearCut { enzyme: enzyme.name.clone(), pos: cut.pos, distance, end, bases, cuts: cuts.len(), vector_cuts });
        }
    }

    // Nearest first on each side, with ties in name order so the list is the same every time.
    for side in [&mut found.before, &mut found.after] {
        side.sort_by(|a, b| a.distance.cmp(&b.distance).then_with(|| a.enzyme.to_ascii_lowercase().cmp(&b.enzyme.to_ascii_lowercase())));
        side.truncate(per_side);
    }
    found.inside.sort_by_key(|name| name.to_ascii_lowercase());
    found
}

/// What cutting `doc` with `left` before `region` and `right` after it gives: the piece that
/// holds the region, with its ends, and with a `vector` whether the two enzymes open it for the
/// piece. The same enzyme may make both ends. An enzyme that does not cut, or cuts inside the
/// region, or a pair that leaves no piece holding all of the region, is an error that says so.
pub fn excision(doc: &Document, region: Span, left: &Enzyme, right: &Enzyme, vector: Option<&Document>) -> Result<Excision> {
    let len = doc.len();
    let circular = doc.is_circular();
    if len == 0 || !doc.molecule.is_nucleic() {
        return Err(invalid("only DNA can be cut"));
    }
    if region.is_empty() || region.len() >= len {
        return Err(invalid("mark a region smaller than the whole molecule"));
    }
    let region = region.normalize(len);
    let same = left.name.eq_ignore_ascii_case(&right.name);
    let pair: Vec<&Enzyme> = if same { vec![left] } else { vec![left, right] };

    // Each enzyme must cut, and must leave the region whole.
    let cuts = cuts_by_enzyme(doc, &pair);
    for enzyme in &pair {
        let Some(positions) = cuts.get(&enzyme.name) else {
            return Err(invalid(format!("{} does not cut {}", enzyme.name, doc.name)));
        };
        if positions.iter().any(|c| distances(region, c.pos, len, circular).is_none()) {
            return Err(invalid(format!("{} cuts inside the region, so it would cut the region too", enzyme.name)));
        }
    }

    // One cut in a circle only opens it: the one piece is the whole molecule.
    let mut positions: Vec<usize> = cuts.values().flatten().map(|c| c.pos).collect();
    positions.sort_unstable();
    positions.dedup();
    if circular && positions.len() == 1 {
        return Err(invalid(format!("{} cuts {} only once, which opens the circle without cutting the region out", left.name, doc.name)));
    }

    // The digest gives the exact piece. It must have enzyme-made ends on both sides and hold the
    // whole region.
    let cut = digest(doc, &pair)?;
    let Some(fragment) =
        cut.fragments.iter().find(|f| f.left_enzyme.is_some() && f.right_enzyme.is_some() && !f.doc.is_circular() && f.span.contains_span(&region, len))
    else {
        return Err(invalid("no piece holds the whole region: choose an enzyme that cuts before it and one that cuts after it"));
    };

    // The cuts that bound the piece, found among the digest's cuts by the boundary each one
    // makes. They may differ from the pair asked for, when one enzyme also cuts nearer the
    // region on the other enzyme's side.
    let bounding = |name: &str, want: usize, at_start: bool| -> Option<usize> {
        cut.sites
            .iter()
            .filter(|s| s.enzyme == name)
            .flat_map(|s| s.cuts.iter())
            .find(|c| {
                let (start, end) = cut_bounds(c, len);
                if at_start { start == want } else { end == want }
            })
            .map(|c| c.pos)
    };
    let (Some(left_name), Some(right_name)) = (fragment.left_enzyme.as_deref(), fragment.right_enzyme.as_deref()) else {
        return Err(invalid("the piece that holds the region has an end no enzyme made"));
    };
    let (Some(left_pos), Some(right_pos)) = (bounding(left_name, fragment.span.start, true), bounding(right_name, fragment.span.end % len, false)) else {
        return Err(invalid("the cuts that bound the piece could not be found"));
    };
    let (Some((Some(extra_left), _)), Some((_, Some(extra_right)))) = (distances(region, left_pos, len, circular), distances(region, right_pos, len, circular))
    else {
        return Err(invalid("the piece's cuts do not lie on either side of the region"));
    };

    // Every distinct cut of the bounding enzymes beyond the two that bound the piece adds a band
    // to the gel.
    let bounding_pair: Vec<&Enzyme> = pair.iter().copied().filter(|e| e.name == left_name || e.name == right_name).collect();
    let mut all_cuts: Vec<Cut> = cut.sites.iter().filter(|s| s.enzyme == left_name || s.enzyme == right_name).flat_map(|s| s.cuts.iter().copied()).collect();
    all_cuts.sort();
    all_cuts.dedup();
    Ok(Excision {
        left: FlankCut { enzyme: left_name.to_string(), pos: left_pos, extra: extra_left },
        right: FlankCut { enzyme: right_name.to_string(), pos: right_pos, extra: extra_right },
        fragment: fragment.span,
        length: fragment.len(),
        left_end: fragment.doc.ends.left,
        right_end: fragment.doc.ends.right,
        left_bases: String::from_utf8_lossy(&overhang_bases(&fragment.doc, true)).into_owned(),
        right_bases: String::from_utf8_lossy(&overhang_bases(&fragment.doc, false)).into_owned(),
        directional: !ends_compatible(&fragment.doc, &fragment.doc),
        other_cuts: all_cuts.len().saturating_sub(2),
        vector: vector.map(|v| vector_fit(v, &bounding_pair, &fragment.doc)),
    })
}

/// How `pair` opens the vector for `insert`: two enzymes must cut once each, and one enzyme
/// once or twice (twice replaces the bases between its sites). The backbone's ends must take
/// the insert one way or the other.
fn vector_fit(vector: &Document, pair: &[&Enzyme], insert: &Document) -> VectorFit {
    let len = vector.len();
    let circular = vector.is_circular();
    let found = cuts_by_enzyme(vector, pair);
    let mut cuts = Vec::new();
    let mut problems = Vec::new();
    let allowed = if pair.len() == 1 { 2 } else { 1 };
    for enzyme in pair {
        let positions: Vec<usize> = found.get(&enzyme.name).map(|c| c.iter().map(|c| c.pos).collect()).unwrap_or_default();
        match positions.len() {
            0 => problems.push(format!("{} does not cut the vector", enzyme.name)),
            n if n > allowed => problems.push(format!("{} cuts the vector {n} times", enzyme.name)),
            _ => cuts.extend(positions.iter().map(|p| (enzyme.name.clone(), *p))),
        }
    }
    // The ends of the backbone (the largest piece of the opened vector) must pair with the
    // insert's. Enzymes with variable sites can leave overhangs that differ between the maps.
    if problems.is_empty()
        && let Ok(opened) = digest(vector, pair)
        && let Some(backbone) = opened.fragments.iter().max_by_key(|f| f.len())
    {
        let mut flipped = insert.clone();
        crate::edit::reverse_complement(&mut flipped);
        let takes = |piece: &Document| ends_compatible(&backbone.doc, piece) && ends_compatible(piece, &backbone.doc);
        if !takes(insert) && !takes(&flipped) {
            problems.push("the ends it leaves on the vector do not pair with the insert's".to_string());
        }
    }
    let dropped = match cuts.as_slice() {
        [(_, a), (_, b)] => {
            let forward = if circular { (b + len - a) % len } else { b.abs_diff(*a) };
            if circular { forward.min(len - forward) } else { forward }
        }
        _ => 0,
    };
    VectorFit { cuts, dropped, problem: (!problems.is_empty()).then(|| problems.join("; ")) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzymes::EnzymeDb;

    /// A background of A and C only, in a pseudo-random order. Every common enzyme's site holds
    /// a G or a T, so only the sites planted into it cut, and a planted site cannot complete
    /// another site with the bases beside it.
    fn clean(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"AC"[(x >> 33) as usize % 2]
            })
            .collect()
    }

    /// The common enzymes that cut a document, sorted, so a fixture can say what it holds.
    fn cutters(doc: &Document) -> Vec<String> {
        let mut names: Vec<String> = crate::enzymes::cut_counts_in(doc, &common()).into_iter().filter(|(_, n)| *n > 0).map(|(e, _)| e.name.clone()).collect();
        names.sort();
        names
    }

    /// A 3,000 bp circular donor with geneX at 1000..1500, an EcoRI cut 4 bp before it, a
    /// BamHI cut 8 bp after it, and a second BamHI site far away.
    fn donor() -> Document {
        let mut seq = clean(3000, 1);
        // EcoRI cuts G^AATTC, so a site at 995 cuts at 996, 4 bases before the region.
        seq[995..1001].copy_from_slice(b"GAATTC");
        // BamHI cuts G^GATCC, so a site at 1507 cuts at 1508, 8 bases after the region.
        seq[1507..1513].copy_from_slice(b"GGATCC");
        seq[2400..2406].copy_from_slice(b"GGATCC");
        let mut doc = Document::new_dna("donor", seq, Topology::Circular);
        doc.features.push(Feature::new("geneX", "CDS", Span::new(1000, 1500), Strand::Forward));
        assert_eq!(cutters(&doc), ["BamHI", "EcoRI"]);
        doc
    }

    fn common() -> Vec<&'static Enzyme> {
        EnzymeDb::builtin().set("common").unwrap()
    }

    fn enzyme(name: &str) -> &'static Enzyme {
        EnzymeDb::builtin().get(name).unwrap()
    }

    fn region() -> Span {
        Span::new(1000, 1500)
    }

    #[test]
    fn distances_count_from_the_boundaries_either_way_round_a_circle() {
        let span = Span::new(100, 200);
        assert_eq!(distances(span, 150, 1000, true), None, "inside");
        assert_eq!(distances(span, 100, 1000, true), Some((Some(0), Some(900))), "at the start");
        assert_eq!(distances(span, 200, 1000, true), Some((Some(900), Some(0))), "at the end");
        assert_eq!(distances(span, 96, 1000, true), Some((Some(4), Some(896))));
        assert_eq!(distances(span, 208, 1000, true), Some((Some(892), Some(8))));
        // A span through the origin.
        let wrapped = Span::new(950, 1050);
        assert_eq!(distances(wrapped, 20, 1000, true), None);
        assert_eq!(distances(wrapped, 946, 1000, true), Some((Some(4), Some(896))));
        assert_eq!(distances(wrapped, 58, 1000, true), Some((Some(892), Some(8))));
        // On a linear molecule a cut lies on one side only.
        assert_eq!(distances(span, 96, 1000, false), Some((Some(4), None)));
        assert_eq!(distances(span, 208, 1000, false), Some((None, Some(8))));
        assert_eq!(distances(span, 150, 1000, false), None);
    }

    #[test]
    fn the_nearest_cuts_on_each_side_come_first() {
        let found = nearest_cuts(&donor(), region(), &common(), None, DEFAULT_PER_SIDE);
        // Before the region: EcoRI 4 bp away, then BamHI's nearest cut that way round the circle
        // (its site at 2,400 lies 1,599 bp before the region's start).
        let before: Vec<(&str, usize)> = found.before.iter().map(|c| (c.enzyme.as_str(), c.distance)).collect();
        assert_eq!(before, vec![("EcoRI", 4), ("BamHI", 1599)]);
        // After the region: BamHI 8 bp away, then EcoRI the long way round.
        let after: Vec<(&str, usize)> = found.after.iter().map(|c| (c.enzyme.as_str(), c.distance)).collect();
        assert_eq!(after, vec![("BamHI", 8), ("EcoRI", 2496)]);
        // Each cut says where it is, the end it leaves and how often its enzyme cuts.
        let ecori = &found.before[0];
        assert_eq!((ecori.pos, ecori.end_text(), ecori.cuts, ecori.vector_cuts), (996, "5' AATT".to_string(), 1, None));
        let bamhi = &found.after[0];
        assert_eq!((bamhi.pos, bamhi.end_text(), bamhi.cuts), (1508, "5' GATC".to_string(), 2));
        assert!(found.inside.is_empty(), "{:?}", found.inside);
        // JSON names the same things.
        let value = found.to_json();
        assert_eq!(value["before"][0]["enzyme"], "EcoRI");
        assert_eq!(value["after"][0]["distance"], 8);
        assert_eq!(value["after"][0]["end"], "5' GATC");
    }

    #[test]
    fn an_enzyme_that_cuts_inside_the_region_is_set_apart() {
        let mut doc = donor();
        // HindIII (A^AGCTT) in the middle of geneX.
        doc.seq[1200..1206].copy_from_slice(b"AAGCTT");
        let found = nearest_cuts(&doc, region(), &common(), None, DEFAULT_PER_SIDE);
        assert_eq!(found.inside, vec!["HindIII".to_string()]);
        assert!(found.before.iter().chain(&found.after).all(|c| c.enzyme != "HindIII"));
        // Choosing it anyway says why it cannot be used.
        let error = excision(&doc, region(), enzyme("HindIII"), enzyme("BamHI"), None).unwrap_err().to_string();
        assert!(error.contains("HindIII cuts inside the region"), "{error}");
    }

    #[test]
    fn the_list_is_as_long_as_asked_and_counts_the_vectors_cuts() {
        let vector = vector();
        let found = nearest_cuts(&donor(), region(), &common(), Some(&vector), 1);
        assert_eq!((found.before.len(), found.after.len()), (1, 1));
        assert_eq!((found.before[0].enzyme.as_str(), found.before[0].vector_cuts), ("EcoRI", Some(1)));
        assert_eq!((found.after[0].enzyme.as_str(), found.after[0].vector_cuts), ("BamHI", Some(1)));
    }

    #[test]
    fn a_pick_gives_the_piece_that_holds_the_region() {
        let best = excision(&donor(), region(), enzyme("EcoRI"), enzyme("BamHI"), None).unwrap();
        assert_eq!((best.left.enzyme.as_str(), best.right.enzyme.as_str()), ("EcoRI", "BamHI"));
        assert_eq!((best.left.pos, best.right.pos), (996, 1508));
        assert_eq!((best.left.extra, best.right.extra), (4, 8));
        assert!(best.directional, "AATT and GATC ends differ");
        assert_eq!(best.ends_text(), "5' AATT / 5' GATC");
        // The fragment runs from the EcoRI cut to the end of the BamHI overhang: 516 bp.
        assert_eq!(best.fragment, Span::new(996, 1512));
        assert_eq!(best.length, 516);
        assert_eq!(best.other_cuts, 1, "the second BamHI site adds a band");
        assert_eq!(best.enzymes(), vec!["EcoRI".to_string(), "BamHI".to_string()]);
        // JSON counts positions from one.
        let value = best.to_json(3000);
        assert_eq!((value["start"].as_u64(), value["end"].as_u64(), value["length"].as_u64()), (Some(997), Some(1512), Some(516)));
        assert_eq!(value["left"]["cut"], 996);
    }

    #[test]
    fn the_same_enzyme_on_both_sides_goes_in_either_way() {
        let mut doc = donor();
        // A second EcoRI site 20 bp after the region, so EcoRI alone cuts the region out too.
        doc.seq[1519..1525].copy_from_slice(b"GAATTC");
        let single = excision(&doc, region(), enzyme("EcoRI"), enzyme("EcoRI"), None).unwrap();
        assert!(!single.directional);
        assert_eq!(single.enzymes(), vec!["EcoRI".to_string()]);
        assert_eq!(single.label(), "EcoRI");
        assert_eq!((single.left.extra, single.right.extra), (4, 20));
        assert_eq!(single.other_cuts, 0);
        // An enzyme that cuts the circle once only opens it.
        let error = excision(&donor(), region(), enzyme("EcoRI"), enzyme("EcoRI"), None).unwrap_err().to_string();
        assert!(error.contains("only once"), "{error}");
    }

    #[test]
    fn a_region_through_the_origin_is_cut_out_too() {
        // The same donor with its origin moved into geneX, so the region wraps.
        let mut doc = donor();
        crate::edit::set_origin(&mut doc, 1250).unwrap();
        let gene = doc.features.iter().find(|f| f.name == "geneX").unwrap().span(doc.len());
        assert!(gene.wraps(doc.len()), "{gene:?}");
        let found = nearest_cuts(&doc, gene, &common(), None, DEFAULT_PER_SIDE);
        assert_eq!((found.before[0].enzyme.as_str(), found.before[0].distance), ("EcoRI", 4));
        assert_eq!((found.after[0].enzyme.as_str(), found.after[0].distance), ("BamHI", 8));
        let best = excision(&doc, gene, enzyme("EcoRI"), enzyme("BamHI"), None).unwrap();
        assert_eq!((best.left.extra, best.right.extra), (4, 8));
        assert_eq!(best.length, 516);
        assert!(best.fragment.contains_span(&gene, doc.len()));
    }

    #[test]
    fn a_linear_donor_has_no_cuts_past_its_ends() {
        let mut doc = donor();
        doc.topology = Topology::Linear;
        let found = nearest_cuts(&doc, region(), &common(), None, DEFAULT_PER_SIDE);
        // Each cut lies on one side only, so each enzyme appears where it really is.
        let names = |side: &[NearCut]| side.iter().map(|c| c.enzyme.clone()).collect::<Vec<_>>();
        assert_eq!((names(&found.before), names(&found.after)), (vec!["EcoRI".to_string()], vec!["BamHI".to_string()]));
        assert!(excision(&doc, region(), enzyme("EcoRI"), enzyme("BamHI"), None).is_ok());
        // A region at the very start has no cut before it.
        let found = nearest_cuts(&doc, Span::new(0, 500), &common(), None, DEFAULT_PER_SIDE);
        assert!(found.before.is_empty());
        assert!(excision(&doc, Span::new(0, 500), enzyme("EcoRI"), enzyme("BamHI"), None).is_err());
    }

    /// A 4,000 bp circular vector with an EcoRI site and a BamHI site 21 bp apart in a short
    /// multiple cloning site at 500, and a KpnI site elsewhere.
    fn vector() -> Document {
        let mut seq = clean(4000, 2);
        seq[500..506].copy_from_slice(b"GAATTC");
        seq[521..527].copy_from_slice(b"GGATCC");
        seq[3000..3006].copy_from_slice(b"GGTACC");
        let mut doc = Document::new_dna("vector", seq, Topology::Circular);
        doc.features.push(Feature::new("MCS", "misc_feature", Span::new(495, 535), Strand::Forward));
        assert_eq!(cutters(&doc), ["BamHI", "EcoRI", "KpnI"]);
        doc
    }

    #[test]
    fn a_pick_is_checked_on_the_vector() {
        let vector = vector();
        let best = excision(&donor(), region(), enzyme("EcoRI"), enzyme("BamHI"), Some(&vector)).unwrap();
        let fit = best.vector.as_ref().unwrap();
        assert_eq!(fit.cuts, vec![("EcoRI".to_string(), 501), ("BamHI".to_string(), 522)]);
        assert_eq!(fit.dropped, 21);
        assert_eq!(fit.problem, None);
        assert_eq!(fit.describe(), "cuts once each, at 501 and 522, dropping 21 bp");
        // A vector with a second BamHI site refuses the pair, and says so.
        let mut twice = vector.clone();
        twice.seq[2000..2006].copy_from_slice(b"GGATCC");
        let pair = excision(&donor(), region(), enzyme("EcoRI"), enzyme("BamHI"), Some(&twice)).unwrap();
        assert_eq!(pair.vector.as_ref().unwrap().problem.as_deref(), Some("BamHI cuts the vector 2 times"));
        // A vector without an EcoRI site says that too.
        let mut none = vector.clone();
        none.seq[500..506].copy_from_slice(b"ACACAC");
        let pair = excision(&donor(), region(), enzyme("EcoRI"), enzyme("BamHI"), Some(&none)).unwrap();
        assert_eq!(pair.vector.as_ref().unwrap().problem.as_deref(), Some("EcoRI does not cut the vector"));
    }

    /// A BsaI site reads GGTCTC and cuts one base after it on the top strand, five on the
    /// bottom, so each site leaves its own four-base overhang. Here the donor's region sits
    /// between two BsaI sites whose overhangs are CATG and TGCA, and the vector carries two
    /// BsaI sites with the overhangs given.
    fn bsai_pair(vector_left: &[u8; 4], vector_right: &[u8; 4]) -> (Document, Document) {
        let mut donor = clean(2000, 3);
        // Forward site at 989: cut at 996 and overhang 996..1000, just before the region.
        donor[989..995].copy_from_slice(b"GGTCTC");
        donor[996..1000].copy_from_slice(b"CATG");
        // Reverse site at 1509 (GAGACC): cut at 1504 and overhang 1504..1508, after the region.
        donor[1504..1508].copy_from_slice(b"TGCA");
        donor[1509..1515].copy_from_slice(b"GAGACC");
        let donor = Document::new_dna("donor", donor, Topology::Circular);
        let mut vector = clean(3000, 4);
        vector[489..495].copy_from_slice(b"GGTCTC");
        vector[496..500].copy_from_slice(vector_left);
        vector[554..558].copy_from_slice(vector_right);
        vector[559..565].copy_from_slice(b"GAGACC");
        let vector = Document::new_dna("vector", vector, Topology::Circular);
        assert_eq!((cutters(&donor), cutters(&vector)), (vec!["BsaI".to_string()], vec!["BsaI".to_string()]));
        (donor, vector)
    }

    #[test]
    fn each_site_of_a_type_iis_enzyme_leaves_its_own_end() {
        let bsai = EnzymeDb::builtin().resolve(&["BsaI"]).unwrap();
        let (donor, vector) = bsai_pair(b"CATG", b"TGCA");
        // The list shows the overhang each site leaves.
        let found = nearest_cuts(&donor, region(), &bsai, Some(&vector), DEFAULT_PER_SIDE);
        assert_eq!((found.before[0].end_text(), found.after[0].end_text()), ("5' CATG".to_string(), "5' TGCA".to_string()));
        let option = excision(&donor, region(), bsai[0], bsai[0], Some(&vector)).unwrap();
        assert_eq!(option.ends_text(), "5' CATG / 5' TGCA");
        assert!(option.directional);
        let fit = option.vector.as_ref().unwrap();
        assert_eq!(fit.problem, None);
        assert_eq!(fit.describe(), "cuts twice, at 496 and 554, dropping 58 bp");
        // A vector whose sites leave other overhangs cannot take the piece.
        let (donor, vector) = bsai_pair(b"GGGG", b"TGCA");
        let option = excision(&donor, region(), bsai[0], bsai[0], Some(&vector)).unwrap();
        assert_eq!(option.vector.unwrap().problem.as_deref(), Some("the ends it leaves on the vector do not pair with the insert's"));
    }

    #[test]
    fn an_end_joins_a_piece_with_a_matching_end_either_way_round() {
        // The EcoRI-BamHI piece of the donor has ends AATT and GATC.
        let piece = excision(&donor(), region(), enzyme("EcoRI"), enzyme("BamHI"), None).unwrap();
        let cut = digest(&donor(), &[enzyme("EcoRI"), enzyme("BamHI")]).unwrap();
        let fragment = cut.fragments.iter().find(|f| f.span == piece.fragment).unwrap();
        assert!(can_join(Overhang::FivePrime(4), "AATT", &fragment.doc));
        assert!(can_join(Overhang::FivePrime(4), "gatc", &fragment.doc));
        assert!(!can_join(Overhang::FivePrime(4), "CTAG", &fragment.doc), "an XbaI end does not pair");
        assert!(!can_join(Overhang::Blunt, "", &fragment.doc));
        assert!(!can_join(Overhang::FivePrime(4), "AATT", &donor()), "a circle has no ends");
    }

    #[test]
    fn nothing_comes_of_an_empty_or_whole_region() {
        let doc = donor();
        assert_eq!(nearest_cuts(&doc, Span::new(100, 100), &common(), None, DEFAULT_PER_SIDE), NearestCuts::default());
        assert_eq!(nearest_cuts(&doc, Span::new(0, 3000), &common(), None, DEFAULT_PER_SIDE), NearestCuts::default());
        assert!(excision(&doc, Span::new(0, 3000), enzyme("EcoRI"), enzyme("BamHI"), None).is_err());
    }
}
