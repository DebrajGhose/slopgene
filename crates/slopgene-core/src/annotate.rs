//! Automatic annotation from a library of known features.
//!
//! The library is built from annotated files the user already trusts (their own plasmid maps,
//! public GenBank records). Each library feature's sequence is looked for in the target on both
//! strands, through the origin of circular targets, allowing a set fraction of mismatched
//! bases. This is the job SnapGene's "Detect Common Features" and pLannotate do, without
//! relying on anyone else's proprietary feature database.

use crate::alphabet;
use crate::model::*;
use crate::span::Span;
use std::collections::{HashMap, HashSet};

const SEED: usize = 16;

/// One known feature.
#[derive(Clone, Debug)]
pub struct LibraryEntry {
    pub name: String,
    pub kind: String,
    /// The feature's bases, 5'->3' along its own strand.
    pub seq: Vec<u8>,
    /// Whether the feature has a direction (drawn as an arrow).
    pub directional: bool,
    pub color: Option<Color>,
    /// Where the entry came from (document name).
    pub source: String,
}

/// A set of known features to annotate with.
#[derive(Clone, Debug, Default)]
pub struct FeatureLibrary {
    pub entries: Vec<LibraryEntry>,
    seen: HashSet<Vec<u8>>,
}

impl FeatureLibrary {
    pub fn new() -> FeatureLibrary {
        FeatureLibrary::default()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Add every feature of `doc` at least `min_len` bases long. Features with the same bases
    /// as one already in the library are skipped (the first name wins).
    pub fn add_document(&mut self, doc: &Document, min_len: usize) -> usize {
        let mut added = 0;
        for feature in &doc.features {
            if feature.kind.eq_ignore_ascii_case("source") || feature.len() < min_len.max(SEED) {
                continue;
            }
            let seq = alphabet::upper(&feature.extract(&doc.seq));
            // The same element seen in the other orientation is the same entry.
            let key = seq.clone().min(alphabet::revcomp(&seq));
            if !self.seen.insert(key) {
                continue;
            }
            self.entries.push(LibraryEntry {
                name: feature.name.clone(),
                kind: feature.kind.clone(),
                seq,
                directional: feature.strand != Strand::None,
                color: feature.color(),
                source: doc.name.clone(),
            });
            added += 1;
        }
        added
    }

    /// The library as one linear document: entries laid end to end, each annotated (the same
    /// layout SnapGene uses for its feature-library files), so it can be saved and reused.
    pub fn to_document(&self, name: &str) -> Document {
        let mut doc = Document::new_dna(name, b"", Topology::Linear);
        for entry in &self.entries {
            let start = doc.seq.len();
            doc.seq.extend_from_slice(&entry.seq);
            let strand = if entry.directional { Strand::Forward } else { Strand::None };
            let mut f = Feature::new(entry.name.clone(), entry.kind.clone(), Span::new(start, doc.seq.len()), strand);
            f.set_color(entry.color);
            if entry.kind.eq_ignore_ascii_case("CDS") {
                f.segments[0].translated = Some(true);
            }
            f.set_qualifier("note", QualValue::Text(format!("from {}", entry.source)));
            doc.features.push(f);
        }
        doc
    }
}

/// One feature found in a target.
#[derive(Clone, Debug)]
pub struct Annotation {
    pub feature: Feature,
    /// Matching bases as a fraction of the feature length (0..1).
    pub identity: f64,
    pub source: String,
}

fn pack(kmer: &[u8]) -> Option<u64> {
    let mut v = 0u64;
    for &b in kmer {
        v = (v << 2)
            | match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' | b'U' => 3,
                _ => return None,
            };
    }
    Some(v)
}

/// Find library features in `doc`. Matches must cover the whole library feature with at least
/// `min_identity` (0..1) of its bases matching, without gaps. Where several entries match the
/// same stretch, the most identical (then the longest) is kept. Stretches already annotated
/// with a feature of the same name are skipped.
pub fn annotate(doc: &Document, library: &FeatureLibrary, min_identity: f64) -> Vec<Annotation> {
    let seq = alphabet::upper(&doc.seq);
    let len = seq.len();
    let circular = doc.is_circular();
    if len < SEED {
        return Vec::new();
    }
    // Seed index of the target.
    let mut index: HashMap<u64, Vec<u32>> = HashMap::new();
    let starts = if circular { len } else { len - SEED + 1 };
    let mut window = Vec::with_capacity(SEED);
    for start in 0..starts {
        window.clear();
        window.extend((0..SEED).map(|ii| seq[(start + ii) % len]));
        if let Some(code) = pack(&window) {
            index.entry(code).or_default().push(start as u32);
        }
    }

    let mut found: Vec<Annotation> = Vec::new();
    for entry in &library.entries {
        let n = entry.seq.len();
        if n < SEED || n > len {
            continue;
        }
        let allowed = ((1.0 - min_identity) * n as f64).floor() as usize;
        // ii = 0 looks for the entry on the top strand, ii = 1 on the bottom strand.
        for ii in 0..2 {
            let probe = if ii == 0 { entry.seq.clone() } else { alphabet::revcomp(&entry.seq) };
            // Seeds at the start, middle and end, so a few mismatches cannot hide a match.
            let mut candidates: HashSet<usize> = HashSet::new();
            for offset in [0, n / 2 - SEED / 2, n - SEED] {
                let Some(code) = pack(&probe[offset..offset + SEED]) else { continue };
                if let Some(hits) = index.get(&code) {
                    for &h in hits.iter().take(64) {
                        let start = (h as i64 - offset as i64).rem_euclid(if circular { len as i64 } else { i64::MAX }) as usize;
                        if circular || start + n <= len {
                            candidates.insert(start);
                        }
                    }
                }
            }
            for start in candidates {
                let mut mismatches = 0;
                let mut ok = true;
                for (k, &b) in probe.iter().enumerate() {
                    if !alphabet::iupac_matches(b, seq[(start + k) % len]) {
                        mismatches += 1;
                        if mismatches > allowed {
                            ok = false;
                            break;
                        }
                    }
                }
                if !ok {
                    continue;
                }
                let strand = match (entry.directional, ii) {
                    (false, _) => Strand::None,
                    (true, 0) => Strand::Forward,
                    (true, _) => Strand::Reverse,
                };
                let span = Span::new(start, start + n).normalize(len);
                let mut feature = Feature::new(entry.name.clone(), entry.kind.clone(), span, strand);
                feature.set_color(entry.color);
                if entry.kind.eq_ignore_ascii_case("CDS") {
                    feature.segments[0].translated = Some(true);
                }
                let identity = 1.0 - mismatches as f64 / n as f64;
                feature.set_qualifier(
                    "note",
                    QualValue::Text(format!("detected by SlopGene: {:.1}% identical to {} in {}", identity * 100.0, entry.name, entry.source)),
                );
                found.push(Annotation { feature, identity, source: entry.source.clone() });
            }
        }
    }

    // One annotation per stretch: the best identity wins, then the longer feature.
    found.sort_by(|a, b| b.identity.partial_cmp(&a.identity).unwrap_or(std::cmp::Ordering::Equal).then(b.feature.len().cmp(&a.feature.len())));
    let mut kept: Vec<Annotation> = Vec::new();
    for a in found {
        let span = a.feature.span(len);
        // One label per element: skip a hit when a feature of the same name (already in the
        // file or found just now) covers most of it.
        let named_twice = doc
            .features
            .iter()
            .chain(kept.iter().map(|k| &k.feature))
            .any(|f| f.name.eq_ignore_ascii_case(&a.feature.name) && mostly_overlaps(f.span(len), span, len));
        if named_twice {
            continue;
        }
        // A synonym (same type, nearly or exactly the same stretch as a better hit) becomes an
        // /also_called name on that hit.
        let synonym_of = kept.iter_mut().find(|k| k.feature.kind.eq_ignore_ascii_case(&a.feature.kind) && nearly_same(k.feature.span(len), span, len));
        if let Some(k) = synonym_of {
            let entry = QualValue::Text(a.feature.name.clone());
            match k.feature.qualifiers.iter_mut().find(|q| q.key == "also_called") {
                Some(q) if q.values.contains(&entry) => {}
                Some(q) => q.values.push(entry),
                None => k.feature.qualifiers.push(Qualifier { key: "also_called".into(), values: vec![entry] }),
            }
            continue;
        }
        // A hit of another type on exactly the stretch of a better hit, or a hit that repeats a
        // feature the file already has, adds nothing.
        let same_stretch = kept.iter().any(|k| k.feature.span(len) == span);
        let already_annotated = doc.features.iter().any(|f| f.kind.eq_ignore_ascii_case(&a.feature.kind) && nearly_same(f.span(len), span, len));
        if same_stretch || already_annotated {
            continue;
        }
        kept.push(a);
    }
    kept.sort_by_key(|a| a.feature.span(len).start);
    kept
}

/// Whether two spans cover nearly the same bases: each shares at least 80% of itself.
fn nearly_same(a: Span, b: Span, len: usize) -> bool {
    let mut shared = 0;
    for (a0, a1) in a.pieces(len) {
        for (b0, b1) in b.pieces(len) {
            shared += a1.min(b1).saturating_sub(a0.max(b0));
        }
    }
    shared * 5 >= a.len() * 4 && shared * 5 >= b.len() * 4
}

/// Whether two spans share more than half of the shorter one.
fn mostly_overlaps(a: Span, b: Span, len: usize) -> bool {
    let mut shared = 0;
    for (a0, a1) in a.pieces(len) {
        for (b0, b1) in b.pieces(len) {
            shared += a1.min(b1).saturating_sub(a0.max(b0));
        }
    }
    shared * 2 > a.len().min(b.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_dna(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn finds_features_on_both_strands_and_through_the_origin() {
        let gene = random_dna(300, 1);
        let promoter = random_dna(80, 2);
        let mut source = Document::new_dna("lib", [promoter.clone(), random_dna(50, 3), gene.clone()].concat(), Topology::Linear);
        source.features.push(Feature::new("pX", "promoter", Span::new(0, 80), Strand::Forward));
        source.features.push(Feature::new("geneY", "CDS", Span::new(130, 430), Strand::Forward));
        let mut library = FeatureLibrary::new();
        assert_eq!(library.add_document(&source, 15), 2);

        // Target: geneY reversed, and the promoter split across the origin with one mismatch.
        let mut target = random_dna(1000, 4);
        target[200..500].copy_from_slice(&alphabet::revcomp(&gene));
        let mut p = promoter.clone();
        p[40] = if p[40] == b'A' { b'C' } else { b'A' };
        target[960..1000].copy_from_slice(&p[..40]);
        target[..40].copy_from_slice(&p[40..]);
        let doc = Document::new_dna("target", target, Topology::Circular);

        let hits = annotate(&doc, &library, 0.95);
        assert_eq!(hits.len(), 2, "{hits:?}");
        let gene_hit = hits.iter().find(|h| h.feature.name == "geneY").unwrap();
        assert_eq!(gene_hit.feature.strand, Strand::Reverse);
        assert_eq!(gene_hit.feature.span(1000), Span::new(200, 500));
        let promoter_hit = hits.iter().find(|h| h.feature.name == "pX").unwrap();
        assert_eq!(promoter_hit.feature.span(1000), Span::new(960, 1040));
        assert!(promoter_hit.identity < 1.0 && promoter_hit.identity > 0.98);
        // Strict identity misses the mutated promoter.
        assert_eq!(annotate(&doc, &library, 1.0).len(), 1);
    }

    #[test]
    fn synonyms_on_the_same_stretch_become_also_called() {
        // Two library maps name the same 200 bp element differently; one copy has a point
        // difference. The exact copy wins and the other name is kept as /also_called.
        let element = random_dna(200, 7);
        let mut variant = element.clone();
        variant[100] = if variant[100] == b'G' { b'T' } else { b'G' };
        let mut first = Document::new_dna("map1", element.clone(), Topology::Linear);
        first.features.push(Feature::new("ori", "rep_origin", Span::new(0, 200), Strand::Forward));
        let mut second = Document::new_dna("map2", variant, Topology::Linear);
        second.features.push(Feature::new("pUC ori", "rep_origin", Span::new(0, 200), Strand::Forward));
        let mut library = FeatureLibrary::new();
        library.add_document(&first, 30);
        library.add_document(&second, 30);
        assert_eq!(library.len(), 2);

        let mut target = random_dna(900, 8);
        target[300..500].copy_from_slice(&element);
        let doc = Document::new_dna("target", target, Topology::Circular);
        let hits = annotate(&doc, &library, 0.95);
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert_eq!(hits[0].feature.name, "ori");
        let also = hits[0].feature.qualifiers.iter().find(|q| q.key == "also_called").expect("also_called");
        assert_eq!(also.values, vec![QualValue::Text("pUC ori".into())]);
    }
}
