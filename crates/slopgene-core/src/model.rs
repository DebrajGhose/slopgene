//! The sequence document: one DNA, RNA or protein molecule with its annotations.
//!
//! All coordinates in the model are 0-based and half-open ([`Span`]). The model is format
//! neutral; the readers in [`crate::io`] fill it in, and anything a format stores that the model
//! does not cover is kept in [`RawData`] so files survive a read-write cycle unchanged.

use crate::span::Span;
use crate::xml::Element;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Molecule {
    #[default]
    Dna,
    Rna,
    Protein,
}

impl Molecule {
    pub fn name(self) -> &'static str {
        match self {
            Molecule::Dna => "dna",
            Molecule::Rna => "rna",
            Molecule::Protein => "protein",
        }
    }

    pub fn is_nucleic(self) -> bool {
        !matches!(self, Molecule::Protein)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Topology {
    #[default]
    Linear,
    Circular,
}

impl Topology {
    pub fn name(self) -> &'static str {
        match self {
            Topology::Linear => "linear",
            Topology::Circular => "circular",
        }
    }

    pub fn parse(text: &str) -> Option<Topology> {
        match text.trim().to_ascii_lowercase().as_str() {
            "linear" | "lin" | "l" => Some(Topology::Linear),
            "circular" | "circ" | "c" => Some(Topology::Circular),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Strandedness {
    #[default]
    Double,
    Single,
}

/// Which way a feature or primer points. SnapGene calls this "directionality".
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Strand {
    /// No direction (drawn as a plain box).
    #[default]
    None,
    /// Reads 5'->3' along the top strand.
    Forward,
    /// Reads 5'->3' along the bottom strand.
    Reverse,
    /// Points both ways.
    Both,
}

impl Strand {
    /// The strand seen from the other side of the molecule (after reverse-complementing).
    pub fn flipped(self) -> Strand {
        match self {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
            other => other,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Strand::None => "none",
            Strand::Forward => "forward",
            Strand::Reverse => "reverse",
            Strand::Both => "both",
        }
    }

    pub fn parse(text: &str) -> Option<Strand> {
        match text.trim().to_ascii_lowercase().as_str() {
            "none" | "0" | "." | "" | "nondirectional" => Some(Strand::None),
            "forward" | "fwd" | "f" | "+" | "1" | "+1" | "top" => Some(Strand::Forward),
            "reverse" | "rev" | "r" | "-" | "-1" | "bottom" | "complement" => Some(Strand::Reverse),
            "both" | "bidirectional" | "+-" | "2" => Some(Strand::Both),
            _ => None,
        }
    }

    pub fn is_reverse(self) -> bool {
        matches!(self, Strand::Reverse)
    }
}

/// An sRGB colour.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl Color {
    pub const fn rgb(r: u8, g: u8, b: u8) -> Color {
        Color { r, g, b }
    }

    /// Parse `#rrggbb`, `rrggbb` or `#rgb`.
    pub fn parse_hex(text: &str) -> Option<Color> {
        let hex = text.trim().trim_start_matches('#');
        let channel = |s: &str| u8::from_str_radix(s, 16).ok();
        match hex.len() {
            6 => Some(Color::rgb(channel(&hex[0..2])?, channel(&hex[2..4])?, channel(&hex[4..6])?)),
            3 => {
                let expand = |s: &str| channel(s).map(|v| v * 17);
                Some(Color::rgb(expand(&hex[0..1])?, expand(&hex[1..2])?, expand(&hex[2..3])?))
            }
            _ => None,
        }
    }

    pub fn to_hex(self) -> String {
        format!("#{:02x}{:02x}{:02x}", self.r, self.g, self.b)
    }
}

/// One value of a feature qualifier (GenBank `/key=value`; SnapGene `<Q><V .../></Q>`).
#[derive(Clone, Debug, PartialEq)]
pub enum QualValue {
    Text(String),
    Int(i64),
    /// SnapGene "predefined" value: a member of a fixed vocabulary, stored as text.
    Predef(String),
    /// A qualifier with no value, such as GenBank `/pseudo`.
    Empty,
    /// A SnapGene value element with several attributes (for example both `predef` and
    /// `text`), kept exactly as written.
    Raw(Vec<(String, String)>),
}

impl QualValue {
    pub fn as_text(&self) -> String {
        match self {
            QualValue::Text(t) | QualValue::Predef(t) => t.clone(),
            QualValue::Int(n) => n.to_string(),
            QualValue::Empty => String::new(),
            QualValue::Raw(attrs) => ["text", "predef", "int"]
                .iter()
                .find_map(|key| attrs.iter().find(|(k, _)| k == key).map(|(_, v)| v.clone()))
                .or_else(|| attrs.first().map(|(_, v)| v.clone()))
                .unwrap_or_default(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Qualifier {
    pub key: String,
    pub values: Vec<QualValue>,
}

/// What a SnapGene element looked like when read, so it can be written back the same way.
///
/// `attrs` holds every attribute in document order, including the ones the model has its own
/// fields for (name, type, ...). On writing, modelled attributes take the model's current
/// values in their original positions; the others are copied as they were. `children` holds
/// child elements the model does not cover.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct XmlExtra {
    pub attrs: Vec<(String, String)>,
    pub children: Vec<Element>,
}

impl XmlExtra {
    pub fn attr(&self, name: &str) -> Option<&str> {
        self.attrs.iter().find(|(k, _)| k == name).map(|(_, v)| v.as_str())
    }
}

/// One contiguous piece of a feature.
#[derive(Clone, Debug, PartialEq)]
pub struct Segment {
    pub span: Span,
    pub color: Option<Color>,
    pub name: Option<String>,
    /// SnapGene "gap" segment: part of the feature's extent that is not part of the feature
    /// itself (for example an intron drawn as a thin line).
    pub gap: bool,
    /// Whether this segment is translated when the feature is shown with its protein.
    pub translated: Option<bool>,
    /// Every attribute as read, in order (see [`XmlExtra`]).
    pub extra: Vec<(String, String)>,
}

impl Segment {
    pub fn new(span: Span) -> Segment {
        Segment { span, color: None, name: None, gap: false, translated: None, extra: Vec::new() }
    }
}

static NEXT_UID: AtomicU64 = AtomicU64::new(1);

/// A process-wide unique id, used to track features and primers across edits (selection,
/// undo). Ids are not saved to files.
pub fn new_uid() -> u64 {
    NEXT_UID.fetch_add(1, Ordering::Relaxed)
}

/// An annotated region of the sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct Feature {
    pub uid: u64,
    pub name: String,
    /// GenBank feature key, such as `CDS`, `promoter` or `misc_feature`.
    pub kind: String,
    pub strand: Strand,
    /// Pieces in top-strand order (left to right), following the feature through the origin
    /// when it wraps. At least one.
    pub segments: Vec<Segment>,
    pub qualifiers: Vec<Qualifier>,
    pub extra: XmlExtra,
}

impl Feature {
    pub fn new(name: impl Into<String>, kind: impl Into<String>, span: Span, strand: Strand) -> Feature {
        Feature {
            uid: new_uid(),
            name: name.into(),
            kind: kind.into(),
            strand,
            segments: vec![Segment::new(span)],
            qualifiers: Vec::new(),
            extra: XmlExtra::default(),
        }
    }

    /// The whole extent, from the first segment's start to the last segment's end.
    pub fn span(&self, seq_len: usize) -> Span {
        let Some(first) = self.segments.first() else {
            return Span::default();
        };
        let last = self.segments.last().unwrap();
        let start = first.span.start;
        let mut end = last.span.end;
        // A later segment that starts before the first one has passed the origin.
        if self.segments.len() > 1 && last.span.start < start {
            end += seq_len;
        }
        Span::new(start, end.max(start))
    }

    /// Bases covered by the non-gap segments.
    pub fn len(&self) -> usize {
        self.segments.iter().filter(|s| !s.gap).map(|s| s.span.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Colour of the first segment that has one.
    pub fn color(&self) -> Option<Color> {
        self.segments.iter().find_map(|s| s.color)
    }

    /// Give every segment the same colour.
    pub fn set_color(&mut self, color: Option<Color>) {
        for segment in &mut self.segments {
            segment.color = color;
        }
    }

    pub fn qualifier(&self, key: &str) -> Option<&Qualifier> {
        self.qualifiers.iter().find(|q| q.key == key)
    }

    /// Text of the first value of qualifier `key`.
    pub fn qualifier_text(&self, key: &str) -> Option<String> {
        self.qualifier(key).and_then(|q| q.values.first()).map(|v| v.as_text())
    }

    /// Set qualifier `key` to a single value, replacing what was there.
    pub fn set_qualifier(&mut self, key: &str, value: QualValue) {
        match self.qualifiers.iter_mut().find(|q| q.key == key) {
            Some(q) => q.values = vec![value],
            None => self.qualifiers.push(Qualifier { key: key.to_string(), values: vec![value] }),
        }
    }

    pub fn remove_qualifier(&mut self, key: &str) {
        self.qualifiers.retain(|q| q.key != key);
    }

    /// Whether the feature codes for protein (a CDS, or any segment marked as translated).
    pub fn is_coding(&self) -> bool {
        self.kind.eq_ignore_ascii_case("CDS") || self.segments.iter().any(|s| s.translated == Some(true))
    }

    /// The feature's own sequence, read 5'->3' along its strand: non-gap segments joined in
    /// order, reverse-complemented for reverse-strand features.
    pub fn extract(&self, seq: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.len());
        for segment in self.segments.iter().filter(|s| !s.gap) {
            out.extend(segment.span.extract(seq));
        }
        if self.strand.is_reverse() {
            crate::alphabet::revcomp_in_place(&mut out);
        }
        out
    }

    /// Put segments into top-strand order starting from the leftmost one, and make sure every
    /// span starts inside the sequence.
    pub fn normalize(&mut self, seq_len: usize, circular: bool) {
        for segment in &mut self.segments {
            segment.span = segment.span.normalize(seq_len);
        }
        if self.segments.len() < 2 {
            return;
        }
        if !circular {
            self.segments.sort_by_key(|s| s.span.start);
            return;
        }
        // On a circular molecule keep the order as written, but rotate the list so it begins at
        // the segment that follows the largest gap. That segment is where the feature starts.
        let count = self.segments.len();
        let mut best_start = 0;
        let mut best_gap = 0usize;
        for ii in 0..count {
            let prev_end = self.segments[(ii + count - 1) % count].span.end % seq_len.max(1);
            let this_start = self.segments[ii].span.start;
            let gap = (this_start + seq_len - prev_end) % seq_len.max(1);
            if ii == 0 || gap > best_gap {
                best_gap = gap;
                best_start = ii;
            }
        }
        self.segments.rotate_left(best_start);
    }
}

/// Where a primer anneals.
#[derive(Clone, Debug, PartialEq)]
pub struct BindingSite {
    /// The annealed bases, in top-strand coordinates.
    pub span: Span,
    /// `Forward`: the primer reads like the top strand (it anneals to the bottom strand).
    /// `Reverse`: the primer reads like the bottom strand.
    pub strand: Strand,
    pub tm: Option<f64>,
    pub extra: XmlExtra,
}

/// An oligonucleotide, usually a PCR or sequencing primer.
#[derive(Clone, Debug, PartialEq)]
pub struct Primer {
    pub uid: u64,
    pub name: String,
    /// 5'->3'.
    pub seq: String,
    pub description: String,
    /// Binding sites as stored in the file. Recompute with [`crate::primers::find_binding_sites`]
    /// after editing the sequence.
    pub sites: Vec<BindingSite>,
    pub extra: XmlExtra,
}

impl Primer {
    pub fn new(name: impl Into<String>, seq: impl Into<String>) -> Primer {
        Primer { uid: new_uid(), name: name.into(), seq: seq.into(), description: String::new(), sites: Vec::new(), extra: XmlExtra::default() }
    }
}

/// What one end of a linear double-stranded molecule looks like.
///
/// The document's sequence always spans the full extent of both strands, so overhanging bases
/// are part of `Document::seq`. `FivePrime(n)` on the left end means the top strand sticks out by
/// `n` bases; `ThreePrime(n)` on the left end means the bottom strand sticks out. On the right
/// end `ThreePrime(n)` means the top strand sticks out and `FivePrime(n)` the bottom strand.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Overhang {
    #[default]
    Blunt,
    FivePrime(usize),
    ThreePrime(usize),
}

impl Overhang {
    /// Length of the single-stranded overhang in bases (0 for a blunt end). An end is never
    /// "empty", so there is no `is_empty` to go with it.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(self) -> usize {
        match self {
            Overhang::Blunt => 0,
            Overhang::FivePrime(n) | Overhang::ThreePrime(n) => n,
        }
    }

    pub fn is_blunt(self) -> bool {
        self.len() == 0
    }

    pub fn describe(self) -> String {
        match self {
            Overhang::Blunt => "blunt".to_string(),
            Overhang::FivePrime(n) => format!("5' overhang ({n} nt)"),
            Overhang::ThreePrime(n) => format!("3' overhang ({n} nt)"),
        }
    }
}

/// Both ends of a linear double-stranded molecule, plus their 5' phosphorylation state.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Ends {
    pub left: Overhang,
    pub right: Overhang,
    /// SnapGene's end modification name, for example `FivePrimePhosphorylated` or `Unmodified`.
    pub left_modification: Option<String>,
    pub right_modification: Option<String>,
}

/// Methylation of the DNA by host enzymes, which can block some restriction enzymes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Methylation {
    pub dam: bool,
    pub dcm: bool,
    pub ecoki: bool,
}

/// A literature reference attached to a sequence.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Reference {
    pub authors: String,
    pub title: String,
    pub journal: String,
    pub pubmed: Option<String>,
}

/// Descriptive metadata. Dates are ISO `YYYY-MM-DD`.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Meta {
    pub description: String,
    pub accession: Option<String>,
    pub version: Option<String>,
    pub organism: Option<String>,
    pub keywords: Option<String>,
    pub comments: Vec<String>,
    pub references: Vec<Reference>,
    pub created: Option<String>,
    pub modified: Option<String>,
    pub created_by: Option<String>,
    /// `Some(true)` for a synthetic (designed) sequence, `Some(false)` for a natural one.
    pub synthetic: Option<bool>,
    /// GenBank division code (`SYN`, `UNA`, `BCT`, ...).
    pub division: Option<String>,
}

/// One packet of a SnapGene file, as read.
#[derive(Clone, Debug, PartialEq)]
pub struct RawPacket {
    pub kind: u8,
    pub data: Vec<u8>,
}

/// What a SnapGene file contained beyond the model, so it can be written back faithfully.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SnapGeneRaw {
    pub seq_type: u16,
    pub export_version: u16,
    pub import_version: u16,
    /// Every packet in file order, including the ones the model was read from.
    pub packets: Vec<RawPacket>,
    /// Fingerprint of the sequence, topology and strandedness at load time. A writer compares it
    /// with the current document to decide whether sequence-dependent packets went stale.
    pub fingerprint: u64,
    /// `<Features>` root without its `<Feature>` children (keeps root attributes).
    pub features_root: Option<Element>,
    /// `<Primers>` root without its `<Primer>` children (keeps `HybridizationParams`).
    pub primers_root: Option<Element>,
    pub notes: Option<Element>,
    pub properties: Option<Element>,
    /// Description as the original HTML, and the plain text it was read as.
    pub description_html: Option<(String, String)>,
}

/// Format-specific leftovers for lossless round trips.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RawData {
    pub snapgene: Option<SnapGeneRaw>,
    /// Qualifiers of the GenBank `source` feature, which is folded into the metadata on read.
    pub genbank_source: Option<Vec<Qualifier>>,
}

/// One molecule with its annotations.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Document {
    pub name: String,
    pub molecule: Molecule,
    pub topology: Topology,
    pub strandedness: Strandedness,
    /// The sequence as ASCII bytes (IUPAC codes for nucleotides, one-letter codes for protein).
    pub seq: Vec<u8>,
    pub features: Vec<Feature>,
    pub primers: Vec<Primer>,
    pub ends: Ends,
    pub methylation: Methylation,
    pub meta: Meta,
    pub raw: RawData,
}

impl Document {
    /// A new DNA document with no annotations.
    pub fn new_dna(name: impl Into<String>, seq: impl AsRef<[u8]>, topology: Topology) -> Document {
        Document { name: name.into(), molecule: Molecule::Dna, topology, strandedness: Strandedness::Double, seq: seq.as_ref().to_vec(), ..Document::default() }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn is_circular(&self) -> bool {
        self.topology == Topology::Circular
    }

    /// The sequence as text. Sequences are ASCII, so this never fails.
    pub fn seq_str(&self) -> String {
        String::from_utf8_lossy(&self.seq).into_owned()
    }

    /// Bases of `span` (following the origin on circular sequences).
    pub fn subseq(&self, span: Span) -> Vec<u8> {
        span.extract(&self.seq)
    }

    /// Index of the feature with this uid.
    pub fn feature_index(&self, uid: u64) -> Option<usize> {
        self.features.iter().position(|f| f.uid == uid)
    }

    /// Features whose name matches `name` (case-insensitive).
    pub fn features_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Feature> + 'a {
        self.features.iter().filter(move |f| f.name.eq_ignore_ascii_case(name))
    }

    /// Fingerprint of everything that sequence-dependent data (enzyme caches, history) depends
    /// on: the bases, topology and strandedness. FNV-1a, 64 bit.
    pub fn sequence_fingerprint(&self) -> u64 {
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |byte: u8| {
            hash ^= byte as u64;
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        };
        for &b in &self.seq {
            feed(b);
        }
        feed(self.topology as u8);
        feed(self.strandedness as u8);
        hash
    }

    /// Whether the sequence differs from what was read from a SnapGene file.
    pub fn sequence_changed_since_load(&self) -> bool {
        match &self.raw.snapgene {
            Some(raw) => raw.fingerprint != self.sequence_fingerprint(),
            None => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colors_parse_and_format() {
        assert_eq!(Color::parse_hex("#31849b"), Some(Color::rgb(0x31, 0x84, 0x9b)));
        assert_eq!(Color::parse_hex("f0a"), Some(Color::rgb(0xff, 0x00, 0xaa)));
        assert_eq!(Color::rgb(1, 2, 255).to_hex(), "#0102ff");
        assert_eq!(Color::parse_hex("#12345"), None);
    }

    #[test]
    fn feature_span_through_origin() {
        let mut feature = Feature::new("x", "CDS", Span::new(90, 100), Strand::Forward);
        feature.segments.push(Segment::new(Span::new(0, 10)));
        assert_eq!(feature.span(100), Span::new(90, 110));
        assert_eq!(feature.len(), 20);
    }

    #[test]
    fn normalize_rotates_segments_to_feature_start() {
        let mut feature = Feature::new("x", "CDS", Span::new(0, 10), Strand::Forward);
        feature.segments.push(Segment::new(Span::new(90, 100)));
        feature.normalize(100, true);
        assert_eq!(feature.segments[0].span, Span::new(90, 100));
        assert_eq!(feature.span(100), Span::new(90, 110));
    }

    #[test]
    fn extract_reverse_feature() {
        let seq = b"AAAACCCGGG";
        let feature = Feature::new("x", "misc_feature", Span::new(4, 10), Strand::Reverse);
        assert_eq!(feature.extract(seq), b"CCCGGG".to_vec());
        let forward = Feature::new("y", "misc_feature", Span::new(0, 5), Strand::Forward);
        assert_eq!(forward.extract(seq), b"AAAAC".to_vec());
    }
}
