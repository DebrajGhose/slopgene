//! Plasmid and sequence maps.
//!
//! Layout code turns a document into a [`Scene`] of simple shapes. Back ends draw the scene:
//! [`svg::to_svg`] for files, and the desktop app paints the same scene with egui, so maps look
//! the same everywhere.

pub mod bench;
pub mod circular;
pub mod linear;
pub mod ring;
pub mod scene;
pub mod svg;
pub mod theme;

pub use scene::*;
pub use theme::Theme;

use crate::enzymes::{EnzymeDb, find_sites_in};
use crate::model::{Document, Feature};
use crate::span::Span;

/// An enzyme site to mark on a map.
#[derive(Clone, Debug, PartialEq)]
pub struct SiteMark {
    pub name: String,
    /// Top-strand cut position (a boundary, 0..=length).
    pub pos: usize,
}

/// Where feature colours come from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ColorMode {
    /// One calm colour per feature type, from the SlopGene palette.
    #[default]
    Palette,
    /// The colours stored in the file (SnapGene's own), falling back to the palette.
    File,
}

impl ColorMode {
    pub fn parse(text: &str) -> Option<ColorMode> {
        match text.trim().to_ascii_lowercase().as_str() {
            "palette" | "type" | "slopgene" => Some(ColorMode::Palette),
            "file" | "snapgene" | "original" => Some(ColorMode::File),
            _ => None,
        }
    }
}

/// What kind of difference between two aligned sequences a mark stands for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DifferenceKind {
    /// Bases that differ: substitutions and inversions.
    Changed,
    /// Bases only the aligned sequence has: insertions, duplications, and the sequence that
    /// stands in place of replaced bases.
    Inserted,
    /// Reference bases the aligned sequence lacks: deletions, and the bases of a replacement.
    Missing,
}

impl DifferenceKind {
    /// The kinds in the order a legend lists them, after the aligned stretches.
    pub const ALL: [DifferenceKind; 3] = [DifferenceKind::Missing, DifferenceKind::Changed, DifferenceKind::Inserted];

    /// What the colour means, for a legend.
    pub fn legend(self) -> &'static str {
        match self {
            DifferenceKind::Changed => "different bases",
            DifferenceKind::Inserted => "bases only the aligned sequence has",
            DifferenceKind::Missing => "bases the aligned sequence lacks",
        }
    }
}

/// What an aligned sequence has over a stretch of the reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SegmentKind {
    /// The same bases in the same order. Changes of a few bases are marks on the stretch.
    Aligned,
    /// The same bases turned round (an inversion).
    Inverted,
    /// None of these bases: they were deleted, or replaced by other sequence.
    Missing,
}

/// One stretch of the reference as an aligned sequence has it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrackSegment {
    pub kind: SegmentKind,
    /// The reference bases of the stretch (0-based, may run through the origin).
    pub span: Span,
}

/// One difference drawn on a map.
#[derive(Clone, Debug, PartialEq)]
pub struct DifferenceMark {
    pub kind: DifferenceKind,
    /// The reference bases involved (0-based). Empty at the point where bases were inserted.
    pub span: Span,
    /// A short label, such as "501 T>A" or "2,001..3,500 missing".
    pub label: String,
    /// One line for a tooltip.
    pub description: String,
    /// A structural change (50 bp or more) rather than a change of a few bases.
    pub structural: bool,
}

/// An aligned sequence drawn on a map: a ring around a circular map, or a line under a linear
/// one. The ring is solid where the sequence aligns, thin where it lacks the reference's bases,
/// and absent where it gives no information (outside a partial read). A mark sits at each
/// difference.
#[derive(Clone, Debug, PartialEq)]
pub struct AlignmentTrack {
    pub name: String,
    /// The reference stretch the sequence covers.
    pub covered: Span,
    /// The covered stretch cut where the sequence lacks the reference's bases or has them turned
    /// round, in order along it. Together the segments make up `covered`.
    pub segments: Vec<TrackSegment>,
    pub marks: Vec<DifferenceMark>,
    /// The track shown in full colour and with labels. Other tracks are drawn muted.
    pub emphasis: bool,
}

impl AlignmentTrack {
    /// A track for a sequence that does not align: nothing to draw.
    pub fn empty(name: &str, emphasis: bool) -> AlignmentTrack {
        AlignmentTrack { name: name.to_string(), covered: Span::new(0, 0), segments: Vec::new(), marks: Vec::new(), emphasis }
    }

    /// A track that aligns over all of `covered`, with these marks.
    pub fn aligned(name: &str, covered: Span, marks: Vec<DifferenceMark>, emphasis: bool) -> AlignmentTrack {
        let segments = if covered.is_empty() { Vec::new() } else { vec![TrackSegment { kind: SegmentKind::Aligned, span: covered }] };
        AlignmentTrack { name: name.to_string(), covered, segments, marks, emphasis }
    }
}

/// What to draw and how.
#[derive(Clone, Debug)]
pub struct MapOptions {
    pub width: f32,
    pub height: f32,
    pub theme: Theme,
    /// Label text size in pixels.
    pub font_size: f32,
    pub show_features: bool,
    pub show_primers: bool,
    pub show_ticks: bool,
    pub enzymes: Vec<SiteMark>,
    /// Enzyme marks drawn in the accent colour and in bold, such as the cuts a protocol chose.
    /// They are drawn whether or not `enzymes` lists them, and their labels stay when space runs
    /// out.
    pub chosen_enzymes: Vec<SiteMark>,
    /// Enzyme marks for sites that methylation blocks, drawn faint at the cut each would make.
    pub blocked_enzymes: Vec<SiteMark>,
    pub selection: Option<Span>,
    pub cursor: Option<usize>,
    /// Feature to outline (by uid).
    pub highlight: Option<u64>,
    /// Title in the middle of a circular map; defaults to the document name.
    pub title: Option<String>,
    /// The base drawn at 12 o'clock on a circular map, so the part being worked on can sit at
    /// the top. 0 keeps the origin there.
    pub top: usize,
    /// At most this many feature labels (the longest features win).
    pub max_feature_labels: usize,
    /// At most this many labels in each label column of a circular map. Primers are dropped
    /// first, then enzymes, then the shortest features.
    pub max_labels_per_side: usize,
    pub colors: ColorMode,
    /// Aligned sequences to draw with their differences (see [`AlignmentTrack`]).
    pub alignments: Vec<AlignmentTrack>,
}

impl Default for MapOptions {
    fn default() -> Self {
        MapOptions {
            width: 900.0,
            height: 700.0,
            theme: Theme::light(),
            font_size: 11.0,
            show_features: true,
            show_primers: true,
            show_ticks: true,
            enzymes: Vec::new(),
            chosen_enzymes: Vec::new(),
            blocked_enzymes: Vec::new(),
            selection: None,
            cursor: None,
            highlight: None,
            title: None,
            top: 0,
            max_feature_labels: 60,
            max_labels_per_side: 30,
            colors: ColorMode::Palette,
            alignments: Vec::new(),
        }
    }
}

/// Unique cutters among `names` (or the common set when empty), as site marks for a map.
pub fn unique_cutter_marks(doc: &Document, enzymes: &[&crate::enzymes::Enzyme]) -> Vec<SiteMark> {
    let sites = find_sites_in(doc, enzymes);
    let groups = crate::enzymes::cutting_sites_by_enzyme(&sites);
    let mut marks = Vec::new();
    for enzyme in enzymes {
        let cutting: &[&crate::enzymes::Site] = groups.get(enzyme.name.as_str()).map_or(&[], |g| g.as_slice());
        if cutting.len() == 1 && cutting[0].cuts.len() == 1 {
            marks.push(SiteMark { name: enzyme.name.clone(), pos: cutting[0].cuts[0].pos });
        }
    }
    marks.sort_by_key(|m| m.pos);
    marks
}

/// Marks for the sites of `enzymes` that methylation blocks, at the cut each would have made.
pub fn blocked_marks(doc: &Document, enzymes: &[&crate::enzymes::Enzyme]) -> Vec<SiteMark> {
    find_sites_in(doc, enzymes)
        .iter()
        .filter(|s| s.blocked.is_some())
        .filter_map(|s| s.blocked_cuts.first().map(|c| SiteMark { name: s.enzyme.clone(), pos: c.pos }))
        .collect()
}

/// Site marks for the default map: unique cutters from the common enzyme set.
pub fn default_site_marks(doc: &Document) -> Vec<SiteMark> {
    let db = EnzymeDb::builtin();
    let enzymes = db.set("common").unwrap_or_default();
    unique_cutter_marks(doc, &enzymes)
}

/// Every enzyme mark a map draws, each with whether it is chosen: the marks of `enzymes`, then
/// the chosen marks that `enzymes` does not list.
pub(crate) fn all_marks(options: &MapOptions) -> Vec<(SiteMark, MarkStyle)> {
    let style = |m: &SiteMark| if options.chosen_enzymes.contains(m) { MarkStyle::Chosen } else { MarkStyle::Plain };
    let mut marks: Vec<(SiteMark, MarkStyle)> = options.enzymes.iter().map(|m| (m.clone(), style(m))).collect();
    for mark in &options.chosen_enzymes {
        if !options.enzymes.contains(mark) {
            marks.push((mark.clone(), MarkStyle::Chosen));
        }
    }
    marks.extend(options.blocked_enzymes.iter().map(|m| (m.clone(), MarkStyle::Blocked)));
    marks
}

/// How a map draws an enzyme mark.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MarkStyle {
    /// A cut.
    Plain,
    /// A cut a protocol chose: the accent colour, bold, and its label outranks the rest.
    Chosen,
    /// A site that methylation blocks: faint, at the position the cut would have.
    Blocked,
}

/// Draw a map of the document: circular for circular molecules, linear otherwise.
pub fn render_map(doc: &Document, options: &MapOptions, measure: &dyn TextMeasure) -> Scene {
    if doc.is_circular() { circular::render(doc, options, measure) } else { linear::render(doc, options, measure) }
}

/// The colour a feature is drawn in.
pub fn feature_rgba(feature: &Feature, mode: ColorMode) -> Rgba {
    match mode {
        ColorMode::Palette => theme::feature_color(&feature.kind),
        ColorMode::File => feature.color().map(Rgba::from).unwrap_or_else(|| theme::feature_color(&feature.kind)),
    }
}

/// A darker edge for colours too light to stand out on a light background.
pub fn edge_for(color: Rgba, theme: &Theme) -> Option<Rgba> {
    (!theme.dark && color.luminance() > 0.6).then(|| color.mix(theme::INK, 0.45))
}

/// Features worth drawing on a map (everything except whole-record `source` features).
pub fn drawable_features(doc: &Document) -> Vec<&Feature> {
    doc.features.iter().filter(|f| !f.kind.eq_ignore_ascii_case("source") && !f.segments.is_empty()).collect()
}

/// A tick spacing that gives roughly `target` ticks along `len` bases: 1, 2 or 5 times a power
/// of ten.
pub fn nice_tick_step(len: usize, target: usize) -> usize {
    let raw = (len as f64 / target.max(1) as f64).max(1.0);
    let magnitude = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0].iter().map(|m| m * magnitude).find(|s| *s >= raw).unwrap_or(10.0 * magnitude);
    step.max(1.0) as usize
}

/// Format a base count with thousands separators: 12345 -> "12,345".
pub fn thousands(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (ii, ch) in digits.chars().enumerate() {
        if ii > 0 && (digits.len() - ii).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

/// Assign each span to the lowest track where it overlaps nothing already placed.
/// `pad` widens spans by that many bases when testing overlap, so arrows do not touch.
pub fn assign_tracks(spans: &[Span], seq_len: usize, pad: usize) -> Vec<usize> {
    // Place long features first so they sit on the backbone; short ones stack around them.
    let mut order: Vec<usize> = (0..spans.len()).collect();
    order.sort_by_key(|&ii| std::cmp::Reverse(spans[ii].len()));
    let mut tracks: Vec<Vec<Span>> = Vec::new();
    let mut assigned = vec![0usize; spans.len()];
    for ii in order {
        let padded = Span::new(spans[ii].start.saturating_sub(pad), spans[ii].end + pad);
        let mut chosen = None;
        for (iii, track) in tracks.iter().enumerate() {
            if !track.iter().any(|s| s.overlaps(&padded, seq_len.max(1))) {
                chosen = Some(iii);
                break;
            }
        }
        let track = chosen.unwrap_or_else(|| {
            tracks.push(Vec::new());
            tracks.len() - 1
        });
        tracks[track].push(padded);
        assigned[ii] = track;
    }
    assigned
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Strand, Topology};

    #[test]
    fn tick_steps_and_separators() {
        assert_eq!(nice_tick_step(2686, 8), 500);
        assert_eq!(nice_tick_step(10_000, 10), 1000);
        assert_eq!(nice_tick_step(40, 8), 5);
        assert_eq!(thousands(2686), "2,686");
        assert_eq!(thousands(1_234_567), "1,234,567");
        assert_eq!(thousands(12), "12");
    }

    #[test]
    fn overlapping_spans_go_to_different_tracks() {
        let spans = [Span::new(0, 100), Span::new(50, 80), Span::new(200, 300)];
        let tracks = assign_tracks(&spans, 1000, 0);
        assert_eq!(tracks[0], 0);
        assert_eq!(tracks[1], 1);
        assert_eq!(tracks[2], 0);
    }

    #[test]
    fn maps_render_to_svg() {
        let mut doc = Document::new_dna("demo", "ATGC".repeat(500), Topology::Circular);
        doc.features.push(Feature::new("AmpR", "CDS", Span::new(100, 900), Strand::Reverse));
        doc.features.push(Feature::new("ori", "rep_origin", Span::new(1900, 2100), Strand::Forward));
        let options = MapOptions { enzymes: vec![SiteMark { name: "EcoRI".into(), pos: 396 }], ..MapOptions::default() };
        let scene = render_map(&doc, &options, &ApproxMeasure);
        let svg = svg::to_svg(&scene);
        assert!(svg.contains("AmpR"));
        assert!(svg.contains("EcoRI (396)"));
        assert!(svg.contains("2,000 bp"));
        doc.topology = Topology::Linear;
        let linear = svg::to_svg(&render_map(&doc, &options, &ApproxMeasure));
        assert!(linear.contains("AmpR"));
    }

    #[test]
    fn aligned_sequences_are_drawn_with_every_difference() {
        // Two tracks: the emphasised one with a substitution, a missing stretch and an insertion
        // point, and a faded one with no differences.
        let mut doc = Document::new_dna("demo", "ATGC".repeat(500), Topology::Circular);
        doc.features.push(Feature::new("AmpR", "CDS", Span::new(100, 900), Strand::Reverse));
        let mark = |kind, span, label: &str, structural| DifferenceMark { kind, span, label: label.into(), description: String::new(), structural };
        let read = AlignmentTrack {
            name: "read".into(),
            covered: Span::new(0, 2000),
            segments: vec![
                TrackSegment { kind: SegmentKind::Aligned, span: Span::new(0, 1000) },
                TrackSegment { kind: SegmentKind::Missing, span: Span::new(1000, 1400) },
                TrackSegment { kind: SegmentKind::Aligned, span: Span::new(1400, 2000) },
            ],
            marks: vec![
                mark(DifferenceKind::Changed, Span::new(150, 151), "151 A>G", false),
                mark(DifferenceKind::Missing, Span::new(1000, 1400), "1,001..1,400 missing", true),
                mark(DifferenceKind::Inserted, Span::new(1700, 1700), "1,700 +GGTA", false),
            ],
            emphasis: true,
        };
        let other = AlignmentTrack::aligned("other", Span::new(1800, 2300).normalize(2000), Vec::new(), false);
        let options = MapOptions { alignments: vec![read, other], ..MapOptions::default() };
        for topology in [Topology::Circular, Topology::Linear] {
            doc.topology = topology;
            let scene = render_map(&doc, &options, &ApproxMeasure);
            let tagged = |tag: Tag| scene.items.iter().filter(|item| item.tag == tag).count();
            assert!(tagged(Tag::Track(0)) >= 1 && tagged(Tag::Track(1)) >= 1, "{topology:?}: both tracks are drawn");
            for m in 0..3 {
                assert!(tagged(Tag::Difference(0, m)) >= 1, "{topology:?}: mark {m} is drawn");
            }
            let svg = svg::to_svg(&scene);
            assert!(svg.contains("151 A&gt;G") || svg.contains("151 A>G"), "{topology:?}: the substitution is labelled");
            assert!(svg.contains(&theme::CHANGED.hex()), "{topology:?}: in the colour for changed bases");
            // The aligned stretches are green and the missing one red.
            assert!(svg.contains(&theme::ALIGNED.hex()), "{topology:?}: aligned stretches are drawn");
            assert!(svg.contains(&theme::MISSING.hex()), "{topology:?}: the missing stretch is drawn");
        }
    }
}
