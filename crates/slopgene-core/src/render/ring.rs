//! Small ring maps for the protocol bench.
//!
//! A ring draws a molecule as a thick circle, or a bar when it is linear, in a few colours: green
//! where the bases go on into the product, pale grey where they are left over, and the backbone's
//! grey where nothing is decided yet. A cut molecule is drawn as its pieces, each a stretch of the ring with a small
//! gap at every cut. A piece can be lifted away from the centre, so the pieces read as separate
//! molecules that came out of the one ring. The part being worked on is turned to 12 o'clock.
//! Features lie just inside the ring as thin arcs in their usual colours.
//!
//! The ring's centre is the middle of the canvas and its radius is [`radius`], so the app can
//! line a ring up with a full map of the same molecule when it zooms from one to the other.

use super::circular::fit_text;
use super::scene::*;
use super::{ColorMode, SiteMark, Theme, assign_tracks, drawable_features, feature_rgba};
use crate::model::{Document, Strand};
use crate::span::Span;
use std::f32::consts::TAU;

/// What a stretch of the ring stands for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    /// Bases that go on into the product: green.
    Kept,
    /// Bases left over: pale grey.
    LeftOver,
    /// Nothing decided yet: the backbone's grey.
    Plain,
}

/// One stretch of the ring.
#[derive(Clone, Debug, PartialEq)]
pub struct RingPiece {
    pub span: Span,
    pub role: Role,
    /// How far the piece is pulled away from the centre, in pixels, along the middle of its arc
    /// (upward for a linear molecule).
    pub lift: f32,
    /// A label beside the piece, such as its length.
    pub label: Option<String>,
}

/// What to draw and how.
#[derive(Clone, Debug)]
pub struct RingOptions {
    pub width: f32,
    pub height: f32,
    pub theme: Theme,
    pub colors: ColorMode,
    /// Label text size in pixels.
    pub font_size: f32,
    /// The base drawn at 12 o'clock on a circle.
    pub top: usize,
    /// The stretches of the ring. Without any, the whole ring is drawn plain.
    pub pieces: Vec<RingPiece>,
    /// Where the ring is cut or joined, each with a name such as an enzyme's.
    pub marks: Vec<SiteMark>,
    /// Whether the marks' names are written beside the ring.
    pub mark_labels: bool,
    /// Whether each piece stops short of its ends, so a cut shows as a gap.
    pub gaps: bool,
    /// Whether the features are drawn inside the ring.
    pub features: bool,
    /// Text in the middle of a circle, or above a bar: a name, then a smaller line.
    pub title: Option<String>,
    pub subtitle: Option<String>,
}

impl Default for RingOptions {
    fn default() -> Self {
        RingOptions {
            width: 220.0,
            height: 200.0,
            theme: Theme::light(),
            colors: ColorMode::Palette,
            font_size: 11.0,
            top: 0,
            pieces: Vec::new(),
            marks: Vec::new(),
            mark_labels: false,
            gaps: false,
            features: true,
            title: None,
            subtitle: None,
        }
    }
}

/// The share of a canvas's smaller side that a ring's radius takes.
pub const RADIUS_SHARE: f32 = 0.34;

/// The radius of the ring on a canvas of this size. A caller that turns a ring into a full map
/// gives the ring a canvas whose radius matches the map's circle.
pub fn radius(width: f32, height: f32) -> f32 {
    (width.min(height) * RADIUS_SHARE).max(14.0)
}

/// The colour of a stretch of this role.
pub fn role_color(role: Role, theme: &Theme) -> Rgba {
    match role {
        Role::Kept => theme.kept,
        Role::LeftOver => theme.left_over,
        Role::Plain => theme.backbone,
    }
}

/// Draw the molecule as a ring: a circle for a circular molecule, a bar for a linear one.
pub fn render_ring(doc: &Document, options: &RingOptions, measure: &dyn TextMeasure) -> Scene {
    if doc.is_circular() { circle(doc, options, measure) } else { bar(doc, options, measure) }
}

/// A label waiting for a free place beside the ring.
struct Pending {
    text: String,
    /// Where the label wants to sit, and the direction it moves in when that place is taken.
    pos: [f32; 2],
    away: [f32; 2],
    anchor: Anchor,
    size: f32,
    color: Rgba,
    bold: bool,
    /// A label that must show moves as far as it needs; another one is dropped after a few tries.
    required: bool,
    tag: Tag,
}

/// Place labels one after another, moving each away from the ring until it overlaps none placed
/// before it. Labels that are not required are dropped when three moves do not free them.
fn place(scene: &mut Scene, labels: Vec<Pending>, measure: &dyn TextMeasure) {
    let mut taken: Vec<[f32; 4]> = Vec::new();
    for label in labels {
        let width = measure.text_width(&label.text, label.size, label.bold);
        let rect_at = |p: [f32; 2]| {
            let x0 = match label.anchor {
                Anchor::Start => p[0],
                Anchor::Middle => p[0] - width / 2.0,
                Anchor::End => p[0] - width,
            };
            [x0 - 2.0, p[1] - label.size * 0.65, x0 + width + 2.0, p[1] + label.size * 0.65]
        };
        let free = |r: &[f32; 4]| !taken.iter().any(|t| r[0] < t[2] && t[0] < r[2] && r[1] < t[3] && t[1] < r[3]);
        let mut pos = label.pos;
        let mut tries = 0;
        while !free(&rect_at(pos)) && (tries < 3 || label.required && tries < 12) {
            pos = [pos[0] + label.away[0] * label.size * 1.2, pos[1] + label.away[1] * label.size * 1.2];
            tries += 1;
        }
        let rect = rect_at(pos);
        if !free(&rect) && !label.required {
            continue;
        }
        taken.push(rect);
        scene.push(Shape::Text { pos, text: label.text, size: label.size, color: label.color, anchor: label.anchor, bold: label.bold }, label.tag);
    }
}

/// How a label beside a circle at angle `a` is anchored, so it reads away from the circle.
fn anchor_at(a: f32) -> Anchor {
    let s = a.sin();
    if s > 0.3 {
        Anchor::Start
    } else if s < -0.3 {
        Anchor::End
    } else {
        Anchor::Middle
    }
}

/// The stretches of `inner` that lie inside `outer`, in bases from `outer`'s start. Either may
/// run through the origin of a molecule of length `len`.
fn within(inner: Span, outer: Span, len: usize) -> Vec<(usize, usize)> {
    let rel = (inner.start % len + len - outer.start % len) % len;
    let mut out = Vec::new();
    // The stretch as it starts after `outer`'s start, and as it was before it, once round.
    for start in [rel as i64, rel as i64 - len as i64] {
        let end = start + inner.len() as i64;
        let (a, b) = (start.max(0), end.min(outer.len() as i64));
        if a < b {
            out.push((a as usize, b as usize));
        }
    }
    out
}

/// A circular molecule as a ring.
fn circle(doc: &Document, options: &RingOptions, measure: &dyn TextMeasure) -> Scene {
    let theme = &options.theme;
    let (w, h) = (options.width, options.height);
    let mut scene = Scene::new(w, h, theme.background);
    let len = doc.len().max(1);
    let center = [w / 2.0, h / 2.0];
    let r = radius(w, h);
    let band = (r * 0.09).clamp(3.5, 14.0);
    let font = options.font_size;
    let top = (options.top % len) as f32;
    let angle_of = |pos: f32| (pos - top) / len as f32 * TAU;
    // Each piece stops this far short of a cut, about 3 pixels of arc on each side.
    let gap = if options.gaps { (band * 0.4 + 2.5) / r } else { 0.0 };
    let pieces: Vec<RingPiece> =
        if options.pieces.is_empty() { vec![RingPiece { span: Span::new(0, len), role: Role::Plain, lift: 0.0, label: None }] } else { options.pieces.clone() };

    // Features take thin tracks just inside the band, at most three, the longest outermost.
    let features = if options.features { drawable_features(doc) } else { Vec::new() };
    let spans: Vec<Span> = features.iter().map(|f| f.span(len).normalize(len)).collect();
    let tracks = assign_tracks(&spans, len, len / 300);
    let fine = (band * 0.55).max(2.0);
    let track_radius = |track: usize| r - band / 2.0 - 2.5 - fine / 2.0 - track as f32 * (fine + 1.5);

    let mut labels: Vec<Pending> = Vec::new();
    for (ii, piece) in pieces.iter().enumerate() {
        let a0 = angle_of(piece.span.start as f32);
        let a1 = angle_of((piece.span.start + piece.span.len()) as f32);
        let mid = (a0 + a1) / 2.0;
        // A lifted piece moves out along the middle of its arc.
        let c = [center[0] + piece.lift * mid.sin(), center[1] - piece.lift * mid.cos()];
        let color = role_color(piece.role, theme);
        let tag = Tag::Piece(ii);
        if piece.span.len() >= len && gap == 0.0 {
            scene.push(Shape::Circle { center: c, radius: r, fill: None, stroke: Some((color, band)) }, tag.clone());
        } else if a1 - gap > a0 + gap {
            scene.push(Shape::Polyline { points: arc_points(c, r, a0 + gap, a1 - gap), color, width: band }, tag.clone());
        }
        // The features on this piece, faded on a piece left over.
        for (iii, feature) in features.iter().enumerate() {
            if tracks[iii] > 2 {
                continue;
            }
            let mut fill = feature_rgba(feature, options.colors);
            if piece.role == Role::LeftOver {
                fill = fill.mix(theme.background, 0.55);
            }
            // Each stretch of the feature on this piece is a thin arc. The stretch that holds the
            // feature's leading end gets an arrowhead, so the direction a gene reads shows, and an
            // insert turned round reads the other way.
            let tr = track_radius(tracks[iii]);
            let head = fine * 1.7 / tr.max(1.0);
            let feature_rel = (spans[iii].start % len + len - piece.span.start % len) % len;
            for (from, to) in within(spans[iii], piece.span, len) {
                let (mut f0, mut f1) = (angle_of((piece.span.start + from) as f32).max(a0 + gap), angle_of((piece.span.start + to) as f32).min(a1 - gap));
                if f1 <= f0 {
                    continue;
                }
                let starts_here = from == feature_rel;
                let ends_here = to == feature_rel + spans[iii].len();
                let forward = matches!(feature.strand, Strand::Forward) && ends_here;
                let reverse = matches!(feature.strand, Strand::Reverse) && starts_here;
                let tag = Tag::Feature(feature.uid);
                if forward && f1 - f0 > head * 1.5 {
                    f1 -= head;
                    let tip = f1 + head;
                    let points = vec![polar(c, tr + fine * 0.95, f1), polar(c, tr - fine * 0.95, f1), polar(c, tr, tip)];
                    scene.push(Shape::Polygon { points, fill, stroke: None }, tag.clone());
                }
                if reverse && f1 - f0 > head * 1.5 {
                    f0 += head;
                    let tip = f0 - head;
                    let points = vec![polar(c, tr + fine * 0.95, f0), polar(c, tr - fine * 0.95, f0), polar(c, tr, tip)];
                    scene.push(Shape::Polygon { points, fill, stroke: None }, tag.clone());
                }
                scene.push(Shape::Polyline { points: arc_points(c, tr, f0, f1), color: fill, width: fine }, tag);
            }
        }
        if let Some(text) = &piece.label {
            let at = polar(c, r + band / 2.0 + font * 0.9, mid);
            let away = [mid.sin(), -mid.cos()];
            let ink = if piece.role == Role::LeftOver { theme.muted } else { theme.ink };
            labels.push(Pending {
                text: text.clone(),
                pos: at,
                away,
                anchor: anchor_at(mid),
                size: font,
                color: ink,
                bold: piece.role == Role::Kept,
                required: true,
                tag,
            });
        }
    }
    // The names of the cuts or joins, just inside the ring where it breaks. Outside, the pieces
    // pulled away from the ring would cross them. Inside the features' tracks, a name reads
    // toward the centre and moves further in when its place is taken.
    if options.mark_labels {
        let inner = track_radius(2) - fine / 2.0 - font * 0.9;
        for mark in &options.marks {
            let a = angle_of(mark.pos as f32);
            let anchor = match anchor_at(a) {
                Anchor::Start => Anchor::End,
                Anchor::End => Anchor::Start,
                Anchor::Middle => Anchor::Middle,
            };
            labels.push(Pending {
                text: mark.name.clone(),
                pos: polar(center, inner, a),
                away: [-a.sin(), a.cos()],
                anchor,
                size: font * 0.92,
                color: theme.enzyme,
                bold: false,
                required: false,
                tag: Tag::Enzyme(mark.name.clone(), mark.pos),
            });
        }
    }
    place(&mut scene, labels, measure);

    // The name and a second line in the middle, shortened to fit inside the ring.
    let room = (r - band) * 1.6;
    if let Some(title) = &options.title {
        let text = fit_text(title, font * 1.05, true, room, measure);
        let y = if options.subtitle.is_some() { center[1] - font * 0.6 } else { center[1] };
        scene.push(Shape::Text { pos: [center[0], y], text, size: font * 1.05, color: theme.ink, anchor: Anchor::Middle, bold: true }, Tag::None);
    }
    if let Some(subtitle) = &options.subtitle {
        let text = fit_text(subtitle, font * 0.92, false, room, measure);
        let y = if options.title.is_some() { center[1] + font * 0.75 } else { center[1] };
        scene.push(Shape::Text { pos: [center[0], y], text, size: font * 0.92, color: theme.muted, anchor: Anchor::Middle, bold: false }, Tag::None);
    }
    scene
}

/// A linear molecule as a bar across the canvas, its pieces lifted upward.
fn bar(doc: &Document, options: &RingOptions, measure: &dyn TextMeasure) -> Scene {
    let theme = &options.theme;
    let (w, h) = (options.width, options.height);
    let mut scene = Scene::new(w, h, theme.background);
    let len = doc.len().max(1);
    let font = options.font_size;
    let band = (h * 0.045).clamp(3.5, 12.0);
    let (x0, x1) = (font, w - font);
    let y = h * 0.58;
    let x_of = |pos: usize| x0 + (x1 - x0) * pos as f32 / len as f32;
    let gap = if options.gaps { band * 0.4 + 2.5 } else { 0.0 };
    let pieces: Vec<RingPiece> =
        if options.pieces.is_empty() { vec![RingPiece { span: Span::new(0, len), role: Role::Plain, lift: 0.0, label: None }] } else { options.pieces.clone() };
    let features = if options.features { drawable_features(doc) } else { Vec::new() };
    let fine = (band * 0.55).max(2.0);

    let mut labels: Vec<Pending> = Vec::new();
    for (ii, piece) in pieces.iter().enumerate() {
        let (start, end) = (piece.span.start.min(len), (piece.span.start + piece.span.len()).min(len));
        let (a, b) = (x_of(start) + gap, x_of(end) - gap);
        let yy = y - piece.lift;
        let color = role_color(piece.role, theme);
        let tag = Tag::Piece(ii);
        if b > a {
            scene.push(Shape::Polyline { points: vec![[a, yy], [b, yy]], color, width: band }, tag.clone());
        }
        // The features on this piece, in one thin line under it.
        for feature in &features {
            let span = feature.span(len);
            let (f0, f1) = (x_of(span.start.max(start)).max(a), x_of(span.end.min(end)).min(b));
            if f1 > f0 {
                let mut fill = feature_rgba(feature, options.colors);
                if piece.role == Role::LeftOver {
                    fill = fill.mix(theme.background, 0.55);
                }
                let fy = yy + band / 2.0 + 2.5 + fine / 2.0;
                scene.push(Shape::Polyline { points: vec![[f0, fy], [f1, fy]], color: fill, width: fine }, Tag::Feature(feature.uid));
            }
        }
        if let Some(text) = &piece.label {
            let ink = if piece.role == Role::LeftOver { theme.muted } else { theme.ink };
            labels.push(Pending {
                text: text.clone(),
                pos: [(a + b) / 2.0, yy - band / 2.0 - font * 0.9],
                away: [0.0, -1.0],
                anchor: Anchor::Middle,
                size: font,
                color: ink,
                bold: piece.role == Role::Kept,
                required: true,
                tag,
            });
        }
    }
    if options.mark_labels {
        for mark in &options.marks {
            labels.push(Pending {
                text: mark.name.clone(),
                pos: [x_of(mark.pos.min(len)), y + band / 2.0 + font * 1.9],
                away: [0.0, 1.0],
                anchor: Anchor::Middle,
                size: font * 0.92,
                color: theme.enzyme,
                bold: false,
                required: false,
                tag: Tag::Enzyme(mark.name.clone(), mark.pos),
            });
        }
    }
    place(&mut scene, labels, measure);
    // The name and a second line above the bar.
    if let Some(title) = &options.title {
        let text = fit_text(title, font * 1.05, true, w - 2.0 * font, measure);
        scene.push(Shape::Text { pos: [w / 2.0, font * 1.1], text, size: font * 1.05, color: theme.ink, anchor: Anchor::Middle, bold: true }, Tag::None);
    }
    if let Some(subtitle) = &options.subtitle {
        let text = fit_text(subtitle, font * 0.92, false, w - 2.0 * font, measure);
        scene.push(Shape::Text { pos: [w / 2.0, font * 2.5], text, size: font * 0.92, color: theme.muted, anchor: Anchor::Middle, bold: false }, Tag::None);
    }
    scene
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Feature, Strand, Topology};
    use crate::render::{ApproxMeasure, svg};

    fn plasmid() -> Document {
        let mut doc = Document::new_dna("pDemo", "ACGT".repeat(1000), Topology::Circular);
        doc.features.push(Feature::new("geneX", "CDS", Span::new(1000, 1500), Strand::Forward));
        doc.features.push(Feature::new("AmpR", "CDS", Span::new(2500, 3300), Strand::Reverse));
        doc
    }

    /// How many shapes carry each tag.
    fn count(scene: &Scene, tag: &Tag) -> usize {
        scene.items.iter().filter(|item| item.tag == *tag).count()
    }

    #[test]
    fn a_plain_ring_draws_the_whole_molecule_and_its_features() {
        let doc = plasmid();
        let options = RingOptions { title: Some("pDemo".into()), subtitle: Some("4,000 bp".into()), ..RingOptions::default() };
        let scene = render_ring(&doc, &options, &ApproxMeasure);
        assert_eq!(count(&scene, &Tag::Piece(0)), 1, "one whole ring");
        assert!(scene.items.iter().any(|i| matches!(&i.shape, Shape::Circle { stroke: Some((c, _)), .. } if *c == options.theme.backbone)));
        assert!(doc.features.iter().all(|f| count(&scene, &Tag::Feature(f.uid)) >= 1), "each feature is drawn inside the ring");
        // A gene on the forward strand ends in an arrowhead.
        let gene = doc.features[0].uid;
        assert!(scene.items.iter().any(|i| i.tag == Tag::Feature(gene) && matches!(i.shape, Shape::Polygon { .. })));
        let svg = svg::to_svg(&scene);
        assert!(svg.contains("pDemo") && svg.contains("4,000 bp"));
    }

    #[test]
    fn cut_pieces_are_coloured_by_role_and_pulled_apart() {
        let doc = plasmid();
        // geneX's piece is kept and lifted; the rest is left over and stays.
        let pieces = vec![
            RingPiece { span: Span::new(950, 1550), role: Role::Kept, lift: 14.0, label: Some("geneX · 600 bp".into()) },
            RingPiece { span: Span::new(1550, 4950), role: Role::LeftOver, lift: 0.0, label: Some("3,400 bp".into()) },
        ];
        let marks = vec![SiteMark { name: "EcoRI".into(), pos: 950 }, SiteMark { name: "BamHI".into(), pos: 1550 }];
        let options = RingOptions { top: 1250, pieces, marks, mark_labels: true, gaps: true, ..RingOptions::default() };
        let scene = render_ring(&doc, &options, &ApproxMeasure);
        let arc_of = |tag: &Tag| {
            scene
                .items
                .iter()
                .find_map(|i| match (&i.shape, &i.tag) {
                    (Shape::Polyline { points, color, .. }, t) if t == tag => Some((points.clone(), *color)),
                    _ => None,
                })
                .unwrap()
        };
        let (kept, kept_color) = arc_of(&Tag::Piece(0));
        let (rest, rest_color) = arc_of(&Tag::Piece(1));
        assert_eq!((kept_color, rest_color), (options.theme.kept, options.theme.left_over));
        // The kept piece sits at 12 o'clock and above where the ring would be: lifted upward.
        let top_y = kept.iter().map(|p| p[1]).fold(f32::MAX, f32::min);
        let ring_top = options.height / 2.0 - radius(options.width, options.height);
        assert!(top_y < ring_top - 10.0, "{top_y} vs {ring_top}");
        // The rest stops short of the cuts on both sides, leaving gaps.
        assert!(rest.len() > 2);
        let svg = svg::to_svg(&scene);
        assert!(svg.contains("geneX · 600 bp") && svg.contains("3,400 bp"));
        assert!(svg.contains("EcoRI") && svg.contains("BamHI"));
    }

    #[test]
    fn a_linear_molecule_is_a_bar() {
        let mut doc = plasmid();
        doc.topology = Topology::Linear;
        let pieces = vec![
            RingPiece { span: Span::new(0, 1000), role: Role::LeftOver, lift: 0.0, label: None },
            RingPiece { span: Span::new(1000, 4000), role: Role::Kept, lift: 8.0, label: Some("kept".into()) },
        ];
        let scene = render_ring(&doc, &RingOptions { pieces, gaps: true, ..RingOptions::default() }, &ApproxMeasure);
        assert_eq!(count(&scene, &Tag::Piece(0)), 1);
        assert!(svg::to_svg(&scene).contains("kept"));
    }

    #[test]
    fn a_feature_is_clipped_to_each_piece_it_lies_on() {
        // A stretch through the origin, and one that runs past the piece's end.
        assert_eq!(within(Span::new(90, 110), Span::new(80, 100), 100), vec![(10, 20)]);
        assert_eq!(within(Span::new(5, 15), Span::new(80, 110), 100), vec![(25, 30)]);
        assert_eq!(within(Span::new(20, 30), Span::new(80, 110), 100), Vec::<(usize, usize)>::new());
        assert_eq!(within(Span::new(10, 20), Span::new(0, 100), 100), vec![(10, 20)]);
    }
}
