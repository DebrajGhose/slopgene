//! Circular plasmid maps.
//!
//! Geometry: base position p sits at angle 2*pi*p/L, clockwise from 12 o'clock. Features are
//! thick arcs on tracks: track 0 is centred on the backbone and further tracks step inward.
//! Labels sit outside the circle in two columns (right and left halves), pushed apart
//! vertically so they never overlap, and joined to their feature by a two-segment leader.
//! Aligned sequences are rings outside the primers: solid where they align, thin where they lack
//! the map's bases, with a mark at each difference.

use super::scene::*;
use super::theme::{difference_color, difference_ink};
use super::{ColorMode, DifferenceKind, MapOptions, SegmentKind, assign_tracks, drawable_features, edge_for, feature_rgba, nice_tick_step, thousands};
use crate::model::{Document, Strand};
use crate::span::Span;
use std::f32::consts::TAU;

/// The most aligned sequences drawn as rings, and the distance between neighbouring rings.
const MAX_RINGS: usize = 6;
const RING_PITCH: f32 = 12.0;
/// At most this many differences are labelled, so the features keep room for their labels.
const MAX_DIFFERENCE_LABELS: usize = 14;

/// A label waiting for a place in one of the two label columns.
struct LabelRequest {
    text: String,
    /// Angle of the thing being labelled.
    angle: f32,
    /// Radius where the leader line starts.
    from_radius: f32,
    color: Rgba,
    tag: Tag,
    /// Higher wins when space runs out.
    priority: f32,
    bold: bool,
}

pub fn render(doc: &Document, options: &MapOptions, measure: &dyn TextMeasure) -> Scene {
    let theme = &options.theme;
    let (w, h) = (options.width, options.height);
    let mut scene = Scene::new(w, h, theme.background);
    let seq_len = doc.len().max(1);
    let center = [w / 2.0, h / 2.0];

    // Sizes scale with the canvas. The backbone leaves room for label columns on both sides.
    // Aligned sequences take rings outside the primers, so the circle shrinks by their width and
    // the label columns stay where they were.
    let rings = options.alignments.len().min(MAX_RINGS);
    let ring_space = if rings == 0 { 0.0 } else { 10.0 + rings as f32 * RING_PITCH };
    let radius = ((h * 0.31).min(w * 0.23) - ring_space).max(40.0);
    let thickness = (radius * 0.085).clamp(7.0, 20.0);
    let track_gap = 3.0;
    let font = options.font_size;
    // Base `top` sits at 12 o'clock, and the rest follow clockwise.
    let top = (options.top % seq_len) as f64;
    let angle_of = |pos: f64| ((pos - top) / seq_len as f64) as f32 * TAU;

    // Feature tracks.
    let features = if options.show_features { drawable_features(doc) } else { Vec::new() };
    let spans: Vec<Span> = features.iter().map(|f| f.span(seq_len)).collect();
    let tracks = assign_tracks(&spans, seq_len, seq_len / 400);
    let track_count = tracks.iter().copied().max().map(|t| t + 1).unwrap_or(1);
    let track_radius = |track: usize| radius - track as f32 * (thickness + track_gap);
    let innermost = track_radius(track_count.saturating_sub(1)) - thickness / 2.0;

    // Selection: a translucent band across every track.
    if let Some(sel) = options.selection.filter(|s| !s.is_empty()) {
        let outer = radius + thickness / 2.0 + 2.0;
        let inner = innermost - 2.0;
        let mid = (outer + inner) / 2.0;
        let a0 = angle_of(sel.start as f64);
        let a1 = angle_of(sel.end as f64);
        scene.push(Shape::Polyline { points: arc_points(center, mid, a0, a1), color: theme.selection, width: outer - inner }, Tag::Selection);
    }

    // Backbone.
    scene.push(Shape::Circle { center, radius, fill: None, stroke: Some((theme.backbone, 1.4)) }, Tag::Backbone);

    // Name and length for the middle. They are measured before the ticks because on a small
    // circle the position numbers inside it would run into them.
    let title = options.title.clone().unwrap_or_else(|| doc.name.clone());
    let title_size = (font * 1.55).min(radius * 0.16);
    let title = fit_text(&title, title_size, true, innermost * 1.5, measure);
    let length_text = format!("{} bp", thousands(doc.len()));
    let length_size = font * 1.05;
    let block_half_width = measure.text_width(&title, title_size, true).max(measure.text_width(&length_text, length_size, false)) / 2.0;
    let block_half_height = title_size * 1.3;

    // Position ticks just inside the innermost track, numbered when the numbers leave the
    // middle clear: the corners of the name block must stay inside the ring of numbers.
    if options.show_ticks {
        let step = nice_tick_step(seq_len, 8);
        let tick_outer = innermost - 4.0;
        let tick_inner = tick_outer - 5.0;
        let number_size = font * 0.82;
        let number_radius = tick_inner - font * 1.3;
        let widest_number = measure.text_width(&thousands(seq_len.saturating_sub(1) / step * step), number_size, false);
        let clear_radius = number_radius - widest_number / 2.0 - 3.0;
        let numbered = block_half_width.hypot(block_half_height) <= clear_radius;
        let mut pos = step;
        while pos < seq_len {
            let a = angle_of(pos as f64);
            scene.push(Shape::Polyline { points: vec![polar(center, tick_outer, a), polar(center, tick_inner, a)], color: theme.faint, width: 1.0 }, Tag::None);
            if numbered {
                let label_pos = polar(center, number_radius, a);
                scene.push(
                    Shape::Text { pos: label_pos, text: thousands(pos), size: number_size, color: theme.muted, anchor: Anchor::Middle, bold: false },
                    Tag::None,
                );
            }
            pos += step;
        }
    }

    // Features, with labels queued for placement.
    let mut labels: Vec<LabelRequest> = Vec::new();
    for (ii, feature) in features.iter().enumerate() {
        let r = track_radius(tracks[ii]);
        let color = feature_rgba(feature, options.colors);
        let tag = Tag::Feature(feature.uid);
        if options.highlight == Some(feature.uid) {
            let span = spans[ii];
            scene.push(
                Shape::Polyline {
                    points: arc_points(center, r, angle_of(span.start as f64), angle_of(span.end as f64)),
                    color: theme.accent,
                    width: thickness + 5.0,
                },
                tag.clone(),
            );
        }
        let segment_count = feature.segments.len();
        for (iii, segment) in feature.segments.iter().enumerate() {
            let a0 = angle_of(segment.span.start as f64);
            let a1 = angle_of(segment.span.end as f64);
            if segment.gap {
                scene.push(Shape::Polyline { points: arc_points(center, r, a0, a1), color, width: 1.5 }, tag.clone());
                continue;
            }
            let seg_color = match options.colors {
                ColorMode::File => segment.color.map(Rgba::from).unwrap_or(color),
                ColorMode::Palette => color,
            };
            // Arrowheads go on the segment at the feature's leading end.
            let head_fwd = matches!(feature.strand, Strand::Forward | Strand::Both) && iii + 1 == segment_count;
            let head_rev = matches!(feature.strand, Strand::Reverse | Strand::Both) && iii == 0;
            let edge = edge_for(seg_color, theme);
            draw_arrow(&mut scene, center, r, thickness, a0, a1, head_fwd, head_rev, seg_color, edge, tag.clone());
        }
        if !feature.name.is_empty() {
            let span = spans[ii];
            let mid = angle_of((span.start as f64 + span.end as f64) / 2.0);
            labels.push(LabelRequest {
                text: feature.name.clone(),
                angle: mid,
                from_radius: r + thickness / 2.0 + 1.0,
                color: theme.ink,
                tag,
                priority: 1000.0 + span.len() as f32 / seq_len as f32,
                bold: false,
            });
        }
    }
    // Keep only the longest features' labels when there are too many.
    if labels.len() > options.max_feature_labels {
        labels.sort_by(|a, b| b.priority.partial_cmp(&a.priority).unwrap_or(std::cmp::Ordering::Equal));
        labels.truncate(options.max_feature_labels);
    }

    // Primers: thin arrows just outside the backbone.
    if options.show_primers {
        let r = radius + thickness / 2.0 + 5.0;
        for primer in &doc.primers {
            for site in primer.sites.iter().filter(|s| s.extra.attr("simplified") != Some("1")) {
                let a0 = angle_of(site.span.start as f64);
                let a1 = angle_of(site.span.end as f64);
                let forward = site.strand != Strand::Reverse;
                draw_arrow(&mut scene, center, r, 3.5, a0, a1, forward, !forward, theme.primer, None, Tag::Primer(primer.uid));
                labels.push(LabelRequest {
                    text: primer.name.clone(),
                    angle: (a0 + a1) / 2.0,
                    from_radius: r + 3.0,
                    color: theme.primer,
                    tag: Tag::Primer(primer.uid),
                    priority: 10.0,
                    bold: false,
                });
            }
        }
    }

    // Aligned sequences: a ring each, outside the map, that shows what the sequence has of it.
    // The ring is a solid green band where the sequence has the map's bases, an orange band where
    // it has them turned round, and a thin red line where it lacks them. Where the sequence gives
    // no information (past the ends of a partial read) there is no ring at all. A change of a few
    // bases is a short radial bar across the band, so even one base shows. Extra bases the map
    // lacks stick out of the ring as a flag. Only the emphasised track is labelled, and the
    // others are drawn faded.
    let first_ring = radius + thickness / 2.0 + 12.0;
    let mut difference_labels: Vec<LabelRequest> = Vec::new();
    for (k, track) in options.alignments.iter().take(rings).enumerate() {
        let r = first_ring + k as f32 * RING_PITCH;
        let fade = |c: Rgba| if track.emphasis { c } else { c.mix(theme.background, 0.55) };
        let band = if track.emphasis { 6.0 } else { 3.5 };
        let thin = 1.5;
        // The arc of a stretch, or a whole circle for a stretch that goes all the way round.
        let stroke = |scene: &mut Scene, span: Span, color: Rgba, width: f32, tag: Tag| {
            if span.len() >= seq_len {
                scene.push(Shape::Circle { center, radius: r, fill: None, stroke: Some((color, width)) }, tag);
            } else if !span.is_empty() {
                let (a0, a1) = (angle_of(span.start as f64), angle_of((span.start + span.len()) as f64));
                scene.push(Shape::Polyline { points: arc_points(center, r, a0, a1), color, width }, tag);
            }
        };
        for segment in &track.segments {
            let (color, width) = match segment.kind {
                SegmentKind::Aligned => (theme.aligned, band),
                SegmentKind::Inverted => (theme.changed, band),
                SegmentKind::Missing => (theme.missing, thin),
            };
            stroke(&mut scene, segment.span, fade(color), width, Tag::Track(k));
        }
        for (m, mark) in track.marks.iter().enumerate() {
            let color = fade(difference_color(mark.kind, theme));
            let tag = Tag::Difference(k, m);
            let (a0, a1) = (angle_of(mark.span.start as f64), angle_of((mark.span.start + mark.span.len()) as f64));
            let a = (a0 + a1) / 2.0;
            let mut label_radius = r + band / 2.0 + 3.0;
            if mark.structural && !mark.span.is_empty() {
                // A structural change's stretch is drawn by its segment. The same arc again carries
                // the change's tag, so pointing at it names the change.
                let width = if mark.kind == DifferenceKind::Missing { thin } else { band };
                stroke(&mut scene, mark.span, color, width, tag.clone());
            } else if mark.structural {
                // Extra sequence: a flag that sticks out of the ring, with a dot at its tip.
                scene.push(Shape::Polyline { points: vec![polar(center, r - band / 2.0, a), polar(center, r + 10.0, a)], color, width: 2.5 }, tag.clone());
                scene.push(Shape::Circle { center: polar(center, r + 12.0, a), radius: 3.0, fill: Some(color), stroke: None }, tag.clone());
                label_radius = r + 16.0;
            } else {
                let reach = band / 2.0 + 3.5;
                scene.push(Shape::Polyline { points: vec![polar(center, r - reach, a), polar(center, r + reach, a)], color, width: 2.5 }, tag.clone());
            }
            if track.emphasis {
                difference_labels.push(LabelRequest {
                    text: mark.label.clone(),
                    angle: a,
                    from_radius: label_radius,
                    color: difference_ink(mark.kind, theme),
                    tag,
                    priority: if mark.structural { 3000.0 } else { 2000.0 },
                    bold: mark.structural,
                });
            }
        }
    }
    // Structural changes are labelled first, then base changes in their order along the map
    // (the sort keeps that order among equals).
    difference_labels.sort_by(|a, b| b.priority.partial_cmp(&a.priority).unwrap_or(std::cmp::Ordering::Equal));
    difference_labels.truncate(MAX_DIFFERENCE_LABELS);
    labels.extend(difference_labels);
    // Label leaders start beyond the outermost ring and the flags that stick out of it.
    let outer_radius = if rings == 0 { radius + thickness / 2.0 } else { first_ring + (rings - 1) as f32 * RING_PITCH + 16.0 };

    // Enzyme sites: a short tick across the backbone and a label. A chosen cut is drawn in the
    // accent colour, longer and bolder, and its label outranks every other.
    for (mark, style) in super::all_marks(options) {
        let a = angle_of(mark.pos as f64);
        let tag = Tag::Enzyme(mark.name.clone(), mark.pos);
        let chosen = style == super::MarkStyle::Chosen;
        let (color, width, reach) = match style {
            super::MarkStyle::Chosen => (theme.accent, 2.4, 7.0),
            super::MarkStyle::Plain => (theme.enzyme, 1.2, 3.0),
            super::MarkStyle::Blocked => (theme.faint, 1.2, 3.0),
        };
        scene.push(
            Shape::Polyline {
                points: vec![polar(center, radius - thickness / 2.0 - 1.0, a), polar(center, radius + thickness / 2.0 + reach, a)],
                color,
                width,
            },
            tag.clone(),
        );
        labels.push(LabelRequest {
            text: format!("{} ({})", mark.name, mark.pos.max(1)),
            angle: a,
            from_radius: radius + thickness / 2.0 + reach,
            color: if style == super::MarkStyle::Blocked { theme.muted } else { color },
            tag,
            priority: match style {
                super::MarkStyle::Chosen => 5000.0,
                super::MarkStyle::Plain => 100.0,
                super::MarkStyle::Blocked => 50.0,
            },
            bold: chosen,
        });
    }

    place_labels(&mut scene, labels, center, outer_radius, font, h, options.max_labels_per_side, measure, theme);

    // Cursor.
    if let Some(cursor) = options.cursor {
        let a = angle_of(cursor as f64);
        scene.push(
            Shape::Polyline {
                points: vec![polar(center, innermost - 2.0, a), polar(center, radius + thickness / 2.0 + 4.0, a)],
                color: theme.accent,
                width: 1.6,
            },
            Tag::None,
        );
    }

    // Name and length in the middle.
    scene.push(
        Shape::Text { pos: [center[0], center[1] - title_size * 0.55], text: title, size: title_size, color: theme.ink, anchor: Anchor::Middle, bold: true },
        Tag::None,
    );
    scene.push(
        Shape::Text {
            pos: [center[0], center[1] + title_size * 0.75],
            text: length_text,
            size: length_size,
            color: theme.muted,
            anchor: Anchor::Middle,
            bold: false,
        },
        Tag::None,
    );
    scene
}

/// A thick arc with optional arrowheads. Short arcs become a single triangle. `edge` outlines
/// colours that would otherwise vanish into the background.
#[allow(clippy::too_many_arguments)]
fn draw_arrow(
    scene: &mut Scene,
    center: [f32; 2],
    r: f32,
    thickness: f32,
    a0: f32,
    a1: f32,
    head_fwd: bool,
    head_rev: bool,
    color: Rgba,
    edge: Option<Rgba>,
    tag: Tag,
) {
    let head_len = thickness * 0.9;
    let head = head_len / r.max(1.0);
    let span = a1 - a0;
    let heads = head_fwd as u8 + head_rev as u8;
    let stroke = edge.map(|e| (e, 1.0));
    if heads > 0 && span <= head * heads as f32 * 1.1 {
        // Too short for a body: draw one triangle pointing the right way.
        let (base, tip) = if head_fwd { (a0, a1) } else { (a1, a0) };
        scene.push(
            Shape::Polygon {
                points: vec![polar(center, r + thickness / 2.0, base), polar(center, r - thickness / 2.0, base), polar(center, r, tip)],
                fill: color,
                stroke,
            },
            tag,
        );
        return;
    }
    let body_start = if head_rev { a0 + head } else { a0 };
    let body_end = if head_fwd { a1 - head } else { a1 };
    if body_end > body_start {
        if let Some(e) = edge {
            // A slightly wider arc underneath reads as a thin outline.
            scene.push(Shape::Polyline { points: arc_points(center, r, body_start, body_end), color: e, width: thickness + 2.0 }, tag.clone());
        }
        scene.push(Shape::Polyline { points: arc_points(center, r, body_start, body_end), color, width: thickness }, tag.clone());
    }
    if head_fwd {
        scene.push(
            Shape::Polygon {
                points: vec![polar(center, r + thickness / 2.0, body_end), polar(center, r - thickness / 2.0, body_end), polar(center, r, a1)],
                fill: color,
                stroke,
            },
            tag.clone(),
        );
    }
    if head_rev {
        scene.push(
            Shape::Polygon {
                points: vec![polar(center, r + thickness / 2.0, body_start), polar(center, r - thickness / 2.0, body_start), polar(center, r, a0)],
                fill: color,
                stroke,
            },
            tag,
        );
    }
}

/// Shorten text with an ellipsis until it fits `max_width`.
pub(crate) fn fit_text(text: &str, size: f32, bold: bool, max_width: f32, measure: &dyn TextMeasure) -> String {
    if measure.text_width(text, size, bold) <= max_width {
        return text.to_string();
    }
    let mut chars: Vec<char> = text.chars().collect();
    while chars.len() > 1 {
        chars.pop();
        let candidate: String = chars.iter().collect::<String>() + "…";
        if measure.text_width(&candidate, size, bold) <= max_width {
            return candidate;
        }
    }
    text.chars().take(1).collect()
}

/// Spread sorted positions so neighbours are at least `gap` apart, moving each as little as
/// possible, within `lo..=hi`.
///
/// Labels that would overlap are merged into a cluster laid out at `gap` spacing and centred on
/// the mean of where its members want to be. Merging repeats until no clusters overlap. This
/// keeps every label near its anchor and preserves their order, so leader lines do not cross.
pub fn spread(desired: &[f32], gap: f32, lo: f32, hi: f32) -> Vec<f32> {
    // Each cluster is (first index, member count, top position).
    let mut clusters: Vec<(usize, usize, f32)> = Vec::new();
    let place = |first: usize, count: usize| -> f32 {
        let mean: f32 = (first..first + count).map(|k| desired[k] - (k - first) as f32 * gap).sum::<f32>() / count as f32;
        let lowest_top = (hi - (count as f32 - 1.0) * gap).max(lo);
        mean.clamp(lo, lowest_top)
    };
    for ii in 0..desired.len() {
        clusters.push((ii, 1, place(ii, 1)));
        // Merge backwards while the newest cluster overlaps the one before it.
        while clusters.len() >= 2 {
            let n = clusters.len();
            let (pf, pc, pt) = clusters[n - 2];
            let (_, cc, ct) = clusters[n - 1];
            if pt + pc as f32 * gap <= ct + 0.01 {
                break;
            }
            clusters.truncate(n - 2);
            clusters.push((pf, pc + cc, place(pf, pc + cc)));
        }
    }
    let mut out = vec![0.0; desired.len()];
    for (first, count, top) in clusters {
        for iii in 0..count {
            out[first + iii] = top + iii as f32 * gap;
        }
    }
    out
}

/// Lay labels out in a right and a left column so none overlap, and draw them with leaders.
#[allow(clippy::too_many_arguments)]
fn place_labels(
    scene: &mut Scene,
    mut labels: Vec<LabelRequest>,
    center: [f32; 2],
    outer_radius: f32,
    font: f32,
    height: f32,
    max_per_side: usize,
    measure: &dyn TextMeasure,
    theme: &super::Theme,
) {
    let line = font * 1.35;
    let column_radius = outer_radius + font * 3.0;
    let knee_radius = outer_radius + font * 1.4;
    // Keep the columns beside the circle: above or below it a column would collapse onto the
    // vertical axis, where its text would run into other labels' leaders.
    let top = (center[1] - column_radius - 2.0 * line).max(font);
    let bottom = (center[1] + column_radius + 2.0 * line).min(height - font);
    let capacity = (((bottom - top) / line).floor().max(1.0) as usize).min(max_per_side.max(1));
    let width = center[0] * 2.0;

    // ii = 0 is the right half (angles 0..pi), ii = 1 the left half.
    for ii in 0..2 {
        let mut side: Vec<LabelRequest> = Vec::new();
        let mut rest = Vec::new();
        for label in labels.drain(..) {
            let a = label.angle.rem_euclid(TAU);
            let on_right = a < std::f32::consts::PI;
            if (ii == 0) == on_right { side.push(label) } else { rest.push(label) }
        }
        labels = rest;
        if side.len() > capacity {
            side.sort_by(|a, b| b.priority.partial_cmp(&a.priority).unwrap_or(std::cmp::Ordering::Equal));
            side.truncate(capacity);
        }

        // Desired heights straight out from each labelled point, sorted top to bottom, then
        // spread so no two overlap while each stays as close as possible to its point.
        let mut placed: Vec<(f32, LabelRequest)> = side.into_iter().map(|l| (center[1] - column_radius * l.angle.cos(), l)).collect();
        placed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let desired: Vec<f32> = placed.iter().map(|(y, _)| *y).collect();
        let spread_ys = spread(&desired, line, top, bottom);
        for (iii, y) in spread_ys.into_iter().enumerate() {
            placed[iii].0 = y;
        }

        for (y, label) in placed {
            // Horizontal position on the label circle at this height; the two rows that may sit
            // just above or below the circle keep a small offset from the axis.
            let dy = y - center[1];
            let dx = (column_radius * column_radius - dy * dy).max(0.0).sqrt().max(font * 3.0);
            let x = if ii == 0 { center[0] + dx } else { center[0] - dx };
            let start = polar(center, label.from_radius, label.angle);
            let knee = polar(center, knee_radius, label.angle);
            let end = [if ii == 0 { x - 2.0 } else { x + 2.0 }, y];
            scene.push(Shape::Polyline { points: vec![start, knee, end], color: theme.faint.mix(theme.muted, 0.35), width: 0.8 }, label.tag.clone());
            let text_x = if ii == 0 { x + 1.0 } else { x - 1.0 };
            let anchor = if ii == 0 { Anchor::Start } else { Anchor::End };
            // Shorten labels that would run off the canvas.
            let room = if ii == 0 { width - text_x - 4.0 } else { text_x - 4.0 };
            let text = fit_text(&label.text, font, label.bold, room.max(font * 3.0), measure);
            scene.push(Shape::Text { pos: [text_x, y], text, size: font, color: label.color, anchor, bold: label.bold }, label.tag);
        }
    }
}
