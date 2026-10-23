//! Linear maps: a horizontal backbone, forward-strand features in lanes above it,
//! reverse-strand features in lanes below, enzyme labels in staggered rows on top and a
//! position ruler at the bottom. Aligned sequences are lines between the features and the
//! ruler: solid where they align, thin where they lack the map's bases, with a mark at each
//! difference.

use super::scene::*;
use super::theme::{difference_color, difference_ink};
use super::{DifferenceKind, MapOptions, SegmentKind, drawable_features, edge_for, feature_rgba, nice_tick_step, thousands};
use crate::model::{Document, Strand};
use crate::span::Span;

/// Height of each aligned sequence's line, and the most rows of difference labels under them.
const TRACK_PITCH: f32 = 16.0;
const DIFFERENCE_LABEL_ROWS: usize = 2;

/// Greedy interval packing: give each `(start, end)` (pixels) the first row where it fits.
fn pack_rows(intervals: &[(f32, f32)], gap: f32) -> Vec<usize> {
    let mut order: Vec<usize> = (0..intervals.len()).collect();
    order.sort_by(|&a, &b| intervals[a].0.partial_cmp(&intervals[b].0).unwrap_or(std::cmp::Ordering::Equal));
    let mut row_ends: Vec<f32> = Vec::new();
    let mut rows = vec![0usize; intervals.len()];
    for ii in order {
        let (start, end) = intervals[ii];
        let row = match row_ends.iter().position(|&e| e + gap <= start) {
            Some(r) => r,
            None => {
                row_ends.push(f32::MIN);
                row_ends.len() - 1
            }
        };
        row_ends[row] = end;
        rows[ii] = row;
    }
    rows
}

pub fn render(doc: &Document, options: &MapOptions, measure: &dyn TextMeasure) -> Scene {
    let theme = &options.theme;
    let font = options.font_size;
    let width = options.width.max(200.0);
    let seq_len = doc.len().max(1);
    let margin = 24.0;
    let scale = (width - 2.0 * margin) / seq_len as f32;
    let x_of = |pos: f64| margin + pos as f32 * scale;
    let thickness = 11.0;
    let lane_pitch = thickness + font * 1.5 + 4.0;
    let line = font * 1.35;

    // Features split by strand; each takes the horizontal room of its arrow or its label,
    // whichever is wider, so labels never collide inside a lane.
    let features = if options.show_features { drawable_features(doc) } else { Vec::new() };
    let mut upper = Vec::new();
    let mut lower = Vec::new();
    for feature in &features {
        let span = feature.span(seq_len);
        let x0 = x_of(span.start as f64);
        let x1 = x_of(span.end.min(seq_len) as f64).max(x0 + 2.0);
        let label_w = measure.text_width(&feature.name, font, false);
        let mid = (x0 + x1) / 2.0;
        let interval = (x0.min(mid - label_w / 2.0), x1.max(mid + label_w / 2.0));
        if feature.strand == Strand::Reverse {
            lower.push((*feature, x0, x1, interval));
        } else {
            upper.push((*feature, x0, x1, interval));
        }
    }
    let upper_rows = pack_rows(&upper.iter().map(|u| u.3).collect::<Vec<_>>(), 6.0);
    let lower_rows = pack_rows(&lower.iter().map(|l| l.3).collect::<Vec<_>>(), 6.0);
    let upper_lanes = upper_rows.iter().copied().max().map(|m| m + 1).unwrap_or(0);
    let lower_lanes = lower_rows.iter().copied().max().map(|m| m + 1).unwrap_or(0);

    // Enzyme labels in staggered rows above everything. A chosen cut is bold.
    let marks = super::all_marks(options);
    let enzyme_labels: Vec<(String, f32)> = marks.iter().map(|(m, _)| (format!("{} ({})", m.name, m.pos.max(1)), x_of(m.pos as f64))).collect();
    let enzyme_intervals: Vec<(f32, f32)> = enzyme_labels
        .iter()
        .zip(&marks)
        .map(|((t, x), (_, chosen))| {
            let w = measure.text_width(t, font, *chosen);
            (x - w / 2.0, x + w / 2.0)
        })
        .collect();
    let enzyme_rows = pack_rows(&enzyme_intervals, 6.0);
    let enzyme_row_count = enzyme_rows.iter().copied().max().map(|m| m + 1).unwrap_or(0);

    // Labels of the emphasised aligned sequence's differences, packed into rows under all the
    // aligned lines. A label that finds no room in the first rows is left out; its mark stays.
    let mut difference_labels: Vec<(usize, usize, String, f32)> = Vec::new();
    for (k, track) in options.alignments.iter().enumerate().filter(|(_, t)| t.emphasis) {
        for (m, mark) in track.marks.iter().enumerate() {
            let x = (x_of(mark.span.start as f64) + x_of((mark.span.start + mark.span.len()) as f64)) / 2.0;
            difference_labels.push((k, m, mark.label.clone(), x));
        }
    }
    let label_intervals: Vec<(f32, f32)> = difference_labels
        .iter()
        .map(|(_, _, text, x)| {
            let w = measure.text_width(text, font * 0.9, false);
            (x - w / 2.0, x + w / 2.0)
        })
        .collect();
    let label_rows = pack_rows(&label_intervals, 6.0);
    let label_row_count = label_rows.iter().copied().filter(|&r| r < DIFFERENCE_LABEL_ROWS).max().map(|m| m + 1).unwrap_or(0);

    // Vertical layout, top to bottom.
    let title_h = font * 2.4;
    let enzymes_top = title_h + font * 0.5;
    let enzymes_h = enzyme_row_count as f32 * line;
    let upper_h = upper_lanes as f32 * lane_pitch;
    let backbone_y = enzymes_top + enzymes_h + upper_h + thickness + 6.0;
    let lower_h = lower_lanes as f32 * lane_pitch;
    let tracks_top = backbone_y + thickness / 2.0 + lower_h + 8.0;
    let tracks_h = if options.alignments.is_empty() { 0.0 } else { 8.0 + options.alignments.len() as f32 * TRACK_PITCH + label_row_count as f32 * line };
    let ruler_y = backbone_y + thickness / 2.0 + lower_h + tracks_h + font * 1.2;
    let height = (ruler_y + font * 2.6).max(options.height.min(ruler_y + font * 2.6));
    let mut scene = Scene::new(width, height, theme.background);

    // Title.
    scene.push(
        Shape::Text { pos: [margin, title_h / 2.0 + 2.0], text: doc.name.clone(), size: font * 1.35, color: theme.ink, anchor: Anchor::Start, bold: true },
        Tag::None,
    );
    let title_w = measure.text_width(&doc.name, font * 1.35, true);
    scene.push(
        Shape::Text {
            pos: [margin + title_w + font, title_h / 2.0 + 2.0],
            text: format!("{} bp, linear", thousands(doc.len())),
            size: font,
            color: theme.muted,
            anchor: Anchor::Start,
            bold: false,
        },
        Tag::None,
    );

    // Selection band behind everything.
    if let Some(sel) = options.selection.filter(|s| !s.is_empty()) {
        let x0 = x_of(sel.start as f64);
        let x1 = x_of(sel.end.min(seq_len) as f64);
        let top = backbone_y - upper_h - thickness;
        let bottom = backbone_y + thickness / 2.0 + lower_h;
        scene.push(Shape::Polygon { points: vec![[x0, top], [x1, top], [x1, bottom], [x0, bottom]], fill: theme.selection, stroke: None }, Tag::Selection);
    }

    // Backbone.
    scene.push(Shape::Polyline { points: vec![[x_of(0.0), backbone_y], [x_of(seq_len as f64), backbone_y]], color: theme.backbone, width: 2.0 }, Tag::Backbone);

    // Features: arrows in lanes, labels just outside each arrow.
    let mut draw_lane_feature = |feature: &crate::model::Feature, x0: f32, x1: f32, lane: usize, above: bool| {
        let y = if above {
            backbone_y - thickness / 2.0 - 4.0 - lane as f32 * lane_pitch - thickness / 2.0
        } else {
            backbone_y + thickness / 2.0 + 4.0 + lane as f32 * lane_pitch + thickness / 2.0
        };
        let color = feature_rgba(feature, options.colors);
        let edge = edge_for(color, theme);
        let tag = Tag::Feature(feature.uid);
        if options.highlight == Some(feature.uid) {
            scene.push(
                Shape::Polygon {
                    points: vec![
                        [x0 - 2.5, y - thickness / 2.0 - 2.5],
                        [x1 + 2.5, y - thickness / 2.0 - 2.5],
                        [x1 + 2.5, y + thickness / 2.0 + 2.5],
                        [x0 - 2.5, y + thickness / 2.0 + 2.5],
                    ],
                    fill: theme.accent,
                    stroke: None,
                },
                tag.clone(),
            );
        }
        let head = (thickness * 0.9).min((x1 - x0) * 0.6);
        let (t, b) = (y - thickness / 2.0, y + thickness / 2.0);
        let points = match feature.strand {
            Strand::Forward => vec![[x0, t], [x1 - head, t], [x1, y], [x1 - head, b], [x0, b]],
            Strand::Reverse => vec![[x0 + head, t], [x1, t], [x1, b], [x0 + head, b], [x0, y]],
            Strand::Both => vec![[x0 + head, t], [x1 - head, t], [x1, y], [x1 - head, b], [x0 + head, b], [x0, y]],
            Strand::None => vec![[x0, t], [x1, t], [x1, b], [x0, b]],
        };
        scene.push(Shape::Polygon { points, fill: color, stroke: edge.map(|e| (e, 1.0)) }, tag.clone());
        let label_y = if above { t - font * 0.75 } else { b + font * 0.8 };
        scene.push(
            Shape::Text { pos: [(x0 + x1) / 2.0, label_y], text: feature.name.clone(), size: font, color: theme.ink, anchor: Anchor::Middle, bold: false },
            tag,
        );
    };
    for (ii, (feature, x0, x1, _)) in upper.iter().enumerate() {
        draw_lane_feature(feature, *x0, *x1, upper_rows[ii], true);
    }
    for (ii, (feature, x0, x1, _)) in lower.iter().enumerate() {
        draw_lane_feature(feature, *x0, *x1, lower_rows[ii], false);
    }

    // Primers as thin arrows hugging the backbone.
    if options.show_primers {
        for primer in &doc.primers {
            for site in primer.sites.iter().filter(|s| s.extra.attr("simplified") != Some("1")) {
                let x0 = x_of(site.span.start as f64);
                let x1 = x_of(site.span.end.min(seq_len) as f64);
                let y = if site.strand == Strand::Reverse { backbone_y + 4.0 } else { backbone_y - 4.0 };
                scene.push(Shape::Polyline { points: vec![[x0, y], [x1, y]], color: theme.primer, width: 2.5 }, Tag::Primer(primer.uid));
            }
        }
    }

    // Enzyme ticks with leaders up to their label rows. A chosen cut is drawn in the accent colour
    // and longer.
    for (ii, (mark, chosen)) in marks.iter().enumerate() {
        let x = enzyme_labels[ii].1;
        let label_y = enzymes_top + enzyme_rows[ii] as f32 * line + line / 2.0;
        let tag = Tag::Enzyme(mark.name.clone(), mark.pos);
        let (color, width, reach) = if *chosen { (theme.accent, 2.4, 8.0) } else { (theme.enzyme, 1.2, 5.0) };
        scene.push(Shape::Polyline { points: vec![[x, label_y + font * 0.6], [x, backbone_y + 4.0]], color: theme.faint, width: 0.8 }, tag.clone());
        scene.push(Shape::Polyline { points: vec![[x, backbone_y - reach], [x, backbone_y + reach]], color, width }, tag.clone());
        scene.push(Shape::Text { pos: [x, label_y], text: enzyme_labels[ii].0.clone(), size: font, color, anchor: Anchor::Middle, bold: *chosen }, tag);
    }

    // Aligned sequences: a line each under the features, drawn as the rings of a circular map
    // are. The line is a solid green band where the sequence has the map's bases, an orange band
    // where it has them turned round, and thin and red where it lacks them. There is no line
    // where the sequence gives no information. A change of a few bases is a short bar across the
    // band, and extra bases hang from it as a flag. Tracks other than the emphasised one are
    // faded.
    for (k, track) in options.alignments.iter().enumerate() {
        let y = tracks_top + k as f32 * TRACK_PITCH + TRACK_PITCH / 2.0;
        let fade = |c: Rgba| if track.emphasis { c } else { c.mix(theme.background, 0.55) };
        let band = if track.emphasis { 6.0 } else { 3.5 };
        let thin = 1.5;
        let ends = |span: Span| (x_of(span.start as f64), x_of((span.start + span.len()).min(seq_len) as f64));
        for segment in &track.segments {
            let (color, width) = match segment.kind {
                SegmentKind::Aligned => (theme.aligned, band),
                SegmentKind::Inverted => (theme.changed, band),
                SegmentKind::Missing => (theme.missing, thin),
            };
            let (x0, x1) = ends(segment.span);
            scene.push(Shape::Polyline { points: vec![[x0, y], [x1, y]], color: fade(color), width }, Tag::Track(k));
        }
        for (m, mark) in track.marks.iter().enumerate() {
            let color = fade(difference_color(mark.kind, theme));
            let (x0, x1) = ends(mark.span);
            let x = (x0 + x1) / 2.0;
            let tag = Tag::Difference(k, m);
            if mark.structural && !mark.span.is_empty() {
                // The segment draws the stretch; the same line again carries the change's tag.
                let width = if mark.kind == DifferenceKind::Missing { thin } else { band };
                scene.push(Shape::Polyline { points: vec![[x0, y], [x1, y]], color, width }, tag);
            } else if mark.structural {
                // Extra sequence: a flag that hangs from the line, with a dot at its tip.
                scene.push(Shape::Polyline { points: vec![[x, y - band / 2.0], [x, y + 8.0]], color, width: 2.5 }, tag.clone());
                scene.push(Shape::Circle { center: [x, y + 10.0], radius: 3.0, fill: Some(color), stroke: None }, tag);
            } else {
                let reach = band / 2.0 + 3.5;
                scene.push(Shape::Polyline { points: vec![[x, y - reach], [x, y + reach]], color, width: 2.5 }, tag);
            }
        }
    }
    let labels_top = tracks_top + options.alignments.len() as f32 * TRACK_PITCH + 4.0;
    for (ii, (k, m, text, x)) in difference_labels.into_iter().enumerate() {
        if label_rows[ii] >= DIFFERENCE_LABEL_ROWS {
            continue;
        }
        let mark = &options.alignments[k].marks[m];
        let y = labels_top + label_rows[ii] as f32 * line + line / 2.0;
        let color = difference_ink(mark.kind, theme);
        scene.push(Shape::Text { pos: [x, y], text, size: font * 0.9, color, anchor: Anchor::Middle, bold: mark.structural }, Tag::Difference(k, m));
    }

    // Ruler.
    if options.show_ticks {
        let step = nice_tick_step(seq_len, 10);
        scene.push(Shape::Polyline { points: vec![[x_of(0.0), ruler_y], [x_of(seq_len as f64), ruler_y]], color: theme.faint, width: 1.0 }, Tag::None);
        let mut pos = 0;
        while pos <= seq_len {
            let x = x_of(pos as f64);
            scene.push(Shape::Polyline { points: vec![[x, ruler_y], [x, ruler_y + 4.0]], color: theme.faint, width: 1.0 }, Tag::None);
            scene.push(
                Shape::Text {
                    pos: [x, ruler_y + font * 1.2],
                    text: thousands(pos.max(1)),
                    size: font * 0.85,
                    color: theme.muted,
                    anchor: Anchor::Middle,
                    bold: false,
                },
                Tag::None,
            );
            pos += step;
        }
    }

    if let Some(cursor) = options.cursor {
        let x = x_of(cursor as f64);
        scene.push(
            Shape::Polyline { points: vec![[x, backbone_y - upper_h - thickness], [x, backbone_y + thickness + lower_h]], color: theme.accent, width: 1.5 },
            Tag::None,
        );
    }
    scene
}
