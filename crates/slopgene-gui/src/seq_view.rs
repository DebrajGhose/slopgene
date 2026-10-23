//! The sequence view: rows of bases with the complementary strand, a ruler, feature bars,
//! translations of coding features and enzyme cuts, with a text-editor-like cursor.
//!
//! Each row holds a multiple of ten bases, with a small gap after every ten so positions are
//! easy to count. Rows are laid out once per document revision and width (feature lanes,
//! translations), and only the rows on screen are painted, so long sequences stay fast.
//!
//! Keys (when the view has focus): letters type bases, Backspace/Delete remove, arrows move
//! (Shift extends the selection), Home/End go to the row's ends, Ctrl+A selects all,
//! Ctrl+C / Ctrl+X / Ctrl+V copy, cut and paste.

use crate::app::ViewSettings;
use crate::doc::DocTab;
use crate::style::{self, Palette};
use eframe::egui::{self, Align2, Color32, CornerRadius, Event, EventFilter, Key, Pos2, Rect, Sense, Stroke, pos2, vec2};
use slopgene_core::render::{feature_rgba, thousands};
use slopgene_core::translate::GeneticCode;
use slopgene_core::*;

const LEFT_MARGIN: f32 = 72.0;
const RIGHT_MARGIN: f32 = 18.0;
const BASE_FONT: f32 = 13.5;
const SMALL_FONT: f32 = 10.0;
const LINE_H: f32 = 18.0;
const ENZYME_LANE_H: f32 = 13.0;
const ENZYME_SUB_LANES: usize = 3;
const RULER_H: f32 = 13.0;
const BAR_H: f32 = 13.0;
const AA_H: f32 = 14.0;
const LANE_GAP: f32 = 4.0;
const ROW_PAD_TOP: f32 = 6.0;
const ROW_PAD_BOTTOM: f32 = 12.0;

/// A piece of a feature inside one row.
#[derive(Clone)]
struct Bar {
    uid: u64,
    /// Base range covered in this row (0-based, half-open).
    start: usize,
    end: usize,
    lane: usize,
    color: Color32,
    label: String,
    /// Arrowhead on the right (forward feature ending here) or left (reverse feature starting here).
    arrow_right: bool,
    arrow_left: bool,
    /// Amino acids of the translation whose middle base falls in this bar: (position, letter).
    amino_acids: Vec<(usize, u8)>,
}

struct Row {
    start: usize,
    end: usize,
    bars: Vec<Bar>,
    lane_heights: Vec<f32>,
    /// Enzyme cuts in this row: (cut position, name, sub-lane).
    enzymes: Vec<(usize, String, usize)>,
    enzyme_lanes: usize,
    height: f32,
}

/// Everything needed to paint and hit-test the view; rebuilt when the key changes.
pub struct Layout {
    key: LayoutKey,
    per_row: usize,
    char_w: f32,
    group_gap: f32,
    rows: Vec<Row>,
    tops: Vec<f32>,
    total_height: f32,
}

impl Layout {
    /// x of the left edge of base `i` of a row (i may equal per_row: the row's right end).
    fn base_x(&self, i: usize) -> f32 {
        LEFT_MARGIN + i as f32 * self.char_w + (i / 10) as f32 * self.group_gap
    }

    /// x of the boundary before base `i` (the cursor position), centred in group gaps.
    fn boundary_x(&self, i: usize) -> f32 {
        if i > 0 && i.is_multiple_of(10) { self.base_x(i) - self.group_gap / 2.0 } else { self.base_x(i) }
    }

    /// The boundary nearest to x within a row (0..=row length).
    fn boundary_at(&self, x: f32, row_len: usize) -> usize {
        let group_w = 10.0 * self.char_w + self.group_gap;
        let rel = (x - LEFT_MARGIN).max(0.0);
        let group = (rel / group_w).floor() as usize;
        let within = ((rel - group as f32 * group_w) / self.char_w).round().clamp(0.0, 10.0) as usize;
        (group * 10 + within).min(row_len)
    }

    fn row_at(&self, y: f32) -> usize {
        match self.tops.binary_search_by(|t| t.partial_cmp(&y).unwrap_or(std::cmp::Ordering::Less)) {
            Ok(ii) => ii,
            Err(ii) => ii.saturating_sub(1),
        }
        .min(self.rows.len().saturating_sub(1))
    }

    fn row_of(&self, pos: usize) -> usize {
        (pos / self.per_row.max(1)).min(self.rows.len().saturating_sub(1))
    }
}

fn strand_block(settings: &ViewSettings) -> f32 {
    LINE_H * if settings.show_complement { 2.0 } else { 1.0 }
}

/// Codon positions and amino acids of a coding feature, as (middle base, letter).
fn translation_marks(feature: &Feature, doc: &Document) -> Vec<(usize, u8)> {
    let len = doc.len();
    // Base positions in the order the feature reads them.
    let mut positions: Vec<usize> = Vec::with_capacity(feature.len());
    for segment in feature.segments.iter().filter(|s| !s.gap) {
        for (a, b) in segment.span.pieces(len) {
            positions.extend(a..b);
        }
    }
    if feature.strand.is_reverse() {
        positions.reverse();
    }
    let codon_start = feature.qualifier_text("codon_start").and_then(|t| t.trim().parse::<usize>().ok()).unwrap_or(1).clamp(1, 3);
    let table =
        feature.qualifier_text("transl_table").and_then(|t| t.trim().parse::<u32>().ok()).and_then(GeneticCode::get).unwrap_or_else(GeneticCode::standard);
    let mut marks = Vec::with_capacity(positions.len() / 3);
    for codon in positions[(codon_start - 1).min(positions.len())..].as_chunks::<3>().0 {
        let mut bases = [doc.seq[codon[0]], doc.seq[codon[1]], doc.seq[codon[2]]];
        if feature.strand.is_reverse() {
            for b in &mut bases {
                *b = alphabet::complement_base(*b);
            }
        }
        marks.push((codon[1], table.translate_codon(&bases)));
    }
    marks
}

/// Lay out rows for the current width and document state.
fn build_layout(ui: &egui::Ui, tab: &mut DocTab, settings: &ViewSettings, width: f32, key: LayoutKey) -> Layout {
    let char_w = ui.ctx().fonts_mut(|f| f.layout_no_wrap("ACGTACGTAC".into(), style::mono_font(BASE_FONT), Color32::WHITE).size().x / 10.0);
    let group_gap = (char_w * 0.6).round().max(4.0);
    let usable = (width - LEFT_MARGIN - RIGHT_MARGIN).max(10.0 * char_w);
    let groups = ((usable + group_gap) / (10.0 * char_w + group_gap)).floor().max(1.0) as usize;
    let per_row = groups * 10;
    let len = tab.len();
    let row_count = len.div_ceil(per_row).max(1);

    let marks = if settings.show_enzymes { tab.sites().0.to_vec() } else { Vec::new() };
    let doc = &tab.doc;
    let mut rows: Vec<Row> = (0..row_count)
        .map(|ii| Row {
            start: ii * per_row,
            end: ((ii + 1) * per_row).min(len),
            bars: Vec::new(),
            lane_heights: Vec::new(),
            enzymes: Vec::new(),
            enzyme_lanes: 0,
            height: 0.0,
        })
        .collect();

    // Feature pieces, cut at row boundaries.
    for feature in doc.features.iter().filter(|f| !f.kind.eq_ignore_ascii_case("source")) {
        let rgba = feature_rgba(feature, settings.colors);
        let color = Color32::from_rgb(rgba.r, rgba.g, rgba.b);
        let extent = feature.span(len);
        let first_base = extent.start % len.max(1);
        let last_base = (extent.end + len - 1) % len.max(1);
        let amino = if settings.show_translations && feature.is_coding() { translation_marks(feature, doc) } else { Vec::new() };
        for segment in feature.segments.iter().filter(|s| !s.gap) {
            for (a, b) in segment.span.pieces(len) {
                let mut pos = a;
                while pos < b {
                    let row = pos / per_row;
                    let row_end = ((row + 1) * per_row).min(b);
                    let aa: Vec<(usize, u8)> = amino.iter().copied().filter(|(p, _)| *p >= pos && *p < row_end).collect();
                    rows[row].bars.push(Bar {
                        uid: feature.uid,
                        start: pos,
                        end: row_end,
                        lane: 0,
                        color,
                        label: feature.name.clone(),
                        arrow_right: feature.strand == Strand::Forward && row_end - 1 == last_base,
                        arrow_left: feature.strand == Strand::Reverse && pos == first_base,
                        amino_acids: aa,
                    });
                    pos = row_end;
                }
            }
        }
    }
    for mark in &marks {
        let row = (mark.pos.saturating_sub(1)) / per_row;
        if row < rows.len() {
            rows[row].enzymes.push((mark.pos, mark.name.clone(), 0));
        }
    }

    // Lanes: greedy packing by start, and two staggered sub-lanes for enzyme names.
    let label_w = |text: &str| ui.ctx().fonts_mut(|f| f.layout_no_wrap(text.to_string(), style::ui_font(SMALL_FONT), Color32::WHITE).size().x);
    let mut tops = Vec::with_capacity(rows.len());
    let mut y = 0.0;
    for row in &mut rows {
        row.bars.sort_by_key(|b| (b.start, std::cmp::Reverse(b.end)));
        let mut lane_ends: Vec<usize> = Vec::new();
        for bar in &mut row.bars {
            let lane = lane_ends.iter().position(|&e| e <= bar.start).unwrap_or_else(|| {
                lane_ends.push(0);
                lane_ends.len() - 1
            });
            lane_ends[lane] = bar.end;
            bar.lane = lane;
        }
        row.lane_heights = (0..lane_ends.len())
            .map(|lane| {
                let has_aa = row.bars.iter().any(|b| b.lane == lane && !b.amino_acids.is_empty());
                BAR_H + if has_aa { AA_H } else { 0.0 }
            })
            .collect();
        // Enzyme names go into up to three staggered rows; a name that fits nowhere is left out
        // (its cut tick is still drawn) rather than drawn over another.
        row.enzymes.sort_by_key(|e| e.0);
        let mut sub_ends = [f32::MIN; ENZYME_SUB_LANES];
        for enzyme in &mut row.enzymes {
            let i = enzyme.0.saturating_sub(row.start);
            let x = LEFT_MARGIN + i as f32 * char_w + (i / 10) as f32 * group_gap;
            let w = label_w(&enzyme.1);
            let lane = (0..ENZYME_SUB_LANES).find(|&l| x - w / 2.0 > sub_ends[l] + 4.0);
            match lane {
                Some(l) => {
                    sub_ends[l] = x + w / 2.0;
                    enzyme.2 = l;
                }
                None => enzyme.2 = usize::MAX,
            }
        }
        row.enzyme_lanes = row.enzymes.iter().filter(|e| e.2 != usize::MAX).map(|e| e.2 + 1).max().unwrap_or(if row.enzymes.is_empty() { 0 } else { 1 });
        let lanes_h: f32 = row.lane_heights.iter().map(|h| h + LANE_GAP).sum();
        row.height = ROW_PAD_TOP + row.enzyme_lanes as f32 * ENZYME_LANE_H + strand_block(settings) + RULER_H + lanes_h + ROW_PAD_BOTTOM;
        tops.push(y);
        y += row.height;
    }
    Layout { key, per_row, char_w, group_gap, rows, tops, total_height: y }
}

/// What the layout depends on: document revision and id, width, display options, enzymes, and
/// the cuts of a cloning plan being made.
type LayoutKey = (u64, u64, u32, bool, bool, bool, String, bool, bool);

fn layout_key(revision: u64, id: u64, width: f32, settings: &ViewSettings, set: &str, unique: bool) -> LayoutKey {
    (
        revision,
        id,
        width.round() as u32,
        settings.show_complement,
        settings.show_translations,
        settings.show_enzymes,
        set.to_string(),
        unique,
        settings.colors == slopgene_core::render::ColorMode::Palette,
    )
}

/// Show the sequence view. The layout is cached in egui's memory per document.
pub fn show(ui: &mut egui::Ui, tab: &mut DocTab, settings: &ViewSettings) {
    let palette = Palette::get(settings.dark);
    let width = ui.available_width();
    let cache_id = egui::Id::unique(("seq-layout", tab.id));
    let key = layout_key(tab.revision, tab.id, width, settings, &tab.enzyme_set, tab.unique_only);
    let cached: Option<std::sync::Arc<Layout>> = ui.ctx().data(|d| d.get_temp(cache_id));
    let layout = match cached {
        Some(l) if l.key == key => l,
        _ => {
            let l = std::sync::Arc::new(build_layout(ui, tab, settings, width, key));
            ui.ctx().data_mut(|d| d.insert_temp(cache_id, l.clone()));
            l
        }
    };

    let view_id = egui::Id::unique(("seq-view", tab.id));
    egui::ScrollArea::vertical().id_salt(("seq-scroll", tab.id)).auto_shrink([false, false]).show_viewport(ui, |ui, viewport| {
        let size = vec2(width, layout.total_height.max(ui.available_height()));
        // A stable id keeps keyboard focus on the view from frame to frame.
        let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
        let response = ui.interact(rect, view_id, Sense::click_and_drag());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, palette.canvas);

        // Scroll a requested position into view.
        if let Some(pos) = tab.scroll_to.take() {
            let row = layout.row_of(pos.min(tab.len().saturating_sub(1)));
            let top = rect.min.y + layout.tops[row];
            ui.scroll_to_rect(Rect::from_min_size(pos2(rect.min.x, top), vec2(10.0, layout.rows[row].height)), Some(egui::Align::Center));
        }

        // Paint the rows on screen.
        let first = layout.row_at(viewport.min.y);
        let last = layout.row_at(viewport.max.y);
        for ii in first..=last.min(layout.rows.len().saturating_sub(1)) {
            paint_row(&painter, &layout, ii, rect.min + vec2(0.0, layout.tops[ii]), tab, settings, &palette);
        }

        handle_pointer(ui, &response, rect, &layout, tab, settings);
        if response.clicked() || response.drag_started() {
            response.request_focus();
        }
        if response.has_focus() {
            ui.memory_mut(|m| m.set_focus_lock_filter(view_id, EventFilter { tab: false, horizontal_arrows: true, vertical_arrows: true, escape: false }));
            handle_keys(ui, &layout, tab);
        }
    });
}

#[allow(clippy::too_many_arguments)]
fn paint_row(painter: &egui::Painter, layout: &Layout, index: usize, origin: Pos2, tab: &DocTab, settings: &ViewSettings, palette: &Palette) {
    let row = &layout.rows[index];
    let seq = &tab.doc.seq;
    let len = seq.len();
    let x = |i: usize| origin.x + layout.base_x(i);
    let mut y = origin.y + ROW_PAD_TOP;

    // Enzyme cuts: name above, tick down to the strand.
    let enzyme_top = y;
    y += row.enzyme_lanes as f32 * ENZYME_LANE_H;
    let strand_top = y;
    for (pos, name, lane) in &row.enzymes {
        let cx = origin.x + layout.boundary_x(pos - row.start);
        if *lane == usize::MAX {
            painter.line_segment([pos2(cx, strand_top - 4.0), pos2(cx, strand_top + LINE_H)], Stroke::new(1.0, palette.faint));
            continue;
        }
        let ty = enzyme_top + (*lane as f32 + 0.5) * ENZYME_LANE_H;
        painter.text(pos2(cx, ty), Align2::CENTER_CENTER, name, style::ui_font(SMALL_FONT), palette.muted);
        painter.line_segment([pos2(cx, ty + 5.0), pos2(cx, strand_top + LINE_H)], Stroke::new(1.0, palette.faint));
    }

    // Selection behind the bases.
    let block = strand_block(settings);
    if let Some(sel) = tab.selection {
        for (a, b) in sel.pieces(len) {
            let s = a.max(row.start);
            let e = b.min(row.end);
            if s < e {
                let r = Rect::from_min_max(pos2(x(s - row.start) - 1.0, strand_top), pos2(x(e - 1 - row.start) + layout.char_w + 1.0, strand_top + block));
                painter.rect_filled(r, CornerRadius::same(2), palette.accent_soft);
            }
        }
    }

    // Position number, then the bases in groups of ten.
    painter.text(
        pos2(origin.x + LEFT_MARGIN - 14.0, strand_top + LINE_H / 2.0),
        Align2::RIGHT_CENTER,
        thousands(row.start + 1),
        style::mono_font(11.0),
        palette.faint,
    );
    let rna = tab.doc.molecule == Molecule::Rna;
    let mut group = row.start;
    while group < row.end {
        let group_end = (group + 10).min(row.end);
        let top: String = seq[group..group_end].iter().map(|&b| if rna && (b == b'T' || b == b't') { (b + 1) as char } else { b as char }).collect();
        painter.text(pos2(x(group - row.start), strand_top + LINE_H / 2.0), Align2::LEFT_CENTER, top, style::mono_font(BASE_FONT), palette.ink);
        if settings.show_complement {
            let bottom: String = seq[group..group_end].iter().map(|&b| alphabet::complement_base(b) as char).collect();
            painter.text(pos2(x(group - row.start), strand_top + LINE_H * 1.5), Align2::LEFT_CENTER, bottom, style::mono_font(BASE_FONT), palette.faint);
        }
        group = group_end;
    }

    // Cursor (only without a selection).
    if tab.selection.is_none() && tab.cursor >= row.start && (tab.cursor < row.end || (tab.cursor == len && row.end == len)) {
        let cx = origin.x + layout.boundary_x(tab.cursor - row.start);
        painter.line_segment([pos2(cx, strand_top - 1.0), pos2(cx, strand_top + block + 1.0)], Stroke::new(1.6, palette.accent));
    }

    // Ruler: a tick after every tenth base, numbered.
    let ruler_y = strand_top + block + 3.0;
    let mut tick = row.start + 10 - (row.start % 10);
    while tick <= row.end {
        let tx = x(tick - 1 - row.start) + layout.char_w / 2.0;
        painter.line_segment([pos2(tx, ruler_y), pos2(tx, ruler_y + 3.0)], Stroke::new(1.0, palette.faint));
        painter.text(pos2(tx, ruler_y + 8.0), Align2::CENTER_CENTER, thousands(tick), style::ui_font(9.0), palette.faint);
        tick += 10;
    }

    // Feature bars, lane by lane, with labels inside when they fit and translations below.
    let mut lane_y = ruler_y + RULER_H;
    let mut lane_tops = Vec::with_capacity(row.lane_heights.len());
    for h in &row.lane_heights {
        lane_tops.push(lane_y);
        lane_y += h + LANE_GAP;
    }
    for bar in &row.bars {
        let top = lane_tops[bar.lane];
        let x0 = x(bar.start - row.start);
        let x1 = x(bar.end - 1 - row.start) + layout.char_w;
        let selected = tab.selected_feature == Some(bar.uid);
        let r = Rect::from_min_max(pos2(x0, top), pos2(x1, top + BAR_H));
        let head = (BAR_H * 0.55).min(r.width() * 0.5);
        let mut points = vec![r.left_top(), r.right_top(), r.right_bottom(), r.left_bottom()];
        if bar.arrow_right {
            points = vec![r.left_top(), pos2(r.right() - head, r.top()), pos2(r.right(), r.center().y), pos2(r.right() - head, r.bottom()), r.left_bottom()];
        }
        if bar.arrow_left {
            points = vec![pos2(r.left() + head, r.top()), r.right_top(), r.right_bottom(), pos2(r.left() + head, r.bottom()), pos2(r.left(), r.center().y)];
        }
        let stroke = if selected { Stroke::new(1.5, palette.accent) } else { Stroke::NONE };
        painter.add(egui::Shape::convex_polygon(points, bar.color, stroke));
        let text_color = if luminance(bar.color) > 0.55 { Color32::from_rgb(0x22, 0x22, 0x22) } else { Color32::WHITE };
        let galley = painter.layout_no_wrap(bar.label.clone(), style::ui_font(SMALL_FONT), text_color);
        if galley.size().x + 8.0 < r.width() {
            painter.galley(pos2(r.center().x - galley.size().x / 2.0, r.center().y - galley.size().y / 2.0), galley, text_color);
        }
        for (pos, aa) in &bar.amino_acids {
            let ax = x(pos - row.start) + layout.char_w / 2.0;
            painter.text(
                pos2(ax, top + BAR_H + AA_H / 2.0),
                Align2::CENTER_CENTER,
                (*aa as char).to_string(),
                style::mono_font(11.0),
                darker(bar.color, palette),
            );
        }
    }
}

fn luminance(c: Color32) -> f32 {
    (0.2126 * c.r() as f32 + 0.7152 * c.g() as f32 + 0.0722 * c.b() as f32) / 255.0
}

/// A colour for text drawn on the canvas in a feature's colour.
fn darker(c: Color32, palette: &Palette) -> Color32 {
    if palette.dark { c } else { Color32::from_rgb((c.r() as f32 * 0.62) as u8, (c.g() as f32 * 0.62) as u8, (c.b() as f32 * 0.62) as u8) }
}

/// Base boundary under a screen point.
fn position_at(layout: &Layout, rect: Rect, p: Pos2) -> usize {
    let row_index = layout.row_at(p.y - rect.min.y);
    let row = &layout.rows[row_index];
    row.start + layout.boundary_at(p.x - rect.min.x, row.end - row.start)
}

/// The feature bar under a screen point, if any.
fn bar_at(layout: &Layout, rect: Rect, p: Pos2, settings_block: f32) -> Option<u64> {
    let row_index = layout.row_at(p.y - rect.min.y);
    let row = &layout.rows[row_index];
    let mut lane_y = rect.min.y + layout.tops[row_index] + ROW_PAD_TOP + row.enzyme_lanes as f32 * ENZYME_LANE_H + settings_block + 3.0 + RULER_H;
    for (lane, h) in row.lane_heights.iter().enumerate() {
        if p.y >= lane_y && p.y <= lane_y + BAR_H {
            for bar in row.bars.iter().filter(|b| b.lane == lane) {
                let x0 = rect.min.x + layout.base_x(bar.start - row.start);
                let x1 = rect.min.x + layout.base_x(bar.end - 1 - row.start) + layout.char_w;
                if p.x >= x0 && p.x <= x1 {
                    return Some(bar.uid);
                }
            }
        }
        lane_y += h + LANE_GAP;
    }
    None
}

fn handle_pointer(ui: &egui::Ui, response: &egui::Response, rect: Rect, layout: &Layout, tab: &mut DocTab, settings: &ViewSettings) {
    let block = strand_block(settings);
    let shift = ui.input(|i| i.modifiers.shift);
    if response.drag_started()
        && let Some(p) = response.interact_pointer_pos()
    {
        let pos = position_at(layout, rect, p);
        if shift {
            let anchor = tab.anchor.unwrap_or(tab.cursor);
            tab.anchor = Some(anchor);
        } else {
            tab.anchor = Some(pos);
            tab.select(None);
        }
        tab.cursor = pos;
    }
    if response.dragged()
        && let (Some(p), Some(anchor)) = (response.interact_pointer_pos(), tab.anchor)
    {
        let pos = position_at(layout, rect, p);
        tab.cursor = pos;
        tab.selection = (pos != anchor).then(|| Span::new(anchor.min(pos), anchor.max(pos)));
        tab.selected_feature = None;
    }
    if response.clicked()
        && let Some(p) = response.interact_pointer_pos()
    {
        if let Some(uid) = bar_at(layout, rect, p, block) {
            tab.select_feature(uid);
            tab.scroll_to = None;
            return;
        }
        let pos = position_at(layout, rect, p);
        if shift {
            let anchor = tab.anchor.unwrap_or(tab.cursor);
            tab.selection = (pos != anchor).then(|| Span::new(anchor.min(pos), anchor.max(pos)));
            tab.cursor = pos;
        } else {
            tab.select(None);
            tab.cursor = pos;
            tab.anchor = Some(pos);
        }
    }
}

/// Typing, deleting and moving the cursor. Copy, cut, paste, select all, undo and redo are the
/// app's (see `app::read_shortcuts`), so they also work from the map and the inspector.
fn handle_keys(ui: &egui::Ui, layout: &Layout, tab: &mut DocTab) {
    let len = tab.len();
    let events = ui.input(|i| i.events.clone());
    let mods = ui.input(|i| i.modifiers);
    let mut moved_to: Option<usize> = None;
    for event in events {
        match event {
            Event::Text(text) => {
                // Spaces and digits are skipped quietly. Other characters that are not sequence
                // codes are refused with a note, so a typo never lands in the sequence.
                let typed: String = text.chars().filter(|c| !c.is_whitespace()).collect();
                if typed.is_empty() {
                    continue;
                }
                let outcome = tab
                    .clean_input(&typed)
                    .map_err(|_| format!("'{typed}' is not a sequence code"))
                    .and_then(|bases| if bases.is_empty() { Ok(()) } else { tab.type_bases(&bases).map_err(|e| e.to_string()) });
                match outcome {
                    Ok(()) => moved_to = Some(tab.cursor),
                    Err(e) => tab.notice = Some((e, std::time::Instant::now())),
                }
            }
            Event::Key { key, pressed: true, modifiers, .. } => {
                let row = layout.row_of(tab.cursor.min(len.saturating_sub(1)));
                let (row_start, row_end) = (layout.rows[row].start, layout.rows[row].end);
                let target = match key {
                    Key::ArrowLeft => Some(tab.cursor.saturating_sub(1)),
                    Key::ArrowRight => Some((tab.cursor + 1).min(len)),
                    Key::ArrowUp => Some(tab.cursor.saturating_sub(layout.per_row)),
                    Key::ArrowDown => Some((tab.cursor + layout.per_row).min(len)),
                    Key::Home if modifiers.command => Some(0),
                    Key::End if modifiers.command => Some(len),
                    Key::Home => Some(row_start),
                    Key::End => Some(row_end),
                    Key::Backspace | Key::Delete => {
                        let before = tab.selection.map(|s| s.len().min(len));
                        let result = if key == Key::Backspace { tab.delete_backward() } else { tab.delete_forward() };
                        match result {
                            // A deleted selection gets a note, because it may be large.
                            Ok(()) => {
                                if let Some(count) = before {
                                    tab.notice = Some((deleted_note(count, tab.unit()), std::time::Instant::now()));
                                }
                            }
                            Err(e) => tab.notice = Some((e.to_string(), std::time::Instant::now())),
                        }
                        moved_to = Some(tab.cursor);
                        None
                    }
                    Key::Escape => {
                        tab.select(None);
                        None
                    }
                    _ => None,
                };
                if let Some(t) = target {
                    if mods.shift {
                        let anchor = tab.anchor.unwrap_or(tab.cursor);
                        tab.anchor = Some(anchor);
                        tab.selection = (t != anchor).then(|| Span::new(anchor.min(t), anchor.max(t)));
                        tab.selected_feature = None;
                    } else {
                        tab.select(None);
                        tab.anchor = Some(t);
                    }
                    tab.cursor = t;
                    moved_to = Some(t);
                }
            }
            _ => {}
        }
    }
    if let Some(pos) = moved_to {
        tab.scroll_to = Some(pos);
    }
}

/// The status-bar note after a selection was deleted, with the way back.
pub fn deleted_note(count: usize, unit: &str) -> String {
    let undo = if cfg!(target_os = "macos") { "⌘Z" } else { "Ctrl+Z" };
    format!("deleted {} {unit} · {undo} to undo", slopgene_core::render::thousands(count))
}
