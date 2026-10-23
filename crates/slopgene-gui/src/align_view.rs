//! The alignment view: sequences aligned to the open document, with the plasmid map and the
//! aligned bases side by side, as in SnapGene's "Align to Reference".
//!
//! The map on the left draws the document as its circle, and each aligned sequence as a ring
//! outside it (a line under a linear map). The ring shows what the sequence has of the map: a
//! solid green band where the two align, a thin red line where the sequence lacks the map's bases,
//! and nothing where it gives no information (past the ends of a partial read). Marks on the ring
//! show differences of a few bases, and flags show bases only the sequence has.
//!
//! The pane on the right lays the alignment out in rows: the document's features over its bases,
//! the aligned sequence's bases below, a coloured band behind every difference and a short label
//! saying what changed. Each structural change stands in one row in place of its bases. Over the
//! rows an overview bar draws the whole map coloured as the ring is, with the part on screen
//! outlined, so the rows and the map can be read together. The inspector beside them lists the
//! aligned sequences, how much of the map aligns, the colour key, the structural changes, the
//! differences and the alignment method.
//!
//! Colours (`slopgene_core::render::theme`): green for aligned bases, red for bases the aligned
//! sequence lacks, orange for different bases and violet for bases only it has. Labels and
//! shapes say the same, so colour is never the only cue.

use crate::app::ViewSettings;
use crate::doc::DocTab;
use crate::map_view;
use crate::panels::{section_title, swatch};
use crate::scene_paint::color;
use crate::style::{self, Palette};
use eframe::egui::text::{LayoutJob, TextFormat};
use eframe::egui::{self, Align2, Color32, CornerRadius, FontId, Pos2, Rect, RichText, Sense, Stroke, StrokeKind, pos2, vec2};
use slopgene_core::align::{AlignedColumn, Alignment, AlignmentBlock, EventKind, StructuralEvent, align_read, method_line, method_summary};
use slopgene_core::model::Feature;
use slopgene_core::render::theme::{difference_color, difference_ink};
use slopgene_core::render::{AlignmentTrack, ColorMode, DifferenceKind, DifferenceMark, SegmentKind, Theme, feature_rgba, thousands};
use slopgene_core::span::Span;
use slopgene_core::{Document, Strand};
use std::ops::Range;
use std::sync::Arc;

/// Font size of the aligned bases, and the height of each sequence's line.
const BASE_FONT: f32 = 13.0;
const LINE_H: f32 = 18.0;
/// The line between the two sequences, which carries a tick for each base that agrees.
const MARK_H: f32 = 12.0;
/// Feature bars over the bases.
const BAR_H: f32 = 13.0;
const LANE_GAP: f32 = 3.0;
/// Difference labels under the bases.
const LABEL_H: f32 = 15.0;
const LABEL_FONT: f32 = 10.5;
const ROW_PAD_TOP: f32 = 8.0;
const ROW_PAD_BOTTOM: f32 = 10.0;
/// The row of a structural change: a headline and a line of detail.
const EVENT_H: f32 = 60.0;
/// The overview of the whole map over the rows.
const OVERVIEW_H: f32 = 30.0;

/// One sequence aligned to the document.
pub struct AlignedEntry {
    pub name: String,
    /// The sequence and whether it is a whole circular molecule, kept to align it again after the
    /// document changes.
    bases: Vec<u8>,
    circular: bool,
    /// The alignment, or `None` when the two sequences share no 15-base stretch.
    pub result: Option<Alignment>,
}

/// The sequences aligned to one document, and what the alignment view shows of them.
#[derive(Default)]
pub struct Alignments {
    pub entries: Vec<AlignedEntry>,
    /// The entry shown in detail.
    pub current: usize,
    /// The document revision the entries were aligned against. Positions go stale when the
    /// document changes, so a new revision aligns them all again.
    revision: u64,
    /// Bumped whenever an entry is added, removed or aligned again, so a cached layout knows.
    generation: u64,
    /// A reference position (0-based) to bring into view among the aligned bases, and for how
    /// many more frames. The request is repeated for a few frames because the rows are laid out
    /// again when the panels beside them settle, which moves the position.
    reveal: Option<usize>,
    reveal_frames: u8,
}

/// Frames for which a request to bring a position into view is repeated.
const REVEAL_FRAMES: u8 = 6;

impl Alignments {
    /// Align a sequence to the document (at `revision`) and show it.
    pub fn add(&mut self, doc: &Document, revision: u64, name: String, bases: Vec<u8>, circular: bool) {
        // Entries aligned before an edit are aligned again first, so all of them fit the document.
        self.refresh(doc, revision);
        let result = align_read(doc, &bases, circular);
        self.entries.push(AlignedEntry { name, bases, circular, result });
        self.current = self.entries.len() - 1;
        self.generation += 1;
        self.reveal_frames = 0;
    }

    /// Align every sequence again if the document changed since they were aligned.
    pub fn refresh(&mut self, doc: &Document, revision: u64) {
        if revision == self.revision {
            return;
        }
        self.revision = revision;
        for entry in &mut self.entries {
            entry.result = align_read(doc, &entry.bases, entry.circular);
        }
        self.generation += 1;
    }

    /// Take entry `ii` off the list. The entry shown stays the same one, unless it was the one
    /// removed.
    pub fn remove(&mut self, ii: usize) {
        if ii >= self.entries.len() {
            return;
        }
        self.entries.remove(ii);
        if ii < self.current {
            self.current -= 1;
        }
        self.current = self.current.min(self.entries.len().saturating_sub(1));
        self.generation += 1;
    }

    /// Show entry `ii` in detail.
    pub fn choose(&mut self, ii: usize) {
        if ii < self.entries.len() && ii != self.current {
            self.current = ii;
            self.reveal_frames = 0;
        }
    }

    /// Bring a reference position (0-based) into view among the aligned bases.
    pub fn reveal(&mut self, position: usize) {
        self.reveal = Some(position);
        self.reveal_frames = REVEAL_FRAMES;
    }

    /// The position to bring into view this frame, if a request is still being carried out.
    fn take_reveal(&mut self) -> Option<usize> {
        if self.reveal_frames == 0 {
            return None;
        }
        self.reveal_frames -= 1;
        self.reveal
    }

    /// One map track per entry, in list order, with the entry shown in detail emphasised. An
    /// entry that does not align gets an empty track, so track numbers match entry numbers.
    pub fn tracks(&self) -> Vec<AlignmentTrack> {
        self.entries
            .iter()
            .enumerate()
            .map(|(ii, entry)| match &entry.result {
                Some(aln) => aln.track(&entry.name, ii == self.current),
                None => AlignmentTrack::empty(&entry.name, ii == self.current),
            })
            .collect()
    }

    /// Mark `mark` of track `track`, as [`Alignments::tracks`] draws it.
    pub fn mark(&self, track: usize, mark: usize) -> Option<DifferenceMark> {
        self.entries.get(track)?.result.as_ref()?.track("", false).marks.get(mark).cloned()
    }
}

/// What the alignment view asks of the app after a frame.
#[derive(Debug, Default, PartialEq)]
pub struct Requests {
    /// Ask for files to align.
    pub pick_files: bool,
    /// Align this open document, by tab id.
    pub add_tab: Option<u64>,
}

impl Requests {
    /// Both sets of requests together.
    pub fn and(self, other: Requests) -> Requests {
        Requests { pick_files: self.pick_files || other.pick_files, add_tab: self.add_tab.or(other.add_tab) }
    }
}

/// The colour of a kind of difference, from the map theme, for this appearance.
fn kind_color(kind: DifferenceKind, dark: bool) -> Color32 {
    color(difference_color(kind, &if dark { Theme::dark() } else { Theme::light() }))
}

/// The same colour for text.
fn kind_ink(kind: DifferenceKind, dark: bool) -> Color32 {
    color(difference_ink(kind, &if dark { Theme::dark() } else { Theme::light() }))
}

/// The colour of aligned bases, for this appearance.
fn aligned_color(dark: bool) -> Color32 {
    color(if dark { Theme::dark() } else { Theme::light() }.aligned)
}

/// The colour key, in the aligned sequence's name.
fn key(read: &str, dark: bool) -> [(Color32, String); 4] {
    [
        (aligned_color(dark), "aligned".to_string()),
        (kind_color(DifferenceKind::Missing, dark), format!("missing from {read}")),
        (kind_color(DifferenceKind::Changed, dark), "different base".to_string()),
        (kind_color(DifferenceKind::Inserted, dark), format!("only in {read}")),
    ]
}

/// The alignment view: the map with the aligned sequences on the left, and the aligned bases on
/// the right. `others` lists the other open DNA documents (tab id and title) that could be added.
pub fn show(ui: &mut egui::Ui, tab: &mut DocTab, settings: &ViewSettings, others: &[(u64, String)]) -> Requests {
    let palette = Palette::get(settings.dark);
    let mut requests = Requests::default();
    tab.alignments.refresh(&tab.doc, tab.revision);
    let width = ui.available_width();
    egui::Panel::left(("align-map", tab.id)).resizable(true).default_size(width * 0.4).min_size(240.0).frame(egui::Frame::new().fill(palette.canvas)).show(
        ui,
        |ui| {
            map_view::show(ui, tab, settings);
        },
    );
    if tab.alignments.entries.is_empty() {
        empty_state(ui, &tab.doc.name, others, &palette, &mut requests);
    } else {
        aligned_pane(ui, tab, settings, &palette);
    }
    requests
}

/// The buttons that add sequences: files, or another open document.
fn add_buttons(ui: &mut egui::Ui, others: &[(u64, String)], requests: &mut Requests) {
    if ui.button("Add files…").on_hover_text("Sequencing results or other sequences (GenBank, FASTA, SnapGene)").clicked() {
        requests.pick_files = true;
    }
    ui.add_enabled_ui(!others.is_empty(), |ui| {
        ui.menu_button("Add an open sequence", |ui| {
            for (id, title) in others {
                if ui.button(title).clicked() {
                    requests.add_tab = Some(*id);
                    ui.close();
                }
            }
        })
        .response
        .on_disabled_hover_text("Open another sequence to compare it with this one");
    });
}

/// What the view shows before anything is aligned: what it does, and the buttons to start.
fn empty_state(ui: &mut egui::Ui, name: &str, others: &[(u64, String)], palette: &Palette, requests: &mut Requests) {
    ui.vertical_centered(|ui| {
        ui.add_space((ui.available_height() * 0.3).max(20.0));
        ui.label(RichText::new(format!("Compare sequences with {name}")).font(style::bold_font(17.0)).color(palette.ink));
        ui.add_space(4.0);
        ui.label(RichText::new("Add sequencing results or other sequences. Each one is aligned to this map,").color(palette.muted));
        ui.label(RichText::new("and every difference is marked here and on the map.").color(palette.muted));
        ui.add_space(12.0);
        add_buttons(ui, others, requests);
    });
}

/// The pane of aligned bases for the entry shown, with a caption saying what is compared and how.
fn aligned_pane(ui: &mut egui::Ui, tab: &mut DocTab, settings: &ViewSettings, palette: &Palette) {
    // A position to bring into view: a click on the map, or a difference chosen in the panel.
    if let Some(position) = tab.scroll_to.take() {
        tab.alignments.reveal(position);
    }
    let reveal = tab.alignments.take_reveal();
    let (current, generation) = (tab.alignments.current, tab.alignments.generation);
    let Some(entry) = tab.alignments.entries.get(current) else { return };
    caption(ui, &tab.doc.name, entry, settings.dark, palette);
    let Some(aln) = &entry.result else {
        let text = format!("{} shares no 15-base stretch with {}, so it does not align.", entry.name, tab.doc.name);
        ui.label(RichText::new(text).color(palette.warn));
        return;
    };
    // The overview's place, kept now and painted once the rows show which part is on screen.
    let (overview_rect, overview) = ui.allocate_exact_size(vec2(ui.available_width(), OVERVIEW_H), Sense::click());
    let track = aln.track(&entry.name, true);
    ui.horizontal_wrapped(|ui| {
        ui.add_space(8.0);
        ui.label(RichText::new(format!("Method: {}. More in the panel on the right.", method_line())).size(11.0).color(palette.faint));
    });

    // The layout is cached per document, rebuilt when the entry, its alignment, the width or the
    // colour mode change.
    let width = ui.available_width();
    let key: LayoutKey = (tab.id, generation, current, width.round() as u32, settings.colors == ColorMode::Palette);
    let cache_id = egui::Id::unique(("align-layout", tab.id));
    let cached: Option<Arc<Layout>> = ui.ctx().data(|d| d.get_temp(cache_id));
    let layout = match cached {
        Some(layout) if layout.key == key => layout,
        _ => {
            let layout = Arc::new(build_layout(ui, aln, &tab.doc, width, settings.colors, key));
            ui.ctx().data_mut(|d| d.insert_temp(cache_id, layout.clone()));
            layout
        }
    };

    let len = tab.len();
    let selection = tab.selection;
    let mut clicked = None;
    let mut on_screen = 0..0;
    egui::ScrollArea::vertical().id_salt(("align-rows", tab.id)).auto_shrink([false, false]).show_viewport(ui, |ui, viewport| {
        let size = vec2(width, layout.total_height.max(ui.available_height()));
        let (rect, response) = ui.allocate_exact_size(size, Sense::click());
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, palette.canvas);
        if let Some(position) = reveal
            && let Some(ii) = aln.block_at(&layout.blocks, position)
        {
            let row = Rect::from_min_size(pos2(rect.min.x, rect.min.y + layout.tops[ii]), vec2(10.0, layout.rows[ii].height));
            ui.scroll_to_rect(row, Some(egui::Align::Center));
        }
        // Only the rows on screen are painted.
        let first = layout.row_at(viewport.min.y);
        let last = layout.row_at(viewport.max.y).min(layout.rows.len().saturating_sub(1));
        on_screen = first..last + 1;
        for ii in first..=last {
            let origin = rect.min + vec2(0.0, layout.tops[ii]);
            match &layout.blocks[ii] {
                AlignmentBlock::Bases(columns) => paint_bases(&painter, &layout, ii, columns.clone(), origin, aln, selection, settings.dark, palette),
                AlignmentBlock::Event(e) => paint_event(&painter, &layout, origin, rect.width(), aln, *e, (&tab.doc.name, &entry.name), settings.dark, palette),
            }
        }
        let response = match response.hover_pos().and_then(|p| describe_at(&layout, aln, rect, p)) {
            Some(text) => response.on_hover_text_at_pointer(text),
            None => response,
        };
        if response.clicked()
            && let Some(p) = response.interact_pointer_pos()
        {
            clicked = click_at(&layout, aln, rect, p, len);
        }
    });
    // The overview, with the reference stretch of the rows on screen outlined. A click on it
    // brings that part of the map into view among the rows.
    let shown: Vec<Span> = on_screen.filter_map(|ii| row_span(aln, &layout.blocks[ii])).collect();
    paint_overview(&ui.painter_at(overview_rect), overview_rect, &track, len, &shown, settings.dark, palette);
    let pointed = |p: Pos2| overview_position(overview_rect, p, len);
    if let Some(position) = overview.hover_pos().map(pointed) {
        let _ = overview.clone().on_hover_text_at_pointer(describe_position(&track, position, len, &tab.doc.name));
    }
    if overview.clicked()
        && let Some(position) = overview.interact_pointer_pos().map(pointed)
    {
        tab.alignments.reveal(position);
    }
    match clicked {
        Some(Click::Bases(span)) => tab.select(Some(span)),
        Some(Click::Feature(uid)) => tab.select_feature(uid),
        None => {}
    }
}

/// Room on each side of the overview bar.
const OVERVIEW_MARGIN: f32 = 12.0;

/// The x of a reference position (0-based, or the length for the end) on the overview.
fn overview_x(rect: Rect, position: usize, len: usize) -> f32 {
    let usable = rect.width() - 2.0 * OVERVIEW_MARGIN;
    rect.left() + OVERVIEW_MARGIN + usable * position as f32 / len.max(1) as f32
}

/// The reference position (0-based) under a point of the overview.
fn overview_position(rect: Rect, p: Pos2, len: usize) -> usize {
    let usable = (rect.width() - 2.0 * OVERVIEW_MARGIN).max(1.0);
    let share = ((p.x - rect.left() - OVERVIEW_MARGIN) / usable).clamp(0.0, 1.0);
    ((share * len as f32) as usize).min(len.saturating_sub(1))
}

/// The x ranges a reference stretch takes on the overview: one, or two for a stretch through
/// the origin.
fn overview_pieces(rect: Rect, span: Span, len: usize) -> Vec<(f32, f32)> {
    span.pieces(len).into_iter().map(|(a, b)| (overview_x(rect, a, len), overview_x(rect, b, len))).collect()
}

/// The reference stretch a row shows: its bases from first to last, or the stretch of a
/// structural change. `None` for a row with no reference bases.
fn row_span(aln: &Alignment, block: &AlignmentBlock) -> Option<Span> {
    let len = aln.reference_len.max(1);
    match block {
        AlignmentBlock::Bases(columns) => {
            let (first, last) = aln.ends(columns.clone()).reference?;
            Some(Span::new(first, first + (last + len - first) % len + 1).normalize(len))
        }
        AlignmentBlock::Event(e) => {
            let span = aln.events[*e].reference;
            Some(if span.is_empty() { Span::new(span.start, span.start + 1).normalize(len) } else { span })
        }
    }
}

/// The whole map as one bar, coloured as the ring on the map: a green band where the sequence
/// aligns, a thin red line where it lacks the map's bases, and nothing where it gives no
/// information. Marks show differences of a few bases and flags extra bases. The stretches in
/// `shown` (the rows on screen) are outlined.
fn paint_overview(painter: &egui::Painter, rect: Rect, track: &AlignmentTrack, len: usize, shown: &[Span], dark: bool, palette: &Palette) {
    let y = rect.top() + 11.0;
    // The whole map as a faint line, so the bar reads as all of it.
    painter.line_segment([pos2(overview_x(rect, 0, len), y), pos2(overview_x(rect, len, len), y)], Stroke::new(1.0, palette.line));
    for segment in &track.segments {
        let (colour, height) = match segment.kind {
            SegmentKind::Aligned => (aligned_color(dark), 8.0),
            SegmentKind::Inverted => (kind_color(DifferenceKind::Changed, dark), 8.0),
            SegmentKind::Missing => (kind_color(DifferenceKind::Missing, dark), 2.5),
        };
        for (x0, x1) in overview_pieces(rect, segment.span, len) {
            painter.rect_filled(Rect::from_min_max(pos2(x0, y - height / 2.0), pos2(x1.max(x0 + 1.0), y + height / 2.0)), 1.0, colour);
        }
    }
    for mark in &track.marks {
        let colour = kind_color(mark.kind, dark);
        let x = overview_x(rect, mark.span.start, len);
        if mark.structural && mark.span.is_empty() {
            // Extra bases: a flag above the bar.
            painter.line_segment([pos2(x, y - 4.0), pos2(x, y - 9.0)], Stroke::new(2.0, colour));
            painter.circle_filled(pos2(x, y - 10.0), 2.5, colour);
        } else if !mark.structural {
            painter.line_segment([pos2(x, y - 6.0), pos2(x, y + 6.0)], Stroke::new(1.5, colour));
        }
    }
    // The part of the map the rows show.
    let mut outlined: Vec<(f32, f32)> = shown.iter().flat_map(|span| overview_pieces(rect, *span, len)).collect();
    outlined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut joined: Vec<(f32, f32)> = Vec::new();
    for (x0, x1) in outlined {
        match joined.last_mut() {
            Some(last) if x0 <= last.1 + 2.0 => last.1 = last.1.max(x1),
            _ => joined.push((x0, x1)),
        }
    }
    for (x0, x1) in joined {
        let frame = Rect::from_min_max(pos2(x0 - 1.0, y - 8.0), pos2(x1.max(x0 + 2.0) + 1.0, y + 8.0));
        painter.rect_stroke(frame, 2.0, Stroke::new(1.5, palette.accent), StrokeKind::Outside);
    }
    // The ends of the map, and what the outline means.
    let font = style::ui_font(10.0);
    let label_y = rect.top() + 24.0;
    painter.text(pos2(overview_x(rect, 0, len), label_y), Align2::LEFT_CENTER, "1", font.clone(), palette.faint);
    painter.text(pos2(overview_x(rect, len, len), label_y), Align2::RIGHT_CENTER, thousands(len), font.clone(), palette.faint);
    painter.text(
        pos2(rect.center().x, label_y),
        Align2::CENTER_CENTER,
        "the whole map; the outline marks the rows below (click to go there)",
        font,
        palette.faint,
    );
}

/// What the sequence has at a reference position, for the overview's tooltip.
fn describe_position(track: &AlignmentTrack, position: usize, len: usize, reference: &str) -> String {
    let base = format!("{reference} base {}", thousands(position + 1));
    let inside = |span: Span| span.contains(position, len.max(1));
    match track.segments.iter().find(|s| inside(s.span)) {
        Some(segment) => {
            let what = match segment.kind {
                SegmentKind::Aligned => format!("aligned in {}", track.name),
                SegmentKind::Inverted => format!("turned round in {}", track.name),
                SegmentKind::Missing => format!("missing from {}", track.name),
            };
            format!("{base}: {what} ({})", segment.span.format(len))
        }
        None => format!("{base}: {} does not reach here", track.name),
    }
}

/// Lines over the aligned bases: which sequence is the map's circle and the top row, which is
/// the ring and the bottom row, and the colour key.
fn caption(ui: &mut egui::Ui, name: &str, entry: &AlignedEntry, dark: bool, palette: &Palette) {
    let small = |text: String, color: Color32| RichText::new(text).size(11.5).color(color);
    let strand = if entry.result.as_ref().is_some_and(|a| a.reverse) { ", read as its reverse complement" } else { "" };
    ui.add_space(4.0);
    ui.horizontal_wrapped(|ui| {
        ui.add_space(8.0);
        ui.label(small(format!("{name}: the circle of the map, and the top row of bases."), palette.ink));
    });
    ui.horizontal_wrapped(|ui| {
        ui.add_space(8.0);
        ui.label(small(format!("{}: the ring outside the map, and the bottom row{strand}.", entry.name), palette.ink));
    });
    ui.horizontal_wrapped(|ui| {
        ui.add_space(8.0);
        for (colour, words) in key(&entry.name, dark) {
            swatch(ui, colour);
            ui.label(small(words, palette.muted));
            ui.add_space(6.0);
        }
    });
}

/// A feature over a run of columns in one row.
struct RowBar {
    uid: u64,
    /// The columns covered, counted from the row's first column, end exclusive.
    start: usize,
    end: usize,
    lane: usize,
    color: Color32,
    label: String,
    arrow_right: bool,
    arrow_left: bool,
}

/// One row of the layout. A row of bases has feature bars over it and difference labels under
/// it. A structural change's row has neither.
struct Row {
    bars: Vec<RowBar>,
    lanes: usize,
    /// Difference labels: the column (from the row's first), the text and the kind.
    labels: Vec<(usize, String, DifferenceKind)>,
    height: f32,
}

/// What the layout depends on: the tab, the alignment generation, the entry shown, the width and
/// whether feature colours come from the palette.
type LayoutKey = (u64, u64, usize, u32, bool);

/// The aligned bases laid out for one width.
struct Layout {
    key: LayoutKey,
    char_w: f32,
    /// Room for the position numbers on each side of the bases.
    margin: f32,
    /// The alignment's blocks, one per row.
    blocks: Vec<AlignmentBlock>,
    rows: Vec<Row>,
    tops: Vec<f32>,
    total_height: f32,
    /// For every column, the reference base it belongs to: its own, or for an inserted base the
    /// one before it.
    anchors: Vec<usize>,
    /// For every column, the difference it shows, as an index into the alignment's differences.
    differences: Vec<Option<usize>>,
}

impl Layout {
    /// The row at a height in the pane.
    fn row_at(&self, y: f32) -> usize {
        match self.tops.binary_search_by(|t| t.partial_cmp(&y).unwrap_or(std::cmp::Ordering::Less)) {
            Ok(ii) => ii,
            Err(ii) => ii.saturating_sub(1),
        }
        .min(self.rows.len().saturating_sub(1))
    }
}

/// Lay the alignment out in rows that fit `width`.
fn build_layout(ui: &egui::Ui, aln: &Alignment, doc: &Document, width: f32, colors: ColorMode, key: LayoutKey) -> Layout {
    let measure = |text: &str, font: FontId| ui.ctx().fonts_mut(|f| f.layout_no_wrap(text.to_string(), font, Color32::WHITE).size().x);
    let char_w = measure("ACGTACGTAC", style::mono_font(BASE_FONT)) / 10.0;
    let margin = measure(&thousands(aln.reference_len.max(aln.read_len)), style::mono_font(11.0)) + 18.0;
    // Whole tens of columns to a row, so positions are easy to count along it. The 12 points
    // leave room for the scroll bar.
    let usable = width - 2.0 * margin - 12.0;
    let per_row = ((usable / char_w).floor().max(10.0) as usize) / 10 * 10;
    let len = doc.len().max(1);
    let n = aln.aligned.len();

    // Anchors: an inserted base belongs to the reference base before it, or to the first one when
    // the alignment starts with an insertion.
    let mut anchors = vec![0usize; n];
    let mut last = aln.aligned.iter().find_map(|c| c.reference).unwrap_or(0);
    for (k, column) in aln.aligned.iter().enumerate() {
        if let Some(r) = column.reference {
            last = r;
        }
        anchors[k] = last;
    }
    let mut differences = vec![None; n];
    for (ii, d) in aln.differences.iter().enumerate() {
        for k in d.columns.clone().filter(|&k| k < n) {
            differences[k] = Some(ii);
        }
    }

    let features: Vec<&Feature> = doc.features.iter().filter(|f| !f.kind.eq_ignore_ascii_case("source") && !f.segments.is_empty()).collect();
    let covers = |f: &Feature, p: usize| f.segments.iter().any(|s| !s.gap && s.span.contains(p, len));
    let blocks = aln.blocks(per_row);
    let mut rows = Vec::with_capacity(blocks.len());
    for block in &blocks {
        let AlignmentBlock::Bases(range) = block else {
            rows.push(Row { bars: Vec::new(), lanes: 0, labels: Vec::new(), height: EVENT_H });
            continue;
        };
        // Feature bars: each run of columns whose bases lie in a feature, with an arrowhead where
        // the run holds the feature's leading end.
        let mut bars = Vec::new();
        for feature in &features {
            let rgba = feature_rgba(feature, colors);
            let extent = feature.span(len);
            let (first_base, last_base) = (extent.start % len, (extent.end + len - 1) % len);
            let mut k = range.start;
            while k < range.end {
                if !covers(feature, anchors[k]) {
                    k += 1;
                    continue;
                }
                let start = k;
                while k < range.end && covers(feature, anchors[k]) {
                    k += 1;
                }
                bars.push(RowBar {
                    uid: feature.uid,
                    start: start - range.start,
                    end: k - range.start,
                    lane: 0,
                    color: Color32::from_rgb(rgba.r, rgba.g, rgba.b),
                    label: feature.name.clone(),
                    arrow_right: feature.strand == Strand::Forward && aln.aligned[k - 1].reference == Some(last_base),
                    arrow_left: feature.strand == Strand::Reverse && aln.aligned[start].reference == Some(first_base),
                });
            }
        }
        // Lanes: each bar goes into the first lane where it overlaps nothing.
        bars.sort_by_key(|b| (b.start, std::cmp::Reverse(b.end)));
        let mut lane_ends: Vec<usize> = Vec::new();
        for bar in &mut bars {
            let lane = lane_ends.iter().position(|&e| e <= bar.start).unwrap_or_else(|| {
                lane_ends.push(0);
                lane_ends.len() - 1
            });
            lane_ends[lane] = bar.end;
            bar.lane = lane;
        }
        // Labels: one for each difference that starts in this row, left out where it would run
        // into the label before it (the band and the tooltip still show the difference).
        let mut labels = Vec::new();
        let mut free_from = f32::MIN;
        for d in aln.differences.iter().filter(|d| range.contains(&d.columns.start)) {
            let text = d.short_change();
            let x = (d.columns.start - range.start) as f32 * char_w;
            if x >= free_from {
                free_from = x + measure(&text, style::ui_font(LABEL_FONT)) + 6.0;
                labels.push((d.columns.start - range.start, text, d.color_kind()));
            }
        }
        let lanes = lane_ends.len();
        let labels_h = if labels.is_empty() { 0.0 } else { LABEL_H };
        let height = ROW_PAD_TOP + lanes as f32 * (BAR_H + LANE_GAP) + 2.0 * LINE_H + MARK_H + labels_h + ROW_PAD_BOTTOM;
        rows.push(Row { bars, lanes, labels, height });
    }
    let mut tops = Vec::with_capacity(rows.len());
    let mut y = 0.0;
    for row in &rows {
        tops.push(y);
        y += row.height;
    }
    Layout { key, char_w, margin, blocks, rows, tops, total_height: y, anchors, differences }
}

/// One row of bases. From the top: the features, the document's bases, a tick for every base that
/// agrees, and the aligned sequence's bases, then the labels. A band in its kind's colour lies
/// behind each difference, and the selected reference bases are shaded.
#[allow(clippy::too_many_arguments)]
fn paint_bases(
    painter: &egui::Painter,
    layout: &Layout,
    index: usize,
    columns: Range<usize>,
    origin: Pos2,
    aln: &Alignment,
    selection: Option<Span>,
    dark: bool,
    palette: &Palette,
) {
    let row = &layout.rows[index];
    let char_w = layout.char_w;
    let x0 = origin.x + layout.margin;
    let x = |k: usize| x0 + k as f32 * char_w;
    let slice = &aln.aligned[columns.clone()];
    let font = style::mono_font(BASE_FONT);

    // Feature bars, lane by lane, with the name inside when it fits.
    let bars_top = origin.y + ROW_PAD_TOP;
    for bar in &row.bars {
        let top = bars_top + bar.lane as f32 * (BAR_H + LANE_GAP);
        let r = Rect::from_min_max(pos2(x(bar.start), top), pos2(x(bar.end), top + BAR_H));
        let head = (BAR_H * 0.55).min(r.width() * 0.5);
        let points = if bar.arrow_right {
            vec![r.left_top(), pos2(r.right() - head, r.top()), pos2(r.right(), r.center().y), pos2(r.right() - head, r.bottom()), r.left_bottom()]
        } else if bar.arrow_left {
            vec![pos2(r.left() + head, r.top()), r.right_top(), r.right_bottom(), pos2(r.left() + head, r.bottom()), pos2(r.left(), r.center().y)]
        } else {
            vec![r.left_top(), r.right_top(), r.right_bottom(), r.left_bottom()]
        };
        painter.add(egui::Shape::convex_polygon(points, bar.color, Stroke::NONE));
        let light = 0.2126 * bar.color.r() as f32 + 0.7152 * bar.color.g() as f32 + 0.0722 * bar.color.b() as f32 > 140.0;
        let text_color = if light { Color32::from_rgb(0x22, 0x22, 0x22) } else { Color32::WHITE };
        let galley = painter.layout_no_wrap(bar.label.clone(), style::ui_font(10.0), text_color);
        if galley.size().x + 8.0 < r.width() {
            painter.galley(pos2(r.center().x - galley.size().x / 2.0, r.center().y - galley.size().y / 2.0), galley, text_color);
        }
    }

    let top_y = bars_top + row.lanes as f32 * (BAR_H + LANE_GAP);
    let mark_y = top_y + LINE_H;
    let bottom_y = mark_y + MARK_H;
    let column_rect = |k: usize| Rect::from_min_max(pos2(x(k), top_y - 1.0), pos2(x(k + 1), bottom_y + LINE_H + 1.0));
    // The selected reference bases, shaded behind both sequences.
    if let Some(sel) = selection {
        let len = aln.reference_len;
        for (k, column) in columns.clone().enumerate() {
            if sel.contains(layout.anchors[column], len) {
                painter.rect_filled(column_rect(k), 0.0, palette.accent_soft);
            }
        }
    }
    // A band behind every difference, in its kind's colour.
    for (k, column) in slice.iter().enumerate() {
        if let Some(kind) = column.kind() {
            painter.rect_filled(column_rect(k), CornerRadius::same(2), kind_color(kind, dark).gamma_multiply(0.28));
        }
    }
    // The two sequences: bases that agree in ink, differing ones in their kind's colour, gaps faint.
    let tint = |column: &AlignedColumn, base: u8| {
        if base == b'-' { palette.faint } else { column.kind().map_or(palette.ink, |kind| kind_ink(kind, dark)) }
    };
    coloured_line(painter, pos2(x0, top_y + 1.0), slice.iter().map(|c| (c.bases[0], tint(c, c.bases[0]))), &font);
    coloured_line(painter, pos2(x0, bottom_y + 1.0), slice.iter().map(|c| (c.bases[1], tint(c, c.bases[1]))), &font);
    // Between them: a tick where the bases agree, and a block in the kind's colour where they differ.
    for (k, column) in slice.iter().enumerate() {
        let cx = x(k) + char_w / 2.0;
        match column.kind() {
            None => {
                painter.line_segment([pos2(cx, mark_y + 2.0), pos2(cx, mark_y + MARK_H - 2.0)], Stroke::new(1.0, palette.faint));
            }
            Some(kind) => {
                let block = Rect::from_center_size(pos2(cx, mark_y + MARK_H / 2.0), vec2(char_w * 0.7, MARK_H - 3.0));
                painter.rect_filled(block, CornerRadius::same(1), kind_color(kind, dark));
            }
        }
    }
    // 1-based positions at both ends of each sequence's line.
    let ends = aln.ends(columns.clone());
    let small = style::mono_font(11.0);
    let x1 = x(slice.len());
    for (ends, line_y) in [(ends.reference, top_y), (ends.read, bottom_y)] {
        if let Some((first, last)) = ends {
            let y = line_y + LINE_H / 2.0;
            painter.text(pos2(x0 - 8.0, y), Align2::RIGHT_CENTER, thousands(first + 1), small.clone(), palette.faint);
            painter.text(pos2(x1 + 8.0, y), Align2::LEFT_CENTER, thousands(last + 1), small.clone(), palette.faint);
        }
    }
    // What changed, under each difference.
    let label_y = bottom_y + LINE_H + LABEL_H / 2.0;
    for (k, text, kind) in &row.labels {
        painter.text(pos2(x(*k), label_y), Align2::LEFT_CENTER, text, style::ui_font(LABEL_FONT), kind_ink(*kind, dark));
    }
}

/// A structural change in place of its bases: a box in its kind's colour, with a headline in the
/// two sequences' names over the change's full description.
#[allow(clippy::too_many_arguments)]
fn paint_event(
    painter: &egui::Painter,
    layout: &Layout,
    origin: Pos2,
    width: f32,
    aln: &Alignment,
    index: usize,
    names: (&str, &str),
    dark: bool,
    palette: &Palette,
) {
    let event = &aln.events[index];
    let kind = event.color_kind();
    let frame = Rect::from_min_max(pos2(origin.x + layout.margin, origin.y + 8.0), pos2(origin.x + width - layout.margin, origin.y + EVENT_H - 8.0));
    painter.rect_filled(frame, 6.0, kind_color(kind, dark).gamma_multiply(0.14));
    painter.rect_stroke(frame, 6.0, Stroke::new(1.5, kind_color(kind, dark)), StrokeKind::Inside);
    let clip = painter.with_clip_rect(frame.shrink(4.0));
    let headline = event_headline(event, names.0, names.1, aln.reference_len);
    clip.text(frame.left_top() + vec2(12.0, 14.0), Align2::LEFT_CENTER, headline, style::bold_font(12.5), kind_ink(kind, dark));
    clip.text(frame.left_top() + vec2(12.0, 32.0), Align2::LEFT_CENTER, &event.description, style::ui_font(11.0), palette.muted);
}

/// A structural change in a few words, in the names of the map (`reference`) and the aligned
/// sequence (`read`).
fn event_headline(event: &StructuralEvent, reference: &str, read: &str, len: usize) -> String {
    let span = event.reference;
    let bp = |n: usize| format!("{} bp", thousands(n));
    match event.kind {
        EventKind::Deletion => format!("Missing from {read}: {} of {reference} ({})", span.format(len), bp(span.len())),
        EventKind::Replacement => {
            format!(
                "Missing from {read}: {} of {reference} ({}), with {} of other sequence in its place",
                span.format(len),
                bp(span.len()),
                bp(event.read_bases)
            )
        }
        EventKind::Insertion => {
            let after = if span.start == 0 { len } else { span.start };
            format!("Only in {read}: {} after base {} of {reference}", bp(event.read_bases), thousands(after))
        }
        EventKind::Duplication => format!("Only in {read}: {} that repeat part of {reference}", bp(event.read_bases)),
        EventKind::Inversion => format!("Turned round in {read}: {} of {reference} ({})", span.format(len), bp(span.len())),
    }
}

/// Paint a line of characters, each with its own colour, as one piece of text. Runs of one colour
/// share a section, so a line costs a handful of sections rather than one per base.
fn coloured_line(painter: &egui::Painter, at: Pos2, characters: impl Iterator<Item = (u8, Color32)>, font: &FontId) {
    let mut job = LayoutJob::default();
    let mut run = String::new();
    let mut run_colour: Option<Color32> = None;
    for (character, colour) in characters {
        if let Some(previous) = run_colour
            && previous != colour
        {
            job.append(&run, 0.0, TextFormat::simple(font.clone(), previous));
            run.clear();
        }
        run_colour = Some(colour);
        run.push(character as char);
    }
    if let Some(colour) = run_colour {
        job.append(&run, 0.0, TextFormat::simple(font.clone(), colour));
    }
    painter.galley(at, painter.layout_job(job), Color32::PLACEHOLDER);
}

/// The row and column under a point of the pane. The column is `None` over a structural change or
/// beside the bases.
fn column_at(layout: &Layout, rect: Rect, p: Pos2) -> (usize, Option<usize>) {
    let ii = layout.row_at(p.y - rect.min.y);
    let AlignmentBlock::Bases(range) = &layout.blocks[ii] else { return (ii, None) };
    let offset = (p.x - rect.min.x - layout.margin) / layout.char_w;
    if offset < 0.0 || offset >= range.len() as f32 {
        return (ii, None);
    }
    (ii, Some(range.start + offset as usize))
}

/// The feature bar under a point, if any.
fn bar_at(layout: &Layout, rect: Rect, p: Pos2, ii: usize) -> Option<u64> {
    let row = &layout.rows[ii];
    let bars_top = rect.min.y + layout.tops[ii] + ROW_PAD_TOP;
    let offset = ((p.x - rect.min.x - layout.margin) / layout.char_w).floor();
    row.bars
        .iter()
        .find(|bar| {
            let top = bars_top + bar.lane as f32 * (BAR_H + LANE_GAP);
            p.y >= top && p.y <= top + BAR_H && offset >= bar.start as f32 && offset < bar.end as f32
        })
        .map(|bar| bar.uid)
}

/// Tooltip text for the point: the structural change, the difference with its features and
/// effect, or the two positions at that column.
fn describe_at(layout: &Layout, aln: &Alignment, rect: Rect, p: Pos2) -> Option<String> {
    let (ii, column) = column_at(layout, rect, p);
    if let AlignmentBlock::Event(e) = &layout.blocks[ii] {
        return Some(aln.events[*e].description.clone());
    }
    let k = column?;
    if let Some(d) = layout.differences[k] {
        return Some(aln.differences[d].describe());
    }
    let c = &aln.aligned[k];
    let position = |p: Option<usize>| p.map(|p| thousands(p + 1)).unwrap_or_else(|| "-".into());
    Some(format!("map {} · aligned sequence {}", position(c.reference), position(c.read)))
}

/// What a click in the pane selects.
enum Click {
    /// These reference bases.
    Bases(Span),
    /// This feature.
    Feature(u64),
}

/// The selection a click asks for: a structural change's stretch, a difference's bases, a
/// feature, or the one reference base under the pointer.
fn click_at(layout: &Layout, aln: &Alignment, rect: Rect, p: Pos2, len: usize) -> Option<Click> {
    let one_base = |at: usize| Span::new(at, at + 1).normalize(len);
    let (ii, column) = column_at(layout, rect, p);
    if let AlignmentBlock::Event(e) = &layout.blocks[ii] {
        let span = aln.events[*e].reference;
        return Some(Click::Bases(if span.is_empty() { one_base(span.start) } else { span }));
    }
    if let Some(uid) = bar_at(layout, rect, p, ii) {
        return Some(Click::Feature(uid));
    }
    let k = column?;
    if let Some(d) = layout.differences[k] {
        let span = aln.differences[d].span(len);
        return Some(Click::Bases(if span.is_empty() { one_base(aln.differences[d].position) } else { span }));
    }
    Some(Click::Bases(one_base(layout.anchors[k])))
}

/// The inspector of the alignment view: the aligned sequences, the result for the one shown, the
/// colour key, the structural changes and differences, and the method.
pub fn inspector(ui: &mut egui::Ui, tab: &mut DocTab, others: &[(u64, String)], settings: &ViewSettings) -> Requests {
    let palette = Palette::get(settings.dark);
    let dark = settings.dark;
    let mut requests = Requests::default();
    tab.alignments.refresh(&tab.doc, tab.revision);
    let len = tab.len();
    // Choices made in the panel, applied once it is drawn: an entry to show, one to remove, and
    // reference bases to select (and bring into view).
    let mut chosen = None;
    let mut removed = None;
    let mut select: Option<Span> = None;
    egui::ScrollArea::vertical().id_salt(("align-inspector", tab.id)).auto_shrink([false, false]).show(ui, |ui| {
        ui.set_width(ui.available_width());
        section_title(ui, "Aligned sequences", &palette);
        ui.horizontal_wrapped(|ui| add_buttons(ui, others, &mut requests));
        if tab.alignments.entries.is_empty() {
            ui.label(RichText::new("Add sequencing results or other sequences to compare with this one.").small().color(palette.faint));
        }
        for (ii, entry) in tab.alignments.entries.iter().enumerate() {
            ui.horizontal(|ui| {
                if ui.selectable_label(ii == tab.alignments.current, &entry.name).clicked() {
                    chosen = Some(ii);
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.small_button("×").on_hover_text("Remove from the list").clicked() {
                        removed = Some(ii);
                    }
                });
            });
            let (verdict, identity) = match &entry.result {
                Some(aln) => (aln.verdict(), format!("{:.2}%", aln.identity())),
                None => ("does not align", String::new()),
            };
            let colour = if verdict.starts_with("exact") || verdict.starts_with("match") { palette.accent } else { palette.warn };
            ui.label(RichText::new(format!("{verdict}  {identity}")).small().color(colour));
        }

        if let Some(aln) = tab.alignments.entries.get(tab.alignments.current).and_then(|e| e.result.as_ref()) {
            let read = &tab.alignments.entries[tab.alignments.current].name;
            let reference = &tab.doc.name;
            ui.separator();
            section_title(ui, "Result", &palette);
            ui.label(RichText::new(aln.verdict()).font(style::bold_font(14.0)).color(palette.ink));
            // How much of the map the sequence has, lacks or has turned round, and what only it
            // has, each beside its colour.
            let tally = aln.tally();
            let share = |n: usize| format!("{:.0}%", n as f64 * 100.0 / len.max(1) as f64);
            let mut lines = vec![(
                aligned_color(dark),
                format!("Aligned: {} bp of {reference} ({}), {:.2}% identical", thousands(tally.aligned), share(tally.aligned), aln.identity()),
            )];
            if tally.missing > 0 {
                lines.push((
                    kind_color(DifferenceKind::Missing, dark),
                    format!("Missing from {read}: {} bp of {reference} ({})", thousands(tally.missing), share(tally.missing)),
                ));
            }
            if tally.inverted > 0 {
                lines.push((kind_color(DifferenceKind::Changed, dark), format!("Turned round in {read}: {} bp", thousands(tally.inverted))));
            }
            if tally.extra > 0 {
                lines.push((kind_color(DifferenceKind::Inserted, dark), format!("Only in {read}: {} bp", thousands(tally.extra))));
            }
            for (colour, text) in lines {
                ui.horizontal_wrapped(|ui| {
                    swatch(ui, colour);
                    ui.label(RichText::new(text).color(palette.ink));
                });
            }
            // A partial read says nothing about the rest of the map.
            if !aln.whole_read {
                let (first, last) = aln.covered().to_1based(len);
                let reach = format!("{read} reaches only {}..{} of {reference}, so the rest of the map is unknown.", thousands(first), thousands(last));
                ui.label(RichText::new(reach).small().color(palette.muted));
            }
            let strand = if aln.reverse { "reverse" } else { "forward" };
            ui.label(RichText::new(format!("{read} runs along the {strand} strand of {reference}.")).small().color(palette.muted));

            ui.separator();
            section_title(ui, "Key", &palette);
            ui.label(RichText::new(format!("On the map, {reference} is the circle and {read} the ring outside it.")).small().color(palette.muted));
            for (colour, words) in key(read, dark) {
                ui.horizontal(|ui| {
                    swatch(ui, colour);
                    ui.label(RichText::new(words).small().color(palette.muted));
                });
            }

            if !aln.events.is_empty() {
                ui.separator();
                section_title(ui, "Structural changes", &palette);
                for event in &aln.events {
                    ui.horizontal_wrapped(|ui| {
                        swatch(ui, kind_color(event.color_kind(), dark));
                        if ui.link(&event.description).clicked() {
                            let span = event.reference;
                            select = Some(if span.is_empty() { Span::new(span.start, span.start + 1).normalize(len) } else { span });
                        }
                    });
                }
            }
            if !aln.differences.is_empty() {
                ui.separator();
                section_title(ui, &format!("Differences ({})", aln.differences.len()), &palette);
                egui::Grid::new(("align-differences", tab.id)).num_columns(3).spacing([8.0, 4.0]).show(ui, |ui| {
                    for d in &aln.differences {
                        if ui.link(thousands(d.position + 1)).clicked() {
                            let span = d.span(len);
                            select = Some(if span.is_empty() { Span::new(d.position, d.position + 1).normalize(len) } else { span });
                        }
                        ui.label(RichText::new(d.short_change()).font(style::mono_font(12.0)).color(kind_ink(d.color_kind(), dark)));
                        let mut note = d.features.join(", ");
                        if !d.effects.is_empty() {
                            note = if note.is_empty() { d.effects.join("; ") } else { format!("{note}: {}", d.effects.join("; ")) };
                        }
                        if d.homopolymer {
                            note.push_str(if note.is_empty() { "homopolymer" } else { ", homopolymer" });
                        }
                        ui.add(egui::Label::new(RichText::new(note).small().color(palette.muted)).wrap());
                        ui.end_row();
                    }
                });
            }
        }

        ui.separator();
        section_title(ui, "Method", &palette);
        ui.add(egui::Label::new(RichText::new(method_summary()).small().color(palette.muted)).wrap());
        ui.add_space(12.0);
    });
    if let Some(ii) = chosen {
        tab.alignments.choose(ii);
    }
    if let Some(ii) = removed {
        tab.alignments.remove(ii);
    }
    if let Some(span) = select {
        tab.select(Some(span));
        tab.alignments.reveal(span.start);
    }
    requests
}

#[cfg(test)]
mod tests {
    use super::*;
    use slopgene_core::model::Topology;

    /// A 2 kb circular reference made of a fixed pseudo-random sequence, with one forward feature
    /// over bases 201..600.
    fn reference() -> DocTab {
        let mut x = 7u64;
        let bases: Vec<u8> = (0..2000)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect();
        let mut doc = Document::new_dna("ref", bases, Topology::Circular);
        doc.features.push(Feature::new("geneA", "CDS", Span::new(200, 600), Strand::Forward));
        DocTab::new(doc, None)
    }

    #[test]
    fn an_edit_to_the_document_aligns_everything_again() {
        let mut tab = reference();
        let read = tab.doc.seq[100..700].to_vec();
        tab.alignments.add(&tab.doc, tab.revision, "read".into(), read, false);
        assert_eq!(tab.alignments.entries[0].result.as_ref().map(|a| a.ref_start), Some(100));
        // Ten bases inserted before the read's stretch move it ten bases on.
        tab.edit(|d| slopgene_core::edit::insert(d, 0, b"ACGTACGTAC")).unwrap();
        tab.alignments.refresh(&tab.doc, tab.revision);
        assert_eq!(tab.alignments.entries[0].result.as_ref().map(|a| a.ref_start), Some(110));
    }

    #[test]
    fn removing_an_entry_keeps_the_one_shown() {
        let mut tab = reference();
        for (ii, start) in [0usize, 500, 1000].iter().enumerate() {
            let read = tab.doc.seq[*start..start + 300].to_vec();
            tab.alignments.add(&tab.doc, tab.revision, format!("read {ii}"), read, false);
        }
        let alignments = &mut tab.alignments;
        assert_eq!(alignments.current, 2, "the newest entry is shown");
        alignments.remove(0);
        assert_eq!(alignments.entries[alignments.current].name, "read 2");
        // Removing the entry shown moves to its neighbour.
        alignments.remove(alignments.current);
        assert_eq!(alignments.entries[alignments.current].name, "read 1");
    }

    #[test]
    fn tracks_follow_the_entries_and_emphasise_the_one_shown() {
        let mut tab = reference();
        let read = tab.doc.seq[100..700].to_vec();
        tab.alignments.add(&tab.doc, tab.revision, "read".into(), read, false);
        tab.alignments.add(&tab.doc, tab.revision, "unrelated".into(), b"ACGT".repeat(100), false);
        let tracks = tab.alignments.tracks();
        assert_eq!(tracks.len(), 2, "one track per entry, even one that does not align");
        assert_eq!((tracks[0].emphasis, tracks[1].emphasis), (false, true));
        assert!(tracks[1].covered.is_empty());
        assert_eq!(tracks[0].covered, Span::new(100, 700));
    }

    #[test]
    fn the_rows_show_features_differences_and_labels() {
        // A read of bases 151..750 with one substitution inside geneA (reference base 300).
        let mut tab = reference();
        let mut read = tab.doc.seq[150..750].to_vec();
        read[150] = if read[150] == b'A' { b'C' } else { b'A' };
        tab.alignments.add(&tab.doc, tab.revision, "read".into(), read, false);
        let aln = tab.alignments.entries[0].result.as_ref().unwrap();
        assert_eq!(aln.differences.len(), 1);
        let ctx = egui::Context::default();
        let mut layout = None;
        let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
            layout = Some(build_layout(ui, aln, &tab.doc, 900.0, ColorMode::Palette, (0, 0, 0, 900, true)));
        });
        output.textures_delta.clear();
        let layout = layout.unwrap();
        // geneA appears as a bar in the rows that hold reference bases 201..600.
        assert!(layout.rows.iter().any(|row| row.bars.iter().any(|b| b.label == "geneA")));
        // The changed base is the one column that points at the difference, and its row labels it.
        let marked: Vec<usize> = (0..aln.aligned.len()).filter(|&k| layout.differences[k].is_some()).collect();
        assert_eq!(marked.len(), 1);
        assert_eq!(aln.aligned[marked[0]].reference, Some(300));
        let labelled: Vec<&(usize, String, DifferenceKind)> = layout.rows.iter().flat_map(|r| r.labels.iter()).collect();
        assert_eq!(labelled.len(), 1);
        assert_eq!(labelled[0].2, DifferenceKind::Changed);
        // Every row is a whole number of tens of columns wide, except the last.
        for block in &layout.blocks[..layout.blocks.len() - 1] {
            if let AlignmentBlock::Bases(range) = block {
                assert_eq!(range.len() % 10, 0);
            }
        }
    }
}
