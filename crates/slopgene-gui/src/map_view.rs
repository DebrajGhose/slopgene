//! The map view: a circular plasmid map or a linear map, drawn from the shared core layout.
//!
//! Interaction: hovering shows what is under the pointer; clicking a feature selects it,
//! clicking an enzyme moves the cursor to its cut, clicking elsewhere moves the cursor there;
//! dragging selects a range (around the circle, through the origin when needed). A click gives
//! the map the keyboard: Delete or Backspace then removes the selected bases and Escape clears
//! the selection. Copy, cut, paste, select all, undo and redo work here through the app. A
//! right click on a selection offers to clone it, to make a feature of it, or to copy it.
//!
//! In the alignment view the map also draws the aligned sequences, as rings with a mark at each
//! difference, and leaves the enzyme sites out. Clicking a mark selects its bases, and clicking a
//! ring shows that sequence in detail.

use crate::app::{Action, ViewSettings};
use crate::doc::{DocTab, View};
use crate::scene_paint::{EguiMeasure, hit, paint};
use crate::style::{self, Palette};
use eframe::egui::{self, Align2, Key, Pos2, Rect, Sense, vec2};
use slopgene_core::render::{MapOptions, Rgba, Scene, Tag, Theme, render_map};
use slopgene_core::*;
use std::f32::consts::TAU;

/// Something the map wants the app to do.
pub enum MapEvent {
    /// The user double-clicked a feature: open its editor.
    EditFeature(u64),
    /// The user chose a command from the map's menu.
    Act(Action),
}

/// Layout margin of linear maps (must match `slopgene_core::render::linear`).
const LINEAR_MARGIN: f32 = 24.0;

fn to_rgba(c: egui::Color32) -> Rgba {
    Rgba::new(c.r(), c.g(), c.b(), c.a())
}

/// Build the scene for the current state of the document.
fn build_scene(ctx: &egui::Context, tab: &mut DocTab, settings: &ViewSettings, size: egui::Vec2, hovered: Option<&Tag>) -> Scene {
    let palette = Palette::get(settings.dark);
    let aligning = tab.view == View::Alignment;
    let marks = if settings.show_enzymes && !aligning { tab.sites().0.to_vec() } else { Vec::new() };
    let alignments = if aligning { tab.alignments.tracks() } else { Vec::new() };
    let mut theme = if settings.dark { Theme::dark() } else { Theme::light() };
    theme.background = to_rgba(palette.canvas);
    let highlight = match hovered {
        Some(Tag::Feature(uid)) => Some(*uid),
        _ => tab.selected_feature,
    };
    let options = MapOptions {
        width: size.x,
        height: size.y,
        theme,
        font_size: 11.5,
        show_primers: settings.show_primers,
        enzymes: marks,
        selection: tab.selection,
        cursor: Some(tab.cursor),
        highlight,
        colors: settings.colors,
        alignments,
        ..MapOptions::default()
    };
    render_map(&tab.doc, &options, &EguiMeasure { ctx })
}

/// Sequence position under a point of a circular map centred in `rect`, whose base `top` sits
/// at 12 o'clock.
pub(crate) fn circular_position(rect: Rect, p: Pos2, len: usize, top: usize) -> usize {
    let center = rect.center();
    let angle = (p.x - center.x).atan2(-(p.y - center.y)).rem_euclid(TAU);
    ((angle / TAU * len as f32).round() as usize + top) % len.max(1)
}

/// How long a circular map takes to turn to a new top, in seconds: quick, but seen.
const TURN_TIME: f64 = 0.3;

/// A circular map on its way round to a new top: the base that stood at the top when the turn
/// began, the base it turns to, and when it began.
#[derive(Clone, Copy)]
struct Turn {
    from: f32,
    to: usize,
    began: f64,
}

/// The base at 12 o'clock this frame of a circular map of `len` bases that is asked to show
/// `target` there. The map starts where it is first asked to. From then on it turns to each new
/// top the short way round over [`TURN_TIME`], fast at first and settling at the end, so the eye
/// follows the part as it moves instead of losing it in a jump. `id` tells one map from another.
pub(crate) fn turned_top(ui: &egui::Ui, id: egui::Id, target: usize, len: usize) -> usize {
    let len = len.max(1);
    let now = ui.input(|i| i.time);
    let mut turn: Turn = ui.data(|d| d.get_temp(id)).unwrap_or(Turn { from: target as f32, to: target, began: now - TURN_TIME });
    // The way from the turn's start to its end, in bases: the short way round, so the map never
    // spins most of a circle to move a little.
    let way = |turn: &Turn| (turn.to as f32 - turn.from + len as f32 / 2.0).rem_euclid(len as f32) - len as f32 / 2.0;
    // Where the turn stands at `now`, eased so that it moves fast first and settles.
    let standing = |turn: &Turn| {
        let part = ((now - turn.began) / TURN_TIME).clamp(0.0, 1.0) as f32;
        (turn.from + way(turn) * (1.0 - (1.0 - part).powi(3))).rem_euclid(len as f32)
    };
    // A new top starts a turn from where the map stands now, even in the middle of a turn.
    if turn.to != target {
        turn = Turn { from: standing(&turn), to: target, began: now };
    }
    let shown = standing(&turn);
    if now - turn.began < TURN_TIME {
        ui.ctx().request_repaint();
    }
    ui.data_mut(|d| d.insert_temp(id, turn));
    (shown.round() as usize) % len
}

/// Sequence position under a point of a linear map drawn in `rect`.
pub(crate) fn linear_position(rect: Rect, p: Pos2, len: usize) -> usize {
    let scale = (rect.width() - 2.0 * LINEAR_MARGIN) / len.max(1) as f32;
    (((p.x - rect.min.x - LINEAR_MARGIN) / scale).round().max(0.0) as usize).min(len)
}

pub fn show(ui: &mut egui::Ui, tab: &mut DocTab, settings: &ViewSettings) -> Option<MapEvent> {
    let palette = Palette::get(settings.dark);
    let available = ui.available_size();
    if tab.doc.is_empty() {
        let (rect, _) = ui.allocate_exact_size(available, Sense::hover());
        ui.painter().text(rect.center(), Align2::CENTER_CENTER, "Empty sequence", style::ui_font(14.0), palette.muted);
        return None;
    }
    let hover_id = egui::Id::unique(("map-hover", tab.id));
    let hovered: Option<Tag> = ui.ctx().data(|d| d.get_temp(hover_id)).flatten();
    let scene = build_scene(ui.ctx(), tab, settings, available, hovered.as_ref());

    // Linear maps grow taller with more feature lanes; scroll when they do not fit.
    let mut event = None;
    if scene.height > available.y + 1.0 {
        egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
            event = interact(ui, tab, &scene, hover_id, vec2(available.x, scene.height));
        });
    } else {
        event = interact(ui, tab, &scene, hover_id, available);
    }
    event
}

fn interact(ui: &mut egui::Ui, tab: &mut DocTab, scene: &Scene, hover_id: egui::Id, size: egui::Vec2) -> Option<MapEvent> {
    // A stable id lets the map keep keyboard focus from frame to frame.
    let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
    let response = ui.interact(rect, egui::Id::unique(("map-view", tab.id)), Sense::click_and_drag());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, crate::scene_paint::color(scene.background));
    paint(&painter, scene, rect.min);

    let measure = EguiMeasure { ctx: ui.ctx() };
    let len = tab.len();
    let circular = tab.doc.is_circular();
    let position_at = |p: Pos2| if circular { circular_position(rect, p, len, 0) } else { linear_position(rect, p, len) };

    // What is under the pointer, remembered for the next frame's highlight.
    let tag = response.hover_pos().and_then(|p| hit(scene, [p.x - rect.min.x, p.y - rect.min.y], &measure));
    ui.ctx().data_mut(|d| d.insert_temp(hover_id, tag.clone()));
    let mut response = response;
    if let Some(t) = &tag {
        let text = describe(tab, t);
        if !text.is_empty() {
            response = response.on_hover_text_at_pointer(text);
        }
    }

    let mut event = None;
    // The menu on a selection: the commands people reach for once they have marked a region.
    if tab.selection.is_some() && tab.view != View::Alignment {
        response.context_menu(|ui| {
            for (label, hint, action) in [
                (
                    "Clone the selection",
                    "Start a cloning protocol that cuts these bases out; the enzymes that cut nearest them are listed to choose from",
                    Action::NewProtocol,
                ),
                ("Add feature from selection…", "Mark these bases as a feature", Action::AddFeatureFromSelection),
                ("Copy", "Copy the bases", Action::CopySelection),
            ] {
                if ui.button(label).on_hover_text(hint).clicked() {
                    event = Some(MapEvent::Act(action));
                    ui.close();
                }
            }
        });
    }
    if response.drag_started()
        && let Some(p) = response.interact_pointer_pos()
    {
        let pos = position_at(p);
        tab.anchor = Some(pos);
        tab.cursor = pos;
        tab.select(None);
    }
    if response.dragged()
        && let (Some(p), Some(anchor)) = (response.interact_pointer_pos(), tab.anchor)
    {
        let pos = position_at(p);
        tab.cursor = pos;
        tab.selection = drag_span(anchor, pos, len, circular);
        tab.selected_feature = None;
    }
    if response.double_clicked()
        && let Some(Tag::Feature(uid)) = &tag
    {
        tab.select_feature(*uid);
        event = Some(MapEvent::EditFeature(*uid));
    } else if response.clicked() {
        match &tag {
            Some(Tag::Feature(uid)) => tab.select_feature(*uid),
            Some(Tag::Enzyme(_, pos)) => {
                tab.select(None);
                tab.cursor = (*pos).min(len);
                tab.scroll_to = Some(tab.cursor);
            }
            Some(Tag::Primer(uid)) => {
                if let Some(site) = tab.doc.primers.iter().find(|p| p.uid == *uid).and_then(|p| p.sites.first()) {
                    let span = site.span;
                    tab.select(Some(span));
                    tab.scroll_to = Some(span.start);
                }
            }
            // A difference: show its sequence, select its bases (the base before an insertion)
            // and bring it into view among the aligned bases.
            Some(Tag::Difference(track, mark)) => {
                if let Some(mark) = tab.alignments.mark(*track, *mark) {
                    let span = if mark.span.is_empty() { Span::new(mark.span.start, mark.span.start + 1).normalize(len) } else { mark.span };
                    tab.alignments.choose(*track);
                    tab.select(Some(span));
                    tab.scroll_to = Some(span.start);
                }
            }
            Some(Tag::Track(track)) => tab.alignments.choose(*track),
            _ => {
                if let Some(p) = response.interact_pointer_pos() {
                    tab.select(None);
                    tab.cursor = position_at(p);
                    tab.scroll_to = Some(tab.cursor);
                }
            }
        }
    }
    if response.clicked() || response.drag_started() {
        response.request_focus();
    }
    if response.has_focus() {
        handle_keys(ui, tab);
    }
    event
}

/// Keys the map handles while it has the keyboard. Typing bases is left to the sequence view,
/// where the cursor is easy to see.
fn handle_keys(ui: &egui::Ui, tab: &mut DocTab) {
    let delete = ui.input(|i| i.key_pressed(Key::Delete) || i.key_pressed(Key::Backspace));
    if delete && tab.selection.is_some() {
        let unit = tab.unit();
        let note = match tab.delete_selection() {
            Ok(count) => crate::seq_view::deleted_note(count, unit),
            Err(e) => e.to_string(),
        };
        tab.notice = Some((note, std::time::Instant::now()));
    }
    if ui.input(|i| i.key_pressed(Key::Escape)) {
        tab.select(None);
    }
}

/// The range swept by a drag from `anchor` to `pos`. On a circle, take whichever way round is
/// shorter, so dragging back across the origin selects through it.
pub(crate) fn drag_span(anchor: usize, pos: usize, len: usize, circular: bool) -> Option<Span> {
    if anchor == pos || len == 0 {
        return None;
    }
    if !circular {
        return Some(Span::new(anchor.min(pos), anchor.max(pos)));
    }
    let forward = (pos + len - anchor) % len;
    let span = if forward <= len / 2 { Span::new(anchor, anchor + forward) } else { Span::new(pos, pos + (len - forward)) };
    Some(span.normalize(len))
}

/// Tooltip text for the thing under the pointer.
fn describe(tab: &DocTab, tag: &Tag) -> String {
    let len = tab.len();
    match tag {
        Tag::Feature(uid) => match tab.doc.feature_index(*uid) {
            Some(index) => {
                let f = &tab.doc.features[index];
                format!("{}\n{} · {} · {} bp", f.name, f.kind, f.span(len).format(len), f.len())
            }
            None => String::new(),
        },
        Tag::Enzyme(name, pos) => {
            let site = enzymes::EnzymeDb::builtin().get(name).map(|e| e.site_with_cuts()).unwrap_or_default();
            format!("{name}  {site}\ncuts after {pos}")
        }
        Tag::Primer(uid) => match tab.doc.primers.iter().find(|p| p.uid == *uid) {
            Some(p) => format!("{} ({} nt)\n{}", p.name, p.seq.len(), p.seq),
            None => String::new(),
        },
        Tag::Difference(track, mark) => tab.alignments.mark(*track, *mark).map(|m| m.description).unwrap_or_default(),
        Tag::Track(track) => match tab.alignments.entries.get(*track) {
            Some(entry) => match &entry.result {
                Some(aln) => format!("{}\n{}, {:.2}% identity", entry.name, aln.verdict(), aln.identity()),
                None => format!("{}\ndoes not align", entry.name),
            },
            None => String::new(),
        },
        _ => String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One frame at `time` of a map of 1000 bases asked to show `target` at the top.
    fn frame(ctx: &egui::Context, time: f64, target: usize) -> usize {
        let raw = egui::RawInput { time: Some(time), ..Default::default() };
        let mut top = 0;
        ctx.run_ui(raw, |ui| top = turned_top(ui, egui::Id::unique("turn-test"), target, 1000)).textures_delta.clear();
        top
    }

    #[test]
    fn a_map_turns_to_a_new_top_the_short_way_round_and_settles() {
        let ctx = egui::Context::default();
        // The first top shows at once, and stays while nothing asks for another.
        assert_eq!(frame(&ctx, 0.0, 100), 100);
        assert_eq!(frame(&ctx, 0.05, 100), 100);
        // A top on the other side of the origin is reached backwards through the origin, the
        // short way, over the turn's time: a third of the way in, the map has already turned
        // most of the way (fast first), and it has settled at the end.
        assert_eq!(frame(&ctx, 0.1, 900), 100, "the turn starts where the map stands");
        let partway = frame(&ctx, 0.2, 900);
        assert!((930..=975).contains(&partway), "partway {partway}");
        assert_eq!(frame(&ctx, 0.4, 900), 900);
        assert_eq!(frame(&ctx, 1.0, 900), 900);
        // A new top in the middle of a turn starts from where the map stands then.
        assert_eq!(frame(&ctx, 1.1, 300), 900);
        let partway = frame(&ctx, 1.2, 300);
        assert!(!(300..=900).contains(&partway), "partway {partway}");
        assert_eq!(frame(&ctx, 1.5, 300), 300);
    }
}
