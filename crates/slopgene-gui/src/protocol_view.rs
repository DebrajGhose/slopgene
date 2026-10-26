//! Protocol tabs: a cloning protocol as a bench of small plasmids, in the order the work goes.
//!
//! Each map the protocol starts from heads a column, drawn as a small ring. Under it sits its
//! cut: the same ring broken where the enzymes cut, with the pieces pulled apart. Green pieces go
//! on into the product and red pieces are left over. Arrows carry the kept pieces down to the
//! join, which draws the new plasmid. The whole protocol is on the one bench. A map whose part is
//! copied out by PCR has a PCR card between it and its cut, which draws the product as a bar: the
//! copy in green with the primers' tails in grey, which the cut then trims off.
//!
//! A click on a ring zooms into its detailed view. A plasmid shows its full map, where a click on
//! a feature marks the part to cut around and a click on an enzyme cuts there. A cut shows its
//! pieces at full size, where a click on a piece keeps it. A join shows the full map of the new
//! plasmid. "All steps" (or Escape) zooms back out. Both zooms animate, so the eye can follow.
//! Under the picture, the molecule's bases are shown as a sequence view with the same stretch
//! highlighted: the part on a plasmid, the piece kept in a cut, the insert in a join. A stretch
//! selected in either view shows in both, and one selected under a plasmid becomes the part.
//!
//! The panel on the right tells the step in plain words (`slopgene_core::protocol::story`) and
//! holds its choices. For a cut, it lists the enzymes that cut nearest the part on each side
//! (`slopgene_core::excise::nearest_cuts`), and the person picks one on each side. Its last
//! section edits the step by hand.
//!
//! The core runs the protocol again whenever it, or an open map it uses, changes. Maps come from
//! an open document of the same file when there is one, so unsaved edits count, and from the
//! file on disk otherwise.

use crate::app::ViewSettings;
use crate::doc::DocTab;
use crate::map_view;
use crate::panels::{ENZYME_SETS, section_title};
use crate::scene_paint::{self, EguiMeasure, hit, paint};
use crate::seq_view;
use crate::style::{self, Palette};
use eframe::egui::{
    self, Align2, Color32, CornerRadius, CursorIcon, DragAndDrop, Key, LayerId, Order, Pos2, Rect, RichText, Sense, Stroke, UiBuilder, pos2, vec2,
};
use slopgene_core::enzymes::{EnzymeDb, cut_counts};
use slopgene_core::excise::{self, NearCut, NearestCuts};
use slopgene_core::gel::{GelOptions, Ladder, Lane, render_gel};
use slopgene_core::protocol::{
    self, Action, AmplifyParams, DigestParams, LigateParams, OpenParams, Orientation, PartRole, PickParams, Piece, Protocol, Run, Status, Step, StepRun,
    split_reference,
};
use slopgene_core::render::ring::{self, RingOptions, RingPiece, Role};
use slopgene_core::render::{MapOptions, Rgba, Scene, SiteMark, Tag, Theme, render_map, thousands};
use slopgene_core::*;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

/// How many changes a protocol can undo.
const UNDO_LIMIT: usize = 100;

/// Changes of one field that follow each other within this time (dragging a number) share one
/// undo step.
const MERGE_WINDOW: Duration = Duration::from_secs(2);

/// How long a zoom between the bench and a detailed view takes, in seconds.
const ZOOM_SECONDS: f32 = 0.35;

/// A card's size on the bench at scale 1, the gaps between columns and rows, and the margin
/// around the bench.
const CARD: egui::Vec2 = vec2(250.0, 232.0);
const COLUMN_GAP: f32 = 70.0;
const ROW_GAP: f32 = 54.0;
const MARGIN: f32 = 22.0;

/// The heights of a card's title line and status line. The ring fills the rest.
const CARD_TITLE: f32 = 28.0;
const CARD_STATUS: f32 = 26.0;

/// The widths at scale 1 of the places on a join's card for one plasmid each: the product or an
/// option, and another plasmid of the tube, which is drawn smaller. A card with one plasmid is as
/// wide as any other card.
const MAIN_SLOT: f32 = CARD.x - 24.0;
const OTHER_SLOT: f32 = 150.0;

/// The most plasmids a join's card shows. Its status line counts the rest.
const SHOWN: usize = 4;

/// The width of the strip beside a join's detailed picture that lists the other plasmids of its
/// tube.
const STRIP: f32 = 300.0;

/// How many enzymes the chooser lists on each side of a part. The detailed map marks the same.
const NEAREST: usize = 8;

/// What a run was made from: the protocol's revision, then each open document's id and revision.
type RunKey = (u64, Vec<(u64, u64)>);

/// What a list of nearest cuts was made from: the run, the cut's step and the enzyme set.
type NearKey = (u64, Option<RunKey>, usize, String);

/// An open protocol.
pub struct ProtocolTab {
    pub id: u64,
    pub protocol: Protocol,
    pub path: Option<PathBuf>,
    pub dirty: bool,
    /// Bumped on every change, so the run knows when to repeat.
    pub revision: u64,
    undo: Vec<Protocol>,
    redo: Vec<Protocol>,
    /// The field and time of the last change that may merge with the next one.
    last_change: Option<(String, Instant)>,
    /// The latest run, and what it was made from.
    pub run: Run,
    run_key: Option<RunKey>,
    /// The step the panel shows, by index.
    pub selected: usize,
    /// The step whose card is shown in detail, or `None` for the bench.
    pub zoom: Option<usize>,
    /// The card the last zoom went into, so zooming out shrinks back into it.
    zoom_last: Option<usize>,
    /// Where each card and its ring sat on the bench in the last frame, by the card's step, for
    /// the zoom: (step, ring, card).
    rings: Vec<(usize, Rect, Rect)>,
    /// The enzyme set the chooser lists.
    pub set: String,
    /// The nearest cuts as last listed, and what they were listed for.
    near: Option<(NearKey, NearestCuts)>,
    /// The enzymes that can go on a PCR's primer tails, as last listed, and what they were
    /// listed for.
    tails: Option<(NearKey, Vec<TailChoice>)>,
    /// The sequence under a zoomed card's picture, kept from frame to frame.
    preview: Option<Preview>,
    /// Whether a zoomed card shows the molecule's bases under its picture.
    pub show_sequence: bool,
    /// The share of a zoomed card's height that the sequence pane takes.
    pub sequence_share: f32,
    /// A short note for the status bar.
    pub notice: Option<(String, Instant)>,
}

/// What a protocol tab asks the app to do.
#[derive(Default)]
pub struct Requests {
    /// Open these molecules, each in a tab of its own.
    pub open: Vec<Document>,
    /// Read the maps from disk again and run every step.
    pub rerun: bool,
}

/// Finds a protocol's maps. An open document of the same file comes first, so its unsaved edits
/// count, and the file on disk otherwise.
struct OpenMaps<'a> {
    folder: PathBuf,
    docs: &'a [DocTab],
}

/// The open document of the file `path` names, read from `folder`.
fn open_tab<'a>(folder: &Path, docs: &'a [DocTab], path: &str) -> Option<&'a DocTab> {
    let file = folder.join(path);
    let wanted = std::fs::canonicalize(&file).ok();
    let same = |other: &Path| match (&wanted, std::fs::canonicalize(other)) {
        (Some(a), Ok(b)) => *a == b,
        _ => file == other,
    };
    docs.iter().find(|t| t.path.as_deref().is_some_and(same))
}

impl protocol::Loader for OpenMaps<'_> {
    fn load(&self, path: &str) -> slopgene_core::Result<Document> {
        match open_tab(&self.folder, self.docs, path) {
            Some(tab) => Ok(tab.doc.clone()),
            None => io::read_one(self.folder.join(path)),
        }
    }
}

impl ProtocolTab {
    pub fn new(protocol: Protocol, path: Option<PathBuf>) -> ProtocolTab {
        ProtocolTab {
            id: new_uid(),
            protocol,
            path,
            dirty: false,
            revision: 0,
            undo: Vec::new(),
            redo: Vec::new(),
            last_change: None,
            run: Run { steps: Vec::new() },
            run_key: None,
            selected: 0,
            zoom: None,
            zoom_last: None,
            rings: Vec::new(),
            set: "common".into(),
            near: None,
            tails: None,
            preview: None,
            show_sequence: true,
            sequence_share: 0.4,
            notice: None,
        }
    }

    pub fn open(path: &Path) -> std::result::Result<ProtocolTab, String> {
        let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let protocol = Protocol::from_json(&text).map_err(|e| e.to_string())?;
        Ok(ProtocolTab::new(protocol, Some(path.to_path_buf())))
    }

    /// The protocol's name, or its file's name.
    pub fn name(&self) -> String {
        if !self.protocol.name.is_empty() {
            return self.protocol.name.clone();
        }
        self.path.as_deref().and_then(Path::file_stem).map(|s| s.to_string_lossy().into_owned()).unwrap_or_else(|| "untitled protocol".into())
    }

    /// The name, with a dot when there are unsaved changes.
    pub fn title(&self) -> String {
        if self.dirty { format!("{} •", self.name()) } else { self.name() }
    }

    /// The folder the protocol's relative paths start from: its file's folder, or the working
    /// folder before it has a file.
    pub fn folder(&self) -> PathBuf {
        match self.path.as_deref().and_then(Path::parent).filter(|p| !p.as_os_str().is_empty()) {
            Some(folder) => folder.to_path_buf(),
            None => std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
        }
    }

    /// Run the protocol again when it, or an open document, changed since the last run. `force`
    /// runs it anyway, which reads the maps on disk again.
    pub fn refresh(&mut self, docs: &[DocTab], force: bool) {
        let key: RunKey = (self.revision, docs.iter().map(|d| (d.id, d.revision)).collect());
        if !force && self.run_key.as_ref() == Some(&key) {
            return;
        }
        let loader = OpenMaps { folder: self.folder(), docs };
        self.run = protocol::run(&self.protocol, &loader, EnzymeDb::builtin());
        self.run_key = Some(key);
        self.selected = self.selected.min(self.protocol.steps.len().saturating_sub(1));
    }

    /// Change the protocol, keeping the old one for undo. Changes with the same `key` in quick
    /// succession share one undo step.
    pub fn change(&mut self, key: Option<String>, change: impl FnOnce(&mut Protocol)) {
        let before = self.protocol.clone();
        change(&mut self.protocol);
        if self.protocol == before {
            return;
        }
        let merge = key.is_some() && self.last_change.as_ref().is_some_and(|(last, at)| Some(last) == key.as_ref() && at.elapsed() < MERGE_WINDOW);
        if !merge {
            self.undo.push(before);
            if self.undo.len() > UNDO_LIMIT {
                self.undo.remove(0);
            }
        }
        self.redo.clear();
        self.last_change = key.map(|k| (k, Instant::now()));
        self.touched();
    }

    /// Mark the protocol changed: unsaved, due a new run, with the step shown still in range.
    fn touched(&mut self) {
        self.dirty = true;
        self.revision += 1;
        let count = self.protocol.steps.len();
        self.selected = self.selected.min(count.saturating_sub(1));
        if self.zoom.is_some_and(|z| z >= count) {
            self.zoom = None;
        }
    }

    pub fn undo(&mut self) -> bool {
        let Some(before) = self.undo.pop() else { return false };
        let replaced = std::mem::replace(&mut self.protocol, before);
        self.restored(replaced.steps.len());
        self.redo.push(replaced);
        true
    }

    pub fn redo(&mut self) -> bool {
        let Some(after) = self.redo.pop() else { return false };
        let replaced = std::mem::replace(&mut self.protocol, after);
        self.restored(replaced.steps.len());
        self.undo.push(replaced);
        true
    }

    /// After undo or redo: the change is new, and a zoom ends when steps came or went, since the
    /// step it showed may now sit elsewhere.
    fn restored(&mut self, steps_before: usize) {
        self.last_change = None;
        if steps_before != self.protocol.steps.len() {
            self.zoom = None;
        }
        self.touched();
    }

    /// Save to the protocol's own file.
    pub fn save(&mut self) -> std::result::Result<PathBuf, String> {
        let Some(path) = self.path.clone() else { return Err("the protocol has no file yet".into()) };
        io::write_atomic(&path, self.protocol.to_json().as_bytes()).map_err(|e| e.to_string())?;
        self.dirty = false;
        Ok(path)
    }

    /// Save to a new file. Map paths are rewritten from the new file's folder, relative where the
    /// maps lie inside it, so they still name the same files.
    pub fn save_as(&mut self, path: &Path) -> std::result::Result<(), String> {
        let mut moved = self.protocol.clone();
        moved.absolute_paths(&self.folder());
        let folder = path.parent().filter(|p| !p.as_os_str().is_empty()).map(Path::to_path_buf).unwrap_or_else(|| PathBuf::from("."));
        moved.relative_paths(&folder);
        io::write_atomic(path, moved.to_json().as_bytes()).map_err(|e| e.to_string())?;
        self.protocol = moved;
        self.path = Some(path.to_path_buf());
        self.dirty = false;
        self.revision += 1;
        Ok(())
    }

    /// Settle step `step`, which waits for a choice, with its option `option` (both 0-based).
    pub fn choose(&mut self, step: usize, option: usize) {
        let Some(result) = self.run.steps.get(step) else { return };
        let Some(choice) = result.choices.get(option) else { return };
        let (id, pin) = (result.id.clone(), choice.pin.clone());
        let mut outcome = Ok(());
        self.change(None, |p| outcome = p.apply(&id, &pin));
        match outcome {
            Ok(()) => self.notice = Some((format!("{id}: option {} chosen", option + 1), Instant::now())),
            Err(e) => self.notice = Some((e.to_string(), Instant::now())),
        }
    }

    /// Add a step of type `kind` after the step shown, and show it. It takes the step before it
    /// as its input, which is what a chain of steps needs most often.
    pub fn add_step(&mut self, kind: &str) {
        let Some(action) = Action::blank(kind) else { return };
        let at = if self.protocol.steps.is_empty() { 0 } else { self.selected + 1 };
        let inputs = match at.checked_sub(1).and_then(|ii| self.protocol.steps.get(ii)) {
            Some(previous) if !matches!(action, Action::Open(_)) => vec![previous.id.clone()],
            _ => Vec::new(),
        };
        let step = Step { id: self.protocol.new_id(kind), inputs, action, note: String::new() };
        self.change(None, |p| p.steps.insert(at, step));
        self.selected = at;
        self.zoom = None;
    }

    /// Delete step `index`. Steps that took its output then say so when they run.
    pub fn remove_step(&mut self, index: usize) {
        self.change(None, |p| {
            if index < p.steps.len() {
                p.steps.remove(index);
            }
        });
        self.zoom = None;
    }

    /// Move step `index` one place up or down, and keep showing it.
    pub fn move_step(&mut self, index: usize, up: bool) {
        let Some(other) = (if up { index.checked_sub(1) } else { Some(index + 1) }).filter(|&ii| ii < self.protocol.steps.len()) else { return };
        self.change(None, |p| p.steps.swap(index, other));
        self.selected = other;
        self.zoom = None;
    }

    /// How the run went, in a few words: "7 of 8 steps done, 1 waits for a choice".
    pub fn progress(&self) -> String {
        let count = |status: Status| self.run.steps.iter().filter(|s| s.status == status).count();
        let mut text = format!("{} of {} steps done", count(Status::Done), self.protocol.steps.len());
        let states = [
            (Status::NeedsChoice, "waits for a choice", "wait for a choice"),
            (Status::Failed, "failed", "failed"),
            (Status::Skipped, "waits for an earlier step", "wait for earlier steps"),
        ];
        for (status, one, several) in states {
            let n = count(status);
            if n > 0 {
                text.push_str(&format!(", {n} {}", if n == 1 { one } else { several }));
            }
        }
        text
    }

    /// The enzymes that cut nearest the part of cut `step`, listed again only when the protocol,
    /// its maps or the enzyme set changed. For the insert's cut, each enzyme also says how often
    /// it cuts the vector, which is usually cut the same way.
    fn nearest(&mut self, step: usize) -> Option<&NearestCuts> {
        let key: NearKey = (self.revision, self.run_key.clone(), step, self.set.clone());
        if self.near.as_ref().is_none_or(|(k, _)| *k != key) {
            let result = self.run.steps.get(step)?;
            let parent = result.inputs.first()?;
            let region = result.findings.region?;
            let enzymes = EnzymeDb::builtin().set(&self.set).unwrap_or_default();
            let insert = self.protocol.steps.get(step).is_some_and(|s| protocol::part_role(&self.protocol, &s.id) == Some(PartRole::Keep));
            let vector = if insert { other_plasmid(&self.protocol, &self.run, step) } else { None };
            let found = excise::nearest_cuts(&parent.doc, region, &enzymes, vector.as_ref(), NEAREST);
            self.near = Some((key, found));
        }
        self.near.as_ref().map(|(_, n)| n)
    }

    /// The enzymes of the set that could go on the primer tails of PCR step `step`: each with
    /// whether it cuts inside the part, which rules it out, and how often it cuts the vector the
    /// copy goes into. Listed again only when the protocol, its maps or the enzyme set changed.
    fn tails(&mut self, step: usize) -> Option<&[TailChoice]> {
        let key: NearKey = (self.revision, self.run_key.clone(), step, self.set.clone());
        if self.tails.as_ref().is_none_or(|(k, _)| *k != key) {
            let result = self.run.steps.get(step)?;
            let template = result.inputs.first()?;
            let region = result.findings.region?;
            let enzymes = EnzymeDb::builtin().set(&self.set).unwrap_or_default();
            let bases = region.extract(&template.doc.seq);
            let inside = cut_counts(&bases, false, &enzymes);
            let vector = other_plasmid(&self.protocol, &self.run, step);
            let in_vector = vector.as_ref().map(|v| cut_counts(&v.seq, v.is_circular(), &enzymes));
            let mut found: Vec<TailChoice> = inside
                .iter()
                .enumerate()
                .map(|(ii, (enzyme, count))| TailChoice {
                    name: enzyme.name.clone(),
                    cuts_part: *count > 0,
                    vector_cuts: in_vector.as_ref().and_then(|v| v.get(ii)).map(|(_, n)| *n),
                })
                .collect();
            // The enzymes that open the vector once come first, as the tails must match the
            // vector's ends, then the rest by name.
            found.sort_by(|a, b| {
                (a.vector_cuts != Some(1)).cmp(&(b.vector_cuts != Some(1))).then_with(|| a.name.to_ascii_lowercase().cmp(&b.name.to_ascii_lowercase()))
            });
            self.tails = Some((key, found));
        }
        self.tails.as_ref().map(|(_, t)| t.as_slice())
    }

    /// Show step `step` in detail, with a zoom from its card.
    pub fn zoom_to(&mut self, step: usize) {
        self.selected = step.min(self.protocol.steps.len().saturating_sub(1));
        self.zoom = Some(step);
        self.zoom_last = Some(step);
    }
}

/// An enzyme of the set, as a choice for a primer's tail.
#[derive(Clone, Debug, PartialEq, Eq)]
struct TailChoice {
    name: String,
    /// Whether it cuts inside the part, in which case it would cut the product too.
    cuts_part: bool,
    /// How often it cuts the vector the copy goes into, when there is one.
    vector_cuts: Option<usize>,
}

/// The join that the pieces of step `step` go into, and which of its inputs they are.
fn feeds(protocol: &Protocol, step: usize) -> Option<(usize, usize)> {
    let mut lineage: Vec<&str> = vec![protocol.steps.get(step)?.id.as_str()];
    for (ii, later) in protocol.steps.iter().enumerate().skip(step + 1) {
        let Some(position) = later.inputs.iter().position(|r| lineage.contains(&split_reference(r).0)) else { continue };
        if matches!(later.action, Action::Ligate(_)) {
            return Some((ii, position));
        }
        lineage.push(&later.id);
    }
    None
}

/// The map a reference's lineage starts from: its first input, followed back to an open step.
fn source_of(protocol: &Protocol, reference: &str) -> Option<usize> {
    let mut id = split_reference(reference).0.to_string();
    // A lineage a few dozen steps long is far longer than any protocol, and stops a loop.
    for _ in 0..32 {
        let index = protocol.steps.iter().position(|s| s.id == id)?;
        let step = &protocol.steps[index];
        if matches!(step.action, Action::Open(_)) {
            return Some(index);
        }
        id = split_reference(step.inputs.first()?).0.to_string();
    }
    None
}

/// The map at the start of the other lane of the join that cut `step` feeds: the vector, for the
/// insert's cut.
fn other_plasmid(protocol: &Protocol, run: &Run, step: usize) -> Option<Document> {
    let (join, mine) = feeds(protocol, step)?;
    protocol.steps[join]
        .inputs
        .iter()
        .enumerate()
        .filter(|(ii, _)| *ii != mine)
        .find_map(|(_, input)| source_of(protocol, input).and_then(|s| run.steps.get(s)?.outputs.first().map(|p| p.doc.clone())))
}

/// The piece the pieces of cut `step` are joined to, as the run made it: the backbone, for the
/// insert's cut.
fn partner(protocol: &Protocol, run: &Run, step: usize) -> Option<Piece> {
    let (join, mine) = feeds(protocol, step)?;
    protocol.steps[join].inputs.iter().enumerate().filter(|(ii, _)| *ii != mine).find_map(|(_, input)| run.resolve(input).first().map(|p| (*p).clone()))
}

/// What a card on the bench stands for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    /// A map the protocol starts from.
    Plasmid,
    /// A cut, with the picks that keep its pieces.
    Cut,
    /// A PCR, which draws its product: the part copied, with the primers' tails on its ends.
    Copy,
    /// A join, which draws the product.
    Join,
    /// Any other step, told in words.
    Other,
    /// A place to start an empty protocol.
    Start,
    /// A place to put the insert into a vector, while no join takes it.
    PutInto,
    /// A place to add another map, once a join exists.
    AddPlasmid,
}

/// One card on the bench.
#[derive(Clone, Debug, PartialEq)]
struct Card {
    kind: Kind,
    /// The step the card stands for. A place to start or to put the insert has none.
    step: Option<usize>,
    /// The picks folded into a cut's card, by step.
    picks: Vec<usize>,
    /// Where the card sits: a column, between two for a join of both, and a row.
    column: f32,
    row: usize,
    /// How many columns wide it is, around its column: 1 for a card, more for a join that shows
    /// several plasmids.
    width: f32,
}

/// An arrow from one card to another, with the steps on the way (a dephosphorylation, a flip).
#[derive(Clone, Debug, PartialEq)]
struct Link {
    from: usize,
    to: usize,
    treatments: Vec<usize>,
}

/// The cards and arrows of a protocol, and the card that shows each step.
#[derive(Clone, Debug, Default, PartialEq)]
struct Board {
    cards: Vec<Card>,
    links: Vec<Link>,
    home: Vec<Option<usize>>,
}

impl Board {
    /// The card that shows step `step`.
    fn card_of(&self, step: usize) -> Option<usize> {
        self.home.get(step).copied().flatten()
    }

    /// The cut card that a card leads into, such as a plasmid's cut.
    fn cut_after(&self, card: usize) -> Option<usize> {
        self.links.iter().find(|l| l.from == card && self.cards[l.to].kind == Kind::Cut).map(|l| l.to)
    }

    /// The plasmid card that leads into cut card `card`, the other way round from `cut_after`.
    fn plasmid_before(&self, card: usize) -> Option<usize> {
        self.links.iter().find(|l| l.to == card && self.cards[l.from].kind == Kind::Plasmid).map(|l| l.from)
    }

    /// The PCR card that a card leads into, such as a plasmid's PCR.
    fn copy_after(&self, card: usize) -> Option<usize> {
        self.links.iter().find(|l| l.from == card && self.cards[l.to].kind == Kind::Copy).map(|l| l.to)
    }

    /// The card that leads into card `card`, whatever it is: the plasmid or the PCR before a cut.
    fn before(&self, card: usize) -> Option<usize> {
        self.links.iter().find(|l| l.to == card).map(|l| l.from)
    }

    /// The step whose part a plasmid card's map sets: the PCR after the plasmid, or else its cut.
    fn part_step_after(&self, card: usize) -> Option<usize> {
        self.copy_after(card).or_else(|| self.cut_after(card)).and_then(|c| self.cards[c].step)
    }

    /// The bench's size at scale 1.
    fn size(&self) -> egui::Vec2 {
        // The right edge of the card that reaches furthest, in columns.
        let columns = self.cards.iter().map(|c| c.column + (c.width + 1.0) / 2.0).fold(1.0, f32::max);
        let rows = self.cards.iter().map(|c| c.row + 1).max().unwrap_or(1) as f32;
        vec2(2.0 * MARGIN + columns * CARD.x + (columns - 1.0) * COLUMN_GAP, 2.0 * MARGIN + rows * CARD.y + (rows - 1.0) * ROW_GAP)
    }
}

/// How many plasmids a join's card shows: its product or its options in full-size places, and
/// the other plasmids of its tube in smaller places after them, [`SHOWN`] in all.
fn join_slots(result: Option<&StepRun>) -> (usize, usize) {
    let Some(result) = result else { return (1, 0) };
    let main = match result.status {
        Status::NeedsChoice => result.choices.len().clamp(1, SHOWN),
        _ => 1,
    };
    let others = match result.status {
        Status::Done | Status::NeedsChoice => result.findings.byproducts.len().min(SHOWN - main),
        _ => 0,
    };
    (main, others)
}

/// How many columns wide a join's card is: one for one plasmid, and more for each one more.
fn join_width(result: Option<&StepRun>) -> f32 {
    let (main, others) = join_slots(result);
    let pixels = (24.0 + main as f32 * MAIN_SLOT + others as f32 * OTHER_SLOT).max(CARD.x);
    (pixels + COLUMN_GAP) / (CARD.x + COLUMN_GAP)
}

/// Add a card for step `step` under the card `source` came from, or in a new column, with an
/// arrow from its source.
fn add_card(board: &mut Board, columns: &mut f32, kind: Kind, step: usize, source: Option<&(usize, Vec<usize>)>) -> usize {
    let (column, row) = match source {
        Some((card, _)) => (board.cards[*card].column, board.cards[*card].row + 1),
        None => {
            *columns += 1.0;
            (*columns - 1.0, 0)
        }
    };
    board.cards.push(Card { kind, step: Some(step), picks: Vec::new(), column, row, width: 1.0 });
    let me = board.cards.len() - 1;
    if let Some((from, treatments)) = source {
        board.links.push(Link { from: *from, to: me, treatments: treatments.clone() });
    }
    me
}

/// Lay a protocol out as a bench. Each map starts a column. A cut sits under what it cuts, and
/// the picks of its pieces fold into it. A dephosphorylation or a flip rides on the arrow that
/// leaves the card. A join sits under its inputs, between their columns, and is as wide as the
/// plasmids `run` says its tube makes.
fn board(protocol: &Protocol, run: &Run) -> Board {
    let steps = &protocol.steps;
    let mut board = Board { home: vec![None; steps.len()], ..Board::default() };
    // Where each step's output is shown: its card, and the treatments it has had since.
    let mut out: Vec<Option<(usize, Vec<usize>)>> = vec![None; steps.len()];
    let mut columns = 0.0;
    let earlier = |ii: usize, reference: &str| steps[..ii].iter().position(|s| s.id == split_reference(reference).0);
    for (ii, step) in steps.iter().enumerate() {
        let source = step.inputs.first().and_then(|r| earlier(ii, r)).and_then(|jj| out[jj].clone());
        let card = match &step.action {
            Action::Open(_) => add_card(&mut board, &mut columns, Kind::Plasmid, ii, None),
            Action::Digest(_) => add_card(&mut board, &mut columns, Kind::Cut, ii, source.as_ref()),
            Action::Amplify(_) => add_card(&mut board, &mut columns, Kind::Copy, ii, source.as_ref()),
            Action::Pick(_) => match &source {
                Some((card, treatments)) if board.cards[*card].kind == Kind::Cut && treatments.is_empty() => {
                    board.cards[*card].picks.push(ii);
                    *card
                }
                _ => add_card(&mut board, &mut columns, Kind::Other, ii, source.as_ref()),
            },
            Action::Dephosphorylate(_) | Action::Flip(_) => match source {
                Some((card, mut treatments)) => {
                    treatments.push(ii);
                    out[ii] = Some((card, treatments));
                    board.home[ii] = Some(card);
                    continue;
                }
                None => add_card(&mut board, &mut columns, Kind::Other, ii, None),
            },
            Action::Ligate(_) => {
                let sources: Vec<(usize, Vec<usize>)> = step.inputs.iter().filter_map(|r| earlier(ii, r)).filter_map(|jj| out[jj].clone()).collect();
                let row = sources.iter().map(|(c, _)| board.cards[*c].row + 1).max().unwrap_or(0);
                let column = if sources.is_empty() {
                    columns += 1.0;
                    columns - 1.0
                } else {
                    sources.iter().map(|(c, _)| board.cards[*c].column).sum::<f32>() / sources.len() as f32
                };
                // The run counts only when it is of this protocol, as it is between an edit and
                // the next run.
                let result = run.steps.get(ii).filter(|r| r.id == step.id);
                board.cards.push(Card { kind: Kind::Join, step: Some(ii), picks: Vec::new(), column, row, width: join_width(result) });
                let me = board.cards.len() - 1;
                for (from, treatments) in sources {
                    board.links.push(Link { from, to: me, treatments });
                }
                me
            }
        };
        out[ii] = Some((card, Vec::new()));
        board.home[ii] = Some(card);
    }
    // An empty protocol starts with a place for a map. A kept piece that no join takes yet gets a
    // place for the plasmid it goes into, and once a join exists a place to add another map sits
    // at the end of the top row.
    let joins = steps.iter().any(|s| matches!(s.action, Action::Ligate(_)));
    let place = |kind: Kind, column: f32| Card { kind, step: None, picks: Vec::new(), column, row: 0, width: 1.0 };
    if board.cards.is_empty() {
        board.cards.push(place(Kind::Start, 0.0));
    } else if !joins && steps.iter().any(|s| matches!(s.action, Action::Pick(_))) {
        board.cards.push(place(Kind::PutInto, columns));
    } else if joins {
        board.cards.push(place(Kind::AddPlasmid, columns));
    }
    // Two cards that overlap in one row: the later one moves down a row.
    for ii in 1..board.cards.len() {
        let overlaps = |c: &Card, d: &Card| c.row == d.row && (c.column - d.column).abs() < (c.width + d.width) / 2.0 - 0.1;
        while board.cards[..ii].iter().any(|c| overlaps(c, &board.cards[ii])) {
            board.cards[ii].row += 1;
        }
    }
    // A wide join can reach left of the first column, and then the whole bench moves right.
    let left = board.cards.iter().map(|c| c.column - (c.width - 1.0) / 2.0).fold(0.0, f32::min);
    for card in &mut board.cards {
        card.column -= left;
    }
    board
}

/// A cut as the bench draws it: the molecule it cuts, its pieces with their roles (and the
/// output each piece is, once cut), its cuts, its part, and the base turned to the top.
struct CutView {
    parent: Document,
    pieces: Vec<(Span, Role, Option<usize>)>,
    cuts: Vec<SiteMark>,
    region: Option<Span>,
    top: usize,
    /// Whether the enzymes have cut, rather than the part being shown before any cut.
    cut: bool,
}

/// The bases of a molecule outside `region`: one stretch round a circle, up to two on a line.
fn complement(region: Span, len: usize, circular: bool) -> Vec<Span> {
    if circular {
        return if region.len() < len { vec![Span::new(region.end, region.start + len).normalize(len)] } else { Vec::new() };
    }
    let mut out = Vec::new();
    if region.start > 0 {
        out.push(Span::new(0, region.start));
    }
    if region.end < len {
        out.push(Span::new(region.end, len));
    }
    out
}

/// The cut a card shows, from the run: its pieces coloured by what the picks kept. Before any cut
/// the part and the rest are coloured by what the cut is for.
fn cut_view(tab: &ProtocolTab, card: &Card) -> Option<CutView> {
    let step = card.step?;
    let result = tab.run.steps.get(step)?;
    let parent = result.inputs.first()?.doc.clone();
    let len = parent.len().max(1);
    let region = result.findings.region;
    // The pieces the picks kept, and those a pick waits to choose among, by where they lay.
    let picks: Vec<&StepRun> = card.picks.iter().filter_map(|&p| tab.run.steps.get(p)).collect();
    let kept: Vec<Span> = picks.iter().filter(|r| r.status == Status::Done).flat_map(|r| r.outputs.iter().filter_map(|o| o.origin.map(|g| g.span))).collect();
    let open: Vec<Span> = picks
        .iter()
        .filter(|r| r.status == Status::NeedsChoice)
        .flat_map(|r| r.choices.iter().flat_map(|c| c.outputs.iter().filter_map(|o| o.origin.map(|g| g.span))))
        .collect();
    let mut view = CutView { parent, pieces: Vec::new(), cuts: Vec::new(), region, top: 0, cut: false };
    if result.status == Status::Done && !result.outputs.is_empty() {
        view.cut = true;
        for (ii, piece) in result.outputs.iter().enumerate() {
            let Some(origin) = piece.origin else { continue };
            let role = if kept.contains(&origin.span) {
                Role::Kept
            } else if open.contains(&origin.span) || card.picks.is_empty() {
                Role::Plain
            } else {
                Role::LeftOver
            };
            view.pieces.push((origin.span, role, Some(ii)));
        }
        for (name, pos) in &result.findings.cuts {
            if !view.cuts.iter().any(|m| m.pos == *pos) {
                view.cuts.push(SiteMark { name: name.clone(), pos: *pos });
            }
        }
    } else if let Some(region) = region {
        let role = tab.protocol.steps.get(step).and_then(|s| protocol::part_role(&tab.protocol, &s.id));
        let (part, rest) = match role {
            Some(PartRole::Keep) => (Role::Kept, Role::LeftOver),
            Some(PartRole::Replace) => (Role::LeftOver, Role::Kept),
            None => (Role::Plain, Role::Plain),
        };
        view.pieces.push((region, part, None));
        view.pieces.extend(complement(region, len, view.parent.is_circular()).into_iter().map(|s| (s, rest, None)));
    }
    // What sits at the top is where the work happens: the part when one is named. Otherwise it is
    // the piece kept, or for a kept backbone (the largest piece) the opening its cuts leave,
    // which is where the insert goes. Failing both, the smallest piece.
    let largest = view.pieces.iter().map(|p| p.0.len()).max().unwrap_or(0);
    let circular = view.parent.is_circular();
    let kept = view.pieces.iter().find(|p| p.1 == Role::Kept).map(|p| p.0);
    let focus = region
        .or_else(|| kept.map(|k| if k.len() == largest && view.pieces.len() > 1 { complement(k, len, circular).first().copied().unwrap_or(k) } else { k }))
        .or_else(|| view.pieces.iter().map(|p| p.0).min_by_key(|s| s.len()));
    view.top = focus.map_or(0, |s| (s.start + s.len() / 2) % len);
    Some(view)
}

/// The part PCR step `step` copies, as found on its template in the run, with the template's
/// length and shape.
fn copy_region(tab: &ProtocolTab, step: usize) -> Option<(Span, usize, bool)> {
    let result = tab.run.steps.get(step)?;
    let template = result.inputs.first()?;
    Some((result.findings.region?, template.doc.len().max(1), template.doc.is_circular()))
}

/// The tails on the product of a PCR, in bases at each end, as the run found them.
fn tail_lengths(result: &StepRun) -> (usize, usize) {
    let [left, right] = result.findings.tails;
    (left, right)
}

/// The middle base of the last part of a product (the insert, in a template's join), to turn to
/// the top.
fn product_top(parts: &[Span], len: usize) -> usize {
    parts.last().map_or(0, |s| (s.start + s.len() / 2) % len.max(1))
}

/// What one plasmid that a join shows stands for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Made {
    /// The join's product.
    Product,
    /// One of the options of a join that waits for a choice, from 0.
    Option(usize),
    /// Another plasmid the join's tube makes, by its place among the by-products, from 0.
    Other(usize),
}

/// One plasmid a join shows, with the stretch of it each piece became.
struct Shown<'a> {
    made: Made,
    piece: &'a Piece,
    parts: &'a [Span],
}

/// The plasmids a join shows: its product, or its options while it waits for a choice, and then
/// the other plasmids its tube makes. A join that is not made shows none.
fn join_products(result: &StepRun) -> Vec<Shown<'_>> {
    let mut shown = Vec::new();
    match result.status {
        Status::Done => {
            if let Some(piece) = result.outputs.first() {
                shown.push(Shown { made: Made::Product, piece, parts: &result.findings.parts });
            }
        }
        Status::NeedsChoice => {
            for (ii, choice) in result.choices.iter().enumerate() {
                if let Some(piece) = choice.outputs.first() {
                    shown.push(Shown { made: Made::Option(ii), piece, parts: &choice.parts });
                }
            }
        }
        Status::Failed | Status::Skipped => return shown,
    }
    for (ii, other) in result.findings.byproducts.iter().enumerate() {
        shown.push(Shown { made: Made::Other(ii), piece: &other.piece, parts: &other.parts });
    }
    shown
}

/// How to draw one plasmid of a join as a ring. The product and the options are green, as they
/// are what the plan makes. The other plasmids of the tube are grey. Each ring turns its last
/// insert to the top, and a backbone closed on itself turns its seam there, where the insert is
/// missing. The joins are marked.
fn made_ring(shown: &Shown, settings: &ViewSettings, background: Color32) -> RingOptions {
    let len = shown.piece.doc.len().max(1);
    let other = matches!(shown.made, Made::Other(_));
    let role = if other { Role::Plain } else { Role::Kept };
    let title = match shown.made {
        Made::Product => Some("New plasmid".to_string()),
        Made::Option(ii) => Some(format!("Option {}", ii + 1)),
        Made::Other(_) => None,
    };
    RingOptions {
        theme: theme_on(settings, background),
        colors: settings.colors,
        gaps: true,
        top: if shown.parts.len() == 1 { 0 } else { product_top(shown.parts, len) },
        pieces: shown.parts.iter().map(|span| RingPiece { span: *span, role, lift: 0.0, label: None }).collect(),
        marks: shown.parts.iter().map(|span| SiteMark { name: String::new(), pos: span.start % len }).collect(),
        title,
        subtitle: Some(format!("{} bp", thousands(shown.piece.doc.len()))),
        ..RingOptions::default()
    }
}

/// What a join's input is, in a word or two for a caption: "geneX piece", "backbone".
fn short_name(tab: &ProtocolTab, reference: &str) -> String {
    let words = protocol::piece_name(&tab.protocol, &tab.run, reference);
    words.strip_prefix("the ").map(str::to_string).unwrap_or(words)
}

/// A few words for the line under one plasmid of join step `step`: an option's number and which
/// way its insert goes, or how another plasmid of the tube differs from the product.
fn made_caption(tab: &ProtocolTab, step: usize, made: Made) -> String {
    let (Some(result), Some(inputs)) = (tab.run.steps.get(step), tab.protocol.steps.get(step).map(|s| &s.inputs)) else { return String::new() };
    let name = |input: usize| inputs.get(input).map(|r| short_name(tab, r)).unwrap_or_default();
    match made {
        Made::Product => String::new(),
        Made::Option(ii) => {
            // The ring already says which option it is, so the words say which way the inserts
            // go: each by name with one insert, and the ones turned round with several.
            let Some(protocol::Pin::Orientations(ways)) = result.choices.get(ii).map(|c| &c.pin) else { return String::new() };
            let turned: Vec<String> = ways.iter().enumerate().filter(|(_, w)| **w == Orientation::Reverse).map(|(jj, _)| name(jj + 1)).collect();
            match (ways.as_slice(), turned.is_empty()) {
                ([way], _) => format!("{} {}", name(1), way.words()),
                (_, true) => "all as given".into(),
                (_, false) => format!("{} turned round", protocol::and_list(&turned)),
            }
        }
        Made::Other(ii) => {
            let Some(other) = result.findings.byproducts.get(ii) else { return String::new() };
            if other.alone() {
                return "no insert".into();
            }
            let mut said: Vec<String> = Vec::new();
            let missing: Vec<String> = other.missing(inputs.len()).into_iter().map(name).collect();
            if !missing.is_empty() {
                said.push(format!("without {}", protocol::and_list(&missing)));
            }
            if other.reordered() {
                said.push("other order".into());
            }
            let turned = other.turned(&protocol::reference_arrangement(result, inputs.len()));
            // With one insert, which piece turned goes without saying.
            match turned.as_slice() {
                [] => {}
                [_] if inputs.len() == 2 => said.push("other way round".into()),
                _ => said.push(format!("{} turned", protocol::and_list(&turned.into_iter().map(name).collect::<Vec<_>>()))),
            }
            said.join(" · ")
        }
    }
}

/// A ring to draw: the molecule and how to draw it, without its size.
struct Picture {
    doc: Document,
    options: RingOptions,
}

/// The theme for drawings on a background of this colour.
fn theme_on(settings: &ViewSettings, background: Color32) -> Theme {
    let mut theme = if settings.dark { Theme::dark() } else { Theme::light() };
    theme.background = Rgba::new(background.r(), background.g(), background.b(), background.a());
    theme
}

/// The ring a card draws: the map whole for a plasmid, with its pieces coloured once its cut is
/// chosen; the pieces pulled apart for a cut; the new plasmid for a join. `detail` asks for the
/// labels of the large view. `None` when there is nothing to draw yet.
fn picture(tab: &ProtocolTab, board: &Board, card: usize, detail: bool, settings: &ViewSettings, background: Color32) -> Option<Picture> {
    let card_ref = &board.cards[card];
    let step = card_ref.step?;
    let result = tab.run.steps.get(step)?;
    let mut options = RingOptions { theme: theme_on(settings, background), colors: settings.colors, gaps: true, ..RingOptions::default() };
    match card_ref.kind {
        Kind::Plasmid => {
            let doc = result.outputs.first()?.doc.clone();
            options.title = Some(doc.name.clone());
            options.subtitle = Some(format!("{} bp", thousands(doc.len())));
            if let Some(view) = board.cut_after(card).and_then(|c| cut_view(tab, &board.cards[c])) {
                options.top = view.top;
                options.gaps = view.cut;
                options.marks = view.cuts;
                options.pieces = view.pieces.iter().map(|(span, role, _)| RingPiece { span: *span, role: *role, lift: 0.0, label: None }).collect();
            } else if let Some((region, len, circular)) = board.copy_after(card).and_then(|c| board.cards[c].step).and_then(|s| copy_region(tab, s)) {
                // A plasmid a PCR copies from shows the part copied. The rest stays as it is, so
                // it is neither kept nor left over.
                options.top = (region.start + region.len() / 2) % len;
                options.gaps = false;
                options.pieces = vec![RingPiece { span: region, role: Role::Kept, lift: 0.0, label: None }];
                options.pieces.extend(complement(region, len, circular).into_iter().map(|s| RingPiece { span: s, role: Role::Plain, lift: 0.0, label: None }));
            }
            Some(Picture { doc, options })
        }
        Kind::Copy => {
            // The product as a bar: the copy of the part in green, with the primers' tails in grey
            // at its ends, where the enzyme sites are.
            let doc = result.outputs.first()?.doc.clone();
            let len = doc.len();
            let (left, right) = tail_lengths(result);
            let (left, right) = (left.min(len), right.min(len.saturating_sub(left)));
            options.title = Some("PCR product".into());
            options.subtitle = Some(format!("{} bp", thousands(len)));
            options.gaps = false;
            if left > 0 {
                options.pieces.push(RingPiece { span: Span::new(0, left), role: Role::Plain, lift: 0.0, label: None });
            }
            options.pieces.push(RingPiece { span: Span::new(left, len - right), role: Role::Kept, lift: 0.0, label: None });
            if right > 0 {
                options.pieces.push(RingPiece { span: Span::new(len - right, len), role: Role::Plain, lift: 0.0, label: None });
            }
            Some(Picture { doc, options })
        }
        Kind::Cut => {
            let view = cut_view(tab, card_ref).filter(|v| v.cut)?;
            let len = view.parent.len();
            let largest = view.pieces.iter().map(|p| p.0.len()).max().unwrap_or(0);
            let part = match &tab.protocol.steps[step].action {
                Action::Digest(d) => d.around.clone(),
                _ => None,
            };
            let few = view.pieces.len() <= 3;
            options.top = view.top;
            options.marks = view.cuts.clone();
            options.mark_labels = detail;
            options.pieces = view
                .pieces
                .iter()
                .map(|(span, role, _)| {
                    // Every piece but the largest pulls away, so the backbone stays where it was.
                    let lift = if span.len() < largest || view.pieces.len() == 1 { if detail { 26.0 } else { 13.0 } } else { 0.0 };
                    let size = format!("{} bp", thousands(span.len()));
                    let holds_part = view.region.is_some_and(|r| span.contains_span(&r, len));
                    let label = match (detail, role) {
                        (true, Role::Kept) => {
                            Some(format!("{}{size} · kept", part.as_ref().filter(|_| holds_part).map(|p| format!("{p} · ")).unwrap_or_default()))
                        }
                        (true, Role::LeftOver) => Some(format!("{size} · left over")),
                        (true, Role::Plain) => Some(size),
                        (false, Role::Kept) => Some(size),
                        (false, _) if few => Some(size),
                        (false, _) => None,
                    };
                    RingPiece { span: *span, role: *role, lift, label }
                })
                .collect();
            Some(Picture { doc: view.parent, options })
        }
        Kind::Join => {
            // The product, or while the direction waits, the first of the options: the plasmid
            // in the card's first place.
            let products = join_products(result);
            let first = products.first()?;
            Some(Picture { doc: first.piece.doc.clone(), options: made_ring(first, settings, background) })
        }
        Kind::Other | Kind::Start | Kind::PutInto | Kind::AddPlasmid => None,
    }
}

/// A change asked for while the tab was drawn, made once drawing is done.
enum Edit {
    /// Show this step in the panel.
    Select(usize),
    /// Zoom into this step's card, or back out to the bench.
    Zoom(Option<usize>),
    /// Settle a step with one of its options (step, option), both 0-based.
    Choose(usize, usize),
    /// Make another plasmid of join step `step`'s tube its product (step, by-product), both
    /// 0-based, by the way each piece of it goes in.
    UseInstead(usize, usize),
    /// Change the protocol. A key lets quick changes of one field share an undo step.
    Change(Option<String>, Box<dyn FnOnce(&mut Protocol)>),
    /// Give a step a new id, and follow every reference to it.
    Rename(String, String),
    AddStep(&'static str),
    RemoveStep(usize),
    /// Move a step one place up (true) or down.
    MoveStep(usize, bool),
    Rerun,
    OpenDocument(Box<Document>),
    /// Choose an enzyme for one side of a cut's part: the cut's step, whether the side is before
    /// the part, and the enzyme.
    PickEnzyme(usize, bool, String),
    /// Give a cut these enzymes.
    SetEnzymes(usize, Vec<String>),
    /// Name the part a cut is about.
    SetPart(usize, Option<String>),
    /// Keep piece `index` of a cut (0-based among its pieces) with pick step `pick`.
    Keep(usize, usize),
    /// Take a step out and pass its input on.
    LeaveOut(usize),
    /// Put the insert into the map these parameters open.
    PutInto(OpenParams),
    /// Put the piece that pick step `pick` keeps into the map that open step `into` starts from.
    PutPiece(usize, usize),
    /// Add the piece that pick step `pick` keeps to join step `join`.
    AddToJoin(usize, usize),
    /// Take piece `index` (0-based) out of join step `join`.
    TakeOut(usize, usize),
    /// Take the map that open step `open` starts from out of the protocol, with its steps.
    RemoveSource(usize),
    /// Remove the 5' phosphates of piece `index` of join step `join` on its way in.
    Dephosphorylate(usize, usize),
    /// Add a map to start from: open it, cut it and keep a piece of it.
    AddSource(OpenParams),
    /// The enzyme set the chooser lists.
    Set(String),
    /// Copy the part of this step's map out by PCR instead of cutting it out.
    CopyOut(usize),
    /// Cut the part of this step's map out instead of copying it by PCR.
    CutOut(usize),
    /// The primers of this PCR step.
    Primers(usize, PrimerEdit),
    /// Choose the enzyme whose site goes on one primer's tail of this PCR step: whether it is the
    /// forward primer's, and the enzyme. The enzyme there already comes off.
    PickTail(usize, bool, String),
    /// Show or hide the molecule's bases under a zoomed card's picture.
    ShowSequence(bool),
    /// The share of a zoomed card's height the sequence pane takes.
    SequenceShare(f32),
}

/// What to do with a PCR's primers.
enum PrimerEdit {
    /// Keep this pair (from 0) of those the run designed.
    Pair(usize),
    /// Use these two, as typed.
    Given(String, String),
    /// Design them anew on each run.
    Clear,
}

/// A change to the step at `index`. A key merges quick changes of one field into one undo step.
fn change_step(index: usize, key: Option<&str>, change: impl FnOnce(&mut Step) + 'static) -> Edit {
    let key = key.map(|k| format!("{index}:{k}"));
    Edit::Change(
        key,
        Box::new(move |p: &mut Protocol| {
            if let Some(step) = p.steps.get_mut(index) {
                change(step);
            }
        }),
    )
}

/// Which of a cut's enzymes makes each end of the piece that holds its part, from the run.
fn flanks(result: Option<&StepRun>, enzymes: &[String]) -> (Option<String>, Option<String>) {
    let (Some(result), Ok(resolved)) = (result, EnzymeDb::builtin().resolve(enzymes)) else { return (None, None) };
    let (Some(parent), Some(region)) = (result.inputs.first(), result.findings.region) else { return (None, None) };
    let near = excise::nearest_cuts(&parent.doc, region, &resolved, None, 1);
    (near.before.first().map(|c| c.enzyme.clone()), near.after.first().map(|c| c.enzyme.clone()))
}

/// A cut's enzymes once `name` is chosen for one side of its part. It takes the place of the
/// enzyme that made that end, unless that enzyme makes the other end too. Choosing the enzyme
/// that already makes that end takes it off, as a second click on a switch would. The enzyme
/// before the part comes first, so the cut reads in the order of its cuts.
fn choose_side(current: &[String], left: Option<String>, right: Option<String>, before: bool, name: &str) -> Vec<String> {
    let same = |a: &str, b: &str| a.eq_ignore_ascii_case(b);
    let (mine, theirs) = if before { (left, right) } else { (right, left) };
    let mut enzymes: Vec<String> = current.to_vec();
    if mine.as_deref().is_some_and(|m| same(m, name)) {
        enzymes.retain(|e| !same(e, name));
        return enzymes;
    }
    if let Some(old) = mine.filter(|m| !theirs.as_deref().is_some_and(|t| same(t, m))) {
        enzymes.retain(|e| !same(e, &old));
    }
    if !enzymes.iter().any(|e| same(e, name)) {
        enzymes.push(name.to_string());
    }
    if before && let Some(at) = enzymes.iter().position(|e| same(e, name)) {
        let first = enzymes.remove(at);
        enzymes.insert(0, first);
    }
    enzymes
}

/// The name of the map that open step `step` starts from: its own name once read, else its
/// file's name, else the step's id.
fn map_name(tab: &ProtocolTab, step: usize) -> String {
    if let Some(made) = tab.run.steps.get(step).and_then(|r| r.outputs.first()) {
        return made.doc.name.clone();
    }
    match tab.protocol.steps.get(step).map(|s| (&s.action, &s.id)) {
        Some((Action::Open(OpenParams { path: Some(path), .. }), _)) => {
            Path::new(path).file_stem().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default()
        }
        Some((_, id)) => id.clone(),
        None => String::new(),
    }
}

/// What the piece that step `step` makes is called in a sentence: "the EGFP piece", "the 944 bp
/// piece", or else "the piece".
fn piece_words(tab: &ProtocolTab, step: usize) -> String {
    if let Some(Action::Pick(params)) = tab.protocol.steps.get(step).map(|s| &s.action)
        && let Some(name) = params.contains.as_ref()
    {
        return format!("the {name} piece");
    }
    match tab.run.steps.get(step).and_then(|r| r.outputs.first()) {
        Some(piece) => format!("the {} bp piece", thousands(piece.doc.len())),
        None => "the piece".into(),
    }
}

/// The name of the map a cut or a PCR works on, from the run, for messages.
fn cut_map_name(tab: &ProtocolTab, cut: usize) -> String {
    let source = tab.protocol.steps.get(cut).and_then(|s| s.inputs.first()).and_then(|r| source_of(&tab.protocol, r));
    source.map(|s| map_name(tab, s)).unwrap_or_else(|| "the map".into())
}

/// The name of the map at the start of step `step`'s lane, whichever step of the lane it is.
fn lane_map_name(tab: &ProtocolTab, step: usize) -> String {
    tab.protocol
        .steps
        .get(step)
        .and_then(|s| tab.protocol.source_of(&s.id))
        .and_then(|id| tab.protocol.steps.iter().position(|s| s.id == id))
        .map(|s| map_name(tab, s))
        .unwrap_or_else(|| "the map".into())
}

/// Tell the person what an edit did, in the status bar, with the way back.
fn tell(tab: &mut ProtocolTab, text: String) {
    tab.notice = Some((format!("{text} Ctrl+Z undoes it."), Instant::now()));
}

/// Make the changes asked for, and pass on what the app has to do. Each change the person makes
/// on the bench or a map is told in the status bar.
fn apply(tab: &mut ProtocolTab, edits: Vec<Edit>) -> Requests {
    let mut requests = Requests::default();
    let report = |tab: &mut ProtocolTab, outcome: slopgene_core::Result<()>| {
        if let Err(e) = outcome {
            tab.notice = Some((e.to_string(), Instant::now()));
        }
    };
    for edit in edits {
        match edit {
            Edit::Select(ii) => tab.selected = ii.min(tab.protocol.steps.len().saturating_sub(1)),
            Edit::Zoom(Some(step)) => tab.zoom_to(step),
            Edit::Zoom(None) => tab.zoom = None,
            Edit::Choose(step, option) => tab.choose(step, option),
            Edit::UseInstead(step, index) => {
                let Some(result) = tab.run.steps.get(step) else { continue };
                let Some(other) = result.findings.byproducts.get(index) else { continue };
                let Some(pin) = other.pin.clone() else { continue };
                // The pieces that now go in the other way round, by what they are.
                let inputs = tab.protocol.steps.get(step).map(|s| s.inputs.clone()).unwrap_or_default();
                let turned: Vec<String> = other
                    .turned(&protocol::reference_arrangement(result, inputs.len()))
                    .into_iter()
                    .filter_map(|input| inputs.get(input).map(|r| protocol::piece_name(&tab.protocol, &tab.run, r)))
                    .collect();
                let id = result.id.clone();
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.apply(&id, &pin));
                if outcome.is_ok() {
                    let verb = if turned.len() == 1 { "goes" } else { "go" };
                    tell(tab, format!("{} now {verb} in the other way round.", uppercase_first(&protocol::and_list(&turned))));
                }
                report(tab, outcome);
            }
            Edit::Change(key, change) => tab.change(key, change),
            Edit::Rename(old, new) => {
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.rename_step(&old, new.trim()));
                report(tab, outcome);
            }
            Edit::AddStep(kind) => tab.add_step(kind),
            Edit::RemoveStep(ii) => tab.remove_step(ii),
            Edit::MoveStep(ii, up) => tab.move_step(ii, up),
            Edit::Rerun => requests.rerun = true,
            Edit::OpenDocument(doc) => requests.open.push(*doc),
            Edit::PickEnzyme(step, before, name) => {
                let Some(Action::Digest(params)) = tab.protocol.steps.get(step).map(|s| &s.action) else { continue };
                let (left, right) = flanks(tab.run.steps.get(step), &params.enzymes);
                let enzymes = choose_side(&params.enzymes, left, right, before, &name);
                set_enzymes(tab, step, enzymes);
            }
            Edit::SetEnzymes(step, enzymes) => set_enzymes(tab, step, enzymes),
            Edit::SetPart(step, part) => {
                let Some(id) = tab.protocol.steps.get(step).map(|s| s.id.clone()) else { continue };
                let mut outcome = Ok(());
                let said = part.clone();
                tab.change(None, |p| outcome = p.set_part(&id, part));
                if outcome.is_ok() {
                    let map = cut_map_name(tab, step);
                    let text = match (said, protocol::part_role(&tab.protocol, &id)) {
                        (Some(part), Some(PartRole::Replace)) => format!("The insert now replaces {part} in {map}."),
                        (Some(part), _) => format!("{part} of {map} is now the part to cut out."),
                        (None, _) => format!("{map} has no part marked now."),
                    };
                    tell(tab, text);
                }
                report(tab, outcome);
            }
            Edit::Keep(pick, piece) => {
                let cut = tab
                    .protocol
                    .steps
                    .get(pick)
                    .and_then(|s| s.inputs.first())
                    .and_then(|r| tab.protocol.steps.iter().position(|s| s.id == split_reference(r).0));
                let params = |cut: usize| tab.protocol.keep_params(&tab.run, &tab.protocol.steps[pick].id, piece).ok().map(|p| (cut, p));
                if let Some((cut, params)) = cut.and_then(params) {
                    let size = tab.run.steps.get(cut).and_then(|r| r.outputs.get(piece)).map(|p| thousands(p.doc.len())).unwrap_or_default();
                    tab.change(None, |p| {
                        if let Some(step) = p.steps.get_mut(pick) {
                            step.action = Action::Pick(params);
                        }
                    });
                    let map = cut_map_name(tab, cut);
                    tell(tab, format!("Keeping the {size} bp piece of {map}."));
                }
            }
            Edit::LeaveOut(step) => {
                let Some(id) = tab.protocol.steps.get(step).map(|s| s.id.clone()) else { continue };
                let what = protocol::story(&tab.protocol, &tab.run, step).headline;
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.leave_out(&id));
                if outcome.is_ok() {
                    tell(tab, format!("Left out: {}.", what.to_lowercase()));
                }
                report(tab, outcome);
                tab.zoom = None;
            }
            Edit::PutInto(vector) => {
                // The last kept piece goes in, cut with the enzymes of its own cut.
                let Some(insert) = tab.protocol.steps.iter().rev().find(|s| matches!(s.action, Action::Pick(_))).map(|s| s.id.clone()) else { continue };
                let enzymes = tab.protocol.enzymes_of_piece(&insert);
                let before = tab.protocol.steps.len();
                tab.change(None, |p| p.put_into(&insert, vector, enzymes, None, true));
                // The vector's steps went in first, so every step moved down by as many.
                let added = tab.protocol.steps.len().saturating_sub(before + 1);
                tab.selected += added;
                tab.zoom = None;
                tell(tab, "The vector is on the bench, cut with the insert's enzymes, and joined to the insert.".into());
            }
            Edit::PutPiece(pick, into) => {
                let (Some(insert), Some(vector)) = (tab.protocol.steps.get(pick).map(|s| s.id.clone()), tab.protocol.steps.get(into).map(|s| s.id.clone()))
                else {
                    continue;
                };
                let (piece, map) = (piece_words(tab, pick), map_name(tab, into));
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.put_piece(&insert, &vector, true).map(|_| ()));
                if outcome.is_ok() {
                    tell(tab, format!("Put {piece} into {map}."));
                }
                report(tab, outcome);
                tab.zoom = None;
            }
            Edit::AddToJoin(pick, join) => {
                let (Some(piece_id), Some(join_id)) = (tab.protocol.steps.get(pick).map(|s| s.id.clone()), tab.protocol.steps.get(join).map(|s| s.id.clone()))
                else {
                    continue;
                };
                let piece = piece_words(tab, pick);
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.add_to_join(&join_id, &piece_id));
                if outcome.is_ok() {
                    tell(tab, format!("Added {piece} to the new plasmid."));
                }
                report(tab, outcome);
                tab.zoom = None;
            }
            Edit::TakeOut(join, index) => {
                let Some(step) = tab.protocol.steps.get(join) else { continue };
                let (join_id, reference) = (step.id.clone(), step.inputs.get(index).cloned().unwrap_or_default());
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.take_out_of_join(&join_id, index));
                if outcome.is_ok() {
                    tell(tab, format!("Took {} out of the new plasmid.", reference.replace('_', " ")));
                }
                report(tab, outcome);
            }
            Edit::RemoveSource(open) => {
                let Some(id) = tab.protocol.steps.get(open).map(|s| s.id.clone()) else { continue };
                let map = map_name(tab, open);
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.remove_source(&id));
                if outcome.is_ok() {
                    tell(tab, format!("Took {map} out of the protocol."));
                }
                report(tab, outcome);
                tab.zoom = None;
            }
            Edit::Dephosphorylate(join, index) => {
                let Some(join_id) = tab.protocol.steps.get(join).map(|s| s.id.clone()) else { continue };
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.dephosphorylate_input(&join_id, index).map(|_| ()));
                if outcome.is_ok() {
                    tell(tab, "That piece now loses its 5' phosphates on the way in, so it cannot close on itself.".into());
                }
                report(tab, outcome);
                tab.zoom = None;
            }
            Edit::AddSource(source) => {
                let mut kept = String::new();
                tab.change(None, |p| kept = p.add_source(source, Vec::new(), None));
                // Zoom straight into the new map, where its part is marked.
                if let Some(open) = tab.protocol.source_of(&kept).and_then(|id| tab.protocol.steps.iter().position(|s| s.id == id)) {
                    tab.zoom_to(open);
                }
                tell(tab, "Added the map. Click a feature or drag across it to mark the part to cut out.".into());
            }
            Edit::Set(name) => tab.set = name,
            Edit::ShowSequence(on) => tab.show_sequence = on,
            Edit::SequenceShare(share) => tab.sequence_share = share.clamp(0.15, 0.75),
            Edit::CopyOut(step) => {
                let Some(id) = tab.protocol.steps.get(step).map(|s| s.id.clone()) else { continue };
                let map = lane_map_name(tab, step);
                let mut outcome = Ok(String::new());
                tab.change(None, |p| outcome = p.copy_out(&id));
                match outcome {
                    Ok(copy) => {
                        // The PCR step is shown, where the tails are chosen.
                        if let Some(index) = tab.protocol.steps.iter().position(|s| s.id == copy) {
                            tab.selected = index;
                        }
                        tab.zoom = None;
                        tell(tab, format!("The part is now copied out of {map} by PCR. Choose the enzymes for the primers' tails in the panel."));
                    }
                    Err(e) => tab.notice = Some((e.to_string(), Instant::now())),
                }
            }
            Edit::CutOut(step) => {
                let Some(id) = tab.protocol.steps.get(step).map(|s| s.id.clone()) else { continue };
                let map = lane_map_name(tab, step);
                let mut outcome = Ok(String::new());
                tab.change(None, |p| outcome = p.cut_out(&id));
                match outcome {
                    Ok(cut) => {
                        if let Some(index) = tab.protocol.steps.iter().position(|s| s.id == cut) {
                            tab.selected = index;
                        }
                        tab.zoom = None;
                        tell(tab, format!("The part is now cut out of {map} again, without a PCR."));
                    }
                    Err(e) => tab.notice = Some((e.to_string(), Instant::now())),
                }
            }
            Edit::Primers(step, what) => {
                let Some(id) = tab.protocol.steps.get(step).map(|s| s.id.clone()) else { continue };
                let (forward, reverse) = match what {
                    PrimerEdit::Pair(ii) => {
                        let Some(pair) = tab.run.steps.get(step).and_then(|r| r.findings.primer_pairs.get(ii)) else { continue };
                        (Some(pair.forward.clone()), Some(pair.reverse.clone()))
                    }
                    PrimerEdit::Given(f, r) => (Some(f), Some(r)),
                    PrimerEdit::Clear => (None, None),
                };
                let kept = forward.is_some();
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.set_primers(&id, forward, reverse));
                if outcome.is_ok() {
                    let text = if kept {
                        "The primers are kept in the protocol, so they stay as they are when the maps change."
                    } else {
                        "The primers are designed anew on each run."
                    };
                    tell(tab, text.into());
                }
                report(tab, outcome);
            }
            Edit::PickTail(step, forward, name) => {
                let Some((id, params)) = tab.protocol.steps.get(step).and_then(|s| match &s.action {
                    Action::Amplify(a) => Some((s.id.clone(), a.clone())),
                    _ => None,
                }) else {
                    continue;
                };
                // A click on the enzyme already there takes it off, as a second click on a switch
                // would.
                let (mut left, mut right) = (params.forward_enzyme.clone(), params.reverse_enzyme.clone());
                let slot = if forward { &mut left } else { &mut right };
                *slot = if slot.as_deref().is_some_and(|c| c.eq_ignore_ascii_case(&name)) { None } else { Some(name) };
                let said = AmplifyParams { forward_enzyme: left.clone(), reverse_enzyme: right.clone(), ..AmplifyParams::default() }.enzymes();
                let mut outcome = Ok(());
                tab.change(None, |p| outcome = p.set_tails(&id, left, right, true));
                if outcome.is_ok() {
                    let map = cut_map_name(tab, step);
                    let text = if said.is_empty() {
                        format!("The primers for {map} have no tails.")
                    } else {
                        format!("The primers for {map} carry {} sites on their tails.", protocol::and_list(&said))
                    };
                    tell(tab, text);
                }
                report(tab, outcome);
            }
        }
    }
    requests
}

/// Give a cut new enzymes. The insert's cut leads: a vector cut with the same enzymes follows it.
fn set_enzymes(tab: &mut ProtocolTab, step: usize, enzymes: Vec<String>) {
    let Some(id) = tab.protocol.steps.get(step).map(|s| s.id.clone()) else { return };
    let follow = protocol::part_role(&tab.protocol, &id) == Some(PartRole::Keep);
    let said = enzymes.clone();
    let mut outcome = Ok(());
    tab.change(None, |p| outcome = p.set_enzymes(&id, enzymes, follow));
    match outcome {
        Ok(()) => {
            let map = cut_map_name(tab, step);
            let text = if said.is_empty() { format!("{map} has no enzymes chosen.") } else { format!("{map} is cut with {}.", protocol::and_list(&said)) };
            tell(tab, text);
        }
        Err(e) => tab.notice = Some((e.to_string(), Instant::now())),
    }
}

fn status_color(status: Status, palette: &Palette) -> Color32 {
    match status {
        Status::Done => palette.ok,
        Status::NeedsChoice => palette.warn,
        Status::Failed => palette.error,
        Status::Skipped => palette.faint,
    }
}

/// The last part of a path, for showing which file a step opens.
fn file_name(path: &str) -> String {
    Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.to_string())
}

/// The cut whose choices the panel and the detailed view offer: the step shown when it is a cut,
/// the cut a pick belongs to, or the cut that follows a map. The cut that trims a PCR's tails has
/// no choices of its own, so it is left out.
fn focused_cut(tab: &ProtocolTab, board: &Board) -> Option<usize> {
    let step = tab.protocol.steps.get(tab.selected)?;
    let cut = match &step.action {
        Action::Digest(_) => Some(tab.selected),
        Action::Pick(_) => board.card_of(tab.selected).and_then(|c| board.cards[c].step).filter(|&s| matches!(tab.protocol.steps[s].action, Action::Digest(_))),
        Action::Open(_) => board.card_of(tab.selected).and_then(|c| board.cut_after(c)).and_then(|c| board.cards[c].step),
        _ => None,
    }?;
    tab.protocol.amplify_before(&tab.protocol.steps[cut].id).is_none().then_some(cut)
}

/// The PCR whose choices the panel and the detailed view offer: the step shown when it is a PCR,
/// the PCR whose product a cut or a pick works on, or the PCR that follows a map.
fn focused_copy(tab: &ProtocolTab, board: &Board) -> Option<usize> {
    let step = tab.protocol.steps.get(tab.selected)?;
    let index_of = |id: &str| tab.protocol.steps.iter().position(|s| s.id == id);
    match &step.action {
        Action::Amplify(_) => Some(tab.selected),
        Action::Digest(_) => tab.protocol.amplify_before(&step.id).and_then(index_of),
        Action::Pick(_) => {
            let cut = board.card_of(tab.selected).and_then(|c| board.cards[c].step)?;
            tab.protocol.amplify_before(&tab.protocol.steps[cut].id).and_then(index_of)
        }
        Action::Open(_) => board.card_of(tab.selected).and_then(|c| board.copy_after(c)).and_then(|c| board.cards[c].step),
        _ => None,
    }
}

/// Draw the protocol tab: the bench (or the card zoomed into) and, when `panel` is true, the
/// panel on the right. `keys` is false while a dialog has the keyboard; otherwise Escape zooms
/// out and the up and down arrows walk through the steps.
pub fn show(ui: &mut egui::Ui, tab: &mut ProtocolTab, docs: &[DocTab], settings: &ViewSettings, keys: bool, panel: bool) -> Requests {
    let palette = Palette::get(settings.dark);
    let mut edits: Vec<Edit> = Vec::new();
    let board = board(&tab.protocol, &tab.run);
    // The nearest cuts for the cut in focus, listed once per change.
    let near = focused_cut(tab, &board).and_then(|cut| tab.nearest(cut).cloned().map(|n| (cut, n)));
    // The enzymes that can go on the tails of the PCR in focus, listed once per change.
    let tails = focused_copy(tab, &board).and_then(|copy| tab.tails(copy).map(|t| (copy, t.to_vec())));
    if panel {
        let width = ui.available_width();
        egui::Panel::right(("protocol-panel", tab.id))
            .resizable(true)
            .default_size((width * 0.34).clamp(300.0, 420.0))
            .min_size(260.0)
            .frame(egui::Frame::new().fill(palette.panel).inner_margin(egui::Margin::symmetric(14, 10)))
            .show(ui, |ui| step_panel(ui, tab, &board, docs, near.as_ref(), tails.as_ref(), &palette, &mut edits));
    }

    // The bench and the detailed view share the rest. While a zoom runs, the bench fades and the
    // card's ring grows (or shrinks) between its place on the bench and the detailed view.
    let area = ui.available_rect_before_wrap();
    ui.allocate_rect(area, Sense::hover());
    let t = ui.ctx().animate_bool_with_time_and_easing(
        egui::Id::unique(("protocol-zoom", tab.id)),
        tab.zoom.is_some(),
        ZOOM_SECONDS,
        egui::emath::easing::cubic_in_out,
    );
    let zoomed = tab.zoom.or(tab.zoom_last).and_then(|step| board.card_of(step));
    let mut rings = Vec::new();
    if t < 1.0 {
        let mut bench_ui = ui.new_child(UiBuilder::new().max_rect(area).id_salt(("protocol-bench", tab.id)));
        bench_ui.multiply_opacity(1.0 - t);
        rings = bench(&mut bench_ui, tab, &board, docs, settings, &palette, &mut edits, t == 0.0);
    }
    if t > 0.0
        && let Some(card) = zoomed
    {
        // The sequence pane under the picture, when the view shows one, takes its share of the
        // area. The pane's document leaves the tab while the view is drawn.
        let mut preview = tab.preview.take();
        let picture_rect = picture_rect(tab, &board, card, area);
        // Partway, the card grows from its place on the bench to fill the view, carrying its ring
        // to where the detailed drawing will put it. The detailed view fades in over the last
        // stretch, and the ring fades out under it.
        if t < 1.0 {
            let step = board.cards[card].step.unwrap_or(usize::MAX);
            let from = rings.iter().chain(tab.rings.iter()).find(|(s, _, _)| *s == step).map(|(_, ring, card)| (*ring, *card));
            if let (Some((ring_from, card_from)), Some(shown)) = (from, picture(tab, &board, card, false, settings, palette.canvas)) {
                let lerp = |a: Rect, b: Rect| Rect::from_min_max(a.min.lerp(b.min, t), a.max.lerp(b.max, t));
                let mut painter = ui.painter_at(area);
                painter.multiply_opacity(1.0 - smoothstep(0.75, 1.0, t));
                let backdrop = lerp(card_from, area);
                let corner = (10.0 * (1.0 - t)).round() as u8;
                painter.rect(backdrop, CornerRadius::same(corner), palette.canvas, Stroke::new(1.0, palette.line), egui::StrokeKind::Inside);
                let target = match board.cards[card].kind {
                    Kind::Join => join_ring_target(tab.run.steps.get(step), picture_rect),
                    kind => detail_ring_rect(kind, picture_rect),
                };
                let rect = lerp(ring_from, target);
                let mut options = shown.options;
                options.width = rect.width();
                options.height = rect.height();
                options.font_size = 10.5 + 1.5 * t;
                let scene = ring::render_ring(&shown.doc, &options, &EguiMeasure { ctx: ui.ctx() });
                paint(&painter, &scene, rect.min);
            }
        }
        let mut detail_ui = ui.new_child(UiBuilder::new().max_rect(area).id_salt(("protocol-detail", tab.id)));
        detail_ui.multiply_opacity(smoothstep(0.55, 1.0, t));
        detail(&mut detail_ui, tab, &board, card, near.as_ref(), settings, &palette, &mut edits, t == 1.0, &mut preview);
        tab.preview = preview;
    }
    if !rings.is_empty() {
        tab.rings = rings;
    }

    // Escape zooms out; the arrows walk through the steps, and a zoom follows them. The sequence
    // pane keeps the arrows while it has the keyboard.
    let sequence_focused = tab.preview.as_ref().is_some_and(|p| ui.ctx().memory(|m| m.has_focus(egui::Id::unique(("seq-view", p.tab.id)))));
    if keys && !ui.ctx().text_edit_focused() && !tab.protocol.steps.is_empty() {
        let (up, down, escape) = ui.input(|i| (i.key_pressed(Key::ArrowUp), i.key_pressed(Key::ArrowDown), i.key_pressed(Key::Escape)));
        if escape && tab.zoom.is_some() {
            edits.push(Edit::Zoom(None));
        }
        let next = if sequence_focused {
            None
        } else if up && tab.selected > 0 {
            Some(tab.selected - 1)
        } else if down && tab.selected + 1 < tab.protocol.steps.len() {
            Some(tab.selected + 1)
        } else {
            None
        };
        if let Some(next) = next {
            edits.push(Edit::Select(next));
            if tab.zoom.is_some()
                && let Some(step) = board.card_of(next).and_then(|c| board.cards[c].step)
            {
                edits.push(Edit::Zoom(Some(step)));
            }
        }
    }
    apply(tab, edits)
}

/// Hermite smoothstep: 0 below `a`, 1 above `b`, and a smooth rise between.
fn smoothstep(a: f32, b: f32, t: f32) -> f32 {
    let x = ((t - a) / (b - a)).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// The height of the detailed view's header: the way back, the step's headline and a hint.
const DETAIL_HEADER: f32 = 74.0;

/// Where the detailed view draws its picture.
fn detail_picture_rect(area: Rect) -> Rect {
    Rect::from_min_max(pos2(area.min.x + 12.0, area.min.y + DETAIL_HEADER), pos2(area.max.x - 12.0, area.max.y - 10.0))
}

/// Where the ring of a zoomed card ends up: a full map's circle for a plasmid or a product, the
/// whole picture for a cut. A ring as large as the map's circle has the canvas this returns.
fn detail_ring_rect(kind: Kind, picture: Rect) -> Rect {
    match kind {
        Kind::Cut => picture,
        _ => {
            // A full circular map.s radius is 0.31 of its height or 0.23 of its width, whichever
            // is less (`render::circular`). A ring.s radius is a share of its canvas.s smaller side.
            let radius = (picture.height() * 0.31).min(picture.width() * 0.23);
            Rect::from_center_size(picture.center(), egui::Vec2::splat(radius / ring::RADIUS_SHARE))
        }
    }
}

/// The room under each option's ring in a join's detailed view, for its words and its button.
const CHOICE_BELOW: f32 = 90.0;

/// The part of a join's detailed picture that its product or its options take: all of it, or all
/// but the strip on the right when the tube makes other plasmids too.
fn join_main_rect(result: &StepRun, picture: Rect) -> Rect {
    let shows_others = matches!(result.status, Status::Done | Status::NeedsChoice) && !result.findings.byproducts.is_empty();
    if shows_others { Rect::from_min_max(picture.min, pos2(picture.max.x - STRIP, picture.max.y)) } else { picture }
}

/// Where the first ring of a zoomed join ends up: the circle of the product's full map, or the
/// first option's ring while the direction waits for a choice.
fn join_ring_target(result: Option<&StepRun>, picture: Rect) -> Rect {
    let Some(result) = result else { return detail_ring_rect(Kind::Join, picture) };
    let main = join_main_rect(result, picture);
    match result.status {
        Status::NeedsChoice => {
            let count = result.choices.len().max(1) as f32;
            Rect::from_min_size(main.min, vec2(main.width() / count, main.height() - CHOICE_BELOW))
        }
        _ => detail_ring_rect(Kind::Join, main),
    }
}

/// Draw the bench in a scrolling area: the arrows, then the cards, scaled to fill the room.
/// Returns where each card and its ring were drawn, for the zoom.
#[allow(clippy::too_many_arguments)]
fn bench(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    board: &Board,
    docs: &[DocTab],
    settings: &ViewSettings,
    palette: &Palette,
    edits: &mut Vec<Edit>,
    interactive: bool,
) -> Vec<(usize, Rect, Rect)> {
    let mut rings = Vec::new();
    let room = ui.available_size();
    let size = board.size();
    let scale = (room.x / size.x).min(room.y / size.y).clamp(0.8, 1.45);
    // Dragging on the bench carries pieces, so it never scrolls the bench, even on a touch screen.
    let scrolling = egui::scroll_area::ScrollSource { drag: egui::scroll_area::DragScroll::Never, ..Default::default() };
    egui::ScrollArea::both().id_salt(("protocol-bench-scroll", tab.id)).auto_shrink([false, false]).scroll_source(scrolling).show(ui, |ui| {
        let width = (size.x * scale).max(ui.available_width());
        let (rect, _) = ui.allocate_exact_size(vec2(width, size.y * scale), Sense::hover());
        // A bench narrower than the room sits in its middle.
        let origin = pos2(rect.min.x + (width - size.x * scale) / 2.0, rect.min.y);
        // A card wider than one column spreads evenly to both sides of its column.
        let card_rect = |card: &Card| {
            let left = card.column - (card.width - 1.0) / 2.0;
            let min = origin + vec2(MARGIN + left * (CARD.x + COLUMN_GAP), MARGIN + card.row as f32 * (CARD.y + ROW_GAP)) * scale;
            Rect::from_min_size(min, vec2(card.width * (CARD.x + COLUMN_GAP) - COLUMN_GAP, CARD.y) * scale)
        };
        for (ii, link) in board.links.iter().enumerate() {
            let incoming: Vec<usize> = board.links.iter().enumerate().filter(|(_, l)| l.to == link.to).map(|(jj, _)| jj).collect();
            let slot = (incoming.iter().position(|&jj| jj == ii).unwrap_or(0), incoming.len());
            draw_link(ui, tab, board, link, card_rect(&board.cards[link.from]), card_rect(&board.cards[link.to]), slot, scale, palette, edits, interactive);
        }
        for (ii, card) in board.cards.iter().enumerate() {
            let rect = card_rect(card);
            if let Some(ring) = draw_card(ui, tab, board, ii, rect, scale, docs, settings, palette, edits, interactive)
                && let Some(step) = card.step
            {
                rings.push((step, ring, rect));
            }
        }
        if interactive {
            draw_carrying(ui, tab, board, card_rect, palette);
        }
    });
    rings
}

/// Shorten text with an ellipsis until it fits `width` at this font.
fn fit(ui: &egui::Ui, text: &str, font: egui::FontId, width: f32) -> String {
    let measure = |s: &str| ui.ctx().fonts_mut(|f| f.layout_no_wrap(s.to_string(), font.clone(), Color32::WHITE).size().x);
    if measure(text) <= width {
        return text.to_string();
    }
    let mut chars: Vec<char> = text.chars().collect();
    while chars.len() > 1 {
        chars.pop();
        let candidate: String = chars.iter().collect::<String>() + "…";
        if measure(&candidate) <= width {
            return candidate;
        }
    }
    String::new()
}

/// A problem in a few words for a card's status line, from its code. The panel gives it whole.
fn short_problem(issue: &protocol::Issue) -> String {
    match issue.code {
        "ends.incompatible" => "The ends do not pair".into(),
        "ends.no_phosphate" => "No 5' phosphate at a join".into(),
        "digest.no_cut" => "The enzymes do not cut it".into(),
        "digest.unknown_enzyme" | "amplify.unknown_enzyme" => "An enzyme is unknown".into(),
        "amplify.no_part" => "The part is not on the map".into(),
        "amplify.no_primers" => "No primers fit the part".into(),
        "amplify.no_product" => "The primers make no product".into(),
        "amplify.stale" => "The primers no longer fit the part".into(),
        "pick.no_match" => "No piece fits".into(),
        "pick.stale" => "The piece kept changed".into(),
        "open.unreadable" | "open.nothing" => "The map cannot be read".into(),
        _ => issue.message.clone(),
    }
}

/// A card's title (a verb, then what it acts on) and its status line with the line's colour.
fn card_text(tab: &ProtocolTab, board: &Board, card: usize, palette: &Palette) -> (String, String, String, Color32) {
    let c = &board.cards[card];
    let result = c.step.and_then(|s| tab.run.steps.get(s));
    let status = result.map(|r| r.status);
    let problem = || result.and_then(|r| r.error.as_ref()).map(short_problem).unwrap_or_default();
    match c.kind {
        Kind::Plasmid => {
            let made = result.and_then(|r| r.outputs.first());
            let name = made.map(|p| p.doc.name.clone()).unwrap_or_else(|| "choose a map".into());
            match (status, made) {
                (Some(Status::Done), Some(p)) => {
                    let shape = if p.doc.is_circular() { "circular" } else { "linear" };
                    ("START".into(), name, format!("{} bp · {shape}", thousands(p.doc.len())), palette.muted)
                }
                _ => ("START".into(), name, problem(), palette.error),
            }
        }
        Kind::Cut => {
            let Some(Action::Digest(params)) = c.step.map(|s| &tab.protocol.steps[s].action) else { return Default::default() };
            let enzymes = if params.enzymes.is_empty() { "choose enzymes".to_string() } else { params.enzymes.join(" + ") };
            let picks: Vec<&StepRun> = c.picks.iter().filter_map(|&p| tab.run.steps.get(p)).collect();
            let (line, color) = match status {
                Some(Status::NeedsChoice) => {
                    let around = params.around.as_deref().map(|p| format!(" around {p}")).unwrap_or_default();
                    (format!("Choose where to cut{around}"), palette.warn)
                }
                Some(Status::Done) if picks.iter().any(|p| p.status == Status::NeedsChoice) => ("Choose the piece to keep".into(), palette.warn),
                Some(Status::Done) => match picks.iter().find(|p| p.status == Status::Done).and_then(|p| p.outputs.first()) {
                    Some(kept) => (format!("keeps {} bp · {} pieces", thousands(kept.doc.len()), result.map_or(0, |r| r.outputs.len())), palette.muted),
                    None => (picks.iter().find_map(|p| p.error.as_ref()).map(short_problem).unwrap_or_default(), palette.error),
                },
                Some(Status::Failed) => (problem(), palette.error),
                _ => ("waits for the step above".into(), palette.faint),
            };
            ("CUT".into(), enzymes, line, color)
        }
        Kind::Copy => {
            let Some(Action::Amplify(params)) = c.step.map(|s| &tab.protocol.steps[s].action) else { return Default::default() };
            let object = match params.region.as_deref() {
                Some(part) => format!("{part} by PCR"),
                None => "by PCR".to_string(),
            };
            let tails = params.enzymes();
            let (line, color) = match status {
                Some(Status::NeedsChoice) => ("Mark the part to copy".into(), palette.warn),
                Some(Status::Done) => {
                    let size = result.and_then(|r| r.outputs.first()).map(|p| format!("{} bp product", thousands(p.doc.len()))).unwrap_or_default();
                    match tails.as_slice() {
                        [] => (format!("{size} · choose the tail enzymes"), palette.warn),
                        list => (format!("{size} · {} tails", list.join(" + ")), palette.muted),
                    }
                }
                Some(Status::Failed) => (problem(), palette.error),
                _ => ("waits for the step above".into(), palette.faint),
            };
            ("COPY".into(), object, line, color)
        }
        Kind::Join => {
            // The other plasmids of the tube are counted, as the card may not have room for all.
            let others = match result.map_or(0, |r| r.findings.byproducts.len()) {
                0 => String::new(),
                n => format!(" · {n} more can form"),
            };
            let (line, color) = match status {
                Some(Status::Done) => {
                    let made = result.and_then(|r| r.outputs.first());
                    let shape = made.map_or("", |p| if p.doc.is_circular() { "circular" } else { "linear" });
                    (made.map(|p| format!("{} bp · {shape}{others}", thousands(p.doc.len()))).unwrap_or_default(), palette.muted)
                }
                Some(Status::NeedsChoice) => (format!("Choose the direction{others}"), palette.warn),
                Some(Status::Failed) => (problem(), palette.error),
                _ => ("waits for the pieces".into(), palette.faint),
            };
            ("JOIN".into(), "new plasmid".into(), line, color)
        }
        Kind::Other => {
            let step = c.step.map(|s| &tab.protocol.steps[s]);
            let title = step.map(|s| s.action.title().to_uppercase()).unwrap_or_default();
            let line = status.map(|s| s.name().to_string()).unwrap_or_default();
            (title, step.map(|s| s.id.clone()).unwrap_or_default(), line, status.map_or(palette.faint, |s| status_color(s, palette)))
        }
        Kind::Start => ("START".into(), "with a plasmid".into(), String::new(), palette.muted),
        Kind::PutInto => {
            // When a plasmid on the bench can already take the kept piece, the card says a drag
            // does the same as choosing here.
            let pick = board.cards.iter().find_map(|other| carried_by(tab, other));
            let droppable = pick.is_some_and(|pick| (0..board.cards.len()).any(|ii| drop_onto(tab, board, ii, pick).is_some()));
            let line = match (droppable, pick.is_some()) {
                (true, _) => "or drag the green piece onto a plasmid",
                (false, true) => "that the green piece goes into",
                (false, false) => "that the kept piece goes into",
            };
            ("PUT IT INTO".into(), "a plasmid".into(), line.into(), palette.faint)
        }
        Kind::AddPlasmid => ("ADD".into(), "another plasmid".into(), "to cut a piece from".into(), palette.faint),
    }
}

/// A kept piece carried across the bench by a drag: the pick step that keeps it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Carried {
    pick: usize,
}

/// The pick whose kept piece a card can carry to another card: the first pick of a cut card that
/// has kept a piece.
fn carried_by(tab: &ProtocolTab, card: &Card) -> Option<usize> {
    if card.kind != Kind::Cut {
        return None;
    }
    card.picks.iter().copied().find(|&p| tab.run.steps.get(p).is_some_and(|r| r.status == Status::Done))
}

/// What letting the piece of pick `pick` go onto card `card` does, and the words for it. A map
/// other than the piece's own takes it in, and so does a join that does not take it yet.
fn drop_onto(tab: &ProtocolTab, board: &Board, card: usize, pick: usize) -> Option<(Edit, String)> {
    let c = &board.cards[card];
    let step = c.step?;
    let pick_id = &tab.protocol.steps.get(pick)?.id;
    let piece = piece_words(tab, pick);
    match c.kind {
        Kind::Plasmid => {
            let own = tab.protocol.source_of(pick_id)?;
            (tab.protocol.steps[step].id != own).then(|| (Edit::PutPiece(pick, step), format!("Put {piece} into {}", map_name(tab, step))))
        }
        Kind::Join => {
            let takes = tab.protocol.steps[step].inputs.iter().any(|r| split_reference(r).0 == pick_id);
            (!takes).then(|| (Edit::AddToJoin(pick, step), format!("Add {piece} to the new plasmid")))
        }
        _ => None,
    }
}

/// A piece in a few words for the label that rides with the pointer: "EGFP · 944 bp".
fn carried_label(tab: &ProtocolTab, pick: usize) -> String {
    let size = tab.run.steps.get(pick).and_then(|r| r.outputs.first()).map(|p| format!("{} bp", thousands(p.doc.len())));
    let name = match tab.protocol.steps.get(pick).map(|s| &s.action) {
        Some(Action::Pick(params)) => params.contains.clone(),
        _ => None,
    };
    match (name, size) {
        (Some(name), Some(size)) => format!("{name} · {size}"),
        (None, Some(size)) => format!("piece · {size}"),
        (Some(name), None) => name,
        (None, None) => "piece".into(),
    }
}

/// What a join's input is, in words for a menu: the piece behind a treatment, as "the EGFP
/// piece" or "the backbone".
fn input_words(tab: &ProtocolTab, reference: &str) -> String {
    let mut id = split_reference(reference).0.to_string();
    // A treatment is named after what it treats, a few steps back at most.
    for _ in 0..8 {
        match tab.protocol.step(&id) {
            Some(Step { action: Action::Dephosphorylate(_) | Action::Flip(_), inputs, .. }) if !inputs.is_empty() => {
                id = split_reference(&inputs[0]).0.to_string()
            }
            _ => break,
        }
    }
    match tab.protocol.steps.iter().position(|s| s.id == id) {
        Some(index) if matches!(&tab.protocol.steps[index].action, Action::Pick(p) if p.contains.is_some()) => piece_words(tab, index),
        _ => format!("the {}", id.replace('_', " ")),
    }
}

/// The things a card's step can do: the card's "•••" button and its right-click both open these.
fn card_menu(ui: &mut egui::Ui, tab: &ProtocolTab, board: &Board, card: usize, docs: &[DocTab], palette: &Palette, edits: &mut Vec<Edit>) {
    let c = &board.cards[card];
    let Some(step) = c.step else { return };
    match c.kind {
        Kind::Plasmid => {
            if ui.button("Zoom into its map").clicked() {
                edits.push(Edit::Zoom(Some(step)));
                ui.close();
            }
            ui.menu_button("Use another map", |ui| {
                if let Some(source) = map_choices(ui, tab, docs, palette) {
                    edits.push(change_step(step, None, move |s| s.action = Action::Open(source)));
                    ui.close();
                }
            });
            if let Some(made) = tab.run.steps.get(step).and_then(|r| r.outputs.first())
                && ui.button("Open it as a tab").clicked()
            {
                edits.push(Edit::OpenDocument(Box::new(made.doc.clone())));
                ui.close();
            }
            lane_switch(ui, tab, board, card, edits);
            ui.separator();
            if ui.button("Take this plasmid out of the protocol").clicked() {
                edits.push(Edit::RemoveSource(step));
                ui.close();
            }
        }
        Kind::Cut => {
            let pieced = tab.run.steps.get(step).is_some_and(|r| !r.outputs.is_empty());
            if pieced && ui.button("Zoom into the pieces").clicked() {
                edits.push(Edit::Zoom(Some(step)));
                ui.close();
            }
            // The enzymes are chosen on the map of the plasmid the cut comes from, or for the cut
            // that trims a PCR's tails, on the PCR step.
            if let Some(plasmid) = board.plasmid_before(card).and_then(|p| board.cards[p].step)
                && ui.button("Choose the enzymes on its map").clicked()
            {
                edits.push(Edit::Zoom(Some(plasmid)));
                ui.close();
            }
            if let Some(copy) = board.before(card).filter(|&b| board.cards[b].kind == Kind::Copy).and_then(|b| board.cards[b].step)
                && ui.button("Choose the tail enzymes on the PCR step").clicked()
            {
                edits.push(Edit::Zoom(Some(copy)));
                ui.close();
            }
            if matches!(&tab.protocol.steps[step].action, Action::Digest(d) if !d.enzymes.is_empty()) && ui.button("Clear the enzymes").clicked() {
                edits.push(Edit::SetEnzymes(step, Vec::new()));
                ui.close();
            }
            // Where the kept piece can go, for those who would rather not drag it.
            if let Some(pick) = carried_by(tab, c) {
                let places: Vec<usize> = (0..board.cards.len()).filter(|&ii| drop_onto(tab, board, ii, pick).is_some()).collect();
                if !places.is_empty() {
                    ui.separator();
                    ui.menu_button("Put the kept piece into", |ui| {
                        for ii in places {
                            if let Some((edit, words)) = drop_onto(tab, board, ii, pick)
                                && ui.button(words).clicked()
                            {
                                edits.push(edit);
                                ui.close();
                            }
                        }
                    });
                }
            }
            lane_switch(ui, tab, board, card, edits);
        }
        Kind::Copy => {
            if ui.button("Zoom into the map").on_hover_text("The map the part is copied from, with the primers drawn on it").clicked() {
                edits.push(Edit::Zoom(Some(step)));
                ui.close();
            }
            let (pinned, tails) = match &tab.protocol.steps[step].action {
                Action::Amplify(a) => (a.pinned(), !a.enzymes().is_empty()),
                _ => (false, false),
            };
            let designed = tab.run.steps.get(step).is_some_and(|r| !r.findings.primer_pairs.is_empty());
            if designed && ui.button("Keep these primers").on_hover_text("Write the primers into the protocol, so they stay as they are").clicked() {
                edits.push(Edit::Primers(step, PrimerEdit::Pair(0)));
                ui.close();
            }
            if pinned && ui.button("Design the primers anew").clicked() {
                edits.push(Edit::Primers(step, PrimerEdit::Clear));
                ui.close();
            }
            if tails && ui.button("Clear the tail enzymes").clicked() {
                edits.push(Edit::SetEnzymes(step, Vec::new()));
                ui.close();
            }
            if let Some(made) = tab.run.steps.get(step).and_then(|r| r.outputs.first())
                && ui.button("Open the PCR product as a tab").clicked()
            {
                edits.push(Edit::OpenDocument(Box::new(made.doc.clone())));
                ui.close();
            }
            lane_switch(ui, tab, board, card, edits);
        }
        Kind::Join => {
            if ui.button("Zoom into the new plasmid").clicked() {
                edits.push(Edit::Zoom(Some(step)));
                ui.close();
            }
            if let Some(made) = tab.run.steps.get(step).and_then(|r| r.outputs.first())
                && ui.button("Open the new plasmid as a tab").clicked()
            {
                let mut doc = made.doc.clone();
                doc.name = tab.name();
                edits.push(Edit::OpenDocument(Box::new(doc)));
                ui.close();
            }
            // Each piece can come out, and can lose or keep its 5' phosphates on the way in.
            let inputs = tab.protocol.steps[step].inputs.clone();
            for (ii, reference) in inputs.iter().enumerate() {
                ui.separator();
                let name = input_words(tab, reference);
                let treated = tab
                    .protocol
                    .steps
                    .iter()
                    .position(|s| s.id == split_reference(reference).0)
                    .filter(|&index| matches!(tab.protocol.steps[index].action, Action::Dephosphorylate(_)));
                match treated {
                    Some(index) => {
                        if ui.button(format!("Keep the 5' phosphates of {name}")).clicked() {
                            edits.push(Edit::LeaveOut(index));
                            ui.close();
                        }
                    }
                    None => {
                        if ui.button(format!("Remove the 5' phosphates of {name}")).clicked() {
                            edits.push(Edit::Dephosphorylate(step, ii));
                            ui.close();
                        }
                    }
                }
                if ui.button(format!("Take {name} out")).clicked() {
                    edits.push(Edit::TakeOut(step, ii));
                    ui.close();
                }
            }
        }
        Kind::Other => {
            if ui.button("Show this step").clicked() {
                edits.push(Edit::Select(step));
                ui.close();
            }
            if ui.button("Delete this step").clicked() {
                edits.push(Edit::RemoveStep(step));
                ui.close();
            }
        }
        Kind::Start | Kind::PutInto | Kind::AddPlasmid => {}
    }
}

/// The menu items that switch a lane between cutting its part out and copying it by PCR, on the
/// plasmid, PCR or cut card of the lane. A vector's cut takes a part out for the insert, so its
/// lane offers neither.
fn lane_switch(ui: &mut egui::Ui, tab: &ProtocolTab, board: &Board, card: usize, edits: &mut Vec<Edit>) {
    let Some(step) = board.cards[card].step else { return };
    let id = tab.protocol.steps[step].id.as_str();
    // The lane's cut, from whichever of its cards the menu is on.
    let cut: Option<String> = match board.cards[card].kind {
        Kind::Plasmid => board.part_step_after(card).and_then(|s| match &tab.protocol.steps[s].action {
            Action::Digest(_) => Some(tab.protocol.steps[s].id.clone()),
            Action::Amplify(_) => tab.protocol.digest_after(&tab.protocol.steps[s].id).map(str::to_string),
            _ => None,
        }),
        Kind::Copy => tab.protocol.digest_after(id).map(str::to_string),
        Kind::Cut => Some(id.to_string()),
        _ => None,
    };
    let Some(cut) = cut else { return };
    ui.separator();
    if tab.protocol.amplify_before(&cut).is_some() {
        if ui
            .button("Cut the part out instead of copying it")
            .on_hover_text("Take the PCR out: the enzymes then cut the part straight out of the map")
            .clicked()
        {
            edits.push(Edit::CutOut(step));
            ui.close();
        }
    } else if protocol::part_role(&tab.protocol, &cut) != Some(PartRole::Replace)
        && ui
            .button("Copy the part out by PCR instead")
            .on_hover_text("Copy the part with primers that carry the enzyme sites on their tails, then cut the tails off")
            .clicked()
    {
        edits.push(Edit::CopyOut(step));
        ui.close();
    }
}

/// Draw one card: its frame, its title, its ring (or a dashed place for one) and its status
/// line. A click zooms into it, a cut's kept piece can be dragged onto another card, and the
/// card's "•••" button (or a right-click) lists what else it can do. Returns where its ring was
/// drawn.
#[allow(clippy::too_many_arguments)]
fn draw_card(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    board: &Board,
    card: usize,
    rect: Rect,
    scale: f32,
    docs: &[DocTab],
    settings: &ViewSettings,
    palette: &Palette,
    edits: &mut Vec<Edit>,
    interactive: bool,
) -> Option<Rect> {
    let c = &board.cards[card];
    let selected = c.step.is_some() && board.card_of(tab.selected) == Some(card);
    let (verb, object, line, line_color) = card_text(tab, board, card, palette);
    let status = c.step.and_then(|s| tab.run.steps.get(s)).map(|r| r.status);
    // The card under the pointer is outlined, so it is plain that it does something.
    let hovered = interactive && c.step.is_some() && ui.rect_contains_pointer(rect) && DragAndDrop::payload::<Carried>(ui.ctx()).is_none();
    let edge = match (selected, status) {
        (true, _) => Stroke::new(2.0, palette.accent),
        (false, Some(Status::NeedsChoice)) => Stroke::new(1.5, palette.warn),
        (false, Some(Status::Failed)) => Stroke::new(1.5, palette.error),
        _ if hovered => Stroke::new(1.5, palette.accent.gamma_multiply(0.6)),
        _ => Stroke::new(1.0, palette.line),
    };
    let painter = ui.painter();
    painter.rect(rect, CornerRadius::same(10), palette.raised, edge, egui::StrokeKind::Inside);

    // Title: the verb in small capitals, then what it acts on.
    let title_y = rect.min.y + CARD_TITLE * scale * 0.55;
    let verb_font = style::bold_font(10.5 * scale.min(1.2));
    let verb_rect = painter.text(pos2(rect.min.x + 12.0 * scale, title_y), Align2::LEFT_CENTER, &verb, verb_font, palette.accent);
    let room = rect.max.x - verb_rect.max.x - 62.0 * scale;
    let object_font = style::ui_font(12.5 * scale.min(1.2));
    painter.text(pos2(verb_rect.max.x + 7.0, title_y), Align2::LEFT_CENTER, fit(ui, &object, object_font.clone(), room), object_font, palette.ink);
    if let Some(status) = status {
        painter.circle_filled(pos2(rect.max.x - 13.0 * scale, title_y), 4.0, status_color(status, palette));
    }

    // The ring, or a dashed place for one with what it waits for. A cut's pieces pull up and
    // out, so its ring sits a little lower to leave room for the top piece's label.
    let drop = if c.kind == Kind::Cut { 12.0 * scale } else { 0.0 };
    let ring_rect = Rect::from_min_max(pos2(rect.min.x, rect.min.y + CARD_TITLE * scale + drop), pos2(rect.max.x, rect.max.y - CARD_STATUS * scale));
    let ctx = ui.ctx().clone();
    let measure = EguiMeasure { ctx: &ctx };
    // The ring as drawn, and the roles of its pieces, to tell what the pointer is over.
    let mut drawn: Option<(Scene, Vec<Role>)> = None;
    // A join whose tube makes several plasmids shows them side by side. The first one's ring is
    // where a zoom starts from.
    let made = c.step.filter(|_| c.kind == Kind::Join).and_then(|s| tab.run.steps.get(s)).map(join_products).unwrap_or_default();
    let mut first_ring = ring_rect;
    let mut places: Vec<(Rect, Made)> = Vec::new();
    match picture(tab, board, card, false, settings, palette.raised) {
        Some(_) if made.len() > 1 => {
            let step = c.step.unwrap_or_default();
            (first_ring, places) = draw_made(ui, tab, step, &made, rect, ring_rect, scale, settings, palette);
        }
        Some(shown) => {
            let mut options = shown.options;
            options.width = ring_rect.width();
            options.height = ring_rect.height();
            options.font_size = 10.5 * scale.min(1.2);
            let scene = ring::render_ring(&shown.doc, &options, &measure);
            paint(&ui.painter_at(rect), &scene, ring_rect.min);
            drawn = Some((scene, options.pieces.iter().map(|p| p.role).collect()));
        }
        None => {
            let radius = ring::radius(ring_rect.width(), ring_rect.height());
            let points: Vec<Pos2> = (0..=96).map(|ii| ring_rect.center() + radius * egui::Vec2::angled(ii as f32 / 96.0 * std::f32::consts::TAU)).collect();
            ui.painter().extend(egui::Shape::dashed_line(&points, Stroke::new(1.5, palette.faint), 6.0, 5.0));
            // What the empty place waits for, or why its step failed, inside it.
            let waiting = match (c.kind, status) {
                (Kind::Cut, Some(Status::NeedsChoice)) => "click to choose\nwhere to cut",
                (Kind::Copy, Some(Status::NeedsChoice)) => "click to mark\nthe part to copy",
                (Kind::Join | Kind::Cut | Kind::Copy, Some(Status::Failed)) => "see the panel\nfor why",
                (Kind::Join | Kind::Cut | Kind::Copy, _) => "waits for\nthe steps above",
                _ => "",
            };
            // Each line sits centred, the block in the middle of the circle.
            let lines: Vec<&str> = waiting.lines().collect();
            let font = style::ui_font(11.0 * scale.min(1.2));
            let pitch = font.size * 1.3;
            for (ii, line) in lines.iter().enumerate() {
                let y = ring_rect.center().y + (ii as f32 - (lines.len() as f32 - 1.0) / 2.0) * pitch;
                ui.painter().text(pos2(ring_rect.center().x, y), Align2::CENTER_CENTER, *line, font.clone(), palette.faint);
            }
        }
    }
    // Places to start, to put the insert into a vector, or to add a map hold their choices: a
    // file, or a map open in the app.
    if matches!(c.kind, Kind::Start | Kind::PutInto | Kind::AddPlasmid) && interactive {
        let place = Rect::from_center_size(ring_rect.center(), vec2(rect.width() * 0.7, 30.0));
        let mut child = ui.new_child(UiBuilder::new().max_rect(place).id_salt(("protocol-place", tab.id, card)));
        child.vertical_centered(|ui| {
            ui.menu_button("Choose a plasmid…", |ui| {
                if let Some(source) = map_choices(ui, tab, docs, palette) {
                    edits.push(if c.kind == Kind::PutInto { Edit::PutInto(source) } else { Edit::AddSource(source) });
                    ui.close();
                }
            })
            .response
            .on_hover_text(match c.kind {
                Kind::PutInto => "The plasmid the kept piece goes into. It is cut with the same enzymes.",
                Kind::AddPlasmid => "Another plasmid to cut a piece from. Drag that piece into the new plasmid afterwards.",
                _ => "The plasmid that carries the part you want to cut out",
            });
        });
    }

    let line_font = style::ui_font(11.5 * scale.min(1.2));
    let line_y = rect.max.y - CARD_STATUS * scale * 0.5;
    let text = fit(ui, &line, line_font.clone(), rect.width() - 20.0 * scale);
    ui.painter().text(pos2(rect.center().x, line_y), Align2::CENTER_CENTER, text, line_font, line_color);

    let Some(step) = c.step.filter(|_| interactive) else { return Some(first_ring) };
    // A click zooms in. A cut that kept a piece can also be dragged, which carries that piece.
    let carries = carried_by(tab, c);
    let sense = if carries.is_some() { Sense::click_and_drag() } else { Sense::click() };
    let response = ui.interact(rect, egui::Id::unique(("protocol-card", tab.id, card)), sense);
    let carrying = DragAndDrop::payload::<Carried>(ui.ctx());
    if let Some(pick) = carries
        && response.drag_started()
    {
        DragAndDrop::set_payload(ui.ctx(), Carried { pick });
    }
    if response.dragged() && carrying.is_some() {
        ui.ctx().set_cursor_icon(CursorIcon::Grabbing);
    }
    // A piece let go over this card goes into it, when it can.
    if let Some(carried) = response.dnd_release_payload::<Carried>()
        && let Some((edit, _)) = drop_onto(tab, board, card, carried.pick)
    {
        edits.push(edit);
    }
    // A cut without pieces yet has nothing to show but its plasmid's map, so a click on it zooms
    // into the card before it (the plasmid, or the PCR whose product it trims), where the part
    // and the cuts are chosen. The map stays put while the cuts are clicked one by one.
    if response.clicked() {
        let unpieced = c.kind == Kind::Cut && status.is_some_and(|_| tab.run.steps[step].outputs.is_empty());
        let target = if unpieced { board.before(card).and_then(|p| board.cards[p].step).unwrap_or(step) } else { step };
        edits.push(Edit::Zoom(Some(target)));
    }
    // Over the kept piece the cursor offers to grab it; elsewhere it offers to zoom.
    if hovered {
        let over_kept = drawn.as_ref().zip(response.hover_pos()).is_some_and(|((scene, roles), p)| {
            matches!(hit(scene, [p.x - ring_rect.min.x, p.y - ring_rect.min.y], &measure), Some(Tag::Piece(ii)) if roles.get(ii) == Some(&Role::Kept))
        });
        ui.ctx().set_cursor_icon(if over_kept && carries.is_some() { CursorIcon::Grab } else { CursorIcon::ZoomIn });
    }
    response.context_menu(|ui| card_menu(ui, tab, board, card, docs, palette, edits));
    // Over one of a join's plasmids, the hint says what that plasmid is.
    let over = response.hover_pos().and_then(|p| places.iter().find(|(place, _)| place.contains(p))).map(|(_, made)| *made);
    let hint = match (c.kind, carries.is_some(), over) {
        (Kind::Join, _, Some(Made::Other(ii))) => {
            let words =
                tab.run.steps.get(step).and_then(|r| r.findings.byproducts.get(ii)).map(|b| protocol::byproduct_words(&tab.protocol, &tab.run, step, b));
            format!("The same tube also makes {}. Click to zoom in and compare.", words.unwrap_or_default())
        }
        (Kind::Join, _, Some(Made::Option(ii))) => format!("Option {}: the insert fits either way round. Click to zoom in and choose.", ii + 1),
        (Kind::Plasmid, _, _) if board.copy_after(card).is_some() => {
            "Click to zoom into the map: mark the part to copy there, and choose the enzymes for the primers' tails in the panel.".into()
        }
        (Kind::Plasmid, _, _) => "Click to zoom into the map: mark the part there, and click the cuts.".into(),
        (Kind::Copy, _, _) => {
            "Click to zoom into the map, where the primers are drawn: mark the part to copy there, and choose the tail enzymes in the panel.".into()
        }
        (Kind::Cut, false, _) if status == Some(Status::NeedsChoice) => "Click to zoom into the plasmid's map: mark the part there, and click the cuts.".into(),
        (Kind::Cut, true, _) => "Click to zoom into the pieces. Drag the green piece onto another plasmid, or onto the new plasmid, to put it in.".into(),
        (Kind::Cut, false, _) => "Click to zoom into the pieces.".into(),
        (Kind::Join, _, _) => "Click to zoom into the new plasmid. Drop a green piece here to add it.".into(),
        _ => "Click to show this step.".to_string(),
    };
    response.on_hover_text(hint);

    // The "•••" button opens the same menu as a right-click, so the menu can be found. It sits
    // left of the status dot.
    // The dots stay on one line: the button grows to fit them rather than wrapping them.
    let dots = Rect::from_center_size(pos2(rect.max.x - 38.0 * scale, title_y), vec2(28.0, 20.0));
    let mut child = ui.new_child(UiBuilder::new().max_rect(dots).id_salt(("protocol-card-menu", tab.id, card)));
    child.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
    child.spacing_mut().button_padding = vec2(5.0, 1.0);
    child
        .menu_button(RichText::new("•••").size(10.0).color(palette.muted), |ui| card_menu(ui, tab, board, card, docs, palette, edits))
        .response
        .on_hover_text("What you can do with this step");
    Some(first_ring)
}

/// Draw the plasmids of join step `step` side by side on its card `rect`, within `ring_rect`, the
/// card's room for rings. The product or the options take full-size places. The other plasmids of
/// the tube follow in smaller places, past a thin line, each over the words "also forms" and how it
/// differs. Returns where the first ring was drawn, for the zoom, and each place with the plasmid
/// it shows, for the pointer.
#[allow(clippy::too_many_arguments)]
fn draw_made(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    step: usize,
    made: &[Shown],
    rect: Rect,
    ring_rect: Rect,
    scale: f32,
    settings: &ViewSettings,
    palette: &Palette,
) -> (Rect, Vec<(Rect, Made)>) {
    let result = tab.run.steps.get(step);
    let (main, others) = join_slots(result);
    let shown: Vec<&Shown> = made.iter().take(main + others).collect();
    // The places sit side by side in the middle of the card.
    let width = |m: &Shown| (if matches!(m.made, Made::Other(_)) { OTHER_SLOT } else { MAIN_SLOT }) * scale;
    let total: f32 = shown.iter().map(|m| width(m)).sum();
    let mut x = rect.center().x - total / 2.0;
    let ctx = ui.ctx().clone();
    let measure = EguiMeasure { ctx: &ctx };
    let painter = ui.painter_at(rect);
    let caption_font = style::ui_font(10.5 * scale.min(1.2));
    let line = caption_font.size * 1.3;
    let mut first = ring_rect;
    let mut places = Vec::new();
    for (ii, m) in shown.iter().enumerate() {
        let place = Rect::from_min_max(pos2(x, ring_rect.min.y), pos2(x + width(m), ring_rect.max.y));
        x += width(m);
        let other = matches!(m.made, Made::Other(_));
        // A thin line parts the plan's plasmids from the others of the tube.
        if other && ii > 0 && !matches!(shown[ii - 1].made, Made::Other(_)) {
            painter.line_segment([pos2(place.min.x, place.min.y + 8.0 * scale), pos2(place.min.x, place.max.y - 4.0 * scale)], Stroke::new(1.0, palette.line));
        }
        let caption = made_caption(tab, step, m.made);
        let lines = if other {
            2.0
        } else if caption.is_empty() {
            0.0
        } else {
            1.0
        };
        let room = Rect::from_min_max(place.min, pos2(place.max.x, place.max.y - lines * line));
        // The others are drawn smaller, so the plan's plasmids lead the eye.
        let canvas = if other {
            let side = room.width().min(room.height()) * 0.82;
            Rect::from_center_size(pos2(room.center().x, room.max.y - side / 2.0), egui::Vec2::splat(side))
        } else {
            room
        };
        let mut options = made_ring(m, settings, palette.raised);
        options.width = canvas.width();
        options.height = canvas.height();
        options.font_size = (if other { 9.5 } else { 10.5 }) * scale.min(1.2);
        let scene = ring::render_ring(&m.piece.doc, &options, &measure);
        paint(&painter, &scene, canvas.min);
        if ii == 0 {
            first = canvas;
        }
        // The words under the ring: an option's number and direction, or what else forms.
        let text_width = place.width() - 8.0 * scale;
        let mut y = place.max.y - lines * line + line / 2.0;
        if other {
            painter.text(pos2(place.center().x, y), Align2::CENTER_CENTER, "also forms", caption_font.clone(), palette.faint);
            y += line;
        }
        if !caption.is_empty() {
            let color = if other { palette.muted } else { palette.warn };
            painter.text(pos2(place.center().x, y), Align2::CENTER_CENTER, fit(ui, &caption, caption_font.clone(), text_width), caption_font.clone(), color);
        }
        places.push((place, m.made));
    }
    (first, places)
}

/// While a kept piece is carried, every card that takes it lights up, the one under the pointer
/// most, with a label saying what letting go there does. The piece rides beside the pointer.
fn draw_carrying(ui: &mut egui::Ui, tab: &ProtocolTab, board: &Board, card_rect: impl Fn(&Card) -> Rect, palette: &Palette) {
    let Some(carried) = DragAndDrop::payload::<Carried>(ui.ctx()) else { return };
    let pointer = ui.ctx().pointer_latest_pos();
    let font = style::ui_font(12.0);
    for (ii, card) in board.cards.iter().enumerate() {
        let Some((_, words)) = drop_onto(tab, board, ii, carried.pick) else { continue };
        let rect = card_rect(card);
        let over = pointer.is_some_and(|p| rect.contains(p));
        let painter = ui.painter();
        if over {
            painter.rect_filled(rect, CornerRadius::same(10), palette.accent_soft);
        }
        let stroke = if over { Stroke::new(3.0, palette.accent) } else { Stroke::new(1.5, palette.accent.gamma_multiply(0.7)) };
        painter.rect_stroke(rect.expand(3.0), CornerRadius::same(12), stroke, egui::StrokeKind::Outside);
        let caption = if over { format!("Let go to {}", lowercase_first(&words)) } else { words };
        let galley = painter.layout_no_wrap(caption, font.clone(), palette.raised);
        let label = Rect::from_center_size(pos2(rect.center().x, rect.min.y - 14.0), galley.size() + vec2(16.0, 8.0));
        painter.rect_filled(label, CornerRadius::same(8), palette.accent);
        painter.galley(label.center() - galley.size() / 2.0, galley, palette.raised);
    }
    // The piece itself, as a green bar and its name, on top of everything.
    if let Some(p) = pointer {
        let painter = ui.ctx().layer_painter(LayerId::new(Order::Tooltip, egui::Id::unique(("protocol-carry", tab.id))));
        let galley = painter.layout_no_wrap(carried_label(tab, carried.pick), style::bold_font(12.5), palette.ink);
        let body = Rect::from_min_size(p + vec2(14.0, 10.0), galley.size() + vec2(42.0, 12.0));
        painter.rect(body, CornerRadius::same(8), palette.raised, Stroke::new(1.5, palette.ok), egui::StrokeKind::Inside);
        let bar = Rect::from_center_size(pos2(body.min.x + 17.0, body.center().y), vec2(18.0, 6.0));
        painter.rect_filled(
            bar,
            CornerRadius::same(3),
            scene_paint::color(theme_on(&ViewSettings { dark: palette.dark, ..ViewSettings::default() }, palette.raised).kept),
        );
        painter.galley(pos2(body.min.x + 32.0, body.center().y - galley.size().y / 2.0), galley, palette.ink);
    }
}

/// Text with its first letter in lower case, to follow other words in a sentence.
fn lowercase_first(text: &str) -> String {
    let mut chars = text.chars();
    match chars.next() {
        Some(first) => first.to_lowercase().collect::<String>() + chars.as_str(),
        None => String::new(),
    }
}

/// Text with its first letter in upper case, to start a sentence.
fn uppercase_first(text: &str) -> String {
    let mut chars = text.chars();
    match chars.next() {
        Some(first) => first.to_uppercase().collect::<String>() + chars.as_str(),
        None => String::new(),
    }
}

/// Draw an arrow between two cards: straight down within a column, curved into a join. A pill
/// on it says what happens on the way: the enzymes of a cut, a dephosphorylation, a flip.
#[allow(clippy::too_many_arguments)]
fn draw_link(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    board: &Board,
    link: &Link,
    from: Rect,
    to: Rect,
    slot: (usize, usize),
    scale: f32,
    palette: &Palette,
    edits: &mut Vec<Edit>,
    interactive: bool,
) {
    let start = from.center_bottom();
    // Several arrows into one card spread across its top edge.
    let end = pos2(to.min.x + to.width() * (slot.0 + 1) as f32 / (slot.1 + 1) as f32, to.min.y);
    let stroke = Stroke::new(1.6, palette.faint);
    let painter = ui.painter();
    let tip = pos2(end.x, end.y - 2.0);
    if (start.x - end.x).abs() < 1.0 {
        painter.line_segment([start, tip], stroke);
    } else {
        let bend = (end.y - start.y) * 0.5;
        let curve = egui::epaint::CubicBezierShape::from_points_stroke(
            [start, start + vec2(0.0, bend), tip - vec2(0.0, bend), tip],
            false,
            Color32::TRANSPARENT,
            stroke,
        );
        painter.add(curve);
    }
    let head = 6.0;
    painter.add(egui::Shape::convex_polygon(vec![tip, tip + vec2(-head * 0.7, -head), tip + vec2(head * 0.7, -head)], palette.faint, Stroke::NONE));

    // What happens on the way, as pills at the middle of the arrow. A cut's pill leads to the
    // plasmid's map, where its enzymes are chosen. A treatment's pill shows the treatment, and its
    // right-click leaves it out.
    enum Pill {
        Cut(usize),
        Copy(usize),
        Treatment(usize),
    }
    let mut pills: Vec<(String, Color32, Pill)> = Vec::new();
    if board.cards[link.to].kind == Kind::Cut
        && let Some(Action::Digest(params)) = board.cards[link.to].step.map(|s| &tab.protocol.steps[s].action)
        && let Some(from) = board.cards[link.from].step
    {
        // The cut after a PCR trims the tails, and its enzymes are chosen on the PCR step.
        let trims = board.cards[link.from].kind == Kind::Copy;
        let text = match (params.enzymes.is_empty(), trims) {
            (true, true) => "trim the tails: choose enzymes".to_string(),
            (true, false) => "cut: choose enzymes".to_string(),
            (false, true) => format!("trim the tails with {}", params.enzymes.join(" + ")),
            (false, false) => format!("cut with {}", params.enzymes.join(" + ")),
        };
        let color = if params.enzymes.is_empty() { palette.warn } else { palette.ink };
        pills.push((text, color, if trims { Pill::Copy(from) } else { Pill::Cut(from) }));
    }
    if board.cards[link.to].kind == Kind::Copy
        && let Some(copy) = board.cards[link.to].step
        && let Action::Amplify(params) = &tab.protocol.steps[copy].action
    {
        let tails = params.enzymes();
        let text = if tails.is_empty() { "PCR: choose the tail enzymes".to_string() } else { format!("PCR with {} tails", tails.join(" + ")) };
        pills.push((text, if tails.is_empty() { palette.warn } else { palette.ink }, Pill::Copy(copy)));
    }
    for &step in &link.treatments {
        let text = match tab.protocol.steps[step].action {
            Action::Dephosphorylate(_) => "no 5' phosphates",
            Action::Flip(_) => "turned round",
            _ => continue,
        };
        pills.push((text.to_string(), palette.muted, Pill::Treatment(step)));
    }
    let middle = if (start.x - end.x).abs() < 1.0 { pos2(start.x, (start.y + end.y) / 2.0) } else { pos2((start.x + end.x) / 2.0, (start.y + end.y) / 2.0) };
    let font = style::ui_font(11.0 * scale.min(1.2));
    for (ii, (text, color, what)) in pills.iter().enumerate() {
        let at = middle + vec2(0.0, (ii as f32 - (pills.len() as f32 - 1.0) / 2.0) * 22.0 * scale.min(1.2));
        let galley = ui.painter().layout_no_wrap(text.clone(), font.clone(), *color);
        let pill = Rect::from_center_size(at, galley.size() + vec2(14.0, 6.0));
        let id = egui::Id::unique(("protocol-pill", tab.id, link.from, link.to, ii));
        let response = if interactive { Some(ui.interact(pill, id, Sense::click())) } else { None };
        let hot = response.as_ref().is_some_and(|r| r.hovered());
        let edge = if hot { Stroke::new(1.5, palette.accent) } else { Stroke::new(1.0, palette.line) };
        ui.painter().rect(pill, CornerRadius::same(10), palette.panel, edge, egui::StrokeKind::Inside);
        ui.painter().galley(pill.center() - galley.size() / 2.0, galley, *color);
        let Some(response) = response else { continue };
        if hot {
            ui.ctx().set_cursor_icon(CursorIcon::PointingHand);
        }
        match what {
            Pill::Cut(plasmid) => {
                if response.clicked() {
                    edits.push(Edit::Zoom(Some(*plasmid)));
                }
                response.on_hover_text("Click to choose the enzymes on the plasmid's map");
            }
            Pill::Copy(copy) => {
                if response.clicked() {
                    edits.push(Edit::Zoom(Some(*copy)));
                }
                response.on_hover_text("Click to see the primers on the map and choose the enzymes for their tails");
            }
            Pill::Treatment(step) => {
                if response.clicked() {
                    edits.push(Edit::Select(*step));
                }
                response.context_menu(|ui| {
                    if ui.button("Leave this step out").clicked() {
                        edits.push(Edit::LeaveOut(*step));
                        ui.close();
                    }
                });
                let told = protocol::story(&tab.protocol, &tab.run, *step).details.join(" ");
                response.on_hover_text(format!("{told}\nClick to show it in the panel. Right-click to leave it out."));
            }
        }
    }
}

/// The molecule a zoomed card's sequence pane shows: a plasmid with the part marked on it, a cut
/// molecule with the piece kept, or a join's product with the insert.
struct PreviewTarget<'a> {
    /// The step whose molecule it is.
    step: usize,
    doc: &'a Document,
    /// The stretch highlighted: the part, the piece kept, or the insert.
    highlight: Option<Span>,
    /// The cut or PCR whose part a stretch selected in the pane becomes, when there is one.
    setting: Option<usize>,
}

/// The sequence shown under a zoomed card: a read-only document of the molecule, kept from frame
/// to frame so it scrolls and selects like any sequence view.
struct Preview {
    /// The step whose molecule it shows.
    step: usize,
    tab: DocTab,
    /// The highlight the picture showed when the selection was last set from it.
    shown: Option<Span>,
    /// A part set from the pane that the run has not taken up yet: the highlight as it was, and
    /// when the part was set. The selection is not pulled back to the old part meanwhile.
    pending: Option<(Option<Span>, Instant)>,
}

/// How long a part set from the sequence pane is waited for before the picture leads again.
const PENDING_FOR: Duration = Duration::from_secs(1);

/// What the sequence pane of card `card`, zoomed into, shows, or `None` for a card without a
/// molecule to show.
fn preview_target<'a>(tab: &'a ProtocolTab, board: &Board, card: usize) -> Option<PreviewTarget<'a>> {
    let c = &board.cards[card];
    match c.kind {
        Kind::Plasmid | Kind::Copy => {
            let (map_step, setting) = map_of(board, card)?;
            let doc = &tab.run.steps.get(map_step)?.outputs.first()?.doc;
            let highlight = setting.and_then(|s| tab.run.steps.get(s)).and_then(|r| r.findings.region);
            Some(PreviewTarget { step: map_step, doc, highlight, setting })
        }
        // A cut without pieces yet shows the map of the card before it, and so does its pane.
        Kind::Cut if tab.run.steps.get(c.step?).is_none_or(|r| r.outputs.is_empty()) => preview_target(tab, board, board.before(card)?),
        Kind::Cut => {
            let step = c.step?;
            let doc = &tab.run.steps.get(step)?.inputs.first()?.doc;
            let view = cut_view(tab, c)?;
            let highlight = view.pieces.iter().find(|p| p.1 == Role::Kept).map(|p| p.0).or(view.region);
            Some(PreviewTarget { step, doc, highlight, setting: None })
        }
        Kind::Join => {
            let step = c.step?;
            let result = tab.run.steps.get(step).filter(|r| r.status == Status::Done)?;
            let doc = &result.outputs.first()?.doc;
            let highlight = result.findings.parts.last().copied().filter(|_| result.findings.parts.len() > 1);
            Some(PreviewTarget { step, doc, highlight, setting: None })
        }
        _ => None,
    }
}

/// The gap between a zoomed card's picture and the sequence pane under it, which holds the
/// handle that moves the split.
const SPLIT_GAP: f32 = 10.0;

/// A zoomed card's area under its header, as the picture and, when the view shows the molecule's
/// bases, the sequence pane under it, which takes `sequence` of the height.
fn detail_rects(area: Rect, sequence: Option<f32>) -> (Rect, Option<Rect>) {
    let top = area.min.y + DETAIL_HEADER;
    let Some(share) = sequence else { return (detail_picture_rect(area), None) };
    let split = (area.max.y - (area.max.y - top) * share).round();
    let picture = Rect::from_min_max(pos2(area.min.x + 12.0, top), pos2(area.max.x - 12.0, split - SPLIT_GAP));
    (picture, Some(Rect::from_min_max(pos2(area.min.x, split), area.max)))
}

/// The share of a zoomed card's height that the sequence pane takes, when the view shows one.
fn sequence_share(tab: &ProtocolTab, board: &Board, card: usize) -> Option<f32> {
    (tab.show_sequence && preview_target(tab, board, card).is_some()).then_some(tab.sequence_share)
}

/// Where the picture of card `card`, zoomed into, is drawn within `area`.
fn picture_rect(tab: &ProtocolTab, board: &Board, card: usize, area: Rect) -> Rect {
    detail_rects(area, sequence_share(tab, board, card)).0
}

/// The stretch being dragged across the map of card `card` (or of the card before a cut) while
/// the button is held, as [`marking_map`] keeps it.
fn map_drag_span(ui: &egui::Ui, tab: &ProtocolTab, board: &Board, card: usize) -> Option<Span> {
    let card = if board.cards[card].kind == Kind::Cut { board.before(card)? } else { card };
    let (step, _) = map_of(board, card)?;
    let marking: Option<(usize, usize)> = ui.data(|d| d.get_temp(egui::Id::unique(("protocol-map-drag", tab.id, step))));
    let doc = &tab.run.steps.get(step)?.outputs.first()?.doc;
    marking.and_then(|(a, b)| map_view::drag_span(a, b, doc.len(), doc.is_circular()))
}

/// The sequence under a zoomed card's picture: a read-only document of the molecule with the
/// highlight as its selection. The selection follows the picture (the part as it is, or the
/// stretch being dragged across the map, `map_live`) and leads it: a stretch selected here shows
/// on the picture, and once the button is up it becomes the part when the picture sets one.
/// Returns the selection made here while it differs from the picture's highlight.
#[allow(clippy::too_many_arguments)]
fn sequence_pane(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    preview: &mut Option<Preview>,
    target: &PreviewTarget,
    pane: Rect,
    map_live: Option<Span>,
    settings: &ViewSettings,
    palette: &Palette,
    edits: &mut Vec<Edit>,
    interactive: bool,
) -> Option<Span> {
    // The document shown: made for the molecule, and read again when the molecule changes. It
    // keeps its identity across changes, so its scroll position stays.
    if preview.as_ref().is_none_or(|p| p.step != target.step) {
        let mut doc_tab = DocTab::new(target.doc.clone(), None);
        doc_tab.read_only = true;
        *preview = Some(Preview { step: target.step, tab: doc_tab, shown: None, pending: None });
    }
    let preview = preview.as_mut()?;
    if preview.tab.doc != *target.doc {
        preview.tab.doc = target.doc.clone();
        preview.tab.revision += 1;
        preview.tab.select(None);
        preview.shown = None;
    }
    if preview.tab.enzyme_set != tab.set {
        preview.tab.enzyme_set = tab.set.clone();
    }
    // What the picture highlights moves the selection here when it changed. A part just set from
    // here is waited for, so the old part does not pull the selection back meanwhile.
    if let Some((old, since)) = preview.pending
        && (target.highlight != old || since.elapsed() > PENDING_FOR)
    {
        preview.pending = None;
    }
    let wanted = map_live.or(target.highlight);
    if preview.pending.is_none() && wanted != preview.shown {
        preview.tab.select(wanted);
        if let Some(span) = wanted {
            preview.tab.scroll_to = Some(span.start);
        }
        preview.shown = wanted;
    }
    // The pane: a line that says what it shows, then the sequence view.
    ui.painter().line_segment([pane.left_top(), pane.right_top()], Stroke::new(1.0, palette.line));
    let inner = Rect::from_min_max(pane.min + vec2(0.0, 1.0), pane.max);
    let mut child = ui.new_child(UiBuilder::new().max_rect(inner).id_salt(("protocol-sequence", tab.id)));
    // Nothing of the pane may reach past its edge into the panel beside it.
    child.set_clip_rect(inner.intersect(child.clip_rect()));
    child.add_enabled_ui(interactive, |ui| {
        // The line: Hide at the right, placed first so that the words cannot push it past the
        // edge, then the name and the hint, which is cut short when the pane is narrow.
        ui.horizontal(|ui| {
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.small_button("Hide").on_hover_text("Hide the bases. The Sequence button above brings them back.").clicked() {
                    edits.push(Edit::ShowSequence(false));
                }
                ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
                    ui.add_space(12.0);
                    ui.label(
                        RichText::new(format!("{} · {} bp", target.doc.name, thousands(target.doc.len()))).font(style::bold_font(12.5)).color(palette.ink),
                    );
                    let what = match target.setting {
                        Some(_) => "Drag across the bases to mark the part. The map follows, and the map's marks show here.",
                        None => "The bases of what the picture shows, with the same stretch highlighted. A selection here shows on the picture too.",
                    };
                    ui.add(egui::Label::new(RichText::new(what).small().color(palette.muted)).truncate());
                });
            });
        });
        seq_view::show(ui, &mut preview.tab, settings);
    });
    // A selection made here: while the button is held it shows on the picture, and once it is let
    // go it becomes the part when the picture sets one. A feature chosen by its bar goes by its
    // name, as a click on the map's feature does.
    let selected = preview.tab.selection;
    if selected == preview.shown {
        return None;
    }
    let pointer_down = ui.input(|i| i.pointer.primary_down());
    match (selected, target.setting) {
        (Some(span), Some(setting)) if !pointer_down && interactive => {
            let len = target.doc.len();
            let named = preview
                .tab
                .selected_feature
                .and_then(|uid| target.doc.features.iter().find(|f| f.uid == uid))
                .filter(|f| !f.name.is_empty() && target.doc.features_named(&f.name).count() == 1 && f.span(len).normalize(len) == span.normalize(len))
                .map(|f| f.name.clone());
            edits.push(Edit::SetPart(setting, Some(named.unwrap_or_else(|| span.format(len)))));
            preview.pending = Some((preview.shown, Instant::now()));
            preview.shown = Some(span);
            Some(span)
        }
        (Some(span), _) => Some(span),
        // A click cleared the selection. The highlight comes back, as there is nothing to mark.
        (None, _) => {
            preview.tab.select(preview.shown);
            None
        }
    }
}

/// The handle between a zoomed card's picture and its sequence pane: dragging it moves the split.
fn split_handle(ui: &mut egui::Ui, tab: &ProtocolTab, area: Rect, pane: Rect, palette: &Palette, edits: &mut Vec<Edit>, interactive: bool) {
    let handle = Rect::from_min_max(pos2(area.min.x, pane.min.y - SPLIT_GAP), pos2(area.max.x, pane.min.y));
    ui.painter().rect_filled(Rect::from_center_size(handle.center(), vec2(40.0, 3.0)), CornerRadius::same(2), palette.line);
    if !interactive {
        return;
    }
    let response = ui.interact(handle, egui::Id::unique(("protocol-split", tab.id)), Sense::drag());
    if response.hovered() || response.dragged() {
        ui.ctx().set_cursor_icon(CursorIcon::ResizeVertical);
    }
    if response.dragged() {
        let room = (area.height() - DETAIL_HEADER).max(1.0);
        edits.push(Edit::SequenceShare(tab.sequence_share - response.drag_delta().y / room));
    }
}

/// The detailed view of a card: a way back, the step's headline, a hint, the large picture, and
/// the molecule's bases under it.
#[allow(clippy::too_many_arguments)]
fn detail(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    board: &Board,
    card: usize,
    near: Option<&(usize, NearestCuts)>,
    settings: &ViewSettings,
    palette: &Palette,
    edits: &mut Vec<Edit>,
    interactive: bool,
    preview: &mut Option<Preview>,
) {
    let area = ui.max_rect();
    ui.painter().rect_filled(area, 0.0, palette.canvas);
    let c = &board.cards[card];
    let Some(step) = c.step else { return };
    let target = preview_target(tab, board, card);
    // A plasmid's map is where its cut or its PCR is chosen, so the map's heading is that step's.
    let headline_step = match c.kind {
        Kind::Plasmid => board.part_step_after(card).unwrap_or(step),
        _ => step,
    };
    let told = protocol::story(&tab.protocol, &tab.run, headline_step);
    let header = Rect::from_min_size(area.min + vec2(14.0, 10.0), vec2(area.width() - 28.0, DETAIL_HEADER - 14.0));
    let mut top = ui.new_child(UiBuilder::new().max_rect(header).id_salt(("protocol-detail-header", tab.id)));
    top.horizontal(|ui| {
        if ui.add_enabled(interactive, egui::Button::new("← All steps")).on_hover_text("Zoom back out to the whole protocol (Escape)").clicked() {
            edits.push(Edit::Zoom(None));
        }
        ui.add_space(8.0);
        // The bases under the picture can be hidden and shown again from here. The switch is
        // placed first, so a long headline is cut short rather than run under it.
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if target.is_some() {
                let shown = tab.show_sequence;
                let hint = "Show the molecule's bases under the picture, with the same stretch highlighted. A stretch selected in either shows in both.";
                ui.add_enabled_ui(interactive, |ui| {
                    if ui.selectable_label(shown, "Sequence").on_hover_text(hint).clicked() {
                        edits.push(Edit::ShowSequence(!shown));
                    }
                });
                ui.add_space(8.0);
            }
            ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
                ui.add(egui::Label::new(RichText::new(&told.headline).font(style::bold_font(18.0)).color(palette.ink)).truncate());
            });
        });
    });
    let waits = tab.run.steps.get(step).is_some_and(|r| r.status == Status::NeedsChoice);
    // A cut without pieces yet shows the plasmid's map, so its hint is the map's. The map's hint
    // follows the work: first mark the part, then click the cuts that the map then shows.
    let unpieced = tab.run.steps.get(step).is_none_or(|r| r.outputs.is_empty());
    let map_cut = match c.kind {
        Kind::Plasmid => board.part_step_after(card),
        Kind::Copy => Some(step),
        Kind::Cut if unpieced => board.before(card).filter(|&b| board.cards[b].kind == Kind::Copy).and_then(|b| board.cards[b].step).or(Some(step)),
        _ => None,
    };
    let marked = map_cut.and_then(|s| tab.run.steps.get(s)).is_some_and(|r| r.findings.region.is_some());
    let copying = map_cut.is_some_and(|s| matches!(tab.protocol.steps[s].action, Action::Amplify(_)));
    let marking = match (copying, marked) {
        (true, true) => {
            "Drag across the map to mark the part to copy, or click a feature. The primers anneal at its ends and are drawn on the map. Choose the enzymes for their tails in the panel."
        }
        (true, false) => "Drag across the map to mark the part to copy, or click a feature. Primers are then designed for it.",
        (false, true) => {
            "Drag across the map to mark the part, or click a feature. Click an enzyme to cut there: the chosen cuts are bold, and a second click takes one off."
        }
        (false, false) => "Drag across the map to mark the part you want, or click a feature. The enzymes that cut nearest it then appear on the map.",
    };
    let others = tab.run.steps.get(step).is_some_and(|r| !r.findings.byproducts.is_empty());
    let hint = match c.kind {
        Kind::Plasmid if map_cut.is_none() => "This plasmid is not cut yet. On the bench, drag a green piece onto it to put the piece in.",
        Kind::Plasmid | Kind::Copy => marking,
        Kind::Cut if unpieced && board.before(card).is_some() => marking,
        Kind::Cut => "Green pieces go on into the new plasmid. Red pieces are left over. Click a piece to keep it instead.",
        Kind::Join if waits && others => {
            "The insert fits either way round. The arrowheads show which way each gene reads. Choose the plasmid you want. The grey plasmids on the right form in the same tube."
        }
        Kind::Join if waits => "The insert fits either way round. The arrowheads show which way each gene reads. Choose the plasmid you want.",
        Kind::Join if others => "The new plasmid, with the other plasmids its tube makes on the right. The shaded stretch is the insert.",
        Kind::Join => "The new plasmid. The shaded stretch is the insert, and the marks show where the pieces were joined.",
        _ => "",
    };
    top.label(RichText::new(hint).color(palette.muted));

    // The picture, with the molecule's bases under it when the view shows them. The sequence pane
    // is drawn first, so a stretch selected in it shows on the picture this frame.
    let (rect, pane) = detail_rects(area, sequence_share(tab, board, card));
    let mut live: Option<Span> = None;
    if let (Some(target), Some(pane)) = (target.as_ref(), pane) {
        let map_live = target.setting.and_then(|_| map_drag_span(ui, tab, board, card));
        live = sequence_pane(ui, tab, preview, target, pane, map_live, settings, palette, edits, interactive);
        split_handle(ui, tab, area, pane, palette, edits, interactive);
    }
    let ctx = ui.ctx().clone();
    let measure = EguiMeasure { ctx: &ctx };
    match c.kind {
        Kind::Plasmid | Kind::Copy => marking_map(ui, tab, board, card, near, settings, palette, rect, edits, interactive, &measure, live),
        Kind::Cut => {
            // Before the cut has pieces, the map of the plasmid it cuts (or that its PCR copies
            // from) stands in, so the part and the cuts are chosen on the map itself.
            let Some(shown) = picture(tab, board, card, true, settings, palette.canvas) else {
                match board.before(card).filter(|&b| matches!(board.cards[b].kind, Kind::Plasmid | Kind::Copy)) {
                    Some(plasmid) => marking_map(ui, tab, board, plasmid, near, settings, palette, rect, edits, interactive, &measure, live),
                    None => {
                        centered_note(ui, rect, "Choose an enzyme on each side of the part in the panel on the right, and the pieces appear here.", palette)
                    }
                }
                return;
            };
            let mut options = shown.options;
            options.width = rect.width();
            options.height = rect.height();
            options.font_size = 12.5;
            let scene = ring::render_ring(&shown.doc, &options, &measure);
            let view = cut_view(tab, c);
            let pieces = view.as_ref().map(|v| v.pieces.clone()).unwrap_or_default();
            let describe = |tag: &Tag| match tag {
                Tag::Piece(ii) => pieces.get(*ii).and_then(|(span, role, output)| {
                    let piece = output.and_then(|o| tab.run.steps.get(step)?.outputs.get(o))?;
                    let what = match role {
                        Role::Kept => "kept: it goes on into the new plasmid",
                        Role::LeftOver => "left over: click to keep it instead",
                        Role::Plain => "not chosen yet: click to keep it",
                    };
                    Some(format!(
                        "{} bp, {} … {}\n{}\nbases {}",
                        thousands(piece.doc.len()),
                        piece.describe_end(true),
                        piece.describe_end(false),
                        what,
                        span.format(shown.doc.len())
                    ))
                }),
                Tag::Feature(uid) => shown.doc.features.iter().find(|f| f.uid == *uid).map(|f| format!("{}\n{} · {} bp", f.name, f.kind, thousands(f.len()))),
                Tag::Enzyme(name, pos) => Some(format!("{name} cuts after base {}", thousands(*pos))),
                _ => None,
            };
            if let Some(Tag::Piece(ii)) =
                draw_scene(ui, egui::Id::unique(("protocol-detail-scene", tab.id, step)), &scene, rect, interactive, &measure, describe)
                && let (Some(&pick), Some(Some(output))) = (c.picks.first(), pieces.get(ii).map(|p| p.2))
            {
                edits.push(Edit::Keep(pick, output));
            }
        }
        Kind::Join => {
            let Some(result) = tab.run.steps.get(step) else { return };
            // The other plasmids of the tube take a strip on the right, beside the plan's.
            let main = join_main_rect(result, rect);
            if main != rect {
                let strip = Rect::from_min_max(pos2(main.max.x, rect.min.y), rect.max);
                other_plasmids(ui, tab, step, result, strip, settings, palette, edits, interactive);
            }
            let rect = main;
            match result.status {
                Status::Done => {
                    let Some(product) = result.outputs.first() else { return };
                    let doc = &product.doc;
                    let marks = junction_marks(&result.findings.junctions);
                    // The insert (the last part) is outlined, so the eye finds where it went, or
                    // the stretch selected in the sequence under the map.
                    let options = MapOptions {
                        width: rect.width(),
                        height: rect.height(),
                        theme: theme_on(settings, palette.canvas),
                        font_size: 12.0,
                        enzymes: marks,
                        selection: live.or(result.findings.parts.last().copied().filter(|_| result.findings.parts.len() > 1)),
                        top: product_top(&result.findings.parts, doc.len()),
                        title: Some("New plasmid".into()),
                        colors: settings.colors,
                        show_primers: false,
                        max_feature_labels: 40,
                        max_labels_per_side: 26,
                        ..MapOptions::default()
                    };
                    let scene = render_map(doc, &options, &measure);
                    let clicked = draw_scene(ui, egui::Id::unique(("protocol-detail-scene", tab.id, step)), &scene, rect, interactive, &measure, |tag| {
                        describe_tag(doc, tag, None, &[], false)
                    });
                    // A feature clicked on the map is selected in the sequence under it.
                    if let (Some(Tag::Feature(uid)), Some(p)) = (clicked, preview.as_mut()) {
                        p.tab.select_feature(uid);
                    }
                }
                Status::NeedsChoice => choice_rings(ui, tab, step, result, rect, settings, palette, edits, interactive),
                _ => {
                    let text = result.error.as_ref().map(|e| e.message.clone()).unwrap_or_else(|| "The join waits for the pieces above it.".into());
                    centered_note(ui, rect, &text, palette);
                }
            }
        }
        _ => centered_note(ui, rect, &told.details.join(" "), palette),
    }
}

/// The plasmid step a card's detailed map shows, and the step whose part the map sets: the cut or
/// the PCR after a plasmid, or the PCR itself.
fn map_of(board: &Board, card: usize) -> Option<(usize, Option<usize>)> {
    let c = &board.cards[card];
    match c.kind {
        Kind::Plasmid => Some((c.step?, board.part_step_after(card))),
        Kind::Copy => Some((board.plasmid_before(card).and_then(|p| board.cards[p].step)?, c.step)),
        _ => None,
    }
}

/// The detailed map of plasmid card `card` (or of the plasmid a PCR card copies from), on which
/// the part of the cut or the PCR after it is chosen: a drag across the map marks the part, a
/// click on a feature makes it the part, and for a cut a click on an enzyme cuts there or takes
/// that cut off.
#[allow(clippy::too_many_arguments)]
fn marking_map(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    board: &Board,
    card: usize,
    near: Option<&(usize, NearestCuts)>,
    settings: &ViewSettings,
    palette: &Palette,
    rect: Rect,
    edits: &mut Vec<Edit>,
    interactive: bool,
    measure: &EguiMeasure,
    live: Option<Span>,
) {
    let Some((step, _)) = map_of(board, card) else { return };
    // A drag across the map marks the part. Where it started and where the pointer is
    // now are kept from frame to frame while the button is held, and the stretch between
    // them shows as the map's selection.
    let drag_id = egui::Id::unique(("protocol-map-drag", tab.id, step));
    let marking: Option<(usize, usize)> = ui.data(|d| d.get_temp(drag_id));
    let Some(doc) = tab.run.steps.get(step).and_then(|r| r.outputs.first()).map(|p| &p.doc) else { return };
    let (len, circular) = (doc.len(), doc.is_circular());
    let span = marking.and_then(|(a, b)| map_view::drag_span(a, b, len, circular));
    // The stretch selected in the sequence under the map shows while nothing is dragged here.
    let Some(map) = plasmid_map(tab, board, card, near, settings, palette, rect.size(), span.or(live), measure) else { return };
    paint(&ui.painter_at(rect), &map.scene, rect.min);
    if !interactive {
        return;
    }
    let response = ui.interact(rect, egui::Id::unique(("protocol-detail-map", tab.id, step)), Sense::click_and_drag());
    let position_at = |p: Pos2| if circular { map_view::circular_position(rect, p, len, map.top) } else { map_view::linear_position(rect, p, len) };
    let tag = response.hover_pos().and_then(|p| hit(&map.scene, [p.x - rect.min.x, p.y - rect.min.y], measure));
    if map.cut.is_some() {
        // The drag starts where the button went down, not where it had moved to by the
        // time it counted as a drag.
        if response.drag_started()
            && let Some(origin) = ui.input(|i| i.pointer.press_origin()).or(response.interact_pointer_pos())
        {
            let at = position_at(origin);
            ui.data_mut(|d| d.insert_temp(drag_id, (at, at)));
        }
        if response.dragged()
            && let (Some((anchor, _)), Some(p)) = (marking, response.interact_pointer_pos())
        {
            ui.data_mut(|d| d.insert_temp(drag_id, (anchor, position_at(p))));
        }
        if response.drag_stopped() {
            ui.data_mut(|d| d.remove::<(usize, usize)>(drag_id));
            if let (Some(cut), Some(span)) = (map.cut, span) {
                edits.push(Edit::SetPart(cut, Some(span.format(len))));
            }
        }
    }
    // While marking, a label beside the pointer gives the stretch and its length.
    if response.dragged()
        && let (Some(span), Some(p)) = (span, response.interact_pointer_pos())
    {
        let text = format!("{} · {} bp", span.format(len), thousands(span.len()));
        let galley = ui.painter().layout_no_wrap(text, style::bold_font(12.5), palette.raised);
        let label = Rect::from_min_size(p + vec2(14.0, 12.0), galley.size() + vec2(14.0, 8.0));
        ui.painter().rect_filled(label, CornerRadius::same(6), palette.accent);
        ui.painter().galley(label.center() - galley.size() / 2.0, galley, palette.raised);
        ui.ctx().set_cursor_icon(CursorIcon::Crosshair);
    }
    if response.clicked()
        && let Some(cut) = map.cut
    {
        match &tag {
            // A feature becomes the part the cut is about. A name that several features
            // share would find the first of them, so such a feature goes by its bases.
            Some(Tag::Feature(uid)) => {
                if let Some(feature) = doc.features.iter().find(|f| f.uid == *uid) {
                    let shared = doc.features_named(&feature.name).count() > 1 || feature.name.is_empty();
                    let part = if shared { feature.span(len).format(len) } else { feature.name.clone() };
                    edits.push(Edit::SetPart(cut, Some(part)));
                }
            }
            // An enzyme cuts on its side of the part. A chosen cut that is not among the
            // nearest comes off the side it makes.
            Some(Tag::Enzyme(name, pos)) => {
                if let Some(before) = enzyme_side(tab, cut, map.listed, name, *pos) {
                    edits.push(Edit::PickEnzyme(cut, before, name.clone()));
                }
            }
            _ => {}
        }
    }
    if !response.dragged() && response.hovered() {
        let clickable = matches!(tag, Some(Tag::Feature(_) | Tag::Enzyme(..)));
        ui.ctx().set_cursor_icon(if clickable { CursorIcon::PointingHand } else { CursorIcon::Crosshair });
        let chosen = match map.cut.map(|c| &tab.protocol.steps[c].action) {
            Some(Action::Digest(d)) => d.enzymes.clone(),
            _ => Vec::new(),
        };
        let copying = matches!(map.cut.map(|c| &tab.protocol.steps[c].action), Some(Action::Amplify(_)));
        let text = match &tag {
            Some(Tag::Primer(uid)) => map.primers.iter().find(|p| p.uid == *uid).map(|p| format!("{}\n{}\n{}", p.name, p.seq, p.description)),
            Some(tag) => describe_tag(doc, tag, map.listed, &chosen, copying),
            None => Some("Drag around the map to mark the part".to_string()),
        };
        if let Some(text) = text {
            let _ = response.on_hover_text_at_pointer(text);
        }
    }
}

/// A note in the middle of an empty picture.
fn centered_note(ui: &mut egui::Ui, rect: Rect, text: &str, palette: &Palette) {
    let galley = ui.painter().layout(text.to_string(), style::ui_font(14.0), palette.muted, (rect.width() * 0.6).max(200.0));
    ui.painter().galley(rect.center() - galley.size() / 2.0, galley, palette.muted);
}

/// What a tag of a map stands for, for a tooltip: a feature's name, type and length, or an enzyme
/// cut with how far it lies from the part. `copying` says the map's part is copied by PCR.
fn describe_tag(doc: &Document, tag: &Tag, listed: Option<&NearestCuts>, chosen: &[String], copying: bool) -> Option<String> {
    match tag {
        Tag::Feature(uid) => {
            let action = if copying { "Click to copy it." } else { "Click to cut around it." };
            doc.features.iter().find(|f| f.uid == *uid).map(|f| format!("{}\n{} · {} bp\n{action}", f.name, f.kind, thousands(f.len())))
        }
        Tag::Enzyme(name, pos) => {
            let side = listed.and_then(|l| nearer_side(l, name, *pos)).map(|(before, c)| {
                let word = if before { "before" } else { "after" };
                format!("{} bp {word} the part, leaving a {} end", thousands(c.distance), c.end_text())
            });
            let is_chosen = chosen.iter().any(|c| c.eq_ignore_ascii_case(name));
            let action = match (is_chosen, side.is_some()) {
                (true, _) => "It is one of the chosen cuts. Click to take it off.",
                (false, true) => "Click to cut here.",
                (false, false) => "",
            };
            let lines: Vec<String> = [Some(format!("{name} cuts after base {}", thousands(*pos))), side, Some(action.to_string())]
                .into_iter()
                .flatten()
                .filter(|l| !l.is_empty())
                .collect();
            Some(lines.join("\n"))
        }
        _ => None,
    }
}

/// The side of the part on which the cut of `name` at `pos` is listed, true for before it, and
/// the listing. On a small circle one cut can be among the nearest on both sides, as it is near
/// one end going one way and the other end going the other way. It then counts on the side it
/// is nearer.
fn nearer_side<'a>(listed: &'a NearestCuts, name: &str, pos: usize) -> Option<(bool, &'a NearCut)> {
    let find = |side: &'a [NearCut]| side.iter().find(|c| c.enzyme == name && c.pos == pos);
    match (find(&listed.before), find(&listed.after)) {
        (Some(before), Some(after)) => Some(if before.distance <= after.distance { (true, before) } else { (false, after) }),
        (Some(before), None) => Some((true, before)),
        (None, Some(after)) => Some((false, after)),
        (None, None) => None,
    }
}

/// Which side of cut `cut`'s part a click on enzyme `name` at `pos` chooses: true for before it.
/// An enzyme among the nearest chooses the side it is nearer. A chosen enzyme elsewhere comes off
/// the side it makes, if it makes one.
fn enzyme_side(tab: &ProtocolTab, cut: usize, listed: Option<&NearestCuts>, name: &str, pos: usize) -> Option<bool> {
    if let Some((before, _)) = listed.and_then(|l| nearer_side(l, name, pos)) {
        return Some(before);
    }
    let Some(Action::Digest(params)) = tab.protocol.steps.get(cut).map(|s| &s.action) else { return None };
    let (left, right) = flanks(tab.run.steps.get(cut), &params.enzymes);
    if left.as_deref() == Some(name) {
        Some(true)
    } else if right.as_deref() == Some(name) {
        Some(false)
    } else {
        None
    }
}

/// A plasmid's detailed map, and what a click on it acts on.
struct PlasmidMap<'a> {
    scene: Scene,
    /// The cut or the PCR that follows the plasmid, whose part (and for a cut, enzymes) the map
    /// sets.
    cut: Option<usize>,
    /// The enzymes listed as nearest the part, the same as the panel's.
    listed: Option<&'a NearestCuts>,
    /// The base at 12 o'clock.
    top: usize,
    /// A PCR's primers, drawn on the map.
    primers: Vec<Primer>,
}

/// The detailed map of plasmid card `card` at `size`: the part (or the stretch `marking` being
/// dragged across) shaded, the enzymes nearest the part on each side marked, and the chosen cuts
/// drawn in the accent colour.
#[allow(clippy::too_many_arguments)]
fn plasmid_map<'a>(
    tab: &ProtocolTab,
    board: &Board,
    card: usize,
    near: Option<&'a (usize, NearestCuts)>,
    settings: &ViewSettings,
    palette: &Palette,
    size: egui::Vec2,
    marking: Option<Span>,
    measure: &EguiMeasure,
) -> Option<PlasmidMap<'a>> {
    let (step, setting) = map_of(board, card)?;
    let template = &tab.run.steps.get(step)?.outputs.first()?.doc;
    // The cut's pieces and cuts, or the PCR's part and primers, whichever follows the plasmid.
    let setting_card = setting.and_then(|s| board.card_of(s));
    let view = setting_card.filter(|&cc| board.cards[cc].kind == Kind::Cut).and_then(|cc| cut_view(tab, &board.cards[cc]));
    let copy = setting.filter(|&s| matches!(tab.protocol.steps[s].action, Action::Amplify(_)));
    let copied = copy.and_then(|s| copy_region(tab, s));
    let primers: Vec<Primer> = copy.and_then(|s| tab.run.steps.get(s)).map(|r| r.findings.primers.clone()).unwrap_or_default();
    let listed = near.filter(|(listed_for, _)| Some(*listed_for) == setting).map(|(_, n)| n);
    let mut candidates: Vec<SiteMark> = Vec::new();
    for near in listed.iter().flat_map(|l| l.before.iter().chain(&l.after)) {
        if !candidates.iter().any(|m| m.pos == near.pos && m.name == near.enzyme) {
            candidates.push(SiteMark { name: near.enzyme.clone(), pos: near.pos });
        }
    }
    let top = view.as_ref().map(|v| v.top).or_else(|| copied.map(|(r, len, _)| (r.start + r.len() / 2) % len)).unwrap_or(0);
    // A PCR's primers are drawn in place of the map's own.
    let mut doc = template.clone();
    if !primers.is_empty() {
        doc.primers = primers.clone();
    }
    let options = MapOptions {
        width: size.x,
        height: size.y,
        theme: theme_on(settings, palette.canvas),
        font_size: 12.0,
        enzymes: candidates,
        chosen_enzymes: view.as_ref().map(|v| v.cuts.clone()).unwrap_or_default(),
        selection: marking.or(view.as_ref().and_then(|v| v.region)).or(copied.map(|(r, _, _)| r)),
        top,
        colors: settings.colors,
        show_primers: !primers.is_empty(),
        max_feature_labels: 40,
        max_labels_per_side: 26,
        ..MapOptions::default()
    };
    Some(PlasmidMap { scene: render_map(&doc, &options, measure), cut: setting, listed, top, primers })
}

/// Paint a scene in `rect`, name what the pointer is over, and return the tag clicked.
fn draw_scene(
    ui: &mut egui::Ui,
    id: egui::Id,
    scene: &Scene,
    rect: Rect,
    interactive: bool,
    measure: &EguiMeasure,
    describe: impl Fn(&Tag) -> Option<String>,
) -> Option<Tag> {
    paint(&ui.painter_at(rect), scene, rect.min);
    if !interactive {
        return None;
    }
    let response = ui.interact(rect, id, Sense::click());
    let tag = response.hover_pos().and_then(|p| hit(scene, [p.x - rect.min.x, p.y - rect.min.y], measure));
    let clickable = matches!(tag, Some(Tag::Feature(_) | Tag::Enzyme(..) | Tag::Piece(_)));
    if clickable {
        ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
    }
    let clicked = response.clicked();
    if let Some(text) = tag.as_ref().and_then(&describe) {
        let _ = response.on_hover_text_at_pointer(text);
    }
    if clicked { tag } else { None }
}

/// The options of a join that waits for a choice, side by side, each as the ring of the plasmid it
/// gives and a button that chooses it.
#[allow(clippy::too_many_arguments)]
fn choice_rings(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    step: usize,
    result: &StepRun,
    rect: Rect,
    settings: &ViewSettings,
    palette: &Palette,
    edits: &mut Vec<Edit>,
    interactive: bool,
) {
    let count = result.choices.len().max(1);
    let width = rect.width() / count as f32;
    let ctx = ui.ctx().clone();
    let measure = EguiMeasure { ctx: &ctx };
    for (ii, choice) in result.choices.iter().enumerate() {
        let column = Rect::from_min_size(pos2(rect.min.x + ii as f32 * width, rect.min.y), vec2(width, rect.height()));
        let Some(product) = choice.outputs.first() else { continue };
        // Each option is drawn as the bench draws it, only larger.
        let shown = Shown { made: Made::Option(ii), piece: product, parts: &choice.parts };
        let mut options = made_ring(&shown, settings, palette.canvas);
        options.width = column.width();
        options.height = column.height() - CHOICE_BELOW;
        options.font_size = 12.0;
        let scene = ring::render_ring(&product.doc, &options, &measure);
        paint(&ui.painter_at(column), &scene, column.min);
        let below = Rect::from_min_max(pos2(column.min.x + 16.0, column.max.y - CHOICE_BELOW + 4.0), pos2(column.max.x - 16.0, column.max.y));
        let mut child = ui.new_child(UiBuilder::new().max_rect(below).id_salt(("protocol-choice", tab.id, ii)));
        child.vertical_centered(|ui| {
            ui.label(RichText::new(&choice.label).small().color(palette.muted));
            if ui.add_enabled(interactive, egui::Button::new("Use this")).clicked() {
                edits.push(Edit::Choose(step, ii));
            }
        });
    }
}

/// The other plasmids the tube of join step `step` makes, one under another in a strip at the
/// right of its detailed view. Each shows its ring and what it is in plain words. One that holds
/// every piece in order can become the product, and a backbone closed on itself can be stopped by
/// taking its phosphates off.
#[allow(clippy::too_many_arguments)]
fn other_plasmids(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    step: usize,
    result: &StepRun,
    rect: Rect,
    settings: &ViewSettings,
    palette: &Palette,
    edits: &mut Vec<Edit>,
    interactive: bool,
) {
    // A faint line parts the strip from the plan's plasmids.
    ui.painter().line_segment([rect.left_top() + vec2(0.0, 8.0), rect.left_bottom() - vec2(0.0, 8.0)], Stroke::new(1.0, palette.line));
    let mut strip = ui.new_child(UiBuilder::new().max_rect(rect.shrink2(vec2(14.0, 4.0))).id_salt(("protocol-others", tab.id, step)));
    strip.add_enabled_ui(interactive, |ui| {
        egui::ScrollArea::vertical().id_salt(("protocol-others-scroll", tab.id, step)).auto_shrink([false, false]).show(ui, |ui| {
            ui.label(RichText::new("The same tube also makes").font(style::bold_font(14.5)).color(palette.ink));
            ui.label(RichText::new("A colony can carry any of these. Check colonies by a digest or by sequencing.").color(palette.muted));
            let ctx = ui.ctx().clone();
            let measure = EguiMeasure { ctx: &ctx };
            let backbone = tab.protocol.steps.get(step).and_then(|s| s.inputs.first()).map(|r| protocol::piece_name(&tab.protocol, &tab.run, r));
            for (ii, other) in result.findings.byproducts.iter().enumerate() {
                ui.add_space(10.0);
                let shown = Shown { made: Made::Other(ii), piece: &other.piece, parts: &other.parts };
                let (ring_rect, _) = ui.allocate_exact_size(vec2(ui.available_width(), 160.0), Sense::hover());
                let mut options = made_ring(&shown, settings, palette.canvas);
                options.width = ring_rect.width();
                options.height = ring_rect.height();
                options.font_size = 11.0;
                let scene = ring::render_ring(&other.piece.doc, &options, &measure);
                paint(&ui.painter_at(ring_rect), &scene, ring_rect.min);
                let words = protocol::byproduct_words(&tab.protocol, &tab.run, step, other);
                ui.label(RichText::new(uppercase_first(&words)).color(palette.ink));
                ui.label(RichText::new(other.piece.summary()).small().color(palette.muted));
                if other.pin.is_some()
                    && ui.button("Use this instead").on_hover_text("Make this the new plasmid, by turning the pieces that go in the other way round").clicked()
                {
                    edits.push(Edit::UseInstead(step, ii));
                }
                if other.alone()
                    && let Some(backbone) = &backbone
                    && ui.button(format!("Remove the 5' phosphates of {backbone}")).on_hover_text("Without them it cannot close on itself").clicked()
                {
                    edits.push(Edit::Dephosphorylate(step, 0));
                }
            }
        });
    });
}

/// Junctions as marks for a map, labelled J1, J2 and so on with their enzymes.
fn junction_marks(junctions: &[protocol::Junction]) -> Vec<SiteMark> {
    junctions
        .iter()
        .enumerate()
        .map(|(ii, j)| {
            let mut names: Vec<&str> = [&j.left_enzyme, &j.right_enzyme].iter().filter_map(|e| e.as_deref()).collect();
            names.dedup();
            let label = if names.is_empty() { format!("J{}", ii + 1) } else { format!("J{} {}", ii + 1, names.join("/")) };
            SiteMark { name: label, pos: j.position }
        })
        .collect()
}

/// A gel of a digest: a ladder, the uncut molecule and its fragments.
fn gel_scene(input: &Piece, fragments: &[Piece], label: &str) -> Scene {
    let mut lanes = Vec::new();
    if let Some(ladder) = Ladder::get("1kb") {
        lanes.push(Lane::ladder(ladder));
    }
    let length = input.doc.len();
    lanes.push(if input.doc.is_circular() { Lane::uncut_plasmid("uncut", length) } else { Lane::digest("uncut", &[length]) });
    let sizes: Vec<usize> = fragments.iter().map(|f| f.doc.len()).collect();
    lanes.push(Lane::digest(label, &sizes));
    render_gel(&lanes, &GelOptions { height: 200.0, lane_width: 40.0, ..GelOptions::default() })
}

/// The panel: the protocol in a sentence, the step shown in plain words with its choices, and an
/// editor for the step by hand.
#[allow(clippy::too_many_arguments)]
fn step_panel(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    board: &Board,
    docs: &[DocTab],
    near: Option<&(usize, NearestCuts)>,
    tails: Option<&(usize, Vec<TailChoice>)>,
    palette: &Palette,
    edits: &mut Vec<Edit>,
) {
    egui::ScrollArea::vertical().id_salt(("protocol-panel-scroll", tab.id)).auto_shrink([false, false]).show(ui, |ui| {
        ui.set_width(ui.available_width());
        protocol_section(ui, tab, palette, edits);
        ui.separator();
        if tab.protocol.steps.is_empty() {
            start_section(ui, tab, docs, palette, edits, None);
            return;
        }
        let index = tab.selected.min(tab.protocol.steps.len() - 1);
        let step = &tab.protocol.steps[index];
        let result = tab.run.steps.get(index);
        let told = protocol::story(&tab.protocol, &tab.run, index);
        // What the step is, in plain words, with the lab's word for it.
        let kind = match step.action.lab_term() {
            Some(term) => format!("{} · {term}", step.action.title()),
            None => step.action.title().to_string(),
        };
        ui.horizontal(|ui| {
            ui.label(RichText::new(kind.to_uppercase()).font(style::bold_font(10.5)).color(palette.accent));
            if let Some(result) = result {
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(RichText::new(result.status.name()).small().color(status_color(result.status, palette)));
                });
            }
        });
        ui.label(RichText::new(&told.headline).font(style::bold_font(15.5)).color(palette.ink));
        for line in &told.details {
            ui.label(RichText::new(line).color(palette.ink));
        }
        if !step.note.is_empty() {
            ui.label(RichText::new(&step.note).italics().color(palette.muted));
        }
        if let Some(result) = result {
            problems(ui, result, palette);
        }
        ui.add_space(6.0);
        // The nearest cuts listed for the cut in focus.
        let listed_for = |cut: usize| near.filter(|(listed, _)| *listed == cut).map(|(_, n)| n);
        match &step.action {
            Action::Open(params) => {
                source_section(ui, tab, index, params, docs, palette, edits);
                // A map's cut or PCR is chosen while looking at the map, so its choices follow here.
                if let Some(copy) = focused_copy(tab, board) {
                    ui.separator();
                    copy_section(ui, tab, copy, tails, docs, palette, edits);
                } else if let Some(cut) = focused_cut(tab, board) {
                    ui.separator();
                    cut_section(ui, tab, board, cut, listed_for(cut), docs, palette, edits);
                }
            }
            Action::Digest(_) if tab.protocol.amplify_before(&step.id).is_some() => trim_section(ui, tab, board, index, palette, edits),
            Action::Digest(_) => cut_section(ui, tab, board, index, listed_for(index), docs, palette, edits),
            Action::Amplify(_) => copy_section(ui, tab, index, tails, docs, palette, edits),
            Action::Pick(_) => pieces_section(ui, tab, board, index, palette, edits),
            Action::Dephosphorylate(_) | Action::Flip(_) => {
                if ui.button("Leave this step out").on_hover_text("Take the step out and pass its input straight on").clicked() {
                    edits.push(Edit::LeaveOut(index));
                }
            }
            Action::Ligate(_) => join_section(ui, tab, index, palette, edits),
        }
        ui.add_space(10.0);
        egui::CollapsingHeader::new("Change this step by hand").id_salt(("protocol-by-hand", tab.id)).show(ui, |ui| step_editor(ui, tab, palette, edits));
        ui.add_space(4.0);
        ui.menu_button("+ Add a step", |ui| {
            for &(kind, what) in Action::KINDS {
                let title = Action::blank(kind).map(|a| a.title()).unwrap_or(kind);
                if ui.button(format!("{title}: {what}")).clicked() {
                    edits.push(Edit::AddStep(kind));
                    ui.close();
                }
            }
        });
        ui.add_space(12.0);
    });
}

/// The protocol's name, the protocol in a sentence, where it is saved and how the run went.
fn protocol_section(ui: &mut egui::Ui, tab: &ProtocolTab, palette: &Palette, edits: &mut Vec<Edit>) {
    section_title(ui, "Protocol", palette);
    if let Some(name) = committed_text(ui, egui::Id::unique(("protocol-name", tab.id)), &tab.protocol.name, "name") {
        edits.push(Edit::Change(None, Box::new(move |p| p.name = name)));
    }
    let said = protocol::summary(&tab.protocol, &tab.run);
    if !said.is_empty() {
        ui.label(RichText::new(said).color(palette.muted));
    }
    let place = tab.path.as_ref().map(|p| p.display().to_string()).unwrap_or_else(|| "not saved yet".into());
    ui.label(RichText::new(place).small().color(palette.faint));
    ui.horizontal(|ui| {
        ui.label(RichText::new(tab.progress()).small().color(palette.muted));
        if ui.small_button("Read the maps again").on_hover_text("Run every step again with the maps as they are on disk now").clicked() {
            edits.push(Edit::Rerun);
        }
    });
}

/// Warnings, and the problem that stopped the step with what would fix it.
fn problems(ui: &mut egui::Ui, result: &StepRun, palette: &Palette) {
    for warning in &result.warnings {
        ui.label(RichText::new(format!("Note: {}", warning.message)).color(palette.warn));
    }
    if let Some(error) = &result.error {
        ui.add_space(4.0);
        egui::Frame::new().stroke(Stroke::new(1.0, palette.error)).corner_radius(CornerRadius::same(6)).inner_margin(egui::Margin::same(8)).show(ui, |ui| {
            ui.set_width(ui.available_width());
            ui.label(RichText::new(&error.message).color(palette.error));
            for fix in &error.fixes {
                ui.label(RichText::new(format!("Try: {fix}")).color(palette.ink));
            }
        });
    }
}

/// The maps open in the app, and a file chooser, to start from or to put an insert into.
fn map_choices(ui: &mut egui::Ui, tab: &ProtocolTab, docs: &[DocTab], palette: &Palette) -> Option<OpenParams> {
    let folder = tab.folder();
    let mut chosen = None;
    if ui.button("Choose a file…").clicked()
        && let Some(file) = rfd::FileDialog::new().add_filter("Sequence files", &["dna", "gb", "gbk", "genbank", "fasta", "fa", "json"]).pick_file()
    {
        chosen = Some(OpenParams { path: Some(protocol::stored_path(&folder, &file)), document: None });
    }
    let open: Vec<&DocTab> = docs.iter().filter(|d| d.doc.molecule == Molecule::Dna).collect();
    if !open.is_empty() {
        ui.label(RichText::new("or a map open in SlopGene:").small().color(palette.muted));
        for doc in open {
            if ui.selectable_label(false, doc.title()).clicked() {
                chosen = Some(match &doc.path {
                    Some(path) => OpenParams { path: Some(protocol::stored_path(&folder, path)), document: None },
                    // A map with no file yet goes into the protocol itself.
                    None => OpenParams { path: None, document: Some(json::document_to_value(&doc.doc, true)) },
                });
            }
        }
    }
    chosen
}

/// The place to start an empty protocol, or to put the insert into a vector.
fn start_section(ui: &mut egui::Ui, tab: &ProtocolTab, docs: &[DocTab], palette: &Palette, edits: &mut Vec<Edit>, insert: Option<&str>) {
    match insert {
        None => {
            ui.label(RichText::new("Start with a plasmid").font(style::bold_font(15.5)).color(palette.ink));
            ui.label(
                RichText::new("Choose the map that carries the part you want to cut out. The protocol then cuts it and keeps the part.").color(palette.ink),
            );
            if let Some(source) = map_choices(ui, tab, docs, palette) {
                edits.push(Edit::AddSource(source));
            }
        }
        Some(what) => {
            ui.label(RichText::new("Put it into a plasmid").font(style::bold_font(15.5)).color(palette.ink));
            ui.label(
                RichText::new(format!("Choose the plasmid that {what} goes into. It is cut with the same enzymes, and the two are joined.")).color(palette.ink),
            );
            if let Some(vector) = map_choices(ui, tab, docs, palette) {
                edits.push(Edit::PutInto(vector));
            }
        }
    }
}

/// Where a map comes from, and a way to change it or to open it in a tab.
fn source_section(ui: &mut egui::Ui, tab: &ProtocolTab, index: usize, params: &OpenParams, docs: &[DocTab], palette: &Palette, edits: &mut Vec<Edit>) {
    let shown = match (&params.path, &params.document) {
        (Some(path), _) => file_name(path),
        (None, Some(_)) => "a map stored in the protocol".into(),
        (None, None) => "no map chosen yet".into(),
    };
    ui.label(RichText::new(format!("From {shown}")).small().color(palette.muted));
    ui.horizontal(|ui| {
        if let Some(made) = tab.run.steps.get(index).and_then(|r| r.outputs.first())
            && ui.small_button("Open as a tab").on_hover_text("Open this map as a sequence of its own").clicked()
        {
            edits.push(Edit::OpenDocument(Box::new(made.doc.clone())));
        }
        ui.menu_button("Use another map", |ui| {
            if let Some(source) = map_choices(ui, tab, docs, palette) {
                edits.push(change_step(index, None, move |s| s.action = Action::Open(source)));
                ui.close();
            }
        });
    });
    // Without a vector yet, the place to put the insert into one.
    if !tab.protocol.steps.iter().any(|s| matches!(s.action, Action::Ligate(_))) && tab.protocol.steps.iter().any(|s| matches!(s.action, Action::Pick(_))) {
        ui.separator();
        start_section(ui, tab, docs, palette, edits, Some("the kept piece"));
    }
}

/// The open document of the map a cut works on, when the map is open in a tab: its selection can
/// become the part.
fn open_tab_of<'a>(tab: &ProtocolTab, docs: &'a [DocTab], cut: usize) -> Option<&'a DocTab> {
    let source = tab.protocol.steps.get(cut)?.inputs.first().and_then(|r| source_of(&tab.protocol, r))?;
    let Action::Open(OpenParams { path: Some(path), .. }) = &tab.protocol.steps[source].action else { return None };
    open_tab(&tab.folder(), docs, path)
}

/// A cut's choices: the part it cuts around, the enzymes nearest the part on each side to pick one
/// from each, the enzymes as typed, and the pieces it made.
#[allow(clippy::too_many_arguments)]
fn cut_section(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    board: &Board,
    cut: usize,
    listed: Option<&NearestCuts>,
    docs: &[DocTab],
    palette: &Palette,
    edits: &mut Vec<Edit>,
) {
    let Some(Action::Digest(params)) = tab.protocol.steps.get(cut).map(|s| &s.action) else { return };
    let result = tab.run.steps.get(cut);
    let parent = result.and_then(|r| r.inputs.first());
    let part = params.around.clone();
    let role = protocol::part_role(&tab.protocol, &tab.protocol.steps[cut].id);

    // The part: a feature, a range, or the selection on the open map.
    let title = if role == Some(PartRole::Replace) { "The part the insert replaces" } else { "The part to cut out" };
    part_section(ui, tab, cut, parent, &part, result.and_then(|r| r.findings.region), title, docs, palette, edits);

    // The enzymes nearest the part on each side. One on each side cuts the part out.
    let partner = partner(&tab.protocol, &tab.run, cut);
    let other = other_plasmid(&tab.protocol, &tab.run, cut).map(|d| d.name);
    let (left, right) = flanks(result, &params.enzymes);
    if let (Some(listed), Some(part)) = (listed, &part) {
        section_title(ui, "Pick one enzyme on each side", palette);
        for (before, side) in [(true, &listed.before), (false, &listed.after)] {
            let word = if before { "Before" } else { "After" };
            ui.label(RichText::new(format!("{word} {part}, nearest first")).small().color(palette.muted));
            if side.is_empty() {
                ui.label(RichText::new("no enzyme of this set cuts there").small().color(palette.faint));
            }
            let chosen = if before { &left } else { &right };
            ui.horizontal_wrapped(|ui| {
                for near in side.iter() {
                    let on = chosen.as_deref().is_some_and(|c| c.eq_ignore_ascii_case(&near.enzyme));
                    let joins = partner.as_ref().map(|p| excise::can_join(near.end, &near.bases, &p.doc));
                    if chip(ui, near, on, joins, other.as_deref(), part, before, palette).clicked() {
                        edits.push(Edit::PickEnzyme(cut, before, near.enzyme.clone()));
                    }
                }
            });
        }
        if !listed.inside.is_empty() {
            ui.label(RichText::new(format!("These cut inside {part}, so they would cut it too: {}", listed.inside.join(", "))).small().color(palette.faint));
        }
    } else if part.is_none() {
        ui.label(
            RichText::new(
                "Mark the part first: drag across the plasmid's map or click a feature on it, or choose one above. The enzymes that cut nearest it are then listed here.",
            )
            .color(palette.muted),
        );
    }

    // The enzymes as typed, the set the list comes from, and the insert's enzymes for a vector.
    section_title(ui, "Enzymes", palette);
    if let Some(text) =
        committed_text(ui, egui::Id::unique(("protocol-enzymes", tab.id, cut)), &params.enzymes.join(", "), "or type them, such as EcoRI, BamHI")
    {
        let enzymes: Vec<String> = text.split([',', ' ', ';', '+']).filter(|s| !s.is_empty()).map(str::to_string).collect();
        edits.push(Edit::SetEnzymes(cut, enzymes));
    }
    set_chooser(ui, tab, palette, edits);
    if role == Some(PartRole::Replace)
        && let Some(insert_enzymes) = insert_cut(tab, cut).and_then(|c| match &tab.protocol.steps[c].action {
            Action::Digest(d) => Some(d.enzymes.clone()),
            _ => None,
        })
        && !insert_enzymes.is_empty()
        && insert_enzymes != params.enzymes
        && ui.small_button(format!("Use the insert's enzymes ({})", protocol::and_list(&insert_enzymes))).clicked()
    {
        edits.push(Edit::SetEnzymes(cut, insert_enzymes));
    }

    pieces_and_gel(ui, tab, board, cut, &params.enzymes.join("+"), palette, edits);
}

/// The part a cut is about or a PCR copies, to choose: a feature of the map, a range of its
/// bases, or the selection on the map when it is open in a tab. `setting` is the cut or the PCR,
/// `parent` the map as the run read it, and `region` the part as the run found it.
#[allow(clippy::too_many_arguments)]
fn part_section(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    setting: usize,
    parent: Option<&Piece>,
    part: &Option<String>,
    region: Option<Span>,
    title: &str,
    docs: &[DocTab],
    palette: &Palette,
    edits: &mut Vec<Edit>,
) {
    section_title(ui, title, palette);
    let Some(parent) = parent else { return };
    let mut names: Vec<String> = parent.doc.features.iter().map(|f| f.name.clone()).filter(|n| !n.is_empty()).collect();
    names.sort_by_key(|n| n.to_ascii_lowercase());
    names.dedup();
    let shown = part.clone().unwrap_or_else(|| "choose a feature…".into());
    egui::ComboBox::from_id_salt(("protocol-part", tab.id, setting)).selected_text(shown).width(ui.available_width() - 8.0).show_ui(ui, |ui| {
        for name in names {
            if ui.selectable_label(part.as_deref() == Some(name.as_str()), &name).clicked() {
                edits.push(Edit::SetPart(setting, Some(name)));
            }
        }
    });
    let range = part.clone().filter(|p| protocol::is_range(p)).unwrap_or_default();
    if let Some(text) = committed_text(ui, egui::Id::unique(("protocol-part-range", tab.id, setting)), &range, "or bases, such as 1201..1950") {
        let text = text.trim().to_string();
        edits.push(Edit::SetPart(setting, (!text.is_empty()).then_some(text)));
    }
    // A selection on the open map that is not the part already can become the part.
    if let Some(open) = open_tab_of(tab, docs, setting)
        && let Some(selection) = open.selection
        && region != Some(selection.normalize(open.len()))
    {
        let range = selection.format(open.len());
        if ui.small_button(format!("Use the selection on the map ({range})")).clicked() {
            edits.push(Edit::SetPart(setting, Some(range)));
        }
    }
}

/// The enzyme set the choosers list from.
fn set_chooser(ui: &mut egui::Ui, tab: &ProtocolTab, palette: &Palette, edits: &mut Vec<Edit>) {
    ui.horizontal(|ui| {
        ui.label(RichText::new("List enzymes from").small().color(palette.muted));
        egui::ComboBox::from_id_salt(("protocol-set", tab.id)).selected_text(tab.set.as_str()).width(110.0).show_ui(ui, |ui| {
            for name in ENZYME_SETS {
                if ui.selectable_label(tab.set == *name, *name).clicked() {
                    edits.push(Edit::Set(name.to_string()));
                }
            }
        });
    });
}

/// The pieces of cut `cut`, with the one kept, and the gel they make. `label` names the lane.
fn pieces_and_gel(ui: &mut egui::Ui, tab: &ProtocolTab, board: &Board, cut: usize, label: &str, palette: &Palette, edits: &mut Vec<Edit>) {
    let Some(result) = tab.run.steps.get(cut).filter(|r| r.status == Status::Done) else { return };
    section_title(ui, "The pieces", palette);
    if let Some(card) = board.card_of(cut) {
        piece_rows(ui, tab, board, card, palette, edits);
    }
    if let Some(input) = result.inputs.first() {
        let gel = gel_scene(input, &result.outputs, label);
        ui.add_space(4.0);
        let (rect, _) = ui.allocate_exact_size(vec2(gel.width, gel.height), Sense::hover());
        ui.painter().rect_filled(rect, CornerRadius::same(4), scene_paint::color(gel.background));
        paint(ui.painter(), &gel, rect.min);
        ui.label(RichText::new("At the bench, run the cut DNA on a gel and cut out the band of the piece kept.").small().color(palette.muted));
    }
}

/// The cut that trims a PCR product's tails. Its enzymes are the tails', chosen on the PCR step,
/// so it shows what it does and its pieces.
fn trim_section(ui: &mut egui::Ui, tab: &ProtocolTab, board: &Board, cut: usize, palette: &Palette, edits: &mut Vec<Edit>) {
    let Some(Action::Digest(params)) = tab.protocol.steps.get(cut).map(|s| &s.action) else { return };
    let copy = tab.protocol.amplify_before(&tab.protocol.steps[cut].id).and_then(|id| tab.protocol.steps.iter().position(|s| s.id == id));
    section_title(ui, "Trimming the tails", palette);
    let enzymes = if params.enzymes.is_empty() { "no enzymes yet".to_string() } else { protocol::and_list(&params.enzymes) };
    ui.label(
        RichText::new(format!("The cut trims the primers' tails with {enzymes}, so the piece gets the sticky ends the vector has. The enzymes are the tails', chosen on the PCR step."))
            .color(palette.ink),
    );
    if let Some(copy) = copy
        && ui.small_button("Choose the tail enzymes").clicked()
    {
        edits.push(Edit::Select(copy));
    }
    pieces_and_gel(ui, tab, board, cut, &params.enzymes.join("+"), palette, edits);
}

/// The enzymes of the vector's cut in the join that step `step`'s lane feeds, when the vector is
/// cut: the enzymes an insert's tails should match.
fn vector_cut_enzymes(tab: &ProtocolTab, step: usize) -> Option<Vec<String>> {
    match &tab.protocol.steps[insert_cut(tab, step)?].action {
        Action::Digest(d) => Some(d.enzymes.clone()),
        _ => None,
    }
}

/// A PCR's choices: the part it copies, the enzyme whose site goes on each primer's tail, the
/// enzymes as typed, the primers with a way to keep them or use another pair, and the product.
#[allow(clippy::too_many_arguments)]
fn copy_section(
    ui: &mut egui::Ui,
    tab: &ProtocolTab,
    copy: usize,
    tails: Option<&(usize, Vec<TailChoice>)>,
    docs: &[DocTab],
    palette: &Palette,
    edits: &mut Vec<Edit>,
) {
    let Some(Action::Amplify(params)) = tab.protocol.steps.get(copy).map(|s| &s.action) else { return };
    let result = tab.run.steps.get(copy);
    let parent = result.and_then(|r| r.inputs.first());
    let part = params.region.clone();
    part_section(ui, tab, copy, parent, &part, result.and_then(|r| r.findings.region), "The part to copy", docs, palette, edits);

    // The enzyme sites on the primers' tails: one enzyme for each primer, from those that do not
    // cut inside the part. The ones that open the vector once come first, as the ends must match.
    section_title(ui, "Enzyme sites on the primers' tails", palette);
    let listed = tails.filter(|(for_step, _)| *for_step == copy).map(|(_, t)| t.as_slice());
    let other = other_plasmid(&tab.protocol, &tab.run, copy).map(|d| d.name);
    match (listed, &part) {
        (Some(listed), Some(part)) => {
            for (forward, chosen) in [(true, &params.forward_enzyme), (false, &params.reverse_enzyme)] {
                let (which, end) = if forward { ("Forward primer", "start") } else { ("Reverse primer", "end") };
                ui.label(RichText::new(format!("{which}, at the {end} of {part}")).small().color(palette.muted));
                ui.horizontal_wrapped(|ui| {
                    for choice in listed.iter().filter(|c| !c.cuts_part) {
                        let on = chosen.as_deref().is_some_and(|c| c.eq_ignore_ascii_case(&choice.name));
                        if tail_chip(ui, choice, on, other.as_deref(), palette).clicked() {
                            edits.push(Edit::PickTail(copy, forward, choice.name.clone()));
                        }
                    }
                });
            }
            let inside: Vec<&str> = listed.iter().filter(|c| c.cuts_part).map(|c| c.name.as_str()).collect();
            if !inside.is_empty() {
                ui.label(RichText::new(format!("These cut inside {part}, so they are left out: {}", inside.join(", "))).small().color(palette.faint));
            }
        }
        _ => {
            ui.label(
                RichText::new(
                    "Mark the part first: drag across the map or click a feature on it, or choose one above. The enzymes that can go on the tails are then listed here.",
                )
                .color(palette.muted),
            );
        }
    }

    // The enzymes as typed, the set the list comes from, and the vector's enzymes.
    section_title(ui, "Enzymes", palette);
    if let Some(text) =
        committed_text(ui, egui::Id::unique(("protocol-enzymes", tab.id, copy)), &params.enzymes().join(", "), "or type them, such as EcoRI, BamHI")
    {
        let enzymes: Vec<String> = text.split([',', ' ', ';', '+']).filter(|s| !s.is_empty()).map(str::to_string).collect();
        edits.push(Edit::SetEnzymes(copy, enzymes));
    }
    set_chooser(ui, tab, palette, edits);
    if let Some(vector_enzymes) = vector_cut_enzymes(tab, copy)
        && !vector_enzymes.is_empty()
        && vector_enzymes != params.enzymes()
        && ui.small_button(format!("Use the vector's enzymes ({})", protocol::and_list(&vector_enzymes))).clicked()
    {
        edits.push(Edit::SetEnzymes(copy, vector_enzymes));
    }

    // The primers: each with its sequence to copy out, and the way to keep them, to take another
    // pair or to design them anew. Then the product.
    let Some(result) = result.filter(|r| r.status == Status::Done) else { return };
    section_title(ui, "Primers", palette);
    for p in &result.findings.primers {
        ui.horizontal(|ui| {
            ui.label(RichText::new(&p.name).font(style::bold_font(12.5)).color(palette.ink));
            ui.label(RichText::new(&p.description).small().color(palette.muted));
            if ui.small_button("Copy").on_hover_text("Copy the name and the sequence").clicked() {
                ui.ctx().copy_text(format!("{}\t{}", p.name, p.seq));
            }
        });
        ui.add(egui::Label::new(RichText::new(&p.seq).monospace().size(11.5).color(palette.ink)).wrap());
    }
    if params.pinned() {
        ui.label(RichText::new("Kept in the protocol, so they stay as they are when the maps change.").small().color(palette.muted));
        if ui.small_button("Design the primers anew").clicked() {
            edits.push(Edit::Primers(copy, PrimerEdit::Clear));
        }
    } else {
        ui.label(RichText::new("Designed anew on each run. Keep them before ordering them.").small().color(palette.muted));
        if ui.small_button("Keep these primers").clicked() {
            edits.push(Edit::Primers(copy, PrimerEdit::Pair(0)));
        }
        let pairs = &result.findings.primer_pairs;
        if pairs.len() > 1 {
            egui::CollapsingHeader::new(format!("{} other pairs", pairs.len() - 1)).id_salt(("protocol-other-pairs", tab.id, copy)).show(ui, |ui| {
                for (ii, pair) in pairs.iter().enumerate().skip(1) {
                    ui.horizontal_wrapped(|ui| {
                        if ui.small_button(format!("Use pair {}", ii + 1)).clicked() {
                            edits.push(Edit::Primers(copy, PrimerEdit::Pair(ii)));
                        }
                        let mut words = format!(
                            "Tm {:.1} and {:.1} C, {} and {} nt anneal, {} bp product",
                            pair.forward_tm,
                            pair.reverse_tm,
                            pair.forward_anneals,
                            pair.reverse_anneals,
                            thousands(pair.product_len)
                        );
                        if !pair.warnings.is_empty() {
                            words.push_str(&format!(" ({})", pair.warnings.join("; ")));
                        }
                        ui.label(RichText::new(words).small().color(palette.muted));
                    });
                }
            });
        }
    }
    if let Some(product) = result.outputs.first() {
        section_title(ui, "The product", palette);
        ui.label(RichText::new(format!("{} bp, with the tails on its ends and no 5' phosphates.", thousands(product.doc.len()))).color(palette.ink));
        if ui.small_button("Open the PCR product as a tab").clicked() {
            edits.push(Edit::OpenDocument(Box::new(product.doc.clone())));
        }
    }
}

/// One enzyme to put on a primer's tail: its name, chosen or not, with a tick when it opens the
/// vector once, so the ends can match. Hovering tells the rest.
fn tail_chip(ui: &mut egui::Ui, choice: &TailChoice, chosen: bool, other: Option<&str>, palette: &Palette) -> egui::Response {
    let mut job = egui::text::LayoutJob::default();
    let format = |color: Color32, font: egui::FontId| egui::TextFormat { font_id: font, color, ..Default::default() };
    job.append(&choice.name, 0.0, format(palette.ink, style::bold_font(12.5)));
    if choice.vector_cuts == Some(1) {
        job.append("  ✓", 0.0, format(palette.ok, style::bold_font(12.0)));
    }
    let mut lines = vec![format!("Puts the {} site on this primer's tail.", choice.name)];
    if let (Some(n), Some(other)) = (choice.vector_cuts, other) {
        lines.push(format!("Cuts {other} {n} time{}.", if n == 1 { "" } else { "s" }));
    }
    if choice.vector_cuts == Some(1) {
        lines.push("It opens the vector once, so the ends can match.".into());
    }
    ui.selectable_label(chosen, job).on_hover_text(lines.join("\n"))
}

/// The cut of the insert that goes into the same join as cut `cut`: for a vector's cut, the
/// insert's.
fn insert_cut(tab: &ProtocolTab, cut: usize) -> Option<usize> {
    let (join, mine) = feeds(&tab.protocol, cut)?;
    tab.protocol.steps[join].inputs.iter().enumerate().filter(|(ii, _)| *ii != mine).find_map(|(_, input)| {
        // Follow the input back to the cut it came from.
        let mut id = split_reference(input).0.to_string();
        for _ in 0..32 {
            let index = tab.protocol.steps.iter().position(|s| s.id == id)?;
            if matches!(tab.protocol.steps[index].action, Action::Digest(_)) {
                return Some(index);
            }
            id = split_reference(tab.protocol.steps[index].inputs.first()?).0.to_string();
        }
        None
    })
}

/// One enzyme to pick: its name and distance, chosen or not, with a tick when its end can join the
/// piece it is ligated to. Hovering tells the rest.
#[allow(clippy::too_many_arguments)]
fn chip(
    ui: &mut egui::Ui,
    near: &NearCut,
    chosen: bool,
    joins: Option<bool>,
    other: Option<&str>,
    part: &str,
    before: bool,
    palette: &Palette,
) -> egui::Response {
    let mut job = egui::text::LayoutJob::default();
    let format = |color: Color32, font: egui::FontId| egui::TextFormat { font_id: font, color, ..Default::default() };
    job.append(&near.enzyme, 0.0, format(palette.ink, style::bold_font(12.5)));
    job.append(&format!("  {} bp", thousands(near.distance)), 0.0, format(palette.muted, style::ui_font(12.0)));
    if joins == Some(true) {
        job.append("  ✓", 0.0, format(palette.ok, style::bold_font(12.0)));
    }
    let mut lines = vec![
        format!("Cuts {} bp {} {part}, after base {}.", thousands(near.distance), if before { "before" } else { "after" }, thousands(near.pos)),
        format!("Leaves a {} end.", near.end_text()),
        format!("Cuts this map {} time{}.", near.cuts, if near.cuts == 1 { "" } else { "s" }),
    ];
    if let (Some(n), Some(other)) = (near.vector_cuts, other) {
        lines.push(format!("Cuts {other} {n} time{}.", if n == 1 { "" } else { "s" }));
    }
    match joins {
        Some(true) => lines.push("Its end can join the piece it is ligated to.".into()),
        Some(false) => lines.push("Its end does not pair with the piece it is ligated to.".into()),
        None => {}
    }
    ui.selectable_label(chosen, job).on_hover_text(lines.join("\n"))
}

/// The pieces of a cut card: each with its length, ends and features, the one kept marked, and a
/// button that keeps another.
fn piece_rows(ui: &mut egui::Ui, tab: &ProtocolTab, board: &Board, card: usize, palette: &Palette, edits: &mut Vec<Edit>) {
    let c = &board.cards[card];
    let Some(view) = cut_view(tab, c) else { return };
    let Some(result) = c.step.and_then(|s| tab.run.steps.get(s)) else { return };
    // Kept pieces first, then the rest from the largest down.
    let mut order: Vec<&(Span, Role, Option<usize>)> = view.pieces.iter().collect();
    order.sort_by_key(|(span, role, _)| (*role != Role::Kept, std::cmp::Reverse(span.len())));
    for (span, role, output) in order {
        let Some(piece) = output.and_then(|o| result.outputs.get(o)) else { continue };
        ui.horizontal(|ui| {
            let (word, color) = match role {
                Role::Kept => ("kept", palette.ok),
                Role::LeftOver => ("left over", palette.error),
                Role::Plain => ("not chosen", palette.warn),
            };
            ui.label(RichText::new(format!("{} bp", thousands(piece.doc.len()))).font(style::bold_font(12.5)).color(palette.ink));
            ui.label(RichText::new(word).small().color(color));
            // Its ends by the enzymes that made them, and the part it holds, or else the largest
            // feature it carries.
            let ends = format!("{}–{}", piece.left_enzyme.as_deref().unwrap_or("end"), piece.right_enzyme.as_deref().unwrap_or("end"));
            ui.label(RichText::new(ends).small().color(palette.muted));
            let part = match (&tab.protocol.steps[c.step.unwrap_or_default()].action, view.region) {
                (Action::Digest(d), Some(region)) if span.contains_span(&region, view.parent.len()) => d.around.clone(),
                _ => None,
            };
            let largest = piece.doc.features.iter().filter(|f| !f.name.is_empty() && f.len() >= 30).max_by_key(|f| f.len()).map(|f| f.name.clone());
            if let Some(name) = part.or(largest) {
                ui.label(RichText::new(name).small().color(palette.faint));
            }
            if *role != Role::Kept
                && let (Some(&pick), Some(output)) = (c.picks.first(), output)
            {
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.small_button("Keep").on_hover_text(format!("Keep bases {} instead", span.format(view.parent.len()))).clicked() {
                        edits.push(Edit::Keep(pick, *output));
                    }
                });
            }
        });
    }
}

/// A pick's choices: the pieces of its cut, as on the cut's card.
fn pieces_section(ui: &mut egui::Ui, tab: &ProtocolTab, board: &Board, index: usize, palette: &Palette, edits: &mut Vec<Edit>) {
    let Some(card) = board.card_of(index).filter(|&c| board.cards[c].kind == Kind::Cut) else { return };
    section_title(ui, "The pieces", palette);
    piece_rows(ui, tab, board, card, palette, edits);
}

/// A join's junctions in plain words, its options when the direction waits for a choice, the
/// other plasmids its tube makes, and a way to open the new plasmid.
fn join_section(ui: &mut egui::Ui, tab: &ProtocolTab, index: usize, palette: &Palette, edits: &mut Vec<Edit>) {
    let Some(result) = tab.run.steps.get(index) else { return };
    if result.status == Status::NeedsChoice {
        section_title(ui, "Choose the direction", palette);
        for (ii, choice) in result.choices.iter().enumerate() {
            ui.horizontal_wrapped(|ui| {
                if ui.button(format!("Use option {}", ii + 1)).clicked() {
                    edits.push(Edit::Choose(index, ii));
                }
                ui.label(RichText::new(&choice.label).small().color(palette.muted));
            });
        }
    }
    // Each other plasmid of the tube, with the way to make it the product or to stop it forming.
    if !result.findings.byproducts.is_empty() {
        section_title(ui, "Also forms in the same tube", palette);
        let backbone = tab.protocol.steps.get(index).and_then(|s| s.inputs.first()).map(|r| protocol::piece_name(&tab.protocol, &tab.run, r));
        for (ii, other) in result.findings.byproducts.iter().enumerate() {
            let words = uppercase_first(&protocol::byproduct_words(&tab.protocol, &tab.run, index, other));
            ui.label(RichText::new(format!("{words} ({})", other.piece.summary())).color(palette.ink));
            ui.horizontal_wrapped(|ui| {
                if other.pin.is_some() && ui.small_button("Use this instead").clicked() {
                    edits.push(Edit::UseInstead(index, ii));
                }
                if other.alone()
                    && let Some(backbone) = &backbone
                    && ui.small_button(format!("Remove the 5' phosphates of {backbone}")).clicked()
                {
                    edits.push(Edit::Dephosphorylate(index, 0));
                }
            });
        }
    }
    if let Some(product) = result.outputs.first()
        && ui.button("Open the new plasmid as a tab").on_hover_text("Open it as a sequence of its own, to edit or save it").clicked()
    {
        let mut doc = product.doc.clone();
        doc.name = tab.name();
        edits.push(Edit::OpenDocument(Box::new(doc)));
    }
}

/// The editor of the step shown, by hand: its id, inputs, parameters and note, and buttons to
/// move or delete it.
fn step_editor(ui: &mut egui::Ui, tab: &ProtocolTab, palette: &Palette, edits: &mut Vec<Edit>) {
    let index = tab.selected;
    let Some(step) = tab.protocol.steps.get(index) else { return };
    // The id names the step's output for the steps after it.
    field_label(ui, "Name", palette);
    if let Some(id) = committed_text(ui, egui::Id::unique(("protocol-step-id", tab.id, index)), &step.id, "a short name") {
        edits.push(Edit::Rename(step.id.clone(), id));
    }
    input_editor(ui, tab, index, palette, edits);
    match &step.action {
        Action::Digest(params) => digest_editor(ui, tab, index, params, palette, edits),
        Action::Amplify(params) => amplify_editor(ui, tab, index, params, palette, edits),
        Action::Pick(params) => pick_editor(ui, tab, index, params, palette, edits),
        Action::Ligate(params) => ligate_editor(ui, index, params, edits),
        Action::Open(_) | Action::Dephosphorylate(_) | Action::Flip(_) => {}
    }
    field_label(ui, "Note", palette);
    if let Some(note) = committed_multiline(ui, egui::Id::unique(("protocol-step-note", tab.id, index)), &step.note) {
        edits.push(change_step(index, None, move |s| s.note = note));
    }
    ui.add_space(6.0);
    ui.horizontal(|ui| {
        if ui.add_enabled(index > 0, egui::Button::new("Move up")).clicked() {
            edits.push(Edit::MoveStep(index, true));
        }
        if ui.add_enabled(index + 1 < tab.protocol.steps.len(), egui::Button::new("Move down")).clicked() {
            edits.push(Edit::MoveStep(index, false));
        }
        if ui.button("Delete step").clicked() {
            edits.push(Edit::RemoveStep(index));
        }
    });
}

fn field_label(ui: &mut egui::Ui, text: &str, palette: &Palette) {
    ui.add_space(4.0);
    ui.label(RichText::new(text).small().color(palette.muted));
}

/// The inputs of the step: one drop-down for a step that works on one molecule, and for a
/// ligation one per piece, with how each piece after the first goes in.
fn input_editor(ui: &mut egui::Ui, tab: &ProtocolTab, index: usize, palette: &Palette, edits: &mut Vec<Edit>) {
    let step = &tab.protocol.steps[index];
    let references = tab.protocol.references_before(index, Some(&tab.run));
    if step.action.takes_one_input() {
        field_label(ui, "Works on", palette);
        let current = step.inputs.first().cloned().unwrap_or_default();
        if let Some(picked) = reference_combo(ui, ("protocol-input", tab.id, index, 0), &current, &references, "choose…") {
            edits.push(change_step(index, None, move |s| s.inputs = vec![picked]));
        }
        return;
    }
    let Action::Ligate(params) = &step.action else { return };
    field_label(ui, "Joins, in this order", palette);
    for (ii, input) in step.inputs.iter().enumerate() {
        ui.horizontal(|ui| {
            if let Some(picked) = reference_combo(ui, ("protocol-input", tab.id, index, ii), input, &references, "choose…") {
                edits.push(change_step(index, None, move |s| {
                    if let Some(slot) = s.inputs.get_mut(ii) {
                        *slot = picked;
                    }
                }));
            }
            // Orientation ii - 1 belongs to input ii: the first piece sets the direction.
            if ii > 0 {
                let current = params.orientations.get(ii - 1).copied().unwrap_or_default();
                let words = |o: Orientation| match o {
                    Orientation::Auto => "either way",
                    Orientation::Forward => "as given",
                    Orientation::Reverse => "turned round",
                };
                egui::ComboBox::from_id_salt(("protocol-orientation", tab.id, index, ii)).selected_text(words(current)).width(96.0).show_ui(ui, |ui| {
                    for option in [Orientation::Auto, Orientation::Forward, Orientation::Reverse] {
                        if ui.selectable_label(option == current, words(option)).clicked() && option != current {
                            edits.push(change_step(index, None, move |s| {
                                if let Action::Ligate(p) = &mut s.action {
                                    if p.orientations.len() < ii {
                                        p.orientations.resize(ii, Orientation::Auto);
                                    }
                                    p.orientations[ii - 1] = option;
                                }
                            }));
                        }
                    }
                });
            }
            if ui.small_button("×").on_hover_text("Take this piece out").clicked() {
                edits.push(change_step(index, None, move |s| {
                    if ii < s.inputs.len() {
                        s.inputs.remove(ii);
                    }
                    // The orientation of the piece after it moves up with it.
                    if let Action::Ligate(p) = &mut s.action {
                        let dropped = ii.saturating_sub(1);
                        if dropped < p.orientations.len() {
                            p.orientations.remove(dropped);
                        }
                    }
                }));
            }
        });
    }
    if let Some(picked) = reference_combo(ui, ("protocol-input-add", tab.id, index), "", &references, "add a piece…") {
        edits.push(change_step(index, None, move |s| s.inputs.push(picked)));
    }
}

/// A drop-down of the references a step can take. Returns the one picked when it differs from
/// `current`.
fn reference_combo(ui: &mut egui::Ui, id: impl std::hash::Hash + std::fmt::Debug, current: &str, references: &[String], empty: &str) -> Option<String> {
    let mut picked = None;
    egui::ComboBox::from_id_salt(id).selected_text(if current.is_empty() { empty } else { current }).width(150.0).show_ui(ui, |ui| {
        if references.is_empty() {
            ui.label("no earlier steps");
        }
        for reference in references {
            if ui.selectable_label(reference == current, reference).clicked() && reference != current {
                picked = Some(reference.clone());
            }
        }
    });
    picked
}

/// The part a digest is about, as typed.
fn digest_editor(ui: &mut egui::Ui, tab: &ProtocolTab, index: usize, params: &DigestParams, palette: &Palette, edits: &mut Vec<Edit>) {
    field_label(ui, "Cut around (a feature or a range)", palette);
    if let Some(text) = committed_text(ui, egui::Id::unique(("protocol-around", tab.id, index)), params.around.as_deref().unwrap_or(""), "EGFP, or 389..1108") {
        let text = text.trim().to_string();
        edits.push(Edit::SetPart(index, (!text.is_empty()).then_some(text)));
    }
}

/// What a PCR copies and with what, as typed: the part, the enzyme on each primer's tail, and the
/// primers themselves, which are designed when left empty.
fn amplify_editor(ui: &mut egui::Ui, tab: &ProtocolTab, index: usize, params: &AmplifyParams, palette: &Palette, edits: &mut Vec<Edit>) {
    field_label(ui, "Copy (a feature or a range)", palette);
    if let Some(text) = committed_text(ui, egui::Id::unique(("protocol-region", tab.id, index)), params.region.as_deref().unwrap_or(""), "EGFP, or 389..1108") {
        let text = text.trim().to_string();
        edits.push(Edit::SetPart(index, (!text.is_empty()).then_some(text)));
    }
    for (forward, label, current) in
        [(true, "Forward primer's tail enzyme", &params.forward_enzyme), (false, "Reverse primer's tail enzyme", &params.reverse_enzyme)]
    {
        field_label(ui, label, palette);
        if let Some(text) = committed_text(ui, egui::Id::unique(("protocol-tail", tab.id, index, forward)), current.as_deref().unwrap_or(""), "none, or EcoRI")
        {
            // An empty name takes the tail off, which choosing the enzyme there again does.
            let name = text.trim().to_string();
            match (name.is_empty(), current) {
                (true, Some(old)) => edits.push(Edit::PickTail(index, forward, old.clone())),
                (true, None) => {}
                (false, _) => edits.push(Edit::PickTail(index, forward, name)),
            }
        }
    }
    field_label(ui, "Primers, 5' to 3' with their tails (empty to design them)", palette);
    let forward =
        committed_text(ui, egui::Id::unique(("protocol-forward-primer", tab.id, index)), params.forward_primer.as_deref().unwrap_or(""), "forward primer");
    let reverse =
        committed_text(ui, egui::Id::unique(("protocol-reverse-primer", tab.id, index)), params.reverse_primer.as_deref().unwrap_or(""), "reverse primer");
    // Each field commits on its own, and the other keeps its value.
    if forward.is_some() || reverse.is_some() {
        let f = forward.unwrap_or_else(|| params.forward_primer.clone().unwrap_or_default()).trim().to_string();
        let r = reverse.unwrap_or_else(|| params.reverse_primer.clone().unwrap_or_default()).trim().to_string();
        edits.push(Edit::Primers(index, if f.is_empty() && r.is_empty() { PrimerEdit::Clear } else { PrimerEdit::Given(f, r) }));
    }
}

/// What a pick looks for. Every criterion set must hold.
fn pick_editor(ui: &mut egui::Ui, tab: &ProtocolTab, index: usize, params: &PickParams, palette: &Palette, edits: &mut Vec<Edit>) {
    let candidates: &[Piece] = tab.run.steps.get(index).map_or(&[], |r| r.inputs.as_slice());
    let mut edited = params.clone();
    field_label(ui, "Keep the fragment that", palette);
    ui.horizontal(|ui| {
        let mut on = edited.index.is_some();
        if ui.checkbox(&mut on, "is number").changed() {
            edited.index = on.then_some(edited.index.unwrap_or(1));
        }
        if let Some(number) = edited.index.as_mut() {
            ui.add(egui::DragValue::new(number).range(1..=candidates.len().max(1)));
        }
    });
    ui.horizontal(|ui| {
        let mut on = edited.size.is_some();
        if ui.checkbox(&mut on, "is about").changed() {
            edited.size = on.then(|| edited.size.unwrap_or_else(|| candidates.first().map_or(1000, |c| c.doc.len())));
        }
        if let Some(size) = edited.size.as_mut() {
            ui.add(egui::DragValue::new(size).range(1..=100_000_000).speed(10.0).suffix(" bp"));
        }
    });
    ui.horizontal(|ui| {
        ui.label("has the ends");
        if let Some(text) =
            committed_text(ui, egui::Id::unique(("protocol-pick-ends", tab.id, index)), edited.ends.as_deref().unwrap_or(""), "any, or EcoRI/BamHI")
        {
            edited.ends = (!text.trim().is_empty()).then(|| text.trim().to_string());
        }
    });
    // The feature names found among the candidates, to choose from.
    let mut names: Vec<&str> = candidates.iter().flat_map(|c| c.doc.features.iter().map(|f| f.name.as_str())).collect();
    names.sort_unstable();
    names.dedup();
    ui.horizontal(|ui| {
        ui.label("holds");
        let current = edited.contains.clone().unwrap_or_default();
        egui::ComboBox::from_id_salt(("protocol-pick-holds", tab.id, index))
            .selected_text(if current.is_empty() { "any feature" } else { &current })
            .width(150.0)
            .show_ui(ui, |ui| {
                if ui.selectable_label(current.is_empty(), "any feature").clicked() {
                    edited.contains = None;
                }
                for name in &names {
                    if ui.selectable_label(current == *name, *name).clicked() {
                        edited.contains = Some(name.to_string());
                    }
                }
            });
    });
    ui.horizontal(|ui| {
        ui.label("holds bases");
        if let Some(text) =
            committed_text(ui, egui::Id::unique(("protocol-pick-bases", tab.id, index)), edited.holds.as_deref().unwrap_or(""), "any, or 389..1108")
        {
            edited.holds = (!text.trim().is_empty()).then(|| text.trim().to_string());
        }
    });
    ui.horizontal(|ui| {
        if ui.checkbox(&mut edited.largest, "is the largest").changed() && edited.largest {
            edited.smallest = false;
        }
        if ui.checkbox(&mut edited.smallest, "is the smallest").changed() && edited.smallest {
            edited.largest = false;
        }
    });
    if edited != *params {
        edits.push(change_step(index, Some("pick"), move |s| s.action = Action::Pick(edited)));
    }
}

/// Whether a ligation closes its product into a circle. The order and direction of its pieces
/// are set with its inputs.
fn ligate_editor(ui: &mut egui::Ui, index: usize, params: &LigateParams, edits: &mut Vec<Edit>) {
    let mut circular = params.circular;
    if ui.checkbox(&mut circular, "Close the product into a circle").changed() {
        edits.push(change_step(index, None, move |s| {
            if let Action::Ligate(p) = &mut s.action {
                p.circular = circular;
            }
        }));
    }
}

/// A one-line text field that edits a copy of `value` and returns the new text once the person
/// presses Enter or leaves the field, so the protocol does not run again at every key.
pub(crate) fn committed_text(ui: &mut egui::Ui, id: egui::Id, value: &str, hint: &str) -> Option<String> {
    let mut text: String = ui.data_mut(|d| d.get_temp::<String>(id)).unwrap_or_else(|| value.to_string());
    let response = ui.add(egui::TextEdit::singleline(&mut text).id(id).hint_text(hint).desired_width(f32::INFINITY));
    if response.has_focus() {
        ui.data_mut(|d| d.insert_temp(id, text));
        return None;
    }
    // Away from the field the copy is dropped, so a change from elsewhere (undo) shows.
    ui.data_mut(|d| d.remove::<String>(id));
    (response.lost_focus() && text != value).then_some(text)
}

/// A text box for several lines that, like [`committed_text`], hands back its text once the
/// person leaves it.
fn committed_multiline(ui: &mut egui::Ui, id: egui::Id, value: &str) -> Option<String> {
    let mut text: String = ui.data_mut(|d| d.get_temp::<String>(id)).unwrap_or_else(|| value.to_string());
    let response = ui.add(egui::TextEdit::multiline(&mut text).id(id).hint_text("what to watch for at the bench").desired_rows(2).desired_width(f32::INFINITY));
    if response.has_focus() {
        ui.data_mut(|d| d.insert_temp(id, text));
        return None;
    }
    ui.data_mut(|d| d.remove::<String>(id));
    (response.lost_focus() && text != value).then_some(text)
}

#[cfg(test)]
mod tests {
    use super::*;
    use slopgene_core::protocol::CloningTemplate;

    /// Bases of A, C and T only, so no G-containing site (EcoRI, BamHI) appears by chance.
    fn background(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACT"[(x >> 33) as usize % 3]
            })
            .collect()
    }

    /// The maps of the tests: a vector with EcoRI and BamHI sites 60 bp apart, and a donor with
    /// geneX between the same two sites.
    fn maps() -> (Document, Document) {
        let mut vector = background(3000, 1);
        vector[1000..1006].copy_from_slice(b"GAATTC");
        vector[1060..1066].copy_from_slice(b"GGATCC");
        let mut donor = background(1500, 2);
        donor[100..106].copy_from_slice(b"GAATTC");
        donor[900..906].copy_from_slice(b"GGATCC");
        let mut donor = Document::new_dna("donor", donor, Topology::Circular);
        donor.features.push(Feature::new("geneX", "CDS", Span::new(200, 800), Strand::Forward));
        (Document::new_dna("vector", vector, Topology::Circular), donor)
    }

    /// A restriction cloning of geneX into the vector, both maps stored in the protocol itself.
    /// With `enzymes` empty the cuts wait for them.
    fn cloning(enzymes: &[&str], part: Option<&str>) -> Protocol {
        let (vector, donor) = maps();
        let stored = |doc: &Document| OpenParams { path: None, document: Some(json::document_to_value(doc, true)) };
        Protocol::restriction_cloning(&CloningTemplate {
            name: "geneX into vector".into(),
            vector: Some(stored(&vector)),
            donor: stored(&donor),
            insert_enzymes: enzymes.iter().map(|s| s.to_string()).collect(),
            insert: part.map(str::to_string),
            dephosphorylate: true,
            ..CloningTemplate::default()
        })
    }

    /// The same cloning with geneX copied out of the donor by PCR, with the enzymes' sites on
    /// the primers' tails.
    fn cloning_by_pcr(enzymes: &[&str], part: Option<&str>) -> Protocol {
        let (vector, donor) = maps();
        let stored = |doc: &Document| OpenParams { path: None, document: Some(json::document_to_value(doc, true)) };
        Protocol::restriction_cloning(&CloningTemplate {
            name: "geneX into vector by PCR".into(),
            vector: Some(stored(&vector)),
            donor: stored(&donor),
            insert_enzymes: enzymes.iter().map(|s| s.to_string()).collect(),
            insert: part.map(str::to_string),
            dephosphorylate: true,
            pcr: true,
            ..CloningTemplate::default()
        })
    }

    #[test]
    fn a_pick_waits_for_a_choice_and_choosing_settles_it() {
        let mut tab = ProtocolTab::new(cloning(&["EcoRI", "BamHI"], None), None);
        tab.refresh(&[], false);
        // Both donor fragments have an EcoRI end and a BamHI end, so the insert waits.
        let insert = tab.protocol.steps.iter().position(|s| s.id == "insert").unwrap();
        assert_eq!(tab.run.steps[insert].status, Status::NeedsChoice);
        assert_eq!(tab.progress(), "6 of 8 steps done, 1 waits for a choice, 1 waits for an earlier step");
        // Choosing the fragment that holds geneX writes it into the protocol, and the product forms.
        let option = tab.run.steps[insert].choices.iter().position(|c| c.outputs[0].doc.features.iter().any(|f| f.name == "geneX")).unwrap();
        tab.choose(insert, option);
        tab.refresh(&[], false);
        let product = tab.run.step("product").unwrap();
        assert_eq!(product.status, Status::Done, "{:?}", product.error);
        assert!(product.outputs[0].doc.features.iter().any(|f| f.name == "geneX"));
        assert!(tab.dirty);
        // Undo takes the choice back.
        assert!(tab.undo());
        tab.refresh(&[], false);
        assert_eq!(tab.run.steps[insert].status, Status::NeedsChoice);
    }

    #[test]
    fn the_bench_puts_each_map_in_a_column_and_the_join_between_them() {
        let protocol = cloning(&["EcoRI", "BamHI"], Some("geneX"));
        let nothing_run = Run { steps: Vec::new() };
        let laid = board(&protocol, &nothing_run);
        let kinds: Vec<(Kind, f32, usize)> = laid.cards.iter().map(|c| (c.kind, c.column, c.row)).collect();
        // Once a join exists, a place to add another map ends the top row.
        let expected =
            vec![(Kind::Plasmid, 0.0, 0), (Kind::Cut, 0.0, 1), (Kind::Plasmid, 1.0, 0), (Kind::Cut, 1.0, 1), (Kind::Join, 0.5, 2), (Kind::AddPlasmid, 2.0, 0)];
        assert_eq!(kinds, expected);
        // The picks fold into their cuts, and the dephosphorylation rides on the vector's arrow.
        let step = |id: &str| protocol.steps.iter().position(|s| s.id == id).unwrap();
        assert_eq!(laid.cards[1].picks, vec![step("backbone")]);
        assert_eq!(laid.cards[3].picks, vec![step("insert")]);
        let into_join: Vec<&Link> = laid.links.iter().filter(|l| l.to == 4).collect();
        assert_eq!(into_join.len(), 2);
        assert_eq!(into_join[0].treatments, vec![step("cip")]);
        assert_eq!(laid.card_of(step("cip")), Some(1));
        // Only the insert's steps: a place for the vector waits beside them.
        let mut only =
            Protocol::restriction_cloning(&CloningTemplate { donor: OpenParams::default(), insert: Some("geneX".into()), ..CloningTemplate::default() });
        assert_eq!(board(&only, &nothing_run).cards.last().map(|c| c.kind), Some(Kind::PutInto));
        only.steps.clear();
        assert_eq!(board(&only, &nothing_run).cards.iter().map(|c| c.kind).collect::<Vec<_>>(), vec![Kind::Start]);
    }

    #[test]
    fn picking_an_enzyme_on_each_side_cuts_the_part_out_and_the_vector_follows() {
        let mut tab = ProtocolTab::new(cloning(&[], Some("geneX")), None);
        tab.refresh(&[], false);
        let cut = tab.protocol.steps.iter().position(|s| s.id == "cut_donor").unwrap();
        assert_eq!(tab.run.steps[cut].status, Status::NeedsChoice, "the cut waits for its enzymes");
        // The chooser lists EcoRI before geneX and BamHI after it.
        let listed = tab.nearest(cut).cloned().unwrap();
        assert_eq!(listed.before[0].enzyme, "EcoRI");
        assert_eq!(listed.after[0].enzyme, "BamHI");
        assert_eq!(listed.before[0].vector_cuts, Some(1), "EcoRI cuts the vector once");
        apply(&mut tab, vec![Edit::PickEnzyme(cut, true, "EcoRI".into())]);
        tab.refresh(&[], false);
        apply(&mut tab, vec![Edit::PickEnzyme(cut, false, "BamHI".into())]);
        tab.refresh(&[], false);
        let enzymes = |tab: &ProtocolTab, id: &str| match &tab.protocol.step(id).unwrap().action {
            Action::Digest(d) => d.enzymes.clone(),
            _ => Vec::new(),
        };
        assert_eq!(enzymes(&tab, "cut_donor"), vec!["EcoRI".to_string(), "BamHI".to_string()]);
        assert_eq!(enzymes(&tab, "cut_vector"), vec!["EcoRI".to_string(), "BamHI".to_string()], "the vector is cut the same way");
        for step in &tab.run.steps {
            assert_eq!(step.status, Status::Done, "{}: {:?}", step.id, step.error);
        }
        // Picking the same enzyme again takes it off that side.
        apply(&mut tab, vec![Edit::PickEnzyme(cut, false, "BamHI".into())]);
        assert_eq!(enzymes(&tab, "cut_donor"), vec!["EcoRI".to_string()]);
    }

    #[test]
    fn choosing_a_side_replaces_the_enzyme_that_made_that_end() {
        let names = |list: &[&str]| list.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        let side = |name: &str| Some(name.to_string());
        // Another enzyme before the part takes the place of the one there, and comes first.
        assert_eq!(choose_side(&names(&["EcoRI", "BamHI"]), side("EcoRI"), side("BamHI"), true, "KpnI"), names(&["KpnI", "BamHI"]));
        // After the part it takes the other place.
        assert_eq!(choose_side(&names(&["EcoRI", "BamHI"]), side("EcoRI"), side("BamHI"), false, "XbaI"), names(&["EcoRI", "XbaI"]));
        // An enzyme that makes both ends stays for the other side when another takes this one.
        assert_eq!(choose_side(&names(&["EcoRI"]), side("EcoRI"), side("EcoRI"), false, "BamHI"), names(&["EcoRI", "BamHI"]));
        // Choosing a chosen enzyme again takes it off, like a switch.
        assert_eq!(choose_side(&names(&["EcoRI"]), side("EcoRI"), side("EcoRI"), true, "EcoRI"), Vec::<String>::new());
        // The first choice of all.
        assert_eq!(choose_side(&[], None, None, true, "EcoRI"), names(&["EcoRI"]));
    }

    #[test]
    fn clicking_a_piece_keeps_it_and_the_part_still_names_it() {
        let mut tab = ProtocolTab::new(cloning(&["EcoRI", "BamHI"], Some("geneX")), None);
        tab.refresh(&[], false);
        let cut = tab.protocol.steps.iter().position(|s| s.id == "cut_donor").unwrap();
        let insert = tab.protocol.steps.iter().position(|s| s.id == "insert").unwrap();
        let holder = tab.run.steps[cut].findings.holder.unwrap();
        // Keeping the other piece pins it, and the pick no longer asks for geneX.
        let other = 1 - holder;
        apply(&mut tab, vec![Edit::Keep(insert, other)]);
        let Action::Pick(params) = &tab.protocol.steps[insert].action else { panic!() };
        assert_eq!((params.index, params.contains.as_deref()), (Some(other + 1), None));
        // Keeping geneX's piece again names geneX, so the pick finds it after a change of enzymes.
        apply(&mut tab, vec![Edit::Keep(insert, holder)]);
        let Action::Pick(params) = &tab.protocol.steps[insert].action else { panic!() };
        assert_eq!((params.index, params.contains.as_deref()), (Some(holder + 1), Some("geneX")));
    }

    #[test]
    fn a_new_step_follows_the_step_shown_and_takes_its_output() {
        let mut tab = ProtocolTab::new(Protocol::default(), None);
        tab.add_step("open");
        tab.add_step("digest");
        tab.add_step("pick");
        let ids: Vec<&str> = tab.protocol.steps.iter().map(|s| s.id.as_str()).collect();
        assert_eq!(ids, vec!["open", "digest", "pick"]);
        assert_eq!(tab.protocol.steps[2].inputs, vec!["digest".to_string()]);
        // Moving a step keeps it shown, and deleting it leaves the rest.
        tab.move_step(2, true);
        assert_eq!(tab.selected, 1);
        assert_eq!(tab.protocol.steps[1].id, "pick");
        tab.remove_step(1);
        assert_eq!(tab.protocol.steps.len(), 2);
    }

    #[test]
    fn quick_changes_to_one_field_share_an_undo_step() {
        let mut tab = ProtocolTab::new(cloning(&["EcoRI", "BamHI"], None), None);
        let original = tab.protocol.clone();
        for size in [100, 200, 300] {
            tab.change(Some("size".into()), |p| {
                if let Action::Pick(params) = &mut p.steps[2].action {
                    params.size = Some(size);
                }
            });
        }
        assert!(tab.undo());
        assert_eq!(tab.protocol, original, "one undo takes back the whole drag");
    }

    /// Run frames of the tab, each with the panel beside it.
    fn frames(tab: &mut ProtocolTab, count: usize) {
        let ctx = egui::Context::default();
        crate::style::install_fonts(&ctx);
        let settings = ViewSettings::default();
        for _ in 0..count {
            let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
                let requests = show(ui, tab, &[], &settings, true, true);
                assert!(requests.open.is_empty());
            });
            output.textures_delta.clear();
        }
    }

    #[test]
    fn a_click_on_a_card_zooms_in_and_escape_zooms_out() {
        let mut tab = ProtocolTab::new(cloning(&["EcoRI", "BamHI"], Some("geneX")), None);
        tab.refresh(&[], false);
        let ctx = egui::Context::default();
        crate::style::install_fonts(&ctx);
        let settings = ViewSettings::default();
        let screen = Rect::from_min_size(Pos2::ZERO, vec2(1300.0, 850.0));
        let mut time = 0.0;
        // One frame of the tab with these events, at 60 frames a second.
        let mut frame = |tab: &mut ProtocolTab, events: Vec<egui::Event>| {
            time += 1.0 / 60.0;
            let raw = egui::RawInput { events, time: Some(time), screen_rect: Some(screen), ..Default::default() };
            let mut output = ctx.run_ui(raw, |ui| {
                show(ui, tab, &[], &settings, true, true);
            });
            output.textures_delta.clear();
        };
        frame(&mut tab, Vec::new());
        frame(&mut tab, Vec::new());
        // A click in the middle of the donor's cut card.
        let cut = tab.protocol.steps.iter().position(|s| s.id == "cut_donor").unwrap();
        let at = tab.rings.iter().find(|(s, _, _)| *s == cut).map(|(_, _, card)| card.center()).expect("the cut's card is on the bench");
        let button = |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: egui::Modifiers::default() };
        frame(&mut tab, vec![egui::Event::PointerMoved(at)]);
        frame(&mut tab, vec![button(true)]);
        frame(&mut tab, vec![button(false)]);
        assert_eq!((tab.zoom, tab.selected), (Some(cut), cut), "the click zooms into the cut and shows it in the panel");
        // Once the zoom has run, Escape zooms back out to the bench.
        for _ in 0..40 {
            frame(&mut tab, Vec::new());
        }
        let escape = |pressed| egui::Event::Key { key: Key::Escape, physical_key: None, pressed, repeat: false, modifiers: egui::Modifiers::default() };
        frame(&mut tab, vec![escape(true)]);
        frame(&mut tab, vec![escape(false)]);
        assert_eq!(tab.zoom, None);
    }

    /// The cloning cut with EcoRI alone, so the insert fits either way round, and without the
    /// dephosphorylation, so the backbone can also close on itself.
    fn one_enzyme_with_phosphates() -> ProtocolTab {
        let mut tab = ProtocolTab::new(cloning(&["EcoRI"], Some("geneX")), None);
        tab.change(None, |p| p.leave_out("cip").unwrap());
        // The strip beside the zoomed join is read whole, so the bases under the map stay hidden.
        tab.show_sequence = false;
        tab.refresh(&[], false);
        tab
    }

    #[test]
    fn a_join_shows_every_plasmid_its_tube_makes() {
        let mut tab = one_enzyme_with_phosphates();
        let product = step_index(&tab, "product");
        let result = &tab.run.steps[product];
        assert_eq!(result.status, Status::NeedsChoice);
        assert_eq!((result.choices.len(), result.findings.byproducts.len()), (2, 1));
        // The join's card widens for both options and the empty backbone, each with its words.
        let laid = board(&tab.protocol, &tab.run);
        let card = laid.card_of(product).unwrap();
        assert!(laid.cards[card].width > 1.5, "{}", laid.cards[card].width);
        assert!(laid.cards.iter().all(|c| c.column - (c.width - 1.0) / 2.0 >= 0.0), "no card reaches left of the bench");
        let mut screen = Screen::new();
        screen.settle(&mut tab);
        for words in ["geneX piece as given", "geneX piece turned round", "also forms", "no insert", "Choose the direction · 1 more can form"] {
            assert!(screen.texts.iter().any(|t| t == words), "{words} is missing from {:?}", screen.texts);
        }
        // Once a direction is chosen, the other direction stays on the card as a plasmid that
        // also forms, and using it instead turns the insert round.
        tab.choose(product, 0);
        tab.refresh(&[], false);
        let result = &tab.run.steps[product];
        assert_eq!(result.status, Status::Done);
        assert_eq!(result.findings.byproducts.len(), 2);
        let other = result.findings.byproducts.iter().position(|b| b.pin.is_some()).unwrap();
        screen.settle(&mut tab);
        assert!(screen.texts.iter().any(|t| t == "other way round"), "{:?}", screen.texts);
        apply(&mut tab, vec![Edit::UseInstead(product, other)]);
        tab.refresh(&[], false);
        assert_eq!(tab.run.steps[product].findings.flipped, vec![true]);
        assert!(tab.notice.as_ref().is_some_and(|(n, _)| n.starts_with("The geneX piece now goes in the other way round.")), "{:?}", tab.notice);
        // Zoomed in, the other plasmids sit in a strip beside the new plasmid's map.
        tab.zoom_to(product);
        screen.settle(&mut tab);
        assert!(screen.texts.iter().any(|t| t == "The same tube also makes"), "{:?}", screen.texts);
        assert!(screen.texts.iter().any(|t| t == "Use this instead"), "{:?}", screen.texts);
        assert!(screen.texts.iter().any(|t| t == "Remove the 5' phosphates of the backbone"), "{:?}", screen.texts);
    }

    #[test]
    fn every_step_draws_on_the_bench_and_zoomed_in() {
        let tabs = [&["EcoRI", "BamHI"][..], &[][..]].map(|enzymes| ProtocolTab::new(cloning(enzymes, Some("geneX")), None));
        let copies = [ProtocolTab::new(cloning_by_pcr(&["EcoRI", "BamHI"], Some("geneX")), None), ProtocolTab::new(cloning_by_pcr(&[], None), None)];
        for mut tab in tabs.into_iter().chain([one_enzyme_with_phosphates()]).chain(copies) {
            tab.refresh(&[], false);
            // The bench with each step shown in the panel, then each card zoomed into.
            for ii in 0..tab.protocol.steps.len() {
                tab.selected = ii;
                tab.zoom = None;
                frames(&mut tab, 2);
                tab.zoom_to(ii);
                frames(&mut tab, 2);
            }
        }
        // An empty protocol draws its place to start.
        let mut tab = ProtocolTab::new(Protocol::default(), None);
        frames(&mut tab, 2);
    }

    /// The tab drawn frame by frame in a window of fixed size, without the panel, and fed with
    /// pointer events. The protocol runs again after every frame, as the app runs it.
    struct Screen {
        ctx: egui::Context,
        time: f64,
        area: Rect,
        /// Every text the last frame painted, to check what a person would read.
        texts: Vec<String>,
    }

    impl Screen {
        fn new() -> Screen {
            let ctx = egui::Context::default();
            crate::style::install_fonts(&ctx);
            Screen { ctx, time: 0.0, area: Rect::from_min_size(Pos2::ZERO, vec2(1200.0, 900.0)), texts: Vec::new() }
        }

        fn frame(&mut self, tab: &mut ProtocolTab, events: Vec<egui::Event>) {
            self.time += 1.0 / 60.0;
            let raw = egui::RawInput { events, time: Some(self.time), screen_rect: Some(self.area), ..Default::default() };
            let settings = ViewSettings::default();
            let mut output = self.ctx.run_ui(raw, |ui| {
                show(ui, tab, &[], &settings, true, false);
            });
            output.textures_delta.clear();
            self.texts = output
                .shapes
                .iter()
                .filter_map(|clipped| match &clipped.shape {
                    egui::Shape::Text(text) => Some(text.galley.text().to_string()),
                    _ => None,
                })
                .collect();
            tab.refresh(&[], false);
        }

        /// Frames enough for a zoom to finish.
        fn settle(&mut self, tab: &mut ProtocolTab) {
            for _ in 0..40 {
                self.frame(tab, Vec::new());
            }
        }

        fn button(at: Pos2, pressed: bool) -> egui::Event {
            egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: egui::Modifiers::default() }
        }

        fn click(&mut self, tab: &mut ProtocolTab, at: Pos2) {
            self.frame(tab, vec![egui::Event::PointerMoved(at)]);
            self.frame(tab, vec![Screen::button(at, true)]);
            self.frame(tab, vec![Screen::button(at, false)]);
            self.frame(tab, Vec::new());
        }

        /// Press at `from` and move to `to` in steps, holding the button.
        fn hold_and_move(&mut self, tab: &mut ProtocolTab, from: Pos2, to: Pos2) {
            self.frame(tab, vec![egui::Event::PointerMoved(from)]);
            self.frame(tab, vec![Screen::button(from, true)]);
            for ii in 1..=12 {
                self.frame(tab, vec![egui::Event::PointerMoved(from.lerp(to, ii as f32 / 12.0))]);
            }
        }

        /// Let go of the button at `at`.
        fn release(&mut self, tab: &mut ProtocolTab, at: Pos2) {
            self.frame(tab, vec![Screen::button(at, false)]);
            self.frame(tab, Vec::new());
        }

        /// Press at `from`, move to `to` in steps, and let go there.
        fn drag(&mut self, tab: &mut ProtocolTab, from: Pos2, to: Pos2) {
            self.hold_and_move(tab, from, to);
            self.release(tab, to);
        }

        /// The detailed map of plasmid step `step` as the zoomed view draws it now, the rectangle it
        /// is drawn in, and the base at its top.
        fn plasmid_map(&self, tab: &mut ProtocolTab, step: usize) -> (Scene, Rect, usize) {
            let laid = board(&tab.protocol, &tab.run);
            let near = focused_cut(tab, &laid).and_then(|cut| tab.nearest(cut).cloned().map(|n| (cut, n)));
            let card = laid.card_of(step).unwrap();
            let rect = picture_rect(tab, &laid, card, self.area);
            let measure = EguiMeasure { ctx: &self.ctx };
            let map = plasmid_map(tab, &laid, card, near.as_ref(), &ViewSettings::default(), &Palette::get(false), rect.size(), None, &measure).unwrap();
            (map.scene, rect, map.top)
        }
    }

    /// A point of the first shape that draws `tag`: the middle of a line, or a polygon's centre.
    fn point_of(scene: &Scene, tag: &Tag) -> Option<Pos2> {
        use slopgene_core::render::Shape;
        scene.items.iter().find(|item| item.tag == *tag).and_then(|item| match &item.shape {
            Shape::Polyline { points, .. } if points.len() == 2 => Some(pos2((points[0][0] + points[1][0]) / 2.0, (points[0][1] + points[1][1]) / 2.0)),
            Shape::Polyline { points, .. } => points.get(points.len() / 2).map(|p| pos2(p[0], p[1])),
            Shape::Polygon { points, .. } => {
                let n = points.len().max(1) as f32;
                Some(pos2(points.iter().map(|p| p[0]).sum::<f32>() / n, points.iter().map(|p| p[1]).sum::<f32>() / n))
            }
            _ => None,
        })
    }

    fn step_index(tab: &ProtocolTab, id: &str) -> usize {
        tab.protocol.steps.iter().position(|s| s.id == id).unwrap()
    }

    fn part_of(tab: &ProtocolTab, cut: &str) -> Option<String> {
        match &tab.protocol.step(cut)?.action {
            Action::Digest(d) => d.around.clone(),
            _ => None,
        }
    }

    fn enzymes_of(tab: &ProtocolTab, cut: &str) -> Vec<String> {
        match &tab.protocol.step(cut).map(|s| &s.action) {
            Some(Action::Digest(d)) => d.enzymes.clone(),
            _ => Vec::new(),
        }
    }

    #[test]
    fn on_the_zoomed_map_a_click_marks_the_part_and_picks_a_cut_and_a_drag_marks_bases() {
        let mut tab = ProtocolTab::new(cloning(&[], None), None);
        tab.refresh(&[], false);
        let mut screen = Screen::new();
        let donor = step_index(&tab, "donor");
        tab.zoom_to(donor);
        screen.settle(&mut tab);
        // A click on geneX makes it the part to cut out.
        let (scene, rect, _) = screen.plasmid_map(&mut tab, donor);
        let gene = tab.run.steps[donor].outputs[0].doc.features[0].uid;
        let at = point_of(&scene, &Tag::Feature(gene)).expect("geneX is drawn");
        screen.click(&mut tab, rect.min + at.to_vec2());
        assert_eq!(part_of(&tab, "cut_donor").as_deref(), Some("geneX"));
        // The enzymes nearest geneX are marked on the map. A click on EcoRI's cut chooses it, and a
        // second click takes it off.
        screen.settle(&mut tab);
        let (scene, rect, _) = screen.plasmid_map(&mut tab, donor);
        let ecori = scene
            .items
            .iter()
            .find_map(|item| match &item.tag {
                Tag::Enzyme(name, pos) if name == "EcoRI" => Some(Tag::Enzyme(name.clone(), *pos)),
                _ => None,
            })
            .expect("EcoRI, the nearest before geneX, is marked");
        let at = rect.min + point_of(&scene, &ecori).unwrap().to_vec2();
        screen.click(&mut tab, at);
        assert_eq!(enzymes_of(&tab, "cut_donor"), vec!["EcoRI".to_string()]);
        screen.settle(&mut tab);
        screen.click(&mut tab, at);
        assert!(enzymes_of(&tab, "cut_donor").is_empty(), "{:?}", enzymes_of(&tab, "cut_donor"));
        // A drag across the map marks the bases it passes as the part.
        screen.settle(&mut tab);
        let (_, rect, top) = screen.plasmid_map(&mut tab, donor);
        let len = tab.run.steps[donor].outputs[0].doc.len();
        let radius = (rect.height() * 0.31).min(rect.width() * 0.23) * 0.7;
        let point = |base: usize| {
            let angle = (base as f32 - top as f32) / len as f32 * std::f32::consts::TAU;
            rect.center() + vec2(radius * angle.sin(), -radius * angle.cos())
        };
        screen.drag(&mut tab, point(300), point(700));
        let part = part_of(&tab, "cut_donor").expect("the drag marks a part");
        let span = span::parse_range(&part, len, true).unwrap();
        assert!(span.start.abs_diff(300) <= 2 && span.end.abs_diff(700) <= 2, "{part}");
    }

    #[test]
    fn a_click_on_a_cut_that_waits_for_its_enzymes_opens_the_map_where_they_are_chosen() {
        let mut tab = ProtocolTab::new(cloning(&[], None), None);
        tab.refresh(&[], false);
        let mut screen = Screen::new();
        screen.settle(&mut tab);
        let (donor, cut) = (step_index(&tab, "donor"), step_index(&tab, "cut_donor"));
        // The cut has no pieces yet, so a click on its card zooms into the donor's map. The map is
        // headed by the cut, and its hint asks for the part first.
        let at = tab.rings.iter().find(|(s, _, _)| *s == cut).map(|(_, _, card)| card.center()).unwrap();
        screen.click(&mut tab, at);
        assert_eq!(tab.zoom, Some(donor));
        screen.settle(&mut tab);
        assert!(screen.texts.iter().any(|t| t.starts_with("Cut donor")), "{:?}", screen.texts);
        assert!(screen.texts.iter().any(|t| t.starts_with("Drag across the map to mark the part you want")), "{:?}", screen.texts);
        // A drag across the map marks the part.
        let (_, rect, top) = screen.plasmid_map(&mut tab, donor);
        let len = tab.run.steps[donor].outputs[0].doc.len();
        let radius = (rect.height() * 0.31).min(rect.width() * 0.23) * 0.7;
        let point = |base: usize| {
            let angle = (base as f32 - top as f32) / len as f32 * std::f32::consts::TAU;
            rect.center() + vec2(radius * angle.sin(), -radius * angle.cos())
        };
        screen.drag(&mut tab, point(300), point(700));
        assert!(part_of(&tab, "cut_donor").is_some());
        // The cuts nearest the part then show on the map. A click on EcoRI before it and on BamHI
        // after it chooses both, and the view stays on the map all the while.
        for (name, chosen) in [("EcoRI", vec!["EcoRI"]), ("BamHI", vec!["EcoRI", "BamHI"])] {
            screen.settle(&mut tab);
            let (scene, rect, _) = screen.plasmid_map(&mut tab, donor);
            let tag = scene
                .items
                .iter()
                .find_map(|item| match &item.tag {
                    Tag::Enzyme(n, pos) if n == name => Some(Tag::Enzyme(n.clone(), *pos)),
                    _ => None,
                })
                .unwrap_or_else(|| panic!("{name} is marked"));
            screen.click(&mut tab, rect.min + point_of(&scene, &tag).unwrap().to_vec2());
            assert_eq!(enzymes_of(&tab, "cut_donor"), chosen);
            assert_eq!(tab.zoom, Some(donor));
        }
        assert_eq!(tab.run.steps[cut].status, Status::Done, "{:?}", tab.run.steps[cut].error);
    }

    #[test]
    fn a_click_on_a_piece_in_the_zoomed_cut_keeps_it() {
        let mut tab = ProtocolTab::new(cloning(&["EcoRI", "BamHI"], Some("geneX")), None);
        tab.refresh(&[], false);
        let mut screen = Screen::new();
        let cut = step_index(&tab, "cut_donor");
        tab.zoom_to(cut);
        screen.settle(&mut tab);
        // The ring as the zoomed view draws it, and the piece that does not hold geneX.
        let laid = board(&tab.protocol, &tab.run);
        let card = laid.card_of(cut).unwrap();
        let shown = picture(&tab, &laid, card, true, &ViewSettings::default(), Palette::get(false).canvas).unwrap();
        let rect = picture_rect(&tab, &laid, card, screen.area);
        let mut options = shown.options;
        options.width = rect.width();
        options.height = rect.height();
        options.font_size = 12.5;
        let scene = ring::render_ring(&shown.doc, &options, &EguiMeasure { ctx: &screen.ctx });
        let other = 1 - tab.run.steps[cut].findings.holder.unwrap();
        let view = cut_view(&tab, &laid.cards[card]).unwrap();
        let piece = view.pieces.iter().position(|p| p.2 == Some(other)).unwrap();
        let at = point_of(&scene, &Tag::Piece(piece)).unwrap();
        screen.click(&mut tab, rect.min + at.to_vec2());
        let Action::Pick(params) = &tab.protocol.step("insert").unwrap().action else { panic!("the insert is picked") };
        assert_eq!((params.index, params.contains.as_deref()), (Some(other + 1), None));
    }

    #[test]
    fn dragging_a_kept_piece_puts_it_into_a_plasmid_or_the_new_plasmid() {
        let (vector, donor) = maps();
        let stored = |doc: &Document| OpenParams { path: None, document: Some(json::document_to_value(doc, true)) };
        let names = |list: &[&str]| list.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        // The insert's steps, and the vector opened beside them with nothing done to it yet.
        let mut protocol = Protocol::default();
        protocol.add_source(stored(&donor), names(&["EcoRI", "BamHI"]), Some("geneX".into()));
        protocol.steps.push(Step { id: "vector".into(), inputs: Vec::new(), action: Action::Open(stored(&vector)), note: String::new() });
        let mut tab = ProtocolTab::new(protocol, None);
        tab.refresh(&[], false);
        let mut screen = Screen::new();
        screen.settle(&mut tab);
        let card_at = |tab: &ProtocolTab, id: &str| {
            let step = step_index(tab, id);
            tab.rings.iter().find(|(s, _, _)| *s == step).map(|(_, _, card)| card.center()).unwrap()
        };
        // The donor's kept piece, dragged onto the vector: while it is held over the vector, the
        // vector says what letting go does, and the piece's name rides with the pointer.
        let (from, to) = (card_at(&tab, "cut_donor"), card_at(&tab, "vector"));
        screen.hold_and_move(&mut tab, from, to);
        assert!(screen.texts.iter().any(|t| t == "Let go to put the geneX piece into vector"), "{:?}", screen.texts);
        assert!(screen.texts.iter().any(|t| t.starts_with("geneX · ")), "{:?}", screen.texts);
        // Letting go there cuts the vector the same way, keeps its backbone without its
        // phosphates, and joins the two.
        screen.release(&mut tab, to);
        let join = tab.protocol.steps.iter().find(|s| matches!(s.action, Action::Ligate(_))).expect("a join").clone();
        assert_eq!(join.inputs, names(&["cip", "insert"]));
        screen.settle(&mut tab);
        assert_eq!(tab.run.step(&join.id).unwrap().status, Status::Done, "{:?}", tab.run.step(&join.id).unwrap().error);
        // A second donor's kept piece, dragged onto the new plasmid, joins it as well.
        tab.change(None, |p| {
            p.add_source(stored(&donor), names(&["EcoRI", "BamHI"]), Some("geneX".into()));
        });
        tab.refresh(&[], false);
        screen.settle(&mut tab);
        let (from, to) = (card_at(&tab, "cut_donor_2"), card_at(&tab, &join.id));
        screen.drag(&mut tab, from, to);
        assert_eq!(tab.protocol.step(&join.id).unwrap().inputs, names(&["cip", "insert", "insert_2"]));
        // The way back is one undo.
        assert!(tab.undo());
        assert_eq!(tab.protocol.step(&join.id).unwrap().inputs, names(&["cip", "insert"]));
    }

    #[test]
    fn a_pcr_lane_has_a_copy_card_and_a_lane_switches_between_cutting_and_copying() {
        let mut tab = ProtocolTab::new(cloning_by_pcr(&["EcoRI", "BamHI"], Some("geneX")), None);
        tab.refresh(&[], false);
        for step in &tab.run.steps {
            assert_eq!(step.status, Status::Done, "{}: {:?}", step.id, step.error);
        }
        // The PCR card sits between the donor and its cut, and the join one row lower than in
        // a cut lane.
        let laid = board(&tab.protocol, &tab.run);
        let kinds: Vec<(Kind, f32, usize)> = laid.cards.iter().map(|c| (c.kind, c.column, c.row)).collect();
        let expected = vec![
            (Kind::Plasmid, 0.0, 0),
            (Kind::Cut, 0.0, 1),
            (Kind::Plasmid, 1.0, 0),
            (Kind::Copy, 1.0, 1),
            (Kind::Cut, 1.0, 2),
            (Kind::Join, 0.5, 3),
            (Kind::AddPlasmid, 2.0, 0),
        ];
        assert_eq!(kinds, expected);
        let copy = step_index(&tab, "pcr_donor");
        let copy_card = laid.card_of(copy).unwrap();
        assert_eq!(
            (laid.copy_after(2), laid.cut_after(copy_card), laid.before(4), laid.part_step_after(2)),
            (Some(copy_card), Some(4), Some(copy_card), Some(copy))
        );
        // The PCR card draws the product as a bar: the copy green between grey tails. The donor's
        // card shows the part copied, with the rest left as it is.
        let settings = ViewSettings::default();
        let palette = Palette::get(false);
        let shown = picture(&tab, &laid, copy_card, false, &settings, palette.raised).unwrap();
        assert!(!shown.doc.is_circular());
        let pieces: Vec<(Span, Role)> = shown.options.pieces.iter().map(|p| (p.span, p.role)).collect();
        assert_eq!(pieces, vec![(Span::new(0, 12), Role::Plain), (Span::new(12, 612), Role::Kept), (Span::new(612, 624), Role::Plain)]);
        let plasmid = picture(&tab, &laid, 2, false, &settings, palette.raised).unwrap();
        assert_eq!(plasmid.options.pieces.first().map(|p| (p.span, p.role)), Some((Span::new(200, 800), Role::Kept)));
        assert!(plasmid.options.pieces.iter().skip(1).all(|p| p.role == Role::Plain));
        let (verb, object, line, _) = card_text(&tab, &laid, copy_card, &palette);
        assert_eq!((verb.as_str(), object.as_str(), line.as_str()), ("COPY", "geneX by PCR", "624 bp product · EcoRI + BamHI tails"));
        // The whole bench, with the words a person reads on it.
        let mut screen = Screen::new();
        screen.settle(&mut tab);
        for words in ["COPY", "geneX by PCR", "PCR with EcoRI + BamHI tails", "trim the tails with EcoRI + BamHI", "PCR product"] {
            assert!(screen.texts.iter().any(|t| t == words), "{words} is missing from {:?}", screen.texts);
        }
        // Back to cutting the part out, and to copying it again, from the bench.
        apply(&mut tab, vec![Edit::CutOut(copy)]);
        assert!(tab.protocol.step("pcr_donor").is_none());
        assert!(tab.notice.as_ref().is_some_and(|(n, _)| n.starts_with("The part is now cut out of donor again")), "{:?}", tab.notice);
        tab.refresh(&[], false);
        assert_eq!(tab.run.step("product").unwrap().status, Status::Done);
        let donor = step_index(&tab, "donor");
        apply(&mut tab, vec![Edit::CopyOut(donor)]);
        assert!(tab.protocol.step("pcr_donor").is_some());
        assert_eq!(tab.selected, step_index(&tab, "pcr_donor"), "the PCR step is shown, where the tails are chosen");
        tab.refresh(&[], false);
        assert_eq!(tab.run.step("product").unwrap().status, Status::Done, "{:?}", tab.run.step("product").unwrap().error);
        // The way back is one undo each.
        assert!(tab.undo());
        assert!(tab.protocol.step("pcr_donor").is_none());
    }

    #[test]
    fn tail_enzymes_are_picked_one_per_primer_and_the_primers_can_be_kept() {
        let mut tab = ProtocolTab::new(cloning_by_pcr(&[], Some("geneX")), None);
        tab.refresh(&[], false);
        let copy = step_index(&tab, "pcr_donor");
        assert_eq!(tab.run.steps[copy].status, Status::Done, "a PCR without tails still makes its product");
        // The chooser lists the enzymes of the set that do not cut inside geneX, those that open
        // the vector once first.
        let listed = tab.tails(copy).unwrap().to_vec();
        assert!(listed.iter().take(2).all(|c| c.vector_cuts == Some(1)), "{:?}", listed.iter().take(3).collect::<Vec<_>>());
        assert!(listed.iter().any(|c| c.name == "EcoRI" && !c.cuts_part && c.vector_cuts == Some(1)), "{listed:?}");
        // One click per primer sets its tail. The cut that trims the tails and the vector follow.
        apply(&mut tab, vec![Edit::PickTail(copy, true, "EcoRI".into())]);
        tab.refresh(&[], false);
        apply(&mut tab, vec![Edit::PickTail(copy, false, "BamHI".into())]);
        tab.refresh(&[], false);
        assert_eq!(enzymes_of(&tab, "cut_pcr"), vec!["EcoRI".to_string(), "BamHI".to_string()]);
        assert_eq!(enzymes_of(&tab, "cut_vector"), vec!["EcoRI".to_string(), "BamHI".to_string()]);
        for step in &tab.run.steps {
            assert_eq!(step.status, Status::Done, "{}: {:?}", step.id, step.error);
        }
        assert!(tab.notice.as_ref().is_some_and(|(n, _)| n.starts_with("The primers for donor carry EcoRI and BamHI sites")), "{:?}", tab.notice);
        // A second click takes a tail off.
        apply(&mut tab, vec![Edit::PickTail(copy, false, "BamHI".into())]);
        let Action::Amplify(params) = &tab.protocol.steps[copy].action else { panic!("pcr_donor is a PCR") };
        assert_eq!((params.forward_enzyme.as_deref(), params.reverse_enzyme.as_deref()), (Some("EcoRI"), None));
        apply(&mut tab, vec![Edit::PickTail(copy, false, "BamHI".into())]);
        tab.refresh(&[], false);
        // Keeping the primers writes them into the step, and designing them anew lets them go.
        let designed = tab.run.steps[copy].findings.primers[0].seq.clone();
        apply(&mut tab, vec![Edit::Primers(copy, PrimerEdit::Pair(0))]);
        let Action::Amplify(params) = &tab.protocol.steps[copy].action else { panic!("pcr_donor is a PCR") };
        assert_eq!(params.forward_primer.as_deref(), Some(designed.as_str()));
        apply(&mut tab, vec![Edit::Primers(copy, PrimerEdit::Clear)]);
        let Action::Amplify(params) = &tab.protocol.steps[copy].action else { panic!("pcr_donor is a PCR") };
        assert!(!params.pinned());
        // Zoomed into the PCR, the view is headed by the PCR and shows the donor's map with the
        // primers drawn on it.
        let mut screen = Screen::new();
        tab.selected = copy;
        screen.settle(&mut tab);
        tab.zoom_to(copy);
        screen.settle(&mut tab);
        assert!(screen.texts.iter().any(|t| t == "Copy geneX out of donor by PCR"), "{:?}", screen.texts);
        let (scene, _, _) = screen.plasmid_map(&mut tab, copy);
        assert!(scene.items.iter().any(|i| matches!(i.tag, Tag::Primer(_))), "the primers are drawn on the map");
        assert!(scene.items.iter().any(|i| i.tag == Tag::Selection), "the part copied is shaded");
    }

    #[test]
    fn the_zoomed_map_has_the_bases_under_it_and_the_two_select_together() {
        let mut tab = ProtocolTab::new(cloning(&[], Some("geneX")), None);
        tab.refresh(&[], false);
        let mut screen = Screen::new();
        let donor = step_index(&tab, "donor");
        tab.zoom_to(donor);
        screen.settle(&mut tab);
        // The pane shows the donor's bases, read-only, with geneX selected as on the map.
        let preview = tab.preview.as_ref().expect("the pane has a document");
        assert_eq!((preview.step, preview.tab.selection, preview.tab.read_only), (donor, Some(Span::new(200, 800)), true));
        assert!(screen.texts.iter().any(|t| t == "donor · 1,500 bp"), "{:?}", screen.texts);
        assert!(screen.texts.iter().any(|t| t == "Sequence"), "{:?}", screen.texts);
        // The map takes the room above the pane.
        let laid = board(&tab.protocol, &tab.run);
        let card = laid.card_of(donor).unwrap();
        let (map, pane) = detail_rects(screen.area, sequence_share(&tab, &laid, card));
        assert!(pane.is_some_and(|p| p.min.y > map.max.y && p.max.y == screen.area.max.y));
        // A stretch selected in the pane becomes the part once the button is up, and a feature
        // chosen by its bar goes by its name.
        tab.preview.as_mut().unwrap().tab.select(Some(Span::new(300, 700)));
        screen.frame(&mut tab, Vec::new());
        screen.frame(&mut tab, Vec::new());
        assert_eq!(part_of(&tab, "cut_donor").as_deref(), Some("301..700"));
        assert_eq!(tab.preview.as_ref().unwrap().tab.selection, Some(Span::new(300, 700)), "the selection stays while the part catches up");
        let gene = tab.run.steps[donor].outputs[0].doc.features[0].uid;
        tab.preview.as_mut().unwrap().tab.select_feature(gene);
        screen.frame(&mut tab, Vec::new());
        screen.frame(&mut tab, Vec::new());
        assert_eq!(part_of(&tab, "cut_donor").as_deref(), Some("geneX"));
        // A click that cleared the selection brings the part's highlight back.
        tab.preview.as_mut().unwrap().tab.select(None);
        screen.frame(&mut tab, Vec::new());
        assert_eq!(tab.preview.as_ref().unwrap().tab.selection, Some(Span::new(200, 800)));
        // Hidden, the pane leaves the whole area to the map, and the button brings it back.
        apply(&mut tab, vec![Edit::ShowSequence(false)]);
        screen.settle(&mut tab);
        assert!(!screen.texts.iter().any(|t| t == "donor · 1,500 bp"), "{:?}", screen.texts);
        assert_eq!(picture_rect(&tab, &laid, card, screen.area), detail_picture_rect(screen.area));
        apply(&mut tab, vec![Edit::ShowSequence(true), Edit::SequenceShare(0.9)]);
        assert_eq!(tab.sequence_share, 0.75, "the pane leaves room for the map");
        // A join's pane shows the new plasmid with the insert highlighted.
        tab.change(None, |p| p.set_enzymes("cut_donor", vec!["EcoRI".into(), "BamHI".into()], true).unwrap());
        tab.refresh(&[], false);
        let product = step_index(&tab, "product");
        tab.zoom_to(product);
        screen.settle(&mut tab);
        let preview = tab.preview.as_ref().unwrap();
        let made = &tab.run.steps[product];
        assert_eq!((preview.step, preview.tab.selection), (product, made.findings.parts.last().copied()));
    }
}
