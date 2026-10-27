//! The application: open documents in tabs, the main views, side panels, the command palette,
//! dialogs and the script console.

use crate::align_view;
use crate::doc::{DocTab, Saved, View, format_label};
use crate::protocol_view::{self, ProtocolTab};
use crate::record;
use crate::style::{self, Palette};
use crate::update_view::Updater;
use crate::{map_view, panels, scene_paint, seq_view};
use eframe::egui::{self, Key, KeyboardShortcut, Modifiers, RichText, Sense, vec2};
use slopgene_core::digest::digest;
use slopgene_core::enzymes::EnzymeDb;
use slopgene_core::gel::{GelOptions, Ladder, Lane, render_gel};
use slopgene_core::io::Format;
use slopgene_core::protocol::{CloningTemplate, OpenParams, Protocol, stored_path};
use slopgene_core::render::{ColorMode, thousands};
use slopgene_core::search::{Hit, find_dna, find_protein};
use slopgene_core::translate::GeneticCode;
use slopgene_core::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// Display options shared by every view.
#[derive(Clone, Debug)]
pub struct ViewSettings {
    pub dark: bool,
    pub colors: ColorMode,
    pub show_complement: bool,
    pub show_translations: bool,
    pub show_enzymes: bool,
    pub show_primers: bool,
}

impl Default for ViewSettings {
    fn default() -> Self {
        ViewSettings { dark: false, colors: ColorMode::Palette, show_complement: true, show_translations: true, show_enzymes: true, show_primers: true }
    }
}

/// Everything a user can ask for, from the command palette, shortcuts or buttons.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Action {
    Open,
    NewSequence,
    NewProtocol,
    Save,
    SaveAs,
    Close,
    Reload,
    Undo,
    Redo,
    CopySelection,
    CopyReverseComplement,
    Cut,
    SelectAll,
    ViewMap,
    ViewSequence,
    ViewSplit,
    ToggleDark,
    ToggleColors,
    ToggleComplement,
    ToggleTranslations,
    ToggleEnzymes,
    TogglePrimers,
    ToggleInspector,
    ToggleConsole,
    Find,
    GoTo,
    ReverseComplement,
    SetOriginAtCursor,
    ToggleTopology,
    AddFeatureFromSelection,
    DeleteSelectedFeature,
    PrimerFromSelection,
    Digest,
    ExportSvg,
    ExportPng,
    ExportBench,
    ExportGenBank,
    ExportFasta,
    AlignRead,
    DetectFeatures,
    Palette,
    Update,
}

impl Action {
    /// The palette's commands in the order it lists them, group by group (see [`Action::group`]).
    pub const ALL: &'static [Action] = &[
        Action::Open,
        Action::NewSequence,
        Action::NewProtocol,
        Action::Save,
        Action::SaveAs,
        Action::Close,
        Action::Reload,
        Action::Update,
        Action::AlignRead,
        Action::Digest,
        Action::DetectFeatures,
        Action::Undo,
        Action::Redo,
        Action::Find,
        Action::GoTo,
        Action::SelectAll,
        Action::CopySelection,
        Action::CopyReverseComplement,
        Action::Cut,
        Action::AddFeatureFromSelection,
        Action::DeleteSelectedFeature,
        Action::PrimerFromSelection,
        Action::ReverseComplement,
        Action::SetOriginAtCursor,
        Action::ToggleTopology,
        Action::ViewMap,
        Action::ViewSequence,
        Action::ViewSplit,
        Action::ToggleInspector,
        Action::ToggleConsole,
        Action::ToggleDark,
        Action::ToggleColors,
        Action::ToggleComplement,
        Action::ToggleTranslations,
        Action::ToggleEnzymes,
        Action::TogglePrimers,
        Action::ExportSvg,
        Action::ExportPng,
        Action::ExportBench,
        Action::ExportGenBank,
        Action::ExportFasta,
    ];

    /// The heading the command sits under in the palette. The tools come straight after the file
    /// commands, so they show without scrolling.
    pub fn group(self) -> &'static str {
        match self {
            Action::Open | Action::NewSequence | Action::NewProtocol | Action::Save | Action::SaveAs | Action::Close | Action::Reload | Action::Update => {
                "File"
            }
            Action::AlignRead | Action::Digest | Action::DetectFeatures => "Tools",
            Action::ViewMap
            | Action::ViewSequence
            | Action::ViewSplit
            | Action::ToggleInspector
            | Action::ToggleConsole
            | Action::ToggleDark
            | Action::ToggleColors
            | Action::ToggleComplement
            | Action::ToggleTranslations
            | Action::ToggleEnzymes
            | Action::TogglePrimers
            | Action::Palette => "View",
            Action::ExportSvg | Action::ExportPng | Action::ExportBench | Action::ExportGenBank | Action::ExportFasta => "Export",
            _ => "Edit",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Action::Open => "Open file…",
            Action::NewSequence => "New sequence…",
            Action::NewProtocol => "New cloning protocol",
            Action::Save => "Save",
            Action::SaveAs => "Save as…",
            Action::Close => "Close document",
            Action::Reload => "Reload from disk",
            Action::Undo => "Undo",
            Action::Redo => "Redo",
            Action::CopySelection => "Copy selection",
            Action::CopyReverseComplement => "Copy reverse complement",
            Action::Cut => "Cut selection",
            Action::SelectAll => "Select all",
            Action::ViewMap => "Show map",
            Action::ViewSequence => "Show sequence",
            Action::ViewSplit => "Show map and sequence",
            Action::ToggleDark => "Toggle dark mode",
            Action::ToggleColors => "Toggle file colours / palette colours",
            Action::ToggleComplement => "Toggle complementary strand",
            Action::ToggleTranslations => "Toggle translations",
            Action::ToggleEnzymes => "Toggle enzyme sites",
            Action::TogglePrimers => "Toggle primers",
            Action::ToggleInspector => "Toggle inspector panel",
            Action::ToggleConsole => "Toggle script console",
            Action::Find => "Find…",
            Action::GoTo => "Go to position…",
            Action::ReverseComplement => "Reverse-complement the whole sequence",
            Action::SetOriginAtCursor => "Set origin at cursor",
            Action::ToggleTopology => "Toggle circular / linear",
            Action::AddFeatureFromSelection => "Add feature from selection…",
            Action::DeleteSelectedFeature => "Delete selected feature",
            Action::PrimerFromSelection => "Save selection as primer",
            Action::Digest => "Restriction digest…",
            Action::ExportSvg => "Export map as SVG…",
            Action::ExportPng => "Export map as PNG…",
            Action::ExportBench => "Export the bench as a picture…",
            Action::ExportGenBank => "Export as GenBank…",
            Action::ExportFasta => "Export as FASTA…",
            Action::AlignRead => "Align reads or other sequences to this one…",
            Action::DetectFeatures => "Detect features from a library folder…",
            Action::Palette => "Command palette",
            Action::Update => "Check for updates…",
        }
    }

    pub fn shortcut(self) -> Option<KeyboardShortcut> {
        let cmd = Modifiers::COMMAND;
        let cmd_shift = Modifiers::COMMAND.plus(Modifiers::SHIFT);
        Some(match self {
            Action::Open => KeyboardShortcut::new(cmd, Key::O),
            Action::NewSequence => KeyboardShortcut::new(cmd, Key::N),
            Action::NewProtocol => KeyboardShortcut::new(cmd_shift, Key::N),
            Action::Save => KeyboardShortcut::new(cmd, Key::S),
            Action::SaveAs => KeyboardShortcut::new(cmd_shift, Key::S),
            Action::Close => KeyboardShortcut::new(cmd, Key::W),
            Action::Undo => KeyboardShortcut::new(cmd, Key::Z),
            Action::Redo => KeyboardShortcut::new(cmd_shift, Key::Z),
            // Copy and cut arrive as clipboard events (see `read_shortcuts`); these entries
            // show the keys in the command palette.
            Action::CopySelection => KeyboardShortcut::new(cmd, Key::C),
            Action::Cut => KeyboardShortcut::new(cmd, Key::X),
            Action::SelectAll => KeyboardShortcut::new(cmd, Key::A),
            Action::Find => KeyboardShortcut::new(cmd, Key::F),
            Action::GoTo => KeyboardShortcut::new(cmd, Key::G),
            Action::ViewMap => KeyboardShortcut::new(cmd, Key::Num1),
            Action::ViewSequence => KeyboardShortcut::new(cmd, Key::Num2),
            Action::ViewSplit => KeyboardShortcut::new(cmd, Key::Num3),
            Action::ToggleConsole => KeyboardShortcut::new(cmd, Key::J),
            Action::ToggleInspector => KeyboardShortcut::new(cmd, Key::I),
            Action::Digest => KeyboardShortcut::new(cmd_shift, Key::D),
            Action::AlignRead => KeyboardShortcut::new(cmd_shift, Key::A),
            Action::Palette => KeyboardShortcut::new(cmd, Key::K),
            _ => return None,
        })
    }

    /// Editing actions whose keys belong to a text field while one has the keyboard (Ctrl+Z
    /// then undoes typing in that field), and that stay off while a dialog is open.
    fn is_editing(self) -> bool {
        matches!(self, Action::Undo | Action::Redo | Action::CopySelection | Action::Cut | Action::SelectAll)
    }
}

/// Every keyboard shortcut, most specific first. egui ignores an extra Shift or Alt when it
/// matches a shortcut, so Ctrl+Shift+Z has to be tried before Ctrl+Z (or redo would undo) and
/// Ctrl+Shift+S before Ctrl+S. Ctrl+Y is a second key for redo.
pub fn shortcut_table() -> Vec<(Action, KeyboardShortcut)> {
    let mut table: Vec<(Action, KeyboardShortcut)> = Action::ALL.iter().chain([Action::Palette].iter()).filter_map(|&a| a.shortcut().map(|s| (a, s))).collect();
    table.push((Action::Redo, KeyboardShortcut::new(Modifiers::COMMAND, Key::Y)));
    table.sort_by_key(|(_, s)| std::cmp::Reverse(s.modifiers.shift as u8 + s.modifiers.alt as u8));
    table
}

/// Whether the editing keys act on the document this frame. They do unless a text field has the
/// keyboard (its own undo and clipboard come first) or a dialog is open. The sequence view and
/// the map take the keyboard too, but they are not text fields, so the keys still reach the
/// document while either has it.
pub fn document_has_keys(ctx: &egui::Context, modal_open: bool) -> bool {
    !ctx.text_edit_focused() && !modal_open
}

/// What the keyboard asked for in one frame.
#[derive(Debug, Default, PartialEq)]
pub struct KeyRequests {
    pub actions: Vec<Action>,
    /// Text pasted with Ctrl+V (or Shift+Insert), in order.
    pub pastes: Vec<String>,
}

/// Read and consume this frame's shortcuts. `document_keys` is false while a text field has the
/// keyboard or a dialog is open: the editing keys (undo, redo, copy, cut, paste, select all) are
/// then left to the text field or kept from the document behind the dialog. Otherwise they act on
/// the open document, whether the sequence view, the map or nothing has the keyboard.
pub fn read_shortcuts(ctx: &egui::Context, document_keys: bool) -> KeyRequests {
    let mut requests = KeyRequests::default();
    for (action, shortcut) in shortcut_table() {
        if action.is_editing() && !document_keys {
            continue;
        }
        if ctx.input_mut(|i| i.consume_shortcut(&shortcut)) {
            requests.actions.push(action);
        }
    }
    // The window system turns Ctrl+C, Ctrl+X and Ctrl+V into clipboard events, not key presses.
    if document_keys {
        let clipboard: Vec<egui::Event> = ctx.input_mut(|i| {
            let mut taken = Vec::new();
            i.events.retain(|event| {
                let ours = matches!(event, egui::Event::Copy | egui::Event::Cut | egui::Event::Paste(_));
                if ours {
                    taken.push(event.clone());
                }
                !ours
            });
            taken
        });
        for event in clipboard {
            match event {
                egui::Event::Copy => requests.actions.push(Action::CopySelection),
                egui::Event::Cut => requests.actions.push(Action::Cut),
                egui::Event::Paste(text) => requests.pastes.push(text),
                _ => {}
            }
        }
    }
    requests
}

/// Options from the command line.
#[derive(Clone, Debug, Default)]
pub struct StartOptions {
    pub files: Vec<PathBuf>,
    pub dark: bool,
    pub view: Option<View>,
    pub select: Option<String>,
    pub feature: Option<String>,
    pub palette: bool,
    pub console: bool,
    pub digest: Option<String>,
    /// Align this read to the first document at start.
    pub align: Option<PathBuf>,
    /// Show this step (from 1) of the protocol that opens.
    pub step: Option<usize>,
    /// Zoom into the card of that step.
    pub zoom: bool,
    /// A drag to play for a screenshot: press at the first point, move to the second, and hold.
    pub drag: Option<(egui::Pos2, egui::Pos2)>,
    /// Start a cloning protocol from the maps open.
    pub new_protocol: bool,
    /// Enzymes that cut the insert out of the new protocol's donor.
    pub enzymes: Option<String>,
    /// Copy the new protocol's insert out by PCR, with the enzymes' sites on the primers' tails.
    pub pcr: bool,
    /// Save a screenshot of the window to this file after `screenshot_after` frames, then quit.
    pub screenshot: Option<PathBuf>,
    pub screenshot_after: u64,
    /// A script of pointer and key actions to play into the window (see `record::Play`).
    pub play: Option<record::Play>,
    /// Record the window into this animated GIF while the script plays, then quit.
    pub record: Option<PathBuf>,
    /// The window size asked for. A picture or a recording sets it on its second frame, since
    /// eframe first restores the size the window had last.
    pub size: [f32; 2],
}

impl StartOptions {
    /// Whether the run makes a picture or a recording, rather than a session with a person.
    pub fn headless(&self) -> bool {
        self.screenshot.is_some() || self.record.is_some()
    }
}

#[derive(Default)]
struct PaletteState {
    open: bool,
    query: String,
    selected: usize,
    /// Bring the selected command into view on the next frame. Opening the palette, the arrow
    /// keys and a new query set it. Scrolling to the selection on every frame would pull the list
    /// back each time the mouse wheel moved the selection out of view.
    reveal: bool,
    /// How far the list is scrolled, in points, as it was last drawn.
    scrolled: f32,
}

impl PaletteState {
    /// A freshly opened palette: no query, the first command selected and the list at its top.
    fn opened() -> PaletteState {
        PaletteState { open: true, reveal: true, ..PaletteState::default() }
    }
}

#[derive(Default)]
struct FindState {
    open: bool,
    query: String,
    mismatches: usize,
    protein: bool,
    hits: Vec<Hit>,
    current: usize,
    /// (tab id, revision, query, mismatches, protein) the hits were computed for.
    key: Option<(u64, u64, String, usize, bool)>,
}

struct NewSequenceDialog {
    name: String,
    sequence: String,
    circular: bool,
    error: Option<String>,
}

struct FeatureDialog {
    name: String,
    kind: String,
    strand: Strand,
}

/// The fragments of a digest (length, description, fragment), or why it failed.
type DigestOutcome = std::result::Result<Vec<(usize, String, Document)>, String>;

struct DigestDialog {
    enzymes: String,
    result: Option<DigestOutcome>,
    lanes: Vec<Lane>,
}

struct Console {
    code: String,
    output: String,
}

/// Which kind of tab is in front.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Front {
    /// The document tab at `active`.
    Document,
    /// The protocol tab with this index.
    Protocol(usize),
}

/// What waits for "save changes?" before it closes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CloseTarget {
    App,
    Document(usize),
    Protocol(usize),
}

pub struct SlopGeneApp {
    tabs: Vec<DocTab>,
    active: usize,
    /// Open protocols, each in a tab after the documents.
    protocols: Vec<ProtocolTab>,
    front: Front,
    settings: ViewSettings,
    show_inspector: bool,
    show_console: bool,
    palette: PaletteState,
    find: FindState,
    goto: Option<String>,
    new_sequence: Option<NewSequenceDialog>,
    feature_dialog: Option<FeatureDialog>,
    digest_dialog: Option<DigestDialog>,
    console: Console,
    /// What waits for "save changes?" before it closes.
    confirm_close: Option<CloseTarget>,
    allow_quit: bool,
    recent: Vec<PathBuf>,
    notice: Option<(String, Instant)>,
    applied_dark: Option<bool>,
    frames: u64,
    start: StartOptions,
    screenshot_done: Arc<AtomicBool>,
    screenshot_requested: bool,
    /// The script being played into the window, and the recording being made of it.
    play: Option<record::Play>,
    recorder: Option<record::Recorder>,
    /// The check for a newer version, and its dialog.
    updater: Updater,
}

impl SlopGeneApp {
    pub fn new(cc: &eframe::CreationContext<'_>, start: StartOptions) -> SlopGeneApp {
        // A picture or a recording starts from the default layout rather than the panel sizes
        // saved by the last session, so the same command always gives the same picture. This must
        // come before the fonts: egui keeps new font definitions in its memory until the next frame.
        if start.headless() {
            cc.egui_ctx.memory_mut(|memory| *memory = egui::Memory::default());
        }
        style::install_fonts(&cc.egui_ctx);
        let mut recent = Vec::new();
        let mut dark = start.dark;
        let mut check_updates = true;
        if let Some(storage) = cc.storage {
            if let Some(list) = storage.get_string("recent") {
                recent = list.lines().filter(|l| !l.is_empty()).map(PathBuf::from).collect();
            }
            if !start.dark && !start.headless() {
                dark = storage.get_string("dark").as_deref() == Some("1");
            }
            check_updates = storage.get_string("check_updates").as_deref() != Some("0");
        }
        // The programs an update renamed aside are removed, and a newer version is looked for.
        // Pictures and recordings stay off the network.
        if let Ok(exe) = std::env::current_exe()
            && let Some(folder) = exe.parent()
        {
            slopgene_update::tidy(folder);
        }
        let mut updater = Updater::new(check_updates);
        if check_updates && !start.headless() {
            updater.check();
        }
        let mut app = SlopGeneApp {
            tabs: Vec::new(),
            active: 0,
            protocols: Vec::new(),
            front: Front::Document,
            settings: ViewSettings { dark, ..ViewSettings::default() },
            show_inspector: true,
            show_console: start.console,
            palette: if start.palette { PaletteState::opened() } else { PaletteState::default() },
            find: FindState::default(),
            goto: None,
            new_sequence: None,
            feature_dialog: None,
            digest_dialog: None,
            console: Console {
                code: "// DOC is the open document. Changes apply when the script finishes.\nprint(`${DOC.name}: ${DOC.length} bp, ${DOC.features.len()} features`);\n".into(),
                output: String::new(),
            },
            confirm_close: None,
            allow_quit: false,
            recent,
            notice: None,
            applied_dark: None,
            frames: 0,
            start: start.clone(),
            screenshot_done: Arc::new(AtomicBool::new(false)),
            screenshot_requested: false,
            play: start.play.clone(),
            recorder: start.record.clone().map(record::Recorder::new),
            updater,
        };
        for path in &start.files {
            app.open_path(path.clone());
        }
        if let Some(tab) = app.tabs.get_mut(app.active) {
            if let Some(view) = start.view {
                tab.view = view;
            }
            if let Some(range) = &start.select
                && let Ok(span) = slopgene_core::span::parse_range(range, tab.doc.len(), tab.doc.is_circular())
            {
                tab.select(Some(span));
                tab.scroll_to = Some(span.start);
            }
            if let Some(name) = &start.feature
                && let Some(f) = tab.doc.features.iter().find(|f| f.name.eq_ignore_ascii_case(name))
            {
                let uid = f.uid;
                tab.select_feature(uid);
            }
        }
        // `--align READ` aligns the read to the first document and shows the alignment view, and
        // `--select` then also brings the start of the range into view among the aligned bases.
        if let Some(read_path) = &start.align
            && let Some(tab) = app.tabs.get_mut(app.active)
        {
            tab.view = View::Alignment;
            let selected = tab.selection.map(|span| span.start);
            app.add_alignment_files(std::slice::from_ref(read_path));
            if let (Some(tab), Some(position)) = (app.tabs.get_mut(app.active), selected) {
                tab.alignments.reveal(position);
            }
        }
        if let Some(enzymes) = &start.digest {
            app.digest_dialog = Some(DigestDialog { enzymes: enzymes.clone(), result: None, lanes: Vec::new() });
            app.run_digest();
        }
        // `--new-protocol` starts a cloning protocol from the maps open, and `--enzymes` cuts its
        // insert out with those enzymes, or with `--pcr` puts them on the primers' tails.
        if start.new_protocol {
            app.new_cloning_protocol(start.pcr);
            if let (Some(enzymes), Front::Protocol(ii)) = (&start.enzymes, app.front) {
                let names: Vec<String> = enzymes.split([',', ' ', ';', '+']).filter(|s| !s.is_empty()).map(str::to_string).collect();
                let tab = &mut app.protocols[ii];
                let step = if start.pcr { "pcr_donor" } else { "cut_donor" };
                tab.change(None, |p| {
                    let _ = p.set_enzymes(step, names, true);
                });
                tab.refresh(&app.tabs, true);
            }
        }
        // `--step N` shows that step of the protocol in front. `--zoom` zooms into its card a few
        // frames in (see `ui`).
        if let (Some(step), Front::Protocol(ii)) = (start.step, app.front)
            && let Some(tab) = app.protocols.get_mut(ii)
        {
            tab.selected = step.saturating_sub(1).min(tab.protocol.steps.len().saturating_sub(1));
        }
        app
    }

    fn note(&mut self, text: impl Into<String>) {
        self.notice = Some((text.into(), Instant::now()));
    }

    /// The document tab in front, when a document is in front.
    fn tab(&mut self) -> Option<&mut DocTab> {
        if self.front != Front::Document {
            return None;
        }
        self.tabs.get_mut(self.active)
    }

    /// Whether a dialog that blocks the window (the palette, Go to, New sequence, Add feature,
    /// or the question before closing) is open.
    fn modal_open(&self) -> bool {
        self.palette.open || self.goto.is_some() || self.new_sequence.is_some() || self.feature_dialog.is_some() || self.confirm_close.is_some()
    }

    /// Paste clipboard text into the open document at the cursor, replacing the selection.
    fn paste(&mut self, text: &str) {
        let Some(tab) = self.tabs.get_mut(self.active) else { return };
        let unit = tab.unit();
        let outcome = tab.clean_input(text).and_then(|bases| tab.paste_bases(&bases).map(|_| bases.len()).map_err(|e| e.to_string()));
        match outcome {
            Ok(0) => {}
            Ok(count) => self.note(format!("pasted {} {unit}", thousands(count))),
            Err(e) => self.note(format!("not pasted: {e}")),
        }
    }

    pub fn open_path(&mut self, path: PathBuf) {
        if let Some(index) = self.tabs.iter().position(|t| t.path.as_deref() == Some(path.as_path())) {
            self.active = index;
            self.front = Front::Document;
            return;
        }
        if let Some(index) = self.protocols.iter().position(|t| t.path.as_deref() == Some(path.as_path())) {
            self.front = Front::Protocol(index);
            return;
        }
        // A protocol file opens in a protocol tab, and anything else as a sequence.
        let opened = if is_protocol_file(&path) {
            ProtocolTab::open(&path).map(|tab| self.add_protocol(tab))
        } else {
            DocTab::open(&path).map_err(|e| e.to_string()).map(|tab| {
                self.tabs.push(tab);
                self.active = self.tabs.len() - 1;
                self.front = Front::Document;
            })
        };
        match opened {
            Ok(()) => {
                self.recent.retain(|p| p != &path);
                self.recent.insert(0, path);
                self.recent.truncate(12);
            }
            Err(e) => self.note(format!("could not open {}: {e}", path.display())),
        }
    }

    /// Show a protocol in a tab of its own, in front.
    fn add_protocol(&mut self, tab: ProtocolTab) {
        self.protocols.push(tab);
        self.front = Front::Protocol(self.protocols.len() - 1);
    }

    fn open_dialog(&mut self) {
        let files = rfd::FileDialog::new()
            .add_filter("Sequence files", &["dna", "rna", "prot", "gb", "gbk", "genbank", "fasta", "fa", "fas", "json"])
            .add_filter("Protocols", &["json"])
            .add_filter("All files", &["*"])
            .pick_files();
        for path in files.unwrap_or_default() {
            self.open_path(path);
        }
    }

    fn save_as_dialog(&mut self, format: Option<Format>) {
        let Some(tab) = self.tabs.get_mut(self.active) else { return };
        let ext = format.unwrap_or(Format::SnapGene).extension(tab.doc.molecule);
        let suggested = format!("{}.{ext}", if tab.doc.name.is_empty() { "untitled" } else { &tab.doc.name });
        let mut dialog = rfd::FileDialog::new().set_file_name(suggested);
        dialog = match format {
            Some(Format::GenBank) => dialog.add_filter("GenBank", &["gb", "gbk"]),
            Some(Format::Fasta) => dialog.add_filter("FASTA", &["fasta", "fa"]),
            _ => dialog
                .add_filter("SnapGene", &["dna", "rna", "prot"])
                .add_filter("GenBank", &["gb", "gbk"])
                .add_filter("FASTA", &["fasta"])
                .add_filter("JSON", &["json"]),
        };
        if let Some(path) = dialog.save_file() {
            // The format the file is written in, and what of the document it cannot keep.
            let written = format.unwrap_or_else(|| io::format_of(&path));
            let lost = io::losses(&tab.doc, written);
            let result = match format {
                // Exports leave the document attached to its own file.
                Some(f) if f != Format::SnapGene => io::write_path(&tab.doc, &path, Some(f)).map(|_| ()),
                _ => tab.save_as(&path),
            };
            match result {
                Ok(()) if lost.is_empty() => self.note(format!("saved {}", path.display())),
                Ok(()) => self.note(format!("saved {} without {}, which {} cannot keep", path.display(), lost.join(", "), format_label(written))),
                Err(e) => self.note(format!("could not save: {e}")),
            }
        }
    }

    fn export_map(&mut self, png: bool) {
        let Some(tab) = self.tabs.get_mut(self.active) else { return };
        let ext = if png { "png" } else { "svg" };
        let Some(path) = rfd::FileDialog::new().set_file_name(format!("{}.{ext}", tab.doc.name)).add_filter(ext, &[ext]).save_file() else { return };
        let marks = if self.settings.show_enzymes { tab.sites().0.to_vec() } else { Vec::new() };
        let options = slopgene_core::render::MapOptions {
            width: 1000.0,
            height: 760.0,
            theme: if self.settings.dark { slopgene_core::render::Theme::dark() } else { slopgene_core::render::Theme::light() },
            enzymes: marks,
            colors: self.settings.colors,
            show_primers: self.settings.show_primers,
            ..Default::default()
        };
        let scene = slopgene_core::render::render_map(&tab.doc, &options, &slopgene_core::render::ApproxMeasure);
        match save_scene(&scene, &path, 2.0) {
            Ok(()) => self.note(format!("exported {}", path.display())),
            Err(e) => self.note(format!("could not export: {e}")),
        }
    }

    /// Write the bench of protocol `index` as an SVG or PNG picture: the picture the command
    /// line's `protocol render` draws, in the app's current colours.
    fn export_bench(&mut self, index: usize) {
        let Some(tab) = self.protocols.get(index) else { return };
        let name = if tab.protocol.name.is_empty() { "protocol".to_string() } else { tab.protocol.name.clone() };
        let Some(path) = rfd::FileDialog::new().set_file_name(format!("{name}.svg")).add_filter("svg", &["svg"]).add_filter("png", &["png"]).save_file() else {
            return;
        };
        let options = slopgene_core::render::bench::BenchOptions {
            colors: self.settings.colors,
            ..if self.settings.dark { slopgene_core::render::bench::BenchOptions::dark() } else { slopgene_core::render::bench::BenchOptions::default() }
        };
        let scene = slopgene_core::render::bench::render_bench(&tab.protocol, &tab.run, &options, &slopgene_core::render::ApproxMeasure);
        match save_scene(&scene, &path, 2.0) {
            Ok(()) => self.note(format!("exported {}", path.display())),
            Err(e) => self.note(format!("could not export: {e}")),
        }
    }

    /// Run an action from the palette, a shortcut or a button.
    fn perform(&mut self, action: Action, ctx: &egui::Context) {
        // With a protocol in front, the file and history commands act on it.
        if let Front::Protocol(index) = self.front
            && self.perform_on_protocol(action, index)
        {
            return;
        }
        match action {
            Action::Open => self.open_dialog(),
            Action::NewProtocol => self.new_cloning_protocol(false),
            Action::NewSequence => {
                self.new_sequence = Some(NewSequenceDialog { name: "untitled".into(), sequence: String::new(), circular: true, error: None })
            }
            Action::Save => {
                let name = |p: &std::path::Path| p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                match self.tab().map(|t| t.save_keeping_everything()) {
                    Some(Ok(Saved::InPlace(path))) => self.note(format!("saved {}", path.display())),
                    Some(Ok(Saved::Beside { path, original, format, lost })) => self.note(format!(
                        "saved {}: {} cannot keep {}, so {} is unchanged",
                        name(&path),
                        format_label(format),
                        lost.join(", "),
                        name(&original)
                    )),
                    Some(Ok(Saved::NeedsLocation { format, lost })) => {
                        self.note(format!("{} cannot keep {}: choose where to save a SnapGene file", format_label(format), lost.join(", ")));
                        self.save_as_dialog(None);
                    }
                    Some(Ok(Saved::NoFile)) => self.save_as_dialog(None),
                    Some(Err(e)) => self.note(format!("could not save: {e}")),
                    None => {}
                }
            }
            Action::SaveAs => self.save_as_dialog(None),
            Action::ExportGenBank => self.save_as_dialog(Some(Format::GenBank)),
            Action::ExportFasta => self.save_as_dialog(Some(Format::Fasta)),
            Action::ExportSvg => self.export_map(false),
            Action::ExportPng => self.export_map(true),
            Action::ExportBench => self.note("this command draws a cloning protocol: show its tab first"),
            Action::Close => self.request_close(self.active),
            Action::Reload => {
                if let Some(t) = self.tab()
                    && let Err(e) = t.reload()
                {
                    self.note(format!("could not reload: {e}"));
                }
            }
            Action::Undo => {
                if self.tab().is_some_and(|t| !t.undo()) {
                    self.note("nothing to undo");
                }
            }
            Action::Redo => {
                if self.tab().is_some_and(|t| !t.redo()) {
                    self.note("nothing to redo");
                }
            }
            Action::CopySelection | Action::CopyReverseComplement => match self.tab().map(|t| (t.selected_bases(), t.unit())) {
                Some((Some(bases), unit)) => {
                    let bases = if action == Action::CopyReverseComplement { alphabet::revcomp(&bases) } else { bases };
                    let count = bases.len();
                    ctx.copy_text(String::from_utf8_lossy(&bases).into_owned());
                    self.note(format!("copied {} {unit}", thousands(count)));
                }
                Some((None, unit)) => self.note(format!("select {unit} to copy")),
                None => {}
            },
            Action::Cut => {
                let outcome = self.tab().map(|t| (t.cut(), t.unit()));
                match outcome {
                    Some((Ok(Some(bases)), unit)) => {
                        ctx.copy_text(String::from_utf8_lossy(&bases).into_owned());
                        self.note(format!("cut {} {unit}", thousands(bases.len())));
                    }
                    Some((Ok(None), unit)) => self.note(format!("select {unit} to cut")),
                    Some((Err(e), _)) => self.note(format!("not cut: {e}")),
                    None => {}
                }
            }
            Action::SelectAll => {
                if let Some(t) = self.tab() {
                    t.select_all();
                }
            }
            Action::ViewMap | Action::ViewSequence | Action::ViewSplit => {
                if let Some(t) = self.tab() {
                    t.view = match action {
                        Action::ViewMap => View::Map,
                        Action::ViewSequence => View::Sequence,
                        _ => View::Split,
                    };
                }
            }
            Action::ToggleDark => self.settings.dark = !self.settings.dark,
            Action::ToggleColors => {
                self.settings.colors = if self.settings.colors == ColorMode::Palette { ColorMode::File } else { ColorMode::Palette };
            }
            Action::ToggleComplement => self.settings.show_complement = !self.settings.show_complement,
            Action::ToggleTranslations => self.settings.show_translations = !self.settings.show_translations,
            Action::ToggleEnzymes => self.settings.show_enzymes = !self.settings.show_enzymes,
            Action::TogglePrimers => self.settings.show_primers = !self.settings.show_primers,
            Action::ToggleInspector => self.show_inspector = !self.show_inspector,
            Action::ToggleConsole => self.show_console = !self.show_console,
            Action::Find => self.find.open = true,
            Action::GoTo => self.goto = Some(String::new()),
            Action::ReverseComplement => {
                if let Some(t) = self.tab() {
                    let _ = t.edit(|d| {
                        edit::reverse_complement(d);
                        Ok(())
                    });
                }
            }
            Action::SetOriginAtCursor => {
                let result = self.tab().map(|t| {
                    let at = t.cursor.min(t.len().saturating_sub(1));
                    t.edit(|d| edit::set_origin(d, at)).map(|_| t.cursor = 0)
                });
                if let Some(Err(e)) = result {
                    self.note(e.to_string());
                }
            }
            Action::ToggleTopology => {
                if let Some(t) = self.tab() {
                    let next = if t.doc.is_circular() { Topology::Linear } else { Topology::Circular };
                    let _ = t.edit(|d| {
                        edit::set_topology(d, next);
                        Ok(())
                    });
                }
            }
            Action::AddFeatureFromSelection => {
                if self.tab().is_some_and(|t| t.selection.is_some()) {
                    self.feature_dialog = Some(FeatureDialog { name: "new feature".into(), kind: "misc_feature".into(), strand: Strand::Forward });
                } else {
                    self.note("select a range first");
                }
            }
            Action::DeleteSelectedFeature => {
                if let Some(t) = self.tab()
                    && let Some(uid) = t.selected_feature
                {
                    let _ = t.edit(|d| {
                        d.features.retain(|f| f.uid != uid);
                        Ok(())
                    });
                    t.selected_feature = None;
                }
            }
            Action::PrimerFromSelection => {
                if let Some(t) = self.tab()
                    && let Some(sel) = t.selection
                {
                    let bases = String::from_utf8_lossy(&sel.extract(&t.doc.seq)).to_uppercase();
                    let name = format!("primer_{}", sel.start + 1);
                    let _ = t.edit(|d| {
                        d.primers.push(Primer::new(name, bases));
                        primers::refresh_binding_sites(d, &primers::BindingRules::default());
                        Ok(())
                    });
                }
            }
            Action::Digest => {
                let default = self
                    .tab()
                    .map(|t| t.sites().1.iter().filter(|c| c.cuts.len() == 1).take(2).map(|c| c.name.clone()).collect::<Vec<_>>().join(", "))
                    .unwrap_or_default();
                self.digest_dialog = Some(DigestDialog { enzymes: default, result: None, lanes: Vec::new() });
            }
            Action::Palette => self.palette = PaletteState::opened(),
            Action::Update => self.updater.open_dialog(),
            Action::AlignRead => self.open_alignment_view(),
            Action::DetectFeatures => self.detect_features_dialog(),
        }
    }

    /// Show the alignment view of the open document. With nothing aligned yet and no other DNA
    /// document open, the view could only add files, so the file dialog opens straight away.
    fn open_alignment_view(&mut self) {
        let Some(tab) = self.tabs.get_mut(self.active) else { return };
        if tab.doc.molecule != Molecule::Dna {
            self.note("alignment works on DNA sequences");
            return;
        }
        tab.view = View::Alignment;
        let (id, empty) = (tab.id, tab.alignments.entries.is_empty());
        let others = self.tabs.iter().any(|t| t.id != id && t.doc.molecule == Molecule::Dna);
        if empty && !others {
            self.pick_alignment_files();
        }
    }

    /// The other open DNA documents (tab id and title), which can be aligned to the open one.
    fn other_dna_documents(&self) -> Vec<(u64, String)> {
        let active = self.tabs.get(self.active).map(|t| t.id);
        self.tabs.iter().filter(|t| Some(t.id) != active && t.doc.molecule == Molecule::Dna).map(|t| (t.id, t.title())).collect()
    }

    /// Ask for sequence files and align every record in them to the open document.
    fn pick_alignment_files(&mut self) {
        let picked = rfd::FileDialog::new()
            .set_title("Sequences to align")
            .add_filter("Sequence files", &["gb", "gbk", "genbank", "fasta", "fa", "fas", "fna", "seq", "dna", "json"])
            .add_filter("All files", &["*"])
            .pick_files();
        if let Some(paths) = picked {
            self.add_alignment_files(&paths);
        }
    }

    /// Align every DNA record in these files to the open document.
    fn add_alignment_files(&mut self, paths: &[PathBuf]) {
        let Some(tab) = self.tabs.get_mut(self.active) else { return };
        let mut problems = Vec::new();
        for path in paths {
            match io::read_path(path) {
                Ok(docs) => {
                    for doc in docs {
                        if doc.molecule == Molecule::Dna {
                            let circular = doc.is_circular();
                            tab.alignments.add(&tab.doc, tab.revision, doc.name, doc.seq, circular);
                        } else {
                            problems.push(format!("{} is not DNA", doc.name));
                        }
                    }
                }
                Err(e) => problems.push(e.to_string()),
            }
        }
        if !problems.is_empty() {
            self.note(format!("not aligned: {}", problems.join("; ")));
        }
    }

    /// Align another open document to the open one.
    fn align_open_document(&mut self, id: u64) {
        let Some(other) = self.tabs.iter().find(|t| t.id == id) else { return };
        let (name, bases, circular) = (other.doc.name.clone(), other.doc.seq.clone(), other.doc.is_circular());
        if let Some(tab) = self.tabs.get_mut(self.active) {
            tab.alignments.add(&tab.doc, tab.revision, name, bases, circular);
        }
    }

    /// Carry out what the alignment view and its inspector asked for.
    fn handle_alignment_requests(&mut self, requests: align_view::Requests) {
        if let Some(id) = requests.add_tab {
            self.align_open_document(id);
        }
        if requests.pick_files {
            self.pick_alignment_files();
        }
    }

    fn detect_features_dialog(&mut self) {
        if self.tabs.is_empty() {
            return;
        }
        let Some(folder) = rfd::FileDialog::new().set_title("Folder of annotated maps to learn features from").pick_folder() else { return };
        let mut library = slopgene_core::annotate::FeatureLibrary::new();
        let mut stack = vec![folder];
        while let Some(dir) = stack.pop() {
            let Ok(entries) = std::fs::read_dir(&dir) else { continue };
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_dir() {
                    stack.push(path);
                } else if Format::from_path(&path).is_some()
                    && let Ok(docs) = io::read_path(&path)
                {
                    for doc in &docs {
                        library.add_document(doc, 30);
                    }
                }
            }
        }
        let tab = &mut self.tabs[self.active];
        let hits = slopgene_core::annotate::annotate(&tab.doc, &library, 0.95);
        let added = hits.len();
        if added > 0 {
            let _ = tab.edit(|d| {
                d.features.extend(hits.into_iter().map(|h| h.feature));
                Ok(())
            });
        }
        self.note(format!("library of {} features: added {added} to {}", library.len(), self.tabs[self.active].doc.name));
    }

    fn request_close(&mut self, index: usize) {
        if index >= self.tabs.len() {
            return;
        }
        if self.tabs[index].dirty {
            self.confirm_close = Some(CloseTarget::Document(index));
        } else {
            self.close_document(index);
        }
    }

    /// Close document `index`. When no document is left, an open protocol comes to the front.
    fn close_document(&mut self, index: usize) {
        self.tabs.remove(index);
        self.active = self.active.min(self.tabs.len().saturating_sub(1));
        if self.tabs.is_empty() && !self.protocols.is_empty() && self.front == Front::Document {
            self.front = Front::Protocol(0);
        }
    }

    fn request_close_protocol(&mut self, index: usize) {
        if index >= self.protocols.len() {
            return;
        }
        if self.protocols[index].dirty {
            self.confirm_close = Some(CloseTarget::Protocol(index));
        } else {
            self.close_protocol(index);
        }
    }

    /// Close protocol `index`. The tab in front stays in front; when it was the one closed, a
    /// document comes forward, or another protocol when no document is open.
    fn close_protocol(&mut self, index: usize) {
        self.protocols.remove(index);
        self.front = match self.front {
            Front::Protocol(front) if front > index => Front::Protocol(front - 1),
            Front::Protocol(front) if front < index => Front::Protocol(front),
            Front::Protocol(_) if self.tabs.is_empty() && !self.protocols.is_empty() => Front::Protocol(index.min(self.protocols.len() - 1)),
            _ => Front::Document,
        };
    }

    /// Do `action` on protocol `index`, which is in front. Returns false for the commands that
    /// do not act on the tab in front (opening files, view settings), which then run as usual.
    fn perform_on_protocol(&mut self, action: Action, index: usize) -> bool {
        if index >= self.protocols.len() {
            return false;
        }
        match action {
            Action::Save if self.protocols[index].path.is_some() => {
                let outcome = self.protocols[index].save();
                match outcome {
                    Ok(path) => self.note(format!("saved {}", path.display())),
                    Err(e) => self.note(format!("could not save: {e}")),
                }
            }
            Action::Save | Action::SaveAs => self.save_protocol_as(index),
            Action::Close => self.request_close_protocol(index),
            Action::Undo => {
                if !self.protocols[index].undo() {
                    self.note("nothing to undo");
                }
            }
            Action::Redo => {
                if !self.protocols[index].redo() {
                    self.note("nothing to redo");
                }
            }
            Action::Reload => {
                self.protocols[index].refresh(&self.tabs, true);
                self.note("read the maps again");
            }
            Action::ExportBench => self.export_bench(index),
            Action::Open
            | Action::NewSequence
            | Action::NewProtocol
            | Action::ToggleDark
            | Action::ToggleColors
            | Action::ToggleComplement
            | Action::ToggleTranslations
            | Action::ToggleEnzymes
            | Action::TogglePrimers
            | Action::ToggleInspector
            | Action::ToggleConsole
            | Action::Palette => return false,
            _ => self.note("this command works on a sequence: show its tab first"),
        }
        true
    }

    /// Ask where to save protocol `index`, and save it there.
    fn save_protocol_as(&mut self, index: usize) {
        let Some(tab) = self.protocols.get(index) else { return };
        // A file name from the protocol's name, without the characters file systems refuse.
        let stem: String = tab.name().chars().map(|c| if "<>:\"/\\|?*".contains(c) { '_' } else { c }).collect();
        let Some(path) = rfd::FileDialog::new().set_file_name(format!("{stem}.json")).add_filter("Protocol", &["json"]).save_file() else { return };
        let outcome = self.protocols[index].save_as(&path);
        match outcome {
            Ok(()) => self.note(format!("saved {}", path.display())),
            Err(e) => self.note(format!("could not save: {e}")),
        }
    }

    /// Start a cloning protocol in a tab of its own, on the bench, with no window in between. A
    /// document in front with a region marked carries the insert, as the region is what the
    /// person wants cut out; otherwise the document in front is the vector, unless it is the only
    /// map. The next open DNA map takes the other part. The region marked on each map becomes the
    /// part to cut out, or the part of the vector the insert replaces. The enzymes are chosen on
    /// the bench, so the protocol opens on the insert's cut.
    fn new_cloning_protocol(&mut self, pcr: bool) {
        let dna: Vec<usize> = (0..self.tabs.len()).filter(|&ii| self.tabs[ii].doc.molecule == Molecule::Dna).collect();
        let front = Some(self.active).filter(|ii| self.front == Front::Document && dna.contains(ii));
        let marked = front.filter(|&ii| self.tabs[ii].selection.is_some());
        let other = |taken: usize| dna.iter().copied().find(|&ii| ii != taken);
        let (donor, vector) = match (marked, front) {
            (Some(ii), _) => (Some(ii), other(ii)),
            (None, Some(ii)) => match other(ii) {
                Some(jj) => (Some(jj), Some(ii)),
                None => (Some(ii), None),
            },
            (None, None) => (dna.first().copied(), dna.get(1).copied()),
        };
        // Maps with a file are read from it; a map with none goes into the protocol itself.
        let here = std::env::current_dir().unwrap_or_default();
        let source = |ii: usize| {
            let tab = &self.tabs[ii];
            match &tab.path {
                Some(path) => OpenParams { path: Some(stored_path(&here, path)), document: None },
                None => OpenParams { path: None, document: Some(json::document_to_value(&tab.doc, true)) },
            }
        };
        // A region marked on a map, by the name of the feature it is when no other feature shares
        // that name, and by its bases otherwise.
        let part = |ii: usize| -> Option<String> {
            let tab = &self.tabs[ii];
            let span = tab.selection?;
            let len = tab.len();
            let feature = tab.selected_feature.and_then(|uid| tab.doc.feature_index(uid)).map(|index| &tab.doc.features[index]);
            match feature {
                Some(f) if !f.name.is_empty() && f.span(len).normalize(len) == span.normalize(len) && tab.doc.features_named(&f.name).count() == 1 => {
                    Some(f.name.clone())
                }
                _ => Some(span.format(len)),
            }
        };
        let protocol = match donor {
            Some(donor) => {
                let insert = part(donor);
                let what = insert.clone().unwrap_or_else(|| self.tabs[donor].doc.name.clone());
                let name = match vector {
                    Some(v) => format!("{what} into {}", self.tabs[v].doc.name),
                    None => what,
                };
                Protocol::restriction_cloning(&CloningTemplate {
                    name,
                    vector: vector.map(source),
                    donor: source(donor),
                    insert,
                    replace: vector.and_then(part),
                    dephosphorylate: true,
                    pcr,
                    ..CloningTemplate::default()
                })
            }
            None => Protocol { name: "New cloning protocol".into(), ..Protocol::default() },
        };
        let mut tab = ProtocolTab::new(protocol, None);
        tab.dirty = true;
        tab.refresh(&self.tabs, true);
        tab.selected = tab.protocol.steps.iter().position(|s| s.id == if pcr { "pcr_donor" } else { "cut_donor" }).unwrap_or(0);
        self.add_protocol(tab);
    }

    /// Carry out what a protocol tab asked for: open molecules in tabs of their own, or read the
    /// maps from disk again.
    fn handle_protocol_requests(&mut self, requests: protocol_view::Requests) {
        for doc in requests.open {
            let mut tab = DocTab::new(doc, None);
            tab.dirty = true;
            self.tabs.push(tab);
            self.active = self.tabs.len() - 1;
            self.front = Front::Document;
        }
        if requests.rerun
            && let Front::Protocol(ii) = self.front
            && ii < self.protocols.len()
        {
            self.protocols[ii].refresh(&self.tabs, true);
            self.note("read the maps again");
        }
    }

    fn run_digest(&mut self) {
        let Some(dialog) = self.digest_dialog.as_mut() else { return };
        let Some(tab) = self.tabs.get(self.active) else { return };
        let db = EnzymeDb::builtin();
        let names: Vec<&str> = dialog.enzymes.split([',', ' ', ';']).filter(|s| !s.is_empty()).collect();
        let result = db.resolve(&names).map_err(|e| e.to_string()).and_then(|enzymes| digest(&tab.doc, &enzymes).map_err(|e| e.to_string()));
        match result {
            Ok(d) => {
                let sizes: Vec<usize> = d.fragments.iter().map(|f| f.len()).collect();
                let label = names.join("+");
                dialog.lanes = vec![
                    Lane::ladder(Ladder::get("1kb").unwrap()),
                    if tab.doc.is_circular() { Lane::uncut_plasmid("uncut", tab.len()) } else { Lane::digest("uncut", &[tab.len()]) },
                    if d.fragments.len() == 1 && d.fragments[0].doc.is_circular() {
                        Lane::uncut_plasmid(label, tab.len())
                    } else {
                        Lane::digest(label, &sizes)
                    },
                ];
                let fragments = d
                    .fragments
                    .into_iter()
                    .map(|f| {
                        let ends = format!("{} / {}", f.doc.ends.left.describe(), f.doc.ends.right.describe());
                        (f.len(), ends, f.doc)
                    })
                    .collect();
                dialog.result = Some(Ok(fragments));
            }
            Err(e) => {
                dialog.result = Some(Err(e));
                dialog.lanes.clear();
            }
        }
    }

    fn run_console(&mut self) {
        if self.front != Front::Document {
            self.console.output.push_str("show a sequence tab first\n");
            return;
        }
        let Some(tab) = self.tabs.get_mut(self.active) else {
            self.console.output.push_str("open a document first\n");
            return;
        };
        let buffer = Arc::new(Mutex::new(String::new()));
        let options = slopgene_script::RunOptions {
            args: Vec::new(),
            output: Some(buffer.clone()),
            document: Some(slopgene_script::Doc { doc: tab.doc.clone(), path: tab.path.clone() }),
        };
        let started = Instant::now();
        let (result, doc_after) = slopgene_script::run_with_document(&self.console.code, options);
        let printed = buffer.lock().map(|b| b.clone()).unwrap_or_default();
        self.console.output.push_str(&printed);
        match result {
            Ok(value) if !value.is_unit() => self.console.output.push_str(&format!("→ {value}\n")),
            Ok(_) => {}
            Err(e) => self.console.output.push_str(&format!("error: {e}\n")),
        }
        if let Some(after) = doc_after
            && after.doc != tab.doc
        {
            tab.replace_document(after.doc);
            self.console.output.push_str("(document updated; undo with Ctrl+Z)\n");
        }
        self.console.output.push_str(&format!("— {:.0} ms\n", started.elapsed().as_secs_f64() * 1000.0));
    }

    /// Refresh find results when the query or the document changed.
    fn refresh_find(&mut self) {
        let Some(tab) = self.tabs.get(self.active) else { return };
        let key = (tab.id, tab.revision, self.find.query.clone(), self.find.mismatches, self.find.protein);
        if self.find.key.as_ref() == Some(&key) {
            return;
        }
        self.find.key = Some(key);
        self.find.current = 0;
        let query = self.find.query.trim();
        self.find.hits = if query.is_empty() {
            Vec::new()
        } else if self.find.protein {
            alphabet::clean_protein(query).map(|m| find_protein(&tab.doc.seq, tab.doc.is_circular(), &m, GeneticCode::standard())).unwrap_or_default()
        } else {
            alphabet::clean_nucleotides(query).map(|p| find_dna(&tab.doc.seq, tab.doc.is_circular(), &p, self.find.mismatches, true)).unwrap_or_default()
        };
    }

    fn jump_to_hit(&mut self, step: isize) {
        if self.find.hits.is_empty() {
            return;
        }
        let n = self.find.hits.len() as isize;
        self.find.current = ((self.find.current as isize + step).rem_euclid(n)) as usize;
        let span = self.find.hits[self.find.current].span;
        if let Some(tab) = self.tabs.get_mut(self.active) {
            tab.select(Some(span));
            tab.scroll_to = Some(span.start);
        }
    }
}

/// `text` shortened from the front, behind a leading "…", so that it fits in `width` points at
/// the interface font of `size`. Text that fits is returned whole.
fn elide_front(ui: &egui::Ui, text: &str, width: f32, size: f32) -> String {
    let measure = |s: &str| ui.ctx().fonts_mut(|f| f.layout_no_wrap(s.to_string(), style::ui_font(size), egui::Color32::WHITE).size().x);
    if width <= 0.0 {
        return String::new();
    }
    if measure(text) <= width {
        return text.to_string();
    }
    // Dropping more leading characters only ever makes the text narrower, so binary-search for
    // the fewest characters to drop. `starts` holds the byte offset of every character.
    let starts: Vec<usize> = text.char_indices().map(|(offset, _)| offset).collect();
    let tail = |dropped: usize| format!("…{}", &text[starts.get(dropped).copied().unwrap_or(text.len())..]);
    let (mut lo, mut hi) = (1, starts.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if measure(&tail(mid)) <= width { hi = mid } else { lo = mid + 1 }
    }
    tail(lo)
}

/// Whether a file is a SlopGene protocol: JSON that names the `slopgene.protocol` schema, or
/// `fukgene.protocol`, as protocols saved before the program was renamed do.
fn is_protocol_file(path: &std::path::Path) -> bool {
    if !path.extension().is_some_and(|e| e.eq_ignore_ascii_case("json")) {
        return false;
    }
    let Ok(text) = std::fs::read_to_string(path) else { return false };
    let named = ["\"slopgene.protocol/", "\"fukgene.protocol/"].iter().any(|schema| text.contains(schema));
    named && Protocol::from_json(&text).is_ok()
}

/// Save a scene as SVG, or PNG when the path ends in .png.
fn save_scene(scene: &slopgene_core::render::Scene, path: &std::path::Path, scale: f32) -> std::result::Result<(), String> {
    slopgene_script::files::save_scene(scene, path, scale)
}

/// Encode an egui screenshot as PNG.
fn save_screenshot(image: &egui::ColorImage, path: &std::path::Path) -> std::result::Result<(), String> {
    let [w, h] = image.size;
    let mut data = Vec::with_capacity(w * h * 4);
    for p in &image.pixels {
        data.extend_from_slice(&[p.r(), p.g(), p.b(), p.a()]);
    }
    let png = slopgene_script::files::rgba_png(w as u32, h as u32, data)?;
    io::write_atomic(path, &png).map_err(|e| e.to_string())
}

/// Score how well `query` matches `text` (all query characters in order); higher is better.
fn fuzzy_score(query: &str, text: &str) -> Option<i32> {
    let q: Vec<char> = query.to_lowercase().chars().filter(|c| !c.is_whitespace()).collect();
    if q.is_empty() {
        return Some(0);
    }
    let t: Vec<char> = text.to_lowercase().chars().collect();
    let mut score = 0;
    let mut ti = 0;
    let mut last = None;
    for qc in q {
        let found = (ti..t.len()).find(|&ii| t[ii] == qc)?;
        // Consecutive matches and matches at word starts score higher.
        score += if last == Some(found.wrapping_sub(1)) { 6 } else { 1 };
        if found == 0 || t[found - 1] == ' ' {
            score += 4;
        }
        last = Some(found);
        ti = found + 1;
    }
    Some(score - t.len() as i32 / 8)
}

impl eframe::App for SlopGeneApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.frames += 1;
        // A picture or a recording keeps the size asked for on the command line. eframe opens the
        // window at the size it had last, so the size is set again once the window is up.
        if self.frames == 2 && self.start.headless() && self.start.size[0] > 0.0 {
            ctx.send_viewport_cmd(egui::ViewportCommand::InnerSize(vec2(self.start.size[0], self.start.size[1])));
        }
        // `--zoom` zooms into the step shown on the third frame, as a click would, so a screenshot
        // can catch the zoom on its way in as well as at its end.
        if self.start.zoom
            && self.frames == 3
            && let Front::Protocol(ii) = self.front
            && let Some(tab) = self.protocols.get_mut(ii)
        {
            let step = tab.selected;
            tab.zoom_to(step);
        }
        if self.applied_dark != Some(self.settings.dark) {
            style::apply(&ctx, self.settings.dark);
            self.applied_dark = Some(self.settings.dark);
        }
        let palette = Palette::get(self.settings.dark);

        // Files dropped onto the window.
        let dropped: Vec<PathBuf> = ctx.input(|i| i.raw.dropped_files.iter().map(|f| f.path().to_path_buf()).collect());
        for path in dropped {
            self.open_path(path);
        }

        // Edits made outside the app show up live.
        for ii in 0..self.tabs.len() {
            self.tabs[ii].poll_disk();
        }
        ctx.request_repaint_after(Duration::from_secs(1));

        // Shortcuts, with the editing keys going to the document unless a text field or a
        // dialog has them (see `document_has_keys`).
        let requests = read_shortcuts(&ctx, document_has_keys(&ctx, self.modal_open()));
        let mut actions: Vec<Action> = requests.actions;
        for text in requests.pastes {
            self.paste(&text);
        }

        // Closing the window with unsaved work asks first.
        let unsaved = self.tabs.iter().any(|t| t.dirty) || self.protocols.iter().any(|t| t.dirty);
        if ctx.input(|i| i.viewport().close_requested()) && !self.allow_quit && unsaved {
            ctx.send_viewport_cmd(egui::ViewportCommand::CancelClose);
            self.confirm_close = Some(CloseTarget::App);
        }

        // The protocol in front runs again when it, or an open map it may use, changed.
        if let Front::Protocol(ii) = self.front
            && let Some(protocol) = self.protocols.get_mut(ii)
        {
            protocol.refresh(&self.tabs, false);
        }

        self.top_bar(ui, &palette, &mut actions);
        self.status_bar(ui, &palette);
        if self.show_console {
            self.console_panel(ui, &palette);
        }
        // What the alignment view and its inspector ask for, done once both are drawn. The other
        // open DNA documents are the ones they can offer for alignment.
        let mut alignment_requests = align_view::Requests::default();
        let mut protocol_requests = protocol_view::Requests::default();
        let others = self.other_dna_documents();
        // A protocol draws its own panel beside its bench, so the inspector is for documents.
        if self.show_inspector && self.front == Front::Document && !self.tabs.is_empty() {
            let settings = self.settings.clone();
            let active = self.active;
            egui::Panel::right("inspector")
                .resizable(true)
                .default_size(300.0)
                .min_size(220.0)
                .frame(egui::Frame::new().fill(palette.panel).inner_margin(egui::Margin::symmetric(14, 10)))
                .show(ui, |ui| {
                    if let Some(tab) = self.tabs.get_mut(active) {
                        if tab.view == View::Alignment {
                            alignment_requests = align_view::inspector(ui, tab, &others, &settings);
                        } else {
                            panels::inspector(ui, tab, &settings, &mut actions);
                        }
                    }
                });
        }
        let keys = !self.modal_open();
        egui::CentralPanel::default().frame(egui::Frame::new().fill(palette.canvas)).show(ui, |ui| {
            if let Front::Protocol(ii) = self.front
                && let Some(protocol) = self.protocols.get_mut(ii)
            {
                protocol_requests = protocol_view::show(ui, protocol, &self.tabs, &self.settings, keys, self.show_inspector);
                return;
            }
            if self.tabs.is_empty() {
                self.welcome(ui, &palette, &mut actions);
                return;
            }
            if self.find.open {
                self.find_bar(ui, &palette);
            }
            let settings = self.settings.clone();
            let tab = &mut self.tabs[self.active];
            if tab.changed_on_disk {
                ui.horizontal(|ui| {
                    ui.add_space(12.0);
                    ui.label(RichText::new("This file changed on disk while you have unsaved edits.").color(palette.warn));
                    if ui.button("Reload theirs").clicked() {
                        let _ = tab.reload();
                    }
                    if ui.button("Keep mine").clicked() {
                        tab.changed_on_disk = false;
                    }
                });
            }
            match tab.view {
                View::Map => match map_view::show(ui, tab, &settings) {
                    Some(map_view::MapEvent::EditFeature(uid)) => {
                        tab.select_feature(uid);
                        self.show_inspector = true;
                    }
                    Some(map_view::MapEvent::Act(action)) => actions.push(action),
                    None => {}
                },
                View::Sequence => seq_view::show(ui, tab, &settings),
                View::Split => {
                    let width = ui.available_width();
                    egui::Panel::left(("split-map", tab.id))
                        .resizable(true)
                        .default_size(width * 0.42)
                        .min_size(240.0)
                        .frame(egui::Frame::new().fill(palette.canvas))
                        .show(ui, |ui| {
                            if let Some(map_view::MapEvent::Act(action)) = map_view::show(ui, tab, &settings) {
                                actions.push(action);
                            }
                        });
                    seq_view::show(ui, tab, &settings);
                }
                View::Alignment => {
                    let requests = align_view::show(ui, tab, &settings, &others);
                    alignment_requests = std::mem::take(&mut alignment_requests).and(requests);
                }
            }
        });
        self.handle_alignment_requests(alignment_requests);
        self.handle_protocol_requests(protocol_requests);

        self.palette_window(&ctx, &palette, &mut actions);
        self.dialogs(&ctx, &palette);
        for action in actions {
            self.perform(action, &ctx);
        }

        // Window title.
        let title = match self.front {
            Front::Protocol(ii) => self.protocols.get(ii).map(|t| t.title()),
            Front::Document => self.tabs.get(self.active).map(|t| t.title()),
        };
        let version = env!("CARGO_PKG_VERSION");
        let title = title.map(|t| format!("{t} · SlopGene {version}")).unwrap_or_else(|| format!("SlopGene {version}"));
        ctx.send_viewport_cmd(egui::ViewportCommand::Title(title));

        // The cursor of a played script, or of a `--drag`, drawn over everything. A script keeps
        // the frames coming while it plays.
        if let Some(play) = &self.play {
            play.draw(&ctx, palette.accent);
            ctx.request_repaint();
        } else if self.start.drag.is_some()
            && let Some(pointer) = ctx.input(|i| i.pointer.latest_pos())
        {
            record::draw_cursor(&ctx, pointer, &[], 0.0, palette.accent);
        }

        // Screenshot mode: capture once the layout has settled and the script, if any, has
        // played, then quit.
        let played = self.play.as_ref().is_none_or(record::Play::finished);
        if let Some(path) = self.start.screenshot.clone() {
            if self.frames >= self.start.screenshot_after && played && !self.screenshot_requested {
                self.screenshot_requested = true;
                let done = self.screenshot_done.clone();
                ctx.request_screenshot(move |image| {
                    if let Err(e) = save_screenshot(&image, &path) {
                        eprintln!("screenshot failed: {e}");
                    }
                    done.store(true, Ordering::SeqCst);
                });
            }
            if self.screenshot_done.load(Ordering::SeqCst) {
                self.allow_quit = true;
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            }
            ctx.request_repaint();
        }

        // Recording: ask for a frame as the window draws, and quit once the file is written.
        if let Some(recorder) = &mut self.recorder {
            let time = ctx.input(|i| i.time);
            let stop = played && self.frames >= self.start.screenshot_after;
            recorder.tick(&ctx, time, self.frames, stop);
            if recorder.finished() && !self.allow_quit {
                eprintln!("{}", recorder.summary());
                self.allow_quit = true;
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            }
            ctx.request_repaint();
        }
    }

    /// Pictures and recordings leave the saved layout of real sessions alone.
    fn persist_egui_memory(&self) -> bool {
        !self.start.headless()
    }

    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        // A picture or a recording (perhaps with --dark) must not change the user's saved settings.
        if self.start.headless() {
            return;
        }
        let recent: Vec<String> = self.recent.iter().map(|p| p.to_string_lossy().into_owned()).collect();
        storage.set_string("recent", recent.join("\n"));
        storage.set_string("dark", if self.settings.dark { "1".into() } else { "0".into() });
        storage.set_string("check_updates", if self.updater.at_start { "1".into() } else { "0".into() });
    }

    fn clear_color(&self, _visuals: &egui::Visuals) -> [f32; 4] {
        Palette::get(self.settings.dark).canvas.to_normalized_gamma_f32()
    }

    /// `--play` feeds its script's pointer and key events into the window frame by frame.
    ///
    /// `--drag` plays a drag into the window for a screenshot: the pointer arrives at the first
    /// point, presses there, moves to the second point over 12 frames, and holds, so the picture
    /// shows the drag on its way. The drag ends a few frames before the screenshot, so anything
    /// the window does first, such as a zoom, has finished before the button goes down.
    fn raw_input_hook(&mut self, _ctx: &egui::Context, raw_input: &mut egui::RawInput) {
        if let Some(play) = &mut self.play {
            play.tick(raw_input);
        }
        let Some((from, to)) = self.start.drag else { return };
        let start = self.start.screenshot_after.saturating_sub(16).max(4);
        let button = egui::Event::PointerButton { pos: from, button: egui::PointerButton::Primary, pressed: true, modifiers: Modifiers::default() };
        match self.frames.checked_sub(start) {
            Some(0) => raw_input.events.push(egui::Event::PointerMoved(from)),
            Some(1) => raw_input.events.push(button),
            Some(step @ 2..=13) => raw_input.events.push(egui::Event::PointerMoved(from.lerp(to, (step - 1) as f32 / 12.0))),
            _ => {}
        }
    }
}

impl SlopGeneApp {
    fn top_bar(&mut self, ui: &mut egui::Ui, palette: &Palette, actions: &mut Vec<Action>) {
        egui::Panel::top("top-bar").exact_size(42.0).frame(egui::Frame::new().fill(palette.panel).inner_margin(egui::Margin::symmetric(12, 6))).show(
            ui,
            |ui| {
                ui.horizontal_centered(|ui| {
                    ui.label(RichText::new("SlopGene").font(style::bold_font(14.0)).color(palette.accent))
                        .on_hover_text(format!("SlopGene {}", env!("CARGO_PKG_VERSION")));
                    ui.add_space(10.0);
                    // Tabs: the documents, then the protocols, which are set in italics.
                    let mut close = None;
                    let mut activate = None;
                    let documents =
                        self.tabs.iter().enumerate().map(|(ii, t)| (Front::Document, ii, t.title(), self.front == Front::Document && ii == self.active));
                    let protocols = self.protocols.iter().enumerate().map(|(ii, t)| (Front::Protocol(ii), ii, t.title(), self.front == Front::Protocol(ii)));
                    for (kind, ii, title, active) in documents.chain(protocols).collect::<Vec<_>>() {
                        let mut text = RichText::new(title).color(if active { palette.ink } else { palette.muted });
                        if kind != Front::Document {
                            text = text.italics();
                        }
                        let response = ui.selectable_label(active, text);
                        let response = if kind == Front::Document { response } else { response.on_hover_text("A cloning protocol") };
                        if response.clicked() {
                            activate = Some((kind, ii));
                        }
                        if response.middle_clicked() {
                            close = Some((kind, ii));
                        }
                        if active && ui.small_button("×").on_hover_text("Close (Ctrl+W)").clicked() {
                            close = Some((kind, ii));
                        }
                    }
                    let ctrl = if cfg!(target_os = "macos") { "⌘" } else { "Ctrl+" };
                    ui.menu_button("+", |ui| {
                        for (label, keys, action) in [
                            ("Open a file…", "O", Action::Open),
                            ("New sequence…", "N", Action::NewSequence),
                            ("New cloning protocol", "Shift+N", Action::NewProtocol),
                        ] {
                            if ui.button(format!("{label}   {ctrl}{keys}")).clicked() {
                                actions.push(action);
                                ui.close();
                            }
                        }
                    });
                    if let Some((kind, ii)) = activate {
                        if kind == Front::Document {
                            self.active = ii;
                        }
                        self.front = kind;
                    }
                    match close {
                        Some((Front::Document, ii)) => self.request_close(ii),
                        Some((Front::Protocol(_), ii)) => self.request_close_protocol(ii),
                        None => {}
                    }

                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        let theme = if self.settings.dark { "Light" } else { "Dark" };
                        if ui.small_button(theme).on_hover_text("Toggle dark mode").clicked() {
                            actions.push(Action::ToggleDark);
                        }
                        let shortcut = if cfg!(target_os = "macos") { "⌘K" } else { "Ctrl K" };
                        if ui.small_button(format!("Commands  {shortcut}")).clicked() {
                            actions.push(Action::Palette);
                        }
                        // Update is always there. It takes the accent colour once the check has
                        // found a newer version, and its dialog says the rest.
                        ui.add_space(4.0);
                        let (label, hint) = match self.updater.available() {
                            Some(version) => {
                                (RichText::new("Update").color(palette.accent), format!("SlopGene {version} is released. Download and install it."))
                            }
                            None => (RichText::new("Update"), format!("Look for a newer version of SlopGene. You have {}.", env!("CARGO_PKG_VERSION"))),
                        };
                        if ui.small_button(label).on_hover_text(hint).clicked() {
                            actions.push(Action::Update);
                        }
                        // The tools sit apart from the views: Clone, then Align to its left.
                        // Clone starts a cloning protocol whatever tab is in front.
                        ui.add_space(8.0);
                        let clone = ui.selectable_label(false, "Clone").on_hover_text(
                            "Start a cloning protocol: cut a part out of one plasmid and put it into another, step by step on a bench (Ctrl+Shift+N)",
                        );
                        if clone.clicked() {
                            actions.push(Action::NewProtocol);
                        }
                        if self.front == Front::Document
                            && let Some(tab) = self.tabs.get(self.active)
                        {
                            ui.add_space(2.0);
                            let align = ui
                                .selectable_label(tab.view == View::Alignment, "Align")
                                .on_hover_text("Compare sequencing results or other sequences with this one (Ctrl+Shift+A)");
                            if align.clicked() {
                                actions.push(Action::AlignRead);
                            }
                            ui.add_space(6.0);
                            for (view, label, action) in [
                                (View::Split, "Split", Action::ViewSplit),
                                (View::Sequence, "Sequence", Action::ViewSequence),
                                (View::Map, "Map", Action::ViewMap),
                            ] {
                                if ui.selectable_label(tab.view == view, label).clicked() {
                                    actions.push(action);
                                }
                            }
                        }
                    });
                });
            },
        );
    }

    fn status_bar(&mut self, ui: &mut egui::Ui, palette: &Palette) {
        egui::Panel::bottom("status-bar").exact_size(26.0).frame(egui::Frame::new().fill(palette.panel).inner_margin(egui::Margin::symmetric(12, 4))).show(
            ui,
            |ui| {
                ui.horizontal_centered(|ui| {
                    let small = |text: String| RichText::new(text).size(11.5).color(palette.muted);
                    // A notice is cut short rather than run under the labels on the right, which
                    // keep this much room.
                    let notice_label = |ui: &mut egui::Ui, text: String| {
                        let room = ui.available_width() - 300.0;
                        if room > 40.0 {
                            ui.add_space(16.0);
                            ui.scope(|ui| {
                                ui.set_max_width(room - 16.0);
                                ui.add(egui::Label::new(RichText::new(text).size(11.5).color(palette.accent)).truncate());
                            });
                        }
                    };
                    if let Front::Protocol(ii) = self.front
                        && let Some(tab) = self.protocols.get(ii)
                    {
                        ui.label(small(tab.progress()));
                        let notice = tab.notice.clone().or_else(|| self.notice.clone());
                        if let Some((text, at)) = notice
                            && at.elapsed() < Duration::from_secs(6)
                        {
                            notice_label(ui, text);
                        }
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            ui.label(small(format!("protocol  ·  {}", if tab.dirty { "unsaved" } else { "saved" })));
                            let place = tab.path.as_ref().map(|p| p.display().to_string()).unwrap_or_else(|| "not saved yet".into());
                            let room = ui.available_width() - 24.0;
                            ui.label(small(format!("{}  ·", elide_front(ui, &place, room, 11.5))));
                        });
                    } else if let Some(tab) = self.tabs.get(self.active) {
                        let len = tab.len();
                        let left = match tab.selection {
                            Some(sel) => {
                                let bases = sel.extract(&tab.doc.seq);
                                format!("{}  ·  {} bp  ·  GC {:.0}%", sel.format(len), thousands(bases.len()), alphabet::gc_fraction(&bases) * 100.0)
                            }
                            None => format!("cursor {}", thousands(tab.cursor.min(len))),
                        };
                        ui.label(small(left));
                        let notice = tab.notice.clone().or_else(|| self.notice.clone());
                        if let Some((text, at)) = notice
                            && at.elapsed() < Duration::from_secs(6)
                        {
                            notice_label(ui, text);
                        }
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            let state = if tab.dirty { "unsaved" } else { "saved" };
                            let unit = if tab.doc.molecule == Molecule::Protein { "aa" } else { "bp" };
                            ui.label(small(format!("{} {unit}  ·  {}  ·  {state}", thousands(len), tab.doc.topology.name())));
                            // The path gets the room that is left, shortened from the front so the
                            // file name stays visible in a narrow window.
                            let place = tab.path.as_ref().map(|p| p.display().to_string()).unwrap_or_else(|| "not saved yet".into());
                            let room = ui.available_width() - 24.0;
                            ui.label(small(format!("{}  ·", elide_front(ui, &place, room, 11.5))));
                        });
                    } else if let Some((text, at)) = self.notice.clone()
                        && at.elapsed() < Duration::from_secs(6)
                    {
                        ui.label(RichText::new(text).size(11.5).color(palette.accent));
                    }
                });
            },
        );
    }

    fn console_panel(&mut self, ui: &mut egui::Ui, palette: &Palette) {
        let mut run = false;
        egui::Panel::bottom("console")
            .resizable(true)
            .default_size(220.0)
            .min_size(120.0)
            .frame(egui::Frame::new().fill(palette.panel).inner_margin(egui::Margin::symmetric(12, 8)))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new("SCRIPT").font(style::bold_font(10.5)).color(palette.muted));
                    ui.label(RichText::new("Rhai · DOC is the open document · Ctrl+Enter runs").small().color(palette.faint));
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.small_button("Clear output").clicked() {
                            self.console.output.clear();
                        }
                        if ui.button("Run").clicked() {
                            run = true;
                        }
                    });
                });
                ui.columns(2, |cols| {
                    let editor = egui::TextEdit::multiline(&mut self.console.code).code_editor().desired_width(f32::INFINITY).desired_rows(8);
                    let response = egui::ScrollArea::vertical().id_salt("console-code").show(&mut cols[0], |ui| ui.add(editor)).inner;
                    if response.has_focus() && cols[0].input_mut(|i| i.consume_shortcut(&KeyboardShortcut::new(Modifiers::COMMAND, Key::Enter))) {
                        run = true;
                    }
                    egui::ScrollArea::vertical().id_salt("console-output").stick_to_bottom(true).auto_shrink([false, false]).show(&mut cols[1], |ui| {
                        ui.label(RichText::new(&self.console.output).font(style::mono_font(12.0)));
                    });
                });
            });
        if run {
            self.run_console();
        }
    }

    fn find_bar(&mut self, ui: &mut egui::Ui, palette: &Palette) {
        let mut close = false;
        let mut step = 0isize;
        egui::Frame::new().fill(palette.panel).inner_margin(egui::Margin::symmetric(12, 6)).show(ui, |ui| {
            ui.horizontal(|ui| {
                let response = ui.add(egui::TextEdit::singleline(&mut self.find.query).hint_text("find bases (IUPAC) or protein").desired_width(260.0));
                if !response.has_focus() && self.find.query.is_empty() {
                    response.request_focus();
                }
                if response.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                    step = if ui.input(|i| i.modifiers.shift) { -1 } else { 1 };
                    response.request_focus();
                }
                ui.checkbox(&mut self.find.protein, "protein");
                if !self.find.protein {
                    ui.add(egui::DragValue::new(&mut self.find.mismatches).range(0..=3).prefix("mismatches "));
                }
                let count = self.find.hits.len();
                let status = if self.find.query.trim().is_empty() {
                    String::new()
                } else if count == 0 {
                    "no matches".into()
                } else {
                    format!("{} of {count}", self.find.current + 1)
                };
                ui.label(RichText::new(status).color(palette.muted));
                if ui.small_button("↑").clicked() {
                    step = -1;
                }
                if ui.small_button("↓").clicked() {
                    step = 1;
                }
                if ui.small_button("×").clicked() || ui.input(|i| i.key_pressed(Key::Escape)) {
                    close = true;
                }
            });
        });
        self.refresh_find();
        if step != 0 {
            self.jump_to_hit(step);
        }
        if close {
            self.find.open = false;
        }
    }

    fn welcome(&mut self, ui: &mut egui::Ui, palette: &Palette, actions: &mut Vec<Action>) {
        let rect = ui.max_rect();
        let mut child = ui.new_child(egui::UiBuilder::new().max_rect(egui::Rect::from_center_size(rect.center(), vec2(420.0, 420.0))));
        child.vertical_centered(|ui| {
            ui.add_space(30.0);
            ui.label(RichText::new("SlopGene").font(style::bold_font(30.0)).color(palette.ink));
            ui.label(RichText::new("Plasmid maps, sequences and cloning, with a scriptable core.").color(palette.muted));
            ui.add_space(22.0);
            if ui.add_sized([220.0, 32.0], egui::Button::new("Open a file…")).clicked() {
                actions.push(Action::Open);
            }
            if ui.add_sized([220.0, 32.0], egui::Button::new("New sequence…")).clicked() {
                actions.push(Action::NewSequence);
            }
            if ui.add_sized([220.0, 32.0], egui::Button::new("New cloning protocol")).clicked() {
                actions.push(Action::NewProtocol);
            }
            ui.add_space(8.0);
            ui.label(RichText::new("or drop .dna, .gb or .fasta files here").small().color(palette.faint));
            if !self.recent.is_empty() {
                ui.add_space(22.0);
                ui.label(RichText::new("RECENT").font(style::bold_font(10.5)).color(palette.muted));
                let mut open = None;
                for path in self.recent.iter().take(8) {
                    let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                    if ui
                        .add(egui::Label::new(RichText::new(name).color(palette.accent)).sense(Sense::click()))
                        .on_hover_text(path.display().to_string())
                        .clicked()
                    {
                        open = Some(path.clone());
                    }
                }
                if let Some(p) = open {
                    self.open_path(p);
                }
            }
        });
    }

    fn palette_window(&mut self, ctx: &egui::Context, palette: &Palette, actions: &mut Vec<Action>) {
        if !self.palette.open {
            return;
        }
        let mut scored: Vec<(i32, Action)> = Action::ALL.iter().filter_map(|&a| fuzzy_score(&self.palette.query, a.label()).map(|s| (s, a))).collect();
        scored.sort_by_key(|s| std::cmp::Reverse(s.0));
        // With no query every command shows, in its group under a heading. A query lists the
        // matches best first instead.
        let grouped = self.palette.query.trim().is_empty();
        let count = scored.len();
        let mut chosen = None;
        let mut close = false;
        // The list may take most of the window's height.
        let list_height = (ctx.content_rect().height() * 0.6).clamp(200.0, 520.0);
        let modal = egui::Modal::new(egui::Id::unique("command-palette")).show(ctx, |ui| {
            ui.set_width(460.0);
            let response =
                ui.add(egui::TextEdit::singleline(&mut self.palette.query).hint_text("type a command").desired_width(f32::INFINITY).font(style::ui_font(15.0)));
            response.request_focus();
            if response.changed() {
                self.palette.selected = 0;
                self.palette.reveal = true;
            }
            ui.input(|i| {
                if i.key_pressed(Key::ArrowDown) {
                    self.palette.selected = (self.palette.selected + 1).min(count.saturating_sub(1));
                    self.palette.reveal = true;
                }
                if i.key_pressed(Key::ArrowUp) {
                    self.palette.selected = self.palette.selected.saturating_sub(1);
                    self.palette.reveal = true;
                }
                if i.key_pressed(Key::Escape) {
                    close = true;
                }
            });
            if ui.input(|i| i.key_pressed(Key::Enter)) {
                chosen = scored.get(self.palette.selected).map(|s| s.1);
            }
            ui.add_space(4.0);
            // The selection is brought into view once per request, and otherwise the list stays
            // wherever the wheel or the scroll bar leaves it. The first command goes back to the
            // very top of the list, so its group heading shows too.
            let reveal = std::mem::take(&mut self.palette.reveal);
            let mut area = egui::ScrollArea::vertical().max_height(list_height);
            if reveal && self.palette.selected == 0 {
                area = area.vertical_scroll_offset(0.0);
            }
            let output = area.show(ui, |ui| {
                let mut heading = None;
                for (ii, (_, action)) in scored.iter().enumerate() {
                    if grouped && heading != Some(action.group()) {
                        heading = Some(action.group());
                        if ii > 0 {
                            ui.add_space(6.0);
                        }
                        ui.label(RichText::new(action.group().to_uppercase()).font(style::bold_font(10.5)).color(palette.muted));
                    }
                    let selected = ii == self.palette.selected;
                    let shortcut = action.shortcut().map(|s| ctx.format_shortcut(&s)).unwrap_or_default();
                    let response = ui
                        .horizontal(|ui| {
                            let r = ui.selectable_label(selected, action.label());
                            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                ui.label(RichText::new(shortcut).small().color(palette.faint));
                            });
                            r
                        })
                        .inner;
                    if selected && reveal && ii > 0 {
                        response.scroll_to_me(None);
                    }
                    if response.clicked() {
                        chosen = Some(*action);
                    }
                }
            });
            self.palette.scrolled = output.state.clamped_offset().y;
        });
        if modal.should_close() || close {
            self.palette.open = false;
        }
        if let Some(action) = chosen {
            self.palette.open = false;
            actions.push(action);
        }
    }

    fn dialogs(&mut self, ctx: &egui::Context, palette: &Palette) {
        // The update dialog, which closes the app once a new version is on its way.
        let unsaved: Vec<String> =
            self.tabs.iter().filter(|t| t.dirty).map(|t| t.title()).chain(self.protocols.iter().filter(|t| t.dirty).map(|t| t.title())).collect();
        self.updater.show(ctx, palette, &unsaved);
        if self.updater.quit && !self.allow_quit {
            self.allow_quit = true;
            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
        }
        // Go to a position or select a range.
        if let Some(text) = self.goto.as_mut() {
            let mut done = false;
            let mut target = None;
            let modal = egui::Modal::new(egui::Id::unique("goto")).show(ctx, |ui| {
                ui.set_width(300.0);
                ui.label("Go to a position, or select a range such as 120..450");
                let r = ui.add(egui::TextEdit::singleline(text).desired_width(f32::INFINITY));
                r.request_focus();
                if ui.input(|i| i.key_pressed(Key::Enter)) {
                    target = Some(text.clone());
                }
                if ui.input(|i| i.key_pressed(Key::Escape)) {
                    done = true;
                }
            });
            if let Some(t) = target
                && let Some(tab) = self.tabs.get_mut(self.active)
            {
                match slopgene_core::span::parse_range(&t, tab.len(), tab.doc.is_circular()) {
                    Ok(span) if span.len() == 1 => {
                        tab.select(None);
                        tab.cursor = span.start;
                        tab.scroll_to = Some(span.start);
                        done = true;
                    }
                    Ok(span) => {
                        tab.select(Some(span));
                        tab.scroll_to = Some(span.start);
                        done = true;
                    }
                    Err(e) => self.notice = Some((e.to_string(), Instant::now())),
                }
            }
            if done || modal.should_close() {
                self.goto = None;
            }
        }

        // New sequence from pasted bases.
        if let Some(dialog) = self.new_sequence.as_mut() {
            let mut create = false;
            let mut cancel = false;
            let modal = egui::Modal::new(egui::Id::unique("new-sequence")).show(ctx, |ui| {
                ui.set_width(460.0);
                ui.heading("New sequence");
                ui.add(egui::TextEdit::singleline(&mut dialog.name).hint_text("name"));
                ui.checkbox(&mut dialog.circular, "circular");
                ui.add(
                    egui::TextEdit::multiline(&mut dialog.sequence)
                        .hint_text("paste bases (spaces and numbers are ignored)")
                        .font(style::mono_font(12.5))
                        .desired_rows(8)
                        .desired_width(f32::INFINITY),
                );
                if let Some(e) = &dialog.error {
                    ui.label(RichText::new(e).color(palette.warn));
                }
                ui.horizontal(|ui| {
                    if ui.button("Create").clicked() {
                        create = true;
                    }
                    if ui.button("Cancel").clicked() {
                        cancel = true;
                    }
                });
            });
            if create {
                match alphabet::clean_nucleotides(&dialog.sequence) {
                    Ok(bases) if !bases.is_empty() => {
                        let topology = if dialog.circular { Topology::Circular } else { Topology::Linear };
                        let mut doc = Document::new_dna(dialog.name.trim(), bases.to_ascii_uppercase(), topology);
                        doc.meta.synthetic = Some(true);
                        doc.meta.created = Some(clock::today_iso());
                        let mut tab = DocTab::new(doc, None);
                        tab.dirty = true;
                        self.tabs.push(tab);
                        self.active = self.tabs.len() - 1;
                        self.new_sequence = None;
                    }
                    Ok(_) => dialog.error = Some("paste some bases first".into()),
                    Err(e) => dialog.error = Some(e),
                }
            } else if cancel || modal.should_close() {
                self.new_sequence = None;
            }
        }

        // Add a feature over the selection.
        if let Some(dialog) = self.feature_dialog.as_mut() {
            let mut add = false;
            let mut cancel = false;
            let modal = egui::Modal::new(egui::Id::unique("add-feature")).show(ctx, |ui| {
                ui.set_width(340.0);
                ui.heading("Add feature");
                let r = ui.add(egui::TextEdit::singleline(&mut dialog.name).desired_width(f32::INFINITY));
                if r.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                    add = true;
                }
                ui.horizontal(|ui| {
                    egui::ComboBox::from_id_salt("new-feature-type").selected_text(&dialog.kind).show_ui(ui, |ui| {
                        for t in panels::FEATURE_TYPES {
                            ui.selectable_value(&mut dialog.kind, t.to_string(), *t);
                        }
                    });
                    egui::ComboBox::from_id_salt("new-feature-strand").selected_text(dialog.strand.name()).show_ui(ui, |ui| {
                        for s in [Strand::Forward, Strand::Reverse, Strand::None, Strand::Both] {
                            ui.selectable_value(&mut dialog.strand, s, s.name());
                        }
                    });
                });
                ui.horizontal(|ui| {
                    if ui.button("Add").clicked() {
                        add = true;
                    }
                    if ui.button("Cancel").clicked() {
                        cancel = true;
                    }
                });
            });
            if add {
                let (name, kind, strand) = (dialog.name.clone(), dialog.kind.clone(), dialog.strand);
                self.feature_dialog = None;
                if let Some(tab) = self.tabs.get_mut(self.active)
                    && let Some(sel) = tab.selection
                {
                    let mut feature = Feature::new(name, kind.clone(), sel, strand);
                    if kind == "CDS" {
                        feature.segments[0].translated = Some(true);
                    }
                    let uid = feature.uid;
                    let _ = tab.edit(|d| {
                        d.features.push(feature);
                        Ok(())
                    });
                    tab.select_feature(uid);
                }
            } else if cancel || modal.should_close() {
                self.feature_dialog = None;
            }
        }

        // Restriction digest with a virtual gel.
        let mut run = false;
        let mut open_fragment = None;
        if let Some(dialog) = self.digest_dialog.as_mut() {
            let mut open = true;
            egui::Window::new("Restriction digest").open(&mut open).default_width(560.0).collapsible(false).show(ctx, |ui| {
                ui.horizontal(|ui| {
                    let r = ui.add(egui::TextEdit::singleline(&mut dialog.enzymes).hint_text("enzymes, e.g. EcoRI, BamHI").desired_width(300.0));
                    if (r.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter))) || ui.button("Digest").clicked() {
                        run = true;
                    }
                });
                match &dialog.result {
                    Some(Ok(fragments)) => {
                        ui.add_space(6.0);
                        ui.horizontal_top(|ui| {
                            ui.vertical(|ui| {
                                for (ii, (len, ends, _)) in fragments.iter().enumerate() {
                                    ui.horizontal(|ui| {
                                        ui.label(RichText::new(format!("{:>2}", ii + 1)).color(palette.faint));
                                        ui.label(RichText::new(format!("{} bp", thousands(*len))).font(style::bold_font(13.0)));
                                        ui.label(RichText::new(ends).small().color(palette.muted));
                                        if ui.small_button("open").clicked() {
                                            open_fragment = Some(ii);
                                        }
                                    });
                                }
                            });
                            if !dialog.lanes.is_empty() {
                                let scene = render_gel(&dialog.lanes, &GelOptions { height: 300.0, lane_width: 46.0, ..GelOptions::default() });
                                let (rect, _) = ui.allocate_exact_size(vec2(scene.width, scene.height), Sense::hover());
                                scene_paint::paint(ui.painter(), &scene, rect.min);
                            }
                        });
                    }
                    Some(Err(e)) => {
                        ui.label(RichText::new(e).color(palette.warn));
                    }
                    None => {}
                }
            });
            if !open {
                self.digest_dialog = None;
            }
        }
        if run {
            self.run_digest();
        }
        if let Some(ii) = open_fragment
            && let Some(Some(Ok(fragments))) = self.digest_dialog.as_ref().map(|d| d.result.as_ref())
            && let Some((_, _, doc)) = fragments.get(ii)
        {
            let mut tab = DocTab::new(doc.clone(), None);
            tab.dirty = true;
            self.tabs.push(tab);
            self.active = self.tabs.len() - 1;
        }

        // Unsaved changes.
        if let Some(target) = self.confirm_close {
            let documents = self.tabs.iter().filter(|t| t.dirty).map(|t| t.doc.name.clone());
            let protocols = self.protocols.iter().filter(|t| t.dirty).map(|t| t.name());
            let names: Vec<String> = match target {
                CloseTarget::App => documents.chain(protocols).collect(),
                CloseTarget::Document(ii) => vec![self.tabs[ii].doc.name.clone()],
                CloseTarget::Protocol(ii) => vec![self.protocols[ii].name()],
            };
            let mut choice = None;
            egui::Modal::new(egui::Id::unique("confirm-close")).show(ctx, |ui| {
                ui.set_width(360.0);
                ui.heading("Save changes?");
                ui.label(format!("{} {} unsaved changes.", names.join(", "), if names.len() == 1 { "has" } else { "have" }));
                ui.horizontal(|ui| {
                    if ui.button("Save").clicked() {
                        choice = Some(0);
                    }
                    if ui.button("Discard").clicked() {
                        choice = Some(1);
                    }
                    if ui.button("Cancel").clicked() {
                        choice = Some(2);
                    }
                });
            });
            match choice {
                Some(0) | Some(1) => {
                    let save = choice == Some(0);
                    // The documents and the protocols to save first.
                    let (documents, protocols): (Vec<usize>, Vec<usize>) = match target {
                        CloseTarget::App => (
                            (0..self.tabs.len()).filter(|&ii| self.tabs[ii].dirty).collect(),
                            (0..self.protocols.len()).filter(|&ii| self.protocols[ii].dirty).collect(),
                        ),
                        CloseTarget::Document(ii) => (vec![ii], Vec::new()),
                        CloseTarget::Protocol(ii) => (Vec::new(), vec![ii]),
                    };
                    let mut failed = false;
                    if save {
                        // Each goes to its own file, or to a file the person chooses when it has none.
                        for ii in documents {
                            if self.tabs[ii].path.is_some() {
                                failed |= self.tabs[ii].save().is_err();
                            } else {
                                self.active = ii;
                                self.save_as_dialog(None);
                                failed |= self.tabs[ii].dirty;
                            }
                        }
                        for ii in protocols {
                            if self.protocols[ii].path.is_some() {
                                failed |= self.protocols[ii].save().is_err();
                            } else {
                                self.save_protocol_as(ii);
                                failed |= self.protocols[ii].dirty;
                            }
                        }
                    }
                    self.confirm_close = None;
                    if !failed {
                        match target {
                            CloseTarget::App => {
                                self.allow_quit = true;
                                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                            }
                            CloseTarget::Document(ii) => self.close_document(ii),
                            CloseTarget::Protocol(ii) => self.close_protocol(ii),
                        }
                    }
                }
                Some(_) => self.confirm_close = None,
                None => {}
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use slopgene_core::protocol::{OpenParams, Status};

    /// A key press with Ctrl held (Cmd on macOS), and Shift if asked, as the window system
    /// reports it.
    fn ctrl(key: Key, shift: bool) -> egui::Event {
        let modifiers = Modifiers { ctrl: true, command: true, shift, ..Modifiers::default() };
        egui::Event::Key { key, physical_key: None, pressed: true, repeat: false, modifiers }
    }

    /// Run one frame with `events`. Returns what `read_shortcuts` made of them and the events it
    /// left for the widgets.
    fn frame(events: Vec<egui::Event>, document_keys: bool) -> (KeyRequests, Vec<egui::Event>) {
        let ctx = egui::Context::default();
        let mut out = None;
        let mut output = ctx.run_ui(egui::RawInput { events, ..Default::default() }, |ui| {
            let requests = read_shortcuts(ui.ctx(), document_keys);
            let left = ui.input(|i| i.events.clone());
            out = Some((requests, left));
        });
        // Without a renderer nothing uploads the font textures; egui insists they are handled.
        output.textures_delta.clear();
        out.expect("the frame ran")
    }

    #[test]
    fn redo_is_not_mistaken_for_undo() {
        assert_eq!(frame(vec![ctrl(Key::Z, true)], true).0.actions, vec![Action::Redo]);
        assert_eq!(frame(vec![ctrl(Key::Z, false)], true).0.actions, vec![Action::Undo]);
        assert_eq!(frame(vec![ctrl(Key::Y, false)], true).0.actions, vec![Action::Redo]);
        assert_eq!(frame(vec![ctrl(Key::S, true)], true).0.actions, vec![Action::SaveAs]);
        assert_eq!(frame(vec![ctrl(Key::A, false)], true).0.actions, vec![Action::SelectAll]);
        assert_eq!(frame(vec![ctrl(Key::A, true)], true).0.actions, vec![Action::AlignRead]);
    }

    #[test]
    fn specific_shortcuts_are_tried_first() {
        let extra = |s: &KeyboardShortcut| s.modifiers.shift as u8 + s.modifiers.alt as u8;
        let table = shortcut_table();
        for (ii, (_, earlier)) in table.iter().enumerate() {
            for (_, later) in &table[ii + 1..] {
                let shadowed = earlier.logical_key == later.logical_key && extra(later) > extra(earlier);
                assert!(!shadowed, "{earlier:?} would catch {later:?}");
            }
        }
    }

    #[test]
    fn clipboard_events_act_on_the_document() {
        let (requests, left) = frame(vec![egui::Event::Copy, egui::Event::Cut, egui::Event::Paste("ACGT".into())], true);
        assert_eq!(requests.actions, vec![Action::CopySelection, Action::Cut]);
        assert_eq!(requests.pastes, vec!["ACGT".to_string()]);
        assert!(left.is_empty(), "the document took them: {left:?}");
    }

    #[test]
    fn a_text_field_keeps_its_editing_keys() {
        let events = vec![ctrl(Key::Z, false), ctrl(Key::A, false), egui::Event::Copy, egui::Event::Paste("x".into())];
        let (requests, left) = frame(events, false);
        assert_eq!(requests, KeyRequests::default());
        assert_eq!(left.len(), 4, "the text field still sees every event");
        // Shortcuts that are not editing keys still work, such as the command palette.
        assert_eq!(frame(vec![ctrl(Key::K, false)], false).0.actions, vec![Action::Palette]);
    }

    #[test]
    fn the_sequence_view_having_the_keyboard_leaves_undo_working() {
        // The sequence view takes keyboard focus when clicked. Undo used to be skipped whenever
        // any widget had focus. Only a text field should keep the editing keys for itself.
        let ctx = egui::Context::default();
        let view = egui::Id::unique("view");
        let field = egui::Id::unique("field");
        let mut text = String::new();
        // seen[pass] = (some widget has the keyboard, the editing keys reach the document).
        let mut seen = [None; 3];
        for (pass, slot) in seen.iter_mut().enumerate() {
            let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
                let area = egui::Rect::from_min_size(egui::pos2(0.0, 0.0), vec2(50.0, 50.0));
                let response = ui.interact(area, view, Sense::click());
                let edit = ui.add(egui::TextEdit::singleline(&mut text).id(field));
                *slot = Some((ui.ctx().egui_wants_keyboard_input(), document_has_keys(ui.ctx(), false)));
                match pass {
                    0 => response.request_focus(),
                    1 => edit.request_focus(),
                    _ => {}
                }
            });
            output.textures_delta.clear();
        }
        // Pass 1: the view has the keyboard, and the editing keys still reach the document.
        assert_eq!(seen[1], Some((true, true)));
        // Pass 2: a text field has it, and they go to the text field instead.
        assert_eq!(seen[2], Some((true, false)));
        // A dialog keeps them away from the document behind it.
        assert!(!document_has_keys(&ctx, true));
    }

    /// The real app run headless, one frame at a time, fed with pointer and key events.
    struct Harness {
        ctx: egui::Context,
        app: SlopGeneApp,
        frame: eframe::Frame,
        time: f64,
        /// Text the app put on the clipboard.
        copied: Vec<String>,
    }

    impl Harness {
        /// The app with one 12 bp circular document open in `view`.
        fn new(view: View) -> Harness {
            let ctx = egui::Context::default();
            let mut app = SlopGeneApp::new(&eframe::CreationContext::_new_kittest(ctx.clone()), StartOptions::default());
            let mut tab = DocTab::new(Document::new_dna("t", b"AAAACCCCGGGG", Topology::Circular), None);
            tab.view = view;
            app.tabs.push(tab);
            app.active = 0;
            let mut harness = Harness { ctx, app, frame: eframe::Frame::_new_kittest(), time: 0.0, copied: Vec::new() };
            // Two quiet frames, so every widget has a place on screen before the first click.
            harness.run(Vec::new());
            harness.run(Vec::new());
            harness
        }

        fn run(&mut self, events: Vec<egui::Event>) {
            use eframe::App as _;
            self.time += 1.0 / 60.0;
            let screen = egui::Rect::from_min_size(egui::Pos2::ZERO, vec2(1200.0, 800.0));
            let raw = egui::RawInput { events, time: Some(self.time), screen_rect: Some(screen), ..Default::default() };
            let (app, frame) = (&mut self.app, &mut self.frame);
            let mut output = self.ctx.run_ui(raw, |ui| app.ui(ui, frame));
            output.textures_delta.clear();
            for command in output.platform_output.commands {
                if let egui::OutputCommand::CopyText(text) = command {
                    self.copied.push(text);
                }
            }
        }

        fn click(&mut self, pos: egui::Pos2) {
            let button = |pressed| egui::Event::PointerButton { pos, button: egui::PointerButton::Primary, pressed, modifiers: Modifiers::default() };
            self.run(vec![egui::Event::PointerMoved(pos)]);
            self.run(vec![button(true)]);
            self.run(vec![button(false)]);
            self.run(Vec::new());
        }

        fn press(&mut self, key: Key, modifiers: Modifiers) {
            let event = |pressed| egui::Event::Key { key, physical_key: None, pressed, repeat: false, modifiers };
            self.run(vec![event(true)]);
            self.run(vec![event(false)]);
        }

        fn ctrl(&mut self, key: Key) {
            self.press(key, Modifiers { ctrl: true, command: true, ..Modifiers::default() });
        }

        fn ctrl_shift(&mut self, key: Key) {
            self.press(key, Modifiers { ctrl: true, command: true, shift: true, ..Modifiers::default() });
        }

        fn seq(&self) -> String {
            String::from_utf8_lossy(&self.app.tabs[0].doc.seq).into_owned()
        }

        fn select(&mut self, start: usize, end: usize) {
            self.app.tabs[0].select(Some(Span::new(start, end)));
        }
    }

    #[test]
    fn delete_in_the_sequence_view_then_ctrl_z_brings_the_bases_back() {
        let mut h = Harness::new(View::Sequence);
        h.click(egui::pos2(200.0, 140.0));
        h.select(4, 8);
        h.press(Key::Delete, Modifiers::default());
        assert_eq!(h.seq(), "AAAAGGGG");
        h.ctrl(Key::Z);
        assert_eq!(h.seq(), "AAAACCCCGGGG", "Ctrl+Z undoes the deletion");
        assert_eq!(h.app.tabs[0].selection, Some(Span::new(4, 8)), "and selects the bases again");
        h.ctrl_shift(Key::Z);
        assert_eq!(h.seq(), "AAAAGGGG", "Ctrl+Shift+Z redoes it");
        h.ctrl(Key::Z);
        h.ctrl(Key::Y);
        assert_eq!(h.seq(), "AAAAGGGG", "Ctrl+Y redoes too");
    }

    #[test]
    fn clipboard_keys_work_in_the_sequence_view() {
        let mut h = Harness::new(View::Sequence);
        h.click(egui::pos2(200.0, 140.0));
        h.select(4, 8);
        h.run(vec![egui::Event::Copy]);
        assert_eq!(h.copied, vec!["CCCC".to_string()]);
        h.run(vec![egui::Event::Cut]);
        assert_eq!(h.seq(), "AAAAGGGG");
        assert_eq!(h.copied.last().map(String::as_str), Some("CCCC"));
        // The cursor sits where the bases were, so pasting puts them straight back.
        h.run(vec![egui::Event::Paste("cccc".into())]);
        assert_eq!(h.seq(), "AAAACCCCGGGG");
        h.ctrl(Key::A);
        assert_eq!(h.app.tabs[0].selection, Some(Span::new(0, 12)));
    }

    #[test]
    fn the_map_takes_delete_copy_and_undo_too() {
        let mut h = Harness::new(View::Map);
        h.click(egui::pos2(450.0, 420.0));
        h.select(0, 4);
        h.run(vec![egui::Event::Copy]);
        assert_eq!(h.copied, vec!["AAAA".to_string()]);
        h.press(Key::Backspace, Modifiers::default());
        assert_eq!(h.seq(), "CCCCGGGG");
        h.ctrl(Key::Z);
        assert_eq!(h.seq(), "AAAACCCCGGGG");
    }

    /// Frames enough for a scroll animation to finish.
    fn settle(h: &mut Harness) {
        for _ in 0..120 {
            h.run(Vec::new());
        }
    }

    #[test]
    fn the_command_list_stays_where_the_wheel_leaves_it() {
        let mut h = Harness::new(View::Sequence);
        h.ctrl(Key::K);
        assert!(h.app.palette.open);
        // Turn the wheel over the list (the palette sits in the middle of the window).
        let over_list = egui::pos2(600.0, 420.0);
        h.run(vec![egui::Event::PointerMoved(over_list)]);
        let wheel = egui::Event::MouseWheel {
            unit: egui::MouseWheelUnit::Line,
            delta: vec2(0.0, -6.0),
            phase: egui::TouchPhase::Move,
            source: egui::MouseWheelSource::Unknown,
            modifiers: Modifiers::default(),
        };
        h.run(vec![wheel]);
        settle(&mut h);
        let scrolled = h.app.palette.scrolled;
        assert!(scrolled > 50.0, "the wheel scrolled the list, by {scrolled}");
        // The list used to snap back to the selected first command on the next frames.
        settle(&mut h);
        assert_eq!(h.app.palette.scrolled, scrolled, "and it stays there");
    }

    #[test]
    fn the_arrow_keys_still_bring_the_selection_into_view() {
        let mut h = Harness::new(View::Sequence);
        h.ctrl(Key::K);
        for _ in 0..Action::ALL.len() {
            h.press(Key::ArrowDown, Modifiers::default());
        }
        settle(&mut h);
        assert_eq!(h.app.palette.selected, Action::ALL.len() - 1);
        assert!(h.app.palette.scrolled > 100.0, "the last command scrolled into view");
        // A query starts again at the top of the list.
        h.run(vec![egui::Event::Text("align".into())]);
        settle(&mut h);
        assert_eq!((h.app.palette.selected, h.app.palette.scrolled), (0, 0.0));
    }

    #[test]
    fn an_open_document_aligns_to_another() {
        let mut h = Harness::new(View::Sequence);
        // A 2 kb circular reference in the first tab, and in a second tab 600 of its bases with
        // reference base 401 (0-based 400) changed.
        let mut x = 3u64;
        let reference: Vec<u8> = (0..2000)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect();
        let mut read = reference[300..900].to_vec();
        read[100] = if read[100] == b'A' { b'C' } else { b'A' };
        h.app.tabs[0] = DocTab::new(Document::new_dna("ref", reference, Topology::Circular), None);
        h.app.tabs.push(DocTab::new(Document::new_dna("read", read, Topology::Linear), None));
        h.run(Vec::new());
        // Another DNA document is open, so Ctrl+Shift+A shows the alignment view without asking
        // for files.
        h.ctrl_shift(Key::A);
        assert_eq!(h.app.tabs[0].view, View::Alignment);
        let read_id = h.app.tabs[1].id;
        h.app.align_open_document(read_id);
        h.run(Vec::new());
        let aln = h.app.tabs[0].alignments.entries[0].result.as_ref().expect("the read aligns");
        let positions: Vec<usize> = aln.differences.iter().map(|d| d.position).collect();
        assert_eq!(positions, vec![400]);
        // The map draws the read as a ring with its one difference.
        let tracks = h.app.tabs[0].alignments.tracks();
        assert_eq!(tracks.len(), 1);
        assert_eq!(tracks[0].marks.len(), 1);
        // A position asked for from the map is brought into view among the aligned bases.
        h.app.tabs[0].scroll_to = Some(400);
        h.run(Vec::new());
        assert_eq!(h.app.tabs[0].scroll_to, None, "the aligned bases took the request");
        // Another view leaves the alignment in place for when the view comes back.
        h.ctrl(Key::Num1);
        assert_eq!(h.app.tabs[0].view, View::Map);
        assert_eq!(h.app.tabs[0].alignments.entries.len(), 1);
    }

    #[test]
    fn each_group_of_commands_is_listed_together() {
        // The palette writes a heading where the group changes, so a group split in two would get
        // two headings.
        let mut seen: Vec<&str> = Vec::new();
        for action in Action::ALL {
            if seen.last() != Some(&action.group()) {
                assert!(!seen.contains(&action.group()), "{} is split", action.group());
                seen.push(action.group());
            }
        }
        assert_eq!(seen, ["File", "Tools", "Edit", "View", "Export"]);
    }

    #[test]
    fn ctrl_z_in_a_text_field_leaves_the_document_alone() {
        let mut h = Harness::new(View::Sequence);
        h.click(egui::pos2(200.0, 140.0));
        h.select(4, 8);
        h.press(Key::Delete, Modifiers::default());
        // Ctrl+F opens the find bar and puts the keyboard in its text field.
        h.ctrl(Key::F);
        h.run(Vec::new());
        h.ctrl(Key::Z);
        assert_eq!(h.seq(), "AAAAGGGG", "the find field had the keyboard, so the document is untouched");
    }

    #[test]
    fn a_protocol_tab_takes_the_file_and_history_commands() {
        let mut h = Harness::new(View::Map);
        // A protocol of one step, which opens a map stored in the protocol itself.
        let stored = Document::new_dna("stored", b"ACGTACGTAC", Topology::Circular);
        let open = OpenParams { path: None, document: Some(json::document_to_value(&stored, true)) };
        let step =
            slopgene_core::protocol::Step { id: "vector".into(), inputs: Vec::new(), action: slopgene_core::protocol::Action::Open(open), note: String::new() };
        h.app.add_protocol(ProtocolTab::new(Protocol { name: "test".into(), steps: vec![step], ..Protocol::default() }, None));
        h.run(Vec::new());
        assert_eq!(h.app.front, Front::Protocol(0));
        assert_eq!(h.app.protocols[0].run.steps[0].status, Status::Done, "the protocol in front runs");
        // Ctrl+Z takes back a change to the protocol and leaves the document behind it alone.
        h.app.protocols[0].add_step("digest");
        h.ctrl(Key::Z);
        assert_eq!(h.app.protocols[0].protocol.steps.len(), 1);
        assert_eq!(h.seq(), "AAAACCCCGGGG");
        // Closing a changed protocol asks first. Closing it brings the document back to the front.
        h.ctrl(Key::W);
        assert_eq!(h.app.confirm_close, Some(CloseTarget::Protocol(0)));
        h.app.confirm_close = None;
        h.app.close_protocol(0);
        assert_eq!(h.app.front, Front::Document);
        assert_eq!(h.app.tabs.len(), 1);
        // Ctrl+Shift+N, like the Clone button, starts a cloning protocol in a tab of its own.
        h.ctrl_shift(Key::N);
        assert_eq!((h.app.protocols.len(), h.app.front), (1, Front::Protocol(0)));
    }

    #[test]
    fn a_marked_region_makes_the_map_in_front_the_inserts() {
        let mut h = Harness::new(View::Map);
        h.app.tabs.push(DocTab::new(Document::new_dna("other", b"ACGTACGTACGT", Topology::Circular), None));
        h.run(Vec::new());
        // The name of the map a protocol's open step reads.
        let opened = |h: &Harness, id: &str| h.app.protocols.last().and_then(|p| p.run.step(id)).and_then(|r| r.outputs.first()).map(|p| p.doc.name.clone());
        // Nothing marked: the map in front is the vector, and the other map carries the insert.
        h.ctrl_shift(Key::N);
        assert_eq!((opened(&h, "vector").as_deref(), opened(&h, "donor").as_deref()), (Some("t"), Some("other")));
        h.app.close_protocol(0);
        // A region marked on the map in front: that map carries the insert, and the region is the
        // part its cut is about. The protocol opens on that cut, where the enzymes are chosen.
        h.select(2, 8);
        h.ctrl_shift(Key::N);
        assert_eq!((opened(&h, "donor").as_deref(), opened(&h, "vector").as_deref()), (Some("t"), Some("other")));
        let tab = h.app.protocols.last().unwrap();
        let Some(slopgene_core::protocol::Action::Digest(cut)) = tab.protocol.step("cut_donor").map(|s| &s.action) else { panic!("the donor is cut") };
        assert_eq!(cut.around.as_deref(), Some("3..8"));
        assert_eq!(tab.protocol.steps[tab.selected].id, "cut_donor");
        // The bench is a tab, not a window over the maps, and it draws frame after frame.
        h.run(Vec::new());
        assert!(!h.app.modal_open());
        assert_eq!(h.app.front, Front::Protocol(0));
    }

    #[test]
    fn a_protocol_file_is_told_from_a_document_file() {
        let folder = std::env::temp_dir().join(format!("slopgene-protocol-test-{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let protocol = folder.join("clone.json");
        let document = folder.join("map.json");
        std::fs::write(&protocol, Protocol { name: "p".into(), ..Protocol::default() }.to_json()).unwrap();
        let map = Document::new_dna("map", b"ACGT", Topology::Linear);
        std::fs::write(&document, json::document_to_value(&map, true).to_string()).unwrap();
        assert!(is_protocol_file(&protocol));
        assert!(!is_protocol_file(&document), "a document in JSON opens as a sequence");
        let _ = std::fs::remove_dir_all(&folder);
    }
}
