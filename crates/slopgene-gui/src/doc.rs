//! One open document: the sequence, where it lives on disk, the selection, undo history and
//! cached analyses.
//!
//! The file on disk is polled about once a second. When another program (a script, an agent
//! using the `slopgene` command line) changes it and there are no unsaved edits here, the
//! document reloads by itself, so edits made outside the app appear live.

use slopgene_core::enzymes::{EnzymeDb, find_sites};
use slopgene_core::primers::{self, BindingRules};
use slopgene_core::render::SiteMark;
use slopgene_core::*;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant, SystemTime};

/// How many undo steps to keep per document.
const UNDO_LIMIT: usize = 200;

/// Edits of one kind that follow each other within this time (a burst of typing, a run of
/// Backspace) share one undo step.
const MERGE_WINDOW: Duration = Duration::from_secs(2);

/// Where the views stood when a snapshot was taken. Undo and redo restore it with the
/// document, so undoing a deletion brings the bases back selected.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Marks {
    cursor: usize,
    selection: Option<Span>,
    anchor: Option<usize>,
    selected_feature: Option<u64>,
}

/// How [`DocTab::save_keeping_everything`] went.
#[derive(Debug, PartialEq)]
pub enum Saved {
    /// Written to its own file.
    InPlace(PathBuf),
    /// The file's format would have lost something, so a SnapGene file beside it was written
    /// and the document now lives there. The original file is unchanged.
    Beside { path: PathBuf, original: PathBuf, format: io::Format, lost: Vec<&'static str> },
    /// The file's format would lose something and a SnapGene file of that name already exists:
    /// ask where to save.
    NeedsLocation { format: io::Format, lost: Vec<&'static str> },
    /// The document has no file yet.
    NoFile,
}

/// A format's name as people write it.
pub fn format_label(format: io::Format) -> &'static str {
    match format {
        io::Format::SnapGene => "SnapGene",
        io::Format::GenBank => "GenBank",
        io::Format::Fasta => "FASTA",
        io::Format::Json => "JSON",
    }
}

/// One undo or redo step: the document as it was, and where the views stood.
struct Snapshot {
    doc: Document,
    marks: Marks,
}

/// Which main view a document shows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum View {
    Map,
    Sequence,
    Split,
    /// Sequences aligned to the document: the map and the aligned bases side by side.
    Alignment,
}

/// An enzyme and how often it cuts, for the enzymes panel.
#[derive(Clone, Debug)]
pub struct EnzymeCount {
    pub name: String,
    pub site: String,
    pub cuts: Vec<usize>,
}

/// Enzyme sites as last computed: the document revision, enzyme set and unique-only choice they
/// were computed for, then the map marks and the per-enzyme counts.
type SitesCache = (u64, String, bool, Vec<SiteMark>, Vec<EnzymeCount>);

pub struct DocTab {
    pub id: u64,
    pub doc: Document,
    pub path: Option<PathBuf>,
    pub dirty: bool,
    /// Bumped on every change, so caches know when to rebuild.
    pub revision: u64,
    undo: Vec<Snapshot>,
    redo: Vec<Snapshot>,
    /// Selected bases, if any.
    pub selection: Option<Span>,
    /// Insertion point, 0..=length (new bases go before base `cursor`).
    pub cursor: usize,
    /// Where a drag-selection started.
    pub anchor: Option<usize>,
    pub selected_feature: Option<u64>,
    pub view: View,
    /// Enzyme set shown on maps (a set name or a comma-separated list).
    pub enzyme_set: String,
    /// Show only enzymes that cut once.
    pub unique_only: bool,
    sites_cache: Option<SitesCache>,
    disk_mtime: Option<SystemTime>,
    last_disk_check: Instant,
    /// The file changed on disk while there were unsaved edits here.
    pub changed_on_disk: bool,
    /// Ask the sequence view to scroll so this base is visible.
    pub scroll_to: Option<usize>,
    /// A short note for the status bar (for example "reloaded from disk").
    pub notice: Option<(String, Instant)>,
    /// Key and time of the last coalescable edit.
    last_edit: Option<(String, Instant)>,
    /// Sequencing results and other sequences aligned to this document.
    pub alignments: crate::align_view::Alignments,
    /// The bases cannot be typed or deleted here, as in the sequence shown under a protocol's
    /// map. Selecting, copying and moving the cursor still work.
    pub read_only: bool,
}

fn mtime(path: &Path) -> Option<SystemTime> {
    std::fs::metadata(path).and_then(|m| m.modified()).ok()
}

impl DocTab {
    pub fn new(doc: Document, path: Option<PathBuf>) -> DocTab {
        let disk_mtime = path.as_deref().and_then(mtime);
        let view = if doc.is_circular() { View::Map } else { View::Split };
        DocTab {
            id: new_uid(),
            doc,
            path,
            dirty: false,
            revision: 0,
            undo: Vec::new(),
            redo: Vec::new(),
            selection: None,
            cursor: 0,
            anchor: None,
            selected_feature: None,
            view,
            enzyme_set: "common".into(),
            unique_only: true,
            sites_cache: None,
            disk_mtime,
            last_disk_check: Instant::now(),
            changed_on_disk: false,
            scroll_to: None,
            notice: None,
            last_edit: None,
            alignments: crate::align_view::Alignments::default(),
            read_only: false,
        }
    }

    pub fn open(path: &Path) -> Result<DocTab> {
        let doc = io::read_one(path)?;
        Ok(DocTab::new(doc, Some(path.to_path_buf())))
    }

    pub fn title(&self) -> String {
        let name = if self.doc.name.is_empty() { "untitled".to_string() } else { self.doc.name.clone() };
        if self.dirty { format!("{name} •") } else { name }
    }

    pub fn len(&self) -> usize {
        self.doc.len()
    }

    fn marks(&self) -> Marks {
        Marks { cursor: self.cursor, selection: self.selection, anchor: self.anchor, selected_feature: self.selected_feature }
    }

    fn restore_marks(&mut self, marks: Marks) {
        self.cursor = marks.cursor;
        self.selection = marks.selection;
        self.anchor = marks.anchor;
        self.selected_feature = marks.selected_feature;
    }

    /// Put a step on the undo stack. A new step makes the redo stack meaningless.
    fn push_undo(&mut self, doc: Document, marks: Marks) {
        self.undo.push(Snapshot { doc, marks });
        if self.undo.len() > UNDO_LIMIT {
            self.undo.remove(0);
        }
        self.redo.clear();
    }

    /// Change the document as one undo step. If the change fails, the document is left as it
    /// was.
    pub fn edit(&mut self, change: impl FnOnce(&mut Document) -> Result<()>) -> Result<()> {
        let before = self.doc.clone();
        let marks = self.marks();
        match change(&mut self.doc) {
            Ok(()) => {
                self.push_undo(before, marks);
                // Whatever comes next starts its own step.
                self.last_edit = None;
                self.changed();
                Ok(())
            }
            Err(e) => {
                self.doc = before;
                Err(e)
            }
        }
    }

    /// Whether the previous edit was of kind `key` and recent enough to share its undo step.
    fn continues(&self, key: &str) -> bool {
        self.last_edit.as_ref().is_some_and(|(k, at)| k == key && at.elapsed() < MERGE_WINDOW)
    }

    /// Change the document as part of the previous undo step (the second keystroke of a burst).
    fn edit_merged(&mut self, change: impl FnOnce(&mut Document) -> Result<()>) -> Result<()> {
        let before = self.doc.clone();
        match change(&mut self.doc) {
            Ok(()) => {
                self.changed();
                Ok(())
            }
            Err(e) => {
                self.doc = before;
                Err(e)
            }
        }
    }

    /// Like [`DocTab::edit`], but a burst of edits with the same `key` (dragging a colour,
    /// typing a name) within two seconds becomes a single undo step.
    pub fn edit_coalesced(&mut self, key: &str, change: impl FnOnce(&mut Document) -> Result<()>) -> Result<()> {
        if self.continues(key) {
            self.edit_merged(change)?;
        } else {
            self.edit(change)?;
        }
        self.last_edit = Some((key.to_string(), Instant::now()));
        Ok(())
    }

    /// Replace the document wholesale (after a script ran), keeping undo.
    pub fn replace_document(&mut self, doc: Document) {
        let marks = self.marks();
        let before = std::mem::replace(&mut self.doc, doc);
        self.push_undo(before, marks);
        self.last_edit = None;
        self.changed();
    }

    /// Bases a user typed or pasted, checked against this molecule: nucleotide codes, or
    /// amino acids in a protein. Spaces, digits and the header line of a pasted FASTA record are
    /// dropped. Anything else is refused, so a stray character never lands in the sequence.
    pub fn clean_input(&self, text: &str) -> std::result::Result<Vec<u8>, String> {
        let body: String = if text.trim_start().starts_with('>') {
            text.lines().filter(|line| !line.trim_start().starts_with('>')).collect::<Vec<_>>().join("\n")
        } else {
            text.to_string()
        };
        match self.doc.molecule {
            Molecule::Protein => alphabet::clean_protein(&body),
            _ => alphabet::clean_nucleotides(&body).map(|bases| bases.to_ascii_uppercase()),
        }
    }

    /// What one unit of this sequence is called in messages.
    pub fn unit(&self) -> &'static str {
        if self.doc.molecule == Molecule::Protein { "residues" } else { "bases" }
    }

    /// Put `bases` at the cursor, replacing the selection if there is one. With `merge` the
    /// change joins the previous undo step. The cursor ends up after the new bases.
    fn put_bases(&mut self, bases: &[u8], merge: bool) -> Result<()> {
        if bases.is_empty() {
            return Ok(());
        }
        let len = self.len();
        let selection = self.selection;
        // Where the first new base lands. A replacement of the same length overwrites in place;
        // other replacements through the origin, or of everything, start again at base 1.
        let at = match selection {
            Some(sel) if sel.len() == bases.len() => sel.start,
            Some(sel) if sel.len() >= len || sel.wraps(len) => 0,
            Some(sel) => sel.start,
            None => self.cursor.min(len),
        };
        let change = |d: &mut Document| {
            match selection {
                Some(sel) => edit::replace(d, sel, bases)?,
                None => edit::insert(d, at, bases)?,
            }
            primers::refresh_binding_sites(d, &BindingRules::default());
            Ok(())
        };
        if merge {
            self.edit_merged(change)?
        } else {
            self.edit(change)?
        }
        let new_len = self.len();
        let end = at + bases.len();
        self.select(None);
        self.cursor = if end > new_len { end - new_len } else { end };
        self.anchor = Some(self.cursor);
        self.scroll_to = Some(self.cursor);
        Ok(())
    }

    /// Type bases at the cursor, replacing the selection if there is one. Keystrokes that carry
    /// on where the previous one ended share one undo step, as in a text editor.
    pub fn type_bases(&mut self, bases: &[u8]) -> Result<()> {
        let merge = self.selection.is_none() && self.continues(&format!("type@{}", self.cursor));
        self.put_bases(bases, merge)?;
        self.last_edit = Some((format!("type@{}", self.cursor), Instant::now()));
        Ok(())
    }

    /// Paste bases at the cursor, replacing the selection if there is one, as one undo step.
    pub fn paste_bases(&mut self, bases: &[u8]) -> Result<()> {
        self.put_bases(bases, false)
    }

    /// Delete the selected bases as one undo step and return how many went.
    pub fn delete_selection(&mut self) -> Result<usize> {
        let Some(sel) = self.selection else { return Ok(0) };
        let len = self.len();
        let count = sel.len().min(len);
        let at = if sel.wraps(len) { 0 } else { sel.start };
        self.edit(|d| {
            edit::delete(d, sel)?;
            primers::refresh_binding_sites(d, &BindingRules::default());
            Ok(())
        })?;
        self.select(None);
        self.cursor = at.min(self.len());
        self.anchor = Some(self.cursor);
        self.scroll_to = Some(self.cursor);
        Ok(count)
    }

    /// Delete one base at `at`. `key` names the run it belongs to, so a run of Backspace or
    /// Delete presses shares one undo step.
    fn delete_base(&mut self, at: usize, key: &str) -> Result<()> {
        let change = |d: &mut Document| {
            edit::delete(d, Span::new(at, at + 1))?;
            primers::refresh_binding_sites(d, &BindingRules::default());
            Ok(())
        };
        if self.continues(key) { self.edit_merged(change) } else { self.edit(change) }
    }

    /// Backspace: delete the selection, or else the base before the cursor.
    pub fn delete_backward(&mut self) -> Result<()> {
        if self.selection.is_some() {
            return self.delete_selection().map(|_| ());
        }
        let cursor = self.cursor.min(self.len());
        if cursor == 0 {
            return Ok(());
        }
        let at = cursor - 1;
        self.delete_base(at, &format!("back@{cursor}"))?;
        self.cursor = at;
        self.anchor = Some(at);
        self.scroll_to = Some(at);
        self.last_edit = Some((format!("back@{at}"), Instant::now()));
        Ok(())
    }

    /// Delete: delete the selection, or else the base after the cursor.
    pub fn delete_forward(&mut self) -> Result<()> {
        if self.selection.is_some() {
            return self.delete_selection().map(|_| ());
        }
        let at = self.cursor;
        if at >= self.len() {
            return Ok(());
        }
        self.delete_base(at, &format!("forward@{at}"))?;
        self.anchor = Some(at);
        self.scroll_to = Some(at);
        self.last_edit = Some((format!("forward@{at}"), Instant::now()));
        Ok(())
    }

    /// Cut: delete the selection and return its bases for the clipboard. Nothing is returned
    /// when nothing is selected or the deletion is refused.
    pub fn cut(&mut self) -> Result<Option<Vec<u8>>> {
        let Some(bases) = self.selected_bases() else { return Ok(None) };
        self.delete_selection()?;
        Ok(Some(bases))
    }

    /// Select the whole sequence.
    pub fn select_all(&mut self) {
        let len = self.len();
        self.selection = (len > 0).then(|| Span::new(0, len));
        self.selected_feature = None;
        self.anchor = Some(0);
        self.cursor = len;
    }

    fn changed(&mut self) {
        self.dirty = true;
        self.revision += 1;
        self.clamp_positions();
    }

    /// Keep the cursor and selection inside the sequence after it shrank.
    fn clamp_positions(&mut self) {
        let len = self.doc.len();
        self.cursor = self.cursor.min(len);
        if let Some(sel) = self.selection
            && (sel.start >= len || sel.len() > len || (!self.doc.is_circular() && sel.end > len))
        {
            self.selection = None;
        }
        if let Some(uid) = self.selected_feature
            && self.doc.feature_index(uid).is_none()
        {
            self.selected_feature = None;
        }
    }

    /// Step back one edit, and put the cursor and selection back where they were before it.
    /// Returns false when there is nothing to undo.
    pub fn undo(&mut self) -> bool {
        let Some(previous) = self.undo.pop() else { return false };
        let marks = self.marks();
        let current = std::mem::replace(&mut self.doc, previous.doc);
        self.redo.push(Snapshot { doc: current, marks });
        self.restore_marks(previous.marks);
        self.after_history_step();
        true
    }

    /// Repeat an undone edit, with the cursor and selection it left. Returns false when there is
    /// nothing to redo.
    pub fn redo(&mut self) -> bool {
        let Some(next) = self.redo.pop() else { return false };
        let marks = self.marks();
        let current = std::mem::replace(&mut self.doc, next.doc);
        self.undo.push(Snapshot { doc: current, marks });
        self.restore_marks(next.marks);
        self.after_history_step();
        true
    }

    /// After undo or redo: the next keystroke starts a new step, and the views show the change.
    fn after_history_step(&mut self) {
        self.last_edit = None;
        self.changed();
        self.scroll_to = Some(self.cursor);
    }

    /// Save to the document's file.
    pub fn save(&mut self) -> Result<PathBuf> {
        let path = self.path.clone().ok_or_else(|| Error::Invalid("choose where to save first".into()))?;
        self.save_as(&path)?;
        Ok(path)
    }

    /// Save without losing anything. When the document's file is in a format that cannot hold
    /// all of it (features in a FASTA file, primers in a GenBank file), a SnapGene file of the
    /// same name is written beside it instead and the document moves there; the original file
    /// stays as it was.
    pub fn save_keeping_everything(&mut self) -> Result<Saved> {
        let Some(path) = self.path.clone() else { return Ok(Saved::NoFile) };
        let format = io::format_of(&path);
        let lost = io::losses(&self.doc, format);
        if lost.is_empty() {
            self.save_as(&path)?;
            return Ok(Saved::InPlace(path));
        }
        let beside = path.with_extension(io::Format::SnapGene.extension(self.doc.molecule));
        if beside.exists() {
            return Ok(Saved::NeedsLocation { format, lost });
        }
        self.save_as(&beside)?;
        Ok(Saved::Beside { path: beside, original: path, format, lost })
    }

    pub fn save_as(&mut self, path: &Path) -> Result<()> {
        self.doc.meta.modified = Some(clock::today_iso());
        io::write_path(&self.doc, path, None)?;
        // A save in another format (GenBank, FASTA) gives the document that file's name.
        if let Some(stem) = path.file_stem().and_then(|s| s.to_str()) {
            self.doc.name = stem.to_string();
        }
        self.path = Some(path.to_path_buf());
        self.dirty = false;
        self.disk_mtime = mtime(path);
        self.changed_on_disk = false;
        Ok(())
    }

    /// Re-read the file from disk, discarding unsaved edits (kept on the undo stack).
    pub fn reload(&mut self) -> Result<()> {
        let path = self.path.clone().ok_or_else(|| Error::Invalid("document has no file".into()))?;
        let doc = io::read_one(&path)?;
        let marks = self.marks();
        let before = std::mem::replace(&mut self.doc, doc);
        self.push_undo(before, marks);
        self.last_edit = None;
        self.revision += 1;
        self.dirty = false;
        self.changed_on_disk = false;
        self.disk_mtime = mtime(&path);
        self.clamp_positions();
        Ok(())
    }

    /// Poll the file's modification time (at most once a second). Reload when it changed and
    /// nothing is unsaved here; otherwise flag the conflict. Returns true when it reloaded.
    pub fn poll_disk(&mut self) -> bool {
        if self.last_disk_check.elapsed() < Duration::from_millis(900) {
            return false;
        }
        self.last_disk_check = Instant::now();
        let Some(path) = self.path.clone() else { return false };
        let now = mtime(&path);
        if now.is_none() || now == self.disk_mtime {
            return false;
        }
        if self.dirty {
            self.changed_on_disk = true;
            self.disk_mtime = now;
            return false;
        }
        match self.reload() {
            Ok(()) => {
                self.notice = Some(("reloaded: the file changed on disk".into(), Instant::now()));
                true
            }
            Err(e) => {
                // The writer may still be busy; try again on the next poll.
                self.notice = Some((format!("could not reload: {e}"), Instant::now()));
                self.disk_mtime = None;
                false
            }
        }
    }

    /// Enzyme sites for maps and the enzymes panel, recomputed when the document or the enzyme
    /// choice changes.
    pub fn sites(&mut self) -> (&[SiteMark], &[EnzymeCount]) {
        let fresh =
            matches!(&self.sites_cache, Some((rev, set, unique, _, _)) if *rev == self.revision && *set == self.enzyme_set && *unique == self.unique_only);
        if !fresh {
            let db = EnzymeDb::builtin();
            let enzymes = db.set(&self.enzyme_set).unwrap_or_default();
            let found = find_sites(&self.doc.seq, self.doc.is_circular(), &enzymes);
            let groups = slopgene_core::enzymes::cutting_sites_by_enzyme(&found);
            let mut counts = Vec::new();
            let mut marks = Vec::new();
            for e in &enzymes {
                let cuts: Vec<usize> = groups.get(e.name.as_str()).map(|g| g.iter().map(|s| s.cuts[0].pos).collect()).unwrap_or_default();
                if cuts.is_empty() {
                    continue;
                }
                if !self.unique_only || cuts.len() == 1 {
                    for &c in &cuts {
                        marks.push(SiteMark { name: e.name.clone(), pos: c });
                    }
                }
                counts.push(EnzymeCount { name: e.name.clone(), site: e.site_with_cuts(), cuts });
            }
            marks.sort_by_key(|m| m.pos);
            counts.sort_by(|a, b| a.cuts.len().cmp(&b.cuts.len()).then_with(|| a.name.cmp(&b.name)));
            self.sites_cache = Some((self.revision, self.enzyme_set.clone(), self.unique_only, marks, counts));
        }
        let (_, _, _, marks, counts) = self.sites_cache.as_ref().unwrap();
        (marks, counts)
    }

    /// Select a feature: its extent becomes the selection and the views scroll to it.
    pub fn select_feature(&mut self, uid: u64) {
        if let Some(index) = self.doc.feature_index(uid) {
            let span = self.doc.features[index].span(self.doc.len());
            self.selection = Some(span);
            self.cursor = span.start;
            self.selected_feature = Some(uid);
            self.scroll_to = Some(span.start);
        }
    }

    /// Select `span` (or clear the selection with `None`).
    pub fn select(&mut self, span: Option<Span>) {
        self.selection = span.filter(|s| !s.is_empty());
        self.selected_feature = None;
        if let Some(s) = self.selection {
            self.cursor = s.start;
        }
    }

    /// Bases of the selection, or nothing.
    pub fn selected_bases(&self) -> Option<Vec<u8>> {
        self.selection.map(|s| s.extract(&self.doc.seq))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 12 bp linear document with one feature over bases 5..8 (0-based 4..8).
    fn tab() -> DocTab {
        let mut doc = Document::new_dna("t", b"AAAACCCCGGGG", Topology::Linear);
        doc.features.push(Feature::new("f", "misc_feature", Span::new(4, 8), Strand::Forward));
        DocTab::new(doc, None)
    }

    fn seq(tab: &DocTab) -> String {
        String::from_utf8_lossy(&tab.doc.seq).into_owned()
    }

    #[test]
    fn undoing_a_deletion_brings_the_bases_back_selected() {
        let mut t = tab();
        t.select(Some(Span::new(4, 8)));
        assert_eq!(t.delete_selection().unwrap(), 4);
        assert_eq!(seq(&t), "AAAAGGGG");
        assert!(t.undo());
        assert_eq!(seq(&t), "AAAACCCCGGGG");
        assert_eq!(t.selection, Some(Span::new(4, 8)));
        assert_eq!(t.doc.features.len(), 1, "the feature comes back too");
        assert!(t.redo());
        assert_eq!(seq(&t), "AAAAGGGG");
        assert_eq!((t.selection, t.cursor), (None, 4));
        assert!(t.undo());
        assert!(!t.undo(), "only one step was made");
    }

    #[test]
    fn a_burst_of_typing_is_one_undo_step() {
        let mut t = tab();
        t.cursor = 2;
        for base in [b"T", b"T", b"T"] {
            t.type_bases(base).unwrap();
        }
        assert_eq!(seq(&t), "AATTTAACCCCGGGG");
        assert_eq!(t.cursor, 5);
        assert!(t.undo());
        assert_eq!(seq(&t), "AAAACCCCGGGG");
        assert_eq!(t.cursor, 2, "the cursor goes back to where the typing began");
        assert!(!t.undo());
    }

    #[test]
    fn typing_somewhere_else_starts_a_new_step() {
        let mut t = tab();
        t.cursor = 2;
        t.type_bases(b"T").unwrap();
        t.cursor = 0;
        t.type_bases(b"G").unwrap();
        assert_eq!(seq(&t), "GAATAACCCCGGGG");
        assert!(t.undo());
        assert_eq!(seq(&t), "AATAACCCCGGGG");
        assert!(t.undo());
        assert_eq!(seq(&t), "AAAACCCCGGGG");
    }

    #[test]
    fn a_run_of_backspaces_is_one_undo_step() {
        let mut t = tab();
        t.cursor = 8;
        for _ in 0..3 {
            t.delete_backward().unwrap();
        }
        assert_eq!(seq(&t), "AAAACGGGG");
        assert_eq!(t.cursor, 5);
        assert!(t.undo());
        assert_eq!(seq(&t), "AAAACCCCGGGG");
        assert_eq!(t.cursor, 8);
        // The Delete key removes the base after the cursor instead.
        t.cursor = 0;
        t.delete_forward().unwrap();
        t.delete_forward().unwrap();
        assert_eq!(seq(&t), "AACCCCGGGG");
        assert!(t.undo());
        assert_eq!(seq(&t), "AAAACCCCGGGG");
    }

    #[test]
    fn typing_over_a_selection_replaces_it() {
        let mut t = tab();
        t.select(Some(Span::new(4, 8)));
        t.type_bases(b"T").unwrap();
        assert_eq!(seq(&t), "AAAATGGGG");
        assert_eq!((t.selection, t.cursor), (None, 5));
    }

    #[test]
    fn cut_takes_the_bases_and_undo_returns_them() {
        let mut t = tab();
        t.select(Some(Span::new(0, 4)));
        assert_eq!(t.cut().unwrap(), Some(b"AAAA".to_vec()));
        assert_eq!(seq(&t), "CCCCGGGG");
        assert!(t.undo());
        assert_eq!(seq(&t), "AAAACCCCGGGG");
        t.select(None);
        assert_eq!(t.cut().unwrap(), None, "nothing selected, nothing cut");
    }

    #[test]
    fn select_all_then_paste_replaces_the_whole_sequence() {
        let mut t = tab();
        t.select_all();
        assert_eq!(t.selection, Some(Span::new(0, 12)));
        let bases = t.clean_input(">new\nttt ggg\n").unwrap();
        t.paste_bases(&bases).unwrap();
        assert_eq!(seq(&t), "TTTGGG");
        assert_eq!(t.cursor, 6);
        assert!(t.undo());
        assert_eq!(seq(&t), "AAAACCCCGGGG");
        assert_eq!(t.doc.features.len(), 1);
    }

    #[test]
    fn saving_never_drops_annotations() {
        // A FASTA file gets a feature: Save writes a SnapGene file beside it and leaves the FASTA
        // file alone, because FASTA cannot hold features.
        let dir = std::env::temp_dir().join(format!("slopgene-save-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("insert.fasta");
        std::fs::write(&fasta, ">insert\nAAAACCCCGGGGTTTT\n").unwrap();
        let mut t = DocTab::open(&fasta).unwrap();
        assert_eq!(t.save_keeping_everything().unwrap(), Saved::InPlace(fasta.clone()), "nothing to lose yet");
        t.edit(|d| {
            d.features.push(Feature::new("f", "misc_feature", Span::new(2, 8), Strand::Forward));
            Ok(())
        })
        .unwrap();
        let beside = dir.join("insert.dna");
        match t.save_keeping_everything().unwrap() {
            Saved::Beside { path, original, format, lost } => {
                assert_eq!((path, original, format, lost), (beside.clone(), fasta.clone(), io::Format::Fasta, vec!["features"]));
            }
            other => panic!("{other:?}"),
        }
        assert_eq!(std::fs::read_to_string(&fasta).unwrap(), ">insert\nAAAACCCCGGGGTTTT\n", "the FASTA file is unchanged");
        assert_eq!(io::read_one(&beside).unwrap().features.len(), 1, "the SnapGene file has the feature");
        assert_eq!(t.path.as_deref(), Some(beside.as_path()), "the document now lives in the SnapGene file");
        // A second FASTA file whose SnapGene name is taken asks where to save instead.
        let other = dir.join("insert.fa");
        std::fs::write(&other, ">insert\nACGTACGT\n").unwrap();
        let mut t2 = DocTab::open(&other).unwrap();
        t2.edit(|d| {
            d.features.push(Feature::new("g", "misc_feature", Span::new(0, 4), Strand::Forward));
            Ok(())
        })
        .unwrap();
        assert!(matches!(t2.save_keeping_everything().unwrap(), Saved::NeedsLocation { .. }));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn pasted_text_is_checked() {
        let t = tab();
        assert_eq!(t.clean_input("acg t\n12 nn").unwrap(), b"ACGTNN".to_vec());
        assert!(t.clean_input("ACGT hello").is_err(), "letters that are not codes are refused");
    }
}
