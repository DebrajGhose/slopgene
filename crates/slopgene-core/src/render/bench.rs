//! The protocol bench as a drawing.
//!
//! A cloning protocol is shown as a bench: each map heads a column as a small ring, its cut or
//! its PCR sits under it, and arrows carry the pieces kept down to the join, which draws the new
//! plasmid. This module lays the bench out ([`board`]), colours the pieces of each cut by what the
//! picks kept ([`cut_view`]), makes the ring picture of each card ([`picture`]) and writes the
//! words on the cards and the arrows ([`card_text`], [`link_pills`]). The app draws the bench
//! from these parts and lets the person work on it, and [`render_bench`] draws the same bench as
//! a still picture, so `slopgene protocol render` and the scripts give the picture the app shows.

use super::ring::{self, RingOptions, RingPiece, Role, render_ring};
use super::scene::*;
use super::{ColorMode, SiteMark, Theme, thousands};
use crate::model::Document;
use crate::protocol::{self, Action, Issue, Orientation, PartRole, Piece, Pin, Protocol, Run, Status, StepRun, split_reference};
use crate::span::Span;
use std::f32::consts::TAU;

/// A card's size on the bench at scale 1, the gaps between columns and rows, and the margin
/// around the bench.
pub const CARD: [f32; 2] = [250.0, 232.0];
pub const COLUMN_GAP: f32 = 70.0;
pub const ROW_GAP: f32 = 54.0;
pub const MARGIN: f32 = 22.0;

/// The heights of a card's title line and status line. The ring fills the rest.
pub const CARD_TITLE: f32 = 28.0;
pub const CARD_STATUS: f32 = 26.0;

/// The widths at scale 1 of the places on a join's card for one plasmid each: the product or an
/// option, and another plasmid of the tube, which is drawn smaller. A card with one plasmid is as
/// wide as any other card.
pub const MAIN_SLOT: f32 = CARD[0] - 24.0;
pub const OTHER_SLOT: f32 = 150.0;

/// The most plasmids a join's card shows. Its status line counts the rest.
pub const SHOWN: usize = 4;

/// What a card on the bench stands for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
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
pub struct Card {
    pub kind: Kind,
    /// The step the card stands for. A place to start or to put the insert has none.
    pub step: Option<usize>,
    /// The picks folded into a cut's card, by step.
    pub picks: Vec<usize>,
    /// Where the card sits: a column, between two for a join of both, and a row.
    pub column: f32,
    pub row: usize,
    /// How many columns wide it is, around its column: 1 for a card, more for a join that shows
    /// several plasmids.
    pub width: f32,
}

/// An arrow from one card to another, with the steps on the way (a dephosphorylation, a flip).
#[derive(Clone, Debug, PartialEq)]
pub struct Link {
    pub from: usize,
    pub to: usize,
    pub treatments: Vec<usize>,
}

/// The cards and arrows of a protocol, and the card that shows each step.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Board {
    pub cards: Vec<Card>,
    pub links: Vec<Link>,
    pub home: Vec<Option<usize>>,
}

impl Board {
    /// The card that shows step `step`.
    pub fn card_of(&self, step: usize) -> Option<usize> {
        self.home.get(step).copied().flatten()
    }

    /// The cut card that a card leads into, such as a plasmid's cut.
    pub fn cut_after(&self, card: usize) -> Option<usize> {
        self.links.iter().find(|l| l.from == card && self.cards[l.to].kind == Kind::Cut).map(|l| l.to)
    }

    /// The plasmid card that leads into cut card `card`, the other way round from `cut_after`.
    pub fn plasmid_before(&self, card: usize) -> Option<usize> {
        self.links.iter().find(|l| l.to == card && self.cards[l.from].kind == Kind::Plasmid).map(|l| l.from)
    }

    /// The PCR card that a card leads into, such as a plasmid's PCR.
    pub fn copy_after(&self, card: usize) -> Option<usize> {
        self.links.iter().find(|l| l.from == card && self.cards[l.to].kind == Kind::Copy).map(|l| l.to)
    }

    /// The card that leads into card `card`, whatever it is: the plasmid or the PCR before a cut.
    pub fn before(&self, card: usize) -> Option<usize> {
        self.links.iter().find(|l| l.to == card).map(|l| l.from)
    }

    /// The step whose part a plasmid card's map sets: the PCR after the plasmid, or else its cut.
    pub fn part_step_after(&self, card: usize) -> Option<usize> {
        self.copy_after(card).or_else(|| self.cut_after(card)).and_then(|c| self.cards[c].step)
    }

    /// The bench's size at scale 1, as width and height.
    pub fn size(&self) -> [f32; 2] {
        size_of(self.cards.iter())
    }

    /// Where a card sits at scale 1, as left, top, width and height. A card wider than one column
    /// spreads evenly to both sides of its column.
    pub fn rect(&self, card: usize) -> [f32; 4] {
        let c = &self.cards[card];
        let left = c.column - (c.width - 1.0) / 2.0;
        [MARGIN + left * (CARD[0] + COLUMN_GAP), MARGIN + c.row as f32 * (CARD[1] + ROW_GAP), c.width * (CARD[0] + COLUMN_GAP) - COLUMN_GAP, CARD[1]]
    }
}

/// The size of a bench that holds these cards, at scale 1.
fn size_of<'a>(cards: impl Iterator<Item = &'a Card>) -> [f32; 2] {
    let (mut columns, mut rows) = (1.0_f32, 1usize);
    // The right edge of the card that reaches furthest, in columns, and the lowest row.
    for card in cards {
        columns = columns.max(card.column + (card.width + 1.0) / 2.0);
        rows = rows.max(card.row + 1);
    }
    let rows = rows as f32;
    [2.0 * MARGIN + columns * CARD[0] + (columns - 1.0) * COLUMN_GAP, 2.0 * MARGIN + rows * CARD[1] + (rows - 1.0) * ROW_GAP]
}

/// How many plasmids a join's card shows: its product or its options in full-size places, and
/// the other plasmids of its tube in smaller places after them, [`SHOWN`] in all.
pub fn join_slots(result: Option<&StepRun>) -> (usize, usize) {
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
pub fn join_width(result: Option<&StepRun>) -> f32 {
    let (main, others) = join_slots(result);
    let pixels = (24.0 + main as f32 * MAIN_SLOT + others as f32 * OTHER_SLOT).max(CARD[0]);
    (pixels + COLUMN_GAP) / (CARD[0] + COLUMN_GAP)
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
pub fn board(protocol: &Protocol, run: &Run) -> Board {
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
pub struct CutView {
    pub parent: Document,
    pub pieces: Vec<(Span, Role, Option<usize>)>,
    pub cuts: Vec<SiteMark>,
    pub region: Option<Span>,
    pub top: usize,
    /// Whether the enzymes have cut, rather than the part being shown before any cut.
    pub cut: bool,
}

/// The bases of a molecule outside `region`: one stretch round a circle, up to two on a line.
pub fn complement(region: Span, len: usize, circular: bool) -> Vec<Span> {
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
pub fn cut_view(protocol: &Protocol, run: &Run, card: &Card) -> Option<CutView> {
    let step = card.step?;
    let result = run.steps.get(step)?;
    let parent = result.inputs.first()?.doc.clone();
    let len = parent.len().max(1);
    let region = result.findings.region;
    // The pieces the picks kept, and those a pick waits to choose among, by where they lay.
    let picks: Vec<&StepRun> = card.picks.iter().filter_map(|&p| run.steps.get(p)).collect();
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
        let role = protocol.steps.get(step).and_then(|s| protocol::part_role(protocol, &s.id));
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
pub fn copy_region(run: &Run, step: usize) -> Option<(Span, usize, bool)> {
    let result = run.steps.get(step)?;
    let template = result.inputs.first()?;
    Some((result.findings.region?, template.doc.len().max(1), template.doc.is_circular()))
}

/// The tails on the product of a PCR, in bases at each end, as the run found them.
pub fn tail_lengths(result: &StepRun) -> (usize, usize) {
    let [left, right] = result.findings.tails;
    (left, right)
}

/// The middle base of the last part of a product (the insert, in a template's join), to turn to
/// the top.
pub fn product_top(parts: &[Span], len: usize) -> usize {
    parts.last().map_or(0, |s| (s.start + s.len() / 2) % len.max(1))
}

/// What one plasmid that a join shows stands for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Made {
    /// The join's product.
    Product,
    /// One of the options of a join that waits for a choice, from 0.
    Option(usize),
    /// Another plasmid the join's tube makes, by its place among the by-products, from 0.
    Other(usize),
}

/// One plasmid a join shows, with the stretch of it each piece became.
pub struct Shown<'a> {
    pub made: Made,
    pub piece: &'a Piece,
    pub parts: &'a [Span],
}

/// The plasmids a join shows: its product, or its options while it waits for a choice, and then
/// the other plasmids its tube makes. A join that is not made shows none.
pub fn join_products(result: &StepRun) -> Vec<Shown<'_>> {
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
pub fn made_ring(shown: &Shown, theme: Theme, colors: ColorMode) -> RingOptions {
    let len = shown.piece.doc.len().max(1);
    let other = matches!(shown.made, Made::Other(_));
    let role = if other { Role::Plain } else { Role::Kept };
    let title = match shown.made {
        Made::Product => Some("New plasmid".to_string()),
        Made::Option(ii) => Some(format!("Option {}", ii + 1)),
        Made::Other(_) => None,
    };
    RingOptions {
        theme,
        colors,
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
pub fn short_name(protocol: &Protocol, run: &Run, reference: &str) -> String {
    let words = protocol::piece_name(protocol, run, reference);
    words.strip_prefix("the ").map(str::to_string).unwrap_or(words)
}

/// A few words for the line under one plasmid of join step `step`: an option's number and which
/// way its insert goes, or how another plasmid of the tube differs from the product.
pub fn made_caption(protocol: &Protocol, run: &Run, step: usize, made: Made) -> String {
    let (Some(result), Some(inputs)) = (run.steps.get(step), protocol.steps.get(step).map(|s| &s.inputs)) else { return String::new() };
    let name = |input: usize| inputs.get(input).map(|r| short_name(protocol, run, r)).unwrap_or_default();
    match made {
        Made::Product => String::new(),
        Made::Option(ii) => {
            // The ring already says which option it is, so the words say which way the inserts
            // go: each by name with one insert, and the ones turned round with several.
            let Some(Pin::Orientations(ways)) = result.choices.get(ii).map(|c| &c.pin) else { return String::new() };
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
pub struct Picture {
    pub doc: Document,
    pub options: RingOptions,
}

/// The ring a card draws: the map whole for a plasmid, with its pieces coloured once its cut is
/// chosen; the pieces pulled apart for a cut; the new plasmid for a join. `detail` asks for the
/// labels of the large view. `None` when there is nothing to draw yet.
pub fn picture(protocol: &Protocol, run: &Run, board: &Board, card: usize, detail: bool, theme: Theme, colors: ColorMode) -> Option<Picture> {
    let card_ref = &board.cards[card];
    let step = card_ref.step?;
    let result = run.steps.get(step)?;
    let mut options = RingOptions { theme, colors, gaps: true, ..RingOptions::default() };
    match card_ref.kind {
        Kind::Plasmid => {
            let doc = result.outputs.first()?.doc.clone();
            options.title = Some(doc.name.clone());
            options.subtitle = Some(format!("{} bp", thousands(doc.len())));
            if let Some(view) = board.cut_after(card).and_then(|c| cut_view(protocol, run, &board.cards[c])) {
                options.top = view.top;
                options.gaps = view.cut;
                options.marks = view.cuts;
                options.pieces = view.pieces.iter().map(|(span, role, _)| RingPiece { span: *span, role: *role, lift: 0.0, label: None }).collect();
            } else if let Some((region, len, circular)) = board.copy_after(card).and_then(|c| board.cards[c].step).and_then(|s| copy_region(run, s)) {
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
            let view = cut_view(protocol, run, card_ref).filter(|v| v.cut)?;
            let len = view.parent.len();
            let largest = view.pieces.iter().map(|p| p.0.len()).max().unwrap_or(0);
            let part = match &protocol.steps[step].action {
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
            Some(Picture { doc: first.piece.doc.clone(), options: made_ring(first, theme, colors) })
        }
        Kind::Other | Kind::Start | Kind::PutInto | Kind::AddPlasmid => None,
    }
}

/// The colour a line of words is drawn in, by what it says: plain text, a quieter note, a faint
/// hint, and the three states of a step.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Tone {
    Ink,
    #[default]
    Muted,
    Faint,
    Ok,
    Warn,
    Error,
}

/// The tone of a step's state: done, waiting for a choice, failed, or skipped.
pub fn status_tone(status: Status) -> Tone {
    match status {
        Status::Done => Tone::Ok,
        Status::NeedsChoice => Tone::Warn,
        Status::Failed => Tone::Error,
        Status::Skipped => Tone::Faint,
    }
}

/// The words on a card: a verb in small capitals, what it acts on, and a status line in its tone.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CardText {
    pub verb: String,
    pub object: String,
    pub line: String,
    pub tone: Tone,
}

/// A problem in a few words for a card's status line, from its code. The panel gives it whole.
pub fn short_problem(issue: &Issue) -> String {
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

/// A card's title (a verb, then what it acts on) and its status line with the line's tone.
pub fn card_text(protocol: &Protocol, run: &Run, board: &Board, card: usize) -> CardText {
    let c = &board.cards[card];
    let result = c.step.and_then(|s| run.steps.get(s));
    let status = result.map(|r| r.status);
    let problem = || result.and_then(|r| r.error.as_ref()).map(short_problem).unwrap_or_default();
    let text = |verb: &str, object: String, line: String, tone: Tone| CardText { verb: verb.into(), object, line, tone };
    match c.kind {
        Kind::Plasmid => {
            let made = result.and_then(|r| r.outputs.first());
            let name = made.map(|p| p.doc.name.clone()).unwrap_or_else(|| "choose a map".into());
            match (status, made) {
                (Some(Status::Done), Some(p)) => {
                    let shape = if p.doc.is_circular() { "circular" } else { "linear" };
                    text("START", name, format!("{} bp · {shape}", thousands(p.doc.len())), Tone::Muted)
                }
                _ => text("START", name, problem(), Tone::Error),
            }
        }
        Kind::Cut => {
            let Some(Action::Digest(params)) = c.step.map(|s| &protocol.steps[s].action) else { return CardText::default() };
            let enzymes = if params.enzymes.is_empty() { "choose enzymes".to_string() } else { params.enzymes.join(" + ") };
            let picks: Vec<&StepRun> = c.picks.iter().filter_map(|&p| run.steps.get(p)).collect();
            let (line, tone) = match status {
                Some(Status::NeedsChoice) => {
                    let around = params.around.as_deref().map(|p| format!(" around {p}")).unwrap_or_default();
                    (format!("Choose where to cut{around}"), Tone::Warn)
                }
                Some(Status::Done) if picks.iter().any(|p| p.status == Status::NeedsChoice) => ("Choose the piece to keep".into(), Tone::Warn),
                Some(Status::Done) => match picks.iter().find(|p| p.status == Status::Done).and_then(|p| p.outputs.first()) {
                    Some(kept) => (format!("keeps {} bp · {} pieces", thousands(kept.doc.len()), result.map_or(0, |r| r.outputs.len())), Tone::Muted),
                    None => (picks.iter().find_map(|p| p.error.as_ref()).map(short_problem).unwrap_or_default(), Tone::Error),
                },
                Some(Status::Failed) => (problem(), Tone::Error),
                _ => ("waits for the step above".into(), Tone::Faint),
            };
            text("CUT", enzymes, line, tone)
        }
        Kind::Copy => {
            let Some(Action::Amplify(params)) = c.step.map(|s| &protocol.steps[s].action) else { return CardText::default() };
            let object = match params.region.as_deref() {
                Some(part) => format!("{part} by PCR"),
                None => "by PCR".to_string(),
            };
            let tails = params.enzymes();
            let (line, tone) = match status {
                Some(Status::NeedsChoice) => ("Mark the part to copy".into(), Tone::Warn),
                Some(Status::Done) => {
                    let size = result.and_then(|r| r.outputs.first()).map(|p| format!("{} bp product", thousands(p.doc.len()))).unwrap_or_default();
                    match tails.as_slice() {
                        [] => (format!("{size} · choose the tail enzymes"), Tone::Warn),
                        list => (format!("{size} · {} tails", list.join(" + ")), Tone::Muted),
                    }
                }
                Some(Status::Failed) => (problem(), Tone::Error),
                _ => ("waits for the step above".into(), Tone::Faint),
            };
            text("COPY", object, line, tone)
        }
        Kind::Join => {
            // The other plasmids of the tube are counted, as the card may not have room for all.
            let others = match result.map_or(0, |r| r.findings.byproducts.len()) {
                0 => String::new(),
                n => format!(" · {n} more can form"),
            };
            let (line, tone) = match status {
                Some(Status::Done) => {
                    let made = result.and_then(|r| r.outputs.first());
                    let shape = made.map_or("", |p| if p.doc.is_circular() { "circular" } else { "linear" });
                    (made.map(|p| format!("{} bp · {shape}{others}", thousands(p.doc.len()))).unwrap_or_default(), Tone::Muted)
                }
                Some(Status::NeedsChoice) => (format!("Choose the direction{others}"), Tone::Warn),
                Some(Status::Failed) => (problem(), Tone::Error),
                _ => ("waits for the pieces".into(), Tone::Faint),
            };
            text("JOIN", "new plasmid".into(), line, tone)
        }
        Kind::Other => {
            let step = c.step.map(|s| &protocol.steps[s]);
            let title = step.map(|s| s.action.title().to_uppercase()).unwrap_or_default();
            let line = status.map(|s| s.name().to_string()).unwrap_or_default();
            CardText { verb: title, object: step.map(|s| s.id.clone()).unwrap_or_default(), line, tone: status.map_or(Tone::Faint, status_tone) }
        }
        Kind::Start => text("START", "with a plasmid".into(), String::new(), Tone::Muted),
        Kind::PutInto => text("PUT IT INTO", "a plasmid".into(), "that the kept piece goes into".into(), Tone::Faint),
        Kind::AddPlasmid => text("ADD", "another plasmid".into(), "to cut a piece from".into(), Tone::Faint),
    }
}

/// What a pill on an arrow leads to: the plasmid whose map chooses a cut's enzymes, the PCR whose
/// tails are chosen, or a treatment step.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PillTarget {
    Cut(usize),
    Copy(usize),
    Treatment(usize),
}

/// A pill on an arrow between two cards: what happens on the way, in its tone.
#[derive(Clone, Debug, PartialEq)]
pub struct LinkPill {
    pub text: String,
    pub tone: Tone,
    pub target: PillTarget,
}

/// The pills on an arrow: the enzymes of the cut it leads into (or of the trim after a PCR), the
/// tails of the PCR it leads into, and the treatments on the way (a dephosphorylation, a flip).
pub fn link_pills(protocol: &Protocol, board: &Board, link: &Link) -> Vec<LinkPill> {
    let mut pills = Vec::new();
    if board.cards[link.to].kind == Kind::Cut
        && let Some(Action::Digest(params)) = board.cards[link.to].step.map(|s| &protocol.steps[s].action)
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
        let tone = if params.enzymes.is_empty() { Tone::Warn } else { Tone::Ink };
        pills.push(LinkPill { text, tone, target: if trims { PillTarget::Copy(from) } else { PillTarget::Cut(from) } });
    }
    if board.cards[link.to].kind == Kind::Copy
        && let Some(copy) = board.cards[link.to].step
        && let Action::Amplify(params) = &protocol.steps[copy].action
    {
        let tails = params.enzymes();
        let text = if tails.is_empty() { "PCR: choose the tail enzymes".to_string() } else { format!("PCR with {} tails", tails.join(" + ")) };
        pills.push(LinkPill { text, tone: if tails.is_empty() { Tone::Warn } else { Tone::Ink }, target: PillTarget::Copy(copy) });
    }
    for &step in &link.treatments {
        let text = match protocol.steps[step].action {
            Action::Dephosphorylate(_) => "no 5' phosphates",
            Action::Flip(_) => "turned round",
            _ => continue,
        };
        pills.push(LinkPill { text: text.to_string(), tone: Tone::Muted, target: PillTarget::Treatment(step) });
    }
    pills
}

/// The colours of the bench beyond the theme's: a card and its edge, a pill, faint words, and
/// the three states of a step. They match the app's interface colours.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BenchColors {
    pub card: Rgba,
    pub line: Rgba,
    pub pill: Rgba,
    pub faint_text: Rgba,
    pub ok: Rgba,
    pub warn: Rgba,
    pub error: Rgba,
}

impl BenchColors {
    pub fn light() -> BenchColors {
        BenchColors {
            card: Rgba::rgb(0xFF, 0xFF, 0xFF),
            line: Rgba::rgb(0xDA, 0xDD, 0xE3),
            pill: Rgba::rgb(0xEF, 0xF1, 0xF4),
            faint_text: Rgba::rgb(0xA0, 0xA4, 0xAB),
            ok: Rgba::rgb(0x2F, 0x80, 0x55),
            warn: Rgba::rgb(0xB0, 0x6A, 0x00),
            error: Rgba::rgb(0xB4, 0x3C, 0x36),
        }
    }

    pub fn dark() -> BenchColors {
        BenchColors {
            card: Rgba::rgb(0x2A, 0x2B, 0x30),
            line: Rgba::rgb(0x34, 0x35, 0x3A),
            pill: Rgba::rgb(0x22, 0x23, 0x27),
            faint_text: Rgba::rgb(0x66, 0x67, 0x6C),
            ok: Rgba::rgb(0x6C, 0xC0, 0x8E),
            warn: Rgba::rgb(0xE0, 0xA0, 0x40),
            error: Rgba::rgb(0xE8, 0x7C, 0x72),
        }
    }
}

/// How to draw a bench as a still picture.
#[derive(Clone, Debug, PartialEq)]
pub struct BenchOptions {
    pub theme: Theme,
    pub colors: ColorMode,
    pub bench: BenchColors,
    /// The size of the words on the arrows, in pixels. The other sizes follow it.
    pub font_size: f32,
}

impl Default for BenchOptions {
    fn default() -> Self {
        BenchOptions { theme: Theme::light(), colors: ColorMode::Palette, bench: BenchColors::light(), font_size: 11.0 }
    }
}

impl BenchOptions {
    /// The options for a dark picture.
    pub fn dark() -> BenchOptions {
        BenchOptions { theme: Theme::dark(), bench: BenchColors::dark(), ..BenchOptions::default() }
    }

    /// The colour of words in this tone.
    pub fn tone(&self, tone: Tone) -> Rgba {
        match tone {
            Tone::Ink => self.theme.ink,
            Tone::Muted => self.theme.muted,
            Tone::Faint => self.bench.faint_text,
            Tone::Ok => self.bench.ok,
            Tone::Warn => self.bench.warn,
            Tone::Error => self.bench.error,
        }
    }
}

/// Draw the bench of a protocol as a still picture, as the app shows it: the cards with their
/// rings and words, and the arrows between them with their pills. The places the app offers to
/// start, to put the insert somewhere or to add a map are left out, as they are for the mouse.
pub fn render_bench(protocol: &Protocol, run: &Run, options: &BenchOptions, measure: &dyn TextMeasure) -> Scene {
    let board = board(protocol, run);
    let shown: Vec<usize> = (0..board.cards.len()).filter(|&ii| board.cards[ii].step.is_some()).collect();
    let [width, height] = size_of(shown.iter().map(|&ii| &board.cards[ii]));
    let mut scene = Scene::new(width, height, options.theme.background);
    // The arrows go first, so the cards lie over their ends.
    for (ii, link) in board.links.iter().enumerate() {
        if !shown.contains(&link.from) || !shown.contains(&link.to) {
            continue;
        }
        // Several arrows into one card spread across its top edge.
        let incoming: Vec<usize> = board.links.iter().enumerate().filter(|(_, l)| l.to == link.to).map(|(jj, _)| jj).collect();
        let slot = (incoming.iter().position(|&jj| jj == ii).unwrap_or(0), incoming.len());
        draw_arrow(&mut scene, protocol, &board, link, slot, options, measure);
    }
    for &ii in &shown {
        draw_card(&mut scene, protocol, run, &board, ii, options, measure);
    }
    scene
}

/// An arrow between two cards: straight down within a column, curved into a join, with the
/// pills that say what happens on the way at its middle.
fn draw_arrow(scene: &mut Scene, protocol: &Protocol, board: &Board, link: &Link, slot: (usize, usize), options: &BenchOptions, measure: &dyn TextMeasure) {
    let (from, to) = (board.rect(link.from), board.rect(link.to));
    let start = [from[0] + from[2] / 2.0, from[1] + from[3]];
    let end = [to[0] + to[2] * (slot.0 + 1) as f32 / (slot.1 + 1) as f32, to[1]];
    let tip = [end[0], end[1] - 2.0];
    let color = options.bench.faint_text;
    let points = if (start[0] - end[0]).abs() < 1.0 {
        vec![start, tip]
    } else {
        let bend = (end[1] - start[1]) * 0.5;
        bezier([start, [start[0], start[1] + bend], [tip[0], tip[1] - bend], tip], 24)
    };
    scene.items.push(Item { shape: Shape::Polyline { points, color, width: 1.6 }, tag: Tag::None });
    let head = 6.0;
    let points = vec![tip, [tip[0] - head * 0.7, tip[1] - head], [tip[0] + head * 0.7, tip[1] - head]];
    scene.items.push(Item { shape: Shape::Polygon { points, fill: color, stroke: None }, tag: Tag::None });
    // The pills, one under the other at the middle of the arrow.
    let pills = link_pills(protocol, board, link);
    let middle = [(start[0] + end[0]) / 2.0, (start[1] + end[1]) / 2.0];
    let size = options.font_size;
    for (ii, pill) in pills.iter().enumerate() {
        let at = [middle[0], middle[1] + (ii as f32 - (pills.len() as f32 - 1.0) / 2.0) * 22.0];
        let (w, h) = (measure.text_width(&pill.text, size, false) + 14.0, size * 1.3 + 6.0);
        let points = rounded_rect([at[0] - w / 2.0, at[1] - h / 2.0, w, h], (h / 2.0).min(10.0));
        scene.items.push(Item { shape: Shape::Polygon { points, fill: options.bench.pill, stroke: Some((options.bench.line, 1.0)) }, tag: Tag::None });
        scene.items.push(words(at, &pill.text, size, options.tone(pill.tone), Anchor::Middle, false));
    }
}

/// One card: its frame, its title with the state dot, its ring (or a dashed place for one with
/// what it waits for) and its status line.
fn draw_card(scene: &mut Scene, protocol: &Protocol, run: &Run, board: &Board, card: usize, options: &BenchOptions, measure: &dyn TextMeasure) {
    let c = &board.cards[card];
    let [x, y, w, h] = board.rect(card);
    let text = card_text(protocol, run, board, card);
    let status = c.step.and_then(|s| run.steps.get(s)).map(|r| r.status);
    let edge = match status {
        Some(Status::NeedsChoice) => (options.bench.warn, 1.5),
        Some(Status::Failed) => (options.bench.error, 1.5),
        _ => (options.bench.line, 1.0),
    };
    scene.items.push(Item { shape: Shape::Polygon { points: rounded_rect([x, y, w, h], 10.0), fill: options.bench.card, stroke: Some(edge) }, tag: Tag::None });

    // The title: the verb in small capitals, then what it acts on, and the state as a dot.
    let title_y = y + CARD_TITLE * 0.55;
    let (verb_size, object_size) = (options.font_size - 0.5, options.font_size + 1.5);
    let verb_width = measure.text_width(&text.verb, verb_size, true);
    scene.items.push(words([x + 12.0, title_y], &text.verb, verb_size, options.theme.accent, Anchor::Start, true));
    let object = fit(measure, &text.object, object_size, false, w - 12.0 - verb_width - 7.0 - 62.0);
    scene.items.push(words([x + 12.0 + verb_width + 7.0, title_y], &object, object_size, options.theme.ink, Anchor::Start, false));
    if let Some(status) = status {
        let fill = Some(options.tone(status_tone(status)));
        scene.items.push(Item { shape: Shape::Circle { center: [x + w - 13.0, title_y], radius: 4.0, fill, stroke: None }, tag: Tag::None });
    }

    // The ring. A cut's pieces pull up and out, so its ring sits a little lower to leave room for
    // the top piece's label.
    let drop = if c.kind == Kind::Cut { 12.0 } else { 0.0 };
    let ring = [x, y + CARD_TITLE + drop, w, h - CARD_TITLE - CARD_STATUS - drop];
    match picture(protocol, run, board, card, false, options.theme, options.colors) {
        Some(shown) => {
            let mut ring_options = shown.options;
            ring_options.width = ring[2];
            ring_options.height = ring[3];
            ring_options.font_size = options.font_size - 0.5;
            let drawn = render_ring(&shown.doc, &ring_options, measure);
            scene.items.extend(shifted(drawn, ring[0], ring[1]));
        }
        None => {
            // A dashed circle where the ring will be, with what it waits for inside.
            let center = [ring[0] + ring[2] / 2.0, ring[1] + ring[3] / 2.0];
            let radius = ring::radius(ring[2], ring[3]);
            let dashes = 28;
            for ii in 0..dashes {
                let a0 = ii as f32 / dashes as f32 * TAU;
                let points = arc_points(center, radius, a0, a0 + 0.55 * TAU / dashes as f32);
                scene.items.push(Item { shape: Shape::Polyline { points, color: options.bench.faint_text, width: 1.5 }, tag: Tag::None });
            }
            let waiting = match (c.kind, status) {
                (Kind::Cut, Some(Status::NeedsChoice)) => "choose\nwhere to cut",
                (Kind::Copy, Some(Status::NeedsChoice)) => "mark\nthe part to copy",
                (Kind::Join | Kind::Cut | Kind::Copy, Some(Status::Failed)) => "failed:\nsee the step",
                (Kind::Join | Kind::Cut | Kind::Copy, _) => "waits for\nthe steps above",
                _ => "",
            };
            let lines: Vec<&str> = waiting.lines().collect();
            let pitch = options.font_size * 1.3;
            for (ii, line) in lines.iter().enumerate() {
                let at = [center[0], center[1] + (ii as f32 - (lines.len() as f32 - 1.0) / 2.0) * pitch];
                scene.items.push(words(at, line, options.font_size, options.bench.faint_text, Anchor::Middle, false));
            }
        }
    }

    // The status line.
    let line_size = options.font_size + 0.5;
    let line = fit(measure, &text.line, line_size, false, w - 20.0);
    scene.items.push(words([x + w / 2.0, y + h - CARD_STATUS * 0.5], &line, line_size, options.tone(text.tone), Anchor::Middle, false));
}

/// A text item.
fn words(pos: [f32; 2], text: &str, size: f32, color: Rgba, anchor: Anchor, bold: bool) -> Item {
    Item { shape: Shape::Text { pos, text: text.to_string(), size, color, anchor, bold }, tag: Tag::None }
}

/// Shorten text with an ellipsis until it fits `width` at this size.
fn fit(measure: &dyn TextMeasure, text: &str, size: f32, bold: bool, width: f32) -> String {
    if measure.text_width(text, size, bold) <= width {
        return text.to_string();
    }
    let mut chars: Vec<char> = text.chars().collect();
    while chars.len() > 1 {
        chars.pop();
        let candidate: String = chars.iter().collect::<String>() + "…";
        if measure.text_width(&candidate, size, bold) <= width {
            return candidate;
        }
    }
    String::new()
}

/// The outline of a rectangle (left, top, width, height) with rounded corners of this radius,
/// as a convex polygon.
fn rounded_rect(rect: [f32; 4], radius: f32) -> Vec<[f32; 2]> {
    let [x, y, w, h] = rect;
    let r = radius.min(w / 2.0).min(h / 2.0);
    // The centres of the four corner arcs, clockwise from the top left, and the angle each arc
    // starts at (clockwise from 12 o'clock).
    let corners = [([x + w - r, y + r], 0.0), ([x + w - r, y + h - r], 0.25), ([x + r, y + h - r], 0.5), ([x + r, y + r], 0.75)];
    let mut points = Vec::with_capacity(28);
    for (center, start) in corners {
        for ii in 0..=6 {
            let angle = (start + ii as f32 / 6.0 * 0.25) * TAU;
            points.push([center[0] + r * angle.sin(), center[1] - r * angle.cos()]);
        }
    }
    points
}

/// Points along a cubic Bezier curve.
fn bezier(p: [[f32; 2]; 4], steps: usize) -> Vec<[f32; 2]> {
    (0..=steps)
        .map(|ii| {
            let t = ii as f32 / steps as f32;
            let u = 1.0 - t;
            let weights = [u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t];
            let mut point = [0.0, 0.0];
            for (w, q) in weights.iter().zip(p.iter()) {
                point[0] += w * q[0];
                point[1] += w * q[1];
            }
            point
        })
        .collect()
}

/// The items of a drawing moved right by `dx` and down by `dy`, to place it inside another.
fn shifted(scene: Scene, dx: f32, dy: f32) -> Vec<Item> {
    let moved = |p: [f32; 2]| [p[0] + dx, p[1] + dy];
    scene
        .items
        .into_iter()
        .map(|item| {
            let shape = match item.shape {
                Shape::Polyline { points, color, width } => Shape::Polyline { points: points.into_iter().map(moved).collect(), color, width },
                Shape::Polygon { points, fill, stroke } => Shape::Polygon { points: points.into_iter().map(moved).collect(), fill, stroke },
                Shape::Circle { center, radius, fill, stroke } => Shape::Circle { center: moved(center), radius, fill, stroke },
                Shape::Text { pos, text, size, color, anchor, bold } => Shape::Text { pos: moved(pos), text, size, color, anchor, bold },
            };
            Item { shape, tag: item.tag }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzymes::EnzymeDb;
    use crate::model::{Feature, Strand, Topology};
    use crate::protocol::{CloningTemplate, Loader, OpenParams, run};
    use crate::render::svg::to_svg;

    /// A loader that has no files: the maps are carried in the steps.
    struct NoFiles;

    impl Loader for NoFiles {
        fn load(&self, path: &str) -> crate::Result<Document> {
            Err(crate::Error::Invalid(format!("{path} is not here")))
        }
    }

    /// Random bases without the sites of the enzymes used here, so only planted sites cut, and
    /// without repeats, so primers anneal where they are designed to.
    fn filler(seed: u64, n: usize) -> Vec<u8> {
        let sites: [&[u8]; 3] = [b"GAATTC", b"GGATCC", b"AGATCT"];
        let mut x = seed;
        let mut seq: Vec<u8> = Vec::with_capacity(n);
        while seq.len() < n {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            seq.push(b"ACGT"[(x >> 33) as usize % 4]);
            // Drop the last base whenever it completes a site.
            if seq.len() >= 6 && sites.iter().any(|site| seq.ends_with(site)) {
                seq.pop();
            }
        }
        seq
    }

    /// A donor with geneX between EcoRI and BamHI sites, and a vector with EcoRI and BglII sites.
    fn protocol_in_memory(pcr: bool) -> Protocol {
        let mut donor = filler(1, 1000);
        donor[200..206].copy_from_slice(b"GAATTC");
        donor[700..706].copy_from_slice(b"GGATCC");
        let mut donor = Document::new_dna("donor", donor, Topology::Circular);
        donor.features.push(Feature::new("geneX", "CDS", Span::new(300, 600), Strand::Forward));
        let mut vector = filler(2, 3000);
        vector[1000..1006].copy_from_slice(b"AGATCT");
        vector[1100..1106].copy_from_slice(b"GAATTC");
        let vector = Document::new_dna("vector", vector, Topology::Circular);
        let carried = |doc: Document| OpenParams { path: None, document: Some(crate::json::document_to_value(&doc, true)) };
        Protocol::restriction_cloning(&CloningTemplate {
            name: "geneX into vector".into(),
            vector: Some(carried(vector)),
            donor: carried(donor),
            insert_enzymes: vec!["EcoRI".into(), "BamHI".into()],
            vector_enzymes: vec!["EcoRI".into(), "BglII".into()],
            insert: Some("geneX".into()),
            dephosphorylate: true,
            pcr,
            ..CloningTemplate::default()
        })
    }

    fn texts(scene: &Scene) -> Vec<String> {
        scene
            .items
            .iter()
            .filter_map(|i| match &i.shape {
                Shape::Text { text, .. } => Some(text.clone()),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn bench_picture_shows_the_steps() {
        let db = EnzymeDb::builtin();
        let protocol = protocol_in_memory(false);
        let ran = run(&protocol, &NoFiles, db);
        assert!(ran.steps.iter().all(|s| s.status == Status::Done), "{}", ran.text(&protocol));
        let laid = board(&protocol, &ran);
        // Two columns of a plasmid and its cut, the join between them, and the place to add a map.
        assert_eq!(
            laid.cards.iter().map(|c| c.kind).collect::<Vec<_>>(),
            vec![Kind::Plasmid, Kind::Cut, Kind::Plasmid, Kind::Cut, Kind::Join, Kind::AddPlasmid]
        );
        let scene = render_bench(&protocol, &ran, &BenchOptions::default(), &ApproxMeasure);
        // The picture holds the five step cards and no place, so it is two columns wide and three
        // rows tall.
        assert_eq!([scene.width, scene.height], size_of(laid.cards.iter().take(5)));
        assert!(scene.width > 2.0 * CARD[0] && scene.height > 2.0 * CARD[1]);
        let words = texts(&scene);
        for expected in ["START", "CUT", "JOIN", "vector", "donor", "New plasmid", "cut with EcoRI + BglII", "cut with EcoRI + BamHI", "no 5' phosphates"] {
            assert!(words.iter().any(|w| w == expected), "{expected:?} is missing from {words:?}");
        }
        // The pieces kept are green, the leftovers grey, in the rings.
        let colours: Vec<Rgba> = scene
            .items
            .iter()
            .filter_map(|i| match &i.shape {
                Shape::Polyline { color, .. } => Some(*color),
                _ => None,
            })
            .collect();
        assert!(colours.contains(&Theme::light().kept) && colours.contains(&Theme::light().left_over));
        let svg = to_svg(&scene);
        assert!(svg.contains("JOIN") && svg.contains("<svg"));
    }

    #[test]
    fn bench_picture_of_a_pcr_lane_and_of_nothing() {
        let db = EnzymeDb::builtin();
        let protocol = protocol_in_memory(true);
        let ran = run(&protocol, &NoFiles, db);
        let scene = render_bench(&protocol, &ran, &BenchOptions::dark(), &ApproxMeasure);
        let words = texts(&scene);
        assert!(words.iter().any(|w| w == "COPY") && words.iter().any(|w| w == "PCR product"), "{words:?}");
        assert_eq!(scene.background, Theme::dark().background);
        // An empty protocol draws an empty bench without complaint.
        let empty = Protocol::default();
        let scene = render_bench(&empty, &run(&empty, &NoFiles, db), &BenchOptions::default(), &ApproxMeasure);
        assert!(scene.items.is_empty() && scene.width > 0.0);
    }
}
