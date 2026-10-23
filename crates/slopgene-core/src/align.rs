//! Aligning sequencing reads to a reference and reporting the differences.
//!
//! The method is seed, chain and align, as in SnapGene's "Align to Reference" and minimap2:
//!
//! 1. **Seeding.** Every 15-mer of the read (both orientations) is looked up in the reference;
//!    each hit votes for a diagonal. The winning orientation and diagonal say where the read
//!    sits, whatever base sequencing started from.
//! 2. **Chaining.** The read's k-mer hits in the reference (anchors) are chained: the best set
//!    of anchors that runs in the same order along both sequences, scored by bases matched minus
//!    a cost for each change in diagonal. The cost grows slowly with the size of the jump, so a
//!    chain carries on across a deletion or insertion of any length.
//! 3. **Alignment.** Runs of chained anchors are exact matches. The stretches between them are
//!    aligned base by base with affine gaps (Gotoh), as are the read's two ends, with free end
//!    gaps so a Sanger read can cover part of a plasmid.
//! 4. **Differences and events.** Substitutions and small indels are listed in reference
//!    coordinates with the features they fall in, the protein change inside coding features and
//!    a flag for indels in homopolymer runs (the usual nanopore miscall). Changes of
//!    [`STRUCTURAL_MIN`] bases or more are structural events: deletions, insertions, duplications
//!    (such as a dimer), inversions and replacements, each with the features it touches.
//!
//! A whole-plasmid read is a circle opened at an arbitrary base. It is aligned twice: the second
//! time opened again in the middle of its longest exact match, so no change is split between the
//! read's two ends, and whatever remains where the ends meet is reported as well.
//!
//! The alignment itself is kept column by column ([`Alignment::aligned`]), so it can be shown
//! base by base: [`Alignment::blocks`] lays it out in lines with each structural change on a line
//! of its own, and [`Alignment::pairwise_text`] writes those lines as text.

use crate::alphabet;
use crate::model::{Document, Feature, Strand};
use crate::render::{AlignmentTrack, DifferenceKind, DifferenceMark, SegmentKind, TrackSegment, thousands};
use crate::span::Span;
use crate::translate::GeneticCode;
use crate::translate::{OrfOptions, find_orfs};
use std::collections::{HashMap, HashSet};
use std::ops::Range;

const K: usize = 15;
/// Alignment scores: match, mismatch, gap open, gap extend.
const MATCH: i32 = 2;
const MISMATCH: i32 = -3;
const GAP_OPEN: i32 = -5;
const GAP_EXTEND: i32 = -2;

/// Where a read was placed on a reference by seeding.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Placement {
    /// True when the read matches the reference's bottom strand.
    pub reverse: bool,
    /// Reference position (0-based) that read position 0 lines up with.
    pub diagonal: i64,
    /// Fraction of the read's k-mers found in the reference (0..1): a fast similarity score.
    pub containment: f64,
    /// Votes for the winning diagonal.
    pub votes: u32,
}

/// Pack a k-mer of unambiguous bases into 2 bits per base; `None` when it has another code.
fn pack(kmer: &[u8]) -> Option<u64> {
    let mut v = 0u64;
    for &b in kmer {
        let code = match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' | b'U' => 3,
            _ => return None,
        };
        v = (v << 2) | code;
    }
    Some(v)
}

/// Reverse complement of a 2-bit packed k-mer (A=0, C=1, G=2, T=3, so complement is 3 - x).
fn revcomp_code(code: u64, k: usize) -> u64 {
    let mut out = 0u64;
    let mut c = code;
    for _ in 0..k {
        out = (out << 2) | (3 - (c & 3));
        c >>= 2;
    }
    out
}

/// The strand-independent form of a k-mer: the smaller of it and its reverse complement.
fn canonical(code: u64) -> u64 {
    code.min(revcomp_code(code, K))
}

/// Distinct canonical k-mers of a sequence (through the origin when circular).
fn canonical_set(seq: &[u8], circular: bool) -> HashSet<u64> {
    let len = seq.len();
    let mut set = HashSet::new();
    if len < K {
        return set;
    }
    let starts = if circular { len } else { len - K + 1 };
    let mut window = Vec::with_capacity(K);
    for start in 0..starts {
        window.clear();
        for ii in 0..K {
            window.push(seq[(start + ii) % len]);
        }
        if let Some(code) = pack(&window) {
            set.insert(canonical(code));
        }
    }
    set
}

/// A k-mer index of a reference (k-mers through the origin included for circular ones).
pub struct KmerIndex {
    map: HashMap<u64, Vec<u32>>,
    /// Distinct canonical k-mers, for the strand-independent similarity score.
    canonical: HashSet<u64>,
    len: usize,
    circular: bool,
}

impl KmerIndex {
    pub fn new(seq: &[u8], circular: bool) -> KmerIndex {
        let len = seq.len();
        let mut map: HashMap<u64, Vec<u32>> = HashMap::new();
        if len >= K {
            let starts = if circular { len } else { len - K + 1 };
            let mut window = Vec::with_capacity(K);
            for start in 0..starts {
                window.clear();
                for ii in 0..K {
                    window.push(seq[(start + ii) % len]);
                }
                if let Some(code) = pack(&window) {
                    map.entry(code).or_default().push(start as u32);
                }
            }
        }
        KmerIndex { map, canonical: canonical_set(seq, circular), len, circular }
    }

    /// Strand-independent similarity between a read and this reference: shared distinct
    /// k-mers divided by the larger of the two k-mer sets. Unlike containment, extra sequence
    /// on either side lowers the score, so a construct that merely contains the read (or that
    /// the read merely contains) ranks below the true match.
    pub fn similarity(&self, read: &[u8]) -> f64 {
        let read_set = canonical_set(read, false);
        let shared = read_set.iter().filter(|k| self.canonical.contains(k)).count();
        let larger = read_set.len().max(self.canonical.len());
        if larger == 0 { 0.0 } else { shared as f64 / larger as f64 }
    }

    /// Place a read: best orientation and diagonal.
    pub fn place(&self, read: &[u8]) -> Option<Placement> {
        if read.len() < K || self.len < K {
            return None;
        }
        let mut best: Option<Placement> = None;
        for reverse in [false, true] {
            let probe = if reverse { alphabet::revcomp(read) } else { alphabet::upper(read) };
            let mut votes: HashMap<i64, u32> = HashMap::new();
            let mut found = 0usize;
            let mut total = 0usize;
            for (q, kmer) in probe.windows(K).enumerate() {
                let Some(code) = pack(kmer) else { continue };
                total += 1;
                let Some(hits) = self.map.get(&code) else { continue };
                found += 1;
                // Repetitive k-mers say little about placement.
                if hits.len() > 8 {
                    continue;
                }
                for &p in hits {
                    let mut diagonal = p as i64 - q as i64;
                    if self.circular {
                        diagonal = diagonal.rem_euclid(self.len as i64);
                    }
                    *votes.entry(diagonal).or_default() += 1;
                }
            }
            let Some((&diagonal, &count)) = votes.iter().max_by_key(|(d, c)| (**c, -(**d))) else { continue };
            let containment = if total == 0 { 0.0 } else { found as f64 / total as f64 };
            let candidate = Placement { reverse, diagonal, containment, votes: count };
            if best.is_none_or(|b| (candidate.votes, (candidate.containment * 1e6) as u64) > (b.votes, (b.containment * 1e6) as u64)) {
                best = Some(candidate);
            }
        }
        best
    }
}

/// One aligned column.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Column {
    /// Reference base i aligned to read base j.
    Pair(usize, usize),
    /// Reference base i has no partner (a deletion in the read).
    Deletion(usize),
    /// Read base j has no partner (an insertion in the read), after reference base `after`.
    Insertion { after: Option<usize>, read: usize },
}

/// Changes at least this long are structural (whole elements gained, lost or turned round) and
/// are reported as events rather than as base changes. 50 bp is the usual cut-off between small
/// indels and structural variants in sequencing analysis.
pub const STRUCTURAL_MIN: usize = 50;
/// The largest dynamic-programming matrix filled for one stretch, in cells (about 50 MB).
const MAX_CELLS: usize = 4_000_000;
/// The most runs of exact matches considered when chaining (the longest are kept).
const MAX_RUNS: usize = 4_000;
/// How far past the first or last anchor the read is still aligned base by base.
const END_REACH: usize = 2_000;
/// A stretch between anchors, long on both sides, that aligns with less identity than this is
/// a replacement (or an inversion), not a diverged copy of the same sequence.
const REPLACED_BELOW: f64 = 0.6;

/// The alignment method and its settings in a few sentences, as the app and the command line
/// state it.
pub fn method_summary() -> String {
    format!(
        "Seed, chain and align, the method of minimap2 and of SnapGene's Align to Reference. \
         Shared {K}-base words place the sequence on the map, on either strand and through the origin. \
         Runs of exact matches are chained with minimap2's gap cost, so a chain carries on across an insertion or deletion of any size. \
         The bases between runs are aligned by affine-gap dynamic programming (Gotoh): match +{MATCH}, mismatch {MISMATCH}, \
         a gap {GAP_OPEN} to open and {GAP_EXTEND} per base, with free end gaps so a partial read fits. \
         Changes of {STRUCTURAL_MIN} bp or more are structural changes, and identity is measured outside them."
    )
}

/// The alignment method and its scores in one line.
pub fn method_line() -> String {
    format!(
        "seed, chain and align ({K}-base seeds, minimap2-style chaining, Gotoh affine-gap alignment: match +{MATCH}, mismatch {MISMATCH}, gap {GAP_OPEN} to open and {GAP_EXTEND} per base)"
    )
}

/// The alignment method and its settings, for JSON output.
pub fn method_json() -> serde_json::Value {
    serde_json::json!({
        "name": "seed, chain and align",
        "seed_length": K,
        "chaining": "minimap2 gap cost, 0.01 x seed length x indel + 0.5 x log2(indel)",
        "alignment": "affine-gap dynamic programming (Gotoh) with free end gaps",
        "scores": { "match": MATCH, "mismatch": MISMATCH, "gap_open": GAP_OPEN, "gap_extend": GAP_EXTEND },
        "structural_min_bp": STRUCTURAL_MIN,
    })
}

fn score(a: u8, b: u8) -> i32 {
    if alphabet::iupac_matches(a, b) || alphabet::iupac_matches(b, a) { MATCH } else { MISMATCH }
}

/// How the ends of a stretch may be left unaligned.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Ends {
    /// Both stretches are aligned from their first base to their last.
    Global,
    /// Leading bases of either stretch may stay unaligned (the read before its first anchor).
    FreeStart,
    /// Trailing bases of either stretch may stay unaligned (the read after its last anchor).
    FreeEnd,
}

/// Affine-gap alignment (Gotoh) of a reference stretch with a read stretch. Column indices are
/// local to the two stretches. Returns `None` when the matrix would be too large.
fn align_stretch(reference: &[u8], read: &[u8], ends: Ends) -> Option<Vec<Column>> {
    let (n, m) = (reference.len(), read.len());
    if (n + 1).saturating_mul(m + 1) > MAX_CELLS {
        return None;
    }
    if n == 0 || m == 0 {
        // Nothing to pair: a global alignment keeps the bases as gaps, free ends drop them.
        if ends != Ends::Global {
            return Some(Vec::new());
        }
        return Some(if n == 0 { (0..m).map(|j| Column::Insertion { after: None, read: j }).collect() } else { (0..n).map(Column::Deletion).collect() });
    }
    let neg = i32::MIN / 4;
    let width = m + 1;
    let at = |i: usize, j: usize| i * width + j;
    // H: best score ending at (i, j); E: ending with a deletion (a reference base unpaired);
    // F: ending with an insertion (a read base unpaired). Traceback bits per cell: bits 0-1
    // where H came from (0 diagonal, 1 E, 2 F, 3 start), bit 2 E extended, bit 3 F extended.
    let cells = (n + 1) * width;
    let mut h = vec![neg; cells];
    let mut e = vec![neg; cells];
    let mut f = vec![neg; cells];
    let mut tb = vec![0u8; cells];
    for i in 0..=n {
        for j in 0..=m {
            let c = at(i, j);
            if (i == 0 && j == 0) || (ends == Ends::FreeStart && (i == 0 || j == 0)) {
                h[c] = 0;
                tb[c] = 3;
                continue;
            }
            if i > 0 {
                let up = at(i - 1, j);
                let (open, extend) = (h[up] + GAP_OPEN + GAP_EXTEND, e[up] + GAP_EXTEND);
                if extend > open {
                    e[c] = extend;
                    tb[c] |= 4;
                } else {
                    e[c] = open;
                }
            }
            if j > 0 {
                let left = at(i, j - 1);
                let (open, extend) = (h[left] + GAP_OPEN + GAP_EXTEND, f[left] + GAP_EXTEND);
                if extend > open {
                    f[c] = extend;
                    tb[c] |= 8;
                } else {
                    f[c] = open;
                }
            }
            let (mut best, mut from) = if i > 0 && j > 0 { (h[at(i - 1, j - 1)] + score(reference[i - 1], read[j - 1]), 0u8) } else { (neg, 3u8) };
            if e[c] > best {
                best = e[c];
                from = 1;
            }
            if f[c] > best {
                best = f[c];
                from = 2;
            }
            h[c] = best;
            tb[c] = (tb[c] & !3) | from;
        }
    }
    // Where the alignment ends: the far corner, or with free trailing ends the best cell of the
    // last row or column.
    let end = if ends == Ends::FreeEnd {
        let last_row = (0..=m).map(|j| (n, j));
        let last_col = (0..=n).map(|i| (i, m));
        last_row.chain(last_col).max_by_key(|&(i, j)| h[at(i, j)]).unwrap_or((n, m))
    } else {
        (n, m)
    };
    // Trace back, switching between the three matrices.
    let (mut i, mut j) = end;
    let mut state = 0u8; // 0 = H, 1 = E, 2 = F
    let mut columns = Vec::new();
    while i > 0 || j > 0 {
        let c = at(i, j);
        match state {
            0 => match tb[c] & 3 {
                0 => {
                    columns.push(Column::Pair(i - 1, j - 1));
                    i -= 1;
                    j -= 1;
                }
                1 => state = 1,
                2 => state = 2,
                _ => break,
            },
            1 => {
                columns.push(Column::Deletion(i - 1));
                let extended = tb[c] & 4 != 0;
                i -= 1;
                if !extended {
                    state = 0;
                }
            }
            _ => {
                columns.push(Column::Insertion { after: i.checked_sub(1), read: j - 1 });
                let extended = tb[c] & 8 != 0;
                j -= 1;
                if !extended {
                    state = 0;
                }
            }
        }
    }
    columns.reverse();
    Some(columns)
}

/// Columns of a stretch moved from stretch-local to window and read coordinates.
fn shifted(columns: Vec<Column>, p0: usize, q0: usize) -> impl Iterator<Item = Column> {
    columns.into_iter().map(move |c| match c {
        Column::Pair(i, j) => Column::Pair(i + p0, j + q0),
        Column::Deletion(i) => Column::Deletion(i + p0),
        Column::Insertion { after, read } => Column::Insertion { after: after.map(|a| a + p0).or(p0.checked_sub(1)), read: read + q0 },
    })
}

/// Fraction of columns that pair matching bases.
fn identity_of(columns: &[Column], reference: &[u8], read: &[u8]) -> f64 {
    let matched = columns.iter().filter(|c| matches!(c, Column::Pair(i, j) if score(reference[*i], read[*j]) == MATCH)).count();
    if columns.is_empty() { 0.0 } else { matched as f64 / columns.len() as f64 }
}

/// A reference stretch replaced wholesale: every reference base deleted, every read base inserted.
fn replaced(p0: usize, p1: usize, q0: usize, q1: usize) -> impl Iterator<Item = Column> {
    let after = p1.checked_sub(1);
    (p0..p1).map(Column::Deletion).chain((q0..q1).map(move |q| Column::Insertion { after, read: q }))
}

/// An exact k-mer match: read position `q` (in the oriented read) and window position `p`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Anchor {
    q: usize,
    p: usize,
}

/// A run of exact matches on one diagonal: read start `q`, window start `p`, length.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Run {
    q: usize,
    p: usize,
    len: usize,
}

/// Anchors merged into runs: anchors on the same diagonal whose k-mers touch or overlap form one
/// stretch of exact matches. Sorted by read position.
fn runs_of(anchors: &[Anchor]) -> Vec<Run> {
    let mut by_diagonal: HashMap<i64, Vec<usize>> = HashMap::new();
    for a in anchors {
        by_diagonal.entry(a.p as i64 - a.q as i64).or_default().push(a.q);
    }
    let mut runs = Vec::new();
    for (diagonal, mut starts) in by_diagonal {
        starts.sort_unstable();
        let mut current: Option<Run> = None;
        for q in starts {
            match current.as_mut() {
                Some(run) if q <= run.q + run.len => run.len = (q + K).max(run.q + run.len) - run.q,
                _ => {
                    runs.extend(current.take());
                    current = Some(Run { q, p: (q as i64 + diagonal) as usize, len: K });
                }
            }
        }
        runs.extend(current);
    }
    // A pathological read can produce very many short runs; the longest carry the alignment.
    if runs.len() > MAX_RUNS {
        runs.sort_by_key(|r| std::cmp::Reverse(r.len));
        runs.truncate(MAX_RUNS);
    }
    runs.sort_by_key(|r| (r.q, r.p));
    runs
}

/// The best chain of colinear runs, by the dynamic programme minimap2 uses for its anchors: each
/// run extends the best earlier run that starts before it in both sequences, gaining the bases
/// it adds (less any overlap) and paying for the difference in distance between the two, which is
/// an indel. The cost grows slowly with the indel's size, so a chain carries on across a deletion
/// or insertion of any length when well-matched sequence follows. Every earlier run is a
/// candidate, so an insertion that itself holds matches (a duplicated stretch) cannot hide the run
/// before it. Overlaps between chained runs are trimmed off the later run.
fn chain_runs(runs: &[Run]) -> Vec<Run> {
    if runs.is_empty() {
        return Vec::new();
    }
    let k = K as f64;
    let mut score = vec![0f64; runs.len()];
    let mut previous = vec![usize::MAX; runs.len()];
    for ii in 0..runs.len() {
        let a = runs[ii];
        score[ii] = a.len as f64;
        // iii = an earlier run that this one might follow.
        for iii in 0..ii {
            let b = runs[iii];
            if b.q >= a.q || b.p >= a.p {
                continue;
            }
            let overlap = (b.q + b.len).saturating_sub(a.q).max((b.p + b.len).saturating_sub(a.p));
            if overlap >= a.len {
                continue;
            }
            let dq = a.q + overlap - (b.q + b.len);
            let dp = a.p + overlap - (b.p + b.len);
            let indel = dq.abs_diff(dp) as f64;
            let cost = if indel == 0.0 { 0.0 } else { 0.01 * k * indel + 0.5 * indel.log2() };
            let candidate = score[iii] + (a.len - overlap) as f64 - cost;
            if candidate > score[ii] {
                score[ii] = candidate;
                previous[ii] = iii;
            }
        }
    }
    let mut at = (0..runs.len()).max_by(|&x, &y| score[x].total_cmp(&score[y])).unwrap_or(0);
    let mut chain = vec![runs[at]];
    while previous[at] != usize::MAX {
        at = previous[at];
        chain.push(runs[at]);
    }
    chain.reverse();
    let mut out: Vec<Run> = Vec::new();
    for run in chain {
        let shift = match out.last() {
            Some(b) => (b.q + b.len).saturating_sub(run.q).max((b.p + b.len).saturating_sub(run.p)),
            None => 0,
        };
        if shift < run.len {
            out.push(Run { q: run.q + shift, p: run.p + shift, len: run.len - shift });
        }
    }
    out
}

/// One alignment of a read, as given, to the reference.
struct Pass {
    reverse: bool,
    /// The read on the reference's top strand, upper case.
    probe: Vec<u8>,
    /// The reference stretch aligned to: the whole circle plus a margin either side for circular
    /// references (starting at reference position `window_start`), the whole sequence otherwise.
    window: Vec<u8>,
    window_start: usize,
    /// Aligned columns, without unaligned ends.
    columns: Vec<Column>,
    containment: f64,
    /// Read start and length of the longest exact run.
    longest_run: (usize, usize),
    /// Where the read as given starts in the read this pass aligned: a whole-plasmid read is
    /// opened again at this base for the second pass (0 otherwise).
    read_offset: usize,
}

/// Seed, chain and align a read in one pass.
fn align_pass(reference: &Document, index: &KmerIndex, read: &[u8]) -> Option<Pass> {
    let placement = index.place(read)?;
    let len = reference.len();
    let probe = alphabet::upper(&if placement.reverse { alphabet::revcomp(read) } else { read.to_vec() });
    let m = probe.len();
    // The window is the whole circle for circular references, starting a margin before the read.
    // It is longer by however much the read is shorter or longer than the reference: an indel
    // before the stretch that won the vote shifts that stretch's diagonal away from where the read
    // begins, and a read longer than the plasmid (a dimer) needs a second turn.
    let margin = 40 + m / 20;
    let (window_start, window_len) = if reference.is_circular() {
        let (shorter, longer) = (len.saturating_sub(m), m.saturating_sub(len));
        ((placement.diagonal - (margin + shorter) as i64).rem_euclid(len as i64) as usize, len + 2 * margin + shorter + longer)
    } else {
        (0, len)
    };
    let window: Vec<u8> = (0..window_len).map(|ii| reference.seq[(window_start + ii) % len].to_ascii_uppercase()).collect();

    // Anchors: read k-mers found in the window. K-mers found many times say little about where
    // the read lies and are skipped.
    let mut table: HashMap<u64, Vec<u32>> = HashMap::new();
    for (p, kmer) in window.windows(K).enumerate() {
        if let Some(code) = pack(kmer) {
            table.entry(code).or_default().push(p as u32);
        }
    }
    let mut anchors = Vec::new();
    for (q, kmer) in probe.windows(K).enumerate() {
        if let Some(hits) = pack(kmer).and_then(|code| table.get(&code))
            && hits.len() <= 8
        {
            anchors.extend(hits.iter().map(|&p| Anchor { q, p: p as usize }));
        }
    }
    let runs = chain_runs(&runs_of(&anchors));
    let (&first, &last) = (runs.first()?, runs.last()?);

    // Columns: the read before the first run, each run, the stretch between each pair of runs,
    // and the read after the last run.
    let mut columns: Vec<Column> = Vec::new();
    let (q0, p0) = (first.q, first.p);
    let lead_q = q0.saturating_sub(END_REACH);
    let lead_p = p0.saturating_sub(q0 - lead_q + 50);
    if let Some(c) = align_stretch(&window[lead_p..p0], &probe[lead_q..q0], Ends::FreeStart) {
        columns.extend(shifted(c, lead_p, lead_q));
    }
    for (ii, &Run { q, p, len: run_len }) in runs.iter().enumerate() {
        if ii > 0 {
            let before = runs[ii - 1];
            let (p_from, q_from) = (before.p + before.len, before.q + before.len);
            let (r, s) = (&window[p_from..p], &probe[q_from..q]);
            let aligned = align_stretch(r, s, Ends::Global);
            // Long on both sides and poorly matched: a replaced stretch, kept whole so it can be
            // reported (and recognised as an inversion) rather than as hundreds of substitutions.
            let replacement = r.len() >= STRUCTURAL_MIN && s.len() >= STRUCTURAL_MIN && aligned.as_ref().is_none_or(|c| identity_of(c, r, s) < REPLACED_BELOW);
            match aligned {
                Some(c) if !replacement => columns.extend(shifted(c, p_from, q_from)),
                _ => columns.extend(replaced(p_from, p, q_from, q)),
            }
        }
        columns.extend((0..run_len).map(|k| Column::Pair(p + k, q + k)));
    }
    let (tail_q, tail_p) = (last.q + last.len, last.p + last.len);
    let end_q = (tail_q + END_REACH).min(m);
    let end_p = (tail_p + (end_q - tail_q) + 50).min(window.len());
    if let Some(c) = align_stretch(&window[tail_p..end_p], &probe[tail_q..end_q], Ends::FreeEnd) {
        columns.extend(shifted(c, tail_p, tail_q));
    }
    // Gaps at the very ends are overhang, not differences.
    let lead = columns.iter().take_while(|c| !matches!(c, Column::Pair(..))).count();
    columns.drain(..lead);
    while matches!(columns.last(), Some(Column::Deletion(_)) | Some(Column::Insertion { .. })) {
        columns.pop();
    }
    let longest = runs.iter().max_by_key(|r| r.len).map(|r| (r.q, r.len)).unwrap_or((0, 0));
    Some(Pass { reverse: placement.reverse, probe, window, window_start, columns, containment: placement.containment, longest_run: longest, read_offset: 0 })
}

/// Align a read (bases, either orientation) to a reference document. `read_circular` says the
/// read is a whole circular molecule, such as a whole-plasmid sequencing consensus. A read on a
/// circular reference that is nearly as long as it is treated as circular too.
pub fn align_read(reference: &Document, read: &[u8], read_circular: bool) -> Option<Alignment> {
    let index = KmerIndex::new(&reference.seq, reference.is_circular());
    align_read_with(reference, &index, read, read_circular)
}

/// Like [`align_read`] with a prebuilt index (for many reads against one reference).
pub fn align_read_with(reference: &Document, index: &KmerIndex, read: &[u8], read_circular: bool) -> Option<Alignment> {
    let len = reference.len();
    if read.is_empty() || len == 0 {
        return None;
    }
    let whole = reference.is_circular() && (read_circular || read.len() * 10 >= len * 9);
    let first = align_pass(reference, index, read)?;
    if !whole {
        return Some(finish(reference, index, first, false));
    }
    // A whole-plasmid read is a circle opened wherever sequencing started. Open it again in the
    // middle of its longest exact match, so no change is split between the read's two ends.
    let m = read.len();
    let middle = first.longest_run.0 + first.longest_run.1 / 2;
    let cut = if first.reverse { (m - middle.min(m)) % m } else { middle % m };
    let rotated: Vec<u8> = read[cut..].iter().chain(&read[..cut]).copied().collect();
    let second = align_pass(reference, index, &rotated).map(|pass| Pass { read_offset: cut, ..pass }).unwrap_or(first);
    Some(finish(reference, index, second, true))
}

/// One raw difference with the alignment columns it spans.
struct Raw {
    difference: Difference,
    columns: Range<usize>,
}

/// Turn a pass into an alignment: differences, structural events and the summary numbers.
fn finish(reference: &Document, index: &KmerIndex, pass: Pass, whole: bool) -> Alignment {
    let len = reference.len();
    let to_ref = |i: usize| (pass.window_start + i) % len;
    let window = &pass.window;
    let probe = &pass.probe;
    let columns = &pass.columns;
    let blank = |kind, position, reference: String, read: String| Difference {
        kind,
        position,
        reference,
        read,
        features: vec![],
        effects: vec![],
        homopolymer: false,
        columns: 0..0,
    };

    // Walk the columns: count matches and collect runs of differences.
    let mut matches = 0;
    let mut raw: Vec<Raw> = Vec::new();
    let mut ii = 0;
    while ii < columns.len() {
        let start = ii;
        match columns[ii] {
            Column::Pair(i, j) => {
                if score(window[i], probe[j]) == MATCH {
                    matches += 1;
                    ii += 1;
                    continue;
                }
                // A run of mismatches becomes one substitution.
                let (mut refs, mut reads) = (String::new(), String::new());
                while let Some(&Column::Pair(i2, j2)) = columns.get(ii) {
                    if score(window[i2], probe[j2]) == MATCH {
                        break;
                    }
                    refs.push(window[i2] as char);
                    reads.push(probe[j2] as char);
                    ii += 1;
                }
                raw.push(Raw { difference: blank(DiffKind::Substitution, to_ref(i), refs, reads), columns: start..ii });
            }
            Column::Deletion(i) => {
                let mut refs = String::new();
                while let Some(&Column::Deletion(i2)) = columns.get(ii) {
                    refs.push(window[i2] as char);
                    ii += 1;
                }
                raw.push(Raw { difference: blank(DiffKind::Deletion, to_ref(i), refs, String::new()), columns: start..ii });
            }
            Column::Insertion { after, .. } => {
                let mut reads = String::new();
                while let Some(&Column::Insertion { read, .. }) = columns.get(ii) {
                    reads.push(probe[read] as char);
                    ii += 1;
                }
                raw.push(Raw { difference: blank(DiffKind::Insertion, to_ref(after.unwrap_or(0)), String::new(), reads), columns: start..ii });
            }
        }
    }

    // Long indels are structural events. A long deletion followed at once by a long insertion is
    // one replacement, which may be an inversion. Each event remembers the columns it covers, so
    // a view of the alignment can show it in place of its bases.
    let mut small: Vec<Difference> = Vec::new();
    let mut events: Vec<StructuralEvent> = Vec::new();
    let mut event_columns = 0;
    let mut ii = 0;
    while ii < raw.len() {
        let d = &raw[ii].difference;
        let size = d.reference.len().max(d.read.len());
        if d.kind == DiffKind::Substitution || size < STRUCTURAL_MIN {
            small.push(Difference { columns: raw[ii].columns.clone(), ..d.clone() });
            ii += 1;
            continue;
        }
        event_columns += raw[ii].columns.len();
        let next =
            raw.get(ii + 1).filter(|n| d.kind == DiffKind::Deletion && n.difference.kind == DiffKind::Insertion && n.difference.read.len() >= STRUCTURAL_MIN);
        match next {
            Some(n) => {
                event_columns += n.columns.len();
                let mut event = replacement_event(reference, index, d, &n.difference);
                event.columns = Some(raw[ii].columns.start..n.columns.end);
                events.push(event);
                ii += 2;
            }
            None => {
                let mut event = indel_event(reference, index, d);
                event.columns = Some(raw[ii].columns.clone());
                events.push(event);
                ii += 1;
            }
        }
    }

    // Which reference bases the alignment accounts for (paired or deleted), and how much of the
    // read it uses.
    let mut seen = vec![false; len];
    for c in columns {
        if let Column::Pair(i, _) | Column::Deletion(i) = *c {
            seen[to_ref(i)] = true;
        }
    }
    let ref_covered = seen.iter().filter(|&&s| s).count();
    let read_aligned = columns.iter().filter(|c| !matches!(c, Column::Deletion(_))).count();
    let first_read = columns.iter().find_map(|c| if let Column::Pair(_, j) = c { Some(*j) } else { None }).unwrap_or(0);
    let last_read = columns.iter().rev().find_map(|c| if let Column::Pair(_, j) = c { Some(*j) } else { None }).unwrap_or(0);
    let first_ref = columns.iter().find_map(|c| if let Column::Pair(i, _) = c { Some(*i) } else { None }).unwrap_or(0);
    let last_ref = columns.iter().rev().find_map(|c| if let Column::Pair(i, _) = c { Some(*i) } else { None }).unwrap_or(0);

    // A whole-plasmid read accounts for every base on both sides. An alignment that goes round
    // the circle more than once means the read holds that much of the plasmid twice (a dimer
    // holds all of it twice). What is left where the read's two ends meet is a change too: read
    // bases the reference lacks, reference bases the read lacks, or both.
    if whole {
        let first_window = columns.iter().find_map(|c| if let Column::Pair(i, _) = c { Some(*i) } else { None }).unwrap_or(0);
        let last_window = columns.iter().rev().find_map(|c| if let Column::Pair(i, _) = c { Some(*i) } else { None }).unwrap_or(0);
        let span = last_window + 1 - first_window.min(last_window);
        if span >= len + STRUCTURAL_MIN {
            let twice = span - len;
            let description = if twice * 4 >= len * 3 {
                "the read holds the whole plasmid twice (a dimer)".to_string()
            } else {
                format!("{} bp of the plasmid appear twice in the read", thousands(twice))
            };
            let point = Span::new((to_ref(last_window) + 1) % len, (to_ref(last_window) + 1) % len);
            events.push(StructuralEvent {
                kind: EventKind::Duplication,
                reference: point,
                read_bases: twice,
                features: Vec::new(),
                description,
                columns: None,
            });
        }
        let leftover: Vec<u8> = probe[last_read + 1..].iter().chain(&probe[..first_read]).copied().collect();
        let missing = len - ref_covered.min(len);
        let junction = to_ref(last_ref);
        if leftover.len() >= STRUCTURAL_MIN || missing >= STRUCTURAL_MIN {
            let lost = Span::new((junction + 1) % len, (junction + 1) % len + missing).normalize(len);
            let lost_bases = lost.extract(&reference.seq);
            let deletion = blank(DiffKind::Deletion, lost.start, String::from_utf8_lossy(&lost_bases).into_owned(), String::new());
            let insertion = blank(DiffKind::Insertion, junction, String::new(), String::from_utf8_lossy(&leftover).into_owned());
            let mut event = match (missing >= STRUCTURAL_MIN, leftover.len() >= STRUCTURAL_MIN) {
                (true, true) => replacement_event(reference, index, &deletion, &insertion),
                (true, false) => indel_event(reference, index, &deletion),
                _ => indel_event(reference, index, &insertion),
            };
            // The change sits where the last column meets the first: no columns of its own, at the
            // end of the alignment.
            event.columns = Some(columns.len()..columns.len());
            events.push(event);
        }
    }

    // Open reading frames of the reference, for changes in coding sequence the map does not
    // annotate as CDS (payloads are often drawn as misc_feature). Only ORFs that overlap an
    // annotated coding element on the same strand (a CDS, signal peptide, mature peptide or
    // gene) count; chance ORFs in backbone sequence would only add noise.
    let circular = reference.is_circular();
    let coding_marks: Vec<(Span, Strand)> = reference
        .features
        .iter()
        .filter(|f| ["CDS", "sig_peptide", "mat_peptide", "gene"].iter().any(|k| f.kind.eq_ignore_ascii_case(k)))
        .map(|f| (f.span(len), f.strand))
        .collect();
    let orfs: Vec<Feature> = find_orfs(&reference.seq, circular, &OrfOptions { min_codons: 100, ..OrfOptions::default() })
        .into_iter()
        .filter(|o| coding_marks.iter().any(|(span, strand)| *strand == o.strand && span.overlaps(&o.span, len)))
        .map(|o| {
            let (a, b) = o.span.to_1based(len);
            let mut f = Feature::new(format!("ORF {a}..{b}"), "CDS", Span::new(o.span.start, o.span.end), o.strand);
            f.segments[0].translated = Some(true);
            f
        })
        .collect();
    for d in &mut small {
        annotate(reference, &orfs, d);
    }
    events.sort_by_key(|e| e.reference.start);

    // The columns for showing the alignment, in reference positions and in positions of the read
    // as given: the read's index in the pass is turned back for a reverse read, and moved back by
    // the base a whole-plasmid read was opened at for the second pass.
    let m = probe.len();
    let read_position = |j: usize| ((if pass.reverse { m - 1 - j } else { j }) + pass.read_offset) % m;
    let mut aligned: Vec<AlignedColumn> = columns
        .iter()
        .map(|c| match *c {
            Column::Pair(i, j) => AlignedColumn { reference: Some(to_ref(i)), read: Some(read_position(j)), bases: [window[i], probe[j]] },
            Column::Deletion(i) => AlignedColumn { reference: Some(to_ref(i)), read: None, bases: [window[i], b'-'] },
            Column::Insertion { read, .. } => AlignedColumn { reference: None, read: Some(read_position(read)), bases: [b'-', probe[read]] },
        })
        .collect();
    // A whole-plasmid read is shown from the reference's first base, and its differences are
    // listed in that order too.
    if whole {
        start_at_the_origin(&mut aligned, &mut events, &mut small);
        small.sort_by_key(|d| d.columns.start);
    }
    Alignment {
        reverse: pass.reverse,
        ref_start: to_ref(first_ref),
        ref_covered,
        read_aligned,
        read_len: probe.len(),
        reference_len: len,
        whole_read: whole,
        matches,
        columns: columns.len() - event_columns,
        differences: small,
        events,
        containment: pass.containment,
        aligned,
    }
}

/// Turn the columns of a whole-plasmid read so they start at the reference's first base, as the
/// map does, rather than wherever the read happened to be opened. The turn comes at the first
/// column whose reference position is lower than the one before (the origin), moved on past any
/// structural change it would split. The column ranges of the events and of the differences turn
/// with the columns.
fn start_at_the_origin(aligned: &mut [AlignedColumn], events: &mut [StructuralEvent], differences: &mut [Difference]) {
    let n = aligned.len();
    let mut previous = None;
    let mut turn = 0;
    for (k, column) in aligned.iter().enumerate() {
        if let Some(r) = column.reference {
            if previous.is_some_and(|p| r < p) {
                turn = k;
                break;
            }
            previous = Some(r);
        }
    }
    for event in events.iter() {
        if let Some(columns) = &event.columns
            && columns.start < turn
            && turn < columns.end
        {
            turn = columns.end;
        }
    }
    if turn == 0 || turn >= n {
        return;
    }
    aligned.rotate_left(turn);
    // A range keeps its length and moves with its first column. No range straddles the turn.
    let shifted = |columns: &Range<usize>| {
        let start = if columns.start >= turn { columns.start - turn } else { columns.start + n - turn };
        start..start + columns.len()
    };
    for event in events.iter_mut() {
        if let Some(columns) = event.columns.as_mut() {
            *columns = shifted(columns);
        }
    }
    for difference in differences.iter_mut() {
        difference.columns = shifted(&difference.columns);
    }
}

/// A structural change between the read and the reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EventKind {
    /// Reference bases the read lacks.
    Deletion,
    /// Read bases the reference lacks, not found elsewhere in the reference.
    Insertion,
    /// Read bases the reference lacks that copy part of the reference (in either orientation).
    Duplication,
    /// A reference stretch the read has in the opposite orientation.
    Inversion,
    /// A reference stretch the read has replaced by other sequence.
    Replacement,
}

impl EventKind {
    pub fn name(self) -> &'static str {
        match self {
            EventKind::Deletion => "deletion",
            EventKind::Insertion => "insertion",
            EventKind::Duplication => "duplication",
            EventKind::Inversion => "inversion",
            EventKind::Replacement => "replacement",
        }
    }
}

/// One structural change, in reference coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct StructuralEvent {
    pub kind: EventKind,
    /// The reference stretch involved (0-based, may run through the origin). For insertions and
    /// duplications it is empty, at the place the new bases go.
    pub reference: Span,
    /// Read bases in place of that stretch (inserted, duplicated or replacing bases).
    pub read_bases: usize,
    /// Reference features the change touches.
    pub features: Vec<String>,
    /// One line for people, such as "1,500 bp deletion of reference 2,001..3,500".
    pub description: String,
    /// The columns of [`Alignment::aligned`] that hold the change's bases. A change found where a
    /// whole-plasmid read's two ends meet has an empty range at that point. `None` for a second
    /// copy of the whole plasmid, which spans the alignment rather than sitting in it.
    pub columns: Option<Range<usize>>,
}

/// Reference features a span touches, with how: "removes X" when the span holds all of X,
/// "cuts X" when it holds part of it.
fn touched_features(reference: &Document, span: Span) -> (Vec<String>, String) {
    let len = reference.len();
    let mut names = Vec::new();
    let mut words = Vec::new();
    for feature in reference.features.iter().filter(|f| !f.kind.eq_ignore_ascii_case("source")) {
        let fs = feature.span(len);
        if span.is_empty() {
            // An insertion point inside a feature.
            if fs.contains(span.start, len) && fs.start != span.start {
                names.push(feature.name.clone());
                words.push(format!("inside {}", feature.name));
            }
        } else if fs.overlaps(&span, len) {
            names.push(feature.name.clone());
            let all = fs.pieces(len).iter().all(|&(a, b)| (a..b).all(|p| span.contains(p, len)));
            words.push(format!("{} {}", if all { "removes" } else { "cuts" }, feature.name));
        }
    }
    let text = if words.is_empty() { String::new() } else { format!(" ({})", words.join(", ")) };
    (names, text)
}

/// Where a stretch of read sequence comes from: a copy of part of the reference (in either
/// orientation), or nowhere in it.
fn origin_of(reference: &Document, index: &KmerIndex, bases: &[u8]) -> Option<(bool, Span)> {
    let len = reference.len();
    let placement = index.place(bases)?;
    if placement.containment < 0.5 {
        return None;
    }
    let start = placement.diagonal.rem_euclid(len as i64) as usize;
    Some((placement.reverse, Span::new(start, start + bases.len().min(len)).normalize(len)))
}

fn indel_event(reference: &Document, index: &KmerIndex, d: &Difference) -> StructuralEvent {
    let len = reference.len();
    match d.kind {
        DiffKind::Deletion => {
            let span = Span::new(d.position, d.position + d.reference.len()).normalize(len);
            let (features, words) = touched_features(reference, span);
            StructuralEvent {
                kind: EventKind::Deletion,
                reference: span,
                read_bases: 0,
                features,
                description: format!("{} bp deletion of reference {}{words}", thousands(d.reference.len()), span.format(len)),
                columns: None,
            }
        }
        _ => {
            let point = Span::new((d.position + 1) % len, (d.position + 1) % len);
            let (features, words) = touched_features(reference, point);
            let n = d.read.len();
            let (kind, what) = match origin_of(reference, index, d.read.as_bytes()) {
                Some((false, _)) if n * 4 >= len * 3 => (EventKind::Duplication, "a second copy of the plasmid (the read is a dimer)".to_string()),
                Some((reverse, from)) => {
                    (EventKind::Duplication, format!("{} bp {}copy of reference {}", thousands(n), if reverse { "inverted " } else { "" }, from.format(len)))
                }
                None => (EventKind::Insertion, format!("{} bp of sequence not in the reference", thousands(n))),
            };
            StructuralEvent {
                kind,
                reference: point,
                read_bases: n,
                features,
                description: format!("{what}, inserted after reference {}{words}", thousands(d.position + 1)),
                columns: None,
            }
        }
    }
}

/// A reference stretch with other read sequence in its place. When that sequence is the stretch
/// itself turned round, the change is an inversion.
fn replacement_event(reference: &Document, index: &KmerIndex, deletion: &Difference, insertion: &Difference) -> StructuralEvent {
    let len = reference.len();
    let span = Span::new(deletion.position, deletion.position + deletion.reference.len()).normalize(len);
    let (features, words) = touched_features(reference, span);
    let flipped = alphabet::revcomp(insertion.read.as_bytes());
    let removed = deletion.reference.as_bytes();
    let inverted = align_stretch(removed, &flipped, Ends::Global).is_some_and(|c| identity_of(&c, removed, &flipped) >= 0.9);
    if inverted {
        return StructuralEvent {
            kind: EventKind::Inversion,
            reference: span,
            read_bases: insertion.read.len(),
            features,
            description: format!("{} bp inversion of reference {}{words}", thousands(span.len()), span.format(len)),
            columns: None,
        };
    }
    let instead = match origin_of(reference, index, insertion.read.as_bytes()) {
        Some((reverse, from)) => format!("a {}copy of reference {}", if reverse { "reversed " } else { "" }, from.format(len)),
        None => "sequence not in the reference".to_string(),
    };
    StructuralEvent {
        kind: EventKind::Replacement,
        reference: span,
        read_bases: insertion.read.len(),
        features,
        description: format!(
            "reference {} ({} bp) replaced by {} bp of {instead}{words}",
            span.format(len),
            thousands(span.len()),
            thousands(insertion.read.len())
        ),
        columns: None,
    }
}

/// The kind of a difference between read and reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiffKind {
    Substitution,
    Insertion,
    Deletion,
}

impl DiffKind {
    pub fn name(self) -> &'static str {
        match self {
            DiffKind::Substitution => "substitution",
            DiffKind::Insertion => "insertion",
            DiffKind::Deletion => "deletion",
        }
    }
}

/// One difference, in reference coordinates (0-based).
#[derive(Clone, Debug, PartialEq)]
pub struct Difference {
    pub kind: DiffKind,
    /// First reference base involved (for insertions: the base after which the read has extra
    /// bases).
    pub position: usize,
    /// Reference bases (empty for insertions).
    pub reference: String,
    /// Read bases (empty for deletions), on the reference's top strand.
    pub read: String,
    /// Features the difference falls in.
    pub features: Vec<String>,
    /// Protein consequence inside coding features, e.g. "AmpR: L191V" or "lacZα: frameshift at codon 12".
    pub effects: Vec<String>,
    /// The indel sits in a run of 4 or more identical bases (often a sequencing artefact).
    pub homopolymer: bool,
    /// The columns of [`Alignment::aligned`] that hold the difference.
    pub columns: Range<usize>,
}

impl Difference {
    /// The change in a few characters: "T>A", "+GGTA" or "−A". Runs longer than six bases are
    /// counted instead, as in "+12 bp".
    pub fn short_change(&self) -> String {
        let bases = |text: &str| if text.len() <= 6 { text.to_string() } else { format!("{} bp", thousands(text.len())) };
        match self.kind {
            DiffKind::Substitution => format!("{}>{}", self.reference, self.read),
            DiffKind::Insertion => format!("+{}", bases(&self.read)),
            DiffKind::Deletion => format!("−{}", bases(&self.reference)),
        }
    }

    /// One line for people: what changed, where, in which features and to what effect.
    pub fn describe(&self) -> String {
        let mut text = format!("{} at {}: {}", self.kind.name(), thousands(self.position + 1), self.short_change());
        if !self.features.is_empty() {
            text.push_str(&format!(" in {}", self.features.join(", ")));
        }
        if !self.effects.is_empty() {
            text.push_str(&format!(" ({})", self.effects.join("; ")));
        }
        if self.homopolymer {
            text.push_str(", in a run of one repeated base (a common sequencing error)");
        }
        text
    }

    /// The reference bases the difference involves (0-based). An insertion's span is empty, at
    /// the base after which the extra bases go.
    pub fn span(&self, reference_len: usize) -> Span {
        let len = reference_len.max(1);
        match self.kind {
            DiffKind::Insertion => Span::new((self.position + 1) % len, (self.position + 1) % len),
            _ => Span::new(self.position, self.position + self.reference.len()).normalize(len),
        }
    }

    /// The kind of difference, as maps colour it.
    pub fn color_kind(&self) -> DifferenceKind {
        match self.kind {
            DiffKind::Substitution => DifferenceKind::Changed,
            DiffKind::Insertion => DifferenceKind::Inserted,
            DiffKind::Deletion => DifferenceKind::Missing,
        }
    }
}

impl StructuralEvent {
    /// The kind of difference, as maps colour it. Deletions and replacements take reference
    /// bases away, insertions and duplications add bases, and an inversion changes them.
    pub fn color_kind(&self) -> DifferenceKind {
        match self.kind {
            EventKind::Deletion | EventKind::Replacement => DifferenceKind::Missing,
            EventKind::Insertion | EventKind::Duplication => DifferenceKind::Inserted,
            EventKind::Inversion => DifferenceKind::Changed,
        }
    }
}

/// How much of each sequence an alignment pairs, and how much only one of them has.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Tally {
    /// Reference bases paired with read bases, alike or not.
    pub aligned: usize,
    /// Reference bases the read lacks: deleted, or replaced by other sequence.
    pub missing: usize,
    /// Reference bases the read has turned round.
    pub inverted: usize,
    /// Read bases the reference lacks: inserted, duplicated, or in place of replaced bases.
    pub extra: usize,
}

/// Result of aligning one read to one reference.
#[derive(Clone, Debug)]
pub struct Alignment {
    pub reverse: bool,
    /// Reference position (0-based) the alignment starts at.
    pub ref_start: usize,
    /// Reference bases the alignment accounts for (paired with the read, or deleted from it).
    pub ref_covered: usize,
    /// Read bases inside the alignment.
    pub read_aligned: usize,
    pub read_len: usize,
    pub reference_len: usize,
    /// The read was treated as a whole circular molecule.
    pub whole_read: bool,
    /// Matching columns, and all columns outside structural events.
    pub matches: usize,
    pub columns: usize,
    /// Substitutions and indels shorter than [`STRUCTURAL_MIN`].
    pub differences: Vec<Difference>,
    /// Deletions, insertions, duplications, inversions and replacements of [`STRUCTURAL_MIN`]
    /// bases or more.
    pub events: Vec<StructuralEvent>,
    pub containment: f64,
    /// The alignment column by column, from its first paired base to its last. A whole-plasmid
    /// read starts at the reference's first base.
    pub aligned: Vec<AlignedColumn>,
}

/// One column of an alignment: a reference base over a read base. Either may be missing, where
/// the read lacks a reference base (a deletion) or holds a base the reference lacks (an
/// insertion).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AlignedColumn {
    /// Reference position (0-based) of the reference base.
    pub reference: Option<usize>,
    /// Position (0-based) of the read base in the read as given. It counts down along the
    /// alignment for a read on the reference's bottom strand.
    pub read: Option<usize>,
    /// The reference base and the read base on the reference's top strand, in upper case, with
    /// `b'-'` for a missing base.
    pub bases: [u8; 2],
}

impl AlignedColumn {
    /// Both bases are present and agree. An ambiguity code agrees with the bases it stands for.
    pub fn matches(&self) -> bool {
        self.reference.is_some() && self.read.is_some() && score(self.bases[0], self.bases[1]) == MATCH
    }

    /// What kind of difference the column shows, or `None` where the bases agree.
    pub fn kind(&self) -> Option<DifferenceKind> {
        if self.matches() {
            None
        } else if self.reference.is_none() {
            Some(DifferenceKind::Inserted)
        } else if self.read.is_none() {
            Some(DifferenceKind::Missing)
        } else {
            Some(DifferenceKind::Changed)
        }
    }

    /// The mark between the two bases when the alignment is written out: `|` for the same base,
    /// `:` for an ambiguity code that agrees, and a space for a difference.
    pub fn mark(&self) -> char {
        if !self.matches() {
            ' '
        } else if self.bases[0] == self.bases[1] {
            '|'
        } else {
            ':'
        }
    }
}

/// A piece of an alignment laid out for reading.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AlignmentBlock {
    /// One line of bases: this range of [`Alignment::aligned`].
    Bases(Range<usize>),
    /// A structural change, by its index in [`Alignment::events`], shown as one line in place of
    /// its bases.
    Event(usize),
}

/// The first and last positions (0-based) of each sequence in a stretch of alignment columns.
/// A sequence has `None` where the stretch holds only gaps in it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LineEnds {
    pub reference: Option<(usize, usize)>,
    pub read: Option<(usize, usize)>,
}

impl Alignment {
    /// The alignment laid out in lines of `width` columns. Each structural change takes a block of
    /// its own in place of its bases. A second copy of the whole plasmid has no columns and comes
    /// last.
    pub fn blocks(&self, width: usize) -> Vec<AlignmentBlock> {
        let width = width.max(1);
        // Lines of bases covering the columns from `from` up to `to`.
        let lines = |blocks: &mut Vec<AlignmentBlock>, from: usize, to: usize| {
            let mut start = from;
            while start < to {
                let end = (start + width).min(to);
                blocks.push(AlignmentBlock::Bases(start..end));
                start = end;
            }
        };
        let mut changes: Vec<(Range<usize>, usize)> = self.events.iter().enumerate().filter_map(|(ii, e)| e.columns.clone().map(|c| (c, ii))).collect();
        changes.sort_by_key(|(columns, _)| columns.start);
        let mut blocks = Vec::new();
        let mut at = 0;
        for (columns, ii) in changes {
            lines(&mut blocks, at, columns.start.max(at));
            blocks.push(AlignmentBlock::Event(ii));
            at = at.max(columns.end);
        }
        lines(&mut blocks, at, self.aligned.len());
        blocks.extend(self.events.iter().enumerate().filter(|(_, e)| e.columns.is_none()).map(|(ii, _)| AlignmentBlock::Event(ii)));
        blocks
    }

    /// The first and last reference and read positions among a range of columns.
    pub fn ends(&self, columns: Range<usize>) -> LineEnds {
        let slice = &self.aligned[columns];
        let first_last = |pick: fn(&AlignedColumn) -> Option<usize>| Some((slice.iter().find_map(pick)?, slice.iter().rev().find_map(pick)?));
        LineEnds { reference: first_last(|c| c.reference), read: first_last(|c| c.read) }
    }

    /// The block that shows a reference position (0-based): the line that holds it, or the
    /// structural change that took it away (or that adds bases right there).
    pub fn block_at(&self, blocks: &[AlignmentBlock], position: usize) -> Option<usize> {
        blocks.iter().position(|block| match block {
            AlignmentBlock::Bases(columns) => self.aligned[columns.clone()].iter().any(|c| c.reference == Some(position)),
            AlignmentBlock::Event(ii) => {
                let span = self.events[*ii].reference;
                if span.is_empty() { span.start == position } else { span.contains(position, self.reference_len) }
            }
        })
    }

    /// The reference stretch the alignment covers: all of it for a whole-plasmid read, otherwise
    /// from its first reference base to its last (through the origin when the read runs through it).
    pub fn covered(&self) -> Span {
        let len = self.reference_len;
        if self.whole_read || len == 0 {
            return Span::new(0, len);
        }
        match self.ends(0..self.aligned.len()).reference {
            Some((first, last)) => Span::new(first, first + (last + len - first) % len + 1).normalize(len),
            None => Span::new(0, 0),
        }
    }

    /// The covered stretch of the reference cut by what the read has there: aligned stretches,
    /// stretches it has turned round (inversions), and stretches it lacks (deletions and replaced
    /// bases), in order from the covered stretch's start.
    pub fn segments(&self) -> Vec<TrackSegment> {
        let len = self.reference_len;
        let covered = self.covered();
        if len == 0 || covered.is_empty() {
            return Vec::new();
        }
        let room = covered.len();
        // The structural changes that take reference bases away or turn them round, measured
        // from the covered stretch's start. On a whole plasmid a change may run past the end of
        // the turn, and its tail then sits at the start.
        let mut cuts: Vec<(usize, usize, SegmentKind)> = Vec::new();
        for event in &self.events {
            let kind = match event.kind {
                EventKind::Deletion | EventKind::Replacement => SegmentKind::Missing,
                EventKind::Inversion => SegmentKind::Inverted,
                EventKind::Insertion | EventKind::Duplication => continue,
            };
            if event.reference.is_empty() {
                continue;
            }
            let from = (event.reference.start + len - covered.start) % len;
            let to = from + event.reference.len();
            if to <= room {
                cuts.push((from, to, kind));
            } else if room >= len {
                cuts.push((from, room, kind));
                cuts.push((0, to - room, kind));
            } else {
                cuts.push((from.min(room), room, kind));
            }
        }
        cuts.sort_by_key(|cut| cut.0);
        // Walk along the covered stretch: aligned between the cuts, and each cut as its kind.
        let piece = |from: usize, to: usize, kind: SegmentKind| TrackSegment { kind, span: Span::new(covered.start + from, covered.start + to).normalize(len) };
        let mut segments = Vec::new();
        let mut at = 0;
        for (from, to, kind) in cuts {
            let from = from.max(at);
            if to <= from {
                continue;
            }
            if from > at {
                segments.push(piece(at, from, SegmentKind::Aligned));
            }
            segments.push(piece(from, to, kind));
            at = to;
        }
        if at < room {
            segments.push(piece(at, room, SegmentKind::Aligned));
        }
        segments
    }

    /// How many reference bases the read pairs, lacks or has turned round, and how many bases
    /// only the read has.
    pub fn tally(&self) -> Tally {
        let aligned = self.aligned.iter().filter(|c| c.reference.is_some() && c.read.is_some()).count();
        let mut tally = Tally { aligned, ..Tally::default() };
        for d in &self.differences {
            match d.kind {
                DiffKind::Deletion => tally.missing += d.reference.len(),
                DiffKind::Insertion => tally.extra += d.read.len(),
                DiffKind::Substitution => {}
            }
        }
        for e in &self.events {
            match e.kind {
                EventKind::Deletion => tally.missing += e.reference.len(),
                EventKind::Replacement => {
                    tally.missing += e.reference.len();
                    tally.extra += e.read_bases;
                }
                EventKind::Inversion => tally.inverted += e.reference.len(),
                EventKind::Insertion | EventKind::Duplication => tally.extra += e.read_bases,
            }
        }
        tally
    }

    /// The alignment as a track for a map: its segments, and a mark for each difference and
    /// structural change, in order along the reference.
    pub fn track(&self, name: &str, emphasis: bool) -> AlignmentTrack {
        let len = self.reference_len.max(1);
        let mut marks = Vec::new();
        for event in &self.events {
            let span = event.reference;
            let label = match event.kind {
                EventKind::Deletion | EventKind::Replacement => format!("{} missing ({} bp)", span.format(len), thousands(span.len())),
                EventKind::Insertion => format!("+{} bp inserted", thousands(event.read_bases)),
                EventKind::Duplication => format!("+{} bp duplicated", thousands(event.read_bases)),
                EventKind::Inversion => format!("{} turned round", span.format(len)),
            };
            marks.push(DifferenceMark { kind: event.color_kind(), span, label, description: event.description.clone(), structural: true });
            // The sequence in place of replaced bases is extra sequence. Its mark sits at the
            // middle of the stretch it stands in for.
            if event.kind == EventKind::Replacement && event.read_bases > 0 {
                let middle = (span.start + span.len() / 2) % len;
                marks.push(DifferenceMark {
                    kind: DifferenceKind::Inserted,
                    span: Span::new(middle, middle),
                    label: format!("+{} bp other sequence", thousands(event.read_bases)),
                    description: event.description.clone(),
                    structural: true,
                });
            }
        }
        for d in &self.differences {
            marks.push(DifferenceMark {
                kind: d.color_kind(),
                span: d.span(len),
                label: format!("{} {}", thousands(d.position + 1), d.short_change()),
                description: d.describe(),
                structural: false,
            });
        }
        marks.sort_by_key(|m| m.span.start);
        AlignmentTrack { name: name.to_string(), covered: self.covered(), segments: self.segments(), marks, emphasis }
    }

    /// The two rows of the alignment, with `-` for each missing base: the reference, then the read
    /// on the reference's top strand.
    pub fn rows(&self) -> (String, String) {
        self.aligned.iter().map(|c| (c.bases[0] as char, c.bases[1] as char)).unzip()
    }

    /// The alignment written out for people, `width` columns to a line. The reference sits over
    /// the read, with the marks of [`AlignedColumn::mark`] between them and 1-based positions at
    /// both ends of each line. Each structural change takes one line, in brackets, in place of
    /// its bases.
    pub fn pairwise_text(&self, width: usize, reference_name: &str, read_name: &str) -> String {
        // Names are cut to 12 characters and positions padded, so every line of bases starts in
        // the same column.
        let short = |name: &str| -> String { name.chars().take(12).collect() };
        let (top, bottom) = (short(reference_name), short(read_name));
        let name_w = top.chars().count().max(bottom.chars().count());
        let number_w = thousands(self.reference_len.max(self.read_len)).len();
        let number = |p: Option<usize>| p.map(|p| thousands(p + 1)).unwrap_or_default();
        let mut out = String::new();
        for block in self.blocks(width) {
            match block {
                AlignmentBlock::Bases(columns) => {
                    let ends = self.ends(columns.clone());
                    let slice = &self.aligned[columns];
                    let upper: String = slice.iter().map(|c| c.bases[0] as char).collect();
                    let marks: String = slice.iter().map(|c| c.mark()).collect();
                    let lower: String = slice.iter().map(|c| c.bases[1] as char).collect();
                    let (first, last) = (number(ends.reference.map(|r| r.0)), number(ends.reference.map(|r| r.1)));
                    out.push_str(&format!("{top:<name_w$} {first:>number_w$} {upper} {last}\n"));
                    out.push_str(format!("{:<name_w$} {:>number_w$} {marks}", "", "").trim_end());
                    out.push('\n');
                    let (first, last) = (number(ends.read.map(|r| r.0)), number(ends.read.map(|r| r.1)));
                    out.push_str(&format!("{bottom:<name_w$} {first:>number_w$} {lower} {last}\n\n"));
                }
                AlignmentBlock::Event(ii) => {
                    out.push_str(&format!("{:<name_w$} {:>number_w$} [{}]\n\n", "", "", self.events[ii].description));
                }
            }
        }
        out
    }

    /// Matching columns as a percentage of aligned columns, outside structural events.
    pub fn identity(&self) -> f64 {
        if self.columns == 0 { 0.0 } else { self.matches as f64 * 100.0 / self.columns as f64 }
    }

    /// Reference bases accounted for, as a fraction of the reference.
    pub fn coverage(&self) -> f64 {
        self.ref_covered as f64 / self.reference_len.max(1) as f64
    }

    /// The alignment as JSON, the same from the command line and from scripts. Positions are
    /// 1-based; a stretch through the origin has its start after its end.
    pub fn to_json(&self) -> serde_json::Value {
        let len = self.reference_len;
        let round = |x: f64, places: i32| (x * 10f64.powi(places)).round() / 10f64.powi(places);
        let differences: Vec<serde_json::Value> = self
            .differences
            .iter()
            .map(|d| {
                serde_json::json!({ "type": d.kind.name(), "position": d.position + 1, "reference": d.reference, "read": d.read,
                    "features": d.features, "effects": d.effects, "homopolymer": d.homopolymer })
            })
            .collect();
        let events: Vec<serde_json::Value> = self
            .events
            .iter()
            .map(|e| {
                let mut v = serde_json::json!({ "type": e.kind.name(), "length": e.reference.len(), "read_bases": e.read_bases,
                    "features": e.features, "description": e.description });
                if e.reference.is_empty() {
                    v["after"] = serde_json::json!(e.reference.start);
                } else {
                    let (start, end) = e.reference.to_1based(len);
                    v["reference_start"] = serde_json::json!(start);
                    v["reference_end"] = serde_json::json!(end);
                }
                v
            })
            .collect();
        serde_json::json!({
            "verdict": self.verdict(),
            "identity": round(self.identity(), 2),
            "orientation": if self.reverse { "reverse" } else { "forward" },
            "reference_start": self.ref_start + 1,
            "reference_covered": self.ref_covered,
            "reference_length": len,
            "coverage": round(self.coverage() * 100.0, 2),
            "read_length": self.read_len,
            "read_aligned": self.read_aligned,
            "whole_read": self.whole_read,
            "kmer_containment": round(self.containment * 100.0, 2),
            "events": events,
            "differences": differences,
            "method": method_json(),
        })
    }

    /// A one-line summary for people.
    pub fn summary(&self) -> String {
        let outside = if self.events.is_empty() { String::new() } else { format!(" outside {} structural change(s)", self.events.len()) };
        format!(
            "{:.2}% identity{outside} over {} of {} bp ({} strand), {} difference(s): {}",
            self.identity(),
            thousands(self.ref_covered),
            thousands(self.reference_len),
            if self.reverse { "reverse" } else { "forward" },
            self.differences.len(),
            self.verdict()
        )
    }

    /// A short verdict for people and scripts, by identity outside structural events:
    /// 99.999% and no differences is a match, 99% a match to inspect, 95% a near match.
    pub fn verdict(&self) -> &'static str {
        let identity = self.identity();
        let real = self.differences.iter().filter(|d| !d.homopolymer).count();
        if identity < 95.0 || self.read_aligned * 2 < self.read_len {
            "not this plasmid"
        } else if !self.events.is_empty() {
            if identity >= 99.0 { "structural change" } else { "structural change in a sibling construct" }
        } else if identity >= 99.999 && self.differences.is_empty() {
            if self.coverage() > 0.98 { "exact match" } else { "exact over the region read" }
        } else if identity >= 99.0 {
            if real == 0 { "match (homopolymer indels only)" } else { "match with differences to inspect" }
        } else {
            "near match: probably a sibling construct"
        }
    }
}

/// Fill in features, protein effects and the homopolymer flag for a difference. When no
/// annotated coding feature covers it, the reference's long ORFs are used for the effect.
fn annotate(reference: &Document, orfs: &[Feature], d: &mut Difference) {
    let len = reference.len();
    let seq = &reference.seq;
    // Bases the difference touches; an insertion touches the bases on both sides of it.
    let touched: Vec<usize> = match d.kind {
        DiffKind::Insertion => vec![d.position, (d.position + 1) % len],
        _ => (0..d.reference.len().max(1)).map(|k| (d.position + k) % len).collect(),
    };
    for feature in reference.features.iter().filter(|f| !f.kind.eq_ignore_ascii_case("source")) {
        let span = feature.span(len);
        if touched.iter().all(|&p| span.contains(p, len)) || (d.kind != DiffKind::Insertion && touched.iter().any(|&p| span.contains(p, len))) {
            d.features.push(feature.name.clone());
            if feature.is_coding()
                && let Some(effect) = coding_effect(reference, feature, d)
            {
                d.effects.push(format!("{}: {effect}", feature.name));
            }
        }
    }
    if d.effects.is_empty() {
        for orf in orfs {
            let span = orf.span(len);
            if touched.iter().all(|&p| span.contains(p, len))
                && let Some(effect) = coding_effect(reference, orf, d)
            {
                d.effects.push(format!("{}: {effect}", orf.name));
            }
        }
    }
    // Homopolymer: an indel of one repeated base inside (or next to) a run of 4+ of that base.
    if d.kind != DiffKind::Substitution {
        let bases = if d.kind == DiffKind::Insertion { &d.read } else { &d.reference };
        if let Some(b) = bases.bytes().next()
            && bases.bytes().all(|x| x == b)
        {
            let mut run = 0;
            let mut p = (d.position + len - 3) % len;
            for _ in 0..10 {
                if seq[p].to_ascii_uppercase() == b {
                    run += 1;
                    if run >= 4 {
                        d.homopolymer = true;
                        break;
                    }
                } else {
                    run = 0;
                }
                p = (p + 1) % len;
            }
        }
    }
}

/// Positions of a coding feature's bases in reading order, and the genetic code it uses.
fn coding_positions(reference: &Document, feature: &Feature) -> (Vec<usize>, &'static GeneticCode) {
    let len = reference.len();
    let mut positions = Vec::new();
    for segment in feature.segments.iter().filter(|s| !s.gap) {
        for (a, b) in segment.span.pieces(len) {
            positions.extend(a..b);
        }
    }
    if feature.strand == Strand::Reverse {
        positions.reverse();
    }
    let codon_start = feature.qualifier_text("codon_start").and_then(|t| t.trim().parse::<usize>().ok()).unwrap_or(1).clamp(1, 3);
    positions.drain(..(codon_start - 1).min(positions.len()));
    let code =
        feature.qualifier_text("transl_table").and_then(|t| t.trim().parse::<u32>().ok()).and_then(GeneticCode::get).unwrap_or_else(GeneticCode::standard);
    (positions, code)
}

/// Protein-level consequence of a difference inside a coding feature.
fn coding_effect(reference: &Document, feature: &Feature, d: &Difference) -> Option<String> {
    let (positions, code) = coding_positions(reference, feature);
    let reverse = feature.strand == Strand::Reverse;
    let base_at = |p: usize| -> u8 {
        let b = reference.seq[p].to_ascii_uppercase();
        if reverse { alphabet::complement_base(b) } else { b }
    };
    let index_of = |p: usize| positions.iter().position(|&x| x == p);
    match d.kind {
        DiffKind::Substitution => {
            // Codons touched, with the read's bases substituted in.
            let mut changed: HashMap<usize, u8> = HashMap::new();
            for (k, b) in d.read.bytes().enumerate() {
                let p = (d.position + k) % reference.len();
                changed.insert(p, if reverse { alphabet::complement_base(b.to_ascii_uppercase()) } else { b.to_ascii_uppercase() });
            }
            let mut codons: Vec<usize> = changed.keys().filter_map(|&p| index_of(p)).map(|i| i / 3).collect();
            codons.sort();
            codons.dedup();
            let mut parts = Vec::new();
            for codon in codons {
                let at = codon * 3;
                if at + 3 > positions.len() {
                    continue;
                }
                let ref_codon: Vec<u8> = positions[at..at + 3].iter().map(|&p| base_at(p)).collect();
                let alt_codon: Vec<u8> = positions[at..at + 3].iter().map(|&p| changed.get(&p).copied().unwrap_or_else(|| base_at(p))).collect();
                let (a, b) = (code.translate_codon(&ref_codon), code.translate_codon(&alt_codon));
                parts.push(if a == b { format!("{}{} silent", a as char, codon + 1) } else { format!("{}{}{}", a as char, codon + 1, b as char) });
            }
            (!parts.is_empty()).then(|| parts.join(", "))
        }
        DiffKind::Insertion | DiffKind::Deletion => {
            let n = d.read.len().max(d.reference.len());
            let codon = index_of(d.position).map(|i| i / 3 + 1)?;
            Some(if n.is_multiple_of(3) {
                format!("in-frame {} of {} aa at codon {codon}", d.kind.name(), n / 3)
            } else {
                format!("frameshift at codon {codon}")
            })
        }
    }
}

/// References ranked by k-mer similarity to a read (fast; no alignment). Best first.
pub fn rank_references(read: &[u8], references: &[(String, KmerIndex)]) -> Vec<(usize, f64)> {
    let mut scores: Vec<(usize, f64)> = references.iter().enumerate().map(|(ii, (_, index))| (ii, index.similarity(read))).collect();
    scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    scores
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Feature, Topology};
    use crate::span::Span;

    fn random_dna(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    fn plasmid() -> Document {
        let mut seq = random_dna(3000, 5);
        // A CDS at 500..800 starting with ATG (forward).
        seq[500..503].copy_from_slice(b"ATG");
        let mut doc = Document::new_dna("ref", seq, Topology::Circular);
        let mut cds = Feature::new("geneX", "CDS", Span::new(500, 800), Strand::Forward);
        cds.segments[0].translated = Some(true);
        doc.features.push(cds);
        doc
    }

    #[test]
    fn rotated_reverse_read_with_a_substitution() {
        let doc = plasmid();
        // The read starts at base 1200, runs once round the circle, and is reverse-complemented.
        let mut read: Vec<u8> = (0..3000).map(|k| doc.seq[(1200 + k) % 3000]).collect();
        // Mutate reference base 650 (inside geneX): read index 3000 - 1200 + 650 = 2450.
        let original = read[2450];
        read[2450] = if original == b'A' { b'C' } else { b'A' };
        let read = alphabet::revcomp(&read);
        let aln = align_read(&doc, &read, false).unwrap();
        assert!(aln.reverse);
        assert_eq!(aln.differences.len(), 1, "{:?}", aln.differences);
        let d = &aln.differences[0];
        assert_eq!(d.kind, DiffKind::Substitution);
        assert_eq!(d.position, 650);
        assert_eq!(d.features, vec!["geneX".to_string()]);
        assert_eq!(d.effects.len(), 1, "codon 51 changes or stays silent: {:?}", d.effects);
        assert!(aln.identity() > 99.9);
    }

    #[test]
    fn indels_and_frameshift() {
        let doc = plasmid();
        let mut read = doc.seq.clone();
        read.remove(700); // single-base deletion inside the CDS
        read.splice(1500..1500, b"GATTACA".iter().copied()); // insertion outside the CDS
        let aln = align_read(&doc, &read, false).unwrap();
        let kinds: Vec<DiffKind> = aln.differences.iter().map(|d| d.kind).collect();
        assert!(kinds.contains(&DiffKind::Deletion), "{:?}", aln.differences);
        assert!(kinds.contains(&DiffKind::Insertion), "{:?}", aln.differences);
        let del = aln.differences.iter().find(|d| d.kind == DiffKind::Deletion).unwrap();
        assert!(del.effects.iter().any(|e| e.contains("frameshift")), "{:?}", del.effects);
        let ins = aln.differences.iter().find(|d| d.kind == DiffKind::Insertion).unwrap();
        // An insertion next to matching bases can be placed a base either way; its size is fixed.
        assert_eq!(ins.read.len(), 7);
    }

    #[test]
    fn partial_sanger_style_read() {
        let doc = plasmid();
        let read = doc.seq[900..1700].to_vec();
        let aln = align_read(&doc, &read, false).unwrap();
        assert!(aln.differences.is_empty(), "{:?}", aln.differences);
        assert_eq!(aln.ref_start, 900);
        assert_eq!(aln.ref_covered, 800);
    }

    #[test]
    fn similarity_prefers_the_exact_construct() {
        let base = random_dna(4000, 9);
        let mut superset = base.clone();
        superset.extend(random_dna(1500, 10)); // the same construct with an extra cassette
        let exact = KmerIndex::new(&base, true);
        let bigger = KmerIndex::new(&superset, true);
        let read = alphabet::revcomp(&base);
        assert!(exact.similarity(&read) > 0.99);
        assert!(bigger.similarity(&read) < exact.similarity(&read) - 0.2);
        assert_eq!(revcomp_code(pack(b"AACGTTTTTTTTTGC").unwrap(), K), pack(b"GCAAAAAAAAACGTT").unwrap());
    }

    /// A 6 kb circular plasmid with a promoter at 1900..1990 and a gene at 2100..2400 (0-based).
    fn big_plasmid() -> Document {
        let mut doc = Document::new_dna("ref", random_dna(6000, 21), Topology::Circular);
        doc.features.push(Feature::new("pX", "promoter", Span::new(1900, 1990), Strand::Forward));
        doc.features.push(Feature::new("geneA", "CDS", Span::new(2100, 2400), Strand::Forward));
        doc
    }

    /// The only structural event of an alignment, checked to start within a few bases of
    /// `start` (an indel between repeated bases can sit a base or two either way).
    fn single_event(aln: &Alignment, kind: EventKind, start: usize) -> &StructuralEvent {
        assert_eq!(aln.events.len(), 1, "{:#?}", aln.events);
        let event = &aln.events[0];
        assert_eq!(event.kind, kind, "{event:#?}");
        assert!(event.reference.start.abs_diff(start) <= 3, "{event:#?}");
        event
    }

    #[test]
    fn a_large_deletion_is_one_event_and_the_rest_still_matches() {
        let doc = big_plasmid();
        let read = [&doc.seq[..2000], &doc.seq[2300..]].concat();
        let aln = align_read(&doc, &read, true).unwrap();
        let event = single_event(&aln, EventKind::Deletion, 2000);
        assert_eq!(event.reference.len(), 300);
        assert!(event.description.contains("cuts geneA"), "{}", event.description);
        assert!(aln.differences.is_empty(), "{:?}", aln.differences);
        assert!(aln.identity() > 99.99);
        assert_eq!(aln.verdict(), "structural change");
        // Bigger still: 1,500 bp, which removes the gene outright.
        let read = [&doc.seq[..2000], &doc.seq[3500..]].concat();
        let aln = align_read(&doc, &read, true).unwrap();
        let event = single_event(&aln, EventKind::Deletion, 2000);
        assert_eq!(event.reference.len(), 1500);
        assert!(event.description.contains("removes geneA"), "{}", event.description);
    }

    #[test]
    fn insertions_are_told_apart_by_where_the_bases_come_from() {
        let doc = big_plasmid();
        // Foreign sequence.
        let read = [&doc.seq[..3000], &random_dna(800, 77)[..], &doc.seq[3000..]].concat();
        let aln = align_read(&doc, &read, true).unwrap();
        let event = single_event(&aln, EventKind::Insertion, 3000);
        assert_eq!(event.read_bases, 800);
        assert!(aln.identity() > 99.99);
        // A tandem copy of 400 bp of the plasmid itself.
        let read = [&doc.seq[..1400], &doc.seq[1000..1400], &doc.seq[1400..]].concat();
        let aln = align_read(&doc, &read, true).unwrap();
        let event = &aln.events[0];
        assert_eq!(aln.events.len(), 1, "{:#?}", aln.events);
        assert_eq!(event.kind, EventKind::Duplication, "{event:#?}");
        assert_eq!(event.read_bases, 400);
    }

    #[test]
    fn an_inverted_stretch_is_an_inversion() {
        let doc = big_plasmid();
        let read = [&doc.seq[..2000], &alphabet::revcomp(&doc.seq[2000..2500])[..], &doc.seq[2500..]].concat();
        let aln = align_read(&doc, &read, true).unwrap();
        let event = single_event(&aln, EventKind::Inversion, 2000);
        assert!(event.reference.len().abs_diff(500) <= 6, "{event:#?}");
        assert_eq!(aln.verdict(), "structural change");
    }

    #[test]
    fn a_dimer_is_a_second_copy() {
        let doc = big_plasmid();
        let read = [&doc.seq[..], &doc.seq[..]].concat();
        let aln = align_read(&doc, &read, true).unwrap();
        assert_eq!(aln.events.len(), 1, "{:#?}", aln.events);
        assert_eq!(aln.events[0].kind, EventKind::Duplication);
        assert!(aln.events[0].description.contains("dimer"), "{}", aln.events[0].description);
    }

    #[test]
    fn a_deletion_across_the_reads_start_is_still_found() {
        // The sequencer happened to start reading right where 200 bp are missing: the read's
        // two ends meet at the deletion.
        let doc = big_plasmid();
        let read = doc.seq[100..5900].to_vec();
        let aln = align_read(&doc, &read, true).unwrap();
        let event = single_event(&aln, EventKind::Deletion, 5900);
        assert_eq!(event.reference.len(), 200);
        assert!(event.reference.wraps(6000), "the deletion runs through the origin: {event:#?}");
    }

    #[test]
    fn a_relative_with_scattered_changes_is_a_near_match_not_a_structural_change() {
        let doc = big_plasmid();
        let mut read = doc.seq.clone();
        // ii = every 20th base changed: 5% different, spread out.
        for ii in (10..read.len()).step_by(20) {
            read[ii] = if read[ii] == b'A' { b'C' } else { b'A' };
        }
        let aln = align_read(&doc, &read, true).unwrap();
        assert!(aln.events.is_empty(), "{:#?}", aln.events);
        assert!((94.0..96.0).contains(&aln.identity()), "{}", aln.identity());
        assert_eq!(aln.verdict(), "near match: probably a sibling construct");
        // And something unrelated does not align at all.
        assert!(align_read(&doc, &random_dna(6000, 99), true).is_none());
    }

    #[test]
    fn the_columns_point_back_at_both_sequences() {
        // A whole-plasmid read, reverse-complemented and opened at base 1200, with one
        // substitution: the second pass opens it again, so both of its turns must be undone.
        let doc = plasmid();
        let mut read: Vec<u8> = (0..3000).map(|k| doc.seq[(1200 + k) % 3000]).collect();
        // Reference bases 100 and 650 changed (read indexes 1900 and 2450).
        for at in [1900, 2450] {
            read[at] = if read[at] == b'A' { b'C' } else { b'A' };
        }
        let read = alphabet::revcomp(&read);
        let aln = align_read(&doc, &read, true).unwrap();
        let listed: Vec<usize> = aln.differences.iter().map(|d| d.position).collect();
        assert_eq!(listed, vec![100, 650], "differences are listed from the reference's first base");
        assert_eq!(aln.aligned.len(), 3000);
        assert_eq!(aln.aligned[0].reference, Some(0), "a whole-plasmid read is shown from the reference's first base");
        let mut previous: Option<usize> = None;
        for column in &aln.aligned {
            let (r, q) = (column.reference.unwrap(), column.read.unwrap());
            assert_eq!(column.bases[0], doc.seq[r].to_ascii_uppercase());
            // The read is on the bottom strand, so its base is the complement of the one shown.
            assert_eq!(column.bases[1], alphabet::complement_base(read[q].to_ascii_uppercase()));
            if let Some(p) = previous {
                assert_eq!(r, (p + 1) % 3000, "reference positions run on through the origin");
            }
            previous = Some(r);
        }
        let differing: Vec<usize> = aln.aligned.iter().filter(|c| !c.matches()).map(|c| c.reference.unwrap()).collect();
        assert_eq!(differing, vec![100, 650]);
    }

    #[test]
    fn a_structural_change_takes_one_block_in_place_of_its_bases() {
        let doc = big_plasmid();
        let read = [&doc.seq[..2000], &doc.seq[3500..]].concat();
        let aln = align_read(&doc, &read, true).unwrap();
        let blocks = aln.blocks(60);
        let events: Vec<&AlignmentBlock> = blocks.iter().filter(|b| matches!(b, AlignmentBlock::Event(_))).collect();
        assert_eq!(events, vec![&AlignmentBlock::Event(0)]);
        // No line of bases shows a deleted base, and the lines and the event cover every column.
        let deleted = aln.events[0].reference;
        let mut covered = aln.events[0].columns.clone().unwrap().len();
        for block in &blocks {
            if let AlignmentBlock::Bases(columns) = block {
                assert!(columns.len() <= 60);
                covered += columns.len();
                assert!(aln.aligned[columns.clone()].iter().all(|c| c.reference.is_none_or(|r| !deleted.contains(r, 6000))));
            }
        }
        assert_eq!(covered, aln.aligned.len());
        // Any deleted base leads to the event's block.
        let event_block = blocks.iter().position(|b| *b == AlignmentBlock::Event(0));
        assert_eq!(aln.block_at(&blocks, 2700), event_block);
        // The written alignment has the change on a line of its own, in brackets.
        let text = aln.pairwise_text(60, "ref", "read");
        let line = text.lines().find(|l| l.contains("bp deletion")).unwrap_or_default().trim();
        assert_eq!(line, format!("[{}]", aln.events[0].description));
        assert!(line.starts_with("[1,500 bp deletion of reference") && line.ends_with("(removes geneA)]"), "{line}");
    }

    #[test]
    fn the_written_alignment_marks_each_difference() {
        let doc = plasmid();
        let mut read = doc.seq[900..1020].to_vec();
        read[5] = if read[5] == b'A' { b'C' } else { b'A' };
        let aln = align_read(&doc, &read, false).unwrap();
        let text = aln.pairwise_text(60, "ref", "read");
        let lines: Vec<&str> = text.lines().collect();
        // Reference line, marks, read line, blank: the first line of bases runs 901..960 on the
        // reference and 1..60 on the read, with a gap in the marks over the changed base.
        assert!(lines[0].starts_with("ref    901 "), "{text}");
        assert!(lines[0].ends_with(" 960"), "{text}");
        assert!(lines[2].starts_with("read     1 "), "{text}");
        assert!(lines[2].ends_with(" 60"), "{text}");
        let marks = &lines[1]["ref    901 ".len()..];
        assert_eq!(marks.find(' '), Some(5), "{text}");
        assert_eq!(marks.matches('|').count(), 59);
        let (upper, lower) = aln.rows();
        assert_eq!(upper.len(), 120);
        assert_eq!(lower, String::from_utf8(read.to_ascii_uppercase()).unwrap());
    }

    #[test]
    fn each_difference_owns_exactly_its_differing_columns() {
        // A whole-plasmid read with a substitution, a one-base deletion and a four-base insertion,
        // so the columns turn to the origin as well.
        let doc = big_plasmid();
        let mut read = doc.seq.clone();
        read[500] = if read[500] == b'A' { b'C' } else { b'A' };
        read.remove(3000);
        read.splice(4500..4500, b"GGTA".iter().copied());
        let aln = align_read(&doc, &read, true).unwrap();
        let mut owned = vec![false; aln.aligned.len()];
        for d in &aln.differences {
            assert!(!d.columns.is_empty(), "{d:?}");
            for k in d.columns.clone() {
                assert!(!aln.aligned[k].matches(), "column {k} of {d:?} matches");
                owned[k] = true;
            }
        }
        // Every differing column belongs to a difference.
        for (k, column) in aln.aligned.iter().enumerate() {
            assert_eq!(column.kind().is_some(), owned[k], "column {k}: {column:?}");
        }
        let kinds: Vec<Option<DifferenceKind>> = aln.differences.iter().map(|d| aln.aligned[d.columns.start].kind()).collect();
        assert_eq!(kinds, vec![Some(DifferenceKind::Changed), Some(DifferenceKind::Missing), Some(DifferenceKind::Inserted)]);
    }

    #[test]
    fn a_track_marks_every_difference_where_it_lies() {
        let doc = big_plasmid();
        // A whole read missing 1,500 bp covers the whole map and has one structural mark.
        let read = [&doc.seq[..2000], &doc.seq[3500..]].concat();
        let track = align_read(&doc, &read, true).unwrap().track("read", true);
        assert_eq!(track.covered, Span::new(0, 6000));
        assert_eq!(track.marks.len(), 1);
        let mark = &track.marks[0];
        assert_eq!((mark.kind, mark.structural, mark.span.len()), (DifferenceKind::Missing, true, 1500));
        // Its ring is aligned up to the missing stretch and after it.
        let kinds: Vec<(SegmentKind, Span)> = track.segments.iter().map(|s| (s.kind, s.span)).collect();
        assert_eq!(kinds.len(), 3, "{kinds:?}");
        assert_eq!(kinds[1].0, SegmentKind::Missing);
        assert_eq!(kinds[1].1.len(), 1500);
        assert!(kinds[0].0 == SegmentKind::Aligned && kinds[2].0 == SegmentKind::Aligned);
        // A Sanger read covers its own stretch, and a substitution is one changed base.
        let mut read = doc.seq[900..1700].to_vec();
        read[100] = if read[100] == b'A' { b'C' } else { b'A' };
        let track = align_read(&doc, &read, false).unwrap().track("sanger", true);
        assert_eq!(track.covered, Span::new(900, 1700));
        assert_eq!(track.segments, vec![TrackSegment { kind: SegmentKind::Aligned, span: Span::new(900, 1700) }]);
        assert_eq!(track.marks.len(), 1);
        assert_eq!((track.marks[0].kind, track.marks[0].span), (DifferenceKind::Changed, Span::new(1000, 1001)));
        assert!(track.marks[0].label.starts_with("1,001 "), "{}", track.marks[0].label);
    }

    #[test]
    fn a_plasmid_that_shares_only_its_backbone_is_mostly_missing() {
        // The other plasmid holds 1,500 bp of the reference (1,001..2,500) and 3,000 bp of its
        // own, so the rest of the reference (4,500 bp, through the origin) is replaced.
        let doc = big_plasmid();
        let read = [&doc.seq[1000..2500], &random_dna(3000, 77)[..]].concat();
        let aln = align_read(&doc, &read, true).unwrap();
        let tally = aln.tally();
        assert_eq!(tally.aligned, 1500, "{tally:?}");
        assert_eq!((tally.missing, tally.extra), (4500, 3000), "{tally:?}");
        // The ring is aligned over the shared stretch and missing everywhere else.
        let track = aln.track("other", true);
        let aligned: usize = track.segments.iter().filter(|s| s.kind == SegmentKind::Aligned).map(|s| s.span.len()).sum();
        let missing: usize = track.segments.iter().filter(|s| s.kind == SegmentKind::Missing).map(|s| s.span.len()).sum();
        assert_eq!((aligned, missing), (1500, 4500), "{:?}", track.segments);
        // A missing stretch in red, and a mark for the other plasmid's own 3,000 bp.
        assert!(track.marks.iter().any(|m| m.kind == DifferenceKind::Missing && m.label.contains("missing (4,500 bp)")), "{:?}", track.marks);
        assert!(track.marks.iter().any(|m| m.kind == DifferenceKind::Inserted && m.label == "+3,000 bp other sequence"), "{:?}", track.marks);
    }

    #[test]
    fn the_method_names_its_settings() {
        let text = method_summary();
        for part in ["Gotoh", "minimap2", "15-base", "match +2", "mismatch -3", "50 bp"] {
            assert!(text.contains(part), "{part} missing from: {text}");
        }
        assert_eq!(method_json()["scores"]["gap_open"], -5);
    }

    #[test]
    fn homopolymer_indel_is_flagged() {
        let mut doc = plasmid();
        doc.seq[1000..1006].copy_from_slice(b"AAAAAA");
        let mut read = doc.seq.clone();
        read.insert(1003, b'A');
        let aln = align_read(&doc, &read, false).unwrap();
        assert_eq!(aln.differences.len(), 1, "{:?}", aln.differences);
        assert!(aln.differences[0].homopolymer);
    }
}
