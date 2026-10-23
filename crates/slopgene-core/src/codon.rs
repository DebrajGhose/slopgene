//! Codon usage tables, reverse translation, codon optimisation and silent restriction-site edits.
//!
//! # Usage tables
//!
//! The built-in tables hold the raw codon counts of the Kazusa Codon Usage Database
//! (<https://www.kazusa.or.jp/codon/>), which tabulated them from NCBI-GenBank Flat File Release
//! 160.0: Nakamura Y, Gojobori T, Ikemura T, "Codon usage tabulated from international DNA
//! sequence databases: status for the year 2000", Nucleic Acids Res. 28:292 (2000). The counts
//! live in `data/codon_usage.tsv` with the organism, its NCBI taxonomy id, the number of coding
//! sequences and the page each row was read from. Frequencies are derived from the counts here,
//! never stored, so the numbers in the file are exactly what the source shows.
//!
//! # Derived measures
//!
//! * Frequency per thousand is a codon's count divided by the table's total, times 1000.
//! * Relative synonymous usage is a codon's share of its amino acid's codons under a given
//!   genetic code. [`CodonTable::rscu`] is the RSCU of Sharp PM, Tuohy TMF, Mosurski KR,
//!   "Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed
//!   genes", Nucleic Acids Res. 14:5125-5143 (1986): the share divided by the share an evenly
//!   used family would have, so 1 means average and 2 means twice as often as its synonyms.
//! * The codon adaptation index is the geometric mean of the relative adaptiveness of a coding
//!   sequence's codons, after Sharp PM, Li WH, "The codon adaptation index: a measure of
//!   directional synonymous codon usage bias, and its potential applications", Nucleic Acids
//!   Res. 15:1281-1295 (1987). Stop codons and the single-codon families (Met and Trp under the
//!   standard code) carry no information and are left out, as that paper does.
//!
//! # Edits
//!
//! [`reverse_translate`] writes a coding sequence for a protein, [`optimize`] rewrites one for a
//! host while keeping the protein, and [`silent_site_changes`] lists the silent substitutions
//! that add or remove a restriction site inside a coding sequence, the way a plasmid editor
//! offers them. Weighted choices use SplitMix64 (Steele GL, Lea D, Flood CH, "Fast splittable
//! pseudorandom number generators", OOPSLA 2014) from an explicit seed, so a run is reproducible
//! on every machine.

use crate::alphabet::{self, iupac_matches};
use crate::enzymes::{Enzyme, EnzymeDb};
use crate::error::{Error, Result, invalid};
use crate::model::Strand;
use crate::span::Span;
use crate::translate::GeneticCode;
use std::collections::HashMap;
use std::sync::OnceLock;

/// The embedded usage tables (see the header of the file for the sources).
const BUILTIN_TSV: &str = include_str!("../data/codon_usage.tsv");

/// Number of codons in the genetic code, and the length of every count array here.
pub const CODON_COUNT: usize = 64;

/// Index of a codon in NCBI order: the first, second and third base each run through T, C, A, G,
/// so the index is `16 * b1 + 4 * b2 + b3` with T = 0, C = 1, A = 2, G = 3. `U` counts as `T`,
/// case is ignored, and an ambiguous or non-nucleotide base gives `None`.
pub fn codon_index(codon: &[u8]) -> Option<usize> {
    if codon.len() != 3 {
        return None;
    }
    let rank = |base: u8| match base.to_ascii_uppercase() {
        b'T' | b'U' => Some(0usize),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    Some(16 * rank(codon[0])? + 4 * rank(codon[1])? + rank(codon[2])?)
}

/// The three bases of the codon with this index (see [`codon_index`]).
pub fn codon_bases(index: usize) -> [u8; 3] {
    const ORDER: [u8; 4] = *b"TCAG";
    [ORDER[(index / 16) % 4], ORDER[(index / 4) % 4], ORDER[index % 4]]
}

/// The codon with this index, as text.
pub fn codon_name(index: usize) -> String {
    String::from_utf8(codon_bases(index).to_vec()).unwrap_or_default()
}

/// Every codon that codes for `aa` under `code`, in NCBI order. `*` means stop.
pub fn synonymous_codons(aa: u8, code: &GeneticCode) -> Vec<usize> {
    let aa = aa.to_ascii_uppercase();
    (0..CODON_COUNT).filter(|&ii| code.amino_acids[ii] == aa).collect()
}

/// How often each codon is used in one organism.
///
/// The counts are the raw numbers from the source table, in NCBI order (see [`codon_index`]).
/// Frequencies are computed from them on demand.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CodonTable {
    /// Short name SlopGene uses, such as `ecoli` or `human`.
    pub key: String,
    /// Other names that resolve to this table, such as `e. coli` or `homo sapiens`.
    pub aliases: Vec<String>,
    /// Organism as the source names it, such as `Escherichia coli W3110`.
    pub organism: String,
    /// NCBI taxonomy id.
    pub taxid: u32,
    /// GenBank division the source tabulated, such as `gbbct`.
    pub division: String,
    /// Number of coding sequences the counts came from.
    pub cds: u64,
    /// Number of codons the source states.
    pub codons: u64,
    /// Page the row was read from.
    pub source: String,
    /// Count of each codon, in NCBI order.
    pub counts: [u64; CODON_COUNT],
}

impl CodonTable {
    /// Count of one codon, or `None` for an ambiguous or malformed codon.
    pub fn count(&self, codon: &[u8]) -> Option<u64> {
        codon_index(codon).map(|ii| self.counts[ii])
    }

    /// Total of every count. This equals [`CodonTable::codons`] for the built-in tables.
    pub fn total(&self) -> u64 {
        self.counts.iter().sum()
    }

    /// Frequency of a codon per thousand codons.
    pub fn per_thousand(&self, index: usize) -> f64 {
        let total = self.total();
        if total == 0 { 0.0 } else { self.counts[index] as f64 * 1000.0 / total as f64 }
    }

    /// Relative synonymous usage: a codon's share of the codons for its amino acid under `code`.
    /// The shares of one amino acid add up to 1, and an amino acid with no counts gives 0.
    pub fn relative_usage(&self, index: usize, code: &GeneticCode) -> f64 {
        let family: u64 = synonymous_codons(code.amino_acids[index], code).into_iter().map(|ii| self.counts[ii]).sum();
        if family == 0 { 0.0 } else { self.counts[index] as f64 / family as f64 }
    }

    /// Relative synonymous codon usage (RSCU): the share divided by the share an evenly used
    /// family would have, so 1 is average usage (Sharp, Tuohy and Mosurski 1986).
    pub fn rscu(&self, index: usize, code: &GeneticCode) -> f64 {
        let family = synonymous_codons(code.amino_acids[index], code);
        let total: u64 = family.iter().map(|&ii| self.counts[ii]).sum();
        if total == 0 { 0.0 } else { self.counts[index] as f64 * family.len() as f64 / total as f64 }
    }

    /// The most used codon for an amino acid under `code`, or `None` when the code has none.
    /// Ties go to the first codon in NCBI order, so the answer never depends on iteration order.
    pub fn best_codon_index(&self, aa: u8, code: &GeneticCode) -> Option<usize> {
        synonymous_codons(aa, code).into_iter().max_by_key(|&ii| (self.counts[ii], std::cmp::Reverse(ii)))
    }

    /// The most used codon for an amino acid, as bases.
    pub fn best_codon(&self, aa: u8, code: &GeneticCode) -> Option<[u8; 3]> {
        self.best_codon_index(aa, code).map(codon_bases)
    }

    /// Codons for an amino acid, most used first (ties in NCBI order).
    pub fn ranked_codons(&self, aa: u8, code: &GeneticCode) -> Vec<usize> {
        let mut family = synonymous_codons(aa, code);
        family.sort_by_key(|&ii| (std::cmp::Reverse(self.counts[ii]), ii));
        family
    }

    /// Relative adaptiveness of every codon: its count divided by the largest count in its
    /// family (Sharp and Li 1987). A codon the table never saw is given half a count, the usual
    /// substitution that keeps the geometric mean of [`cai`] finite for sparse tables.
    pub fn cai_weights(&self, code: &GeneticCode) -> [f64; CODON_COUNT] {
        let mut weights = [0.0; CODON_COUNT];
        let count_of = |ii: usize| if self.counts[ii] == 0 { 0.5 } else { self.counts[ii] as f64 };
        for (ii, weight) in weights.iter_mut().enumerate() {
            let family = synonymous_codons(code.amino_acids[ii], code);
            let best = family.iter().map(|&iii| count_of(iii)).fold(0.0f64, f64::max);
            *weight = if best > 0.0 { count_of(ii) / best } else { 0.0 };
        }
        weights
    }

    /// Whether `name` is this table's key, one of its aliases, its organism or its taxonomy id.
    pub fn matches_name(&self, name: &str) -> bool {
        let name = name.trim().to_ascii_lowercase();
        name == self.key || name == self.taxid.to_string() || name == self.organism.to_ascii_lowercase() || self.aliases.contains(&name)
    }

    /// Number of codons in the amino acid's family under `code`.
    fn family_size(&self, index: usize, code: &GeneticCode) -> usize {
        synonymous_codons(code.amino_acids[index], code).len()
    }
}

/// The built-in usage tables, parsed once on first use.
pub fn builtin() -> &'static [CodonTable] {
    static TABLES: OnceLock<Vec<CodonTable>> = OnceLock::new();
    TABLES.get_or_init(|| parse_tsv(BUILTIN_TSV).expect("built-in codon usage table is valid")).as_slice()
}

/// The built-in table with this key, alias, organism name or taxonomy id (case-insensitive).
pub fn table(name: &str) -> Option<&'static CodonTable> {
    builtin().iter().find(|t| t.matches_name(name))
}

/// Parse SlopGene's codon usage TSV: `#` lines are comments, one header row names the columns,
/// and each later row is one organism (see `data/codon_usage.tsv`).
pub fn parse_tsv(text: &str) -> Result<Vec<CodonTable>> {
    let mut tables = Vec::new();
    for (line_no, line) in text.lines().enumerate() {
        // Comments carry the sources; the header row names the 64 codon columns.
        if line.starts_with('#') || line.trim().is_empty() || line.starts_with("key\t") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let where_ = || format!("codon usage table line {}", line_no + 1);
        if fields.len() != 8 + CODON_COUNT {
            return Err(Error::Parse(format!("{}: {} columns, expected {}", where_(), fields.len(), 8 + CODON_COUNT)));
        }
        let number = |text: &str| -> Result<u64> { text.trim().parse::<u64>().map_err(|_| Error::Parse(format!("{}: '{text}' is not a number", where_()))) };
        let mut counts = [0u64; CODON_COUNT];
        for (ii, field) in fields[8..].iter().enumerate() {
            counts[ii] = number(field)?;
        }
        tables.push(CodonTable {
            key: fields[0].trim().to_ascii_lowercase(),
            aliases: fields[1].split(',').map(|a| a.trim().to_ascii_lowercase()).filter(|a| !a.is_empty()).collect(),
            organism: fields[2].trim().to_string(),
            taxid: number(fields[3])? as u32,
            division: fields[4].trim().to_string(),
            cds: number(fields[5])?,
            codons: number(fields[6])?,
            source: fields[7].trim().to_string(),
            counts,
        });
    }
    if tables.is_empty() {
        return Err(Error::Parse("codon usage table holds no organisms".into()));
    }
    Ok(tables)
}

/// How [`reverse_translate`] picks among the codons for one amino acid.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ReverseMode {
    /// Always the host's most used codon. The same protein always gives the same sequence.
    #[default]
    MostFrequent,
    /// Pick at random, weighted by how often the host uses each codon, from this seed. The same
    /// seed always gives the same sequence.
    Weighted { seed: u64 },
}

/// One step of SplitMix64: advance the state and return the next value.
fn next_random(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Write a coding sequence for `protein`, choosing codons as `mode` says. `*` is a stop codon.
///
/// An amino acid the genetic code has no codon for becomes `NNN`, so the result always has three
/// bases per residue and the caller can see what could not be written.
pub fn reverse_translate(protein: &[u8], table: &CodonTable, code: &GeneticCode, mode: ReverseMode) -> Vec<u8> {
    let mut out = Vec::with_capacity(protein.len() * 3);
    let mut state = match mode {
        ReverseMode::Weighted { seed } => seed,
        ReverseMode::MostFrequent => 0,
    };
    // Cache each amino acid's codons, so a long protein does not rescan the code per residue.
    let mut families: HashMap<u8, Vec<usize>> = HashMap::new();
    for &aa in protein {
        let aa = aa.to_ascii_uppercase();
        let family = families.entry(aa).or_insert_with(|| table.ranked_codons(aa, code));
        if family.is_empty() {
            out.extend_from_slice(b"NNN");
            continue;
        }
        let chosen = match mode {
            // `ranked_codons` puts the most used codon first.
            ReverseMode::MostFrequent => family[0],
            ReverseMode::Weighted { .. } => {
                // Walk the family until the weighted draw is used up. A family the table never
                // saw has no weights to draw from, so fall back to its most used codon.
                let total: u64 = family.iter().map(|&ii| table.counts[ii]).sum();
                if total == 0 {
                    family[0]
                } else {
                    let mut draw = next_random(&mut state) % total;
                    let mut pick = family[0];
                    for &ii in family.iter() {
                        if draw < table.counts[ii] {
                            pick = ii;
                            break;
                        }
                        draw -= table.counts[ii];
                    }
                    pick
                }
            }
        };
        out.extend_from_slice(&codon_bases(chosen));
    }
    out
}

/// Write a degenerate coding sequence for `protein`: at every codon position the IUPAC code for
/// every base the amino acid's codons allow there.
///
/// This is what a degenerate oligo needs. Such an oligo must anneal whatever synonymous codon
/// the target happens to use, so the code at each position is the union over the whole family.
/// For amino acids whose codons do not fill a neat box (Leu, Arg, Ser and the stops under the
/// standard code) the union also matches codons of other amino acids: Leu becomes `YTN`, which
/// matches the Phe codons too. An amino acid the code has no codon for becomes `NNN`.
pub fn degenerate_reverse_translate(protein: &[u8], code: &GeneticCode) -> Vec<u8> {
    // IUPAC code for a set of bases, as the bit mask A = 1, C = 2, G = 4, T = 8.
    const CODES: [(u8, u8); 15] = [
        (1, b'A'),
        (2, b'C'),
        (4, b'G'),
        (8, b'T'),
        (1 | 4, b'R'),
        (2 | 8, b'Y'),
        (2 | 4, b'S'),
        (1 | 8, b'W'),
        (4 | 8, b'K'),
        (1 | 2, b'M'),
        (2 | 4 | 8, b'B'),
        (1 | 4 | 8, b'D'),
        (1 | 2 | 8, b'H'),
        (1 | 2 | 4, b'V'),
        (1 | 2 | 4 | 8, b'N'),
    ];
    let mut out = Vec::with_capacity(protein.len() * 3);
    let mut families: HashMap<u8, [u8; 3]> = HashMap::new();
    for &aa in protein {
        let aa = aa.to_ascii_uppercase();
        let codon = *families.entry(aa).or_insert_with(|| {
            let family = synonymous_codons(aa, code);
            if family.is_empty() {
                return *b"NNN";
            }
            // Collect the bases each position takes across the family, then name the union.
            let mut masks = [0u8; 3];
            for &index in &family {
                let bases = codon_bases(index);
                for ii in 0..3 {
                    masks[ii] |= alphabet::base_mask(bases[ii]);
                }
            }
            let mut codon = *b"NNN";
            for ii in 0..3 {
                codon[ii] = CODES.iter().find(|&&(mask, _)| mask == masks[ii]).map(|&(_, letter)| letter).unwrap_or(b'N');
            }
            codon
        });
        out.extend_from_slice(&codon);
    }
    out
}

/// Codon adaptation index of a coding sequence against a host's usage (Sharp and Li 1987).
///
/// The result is the geometric mean of the codons' relative adaptiveness, between 0 and 1, where
/// 1 means every codon is the host's favourite. Stop codons, single-codon families (Met and Trp
/// under the standard code), ambiguous codons and a trailing partial codon are left out. A
/// sequence with nothing left to count gives 0.
pub fn cai(cds: &[u8], table: &CodonTable, code: &GeneticCode) -> f64 {
    let weights = table.cai_weights(code);
    let mut sum_log = 0.0;
    let mut counted = 0usize;
    for codon in cds.as_chunks::<3>().0 {
        let Some(index) = codon_index(codon) else { continue };
        // A stop carries no choice of its own, and a one-codon family always scores 1.
        if code.amino_acids[index] == b'*' || table.family_size(index, code) < 2 {
            continue;
        }
        if weights[index] <= 0.0 {
            continue;
        }
        sum_log += weights[index].ln();
        counted += 1;
    }
    if counted == 0 { 0.0 } else { (sum_log / counted as f64).exp() }
}

/// Fraction of G + C in a sequence (see [`crate::alphabet::gc_fraction`]).
pub fn gc_fraction(seq: &[u8]) -> f64 {
    alphabet::gc_fraction(seq)
}

/// Fraction of G + C at the third base of each complete codon (GC3). A trailing partial codon is
/// left out.
pub fn gc3_fraction(cds: &[u8]) -> f64 {
    let thirds: Vec<u8> = cds.as_chunks::<3>().0.iter().map(|codon| codon[2]).collect();
    alphabet::gc_fraction(&thirds)
}

/// How [`optimize`] chooses a codon before the constraints are repaired.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum CodonChoice {
    /// The host's most used codon for every amino acid.
    #[default]
    MostFrequent,
    /// Codons drawn at random, weighted by the host's usage, so the new sequence has the host's
    /// mix of codons rather than one codon per amino acid. Codons below
    /// [`OptimizeOptions::min_relative_usage`] are left out of the draw.
    Matched,
}

/// What [`optimize`] must respect while it rewrites a coding sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct OptimizeOptions {
    pub choice: CodonChoice,
    /// Seed for [`CodonChoice::Matched`], so a run can be repeated.
    pub seed: u64,
    /// Enzymes whose sites the new sequence must not contain, by name (resolved through
    /// [`EnzymeDb::builtin`]).
    pub avoid_enzymes: Vec<String>,
    /// Longest run of one base allowed. 0 means no limit.
    pub max_run: usize,
    /// Allowed fraction of G + C in every window, as `(low, high)`. `None` means no limit.
    pub gc_range: Option<(f64, f64)>,
    /// Width of that sliding window in bases.
    pub gc_window: usize,
    /// Smallest relative synonymous usage a chosen codon may have, so the host's rare codons are
    /// avoided. 0 allows any codon.
    pub min_relative_usage: f64,
    /// Upper bound on the repair passes, so the search always ends.
    pub max_passes: usize,
}

impl Default for OptimizeOptions {
    fn default() -> Self {
        // Defaults follow what DNA synthesis vendors ask for: runs of more than nine identical
        // bases and GC outside 30-70% make a fragment hard to synthesise and to sequence, and
        // codons used in under a tenth of their family are the ones hosts translate slowly.
        OptimizeOptions {
            choice: CodonChoice::MostFrequent,
            seed: 1,
            avoid_enzymes: Vec::new(),
            max_run: 9,
            gc_range: Some((0.30, 0.70)),
            gc_window: 50,
            min_relative_usage: 0.10,
            max_passes: 20,
        }
    }
}

/// One codon [`optimize`] rewrote.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CodonChange {
    /// 0-based codon number inside the coding sequence.
    pub codon: usize,
    /// The codon's bases inside the coding sequence, 0-based and half-open.
    pub span: Span,
    pub old: [u8; 3],
    pub new: [u8; 3],
    pub amino_acid: u8,
}

/// What [`optimize`] produced.
#[derive(Clone, Debug, PartialEq)]
pub struct Optimized {
    /// The rewritten coding sequence.
    pub seq: Vec<u8>,
    /// The protein, which is the same before and after.
    pub protein: Vec<u8>,
    pub cai_before: f64,
    pub cai_after: f64,
    pub gc_before: f64,
    pub gc_after: f64,
    pub gc3_before: f64,
    pub gc3_after: f64,
    pub changes: Vec<CodonChange>,
    /// Constraints still broken when the search stopped, in words for a person to read, with
    /// 1-based positions.
    pub unmet: Vec<String>,
}

/// Something about a candidate sequence that breaks a constraint.
#[derive(Clone, Debug, PartialEq)]
enum Problem {
    /// A site of an enzyme that should not be there.
    Site { enzyme: String, start: usize, len: usize },
    /// A run of one base longer than allowed.
    Run { start: usize, len: usize, base: u8 },
    /// A window whose GC lies outside the range.
    Gc { start: usize, end: usize, fraction: f64, above: bool },
    /// A codon the host hardly uses.
    Rare { codon: usize, bases: [u8; 3], usage: f64 },
}

impl Problem {
    /// The bases the problem covers, as `start..end` inside the coding sequence.
    fn range(&self) -> (usize, usize) {
        match self {
            Problem::Site { start, len, .. } => (*start, start + len),
            Problem::Run { start, len, .. } => (*start, start + len),
            Problem::Gc { start, end, .. } => (*start, *end),
            Problem::Rare { codon, .. } => (codon * 3, codon * 3 + 3),
        }
    }

    /// The problem in words, with 1-based positions as people read them.
    fn describe(&self) -> String {
        match self {
            Problem::Site { enzyme, start, len } => format!("{enzyme} site at {}..{} could not be removed", start + 1, start + len),
            Problem::Run { start, len, base } => format!("run of {len} {} at {}..{} could not be shortened", *base as char, start + 1, start + len),
            Problem::Gc { start, end, fraction, above } => {
                let side = if *above { "above" } else { "below" };
                format!("GC {:.0}% in {}..{} stays {side} the range", fraction * 100.0, start + 1, end)
            }
            Problem::Rare { codon, bases, usage } => {
                format!("codon {} ({}) stays at {:.0}% of its family", codon + 1, String::from_utf8_lossy(bases), usage * 100.0)
            }
        }
    }
}

/// Everything the problem search needs, gathered once per [`optimize`] call.
struct Constraints<'a> {
    enzymes: Vec<&'a Enzyme>,
    max_run: usize,
    gc_range: Option<(f64, f64)>,
    gc_window: usize,
    min_usage: f64,
    usage: [f64; CODON_COUNT],
}

/// How badly a candidate sequence breaks the constraints, as a cost to drive down.
///
/// The parts are compared in order, so a site or a long run is always worth clearing even when
/// that pushes GC further out. The GC term counts every window on its own, unlike the stretches
/// [`find_problems`] reports: merging two bad stretches into one would otherwise look like an
/// improvement and stop the search early.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Score {
    /// Unwanted enzyme sites plus runs that are too long: what a synthesis vendor refuses.
    hard: usize,
    /// How far the GC windows stray outside the range, added up.
    gc_excess: f64,
    /// Codons the host hardly uses.
    rare: usize,
}

impl Score {
    /// Whether this candidate is better than `other`: fewer hard problems first, then less GC
    /// strayed, then fewer rare codons.
    fn better_than(&self, other: &Score) -> bool {
        if self.hard != other.hard {
            return self.hard < other.hard;
        }
        if (self.gc_excess - other.gc_excess).abs() > 1e-12 {
            return self.gc_excess < other.gc_excess;
        }
        self.rare < other.rare
    }
}

/// Score a candidate sequence against the constraints.
fn cost(seq: &[u8], limits: &Constraints) -> Score {
    let mut hard = 0usize;
    if !limits.enzymes.is_empty() {
        hard += crate::enzymes::find_sites(seq, false, &limits.enzymes).len();
    }
    if limits.max_run > 0 {
        let mut start = 0usize;
        for ii in 1..=seq.len() {
            if ii == seq.len() || !seq[ii].eq_ignore_ascii_case(&seq[start]) {
                hard += usize::from(ii - start > limits.max_run);
                start = ii;
            }
        }
    }
    let mut gc_excess = 0.0;
    if let Some((low, high)) = limits.gc_range {
        let width = limits.gc_window.max(1).min(seq.len().max(1));
        for start in 0..=seq.len().saturating_sub(width) {
            let fraction = alphabet::gc_fraction(&seq[start..start + width]);
            gc_excess += (fraction - high).max(0.0) + (low - fraction).max(0.0);
        }
    }
    let mut rare = 0usize;
    if limits.min_usage > 0.0 {
        for bases in seq.as_chunks::<3>().0 {
            if let Some(index) = codon_index(bases) {
                rare += usize::from(limits.usage[index] < limits.min_usage);
            }
        }
    }
    Score { hard, gc_excess, rare }
}

/// Every constraint `seq` breaks, in position order.
fn find_problems(seq: &[u8], limits: &Constraints) -> Vec<Problem> {
    let mut problems = Vec::new();
    // Unwanted enzyme sites. The coding sequence is handled on its own, so it is linear here.
    if !limits.enzymes.is_empty() {
        for site in crate::enzymes::find_sites(seq, false, &limits.enzymes) {
            problems.push(Problem::Site { enzyme: site.enzyme, start: site.start, len: site.len });
        }
    }
    // Runs of one base: walk the sequence and close a run when the base changes.
    if limits.max_run > 0 {
        let mut start = 0usize;
        for ii in 1..=seq.len() {
            if ii == seq.len() || !seq[ii].eq_ignore_ascii_case(&seq[start]) {
                if ii - start > limits.max_run {
                    problems.push(Problem::Run { start, len: ii - start, base: seq[start].to_ascii_uppercase() });
                }
                start = ii;
            }
        }
    }
    // GC in a sliding window. Neighbouring windows that break the same side of the range are
    // merged, so one GC-rich stretch is one problem rather than dozens.
    if let Some((low, high)) = limits.gc_range {
        let width = limits.gc_window.max(1).min(seq.len().max(1));
        let mut merged: Option<Problem> = None;
        for start in 0..=seq.len().saturating_sub(width) {
            let end = start + width;
            let fraction = alphabet::gc_fraction(&seq[start..end]);
            let above = fraction > high;
            if !above && fraction >= low {
                if let Some(problem) = merged.take() {
                    problems.push(problem);
                }
                continue;
            }
            match &mut merged {
                // Extend the stretch while it breaks the same side of the range.
                Some(Problem::Gc { end: last_end, fraction: worst, above: last_above, .. }) if *last_above == above && start <= *last_end => {
                    *last_end = end;
                    if (above && fraction > *worst) || (!above && fraction < *worst) {
                        *worst = fraction;
                    }
                }
                _ => {
                    if let Some(problem) = merged.take() {
                        problems.push(problem);
                    }
                    merged = Some(Problem::Gc { start, end, fraction, above });
                }
            }
        }
        if let Some(problem) = merged.take() {
            problems.push(problem);
        }
    }
    // Codons the host hardly uses.
    if limits.min_usage > 0.0 {
        for (codon, bases) in seq.as_chunks::<3>().0.iter().enumerate() {
            if let Some(index) = codon_index(bases)
                && limits.usage[index] < limits.min_usage
            {
                problems.push(Problem::Rare { codon, bases: *bases, usage: limits.usage[index] });
            }
        }
    }
    problems.sort_by_key(|problem| problem.range());
    problems
}

/// Rewrite a coding sequence for a host, keeping the protein exactly.
///
/// Codons are chosen as `options.choice` says, then the constraints are repaired codon by codon:
/// each broken constraint is attacked by trying the synonymous codons of the codons it covers,
/// most used first, and a change is kept only when it leaves fewer constraints broken than
/// before. The repair stops when nothing is broken, when no single codon change improves matters
/// or after [`OptimizeOptions::max_passes`] passes, so it always ends. Whatever is still broken
/// is listed in [`Optimized::unmet`].
///
/// The sequence must be a whole number of codons. Codons with ambiguous bases are kept as they
/// are and reported in `unmet`.
pub fn optimize(cds: &[u8], table: &CodonTable, code: &GeneticCode, options: &OptimizeOptions) -> Result<Optimized> {
    if cds.is_empty() {
        return Err(invalid("cannot optimise an empty coding sequence"));
    }
    if !cds.len().is_multiple_of(3) {
        return Err(invalid(format!("coding sequence of {} bases is not a whole number of codons", cds.len())));
    }
    let original = alphabet::upper(cds);
    let codon_count = original.len() / 3;
    let avoid = EnzymeDb::builtin().resolve(&options.avoid_enzymes)?;

    // Relative usage of every codon, and the candidates for every codon of the sequence.
    let mut usage = [0.0f64; CODON_COUNT];
    for (ii, share) in usage.iter_mut().enumerate() {
        *share = table.relative_usage(ii, code);
    }
    let limits = Constraints {
        enzymes: avoid,
        max_run: options.max_run,
        gc_range: options.gc_range,
        gc_window: options.gc_window,
        min_usage: options.min_relative_usage,
        usage,
    };

    // Protein, and the candidate codons of each position: synonymous, most used first, and above
    // the rare-codon threshold unless nothing is.
    let mut protein = Vec::with_capacity(codon_count);
    let mut candidates: Vec<Vec<usize>> = Vec::with_capacity(codon_count);
    let mut unmet: Vec<String> = Vec::new();
    for (codon, bases) in original.as_chunks::<3>().0.iter().enumerate() {
        match codon_index(bases) {
            Some(index) => {
                let aa = code.amino_acids[index];
                protein.push(aa);
                let ranked = table.ranked_codons(aa, code);
                let allowed: Vec<usize> = ranked.iter().copied().filter(|&ii| usage[ii] >= options.min_relative_usage).collect();
                candidates.push(if allowed.is_empty() { ranked } else { allowed });
            }
            None => {
                // An ambiguous codon has no known amino acid, so it cannot be rewritten safely.
                protein.push(b'X');
                candidates.push(Vec::new());
                unmet.push(format!("codon {} ({}) has ambiguous bases and was left alone", codon + 1, String::from_utf8_lossy(bases)));
            }
        }
    }

    // First pass: the preferred codon for every amino acid.
    let mut seq = original.clone();
    let mut state = options.seed;
    for codon in 0..codon_count {
        let choices = &candidates[codon];
        if choices.is_empty() {
            continue;
        }
        let chosen = match options.choice {
            CodonChoice::MostFrequent => choices[0],
            CodonChoice::Matched => {
                // Weighted draw over the allowed codons, as in `reverse_translate`.
                let total: u64 = choices.iter().map(|&ii| table.counts[ii]).sum();
                if total == 0 {
                    choices[0]
                } else {
                    let mut draw = next_random(&mut state) % total;
                    let mut pick = choices[0];
                    for &ii in choices {
                        if draw < table.counts[ii] {
                            pick = ii;
                            break;
                        }
                        draw -= table.counts[ii];
                    }
                    pick
                }
            }
        };
        seq[codon * 3..codon * 3 + 3].copy_from_slice(&codon_bases(chosen));
    }

    // Repair: each pass walks the broken constraints and tries the synonymous codons of the
    // codons each one covers, keeping a substitution only when the whole sequence scores better.
    // A pass that improves nothing ends the search, so it always ends.
    let mut problems = find_problems(&seq, &limits);
    let mut score = cost(&seq, &limits);
    for _ in 0..options.max_passes {
        if problems.is_empty() {
            break;
        }
        let mut improved = false;
        for problem in &problems {
            let (lo, hi) = problem.range();
            if hi == 0 {
                continue;
            }
            // The codons the problem covers, from its start to its last base.
            let first_codon = lo / 3;
            let last_codon = ((hi - 1) / 3).min(codon_count.saturating_sub(1));
            let mut fixed = false;
            for (codon, choices) in candidates.iter().enumerate().take(last_codon + 1).skip(first_codon) {
                let at = codon * 3;
                let current: [u8; 3] = [seq[at], seq[at + 1], seq[at + 2]];
                for &candidate in choices {
                    let bases = codon_bases(candidate);
                    if bases == current {
                        continue;
                    }
                    seq[at..at + 3].copy_from_slice(&bases);
                    let next = cost(&seq, &limits);
                    if next.better_than(&score) {
                        score = next;
                        improved = true;
                        fixed = true;
                        break;
                    }
                    seq[at..at + 3].copy_from_slice(&current);
                }
                if fixed {
                    break;
                }
            }
        }
        if !improved {
            break;
        }
        problems = find_problems(&seq, &limits);
    }
    unmet.extend(problems.iter().map(Problem::describe));

    // What changed, codon by codon, against the sequence that came in.
    let mut changes = Vec::new();
    for (codon, &amino_acid) in protein.iter().enumerate() {
        let at = codon * 3;
        let old: [u8; 3] = [original[at], original[at + 1], original[at + 2]];
        let new: [u8; 3] = [seq[at], seq[at + 1], seq[at + 2]];
        if old != new {
            changes.push(CodonChange { codon, span: Span::new(at, at + 3), old, new, amino_acid });
        }
    }

    Ok(Optimized {
        cai_before: cai(&original, table, code),
        cai_after: cai(&seq, table, code),
        gc_before: gc_fraction(&original),
        gc_after: gc_fraction(&seq),
        gc3_before: gc3_fraction(&original),
        gc3_after: gc3_fraction(&seq),
        seq,
        protein,
        changes,
        unmet,
    })
}

/// Whether a silent change adds a restriction site or takes one away.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SilentKind {
    Add,
    Remove,
}

impl SilentKind {
    pub fn name(self) -> &'static str {
        match self {
            SilentKind::Add => "add",
            SilentKind::Remove => "remove",
        }
    }
}

/// One base swapped for another, in top-strand coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Substitution {
    /// 0-based position on the top strand, whichever strand the coding sequence reads from.
    pub pos: usize,
    pub from: u8,
    pub to: u8,
}

/// One codon a silent change rewrites, read 5'->3' along the coding sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CodonEdit {
    /// 0-based codon number inside the coding sequence.
    pub codon: usize,
    /// Top-strand positions of the codon's three bases, in reading order. For a reverse-strand
    /// coding sequence these count down.
    pub positions: [usize; 3],
    pub old: [u8; 3],
    pub new: [u8; 3],
    pub amino_acid: u8,
}

/// A restriction site, named by where it is rather than by what cuts there.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiteRef {
    pub enzyme: String,
    /// First base of the recognition site on the top strand (0-based).
    pub start: usize,
    pub len: usize,
    /// `Forward` when the site reads 5'->3' on the top strand.
    pub strand: Strand,
}

/// One way to add or remove a restriction site without changing the protein.
#[derive(Clone, Debug, PartialEq)]
pub struct SilentChange {
    pub enzyme: String,
    pub kind: SilentKind,
    /// The bases to swap, in position order.
    pub subs: Vec<Substitution>,
    /// The codons those bases sit in, with the amino acid they keep.
    pub codons: Vec<CodonEdit>,
    /// The site that appears or disappears, in top-strand coordinates. The span runs through the
    /// origin when the site does.
    pub site: Span,
    /// Orientation of that site.
    pub strand: Strand,
    /// Where the enzyme cuts the top strand of that site, as a boundary (the strand is cut
    /// between bases `cut - 1` and `cut`). `None` when the cut falls off a linear molecule.
    pub cut: Option<usize>,
    /// Other sites of the given enzymes the same change creates.
    pub also_created: Vec<SiteRef>,
    /// Other sites of the given enzymes the same change destroys.
    pub also_destroyed: Vec<SiteRef>,
}

impl SilentChange {
    /// How many other sites the change disturbs.
    pub fn side_effects(&self) -> usize {
        self.also_created.len() + self.also_destroyed.len()
    }
}

/// Length of the site prefix used to index enzymes, as in [`crate::enzymes::find_sites`].
const PREFIX: usize = 4;
/// Prefixes that stand for more concrete 4-mers than this are checked at every position.
const MAX_PREFIX_FORMS: usize = 16;

/// One orientation of one enzyme's recognition site.
struct Probe<'a> {
    enzyme: &'a Enzyme,
    pattern: Vec<u8>,
    strand: Strand,
}

/// One match of a probe.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Found {
    probe: u32,
    start: usize,
}

/// Enzyme sites indexed by the first four bases of each pattern, so a short window can be
/// rescanned cheaply.
///
/// [`crate::enzymes::find_sites`] indexes the same way but rebuilds its index on every call,
/// which is far too slow here: a silent-change search rescans a window for each of thousands of
/// candidate edits, so the index is built once and kept.
struct SiteIndex<'a> {
    probes: Vec<Probe<'a>>,
    by_prefix: Vec<Vec<u32>>,
    everywhere: Vec<u32>,
    /// Longest pattern, which is how far a substitution can reach into a site.
    max_site: usize,
}

/// Rank of a concrete base for the prefix index (A = 0, C = 1, G = 2, T = 3).
fn base_rank(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// Every concrete 4-mer code an IUPAC prefix stands for, or `None` when there are too many.
fn prefix_codes(prefix: &[u8]) -> Option<Vec<usize>> {
    let mut codes = vec![0usize];
    for &base in prefix {
        let mask = alphabet::base_mask(base);
        let mut next = Vec::with_capacity(codes.len() * 4);
        for code in &codes {
            for (bit, rank) in [(1u8, 0usize), (2, 1), (4, 2), (8, 3)] {
                if mask & bit != 0 {
                    next.push(code * 4 + rank);
                }
            }
        }
        if next.len() > MAX_PREFIX_FORMS {
            return None;
        }
        codes = next;
    }
    Some(codes)
}

impl<'a> SiteIndex<'a> {
    /// Index both orientations of every enzyme that can fit in the sequence.
    fn new(enzymes: &[&'a Enzyme], seq_len: usize, circular: bool) -> SiteIndex<'a> {
        let mut probes: Vec<Probe<'a>> = Vec::new();
        for &enzyme in enzymes {
            let n = enzyme.site.len();
            if n == 0 || (n > seq_len && !circular) {
                continue;
            }
            let reverse = alphabet::revcomp(&enzyme.site);
            let palindromic = reverse == enzyme.site;
            probes.push(Probe { enzyme, pattern: enzyme.site.clone(), strand: Strand::Forward });
            if !palindromic {
                probes.push(Probe { enzyme, pattern: reverse, strand: Strand::Reverse });
            }
        }
        let mut by_prefix: Vec<Vec<u32>> = vec![Vec::new(); 256];
        let mut everywhere: Vec<u32> = Vec::new();
        for (ii, probe) in probes.iter().enumerate() {
            match (probe.pattern.len() >= PREFIX).then(|| prefix_codes(&probe.pattern[..PREFIX])).flatten() {
                Some(codes) => codes.into_iter().for_each(|code| by_prefix[code].push(ii as u32)),
                None => everywhere.push(ii as u32),
            }
        }
        let max_site = probes.iter().map(|probe| probe.pattern.len()).max().unwrap_or(0);
        SiteIndex { probes, by_prefix, everywhere, max_site }
    }

    /// Match one probe at one position.
    fn try_probe(&self, probe_index: u32, start: usize, seq: &[u8], circular: bool, out: &mut Vec<Found>) {
        let probe = &self.probes[probe_index as usize];
        let n = probe.pattern.len();
        let seq_len = seq.len();
        if !circular && start + n > seq_len {
            return;
        }
        if !(0..n).all(|ii| iupac_matches(probe.pattern[ii], seq[(start + ii) % seq_len])) {
            return;
        }
        out.push(Found { probe: probe_index, start });
    }

    /// Every site that starts at one of `count` positions from `from` (positions wrap on a
    /// circular sequence). The answer is what [`crate::enzymes::find_sites`] reports for the same
    /// positions, including its fallback to every probe where the sequence is ambiguous.
    fn scan(&self, seq: &[u8], circular: bool, from: usize, count: usize, out: &mut Vec<Found>) {
        let seq_len = seq.len();
        if seq_len == 0 {
            return;
        }
        for ii in 0..count.min(seq_len) {
            let start = (from + ii) % seq_len;
            let window_fits = circular || start + PREFIX <= seq_len;
            let code =
                if window_fits { (0..PREFIX).try_fold(0usize, |acc, iii| base_rank(seq[(start + iii) % seq_len]).map(|rank| acc * 4 + rank)) } else { None };
            match code {
                Some(code) => {
                    for &probe in &self.by_prefix[code] {
                        self.try_probe(probe, start, seq, circular, out);
                    }
                    for &probe in &self.everywhere {
                        self.try_probe(probe, start, seq, circular, out);
                    }
                }
                // An ambiguous base in the window: compare every probe the slow way.
                None if window_fits => {
                    for probe in 0..self.probes.len() as u32 {
                        self.try_probe(probe, start, seq, circular, out);
                    }
                }
                // Near the end of a linear sequence only probes shorter than the prefix can fit.
                None => {
                    for &probe in &self.everywhere {
                        self.try_probe(probe, start, seq, circular, out);
                    }
                }
            }
        }
    }

    /// Describe a match.
    fn site_ref(&self, found: Found) -> SiteRef {
        let probe = &self.probes[found.probe as usize];
        SiteRef { enzyme: probe.enzyme.name.clone(), start: found.start, len: probe.pattern.len(), strand: probe.strand }
    }

    /// Where the enzyme cuts the top strand for this match, as a boundary.
    ///
    /// The reasoning is that of `place_cut` in [`crate::enzymes`]: a site on the bottom strand is
    /// the enzyme seen from the other side, so its boundaries mirror around the site. An enzyme
    /// with two cuts gives the first one that lands inside a linear molecule, which is what
    /// [`crate::enzymes::find_sites`] reports first.
    fn cut_position(&self, found: Found, seq_len: usize, circular: bool) -> Option<usize> {
        let probe = &self.probes[found.probe as usize];
        let n = probe.pattern.len() as i64;
        let start = found.start as i64;
        let len = seq_len as i64;
        probe.enzyme.cuts.iter().find_map(|&(top, bottom)| {
            let (top_abs, bottom_abs) = match probe.strand {
                Strand::Reverse => (start + n - bottom as i64, start + n - top as i64),
                _ => (start + top as i64, start + bottom as i64),
            };
            if circular {
                Some(top_abs.rem_euclid(len) as usize)
            } else {
                // Both strands must be cut inside the molecule, and not at its very ends.
                (top_abs.min(bottom_abs) > 0 && top_abs.max(bottom_abs) < len).then_some(top_abs as usize)
            }
        })
    }
}

/// One codon of the coding sequence, laid out on the top strand.
struct CdsCodon {
    /// Top-strand positions of the codon's bases, in reading order.
    positions: [usize; 3],
    /// The codon as read along the coding sequence.
    bases: [u8; 3],
    amino_acid: u8,
}

/// One synonymous rewrite of one codon.
struct Alternative {
    codon: usize,
    /// Bases to swap on the top strand, in position order.
    subs: Vec<Substitution>,
    /// The new codon, as read along the coding sequence.
    bases: [u8; 3],
}

/// Silent substitutions that add or remove a restriction site inside a coding sequence.
///
/// `cds` is the coding sequence's extent in top-strand coordinates, which may run through the
/// origin of a circular molecule, and `reverse` says it is read on the bottom strand. Only bases
/// inside the coding sequence are changed, and every change keeps the protein exactly. For each
/// of `enzymes` the result holds the changes of at most `max_changes` bases (1 or 2) that create
/// a site of it, and the changes that destroy an existing site of it that overlaps the coding
/// sequence. Each change also lists the other sites of the same enzymes it creates or destroys,
/// so a side effect is visible before the edit is made.
///
/// Changes are sorted by the number of substitutions, then by the number of side effects, then by
/// position. Where one site can be created or destroyed both by one substitution and by two, only
/// the one-substitution ways are kept. A trailing partial codon, a codon with ambiguous bases and
/// a coding sequence that runs off the end of a linear molecule yield no changes.
pub fn silent_site_changes(
    seq: &[u8],
    circular: bool,
    cds: Span,
    reverse: bool,
    enzymes: &[&Enzyme],
    code: &GeneticCode,
    max_changes: usize,
) -> Vec<SilentChange> {
    let seq_len = seq.len();
    let max_changes = max_changes.clamp(1, 2);
    if seq_len == 0 || cds.len() < 3 || enzymes.is_empty() || (!circular && cds.end > seq_len) {
        return Vec::new();
    }
    let work = alphabet::upper(seq);

    // Lay the codons out on the top strand. A forward coding sequence reads along the span; a
    // reverse one reads backwards from its last base, with every base complemented.
    let codon_count = cds.len() / 3;
    let mut codons: Vec<CdsCodon> = Vec::with_capacity(codon_count);
    let mut position_of: HashMap<usize, (usize, usize)> = HashMap::with_capacity(codon_count * 3);
    for codon in 0..codon_count {
        let mut positions = [0usize; 3];
        let mut bases = [0u8; 3];
        for offset in 0..3 {
            let along = codon * 3 + offset;
            let pos = if reverse { (cds.end - 1 - along) % seq_len } else { (cds.start + along) % seq_len };
            positions[offset] = pos;
            bases[offset] = if reverse { alphabet::complement_base(work[pos]) } else { work[pos] };
            position_of.insert(pos, (codon, offset));
        }
        let amino_acid = code.translate_codon(&bases);
        codons.push(CdsCodon { positions, bases, amino_acid });
    }

    // Every synonymous rewrite of every codon, within the substitution budget. A codon with
    // ambiguous bases has no known amino acid, so it is left alone.
    let mut alternatives: Vec<Vec<Alternative>> = Vec::with_capacity(codon_count);
    for (codon, layout) in codons.iter().enumerate() {
        let mut alts = Vec::new();
        if let Some(index) = codon_index(&layout.bases) {
            for other in synonymous_codons(code.amino_acids[index], code) {
                if other == index {
                    continue;
                }
                let bases = codon_bases(other);
                // Each differing base becomes one top-strand substitution.
                let mut subs: Vec<Substitution> = (0..3)
                    .filter(|&offset| bases[offset] != layout.bases[offset])
                    .map(|offset| {
                        let to = if reverse { alphabet::complement_base(bases[offset]) } else { bases[offset] };
                        Substitution { pos: layout.positions[offset], from: work[layout.positions[offset]], to }
                    })
                    .collect();
                if subs.is_empty() || subs.len() > max_changes {
                    continue;
                }
                subs.sort_by_key(|sub| sub.pos);
                alts.push(Alternative { codon, subs, bases });
            }
        }
        alternatives.push(alts);
    }

    let index = SiteIndex::new(enzymes, seq_len, circular);
    if index.max_site == 0 {
        return Vec::new();
    }

    // Candidate edits: one rewritten codon, or, when two substitutions are allowed, two codons
    // one substitution each. Two codons only matter when one site can span both substitutions,
    // so distant pairs are skipped: a site that needs only one of them is found on its own.
    let mut candidates: Vec<Vec<&Alternative>> = Vec::new();
    for alts in &alternatives {
        for alt in alts {
            candidates.push(vec![alt]);
        }
    }
    if max_changes == 2 {
        let reach = index.max_site - 1;
        // Codons are 3 bases, so a site of `max_site` bases spans at most this many codons.
        let span_codons = reach / 3 + 2;
        for first in 0..codon_count {
            for second in first + 1..(first + span_codons + 1).min(codon_count) {
                for left in alternatives[first].iter().filter(|alt| alt.subs.len() == 1) {
                    for right in alternatives[second].iter().filter(|alt| alt.subs.len() == 1) {
                        if circular_gap(left.subs[0].pos, right.subs[0].pos, seq_len) <= reach {
                            candidates.push(vec![left, right]);
                        }
                    }
                }
            }
        }
    }

    // Try every candidate: scan the window around it, apply it, scan again, and compare.
    let mut found = Vec::new();
    let mut before = Vec::new();
    let mut after = Vec::new();
    let mut work = work;
    for candidate in &candidates {
        let subs: Vec<Substitution> = {
            let mut subs: Vec<Substitution> = candidate.iter().flat_map(|alt| alt.subs.iter().copied()).collect();
            subs.sort_by_key(|sub| sub.pos);
            subs
        };
        // A site can only change when it covers a substituted base, so a window reaching one
        // site length either side of the substitutions holds every site that can change.
        let reach = index.max_site - 1;
        let (window_start, window_count) = window_around(&subs, reach, seq_len, circular);
        before.clear();
        after.clear();
        index.scan(&work, circular, window_start, window_count, &mut before);
        for sub in &subs {
            work[sub.pos] = sub.to;
        }
        index.scan(&work, circular, window_start, window_count, &mut after);
        for sub in &subs {
            work[sub.pos] = sub.from;
        }
        if before == after {
            continue;
        }
        let created: Vec<Found> = after.iter().copied().filter(|site| !before.contains(site)).collect();
        let destroyed: Vec<Found> = before.iter().copied().filter(|site| !after.contains(site)).collect();
        if created.is_empty() && destroyed.is_empty() {
            continue;
        }
        // The codons this candidate rewrites, with the amino acid each keeps.
        let edits: Vec<CodonEdit> = candidate
            .iter()
            .map(|alt| {
                let layout = &codons[alt.codon];
                CodonEdit { codon: alt.codon, positions: layout.positions, old: layout.bases, new: alt.bases, amino_acid: layout.amino_acid }
            })
            .collect();
        // One entry per site that appears, and one per site that goes away.
        for (kind, sites) in [(SilentKind::Add, &created), (SilentKind::Remove, &destroyed)] {
            for &site in sites {
                let described = index.site_ref(site);
                let others = |list: &[Found]| -> Vec<SiteRef> { list.iter().filter(|&&other| other != site).map(|&other| index.site_ref(other)).collect() };
                found.push(SilentChange {
                    enzyme: described.enzyme.clone(),
                    kind,
                    subs: subs.clone(),
                    codons: edits.clone(),
                    site: Span::new(described.start, described.start + described.len),
                    strand: described.strand,
                    cut: index.cut_position(site, seq_len, circular),
                    also_created: if kind == SilentKind::Add { others(&created) } else { created.iter().map(|&other| index.site_ref(other)).collect() },
                    also_destroyed: if kind == SilentKind::Remove {
                        others(&destroyed)
                    } else {
                        destroyed.iter().map(|&other| index.site_ref(other)).collect()
                    },
                });
            }
        }
    }
    // Keep, for each site, only the ways that need the fewest substitutions, and drop duplicates.
    let mut fewest: HashMap<(String, SilentKind, usize, u8), usize> = HashMap::new();
    for change in &found {
        let key = (change.enzyme.clone(), change.kind, change.site.start, change.strand as u8);
        let entry = fewest.entry(key).or_insert(usize::MAX);
        *entry = (*entry).min(change.subs.len());
    }
    found.retain(|change| fewest[&(change.enzyme.clone(), change.kind, change.site.start, change.strand as u8)] == change.subs.len());
    found.sort_by(|a, b| {
        a.subs
            .len()
            .cmp(&b.subs.len())
            .then_with(|| a.side_effects().cmp(&b.side_effects()))
            .then_with(|| a.subs[0].pos.cmp(&b.subs[0].pos))
            .then_with(|| a.site.start.cmp(&b.site.start))
            .then_with(|| a.enzyme.cmp(&b.enzyme))
            .then_with(|| a.kind.cmp(&b.kind))
            .then_with(|| a.subs.last().map(|s| s.pos).cmp(&b.subs.last().map(|s| s.pos)))
    });
    found.dedup_by(|a, b| a.enzyme == b.enzyme && a.kind == b.kind && a.site == b.site && a.strand == b.strand && a.subs == b.subs);
    found
}

/// Smallest number of bases between two positions, either way round a circular sequence.
fn circular_gap(a: usize, b: usize, seq_len: usize) -> usize {
    let straight = a.abs_diff(b);
    straight.min(seq_len - straight)
}

/// The positions to rescan around a set of substitutions: `reach` bases either side, so every
/// site that can cover a substituted base starts inside the window.
fn window_around(subs: &[Substitution], reach: usize, seq_len: usize, circular: bool) -> (usize, usize) {
    if circular {
        // On a circle the substitutions may straddle the origin, so measure from the first one
        // and stretch the window to cover the gap to the last.
        let first = subs.iter().map(|sub| sub.pos).min().unwrap_or(0);
        let last = subs.iter().map(|sub| sub.pos).max().unwrap_or(0);
        let straight = last - first;
        let (anchor, spread) = if straight * 2 <= seq_len { (first, straight) } else { (last, seq_len - straight) };
        let start = (anchor + seq_len - reach) % seq_len;
        ((start), (spread + 2 * reach + 1).min(seq_len))
    } else {
        let first = subs.iter().map(|sub| sub.pos).min().unwrap_or(0);
        let last = subs.iter().map(|sub| sub.pos).max().unwrap_or(0);
        let start = first.saturating_sub(reach);
        let end = (last + reach + 1).min(seq_len);
        (start, end - start)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzymes::find_sites;
    use crate::translate::translate;

    fn standard() -> &'static GeneticCode {
        GeneticCode::standard()
    }

    #[test]
    fn every_table_is_complete_and_adds_up() {
        let tables = builtin();
        assert_eq!(tables.len(), 13, "{} tables", tables.len());
        for table in tables {
            // The file stores 64 counts, and they add up to the codon total the source states.
            assert_eq!(table.counts.len(), CODON_COUNT);
            assert_eq!(table.total(), table.codons, "{} totals differ", table.key);
            assert!(table.cds > 0 && table.taxid > 0 && table.source.starts_with("https://www.kazusa.or.jp/codon/"), "{}", table.key);
            // Every amino acid's synonymous fractions add up to 1, and the per-thousand
            // frequencies add up to 1000.
            let mut seen = std::collections::HashSet::new();
            for ii in 0..CODON_COUNT {
                let aa = standard().amino_acids[ii];
                if !seen.insert(aa) {
                    continue;
                }
                let sum: f64 = synonymous_codons(aa, standard()).into_iter().map(|iii| table.relative_usage(iii, standard())).sum();
                assert!((sum - 1.0).abs() < 1e-12, "{} {}: fractions add up to {sum}", table.key, aa as char);
            }
            let per_thousand: f64 = (0..CODON_COUNT).map(|ii| table.per_thousand(ii)).sum();
            assert!((per_thousand - 1000.0).abs() < 1e-9, "{}: {per_thousand} per thousand", table.key);
        }
    }

    #[test]
    fn tables_match_the_source_pages() {
        // Values copied from the Kazusa pages named in each table's source column, read
        // 2026-10-23. Each line is [codon] [fraction] [per thousand] ([count]) as the page prints
        // it, so both the stored count and the derived frequency are checked.
        let checks: [(&str, u64, u64, &str, f64, f64, u64); 6] = [
            // https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=316407 (E. coli W3110)
            // "CUG L 0.50 53.1 ( 72898)", 4332 CDS's (1372057 codons).
            ("ecoli", 4332, 1372057, "CTG", 0.50, 53.1, 72898),
            // ?species=9606 (Homo sapiens) "CUG L 0.40 39.6 (1611801)", 93487 CDS's.
            ("human", 93487, 40662582, "CTG", 0.40, 39.6, 1611801),
            // ?species=9606 "AGA R 0.21 12.2 (494682)".
            ("human", 93487, 40662582, "AGA", 0.21, 12.2, 494682),
            // ?species=4932 (S. cerevisiae) "AGA R 0.48 21.3 (139081)", 14411 CDS's.
            ("yeast", 14411, 6534504, "AGA", 0.48, 21.3, 139081),
            // ?species=4922 (Pichia pastoris) "AAA K 0.47 29.9 ( 2433)", 137 CDS's.
            ("pichia", 137, 81301, "AAA", 0.47, 29.9, 2433),
            // ?species=7108 (Spodoptera frugiperda) "UUU F 0.27 10.1 ( 392)", 142 CDS's.
            ("insect", 142, 38777, "TTT", 0.27, 10.1, 392),
        ];
        for (key, cds, codons, codon, fraction, per_thousand, count) in checks {
            let table = table(key).unwrap_or_else(|| panic!("no table {key}"));
            assert_eq!(table.cds, cds, "{key} CDS count");
            assert_eq!(table.codons, codons, "{key} codon count");
            let index = codon_index(codon.as_bytes()).unwrap();
            assert_eq!(table.counts[index], count, "{key} {codon} count");
            assert!((table.per_thousand(index) - per_thousand).abs() <= 0.05, "{key} {codon}: {} per thousand", table.per_thousand(index));
            assert!((table.relative_usage(index, standard()) - fraction).abs() <= 0.005, "{key} {codon}: fraction {}", table.relative_usage(index, standard()));
        }
        // The page for taxid 9606 states 93487 CDS's, which is what the row says it read.
        assert_eq!(table("human").unwrap().taxid, 9606);
        // E. coli uses the K-12 substrain W3110, because Kazusa's taxid 83333 table has 14 CDSs.
        assert_eq!(table("ecoli").unwrap().taxid, 316407);
        assert!(table("ecoli").unwrap().organism.contains("W3110"));
    }

    #[test]
    fn names_aliases_and_taxids_find_tables() {
        for (name, key) in [
            ("human", "human"),
            ("Homo sapiens", "human"),
            ("9606", "human"),
            ("E. coli", "ecoli"),
            ("ecoli", "ecoli"),
            ("CHO", "cho"),
            ("chinese hamster", "cho"),
            ("sf9", "insect"),
            ("komagataella", "pichia"),
            ("S. cerevisiae", "yeast"),
            ("fruit fly", "fly"),
            ("c. elegans", "worm"),
            ("danio rerio", "zebrafish"),
            ("b. subtilis", "bsubtilis"),
            ("thale cress", "arabidopsis"),
            ("n. tabacum", "tobacco"),
            ("mus musculus", "mouse"),
        ] {
            assert_eq!(table(name).map(|t| t.key.as_str()), Some(key), "looking up {name}");
        }
        assert!(table("Tyrannosaurus rex").is_none());
    }

    #[test]
    fn best_codons_and_usage_follow_the_counts() {
        let ecoli = table("ecoli").unwrap();
        // CTG is E. coli's leucine codon and GAA its glutamate codon.
        assert_eq!(ecoli.best_codon(b'L', standard()), Some(*b"CTG"));
        assert_eq!(ecoli.best_codon(b'E', standard()), Some(*b"GAA"));
        // Met and Trp have one codon each, so their usage and RSCU are both 1.
        assert_eq!(ecoli.relative_usage(codon_index(b"ATG").unwrap(), standard()), 1.0);
        assert_eq!(ecoli.rscu(codon_index(b"TGG").unwrap(), standard()), 1.0);
        // An RSCU above 1 means a codon is used more than its share of an even family.
        assert!(ecoli.rscu(codon_index(b"CTG").unwrap(), standard()) > 1.0);
        assert!(ecoli.rscu(codon_index(b"AGA").unwrap(), standard()) < 1.0);
        // Yeast prefers AGA for arginine where E. coli prefers CGC.
        assert_eq!(table("yeast").unwrap().best_codon(b'R', standard()), Some(*b"AGA"));
        assert_eq!(ecoli.best_codon(b'R', standard()), Some(*b"CGC"));
        // The CAI weight of the favourite codon is 1 by definition.
        let weights = ecoli.cai_weights(standard());
        assert_eq!(weights[codon_index(b"CTG").unwrap()], 1.0);
        assert!(weights[codon_index(b"CTA").unwrap()] < 0.2);
    }

    #[test]
    fn reverse_translation_gives_back_the_protein() {
        let protein = b"MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGQIEWFKK*";
        for key in ["ecoli", "human", "yeast", "pichia", "insect"] {
            let table = table(key).unwrap();
            // Several genetic codes: the standard one, the bacterial one and two with reassigned
            // stops, which change which codons each amino acid has.
            for code_id in [1u32, 11, 2, 6] {
                let code = GeneticCode::get(code_id).unwrap();
                let dna = reverse_translate(protein, table, code, ReverseMode::MostFrequent);
                assert_eq!(dna.len(), protein.len() * 3);
                assert_eq!(translate(&dna, code, false), protein.to_vec(), "{key} code {code_id} most frequent");
                for seed in [1u64, 7, 12345] {
                    let dna = reverse_translate(protein, table, code, ReverseMode::Weighted { seed });
                    assert_eq!(translate(&dna, code, false), protein.to_vec(), "{key} code {code_id} seed {seed}");
                }
            }
        }
        // The same seed always gives the same sequence, and different seeds usually differ.
        let table = table("human").unwrap();
        let once = reverse_translate(protein, table, standard(), ReverseMode::Weighted { seed: 99 });
        let twice = reverse_translate(protein, table, standard(), ReverseMode::Weighted { seed: 99 });
        assert_eq!(once, twice);
        assert_ne!(once, reverse_translate(protein, table, standard(), ReverseMode::Weighted { seed: 100 }));
        // Most-frequent mode uses the host's favourite codon everywhere.
        let dna = reverse_translate(b"LLL", table, standard(), ReverseMode::MostFrequent);
        assert_eq!(dna, b"CTGCTGCTG".to_vec());
        // An amino acid the code cannot write becomes NNN.
        assert_eq!(reverse_translate(b"X", table, standard(), ReverseMode::MostFrequent), b"NNN".to_vec());
    }

    #[test]
    fn degenerate_codons_match_every_codon_of_their_amino_acid() {
        for code_id in [1u32, 2, 11] {
            let code = GeneticCode::get(code_id).unwrap();
            let mut seen = std::collections::HashSet::new();
            for ii in 0..CODON_COUNT {
                let aa = code.amino_acids[ii];
                if !seen.insert(aa) {
                    continue;
                }
                let degenerate = degenerate_reverse_translate(&[aa], code);
                assert_eq!(degenerate.len(), 3);
                // Every codon of this amino acid must match the degenerate codon position by
                // position, which is what a degenerate oligo has to anneal to.
                for index in synonymous_codons(aa, code) {
                    let codon = codon_bases(index);
                    for offset in 0..3 {
                        assert!(
                            iupac_matches(degenerate[offset], codon[offset]),
                            "code {code_id} {}: {} does not match {} at {offset}",
                            aa as char,
                            String::from_utf8_lossy(&degenerate),
                            String::from_utf8_lossy(&codon)
                        );
                    }
                }
            }
        }
        // Spot checks under the standard code: the clean families collapse to one code each.
        assert_eq!(degenerate_reverse_translate(b"AGKFLSR*", standard()), b"GCNGGNAARTTYYTNWSNMGNTRR".to_vec());
        assert_eq!(degenerate_reverse_translate(b"MW", standard()), b"ATGTGG".to_vec());
        assert_eq!(degenerate_reverse_translate(b"X", standard()), b"NNN".to_vec());
    }

    #[test]
    fn cai_rates_the_best_and_worst_codons() {
        let ecoli = table("ecoli").unwrap();
        let protein = b"MKVLATTGSLNDIEWQRPHFYACG";
        // A sequence of the host's favourite codons scores 1; the rarest codons score far lower.
        let best = reverse_translate(protein, ecoli, standard(), ReverseMode::MostFrequent);
        assert!((cai(&best, ecoli, standard()) - 1.0).abs() < 1e-12, "{}", cai(&best, ecoli, standard()));
        let mut worst = Vec::new();
        for &aa in protein {
            let rarest = *ecoli.ranked_codons(aa, standard()).last().unwrap();
            worst.extend_from_slice(&codon_bases(rarest));
        }
        let low = cai(&worst, ecoli, standard());
        assert!(low < 0.4, "worst-codon CAI {low}");
        // A host that likes different codons rates the same sequence differently.
        assert!(cai(&best, table("yeast").unwrap(), standard()) < 1.0);
        // GC and GC3 come out of the sequence itself.
        assert!((gc_fraction(b"GGCCAATT") - 0.5).abs() < 1e-12);
        assert!((gc3_fraction(b"AAGAAGAAT") - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn optimisation_keeps_the_protein_and_raises_cai() {
        let ecoli = table("ecoli").unwrap();
        let protein = b"MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQ*";
        // Start from the host's rarest codons, then plant an EcoRI site that a silent change can
        // clear: GAA TTC codes for Glu-Phe, and both have another codon.
        let mut starting = Vec::new();
        for &aa in protein {
            let rarest = *ecoli.ranked_codons(aa, standard()).last().unwrap();
            starting.extend_from_slice(&codon_bases(rarest));
        }
        let at = 30;
        starting[at..at + 6].copy_from_slice(b"GAATTC");
        let ecori = EnzymeDb::builtin().get("EcoRI").unwrap();
        assert_eq!(find_sites(&starting, false, &[ecori]).len(), 1);

        let options = OptimizeOptions { avoid_enzymes: vec!["EcoRI".to_string()], max_run: 5, ..OptimizeOptions::default() };
        let result = optimize(&starting, ecoli, standard(), &options).unwrap();
        // The protein is untouched, and the sequence is better adapted than it was.
        assert_eq!(translate(&result.seq, standard(), false), translate(&starting, standard(), false));
        assert_eq!(result.protein, translate(&starting, standard(), false));
        assert!(result.cai_after > result.cai_before, "CAI {} -> {}", result.cai_before, result.cai_after);
        assert!(result.cai_after > 0.8, "CAI after {}", result.cai_after);
        assert!(!result.changes.is_empty());
        // Every reported change is a real synonymous swap at its own codon.
        for change in &result.changes {
            assert_eq!(change.span, Span::new(change.codon * 3, change.codon * 3 + 3));
            assert_eq!(&result.seq[change.span.start..change.span.end], &change.new);
            assert_eq!(&starting[change.span.start..change.span.end], &change.old);
            assert_eq!(standard().translate_codon(&change.old), standard().translate_codon(&change.new));
        }
        // The constraints hold: no EcoRI site, no long runs, GC inside the window range.
        assert!(find_sites(&result.seq, false, &[ecori]).is_empty(), "EcoRI survived");
        let mut run = 1;
        for ii in 1..result.seq.len() {
            run = if result.seq[ii] == result.seq[ii - 1] { run + 1 } else { 1 };
            assert!(run <= 5, "run of {run} at {ii}");
        }
        let (low, high) = options.gc_range.unwrap();
        for window in result.seq.windows(options.gc_window) {
            let gc = alphabet::gc_fraction(window);
            assert!((low..=high).contains(&gc), "window GC {gc}");
        }
        assert!(result.unmet.is_empty(), "{:?}", result.unmet);
        // Matched choice keeps the protein too, and uses more than one codon per amino acid.
        let matched = optimize(&starting, ecoli, standard(), &OptimizeOptions { choice: CodonChoice::Matched, seed: 4, ..options.clone() }).unwrap();
        assert_eq!(translate(&matched.seq, standard(), false), result.protein);
        assert_ne!(matched.seq, result.seq);
    }

    #[test]
    fn optimisation_reports_what_it_cannot_do_and_always_stops() {
        let ecoli = table("ecoli").unwrap();
        // Lysine has only AAA and AAG, so a GC-rich window is out of reach: the search must give
        // up and say so rather than run on.
        let starting = reverse_translate(b"KKKKKKKKKKKKKKKKKKKK", ecoli, standard(), ReverseMode::MostFrequent);
        let options = OptimizeOptions { gc_range: Some((0.90, 1.0)), gc_window: 30, max_run: 3, ..OptimizeOptions::default() };
        let result = optimize(&starting, ecoli, standard(), &options).unwrap();
        assert_eq!(translate(&result.seq, standard(), false), b"KKKKKKKKKKKKKKKKKKKK".to_vec());
        assert!(result.unmet.iter().any(|note| note.contains("GC")), "{:?}", result.unmet);
        // An unknown enzyme name is an error, and so is a sequence that is not whole codons.
        assert!(optimize(&starting, ecoli, standard(), &OptimizeOptions { avoid_enzymes: vec!["NotAnEnzyme".into()], ..Default::default() }).is_err());
        assert!(optimize(b"ATGAA", ecoli, standard(), &OptimizeOptions::default()).is_err());
        assert!(optimize(b"", ecoli, standard(), &OptimizeOptions::default()).is_err());
        // An ambiguous codon is kept and reported.
        let mut ambiguous = starting.clone();
        ambiguous[3..6].copy_from_slice(b"ANA");
        let result = optimize(&ambiguous, ecoli, standard(), &OptimizeOptions::default()).unwrap();
        assert_eq!(&result.seq[3..6], b"ANA");
        assert!(result.unmet.iter().any(|note| note.contains("ambiguous")), "{:?}", result.unmet);
    }

    /// A sequence with a coding region, built from a protein so every codon is known.
    fn coding_sequence(protein: &[u8]) -> Vec<u8> {
        reverse_translate(protein, table("human").unwrap(), standard(), ReverseMode::Weighted { seed: 2024 })
    }

    /// Apply a change's substitutions to a copy of the sequence.
    fn apply(seq: &[u8], change: &SilentChange) -> Vec<u8> {
        let mut out = seq.to_vec();
        for sub in &change.subs {
            assert_eq!(out[sub.pos], sub.from, "substitution at {} expected {}", sub.pos, sub.from as char);
            out[sub.pos] = sub.to;
        }
        out
    }

    /// The protein a coding region translates to.
    fn protein_of(seq: &[u8], cds: Span, reverse: bool) -> Vec<u8> {
        let mut bases = cds.extract(seq);
        if reverse {
            alphabet::revcomp_in_place(&mut bases);
        }
        translate(&bases, standard(), false)
    }

    /// Where one enzyme's sites are, by start and orientation. A non-palindromic site is reported
    /// in both orientations at the same position, so the orientation is part of its identity.
    fn sites_of(seq: &[u8], circular: bool, enzyme: &Enzyme) -> Vec<(usize, u8)> {
        find_sites(seq, circular, &[enzyme]).iter().map(|site| (site.start, site.strand as u8)).collect()
    }

    #[test]
    fn silent_changes_add_and_remove_sites_on_a_forward_coding_sequence() {
        let db = EnzymeDb::builtin();
        let enzymes = db.set("common").unwrap();
        let coding = coding_sequence(b"MGEFLKVTTANDSRQPWYHCIG");
        let mut seq = b"TTGACATAGGCTCA".to_vec();
        let cds = Span::new(seq.len(), seq.len() + coding.len());
        seq.extend_from_slice(&coding);
        seq.extend_from_slice(b"GGTACCAATTGTT");
        let protein = protein_of(&seq, cds, false);

        let changes = silent_site_changes(&seq, false, cds, false, &enzymes, standard(), 2);
        assert!(!changes.is_empty());
        assert!(changes.iter().any(|change| change.kind == SilentKind::Add));
        assert!(changes.iter().any(|change| change.subs.len() == 1));
        // Sorted by substitutions, then side effects.
        for pair in changes.windows(2) {
            let (a, b) = (&pair[0], &pair[1]);
            assert!((a.subs.len(), a.side_effects()) <= (b.subs.len(), b.side_effects()));
        }
        for change in &changes {
            // Every substitution lies inside the coding sequence and changes the base it names.
            assert!(!change.subs.is_empty() && change.subs.len() <= 2);
            for sub in &change.subs {
                assert!(cds.contains(sub.pos, seq.len()), "substitution at {} is outside the CDS", sub.pos);
                assert_ne!(sub.from, sub.to);
                assert_eq!(seq[sub.pos], sub.from);
            }
            // The protein is unchanged, and each rewritten codon keeps its amino acid.
            let edited = apply(&seq, change);
            assert_eq!(protein_of(&edited, cds, false), protein, "{change:?}");
            for codon in &change.codons {
                assert_eq!(standard().translate_codon(&codon.old), codon.amino_acid);
                assert_eq!(standard().translate_codon(&codon.new), codon.amino_acid);
                for offset in 0..3 {
                    assert_eq!(edited[codon.positions[offset]], codon.new[offset]);
                }
            }
            // The site really appears or disappears, where the change says it does.
            let enzyme = db.get(&change.enzyme).unwrap();
            let named = (change.site.start, change.strand as u8);
            let before = sites_of(&seq, false, enzyme);
            let after = sites_of(&edited, false, enzyme);
            match change.kind {
                SilentKind::Add => {
                    assert!(!before.contains(&named) && after.contains(&named), "{change:?}");
                    // The cut position is the one the enzyme database places for that site.
                    let site = find_sites(&edited, false, &[enzyme]).into_iter().find(|site| (site.start, site.strand as u8) == named).unwrap();
                    assert_eq!(change.cut, site.cuts.first().map(|cut| cut.pos), "{change:?}");
                }
                SilentKind::Remove => assert!(before.contains(&named) && !after.contains(&named), "{change:?}"),
            }
        }
    }

    #[test]
    fn silent_changes_report_every_side_effect() {
        // Applying a change must create and destroy exactly the sites it says it does.
        let db = EnzymeDb::builtin();
        let enzymes = db.set("common").unwrap();
        // GCTAGC is both an NheI and a BmtI site, and CCCGGG both a SmaI and an XmaI site, so a
        // change that clears one clears the other as well: exactly the side effect to report.
        let coding = {
            let mut coding = coding_sequence(b"MEFKT");
            coding.extend_from_slice(b"GCTAGC");
            coding.extend_from_slice(&coding_sequence(b"SGGRLL"));
            coding.extend_from_slice(b"CCCGGG");
            coding.extend_from_slice(&coding_sequence(b"DDWQAHNPVCYM"));
            coding
        };
        let mut seq = b"CCGGAATTCGATATCGGTACCGC".to_vec();
        let cds = Span::new(seq.len(), seq.len() + coding.len());
        seq.extend_from_slice(&coding);
        seq.extend_from_slice(b"AAGCTTGCGGCCGCTT");
        let key = |site: &crate::enzymes::Site| (site.enzyme.clone(), site.start, site.strand as u8);
        let before: std::collections::HashSet<_> = find_sites(&seq, false, &enzymes).iter().map(key).collect();
        let changes = silent_site_changes(&seq, false, cds, false, &enzymes, standard(), 2);
        assert!(!changes.is_empty());
        let mut with_effects = 0;
        for change in &changes {
            let edited = apply(&seq, change);
            let after: std::collections::HashSet<_> = find_sites(&edited, false, &enzymes).iter().map(key).collect();
            // Everything the change claims about the main site and the side effects, checked
            // against a full rescan of the whole sequence.
            let mut created: Vec<_> = after.difference(&before).cloned().collect();
            let mut destroyed: Vec<_> = before.difference(&after).cloned().collect();
            let mut claimed_created: Vec<_> = change.also_created.iter().map(|site| (site.enzyme.clone(), site.start, site.strand as u8)).collect();
            let mut claimed_destroyed: Vec<_> = change.also_destroyed.iter().map(|site| (site.enzyme.clone(), site.start, site.strand as u8)).collect();
            match change.kind {
                SilentKind::Add => claimed_created.push((change.enzyme.clone(), change.site.start, change.strand as u8)),
                SilentKind::Remove => claimed_destroyed.push((change.enzyme.clone(), change.site.start, change.strand as u8)),
            }
            created.sort();
            destroyed.sort();
            claimed_created.sort();
            claimed_destroyed.sort();
            assert_eq!(created, claimed_created, "created sites of {change:?}");
            assert_eq!(destroyed, claimed_destroyed, "destroyed sites of {change:?}");
            with_effects += usize::from(change.side_effects() > 0);
        }
        // Clearing one of the shared sites disturbs its twin, which is what the side-effect lists
        // exist for.
        assert!(with_effects > 0, "{} changes, none with side effects", changes.len());
    }

    #[test]
    fn silent_changes_work_on_a_reverse_coding_sequence() {
        let db = EnzymeDb::builtin();
        let enzymes = db.set("common").unwrap();
        let protein = b"MGEFLKVTTANDSRQPWYHCIG";
        let coding = coding_sequence(protein);
        // Put the coding sequence on the bottom strand: the top strand holds its complement.
        let mut seq = b"GTTAACCCTGCA".to_vec();
        let cds = Span::new(seq.len(), seq.len() + coding.len());
        seq.extend_from_slice(&alphabet::revcomp(&coding));
        seq.extend_from_slice(b"TTCGAACCGGTT");
        assert_eq!(protein_of(&seq, cds, true), translate(&coding, standard(), false));
        let expected = protein_of(&seq, cds, true);

        let changes = silent_site_changes(&seq, false, cds, true, &enzymes, standard(), 2);
        assert!(!changes.is_empty());
        for change in &changes {
            let edited = apply(&seq, change);
            // Positions are top-strand even though the coding sequence reads the other way.
            for sub in &change.subs {
                assert!(cds.contains(sub.pos, seq.len()));
            }
            // Each codon's positions count down the top strand for a reverse coding sequence.
            for codon in &change.codons {
                assert!(codon.positions[0] > codon.positions[1] && codon.positions[1] > codon.positions[2]);
                for offset in 0..3 {
                    assert_eq!(alphabet::complement_base(edited[codon.positions[offset]]), codon.new[offset]);
                }
            }
            assert_eq!(protein_of(&edited, cds, true), expected, "{change:?}");
            let enzyme = db.get(&change.enzyme).unwrap();
            let after = sites_of(&edited, false, enzyme);
            let named = (change.site.start, change.strand as u8);
            match change.kind {
                SilentKind::Add => assert!(after.contains(&named), "{change:?}"),
                SilentKind::Remove => assert!(!after.contains(&named), "{change:?}"),
            }
        }
        assert!(changes.iter().any(|change| change.kind == SilentKind::Add));
    }

    #[test]
    fn silent_changes_work_through_the_origin_of_a_circular_sequence() {
        let db = EnzymeDb::builtin();
        // The whole table, so sites through the origin are certain to turn up.
        let enzymes: Vec<&Enzyme> = db.enzymes.iter().collect();
        let protein = b"MGEFLKVTTANDSRQPWYHCIGEM";
        let coding = coding_sequence(protein);
        // Rotate the plasmid so the coding sequence straddles the origin.
        let tail = 15usize;
        let mut seq = coding[coding.len() - tail..].to_vec();
        seq.extend_from_slice(b"GACCTGAGGATCCTTAGGCA");
        seq.extend_from_slice(&coding[..coding.len() - tail]);
        let cds = Span::new(seq.len() - (coding.len() - tail), seq.len() + tail);
        assert!(cds.wraps(seq.len()));
        assert_eq!(cds.extract(&seq), coding);
        let expected = protein_of(&seq, cds, false);

        let changes = silent_site_changes(&seq, true, cds, false, &enzymes, standard(), 2);
        assert!(!changes.is_empty());
        let mut wrapping = 0;
        // Every change through the origin, plus a sample of the rest: the whole enzyme table
        // yields thousands of changes, and rescanning for each would dominate the run time.
        for (ii, change) in changes.iter().enumerate() {
            let wraps = change.site.wraps(seq.len());
            wrapping += usize::from(wraps);
            if !wraps && !ii.is_multiple_of(29) {
                continue;
            }
            let edited = apply(&seq, change);
            assert_eq!(protein_of(&edited, cds, false), expected, "{change:?}");
            let enzyme = db.get(&change.enzyme).unwrap();
            let sites = find_sites(&edited, true, &[enzyme]);
            let named = (change.site.start, change.strand as u8);
            let present: Vec<(usize, u8)> = sites.iter().map(|site| (site.start, site.strand as u8)).collect();
            match change.kind {
                SilentKind::Add => {
                    assert!(present.contains(&named), "{change:?}");
                    let site = sites.iter().find(|site| (site.start, site.strand as u8) == named).unwrap();
                    assert_eq!(change.cut, site.cuts.first().map(|cut| cut.pos), "{change:?}");
                }
                SilentKind::Remove => assert!(!present.contains(&named), "{change:?}"),
            }
        }
        // Changes near the origin create sites that run through it.
        assert!(wrapping > 0, "no site through the origin among {} changes", changes.len());
    }

    #[test]
    fn silent_changes_remove_a_planted_site() {
        let db = EnzymeDb::builtin();
        let ecori = db.get("EcoRI").unwrap();
        let bamhi = db.get("BamHI").unwrap();
        let enzymes = vec![ecori, bamhi];
        // GAA TTC codes for Glu-Phe and GGA TCC for Gly-Ser, so both sites can go silently.
        let coding = {
            let mut coding = coding_sequence(b"MKVTT");
            coding.extend_from_slice(b"GAATTC");
            coding.extend_from_slice(&coding_sequence(b"ADLNP"));
            coding.extend_from_slice(b"GGATCC");
            coding.extend_from_slice(&coding_sequence(b"QHRGM"));
            coding
        };
        let mut seq = b"CCTTAATTAACC".to_vec();
        let cds = Span::new(seq.len(), seq.len() + coding.len());
        seq.extend_from_slice(&coding);
        seq.extend_from_slice(b"TTAATTAAGG");
        let expected = protein_of(&seq, cds, false);
        assert_eq!(find_sites(&seq, false, &enzymes).len(), 2);

        let changes = silent_site_changes(&seq, false, cds, false, &enzymes, standard(), 1);
        let removals: Vec<&SilentChange> = changes.iter().filter(|change| change.kind == SilentKind::Remove).collect();
        assert!(removals.iter().any(|change| change.enzyme == "EcoRI"), "{changes:?}");
        assert!(removals.iter().any(|change| change.enzyme == "BamHI"), "{changes:?}");
        for change in removals {
            assert_eq!(change.subs.len(), 1);
            let edited = apply(&seq, change);
            assert_eq!(protein_of(&edited, cds, false), expected);
            let enzyme = db.get(&change.enzyme).unwrap();
            assert!(!sites_of(&edited, false, enzyme).contains(&(change.site.start, change.strand as u8)));
        }
        // With one substitution allowed, no change asks for two.
        assert!(changes.iter().all(|change| change.subs.len() == 1));
        // Two substitutions open up sites one cannot reach.
        let one = silent_site_changes(&seq, false, cds, false, &db.set("common").unwrap(), standard(), 1).len();
        let two = silent_site_changes(&seq, false, cds, false, &db.set("common").unwrap(), standard(), 2).len();
        assert!(two > one, "{one} changes with one substitution, {two} with two");
    }

    #[test]
    fn the_window_scanner_agrees_with_find_sites() {
        // The silent-change search uses its own scanner so it can rescan a window cheaply. It
        // must see exactly what the enzyme module's search sees.
        let enzymes: Vec<&Enzyme> = EnzymeDb::builtin().enzymes.iter().collect();
        let mut state = 7u64;
        for round in 0..4 {
            let len = 300 + round * 200;
            let seq: Vec<u8> = (0..len)
                .map(|_| {
                    let value = (next_random(&mut state) >> 20) as usize % 100;
                    if value == 0 {
                        b'N'
                    } else if value == 1 {
                        b'R'
                    } else {
                        b"ACGT"[value % 4]
                    }
                })
                .collect();
            for circular in [true, false] {
                let index = SiteIndex::new(&enzymes, seq.len(), circular);
                let mut found = Vec::new();
                index.scan(&seq, circular, 0, seq.len(), &mut found);
                let mut mine: Vec<(String, usize, u8, Option<usize>)> = found
                    .iter()
                    .map(|&site| {
                        let described = index.site_ref(site);
                        (described.enzyme, described.start, described.strand as u8, index.cut_position(site, seq.len(), circular))
                    })
                    .collect();
                let mut theirs: Vec<(String, usize, u8, Option<usize>)> = find_sites(&seq, circular, &enzymes)
                    .into_iter()
                    .map(|site| (site.enzyme, site.start, site.strand as u8, site.cuts.first().map(|cut| cut.pos)))
                    .collect();
                mine.sort();
                theirs.sort();
                assert_eq!(mine, theirs, "length {len}, circular {circular}");
            }
        }
    }

    #[test]
    fn silent_changes_handle_a_kilobase_against_every_enzyme() {
        // The search a person waits for: a 1 kb coding sequence against the whole enzyme table.
        let enzymes: Vec<&Enzyme> = EnzymeDb::builtin().enzymes.iter().collect();
        let protein: Vec<u8> = b"MKVLATTGSLNDIEWQRPHFYACGSTVEDKLMNQR".iter().cycle().take(333).copied().collect();
        let coding = coding_sequence(&protein);
        assert_eq!(coding.len(), 999);
        let mut seq = b"TTGACAGCTAGCATG".to_vec();
        let cds = Span::new(seq.len(), seq.len() + coding.len());
        seq.extend_from_slice(&coding);
        seq.extend_from_slice(b"TAAGGATCCTTGCA");
        let expected = protein_of(&seq, cds, false);
        // The three searches a person waits for, timed: one substitution and two against the
        // whole table, and two against the enzymes a cloning window usually offers.
        let common = EnzymeDb::builtin().set("common").unwrap();
        for (label, set, max_changes) in [("all/1", &enzymes, 1usize), ("all/2", &enzymes, 2), ("common/2", &common, 2)] {
            let started = std::time::Instant::now();
            let found = silent_site_changes(&seq, false, cds, false, set, standard(), max_changes);
            eprintln!("1 kb CDS, {} enzymes, up to {max_changes} substitutions ({label}): {} changes in {:?}", set.len(), found.len(), started.elapsed());
        }
        let changes = silent_site_changes(&seq, false, cds, false, &enzymes, standard(), 2);
        assert!(changes.len() > 100, "only {} changes", changes.len());
        // Spot-check a sample, since checking every change here would dominate the run time.
        for change in changes.iter().step_by(37) {
            let edited = apply(&seq, change);
            assert_eq!(protein_of(&edited, cds, false), expected, "{change:?}");
            let enzyme = EnzymeDb::builtin().get(&change.enzyme).unwrap();
            let present = sites_of(&edited, false, enzyme);
            let named = (change.site.start, change.strand as u8);
            match change.kind {
                SilentKind::Add => assert!(present.contains(&named), "{change:?}"),
                SilentKind::Remove => assert!(!present.contains(&named), "{change:?}"),
            }
        }
    }
}
