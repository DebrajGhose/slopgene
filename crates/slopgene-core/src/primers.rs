//! Primers: melting temperature, binding sites and quality checks, by Primer3's methods.
//!
//! Melting temperatures use the nearest-neighbour model with the unified parameters of
//! SantaLucia (1998), the entropy salt correction of the same paper and the sodium equivalent of
//! magnesium from von Ahsen et al. (2001), as Primer3 computes them (its `PRIMER_TM_FORMULA` 1
//! and `PRIMER_SALT_CORRECTIONS` 1). The default conditions are Primer3's: 50 mM monovalent
//! salt, 1.5 mM Mg2+, 0.6 mM dNTPs and 50 nM oligo. The quality checks are Primer3's too: its
//! alignment scores for a primer pairing with itself or with another primer, and its stability
//! of the last five bases. `primer_design` builds the designer on them.

use crate::alphabet::{self, iupac_matches};
use crate::model::{BindingSite, Document, Strand, XmlExtra};
use crate::span::Span;

/// Gas constant in cal/(K mol), rounded as Primer3 rounds it.
const R: f64 = 1.987;

/// Conditions for a melting-temperature estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TmConditions {
    /// Monovalent cations (Na+, K+), millimolar.
    pub monovalent_mm: f64,
    /// Mg2+, millimolar. Converted to a sodium equivalent (von Ahsen et al. 2001).
    pub magnesium_mm: f64,
    /// dNTPs, millimolar (they bind Mg2+).
    pub dntp_mm: f64,
    /// Oligo concentration, nanomolar (Primer3's PRIMER_DNA_CONC).
    pub oligo_nm: f64,
}

impl Default for TmConditions {
    /// Primer3's defaults: 50 mM monovalent salt, 1.5 mM Mg2+, 0.6 mM dNTPs and 50 nM oligo.
    fn default() -> Self {
        TmConditions { monovalent_mm: 50.0, magnesium_mm: 1.5, dntp_mm: 0.6, oligo_nm: 50.0 }
    }
}

/// Nearest-neighbour enthalpy (kcal/mol) and entropy (cal/K/mol) for the dinucleotide `a b`
/// read 5'->3' on one strand (SantaLucia 1998, Table 1).
fn nearest_neighbour(a: u8, b: u8) -> Option<(f64, f64)> {
    let pair = [a.to_ascii_uppercase(), b.to_ascii_uppercase()];
    // Each stack appears twice, once per strand: AA/TT is the same stack as TT/AA, and so on.
    let value = match &pair {
        b"AA" | b"TT" => (-7.9, -22.2),
        b"AT" => (-7.2, -20.4),
        b"TA" => (-7.2, -21.3),
        b"CA" | b"TG" => (-8.5, -22.7),
        b"GT" | b"AC" => (-8.4, -22.4),
        b"CT" | b"AG" => (-7.8, -21.0),
        b"GA" | b"TC" => (-8.2, -22.2),
        b"CG" => (-10.6, -27.2),
        b"GC" => (-9.8, -24.4),
        b"GG" | b"CC" => (-8.0, -19.9),
        _ => return None,
    };
    Some(value)
}

/// Duplex enthalpy and entropy (1 M NaCl) for a perfectly matched duplex of `seq`.
/// Returns `None` when the sequence has ambiguous bases or is shorter than 2.
pub fn duplex_thermodynamics(seq: &[u8]) -> Option<(f64, f64)> {
    if seq.len() < 2 {
        return None;
    }
    let mut dh = 0.0;
    let mut ds = 0.0;
    for pair in seq.windows(2) {
        let (h, s) = nearest_neighbour(pair[0], pair[1])?;
        dh += h;
        ds += s;
    }
    // Initiation: one term per duplex end, depending on whether the end is G/C or A/T.
    for &end in [seq[0], seq[seq.len() - 1]].iter() {
        match end.to_ascii_uppercase() {
            b'G' | b'C' => {
                dh += 0.1;
                ds += -2.8;
            }
            _ => {
                dh += 2.3;
                ds += 4.1;
            }
        }
    }
    // Symmetry correction for self-complementary duplexes.
    if alphabet::revcomp(seq).eq_ignore_ascii_case(seq) {
        ds += -1.4;
    }
    Some((dh, ds))
}

/// Melting temperature (deg C) of `seq` annealed to its perfect complement, as Primer3
/// computes it.
///
/// Falls back to the Wallace rule (2 C per A/T, 4 C per G/C) for sequences with ambiguous bases
/// or shorter than 8 bases, where nearest-neighbour values are unreliable.
pub fn melting_temperature(seq: &[u8], conditions: &TmConditions) -> f64 {
    let wallace = || {
        let gc = seq.iter().filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C' | b'S')).count() as f64;
        let at = seq.len() as f64 - gc;
        2.0 * at + 4.0 * gc
    };
    if seq.len() < 8 {
        return wallace();
    }
    let Some((dh, ds)) = duplex_thermodynamics(seq) else {
        return wallace();
    };
    // The sodium equivalent (mM) of the salts: the monovalent cations plus 120 times the square
    // root of the Mg2+ left after the dNTPs have bound their share (von Ahsen et al. 2001).
    let free_mg = (conditions.magnesium_mm - conditions.dntp_mm).max(0.0);
    let k_mm = conditions.monovalent_mm + 120.0 * free_mg.sqrt();
    let ds_salt = ds + 0.368 * (seq.len() as f64 - 1.0) * (k_mm / 1000.0).max(1e-6).ln();
    // The concentration term as Primer3 takes it: the oligo concentration for an oligo that is
    // its own complement, and a quarter of it otherwise, for two strands at the same
    // concentration.
    let oligo = conditions.oligo_nm.max(1e-3) * 1e-9;
    let c = if alphabet::revcomp(seq).eq_ignore_ascii_case(seq) { oligo } else { oligo / 4.0 };
    dh * 1000.0 / (ds_salt + R * c.ln()) - 273.15
}

/// The stability of the stack `a b` (read 5'->3' on one strand) at 37 °C, as -ΔG° in kcal/mol
/// from SantaLucia (1998), the table Primer3 uses. A stack with an ambiguous base takes the
/// value Primer3 gives it.
fn stacking_stability(a: u8, b: u8) -> f64 {
    let code = |x: u8| match x.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    };
    const TABLE: [[f64; 5]; 5] = [
        [1.00, 1.44, 1.28, 0.88, 0.88], // AA, AC, AG, AT, AN
        [1.45, 1.84, 2.17, 1.28, 1.45], // CA, CC, CG, CT, CN
        [1.30, 2.24, 1.84, 1.44, 1.30], // GA, GC, GG, GT, GN
        [0.58, 1.30, 1.45, 1.00, 0.58], // TA, TC, TG, TT, TN
        [0.58, 1.30, 1.28, 0.88, 0.58], // NA, NC, NG, NT, NN
    ];
    TABLE[code(a)][code(b)]
}

/// Primer3's stability of a duplex of `seq`, in kcal/mol: -ΔG°37 from SantaLucia's (1998)
/// nearest-neighbour free energies, that is the stacks, less 1.96 for starting the duplex, 0.05
/// for each A or T at an end and 0.43 for an oligo that is its own complement. Bigger is more
/// stable.
pub fn duplex_stability(seq: &[u8]) -> f64 {
    if seq.len() < 2 {
        return 0.0;
    }
    let mut stability: f64 = seq.windows(2).map(|pair| stacking_stability(pair[0], pair[1])).sum();
    stability -= 1.96;
    for &end in &[seq[0], seq[seq.len() - 1]] {
        if matches!(end.to_ascii_uppercase(), b'A' | b'T') {
            stability -= 0.05;
        }
    }
    if alphabet::revcomp(seq).eq_ignore_ascii_case(seq) {
        stability -= 0.43;
    }
    stability
}

/// The stability of the last five bases, Primer3's end stability: how firmly a primer's 3' end
/// holds on, and so how readily a polymerase starts from a place where only the end pairs.
/// GCGCG gives 6.86 kcal/mol, the most possible, and TATAT 0.86.
pub fn end_stability(seq: &[u8]) -> f64 {
    duplex_stability(&seq[seq.len().saturating_sub(5)..])
}

/// How Primer3's alignment of two sequences is anchored.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Anchor {
    /// The best stretch anywhere, as in a local alignment.
    Anywhere,
    /// An alignment that reaches the last base of the first sequence, Primer3's test for a 3' end
    /// that a polymerase could extend.
    End,
}

/// Primer3's score of one base of a sequence against one base of another, in hundredths: 100 for
/// the same base, -100 for different bases and -25 when either is N or another ambiguity code.
/// The second sequence is the reverse complement of the first one's partner, so the same base
/// means that the two pair.
fn align_score(a: u8, b: u8) -> i32 {
    let (a, b) = (a.to_ascii_uppercase(), b.to_ascii_uppercase());
    if !matches!(a, b'A' | b'C' | b'G' | b'T') || !matches!(b, b'A' | b'C' | b'G' | b'T') {
        -25
    } else if a == b {
        100
    } else {
        -100
    }
}

/// Primer3's alignment score between `s1` and `s2`, as its `dpal` routine computes it with the
/// settings Primer3 uses for primers: a match scores 1.00, a mismatch -1.00, a base against an N
/// -0.25 and a gap of one base -2.00, with no longer gaps. With [`Anchor::Anywhere`] the best
/// stretch anywhere counts and stretches below zero are dropped, as in a local alignment. With
/// [`Anchor::End`] the alignment must reach the last base of `s1` and runs back to the start of
/// one of the sequences, mismatches and all. A score below zero is reported as zero.
pub fn primer3_alignment(s1: &[u8], s2: &[u8], anchor: Anchor) -> f64 {
    let x = alphabet::upper(s1);
    let y = alphabet::upper(s2);
    let (xlen, ylen) = (x.len(), y.len());
    if xlen == 0 || ylen == 0 {
        return 0.0;
    }
    const GAP: i32 = -200;
    let local = anchor == Anchor::Anywhere;
    let mut s = vec![vec![0i32; ylen]; xlen];
    let mut best_score = i32::MIN;
    // The first column and the first row: one pair each, with nothing before it. A local
    // alignment may end anywhere, so each cell is a candidate; an anchored one ends in the last
    // row only, and starts from the last row's first cell.
    for ii in 0..xlen {
        let mut score = align_score(x[ii], y[0]);
        if local {
            score = score.max(0);
            best_score = best_score.max(score);
        }
        s[ii][0] = score;
    }
    if !local {
        best_score = s[xlen - 1][0];
    }
    for jj in 0..ylen {
        let mut score = align_score(x[0], y[jj]);
        if local {
            score = score.max(0);
            best_score = best_score.max(score);
        }
        s[0][jj] = score;
    }
    // Each further cell continues the best of the cell diagonally before it and the two cells a
    // one-base gap away, and adds the score of its own pair.
    for ii in 1..xlen {
        for jj in 1..ylen {
            let mut before = s[ii - 1][jj - 1];
            if ii > 1 {
                before = before.max(s[ii - 2][jj - 1] + GAP);
            }
            if jj > 1 {
                before = before.max(s[ii - 1][jj - 2] + GAP);
            }
            let mut score = before + align_score(x[ii], y[jj]);
            if local || ii == xlen - 1 {
                best_score = best_score.max(score);
            }
            if local {
                score = score.max(0);
            }
            s[ii][jj] = score;
        }
    }
    f64::from(best_score.max(0)) / 100.0
}

/// Primer3's score for an oligo pairing with another copy of itself anywhere (its SELF_ANY).
/// Primer3 accepts up to 8.
pub fn self_any(oligo: &[u8]) -> f64 {
    primer3_alignment(oligo, &alphabet::revcomp(oligo), Anchor::Anywhere)
}

/// Primer3's score for an oligo's 3' end pairing with another copy of itself (its SELF_END), the
/// primer-dimer test. Primer3 accepts up to 3.
pub fn self_end(oligo: &[u8]) -> f64 {
    primer3_alignment(oligo, &alphabet::revcomp(oligo), Anchor::End)
}

/// Primer3's score for two primers pairing with each other anywhere (its PAIR_COMPL_ANY).
pub fn pair_any(a: &[u8], b: &[u8]) -> f64 {
    primer3_alignment(a, &alphabet::revcomp(b), Anchor::Anywhere)
}

/// Primer3's score for the 3' end of either primer pairing with the other (its PAIR_COMPL_END).
pub fn pair_end(a: &[u8], b: &[u8]) -> f64 {
    primer3_alignment(a, &alphabet::revcomp(b), Anchor::End).max(primer3_alignment(b, &alphabet::revcomp(a), Anchor::End))
}

/// Rules for finding where a primer anneals, after SnapGene's hybridization parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BindingRules {
    /// The primer's 3' end must match the template for at least this many contiguous bases.
    pub min_continuous_match: usize,
    /// Whether mismatches are tolerated upstream of that 3' stretch.
    pub allow_mismatch: bool,
    /// The annealed part must reach this melting temperature (deg C).
    pub min_tm: f64,
    pub conditions: TmConditions,
}

impl Default for BindingRules {
    fn default() -> Self {
        // SnapGene's defaults: 10 contiguous 3' bases, mismatches allowed, Tm at least 40 C.
        BindingRules { min_continuous_match: 10, allow_mismatch: true, min_tm: 40.0, conditions: TmConditions::default() }
    }
}

/// A primer binding site found on a template.
#[derive(Clone, Debug, PartialEq)]
pub struct Binding {
    /// Annealed bases on the template, top-strand coordinates.
    pub span: Span,
    /// `Forward`: the primer reads like the top strand. `Reverse`: like the bottom strand.
    pub strand: Strand,
    /// Number of primer bases annealed, counted from the 3' end.
    pub annealed: usize,
    pub mismatches: usize,
    pub tm: f64,
}

/// Find every site where `primer` (5'->3') anneals to `template`.
pub fn find_binding_sites(template: &[u8], circular: bool, primer: &[u8], rules: &BindingRules) -> Vec<Binding> {
    let n = primer.len();
    let seq_len = template.len();
    let mut found = Vec::new();
    if n < rules.min_continuous_match.max(1) || seq_len == 0 {
        return found;
    }
    let primer = alphabet::upper(primer);
    let top = alphabet::upper(template);
    let base_at = |pos: i64| -> Option<u8> {
        if circular {
            Some(top[pos.rem_euclid(seq_len as i64) as usize])
        } else if pos >= 0 && (pos as usize) < seq_len {
            Some(top[pos as usize])
        } else {
            None
        }
    };

    // Strand ii = 0: the primer reads like the top strand, so primer base k pairs with template
    // position (end - n + k), and the 3' end sits at the right. Strand ii = 1: the primer reads
    // like the bottom strand; compare its reverse complement with the top strand, 3' end at
    // the left.
    let rc = alphabet::revcomp(&primer);
    for ii in 0..2 {
        let strand = if ii == 0 { Strand::Forward } else { Strand::Reverse };
        let probe = if ii == 0 { &primer } else { &rc };
        let positions = if circular { seq_len } else { seq_len + n };
        for anchor in 0..positions as i64 {
            // `anchor` is the template position paired with the primer's 3'-most base.
            let anchored = (0..rules.min_continuous_match).all(|iii| {
                let (probe_index, template_pos) = if ii == 0 { (n - 1 - iii, anchor - iii as i64) } else { (iii, anchor + iii as i64) };
                base_at(template_pos).is_some_and(|b| iupac_matches(b, probe[probe_index]) || iupac_matches(probe[probe_index], b))
            });
            if !anchored {
                continue;
            }
            // Extend toward the primer's 5' end. A mismatch is tolerated only when the next
            // three bases match again; otherwise annealing stops there.
            let mut annealed = rules.min_continuous_match;
            let mut mismatches = 0;
            let mut iii = rules.min_continuous_match;
            while iii < n {
                let pairs = |offset: usize| -> bool {
                    let k = iii + offset;
                    if k >= n {
                        return false;
                    }
                    let (probe_index, template_pos) = if ii == 0 { (n - 1 - k, anchor - k as i64) } else { (k, anchor + k as i64) };
                    base_at(template_pos).is_some_and(|b| iupac_matches(b, probe[probe_index]) || iupac_matches(probe[probe_index], b))
                };
                if pairs(0) {
                    iii += 1;
                    annealed = iii;
                } else if rules.allow_mismatch && pairs(1) && pairs(2) && pairs(3) {
                    mismatches += 1;
                    iii += 1;
                } else {
                    break;
                }
            }
            // Template span of the annealed bases.
            let (start, end) = if ii == 0 { (anchor - annealed as i64 + 1, anchor + 1) } else { (anchor, anchor + annealed as i64) };
            if !circular && (start < 0 || end > seq_len as i64) {
                continue;
            }
            if circular && annealed > seq_len {
                continue;
            }
            let annealed_part: Vec<u8> = if ii == 0 { primer[n - annealed..].to_vec() } else { rc[..annealed].to_vec() };
            let tm = melting_temperature(&annealed_part, &rules.conditions);
            if tm < rules.min_tm {
                continue;
            }
            let start = start.rem_euclid(seq_len as i64) as usize;
            found.push(Binding { span: Span::new(start, start + annealed), strand, annealed, mismatches, tm });
        }
    }
    found.sort_by_key(|b| (b.span.start, b.strand as u8));
    found
}

/// Recompute the binding sites of every primer in a document (after editing the sequence).
pub fn refresh_binding_sites(doc: &mut Document, rules: &BindingRules) {
    let circular = doc.is_circular();
    for primer in &mut doc.primers {
        primer.sites = find_binding_sites(&doc.seq, circular, primer.seq.as_bytes(), rules)
            .into_iter()
            .map(|b| BindingSite { span: b.span, strand: b.strand, tm: Some(b.tm), extra: XmlExtra::default() })
            .collect();
    }
}

/// Basic properties and warning signs of a primer.
#[derive(Clone, Debug, PartialEq)]
pub struct PrimerReport {
    pub length: usize,
    pub gc_percent: f64,
    pub tm: f64,
    /// Primer3's alignment scores of the oligo against its own reverse complement: anywhere, and
    /// ending at its 3' end. Primer3 accepts up to 8 and 3.
    pub self_any: f64,
    pub self_end: f64,
    /// Primer3's stability of the last five bases, in kcal/mol.
    pub end_stability: f64,
    /// Longest stretch of the primer that can pair with another copy of itself.
    pub self_dimer: usize,
    /// Longest such stretch that includes the primer's 3'-terminal base (the risky kind).
    pub three_prime_dimer: usize,
    /// Longest hairpin stem (at least 3 unpaired bases in the loop).
    pub hairpin_stem: usize,
    /// Longest run of a single base.
    pub longest_run: usize,
    /// G or C bases among the last five at the 3' end.
    pub three_prime_gc: usize,
    pub warnings: Vec<String>,
}

/// Analyse a primer: Primer3's numbers and the limits of its that the primer breaks, with the
/// plain stretch counts (self-dimer, hairpin stem, run) as further information.
pub fn analyze(primer: &[u8], conditions: &TmConditions) -> PrimerReport {
    let seq = alphabet::upper(primer);
    let n = seq.len();
    let rc = alphabet::revcomp(&seq);
    let gc_percent = alphabet::gc_fraction(&seq) * 100.0;
    let tm = melting_temperature(&seq, conditions);

    // Self-dimer: longest common stretch between the primer and its reverse complement, at any
    // offset. A stretch that reaches the primer's last base can be extended by a polymerase.
    let mut self_dimer = 0;
    let mut three_prime_dimer = 0;
    for ii in 0..n {
        for iii in 0..n {
            let mut run = 0;
            while ii + run < n && iii + run < n && seq[ii + run] == rc[iii + run] {
                run += 1;
            }
            if run > 0 {
                self_dimer = self_dimer.max(run);
                if ii + run == n {
                    three_prime_dimer = three_prime_dimer.max(run);
                }
            }
        }
    }

    // Hairpin: a stem of paired bases i..i+s pairing with j-s+1..j, loop of at least 3.
    let mut hairpin_stem = 0;
    for ii in 0..n {
        for iii in (ii + 4..n).rev() {
            let mut stem = 0;
            while ii + stem < iii.saturating_sub(stem) && iii >= stem && iii - stem > ii + stem + 2 {
                let a = seq[ii + stem];
                let b = seq[iii - stem];
                if alphabet::complement_base(a) == b {
                    stem += 1;
                } else {
                    break;
                }
            }
            hairpin_stem = hairpin_stem.max(stem);
        }
    }

    let mut longest_run = 0;
    let mut run = 0;
    for ii in 0..n {
        run = if ii > 0 && seq[ii] == seq[ii - 1] { run + 1 } else { 1 };
        longest_run = longest_run.max(run);
    }
    let three_prime_gc = seq[n.saturating_sub(5)..].iter().filter(|b| matches!(b, b'G' | b'C')).count();

    let self_any = self_any(&seq);
    let self_end = self_end(&seq);
    let end_stability = end_stability(&seq);

    // The warnings are Primer3's limits, broken.
    let mut warnings = Vec::new();
    if n < 15 {
        warnings.push(format!("short primer ({n} nt)"));
    }
    if !(20.0..=80.0).contains(&gc_percent) {
        warnings.push(format!("GC content {gc_percent:.0}% is outside 20-80%"));
    }
    if self_any > 8.0 {
        warnings.push(format!("pairs with itself (score {self_any:.1}, limit 8)"));
    }
    if self_end > 3.0 {
        warnings.push(format!("its 3' end pairs with itself (score {self_end:.1}, limit 3)"));
    }
    if longest_run > 5 {
        warnings.push(format!("run of {longest_run} identical bases, more than 5"));
    }
    PrimerReport {
        length: n,
        gc_percent,
        tm,
        self_any,
        self_end,
        end_stability,
        self_dimer,
        three_prime_dimer,
        hairpin_stem,
        longest_run,
        three_prime_gc,
        warnings,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tm_is_in_the_expected_range() {
        let c = TmConditions::default();
        // M13 forward (-20); SnapGene reports about 52-53 C for it under its defaults.
        let tm = melting_temperature(b"GTAAAACGACGGCCAGT", &c);
        assert!((48.0..58.0).contains(&tm), "M13 fwd Tm {tm}");
        // A GC-rich 20-mer melts higher than an AT-rich one.
        let gc = melting_temperature(b"GCGGCCGCGGCGCCGCGGCC", &c);
        let at = melting_temperature(b"ATATTAAATATTTAATATAA", &c);
        assert!(gc > at + 20.0, "gc {gc} at {at}");
        // More salt stabilises the duplex.
        let salty = melting_temperature(b"GTAAAACGACGGCCAGT", &TmConditions { monovalent_mm: 200.0, ..c });
        assert!(salty > tm);
        // Short oligos fall back to the Wallace rule.
        assert_eq!(melting_temperature(b"ACGT", &c), 12.0);
    }

    #[test]
    fn finds_forward_and_reverse_sites_with_tails() {
        let template = b"TTTTTTTTTTGACTGACTGACTGGCATGCATGCATGCAAAAAAAAAACAGTCGATTGCTAACGTGTCAGTACCCC";
        let rules = BindingRules { min_tm: 20.0, ..BindingRules::default() };
        // Forward primer with a 5' tail of G's. The last tail G happens to pair with the G at
        // template position 10, so 18 bases anneal and the other six hang off.
        let fwd = b"GGGGGGGACTGACTGACTGGCATG";
        let sites = find_binding_sites(template, false, fwd, &rules);
        assert_eq!(sites.len(), 1, "{sites:?}");
        assert_eq!(sites[0].strand, Strand::Forward);
        assert_eq!(sites[0].annealed, 18);
        assert_eq!(sites[0].span, Span::new(10, 28));
        // Reverse primer: reverse complement of a template stretch.
        let rev = alphabet::revcomp(b"CAGTCGATTGCTAACGTGTCAGTA");
        let sites = find_binding_sites(template, false, &rev, &rules);
        assert_eq!(sites.len(), 1, "{sites:?}");
        assert_eq!(sites[0].strand, Strand::Reverse);
        assert_eq!(sites[0].span.len(), 24);
    }

    #[test]
    fn single_mismatch_is_tolerated() {
        let template = b"CCCCCCCCCCATGGCTAGCTAGGATCCGATCGATCGTTTTTTTTTT";
        let mut primer = b"ATGGCTAGCTAGGATCCGATCGATCG".to_vec();
        primer[3] = b'C'; // mismatch near the 5' end
        let sites = find_binding_sites(template, false, &primer, &BindingRules::default());
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].mismatches, 1);
        assert_eq!(sites[0].annealed, primer.len());
    }

    #[test]
    fn primer_report_flags_problems() {
        let c = TmConditions::default();
        let good = analyze(b"AGCGGATAACAATTTCACACAGGA", &c);
        assert!(good.warnings.is_empty(), "{:?}", good.warnings);
        let dimer = analyze(b"ATGCGAATTCGCAT", &c);
        assert!(dimer.three_prime_dimer >= 4);
        assert!(dimer.warnings.iter().any(|w| w.starts_with("its 3' end pairs with itself (score 14.0")), "{:?}", dimer.warnings);
    }

    // The reference values below come from Primer3 itself (primer3_core 2.6 through primer3-py
    // 2.4.0), at Primer3's default conditions unless said otherwise.

    #[test]
    fn melting_temperatures_match_primer3() {
        let c = TmConditions::default();
        for (seq, tm) in [
            ("CGTGACGTGACGGACT", 55.973),
            ("GTAAAACGACGGCCAGT", 54.695),
            ("AGCGGATAACAATTTCACACAGGA", 60.806),
            ("ATGCGAATTCGCAT", 47.346),
            ("GCGCGCGCGCGCGCGCGCGC", 83.405),
            ("TTTTAAAATTTTAAAATTTT", 39.100),
            ("CAGGAAACAGCTATGACCATG", 56.676),
            ("AAAAAGCCCCCTTTTTGGGGG", 61.610),
        ] {
            let ours = melting_temperature(seq.as_bytes(), &c);
            assert!((ours - tm).abs() < 0.01, "{seq}: {ours} vs Primer3's {tm}");
        }
        // The example in Primer3's manual: 50 mM salt alone and 50 nM oligo give 50.55 °C.
        let plain = TmConditions { magnesium_mm: 0.0, dntp_mm: 0.0, ..c };
        assert!((melting_temperature(b"CGTGACGTGACGGACT", &plain) - 50.55).abs() < 0.01);
    }

    #[test]
    fn pairing_scores_match_primer3() {
        for (seq, any, end) in [
            ("GTAAAACGACGGCCAGT", 4.0, 4.0),
            ("AGCGGATAACAATTTCACACAGGA", 4.0, 0.0),
            ("ATGCGAATTCGCAT", 14.0, 14.0),
            ("CAGGAAACAGCTATGACCATG", 4.0, 4.0),
            ("GCGCGCGCGCGCGCGCGCGC", 20.0, 20.0),
            ("TTTTAAAATTTTAAAATTTT", 16.0, 16.0),
            ("AAAAAGCCCCCTTTTTGGGGG", 8.0, 5.0),
            ("CGTGACGTGACGGACT", 8.0, 1.0),
            ("ACGTACGTTGCAGGCATGCAAGCT", 10.0, 8.0),
            ("AGCTTGCATGCCTGCAGGTCGACT", 8.0, 4.0),
            // With tails, which Primer3 calls overhangs.
            ("TAAGCAGAATTCGTAAAACGACGGCCAGT", 6.0, 4.0),
            ("TAAGCAGGATCCAGCGGATAACAATTTCACACAGGA", 8.0, 3.0),
        ] {
            assert_eq!(self_any(seq.as_bytes()), any, "{seq} self_any");
            assert_eq!(self_end(seq.as_bytes()), end, "{seq} self_end");
        }
        for (a, b, any, end) in [
            ("GTAAAACGACGGCCAGT", "AGCGGATAACAATTTCACACAGGA", 4.0, 2.0),
            ("ATGCGAATTCGCAT", "CAGGAAACAGCTATGACCATG", 3.0, 3.0),
            ("GCGCGCGCGCGCGCGCGCGC", "TTTTAAAATTTTAAAATTTT", 0.0, 0.0),
            ("AAAAAGCCCCCTTTTTGGGGG", "CGTGACGTGACGGACT", 3.0, 0.0),
            ("ACGTACGTTGCAGGCATGCAAGCT", "AGCTTGCATGCCTGCAGGTCGACT", 16.0, 14.0),
            ("TAAGCAGAATTCGTAAAACGACGGCCAGT", "TAAGCAGGATCCAGCGGATAACAATTTCACACAGGA", 4.0, 2.0),
        ] {
            assert_eq!(pair_any(a.as_bytes(), b.as_bytes()), any, "{a} {b} pair_any");
            assert_eq!(pair_end(a.as_bytes(), b.as_bytes()), end, "{a} {b} pair_end");
        }
    }

    #[test]
    fn end_stability_matches_primer3() {
        for (seq, stability) in [
            ("GCGCG", 6.86),
            ("TATAT", 0.86),
            ("GTAAAACGACGGCCAGT", 4.00),
            ("AGCGGATAACAATTTCACACAGGA", 3.86),
            ("ATGCGAATTCGCAT", 4.73),
            ("CAGGAAACAGCTATGACCATG", 3.66),
            ("TTTTAAAATTTTAAAATTTT", 1.82),
            ("AAAAAGCCCCCTTTTTGGGGG", 5.40),
            ("CGTGACGTGACGGACT", 3.85),
            ("ACGTACGTTGCAGGCATGCAAGCT", 3.74),
            ("AGCTTGCATGCCTGCAGGTCGACT", 4.18),
        ] {
            let ours = end_stability(seq.as_bytes());
            assert!((ours - stability).abs() < 0.005, "{seq}: {ours} vs Primer3's {stability}");
        }
    }
}
