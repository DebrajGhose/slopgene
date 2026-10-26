//! Primers: melting temperature, binding sites, quality checks and simple design.
//!
//! Melting temperatures use the nearest-neighbour model with the unified parameters of
//! SantaLucia (1998) and SantaLucia & Hicks (2004), with the entropy salt correction from the
//! same papers. The defaults match the conditions SnapGene states it assumes: 50 mM Na+ and
//! 0.25 uM oligo, no Mg2+.

use crate::alphabet::{self, iupac_matches};
use crate::model::{BindingSite, Document, Primer, Strand, XmlExtra};
use crate::span::Span;

/// Gas constant in cal/(K mol).
const R: f64 = 1.9872;

/// Conditions for a melting-temperature estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TmConditions {
    /// Monovalent cations (Na+, K+), millimolar.
    pub monovalent_mm: f64,
    /// Mg2+, millimolar. Converted to a sodium equivalent (von Ahsen et al. 2001).
    pub magnesium_mm: f64,
    /// dNTPs, millimolar (they bind Mg2+).
    pub dntp_mm: f64,
    /// Oligo concentration, nanomolar. For a primer in excess over its template this is the
    /// primer concentration.
    pub oligo_nm: f64,
}

impl Default for TmConditions {
    fn default() -> Self {
        TmConditions { monovalent_mm: 50.0, magnesium_mm: 0.0, dntp_mm: 0.0, oligo_nm: 250.0 }
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

/// Melting temperature (deg C) of `seq` annealed to its perfect complement.
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
    // Sodium-equivalent concentration (M): Na+ plus 120 * sqrt(free Mg2+), in mM.
    let free_mg = (conditions.magnesium_mm - conditions.dntp_mm).max(0.0);
    let na_eq = (conditions.monovalent_mm + 120.0 * free_mg.sqrt()) / 1000.0;
    let ds_salt = ds + 0.368 * (seq.len() as f64 - 1.0) * na_eq.max(1e-6).ln();
    // A primer anneals to a template it is in excess over, so the concentration term is the
    // oligo concentration itself (x = 1 in SantaLucia's notation). SnapGene's "0.5 uM each for
    // two annealed oligos" gives the same value, 1 uM total divided by 4.
    let oligo = conditions.oligo_nm.max(1e-3) * 1e-9;
    dh * 1000.0 / (ds_salt + R * oligo.ln()) - 273.15
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

/// Analyse a primer.
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

    let mut warnings = Vec::new();
    if n < 15 {
        warnings.push(format!("short primer ({n} nt)"));
    }
    if !(40.0..=60.0).contains(&gc_percent) {
        warnings.push(format!("GC content {gc_percent:.0}% is outside 40-60%"));
    }
    if three_prime_dimer >= 4 {
        warnings.push(format!("3' end can pair with itself over {three_prime_dimer} bases (primer dimer)"));
    }
    if hairpin_stem >= 4 {
        warnings.push(format!("possible hairpin with a {hairpin_stem}-base stem"));
    }
    if longest_run >= 5 {
        warnings.push(format!("run of {longest_run} identical bases"));
    }
    if three_prime_gc > 3 {
        warnings.push(format!("{three_prime_gc} G/C among the last 5 bases (3' end too sticky)"));
    }
    PrimerReport { length: n, gc_percent, tm, self_dimer, three_prime_dimer, hairpin_stem, longest_run, three_prime_gc, warnings }
}

/// Pick the length (within `min_len..=max_len`) of a primer starting at the 5' end of `region`
/// so its Tm is closest to `target_tm`. `region` is read 5'->3' along the primer's strand.
pub fn choose_primer_length(region: &[u8], target_tm: f64, min_len: usize, max_len: usize, conditions: &TmConditions) -> Option<usize> {
    let max_len = max_len.min(region.len());
    if min_len > max_len {
        return None;
    }
    (min_len..=max_len).min_by(|&a, &b| {
        let da = (melting_temperature(&region[..a], conditions) - target_tm).abs();
        let db = (melting_temperature(&region[..b], conditions) - target_tm).abs();
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    })
}

/// Design a primer pair that amplifies exactly `span` of `doc`, optionally with 5' tails
/// (restriction sites, overlaps). Returns (forward, reverse) primers.
pub fn design_pcr_primers(
    doc: &Document,
    span: Span,
    target_tm: f64,
    forward_tail: &[u8],
    reverse_tail: &[u8],
    conditions: &TmConditions,
) -> crate::error::Result<(Primer, Primer)> {
    let region = doc.subseq(span);
    if region.len() < 30 {
        return Err(crate::error::Error::Invalid("region is too short to amplify (needs at least 30 bp)".into()));
    }
    let fwd_len = choose_primer_length(&region, target_tm, 18, 32, conditions).unwrap_or(20);
    let rc_region = alphabet::revcomp(&region);
    let rev_len = choose_primer_length(&rc_region, target_tm, 18, 32, conditions).unwrap_or(20);
    let mut fwd = forward_tail.to_vec();
    fwd.extend_from_slice(&region[..fwd_len]);
    let mut rev = reverse_tail.to_vec();
    rev.extend_from_slice(&rc_region[..rev_len]);
    let (a, b) = span.to_1based(doc.len());
    let mut forward = Primer::new(format!("{}_{a}_F", doc.name), String::from_utf8_lossy(&fwd).to_uppercase());
    let mut reverse = Primer::new(format!("{}_{b}_R", doc.name), String::from_utf8_lossy(&rev).to_uppercase());
    forward.description = format!("anneals {fwd_len} nt, Tm {:.1} °C", melting_temperature(&region[..fwd_len], conditions));
    reverse.description = format!("anneals {rev_len} nt, Tm {:.1} °C", melting_temperature(&rc_region[..rev_len], conditions));
    Ok((forward, reverse))
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
        assert!(good.warnings.iter().all(|w| !w.contains("dimer")), "{:?}", good.warnings);
        let dimer = analyze(b"ATGCGAATTCGCAT", &c);
        assert!(dimer.three_prime_dimer >= 4);
    }
}
