//! Primer design after Primer3: primer pairs that amplify a region, judged and ranked as Primer3
//! judges and ranks them, and 5' tails for cloning.
//!
//! Each primer anneals with its 3' end toward the region it copies. The designer tries every
//! start and length allowed on each side and holds each candidate to Primer3's rules
//! (Untergasser et al. 2012, Koressaar and Remm 2007), with Primer3's default limits and weights
//! unless [`DesignOptions`] says otherwise:
//!
//! - its melting temperature, by SantaLucia's (1998) nearest-neighbour model with Primer3's salt
//!   correction and conditions, must lie within a window around the target, 57 to 63 °C by
//!   default;
//! - its GC content must lie between 20 and 80%;
//! - no base may repeat more than five times in a row;
//! - it must not pair with itself: Primer3's alignment score of the primer against its own
//!   reverse complement may not exceed 8 anywhere, and may not exceed 3 for an alignment that
//!   ends at its 3' end, which is the primer-dimer test;
//! - it may anneal nowhere else on the template, SlopGene's own check in place of Primer3's
//!   mispriming libraries.
//!
//! Among the candidates that pass, Primer3's penalty ranks them: one point per degree from the
//! target temperature and one per base from the optimum length of 20, the other weights being
//! zero by default. A pair must also pass the same two alignment tests between its primers, and
//! its penalty is the sum of its primers' penalties. When no candidate passes every test, the
//! designer still answers, as Primer3 does when told to pick anyway: it takes the candidates that
//! break the fewest rules, with the rules they break written out as warnings.
//!
//! Temperatures and GC content are those of the annealing part only, as a tail does not pair with
//! the template in the first cycles. The pairing tests use the whole primer, tail included, as
//! Primer3 does with its overhangs.
//!
//! Primer3's thermodynamic model of hairpins and dimers (its `_TH` settings) is not part of this.
//! The pairing tests are Primer3's alignment scores, which it used before that model and still
//! offers.
//!
//! Tails for cloning are made here too: a restriction site with a few spare bases in front of it
//! (so the enzyme can cut near the end of the product), a Golden Gate site that leaves a chosen
//! four-base overhang, and the stretches of a vector that an insert must repeat to be joined by
//! homology (Gibson, NEBuilder HiFi, In-Fusion).

use crate::alphabet::{self, revcomp};
use crate::enzymes::{Enzyme, find_sites};
use crate::error::{Error, Result};
use crate::model::{Document, Primer, Strand};
use crate::primers::{self, BindingRules, TmConditions, find_binding_sites, melting_temperature};
use crate::span::Span;

/// Spare bases put in front of a restriction site so the enzyme can cut near the end of a PCR
/// product. Six bases suit most enzymes. The filler holds no common recognition site and does
/// not pair with itself.
pub const SITE_PADDING: &[u8] = b"TAAGCA";

/// The limits a primer must keep, Primer3's by default.
#[derive(Clone, Debug, PartialEq)]
pub struct Limits {
    /// How far the melting temperature may lie from the target, in degrees. Primer3 allows 57 to
    /// 63 °C around its 60.
    pub tm_window: f64,
    /// GC content, in percent.
    pub min_gc: f64,
    pub max_gc: f64,
    /// The longest run of one base.
    pub max_poly_x: usize,
    /// G or C bases required at the very 3' end.
    pub gc_clamp: usize,
    /// The most G or C among the last five bases.
    pub max_end_gc: usize,
    /// Primer3's alignment score of a primer against its own reverse complement, anywhere.
    pub max_self_any: f64,
    /// The same score for an alignment that ends at the primer's 3' end.
    pub max_self_end: f64,
    /// The same two scores between the two primers of a pair.
    pub max_pair_any: f64,
    pub max_pair_end: f64,
    /// The stability of the last five bases, in kcal/mol. Primer3 does not limit it (100).
    pub max_end_stability: f64,
    /// How far apart the two primers' melting temperatures may be. Primer3 does not limit it
    /// (100).
    pub max_tm_difference: f64,
    /// How many other places on the template a primer may anneal at.
    pub max_other_sites: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Limits {
            tm_window: 3.0,
            min_gc: 20.0,
            max_gc: 80.0,
            max_poly_x: 5,
            gc_clamp: 0,
            max_end_gc: 5,
            max_self_any: 8.0,
            max_self_end: 3.0,
            max_pair_any: 8.0,
            max_pair_end: 3.0,
            max_end_stability: 100.0,
            max_tm_difference: 100.0,
            max_other_sites: 0,
        }
    }
}

/// Primer3's penalty weights: what each fault costs. Only the temperature and the length count by
/// default.
#[derive(Clone, Debug, PartialEq)]
pub struct Weights {
    /// Per degree below and above the target melting temperature.
    pub tm_below: f64,
    pub tm_above: f64,
    /// Per base below and above the optimum length.
    pub size_below: f64,
    pub size_above: f64,
    /// Per percent of GC below and above the optimum.
    pub gc_below: f64,
    pub gc_above: f64,
    /// Per point of the two self-pairing scores and per kcal/mol of 3' end stability.
    pub self_any: f64,
    pub self_end: f64,
    pub end_stability: f64,
    /// Per degree between the two primers of a pair, and per point of their pairing scores.
    pub pair_tm_difference: f64,
    pub pair_any: f64,
    pub pair_end: f64,
}

impl Default for Weights {
    fn default() -> Self {
        Weights {
            tm_below: 1.0,
            tm_above: 1.0,
            size_below: 1.0,
            size_above: 1.0,
            gc_below: 0.0,
            gc_above: 0.0,
            self_any: 0.0,
            self_end: 0.0,
            end_stability: 0.0,
            pair_tm_difference: 0.0,
            pair_any: 0.0,
            pair_end: 0.0,
        }
    }
}

/// How the designer chooses primers.
#[derive(Clone, Debug, PartialEq)]
pub struct DesignOptions {
    /// The melting temperature to aim for (deg C), for the annealing part of each primer.
    pub target_tm: f64,
    /// Shortest, best and longest annealing part, in bases.
    pub min_len: usize,
    pub opt_len: usize,
    pub max_len: usize,
    /// The GC content to aim for, in percent. It counts only with a weight on it.
    pub opt_gc: f64,
    /// How far outside the region a primer may start. 0 amplifies exactly the region; more lets
    /// the designer move each primer outward by up to this many bases to find a better one.
    pub flank: usize,
    /// 5' tails, written 5'->3' as they sit on each primer.
    pub forward_tail: Vec<u8>,
    pub reverse_tail: Vec<u8>,
    pub conditions: TmConditions,
    /// How many pairs to return, best first.
    pub count: usize,
    pub limits: Limits,
    pub weights: Weights,
}

/// The melting temperature a design aims for unless told another, in degrees Celsius.
pub const DEFAULT_TARGET_TM: f64 = 60.0;

impl Default for DesignOptions {
    fn default() -> Self {
        DesignOptions {
            target_tm: DEFAULT_TARGET_TM,
            min_len: 18,
            opt_len: 20,
            max_len: 27,
            opt_gc: 50.0,
            flank: 0,
            forward_tail: Vec::new(),
            reverse_tail: Vec::new(),
            conditions: TmConditions::default(),
            count: 5,
            limits: Limits::default(),
            weights: Weights::default(),
        }
    }
}

/// One designed primer.
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    /// The whole primer, tail first, 5'->3'.
    pub seq: Vec<u8>,
    /// Bases of tail at its 5' end; the rest anneals.
    pub tail_len: usize,
    /// Where the annealing part sits on the template, in top-strand coordinates.
    pub site: Span,
    /// `Forward` reads like the top strand; `Reverse` like the bottom strand.
    pub strand: Strand,
    /// Melting temperature of the annealing part.
    pub tm: f64,
    pub gc_percent: f64,
    /// Primer3's alignment scores of the whole primer against its own reverse complement:
    /// anywhere, and ending at its 3' end.
    pub self_any: f64,
    pub self_end: f64,
    /// Primer3's stability of the last five bases, in kcal/mol.
    pub end_stability: f64,
    /// Other places on the template where the primer would also anneal.
    pub other_sites: usize,
    /// The rules it breaks, in a few words each. Empty for a primer Primer3 would accept.
    pub warnings: Vec<String>,
    /// Primer3's penalty: the sum of its faults, weighted. Lower is better.
    pub penalty: f64,
}

impl Candidate {
    /// The part that anneals to the template.
    pub fn annealing(&self) -> &[u8] {
        &self.seq[self.tail_len..]
    }

    /// The primer as text.
    pub fn text(&self) -> String {
        String::from_utf8_lossy(&self.seq).into_owned()
    }
}

/// A forward and a reverse primer that work together.
#[derive(Clone, Debug, PartialEq)]
pub struct PrimerPair {
    pub forward: Candidate,
    pub reverse: Candidate,
    /// Template bases the product copies, from the forward primer's annealing start to the
    /// reverse primer's annealing end.
    pub amplified: Span,
    /// Length of the product, tails included.
    pub product_len: usize,
    /// Primer3's alignment scores between the two primers: anywhere, and ending at either 3' end.
    pub compl_any: f64,
    pub compl_end: f64,
    /// The rules the pair breaks between its primers. Each primer's own are in it.
    pub warnings: Vec<String>,
    /// Primer3's pair penalty: both primers' penalties, plus the pair's weighted faults.
    pub penalty: f64,
}

impl PrimerPair {
    /// The two primers as they are stored in a document, named after `name` and the 1-based
    /// positions where they start, with what they anneal to in their descriptions.
    pub fn primers(&self, name: &str, seq_len: usize) -> (Primer, Primer) {
        let (start, end) = self.amplified.to_1based(seq_len);
        let describe = |c: &Candidate| {
            let mut text = format!("anneals {} nt, Tm {:.1} °C", c.annealing().len(), c.tm);
            if c.tail_len > 0 {
                text.push_str(&format!(", {}-nt 5' tail", c.tail_len));
            }
            text
        };
        let mut forward = Primer::new(format!("{name}_{start}_F"), self.forward.text());
        forward.description = describe(&self.forward);
        let mut reverse = Primer::new(format!("{name}_{end}_R"), self.reverse.text());
        reverse.description = describe(&self.reverse);
        (forward, reverse)
    }

    /// How many rules the pair and its primers break together.
    fn faults(&self) -> usize {
        self.forward.warnings.len() + self.reverse.warnings.len() + self.warnings.len()
    }
}

/// How many of the best candidates on each side are checked for other binding sites and paired.
const SHORTLIST: usize = 40;

/// What the judge finds about one primer.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Judged {
    pub tm: f64,
    pub gc_percent: f64,
    pub self_any: f64,
    pub self_end: f64,
    pub end_stability: f64,
    pub penalty: f64,
    pub warnings: Vec<String>,
}

/// The longest run of one base in `seq`.
fn longest_run(seq: &[u8]) -> usize {
    let mut longest = 0;
    let mut run = 0;
    for ii in 0..seq.len() {
        run = if ii > 0 && seq[ii].eq_ignore_ascii_case(&seq[ii - 1]) { run + 1 } else { 1 };
        longest = longest.max(run);
    }
    longest
}

/// Judge one primer by Primer3's rules. `annealing` is the part that pairs with the template,
/// which sets the temperature and the GC content. `full` includes the tail, since a tail can
/// pair with the primer or with the other primer too. The penalty is Primer3's, and the warnings
/// name every limit the primer breaks.
pub(crate) fn judge(full: &[u8], annealing: &[u8], options: &DesignOptions) -> Judged {
    let limits = &options.limits;
    let weights = &options.weights;
    let tm = melting_temperature(annealing, &options.conditions);
    let gc_percent = alphabet::gc_fraction(annealing) * 100.0;
    let self_any = primers::self_any(full);
    let self_end = primers::self_end(full);
    let end_stability = primers::end_stability(annealing);
    let len = annealing.len();

    // Primer3's penalty: each fault times its weight, above or below the optimum.
    let mut penalty = if tm < options.target_tm { weights.tm_below * (options.target_tm - tm) } else { weights.tm_above * (tm - options.target_tm) };
    penalty += if len < options.opt_len { weights.size_below * (options.opt_len - len) as f64 } else { weights.size_above * (len - options.opt_len) as f64 };
    penalty += if gc_percent < options.opt_gc { weights.gc_below * (options.opt_gc - gc_percent) } else { weights.gc_above * (gc_percent - options.opt_gc) };
    penalty += weights.self_any * self_any + weights.self_end * self_end + weights.end_stability * end_stability;

    // The limits, each written out when broken.
    let mut warnings = Vec::new();
    let (min_tm, max_tm) = (options.target_tm - limits.tm_window, options.target_tm + limits.tm_window);
    if tm < min_tm {
        warnings.push(format!("Tm {tm:.1} °C is below {min_tm:.0}"));
    } else if tm > max_tm {
        warnings.push(format!("Tm {tm:.1} °C is above {max_tm:.0}"));
    }
    if gc_percent < limits.min_gc {
        warnings.push(format!("GC {gc_percent:.0}% is below {:.0}%", limits.min_gc));
    } else if gc_percent > limits.max_gc {
        warnings.push(format!("GC {gc_percent:.0}% is above {:.0}%", limits.max_gc));
    }
    let run = longest_run(annealing);
    if run > limits.max_poly_x {
        warnings.push(format!("a run of {run} identical bases, more than {}", limits.max_poly_x));
    }
    let last_five = &annealing[len.saturating_sub(5)..];
    let is_gc = |b: &u8| matches!(b.to_ascii_uppercase(), b'G' | b'C');
    if limits.gc_clamp > 0 && !annealing[len.saturating_sub(limits.gc_clamp)..].iter().all(is_gc) {
        warnings.push(format!("fewer than {} G or C at the 3' end", limits.gc_clamp));
    }
    let end_gc = last_five.iter().filter(|b| is_gc(b)).count();
    if end_gc > limits.max_end_gc {
        warnings.push(format!("{end_gc} G or C among the last five bases, more than {}", limits.max_end_gc));
    }
    if self_any > limits.max_self_any {
        warnings.push(format!("pairs with itself (score {self_any:.1}, limit {:.0})", limits.max_self_any));
    }
    if self_end > limits.max_self_end {
        warnings.push(format!("its 3' end pairs with itself (score {self_end:.1}, limit {:.0})", limits.max_self_end));
    }
    if end_stability > limits.max_end_stability {
        warnings.push(format!("3' end stability {end_stability:.1} kcal/mol, more than {:.1}", limits.max_end_stability));
    }
    Judged { tm, gc_percent, self_any, self_end, end_stability, penalty, warnings }
}

/// Every forward-primer candidate for `region`: starts from `flank` bases before the region to
/// its first base, each length allowed, judged but not yet checked for other sites. Reverse
/// primers likewise, from the region's last base outward.
fn side_candidates(doc: &Document, region: Span, forward: bool, options: &DesignOptions) -> Vec<Candidate> {
    let len = doc.len();
    let circular = doc.is_circular();
    let tail = if forward { &options.forward_tail } else { &options.reverse_tail };
    let mut found = Vec::new();
    for shift in 0..=options.flank {
        for size in options.min_len..=options.max_len {
            // The annealing part on the top strand. A forward primer starts `shift` bases before
            // the region; a reverse primer ends `shift` bases after it.
            let site = if forward {
                if circular {
                    Span::new(region.start + len - shift, region.start + len - shift + size).normalize(len)
                } else if shift > region.start {
                    continue;
                } else {
                    Span::new(region.start - shift, region.start - shift + size)
                }
            } else {
                let end = region.end + shift;
                if !circular && end > len {
                    continue;
                }
                if end < size {
                    if !circular {
                        continue;
                    }
                    Span::new(end + len - size, end + len).normalize(len)
                } else {
                    Span::new(end - size, end).normalize(len.max(1))
                }
            };
            // The primer must stay inside the stretch it amplifies.
            if size > region.len() + 2 * shift || (!circular && site.end > len) {
                continue;
            }
            let top = site.extract(&doc.seq);
            let annealing = if forward { alphabet::upper(&top) } else { revcomp(&alphabet::upper(&top)) };
            if annealing.iter().any(|&b| !alphabet::is_unambiguous(b)) {
                continue;
            }
            let mut full = tail.clone();
            full.extend_from_slice(&annealing);
            let full = alphabet::upper(&full);
            let judged = judge(&full, &annealing, options);
            found.push(Candidate {
                seq: full,
                tail_len: tail.len(),
                site,
                strand: if forward { Strand::Forward } else { Strand::Reverse },
                tm: judged.tm,
                gc_percent: judged.gc_percent,
                self_any: judged.self_any,
                self_end: judged.self_end,
                end_stability: judged.end_stability,
                other_sites: 0,
                warnings: judged.warnings,
                // A primer moved outward costs a little, so the tightest product wins ties.
                penalty: judged.penalty + 0.05 * shift as f64,
            });
        }
    }
    found.sort_by(|a, b| (a.warnings.len(), a.penalty).partial_cmp(&(b.warnings.len(), b.penalty)).unwrap_or(std::cmp::Ordering::Equal));
    found
}

/// Count the places besides its own where a candidate's annealing part binds the template.
/// Primers that start in the wrong place give extra products, so more than the limit allows is
/// a fault.
fn check_other_sites(doc: &Document, candidate: &mut Candidate, limits: &Limits) {
    let rules = BindingRules { min_tm: 35.0, ..BindingRules::default() };
    let sites = find_binding_sites(&doc.seq, doc.is_circular(), candidate.annealing(), &rules);
    let others = sites.iter().filter(|b| !(b.strand == candidate.strand && b.span == candidate.site)).count();
    candidate.other_sites = others;
    if others > limits.max_other_sites {
        candidate.warnings.push(format!("also anneals at {others} other site{}", if others == 1 { "" } else { "s" }));
    }
}

/// Design primer pairs that amplify `region` of `doc`, best first.
///
/// The product holds the whole region; with `flank` above 0 it may also hold up to that many
/// bases on either side. Tails are added to the 5' ends as given. Pairs whose primers pass every
/// limit come first, and when none do, the pairs that break the fewest rules are returned with
/// those rules as warnings.
pub fn design_pairs(doc: &Document, region: Span, options: &DesignOptions) -> Result<Vec<PrimerPair>> {
    if !doc.molecule.is_nucleic() {
        return Err(Error::Invalid("primers need a DNA template".into()));
    }
    if options.min_len < 10 || options.min_len > options.max_len {
        return Err(Error::Invalid(format!("primer lengths {}..{} make no sense; use at least 10 bases", options.min_len, options.max_len)));
    }
    if region.len() < options.min_len + 2 {
        return Err(Error::Invalid(format!("the region is {} bp, too short to amplify with {}-base primers", region.len(), options.min_len)));
    }
    let len = doc.len();
    let region = region.normalize(len);
    let limits = &options.limits;
    let weights = &options.weights;

    // The best candidates on each side, checked for other binding sites.
    let shortlist = |forward: bool| {
        let mut side = side_candidates(doc, region, forward, options);
        side.truncate(SHORTLIST);
        for candidate in &mut side {
            check_other_sites(doc, candidate, limits);
        }
        side.sort_by(|a, b| (a.warnings.len(), a.penalty).partial_cmp(&(b.warnings.len(), b.penalty)).unwrap_or(std::cmp::Ordering::Equal));
        side
    };
    let forwards = shortlist(true);
    let reverses = shortlist(false);
    if forwards.is_empty() || reverses.is_empty() {
        return Err(Error::Invalid("no primer fits there: the ends of the region hold ambiguous bases or run off the sequence".into()));
    }

    // Every forward with every reverse: their own penalties and faults, plus how the two pair
    // with each other and how far apart their melting temperatures are.
    let mut pairs = Vec::new();
    for forward in &forwards {
        for reverse in &reverses {
            let start = forward.site.start;
            let mut end = reverse.site.end;
            if end <= start {
                end += len;
            }
            let amplified = Span::new(start, end);
            if amplified.len() > len || amplified.len() < forward.site.len().max(reverse.site.len()) {
                continue;
            }
            let compl_any = primers::pair_any(&forward.seq, &reverse.seq);
            let compl_end = primers::pair_end(&forward.seq, &reverse.seq);
            let tm_gap = (forward.tm - reverse.tm).abs();
            let penalty = forward.penalty + reverse.penalty + weights.pair_tm_difference * tm_gap + weights.pair_any * compl_any + weights.pair_end * compl_end;
            let mut warnings = Vec::new();
            if compl_any > limits.max_pair_any {
                warnings.push(format!("the primers pair with each other (score {compl_any:.1}, limit {:.0})", limits.max_pair_any));
            }
            if compl_end > limits.max_pair_end {
                warnings.push(format!("the primers' 3' ends pair with each other (score {compl_end:.1}, limit {:.0})", limits.max_pair_end));
            }
            if tm_gap > limits.max_tm_difference {
                warnings.push(format!("melting temperatures {tm_gap:.1} °C apart, more than {:.0}", limits.max_tm_difference));
            }
            pairs.push(PrimerPair {
                forward: forward.clone(),
                reverse: reverse.clone(),
                amplified,
                product_len: amplified.len() + forward.tail_len + reverse.tail_len,
                compl_any,
                compl_end,
                warnings,
                penalty,
            });
        }
    }
    pairs.sort_by(|a, b| (a.faults(), a.penalty).partial_cmp(&(b.faults(), b.penalty)).unwrap_or(std::cmp::Ordering::Equal));

    // Keep pairs that differ: no two with the same forward or the same reverse primer, so the
    // list offers real alternatives rather than one primer with many partners.
    let mut chosen: Vec<PrimerPair> = Vec::new();
    for pair in pairs {
        if chosen.len() >= options.count.max(1) {
            break;
        }
        if chosen.iter().any(|c| c.forward.seq == pair.forward.seq || c.reverse.seq == pair.reverse.seq) {
            continue;
        }
        chosen.push(pair);
    }
    if chosen.is_empty() {
        return Err(Error::Invalid("the primers found do not face each other across the region".into()));
    }
    Ok(chosen)
}

/// Design the best primer pair that amplifies exactly `span` of `doc`, with 5' tails
/// (restriction sites, overlaps) when given, as primers named after the document. The pair is
/// the first that [`design_pairs`] returns at `target_tm`.
pub fn design_pcr_primers(
    doc: &Document,
    span: Span,
    target_tm: f64,
    forward_tail: &[u8],
    reverse_tail: &[u8],
    conditions: &TmConditions,
) -> Result<(Primer, Primer)> {
    let options = DesignOptions {
        target_tm,
        forward_tail: forward_tail.to_vec(),
        reverse_tail: reverse_tail.to_vec(),
        conditions: *conditions,
        count: 1,
        ..DesignOptions::default()
    };
    let pairs = design_pairs(doc, span, &options)?;
    Ok(pairs[0].primers(&doc.name, doc.len()))
}

/// A concrete base for each position of an enzyme site, so it can be written into a primer:
/// an ambiguity code becomes the first base it stands for among G, C, A and T.
fn concrete(site: &[u8]) -> Vec<u8> {
    site.iter()
        .map(|&code| {
            if alphabet::is_unambiguous(code) {
                code.to_ascii_uppercase()
            } else {
                *b"GCAT".iter().find(|&&base| alphabet::iupac_matches(code, base)).unwrap_or(&b'A')
            }
        })
        .collect()
}

/// A 5' tail that puts a restriction site on a primer, with `padding` spare bases in front.
///
/// The site is written so the enzyme cuts toward the primer's 3' end, which is where the insert
/// lies for forward and reverse primers alike. That matters only for enzymes that cut outside a
/// site that is not palindromic (Type IIS); for those the bases between the site and the cut
/// are filled with A.
pub fn site_tail(enzyme: &Enzyme, padding: usize) -> Vec<u8> {
    let site_len = enzyme.site.len() as i32;
    let upstream = enzyme.cuts.iter().all(|&(top, bottom)| top < 0 && bottom < 0);
    let mut tail: Vec<u8> = SITE_PADDING.iter().copied().cycle().take(padding).collect();
    if upstream {
        // A site that cuts before itself is turned round, so the cut falls after it.
        tail.extend(revcomp(&concrete(&enzyme.site)));
    } else {
        tail.extend(concrete(&enzyme.site));
    }
    // Type IIS: room for the cut beyond the site.
    let reach = enzyme.cuts.iter().map(|&(top, bottom)| if upstream { -top.min(bottom) } else { top.max(bottom) - site_len }).max().unwrap_or(0);
    tail.extend(std::iter::repeat_n(b'A', reach.max(0) as usize));
    tail
}

/// A 5' tail for Golden Gate assembly: spare bases, the Type IIS site, the spacer up to the cut,
/// then `overhang`, the bases the cut leaves single-stranded at the end of the product.
pub fn golden_gate_tail(enzyme: &Enzyme, overhang: &[u8], padding: usize) -> Result<Vec<u8>> {
    let site_len = enzyme.site.len() as i32;
    let Some(&(top, bottom)) = enzyme.cuts.first() else {
        return Err(Error::Invalid(format!("{} has no cut positions", enzyme.name)));
    };
    if enzyme.cuts.len() != 1 || top <= site_len || bottom <= top {
        return Err(Error::Invalid(format!("{} does not cut after its site leaving a 5' overhang, so it cannot make Golden Gate ends", enzyme.name)));
    }
    let width = (bottom - top) as usize;
    let overhang = alphabet::upper(overhang);
    if overhang.len() != width || overhang.iter().any(|&b| !alphabet::is_unambiguous(b)) {
        return Err(Error::Invalid(format!("{} leaves {width}-base overhangs; give {width} plain bases", enzyme.name)));
    }
    let mut tail: Vec<u8> = SITE_PADDING.iter().copied().cycle().take(padding).collect();
    tail.extend(concrete(&enzyme.site));
    // The spacer between the site and the top-strand cut; A is as good as any base there.
    tail.extend(std::iter::repeat_n(b'A', (top - site_len) as usize));
    tail.extend_from_slice(&overhang);
    Ok(tail)
}

/// Tails that let a PCR product be joined into a vector by homology (Gibson, NEBuilder HiFi,
/// In-Fusion). The insert goes in place of `opening` (an empty span is a point between two
/// bases). The forward tail repeats the `overlap` vector bases before the opening; the reverse
/// tail repeats, on the other strand, the `overlap` bases after it.
pub fn homology_tails(vector: &Document, opening: Span, overlap: usize) -> Result<(Vec<u8>, Vec<u8>)> {
    let len = vector.len();
    if overlap == 0 || overlap * 2 > len {
        return Err(Error::Invalid(format!("an overlap of {overlap} bases does not fit a {len} bp vector")));
    }
    let circular = vector.is_circular();
    if !circular && (opening.start < overlap || opening.end + overlap > len) {
        return Err(Error::Invalid(format!("the vector has fewer than {overlap} bases on one side of the opening")));
    }
    // Bases just before the opening, and just after it, read along the top strand.
    let before = Span::new(opening.start + len - overlap, opening.start + len).normalize(len);
    let after = Span::new(opening.end, opening.end + overlap).normalize(len);
    let forward = alphabet::upper(&before.extract(&vector.seq));
    let reverse = revcomp(&alphabet::upper(&after.extract(&vector.seq)));
    Ok((forward, reverse))
}

/// Whether a primer's tail holds a site of any of `enzymes`, so a person can see that a tail
/// will be cut where it should be.
pub fn tail_sites(candidate: &Candidate, enzymes: &[&Enzyme]) -> Vec<String> {
    let mut names: Vec<String> = find_sites(&candidate.seq[..candidate.tail_len], false, enzymes).into_iter().map(|s| s.enzyme.to_string()).collect();
    names.sort();
    names.dedup();
    names
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzymes::EnzymeDb;
    use crate::model::Topology;
    use crate::pcr::pcr;

    fn random_dna(n: usize, seed: u64) -> Vec<u8> {
        // A small linear congruential generator keeps the tests deterministic.
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn the_best_pair_amplifies_exactly_the_region_and_passes_every_rule() {
        let doc = Document::new_dna("t", random_dna(1500, 3), Topology::Circular);
        let region = Span::new(400, 1100);
        let pairs = design_pairs(&doc, region, &DesignOptions::default()).unwrap();
        assert!(!pairs.is_empty());
        let best = &pairs[0];
        assert_eq!(best.amplified, region, "with no flank the product is the region");
        // Random DNA offers primers Primer3 would accept: within 57 to 63 °C and clean.
        assert!(best.forward.warnings.is_empty() && best.reverse.warnings.is_empty() && best.warnings.is_empty(), "{best:?}");
        assert!((57.0..=63.0).contains(&best.forward.tm) && (57.0..=63.0).contains(&best.reverse.tm), "{} {}", best.forward.tm, best.reverse.tm);
        assert!(best.forward.self_any <= 8.0 && best.forward.self_end <= 3.0 && best.compl_end <= 3.0);
        // The pair amplifies the region in a simulated PCR.
        let (f, r) = best.primers("t", doc.len());
        let product = pcr(&doc, &f, &r, &BindingRules::default()).unwrap();
        assert_eq!(product.doc.seq, region.extract(&doc.seq));
        // Pairs are listed best first and differ in both primers.
        for window in pairs.windows(2) {
            assert!((window[0].faults(), window[0].penalty) <= (window[1].faults(), window[1].penalty));
            assert_ne!(window[0].forward.seq, window[1].forward.seq);
        }
    }

    #[test]
    fn the_penalty_is_primer3s() {
        // Primer3 (primer3-py 2.4.0) gives these penalties at its defaults: one point per degree
        // from 60 °C and one per base from 20.
        let options = DesignOptions::default();
        let m13_forward = judge(b"GTAAAACGACGGCCAGT", b"GTAAAACGACGGCCAGT", &options);
        assert!((m13_forward.penalty - 8.3045).abs() < 0.001, "{}", m13_forward.penalty);
        let m13_reverse = judge(b"AGCGGATAACAATTTCACACAGGA", b"AGCGGATAACAATTTCACACAGGA", &options);
        assert!((m13_reverse.penalty - 4.8059).abs() < 0.001, "{}", m13_reverse.penalty);
        // The warnings name the limits broken: M13 forward melts below 57 °C, and its 3' end
        // pairs with another copy of itself over Primer3's limit of 3.
        assert_eq!(m13_forward.warnings, vec!["Tm 54.7 °C is below 57".to_string(), "its 3' end pairs with itself (score 4.0, limit 3)".to_string()]);
        assert!(m13_reverse.warnings.is_empty(), "{:?}", m13_reverse.warnings);
        // A palindrome pairs with itself far over both limits.
        let palindrome = judge(b"ATGCGAATTCGCAT", b"ATGCGAATTCGCAT", &options);
        assert!(palindrome.warnings.iter().any(|w| w.starts_with("pairs with itself (score 14.0")), "{:?}", palindrome.warnings);
        assert!(palindrome.warnings.iter().any(|w| w.starts_with("its 3' end pairs with itself (score 14.0")), "{:?}", palindrome.warnings);
        // A run of six breaks the poly-X limit, and the GC limits are Primer3's.
        let run = judge(b"GACTGACTGAAAAAATCAGTCG", b"GACTGACTGAAAAAATCAGTCG", &options);
        assert!(run.warnings.iter().any(|w| w == "a run of 6 identical bases, more than 5"), "{:?}", run.warnings);
        let at = judge(b"TTTTAAAATTTTAAAATTTT", b"TTTTAAAATTTTAAAATTTT", &options);
        assert!(at.warnings.iter().any(|w| w == "GC 0% is below 20%"), "{:?}", at.warnings);
    }

    #[test]
    fn when_no_primer_passes_the_least_faulty_is_picked_anyway() {
        // A region whose ends keep only every fourth G or C: no length up to 27 reaches 57 °C.
        let mut bases = random_dna(600, 11);
        for (ii, b) in bases.iter_mut().enumerate() {
            if ((200..230).contains(&ii) || (370..400).contains(&ii)) && ii % 4 != 0 {
                *b = match *b {
                    b'G' => b'A',
                    b'C' => b'T',
                    other => other,
                };
            }
        }
        let doc = Document::new_dna("t", bases, Topology::Linear);
        let pairs = design_pairs(&doc, Span::new(200, 400), &DesignOptions::default()).unwrap();
        let best = &pairs[0];
        assert!(best.forward.warnings.iter().any(|w| w.contains("is below 57")), "{:?}", best.forward.warnings);
        assert_eq!(best.forward.warnings.len(), 1, "the temperature is its only fault: {:?}", best.forward.warnings);
        assert_eq!(best.amplified, Span::new(200, 400));
        // The pairs still come fewest faults first.
        for window in pairs.windows(2) {
            assert!(window[0].faults() <= window[1].faults());
        }
    }

    #[test]
    fn a_flank_lets_the_primers_move_outward_but_never_inward() {
        let doc = Document::new_dna("t", random_dna(1200, 9), Topology::Linear);
        let region = Span::new(300, 800);
        let options = DesignOptions { flank: 40, count: 8, ..DesignOptions::default() };
        for pair in design_pairs(&doc, region, &options).unwrap() {
            assert!(pair.amplified.start <= region.start && pair.amplified.start + 40 >= region.start);
            assert!(pair.amplified.end >= region.end && pair.amplified.end <= region.end + 40);
        }
    }

    #[test]
    fn tails_end_up_on_the_product_and_sites_cut_toward_the_insert() {
        let db = EnzymeDb::builtin();
        let ecori = db.get("EcoRI").unwrap();
        let bamhi = db.get("BamHI").unwrap();
        let doc = Document::new_dna("t", random_dna(900, 5), Topology::Linear);
        let options = DesignOptions { forward_tail: site_tail(ecori, 6), reverse_tail: site_tail(bamhi, 6), ..DesignOptions::default() };
        let best = &design_pairs(&doc, Span::new(200, 700), &options).unwrap()[0];
        assert_eq!(best.forward.tail_len, 12);
        assert_eq!(tail_sites(&best.forward, &[ecori, bamhi]), vec!["EcoRI".to_string()]);
        let (f, r) = best.primers("t", doc.len());
        let product = pcr(&doc, &f, &r, &BindingRules::default()).unwrap();
        assert_eq!(product.doc.len(), 500 + 24);
        assert!(product.doc.seq.starts_with(b"TAAGCAGAATTC"));
        assert!(product.doc.seq.ends_with(&revcomp(b"TAAGCAGGATCC")));
        // The tails count in the pairing tests but not in the temperature.
        let bare = &design_pairs(&doc, Span::new(200, 700), &DesignOptions::default()).unwrap()[0];
        assert!((best.forward.tm - bare.forward.tm).abs() < 1e-9 || best.forward.annealing() != bare.forward.annealing());
        // The named primers carry the tails and say what anneals.
        let (f, _) = design_pcr_primers(&doc, Span::new(200, 700), 60.0, &site_tail(ecori, 6), &[], &TmConditions::default()).unwrap();
        assert!(f.seq.starts_with("TAAGCAGAATTC"));
        assert!(f.description.contains("12-nt 5' tail"), "{}", f.description);
    }

    #[test]
    fn a_golden_gate_tail_leaves_the_chosen_overhang() {
        let db = EnzymeDb::builtin();
        let bsai = db.get("BsaI").unwrap();
        let tail = golden_gate_tail(bsai, b"AATG", 4).unwrap();
        assert_eq!(tail, b"TAAGGGTCTCAAATG".to_vec());
        // BsaI cuts the top strand one base after its site, so the 4 bases that follow the
        // spacer are the single-stranded end.
        let doc = Document::new_dna("p", [tail.clone(), random_dna(60, 2)].concat(), Topology::Linear);
        let cut = crate::digest::digest(&doc, &[bsai]).unwrap();
        let piece = cut.fragments.iter().max_by_key(|f| f.len()).unwrap();
        assert_eq!(crate::digest::overhang_bases(&piece.doc, true), b"AATG".to_vec());
        assert!(golden_gate_tail(bsai, b"AAT", 4).is_err(), "BsaI overhangs are four bases");
        assert!(golden_gate_tail(db.get("EcoRI").unwrap(), b"AATT", 4).is_err(), "EcoRI cuts inside its site");
    }

    #[test]
    fn homology_tails_repeat_the_vector_on_both_sides_of_the_opening() {
        let vector = Document::new_dna("v", random_dna(300, 4), Topology::Circular);
        let (f, r) = homology_tails(&vector, Span::new(100, 120), 20).unwrap();
        assert_eq!(f, vector.seq[80..100].to_vec());
        assert_eq!(r, revcomp(&vector.seq[120..140]));
        // An opening at the origin takes the bases before it from the end of the circle.
        let (f, _) = homology_tails(&vector, Span::new(0, 0), 15).unwrap();
        assert_eq!(f, vector.seq[285..300].to_vec());
    }
}
