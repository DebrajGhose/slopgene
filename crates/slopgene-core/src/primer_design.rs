//! Primer design: primer pairs that amplify a region, ranked, and 5' tails for cloning.
//!
//! Each primer anneals with its 3' end toward the region it copies. The designer tries every
//! start and length allowed on each side, scores each candidate, and pairs the best of them.
//! A candidate's penalty adds up, in degrees or in comparable steps:
//!
//! - the distance of its melting temperature from the target,
//! - GC content outside 40 to 60%,
//! - a 3' end without a G or C clamp, or with more than three G/C among its last five bases,
//! - a 3' end that pairs with another copy of the primer, a hairpin, a run of one base,
//! - each other place on the template where it would also anneal.
//!
//! A pair adds the difference between its two melting temperatures and how well the two
//! primers' 3' ends pair with each other. Temperatures are those of the annealing part only, as
//! a tail does not pair with the template in the first cycles. The weights follow the usual
//! rules of thumb (Primer3's defaults and the guidance of polymerase suppliers): close Tm first,
//! then a clean 3' end.
//!
//! Tails for cloning are made here too: a restriction site with a few spare bases in front of it
//! (so the enzyme can cut near the end of the product), a Golden Gate site that leaves a chosen
//! four-base overhang, and the stretches of a vector that an insert must repeat to be joined by
//! homology (Gibson, NEBuilder HiFi, In-Fusion).

use crate::alphabet::{self, revcomp};
use crate::enzymes::{Enzyme, find_sites};
use crate::error::{Error, Result};
use crate::model::{Document, Primer, Strand};
use crate::primers::{BindingRules, TmConditions, analyze, find_binding_sites, melting_temperature};
use crate::span::Span;

/// Spare bases put in front of a restriction site so the enzyme can cut near the end of a PCR
/// product. Six bases suit most enzymes. The filler holds no common recognition site and does
/// not pair with itself.
pub const SITE_PADDING: &[u8] = b"TAAGCA";

/// How the designer chooses primers.
#[derive(Clone, Debug, PartialEq)]
pub struct DesignOptions {
    /// The melting temperature to aim for (deg C), for the annealing part of each primer.
    pub target_tm: f64,
    /// Shortest and longest annealing part, in bases.
    pub min_len: usize,
    pub max_len: usize,
    /// How far outside the region a primer may start. 0 amplifies exactly the region; more lets
    /// the designer move each primer outward by up to this many bases to find a better one.
    pub flank: usize,
    /// 5' tails, written 5'->3' as they sit on each primer.
    pub forward_tail: Vec<u8>,
    pub reverse_tail: Vec<u8>,
    pub conditions: TmConditions,
    /// How many pairs to return, best first.
    pub count: usize,
}

/// The melting temperature a design aims for unless told another, in degrees Celsius.
pub const DEFAULT_TARGET_TM: f64 = 60.0;

impl Default for DesignOptions {
    fn default() -> Self {
        DesignOptions {
            target_tm: DEFAULT_TARGET_TM,
            min_len: 18,
            max_len: 30,
            flank: 0,
            forward_tail: Vec::new(),
            reverse_tail: Vec::new(),
            conditions: TmConditions::default(),
            count: 5,
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
    /// Other places on the template where the primer would also anneal.
    pub other_sites: usize,
    /// What a person should know about it, in a few words each.
    pub warnings: Vec<String>,
    /// The sum of its faults; lower is better.
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
    /// Longest stretch where one primer's 3' end pairs with the other primer.
    pub cross_dimer: usize,
    pub warnings: Vec<String>,
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
}

/// How many of the best candidates on each side are checked for other binding sites and paired.
const SHORTLIST: usize = 40;

/// The fault score of one primer, with its warnings. `annealing` is the part that pairs with the
/// template; `full` includes the tail, since a tail can fold back or pair with itself too.
fn judge(full: &[u8], annealing: &[u8], options: &DesignOptions) -> (f64, f64, f64, Vec<String>) {
    let report = analyze(full, &options.conditions);
    let tm = melting_temperature(annealing, &options.conditions);
    let gc_percent = alphabet::gc_fraction(annealing) * 100.0;
    let mut penalty = (tm - options.target_tm).abs();
    let mut warnings = Vec::new();

    // GC content outside the comfortable range.
    if gc_percent < 40.0 {
        penalty += (40.0 - gc_percent) * 0.2;
    } else if gc_percent > 60.0 {
        penalty += (gc_percent - 60.0) * 0.2;
    }
    if !(30.0..=70.0).contains(&gc_percent) {
        warnings.push(format!("GC {gc_percent:.0}%"));
    }

    // The 3' end: one G or C as the last base helps the polymerase start; more than three among
    // the last five makes the end sticky enough to start in the wrong places.
    let last = annealing.last().copied().unwrap_or(b'A').to_ascii_uppercase();
    if !matches!(last, b'G' | b'C') {
        penalty += 1.0;
    }
    let tail_gc = annealing[annealing.len().saturating_sub(5)..].iter().filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C')).count();
    if tail_gc > 3 {
        penalty += 2.0 * (tail_gc - 3) as f64;
        warnings.push(format!("{tail_gc} G/C in the last 5 bases"));
    }

    // Pairing with itself: a 3' end that pairs can be extended into a primer dimer, and a hairpin
    // hides the primer from its template.
    if report.three_prime_dimer >= 4 {
        penalty += 2.0 * (report.three_prime_dimer - 3) as f64;
        warnings.push(format!("3' end pairs with itself over {} bases", report.three_prime_dimer));
    }
    if report.hairpin_stem >= 4 {
        penalty += 2.0 * (report.hairpin_stem - 3) as f64;
        warnings.push(format!("hairpin with a {}-base stem", report.hairpin_stem));
    }
    if report.longest_run >= 5 {
        penalty += 2.0 * (report.longest_run - 4) as f64;
        warnings.push(format!("run of {} identical bases", report.longest_run));
    }
    (penalty, tm, gc_percent, warnings)
}

/// The longest stretch ending at `a`'s 3' base that pairs with any part of `b`. Such a stretch
/// lets the polymerase extend `a` on `b`, which makes primer dimers.
pub fn three_prime_pairing(a: &[u8], b: &[u8]) -> usize {
    let a = alphabet::upper(a);
    let b_rc = revcomp(&alphabet::upper(b));
    let mut best = 0;
    // ii = where in b's reverse complement the run that ends at a's last base ends.
    for ii in 0..b_rc.len() {
        let mut run = 0;
        while run < a.len() && run <= ii && a[a.len() - 1 - run] == b_rc[ii - run] {
            run += 1;
        }
        best = best.max(run);
    }
    best
}

/// The longest stretch at either primer's 3' end that pairs with the other primer.
pub fn cross_dimer(a: &[u8], b: &[u8]) -> usize {
    three_prime_pairing(a, b).max(three_prime_pairing(b, a))
}

/// Every forward-primer candidate for `region`: starts from `flank` bases before the region to
/// its first base, each length allowed, scored but not yet checked for other sites.
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
            let (penalty, tm, gc_percent, warnings) = judge(&full, &annealing, options);
            found.push(Candidate {
                seq: full,
                tail_len: tail.len(),
                site,
                strand: if forward { Strand::Forward } else { Strand::Reverse },
                tm,
                gc_percent,
                other_sites: 0,
                warnings,
                penalty: penalty + 0.05 * shift as f64,
            });
        }
    }
    found.sort_by(|a, b| a.penalty.total_cmp(&b.penalty));
    found
}

/// Count the places besides its own where a candidate's annealing part binds the template, and
/// add them to its penalty. Primers that start in the wrong place give extra products.
fn check_other_sites(doc: &Document, candidate: &mut Candidate) {
    let rules = BindingRules { min_tm: 35.0, ..BindingRules::default() };
    let sites = find_binding_sites(&doc.seq, doc.is_circular(), candidate.annealing(), &rules);
    let others = sites.iter().filter(|b| !(b.strand == candidate.strand && b.span == candidate.site)).count();
    candidate.other_sites = others;
    if others > 0 {
        candidate.penalty += 5.0 * others as f64;
        candidate.warnings.push(format!("also anneals at {others} other site{}", if others == 1 { "" } else { "s" }));
    }
}

/// Design primer pairs that amplify `region` of `doc`, best first.
///
/// The product holds the whole region; with `flank` above 0 it may also hold up to that many
/// bases on either side. Tails are added to the 5' ends as given.
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

    // The best candidates on each side, checked for other binding sites.
    let shortlist = |forward: bool| {
        let mut side = side_candidates(doc, region, forward, options);
        side.truncate(SHORTLIST);
        for candidate in &mut side {
            check_other_sites(doc, candidate);
        }
        side.sort_by(|a, b| a.penalty.total_cmp(&b.penalty));
        side
    };
    let forwards = shortlist(true);
    let reverses = shortlist(false);
    if forwards.is_empty() || reverses.is_empty() {
        return Err(Error::Invalid("no primer fits there: the ends of the region hold ambiguous bases or run off the sequence".into()));
    }

    // Every forward with every reverse: their own penalties, plus how far apart their melting
    // temperatures are and whether their 3' ends pair with each other.
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
            let cross = cross_dimer(&forward.seq, &reverse.seq);
            let tm_gap = (forward.tm - reverse.tm).abs();
            let mut penalty = forward.penalty + reverse.penalty + tm_gap;
            let mut warnings = Vec::new();
            if cross >= 4 {
                penalty += 2.0 * (cross - 3) as f64;
                warnings.push(format!("the primers' 3' ends pair with each other over {cross} bases"));
            }
            if tm_gap > 5.0 {
                warnings.push(format!("melting temperatures {tm_gap:.1} C apart"));
            }
            pairs.push(PrimerPair {
                forward: forward.clone(),
                reverse: reverse.clone(),
                amplified,
                product_len: amplified.len() + forward.tail_len + reverse.tail_len,
                cross_dimer: cross,
                warnings,
                penalty,
            });
        }
    }
    pairs.sort_by(|a, b| a.penalty.total_cmp(&b.penalty));

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
    fn the_best_pair_amplifies_exactly_the_region_with_matched_temperatures() {
        let doc = Document::new_dna("t", random_dna(1500, 3), Topology::Circular);
        let region = Span::new(400, 1100);
        let pairs = design_pairs(&doc, region, &DesignOptions::default()).unwrap();
        assert!(!pairs.is_empty());
        let best = &pairs[0];
        assert_eq!(best.amplified, region, "with no flank the product is the region");
        assert!((best.forward.tm - 60.0).abs() < 4.0 && (best.reverse.tm - 60.0).abs() < 4.0, "{} {}", best.forward.tm, best.reverse.tm);
        // The pair amplifies the region in a simulated PCR.
        let (f, r) = best.primers("t", doc.len());
        let product = pcr(&doc, &f, &r, &BindingRules::default()).unwrap();
        assert_eq!(product.doc.seq, region.extract(&doc.seq));
        // Pairs are listed best first and differ in both primers.
        for window in pairs.windows(2) {
            assert!(window[0].penalty <= window[1].penalty);
            assert_ne!(window[0].forward.seq, window[1].forward.seq);
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

    #[test]
    fn three_prime_pairing_finds_ends_that_extend_on_each_other() {
        // GAATTC pairs with itself, so a primer ending in it pairs with a copy of itself over
        // those six bases.
        assert_eq!(three_prime_pairing(b"AAAAGAATTC", b"AAAAGAATTC"), 6);
        // The last eight bases of a, AAAAGCGC, pair with GCGCTTTT inside b.
        assert_eq!(three_prime_pairing(b"AAAAAAAAGCGC", b"TTTTGCGCTTTT"), 8);
        assert_eq!(three_prime_pairing(b"AAAAAAAAAAAA", b"CCCCCCCCCCCC"), 0);
        // A pairing stretch away from the 3' ends cannot be extended, so it does not count.
        assert!(cross_dimer(b"GAATTCAAAAAA", b"GAATTCAAAAAA") < 4);
    }
}
