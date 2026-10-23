//! Restriction enzymes: the built-in database (commercially available enzymes from REBASE),
//! named enzyme sets, and recognition-site search.
//!
//! Cut positions are stored as *boundaries* relative to the first base of the recognition
//! site: boundary 0 lies just before the site, boundary `n` just after its `n`th base, and
//! negative boundaries lie upstream. Each double-strand cut has a top-strand boundary and a
//! bottom-strand boundary, both in top-strand coordinates. EcoRI (G^AATTC) cuts the top strand
//! at 1 and the bottom strand at 5, leaving a 4-base 5' overhang.

use crate::alphabet::{self, iupac_matches};
use crate::error::{Error, Result};
use crate::model::Strand;
use std::collections::HashMap;
use std::sync::OnceLock;

/// The embedded table, generated from REBASE (see the header of the file).
const BUILTIN_TSV: &str = include_str!("../data/enzymes.tsv");

/// REBASE supplier codes.
pub const SUPPLIERS: &[(char, &str)] = &[
    ('B', "Thermo Fisher Scientific"),
    ('C', "Minotech Biotechnology"),
    ('E', "Agilent Technologies"),
    ('I', "SibEnzyme"),
    ('J', "Nippon Gene"),
    ('K', "Takara Bio"),
    ('M', "Roche"),
    ('N', "New England Biolabs"),
    ('O', "Toyobo"),
    ('Q', "Molecular Biology Resources (CHIMERx)"),
    ('R', "Promega"),
    ('S', "Sigma-Aldrich"),
    ('V', "Vivantis"),
    ('X', "EURx"),
];

/// One restriction enzyme.
#[derive(Clone, Debug, PartialEq)]
pub struct Enzyme {
    pub name: String,
    /// Recognition site, uppercase IUPAC, 5'->3' on the top strand.
    pub site: Vec<u8>,
    /// One or two double-strand cuts as `(top, bottom)` boundaries relative to the site.
    pub cuts: Vec<(i32, i32)>,
    /// REBASE supplier codes (see [`SUPPLIERS`]).
    pub suppliers: String,
}

/// How a cut leaves the ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EndType {
    Blunt,
    FivePrime(u32),
    ThreePrime(u32),
}

impl Enzyme {
    pub fn site_str(&self) -> &str {
        std::str::from_utf8(&self.site).unwrap_or("")
    }

    pub fn site_len(&self) -> usize {
        self.site.len()
    }

    /// Number of specified (non-N) bases in the site: "6-cutter", "8-cutter".
    pub fn specificity(&self) -> usize {
        self.site.iter().filter(|&&b| b != b'N').count()
    }

    /// Whether the site reads the same on both strands.
    pub fn is_palindromic(&self) -> bool {
        alphabet::revcomp(&self.site) == self.site
    }

    /// Whether any cut lies outside the recognition site (Type IIS enzymes such as BsaI).
    pub fn cuts_outside_site(&self) -> bool {
        let n = self.site.len() as i32;
        self.cuts.iter().any(|&(t, b)| t < 0 || b < 0 || t > n || b > n)
    }

    /// Ends left by the (first) cut.
    pub fn end_type(&self) -> EndType {
        let (top, bottom) = self.cuts[0];
        match bottom - top {
            0 => EndType::Blunt,
            d if d > 0 => EndType::FivePrime(d as u32),
            d => EndType::ThreePrime((-d) as u32),
        }
    }

    /// Describe the site with its cuts, REBASE style: `G^AATTC`, `GGTCTC(1/5)`,
    /// `(10/15)ACNNNNGTAYC(12/7)`.
    pub fn site_with_cuts(&self) -> String {
        let n = self.site.len() as i32;
        let site = self.site_str();
        let (top, bottom) = self.cuts[0];
        if self.cuts.len() == 2 {
            let (top2, bottom2) = self.cuts[1];
            return format!("({}/{}){site}({}/{})", -top, -bottom, top2 - n, bottom2 - n);
        }
        if (0..=n).contains(&top) && (0..=n).contains(&bottom) {
            let mut out = String::new();
            out.push_str(&site[..top as usize]);
            out.push('^');
            out.push_str(&site[top as usize..]);
            out
        } else if top > n || bottom > n {
            format!("{site}({}/{})", top - n, bottom - n)
        } else {
            format!("({}/{}){site}", -top, -bottom)
        }
    }

    /// Supplier names.
    pub fn supplier_names(&self) -> Vec<&'static str> {
        self.suppliers.chars().filter_map(|c| SUPPLIERS.iter().find(|(code, _)| *code == c).map(|(_, name)| *name)).collect()
    }
}

/// A searchable set of enzymes.
#[derive(Clone, Debug, Default)]
pub struct EnzymeDb {
    pub enzymes: Vec<Enzyme>,
    /// Lower-case name -> index.
    index: HashMap<String, usize>,
    /// Where the data came from (for example "REBASE version 610").
    pub source: String,
}

impl EnzymeDb {
    /// The built-in database (parsed once, on first use).
    pub fn builtin() -> &'static EnzymeDb {
        static DB: OnceLock<EnzymeDb> = OnceLock::new();
        DB.get_or_init(|| EnzymeDb::from_tsv(BUILTIN_TSV).expect("built-in enzyme table is valid"))
    }

    fn rebuild_index(&mut self) {
        self.index = self.enzymes.iter().enumerate().map(|(ii, e)| (normalize_name(&e.name), ii)).collect();
    }

    pub fn from_enzymes(enzymes: Vec<Enzyme>, source: impl Into<String>) -> EnzymeDb {
        let mut db = EnzymeDb { enzymes, index: HashMap::new(), source: source.into() };
        db.rebuild_index();
        db
    }

    /// Parse SlopGene's TSV enzyme table.
    pub fn from_tsv(text: &str) -> Result<EnzymeDb> {
        let mut enzymes = Vec::new();
        let mut source = String::new();
        for (ii, line) in text.lines().enumerate() {
            if let Some(comment) = line.strip_prefix('#') {
                if source.is_empty()
                    && let Some(at) = comment.find("REBASE version")
                {
                    source = comment[at..].split_whitespace().take(3).collect::<Vec<_>>().join(" ");
                }
                continue;
            }
            if line.trim().is_empty() || line.starts_with("name\t") {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 4 {
                return Err(Error::Parse(format!("enzyme table line {}: too few columns", ii + 1)));
            }
            let int =
                |s: &str| -> Result<i32> { s.trim().parse::<i32>().map_err(|_| Error::Parse(format!("enzyme table line {}: '{s}' is not a number", ii + 1))) };
            let mut cuts = vec![(int(fields[2])?, int(fields[3])?)];
            if fields.len() > 5 && !fields[4].trim().is_empty() && !fields[5].trim().is_empty() {
                cuts.push((int(fields[4])?, int(fields[5])?));
            }
            enzymes.push(Enzyme {
                name: fields[0].trim().to_string(),
                site: fields[1].trim().to_ascii_uppercase().into_bytes(),
                cuts,
                suppliers: fields.get(6).map(|s| s.trim().to_string()).unwrap_or_default(),
            });
        }
        Ok(EnzymeDb::from_enzymes(enzymes, source))
    }

    /// Write the database as SlopGene's TSV table.
    pub fn to_tsv(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!("# SlopGene restriction enzymes: commercially available enzymes from {}\n", self.source));
        out.push_str("# REBASE, the Restriction Enzyme Database, http://rebase.neb.com (Roberts RJ, Vincze T, Posfai J, Macelis D. Nucleic Acids Res.)\n");
        out.push_str("# Cut positions are boundaries counted from the first base of the site (0 = just before it,\n");
        out.push_str("# negative = upstream), on the top strand (top) and bottom strand (bottom), in top-strand coordinates.\n");
        out.push_str("# Regenerate with: slopgene enzymes import-rebase --emboss-e link_emboss_e --emboss-r link_emboss_r\n");
        out.push_str("name\tsite\ttop1\tbottom1\ttop2\tbottom2\tsuppliers\n");
        for e in &self.enzymes {
            let (t1, b1) = e.cuts[0];
            let (t2, b2) = e.cuts.get(1).map(|&(t, b)| (t.to_string(), b.to_string())).unwrap_or_default();
            out.push_str(&format!("{}\t{}\t{t1}\t{b1}\t{t2}\t{b2}\t{}\n", e.name, e.site_str(), e.suppliers));
        }
        out
    }

    /// Build a database from REBASE's EMBOSS files (`link_emboss_e` with sites and cuts,
    /// `link_emboss_r` with suppliers). Only commercially available enzymes with known cuts are
    /// kept, which is what the built-in table holds.
    pub fn from_rebase_emboss(emboss_e: &str, emboss_r: &str) -> Result<EnzymeDb> {
        // REBASE numbers residues ... -2 -1 1 2 ... (there is no zero) and says a cut falls
        // "after residue c". As a boundary that is c itself for c >= 1 and c + 1 below zero.
        let boundary = |c: i32| if c >= 1 { c } else { c + 1 };
        let mut version = String::new();
        // Enzyme name, recognition site, and (top, bottom) cut boundaries.
        type NamedCuts = (String, Vec<u8>, Vec<(i32, i32)>);
        let mut cuts_by_name: Vec<NamedCuts> = Vec::new();
        for line in emboss_e.lines() {
            if line.starts_with('#') {
                if let Some(at) = line.find("REBASE version") {
                    version = line[at..].split_whitespace().take(3).collect::<Vec<_>>().join(" ");
                }
                continue;
            }
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 9 {
                continue;
            }
            let num = |s: &str| s.parse::<i32>().unwrap_or(0);
            let ncuts = num(f[3]);
            if ncuts != 2 && ncuts != 4 {
                continue;
            }
            let mut cuts = vec![(boundary(num(f[5])), boundary(num(f[6])))];
            if ncuts == 4 {
                cuts.push((boundary(num(f[7])), boundary(num(f[8]))));
            }
            cuts_by_name.push((f[0].to_string(), f[1].to_ascii_uppercase().into_bytes(), cuts));
        }
        // emboss_r: blocks end with "//"; line 1 is the name, line 6 the supplier codes.
        let mut suppliers: HashMap<String, String> = HashMap::new();
        let mut block: Vec<&str> = Vec::new();
        for line in emboss_r.lines() {
            if line.starts_with('#') {
                continue;
            }
            if line == "//" {
                if let Some(name) = block.first() {
                    suppliers.insert(name.trim().to_string(), block.get(5).map(|s| s.trim().to_string()).unwrap_or_default());
                }
                block.clear();
            } else {
                block.push(line);
            }
        }
        let mut enzymes: Vec<Enzyme> = cuts_by_name
            .into_iter()
            .filter_map(|(name, site, cuts)| {
                let sup = suppliers.get(&name).cloned().unwrap_or_default();
                (!sup.is_empty()).then_some(Enzyme { name, site, cuts, suppliers: sup })
            })
            .collect();
        enzymes.sort_by_key(|e| e.name.to_ascii_lowercase());
        if enzymes.is_empty() {
            return Err(Error::Parse("no commercially available enzymes found in the REBASE files".into()));
        }
        Ok(EnzymeDb::from_enzymes(enzymes, version))
    }

    pub fn len(&self) -> usize {
        self.enzymes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.enzymes.is_empty()
    }

    /// Look up an enzyme by name, ignoring case and Greek-letter spelling (`BtsαI` = `BtsaI`).
    pub fn get(&self, name: &str) -> Option<&Enzyme> {
        self.index.get(&normalize_name(name)).map(|&ii| &self.enzymes[ii])
    }

    /// Look up several names; unknown names are an error that suggests close matches.
    pub fn resolve(&self, names: &[impl AsRef<str>]) -> Result<Vec<&Enzyme>> {
        let mut out = Vec::with_capacity(names.len());
        for name in names {
            let name = name.as_ref().trim();
            if name.is_empty() {
                continue;
            }
            match self.get(name) {
                Some(e) => out.push(e),
                None => {
                    let suggestions = self.suggest(name, 3);
                    let hint = if suggestions.is_empty() { String::new() } else { format!(" (did you mean {}?)", suggestions.join(", ")) };
                    return Err(Error::NotFound(format!("enzyme '{name}'{hint}")));
                }
            }
        }
        Ok(out)
    }

    /// Names closest to `name` by edit distance.
    pub fn suggest(&self, name: &str, count: usize) -> Vec<String> {
        let target = normalize_name(name);
        let mut scored: Vec<(usize, &str)> =
            self.enzymes.iter().map(|e| (edit_distance(&target, &normalize_name(&e.name)), e.name.as_str())).filter(|(d, _)| *d <= 2).collect();
        scored.sort();
        scored.into_iter().take(count).map(|(_, n)| n.to_string()).collect()
    }

    /// Enzymes in a named set (see [`ENZYME_SETS`]), or a comma-separated list of names.
    pub fn set(&self, spec: &str) -> Result<Vec<&Enzyme>> {
        let key = spec.trim().to_ascii_lowercase();
        let named = ENZYME_SETS.iter().find(|s| s.name == key);
        match named {
            Some(set) => Ok(self.enzymes.iter().filter(|e| (set.filter)(e)).collect()),
            None => {
                let names: Vec<&str> = spec.split([',', ' ', ';']).filter(|s| !s.is_empty()).collect();
                self.resolve(&names)
            }
        }
    }
}

/// Lower-case a name and spell Greek letters out, so `BtsαI` and `BtsaI` match.
fn normalize_name(name: &str) -> String {
    let mut out = String::with_capacity(name.len());
    for ch in name.trim().chars() {
        match ch {
            'α' => out.push('a'),
            'β' => out.push('b'),
            'γ' => out.push('g'),
            'δ' => out.push('d'),
            _ => out.extend(ch.to_lowercase()),
        }
    }
    out
}

fn edit_distance(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for ii in 1..=a.len() {
        let mut row = vec![ii; b.len() + 1];
        for iii in 1..=b.len() {
            let cost = usize::from(a[ii - 1] != b[iii - 1]);
            row[iii] = (prev[iii] + 1).min(row[iii - 1] + 1).min(prev[iii - 1] + cost);
        }
        prev = row;
    }
    prev[b.len()]
}

/// A named group of enzymes.
pub struct EnzymeSet {
    pub name: &'static str,
    pub description: &'static str,
    pub filter: fn(&Enzyme) -> bool,
}

/// Enzymes found in most multiple cloning sites and plasmid workflows.
pub const COMMON_ENZYMES: &[&str] = &[
    "AatII", "AflII", "AgeI", "ApaI", "ApaLI", "AscI", "AvrII", "BamHI", "BbsI", "BglII", "BmtI", "BsaI", "BsiWI", "BsmBI", "BspEI", "BsrGI", "BstBI", "ClaI",
    "DraI", "EagI", "EcoRI", "EcoRV", "Esp3I", "FseI", "HindIII", "HpaI", "KpnI", "MfeI", "MluI", "NcoI", "NdeI", "NheI", "NotI", "NruI", "NsiI", "PacI",
    "PaqCI", "PmeI", "PstI", "PvuI", "PvuII", "SacI", "SacII", "SalI", "SapI", "SbfI", "ScaI", "SfiI", "SmaI", "SpeI", "SphI", "StuI", "SwaI", "XbaI", "XhoI",
    "XmaI",
];

/// Type IIS enzymes used for Golden Gate and related assembly methods.
pub const GOLDEN_GATE_ENZYMES: &[&str] = &["BsaI", "BsmBI", "Esp3I", "BbsI", "BpiI", "SapI", "BspQI", "PaqCI", "AarI", "BtgZI", "BsmAI"];

/// The built-in enzyme sets.
pub static ENZYME_SETS: &[EnzymeSet] = &[
    EnzymeSet { name: "all", description: "every commercially available enzyme", filter: |_| true },
    EnzymeSet { name: "common", description: "enzymes found in common multiple cloning sites", filter: |e| COMMON_ENZYMES.contains(&e.name.as_str()) },
    EnzymeSet { name: "6+", description: "enzymes with at least 6 specified bases in their site", filter: |e| e.specificity() >= 6 },
    EnzymeSet { name: "8+", description: "enzymes with at least 8 specified bases in their site", filter: |e| e.specificity() >= 8 },
    EnzymeSet { name: "neb", description: "enzymes sold by New England Biolabs", filter: |e| e.suppliers.contains('N') },
    EnzymeSet { name: "thermo", description: "enzymes sold by Thermo Fisher Scientific", filter: |e| e.suppliers.contains('B') },
    EnzymeSet {
        name: "golden-gate",
        description: "Type IIS enzymes used for Golden Gate assembly",
        filter: |e| GOLDEN_GATE_ENZYMES.contains(&e.name.as_str()),
    },
    EnzymeSet { name: "type-iis", description: "enzymes that cut outside their recognition site", filter: |e| e.cuts_outside_site() },
    EnzymeSet { name: "blunt", description: "enzymes that leave blunt ends", filter: |e| e.end_type() == EndType::Blunt },
];

/// One double-strand cut in a sequence.
///
/// `pos` is the top-strand cut boundary (the top strand is cut between bases `pos - 1` and
/// `pos`). The bottom strand is cut at `pos + overhang`: a positive overhang leaves a 5'
/// overhang, a negative one a 3' overhang, zero a blunt end. On circular sequences positions
/// are taken modulo the length.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Cut {
    pub pos: usize,
    pub overhang: i32,
}

/// One occurrence of a recognition site.
#[derive(Clone, Debug, PartialEq)]
pub struct Site {
    pub enzyme: String,
    /// First base of the recognition site on the top strand (0-based).
    pub start: usize,
    pub len: usize,
    /// `Forward` when the site reads 5'->3' on the top strand, `Reverse` for the other
    /// orientation of a non-palindromic site.
    pub strand: Strand,
    /// The cuts this occurrence makes. Empty when the cut would fall outside a linear sequence.
    pub cuts: Vec<Cut>,
}

/// Length of the site prefix used to index enzymes for scanning.
const PREFIX: usize = 4;
/// Prefixes that stand for more concrete 4-mers than this are checked at every position.
const MAX_PREFIX_FORMS: usize = 16;

/// One orientation of one enzyme's site.
struct Probe<'a> {
    enzyme: &'a Enzyme,
    pattern: Vec<u8>,
    strand: Strand,
}

/// Every concrete A/C/G/T 4-mer (as a code 0..256) that an IUPAC prefix stands for, or `None`
/// when there are more than [`MAX_PREFIX_FORMS`].
fn prefix_codes(prefix: &[u8]) -> Option<Vec<usize>> {
    let mut codes = vec![0usize];
    for &p in prefix {
        let mask = alphabet::base_mask(p);
        let mut next = Vec::with_capacity(codes.len() * 4);
        for code in &codes {
            // Bit order of the IUPAC mask is A=1, C=2, G=4, T=8; codes use A=0, C=1, G=2, T=3.
            for (bit, value) in [(1u8, 0usize), (2, 1), (4, 2), (8, 3)] {
                if mask & bit != 0 {
                    next.push(code * 4 + value);
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

/// Find every occurrence of the given enzymes' sites. Results are sorted by position, then
/// enzyme name, then orientation.
///
/// Enzymes are indexed by the first four bases of their site (both orientations), so each
/// position of the sequence is compared only with the few enzymes whose site can start there.
/// Positions where the sequence has an ambiguous base fall back to checking every enzyme, so the
/// result is exactly that of comparing every enzyme everywhere ([`find_sites_exhaustive`]).
pub fn find_sites(seq: &[u8], circular: bool, enzymes: &[&Enzyme]) -> Vec<Site> {
    let seq_len = seq.len();
    let mut sites = Vec::new();
    if seq_len == 0 {
        return sites;
    }
    // Probes: each enzyme as written, plus its reverse complement unless it is palindromic.
    let mut probes: Vec<Probe> = Vec::new();
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
    // Index by concrete 4-base prefix; short or very degenerate probes are checked everywhere.
    let mut table: Vec<Vec<usize>> = vec![Vec::new(); 256];
    let mut everywhere: Vec<usize> = Vec::new();
    for (ii, probe) in probes.iter().enumerate() {
        match (probe.pattern.len() >= PREFIX).then(|| prefix_codes(&probe.pattern[..PREFIX])).flatten() {
            Some(codes) => codes.into_iter().for_each(|c| table[c].push(ii)),
            None => everywhere.push(ii),
        }
    }
    let base_at = |pos: usize| seq[pos % seq_len];
    let code_of = |b: u8| match b.to_ascii_uppercase() {
        b'A' => Some(0usize),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    };
    let mut try_probe = |ii: usize, start: usize| {
        let probe = &probes[ii];
        let n = probe.pattern.len();
        if !circular && start + n > seq_len {
            return;
        }
        if !(0..n).all(|k| iupac_matches(probe.pattern[k], base_at(start + k))) {
            return;
        }
        let cuts = probe.enzyme.cuts.iter().filter_map(|&(top, bottom)| place_cut(start, n, top, bottom, probe.strand, seq_len, circular)).collect();
        sites.push(Site { enzyme: probe.enzyme.name.clone(), start, len: n, strand: probe.strand, cuts });
    };
    for start in 0..seq_len {
        let window_fits = circular || start + PREFIX <= seq_len;
        let code = if window_fits { (0..PREFIX).try_fold(0usize, |acc, k| code_of(base_at(start + k)).map(|c| acc * 4 + c)) } else { None };
        match code {
            Some(c) => {
                for &ii in &table[c] {
                    try_probe(ii, start);
                }
                for &ii in &everywhere {
                    try_probe(ii, start);
                }
            }
            // An ambiguous base in the window: compare every probe the slow way.
            None if window_fits => (0..probes.len()).for_each(|ii| try_probe(ii, start)),
            // Near the end of a linear sequence only probes shorter than the prefix can fit.
            None => everywhere.iter().for_each(|&ii| try_probe(ii, start)),
        }
    }
    sites.sort_by(|a, b| a.start.cmp(&b.start).then_with(|| a.enzyme.cmp(&b.enzyme)).then_with(|| (a.strand as u8).cmp(&(b.strand as u8))));
    sites
}

/// Reference implementation of [`find_sites`]: every enzyme compared at every position.
pub fn find_sites_exhaustive(seq: &[u8], circular: bool, enzymes: &[&Enzyme]) -> Vec<Site> {
    let seq_len = seq.len();
    let mut sites = Vec::new();
    if seq_len == 0 {
        return sites;
    }
    for enzyme in enzymes {
        let n = enzyme.site.len();
        if n == 0 || (n > seq_len && !circular) {
            continue;
        }
        let reverse_site = alphabet::revcomp(&enzyme.site);
        let palindromic = reverse_site == enzyme.site;
        // On a circular molecule a site may start near the end and continue at the start.
        let last_start = if circular { seq_len } else { seq_len - n + 1 };
        let base_at = |pos: usize| seq[pos % seq_len];
        for start in 0..last_start {
            // ii = orientation: 0 = the site as written, 1 = its reverse complement.
            for ii in 0..if palindromic { 1 } else { 2 } {
                let pattern = if ii == 0 { &enzyme.site } else { &reverse_site };
                if !(0..n).all(|iii| iupac_matches(pattern[iii], base_at(start + iii))) {
                    continue;
                }
                let strand = if ii == 0 { Strand::Forward } else { Strand::Reverse };
                let cuts = enzyme.cuts.iter().filter_map(|&(top, bottom)| place_cut(start, n, top, bottom, strand, seq_len, circular)).collect();
                sites.push(Site { enzyme: enzyme.name.clone(), start, len: n, strand, cuts });
            }
        }
    }
    sites.sort_by(|a, b| a.start.cmp(&b.start).then_with(|| a.enzyme.cmp(&b.enzyme)).then_with(|| (a.strand as u8).cmp(&(b.strand as u8))));
    sites
}

/// Turn an enzyme's relative cut into a position in the sequence.
fn place_cut(start: usize, n: usize, top: i32, bottom: i32, strand: Strand, seq_len: usize, circular: bool) -> Option<Cut> {
    // A site on the bottom strand is the enzyme seen from the other side: its top-strand cut
    // lands on our bottom strand and the boundaries mirror around the site.
    let (top_abs, bottom_abs) = match strand {
        Strand::Reverse => (start as i64 + n as i64 - bottom as i64, start as i64 + n as i64 - top as i64),
        _ => (start as i64 + top as i64, start as i64 + bottom as i64),
    };
    let overhang = (bottom_abs - top_abs) as i32;
    let len = seq_len as i64;
    if circular {
        Some(Cut { pos: top_abs.rem_euclid(len) as usize, overhang })
    } else {
        // Both strands must be cut inside the molecule, and not at its very ends.
        let lo = top_abs.min(bottom_abs);
        let hi = top_abs.max(bottom_abs);
        (lo > 0 && hi < len).then_some(Cut { pos: top_abs as usize, overhang })
    }
}

/// The sites that cut, grouped by enzyme name and kept in position order. Sites whose cuts fall
/// outside a linear molecule are left out.
pub fn cutting_sites_by_enzyme(sites: &[Site]) -> HashMap<&str, Vec<&Site>> {
    let mut groups: HashMap<&str, Vec<&Site>> = HashMap::new();
    for site in sites.iter().filter(|s| !s.cuts.is_empty()) {
        groups.entry(site.enzyme.as_str()).or_default().push(site);
    }
    groups
}

/// How many times each enzyme cuts (enzymes that do not cut are included with 0).
pub fn cut_counts<'a>(seq: &[u8], circular: bool, enzymes: &[&'a Enzyme]) -> Vec<(&'a Enzyme, usize)> {
    let sites = find_sites(seq, circular, enzymes);
    let groups = cutting_sites_by_enzyme(&sites);
    enzymes.iter().map(|e| (*e, groups.get(e.name.as_str()).map_or(0, |g| g.len()))).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn db() -> &'static EnzymeDb {
        EnzymeDb::builtin()
    }

    #[test]
    fn builtin_table_loads() {
        let db = db();
        assert!(db.len() > 500, "only {} enzymes", db.len());
        assert!(db.source.starts_with("REBASE version"), "{}", db.source);
        let ecori = db.get("ecori").unwrap();
        assert_eq!(ecori.site_str(), "GAATTC");
        assert_eq!(ecori.cuts, vec![(1, 5)]);
        assert_eq!(ecori.end_type(), EndType::FivePrime(4));
        assert_eq!(ecori.site_with_cuts(), "G^AATTC");
        assert_eq!(db.get("KpnI").unwrap().end_type(), EndType::ThreePrime(4));
        assert_eq!(db.get("BsaI").unwrap().site_with_cuts(), "GGTCTC(1/5)");
        assert_eq!(db.get("BaeI").unwrap().site_with_cuts(), "(10/15)ACNNNNGTAYC(12/7)");
        assert!(db.get("BsaI").unwrap().cuts_outside_site());
        assert!(db.resolve(&["EcoRl"]).unwrap_err().to_string().contains("EcoRI"));
    }

    #[test]
    fn tsv_round_trip() {
        let db = db();
        let again = EnzymeDb::from_tsv(&db.to_tsv()).unwrap();
        assert_eq!(again.enzymes, db.enzymes);
    }

    #[test]
    fn finds_palindromic_and_type_iis_sites() {
        let db = db();
        let ecori = db.get("EcoRI").unwrap();
        let bsai = db.get("BsaI").unwrap();
        // EcoRI at 10; BsaI forward at 30 (cuts 1/5 downstream) and reverse at 50.
        let mut seq = vec![b'A'; 80];
        seq[10..16].copy_from_slice(b"GAATTC");
        seq[30..36].copy_from_slice(b"GGTCTC");
        seq[50..56].copy_from_slice(b"GAGACC");
        let sites = find_sites(&seq, false, &[ecori, bsai]);
        assert_eq!(sites.len(), 3);
        assert_eq!(sites[0].cuts, vec![Cut { pos: 11, overhang: 4 }]);
        assert_eq!(sites[1].strand, Strand::Forward);
        assert_eq!(sites[1].cuts, vec![Cut { pos: 37, overhang: 4 }]);
        // Reverse BsaI: cuts upstream of the site, top at 50-5, bottom at 50-1.
        assert_eq!(sites[2].strand, Strand::Reverse);
        assert_eq!(sites[2].cuts, vec![Cut { pos: 45, overhang: 4 }]);
    }

    #[test]
    fn site_across_origin() {
        let ecori = db().get("EcoRI").unwrap();
        let mut seq = b"ATTCAAAAAAAAAAGA".to_vec(); // GA | ATTC split by the origin
        let sites = find_sites(&seq, true, &[ecori]);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].start, 14);
        assert_eq!(sites[0].cuts, vec![Cut { pos: 15, overhang: 4 }]);
        assert!(find_sites(&seq, false, &[ecori]).is_empty());
        seq.clear();
    }

    #[test]
    fn indexed_search_matches_exhaustive_search() {
        // Random sequences with a sprinkling of ambiguous bases, circular and linear, against
        // every built-in enzyme: the fast scan must find exactly what the slow scan finds.
        let all: Vec<&Enzyme> = db().enzymes.iter().collect();
        let mut x = 12345u64;
        for ii in 0..6 {
            let len = 400 + ii * 350;
            let seq: Vec<u8> = (0..len)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    let r = (x >> 33) as usize % 100;
                    if r == 0 {
                        b'N'
                    } else if r == 1 {
                        b'R'
                    } else {
                        b"ACGT"[r % 4]
                    }
                })
                .collect();
            for circular in [true, false] {
                let fast = find_sites(&seq, circular, &all);
                let slow = find_sites_exhaustive(&seq, circular, &all);
                assert_eq!(fast, slow, "length {len}, circular {circular}");
            }
        }
    }

    #[test]
    fn named_sets() {
        let db = db();
        assert!(db.set("common").unwrap().iter().any(|e| e.name == "EcoRI"));
        assert!(db.set("golden-gate").unwrap().iter().any(|e| e.name == "BsaI"));
        assert_eq!(db.set("EcoRI, BamHI").unwrap().len(), 2);
        assert!(db.set("6+").unwrap().iter().all(|e| e.specificity() >= 6));
    }
}
