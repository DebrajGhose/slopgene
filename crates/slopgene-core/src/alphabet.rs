//! Nucleotide alphabets, including the IUPAC ambiguity codes.
//!
//! Sequences are stored as ASCII bytes. Every function here keeps the case of its input, so a
//! lowercase stretch stays lowercase after complementing.

/// Build the complement lookup table at compile time.
///
/// Every byte maps to itself unless it is a nucleotide code, so complementing text that contains
/// gaps or unknown symbols leaves those symbols alone.
const fn build_complement_table() -> [u8; 256] {
    let mut table = [0u8; 256];
    let mut ii = 0;
    while ii < 256 {
        table[ii] = ii as u8;
        ii += 1;
    }

    // Uppercase pairs; the lowercase entries are filled in from the same list (ASCII + 32).
    // U complements to A because RNA is read here as DNA. Use `revcomp_rna` when U is wanted back.
    let pairs: [(u8, u8); 16] = [
        (b'A', b'T'),
        (b'T', b'A'),
        (b'U', b'A'),
        (b'G', b'C'),
        (b'C', b'G'),
        (b'R', b'Y'),
        (b'Y', b'R'),
        (b'S', b'S'),
        (b'W', b'W'),
        (b'K', b'M'),
        (b'M', b'K'),
        (b'B', b'V'),
        (b'V', b'B'),
        (b'D', b'H'),
        (b'H', b'D'),
        (b'N', b'N'),
    ];
    let mut ii = 0;
    while ii < pairs.len() {
        let (base, comp) = pairs[ii];
        table[base as usize] = comp;
        table[(base + 32) as usize] = comp + 32;
        ii += 1;
    }
    table
}

/// Build the IUPAC bit-mask table: A = 1, C = 2, G = 4, T/U = 8, and ambiguity codes are the
/// union of the bases they stand for. Anything that is not a nucleotide code maps to 0.
const fn build_mask_table() -> [u8; 256] {
    let mut table = [0u8; 256];
    let codes: [(u8, u8); 16] = [
        (b'A', 1),
        (b'C', 2),
        (b'G', 4),
        (b'T', 8),
        (b'U', 8),
        (b'R', 1 | 4),
        (b'Y', 2 | 8),
        (b'S', 2 | 4),
        (b'W', 1 | 8),
        (b'K', 4 | 8),
        (b'M', 1 | 2),
        (b'B', 2 | 4 | 8),
        (b'D', 1 | 4 | 8),
        (b'H', 1 | 2 | 8),
        (b'V', 1 | 2 | 4),
        (b'N', 1 | 2 | 4 | 8),
    ];
    let mut ii = 0;
    while ii < codes.len() {
        let (code, mask) = codes[ii];
        table[code as usize] = mask;
        table[(code + 32) as usize] = mask;
        ii += 1;
    }
    table
}

static COMPLEMENT: [u8; 256] = build_complement_table();
static MASK: [u8; 256] = build_mask_table();

/// Complement of one nucleotide code (case kept). Non-nucleotide bytes come back unchanged.
#[inline]
pub fn complement_base(base: u8) -> u8 {
    COMPLEMENT[base as usize]
}

/// Complement without reversing.
pub fn complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| complement_base(b)).collect()
}

/// Reverse complement of a DNA sequence (U is complemented to A).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// Whether two circular molecules are the same sequence, whatever base each was opened at and
/// whichever strand each was written from.
pub fn same_circle(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let a = a.to_ascii_uppercase();
    let doubled = [b.to_ascii_uppercase(), b.to_ascii_uppercase()].concat();
    let reverse = revcomp(&a);
    doubled.windows(a.len()).any(|w| w == a.as_slice() || w == reverse.as_slice())
}

/// Reverse complement of an RNA sequence: like [`revcomp`], but writes U where DNA would have T.
pub fn revcomp_rna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match complement_base(b) {
            b'T' => b'U',
            b't' => b'u',
            other => other,
        })
        .collect()
}

/// Reverse complement in place.
pub fn revcomp_in_place(seq: &mut [u8]) {
    seq.reverse();
    for base in seq.iter_mut() {
        *base = complement_base(*base);
    }
}

/// IUPAC bit mask of a nucleotide code (0 when the byte is not a nucleotide code).
#[inline]
pub fn base_mask(base: u8) -> u8 {
    MASK[base as usize]
}

/// Whether sequence base `base` is matched by pattern code `pattern`.
///
/// A base matches when every nucleotide it can stand for is allowed by the pattern. So the
/// pattern `N` matches anything, but the sequence base `N` only matches a pattern `N`: an unknown
/// base in a sequence is never reported as, say, an EcoRI site.
#[inline]
pub fn iupac_matches(pattern: u8, base: u8) -> bool {
    let seq_mask = MASK[base as usize];
    seq_mask != 0 && (seq_mask & !MASK[pattern as usize]) == 0
}

/// Whether the byte is one of the 16 IUPAC nucleotide codes (either case).
#[inline]
pub fn is_nucleotide(base: u8) -> bool {
    MASK[base as usize] != 0
}

/// Whether the byte is A, C, G, T or U (either case).
#[inline]
pub fn is_unambiguous(base: u8) -> bool {
    matches!(MASK[base as usize], 1 | 2 | 4 | 8)
}

/// Whether two codes could be the same base (their masks overlap).
#[inline]
pub fn could_pair_equal(a: u8, b: u8) -> bool {
    (MASK[a as usize] & MASK[b as usize]) != 0
}

/// Turn free text into a clean nucleotide sequence.
///
/// Whitespace, digits and the separators people paste in (`-`, `.`, `/`, `|`, `*`) are dropped.
/// Anything else that is not an IUPAC code is an error, so typos are caught instead of being
/// silently stored. Case is kept.
pub fn clean_nucleotides(text: &str) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(text.len());
    for (ii, ch) in text.char_indices() {
        if ch.is_whitespace() || ch.is_ascii_digit() || matches!(ch, '-' | '.' | '/' | '|' | '*') {
            continue;
        }
        if ch.is_ascii() && is_nucleotide(ch as u8) {
            out.push(ch as u8);
        } else {
            return Err(format!("'{ch}' at character {} is not a nucleotide code", ii + 1));
        }
    }
    Ok(out)
}

/// Turn free text into a clean protein sequence: letters and `*` (stop) are kept, case is
/// normalised to uppercase, whitespace and digits are dropped.
pub fn clean_protein(text: &str) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(text.len());
    for (ii, ch) in text.char_indices() {
        if ch.is_whitespace() || ch.is_ascii_digit() {
            continue;
        }
        if ch.is_ascii_alphabetic() || ch == '*' {
            out.push(ch.to_ascii_uppercase() as u8);
        } else {
            return Err(format!("'{ch}' at character {} is not an amino-acid code", ii + 1));
        }
    }
    Ok(out)
}

/// Fraction of G + C among the unambiguous bases (S counts as G/C, W as A/T).
/// Returns 0 for a sequence with no countable bases.
pub fn gc_fraction(seq: &[u8]) -> f64 {
    let mut gc = 0usize;
    let mut total = 0usize;
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'G' | b'C' | b'S' => {
                gc += 1;
                total += 1;
            }
            b'A' | b'T' | b'U' | b'W' => total += 1,
            _ => {}
        }
    }
    if total == 0 { 0.0 } else { gc as f64 / total as f64 }
}

/// Count of each of A, C, G, T (U counted as T); anything else is counted in the fifth slot.
pub fn base_counts(seq: &[u8]) -> [usize; 5] {
    let mut counts = [0usize; 5];
    for &b in seq {
        let slot = match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' | b'U' => 3,
            _ => 4,
        };
        counts[slot] += 1;
    }
    counts
}

/// Uppercase copy of a sequence.
pub fn upper(seq: &[u8]) -> Vec<u8> {
    seq.to_ascii_uppercase()
}

/// Three-letter name of an amino acid one-letter code ("Ter" for the stop `*`).
pub fn amino_acid_three_letter(aa: u8) -> &'static str {
    match aa.to_ascii_uppercase() {
        b'A' => "Ala",
        b'R' => "Arg",
        b'N' => "Asn",
        b'D' => "Asp",
        b'C' => "Cys",
        b'E' => "Glu",
        b'Q' => "Gln",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'L' => "Leu",
        b'K' => "Lys",
        b'M' => "Met",
        b'F' => "Phe",
        b'P' => "Pro",
        b'S' => "Ser",
        b'T' => "Thr",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'V' => "Val",
        b'U' => "Sec",
        b'O' => "Pyl",
        b'B' => "Asx",
        b'Z' => "Glx",
        b'J' => "Xle",
        b'*' => "Ter",
        _ => "Xaa",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_keeps_case_and_ambiguity() {
        assert_eq!(revcomp(b"GAATTCnRY"), b"RYnGAATTC".to_vec());
        assert_eq!(revcomp(b"ACGTacgt"), b"acgtACGT".to_vec());
        assert_eq!(revcomp_rna(b"AUGC"), b"GCAU".to_vec());
    }

    #[test]
    fn iupac_matching_rules() {
        assert!(iupac_matches(b'N', b'A'));
        assert!(iupac_matches(b'R', b'g'));
        assert!(!iupac_matches(b'R', b'C'));
        // An unknown base in the sequence only matches N in the pattern.
        assert!(!iupac_matches(b'A', b'N'));
        assert!(iupac_matches(b'N', b'N'));
        assert!(!iupac_matches(b'A', b'-'));
    }

    #[test]
    fn cleaning_rejects_typos() {
        assert_eq!(clean_nucleotides("  gaa ttc 12\nAC").unwrap(), b"gaattcAC".to_vec());
        assert!(clean_nucleotides("GAXTTC").is_err());
    }

    #[test]
    fn gc_content() {
        assert!((gc_fraction(b"GGCCAATT") - 0.5).abs() < 1e-12);
        assert_eq!(gc_fraction(b"NNNN"), 0.0);
    }
}
