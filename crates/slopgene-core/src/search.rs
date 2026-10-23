//! Searching sequences: nucleotide motifs (IUPAC codes, mismatches, both strands) and protein
//! motifs in the six reading frames.

use crate::alphabet::{self, iupac_matches};
use crate::model::Strand;
use crate::span::Span;
use crate::translate::GeneticCode;

/// One match of a search.
#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    /// Matched bases in top-strand coordinates.
    pub span: Span,
    pub strand: Strand,
    pub mismatches: usize,
}

/// Find `pattern` (IUPAC codes allowed in the pattern) in `seq`, allowing up to
/// `max_mismatches` substitutions. Matches on both strands are reported when `both_strands`.
/// On circular sequences matches may run through the origin.
pub fn find_dna(seq: &[u8], circular: bool, pattern: &[u8], max_mismatches: usize, both_strands: bool) -> Vec<Hit> {
    let n = pattern.len();
    let seq_len = seq.len();
    let mut hits = Vec::new();
    if n == 0 || seq_len == 0 || (n > seq_len && !circular) || n > seq_len.max(1) * 2 {
        return hits;
    }
    let forward = alphabet::upper(pattern);
    let reverse = alphabet::revcomp(&forward);
    let palindromic = forward == reverse;
    let last_start = if circular { seq_len } else { seq_len - n + 1 };
    // ii = 0 searches the pattern as given (forward strand), ii = 1 its reverse complement.
    for ii in 0..if both_strands && !palindromic { 2 } else { 1 } {
        let probe = if ii == 0 { &forward } else { &reverse };
        for start in 0..last_start {
            let mut mismatches = 0;
            let mut ok = true;
            for iii in 0..n {
                if !iupac_matches(probe[iii], seq[(start + iii) % seq_len]) {
                    mismatches += 1;
                    if mismatches > max_mismatches {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                let strand = if ii == 0 { Strand::Forward } else { Strand::Reverse };
                hits.push(Hit { span: Span::new(start, start + n), strand, mismatches });
            }
        }
    }
    hits.sort_by_key(|h| (h.span.start, h.strand as u8));
    hits
}

/// Find a protein motif (one-letter codes, `X` matches anything) in all six reading frames.
/// Spans are the coding bases in top-strand coordinates.
pub fn find_protein(seq: &[u8], circular: bool, motif: &[u8], code: &GeneticCode) -> Vec<Hit> {
    let motif = motif.to_ascii_uppercase();
    let m = motif.len();
    let seq_len = seq.len();
    let mut hits = Vec::new();
    if m == 0 || seq_len < 3 * m {
        return hits;
    }
    let upper = alphabet::upper(seq);
    // On circular sequences, append enough of the start to see codons through the origin.
    let extended = |s: &[u8]| -> Vec<u8> {
        let mut v = s.to_vec();
        if circular {
            v.extend_from_slice(&s[..(3 * m + 2).min(s.len())]);
        }
        v
    };
    for (strand, strand_seq) in [(Strand::Forward, upper.clone()), (Strand::Reverse, alphabet::revcomp(&upper))] {
        let ext = extended(&strand_seq);
        for frame in 0..3 {
            let protein: Vec<u8> = ext[frame..].as_chunks::<3>().0.iter().map(|c| code.translate_codon(c)).collect();
            if protein.len() < m {
                continue;
            }
            for ii in 0..=protein.len() - m {
                let matches = (0..m).all(|iii| motif[iii] == b'X' || motif[iii] == protein[ii + iii]);
                if !matches {
                    continue;
                }
                let start_on_strand = frame + 3 * ii;
                if start_on_strand >= seq_len {
                    continue;
                }
                let len = 3 * m;
                let span = match strand {
                    Strand::Reverse => {
                        // Position p on the reverse complement is top-strand base L - 1 - p, so
                        // the stretch [s, s + len) there is [L - s - len, L - s) on the top strand.
                        let end_on_strand = (start_on_strand + len) as i64;
                        let start = (seq_len as i64 - end_on_strand).rem_euclid(seq_len as i64) as usize;
                        Span::new(start, start + len)
                    }
                    _ => Span::new(start_on_strand, start_on_strand + len),
                };
                if !circular && span.end > seq_len {
                    continue;
                }
                hits.push(Hit { span, strand, mismatches: 0 });
            }
        }
    }
    hits.sort_by_key(|h| (h.span.start, h.strand as u8));
    hits.dedup();
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_search_with_mismatches_and_strands() {
        let seq = b"AAAAGAATTCAAAAGTATTCAAAACCTTAAAA";
        let exact = find_dna(seq, false, b"GAATTC", 0, true);
        assert_eq!(exact.len(), 1);
        assert_eq!(exact[0].span, Span::new(4, 10));
        let fuzzy = find_dna(seq, false, b"GAATTC", 1, true);
        assert!(fuzzy.iter().any(|h| h.span.start == 14 && h.mismatches == 1));
        // Non-palindromic pattern found on the reverse strand.
        let rev = find_dna(seq, false, b"AAGGTT", 0, true);
        assert_eq!(rev.len(), 1);
        assert_eq!(rev[0].strand, Strand::Reverse);
        assert_eq!(rev[0].span, Span::new(22, 28));
        // IUPAC pattern.
        assert_eq!(find_dna(seq, false, b"GWATTC", 0, false).len(), 2);
    }

    #[test]
    fn protein_search_forward_and_reverse() {
        let code = GeneticCode::standard();
        let coding = b"ATGGCTAGCAAGGGC"; // M A S K G
        let mut seq = b"TT".to_vec();
        seq.extend_from_slice(coding);
        seq.extend_from_slice(b"TTT");
        seq.extend(alphabet::revcomp(coding));
        let hits = find_protein(&seq, false, b"ASK", code);
        assert_eq!(hits.len(), 2, "{hits:?}");
        assert_eq!(hits[0].span, Span::new(5, 14));
        assert_eq!(hits[1].strand, Strand::Reverse);
        assert_eq!(alphabet::revcomp(&hits[1].span.extract(&seq)), b"GCTAGCAAG".to_vec());
    }
}
