//! Translation with the NCBI genetic codes, and open-reading-frame detection.

use crate::alphabet::{self, base_mask};
use crate::model::{Document, Feature, Strand};
use crate::span::Span;

/// One NCBI genetic code. Both strings list the 64 codons in NCBI order: the first, second and
/// third bases each run through T, C, A, G.
pub struct GeneticCode {
    pub id: u32,
    pub name: &'static str,
    /// Amino acid for each codon (`*` = stop).
    pub amino_acids: &'static [u8; 64],
    /// `M` marks codons that may start translation, `*` stops, `-` neither.
    pub starts: &'static [u8; 64],
}

/// The NCBI genetic codes, transcribed from NCBI's `gc.prt`
/// (https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt).
pub static GENETIC_CODES: &[GeneticCode] = &[
    GeneticCode {
        id: 1,
        name: "Standard",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M---------------M----------------------------",
    },
    GeneticCode {
        id: 2,
        name: "Vertebrate Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        starts: b"----------**--------------------MMMM----------**---M------------",
    },
    GeneticCode {
        id: 3,
        name: "Yeast Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**----------------------MM---------------M------------",
    },
    GeneticCode {
        id: 4,
        name: "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--MM------**-------M------------MMMM---------------M------------",
    },
    GeneticCode {
        id: 5,
        name: "Invertebrate Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        starts: b"---M------**--------------------MMMM---------------M------------",
    },
    GeneticCode {
        id: 6,
        name: "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
        amino_acids: b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    GeneticCode {
        id: 9,
        name: "Echinoderm Mitochondrial; Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    GeneticCode {
        id: 10,
        name: "Euplotid Nuclear",
        amino_acids: b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    GeneticCode {
        id: 11,
        name: "Bacterial, Archaeal and Plant Plastid",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M------------MMMM---------------M------------",
    },
    GeneticCode {
        id: 12,
        name: "Alternative Yeast Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    GeneticCode {
        id: 13,
        name: "Ascidian Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        starts: b"---M------**----------------------MM---------------M------------",
    },
    GeneticCode {
        id: 14,
        name: "Alternative Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"-----------*-----------------------M----------------------------",
    },
    GeneticCode {
        id: 15,
        name: "Blepharisma Macronuclear",
        amino_acids: b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    GeneticCode {
        id: 16,
        name: "Chlorophycean Mitochondrial",
        amino_acids: b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    GeneticCode {
        id: 21,
        name: "Trematode Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    GeneticCode {
        id: 22,
        name: "Scenedesmus obliquus Mitochondrial",
        amino_acids: b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"------*---*---*--------------------M----------------------------",
    },
    GeneticCode {
        id: 23,
        name: "Thraustochytrium Mitochondrial",
        amino_acids: b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--*-------**--*-----------------M--M---------------M------------",
    },
    GeneticCode {
        id: 24,
        name: "Rhabdopleuridae Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M------**-------M---------------M---------------M------------",
    },
    GeneticCode {
        id: 25,
        name: "Candidate Division SR1 and Gracilibacteria",
        amino_acids: b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**-----------------------M---------------M------------",
    },
    GeneticCode {
        id: 26,
        name: "Pachysolen tannophilus Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    GeneticCode {
        id: 27,
        name: "Karyorelict Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    GeneticCode {
        id: 28,
        name: "Condylostoma Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*--------------------M----------------------------",
    },
    GeneticCode {
        id: 29,
        name: "Mesodinium Nuclear",
        amino_acids: b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    GeneticCode {
        id: 30,
        name: "Peritrich Nuclear",
        amino_acids: b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    GeneticCode {
        id: 31,
        name: "Blastocrithidia Nuclear",
        amino_acids: b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    GeneticCode {
        id: 32,
        name: "Balanophoraceae Plastid",
        amino_acids: b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------*---*----M------------MMMM---------------M------------",
    },
    GeneticCode {
        id: 33,
        name: "Cephalodiscidae Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M-------*-------M---------------M---------------M------------",
    },
];

impl GeneticCode {
    /// The genetic code with this NCBI id.
    pub fn get(id: u32) -> Option<&'static GeneticCode> {
        GENETIC_CODES.iter().find(|c| c.id == id)
    }

    /// The standard code (NCBI table 1).
    pub fn standard() -> &'static GeneticCode {
        &GENETIC_CODES[0]
    }

    /// Amino acid for a codon. Ambiguous codons are translated when every base they could
    /// stand for gives the same amino acid (GCN is alanine); otherwise the result is `X`.
    pub fn translate_codon(&self, codon: &[u8]) -> u8 {
        self.resolve(codon, |index| self.amino_acids[index]).unwrap_or(b'X')
    }

    /// Whether the codon is a start codon in this code (alternative starts included).
    pub fn is_start(&self, codon: &[u8]) -> bool {
        self.resolve(codon, |index| self.starts[index]) == Some(b'M')
    }

    /// Whether the codon is a stop codon (all of its possible readings are stops).
    pub fn is_stop(&self, codon: &[u8]) -> bool {
        self.resolve(codon, |index| self.amino_acids[index]) == Some(b'*')
    }

    /// Apply `lookup` to every unambiguous codon the (possibly ambiguous) codon stands for, and
    /// return the answer if they all agree.
    fn resolve(&self, codon: &[u8], lookup: impl Fn(usize) -> u8) -> Option<u8> {
        if codon.len() != 3 {
            return None;
        }
        // Position of each base in NCBI's T, C, A, G order, as bits of the IUPAC mask.
        const ORDER: [(u8, usize); 4] = [(8, 0), (2, 1), (1, 2), (4, 3)];
        let masks = [base_mask(codon[0]), base_mask(codon[1]), base_mask(codon[2])];
        if masks.contains(&0) {
            return None;
        }
        let mut answer: Option<u8> = None;
        for &(bit1, i1) in &ORDER {
            if masks[0] & bit1 == 0 {
                continue;
            }
            for &(bit2, i2) in &ORDER {
                if masks[1] & bit2 == 0 {
                    continue;
                }
                for &(bit3, i3) in &ORDER {
                    if masks[2] & bit3 == 0 {
                        continue;
                    }
                    let value = lookup(16 * i1 + 4 * i2 + i3);
                    match answer {
                        None => answer = Some(value),
                        Some(prev) if prev != value => return None,
                        _ => {}
                    }
                }
            }
        }
        answer
    }
}

/// Translate `seq` from its first base, codon by codon. A trailing partial codon is dropped.
/// With `to_stop` the translation ends before the first stop codon.
pub fn translate(seq: &[u8], code: &GeneticCode, to_stop: bool) -> Vec<u8> {
    let mut protein = Vec::with_capacity(seq.len() / 3);
    for codon in seq.as_chunks::<3>().0 {
        let aa = code.translate_codon(codon);
        if to_stop && aa == b'*' {
            break;
        }
        protein.push(aa);
    }
    protein
}

/// Translate a coding feature, honouring its `/codon_start` and `/transl_table` qualifiers.
/// The first codon is translated as M when it is an alternative start codon, as ribosomes do.
pub fn translate_feature(feature: &Feature, doc: &Document) -> Vec<u8> {
    let table =
        feature.qualifier_text("transl_table").and_then(|t| t.trim().parse::<u32>().ok()).and_then(GeneticCode::get).unwrap_or_else(GeneticCode::standard);
    let codon_start = feature.qualifier_text("codon_start").and_then(|t| t.trim().parse::<usize>().ok()).unwrap_or(1).clamp(1, 3);
    let coding = feature.extract(&doc.seq);
    let coding = if coding.len() >= codon_start - 1 { &coding[codon_start - 1..] } else { &coding[..0] };
    let mut protein = translate(coding, table, false);
    if codon_start == 1 && coding.len() >= 3 && table.is_start(&coding[..3]) {
        protein[0] = b'M';
    }
    protein
}

/// Which codons may start an open reading frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum StartCodons {
    /// ATG only (SnapGene's default).
    #[default]
    Atg,
    /// Every start codon of the genetic code (for example GTG and TTG in bacteria).
    Alternative,
    /// Any codon after a stop: stop-to-stop frames.
    Any,
}

/// Options for [`find_orfs`].
#[derive(Clone, Copy, Debug)]
pub struct OrfOptions {
    /// Minimum length in codons, not counting the stop codon.
    pub min_codons: usize,
    pub starts: StartCodons,
    pub genetic_code: u32,
    pub both_strands: bool,
}

impl Default for OrfOptions {
    fn default() -> Self {
        // SnapGene shows ORFs of at least 75 codons by default.
        OrfOptions { min_codons: 75, starts: StartCodons::Atg, genetic_code: 1, both_strands: true }
    }
}

/// An open reading frame: start codon through stop codon (inclusive).
#[derive(Clone, Debug, PartialEq)]
pub struct Orf {
    pub span: Span,
    pub strand: Strand,
    /// Reading frame 1-3 on its strand, counted from the strand's first base.
    pub frame: u8,
    /// Number of codons before the stop codon.
    pub codons: usize,
}

/// Find open reading frames on both strands. On circular sequences ORFs may run through the
/// origin. For nested starts in one frame only the longest ORF (first start after the
/// previous stop) is reported. Results are sorted by start position.
pub fn find_orfs(seq: &[u8], circular: bool, options: &OrfOptions) -> Vec<Orf> {
    let code = GeneticCode::get(options.genetic_code).unwrap_or_else(GeneticCode::standard);
    let seq_len = seq.len();
    let mut orfs = Vec::new();
    if seq_len < 3 {
        return orfs;
    }
    let upper = alphabet::upper(seq);
    let strands: Vec<(Strand, Vec<u8>)> = if options.both_strands {
        vec![(Strand::Forward, upper.clone()), (Strand::Reverse, alphabet::revcomp(&upper))]
    } else {
        vec![(Strand::Forward, upper)]
    };

    for (strand, strand_seq) in strands {
        // On a circular molecule scan three turns and keep ORFs whose start lies in the middle
        // turn: by then the upstream context (the previous stop) is known, and no ORF is
        // longer than one turn.
        let (scan, window) = if circular {
            let mut tripled = Vec::with_capacity(seq_len * 3);
            for _ in 0..3 {
                tripled.extend_from_slice(&strand_seq);
            }
            (tripled, seq_len..2 * seq_len)
        } else {
            (strand_seq.clone(), 0..seq_len)
        };
        let mut seen = std::collections::HashSet::new();
        for frame in 0..3usize {
            // `open` is the position of the first start codon since the last stop.
            let mut open: Option<usize> = if options.starts == StartCodons::Any { Some(frame) } else { None };
            let mut pos = frame;
            while pos + 3 <= scan.len() {
                let codon = &scan[pos..pos + 3];
                if code.is_stop(codon) {
                    if let Some(start) = open.take() {
                        let end = pos + 3;
                        let codons = (end - start) / 3 - 1;
                        if window.contains(&start) && end - start <= seq_len && codons >= options.min_codons {
                            let local = start % seq_len;
                            if seen.insert(local) {
                                orfs.push(make_orf(local, end - start, strand, seq_len));
                            }
                        }
                    }
                    if options.starts == StartCodons::Any {
                        open = Some(pos + 3);
                    }
                } else if open.is_none() {
                    let is_start = match options.starts {
                        StartCodons::Atg => codon == b"ATG",
                        StartCodons::Alternative => code.is_start(codon),
                        StartCodons::Any => true,
                    };
                    if is_start {
                        open = Some(pos);
                    }
                }
                pos += 3;
            }
        }
    }
    orfs.sort_by_key(|o| (o.span.start, o.span.end));
    orfs
}

/// Build an ORF record, converting reverse-strand positions back to top-strand coordinates.
fn make_orf(start_on_strand: usize, len: usize, strand: Strand, seq_len: usize) -> Orf {
    let span_on_strand = Span::new(start_on_strand, start_on_strand + len);
    let span = if strand == Strand::Reverse {
        // Position p on the reverse complement is base L - 1 - p on the top strand.
        let start = (2 * seq_len - span_on_strand.end) % seq_len;
        Span::new(start, start + len)
    } else {
        span_on_strand
    };
    let frame = (start_on_strand % 3) as u8 + 1;
    Orf { span, strand, frame, codons: len / 3 - 1 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_code_basics() {
        let code = GeneticCode::standard();
        assert_eq!(translate(b"ATGGCCTAA", code, false), b"MA*".to_vec());
        assert_eq!(translate(b"ATGGCCTAA", code, true), b"MA".to_vec());
        assert_eq!(code.translate_codon(b"GCN"), b'A');
        assert_eq!(code.translate_codon(b"NNN"), b'X');
        assert_eq!(code.translate_codon(b"TAR"), b'*');
        assert!(code.is_stop(b"TGA"));
        assert!(!code.is_stop(b"TGG"));
        assert!(code.is_start(b"ATG"));
        // NCBI table 11 lets GTG start; the standard code lists TTG, CTG and ATG.
        assert!(GeneticCode::get(11).unwrap().is_start(b"GTG"));
        assert_eq!(GeneticCode::get(2).unwrap().translate_codon(b"TGA"), b'W');
    }

    #[test]
    fn orfs_on_both_strands() {
        // Forward ORF: ATG + 4 codons + stop. Reverse ORF placed as reverse complement.
        let fwd = b"ATGAAACCCGGGTTTTAA";
        let rev_orf = crate::alphabet::revcomp(b"ATGCCCCCCCCCTGA");
        let mut seq = b"CC".to_vec();
        seq.extend_from_slice(fwd);
        seq.extend_from_slice(b"GG");
        seq.extend_from_slice(&rev_orf);
        let options = OrfOptions { min_codons: 3, ..OrfOptions::default() };
        let orfs = find_orfs(&seq, false, &options);
        assert_eq!(orfs.len(), 2, "{orfs:?}");
        assert_eq!(orfs[0].span, Span::new(2, 20));
        assert_eq!(orfs[0].strand, Strand::Forward);
        assert_eq!(orfs[0].codons, 5);
        assert_eq!(orfs[1].strand, Strand::Reverse);
        assert_eq!(orfs[1].span, Span::new(22, 37));
    }

    #[test]
    fn orf_through_origin_of_circular_sequence() {
        // Rotate an ORF so it straddles the origin.
        let orf = b"ATGAAAAAAAAAAAATAG";
        let mut seq = orf[9..].to_vec();
        seq.extend_from_slice(b"CCCCCC");
        seq.extend_from_slice(&orf[..9]);
        let options = OrfOptions { min_codons: 3, both_strands: false, ..OrfOptions::default() };
        let found = find_orfs(&seq, true, &options);
        assert_eq!(found.len(), 1, "{found:?}");
        let span = found[0].span;
        assert_eq!(span.extract(&seq), orf.to_vec());
        assert!(span.wraps(seq.len()));
        assert!(find_orfs(&seq, false, &options).is_empty());
    }
}
