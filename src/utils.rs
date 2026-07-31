use bio::alphabets::dna;
use bio::io::fastq;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhredRange {
    Solexa,
    Illumina1_3,
    Illumina1_5,
    Default,
}

impl PhredRange {
    pub fn range(&self) -> Option<(u8, u8)> {
        match self {
            Self::Default => Some((33, 126)),
            Self::Solexa => Some((59, 104)),
            Self::Illumina1_3 => Some((64, 126)),
            Self::Illumina1_5 => Some((64, 126)),
        }
    }

    pub fn all() -> &'static [Self] {
        &[
            Self::Solexa,
            Self::Illumina1_3,
            Self::Illumina1_5,
            Self::Default,
        ]
    }

    pub fn from_min_max(min: u8, max: u8) -> Self {
        for variant in Self::all() {
            if let Some((rmin, rmax)) = variant.range() {
                if min >= rmin && max <= rmax {
                    return *variant;
                }
            }
        }
        Self::Default
    }

    pub fn base_phred(&self) -> u8 {
        self.range().map(|(min, _)| min).unwrap_or(0)
    }

    pub fn top_phred(&self) -> u8 {
        self.range().map(|(_, max)| max).unwrap_or(0)
    }

    pub fn name(&self) -> &'static str {
        match self {
            Self::Solexa => "Solexa",
            Self::Illumina1_3 => "Illumina 1.3",
            Self::Illumina1_5 => "Illumina 1.8",
            Self::Default => "Default",
        }
    }
}

/// Orientation in which reads are displayed. Quality scores always stay
/// paired with their base, so the quality string is reversed together with
/// the sequence. FASTQ stores every read in 5'->3' order, so only the
/// reverse complement (an unoriented read's other strand) is meaningful.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadOrientation {
    AsIs,
    ReverseComplement,
}

impl ReadOrientation {
    pub fn next(&self) -> Self {
        match self {
            Self::AsIs => Self::ReverseComplement,
            Self::ReverseComplement => Self::AsIs,
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            Self::AsIs => "5'->3'",
            Self::ReverseComplement => "RC",
        }
    }
}

// Standard genetic code (NCBI transl_table=1), cross-checked against
// Biopython's Bio/Data/CodonTable.py (Table 1). Index:
// (first*16) + (second*4) + third, where A=0, C=1, G=2, T=3.
// '*' marks a stop codon.
const CODON_TABLE: [char; 64] = [
    'K', 'N', 'K', 'N', // AAA AAC AAG AAT
    'T', 'T', 'T', 'T', // ACA ACC ACG ACT
    'R', 'S', 'R', 'S', // AGA AGC AGG AGT
    'I', 'I', 'M', 'I', // ATA ATC ATG ATT
    'Q', 'H', 'Q', 'H', // CAA CAC CAG CAT
    'P', 'P', 'P', 'P', // CCA CCC CCG CCT
    'R', 'R', 'R', 'R', // CGA CGC CGG CGT
    'L', 'L', 'L', 'L', // CTA CTC CTG CTT
    'E', 'D', 'E', 'D', // GAA GAC GAG GAT
    'A', 'A', 'A', 'A', // GCA GCC GCG GCT
    'G', 'G', 'G', 'G', // GGA GGC GGG GGT
    'V', 'V', 'V', 'V', // GTA GTC GTG GTT
    '*', 'Y', '*', 'Y', // TAA TAC TAG TAT
    'S', 'S', 'S', 'S', // TCA TCC TCG TCT
    '*', 'C', 'W', 'C', // TGA TGC TGG TGT
    'L', 'F', 'L', 'F', // TTA TTC TTG TTT
];

fn base_index(base: u8) -> usize {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 4,
    }
}

fn codon_index(first: u8, second: u8, third: u8) -> Option<usize> {
    let first = base_index(first);
    let second = base_index(second);
    let third = base_index(third);
    if first == 4 || second == 4 || third == 4 {
        None
    } else {
        Some((first << 4) | (second << 2) | third)
    }
}

/// Per-codon position weights, packed two bits per position
/// (w0<<4 | w1<<2 | w2). Each weight is the number of the three possible
/// base substitutions at that position (out of 3) that change the amino
/// acid. A weight of 0 marks a position that never affects the amino acid
/// (e.g. the third base of four-fold degenerate codons). Indexed exactly
/// like CODON_TABLE.
const CODON_WEIGHTS: [u8; 64] = [
    0x3e, 0x3e, 0x3e, 0x3e, 0x3c, 0x3c, 0x3c, 0x3c, 0x2e, 0x3e, 0x2e, 0x3e, 0x3d, 0x3d, 0x3f, 0x3d,
    0x3e, 0x3e, 0x3e, 0x3e, 0x3c, 0x3c, 0x3c, 0x3c, 0x2c, 0x3c, 0x2c, 0x3c, 0x2c, 0x3c, 0x2c, 0x3c,
    0x3e, 0x3e, 0x3e, 0x3e, 0x3c, 0x3c, 0x3c, 0x3c, 0x3c, 0x3c, 0x3c, 0x3c, 0x3c, 0x3c, 0x3c, 0x3c,
    0x3a, 0x3e, 0x3e, 0x3e, 0x3c, 0x3c, 0x3c, 0x3c, 0x3b, 0x3e, 0x3f, 0x3e, 0x2e, 0x3e, 0x2e, 0x3e,
];

/// Phred boost applied to a position based on how many of the three
/// possible substitutions change the amino acid: round(10*log10(3/w)).
/// A position that only rarely changes the amino acid contributes a higher
/// effective quality. w=0 is never looked up.
const WEIGHT_BOOST: [u8; 4] = [0, 5, 2, 0];

/// Translate one codon to an amino acid and an effective quality score
/// (clamped to 0..=40). The score is the minimum over the positions that
/// can change the amino acid of the position's phred score boosted by its
/// weight, so a low-quality base that determines the amino acid flags the
/// whole amino acid as uncertain while a redundant third base is ignored.
/// Ambiguous codons translate to 'X' at quality 0; stop codons always
/// return quality 0 so they render in red.
fn codon_aa_quality(
    first: u8,
    second: u8,
    third: u8,
    qual_first: u8,
    qual_second: u8,
    qual_third: u8,
    base_phred: u8,
) -> (char, u8) {
    let Some(idx) = codon_index(first, second, third) else {
        return ('X', 0);
    };
    let aa = CODON_TABLE[idx];
    if aa == '*' {
        return (aa, 0);
    }
    let weights = CODON_WEIGHTS[idx];
    let quals = [
        qual_first.saturating_sub(base_phred).min(40),
        qual_second.saturating_sub(base_phred).min(40),
        qual_third.saturating_sub(base_phred).min(40),
    ];
    let mut best = u8::MAX;
    for (weight, qual) in [weights >> 4, (weights >> 2) & 0b11, weights & 0b11]
        .into_iter()
        .zip(quals)
    {
        if weight != 0 {
            best = best.min(qual + WEIGHT_BOOST[weight as usize]).min(40);
        }
    }
    (aa, best)
}

/// Translate one reading frame of a nucleotide sequence into amino acids,
/// paired with an effective quality score per amino acid. Only complete
/// codons are translated; a trailing partial codon is ignored. If the
/// quality string is shorter than the sequence, missing scores are treated
/// as quality 0.
fn translate_frame_qual(
    strand: &[u8],
    qual: &[u8],
    frame: usize,
    base_phred: u8,
) -> (String, Vec<u8>) {
    let mut aa = String::new();
    let mut scores = Vec::new();
    let mut i = frame;
    while i + 2 < strand.len() {
        let (ch, q) = codon_aa_quality(
            strand[i],
            strand[i + 1],
            strand[i + 2],
            qual.get(i).copied().unwrap_or(0),
            qual.get(i + 1).copied().unwrap_or(0),
            qual.get(i + 2).copied().unwrap_or(0),
            base_phred,
        );
        aa.push(ch);
        scores.push(q);
        i += 3;
    }
    (aa, scores)
}

/// Translate a nucleotide sequence in all six reading frames, pairing each
/// amino acid with an effective quality score derived from the phred scores
/// of the codon positions that determine the amino acid. Frames 1-3 are
/// read from the given strand, frames 4-6 from its reverse complement; the
/// quality string is reversed along with the sequence so scores stay paired
/// with their bases.
pub fn translate_frames_with_quality(
    seq: &[u8],
    qual: &[u8],
    base_phred: u8,
) -> Vec<(String, Vec<u8>)> {
    let revcomp = dna::revcomp(seq);
    let qual_rev: Vec<u8> = qual.iter().rev().copied().collect();
    (0..3)
        .map(|frame| translate_frame_qual(seq, qual, frame, base_phred))
        .chain((0..3).map(|frame| translate_frame_qual(&revcomp, &qual_rev, frame, base_phred)))
        .collect()
}

/// Translate a nucleotide sequence in all six reading frames. Frames 1-3
/// are read from the given strand, frames 4-6 from its reverse complement.
pub fn translate_frames(seq: &[u8]) -> Vec<String> {
    translate_frames_with_quality(seq, &[], 0)
        .into_iter()
        .map(|(frame, _)| frame)
        .collect()
}

pub fn determine_min_max_phred(records: &[fastq::Record]) -> (u8, u8) {
    let mut min_phred = u8::MAX;
    let mut max_phred = u8::MIN;

    for record in records {
        for &q in record.qual() {
            min_phred = min_phred.min(q);
            max_phred = max_phred.max(q);
        }
    }

    (min_phred, max_phred)
}

/// Calculate how many lines a FASTQ record takes on screen
pub fn calculate_record_lines(
    record: &fastq::Record,
    terminal_width: usize,
    no_wrap: bool,
    show_quality: bool,
    translate: bool,
) -> usize {
    let header_lines = 1; // Header line (ID + description)

    if translate {
        // One line for the header plus one line per reading frame
        let frames = translate_frames(record.seq());
        return header_lines
            + frames
                .iter()
                .map(|frame| {
                    if no_wrap {
                        1
                    } else {
                        frame.len().max(1).div_ceil(terminal_width)
                    }
                })
                .sum::<usize>();
    }

    let sequence_lines = if no_wrap {
        1 // No wrapping, always 1 line
    } else {
        // Calculate wrapped lines: ceil(sequence_length / terminal_width)
        let seq_len = record.seq().len();
        if seq_len == 0 {
            1
        } else {
            seq_len.div_ceil(terminal_width) // Ceiling division
        }
    };

    let quality_lines = if show_quality {
        // same as sequence_lines
        sequence_lines
    } else {
        0
    };

    header_lines + sequence_lines + quality_lines
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phred_range_detection() {
        // Test Default/Sanger range (33-126)
        let range = PhredRange::from_min_max(33, 75);
        assert_eq!(range, PhredRange::Default);

        // Test Solexa range (59-104)
        let range = PhredRange::from_min_max(59, 100);
        assert_eq!(range, PhredRange::Solexa);

        // Test Illumina 1.3+ range (64-126) - needs higher max to distinguish from Solexa
        let range = PhredRange::from_min_max(64, 110);
        assert_eq!(range, PhredRange::Illumina1_3);
    }

    #[test]
    fn test_phred_range_methods() {
        let default_range = PhredRange::Default;
        assert_eq!(default_range.base_phred(), 33);
        assert_eq!(default_range.top_phred(), 126);
        assert_eq!(default_range.name(), "Default");
    }

    #[test]
    fn test_determine_min_max_phred() {
        use bio::io::fastq::Record;

        // Create a record with known quality scores
        let record = Record::with_attrs("test", None, b"ACGT", b"!IJK"); // ASCII 33, 73, 74, 75

        let records = vec![record];
        let (min, max) = determine_min_max_phred(&records);

        assert_eq!(min, 33);
        assert_eq!(max, 75);
    }

    #[test]
    fn test_calculate_record_lines() {
        use bio::io::fastq::Record;

        // Create a test record with 100bp sequence
        let record = Record::with_attrs("test", None, b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");

        // Test with wrapping, no quality - 100bp in 80 width terminal = 2 lines + 1 header = 3 lines
        let lines = calculate_record_lines(&record, 80, false, false, false);
        assert_eq!(lines, 3); // 1 header + 2 sequence lines

        // Test with no wrapping - always 2 lines (header + sequence)
        let lines = calculate_record_lines(&record, 80, true, false, false);
        assert_eq!(lines, 2);

        // Test with quality shown and wrapping - doubles sequence lines
        let lines = calculate_record_lines(&record, 80, false, true, false);
        assert_eq!(lines, 5); // 1 header + 2 sequence + 2 quality
    }

    #[test]
    fn test_calculate_record_lines_translation() {
        use bio::io::fastq::Record;

        // 9 bases -> 6 frames, each 1 line (fits in width) -> 7 lines total
        let record = Record::with_attrs("test", None, b"ATGAAATAA", b"IIIIIIIII");

        let lines = calculate_record_lines(&record, 80, false, false, true);
        assert_eq!(lines, 7); // 1 header + 6 frames

        let lines = calculate_record_lines(&record, 80, true, false, true);
        assert_eq!(lines, 7);
    }

    #[test]
    fn test_read_orientation_cycle() {
        assert_eq!(
            ReadOrientation::AsIs.next(),
            ReadOrientation::ReverseComplement
        );
        assert_eq!(
            ReadOrientation::ReverseComplement.next(),
            ReadOrientation::AsIs
        );
        assert_eq!(ReadOrientation::AsIs.name(), "5'->3'");
        assert_eq!(ReadOrientation::ReverseComplement.name(), "RC");
    }

    #[test]
    fn test_translate_frames_simple() {
        // ATG AAA TAA -> M K *
        let frames = translate_frames(b"ATGAAATAA");
        assert_eq!(frames.len(), 6);
        assert_eq!(frames[0], "MK*");

        // Frame 2 starts at index 1: TGA AAC -> * N
        assert_eq!(frames[1], "*N");
        // Frame 3 starts at index 2: GAA ATA -> E I
        assert_eq!(frames[2], "EI");
    }

    #[test]
    fn test_translate_frames_reverse_complement() {
        // Reverse complement of AAAACCCC is GGGGTTTT, so frame 4 reads
        // GGG GTT -> G V
        let frames = translate_frames(b"AAAACCCC");
        assert_eq!(frames[3], "GV");
        // Frame 6 starts at index 2 of revcomp: GGT TTT -> G F
        assert_eq!(frames[5], "GF");
    }

    #[test]
    fn test_translate_frames_ambiguous_and_stops() {
        // N codons translate to X, stop codons to '*'
        assert_eq!(translate_frames(b"TGA")[0], "*");
        assert_eq!(translate_frames(b"NNN")[0], "X");
    }

    #[test]
    fn test_codon_table_matches_standard_genetic_code() {
        // Full-table regression check against the NCBI standard genetic code
        // (transl_table=1), cross-checked with Biopython CodonTable.py Table 1.
        // Listed in A/C/G/T order to match CODON_TABLE's index layout.
        assert_eq!(
            CODON_TABLE.iter().collect::<String>(),
            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
        );
    }

    #[test]
    fn test_codon_weights_consistent_with_table() {
        // For every codon, recompute how many of the three possible base
        // substitutions at each position change the amino acid and check it
        // matches the packed CODON_WEIGHTS constant.
        let bases = b"ACGT";
        let mut idx = 0;
        for &first in bases {
            for &second in bases {
                for &third in bases {
                    let ref_aa = CODON_TABLE
                        [(base_index(first) << 4) | (base_index(second) << 2) | base_index(third)];
                    let mut weight = [0u8; 3];
                    for (pos, base) in [first, second, third].iter().enumerate() {
                        for &alt in bases {
                            if alt == *base {
                                continue;
                            }
                            let mut codon = [first, second, third];
                            codon[pos] = alt;
                            let aa = CODON_TABLE[(base_index(codon[0]) << 4)
                                | (base_index(codon[1]) << 2)
                                | base_index(codon[2])];
                            if aa != ref_aa {
                                weight[pos] += 1;
                            }
                        }
                    }
                    let packed = (weight[0] << 4) | (weight[1] << 2) | weight[2];
                    assert_eq!(CODON_WEIGHTS[idx], packed, "codon index {idx}");
                    idx += 1;
                }
            }
        }
    }

    #[test]
    fn test_codon_quality_fourfold_excludes_third_base() {
        // GGG is Gly, a four-fold degenerate codon: positions 1-2 determine
        // the amino acid, so a bad third base must not drag the score down.
        let (aa, q) = codon_aa_quality(b'G', b'G', b'G', 73, 73, 33, 33);
        assert_eq!(aa, 'G');
        assert_eq!(q, 40);
    }

    #[test]
    fn test_codon_quality_twofolds_keep_third_base() {
        // GAA is Glu; the third base still decides Glu vs Asp, so it counts
        // (weight 2 -> +2 boost), turning a score-0 base into effective 2.
        let (aa, q) = codon_aa_quality(b'G', b'A', b'A', 73, 73, 33, 33);
        assert_eq!(aa, 'E');
        assert_eq!(q, 2);
    }

    #[test]
    fn test_codon_quality_threefold_keeps_third_base() {
        // ATA is Ile; only one of three substitutions at the third base
        // (to ATG/Met) changes the amino acid, so it is kept with a +5 boost.
        let (aa, q) = codon_aa_quality(b'A', b'T', b'A', 73, 73, 33, 33);
        assert_eq!(aa, 'I');
        assert_eq!(q, 5);
    }

    #[test]
    fn test_codon_quality_min_of_informative_positions() {
        // ATG is Met, unique codon: all three positions matter, so the score
        // is the worst of the three.
        let (aa, q) = codon_aa_quality(b'A', b'T', b'G', 73, 43, 73, 33);
        assert_eq!(aa, 'M');
        assert_eq!(q, 10);
    }

    #[test]
    fn test_codon_quality_boost_clamped() {
        // All-high-quality GAA: the boosted third base clamps back to 40.
        let (aa, q) = codon_aa_quality(b'G', b'A', b'A', 73, 73, 73, 33);
        assert_eq!(aa, 'E');
        assert_eq!(q, 40);
    }

    #[test]
    fn test_codon_quality_ambiguous_and_stops() {
        // N codons translate to X at quality 0, stop codons stay at quality 0.
        assert_eq!(codon_aa_quality(b'N', b'N', b'N', 73, 73, 73, 33), ('X', 0));
        assert_eq!(codon_aa_quality(b'T', b'G', b'A', 73, 73, 73, 33), ('*', 0));
    }

    #[test]
    fn test_translate_frames_with_quality() {
        // ATG GGG TAA -> M G *, with the third base of the four-fold GGG
        // codon at quality 0 (ASCII 33).
        let seq = b"ATGGGGTAA";
        let qual = [73u8, 73, 73, 73, 73, 33, 73, 73, 73];
        let frames = translate_frames_with_quality(seq, &qual, 33);
        assert_eq!(frames[0].0, "MG*");
        assert_eq!(frames[0].1, vec![40, 40, 0]);
    }

    #[test]
    fn test_translate_known_protein_hbb() {
        // Human beta-globin CDS from NCBI NM_000518.5 (positions 51..494,
        // 444 nt, ends in the TAA stop codon) and its translation
        // NP_000509.1 (hemoglobin subunit beta, 147 aa).
        let cds = concat!(
            "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAG",
            "TTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGG",
            "GGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGT",
            "GCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACT",
            "GTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCA",
            "TCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAAT",
            "GCCCTGGCCCACAAGTATCACTAA",
        );
        let protein = concat!(
            "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLG",
            "AFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVAN",
            "ALAHKYH*",
        );

        let frames = translate_frames(cds.as_bytes());
        assert_eq!(frames.len(), 6);
        assert_eq!(frames[0], protein);
        assert_eq!(frames[0].len(), 148);
    }
}
