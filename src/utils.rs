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
/// the sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadOrientation {
    AsIs,
    Reverse,
    ReverseComplement,
}

impl ReadOrientation {
    pub fn next(&self) -> Self {
        match self {
            Self::AsIs => Self::Reverse,
            Self::Reverse => Self::ReverseComplement,
            Self::ReverseComplement => Self::AsIs,
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            Self::AsIs => "5'->3'",
            Self::Reverse => "3'->5'",
            Self::ReverseComplement => "RC",
        }
    }
}

// Standard genetic code. Index: (first*16) + (second*4) + third, where
// A=0, C=1, G=2, T=3. '*' marks a stop codon.
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

fn translate_codon(first: u8, second: u8, third: u8) -> char {
    let first = base_index(first);
    let second = base_index(second);
    let third = base_index(third);
    if first == 4 || second == 4 || third == 4 {
        'X'
    } else {
        CODON_TABLE[(first << 4) | (second << 2) | third]
    }
}

/// Translate one reading frame of a nucleotide sequence into amino acids.
/// Only complete codons are translated; a trailing partial codon is ignored.
fn translate_frame(strand: &[u8], frame: usize) -> String {
    let mut aa = String::new();
    let mut i = frame;
    while i + 2 < strand.len() {
        aa.push(translate_codon(strand[i], strand[i + 1], strand[i + 2]));
        i += 3;
    }
    aa
}

/// Translate a nucleotide sequence in all six reading frames. Frames 1-3
/// are read from the given strand, frames 4-6 from its reverse complement.
pub fn translate_frames(seq: &[u8]) -> Vec<String> {
    let revcomp = dna::revcomp(seq);
    (0..3)
        .map(|frame| translate_frame(seq, frame))
        .chain((0..3).map(|frame| translate_frame(&revcomp, frame)))
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
        assert_eq!(ReadOrientation::AsIs.next(), ReadOrientation::Reverse);
        assert_eq!(
            ReadOrientation::Reverse.next(),
            ReadOrientation::ReverseComplement
        );
        assert_eq!(
            ReadOrientation::ReverseComplement.next(),
            ReadOrientation::AsIs
        );
        assert_eq!(ReadOrientation::AsIs.name(), "5'->3'");
        assert_eq!(ReadOrientation::Reverse.name(), "3'->5'");
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
}
