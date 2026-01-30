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
) -> usize {
    let header_lines = 1; // Header line (ID + description)

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
        let lines = calculate_record_lines(&record, 80, false, false);
        assert_eq!(lines, 3); // 1 header + 2 sequence lines

        // Test with no wrapping - always 2 lines (header + sequence)
        let lines = calculate_record_lines(&record, 80, true, false);
        assert_eq!(lines, 2);

        // Test with quality shown and wrapping - doubles sequence lines
        let lines = calculate_record_lines(&record, 80, false, true);
        assert_eq!(lines, 5); // 1 header + 2 sequence + 2 quality
    }
}
