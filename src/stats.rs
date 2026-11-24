//! Read statistics shown on the stats page.
//!
//! [`StatsAccumulator`] is fed one record at a time by a background scan (or,
//! for stdin, by the main thread as reads are loaded) and turned into a
//! [`ReadStats`] snapshot for drawing.

use crate::adapter::{AdapterDetector, AdapterStats};
use crate::record::Record;
use crate::utils::PhredRange;

// Bins for the average-quality-per-read histogram: width 5, capped at Q40.
// The accumulator increments these bins directly instead of keeping one
// value per read, which would grow unboundedly for large files.
pub const AVG_QUAL_BIN_WIDTH: usize = 5;
const AVG_QUAL_MAX: usize = 40;
pub const AVG_QUAL_NUM_BINS: usize = (AVG_QUAL_MAX / AVG_QUAL_BIN_WIDTH) + 1;

/// Bin index in the average-quality histogram for a per-read average quality.
fn average_quality_bin(avg_qual: f64) -> usize {
    let idx = (avg_qual.max(0.0) as usize) / AVG_QUAL_BIN_WIDTH;
    idx.min(AVG_QUAL_NUM_BINS - 1)
}

#[derive(Debug, Clone)]
pub struct ReadStats {
    pub processed_reads: u64,
    pub avg_length: f64,
    pub min_length: usize,
    pub max_length: usize,
    /// Histogram of individual base quality scores, indexed by score.
    pub quality_histogram: Vec<u64>,
    /// Average quality at each position, over the reads that reach it.
    pub position_quality: Vec<f64>,
    pub gc_content: f64,
    pub n_content: f64,
    /// True once every read of the input has been processed.
    pub scanned_all: bool,
    /// Histogram of per-read average quality scores, binned with
    /// [`AVG_QUAL_BIN_WIDTH`]. One entry per bin, so its size is
    /// independent of the number of reads.
    pub average_read_quality_histogram: Vec<u64>,
    pub adapter_stats: AdapterStats,
    /// False when no read had base qualities (FASTA), so the quality
    /// charts are hidden instead of showing misleading empty data.
    pub quality_available: bool,
}

impl Default for ReadStats {
    fn default() -> Self {
        Self {
            processed_reads: 0,
            avg_length: 0.0,
            min_length: 0,
            max_length: 0,
            quality_histogram: vec![0; 94], // Quality scores 0-93
            position_quality: Vec::new(),
            gc_content: 0.0,
            n_content: 0.0,
            scanned_all: false,
            average_read_quality_histogram: vec![0; AVG_QUAL_NUM_BINS],
            adapter_stats: AdapterStats::default(),
            quality_available: false,
        }
    }
}

/// Running totals for [`ReadStats`].
///
/// The per-base work is the hot loop of the stats scan, so it is kept free
/// of data dependencies: quality characters are counted raw in four
/// interleaved histograms (most bases share one character, and a single
/// histogram would make every increment wait for the previous one), and the
/// number of reads reaching each position comes from a read-length histogram
/// instead of a per-base counter. Scores are derived at snapshot time.
pub struct StatsAccumulator {
    phred_range: PhredRange,
    reads: u64,
    total_length: u64,
    min_length: usize,
    max_length: usize,
    total_gc: u64,
    total_n: u64,
    /// Counts of raw quality characters, split four ways (see above).
    quality_chars: [[u64; 256]; 4],
    /// Sum of quality scores at each position.
    position_quality_sums: Vec<u64>,
    /// Number of reads whose quality string has exactly this length.
    quality_lengths: Vec<u64>,
    average_read_quality_histogram: Vec<u64>,
    adapter_detector: AdapterDetector,
    adapter_stats: AdapterStats,
    /// True once any read carried qualities.
    saw_quality: bool,
}

impl StatsAccumulator {
    pub fn new(phred_range: PhredRange) -> Self {
        Self {
            phred_range,
            reads: 0,
            total_length: 0,
            min_length: usize::MAX,
            max_length: 0,
            total_gc: 0,
            total_n: 0,
            quality_chars: [[0; 256]; 4],
            position_quality_sums: Vec::new(),
            quality_lengths: Vec::new(),
            average_read_quality_histogram: vec![0; AVG_QUAL_NUM_BINS],
            adapter_detector: AdapterDetector::new(),
            adapter_stats: AdapterStats::default(),
            saw_quality: false,
        }
    }

    /// Number of reads added so far.
    pub fn reads(&self) -> u64 {
        self.reads
    }

    pub fn add(&mut self, record: &Record) {
        let seq = record.seq();
        let qual = record.qual();

        self.reads += 1;
        self.total_length += seq.len() as u64;
        self.min_length = self.min_length.min(seq.len());
        self.max_length = self.max_length.max(seq.len());

        // Counting with filter().count() lets the compiler vectorize; `| 0x20`
        // folds upper case onto lower case.
        self.total_gc += seq
            .iter()
            .filter(|&&b| matches!(b | 0x20, b'g' | b'c'))
            .count() as u64;
        self.total_n += seq.iter().filter(|&&b| b | 0x20 == b'n').count() as u64;

        let detections = self.adapter_detector.detect_adapters(seq);
        self.adapter_detector
            .update_stats(&mut self.adapter_stats, &detections);

        // FASTA records have no qualities; an empty read has none to count.
        let Some(qual) = qual else {
            return;
        };
        self.saw_quality = true;
        if qual.is_empty() {
            return;
        }

        let (chunks, remainder) = qual.as_chunks::<4>();
        for chunk in chunks {
            self.quality_chars[0][chunk[0] as usize] += 1;
            self.quality_chars[1][chunk[1] as usize] += 1;
            self.quality_chars[2][chunk[2] as usize] += 1;
            self.quality_chars[3][chunk[3] as usize] += 1;
        }
        for &q in remainder {
            self.quality_chars[0][q as usize] += 1;
        }

        // Scores below the encoding's offset count as Q0.
        let base = self.phred_range.base_phred();
        if qual.len() > self.position_quality_sums.len() {
            self.position_quality_sums.resize(qual.len(), 0);
        }
        let mut total_qual = 0u64;
        for (sum, &q) in self.position_quality_sums.iter_mut().zip(qual) {
            let score = q.max(base) - base;
            *sum += score as u64;
            total_qual += score as u64;
        }

        if qual.len() >= self.quality_lengths.len() {
            self.quality_lengths.resize(qual.len() + 1, 0);
        }
        self.quality_lengths[qual.len()] += 1;

        let average = total_qual as f64 / qual.len() as f64;
        self.average_read_quality_histogram[average_quality_bin(average)] += 1;
    }

    /// Write the current totals into `stats`, reusing its buffers.
    pub fn snapshot_into(&self, stats: &mut ReadStats, scanned_all: bool) {
        let percent = |count: u64| {
            if self.total_length > 0 {
                count as f64 / self.total_length as f64 * 100.0
            } else {
                0.0
            }
        };

        stats.processed_reads = self.reads;
        stats.avg_length = if self.reads > 0 {
            self.total_length as f64 / self.reads as f64
        } else {
            0.0
        };
        stats.min_length = if self.reads > 0 { self.min_length } else { 0 };
        stats.max_length = self.max_length;
        stats.gc_content = percent(self.total_gc);
        stats.n_content = percent(self.total_n);

        // Map raw characters to scores. Characters below the offset count as
        // Q0, characters above the encoding's top are not shown.
        let base = self.phred_range.base_phred() as usize;
        let top = self.phred_range.top_phred() as usize;
        stats.quality_histogram.clear();
        stats.quality_histogram.resize(top - base + 1, 0);
        for ch in 0..=top {
            let count: u64 = self.quality_chars.iter().map(|h| h[ch]).sum();
            stats.quality_histogram[ch.saturating_sub(base)] += count;
        }

        stats
            .average_read_quality_histogram
            .clone_from(&self.average_read_quality_histogram);
        stats.adapter_stats.clone_from(&self.adapter_stats);

        // Reads reaching position p are those with a quality string longer
        // than p: a running sum of the length histogram from the far end.
        stats.position_quality.clear();
        stats
            .position_quality
            .resize(self.position_quality_sums.len(), 0.0);
        let mut reaching = 0u64;
        for pos in (0..self.position_quality_sums.len()).rev() {
            reaching += self.quality_lengths.get(pos + 1).copied().unwrap_or(0);
            stats.position_quality[pos] =
                self.position_quality_sums[pos] as f64 / reaching.max(1) as f64;
        }
        stats.quality_available = self.saw_quality;
        stats.scanned_all = scanned_all;
    }

    #[cfg(test)]
    pub fn snapshot(&self, scanned_all: bool) -> ReadStats {
        let mut stats = ReadStats::default();
        self.snapshot_into(&mut stats, scanned_all);
        stats
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Record;

    fn record(seq: &[u8], qual: &[u8]) -> Record {
        Record::with_attrs("r", None, seq, Some(qual))
    }

    #[test]
    fn test_records_without_quality() {
        let mut acc = StatsAccumulator::new(PhredRange::Default);
        acc.add(&Record::with_attrs("r", None, b"GGAT", None));
        let stats = acc.snapshot(true);
        assert_eq!(stats.processed_reads, 1);
        assert!((stats.gc_content - 50.0).abs() < 1e-9);
        assert!(!stats.quality_available);
        assert!(stats.position_quality.is_empty());
        assert_eq!(stats.average_read_quality_histogram.iter().sum::<u64>(), 0);

        acc.add(&record(b"AC", b"II"));
        assert!(acc.snapshot(true).quality_available);
    }

    #[test]
    fn test_empty_input() {
        let acc = StatsAccumulator::new(PhredRange::Default);
        let stats = acc.snapshot(true);
        assert_eq!(stats.processed_reads, 0);
        assert_eq!(stats.min_length, 0);
        assert_eq!(stats.max_length, 0);
        assert_eq!(stats.avg_length, 0.0);
        assert_eq!(stats.gc_content, 0.0);
        assert!(stats.scanned_all);
    }

    #[test]
    fn test_lengths_and_composition() {
        let mut acc = StatsAccumulator::new(PhredRange::Default);
        acc.add(&record(b"GGCC", b"IIII"));
        acc.add(&record(b"ATNNAT", b"IIIIII"));
        let stats = acc.snapshot(false);

        assert_eq!(stats.processed_reads, 2);
        assert_eq!(stats.min_length, 4);
        assert_eq!(stats.max_length, 6);
        assert_eq!(stats.avg_length, 5.0);
        assert!((stats.gc_content - 40.0).abs() < 1e-9); // 4 of 10
        assert!((stats.n_content - 20.0).abs() < 1e-9); // 2 of 10
        assert!(!stats.scanned_all);
    }

    #[test]
    fn test_quality_histograms_use_phred_offset() {
        let mut acc = StatsAccumulator::new(PhredRange::Default);
        // '!' = Q0, '+' = Q10, 'I' = Q40
        acc.add(&record(b"ACG", b"!+I"));
        let stats = acc.snapshot(true);

        assert_eq!(stats.quality_histogram[0], 1);
        assert_eq!(stats.quality_histogram[10], 1);
        assert_eq!(stats.quality_histogram[40], 1);
        // Average (0 + 10 + 40) / 3 = 16.7 -> bin Q15
        assert_eq!(stats.average_read_quality_histogram[3], 1);
        assert_eq!(stats.position_quality, vec![0.0, 10.0, 40.0]);
    }

    #[test]
    fn test_position_quality_averages_only_reads_reaching_the_position() {
        let mut acc = StatsAccumulator::new(PhredRange::Default);
        acc.add(&record(b"AC", b"++")); // Q10 Q10
        acc.add(&record(b"ACGT", b"5555")); // Q20 x4
        let stats = acc.snapshot(true);
        assert_eq!(stats.position_quality, vec![15.0, 15.0, 20.0, 20.0]);
    }

    #[test]
    fn test_short_quality_string_does_not_panic() {
        let mut acc = StatsAccumulator::new(PhredRange::Default);
        acc.add(&record(b"ACGTACGT", b"III"));
        let stats = acc.snapshot(true);
        assert_eq!(stats.position_quality.len(), 3);
        assert_eq!(stats.average_read_quality_histogram[8], 1); // Q40
    }

    /// The accumulator's split histograms and length-based position counts
    /// must give the same numbers as the obvious per-base computation.
    #[test]
    fn test_matches_naive_reference() {
        // Small deterministic generator, so the test needs no extra crate.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = move |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };

        for phred_range in [PhredRange::Default, PhredRange::Illumina1_3] {
            let base = phred_range.base_phred();
            let bins = (phred_range.top_phred() - base) as usize + 1;
            let mut acc = StatsAccumulator::new(phred_range);
            let mut hist = vec![0u64; bins];
            let mut avg_hist = vec![0u64; AVG_QUAL_NUM_BINS];
            let mut pos_sum: Vec<u64> = Vec::new();
            let mut pos_count: Vec<u64> = Vec::new();

            for _ in 0..500 {
                let len = next(40) as usize; // includes 0 and lengths not divisible by 4
                let seq: Vec<u8> = (0..len).map(|_| b"ACGTNacgtn"[next(10) as usize]).collect();
                // Every quality character from `!` to `~`; for Phred+64 that
                // includes characters below the offset.
                let qual: Vec<u8> = (0..len).map(|_| b'!' + next(94) as u8).collect();
                acc.add(&record(&seq, &qual));

                if pos_sum.len() < len {
                    pos_sum.resize(len, 0);
                    pos_count.resize(len, 0);
                }
                let mut total = 0u64;
                for (i, &q) in qual.iter().enumerate() {
                    let score = q.saturating_sub(base) as u64;
                    if let Some(bin) = hist.get_mut(score as usize) {
                        *bin += 1;
                    }
                    pos_sum[i] += score;
                    pos_count[i] += 1;
                    total += score;
                }
                if len > 0 {
                    avg_hist[average_quality_bin(total as f64 / len as f64)] += 1;
                }
            }

            let stats = acc.snapshot(true);
            assert_eq!(stats.quality_histogram, hist, "{phred_range:?}");
            assert_eq!(stats.average_read_quality_histogram, avg_hist);
            let expected: Vec<f64> = pos_sum
                .iter()
                .zip(&pos_count)
                .map(|(&s, &c)| s as f64 / c.max(1) as f64)
                .collect();
            assert_eq!(stats.position_quality, expected);
        }
    }

    #[test]
    fn test_average_quality_bin_caps_at_q40() {
        assert_eq!(average_quality_bin(0.0), 0);
        assert_eq!(average_quality_bin(4.9), 0);
        assert_eq!(average_quality_bin(5.0), 1);
        assert_eq!(average_quality_bin(41.0), AVG_QUAL_NUM_BINS - 1);
        assert_eq!(average_quality_bin(-3.0), 0);
    }
}
