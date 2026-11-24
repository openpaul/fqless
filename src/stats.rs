//! Read statistics shown on the stats page.
//!
//! [`StatsAccumulator`] is fed one record at a time by a background scan (or,
//! for stdin, by the main thread as reads are loaded) and turned into a
//! [`FastqStats`] snapshot for drawing.

use crate::adapter::{AdapterDetector, AdapterStats};
use crate::utils::PhredRange;
use bio::io::fastq;

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
pub struct FastqStats {
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
}

impl Default for FastqStats {
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
        }
    }
}

/// Running totals for [`FastqStats`].
pub struct StatsAccumulator {
    base_phred: u8,
    reads: u64,
    total_length: u64,
    min_length: usize,
    max_length: usize,
    total_gc: u64,
    total_n: u64,
    quality_histogram: Vec<u64>,
    position_quality_sums: Vec<u64>,
    position_counts: Vec<u64>,
    average_read_quality_histogram: Vec<u64>,
    adapter_detector: AdapterDetector,
    adapter_stats: AdapterStats,
}

impl StatsAccumulator {
    pub fn new(phred_range: PhredRange) -> Self {
        let bins = (phred_range.top_phred() - phred_range.base_phred()) as usize + 1;
        Self {
            base_phred: phred_range.base_phred(),
            reads: 0,
            total_length: 0,
            min_length: usize::MAX,
            max_length: 0,
            total_gc: 0,
            total_n: 0,
            quality_histogram: vec![0; bins],
            position_quality_sums: Vec::new(),
            position_counts: Vec::new(),
            average_read_quality_histogram: vec![0; AVG_QUAL_NUM_BINS],
            adapter_detector: AdapterDetector::new(),
            adapter_stats: AdapterStats::default(),
        }
    }

    /// Number of reads added so far.
    pub fn reads(&self) -> u64 {
        self.reads
    }

    pub fn add(&mut self, record: &fastq::Record) {
        let seq = record.seq();
        let qual = record.qual();

        self.reads += 1;
        self.total_length += seq.len() as u64;
        self.min_length = self.min_length.min(seq.len());
        self.max_length = self.max_length.max(seq.len());

        for &base in seq {
            match base {
                b'G' | b'g' | b'C' | b'c' => self.total_gc += 1,
                b'N' | b'n' => self.total_n += 1,
                _ => {}
            }
        }

        let detections = self.adapter_detector.detect_adapters(seq);
        self.adapter_detector
            .update_stats(&mut self.adapter_stats, &detections);

        if qual.len() > self.position_quality_sums.len() {
            self.position_quality_sums.resize(qual.len(), 0);
            self.position_counts.resize(qual.len(), 0);
        }

        let mut total_qual = 0u64;
        for (pos, &q) in qual.iter().enumerate() {
            let score = q.saturating_sub(self.base_phred);
            total_qual += score as u64;
            if let Some(bin) = self.quality_histogram.get_mut(score as usize) {
                *bin += 1;
            }
            self.position_quality_sums[pos] += score as u64;
            self.position_counts[pos] += 1;
        }

        // A record without qualities has no average to bin.
        if !qual.is_empty() {
            let average = total_qual as f64 / qual.len() as f64;
            self.average_read_quality_histogram[average_quality_bin(average)] += 1;
        }
    }

    /// Write the current totals into `stats`, reusing its buffers.
    pub fn snapshot_into(&self, stats: &mut FastqStats, scanned_all: bool) {
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
        stats.quality_histogram.clone_from(&self.quality_histogram);
        stats
            .average_read_quality_histogram
            .clone_from(&self.average_read_quality_histogram);
        stats.adapter_stats.clone_from(&self.adapter_stats);
        stats.position_quality.clear();
        stats.position_quality.extend(
            self.position_quality_sums
                .iter()
                .zip(&self.position_counts)
                .map(|(&sum, &count)| sum as f64 / count.max(1) as f64),
        );
        stats.scanned_all = scanned_all;
    }

    #[cfg(test)]
    pub fn snapshot(&self, scanned_all: bool) -> FastqStats {
        let mut stats = FastqStats::default();
        self.snapshot_into(&mut stats, scanned_all);
        stats
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::Record;

    fn record(seq: &[u8], qual: &[u8]) -> Record {
        Record::with_attrs("r", None, seq, qual)
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

    #[test]
    fn test_average_quality_bin_caps_at_q40() {
        assert_eq!(average_quality_bin(0.0), 0);
        assert_eq!(average_quality_bin(4.9), 0);
        assert_eq!(average_quality_bin(5.0), 1);
        assert_eq!(average_quality_bin(41.0), AVG_QUAL_NUM_BINS - 1);
        assert_eq!(average_quality_bin(-3.0), 0);
    }
}
