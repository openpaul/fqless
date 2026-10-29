use aho_corasick::{AhoCorasick, AhoCorasickBuilder, MatchKind};
use std::collections::HashMap;

/// Maximum number of adapter families we can track per read in a `u64` bitmask.
const MAX_ADAPTER_FAMILIES: usize = 64;

#[derive(Debug, Clone, Default)]
pub struct AdapterStats {
    pub total_adapters_found: u64,
    pub contaminated_reads: u64,
    pub adapters_detected: HashMap<String, AdapterMatch>,
}

#[derive(Debug, Clone)]
pub struct AdapterMatch {
    pub name: String,
    pub count: u64,
    pub avg_position: f64,
    pub min_position: usize,
    pub max_position: usize,
}

pub struct AdapterDetector {
    matcher: AhoCorasick,
    /// Family id and display name, indexed by the matcher's pattern id.
    pattern_families: Vec<(u8, &'static str)>,
    family_count: usize,
}

#[derive(Debug, Clone)]
pub struct AdapterDetection {
    pub adapter_name: String,
    pub position: usize,
}

/// Adapter seeds used for detection.
///
/// These are deliberately shorter than the full adapter sequences: sequencing
/// reads often run off the end of an insert and contain only a truncated
/// adapter. Seeds are long enough (12-15 bp) to be specific in random DNA
/// while still catching a partially sequenced adapter at the 3' end. The
/// reverse complement of every seed is searched as well, so both read
/// orientations are covered.
///
/// Entries sharing a name are treated as one adapter family: a read is only
/// reported once per family, even if several seeds match or the adapter is
/// present more than once.
const ADAPTER_SEEDS: &[(&str, &[u8])] = &[
    // Illumina
    ("Illumina_TruSeq", b"AGATCGGAAGAGC"),
    ("Illumina_Nextera", b"TCGTCGGCAGCGT"),
    ("Illumina_Nextera", b"CTGTCTCTTATAC"),
    ("Illumina_Small_RNA", b"TGGAATTCTCGG"),
    // Oxford Nanopore
    ("ONT_Rapid", b"TTTCTGTTGGTGC"),
    ("ONT_Ligation", b"GGCGTCTGCTTGG"),
    // PacBio
    ("PacBio_SMRTbell", b"ATCTCTCTCTTTTC"),
    ("PacBio_IsoSeq", b"AAGCAGTGGTATC"),
    // Homopolymer tails (not strictly adapters, but useful contamination signal)
    ("PolyA_tail", b"AAAAAAAAAAAA"),
];

impl AdapterDetector {
    pub fn new() -> Self {
        let mut patterns = Vec::with_capacity(ADAPTER_SEEDS.len() * 2);
        let mut pattern_families = Vec::with_capacity(ADAPTER_SEEDS.len() * 2);
        let mut family_names: Vec<&'static str> = Vec::new();

        for &(name, seed) in ADAPTER_SEEDS {
            let family_id = match family_names.iter().position(|known| *known == name) {
                Some(id) => id,
                None => {
                    family_names.push(name);
                    family_names.len() - 1
                }
            };
            assert!(
                family_names.len() <= MAX_ADAPTER_FAMILIES,
                "too many adapter families for the per-read bitmask"
            );

            patterns.push(seed.to_vec());
            pattern_families.push((family_id as u8, name));

            let rc = Self::revcomp(seed);
            if rc.as_slice() != seed {
                patterns.push(rc);
                pattern_families.push((family_id as u8, name));
            }
        }

        let matcher = AhoCorasickBuilder::new()
            .ascii_case_insensitive(true)
            .match_kind(MatchKind::LeftmostFirst)
            .build(patterns)
            .expect("adapter seeds are valid patterns");

        Self {
            matcher,
            pattern_families,
            family_count: family_names.len(),
        }
    }

    fn rc_base(b: u8) -> u8 {
        match b {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            _ => b'N',
        }
    }

    fn revcomp(seq: &[u8]) -> Vec<u8> {
        let mut rc = Vec::with_capacity(seq.len());
        for &b in seq.iter().rev() {
            rc.push(Self::rc_base(b));
        }
        rc
    }

    /// Find adapter contamination in one read.
    ///
    /// Matching is case-insensitive and reports at most one hit per adapter
    /// family. The reported position is the start of the first matching seed.
    pub fn detect_adapters(&self, sequence: &[u8]) -> Vec<AdapterDetection> {
        if sequence.is_empty() {
            return Vec::new();
        }

        let mut seen_families = 0u64;
        let mut detections = Vec::new();

        // All seeds (forward and reverse complement) are matched in a single
        // case-insensitive pass. Leftmost-first keeps the results
        // non-overlapping.
        for found in self.matcher.find_iter(sequence) {
            let (family_id, name) = self.pattern_families[found.pattern().as_usize()];
            let family_bit = 1u64 << family_id;
            if seen_families & family_bit != 0 {
                continue;
            }

            seen_families |= family_bit;
            detections.push(AdapterDetection {
                adapter_name: name.to_string(),
                position: found.start(),
            });

            if detections.len() == self.family_count {
                break;
            }
        }

        detections
    }

    pub fn update_stats(&self, stats: &mut AdapterStats, detections: &[AdapterDetection]) {
        if detections.is_empty() {
            return;
        }

        stats.contaminated_reads += 1;
        stats.total_adapters_found += detections.len() as u64;

        for detection in detections {
            let entry = stats
                .adapters_detected
                .entry(detection.adapter_name.clone())
                .or_insert_with(|| AdapterMatch {
                    name: detection.adapter_name.clone(),
                    count: 0,
                    avg_position: 0.0,
                    min_position: usize::MAX,
                    max_position: 0,
                });

            let old_count = entry.count;
            entry.count += 1;

            entry.avg_position = ((entry.avg_position * old_count as f64)
                + detection.position as f64)
                / entry.count as f64;

            entry.min_position = entry.min_position.min(detection.position);
            entry.max_position = entry.max_position.max(detection.position);
        }
    }
}

impl Default for AdapterDetector {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let seq = b"ACGT";
        let rc = AdapterDetector::revcomp(seq);
        assert_eq!(rc, b"ACGT", "ACGT should be its own reverse complement");

        let seq2 = b"AAAA";
        let rc2 = AdapterDetector::revcomp(seq2);
        assert_eq!(rc2, b"TTTT", "AAAA reverse complement should be TTTT");

        let seq3 = b"GATTACA";
        let rc3 = AdapterDetector::revcomp(seq3);
        assert_eq!(
            rc3, b"TGTAATC",
            "GATTACA reverse complement should be TGTAATC"
        );
    }

    #[test]
    fn test_adapter_detection() {
        let detector = AdapterDetector::new();

        // Test with a known Illumina adapter sequence
        let sequence = b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAAGATCGGAAGAGCACAC";
        let detections = detector.detect_adapters(sequence);

        assert!(
            !detections.is_empty(),
            "Should detect Illumina_TruSeq adapter"
        );
        assert_eq!(detections[0].adapter_name, "Illumina_TruSeq");
        assert_eq!(detections[0].position, 32);
    }

    #[test]
    fn test_reverse_complement_adapter_detection() {
        let detector = AdapterDetector::new();
        let adapter = AdapterDetector::revcomp(b"AGATCGGAAGAGC");
        let mut sequence = b"GCTAGCTAGCTAGCTAGCTA".to_vec();
        sequence.extend_from_slice(&adapter);

        let detections = detector.detect_adapters(&sequence);
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].adapter_name, "Illumina_TruSeq");
        assert_eq!(detections[0].position, 20);
    }

    #[test]
    fn test_case_insensitive_detection() {
        let detector = AdapterDetector::new();
        let mut sequence = b"GCTAGCTAGCTA".to_vec();
        sequence.extend_from_slice(b"agatcggaagagc");

        let detections = detector.detect_adapters(&sequence);
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].adapter_name, "Illumina_TruSeq");
        assert_eq!(detections[0].position, 12);
    }

    #[test]
    fn test_truncated_adapter_at_read_end() {
        let detector = AdapterDetector::new();
        // Only the first 13 bases of the adapter were sequenced.
        let mut sequence = b"GCTAGCTAGCTAGCTAGCTAGCTA".to_vec();
        sequence.extend_from_slice(b"AGATCGGAAGAGC");

        let detections = detector.detect_adapters(&sequence);
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].adapter_name, "Illumina_TruSeq");
        assert_eq!(detections[0].position, 24);
    }

    #[test]
    fn test_shared_prefix_counts_once_per_family() {
        let detector = AdapterDetector::new();
        // TruSeq Universal followed by TruSeq Index: both share the seed,
        // but this is still a single adapter family in one read.
        let sequence = b"AGATCGGAAGAGCACACAGATCGGAAGAGCGTCGTGT";

        let detections = detector.detect_adapters(sequence);
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].adapter_name, "Illumina_TruSeq");
    }

    #[test]
    fn test_multiple_seeds_for_one_family_count_once() {
        let detector = AdapterDetector::new();
        // Two different Nextera seeds in the same read.
        let sequence = b"TCGTCGGCAGCGTACGTACGTACGTCTGTCTCTTATAC";

        let detections = detector.detect_adapters(sequence);
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].adapter_name, "Illumina_Nextera");
    }

    #[test]
    fn test_no_false_positive_on_random_sequence() {
        let detector = AdapterDetector::new();
        let sequence = b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";

        assert!(detector.detect_adapters(sequence).is_empty());
    }

    #[test]
    fn test_poly_t_tail_detected_as_poly_a_family() {
        let detector = AdapterDetector::new();
        let sequence = b"GCTAGCTAGCTAGCTTTTTTTTTTTTTTTT";

        let detections = detector.detect_adapters(sequence);
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].adapter_name, "PolyA_tail");
    }

    #[test]
    fn test_stats_update() {
        let detector = AdapterDetector::new();
        let mut stats = AdapterStats::default();

        let detections = vec![AdapterDetection {
            adapter_name: "Illumina_TruSeq".to_string(),
            position: 50,
        }];

        detector.update_stats(&mut stats, &detections);

        assert_eq!(stats.contaminated_reads, 1);
        assert_eq!(stats.total_adapters_found, 1);
        assert!(stats.adapters_detected.contains_key("Illumina_TruSeq"));

        let adapter_match = stats.adapters_detected.get("Illumina_TruSeq").unwrap();
        assert_eq!(adapter_match.count, 1);
        assert_eq!(adapter_match.avg_position, 50.0);
        assert_eq!(adapter_match.min_position, 50);
        assert_eq!(adapter_match.max_position, 50);
    }

    #[test]
    fn test_stats_average_position_and_counts() {
        let detector = AdapterDetector::new();
        let mut stats = AdapterStats::default();

        for position in [10, 20, 60] {
            detector.update_stats(
                &mut stats,
                &[AdapterDetection {
                    adapter_name: "Illumina_TruSeq".to_string(),
                    position,
                }],
            );
        }

        let adapter_match = stats.adapters_detected.get("Illumina_TruSeq").unwrap();
        assert_eq!(stats.contaminated_reads, 3);
        assert_eq!(stats.total_adapters_found, 3);
        assert_eq!(adapter_match.count, 3);
        assert_eq!(adapter_match.avg_position, 30.0);
        assert_eq!(adapter_match.min_position, 10);
        assert_eq!(adapter_match.max_position, 60);
    }
}
