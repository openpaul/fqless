use bio::pattern_matching::bndm::BNDM;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct AdapterStats {
    pub total_adapters_found: u64,
    pub adapter_positions: Vec<u64>, // Histogram of positions where adapters are found (0-150 positions)
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

impl Default for AdapterStats {
    fn default() -> Self {
        Self {
            total_adapters_found: 0,
            adapter_positions: vec![0; 151],
            contaminated_reads: 0,
            adapters_detected: HashMap::new(),
        }
    }
}

pub struct AdapterDetector {
    adapters: Vec<KnownAdapter>,
}

#[derive(Debug, Clone)]
struct KnownAdapter {
    name: String,
    sequence: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct AdapterDetection {
    pub adapter_name: String,
    pub position: usize,
}

impl AdapterDetector {
    pub fn new() -> Self {
        Self {
            adapters: Self::default_adapters(),
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

    fn default_adapters() -> Vec<KnownAdapter> {
        macro_rules! add_adapter {
            ($vec:ident, $name:expr, $seq:expr) => {{
                let seq = $seq.to_vec();
                let rc = Self::revcomp($seq);
                $vec.push(KnownAdapter {
                    name: $name.to_string(),
                    sequence: seq,
                });
                $vec.push(KnownAdapter {
                    name: format!("{}_RC", $name),
                    sequence: rc,
                });
            }};
        }

        let mut adapters = Vec::new();

        // Illumina
        add_adapter!(adapters, "TruSeq_Universal", b"AGATCGGAAGAGCACAC");
        add_adapter!(adapters, "TruSeq_Index", b"AGATCGGAAGAGCGTCG");
        add_adapter!(adapters, "Nextera", b"TCGTCGGCAGCGTCAGA");
        add_adapter!(adapters, "PolyA", b"AAAAAAAAAAAAAAA");

        // Oxford Nanopore (ONT) - rapid adapters
        add_adapter!(adapters, "ONT_Rapid_Adapter_1", b"TTTCTGTTGGTGCTGATATTGC");
        add_adapter!(adapters, "ONT_Rapid_Adapter_2", b"ACTTGCCTGTCGCTCTATCTTC");

        // ONT ligation adapters (common prefix)
        add_adapter!(adapters, "ONT_Ligation_Adapter_1", b"AGATCGGAAGAGCACACGTC");
        add_adapter!(adapters, "ONT_Ligation_Adapter_2", b"AGATCGGAAGAGCGTCGTGT");

        // PacBio SMRTbell hairpin adapters (short prefix of forward)
        add_adapter!(adapters, "PacBio_SMRTbell_FWD", b"ATCTCTCTCTTTTCCTCCTCC");
        add_adapter!(adapters, "PacBio_SMRTbell_REV", b"ATCTCTCTCAACAACAACAA");

        // PacBio Iso-Seq universal primer
        add_adapter!(adapters, "PacBio_IsoSeq_Primer", b"AAGCAGTGGTATCAACGCAG");

        adapters
    }

    pub fn detect_adapters(&self, sequence: &[u8]) -> Vec<AdapterDetection> {
        let mut detections = Vec::new();

        for adapter in &self.adapters {
            if adapter.sequence.is_empty() || adapter.sequence.len() > sequence.len() {
                continue;
            }
            let matcher = BNDM::new(&adapter.sequence);
            for pos in matcher.find_all(sequence) {
                detections.push(AdapterDetection {
                    adapter_name: adapter.name.clone(),
                    position: pos,
                });
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
            if detection.position < stats.adapter_positions.len() {
                stats.adapter_positions[detection.position] += 1;
            }

            // Normalize adapter name by stripping "_RC" suffix if present
            let base_name = detection
                .adapter_name
                .strip_suffix("_RC")
                .unwrap_or(&detection.adapter_name);

            let entry = stats
                .adapters_detected
                .entry(base_name.to_string())
                .or_insert_with(|| AdapterMatch {
                    name: base_name.to_string(),
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
