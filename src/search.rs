use bio::alphabets::dna;
use bio::pattern_matching::bom::BOM;

/// Bitmap of which records contain at least one search match. One bit per
/// record, so memory use is independent of how frequent the pattern is.
pub struct SearchIndex {
    bits: Vec<u64>,
    count: u64,
}

impl SearchIndex {
    pub fn new() -> Self {
        Self {
            bits: Vec::new(),
            count: 0,
        }
    }

    /// Clear all bits and the match count.
    pub fn reset(&mut self) {
        self.bits.clear();
        self.count = 0;
    }

    /// Mark a record as containing a match. Returns true if the bit was
    /// newly set (i.e. the total match count increased).
    pub fn set(&mut self, record: u64) -> bool {
        let word = (record / 64) as usize;
        let bit = record % 64;
        if word >= self.bits.len() {
            self.bits.resize(word + 1, 0);
        }
        let mask = 1u64 << bit;
        if self.bits[word] & mask == 0 {
            self.bits[word] |= mask;
            self.count += 1;
            true
        } else {
            false
        }
    }

    #[cfg(test)]
    pub fn test(&self, record: u64) -> bool {
        let word = (record / 64) as usize;
        let bit = record % 64;
        self.bits.get(word).is_some_and(|w| (w >> bit) & 1 == 1)
    }

    pub fn count(&self) -> u64 {
        self.count
    }

    /// Number of set bits with index <= `record`.
    pub fn count_up_to(&self, record: u64) -> u64 {
        let word = (record / 64) as usize;
        let offset = record % 64;
        let mut total = 0u64;
        let full = word.min(self.bits.len());
        for w in &self.bits[..full] {
            total += w.count_ones() as u64;
        }
        if word < self.bits.len() {
            let w = self.bits[word];
            let keep = if offset == 63 {
                !0u64
            } else {
                (1u64 << (offset + 1)) - 1
            };
            total += (w & keep).count_ones() as u64;
        }
        total
    }

    /// First set bit at or after `from`.
    pub fn next(&self, from: u64) -> Option<u64> {
        let start_word = (from / 64) as usize;
        let offset = from % 64;
        for (i, word) in self.bits.iter().enumerate().skip(start_word) {
            let mut w = *word;
            if i == start_word {
                w &= !0u64 << offset;
            }
            if w != 0 {
                return Some((i as u64) * 64 + w.trailing_zeros() as u64);
            }
        }
        None
    }

    /// Last set bit strictly below `below`.
    pub fn prev(&self, below: u64) -> Option<u64> {
        if below == 0 {
            return None;
        }
        let last = below - 1;
        let last_word = (last / 64) as usize;
        let offset = last % 64;
        for i in (0..self.bits.len().min(last_word + 1)).rev() {
            let mut w = self.bits[i];
            if i == last_word {
                w &= if offset == 63 {
                    !0u64
                } else {
                    (1u64 << (offset + 1)) - 1
                };
            }
            if w != 0 {
                return Some((i as u64) * 64 + (63 - w.leading_zeros()) as u64);
            }
        }
        None
    }
}

impl Default for SearchIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// A case-insensitive exact-substring pattern, optionally searched together
/// with its reverse complement. Reads are unoriented, so a motif may appear
/// as either the literal pattern or its reverse complement in a read.
#[derive(Clone)]
pub struct Pattern {
    literal: Vec<u8>,
    literal_bom: BOM,
    rc_bom: Option<BOM>,
}

impl Pattern {
    /// Build a pattern from user input (leading/trailing whitespace is
    /// trimmed). Returns None for empty input.
    pub fn new(input: &str, include_rc: bool) -> Option<Self> {
        let input = input.trim();
        if input.is_empty() {
            return None;
        }
        let literal: Vec<u8> = input
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_lowercase())
            .collect();
        let rc_valid = literal
            .iter()
            .all(|b| matches!(b, b'a' | b'c' | b'g' | b't' | b'u' | b'n'));
        let rc = if rc_valid {
            dna::revcomp(&literal)
        } else {
            literal.clone()
        };
        let literal_bom = BOM::new(&literal);
        let rc_bom = if rc_valid && include_rc {
            Some(BOM::new(&rc))
        } else {
            None
        };
        Some(Self {
            literal,
            literal_bom,
            rc_bom,
        })
    }

    /// True if any occurrence of the pattern (or its reverse complement)
    /// appears in `seq`.
    pub fn has_match(&self, seq: &[u8]) -> bool {
        let lc: Vec<u8> = seq.iter().map(|b| b.to_ascii_lowercase()).collect();
        if self.literal_bom.find_all(&lc).next().is_some() {
            return true;
        }
        self.rc_bom
            .as_ref()
            .is_some_and(|bom| bom.find_all(&lc).next().is_some())
    }

    /// Find all matches in `seq`, merging overlapping and adjacent hits into
    /// contiguous ranges. Ranges are [start, end) byte offsets into `seq`.
    pub fn find_merged_matches(&self, seq: &[u8]) -> Vec<(usize, usize)> {
        let lc: Vec<u8> = seq.iter().map(|b| b.to_ascii_lowercase()).collect();
        let mut starts: Vec<usize> = self.literal_bom.find_all(&lc).collect();
        if let Some(rc) = &self.rc_bom {
            starts.extend(rc.find_all(&lc));
        }
        starts.sort_unstable();

        let mut ranges: Vec<(usize, usize)> = Vec::new();
        for start in starts {
            let end = start + self.literal.len();
            if let Some((_, last_end)) = ranges.last_mut() {
                if start <= *last_end {
                    *last_end = (*last_end).max(end);
                    continue;
                }
            }
            ranges.push((start, end));
        }
        ranges
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index_with(matches: &[u64]) -> SearchIndex {
        let mut idx = SearchIndex::new();
        for &m in matches {
            idx.set(m);
        }
        idx
    }

    #[test]
    fn test_set_test_count() {
        let mut idx = SearchIndex::new();
        assert_eq!(idx.count(), 0);
        assert!(idx.set(0));
        assert!(!idx.set(0)); // duplicate set does not change count
        assert!(idx.set(64));
        assert!(idx.set(1000));
        assert_eq!(idx.count(), 3);
        assert!(idx.test(0));
        assert!(idx.test(64));
        assert!(idx.test(1000));
        assert!(!idx.test(63));
        assert!(!idx.test(65));
    }

    #[test]
    fn test_next() {
        let idx = index_with(&[5, 100, 101]);
        assert_eq!(idx.next(0), Some(5));
        assert_eq!(idx.next(5), Some(5));
        assert_eq!(idx.next(6), Some(100));
        assert_eq!(idx.next(101), Some(101));
        assert_eq!(idx.next(102), None);
    }

    #[test]
    fn test_next_across_words() {
        let idx = index_with(&[63, 64, 128]);
        assert_eq!(idx.next(0), Some(63));
        assert_eq!(idx.next(64), Some(64));
        assert_eq!(idx.next(65), Some(128));
        assert_eq!(idx.next(129), None);
    }

    #[test]
    fn test_prev() {
        let idx = index_with(&[5, 100, 101]);
        assert_eq!(idx.prev(102), Some(101));
        assert_eq!(idx.prev(101), Some(100));
        assert_eq!(idx.prev(6), Some(5));
        assert_eq!(idx.prev(5), None); // strictly below
        assert_eq!(idx.prev(0), None);
    }

    #[test]
    fn test_prev_empty() {
        let idx = SearchIndex::new();
        assert_eq!(idx.prev(1000), None);
        assert_eq!(idx.next(0), None);
    }

    #[test]
    fn test_count_up_to() {
        let idx = index_with(&[1, 64, 65, 1000]);
        assert_eq!(idx.count_up_to(0), 0);
        assert_eq!(idx.count_up_to(1), 1);
        assert_eq!(idx.count_up_to(63), 1);
        assert_eq!(idx.count_up_to(64), 2);
        assert_eq!(idx.count_up_to(65), 3);
        assert_eq!(idx.count_up_to(1000), 4);
        assert_eq!(idx.count_up_to(9999), 4);
    }

    #[test]
    fn test_reset() {
        let mut idx = index_with(&[1, 2, 3]);
        idx.reset();
        assert_eq!(idx.count(), 0);
        assert_eq!(idx.next(0), None);
    }

    #[test]
    fn test_pattern_case_insensitive() {
        let p = Pattern::new("atg", false).unwrap();
        assert!(p.has_match(b"CCATGCC"));
        assert!(p.has_match(b"ccatgcc"));
        assert!(p.has_match(b"ATG"));
        assert!(!p.has_match(b"GCA"));
    }

    #[test]
    fn test_pattern_reverse_complement() {
        let p = Pattern::new("atgc", true).unwrap();
        assert!(p.has_match(b"ATGC"));
        assert!(p.has_match(b"GCAT")); // rc(ATGC) = GCAT
        assert!(!p.has_match(b"TGCA"));
    }

    #[test]
    fn test_pattern_rc_disabled() {
        let p = Pattern::new("atgc", false).unwrap();
        assert!(p.has_match(b"ATGC"));
        assert!(!p.has_match(b"GCAT"));
    }

    #[test]
    fn test_pattern_invalid_bases_no_rc() {
        // Pattern with non-DNA chars still matches literally without panicking.
        let p = Pattern::new("acgtx", true).unwrap();
        assert!(p.has_match(b"ACGTX"));
        assert!(!p.has_match(b"XTGCA"));
    }

    #[test]
    fn test_pattern_empty() {
        assert!(Pattern::new("", true).is_none());
        assert!(Pattern::new("   ", true).is_none());
    }

    #[test]
    fn test_merged_matches_overlap() {
        // "AAA" occurs overlapping in "AAAA"; hits merge into one span.
        let p = Pattern::new("aaa", false).unwrap();
        assert_eq!(p.find_merged_matches(b"AAAA"), vec![(0, 4)]);
        assert_eq!(p.find_merged_matches(b"CCAAAC"), vec![(2, 5)]);
    }

    #[test]
    fn test_merged_matches_rc_merge() {
        // Literal "ATGC" at 0 and its rc "GCAT" overlapping merge together.
        let p = Pattern::new("atgc", true).unwrap();
        assert_eq!(p.find_merged_matches(b"ATGCAT"), vec![(0, 6)]);
    }

    #[test]
    fn test_merged_matches_multiple_ranges() {
        let p = Pattern::new("tga", false).unwrap();
        assert_eq!(p.find_merged_matches(b"TGACCTTGAC"), vec![(0, 3), (6, 9)]);
    }

    #[test]
    fn test_merged_matches_no_hits() {
        let p = Pattern::new("atgc", true).unwrap();
        assert!(p.find_merged_matches(b"ACGTACGT").is_empty());
    }
}
