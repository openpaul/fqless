//! Sparse read-number index for random access.
//!
//! `.bai`/`.csi`/`.crai` map genomic coordinates to file offsets; they cannot
//! answer "where does read number N start", which is how fqless navigates.
//! This index answers that: while the read counter scans a seekable file
//! (BAM, CRAM) it records a checkpoint about every [`SPACING`] reads. Jumping
//! to read N then means seeking to the closest checkpoint before it and
//! decoding forward at most a few thousand reads.

/// Reads between two checkpoints. 100M reads need 10,000 entries (160 KB).
pub const SPACING: u64 = 10_000;

#[derive(Debug, Default)]
pub struct OrdinalIndex {
    /// (read number, checkpoint of that read), sorted by read number.
    entries: Vec<(u64, u64)>,
}

impl OrdinalIndex {
    pub fn push(&mut self, ordinal: u64, checkpoint: u64) {
        debug_assert!(self.entries.last().is_none_or(|&(last, _)| last < ordinal));
        self.entries.push((ordinal, checkpoint));
    }

    /// The closest checkpoint at or before read `ordinal`.
    pub fn floor(&self, ordinal: u64) -> Option<(u64, u64)> {
        let after = self.entries.partition_point(|&(o, _)| o <= ordinal);
        after.checked_sub(1).map(|i| self.entries[i])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_floor() {
        let mut index = OrdinalIndex::default();
        assert_eq!(index.floor(5), None);
        index.push(10, 1000);
        index.push(20, 2000);
        assert_eq!(index.floor(9), None);
        assert_eq!(index.floor(10), Some((10, 1000)));
        assert_eq!(index.floor(19), Some((10, 1000)));
        assert_eq!(index.floor(25), Some((20, 2000)));
    }
}
