use crate::ordinal::{OrdinalIndex, SPACING};
use crate::reader::{Input, RecordSource};
use crate::record::Record;
use anyhow::Result;
use std::collections::VecDeque;
use std::sync::{Arc, RwLock};

/// Reads kept around the requested range, at most.
const MAX_RECORDS: usize = 200_000;
/// Bytes of sequence, quality and names kept, at most. Long reads are large,
/// so a record count alone does not bound memory.
const MAX_BYTES: usize = 256 << 20;
/// Reads before the requested position that are worth keeping for scrolling
/// back; anything further back is skipped while streaming forward.
const KEEP_BEHIND: u64 = 50_000;
/// Seek forward instead of decoding when a checkpoint skips at least this
/// many reads.
const SEEK_AHEAD_MIN: u64 = SPACING / 2;

/// The decoded reads around the current view position.
///
/// Reads are numbered from 0 in file order. The buffer holds a contiguous
/// window `[start, next)` of them, bounded by [`MAX_RECORDS`] and
/// [`MAX_BYTES`]; reads before the requested range are dropped first.
/// Going back before the window seeks to the closest checkpoint of the
/// ordinal index (BAM, CRAM) or reads the file again from the start, and a
/// far jump forward seeks instead of decoding everything in between.
///
/// Stdin cannot be read twice, so for stdin every read stays loaded.
pub struct DisplayBuffer {
    input: Input,
    source: Box<dyn RecordSource>,
    index: Arc<RwLock<OrdinalIndex>>,
    /// Number of the read `source` returns next.
    next: u64,
    /// Number of `reads[0]`.
    start: u64,
    reads: VecDeque<Record>,
    bytes: usize,
    /// Number of reads in the input, once its end has been reached.
    end: Option<u64>,
    keep_all: bool,
    max_records: usize,
    max_bytes: usize,
}

impl DisplayBuffer {
    pub fn new(input: Input, index: Arc<RwLock<OrdinalIndex>>) -> Result<Self> {
        let source = input.open()?;
        let keep_all = input.is_stdin();
        Ok(DisplayBuffer {
            input,
            source,
            index,
            next: 0,
            start: 0,
            reads: VecDeque::new(),
            bytes: 0,
            end: None,
            keep_all,
            max_records: MAX_RECORDS,
            max_bytes: MAX_BYTES,
        })
    }

    /// Number of reads known to exist: every read before this one has been
    /// decoded at some point.
    pub fn available(&self) -> u64 {
        self.next
    }

    /// Number of reads in the input, once the end has been reached.
    pub fn end(&self) -> Option<u64> {
        self.end
    }

    /// True once the end of the input has been reached.
    pub fn is_fully_loaded(&self) -> bool {
        self.end.is_some()
    }

    /// Numbers of the reads currently held.
    pub fn loaded(&self) -> std::ops::Range<u64> {
        self.start..self.next
    }

    /// A held read.
    pub fn get(&self, ordinal: u64) -> Option<&Record> {
        let offset = ordinal.checked_sub(self.start)?;
        self.reads.get(usize::try_from(offset).ok()?)
    }

    /// Held reads from `ordinal` on (none if it is outside the window).
    pub fn iter_from(&self, ordinal: u64) -> impl Iterator<Item = &Record> {
        let skip = if self.loaded().contains(&ordinal) {
            (ordinal - self.start) as usize
        } else {
            self.reads.len()
        };
        self.reads.iter().skip(skip)
    }

    /// Make reads `[from, from + n)` available (as far as the input goes).
    pub fn ensure(&mut self, from: u64, n: usize) -> Result<()> {
        if !self.keep_all {
            self.reposition_for(from)?;
        }
        let until = from.saturating_add(n as u64);
        // The end of the input stays known after moving back from it.
        while self.next < until && self.end != Some(self.next) {
            let Some(record) = self.source.next_record()? else {
                self.end = Some(self.next);
                break;
            };
            let ordinal = self.next;
            self.next += 1;
            if !self.keep_all && ordinal + KEEP_BEHIND < from {
                // Far before the requested range: not worth keeping.
                self.clear_window();
                continue;
            }
            self.bytes += record_bytes(record);
            self.reads.push_back(record.clone());
            self.trim_before(from);
        }
        Ok(())
    }

    /// Reads `[position, position + n)` that exist, loading them if needed.
    pub fn get_window(&mut self, position: u64, n: usize) -> Result<Vec<Record>> {
        self.ensure(position, n)?;
        Ok(self.iter_from(position).take(n).cloned().collect())
    }

    /// Move the source so that reading forward reaches `from` cheaply.
    fn reposition_for(&mut self, from: u64) -> Result<()> {
        let behind = from < self.start;
        let checkpoint = self.index.read().unwrap().floor(from);
        let far_ahead =
            checkpoint.is_some_and(|(ordinal, _)| ordinal >= self.next + SEEK_AHEAD_MIN);
        if !behind && !far_ahead {
            return Ok(());
        }

        if let Some((ordinal, checkpoint)) = checkpoint {
            // A checkpoint inside the current window is no help going
            // forward; going back, any checkpoint before `from` is.
            if (behind || ordinal > self.next) && self.source.seek(checkpoint).is_ok() {
                self.next = ordinal;
                self.clear_window();
                return Ok(());
            }
        }

        if behind {
            // No way back: read the input again from the start.
            self.source = self.input.open()?;
            self.next = 0;
            self.clear_window();
        }
        Ok(())
    }

    /// Empty the window; it starts again at `next`.
    fn clear_window(&mut self) {
        self.reads.clear();
        self.bytes = 0;
        self.start = self.next;
    }

    /// Drop the oldest reads while over budget, but never reads at or after
    /// `from`: the requested range always stays.
    fn trim_before(&mut self, from: u64) {
        if self.keep_all {
            return;
        }
        while (self.reads.len() > self.max_records || self.bytes > self.max_bytes)
            && self.start < from
        {
            let Some(dropped) = self.reads.pop_front() else {
                break;
            };
            self.bytes -= record_bytes(&dropped);
            self.start += 1;
        }
    }
}

fn record_bytes(record: &Record) -> usize {
    record.id().len()
        + record.desc().map_or(0, str::len)
        + record.seq().len()
        + record.qual().map_or(0, <[u8]>::len)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use tempfile::NamedTempFile;

    fn fastq_file(n: usize) -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("Failed to create temp file");
        for i in 0..n {
            writeln!(file, "@read{i}\nACGT\n+\nIIII").unwrap();
        }
        file.flush().unwrap();
        file
    }

    fn buffer(file: &NamedTempFile) -> DisplayBuffer {
        let input = Input::new(file.path().to_str().unwrap(), None);
        DisplayBuffer::new(input, Arc::default()).expect("Failed to open buffer")
    }

    fn ids(records: &[Record]) -> Vec<&str> {
        records.iter().map(Record::id).collect()
    }

    #[test]
    fn test_window_and_end() {
        let file = fastq_file(3);
        let mut buffer = buffer(&file);
        assert_eq!(buffer.available(), 0);

        let window = buffer.get_window(0, 2).unwrap();
        assert_eq!(ids(&window), ["read0", "read1"]);
        assert_eq!(buffer.end(), None);

        let window = buffer.get_window(1, 10).unwrap();
        assert_eq!(ids(&window), ["read1", "read2"]);
        assert_eq!(buffer.end(), Some(3));
        assert!(buffer.get_window(5, 2).unwrap().is_empty());

        // Going back to the start after reaching the end still reads.
        buffer.max_records = 1;
        buffer.ensure(2, 1).unwrap();
        assert_eq!(ids(&buffer.get_window(0, 2).unwrap()), ["read0", "read1"]);
    }

    #[test]
    fn test_memory_is_bounded_and_going_back_rereads() {
        let file = fastq_file(1_000);
        let mut buffer = buffer(&file);
        buffer.max_records = 100;

        let window = buffer.get_window(900, 5).unwrap();
        assert_eq!(
            ids(&window),
            ["read900", "read901", "read902", "read903", "read904"]
        );
        assert!(
            buffer.reads.len() <= 100,
            "kept {} reads",
            buffer.reads.len()
        );
        assert!(buffer.loaded().start > 0);

        // Before the window: the file is read again, with the same result.
        let window = buffer.get_window(3, 2).unwrap();
        assert_eq!(ids(&window), ["read3", "read4"]);
        assert_eq!(buffer.get(3).map(Record::id), Some("read3"));
    }

    #[test]
    fn test_requested_range_is_never_trimmed() {
        let file = fastq_file(50);
        let mut buffer = buffer(&file);
        buffer.max_records = 5;
        let window = buffer.get_window(10, 20).unwrap();
        assert_eq!(window.len(), 20);
        assert_eq!(window[0].id(), "read10");
    }

    /// Numbered reads; a checkpoint is the read number. Counts decoded reads.
    struct Numbers {
        next: u64,
        total: u64,
        decoded: Arc<AtomicUsize>,
        record: Record,
    }

    impl RecordSource for Numbers {
        fn next_record(&mut self) -> Result<Option<&Record>> {
            if self.next == self.total {
                return Ok(None);
            }
            self.decoded.fetch_add(1, Ordering::Relaxed);
            self.record = Record::with_attrs(&self.next.to_string(), None, b"A", None);
            self.next += 1;
            Ok(Some(&self.record))
        }

        fn checkpoint(&mut self) -> Option<u64> {
            Some(self.next)
        }

        fn seek(&mut self, checkpoint: u64) -> Result<()> {
            self.next = checkpoint;
            Ok(())
        }
    }

    #[test]
    fn test_seeks_with_the_ordinal_index() {
        const TOTAL: u64 = 1_000_000;
        let decoded = Arc::new(AtomicUsize::new(0));
        let mut index = OrdinalIndex::default();
        for ordinal in (SPACING..TOTAL).step_by(SPACING as usize) {
            index.push(ordinal, ordinal);
        }
        let mut buffer = DisplayBuffer {
            input: Input::new("unused", None),
            source: Box::new(Numbers {
                next: 0,
                total: TOTAL,
                decoded: Arc::clone(&decoded),
                record: Record::new(),
            }),
            index: Arc::new(RwLock::new(index)),
            next: 0,
            start: 0,
            reads: VecDeque::new(),
            bytes: 0,
            end: None,
            keep_all: false,
            max_records: MAX_RECORDS,
            max_bytes: MAX_BYTES,
        };

        // Far ahead: seek to the checkpoint at 990,000 and decode from there.
        let window = buffer.get_window(995_000, 3).unwrap();
        assert_eq!(ids(&window), ["995000", "995001", "995002"]);
        assert!(decoded.load(Ordering::Relaxed) < 2 * SPACING as usize);

        // Back before the window: seek to the checkpoint before read 12,345.
        decoded.store(0, Ordering::Relaxed);
        let window = buffer.get_window(12_345, 1).unwrap();
        assert_eq!(ids(&window), ["12345"]);
        assert!(decoded.load(Ordering::Relaxed) < 2 * SPACING as usize);

        // Reaching the end records the total.
        buffer.ensure(TOTAL - 1, 10).unwrap();
        assert_eq!(buffer.end(), Some(TOTAL));
    }

    #[test]
    fn test_stdin_keeps_everything() {
        let file = fastq_file(10);
        let mut buffer = buffer(&file);
        buffer.keep_all = true;
        buffer.max_records = 2;
        buffer.ensure(8, 2).unwrap();
        assert_eq!(buffer.loaded(), 0..10);
    }
}
