use crate::reader::{open_source, RecordSource};
use crate::record::Record;
use anyhow::Result;
use std::sync::{Arc, RwLock};

pub struct DisplayBuffer {
    source: Box<dyn RecordSource>,
    pub reads: Arc<RwLock<Vec<Record>>>,
    buffer_end: u64,
    fully_loaded: bool,
}

/// New buffer from file or stdin
impl DisplayBuffer {
    pub fn new(file_path: &str, reference: Option<&str>) -> Result<Self> {
        let source = open_source(file_path, reference)?;

        let reads = Arc::new(RwLock::new(Vec::new()));

        let buffer_end = 0;
        let fully_loaded = false;

        Ok(DisplayBuffer {
            source,
            reads,
            buffer_end,
            fully_loaded,
        })
    }
    pub fn load_window(&mut self, position: u64, n: usize) -> Result<()> {
        while self.buffer_end < position + n as u64 {
            match self.source.next_record()? {
                Some(record) => {
                    let record = record.clone();
                    self.reads.write().unwrap().push(record);
                    self.buffer_end += 1;
                }
                None => {
                    self.fully_loaded = true;
                    break;
                }
            }
        }
        Ok(())
    }

    /// Implement get_window, which gets the record requested reades plus n trailing it
    pub fn get_window(&mut self, position: u64, n: usize) -> Result<Vec<Record>> {
        // Ensure we have enough records loaded
        if !self.fully_loaded || position + n as u64 > self.buffer_end {
            if position >= self.buffer_end {
                return Ok(vec![]);
            }
            self.load_window(position, n)?;
        }

        // Collect records in the requested range
        let return_until = if self.buffer_end < position + n as u64 {
            self.buffer_end
        } else {
            position + n as u64
        };
        // return reads from the position to the end of the buffer
        Ok(self.reads.read().unwrap()[position as usize..return_until as usize].to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_fastq_file() -> NamedTempFile {
        let content = "@SEQ1\nACGT\n+\nIIII\n@SEQ2\nGGGG\n+\nJJJJ\n@SEQ3\nTTTT\n+\nKKKK\n";
        let mut file = NamedTempFile::new().expect("Failed to create temp file");
        file.write_all(content.as_bytes())
            .expect("Failed to write to temp file");
        file.flush().expect("Failed to flush temp file");
        file
    }

    #[test]
    fn test_buffer_creation() {
        let temp_file = create_test_fastq_file();
        let path = temp_file.path().to_str().unwrap();

        let buffer = DisplayBuffer::new(path, None);
        assert!(buffer.is_ok(), "Should successfully create DisplayBuffer");

        let buffer = buffer.unwrap();
        assert_eq!(
            buffer.reads.read().unwrap().len(),
            0,
            "Buffer should start empty"
        );
    }

    #[test]
    fn test_load_window_and_get_window() {
        let temp_file = create_test_fastq_file();
        let path = temp_file.path().to_str().unwrap();

        let mut buffer = DisplayBuffer::new(path, None).unwrap();
        buffer.load_window(0, 2).unwrap();
        assert_eq!(buffer.reads.read().unwrap().len(), 2);

        let window = buffer.get_window(0, 2).unwrap();
        assert_eq!(window.len(), 2);
        assert_eq!(window[0].id(), "SEQ1");
        assert_eq!(window[0].seq(), b"ACGT");
        assert_eq!(window[1].id(), "SEQ2");

        // Asking past the end returns what is available
        let window = buffer.get_window(1, 10).unwrap();
        assert_eq!(window.len(), 2);
        assert_eq!(window[1].id(), "SEQ3");
    }
}
