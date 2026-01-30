use crate::reader::FastqReader;
use anyhow::Result;
use bio::io::fastq;
use std::{
    io::{BufRead, BufReader},
    sync::{Arc, RwLock},
};

pub struct DisplayBuffer {
    records: fastq::Records<BufReader<Box<dyn BufRead>>>,
    pub reads: Arc<RwLock<Vec<fastq::Record>>>,
    buffer_end: u64,
    fully_loaded: bool,
}

/// New buffer from file or stdin
impl DisplayBuffer {
    pub fn new(file_path: &str) -> Result<Self> {
        let reader = FastqReader::new(file_path)?;
        let fastq_reader = reader.into_fastq_reader();
        let records = fastq_reader.records();

        let reads = Arc::new(RwLock::new(Vec::new()));

        let buffer_end = 0;
        let fully_loaded = false;

        Ok(DisplayBuffer {
            records,
            reads,
            buffer_end,
            fully_loaded,
        })
    }
    pub fn load_window(&mut self, position: u64, n: usize) -> Result<()> {
        while self.buffer_end < position + n as u64 {
            if let Some(record) = self.records.next() {
                let record = record?;
                self.reads.write().unwrap().push(record);
                self.buffer_end += 1;
            } else {
                self.fully_loaded = true;
                break;
            }
        }
        Ok(())
    }

    /// Implement get_window, which gets the record requested reades plus n trailing it
    pub fn get_window(&mut self, position: u64, n: usize) -> Result<Vec<fastq::Record>> {
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

        let buffer = DisplayBuffer::new(path);
        assert!(buffer.is_ok(), "Should successfully create DisplayBuffer");

        let buffer = buffer.unwrap();
        assert_eq!(
            buffer.reads.read().unwrap().len(),
            0,
            "Buffer should start empty"
        );
    }
}
