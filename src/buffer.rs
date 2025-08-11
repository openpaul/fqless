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

        return Ok(DisplayBuffer {
            records,
            reads,
            buffer_end,
            fully_loaded,
        });
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
