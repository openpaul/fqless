use anyhow::{Context, Result};
use bio::io::fastq;
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

/// Unified reader that handles both files and stdin, with automatic compression detection
pub struct FastqReader {
    inner: Box<dyn BufRead>,
}

impl FastqReader {
    /// Create a new FastqReader from a file path or "-" for stdin
    pub fn new(input: &str) -> Result<Self> {
        if input == "-" {
            let stdin = io::stdin();
            let mut buf_reader = BufReader::new(stdin.lock());

            let buffer = buf_reader.fill_buf()?;
            let reader: Box<dyn BufRead> = if buffer.starts_with(&[0x1F, 0x8B]) {
                Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
            } else {
                Box::new(buf_reader)
            };
            Ok(FastqReader { inner: reader })
        } else {
            // Handle file
            let mut file =
                File::open(input).with_context(|| format!("Failed to open file: {}", input))?;

            // Check if file is gzipped by reading magic bytes
            let mut magic = [0; 2];
            let reader: Box<dyn BufRead> = match file.read_exact(&mut magic) {
                Ok(_) => {
                    if magic == [0x1f, 0x8b] {
                        // Gzipped file - reopen to get a clean handle
                        let file = File::open(input)?;
                        let gz_decoder = MultiGzDecoder::new(BufReader::new(file));
                        Box::new(BufReader::new(gz_decoder))
                    } else {
                        // Plain text file - reopen to get a clean handle
                        let file = File::open(input)?;
                        Box::new(BufReader::new(file))
                    }
                }
                Err(_) => {
                    // Empty file or error - treat as plain text
                    let file = File::open(input)?;
                    Box::new(BufReader::new(file))
                }
            };

            Ok(FastqReader { inner: reader })
        }
    }

    /// Get a new FASTQ reader for this input (consumes self)
    pub fn into_fastq_reader(self) -> fastq::Reader<BufReader<Box<dyn BufRead>>> {
        fastq::Reader::new(self.inner)
    }
}
