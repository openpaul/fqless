use anyhow::{Context, Result};
use bio::io::fastq;
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek};

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
                    // Rewind to start of file
                    file.seek(std::io::SeekFrom::Start(0))?;
                    if magic == [0x1f, 0x8b] {
                        // Gzipped file
                        let gz_decoder = MultiGzDecoder::new(BufReader::new(file));
                        Box::new(BufReader::new(gz_decoder))
                    } else {
                        // Plain text file
                        Box::new(BufReader::new(file))
                    }
                }
                Err(_) => {
                    // Empty file or error - treat as plain text
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

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::FastqRead;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Sample FASTQ record for testing
    const SAMPLE_FASTQ: &str = "@SEQ_ID\nGATTACA\n+\nIIIIIII\n";
    const MULTI_RECORD_FASTQ: &str = "@SEQ1\nGATTACA\n+\nIIIIIII\n@SEQ2\nACGTACGT\n+\n!!!!!!!!";

    /// Helper function to create a temporary FASTQ file
    fn create_temp_fastq(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("Failed to create temp file");
        file.write_all(content.as_bytes())
            .expect("Failed to write to temp file");
        file.flush().expect("Failed to flush temp file");
        file
    }

    /// Helper function to create a temporary gzipped FASTQ file
    fn create_temp_gzipped_fastq(content: &str) -> NamedTempFile {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let mut file = NamedTempFile::new().expect("Failed to create temp file");
        {
            let mut encoder = GzEncoder::new(&mut file, Compression::default());
            encoder
                .write_all(content.as_bytes())
                .expect("Failed to write to encoder");
            encoder.finish().expect("Failed to finish encoding");
        }
        file.flush().expect("Failed to flush temp file");
        file
    }

    #[test]
    fn test_new_with_plain_text_file() {
        let temp_file = create_temp_fastq(SAMPLE_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let reader = FastqReader::new(path);
        assert!(
            reader.is_ok(),
            "Failed to create reader for plain text file"
        );
    }

    #[test]
    fn test_new_with_gzipped_file() {
        let temp_file = create_temp_gzipped_fastq(SAMPLE_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let reader = FastqReader::new(path);
        assert!(
            reader.is_ok(),
            "Failed to create reader for gzipped file: {:?}",
            reader.err()
        );
    }

    #[test]
    fn test_new_with_nonexistent_file() {
        let result = FastqReader::new("/nonexistent/path/to/file.fastq");
        assert!(result.is_err(), "Should fail when opening nonexistent file");

        let error_msg = format!("{}", result.err().unwrap());
        assert!(
            error_msg.contains("Failed to open file"),
            "Error message should mention file opening failure"
        );
    }

    #[test]
    fn test_new_with_empty_file() {
        let temp_file = create_temp_fastq("");
        let path = temp_file.path().to_str().unwrap();

        let reader = FastqReader::new(path);
        assert!(reader.is_ok(), "Should handle empty file gracefully");
    }

    #[test]
    fn test_into_fastq_reader_plain_text() {
        let temp_file = create_temp_fastq(SAMPLE_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let reader = FastqReader::new(path).expect("Failed to create FastqReader");
        let mut fastq_reader = reader.into_fastq_reader();

        let mut record = fastq::Record::new();
        let result = fastq_reader.read(&mut record);

        assert!(result.is_ok(), "Failed to read FASTQ record");
        assert_eq!(record.id(), "SEQ_ID");
        assert_eq!(record.seq(), b"GATTACA");
        assert_eq!(record.qual(), b"IIIIIII");
    }

    #[test]
    fn test_into_fastq_reader_gzipped() {
        let temp_file = create_temp_gzipped_fastq(SAMPLE_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let reader = FastqReader::new(path).expect("Failed to create FastqReader");
        let mut fastq_reader = reader.into_fastq_reader();

        let mut record = fastq::Record::new();
        let result = fastq_reader.read(&mut record);

        assert!(
            result.is_ok(),
            "Failed to read FASTQ record from gzipped file"
        );
        assert_eq!(record.id(), "SEQ_ID");
        assert_eq!(record.seq(), b"GATTACA");
        assert_eq!(record.qual(), b"IIIIIII");
    }

    #[test]
    fn test_read_multiple_records() {
        let temp_file = create_temp_fastq(MULTI_RECORD_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let reader = FastqReader::new(path).expect("Failed to create FastqReader");
        let mut fastq_reader = reader.into_fastq_reader();

        // Read first record
        let mut record1 = fastq::Record::new();
        fastq_reader
            .read(&mut record1)
            .expect("Failed to read first record");
        assert_eq!(record1.id(), "SEQ1");
        assert_eq!(record1.seq(), b"GATTACA");

        // Read second record
        let mut record2 = fastq::Record::new();
        fastq_reader
            .read(&mut record2)
            .expect("Failed to read second record");
        assert_eq!(record2.id(), "SEQ2");
        assert_eq!(record2.seq(), b"ACGTACGT");
    }

    #[test]
    fn test_read_multiple_records_gzipped() {
        let temp_file = create_temp_gzipped_fastq(MULTI_RECORD_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let reader = FastqReader::new(path).expect("Failed to create FastqReader");
        let mut fastq_reader = reader.into_fastq_reader();

        // Read first record
        let mut record1 = fastq::Record::new();
        fastq_reader
            .read(&mut record1)
            .expect("Failed to read first record from gzipped file");
        assert_eq!(record1.id(), "SEQ1");

        // Read second record
        let mut record2 = fastq::Record::new();
        fastq_reader
            .read(&mut record2)
            .expect("Failed to read second record from gzipped file");
        assert_eq!(record2.id(), "SEQ2");
    }

    #[test]
    fn test_magic_bytes_detection_gzip() {
        let temp_file = create_temp_gzipped_fastq(SAMPLE_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        // Read the magic bytes to verify they match gzip signature
        let mut file = File::open(path).expect("Failed to open temp file");
        let mut magic = [0u8; 2];
        file.read_exact(&mut magic)
            .expect("Failed to read magic bytes");

        assert_eq!(
            magic,
            [0x1f, 0x8b],
            "Gzipped file should have correct magic bytes"
        );

        // Now test that reader handles it correctly
        let reader = FastqReader::new(path);
        assert!(
            reader.is_ok(),
            "Should successfully create reader for gzipped file"
        );
    }

    #[test]
    fn test_large_file_handling() {
        // Create a FASTQ file with 10 records (reduced for faster testing)
        let mut large_content = String::new();
        for i in 0..100000 {
            large_content.push_str(&format!(
                "@SEQ_{}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n",
                i
            ));
        }

        let temp_file = create_temp_fastq(&large_content);
        let path = temp_file.path().to_str().unwrap();

        let reader = FastqReader::new(path).expect("Failed to create FastqReader");
        let mut fastq_reader = reader.into_fastq_reader();

        let mut count = 0;
        let mut record = fastq::Record::new();
        // Check if record has actual data (non-empty ID) to detect EOF
        while let Ok(()) = fastq_reader.read(&mut record) {
            if record.is_empty() {
                break;
            }
            count += 1;
        }

        assert_eq!(count, 100000, "Should read all records");
    }

    #[test]
    fn test_file_path_with_context() {
        // Test that error messages include the file path
        let bad_path = "/tmp/definitely_does_not_exist_12345.fastq";
        let result = FastqReader::new(bad_path);

        assert!(result.is_err());
        let error_msg = format!("{}", result.err().unwrap());
        assert!(
            error_msg.contains(bad_path),
            "Error message should include the file path that failed"
        );
    }
}
