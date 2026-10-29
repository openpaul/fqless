use crate::record::Record;
use anyhow::{Context, Result};
use bio::alphabets::dna;
use bio::io::fastq::{self, FastqRead};
use flate2::bufread::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek};

/// A source of sequencing records, independent of the file format.
///
/// The returned reference is only valid until the next call; callers that need
/// to keep a record around must clone it. Reusing one record keeps whole-file
/// scans allocation-free.
pub trait RecordSource {
    fn next_record(&mut self) -> Result<Option<&Record>>;
}

/// Open `path` (or stdin for `-`) and choose a parser from its contents.
///
/// The format is detected from the first non-blank byte of the (decompressed)
/// stream, so `.fa`/`.fastq` naming and gzip suffixes do not matter. SAM, BAM
/// and CRAM are detected here as well once they are implemented.
pub fn open_source(path: &str) -> Result<Box<dyn RecordSource>> {
    let mut input = InputReader::new(path)?;
    let is_fasta = {
        let prefix = input.peek(4096)?;
        // A UTF-8 BOM is legal at the start of a text file.
        let prefix = prefix.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(prefix);
        prefix
            .iter()
            .find(|b| !b.is_ascii_whitespace())
            .is_some_and(|&b| b == b'>')
    };

    if is_fasta {
        return Ok(Box::new(FastaSource::new(input.into_inner())));
    }

    let is_sam = {
        let prefix = input.peek(4096)?;
        let prefix = prefix.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(prefix);
        looks_like_sam(prefix)
    };

    if is_sam {
        Ok(Box::new(SamSource::new(input.into_inner())))
    } else {
        Ok(Box::new(FastqSource::new(input.into_fastq_reader())))
    }
}

/// True when the first line looks like a SAM header or alignment record.
///
/// SAM headers start with `@HD`, `@SQ`, `@RG`, `@PG` or `@CO` followed by a
/// tab, and alignment records have at least 11 tab-separated fields. A FASTQ
/// header may also start with `@`, so the tag/field check decides.
fn looks_like_sam(prefix: &[u8]) -> bool {
    let line = match prefix.split(|&b| b == b'\n').next() {
        Some(line) => trim_ascii_end(line),
        None => return false,
    };

    const HEADER_TAGS: [&[u8]; 5] = [b"@HD", b"@SQ", b"@RG", b"@PG", b"@CO"];
    let is_header = HEADER_TAGS
        .iter()
        .any(|tag| line.starts_with(tag) && (line.len() == tag.len() || line[tag.len()] == b'\t'));

    is_header || line.iter().filter(|&&b| b == b'\t').count() >= 10
}

/// Trim trailing ASCII whitespace (`\n`, `\r`, spaces).
fn trim_ascii_end(line: &[u8]) -> &[u8] {
    let end = line
        .iter()
        .rposition(|b| !b.is_ascii_whitespace())
        .map_or(0, |pos| pos + 1);
    &line[..end]
}

/// Unified reader that handles both files and stdin, with automatic compression detection.
///
/// This only produces decompressed bytes; interpreting them as a sequencing
/// format is the job of [`RecordSource`] implementations.
pub struct InputReader {
    inner: Box<dyn BufRead>,
}

impl InputReader {
    /// Create a new reader from a file path or "-" for stdin
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
            Ok(InputReader { inner: reader })
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

            Ok(InputReader { inner: reader })
        }
    }

    /// Wrap the decompressed stream in a FASTQ reader.
    pub fn into_fastq_reader(self) -> fastq::Reader<BufReader<Box<dyn BufRead>>> {
        fastq::Reader::new(self.inner)
    }

    /// Take the decompressed stream out, for parsers that read lines directly.
    pub fn into_inner(self) -> Box<dyn BufRead> {
        self.inner
    }

    /// Look at the next `n` bytes without consuming them.
    ///
    /// Works on the decompressed stream, so gzipped input is sniffed by its
    /// actual content. The returned slice may be shorter than `n` near EOF.
    pub fn peek(&mut self, n: usize) -> Result<&[u8]> {
        let buffer = self.inner.fill_buf()?;
        Ok(&buffer[..buffer.len().min(n)])
    }
}

/// FASTQ parser built on `bio`'s reader, decoding into [`Record`].
pub struct FastqSource {
    reader: fastq::Reader<BufReader<Box<dyn BufRead>>>,
    scratch: fastq::Record,
    record: Record,
    done: bool,
}

impl FastqSource {
    pub fn new(reader: fastq::Reader<BufReader<Box<dyn BufRead>>>) -> Self {
        Self {
            reader,
            scratch: fastq::Record::new(),
            record: Record::new(),
            done: false,
        }
    }
}

impl RecordSource for FastqSource {
    fn next_record(&mut self) -> Result<Option<&Record>> {
        if self.done {
            return Ok(None);
        }
        self.reader.read(&mut self.scratch)?;
        if self.scratch.is_empty() {
            self.done = true;
            return Ok(None);
        }
        self.record.set_from_fastq(&self.scratch);
        Ok(Some(&self.record))
    }
}

/// FASTA parser: `>name [description]` followed by one or more sequence lines.
///
/// Sequences may be wrapped over any number of lines; they are concatenated
/// and internal whitespace is dropped. Case is preserved so soft-masked
/// (lowercase) sequence keeps working with the adapter detector.
pub struct FastaSource {
    inner: Box<dyn BufRead>,
    /// Header of the record after the current one, found while reading ahead.
    pending_header: Option<Vec<u8>>,
    line: Vec<u8>,
    record: Record,
    done: bool,
}

impl FastaSource {
    pub fn new(inner: Box<dyn BufRead>) -> Self {
        Self {
            inner,
            pending_header: None,
            line: Vec::new(),
            record: Record::new(),
            done: false,
        }
    }

    /// Read one line into `self.line` and return its length without the line
    /// ending and trailing whitespace, or `None` at EOF.
    fn read_line_len(&mut self) -> Result<Option<usize>> {
        self.line.clear();
        if self.inner.read_until(b'\n', &mut self.line)? == 0 {
            return Ok(None);
        }
        let end = self
            .line
            .iter()
            .rposition(|b| !b.is_ascii_whitespace())
            .map_or(0, |pos| pos + 1);
        Ok(Some(end))
    }

    /// Skip blank lines and return the next header line without its `>`.
    fn read_header(&mut self) -> Result<Option<Vec<u8>>> {
        while let Some(len) = self.read_line_len()? {
            let line = &self.line[..len];
            if line.is_empty() {
                continue;
            }
            let line = line.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(line);
            match line.strip_prefix(b">") {
                Some(header) => return Ok(Some(header.to_vec())),
                None => anyhow::bail!(
                    "expected '>' at FASTA record start, found {:?}",
                    String::from_utf8_lossy(&line[..line.len().min(40)])
                ),
            }
        }
        Ok(None)
    }
}

/// Split a FASTA header into ID and description at the first whitespace.
fn split_fasta_header(header: &[u8]) -> Result<(&str, Option<&str>)> {
    let header = std::str::from_utf8(header).context("FASTA header is not valid UTF-8")?;
    let header = header.trim();
    let (id, desc) = match header.find(|c: char| c.is_ascii_whitespace()) {
        Some(pos) => (&header[..pos], Some(header[pos..].trim())),
        None => (header, None),
    };
    if id.is_empty() {
        anyhow::bail!("FASTA header is missing an identifier");
    }
    Ok((id, desc.filter(|desc| !desc.is_empty())))
}

impl RecordSource for FastaSource {
    fn next_record(&mut self) -> Result<Option<&Record>> {
        if self.done {
            return Ok(None);
        }

        let header = match self.pending_header.take() {
            Some(header) => header,
            None => match self.read_header()? {
                Some(header) => header,
                None => {
                    self.done = true;
                    return Ok(None);
                }
            },
        };

        let (id, desc) = split_fasta_header(&header)?;
        self.record.set_name(id, desc);
        self.record.seq_mut().clear();
        self.record.set_qual(None);

        while let Some(len) = self.read_line_len()? {
            let line = &self.line[..len];
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix(b">") {
                self.pending_header = Some(rest.to_vec());
                break;
            }
            self.record
                .seq_mut()
                .extend(line.iter().copied().filter(|b| !b.is_ascii_whitespace()));
        }

        Ok(Some(&self.record))
    }
}

/// SAM flag: the sequence is stored as the reverse complement of the read.
const SAM_FLAG_REVERSE: u16 = 0x10;
/// SAM flag: secondary alignment (alternative position of the same read).
const SAM_FLAG_SECONDARY: u16 = 0x100;
/// SAM flag: supplementary alignment (chimeric segment of the same read).
const SAM_FLAG_SUPPLEMENTARY: u16 = 0x800;

/// SAM parser: tab-separated alignment records after an optional header.
///
/// Only sequence-level fields are used (QNAME, FLAG, SEQ, QUAL). Secondary and
/// supplementary alignments are skipped because they repeat a read that is
/// already shown, and reads mapped to the reverse strand are reverse
/// complemented back to their original orientation.
pub struct SamSource {
    inner: Box<dyn BufRead>,
    line: Vec<u8>,
    record: Record,
    line_no: u64,
    done: bool,
}

impl SamSource {
    pub fn new(inner: Box<dyn BufRead>) -> Self {
        Self {
            inner,
            line: Vec::new(),
            record: Record::new(),
            line_no: 0,
            done: false,
        }
    }

    fn malformed(&self, what: &str) -> anyhow::Error {
        anyhow::anyhow!("SAM line {}: {}", self.line_no, what)
    }
}

impl RecordSource for SamSource {
    fn next_record(&mut self) -> Result<Option<&Record>> {
        if self.done {
            return Ok(None);
        }

        loop {
            self.line.clear();
            if self.inner.read_until(b'\n', &mut self.line)? == 0 {
                self.done = true;
                return Ok(None);
            }
            self.line_no += 1;

            let line = trim_ascii_end(&self.line);
            if line.is_empty() {
                continue;
            }
            // Headers and comments only occur before the alignments, but
            // skipping them anywhere is harmless.
            if line[0] == b'@' {
                continue;
            }

            let mut fields = line.split(|&b| b == b'\t');
            let qname = fields
                .next()
                .filter(|qname| !qname.is_empty())
                .ok_or_else(|| self.malformed("missing read name"))?;
            let flag = fields
                .next()
                .ok_or_else(|| self.malformed("missing FLAG field"))?;
            let flag: u16 = std::str::from_utf8(flag)
                .ok()
                .and_then(|flag| flag.parse().ok())
                .ok_or_else(|| self.malformed("FLAG is not a number"))?;

            // RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN are not displayed.
            for _ in 0..7 {
                fields
                    .next()
                    .ok_or_else(|| self.malformed("expected at least 11 fields"))?;
            }

            let seq = fields
                .next()
                .ok_or_else(|| self.malformed("missing SEQ field"))?;
            let qual = fields
                .next()
                .ok_or_else(|| self.malformed("missing QUAL field"))?;

            if flag & (SAM_FLAG_SECONDARY | SAM_FLAG_SUPPLEMENTARY) != 0 {
                continue;
            }
            if seq == b"*" {
                // An alignment without a sequence carries no read.
                continue;
            }

            let id = std::str::from_utf8(qname)
                .map_err(|_| self.malformed("read name is not valid UTF-8"))?;
            self.record.set_name(id, None);
            self.record.seq_mut().clear();
            self.record.seq_mut().extend_from_slice(seq);
            self.record.set_qual((qual != b"*").then_some(qual));

            if flag & SAM_FLAG_REVERSE != 0 {
                let revcomp = dna::revcomp(self.record.seq());
                let seq = self.record.seq_mut();
                seq.clear();
                seq.extend_from_slice(&revcomp);
                if let Some(qual) = self.record.qual_mut() {
                    qual.reverse();
                }
            }

            return Ok(Some(&self.record));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

    fn next_owned(source: &mut dyn RecordSource) -> Result<Option<Record>> {
        Ok(source.next_record()?.cloned())
    }

    #[test]
    fn test_new_with_plain_text_file() {
        let temp_file = create_temp_fastq(SAMPLE_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let reader = InputReader::new(path);
        assert!(
            reader.is_ok(),
            "Failed to create reader for plain text file"
        );
    }

    #[test]
    fn test_new_with_gzipped_file() {
        let temp_file = create_temp_gzipped_fastq(SAMPLE_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let reader = InputReader::new(path);
        assert!(
            reader.is_ok(),
            "Failed to create reader for gzipped file: {:?}",
            reader.err()
        );
    }

    #[test]
    fn test_new_with_nonexistent_file() {
        let result = InputReader::new("/nonexistent/path/to/file.fastq");
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

        let reader = InputReader::new(path);
        assert!(reader.is_ok(), "Should handle empty file gracefully");
    }

    #[test]
    fn test_open_source_plain_text() {
        let temp_file = create_temp_fastq(SAMPLE_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open source");
        let record = source.next_record().expect("Failed to read record");

        let record = record.expect("Expected one record");
        assert_eq!(record.id(), "SEQ_ID");
        assert_eq!(record.seq(), b"GATTACA");
        assert_eq!(record.qual(), Some(b"IIIIIII".as_slice()));
        assert!(source.next_record().unwrap().is_none());
    }

    #[test]
    fn test_open_source_gzipped() {
        let temp_file = create_temp_gzipped_fastq(SAMPLE_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open source");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.id(), "SEQ_ID");
        assert_eq!(record.seq(), b"GATTACA");
        assert_eq!(record.qual(), Some(b"IIIIIII".as_slice()));
    }

    #[test]
    fn test_open_source_multiple_records() {
        let temp_file = create_temp_fastq(MULTI_RECORD_FASTQ);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open source");

        let first = next_owned(source.as_mut())
            .expect("Failed to read first record")
            .expect("Expected a first record");
        assert_eq!(first.id(), "SEQ1");
        assert_eq!(first.seq(), b"GATTACA");

        let second = next_owned(source.as_mut())
            .expect("Failed to read second record")
            .expect("Expected a second record");
        assert_eq!(second.id(), "SEQ2");
        assert_eq!(second.seq(), b"ACGTACGT");

        assert!(source.next_record().unwrap().is_none());
    }

    #[test]
    fn test_open_source_parses_description() {
        let temp_file = create_temp_fastq("@read1 some description\nACGT\n+\nIIII\n");
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open source");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.id(), "read1");
        assert_eq!(record.desc(), Some("some description"));
    }

    #[test]
    fn test_open_source_multiline_sequence() {
        let temp_file = create_temp_fastq("@read1\nACGT\nACGT\n+\nIIIIIIII\n");
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open source");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.seq(), b"ACGTACGT");
        assert_eq!(record.qual(), Some(b"IIIIIIII".as_slice()));
    }

    #[test]
    fn test_record_source_reuses_record_buffers() {
        let content = "@A\nAAAA\n+\nIIII\n@B\nCCCC\n+\nIIII\n";
        let temp_file = create_temp_fastq(content);
        let path = temp_file.path().to_str().unwrap();
        let mut source = open_source(path).expect("Failed to open source");

        let (id1, seq_ptr1) = {
            let record = source.next_record().unwrap().unwrap();
            (record.id().to_string(), record.seq().as_ptr())
        };
        let (id2, seq_ptr2) = {
            let record = source.next_record().unwrap().unwrap();
            (record.id().to_string(), record.seq().as_ptr())
        };

        assert_eq!(id1, "A");
        assert_eq!(id2, "B");
        assert_eq!(seq_ptr1, seq_ptr2, "record buffer should be reused");
        assert!(source.next_record().unwrap().is_none());
    }

    #[test]
    fn test_open_source_fasta_plain() {
        let content = ">seq1 first record\nACGT\nACGT\n\n>seq2\nTTTT\n";
        let temp_file = create_temp_fastq(content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open FASTA");

        let first = source
            .next_record()
            .expect("Failed to read first record")
            .expect("Expected a first record");
        assert_eq!(first.id(), "seq1");
        assert_eq!(first.desc(), Some("first record"));
        assert_eq!(first.seq(), b"ACGTACGT");
        assert_eq!(first.qual(), None);
        assert!(!first.has_quality());

        let second = source
            .next_record()
            .expect("Failed to read second record")
            .expect("Expected a second record");
        assert_eq!(second.id(), "seq2");
        assert_eq!(second.desc(), None);
        assert_eq!(second.seq(), b"TTTT");

        assert!(source.next_record().unwrap().is_none());
    }

    #[test]
    fn test_open_source_fasta_gzipped() {
        let content = ">seq1\nACGT\n";
        let temp_file = create_temp_gzipped_fastq(content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open gzipped FASTA");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.id(), "seq1");
        assert_eq!(record.seq(), b"ACGT");
        assert_eq!(record.qual(), None);
    }

    #[test]
    fn test_open_source_fasta_without_trailing_newline() {
        let temp_file = create_temp_fastq(">seq1\nACGT");
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open FASTA");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.id(), "seq1");
        assert_eq!(record.seq(), b"ACGT");
        assert!(source.next_record().unwrap().is_none());
    }

    #[test]
    fn test_open_source_fasta_crlf_and_blank_header_lines() {
        let content = "\r\n>seq1 desc\r\nACGT\r\nACGT\r\n";
        let temp_file = create_temp_fastq(content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open FASTA");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.id(), "seq1");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ACGTACGT");
    }

    #[test]
    fn test_open_source_fasta_with_utf8_bom() {
        let temp_file = create_temp_fastq("\u{feff}>seq1 desc\nACGT\n");
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open FASTA with BOM");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.id(), "seq1");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ACGT");
    }

    #[test]
    fn test_open_source_fasta_drops_internal_whitespace() {
        let temp_file = create_temp_fastq(">seq1\nAC GT\nAC\tGT\n");
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open FASTA");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.seq(), b"ACGTACGT");
    }

    #[test]
    fn test_fasta_source_rejects_sequence_before_header() {
        let input: Box<dyn BufRead> = Box::new(BufReader::new(&b"ACGT\n>seq1\nACGT\n"[..]));
        let mut source = FastaSource::new(input);
        let error = source
            .next_record()
            .expect_err("Sequence before the first header should fail");

        assert!(
            format!("{error}").contains("expected '>'"),
            "Unexpected error: {error}"
        );
    }

    #[test]
    fn test_open_source_fasta_empty_file() {
        let temp_file = create_temp_fastq("");
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open empty file");
        assert!(source.next_record().unwrap().is_none());
    }

    const SAM_HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";

    #[test]
    fn test_open_source_sam_with_header() {
        let content = format!(
            "{SAM_HEADER}read1\t0\tchr1\t7\t60\t4M\t*\t0\t0\tACGT\tIIII\nread2\t0\tchr1\t8\t60\t3M\t*\t0\t0\tGGG\t###\n"
        );
        let temp_file = create_temp_fastq(&content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open SAM");

        let first = source
            .next_record()
            .expect("Failed to read first record")
            .expect("Expected a first record");
        assert_eq!(first.id(), "read1");
        assert_eq!(first.seq(), b"ACGT");
        assert_eq!(first.qual(), Some(b"IIII".as_slice()));

        let second = source
            .next_record()
            .expect("Failed to read second record")
            .expect("Expected a second record");
        assert_eq!(second.id(), "read2");
        assert_eq!(second.seq(), b"GGG");

        assert!(source.next_record().unwrap().is_none());
    }

    #[test]
    fn test_open_source_sam_without_header() {
        let content = "read1\t0\tchr1\t7\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let temp_file = create_temp_fastq(content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open headerless SAM");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.id(), "read1");
        assert_eq!(record.seq(), b"ACGT");
    }

    #[test]
    fn test_open_source_fastq_starting_with_hd_is_not_sam() {
        let temp_file = create_temp_fastq("@HDREAD\nACGT\n+\nIIII\n");
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open FASTQ");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.id(), "HDREAD");
        assert_eq!(record.qual(), Some(b"IIII".as_slice()));
    }

    #[test]
    fn test_sam_skips_secondary_and_supplementary() {
        let content = format!(
            "{SAM_HEADER}primary\t0\tchr1\t7\t60\t4M\t*\t0\t0\tACGT\tIIII\nsecondary\t256\tchr1\t7\t60\t4M\t*\t0\t0\tACGT\tIIII\nsupplementary\t2048\tchr1\t9\t60\t2M\t*\t0\t0\tAC\tII\n"
        );
        let temp_file = create_temp_fastq(&content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open SAM");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected the primary record");

        assert_eq!(record.id(), "primary");
        assert!(source.next_record().unwrap().is_none());
    }

    #[test]
    fn test_sam_reverse_strand_restores_original_orientation() {
        // SAM stores reverse-strand reads reverse complemented, so the viewer
        // should show the original read again after flipping it back.
        let original = b"ACGTACGTTT";
        let stored = dna::revcomp(original);
        let stored = std::str::from_utf8(&stored).unwrap();
        let content =
            format!("{SAM_HEADER}rev\t16\tchr1\t7\t60\t10M\t*\t0\t0\t{stored}\tABCDEFGHIJ\n");
        let temp_file = create_temp_fastq(&content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open SAM");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.seq(), original);
        assert_eq!(record.qual(), Some(b"JIHGFEDCBA".as_slice()));
    }

    #[test]
    fn test_sam_missing_quality_is_none() {
        let content = format!("{SAM_HEADER}noqual\t0\tchr1\t7\t60\t4M\t*\t0\t0\tACGT\t*\n");
        let temp_file = create_temp_fastq(&content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open SAM");
        let record = source
            .next_record()
            .expect("Failed to read record")
            .expect("Expected one record");

        assert_eq!(record.seq(), b"ACGT");
        assert_eq!(record.qual(), None);
        assert!(!record.has_quality());
    }

    #[test]
    fn test_sam_rejects_line_with_too_few_fields() {
        let content = format!("{SAM_HEADER}read1\t0\tchr1\t7\n");
        let temp_file = create_temp_fastq(&content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open SAM");
        let error = source
            .next_record()
            .expect_err("A truncated SAM line should fail");

        assert!(
            format!("{error}").contains("SAM line 3"),
            "Unexpected error: {error}"
        );
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
        let reader = InputReader::new(path);
        assert!(
            reader.is_ok(),
            "Should successfully create reader for gzipped file"
        );
    }

    #[test]
    fn test_large_file_handling() {
        // Create a FASTQ file with 100000 records
        let mut large_content = String::new();
        for i in 0..100000 {
            large_content.push_str(&format!(
                "@SEQ_{}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n",
                i
            ));
        }

        let temp_file = create_temp_fastq(&large_content);
        let path = temp_file.path().to_str().unwrap();

        let mut source = open_source(path).expect("Failed to open source");

        let mut count = 0;
        while source.next_record().unwrap().is_some() {
            count += 1;
        }

        assert_eq!(count, 100000, "Should read all records");
    }

    #[test]
    fn test_file_path_with_context() {
        // Test that error messages include the file path
        let bad_path = "/tmp/definitely_does_not_exist_12345.fastq";
        let result = InputReader::new(bad_path);

        assert!(result.is_err());
        let error_msg = format!("{}", result.err().unwrap());
        assert!(
            error_msg.contains(bad_path),
            "Error message should include the file path that failed"
        );
    }
}
