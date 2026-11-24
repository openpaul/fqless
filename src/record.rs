//! Format-independent sequencing record.
//!
//! Every input format fqless supports is reduced to this type, so the viewer,
//! search and statistics code never has to know whether a read came from
//! FASTQ, FASTA, SAM, BAM or CRAM.
//!
//! Quality is optional because not every format stores it: FASTA has none and
//! BAM/CRAM may omit it per record. `Some` means the qualities are present and
//! exactly as long as the sequence.

use bio::io::fastq;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
}

impl Record {
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a record from already-owned data.
    ///
    /// Test helper for now; the FASTA and alignment parsers will use it once
    /// they land, at which point it stops being test-only.
    #[cfg(test)]
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: &[u8], qual: Option<&[u8]>) -> Self {
        Self {
            id: id.to_owned(),
            desc: desc.map(str::to_owned),
            seq: seq.to_vec(),
            qual: qual.map(<[u8]>::to_vec),
        }
    }

    /// Record name (FASTQ/FASTA ID, SAM QNAME).
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Description following the ID, if the format has one.
    pub fn desc(&self) -> Option<&str> {
        self.desc.as_deref()
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    /// Base qualities, or `None` when the format does not provide them.
    pub fn qual(&self) -> Option<&[u8]> {
        self.qual.as_deref()
    }

    pub fn has_quality(&self) -> bool {
        self.qual.is_some()
    }

    /// Set the name fields, reusing the existing buffers.
    pub(crate) fn set_name(&mut self, id: &str, desc: Option<&str>) {
        self.id.clear();
        self.id.push_str(id);

        match (desc, &mut self.desc) {
            (Some(desc), Some(existing)) => {
                existing.clear();
                existing.push_str(desc);
            }
            (Some(desc), slot) => *slot = Some(desc.to_owned()),
            (None, slot) => *slot = None,
        }
    }

    /// Mutable access to the sequence, so parsers can append without copying.
    pub(crate) fn seq_mut(&mut self) -> &mut Vec<u8> {
        &mut self.seq
    }

    /// Mutable access to the qualities, if the record has any.
    pub(crate) fn qual_mut(&mut self) -> Option<&mut Vec<u8>> {
        self.qual.as_mut()
    }

    /// Replace the qualities with `chars`, reusing the existing buffer.
    pub(crate) fn set_qual_from(&mut self, chars: impl IntoIterator<Item = u8>) {
        let qual = self.qual.get_or_insert_with(Vec::new);
        qual.clear();
        qual.extend(chars);
    }

    /// Set or clear the qualities, reusing the existing buffer when possible.
    pub(crate) fn set_qual(&mut self, qual: Option<&[u8]>) {
        match (qual, &mut self.qual) {
            (Some(qual), Some(existing)) => {
                existing.clear();
                existing.extend_from_slice(qual);
            }
            (Some(qual), slot) => *slot = Some(qual.to_vec()),
            (None, slot) => *slot = None,
        }
    }

    /// Replace the contents with `src`, reusing the existing buffers so a
    /// whole-file scan does not allocate per read.
    pub(crate) fn set_from_fastq(&mut self, src: &fastq::Record) {
        self.set_name(src.id(), src.desc());

        self.seq.clear();
        self.seq.extend_from_slice(src.seq());

        self.set_qual(Some(src.qual()));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_with_attrs_roundtrip() {
        let record = Record::with_attrs("read1", Some("extra info"), b"ACGT", Some(b"IIII"));

        assert_eq!(record.id(), "read1");
        assert_eq!(record.desc(), Some("extra info"));
        assert_eq!(record.seq(), b"ACGT");
        assert_eq!(record.qual(), Some(b"IIII".as_slice()));
    }

    #[test]
    fn test_quality_is_optional() {
        let record = Record::with_attrs("read1", None, b"ACGT", None);

        assert_eq!(record.qual(), None);
        assert_eq!(record.seq(), b"ACGT");
    }

    #[test]
    fn test_new_record_is_empty() {
        let record = Record::new();
        assert_eq!(record.id(), "");
        assert_eq!(record.desc(), None);
        assert_eq!(record.seq(), b"");
        assert_eq!(record.qual(), None);
    }

    #[test]
    fn test_set_from_fastq_reuses_buffers() {
        let src1 = fastq::Record::with_attrs("read1", Some("first"), b"ACGT", b"IIII");
        let src2 = fastq::Record::with_attrs("read2", Some("second"), b"ACGTACGT", b"IIIIIIII");

        let mut record = Record::new();
        record.set_from_fastq(&src1);
        let seq_ptr = record.seq().as_ptr();
        assert_eq!(record.id(), "read1");
        assert_eq!(record.desc(), Some("first"));
        assert_eq!(record.qual(), Some(b"IIII".as_slice()));

        record.set_from_fastq(&src2);
        assert_eq!(record.id(), "read2");
        assert_eq!(record.desc(), Some("second"));
        assert_eq!(record.seq(), b"ACGTACGT");
        assert_eq!(record.qual(), Some(b"IIIIIIII".as_slice()));
        assert_eq!(
            seq_ptr,
            record.seq().as_ptr(),
            "sequence buffer should be reused"
        );
    }

    #[test]
    fn test_set_from_fastq_drops_absent_description() {
        let with_desc = fastq::Record::with_attrs("read1", Some("note"), b"ACGT", b"IIII");
        let without_desc = fastq::Record::with_attrs("read2", None, b"ACGT", b"IIII");

        let mut record = Record::new();
        record.set_from_fastq(&with_desc);
        assert_eq!(record.desc(), Some("note"));

        record.set_from_fastq(&without_desc);
        assert_eq!(record.desc(), None);
    }
}
