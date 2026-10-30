//! Generate the CRAM examples from the SAM examples and the fake reference.
//!
//! CRAM stores sequences as differences against a reference, so unlike the
//! other formats it cannot be produced by a small script. This example reads a
//! SAM file and a reference FASTA with noodles and writes a CRAM file plus its
//! `.crai` index.
//!
//! Run from the repository root:
//!
//! ```sh
//! cargo run --example make_cram -- \
//!     example/simple.sam example/fake_reference.fasta example/simple.cram
//! ```

use std::{
    env,
    error::Error,
    fs::File,
    io::{self, BufReader, BufWriter},
};

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, alignment::io::Write as _};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args().skip(1);
    let sam_path = args.next().ok_or_else(usage)?;
    let reference_path = args.next().ok_or_else(usage)?;
    let cram_path = args.next().ok_or_else(usage)?;

    let mut sam_reader = sam::io::Reader::new(BufReader::new(File::open(&sam_path)?));
    let header = sam_reader.read_header()?;
    let records = sam_reader.records().collect::<io::Result<Vec<_>>>()?;

    let mut reference_reader = fasta::io::Reader::new(BufReader::new(File::open(&reference_path)?));
    let references = reference_reader.records().collect::<io::Result<Vec<_>>>()?;
    let repository = fasta::Repository::new(references);

    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_writer(BufWriter::new(File::create(&cram_path)?));

    writer.write_header(&header)?;
    for record in &records {
        writer.write_alignment_record(&header, record)?;
    }
    writer.try_finish(&header)?;
    // try_finish writes through the writer; make sure the file is on disk
    // before indexing it.
    std::io::Write::flush(writer.get_mut())?;

    let index = cram::fs::index(&cram_path)?;
    cram::crai::fs::write(format!("{cram_path}.crai"), &index)?;

    eprintln!(
        "wrote {cram_path} ({} records) and {cram_path}.crai",
        records.len()
    );
    Ok(())
}

fn usage() -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        "usage: make_cram <input.sam> <reference.fa> <output.cram>",
    )
}
