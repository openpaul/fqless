mod adapter;
mod app;
mod buffer;
mod color;
mod reader;
mod record;
mod search;
mod utils;
mod viewer;
use app::run;
use std::env;

fn print_usage() {
    println!(
        "fqless v{} - FastQ/FASTA/SAM/BAM File Viewer\n",
        env!("CARGO_PKG_VERSION")
    );
    println!("Usage: fqless <FASTQ|FASTA|SAM|BAM file> [<R2 file>]\n");
    println!("Options:");
    println!("  -h, --help     Show this help message");
    println!("\nThe format (FASTQ/FASTA/SAM/BAM, plain or gzipped) is detected automatically.");
    println!("If the input is '-', reads from stdin.");
    println!();
    println!("Examples:");
    println!("  fqless myfile.fastq.gz                      # View a gzipped FASTQ file");
    println!("  fqless myfile.fastq                         # View a regular FASTQ file");
    println!("  fqless myfile.fasta                         # View a FASTA file (no quality)");
    println!("  fqless myfile.sam                           # View reads from a SAM file");
    println!("  fqless myfile.bam                           # View reads from a BAM file");
    println!("  fqless R1.fastq.gz R2.fastq.gz              # View paired-end FASTQ files");
    println!("  cat myfile.fastq.gz | fqless -              # View gzipped FASTQ from stdin");
    println!("  cat myfile.fastq | fqless -                 # View regular FASTQ from stdin");
    println!();
    println!("Note: Statistics are limited to the loaded reads when using stdin.");
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        print_usage();
        return;
    }
    let path = &args[1];
    let path2 = if args.len() > 2 && args[2] != "-h" && args[2] != "--help" {
        Some(args[2].clone())
    } else {
        None
    };
    if let Err(e) = run(path, path2.as_deref()) {
        eprintln!("Application error: {:?}", e);
        std::process::exit(1);
    }
}
