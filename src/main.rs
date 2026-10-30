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
        "fqless v{} - FastQ/FASTA/SAM/BAM/CRAM File Viewer\n",
        env!("CARGO_PKG_VERSION")
    );
    println!("Usage: fqless [OPTIONS] <FASTQ|FASTA|SAM|BAM|CRAM file> [<R2 file>]\n");
    println!("Options:");
    println!("  -h, --help                Show this help message");
    println!("  -r, --reference <FASTA>   Reference FASTA, required to decode CRAM input");
    println!(
        "\nThe format (FASTQ/FASTA/SAM/BAM/CRAM, plain or gzipped) is detected automatically."
    );
    println!("If the input is '-', reads from stdin.");
    println!();
    println!("Examples:");
    println!("  fqless myfile.fastq.gz                      # View a gzipped FASTQ file");
    println!("  fqless myfile.fastq                         # View a regular FASTQ file");
    println!("  fqless myfile.fasta                         # View a FASTA file (no quality)");
    println!("  fqless myfile.sam                           # View reads from a SAM file");
    println!("  fqless myfile.bam                           # View reads from a BAM file");
    println!("  fqless -r reference.fa myfile.cram          # View a CRAM file");
    println!("  fqless R1.fastq.gz R2.fastq.gz              # View paired-end FASTQ files");
    println!("  cat myfile.fastq.gz | fqless -              # View gzipped FASTQ from stdin");
    println!("  cat myfile.fastq | fqless -                 # View regular FASTQ from stdin");
    println!();
    println!("Note: Statistics are limited to the loaded reads when using stdin.");
}

/// Parsed command line: positional input file(s) plus the CRAM reference.
struct Args {
    path: String,
    path2: Option<String>,
    reference: Option<String>,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut paths: Vec<String> = Vec::new();
    let mut reference = None;

    let mut iter = argv.iter().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print_usage();
                std::process::exit(0);
            }
            "-r" | "--reference" => {
                let value = iter
                    .next()
                    .ok_or_else(|| format!("{arg} requires a FASTA path"))?;
                reference = Some(value.clone());
            }
            other if other.starts_with("--reference=") => {
                reference = Some(other["--reference=".len()..].to_string());
            }
            other if other.starts_with('-') && other != "-" => {
                return Err(format!("unknown option: {other}"));
            }
            other => paths.push(other.to_string()),
        }
    }

    let mut paths = paths.into_iter();
    let path = paths
        .next()
        .ok_or_else(|| "missing input file".to_string())?;
    let path2 = paths.next();
    if paths.next().is_some() {
        return Err("too many input files (expected at most two: R1 and R2)".to_string());
    }

    Ok(Args {
        path,
        path2,
        reference,
    })
}

fn main() {
    let argv: Vec<String> = env::args().collect();
    if argv.len() < 2 {
        print_usage();
        return;
    }

    let args = match parse_args(&argv) {
        Ok(args) => args,
        Err(message) => {
            eprintln!("fqless: {message}\n");
            print_usage();
            std::process::exit(2);
        }
    };

    if let Err(e) = run(&args.path, args.path2.as_deref(), args.reference.as_deref()) {
        eprintln!("Application error: {:?}", e);
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Result<Args, String> {
        let argv: Vec<String> = std::iter::once("fqless")
            .chain(args.iter().copied())
            .map(str::to_string)
            .collect();
        parse_args(&argv)
    }

    #[test]
    fn test_parse_plain_input() {
        let args = parse(&["reads.fastq"]).unwrap();
        assert_eq!(args.path, "reads.fastq");
        assert!(args.path2.is_none());
        assert!(args.reference.is_none());
    }

    #[test]
    fn test_parse_paired_input() {
        let args = parse(&["R1.fastq", "R2.fastq"]).unwrap();
        assert_eq!(args.path, "R1.fastq");
        assert_eq!(args.path2.as_deref(), Some("R2.fastq"));
    }

    #[test]
    fn test_parse_reference_forms() {
        let args = parse(&["-r", "ref.fa", "reads.cram"]).unwrap();
        assert_eq!(args.reference.as_deref(), Some("ref.fa"));
        assert_eq!(args.path, "reads.cram");

        let args = parse(&["reads.cram", "--reference", "ref.fa"]).unwrap();
        assert_eq!(args.reference.as_deref(), Some("ref.fa"));

        let args = parse(&["--reference=ref.fa", "reads.cram"]).unwrap();
        assert_eq!(args.reference.as_deref(), Some("ref.fa"));
    }

    #[test]
    fn test_parse_stdin_is_not_an_option() {
        let args = parse(&["-"]).unwrap();
        assert_eq!(args.path, "-");
    }

    #[test]
    fn test_parse_rejects_bad_input() {
        assert!(parse(&[]).is_err());
        assert!(parse(&["a", "b", "c"]).is_err());
        assert!(parse(&["--index", "x"]).is_err());
        assert!(parse(&["--reference"]).is_err());
    }
}
