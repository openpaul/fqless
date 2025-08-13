mod adapter;
mod app;
mod buffer;
mod color;
mod reader;
mod viewer;
use app::run;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        print_usage();
        return;
    }
    let path = &args[1];
    if let Err(e) = run(path) {
        eprintln!("Application error: {:?}", e);
        std::process::exit(1);
    }

    fn print_usage() {
        println!("Usage: fqless <FASTQ file>\n");
        println!("Options:");
        println!("  -h, --help     Show this help message");
        println!("\nIf <FASTQ file> is '-', reads from stdin.");
        println!("");
        println!("Examples:");
        println!("  fqless myfile.fastq.gz           # View a gzipped FASTQ file");
        println!("  fqless myfile.fastq              # View a regular FASTQ file");
        println!("  cat myfile.fastq.gz | fqless -   # View gzipped FASTQ from stdin");
        println!("  cat myfile.fastq | fqless -      # View regular FASTQ from stdin");
        println!("");
        println!("Note: Statistics are limited to the loaded reads when using stdin.");
    }
}
