mod adapter;
mod app;
mod buffer;
mod color;
mod reader;
mod viewer;
use app::run;
use std::env;

fn main() {
    let path = env::args().nth(1).unwrap_or_else(|| "-".to_string());
    if let Err(e) = run(&path) {
        eprintln!("Application error: {:?}", e);
        std::process::exit(1);
    }
}
