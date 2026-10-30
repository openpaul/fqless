use crate::viewer::TuiViewer;

// run needs to support stdin or a file path (or two for paired-end), plus the
// reference FASTA used for CRAM input.
pub fn run(path: &str, path2: Option<&str>, reference: Option<&str>) -> anyhow::Result<()> {
    let mut viewer = TuiViewer::new(
        String::from(path),
        path2.map(String::from),
        reference.map(String::from),
    )?;
    viewer.run()?;
    Ok(())
}
