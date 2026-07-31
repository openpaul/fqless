use crate::viewer::TuiViewer;

// run needs to support stdin or a file path (or two for paired-end)
pub fn run(path: &str, path2: Option<&str>) -> anyhow::Result<()> {
    let mut viewer = TuiViewer::new(String::from(path), path2.map(String::from))?;
    viewer.run()?;
    Ok(())
}
