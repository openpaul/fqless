use crate::viewer::TuiViewer;

// run needs to support stdin or a file path
pub fn run(path: &str) -> anyhow::Result<()> {
    let mut viewer = TuiViewer::new(String::from(path))?;
    viewer.run()?;
    Ok(())
}
