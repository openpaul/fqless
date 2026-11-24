//! A background thread that can be stopped, is joined when dropped, and keeps
//! its error so the UI can show it (printing from a worker would land on the
//! alternate screen and be wiped by the next frame).

use anyhow::Result;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

#[derive(Default)]
pub struct Worker {
    handle: Option<JoinHandle<()>>,
    stop: Arc<AtomicBool>,
    error: Arc<Mutex<Option<String>>>,
}

impl Worker {
    /// Run `job` on a new thread. The job should return early once the stop
    /// flag it is given becomes true.
    pub fn spawn<F>(job: F) -> Self
    where
        F: FnOnce(&AtomicBool) -> Result<()> + Send + 'static,
    {
        let stop = Arc::new(AtomicBool::new(false));
        let error = Arc::new(Mutex::new(None));
        let handle = {
            let stop = Arc::clone(&stop);
            let error = Arc::clone(&error);
            thread::spawn(move || {
                if let Err(e) = job(&stop) {
                    *error.lock().unwrap() = Some(format!("{e:#}"));
                }
            })
        };
        Self {
            handle: Some(handle),
            stop,
            error,
        }
    }

    /// Ask the job to stop and wait for it.
    pub fn stop(&mut self) {
        if let Some(handle) = self.handle.take() {
            self.stop.store(true, Ordering::SeqCst);
            let _ = handle.join();
        }
    }

    pub fn is_running(&self) -> bool {
        self.handle.as_ref().is_some_and(|h| !h.is_finished())
    }

    /// The error the job returned, if it failed.
    pub fn error(&self) -> Option<String> {
        self.error.lock().unwrap().clone()
    }
}

impl Drop for Worker {
    fn drop(&mut self) {
        self.stop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn test_stop_ends_a_running_job() {
        let mut worker = Worker::spawn(|stop| {
            while !stop.load(Ordering::Relaxed) {
                thread::sleep(Duration::from_millis(1));
            }
            Ok(())
        });
        assert!(worker.is_running());
        worker.stop();
        assert!(!worker.is_running());
        assert!(worker.error().is_none());
    }

    #[test]
    fn test_error_is_kept() {
        let mut worker = Worker::spawn(|_| anyhow::bail!("broken input"));
        while worker.is_running() {
            thread::sleep(Duration::from_millis(1));
        }
        worker.stop();
        assert_eq!(worker.error().as_deref(), Some("broken input"));
    }

    #[test]
    fn test_default_worker_is_idle() {
        let mut worker = Worker::default();
        assert!(!worker.is_running());
        worker.stop();
    }
}
