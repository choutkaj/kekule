//! Publish complete report snapshots without destroying the previous snapshot.
use std::{
    error::Error,
    fs, io,
    path::{Path, PathBuf},
    process,
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};

pub(super) struct ReportFile {
    path: PathBuf,
    published: bool,
}

impl ReportFile {
    pub(super) fn new(path: &Path) -> Self {
        Self {
            path: path.to_owned(),
            published: false,
        }
    }

    pub(super) fn write(
        &mut self,
        write: impl FnOnce(&mut fs::File) -> Result<(), Box<dyn Error>>,
    ) -> Result<(), Box<dyn Error>> {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let temporary = self.path.with_extension(format!(
            "{}-{}-{}.tmp",
            process::id(),
            SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos(),
            NEXT.fetch_add(1, Ordering::Relaxed),
        ));
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|error| report_io_error("create staging file for", &self.path, error))?;
        let result = (|| {
            write(&mut file)?;
            file.sync_all()
                .map_err(|error| report_io_error("sync staging file for", &self.path, error))?;
            drop(file);
            let publish = || {
                if self.published {
                    fs::rename(&temporary, &self.path)
                } else {
                    // First publication must not replace another run's report.
                    fs::hard_link(&temporary, &self.path)
                }
            };
            #[cfg(windows)]
            let publication = retry_publication(publish);
            #[cfg(not(windows))]
            let publication = publish();
            publication.map_err(|error| {
                report_io_error(
                    if self.published { "replace" } else { "publish" },
                    &self.path,
                    error,
                )
            })?;
            Ok(())
        })();
        let _ = fs::remove_file(&temporary);
        if result.is_ok() {
            self.published = true;
        }
        result
    }
}

fn report_io_error(operation: &str, path: &Path, error: io::Error) -> io::Error {
    io::Error::new(
        error.kind(),
        format!("cannot {operation} report '{}': {error}", path.display()),
    )
}

#[cfg(windows)]
fn retry_publication(mut publish: impl FnMut() -> io::Result<()>) -> io::Result<()> {
    use std::time::{Duration, Instant};
    let start = Instant::now();
    let timeout = Duration::from_secs(2);
    let mut delay = Duration::from_millis(10);
    loop {
        match publish() {
            // A reader without FILE_SHARE_DELETE can make rename return access
            // denied (5), sharing violation (32), or lock violation (33). Retry
            // only publication, never serialization or benchmark computation.
            Err(error)
                if matches!(error.raw_os_error(), Some(5 | 32 | 33))
                    && start.elapsed() < timeout =>
            {
                std::thread::sleep(delay.min(timeout.saturating_sub(start.elapsed())));
                delay = (delay * 2).min(Duration::from_millis(100));
            }
            result => return result,
        }
    }
}
