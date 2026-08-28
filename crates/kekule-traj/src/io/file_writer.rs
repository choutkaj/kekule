use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::{
    atomic::{AtomicU64, Ordering},
    Arc,
};

use kekule::topology::Topology;

use crate::{
    TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat, TrajectoryFrameView,
    TrajectoryIoOperation, TrajectoryWriter,
};

use super::{codec_context, dcd, io_context, trr, xtc, xyz};

static TEMPORARY_FILE_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// Existing-destination policy for a path writer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum OverwritePolicy {
    Forbid,
    Replace,
}

/// Format-agnostic path-writer configuration.
#[derive(Debug, Clone)]
pub struct TrajectoryWriteOptions {
    format: TrajectoryFormat,
    overwrite: OverwritePolicy,
    xyz: xyz::XyzWriteOptions,
    dcd: dcd::DcdWriteOptions,
    trr: trr::TrrWriteOptions,
    xtc: xtc::XtcWriteOptions,
}

impl TrajectoryWriteOptions {
    pub fn new(format: TrajectoryFormat) -> Self {
        Self {
            format,
            overwrite: OverwritePolicy::Forbid,
            xyz: xyz::XyzWriteOptions::default(),
            dcd: dcd::DcdWriteOptions::default(),
            trr: trr::TrrWriteOptions::default(),
            xtc: xtc::XtcWriteOptions::default(),
        }
    }

    pub fn with_overwrite_policy(mut self, overwrite: OverwritePolicy) -> Self {
        self.overwrite = overwrite;
        self
    }

    pub fn with_xyz_options(mut self, options: xyz::XyzWriteOptions) -> Self {
        self.xyz = options;
        self
    }

    pub fn with_dcd_options(mut self, options: dcd::DcdWriteOptions) -> Self {
        self.dcd = options;
        self
    }

    pub fn with_trr_options(mut self, options: trr::TrrWriteOptions) -> Self {
        self.trr = options;
        self
    }

    pub fn with_xtc_options(mut self, options: xtc::XtcWriteOptions) -> Self {
        self.xtc = options;
        self
    }
}

enum FileWriterInner {
    Xyz(xyz::XyzWriter<BufWriter<File>>),
    Dcd(dcd::DcdWriter<BufWriter<File>>),
    Trr(trr::TrrWriter<BufWriter<File>>),
    Xtc(xtc::XtcWriter<BufWriter<File>>),
}

impl FileWriterInner {
    fn topology(&self) -> &Topology {
        match self {
            Self::Xyz(writer) => writer.topology(),
            Self::Dcd(writer) => writer.topology(),
            Self::Trr(writer) => writer.topology(),
            Self::Xtc(writer) => writer.topology(),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match self {
            Self::Xyz(writer) => writer.shared_topology(),
            Self::Dcd(writer) => writer.shared_topology(),
            Self::Trr(writer) => writer.shared_topology(),
            Self::Xtc(writer) => writer.shared_topology(),
        }
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        match self {
            Self::Xyz(writer) => writer.write_frame(frame),
            Self::Dcd(writer) => writer.write_frame(frame),
            Self::Trr(writer) => writer.write_frame(frame),
            Self::Xtc(writer) => writer.write_frame(frame),
        }
    }

    fn flush_and_sync(&mut self, label: &str) -> Result<(), TrajectoryError> {
        match self {
            Self::Xyz(writer) => {
                writer.validate_finish()?;
                writer.flush().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Xyz),
                        label,
                        error,
                    )
                })?;
                writer.writer().get_ref().sync_all().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Xyz),
                        label,
                        error,
                    )
                })
            }
            Self::Dcd(writer) => {
                writer.finalize()?;
                writer.flush().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Dcd),
                        label,
                        error,
                    )
                })?;
                writer.writer().get_ref().sync_all().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Dcd),
                        label,
                        error,
                    )
                })
            }
            Self::Trr(writer) => {
                writer.validate_finish()?;
                writer.flush().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Trr),
                        label,
                        error,
                    )
                })?;
                writer.writer().get_ref().sync_all().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Trr),
                        label,
                        error,
                    )
                })
            }
            Self::Xtc(writer) => {
                writer.validate_finish()?;
                writer.flush().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Xtc),
                        label,
                        error,
                    )
                })?;
                writer.writer().get_ref().sync_all().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Xtc),
                        label,
                        error,
                    )
                })
            }
        }
    }
}

/// Strict atomic path writer.
///
/// A nonempty file is published only by successful [`Self::finish`].
pub struct FileTrajectoryWriter {
    inner: Option<FileWriterInner>,
    format: TrajectoryFormat,
    destination: PathBuf,
    temporary: PathBuf,
    overwrite: OverwritePolicy,
    failed: bool,
    published: bool,
}

impl FileTrajectoryWriter {
    /// Flushes, synchronizes, and atomically publishes a nonempty trajectory.
    ///
    /// Finishing before any successful frame write returns a structured error
    /// and removes the unpublished temporary sibling.
    pub fn finish(mut self) -> Result<(), TrajectoryError> {
        let label = self.destination.display().to_string();
        if self.failed {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::Finish,
                Some(self.format()),
                &label,
                "cannot publish a trajectory after an earlier frame-write failure",
            ));
        }
        if let Some(inner) = &mut self.inner {
            inner.flush_and_sync(&label)?;
        }
        self.inner.take();
        match self.overwrite {
            OverwritePolicy::Forbid => {
                std::fs::hard_link(&self.temporary, &self.destination).map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(self.format()),
                        &label,
                        error,
                    )
                })?;
                self.published = true;
                let _ = std::fs::remove_file(&self.temporary);
            }
            OverwritePolicy::Replace => {
                std::fs::rename(&self.temporary, &self.destination).map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(self.format()),
                        &label,
                        error,
                    )
                })?;
                self.published = true;
            }
        }
        Ok(())
    }

    pub const fn format(&self) -> TrajectoryFormat {
        self.format
    }
}

impl TrajectoryWriter for FileTrajectoryWriter {
    fn topology(&self) -> &Topology {
        match &self.inner {
            Some(inner) => inner.topology(),
            None => unreachable!("finished path writers are consumed"),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match &self.inner {
            Some(inner) => inner.shared_topology(),
            None => unreachable!("finished path writers are consumed"),
        }
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if self.failed {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::WriteFrame,
                Some(self.format),
                &self.destination.display().to_string(),
                "trajectory writer is poisoned by an earlier frame-write failure",
            ));
        }
        let result = match &mut self.inner {
            Some(inner) => inner.write_frame(frame),
            None => Err(codec_context(
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::WriteFrame,
                Some(self.format),
                &self.destination.display().to_string(),
                "trajectory writer has already finished",
            )),
        };
        if result.is_err() {
            self.failed = true;
        }
        result
    }
}

impl Drop for FileTrajectoryWriter {
    fn drop(&mut self) {
        if !self.published {
            self.inner.take();
            let _ = std::fs::remove_file(&self.temporary);
        }
    }
}

/// Creates a strict path writer backed by a temporary sibling file.
pub fn create_trajectory_writer(
    path: impl AsRef<Path>,
    topology: Arc<Topology>,
    options: TrajectoryWriteOptions,
) -> Result<FileTrajectoryWriter, TrajectoryError> {
    let destination = path.as_ref().to_path_buf();
    let label = destination.display().to_string();
    if options.overwrite == OverwritePolicy::Forbid && destination.exists() {
        return Err(io_context(
            TrajectoryIoOperation::Open,
            Some(options.format),
            &label,
            io::Error::new(io::ErrorKind::AlreadyExists, "destination already exists"),
        ));
    }
    #[cfg(windows)]
    if options.overwrite == OverwritePolicy::Replace && destination.exists() {
        return Err(codec_context(
            TrajectoryCodecErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::Open,
            Some(options.format),
            &label,
            "atomic replacement of an existing destination is unavailable on this platform",
        ));
    }
    let temporary = create_temporary_sibling(&destination, options.format)?;
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(options.format),
                &label,
                error,
            )
        })?;
    let inner_result = match options.format {
        TrajectoryFormat::Xyz => {
            xyz::XyzWriter::new(BufWriter::new(file), topology, options.xyz, label)
                .map(FileWriterInner::Xyz)
        }
        TrajectoryFormat::Dcd => {
            dcd::DcdWriter::new(BufWriter::new(file), topology, options.dcd, label)
                .map(FileWriterInner::Dcd)
        }
        TrajectoryFormat::Trr => {
            trr::TrrWriter::new(BufWriter::new(file), topology, options.trr, label)
                .map(FileWriterInner::Trr)
        }
        TrajectoryFormat::Xtc => {
            xtc::XtcWriter::new(BufWriter::new(file), topology, options.xtc, label)
                .map(FileWriterInner::Xtc)
        }
    };
    let inner = match inner_result {
        Ok(inner) => inner,
        Err(error) => {
            let _ = std::fs::remove_file(&temporary);
            return Err(error);
        }
    };
    Ok(FileTrajectoryWriter {
        inner: Some(inner),
        format: options.format,
        destination,
        temporary,
        overwrite: options.overwrite,
        failed: false,
        published: false,
    })
}

fn create_temporary_sibling(
    destination: &Path,
    format: TrajectoryFormat,
) -> Result<PathBuf, TrajectoryError> {
    let parent = destination.parent().unwrap_or_else(|| Path::new("."));
    let name = destination
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("trajectory");
    for _ in 0..128 {
        let sequence = TEMPORARY_FILE_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let candidate = parent.join(format!(".{name}.{}.{}.tmp", std::process::id(), sequence));
        if !candidate.exists() {
            return Ok(candidate);
        }
    }
    Err(codec_context(
        TrajectoryCodecErrorKind::ResourceLimitExceeded,
        TrajectoryIoOperation::Open,
        Some(format),
        &destination.display().to_string(),
        "could not reserve a unique temporary sibling name",
    ))
}
