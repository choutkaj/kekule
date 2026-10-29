use crate::TrajectoryField;
use std::fmt;
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Seek, Write};
use std::path::{Path, PathBuf};
use std::sync::{
    atomic::{AtomicU64, Ordering},
    Arc,
};

use kekule::diagnostics::{UnrepresentablePolicy, WriteReport};
use kekule::topology::Topology;

use crate::{
    TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
    TrajectoryWriter,
};
use kekule::structure::{Trajectory, TrajectoryFrameView};

use super::{codec_context, dcd, detect, io_context, trr, xtc, xyz, FrameTimingPolicy};

static TEMPORARY_FILE_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// Existing-destination policy for a path writer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum OverwritePolicy {
    /// Never replace a file. Publication hard-links the completed temporary
    /// file to the destination, which fails atomically if the destination
    /// exists; no portable rename refuses to replace. Hard-link support is
    /// checked when the writer is created, so a filesystem without it (FAT,
    /// exFAT, some network shares) is rejected before any frame is written;
    /// use [`Self::Replace`] there.
    Forbid,
    /// Atomically rename the completed temporary file over any existing
    /// destination.
    Replace,
}

/// Format-agnostic writer configuration.
///
/// [`Default`] takes the format from the destination extension (`.xyz`,
/// `.dcd`, `.trr`, or `.xtc`, case-insensitively); [`Self::new`] names it
/// explicitly, as a writer over a host's stream requires. Both forbid
/// overwriting and use each codec's default options.
#[derive(Debug, Clone)]
pub struct TrajectoryWriteOptions {
    format: Option<TrajectoryFormat>,
    overwrite: OverwritePolicy,
    xyz: xyz::XyzWriteOptions,
    dcd: dcd::DcdWriteOptions,
    trr: trr::TrrWriteOptions,
    xtc: xtc::XtcWriteOptions,
}

impl Default for TrajectoryWriteOptions {
    fn default() -> Self {
        Self {
            format: None,
            overwrite: OverwritePolicy::Forbid,
            xyz: xyz::XyzWriteOptions::default(),
            dcd: dcd::DcdWriteOptions::default(),
            trr: trr::TrrWriteOptions::default(),
            xtc: xtc::XtcWriteOptions::default(),
        }
    }
}

impl TrajectoryWriteOptions {
    /// Options with an explicit format, which takes precedence over the
    /// destination extension.
    pub fn new(format: TrajectoryFormat) -> Self {
        Self {
            format: Some(format),
            overwrite: OverwritePolicy::Forbid,
            xyz: xyz::XyzWriteOptions::default(),
            dcd: dcd::DcdWriteOptions::default(),
            trr: trr::TrrWriteOptions::default(),
            xtc: xtc::XtcWriteOptions::default(),
        }
    }

    /// Sets [`Self::overwrite_policy`].
    #[must_use]
    pub fn with_overwrite_policy(mut self, overwrite: OverwritePolicy) -> Self {
        self.overwrite = overwrite;
        self
    }

    /// Sets [`Self::xyz_options`].
    #[must_use]
    pub fn with_xyz_options(mut self, options: xyz::XyzWriteOptions) -> Self {
        self.xyz = options;
        self
    }

    /// Sets [`Self::dcd_options`].
    #[must_use]
    pub fn with_dcd_options(mut self, options: dcd::DcdWriteOptions) -> Self {
        self.dcd = options;
        self
    }

    /// Sets [`Self::trr_options`].
    #[must_use]
    pub fn with_trr_options(mut self, options: trr::TrrWriteOptions) -> Self {
        self.trr = options;
        self
    }

    /// Sets [`Self::xtc_options`].
    #[must_use]
    pub fn with_xtc_options(mut self, options: xtc::XtcWriteOptions) -> Self {
        self.xtc = options;
        self
    }

    /// Sets the `frame_metadata` policy of the codec options whose formats
    /// cannot carry some frame velocities, forces, times, or steps, so that
    /// whichever is written leaves them out, or fails on them:
    /// [`xyz::XyzWriteOptions::frame_metadata`],
    /// [`dcd::DcdWriteOptions::frame_metadata`], and
    /// [`xtc::XtcWriteOptions::frame_metadata`]; TRR carries them all. The
    /// writer's report names what was left out. A later `with_*_options`
    /// call replaces the policy this sets for its codec.
    #[must_use]
    pub fn with_frame_metadata(mut self, frame_metadata: UnrepresentablePolicy) -> Self {
        self.xyz = self.xyz.with_frame_metadata(frame_metadata);
        self.dcd = self.dcd.with_frame_metadata(frame_metadata);
        self.xtc = self.xtc.with_frame_metadata(frame_metadata);
        self
    }

    /// Sets the `frame_properties` policy of every codec's options, for the
    /// atom, bond, and frame properties none of the formats carry:
    /// [`xyz::XyzWriteOptions::frame_properties`],
    /// [`dcd::DcdWriteOptions::frame_properties`],
    /// [`trr::TrrWriteOptions::frame_properties`], and
    /// [`xtc::XtcWriteOptions::frame_properties`]. A later `with_*_options`
    /// call replaces the policy this sets for its codec.
    #[must_use]
    pub fn with_frame_properties(mut self, frame_properties: UnrepresentablePolicy) -> Self {
        self.xyz = self.xyz.with_frame_properties(frame_properties);
        self.dcd = self.dcd.with_frame_properties(frame_properties);
        self.trr = self.trr.with_frame_properties(frame_properties);
        self.xtc = self.xtc.with_frame_properties(frame_properties);
        self
    }

    /// Sets the `frame_timing` policy of the TRR and XTC options, the
    /// formats that record a time and a step for every frame:
    /// [`trr::TrrWriteOptions::frame_timing`] and
    /// [`xtc::XtcWriteOptions::frame_timing`]. A later `with_*_options` call
    /// replaces the policy this sets for its codec. DCD numbers frames on its
    /// header's step sequence instead, under its own
    /// [`dcd::DcdWriteOptions::frame_timing`].
    #[must_use]
    pub fn with_frame_timing(mut self, frame_timing: FrameTimingPolicy) -> Self {
        self.trr = self.trr.with_frame_timing(frame_timing);
        self.xtc = self.xtc.with_frame_timing(frame_timing);
        self
    }

    /// The explicit format, or `None` to take it from the destination extension.
    pub const fn format(&self) -> Option<TrajectoryFormat> {
        self.format
    }

    /// What to do when the destination exists; [`OverwritePolicy::Forbid`] by
    /// default. Writers over a host's stream ignore it.
    pub const fn overwrite_policy(&self) -> OverwritePolicy {
        self.overwrite
    }

    /// Options for the XYZ writer, used when that format is selected.
    pub const fn xyz_options(&self) -> &xyz::XyzWriteOptions {
        &self.xyz
    }

    /// Options for the DCD writer, used when that format is selected.
    pub const fn dcd_options(&self) -> &dcd::DcdWriteOptions {
        &self.dcd
    }

    /// Options for the TRR writer, used when that format is selected.
    pub const fn trr_options(&self) -> &trr::TrrWriteOptions {
        &self.trr
    }

    /// Options for the XTC writer, used when that format is selected.
    pub const fn xtc_options(&self) -> &xtc::XtcWriteOptions {
        &self.xtc
    }
}

/// The codec writer of the selected format.
enum CodecWriter<W> {
    Xyz(xyz::XyzWriter<W>),
    Dcd(dcd::DcdWriter<W>),
    Trr(trr::TrrWriter<W>),
    Xtc(xtc::XtcWriter<W>),
}

impl<W: Write + Seek> CodecWriter<W> {
    fn new(
        format: TrajectoryFormat,
        writer: W,
        topology: Arc<Topology>,
        options: TrajectoryWriteOptions,
        label: String,
    ) -> Result<Self, TrajectoryIoError> {
        match format {
            TrajectoryFormat::Xyz => {
                xyz::XyzWriter::new(writer, topology, options.xyz, label).map(Self::Xyz)
            }
            TrajectoryFormat::Dcd => {
                dcd::DcdWriter::new(writer, topology, options.dcd, label).map(Self::Dcd)
            }
            TrajectoryFormat::Trr => {
                trr::TrrWriter::new(writer, topology, options.trr, label).map(Self::Trr)
            }
            TrajectoryFormat::Xtc => {
                xtc::XtcWriter::new(writer, topology, options.xtc, label).map(Self::Xtc)
            }
        }
    }

    fn report(&self) -> &WriteReport {
        match self {
            Self::Xyz(writer) => writer.report(),
            Self::Dcd(writer) => writer.report(),
            Self::Trr(writer) => writer.report(),
            Self::Xtc(writer) => writer.report(),
        }
    }

    fn writer(&self) -> &W {
        match self {
            Self::Xyz(writer) => writer.writer(),
            Self::Dcd(writer) => writer.writer(),
            Self::Trr(writer) => writer.writer(),
            Self::Xtc(writer) => writer.writer(),
        }
    }

    fn finish(self) -> Result<W, TrajectoryIoError> {
        match self {
            Self::Xyz(writer) => writer.finish(),
            Self::Dcd(writer) => writer.finish(),
            Self::Trr(writer) => writer.finish(),
            Self::Xtc(writer) => writer.finish(),
        }
    }

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

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryIoError> {
        match self {
            Self::Xyz(writer) => writer.write_frame(frame),
            Self::Dcd(writer) => writer.write_frame(frame),
            Self::Trr(writer) => writer.write_frame(frame),
            Self::Xtc(writer) => writer.write_frame(frame),
        }
    }

    fn flush(&mut self) -> Result<(), TrajectoryIoError> {
        match self {
            Self::Xyz(writer) => writer.flush(),
            Self::Dcd(writer) => writer.flush(),
            Self::Trr(writer) => writer.flush(),
            Self::Xtc(writer) => writer.flush(),
        }
    }
}

/// Strict writer of one format over a host's stream, from
/// [`create_trajectory_stream_writer`].
///
/// The stream holds a complete trajectory only after successful
/// [`Self::finish`]. A failed frame write poisons the writer, as it does a
/// [`FileTrajectoryWriter`], so that a host never takes a trajectory missing
/// a frame for a complete one; the codec writers of the [`xyz`], [`dcd`],
/// [`trr`], and [`xtc`] modules leave the stream as it was when they reject
/// a frame and accept the next.
pub struct TrajectoryStreamWriter<W> {
    inner: CodecWriter<W>,
    format: TrajectoryFormat,
    source_label: String,
    failed: bool,
}

impl<W> fmt::Debug for TrajectoryStreamWriter<W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("TrajectoryStreamWriter")
            .field("format", &self.format)
            .field("source_label", &self.source_label)
            .field("failed", &self.failed)
            .finish_non_exhaustive()
    }
}

impl<W: Write + Seek> TrajectoryStreamWriter<W> {
    /// Finalizes format metadata, such as the DCD frame count, flushes, and
    /// returns the stream, positioned at the end of the trajectory. The
    /// host synchronizes and places the file itself.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidFrame`] after an earlier failed
    /// frame write, when no frame was written, or when the stream ends in a
    /// partial frame; the codec's finalization errors; and
    /// [`TrajectoryIoErrorKind::Io`] when finalizing or flushing fails.
    pub fn finish(self) -> Result<W, TrajectoryIoError> {
        if self.failed {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::Finish,
                Some(self.format),
                &self.source_label,
                "cannot finish a trajectory after an earlier frame-write failure",
            ));
        }
        self.inner.finish()
    }

    /// The format being written.
    pub const fn format(&self) -> TrajectoryFormat {
        self.format
    }

    /// The stream being written.
    pub fn writer(&self) -> &W {
        self.inner.writer()
    }

    /// The content the frames written so far leave out: for XYZ,
    /// [`xyz::XyzWriter::report`]; for DCD, TRR, and XTC, the frame state
    /// their `frame_metadata` options omit ([`dcd::DcdWriter::report`],
    /// [`trr::TrrWriter::report`], [`xtc::XtcWriter::report`]).
    pub fn report(&self) -> &WriteReport {
        self.inner.report()
    }

    fn poisoned(&self, operation: TrajectoryIoOperation) -> TrajectoryIoError {
        codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            operation,
            Some(self.format),
            &self.source_label,
            "trajectory writer is poisoned by an earlier frame-write failure",
        )
    }
}

impl<W: Write + Seek> TrajectoryWriter for TrajectoryStreamWriter<W> {
    fn topology(&self) -> &Topology {
        self.inner.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.inner.shared_topology()
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryIoError> {
        if self.failed {
            return Err(self.poisoned(TrajectoryIoOperation::WriteFrame));
        }
        let result = self.inner.write_frame(frame);
        self.failed = result.is_err();
        result
    }

    /// Flushes the codec writer, so the host's stream holds every frame
    /// written so far, whole; a DCD header then counts them. Use it to let
    /// another process read a trajectory while it is written.
    fn flush(&mut self) -> Result<(), TrajectoryIoError> {
        if self.failed {
            return Err(self.poisoned(TrajectoryIoOperation::WriteFrame));
        }
        self.inner.flush()
    }
}

/// Strict atomic path writer.
///
/// A nonempty file is published only by successful [`Self::finish`].
pub struct FileTrajectoryWriter {
    inner: Option<TrajectoryStreamWriter<BufWriter<File>>>,
    format: TrajectoryFormat,
    destination: PathBuf,
    temporary: PathBuf,
    overwrite: OverwritePolicy,
    published: bool,
}

impl fmt::Debug for FileTrajectoryWriter {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FileTrajectoryWriter")
            .field("format", &self.format)
            .field("destination", &self.destination)
            .field("overwrite", &self.overwrite)
            .field("published", &self.published)
            .finish_non_exhaustive()
    }
}

impl FileTrajectoryWriter {
    /// Flushes, synchronizes, and atomically publishes a nonempty trajectory,
    /// then synchronizes the destination directory on Unix so the publication
    /// survives a crash. An error from that last step means the trajectory is
    /// published but its directory entry may not yet be durable.
    ///
    /// Finishing before any successful frame write returns a structured error
    /// and removes the unpublished temporary sibling.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidFrame`] after an earlier failed
    /// frame write, when no frame was written, or when the stream ends in a
    /// partial frame; the codec's finalization errors; and
    /// [`TrajectoryIoErrorKind::Io`] when flushing, synchronizing or
    /// publishing fails, including when [`OverwritePolicy::Forbid`] finds the
    /// destination created meanwhile.
    pub fn finish(mut self) -> Result<(), TrajectoryIoError> {
        let label = self.destination.display().to_string();
        if self.stream().failed {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::Finish,
                Some(self.format),
                &label,
                "cannot publish a trajectory after an earlier frame-write failure",
            ));
        }
        if let Some(inner) = self.inner.take() {
            let file = inner.finish()?;
            file.get_ref().sync_all().map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Finish,
                    Some(self.format),
                    &label,
                    error,
                )
            })?;
        }
        let finish_error = |error| {
            io_context(
                TrajectoryIoOperation::Finish,
                Some(self.format),
                &label,
                error,
            )
        };
        match self.overwrite {
            OverwritePolicy::Forbid => {
                std::fs::hard_link(&self.temporary, &self.destination).map(|()| {
                    let _ = std::fs::remove_file(&self.temporary);
                })
            }
            OverwritePolicy::Replace => std::fs::rename(&self.temporary, &self.destination),
        }
        .map_err(finish_error)?;
        self.published = true;
        sync_parent_directory(&self.destination).map_err(finish_error)
    }

    /// The format being written.
    pub const fn format(&self) -> TrajectoryFormat {
        self.format
    }

    /// The content the frames written so far leave out, as
    /// [`TrajectoryStreamWriter::report`] gives it.
    pub fn report(&self) -> &WriteReport {
        self.stream().report()
    }

    fn stream(&self) -> &TrajectoryStreamWriter<BufWriter<File>> {
        match &self.inner {
            Some(inner) => inner,
            None => unreachable!("finished path writers are consumed"),
        }
    }
}

impl TrajectoryWriter for FileTrajectoryWriter {
    fn topology(&self) -> &Topology {
        self.stream().topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.stream().shared_topology()
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryIoError> {
        match &mut self.inner {
            Some(inner) => inner.write_frame(frame),
            None => unreachable!("finished path writers are consumed"),
        }
    }

    /// Flushes the frames written so far to the temporary sibling; the
    /// destination appears only when the writer finishes.
    fn flush(&mut self) -> Result<(), TrajectoryIoError> {
        match &mut self.inner {
            Some(inner) => inner.flush(),
            None => unreachable!("finished path writers are consumed"),
        }
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

/// Writes a loaded trajectory.
///
/// The codec is the format of `options` or, by default, the one its destination
/// extension names: `.xyz`, `.dcd`, `.trr`, or `.xtc`, case-insensitively.
/// Unknown extensions require an explicit format. `options` also sets the
/// precision, field, and overwrite policies; by default existing files are
/// rejected. Unsupported metadata is rejected; no fields are silently discarded. Collection
/// properties are unsupported by these codecs and must be explicitly cleared first.
/// The destination is published atomically only after every frame and finalization
/// succeeds. Empty trajectories are rejected by the writer.
///
/// Returns the [`FileTrajectoryWriter::report`] of the content the file
/// leaves out, such as the bonds of an XYZ file or the cells and frame state
/// its options omit.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::UnknownFormat`] when no format is given and
/// the extension names none, [`TrajectoryIoErrorKind::UnsupportedField`] for
/// collection properties or a frame field the format cannot store, and the
/// errors of [`create_trajectory_writer`], the codec's
/// [`TrajectoryWriter::write_frame`] and [`FileTrajectoryWriter::finish`]. No
/// file is published on error.
///
/// ```no_run
/// # use kekule::structure::Trajectory;
/// # fn save(trajectory: &Trajectory) -> Result<(), kekule_traj::TrajectoryIoError> {
/// let report = kekule_traj::io::write_trajectory("aligned.xyz", trajectory, Default::default())?;
/// for omission in report.omissions() {
///     eprintln!("left out {omission:?}");
/// }
/// # Ok(()) }
/// ```
pub fn write_trajectory(
    path: impl AsRef<Path>,
    trajectory: &Trajectory,
    options: TrajectoryWriteOptions,
) -> Result<WriteReport, TrajectoryIoError> {
    let path = path.as_ref();
    resolve_format(path, options.format)?;
    if !trajectory.properties().is_empty() {
        return Err(TrajectoryIoError::unsupported_field(
            TrajectoryField::CollectionProperties,
        ));
    }
    let mut writer = create_trajectory_writer(path, trajectory.shared_topology(), options)?;
    for frame in trajectory {
        writer.write_frame(frame)?;
    }
    let report = writer.report().clone();
    writer.finish()?;
    Ok(report)
}

/// Writes a loaded trajectory to a host's stream, as [`write_trajectory`]
/// writes a file, so that the host places the file itself, for example
/// under a temporary name of its own that it then moves over the
/// destination.
///
/// `options` must name the format ([`TrajectoryWriteOptions::new`]); its
/// overwrite policy does not apply. The trajectory is written from the
/// stream's current position and the stream is left at its end. Wrap an
/// unbuffered sink such as a `File` in a [`BufWriter`]. Returns the
/// [`TrajectoryStreamWriter::report`] of the content the stream leaves out.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::UnknownFormat`] when `options` names no
/// format, [`TrajectoryIoErrorKind::UnsupportedField`] for collection
/// properties or a frame field the format cannot store, and the errors of
/// [`create_trajectory_stream_writer`], the codec's
/// [`TrajectoryWriter::write_frame`] and [`TrajectoryStreamWriter::finish`].
/// Output before a failure has already reached `writer`.
///
/// ```
/// # use std::io::Cursor;
/// # use std::sync::Arc;
/// # use kekule::geometry::Point3;
/// # use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
/// # use kekule::units::{Quantity, ANGSTROM, PICOSECOND};
/// # use kekule::{smiles, topology::Topology};
/// use kekule::diagnostics::{StructureField, UnrepresentablePolicy, WriteOmission};
/// use kekule_traj::io::{write_trajectory_to, TrajectoryWriteOptions};
/// use kekule_traj::TrajectoryFormat;
/// # let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
/// # let topology = Arc::new(Topology::from_molecule(molecule)?);
/// # let points = vec![Point3::new(1.5, 0.0, 0.0)];
/// # let mut frame = TrajectoryFrame::new(Positions::from_vec(Quantity::new(points, ANGSTROM))?);
/// # frame.set_step(Some(0));
/// # frame.set_time(Some(Quantity::new(0.0, PICOSECOND)))?;
/// # let trajectory = Trajectory::from_items(topology, [frame])?;
///
/// // A frame with a time, as read from XTC, written as DCD without it.
/// let mut stream = Cursor::new(Vec::new());
/// let options = TrajectoryWriteOptions::new(TrajectoryFormat::Dcd)
///     .with_frame_metadata(UnrepresentablePolicy::Omit);
/// let report = write_trajectory_to(&mut stream, &trajectory, options)?;
/// let omitted = WriteOmission::StructureState { field: StructureField::Time, structures: 1 };
/// assert_eq!(report.omissions(), [omitted]);
/// assert!(!stream.into_inner().is_empty());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn write_trajectory_to<W: Write + Seek + ?Sized>(
    writer: &mut W,
    trajectory: &Trajectory,
    options: TrajectoryWriteOptions,
) -> Result<WriteReport, TrajectoryIoError> {
    if !trajectory.properties().is_empty() {
        return Err(TrajectoryIoError::unsupported_field(
            TrajectoryField::CollectionProperties,
        ));
    }
    let mut writer = create_trajectory_stream_writer(
        writer,
        trajectory.shared_topology(),
        options,
        "trajectory stream",
    )?;
    for frame in trajectory {
        writer.write_frame(frame)?;
    }
    let report = writer.report().clone();
    writer.finish()?;
    Ok(report)
}

/// Creates a strict writer of the format `options` names over a host's
/// stream; `source_label` names the destination in errors.
///
/// The writer writes from the stream's current position: the DCD header
/// there, patched when the writer finishes, and the other formats' frames
/// as they are written. `options` must name the format
/// ([`TrajectoryWriteOptions::new`]); its overwrite policy does not apply.
/// Wrap an unbuffered sink such as a `File` in a [`BufWriter`].
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::UnknownFormat`] when `options` names no
/// format, and the option and stream errors of the selected codec writer's
/// constructor.
pub fn create_trajectory_stream_writer<W: Write + Seek>(
    writer: W,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryWriteOptions,
    source_label: impl Into<String>,
) -> Result<TrajectoryStreamWriter<W>, TrajectoryIoError> {
    let source_label = source_label.into();
    let format = options.format.ok_or_else(|| {
        codec_context(
            TrajectoryIoErrorKind::UnknownFormat,
            TrajectoryIoOperation::Open,
            None,
            &source_label,
            "a stream writer needs an explicit format; use TrajectoryWriteOptions::new",
        )
    })?;
    stream_writer(format, writer, topology.into(), options, source_label)
}

fn stream_writer<W: Write + Seek>(
    format: TrajectoryFormat,
    writer: W,
    topology: Arc<Topology>,
    options: TrajectoryWriteOptions,
    source_label: String,
) -> Result<TrajectoryStreamWriter<W>, TrajectoryIoError> {
    Ok(TrajectoryStreamWriter {
        inner: CodecWriter::new(format, writer, topology, options, source_label.clone())?,
        format,
        source_label,
        failed: false,
    })
}

/// The explicit format, or the one the destination extension names.
fn resolve_format(
    path: &Path,
    format: Option<TrajectoryFormat>,
) -> Result<TrajectoryFormat, TrajectoryIoError> {
    match format {
        Some(format) => Ok(format),
        None => detect::extension_format(path).ok_or_else(|| codec_context(
            TrajectoryIoErrorKind::UnknownFormat, TrajectoryIoOperation::Open, None,
            &path.display().to_string(), "destination extension does not identify a supported format; supply explicit write options",
        )),
    }
}

/// Creates a strict path writer backed by a temporary sibling file.
///
/// The codec is the format of `options` or, by default, the one the
/// destination extension names, as for [`write_trajectory`].
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::UnknownFormat`] when no format is given and
/// the extension names none, [`TrajectoryIoErrorKind::Io`] with
/// [`std::io::ErrorKind::AlreadyExists`] when the destination exists under
/// [`OverwritePolicy::Forbid`] or another I/O error when the temporary file
/// cannot be created, [`TrajectoryIoErrorKind::UnsupportedVariant`] when
/// [`OverwritePolicy::Forbid`] is used where hard links are unavailable, and
/// the option errors of the selected codec writer's constructor.
pub fn create_trajectory_writer(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryWriteOptions,
) -> Result<FileTrajectoryWriter, TrajectoryIoError> {
    let topology = topology.into();
    let destination = path.as_ref().to_path_buf();
    let format = resolve_format(&destination, options.format)?;
    let label = destination.display().to_string();
    if options.overwrite == OverwritePolicy::Forbid && destination.exists() {
        return Err(io_context(
            TrajectoryIoOperation::Open,
            Some(format),
            &label,
            io::Error::new(io::ErrorKind::AlreadyExists, "destination already exists"),
        ));
    }
    let temporary = create_temporary_sibling(&destination, format)?;
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .map_err(|error| io_context(TrajectoryIoOperation::Open, Some(format), &label, error))?;
    if options.overwrite == OverwritePolicy::Forbid {
        if let Err(error) = require_hard_links(&temporary, |from, to| std::fs::hard_link(from, to))
        {
            let _ = std::fs::remove_file(&temporary);
            return Err(codec_context(
                TrajectoryIoErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::Open,
                Some(format),
                &label,
                format!(
                    "OverwritePolicy::Forbid publishes through a hard link, which this destination cannot use ({error}); OverwritePolicy::Replace publishes by rename"
                ),
            ));
        }
    }
    let overwrite = options.overwrite;
    let inner = match stream_writer(format, BufWriter::new(file), topology, options, label) {
        Ok(inner) => inner,
        Err(error) => {
            let _ = std::fs::remove_file(&temporary);
            return Err(error);
        }
    };
    Ok(FileTrajectoryWriter {
        inner: Some(inner),
        format,
        destination,
        temporary,
        overwrite,
        published: false,
    })
}

fn create_temporary_sibling(
    destination: &Path,
    format: TrajectoryFormat,
) -> Result<PathBuf, TrajectoryIoError> {
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
        TrajectoryIoErrorKind::LimitExceeded,
        TrajectoryIoOperation::Open,
        Some(format),
        &destination.display().to_string(),
        "could not reserve a unique temporary sibling name",
    ))
}

/// Checks that the temporary file can be hard-linked beside itself, as
/// [`OverwritePolicy::Forbid`] publication will need. `link` is the hard-link
/// primitive, injectable so an unsupported filesystem can be tested.
fn require_hard_links(
    temporary: &Path,
    link: impl FnOnce(&Path, &Path) -> io::Result<()>,
) -> io::Result<()> {
    let probe = temporary.with_extension("link-probe");
    link(temporary, &probe)?;
    let _ = std::fs::remove_file(&probe);
    Ok(())
}

/// Makes a publication durable: a rename or link survives a crash only once
/// its directory is synchronized. Filesystems that cannot synchronize a
/// directory report `InvalidInput` or `Unsupported`, which is accepted.
#[cfg(unix)]
fn sync_parent_directory(destination: &Path) -> io::Result<()> {
    let parent = destination
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    match File::open(parent).and_then(|directory| directory.sync_all()) {
        Err(error)
            if matches!(
                error.kind(),
                io::ErrorKind::InvalidInput | io::ErrorKind::Unsupported
            ) =>
        {
            Ok(())
        }
        result => result,
    }
}

/// The standard library cannot open a directory handle on this platform, so
/// directory synchronization is left to the filesystem.
#[cfg(not(unix))]
fn sync_parent_directory(_destination: &Path) -> io::Result<()> {
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    #[test]
    fn forbid_requires_hard_links_and_leaves_no_probe_behind() {
        let directory = std::env::temp_dir().join(format!(
            "kekule-publish-{}-{}",
            std::process::id(),
            TEMPORARY_FILE_SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&directory).unwrap();
        let temporary = directory.join(".out.xyz.1.2.tmp");
        fs::write(&temporary, b"trajectory").unwrap();

        require_hard_links(&temporary, |from, to| fs::hard_link(from, to)).unwrap();
        let error = require_hard_links(&temporary, |_, _| {
            Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "filesystem has no hard links",
            ))
        })
        .unwrap_err();
        assert_eq!(
            (error.kind(), error.to_string()),
            (
                io::ErrorKind::Unsupported,
                "filesystem has no hard links".to_owned()
            )
        );
        let mut names = fs::read_dir(&directory)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect::<Vec<_>>();
        names.sort();
        assert_eq!(names, [".out.xyz.1.2.tmp"]);
        fs::remove_dir_all(directory).unwrap();
    }
}
