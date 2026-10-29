use std::fmt;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use kekule::topology::Topology;

use crate::{
    FrameBuffer, MemoryTrajectoryWriter, SeekableTrajectoryReader, TrajectoryFormat,
    TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation, TrajectoryReader,
    TrajectoryWriter,
};
use kekule::diagnostics::{LimitExceeded, Resource};
use kekule::structure::Trajectory;

use super::limits::index_hard_capacity;
use super::{
    dcd, detect, index_mismatch, io_context, trr, xtc, xyz, FileTrajectoryMetadata, FrameOffsets,
    IndexBuildControl, TrajectoryFileFingerprint, TrajectoryFrameIndex, TrajectoryOpenOptions,
    TrajectoryOpenReport,
};

enum SequentialReaderInner {
    Xyz(xyz::XyzReader<BufReader<File>>),
    Dcd(dcd::DcdReader<BufReader<File>>),
    Trr(Box<trr::TrrReader<BufReader<File>>>),
    Xtc(xtc::XtcReader<BufReader<File>>),
}

/// Format-agnostic path-backed sequential reader retaining one file handle.
pub struct SequentialFileTrajectoryReader {
    inner: SequentialReaderInner,
    metadata: FileTrajectoryMetadata,
    open_report: TrajectoryOpenReport,
}

impl fmt::Debug for SequentialFileTrajectoryReader {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SequentialFileTrajectoryReader")
            .field("metadata", &self.metadata)
            .field("open_report", &self.open_report)
            .finish_non_exhaustive()
    }
}

impl SequentialFileTrajectoryReader {
    /// Format detection and non-fatal facts recorded when this reader was opened.
    pub fn open_report(&self) -> &TrajectoryOpenReport {
        &self.open_report
    }

    /// Returns metadata verified through the most recent successful read.
    ///
    /// In particular, mixed-width TRR input is reported as mixed as soon as
    /// the second scalar width has been observed.
    pub fn metadata(&self) -> &FileTrajectoryMetadata {
        &self.metadata
    }
}

impl TrajectoryReader for SequentialFileTrajectoryReader {
    fn topology(&self) -> &Topology {
        match &self.inner {
            SequentialReaderInner::Xyz(reader) => reader.topology(),
            SequentialReaderInner::Dcd(reader) => reader.topology(),
            SequentialReaderInner::Trr(reader) => reader.topology(),
            SequentialReaderInner::Xtc(reader) => reader.topology(),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match &self.inner {
            SequentialReaderInner::Xyz(reader) => reader.shared_topology(),
            SequentialReaderInner::Dcd(reader) => reader.shared_topology(),
            SequentialReaderInner::Trr(reader) => reader.shared_topology(),
            SequentialReaderInner::Xtc(reader) => reader.shared_topology(),
        }
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        match &mut self.inner {
            SequentialReaderInner::Xyz(reader) => reader.read_next(destination),
            SequentialReaderInner::Dcd(reader) => reader.read_next(destination),
            SequentialReaderInner::Trr(reader) => {
                let read = reader.read_next(destination)?;
                self.metadata
                    .update_trr_precision(reader.first_header(), reader.precision_mixed());
                Ok(read)
            }
            SequentialReaderInner::Xtc(reader) => reader.read_next(destination),
        }
    }
}

enum IndexedReaderInner {
    Xyz(xyz::IndexedXyzReader<BufReader<File>>),
    Dcd(dcd::IndexedDcdReader<BufReader<File>>),
    Trr(Box<trr::IndexedTrrReader<BufReader<File>>>),
    Xtc(xtc::IndexedXtcReader<BufReader<File>>),
}

/// What an index extension may have to undo.
#[derive(Debug, Clone, Copy)]
struct IndexSnapshot {
    mark: (usize, u64),
    pending_bytes: u64,
    precision_mixed: bool,
    declared_frames: Option<usize>,
}

impl IndexedReaderInner {
    fn offsets(&self) -> &FrameOffsets {
        match self {
            Self::Xyz(reader) => reader.offsets(),
            Self::Dcd(reader) => reader.offsets(),
            Self::Trr(reader) => reader.offsets(),
            Self::Xtc(reader) => reader.offsets(),
        }
    }

    fn pending_bytes(&self) -> u64 {
        match self {
            Self::Xyz(reader) => reader.pending_bytes(),
            Self::Dcd(reader) => reader.pending_bytes(),
            Self::Trr(reader) => reader.pending_bytes(),
            Self::Xtc(reader) => reader.pending_bytes(),
        }
    }

    fn precision_mixed(&self) -> bool {
        match self {
            Self::Trr(reader) => reader.precision_mixed(),
            _ => false,
        }
    }

    fn declared_frames(&self) -> Option<usize> {
        match self {
            Self::Dcd(reader) => Some(reader.header().declared_frames()),
            _ => None,
        }
    }

    fn snapshot(&self) -> IndexSnapshot {
        IndexSnapshot {
            mark: self.offsets().mark(),
            pending_bytes: self.pending_bytes(),
            precision_mixed: self.precision_mixed(),
            declared_frames: self.declared_frames(),
        }
    }

    fn restore(&mut self, snapshot: IndexSnapshot) {
        let IndexSnapshot {
            mark,
            pending_bytes,
            precision_mixed,
            declared_frames,
        } = snapshot;
        match self {
            Self::Xyz(reader) => reader.restore_index(mark, pending_bytes),
            Self::Dcd(reader) => {
                reader.restore_index(mark, declared_frames.unwrap_or(0), pending_bytes);
            }
            Self::Trr(reader) => reader.restore_index(mark, precision_mixed, pending_bytes),
            Self::Xtc(reader) => reader.restore_index(mark, pending_bytes),
        }
    }

    /// Indexes the complete frames after the indexed ones; see the codecs'
    /// `extend`.
    fn extend(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        match self {
            Self::Xyz(reader) => reader.extend(growing, control),
            Self::Dcd(reader) => reader.extend(growing, control),
            Self::Trr(reader) => reader.extend(growing, control),
            Self::Xtc(reader) => reader.extend(growing, control),
        }
    }
}

/// Format-agnostic path-backed indexed reader retaining one file handle.
///
/// It indexes the complete frames of the file and verifies each before it
/// reads it. A file that keeps growing after opening, such as the
/// trajectory of a running simulation, is followed with [`Self::refresh`],
/// which indexes the frames appended since without scanning the indexed
/// ones again; open such a file with
/// [`TrajectoryOpenOptions::with_growing`], so a frame still being written
/// is not yet an error. Sequential reads end after the last indexed frame.
pub struct IndexedFileTrajectoryReader {
    inner: IndexedReaderInner,
    metadata: FileTrajectoryMetadata,
    open_report: TrajectoryOpenReport,
    /// The fingerprint of the indexed bytes: the file up to the end of its
    /// last indexed frame.
    fingerprint: TrajectoryFileFingerprint,
    path: PathBuf,
    growing: bool,
}

impl fmt::Debug for IndexedFileTrajectoryReader {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedFileTrajectoryReader")
            .field("metadata", &self.metadata)
            .field("open_report", &self.open_report)
            .finish_non_exhaustive()
    }
}

impl IndexedFileTrajectoryReader {
    /// Format detection and non-fatal facts recorded when this reader was opened.
    pub fn open_report(&self) -> &TrajectoryOpenReport {
        &self.open_report
    }

    /// Metadata verified for the indexed frames while building the index,
    /// or recorded in the imported index whose bytes the file begins with,
    /// and for the frames indexed since.
    pub fn metadata(&self) -> &FileTrajectoryMetadata {
        &self.metadata
    }

    /// Exports the frame offsets, times, and steps with the fingerprint of
    /// the indexed bytes, copying them, for
    /// [`open_indexed_trajectory_from_index`].
    pub fn to_frame_index(&self) -> TrajectoryFrameIndex {
        TrajectoryFrameIndex::new(
            self.metadata.format(),
            self.fingerprint.clone(),
            self.inner.offsets().clone(),
            self.inner.precision_mixed(),
        )
    }

    /// The bytes the file held after the last indexed frame when the reader
    /// last opened or [refreshed](Self::refresh): a frame still being
    /// written, or more that arrived during the scan. Zero when the indexed
    /// frames reached the end of the file. Only a reader opened with
    /// [`TrajectoryOpenOptions::with_growing`] leaves a partial frame
    /// pending; for others it is an error.
    pub fn pending_byte_count(&self) -> u64 {
        self.inner.pending_bytes()
    }

    /// Indexes the frames appended to the file since the reader was opened
    /// or last refreshed, and returns how many there are.
    ///
    /// The file must still begin with the indexed bytes, as a file that only
    /// grew by appending does; the frame-count fields of a DCD header, which
    /// writers update as they append, are not compared (see
    /// [`TrajectoryFileFingerprint`]). The indexed frames are not read
    /// again: the scan starts after the last of them, and verifies every new
    /// frame as opening does. Under [`TrajectoryOpenOptions::with_growing`] a
    /// final frame still being written stays pending
    /// ([`Self::pending_byte_count`]) until a later refresh, and DCD frames
    /// are counted in the file rather than checked against the header's
    /// `NSET`; otherwise a partial frame is an error, and a DCD must hold
    /// exactly the frames its header now counts. Sequential reads resume
    /// where they stopped and continue into the new frames, and
    /// [`Self::to_frame_index`] then exports them too.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::IndexMismatch`] when the indexed bytes
    /// changed, including while the new frames were read, the decoding error
    /// of the first invalid new frame,
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when the index outgrows the
    /// limits, and [`TrajectoryIoErrorKind::Io`] when the file cannot be
    /// read. Failure changes nothing: the reader keeps its frames and can be
    /// refreshed again.
    ///
    /// ```no_run
    /// use kekule::{smiles, topology::Topology};
    /// use kekule_traj::io::{open_indexed_trajectory, TrajectoryOpenOptions};
    /// use kekule_traj::{SeekableTrajectoryReader, TrajectoryReader};
    ///
    /// let molecule = smiles::to_molecules("CC")?.pop().unwrap();
    /// let topology = Topology::from_molecule(molecule)?;
    /// // The trajectory of a running simulation.
    /// let options = TrajectoryOpenOptions::default().with_growing(true);
    /// let mut reader = open_indexed_trajectory("run.xtc", topology, options)?;
    /// let mut frame = reader.frame_buffer();
    /// loop {
    ///     let added = reader.refresh()?;
    ///     let count = reader.frame_count().unwrap_or(0);
    ///     for index in count - added..count {
    ///         reader.read_frame(index, &mut frame)?;
    ///         // Show the new frame.
    ///     }
    ///     # break;
    ///     // Wait for the simulation to write more.
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// [`TrajectoryIoErrorKind::IndexMismatch`]: crate::TrajectoryIoErrorKind::IndexMismatch
    /// [`TrajectoryIoErrorKind::LimitExceeded`]: crate::TrajectoryIoErrorKind::LimitExceeded
    /// [`TrajectoryIoErrorKind::Io`]: crate::TrajectoryIoErrorKind::Io
    pub fn refresh(&mut self) -> Result<usize, TrajectoryIoError> {
        let label = self.path.display().to_string();
        let format = self.metadata.format();
        let mut file = File::open(&self.path).map_err(|error| {
            io_context(TrajectoryIoOperation::Index, Some(format), &label, error)
        })?;
        let changed = || {
            index_mismatch(
                TrajectoryIoOperation::Index,
                format,
                &label,
                "the file's indexed bytes changed since they were indexed",
            )
        };
        if !self
            .fingerprint
            .is_prefix_of_labelled_file(&mut file, &label)?
        {
            return Err(changed());
        }
        let length = file
            .metadata()
            .map_err(|error| io_context(TrajectoryIoOperation::Index, Some(format), &label, error))?
            .len();
        let indexed = self.inner.offsets().end();
        if self.growing && length == indexed + self.inner.pending_bytes() {
            // Nothing arrived; a pending frame is as incomplete as before.
            return Ok(0);
        }
        let before = self.inner.offsets().len();
        let snapshot = self.inner.snapshot();
        self.inner
            .extend(self.growing, &IndexBuildControl::default())?;
        let end = self.inner.offsets().end();
        let checked = self
            .fingerprint
            .is_prefix_of_labelled_file(&mut file, &label)
            .and_then(|unchanged| {
                if unchanged {
                    TrajectoryFileFingerprint::of_labelled_file_prefix(&mut file, end, &label)
                } else {
                    Err(changed())
                }
            });
        match checked {
            Ok(fingerprint) => {
                self.fingerprint = fingerprint;
                self.record_counts();
                Ok(self.inner.offsets().len() - before)
            }
            Err(error) => {
                self.inner.restore(snapshot);
                Err(error)
            }
        }
    }

    /// Updates the metadata with the indexed frames.
    fn record_counts(&mut self) {
        let count = self.inner.offsets().len();
        self.metadata
            .update_frame_counts(count, self.inner.declared_frames());
        if let IndexedReaderInner::Trr(reader) = &self.inner {
            self.metadata
                .update_trr_precision(reader.first_header(), reader.precision_mixed());
        }
    }
}

impl TrajectoryReader for IndexedFileTrajectoryReader {
    fn topology(&self) -> &Topology {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.topology(),
            IndexedReaderInner::Dcd(reader) => reader.topology(),
            IndexedReaderInner::Trr(reader) => reader.topology(),
            IndexedReaderInner::Xtc(reader) => reader.topology(),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.shared_topology(),
            IndexedReaderInner::Dcd(reader) => reader.shared_topology(),
            IndexedReaderInner::Trr(reader) => reader.shared_topology(),
            IndexedReaderInner::Xtc(reader) => reader.shared_topology(),
        }
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        match &mut self.inner {
            IndexedReaderInner::Xyz(reader) => reader.read_next(destination),
            IndexedReaderInner::Dcd(reader) => reader.read_next(destination),
            IndexedReaderInner::Trr(reader) => reader.read_next(destination),
            IndexedReaderInner::Xtc(reader) => reader.read_next(destination),
        }
    }
}

impl SeekableTrajectoryReader for IndexedFileTrajectoryReader {
    fn frame_count(&self) -> Option<usize> {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.frame_count(),
            IndexedReaderInner::Dcd(reader) => reader.frame_count(),
            IndexedReaderInner::Trr(reader) => reader.frame_count(),
            IndexedReaderInner::Xtc(reader) => reader.frame_count(),
        }
    }

    fn read_frame(
        &mut self,
        index: usize,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        match &mut self.inner {
            IndexedReaderInner::Xyz(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Dcd(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Trr(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Xtc(reader) => reader.read_frame(index, destination),
        }
    }
}

/// Reads an entire trajectory file into memory in file frame order.
///
/// Accepts an owned topology or a shared `Arc<Topology>` and detects the format
/// automatically. File coordinates must follow the topology's dense atom order;
/// counts and available format metadata are checked, but matching counts alone
/// cannot establish atom identity. All decoded frame fields are retained, and
/// a supplied shared topology retains its exact allocation.
///
/// This consumes the sequential reader through clean EOF. Any open or decoding
/// error is returned without publishing a partial trajectory. Memory usage grows
/// with the number of frames and atoms; use [`open_trajectory`] to process large
/// files one frame at a time or inspect file metadata and opening diagnostics.
///
/// `options` sets the format policies and limits; frame limits cause an error
/// rather than truncating the result. These codec limits do not impose a total
/// resident-memory budget for the loaded trajectory. No random-access index is
/// built.
///
/// # Errors
///
/// Returns the errors of [`open_trajectory`] and the first decoding,
/// validation or [`TrajectoryIoErrorKind::Io`] error of any frame.
///
/// [`TrajectoryIoErrorKind::Io`]: crate::TrajectoryIoErrorKind::Io
///
/// ```no_run
/// use kekule::mmcif;
/// use kekule_traj::io::read_trajectory;
///
/// let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?, Default::default())?;
/// let topology = document.interpret(Default::default())?.into_topology();
/// let trajectory = read_trajectory("trajectory.xtc", topology, Default::default())?;
/// println!("{} frames, {} atoms", trajectory.len(), trajectory.topology().atom_count());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn read_trajectory(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryOpenOptions,
) -> Result<Trajectory, TrajectoryIoError> {
    let mut reader = open_trajectory(path, topology, options)?;
    let mut buffer = reader.frame_buffer();
    let mut writer = MemoryTrajectoryWriter::new(reader.shared_topology());
    while reader.read_next(&mut buffer)? {
        writer.write_frame(buffer.frame_view())?;
    }
    Ok(writer.into_trajectory())
}

/// Opens one fast sequential path-backed trajectory reader without loading all frames.
///
/// Accepts an owned topology or a shared `Arc<Topology>`. File coordinates must
/// follow its dense atom order; counts and available format metadata are checked
/// automatically. Matching counts alone cannot establish atom identity.
/// `options` sets the format policies and limits. For a fully loaded in-memory
/// trajectory, use [`read_trajectory`] instead.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::Io`] when the file cannot be opened or
/// read, the detection errors of [`super::detect_trajectory_format`], and the
/// header errors of the selected codec's reader, such as
/// [`TrajectoryIoErrorKind::InconsistentAtomCount`] when the file's atom count
/// differs from the topology's.
///
/// [`TrajectoryIoErrorKind::Io`]: crate::TrajectoryIoErrorKind::Io
/// [`TrajectoryIoErrorKind::InconsistentAtomCount`]:
///     crate::TrajectoryIoErrorKind::InconsistentAtomCount
///
/// ```no_run
/// use kekule::{smiles, topology::Topology};
/// use kekule_traj::{io::open_trajectory, TrajectoryReader};
///
/// let molecule = smiles::to_molecules("CC")?.pop().unwrap();
/// let topology = Topology::from_molecule(molecule)?;
/// let mut reader = open_trajectory("ethane.xyz", topology, Default::default())?;
/// let mut frame = reader.frame_buffer();
/// while reader.read_next(&mut frame)? {
///     let model = frame.as_model_view();
///     assert_eq!(model.atom_count(), 2);
/// }
/// println!("{:?}", reader.open_report().selected_format());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn open_trajectory(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryOpenOptions,
) -> Result<SequentialFileTrajectoryReader, TrajectoryIoError> {
    let topology = topology.into();
    let path = path.as_ref();
    let label = path.display().to_string();
    let file = File::open(path)
        .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &label, error))?;
    let mut reader = BufReader::new(file);
    let detection = detect::select_format(
        &mut reader,
        path,
        options.format_hint,
        &options.limits,
        &label,
        false,
    )?;
    let (inner, metadata, notes) = match detection.format {
        TrajectoryFormat::Xyz => {
            let mut reader = xyz::XyzReader::new(
                reader,
                topology,
                options
                    .xyz
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            reader.validate_first_frame()?;
            let atom_count = reader.topology().atom_count();
            (
                SequentialReaderInner::Xyz(reader),
                FileTrajectoryMetadata::xyz(atom_count, None),
                vec!["XYZ coordinates use the configured explicit/default length unit".into()],
            )
        }
        TrajectoryFormat::Dcd => {
            let reader = dcd::DcdReader::new(
                reader,
                topology,
                options
                    .dcd
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let metadata = FileTrajectoryMetadata::dcd(reader.header(), None);
            (
                SequentialReaderInner::Dcd(reader),
                metadata,
                vec![
                    "DCD steps are preserved from ISTART/NSAVC; time follows the explicit policy"
                        .into(),
                ],
            )
        }
        TrajectoryFormat::Trr => {
            let reader = trr::TrrReader::new(
                reader,
                topology,
                options
                    .trr
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let metadata = FileTrajectoryMetadata::trr(reader.first_header(), None, false);
            (
                SequentialReaderInner::Trr(Box::new(reader)),
                metadata,
                vec![
                    "TRR uses XDR big-endian scalars and explicit lambda preservation policy"
                        .into(),
                ],
            )
        }
        TrajectoryFormat::Xtc => {
            let cell_policy = options.xtc.cell_policy();
            let xtc_options = options.xtc;
            let reader = xtc::XtcReader::new(
                reader,
                topology,
                xtc_options
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let metadata = FileTrajectoryMetadata::xtc(reader.first_info(), cell_policy, None);
            (
                SequentialReaderInner::Xtc(reader),
                metadata,
                vec!["XTC decoding uses a bounded checked reader".into()],
            )
        }
    };
    let report = TrajectoryOpenReport {
        selected_format: detection.format,
        detection_evidence: detection.evidence,
        notes,
    };
    Ok(SequentialFileTrajectoryReader {
        inner,
        metadata,
        open_report: report,
    })
}

/// Opens one fully verified indexed path-backed trajectory reader.
///
/// Uses the topology order contract of [`open_trajectory`]: file coordinate
/// index `i` is topology dense atom index `i`, and readers check counts and
/// available metadata, including XYZ element order. Opening scans and validates
/// the entire file and stores bounded frame offsets for random access.
/// `options` sets the format policies and limits, and may set a cancellation
/// and a progress callback for the scan
/// ([`TrajectoryOpenOptions::with_index_cancellation`],
/// [`TrajectoryOpenOptions::with_index_progress`]). A file still being
/// written opens with [`TrajectoryOpenOptions::with_growing`]: its complete
/// frames are indexed, a final frame still being written is
/// [pending](IndexedFileTrajectoryReader::pending_byte_count), and
/// [`IndexedFileTrajectoryReader::refresh`] indexes the frames written
/// later.
///
/// The bytes the file holds when it is opened must not change during the
/// scan, though more may be appended. The reader fingerprints the indexed
/// bytes, so [`IndexedFileTrajectoryReader::to_frame_index`] can export the
/// offsets for [`open_indexed_trajectory_from_index`] to reopen the file
/// without scanning it again.
///
/// # Errors
///
/// Returns the errors of [`open_trajectory`], the decoding error of the first
/// invalid frame, [`TrajectoryIoErrorKind::LimitExceeded`] when the index
/// outgrows the limits, [`TrajectoryIoErrorKind::Cancelled`] when the
/// cancellation stops the scan, [`TrajectoryIoErrorKind::IndexMismatch`]
/// when bytes of the file change while it is scanned, and, for a growing
/// file, [`TrajectoryIoErrorKind::Incomplete`] when it does not yet hold a
/// complete frame, or for DCD a complete header.
///
/// [`TrajectoryIoErrorKind::LimitExceeded`]: crate::TrajectoryIoErrorKind::LimitExceeded
/// [`TrajectoryIoErrorKind::Cancelled`]: crate::TrajectoryIoErrorKind::Cancelled
/// [`TrajectoryIoErrorKind::IndexMismatch`]: crate::TrajectoryIoErrorKind::IndexMismatch
/// [`TrajectoryIoErrorKind::Incomplete`]: crate::TrajectoryIoErrorKind::Incomplete
///
/// ```no_run
/// use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
/// use std::sync::Arc;
///
/// use kekule::diagnostics::Cancellation;
/// use kekule::{smiles, topology::Topology};
/// use kekule_traj::io::{open_indexed_trajectory, TrajectoryOpenOptions};
///
/// let molecule = smiles::to_molecules("CC")?.pop().unwrap();
/// let topology = Topology::from_molecule(molecule)?;
/// let cancel = Arc::new(AtomicBool::new(false));
/// let done = Arc::new(AtomicU64::new(0));
/// let progress = Arc::clone(&done);
/// let options = TrajectoryOpenOptions::default()
///     .with_index_cancellation(Cancellation::from_flag(Arc::clone(&cancel)))
///     .with_index_progress(move |step| progress.store(step.byte_count(), Ordering::Relaxed));
/// // Another thread may read `done` for a progress bar and set `cancel`.
/// let reader = open_indexed_trajectory("ethane.xtc", topology, options)?;
/// let index = reader.to_frame_index();
/// println!("{} frames", index.frame_count());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn open_indexed_trajectory(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryOpenOptions,
) -> Result<IndexedFileTrajectoryReader, TrajectoryIoError> {
    open_indexed(path.as_ref(), topology.into(), None, options)
}

/// Opens an indexed path-backed trajectory reader from a previously exported
/// [`TrajectoryFrameIndex`] without scanning the frames it indexes.
///
/// The topology order contract and `options` are those of
/// [`open_indexed_trajectory`]. Opening checks that the file begins with the
/// bytes the index covers ([`TrajectoryFileFingerprint::is_prefix_of_file`]),
/// that its format, header, and first frame agree with the index, and reads
/// the frame offsets, times, and steps from the index alone. A file that
/// grew by appending since the index was exported is the same file with
/// more frames: the reader scans the appended bytes alone, verifying and
/// indexing every complete frame there as [`IndexedFileTrajectoryReader::refresh`]
/// does, with the index cancellation and progress callback of `options`.
/// Every frame read later is fully verified as it is decoded, and a random
/// read checks that the frame ends exactly where the index places the next
/// frame, or at the end of the indexed bytes, and records the index's time
/// and step, so a wrong index yields an error rather than a wrong frame.
/// DCD times follow the reader's time policy, so a DCD index opens only
/// under the policy it was built with. Metadata reports the frame count
/// and, for TRR, whether the frames mix scalar widths.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::IndexMismatch`] when the index was built
/// for a different file or the file changed since: the indexed bytes, the
/// format, the first frame offset, or DCD steps and times under the reader's
/// time policy disagree, or a DCD opened without
/// [`TrajectoryOpenOptions::with_growing`] does not hold the frames its
/// header counts. Returns the open and header errors of [`open_trajectory`],
/// the decoding error of the first invalid appended frame, and the frame
/// index limits of `options` bound the index as
/// [`TrajectoryIoErrorKind::LimitExceeded`]. Later reads return
/// [`TrajectoryIoErrorKind::IndexMismatch`] for a frame that does not end
/// where the index says or records another time or step, besides the
/// codec's decoding errors.
///
/// [`TrajectoryIoErrorKind::IndexMismatch`]: crate::TrajectoryIoErrorKind::IndexMismatch
/// [`TrajectoryIoErrorKind::LimitExceeded`]: crate::TrajectoryIoErrorKind::LimitExceeded
pub fn open_indexed_trajectory_from_index(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    index: TrajectoryFrameIndex,
    options: TrajectoryOpenOptions,
) -> Result<IndexedFileTrajectoryReader, TrajectoryIoError> {
    open_indexed(path.as_ref(), topology.into(), Some(index), options)
}

/// Opens an indexed reader, scanning the file unless `index` is given, and
/// then the bytes appended after the frames `index` covers.
fn open_indexed(
    path: &Path,
    topology: Arc<Topology>,
    index: Option<TrajectoryFrameIndex>,
    options: TrajectoryOpenOptions,
) -> Result<IndexedFileTrajectoryReader, TrajectoryIoError> {
    let label = path.display().to_string();
    let growing = options.growing;
    let mut file = File::open(path)
        .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &label, error))?;
    let length = file
        .metadata()
        .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &label, error))?
        .len();
    if growing && length == 0 {
        return Err(incomplete(None, &label));
    }
    // A scan must not see the bytes the file holds now change.
    let scanned_from = match index {
        Some(_) => None,
        None => Some(TrajectoryFileFingerprint::of_labelled_file(
            &mut file, &label,
        )?),
    };
    let mut reader = BufReader::new(file);
    let detection = detect::select_format(
        &mut reader,
        path,
        options.format_hint,
        &options.limits,
        &label,
        growing,
    )?;
    let format = detection.format;
    let imported = match index {
        Some(index) => {
            if index.format() != format {
                return Err(index_mismatch(
                    TrajectoryIoOperation::Open,
                    format,
                    &label,
                    format!("the index was built for a {} file", index.format()),
                ));
            }
            // Detection left the reader at the start with an empty buffer.
            if !index
                .fingerprint()
                .is_prefix_of_labelled_file(reader.get_mut(), &label)?
            {
                return Err(index_mismatch(
                    TrajectoryIoOperation::Open,
                    format,
                    &label,
                    "the file does not begin with the bytes the index covers",
                ));
            }
            if index.frame_count() > index_hard_capacity(&options.limits) {
                return Err(TrajectoryIoError::codec(
                    TrajectoryIoErrorKind::LimitExceeded,
                    TrajectoryIoOperation::Open,
                    Some(format),
                )
                .with_source_label(&label)
                .with_detail("the index frame count exceeds the configured index limits")
                .with_limit(LimitExceeded::new(
                    Resource::Frames,
                    index_hard_capacity(&options.limits),
                    index.frame_count(),
                )));
            }
            Some(index)
        }
        None => None,
    };
    let imported_count = imported.as_ref().map(TrajectoryFrameIndex::frame_count);
    let opened = match (&imported, scanned_from) {
        (Some(index), _) => index.fingerprint().clone(),
        (None, Some(fingerprint)) => fingerprint,
        (None, None) => unreachable!("an index or a scan fingerprint is taken"),
    };
    let imported = imported.map(TrajectoryFrameIndex::into_parts);
    let control = IndexBuildControl::new(
        options.index_cancellation.as_ref(),
        options.index_progress.as_ref(),
        length,
    );
    // A growing file whose first frame or header is still being written.
    let unfinished = |error: TrajectoryIoError| {
        if growing && error.kind() == TrajectoryIoErrorKind::TruncatedRecord {
            incomplete(Some(format), &label).with_source(error)
        } else {
            error
        }
    };
    let (inner, metadata, note) = match format {
        TrajectoryFormat::Xyz => {
            let xyz_options = options
                .xyz
                .with_limits(options.limits)
                .with_source_label(label.clone());
            let reader = match imported {
                Some((offsets, _)) => {
                    let mut reader = xyz::IndexedXyzReader::from_offsets(
                        reader,
                        topology,
                        xyz_options,
                        offsets,
                        growing,
                    )?;
                    reader.extend(growing, &control)?;
                    reader
                }
                None => {
                    xyz::IndexedXyzReader::build(reader, topology, xyz_options, growing, &control)
                        .map_err(unfinished)?
                }
            };
            let count = reader.frame_count().unwrap_or(0);
            let atom_count = reader.topology().atom_count();
            (
                IndexedReaderInner::Xyz(reader),
                FileTrajectoryMetadata::xyz(atom_count, Some(count)),
                "XYZ index verified every complete frame",
            )
        }
        TrajectoryFormat::Dcd => {
            let dcd_options = options
                .dcd
                .with_limits(options.limits)
                .with_source_label(label.clone());
            let reader = match imported {
                Some((offsets, _)) => {
                    let mut reader = dcd::IndexedDcdReader::from_offsets(
                        reader,
                        topology,
                        dcd_options,
                        offsets,
                        growing,
                    )?;
                    reader.extend(growing, &control).map_err(|error| {
                        if error.kind() == TrajectoryIoErrorKind::InconsistentMetadata
                            && error.expected().is_some()
                        {
                            // The frames do not add up to the header's count.
                            index_mismatch(
                                TrajectoryIoOperation::Open,
                                TrajectoryFormat::Dcd,
                                &label,
                                "DCD index frame count differs from the header NSET count",
                            )
                            .with_counts(error.expected().unwrap_or(0), error.actual().unwrap_or(0))
                            .with_source(error)
                        } else {
                            error
                        }
                    })?;
                    reader
                }
                None => {
                    dcd::IndexedDcdReader::build(reader, topology, dcd_options, growing, &control)
                        .map_err(unfinished)?
                }
            };
            let count = reader.frame_count().unwrap_or(0);
            let metadata = FileTrajectoryMetadata::dcd(reader.header(), Some(count));
            (
                IndexedReaderInner::Dcd(reader),
                metadata,
                "DCD index verified record markers and the declared frame count",
            )
        }
        TrajectoryFormat::Trr => {
            let trr_options = options
                .trr
                .with_limits(options.limits)
                .with_source_label(label.clone());
            let reader = match imported {
                Some((offsets, precision_mixed)) => {
                    let mut reader = trr::IndexedTrrReader::from_offsets(
                        reader,
                        topology,
                        trr_options,
                        offsets,
                        precision_mixed,
                    )?;
                    reader.extend(growing, &control)?;
                    reader
                }
                None => {
                    trr::IndexedTrrReader::build(reader, topology, trr_options, growing, &control)
                        .map_err(unfinished)?
                }
            };
            let count = reader.frame_count().unwrap_or(0);
            let metadata = FileTrajectoryMetadata::trr(
                reader.first_header(),
                Some(count),
                reader.precision_mixed(),
            );
            (
                IndexedReaderInner::Trr(Box::new(reader)),
                metadata,
                "TRR index verified every XDR frame and payload block",
            )
        }
        TrajectoryFormat::Xtc => {
            let cell_policy = options.xtc.cell_policy();
            let xtc_options = options
                .xtc
                .with_limits(options.limits)
                .with_source_label(label.clone());
            let reader = match imported {
                Some((offsets, _)) => {
                    let mut reader = xtc::IndexedXtcReader::from_offsets(
                        reader,
                        topology,
                        xtc_options,
                        offsets,
                    )?;
                    reader.extend(growing, &control)?;
                    reader
                }
                None => {
                    xtc::IndexedXtcReader::build(reader, topology, xtc_options, growing, &control)
                        .map_err(unfinished)?
                }
            };
            let count = reader.frame_count().unwrap_or(0);
            let metadata =
                FileTrajectoryMetadata::xtc(reader.first_info(), cell_policy, Some(count));
            (
                IndexedReaderInner::Xtc(reader),
                metadata,
                "XTC index fully decoded and verified every compressed frame",
            )
        }
    };
    if growing && format != TrajectoryFormat::Dcd && inner.offsets().is_empty() {
        // Only DCD, whose header declares the frame count, opens without a
        // frame; the others need their first frame for its atoms and layout.
        return Err(incomplete(Some(format), &label));
    }
    let count = inner.offsets().len();
    let fingerprint = if imported_count == Some(count) {
        // Nothing was scanned, so nothing scanned can have changed.
        opened
    } else {
        let mut file = File::open(path)
            .map_err(|error| io_context(TrajectoryIoOperation::Index, None, &label, error))?;
        if !opened.is_prefix_of_labelled_file(&mut file, &label)? {
            return Err(index_mismatch(
                TrajectoryIoOperation::Index,
                format,
                &label,
                "the file changed while it was indexed",
            ));
        }
        let end = inner.offsets().end();
        if end == opened.byte_count() {
            opened
        } else {
            TrajectoryFileFingerprint::of_labelled_file_prefix(&mut file, end, &label)?
        }
    };
    let note = match imported_count {
        Some(imported) if imported == count => {
            format!("{format} index imported; every frame is verified as it is read")
        }
        Some(imported) => format!(
            "{format} index imported and {} appended frames verified; every frame is verified as it is read",
            count - imported
        ),
        None => note.to_owned(),
    };
    let report = TrajectoryOpenReport {
        selected_format: format,
        detection_evidence: detection.evidence,
        notes: vec![note],
    };
    Ok(IndexedFileTrajectoryReader {
        inner,
        metadata,
        open_report: report,
        fingerprint,
        path: path.to_path_buf(),
        growing,
    })
}

/// A growing file that does not yet hold its first frame or header.
fn incomplete(format: Option<TrajectoryFormat>, source_label: &str) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        TrajectoryIoErrorKind::Incomplete,
        TrajectoryIoOperation::Open,
        format,
    )
    .with_source_label(source_label)
    .with_detail("the file does not yet hold a complete frame; open it again once more is written")
}
