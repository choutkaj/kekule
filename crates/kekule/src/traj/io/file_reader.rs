use std::fmt;
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::topology::Topology;

use crate::diagnostics::{LimitExceeded, Resource};
use crate::structure::{DecodedGeometry, Geometry, StructureSet};
use crate::traj::{
    FrameBuffer, MemoryTrajectoryWriter, SeekableTrajectoryReader, TrajectoryFormat,
    TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation, TrajectoryReader,
};

use crate::traj::trajectory::COLUMN_MISMATCH;

use super::limits::index_hard_capacity;
use super::{
    dcd, detect, index_mismatch, io_context, trr, xtc, xyz, FileTrajectoryMetadata, FrameOffsets,
    IndexBuildControl, SharedFile, TrajectoryFileFingerprint, TrajectoryFrameIndex,
    TrajectoryOpenOptions, TrajectoryOpenReport, TrajectoryWriteOptions,
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

/// Every indexed reader reads its file through a shared position, so that
/// readers forked from it read the same file object.
enum IndexedReaderInner {
    Xyz(xyz::IndexedXyzReader<BufReader<SharedFile>>),
    Dcd(dcd::IndexedDcdReader<BufReader<SharedFile>>),
    Trr(Box<trr::IndexedTrrReader<BufReader<SharedFile>>>),
    Xtc(xtc::IndexedXtcReader<BufReader<SharedFile>>),
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

    fn shared_offsets(&self) -> &Arc<FrameOffsets> {
        match self {
            Self::Xyz(reader) => reader.shared_offsets(),
            Self::Dcd(reader) => reader.shared_offsets(),
            Self::Trr(reader) => reader.shared_offsets(),
            Self::Xtc(reader) => reader.shared_offsets(),
        }
    }

    /// See the codecs' `release_index`.
    fn release_index(&mut self) {
        match self {
            Self::Xyz(reader) => reader.release_index(),
            Self::Dcd(reader) => reader.release_index(),
            Self::Trr(reader) => reader.release_index(),
            Self::Xtc(reader) => reader.release_index(),
        }
    }

    /// Reads by the index of `other`, a reader of the same file, and the
    /// file object it reads, from now on.
    fn share_index(&mut self, other: &Self) {
        match (&mut *self, other) {
            (Self::Xyz(reader), Self::Xyz(other)) => reader.share_index(other),
            (Self::Dcd(reader), Self::Dcd(other)) => reader.share_index(other),
            (Self::Trr(reader), Self::Trr(other)) => reader.share_index(other),
            (Self::Xtc(reader), Self::Xtc(other)) => reader.share_index(other),
            _ => unreachable!("readers of one file read one format"),
        }
        let file = other.file();
        let stream = self.stream_mut();
        if !stream.get_ref().same_file(file) {
            // The buffered bytes were read but not yet consumed.
            let position = stream.get_ref().position() - stream.buffer().len() as u64;
            *stream = BufReader::new(file.at(position));
        }
    }

    /// Another reader of the same frames and file object, sharing the
    /// index; see the codecs' `fork`.
    fn fork(&self) -> Result<Self, TrajectoryIoError> {
        let stream = BufReader::new(self.file().at(0));
        Ok(match self {
            Self::Xyz(reader) => Self::Xyz(reader.fork(stream)?),
            Self::Dcd(reader) => Self::Dcd(reader.fork(stream)?),
            Self::Trr(reader) => Self::Trr(Box::new(reader.fork(stream)?)),
            Self::Xtc(reader) => Self::Xtc(reader.fork(stream)?),
        })
    }

    /// See the codecs' `read_decoded`.
    fn read_decoded(&mut self, index: usize) -> Result<DecodedGeometry, TrajectoryIoError> {
        match self {
            Self::Xyz(reader) => reader.read_decoded(index),
            Self::Dcd(reader) => reader.read_decoded(index),
            Self::Trr(reader) => reader.read_decoded(index),
            Self::Xtc(reader) => reader.read_decoded(index),
        }
    }

    /// The file every read goes through.
    fn file(&self) -> &SharedFile {
        match self {
            Self::Xyz(reader) => reader.stream().get_ref(),
            Self::Dcd(reader) => reader.stream().get_ref(),
            Self::Trr(reader) => reader.stream().get_ref(),
            Self::Xtc(reader) => reader.stream().get_ref(),
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

    /// The stream every read goes through.
    fn stream_mut(&mut self) -> &mut BufReader<SharedFile> {
        match self {
            Self::Xyz(reader) => reader.stream_mut(),
            Self::Dcd(reader) => reader.stream_mut(),
            Self::Trr(reader) => reader.stream_mut(),
            Self::Xtc(reader) => reader.stream_mut(),
        }
    }

    /// Runs `check` on the file the reader reads, then puts the stream back
    /// where the reader left it.
    fn check_file<T>(
        &mut self,
        format: TrajectoryFormat,
        label: &str,
        check: impl FnOnce(&mut SharedFile) -> Result<T, TrajectoryIoError>,
    ) -> Result<T, TrajectoryIoError> {
        let io = |error| io_context(TrajectoryIoOperation::Index, Some(format), label, error);
        let stream = self.stream_mut();
        let position = stream.stream_position().map_err(io)?;
        let checked = check(stream.get_mut());
        // Seeking also drops the buffer, which the check left behind.
        let restored = stream.seek(SeekFrom::Start(position));
        let value = checked?;
        restored.map_err(io)?;
        Ok(value)
    }

    /// Reads on from `file`, which must begin with the indexed bytes, where
    /// the reader left its stream, and returns the stream it replaces.
    fn rebind(
        &mut self,
        file: File,
        format: TrajectoryFormat,
        label: &str,
    ) -> Result<BufReader<SharedFile>, TrajectoryIoError> {
        let io = |error| io_context(TrajectoryIoOperation::Index, Some(format), label, error);
        let stream = self.stream_mut();
        let position = stream.stream_position().map_err(io)?;
        let mut replacement = BufReader::new(SharedFile::new(file));
        replacement.seek(SeekFrom::Start(position)).map_err(io)?;
        Ok(std::mem::replace(stream, replacement))
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
    /// the indexed bytes, for [`open_indexed_trajectory_from_index`]. The
    /// index shares the reader's offsets until either changes them.
    pub fn to_frame_index(&self) -> TrajectoryFrameIndex {
        TrajectoryFrameIndex::new(
            self.metadata.format(),
            self.fingerprint.clone(),
            Arc::clone(self.inner.shared_offsets()),
            self.inner.precision_mixed(),
        )
    }

    /// The path the reader opened.
    pub(crate) fn path(&self) -> &Path {
        &self.path
    }

    /// The fingerprint of the indexed bytes.
    pub(crate) const fn fingerprint(&self) -> &TrajectoryFileFingerprint {
        &self.fingerprint
    }

    /// Another reader of the indexed frames that reads the same file object
    /// and shares the index. Opening it reads the start of the file alone:
    /// the header, and the first frame where the format checks it.
    pub(crate) fn fork(&self) -> Result<Self, TrajectoryIoError> {
        Ok(Self {
            inner: self.inner.fork()?,
            metadata: self.metadata.clone(),
            open_report: self.open_report.clone(),
            fingerprint: self.fingerprint.clone(),
            path: self.path.clone(),
            growing: self.growing,
        })
    }

    /// Lets go of the index until [`Self::share_index`], so that refreshing
    /// another reader extends the index without copying it.
    pub(crate) fn release_index(&mut self) {
        self.inner.release_index();
    }

    /// Reads by the index of `other`, a reader of the same file, and the
    /// file object it reads, from now on.
    pub(crate) fn share_index(&mut self, other: &Self) {
        self.inner.share_index(&other.inner);
        self.fingerprint = other.fingerprint.clone();
        self.metadata = other.metadata.clone();
    }

    /// Reads frame `index` as a geometry, as [`SeekableTrajectoryReader::read_frame`]
    /// verifies it, moving the decoded arrays into the geometry.
    pub(crate) fn read_geometry(&mut self, index: usize) -> Result<Geometry, TrajectoryIoError> {
        let decoded = self.inner.read_decoded(index)?;
        Ok(Geometry::from_decoded(decoded, self.topology())?)
    }

    /// The error for a read of frame `index` that failed with `error`: one
    /// naming the frame as lost when the file no longer holds all of it,
    /// whatever the read reported, otherwise `error`.
    pub(crate) fn lost_frame(&self, index: usize, error: TrajectoryIoError) -> TrajectoryIoError {
        let Some((_, end)) = self.inner.offsets().frame(index) else {
            return error;
        };
        match self.inner.file().len() {
            Ok(length) if length < end => TrajectoryIoError::codec(
                TrajectoryIoErrorKind::IndexMismatch,
                TrajectoryIoOperation::ReadFrame,
                Some(self.metadata.format()),
            )
            .with_source_label(self.path.display().to_string())
            .with_frame(index)
            .with_detail(format!(
                "The file no longer holds frame {index}: it is {length} bytes, and the frame ends at byte {end}. Open the file again."
            )),
            _ => error,
        }
    }

    /// The bytes the file held after the last indexed frame when the reader
    /// last opened or [refreshed](Self::refresh): a frame still being
    /// written, or more that arrived during the scan. Zero when the indexed
    /// frames reached the end of the file. Every refresh reads them again.
    /// Only a reader opened with [`TrajectoryOpenOptions::with_growing`]
    /// leaves a partial frame pending; for others it is an error.
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
    /// A refresh reads the file the path names when it runs: a file moved
    /// to the path in place of the one read so far is read from then on
    /// when it begins with the indexed bytes, as a copy of the file grown
    /// further does. The refresh checks, scans, and fingerprints that one
    /// file, so the index and its fingerprint always describe the same bytes.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::IndexMismatch`] when the indexed bytes
    /// changed, including while the new frames were read, the decoding error
    /// of the first invalid new frame,
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when the index outgrows the
    /// limits, and [`TrajectoryIoErrorKind::Io`] when the file cannot be
    /// read. Failure changes nothing: the reader keeps its frames and the
    /// file it read them from, and can be refreshed again. A file whose
    /// indexed bytes changed, such as a truncated or rewritten one, never
    /// refreshes; opening it again is the way back.
    ///
    /// ```no_run
    /// use kekule::{smiles, topology::Topology};
    /// use kekule::traj::io::{open_indexed_trajectory, TrajectoryOpenOptions};
    /// use kekule::traj::{SeekableTrajectoryReader, TrajectoryReader};
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
    /// [`TrajectoryIoErrorKind::IndexMismatch`]: crate::traj::TrajectoryIoErrorKind::IndexMismatch
    /// [`TrajectoryIoErrorKind::LimitExceeded`]: crate::traj::TrajectoryIoErrorKind::LimitExceeded
    /// [`TrajectoryIoErrorKind::Io`]: crate::traj::TrajectoryIoErrorKind::Io
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
                "The file's indexed bytes changed since they were indexed. Open the file again.",
            )
        };
        if !self
            .fingerprint
            .is_prefix_of_labelled_file(&mut file, &label)?
        {
            return Err(changed());
        }
        // The scan, the checks, and the fingerprint all read the file
        // checked here, which may not be the one read so far.
        let previous = self.inner.rebind(file, format, &label)?;
        let before = self.inner.offsets().len();
        let snapshot = self.inner.snapshot();
        let indexed = &self.fingerprint;
        let checked = self
            .inner
            .extend(self.growing, &IndexBuildControl::default())
            .and_then(|_| {
                let end = self.inner.offsets().end();
                self.inner.check_file(format, &label, |file| {
                    if indexed.is_prefix_of_labelled_file(file, &label)? {
                        TrajectoryFileFingerprint::of_labelled_file_prefix(file, end, &label)
                    } else {
                        Err(changed())
                    }
                })
            });
        match checked {
            Ok(fingerprint) => {
                self.fingerprint = fingerprint;
                self.record_counts();
                Ok(self.inner.offsets().len() - before)
            }
            Err(error) => {
                self.inner.restore(snapshot);
                *self.inner.stream_mut() = previous;
                Err(error)
            }
        }
    }

    /// Fails unless a writer with `options` continues the indexed frames as
    /// the file's writer would have: under the same DCD header, or with the
    /// magic and precision of the first XTC frame or the scalar width of the
    /// TRR frames. XYZ records none of its writer's settings.
    pub(crate) fn check_append(
        &mut self,
        options: &TrajectoryWriteOptions,
    ) -> Result<(), TrajectoryIoError> {
        let label = self.path.display().to_string();
        match &self.inner {
            IndexedReaderInner::Xyz(_) => Ok(()),
            IndexedReaderInner::Xtc(reader) => {
                xtc::check_append_frames(reader.first_info(), options.xtc_options(), &label)
            }
            IndexedReaderInner::Trr(reader) => trr::check_append_frames(
                reader.first_header(),
                reader.precision_mixed(),
                options.trr_options(),
                &label,
            ),
            IndexedReaderInner::Dcd(_) => {
                let atoms = self.metadata.atom_count();
                let dcd = *options.dcd_options();
                self.inner
                    .check_file(TrajectoryFormat::Dcd, &label, |file| {
                        file.seek(SeekFrom::Start(0)).map_err(|error| {
                            io_context(
                                TrajectoryIoOperation::Open,
                                Some(TrajectoryFormat::Dcd),
                                &label,
                                error,
                            )
                        })?;
                        dcd::check_append_header(file, dcd, atoms, &label)
                    })
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

/// Reads an entire trajectory file into a structure set in memory, one geometry
/// per frame in file frame order: the frames' times and steps become a time
/// order, and each frame's positions, cell, properties, velocities, and
/// forces its geometry's. The set has no distribution.
///
/// Accepts an owned topology or a shared `Arc<Topology>` and detects the format
/// automatically. File coordinates must follow the topology's dense atom order;
/// counts and available format metadata are checked, but matching counts alone
/// cannot establish atom identity. All decoded frame fields are retained, and
/// a supplied shared topology retains its exact allocation.
///
/// This consumes the sequential reader through clean EOF. Any open or decoding
/// error is returned without publishing a partial set. Memory usage grows
/// with the number of frames and atoms; use [`open_trajectory`] to process
/// large files one frame at a time or inspect file metadata and opening
/// diagnostics, and [`open_indexed_structure_set`](super::open_indexed_structure_set)
/// to read geometries on demand.
///
/// `options` sets the format policies and limits; frame limits cause an error
/// rather than truncating the result. These codec limits do not impose a total
/// resident-memory budget for the loaded set. No random-access index is
/// built.
///
/// # Errors
///
/// Returns the errors of [`open_trajectory`], the first decoding, validation
/// or [`TrajectoryIoErrorKind::Io`] error of any frame, and
/// [`TrajectoryIoErrorKind::InconsistentMetadata`] for the first frame
/// without a time or step that the first frame has, or with one the first
/// frame lacks, with [`TrajectoryIoError::frame`] its index; the
/// [`StructureSetError`](crate::structure::StructureSetError) is the source.
///
/// [`TrajectoryIoErrorKind::Io`]: crate::traj::TrajectoryIoErrorKind::Io
/// [`TrajectoryIoErrorKind::InconsistentMetadata`]:
///     crate::traj::TrajectoryIoErrorKind::InconsistentMetadata
///
/// ```no_run
/// use kekule::mmcif;
/// use kekule::traj::io::read_structure_set;
///
/// let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?, Default::default())?;
/// let topology = document.interpret(Default::default())?.into_topology();
/// let frames = read_structure_set("trajectory.xtc", topology, Default::default())?;
/// println!("{} frames, {} atoms", frames.len(), frames.topology().atom_count());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn read_structure_set(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryOpenOptions,
) -> Result<StructureSet, TrajectoryIoError> {
    let path = path.as_ref();
    let mut reader = open_trajectory(path, topology, options)?;
    let mut buffer = reader.frame_buffer();
    let mut writer = MemoryTrajectoryWriter::new(reader.shared_topology());
    let format = reader.metadata().format();
    let mut index = 0;
    while reader.read_next(&mut buffer)? {
        writer.collect(buffer.frame_view()).map_err(|error| {
            TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::ReadFrame,
                Some(format),
            )
            .with_source_label(path.display().to_string())
            .with_frame(index)
            .with_detail(COLUMN_MISMATCH)
            .with_source(error)
        })?;
        index += 1;
    }
    Ok(writer.into_structure_set())
}

/// Opens one fast sequential path-backed trajectory reader without loading all frames.
///
/// Accepts an owned topology or a shared `Arc<Topology>`. File coordinates must
/// follow its dense atom order; counts and available format metadata are checked
/// automatically. Matching counts alone cannot establish atom identity.
/// `options` sets the format policies and limits. For a fully loaded in-memory
/// structure set, use [`read_structure_set`] instead.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::Io`] when the file cannot be opened or
/// read, the detection errors of [`super::detect_trajectory_format`], and the
/// header errors of the selected codec's reader, such as
/// [`TrajectoryIoErrorKind::InconsistentAtomCount`] when the file's atom count
/// differs from the topology's.
///
/// [`TrajectoryIoErrorKind::Io`]: crate::traj::TrajectoryIoErrorKind::Io
/// [`TrajectoryIoErrorKind::InconsistentAtomCount`]:
///     crate::traj::TrajectoryIoErrorKind::InconsistentAtomCount
///
/// ```no_run
/// use kekule::{smiles, topology::Topology};
/// use kekule::traj::{io::open_trajectory, TrajectoryReader};
///
/// let molecule = smiles::to_molecules("CC")?.pop().unwrap();
/// let topology = Topology::from_molecule(molecule)?;
/// let mut reader = open_trajectory("ethane.xyz", topology, Default::default())?;
/// let mut frame = reader.frame_buffer();
/// while reader.read_next(&mut frame)? {
///     let structure = frame.as_structure_view();
///     assert_eq!(structure.atom_count(), 2);
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
/// [`TrajectoryIoErrorKind::LimitExceeded`]: crate::traj::TrajectoryIoErrorKind::LimitExceeded
/// [`TrajectoryIoErrorKind::Cancelled`]: crate::traj::TrajectoryIoErrorKind::Cancelled
/// [`TrajectoryIoErrorKind::IndexMismatch`]: crate::traj::TrajectoryIoErrorKind::IndexMismatch
/// [`TrajectoryIoErrorKind::Incomplete`]: crate::traj::TrajectoryIoErrorKind::Incomplete
///
/// ```no_run
/// use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
/// use std::sync::Arc;
///
/// use kekule::diagnostics::Cancellation;
/// use kekule::{smiles, topology::Topology};
/// use kekule::traj::io::{open_indexed_trajectory, TrajectoryOpenOptions};
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
/// [`TrajectoryIoErrorKind::IndexMismatch`]: crate::traj::TrajectoryIoErrorKind::IndexMismatch
/// [`TrajectoryIoErrorKind::LimitExceeded`]: crate::traj::TrajectoryIoErrorKind::LimitExceeded
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
    let mut reader = BufReader::new(SharedFile::new(file));
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
                    format!("The index was built for a {} file.", index.format()),
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
                    "The file does not begin with the bytes the index covers.",
                ));
            }
            if index.frame_count() > index_hard_capacity(&options.limits) {
                return Err(TrajectoryIoError::codec(
                    TrajectoryIoErrorKind::LimitExceeded,
                    TrajectoryIoOperation::Open,
                    Some(format),
                )
                .with_source_label(&label)
                .exceeding(
                    "The index",
                    LimitExceeded::new(
                        Resource::Frames,
                        index_hard_capacity(&options.limits),
                        index.frame_count(),
                    ),
                ));
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
    let (mut inner, metadata, note) = match format {
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
                                "The index's frame count differs from the frame count NSET in the header.",
                            )
                            .with_written_counts(
                                error.expected().unwrap_or(0),
                                error.actual().unwrap_or(0),
                            )
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
        // In the file that was scanned, whatever the path names now.
        let end = inner.offsets().end();
        inner.check_file(format, &label, |file| {
            if !opened.is_prefix_of_labelled_file(file, &label)? {
                return Err(index_mismatch(
                    TrajectoryIoOperation::Index,
                    format,
                    &label,
                    "The file changed while it was indexed.",
                ));
            }
            if end == opened.byte_count() {
                Ok(opened)
            } else {
                TrajectoryFileFingerprint::of_labelled_file_prefix(file, end, &label)
            }
        })?
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
    .with_detail("The file does not yet hold a complete frame. Open it again once more is written.")
}
