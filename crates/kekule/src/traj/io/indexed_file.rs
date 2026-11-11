//! A sample set over an indexed trajectory file: [`IndexedFile`].

use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, PoisonError};

use crate::structure::{
    sample_storage, Distribution, Order, Sample, SampleRead, SampleReadError, SampleSet,
    SampleStorage, SampleView, TimeAxis,
};
use crate::topology::Topology;
use crate::traj::{FrameBuffer, SeekableTrajectoryReader, TrajectoryIoError, TrajectoryReader};

use super::{
    open_indexed_trajectory, open_indexed_trajectory_from_index, IndexedFileTrajectoryReader,
    TrajectoryFileFingerprint, TrajectoryFrameIndex, TrajectoryOpenOptions,
};

/// The storage of a [`SampleSet`] read from an indexed trajectory file: the
/// file's path and frame index, and a pool of readers.
///
/// The set holds no sample in memory. [`SampleRead::read`] decodes frame
/// `i` into the returned [`Sample`], positions, cell, velocities, forces,
/// and properties, without a second copy of the coordinates;
/// [`SampleRead::peek`] lends none. Reads from several threads run at once:
/// each takes a reader of the pool, or opens one from the frame index when
/// all are busy, so the pool grows to the most reads that ran at once.
/// The order is a time order of the index's times and steps, as
/// [`read_sample_set`](super::read_sample_set) gives them.
///
/// Open one with [`open_indexed_sample_set`] or
/// [`open_indexed_sample_set_from_index`], follow a growing file with
/// [`SampleSet::refresh`], and load samples into memory with
/// [`SampleSet::load`].
pub struct IndexedFile {
    path: PathBuf,
    /// Opens the pool's readers; without index cancellation or progress.
    options: TrajectoryOpenOptions,
    index: TrajectoryFrameIndex,
    idle: Mutex<Vec<IndexedFileTrajectoryReader>>,
}

impl fmt::Debug for IndexedFile {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedFile")
            .field("path", &self.path)
            .field("index", &self.index)
            .finish_non_exhaustive()
    }
}

impl IndexedFile {
    /// The trajectory file.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The frame index of the samples: one frame per sample. A host keeps it
    /// with [`TrajectoryFrameIndex::to_bytes`] and opens the file again from
    /// it without scanning ([`open_indexed_sample_set_from_index`]).
    pub fn frame_index(&self) -> &TrajectoryFrameIndex {
        &self.index
    }

    /// The fingerprint of the bytes the frame index covers.
    pub fn fingerprint(&self) -> &TrajectoryFileFingerprint {
        self.index.fingerprint()
    }

    /// The pool's readers waiting for a read.
    fn idle(&self) -> std::sync::MutexGuard<'_, Vec<IndexedFileTrajectoryReader>> {
        self.idle.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// A reader of the pool, or a new one opened from the index.
    fn take(
        &self,
        topology: &Arc<Topology>,
    ) -> Result<IndexedFileTrajectoryReader, TrajectoryIoError> {
        if let Some(reader) = self.idle().pop() {
            return Ok(reader);
        }
        open_indexed_trajectory_from_index(
            &self.path,
            Arc::clone(topology),
            self.index.clone(),
            self.options.clone(),
        )
    }

    /// Returns a reader to the pool.
    fn give_back(&self, reader: IndexedFileTrajectoryReader) {
        self.idle().push(reader);
    }
}

impl sample_storage::Storage for IndexedFile {
    fn sample_count(&self) -> usize {
        self.index.frame_count()
    }
}

impl SampleStorage for IndexedFile {}

/// Opens a trajectory file as a [`SampleSet`] over its frames, indexing
/// them through [`open_indexed_trajectory`], which verifies every frame.
///
/// The set's order is a time order of the frames' times and steps, and it
/// has no distribution; give it one with [`SampleSet::set_distribution`].
/// `options` set the codec, its limits, and the indexing's cancellation and
/// progress; readers the set opens later use the same options without
/// them.
///
/// ```no_run
/// use kekule::structure::SampleRead;
/// use kekule::traj::io::open_indexed_sample_set;
///
/// let topology = kekule::smiles::to_topology("CC")?;
/// let set = open_indexed_sample_set("run.xtc", topology, Default::default())?;
/// // Keep the index to open the file again without scanning it.
/// let index = set.storage().frame_index().to_bytes();
/// let sample = set.read(set.len() - 1)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
///
/// The errors of [`open_indexed_trajectory`].
pub fn open_indexed_sample_set(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryOpenOptions,
) -> Result<SampleSet<IndexedFile>, TrajectoryIoError> {
    let path = path.as_ref();
    let reader = open_indexed_trajectory(path, topology, options.clone())?;
    Ok(sample_set(path, reader, options))
}

/// Opens a trajectory file as a [`SampleSet`] over its frames from a frame
/// index the host kept, through [`open_indexed_trajectory_from_index`],
/// without scanning the frames it covers; as
/// [`open_indexed_sample_set`] does otherwise. Frames appended to the file
/// after the index are indexed too.
///
/// # Errors
///
/// The errors of [`open_indexed_trajectory_from_index`].
pub fn open_indexed_sample_set_from_index(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    index: TrajectoryFrameIndex,
    options: TrajectoryOpenOptions,
) -> Result<SampleSet<IndexedFile>, TrajectoryIoError> {
    let path = path.as_ref();
    let reader = open_indexed_trajectory_from_index(path, topology, index, options.clone())?;
    Ok(sample_set(path, reader, options))
}

/// The set over `reader`'s frames, `reader` the first of its pool.
fn sample_set(
    path: &Path,
    reader: IndexedFileTrajectoryReader,
    mut options: TrajectoryOpenOptions,
) -> SampleSet<IndexedFile> {
    options.index_cancellation = None;
    options.index_progress = None;
    let index = reader.to_frame_index();
    let topology = reader.shared_topology();
    let order = time_order(&index);
    let storage = IndexedFile {
        path: path.to_path_buf(),
        options,
        index,
        idle: Mutex::new(vec![reader]),
    };
    SampleSet::from_storage(topology, storage, Some(order), None)
}

/// The time order of an index's frames.
fn time_order(index: &TrajectoryFrameIndex) -> Order {
    Order::Time(TimeAxis::from_columns(
        index.times().map(|times| times.value().to_vec()),
        index.steps().map(<[u64]>::to_vec),
    ))
}

impl SampleSet<IndexedFile> {
    /// Indexes the frames appended to the file since the set was opened or
    /// last refreshed, as [`IndexedFileTrajectoryReader::refresh`] does, and
    /// takes them as new samples; returns how many there are.
    ///
    /// The new samples extend the time order with their times and steps. A
    /// distribution over every sample becomes one over the samples it
    /// covered, so the new samples are outside it, as
    /// [`SampleSet::push`] leaves an appended sample; a distribution over a
    /// subset keeps its members. Readers opened from the old index are
    /// dropped. Failure changes nothing.
    ///
    /// ```no_run
    /// use kekule::structure::SampleRead;
    /// use kekule::traj::io::{open_indexed_sample_set, TrajectoryOpenOptions};
    ///
    /// let topology = kekule::smiles::to_topology("CC")?;
    /// // The trajectory of a running simulation.
    /// let options = TrajectoryOpenOptions::default().with_growing(true);
    /// let mut set = open_indexed_sample_set("run.xtc", topology, options)?;
    /// loop {
    ///     let added = set.refresh()?;
    ///     for index in set.len() - added..set.len() {
    ///         let sample = set.read(index)?;
    ///         // Show the new sample.
    ///     }
    ///     # break;
    ///     // Wait for the simulation to write more.
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// The errors of [`IndexedFileTrajectoryReader::refresh`], and of
    /// [`open_indexed_trajectory_from_index`] when the set holds no reader.
    pub fn refresh(&mut self) -> Result<usize, TrajectoryIoError> {
        let topology = self.shared_topology();
        let storage = self.storage();
        let mut reader = storage.take(&topology)?;
        if let Err(error) = reader.refresh() {
            storage.give_back(reader);
            return Err(error);
        }
        // A reader the pool opened after the file grew indexed more frames
        // than the set holds, so the new samples are counted from the set.
        let previous = self.len();
        let count = reader.frame_count().unwrap_or(previous);
        if count <= previous {
            storage.give_back(reader);
            return Ok(0);
        }
        let added = count - previous;
        let index = reader.to_frame_index();
        let order = time_order(&index);
        let storage = self.storage_mut();
        storage.index = index;
        *storage
            .idle
            .get_mut()
            .unwrap_or_else(PoisonError::into_inner) = vec![reader];
        self.grown(previous, Some(order));
        Ok(added)
    }

    /// Reads the samples `indices` into a set in memory, in the requested
    /// order: each with its order value, and the distribution restricted to
    /// them, as [`SampleSet::select`] restricts it. Ranges and strides can be
    /// written as `(start..end).step_by(stride)`, and every sample as
    /// `0..set.len()`; duplicates and reordering are allowed. Set
    /// properties are kept.
    ///
    /// # Errors
    ///
    /// [`SampleReadError::IndexOutOfRange`] for an index past the last
    /// sample, and the first failure to read a sample.
    pub fn load(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<SampleSet, SampleReadError> {
        let indices = indices.into_iter().collect::<Vec<_>>();
        let samples = indices
            .iter()
            .map(|index| self.read(*index))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(self.loaded(&indices, samples))
    }
}

impl SampleRead for SampleSet<IndexedFile> {
    fn topology(&self) -> &Arc<Topology> {
        self.topology_handle()
    }

    fn len(&self) -> usize {
        SampleSet::len(self)
    }

    fn order(&self) -> Option<&Order> {
        SampleSet::order(self)
    }

    fn distribution(&self) -> Option<&Distribution> {
        SampleSet::distribution(self)
    }

    fn peek(&self, _index: usize) -> Option<SampleView<'_>> {
        None
    }

    fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
        let len = SampleSet::len(self);
        if index >= len {
            return Err(SampleReadError::IndexOutOfRange { index, len });
        }
        let topology = self.topology_handle();
        let storage = self.storage();
        let mut reader = storage.take(topology)?;
        let mut buffer = FrameBuffer::new(Arc::clone(topology));
        let read = reader.read_frame(index, &mut buffer);
        storage.give_back(reader);
        read?;
        let (_, _, mut frames) = buffer.into_trajectory().into_parts();
        let frame = frames.pop().expect("a frame buffer holds one frame");
        Ok(Sample::from_frame(frame))
    }
}
