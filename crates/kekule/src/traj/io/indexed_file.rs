//! A structure set over an indexed trajectory file: [`IndexedFile`].

use std::fmt;
use std::path::Path;
use std::sync::{Arc, Mutex, PoisonError};

use crate::properties::OwnerProperties;
use crate::structure::{
    geometry_storage, Distribution, Geometry, GeometryStorage, GeometryView, Order, StructureSet,
    StructureSetRead, StructureSetReadError, TimeAxis,
};
use crate::topology::Topology;
use crate::traj::{SeekableTrajectoryReader, TrajectoryIoError, TrajectoryReader};

use super::{
    open_indexed_trajectory, open_indexed_trajectory_from_index, IndexedFileTrajectoryReader,
    TrajectoryFileFingerprint, TrajectoryFrameIndex, TrajectoryOpenOptions,
};

/// The storage of a [`StructureSet`] read from an indexed trajectory file: the
/// open file, its frame index, and a pool of readers.
///
/// The set holds no geometry in memory. [`StructureSetRead::read`] decodes frame
/// `i` into the returned [`Geometry`]: positions, cell, velocities, forces,
/// and properties, written once into the arrays the geometry keeps.
/// [`StructureSetRead::peek`] lends none. Reads from several threads run at once:
/// each takes a reader of the pool, or opens one when all are busy, so the
/// pool grows to the most reads that ran at once.
///
/// Every reader reads the file object the set opened and shares the set's
/// frame index, so opening one reads the start of the file alone and a read
/// gives the same geometry whichever reader serves it. The readers share one
/// file handle. On Windows, reads of one handle run one after another, so
/// reads from several threads take turns reading the file and decode at
/// once.
///
/// A read decodes whatever bytes the file holds at the indexed offsets;
/// kekule does not watch the file.
/// - A file truncated while the set is open fails every read of the frames
///   it no longer holds with one error,
///   [`TrajectoryIoErrorKind::IndexMismatch`] naming the frame. The frames
///   before them read as before.
/// - Bytes rewritten in place are read as they now are. A rewritten frame
///   of the same size reads without error, until [`StructureSet::refresh`]
///   finds that the indexed bytes changed and fails.
/// - A file replaced at its path keeps reading through the open file. A
///   refresh opens the path again. When the new file begins with the
///   indexed bytes, the refresh takes it, and every reader reads it from
///   then on. Otherwise the refresh fails.
///
/// After a failed refresh, opening the file again is the way back. A host
/// that must know of changes between refreshes compares
/// [`Self::fingerprint`] with [`TrajectoryFileFingerprint::is_prefix_of_path`].
///
/// The order is a time order of the index's times and steps, as
/// [`read_structure_set`](super::read_structure_set) gives them.
///
/// Open one with [`open_indexed_structure_set`] or
/// [`open_indexed_structure_set_from_index`], or make one from a reader already
/// open with [`IndexedFileTrajectoryReader::into_structure_set`]. Follow a
/// growing file with [`StructureSet::refresh`], and load geometries into memory
/// with [`StructureSet::load`].
///
/// Cloning a set over a file reads nothing: the clone reads the file
/// object this set reads, shares its frame index, and keeps the set
/// properties, the order, and the distribution, so a host that changes only
/// the distribution or the properties copies the set as cheaply as one in
/// memory. Its pool starts empty. The clone keeps its geometries while the
/// set follows the file: refreshing either adds frames to it alone. A
/// host hands a clone to an analysis and refreshes the set meanwhile. While
/// both share the index, a refresh that finds new frames copies its whole
/// offset table first, one offset per frame, and the first refresh of
/// either reads the start of the file once more, as opening a reader of
/// the pool does.
///
/// A read allocates the geometry's arrays on the thread that reads it, and
/// the pool's readers keep none for their next read, so that the arrays are
/// allocated and freed on one thread wherever the readers move.
///
/// [`TrajectoryIoErrorKind::IndexMismatch`]: crate::traj::TrajectoryIoErrorKind::IndexMismatch
pub struct IndexedFile {
    /// The reader the set opened: the pool's readers are forked from it,
    /// and it alone refreshes. It serves no read, and clones of the set
    /// share it until one of them refreshes.
    reader: Arc<IndexedFileTrajectoryReader>,
    idle: Mutex<Vec<IndexedFileTrajectoryReader>>,
}

/// A storage of the same file object and frame index, with an empty pool;
/// it reads nothing.
impl Clone for IndexedFile {
    fn clone(&self) -> Self {
        Self {
            reader: Arc::clone(&self.reader),
            idle: Mutex::new(Vec::new()),
        }
    }
}

impl fmt::Debug for IndexedFile {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedFile")
            .field("path", &self.path())
            .field("frame_count", &self.reader.frame_count())
            .finish_non_exhaustive()
    }
}

impl IndexedFile {
    /// The trajectory file.
    pub fn path(&self) -> &Path {
        self.reader.path()
    }

    /// The frame index of the geometries: one frame per geometry, sharing the
    /// set's offsets. A host keeps it with [`TrajectoryFrameIndex::to_bytes`]
    /// and opens the file again from it without scanning
    /// ([`open_indexed_structure_set_from_index`]).
    ///
    /// While a host keeps the index, each [`StructureSet::refresh`] that adds
    /// frames copies the offsets, so the host's index stays as it was.
    pub fn frame_index(&self) -> TrajectoryFrameIndex {
        self.reader.to_frame_index()
    }

    /// The fingerprint of the bytes the frame index covers.
    pub fn fingerprint(&self) -> &TrajectoryFileFingerprint {
        self.reader.fingerprint()
    }

    /// The pool's readers waiting for a read.
    fn idle(&self) -> std::sync::MutexGuard<'_, Vec<IndexedFileTrajectoryReader>> {
        self.idle.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// A reader of the pool, or a new one forked from the set's.
    fn take(&self) -> Result<IndexedFileTrajectoryReader, TrajectoryIoError> {
        if let Some(reader) = self.idle().pop() {
            return Ok(reader);
        }
        self.reader.fork()
    }

    /// Returns a reader to the pool.
    fn give_back(&self, reader: IndexedFileTrajectoryReader) {
        self.idle().push(reader);
    }

    /// Reads geometry `index`, which the index holds.
    fn read(&self, index: usize) -> Result<Geometry, TrajectoryIoError> {
        let lost = |error| self.reader.lost_frame(index, error);
        let mut reader = self.take().map_err(lost)?;
        let read = reader.read_sample(index);
        // A failed read puts the stream back, so the reader stays usable.
        self.give_back(reader);
        read.map_err(lost)
    }
}

impl geometry_storage::Storage for IndexedFile {
    fn geometry_count(&self) -> usize {
        self.reader.frame_count().unwrap_or(0)
    }
}

impl GeometryStorage for IndexedFile {}

/// Opens a trajectory file as a [`StructureSet`] over its frames, indexing
/// them through [`open_indexed_trajectory`], which verifies every frame.
///
/// The set's order is a time order of the frames' times and steps, and it
/// has no distribution; give it one with [`StructureSet::set_distribution`].
/// `options` set the codec, its limits, and the indexing's cancellation and
/// progress.
///
/// ```no_run
/// use kekule::structure::StructureSetRead;
/// use kekule::traj::io::open_indexed_structure_set;
///
/// let topology = kekule::smiles::to_topology("CC")?;
/// let set = open_indexed_structure_set("run.xtc", topology, Default::default())?;
/// // Keep the index to open the file again without scanning it.
/// let index = set.storage().frame_index().to_bytes();
/// let geometry = set.read(set.len() - 1)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
///
/// The errors of [`open_indexed_trajectory`].
pub fn open_indexed_structure_set(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryOpenOptions,
) -> Result<StructureSet<IndexedFile>, TrajectoryIoError> {
    Ok(open_indexed_trajectory(path, topology, options)?.into_structure_set())
}

/// Opens a trajectory file as a [`StructureSet`] over its frames from a frame
/// index the host kept, through [`open_indexed_trajectory_from_index`],
/// without scanning the frames it covers; as
/// [`open_indexed_structure_set`] does otherwise. Frames appended to the file
/// after the index are indexed too.
///
/// # Errors
///
/// The errors of [`open_indexed_trajectory_from_index`].
pub fn open_indexed_structure_set_from_index(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    index: TrajectoryFrameIndex,
    options: TrajectoryOpenOptions,
) -> Result<StructureSet<IndexedFile>, TrajectoryIoError> {
    Ok(open_indexed_trajectory_from_index(path, topology, index, options)?.into_structure_set())
}

impl IndexedFileTrajectoryReader {
    /// A [`StructureSet`] over the reader's indexed frames, as
    /// [`open_indexed_structure_set`] opens it, without opening or scanning the
    /// file again: for a host that opened the reader itself, such as to show
    /// the scan's progress or to follow a file still being written.
    ///
    /// The set reads the file object the reader opened. The reader's
    /// sequential position does not matter, since the set reads by index.
    ///
    /// ```no_run
    /// use kekule::structure::StructureSetRead;
    /// use kekule::traj::io::{open_indexed_trajectory, TrajectoryOpenOptions};
    ///
    /// let topology = kekule::smiles::to_topology("CC")?;
    /// let options = TrajectoryOpenOptions::default().with_growing(true);
    /// let reader = open_indexed_trajectory("run.xtc", topology, options)?;
    /// println!("{} bytes still being written", reader.pending_byte_count());
    /// let set = reader.into_structure_set();
    /// let geometry = set.read(0)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_structure_set(self) -> StructureSet<IndexedFile> {
        let topology = self.shared_topology();
        let order = time_order(&self.to_frame_index());
        let storage = IndexedFile {
            reader: Arc::new(self),
            idle: Mutex::new(Vec::new()),
        };
        StructureSet::from_storage(topology, storage, Some(order), None)
    }
}

/// The time order of an index's frames.
fn time_order(index: &TrajectoryFrameIndex) -> Order {
    Order::Time(TimeAxis::from_columns(
        index.times().map(|times| times.value().to_vec()),
        index.steps().map(<[u64]>::to_vec),
    ))
}

impl StructureSet<IndexedFile> {
    /// Indexes the frames appended to the file since the set was opened or
    /// last refreshed, as [`IndexedFileTrajectoryReader::refresh`] does, and
    /// takes them as new geometries; returns how many there are.
    ///
    /// The new geometries extend the time order with their times and steps. A
    /// distribution over every geometry becomes one over the geometries it
    /// covered, so the new geometries are outside it, as
    /// [`StructureSet::push`] leaves an appended geometry; a distribution over a
    /// subset keeps its members, as [`Distribution::extended_to`] gives it.
    /// The pool's readers are kept and read by the new index, which they
    /// share. A set whose clone shares its reader first opens a reader of
    /// its own, as the pool does. Failure changes nothing.
    ///
    /// ```no_run
    /// use kekule::structure::StructureSetRead;
    /// use kekule::traj::io::{open_indexed_structure_set, TrajectoryOpenOptions};
    ///
    /// let topology = kekule::smiles::to_topology("CC")?;
    /// // The trajectory of a running simulation.
    /// let options = TrajectoryOpenOptions::default().with_growing(true);
    /// let mut set = open_indexed_structure_set("run.xtc", topology, options)?;
    /// loop {
    ///     let added = set.refresh()?;
    ///     for index in set.len() - added..set.len() {
    ///         let geometry = set.read(index)?;
    ///         // Show the new geometry.
    ///     }
    ///     # break;
    ///     // Wait for the simulation to write more.
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// The errors of [`IndexedFileTrajectoryReader::refresh`]. A truncated
    /// or rewritten file fails with
    /// [`TrajectoryIoErrorKind::IndexMismatch`]: open it again.
    ///
    /// [`TrajectoryIoErrorKind::IndexMismatch`]: crate::traj::TrajectoryIoErrorKind::IndexMismatch
    pub fn refresh(&mut self) -> Result<usize, TrajectoryIoError> {
        let previous = self.len();
        let storage = self.storage_mut();
        if Arc::get_mut(&mut storage.reader).is_none() {
            // A clone shares the reader, so this set refreshes one of its own.
            storage.reader = Arc::new(storage.reader.fork()?);
        }
        let reader = Arc::get_mut(&mut storage.reader).expect("the set alone holds its reader");
        let idle = storage
            .idle
            .get_mut()
            .unwrap_or_else(PoisonError::into_inner);
        // The pool lets go of the index, so the refresh extends it in place.
        for pooled in idle.iter_mut() {
            pooled.release_index();
        }
        let refreshed = reader.refresh();
        for pooled in idle.iter_mut() {
            pooled.share_index(reader);
        }
        refreshed?;
        let count = reader.frame_count().unwrap_or(previous);
        if count <= previous {
            return Ok(0);
        }
        let order = time_order(&reader.to_frame_index());
        self.grown(previous, Some(order));
        Ok(count - previous)
    }

    /// Reads the geometries `indices` into a set in memory, as
    /// [`StructureSet::from_reader`] does, keeping the set properties.
    ///
    /// # Errors
    ///
    /// As [`StructureSet::from_reader`].
    pub fn load(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<StructureSet, StructureSetReadError> {
        let mut set = StructureSet::from_reader(self, indices)?;
        set.properties_mut().clone_from(self.properties());
        Ok(set)
    }
}

impl StructureSetRead for StructureSet<IndexedFile> {
    fn topology(&self) -> &Arc<Topology> {
        self.topology_handle()
    }

    fn len(&self) -> usize {
        StructureSet::len(self)
    }

    fn order(&self) -> Option<&Order> {
        StructureSet::order(self)
    }

    fn distribution(&self) -> Option<&Distribution> {
        StructureSet::distribution(self)
    }

    fn properties(&self) -> Option<&OwnerProperties> {
        Some(StructureSet::properties(self))
    }

    fn peek(&self, _index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        let len = StructureSet::len(self);
        if index >= len {
            return Err(StructureSetReadError::IndexOutOfRange { index, len });
        }
        self.storage()
            .read(index)
            .map_err(|source| StructureSetReadError::Trajectory {
                geometry: index,
                source,
            })
    }
}
