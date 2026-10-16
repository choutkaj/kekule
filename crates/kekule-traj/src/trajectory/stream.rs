use std::fmt;
use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::structure::{Positions, StructureError, Trajectory, TrajectoryFrameView};
use kekule::topology::{InstanceAtomId, Topology};
use kekule::units::Quantity;

use super::buffer::FrameBuffer;
use super::TrajectoryIoError;

/// Sequential reader that publishes complete frames into reusable storage.
///
/// `read_next` returns `Ok(false)` only at clean end-of-stream. Implementations
/// validate a frame completely before replacing the destination buffer.
pub trait TrajectoryReader {
    fn topology(&self) -> &Topology;

    fn shared_topology(&self) -> Arc<Topology>;

    /// Creates reusable frame storage sharing this reader's exact topology.
    fn frame_buffer(&self) -> FrameBuffer {
        FrameBuffer::new(self.shared_topology())
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError>;
}

/// Trajectory reader with random frame access.
///
/// Implementations may build an index eagerly; callers should inspect
/// format-specific metadata when index construction cost matters.
/// Random reads preserve the sequential cursor: the next `read_next` call reads
/// the same frame it would have read before `read_frame`. Invalid indices,
/// topology mismatches, and recoverable decoding errors also preserve the cursor.
/// An I/O failure while seeking or restoring the stream can leave its position
/// unspecified; callers should reopen the reader before continuing in that case.
pub trait SeekableTrajectoryReader: TrajectoryReader {
    fn frame_count(&self) -> Option<u64>;

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError>;
}

/// Streaming writer for topology-bound frame views.
///
/// Writers reject topology mismatches and unsupported frame state rather than
/// silently dropping data. A rejected frame writes nothing, so the writer
/// remains usable. A failure while a frame's bytes are being written, such as
/// an I/O error, may leave a partial frame in the stream; the writer then
/// refuses every further frame and finishing, rather than extend or report a
/// corrupt stream as complete.
pub trait TrajectoryWriter {
    fn topology(&self) -> &Topology;

    fn shared_topology(&self) -> Arc<Topology>;

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryIoError>;
}

/// Sequential and seekable reader over an in-memory trajectory.
pub struct MemoryTrajectoryReader<'a> {
    trajectory: &'a Trajectory,
    cursor: usize,
}

impl fmt::Debug for MemoryTrajectoryReader<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MemoryTrajectoryReader")
            .field("frame_count", &self.trajectory.len())
            .field("cursor", &self.cursor)
            .finish_non_exhaustive()
    }
}

impl<'a> MemoryTrajectoryReader<'a> {
    pub fn new(trajectory: &'a Trajectory) -> Self {
        Self {
            trajectory,
            cursor: 0,
        }
    }
}

impl TrajectoryReader for MemoryTrajectoryReader<'_> {
    fn topology(&self) -> &Topology {
        self.trajectory.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.trajectory.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        let Some(frame) = self.trajectory.get(self.cursor) else {
            return Ok(false);
        };
        destination.copy_from(frame)?;
        self.cursor += 1;
        Ok(true)
    }
}

impl SeekableTrajectoryReader for MemoryTrajectoryReader<'_> {
    fn frame_count(&self) -> Option<u64> {
        u64::try_from(self.trajectory.len()).ok()
    }

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        let index =
            usize::try_from(index).map_err(|_| TrajectoryIoError::FrameIndexOutOfRange(index))?;
        let frame = self
            .trajectory
            .get(index)
            .ok_or(TrajectoryIoError::FrameIndexOutOfRange(index as u64))?;
        destination.copy_from(frame)?;
        Ok(())
    }
}

/// In-memory writer that preserves every core frame field.
///
/// Useful as a [`TrajectoryWriter`] destination for streaming pipelines. To load
/// a complete file directly, use [`crate::io::read_trajectory`].
pub struct MemoryTrajectoryWriter {
    trajectory: Trajectory,
}

impl fmt::Debug for MemoryTrajectoryWriter {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MemoryTrajectoryWriter")
            .field("frame_count", &self.trajectory.len())
            .finish_non_exhaustive()
    }
}

impl MemoryTrajectoryWriter {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            trajectory: Trajectory::new(topology),
        }
    }

    pub fn into_trajectory(self) -> Trajectory {
        self.trajectory
    }
}

impl TrajectoryWriter for MemoryTrajectoryWriter {
    fn topology(&self) -> &Topology {
        self.trajectory.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.trajectory.shared_topology()
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryIoError> {
        if !self.trajectory.topology().shares_layout(frame.topology()) {
            return Err(TrajectoryIoError::TopologyMismatch);
        }
        Ok(self.trajectory.push(frame.payload().clone())?)
    }
}

/// Checks an independently supplied semantic atom sequence against topology order.
///
/// This compares the supplied IDs only; it does not inspect a coordinate file or
/// infer atom identity from matching counts. Readers interpret coordinate index
/// `i` as topology dense atom index `i` without requiring this optional check.
pub fn validate_atom_order(
    topology: &Topology,
    atom_order: &[InstanceAtomId],
) -> Result<(), TrajectoryIoError> {
    if topology.atom_ids() != atom_order {
        return Err(TrajectoryIoError::AtomOrderMismatch);
    }
    Ok(())
}

/// Reference reader for a topology-free coordinate source.
///
/// Construction collects the source and validates every frame's units and atom
/// count. Coordinates follow the supplied topology's dense atom order.
pub struct CoordinateFrameReader {
    topology: Arc<Topology>,
    frames: Vec<Vec<Point3>>,
    cursor: usize,
}

impl fmt::Debug for CoordinateFrameReader {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CoordinateFrameReader")
            .field("atom_count", &self.topology.atom_count())
            .field("frame_count", &self.frames.len())
            .field("cursor", &self.cursor)
            .finish_non_exhaustive()
    }
}

impl CoordinateFrameReader {
    pub fn new(
        topology: impl Into<Arc<Topology>>,
        frames: impl IntoIterator<Item = Quantity<Vec<Point3>>>,
    ) -> Result<Self, TrajectoryIoError> {
        let topology = topology.into();
        let frames = frames
            .into_iter()
            .map(|frame| {
                let positions = Positions::new(frame)?;
                if positions.len() != topology.atom_count() {
                    return Err(StructureError::AtomCountMismatch {
                        expected: topology.atom_count(),
                        actual: positions.len(),
                    }
                    .into());
                }
                Ok::<_, TrajectoryIoError>(positions.values().value().to_vec())
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology,
            frames,
            cursor: 0,
        })
    }
}

impl TrajectoryReader for CoordinateFrameReader {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        if !self.topology.shares_layout(destination.topology()) {
            return Err(TrajectoryIoError::TopologyMismatch);
        }
        let Some(frame) = self.frames.get(self.cursor) else {
            return Ok(false);
        };
        destination.replace_from_data(crate::FrameBufferData::new(Quantity::new(
            frame.as_slice(),
            kekule::units::CANONICAL_LENGTH_UNIT,
        )))?;
        self.cursor += 1;
        Ok(true)
    }
}
