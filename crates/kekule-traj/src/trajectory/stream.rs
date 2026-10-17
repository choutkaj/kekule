use std::fmt;
use std::sync::Arc;

use kekule::structure::{Trajectory, TrajectoryFrameView};
use kekule::topology::Topology;

use super::buffer::FrameBuffer;
use super::{TrajectoryIoError, TrajectoryIoErrorKind};

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
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        Ok(self.trajectory.push(frame.payload().clone())?)
    }
}

#[cfg(test)]
mod tests;
