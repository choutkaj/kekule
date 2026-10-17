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
///
/// # Implementing this trait
///
/// This trait is an extension point for codecs outside this crate. An
/// implementation must keep this contract:
///
/// - [`Self::topology`] and [`Self::shared_topology`] return the same snapshot
///   for the reader's whole life, and every published frame belongs to it.
/// - [`Self::read_next`] rejects a destination whose topology does not share
///   that layout, and on any error leaves the destination unchanged.
/// - A frame is published only after every field it carries has been decoded
///   and validated; coordinates are stored in canonical units.
/// - `Ok(false)` means a clean end of stream with no partial frame; a
///   truncated frame is an error.
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
///
/// # Implementing this trait
///
/// Besides the [`TrajectoryReader`] contract, an implementation must:
///
/// - number frames from zero in stream order, matching the order
///   [`TrajectoryReader::read_next`] publishes them;
/// - return the same frame from [`Self::read_frame`] that sequential reading
///   would publish at that index;
/// - preserve the sequential cursor across random reads, as described above;
/// - report [`Self::frame_count`] as `None` only when the count is unknown
///   without reading the whole stream.
pub trait SeekableTrajectoryReader: TrajectoryReader {
    fn frame_count(&self) -> Option<usize>;

    fn read_frame(
        &mut self,
        index: usize,
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
///
/// # Implementing this trait
///
/// This trait is an extension point for codecs outside this crate. An
/// implementation must keep this contract:
///
/// - [`Self::topology`] and [`Self::shared_topology`] return the snapshot the
///   writer was opened for, and [`Self::write_frame`] rejects a frame whose
///   topology does not share its layout.
/// - A frame field the format cannot represent is an error, never silently
///   dropped or approximated beyond the format's documented precision.
/// - A rejected frame writes nothing; only a failure part-way through writing
///   may poison the writer, which then refuses every further frame.
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
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        let topology = topology.into();
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
