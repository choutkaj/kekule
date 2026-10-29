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
    /// The topology every frame of this reader belongs to.
    fn topology(&self) -> &Topology;

    /// A shared handle to [`Self::topology`].
    fn shared_topology(&self) -> Arc<Topology>;

    /// Creates reusable frame storage sharing this reader's exact topology.
    fn frame_buffer(&self) -> FrameBuffer {
        FrameBuffer::new(self.shared_topology())
    }

    /// Reads the next frame into `destination`.
    ///
    /// Returns `Ok(true)` when a frame was published and `Ok(false)` at a
    /// clean end of stream, leaving `destination` unchanged.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::TopologyMismatch`] when `destination`
    /// does not share the reader's layout, and a decoding, validation or
    /// [`TrajectoryIoErrorKind::Io`] error when the frame cannot be read; on
    /// error `destination` is unchanged.
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
    /// The number of frames, or `None` when it is unknown without reading the
    /// whole stream.
    fn frame_count(&self) -> Option<usize>;

    /// Reads the zero-based frame `index` into `destination` without moving
    /// the sequential cursor.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::FrameIndexOutOfRange`] when no frame
    /// has `index`, [`TrajectoryIoErrorKind::TopologyMismatch`] when
    /// `destination` does not share the reader's layout, and a decoding or
    /// [`TrajectoryIoErrorKind::Io`] error when the frame cannot be read; on
    /// error `destination` is unchanged.
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
/// # Writing a file another process reads
///
/// The file codec writers of [`crate::io`] append each frame's bytes after
/// the previous frame's, in order, and never rewrite a byte of a frame once
/// it is written; DCD writes its header first and rewrites only the
/// header's frame count. After [`Self::flush`] the stream therefore holds a
/// complete trajectory of every frame written so far, which a reader in
/// another process sees as whole frames; while a frame is being written, or
/// between flushes of a buffered stream, the file may end in part of a
/// frame, which an indexed reader opened with
/// [`TrajectoryOpenOptions::with_growing`](crate::io::TrajectoryOpenOptions::with_growing)
/// takes as not yet written. A simulation whose trajectory is viewed while
/// it runs writes through a stream writer over the destination file
/// ([`crate::io::create_trajectory_stream_writer`]) and flushes after each
/// frame, or as often as readers should see frames; a path writer
/// ([`crate::io::create_trajectory_writer`]) publishes its file only when it
/// finishes.
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
/// - [`Self::flush`] leaves the destination holding every frame written,
///   whole; writers without a destination of their own may keep the default,
///   which does nothing.
pub trait TrajectoryWriter {
    /// The topology this writer was opened for.
    fn topology(&self) -> &Topology;

    /// A shared handle to [`Self::topology`].
    fn shared_topology(&self) -> Arc<Topology>;

    /// Appends one frame.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::TopologyMismatch`] when the frame's
    /// topology does not share the writer's layout,
    /// [`TrajectoryIoErrorKind::UnsupportedField`] when the format cannot store
    /// a field the frame carries, another validation kind when the frame
    /// cannot be encoded, and [`TrajectoryIoErrorKind::Io`] when writing fails.
    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryIoError>;

    /// Makes every frame written so far reach the destination whole, so
    /// that a reader of it sees them; see the trait documentation. More
    /// frames can be written afterwards. The default does nothing.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidFrame`] after an earlier
    /// write left a partial frame, and [`TrajectoryIoErrorKind::Io`] when
    /// the stream fails.
    fn flush(&mut self) -> Result<(), TrajectoryIoError> {
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
    /// An empty writer collecting frames of `topology`.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        let topology = topology.into();
        Self {
            trajectory: Trajectory::new(topology),
        }
    }

    /// The frames written so far, in write order.
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
