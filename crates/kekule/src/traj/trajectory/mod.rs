//! Frames, reusable frame buffers, streaming reader/writer contracts, and
//! the trajectory I/O error type.
//!
//! A frame is a [`crate::structure::Geometry`] with its time and MD step:
//! [`FrameView`] borrows one, and [`FrameBuffer`] provides reusable
//! validated storage for streaming decoders. The reader and writer traits
//! publish complete frames transactionally. [`TrajectoryIoError`] reports
//! every codec and file failure.

mod buffer;
mod error;
mod format;
mod frame;
mod stream;

pub use buffer::{FrameBuffer, FrameData};
pub use error::{TrajectoryField, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation};
pub use format::TrajectoryFormat;
pub use frame::FrameView;
pub(crate) use stream::COLUMN_MISMATCH;
pub use stream::{
    MemoryTrajectoryWriter, SeekableTrajectoryReader, TrajectoryReader, TrajectoryWriter,
};
