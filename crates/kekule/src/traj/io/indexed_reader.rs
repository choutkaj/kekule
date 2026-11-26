//! The indexed reader every codec shares.
//!
//! A codec's indexed reader wraps an [`IndexedReader`] over the codec's
//! sequential reader, which implements [`IndexedCodec`]. The index lists the
//! verified byte range, time, and step of every frame from the stream's
//! frame 0. A scan or a random read moves the stream and puts back the
//! stream position, frame cursor, and decode state of the sequential reads
//! afterwards, whether it succeeds or not, and a failed scan leaves the index
//! as it was.

use crate::traj::trajectory::reason::Part;
use std::fmt;
use std::io::{Read, Seek, SeekFrom};
use std::sync::Arc;

use crate::diagnostics::LimitExceeded;
use crate::topology::Topology;
use crate::traj::io::DecodedGeometry;
use crate::traj::{
    FrameBuffer, SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoError,
    TrajectoryIoErrorKind, TrajectoryIoOperation, TrajectoryReader,
};

use super::codec::limit_error;
use super::options::ReadSource;
use super::{
    frame_offset_context, io_context, probe_seekable_eof, projected_index_limit, FrameOffsets,
    FrameStamp, IndexBuildControl,
};

/// What an index records about its file beyond the frames.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct IndexFacts {
    /// The DCD header's frame count `NSET` when the index last grew.
    pub(crate) declared_frames: Option<usize>,
    /// Whether the indexed TRR frames mix scalar widths.
    pub(crate) precision_mixed: bool,
    /// The TRR frames without positions before indexed frames, which the
    /// scans skipped.
    pub(crate) skipped_frames: usize,
    /// The TRR frames without positions after the last indexed frame, which
    /// the last scan skipped.
    pub(crate) trailing_frames: usize,
}

/// A sequential reader that an [`IndexedReader`] scans and seeks.
///
/// `State` is the decode state beyond the stream position and frame cursor,
/// such as a frame header read ahead, that a scan or random read saves and
/// restores. The hooks with defaults suit a codec that reads no header ahead
/// and records no [`IndexFacts`].
pub(crate) trait IndexedCodec: TrajectoryReader {
    type Stream: Read + Seek;
    /// What a frame decoded into the reader's buffers carries besides them.
    type Frame;
    type State;
    const FORMAT: TrajectoryFormat;
    /// The operation that a random read's seeks report.
    const SEEK_OPERATION: TrajectoryIoOperation;

    fn source(&self) -> &ReadSource;
    fn stream(&self) -> &Self::Stream;
    fn stream_and_source(&mut self) -> (&mut Self::Stream, &ReadSource);
    fn frame_cursor(&self) -> usize;
    fn set_frame_cursor(&mut self, frame: usize);
    fn save_state(&self) -> Self::State;
    fn restore_state(&mut self, state: Self::State);

    fn index_facts(&self) -> IndexFacts {
        IndexFacts::default()
    }

    fn set_index_facts(&mut self, _facts: IndexFacts) {}

    /// Whether the next frame's header was already read.
    fn frame_pending(&self) -> bool {
        false
    }

    /// Prepares to scan from `end`, where the last indexed frame ends.
    fn start_scan(&mut self, end: u64, growing: bool) -> Result<(), TrajectoryIoError>;

    /// Records in the index facts what a scan found.
    fn end_scan(&mut self) {}

    /// Ends the scan before the frame after `indexed` ones, returning the
    /// pending bytes, or lets it go on with `None`.
    fn before_frame(
        &mut self,
        _indexed: usize,
        _growing: bool,
    ) -> Result<Option<u64>, TrajectoryIoError> {
        Ok(None)
    }

    /// Ends a scan whose index has no room for another frame: cleanly at the
    /// end of the stream, otherwise with an error.
    fn index_full(
        &mut self,
        indexed: usize,
        _growing: bool,
        limit: LimitExceeded,
    ) -> Result<u64, TrajectoryIoError> {
        let clean_eof = !self.frame_pending() && {
            let (stream, source) = self.stream_and_source();
            probe_seekable_eof(
                stream,
                TrajectoryIoOperation::Index,
                Self::FORMAT,
                &source.label,
            )?
        };
        if clean_eof {
            return Ok(0);
        }
        Err(limit_error(
            Self::FORMAT,
            TrajectoryIoOperation::ReadFrame,
            &self.source().label,
            Some(indexed),
            None,
        )
        .exceeding(Part::FrameIndex, limit))
    }

    /// Where the next frame starts.
    fn frame_offset(&mut self) -> Result<u64, TrajectoryIoError> {
        position(self, TrajectoryIoOperation::Index)
    }

    /// Decodes the frame after `indexed` ones, which starts at `offset`, and
    /// returns its stamp, or `None` at the end of the stream.
    fn index_next(
        &mut self,
        indexed: usize,
        offset: u64,
        growing: bool,
    ) -> Result<Option<FrameStamp>, TrajectoryIoError>;

    /// Prepares to decode the frame at `offset`, where the stream now is.
    fn start_read(&mut self, _offset: u64) {}

    /// Decodes the frame at the stream position, or returns `None` at the
    /// end of the stream.
    fn read_frame(&mut self) -> Result<Option<Self::Frame>, TrajectoryIoError>;

    /// Adds the frame and its offset to an error of a random read.
    fn frame_context(error: TrajectoryIoError, frame: usize, offset: u64) -> TrajectoryIoError {
        frame_offset_context(error, frame, offset)
    }

    /// Checks a randomly read frame against the index beyond its extent.
    fn check_frame(
        &self,
        _frame: &Self::Frame,
        _index: usize,
        _offsets: &FrameOffsets,
    ) -> Result<(), TrajectoryIoError> {
        Ok(())
    }

    fn publish_frame(
        &self,
        frame: Self::Frame,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError>;

    /// Hands over the decoded arrays of a random read. They are taken
    /// whether the read succeeded or failed, so that the reader keeps no
    /// array for a read on another thread.
    fn take_decoded(
        &mut self,
        read: Result<Self::Frame, TrajectoryIoError>,
    ) -> Result<DecodedGeometry, TrajectoryIoError>;
}

fn position<C: IndexedCodec + ?Sized>(
    codec: &mut C,
    operation: TrajectoryIoOperation,
) -> Result<u64, TrajectoryIoError> {
    let (stream, source) = codec.stream_and_source();
    stream
        .stream_position()
        .map_err(|error| io_context(operation, Some(C::FORMAT), &source.label, error))
}

fn seek<C: IndexedCodec>(
    codec: &mut C,
    to: SeekFrom,
    operation: TrajectoryIoOperation,
) -> Result<u64, TrajectoryIoError> {
    let (stream, source) = codec.stream_and_source();
    stream
        .seek(to)
        .map_err(|error| io_context(operation, Some(C::FORMAT), &source.label, error))
}

/// A codec reader with verified offsets for its frames.
///
/// Sequential reads go through the codec reader and end after the last
/// indexed frame, even when the stream holds more. Random reads verify each
/// frame against the index.
pub(crate) struct IndexedReader<C> {
    pub(crate) codec: C,
    offsets: Arc<FrameOffsets>,
    /// Bytes after the last indexed frame when the index last grew.
    pending_bytes: u64,
}

impl<C: fmt::Debug> IndexedReader<C> {
    pub(crate) fn fmt_named(&self, name: &str, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct(name)
            .field("reader", &self.codec)
            .field("frame_count", &self.offsets.len())
            .finish_non_exhaustive()
    }
}

impl<C: IndexedCodec> IndexedReader<C> {
    /// Reads by previously verified `offsets`.
    pub(crate) fn new(codec: C, offsets: Arc<FrameOffsets>) -> Self {
        Self {
            codec,
            offsets,
            pending_bytes: 0,
        }
    }

    /// Indexes the frames from `start`, frame 0's offset, reporting to
    /// `control`; see [`IndexedSource::extend`].
    pub(crate) fn build(
        codec: C,
        start: u64,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        let mut indexed = Self::new(codec, Arc::new(FrameOffsets::empty(start)));
        indexed.extend(growing, control)?;
        Ok(indexed)
    }

    /// Reads by the index and index facts of `other`, a reader of the same
    /// file, from now on.
    pub(crate) fn share_index(&mut self, other: &Self) {
        self.offsets = Arc::clone(&other.offsets);
        self.codec.set_index_facts(other.codec.index_facts());
    }

    /// Indexes frames from the current position.
    fn scan(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        let label = self.codec.source().label.clone();
        control.checkpoint(self.offsets.len(), self.offsets.end(), C::FORMAT, &label)?;
        loop {
            let indexed = self.offsets.len();
            if let Some(pending) = self.codec.before_frame(indexed, growing)? {
                return Ok(pending);
            }
            if let Some(limit) = projected_index_limit(indexed, &self.codec.source().limits) {
                return self.codec.index_full(indexed, growing, limit);
            }
            let offset = self.codec.frame_offset()?;
            let stamp = match self.codec.index_next(indexed, offset, growing) {
                Ok(Some(stamp)) => stamp,
                Ok(None) => return Ok(0),
                Err(error) if growing && error.kind() == TrajectoryIoErrorKind::TruncatedRecord => {
                    let file_end = seek(
                        &mut self.codec,
                        SeekFrom::End(0),
                        TrajectoryIoOperation::Index,
                    )?;
                    return Ok(file_end.saturating_sub(offset));
                }
                Err(error) => return Err(error),
            };
            let end = position(&mut self.codec, TrajectoryIoOperation::Index)?;
            Arc::make_mut(&mut self.offsets).push(
                offset,
                end,
                stamp,
                &self.codec.source().limits,
                C::FORMAT,
                &label,
            )?;
            control.checkpoint(self.offsets.len(), end, C::FORMAT, &label)?;
        }
    }

    /// Decodes and verifies frame `index` into the codec reader's buffers.
    fn read_indexed(&mut self, index: usize) -> Result<C::Frame, TrajectoryIoError> {
        let (offset, _) = self
            .offsets
            .frame(index)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        let saved_offset = position(&mut self.codec, C::SEEK_OPERATION)?;
        let saved_cursor = self.codec.frame_cursor();
        let saved_state = self.codec.save_state();
        seek(&mut self.codec, SeekFrom::Start(offset), C::SEEK_OPERATION)?;
        self.codec.start_read(offset);
        self.codec.set_frame_cursor(index);
        let result = self
            .codec
            .read_frame()
            .map_err(|error| C::frame_context(error, index, offset))
            .and_then(|frame| frame.ok_or(TrajectoryIoError::frame_index_out_of_range(index)))
            .and_then(|frame| {
                let (stream, source) = self.codec.stream_and_source();
                self.offsets
                    .verify_stream_frame_end(index, stream, C::FORMAT, &source.label)?;
                self.codec.check_frame(&frame, index, &self.offsets)?;
                Ok(frame)
            });
        let restore = seek(
            &mut self.codec,
            SeekFrom::Start(saved_offset),
            C::SEEK_OPERATION,
        );
        self.codec.set_frame_cursor(saved_cursor);
        self.codec.restore_state(saved_state);
        let frame = result?;
        restore?;
        Ok(frame)
    }
}

/// What the path-backed indexed reader uses of a codec's indexed reader over
/// a stream `S`, for every format alike.
pub(crate) trait IndexedSource<S>: SeekableTrajectoryReader {
    fn offsets(&self) -> &FrameOffsets;
    fn shared_offsets(&self) -> &Arc<FrameOffsets>;

    /// Lets go of the index until [`IndexedReader::share_index`], so that a
    /// reader extending it does not copy it.
    fn release_index(&mut self);

    fn facts(&self) -> IndexFacts;

    /// Bytes after the last indexed frame when the index last grew.
    fn pending_bytes(&self) -> u64;

    /// Drops the frames indexed after `mark`, with the index facts and
    /// pending bytes the reader had then; see [`FrameOffsets::mark`].
    fn restore_index(&mut self, mark: (usize, u64), facts: IndexFacts, pending_bytes: u64);

    /// Verifies and indexes every complete frame after the last indexed
    /// one, up to the end of the stream, and returns the bytes left after
    /// them.
    ///
    /// A `growing` stream may end in a frame still being written: the scan
    /// stops before a frame the stream ends inside, which is not yet an
    /// error. Otherwise such a frame is an error, as is any frame that is
    /// invalid with all its bytes present. Failure leaves the index as it
    /// was; the sequential cursor is kept either way.
    fn extend(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError>;

    /// Reads frame `index` as [`SeekableTrajectoryReader::read_frame`]
    /// does and hands over its arrays, each written once as it was decoded.
    /// The read allocates them on the calling thread and the reader keeps no
    /// buffer for its next read, which may run on another.
    fn read_decoded(&mut self, index: usize) -> Result<DecodedGeometry, TrajectoryIoError>;

    /// The stream every read goes through.
    fn stream(&self) -> &S;
    fn stream_mut(&mut self) -> &mut S;
}

impl<C: IndexedCodec> IndexedSource<C::Stream> for IndexedReader<C> {
    fn offsets(&self) -> &FrameOffsets {
        &self.offsets
    }

    fn shared_offsets(&self) -> &Arc<FrameOffsets> {
        &self.offsets
    }

    fn release_index(&mut self) {
        self.offsets = Arc::new(FrameOffsets::empty(self.offsets.end()));
    }

    fn facts(&self) -> IndexFacts {
        self.codec.index_facts()
    }

    fn pending_bytes(&self) -> u64 {
        self.pending_bytes
    }

    fn restore_index(&mut self, mark: (usize, u64), facts: IndexFacts, pending_bytes: u64) {
        Arc::make_mut(&mut self.offsets).restore(mark);
        self.codec.set_index_facts(facts);
        self.pending_bytes = pending_bytes;
    }

    fn extend(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        let mark = self.offsets.mark();
        let facts = self.codec.index_facts();
        let saved_offset = position(&mut self.codec, TrajectoryIoOperation::Index)?;
        let saved_cursor = self.codec.frame_cursor();
        let saved_state = self.codec.save_state();
        let scanned = self
            .codec
            .start_scan(self.offsets.end(), growing)
            .and_then(|()| {
                self.codec.set_frame_cursor(self.offsets.len());
                self.scan(growing, control)
            });
        self.codec.end_scan();
        let restored = seek(
            &mut self.codec,
            SeekFrom::Start(saved_offset),
            TrajectoryIoOperation::Index,
        );
        self.codec.set_frame_cursor(saved_cursor);
        self.codec.restore_state(saved_state);
        match scanned.and_then(|pending| restored.map(|_| pending)) {
            Ok(pending) => {
                self.pending_bytes = pending;
                Ok(pending)
            }
            Err(error) => {
                Arc::make_mut(&mut self.offsets).restore(mark);
                self.codec.set_index_facts(facts);
                Err(error)
            }
        }
    }

    fn read_decoded(&mut self, index: usize) -> Result<DecodedGeometry, TrajectoryIoError> {
        let read = self.read_indexed(index);
        self.codec.take_decoded(read)
    }

    fn stream(&self) -> &C::Stream {
        self.codec.stream()
    }

    fn stream_mut(&mut self) -> &mut C::Stream {
        self.codec.stream_and_source().0
    }
}

impl<C: IndexedCodec> TrajectoryReader for IndexedReader<C> {
    fn topology(&self) -> &Topology {
        self.codec.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.codec.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        if self.topology().shares_layout(destination.topology())
            && self.codec.frame_cursor() >= self.offsets.len()
        {
            return Ok(false);
        }
        self.codec.read_next(destination)
    }
}

impl<C: IndexedCodec> SeekableTrajectoryReader for IndexedReader<C> {
    fn frame_count(&self) -> Option<usize> {
        Some(self.offsets.len())
    }

    fn read_frame(
        &mut self,
        index: usize,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        if !self.topology().shares_layout(destination.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        let frame = self.read_indexed(index)?;
        self.codec.publish_frame(frame, destination)
    }
}

/// Implements `Debug`, [`TrajectoryReader`], and [`SeekableTrajectoryReader`]
/// for a codec's indexed reader `$name<R>` by its `inner` [`IndexedReader`].
macro_rules! indexed_reader_traits {
    ($name:ident, $($bound:tt)+) => {
        impl<R> std::fmt::Debug for $name<R> {
            fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                self.inner.fmt_named(stringify!($name), formatter)
            }
        }

        impl<R: $($bound)+> crate::traj::TrajectoryReader for $name<R> {
            fn topology(&self) -> &crate::topology::Topology {
                crate::traj::TrajectoryReader::topology(&self.inner)
            }

            fn shared_topology(&self) -> std::sync::Arc<crate::topology::Topology> {
                crate::traj::TrajectoryReader::shared_topology(&self.inner)
            }

            /// Reads the indexed frames in order and ends after the last of
            /// them, even when the stream holds more.
            fn read_next(
                &mut self,
                destination: &mut crate::traj::FrameBuffer,
            ) -> Result<bool, crate::traj::TrajectoryIoError> {
                crate::traj::TrajectoryReader::read_next(&mut self.inner, destination)
            }
        }

        impl<R: $($bound)+> crate::traj::SeekableTrajectoryReader for $name<R> {
            fn frame_count(&self) -> Option<usize> {
                crate::traj::SeekableTrajectoryReader::frame_count(&self.inner)
            }

            fn read_frame(
                &mut self,
                index: usize,
                destination: &mut crate::traj::FrameBuffer,
            ) -> Result<(), crate::traj::TrajectoryIoError> {
                crate::traj::SeekableTrajectoryReader::read_frame(
                    &mut self.inner,
                    index,
                    destination,
                )
            }
        }
    };
}

pub(crate) use indexed_reader_traits;
