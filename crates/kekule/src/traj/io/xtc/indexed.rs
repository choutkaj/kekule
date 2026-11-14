use super::*;

/// Fully decoded-and-verified indexed XTC reader.
///
/// Random and sequential reads share decode scratch. A random read restores
/// stream/header state and invalidates cached coordinates, so the next sequential
/// frame is decoded again when it was already pending.
pub struct IndexedXtcReader<R> {
    pub(super) inner: XtcReader<R>,
    pub(super) offsets: Arc<FrameOffsets>,
    /// Bytes after the last indexed frame when the index last grew.
    pub(super) pending_bytes: u64,
}

impl<R> fmt::Debug for IndexedXtcReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedXtcReader")
            .field("reader", &self.inner)
            .field("frame_count", &self.offsets.len())
            .finish_non_exhaustive()
    }
}

impl<R: Read + Seek> IndexedXtcReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order
    /// and indexes it for random access.
    ///
    /// The stream's position at construction is frame 0. Every frame from there to
    /// the end is decoded and verified before the reader is returned, and both
    /// sequential and indexed reads then address frames from that start. Topology,
    /// limits, and diagnostic labels follow [`XtcReader::new`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`XtcReader::new`], the decoding error of the
    /// first invalid frame, [`TrajectoryIoErrorKind::LimitExceeded`] when the
    /// index outgrows the limits, and [`TrajectoryIoErrorKind::Io`] when
    /// reading or seeking fails.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: XtcReadOptions,
    ) -> Result<Self, TrajectoryIoError> {
        Self::build(
            reader,
            topology.into(),
            options,
            false,
            &IndexBuildControl::default(),
        )
    }

    /// Indexes the stream like [`Self::new`], reporting to `control`; a
    /// `growing` stream may end in a partial frame (see [`Self::extend`]).
    pub(crate) fn build(
        reader: R,
        topology: Arc<Topology>,
        options: XtcReadOptions,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        XtcReader::new(reader, topology, options)?.index_frames(growing, control)
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; random reads verify each frame against them.
    pub(crate) fn from_offsets(
        reader: R,
        topology: Arc<Topology>,
        options: XtcReadOptions,
        offsets: Arc<FrameOffsets>,
    ) -> Result<Self, TrajectoryIoError> {
        let inner = XtcReader::new(reader, topology, options)?;
        offsets.verify_start(
            inner.stream_start,
            TrajectoryFormat::Xtc,
            &inner.options.source_label,
        )?;
        Ok(Self {
            inner,
            offsets,
            pending_bytes: 0,
        })
    }

    pub(crate) fn offsets(&self) -> &FrameOffsets {
        &self.offsets
    }

    pub(crate) fn shared_offsets(&self) -> &Arc<FrameOffsets> {
        &self.offsets
    }

    /// Lets go of the index until [`Self::share_index`], so that a reader
    /// extending it does not copy it.
    pub(crate) fn release_index(&mut self) {
        self.offsets = Arc::new(FrameOffsets::empty(self.offsets.end()));
    }

    /// Reads by the index of `other`, a reader of the same file, from now
    /// on.
    pub(crate) fn share_index(&mut self, other: &Self) {
        self.offsets = Arc::clone(&other.offsets);
    }

    /// Another reader of the same frames over `stream`, which must read the
    /// same file from its start; it shares the index and reads the first
    /// frame alone.
    pub(crate) fn fork(&self, stream: R) -> Result<Self, TrajectoryIoError> {
        Self::from_offsets(
            stream,
            Arc::clone(&self.inner.topology),
            self.inner.options.clone(),
            Arc::clone(&self.offsets),
        )
    }

    /// The stream every read goes through.
    pub(crate) fn stream(&self) -> &R {
        &self.inner.adapter.reader
    }

    /// The stream every read goes through.
    pub(crate) fn stream_mut(&mut self) -> &mut R {
        &mut self.inner.adapter.reader
    }

    /// Reads frame `index` as [`SeekableTrajectoryReader::read_frame`]
    /// does and hands over its positions, written once as they were
    /// decoded. The read allocates them on the calling thread and the
    /// reader keeps no buffer for its next read, which may run on another.
    pub(crate) fn read_decoded(
        &mut self,
        index: usize,
    ) -> Result<DecodedGeometry, TrajectoryIoError> {
        let read = self.read_indexed(index);
        // Taken whether the read succeeds or fails, so that the reader keeps
        // no array for a read on another thread.
        let positions = std::mem::take(&mut self.inner.positions);
        let info = read?;
        Ok(DecodedGeometry {
            positions: Quantity::new(positions, NANOMETER),
            cell: info.cell,
            velocities: None,
            forces: None,
            owner: None,
        })
    }

    /// Bytes after the last indexed frame when the index last grew.
    pub(crate) const fn pending_bytes(&self) -> u64 {
        self.pending_bytes
    }

    /// Drops the frames indexed after `mark`; see [`FrameOffsets::mark`].
    pub(crate) fn restore_index(&mut self, mark: (usize, u64), pending_bytes: u64) {
        Arc::make_mut(&mut self.offsets).restore(mark);
        self.pending_bytes = pending_bytes;
    }

    /// Verifies and indexes every complete frame after the last indexed one,
    /// up to the end of the stream, and returns the bytes left after them.
    ///
    /// A `growing` stream may end in a frame still being written: the scan
    /// stops before a frame the stream ends inside, which is not yet an
    /// error. Otherwise such a frame is an error, as is any frame that is
    /// invalid with all its bytes present. Failure leaves the index as it
    /// was; the sequential cursor is kept either way.
    pub(crate) fn extend(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        let mark = self.offsets.mark();
        let saved_offset = self.inner.position(TrajectoryIoOperation::Index)?;
        let saved_cursor = self.inner.frame_cursor;
        let saved_pending = self.inner.pending_info.clone();
        let end = self.offsets.end();
        // A freshly opened reader has already checked frame 0 there.
        if saved_pending.as_ref().is_none_or(|info| info.start != end) {
            self.inner.pending_info = None;
            self.inner.adapter.decoded_start = None;
            self.inner
                .seek_to(SeekFrom::Start(end), TrajectoryIoOperation::Index)?;
        }
        self.inner.frame_cursor = self.offsets.len();
        let scanned = self.scan(growing, control);
        let restored = self
            .inner
            .seek_to(SeekFrom::Start(saved_offset), TrajectoryIoOperation::Index);
        self.inner.frame_cursor = saved_cursor;
        self.inner.pending_info = saved_pending;
        self.inner.adapter.decoded_start = None;
        match scanned.and_then(|pending| restored.map(|_| pending)) {
            Ok(pending) => {
                self.pending_bytes = pending;
                Ok(pending)
            }
            Err(error) => {
                Arc::make_mut(&mut self.offsets).restore(mark);
                Err(error)
            }
        }
    }

    /// Indexes frames from the current position; see [`Self::extend`].
    pub(super) fn scan(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        let label = self.inner.options.source_label.clone();
        control.checkpoint(
            self.offsets.len(),
            self.offsets.end(),
            TrajectoryFormat::Xtc,
            &label,
        )?;
        loop {
            if let Some(limit) =
                projected_index_limit(self.offsets.len(), &self.inner.options.limits)
            {
                let clean_eof = self.inner.pending_info.is_none()
                    && probe_seekable_eof(
                        &mut self.inner.adapter.reader,
                        TrajectoryIoOperation::Index,
                        TrajectoryFormat::Xtc,
                        &label,
                    )?;
                if clean_eof {
                    return Ok(0);
                }
                return Err(resource_error(&label, Some(self.offsets.len()), "")
                    .exceeding("The frame index", limit));
            }
            let offset = match &self.inner.pending_info {
                Some(info) => info.start,
                None => self.inner.position(TrajectoryIoOperation::Index)?,
            };
            let info = match self.inner.parse_next() {
                Ok(Some(info)) => info,
                Ok(None) => return Ok(0),
                Err(error) if growing && error.kind() == TrajectoryIoErrorKind::TruncatedRecord => {
                    let file_end = self
                        .inner
                        .seek_to(SeekFrom::End(0), TrajectoryIoOperation::Index)?;
                    return Ok(file_end.saturating_sub(offset));
                }
                Err(error) => {
                    return Err(frame_offset_context(error, self.inner.frame_cursor, offset))
                }
            };
            let end = self.inner.position(TrajectoryIoOperation::Index)?;
            Arc::make_mut(&mut self.offsets).push(
                offset,
                end,
                info.stamp(),
                &self.inner.options.limits,
                TrajectoryFormat::Xtc,
                &label,
            )?;
            control.checkpoint(self.offsets.len(), end, TrajectoryFormat::Xtc, &label)?;
        }
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.inner.topology()
    }

    pub(crate) fn first_info(&self) -> &XtcFrameInfo {
        self.inner.first_info()
    }
}

impl<R: Read + Seek> TrajectoryReader for IndexedXtcReader<R> {
    fn topology(&self) -> &Topology {
        self.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.inner.shared_topology()
    }

    /// Reads the indexed frames in order and ends after the last of them,
    /// even when the stream holds more.
    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        if self.topology().shares_layout(destination.topology())
            && self.inner.frame_cursor >= self.offsets.len()
        {
            return Ok(false);
        }
        self.inner.read_next(destination)
    }
}

impl<R: Read + Seek> SeekableTrajectoryReader for IndexedXtcReader<R> {
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
        let info = self.read_indexed(index)?;
        self.inner
            .publish(&self.inner.positions, &info, destination)
    }
}

impl<R: Read + Seek> IndexedXtcReader<R> {
    /// Decodes and verifies frame `index` into the reader's buffers, then
    /// puts the stream back where it was.
    pub(super) fn read_indexed(&mut self, index: usize) -> Result<XtcFrameInfo, TrajectoryIoError> {
        let (offset, _) = self
            .offsets
            .frame(index)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        let saved_offset = self
            .inner
            .adapter
            .reader
            .stream_position()
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xtc),
                    &self.inner.options.source_label,
                    error,
                )
            })?;
        let saved_cursor = self.inner.frame_cursor;
        let saved_pending = self.inner.pending_info.clone();
        self.inner
            .adapter
            .reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xtc),
                    &self.inner.options.source_label,
                    error,
                )
            })?;
        self.inner.pending_info = None;
        self.inner.adapter.decoded_start = None;
        self.inner.frame_cursor = index;
        let result = self
            .inner
            .parse_next()
            .map_err(|error| frame_offset_context(error, index, offset))
            .and_then(|info| info.ok_or(TrajectoryIoError::frame_index_out_of_range(index)))
            .and_then(|info| {
                self.offsets.verify_stream_frame_end(
                    index,
                    &mut self.inner.adapter.reader,
                    TrajectoryFormat::Xtc,
                    &self.inner.options.source_label,
                )?;
                self.offsets.verify_stamp(
                    index,
                    info.stamp(),
                    TrajectoryIoOperation::ReadFrame,
                    TrajectoryFormat::Xtc,
                    &self.inner.options.source_label,
                )?;
                Ok(info)
            });
        let restore = self
            .inner
            .adapter
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xtc),
                    &self.inner.options.source_label,
                    error,
                )
            });
        self.inner.frame_cursor = saved_cursor;
        self.inner.pending_info = saved_pending;
        self.inner.adapter.decoded_start = None;
        let info = result?;
        restore?;
        Ok(info)
    }
}
