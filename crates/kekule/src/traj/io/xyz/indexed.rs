use super::*;

/// Fully scanned XYZ reader with verified frame offsets.
pub struct IndexedXyzReader<R> {
    pub(super) reader: XyzReader<R>,
    pub(super) offsets: Arc<FrameOffsets>,
    /// Bytes after the last indexed frame when the index last grew.
    pub(super) pending_bytes: u64,
}

impl<R> fmt::Debug for IndexedXyzReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedXyzReader")
            .field("reader", &self.reader)
            .field("frame_count", &self.offsets.len())
            .finish_non_exhaustive()
    }
}

impl<R: BufRead + Seek> IndexedXyzReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order
    /// and indexes it for random access.
    ///
    /// The stream's position at construction is frame 0. Every frame from there to
    /// the end is decoded and verified before the reader is returned, and both
    /// sequential and indexed reads then address frames from that start. Topology,
    /// limits, and diagnostic labels follow [`XyzReader::new`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`XyzReader::new`], the decoding error of the
    /// first invalid frame, [`TrajectoryIoErrorKind::LimitExceeded`] when the
    /// index outgrows the limits, and [`TrajectoryIoErrorKind::Io`] when
    /// reading or seeking fails.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: XyzReadOptions,
    ) -> Result<Self, TrajectoryIoError> {
        Self::build(
            reader,
            topology.into(),
            options,
            false,
            &IndexBuildControl::default(),
        )
    }

    /// Indexes the stream like [`Self::new`], reporting to `control`. A
    /// `growing` stream may end in a partial frame, and only lines with a
    /// line break count (see [`Self::extend`]).
    pub(crate) fn build(
        mut reader: R,
        topology: Arc<Topology>,
        options: XyzReadOptions,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        let stream_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &options.source_label,
                error,
            )
        })?;
        XyzReader::new(reader, topology, options)?.index_frames(stream_start, growing, control)
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; frame 0 is verified here, and random reads verify each
    /// frame against the offsets. A `growing` stream is extended as in
    /// [`Self::build`].
    pub(crate) fn from_offsets(
        mut reader: R,
        topology: Arc<Topology>,
        options: XyzReadOptions,
        offsets: Arc<FrameOffsets>,
        growing: bool,
    ) -> Result<Self, TrajectoryIoError> {
        let stream_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &options.source_label,
                error,
            )
        })?;
        let mut reader = XyzReader::new(reader, topology, options)?;
        reader.growing = growing;
        offsets.verify_start(
            stream_start,
            TrajectoryFormat::Xyz,
            &reader.options.source_label,
        )?;
        reader.validate_first_frame()?;
        Ok(Self {
            reader,
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
            Arc::clone(&self.reader.topology),
            self.reader.options.clone(),
            Arc::clone(&self.offsets),
            self.reader.growing,
        )
    }

    /// The stream every read goes through.
    pub(crate) fn stream(&self) -> &R {
        &self.reader.reader
    }

    /// The stream every read goes through.
    pub(crate) fn stream_mut(&mut self) -> &mut R {
        &mut self.reader.reader
    }

    /// Reads frame `index` as [`SeekableTrajectoryReader::read_frame`]
    /// does and hands over its positions, written once as they were parsed.
    /// The read allocates them on the calling thread and the reader keeps
    /// no buffer for its next read, which may run on another.
    pub(crate) fn read_decoded(
        &mut self,
        index: usize,
    ) -> Result<DecodedGeometry, TrajectoryIoError> {
        let read = self.read_indexed(index);
        // Taken whether the read succeeds or fails, so that the reader keeps
        // no array for a read on another thread.
        let positions = std::mem::take(&mut self.reader.positions);
        read?;
        Ok(DecodedGeometry {
            positions: Quantity::new(positions, self.reader.options.length_unit),
            cell: None,
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
    /// A growing stream (see [`Self::build`]) may end in a frame still being
    /// written, such as one whose last line has no line break yet: the scan
    /// stops before it, which is not yet an error. Otherwise such a frame is
    /// an error, as is any frame that is invalid with all its lines present.
    /// Failure leaves the index as it was; the sequential cursor is kept
    /// either way.
    pub(crate) fn extend(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        self.reader.growing = growing;
        let mark = self.offsets.mark();
        let saved_offset = self.reader.position(TrajectoryIoOperation::Index)?;
        let saved_cursor = self.reader.frame_cursor;
        let scanned = self
            .reader
            .seek_to(
                SeekFrom::Start(self.offsets.end()),
                TrajectoryIoOperation::Index,
            )
            .and_then(|_| {
                self.reader.frame_cursor = self.offsets.len();
                self.scan(control)
            });
        let restored = self
            .reader
            .seek_to(SeekFrom::Start(saved_offset), TrajectoryIoOperation::Index);
        self.reader.frame_cursor = saved_cursor;
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
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        let label = self.reader.options.source_label.clone();
        control.checkpoint(
            self.offsets.len(),
            self.offsets.end(),
            TrajectoryFormat::Xyz,
            &label,
        )?;
        loop {
            if let Some(limit) =
                projected_index_limit(self.offsets.len(), &self.reader.options.limits)
            {
                if self
                    .reader
                    .reader
                    .fill_buf()
                    .map_err(|error| {
                        io_context(TrajectoryIoOperation::Index, Some(FORMAT), &label, error)
                    })?
                    .is_empty()
                {
                    return Ok(0);
                }
                return Err(codec_context(
                    TrajectoryIoErrorKind::LimitExceeded,
                    TrajectoryIoOperation::Index,
                    Some(FORMAT),
                    &label,
                    "",
                )
                .exceeding("The frame index", limit));
            }
            let offset = self.reader.position(TrajectoryIoOperation::Index)?;
            match self.reader.parse_next(false) {
                Ok(true) => {}
                Ok(false) => return Ok(0),
                Err(error)
                    if self.reader.growing
                        && error.kind() == TrajectoryIoErrorKind::TruncatedRecord =>
                {
                    let file_end = self
                        .reader
                        .seek_to(SeekFrom::End(0), TrajectoryIoOperation::Index)?;
                    return Ok(file_end.saturating_sub(offset));
                }
                Err(error) => return Err(error),
            }
            let end = self.reader.position(TrajectoryIoOperation::Index)?;
            Arc::make_mut(&mut self.offsets).push(
                offset,
                end,
                FrameStamp::default(),
                &self.reader.options.limits,
                TrajectoryFormat::Xyz,
                &label,
            )?;
            control.checkpoint(self.offsets.len(), end, TrajectoryFormat::Xyz, &label)?;
            self.reader.frame_cursor =
                self.reader.frame_cursor.checked_add(1).ok_or_else(|| {
                    frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::LimitExceeded,
                        TrajectoryIoOperation::Index,
                        &label,
                        self.reader.frame_cursor,
                        "The frame number exceeds this platform's address range.",
                    )
                })?;
        }
    }
}

impl<R: BufRead + Seek> TrajectoryReader for IndexedXyzReader<R> {
    fn topology(&self) -> &Topology {
        self.reader.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.reader.shared_topology()
    }

    /// Reads the indexed frames in order and ends after the last of them,
    /// even when the stream holds more.
    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        if self.topology().shares_layout(destination.topology())
            && self.reader.frame_cursor >= self.offsets.len()
        {
            return Ok(false);
        }
        self.reader.read_next(destination)
    }
}

impl<R: BufRead + Seek> SeekableTrajectoryReader for IndexedXyzReader<R> {
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
        self.read_indexed(index)?;
        self.reader.publish(destination)
    }
}

impl<R: BufRead + Seek> IndexedXyzReader<R> {
    /// Parses and verifies frame `index` into the reader's buffer, then
    /// puts the stream back where it was.
    pub(super) fn read_indexed(&mut self, index: usize) -> Result<(), TrajectoryIoError> {
        let (offset, _) = self
            .offsets
            .frame(index)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        let saved_offset = self.reader.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(FORMAT),
                &self.reader.options.source_label,
                error,
            )
        })?;
        let saved_cursor = self.reader.frame_cursor;
        self.reader
            .reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadFrame,
                    Some(FORMAT),
                    &self.reader.options.source_label,
                    error,
                )
            })?;
        self.reader.frame_cursor = index;
        let parsed = self.reader.parse_next(true).and_then(|parsed| {
            if parsed {
                self.offsets.verify_stream_frame_end(
                    index,
                    &mut self.reader.reader,
                    TrajectoryFormat::Xyz,
                    &self.reader.options.source_label,
                )?;
            }
            Ok(parsed)
        });
        let restore = self
            .reader
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadFrame,
                    Some(FORMAT),
                    &self.reader.options.source_label,
                    error,
                )
            });
        self.reader.frame_cursor = saved_cursor;
        if !parsed? {
            return Err(TrajectoryIoError::frame_index_out_of_range(index));
        }
        restore.map(|_| ())
    }
}
