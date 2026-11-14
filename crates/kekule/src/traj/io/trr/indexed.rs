use super::*;

/// Fully verified indexed TRR reader.
///
/// Sequential and random reads share one set of decode scratch arrays. Random
/// reads restore stream position and header state before publishing a frame.
pub struct IndexedTrrReader<R> {
    pub(super) inner: TrrReader<R>,
    pub(super) offsets: Arc<FrameOffsets>,
    /// Whether the indexed frames mix scalar widths.
    pub(super) precision_mixed: bool,
    /// Bytes after the last indexed frame when the index last grew.
    pub(super) pending_bytes: u64,
}

impl<R> fmt::Debug for IndexedTrrReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedTrrReader")
            .field("reader", &self.inner)
            .field("frame_count", &self.offsets.len())
            .finish_non_exhaustive()
    }
}

impl<R: Read + Seek> IndexedTrrReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order
    /// and indexes it for random access.
    ///
    /// The stream's position at construction is frame 0. Every frame from there to
    /// the end is decoded and verified before the reader is returned, and both
    /// sequential and indexed reads then address frames from that start. Topology,
    /// limits, and diagnostic labels follow [`TrrReader::new`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`TrrReader::new`], the decoding error of the
    /// first invalid frame, [`TrajectoryIoErrorKind::LimitExceeded`] when the
    /// index outgrows the limits, and [`TrajectoryIoErrorKind::Io`] when
    /// reading or seeking fails.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: TrrReadOptions,
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
        options: TrrReadOptions,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        TrrReader::new(reader, topology, options)?.index_frames(growing, control)
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; random reads verify each frame against them, including
    /// that its scalar width is the first frame's unless `precision_mixed`.
    pub(crate) fn from_offsets(
        reader: R,
        topology: Arc<Topology>,
        options: TrrReadOptions,
        offsets: Arc<FrameOffsets>,
        precision_mixed: bool,
    ) -> Result<Self, TrajectoryIoError> {
        let mut inner = TrrReader::new(reader, topology, options)?;
        offsets.verify_start(
            inner.stream_start,
            TrajectoryFormat::Trr,
            &inner.options.source_label,
        )?;
        inner.precision_mixed = precision_mixed;
        Ok(Self {
            inner,
            offsets,
            precision_mixed,
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
        self.precision_mixed = other.precision_mixed;
    }

    /// Another reader of the same frames over `stream`, which must read the
    /// same file from its start; it shares the index and reads the first
    /// frame's header alone.
    pub(crate) fn fork(&self, stream: R) -> Result<Self, TrajectoryIoError> {
        Self::from_offsets(
            stream,
            Arc::clone(&self.inner.topology),
            self.inner.options.clone(),
            Arc::clone(&self.offsets),
            self.precision_mixed,
        )
    }

    /// The stream every read goes through.
    pub(crate) fn stream(&self) -> &R {
        &self.inner.reader
    }

    /// The stream every read goes through.
    pub(crate) fn stream_mut(&mut self) -> &mut R {
        &mut self.inner.reader
    }

    /// Reads frame `index` as [`SeekableTrajectoryReader::read_frame`]
    /// does and hands over its positions, velocities, and forces, each
    /// written once as it was decoded. The read allocates them on the
    /// calling thread and the reader keeps no buffer for its next read,
    /// which may run on another.
    pub(crate) fn read_decoded(
        &mut self,
        index: usize,
    ) -> Result<DecodedGeometry, TrajectoryIoError> {
        let read = self.read_indexed(index);
        // Taken whether the read succeeds or fails, and whether the frame
        // holds velocities and forces or not, so that the reader keeps no
        // array for a read on another thread.
        let positions = std::mem::take(&mut self.inner.positions);
        let velocities = std::mem::take(&mut self.inner.velocities);
        let forces = std::mem::take(&mut self.inner.forces);
        let decoded = read?;
        Ok(DecodedGeometry {
            positions: Quantity::new(positions, NANOMETER),
            cell: decoded.cell,
            velocities: decoded
                .header
                .has_velocities()
                .then(|| Quantity::new(velocities, NANOMETER_PER_PICOSECOND)),
            forces: decoded
                .header
                .has_forces()
                .then(|| Quantity::new(forces, KILOJOULE_PER_MOLE_PER_NANOMETER)),
            owner: Some(self.inner.properties.clone()),
        })
    }

    /// Bytes after the last indexed frame when the index last grew.
    pub(crate) const fn pending_bytes(&self) -> u64 {
        self.pending_bytes
    }

    /// Drops the frames indexed after `mark`, with the mixed-precision flag
    /// and pending bytes the index had then; see [`FrameOffsets::mark`].
    pub(crate) fn restore_index(
        &mut self,
        mark: (usize, u64),
        precision_mixed: bool,
        pending_bytes: u64,
    ) {
        Arc::make_mut(&mut self.offsets).restore(mark);
        self.precision_mixed = precision_mixed;
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
        let saved_pending = self.inner.pending_header.clone();
        let saved_header_offset = self.inner.current_header_offset;
        let saved_precision_mixed = self.inner.precision_mixed;
        let end = self.offsets.end();
        // A freshly opened reader has already read frame 0's header there.
        if saved_pending.is_none() || saved_header_offset != end {
            self.inner.pending_header = None;
            self.inner
                .seek_to(SeekFrom::Start(end), TrajectoryIoOperation::Index)?;
        }
        self.inner.current_header_offset = end;
        self.inner.frame_cursor = self.offsets.len();
        self.inner.precision_mixed = self.precision_mixed;
        let scanned = self.scan(growing, control);
        let scanned_mixed = self.inner.precision_mixed;
        let restored = self
            .inner
            .seek_to(SeekFrom::Start(saved_offset), TrajectoryIoOperation::Index);
        self.inner.frame_cursor = saved_cursor;
        self.inner.pending_header = saved_pending;
        self.inner.current_header_offset = saved_header_offset;
        self.inner.precision_mixed = saved_precision_mixed;
        match scanned.and_then(|pending| restored.map(|_| pending)) {
            Ok(pending) => {
                self.precision_mixed = scanned_mixed;
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
            TrajectoryFormat::Trr,
            &label,
        )?;
        loop {
            if let Some(limit) =
                projected_index_limit(self.offsets.len(), &self.inner.options.limits)
            {
                let clean_eof = self.inner.pending_header.is_none()
                    && probe_seekable_eof(
                        &mut self.inner.reader,
                        TrajectoryIoOperation::Index,
                        TrajectoryFormat::Trr,
                        &label,
                    )?;
                if clean_eof {
                    return Ok(0);
                }
                return Err(resource_error(&label, Some(self.offsets.len()), "")
                    .exceeding("The frame index", limit));
            }
            let offset = if self.inner.pending_header.is_some() {
                self.inner.current_header_offset
            } else {
                self.inner.position(TrajectoryIoOperation::Index)?
            };
            let decoded = match self.inner.parse_next() {
                Ok(Some(decoded)) => decoded,
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
                decoded.header.stamp(),
                &self.inner.options.limits,
                TrajectoryFormat::Trr,
                &label,
            )?;
            control.checkpoint(self.offsets.len(), end, TrajectoryFormat::Trr, &label)?;
        }
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.inner.topology()
    }

    pub(crate) fn first_header(&self) -> &TrrFrameHeader {
        self.inner.first_header()
    }

    pub(crate) const fn precision_mixed(&self) -> bool {
        self.precision_mixed
    }
}

impl<R: Read + Seek> TrajectoryReader for IndexedTrrReader<R> {
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

impl<R: Read + Seek> SeekableTrajectoryReader for IndexedTrrReader<R> {
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
        let decoded = self.read_indexed(index)?;
        self.inner.publish(
            &self.inner.positions,
            &self.inner.velocities,
            &self.inner.forces,
            &self.inner.properties,
            &decoded,
            destination,
        )
    }
}

impl<R: Read + Seek> IndexedTrrReader<R> {
    /// Decodes and verifies frame `index` into the reader's buffers, then
    /// puts the stream back where it was.
    pub(super) fn read_indexed(
        &mut self,
        index: usize,
    ) -> Result<TrrDecodedFrame, TrajectoryIoError> {
        let (offset, _) = self
            .offsets
            .frame(index)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        let saved_offset = self.inner.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Trr),
                &self.inner.options.source_label,
                error,
            )
        })?;
        let saved_cursor = self.inner.frame_cursor;
        let saved_pending = self.inner.pending_header.clone();
        let saved_header_offset = self.inner.current_header_offset;
        let saved_precision_mixed = self.inner.precision_mixed;
        self.inner
            .reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Trr),
                    &self.inner.options.source_label,
                    error,
                )
            })?;
        self.inner.pending_header = None;
        self.inner.current_header_offset = offset;
        self.inner.frame_cursor = index;
        let result = self
            .inner
            .parse_next()
            .map_err(|error| frame_offset_context(error, index, offset))
            .and_then(|decoded| decoded.ok_or(TrajectoryIoError::frame_index_out_of_range(index)))
            .and_then(|decoded| {
                self.offsets.verify_stream_frame_end(
                    index,
                    &mut self.inner.reader,
                    TrajectoryFormat::Trr,
                    &self.inner.options.source_label,
                )?;
                self.offsets.verify_stamp(
                    index,
                    decoded.header.stamp(),
                    TrajectoryIoOperation::ReadFrame,
                    TrajectoryFormat::Trr,
                    &self.inner.options.source_label,
                )?;
                if !self.precision_mixed
                    && decoded.header.precision != self.inner.first_header.precision
                {
                    return Err(index_mismatch(
                        TrajectoryIoOperation::ReadFrame,
                        TrajectoryFormat::Trr,
                        &self.inner.options.source_label,
                        "The frame's scalar width differs from the first frame's.",
                    )
                    .with_frame(index));
                }
                Ok(decoded)
            });
        let restore = self
            .inner
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Trr),
                    &self.inner.options.source_label,
                    error,
                )
            });
        self.inner.frame_cursor = saved_cursor;
        self.inner.pending_header = saved_pending;
        self.inner.current_header_offset = saved_header_offset;
        self.inner.precision_mixed = saved_precision_mixed;
        let decoded = result?;
        restore?;
        Ok(decoded)
    }
}
