use super::*;

/// Fully verified DCD indexed reader.
///
/// Random and sequential reads share one set of decode scratch. Fixed-atom
/// reference coordinates remain immutable during random reads.
pub struct IndexedDcdReader<R> {
    pub(super) inner: DcdReader<R>,
    pub(super) offsets: Arc<FrameOffsets>,
    /// Bytes after the last indexed frame when the index last grew.
    pub(super) pending_bytes: u64,
}

impl<R> fmt::Debug for IndexedDcdReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedDcdReader")
            .field("reader", &self.inner)
            .field("frame_count", &self.offsets.len())
            .finish_non_exhaustive()
    }
}

impl<R: Read + Seek> IndexedDcdReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order
    /// and indexes it for random access.
    ///
    /// The stream's position at construction is frame 0. Every frame from there to
    /// the end is decoded and verified before the reader is returned, and both
    /// sequential and indexed reads then address frames from that start. Topology,
    /// limits, and diagnostic labels follow [`DcdReader::new`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`DcdReader::new`], the decoding error of the
    /// first invalid frame, [`TrajectoryIoErrorKind::InconsistentMetadata`]
    /// when the frame count differs from the header `NSET` count,
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when the index outgrows the
    /// limits, and [`TrajectoryIoErrorKind::Io`] when reading or seeking
    /// fails.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: DcdReadOptions,
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
    /// `growing` stream may end in a partial frame, and its frames are
    /// counted in the stream rather than checked against `NSET` (see
    /// [`Self::extend`]).
    pub(crate) fn build(
        reader: R,
        topology: Arc<Topology>,
        options: DcdReadOptions,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        DcdReader::new(reader, topology, options)?.index_frames(growing, control)
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; random reads verify each frame against them. Unless the
    /// stream is `growing`, the offsets must count at most the header's
    /// `NSET` frames, and [`Self::extend`] then requires the rest. With
    /// fixed atoms, frame 0 is decoded and verified here for their reference
    /// coordinates.
    pub(crate) fn from_offsets(
        reader: R,
        topology: Arc<Topology>,
        options: DcdReadOptions,
        offsets: Arc<FrameOffsets>,
        growing: bool,
    ) -> Result<Self, TrajectoryIoError> {
        let mut inner = DcdReader::new(reader, topology, options)?;
        inner.growing = growing;
        let label = inner.options.source_label.clone();
        if !growing && offsets.len() > inner.header.declared_frames {
            return Err(index_mismatch(
                TrajectoryIoOperation::Open,
                TrajectoryFormat::Dcd,
                &label,
                "The index counts more frames than the frame count NSET in the header.",
            )
            .with_counts(inner.header.declared_frames as u64, offsets.len() as u64));
        }
        offsets.verify_start(inner.header.data_start, TrajectoryFormat::Dcd, &label)?;
        for frame in 0..offsets.len() {
            offsets.verify_stamp(
                frame,
                inner.stamp(frame)?,
                TrajectoryIoOperation::Open,
                TrajectoryFormat::Dcd,
                &label,
            )?;
        }
        if inner.header.fixed_count > 0 && !offsets.is_empty() {
            let data_start = inner.header.data_start;
            inner
                .parse_next(true)
                .map_err(|error| frame_offset_context(error, 0, data_start))?
                .ok_or(TrajectoryIoError::frame_index_out_of_range(0))?;
            offsets.verify_stream_frame_end(0, &mut inner.reader, TrajectoryFormat::Dcd, &label)?;
            inner
                .reader
                .seek(SeekFrom::Start(data_start))
                .map_err(|error| {
                    io_context(TrajectoryIoOperation::Open, Some(FORMAT), &label, error)
                })?;
            inner.frame_cursor = 0;
        }
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

    /// Reads by the index and frame count `NSET` of `other`, a reader of
    /// the same file, from now on.
    pub(crate) fn share_index(&mut self, other: &Self) {
        self.offsets = Arc::clone(&other.offsets);
        self.inner.header.declared_frames = other.inner.header.declared_frames;
    }

    /// Another reader of the same frames over `stream`, which must read the
    /// same file from its start; it shares the index and reads the header,
    /// and with fixed atoms frame 0.
    pub(crate) fn fork(&self, stream: R) -> Result<Self, TrajectoryIoError> {
        let mut fork = Self::from_offsets(
            stream,
            Arc::clone(&self.inner.topology),
            self.inner.options.clone(),
            Arc::clone(&self.offsets),
            self.inner.growing,
        )?;
        fork.inner.header.declared_frames = self.inner.header.declared_frames;
        Ok(fork)
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
    /// does and hands over its positions. Without fixed atoms they are
    /// written once as they are decoded; with them, a frame after the first
    /// fills a fresh buffer with the fixed atoms first. The read allocates
    /// them on the calling thread and the reader keeps no buffer for its
    /// next read, which may run on another.
    pub(crate) fn read_decoded(
        &mut self,
        index: usize,
    ) -> Result<DecodedGeometry, TrajectoryIoError> {
        let read = self.read_indexed(index);
        // Taken whether the read succeeds or fails, so that the reader keeps
        // no array for a read on another thread.
        let positions = std::mem::take(&mut self.inner.positions);
        let decoded = read?;
        Ok(DecodedGeometry {
            positions: Quantity::new(positions, ANGSTROM),
            cell: decoded.cell,
            velocities: None,
            forces: None,
            owner: None,
        })
    }

    /// Bytes after the last indexed frame when the index last grew.
    pub(crate) const fn pending_bytes(&self) -> u64 {
        self.pending_bytes
    }

    /// Drops the frames indexed after `mark`, with the declared frame count
    /// and pending bytes the reader had then; see [`FrameOffsets::mark`].
    pub(crate) fn restore_index(
        &mut self,
        mark: (usize, u64),
        declared_frames: usize,
        pending_bytes: u64,
    ) {
        Arc::make_mut(&mut self.offsets).restore(mark);
        self.inner.header.declared_frames = declared_frames;
        self.pending_bytes = pending_bytes;
    }

    /// Verifies and indexes every complete frame after the last indexed one,
    /// up to the end of the stream, and returns the bytes left after them.
    ///
    /// The header's `NSET` is read again first, since a writer appending
    /// frames updates it. A growing stream (see [`Self::build`]) may end in
    /// a frame still being written: the scan stops before a frame the stream
    /// ends inside, which is not yet an error, and `NSET` is not checked.
    /// Otherwise such a frame is an error, as is any frame that is invalid
    /// with all its bytes present, and the stream must hold exactly `NSET`
    /// frames. Failure leaves the index as it was; the sequential cursor is
    /// kept either way.
    pub(crate) fn extend(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        self.inner.growing = growing;
        let mark = self.offsets.mark();
        let declared = self.inner.header.declared_frames;
        let saved_offset = self.inner.position(TrajectoryIoOperation::Index)?;
        let saved_cursor = self.inner.frame_cursor;
        let scanned = self
            .inner
            .reread_declared_frames()
            .and_then(|declared| {
                self.inner.header.declared_frames = declared;
                self.inner.seek_to(
                    SeekFrom::Start(self.offsets.end()),
                    TrajectoryIoOperation::Index,
                )
            })
            .and_then(|_| {
                self.inner.frame_cursor = self.offsets.len();
                self.scan(control)
            });
        let restored = self
            .inner
            .seek_to(SeekFrom::Start(saved_offset), TrajectoryIoOperation::Index);
        self.inner.frame_cursor = saved_cursor;
        match scanned.and_then(|pending| restored.map(|_| pending)) {
            Ok(pending) => {
                self.pending_bytes = pending;
                Ok(pending)
            }
            Err(error) => {
                Arc::make_mut(&mut self.offsets).restore(mark);
                self.inner.header.declared_frames = declared;
                Err(error)
            }
        }
    }

    /// Indexes frames from the current position; see [`Self::extend`].
    pub(super) fn scan(
        &mut self,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        let label = self.inner.options.source_label.clone();
        let growing = self.inner.growing;
        let declared = self.inner.header.declared_frames;
        control.checkpoint(
            self.offsets.len(),
            self.offsets.end(),
            TrajectoryFormat::Dcd,
            &label,
        )?;
        let ended_early = |indexed: usize| {
            TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Index,
                Some(FORMAT),
            )
            .with_source_label(&label)
            .with_counts(declared as u64, indexed as u64)
            .with_detail("The file ends before the frame count NSET in its header.")
        };
        loop {
            if !growing && self.offsets.len() >= declared {
                if probe_seekable_eof(
                    &mut self.inner.reader,
                    TrajectoryIoOperation::Index,
                    TrajectoryFormat::Dcd,
                    &label,
                )? {
                    return Ok(0);
                }
                return Err(TrajectoryIoError::codec(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Index,
                    Some(FORMAT),
                )
                .with_source_label(&label)
                .with_counts(
                    declared as u64,
                    (self.offsets.len() as u64).saturating_add(1),
                )
                .with_detail(
                    "The file holds more frames than the frame count NSET in its header.",
                ));
            }
            if let Some(limit) =
                projected_index_limit(self.offsets.len(), &self.inner.options.limits)
            {
                if probe_seekable_eof(
                    &mut self.inner.reader,
                    TrajectoryIoOperation::Index,
                    TrajectoryFormat::Dcd,
                    &label,
                )? {
                    if growing {
                        return Ok(0);
                    }
                    return Err(ended_early(self.offsets.len()));
                }
                return Err(limit_error(
                    FORMAT,
                    TrajectoryIoOperation::Index,
                    &label,
                    Some(self.offsets.len()),
                    "",
                )
                .exceeding("The frame index", limit));
            }
            let offset = self.inner.position(TrajectoryIoOperation::Index)?;
            let frame = self.inner.frame_cursor;
            let stamp = self
                .inner
                .stamp(frame)
                .map_err(|error| frame_offset_context(error, frame, offset))?;
            match self.inner.parse_next(true) {
                Ok(Some(_)) => {}
                Ok(None) if growing => return Ok(0),
                Ok(None) => return Err(ended_early(self.offsets.len())),
                Err(error) if growing && error.kind() == TrajectoryIoErrorKind::TruncatedRecord => {
                    let file_end = self
                        .inner
                        .seek_to(SeekFrom::End(0), TrajectoryIoOperation::Index)?;
                    return Ok(file_end.saturating_sub(offset));
                }
                Err(error) => return Err(frame_offset_context(error, frame, offset)),
            }
            let end = self.inner.position(TrajectoryIoOperation::Index)?;
            Arc::make_mut(&mut self.offsets).push(
                offset,
                end,
                stamp,
                &self.inner.options.limits,
                TrajectoryFormat::Dcd,
                &label,
            )?;
            control.checkpoint(self.offsets.len(), end, TrajectoryFormat::Dcd, &label)?;
        }
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.inner.topology()
    }

    pub(crate) fn header(&self) -> &DcdHeader {
        self.inner.header()
    }
}

impl<R: Read + Seek> TrajectoryReader for IndexedDcdReader<R> {
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

impl<R: Read + Seek> SeekableTrajectoryReader for IndexedDcdReader<R> {
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
        self.inner
            .publish(&self.inner.positions, decoded, destination)
    }
}

impl<R: Read + Seek> IndexedDcdReader<R> {
    /// Decodes and verifies frame `index` into the reader's buffers, then
    /// puts the stream back where it was.
    pub(super) fn read_indexed(
        &mut self,
        index: usize,
    ) -> Result<DcdDecodedFrame, TrajectoryIoError> {
        let (offset, _) = self
            .offsets
            .frame(index)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        let saved_offset = self.inner.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Index,
                Some(FORMAT),
                &self.inner.options.source_label,
                error,
            )
        })?;
        let saved_cursor = self.inner.frame_cursor;
        self.inner
            .reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(FORMAT),
                    &self.inner.options.source_label,
                    error,
                )
            })?;
        self.inner.frame_cursor = index;
        let result = self
            .inner
            .parse_next(false)
            .map_err(|error| frame_offset_context(error, index, offset))
            .and_then(|decoded| decoded.ok_or(TrajectoryIoError::frame_index_out_of_range(index)))
            .and_then(|decoded| {
                self.offsets.verify_stream_frame_end(
                    index,
                    &mut self.inner.reader,
                    TrajectoryFormat::Dcd,
                    &self.inner.options.source_label,
                )?;
                Ok(decoded)
            });
        let restore = self
            .inner
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(FORMAT),
                    &self.inner.options.source_label,
                    error,
                )
            });
        self.inner.frame_cursor = saved_cursor;
        let decoded = result?;
        restore?;
        Ok(decoded)
    }
}
