use super::*;

/// Fully verified indexed TRR reader.
///
/// Sequential and random reads share one set of decode scratch arrays. Random
/// reads restore stream position and header state before publishing a frame.
pub struct IndexedTrrReader<R> {
    pub(crate) inner: IndexedReader<TrrReader<R>>,
}

indexed_reader_traits!(IndexedTrrReader, Read + Seek);

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
        options: impl Into<TrajectoryOpenOptions>,
    ) -> Result<Self, TrajectoryIoError> {
        let options = options.into();
        let source = options.read_source(FORMAT);
        Self::build(
            reader,
            topology.into(),
            options.trr,
            source,
            false,
            &IndexBuildControl::default(),
        )
    }

    /// Indexes the stream like [`Self::new`], reporting to `control`; a
    /// `growing` stream may end in a partial frame.
    pub(crate) fn build(
        reader: R,
        topology: Arc<Topology>,
        options: TrrReadOptions,
        source: ReadSource,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        let reader = TrrReader::open(reader, topology, options, source)?;
        let stream_start = reader.stream_start;
        let inner = IndexedReader::build(reader, stream_start, growing, control)?;
        Ok(Self { inner })
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; random reads verify each frame against them, including
    /// that its scalar width is the first frame's unless `precision_mixed`.
    pub(crate) fn from_offsets(
        reader: R,
        topology: Arc<Topology>,
        options: TrrReadOptions,
        source: ReadSource,
        offsets: Arc<FrameOffsets>,
        precision_mixed: bool,
    ) -> Result<Self, TrajectoryIoError> {
        let mut inner = TrrReader::open(reader, topology, options, source)?;
        offsets.verify_start(
            inner.stream_start,
            TrajectoryFormat::Trr,
            &inner.source.label,
        )?;
        inner.precision_mixed = precision_mixed;
        inner.index_mixed = precision_mixed;
        Ok(Self {
            inner: IndexedReader::new(inner, offsets),
        })
    }

    /// Another reader of the same frames over `stream`, which must read the
    /// same file from its start; it shares the index and reads the first
    /// frame's header alone.
    pub(crate) fn fork(&self, stream: R) -> Result<Self, TrajectoryIoError> {
        let reader = &self.inner.codec;
        Self::from_offsets(
            stream,
            Arc::clone(&reader.topology),
            reader.options.clone(),
            reader.source.clone(),
            Arc::clone(self.inner.shared_offsets()),
            reader.index_mixed,
        )
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.inner.codec.topology()
    }

    pub(crate) fn first_header(&self) -> &TrrFrameHeader {
        self.inner.codec.first_header()
    }
}

/// The index facts record whether the indexed frames mix scalar widths,
/// apart from what sequential reads have seen. A random read checks the
/// frame's stamp, and its scalar width against the first frame's unless the
/// index is mixed.
impl<R: Read + Seek> IndexedCodec for TrrReader<R> {
    type Stream = R;
    type Frame = TrrDecodedFrame;
    /// The header read ahead, where it starts, and whether sequential reads
    /// have seen mixed widths.
    type State = (Option<TrrFrameHeader>, u64, bool);
    const FORMAT: TrajectoryFormat = FORMAT;
    const SEEK_OPERATION: TrajectoryIoOperation = TrajectoryIoOperation::Index;

    fn source(&self) -> &ReadSource {
        &self.source
    }

    fn stream(&self) -> &R {
        &self.reader
    }

    fn stream_and_source(&mut self) -> (&mut R, &ReadSource) {
        (&mut self.reader, &self.source)
    }

    fn frame_cursor(&self) -> usize {
        self.frame_cursor
    }

    fn set_frame_cursor(&mut self, frame: usize) {
        self.frame_cursor = frame;
    }

    fn save_state(&self) -> Self::State {
        (
            self.pending_header.clone(),
            self.current_header_offset,
            self.precision_mixed,
        )
    }

    fn restore_state(&mut self, (pending, header_offset, precision_mixed): Self::State) {
        self.pending_header = pending;
        self.current_header_offset = header_offset;
        self.precision_mixed = precision_mixed;
    }

    fn index_facts(&self) -> IndexFacts {
        IndexFacts {
            declared_frames: None,
            precision_mixed: self.index_mixed,
        }
    }

    fn set_index_facts(&mut self, facts: IndexFacts) {
        self.index_mixed = facts.precision_mixed;
    }

    fn frame_pending(&self) -> bool {
        self.pending_header.is_some()
    }

    fn start_scan(&mut self, end: u64, _growing: bool) -> Result<(), TrajectoryIoError> {
        // A freshly opened reader has already read frame 0's header there.
        if self.pending_header.is_none() || self.current_header_offset != end {
            self.pending_header = None;
            self.seek_to(SeekFrom::Start(end), TrajectoryIoOperation::Index)?;
        }
        self.current_header_offset = end;
        self.precision_mixed = self.index_mixed;
        Ok(())
    }

    fn end_scan(&mut self) {
        self.index_mixed = self.precision_mixed;
    }

    fn frame_offset(&mut self) -> Result<u64, TrajectoryIoError> {
        if self.pending_header.is_some() {
            Ok(self.current_header_offset)
        } else {
            self.position(TrajectoryIoOperation::Index)
        }
    }

    fn index_next(
        &mut self,
        _indexed: usize,
        offset: u64,
        _growing: bool,
    ) -> Result<Option<FrameStamp>, TrajectoryIoError> {
        match self.parse_next() {
            Ok(decoded) => Ok(decoded.map(|decoded| decoded.header.stamp())),
            Err(error) => Err(frame_offset_context(error, self.frame_cursor, offset)),
        }
    }

    fn start_read(&mut self, offset: u64) {
        self.pending_header = None;
        self.current_header_offset = offset;
    }

    fn read_frame(&mut self) -> Result<Option<TrrDecodedFrame>, TrajectoryIoError> {
        self.parse_next()
    }

    fn check_frame(
        &self,
        frame: &TrrDecodedFrame,
        index: usize,
        offsets: &FrameOffsets,
    ) -> Result<(), TrajectoryIoError> {
        offsets.verify_stamp(
            index,
            frame.header.stamp(),
            TrajectoryIoOperation::ReadFrame,
            TrajectoryFormat::Trr,
            &self.source.label,
        )?;
        if !self.index_mixed && frame.header.precision != self.first_header.precision {
            return Err(index_mismatch(
                TrajectoryIoOperation::ReadFrame,
                TrajectoryFormat::Trr,
                &self.source.label,
                "The frame's scalar width differs from the first frame's.",
            )
            .with_frame(index));
        }
        Ok(())
    }

    fn publish_frame(
        &self,
        frame: TrrDecodedFrame,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        self.publish(
            &self.positions,
            &self.velocities,
            &self.forces,
            &self.properties,
            &frame,
            destination,
        )
    }

    /// Velocities and forces are taken whether the frame holds them or not.
    fn take_decoded(
        &mut self,
        read: Result<TrrDecodedFrame, TrajectoryIoError>,
    ) -> Result<DecodedGeometry, TrajectoryIoError> {
        let positions = std::mem::take(&mut self.positions);
        let velocities = std::mem::take(&mut self.velocities);
        let forces = std::mem::take(&mut self.forces);
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
            owner: Some(self.properties.clone()),
        })
    }
}
