use super::*;

/// Fully verified DCD indexed reader. Random reads leave the fixed-atom
/// reference coordinates as they are.
pub struct IndexedDcdReader<R> {
    pub(crate) inner: IndexedReader<DcdReader<R>>,
}

indexed_reader_traits!(IndexedDcdReader, Read + Seek);

impl<R: Read + Seek> IndexedDcdReader<R> {
    /// Opens a stream as [`DcdReader::new`] does and indexes it for random
    /// access from its current position, frame 0, verifying every frame (see
    /// [`traj::io`](crate::traj::io)).
    ///
    /// # Errors
    ///
    /// Returns the errors of [`DcdReader::new`], the decoding error of the
    /// first invalid frame, [`TrajectoryIoErrorKind::LimitExceeded`] when
    /// the index outgrows the limits, and [`TrajectoryIoErrorKind::Io`] when
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
            options.dcd,
            source,
            false,
            &IndexBuildControl::default(),
        )
    }

    /// Indexes the stream like [`Self::new`], reporting to `control`. A
    /// `growing` stream may end in a partial frame.
    pub(crate) fn build(
        reader: R,
        topology: Arc<Topology>,
        options: DcdReadOptions,
        source: ReadSource,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        let reader = DcdReader::open(reader, topology, options, source)?;
        let data_start = reader.header.data_start;
        let inner = IndexedReader::build(reader, data_start, growing, control)?;
        Ok(Self { inner })
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; random reads verify each frame against them. With fixed
    /// atoms, frame 0 is decoded and verified here for their reference
    /// coordinates.
    pub(crate) fn from_offsets(
        reader: R,
        topology: Arc<Topology>,
        options: DcdReadOptions,
        source: ReadSource,
        offsets: Arc<FrameOffsets>,
    ) -> Result<Self, TrajectoryIoError> {
        let mut inner = DcdReader::open(reader, topology, options, source)?;
        let label = inner.source.label.clone();
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
            inner: IndexedReader::new(inner, offsets),
        })
    }

    /// Another reader of the same frames over `stream`, which must read the
    /// same file from its start; it shares the index and reads the header,
    /// and with fixed atoms frame 0.
    pub(crate) fn fork(&self, stream: R) -> Result<Self, TrajectoryIoError> {
        let reader = &self.inner.codec;
        let mut fork = Self::from_offsets(
            stream,
            Arc::clone(&reader.topology),
            reader.options.clone(),
            reader.source.clone(),
            Arc::clone(self.inner.shared_offsets()),
        )?;
        fork.inner.codec.header.declared_frames = reader.header.declared_frames;
        Ok(fork)
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.inner.codec.topology()
    }

    pub(crate) fn header(&self) -> &DcdHeader {
        self.inner.codec.header()
    }
}

/// A scan reads the header's `NSET` again first, since a writer appending
/// frames updates it, and records it in the index facts. The frames are
/// counted in the stream, whatever `NSET` says.
impl<R: Read + Seek> IndexedCodec for DcdReader<R> {
    type Stream = R;
    type Frame = DcdDecodedFrame;
    type State = ();
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

    fn save_state(&self) {}

    fn restore_state(&mut self, (): ()) {}

    fn index_facts(&self) -> IndexFacts {
        IndexFacts {
            declared_frames: Some(self.header.declared_frames),
            ..IndexFacts::default()
        }
    }

    fn set_index_facts(&mut self, facts: IndexFacts) {
        if let Some(declared) = facts.declared_frames {
            self.header.declared_frames = declared;
        }
    }

    fn start_scan(&mut self, end: u64, _growing: bool) -> Result<(), TrajectoryIoError> {
        self.header.declared_frames = self.reread_declared_frames()?;
        self.seek_to(SeekFrom::Start(end), TrajectoryIoOperation::Index)
            .map(drop)
    }

    fn index_next(
        &mut self,
        _indexed: usize,
        offset: u64,
        _growing: bool,
        _control: &IndexBuildControl<'_>,
    ) -> Result<Option<FrameStamp>, TrajectoryIoError> {
        let frame = self.frame_cursor;
        let stamp = self
            .stamp(frame)
            .map_err(|error| frame_offset_context(error, frame, offset))?;
        match self.parse_next(true) {
            Ok(decoded) => Ok(decoded.map(|_| stamp)),
            Err(error) => Err(frame_offset_context(error, frame, offset)),
        }
    }

    fn read_frame(&mut self) -> Result<Option<DcdDecodedFrame>, TrajectoryIoError> {
        self.parse_next(false)
    }

    fn publish_frame(
        &self,
        frame: DcdDecodedFrame,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        self.publish(&self.positions, frame, destination)
    }

    /// Without fixed atoms the positions are written once as they are
    /// decoded; with them, a frame after the first fills a fresh buffer with
    /// the fixed atoms first.
    fn take_decoded(
        &mut self,
        read: Result<DcdDecodedFrame, TrajectoryIoError>,
    ) -> Result<DecodedGeometry, TrajectoryIoError> {
        let positions = std::mem::take(&mut self.positions);
        let decoded = read?;
        Ok(DecodedGeometry {
            positions: Quantity::new(positions, ANGSTROM),
            cell: decoded.cell,
            velocities: None,
            forces: None,
            owner: None,
        })
    }
}
