use super::*;

/// Fully decoded-and-verified indexed XTC reader.
///
/// Random and sequential reads share decode scratch. A random read restores
/// stream/header state and invalidates cached coordinates, so the next sequential
/// frame is decoded again when it was already pending.
pub struct IndexedXtcReader<R> {
    pub(crate) inner: IndexedReader<XtcReader<R>>,
}

indexed_reader_traits!(IndexedXtcReader, Read + Seek);

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
        options: impl Into<TrajectoryOpenOptions>,
    ) -> Result<Self, TrajectoryIoError> {
        let options = options.into();
        let source = options.read_source(FORMAT);
        Self::build(
            reader,
            topology.into(),
            options.xtc,
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
        options: XtcReadOptions,
        source: ReadSource,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        let reader = XtcReader::open(reader, topology, options, source)?;
        let stream_start = reader.stream_start;
        let inner = IndexedReader::build(reader, stream_start, growing, control)?;
        Ok(Self { inner })
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; random reads verify each frame against them.
    pub(crate) fn from_offsets(
        reader: R,
        topology: Arc<Topology>,
        options: XtcReadOptions,
        source: ReadSource,
        offsets: Arc<FrameOffsets>,
    ) -> Result<Self, TrajectoryIoError> {
        let inner = XtcReader::open(reader, topology, options, source)?;
        offsets.verify_start(
            inner.stream_start,
            TrajectoryFormat::Xtc,
            &inner.source.label,
        )?;
        Ok(Self {
            inner: IndexedReader::new(inner, offsets),
        })
    }

    /// Another reader of the same frames over `stream`, which must read the
    /// same file from its start; it shares the index and reads the first
    /// frame alone.
    pub(crate) fn fork(&self, stream: R) -> Result<Self, TrajectoryIoError> {
        let reader = &self.inner.codec;
        Self::from_offsets(
            stream,
            Arc::clone(&reader.topology),
            reader.options.clone(),
            reader.source.clone(),
            Arc::clone(self.inner.shared_offsets()),
        )
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.inner.codec.topology()
    }

    pub(crate) fn first_info(&self) -> &XtcFrameInfo {
        self.inner.codec.first_info()
    }
}

/// A random read checks the frame's stamp against the index, and moving the
/// stream drops the coordinates decoded for a frame read ahead.
impl<R: Read + Seek> IndexedCodec for XtcReader<R> {
    type Stream = R;
    type Frame = XtcFrameInfo;
    /// The frame read ahead.
    type State = Option<XtcFrameInfo>;
    const FORMAT: TrajectoryFormat = FORMAT;
    const SEEK_OPERATION: TrajectoryIoOperation = TrajectoryIoOperation::Index;

    fn source(&self) -> &ReadSource {
        &self.source
    }

    fn stream(&self) -> &R {
        &self.adapter.reader
    }

    fn stream_and_source(&mut self) -> (&mut R, &ReadSource) {
        (&mut self.adapter.reader, &self.source)
    }

    fn frame_cursor(&self) -> usize {
        self.frame_cursor
    }

    fn set_frame_cursor(&mut self, frame: usize) {
        self.frame_cursor = frame;
    }

    fn save_state(&self) -> Option<XtcFrameInfo> {
        self.pending_info.clone()
    }

    fn restore_state(&mut self, pending: Option<XtcFrameInfo>) {
        self.pending_info = pending;
        self.adapter.decoded_start = None;
    }

    fn frame_pending(&self) -> bool {
        self.pending_info.is_some()
    }

    fn start_scan(&mut self, end: u64, _growing: bool) -> Result<(), TrajectoryIoError> {
        // A freshly opened reader has already checked frame 0 there.
        if self
            .pending_info
            .as_ref()
            .is_none_or(|info| info.start != end)
        {
            self.pending_info = None;
            self.adapter.decoded_start = None;
            self.seek_to(SeekFrom::Start(end), TrajectoryIoOperation::Index)?;
        }
        Ok(())
    }

    fn frame_offset(&mut self) -> Result<u64, TrajectoryIoError> {
        match &self.pending_info {
            Some(info) => Ok(info.start),
            None => self.position(TrajectoryIoOperation::Index),
        }
    }

    fn index_next(
        &mut self,
        _indexed: usize,
        offset: u64,
        _growing: bool,
    ) -> Result<Option<FrameStamp>, TrajectoryIoError> {
        match self.parse_next() {
            Ok(info) => Ok(info.map(|info| info.stamp())),
            Err(error) => Err(frame_offset_context(error, self.frame_cursor, offset)),
        }
    }

    fn start_read(&mut self, _offset: u64) {
        self.pending_info = None;
        self.adapter.decoded_start = None;
    }

    fn read_frame(&mut self) -> Result<Option<XtcFrameInfo>, TrajectoryIoError> {
        self.parse_next()
    }

    fn check_frame(
        &self,
        frame: &XtcFrameInfo,
        index: usize,
        offsets: &FrameOffsets,
    ) -> Result<(), TrajectoryIoError> {
        offsets.verify_stamp(
            index,
            frame.stamp(),
            TrajectoryIoOperation::ReadFrame,
            TrajectoryFormat::Xtc,
            &self.source.label,
        )
    }

    fn publish_frame(
        &self,
        frame: XtcFrameInfo,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        self.publish(&self.positions, &frame, destination)
    }

    fn take_decoded(
        &mut self,
        read: Result<XtcFrameInfo, TrajectoryIoError>,
    ) -> Result<DecodedGeometry, TrajectoryIoError> {
        let positions = std::mem::take(&mut self.positions);
        let info = read?;
        Ok(DecodedGeometry {
            positions: Quantity::new(positions, NANOMETER),
            cell: info.cell,
            velocities: None,
            forces: None,
            owner: None,
        })
    }
}
