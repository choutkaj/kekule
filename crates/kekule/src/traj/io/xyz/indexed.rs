use super::*;
use crate::traj::trajectory::reason::{Part, Reason};

/// Fully scanned XYZ reader with verified frame offsets.
pub struct IndexedXyzReader<R> {
    pub(crate) inner: IndexedReader<XyzReader<R>>,
}

indexed_reader_traits!(IndexedXyzReader, BufRead + Seek);

impl<R: BufRead + Seek> IndexedXyzReader<R> {
    /// Opens a stream as [`XyzReader::new`] does and indexes it for random
    /// access from its current position, frame 0, verifying every frame (see
    /// [`traj::io`](crate::traj::io)).
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
        options: impl Into<TrajectoryOpenOptions>,
    ) -> Result<Self, TrajectoryIoError> {
        let options = options.into();
        let source = options.read_source(FORMAT);
        Self::build(
            reader,
            topology.into(),
            options.xyz,
            source,
            false,
            &IndexBuildControl::default(),
        )
    }

    /// Indexes the stream like [`Self::new`], reporting to `control`. A
    /// `growing` stream may end in a partial frame, need not hold a complete
    /// one yet, and only lines with a line break count.
    pub(crate) fn build(
        mut reader: R,
        topology: Arc<Topology>,
        options: XyzReadOptions,
        source: ReadSource,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        let stream_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &source.label,
                error,
            )
        })?;
        let reader = XyzReader::open(reader, topology, options, source)?;
        let inner = IndexedReader::build(reader, stream_start, growing, control)?;
        if inner.offsets().is_empty() && !growing {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidHeader,
                TrajectoryIoOperation::Index,
                Some(FORMAT),
                &inner.codec.source.label,
                Reason::FileHoldsNoFrame,
            ));
        }
        Ok(Self { inner })
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; frame 0 is verified here, and random reads verify each
    /// frame against the offsets. A `growing` stream is extended as in
    /// [`Self::build`].
    pub(crate) fn from_offsets(
        mut reader: R,
        topology: Arc<Topology>,
        options: XyzReadOptions,
        source: ReadSource,
        offsets: Arc<FrameOffsets>,
        growing: bool,
    ) -> Result<Self, TrajectoryIoError> {
        let stream_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &source.label,
                error,
            )
        })?;
        let mut reader = XyzReader::open(reader, topology, options, source)?;
        reader.growing = growing;
        offsets.verify_start(stream_start, TrajectoryFormat::Xyz, &reader.source.label)?;
        reader.validate_first_frame()?;
        Ok(Self {
            inner: IndexedReader::new(reader, offsets),
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
            reader.growing,
        )
    }
}

/// Random reads seek with [`TrajectoryIoOperation::ReadFrame`] and report
/// errors without frame context, which XYZ frames carry themselves.
impl<R: BufRead + Seek> IndexedCodec for XyzReader<R> {
    type Stream = R;
    type Frame = ();
    type State = ();
    const FORMAT: TrajectoryFormat = FORMAT;
    const SEEK_OPERATION: TrajectoryIoOperation = TrajectoryIoOperation::ReadFrame;

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
            extra_column_frames: self.index_extra_columns,
            ..IndexFacts::default()
        }
    }

    fn set_index_facts(&mut self, facts: IndexFacts) {
        self.index_extra_columns = facts.extra_column_frames;
    }

    fn start_scan(&mut self, end: u64, growing: bool) -> Result<(), TrajectoryIoError> {
        self.growing = growing;
        self.seek_to(SeekFrom::Start(end), TrajectoryIoOperation::Index)
            .map(drop)
    }

    fn index_full(
        &mut self,
        _indexed: usize,
        _growing: bool,
        limit: LimitExceeded,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        let offset = self.frame_offset()?;
        let scan = ScanPoll { control, offset };
        if self.at_end(TrajectoryIoOperation::Index, Some(scan))? {
            return Ok(0);
        }
        let label = &self.source.label;
        Err(codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::Index,
            Some(FORMAT),
            label,
            None,
        )
        .exceeding(Part::FrameIndex, limit))
    }

    fn index_next(
        &mut self,
        _indexed: usize,
        offset: u64,
        _growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Option<FrameStamp>, TrajectoryIoError> {
        if !self.parse_next(false, Some(ScanPoll { control, offset }))? {
            return Ok(None);
        }
        if self.extra_columns {
            self.index_extra_columns = self.index_extra_columns.saturating_add(1);
        }
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            frame_error(
                FORMAT,
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Index,
                &self.source.label,
                self.frame_cursor,
                Reason::AddressRange(Part::FrameNumber),
            )
        })?;
        Ok(Some(FrameStamp::default()))
    }

    fn read_frame(&mut self) -> Result<Option<()>, TrajectoryIoError> {
        self.parse_next(true, None)
            .map(|parsed| parsed.then_some(()))
    }

    fn frame_context(error: TrajectoryIoError, _frame: usize, _offset: u64) -> TrajectoryIoError {
        error
    }

    fn publish_frame(
        &self,
        (): (),
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        self.publish(destination)
    }

    fn take_decoded(
        &mut self,
        read: Result<(), TrajectoryIoError>,
    ) -> Result<DecodedGeometry, TrajectoryIoError> {
        let positions = std::mem::take(&mut self.positions);
        read?;
        Ok(DecodedGeometry {
            positions: Quantity::new(positions, self.options.length_unit),
            cell: None,
            velocities: None,
            forces: None,
            owner: None,
        })
    }
}
