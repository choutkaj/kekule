use super::*;
use crate::traj::trajectory::reason::{Compressed, Part, Reason};

pub(super) struct CheckedXtcReaderAdapter<R> {
    pub(super) reader: R,
    pub(super) scratch: Vec<u8>,
    pub(super) decoded_positions: Vec<f32>,
    pub(super) decoded_start: Option<u64>,
}

impl<R: Read + Seek> CheckedXtcReaderAdapter<R> {
    pub(super) fn new(reader: R) -> Self {
        Self {
            reader,
            scratch: Vec::new(),
            decoded_positions: Vec::new(),
            decoded_start: None,
        }
    }

    pub(super) fn preflight(
        &mut self,
        topology: &Topology,
        options: &XtcReadOptions,
        source: &ReadSource,
        frame_index: usize,
        clean_eof: bool,
        positions_capacity: usize,
    ) -> Result<Option<XtcFrameInfo>, TrajectoryIoError> {
        let limits = &source.limits;
        let source_label = &source.label;
        self.decoded_start = None;
        let start = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                error,
            )
        })?;
        let mut header = [0_u8; XTC_HEADER_BYTES];
        match self.reader.read(&mut header[..1]) {
            Ok(0) if clean_eof => return Ok(None),
            Ok(0) => {
                return Err(header_error(
                    FORMAT,
                    source_label,
                    Reason::Missing(Part::FrameHeader),
                ))
            }
            Ok(_) => {}
            Err(error) => {
                return Err(io_context(
                    TrajectoryIoOperation::ReadHeader,
                    Some(FORMAT),
                    source_label,
                    error,
                ))
            }
        }
        // The magic is checked as soon as it is read, so a frame that a
        // stream ends inside is rejected when what it holds is already wrong.
        read_exact(
            FORMAT,
            &mut self.reader,
            &mut header[1..4],
            TrajectoryIoOperation::ReadFrame,
            source_label,
            Some(frame_index),
            Reason::Truncated(Part::FrameHeader),
        )?;
        let magic = match i32::from_be_bytes(header[..4].try_into().expect("magic")) {
            1995 => XtcMagic::Xtc1995,
            2023 => XtcMagic::Xtc2023,
            _ => return Err(header_error(FORMAT, source_label, Reason::XtcMagic)),
        };
        read_exact(
            FORMAT,
            &mut self.reader,
            &mut header[4..],
            TrajectoryIoOperation::ReadFrame,
            source_label,
            Some(frame_index),
            Reason::Truncated(Part::FrameHeader),
        )?;
        let atom_count_raw = i32::from_be_bytes(header[4..8].try_into().expect("natoms"));
        let atom_count = nonnegative_xdr_count(
            atom_count_raw,
            Part::AtomCount,
            source_label,
            frame_index,
            start + 4,
        )?;
        validate_atom_count(FORMAT, atom_count, u32::MAX as usize, limits, source_label)?;
        let expected_atoms = topology.atom_count();
        if atom_count != expected_atoms {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
            )
            .with_source_label(source_label)
            .with_frame(frame_index)
            .with_counts(expected_atoms as u64, atom_count as u64));
        }
        let repeated_raw = i32::from_be_bytes(header[52..56].try_into().expect("repeat"));
        let repeated = nonnegative_xdr_count(
            repeated_raw,
            Part::RepeatedAtomCount,
            source_label,
            frame_index,
            start + 52,
        )?;
        if repeated != atom_count {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
            )
            .with_source_label(source_label)
            .with_frame(frame_index)
            .with_counts(atom_count as u64, repeated as u64)
            .with_reason(Reason::RepeatedAtomCountDiffers));
        }
        let step_raw = i32::from_be_bytes(header[8..12].try_into().expect("step"));
        let step = u64::try_from(step_raw).map_err(|_| {
            TrajectoryIoError::codec(
                TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
            )
            .with_source_label(source_label)
            .with_frame(frame_index)
            .with_byte_offset(start + 8)
            .with_reason(Reason::Negative(Part::Step))
        })?;
        let time = f64::from(f32::from_be_bytes(header[12..16].try_into().expect("time")));
        if !time.is_finite() {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame_index,
                Reason::NotFinite(Part::FrameTime),
            ));
        }
        let mut box_values = [0_f64; 9];
        for (value, chunk) in box_values
            .iter_mut()
            .zip(header[16..52].as_chunks::<4>().0.iter())
        {
            *value = f64::from(f32::from_be_bytes(*chunk));
        }
        let cell = decode_cell(box_values, options.cell_policy, source_label, frame_index)?;
        let coordinate_scalars = atom_count.checked_mul(3).ok_or_else(|| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                Reason::AddressRange(Part::CoordinateCount),
            )
        })?;
        let (precision, compressed_bytes, end) = if atom_count <= 9 {
            let coordinate_bytes = atom_count.checked_mul(12).ok_or_else(|| {
                limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadHeader,
                    source_label,
                    None,
                    Reason::AddressRange(Part::FrameSize),
                )
            })?;
            self.reserve_decode_scratch(
                coordinate_scalars,
                coordinate_bytes,
                positions_capacity,
                limits,
                source_label,
                frame_index,
            )?;
            self.scratch.resize(coordinate_bytes, 0);
            read_exact(
                FORMAT,
                &mut self.reader,
                &mut self.scratch,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                Some(frame_index),
                Reason::Truncated(Part::Coordinates),
            )?;
            for chunk in self.scratch.as_chunks::<4>().0 {
                let value = f32::from_be_bytes(*chunk);
                if !value.is_finite() {
                    return Err(frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::CorruptCompressedData,
                        TrajectoryIoOperation::ReadFrame,
                        source_label,
                        frame_index,
                        Reason::NotFinite(Part::Coordinate),
                    ));
                }
                self.decoded_positions.push(value);
            }
            let end = start
                .checked_add(XTC_HEADER_BYTES as u64)
                .and_then(|offset| offset.checked_add(coordinate_bytes as u64))
                .ok_or_else(|| {
                    limit_error(
                        FORMAT,
                        TrajectoryIoOperation::ReadHeader,
                        source_label,
                        None,
                        Reason::Range64(Part::FrameEnd),
                    )
                })?;
            (None, coordinate_bytes, end)
        } else {
            let mut prelude = [0_u8; XTC_COMPRESSED_PRELUDE_BYTES];
            read_exact(
                FORMAT,
                &mut self.reader,
                &mut prelude,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                Some(frame_index),
                Reason::Truncated(Part::CompressedCoordinateHeader),
            )?;
            let precision = f32::from_be_bytes(prelude[..4].try_into().expect("precision"));
            validate_precision(precision, source_label)?;
            let mut minimum = [0_i32; 3];
            let mut maximum = [0_i32; 3];
            for axis in 0..3 {
                let offset = 4 + axis * 4;
                let min =
                    i32::from_be_bytes(prelude[offset..offset + 4].try_into().expect("minimum"));
                let max_offset = 16 + axis * 4;
                let max = i32::from_be_bytes(
                    prelude[max_offset..max_offset + 4]
                        .try_into()
                        .expect("maximum"),
                );
                minimum[axis] = min;
                maximum[axis] = max;
                if min > max {
                    return Err(frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::CorruptCompressedData,
                        TrajectoryIoOperation::ReadFrame,
                        source_label,
                        frame_index,
                        Reason::Compressed(Compressed::BoundsReversed),
                    ));
                }
            }
            let small_index = usize::try_from(u32::from_be_bytes(
                prelude[28..32].try_into().expect("small index"),
            ))
            .map_err(|_| {
                limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    Some(frame_index),
                    Reason::Compressed(Compressed::SmallIndexRange),
                )
            })?;
            if !(XTC_FIRST_SMALL_INDEX..XTC_MAGIC_INTS.len()).contains(&small_index) {
                return Err(frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::CorruptCompressedData,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    frame_index,
                    Reason::Compressed(Compressed::SmallIndexOutsideTable),
                ));
            }
            let compressed_bytes = read_nbytes(&mut self.reader, magic, source_label, frame_index)?;
            if compressed_bytes == 0 {
                return Err(limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    Some(frame_index),
                    Reason::Compressed(Compressed::NoBytes),
                ));
            }
            if compressed_bytes as u64 > limits.max_record_bytes() {
                return Err(limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    Some(frame_index),
                    None,
                )
                .exceeding(
                    Part::CompressedCoordinates,
                    exceeded_bytes(limits.max_record_bytes(), compressed_bytes as u64),
                ));
            }
            if compressed_bytes > limits.max_scratch_bytes() {
                return Err(limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    Some(frame_index),
                    None,
                )
                .exceeding(
                    Part::CompressedCoordinates,
                    LimitExceeded::new(
                        Resource::BufferBytes,
                        limits.max_scratch_bytes(),
                        compressed_bytes,
                    ),
                ));
            }
            self.reserve_decode_scratch(
                coordinate_scalars,
                compressed_bytes,
                positions_capacity,
                limits,
                source_label,
                frame_index,
            )?;
            self.scratch.resize(compressed_bytes, 0);
            read_exact(
                FORMAT,
                &mut self.reader,
                &mut self.scratch,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                Some(frame_index),
                Reason::Truncated(Part::CompressedCoordinates),
            )?;
            decode_compressed_payload(
                &self.scratch,
                CompressedLayout {
                    atom_count,
                    minimum,
                    maximum,
                    initial_small_index: small_index,
                    precision,
                },
                &mut self.decoded_positions,
                source_label,
                frame_index,
            )?;
            let padded = compressed_bytes
                .checked_add((4 - compressed_bytes % 4) % 4)
                .ok_or_else(|| {
                    limit_error(
                        FORMAT,
                        TrajectoryIoOperation::ReadHeader,
                        source_label,
                        None,
                        Reason::AddressRange(Part::Padding),
                    )
                })?;
            let end = start
                .checked_add(XTC_HEADER_BYTES as u64)
                .and_then(|offset| offset.checked_add(XTC_COMPRESSED_PRELUDE_BYTES as u64))
                .and_then(|offset| offset.checked_add(magic.nbytes_width() as u64))
                .and_then(|offset| offset.checked_add(padded as u64))
                .ok_or_else(|| {
                    limit_error(
                        FORMAT,
                        TrajectoryIoOperation::ReadHeader,
                        source_label,
                        None,
                        Reason::Range64(Part::FrameEnd),
                    )
                })?;
            (Some(precision), compressed_bytes, end)
        };
        let frame_bytes = end.checked_sub(start).ok_or_else(|| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                Reason::FrameEndsBeforeStart,
            )
        })?;
        if frame_bytes > limits.max_frame_bytes() {
            return Err(limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                Some(frame_index),
                None,
            )
            .exceeding(
                Part::Frame,
                exceeded_bytes(limits.max_frame_bytes(), frame_bytes),
            ));
        }
        let file_end = self.reader.seek(SeekFrom::End(0)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                error,
            )
        })?;
        if end > file_end {
            self.reader.seek(SeekFrom::Start(start)).map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadHeader,
                    Some(FORMAT),
                    source_label,
                    error,
                )
            })?;
            return Err(stream_error(
                FORMAT,
                TrajectoryIoErrorKind::TruncatedRecord,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                Some(frame_index),
                Reason::FrameEndsAfterFile,
            ));
        }
        if let Some(precision) = precision {
            let padding = (4 - compressed_bytes % 4) % 4;
            if padding > 0 {
                self.reader
                    .seek(SeekFrom::Start(end - padding as u64))
                    .map_err(|error| {
                        io_context(
                            TrajectoryIoOperation::ReadHeader,
                            Some(FORMAT),
                            source_label,
                            error,
                        )
                    })?;
                let mut pad = [0_u8; 3];
                read_exact(
                    FORMAT,
                    &mut self.reader,
                    &mut pad[..padding],
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    Some(frame_index),
                    Reason::Truncated(Part::Padding),
                )?;
                if pad[..padding].iter().any(|byte| *byte != 0) {
                    return Err(frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::CorruptCompressedData,
                        TrajectoryIoOperation::ReadFrame,
                        source_label,
                        frame_index,
                        Reason::NotZero(Part::Padding),
                    ));
                }
            }
            let _ = precision;
        }
        self.reader.seek(SeekFrom::Start(start)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                error,
            )
        })?;
        self.decoded_start = Some(start);
        Ok(Some(XtcFrameInfo {
            start,
            end,
            magic,
            atom_count,
            step,
            time,
            cell,
            precision,
            compressed_bytes,
        }))
    }

    pub(super) fn reserve_decode_scratch(
        &mut self,
        coordinate_scalars: usize,
        raw_bytes: usize,
        positions_capacity: usize,
        limits: &TrajectoryIoLimits,
        source_label: &str,
        frame: usize,
    ) -> Result<(), TrajectoryIoError> {
        // The position array holds one point per atom, whether it is kept or
        // allocated by this read.
        let position_bytes = positions_capacity
            .max(coordinate_scalars / 3)
            .checked_mul(std::mem::size_of::<Point3>());
        let error = |reason: Reason| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                Some(frame),
                reason,
            )
        };
        // Include the forthcoming byte buffer before reserving scalar storage;
        // neither allocation may temporarily escape the aggregate limit.
        reserve_scratch(
            &mut self.decoded_positions,
            coordinate_scalars,
            position_bytes
                .and_then(|bytes| bytes.checked_add(self.scratch.capacity().max(raw_bytes))),
            limits.max_scratch_bytes(),
            error,
        )?;
        reserve_scratch(
            &mut self.scratch,
            raw_bytes,
            position_bytes
                .and_then(|bytes| bytes.checked_add(scratch_bytes(&self.decoded_positions)?)),
            limits.max_scratch_bytes(),
            error,
        )?;
        self.decoded_positions.clear();
        self.scratch.clear();
        Ok(())
    }

    pub(super) fn decode(
        &mut self,
        info: &XtcFrameInfo,
        source_label: &str,
        frame_index: usize,
    ) -> Result<(), TrajectoryIoError> {
        if self.decoded_start != Some(info.start) {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::CorruptCompressedData,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame_index,
                Reason::FrameChanged,
            ));
        }
        self.reader
            .seek(SeekFrom::Start(info.end))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadFrame,
                    Some(FORMAT),
                    source_label,
                    error,
                )
            })?;
        let position = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(FORMAT),
                source_label,
                error,
            )
        })?;
        if position != info.end || self.decoded_positions.len() != info.atom_count.saturating_mul(3)
        {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::CorruptCompressedData,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame_index,
                Reason::DecodedFrameDiffers,
            ));
        }
        if self
            .decoded_positions
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::CorruptCompressedData,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame_index,
                Reason::NotFinite(Part::DecodedCoordinate),
            ));
        }
        self.decoded_start = None;
        Ok(())
    }
}

/// Sequential XTC reader using a private checked decoder for untrusted input.
pub struct XtcReader<R> {
    pub(super) adapter: CheckedXtcReaderAdapter<R>,
    pub(super) topology: Arc<Topology>,
    pub(super) options: XtcReadOptions,
    pub(super) source: ReadSource,
    pub(super) stream_start: u64,
    pub(super) first_info: XtcFrameInfo,
    pub(super) pending_info: Option<XtcFrameInfo>,
    pub(super) positions: Vec<Point3>,
    pub(super) frame_cursor: usize,
}

impl<R> fmt::Debug for XtcReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("XtcReader")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("source_label", &self.source.label)
            .field("frame_cursor", &self.frame_cursor)
            .finish_non_exhaustive()
    }
}

impl<R: Read + Seek> XtcReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order.
    ///
    /// Accepts an owned topology or a shared `Arc<Topology>`. Counts and available
    /// format metadata are checked automatically; matching counts alone cannot
    /// establish atom identity. Limits and diagnostic labels are set in `options`.
    /// The first frame is checked here.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidHeader`] for an empty stream or
    /// malformed header, [`TrajectoryIoErrorKind::UnsupportedVariant`] for an
    /// unsupported magic, [`TrajectoryIoErrorKind::InconsistentAtomCount`] when
    /// its atom count differs from the topology's, a decoding error when the
    /// first frame is invalid, [`TrajectoryIoErrorKind::LimitExceeded`] when
    /// the limits are exceeded, and [`TrajectoryIoErrorKind::Io`] when reading
    /// fails.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: impl Into<TrajectoryOpenOptions>,
    ) -> Result<Self, TrajectoryIoError> {
        let options = options.into();
        let source = options.read_source(FORMAT);
        Self::open(reader, topology.into(), options.xtc, source)
    }

    /// Opens the stream like [`Self::new`] with the format options and
    /// source taken apart.
    pub(crate) fn open(
        reader: R,
        topology: Arc<Topology>,
        options: XtcReadOptions,
        source: ReadSource,
    ) -> Result<Self, TrajectoryIoError> {
        let source_label = &source.label;
        let limits = &source.limits;
        validate_atom_count(
            FORMAT,
            topology.atom_count(),
            u32::MAX as usize,
            limits,
            source_label,
        )?;
        let atom_count = topology.atom_count();
        let mut positions = Vec::new();
        reserve_scratch(
            &mut positions,
            atom_count,
            Some(0),
            limits.max_scratch_bytes(),
            |reason| {
                limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadHeader,
                    source_label,
                    None,
                    reason,
                )
            },
        )?;
        let mut adapter = CheckedXtcReaderAdapter::new(reader);
        let stream_start = adapter.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                source_label,
                error,
            )
        })?;
        let first_info = adapter
            .preflight(&topology, &options, &source, 0, false, positions.capacity())
            .map_err(|error| frame_offset_context(error, 0, stream_start))?
            .ok_or_else(|| header_error(FORMAT, source_label, Reason::EmptyStream))?;
        Ok(Self {
            adapter,
            topology,
            options,
            source,
            stream_start,
            first_info: first_info.clone(),
            pending_info: Some(first_info),
            positions,
            frame_cursor: 0,
        })
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.topology.as_ref()
    }

    pub(crate) fn first_info(&self) -> &XtcFrameInfo {
        &self.first_info
    }

    pub(super) fn next_info(&mut self) -> Result<Option<XtcFrameInfo>, TrajectoryIoError> {
        if let Some(info) = self.pending_info.take() {
            return Ok(Some(info));
        }
        let offset = self.adapter.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                &self.source.label,
                error,
            )
        })?;
        self.adapter
            .preflight(
                &self.topology,
                &self.options,
                &self.source,
                self.frame_cursor,
                true,
                self.positions.capacity(),
            )
            .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))
    }

    pub(super) fn parse_next(&mut self) -> Result<Option<XtcFrameInfo>, TrajectoryIoError> {
        if self.frame_cursor >= self.source.limits.max_frames() {
            let clean_eof = if self.pending_info.is_some() {
                false
            } else {
                probe_seekable_eof(
                    &mut self.adapter.reader,
                    TrajectoryIoOperation::ReadFrame,
                    TrajectoryFormat::Xtc,
                    &self.source.label,
                )?
            };
            if clean_eof {
                return Ok(None);
            }
            return Err(limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.source.label,
                Some(self.frame_cursor),
                None,
            )
            .exceeding(
                Part::Trajectory,
                LimitExceeded::first_beyond(Resource::Frames, self.source.limits.max_frames()),
            ));
        }
        let Some(info) = self.next_info()? else {
            return Ok(None);
        };
        if info.magic != self.first_info.magic || info.precision != self.first_info.precision {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::ReadFrame,
                Some(FORMAT),
            )
            .with_source_label(&self.source.label)
            .with_frame(self.frame_cursor)
            .with_byte_offset(info.start)
            .with_reason(Reason::MagicOrPrecisionDiffers));
        }
        if self.adapter.decoded_start != Some(info.start) {
            self.adapter
                .reader
                .seek(SeekFrom::Start(info.start))
                .map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::ReadFrame,
                        Some(FORMAT),
                        &self.source.label,
                        error,
                    )
                })?;
            let refreshed = self
                .adapter
                .preflight(
                    &self.topology,
                    &self.options,
                    &self.source,
                    self.frame_cursor,
                    false,
                    self.positions.capacity(),
                )?
                .ok_or_else(|| {
                    stream_error(
                        FORMAT,
                        TrajectoryIoErrorKind::TruncatedRecord,
                        TrajectoryIoOperation::ReadFrame,
                        &self.source.label,
                        Some(self.frame_cursor),
                        Reason::IndexedFrameGone,
                    )
                })?;
            if refreshed != info {
                return Err(frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::CorruptCompressedData,
                    TrajectoryIoOperation::ReadFrame,
                    &self.source.label,
                    self.frame_cursor,
                    Reason::FrameChanged,
                ));
            }
        }
        self.adapter
            .decode(&info, &self.source.label, self.frame_cursor)?;
        // Written once. The indexed reader hands the buffer over, so the
        // next read allocates its own on the thread that reads.
        self.positions.clear();
        self.positions
            .reserve_exact(self.adapter.decoded_positions.len() / 3);
        self.positions.extend(
            self.adapter
                .decoded_positions
                .as_chunks::<3>()
                .0
                .iter()
                .map(|values| {
                    Point3::new(
                        f64::from(values[0]),
                        f64::from(values[1]),
                        f64::from(values[2]),
                    )
                }),
        );
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.source.label,
                Some(self.frame_cursor),
                Reason::AddressRange(Part::FrameNumber),
            )
        })?;
        Ok(Some(info))
    }

    pub(super) fn publish(
        &self,
        positions: &[Point3],
        info: &XtcFrameInfo,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        let mut data = FrameData::new(Quantity::new(positions, NANOMETER))
            .with_time(Quantity::new(info.time, PICOSECOND))
            .with_step(info.step);
        if let Some(cell) = info.cell {
            data = data.with_cell(cell);
        }
        destination.replace_from_data(data).map_err(Into::into)
    }

    pub(super) fn position(
        &mut self,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.adapter
            .reader
            .stream_position()
            .map_err(|error| io_context(operation, Some(FORMAT), &self.source.label, error))
    }

    pub(super) fn seek_to(
        &mut self,
        offset: SeekFrom,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.adapter
            .reader
            .seek(offset)
            .map_err(|error| io_context(operation, Some(FORMAT), &self.source.label, error))
    }
}

impl<R: Read + Seek> TrajectoryReader for XtcReader<R> {
    fn topology(&self) -> &Topology {
        self.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        if !self.topology().shares_layout(destination.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        let offset = self.pending_info.as_ref().map_or_else(
            || {
                self.adapter.reader.stream_position().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::ReadFrame,
                        Some(FORMAT),
                        &self.source.label,
                        error,
                    )
                })
            },
            |info| Ok(info.start),
        )?;
        let info = self
            .parse_next()
            .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))?;
        let Some(info) = info else {
            return Ok(false);
        };
        self.publish(&self.positions, &info, destination)?;
        Ok(true)
    }
}

pub(super) struct CheckedByteBuffer<'a> {
    pub(super) bytes: &'a [u8],
    pub(super) index: usize,
}

impl<'a> CheckedByteBuffer<'a> {
    pub(super) const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, index: 0 }
    }

    pub(super) fn pop(
        &mut self,
        source_label: &str,
        frame_index: usize,
    ) -> Result<u8, TrajectoryIoError> {
        let value = self.bytes.get(self.index).copied().ok_or_else(|| {
            frame_error(
                FORMAT,
                TrajectoryIoErrorKind::CorruptCompressedData,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame_index,
                Reason::Compressed(Compressed::EndsBeforeLastAtom),
            )
        })?;
        self.index += 1;
        Ok(value)
    }
}

#[derive(Default)]
pub(super) struct CheckedDecodeState {
    pub(super) last_bits: usize,
    pub(super) last_byte: u8,
}

pub(super) struct CompressedLayout {
    pub(super) atom_count: usize,
    pub(super) minimum: [i32; 3],
    pub(super) maximum: [i32; 3],
    pub(super) initial_small_index: usize,
    pub(super) precision: f32,
}

pub(super) fn decode_compressed_payload(
    payload: &[u8],
    layout: CompressedLayout,
    output: &mut Vec<f32>,
    source_label: &str,
    frame_index: usize,
) -> Result<(), TrajectoryIoError> {
    let CompressedLayout {
        atom_count,
        minimum,
        maximum,
        initial_small_index,
        precision,
    } = layout;
    let mut sizes = [0_u32; 3];
    let mut axis_bits = [0_u32; 3];
    for axis in 0..3 {
        let span = maximum[axis].checked_sub(minimum[axis]).ok_or_else(|| {
            frame_error(
                FORMAT,
                TrajectoryIoErrorKind::CorruptCompressedData,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame_index,
                Reason::Compressed(Compressed::SpanTooWide),
            )
        })?;
        sizes[axis] = u32::try_from(span)
            .ok()
            .and_then(|span| span.checked_add(1))
            .ok_or_else(|| {
                frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::CorruptCompressedData,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    frame_index,
                    Reason::AddressRange(Part::CompressedSize),
                )
            })?;
    }
    let combined_bits = if sizes.iter().fold(0_u32, |bits, size| bits | size) > 0x00ff_ffff {
        for (target, size) in axis_bits.iter_mut().zip(sizes) {
            *target = 32 - size.leading_zeros();
            if *target == 32 {
                return Err(frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::CorruptCompressedData,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    frame_index,
                    Reason::Compressed(Compressed::Axis32Bit),
                ));
            }
        }
        0
    } else {
        combined_bit_count(sizes)
    };
    let mut buffer = CheckedByteBuffer::new(payload);
    let mut state = CheckedDecodeState::default();
    let mut small_index = initial_small_index;
    let mut smaller = XTC_MAGIC_INTS[small_index.saturating_sub(1).max(XTC_FIRST_SMALL_INDEX)] / 2;
    let mut small_number = XTC_MAGIC_INTS[small_index] / 2;
    let mut read_atoms = 0_usize;
    let inverse_precision = precision.recip();
    // A zero update flag preserves the preceding run length. Resetting it for
    // every absolute coordinate desynchronizes ordinary compressed XTC streams.
    let mut run = 0_u32;

    while read_atoms < atom_count {
        let mut coordinate = if combined_bits == 0 {
            let mut decoded = [0_i32; 3];
            for axis in 0..3 {
                let value = decode_checked_bits(
                    &mut buffer,
                    &mut state,
                    axis_bits[axis],
                    source_label,
                    frame_index,
                )?;
                if value >= sizes[axis] || value > i32::MAX as u32 {
                    return Err(frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::CorruptCompressedData,
                        TrajectoryIoOperation::ReadFrame,
                        source_label,
                        frame_index,
                        Reason::Compressed(Compressed::OutsideDeclaredBounds),
                    ));
                }
                decoded[axis] = value as i32;
            }
            decoded
        } else {
            decode_checked_triplet(
                &mut buffer,
                &mut state,
                combined_bits,
                sizes,
                source_label,
                frame_index,
            )?
        };
        for axis in 0..3 {
            coordinate[axis] = minimum[axis]
                .checked_add(coordinate[axis])
                .filter(|value| *value <= maximum[axis])
                .ok_or_else(|| {
                    frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::CorruptCompressedData,
                        TrajectoryIoOperation::ReadFrame,
                        source_label,
                        frame_index,
                        Reason::Compressed(Compressed::ExceedsDeclaredBounds),
                    )
                })?;
        }
        let mut previous = coordinate;
        let flag = decode_checked_bits(&mut buffer, &mut state, 1, source_label, frame_index)? != 0;
        let mut smaller_change = 0_i32;
        if flag {
            run = decode_checked_bits(&mut buffer, &mut state, 5, source_label, frame_index)?;
            let remainder = run % 3;
            run -= remainder;
            smaller_change = remainder as i32 - 1;
        }
        if run == 0 {
            push_checked_coordinate(
                output,
                coordinate,
                inverse_precision,
                source_label,
                frame_index,
            )?;
            read_atoms += 1;
        } else {
            let emitted_atoms = usize::try_from(run / 3)
                .ok()
                .and_then(|count| count.checked_add(1))
                .ok_or_else(|| {
                    frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::CorruptCompressedData,
                        TrajectoryIoOperation::ReadFrame,
                        source_label,
                        frame_index,
                        Reason::AddressRange(Part::CompressedRunAtomCount),
                    )
                })?;
            if emitted_atoms > atom_count - read_atoms {
                return Err(frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::CorruptCompressedData,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    frame_index,
                    Reason::Compressed(Compressed::RunPastAtomCount),
                ));
            }
            let small_size = u32::try_from(XTC_MAGIC_INTS[small_index])
                .expect("audited XTC magic integers are positive");
            let small_sizes = [small_size; 3];
            for offset in (0..run).step_by(3) {
                let delta = decode_checked_triplet(
                    &mut buffer,
                    &mut state,
                    small_index as u32,
                    small_sizes,
                    source_label,
                    frame_index,
                )?;
                for axis in 0..3 {
                    coordinate[axis] = previous[axis]
                        .checked_sub(small_number)
                        .and_then(|base| base.checked_add(delta[axis]))
                        .filter(|value| *value >= minimum[axis] && *value <= maximum[axis])
                        .ok_or_else(|| {
                            frame_error(
                                FORMAT,
                                TrajectoryIoErrorKind::CorruptCompressedData,
                                TrajectoryIoOperation::ReadFrame,
                                source_label,
                                frame_index,
                                Reason::Compressed(Compressed::RunExceedsBounds),
                            )
                        })?;
                }
                if offset == 0 {
                    std::mem::swap(&mut coordinate, &mut previous);
                    push_checked_coordinate(
                        output,
                        previous,
                        inverse_precision,
                        source_label,
                        frame_index,
                    )?;
                    read_atoms += 1;
                } else {
                    previous = coordinate;
                }
                push_checked_coordinate(
                    output,
                    coordinate,
                    inverse_precision,
                    source_label,
                    frame_index,
                )?;
                read_atoms += 1;
            }
        }

        match smaller_change.cmp(&0) {
            std::cmp::Ordering::Less => {
                if small_index <= XTC_FIRST_SMALL_INDEX {
                    return Err(frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::CorruptCompressedData,
                        TrajectoryIoOperation::ReadFrame,
                        source_label,
                        frame_index,
                        Reason::Compressed(Compressed::SmallIndexBelowTable),
                    ));
                }
                small_index -= 1;
                small_number = smaller;
                smaller = if small_index > XTC_FIRST_SMALL_INDEX {
                    XTC_MAGIC_INTS[small_index - 1] / 2
                } else {
                    0
                };
            }
            std::cmp::Ordering::Greater => {
                if small_index + 1 >= XTC_MAGIC_INTS.len() {
                    return Err(frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::CorruptCompressedData,
                        TrajectoryIoOperation::ReadFrame,
                        source_label,
                        frame_index,
                        Reason::Compressed(Compressed::SmallIndexAboveTable),
                    ));
                }
                small_index += 1;
                smaller = small_number;
                small_number = XTC_MAGIC_INTS[small_index] / 2;
            }
            std::cmp::Ordering::Equal => {}
        }
    }
    validate_decoded_count(output, atom_count, source_label, frame_index)?;
    if buffer.index != payload.len() {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::CorruptCompressedData,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame_index,
            Reason::Compressed(Compressed::TrailingBytes),
        ));
    }
    let trailing_mask = if state.last_bits == 0 {
        0
    } else {
        (1_u16 << state.last_bits) - 1
    };
    if u16::from(state.last_byte) & trailing_mask != 0 {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::CorruptCompressedData,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame_index,
            Reason::Compressed(Compressed::TrailingBits),
        ));
    }
    Ok(())
}

pub(super) fn decode_checked_bits(
    buffer: &mut CheckedByteBuffer<'_>,
    state: &mut CheckedDecodeState,
    count: u32,
    source_label: &str,
    frame_index: usize,
) -> Result<u32, TrajectoryIoError> {
    if count > 32 {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::CorruptCompressedData,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame_index,
            Reason::Compressed(Compressed::ValueOver32Bits),
        ));
    }
    let mask = if count == 32 {
        u32::MAX
    } else {
        (1_u32 << count) - 1
    };
    let mut remaining = count as usize;
    let mut last_bits = state.last_bits;
    let mut last_byte = u32::from(state.last_byte);
    let mut value = 0_u32;
    while remaining >= 8 {
        last_byte = last_byte.wrapping_shl(8) | u32::from(buffer.pop(source_label, frame_index)?);
        value |= (last_byte >> last_bits) << (remaining - 8);
        remaining -= 8;
    }
    if remaining > 0 {
        if last_bits < remaining {
            last_bits += 8;
            last_byte =
                last_byte.wrapping_shl(8) | u32::from(buffer.pop(source_label, frame_index)?);
        }
        last_bits -= remaining;
        value |= (last_byte >> last_bits) & mask;
    }
    state.last_bits = last_bits;
    state.last_byte = (last_byte & 0xff) as u8;
    Ok(value & mask)
}

pub(super) fn push_checked_coordinate(
    output: &mut Vec<f32>,
    coordinate: [i32; 3],
    inverse_precision: f32,
    source_label: &str,
    frame_index: usize,
) -> Result<(), TrajectoryIoError> {
    for value in coordinate {
        let value = value as f32 * inverse_precision;
        if !value.is_finite() {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::CorruptCompressedData,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame_index,
                Reason::NotFinite(Part::DecodedCoordinate),
            ));
        }
        output.push(value);
    }
    Ok(())
}

pub(super) fn validate_decoded_count(
    output: &[f32],
    atom_count: usize,
    source_label: &str,
    frame_index: usize,
) -> Result<(), TrajectoryIoError> {
    if output.len() != atom_count.saturating_mul(3) {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::CorruptCompressedData,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame_index,
            Reason::Compressed(Compressed::AtomCountDiffers),
        ));
    }
    Ok(())
}

pub(super) fn combined_bit_count(sizes: [u32; 3]) -> u32 {
    let mut byte_count = 1_usize;
    let mut bytes = [0_u8; 32];
    bytes[0] = 1;
    for size in sizes {
        let mut carry = 0_u32;
        for byte in &mut bytes[..byte_count] {
            carry += u32::from(*byte) * size;
            *byte = (carry & 0xff) as u8;
            carry >>= 8;
        }
        while carry != 0 {
            bytes[byte_count] = (carry & 0xff) as u8;
            byte_count += 1;
            carry >>= 8;
        }
    }
    let high = bytes[byte_count - 1];
    (byte_count as u32 - 1) * 8 + (8 - high.leading_zeros())
}

pub(super) fn decode_checked_triplet(
    buffer: &mut CheckedByteBuffer<'_>,
    state: &mut CheckedDecodeState,
    mut count: u32,
    sizes: [u32; 3],
    source_label: &str,
    frame_index: usize,
) -> Result<[i32; 3], TrajectoryIoError> {
    let mut encoded = [0_u8; 32];
    let mut bytes = 0_usize;
    while count >= 8 {
        encoded[bytes] = decode_checked_bits(buffer, state, 8, source_label, frame_index)? as u8;
        bytes += 1;
        count -= 8;
    }
    if count > 0 {
        encoded[bytes] =
            decode_checked_bits(buffer, state, count, source_label, frame_index)? as u8;
        bytes += 1;
    }
    let mut decoded = [0_u32; 3];
    for axis in (1..=2).rev() {
        let mut remainder = 0_u64;
        for index in (0..bytes).rev() {
            remainder = (remainder << 8) | u64::from(encoded[index]);
            let quotient = remainder / u64::from(sizes[axis]);
            if quotient > u64::from(u8::MAX) {
                return Err(frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::CorruptCompressedData,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    frame_index,
                    Reason::Compressed(Compressed::TripleExceedsBounds),
                ));
            }
            encoded[index] = quotient as u8;
            remainder -= quotient * u64::from(sizes[axis]);
        }
        decoded[axis] = remainder as u32;
    }
    if bytes > 4 && encoded[4..bytes].iter().any(|byte| *byte != 0) {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::CorruptCompressedData,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame_index,
            Reason::Compressed(Compressed::XExceedsBounds),
        ));
    }
    let mut x = [0_u8; 4];
    let copied = bytes.min(x.len());
    x[..copied].copy_from_slice(&encoded[..copied]);
    decoded[0] = u32::from_le_bytes(x);
    for axis in 0..3 {
        if decoded[axis] >= sizes[axis] || decoded[axis] > i32::MAX as u32 {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::CorruptCompressedData,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame_index,
                Reason::Compressed(Compressed::OutsideBounds),
            ));
        }
    }
    Ok(decoded.map(|value| value as i32))
}

pub(super) fn decode_cell(
    values: [f64; 9],
    policy: XtcCellPolicy,
    source_label: &str,
    frame: usize,
) -> Result<Option<PeriodicCell>, TrajectoryIoError> {
    if values.iter().any(|value| !value.is_finite()) {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame,
            Reason::BoxNotFinite,
        ));
    }
    if values.iter().all(|value| *value == 0.0) {
        return match policy {
            XtcCellPolicy::ZeroMatrixAsAbsent => Ok(None),
            XtcCellPolicy::RequirePeriodic => Err(codec_context(
                TrajectoryIoErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::ReadFrame,
                Some(FORMAT),
                source_label,
                Reason::BoxZero,
            )),
        };
    }
    let vectors = [
        Vector3::new(values[0], values[1], values[2]),
        Vector3::new(values[3], values[4], values[5]),
        Vector3::new(values[6], values[7], values[8]),
    ];
    PeriodicCell::new(Quantity::new(vectors, NANOMETER), [true; 3])
        .map(Some)
        .map_err(|error| {
            frame_error(
                FORMAT,
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame,
                Reason::InvalidPeriodicCell,
            )
            .with_source(error)
        })
}

pub(super) fn read_nbytes<R: Read>(
    reader: &mut R,
    magic: XtcMagic,
    source_label: &str,
    frame: usize,
) -> Result<usize, TrajectoryIoError> {
    let mut bytes = [0_u8; 8];
    read_exact(
        FORMAT,
        reader,
        &mut bytes[..magic.nbytes_width()],
        TrajectoryIoOperation::ReadFrame,
        source_label,
        Some(frame),
        Reason::Truncated(Part::CompressedByteCount),
    )?;
    let value = match magic {
        XtcMagic::Xtc1995 => u64::from(u32::from_be_bytes(bytes[..4].try_into().expect("u32"))),
        XtcMagic::Xtc2023 => u64::from_be_bytes(bytes),
    };
    usize::try_from(value).map_err(|_| {
        limit_error(
            FORMAT,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            Some(frame),
            Reason::AddressRange(Part::CompressedByteCount),
        )
    })
}

#[cfg(test)]
mod compressed_run_tests {
    use super::*;

    #[test]
    fn unchanged_run_flag_reuses_the_previous_nonzero_run_length() {
        // Five pairs. Only the first pair supplies a run-length update; later
        // zero flags mean "reuse", not "no run" (xdrfile decompression contract).
        // Pack this focused bitstream independently of the production encoder.
        let mut bits = Vec::new();
        let mut write = |value: u32, count: u32| {
            for bit in (0..count).rev() {
                bits.push((value >> bit) & 1);
            }
        };
        for pair in 0..5 {
            // Absolute mixed-radix [pair, 0, 0], radices [8, 8, 8].
            let packed = pair * 64;
            write(packed & 255, 8);
            write(packed >> 8, 2);
            write(u32::from(pair == 0), 1);
            if pair == 0 {
                write(4, 5);
            } // Run 3; small-index change 0.
              // Delta [5, 4, 4] relative to the absolute coordinate, centered at 4.
            write(356 & 255, 8);
            write(356 >> 8, 1);
        }
        let mut payload = vec![0_u8; bits.len().div_ceil(8)];
        for (index, bit) in bits.into_iter().enumerate() {
            payload[index / 8] |= (bit as u8) << (7 - index % 8);
        }
        let mut output = Vec::new();
        decode_compressed_payload(
            &payload,
            CompressedLayout {
                atom_count: 10,
                minimum: [0; 3],
                maximum: [7; 3],
                initial_small_index: XTC_FIRST_SMALL_INDEX,
                precision: 1.0,
            },
            &mut output,
            "repeated-run",
            0,
        )
        .unwrap();
        let expected = (0..5)
            .flat_map(|pair| [(pair + 1) as f32, 0.0, 0.0, pair as f32, 0.0, 0.0])
            .collect::<Vec<_>>();
        assert_eq!(output, expected);
    }
}
