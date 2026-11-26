use super::*;
use crate::traj::trajectory::reason::{AtomRange, Datum, Field32, Need, Part, Reason, UnitOf};

/// Canonical DCD writer over one seekable stream.
pub struct DcdWriter<W> {
    pub(super) writer: W,
    pub(super) topology: Arc<Topology>,
    pub(super) options: DcdWriteOptions,
    pub(super) source_label: String,
    pub(super) header_start: u64,
    pub(super) frame_count: usize,
    /// Set while a frame's bytes are being written and left set if writing
    /// fails, so a stream that may end in a partial frame is never extended
    /// or reported as finished.
    pub(super) partial: bool,
    pub(super) axis: Vec<f32>,
    pub(super) report: WriteReport,
}

impl<W> fmt::Debug for DcdWriter<W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DcdWriter")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("source_label", &self.source_label)
            .field("frame_count", &self.frame_count)
            .field("partial", &self.partial)
            .finish_non_exhaustive()
    }
}

impl<W: Write + Seek> DcdWriter<W> {
    /// A writer of `topology` frames that writes the header at the stream's
    /// current position.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the step
    /// sequence or `DELTA` does not fit the i32/f32 header or the time-policy
    /// unit is not a time, [`TrajectoryIoErrorKind::LimitExceeded`] when the
    /// atom count is zero or exceeds `i32::MAX` or scratch cannot be reserved,
    /// and [`TrajectoryIoErrorKind::Io`] when writing the header fails.
    pub fn new(
        writer: W,
        topology: impl Into<Arc<Topology>>,
        options: impl Into<TrajectoryWriteOptions>,
    ) -> Result<Self, TrajectoryIoError> {
        let options = options.into();
        let source_label = options.stream_label(FORMAT);
        Self::open(writer, topology.into(), options.dcd, source_label)
    }

    /// Opens the writer like [`Self::new`] with the format options and
    /// label taken apart.
    pub(crate) fn open(
        mut writer: W,
        topology: Arc<Topology>,
        options: DcdWriteOptions,
        source_label: String,
    ) -> Result<Self, TrajectoryIoError> {
        let header = header_records(options, topology.atom_count(), &source_label)?;
        let io = |error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &source_label,
                error,
            )
        };
        let header_start = writer.stream_position().map_err(io)?;
        writer.write_all(&header).map_err(io)?;
        Self::with_frames(writer, topology, options, source_label, header_start, 0)
    }

    /// A writer that appends to a DCD whose `frame_count` frames end at the
    /// stream's position, after the header `options` write starting at
    /// `header_start` ([`check_append_header`]). Nothing is written.
    pub(crate) fn resume(
        writer: W,
        topology: Arc<Topology>,
        options: DcdWriteOptions,
        source_label: String,
        header_start: u64,
        frame_count: usize,
    ) -> Result<Self, TrajectoryIoError> {
        header_records(options, topology.atom_count(), &source_label)?;
        Self::with_frames(
            writer,
            topology,
            options,
            source_label,
            header_start,
            frame_count,
        )
    }

    pub(super) fn with_frames(
        writer: W,
        topology: Arc<Topology>,
        options: DcdWriteOptions,
        source_label: String,
        header_start: u64,
        frame_count: usize,
    ) -> Result<Self, TrajectoryIoError> {
        let mut axis = Vec::new();
        axis.try_reserve_exact(topology.atom_count()).map_err(|_| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::Open,
                &source_label,
                None,
                Reason::Allocation(Part::WriterBuffers),
            )
        })?;
        Ok(Self {
            writer,
            topology,
            options,
            source_label,
            header_start,
            frame_count,
            partial: false,
            axis,
            report: WriteReport::empty(),
        })
    }

    /// The frame state the frames written so far leave out under
    /// [`DcdWriteOptions::geometry_metadata`], [`DcdWriteOptions::cell`], and
    /// [`DcdWriteOptions::properties`], and the steps and times
    /// [`FrameTimingPolicy::Number`] replaced, counted over every frame
    /// written. DCD stores coordinates for a topology its readers receive
    /// separately, so nothing else is reported.
    pub fn report(&self) -> &WriteReport {
        &self.report
    }

    /// The topology every written frame must share.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// The underlying stream.
    pub fn writer(&self) -> &W {
        &self.writer
    }

    /// Writes the frame count to the header's `NSET`, after every frame
    /// written so far has reached the stream, and leaves the stream at its
    /// end.
    pub(super) fn write_frame_count(
        &mut self,
        operation: TrajectoryIoOperation,
    ) -> Result<(), TrajectoryIoError> {
        if self.frame_count > i32::MAX as usize {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                operation,
                Some(FORMAT),
                &self.source_label,
                Reason::FrameCountField,
            ));
        }
        let io = |error| io_context(operation, Some(FORMAT), &self.source_label, error);
        let end = self.writer.stream_position().map_err(io)?;
        let count_offset = self.header_start.checked_add(8).ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                operation,
                Some(FORMAT),
                &self.source_label,
                Reason::Range64(Part::HeaderOffset),
            )
        })?;
        let count = self
            .options
            .endian
            .encode_i32(i32::try_from(self.frame_count).expect("checked frame count"));
        // Frames reach the stream before the count that includes them, so a
        // reader never sees the header count a frame that is not there.
        self.partial = true;
        self.writer
            .flush()
            .and_then(|()| self.writer.seek(SeekFrom::Start(count_offset)))
            .and_then(|_| self.writer.write_all(&count))
            .and_then(|()| self.writer.seek(SeekFrom::Start(end)))
            .map_err(|error| io_context(operation, Some(FORMAT), &self.source_label, error))?;
        self.partial = false;
        Ok(())
    }

    /// Patches the header frame count, flushes, and returns the nonempty DCD
    /// stream, positioned at its end.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidFrame`] when no frame was
    /// written or an earlier write left a partial frame,
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when the frame count exceeds
    /// `i32::MAX`, and [`TrajectoryIoErrorKind::Io`] when seeking, writing,
    /// or flushing fails.
    pub fn finish(mut self) -> Result<W, TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::Finish,
            TrajectoryFormat::Dcd,
            &self.source_label,
        )?;
        require_nonempty_writer(self.frame_count, TrajectoryFormat::Dcd, &self.source_label)?;
        self.write_frame_count(TrajectoryIoOperation::Finish)?;
        self.writer.flush().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Finish,
                Some(FORMAT),
                &self.source_label,
                error,
            )
        })?;
        Ok(self.writer)
    }
}

impl<W: Write + Seek> TrajectoryWriter for DcdWriter<W> {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Writes the header's `NSET` for the frames written so far, after
    /// them, and flushes the stream, which then holds a complete DCD of
    /// those frames.
    fn flush(&mut self) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Dcd,
            &self.source_label,
        )?;
        self.write_frame_count(TrajectoryIoOperation::WriteFrame)?;
        self.writer.flush().map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                error,
            )
        })
    }

    fn write_frame(&mut self, frame: FrameView<'_>) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Dcd,
            &self.source_label,
        )?;
        if !frame.topology().shares_layout(self.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        let omitted = check_frame_state(
            TrajectoryFormat::Dcd,
            &self.source_label,
            self.frame_count,
            frame,
            self.options.state_policies(),
        )?;
        if self.frame_count >= i32::MAX as usize {
            return Err(writer_limit(
                FORMAT,
                &self.source_label,
                Reason::FrameCountField,
            ));
        }
        let expected_step = self
            .options
            .start_step
            .checked_add(
                (self.frame_count as u64)
                    .checked_mul(self.options.step_interval)
                    .ok_or_else(|| writer_limit(FORMAT, &self.source_label, Reason::StepRange64))?,
            )
            .ok_or_else(|| writer_limit(FORMAT, &self.source_label, Reason::StepRange64))?;
        let numbers = self.options.frame_timing == FrameTimingPolicy::Number;
        let replaced_step = match frame.step() {
            Some(step) if step == expected_step => false,
            None if numbers => false,
            Some(_) if numbers => true,
            _ => {
                return Err(codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    &self.source_label,
                    Reason::StepNotInSequence(expected_step),
                ))
            }
        };
        let mut replaced_time = false;
        match (self.options.time_policy, frame.time()) {
            // A time was omitted above or failed there.
            (DcdTimePolicy::Absent, _) => {}
            (DcdTimePolicy::HeaderDelta { unit }, Some(time)) => {
                let actual = time.value_in(unit).map_err(|error| {
                    codec_context(
                        TrajectoryIoErrorKind::InconsistentMetadata,
                        TrajectoryIoOperation::WriteFrame,
                        Some(FORMAT),
                        &self.source_label,
                        Reason::IncompatibleUnit(Some(FORMAT), UnitOf::Of(Datum::Time)),
                    )
                    .with_source(error)
                })?;
                // Readers derive time as step × DELTA, and DELTA is stored
                // as f32. A time therefore matches when it agrees to f32
                // precision; rounding a timestep such as 0.002 to f32
                // already moves it by up to half an f32 epsilon.
                let expected = (expected_step as f64) * f64::from(self.options.delta);
                if !actual.is_finite()
                    || (actual - expected).abs() > f64::from(f32::EPSILON) * expected.abs()
                {
                    if !numbers {
                        return Err(codec_context(
                            TrajectoryIoErrorKind::InconsistentMetadata,
                            TrajectoryIoOperation::WriteFrame,
                            Some(FORMAT),
                            &self.source_label,
                            Reason::TimeNotFromDelta,
                        ));
                    }
                    replaced_time = true;
                }
            }
            (DcdTimePolicy::HeaderDelta { .. }, None) if numbers => {}
            (DcdTimePolicy::HeaderDelta { .. }, None) => {
                return Err(codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    &self.source_label,
                    Reason::FrameLacks(Datum::Time, Need::TimePolicy),
                ))
            }
        }
        // Without unit cells, a cell was omitted above or failed there.
        let cell = frame.cell().copied().filter(|_| self.options.unit_cells);
        if self.options.unit_cells && cell.is_none() {
            return Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                Reason::FrameLacks(Datum::Cell, Need::CellFlag),
            ));
        }
        let factor = CANONICAL_LENGTH_UNIT
            .conversion_factor_to(ANGSTROM)
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    &self.source_label,
                    Reason::IncompatibleUnit(Some(FORMAT), UnitOf::Of(Datum::Position)),
                )
                .with_source(error)
            })?;
        let positions = frame.positions().values();
        for point in *positions.value() {
            for value in [point.x * factor, point.y * factor, point.z * factor] {
                if !value.is_finite() || !(value as f32).is_finite() {
                    return Err(codec_context(
                        TrajectoryIoErrorKind::InvalidFrame,
                        TrajectoryIoOperation::WriteFrame,
                        Some(FORMAT),
                        &self.source_label,
                        Reason::DoesNotFitSingle(Part::Coordinate),
                    ));
                }
            }
        }
        // Everything that can reject the frame is checked above or here, so a
        // rejected frame leaves the stream untouched.
        let cell_payload = cell
            .map(|cell| -> Result<[u8; CELL_BYTES], TrajectoryIoError> {
                let values = encode_cell(cell, &self.source_label)?;
                let mut payload = [0_u8; CELL_BYTES];
                for (chunk, value) in payload.as_chunks_mut::<8>().0.iter_mut().zip(values) {
                    chunk.copy_from_slice(&self.options.endian.encode_f64(value));
                }
                Ok(payload)
            })
            .transpose()?;
        let byte_len = positions
            .value()
            .len()
            .checked_mul(4)
            .and_then(|len| u32::try_from(len).ok())
            .ok_or_else(|| {
                writer_limit(
                    FORMAT,
                    &self.source_label,
                    Reason::Exceeds32Bit(Part::CoordinateRecord, Field32::RecordLength),
                )
            })?;
        self.partial = true;
        if let Some(payload) = cell_payload {
            write_record(
                &mut self.writer,
                self.options.endian,
                &payload,
                &self.source_label,
                TrajectoryIoOperation::WriteFrame,
            )?;
        }
        for axis in 0..3 {
            self.axis.clear();
            for point in *positions.value() {
                let value = match axis {
                    0 => point.x,
                    1 => point.y,
                    _ => point.z,
                } * factor;
                self.axis.push(value as f32);
            }
            self.writer
                .write_all(&self.options.endian.encode_u32(byte_len))
                .and_then(|_| {
                    for value in &self.axis {
                        self.writer
                            .write_all(&self.options.endian.encode_f32(*value))?;
                    }
                    self.writer
                        .write_all(&self.options.endian.encode_u32(byte_len))
                })
                .map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::WriteFrame,
                        Some(FORMAT),
                        &self.source_label,
                        error,
                    )
                })?;
        }
        self.partial = false;
        self.frame_count = self.frame_count.checked_add(1).ok_or_else(|| {
            writer_limit(
                FORMAT,
                &self.source_label,
                Reason::AddressRange(Part::FrameCount),
            )
        })?;
        record(&mut self.report, omitted);
        record(
            &mut self.report,
            [
                (replaced_step, GeometryField::Step),
                (replaced_time, GeometryField::Time),
            ]
            .map(|(replaced, field)| WriteOmission::GeometryState {
                field,
                geometries: usize::from(replaced),
            }),
        );
        Ok(())
    }
}

/// The header, title, and atom-count records a writer with `options` writes
/// for `atom_count` atoms, counting no frames.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the step
/// sequence or `DELTA` does not fit the i32/f32 header or the time-policy
/// unit is not a time, and [`TrajectoryIoErrorKind::LimitExceeded`] when the
/// atom count is zero or exceeds `i32::MAX`.
pub(super) fn header_records(
    options: DcdWriteOptions,
    atom_count: usize,
    source_label: &str,
) -> Result<Vec<u8>, TrajectoryIoError> {
    if atom_count == 0 || atom_count > i32::MAX as usize {
        return Err(codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::Open,
            Some(FORMAT),
            source_label,
            Reason::AtomCountOutOfRange(atom_count, AtomRange::DcdHeader),
        ));
    }
    let mut records = Vec::new();
    let mut title = vec![0_u8; 84];
    title[..4].copy_from_slice(&options.endian.encode_i32(1));
    let text = b"Created by kekule";
    title[4..4 + text.len()].copy_from_slice(text);
    let atoms = options
        .endian
        .encode_i32(i32::try_from(atom_count).expect("checked atom count"));
    for payload in [&encode_header(options, 0)[..], &title, &atoms] {
        write_record(
            &mut records,
            options.endian,
            payload,
            source_label,
            TrajectoryIoOperation::Open,
        )?;
    }
    Ok(records)
}

/// Fails unless `stream`, at a DCD's start, holds the header, title, and
/// atom-count records a writer with `options` writes for `atom_count` atoms,
/// whatever frame count `NSET` they hold, so that appending with `options`
/// continues the file as one writer would have written it.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the records
/// differ: another endianness, step sequence, `DELTA`, unit-cell flag,
/// title, or atom count; the errors of the options' validation; and
/// [`TrajectoryIoErrorKind::Io`] when the records cannot be read.
pub(crate) fn check_append_header(
    stream: &mut impl Read,
    options: DcdWriteOptions,
    atom_count: usize,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    let expected = header_records(options, atom_count, source_label)?;
    let mut found = vec![0_u8; expected.len()];
    stream.read_exact(&mut found).map_err(|error| {
        io_context(
            TrajectoryIoOperation::Open,
            Some(FORMAT),
            source_label,
            error,
        )
    })?;
    // NSET, after the header record's marker and signature, counts frames.
    found[8..12].copy_from_slice(&expected[8..12]);
    if found != expected {
        return Err(codec_context(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::Open,
            Some(FORMAT),
            source_label,
            Reason::AppendHeaderDiffers,
        ));
    }
    Ok(())
}

pub(super) fn encode_header(options: DcdWriteOptions, frames: i32) -> [u8; HEADER_BYTES] {
    let mut header = [0_u8; HEADER_BYTES];
    header[..4].copy_from_slice(b"CORD");
    let controls = [
        frames,
        options.start_step as i32,
        options.step_interval as i32,
        0,
        0,
        0,
        0,
        0,
        0,
        i32::from_ne_bytes(options.delta.to_bits().to_ne_bytes()),
        i32::from(options.unit_cells),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        24,
    ];
    for (chunk, control) in header[4..].as_chunks_mut::<4>().0.iter_mut().zip(controls) {
        chunk.copy_from_slice(&options.endian.encode_i32(control));
    }
    header[40..44].copy_from_slice(&options.endian.encode_f32(options.delta));
    header
}

pub(super) fn encode_cell(
    cell: PeriodicCell,
    source_label: &str,
) -> Result<[f64; 6], TrajectoryIoError> {
    if cell.periodic_axes() != [true; 3] {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::WriteFrame,
            Some(FORMAT),
            source_label,
            Reason::NotFullyPeriodic(FORMAT),
        ));
    }
    let [a, b, c] = cell.vectors().value_in(ANGSTROM).map_err(|error| {
        codec_context(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::WriteFrame,
            Some(FORMAT),
            source_label,
            Reason::IncompatibleUnit(Some(FORMAT), UnitOf::Of(Datum::Cell)),
        )
        .with_source(error)
    })?;
    let [a, b, c] = standard_orientation([a, b, c], FORMAT, source_label)?;
    let length =
        |vector: Vector3| (vector.x * vector.x + vector.y * vector.y + vector.z * vector.z).sqrt();
    let dot =
        |left: Vector3, right: Vector3| left.x * right.x + left.y * right.y + left.z * right.z;
    let la = length(a);
    let lb = length(b);
    let lc = length(c);
    let cosine = |left: Vector3, right: Vector3, ll: f64, lr: f64| {
        (dot(left, right) / (ll * lr)).clamp(-1.0, 1.0)
    };
    let alpha = cosine(b, c, lb, lc);
    let beta = cosine(a, c, la, lc);
    let gamma = cosine(a, b, la, lb);
    let values = [la, gamma, lb, beta, alpha, lc];
    if values.iter().any(|value| !value.is_finite()) {
        return Err(codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(FORMAT),
            source_label,
            Reason::CellNoParameters,
        ));
    }
    let mut record = [0; CELL_BYTES];
    for (chunk, value) in record.as_chunks_mut::<8>().0.iter_mut().zip(values) {
        *chunk = value.to_le_bytes();
    }
    decode_cell(&record, DcdEndian::Little, source_label, 0).map_err(|error| {
        codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(FORMAT),
            source_label,
            Reason::CellDegenerateParameters,
        )
        .with_source(validation_cause(error))
    })?;
    Ok(values)
}

pub(super) fn write_record<W: Write>(
    writer: &mut W,
    endian: DcdEndian,
    payload: &[u8],
    source_label: &str,
    operation: TrajectoryIoOperation,
) -> Result<(), TrajectoryIoError> {
    let size = u32::try_from(payload.len()).map_err(|_| {
        writer_limit(
            FORMAT,
            source_label,
            Reason::Exceeds32Bit(Part::Record, Field32::RecordLength),
        )
    })?;
    writer
        .write_all(&endian.encode_u32(size))
        .and_then(|_| writer.write_all(payload))
        .and_then(|_| writer.write_all(&endian.encode_u32(size)))
        .map_err(|error| io_context(operation, Some(FORMAT), source_label, error))
}
