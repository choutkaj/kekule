use super::*;
use crate::traj::trajectory::reason::{AtomRange, Datum, Field32, Need, Part, Reason, UnitOf};

/// Pure-Rust TRR writer.
pub struct TrrWriter<W> {
    pub(super) writer: W,
    pub(super) topology: Arc<Topology>,
    pub(super) options: TrrWriteOptions,
    pub(super) source_label: String,
    pub(super) frame_count: usize,
    /// Set while a frame's bytes are being written and left set if writing
    /// fails, so a stream that may end in a partial frame is never extended
    /// or reported as finished.
    pub(super) partial: bool,
    pub(super) raw: Vec<u8>,
    pub(super) report: WriteReport,
}

impl<W> fmt::Debug for TrrWriter<W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("TrrWriter")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("source_label", &self.source_label)
            .field("frame_count", &self.frame_count)
            .field("partial", &self.partial)
            .finish_non_exhaustive()
    }
}

impl<W: Write> TrrWriter<W> {
    /// A writer of `topology` frames to `writer`. Nothing is written until the first frame.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::LimitExceeded`] when the topology has
    /// no atoms or more than `i32::MAX / 3`.
    pub fn new(
        writer: W,
        topology: impl Into<Arc<Topology>>,
        options: impl Into<TrajectoryWriteOptions>,
    ) -> Result<Self, TrajectoryIoError> {
        let options = options.into();
        let source_label = options.stream_label(FORMAT);
        Self::open(writer, topology.into(), options.trr, source_label)
    }

    /// Opens the writer like [`Self::new`] with the format options and
    /// label taken apart.
    pub(crate) fn open(
        writer: W,
        topology: Arc<Topology>,
        options: TrrWriteOptions,
        source_label: String,
    ) -> Result<Self, TrajectoryIoError> {
        if topology.atom_count() == 0 || topology.atom_count() > i32::MAX as usize / 3 {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &source_label,
                Reason::AtomCountOutOfRange(topology.atom_count(), AtomRange::Trr),
            ));
        }
        Ok(Self {
            writer,
            topology,
            options,
            source_label,
            frame_count: 0,
            partial: false,
            raw: Vec::new(),
            report: WriteReport::empty(),
        })
    }

    /// The frame state the frames written so far leave out under
    /// [`TrrWriteOptions::geometry_metadata`] and
    /// [`TrrWriteOptions::properties`], counted over every frame
    /// written. TRR stores coordinates for a topology its readers receive
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

    /// Takes the stream as holding `frame_count` frames already, which the
    /// next frame follows, as when appending to a file.
    pub(crate) const fn resumed(mut self, frame_count: usize) -> Self {
        self.frame_count = frame_count;
        self
    }

    pub(crate) fn validate_finish(&self) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::Finish,
            TrajectoryFormat::Trr,
            &self.source_label,
        )?;
        require_nonempty_writer(self.frame_count, TrajectoryFormat::Trr, &self.source_label)
    }

    /// Flushes and returns the completed nonempty TRR stream.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidFrame`] when no frame was
    /// written or an earlier write left a partial frame, and
    /// [`TrajectoryIoErrorKind::Io`] when flushing fails.
    pub fn finish(mut self) -> Result<W, TrajectoryIoError> {
        self.validate_finish()?;
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

impl<W: Write> TrajectoryWriter for TrrWriter<W> {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn write_frame(&mut self, frame: FrameView<'_>) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Trr,
            &self.source_label,
        )?;
        if !frame.topology().shares_layout(self.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        let properties = frame.properties();
        let mut omitted = check_frame_state(
            TrajectoryFormat::Trr,
            &self.source_label,
            self.frame_count,
            frame,
            self.options.state_policies(),
        )?;
        if self.options.properties == UnrepresentablePolicy::Omit
            && has_unwritten_owner_properties(properties.owner(), self.options.lambda_policy)
        {
            omitted.push(WriteOmission::GeometryState {
                field: GeometryField::Properties,
                geometries: 1,
            });
        }
        if self.frame_count == usize::MAX {
            return Err(writer_limit(
                FORMAT,
                &self.source_label,
                Reason::AddressRange(Part::FrameCount),
            ));
        }
        let (step, time) = self.options.frame_timing.timing(
            self.frame_count,
            self.options.frame_interval,
            frame.step(),
            frame.time(),
        );
        let step = step.ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                Reason::FrameLacks(Datum::Step, Need::Format(FORMAT)),
            )
        })?;
        let step = i32::try_from(step).map_err(|_| {
            codec_context(
                TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                Reason::Exceeds32Bit(Part::FrameStep, Field32::StepField),
            )
        })?;
        let time = time.ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                Reason::FrameLacks(Datum::Time, Need::Format(FORMAT)),
            )
        })?;
        let time = time.value_in(PICOSECOND).map_err(|error| {
            codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                Reason::IncompatibleUnit(Some(FORMAT), UnitOf::Of(Datum::Time)),
            )
            .with_source(error)
        })?;
        let lambda = writer_lambda(
            properties,
            self.options.lambda_policy,
            self.options.properties,
            &self.source_label,
            self.frame_count,
        )?;
        let scalar_bytes = self.options.precision.bytes();
        validate_representable(
            time,
            self.options.precision,
            &self.source_label,
            Datum::Time,
        )?;
        validate_representable(
            lambda,
            self.options.precision,
            &self.source_label,
            Datum::Lambda,
        )?;
        let atom_scalars = self.topology.atom_count().checked_mul(3).ok_or_else(|| {
            writer_limit(
                FORMAT,
                &self.source_label,
                Reason::AddressRange(Part::FrameValues),
            )
        })?;
        let vector_bytes = atom_scalars.checked_mul(scalar_bytes).ok_or_else(|| {
            writer_limit(
                FORMAT,
                &self.source_label,
                Reason::AddressRange(Part::FrameVectors),
            )
        })?;
        let box_size = if frame.cell().is_some() {
            9 * scalar_bytes
        } else {
            0
        };
        let v_size = if frame.velocities().is_some() {
            vector_bytes
        } else {
            0
        };
        let f_size = if frame.forces().is_some() {
            vector_bytes
        } else {
            0
        };
        let cell = frame.cell().copied();
        let cell_data = if let Some(cell) = cell {
            if cell.periodic_axes() != [true; 3] {
                return Err(codec_context(
                    TrajectoryIoErrorKind::UnsupportedVariant,
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    &self.source_label,
                    Reason::NotFullyPeriodic(FORMAT),
                ));
            }
            let factor = CANONICAL_LENGTH_UNIT
                .conversion_factor_to(NANOMETER)
                .map_err(|error| writer_unit(&self.source_label, UnitOf::Of(Datum::Cell), error))?;
            Some((cell.vectors().into_value(), factor))
        } else {
            None
        };
        let positions = frame.positions().values();
        let position_factor = CANONICAL_LENGTH_UNIT
            .conversion_factor_to(NANOMETER)
            .map_err(|error| {
                writer_unit(&self.source_label, UnitOf::AllOf(Datum::Position), error)
            })?;
        let velocity_factor = frame
            .velocities()
            .map(|velocities| {
                velocities
                    .values()
                    .unit()
                    .conversion_factor_to(NANOMETER_PER_PICOSECOND)
                    .map_err(|error| {
                        writer_unit(&self.source_label, UnitOf::AllOf(Datum::Velocity), error)
                    })
            })
            .transpose()?;
        let force_factor = frame
            .forces()
            .map(|forces| {
                forces
                    .values()
                    .unit()
                    .conversion_factor_to(KILOJOULE_PER_MOLE_PER_NANOMETER)
                    .map_err(|error| {
                        writer_unit(&self.source_label, UnitOf::AllOf(Datum::Force), error)
                    })
            })
            .transpose()?;

        if let Some((vectors, factor)) = &cell_data {
            encode_vectors_to_raw(
                &mut self.raw,
                vectors,
                *factor,
                self.options.precision,
                &self.source_label,
                Datum::Cell,
            )?;
            decode_cell(
                &self.raw,
                self.options.precision,
                &self.source_label,
                self.frame_count,
            )
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InvalidFrame,
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    &self.source_label,
                    Reason::DoesNotFitSelected(Part::Cell),
                )
                .with_source(validation_cause(error))
            })?;
        }
        encode_points_to_raw(
            &mut self.raw,
            positions.value(),
            position_factor,
            self.options.precision,
            &self.source_label,
        )?;
        if let (Some(velocities), Some(factor)) = (
            frame.velocities().map(|values| values.values()),
            velocity_factor,
        ) {
            encode_vectors_to_raw(
                &mut self.raw,
                velocities.value(),
                factor,
                self.options.precision,
                &self.source_label,
                Datum::Velocity,
            )?;
        }
        if let (Some(forces), Some(factor)) =
            (frame.forces().map(|values| values.values()), force_factor)
        {
            encode_vectors_to_raw(
                &mut self.raw,
                forces.value(),
                factor,
                self.options.precision,
                &self.source_label,
                Datum::Force,
            )?;
        }

        let mut block_sizes = [0_i32; 10];
        for (target, size) in
            block_sizes
                .iter_mut()
                .zip([0, 0, box_size, 0, 0, 0, 0, vector_bytes, v_size, f_size])
        {
            *target = i32::try_from(size).map_err(|_| {
                writer_limit(
                    FORMAT,
                    &self.source_label,
                    Reason::Exceeds32Bit(Part::Block, Field32::BlockSize),
                )
            })?;
        }
        // Every scalar and size was checked above, so a rejected frame leaves
        // the stream untouched; from here only I/O can fail.
        self.partial = true;
        write_i32(&mut self.writer, TRR_MAGIC, &self.source_label)?;
        write_i32(
            &mut self.writer,
            i32::try_from(TRR_VERSION.len() + 1).expect("TRR version length fits i32"),
            &self.source_label,
        )?;
        write_xdr_string(&mut self.writer, TRR_VERSION, &self.source_label)?;
        for size in block_sizes {
            write_i32(&mut self.writer, size, &self.source_label)?;
        }
        write_i32(
            &mut self.writer,
            i32::try_from(self.topology.atom_count()).expect("checked atom count"),
            &self.source_label,
        )?;
        write_i32(&mut self.writer, step, &self.source_label)?;
        write_i32(&mut self.writer, 0, &self.source_label)?;
        write_scalar(
            &mut self.writer,
            time,
            self.options.precision,
            &self.source_label,
        )?;
        write_scalar(
            &mut self.writer,
            lambda,
            self.options.precision,
            &self.source_label,
        )?;
        if let Some((vectors, factor)) = cell_data {
            encode_vectors_to_raw(
                &mut self.raw,
                &vectors,
                factor,
                self.options.precision,
                &self.source_label,
                Datum::Cell,
            )?;
            write_bytes(&mut self.writer, &self.raw, &self.source_label)?;
        }
        encode_points_to_raw(
            &mut self.raw,
            positions.value(),
            position_factor,
            self.options.precision,
            &self.source_label,
        )?;
        write_bytes(&mut self.writer, &self.raw, &self.source_label)?;
        if let (Some(velocities), Some(factor)) = (
            frame.velocities().map(|values| values.values()),
            velocity_factor,
        ) {
            encode_vectors_to_raw(
                &mut self.raw,
                velocities.value(),
                factor,
                self.options.precision,
                &self.source_label,
                Datum::Velocity,
            )?;
            write_bytes(&mut self.writer, &self.raw, &self.source_label)?;
        }
        if let (Some(forces), Some(factor)) =
            (frame.forces().map(|values| values.values()), force_factor)
        {
            encode_vectors_to_raw(
                &mut self.raw,
                forces.value(),
                factor,
                self.options.precision,
                &self.source_label,
                Datum::Force,
            )?;
            write_bytes(&mut self.writer, &self.raw, &self.source_label)?;
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
        Ok(())
    }

    /// Flushes the stream, which then ends after the last whole frame.
    fn flush(&mut self) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Trr,
            &self.source_label,
        )?;
        self.writer.flush().map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                error,
            )
        })
    }
}

/// Fails unless a writer with `options` writes frames of the scalar width of
/// the TRR frames whose first is `first`, which must not mix widths
/// (`mixed`).
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the frames
/// mix widths or have another width than `options`.
pub(crate) fn check_append_frames(
    first: &TrrFrameHeader,
    mixed: bool,
    options: &TrrWriteOptions,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    let reason = if mixed {
        Reason::AppendMixedTrrWidths
    } else if first.precision != options.precision {
        Reason::AppendTrrPrecision {
            file: first.precision.name(),
            options: options.precision.name(),
        }
    } else {
        return Ok(());
    };
    Err(codec_context(
        TrajectoryIoErrorKind::InconsistentMetadata,
        TrajectoryIoOperation::Open,
        Some(FORMAT),
        source_label,
        reason,
    ))
}

pub(super) fn encode_points_to_raw(
    raw: &mut Vec<u8>,
    points: &[Point3],
    factor: f64,
    precision: TrrScalarPrecision,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    encode_triplets_to_raw(
        raw,
        points.iter().map(|p| [p.x, p.y, p.z]),
        factor,
        precision,
        source_label,
        Datum::Position,
        Reason::AddressRange(Part::Positions),
    )
}

pub(super) fn encode_vectors_to_raw(
    raw: &mut Vec<u8>,
    vectors: &[Vector3],
    factor: f64,
    precision: TrrScalarPrecision,
    source_label: &str,
    field: Datum,
) -> Result<(), TrajectoryIoError> {
    encode_triplets_to_raw(
        raw,
        vectors.iter().map(|v| [v.x, v.y, v.z]),
        factor,
        precision,
        source_label,
        field,
        Reason::AddressRange(Part::Vectors),
    )
}

pub(super) fn encode_triplets_to_raw(
    raw: &mut Vec<u8>,
    values: impl ExactSizeIterator<Item = [f64; 3]>,
    factor: f64,
    precision: TrrScalarPrecision,
    source_label: &str,
    field: Datum,
    overflow: Reason,
) -> Result<(), TrajectoryIoError> {
    raw.clear();
    let bytes = values
        .len()
        .checked_mul(3 * precision.bytes())
        .ok_or_else(|| writer_limit(FORMAT, source_label, overflow))?;
    raw.try_reserve(bytes.saturating_sub(raw.len()))
        .map_err(|_| {
            writer_limit(
                FORMAT,
                source_label,
                Reason::Allocation(Part::WriterBuffers),
            )
        })?;
    for triplet in values {
        for value in triplet {
            push_scalar(raw, value * factor, precision, source_label, field)?;
        }
    }
    Ok(())
}

pub(super) fn push_scalar(
    bytes: &mut Vec<u8>,
    value: f64,
    precision: TrrScalarPrecision,
    source_label: &str,
    field: Datum,
) -> Result<(), TrajectoryIoError> {
    validate_representable(value, precision, source_label, field)?;
    match precision {
        TrrScalarPrecision::Float32 => bytes.extend_from_slice(&(value as f32).to_be_bytes()),
        TrrScalarPrecision::Float64 => bytes.extend_from_slice(&value.to_be_bytes()),
    }
    Ok(())
}

pub(super) fn validate_representable(
    value: f64,
    precision: TrrScalarPrecision,
    source_label: &str,
    field: Datum,
) -> Result<(), TrajectoryIoError> {
    if !value.is_finite()
        || (precision == TrrScalarPrecision::Float32 && !(value as f32).is_finite())
    {
        return Err(codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(FORMAT),
            source_label,
            Reason::DoesNotFitSelected(Part::ValueOf(field)),
        ));
    }
    Ok(())
}

/// Whether `owner` holds frame properties a frame written under `policy`
/// leaves out: any under [`TrrLambdaPolicy::RequireZero`], and any but the
/// lambda under [`TrrLambdaPolicy::FrameProperty`].
pub(super) fn has_unwritten_owner_properties(
    owner: &OwnerProperties,
    policy: TrrLambdaPolicy,
) -> bool {
    match policy {
        TrrLambdaPolicy::FrameProperty => {
            owner.len() > usize::from(owner.get(&trr_lambda_key()).is_some())
        }
        TrrLambdaPolicy::RequireZero => !owner.is_empty(),
    }
}

/// The lambda of a frame with `properties`. Under
/// [`UnrepresentablePolicy::Omit`] for `properties`, properties
/// [`has_unwritten_owner_properties`] names are ignored, and atom and bond
/// properties are left to the caller.
pub(super) fn writer_lambda(
    properties: &GeometryProperties,
    policy: TrrLambdaPolicy,
    properties_policy: UnrepresentablePolicy,
    source_label: &str,
    frame: usize,
) -> Result<f64, TrajectoryIoError> {
    let omit = properties_policy == UnrepresentablePolicy::Omit;
    match policy {
        TrrLambdaPolicy::FrameProperty => {
            let other = if omit {
                properties.owner().get(&trr_lambda_key()).is_none()
            } else {
                properties.owner().len() != 1
                    || properties.atoms().has_data()
                    || properties.bonds().has_data()
            };
            if other {
                return Err(codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    source_label,
                    Reason::LambdaPolicyKey,
                ));
            }
            match properties.owner().get(&trr_lambda_key()) {
                Some(PropertyValue::Real { value, unit })
                    if value.is_finite() && *unit == DIMENSIONLESS =>
                {
                    Ok(*value)
                }
                _ => Err(codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    source_label,
                    Reason::LambdaPropertyValue,
                )),
            }
        }
        TrrLambdaPolicy::RequireZero => {
            if !omit && !properties.is_empty() {
                return Err(writer_field_error(
                    TrajectoryFormat::Trr,
                    source_label,
                    frame,
                    GeometryField::Properties,
                ));
            }
            Ok(0.0)
        }
    }
}

pub(super) fn write_i32<W: Write>(
    writer: &mut W,
    value: i32,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    write_bytes(writer, &value.to_be_bytes(), source_label)
}

pub(super) fn write_xdr_string<W: Write>(
    writer: &mut W,
    value: &[u8],
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    write_i32(
        writer,
        i32::try_from(value.len()).map_err(|_| {
            writer_limit(
                FORMAT,
                source_label,
                Reason::Exceeds32Bit(Part::VersionString, Field32::LengthField),
            )
        })?,
        source_label,
    )?;
    write_bytes(writer, value, source_label)?;
    let padding = (4 - value.len() % 4) % 4;
    write_bytes(writer, &[0_u8; 3][..padding], source_label)
}

pub(super) fn write_scalar<W: Write>(
    writer: &mut W,
    value: f64,
    precision: TrrScalarPrecision,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    match precision {
        TrrScalarPrecision::Float32 => {
            write_bytes(writer, &(value as f32).to_be_bytes(), source_label)
        }
        TrrScalarPrecision::Float64 => write_bytes(writer, &value.to_be_bytes(), source_label),
    }
}

pub(super) fn write_bytes<W: Write>(
    writer: &mut W,
    bytes: &[u8],
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    writer.write_all(bytes).map_err(|error| {
        io_context(
            TrajectoryIoOperation::WriteFrame,
            Some(FORMAT),
            source_label,
            error,
        )
    })
}
