use super::*;

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
    /// A writer of `topology` frames to `writer`; `source_label` names the
    /// destination in errors. Nothing is written until the first frame.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::LimitExceeded`] when the topology has
    /// no atoms or more than `i32::MAX / 3`.
    pub fn new(
        writer: W,
        topology: impl Into<Arc<Topology>>,
        options: TrrWriteOptions,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryIoError> {
        let topology = topology.into();
        let source_label = source_label.into();
        if topology.atom_count() == 0 || topology.atom_count() > i32::MAX as usize / 3 {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Trr),
                &source_label,
                format!(
                    "The topology has {} atoms. TRR holds 1 to 715,827,882.",
                    topology.atom_count()
                ),
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
                Some(TrajectoryFormat::Trr),
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
                &self.source_label,
                "The frame count exceeds this platform's address range.",
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
                Some(TrajectoryFormat::Trr),
                &self.source_label,
                "The frame has no step. TRR frames need one.",
            )
        })?;
        let step = i32::try_from(step).map_err(|_| {
            codec_context(
                TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Trr),
                &self.source_label,
                "The frame's step exceeds the 32-bit step field.",
            )
        })?;
        let time = time.ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Trr),
                &self.source_label,
                "The frame has no time. TRR frames need one.",
            )
        })?;
        let time = time.value_in(PICOSECOND).map_err(|error| {
            codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Trr),
                &self.source_label,
                "The TRR time unit is incompatible.",
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
        validate_representable(time, self.options.precision, &self.source_label, "time")?;
        validate_representable(lambda, self.options.precision, &self.source_label, "lambda")?;
        let atom_scalars = self.topology.atom_count().checked_mul(3).ok_or_else(|| {
            writer_limit(
                &self.source_label,
                "The frame's values exceed this platform's address range.",
            )
        })?;
        let vector_bytes = atom_scalars.checked_mul(scalar_bytes).ok_or_else(|| {
            writer_limit(
                &self.source_label,
                "The frame's vectors exceed this platform's address range.",
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
                    Some(TrajectoryFormat::Trr),
                    &self.source_label,
                    "The cell is not periodic along all three axes. TRR stores only fully periodic cells.",
                ));
            }
            let factor = CANONICAL_LENGTH_UNIT
                .conversion_factor_to(NANOMETER)
                .map_err(|error| writer_unit(&self.source_label, "cell", error))?;
            Some((cell.vectors().into_value(), factor))
        } else {
            None
        };
        let positions = frame.positions().values();
        let position_factor = CANONICAL_LENGTH_UNIT
            .conversion_factor_to(NANOMETER)
            .map_err(|error| writer_unit(&self.source_label, "positions", error))?;
        let velocity_factor = frame
            .velocities()
            .map(|velocities| {
                velocities
                    .values()
                    .unit()
                    .conversion_factor_to(NANOMETER_PER_PICOSECOND)
                    .map_err(|error| writer_unit(&self.source_label, "velocities", error))
            })
            .transpose()?;
        let force_factor = frame
            .forces()
            .map(|forces| {
                forces
                    .values()
                    .unit()
                    .conversion_factor_to(KILOJOULE_PER_MOLE_PER_NANOMETER)
                    .map_err(|error| writer_unit(&self.source_label, "forces", error))
            })
            .transpose()?;

        if let Some((vectors, factor)) = &cell_data {
            encode_vectors_to_raw(
                &mut self.raw,
                vectors,
                *factor,
                self.options.precision,
                &self.source_label,
                "cell",
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
                    Some(TrajectoryFormat::Trr),
                    &self.source_label,
                    "The cell does not fit the selected precision.",
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
                "velocity",
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
                "force",
            )?;
        }

        let mut block_sizes = [0_i32; 10];
        for (target, size) in
            block_sizes
                .iter_mut()
                .zip([0, 0, box_size, 0, 0, 0, 0, vector_bytes, v_size, f_size])
        {
            *target = i32::try_from(size).map_err(|_| {
                writer_limit(&self.source_label, "A block exceeds the 32-bit block size.")
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
                "cell",
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
                "velocity",
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
                "force",
            )?;
            write_bytes(&mut self.writer, &self.raw, &self.source_label)?;
        }
        self.partial = false;
        self.frame_count = self.frame_count.checked_add(1).ok_or_else(|| {
            writer_limit(
                &self.source_label,
                "The frame count exceeds this platform's address range.",
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
                Some(TrajectoryFormat::Trr),
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
    let detail = if mixed {
        "The file's TRR frames mix scalar widths. One writer cannot continue them.".to_owned()
    } else if first.precision != options.precision {
        format!(
            "The file's TRR frames are in {}, the write options in {}. Append with the options the file was written with.",
            first.precision.name(),
            options.precision.name()
        )
    } else {
        return Ok(());
    };
    Err(codec_context(
        TrajectoryIoErrorKind::InconsistentMetadata,
        TrajectoryIoOperation::Open,
        Some(TrajectoryFormat::Trr),
        source_label,
        detail,
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
        "position",
        "The positions exceed this platform's address range.",
    )
}

pub(super) fn encode_vectors_to_raw(
    raw: &mut Vec<u8>,
    vectors: &[Vector3],
    factor: f64,
    precision: TrrScalarPrecision,
    source_label: &str,
    field: &str,
) -> Result<(), TrajectoryIoError> {
    encode_triplets_to_raw(
        raw,
        vectors.iter().map(|v| [v.x, v.y, v.z]),
        factor,
        precision,
        source_label,
        field,
        "The vectors exceed this platform's address range.",
    )
}

pub(super) fn encode_triplets_to_raw(
    raw: &mut Vec<u8>,
    values: impl ExactSizeIterator<Item = [f64; 3]>,
    factor: f64,
    precision: TrrScalarPrecision,
    source_label: &str,
    field: &str,
    overflow_message: &str,
) -> Result<(), TrajectoryIoError> {
    raw.clear();
    let bytes = values
        .len()
        .checked_mul(3 * precision.bytes())
        .ok_or_else(|| writer_limit(source_label, overflow_message))?;
    raw.try_reserve(bytes.saturating_sub(raw.len()))
        .map_err(|_| writer_limit(source_label, "Cannot allocate the writer's buffers."))?;
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
    field: &str,
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
    field: &str,
) -> Result<(), TrajectoryIoError> {
    if !value.is_finite()
        || (precision == TrrScalarPrecision::Float32 && !(value as f32).is_finite())
    {
        return Err(codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Trr),
            source_label,
            format!("A {field} value does not fit the selected precision."),
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
                    Some(TrajectoryFormat::Trr),
                    source_label,
                    "The lambda policy names a property other than gromacs.trr.lambda.",
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
                    Some(TrajectoryFormat::Trr),
                    source_label,
                    "The lambda property is not one finite, dimensionless number.",
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
                source_label,
                "The version string exceeds the 32-bit length field.",
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
            Some(TrajectoryFormat::Trr),
            source_label,
            error,
        )
    })
}
