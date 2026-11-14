use super::*;

/// Sequential TRR reader retaining one seekable stream and reusable scratch.
pub struct TrrReader<R> {
    pub(super) reader: R,
    pub(super) topology: Arc<Topology>,
    pub(super) options: TrrReadOptions,
    pub(super) stream_start: u64,
    pub(super) current_header_offset: u64,
    pub(super) first_header: TrrFrameHeader,
    pub(super) pending_header: Option<TrrFrameHeader>,
    pub(super) positions: Vec<Point3>,
    pub(super) velocities: Vec<Vector3>,
    pub(super) forces: Vec<Vector3>,
    pub(super) raw: Vec<u8>,
    pub(super) properties: OwnerProperties,
    pub(super) frame_cursor: usize,
    pub(super) precision_mixed: bool,
}

impl<R> fmt::Debug for TrrReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("TrrReader")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("frame_cursor", &self.frame_cursor)
            .finish_non_exhaustive()
    }
}

#[derive(Debug, Clone)]
pub(super) struct TrrDecodedFrame {
    pub(super) header: TrrFrameHeader,
    pub(super) cell: Option<PeriodicCell>,
}

impl<R: Read + Seek> TrrReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order.
    ///
    /// Accepts an owned topology or a shared `Arc<Topology>`. Counts and available
    /// format metadata are checked automatically; matching counts alone cannot
    /// establish atom identity. Limits and diagnostic labels are set in `options`.
    /// The first frame header is read and validated here.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidHeader`] for an empty stream or
    /// malformed header, [`TrajectoryIoErrorKind::UnsupportedVariant`] for an
    /// unsupported header layout, [`TrajectoryIoErrorKind::InconsistentAtomCount`]
    /// when its atom count differs from the topology's,
    /// [`TrajectoryIoErrorKind::InconsistentMetadata`] for a nonzero lambda
    /// under [`TrrLambdaPolicy::RequireZero`],
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when the limits are exceeded,
    /// and [`TrajectoryIoErrorKind::Io`] when reading fails.
    pub fn new(
        mut reader: R,
        topology: impl Into<Arc<Topology>>,
        options: TrrReadOptions,
    ) -> Result<Self, TrajectoryIoError> {
        let topology = topology.into();
        let source_label = &options.source_label;
        let limits = &options.limits;
        let atom_count = topology.atom_count();
        validate_atom_count(
            FORMAT,
            atom_count,
            i32::MAX as usize / 3,
            limits,
            source_label,
        )?;
        let stream_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                source_label,
                error,
            )
        })?;
        let first_header = read_header(&mut reader, limits, source_label, Some(atom_count), false)?
            .ok_or_else(|| header_error(FORMAT, source_label, "The stream is empty."))?;
        validate_lambda(first_header.lambda, options.lambda_policy, source_label, 0)?;
        let mut positions = Vec::new();
        let mut velocities = Vec::new();
        let mut forces = Vec::new();
        let vector_bytes = atom_count
            .checked_mul(std::mem::size_of::<Vector3>())
            .ok_or_else(|| {
                limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadHeader,
                    source_label,
                    None,
                    "The vector buffers exceed this platform's address range.",
                )
            })?;
        let scratch_bytes = vector_bytes.checked_mul(3).ok_or_else(|| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                "The vector buffers exceed this platform's address range.",
            )
        })?;
        if scratch_bytes > limits.max_scratch_bytes() {
            return Err(limit_error(FORMAT, TrajectoryIoOperation::ReadHeader, source_label, None, format!(
                    "The vector buffers need {scratch_bytes} bytes, more than the buffer limit of {} bytes.",
                    limits.max_scratch_bytes()
                )));
        }
        positions.try_reserve_exact(atom_count).map_err(|_| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                "Cannot allocate the position buffer.",
            )
        })?;
        velocities.try_reserve_exact(atom_count).map_err(|_| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                "Cannot allocate the velocity buffer.",
            )
        })?;
        forces.try_reserve_exact(atom_count).map_err(|_| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                "Cannot allocate the force buffer.",
            )
        })?;
        let mut properties = OwnerProperties::new();
        if options.lambda_policy == TrrLambdaPolicy::FrameProperty {
            properties
                .insert(
                    trr_lambda_key(),
                    PropertyValue::Real {
                        value: 0.0,
                        unit: DIMENSIONLESS,
                    },
                )
                .expect("TRR lambda scratch property is valid");
        }
        Ok(Self {
            reader,
            topology,
            options,
            stream_start,
            current_header_offset: stream_start,
            first_header: first_header.clone(),
            pending_header: Some(first_header),
            positions,
            velocities,
            forces,
            raw: Vec::new(),
            properties,
            frame_cursor: 0,
            precision_mixed: false,
        })
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.topology.as_ref()
    }

    pub(crate) fn first_header(&self) -> &TrrFrameHeader {
        &self.first_header
    }

    pub(crate) const fn precision_mixed(&self) -> bool {
        self.precision_mixed
    }

    pub(super) fn next_header(&mut self) -> Result<Option<TrrFrameHeader>, TrajectoryIoError> {
        if let Some(header) = self.pending_header.take() {
            return Ok(Some(header));
        }
        self.current_header_offset = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                &self.options.source_label,
                error,
            )
        })?;
        let atom_count = self.topology.as_ref().atom_count();
        read_header(
            &mut self.reader,
            &self.options.limits,
            &self.options.source_label,
            Some(atom_count),
            true,
        )
        .map_err(|error| frame_offset_context(error, self.frame_cursor, self.current_header_offset))
    }

    pub(super) fn parse_next(&mut self) -> Result<Option<TrrDecodedFrame>, TrajectoryIoError> {
        if self.frame_cursor >= self.options.limits.max_frames() {
            let clean_eof = if self.pending_header.is_some() {
                false
            } else {
                probe_seekable_eof(
                    &mut self.reader,
                    TrajectoryIoOperation::ReadFrame,
                    TrajectoryFormat::Trr,
                    &self.options.source_label,
                )?
            };
            if clean_eof {
                return Ok(None);
            }
            return Err(limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.options.source_label,
                Some(self.frame_cursor),
                "",
            )
            .exceeding(
                "The trajectory",
                LimitExceeded::first_beyond(Resource::Frames, self.options.limits.max_frames()),
            ));
        }
        let Some(header) = self.next_header()? else {
            return Ok(None);
        };
        self.precision_mixed |= header.precision != self.first_header.precision;
        validate_lambda(
            header.lambda,
            self.options.lambda_policy,
            &self.options.source_label,
            self.frame_cursor,
        )?;
        let payload_bytes = header
            .box_size
            .checked_add(header.x_size)
            .and_then(|bytes| bytes.checked_add(header.v_size))
            .and_then(|bytes| bytes.checked_add(header.f_size))
            .ok_or_else(|| {
                limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadFrame,
                    &self.options.source_label,
                    Some(self.frame_cursor),
                    "The frame's size exceeds this platform's address range.",
                )
            })?;
        let frame_bytes = header
            .header_bytes
            .checked_add(payload_bytes as u64)
            .ok_or_else(|| {
                limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadFrame,
                    &self.options.source_label,
                    Some(self.frame_cursor),
                    "The frame's size exceeds the 64-bit range.",
                )
            })?;
        if frame_bytes > self.options.limits.max_frame_bytes() {
            return Err(limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.options.source_label,
                Some(self.frame_cursor),
                "",
            )
            .exceeding(
                "The frame",
                exceeded_bytes(self.options.limits.max_frame_bytes(), frame_bytes),
            ));
        }

        // The reusable raw record and all three dense arrays coexist. Account
        // for retained capacity as well as any growth needed by this frame;
        // each dense array holds one vector per atom, whether it is kept or
        // allocated by this read, as the check at opening counts them.
        let atom_count = self.topology.atom_count();
        let raw_bytes = self
            .raw
            .capacity()
            .max(header.box_size)
            .max(header.x_size)
            .max(header.v_size)
            .max(header.f_size);
        let scratch_bytes = filled_bytes(&self.positions, atom_count)
            .and_then(|bytes| bytes.checked_add(filled_bytes(&self.velocities, atom_count)?))
            .and_then(|bytes| bytes.checked_add(filled_bytes(&self.forces, atom_count)?))
            .and_then(|bytes| bytes.checked_add(raw_bytes));
        if scratch_bytes.is_none_or(|bytes| bytes > self.options.limits.max_scratch_bytes()) {
            return Err(limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.options.source_label,
                Some(self.frame_cursor),
                format!(
                    "The frame's buffers exceed the buffer limit of {} bytes.",
                    self.options.limits.max_scratch_bytes()
                ),
            ));
        }

        let cell = if header.has_cell() {
            read_raw(
                &mut self.reader,
                &mut self.raw,
                header.box_size,
                &self.options.limits,
                &self.options.source_label,
                self.frame_cursor,
                "box",
            )?;
            Some(decode_cell(
                &self.raw,
                header.precision,
                &self.options.source_label,
                self.frame_cursor,
            )?)
        } else {
            None
        };
        read_raw(
            &mut self.reader,
            &mut self.raw,
            header.x_size,
            &self.options.limits,
            &self.options.source_label,
            self.frame_cursor,
            "position",
        )?;
        decode_triplets_into(
            &self.raw,
            header.precision,
            &mut self.positions,
            Point3::new,
            &self.options.source_label,
            self.frame_cursor,
            "position",
        )?;
        if header.has_velocities() {
            read_raw(
                &mut self.reader,
                &mut self.raw,
                header.v_size,
                &self.options.limits,
                &self.options.source_label,
                self.frame_cursor,
                "velocity",
            )?;
            decode_triplets_into(
                &self.raw,
                header.precision,
                &mut self.velocities,
                Vector3::new,
                &self.options.source_label,
                self.frame_cursor,
                "velocity",
            )?;
        }
        if header.has_forces() {
            read_raw(
                &mut self.reader,
                &mut self.raw,
                header.f_size,
                &self.options.limits,
                &self.options.source_label,
                self.frame_cursor,
                "force",
            )?;
            decode_triplets_into(
                &self.raw,
                header.precision,
                &mut self.forces,
                Vector3::new,
                &self.options.source_label,
                self.frame_cursor,
                "force",
            )?;
        }
        if self.options.lambda_policy == TrrLambdaPolicy::FrameProperty {
            self.properties
                .insert(
                    trr_lambda_key(),
                    PropertyValue::Real {
                        value: header.lambda,
                        unit: DIMENSIONLESS,
                    },
                )
                .map_err(|_| {
                    frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::InvalidFrame,
                        TrajectoryIoOperation::ReadFrame,
                        &self.options.source_label,
                        self.frame_cursor,
                        "The lambda value is not a valid frame property.",
                    )
                })?;
        }
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.options.source_label,
                Some(self.frame_cursor),
                "The frame number exceeds this platform's address range.",
            )
        })?;
        Ok(Some(TrrDecodedFrame { header, cell }))
    }

    pub(super) fn publish(
        &self,
        positions: &[Point3],
        velocities: &[Vector3],
        forces: &[Vector3],
        properties: &OwnerProperties,
        decoded: &TrrDecodedFrame,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        let mut data = FrameData::new(Quantity::new(positions, NANOMETER))
            .with_time(Quantity::new(decoded.header.time, PICOSECOND))
            .with_step(decoded.header.step)
            .with_owner_properties(properties);
        if let Some(cell) = decoded.cell {
            data = data.with_cell(cell);
        }
        if decoded.header.has_velocities() {
            data = data.with_velocities(Quantity::new(velocities, NANOMETER_PER_PICOSECOND));
        }
        if decoded.header.has_forces() {
            data = data.with_forces(Quantity::new(forces, KILOJOULE_PER_MOLE_PER_NANOMETER));
        }
        destination.replace_from_data(data).map_err(Into::into)
    }

    /// Indexes a freshly opened reader, whose next frame is frame 0.
    pub(super) fn index_frames(
        self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<IndexedTrrReader<R>, TrajectoryIoError> {
        let mut indexed = IndexedTrrReader {
            offsets: Arc::new(FrameOffsets::empty(self.stream_start)),
            precision_mixed: false,
            inner: self,
            pending_bytes: 0,
        };
        indexed.extend(growing, control)?;
        Ok(indexed)
    }

    pub(super) fn position(
        &mut self,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.reader
            .stream_position()
            .map_err(|error| io_context(operation, Some(FORMAT), &self.options.source_label, error))
    }

    pub(super) fn seek_to(
        &mut self,
        offset: SeekFrom,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.reader
            .seek(offset)
            .map_err(|error| io_context(operation, Some(FORMAT), &self.options.source_label, error))
    }
}

impl<R: Read + Seek> TrajectoryReader for TrrReader<R> {
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
        let offset = if self.pending_header.is_some() {
            self.current_header_offset
        } else {
            self.reader.stream_position().map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadFrame,
                    Some(FORMAT),
                    &self.options.source_label,
                    error,
                )
            })?
        };
        let decoded = self
            .parse_next()
            .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))?;
        let Some(decoded) = decoded else {
            return Ok(false);
        };
        self.publish(
            &self.positions,
            &self.velocities,
            &self.forces,
            &self.properties,
            &decoded,
            destination,
        )?;
        Ok(true)
    }
}

pub(super) fn read_header<R: Read>(
    reader: &mut R,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    expected_atoms: Option<usize>,
    clean_eof: bool,
) -> Result<Option<TrrFrameHeader>, TrajectoryIoError> {
    let mut magic = [0_u8; 4];
    match reader.read(&mut magic[..1]) {
        Ok(0) if clean_eof => return Ok(None),
        Ok(0) => {
            return Err(header_error(
                FORMAT,
                source_label,
                "The magic number is missing.",
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
    read_exact(
        FORMAT,
        reader,
        &mut magic[1..],
        TrajectoryIoOperation::ReadHeader,
        source_label,
        None,
        "The magic number is incomplete.",
    )?;
    if i32::from_be_bytes(magic) != TRR_MAGIC {
        return Err(header_error(
            FORMAT,
            source_label,
            "The magic number is not 1993.",
        ));
    }
    let declared_version_bytes = read_i32(
        reader,
        source_label,
        None,
        "The version length is truncated.",
    )?;
    if declared_version_bytes
        != i32::try_from(TRR_VERSION.len() + 1).expect("TRR version length fits i32")
    {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
            source_label,
            "The version length is not that of GMX_trn_file.",
        ));
    }
    let (version_matches, version_bytes) = read_xdr_string(reader, limits, source_label)?;
    if !version_matches {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
            source_label,
            "The version string is not GMX_trn_file.",
        ));
    }
    let mut sizes = [0_usize; 10];
    for size in &mut sizes {
        let value = read_i32(reader, source_label, None, "A block size is truncated.")?;
        *size = usize::try_from(value)
            .map_err(|_| header_error(FORMAT, source_label, "A block size is negative."))?;
        if *size as u64 > limits.max_record_bytes() {
            return Err(limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                "",
            )
            .exceeding(
                "A block",
                exceeded_bytes(limits.max_record_bytes(), *size as u64),
            ));
        }
    }
    let atom_count = usize::try_from(read_i32(
        reader,
        source_label,
        None,
        "The atom count is truncated.",
    )?)
    .map_err(|_| header_error(FORMAT, source_label, "The atom count is negative."))?;
    validate_atom_count(
        FORMAT,
        atom_count,
        i32::MAX as usize / 3,
        limits,
        source_label,
    )?;
    if let Some(expected) = expected_atoms {
        if atom_count != expected {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
            )
            .with_source_label(source_label)
            .with_counts(expected as u64, atom_count as u64));
        }
    }
    let [ir, energy, box_size, virial, pressure, topology, symbols, x_size, v_size, f_size] = sizes;
    if [ir, energy, virial, pressure, topology, symbols]
        .into_iter()
        .any(|size| size != 0)
    {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedField,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
            source_label,
            "The frame has input record, energy, virial, pressure, topology or symbol blocks. Kekule cannot keep them.",
        ));
    }
    if x_size == 0 {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
            source_label,
            "The frame has no positions.",
        ));
    }
    let precision = infer_precision(atom_count, box_size, x_size, v_size, f_size, source_label)?;
    let expected_vectors = atom_count
        .checked_mul(3)
        .and_then(|count| count.checked_mul(precision.bytes()))
        .ok_or_else(|| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                "A vector block exceeds this platform's address range.",
            )
        })?;
    if x_size != expected_vectors
        || (v_size != 0 && v_size != expected_vectors)
        || (f_size != 0 && f_size != expected_vectors)
        || (box_size != 0 && box_size != 9 * precision.bytes())
    {
        return Err(codec_context(
            TrajectoryIoErrorKind::InvalidRecordLength,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
            source_label,
            "The block sizes do not match the atom count and the scalar precision.",
        ));
    }
    let step_raw = read_i32(reader, source_label, None, "The step is truncated.")?;
    let step = u64::try_from(step_raw).map_err(|_| {
        codec_context(
            TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
            source_label,
            "The step is negative.",
        )
    })?;
    let nre = read_i32(
        reader,
        source_label,
        None,
        "The energy count nre is truncated.",
    )?;
    if nre != 0 {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedField,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
            source_label,
            "The frame has energy terms. Kekule cannot keep them.",
        ));
    }
    let time = read_scalar(
        reader,
        precision,
        source_label,
        None,
        "The time is truncated.",
    )?;
    let lambda = read_scalar(
        reader,
        precision,
        source_label,
        None,
        "The lambda value is truncated.",
    )?;
    if !time.is_finite() || !lambda.is_finite() {
        return Err(header_error(
            FORMAT,
            source_label,
            format!("The time {time} or the lambda value {lambda} is not finite."),
        ));
    }
    let header_bytes = 4_u64
        .checked_add(4)
        .and_then(|bytes| bytes.checked_add(version_bytes))
        .and_then(|bytes| bytes.checked_add(11 * 4 + 2 * 4))
        .and_then(|bytes| bytes.checked_add((2 * precision.bytes()) as u64))
        .ok_or_else(|| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                "The frame header's size exceeds the 64-bit range.",
            )
        })?;
    Ok(Some(TrrFrameHeader {
        precision,
        atom_count,
        box_size,
        x_size,
        v_size,
        f_size,
        step,
        time,
        lambda,
        header_bytes,
    }))
}

pub(super) fn infer_precision(
    atom_count: usize,
    box_size: usize,
    x_size: usize,
    v_size: usize,
    f_size: usize,
    source_label: &str,
) -> Result<TrrScalarPrecision, TrajectoryIoError> {
    let candidates = [
        (box_size, 9_usize),
        (x_size, atom_count.saturating_mul(3)),
        (v_size, atom_count.saturating_mul(3)),
        (f_size, atom_count.saturating_mul(3)),
    ];
    let (size, scalars) = candidates
        .into_iter()
        .find(|(size, _)| *size != 0)
        .ok_or_else(|| {
            header_error(
                FORMAT,
                source_label,
                "The frame has no block that shows its precision.",
            )
        })?;
    let width = size
        .checked_div(scalars)
        .filter(|_| size % scalars == 0)
        .ok_or_else(|| {
            header_error(
                FORMAT,
                source_label,
                "The block sizes give no whole scalar width.",
            )
        })?;
    match width {
        4 => Ok(TrrScalarPrecision::Float32),
        8 => Ok(TrrScalarPrecision::Float64),
        _ => Err(codec_context(
            TrajectoryIoErrorKind::InvalidPrecision,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
            source_label,
            format!("The scalar width is {width} bytes, not 4 or 8."),
        )),
    }
}

pub(super) fn read_xdr_string<R: Read>(
    reader: &mut R,
    limits: &TrajectoryIoLimits,
    source_label: &str,
) -> Result<(bool, u64), TrajectoryIoError> {
    let length = usize::try_from(read_i32(
        reader,
        source_label,
        None,
        "The version string's length is truncated.",
    )?)
    .map_err(|_| {
        header_error(
            FORMAT,
            source_label,
            "The version string's length is negative.",
        )
    })?;
    if length > 256 || length > limits.max_text_line_bytes() || length > limits.max_scratch_bytes()
    {
        return Err(limit_error(
            FORMAT,
            TrajectoryIoOperation::ReadHeader,
            source_label,
            None,
            format!("The version string has {length} bytes, more than 256 or the line limit."),
        ));
    }
    let padded = length
        .checked_add(3)
        .map(|value| value & !3)
        .ok_or_else(|| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadHeader,
                source_label,
                None,
                "The version string's padding exceeds this platform's address range.",
            )
        })?;
    let mut bytes = [0_u8; 256];
    read_exact(
        FORMAT,
        reader,
        &mut bytes[..padded],
        TrajectoryIoOperation::ReadHeader,
        source_label,
        None,
        "The version string is truncated.",
    )?;
    if bytes[length..].iter().any(|byte| *byte != 0) {
        return Err(header_error(
            FORMAT,
            source_label,
            "The version string's padding is not zero.",
        ));
    }
    Ok((&bytes[..length] == TRR_VERSION, 4 + padded as u64))
}

pub(super) fn read_raw<R: Read>(
    reader: &mut R,
    raw: &mut Vec<u8>,
    length: usize,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    frame: usize,
    field: &str,
) -> Result<(), TrajectoryIoError> {
    if length > limits.max_scratch_bytes() || length as u64 > limits.max_record_bytes() {
        return Err(limit_error(FORMAT, TrajectoryIoOperation::ReadFrame, source_label, Some(frame), format!(
                "The {field} block has {length} bytes, more than the record limit of {} bytes or the buffer limit of {} bytes.",
                limits.max_record_bytes(),
                limits.max_scratch_bytes()
            )));
    }
    if raw.capacity() < length {
        raw.try_reserve_exact(length.saturating_sub(raw.len()))
            .map_err(|_| {
                limit_error(
                    FORMAT,
                    TrajectoryIoOperation::ReadFrame,
                    source_label,
                    Some(frame),
                    format!("Cannot allocate the {field} buffer."),
                )
            })?;
    }
    raw.resize(length, 0);
    read_exact(
        FORMAT,
        reader,
        raw,
        TrajectoryIoOperation::ReadFrame,
        source_label,
        Some(frame),
        &format!("truncated {field}"),
    )
}

pub(super) fn decode_triplets<T>(
    raw: &[u8],
    precision: TrrScalarPrecision,
    destination: &mut [T],
    construct: impl Fn(f64, f64, f64) -> T,
    source_label: &str,
    frame: usize,
    field: &str,
) -> Result<(), TrajectoryIoError> {
    for (destination, values) in destination
        .iter_mut()
        .zip(raw.chunks_exact(3 * precision.bytes()))
    {
        let x = scalar_at(values, 0, precision);
        let y = scalar_at(values, 1, precision);
        let z = scalar_at(values, 2, precision);
        if ![x, y, z].into_iter().all(f64::is_finite) {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame,
                format!("A {field} value is not finite."),
            ));
        }
        *destination = construct(x, y, z);
    }
    Ok(())
}

/// As [`decode_triplets`], writing each triplet once into `destination`
/// in place of what it held, so the buffer can be handed over.
pub(super) fn decode_triplets_into<T>(
    raw: &[u8],
    precision: TrrScalarPrecision,
    destination: &mut Vec<T>,
    construct: impl Fn(f64, f64, f64) -> T,
    source_label: &str,
    frame: usize,
    field: &str,
) -> Result<(), TrajectoryIoError> {
    destination.clear();
    // The indexed reader hands the buffers over, so the next read
    // allocates its own on the thread that reads.
    destination.reserve_exact(raw.len() / (3 * precision.bytes()));
    for values in raw.chunks_exact(3 * precision.bytes()) {
        let x = scalar_at(values, 0, precision);
        let y = scalar_at(values, 1, precision);
        let z = scalar_at(values, 2, precision);
        if ![x, y, z].into_iter().all(f64::is_finite) {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame,
                format!("A {field} value is not finite."),
            ));
        }
        destination.push(construct(x, y, z));
    }
    Ok(())
}

pub(super) fn decode_cell(
    raw: &[u8],
    precision: TrrScalarPrecision,
    source_label: &str,
    frame: usize,
) -> Result<PeriodicCell, TrajectoryIoError> {
    let mut vectors = [Vector3::zero(); 3];
    decode_triplets(
        raw,
        precision,
        &mut vectors,
        Vector3::new,
        source_label,
        frame,
        "box",
    )?;
    PeriodicCell::new(Quantity::new(vectors, NANOMETER), [true; 3]).map_err(|error| {
        frame_error(
            FORMAT,
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame,
            "The periodic cell is invalid.",
        )
        .with_source(error)
    })
}

pub(super) fn scalar_at(bytes: &[u8], index: usize, precision: TrrScalarPrecision) -> f64 {
    let start = index * precision.bytes();
    match precision {
        TrrScalarPrecision::Float32 => f64::from(f32::from_be_bytes(
            bytes[start..start + 4].try_into().expect("f32"),
        )),
        TrrScalarPrecision::Float64 => {
            f64::from_be_bytes(bytes[start..start + 8].try_into().expect("f64"))
        }
    }
}

pub(super) fn validate_lambda(
    lambda: f64,
    policy: TrrLambdaPolicy,
    source_label: &str,
    frame: usize,
) -> Result<(), TrajectoryIoError> {
    if policy == TrrLambdaPolicy::RequireZero && lambda != 0.0 {
        return Err(TrajectoryIoError::codec(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::ReadFrame,
            Some(FORMAT),
        )
        .with_source_label(source_label)
        .with_frame(frame)
        .with_detail("The lambda value is not zero."));
    }
    Ok(())
}

pub(super) fn read_i32<R: Read>(
    reader: &mut R,
    source_label: &str,
    frame: Option<usize>,
    detail: &str,
) -> Result<i32, TrajectoryIoError> {
    let mut bytes = [0_u8; 4];
    read_exact(
        FORMAT,
        reader,
        &mut bytes,
        read_operation(frame),
        source_label,
        frame,
        detail,
    )?;
    Ok(i32::from_be_bytes(bytes))
}

pub(super) fn read_scalar<R: Read>(
    reader: &mut R,
    precision: TrrScalarPrecision,
    source_label: &str,
    frame: Option<usize>,
    detail: &str,
) -> Result<f64, TrajectoryIoError> {
    let mut bytes = [0_u8; 8];
    read_exact(
        FORMAT,
        reader,
        &mut bytes[..precision.bytes()],
        read_operation(frame),
        source_label,
        frame,
        detail,
    )?;
    Ok(match precision {
        TrrScalarPrecision::Float32 => {
            f64::from(f32::from_be_bytes(bytes[..4].try_into().expect("f32")))
        }
        TrrScalarPrecision::Float64 => f64::from_be_bytes(bytes),
    })
}
