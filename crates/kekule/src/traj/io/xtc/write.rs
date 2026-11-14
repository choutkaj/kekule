use super::*;

/// Strict XTC writer using explicit lossy precision.
pub struct XtcWriter<W> {
    pub(super) writer: W,
    pub(super) topology: Arc<Topology>,
    pub(super) options: XtcWriteOptions,
    pub(super) source_label: String,
    pub(super) frame_count: usize,
    /// Set while a frame's bytes are being written and left set if writing
    /// fails, so a stream that may end in a partial frame is never extended
    /// or reported as finished.
    pub(super) partial: bool,
    /// Nanometre coordinates of the frame being written.
    pub(super) positions: Vec<f32>,
    /// Scaled integer coordinates of a compressed frame.
    pub(super) scaled: Vec<[i32; 3]>,
    /// Compressed bit stream of a frame.
    pub(super) stream: Vec<u8>,
    /// Complete encoded frame, written to the stream in one call.
    pub(super) bytes: Vec<u8>,
    pub(super) report: WriteReport,
}

impl<W> fmt::Debug for XtcWriter<W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("XtcWriter")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("source_label", &self.source_label)
            .field("frame_count", &self.frame_count)
            .field("partial", &self.partial)
            .finish_non_exhaustive()
    }
}

impl<W: Write> XtcWriter<W> {
    /// A writer of `topology` frames to `writer`; `source_label` names the
    /// destination in errors. Nothing is written until the first frame.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidPrecision`] for a precision that
    /// is not finite and positive, and [`TrajectoryIoErrorKind::LimitExceeded`]
    /// when the topology has no atoms or more than the magic's byte count
    /// allows.
    pub fn new(
        writer: W,
        topology: impl Into<Arc<Topology>>,
        options: XtcWriteOptions,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryIoError> {
        let topology = topology.into();
        let source_label = source_label.into();
        validate_precision(options.precision, &source_label)?;
        if topology.atom_count() == 0 || topology.atom_count() > i32::MAX as usize {
            return Err(resource_error(
                &source_label,
                None,
                format!(
                    "The topology has {} atoms. XTC holds 1 to 2,147,483,647.",
                    topology.atom_count()
                ),
            ));
        }
        if options.magic == XtcMagic::Xtc1995 && topology.atom_count() > XTC_1995_MAX_ATOMS {
            return Err(resource_error(
                &source_label,
                None,
                format!(
                    "The topology has {} atoms. XTC 1995 holds at most 298,261,617. Use XTC 2023.",
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
            positions: Vec::new(),
            scaled: Vec::new(),
            stream: Vec::new(),
            bytes: Vec::new(),
            report: WriteReport::empty(),
        })
    }

    /// The frame state the frames written so far leave out under
    /// [`XtcWriteOptions::geometry_metadata`] and
    /// [`XtcWriteOptions::properties`], counted over every frame written. XTC stores coordinates for a topology its readers receive
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
            TrajectoryFormat::Xtc,
            &self.source_label,
        )?;
        require_nonempty_writer(self.frame_count, TrajectoryFormat::Xtc, &self.source_label)
    }

    /// Flushes and returns the completed nonempty XTC stream.
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
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                error,
            )
        })?;
        Ok(self.writer)
    }
}

impl<W: Write> TrajectoryWriter for XtcWriter<W> {
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
            TrajectoryFormat::Xtc,
            &self.source_label,
        )?;
        if !frame.topology().shares_layout(self.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        let omitted = check_frame_state(
            TrajectoryFormat::Xtc,
            &self.source_label,
            self.frame_count,
            frame,
            self.options.state_policies(),
        )?;
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
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "The frame has no step. XTC frames need one.",
            )
        })?;
        let step = i32::try_from(step).map_err(|_| {
            codec_context(
                TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "The frame's step exceeds the 32-bit step field.",
            )
        })?;
        let time = time.ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "The frame has no time. XTC frames need one.",
            )
        })?;
        let time = time.value_in(PICOSECOND).map_err(|error| {
            codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "The XTC time unit is incompatible.",
            )
            .with_source(error)
        })?;
        let time = finite_f32(time, &self.source_label, "time")?;
        let mut boxvec = [0.0_f32; 9];
        match (frame.cell().copied(), self.options.cell_policy) {
            (None, XtcCellPolicy::ZeroMatrixAsAbsent) => {}
            (None, XtcCellPolicy::RequirePeriodic) => {
                return Err(codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Xtc),
                    &self.source_label,
                    "The frame has no periodic cell. XTC frames need one.",
                ))
            }
            (Some(cell), _) => {
                if cell.periodic_axes() != [true; 3] {
                    return Err(codec_context(
                        TrajectoryIoErrorKind::UnsupportedVariant,
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Xtc),
                        &self.source_label,
                        "The cell is not periodic along all three axes. XTC stores only fully periodic cells.",
                    ));
                }
                let vectors = cell.vectors().value_in(NANOMETER).map_err(|error| {
                    codec_context(
                        TrajectoryIoErrorKind::InconsistentMetadata,
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Xtc),
                        &self.source_label,
                        "The XTC cell unit is incompatible.",
                    )
                    .with_source(error)
                })?;
                for (target, value) in boxvec.iter_mut().zip(
                    vectors
                        .into_iter()
                        .flat_map(|vector| [vector.x, vector.y, vector.z]),
                ) {
                    *target = finite_f32(value, &self.source_label, "box")?;
                }
                // Scalar finiteness alone cannot ensure the narrowed matrix still spans
                // three dimensions. Validate the representation the reader will receive.
                decode_cell(
                    boxvec.map(f64::from),
                    XtcCellPolicy::RequirePeriodic,
                    &self.source_label,
                    self.frame_count,
                )
                .map_err(|error| {
                    codec_context(
                        TrajectoryIoErrorKind::InvalidFrame,
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Xtc),
                        &self.source_label,
                        "The cell does not fit single precision.",
                    )
                    .with_source(validation_cause(error))
                })?;
            }
        }
        let positions = frame.positions().values();
        let factor = positions
            .unit()
            .conversion_factor_to(NANOMETER)
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Xtc),
                    &self.source_label,
                    "The XTC position unit is incompatible.",
                )
                .with_source(error)
            })?;
        let atom_count = positions.value().len();
        // Compressed frames store coordinates as integers relative to each
        // axis's minimum, so both the scaled values and their per-axis spans
        // must stay within GROMACS's range.
        let compressed = self.topology.atom_count() > 9;
        let reserve_error = |_| {
            resource_error(
                &self.source_label,
                Some(self.frame_count),
                "Cannot allocate the coordinate buffer.",
            )
        };
        self.positions.clear();
        self.scaled.clear();
        self.positions
            .try_reserve(atom_count.checked_mul(3).ok_or_else(|| {
                resource_error(
                    &self.source_label,
                    Some(self.frame_count),
                    "The coordinate count exceeds this platform's address range.",
                )
            })?)
            .map_err(reserve_error)?;
        if compressed {
            self.scaled.try_reserve(atom_count).map_err(reserve_error)?;
        }
        let mut minimum = [i64::MAX; 3];
        let mut maximum = [i64::MIN; 3];
        for point in *positions.value() {
            for (axis, value) in [point.x * factor, point.y * factor, point.z * factor]
                .into_iter()
                .enumerate()
            {
                let value = finite_f32(value, &self.source_label, "position")?;
                if compressed {
                    // GROMACS scales in f32 and rounds half away from zero.
                    let scaled = (value * self.options.precision).round();
                    if !scaled.is_finite() || f64::from(scaled).abs() > XTC_MAX_SCALED as f64 {
                        return Err(codec_context(
                            TrajectoryIoErrorKind::InvalidFrame,
                            TrajectoryIoOperation::WriteFrame,
                            Some(TrajectoryFormat::Xtc),
                            &self.source_label,
                            "A coordinate exceeds the integer range of the selected precision.",
                        ));
                    }
                    let scaled = scaled as i32;
                    minimum[axis] = minimum[axis].min(i64::from(scaled));
                    maximum[axis] = maximum[axis].max(i64::from(scaled));
                    if axis == 0 {
                        self.scaled.push([scaled, 0, 0]);
                    } else if let Some(last) = self.scaled.last_mut() {
                        last[axis] = scaled;
                    }
                }
                self.positions.push(value);
            }
        }
        if compressed && (0..3).any(|axis| maximum[axis] - minimum[axis] >= XTC_MAX_SCALED) {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "The coordinates span more along one axis than the integer range of the selected precision holds.",
            ));
        }
        let atom_count = i32::try_from(atom_count).map_err(|_| {
            resource_error(
                &self.source_label,
                Some(self.frame_count),
                "The atom count exceeds the 32-bit atom-count field.",
            )
        })?;
        self.bytes.clear();
        for value in [self.options.magic.number(), atom_count, step] {
            self.bytes.extend_from_slice(&value.to_be_bytes());
        }
        for value in std::iter::once(time).chain(boxvec) {
            self.bytes.extend_from_slice(&value.to_be_bytes());
        }
        self.bytes.extend_from_slice(&atom_count.to_be_bytes());
        if compressed {
            self.bytes
                .extend_from_slice(&self.options.precision.to_be_bytes());
            encode::encode_coordinates(
                &mut self.scaled,
                self.options.magic,
                &mut self.stream,
                &mut self.bytes,
            )
            .map_err(|detail| resource_error(&self.source_label, Some(self.frame_count), detail))?;
        } else {
            for value in &self.positions {
                self.bytes.extend_from_slice(&value.to_be_bytes());
            }
        }
        // One call writes the whole frame, but a failing writer may still
        // have taken part of it.
        self.partial = true;
        self.writer.write_all(&self.bytes).map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                error,
            )
        })?;
        self.partial = false;
        self.frame_count = self.frame_count.checked_add(1).ok_or_else(|| {
            resource_error(
                &self.source_label,
                None,
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
            TrajectoryFormat::Xtc,
            &self.source_label,
        )?;
        self.writer.flush().map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                error,
            )
        })
    }
}

/// Fails unless a writer with `options` writes frames as the XTC frames whose
/// first `first` describes were written: with the same magic and, for
/// compressed frames, the same precision.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] naming the setting
/// that differs.
pub(crate) fn check_append_frames(
    first: &XtcFrameInfo,
    options: &XtcWriteOptions,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    let differs = |detail: String| {
        codec_context(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::Open,
            Some(TrajectoryFormat::Xtc),
            source_label,
            format!("{detail} Append with the options the file was written with."),
        )
    };
    if first.magic != options.magic {
        return Err(differs(format!(
            "The file's XTC frames use magic {} and the write options {}.",
            first.magic.number(),
            options.magic.number()
        )));
    }
    match first.precision {
        Some(precision) if precision != options.precision => Err(differs(format!(
            "The file's XTC frames have precision {precision} per nm and the write options {} per nm.",
            options.precision
        ))),
        _ => Ok(()),
    }
}
