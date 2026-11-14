use super::*;

/// Sequential XYZ reader over a caller-supplied buffered stream.
pub struct XyzReader<R> {
    pub(super) reader: R,
    pub(super) topology: Arc<Topology>,
    pub(super) options: XyzReadOptions,
    pub(super) frame_cursor: usize,
    pub(super) line: Vec<u8>,
    pub(super) positions: Vec<Point3>,
    /// Whether a line counts only once its line break is written, as in a
    /// file still being written.
    pub(super) growing: bool,
}

impl<R> fmt::Debug for XyzReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("XyzReader")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("frame_cursor", &self.frame_cursor)
            .finish_non_exhaustive()
    }
}

impl<R: BufRead> XyzReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order.
    ///
    /// Accepts an owned topology or a shared `Arc<Topology>`. Counts and available
    /// format metadata are checked automatically; matching counts alone cannot
    /// establish atom identity. Limits and diagnostic labels are set in `options`.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the length
    /// unit is not a length, and [`TrajectoryIoErrorKind::LimitExceeded`] when
    /// the atom count exceeds the limits or position scratch cannot be
    /// reserved.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: XyzReadOptions,
    ) -> Result<Self, TrajectoryIoError> {
        let topology = topology.into();
        let source_label = &options.source_label;
        let limits = &options.limits;
        options
            .length_unit
            .conversion_factor_to(CANONICAL_LENGTH_UNIT)
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Open,
                    Some(TrajectoryFormat::Xyz),
                    source_label,
                    "The XYZ length unit is not a length.",
                )
                .with_source(error)
            })?;
        let atom_count = topology.atom_count();
        validate_atom_limit(atom_count, limits, source_label)?;
        let mut positions = Vec::new();
        positions.try_reserve_exact(atom_count).map_err(|_| {
            codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                source_label,
                "Cannot allocate the position buffer.",
            )
        })?;
        Ok(Self {
            reader,
            topology,
            options,
            frame_cursor: 0,
            line: Vec::new(),
            positions,
            growing: false,
        })
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.topology.as_ref()
    }

    pub(super) fn parse_next(&mut self, retain_positions: bool) -> Result<bool, TrajectoryIoError> {
        let operation = if self.frame_cursor == 0 {
            TrajectoryIoOperation::ReadHeader
        } else {
            TrajectoryIoOperation::ReadFrame
        };
        if self.frame_cursor >= self.options.limits.max_frames() {
            if self
                .reader
                .fill_buf()
                .map_err(|error| {
                    io_context(
                        operation,
                        Some(TrajectoryFormat::Xyz),
                        &self.options.source_label,
                        error,
                    )
                })?
                .is_empty()
            {
                return Ok(false);
            }
            return Err(frame_codec_error(
                TrajectoryIoErrorKind::LimitExceeded,
                operation,
                &self.options.source_label,
                self.frame_cursor,
                "",
            )
            .exceeding(
                "The trajectory",
                LimitExceeded::first_beyond(Resource::Frames, self.options.limits.max_frames()),
            ));
        }
        // The position array holds one point per atom, whether it is kept or
        // allocated by this read; the line buffer may grow only within what
        // it leaves of the scratch budget.
        let position_bytes = filled_bytes(&self.positions, self.topology().atom_count());
        let growing = self.growing;
        let Some((count_line, consumed)) = read_line(
            &mut self.reader,
            &mut self.line,
            self.options.limits.max_text_line_bytes(),
            position_bytes,
            &self.options.limits,
            &self.options.source_label,
            self.frame_cursor,
            operation,
            growing,
        )?
        else {
            return Ok(false);
        };
        let mut frame_bytes = checked_frame_add(
            0,
            consumed,
            &self.options.source_label,
            self.frame_cursor,
            operation,
        )?;
        validate_frame_bytes(
            frame_bytes,
            &self.options.limits,
            &self.options.source_label,
            self.frame_cursor,
            operation,
        )?;
        let atom_count = count_line.trim().parse::<usize>().map_err(|_| {
            frame_codec_error(
                TrajectoryIoErrorKind::InvalidFrame,
                operation,
                &self.options.source_label,
                self.frame_cursor,
                "The atom-count line is not one nonnegative whole number.",
            )
        })?;
        if atom_count == 0 {
            return Err(frame_codec_error(
                TrajectoryIoErrorKind::InvalidFrame,
                operation,
                &self.options.source_label,
                self.frame_cursor,
                "The frame has zero atoms.",
            ));
        }
        validate_atom_limit(atom_count, &self.options.limits, &self.options.source_label)?;
        if atom_count != self.topology().atom_count() {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentAtomCount,
                operation,
                Some(TrajectoryFormat::Xyz),
            )
            .with_source_label(&self.options.source_label)
            .with_frame(self.frame_cursor)
            .with_counts(self.topology().atom_count() as u64, atom_count as u64));
        }
        let Some((_, consumed)) = read_line(
            &mut self.reader,
            &mut self.line,
            self.options.limits.max_comment_bytes(),
            position_bytes,
            &self.options.limits,
            &self.options.source_label,
            self.frame_cursor,
            operation,
            growing,
        )?
        else {
            return Err(frame_codec_error(
                TrajectoryIoErrorKind::TruncatedRecord,
                operation,
                &self.options.source_label,
                self.frame_cursor,
                "The frame has no comment line.",
            ));
        };
        frame_bytes = checked_frame_add(
            frame_bytes,
            consumed,
            &self.options.source_label,
            self.frame_cursor,
            operation,
        )?;
        validate_frame_bytes(
            frame_bytes,
            &self.options.limits,
            &self.options.source_label,
            self.frame_cursor,
            operation,
        )?;

        // The indexed reader hands the buffer over, so the next read
        // allocates its own on the thread that reads.
        self.positions.clear();
        self.positions.reserve_exact(atom_count);
        for atom_index in 0..atom_count {
            let atom_id = self.topology().atom_ids()[atom_index];
            let expected_element = self
                .topology()
                .atom(atom_id)
                .expect("dense topology atom IDs name topology atoms")
                .element();
            let Some((line, consumed)) = read_line(
                &mut self.reader,
                &mut self.line,
                self.options.limits.max_text_line_bytes(),
                position_bytes,
                &self.options.limits,
                &self.options.source_label,
                self.frame_cursor,
                operation,
                growing,
            )?
            else {
                return Err(frame_codec_error(
                    TrajectoryIoErrorKind::TruncatedRecord,
                    operation,
                    &self.options.source_label,
                    self.frame_cursor,
                    format!("The frame ends before atom line {atom_index}."),
                ));
            };
            frame_bytes = checked_frame_add(
                frame_bytes,
                consumed,
                &self.options.source_label,
                self.frame_cursor,
                operation,
            )?;
            validate_frame_bytes(
                frame_bytes,
                &self.options.limits,
                &self.options.source_label,
                self.frame_cursor,
                operation,
            )?;
            let mut fields = line.split_whitespace();
            let (Some(symbol), Some(x), Some(y), Some(z), None) = (
                fields.next(),
                fields.next(),
                fields.next(),
                fields.next(),
                fields.next(),
            ) else {
                return Err(frame_codec_error(
                    TrajectoryIoErrorKind::InvalidFrame,
                    operation,
                    &self.options.source_label,
                    self.frame_cursor,
                    format!("Atom line {atom_index} does not have four fields."),
                ));
            };
            let element = Element::from_symbol(symbol).ok_or_else(|| {
                frame_codec_error(
                    TrajectoryIoErrorKind::InvalidFrame,
                    operation,
                    &self.options.source_label,
                    self.frame_cursor,
                    format!("Atom line {atom_index} has the unknown element symbol {symbol:?}."),
                )
            })?;
            if element != expected_element {
                return Err(frame_codec_error(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    operation,
                    &self.options.source_label,
                    self.frame_cursor,
                    format!(
                        "Atom {atom_index} is {}. The topology has {} there.",
                        element.symbol(),
                        expected_element.symbol()
                    ),
                ));
            }
            let parse_coordinate = |value: &str| {
                value.parse::<f64>().map_err(|_| {
                    frame_codec_error(
                        TrajectoryIoErrorKind::InvalidFrame,
                        operation,
                        &self.options.source_label,
                        self.frame_cursor,
                        format!("Atom line {atom_index} has a coordinate that is not a number."),
                    )
                })
            };
            let x = parse_coordinate(x)?;
            let y = parse_coordinate(y)?;
            let z = parse_coordinate(z)?;
            if ![x, y, z].into_iter().all(f64::is_finite) {
                return Err(frame_codec_error(
                    TrajectoryIoErrorKind::InvalidFrame,
                    operation,
                    &self.options.source_label,
                    self.frame_cursor,
                    format!("Atom line {atom_index} has a coordinate that is not finite."),
                ));
            }
            if retain_positions {
                self.positions.push(Point3::new(x, y, z));
            }
        }
        Ok(true)
    }

    pub(super) fn publish(&self, destination: &mut FrameBuffer) -> Result<(), TrajectoryIoError> {
        destination
            .replace_from_data(FrameData::new(Quantity::new(
                self.positions.as_slice(),
                self.options.length_unit,
            )))
            .map_err(Into::into)
    }
}

impl<R: BufRead> TrajectoryReader for XyzReader<R> {
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
        if !self.parse_next(true)? {
            return Ok(false);
        }
        self.publish(destination)?;
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            frame_codec_error(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::ReadFrame,
                &self.options.source_label,
                self.frame_cursor,
                "The frame number exceeds this platform's address range.",
            )
        })?;
        Ok(true)
    }
}

impl<R: BufRead + Seek> XyzReader<R> {
    /// Parses and validates the next frame without publishing it, then
    /// restores the stream position and frame cursor.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidHeader`] when the stream holds
    /// no frame, the decoding or [`TrajectoryIoErrorKind::Io`] error of an
    /// invalid frame, and an [`TrajectoryIoErrorKind::Io`] error when the
    /// position cannot be restored.
    pub fn validate_first_frame(&mut self) -> Result<(), TrajectoryIoError> {
        let start = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xyz),
                &self.options.source_label,
                error,
            )
        })?;
        let result = self.parse_next(false);
        self.reader.seek(SeekFrom::Start(start)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xyz),
                &self.options.source_label,
                error,
            )
        })?;
        self.frame_cursor = 0;
        if !result? {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidHeader,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xyz),
                &self.options.source_label,
                "The file holds no frame.",
            ));
        }
        Ok(())
    }

    /// Indexes a freshly opened reader, whose next frame is frame 0. A
    /// `growing` stream may end in a partial frame and need not hold a
    /// complete one yet (see [`IndexedXyzReader::extend`]).
    pub(super) fn index_frames(
        self,
        stream_start: u64,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<IndexedXyzReader<R>, TrajectoryIoError> {
        let mut indexed = IndexedXyzReader {
            reader: self,
            offsets: Arc::new(FrameOffsets::empty(stream_start)),
            pending_bytes: 0,
        };
        indexed.extend(growing, control)?;
        if indexed.offsets.is_empty() && !growing {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidHeader,
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Xyz),
                &indexed.reader.options.source_label,
                "The file holds no frame.",
            ));
        }
        Ok(indexed)
    }

    pub(super) fn position(
        &mut self,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.reader.stream_position().map_err(|error| {
            io_context(
                operation,
                Some(TrajectoryFormat::Xyz),
                &self.options.source_label,
                error,
            )
        })
    }

    pub(super) fn seek_to(
        &mut self,
        offset: SeekFrom,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.reader.seek(offset).map_err(|error| {
            io_context(
                operation,
                Some(TrajectoryFormat::Xyz),
                &self.options.source_label,
                error,
            )
        })
    }
}

pub(super) fn validate_atom_limit(
    atom_count: usize,
    limits: &TrajectoryIoLimits,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    if atom_count > limits.max_atoms() {
        return Err(codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Xyz),
            source_label,
            "",
        )
        .exceeding(
            "The atom count",
            LimitExceeded::new(Resource::Atoms, limits.max_atoms(), atom_count),
        ));
    }
    let scratch_bytes = atom_count
        .checked_mul(std::mem::size_of::<Point3>())
        .ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xyz),
                source_label,
                "The position buffer exceeds this platform's address range.",
            )
        })?;
    if scratch_bytes > limits.max_scratch_bytes() {
        return Err(codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Xyz),
            source_label,
            format!(
                "The position buffer needs {scratch_bytes} bytes, more than the buffer limit of {} bytes.",
                limits.max_scratch_bytes()
            ),
        )
        .with_limit(LimitExceeded::new(
            Resource::RecordBytes,
            limits.max_scratch_bytes(),
            scratch_bytes,
        )));
    }
    Ok(())
}

pub(super) fn validate_frame_bytes(
    frame_bytes: u64,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    frame: usize,
    operation: TrajectoryIoOperation,
) -> Result<(), TrajectoryIoError> {
    if frame_bytes > limits.max_frame_bytes() {
        return Err(frame_codec_error(
            TrajectoryIoErrorKind::LimitExceeded,
            operation,
            source_label,
            frame,
            "",
        )
        .exceeding(
            "The frame",
            exceeded_bytes(limits.max_frame_bytes(), frame_bytes),
        ));
    }
    Ok(())
}

/// Adds the bytes one line consumed from the source, terminator included.
pub(super) fn checked_frame_add(
    frame_bytes: u64,
    consumed_bytes: usize,
    source_label: &str,
    frame: usize,
    operation: TrajectoryIoOperation,
) -> Result<u64, TrajectoryIoError> {
    u64::try_from(consumed_bytes)
        .ok()
        .and_then(|bytes| frame_bytes.checked_add(bytes))
        .ok_or_else(|| {
            frame_codec_error(
                TrajectoryIoErrorKind::LimitExceeded,
                operation,
                source_label,
                frame,
                "The frame's size exceeds the 64-bit range.",
            )
        })
}

/// Reads one line and returns it without its LF or CRLF terminator, together
/// with the exact number of source bytes consumed, terminator included.
///
/// The line, terminator included, must fit `limit`. The retained `line` buffer
/// grows only while its capacity plus `other_scratch_bytes` stays within
/// `max_scratch_bytes`, so whether a line fits depends only on its length.
/// With `terminated`, a line the stream ends inside, without its line break,
/// is truncated: its writer may not have finished it.
#[allow(clippy::too_many_arguments)]
pub(super) fn read_line<'a, R: BufRead>(
    reader: &mut R,
    line: &'a mut Vec<u8>,
    limit: usize,
    other_scratch_bytes: Option<usize>,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    frame: usize,
    operation: TrajectoryIoOperation,
    terminated: bool,
) -> Result<Option<(&'a str, usize)>, TrajectoryIoError> {
    line.clear();
    loop {
        let available = match reader.fill_buf() {
            Ok(available) => available,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => {
                return Err(io_context(
                    operation,
                    Some(TrajectoryFormat::Xyz),
                    source_label,
                    error,
                ))
            }
        };
        if available.is_empty() {
            break;
        }
        // Scan at most one byte past the limit, enough to detect an overlong line.
        let window = &available[..available.len().min((limit - line.len()).saturating_add(1))];
        let (chunk, complete) = match window.iter().position(|byte| *byte == b'\n') {
            Some(end) => (&window[..=end], true),
            None => (window, false),
        };
        let needed = line.len() + chunk.len();
        if needed > limit {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                operation,
                Some(TrajectoryFormat::Xyz),
                source_label,
                "",
            )
            .exceeding(
                "A line",
                LimitExceeded::new(Resource::RecordBytes, limit, needed),
            ));
        }
        if needed > line.capacity() {
            // Grow geometrically so a long line is not copied once per buffered
            // chunk, but never past the line limit or the remaining budget.
            let remaining = other_scratch_bytes
                .and_then(|bytes| limits.max_scratch_bytes().checked_sub(bytes))
                .unwrap_or(0);
            let target = line
                .capacity()
                .saturating_mul(2)
                .min(limit)
                .min(remaining)
                .max(needed);
            reserve_scratch(
                line,
                target,
                other_scratch_bytes,
                limits.max_scratch_bytes(),
                |detail| {
                    frame_codec_error(
                        TrajectoryIoErrorKind::LimitExceeded,
                        operation,
                        source_label,
                        frame,
                        format!("The line buffer cannot grow. {detail}"),
                    )
                },
            )?;
        }
        line.extend_from_slice(chunk);
        let consumed = chunk.len();
        reader.consume(consumed);
        if complete {
            break;
        }
    }
    let read = line.len();
    if read == 0 {
        return Ok(None);
    }
    if terminated && !line.ends_with(b"\n") {
        return Err(frame_codec_error(
            TrajectoryIoErrorKind::TruncatedRecord,
            operation,
            source_label,
            frame,
            "The last line has no line break yet.",
        ));
    }
    let mut content = line.as_slice();
    if let Some(stripped) = content.strip_suffix(b"\n") {
        content = stripped.strip_suffix(b"\r").unwrap_or(stripped);
    }
    let content = std::str::from_utf8(content).map_err(|_| {
        io_context(
            operation,
            Some(TrajectoryFormat::Xyz),
            source_label,
            io::Error::new(io::ErrorKind::InvalidData, "The stream is not valid UTF-8."),
        )
    })?;
    Ok(Some((content, read)))
}

pub(super) fn frame_codec_error(
    kind: TrajectoryIoErrorKind,
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: usize,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    TrajectoryIoError::codec(kind, operation, Some(TrajectoryFormat::Xyz))
        .with_source_label(source_label)
        .with_frame(frame)
        .with_detail(detail)
}
