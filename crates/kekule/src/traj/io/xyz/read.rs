use super::*;
use crate::traj::trajectory::reason::{Part, Reason};

/// Sequential XYZ reader over a caller-supplied buffered stream.
pub struct XyzReader<R> {
    pub(super) reader: R,
    pub(super) topology: Arc<Topology>,
    pub(super) options: XyzReadOptions,
    pub(super) source: ReadSource,
    pub(super) frame_cursor: usize,
    pub(super) line: Vec<u8>,
    pub(super) positions: Vec<Point3>,
    /// Whether a line counts only once its line break is written, as in a
    /// file still being written.
    pub(super) growing: bool,
    /// Whether the last frame read has columns after the coordinates.
    pub(super) extra_columns: bool,
    /// The frames with columns after the coordinates that the scans of an
    /// indexed reader indexed.
    pub(super) index_extra_columns: usize,
}

impl<R> fmt::Debug for XyzReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("XyzReader")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("source_label", &self.source.label)
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
    /// Returns [`TrajectoryIoErrorKind::LimitExceeded`] when the atom count
    /// exceeds the limits or position scratch cannot be reserved.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: impl Into<TrajectoryOpenOptions>,
    ) -> Result<Self, TrajectoryIoError> {
        let options = options.into();
        let source = options.read_source(FORMAT);
        Self::open(reader, topology.into(), options.xyz, source)
    }

    /// Opens the stream like [`Self::new`] with the format options and
    /// source taken apart.
    pub(crate) fn open(
        reader: R,
        topology: Arc<Topology>,
        options: XyzReadOptions,
        source: ReadSource,
    ) -> Result<Self, TrajectoryIoError> {
        let source_label = &source.label;
        let limits = &source.limits;
        let atom_count = topology.atom_count();
        validate_atom_limit(atom_count, limits, source_label)?;
        let mut positions = Vec::new();
        positions.try_reserve_exact(atom_count).map_err(|_| {
            codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                source_label,
                Reason::Allocation(Part::PositionBuffer),
            )
        })?;
        Ok(Self {
            reader,
            topology,
            options,
            source,
            frame_cursor: 0,
            line: Vec::new(),
            positions,
            growing: false,
            extra_columns: false,
            index_extra_columns: 0,
        })
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.topology.as_ref()
    }

    /// Reads the next frame, or returns `false` at the end of the stream. A
    /// scan passes its control and the offset the frame starts at, which
    /// the blank lines after the last frame poll.
    pub(super) fn parse_next(
        &mut self,
        retain_positions: bool,
        scan: Option<ScanPoll<'_, '_>>,
    ) -> Result<bool, TrajectoryIoError> {
        let operation = if self.frame_cursor == 0 {
            TrajectoryIoOperation::ReadHeader
        } else {
            TrajectoryIoOperation::ReadFrame
        };
        if self.frame_cursor >= self.source.limits.max_frames() {
            if self.at_end(operation, scan)? {
                return Ok(false);
            }
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::LimitExceeded,
                operation,
                &self.source.label,
                self.frame_cursor,
                None,
            )
            .exceeding(
                Part::Trajectory,
                LimitExceeded::first_beyond(Resource::Frames, self.source.limits.max_frames()),
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
            self.source.limits.max_text_line_bytes(),
            position_bytes,
            &self.source.limits,
            &self.source.label,
            self.frame_cursor,
            operation,
            growing,
        )?
        else {
            return Ok(false);
        };
        if count_line.trim().is_empty() {
            // Blank lines may only trail the last frame.
            self.read_blank_tail(operation, scan.map(|scan| scan.after(consumed)))?;
            return Ok(false);
        }
        let mut frame_bytes = checked_frame_add(
            0,
            consumed,
            &self.source.label,
            self.frame_cursor,
            operation,
        )?;
        validate_frame_bytes(
            frame_bytes,
            &self.source.limits,
            &self.source.label,
            self.frame_cursor,
            operation,
        )?;
        let atom_count = count_line.trim().parse::<usize>().map_err(|_| {
            frame_error(
                FORMAT,
                TrajectoryIoErrorKind::InvalidFrame,
                operation,
                &self.source.label,
                self.frame_cursor,
                Reason::AtomCountLine,
            )
        })?;
        if atom_count == 0 {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::InvalidFrame,
                operation,
                &self.source.label,
                self.frame_cursor,
                Reason::ZeroAtoms,
            ));
        }
        validate_atom_limit(atom_count, &self.source.limits, &self.source.label)?;
        if atom_count != self.topology().atom_count() {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentAtomCount,
                operation,
                Some(FORMAT),
            )
            .with_source_label(&self.source.label)
            .with_frame(self.frame_cursor)
            .with_counts(self.topology().atom_count() as u64, atom_count as u64));
        }
        let Some((_, consumed)) = read_line(
            &mut self.reader,
            &mut self.line,
            self.source.limits.max_comment_bytes(),
            position_bytes,
            &self.source.limits,
            &self.source.label,
            self.frame_cursor,
            operation,
            growing,
        )?
        else {
            return Err(frame_error(
                FORMAT,
                TrajectoryIoErrorKind::TruncatedRecord,
                operation,
                &self.source.label,
                self.frame_cursor,
                Reason::NoCommentLine,
            ));
        };
        frame_bytes = checked_frame_add(
            frame_bytes,
            consumed,
            &self.source.label,
            self.frame_cursor,
            operation,
        )?;
        validate_frame_bytes(
            frame_bytes,
            &self.source.limits,
            &self.source.label,
            self.frame_cursor,
            operation,
        )?;

        // The indexed reader hands the buffer over, so the next read
        // allocates its own on the thread that reads.
        self.positions.clear();
        self.positions.reserve_exact(atom_count);
        self.extra_columns = false;
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
                self.source.limits.max_text_line_bytes(),
                position_bytes,
                &self.source.limits,
                &self.source.label,
                self.frame_cursor,
                operation,
                growing,
            )?
            else {
                return Err(frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::TruncatedRecord,
                    operation,
                    &self.source.label,
                    self.frame_cursor,
                    Reason::EndsBeforeAtomLine(atom_index),
                ));
            };
            frame_bytes = checked_frame_add(
                frame_bytes,
                consumed,
                &self.source.label,
                self.frame_cursor,
                operation,
            )?;
            validate_frame_bytes(
                frame_bytes,
                &self.source.limits,
                &self.source.label,
                self.frame_cursor,
                operation,
            )?;
            let atom = parse_atom_line(line).map_err(|error| {
                frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::InvalidFrame,
                    operation,
                    &self.source.label,
                    self.frame_cursor,
                    match error {
                        AtomLineError::MissingField => Reason::AtomLineFields(atom_index),
                        AtomLineError::NotNumber(_) => Reason::CoordinateNotNumber(atom_index),
                        AtomLineError::NotFinite(_) => {
                            Reason::CoordinateNotFiniteOnLine(atom_index)
                        }
                    },
                )
            })?;
            let element = element_of_token(atom.element).ok_or_else(|| {
                frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::InvalidFrame,
                    operation,
                    &self.source.label,
                    self.frame_cursor,
                    Reason::UnknownElement {
                        line: atom_index,
                        symbol: atom.element.into(),
                    },
                )
            })?;
            if element != expected_element {
                return Err(frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    operation,
                    &self.source.label,
                    self.frame_cursor,
                    Reason::ElementMismatch {
                        atom: atom_index,
                        found: element,
                        expected: expected_element,
                    },
                ));
            }
            if atom.extra.is_some() {
                if self.options.extra_columns != UnsupportedContent::Report {
                    return Err(frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::UnsupportedVariant,
                        operation,
                        &self.source.label,
                        self.frame_cursor,
                        Reason::ExtraColumns(atom_index),
                    ));
                }
                self.extra_columns = true;
            }
            if retain_positions {
                self.positions.push(atom.position);
            }
        }
        Ok(true)
    }

    /// Whether the stream holds no frame after the last one read: it ends,
    /// or holds only blank lines. A line that starts with no whitespace
    /// starts a frame and is not consumed. A scan's blank lines poll `scan`.
    pub(super) fn at_end(
        &mut self,
        operation: TrajectoryIoOperation,
        scan: Option<ScanPoll<'_, '_>>,
    ) -> Result<bool, TrajectoryIoError> {
        let next = loop {
            match self.reader.fill_buf() {
                Ok(available) => break available.first().copied(),
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                Err(error) => {
                    return Err(io_context(
                        operation,
                        Some(FORMAT),
                        &self.source.label,
                        error,
                    ))
                }
            }
        };
        let Some(next) = next else {
            return Ok(true);
        };
        if !next.is_ascii_whitespace() {
            return Ok(false);
        }
        // An atom-count line may begin with spaces; reading it tells.
        let position_bytes = filled_bytes(&self.positions, self.topology().atom_count());
        let Some((line, consumed)) = read_line(
            &mut self.reader,
            &mut self.line,
            self.source.limits.max_text_line_bytes(),
            position_bytes,
            &self.source.limits,
            &self.source.label,
            self.frame_cursor,
            operation,
            self.growing,
        )?
        else {
            return Ok(true);
        };
        if !line.trim().is_empty() {
            return Ok(false);
        }
        self.read_blank_tail(operation, scan.map(|scan| scan.after(consumed)))?;
        Ok(true)
    }

    /// Reads the lines after a blank line where a frame would start, which
    /// must all be blank too. After each line, a scan's control reports
    /// progress and fails a cancelled scan.
    fn read_blank_tail(
        &mut self,
        operation: TrajectoryIoOperation,
        mut scan: Option<ScanPoll<'_, '_>>,
    ) -> Result<(), TrajectoryIoError> {
        let position_bytes = filled_bytes(&self.positions, self.topology().atom_count());
        while let Some((line, consumed)) = read_line(
            &mut self.reader,
            &mut self.line,
            self.source.limits.max_text_line_bytes(),
            position_bytes,
            &self.source.limits,
            &self.source.label,
            self.frame_cursor,
            operation,
            self.growing,
        )? {
            if !line.trim().is_empty() {
                return Err(frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::InvalidFrame,
                    operation,
                    &self.source.label,
                    self.frame_cursor,
                    Reason::ContentAfterBlankLine,
                ));
            }
            if let Some(poll) = &mut scan {
                *poll = poll.after(consumed);
                poll.control.checkpoint(
                    self.frame_cursor,
                    poll.offset,
                    FORMAT,
                    &self.source.label,
                )?;
            }
        }
        Ok(())
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

/// The control of a scan, which the blank lines after the last frame poll,
/// and the offset the stream is at.
#[derive(Clone, Copy)]
pub(super) struct ScanPoll<'a, 'b> {
    pub(super) control: &'a IndexBuildControl<'b>,
    pub(super) offset: u64,
}

impl ScanPoll<'_, '_> {
    /// The poll at `consumed` bytes further on.
    fn after(self, consumed: usize) -> Self {
        Self {
            offset: self
                .offset
                .saturating_add(u64::try_from(consumed).unwrap_or(u64::MAX)),
            ..self
        }
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
        if !self.parse_next(true, None)? {
            return Ok(false);
        }
        self.publish(destination)?;
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            frame_error(
                FORMAT,
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::ReadFrame,
                &self.source.label,
                self.frame_cursor,
                Reason::AddressRange(Part::FrameNumber),
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
                Some(FORMAT),
                &self.source.label,
                error,
            )
        })?;
        let result = self.parse_next(false, None);
        self.reader.seek(SeekFrom::Start(start)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                &self.source.label,
                error,
            )
        })?;
        self.frame_cursor = 0;
        if !result? {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidHeader,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                &self.source.label,
                Reason::FileHoldsNoFrame,
            ));
        }
        Ok(())
    }

    pub(super) fn seek_to(
        &mut self,
        offset: SeekFrom,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.reader
            .seek(offset)
            .map_err(|error| io_context(operation, Some(FORMAT), &self.source.label, error))
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
            Some(FORMAT),
            source_label,
            None,
        )
        .exceeding(
            Part::AtomCount,
            LimitExceeded::new(Resource::Atoms, limits.max_atoms(), atom_count),
        ));
    }
    let scratch_bytes = atom_count
        .checked_mul(std::mem::size_of::<Point3>())
        .ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                Reason::AddressRange(Part::PositionBuffer),
            )
        })?;
    if scratch_bytes > limits.max_scratch_bytes() {
        return Err(codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
            source_label,
            None,
        )
        .exceeding(
            Part::PositionBuffer,
            LimitExceeded::new(
                Resource::BufferBytes,
                limits.max_scratch_bytes(),
                scratch_bytes,
            ),
        ));
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
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::LimitExceeded,
            operation,
            source_label,
            frame,
            None,
        )
        .exceeding(
            Part::Frame,
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
            frame_error(
                FORMAT,
                TrajectoryIoErrorKind::LimitExceeded,
                operation,
                source_label,
                frame,
                Reason::Range64(Part::FrameSize),
            )
        })
}

/// Reads one line and returns it without its LF, CRLF, or CR terminator,
/// together with the exact number of source bytes consumed, terminator
/// included.
///
/// The line, terminator included, must fit `limit`. The retained `line` buffer
/// grows only while its capacity plus `other_scratch_bytes` stays within
/// `max_scratch_bytes`, so whether a line fits depends only on its length.
/// With `terminated`, a line the stream ends inside, without its line break,
/// is truncated: its writer may not have finished it. A CR the stream ends
/// with is such a line too, since its LF may follow.
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
    let too_long = |needed: usize| {
        codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            operation,
            Some(FORMAT),
            source_label,
            None,
        )
        .exceeding(
            Part::Line,
            LimitExceeded::new(Resource::RecordBytes, limit, needed),
        )
    };
    line.clear();
    // The LF of a CRLF that the buffered chunk ends between, consumed but
    // not kept.
    let mut split_lf = false;
    // Whether the stream ends right after a CR.
    let mut cr_at_end = false;
    loop {
        let available = match reader.fill_buf() {
            Ok(available) => available,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(io_context(operation, Some(FORMAT), source_label, error)),
        };
        if available.is_empty() {
            break;
        }
        // Scan at most one byte past the limit, enough to detect an overlong line.
        let window = &available[..available.len().min((limit - line.len()).saturating_add(1))];
        let (chunk, complete) = match window
            .iter()
            .position(|byte| *byte == b'\n' || *byte == b'\r')
        {
            Some(end) if window[end] == b'\r' && window.get(end + 1) == Some(&b'\n') => {
                (&window[..end + 2], true)
            }
            Some(end) => (&window[..=end], true),
            None => (window, false),
        };
        let needed = line.len() + chunk.len();
        if needed > limit {
            return Err(too_long(needed));
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
                |reason| {
                    frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::LimitExceeded,
                        operation,
                        source_label,
                        frame,
                        Reason::LineBufferGrowth(Box::new(reason)),
                    )
                },
            )?;
        }
        line.extend_from_slice(chunk);
        let consumed = chunk.len();
        reader.consume(consumed);
        if complete {
            if line.ends_with(b"\r") {
                // The LF of a CRLF may start the next buffered chunk.
                let next = loop {
                    match reader.fill_buf() {
                        Ok(next) => break next,
                        Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                        Err(error) => {
                            return Err(io_context(operation, Some(FORMAT), source_label, error))
                        }
                    }
                };
                match next.first() {
                    Some(b'\n') => {
                        if line.len() + 1 > limit {
                            return Err(too_long(line.len() + 1));
                        }
                        reader.consume(1);
                        split_lf = true;
                    }
                    Some(_) => {}
                    None => cr_at_end = true,
                }
            }
            break;
        }
    }
    let read = line.len() + usize::from(split_lf);
    if read == 0 {
        return Ok(None);
    }
    let broken = line.ends_with(b"\n") || (line.ends_with(b"\r") && !cr_at_end);
    if terminated && !broken {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::TruncatedRecord,
            operation,
            source_label,
            frame,
            Reason::NoLineBreak,
        ));
    }
    let mut content = line.as_slice();
    content = content.strip_suffix(b"\n").unwrap_or(content);
    content = content.strip_suffix(b"\r").unwrap_or(content);
    let content = std::str::from_utf8(content).map_err(|_| {
        io_context(
            operation,
            Some(FORMAT),
            source_label,
            io::Error::new(io::ErrorKind::InvalidData, "The stream is not valid UTF-8."),
        )
    })?;
    Ok(Some((content, read)))
}
