//! Strict bounded multi-frame XYZ trajectory I/O.

use std::io::{self, BufRead, Read, Seek, SeekFrom, Write};
use std::sync::Arc;

use crate::{
    FrameBuffer, FrameBufferData, SeekableTrajectoryReader, TrajectoryCodecErrorContext,
    TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat, TrajectoryFrameView,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use kekule::core::Element;
use kekule::geometry::Point3;
use kekule::topology::Topology;
use kekule::units::{Quantity, Unit, ANGSTROM, MODEL_LENGTH_UNIT};

use super::{
    codec_context, io_context, projected_index_limit, require_nonempty_writer,
    reserve_index_for_push, TrajectoryIoLimits, TrajectoryTopologyBinding,
};

const MAX_WRITER_COMMENT_BYTES: usize = 1_048_576;

/// XYZ length-unit interpretation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XyzReadOptions {
    length_unit: Unit,
}

impl Default for XyzReadOptions {
    fn default() -> Self {
        Self {
            length_unit: ANGSTROM,
        }
    }
}

impl XyzReadOptions {
    pub const fn with_length_unit(mut self, length_unit: Unit) -> Self {
        self.length_unit = length_unit;
        self
    }

    pub const fn length_unit(&self) -> Unit {
        self.length_unit
    }
}

/// Deterministic strict XYZ writer options.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzWriteOptions {
    length_unit: Unit,
    decimal_places: usize,
    comment: String,
}

impl Default for XyzWriteOptions {
    fn default() -> Self {
        Self {
            length_unit: ANGSTROM,
            decimal_places: 8,
            comment: "written by kekule-traj".into(),
        }
    }
}

impl XyzWriteOptions {
    pub const fn with_length_unit(mut self, length_unit: Unit) -> Self {
        self.length_unit = length_unit;
        self
    }

    pub const fn with_decimal_places(mut self, decimal_places: usize) -> Self {
        self.decimal_places = decimal_places;
        self
    }

    pub fn with_comment(mut self, comment: impl Into<String>) -> Self {
        self.comment = comment.into();
        self
    }
}

/// Sequential XYZ reader over a caller-supplied buffered stream.
pub struct XyzReader<R> {
    reader: R,
    binding: TrajectoryTopologyBinding,
    options: XyzReadOptions,
    limits: TrajectoryIoLimits,
    source_label: String,
    frame_cursor: u64,
    line: String,
    positions: Vec<Point3>,
}

impl<R: BufRead> XyzReader<R> {
    pub fn new(
        reader: R,
        binding: TrajectoryTopologyBinding,
        options: XyzReadOptions,
        limits: TrajectoryIoLimits,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryError> {
        let source_label = source_label.into();
        options
            .length_unit
            .conversion_factor_to(MODEL_LENGTH_UNIT)
            .map_err(|error| {
                codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Open,
                    Some(TrajectoryFormat::Xyz),
                    &source_label,
                    format!("XYZ length unit is incompatible: {error}"),
                )
            })?;
        let atom_count = binding.topology().atom_count();
        validate_atom_limit(atom_count, &limits, &source_label)?;
        let mut positions = Vec::new();
        positions.try_reserve_exact(atom_count).map_err(|_| {
            codec_context(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                &source_label,
                "could not reserve XYZ position scratch",
            )
        })?;
        Ok(Self {
            reader,
            binding,
            options,
            limits,
            source_label,
            frame_cursor: 0,
            line: String::new(),
            positions,
        })
    }

    pub fn topology(&self) -> &Topology {
        self.binding.topology()
    }

    fn parse_next(&mut self, retain_positions: bool) -> Result<bool, TrajectoryError> {
        let operation = if self.frame_cursor == 0 {
            TrajectoryIoOperation::ReadHeader
        } else {
            TrajectoryIoOperation::ReadFrame
        };
        if self.frame_cursor >= self.limits.max_frames {
            if self
                .reader
                .fill_buf()
                .map_err(|error| {
                    io_context(
                        operation,
                        Some(TrajectoryFormat::Xyz),
                        &self.source_label,
                        error,
                    )
                })?
                .is_empty()
            {
                return Ok(false);
            }
            return Err(frame_codec_error(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                operation,
                &self.source_label,
                self.frame_cursor,
                "XYZ frame count exceeds the configured limit",
            ));
        }
        let Some(count_line) = read_line(
            &mut self.reader,
            &mut self.line,
            self.limits.max_text_line_bytes,
            &self.source_label,
            operation,
        )?
        else {
            return Ok(false);
        };
        let mut frame_bytes = u64::try_from(count_line.len())
            .ok()
            .and_then(|bytes| bytes.checked_add(1))
            .ok_or_else(|| {
                frame_codec_error(
                    TrajectoryCodecErrorKind::ResourceLimitExceeded,
                    operation,
                    &self.source_label,
                    self.frame_cursor,
                    "XYZ frame byte count overflows",
                )
            })?;
        validate_frame_bytes(
            frame_bytes,
            &self.limits,
            &self.source_label,
            self.frame_cursor,
            operation,
        )?;
        let atom_count = count_line.trim().parse::<usize>().map_err(|_| {
            frame_codec_error(
                TrajectoryCodecErrorKind::InvalidFrame,
                operation,
                &self.source_label,
                self.frame_cursor,
                "XYZ atom-count line is not one unsigned integer",
            )
        })?;
        if atom_count == 0 {
            return Err(frame_codec_error(
                TrajectoryCodecErrorKind::InvalidFrame,
                operation,
                &self.source_label,
                self.frame_cursor,
                "XYZ frames cannot contain zero atoms",
            ));
        }
        validate_atom_limit(atom_count, &self.limits, &self.source_label)?;
        if atom_count != self.topology().atom_count() {
            return Err(TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::InconsistentAtomCount,
                operation,
                Some(TrajectoryFormat::Xyz),
            )
            .with_source_label(&self.source_label)
            .with_frame(self.frame_cursor)
            .with_counts(self.topology().atom_count() as u64, atom_count as u64)
            .into());
        }
        if read_line(
            &mut self.reader,
            &mut self.line,
            self.limits.max_comment_bytes,
            &self.source_label,
            operation,
        )?
        .is_none()
        {
            return Err(frame_codec_error(
                TrajectoryCodecErrorKind::TruncatedRecord,
                operation,
                &self.source_label,
                self.frame_cursor,
                "XYZ frame is missing its comment line",
            ));
        }
        frame_bytes = checked_frame_add(
            frame_bytes,
            self.line.len(),
            &self.source_label,
            self.frame_cursor,
            operation,
        )?;
        validate_frame_bytes(
            frame_bytes,
            &self.limits,
            &self.source_label,
            self.frame_cursor,
            operation,
        )?;

        self.positions.clear();
        for atom_index in 0..atom_count {
            let atom_id = self.topology().atom_ids()[atom_index];
            let expected_element = self
                .topology()
                .atom(atom_id)
                .map_err(|error| {
                    frame_codec_error(
                        TrajectoryCodecErrorKind::InconsistentMetadata,
                        operation,
                        &self.source_label,
                        self.frame_cursor,
                        format!("topology atom lookup failed: {error}"),
                    )
                })?
                .element;
            let Some(line) = read_line(
                &mut self.reader,
                &mut self.line,
                self.limits.max_text_line_bytes,
                &self.source_label,
                operation,
            )?
            else {
                return Err(frame_codec_error(
                    TrajectoryCodecErrorKind::TruncatedRecord,
                    operation,
                    &self.source_label,
                    self.frame_cursor,
                    format!("XYZ frame ends before atom line {atom_index}"),
                ));
            };
            frame_bytes = checked_frame_add(
                frame_bytes,
                line.len(),
                &self.source_label,
                self.frame_cursor,
                operation,
            )?;
            validate_frame_bytes(
                frame_bytes,
                &self.limits,
                &self.source_label,
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
                    TrajectoryCodecErrorKind::InvalidFrame,
                    operation,
                    &self.source_label,
                    self.frame_cursor,
                    format!("XYZ atom line {atom_index} must contain exactly four fields"),
                ));
            };
            let element = Element::from_symbol(symbol).ok_or_else(|| {
                frame_codec_error(
                    TrajectoryCodecErrorKind::InvalidFrame,
                    operation,
                    &self.source_label,
                    self.frame_cursor,
                    format!("XYZ atom line {atom_index} has invalid element symbol {symbol:?}"),
                )
            })?;
            if element != expected_element {
                return Err(frame_codec_error(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    operation,
                    &self.source_label,
                    self.frame_cursor,
                    format!(
                        "XYZ atom {atom_index} is {}, but topology order requires {}",
                        element.symbol(),
                        expected_element.symbol()
                    ),
                ));
            }
            let parse_coordinate = |value: &str| {
                value.parse::<f64>().map_err(|_| {
                    frame_codec_error(
                        TrajectoryCodecErrorKind::InvalidFrame,
                        operation,
                        &self.source_label,
                        self.frame_cursor,
                        format!("XYZ atom line {atom_index} has an invalid coordinate"),
                    )
                })
            };
            let x = parse_coordinate(x)?;
            let y = parse_coordinate(y)?;
            let z = parse_coordinate(z)?;
            if ![x, y, z].into_iter().all(f64::is_finite) {
                return Err(frame_codec_error(
                    TrajectoryCodecErrorKind::InvalidFrame,
                    operation,
                    &self.source_label,
                    self.frame_cursor,
                    format!("XYZ atom line {atom_index} has a non-finite coordinate"),
                ));
            }
            if retain_positions {
                self.positions.push(Point3::new(x, y, z));
            }
        }
        Ok(true)
    }

    fn publish(&self, destination: &mut FrameBuffer) -> Result<(), TrajectoryError> {
        destination
            .replace_from_data(FrameBufferData::new(
                self.binding.topology_arc(),
                Quantity::new(self.positions.as_slice(), self.options.length_unit),
            ))
            .map_err(Into::into)
    }
}

impl<R: BufRead> TrajectoryReader for XyzReader<R> {
    fn topology(&self) -> &Topology {
        self.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.binding.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        if !std::ptr::eq(self.topology(), destination.topology()) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        if !self.parse_next(true)? {
            return Ok(false);
        }
        self.publish(destination)?;
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            frame_codec_error(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::ReadFrame,
                &self.source_label,
                self.frame_cursor,
                "XYZ frame cursor overflow",
            )
        })?;
        Ok(true)
    }
}

impl<R: BufRead + Seek> XyzReader<R> {
    pub fn validate_first_frame(&mut self) -> Result<(), TrajectoryError> {
        let start = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                error,
            )
        })?;
        let result = self.parse_next(false);
        self.reader.seek(SeekFrom::Start(start)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                error,
            )
        })?;
        self.frame_cursor = 0;
        if !result? {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidHeader,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                "XYZ trajectory contains no frame",
            ));
        }
        Ok(())
    }

    pub fn into_indexed(mut self) -> Result<IndexedXyzReader<R>, TrajectoryError> {
        self.reader.seek(SeekFrom::Start(0)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                error,
            )
        })?;
        self.frame_cursor = 0;
        let mut offsets = Vec::new();
        loop {
            if let Some(limit) = projected_index_limit(offsets.len(), &self.limits) {
                if self
                    .reader
                    .fill_buf()
                    .map_err(|error| {
                        io_context(
                            TrajectoryIoOperation::Index,
                            Some(TrajectoryFormat::Xyz),
                            &self.source_label,
                            error,
                        )
                    })?
                    .is_empty()
                {
                    break;
                }
                return Err(codec_context(
                    TrajectoryCodecErrorKind::ResourceLimitExceeded,
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xyz),
                    &self.source_label,
                    format!("XYZ index {limit} exceeds the configured limit"),
                ));
            }
            let offset = self.reader.stream_position().map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xyz),
                    &self.source_label,
                    error,
                )
            })?;
            if !self.parse_next(false)? {
                break;
            }
            reserve_index_for_push(
                &mut offsets,
                &self.limits,
                TrajectoryFormat::Xyz,
                &self.source_label,
                self.frame_cursor,
            )?;
            offsets.push(offset);
            self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
                frame_codec_error(
                    TrajectoryCodecErrorKind::ResourceLimitExceeded,
                    TrajectoryIoOperation::Index,
                    &self.source_label,
                    self.frame_cursor,
                    "XYZ frame cursor overflows",
                )
            })?;
        }
        if offsets.is_empty() {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidHeader,
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                "XYZ trajectory contains no frame",
            ));
        }
        self.reader.seek(SeekFrom::Start(0)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                error,
            )
        })?;
        self.frame_cursor = 0;
        Ok(IndexedXyzReader {
            reader: self,
            offsets,
        })
    }
}

/// Fully scanned XYZ reader with verified frame offsets.
pub struct IndexedXyzReader<R> {
    reader: XyzReader<R>,
    offsets: Vec<u64>,
}

impl<R: BufRead + Seek> TrajectoryReader for IndexedXyzReader<R> {
    fn topology(&self) -> &Topology {
        self.reader.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.reader.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        self.reader.read_next(destination)
    }
}

impl<R: BufRead + Seek> SeekableTrajectoryReader for IndexedXyzReader<R> {
    fn frame_count(&self) -> Option<u64> {
        u64::try_from(self.offsets.len()).ok()
    }

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        if !std::ptr::eq(self.topology(), destination.topology()) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let index_usize =
            usize::try_from(index).map_err(|_| TrajectoryError::FrameIndexOutOfRange(index))?;
        let offset = *self
            .offsets
            .get(index_usize)
            .ok_or(TrajectoryError::FrameIndexOutOfRange(index))?;
        let saved_offset = self.reader.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Xyz),
                &self.reader.source_label,
                error,
            )
        })?;
        let saved_cursor = self.reader.frame_cursor;
        self.reader
            .reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadFrame,
                    Some(TrajectoryFormat::Xyz),
                    &self.reader.source_label,
                    error,
                )
            })?;
        self.reader.frame_cursor = index;
        let parsed = self.reader.parse_next(true);
        let restore = self
            .reader
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadFrame,
                    Some(TrajectoryFormat::Xyz),
                    &self.reader.source_label,
                    error,
                )
            });
        self.reader.frame_cursor = saved_cursor;
        if !parsed? {
            return Err(TrajectoryError::FrameIndexOutOfRange(index));
        }
        restore?;
        self.reader.publish(destination)
    }
}

/// Strict multi-frame XYZ writer over a caller-supplied stream.
pub struct XyzWriter<W> {
    writer: W,
    topology: Arc<Topology>,
    options: XyzWriteOptions,
    source_label: String,
    frame_count: u64,
}

impl<W: Write> XyzWriter<W> {
    pub fn new(
        writer: W,
        topology: Arc<Topology>,
        options: XyzWriteOptions,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryError> {
        let source_label = source_label.into();
        options
            .length_unit
            .conversion_factor_to(MODEL_LENGTH_UNIT)
            .map_err(|error| {
                codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Open,
                    Some(TrajectoryFormat::Xyz),
                    &source_label,
                    format!("XYZ length unit is incompatible: {error}"),
                )
            })?;
        if options.decimal_places > 17 {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidPrecision,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                &source_label,
                "XYZ decimal precision must be at most 17 places",
            ));
        }
        if options.comment.contains(['\r', '\n']) {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                &source_label,
                "XYZ writer comment cannot contain a line break",
            ));
        }
        if options.comment.len() > MAX_WRITER_COMMENT_BYTES {
            return Err(codec_context(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                &source_label,
                format!("XYZ writer comment exceeds the {MAX_WRITER_COMMENT_BYTES}-byte limit"),
            ));
        }
        if topology.atom_count() == 0 {
            return Err(codec_context(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                &source_label,
                "XYZ writer atom count must be positive",
            ));
        }
        Ok(Self {
            writer,
            topology,
            options,
            source_label,
            frame_count: 0,
        })
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn writer(&self) -> &W {
        &self.writer
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }

    pub(crate) fn validate_finish(&self) -> Result<(), TrajectoryError> {
        require_nonempty_writer(self.frame_count, TrajectoryFormat::Xyz, &self.source_label)
    }

    /// Flushes and returns the completed nonempty XYZ stream.
    pub fn finish(mut self) -> Result<W, TrajectoryError> {
        self.validate_finish()?;
        self.flush().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Finish,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                error,
            )
        })?;
        Ok(self.writer)
    }
}

impl<W: Write> TrajectoryWriter for XyzWriter<W> {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if !std::ptr::eq(self.topology(), frame.topology()) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let configuration = frame.configuration();
        let unsupported = [
            (configuration.cell().is_some(), "periodic cell"),
            (frame.velocities().is_some(), "velocities"),
            (frame.forces().is_some(), "forces"),
            (frame.time().is_some(), "time"),
            (frame.step().is_some(), "step"),
            (frame.observation().is_some(), "structure observation"),
            (!frame.props().is_empty(), "frame properties"),
        ];
        if let Some((_, field)) = unsupported.into_iter().find(|(present, _)| *present) {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                field,
            ));
        }
        let factor = MODEL_LENGTH_UNIT
            .conversion_factor_to(self.options.length_unit)
            .map_err(|error| {
                codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Xyz),
                    &self.source_label,
                    format!("XYZ output unit is incompatible: {error}"),
                )
            })?;
        let positions = configuration.positions().values();
        for point in positions.value().iter() {
            if !Point3::new(point.x * factor, point.y * factor, point.z * factor).is_finite() {
                return Err(codec_context(
                    TrajectoryCodecErrorKind::InvalidFrame,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Xyz),
                    &self.source_label,
                    "XYZ output coordinate is not finite after conversion",
                ));
            }
        }

        writeln!(self.writer, "{}", self.topology.atom_count()).map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                error,
            )
        })?;
        writeln!(self.writer, "{}", self.options.comment).map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                error,
            )
        })?;
        for (atom_index, point) in positions.value().iter().enumerate() {
            let atom_id = self.topology.atom_ids()[atom_index];
            let element = self.topology.atom(atom_id).map_err(|error| {
                codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Xyz),
                    &self.source_label,
                    format!("topology atom lookup failed: {error}"),
                )
            })?;
            writeln!(
                self.writer,
                "{} {x:.precision$} {y:.precision$} {z:.precision$}",
                element.element.symbol(),
                x = point.x * factor,
                y = point.y * factor,
                z = point.z * factor,
                precision = self.options.decimal_places,
            )
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Xyz),
                    &self.source_label,
                    error,
                )
            })?;
        }
        self.frame_count = self.frame_count.checked_add(1).ok_or_else(|| {
            codec_context(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xyz),
                &self.source_label,
                "XYZ writer frame count overflow",
            )
        })?;
        Ok(())
    }
}

fn validate_atom_limit(
    atom_count: usize,
    limits: &TrajectoryIoLimits,
    source_label: &str,
) -> Result<(), TrajectoryError> {
    if atom_count > limits.max_atoms {
        return Err(codec_context(
            TrajectoryCodecErrorKind::ResourceLimitExceeded,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Xyz),
            source_label,
            format!(
                "XYZ atom count {atom_count} exceeds configured maximum {}",
                limits.max_atoms
            ),
        ));
    }
    let scratch_bytes = atom_count
        .checked_mul(std::mem::size_of::<Point3>())
        .ok_or_else(|| {
            codec_context(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xyz),
                source_label,
                "XYZ position scratch size overflows",
            )
        })?;
    if scratch_bytes > limits.max_scratch_bytes {
        return Err(codec_context(
            TrajectoryCodecErrorKind::ResourceLimitExceeded,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Xyz),
            source_label,
            format!(
                "XYZ position scratch requires {scratch_bytes} bytes, exceeding configured maximum {}",
                limits.max_scratch_bytes
            ),
        ));
    }
    Ok(())
}

fn validate_frame_bytes(
    frame_bytes: u64,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    frame: u64,
    operation: TrajectoryIoOperation,
) -> Result<(), TrajectoryError> {
    if frame_bytes > limits.max_frame_bytes {
        return Err(frame_codec_error(
            TrajectoryCodecErrorKind::ResourceLimitExceeded,
            operation,
            source_label,
            frame,
            format!(
                "XYZ frame exceeds configured {}-byte limit",
                limits.max_frame_bytes
            ),
        ));
    }
    Ok(())
}

fn checked_frame_add(
    frame_bytes: u64,
    line_bytes: usize,
    source_label: &str,
    frame: u64,
    operation: TrajectoryIoOperation,
) -> Result<u64, TrajectoryError> {
    u64::try_from(line_bytes)
        .ok()
        .and_then(|bytes| bytes.checked_add(1))
        .and_then(|bytes| frame_bytes.checked_add(bytes))
        .ok_or_else(|| {
            frame_codec_error(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                operation,
                source_label,
                frame,
                "XYZ frame byte count overflows",
            )
        })
}

fn read_line<'a, R: BufRead>(
    reader: &mut R,
    line: &'a mut String,
    limit: usize,
    source_label: &str,
    operation: TrajectoryIoOperation,
) -> Result<Option<&'a str>, TrajectoryError> {
    line.clear();
    let bounded = limit.checked_add(1).ok_or_else(|| {
        codec_context(
            TrajectoryCodecErrorKind::ResourceLimitExceeded,
            operation,
            Some(TrajectoryFormat::Xyz),
            source_label,
            "XYZ line limit overflows",
        )
    })?;
    let read = Read::by_ref(reader)
        .take(bounded as u64)
        .read_line(line)
        .map_err(|error| io_context(operation, Some(TrajectoryFormat::Xyz), source_label, error))?;
    if read == 0 {
        return Ok(None);
    }
    if read > limit {
        return Err(codec_context(
            TrajectoryCodecErrorKind::ResourceLimitExceeded,
            operation,
            Some(TrajectoryFormat::Xyz),
            source_label,
            format!("XYZ line exceeds configured {limit}-byte limit"),
        ));
    }
    if line.ends_with('\n') {
        line.pop();
        if line.ends_with('\r') {
            line.pop();
        }
    }
    Ok(Some(line.as_str()))
}

fn frame_codec_error(
    kind: TrajectoryCodecErrorKind,
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: u64,
    detail: impl Into<String>,
) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(kind, operation, Some(TrajectoryFormat::Xyz))
        .with_source_label(source_label)
        .with_frame(frame)
        .with_detail(detail)
        .into()
}

fn writer_field_error(source_label: &str, frame: u64, field: &str) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::UnsupportedField,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Xyz),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_detail(format!("XYZ cannot preserve {field}"))
    .into()
}
