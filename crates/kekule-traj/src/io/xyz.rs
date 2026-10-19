//! Strict bounded multi-frame XYZ trajectory I/O.
//!
//! [`XyzReader`] and [`IndexedXyzReader`] read constant-count
//! element/x/y/z frames whose atom lines follow the topology's dense order and
//! elements; [`XyzWriter`] writes them deterministically. XYZ stores
//! positions only: the writer rejects cells, velocities, forces, time, step
//! and properties instead of dropping them. To obtain a topology from an XYZ
//! file itself, interpret one of its frames with `kekule::xyz`, which keeps
//! atom-line order as dense order and infers bonds from distances. This
//! reader accepts element symbols in their standard letter case only, while
//! interpretation also accepts other cases and atomic numbers.
//!
//! # Examples
//!
//! ```
//! use std::io::Cursor;
//! use std::sync::Arc;
//!
//! use kekule::geometry::Point3;
//! use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
//! use kekule::units::{Quantity, ANGSTROM};
//! use kekule::{smiles, topology::Topology};
//! use kekule_traj::io::xyz::{XyzReadOptions, XyzReader, XyzWriteOptions, XyzWriter};
//! use kekule_traj::{TrajectoryReader, TrajectoryWriter};
//!
//! let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//! let positions = Positions::from_vec(Quantity::new(vec![Point3::new(1.5, 0.0, 0.0)], ANGSTROM))?;
//! let frame = TrajectoryFrame::new(positions);
//! let trajectory = Trajectory::from_items(Arc::clone(&topology), [frame])?;
//!
//! let options = XyzWriteOptions::default().with_decimal_places(3).with_comment("helium");
//! let mut writer = XyzWriter::new(Vec::new(), Arc::clone(&topology), options, "memory")?;
//! writer.write_frame(trajectory.get(0).unwrap())?;
//! let bytes = writer.finish()?;
//! assert_eq!(String::from_utf8(bytes.clone())?, "1\nhelium\nHe 1.500 0.000 0.000\n");
//!
//! let mut reader = XyzReader::new(Cursor::new(bytes), topology, XyzReadOptions::default())?;
//! let mut frame = reader.frame_buffer();
//! assert!(reader.read_next(&mut frame)?);
//! let position = frame.positions().get(0).unwrap().value_in(ANGSTROM)?;
//! assert_eq!(position, Point3::new(1.5, 0.0, 0.0));
//! assert!(!reader.read_next(&mut frame)?);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::TrajectoryField;
use kekule::diagnostics::{LimitExceeded, Resource};
use std::fmt;
use std::io::{self, BufRead, Seek, SeekFrom, Write};
use std::sync::Arc;

use crate::{
    FrameBuffer, FrameData, SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoError,
    TrajectoryIoErrorKind, TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use kekule::geometry::Point3;
use kekule::molecule::Element;
use kekule::structure::TrajectoryFrameView;
use kekule::topology::Topology;
use kekule::units::{Quantity, Unit, ANGSTROM, CANONICAL_LENGTH_UNIT};

use super::limits::{reserve_scratch, scratch_bytes};
use super::{
    codec_context, io_context, projected_index_limit, require_no_partial_frame,
    require_nonempty_writer, reserve_index_for_push, TrajectoryIoLimits,
};

const MAX_WRITER_COMMENT_BYTES: usize = 1_048_576;

/// XYZ length-unit interpretation.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzReadOptions {
    length_unit: Unit,
    limits: TrajectoryIoLimits,
    source_label: String,
}

impl Default for XyzReadOptions {
    fn default() -> Self {
        Self {
            limits: TrajectoryIoLimits::default(),
            source_label: "xyz stream".into(),
            length_unit: ANGSTROM,
        }
    }
}

impl XyzReadOptions {
    /// Sets limits applied before allocation, scanning, and seeking.
    /// `max_scratch_bytes` bounds the combined capacities of the dense position
    /// array and the reusable line buffer, which grows to hold the longest
    /// line read so far, terminator included, in both sequential and indexed
    /// modes. The frame index has its own `max_index_bytes` limit.
    #[must_use]
    pub fn with_limits(mut self, limits: TrajectoryIoLimits) -> Self {
        self.limits = limits;
        self
    }

    /// Sets the source label used in decoding errors.
    #[must_use]
    pub fn with_source_label(mut self, source_label: impl Into<String>) -> Self {
        self.source_label = source_label.into();
        self
    }

    /// Sets [`Self::length_unit`].
    #[must_use]
    pub const fn with_length_unit(mut self, length_unit: Unit) -> Self {
        self.length_unit = length_unit;
        self
    }

    /// The length unit of file coordinates; angstrom by default. Positions are
    /// converted to the canonical length unit on publication.
    pub const fn length_unit(&self) -> Unit {
        self.length_unit
    }

    /// Limits applied before allocation, scanning, and seeking.
    pub const fn limits(&self) -> &TrajectoryIoLimits {
        &self.limits
    }

    /// The source label used in decoding errors.
    pub fn source_label(&self) -> &str {
        &self.source_label
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
    /// Sets [`Self::length_unit`].
    #[must_use]
    pub const fn with_length_unit(mut self, length_unit: Unit) -> Self {
        self.length_unit = length_unit;
        self
    }

    /// Sets [`Self::decimal_places`].
    #[must_use]
    pub const fn with_decimal_places(mut self, decimal_places: usize) -> Self {
        self.decimal_places = decimal_places;
        self
    }

    /// Sets [`Self::comment`].
    #[must_use]
    pub fn with_comment(mut self, comment: impl Into<String>) -> Self {
        self.comment = comment.into();
        self
    }

    /// The length unit coordinates are written in; angstrom by default.
    pub const fn length_unit(&self) -> Unit {
        self.length_unit
    }

    /// Fixed decimal places per coordinate, at most 17; 8 by default.
    pub const fn decimal_places(&self) -> usize {
        self.decimal_places
    }

    /// The comment line written in every frame, without line breaks and at
    /// most 1 MiB; `"written by kekule-traj"` by default.
    pub fn comment(&self) -> &str {
        &self.comment
    }
}

/// Sequential XYZ reader over a caller-supplied buffered stream.
pub struct XyzReader<R> {
    reader: R,
    topology: Arc<Topology>,
    options: XyzReadOptions,
    frame_cursor: usize,
    line: Vec<u8>,
    positions: Vec<Point3>,
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
                    format!("XYZ length unit is incompatible: {error}"),
                )
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
                "could not reserve XYZ position scratch",
            )
        })?;
        Ok(Self {
            reader,
            topology,
            options,
            frame_cursor: 0,
            line: Vec::new(),
            positions,
        })
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.topology.as_ref()
    }

    fn parse_next(&mut self, retain_positions: bool) -> Result<bool, TrajectoryIoError> {
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
                "XYZ frame count exceeds the configured limit",
            ));
        }
        // The position array never grows after construction; the line buffer
        // may grow only within what it leaves of the scratch budget.
        let position_bytes = scratch_bytes(&self.positions);
        let Some((count_line, consumed)) = read_line(
            &mut self.reader,
            &mut self.line,
            self.options.limits.max_text_line_bytes(),
            position_bytes,
            &self.options.limits,
            &self.options.source_label,
            self.frame_cursor,
            operation,
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
                "XYZ atom-count line is not one unsigned integer",
            )
        })?;
        if atom_count == 0 {
            return Err(frame_codec_error(
                TrajectoryIoErrorKind::InvalidFrame,
                operation,
                &self.options.source_label,
                self.frame_cursor,
                "XYZ frames cannot contain zero atoms",
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
        )?
        else {
            return Err(frame_codec_error(
                TrajectoryIoErrorKind::TruncatedRecord,
                operation,
                &self.options.source_label,
                self.frame_cursor,
                "XYZ frame is missing its comment line",
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

        self.positions.clear();
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
            )?
            else {
                return Err(frame_codec_error(
                    TrajectoryIoErrorKind::TruncatedRecord,
                    operation,
                    &self.options.source_label,
                    self.frame_cursor,
                    format!("XYZ frame ends before atom line {atom_index}"),
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
                    format!("XYZ atom line {atom_index} must contain exactly four fields"),
                ));
            };
            let element = Element::from_symbol(symbol).ok_or_else(|| {
                frame_codec_error(
                    TrajectoryIoErrorKind::InvalidFrame,
                    operation,
                    &self.options.source_label,
                    self.frame_cursor,
                    format!("XYZ atom line {atom_index} has invalid element symbol {symbol:?}"),
                )
            })?;
            if element != expected_element {
                return Err(frame_codec_error(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    operation,
                    &self.options.source_label,
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
                        TrajectoryIoErrorKind::InvalidFrame,
                        operation,
                        &self.options.source_label,
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
                    TrajectoryIoErrorKind::InvalidFrame,
                    operation,
                    &self.options.source_label,
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

    fn publish(&self, destination: &mut FrameBuffer) -> Result<(), TrajectoryIoError> {
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
                "XYZ frame cursor overflow",
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
                "XYZ trajectory contains no frame",
            ));
        }
        Ok(())
    }

    /// Indexes a freshly opened reader, whose next frame is frame 0.
    fn index_frames(mut self, stream_start: u64) -> Result<IndexedXyzReader<R>, TrajectoryIoError> {
        let mut offsets = Vec::new();
        loop {
            if let Some(limit) = projected_index_limit(offsets.len(), &self.options.limits) {
                if self
                    .reader
                    .fill_buf()
                    .map_err(|error| {
                        io_context(
                            TrajectoryIoOperation::Index,
                            Some(TrajectoryFormat::Xyz),
                            &self.options.source_label,
                            error,
                        )
                    })?
                    .is_empty()
                {
                    break;
                }
                return Err(codec_context(
                    TrajectoryIoErrorKind::LimitExceeded,
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xyz),
                    &self.options.source_label,
                    format!("XYZ index {limit} exceeds the configured limit"),
                ));
            }
            let offset = self.reader.stream_position().map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xyz),
                    &self.options.source_label,
                    error,
                )
            })?;
            if !self.parse_next(false)? {
                break;
            }
            reserve_index_for_push(
                &mut offsets,
                &self.options.limits,
                TrajectoryFormat::Xyz,
                &self.options.source_label,
                self.frame_cursor,
            )?;
            offsets.push(offset);
            self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
                frame_codec_error(
                    TrajectoryIoErrorKind::LimitExceeded,
                    TrajectoryIoOperation::Index,
                    &self.options.source_label,
                    self.frame_cursor,
                    "XYZ frame cursor overflows",
                )
            })?;
        }
        if offsets.is_empty() {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidHeader,
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Xyz),
                &self.options.source_label,
                "XYZ trajectory contains no frame",
            ));
        }
        self.reader
            .seek(SeekFrom::Start(stream_start))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xyz),
                    &self.options.source_label,
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

impl<R> fmt::Debug for IndexedXyzReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedXyzReader")
            .field("reader", &self.reader)
            .field("frame_count", &self.offsets.len())
            .finish_non_exhaustive()
    }
}

impl<R: BufRead + Seek> IndexedXyzReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order
    /// and indexes it for random access.
    ///
    /// The stream's position at construction is frame 0. Every frame from there to
    /// the end is decoded and verified before the reader is returned, and both
    /// sequential and indexed reads then address frames from that start. Topology,
    /// limits, and diagnostic labels follow [`XyzReader::new`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`XyzReader::new`], the decoding error of the
    /// first invalid frame, [`TrajectoryIoErrorKind::LimitExceeded`] when the
    /// index outgrows the limits, and [`TrajectoryIoErrorKind::Io`] when
    /// reading or seeking fails.
    pub fn new(
        mut reader: R,
        topology: impl Into<Arc<Topology>>,
        options: XyzReadOptions,
    ) -> Result<Self, TrajectoryIoError> {
        let stream_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                &options.source_label,
                error,
            )
        })?;
        XyzReader::new(reader, topology, options)?.index_frames(stream_start)
    }
}

impl<R: BufRead + Seek> TrajectoryReader for IndexedXyzReader<R> {
    fn topology(&self) -> &Topology {
        self.reader.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.reader.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        self.reader.read_next(destination)
    }
}

impl<R: BufRead + Seek> SeekableTrajectoryReader for IndexedXyzReader<R> {
    fn frame_count(&self) -> Option<usize> {
        Some(self.offsets.len())
    }

    fn read_frame(
        &mut self,
        index: usize,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        if !self.topology().shares_layout(destination.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        let index_usize = index;
        let offset = *self
            .offsets
            .get(index_usize)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        let saved_offset = self.reader.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Xyz),
                &self.reader.options.source_label,
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
                    &self.reader.options.source_label,
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
                    &self.reader.options.source_label,
                    error,
                )
            });
        self.reader.frame_cursor = saved_cursor;
        if !parsed? {
            return Err(TrajectoryIoError::frame_index_out_of_range(index));
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
    frame_count: usize,
    /// Set while a frame's bytes are being written and left set if writing
    /// fails, so a stream that may end in a partial frame is never extended
    /// or reported as finished.
    partial: bool,
}

impl<W> fmt::Debug for XyzWriter<W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("XyzWriter")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("source_label", &self.source_label)
            .field("frame_count", &self.frame_count)
            .field("partial", &self.partial)
            .finish_non_exhaustive()
    }
}

impl<W: Write> XyzWriter<W> {
    /// A writer of `topology` frames to `writer`; `source_label` names the
    /// destination in errors. Nothing is written until the first frame.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the length
    /// unit is not a length or the comment contains a line break,
    /// [`TrajectoryIoErrorKind::InvalidPrecision`] for more than 17 decimal
    /// places, and [`TrajectoryIoErrorKind::LimitExceeded`] for a comment over
    /// 1 MiB or a topology without atoms.
    pub fn new(
        writer: W,
        topology: impl Into<Arc<Topology>>,
        options: XyzWriteOptions,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryIoError> {
        let topology = topology.into();
        let source_label = source_label.into();
        options
            .length_unit
            .conversion_factor_to(CANONICAL_LENGTH_UNIT)
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Open,
                    Some(TrajectoryFormat::Xyz),
                    &source_label,
                    format!("XYZ length unit is incompatible: {error}"),
                )
            })?;
        if options.decimal_places > 17 {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidPrecision,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                &source_label,
                "XYZ decimal precision must be at most 17 places",
            ));
        }
        if options.comment.contains(['\r', '\n']) {
            return Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                &source_label,
                "XYZ writer comment cannot contain a line break",
            ));
        }
        if options.comment.len() > MAX_WRITER_COMMENT_BYTES {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xyz),
                &source_label,
                format!("XYZ writer comment exceeds the {MAX_WRITER_COMMENT_BYTES}-byte limit"),
            ));
        }
        if topology.atom_count() == 0 {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
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
            partial: false,
        })
    }

    /// The topology every written frame must share.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// The underlying stream.
    pub fn writer(&self) -> &W {
        &self.writer
    }

    /// Flushes the underlying stream.
    ///
    /// # Errors
    ///
    /// Returns the stream's flush error.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }

    pub(crate) fn validate_finish(&self) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::Finish,
            TrajectoryFormat::Xyz,
            &self.source_label,
        )?;
        require_nonempty_writer(self.frame_count, TrajectoryFormat::Xyz, &self.source_label)
    }

    /// Flushes and returns the completed nonempty XYZ stream.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidFrame`] when no frame was
    /// written or an earlier write left a partial frame, and
    /// [`TrajectoryIoErrorKind::Io`] when flushing fails.
    pub fn finish(mut self) -> Result<W, TrajectoryIoError> {
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

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Xyz,
            &self.source_label,
        )?;
        if !self.topology().shares_layout(frame.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        let unsupported = [
            (frame.cell().is_some(), TrajectoryField::PeriodicCell),
            (frame.velocities().is_some(), TrajectoryField::Velocities),
            (frame.forces().is_some(), TrajectoryField::Forces),
            (frame.time().is_some(), TrajectoryField::Time),
            (frame.step().is_some(), TrajectoryField::Step),
            (
                frame.properties().atoms().has_data(),
                TrajectoryField::AtomProperties,
            ),
            (
                frame.properties().bonds().has_data(),
                TrajectoryField::BondProperties,
            ),
            (
                !frame.properties().owner().is_empty(),
                TrajectoryField::Properties,
            ),
        ];
        if let Some((_, field)) = unsupported.into_iter().find(|(present, _)| *present) {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                field,
            ));
        }
        let factor = CANONICAL_LENGTH_UNIT
            .conversion_factor_to(self.options.length_unit)
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Xyz),
                    &self.source_label,
                    format!("XYZ output unit is incompatible: {error}"),
                )
            })?;
        let positions = frame.positions().values();
        for point in positions.value().iter() {
            if !Point3::new(point.x * factor, point.y * factor, point.z * factor).is_finite() {
                return Err(codec_context(
                    TrajectoryIoErrorKind::InvalidFrame,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Xyz),
                    &self.source_label,
                    "XYZ output coordinate is not finite after conversion",
                ));
            }
        }

        self.partial = true;
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
            let element = self
                .topology
                .atom(atom_id)
                .expect("dense topology atom IDs name topology atoms");
            writeln!(
                self.writer,
                "{} {x:.precision$} {y:.precision$} {z:.precision$}",
                element.element().symbol(),
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
        self.partial = false;
        self.frame_count = self.frame_count.checked_add(1).ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
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
) -> Result<(), TrajectoryIoError> {
    if atom_count > limits.max_atoms() {
        return Err(codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Xyz),
            source_label,
            format!(
                "XYZ atom count {atom_count} exceeds configured maximum {}",
                limits.max_atoms()
            ),
        )
        .with_limit(LimitExceeded::new(
            Resource::Atoms,
            limits.max_atoms(),
            atom_count,
        )));
    }
    let scratch_bytes = atom_count
        .checked_mul(std::mem::size_of::<Point3>())
        .ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xyz),
                source_label,
                "XYZ position scratch size overflows",
            )
        })?;
    if scratch_bytes > limits.max_scratch_bytes() {
        return Err(codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Xyz),
            source_label,
            format!(
                "XYZ position scratch requires {scratch_bytes} bytes, exceeding configured maximum {}",
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

fn validate_frame_bytes(
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
            format!(
                "XYZ frame exceeds configured {}-byte limit",
                limits.max_frame_bytes()
            ),
        )
        .with_limit(LimitExceeded::new(
            Resource::RecordBytes,
            usize::try_from(limits.max_frame_bytes()).unwrap_or(usize::MAX),
            usize::try_from(frame_bytes).unwrap_or(usize::MAX),
        )));
    }
    Ok(())
}

/// Adds the bytes one line consumed from the source, terminator included.
fn checked_frame_add(
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
                "XYZ frame byte count overflows",
            )
        })
}

/// Reads one line and returns it without its LF or CRLF terminator, together
/// with the exact number of source bytes consumed, terminator included.
///
/// The line, terminator included, must fit `limit`. The retained `line` buffer
/// grows only while its capacity plus `other_scratch_bytes` stays within
/// `max_scratch_bytes`, so whether a line fits depends only on its length.
#[allow(clippy::too_many_arguments)]
fn read_line<'a, R: BufRead>(
    reader: &mut R,
    line: &'a mut Vec<u8>,
    limit: usize,
    other_scratch_bytes: Option<usize>,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    frame: usize,
    operation: TrajectoryIoOperation,
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
                format!("XYZ line exceeds configured {limit}-byte limit"),
            )
            .with_limit(LimitExceeded::new(Resource::RecordBytes, limit, needed)));
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
                        format!("XYZ line buffer: {detail}"),
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
    let mut content = line.as_slice();
    if let Some(stripped) = content.strip_suffix(b"\n") {
        content = stripped.strip_suffix(b"\r").unwrap_or(stripped);
    }
    let content = std::str::from_utf8(content).map_err(|_| {
        io_context(
            operation,
            Some(TrajectoryFormat::Xyz),
            source_label,
            io::Error::new(
                io::ErrorKind::InvalidData,
                "stream did not contain valid UTF-8",
            ),
        )
    })?;
    Ok(Some((content, read)))
}

fn frame_codec_error(
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

fn writer_field_error(
    source_label: &str,
    frame: usize,
    field: TrajectoryField,
) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        TrajectoryIoErrorKind::UnsupportedField,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Xyz),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_field(field)
    .with_detail(format!("XYZ cannot preserve {field}"))
}
