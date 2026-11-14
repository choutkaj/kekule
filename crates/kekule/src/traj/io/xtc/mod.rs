//! Bounded GROMACS XTC trajectories with a checked reader and writer.
//!
//! [`XtcReader`] and [`IndexedXtcReader`] read GROMACS magic 1995 and 2023
//! frames; [`XtcWriter`] writes them at an explicit lossy precision. XTC
//! stores positions, a box, a step and a time in picoseconds; frames of more
//! than nine atoms are compressed with nominal resolution `1 / precision`
//! nanometres, smaller frames are stored as plain f32. The box follows an
//! [`XtcCellPolicy`]. A frame without a time or step can be numbered under
//! [`XtcWriteOptions::frame_timing`]; the velocities, forces, occupancies,
//! B-factors, and properties XTC cannot carry are rejected unless
//! [`XtcWriteOptions::geometry_metadata`] and
//! [`XtcWriteOptions::properties`] omit them, as [`XtcWriter::report`] then
//! reports.
//!
//! The scratch limit bounds the combined capacities of the encoded bytes,
//! the decoded scalar coordinates, and the dense positions.
//!
//! # Examples
//!
//! ```
//! use std::io::Cursor;
//! use std::sync::Arc;
//!
//! use kekule::math::Point3;
//! use kekule::units::{Quantity, NANOMETER, PICOSECOND};
//! use kekule::{smiles, topology::Topology};
//! use kekule::traj::io::xtc::{
//!     XtcCellPolicy, XtcReadOptions, XtcReader, XtcWriteOptions, XtcWriter,
//! };
//! use kekule::traj::{FrameBuffer, TrajectoryReader, TrajectoryWriter};
//!
//! let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//! let mut frame = FrameBuffer::new(Arc::clone(&topology));
//! let points = vec![Point3::new(0.25, 0.5, 0.0)];
//! frame.geometry_mut().set_positions(Quantity::new(points, NANOMETER))?;
//! frame.set_step(Some(100));
//! frame.set_time(Some(Quantity::new(0.5, PICOSECOND)))?;
//!
//! // A frame without a cell is written as the GROMACS zero box.
//! let options = XtcWriteOptions::default().with_cell_policy(XtcCellPolicy::ZeroMatrixAsAbsent);
//! let mut writer = XtcWriter::new(Vec::new(), Arc::clone(&topology), options)?;
//! writer.write_frame(frame.frame_view())?;
//! let bytes = writer.finish()?;
//!
//! let read_options =
//!     XtcReadOptions::default().with_cell_policy(XtcCellPolicy::ZeroMatrixAsAbsent);
//! let mut reader = XtcReader::new(Cursor::new(bytes), topology, read_options)?;
//! let mut buffer = reader.frame_buffer();
//! assert!(reader.read_next(&mut buffer)?);
//! assert_eq!(buffer.step(), Some(100));
//! assert_eq!(buffer.time().unwrap().value_in(PICOSECOND)?, 0.5);
//! assert!(buffer.cell().is_none());
//! let position = buffer.positions().get(0).unwrap().value_in(NANOMETER)?;
//! assert_eq!(position, Point3::new(0.25, 0.5, 0.0));
//! assert!(!reader.read_next(&mut buffer)?);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::traj::trajectory::reason::{Datum, Part, Reason};
use crate::traj::FrameView;
use std::fmt;
use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::Arc;

use crate::diagnostics::{LimitExceeded, Resource, UnrepresentablePolicy, WriteReport};
use crate::io::write_policy::StatePolicies;
use crate::math::{PeriodicCell, Point3, Vector3};
use crate::topology::Topology;
use crate::traj::io::DecodedGeometry;
use crate::traj::{
    FrameBuffer, FrameData, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use crate::units::{Quantity, NANOMETER, PER_NANOMETER, PICOSECOND};

use super::codec::{
    frame_error, header_error, limit_error, read_exact, stream_error, validate_atom_count,
};
use super::frame_state::{check_frame_state, frame_interval, record};
use super::indexed_reader::{indexed_reader_traits, IndexedCodec, IndexedReader, IndexedSource};
use super::limits::{exceeded_bytes, reserve_scratch, scratch_bytes};
use super::metadata::{FieldAvailability, TrajectoryFieldAvailability};
use super::options::ReadSource;
use super::TrajectoryWriteOptions;
use super::{
    codec_context, frame_offset_context, io_context, probe_seekable_eof, require_no_partial_frame,
    require_nonempty_writer, validation_cause, FrameOffsets, FrameStamp, FrameTimingPolicy,
    IndexBuildControl, TrajectoryIoLimits, TrajectoryOpenOptions,
};

/// The format of this codec's streams.
const FORMAT: TrajectoryFormat = TrajectoryFormat::Xtc;

mod indexed;
mod read;
mod write;

pub use indexed::IndexedXtcReader;
pub use read::XtcReader;
use read::*;
pub(crate) use write::check_append_frames;
pub use write::XtcWriter;

/// What an XTC frame carries: a cell, which a zero box makes absent
/// under [`XtcCellPolicy::ZeroMatrixAsAbsent`].
pub(super) const FIELDS: TrajectoryFieldAvailability = TrajectoryFieldAvailability {
    positions: FieldAvailability::Required,
    cell: FieldAvailability::Required,
    velocities: FieldAvailability::Absent,
    forces: FieldAvailability::Absent,
    time: FieldAvailability::Required,
    step: FieldAvailability::Required,
    properties: FieldAvailability::Absent,
};

const XTC_HEADER_BYTES: usize = 56;
const XTC_COMPRESSED_PRELUDE_BYTES: usize = 32;
const XTC_MAGIC_INTS: [i32; 73] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64, 80, 101, 128, 161, 203, 256,
    322, 406, 512, 645, 812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321,
    13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031, 131072, 165140, 208063,
    262144, 330280, 416127, 524287, 660561, 832255, 1048576, 1321122, 1664510, 2097152, 2642245,
    3329021, 4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216,
];
const XTC_FIRST_SMALL_INDEX: usize = 9;
/// Largest scaled coordinate magnitude, and bound on each axis's span, that
/// GROMACS writes (its `MAXABS`). Wider spans overflow the encoder's sizes.
const XTC_MAX_SCALED: i64 = i32::MAX as i64 - 2;
/// Most atoms a 1995-magic frame holds before its signed 32-bit compressed
/// byte count can overflow (GROMACS `xtc_1995_max_natoms`).
const XTC_1995_MAX_ATOMS: usize = 298_261_617;

mod encode;

/// Supported XTC magic and compressed-byte-count profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum XtcMagic {
    /// Magic 1995 with a 32-bit compressed byte count, readable by every
    /// GROMACS version.
    Xtc1995,
    /// Magic 2023 with a 64-bit compressed byte count, for frames too large
    /// for 1995.
    Xtc2023,
}

impl XtcMagic {
    const fn nbytes_width(self) -> usize {
        match self {
            Self::Xtc1995 => 4,
            Self::Xtc2023 => 8,
        }
    }

    const fn number(self) -> i32 {
        match self {
            Self::Xtc1995 => 1995,
            Self::Xtc2023 => 2023,
        }
    }
}

/// Explicit handling for the mandatory XTC box matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum XtcCellPolicy {
    /// Require a finite, nondegenerate, fully periodic cell.
    RequirePeriodic,
    /// Treat an exactly zero matrix as explicitly absent; validate all others.
    ZeroMatrixAsAbsent,
}

/// XTC reader policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct XtcReadOptions {
    cell_policy: XtcCellPolicy,
}

impl Default for XtcReadOptions {
    fn default() -> Self {
        Self {
            cell_policy: XtcCellPolicy::RequirePeriodic,
        }
    }
}

impl XtcReadOptions {
    /// Sets [`Self::cell_policy`].
    #[must_use]
    pub const fn with_cell_policy(mut self, cell_policy: XtcCellPolicy) -> Self {
        self.cell_policy = cell_policy;
        self
    }

    /// How the box matrix is interpreted; [`XtcCellPolicy::RequirePeriodic`]
    /// by default.
    pub const fn cell_policy(&self) -> XtcCellPolicy {
        self.cell_policy
    }
}

/// XTC writer policy with explicit lossy coordinate precision.
///
/// Box components always use f32. Cells that become degenerate at that precision
/// are rejected before writing the frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XtcWriteOptions {
    magic: XtcMagic,
    precision: f32,
    cell_policy: XtcCellPolicy,
    geometry_metadata: UnrepresentablePolicy,
    properties: UnrepresentablePolicy,
    frame_timing: FrameTimingPolicy,
    frame_interval: Quantity<f64>,
}

impl Default for XtcWriteOptions {
    fn default() -> Self {
        Self {
            magic: XtcMagic::Xtc1995,
            precision: 1000.0,
            cell_policy: XtcCellPolicy::RequirePeriodic,
            geometry_metadata: UnrepresentablePolicy::Fail,
            properties: UnrepresentablePolicy::Fail,
            frame_timing: FrameTimingPolicy::Require,
            frame_interval: Quantity::new(1.0, PICOSECOND),
        }
    }
}

impl XtcWriteOptions {
    /// Sets [`Self::magic`].
    #[must_use]
    pub const fn with_magic(mut self, magic: XtcMagic) -> Self {
        self.magic = magic;
        self
    }

    /// Sets the coordinate scale, an inverse length such as `1000 1/nm`.
    ///
    /// Compressed coordinates have nominal resolution `1 / precision`. The
    /// scale is stored as an `f32` in inverse nanometres, as XTC records it.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidPrecision`] when `precision` is
    /// not an inverse length or is not finite and positive as an `f32`.
    pub fn with_precision(mut self, precision: Quantity<f64>) -> Result<Self, TrajectoryIoError> {
        let invalid = || {
            codec_context(
                TrajectoryIoErrorKind::InvalidPrecision,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                "XTC writer options",
                Reason::PrecisionNotInverseLength,
            )
        };
        let precision = precision.value_in(PER_NANOMETER).map_err(|_| invalid())? as f32;
        validate_precision(precision, "XTC writer options")?;
        self.precision = precision;
        Ok(self)
    }

    /// Chooses how frames without a cell are written, mirroring the reader's
    /// [`XtcCellPolicy`]. [`XtcCellPolicy::ZeroMatrixAsAbsent`] writes the
    /// zero box that GROMACS uses for non-periodic systems; the default,
    /// [`XtcCellPolicy::RequirePeriodic`], rejects such frames.
    #[must_use]
    pub const fn with_cell_policy(mut self, cell_policy: XtcCellPolicy) -> Self {
        self.cell_policy = cell_policy;
        self
    }

    /// Sets [`Self::geometry_metadata`].
    #[must_use]
    pub const fn with_geometry_metadata(mut self, policy: UnrepresentablePolicy) -> Self {
        self.geometry_metadata = policy;
        self
    }

    /// Sets [`Self::properties`].
    #[must_use]
    pub const fn with_properties(mut self, policy: UnrepresentablePolicy) -> Self {
        self.properties = policy;
        self
    }

    /// Sets [`Self::frame_timing`].
    #[must_use]
    pub const fn with_frame_timing(mut self, frame_timing: FrameTimingPolicy) -> Self {
        self.frame_timing = frame_timing;
        self
    }

    /// Sets [`Self::frame_interval`].
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when
    /// `interval` is not a finite, nonnegative time.
    pub fn with_frame_interval(
        mut self,
        interval: Quantity<f64>,
    ) -> Result<Self, TrajectoryIoError> {
        self.frame_interval = frame_interval(TrajectoryFormat::Xtc, interval)?;
        Ok(self)
    }

    /// The magic written; [`XtcMagic::Xtc1995`] by default.
    pub const fn magic(self) -> XtcMagic {
        self.magic
    }

    /// The coordinate scale, in inverse nanometres; 1000 by default.
    pub fn precision(self) -> Quantity<f64> {
        Quantity::new(f64::from(self.precision), PER_NANOMETER)
    }

    /// How frames without a cell are written;
    /// [`XtcCellPolicy::RequirePeriodic`] by default.
    pub const fn cell_policy(self) -> XtcCellPolicy {
        self.cell_policy
    }

    /// A frame's velocities, forces, occupancies, and B-factors, which XTC
    /// cannot carry; by default a frame with any of them fails to write.
    pub const fn geometry_metadata(self) -> UnrepresentablePolicy {
        self.geometry_metadata
    }

    /// A frame's atom, bond, and own properties, which XTC cannot carry; by
    /// default a frame with any of them fails to write.
    pub const fn properties(self) -> UnrepresentablePolicy {
        self.properties
    }

    /// The time and step of a frame without them, which XTC records for
    /// every frame; [`FrameTimingPolicy::Require`] by default.
    pub const fn frame_timing(self) -> FrameTimingPolicy {
        self.frame_timing
    }

    /// The time between frames [`FrameTimingPolicy::Number`] gives; 1 ps by
    /// default.
    pub const fn frame_interval(self) -> Quantity<f64> {
        self.frame_interval
    }

    const fn state_policies(self) -> StatePolicies {
        StatePolicies {
            carries: FIELDS.carries(),
            geometry_metadata: self.geometry_metadata,
            properties: self.properties,
            cell: UnrepresentablePolicy::Fail,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct XtcFrameInfo {
    start: u64,
    end: u64,
    magic: XtcMagic,
    atom_count: usize,
    step: u64,
    time: f64,
    cell: Option<PeriodicCell>,
    precision: Option<f32>,
    compressed_bytes: usize,
}

impl XtcFrameInfo {
    pub(crate) const fn atom_count(&self) -> usize {
        self.atom_count
    }

    pub(crate) const fn magic(&self) -> XtcMagic {
        self.magic
    }

    pub(crate) const fn precision(&self) -> Option<f32> {
        self.precision
    }

    /// The time, in picoseconds, and step the frame records.
    const fn stamp(&self) -> FrameStamp {
        FrameStamp {
            time: Some(self.time),
            step: Some(self.step),
        }
    }

    pub(crate) const fn compressed_bytes(&self) -> usize {
        self.compressed_bytes
    }
}

fn nonnegative_xdr_count(
    value: i32,
    field: Part,
    source_label: &str,
    frame: usize,
    byte_offset: u64,
) -> Result<usize, TrajectoryIoError> {
    usize::try_from(value).map_err(|_| {
        TrajectoryIoError::codec(
            TrajectoryIoErrorKind::InvalidHeader,
            TrajectoryIoOperation::ReadHeader,
            Some(FORMAT),
        )
        .with_source_label(source_label)
        .with_frame(frame)
        .with_byte_offset(byte_offset)
        .with_reason(Reason::Negative(field))
    })
}

fn validate_precision(precision: f32, source_label: &str) -> Result<(), TrajectoryIoError> {
    if !precision.is_finite() || precision <= 0.0 {
        return Err(codec_context(
            TrajectoryIoErrorKind::InvalidPrecision,
            TrajectoryIoOperation::Open,
            Some(FORMAT),
            source_label,
            Reason::PrecisionNotPositive,
        ));
    }
    Ok(())
}

fn finite_f32(value: f64, source_label: &str, field: Datum) -> Result<f32, TrajectoryIoError> {
    let narrowed = value as f32;
    if !value.is_finite() || !narrowed.is_finite() {
        return Err(codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(FORMAT),
            source_label,
            Reason::DoesNotFitSingle(Part::ValueOf(field)),
        ));
    }
    Ok(narrowed)
}
