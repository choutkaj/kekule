//! Strict common-profile CHARMM/NAMD/OpenMM DCD trajectory I/O.
//!
//! [`DcdReader`] and [`IndexedDcdReader`] read 32-bit-record `CORD` files in
//! either byte order, with optional unit cells and fixed atoms. They read
//! every complete frame in the file, as VMD and MDAnalysis do, whatever
//! frame count `NSET` the header declares; an indexed reader's open report
//! notes a count that differs ([`TrajectoryOpenNote::DcdFrameCount`]).
//! [`DcdWriter`]
//! writes canonical all-atom files. DCD stores f32 coordinates in angstrom and
//! one step per frame on the `ISTART`/`NSAVC` sequence; time is derived from
//! the header `DELTA` only under an explicit [`DcdTimePolicy`]. The writer
//! rejects frame state DCD cannot carry (velocities, forces, occupancies,
//! B-factors, times under [`DcdTimePolicy::Absent`], cells without
//! [`DcdWriteOptions::unit_cells`], and properties) unless
//! [`DcdWriteOptions::geometry_metadata`], [`DcdWriteOptions::cell`], and
//! [`DcdWriteOptions::properties`] omit it, and reports what it leaves out
//! in [`DcdWriter::report`]. Under [`FrameTimingPolicy::Number`] it numbers
//! frames whose steps or times are missing or off the header's sequence, and
//! reports the steps and times it replaced.
//!
//! The scratch limit bounds the combined capacities of the coordinates, the
//! record bytes, and the fixed-atom buffers.
//!
//! [`TrajectoryOpenNote::DcdFrameCount`]: super::TrajectoryOpenNote::DcdFrameCount
//!
//! # Examples
//!
//! ```
//! use std::io::Cursor;
//! use std::sync::Arc;
//!
//! use kekule::math::Point3;
//! use kekule::units::{Quantity, ANGSTROM};
//! use kekule::{smiles, topology::Topology};
//! use kekule::traj::io::dcd::{DcdReadOptions, DcdReader, DcdWriteOptions, DcdWriter};
//! use kekule::traj::{FrameBuffer, TrajectoryReader, TrajectoryWriter};
//!
//! let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//!
//! // Frame steps must follow the writer's step sequence exactly.
//! let options = DcdWriteOptions::default().with_step_sequence(10, 10)?;
//! let stream = Cursor::new(Vec::new());
//! let mut writer = DcdWriter::new(stream, Arc::clone(&topology), options)?;
//! let mut frame = FrameBuffer::new(Arc::clone(&topology));
//! for step in [10, 20] {
//!     let x = step as f64;
//!     let points = vec![Point3::new(x, 0.0, 0.0)];
//!     frame.geometry_mut().set_positions(Quantity::new(points, ANGSTROM))?;
//!     frame.set_step(Some(step));
//!     writer.write_frame(frame.frame_view())?;
//! }
//! let mut stream = writer.finish()?;
//! stream.set_position(0);
//!
//! let mut reader = DcdReader::new(stream, topology, DcdReadOptions::default())?;
//! let mut frame = reader.frame_buffer();
//! assert!(reader.read_next(&mut frame)?);
//! assert_eq!(frame.step(), Some(10));
//! assert_eq!(frame.time(), None);
//! assert!(reader.read_next(&mut frame)?);
//! assert_eq!(frame.step(), Some(20));
//! // DCD stores f32 angstrom; compare at that precision.
//! let position = frame.positions().get(0).unwrap().value_in(ANGSTROM)?;
//! assert_eq!([position.x as f32, position.y as f32, position.z as f32], [20.0, 0.0, 0.0]);
//! assert!(!reader.read_next(&mut frame)?);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::traj::trajectory::reason::{DcdField, Part, Reason, UnitOf};
use crate::traj::FrameView;
use std::fmt;
use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::Arc;

use crate::diagnostics::{
    GeometryField, LimitExceeded, Resource, UnrepresentablePolicy, WriteOmission, WriteReport,
};
use crate::io::write_policy::{Carries, StatePolicies};
use crate::math::{PeriodicCell, Point3, Vector3};
use crate::topology::Topology;
use crate::traj::io::DecodedGeometry;
use crate::traj::{
    FrameBuffer, FrameData, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use crate::units::{Quantity, Unit, ANGSTROM, CANONICAL_LENGTH_UNIT, CANONICAL_TIME_UNIT};

use super::cell::standard_orientation;
use super::codec::{
    frame_error, header_error, limit_error, read_exact, stream_error, validate_atom_count,
    writer_limit,
};
use super::frame_state::{check_frame_state, record};
use super::indexed_reader::{
    indexed_reader_traits, IndexFacts, IndexedCodec, IndexedReader, IndexedSource,
};
use super::limits::{exceeded_bytes, filled_bytes, reserve_scratch, scratch_bytes};
use super::metadata::{FieldAvailability, TrajectoryFieldAvailability};
use super::options::ReadSource;
use super::TrajectoryWriteOptions;
use super::{
    codec_context, frame_offset_context, io_context, probe_seekable_eof, require_no_partial_frame,
    require_nonempty_writer, validation_cause, FrameOffsets, FrameStamp, FrameTimingPolicy,
    IndexBuildControl, TrajectoryIoLimits, TrajectoryOpenOptions,
};

/// The format of this codec's streams.
const FORMAT: TrajectoryFormat = TrajectoryFormat::Dcd;

mod indexed;
mod read;
mod write;

pub use indexed::IndexedDcdReader;
pub use read::DcdReader;
use read::*;
pub(crate) use write::check_append_header;
pub use write::DcdWriter;

/// What a DCD frame carries: a cell when the header's unit-cell flag is
/// set, and a time under [`DcdTimePolicy::HeaderDelta`].
pub(super) const FIELDS: TrajectoryFieldAvailability = TrajectoryFieldAvailability {
    positions: FieldAvailability::Required,
    cell: FieldAvailability::Optional,
    velocities: FieldAvailability::Absent,
    forces: FieldAvailability::Absent,
    time: FieldAvailability::Optional,
    step: FieldAvailability::Required,
    properties: FieldAvailability::Absent,
};

const HEADER_BYTES: usize = 84;
const CELL_BYTES: usize = 48;

/// Byte order used by DCD record markers and numeric payloads.
///
/// This enum is closed: a DCD file is either little- or big-endian.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DcdEndian {
    /// Little-endian, the common default.
    Little,
    /// Big-endian.
    Big,
}

impl DcdEndian {
    fn u32(self, bytes: [u8; 4]) -> u32 {
        match self {
            Self::Little => u32::from_le_bytes(bytes),
            Self::Big => u32::from_be_bytes(bytes),
        }
    }

    fn i32(self, bytes: [u8; 4]) -> i32 {
        match self {
            Self::Little => i32::from_le_bytes(bytes),
            Self::Big => i32::from_be_bytes(bytes),
        }
    }

    fn f32(self, bytes: [u8; 4]) -> f32 {
        f32::from_bits(self.u32(bytes))
    }

    fn f64(self, bytes: [u8; 8]) -> f64 {
        match self {
            Self::Little => f64::from_le_bytes(bytes),
            Self::Big => f64::from_be_bytes(bytes),
        }
    }

    fn encode_u32(self, value: u32) -> [u8; 4] {
        match self {
            Self::Little => value.to_le_bytes(),
            Self::Big => value.to_be_bytes(),
        }
    }

    fn encode_i32(self, value: i32) -> [u8; 4] {
        match self {
            Self::Little => value.to_le_bytes(),
            Self::Big => value.to_be_bytes(),
        }
    }

    fn encode_f32(self, value: f32) -> [u8; 4] {
        self.encode_u32(value.to_bits())
    }

    fn encode_f64(self, value: f64) -> [u8; 8] {
        match self {
            Self::Little => value.to_le_bytes(),
            Self::Big => value.to_be_bytes(),
        }
    }
}

/// Explicit interpretation of the DCD header's unitless `DELTA` value.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum DcdTimePolicy {
    /// Preserve steps but do not synthesize per-frame time.
    Absent,
    /// Treat `DELTA` as the duration of one integration step in this unit.
    HeaderDelta {
        /// The time unit of `DELTA`.
        unit: Unit,
    },
}

/// DCD reader policy.
#[derive(Debug, Clone, PartialEq)]
pub struct DcdReadOptions {
    time_policy: DcdTimePolicy,
}

impl Default for DcdReadOptions {
    fn default() -> Self {
        Self {
            time_policy: DcdTimePolicy::Absent,
        }
    }
}

impl DcdReadOptions {
    /// Sets [`Self::time_policy`].
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the
    /// policy's unit is not a time.
    pub fn with_time_policy(
        mut self,
        time_policy: DcdTimePolicy,
    ) -> Result<Self, TrajectoryIoError> {
        validate_time_policy(time_policy, "DCD reader options")?;
        self.time_policy = time_policy;
        Ok(self)
    }

    /// How frame times are derived; [`DcdTimePolicy::Absent`] by default.
    pub const fn time_policy(&self) -> DcdTimePolicy {
        self.time_policy
    }
}

/// Canonical DCD writer policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DcdWriteOptions {
    endian: DcdEndian,
    unit_cells: bool,
    start_step: u64,
    step_interval: u64,
    delta: f32,
    time_policy: DcdTimePolicy,
    frame_timing: FrameTimingPolicy,
    geometry_metadata: UnrepresentablePolicy,
    properties: UnrepresentablePolicy,
    cell: UnrepresentablePolicy,
}

impl Default for DcdWriteOptions {
    fn default() -> Self {
        Self {
            endian: DcdEndian::Little,
            unit_cells: false,
            start_step: 0,
            step_interval: 1,
            delta: 1.0,
            time_policy: DcdTimePolicy::Absent,
            frame_timing: FrameTimingPolicy::Require,
            geometry_metadata: UnrepresentablePolicy::Fail,
            properties: UnrepresentablePolicy::Fail,
            cell: UnrepresentablePolicy::Fail,
        }
    }
}

impl DcdWriteOptions {
    /// Sets [`Self::endian`].
    #[must_use]
    pub const fn with_endian(mut self, endian: DcdEndian) -> Self {
        self.endian = endian;
        self
    }

    /// Sets [`Self::unit_cells`].
    #[must_use]
    pub const fn with_unit_cells(mut self, unit_cells: bool) -> Self {
        self.unit_cells = unit_cells;
        self
    }

    /// Sets [`Self::cell`].
    #[must_use]
    pub const fn with_cell(mut self, cell: UnrepresentablePolicy) -> Self {
        self.cell = cell;
        self
    }

    /// Sets the header step sequence: frame `n` must carry step
    /// `start_step + n * step_interval`, or is numbered with it under
    /// [`FrameTimingPolicy::Number`].
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the
    /// interval is zero or either exceeds the header's `i32`.
    pub fn with_step_sequence(
        mut self,
        start_step: u64,
        step_interval: u64,
    ) -> Result<Self, TrajectoryIoError> {
        if step_interval == 0 {
            return Err(writer_option_error(Reason::StepIntervalZero));
        }
        if start_step > i32::MAX as u64 || step_interval > i32::MAX as u64 {
            return Err(writer_option_error(Reason::StepSequenceRange));
        }
        self.start_step = start_step;
        self.step_interval = step_interval;
        Ok(self)
    }

    /// Stores `DELTA`, the duration of one step, as an `f32` in the time
    /// unit of [`DcdTimePolicy::HeaderDelta`], or in the canonical time unit
    /// for [`DcdTimePolicy::Absent`]. With `HeaderDelta`, every frame time
    /// must equal its step times `DELTA` to f32 precision, the precision
    /// `DELTA` has in the file.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when `delta`
    /// is not a time, or not finite and nonnegative as an `f32`, or the
    /// policy's unit is not a time.
    pub fn with_header_delta(
        mut self,
        delta: Quantity<f64>,
        time_policy: DcdTimePolicy,
    ) -> Result<Self, TrajectoryIoError> {
        validate_time_policy(time_policy, "DCD writer options")?;
        let unit = header_delta_unit(time_policy);
        let delta = delta
            .value_in(unit)
            .map_err(|error| writer_option_error(Reason::DeltaNotTime).with_source(error))?
            as f32;
        if !delta.is_finite() || delta < 0.0 {
            return Err(writer_option_error(Reason::InvalidDelta(delta)));
        }
        self.delta = delta;
        self.time_policy = time_policy;
        Ok(self)
    }

    /// Sets [`Self::frame_timing`].
    #[must_use]
    pub const fn with_frame_timing(mut self, frame_timing: FrameTimingPolicy) -> Self {
        self.frame_timing = frame_timing;
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

    /// A frame's velocities, forces, occupancies, and B-factors, and its
    /// time under [`DcdTimePolicy::Absent`], which DCD cannot carry; by
    /// default a frame with any of them fails to write.
    pub const fn geometry_metadata(&self) -> UnrepresentablePolicy {
        self.geometry_metadata
    }

    /// A frame's atom, bond, and own properties, which DCD cannot carry; by
    /// default a frame with any of them fails to write.
    pub const fn properties(&self) -> UnrepresentablePolicy {
        self.properties
    }

    /// The periodic cell of a frame written without [`Self::unit_cells`];
    /// by default such a frame fails to write.
    pub const fn cell(&self) -> UnrepresentablePolicy {
        self.cell
    }

    /// How a frame whose step or time is missing or off the header's
    /// sequence is written; [`FrameTimingPolicy::Require`], under which
    /// such a frame fails, by default.
    pub const fn frame_timing(&self) -> FrameTimingPolicy {
        self.frame_timing
    }

    /// The byte order written; [`DcdEndian::Little`] by default.
    pub const fn endian(&self) -> DcdEndian {
        self.endian
    }

    /// Whether the header's unit-cell flag is set, so that every frame
    /// carries a cell; `false` by default. Cells are written in the
    /// canonical DCD orientation: `a` along positive x, `b` in the xy plane
    /// with positive y, and `c` with positive z. Other orientations are
    /// rejected because DCD lengths and angles cannot preserve their
    /// relation to the stored Cartesian coordinates. Angles are written
    /// using the unambiguous supported cosine convention.
    pub const fn unit_cells(&self) -> bool {
        self.unit_cells
    }

    /// The step of the first frame (`ISTART`); 0 by default.
    pub const fn start_step(&self) -> u64 {
        self.start_step
    }

    /// The step increment between frames (`NSAVC`); 1 by default.
    pub const fn step_interval(&self) -> u64 {
        self.step_interval
    }

    /// `DELTA` in the unit [`Self::with_header_delta`] stores it in.
    pub fn delta(&self) -> Quantity<f64> {
        Quantity::new(f64::from(self.delta), header_delta_unit(self.time_policy))
    }

    /// How frame times relate to `DELTA`; [`DcdTimePolicy::Absent`] by
    /// default, under which a frame's time follows
    /// [`Self::geometry_metadata`].
    pub const fn time_policy(&self) -> DcdTimePolicy {
        self.time_policy
    }

    /// The policies of the shared check: a frame's cell is carried with
    /// unit cells, and its time under [`DcdTimePolicy::HeaderDelta`].
    const fn state_policies(&self) -> StatePolicies {
        StatePolicies {
            carries: Carries {
                cell: self.unit_cells,
                time: matches!(self.time_policy, DcdTimePolicy::HeaderDelta { .. }),
                ..FIELDS.carries()
            },
            geometry_metadata: self.geometry_metadata,
            properties: self.properties,
            cell: self.cell,
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct DcdHeader {
    endian: DcdEndian,
    declared_frames: usize,
    start_step: u64,
    step_interval: u64,
    delta: f32,
    fixed_count: usize,
    atom_count: usize,
    has_cell: bool,
    charmm_version: i32,
    /// Where the header record starts.
    header_start: u64,
    data_start: u64,
}

impl DcdHeader {
    pub(crate) const fn declared_frames(&self) -> usize {
        self.declared_frames
    }

    pub(crate) const fn atom_count(&self) -> usize {
        self.atom_count
    }

    pub(crate) const fn has_cell(&self) -> bool {
        self.has_cell
    }

    pub(crate) fn variant(&self) -> String {
        format!(
            "CHARMM-compatible v{} {:?}{}{}",
            self.charmm_version,
            self.endian,
            if self.has_cell { " cell" } else { "" },
            if self.fixed_count > 0 {
                " fixed-atoms"
            } else {
                ""
            }
        )
    }
}

/// The unit `DELTA` is written in under a time policy.
const fn header_delta_unit(policy: DcdTimePolicy) -> Unit {
    match policy {
        DcdTimePolicy::HeaderDelta { unit } => unit,
        DcdTimePolicy::Absent => CANONICAL_TIME_UNIT,
    }
}

fn validate_time_policy(
    policy: DcdTimePolicy,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    if let DcdTimePolicy::HeaderDelta { unit } = policy {
        unit.conversion_factor_to(CANONICAL_TIME_UNIT)
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Open,
                    Some(FORMAT),
                    source_label,
                    Reason::IncompatibleUnit(Some(FORMAT), UnitOf::TimePolicy),
                )
                .with_source(error)
            })?;
    }
    Ok(())
}

/// An invalid DCD writer option.
fn writer_option_error(reason: impl Into<Option<Reason>>) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::InconsistentMetadata,
        TrajectoryIoOperation::Open,
        Some(FORMAT),
        "DCD writer options",
        reason,
    )
}

fn checked_coordinate_bytes(count: usize, source_label: &str) -> Result<usize, TrajectoryIoError> {
    count.checked_mul(4).ok_or_else(|| {
        limit_error(
            FORMAT,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            None,
            Reason::AddressRange(Part::CoordinateRecordSize),
        )
    })
}

fn nonnegative_u64(
    value: i32,
    field: DcdField,
    source_label: &str,
) -> Result<u64, TrajectoryIoError> {
    u64::try_from(value).map_err(|_| {
        header_error(
            FORMAT,
            source_label,
            Reason::HeaderFieldNegative(field, value),
        )
    })
}

fn positive_u64(value: i32, field: DcdField, source_label: &str) -> Result<u64, TrajectoryIoError> {
    if value <= 0 {
        return Err(header_error(
            FORMAT,
            source_label,
            Reason::HeaderFieldNotPositive(field, value),
        ));
    }
    Ok(value as u64)
}

fn nonnegative_usize(
    value: i32,
    field: DcdField,
    source_label: &str,
) -> Result<usize, TrajectoryIoError> {
    usize::try_from(value).map_err(|_| {
        header_error(
            FORMAT,
            source_label,
            Reason::HeaderFieldNegative(field, value),
        )
    })
}
