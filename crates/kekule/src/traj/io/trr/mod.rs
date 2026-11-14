//! Pure-Rust GROMACS TRR/XDR trajectory I/O.
//!
//! [`TrrReader`] and [`IndexedTrrReader`] read big-endian XDR frames with f32
//! or f64 position, box, velocity and force blocks; [`TrrWriter`] writes them
//! at one explicit [`TrrScalarPrecision`]. Every frame carries a step and a
//! time; GROMACS nanometre, picosecond and kJ/(mol nm) units are converted at
//! the codec boundary. The lambda scalar follows a [`TrrLambdaPolicy`]. A
//! frame without a time or step can be numbered under
//! [`TrrWriteOptions::frame_timing`]; the occupancies, B-factors, and
//! properties TRR cannot carry are rejected unless
//! [`TrrWriteOptions::geometry_metadata`] and
//! [`TrrWriteOptions::properties`] omit them, as [`TrrWriter::report`] then
//! reports.
//!
//! The scratch limit bounds the three dense decode arrays, one vector per
//! atom each, and the reusable raw record buffer.
//!
//! # Examples
//!
//! ```
//! use std::io::Cursor;
//! use std::sync::Arc;
//!
//! use kekule::math::{Point3, Vector3};
//! use kekule::units::{Quantity, NANOMETER, NANOMETER_PER_PICOSECOND, PICOSECOND};
//! use kekule::{smiles, topology::Topology};
//! use kekule::traj::io::trr::{
//!     TrrLambdaPolicy, TrrReadOptions, TrrReader, TrrWriteOptions, TrrWriter,
//! };
//! use kekule::traj::{FrameBuffer, TrajectoryReader, TrajectoryWriter};
//!
//! let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//! let mut frame = FrameBuffer::new(Arc::clone(&topology));
//! let points = vec![Point3::new(0.25, 0.5, 0.0)];
//! frame.geometry_mut().set_positions(Quantity::new(points, NANOMETER))?;
//! let velocity = Vector3::new(1.0, 0.0, 0.0);
//! frame.set_velocities(Some(Quantity::new([velocity], NANOMETER_PER_PICOSECOND)))?;
//! frame.set_step(Some(500));
//! frame.set_time(Some(Quantity::new(1.0, PICOSECOND)))?;
//!
//! // Frames without a lambda property need the RequireZero policy.
//! let options = TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero);
//! let mut writer = TrrWriter::new(Vec::new(), Arc::clone(&topology), options)?;
//! writer.write_frame(frame.frame_view())?;
//! let bytes = writer.finish()?;
//!
//! let read_options = TrrReadOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero);
//! let mut reader = TrrReader::new(Cursor::new(bytes), topology, read_options)?;
//! let mut buffer = reader.frame_buffer();
//! assert!(reader.read_next(&mut buffer)?);
//! assert_eq!(buffer.step(), Some(500));
//! assert_eq!(buffer.time().unwrap().value_in(PICOSECOND)?, 1.0);
//! let position = buffer.positions().get(0).unwrap().value_in(NANOMETER)?;
//! assert_eq!(position, Point3::new(0.25, 0.5, 0.0));
//! let read_velocities = buffer.velocities().unwrap().clone().into_values();
//! let read_velocities = read_velocities.into_unit(NANOMETER_PER_PICOSECOND)?.into_value();
//! assert_eq!(read_velocities, [velocity]);
//! assert!(buffer.forces().is_none());
//! assert!(!reader.read_next(&mut buffer)?);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::traj::FrameView;
use std::fmt;
use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::Arc;

use crate::diagnostics::{
    GeometryField, LimitExceeded, Resource, UnrepresentablePolicy, WriteOmission, WriteReport,
};
use crate::io::write_policy::StatePolicies;
use crate::math::{PeriodicCell, Point3, Vector3};
use crate::properties::{GeometryProperties, OwnerProperties, PropertyKey, PropertyValue};
use crate::structure::DecodedGeometry;
use crate::topology::Topology;
use crate::traj::{
    FrameBuffer, FrameData, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use crate::units::{
    Quantity, CANONICAL_LENGTH_UNIT, DIMENSIONLESS, KILOJOULE_PER_MOLE_PER_NANOMETER, NANOMETER,
    NANOMETER_PER_PICOSECOND, PICOSECOND,
};

use super::codec::{
    frame_error, header_error, limit_error, read_exact, read_operation, validate_atom_count,
    writer_limit,
};
use super::frame_state::{check_frame_state, frame_interval, record, writer_field_error};
use super::indexed_reader::{
    indexed_reader_traits, IndexFacts, IndexedCodec, IndexedReader, IndexedSource,
};
use super::limits::{exceeded_bytes, filled_bytes};
use super::metadata::{FieldAvailability, TrajectoryFieldAvailability};
use super::options::ReadSource;
use super::TrajectoryWriteOptions;
use super::{
    codec_context, frame_offset_context, index_mismatch, io_context, probe_seekable_eof,
    require_no_partial_frame, require_nonempty_writer, validation_cause, FrameOffsets, FrameStamp,
    FrameTimingPolicy, IndexBuildControl, TrajectoryIoLimits, TrajectoryOpenOptions,
};

/// The format of this codec's streams.
const FORMAT: TrajectoryFormat = TrajectoryFormat::Trr;

mod indexed;
mod read;
mod write;

pub use indexed::IndexedTrrReader;
pub use read::TrrReader;
use read::*;
pub(crate) use write::check_append_frames;
pub use write::TrrWriter;

/// What a TRR frame carries; its only property is the lambda.
pub(super) const FIELDS: TrajectoryFieldAvailability = TrajectoryFieldAvailability {
    positions: FieldAvailability::Required,
    cell: FieldAvailability::Optional,
    velocities: FieldAvailability::Optional,
    forces: FieldAvailability::Optional,
    time: FieldAvailability::Required,
    step: FieldAvailability::Required,
    properties: FieldAvailability::Optional,
};

const TRR_MAGIC: i32 = 1993;
const TRR_VERSION: &[u8] = b"GMX_trn_file";
/// Stable frame-property key used to preserve the GROMACS TRR lambda scalar.
pub const TRR_LAMBDA_PROPERTY: &str = "gromacs.trr.lambda";

/// Native scalar width used by one TRR writer or frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum TrrScalarPrecision {
    /// IEEE 754 single precision (GROMACS mixed-precision builds).
    Float32,
    /// IEEE 754 double precision.
    Float64,
}

impl TrrScalarPrecision {
    const fn bytes(self) -> usize {
        match self {
            Self::Float32 => 4,
            Self::Float64 => 8,
        }
    }

    /// The precision's name in messages.
    const fn name(self) -> &'static str {
        match self {
            Self::Float32 => "single precision",
            Self::Float64 => "double precision",
        }
    }
}

/// Explicit policy for the TRR lambda scalar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum TrrLambdaPolicy {
    /// Preserve lambda in the `gromacs.trr.lambda` frame property.
    FrameProperty,
    /// Require lambda to be exactly zero and do not publish a property.
    RequireZero,
}

/// TRR reader policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrrReadOptions {
    lambda_policy: TrrLambdaPolicy,
}

impl Default for TrrReadOptions {
    fn default() -> Self {
        Self {
            lambda_policy: TrrLambdaPolicy::FrameProperty,
        }
    }
}

impl TrrReadOptions {
    /// Sets [`Self::lambda_policy`].
    #[must_use]
    pub const fn with_lambda_policy(mut self, lambda_policy: TrrLambdaPolicy) -> Self {
        self.lambda_policy = lambda_policy;
        self
    }

    /// How lambda is published; [`TrrLambdaPolicy::FrameProperty`] by default.
    pub const fn lambda_policy(&self) -> TrrLambdaPolicy {
        self.lambda_policy
    }
}

/// TRR writer policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrrWriteOptions {
    precision: TrrScalarPrecision,
    lambda_policy: TrrLambdaPolicy,
    geometry_metadata: UnrepresentablePolicy,
    properties: UnrepresentablePolicy,
    frame_timing: FrameTimingPolicy,
    frame_interval: Quantity<f64>,
}

impl Default for TrrWriteOptions {
    fn default() -> Self {
        Self {
            precision: TrrScalarPrecision::Float32,
            lambda_policy: TrrLambdaPolicy::FrameProperty,
            geometry_metadata: UnrepresentablePolicy::Fail,
            properties: UnrepresentablePolicy::Fail,
            frame_timing: FrameTimingPolicy::Require,
            frame_interval: Quantity::new(1.0, PICOSECOND),
        }
    }
}

impl TrrWriteOptions {
    /// Selects the precision for every scalar, including periodic-cell vectors.
    /// Cells that become degenerate at that precision are rejected before any
    /// frame bytes are written; use Float64 when Float32 cannot preserve a cell.
    #[must_use]
    pub const fn with_precision(mut self, precision: TrrScalarPrecision) -> Self {
        self.precision = precision;
        self
    }

    /// Sets [`Self::lambda_policy`].
    #[must_use]
    pub const fn with_lambda_policy(mut self, lambda_policy: TrrLambdaPolicy) -> Self {
        self.lambda_policy = lambda_policy;
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
        self.frame_interval = frame_interval(TrajectoryFormat::Trr, interval)?;
        Ok(self)
    }

    /// The scalar width written; [`TrrScalarPrecision::Float32`] by default.
    pub const fn precision(&self) -> TrrScalarPrecision {
        self.precision
    }

    /// Where lambda comes from; [`TrrLambdaPolicy::FrameProperty`] by default,
    /// which requires every frame to carry exactly the
    /// [`TRR_LAMBDA_PROPERTY`] frame property.
    pub const fn lambda_policy(&self) -> TrrLambdaPolicy {
        self.lambda_policy
    }

    /// A frame's occupancies and B-factors, which TRR cannot carry; by
    /// default a frame with either fails to write.
    pub const fn geometry_metadata(&self) -> UnrepresentablePolicy {
        self.geometry_metadata
    }

    /// A frame's atom and bond properties, and its own properties other than
    /// the [`TRR_LAMBDA_PROPERTY`] the [`TrrLambdaPolicy::FrameProperty`]
    /// lambda policy writes, which TRR cannot carry; by default a frame with
    /// any of them fails to write. The lambda policy's own requirements are
    /// never omitted.
    pub const fn properties(&self) -> UnrepresentablePolicy {
        self.properties
    }

    /// The time and step of a frame without them, which TRR records for
    /// every frame; [`FrameTimingPolicy::Require`] by default.
    pub const fn frame_timing(&self) -> FrameTimingPolicy {
        self.frame_timing
    }

    /// The time between frames [`FrameTimingPolicy::Number`] gives; 1 ps by
    /// default.
    pub const fn frame_interval(&self) -> Quantity<f64> {
        self.frame_interval
    }

    /// The policies of the shared check; the lambda policy decides about
    /// frame properties.
    const fn state_policies(&self) -> StatePolicies {
        StatePolicies {
            carries: FIELDS.carries(),
            geometry_metadata: self.geometry_metadata,
            properties: self.properties,
            cell: UnrepresentablePolicy::Fail,
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct TrrFrameHeader {
    precision: TrrScalarPrecision,
    atom_count: usize,
    box_size: usize,
    x_size: usize,
    v_size: usize,
    f_size: usize,
    step: u64,
    time: f64,
    lambda: f64,
    header_bytes: u64,
}

impl TrrFrameHeader {
    pub(crate) const fn precision(&self) -> TrrScalarPrecision {
        self.precision
    }

    pub(crate) const fn atom_count(&self) -> usize {
        self.atom_count
    }

    pub(crate) const fn has_cell(&self) -> bool {
        self.box_size != 0
    }

    pub(crate) const fn has_velocities(&self) -> bool {
        self.v_size != 0
    }

    pub(crate) const fn has_forces(&self) -> bool {
        self.f_size != 0
    }

    /// The time, in picoseconds, and step the frame records.
    const fn stamp(&self) -> FrameStamp {
        FrameStamp {
            time: Some(self.time),
            step: Some(self.step),
        }
    }
}

fn trr_lambda_key() -> PropertyKey {
    PropertyKey::new(TRR_LAMBDA_PROPERTY).expect("TRR lambda property key is valid")
}

fn writer_unit(
    source_label: &str,
    field: &str,
    error: impl std::error::Error + Send + Sync + 'static,
) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::InconsistentMetadata,
        TrajectoryIoOperation::WriteFrame,
        Some(FORMAT),
        source_label,
        format!("The {field} unit is incompatible."),
    )
    .with_source(error)
}
