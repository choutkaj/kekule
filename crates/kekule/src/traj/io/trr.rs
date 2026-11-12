//! Pure-Rust GROMACS TRR/XDR trajectory I/O.
//!
//! [`TrrReader`] and [`IndexedTrrReader`] read big-endian XDR frames with f32
//! or f64 position, box, velocity and force blocks; [`TrrWriter`] writes them
//! at one explicit [`TrrScalarPrecision`]. Every frame carries a step and a
//! time; GROMACS nanometre, picosecond and kJ/(mol nm) units are converted at
//! the codec boundary. The lambda scalar follows a [`TrrLambdaPolicy`]. A
//! frame without a time or step can be numbered under
//! [`TrrWriteOptions::frame_timing`]; the properties TRR cannot carry are
//! rejected unless [`TrrWriteOptions::frame_properties`] omits them, as
//! [`TrrWriter::report`] then reports.
//!
//! # Examples
//!
//! ```
//! use std::io::Cursor;
//! use std::sync::Arc;
//!
//! use kekule::geometry::{Point3, Vector3};
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
//! frame.sample_mut().set_positions(Quantity::new(points, NANOMETER))?;
//! let velocity = Vector3::new(1.0, 0.0, 0.0);
//! frame.set_velocities(Some(Quantity::new([velocity], NANOMETER_PER_PICOSECOND)))?;
//! frame.set_step(Some(500));
//! frame.set_time(Some(Quantity::new(1.0, PICOSECOND)))?;
//!
//! // Frames without a lambda property need the RequireZero policy.
//! let options = TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero);
//! let mut writer = TrrWriter::new(Vec::new(), Arc::clone(&topology), options, "memory")?;
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

use crate::traj::{FrameView, TrajectoryField};
use std::fmt;
use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::Arc;

use crate::diagnostics::{LimitExceeded, Resource, UnrepresentablePolicy, WriteReport};
use crate::geometry::{PeriodicCell, Point3, Vector3};
use crate::properties::{OwnerProperties, PropertyKey, PropertyValue, StructureProperties};
use crate::structure::DecodedSample;
use crate::topology::Topology;
use crate::traj::{
    FrameBuffer, FrameData, SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoError,
    TrajectoryIoErrorKind, TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use crate::units::{
    Quantity, CANONICAL_LENGTH_UNIT, DIMENSIONLESS, KILOJOULE_PER_MOLE_PER_NANOMETER, NANOMETER,
    NANOMETER_PER_PICOSECOND, PICOSECOND,
};

use super::frame_state::{record, unrepresentable_frame_state};
use super::limits::{exceeded_bytes, filled_bytes};
use super::{
    codec_context, frame_offset_context, index_mismatch, io_context, probe_seekable_eof,
    projected_index_limit, require_no_partial_frame, require_nonempty_writer, validation_cause,
    FrameOffsets, FrameStamp, FrameTimingPolicy, IndexBuildControl, TrajectoryIoLimits,
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
    limits: TrajectoryIoLimits,
    source_label: String,
}

impl Default for TrrReadOptions {
    fn default() -> Self {
        Self {
            limits: TrajectoryIoLimits::default(),
            source_label: "trr stream".into(),
            lambda_policy: TrrLambdaPolicy::FrameProperty,
        }
    }
}

impl TrrReadOptions {
    /// Sets limits applied before allocation, scanning, and seeking.
    /// `max_scratch_bytes` bounds the three dense decode arrays, one vector
    /// per atom each, and the reusable raw record buffer, in both sequential
    /// and indexed modes. The frame index has its own `max_index_bytes` limit.
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

    /// Sets [`Self::lambda_policy`].
    #[must_use]
    pub const fn with_lambda_policy(mut self, lambda_policy: TrrLambdaPolicy) -> Self {
        self.lambda_policy = lambda_policy;
        self
    }

    /// Limits applied before allocation, scanning, and seeking.
    pub const fn limits(&self) -> &TrajectoryIoLimits {
        &self.limits
    }

    /// The source label used in decoding errors.
    pub fn source_label(&self) -> &str {
        &self.source_label
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
    frame_properties: UnrepresentablePolicy,
    frame_timing: FrameTimingPolicy,
}

impl Default for TrrWriteOptions {
    fn default() -> Self {
        Self {
            precision: TrrScalarPrecision::Float32,
            lambda_policy: TrrLambdaPolicy::FrameProperty,
            frame_properties: UnrepresentablePolicy::Fail,
            frame_timing: FrameTimingPolicy::Require,
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

    /// Sets [`Self::frame_properties`].
    #[must_use]
    pub const fn with_frame_properties(mut self, frame_properties: UnrepresentablePolicy) -> Self {
        self.frame_properties = frame_properties;
        self
    }

    /// Sets [`Self::frame_timing`]; [`TrrWriter::new`] checks its interval.
    #[must_use]
    pub const fn with_frame_timing(mut self, frame_timing: FrameTimingPolicy) -> Self {
        self.frame_timing = frame_timing;
        self
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

    /// A trajectory frame's atom and bond properties, and its frame
    /// properties other than the [`TRR_LAMBDA_PROPERTY`] the
    /// [`TrrLambdaPolicy::FrameProperty`] lambda policy writes, which TRR
    /// cannot carry: by default a frame with any of them fails to write;
    /// [`UnrepresentablePolicy::Omit`] writes the frame without them and
    /// reports a
    /// [`WriteOmission::StructureState`](crate::diagnostics::WriteOmission::StructureState)
    /// for each kind left out (see [`TrrWriter::report`]). The lambda
    /// policy's own requirements are never omitted.
    pub const fn frame_properties(&self) -> UnrepresentablePolicy {
        self.frame_properties
    }

    /// The time and step of a frame without them, which TRR records for
    /// every frame; [`FrameTimingPolicy::Require`], under which such a frame
    /// fails to write, by default.
    pub const fn frame_timing(&self) -> FrameTimingPolicy {
        self.frame_timing
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

/// Sequential TRR reader retaining one seekable stream and reusable scratch.
pub struct TrrReader<R> {
    reader: R,
    topology: Arc<Topology>,
    options: TrrReadOptions,
    stream_start: u64,
    current_header_offset: u64,
    first_header: TrrFrameHeader,
    pending_header: Option<TrrFrameHeader>,
    positions: Vec<Point3>,
    velocities: Vec<Vector3>,
    forces: Vec<Vector3>,
    raw: Vec<u8>,
    properties: OwnerProperties,
    frame_cursor: usize,
    precision_mixed: bool,
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
struct TrrDecodedFrame {
    header: TrrFrameHeader,
    cell: Option<PeriodicCell>,
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
        validate_atom_count(atom_count, limits, source_label)?;
        let stream_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Trr),
                source_label,
                error,
            )
        })?;
        let first_header = read_header(&mut reader, limits, source_label, Some(atom_count), false)?
            .ok_or_else(|| header_error(source_label, "The stream is empty."))?;
        validate_lambda(first_header.lambda, options.lambda_policy, source_label, 0)?;
        let mut positions = Vec::new();
        let mut velocities = Vec::new();
        let mut forces = Vec::new();
        let vector_bytes = atom_count
            .checked_mul(std::mem::size_of::<Vector3>())
            .ok_or_else(|| {
                resource_error(
                    source_label,
                    None,
                    "The vector buffers exceed this platform's address range.",
                )
            })?;
        let scratch_bytes = vector_bytes.checked_mul(3).ok_or_else(|| {
            resource_error(
                source_label,
                None,
                "The vector buffers exceed this platform's address range.",
            )
        })?;
        if scratch_bytes > limits.max_scratch_bytes() {
            return Err(resource_error(
                source_label,
                None,
                format!(
                    "The vector buffers need {scratch_bytes} bytes, more than the buffer limit of {} bytes.",
                    limits.max_scratch_bytes()
                ),
            ));
        }
        positions.try_reserve_exact(atom_count).map_err(|_| {
            resource_error(source_label, None, "Cannot allocate the position buffer.")
        })?;
        velocities.try_reserve_exact(atom_count).map_err(|_| {
            resource_error(source_label, None, "Cannot allocate the velocity buffer.")
        })?;
        forces
            .try_reserve_exact(atom_count)
            .map_err(|_| resource_error(source_label, None, "Cannot allocate the force buffer."))?;
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

    fn next_header(&mut self) -> Result<Option<TrrFrameHeader>, TrajectoryIoError> {
        if let Some(header) = self.pending_header.take() {
            return Ok(Some(header));
        }
        self.current_header_offset = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Trr),
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

    fn parse_next(&mut self) -> Result<Option<TrrDecodedFrame>, TrajectoryIoError> {
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
            return Err(
                resource_error(&self.options.source_label, Some(self.frame_cursor), "").exceeding(
                    "The trajectory",
                    LimitExceeded::first_beyond(Resource::Frames, self.options.limits.max_frames()),
                ),
            );
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
                resource_error(
                    &self.options.source_label,
                    Some(self.frame_cursor),
                    "The frame's size exceeds this platform's address range.",
                )
            })?;
        let frame_bytes = header
            .header_bytes
            .checked_add(payload_bytes as u64)
            .ok_or_else(|| {
                resource_error(
                    &self.options.source_label,
                    Some(self.frame_cursor),
                    "The frame's size exceeds the 64-bit range.",
                )
            })?;
        if frame_bytes > self.options.limits.max_frame_bytes() {
            return Err(
                resource_error(&self.options.source_label, Some(self.frame_cursor), "").exceeding(
                    "The frame",
                    exceeded_bytes(self.options.limits.max_frame_bytes(), frame_bytes),
                ),
            );
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
            return Err(resource_error(
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
                        &self.options.source_label,
                        self.frame_cursor,
                        "The lambda value is not a valid frame property.",
                    )
                })?;
        }
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            resource_error(
                &self.options.source_label,
                Some(self.frame_cursor),
                "The frame number exceeds this platform's address range.",
            )
        })?;
        Ok(Some(TrrDecodedFrame { header, cell }))
    }

    fn publish(
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
    fn index_frames(
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

    fn position(&mut self, operation: TrajectoryIoOperation) -> Result<u64, TrajectoryIoError> {
        self.reader.stream_position().map_err(|error| {
            io_context(
                operation,
                Some(TrajectoryFormat::Trr),
                &self.options.source_label,
                error,
            )
        })
    }

    fn seek_to(
        &mut self,
        offset: SeekFrom,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.reader.seek(offset).map_err(|error| {
            io_context(
                operation,
                Some(TrajectoryFormat::Trr),
                &self.options.source_label,
                error,
            )
        })
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
                    Some(TrajectoryFormat::Trr),
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

/// Fully verified indexed TRR reader.
///
/// Sequential and random reads share one set of decode scratch arrays. Random
/// reads restore stream position and header state before publishing a frame.
pub struct IndexedTrrReader<R> {
    inner: TrrReader<R>,
    offsets: Arc<FrameOffsets>,
    /// Whether the indexed frames mix scalar widths.
    precision_mixed: bool,
    /// Bytes after the last indexed frame when the index last grew.
    pending_bytes: u64,
}

impl<R> fmt::Debug for IndexedTrrReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedTrrReader")
            .field("reader", &self.inner)
            .field("frame_count", &self.offsets.len())
            .finish_non_exhaustive()
    }
}

impl<R: Read + Seek> IndexedTrrReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order
    /// and indexes it for random access.
    ///
    /// The stream's position at construction is frame 0. Every frame from there to
    /// the end is decoded and verified before the reader is returned, and both
    /// sequential and indexed reads then address frames from that start. Topology,
    /// limits, and diagnostic labels follow [`TrrReader::new`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`TrrReader::new`], the decoding error of the
    /// first invalid frame, [`TrajectoryIoErrorKind::LimitExceeded`] when the
    /// index outgrows the limits, and [`TrajectoryIoErrorKind::Io`] when
    /// reading or seeking fails.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: TrrReadOptions,
    ) -> Result<Self, TrajectoryIoError> {
        Self::build(
            reader,
            topology.into(),
            options,
            false,
            &IndexBuildControl::default(),
        )
    }

    /// Indexes the stream like [`Self::new`], reporting to `control`; a
    /// `growing` stream may end in a partial frame (see [`Self::extend`]).
    pub(crate) fn build(
        reader: R,
        topology: Arc<Topology>,
        options: TrrReadOptions,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        TrrReader::new(reader, topology, options)?.index_frames(growing, control)
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; random reads verify each frame against them, including
    /// that its scalar width is the first frame's unless `precision_mixed`.
    pub(crate) fn from_offsets(
        reader: R,
        topology: Arc<Topology>,
        options: TrrReadOptions,
        offsets: Arc<FrameOffsets>,
        precision_mixed: bool,
    ) -> Result<Self, TrajectoryIoError> {
        let mut inner = TrrReader::new(reader, topology, options)?;
        offsets.verify_start(
            inner.stream_start,
            TrajectoryFormat::Trr,
            &inner.options.source_label,
        )?;
        inner.precision_mixed = precision_mixed;
        Ok(Self {
            inner,
            offsets,
            precision_mixed,
            pending_bytes: 0,
        })
    }

    pub(crate) fn offsets(&self) -> &FrameOffsets {
        &self.offsets
    }

    pub(crate) fn shared_offsets(&self) -> &Arc<FrameOffsets> {
        &self.offsets
    }

    /// Lets go of the index until [`Self::share_index`], so that a reader
    /// extending it does not copy it.
    pub(crate) fn release_index(&mut self) {
        self.offsets = Arc::new(FrameOffsets::empty(self.offsets.end()));
    }

    /// Reads by the index of `other`, a reader of the same file, from now
    /// on.
    pub(crate) fn share_index(&mut self, other: &Self) {
        self.offsets = Arc::clone(&other.offsets);
        self.precision_mixed = other.precision_mixed;
    }

    /// Another reader of the same frames over `stream`, which must read the
    /// same file from its start; it shares the index and reads the first
    /// frame's header alone.
    pub(crate) fn fork(&self, stream: R) -> Result<Self, TrajectoryIoError> {
        Self::from_offsets(
            stream,
            Arc::clone(&self.inner.topology),
            self.inner.options.clone(),
            Arc::clone(&self.offsets),
            self.precision_mixed,
        )
    }

    /// The stream every read goes through.
    pub(crate) fn stream(&self) -> &R {
        &self.inner.reader
    }

    /// The stream every read goes through.
    pub(crate) fn stream_mut(&mut self) -> &mut R {
        &mut self.inner.reader
    }

    /// Reads frame `index` as [`SeekableTrajectoryReader::read_frame`]
    /// does and hands over its positions, velocities, and forces, each
    /// written once as it was decoded. The read allocates them on the
    /// calling thread and the reader keeps no buffer for its next read,
    /// which may run on another.
    pub(crate) fn read_decoded(
        &mut self,
        index: usize,
    ) -> Result<DecodedSample, TrajectoryIoError> {
        let read = self.read_indexed(index);
        // Taken whether the read succeeds or fails, and whether the frame
        // holds velocities and forces or not, so that the reader keeps no
        // array for a read on another thread.
        let positions = std::mem::take(&mut self.inner.positions);
        let velocities = std::mem::take(&mut self.inner.velocities);
        let forces = std::mem::take(&mut self.inner.forces);
        let decoded = read?;
        Ok(DecodedSample {
            positions: Quantity::new(positions, NANOMETER),
            cell: decoded.cell,
            velocities: decoded
                .header
                .has_velocities()
                .then(|| Quantity::new(velocities, NANOMETER_PER_PICOSECOND)),
            forces: decoded
                .header
                .has_forces()
                .then(|| Quantity::new(forces, KILOJOULE_PER_MOLE_PER_NANOMETER)),
            owner: Some(self.inner.properties.clone()),
        })
    }

    /// Bytes after the last indexed frame when the index last grew.
    pub(crate) const fn pending_bytes(&self) -> u64 {
        self.pending_bytes
    }

    /// Drops the frames indexed after `mark`, with the mixed-precision flag
    /// and pending bytes the index had then; see [`FrameOffsets::mark`].
    pub(crate) fn restore_index(
        &mut self,
        mark: (usize, u64),
        precision_mixed: bool,
        pending_bytes: u64,
    ) {
        Arc::make_mut(&mut self.offsets).restore(mark);
        self.precision_mixed = precision_mixed;
        self.pending_bytes = pending_bytes;
    }

    /// Verifies and indexes every complete frame after the last indexed one,
    /// up to the end of the stream, and returns the bytes left after them.
    ///
    /// A `growing` stream may end in a frame still being written: the scan
    /// stops before a frame the stream ends inside, which is not yet an
    /// error. Otherwise such a frame is an error, as is any frame that is
    /// invalid with all its bytes present. Failure leaves the index as it
    /// was; the sequential cursor is kept either way.
    pub(crate) fn extend(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        let mark = self.offsets.mark();
        let saved_offset = self.inner.position(TrajectoryIoOperation::Index)?;
        let saved_cursor = self.inner.frame_cursor;
        let saved_pending = self.inner.pending_header.clone();
        let saved_header_offset = self.inner.current_header_offset;
        let saved_precision_mixed = self.inner.precision_mixed;
        let end = self.offsets.end();
        // A freshly opened reader has already read frame 0's header there.
        if saved_pending.is_none() || saved_header_offset != end {
            self.inner.pending_header = None;
            self.inner
                .seek_to(SeekFrom::Start(end), TrajectoryIoOperation::Index)?;
        }
        self.inner.current_header_offset = end;
        self.inner.frame_cursor = self.offsets.len();
        self.inner.precision_mixed = self.precision_mixed;
        let scanned = self.scan(growing, control);
        let scanned_mixed = self.inner.precision_mixed;
        let restored = self
            .inner
            .seek_to(SeekFrom::Start(saved_offset), TrajectoryIoOperation::Index);
        self.inner.frame_cursor = saved_cursor;
        self.inner.pending_header = saved_pending;
        self.inner.current_header_offset = saved_header_offset;
        self.inner.precision_mixed = saved_precision_mixed;
        match scanned.and_then(|pending| restored.map(|_| pending)) {
            Ok(pending) => {
                self.precision_mixed = scanned_mixed;
                self.pending_bytes = pending;
                Ok(pending)
            }
            Err(error) => {
                Arc::make_mut(&mut self.offsets).restore(mark);
                Err(error)
            }
        }
    }

    /// Indexes frames from the current position; see [`Self::extend`].
    fn scan(
        &mut self,
        growing: bool,
        control: &IndexBuildControl<'_>,
    ) -> Result<u64, TrajectoryIoError> {
        let label = self.inner.options.source_label.clone();
        control.checkpoint(
            self.offsets.len(),
            self.offsets.end(),
            TrajectoryFormat::Trr,
            &label,
        )?;
        loop {
            if let Some(limit) =
                projected_index_limit(self.offsets.len(), &self.inner.options.limits)
            {
                let clean_eof = self.inner.pending_header.is_none()
                    && probe_seekable_eof(
                        &mut self.inner.reader,
                        TrajectoryIoOperation::Index,
                        TrajectoryFormat::Trr,
                        &label,
                    )?;
                if clean_eof {
                    return Ok(0);
                }
                return Err(resource_error(&label, Some(self.offsets.len()), "")
                    .exceeding("The frame index", limit));
            }
            let offset = if self.inner.pending_header.is_some() {
                self.inner.current_header_offset
            } else {
                self.inner.position(TrajectoryIoOperation::Index)?
            };
            let decoded = match self.inner.parse_next() {
                Ok(Some(decoded)) => decoded,
                Ok(None) => return Ok(0),
                Err(error) if growing && error.kind() == TrajectoryIoErrorKind::TruncatedRecord => {
                    let file_end = self
                        .inner
                        .seek_to(SeekFrom::End(0), TrajectoryIoOperation::Index)?;
                    return Ok(file_end.saturating_sub(offset));
                }
                Err(error) => {
                    return Err(frame_offset_context(error, self.inner.frame_cursor, offset))
                }
            };
            let end = self.inner.position(TrajectoryIoOperation::Index)?;
            Arc::make_mut(&mut self.offsets).push(
                offset,
                end,
                decoded.header.stamp(),
                &self.inner.options.limits,
                TrajectoryFormat::Trr,
                &label,
            )?;
            control.checkpoint(self.offsets.len(), end, TrajectoryFormat::Trr, &label)?;
        }
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.inner.topology()
    }

    pub(crate) fn first_header(&self) -> &TrrFrameHeader {
        self.inner.first_header()
    }

    pub(crate) const fn precision_mixed(&self) -> bool {
        self.precision_mixed
    }
}

impl<R: Read + Seek> TrajectoryReader for IndexedTrrReader<R> {
    fn topology(&self) -> &Topology {
        self.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.inner.shared_topology()
    }

    /// Reads the indexed frames in order and ends after the last of them,
    /// even when the stream holds more.
    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        if self.topology().shares_layout(destination.topology())
            && self.inner.frame_cursor >= self.offsets.len()
        {
            return Ok(false);
        }
        self.inner.read_next(destination)
    }
}

impl<R: Read + Seek> SeekableTrajectoryReader for IndexedTrrReader<R> {
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
        let decoded = self.read_indexed(index)?;
        self.inner.publish(
            &self.inner.positions,
            &self.inner.velocities,
            &self.inner.forces,
            &self.inner.properties,
            &decoded,
            destination,
        )
    }
}

impl<R: Read + Seek> IndexedTrrReader<R> {
    /// Decodes and verifies frame `index` into the reader's buffers, then
    /// puts the stream back where it was.
    fn read_indexed(&mut self, index: usize) -> Result<TrrDecodedFrame, TrajectoryIoError> {
        let (offset, _) = self
            .offsets
            .frame(index)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        let saved_offset = self.inner.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Trr),
                &self.inner.options.source_label,
                error,
            )
        })?;
        let saved_cursor = self.inner.frame_cursor;
        let saved_pending = self.inner.pending_header.clone();
        let saved_header_offset = self.inner.current_header_offset;
        let saved_precision_mixed = self.inner.precision_mixed;
        self.inner
            .reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Trr),
                    &self.inner.options.source_label,
                    error,
                )
            })?;
        self.inner.pending_header = None;
        self.inner.current_header_offset = offset;
        self.inner.frame_cursor = index;
        let result = self
            .inner
            .parse_next()
            .map_err(|error| frame_offset_context(error, index, offset))
            .and_then(|decoded| decoded.ok_or(TrajectoryIoError::frame_index_out_of_range(index)))
            .and_then(|decoded| {
                self.offsets.verify_stream_frame_end(
                    index,
                    &mut self.inner.reader,
                    TrajectoryFormat::Trr,
                    &self.inner.options.source_label,
                )?;
                self.offsets.verify_stamp(
                    index,
                    decoded.header.stamp(),
                    TrajectoryIoOperation::ReadFrame,
                    TrajectoryFormat::Trr,
                    &self.inner.options.source_label,
                )?;
                if !self.precision_mixed
                    && decoded.header.precision != self.inner.first_header.precision
                {
                    return Err(index_mismatch(
                        TrajectoryIoOperation::ReadFrame,
                        TrajectoryFormat::Trr,
                        &self.inner.options.source_label,
                        "The frame's scalar width differs from the first frame's.",
                    )
                    .with_frame(index));
                }
                Ok(decoded)
            });
        let restore = self
            .inner
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Trr),
                    &self.inner.options.source_label,
                    error,
                )
            });
        self.inner.frame_cursor = saved_cursor;
        self.inner.pending_header = saved_pending;
        self.inner.current_header_offset = saved_header_offset;
        self.inner.precision_mixed = saved_precision_mixed;
        let decoded = result?;
        restore?;
        Ok(decoded)
    }
}

/// Pure-Rust TRR writer.
pub struct TrrWriter<W> {
    writer: W,
    topology: Arc<Topology>,
    options: TrrWriteOptions,
    source_label: String,
    frame_count: usize,
    /// Set while a frame's bytes are being written and left set if writing
    /// fails, so a stream that may end in a partial frame is never extended
    /// or reported as finished.
    partial: bool,
    raw: Vec<u8>,
    report: WriteReport,
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
    /// no atoms or more than `i32::MAX / 3`, and
    /// [`TrajectoryIoErrorKind::InconsistentMetadata`] for a
    /// [`FrameTimingPolicy::Number`] interval that is not a finite,
    /// nonnegative time.
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
        options
            .frame_timing
            .validate(TrajectoryFormat::Trr, &source_label)?;
        Ok(Self {
            writer,
            topology,
            options,
            source_label,
            frame_count: 0,
            partial: false,
            raw: Vec::new(),
            report: WriteReport::default(),
        })
    }

    /// The frame state the frames written so far leave out under
    /// [`TrrWriteOptions::frame_properties`], counted over every frame
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
        let omitted = unrepresentable_frame_state(
            [
                (
                    properties.atoms().has_data(),
                    TrajectoryField::AtomProperties,
                    self.options.frame_properties,
                ),
                (
                    properties.bonds().has_data(),
                    TrajectoryField::BondProperties,
                    self.options.frame_properties,
                ),
                (
                    self.options.frame_properties == UnrepresentablePolicy::Omit
                        && has_unwritten_owner_properties(
                            properties.owner(),
                            self.options.lambda_policy,
                        ),
                    TrajectoryField::Properties,
                    self.options.frame_properties,
                ),
            ],
            |field| writer_field_error(&self.source_label, self.frame_count, field),
        )?;
        if self.frame_count == usize::MAX {
            return Err(writer_limit(
                &self.source_label,
                "The frame count exceeds this platform's address range.",
            ));
        }
        let (step, time) =
            self.options
                .frame_timing
                .timing(self.frame_count, frame.step(), frame.time());
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
            self.options.frame_properties,
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

fn read_header<R: Read>(
    reader: &mut R,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    expected_atoms: Option<usize>,
    clean_eof: bool,
) -> Result<Option<TrrFrameHeader>, TrajectoryIoError> {
    let mut magic = [0_u8; 4];
    match reader.read(&mut magic[..1]) {
        Ok(0) if clean_eof => return Ok(None),
        Ok(0) => return Err(header_error(source_label, "The magic number is missing.")),
        Ok(_) => {}
        Err(error) => {
            return Err(io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Trr),
                source_label,
                error,
            ))
        }
    }
    read_exact(
        reader,
        &mut magic[1..],
        source_label,
        None,
        "The magic number is incomplete.",
    )?;
    if i32::from_be_bytes(magic) != TRR_MAGIC {
        return Err(header_error(source_label, "The magic number is not 1993."));
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
            Some(TrajectoryFormat::Trr),
            source_label,
            "The version length is not that of GMX_trn_file.",
        ));
    }
    let (version_matches, version_bytes) = read_xdr_string(reader, limits, source_label)?;
    if !version_matches {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "The version string is not GMX_trn_file.",
        ));
    }
    let mut sizes = [0_usize; 10];
    for size in &mut sizes {
        let value = read_i32(reader, source_label, None, "A block size is truncated.")?;
        *size = usize::try_from(value)
            .map_err(|_| header_error(source_label, "A block size is negative."))?;
        if *size as u64 > limits.max_record_bytes() {
            return Err(resource_error(source_label, None, "").exceeding(
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
    .map_err(|_| header_error(source_label, "The atom count is negative."))?;
    validate_atom_count(atom_count, limits, source_label)?;
    if let Some(expected) = expected_atoms {
        if atom_count != expected {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Trr),
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
            Some(TrajectoryFormat::Trr),
            source_label,
            "The frame has input record, energy, virial, pressure, topology or symbol blocks. Kekule cannot keep them.",
        ));
    }
    if x_size == 0 {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "The frame has no positions.",
        ));
    }
    let precision = infer_precision(atom_count, box_size, x_size, v_size, f_size, source_label)?;
    let expected_vectors = atom_count
        .checked_mul(3)
        .and_then(|count| count.checked_mul(precision.bytes()))
        .ok_or_else(|| {
            resource_error(
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
            Some(TrajectoryFormat::Trr),
            source_label,
            "The block sizes do not match the atom count and the scalar precision.",
        ));
    }
    let step_raw = read_i32(reader, source_label, None, "The step is truncated.")?;
    let step = u64::try_from(step_raw).map_err(|_| {
        codec_context(
            TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
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
            Some(TrajectoryFormat::Trr),
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
            resource_error(
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

fn infer_precision(
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
                source_label,
                "The frame has no block that shows its precision.",
            )
        })?;
    let width = size
        .checked_div(scalars)
        .filter(|_| size % scalars == 0)
        .ok_or_else(|| header_error(source_label, "The block sizes give no whole scalar width."))?;
    match width {
        4 => Ok(TrrScalarPrecision::Float32),
        8 => Ok(TrrScalarPrecision::Float64),
        _ => Err(codec_context(
            TrajectoryIoErrorKind::InvalidPrecision,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            format!("The scalar width is {width} bytes, not 4 or 8."),
        )),
    }
}

fn read_xdr_string<R: Read>(
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
    .map_err(|_| header_error(source_label, "The version string's length is negative."))?;
    if length > 256 || length > limits.max_text_line_bytes() || length > limits.max_scratch_bytes()
    {
        return Err(resource_error(
            source_label,
            None,
            format!("The version string has {length} bytes, more than 256 or the line limit."),
        ));
    }
    let padded = length
        .checked_add(3)
        .map(|value| value & !3)
        .ok_or_else(|| {
            resource_error(
                source_label,
                None,
                "The version string's padding exceeds this platform's address range.",
            )
        })?;
    let mut bytes = [0_u8; 256];
    read_exact(
        reader,
        &mut bytes[..padded],
        source_label,
        None,
        "The version string is truncated.",
    )?;
    if bytes[length..].iter().any(|byte| *byte != 0) {
        return Err(header_error(
            source_label,
            "The version string's padding is not zero.",
        ));
    }
    Ok((&bytes[..length] == TRR_VERSION, 4 + padded as u64))
}

fn read_raw<R: Read>(
    reader: &mut R,
    raw: &mut Vec<u8>,
    length: usize,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    frame: usize,
    field: &str,
) -> Result<(), TrajectoryIoError> {
    if length > limits.max_scratch_bytes() || length as u64 > limits.max_record_bytes() {
        return Err(resource_error(
            source_label,
            Some(frame),
            format!(
                "The {field} block has {length} bytes, more than the record limit of {} bytes or the buffer limit of {} bytes.",
                limits.max_record_bytes(),
                limits.max_scratch_bytes()
            ),
        ));
    }
    if raw.capacity() < length {
        raw.try_reserve_exact(length.saturating_sub(raw.len()))
            .map_err(|_| {
                resource_error(
                    source_label,
                    Some(frame),
                    format!("Cannot allocate the {field} buffer."),
                )
            })?;
    }
    raw.resize(length, 0);
    read_exact(
        reader,
        raw,
        source_label,
        Some(frame),
        &format!("truncated {field}"),
    )
}

fn decode_triplets<T>(
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
fn decode_triplets_into<T>(
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
                source_label,
                frame,
                format!("A {field} value is not finite."),
            ));
        }
        destination.push(construct(x, y, z));
    }
    Ok(())
}

fn decode_cell(
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
        frame_error(source_label, frame, "The periodic cell is invalid.").with_source(error)
    })
}

fn scalar_at(bytes: &[u8], index: usize, precision: TrrScalarPrecision) -> f64 {
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

fn encode_points_to_raw(
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

fn encode_vectors_to_raw(
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

fn encode_triplets_to_raw(
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

fn push_scalar(
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

fn validate_representable(
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

fn trr_lambda_key() -> PropertyKey {
    PropertyKey::new(TRR_LAMBDA_PROPERTY).expect("TRR lambda property key is valid")
}

/// Whether `owner` holds frame properties a frame written under `policy`
/// leaves out: any under [`TrrLambdaPolicy::RequireZero`], and any but the
/// lambda under [`TrrLambdaPolicy::FrameProperty`].
fn has_unwritten_owner_properties(owner: &OwnerProperties, policy: TrrLambdaPolicy) -> bool {
    match policy {
        TrrLambdaPolicy::FrameProperty => {
            owner.len() > usize::from(owner.get(&trr_lambda_key()).is_some())
        }
        TrrLambdaPolicy::RequireZero => !owner.is_empty(),
    }
}

/// The lambda of a frame with `properties`. Under
/// [`UnrepresentablePolicy::Omit`] for `frame_properties`, properties
/// [`has_unwritten_owner_properties`] names are ignored, and atom and bond
/// properties are left to the caller.
fn writer_lambda(
    properties: &StructureProperties,
    policy: TrrLambdaPolicy,
    frame_properties: UnrepresentablePolicy,
    source_label: &str,
    frame: usize,
) -> Result<f64, TrajectoryIoError> {
    let omit = frame_properties == UnrepresentablePolicy::Omit;
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
                    source_label,
                    frame,
                    TrajectoryField::Properties,
                ));
            }
            Ok(0.0)
        }
    }
}

fn validate_lambda(
    lambda: f64,
    policy: TrrLambdaPolicy,
    source_label: &str,
    frame: usize,
) -> Result<(), TrajectoryIoError> {
    if policy == TrrLambdaPolicy::RequireZero && lambda != 0.0 {
        return Err(TrajectoryIoError::codec(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::ReadFrame,
            Some(TrajectoryFormat::Trr),
        )
        .with_source_label(source_label)
        .with_frame(frame)
        .with_detail("The lambda value is not zero."));
    }
    Ok(())
}

fn read_i32<R: Read>(
    reader: &mut R,
    source_label: &str,
    frame: Option<usize>,
    detail: &str,
) -> Result<i32, TrajectoryIoError> {
    let mut bytes = [0_u8; 4];
    read_exact(reader, &mut bytes, source_label, frame, detail)?;
    Ok(i32::from_be_bytes(bytes))
}

fn read_scalar<R: Read>(
    reader: &mut R,
    precision: TrrScalarPrecision,
    source_label: &str,
    frame: Option<usize>,
    detail: &str,
) -> Result<f64, TrajectoryIoError> {
    let mut bytes = [0_u8; 8];
    read_exact(
        reader,
        &mut bytes[..precision.bytes()],
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

fn read_exact<R: Read>(
    reader: &mut R,
    bytes: &mut [u8],
    source_label: &str,
    frame: Option<usize>,
    detail: &str,
) -> Result<(), TrajectoryIoError> {
    reader.read_exact(bytes).map_err(|error| {
        if error.kind() == std::io::ErrorKind::UnexpectedEof {
            let mut context = TrajectoryIoError::codec(
                TrajectoryIoErrorKind::TruncatedRecord,
                if frame.is_some() {
                    TrajectoryIoOperation::ReadFrame
                } else {
                    TrajectoryIoOperation::ReadHeader
                },
                Some(TrajectoryFormat::Trr),
            )
            .with_source_label(source_label)
            .with_detail(detail);
            if let Some(frame) = frame {
                context = context.with_frame(frame);
            }
            context
        } else {
            io_context(
                if frame.is_some() {
                    TrajectoryIoOperation::ReadFrame
                } else {
                    TrajectoryIoOperation::ReadHeader
                },
                Some(TrajectoryFormat::Trr),
                source_label,
                error,
            )
        }
    })
}

fn write_i32<W: Write>(
    writer: &mut W,
    value: i32,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    write_bytes(writer, &value.to_be_bytes(), source_label)
}

fn write_xdr_string<W: Write>(
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

fn write_scalar<W: Write>(
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

fn write_bytes<W: Write>(
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

fn validate_atom_count(
    atom_count: usize,
    limits: &TrajectoryIoLimits,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    if atom_count == 0 {
        return Err(resource_error(
            source_label,
            None,
            "The header declares no atoms.",
        ));
    }
    let limit = limits.max_atoms().min(i32::MAX as usize / 3);
    if atom_count > limit {
        return Err(resource_error(source_label, None, "").exceeding(
            "The atom count",
            LimitExceeded::new(Resource::Atoms, limit, atom_count),
        ));
    }
    Ok(())
}

fn header_error(source_label: &str, detail: impl Into<String>) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::InvalidHeader,
        TrajectoryIoOperation::ReadHeader,
        Some(TrajectoryFormat::Trr),
        source_label,
        detail,
    )
}

fn frame_error(source_label: &str, frame: usize, detail: impl Into<String>) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        TrajectoryIoErrorKind::InvalidFrame,
        TrajectoryIoOperation::ReadFrame,
        Some(TrajectoryFormat::Trr),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_detail(detail)
}

fn resource_error(
    source_label: &str,
    frame: Option<usize>,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    let mut context = TrajectoryIoError::codec(
        TrajectoryIoErrorKind::LimitExceeded,
        if frame.is_some() {
            TrajectoryIoOperation::ReadFrame
        } else {
            TrajectoryIoOperation::ReadHeader
        },
        Some(TrajectoryFormat::Trr),
    )
    .with_source_label(source_label)
    .with_detail(detail);
    if let Some(frame) = frame {
        context = context.with_frame(frame);
    }
    context
}

fn writer_limit(source_label: &str, detail: impl Into<String>) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::LimitExceeded,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Trr),
        source_label,
        detail,
    )
}

fn writer_unit(
    source_label: &str,
    field: &str,
    error: impl std::error::Error + Send + Sync + 'static,
) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::InconsistentMetadata,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Trr),
        source_label,
        format!("The {field} unit is incompatible."),
    )
    .with_source(error)
}

fn writer_field_error(
    source_label: &str,
    frame: usize,
    field: TrajectoryField,
) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        TrajectoryIoErrorKind::UnsupportedField,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Trr),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_field(field)
    .with_detail(format!("TRR cannot store the {field}."))
}
