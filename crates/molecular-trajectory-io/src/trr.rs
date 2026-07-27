//! Pure-Rust GROMACS TRR/XDR trajectory I/O.

use std::collections::BTreeMap;
use std::io::{Read, Seek, SeekFrom, Write};

use molecular::core::{PropMap, PropValue};
use molecular::geometry::{PeriodicCell, Point3, Vector3};
use molecular::topology::Topology;
use molecular::trajectory::{
    FrameBuffer, FrameBufferData, SeekableTrajectoryReader, TrajectoryCodecErrorContext,
    TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat, TrajectoryFrameView,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use molecular::units::{Quantity, KILOJOULE_PER_MOLE, MODEL_LENGTH_UNIT, NANOMETER, PICOSECOND};

use crate::{codec_context, io_context, TrajectoryIoLimits, TrajectoryTopologyBinding};

const TRR_MAGIC: i32 = 1993;
const TRR_VERSION: &[u8] = b"GMX_trn_file";
const TRR_LAMBDA_PROPERTY: &str = "gromacs.trr.lambda";

/// Native scalar width used by one TRR writer or frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum TrrScalarPrecision {
    Float32,
    Float64,
}

impl TrrScalarPrecision {
    const fn bytes(self) -> usize {
        match self {
            Self::Float32 => 4,
            Self::Float64 => 8,
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
    pub const fn with_lambda_policy(mut self, lambda_policy: TrrLambdaPolicy) -> Self {
        self.lambda_policy = lambda_policy;
        self
    }
}

/// TRR writer policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrrWriteOptions {
    precision: TrrScalarPrecision,
    lambda_policy: TrrLambdaPolicy,
}

impl Default for TrrWriteOptions {
    fn default() -> Self {
        Self {
            precision: TrrScalarPrecision::Float32,
            lambda_policy: TrrLambdaPolicy::FrameProperty,
        }
    }
}

impl TrrWriteOptions {
    pub const fn with_precision(mut self, precision: TrrScalarPrecision) -> Self {
        self.precision = precision;
        self
    }

    pub const fn with_lambda_policy(mut self, lambda_policy: TrrLambdaPolicy) -> Self {
        self.lambda_policy = lambda_policy;
        self
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
}

/// Sequential TRR reader retaining one seekable stream and reusable scratch.
pub struct TrrReader<R> {
    reader: R,
    binding: TrajectoryTopologyBinding,
    options: TrrReadOptions,
    limits: TrajectoryIoLimits,
    source_label: String,
    stream_start: u64,
    current_header_offset: u64,
    first_header: TrrFrameHeader,
    pending_header: Option<TrrFrameHeader>,
    positions: Vec<Point3>,
    velocities: Vec<Vector3>,
    forces: Vec<Vector3>,
    raw: Vec<u8>,
    props: PropMap,
    frame_cursor: u64,
    precision_mixed: bool,
}

impl<R: Read + Seek> TrrReader<R> {
    pub fn new(
        mut reader: R,
        binding: TrajectoryTopologyBinding,
        options: TrrReadOptions,
        limits: TrajectoryIoLimits,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryError> {
        let source_label = source_label.into();
        let atom_count = binding.topology().atom_count();
        validate_atom_count(atom_count, &limits, &source_label)?;
        let stream_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Trr),
                &source_label,
                error,
            )
        })?;
        let first_header =
            read_header(&mut reader, &limits, &source_label, Some(atom_count), false)?
                .ok_or_else(|| header_error(&source_label, "TRR stream is empty"))?;
        validate_lambda(first_header.lambda, options.lambda_policy, &source_label, 0)?;
        let mut positions = Vec::new();
        let mut velocities = Vec::new();
        let mut forces = Vec::new();
        let vector_bytes = atom_count
            .checked_mul(std::mem::size_of::<Vector3>())
            .ok_or_else(|| resource_error(&source_label, None, "TRR vector scratch overflows"))?;
        let scratch_bytes = vector_bytes.checked_mul(3).ok_or_else(|| {
            resource_error(
                &source_label,
                None,
                "TRR aggregate vector scratch overflows",
            )
        })?;
        if scratch_bytes > limits.max_scratch_bytes {
            return Err(resource_error(
                &source_label,
                None,
                "TRR vector scratch exceeds the configured limit",
            ));
        }
        positions.try_reserve_exact(atom_count).map_err(|_| {
            resource_error(
                &source_label,
                None,
                "could not reserve TRR position scratch",
            )
        })?;
        velocities.try_reserve_exact(atom_count).map_err(|_| {
            resource_error(
                &source_label,
                None,
                "could not reserve TRR velocity scratch",
            )
        })?;
        forces.try_reserve_exact(atom_count).map_err(|_| {
            resource_error(&source_label, None, "could not reserve TRR force scratch")
        })?;
        positions.resize(atom_count, Point3::new(0.0, 0.0, 0.0));
        velocities.resize(atom_count, Vector3::zero());
        forces.resize(atom_count, Vector3::zero());
        Ok(Self {
            reader,
            binding,
            options,
            limits,
            source_label,
            stream_start,
            current_header_offset: stream_start,
            first_header: first_header.clone(),
            pending_header: Some(first_header),
            positions,
            velocities,
            forces,
            raw: Vec::new(),
            props: BTreeMap::new(),
            frame_cursor: 0,
            precision_mixed: false,
        })
    }

    pub fn topology(&self) -> &Topology {
        self.binding.topology()
    }

    pub(crate) fn first_header(&self) -> &TrrFrameHeader {
        &self.first_header
    }

    pub(crate) const fn precision_mixed(&self) -> bool {
        self.precision_mixed
    }

    fn next_header(&mut self) -> Result<Option<TrrFrameHeader>, TrajectoryError> {
        if let Some(header) = self.pending_header.take() {
            return Ok(Some(header));
        }
        self.current_header_offset = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Trr),
                &self.source_label,
                error,
            )
        })?;
        let atom_count = self.binding.topology().atom_count();
        read_header(
            &mut self.reader,
            &self.limits,
            &self.source_label,
            Some(atom_count),
            true,
        )
    }

    fn parse_next(&mut self, publish: Option<&mut FrameBuffer>) -> Result<bool, TrajectoryError> {
        if self.frame_cursor >= self.limits.max_frames {
            return Err(resource_error(
                &self.source_label,
                Some(self.frame_cursor),
                "TRR frame count exceeds the configured limit",
            ));
        }
        let Some(header) = self.next_header()? else {
            return Ok(false);
        };
        self.precision_mixed |= header.precision != self.first_header.precision;
        validate_lambda(
            header.lambda,
            self.options.lambda_policy,
            &self.source_label,
            self.frame_cursor,
        )?;
        let payload_bytes = header
            .box_size
            .checked_add(header.x_size)
            .and_then(|bytes| bytes.checked_add(header.v_size))
            .and_then(|bytes| bytes.checked_add(header.f_size))
            .ok_or_else(|| {
                resource_error(
                    &self.source_label,
                    Some(self.frame_cursor),
                    "TRR payload size overflows",
                )
            })?;
        let frame_bytes = header
            .header_bytes
            .checked_add(payload_bytes as u64)
            .ok_or_else(|| {
                resource_error(
                    &self.source_label,
                    Some(self.frame_cursor),
                    "TRR frame size overflows",
                )
            })?;
        if frame_bytes > self.limits.max_frame_bytes {
            return Err(resource_error(
                &self.source_label,
                Some(self.frame_cursor),
                "TRR frame exceeds the configured byte limit",
            ));
        }

        let cell = if header.has_cell() {
            read_raw(
                &mut self.reader,
                &mut self.raw,
                header.box_size,
                &self.limits,
                &self.source_label,
                self.frame_cursor,
                "TRR box",
            )?;
            Some(decode_cell(
                &self.raw,
                header.precision,
                &self.source_label,
                self.frame_cursor,
            )?)
        } else {
            None
        };
        read_raw(
            &mut self.reader,
            &mut self.raw,
            header.x_size,
            &self.limits,
            &self.source_label,
            self.frame_cursor,
            "TRR positions",
        )?;
        decode_points(
            &self.raw,
            header.precision,
            &mut self.positions,
            &self.source_label,
            self.frame_cursor,
            "position",
        )?;
        if header.has_velocities() {
            read_raw(
                &mut self.reader,
                &mut self.raw,
                header.v_size,
                &self.limits,
                &self.source_label,
                self.frame_cursor,
                "TRR velocities",
            )?;
            decode_vectors(
                &self.raw,
                header.precision,
                &mut self.velocities,
                &self.source_label,
                self.frame_cursor,
                "velocity",
            )?;
        }
        if header.has_forces() {
            read_raw(
                &mut self.reader,
                &mut self.raw,
                header.f_size,
                &self.limits,
                &self.source_label,
                self.frame_cursor,
                "TRR forces",
            )?;
            decode_vectors(
                &self.raw,
                header.precision,
                &mut self.forces,
                &self.source_label,
                self.frame_cursor,
                "force",
            )?;
        }
        if let Some(destination) = publish {
            self.props.clear();
            if self.options.lambda_policy == TrrLambdaPolicy::FrameProperty {
                self.props.insert(
                    TRR_LAMBDA_PROPERTY.to_owned(),
                    PropValue::Float(header.lambda),
                );
            }
            let mut data = FrameBufferData::new(
                self.topology(),
                Quantity::new(self.positions.as_slice(), NANOMETER),
            )
            .with_time(Quantity::new(header.time, PICOSECOND))
            .with_step(header.step)
            .with_props(&self.props);
            if let Some(cell) = cell {
                data = data.with_cell(cell);
            }
            if header.has_velocities() {
                data = data.with_velocities(Quantity::new(
                    self.velocities.as_slice(),
                    NANOMETER / PICOSECOND,
                ));
            }
            if header.has_forces() {
                data = data.with_forces(Quantity::new(
                    self.forces.as_slice(),
                    KILOJOULE_PER_MOLE / NANOMETER,
                ));
            }
            destination.replace_from_data(data)?;
        }
        self.frame_cursor += 1;
        Ok(true)
    }

    pub fn into_indexed(mut self) -> Result<IndexedTrrReader<R>, TrajectoryError> {
        let mut offsets = Vec::new();
        loop {
            if offsets.len() >= self.limits.max_index_entries {
                return Err(resource_error(
                    &self.source_label,
                    Some(offsets.len() as u64),
                    "TRR index entry limit exceeded",
                ));
            }
            let offset = if self.pending_header.is_some() {
                self.current_header_offset
            } else {
                self.reader.stream_position().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Index,
                        Some(TrajectoryFormat::Trr),
                        &self.source_label,
                        error,
                    )
                })?
            };
            if !self.parse_next(None)? {
                break;
            }
            offsets.try_reserve(1).map_err(|_| {
                resource_error(
                    &self.source_label,
                    Some(offsets.len() as u64),
                    "could not grow TRR index",
                )
            })?;
            offsets.push(offset);
            if offsets
                .len()
                .checked_mul(std::mem::size_of::<u64>())
                .is_none_or(|bytes| bytes > self.limits.max_index_bytes)
            {
                return Err(resource_error(
                    &self.source_label,
                    Some(offsets.len() as u64),
                    "TRR index byte limit exceeded",
                ));
            }
        }
        self.rewind()?;
        Ok(IndexedTrrReader {
            inner: self,
            offsets,
        })
    }

    fn rewind(&mut self) -> Result<(), TrajectoryError> {
        self.reader
            .seek(SeekFrom::Start(self.stream_start))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Trr),
                    &self.source_label,
                    error,
                )
            })?;
        self.current_header_offset = self.stream_start;
        let atom_count = self.binding.topology().atom_count();
        self.pending_header = read_header(
            &mut self.reader,
            &self.limits,
            &self.source_label,
            Some(atom_count),
            false,
        )?;
        self.frame_cursor = 0;
        Ok(())
    }
}

impl<R: Read + Seek> TrajectoryReader for TrrReader<R> {
    fn topology(&self) -> &Topology {
        self.topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        self.parse_next(Some(destination))
    }
}

/// Fully verified indexed TRR reader.
pub struct IndexedTrrReader<R> {
    inner: TrrReader<R>,
    offsets: Vec<u64>,
}

impl<R: Read + Seek> IndexedTrrReader<R> {
    pub fn topology(&self) -> &Topology {
        self.inner.topology()
    }

    pub(crate) fn first_header(&self) -> &TrrFrameHeader {
        self.inner.first_header()
    }

    pub(crate) const fn precision_mixed(&self) -> bool {
        self.inner.precision_mixed()
    }
}

impl<R: Read + Seek> TrajectoryReader for IndexedTrrReader<R> {
    fn topology(&self) -> &Topology {
        self.topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        self.inner.read_next(destination)
    }
}

impl<R: Read + Seek> SeekableTrajectoryReader for IndexedTrrReader<R> {
    fn frame_count(&self) -> Option<u64> {
        Some(self.offsets.len() as u64)
    }

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        let offset = self
            .offsets
            .get(usize::try_from(index).map_err(|_| TrajectoryError::FrameIndexOutOfRange(index))?)
            .copied()
            .ok_or(TrajectoryError::FrameIndexOutOfRange(index))?;
        let saved_offset = self.inner.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Trr),
                &self.inner.source_label,
                error,
            )
        })?;
        let saved_cursor = self.inner.frame_cursor;
        let saved_pending = self.inner.pending_header.take();
        let saved_header_offset = self.inner.current_header_offset;
        self.inner
            .reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Trr),
                    &self.inner.source_label,
                    error,
                )
            })?;
        self.inner.current_header_offset = offset;
        self.inner.frame_cursor = index;
        let result = self.inner.parse_next(Some(destination)).and_then(|read| {
            if read {
                Ok(())
            } else {
                Err(TrajectoryError::FrameIndexOutOfRange(index))
            }
        });
        let restore = self
            .inner
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Trr),
                    &self.inner.source_label,
                    error,
                )
            });
        self.inner.frame_cursor = saved_cursor;
        self.inner.pending_header = saved_pending;
        self.inner.current_header_offset = saved_header_offset;
        result.and(restore.map(|_| ()))
    }
}

/// Pure-Rust TRR writer.
pub struct TrrWriter<W> {
    writer: W,
    topology: Topology,
    options: TrrWriteOptions,
    source_label: String,
    frame_count: u64,
    raw: Vec<u8>,
}

impl<W: Write> TrrWriter<W> {
    pub fn new(
        writer: W,
        topology: Topology,
        options: TrrWriteOptions,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryError> {
        let source_label = source_label.into();
        if topology.atom_count() == 0 || topology.atom_count() > i32::MAX as usize / 3 {
            return Err(codec_context(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Trr),
                &source_label,
                "TRR atom count must fit the signed 32-bit format",
            ));
        }
        Ok(Self {
            writer,
            topology,
            options,
            source_label,
            frame_count: 0,
            raw: Vec::new(),
        })
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn writer(&self) -> &W {
        &self.writer
    }

    pub fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }

    pub fn finish(mut self) -> Result<W, TrajectoryError> {
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

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if frame.topology().identity() != self.topology.identity() {
            return Err(TrajectoryError::TopologyIdentityMismatch);
        }
        if frame.observation().is_some() {
            return Err(TrajectoryError::UnsupportedField("observation"));
        }
        let step = frame.step().ok_or_else(|| {
            codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Trr),
                &self.source_label,
                "TRR requires an explicit frame step",
            )
        })?;
        let step = i32::try_from(step).map_err(|_| {
            codec_context(
                TrajectoryCodecErrorKind::NegativeOrUnrepresentableStep,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Trr),
                &self.source_label,
                "TRR step exceeds signed 32-bit capacity",
            )
        })?;
        let time = frame.time().ok_or_else(|| {
            codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Trr),
                &self.source_label,
                "TRR requires explicit frame time",
            )
        })?;
        let time = time.value_in(PICOSECOND).map_err(|error| {
            codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Trr),
                &self.source_label,
                format!("TRR time unit is incompatible: {error}"),
            )
        })?;
        let lambda = writer_lambda(
            frame.props(),
            self.options.lambda_policy,
            &self.source_label,
        )?;
        let scalar_bytes = self.options.precision.bytes();
        validate_representable(time, self.options.precision, &self.source_label, "time")?;
        validate_representable(lambda, self.options.precision, &self.source_label, "lambda")?;
        let atom_scalars =
            self.topology.atom_count().checked_mul(3).ok_or_else(|| {
                writer_limit(&self.source_label, "TRR atom scalar count overflows")
            })?;
        let vector_bytes = atom_scalars
            .checked_mul(scalar_bytes)
            .ok_or_else(|| writer_limit(&self.source_label, "TRR vector byte size overflows"))?;
        let box_size = if frame.configuration().cell().is_some() {
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
        write_i32(&mut self.writer, TRR_MAGIC, &self.source_label)?;
        write_i32(
            &mut self.writer,
            i32::try_from(TRR_VERSION.len() + 1).expect("TRR version length fits i32"),
            &self.source_label,
        )?;
        write_xdr_string(&mut self.writer, TRR_VERSION, &self.source_label)?;
        for size in [0, 0, box_size, 0, 0, 0, 0, vector_bytes, v_size, f_size] {
            write_i32(
                &mut self.writer,
                i32::try_from(size)
                    .map_err(|_| writer_limit(&self.source_label, "TRR block size exceeds i32"))?,
                &self.source_label,
            )?;
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
        if let Some(cell) = frame.configuration().cell().copied() {
            if cell.periodic_axes() != [true; 3] {
                return Err(codec_context(
                    TrajectoryCodecErrorKind::UnsupportedVariant,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Trr),
                    &self.source_label,
                    "TRR cell requires all three periodic axes",
                ));
            }
            let factor = MODEL_LENGTH_UNIT
                .conversion_factor_to(NANOMETER)
                .map_err(|error| writer_unit(&self.source_label, "cell", error))?;
            let vectors = cell.vectors().into_value();
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
        let positions = frame.configuration().positions().values();
        let factor = MODEL_LENGTH_UNIT
            .conversion_factor_to(NANOMETER)
            .map_err(|error| writer_unit(&self.source_label, "positions", error))?;
        encode_points_to_raw(
            &mut self.raw,
            positions.value(),
            factor,
            self.options.precision,
            &self.source_label,
        )?;
        write_bytes(&mut self.writer, &self.raw, &self.source_label)?;
        if let Some(velocities) = frame.velocities() {
            let factor = velocities
                .unit()
                .conversion_factor_to(NANOMETER / PICOSECOND)
                .map_err(|error| writer_unit(&self.source_label, "velocities", error))?;
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
        if let Some(forces) = frame.forces() {
            let factor = forces
                .unit()
                .conversion_factor_to(KILOJOULE_PER_MOLE / NANOMETER)
                .map_err(|error| writer_unit(&self.source_label, "forces", error))?;
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
        self.frame_count = self
            .frame_count
            .checked_add(1)
            .ok_or_else(|| writer_limit(&self.source_label, "TRR frame count overflows"))?;
        Ok(())
    }
}

fn read_header<R: Read>(
    reader: &mut R,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    expected_atoms: Option<usize>,
    clean_eof: bool,
) -> Result<Option<TrrFrameHeader>, TrajectoryError> {
    let mut magic = [0_u8; 4];
    match reader.read(&mut magic[..1]) {
        Ok(0) if clean_eof => return Ok(None),
        Ok(0) => return Err(header_error(source_label, "missing TRR magic")),
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
        "partial TRR magic",
    )?;
    if i32::from_be_bytes(magic) != TRR_MAGIC {
        return Err(header_error(source_label, "TRR magic is not 1993"));
    }
    let declared_version_bytes =
        read_i32(reader, source_label, None, "TRR declared version storage")?;
    if declared_version_bytes
        != i32::try_from(TRR_VERSION.len() + 1).expect("TRR version length fits i32")
    {
        return Err(codec_context(
            TrajectoryCodecErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "TRR declared version storage is not the GMX_trn_file profile",
        ));
    }
    let (version, version_bytes) = read_xdr_string(reader, limits, source_label)?;
    if version != TRR_VERSION {
        return Err(codec_context(
            TrajectoryCodecErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "TRR version string is not GMX_trn_file",
        ));
    }
    let mut sizes = [0_usize; 10];
    for size in &mut sizes {
        let value = read_i32(reader, source_label, None, "TRR size field")?;
        *size = usize::try_from(value)
            .map_err(|_| header_error(source_label, "TRR block size is negative"))?;
        if *size as u64 > limits.max_record_bytes {
            return Err(resource_error(
                source_label,
                None,
                "TRR block exceeds the configured record limit",
            ));
        }
    }
    let atom_count = usize::try_from(read_i32(reader, source_label, None, "TRR atom count")?)
        .map_err(|_| header_error(source_label, "TRR atom count is negative"))?;
    validate_atom_count(atom_count, limits, source_label)?;
    if let Some(expected) = expected_atoms {
        if atom_count != expected {
            return Err(TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Trr),
            )
            .with_source_label(source_label)
            .with_counts(expected as u64, atom_count as u64)
            .into());
        }
    }
    let [ir, energy, box_size, virial, pressure, topology, symbols, x_size, v_size, f_size] = sizes;
    if [ir, energy, virial, pressure, topology, symbols]
        .into_iter()
        .any(|size| size != 0)
    {
        return Err(codec_context(
            TrajectoryCodecErrorKind::UnsupportedField,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "TRR inputrec/energy/virial/pressure/topology/symbol blocks cannot be preserved",
        ));
    }
    if x_size == 0 {
        return Err(codec_context(
            TrajectoryCodecErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "TRR frames without positions cannot populate Molecular FrameBuffer",
        ));
    }
    let precision = infer_precision(atom_count, box_size, x_size, v_size, f_size, source_label)?;
    let expected_vectors = atom_count
        .checked_mul(3)
        .and_then(|count| count.checked_mul(precision.bytes()))
        .ok_or_else(|| resource_error(source_label, None, "TRR vector size overflows"))?;
    if x_size != expected_vectors
        || (v_size != 0 && v_size != expected_vectors)
        || (f_size != 0 && f_size != expected_vectors)
        || (box_size != 0 && box_size != 9 * precision.bytes())
    {
        return Err(codec_context(
            TrajectoryCodecErrorKind::InvalidRecordLength,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "TRR block sizes do not match atom count and scalar precision",
        ));
    }
    let step_raw = read_i32(reader, source_label, None, "TRR step")?;
    let step = u64::try_from(step_raw).map_err(|_| {
        codec_context(
            TrajectoryCodecErrorKind::NegativeOrUnrepresentableStep,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "TRR step is negative",
        )
    })?;
    let nre = read_i32(reader, source_label, None, "TRR nre")?;
    if nre != 0 {
        return Err(codec_context(
            TrajectoryCodecErrorKind::UnsupportedField,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "TRR nonzero nre metadata cannot be preserved",
        ));
    }
    let time = read_scalar(reader, precision, source_label, None, "TRR time")?;
    let lambda = read_scalar(reader, precision, source_label, None, "TRR lambda")?;
    if !time.is_finite() || !lambda.is_finite() {
        return Err(header_error(
            source_label,
            "TRR time and lambda must be finite",
        ));
    }
    let header_bytes = 4_u64
        .checked_add(4)
        .and_then(|bytes| bytes.checked_add(version_bytes))
        .and_then(|bytes| bytes.checked_add(11 * 4 + 2 * 4))
        .and_then(|bytes| bytes.checked_add((2 * precision.bytes()) as u64))
        .ok_or_else(|| resource_error(source_label, None, "TRR header size overflows"))?;
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
) -> Result<TrrScalarPrecision, TrajectoryError> {
    let candidates = [
        (box_size, 9_usize),
        (x_size, atom_count.saturating_mul(3)),
        (v_size, atom_count.saturating_mul(3)),
        (f_size, atom_count.saturating_mul(3)),
    ];
    let (size, scalars) = candidates
        .into_iter()
        .find(|(size, _)| *size != 0)
        .ok_or_else(|| header_error(source_label, "TRR has no precision-bearing block"))?;
    let width = size
        .checked_div(scalars)
        .filter(|_| size % scalars == 0)
        .ok_or_else(|| header_error(source_label, "TRR scalar width is inconsistent"))?;
    match width {
        4 => Ok(TrrScalarPrecision::Float32),
        8 => Ok(TrrScalarPrecision::Float64),
        _ => Err(codec_context(
            TrajectoryCodecErrorKind::InvalidPrecision,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Trr),
            source_label,
            "TRR scalar width must be f32 or f64",
        )),
    }
}

fn read_xdr_string<R: Read>(
    reader: &mut R,
    limits: &TrajectoryIoLimits,
    source_label: &str,
) -> Result<(Vec<u8>, u64), TrajectoryError> {
    let length = usize::try_from(read_i32(reader, source_label, None, "TRR version length")?)
        .map_err(|_| header_error(source_label, "TRR version length is negative"))?;
    if length > 256 || length > limits.max_text_line_bytes || length > limits.max_scratch_bytes {
        return Err(resource_error(
            source_label,
            None,
            "TRR version string exceeds configured limits",
        ));
    }
    let padded = length
        .checked_add(3)
        .map(|value| value & !3)
        .ok_or_else(|| resource_error(source_label, None, "TRR string padding overflows"))?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(padded)
        .map_err(|_| resource_error(source_label, None, "could not reserve TRR version scratch"))?;
    bytes.resize(padded, 0);
    read_exact(
        reader,
        &mut bytes,
        source_label,
        None,
        "truncated TRR version string",
    )?;
    if bytes[length..].iter().any(|byte| *byte != 0) {
        return Err(header_error(
            source_label,
            "TRR XDR string padding is not zero",
        ));
    }
    bytes.truncate(length);
    Ok((bytes, 4 + padded as u64))
}

fn read_raw<R: Read>(
    reader: &mut R,
    raw: &mut Vec<u8>,
    length: usize,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    frame: u64,
    field: &str,
) -> Result<(), TrajectoryError> {
    if length > limits.max_scratch_bytes || length as u64 > limits.max_record_bytes {
        return Err(resource_error(
            source_label,
            Some(frame),
            format!("{field} exceeds configured scratch or record limits"),
        ));
    }
    if raw.capacity() < length {
        raw.try_reserve_exact(length - raw.capacity())
            .map_err(|_| {
                resource_error(
                    source_label,
                    Some(frame),
                    format!("could not reserve {field} scratch"),
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

fn decode_points(
    raw: &[u8],
    precision: TrrScalarPrecision,
    destination: &mut [Point3],
    source_label: &str,
    frame: u64,
    field: &str,
) -> Result<(), TrajectoryError> {
    for (point, values) in destination
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
                format!("TRR {field} contains a non-finite value"),
            ));
        }
        *point = Point3::new(x, y, z);
    }
    Ok(())
}

fn decode_vectors(
    raw: &[u8],
    precision: TrrScalarPrecision,
    destination: &mut [Vector3],
    source_label: &str,
    frame: u64,
    field: &str,
) -> Result<(), TrajectoryError> {
    for (vector, values) in destination
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
                format!("TRR {field} contains a non-finite value"),
            ));
        }
        *vector = Vector3::new(x, y, z);
    }
    Ok(())
}

fn decode_cell(
    raw: &[u8],
    precision: TrrScalarPrecision,
    source_label: &str,
    frame: u64,
) -> Result<PeriodicCell, TrajectoryError> {
    let mut vectors = [Vector3::zero(); 3];
    decode_vectors(raw, precision, &mut vectors, source_label, frame, "box")?;
    PeriodicCell::new(Quantity::new(vectors, NANOMETER), [true; 3]).map_err(|error| {
        frame_error(
            source_label,
            frame,
            format!("TRR periodic cell is invalid: {error}"),
        )
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
) -> Result<(), TrajectoryError> {
    raw.clear();
    let bytes = points
        .len()
        .checked_mul(3 * precision.bytes())
        .ok_or_else(|| writer_limit(source_label, "TRR position bytes overflow"))?;
    raw.try_reserve(bytes.saturating_sub(raw.capacity()))
        .map_err(|_| writer_limit(source_label, "could not reserve TRR writer scratch"))?;
    for point in points {
        for value in [point.x * factor, point.y * factor, point.z * factor] {
            push_scalar(raw, value, precision, source_label, "position")?;
        }
    }
    Ok(())
}

fn encode_vectors_to_raw(
    raw: &mut Vec<u8>,
    vectors: &[Vector3],
    factor: f64,
    precision: TrrScalarPrecision,
    source_label: &str,
    field: &str,
) -> Result<(), TrajectoryError> {
    raw.clear();
    let bytes = vectors
        .len()
        .checked_mul(3 * precision.bytes())
        .ok_or_else(|| writer_limit(source_label, "TRR vector bytes overflow"))?;
    raw.try_reserve(bytes.saturating_sub(raw.capacity()))
        .map_err(|_| writer_limit(source_label, "could not reserve TRR writer scratch"))?;
    for vector in vectors {
        for value in [vector.x * factor, vector.y * factor, vector.z * factor] {
            push_scalar(raw, value, precision, source_label, field)?;
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
) -> Result<(), TrajectoryError> {
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
) -> Result<(), TrajectoryError> {
    if !value.is_finite()
        || (precision == TrrScalarPrecision::Float32 && !(value as f32).is_finite())
    {
        return Err(codec_context(
            TrajectoryCodecErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Trr),
            source_label,
            format!("TRR {field} is not representable at the selected precision"),
        ));
    }
    Ok(())
}

fn writer_lambda(
    props: &PropMap,
    policy: TrrLambdaPolicy,
    source_label: &str,
) -> Result<f64, TrajectoryError> {
    match policy {
        TrrLambdaPolicy::FrameProperty => {
            if props.len() != 1 {
                return Err(codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Trr),
                    source_label,
                    "TRR FrameProperty policy requires exactly gromacs.trr.lambda",
                ));
            }
            match props.get(TRR_LAMBDA_PROPERTY) {
                Some(PropValue::Float(value)) if value.is_finite() => Ok(*value),
                _ => Err(codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Trr),
                    source_label,
                    "TRR lambda property must be one finite float",
                )),
            }
        }
        TrrLambdaPolicy::RequireZero => {
            if !props.is_empty() {
                return Err(TrajectoryError::UnsupportedField("properties"));
            }
            Ok(0.0)
        }
    }
}

fn validate_lambda(
    lambda: f64,
    policy: TrrLambdaPolicy,
    source_label: &str,
    frame: u64,
) -> Result<(), TrajectoryError> {
    if policy == TrrLambdaPolicy::RequireZero && lambda != 0.0 {
        return Err(TrajectoryCodecErrorContext::new(
            TrajectoryCodecErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::ReadFrame,
            Some(TrajectoryFormat::Trr),
        )
        .with_source_label(source_label)
        .with_frame(frame)
        .with_detail("TRR lambda is nonzero under RequireZero policy")
        .into());
    }
    Ok(())
}

fn read_i32<R: Read>(
    reader: &mut R,
    source_label: &str,
    frame: Option<u64>,
    detail: &str,
) -> Result<i32, TrajectoryError> {
    let mut bytes = [0_u8; 4];
    read_exact(reader, &mut bytes, source_label, frame, detail)?;
    Ok(i32::from_be_bytes(bytes))
}

fn read_scalar<R: Read>(
    reader: &mut R,
    precision: TrrScalarPrecision,
    source_label: &str,
    frame: Option<u64>,
    detail: &str,
) -> Result<f64, TrajectoryError> {
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
    frame: Option<u64>,
    detail: &str,
) -> Result<(), TrajectoryError> {
    reader.read_exact(bytes).map_err(|error| {
        if error.kind() == std::io::ErrorKind::UnexpectedEof {
            let mut context = TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::TruncatedRecord,
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
            context.into()
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
) -> Result<(), TrajectoryError> {
    write_bytes(writer, &value.to_be_bytes(), source_label)
}

fn write_xdr_string<W: Write>(
    writer: &mut W,
    value: &[u8],
    source_label: &str,
) -> Result<(), TrajectoryError> {
    write_i32(
        writer,
        i32::try_from(value.len())
            .map_err(|_| writer_limit(source_label, "TRR version is too long"))?,
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
) -> Result<(), TrajectoryError> {
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
) -> Result<(), TrajectoryError> {
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
) -> Result<(), TrajectoryError> {
    if atom_count == 0 || atom_count > limits.max_atoms || atom_count > i32::MAX as usize / 3 {
        return Err(resource_error(
            source_label,
            None,
            "TRR atom count is zero or exceeds configured/profile limits",
        ));
    }
    Ok(())
}

fn header_error(source_label: &str, detail: impl Into<String>) -> TrajectoryError {
    codec_context(
        TrajectoryCodecErrorKind::InvalidHeader,
        TrajectoryIoOperation::ReadHeader,
        Some(TrajectoryFormat::Trr),
        source_label,
        detail,
    )
}

fn frame_error(source_label: &str, frame: u64, detail: impl Into<String>) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::InvalidFrame,
        TrajectoryIoOperation::ReadFrame,
        Some(TrajectoryFormat::Trr),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_detail(detail)
    .into()
}

fn resource_error(
    source_label: &str,
    frame: Option<u64>,
    detail: impl Into<String>,
) -> TrajectoryError {
    let mut context = TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::ResourceLimitExceeded,
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
    context.into()
}

fn writer_limit(source_label: &str, detail: impl Into<String>) -> TrajectoryError {
    codec_context(
        TrajectoryCodecErrorKind::ResourceLimitExceeded,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Trr),
        source_label,
        detail,
    )
}

fn writer_unit(source_label: &str, field: &str, error: impl std::fmt::Display) -> TrajectoryError {
    codec_context(
        TrajectoryCodecErrorKind::InconsistentMetadata,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Trr),
        source_label,
        format!("TRR {field} unit is incompatible: {error}"),
    )
}
