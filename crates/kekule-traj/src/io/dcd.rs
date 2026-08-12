//! Strict common-profile CHARMM/NAMD/OpenMM DCD trajectory I/O.

use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::Arc;

use crate::{
    FrameBuffer, FrameBufferData, SeekableTrajectoryReader, TrajectoryCodecErrorContext,
    TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat, TrajectoryFrameView,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::topology::Topology;
use kekule::units::{Quantity, Unit, ANGSTROM, MODEL_LENGTH_UNIT, MODEL_TIME_UNIT};

use super::{
    codec_context, frame_offset_context, io_context, probe_seekable_eof, projected_index_limit,
    require_nonempty_writer, reserve_index_for_push, TrajectoryIoLimits, TrajectoryTopologyBinding,
};

const HEADER_BYTES: usize = 84;
const CELL_BYTES: usize = 48;

/// Byte order used by DCD record markers and numeric payloads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DcdEndian {
    Little,
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
    HeaderDelta { unit: Unit },
}

/// DCD reader policy.
#[derive(Debug, Clone, Copy, PartialEq)]
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
    pub const fn with_time_policy(mut self, time_policy: DcdTimePolicy) -> Self {
        self.time_policy = time_policy;
        self
    }

    pub const fn time_policy(self) -> DcdTimePolicy {
        self.time_policy
    }
}

/// Canonical DCD writer policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DcdWriteOptions {
    endian: DcdEndian,
    write_cell: bool,
    start_step: u64,
    step_interval: u64,
    delta: f32,
    time_policy: DcdTimePolicy,
}

impl Default for DcdWriteOptions {
    fn default() -> Self {
        Self {
            endian: DcdEndian::Little,
            write_cell: false,
            start_step: 0,
            step_interval: 1,
            delta: 1.0,
            time_policy: DcdTimePolicy::Absent,
        }
    }
}

impl DcdWriteOptions {
    pub const fn with_endian(mut self, endian: DcdEndian) -> Self {
        self.endian = endian;
        self
    }

    pub const fn with_cells(mut self, write_cell: bool) -> Self {
        self.write_cell = write_cell;
        self
    }

    pub const fn with_step_sequence(mut self, start_step: u64, step_interval: u64) -> Self {
        self.start_step = start_step;
        self.step_interval = step_interval;
        self
    }

    pub const fn with_header_delta(mut self, delta: f32, time_policy: DcdTimePolicy) -> Self {
        self.delta = delta;
        self.time_policy = time_policy;
        self
    }
}

#[derive(Debug, Clone)]
pub(crate) struct DcdHeader {
    endian: DcdEndian,
    declared_frames: u64,
    start_step: u64,
    step_interval: u64,
    delta: f32,
    fixed_count: usize,
    atom_count: usize,
    has_cell: bool,
    charmm_version: i32,
    data_start: u64,
}

impl DcdHeader {
    pub(crate) const fn declared_frames(&self) -> u64 {
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

/// Sequential DCD reader over one seekable stream.
pub struct DcdReader<R> {
    reader: R,
    binding: TrajectoryTopologyBinding,
    options: DcdReadOptions,
    limits: TrajectoryIoLimits,
    source_label: String,
    header: DcdHeader,
    free_indices: Vec<usize>,
    fixed_indices: Vec<usize>,
    fixed_reference: Vec<Point3>,
    positions: Vec<Point3>,
    record: Vec<u8>,
    frame_cursor: u64,
}

#[derive(Debug, Clone, Copy)]
struct DcdDecodedFrame {
    cell: Option<PeriodicCell>,
    step: u64,
    time: Option<f64>,
}

impl<R: Read + Seek> DcdReader<R> {
    pub fn new(
        mut reader: R,
        binding: TrajectoryTopologyBinding,
        options: DcdReadOptions,
        limits: TrajectoryIoLimits,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryError> {
        let source_label = source_label.into();
        validate_time_policy(options.time_policy, &source_label)?;
        let atom_count = binding.topology().atom_count();
        validate_atom_count(atom_count, &limits, &source_label)?;
        let start = reader
            .stream_position()
            .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &source_label, error))?;
        let mut marker = [0_u8; 4];
        read_exact_required(
            &mut reader,
            &mut marker,
            TrajectoryIoOperation::ReadHeader,
            &source_label,
            None,
            "DCD header marker",
        )?;
        let endian = if u32::from_le_bytes(marker) == HEADER_BYTES as u32 {
            DcdEndian::Little
        } else if u32::from_be_bytes(marker) == HEADER_BYTES as u32 {
            DcdEndian::Big
        } else {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidHeader,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                &source_label,
                "DCD first record marker is not 84 in either byte order",
            ));
        };
        reader
            .seek(SeekFrom::Start(start))
            .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &source_label, error))?;

        let mut record = Vec::new();
        read_record(
            &mut reader,
            endian,
            &mut record,
            &limits,
            &source_label,
            TrajectoryIoOperation::ReadHeader,
            None,
            false,
        )?;
        if record.len() != HEADER_BYTES || &record[..4] != b"CORD" {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidHeader,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                &source_label,
                "DCD header must be an 84-byte CORD record",
            ));
        }
        let controls = parse_controls(&record[4..], endian, &source_label)?;
        if controls[11] != 0 || controls[12] != 0 {
            return Err(codec_context(
                TrajectoryCodecErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                &source_label,
                "DCD 4D or charge coordinate records are not supported",
            ));
        }
        let declared_frames = nonnegative_u64(controls[0], "NSET", &source_label)?;
        let start_step = nonnegative_u64(controls[1], "ISTART", &source_label)?;
        let step_interval = positive_u64(controls[2], "NSAVC", &source_label)?;
        let fixed_count = nonnegative_usize(controls[8], "NAMNF", &source_label)?;
        let delta = endian.f32(record[40..44].try_into().expect("header slice"));
        if !delta.is_finite() || delta < 0.0 {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                &source_label,
                "DCD DELTA must be finite and nonnegative",
            ));
        }
        let has_cell = controls[10] != 0;
        let charmm_version = controls[19];
        if charmm_version <= 0 {
            return Err(codec_context(
                TrajectoryCodecErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                &source_label,
                "only the common CHARMM/NAMD/OpenMM DCD profile is supported",
            ));
        }

        read_record(
            &mut reader,
            endian,
            &mut record,
            &limits,
            &source_label,
            TrajectoryIoOperation::ReadHeader,
            None,
            false,
        )?;
        validate_title_record(&record, endian, &source_label)?;
        read_record(
            &mut reader,
            endian,
            &mut record,
            &limits,
            &source_label,
            TrajectoryIoOperation::ReadHeader,
            None,
            false,
        )?;
        if record.len() != 4 {
            return Err(header_error(
                &source_label,
                "DCD atom-count record is not 4 bytes",
            ));
        }
        let file_atoms = nonnegative_usize(
            endian.i32(record[..4].try_into().expect("atom-count slice")),
            "NATOM",
            &source_label,
        )?;
        validate_atom_count(file_atoms, &limits, &source_label)?;
        if file_atoms != atom_count {
            return Err(TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
            )
            .with_source_label(&source_label)
            .with_counts(atom_count as u64, file_atoms as u64)
            .into());
        }
        if fixed_count > atom_count {
            return Err(header_error(
                &source_label,
                "DCD fixed-atom count exceeds atom count",
            ));
        }
        let free_count = atom_count - fixed_count;
        let mut free_indices = Vec::new();
        if fixed_count > 0 {
            read_record(
                &mut reader,
                endian,
                &mut record,
                &limits,
                &source_label,
                TrajectoryIoOperation::ReadHeader,
                None,
                false,
            )?;
            let expected = checked_coordinate_bytes(free_count, &source_label)?;
            if record.len() != expected {
                return Err(header_error(
                    &source_label,
                    "DCD free-atom index record has the wrong size",
                ));
            }
            free_indices.try_reserve_exact(free_count).map_err(|_| {
                resource_error(
                    TrajectoryIoOperation::ReadHeader,
                    &source_label,
                    None,
                    "could not reserve DCD free-atom indices",
                )
            })?;
            let mut seen = Vec::new();
            seen.try_reserve_exact(atom_count).map_err(|_| {
                resource_error(
                    TrajectoryIoOperation::ReadHeader,
                    &source_label,
                    None,
                    "could not reserve DCD index-validation scratch",
                )
            })?;
            seen.resize(atom_count, false);
            for chunk in record.chunks_exact(4) {
                let one_based = endian.i32(chunk.try_into().expect("index chunk"));
                if one_based <= 0 {
                    return Err(header_error(
                        &source_label,
                        "DCD free-atom indices must be positive and one-based",
                    ));
                }
                let index = usize::try_from(one_based - 1).map_err(|_| {
                    header_error(&source_label, "DCD free-atom index does not fit usize")
                })?;
                if index >= atom_count || seen[index] {
                    return Err(header_error(
                        &source_label,
                        "DCD free-atom indices are out of range or duplicated",
                    ));
                }
                seen[index] = true;
                free_indices.push(index);
            }
        } else {
            free_indices.try_reserve_exact(atom_count).map_err(|_| {
                resource_error(
                    TrajectoryIoOperation::ReadHeader,
                    &source_label,
                    None,
                    "could not reserve DCD atom indices",
                )
            })?;
            free_indices.extend(0..atom_count);
        }
        let mut free_mask = Vec::new();
        free_mask.try_reserve_exact(atom_count).map_err(|_| {
            resource_error(
                TrajectoryIoOperation::ReadHeader,
                &source_label,
                None,
                "could not reserve DCD fixed-atom validation scratch",
            )
        })?;
        free_mask.resize(atom_count, false);
        for &index in &free_indices {
            free_mask[index] = true;
        }
        let mut fixed_indices = Vec::new();
        fixed_indices.try_reserve_exact(fixed_count).map_err(|_| {
            resource_error(
                TrajectoryIoOperation::ReadHeader,
                &source_label,
                None,
                "could not reserve DCD fixed-atom indices",
            )
        })?;
        fixed_indices.extend(
            free_mask
                .iter()
                .enumerate()
                .filter_map(|(index, free)| (!free).then_some(index)),
        );
        let scratch_bytes = atom_count
            .checked_mul(std::mem::size_of::<Point3>())
            .ok_or_else(|| {
                resource_error(
                    TrajectoryIoOperation::Open,
                    &source_label,
                    None,
                    "DCD coordinate scratch size overflows",
                )
            })?;
        if scratch_bytes > limits.max_scratch_bytes {
            return Err(resource_error(
                TrajectoryIoOperation::Open,
                &source_label,
                None,
                "DCD coordinate scratch exceeds the configured limit",
            ));
        }
        let mut positions = Vec::new();
        positions.try_reserve_exact(atom_count).map_err(|_| {
            resource_error(
                TrajectoryIoOperation::Open,
                &source_label,
                None,
                "could not reserve DCD coordinate scratch",
            )
        })?;
        positions.resize(atom_count, Point3::new(0.0, 0.0, 0.0));
        let data_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                &source_label,
                error,
            )
        })?;
        let header = DcdHeader {
            endian,
            declared_frames,
            start_step,
            step_interval,
            delta,
            fixed_count,
            atom_count,
            has_cell,
            charmm_version,
            data_start,
        };
        Ok(Self {
            reader,
            binding,
            options,
            limits,
            source_label,
            header,
            free_indices,
            fixed_indices,
            fixed_reference: Vec::new(),
            positions,
            record,
            frame_cursor: 0,
        })
    }

    pub fn topology(&self) -> &Topology {
        self.binding.topology()
    }

    pub(crate) fn header(&self) -> &DcdHeader {
        &self.header
    }

    fn early_eof_error(&self, frame_start: u64) -> TrajectoryError {
        TrajectoryCodecErrorContext::new(
            TrajectoryCodecErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::ReadFrame,
            Some(TrajectoryFormat::Dcd),
        )
        .with_source_label(&self.source_label)
        .with_frame(self.frame_cursor)
        .with_byte_offset(frame_start)
        .with_counts(self.header.declared_frames, self.frame_cursor)
        .with_detail("DCD ended before its declared NSET count")
        .into()
    }

    fn parse_next(
        &mut self,
        capture_fixed_reference: bool,
    ) -> Result<Option<DcdDecodedFrame>, TrajectoryError> {
        let frame_start = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                error,
            )
        })?;
        if self.frame_cursor >= self.header.declared_frames {
            if probe_seekable_eof(
                &mut self.reader,
                TrajectoryIoOperation::ReadFrame,
                TrajectoryFormat::Dcd,
                &self.source_label,
            )? {
                return Ok(None);
            }
            return Err(TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Dcd),
            )
            .with_source_label(&self.source_label)
            .with_frame(self.frame_cursor)
            .with_byte_offset(frame_start)
            .with_counts(
                self.header.declared_frames,
                self.frame_cursor.saturating_add(1),
            )
            .with_detail("DCD contains a frame beyond its declared NSET count")
            .into());
        }
        if self.frame_cursor >= self.limits.max_frames {
            if probe_seekable_eof(
                &mut self.reader,
                TrajectoryIoOperation::ReadFrame,
                TrajectoryFormat::Dcd,
                &self.source_label,
            )? {
                return Err(self.early_eof_error(frame_start));
            }
            return Err(resource_error(
                TrajectoryIoOperation::ReadFrame,
                &self.source_label,
                Some(self.frame_cursor),
                "DCD frame count exceeds the configured limit",
            ));
        }
        let mut cell = None;
        if self.header.has_cell {
            if !read_record(
                &mut self.reader,
                self.header.endian,
                &mut self.record,
                &self.limits,
                &self.source_label,
                TrajectoryIoOperation::ReadFrame,
                Some(self.frame_cursor),
                true,
            )? {
                return Err(self.early_eof_error(frame_start));
            }
            if self.record.len() != CELL_BYTES {
                return Err(frame_error(
                    TrajectoryCodecErrorKind::InvalidFrame,
                    &self.source_label,
                    self.frame_cursor,
                    "DCD unit-cell record must contain six f64 values",
                ));
            }
            cell = Some(decode_cell(
                &self.record,
                self.header.endian,
                &self.source_label,
                self.frame_cursor,
            )?);
        }

        let coordinate_indices = if self.frame_cursor == 0 || self.header.fixed_count == 0 {
            None
        } else {
            Some(self.free_indices.as_slice())
        };
        let coordinate_count = coordinate_indices.map_or(self.header.atom_count, <[usize]>::len);
        let coordinate_bytes = checked_coordinate_bytes(coordinate_count, &self.source_label)?;
        if coordinate_bytes as u64 > self.limits.max_frame_bytes {
            return Err(resource_error(
                TrajectoryIoOperation::ReadFrame,
                &self.source_label,
                Some(self.frame_cursor),
                "DCD coordinate record exceeds the configured frame limit",
            ));
        }
        for axis in 0..3 {
            let clean_frame_eof = !self.header.has_cell && axis == 0;
            if !read_record(
                &mut self.reader,
                self.header.endian,
                &mut self.record,
                &self.limits,
                &self.source_label,
                TrajectoryIoOperation::ReadFrame,
                Some(self.frame_cursor),
                clean_frame_eof,
            )? {
                return Err(self.early_eof_error(frame_start));
            }
            if self.record.len() != coordinate_bytes {
                return Err(frame_error(
                    TrajectoryCodecErrorKind::InvalidFrame,
                    &self.source_label,
                    self.frame_cursor,
                    "DCD coordinate record has the wrong size",
                ));
            }
            for (value_index, chunk) in self.record.chunks_exact(4).enumerate() {
                let value = f64::from(
                    self.header
                        .endian
                        .f32(chunk.try_into().expect("coordinate chunk")),
                );
                if !value.is_finite() {
                    return Err(frame_error(
                        TrajectoryCodecErrorKind::InvalidFrame,
                        &self.source_label,
                        self.frame_cursor,
                        "DCD coordinate is not finite",
                    ));
                }
                let atom_index =
                    coordinate_indices.map_or(value_index, |indices| indices[value_index]);
                match axis {
                    0 => self.positions[atom_index].x = value,
                    1 => self.positions[atom_index].y = value,
                    _ => self.positions[atom_index].z = value,
                }
            }
        }
        let frame_end = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                error,
            )
        })?;
        let frame_bytes = frame_end.checked_sub(frame_start).ok_or_else(|| {
            frame_error(
                TrajectoryCodecErrorKind::InvalidFrame,
                &self.source_label,
                self.frame_cursor,
                "DCD frame end precedes its start",
            )
        })?;
        if frame_bytes > self.limits.max_frame_bytes {
            return Err(resource_error(
                TrajectoryIoOperation::ReadFrame,
                &self.source_label,
                Some(self.frame_cursor),
                "DCD frame exceeds the configured byte limit",
            ));
        }
        if capture_fixed_reference && self.frame_cursor == 0 && self.header.fixed_count > 0 {
            self.fixed_reference.clear();
            self.fixed_reference
                .try_reserve_exact(self.fixed_indices.len())
                .map_err(|_| {
                    resource_error(
                        TrajectoryIoOperation::ReadFrame,
                        &self.source_label,
                        Some(0),
                        "could not reserve DCD fixed-coordinate scratch",
                    )
                })?;
            self.fixed_reference.extend(
                self.fixed_indices
                    .iter()
                    .map(|&index| self.positions[index]),
            );
        }
        let step = self
            .header
            .start_step
            .checked_add(
                self.frame_cursor
                    .checked_mul(self.header.step_interval)
                    .ok_or_else(|| {
                        frame_error(
                            TrajectoryCodecErrorKind::NegativeOrUnrepresentableStep,
                            &self.source_label,
                            self.frame_cursor,
                            "DCD frame step multiplication overflows",
                        )
                    })?,
            )
            .ok_or_else(|| {
                frame_error(
                    TrajectoryCodecErrorKind::NegativeOrUnrepresentableStep,
                    &self.source_label,
                    self.frame_cursor,
                    "DCD frame step addition overflows",
                )
            })?;
        let time = if let DcdTimePolicy::HeaderDelta { .. } = self.options.time_policy {
            let time = (step as f64) * f64::from(self.header.delta);
            if !time.is_finite() {
                return Err(frame_error(
                    TrajectoryCodecErrorKind::InvalidFrame,
                    &self.source_label,
                    self.frame_cursor,
                    "DCD derived frame time is not finite",
                ));
            }
            Some(time)
        } else {
            None
        };
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            resource_error(
                TrajectoryIoOperation::ReadFrame,
                &self.source_label,
                Some(self.frame_cursor),
                "DCD frame cursor overflows",
            )
        })?;
        Ok(Some(DcdDecodedFrame { cell, step, time }))
    }

    fn publish(
        &self,
        positions: &[Point3],
        decoded: DcdDecodedFrame,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        let mut data = FrameBufferData::new(
            self.binding.topology_arc(),
            Quantity::new(positions, ANGSTROM),
        )
        .with_step(decoded.step);
        if let Some(cell) = decoded.cell {
            data = data.with_cell(cell);
        }
        if let (Some(time), DcdTimePolicy::HeaderDelta { unit }) =
            (decoded.time, self.options.time_policy)
        {
            data = data.with_time(Quantity::new(time, unit));
        }
        destination.replace_from_data(data).map_err(Into::into)
    }

    pub fn into_indexed(mut self) -> Result<IndexedDcdReader<R>, TrajectoryError> {
        let mut offsets = Vec::new();
        loop {
            if offsets.len() as u64 == self.header.declared_frames {
                if probe_seekable_eof(
                    &mut self.reader,
                    TrajectoryIoOperation::Index,
                    TrajectoryFormat::Dcd,
                    &self.source_label,
                )? {
                    break;
                }
                return Err(TrajectoryCodecErrorContext::new(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                )
                .with_source_label(&self.source_label)
                .with_counts(
                    self.header.declared_frames,
                    (offsets.len() as u64).saturating_add(1),
                )
                .with_detail("DCD contains frames beyond its declared NSET count")
                .into());
            }
            if let Some(limit) = projected_index_limit(offsets.len(), &self.limits) {
                if probe_seekable_eof(
                    &mut self.reader,
                    TrajectoryIoOperation::Index,
                    TrajectoryFormat::Dcd,
                    &self.source_label,
                )? {
                    return Err(TrajectoryCodecErrorContext::new(
                        TrajectoryCodecErrorKind::InconsistentMetadata,
                        TrajectoryIoOperation::Index,
                        Some(TrajectoryFormat::Dcd),
                    )
                    .with_source_label(&self.source_label)
                    .with_counts(self.header.declared_frames, offsets.len() as u64)
                    .with_detail("DCD ended before its declared NSET count")
                    .into());
                }
                return Err(resource_error(
                    TrajectoryIoOperation::Index,
                    &self.source_label,
                    Some(offsets.len() as u64),
                    format!("DCD index {limit} exceeds the configured limit"),
                ));
            }
            let offset = self.reader.stream_position().map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                    &self.source_label,
                    error,
                )
            })?;
            if self
                .parse_next(true)
                .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))?
                .is_none()
            {
                return Err(TrajectoryCodecErrorContext::new(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                )
                .with_source_label(&self.source_label)
                .with_counts(self.header.declared_frames, offsets.len() as u64)
                .with_detail("DCD ended before its declared NSET count")
                .into());
            }
            reserve_index_for_push(
                &mut offsets,
                &self.limits,
                TrajectoryFormat::Dcd,
                &self.source_label,
                self.frame_cursor.saturating_sub(1),
            )?;
            offsets.push(offset);
        }
        self.reader
            .seek(SeekFrom::Start(self.header.data_start))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                    &self.source_label,
                    error,
                )
            })?;
        self.frame_cursor = 0;
        let mut random_positions = Vec::new();
        random_positions
            .try_reserve_exact(self.header.atom_count)
            .map_err(|_| {
                resource_error(
                    TrajectoryIoOperation::Index,
                    &self.source_label,
                    None,
                    "could not reserve DCD indexed position scratch",
                )
            })?;
        random_positions.resize(self.header.atom_count, Point3::new(0.0, 0.0, 0.0));
        Ok(IndexedDcdReader {
            inner: self,
            offsets,
            random_positions,
            random_record: Vec::new(),
        })
    }
}

impl<R: Read + Seek> TrajectoryReader for DcdReader<R> {
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
        let offset = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                error,
            )
        })?;
        let decoded = self
            .parse_next(true)
            .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))?;
        let Some(decoded) = decoded else {
            return Ok(false);
        };
        self.publish(&self.positions, decoded, destination)?;
        Ok(true)
    }
}

/// Fully verified DCD indexed reader.
pub struct IndexedDcdReader<R> {
    inner: DcdReader<R>,
    offsets: Vec<u64>,
    random_positions: Vec<Point3>,
    random_record: Vec<u8>,
}

impl<R: Read + Seek> IndexedDcdReader<R> {
    pub fn topology(&self) -> &Topology {
        self.inner.topology()
    }

    pub(crate) fn header(&self) -> &DcdHeader {
        self.inner.header()
    }
}

impl<R: Read + Seek> TrajectoryReader for IndexedDcdReader<R> {
    fn topology(&self) -> &Topology {
        self.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.inner.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        self.inner.read_next(destination)
    }
}

impl<R: Read + Seek> SeekableTrajectoryReader for IndexedDcdReader<R> {
    fn frame_count(&self) -> Option<u64> {
        Some(self.offsets.len() as u64)
    }

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        if !std::ptr::eq(self.topology(), destination.topology()) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let offset = self
            .offsets
            .get(usize::try_from(index).map_err(|_| TrajectoryError::FrameIndexOutOfRange(index))?)
            .copied()
            .ok_or(TrajectoryError::FrameIndexOutOfRange(index))?;
        let saved_offset = self.inner.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Dcd),
                &self.inner.source_label,
                error,
            )
        })?;
        let saved_cursor = self.inner.frame_cursor;
        self.inner
            .reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                    &self.inner.source_label,
                    error,
                )
            })?;
        std::mem::swap(&mut self.inner.positions, &mut self.random_positions);
        std::mem::swap(&mut self.inner.record, &mut self.random_record);
        self.inner.frame_cursor = index;
        if index > 0 && self.inner.header.fixed_count > 0 {
            for (&atom_index, &position) in self
                .inner
                .fixed_indices
                .iter()
                .zip(&self.inner.fixed_reference)
            {
                self.inner.positions[atom_index] = position;
            }
        }
        let result = self
            .inner
            .parse_next(false)
            .map_err(|error| frame_offset_context(error, index, offset))
            .and_then(|decoded| decoded.ok_or(TrajectoryError::FrameIndexOutOfRange(index)));
        let restore = self
            .inner
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                    &self.inner.source_label,
                    error,
                )
            });
        self.inner.frame_cursor = saved_cursor;
        std::mem::swap(&mut self.inner.positions, &mut self.random_positions);
        std::mem::swap(&mut self.inner.record, &mut self.random_record);
        let decoded = result?;
        restore?;
        self.inner
            .publish(&self.random_positions, decoded, destination)
    }
}

/// Canonical DCD writer over one seekable stream.
pub struct DcdWriter<W> {
    writer: W,
    topology: Arc<Topology>,
    options: DcdWriteOptions,
    source_label: String,
    header_start: u64,
    frame_count: u64,
    finalized: bool,
    axis: Vec<f32>,
}

impl<W: Write + Seek> DcdWriter<W> {
    pub fn new(
        mut writer: W,
        topology: Arc<Topology>,
        options: DcdWriteOptions,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryError> {
        let source_label = source_label.into();
        if options.step_interval == 0
            || options.start_step > i32::MAX as u64
            || options.step_interval > i32::MAX as u64
            || !options.delta.is_finite()
            || options.delta < 0.0
        {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Dcd),
                &source_label,
                "DCD writer step sequence and DELTA must fit the common i32/f32 profile",
            ));
        }
        validate_time_policy(options.time_policy, &source_label)?;
        let atom_count = topology.atom_count();
        if atom_count == 0 || atom_count > i32::MAX as usize {
            return Err(codec_context(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Dcd),
                &source_label,
                "DCD writer atom count must fit a positive i32",
            ));
        }
        let header_start = writer.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Dcd),
                &source_label,
                error,
            )
        })?;
        let header = encode_header(options, 0);
        write_record(
            &mut writer,
            options.endian,
            &header,
            &source_label,
            TrajectoryIoOperation::Open,
        )?;
        let mut title = vec![0_u8; 84];
        title[..4].copy_from_slice(&options.endian.encode_i32(1));
        let text = b"Created by kekule-traj";
        title[4..4 + text.len()].copy_from_slice(text);
        write_record(
            &mut writer,
            options.endian,
            &title,
            &source_label,
            TrajectoryIoOperation::Open,
        )?;
        write_record(
            &mut writer,
            options.endian,
            &options
                .endian
                .encode_i32(i32::try_from(atom_count).expect("checked atom count")),
            &source_label,
            TrajectoryIoOperation::Open,
        )?;
        let mut axis = Vec::new();
        axis.try_reserve_exact(atom_count).map_err(|_| {
            resource_error(
                TrajectoryIoOperation::Open,
                &source_label,
                None,
                "could not reserve DCD writer scratch",
            )
        })?;
        Ok(Self {
            writer,
            topology,
            options,
            source_label,
            header_start,
            frame_count: 0,
            finalized: false,
            axis,
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

    pub fn finalize(&mut self) -> Result<(), TrajectoryError> {
        if self.finalized {
            return Ok(());
        }
        require_nonempty_writer(self.frame_count, TrajectoryFormat::Dcd, &self.source_label)?;
        if self.frame_count > i32::MAX as u64 {
            return Err(codec_context(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::Finish,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                "DCD frame count exceeds i32 header capacity",
            ));
        }
        let end = self.writer.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Finish,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                error,
            )
        })?;
        let count_offset = self.header_start.checked_add(8).ok_or_else(|| {
            codec_context(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::Finish,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                "DCD header count offset overflows",
            )
        })?;
        self.writer
            .seek(SeekFrom::Start(count_offset))
            .and_then(|_| {
                self.writer.write_all(
                    &self
                        .options
                        .endian
                        .encode_i32(i32::try_from(self.frame_count).expect("checked frame count")),
                )
            })
            .and_then(|_| self.writer.seek(SeekFrom::Start(end)).map(|_| ()))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Finish,
                    Some(TrajectoryFormat::Dcd),
                    &self.source_label,
                    error,
                )
            })?;
        self.finalized = true;
        Ok(())
    }

    /// Finalizes the frame count, flushes, and returns the nonempty DCD stream.
    pub fn finish(mut self) -> Result<W, TrajectoryError> {
        self.finalize()?;
        self.writer.flush().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Finish,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                error,
            )
        })?;
        Ok(self.writer)
    }
}

impl<W: Write + Seek> TrajectoryWriter for DcdWriter<W> {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if self.finalized {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                "cannot write a DCD frame after finalization",
            ));
        }
        if !std::ptr::eq(frame.topology(), self.topology()) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        if frame.velocities().is_some() {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                "velocities",
            ));
        }
        if frame.forces().is_some() {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                "forces",
            ));
        }
        if frame.observation().is_some() {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                "observation",
            ));
        }
        if !frame.props().is_empty() {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                "properties",
            ));
        }
        if self.frame_count >= i32::MAX as u64 {
            return Err(writer_overflow(
                &self.source_label,
                "DCD frame count exceeds i32 header capacity",
            ));
        }
        let expected_step = self
            .options
            .start_step
            .checked_add(
                self.frame_count
                    .checked_mul(self.options.step_interval)
                    .ok_or_else(|| writer_overflow(&self.source_label, "DCD step overflow"))?,
            )
            .ok_or_else(|| writer_overflow(&self.source_label, "DCD step overflow"))?;
        if frame.step() != Some(expected_step) {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                format!("DCD frame step must be exactly {expected_step}"),
            ));
        }
        match (self.options.time_policy, frame.time()) {
            (DcdTimePolicy::Absent, None) => {}
            (DcdTimePolicy::Absent, Some(_)) => {
                return Err(TrajectoryError::UnsupportedField("time"))
            }
            (DcdTimePolicy::HeaderDelta { unit }, Some(time)) => {
                let actual = time.value_in(unit).map_err(|error| {
                    codec_context(
                        TrajectoryCodecErrorKind::InconsistentMetadata,
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Dcd),
                        &self.source_label,
                        format!("DCD time unit is incompatible: {error}"),
                    )
                })?;
                let expected = (expected_step as f64) * f64::from(self.options.delta);
                if !actual.is_finite()
                    || (actual - expected).abs() > 1.0e-9 * expected.abs().max(1.0)
                {
                    return Err(codec_context(
                        TrajectoryCodecErrorKind::InconsistentMetadata,
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Dcd),
                        &self.source_label,
                        "DCD frame time does not match the explicit DELTA policy",
                    ));
                }
            }
            (DcdTimePolicy::HeaderDelta { .. }, None) => {
                return Err(codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Dcd),
                    &self.source_label,
                    "DCD explicit time policy requires frame time",
                ))
            }
        }
        let cell = frame.configuration().cell().copied();
        if self.options.write_cell != cell.is_some() {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                "DCD cell presence must match the writer header policy",
            ));
        }
        let factor = MODEL_LENGTH_UNIT
            .conversion_factor_to(ANGSTROM)
            .map_err(|error| {
                codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Dcd),
                    &self.source_label,
                    format!("DCD position unit is incompatible: {error}"),
                )
            })?;
        let positions = frame.configuration().positions().values();
        for point in *positions.value() {
            for value in [point.x * factor, point.y * factor, point.z * factor] {
                if !value.is_finite() || !(value as f32).is_finite() {
                    return Err(codec_context(
                        TrajectoryCodecErrorKind::InvalidFrame,
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Dcd),
                        &self.source_label,
                        "DCD coordinate cannot be represented as finite f32",
                    ));
                }
            }
        }
        if let Some(cell) = cell {
            let values = encode_cell(cell, &self.source_label)?;
            let mut payload = [0_u8; CELL_BYTES];
            for (chunk, value) in payload.chunks_exact_mut(8).zip(values) {
                chunk.copy_from_slice(&self.options.endian.encode_f64(value));
            }
            write_record(
                &mut self.writer,
                self.options.endian,
                &payload,
                &self.source_label,
                TrajectoryIoOperation::WriteFrame,
            )?;
        }
        for axis in 0..3 {
            self.axis.clear();
            for point in *positions.value() {
                let value = match axis {
                    0 => point.x,
                    1 => point.y,
                    _ => point.z,
                } * factor;
                let narrowed = value as f32;
                self.axis.push(narrowed);
            }
            let byte_len =
                self.axis.len().checked_mul(4).ok_or_else(|| {
                    writer_overflow(&self.source_label, "DCD record size overflows")
                })?;
            if byte_len > u32::MAX as usize {
                return Err(writer_overflow(
                    &self.source_label,
                    "DCD coordinate record exceeds u32",
                ));
            }
            self.writer
                .write_all(&self.options.endian.encode_u32(byte_len as u32))
                .and_then(|_| {
                    for value in &self.axis {
                        self.writer
                            .write_all(&self.options.endian.encode_f32(*value))?;
                    }
                    self.writer
                        .write_all(&self.options.endian.encode_u32(byte_len as u32))
                })
                .map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Dcd),
                        &self.source_label,
                        error,
                    )
                })?;
        }
        self.frame_count = self
            .frame_count
            .checked_add(1)
            .ok_or_else(|| writer_overflow(&self.source_label, "DCD frame count overflows"))?;
        Ok(())
    }
}

fn encode_header(options: DcdWriteOptions, frames: i32) -> [u8; HEADER_BYTES] {
    let mut header = [0_u8; HEADER_BYTES];
    header[..4].copy_from_slice(b"CORD");
    let controls = [
        frames,
        options.start_step as i32,
        options.step_interval as i32,
        0,
        0,
        0,
        0,
        0,
        0,
        i32::from_ne_bytes(options.delta.to_bits().to_ne_bytes()),
        i32::from(options.write_cell),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        24,
    ];
    for (chunk, control) in header[4..].chunks_exact_mut(4).zip(controls) {
        chunk.copy_from_slice(&options.endian.encode_i32(control));
    }
    header[40..44].copy_from_slice(&options.endian.encode_f32(options.delta));
    header
}

fn parse_controls(
    bytes: &[u8],
    endian: DcdEndian,
    source_label: &str,
) -> Result<[i32; 20], TrajectoryError> {
    if bytes.len() != 80 {
        return Err(header_error(
            source_label,
            "DCD control block is not 80 bytes",
        ));
    }
    let mut controls = [0_i32; 20];
    for (target, chunk) in controls.iter_mut().zip(bytes.chunks_exact(4)) {
        *target = endian.i32(chunk.try_into().expect("control chunk"));
    }
    Ok(controls)
}

fn validate_title_record(
    record: &[u8],
    endian: DcdEndian,
    source_label: &str,
) -> Result<(), TrajectoryError> {
    if record.len() < 4 {
        return Err(header_error(source_label, "DCD title record is too short"));
    }
    let count = nonnegative_usize(
        endian.i32(record[..4].try_into().expect("title count")),
        "NTITLE",
        source_label,
    )?;
    let expected = count
        .checked_mul(80)
        .and_then(|bytes| bytes.checked_add(4))
        .ok_or_else(|| header_error(source_label, "DCD title record size overflows"))?;
    if record.len() != expected {
        return Err(header_error(
            source_label,
            "DCD title count does not match record size",
        ));
    }
    Ok(())
}

fn validate_time_policy(policy: DcdTimePolicy, source_label: &str) -> Result<(), TrajectoryError> {
    if let DcdTimePolicy::HeaderDelta { unit } = policy {
        unit.conversion_factor_to(MODEL_TIME_UNIT)
            .map_err(|error| {
                codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Open,
                    Some(TrajectoryFormat::Dcd),
                    source_label,
                    format!("DCD time-policy unit is incompatible: {error}"),
                )
            })?;
    }
    Ok(())
}

fn decode_cell(
    record: &[u8],
    endian: DcdEndian,
    source_label: &str,
    frame: u64,
) -> Result<PeriodicCell, TrajectoryError> {
    let mut values = [0_f64; 6];
    for (target, chunk) in values.iter_mut().zip(record.chunks_exact(8)) {
        *target = endian.f64(chunk.try_into().expect("cell chunk"));
    }
    if values.iter().any(|value| !value.is_finite()) {
        return Err(frame_error(
            TrajectoryCodecErrorKind::InvalidFrame,
            source_label,
            frame,
            "DCD cell contains a non-finite value",
        ));
    }
    let [a, gamma_raw, b, beta_raw, alpha_raw, c] = values;
    if a <= 0.0 || b <= 0.0 || c <= 0.0 {
        return Err(frame_error(
            TrajectoryCodecErrorKind::InvalidFrame,
            source_label,
            frame,
            "DCD cell lengths must be positive",
        ));
    }
    let cosine_variant = [alpha_raw, beta_raw, gamma_raw]
        .into_iter()
        .all(|value| (-1.0..=1.0).contains(&value));
    let angle = |value: f64| {
        if cosine_variant {
            value.acos()
        } else {
            value.to_radians()
        }
    };
    let alpha = angle(alpha_raw);
    let beta = angle(beta_raw);
    let gamma = angle(gamma_raw);
    let sin_gamma = gamma.sin();
    if !alpha.is_finite()
        || !beta.is_finite()
        || !gamma.is_finite()
        || sin_gamma.abs() <= f64::EPSILON
    {
        return Err(frame_error(
            TrajectoryCodecErrorKind::InvalidFrame,
            source_label,
            frame,
            "DCD cell angles are degenerate",
        ));
    }
    let ax = Vector3::new(a, 0.0, 0.0);
    let bx = Vector3::new(b * gamma.cos(), b * sin_gamma, 0.0);
    let cx = c * beta.cos();
    let cy = c * (alpha.cos() - beta.cos() * gamma.cos()) / sin_gamma;
    let cz_squared = c * c - cx * cx - cy * cy;
    if !cz_squared.is_finite() || cz_squared <= 0.0 {
        return Err(frame_error(
            TrajectoryCodecErrorKind::InvalidFrame,
            source_label,
            frame,
            "DCD cell vectors are degenerate",
        ));
    }
    PeriodicCell::new(
        Quantity::new([ax, bx, Vector3::new(cx, cy, cz_squared.sqrt())], ANGSTROM),
        [true; 3],
    )
    .map_err(|error| {
        frame_error(
            TrajectoryCodecErrorKind::InvalidFrame,
            source_label,
            frame,
            format!("invalid DCD cell: {error}"),
        )
    })
}

fn encode_cell(cell: PeriodicCell, source_label: &str) -> Result<[f64; 6], TrajectoryError> {
    if cell.periodic_axes() != [true; 3] {
        return Err(codec_context(
            TrajectoryCodecErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Dcd),
            source_label,
            "DCD writer requires all three periodic axes",
        ));
    }
    let [a, b, c] = cell.vectors().value_in(ANGSTROM).map_err(|error| {
        codec_context(
            TrajectoryCodecErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Dcd),
            source_label,
            format!("DCD cell unit is incompatible: {error}"),
        )
    })?;
    let length =
        |vector: Vector3| (vector.x * vector.x + vector.y * vector.y + vector.z * vector.z).sqrt();
    let dot =
        |left: Vector3, right: Vector3| left.x * right.x + left.y * right.y + left.z * right.z;
    let la = length(a);
    let lb = length(b);
    let lc = length(c);
    let degrees = |left: Vector3, right: Vector3, ll: f64, lr: f64| {
        (dot(left, right) / (ll * lr))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    };
    let alpha = degrees(b, c, lb, lc);
    let beta = degrees(a, c, la, lc);
    let gamma = degrees(a, b, la, lb);
    let values = [la, gamma, lb, beta, alpha, lc];
    if values.iter().any(|value| !value.is_finite()) {
        return Err(codec_context(
            TrajectoryCodecErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Dcd),
            source_label,
            "DCD cell cannot be represented as finite lengths and angles",
        ));
    }
    Ok(values)
}

#[allow(clippy::too_many_arguments)]
fn read_record<R: Read>(
    reader: &mut R,
    endian: DcdEndian,
    scratch: &mut Vec<u8>,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    operation: TrajectoryIoOperation,
    frame: Option<u64>,
    clean_eof: bool,
) -> Result<bool, TrajectoryError> {
    let mut marker = [0_u8; 4];
    match reader.read(&mut marker[..1]) {
        Ok(0) if clean_eof => return Ok(false),
        Ok(0) => {
            return Err(record_error(
                source_label,
                operation,
                frame,
                "missing DCD record marker",
            ))
        }
        Ok(_) => {}
        Err(error) => {
            return Err(io_context(
                operation,
                Some(TrajectoryFormat::Dcd),
                source_label,
                error,
            ))
        }
    }
    read_exact_required(
        reader,
        &mut marker[1..],
        operation,
        source_label,
        frame,
        "partial DCD record marker",
    )?;
    let size = usize::try_from(endian.u32(marker)).map_err(|_| {
        resource_error(
            operation,
            source_label,
            frame,
            "DCD record size does not fit usize",
        )
    })?;
    if size as u64 > limits.max_record_bytes || size > limits.max_scratch_bytes {
        return Err(resource_error(
            operation,
            source_label,
            frame,
            "DCD record exceeds configured record or scratch limit",
        ));
    }
    if scratch.capacity() < size {
        scratch
            .try_reserve_exact(size.saturating_sub(scratch.len()))
            .map_err(|_| {
                resource_error(
                    operation,
                    source_label,
                    frame,
                    "could not reserve DCD record scratch",
                )
            })?;
    }
    scratch.resize(size, 0);
    read_exact_required(
        reader,
        scratch,
        operation,
        source_label,
        frame,
        "truncated DCD record payload",
    )?;
    let mut trailing = [0_u8; 4];
    read_exact_required(
        reader,
        &mut trailing,
        operation,
        source_label,
        frame,
        "missing DCD trailing record marker",
    )?;
    if endian.u32(trailing) != size as u32 {
        let mut context = TrajectoryCodecErrorContext::new(
            TrajectoryCodecErrorKind::RecordMarkerMismatch,
            operation,
            Some(TrajectoryFormat::Dcd),
        )
        .with_source_label(source_label)
        .with_counts(size as u64, u64::from(endian.u32(trailing)))
        .with_detail("DCD leading and trailing record markers differ");
        if let Some(frame) = frame {
            context = context.with_frame(frame);
        }
        return Err(context.into());
    }
    Ok(true)
}

fn read_exact_required<R: Read>(
    reader: &mut R,
    bytes: &mut [u8],
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: Option<u64>,
    detail: &str,
) -> Result<(), TrajectoryError> {
    reader.read_exact(bytes).map_err(|error| {
        if error.kind() == std::io::ErrorKind::UnexpectedEof {
            record_error(source_label, operation, frame, detail)
        } else {
            io_context(operation, Some(TrajectoryFormat::Dcd), source_label, error)
        }
    })
}

fn write_record<W: Write>(
    writer: &mut W,
    endian: DcdEndian,
    payload: &[u8],
    source_label: &str,
    operation: TrajectoryIoOperation,
) -> Result<(), TrajectoryError> {
    let size = u32::try_from(payload.len())
        .map_err(|_| writer_overflow(source_label, "DCD record length exceeds u32 capacity"))?;
    writer
        .write_all(&endian.encode_u32(size))
        .and_then(|_| writer.write_all(payload))
        .and_then(|_| writer.write_all(&endian.encode_u32(size)))
        .map_err(|error| io_context(operation, Some(TrajectoryFormat::Dcd), source_label, error))
}

fn validate_atom_count(
    count: usize,
    limits: &TrajectoryIoLimits,
    source_label: &str,
) -> Result<(), TrajectoryError> {
    if count == 0 || count > limits.max_atoms || count > i32::MAX as usize {
        return Err(resource_error(
            TrajectoryIoOperation::ReadHeader,
            source_label,
            None,
            "DCD atom count is zero or exceeds configured/profile limits",
        ));
    }
    Ok(())
}

fn checked_coordinate_bytes(count: usize, source_label: &str) -> Result<usize, TrajectoryError> {
    count.checked_mul(4).ok_or_else(|| {
        resource_error(
            TrajectoryIoOperation::ReadFrame,
            source_label,
            None,
            "DCD coordinate record size overflows",
        )
    })
}

fn nonnegative_u64(value: i32, field: &str, source_label: &str) -> Result<u64, TrajectoryError> {
    u64::try_from(value)
        .map_err(|_| header_error(source_label, format!("DCD {field} must be nonnegative")))
}

fn positive_u64(value: i32, field: &str, source_label: &str) -> Result<u64, TrajectoryError> {
    if value <= 0 {
        return Err(header_error(
            source_label,
            format!("DCD {field} must be positive"),
        ));
    }
    Ok(value as u64)
}

fn nonnegative_usize(
    value: i32,
    field: &str,
    source_label: &str,
) -> Result<usize, TrajectoryError> {
    usize::try_from(value)
        .map_err(|_| header_error(source_label, format!("DCD {field} must be nonnegative")))
}

fn header_error(source_label: &str, detail: impl Into<String>) -> TrajectoryError {
    codec_context(
        TrajectoryCodecErrorKind::InvalidHeader,
        TrajectoryIoOperation::ReadHeader,
        Some(TrajectoryFormat::Dcd),
        source_label,
        detail,
    )
}

fn record_error(
    source_label: &str,
    operation: TrajectoryIoOperation,
    frame: Option<u64>,
    detail: impl Into<String>,
) -> TrajectoryError {
    let mut context = TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::TruncatedRecord,
        operation,
        Some(TrajectoryFormat::Dcd),
    )
    .with_source_label(source_label)
    .with_detail(detail);
    if let Some(frame) = frame {
        context = context.with_frame(frame);
    }
    context.into()
}

fn frame_error(
    kind: TrajectoryCodecErrorKind,
    source_label: &str,
    frame: u64,
    detail: impl Into<String>,
) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(
        kind,
        TrajectoryIoOperation::ReadFrame,
        Some(TrajectoryFormat::Dcd),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_detail(detail)
    .into()
}

fn resource_error(
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: Option<u64>,
    detail: impl Into<String>,
) -> TrajectoryError {
    let mut context = TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::ResourceLimitExceeded,
        operation,
        Some(TrajectoryFormat::Dcd),
    )
    .with_source_label(source_label)
    .with_detail(detail);
    if let Some(frame) = frame {
        context = context.with_frame(frame);
    }
    context.into()
}

fn writer_overflow(source_label: &str, detail: impl Into<String>) -> TrajectoryError {
    codec_context(
        TrajectoryCodecErrorKind::ResourceLimitExceeded,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Dcd),
        source_label,
        detail,
    )
}

fn writer_field_error(source_label: &str, frame: u64, field: &str) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::UnsupportedField,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Dcd),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_detail(format!("DCD cannot preserve {field}"))
    .into()
}
