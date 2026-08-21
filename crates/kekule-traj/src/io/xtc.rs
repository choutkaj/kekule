//! Bounded GROMACS XTC trajectories with a checked reader and `molly` writer.

use std::io::{Read, Seek, SeekFrom, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::Arc;

use crate::{
    FrameBuffer, FrameBufferData, SeekableTrajectoryReader, TrajectoryCodecErrorContext,
    TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat, TrajectoryFrameView,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::topology::Topology;
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

use super::{
    codec_context, frame_offset_context, io_context, probe_seekable_eof, projected_index_limit,
    require_nonempty_writer, reserve_index_for_push, TrajectoryIoLimits, TrajectoryTopologyBinding,
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

/// Supported XTC magic and compressed-byte-count profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum XtcMagic {
    Xtc1995,
    Xtc2023,
}

impl XtcMagic {
    const fn nbytes_width(self) -> usize {
        match self {
            Self::Xtc1995 => 4,
            Self::Xtc2023 => 8,
        }
    }

    fn to_molly(self) -> molly::Magic {
        match self {
            Self::Xtc1995 => molly::Magic::Xtc1995,
            Self::Xtc2023 => molly::Magic::Xtc2023,
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
    pub const fn with_cell_policy(mut self, cell_policy: XtcCellPolicy) -> Self {
        self.cell_policy = cell_policy;
        self
    }

    pub const fn cell_policy(self) -> XtcCellPolicy {
        self.cell_policy
    }
}

/// XTC writer policy with explicit lossy coordinate precision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XtcWriteOptions {
    magic: XtcMagic,
    precision: f32,
}

impl Default for XtcWriteOptions {
    fn default() -> Self {
        Self {
            magic: XtcMagic::Xtc1995,
            precision: 1000.0,
        }
    }
}

impl XtcWriteOptions {
    pub const fn with_magic(mut self, magic: XtcMagic) -> Self {
        self.magic = magic;
        self
    }

    /// Sets the coordinate scale in inverse nanometers.
    ///
    /// Compressed coordinates have nominal resolution `1 / precision` nm.
    pub fn with_precision(mut self, precision: f32) -> Result<Self, TrajectoryError> {
        validate_precision(precision, "XTC writer options")?;
        self.precision = precision;
        Ok(self)
    }

    pub const fn magic(self) -> XtcMagic {
        self.magic
    }

    pub const fn precision(self) -> f32 {
        self.precision
    }
}

#[derive(Debug, Clone)]
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

    pub(crate) const fn compressed_bytes(&self) -> usize {
        self.compressed_bytes
    }
}

struct CheckedXtcReaderAdapter<R> {
    reader: R,
    scratch: Vec<u8>,
    decoded_positions: Vec<f32>,
    decoded_start: Option<u64>,
}

impl<R: Read + Seek> CheckedXtcReaderAdapter<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            scratch: Vec::new(),
            decoded_positions: Vec::new(),
            decoded_start: None,
        }
    }

    fn preflight(
        &mut self,
        binding: &TrajectoryTopologyBinding,
        options: XtcReadOptions,
        limits: &TrajectoryIoLimits,
        source_label: &str,
        frame_index: u64,
        clean_eof: bool,
    ) -> Result<Option<XtcFrameInfo>, TrajectoryError> {
        self.decoded_start = None;
        let start = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )
        })?;
        let mut header = [0_u8; XTC_HEADER_BYTES];
        match self.reader.read(&mut header[..1]) {
            Ok(0) if clean_eof => return Ok(None),
            Ok(0) => return Err(header_error(source_label, "missing XTC frame header")),
            Ok(_) => {}
            Err(error) => {
                return Err(io_context(
                    TrajectoryIoOperation::ReadHeader,
                    Some(TrajectoryFormat::Xtc),
                    source_label,
                    error,
                ))
            }
        }
        read_exact(
            &mut self.reader,
            &mut header[1..],
            source_label,
            Some(frame_index),
            "truncated XTC frame header",
        )?;
        let magic = match i32::from_be_bytes(header[..4].try_into().expect("magic")) {
            1995 => XtcMagic::Xtc1995,
            2023 => XtcMagic::Xtc2023,
            _ => return Err(header_error(source_label, "XTC magic is not 1995 or 2023")),
        };
        let atom_count_raw = i32::from_be_bytes(header[4..8].try_into().expect("natoms"));
        let atom_count = nonnegative_xdr_count(
            atom_count_raw,
            "atom count",
            source_label,
            frame_index,
            start + 4,
        )?;
        validate_atom_count(atom_count, limits, source_label)?;
        let expected_atoms = binding.topology().atom_count();
        if atom_count != expected_atoms {
            return Err(TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
            )
            .with_source_label(source_label)
            .with_frame(frame_index)
            .with_counts(expected_atoms as u64, atom_count as u64)
            .into());
        }
        let repeated_raw = i32::from_be_bytes(header[52..56].try_into().expect("repeat"));
        let repeated = nonnegative_xdr_count(
            repeated_raw,
            "repeated atom count",
            source_label,
            frame_index,
            start + 52,
        )?;
        if repeated != atom_count {
            return Err(TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
            )
            .with_source_label(source_label)
            .with_frame(frame_index)
            .with_counts(atom_count as u64, repeated as u64)
            .with_detail("XTC repeated atom count differs from its header")
            .into());
        }
        let step_raw = i32::from_be_bytes(header[8..12].try_into().expect("step"));
        let step = u64::try_from(step_raw).map_err(|_| {
            TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::NegativeOrUnrepresentableStep,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
            )
            .with_source_label(source_label)
            .with_frame(frame_index)
            .with_byte_offset(start + 8)
            .with_detail("XTC step is a negative signed XDR integer")
        })?;
        let time = f64::from(f32::from_be_bytes(header[12..16].try_into().expect("time")));
        if !time.is_finite() {
            return Err(frame_error(
                source_label,
                frame_index,
                "XTC frame time is not finite",
            ));
        }
        let mut box_values = [0_f64; 9];
        for (value, chunk) in box_values.iter_mut().zip(header[16..52].chunks_exact(4)) {
            *value = f64::from(f32::from_be_bytes(chunk.try_into().expect("box")));
        }
        let cell = decode_cell(box_values, options.cell_policy, source_label, frame_index)?;
        let base_scratch = atom_count
            .checked_mul(3 * std::mem::size_of::<f32>() + std::mem::size_of::<Point3>())
            .ok_or_else(|| resource_error(source_label, None, "XTC scratch size overflows"))?;
        let coordinate_scalars = atom_count
            .checked_mul(3)
            .ok_or_else(|| resource_error(source_label, None, "XTC coordinate count overflows"))?;
        self.decoded_positions.clear();
        self.decoded_positions
            .try_reserve_exact(coordinate_scalars)
            .map_err(|_| {
                resource_error(
                    source_label,
                    Some(frame_index),
                    "could not reserve XTC decoded coordinate scratch",
                )
            })?;
        let (precision, compressed_bytes, end) = if atom_count <= 9 {
            let coordinate_bytes = atom_count
                .checked_mul(12)
                .ok_or_else(|| resource_error(source_label, None, "XTC small frame overflows"))?;
            self.scratch.clear();
            self.scratch
                .try_reserve_exact(coordinate_bytes)
                .map_err(|_| {
                    resource_error(
                        source_label,
                        Some(frame_index),
                        "could not reserve XTC small-frame scratch",
                    )
                })?;
            self.scratch.resize(coordinate_bytes, 0);
            read_exact(
                &mut self.reader,
                &mut self.scratch,
                source_label,
                Some(frame_index),
                "truncated XTC small-frame coordinates",
            )?;
            for chunk in self.scratch.chunks_exact(4) {
                let value = f32::from_be_bytes(chunk.try_into().expect("small coordinate"));
                if !value.is_finite() {
                    return Err(corrupt_error(
                        source_label,
                        frame_index,
                        "XTC small-frame coordinate is not finite",
                    ));
                }
                self.decoded_positions.push(value);
            }
            let end = start
                .checked_add(XTC_HEADER_BYTES as u64)
                .and_then(|offset| offset.checked_add(coordinate_bytes as u64))
                .ok_or_else(|| resource_error(source_label, None, "XTC frame end overflows"))?;
            (None, coordinate_bytes, end)
        } else {
            let mut prelude = [0_u8; XTC_COMPRESSED_PRELUDE_BYTES];
            read_exact(
                &mut self.reader,
                &mut prelude,
                source_label,
                Some(frame_index),
                "truncated XTC compressed prelude",
            )?;
            let precision = f32::from_be_bytes(prelude[..4].try_into().expect("precision"));
            validate_precision(precision, source_label)?;
            let mut minimum = [0_i32; 3];
            let mut maximum = [0_i32; 3];
            for axis in 0..3 {
                let offset = 4 + axis * 4;
                let min =
                    i32::from_be_bytes(prelude[offset..offset + 4].try_into().expect("minimum"));
                let max_offset = 16 + axis * 4;
                let max = i32::from_be_bytes(
                    prelude[max_offset..max_offset + 4]
                        .try_into()
                        .expect("maximum"),
                );
                minimum[axis] = min;
                maximum[axis] = max;
                if min > max {
                    return Err(corrupt_error(
                        source_label,
                        frame_index,
                        "XTC compressed coordinate bounds are reversed",
                    ));
                }
            }
            let small_index = usize::try_from(u32::from_be_bytes(
                prelude[28..32].try_into().expect("small index"),
            ))
            .map_err(|_| {
                resource_error(
                    source_label,
                    Some(frame_index),
                    "XTC small-index does not fit",
                )
            })?;
            if !(XTC_FIRST_SMALL_INDEX..XTC_MAGIC_INTS.len()).contains(&small_index) {
                return Err(corrupt_error(
                    source_label,
                    frame_index,
                    "XTC compressed small-index is outside the audited table",
                ));
            }
            let compressed_bytes = read_nbytes(&mut self.reader, magic, source_label, frame_index)?;
            if compressed_bytes == 0
                || compressed_bytes as u64 > limits.max_record_bytes
                || compressed_bytes > limits.max_scratch_bytes
            {
                return Err(resource_error(
                    source_label,
                    Some(frame_index),
                    "XTC compressed payload is zero or exceeds configured limits",
                ));
            }
            if base_scratch
                .checked_add(compressed_bytes)
                .is_none_or(|bytes| bytes > limits.max_scratch_bytes)
            {
                return Err(resource_error(
                    source_label,
                    Some(frame_index),
                    "XTC aggregate decode scratch exceeds the configured limit",
                ));
            }
            self.scratch.clear();
            self.scratch
                .try_reserve_exact(compressed_bytes)
                .map_err(|_| {
                    resource_error(
                        source_label,
                        Some(frame_index),
                        "could not reserve XTC compressed validation scratch",
                    )
                })?;
            self.scratch.resize(compressed_bytes, 0);
            read_exact(
                &mut self.reader,
                &mut self.scratch,
                source_label,
                Some(frame_index),
                "truncated XTC compressed payload",
            )?;
            decode_compressed_payload(
                &self.scratch,
                CompressedLayout {
                    atom_count,
                    minimum,
                    maximum,
                    initial_small_index: small_index,
                    precision,
                },
                &mut self.decoded_positions,
                source_label,
                frame_index,
            )?;
            let padded = compressed_bytes
                .checked_add((4 - compressed_bytes % 4) % 4)
                .ok_or_else(|| resource_error(source_label, None, "XTC padding overflows"))?;
            let end = start
                .checked_add(XTC_HEADER_BYTES as u64)
                .and_then(|offset| offset.checked_add(XTC_COMPRESSED_PRELUDE_BYTES as u64))
                .and_then(|offset| offset.checked_add(magic.nbytes_width() as u64))
                .and_then(|offset| offset.checked_add(padded as u64))
                .ok_or_else(|| resource_error(source_label, None, "XTC frame end overflows"))?;
            (Some(precision), compressed_bytes, end)
        };
        let frame_bytes = end
            .checked_sub(start)
            .ok_or_else(|| resource_error(source_label, None, "XTC frame size underflows"))?;
        if frame_bytes > limits.max_frame_bytes {
            return Err(resource_error(
                source_label,
                Some(frame_index),
                "XTC frame exceeds the configured byte limit",
            ));
        }
        if base_scratch
            .checked_add(compressed_bytes)
            .is_none_or(|bytes| bytes > limits.max_scratch_bytes)
        {
            return Err(resource_error(
                source_label,
                Some(frame_index),
                "XTC aggregate decode scratch exceeds the configured limit",
            ));
        }
        let file_end = self.reader.seek(SeekFrom::End(0)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )
        })?;
        if end > file_end {
            self.reader.seek(SeekFrom::Start(start)).map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadHeader,
                    Some(TrajectoryFormat::Xtc),
                    source_label,
                    error,
                )
            })?;
            return Err(truncated_error(
                source_label,
                Some(frame_index),
                "XTC frame payload ends beyond the stream",
            ));
        }
        if let Some(precision) = precision {
            let padding = (4 - compressed_bytes % 4) % 4;
            if padding > 0 {
                self.reader
                    .seek(SeekFrom::Start(end - padding as u64))
                    .map_err(|error| {
                        io_context(
                            TrajectoryIoOperation::ReadHeader,
                            Some(TrajectoryFormat::Xtc),
                            source_label,
                            error,
                        )
                    })?;
                let mut pad = [0_u8; 3];
                read_exact(
                    &mut self.reader,
                    &mut pad[..padding],
                    source_label,
                    Some(frame_index),
                    "truncated XTC padding",
                )?;
                if pad[..padding].iter().any(|byte| *byte != 0) {
                    return Err(corrupt_error(
                        source_label,
                        frame_index,
                        "XTC XDR padding is not zero",
                    ));
                }
            }
            let _ = precision;
        }
        self.reader.seek(SeekFrom::Start(start)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )
        })?;
        self.decoded_start = Some(start);
        Ok(Some(XtcFrameInfo {
            start,
            end,
            magic,
            atom_count,
            step,
            time,
            cell,
            precision,
            compressed_bytes,
        }))
    }

    fn decode(
        &mut self,
        info: &XtcFrameInfo,
        source_label: &str,
        frame_index: u64,
    ) -> Result<(), TrajectoryError> {
        if self.decoded_start != Some(info.start) {
            return Err(corrupt_error(
                source_label,
                frame_index,
                "XTC preflight cache does not match the requested frame",
            ));
        }
        self.reader
            .seek(SeekFrom::Start(info.end))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadFrame,
                    Some(TrajectoryFormat::Xtc),
                    source_label,
                    error,
                )
            })?;
        let position = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )
        })?;
        if position != info.end || self.decoded_positions.len() != info.atom_count.saturating_mul(3)
        {
            return Err(corrupt_error(
                source_label,
                frame_index,
                "XTC adapter output disagrees with bounded preflight metadata",
            ));
        }
        if self
            .decoded_positions
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(corrupt_error(
                source_label,
                frame_index,
                "XTC decoder produced a non-finite coordinate",
            ));
        }
        self.decoded_start = None;
        Ok(())
    }
}

/// Sequential XTC reader using a private checked decoder for untrusted input.
pub struct XtcReader<R> {
    adapter: CheckedXtcReaderAdapter<R>,
    binding: TrajectoryTopologyBinding,
    options: XtcReadOptions,
    limits: TrajectoryIoLimits,
    source_label: String,
    stream_start: u64,
    first_info: XtcFrameInfo,
    pending_info: Option<XtcFrameInfo>,
    positions: Vec<Point3>,
    frame_cursor: u64,
}

impl<R: Read + Seek> XtcReader<R> {
    pub fn new(
        reader: R,
        binding: TrajectoryTopologyBinding,
        options: XtcReadOptions,
        limits: TrajectoryIoLimits,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryError> {
        let source_label = source_label.into();
        validate_atom_count(binding.topology().atom_count(), &limits, &source_label)?;
        let mut adapter = CheckedXtcReaderAdapter::new(reader);
        let stream_start = adapter.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xtc),
                &source_label,
                error,
            )
        })?;
        let first_info = adapter
            .preflight(&binding, options, &limits, &source_label, 0, false)
            .map_err(|error| frame_offset_context(error, 0, stream_start))?
            .ok_or_else(|| header_error(&source_label, "XTC stream is empty"))?;
        let atom_count = binding.topology().atom_count();
        let mut positions = Vec::new();
        positions.try_reserve_exact(atom_count).map_err(|_| {
            resource_error(
                &source_label,
                None,
                "could not reserve XTC position scratch",
            )
        })?;
        positions.resize(atom_count, Point3::new(0.0, 0.0, 0.0));
        Ok(Self {
            adapter,
            binding,
            options,
            limits,
            source_label,
            stream_start,
            first_info: first_info.clone(),
            pending_info: Some(first_info),
            positions,
            frame_cursor: 0,
        })
    }

    pub fn topology(&self) -> &Topology {
        self.binding.topology()
    }

    pub(crate) fn first_info(&self) -> &XtcFrameInfo {
        &self.first_info
    }

    fn next_info(&mut self) -> Result<Option<XtcFrameInfo>, TrajectoryError> {
        if let Some(info) = self.pending_info.take() {
            return Ok(Some(info));
        }
        let offset = self.adapter.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                error,
            )
        })?;
        self.adapter
            .preflight(
                &self.binding,
                self.options,
                &self.limits,
                &self.source_label,
                self.frame_cursor,
                true,
            )
            .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))
    }

    fn parse_next(&mut self) -> Result<Option<XtcFrameInfo>, TrajectoryError> {
        if self.frame_cursor >= self.limits.max_frames {
            let clean_eof = if self.pending_info.is_some() {
                false
            } else {
                probe_seekable_eof(
                    &mut self.adapter.reader,
                    TrajectoryIoOperation::ReadFrame,
                    TrajectoryFormat::Xtc,
                    &self.source_label,
                )?
            };
            if clean_eof {
                return Ok(None);
            }
            return Err(resource_error(
                &self.source_label,
                Some(self.frame_cursor),
                "XTC frame count exceeds the configured limit",
            ));
        }
        let Some(info) = self.next_info()? else {
            return Ok(None);
        };
        if info.magic != self.first_info.magic || info.precision != self.first_info.precision {
            return Err(TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Xtc),
            )
            .with_source_label(&self.source_label)
            .with_frame(self.frame_cursor)
            .with_byte_offset(info.start)
            .with_detail("XTC magic and coordinate precision must remain constant across frames")
            .into());
        }
        if self.adapter.decoded_start != Some(info.start) {
            self.adapter
                .reader
                .seek(SeekFrom::Start(info.start))
                .map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::ReadFrame,
                        Some(TrajectoryFormat::Xtc),
                        &self.source_label,
                        error,
                    )
                })?;
            let refreshed = self
                .adapter
                .preflight(
                    &self.binding,
                    self.options,
                    &self.limits,
                    &self.source_label,
                    self.frame_cursor,
                    false,
                )?
                .ok_or_else(|| {
                    truncated_error(
                        &self.source_label,
                        Some(self.frame_cursor),
                        "XTC indexed frame disappeared during preflight",
                    )
                })?;
            if refreshed.start != info.start
                || refreshed.end != info.end
                || refreshed.magic != info.magic
                || refreshed.atom_count != info.atom_count
                || refreshed.step != info.step
                || refreshed.time != info.time
            {
                return Err(corrupt_error(
                    &self.source_label,
                    self.frame_cursor,
                    "XTC repeated preflight metadata changed",
                ));
            }
        }
        self.adapter
            .decode(&info, &self.source_label, self.frame_cursor)?;
        for (point, values) in self
            .positions
            .iter_mut()
            .zip(self.adapter.decoded_positions.chunks_exact(3))
        {
            *point = Point3::new(
                f64::from(values[0]),
                f64::from(values[1]),
                f64::from(values[2]),
            );
        }
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            resource_error(
                &self.source_label,
                Some(self.frame_cursor),
                "XTC frame cursor overflows",
            )
        })?;
        Ok(Some(info))
    }

    fn publish(
        &self,
        positions: &[Point3],
        info: &XtcFrameInfo,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        let mut data = FrameBufferData::new(
            self.binding.topology_arc(),
            Quantity::new(positions, NANOMETER),
        )
        .with_time(Quantity::new(info.time, PICOSECOND))
        .with_step(info.step);
        if let Some(cell) = info.cell {
            data = data.with_cell(cell);
        }
        destination.replace_from_data(data).map_err(Into::into)
    }

    pub fn to_indexed(mut self) -> Result<IndexedXtcReader<R>, TrajectoryError> {
        let mut offsets = Vec::new();
        loop {
            if let Some(limit) = projected_index_limit(offsets.len(), &self.limits) {
                let clean_eof = if self.pending_info.is_some() {
                    false
                } else {
                    probe_seekable_eof(
                        &mut self.adapter.reader,
                        TrajectoryIoOperation::Index,
                        TrajectoryFormat::Xtc,
                        &self.source_label,
                    )?
                };
                if clean_eof {
                    break;
                }
                return Err(resource_error(
                    &self.source_label,
                    Some(offsets.len() as u64),
                    format!("XTC index {limit} exceeds the configured limit"),
                ));
            }
            let offset = self.pending_info.as_ref().map_or_else(
                || {
                    self.adapter.reader.stream_position().map_err(|error| {
                        io_context(
                            TrajectoryIoOperation::Index,
                            Some(TrajectoryFormat::Xtc),
                            &self.source_label,
                            error,
                        )
                    })
                },
                |info| Ok(info.start),
            )?;
            if self
                .parse_next()
                .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))?
                .is_none()
            {
                break;
            }
            reserve_index_for_push(
                &mut offsets,
                &self.limits,
                TrajectoryFormat::Xtc,
                &self.source_label,
                self.frame_cursor.saturating_sub(1),
            )?;
            offsets.push(offset);
        }
        self.rewind()?;
        let atom_count = self.binding.topology().atom_count();
        let mut random_positions = Vec::new();
        random_positions
            .try_reserve_exact(atom_count)
            .map_err(|_| {
                resource_error(
                    &self.source_label,
                    None,
                    "could not reserve indexed XTC position scratch",
                )
            })?;
        random_positions.resize(atom_count, Point3::new(0.0, 0.0, 0.0));
        Ok(IndexedXtcReader {
            inner: self,
            offsets,
            random_positions,
            random_adapter_scratch: Vec::new(),
            random_decoded_positions: Vec::new(),
        })
    }

    fn rewind(&mut self) -> Result<(), TrajectoryError> {
        self.adapter
            .reader
            .seek(SeekFrom::Start(self.stream_start))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xtc),
                    &self.source_label,
                    error,
                )
            })?;
        self.frame_cursor = 0;
        self.pending_info = self.adapter.preflight(
            &self.binding,
            self.options,
            &self.limits,
            &self.source_label,
            0,
            false,
        )?;
        Ok(())
    }
}

impl<R: Read + Seek> TrajectoryReader for XtcReader<R> {
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
        let offset = self.pending_info.as_ref().map_or_else(
            || {
                self.adapter.reader.stream_position().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::ReadFrame,
                        Some(TrajectoryFormat::Xtc),
                        &self.source_label,
                        error,
                    )
                })
            },
            |info| Ok(info.start),
        )?;
        let info = self
            .parse_next()
            .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))?;
        let Some(info) = info else {
            return Ok(false);
        };
        self.publish(&self.positions, &info, destination)?;
        Ok(true)
    }
}

/// Fully decoded-and-verified indexed XTC reader.
pub struct IndexedXtcReader<R> {
    inner: XtcReader<R>,
    offsets: Vec<u64>,
    random_positions: Vec<Point3>,
    random_adapter_scratch: Vec<u8>,
    random_decoded_positions: Vec<f32>,
}

impl<R: Read + Seek> IndexedXtcReader<R> {
    pub fn topology(&self) -> &Topology {
        self.inner.topology()
    }

    pub(crate) fn first_info(&self) -> &XtcFrameInfo {
        self.inner.first_info()
    }
}

impl<R: Read + Seek> TrajectoryReader for IndexedXtcReader<R> {
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

impl<R: Read + Seek> SeekableTrajectoryReader for IndexedXtcReader<R> {
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
        let saved_offset = self
            .inner
            .adapter
            .reader
            .stream_position()
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xtc),
                    &self.inner.source_label,
                    error,
                )
            })?;
        let saved_cursor = self.inner.frame_cursor;
        let saved_pending = self.inner.pending_info.clone();
        let saved_decoded_start = self.inner.adapter.decoded_start;
        self.inner
            .adapter
            .reader
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xtc),
                    &self.inner.source_label,
                    error,
                )
            })?;
        self.inner.pending_info = None;
        self.inner.adapter.decoded_start = None;
        std::mem::swap(
            &mut self.inner.adapter.scratch,
            &mut self.random_adapter_scratch,
        );
        std::mem::swap(
            &mut self.inner.adapter.decoded_positions,
            &mut self.random_decoded_positions,
        );
        std::mem::swap(&mut self.inner.positions, &mut self.random_positions);
        self.inner.frame_cursor = index;
        let result = self
            .inner
            .parse_next()
            .map_err(|error| frame_offset_context(error, index, offset))
            .and_then(|info| info.ok_or(TrajectoryError::FrameIndexOutOfRange(index)));
        let restore = self
            .inner
            .adapter
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xtc),
                    &self.inner.source_label,
                    error,
                )
            });
        self.inner.frame_cursor = saved_cursor;
        self.inner.pending_info = saved_pending;
        self.inner.adapter.decoded_start = saved_decoded_start;
        std::mem::swap(
            &mut self.inner.adapter.scratch,
            &mut self.random_adapter_scratch,
        );
        std::mem::swap(
            &mut self.inner.adapter.decoded_positions,
            &mut self.random_decoded_positions,
        );
        std::mem::swap(&mut self.inner.positions, &mut self.random_positions);
        let info = result?;
        restore?;
        self.inner
            .publish(&self.random_positions, &info, destination)
    }
}

struct MollyWriterAdapter<W> {
    inner: molly::XTCWriter<W>,
    frame: molly::Frame,
}

impl<W: Write> MollyWriterAdapter<W> {
    fn new(writer: W, options: XtcWriteOptions) -> Self {
        Self {
            inner: molly::XTCWriter::new_with_magic(writer, options.magic.to_molly()),
            frame: molly::Frame {
                precision: options.precision,
                ..molly::Frame::default()
            },
        }
    }

    fn write(&mut self, source_label: &str, frame_index: u64) -> Result<(), TrajectoryError> {
        match catch_unwind(AssertUnwindSafe(|| self.inner.write_frame(&self.frame))) {
            Ok(Ok(())) => Ok(()),
            Ok(Err(error)) => Err(io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )),
            Err(_) => Err(TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::CorruptCompressedData,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
            )
            .with_source_label(source_label)
            .with_frame(frame_index)
            .with_detail("molly panicked while encoding prevalidated XTC data")
            .into()),
        }
    }
}

/// Strict XTC writer using explicit lossy precision.
pub struct XtcWriter<W> {
    adapter: MollyWriterAdapter<W>,
    topology: Arc<Topology>,
    options: XtcWriteOptions,
    source_label: String,
    frame_count: u64,
}

impl<W: Write> XtcWriter<W> {
    pub fn new(
        writer: W,
        topology: Arc<Topology>,
        options: XtcWriteOptions,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryError> {
        let source_label = source_label.into();
        validate_precision(options.precision, &source_label)?;
        if topology.atom_count() == 0 || topology.atom_count() > i32::MAX as usize {
            return Err(resource_error(
                &source_label,
                None,
                "XTC writer atom count must fit a positive signed 32-bit XDR integer",
            ));
        }
        if options.magic == XtcMagic::Xtc1995 && topology.atom_count() > molly::XTC_1995_MAX_NATOMS
        {
            return Err(resource_error(
                &source_label,
                None,
                "XTC 1995 atom count exceeds the audited molly limit",
            ));
        }
        let adapter = MollyWriterAdapter::new(writer, options);
        Ok(Self {
            adapter,
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
        &self.adapter.inner.file
    }

    pub fn flush(&mut self) -> std::io::Result<()> {
        self.adapter.inner.file.flush()
    }

    pub(crate) fn validate_finish(&self) -> Result<(), TrajectoryError> {
        require_nonempty_writer(self.frame_count, TrajectoryFormat::Xtc, &self.source_label)
    }

    /// Flushes and returns the completed nonempty XTC stream.
    pub fn finish(mut self) -> Result<W, TrajectoryError> {
        self.validate_finish()?;
        self.adapter.inner.file.flush().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Finish,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                error,
            )
        })?;
        Ok(self.adapter.inner.file)
    }
}

impl<W: Write> TrajectoryWriter for XtcWriter<W> {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if !std::ptr::eq(frame.topology(), self.topology()) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        for (present, field) in [
            (frame.velocities().is_some(), "velocities"),
            (frame.forces().is_some(), "forces"),
            (!frame.atom_data().is_empty(), "atom data"),
            (!frame.bond_data().is_empty(), "bond data"),
            (!frame.props().is_empty(), "properties"),
        ] {
            if present {
                return Err(writer_field_error(
                    &self.source_label,
                    self.frame_count,
                    field,
                ));
            }
        }
        let step = frame.step().ok_or_else(|| {
            codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "XTC requires explicit frame step",
            )
        })?;
        let step = i32::try_from(step).map_err(|_| {
            codec_context(
                TrajectoryCodecErrorKind::NegativeOrUnrepresentableStep,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "XTC step exceeds nonnegative signed 32-bit XDR capacity",
            )
        })?;
        self.adapter.frame.step = step as u32;
        let time = frame.time().ok_or_else(|| {
            codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "XTC requires explicit frame time",
            )
        })?;
        let time = time.value_in(PICOSECOND).map_err(|error| {
            codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                format!("XTC time unit is incompatible: {error}"),
            )
        })?;
        self.adapter.frame.time = finite_f32(time, &self.source_label, "time")?;
        let cell = frame.cell().copied().ok_or_else(|| {
            codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "XTC writer requires a periodic cell",
            )
        })?;
        if cell.periodic_axes() != [true; 3] {
            return Err(codec_context(
                TrajectoryCodecErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "XTC writer requires all three periodic axes",
            ));
        }
        let vectors = cell.vectors().value_in(NANOMETER).map_err(|error| {
            codec_context(
                TrajectoryCodecErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                format!("XTC cell unit is incompatible: {error}"),
            )
        })?;
        for (target, value) in self.adapter.frame.boxvec.iter_mut().zip(
            vectors
                .into_iter()
                .flat_map(|vector| [vector.x, vector.y, vector.z]),
        ) {
            *target = finite_f32(value, &self.source_label, "box")?;
        }
        let positions = frame.positions().values();
        let factor = positions
            .unit()
            .conversion_factor_to(NANOMETER)
            .map_err(|error| {
                codec_context(
                    TrajectoryCodecErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Xtc),
                    &self.source_label,
                    format!("XTC position unit is incompatible: {error}"),
                )
            })?;
        self.adapter.frame.positions.clear();
        self.adapter
            .frame
            .positions
            .try_reserve(
                positions
                    .value()
                    .len()
                    .checked_mul(3)
                    .ok_or_else(|| {
                        resource_error(
                            &self.source_label,
                            Some(self.frame_count),
                            "XTC writer coordinate count overflows",
                        )
                    })?
                    .saturating_sub(self.adapter.frame.positions.len()),
            )
            .map_err(|_| {
                resource_error(
                    &self.source_label,
                    Some(self.frame_count),
                    "could not reserve XTC writer coordinate scratch",
                )
            })?;
        for point in *positions.value() {
            for value in [point.x * factor, point.y * factor, point.z * factor] {
                let value = finite_f32(value, &self.source_label, "position")?;
                if self.topology.atom_count() > 9 {
                    let scaled = f64::from(value) * f64::from(self.options.precision);
                    if !scaled.is_finite()
                        || scaled.round() < f64::from(i32::MIN)
                        || scaled.round() > f64::from(i32::MAX)
                    {
                        return Err(codec_context(
                            TrajectoryCodecErrorKind::InvalidFrame,
                            TrajectoryIoOperation::WriteFrame,
                            Some(TrajectoryFormat::Xtc),
                            &self.source_label,
                            "XTC coordinate exceeds the selected precision's integer range",
                        ));
                    }
                }
                self.adapter.frame.positions.push(value);
            }
        }
        self.adapter.frame.precision = self.options.precision;
        self.adapter.write(&self.source_label, self.frame_count)?;
        self.frame_count = self
            .frame_count
            .checked_add(1)
            .ok_or_else(|| resource_error(&self.source_label, None, "XTC frame count overflows"))?;
        Ok(())
    }
}

struct CheckedByteBuffer<'a> {
    bytes: &'a [u8],
    index: usize,
}

impl<'a> CheckedByteBuffer<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, index: 0 }
    }

    fn pop(&mut self, source_label: &str, frame_index: u64) -> Result<u8, TrajectoryError> {
        let value = self.bytes.get(self.index).copied().ok_or_else(|| {
            corrupt_error(
                source_label,
                frame_index,
                "XTC compressed bitstream ends before all atoms are decoded",
            )
        })?;
        self.index += 1;
        Ok(value)
    }
}

#[derive(Default)]
struct CheckedDecodeState {
    last_bits: usize,
    last_byte: u8,
}

struct CompressedLayout {
    atom_count: usize,
    minimum: [i32; 3],
    maximum: [i32; 3],
    initial_small_index: usize,
    precision: f32,
}

fn decode_compressed_payload(
    payload: &[u8],
    layout: CompressedLayout,
    output: &mut Vec<f32>,
    source_label: &str,
    frame_index: u64,
) -> Result<(), TrajectoryError> {
    let CompressedLayout {
        atom_count,
        minimum,
        maximum,
        initial_small_index,
        precision,
    } = layout;
    let mut sizes = [0_u32; 3];
    let mut axis_bits = [0_u32; 3];
    for axis in 0..3 {
        let span = maximum[axis].checked_sub(minimum[axis]).ok_or_else(|| {
            corrupt_error(
                source_label,
                frame_index,
                "XTC compressed coordinate span exceeds the audited decoder profile",
            )
        })?;
        sizes[axis] = u32::try_from(span)
            .ok()
            .and_then(|span| span.checked_add(1))
            .ok_or_else(|| {
                corrupt_error(
                    source_label,
                    frame_index,
                    "XTC compressed coordinate size overflows",
                )
            })?;
    }
    let combined_bits = if sizes.iter().fold(0_u32, |bits, size| bits | size) > 0x00ff_ffff {
        for (target, size) in axis_bits.iter_mut().zip(sizes) {
            *target = 32 - size.leading_zeros();
            if *target == 32 {
                return Err(corrupt_error(
                    source_label,
                    frame_index,
                    "XTC compressed axis requires a 32-bit scalar unsupported by the audited adapter",
                ));
            }
        }
        0
    } else {
        combined_bit_count(sizes)
    };
    let mut buffer = CheckedByteBuffer::new(payload);
    let mut state = CheckedDecodeState::default();
    let mut small_index = initial_small_index;
    let mut smaller = XTC_MAGIC_INTS[small_index.saturating_sub(1).max(XTC_FIRST_SMALL_INDEX)] / 2;
    let mut small_number = XTC_MAGIC_INTS[small_index] / 2;
    let mut read_atoms = 0_usize;
    let inverse_precision = precision.recip();

    while read_atoms < atom_count {
        let mut coordinate = if combined_bits == 0 {
            let mut decoded = [0_i32; 3];
            for axis in 0..3 {
                let value = decode_checked_bits(
                    &mut buffer,
                    &mut state,
                    axis_bits[axis],
                    source_label,
                    frame_index,
                )?;
                if value >= sizes[axis] || value > i32::MAX as u32 {
                    return Err(corrupt_error(
                        source_label,
                        frame_index,
                        "XTC compressed absolute coordinate is outside declared bounds",
                    ));
                }
                decoded[axis] = value as i32;
            }
            decoded
        } else {
            decode_checked_triplet(
                &mut buffer,
                &mut state,
                combined_bits,
                sizes,
                source_label,
                frame_index,
            )?
        };
        for axis in 0..3 {
            coordinate[axis] = minimum[axis]
                .checked_add(coordinate[axis])
                .filter(|value| *value <= maximum[axis])
                .ok_or_else(|| {
                    corrupt_error(
                        source_label,
                        frame_index,
                        "XTC compressed absolute coordinate overflows declared bounds",
                    )
                })?;
        }
        let mut previous = coordinate;
        let flag = decode_checked_bits(&mut buffer, &mut state, 1, source_label, frame_index)? != 0;
        let mut smaller_change = 0_i32;
        let mut run = 0_u32;
        if flag {
            run = decode_checked_bits(&mut buffer, &mut state, 5, source_label, frame_index)?;
            let remainder = run % 3;
            run -= remainder;
            smaller_change = remainder as i32 - 1;
        }
        if run == 0 {
            push_checked_coordinate(
                output,
                coordinate,
                inverse_precision,
                source_label,
                frame_index,
            )?;
            read_atoms += 1;
        } else {
            let emitted_atoms = usize::try_from(run / 3)
                .ok()
                .and_then(|count| count.checked_add(1))
                .ok_or_else(|| {
                    corrupt_error(
                        source_label,
                        frame_index,
                        "XTC compressed run atom count overflows",
                    )
                })?;
            if emitted_atoms > atom_count - read_atoms {
                return Err(corrupt_error(
                    source_label,
                    frame_index,
                    "XTC compressed run exceeds the declared atom count",
                ));
            }
            let small_size = u32::try_from(XTC_MAGIC_INTS[small_index])
                .expect("audited XTC magic integers are positive");
            let small_sizes = [small_size; 3];
            for offset in (0..run).step_by(3) {
                let delta = decode_checked_triplet(
                    &mut buffer,
                    &mut state,
                    small_index as u32,
                    small_sizes,
                    source_label,
                    frame_index,
                )?;
                for axis in 0..3 {
                    coordinate[axis] = previous[axis]
                        .checked_sub(small_number)
                        .and_then(|base| base.checked_add(delta[axis]))
                        .filter(|value| *value >= minimum[axis] && *value <= maximum[axis])
                        .ok_or_else(|| {
                            corrupt_error(
                                source_label,
                                frame_index,
                                "XTC compressed run coordinate overflows declared bounds",
                            )
                        })?;
                }
                if offset == 0 {
                    std::mem::swap(&mut coordinate, &mut previous);
                    push_checked_coordinate(
                        output,
                        previous,
                        inverse_precision,
                        source_label,
                        frame_index,
                    )?;
                    read_atoms += 1;
                } else {
                    previous = coordinate;
                }
                push_checked_coordinate(
                    output,
                    coordinate,
                    inverse_precision,
                    source_label,
                    frame_index,
                )?;
                read_atoms += 1;
            }
        }

        match smaller_change.cmp(&0) {
            std::cmp::Ordering::Less => {
                if small_index <= XTC_FIRST_SMALL_INDEX {
                    return Err(corrupt_error(
                        source_label,
                        frame_index,
                        "XTC compressed small-index would underflow",
                    ));
                }
                small_index -= 1;
                small_number = smaller;
                smaller = if small_index > XTC_FIRST_SMALL_INDEX {
                    XTC_MAGIC_INTS[small_index - 1] / 2
                } else {
                    0
                };
            }
            std::cmp::Ordering::Greater => {
                if small_index + 1 >= XTC_MAGIC_INTS.len() {
                    return Err(corrupt_error(
                        source_label,
                        frame_index,
                        "XTC compressed small-index would overflow",
                    ));
                }
                small_index += 1;
                smaller = small_number;
                small_number = XTC_MAGIC_INTS[small_index] / 2;
            }
            std::cmp::Ordering::Equal => {}
        }
    }
    validate_decoded_count(output, atom_count, source_label, frame_index)?;
    if buffer.index != payload.len() {
        return Err(corrupt_error(
            source_label,
            frame_index,
            "XTC compressed payload contains unused trailing bytes",
        ));
    }
    let trailing_mask = if state.last_bits == 0 {
        0
    } else {
        (1_u16 << state.last_bits) - 1
    };
    if u16::from(state.last_byte) & trailing_mask != 0 {
        return Err(corrupt_error(
            source_label,
            frame_index,
            "XTC compressed payload contains nonzero trailing bits",
        ));
    }
    Ok(())
}

fn decode_checked_bits(
    buffer: &mut CheckedByteBuffer<'_>,
    state: &mut CheckedDecodeState,
    count: u32,
    source_label: &str,
    frame_index: u64,
) -> Result<u32, TrajectoryError> {
    if count > 32 {
        return Err(corrupt_error(
            source_label,
            frame_index,
            "XTC compressed scalar requests more than 32 bits",
        ));
    }
    let mask = if count == 32 {
        u32::MAX
    } else {
        (1_u32 << count) - 1
    };
    let mut remaining = count as usize;
    let mut last_bits = state.last_bits;
    let mut last_byte = u32::from(state.last_byte);
    let mut value = 0_u32;
    while remaining >= 8 {
        last_byte = last_byte.wrapping_shl(8) | u32::from(buffer.pop(source_label, frame_index)?);
        value |= (last_byte >> last_bits) << (remaining - 8);
        remaining -= 8;
    }
    if remaining > 0 {
        if last_bits < remaining {
            last_bits += 8;
            last_byte =
                last_byte.wrapping_shl(8) | u32::from(buffer.pop(source_label, frame_index)?);
        }
        last_bits -= remaining;
        value |= (last_byte >> last_bits) & mask;
    }
    state.last_bits = last_bits;
    state.last_byte = (last_byte & 0xff) as u8;
    Ok(value & mask)
}

fn push_checked_coordinate(
    output: &mut Vec<f32>,
    coordinate: [i32; 3],
    inverse_precision: f32,
    source_label: &str,
    frame_index: u64,
) -> Result<(), TrajectoryError> {
    for value in coordinate {
        let value = value as f32 * inverse_precision;
        if !value.is_finite() {
            return Err(corrupt_error(
                source_label,
                frame_index,
                "XTC decoded coordinate is not finite",
            ));
        }
        output.push(value);
    }
    Ok(())
}

fn validate_decoded_count(
    output: &[f32],
    atom_count: usize,
    source_label: &str,
    frame_index: u64,
) -> Result<(), TrajectoryError> {
    if output.len() != atom_count.saturating_mul(3) {
        return Err(corrupt_error(
            source_label,
            frame_index,
            "XTC compressed bitstream did not produce the declared atom count",
        ));
    }
    Ok(())
}

fn combined_bit_count(sizes: [u32; 3]) -> u32 {
    let mut byte_count = 1_usize;
    let mut bytes = [0_u8; 32];
    bytes[0] = 1;
    for size in sizes {
        let mut carry = 0_u32;
        for byte in &mut bytes[..byte_count] {
            carry += u32::from(*byte) * size;
            *byte = (carry & 0xff) as u8;
            carry >>= 8;
        }
        while carry != 0 {
            bytes[byte_count] = (carry & 0xff) as u8;
            byte_count += 1;
            carry >>= 8;
        }
    }
    let high = bytes[byte_count - 1];
    (byte_count as u32 - 1) * 8 + (8 - high.leading_zeros())
}

fn decode_checked_triplet(
    buffer: &mut CheckedByteBuffer<'_>,
    state: &mut CheckedDecodeState,
    mut count: u32,
    sizes: [u32; 3],
    source_label: &str,
    frame_index: u64,
) -> Result<[i32; 3], TrajectoryError> {
    let mut encoded = [0_u8; 32];
    let mut bytes = 0_usize;
    while count >= 8 {
        encoded[bytes] = decode_checked_bits(buffer, state, 8, source_label, frame_index)? as u8;
        bytes += 1;
        count -= 8;
    }
    if count > 0 {
        encoded[bytes] =
            decode_checked_bits(buffer, state, count, source_label, frame_index)? as u8;
        bytes += 1;
    }
    let mut decoded = [0_u32; 3];
    for axis in (1..=2).rev() {
        let mut remainder = 0_u64;
        for index in (0..bytes).rev() {
            remainder = (remainder << 8) | u64::from(encoded[index]);
            let quotient = remainder / u64::from(sizes[axis]);
            if quotient > u64::from(u8::MAX) {
                return Err(corrupt_error(
                    source_label,
                    frame_index,
                    "XTC compressed mixed-radix quotient overflows",
                ));
            }
            encoded[index] = quotient as u8;
            remainder -= quotient * u64::from(sizes[axis]);
        }
        decoded[axis] = remainder as u32;
    }
    if bytes > 4 && encoded[4..bytes].iter().any(|byte| *byte != 0) {
        return Err(corrupt_error(
            source_label,
            frame_index,
            "XTC compressed x coordinate overflows",
        ));
    }
    let mut x = [0_u8; 4];
    let copied = bytes.min(x.len());
    x[..copied].copy_from_slice(&encoded[..copied]);
    decoded[0] = u32::from_le_bytes(x);
    for axis in 0..3 {
        if decoded[axis] >= sizes[axis] || decoded[axis] > i32::MAX as u32 {
            return Err(corrupt_error(
                source_label,
                frame_index,
                "XTC compressed coordinate is outside its mixed-radix bounds",
            ));
        }
    }
    Ok(decoded.map(|value| value as i32))
}

fn decode_cell(
    values: [f64; 9],
    policy: XtcCellPolicy,
    source_label: &str,
    frame: u64,
) -> Result<Option<PeriodicCell>, TrajectoryError> {
    if values.iter().any(|value| !value.is_finite()) {
        return Err(frame_error(
            source_label,
            frame,
            "XTC box contains a non-finite value",
        ));
    }
    if values.iter().all(|value| *value == 0.0) {
        return match policy {
            XtcCellPolicy::ZeroMatrixAsAbsent => Ok(None),
            XtcCellPolicy::RequirePeriodic => Err(codec_context(
                TrajectoryCodecErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Xtc),
                source_label,
                "XTC zero box requires the explicit ZeroMatrixAsAbsent policy",
            )),
        };
    }
    let vectors = [
        Vector3::new(values[0], values[1], values[2]),
        Vector3::new(values[3], values[4], values[5]),
        Vector3::new(values[6], values[7], values[8]),
    ];
    PeriodicCell::new(Quantity::new(vectors, NANOMETER), [true; 3])
        .map(Some)
        .map_err(|error| {
            frame_error(
                source_label,
                frame,
                format!("XTC periodic cell is invalid: {error}"),
            )
        })
}

fn read_nbytes<R: Read>(
    reader: &mut R,
    magic: XtcMagic,
    source_label: &str,
    frame: u64,
) -> Result<usize, TrajectoryError> {
    let mut bytes = [0_u8; 8];
    read_exact(
        reader,
        &mut bytes[..magic.nbytes_width()],
        source_label,
        Some(frame),
        "truncated XTC compressed-byte count",
    )?;
    let value = match magic {
        XtcMagic::Xtc1995 => u64::from(u32::from_be_bytes(bytes[..4].try_into().expect("u32"))),
        XtcMagic::Xtc2023 => u64::from_be_bytes(bytes),
    };
    usize::try_from(value).map_err(|_| {
        resource_error(
            source_label,
            Some(frame),
            "XTC byte count does not fit usize",
        )
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
            truncated_error(source_label, frame, detail)
        } else {
            io_context(
                if frame.is_some() {
                    TrajectoryIoOperation::ReadFrame
                } else {
                    TrajectoryIoOperation::ReadHeader
                },
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )
        }
    })
}

fn validate_atom_count(
    atom_count: usize,
    limits: &TrajectoryIoLimits,
    source_label: &str,
) -> Result<(), TrajectoryError> {
    if atom_count == 0 || atom_count > limits.max_atoms || atom_count > u32::MAX as usize {
        return Err(resource_error(
            source_label,
            None,
            "XTC atom count is zero or exceeds configured/profile limits",
        ));
    }
    Ok(())
}

fn nonnegative_xdr_count(
    value: i32,
    field: &str,
    source_label: &str,
    frame: u64,
    byte_offset: u64,
) -> Result<usize, TrajectoryError> {
    usize::try_from(value).map_err(|_| {
        TrajectoryCodecErrorContext::new(
            TrajectoryCodecErrorKind::InvalidHeader,
            TrajectoryIoOperation::ReadHeader,
            Some(TrajectoryFormat::Xtc),
        )
        .with_source_label(source_label)
        .with_frame(frame)
        .with_byte_offset(byte_offset)
        .with_detail(format!("XTC {field} is a negative signed XDR integer"))
        .into()
    })
}

fn validate_precision(precision: f32, source_label: &str) -> Result<(), TrajectoryError> {
    if !precision.is_finite() || precision <= 0.0 {
        return Err(codec_context(
            TrajectoryCodecErrorKind::InvalidPrecision,
            TrajectoryIoOperation::Open,
            Some(TrajectoryFormat::Xtc),
            source_label,
            "XTC precision must be finite and positive",
        ));
    }
    Ok(())
}

fn finite_f32(value: f64, source_label: &str, field: &str) -> Result<f32, TrajectoryError> {
    let narrowed = value as f32;
    if !value.is_finite() || !narrowed.is_finite() {
        return Err(codec_context(
            TrajectoryCodecErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Xtc),
            source_label,
            format!("XTC {field} cannot be represented as finite f32"),
        ));
    }
    Ok(narrowed)
}

fn header_error(source_label: &str, detail: impl Into<String>) -> TrajectoryError {
    codec_context(
        TrajectoryCodecErrorKind::InvalidHeader,
        TrajectoryIoOperation::ReadHeader,
        Some(TrajectoryFormat::Xtc),
        source_label,
        detail,
    )
}

fn frame_error(source_label: &str, frame: u64, detail: impl Into<String>) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::InvalidFrame,
        TrajectoryIoOperation::ReadFrame,
        Some(TrajectoryFormat::Xtc),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_detail(detail)
    .into()
}

fn corrupt_error(source_label: &str, frame: u64, detail: impl Into<String>) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::CorruptCompressedData,
        TrajectoryIoOperation::ReadFrame,
        Some(TrajectoryFormat::Xtc),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_detail(detail)
    .into()
}

fn truncated_error(
    source_label: &str,
    frame: Option<u64>,
    detail: impl Into<String>,
) -> TrajectoryError {
    let mut context = TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::TruncatedRecord,
        if frame.is_some() {
            TrajectoryIoOperation::ReadFrame
        } else {
            TrajectoryIoOperation::ReadHeader
        },
        Some(TrajectoryFormat::Xtc),
    )
    .with_source_label(source_label)
    .with_detail(detail);
    if let Some(frame) = frame {
        context = context.with_frame(frame);
    }
    context.into()
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
        Some(TrajectoryFormat::Xtc),
    )
    .with_source_label(source_label)
    .with_detail(detail);
    if let Some(frame) = frame {
        context = context.with_frame(frame);
    }
    context.into()
}

fn writer_field_error(source_label: &str, frame: u64, field: &str) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(
        TrajectoryCodecErrorKind::UnsupportedField,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Xtc),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_detail(format!("XTC cannot preserve {field}"))
    .into()
}
