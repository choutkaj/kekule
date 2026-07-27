//! Defensive private adapter over `molly` for GROMACS XTC trajectories.

use std::io::{Read, Seek, SeekFrom, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};

use molecular::geometry::{PeriodicCell, Point3, Vector3};
use molecular::topology::Topology;
use molecular::trajectory::{
    FrameBuffer, FrameBufferData, SeekableTrajectoryReader, TrajectoryCodecErrorContext,
    TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat, TrajectoryFrameView,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use molecular::units::{Quantity, NANOMETER, PICOSECOND};

use crate::{codec_context, io_context, TrajectoryIoLimits, TrajectoryTopologyBinding};

const XTC_HEADER_BYTES: usize = 56;
const XTC_COMPRESSED_PRELUDE_BYTES: usize = 32;

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

    pub(crate) const fn has_cell(&self) -> bool {
        self.cell.is_some()
    }

    pub(crate) const fn compressed_bytes(&self) -> usize {
        self.compressed_bytes
    }
}

struct MollyReaderAdapter<R> {
    inner: molly::XTCReader<R>,
    frame: molly::Frame,
    scratch: Vec<u8>,
}

impl<R: Read + Seek> MollyReaderAdapter<R> {
    fn new(reader: R) -> Self {
        Self {
            inner: molly::XTCReader::new(reader),
            frame: molly::Frame::default(),
            scratch: Vec::new(),
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
        let start = self.inner.file.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )
        })?;
        let mut header = [0_u8; XTC_HEADER_BYTES];
        match self.inner.file.read(&mut header[..1]) {
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
            &mut self.inner.file,
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
        let atom_count =
            usize::try_from(u32::from_be_bytes(header[4..8].try_into().expect("natoms")))
                .map_err(|_| resource_error(source_label, None, "XTC atom count does not fit"))?;
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
        let repeated = usize::try_from(u32::from_be_bytes(
            header[52..56].try_into().expect("repeat"),
        ))
        .map_err(|_| resource_error(source_label, None, "XTC repeated count does not fit"))?;
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
        let step = u64::from(u32::from_be_bytes(header[8..12].try_into().expect("step")));
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
        let (precision, compressed_bytes, end) = if atom_count <= 9 {
            let coordinate_bytes = atom_count
                .checked_mul(12)
                .ok_or_else(|| resource_error(source_label, None, "XTC small frame overflows"))?;
            let end = start
                .checked_add(XTC_HEADER_BYTES as u64)
                .and_then(|offset| offset.checked_add(coordinate_bytes as u64))
                .ok_or_else(|| resource_error(source_label, None, "XTC frame end overflows"))?;
            (None, coordinate_bytes, end)
        } else {
            let mut prelude = [0_u8; XTC_COMPRESSED_PRELUDE_BYTES];
            read_exact(
                &mut self.inner.file,
                &mut prelude,
                source_label,
                Some(frame_index),
                "truncated XTC compressed prelude",
            )?;
            let precision = f32::from_be_bytes(prelude[..4].try_into().expect("precision"));
            validate_precision(precision, source_label)?;
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
                if min > max {
                    return Err(corrupt_error(
                        source_label,
                        frame_index,
                        "XTC compressed coordinate bounds are reversed",
                    ));
                }
            }
            let small_index = u32::from_be_bytes(prelude[28..32].try_into().expect("small index"));
            if !(9..73).contains(&small_index) {
                return Err(corrupt_error(
                    source_label,
                    frame_index,
                    "XTC compressed small-index is outside the audited table",
                ));
            }
            let compressed_bytes =
                read_nbytes(&mut self.inner.file, magic, source_label, frame_index)?;
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
        let base_scratch = atom_count
            .checked_mul(3 * std::mem::size_of::<f32>() + std::mem::size_of::<Point3>())
            .ok_or_else(|| resource_error(source_label, None, "XTC scratch size overflows"))?;
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
        let after_header = self.inner.file.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )
        })?;
        let file_end = self.inner.file.seek(SeekFrom::End(0)).map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )
        })?;
        if end > file_end {
            self.inner
                .file
                .seek(SeekFrom::Start(start))
                .map_err(|error| {
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
                self.inner
                    .file
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
                    &mut self.inner.file,
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
        let _ = after_header;
        self.inner
            .file
            .seek(SeekFrom::Start(start))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadHeader,
                    Some(TrajectoryFormat::Xtc),
                    source_label,
                    error,
                )
            })?;
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
        self.inner
            .file
            .seek(SeekFrom::Start(info.start))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::ReadFrame,
                    Some(TrajectoryFormat::Xtc),
                    source_label,
                    error,
                )
            })?;
        let result = catch_unwind(AssertUnwindSafe(|| {
            self.inner.read_frame_with_scratch(
                &mut self.frame,
                &mut self.scratch,
                &molly::selection::AtomSelection::All,
            )
        }));
        match result {
            Ok(Ok(())) => {}
            Ok(Err(error)) => {
                return Err(corrupt_error(
                    source_label,
                    frame_index,
                    format!("molly rejected XTC compressed data: {error}"),
                ))
            }
            Err(_) => {
                return Err(corrupt_error(
                    source_label,
                    frame_index,
                    "molly panicked while decoding bounded XTC data",
                ))
            }
        }
        let position = self.inner.file.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Xtc),
                source_label,
                error,
            )
        })?;
        if position != info.end
            || self.frame.positions.len() != info.atom_count.saturating_mul(3)
            || u64::from(self.frame.step) != info.step
            || f64::from(self.frame.time) != info.time
        {
            return Err(corrupt_error(
                source_label,
                frame_index,
                "XTC adapter output disagrees with bounded preflight metadata",
            ));
        }
        if self.frame.positions.iter().any(|value| !value.is_finite()) {
            return Err(corrupt_error(
                source_label,
                frame_index,
                "XTC decoder produced a non-finite coordinate",
            ));
        }
        if let Some(precision) = info.precision {
            if self.frame.precision != precision {
                return Err(corrupt_error(
                    source_label,
                    frame_index,
                    "XTC decoded precision disagrees with preflight",
                ));
            }
        }
        Ok(())
    }
}

/// Sequential XTC reader using the private defensive molly adapter.
pub struct XtcReader<R> {
    adapter: MollyReaderAdapter<R>,
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
        let mut adapter = MollyReaderAdapter::new(reader);
        let stream_start = adapter.inner.file.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Xtc),
                &source_label,
                error,
            )
        })?;
        let first_info = adapter
            .preflight(&binding, options, &limits, &source_label, 0, false)?
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
        self.adapter.preflight(
            &self.binding,
            self.options,
            &self.limits,
            &self.source_label,
            self.frame_cursor,
            true,
        )
    }

    fn parse_next(&mut self, publish: Option<&mut FrameBuffer>) -> Result<bool, TrajectoryError> {
        let Some(info) = self.next_info()? else {
            return Ok(false);
        };
        if self.frame_cursor >= self.limits.max_frames {
            return Err(resource_error(
                &self.source_label,
                Some(self.frame_cursor),
                "XTC frame count exceeds the configured limit",
            ));
        }
        self.adapter
            .decode(&info, &self.source_label, self.frame_cursor)?;
        for (point, values) in self
            .positions
            .iter_mut()
            .zip(self.adapter.frame.positions.chunks_exact(3))
        {
            *point = Point3::new(
                f64::from(values[0]),
                f64::from(values[1]),
                f64::from(values[2]),
            );
        }
        if let Some(destination) = publish {
            let mut data = FrameBufferData::new(
                self.topology(),
                Quantity::new(self.positions.as_slice(), NANOMETER),
            )
            .with_time(Quantity::new(info.time, PICOSECOND))
            .with_step(info.step);
            if let Some(cell) = info.cell {
                data = data.with_cell(cell);
            }
            destination.replace_from_data(data)?;
        }
        self.frame_cursor += 1;
        Ok(true)
    }

    pub fn into_indexed(mut self) -> Result<IndexedXtcReader<R>, TrajectoryError> {
        let mut offsets = Vec::new();
        loop {
            let offset = self.pending_info.as_ref().map_or_else(
                || {
                    self.adapter.inner.file.stream_position().map_err(|error| {
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
            if !self.parse_next(None)? {
                break;
            }
            if offsets.len() >= self.limits.max_index_entries {
                return Err(resource_error(
                    &self.source_label,
                    Some(offsets.len() as u64),
                    "XTC index entry limit exceeded",
                ));
            }
            offsets.try_reserve(1).map_err(|_| {
                resource_error(
                    &self.source_label,
                    Some(offsets.len() as u64),
                    "could not grow XTC index",
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
                    "XTC index byte limit exceeded",
                ));
            }
        }
        self.rewind()?;
        Ok(IndexedXtcReader {
            inner: self,
            offsets,
        })
    }

    fn rewind(&mut self) -> Result<(), TrajectoryError> {
        self.adapter
            .inner
            .file
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

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        self.parse_next(Some(destination))
    }
}

/// Fully decoded-and-verified indexed XTC reader.
pub struct IndexedXtcReader<R> {
    inner: XtcReader<R>,
    offsets: Vec<u64>,
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
        let offset = self
            .offsets
            .get(usize::try_from(index).map_err(|_| TrajectoryError::FrameIndexOutOfRange(index))?)
            .copied()
            .ok_or(TrajectoryError::FrameIndexOutOfRange(index))?;
        let saved_offset = self
            .inner
            .adapter
            .inner
            .file
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
        let saved_pending = self.inner.pending_info.take();
        let saved_adapter_step = self.inner.adapter.inner.step;
        self.inner
            .adapter
            .inner
            .file
            .seek(SeekFrom::Start(offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Xtc),
                    &self.inner.source_label,
                    error,
                )
            })?;
        self.inner.frame_cursor = index;
        self.inner.adapter.inner.step = usize::try_from(index).unwrap_or(usize::MAX);
        let result = self.inner.parse_next(Some(destination)).and_then(|read| {
            if read {
                Ok(())
            } else {
                Err(TrajectoryError::FrameIndexOutOfRange(index))
            }
        });
        let restore = self
            .inner
            .adapter
            .inner
            .file
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
        self.inner.adapter.inner.step = saved_adapter_step;
        result.and(restore.map(|_| ()))
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
    topology: Topology,
    options: XtcWriteOptions,
    source_label: String,
    frame_count: u64,
}

impl<W: Write> XtcWriter<W> {
    pub fn new(
        writer: W,
        topology: Topology,
        options: XtcWriteOptions,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryError> {
        let source_label = source_label.into();
        validate_precision(options.precision, &source_label)?;
        if topology.atom_count() == 0 || topology.atom_count() > u32::MAX as usize {
            return Err(resource_error(
                &source_label,
                None,
                "XTC writer atom count must fit a positive u32",
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

    pub fn finish(mut self) -> Result<W, TrajectoryError> {
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

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if frame.topology().identity() != self.topology.identity() {
            return Err(TrajectoryError::TopologyIdentityMismatch);
        }
        for (present, field) in [
            (frame.velocities().is_some(), "velocities"),
            (frame.forces().is_some(), "forces"),
            (frame.observation().is_some(), "observation"),
            (!frame.props().is_empty(), "properties"),
        ] {
            if present {
                return Err(TrajectoryError::UnsupportedField(field));
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
        self.adapter.frame.step = u32::try_from(step).map_err(|_| {
            codec_context(
                TrajectoryCodecErrorKind::NegativeOrUnrepresentableStep,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Xtc),
                &self.source_label,
                "XTC step exceeds u32 capacity",
            )
        })?;
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
        let cell = frame.configuration().cell().copied().ok_or_else(|| {
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
        let positions = frame.configuration().positions().values();
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
                    .saturating_mul(3)
                    .saturating_sub(self.adapter.frame.positions.capacity()),
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
