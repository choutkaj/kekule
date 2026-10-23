//! Strict common-profile CHARMM/NAMD/OpenMM DCD trajectory I/O.
//!
//! [`DcdReader`] and [`IndexedDcdReader`] read 32-bit-record `CORD` files in
//! either byte order, with optional unit cells and fixed atoms. [`DcdWriter`]
//! writes canonical all-atom files. DCD stores f32 coordinates in angstrom and
//! one step per frame on the `ISTART`/`NSAVC` sequence; time is derived from
//! the header `DELTA` only under an explicit [`DcdTimePolicy`].
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
//! use kekule_traj::io::dcd::{DcdReadOptions, DcdReader, DcdWriteOptions, DcdWriter};
//! use kekule_traj::{TrajectoryReader, TrajectoryWriter};
//!
//! let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//! let mut trajectory = Trajectory::new(Arc::clone(&topology));
//! for step in [10, 20] {
//!     let x = step as f64;
//!     let points = vec![Point3::new(x, 0.0, 0.0)];
//!     let positions = Positions::from_vec(Quantity::new(points, ANGSTROM))?;
//!     let mut frame = TrajectoryFrame::new(positions);
//!     frame.set_step(Some(step));
//!     trajectory.push(frame)?;
//! }
//!
//! // Frame steps must follow the writer's step sequence exactly.
//! let options = DcdWriteOptions::default().with_step_sequence(10, 10);
//! let stream = Cursor::new(Vec::new());
//! let mut writer = DcdWriter::new(stream, Arc::clone(&topology), options, "memory")?;
//! for frame in &trajectory {
//!     writer.write_frame(frame)?;
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

use crate::TrajectoryField;
use std::fmt;
use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::Arc;

use crate::{
    FrameBuffer, FrameData, SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoError,
    TrajectoryIoErrorKind, TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::structure::TrajectoryFrameView;
use kekule::topology::Topology;
use kekule::units::{Quantity, Unit, ANGSTROM, CANONICAL_LENGTH_UNIT, CANONICAL_TIME_UNIT};

use super::limits::{reserve_scratch, scratch_bytes};
use super::{
    codec_context, frame_offset_context, index_mismatch, io_context, probe_seekable_eof,
    projected_index_limit, require_no_partial_frame, require_nonempty_writer, validation_detail,
    FrameOffsets, IndexBuildControl, TrajectoryIoLimits,
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

/// Whether a DCD writer stores periodic cells.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum DcdCellPolicy {
    /// Write no unit-cell records; a frame with a cell is rejected.
    #[default]
    Omit,
    /// Write one unit-cell record per frame; a frame without a cell is rejected.
    Write,
}

/// DCD reader policy.
#[derive(Debug, Clone, PartialEq)]
pub struct DcdReadOptions {
    time_policy: DcdTimePolicy,
    limits: TrajectoryIoLimits,
    source_label: String,
}

impl Default for DcdReadOptions {
    fn default() -> Self {
        Self {
            limits: TrajectoryIoLimits::default(),
            source_label: "dcd stream".into(),
            time_policy: DcdTimePolicy::Absent,
        }
    }
}

impl DcdReadOptions {
    /// Sets limits applied before allocation, scanning, and seeking.
    ///
    /// The scratch limit covers the combined retained capacities of coordinates,
    /// record bytes, and fixed-atom reconstruction/validation buffers. Indexed
    /// and sequential reads share those buffers; index offsets have a separate limit.
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

    /// Sets [`Self::time_policy`].
    #[must_use]
    pub const fn with_time_policy(mut self, time_policy: DcdTimePolicy) -> Self {
        self.time_policy = time_policy;
        self
    }

    /// How frame times are derived; [`DcdTimePolicy::Absent`] by default.
    pub const fn time_policy(&self) -> DcdTimePolicy {
        self.time_policy
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

/// Canonical DCD writer policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DcdWriteOptions {
    endian: DcdEndian,
    cell_policy: DcdCellPolicy,
    start_step: u64,
    step_interval: u64,
    delta: f32,
    time_policy: DcdTimePolicy,
}

impl Default for DcdWriteOptions {
    fn default() -> Self {
        Self {
            endian: DcdEndian::Little,
            cell_policy: DcdCellPolicy::Omit,
            start_step: 0,
            step_interval: 1,
            delta: 1.0,
            time_policy: DcdTimePolicy::Absent,
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

    /// Includes cells in the canonical DCD orientation: `a` along positive x,
    /// `b` in the xy plane with positive y, and `c` with positive z.
    /// Other orientations are rejected because DCD lengths and angles cannot
    /// preserve their relation to the stored Cartesian coordinates.
    /// Angles are written using the unambiguous supported cosine convention.
    #[must_use]
    pub const fn with_cell_policy(mut self, cell_policy: DcdCellPolicy) -> Self {
        self.cell_policy = cell_policy;
        self
    }

    /// Sets the header step sequence: frame `n` must carry step
    /// `start_step + n * step_interval`. Both must fit an `i32`, and the
    /// interval must be positive; [`DcdWriter::new`] checks them.
    #[must_use]
    pub const fn with_step_sequence(mut self, start_step: u64, step_interval: u64) -> Self {
        self.start_step = start_step;
        self.step_interval = step_interval;
        self
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
    /// is not a time.
    pub fn with_header_delta(
        mut self,
        delta: Quantity<f64>,
        time_policy: DcdTimePolicy,
    ) -> Result<Self, TrajectoryIoError> {
        let unit = header_delta_unit(time_policy);
        let delta = delta.value_in(unit).map_err(|error| {
            codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Open,
                Some(TrajectoryFormat::Dcd),
                "DCD writer options",
                format!("DCD DELTA must be a time: {error}"),
            )
        })?;
        self.delta = delta as f32;
        self.time_policy = time_policy;
        Ok(self)
    }

    /// The byte order written; [`DcdEndian::Little`] by default.
    pub const fn endian(&self) -> DcdEndian {
        self.endian
    }

    /// Whether unit cells are written; [`DcdCellPolicy::Omit`] by default.
    pub const fn cell_policy(&self) -> DcdCellPolicy {
        self.cell_policy
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

    /// How frame times relate to `DELTA`; [`DcdTimePolicy::Absent`], which
    /// rejects frames carrying a time, by default.
    pub const fn time_policy(&self) -> DcdTimePolicy {
        self.time_policy
    }

    const fn writes_cells(&self) -> bool {
        match self.cell_policy {
            DcdCellPolicy::Omit => false,
            DcdCellPolicy::Write => true,
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

/// Sequential DCD reader over one seekable stream.
pub struct DcdReader<R> {
    reader: R,
    topology: Arc<Topology>,
    options: DcdReadOptions,
    header: DcdHeader,
    free_indices: Vec<usize>,
    fixed_indices: Vec<usize>,
    fixed_reference: Vec<Point3>,
    positions: Vec<Point3>,
    record: Vec<u8>,
    frame_cursor: usize,
}

impl<R> fmt::Debug for DcdReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DcdReader")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("frame_cursor", &self.frame_cursor)
            .finish_non_exhaustive()
    }
}

#[derive(Debug, Clone, Copy)]
struct DcdDecodedFrame {
    cell: Option<PeriodicCell>,
    step: u64,
    time: Option<f64>,
}

impl<R: Read + Seek> DcdReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order.
    ///
    /// Accepts an owned topology or a shared `Arc<Topology>`. Counts and available
    /// format metadata are checked automatically; matching counts alone cannot
    /// establish atom identity. Limits and diagnostic labels are set in `options`.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidHeader`],
    /// [`TrajectoryIoErrorKind::UnsupportedVariant`],
    /// [`TrajectoryIoErrorKind::TruncatedRecord`] or
    /// [`TrajectoryIoErrorKind::RecordMarkerMismatch`] for a malformed or
    /// unsupported header, [`TrajectoryIoErrorKind::InconsistentAtomCount`]
    /// when its atom count differs from the topology's,
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when it exceeds the limits, and
    /// [`TrajectoryIoErrorKind::Io`] when reading fails.
    pub fn new(
        mut reader: R,
        topology: impl Into<Arc<Topology>>,
        options: DcdReadOptions,
    ) -> Result<Self, TrajectoryIoError> {
        let topology = topology.into();
        let source_label = &options.source_label;
        let limits = &options.limits;
        validate_time_policy(options.time_policy, source_label)?;
        let atom_count = topology.atom_count();
        validate_atom_count(atom_count, limits, source_label)?;
        let start = reader
            .stream_position()
            .map_err(|error| io_context(TrajectoryIoOperation::Open, None, source_label, error))?;
        let mut marker = [0_u8; 4];
        read_exact_required(
            &mut reader,
            &mut marker,
            TrajectoryIoOperation::ReadHeader,
            source_label,
            None,
            "DCD header marker",
        )?;
        let endian = if u32::from_le_bytes(marker) == HEADER_BYTES as u32 {
            DcdEndian::Little
        } else if u32::from_be_bytes(marker) == HEADER_BYTES as u32 {
            DcdEndian::Big
        } else {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidHeader,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                source_label,
                "DCD first record marker is not 84 in either byte order",
            ));
        };
        reader
            .seek(SeekFrom::Start(start))
            .map_err(|error| io_context(TrajectoryIoOperation::Open, None, source_label, error))?;

        let mut record = Vec::new();
        read_record(
            &mut reader,
            endian,
            &mut record,
            Some(0),
            limits,
            source_label,
            TrajectoryIoOperation::ReadHeader,
            None,
            false,
        )?;
        if record.len() != HEADER_BYTES || &record[..4] != b"CORD" {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidHeader,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                source_label,
                "DCD header must be an 84-byte CORD record",
            ));
        }
        let controls = parse_controls(&record[4..], endian, source_label)?;
        if controls[11] != 0 || controls[12] != 0 {
            return Err(codec_context(
                TrajectoryIoErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                source_label,
                "DCD 4D or charge coordinate records are not supported",
            ));
        }
        let declared_frames = nonnegative_u64(controls[0], "NSET", source_label)? as usize;
        let start_step = nonnegative_u64(controls[1], "ISTART", source_label)?;
        let step_interval = positive_u64(controls[2], "NSAVC", source_label)?;
        let fixed_count = nonnegative_usize(controls[8], "NAMNF", source_label)?;
        let delta = endian.f32(record[40..44].try_into().expect("header slice"));
        if !delta.is_finite() || delta < 0.0 {
            return Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                source_label,
                "DCD DELTA must be finite and nonnegative",
            ));
        }
        let has_cell = controls[10] != 0;
        let charmm_version = controls[19];
        if charmm_version <= 0 {
            return Err(codec_context(
                TrajectoryIoErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                source_label,
                "only the common CHARMM/NAMD/OpenMM DCD profile is supported",
            ));
        }

        read_record(
            &mut reader,
            endian,
            &mut record,
            Some(0),
            limits,
            source_label,
            TrajectoryIoOperation::ReadHeader,
            None,
            false,
        )?;
        validate_title_record(&record, endian, source_label)?;
        read_record(
            &mut reader,
            endian,
            &mut record,
            Some(0),
            limits,
            source_label,
            TrajectoryIoOperation::ReadHeader,
            None,
            false,
        )?;
        if record.len() != 4 {
            return Err(header_error(
                source_label,
                "DCD atom-count record is not 4 bytes",
            ));
        }
        let file_atoms = nonnegative_usize(
            endian.i32(record[..4].try_into().expect("atom-count slice")),
            "NATOM",
            source_label,
        )?;
        validate_atom_count(file_atoms, limits, source_label)?;
        if file_atoms != atom_count {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
            )
            .with_source_label(source_label)
            .with_counts(atom_count as u64, file_atoms as u64));
        }
        if fixed_count > atom_count {
            return Err(header_error(
                source_label,
                "DCD fixed-atom count exceeds atom count",
            ));
        }
        let free_count = atom_count - fixed_count;
        let mut free_indices = Vec::new();
        let mut fixed_indices = Vec::new();
        let reserve_error =
            |detail: &str| resource_error(TrajectoryIoOperation::Open, source_label, None, detail);
        if fixed_count > 0 {
            read_record(
                &mut reader,
                endian,
                &mut record,
                Some(0),
                limits,
                source_label,
                TrajectoryIoOperation::ReadHeader,
                None,
                false,
            )?;
            let expected = checked_coordinate_bytes(free_count, source_label)?;
            if record.len() != expected {
                return Err(header_error(
                    source_label,
                    "DCD free-atom index record has the wrong size",
                ));
            }
            reserve_scratch(
                &mut free_indices,
                free_count,
                scratch_bytes(&record),
                limits.max_scratch_bytes(),
                reserve_error,
            )?;
            let mut seen = Vec::new();
            reserve_scratch(
                &mut seen,
                atom_count,
                scratch_bytes(&record)
                    .and_then(|bytes| bytes.checked_add(scratch_bytes(&free_indices)?)),
                limits.max_scratch_bytes(),
                reserve_error,
            )?;
            seen.resize(atom_count, false);
            for chunk in record.as_chunks::<4>().0 {
                let one_based = endian.i32(*chunk);
                if one_based <= 0 {
                    return Err(header_error(
                        source_label,
                        "DCD free-atom indices must be positive and one-based",
                    ));
                }
                let index = usize::try_from(one_based - 1).map_err(|_| {
                    header_error(source_label, "DCD free-atom index does not fit usize")
                })?;
                if index >= atom_count || seen[index] {
                    return Err(header_error(
                        source_label,
                        "DCD free-atom indices are out of range or duplicated",
                    ));
                }
                seen[index] = true;
                free_indices.push(index);
            }
            reserve_scratch(
                &mut fixed_indices,
                fixed_count,
                scratch_bytes(&record)
                    .and_then(|bytes| bytes.checked_add(scratch_bytes(&free_indices)?))
                    .and_then(|bytes| bytes.checked_add(scratch_bytes(&seen)?)),
                limits.max_scratch_bytes(),
                reserve_error,
            )?;
            fixed_indices.extend(
                seen.iter()
                    .enumerate()
                    .filter_map(|(index, free)| (!free).then_some(index)),
            );
        }
        let index_bytes = scratch_bytes(&free_indices)
            .and_then(|bytes| bytes.checked_add(scratch_bytes(&fixed_indices)?));
        let mut fixed_reference = Vec::new();
        reserve_scratch(
            &mut fixed_reference,
            fixed_count,
            index_bytes.and_then(|bytes| bytes.checked_add(scratch_bytes(&record)?)),
            limits.max_scratch_bytes(),
            reserve_error,
        )?;
        let mut positions = Vec::new();
        reserve_scratch(
            &mut positions,
            atom_count,
            index_bytes
                .and_then(|bytes| bytes.checked_add(scratch_bytes(&fixed_reference)?))
                .and_then(|bytes| bytes.checked_add(scratch_bytes(&record)?)),
            limits.max_scratch_bytes(),
            reserve_error,
        )?;
        positions.resize(atom_count, Point3::new(0.0, 0.0, 0.0));
        let data_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(TrajectoryFormat::Dcd),
                source_label,
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
            topology,
            options,
            header,
            free_indices,
            fixed_indices,
            fixed_reference,
            positions,
            record,
            frame_cursor: 0,
        })
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.topology.as_ref()
    }

    pub(crate) fn header(&self) -> &DcdHeader {
        &self.header
    }

    fn early_eof_error(&self, frame_start: u64) -> TrajectoryIoError {
        TrajectoryIoError::codec(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::ReadFrame,
            Some(TrajectoryFormat::Dcd),
        )
        .with_source_label(&self.options.source_label)
        .with_frame(self.frame_cursor)
        .with_byte_offset(frame_start)
        .with_counts(self.header.declared_frames as u64, self.frame_cursor as u64)
        .with_detail("DCD ended before its declared NSET count")
    }

    fn parse_next(
        &mut self,
        capture_fixed_reference: bool,
    ) -> Result<Option<DcdDecodedFrame>, TrajectoryIoError> {
        let other_scratch_bytes = scratch_bytes(&self.positions)
            .and_then(|bytes| bytes.checked_add(scratch_bytes(&self.free_indices)?))
            .and_then(|bytes| bytes.checked_add(scratch_bytes(&self.fixed_indices)?))
            .and_then(|bytes| bytes.checked_add(scratch_bytes(&self.fixed_reference)?));
        let frame_start = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Dcd),
                &self.options.source_label,
                error,
            )
        })?;
        if self.frame_cursor >= self.header.declared_frames {
            if probe_seekable_eof(
                &mut self.reader,
                TrajectoryIoOperation::ReadFrame,
                TrajectoryFormat::Dcd,
                &self.options.source_label,
            )? {
                return Ok(None);
            }
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Dcd),
            )
            .with_source_label(&self.options.source_label)
            .with_frame(self.frame_cursor)
            .with_byte_offset(frame_start)
            .with_counts(
                self.header.declared_frames as u64,
                (self.frame_cursor as u64).saturating_add(1),
            )
            .with_detail("DCD contains a frame beyond its declared NSET count"));
        }
        if self.frame_cursor >= self.options.limits.max_frames() {
            if probe_seekable_eof(
                &mut self.reader,
                TrajectoryIoOperation::ReadFrame,
                TrajectoryFormat::Dcd,
                &self.options.source_label,
            )? {
                return Err(self.early_eof_error(frame_start));
            }
            return Err(resource_error(
                TrajectoryIoOperation::ReadFrame,
                &self.options.source_label,
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
                other_scratch_bytes,
                &self.options.limits,
                &self.options.source_label,
                TrajectoryIoOperation::ReadFrame,
                Some(self.frame_cursor),
                true,
            )? {
                return Err(self.early_eof_error(frame_start));
            }
            if self.record.len() != CELL_BYTES {
                return Err(frame_error(
                    TrajectoryIoErrorKind::InvalidFrame,
                    &self.options.source_label,
                    self.frame_cursor,
                    "DCD unit-cell record must contain six f64 values",
                ));
            }
            cell = Some(decode_cell(
                &self.record,
                self.header.endian,
                &self.options.source_label,
                self.frame_cursor,
            )?);
        }

        let coordinate_indices = if self.frame_cursor == 0 || self.header.fixed_count == 0 {
            None
        } else {
            // A prior random read may have overwritten the shared position
            // scratch. Every partial frame restores its immutable fixed atoms.
            for (&index, &point) in self.fixed_indices.iter().zip(&self.fixed_reference) {
                self.positions[index] = point;
            }
            Some(self.free_indices.as_slice())
        };
        let coordinate_count = coordinate_indices.map_or(self.header.atom_count, <[usize]>::len);
        let coordinate_bytes =
            checked_coordinate_bytes(coordinate_count, &self.options.source_label)?;
        if coordinate_bytes as u64 > self.options.limits.max_frame_bytes() {
            return Err(resource_error(
                TrajectoryIoOperation::ReadFrame,
                &self.options.source_label,
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
                other_scratch_bytes,
                &self.options.limits,
                &self.options.source_label,
                TrajectoryIoOperation::ReadFrame,
                Some(self.frame_cursor),
                clean_frame_eof,
            )? {
                return Err(self.early_eof_error(frame_start));
            }
            if self.record.len() != coordinate_bytes {
                return Err(frame_error(
                    TrajectoryIoErrorKind::InvalidFrame,
                    &self.options.source_label,
                    self.frame_cursor,
                    "DCD coordinate record has the wrong size",
                ));
            }
            for (value_index, chunk) in self.record.as_chunks::<4>().0.iter().enumerate() {
                let value = f64::from(self.header.endian.f32(*chunk));
                if !value.is_finite() {
                    return Err(frame_error(
                        TrajectoryIoErrorKind::InvalidFrame,
                        &self.options.source_label,
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
                &self.options.source_label,
                error,
            )
        })?;
        let frame_bytes = frame_end.checked_sub(frame_start).ok_or_else(|| {
            frame_error(
                TrajectoryIoErrorKind::InvalidFrame,
                &self.options.source_label,
                self.frame_cursor,
                "DCD frame end precedes its start",
            )
        })?;
        if frame_bytes > self.options.limits.max_frame_bytes() {
            return Err(resource_error(
                TrajectoryIoOperation::ReadFrame,
                &self.options.source_label,
                Some(self.frame_cursor),
                "DCD frame exceeds the configured byte limit",
            ));
        }
        if capture_fixed_reference && self.frame_cursor == 0 && self.header.fixed_count > 0 {
            self.fixed_reference.clear();
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
                (self.frame_cursor as u64)
                    .checked_mul(self.header.step_interval)
                    .ok_or_else(|| {
                        frame_error(
                            TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
                            &self.options.source_label,
                            self.frame_cursor,
                            "DCD frame step multiplication overflows",
                        )
                    })?,
            )
            .ok_or_else(|| {
                frame_error(
                    TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
                    &self.options.source_label,
                    self.frame_cursor,
                    "DCD frame step addition overflows",
                )
            })?;
        let time = if let DcdTimePolicy::HeaderDelta { .. } = self.options.time_policy {
            let time = (step as f64) * f64::from(self.header.delta);
            if !time.is_finite() {
                return Err(frame_error(
                    TrajectoryIoErrorKind::InvalidFrame,
                    &self.options.source_label,
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
                &self.options.source_label,
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
    ) -> Result<(), TrajectoryIoError> {
        let mut data = FrameData::new(Quantity::new(positions, ANGSTROM)).with_step(decoded.step);
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

    /// Indexes a freshly opened reader, whose next frame is frame 0.
    fn index_frames(
        mut self,
        control: &IndexBuildControl<'_>,
    ) -> Result<IndexedDcdReader<R>, TrajectoryIoError> {
        let mut offsets = FrameOffsets::empty(self.header.data_start);
        control.checkpoint(
            0,
            self.header.data_start,
            TrajectoryFormat::Dcd,
            &self.options.source_label,
        )?;
        loop {
            if offsets.len() == self.header.declared_frames {
                if probe_seekable_eof(
                    &mut self.reader,
                    TrajectoryIoOperation::Index,
                    TrajectoryFormat::Dcd,
                    &self.options.source_label,
                )? {
                    break;
                }
                return Err(TrajectoryIoError::codec(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                )
                .with_source_label(&self.options.source_label)
                .with_counts(
                    self.header.declared_frames as u64,
                    (offsets.len() as u64).saturating_add(1),
                )
                .with_detail("DCD contains frames beyond its declared NSET count"));
            }
            if let Some(limit) = projected_index_limit(offsets.len(), &self.options.limits) {
                if probe_seekable_eof(
                    &mut self.reader,
                    TrajectoryIoOperation::Index,
                    TrajectoryFormat::Dcd,
                    &self.options.source_label,
                )? {
                    return Err(TrajectoryIoError::codec(
                        TrajectoryIoErrorKind::InconsistentMetadata,
                        TrajectoryIoOperation::Index,
                        Some(TrajectoryFormat::Dcd),
                    )
                    .with_source_label(&self.options.source_label)
                    .with_counts(self.header.declared_frames as u64, offsets.len() as u64)
                    .with_detail("DCD ended before its declared NSET count"));
                }
                return Err(resource_error(
                    TrajectoryIoOperation::Index,
                    &self.options.source_label,
                    Some(offsets.len()),
                    format!("DCD index {limit} exceeds the configured limit"),
                ));
            }
            let offset = self.reader.stream_position().map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                    &self.options.source_label,
                    error,
                )
            })?;
            if self
                .parse_next(true)
                .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))?
                .is_none()
            {
                return Err(TrajectoryIoError::codec(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                )
                .with_source_label(&self.options.source_label)
                .with_counts(self.header.declared_frames as u64, offsets.len() as u64)
                .with_detail("DCD ended before its declared NSET count"));
            }
            let end = self.reader.stream_position().map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                    &self.options.source_label,
                    error,
                )
            })?;
            offsets.push(
                offset,
                end,
                &self.options.limits,
                TrajectoryFormat::Dcd,
                &self.options.source_label,
            )?;
            control.checkpoint(
                offsets.len(),
                end,
                TrajectoryFormat::Dcd,
                &self.options.source_label,
            )?;
        }
        self.reader
            .seek(SeekFrom::Start(self.header.data_start))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                    &self.options.source_label,
                    error,
                )
            })?;
        self.frame_cursor = 0;
        Ok(IndexedDcdReader {
            inner: self,
            offsets,
        })
    }
}

impl<R: Read + Seek> TrajectoryReader for DcdReader<R> {
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
        let offset = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(TrajectoryFormat::Dcd),
                &self.options.source_label,
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
///
/// Random and sequential reads share one set of decode scratch. Fixed-atom
/// reference coordinates remain immutable during random reads.
pub struct IndexedDcdReader<R> {
    inner: DcdReader<R>,
    offsets: FrameOffsets,
}

impl<R> fmt::Debug for IndexedDcdReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedDcdReader")
            .field("reader", &self.inner)
            .field("frame_count", &self.offsets.len())
            .finish_non_exhaustive()
    }
}

impl<R: Read + Seek> IndexedDcdReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order
    /// and indexes it for random access.
    ///
    /// The stream's position at construction is frame 0. Every frame from there to
    /// the end is decoded and verified before the reader is returned, and both
    /// sequential and indexed reads then address frames from that start. Topology,
    /// limits, and diagnostic labels follow [`DcdReader::new`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`DcdReader::new`], the decoding error of the
    /// first invalid frame, [`TrajectoryIoErrorKind::InconsistentMetadata`]
    /// when the frame count differs from the header `NSET` count,
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when the index outgrows the
    /// limits, and [`TrajectoryIoErrorKind::Io`] when reading or seeking
    /// fails.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: DcdReadOptions,
    ) -> Result<Self, TrajectoryIoError> {
        Self::build(
            reader,
            topology.into(),
            options,
            &IndexBuildControl::default(),
        )
    }

    /// Indexes the stream like [`Self::new`], reporting to `control`.
    pub(crate) fn build(
        reader: R,
        topology: Arc<Topology>,
        options: DcdReadOptions,
        control: &IndexBuildControl<'_>,
    ) -> Result<Self, TrajectoryIoError> {
        DcdReader::new(reader, topology, options)?.index_frames(control)
    }

    /// Opens the stream with previously verified `offsets` instead of
    /// scanning it; random reads verify each frame against them. The offsets
    /// must count the header's `NSET` frames. With fixed atoms, frame 0 is
    /// decoded and verified here for their reference coordinates.
    pub(crate) fn from_offsets(
        reader: R,
        topology: Arc<Topology>,
        options: DcdReadOptions,
        offsets: FrameOffsets,
    ) -> Result<Self, TrajectoryIoError> {
        let mut inner = DcdReader::new(reader, topology, options)?;
        let label = inner.options.source_label.clone();
        if offsets.len() != inner.header.declared_frames {
            return Err(index_mismatch(
                TrajectoryIoOperation::Open,
                TrajectoryFormat::Dcd,
                &label,
                "DCD index frame count differs from the header NSET count",
            )
            .with_counts(inner.header.declared_frames as u64, offsets.len() as u64));
        }
        offsets.verify_start(inner.header.data_start, TrajectoryFormat::Dcd, &label)?;
        if inner.header.fixed_count > 0 && !offsets.is_empty() {
            let data_start = inner.header.data_start;
            inner
                .parse_next(true)
                .map_err(|error| frame_offset_context(error, 0, data_start))?
                .ok_or(TrajectoryIoError::frame_index_out_of_range(0))?;
            offsets.verify_stream_frame_end(0, &mut inner.reader, TrajectoryFormat::Dcd, &label)?;
            inner
                .reader
                .seek(SeekFrom::Start(data_start))
                .map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Open,
                        Some(TrajectoryFormat::Dcd),
                        &label,
                        error,
                    )
                })?;
            inner.frame_cursor = 0;
        }
        Ok(Self { inner, offsets })
    }

    pub(crate) fn offsets(&self) -> &FrameOffsets {
        &self.offsets
    }

    /// The topology every frame belongs to.
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

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        self.inner.read_next(destination)
    }
}

impl<R: Read + Seek> SeekableTrajectoryReader for IndexedDcdReader<R> {
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
        let (offset, _) = self
            .offsets
            .frame(index)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        let saved_offset = self.inner.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Index,
                Some(TrajectoryFormat::Dcd),
                &self.inner.options.source_label,
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
                    &self.inner.options.source_label,
                    error,
                )
            })?;
        self.inner.frame_cursor = index;
        let result = self
            .inner
            .parse_next(false)
            .map_err(|error| frame_offset_context(error, index, offset))
            .and_then(|decoded| decoded.ok_or(TrajectoryIoError::frame_index_out_of_range(index)))
            .and_then(|decoded| {
                self.offsets.verify_stream_frame_end(
                    index,
                    &mut self.inner.reader,
                    TrajectoryFormat::Dcd,
                    &self.inner.options.source_label,
                )?;
                Ok(decoded)
            });
        let restore = self
            .inner
            .reader
            .seek(SeekFrom::Start(saved_offset))
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::Index,
                    Some(TrajectoryFormat::Dcd),
                    &self.inner.options.source_label,
                    error,
                )
            });
        self.inner.frame_cursor = saved_cursor;
        let decoded = result?;
        restore?;
        self.inner
            .publish(&self.inner.positions, decoded, destination)
    }
}

/// Canonical DCD writer over one seekable stream.
pub struct DcdWriter<W> {
    writer: W,
    topology: Arc<Topology>,
    options: DcdWriteOptions,
    source_label: String,
    header_start: u64,
    frame_count: usize,
    /// Set while a frame's bytes are being written and left set if writing
    /// fails, so a stream that may end in a partial frame is never extended
    /// or reported as finished.
    partial: bool,
    finalized: bool,
    axis: Vec<f32>,
}

impl<W> fmt::Debug for DcdWriter<W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DcdWriter")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("source_label", &self.source_label)
            .field("frame_count", &self.frame_count)
            .field("partial", &self.partial)
            .field("finalized", &self.finalized)
            .finish_non_exhaustive()
    }
}

impl<W: Write + Seek> DcdWriter<W> {
    /// A writer of `topology` frames that writes the header at the stream's
    /// current position; `source_label` names the destination in errors.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the step
    /// sequence or `DELTA` does not fit the i32/f32 header or the time-policy
    /// unit is not a time, [`TrajectoryIoErrorKind::LimitExceeded`] when the
    /// atom count is zero or exceeds `i32::MAX` or scratch cannot be reserved,
    /// and [`TrajectoryIoErrorKind::Io`] when writing the header fails.
    pub fn new(
        mut writer: W,
        topology: impl Into<Arc<Topology>>,
        options: DcdWriteOptions,
        source_label: impl Into<String>,
    ) -> Result<Self, TrajectoryIoError> {
        let topology = topology.into();
        let source_label = source_label.into();
        if options.step_interval == 0
            || options.start_step > i32::MAX as u64
            || options.step_interval > i32::MAX as u64
            || !options.delta.is_finite()
            || options.delta < 0.0
        {
            return Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
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
                TrajectoryIoErrorKind::LimitExceeded,
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
            partial: false,
            finalized: false,
            axis,
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
    pub fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }

    /// Patches the header frame count, leaving the stream at its end; later
    /// calls do nothing. No frame can be written afterwards.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidFrame`] when no frame was
    /// written or an earlier write left a partial frame,
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when the frame count exceeds
    /// `i32::MAX`, and [`TrajectoryIoErrorKind::Io`] when seeking or writing
    /// fails.
    pub fn finalize(&mut self) -> Result<(), TrajectoryIoError> {
        if self.finalized {
            return Ok(());
        }
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::Finish,
            TrajectoryFormat::Dcd,
            &self.source_label,
        )?;
        require_nonempty_writer(self.frame_count, TrajectoryFormat::Dcd, &self.source_label)?;
        if self.frame_count > i32::MAX as usize {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
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
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Finish,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                "DCD header count offset overflows",
            )
        })?;
        self.partial = true;
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
        self.partial = false;
        self.finalized = true;
        Ok(())
    }

    /// Finalizes the frame count, flushes, and returns the nonempty DCD stream.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Self::finalize`] and
    /// [`TrajectoryIoErrorKind::Io`] when flushing fails.
    pub fn finish(mut self) -> Result<W, TrajectoryIoError> {
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

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Dcd,
            &self.source_label,
        )?;
        if self.finalized {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                "cannot write a DCD frame after finalization",
            ));
        }
        if !frame.topology().shares_layout(self.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        if frame.velocities().is_some() {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                TrajectoryField::Velocities,
            ));
        }
        if frame.forces().is_some() {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                TrajectoryField::Forces,
            ));
        }
        if frame.properties().atoms().has_data() {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                TrajectoryField::AtomProperties,
            ));
        }
        if frame.properties().bonds().has_data() {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                TrajectoryField::BondProperties,
            ));
        }
        if !frame.properties().owner().is_empty() {
            return Err(writer_field_error(
                &self.source_label,
                self.frame_count,
                TrajectoryField::Properties,
            ));
        }
        if self.frame_count >= i32::MAX as usize {
            return Err(writer_overflow(
                &self.source_label,
                "DCD frame count exceeds i32 header capacity",
            ));
        }
        let expected_step = self
            .options
            .start_step
            .checked_add(
                (self.frame_count as u64)
                    .checked_mul(self.options.step_interval)
                    .ok_or_else(|| writer_overflow(&self.source_label, "DCD step overflow"))?,
            )
            .ok_or_else(|| writer_overflow(&self.source_label, "DCD step overflow"))?;
        if frame.step() != Some(expected_step) {
            return Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                format!("DCD frame step must be exactly {expected_step}"),
            ));
        }
        match (self.options.time_policy, frame.time()) {
            (DcdTimePolicy::Absent, None) => {}
            (DcdTimePolicy::Absent, Some(_)) => {
                return Err(TrajectoryIoError::unsupported_field(TrajectoryField::Time))
            }
            (DcdTimePolicy::HeaderDelta { unit }, Some(time)) => {
                let actual = time.value_in(unit).map_err(|error| {
                    codec_context(
                        TrajectoryIoErrorKind::InconsistentMetadata,
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Dcd),
                        &self.source_label,
                        format!("DCD time unit is incompatible: {error}"),
                    )
                })?;
                // Readers derive time as step × DELTA, and DELTA is stored
                // as f32. A time therefore matches when it agrees to f32
                // precision; rounding a timestep such as 0.002 to f32
                // already moves it by up to half an f32 epsilon.
                let expected = (expected_step as f64) * f64::from(self.options.delta);
                if !actual.is_finite()
                    || (actual - expected).abs() > f64::from(f32::EPSILON) * expected.abs()
                {
                    return Err(codec_context(
                        TrajectoryIoErrorKind::InconsistentMetadata,
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Dcd),
                        &self.source_label,
                        "DCD frame time does not match the explicit DELTA policy",
                    ));
                }
            }
            (DcdTimePolicy::HeaderDelta { .. }, None) => {
                return Err(codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Dcd),
                    &self.source_label,
                    "DCD explicit time policy requires frame time",
                ))
            }
        }
        let cell = frame.cell().copied();
        if self.options.writes_cells() != cell.is_some() {
            return Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::WriteFrame,
                Some(TrajectoryFormat::Dcd),
                &self.source_label,
                "DCD cell presence must match the writer header policy",
            ));
        }
        let factor = CANONICAL_LENGTH_UNIT
            .conversion_factor_to(ANGSTROM)
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(TrajectoryFormat::Dcd),
                    &self.source_label,
                    format!("DCD position unit is incompatible: {error}"),
                )
            })?;
        let positions = frame.positions().values();
        for point in *positions.value() {
            for value in [point.x * factor, point.y * factor, point.z * factor] {
                if !value.is_finite() || !(value as f32).is_finite() {
                    return Err(codec_context(
                        TrajectoryIoErrorKind::InvalidFrame,
                        TrajectoryIoOperation::WriteFrame,
                        Some(TrajectoryFormat::Dcd),
                        &self.source_label,
                        "DCD coordinate cannot be represented as finite f32",
                    ));
                }
            }
        }
        // Everything that can reject the frame is checked above or here, so a
        // rejected frame leaves the stream untouched.
        let cell_payload = cell
            .map(|cell| -> Result<[u8; CELL_BYTES], TrajectoryIoError> {
                let values = encode_cell(cell, &self.source_label)?;
                let mut payload = [0_u8; CELL_BYTES];
                for (chunk, value) in payload.as_chunks_mut::<8>().0.iter_mut().zip(values) {
                    chunk.copy_from_slice(&self.options.endian.encode_f64(value));
                }
                Ok(payload)
            })
            .transpose()?;
        let byte_len = positions
            .value()
            .len()
            .checked_mul(4)
            .and_then(|len| u32::try_from(len).ok())
            .ok_or_else(|| {
                writer_overflow(&self.source_label, "DCD coordinate record exceeds u32")
            })?;
        self.partial = true;
        if let Some(payload) = cell_payload {
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
                self.axis.push(value as f32);
            }
            self.writer
                .write_all(&self.options.endian.encode_u32(byte_len))
                .and_then(|_| {
                    for value in &self.axis {
                        self.writer
                            .write_all(&self.options.endian.encode_f32(*value))?;
                    }
                    self.writer
                        .write_all(&self.options.endian.encode_u32(byte_len))
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
        self.partial = false;
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
        i32::from(options.writes_cells()),
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
    for (chunk, control) in header[4..].as_chunks_mut::<4>().0.iter_mut().zip(controls) {
        chunk.copy_from_slice(&options.endian.encode_i32(control));
    }
    header[40..44].copy_from_slice(&options.endian.encode_f32(options.delta));
    header
}

fn parse_controls(
    bytes: &[u8],
    endian: DcdEndian,
    source_label: &str,
) -> Result<[i32; 20], TrajectoryIoError> {
    if bytes.len() != 80 {
        return Err(header_error(
            source_label,
            "DCD control block is not 80 bytes",
        ));
    }
    let mut controls = [0_i32; 20];
    for (target, chunk) in controls.iter_mut().zip(bytes.as_chunks::<4>().0.iter()) {
        *target = endian.i32(*chunk);
    }
    Ok(controls)
}

fn validate_title_record(
    record: &[u8],
    endian: DcdEndian,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
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
    frame: usize,
) -> Result<PeriodicCell, TrajectoryIoError> {
    let mut values = [0_f64; 6];
    for (target, chunk) in values.iter_mut().zip(record.as_chunks::<8>().0.iter()) {
        *target = endian.f64(*chunk);
    }
    if values.iter().any(|value| !value.is_finite()) {
        return Err(frame_error(
            TrajectoryIoErrorKind::InvalidFrame,
            source_label,
            frame,
            "DCD cell contains a non-finite value",
        ));
    }
    let [a, gamma_raw, b, beta_raw, alpha_raw, c] = values;
    if a <= 0.0 || b <= 0.0 || c <= 0.0 {
        return Err(frame_error(
            TrajectoryIoErrorKind::InvalidFrame,
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
            TrajectoryIoErrorKind::InvalidFrame,
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
            TrajectoryIoErrorKind::InvalidFrame,
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
            TrajectoryIoErrorKind::InvalidFrame,
            source_label,
            frame,
            format!("invalid DCD cell: {error}"),
        )
    })
}

fn encode_cell(cell: PeriodicCell, source_label: &str) -> Result<[f64; 6], TrajectoryIoError> {
    if cell.periodic_axes() != [true; 3] {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Dcd),
            source_label,
            "DCD writer requires all three periodic axes",
        ));
    }
    let [a, b, c] = cell.vectors().value_in(ANGSTROM).map_err(|error| {
        codec_context(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Dcd),
            source_label,
            format!("DCD cell unit is incompatible: {error}"),
        )
    })?;
    if a.x <= 0.0 || a.y != 0.0 || a.z != 0.0 || b.y <= 0.0 || b.z != 0.0 || c.z <= 0.0 {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Dcd),
            source_label,
            "DCD cells require a along positive x, b in the xy plane with positive y, and c with positive z",
        ));
    }
    let length =
        |vector: Vector3| (vector.x * vector.x + vector.y * vector.y + vector.z * vector.z).sqrt();
    let dot =
        |left: Vector3, right: Vector3| left.x * right.x + left.y * right.y + left.z * right.z;
    let la = length(a);
    let lb = length(b);
    let lc = length(c);
    let cosine = |left: Vector3, right: Vector3, ll: f64, lr: f64| {
        (dot(left, right) / (ll * lr)).clamp(-1.0, 1.0)
    };
    let alpha = cosine(b, c, lb, lc);
    let beta = cosine(a, c, la, lc);
    let gamma = cosine(a, b, la, lb);
    let values = [la, gamma, lb, beta, alpha, lc];
    if values.iter().any(|value| !value.is_finite()) {
        return Err(codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Dcd),
            source_label,
            "DCD cell cannot be represented as finite lengths and angles",
        ));
    }
    let mut record = [0; CELL_BYTES];
    for (chunk, value) in record.as_chunks_mut::<8>().0.iter_mut().zip(values) {
        *chunk = value.to_le_bytes();
    }
    decode_cell(&record, DcdEndian::Little, source_label, 0).map_err(|error| {
        codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(TrajectoryFormat::Dcd),
            source_label,
            format!(
                "DCD cell becomes degenerate in its lengths/cosines representation: {}",
                validation_detail(&error)
            ),
        )
    })?;
    Ok(values)
}

#[allow(clippy::too_many_arguments)]
fn read_record<R: Read>(
    reader: &mut R,
    endian: DcdEndian,
    scratch: &mut Vec<u8>,
    other_scratch_bytes: Option<usize>,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    operation: TrajectoryIoOperation,
    frame: Option<usize>,
    clean_eof: bool,
) -> Result<bool, TrajectoryIoError> {
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
    if size as u64 > limits.max_record_bytes() || size > limits.max_scratch_bytes() {
        return Err(resource_error(
            operation,
            source_label,
            frame,
            "DCD record exceeds configured record or scratch limit",
        ));
    }
    reserve_scratch(
        scratch,
        size,
        other_scratch_bytes,
        limits.max_scratch_bytes(),
        |detail| resource_error(operation, source_label, frame, detail),
    )?;
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
        let mut context = TrajectoryIoError::codec(
            TrajectoryIoErrorKind::RecordMarkerMismatch,
            operation,
            Some(TrajectoryFormat::Dcd),
        )
        .with_source_label(source_label)
        .with_counts(size as u64, u64::from(endian.u32(trailing)))
        .with_detail("DCD leading and trailing record markers differ");
        if let Some(frame) = frame {
            context = context.with_frame(frame);
        }
        return Err(context);
    }
    Ok(true)
}

fn read_exact_required<R: Read>(
    reader: &mut R,
    bytes: &mut [u8],
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: Option<usize>,
    detail: &str,
) -> Result<(), TrajectoryIoError> {
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
) -> Result<(), TrajectoryIoError> {
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
) -> Result<(), TrajectoryIoError> {
    if count == 0 || count > limits.max_atoms() || count > i32::MAX as usize {
        return Err(resource_error(
            TrajectoryIoOperation::ReadHeader,
            source_label,
            None,
            "DCD atom count is zero or exceeds configured/profile limits",
        ));
    }
    Ok(())
}

fn checked_coordinate_bytes(count: usize, source_label: &str) -> Result<usize, TrajectoryIoError> {
    count.checked_mul(4).ok_or_else(|| {
        resource_error(
            TrajectoryIoOperation::ReadFrame,
            source_label,
            None,
            "DCD coordinate record size overflows",
        )
    })
}

fn nonnegative_u64(value: i32, field: &str, source_label: &str) -> Result<u64, TrajectoryIoError> {
    u64::try_from(value)
        .map_err(|_| header_error(source_label, format!("DCD {field} must be nonnegative")))
}

fn positive_u64(value: i32, field: &str, source_label: &str) -> Result<u64, TrajectoryIoError> {
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
) -> Result<usize, TrajectoryIoError> {
    usize::try_from(value)
        .map_err(|_| header_error(source_label, format!("DCD {field} must be nonnegative")))
}

fn header_error(source_label: &str, detail: impl Into<String>) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::InvalidHeader,
        TrajectoryIoOperation::ReadHeader,
        Some(TrajectoryFormat::Dcd),
        source_label,
        detail,
    )
}

fn record_error(
    source_label: &str,
    operation: TrajectoryIoOperation,
    frame: Option<usize>,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    let mut context = TrajectoryIoError::codec(
        TrajectoryIoErrorKind::TruncatedRecord,
        operation,
        Some(TrajectoryFormat::Dcd),
    )
    .with_source_label(source_label)
    .with_detail(detail);
    if let Some(frame) = frame {
        context = context.with_frame(frame);
    }
    context
}

fn frame_error(
    kind: TrajectoryIoErrorKind,
    source_label: &str,
    frame: usize,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        kind,
        TrajectoryIoOperation::ReadFrame,
        Some(TrajectoryFormat::Dcd),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_detail(detail)
}

fn resource_error(
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: Option<usize>,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    let mut context = TrajectoryIoError::codec(
        TrajectoryIoErrorKind::LimitExceeded,
        operation,
        Some(TrajectoryFormat::Dcd),
    )
    .with_source_label(source_label)
    .with_detail(detail);
    if let Some(frame) = frame {
        context = context.with_frame(frame);
    }
    context
}

fn writer_overflow(source_label: &str, detail: impl Into<String>) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::LimitExceeded,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Dcd),
        source_label,
        detail,
    )
}

fn writer_field_error(
    source_label: &str,
    frame: usize,
    field: TrajectoryField,
) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        TrajectoryIoErrorKind::UnsupportedField,
        TrajectoryIoOperation::WriteFrame,
        Some(TrajectoryFormat::Dcd),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_field(field)
    .with_detail(format!("DCD cannot preserve {field}"))
}
