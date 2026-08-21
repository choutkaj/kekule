//! Bounded pure-Rust file codecs for fixed-topology trajectories.
//!
//! The codecs implement this crate's reusable frame-buffer streaming
//! contracts and depend on [`kekule`] only for topology, structural state,
//! geometry, and units.
//!
//! # Supported profiles
//!
//! | Format | Reader | Writer |
//! |---|---|---|
//! | XYZ | strict constant-count multi-frame element/x/y/z text; configured length unit (angstrom by default) | deterministic strict text; optional frame state is rejected |
//! | DCD | common 32-bit-record `CORD` files in either byte order, common unit cells, fixed-atom reconstruction, and strict `NSET` | canonical all-atom `CORD` in either byte order with optional unit cells |
//! | TRR | GROMACS XDR frames with f32 or f64 position, box, velocity, and force blocks | one explicit f32 or f64 precision with per-frame optional blocks |
//! | XTC | GROMACS magic 1995/2023, signed nonnegative i32 counts/steps, small uncompressed and ordinary compressed coordinates | magic 1995/2023 through the private audited encoder adapter at explicit lossy precision |
//!
//! Every open requires a [`TrajectoryTopologyBinding`], which couples one shared
//! `Arc<Topology>` with caller-supplied atom-order evidence. Equal
//! atom count is never correspondence evidence. Native units are converted once
//! at the codec boundary: DCD and default XYZ use angstrom, while TRR/XTC use
//! GROMACS nanometre/picosecond conventions. XTC coordinate resolution is
//! nominally `1 / precision` nanometres.
//!
//! Sequential readers retain one file handle and avoid a whole-file scan.
//! Indexed readers retain one handle, fully verify every frame during an
//! O(file-size) index build, store bounded checked offsets with capped geometric
//! growth, and then decode one complete frame per random read. Decoding
//! validates into reusable private scratch and publishes transactionally into
//! the caller's [`FrameBuffer`]. Random reads restore all sequential reader
//! state before publication. Clean EOF is accepted only between frames,
//! including through a bounded probe at exact frame/index limits.
//!
//! Path writers stage a temporary sibling. Only consuming
//! [`FileTrajectoryWriter::finish`] flushes, synchronizes, finalizes format
//! metadata, and publishes a nonempty trajectory. Any failed frame write or an
//! empty finish prevents publication.
//!
//! # Limits and unsupported formats
//!
//! [`TrajectoryIoLimits`] bounds attacker-controlled atoms, frames, records,
//! scratch, index storage, text, comments, and detection before allocation or
//! seeking. Amber NetCDF, PDB, GRO/G96, Amber ASCII, LAMMPS dump, reactive
//! trajectories, and compressed wrappers are outside this initial profile and
//! return structured unsupported-format or unsupported-variant errors.
pub mod dcd;
mod detect;
pub mod trr;
pub mod xtc;
pub mod xyz;

use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::{
    atomic::{AtomicU64, Ordering},
    Arc,
};

use crate::{
    AtomOrderAssertion, AtomOrderAssertionKind, FrameBuffer, SeekableTrajectoryReader,
    TrajectoryCodecErrorContext, TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat,
    TrajectoryFrameView, TrajectoryIoErrorContext, TrajectoryIoOperation, TrajectoryReader,
    TrajectoryWriter,
};
use kekule::topology::Topology;
use kekule::units::Unit;

static TEMPORARY_FILE_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// Explicit selection or bounded automatic trajectory format detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum TrajectoryFormatHint {
    Auto,
    Explicit(TrajectoryFormat),
}

/// Limits applied before attacker-controlled allocation, scanning, or seeking.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryIoLimits {
    pub max_atoms: usize,
    pub max_frames: u64,
    pub max_frame_bytes: u64,
    pub max_record_bytes: u64,
    pub max_scratch_bytes: usize,
    pub max_index_entries: usize,
    pub max_index_bytes: usize,
    pub max_text_line_bytes: usize,
    pub max_comment_bytes: usize,
    pub max_detection_bytes: usize,
}

impl Default for TrajectoryIoLimits {
    fn default() -> Self {
        Self {
            max_atoms: 10_000_000,
            max_frames: 100_000_000,
            max_frame_bytes: 4 * 1024 * 1024 * 1024,
            max_record_bytes: 4 * 1024 * 1024 * 1024,
            max_scratch_bytes: usize::try_from(4_u64 * 1024 * 1024 * 1024).unwrap_or(usize::MAX),
            max_index_entries: 100_000_000,
            max_index_bytes: 800_000_000,
            max_text_line_bytes: 1_048_576,
            max_comment_bytes: 1_048_576,
            max_detection_bytes: 4096,
        }
    }
}

/// Exact topology and caller-supplied atom-order evidence for a topology-free file.
#[derive(Debug, Clone)]
pub struct TrajectoryTopologyBinding {
    topology: Arc<Topology>,
    atom_order: AtomOrderAssertion,
}

impl TrajectoryTopologyBinding {
    pub fn new(
        topology: Arc<Topology>,
        atom_order: AtomOrderAssertion,
    ) -> Result<Self, TrajectoryError> {
        if !atom_order.is_compatible(&topology) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        Ok(Self {
            topology,
            atom_order,
        })
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub(crate) fn topology_arc(&self) -> &Arc<Topology> {
        &self.topology
    }

    pub fn atom_order(&self) -> &AtomOrderAssertion {
        &self.atom_order
    }
}

/// Per-frame availability of one trajectory field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum FieldAvailability {
    Required,
    Optional,
    Absent,
}

/// Field-presence contract reported by a selected codec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrajectoryFieldAvailability {
    pub positions: FieldAvailability,
    pub cell: FieldAvailability,
    pub velocities: FieldAvailability,
    pub forces: FieldAvailability,
    pub time: FieldAvailability,
    pub step: FieldAvailability,
    pub properties: FieldAvailability,
}

/// Scalar representation used by native coordinate payloads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ScalarPrecision {
    DecimalText,
    Float32,
    Float64,
    Mixed,
}

/// Scientific coordinate encoding reported by a codec.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum CoordinateEncoding {
    Lossless {
        precision: ScalarPrecision,
    },
    Lossy {
        precision: ScalarPrecision,
        resolution: f64,
        unit: Unit,
    },
}

/// Random-access behavior available from an opened reader.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RandomAccessCapability {
    SequentialOnly,
    Indexed,
}

/// Capabilities and structural facts verified for an opened trajectory.
///
/// Indexed metadata describes the fully verified file. Sequential TRR metadata
/// initially reports the first frame's precision and changes to
/// [`ScalarPrecision::Mixed`] after a frame with the other scalar width is
/// successfully read.
#[derive(Debug, Clone, PartialEq)]
pub struct FileTrajectoryMetadata {
    format: TrajectoryFormat,
    atom_count: usize,
    declared_frame_count: Option<u64>,
    indexed_frame_count: Option<u64>,
    fields: TrajectoryFieldAvailability,
    coordinate_encoding: CoordinateEncoding,
    random_access: RandomAccessCapability,
    variant: Option<String>,
}

impl FileTrajectoryMetadata {
    pub const fn format(&self) -> TrajectoryFormat {
        self.format
    }

    pub const fn atom_count(&self) -> usize {
        self.atom_count
    }

    pub const fn declared_frame_count(&self) -> Option<u64> {
        self.declared_frame_count
    }

    pub const fn indexed_frame_count(&self) -> Option<u64> {
        self.indexed_frame_count
    }

    pub const fn fields(&self) -> TrajectoryFieldAvailability {
        self.fields
    }

    pub const fn coordinate_encoding(&self) -> CoordinateEncoding {
        self.coordinate_encoding
    }

    pub const fn random_access(&self) -> RandomAccessCapability {
        self.random_access
    }

    pub fn variant(&self) -> Option<&str> {
        self.variant.as_deref()
    }

    fn xyz(atom_count: usize, indexed_frame_count: Option<u64>) -> Self {
        Self {
            format: TrajectoryFormat::Xyz,
            atom_count,
            declared_frame_count: None,
            indexed_frame_count,
            fields: TrajectoryFieldAvailability {
                positions: FieldAvailability::Required,
                cell: FieldAvailability::Absent,
                velocities: FieldAvailability::Absent,
                forces: FieldAvailability::Absent,
                time: FieldAvailability::Absent,
                step: FieldAvailability::Absent,
                properties: FieldAvailability::Absent,
            },
            coordinate_encoding: CoordinateEncoding::Lossless {
                precision: ScalarPrecision::DecimalText,
            },
            random_access: if indexed_frame_count.is_some() {
                RandomAccessCapability::Indexed
            } else {
                RandomAccessCapability::SequentialOnly
            },
            variant: Some("strict multi-frame XYZ".into()),
        }
    }

    fn dcd(header: &dcd::DcdHeader, indexed_frame_count: Option<u64>) -> Self {
        Self {
            format: TrajectoryFormat::Dcd,
            atom_count: header.atom_count(),
            declared_frame_count: Some(header.declared_frames()),
            indexed_frame_count,
            fields: TrajectoryFieldAvailability {
                positions: FieldAvailability::Required,
                cell: if header.has_cell() {
                    FieldAvailability::Required
                } else {
                    FieldAvailability::Absent
                },
                velocities: FieldAvailability::Absent,
                forces: FieldAvailability::Absent,
                time: FieldAvailability::Optional,
                step: FieldAvailability::Required,
                properties: FieldAvailability::Absent,
            },
            coordinate_encoding: CoordinateEncoding::Lossless {
                precision: ScalarPrecision::Float32,
            },
            random_access: if indexed_frame_count.is_some() {
                RandomAccessCapability::Indexed
            } else {
                RandomAccessCapability::SequentialOnly
            },
            variant: Some(header.variant()),
        }
    }

    fn trr(
        header: &trr::TrrFrameHeader,
        indexed_frame_count: Option<u64>,
        precision_mixed: bool,
    ) -> Self {
        let precision = if precision_mixed {
            ScalarPrecision::Mixed
        } else {
            match header.precision() {
                trr::TrrScalarPrecision::Float32 => ScalarPrecision::Float32,
                trr::TrrScalarPrecision::Float64 => ScalarPrecision::Float64,
            }
        };
        Self {
            format: TrajectoryFormat::Trr,
            atom_count: header.atom_count(),
            declared_frame_count: None,
            indexed_frame_count,
            fields: TrajectoryFieldAvailability {
                positions: FieldAvailability::Required,
                cell: FieldAvailability::Optional,
                velocities: FieldAvailability::Optional,
                forces: FieldAvailability::Optional,
                time: FieldAvailability::Required,
                step: FieldAvailability::Required,
                properties: FieldAvailability::Optional,
            },
            coordinate_encoding: CoordinateEncoding::Lossless { precision },
            random_access: if indexed_frame_count.is_some() {
                RandomAccessCapability::Indexed
            } else {
                RandomAccessCapability::SequentialOnly
            },
            variant: Some(format!(
                "GROMACS TRR/XDR first-frame {:?}{}{}{}",
                header.precision(),
                if header.has_cell() { " cell" } else { "" },
                if header.has_velocities() {
                    " velocities"
                } else {
                    ""
                },
                if header.has_forces() { " forces" } else { "" },
            )),
        }
    }

    fn update_trr_precision(&mut self, header: &trr::TrrFrameHeader, precision_mixed: bool) {
        self.coordinate_encoding = CoordinateEncoding::Lossless {
            precision: if precision_mixed {
                ScalarPrecision::Mixed
            } else {
                match header.precision() {
                    trr::TrrScalarPrecision::Float32 => ScalarPrecision::Float32,
                    trr::TrrScalarPrecision::Float64 => ScalarPrecision::Float64,
                }
            },
        };
    }

    fn xtc(
        info: &xtc::XtcFrameInfo,
        cell_policy: xtc::XtcCellPolicy,
        indexed_frame_count: Option<u64>,
    ) -> Self {
        let coordinate_encoding = match info.precision() {
            Some(precision) => CoordinateEncoding::Lossy {
                precision: ScalarPrecision::Float32,
                resolution: 1.0 / f64::from(precision),
                unit: kekule::units::NANOMETER,
            },
            None => CoordinateEncoding::Lossless {
                precision: ScalarPrecision::Float32,
            },
        };
        Self {
            format: TrajectoryFormat::Xtc,
            atom_count: info.atom_count(),
            declared_frame_count: None,
            indexed_frame_count,
            fields: TrajectoryFieldAvailability {
                positions: FieldAvailability::Required,
                cell: match cell_policy {
                    xtc::XtcCellPolicy::RequirePeriodic => FieldAvailability::Required,
                    xtc::XtcCellPolicy::ZeroMatrixAsAbsent => FieldAvailability::Optional,
                },
                velocities: FieldAvailability::Absent,
                forces: FieldAvailability::Absent,
                time: FieldAvailability::Required,
                step: FieldAvailability::Required,
                properties: FieldAvailability::Absent,
            },
            coordinate_encoding,
            random_access: if indexed_frame_count.is_some() {
                RandomAccessCapability::Indexed
            } else {
                RandomAccessCapability::SequentialOnly
            },
            variant: Some(format!(
                "{:?} checked reader / audited molly 0.6.1 writer; first payload {} bytes",
                info.magic(),
                info.compressed_bytes()
            )),
        }
    }
}

/// Evidence used to select a trajectory format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum FormatDetectionEvidence {
    ExplicitHint,
    Signature,
    Extension,
    ExtensionSignatureAgreement,
    MissingExtension,
}

/// Bounded result of format detection without opening a codec.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryFormatDetection {
    format: TrajectoryFormat,
    evidence: Vec<FormatDetectionEvidence>,
}

impl TrajectoryFormatDetection {
    pub const fn format(&self) -> TrajectoryFormat {
        self.format
    }

    pub fn evidence(&self) -> &[FormatDetectionEvidence] {
        &self.evidence
    }
}

/// Detects a trajectory format from a bounded prefix and an optional filename.
///
/// The reader is restored to its original stream position before this function
/// returns. Automatic detection never trusts a known extension without a
/// conclusive signature.
pub fn detect_trajectory_format<R: io::Read + io::Seek>(
    reader: &mut R,
    source_name: impl AsRef<Path>,
    hint: TrajectoryFormatHint,
    limits: &TrajectoryIoLimits,
) -> Result<TrajectoryFormatDetection, TrajectoryError> {
    let source_name = source_name.as_ref();
    let source_label = source_name.display().to_string();
    let result = detect::select_format(reader, source_name, hint, limits, &source_label)?;
    Ok(TrajectoryFormatDetection {
        format: result.format,
        evidence: result.evidence,
    })
}

/// Non-fatal facts recorded while opening a trajectory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryOpenReport {
    selected_format: TrajectoryFormat,
    detection_evidence: Vec<FormatDetectionEvidence>,
    atom_order_evidence: AtomOrderAssertionKind,
    notes: Vec<String>,
}

impl TrajectoryOpenReport {
    pub const fn selected_format(&self) -> TrajectoryFormat {
        self.selected_format
    }

    pub fn detection_evidence(&self) -> &[FormatDetectionEvidence] {
        &self.detection_evidence
    }

    pub const fn atom_order_evidence(&self) -> AtomOrderAssertionKind {
        self.atom_order_evidence
    }

    pub fn notes(&self) -> &[String] {
        &self.notes
    }
}

/// Format-agnostic open configuration.
#[derive(Debug, Clone)]
pub struct TrajectoryOpenOptions {
    format_hint: TrajectoryFormatHint,
    limits: TrajectoryIoLimits,
    xyz: xyz::XyzReadOptions,
    dcd: dcd::DcdReadOptions,
    trr: trr::TrrReadOptions,
    xtc: xtc::XtcReadOptions,
}

impl Default for TrajectoryOpenOptions {
    fn default() -> Self {
        Self {
            format_hint: TrajectoryFormatHint::Auto,
            limits: TrajectoryIoLimits::default(),
            xyz: xyz::XyzReadOptions::default(),
            dcd: dcd::DcdReadOptions::default(),
            trr: trr::TrrReadOptions::default(),
            xtc: xtc::XtcReadOptions::default(),
        }
    }
}

impl TrajectoryOpenOptions {
    pub fn with_format_hint(mut self, format_hint: TrajectoryFormatHint) -> Self {
        self.format_hint = format_hint;
        self
    }

    pub fn with_limits(mut self, limits: TrajectoryIoLimits) -> Self {
        self.limits = limits;
        self
    }

    pub fn with_xyz_options(mut self, options: xyz::XyzReadOptions) -> Self {
        self.xyz = options;
        self
    }

    pub fn with_dcd_options(mut self, options: dcd::DcdReadOptions) -> Self {
        self.dcd = options;
        self
    }

    pub fn with_trr_options(mut self, options: trr::TrrReadOptions) -> Self {
        self.trr = options;
        self
    }

    pub fn with_xtc_options(mut self, options: xtc::XtcReadOptions) -> Self {
        self.xtc = options;
        self
    }
}

/// Existing-destination policy for a path writer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum OverwritePolicy {
    Forbid,
    Replace,
}

/// Format-agnostic path-writer configuration.
#[derive(Debug, Clone)]
pub struct TrajectoryWriteOptions {
    format: TrajectoryFormat,
    overwrite: OverwritePolicy,
    xyz: xyz::XyzWriteOptions,
    dcd: dcd::DcdWriteOptions,
    trr: trr::TrrWriteOptions,
    xtc: xtc::XtcWriteOptions,
}

impl TrajectoryWriteOptions {
    pub fn new(format: TrajectoryFormat) -> Self {
        Self {
            format,
            overwrite: OverwritePolicy::Forbid,
            xyz: xyz::XyzWriteOptions::default(),
            dcd: dcd::DcdWriteOptions::default(),
            trr: trr::TrrWriteOptions::default(),
            xtc: xtc::XtcWriteOptions::default(),
        }
    }

    pub fn with_overwrite_policy(mut self, overwrite: OverwritePolicy) -> Self {
        self.overwrite = overwrite;
        self
    }

    pub fn with_xyz_options(mut self, options: xyz::XyzWriteOptions) -> Self {
        self.xyz = options;
        self
    }

    pub fn with_dcd_options(mut self, options: dcd::DcdWriteOptions) -> Self {
        self.dcd = options;
        self
    }

    pub fn with_trr_options(mut self, options: trr::TrrWriteOptions) -> Self {
        self.trr = options;
        self
    }

    pub fn with_xtc_options(mut self, options: xtc::XtcWriteOptions) -> Self {
        self.xtc = options;
        self
    }
}

enum SequentialReaderInner {
    Xyz(xyz::XyzReader<BufReader<File>>),
    Dcd(dcd::DcdReader<BufReader<File>>),
    Trr(trr::TrrReader<BufReader<File>>),
    Xtc(xtc::XtcReader<BufReader<File>>),
}

/// Format-agnostic path-backed sequential reader retaining one file handle.
pub struct SequentialFileTrajectoryReader {
    inner: SequentialReaderInner,
    metadata: FileTrajectoryMetadata,
}

impl SequentialFileTrajectoryReader {
    /// Returns metadata verified through the most recent successful read.
    ///
    /// In particular, mixed-width TRR input is reported as mixed as soon as
    /// the second scalar width has been observed.
    pub fn metadata(&self) -> &FileTrajectoryMetadata {
        &self.metadata
    }
}

impl TrajectoryReader for SequentialFileTrajectoryReader {
    fn topology(&self) -> &Topology {
        match &self.inner {
            SequentialReaderInner::Xyz(reader) => reader.topology(),
            SequentialReaderInner::Dcd(reader) => reader.topology(),
            SequentialReaderInner::Trr(reader) => reader.topology(),
            SequentialReaderInner::Xtc(reader) => reader.topology(),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match &self.inner {
            SequentialReaderInner::Xyz(reader) => reader.shared_topology(),
            SequentialReaderInner::Dcd(reader) => reader.shared_topology(),
            SequentialReaderInner::Trr(reader) => reader.shared_topology(),
            SequentialReaderInner::Xtc(reader) => reader.shared_topology(),
        }
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        match &mut self.inner {
            SequentialReaderInner::Xyz(reader) => reader.read_next(destination),
            SequentialReaderInner::Dcd(reader) => reader.read_next(destination),
            SequentialReaderInner::Trr(reader) => {
                let read = reader.read_next(destination)?;
                self.metadata
                    .update_trr_precision(reader.first_header(), reader.precision_mixed());
                Ok(read)
            }
            SequentialReaderInner::Xtc(reader) => reader.read_next(destination),
        }
    }
}

enum IndexedReaderInner {
    Xyz(xyz::IndexedXyzReader<BufReader<File>>),
    Dcd(dcd::IndexedDcdReader<BufReader<File>>),
    Trr(trr::IndexedTrrReader<BufReader<File>>),
    Xtc(xtc::IndexedXtcReader<BufReader<File>>),
}

/// Format-agnostic path-backed indexed reader retaining one file handle.
pub struct IndexedFileTrajectoryReader {
    inner: IndexedReaderInner,
    metadata: FileTrajectoryMetadata,
}

impl IndexedFileTrajectoryReader {
    pub fn metadata(&self) -> &FileTrajectoryMetadata {
        &self.metadata
    }
}

impl TrajectoryReader for IndexedFileTrajectoryReader {
    fn topology(&self) -> &Topology {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.topology(),
            IndexedReaderInner::Dcd(reader) => reader.topology(),
            IndexedReaderInner::Trr(reader) => reader.topology(),
            IndexedReaderInner::Xtc(reader) => reader.topology(),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.shared_topology(),
            IndexedReaderInner::Dcd(reader) => reader.shared_topology(),
            IndexedReaderInner::Trr(reader) => reader.shared_topology(),
            IndexedReaderInner::Xtc(reader) => reader.shared_topology(),
        }
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        match &mut self.inner {
            IndexedReaderInner::Xyz(reader) => reader.read_next(destination),
            IndexedReaderInner::Dcd(reader) => reader.read_next(destination),
            IndexedReaderInner::Trr(reader) => reader.read_next(destination),
            IndexedReaderInner::Xtc(reader) => reader.read_next(destination),
        }
    }
}

impl SeekableTrajectoryReader for IndexedFileTrajectoryReader {
    fn frame_count(&self) -> Option<u64> {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.frame_count(),
            IndexedReaderInner::Dcd(reader) => reader.frame_count(),
            IndexedReaderInner::Trr(reader) => reader.frame_count(),
            IndexedReaderInner::Xtc(reader) => reader.frame_count(),
        }
    }

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        match &mut self.inner {
            IndexedReaderInner::Xyz(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Dcd(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Trr(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Xtc(reader) => reader.read_frame(index, destination),
        }
    }
}

/// Opens one fast sequential path-backed trajectory reader.
pub fn open_trajectory(
    path: impl AsRef<Path>,
    binding: TrajectoryTopologyBinding,
    options: TrajectoryOpenOptions,
) -> Result<(SequentialFileTrajectoryReader, TrajectoryOpenReport), TrajectoryError> {
    let path = path.as_ref();
    let label = path.display().to_string();
    let file = File::open(path)
        .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &label, error))?;
    let mut reader = BufReader::new(file);
    let detection = detect::select_format(
        &mut reader,
        path,
        options.format_hint,
        &options.limits,
        &label,
    )?;
    let order_kind = binding.atom_order.kind();
    let (inner, metadata, notes) = match detection.format {
        TrajectoryFormat::Xyz => {
            let mut reader =
                xyz::XyzReader::new(reader, binding, options.xyz, options.limits, label)?;
            reader.validate_first_frame()?;
            let atom_count = reader.topology().atom_count();
            (
                SequentialReaderInner::Xyz(reader),
                FileTrajectoryMetadata::xyz(atom_count, None),
                vec!["XYZ coordinates use the configured explicit/default length unit".into()],
            )
        }
        TrajectoryFormat::Dcd => {
            let reader = dcd::DcdReader::new(reader, binding, options.dcd, options.limits, label)?;
            let metadata = FileTrajectoryMetadata::dcd(reader.header(), None);
            (
                SequentialReaderInner::Dcd(reader),
                metadata,
                vec![
                    "DCD steps are preserved from ISTART/NSAVC; time follows the explicit policy"
                        .into(),
                ],
            )
        }
        TrajectoryFormat::Trr => {
            let reader = trr::TrrReader::new(reader, binding, options.trr, options.limits, label)?;
            let metadata = FileTrajectoryMetadata::trr(reader.first_header(), None, false);
            (
                SequentialReaderInner::Trr(reader),
                metadata,
                vec![
                    "TRR uses XDR big-endian scalars and explicit lambda preservation policy"
                        .into(),
                ],
            )
        }
        TrajectoryFormat::Xtc => {
            let xtc_options = options.xtc;
            let reader = xtc::XtcReader::new(reader, binding, xtc_options, options.limits, label)?;
            let metadata =
                FileTrajectoryMetadata::xtc(reader.first_info(), xtc_options.cell_policy(), None);
            (
                SequentialReaderInner::Xtc(reader),
                metadata,
                vec![
                    "XTC decoding uses a bounded checked reader; molly is confined to writing"
                        .into(),
                ],
            )
        }
    };
    let report = TrajectoryOpenReport {
        selected_format: detection.format,
        detection_evidence: detection.evidence,
        atom_order_evidence: order_kind,
        notes,
    };
    Ok((SequentialFileTrajectoryReader { inner, metadata }, report))
}

/// Opens one fully verified indexed path-backed trajectory reader.
pub fn open_indexed_trajectory(
    path: impl AsRef<Path>,
    binding: TrajectoryTopologyBinding,
    options: TrajectoryOpenOptions,
) -> Result<(IndexedFileTrajectoryReader, TrajectoryOpenReport), TrajectoryError> {
    let path = path.as_ref();
    let label = path.display().to_string();
    let file = File::open(path)
        .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &label, error))?;
    let mut reader = BufReader::new(file);
    let detection = detect::select_format(
        &mut reader,
        path,
        options.format_hint,
        &options.limits,
        &label,
    )?;
    let order_kind = binding.atom_order.kind();
    let (inner, metadata, notes) = match detection.format {
        TrajectoryFormat::Xyz => {
            let reader = xyz::XyzReader::new(reader, binding, options.xyz, options.limits, label)?;
            let reader = reader.to_indexed()?;
            let count = reader.frame_count().unwrap_or(0);
            let atom_count = reader.topology().atom_count();
            (
                IndexedReaderInner::Xyz(reader),
                FileTrajectoryMetadata::xyz(atom_count, Some(count)),
                vec!["XYZ index verified every complete frame".into()],
            )
        }
        TrajectoryFormat::Dcd => {
            let reader = dcd::DcdReader::new(reader, binding, options.dcd, options.limits, label)?;
            let reader = reader.to_indexed()?;
            let count = reader.frame_count().unwrap_or(0);
            let metadata = FileTrajectoryMetadata::dcd(reader.header(), Some(count));
            (
                IndexedReaderInner::Dcd(reader),
                metadata,
                vec!["DCD index verified record markers and the declared frame count".into()],
            )
        }
        TrajectoryFormat::Trr => {
            let reader = trr::TrrReader::new(reader, binding, options.trr, options.limits, label)?;
            let reader = reader.to_indexed()?;
            let count = reader.frame_count().unwrap_or(0);
            let metadata = FileTrajectoryMetadata::trr(
                reader.first_header(),
                Some(count),
                reader.precision_mixed(),
            );
            (
                IndexedReaderInner::Trr(reader),
                metadata,
                vec!["TRR index verified every XDR frame and payload block".into()],
            )
        }
        TrajectoryFormat::Xtc => {
            let xtc_options = options.xtc;
            let reader = xtc::XtcReader::new(reader, binding, xtc_options, options.limits, label)?;
            let reader = reader.to_indexed()?;
            let count = reader.frame_count().unwrap_or(0);
            let metadata = FileTrajectoryMetadata::xtc(
                reader.first_info(),
                xtc_options.cell_policy(),
                Some(count),
            );
            (
                IndexedReaderInner::Xtc(reader),
                metadata,
                vec!["XTC index fully decoded and verified every compressed frame".into()],
            )
        }
    };
    let report = TrajectoryOpenReport {
        selected_format: detection.format,
        detection_evidence: detection.evidence,
        atom_order_evidence: order_kind,
        notes,
    };
    Ok((IndexedFileTrajectoryReader { inner, metadata }, report))
}

enum FileWriterInner {
    Xyz(xyz::XyzWriter<BufWriter<File>>),
    Dcd(dcd::DcdWriter<BufWriter<File>>),
    Trr(trr::TrrWriter<BufWriter<File>>),
    Xtc(xtc::XtcWriter<BufWriter<File>>),
}

impl FileWriterInner {
    fn topology(&self) -> &Topology {
        match self {
            Self::Xyz(writer) => writer.topology(),
            Self::Dcd(writer) => writer.topology(),
            Self::Trr(writer) => writer.topology(),
            Self::Xtc(writer) => writer.topology(),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match self {
            Self::Xyz(writer) => writer.shared_topology(),
            Self::Dcd(writer) => writer.shared_topology(),
            Self::Trr(writer) => writer.shared_topology(),
            Self::Xtc(writer) => writer.shared_topology(),
        }
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        match self {
            Self::Xyz(writer) => writer.write_frame(frame),
            Self::Dcd(writer) => writer.write_frame(frame),
            Self::Trr(writer) => writer.write_frame(frame),
            Self::Xtc(writer) => writer.write_frame(frame),
        }
    }

    fn flush_and_sync(&mut self, label: &str) -> Result<(), TrajectoryError> {
        match self {
            Self::Xyz(writer) => {
                writer.validate_finish()?;
                writer.flush().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Xyz),
                        label,
                        error,
                    )
                })?;
                writer.writer().get_ref().sync_all().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Xyz),
                        label,
                        error,
                    )
                })
            }
            Self::Dcd(writer) => {
                writer.finalize()?;
                writer.flush().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Dcd),
                        label,
                        error,
                    )
                })?;
                writer.writer().get_ref().sync_all().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Dcd),
                        label,
                        error,
                    )
                })
            }
            Self::Trr(writer) => {
                writer.validate_finish()?;
                writer.flush().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Trr),
                        label,
                        error,
                    )
                })?;
                writer.writer().get_ref().sync_all().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Trr),
                        label,
                        error,
                    )
                })
            }
            Self::Xtc(writer) => {
                writer.validate_finish()?;
                writer.flush().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Xtc),
                        label,
                        error,
                    )
                })?;
                writer.writer().get_ref().sync_all().map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(TrajectoryFormat::Xtc),
                        label,
                        error,
                    )
                })
            }
        }
    }
}

/// Strict atomic path writer.
///
/// A nonempty file is published only by successful [`Self::finish`].
pub struct FileTrajectoryWriter {
    inner: Option<FileWriterInner>,
    format: TrajectoryFormat,
    destination: PathBuf,
    temporary: PathBuf,
    overwrite: OverwritePolicy,
    failed: bool,
    published: bool,
}

impl FileTrajectoryWriter {
    /// Flushes, synchronizes, and atomically publishes a nonempty trajectory.
    ///
    /// Finishing before any successful frame write returns a structured error
    /// and removes the unpublished temporary sibling.
    pub fn finish(mut self) -> Result<(), TrajectoryError> {
        let label = self.destination.display().to_string();
        if self.failed {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::Finish,
                Some(self.format()),
                &label,
                "cannot publish a trajectory after an earlier frame-write failure",
            ));
        }
        if let Some(inner) = &mut self.inner {
            inner.flush_and_sync(&label)?;
        }
        self.inner.take();
        match self.overwrite {
            OverwritePolicy::Forbid => {
                std::fs::hard_link(&self.temporary, &self.destination).map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(self.format()),
                        &label,
                        error,
                    )
                })?;
                self.published = true;
                let _ = std::fs::remove_file(&self.temporary);
            }
            OverwritePolicy::Replace => {
                std::fs::rename(&self.temporary, &self.destination).map_err(|error| {
                    io_context(
                        TrajectoryIoOperation::Finish,
                        Some(self.format()),
                        &label,
                        error,
                    )
                })?;
                self.published = true;
            }
        }
        Ok(())
    }

    pub const fn format(&self) -> TrajectoryFormat {
        self.format
    }
}

impl TrajectoryWriter for FileTrajectoryWriter {
    fn topology(&self) -> &Topology {
        match &self.inner {
            Some(inner) => inner.topology(),
            None => unreachable!("finished path writers are consumed"),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match &self.inner {
            Some(inner) => inner.shared_topology(),
            None => unreachable!("finished path writers are consumed"),
        }
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if self.failed {
            return Err(codec_context(
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::WriteFrame,
                Some(self.format),
                &self.destination.display().to_string(),
                "trajectory writer is poisoned by an earlier frame-write failure",
            ));
        }
        let result = match &mut self.inner {
            Some(inner) => inner.write_frame(frame),
            None => Err(codec_context(
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::WriteFrame,
                Some(self.format),
                &self.destination.display().to_string(),
                "trajectory writer has already finished",
            )),
        };
        if result.is_err() {
            self.failed = true;
        }
        result
    }
}

impl Drop for FileTrajectoryWriter {
    fn drop(&mut self) {
        if !self.published {
            self.inner.take();
            let _ = std::fs::remove_file(&self.temporary);
        }
    }
}

/// Creates a strict path writer backed by a temporary sibling file.
pub fn create_trajectory_writer(
    path: impl AsRef<Path>,
    topology: Arc<Topology>,
    options: TrajectoryWriteOptions,
) -> Result<FileTrajectoryWriter, TrajectoryError> {
    let destination = path.as_ref().to_path_buf();
    let label = destination.display().to_string();
    if options.overwrite == OverwritePolicy::Forbid && destination.exists() {
        return Err(io_context(
            TrajectoryIoOperation::Open,
            Some(options.format),
            &label,
            io::Error::new(io::ErrorKind::AlreadyExists, "destination already exists"),
        ));
    }
    #[cfg(windows)]
    if options.overwrite == OverwritePolicy::Replace && destination.exists() {
        return Err(codec_context(
            TrajectoryCodecErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::Open,
            Some(options.format),
            &label,
            "atomic replacement of an existing destination is unavailable on this platform",
        ));
    }
    let temporary = create_temporary_sibling(&destination, options.format)?;
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .map_err(|error| {
            io_context(
                TrajectoryIoOperation::Open,
                Some(options.format),
                &label,
                error,
            )
        })?;
    let inner_result = match options.format {
        TrajectoryFormat::Xyz => {
            xyz::XyzWriter::new(BufWriter::new(file), topology, options.xyz, label)
                .map(FileWriterInner::Xyz)
        }
        TrajectoryFormat::Dcd => {
            dcd::DcdWriter::new(BufWriter::new(file), topology, options.dcd, label)
                .map(FileWriterInner::Dcd)
        }
        TrajectoryFormat::Trr => {
            trr::TrrWriter::new(BufWriter::new(file), topology, options.trr, label)
                .map(FileWriterInner::Trr)
        }
        TrajectoryFormat::Xtc => {
            xtc::XtcWriter::new(BufWriter::new(file), topology, options.xtc, label)
                .map(FileWriterInner::Xtc)
        }
    };
    let inner = match inner_result {
        Ok(inner) => inner,
        Err(error) => {
            let _ = std::fs::remove_file(&temporary);
            return Err(error);
        }
    };
    Ok(FileTrajectoryWriter {
        inner: Some(inner),
        format: options.format,
        destination,
        temporary,
        overwrite: options.overwrite,
        failed: false,
        published: false,
    })
}

pub(crate) fn frame_offset_context(
    error: TrajectoryError,
    frame: u64,
    byte_offset: u64,
) -> TrajectoryError {
    match error {
        TrajectoryError::Io(context) => {
            let mut context = *context;
            if context.frame().is_none() {
                context = context.with_frame(frame);
            }
            if context.byte_offset().is_none() {
                context = context.with_byte_offset(byte_offset);
            }
            TrajectoryError::Io(Box::new(context))
        }
        TrajectoryError::Codec(context) => {
            let mut context = *context;
            if context.frame().is_none() {
                context = context.with_frame(frame);
            }
            if context.byte_offset().is_none() {
                context = context.with_byte_offset(byte_offset);
            }
            TrajectoryError::Codec(Box::new(context))
        }
        error => error,
    }
}

pub(crate) fn probe_seekable_eof<R: io::Read + io::Seek>(
    reader: &mut R,
    operation: TrajectoryIoOperation,
    format: TrajectoryFormat,
    source_label: &str,
) -> Result<bool, TrajectoryError> {
    let start = reader
        .stream_position()
        .map_err(|error| io_context(operation, Some(format), source_label, error))?;
    let mut byte = [0_u8; 1];
    let read = reader.read(&mut byte);
    let restore = reader
        .seek(io::SeekFrom::Start(start))
        .map_err(|error| io_context(operation, Some(format), source_label, error));
    match (read, restore) {
        (_, Err(error)) => Err(error),
        (Err(error), Ok(_)) => Err(io_context(operation, Some(format), source_label, error)),
        (Ok(0), Ok(_)) => Ok(true),
        (Ok(_), Ok(_)) => Ok(false),
    }
}

pub(crate) fn projected_index_limit(
    current_entries: usize,
    limits: &TrajectoryIoLimits,
) -> Option<&'static str> {
    if u64::try_from(current_entries).map_or(true, |entries| entries >= limits.max_frames) {
        return Some("frame count");
    }
    let Some(projected_entries) = current_entries.checked_add(1) else {
        return Some("entry count");
    };
    if projected_entries > limits.max_index_entries {
        return Some("entry count");
    }
    if projected_entries
        .checked_mul(std::mem::size_of::<u64>())
        .is_none_or(|bytes| bytes > limits.max_index_bytes)
    {
        return Some("byte count");
    }
    None
}

fn index_hard_capacity(limits: &TrajectoryIoLimits) -> usize {
    usize::try_from(limits.max_frames)
        .unwrap_or(usize::MAX)
        .min(limits.max_index_entries)
        .min(limits.max_index_bytes / std::mem::size_of::<u64>())
}

fn next_index_capacity(
    current_entries: usize,
    current_capacity: usize,
    hard_capacity: usize,
) -> Option<usize> {
    if current_entries < current_capacity {
        return None;
    }
    let minimum = current_entries.checked_add(1)?;
    if minimum > hard_capacity {
        return None;
    }
    let geometric = if current_capacity == 0 {
        8
    } else {
        current_capacity.saturating_mul(2)
    };
    Some(geometric.max(minimum).min(hard_capacity))
}

pub(crate) fn reserve_index_for_push(
    offsets: &mut Vec<u64>,
    limits: &TrajectoryIoLimits,
    format: TrajectoryFormat,
    source_label: &str,
    frame: u64,
) -> Result<(), TrajectoryError> {
    if offsets.len() < offsets.capacity() {
        return Ok(());
    }
    let hard_capacity = index_hard_capacity(limits);
    let Some(target_capacity) =
        next_index_capacity(offsets.len(), offsets.capacity(), hard_capacity)
    else {
        return Err(TrajectoryCodecErrorContext::new(
            TrajectoryCodecErrorKind::ResourceLimitExceeded,
            TrajectoryIoOperation::Index,
            Some(format),
        )
        .with_source_label(source_label)
        .with_frame(frame)
        .with_detail(format!(
            "{format} index reached its configured hard capacity"
        ))
        .into());
    };
    offsets
        .try_reserve_exact(target_capacity.saturating_sub(offsets.len()))
        .map_err(|_| {
            TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::Index,
                Some(format),
            )
            .with_source_label(source_label)
            .with_frame(frame)
            .with_detail(format!(
                "could not grow {format} index toward its bounded capacity"
            ))
            .into()
        })
}

pub(crate) fn require_nonempty_writer(
    frame_count: u64,
    format: TrajectoryFormat,
    source_label: &str,
) -> Result<(), TrajectoryError> {
    if frame_count == 0 {
        return Err(codec_context(
            TrajectoryCodecErrorKind::InvalidFrame,
            TrajectoryIoOperation::Finish,
            Some(format),
            source_label,
            format!("{format} production trajectories must contain at least one frame"),
        ));
    }
    Ok(())
}

fn create_temporary_sibling(
    destination: &Path,
    format: TrajectoryFormat,
) -> Result<PathBuf, TrajectoryError> {
    let parent = destination.parent().unwrap_or_else(|| Path::new("."));
    let name = destination
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("trajectory");
    for _ in 0..128 {
        let sequence = TEMPORARY_FILE_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let candidate = parent.join(format!(".{name}.{}.{}.tmp", std::process::id(), sequence));
        if !candidate.exists() {
            return Ok(candidate);
        }
    }
    Err(codec_context(
        TrajectoryCodecErrorKind::ResourceLimitExceeded,
        TrajectoryIoOperation::Open,
        Some(format),
        &destination.display().to_string(),
        "could not reserve a unique temporary sibling name",
    ))
}

pub(crate) fn io_context(
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: &str,
    error: io::Error,
) -> TrajectoryError {
    let mut context = TrajectoryIoErrorContext::new(operation, error.kind(), error.to_string())
        .with_source_label(source_label);
    if let Some(format) = format {
        context = context.with_format(format);
    }
    context.into()
}

pub(crate) fn codec_context(
    kind: TrajectoryCodecErrorKind,
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: &str,
    detail: impl Into<String>,
) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(kind, operation, format)
        .with_source_label(source_label)
        .with_detail(detail)
        .into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounded_index_reservation_grows_logarithmically() {
        let entry_count = 100_000;
        let limits = TrajectoryIoLimits {
            max_frames: entry_count as u64,
            max_index_entries: entry_count,
            max_index_bytes: entry_count * std::mem::size_of::<u64>(),
            ..TrajectoryIoLimits::default()
        };
        let mut offsets = Vec::new();
        let mut growth_events = 0;
        for frame in 0..entry_count {
            assert_eq!(projected_index_limit(offsets.len(), &limits), None);
            let previous_capacity = offsets.capacity();
            reserve_index_for_push(
                &mut offsets,
                &limits,
                TrajectoryFormat::Xyz,
                "capacity-test.xyz",
                frame as u64,
            )
            .unwrap();
            growth_events += usize::from(offsets.capacity() != previous_capacity);
            offsets.push(frame as u64);
        }
        assert_eq!(offsets.len(), entry_count);
        assert!(
            growth_events <= 16,
            "{growth_events} growth events are not logarithmic"
        );
    }

    #[test]
    fn index_hard_capacity_uses_the_smallest_configured_bound() {
        let limits = TrajectoryIoLimits {
            max_frames: 200,
            max_index_entries: 150,
            max_index_bytes: 125 * std::mem::size_of::<u64>(),
            ..TrajectoryIoLimits::default()
        };
        assert_eq!(index_hard_capacity(&limits), 125);
        assert_eq!(next_index_capacity(64, 64, 125), Some(125));
        assert_eq!(next_index_capacity(125, 125, 125), None);
    }
}
