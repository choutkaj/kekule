use std::io;
use std::path::Path;

use kekule::units::Unit;

use crate::{TrajectoryError, TrajectoryFormat};

use super::{dcd, detect, trr, xtc, xyz, TrajectoryIoLimits};

/// Explicit selection or bounded automatic trajectory format detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum TrajectoryFormatHint {
    Auto,
    Explicit(TrajectoryFormat),
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

    pub(super) fn xyz(atom_count: usize, indexed_frame_count: Option<u64>) -> Self {
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

    pub(super) fn dcd(header: &dcd::DcdHeader, indexed_frame_count: Option<u64>) -> Self {
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

    pub(super) fn trr(
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

    pub(super) fn update_trr_precision(
        &mut self,
        header: &trr::TrrFrameHeader,
        precision_mixed: bool,
    ) {
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

    pub(super) fn xtc(
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
    pub(super) selected_format: TrajectoryFormat,
    pub(super) detection_evidence: Vec<FormatDetectionEvidence>,
    pub(super) notes: Vec<String>,
}

impl TrajectoryOpenReport {
    pub const fn selected_format(&self) -> TrajectoryFormat {
        self.selected_format
    }

    pub fn detection_evidence(&self) -> &[FormatDetectionEvidence] {
        &self.detection_evidence
    }

    pub fn notes(&self) -> &[String] {
        &self.notes
    }
}

/// Format-agnostic configuration for opening readers or loading entire trajectories.
///
/// Path readers use these top-level limits and the path as their diagnostic
/// source label, overriding limits and source labels in format-specific options.
#[derive(Debug, Clone)]
pub struct TrajectoryOpenOptions {
    pub(super) format_hint: TrajectoryFormatHint,
    pub(super) limits: TrajectoryIoLimits,
    pub(super) xyz: xyz::XyzReadOptions,
    pub(super) dcd: dcd::DcdReadOptions,
    pub(super) trr: trr::TrrReadOptions,
    pub(super) xtc: xtc::XtcReadOptions,
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
