use crate::io::write_policy::Carries;
use crate::units::Quantity;

use crate::traj::TrajectoryFormat;

use super::{dcd, trr, xtc, xyz, FormatDetectionEvidence};

/// Per-frame availability of one trajectory field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum FieldAvailability {
    /// Every frame carries the field.
    Required,
    /// Each frame may or may not carry the field.
    Optional,
    /// No frame carries the field.
    Absent,
}

/// Field-presence contract reported by a selected codec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrajectoryFieldAvailability {
    pub(crate) positions: FieldAvailability,
    pub(crate) cell: FieldAvailability,
    pub(crate) velocities: FieldAvailability,
    pub(crate) forces: FieldAvailability,
    pub(crate) time: FieldAvailability,
    pub(crate) step: FieldAvailability,
    pub(crate) properties: FieldAvailability,
}

impl TrajectoryFieldAvailability {
    /// Availability of atom positions.
    pub const fn positions(&self) -> FieldAvailability {
        self.positions
    }

    /// Availability of the periodic cell.
    pub const fn cell(&self) -> FieldAvailability {
        self.cell
    }

    /// Availability of atom velocities.
    pub const fn velocities(&self) -> FieldAvailability {
        self.velocities
    }

    /// Availability of atom forces.
    pub const fn forces(&self) -> FieldAvailability {
        self.forces
    }

    /// Availability of the frame time.
    pub const fn time(&self) -> FieldAvailability {
        self.time
    }

    /// Availability of the integration step.
    pub const fn step(&self) -> FieldAvailability {
        self.step
    }

    /// Availability of frame-level properties, such as the TRR lambda.
    pub const fn properties(&self) -> FieldAvailability {
        self.properties
    }

    /// The frame state a writer of a format with these fields carries: every
    /// field some frame may carry.
    pub(crate) const fn carries(self) -> Carries {
        const fn carried(field: FieldAvailability) -> bool {
            !matches!(field, FieldAvailability::Absent)
        }
        Carries {
            cell: carried(self.cell),
            time: carried(self.time),
            step: carried(self.step),
            velocities: carried(self.velocities),
            forces: carried(self.forces),
            properties: carried(self.properties),
            ..Carries::POSITIONS
        }
    }
}

/// Scalar representation used by native coordinate payloads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ScalarPrecision {
    /// Decimal text numbers.
    DecimalText,
    /// IEEE 754 single precision.
    Float32,
    /// IEEE 754 double precision.
    Float64,
    /// Frames of one file use different binary precisions.
    Mixed,
}

/// Scientific coordinate encoding reported by a codec.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum CoordinateEncoding {
    /// Coordinates are stored as plain scalars.
    Lossless {
        /// The stored scalar representation.
        precision: ScalarPrecision,
    },
    /// Coordinates are quantised, as in compressed XTC frames.
    Lossy {
        /// The scalar representation of the decoded coordinates.
        precision: ScalarPrecision,
        /// The nominal coordinate resolution, a length.
        resolution: Quantity<f64>,
    },
}

/// Random-access behavior available from an opened reader.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RandomAccessCapability {
    /// Frames are read in stream order only.
    SequentialOnly,
    /// A verified frame index supports random reads.
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
    declared_frame_count: Option<usize>,
    indexed_frame_count: Option<usize>,
    fields: TrajectoryFieldAvailability,
    coordinate_encoding: CoordinateEncoding,
    random_access: RandomAccessCapability,
    variant: Option<String>,
}

impl FileTrajectoryMetadata {
    /// The opened file's format.
    pub const fn format(&self) -> TrajectoryFormat {
        self.format
    }

    /// The atom count of every frame, equal to the topology's.
    pub const fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// The frame count the header declares, for formats with one (DCD);
    /// `None` otherwise.
    pub const fn declared_frame_count(&self) -> Option<usize> {
        self.declared_frame_count
    }

    /// The verified frame count of an indexed reader; `None` for a
    /// sequential reader.
    pub const fn indexed_frame_count(&self) -> Option<usize> {
        self.indexed_frame_count
    }

    /// Which frame fields the file carries.
    pub const fn fields(&self) -> TrajectoryFieldAvailability {
        self.fields
    }

    /// How coordinates are stored.
    pub const fn coordinate_encoding(&self) -> CoordinateEncoding {
        self.coordinate_encoding
    }

    /// Whether the reader supports random frame access.
    pub const fn random_access(&self) -> RandomAccessCapability {
        self.random_access
    }

    /// A human-readable description of the format variant, such as DCD byte
    /// order or XTC magic; diagnostic text, not a stable identifier.
    pub fn variant(&self) -> Option<&str> {
        self.variant.as_deref()
    }

    pub(super) fn xyz(atom_count: usize, indexed_frame_count: Option<usize>) -> Self {
        Self {
            format: TrajectoryFormat::Xyz,
            atom_count,
            declared_frame_count: None,
            indexed_frame_count,
            fields: xyz::FIELDS,
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

    pub(super) fn dcd(header: &dcd::DcdHeader, indexed_frame_count: Option<usize>) -> Self {
        Self {
            format: TrajectoryFormat::Dcd,
            atom_count: header.atom_count(),
            declared_frame_count: Some(header.declared_frames()),
            indexed_frame_count,
            fields: TrajectoryFieldAvailability {
                cell: if header.has_cell() {
                    FieldAvailability::Required
                } else {
                    FieldAvailability::Absent
                },
                ..dcd::FIELDS
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
        indexed_frame_count: Option<usize>,
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
            fields: trr::FIELDS,
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

    /// Records the frame counts of a refreshed indexed reader.
    pub(super) fn update_frame_counts(
        &mut self,
        indexed_frame_count: usize,
        declared_frame_count: Option<usize>,
    ) {
        self.indexed_frame_count = Some(indexed_frame_count);
        if self.declared_frame_count.is_some() {
            self.declared_frame_count = declared_frame_count;
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
        indexed_frame_count: Option<usize>,
    ) -> Self {
        let coordinate_encoding = match info.precision() {
            Some(precision) => CoordinateEncoding::Lossy {
                precision: ScalarPrecision::Float32,
                resolution: Quantity::new(1.0 / f64::from(precision), crate::units::NANOMETER),
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
                cell: match cell_policy {
                    xtc::XtcCellPolicy::RequirePeriodic => FieldAvailability::Required,
                    xtc::XtcCellPolicy::ZeroMatrixAsAbsent => FieldAvailability::Optional,
                },
                ..xtc::FIELDS
            },
            coordinate_encoding,
            random_access: if indexed_frame_count.is_some() {
                RandomAccessCapability::Indexed
            } else {
                RandomAccessCapability::SequentialOnly
            },
            variant: Some(format!(
                "{:?}; first payload {} bytes",
                info.magic(),
                info.compressed_bytes()
            )),
        }
    }
}

/// Non-fatal facts recorded while opening a trajectory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryOpenReport {
    pub(super) selected_format: TrajectoryFormat,
    pub(super) detection_evidence: Vec<FormatDetectionEvidence>,
    pub(super) notes: Vec<String>,
}

impl TrajectoryOpenReport {
    /// The format the reader was opened with.
    pub const fn selected_format(&self) -> TrajectoryFormat {
        self.selected_format
    }

    /// Why the format was selected; see
    /// [`TrajectoryFormatDetection::evidence`](super::TrajectoryFormatDetection::evidence).
    pub fn detection_evidence(&self) -> &[FormatDetectionEvidence] {
        &self.detection_evidence
    }

    /// Human-readable notes on format conventions applied while opening.
    pub fn notes(&self) -> &[String] {
        &self.notes
    }
}
