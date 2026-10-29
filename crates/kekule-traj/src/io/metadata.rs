use std::io;
use std::path::Path;

use kekule::diagnostics::Cancellation;
use kekule::units::Quantity;

use crate::{TrajectoryFormat, TrajectoryIoError};

use super::frame_index::IndexProgressCallback;
use super::{dcd, detect, trr, xtc, xyz, TrajectoryIndexProgress, TrajectoryIoLimits};

/// Explicit selection or bounded automatic trajectory format detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum TrajectoryFormatHint {
    /// Detect the format from a conclusive content signature, checked against
    /// the file extension when it names a known format.
    Auto,
    /// Use this format without inspecting the signature or extension.
    Explicit(TrajectoryFormat),
}

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

    pub(super) fn dcd(header: &dcd::DcdHeader, indexed_frame_count: Option<usize>) -> Self {
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
                resolution: Quantity::new(1.0 / f64::from(precision), kekule::units::NANOMETER),
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
                "{:?}; first payload {} bytes",
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
    /// The caller named the format with [`TrajectoryFormatHint::Explicit`].
    ExplicitHint,
    /// The content prefix carries a conclusive format signature.
    Signature,
    /// The file extension names the format.
    Extension,
    /// The signature and the extension select the same format.
    ExtensionSignatureAgreement,
    /// The path has no extension naming a supported format, so the signature
    /// alone selected it.
    MissingExtension,
}

/// Bounded result of format detection without opening a codec.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryFormatDetection {
    format: TrajectoryFormat,
    evidence: Vec<FormatDetectionEvidence>,
}

impl TrajectoryFormatDetection {
    /// The selected format.
    pub const fn format(&self) -> TrajectoryFormat {
        self.format
    }

    /// Why the format was selected, in the order the checks ran.
    pub fn evidence(&self) -> &[FormatDetectionEvidence] {
        &self.evidence
    }
}

/// Detects a trajectory format from a bounded prefix and an optional filename.
///
/// The reader is restored to its original stream position before this function
/// returns. Automatic detection never trusts a known extension without a
/// conclusive signature. `source_name` supplies the extension and the
/// diagnostic source label; nothing is opened.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::UnknownFormat`] when automatic detection
/// finds no conclusive signature, [`TrajectoryIoErrorKind::FormatMismatch`]
/// when the signature and a known extension disagree,
/// [`TrajectoryIoErrorKind::UnsupportedVariant`] for a gzip, xz or bzip2
/// wrapper, [`TrajectoryIoErrorKind::LimitExceeded`] when
/// [`TrajectoryIoLimits::max_detection_bytes`] is zero or exceeds
/// [`TrajectoryIoLimits::max_scratch_bytes`], and
/// [`TrajectoryIoErrorKind::Io`] when reading or seeking fails. The wrapper and
/// limit checks apply to explicit hints too.
///
/// [`TrajectoryIoErrorKind::UnknownFormat`]: crate::TrajectoryIoErrorKind::UnknownFormat
/// [`TrajectoryIoErrorKind::FormatMismatch`]: crate::TrajectoryIoErrorKind::FormatMismatch
/// [`TrajectoryIoErrorKind::UnsupportedVariant`]: crate::TrajectoryIoErrorKind::UnsupportedVariant
/// [`TrajectoryIoErrorKind::LimitExceeded`]: crate::TrajectoryIoErrorKind::LimitExceeded
/// [`TrajectoryIoErrorKind::Io`]: crate::TrajectoryIoErrorKind::Io
///
/// # Examples
///
/// ```
/// use std::io::Cursor;
///
/// use kekule_traj::io::{
///     detect_trajectory_format, FormatDetectionEvidence, TrajectoryFormatHint,
///     TrajectoryIoLimits,
/// };
/// use kekule_traj::TrajectoryFormat;
///
/// let mut content = Cursor::new(b"1\ncomment\nO 0.0 0.0 0.0\n".to_vec());
/// let detection = detect_trajectory_format(
///     &mut content,
///     "water.xyz",
///     TrajectoryFormatHint::Auto,
///     &TrajectoryIoLimits::default(),
/// )?;
/// assert_eq!(detection.format(), TrajectoryFormat::Xyz);
/// assert_eq!(
///     detection.evidence(),
///     [
///         FormatDetectionEvidence::Signature,
///         FormatDetectionEvidence::Extension,
///         FormatDetectionEvidence::ExtensionSignatureAgreement,
///     ]
/// );
/// assert_eq!(content.position(), 0);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn detect_trajectory_format<R: io::Read + io::Seek>(
    reader: &mut R,
    source_name: impl AsRef<Path>,
    hint: TrajectoryFormatHint,
    limits: &TrajectoryIoLimits,
) -> Result<TrajectoryFormatDetection, TrajectoryIoError> {
    let source_name = source_name.as_ref();
    let source_label = source_name.display().to_string();
    let result = detect::select_format(reader, source_name, hint, limits, &source_label, false)?;
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
    /// The format the reader was opened with.
    pub const fn selected_format(&self) -> TrajectoryFormat {
        self.selected_format
    }

    /// Why the format was selected; see
    /// [`TrajectoryFormatDetection::evidence`].
    pub fn detection_evidence(&self) -> &[FormatDetectionEvidence] {
        &self.detection_evidence
    }

    /// Human-readable notes on format conventions applied while opening.
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
    pub(super) index_cancellation: Option<Cancellation>,
    pub(super) index_progress: Option<IndexProgressCallback>,
    pub(super) growing: bool,
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
            index_cancellation: None,
            index_progress: None,
            growing: false,
        }
    }
}

impl TrajectoryOpenOptions {
    /// Sets [`Self::format_hint`].
    #[must_use]
    pub fn with_format_hint(mut self, format_hint: TrajectoryFormatHint) -> Self {
        self.format_hint = format_hint;
        self
    }

    /// Sets [`Self::limits`].
    #[must_use]
    pub fn with_limits(mut self, limits: TrajectoryIoLimits) -> Self {
        self.limits = limits;
        self
    }

    /// Sets [`Self::xyz_options`].
    #[must_use]
    pub fn with_xyz_options(mut self, options: xyz::XyzReadOptions) -> Self {
        self.xyz = options;
        self
    }

    /// Sets [`Self::dcd_options`].
    #[must_use]
    pub fn with_dcd_options(mut self, options: dcd::DcdReadOptions) -> Self {
        self.dcd = options;
        self
    }

    /// Sets [`Self::trr_options`].
    #[must_use]
    pub fn with_trr_options(mut self, options: trr::TrrReadOptions) -> Self {
        self.trr = options;
        self
    }

    /// Sets [`Self::xtc_options`].
    #[must_use]
    pub fn with_xtc_options(mut self, options: xtc::XtcReadOptions) -> Self {
        self.xtc = options;
        self
    }

    /// Sets [`Self::index_cancellation`].
    #[must_use]
    pub fn with_index_cancellation(mut self, cancellation: Cancellation) -> Self {
        self.index_cancellation = Some(cancellation);
        self
    }

    /// Sets the callback that [`open_indexed_trajectory`] calls with a
    /// [`TrajectoryIndexProgress`] before the first frame and after every
    /// verified frame while it builds its index, on the opening thread. It
    /// should be cheap, such as storing atomics a progress bar reads. Clones
    /// of these options share the callback. None by default.
    ///
    /// [`open_indexed_trajectory`]: super::open_indexed_trajectory
    #[must_use]
    pub fn with_index_progress(
        mut self,
        progress: impl Fn(TrajectoryIndexProgress) + Send + Sync + 'static,
    ) -> Self {
        self.index_progress = Some(IndexProgressCallback::new(progress));
        self
    }

    /// The cancellation that [`open_indexed_trajectory`] polls before the
    /// first frame and after every verified frame while it builds its index;
    /// once requested, opening fails with
    /// [`TrajectoryIoErrorKind::Cancelled`] and returns no reader. Other
    /// operations do not poll it. None by default.
    ///
    /// [`open_indexed_trajectory`]: super::open_indexed_trajectory
    /// [`TrajectoryIoErrorKind::Cancelled`]: crate::TrajectoryIoErrorKind::Cancelled
    pub const fn index_cancellation(&self) -> Option<&Cancellation> {
        self.index_cancellation.as_ref()
    }

    /// Sets [`Self::growing`].
    #[must_use]
    pub const fn with_growing(mut self, growing: bool) -> Self {
        self.growing = growing;
        self
    }

    /// Whether the file may still be growing: a writer, such as a running
    /// simulation in another process, is appending frames to it. False by
    /// default.
    ///
    /// Indexed readers of a growing file take a partially written final
    /// frame as not yet written rather than as corruption: they index the
    /// complete frames before it and report the bytes after them as
    /// [`pending`](super::IndexedFileTrajectoryReader::pending_byte_count),
    /// and [`refresh`](super::IndexedFileTrajectoryReader::refresh) indexes
    /// it once it is complete. Corrupt data is still an error wherever it
    /// lies. DCD readers count the frames in the file rather than check
    /// them against the header's `NSET`, which writers update only as they
    /// flush or when they finish. XYZ lines count only once their line
    /// break is written, so a growing XYZ file's last line must end in one.
    /// A file that holds no complete frame yet, or for DCD no complete
    /// header, fails to open with [`TrajectoryIoErrorKind::Incomplete`].
    ///
    /// Without it, a partial final frame is an error, as for a complete
    /// file. Sequential readers and [`read_trajectory`](super::read_trajectory)
    /// ignore it.
    ///
    /// [`TrajectoryIoErrorKind::Incomplete`]: crate::TrajectoryIoErrorKind::Incomplete
    pub const fn growing(&self) -> bool {
        self.growing
    }

    /// How the format is selected; [`TrajectoryFormatHint::Auto`] by default.
    pub const fn format_hint(&self) -> TrajectoryFormatHint {
        self.format_hint
    }

    /// Limits for path readers, which override those of the format options.
    pub const fn limits(&self) -> &TrajectoryIoLimits {
        &self.limits
    }

    /// Options for the XYZ reader; path readers override its limits
    /// and source label.
    pub const fn xyz_options(&self) -> &xyz::XyzReadOptions {
        &self.xyz
    }

    /// Options for the DCD reader; path readers override its limits
    /// and source label.
    pub const fn dcd_options(&self) -> &dcd::DcdReadOptions {
        &self.dcd
    }

    /// Options for the TRR reader; path readers override its limits
    /// and source label.
    pub const fn trr_options(&self) -> &trr::TrrReadOptions {
        &self.trr
    }

    /// Options for the XTC reader; path readers override its limits
    /// and source label.
    pub const fn xtc_options(&self) -> &xtc::XtcReadOptions {
        &self.xtc
    }
}
