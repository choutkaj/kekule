//! Reusable frame buffers, streaming reader/writer contracts, and the
//! trajectory I/O error type.
//!
//! In-memory frames and trajectories live in [`crate::structure`]
//! ([`crate::structure::TrajectoryFrame`], [`crate::structure::Trajectory`]).
//! [`FrameBuffer`] provides reusable validated storage for streaming decoders.
//! The reader and writer traits publish complete frames transactionally.
//! [`TrajectoryIoError`] reports every codec and file failure.

mod buffer;
mod sample;
mod stream;

pub use buffer::{FrameBuffer, FrameData};
pub use sample::{
    read_ensemble, read_ensemble_seekable, EnsembleReadOptions, EnsembleReadProgress,
};
pub use stream::{
    MemoryTrajectoryWriter, SeekableTrajectoryReader, TrajectoryReader, TrajectoryWriter,
};

use std::{fmt, io};

use crate::structure::{StructureError, TrajectoryError};

/// Stable identity for a trajectory file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryFormat {
    /// Multi-frame extended XYZ text.
    Xyz,
    /// CHARMM/NAMD/LAMMPS binary DCD.
    Dcd,
    /// GROMACS compressed XTC.
    Xtc,
    /// GROMACS full-precision TRR.
    Trr,
}

impl TrajectoryFormat {
    /// The conventional upper-case format name, such as `"XYZ"`; also its
    /// `Display` form.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Xyz => "XYZ",
            Self::Dcd => "DCD",
            Self::Xtc => "XTC",
            Self::Trr => "TRR",
        }
    }
}

impl fmt::Display for TrajectoryFormat {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

/// File or stream operation active when trajectory I/O failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryIoOperation {
    /// Detecting the format from a path or content.
    Detect,
    /// Opening a file or stream.
    Open,
    /// Building the frame index for random access.
    Index,
    /// Reading the file header.
    ReadHeader,
    /// Reading a frame.
    ReadFrame,
    /// Writing a frame.
    WriteFrame,
    /// Finishing a writer: patching headers and flushing.
    Finish,
}

impl fmt::Display for TrajectoryIoOperation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::Detect => "detect",
            Self::Open => "open",
            Self::Index => "index",
            Self::ReadHeader => "read header",
            Self::ReadFrame => "read frame",
            Self::WriteFrame => "write frame",
            Self::Finish => "finish",
        };
        formatter.write_str(name)
    }
}

/// What kind of trajectory failure a [`TrajectoryIoError`] reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryIoErrorKind {
    /// Reading or writing the file or stream failed; the [`std::io::Error`]
    /// is the source.
    Io,
    /// The content matches no supported trajectory format.
    UnknownFormat,
    /// The content does not match the requested format.
    FormatMismatch,
    /// The file header is malformed.
    InvalidHeader,
    /// The file uses a variant of its format that is not supported.
    UnsupportedVariant,
    /// A record ends before its declared length.
    TruncatedRecord,
    /// A record declares an invalid length.
    InvalidRecordLength,
    /// A record's leading and trailing length markers differ.
    RecordMarkerMismatch,
    /// A frame is malformed or cannot be written as given.
    InvalidFrame,
    /// A frame's atom count differs from the trajectory's.
    InconsistentAtomCount,
    /// Header or frame metadata contradict each other.
    InconsistentMetadata,
    /// A precision setting is invalid.
    InvalidPrecision,
    /// A configured or format bound was exceeded;
    /// [`TrajectoryIoError::limit_exceeded`] names it when it is counted.
    LimitExceeded,
    /// The format cannot store a field of the frame;
    /// [`TrajectoryIoError::field`] names it.
    UnsupportedField,
    /// A step is negative or does not fit the format.
    NegativeOrUnrepresentableStep,
    /// Compressed coordinate data is corrupt.
    CorruptCompressedData,
    /// The object belongs to a different topology.
    TopologyMismatch,
    /// No frame has the requested index; [`TrajectoryIoError::frame`] is it.
    FrameIndexOutOfRange,
    /// A decoded or supplied frame is invalid for its topology; the
    /// [`StructureError`] is the source.
    InvalidStructure,
    /// Frames cannot form an in-memory trajectory; the [`TrajectoryError`]
    /// is the source.
    InvalidTrajectory,
    /// The caller's [`crate::diagnostics::Cancellation`] stopped the
    /// operation; [`crate::diagnostics::Cancelled`] is the source.
    Cancelled,
    /// An encoded frame index or file fingerprint is malformed or of an
    /// unsupported version.
    InvalidIndex,
    /// A frame index does not describe the file: it was built for another
    /// file, the file changed since, or a frame does not lie where the index
    /// places it.
    IndexMismatch,
    /// A growing file does not yet hold what opening it needs: it is empty
    /// or ends inside its header or, except for DCD, its first frame (see
    /// [`crate::traj::io::TrajectoryOpenOptions::with_growing`]). Open it again
    /// once its writer has written more.
    Incomplete,
}

impl fmt::Display for TrajectoryIoErrorKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let description = match self {
            Self::Io => "trajectory I/O failed",
            Self::UnknownFormat => "unknown trajectory format",
            Self::FormatMismatch => "trajectory format mismatch",
            Self::InvalidHeader => "invalid trajectory header",
            Self::UnsupportedVariant => "unsupported trajectory variant",
            Self::TruncatedRecord => "truncated trajectory record",
            Self::InvalidRecordLength => "invalid trajectory record length",
            Self::RecordMarkerMismatch => "trajectory record markers do not match",
            Self::InvalidFrame => "invalid trajectory frame",
            Self::InconsistentAtomCount => "inconsistent trajectory atom count",
            Self::InconsistentMetadata => "inconsistent trajectory metadata",
            Self::InvalidPrecision => "invalid trajectory precision",
            Self::LimitExceeded => "trajectory resource limit exceeded",
            Self::UnsupportedField => "unsupported trajectory field",
            Self::NegativeOrUnrepresentableStep => "negative or unrepresentable trajectory step",
            Self::CorruptCompressedData => "corrupt compressed trajectory data",
            Self::TopologyMismatch => "trajectory object belongs to a different topology",
            Self::FrameIndexOutOfRange => "trajectory frame index is out of range",
            Self::InvalidStructure => "invalid trajectory frame",
            Self::InvalidTrajectory => "invalid trajectory",
            Self::Cancelled => "trajectory operation cancelled",
            Self::InvalidIndex => "invalid trajectory frame index",
            Self::IndexMismatch => "trajectory frame index does not match the file",
            Self::Incomplete => "trajectory file holds no complete frame yet",
        };
        formatter.write_str(description)
    }
}

/// Frame, ensemble, or series state a trajectory format cannot store.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryField {
    /// The periodic cell.
    PeriodicCell,
    /// Atom velocities.
    Velocities,
    /// Atom forces.
    Forces,
    /// The frame time.
    Time,
    /// The integration step.
    Step,
    /// Per-atom properties.
    AtomProperties,
    /// Per-bond properties.
    BondProperties,
    /// Properties of the frame itself.
    Properties,
    /// Properties of the frame collection or ensemble.
    CollectionProperties,
    /// An ensemble sample's relative weight.
    Weight,
    /// The energy of an ensemble sample's state.
    Energy,
    /// The grouping of ensemble samples into states.
    State,
    /// The temperature an ensemble's populations describe.
    Temperature,
    /// Properties of an ensemble state.
    StateProperties,
    /// A series item's step on its series' axis.
    ItemStep,
    /// A series item's coordinate on its series' axis.
    ItemCoordinate,
}

impl fmt::Display for TrajectoryField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::PeriodicCell => "periodic cell",
            Self::Velocities => "velocities",
            Self::Forces => "forces",
            Self::Time => "time",
            Self::Step => "step",
            Self::AtomProperties => "atom properties",
            Self::BondProperties => "bond properties",
            Self::Properties => "frame properties",
            Self::CollectionProperties => "collection properties",
            Self::Weight => "ensemble sample weight",
            Self::Energy => "ensemble state energy",
            Self::State => "ensemble state",
            Self::Temperature => "ensemble temperature",
            Self::StateProperties => "ensemble state properties",
            Self::ItemStep => "series item step",
            Self::ItemCoordinate => "series item coordinate",
        })
    }
}

/// A trajectory read, write, or detection failure.
///
/// [`kind`](Self::kind) classifies it; the accessors give the operation,
/// format, source label, frame, byte offset and counts when known, and
/// [`source`](std::error::Error::source) the underlying [`std::io::Error`],
/// [`StructureError`] or [`TrajectoryError`]. `Display` gives the context
/// and never repeats the source.
///
/// Two errors are equal when every field agrees and their sources have the
/// same message.
#[derive(Debug, Clone)]
pub struct TrajectoryIoError(Box<TrajectoryIoErrorData>);

impl PartialEq for TrajectoryIoError {
    fn eq(&self, other: &Self) -> bool {
        let (left, right) = (&self.0, &other.0);
        left.kind == right.kind
            && left.operation == right.operation
            && left.format == right.format
            && left.source_label == right.source_label
            && left.frame == right.frame
            && left.byte_offset == right.byte_offset
            && left.expected == right.expected
            && left.actual == right.actual
            && left.limit == right.limit
            && left.field == right.field
            && left.detail == right.detail
            && left.source.as_ref().map(ToString::to_string)
                == right.source.as_ref().map(ToString::to_string)
    }
}

#[derive(Debug, Clone)]
struct TrajectoryIoErrorData {
    kind: TrajectoryIoErrorKind,
    operation: Option<TrajectoryIoOperation>,
    format: Option<TrajectoryFormat>,
    source_label: Option<String>,
    frame: Option<usize>,
    byte_offset: Option<u64>,
    expected: Option<u64>,
    actual: Option<u64>,
    /// Whether the detail or the source already writes the counts.
    counts_written: bool,
    limit: Option<crate::diagnostics::LimitExceeded>,
    field: Option<TrajectoryField>,
    detail: Option<String>,
    source: Option<std::sync::Arc<dyn std::error::Error + Send + Sync>>,
}

impl TrajectoryIoError {
    /// An error of `kind` without context; add it with the `with_*` methods.
    pub fn new(kind: TrajectoryIoErrorKind) -> Self {
        Self(Box::new(TrajectoryIoErrorData {
            kind,
            operation: None,
            format: None,
            source_label: None,
            frame: None,
            byte_offset: None,
            expected: None,
            actual: None,
            counts_written: false,
            limit: None,
            field: None,
            detail: None,
            source: None,
        }))
    }

    /// A codec failure of `kind` during `operation` on `format`.
    pub fn codec(
        kind: TrajectoryIoErrorKind,
        operation: TrajectoryIoOperation,
        format: Option<TrajectoryFormat>,
    ) -> Self {
        let mut error = Self::new(kind);
        error.0.operation = Some(operation);
        error.0.format = format;
        error
    }

    /// An I/O failure during `operation`, keeping `error` as the source.
    pub fn io(operation: TrajectoryIoOperation, error: io::Error) -> Self {
        Self::codec(TrajectoryIoErrorKind::Io, operation, None).with_source(error)
    }

    /// No frame has index `index`.
    pub fn frame_index_out_of_range(index: usize) -> Self {
        let mut error = Self::new(TrajectoryIoErrorKind::FrameIndexOutOfRange);
        error.0.frame = Some(index);
        error
    }

    /// The format cannot store `field`.
    pub fn unsupported_field(field: TrajectoryField) -> Self {
        Self::new(TrajectoryIoErrorKind::UnsupportedField).with_field(field)
    }

    #[must_use]
    /// Names the trajectory format.
    pub fn with_format(mut self, format: TrajectoryFormat) -> Self {
        self.0.format = Some(format);
        self
    }

    #[must_use]
    /// Names the file or stream, such as its path.
    pub fn with_source_label(mut self, source_label: impl Into<String>) -> Self {
        self.0.source_label = Some(source_label.into());
        self
    }

    #[must_use]
    /// Names the zero-based frame.
    pub fn with_frame(mut self, frame: usize) -> Self {
        self.0.frame = Some(frame);
        self
    }

    #[must_use]
    /// Names the byte offset in the file or stream.
    pub fn with_byte_offset(mut self, byte_offset: u64) -> Self {
        self.0.byte_offset = Some(byte_offset);
        self
    }

    #[must_use]
    /// Records the expected and actual count of a mismatch.
    pub fn with_counts(mut self, expected: u64, actual: u64) -> Self {
        self.0.expected = Some(expected);
        self.0.actual = Some(actual);
        self
    }

    /// Records the expected and actual count of a mismatch that the detail
    /// or the source already writes, so `Display` does not repeat them.
    #[must_use]
    pub(crate) fn with_written_counts(mut self, expected: u64, actual: u64) -> Self {
        self.0.counts_written = true;
        self.with_counts(expected, actual)
    }

    #[must_use]
    /// Adds human-readable detail.
    pub fn with_detail(mut self, detail: impl Into<String>) -> Self {
        self.0.detail = Some(detail.into());
        self
    }

    #[must_use]
    /// Records the exceeded bound.
    pub fn with_limit(mut self, limit: crate::diagnostics::LimitExceeded) -> Self {
        self.0.limit = Some(limit);
        self
    }

    #[must_use]
    /// Records the exceeded bound, with a detail naming what exceeded it:
    /// `{subject} exceeds its frame limit: 1,001 frames of 1,000.`
    pub(crate) fn exceeding(self, subject: &str, limit: crate::diagnostics::LimitExceeded) -> Self {
        self.with_detail(format!("{subject} exceeds its {}.", limit.describe()))
            .with_limit(limit)
    }

    #[must_use]
    /// Names the unsupported field.
    pub fn with_field(mut self, field: TrajectoryField) -> Self {
        self.0.field = Some(field);
        self
    }

    #[must_use]
    /// Records the underlying cause, returned by `source()`.
    pub fn with_source(mut self, source: impl std::error::Error + Send + Sync + 'static) -> Self {
        self.0.source = Some(std::sync::Arc::new(source));
        self
    }

    /// The error without its operation, format, source label, frame, and
    /// byte offset, for a reused validation failure whose context would
    /// misdescribe the operation that reuses it.
    pub(crate) fn without_context(mut self) -> Self {
        let data = &mut self.0;
        data.operation = None;
        data.format = None;
        data.source_label = None;
        data.frame = None;
        data.byte_offset = None;
        self
    }

    /// The stable classification of the failure.
    pub fn kind(&self) -> TrajectoryIoErrorKind {
        self.0.kind
    }

    /// The file or stream operation active when the failure occurred.
    pub fn operation(&self) -> Option<TrajectoryIoOperation> {
        self.0.operation
    }

    /// The trajectory format, when known.
    pub fn format(&self) -> Option<TrajectoryFormat> {
        self.0.format
    }

    /// The caller's label for the file or stream, such as its path.
    pub fn source_label(&self) -> Option<&str> {
        self.0.source_label.as_deref()
    }

    /// The zero-based frame the failure concerns.
    pub fn frame(&self) -> Option<usize> {
        self.0.frame
    }

    /// The byte offset in the file or stream the failure concerns.
    pub fn byte_offset(&self) -> Option<u64> {
        self.0.byte_offset
    }

    /// The expected count of a mismatch, such as an atom count.
    pub fn expected(&self) -> Option<u64> {
        self.0.expected
    }

    /// The actual count of a mismatch.
    pub fn actual(&self) -> Option<u64> {
        self.0.actual
    }

    /// The exceeded bound of a [`TrajectoryIoErrorKind::LimitExceeded`]
    /// failure, when the bound and the observed value are counted.
    pub fn limit_exceeded(&self) -> Option<crate::diagnostics::LimitExceeded> {
        self.0.limit
    }

    /// The field of a [`TrajectoryIoErrorKind::UnsupportedField`] failure.
    pub fn field(&self) -> Option<TrajectoryField> {
        self.0.field
    }

    /// Human-readable detail of the failure, without its context or source.
    pub fn detail(&self) -> Option<&str> {
        self.0.detail.as_deref()
    }

    /// The kind of the underlying [`std::io::Error`] of an
    /// [`TrajectoryIoErrorKind::Io`] failure.
    pub fn io_error_kind(&self) -> Option<io::ErrorKind> {
        self.0
            .source
            .as_deref()
            .and_then(|source| source.downcast_ref::<io::Error>())
            .map(io::Error::kind)
    }
}

/// Writes the operation that failed and its place, such as `Cannot read frame
/// 3 of the XTC trajectory "run.xtc" at byte 1200.`, then what went wrong:
/// the detail, the exceeded limit, or a sentence for the kind.
impl fmt::Display for TrajectoryIoError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        use TrajectoryIoErrorKind as Kind;
        let data = &self.0;
        let indexed = data.kind == Kind::FrameIndexOutOfRange;
        let mut frame_named = indexed;
        let head = data.operation.is_some();
        if let Some(operation) = data.operation {
            formatter.write_str(match operation {
                TrajectoryIoOperation::Detect => "Cannot detect the format of",
                TrajectoryIoOperation::Open => "Cannot open",
                TrajectoryIoOperation::Index => "Cannot index",
                TrajectoryIoOperation::ReadHeader => "Cannot read the header of",
                TrajectoryIoOperation::ReadFrame => "Cannot read",
                TrajectoryIoOperation::WriteFrame => "Cannot write",
                TrajectoryIoOperation::Finish => "Cannot finish",
            })?;
            let frame_object = matches!(
                operation,
                TrajectoryIoOperation::ReadFrame | TrajectoryIoOperation::WriteFrame
            );
            if let (true, false, Some(frame)) = (frame_object, frame_named, data.frame) {
                write!(formatter, " frame {frame} of")?;
                frame_named = true;
            }
            formatter.write_str(" the")?;
            if let Some(format) = data.format {
                write!(formatter, " {format}")?;
            }
            formatter.write_str(" trajectory")?;
            if let Some(label) = &data.source_label {
                write!(formatter, " \"{label}\"")?;
            }
            if let (false, Some(frame)) = (frame_named, data.frame) {
                write!(formatter, " in frame {frame}")?;
                frame_named = true;
            }
            if let Some(offset) = data.byte_offset {
                write!(formatter, " at byte {offset}")?;
            }
            formatter.write_str(".")?;
        }
        let separator = if head { " " } else { "" };
        let counts = (data.expected, data.actual);
        let mut counted = data.counts_written;
        let frame = data.frame.filter(|_| !frame_named);
        let body = match (&data.detail, data.limit, data.kind) {
            (Some(detail), _, _) => detail.clone(),
            (None, Some(limit), _) => limit.to_string(),
            (None, None, kind) => match kind {
                Kind::Io | Kind::Cancelled if head => return Ok(()),
                Kind::Io => "Reading or writing the trajectory failed.".to_owned(),
                Kind::Cancelled => "The trajectory operation stopped.".to_owned(),
                Kind::UnknownFormat => {
                    "The content matches no supported trajectory format.".to_owned()
                }
                Kind::FormatMismatch => "The content is not in the requested format.".to_owned(),
                Kind::InvalidHeader => "The header is malformed.".to_owned(),
                Kind::UnsupportedVariant => {
                    "This variant of the format is not supported.".to_owned()
                }
                Kind::TruncatedRecord => "A record ends before its declared length.".to_owned(),
                Kind::InvalidRecordLength => "A record declares an invalid length.".to_owned(),
                Kind::RecordMarkerMismatch => {
                    "A record's leading and trailing length markers differ.".to_owned()
                }
                Kind::InvalidFrame => "The frame is malformed.".to_owned(),
                Kind::InconsistentAtomCount => match counts {
                    (Some(expected), Some(actual)) => {
                        counted = true;
                        format!("The frame has {actual} atoms, not {expected}.")
                    }
                    _ => "The frame's atom count differs from the trajectory's.".to_owned(),
                },
                Kind::InconsistentMetadata => {
                    "The header and frame metadata contradict each other.".to_owned()
                }
                Kind::InvalidPrecision => "The precision is invalid.".to_owned(),
                Kind::LimitExceeded => "A format limit was exceeded.".to_owned(),
                Kind::UnsupportedField => match (data.format, data.field) {
                    (Some(format), Some(field)) => format!("{format} cannot store the {field}."),
                    (None, Some(field)) => format!("The format cannot store the {field}."),
                    _ => "The format cannot store a field of the frame.".to_owned(),
                },
                Kind::NegativeOrUnrepresentableStep => {
                    "A step is negative or does not fit the format.".to_owned()
                }
                Kind::CorruptCompressedData => "The compressed coordinates are corrupt.".to_owned(),
                Kind::TopologyMismatch => "The object belongs to a different topology.".to_owned(),
                Kind::FrameIndexOutOfRange => match data.frame {
                    Some(frame) => format!("Frame {frame} does not exist."),
                    None => "The frame does not exist.".to_owned(),
                },
                Kind::InvalidStructure => match frame {
                    Some(frame) => format!("Frame {frame} is invalid for its topology."),
                    None => "The frame is invalid for its topology.".to_owned(),
                },
                Kind::InvalidTrajectory => "The frames cannot form a trajectory.".to_owned(),
                Kind::InvalidIndex => {
                    "The frame index is malformed or of an unsupported version.".to_owned()
                }
                Kind::IndexMismatch => "The frame index does not match the file.".to_owned(),
                Kind::Incomplete => "The file holds no complete frame yet.".to_owned(),
            },
        };
        write!(formatter, "{separator}{body}")?;
        if let (false, Some(expected), Some(actual)) = (counted, counts.0, counts.1) {
            write!(formatter, " Expected {expected}, found {actual}.")?;
        }
        Ok(())
    }
}

impl std::error::Error for TrajectoryIoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.0
            .source
            .as_deref()
            .map(|source| source as &(dyn std::error::Error + 'static))
    }
}

impl From<StructureError> for TrajectoryIoError {
    fn from(error: StructureError) -> Self {
        Self::new(TrajectoryIoErrorKind::InvalidStructure).with_source(error)
    }
}

impl From<crate::structure::DenseArrayError> for TrajectoryIoError {
    fn from(error: crate::structure::DenseArrayError) -> Self {
        StructureError::from(error).into()
    }
}

impl From<TrajectoryError> for TrajectoryIoError {
    fn from(error: TrajectoryError) -> Self {
        match error {
            TrajectoryError::IndexOutOfRange { index, .. } => Self::frame_index_out_of_range(index),
            TrajectoryError::Structure(error) => error.into(),
            error => Self::new(TrajectoryIoErrorKind::InvalidTrajectory).with_source(error),
        }
    }
}
