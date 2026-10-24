//! Reusable frame buffers, streaming reader/writer contracts, and the
//! trajectory I/O error type.
//!
//! In-memory frames and trajectories live in [`kekule::structure`]
//! ([`kekule::structure::TrajectoryFrame`], [`kekule::structure::Trajectory`]).
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

use kekule::structure::{StructureError, TrajectoryError};

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
    /// The caller's [`kekule::diagnostics::Cancellation`] stopped the
    /// operation; [`kekule::diagnostics::Cancelled`] is the source.
    Cancelled,
    /// An encoded frame index is malformed or of an unsupported version.
    InvalidIndex,
    /// A frame index does not describe the file: it was built for another
    /// file, the file changed since, or a frame does not lie where the index
    /// places it.
    IndexMismatch,
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
        };
        formatter.write_str(description)
    }
}

/// Frame state a trajectory format cannot store.
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
    /// Properties of the frame collection.
    CollectionProperties,
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
    limit: Option<kekule::diagnostics::LimitExceeded>,
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

    #[must_use]
    /// Adds human-readable detail.
    pub fn with_detail(mut self, detail: impl Into<String>) -> Self {
        self.0.detail = Some(detail.into());
        self
    }

    #[must_use]
    /// Records the exceeded bound.
    pub fn with_limit(mut self, limit: kekule::diagnostics::LimitExceeded) -> Self {
        self.0.limit = Some(limit);
        self
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
    pub fn limit_exceeded(&self) -> Option<kekule::diagnostics::LimitExceeded> {
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

impl fmt::Display for TrajectoryIoError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let data = &self.0;
        match (data.kind, data.field) {
            (TrajectoryIoErrorKind::UnsupportedField, Some(field)) => {
                write!(formatter, "trajectory writer does not support {field}")?;
            }
            (kind, _) => write!(formatter, "{kind}")?,
        }
        if let Some(format) = data.format {
            write!(formatter, " for {format}")?;
        }
        if let Some(operation) = data.operation {
            write!(formatter, " while attempting to {operation}")?;
        }
        if let Some(source) = &data.source_label {
            write!(formatter, " at {source}")?;
        }
        if let Some(frame) = data.frame {
            write!(formatter, " in frame {frame}")?;
        }
        if let Some(offset) = data.byte_offset {
            write!(formatter, " at byte {offset}")?;
        }
        if let (Some(expected), Some(actual)) = (data.expected, data.actual) {
            write!(formatter, " (expected {expected}, actual {actual})")?;
        }
        if let Some(limit) = data.limit {
            write!(formatter, ": {limit}")?;
        } else if let Some(detail) = &data.detail {
            write!(formatter, ": {detail}")?;
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

impl From<kekule::structure::DenseArrayError> for TrajectoryIoError {
    fn from(error: kekule::structure::DenseArrayError) -> Self {
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
