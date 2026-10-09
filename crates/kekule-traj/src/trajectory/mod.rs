//! Reusable frame buffers and streaming reader/writer contracts.
//!
//! In-memory frames and trajectories live in [`kekule::structure`]
//! ([`kekule::structure::TrajectoryFrame`], [`kekule::structure::Trajectory`]).
//! [`FrameBuffer`] provides reusable validated storage for streaming decoders.
//! The reader and writer traits publish complete frames transactionally.

mod buffer;
mod stream;

pub use buffer::{FrameBuffer, FrameBufferData};
pub use stream::{
    validate_atom_order, CoordinateFrameReader, MemoryTrajectoryReader, MemoryTrajectoryWriter,
    SeekableTrajectoryReader, TrajectoryReader, TrajectoryWriter,
};

use std::{fmt, io};

use kekule::structure::{ConformationError, RealizationError};

/// Stable identity for a trajectory file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryFormat {
    Xyz,
    Dcd,
    Xtc,
    Trr,
}

impl TrajectoryFormat {
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
    Detect,
    Open,
    Index,
    ReadHeader,
    ReadFrame,
    WriteHeader,
    WriteFrame,
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
            Self::WriteHeader => "write header",
            Self::WriteFrame => "write frame",
            Self::Finish => "finish",
        };
        formatter.write_str(name)
    }
}

/// Typed classification for malformed, unsupported, or unsafe codec input.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryCodecErrorKind {
    UnknownFormat,
    FormatMismatch,
    InvalidHeader,
    UnsupportedVariant,
    TruncatedRecord,
    InvalidRecordLength,
    RecordMarkerMismatch,
    InvalidFrame,
    InconsistentAtomCount,
    InconsistentMetadata,
    InvalidPrecision,
    ResourceLimitExceeded,
    UnsupportedField,
    NegativeOrUnrepresentableStep,
    CorruptCompressedData,
}

impl fmt::Display for TrajectoryCodecErrorKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let description = match self {
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
            Self::ResourceLimitExceeded => "trajectory resource limit exceeded",
            Self::UnsupportedField => "unsupported trajectory field",
            Self::NegativeOrUnrepresentableStep => "negative or unrepresentable trajectory step",
            Self::CorruptCompressedData => "corrupt compressed trajectory data",
        };
        formatter.write_str(description)
    }
}

/// Cloneable typed context for an underlying file or stream error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryIoErrorContext {
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: Option<String>,
    frame: Option<u64>,
    byte_offset: Option<u64>,
    error_kind: io::ErrorKind,
    message: String,
}

impl TrajectoryIoErrorContext {
    pub fn new(
        operation: TrajectoryIoOperation,
        error_kind: io::ErrorKind,
        message: impl Into<String>,
    ) -> Self {
        Self {
            operation,
            format: None,
            source_label: None,
            frame: None,
            byte_offset: None,
            error_kind,
            message: message.into(),
        }
    }

    pub fn with_format(mut self, format: TrajectoryFormat) -> Self {
        self.format = Some(format);
        self
    }

    pub fn with_source_label(mut self, source_label: impl Into<String>) -> Self {
        self.source_label = Some(source_label.into());
        self
    }

    pub const fn with_frame(mut self, frame: u64) -> Self {
        self.frame = Some(frame);
        self
    }

    pub const fn with_byte_offset(mut self, byte_offset: u64) -> Self {
        self.byte_offset = Some(byte_offset);
        self
    }

    pub const fn operation(&self) -> TrajectoryIoOperation {
        self.operation
    }

    pub const fn format(&self) -> Option<TrajectoryFormat> {
        self.format
    }

    pub fn source_label(&self) -> Option<&str> {
        self.source_label.as_deref()
    }

    pub const fn frame(&self) -> Option<u64> {
        self.frame
    }

    pub const fn byte_offset(&self) -> Option<u64> {
        self.byte_offset
    }

    pub const fn error_kind(&self) -> io::ErrorKind {
        self.error_kind
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

/// Typed context for a codec validation or capability error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryCodecErrorContext {
    kind: TrajectoryCodecErrorKind,
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: Option<String>,
    frame: Option<u64>,
    byte_offset: Option<u64>,
    expected: Option<u64>,
    actual: Option<u64>,
    detail: Option<String>,
}

impl TrajectoryCodecErrorContext {
    pub const fn new(
        kind: TrajectoryCodecErrorKind,
        operation: TrajectoryIoOperation,
        format: Option<TrajectoryFormat>,
    ) -> Self {
        Self {
            kind,
            operation,
            format,
            source_label: None,
            frame: None,
            byte_offset: None,
            expected: None,
            actual: None,
            detail: None,
        }
    }

    pub fn with_source_label(mut self, source_label: impl Into<String>) -> Self {
        self.source_label = Some(source_label.into());
        self
    }

    pub const fn with_frame(mut self, frame: u64) -> Self {
        self.frame = Some(frame);
        self
    }

    pub const fn with_byte_offset(mut self, byte_offset: u64) -> Self {
        self.byte_offset = Some(byte_offset);
        self
    }

    pub const fn with_counts(mut self, expected: u64, actual: u64) -> Self {
        self.expected = Some(expected);
        self.actual = Some(actual);
        self
    }

    pub fn with_detail(mut self, detail: impl Into<String>) -> Self {
        self.detail = Some(detail.into());
        self
    }

    pub const fn kind(&self) -> TrajectoryCodecErrorKind {
        self.kind
    }

    pub const fn operation(&self) -> TrajectoryIoOperation {
        self.operation
    }

    pub const fn format(&self) -> Option<TrajectoryFormat> {
        self.format
    }

    pub fn source_label(&self) -> Option<&str> {
        self.source_label.as_deref()
    }

    pub const fn frame(&self) -> Option<u64> {
        self.frame
    }

    pub const fn byte_offset(&self) -> Option<u64> {
        self.byte_offset
    }

    pub const fn expected(&self) -> Option<u64> {
        self.expected
    }

    pub const fn actual(&self) -> Option<u64> {
        self.actual
    }

    pub fn detail(&self) -> Option<&str> {
        self.detail.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TrajectoryError {
    TopologyMismatch,
    AtomOrderMismatch,
    FrameIndexOutOfRange(u64),
    UnsupportedRandomAccess,
    MissingRequiredTopology,
    UnsupportedField(&'static str),
    /// Decoded or supplied frame state is invalid for its topology.
    Conformation(Box<ConformationError>),
    /// Frames could not be collected into an in-memory trajectory.
    Realization(Box<RealizationError>),
    Io(Box<TrajectoryIoErrorContext>),
    Codec(Box<TrajectoryCodecErrorContext>),
}

impl fmt::Display for TrajectoryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("trajectory object belongs to a different topology")
            }
            Self::AtomOrderMismatch => {
                formatter.write_str("coordinate-source atom order does not match topology order")
            }
            Self::FrameIndexOutOfRange(index) => {
                write!(formatter, "trajectory frame index {index} is out of range")
            }
            Self::UnsupportedRandomAccess => {
                formatter.write_str("trajectory source does not support random access")
            }
            Self::MissingRequiredTopology => {
                formatter.write_str("coordinate-only trajectory source requires a topology")
            }
            Self::UnsupportedField(field) => {
                write!(formatter, "trajectory writer does not support {field}")
            }
            Self::Conformation(error) => write!(formatter, "invalid trajectory frame: {error}"),
            Self::Realization(error) => write!(formatter, "invalid trajectory: {error}"),
            Self::Io(context) => {
                write!(formatter, "trajectory {} I/O failed", context.operation)?;
                if let Some(format) = context.format {
                    write!(formatter, " for {format}")?;
                }
                if let Some(source) = &context.source_label {
                    write!(formatter, " at {source}")?;
                }
                if let Some(frame) = context.frame {
                    write!(formatter, " in frame {frame}")?;
                }
                if let Some(offset) = context.byte_offset {
                    write!(formatter, " at byte {offset}")?;
                }
                write!(
                    formatter,
                    ": {} ({:?})",
                    context.message, context.error_kind
                )
            }
            Self::Codec(context) => {
                write!(formatter, "{}", context.kind)?;
                if let Some(format) = context.format {
                    write!(formatter, " for {format}")?;
                }
                write!(formatter, " while attempting to {}", context.operation)?;
                if let Some(source) = &context.source_label {
                    write!(formatter, " at {source}")?;
                }
                if let Some(frame) = context.frame {
                    write!(formatter, " in frame {frame}")?;
                }
                if let Some(offset) = context.byte_offset {
                    write!(formatter, " at byte {offset}")?;
                }
                if let (Some(expected), Some(actual)) = (context.expected, context.actual) {
                    write!(formatter, " (expected {expected}, actual {actual})")?;
                }
                if let Some(detail) = &context.detail {
                    write!(formatter, ": {detail}")?;
                }
                Ok(())
            }
        }
    }
}

impl std::error::Error for TrajectoryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Conformation(source) => Some(source.as_ref()),
            Self::Realization(source) => Some(source.as_ref()),
            // Codec and I/O contexts already hold the terminal diagnostic facts.
            _ => None,
        }
    }
}

impl From<ConformationError> for TrajectoryError {
    fn from(error: ConformationError) -> Self {
        Self::Conformation(Box::new(error))
    }
}

impl From<kekule::structure::PositionError> for TrajectoryError {
    fn from(error: kekule::structure::PositionError) -> Self {
        Self::Conformation(Box::new(error.into()))
    }
}

impl From<RealizationError> for TrajectoryError {
    fn from(error: RealizationError) -> Self {
        match error {
            RealizationError::TopologyMismatch => Self::TopologyMismatch,
            RealizationError::IndexOutOfRange { index, .. } => {
                Self::FrameIndexOutOfRange(index as u64)
            }
            RealizationError::Conformation(error) => Self::Conformation(Box::new(error)),
            error => Self::Realization(Box::new(error)),
        }
    }
}

impl From<TrajectoryIoErrorContext> for TrajectoryError {
    fn from(context: TrajectoryIoErrorContext) -> Self {
        Self::Io(Box::new(context))
    }
}

impl From<TrajectoryCodecErrorContext> for TrajectoryError {
    fn from(context: TrajectoryCodecErrorContext) -> Self {
        Self::Codec(Box::new(context))
    }
}
