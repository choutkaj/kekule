//! Fixed-topology trajectory frames, reusable buffers, in-memory storage, and
//! streaming reader/writer contracts.
//!
//! [`TrajectoryFrame`] is topology-free realization payload, [`Trajectory`]
//! supplies the shared topology and temporal ordering, and [`FrameBuffer`]
//! provides reusable validated storage for streaming decoders. The reader and
//! writer traits publish complete frames transactionally.

macro_rules! realization_property_api {
    () => {
        pub const fn atom_properties(&self) -> &PropertyTable {
            self.properties.realization_atom_properties()
        }

        pub const fn bond_properties(&self) -> &PropertyTable {
            self.properties.realization_bond_properties()
        }

        pub fn atom_property(
            &self,
            index: usize,
            key: &PropertyKey,
        ) -> Result<Option<PropertyValue>, FrameError> {
            Ok(self.atom_properties().value(key, index)?)
        }

        pub fn set_atom_property(
            &mut self,
            index: usize,
            key: PropertyKey,
            value: Option<PropertyValue>,
        ) -> Result<(), FrameError> {
            Ok(self
                .properties
                .set_realization_atom_value(key, index, value)?)
        }

        pub fn insert_atom_property_column(
            &mut self,
            key: PropertyKey,
            column: PropertyColumn,
        ) -> Result<Option<PropertyColumn>, FrameError> {
            Ok(self
                .properties
                .insert_realization_atom_column(key, column)?)
        }

        pub fn remove_atom_property_column(
            &mut self,
            key: &PropertyKey,
        ) -> Result<Option<PropertyColumn>, FrameError> {
            Ok(self.properties.remove_realization_atom_column(key)?)
        }

        pub fn bond_property(
            &self,
            index: usize,
            key: &PropertyKey,
        ) -> Result<Option<PropertyValue>, FrameError> {
            Ok(self.bond_properties().value(key, index)?)
        }

        pub fn set_bond_property(
            &mut self,
            index: usize,
            key: PropertyKey,
            value: Option<PropertyValue>,
        ) -> Result<(), FrameError> {
            Ok(self
                .properties
                .set_realization_bond_value(key, index, value)?)
        }

        pub fn insert_bond_property_column(
            &mut self,
            key: PropertyKey,
            column: PropertyColumn,
        ) -> Result<Option<PropertyColumn>, FrameError> {
            Ok(self
                .properties
                .insert_realization_bond_column(key, column)?)
        }

        pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
            self.properties.remove_realization_bond_column(key)
        }

        pub fn occupancy_at(&self, index: usize) -> Result<Option<f64>, FrameError> {
            Ok(self.properties.occupancy_at(index)?)
        }

        pub fn set_occupancy_at(
            &mut self,
            index: usize,
            value: Option<f64>,
        ) -> Result<(), FrameError> {
            Ok(self.properties.set_occupancy_at(index, value)?)
        }

        pub fn b_factor_at(&self, index: usize) -> Result<Option<Quantity<f64>>, FrameError> {
            Ok(self.properties.b_factor_at(index)?)
        }

        pub fn set_b_factor_at(
            &mut self,
            index: usize,
            value: Option<Quantity<f64>>,
        ) -> Result<(), FrameError> {
            Ok(self.properties.set_b_factor_at(index, value)?)
        }
    };
}

mod buffer;
mod collection;
mod frame;
mod stream;

#[cfg(test)]
mod tests;

pub use buffer::{FrameBuffer, FrameBufferData};
pub use collection::Trajectory;
pub use frame::{Forces, TrajectoryFrame, TrajectoryFrameView, Velocities};
pub use stream::{
    validate_atom_order, CoordinateFrameReader, MemoryTrajectoryReader, MemoryTrajectoryWriter,
    SeekableTrajectoryReader, TrajectoryReader, TrajectoryWriter,
};

use std::{fmt, io};

use kekule::properties::PropertyError;
use kekule::structure::{ModelError, PositionError};
use kekule::topology::transform::TopologySubsetError;
use kekule::units::UnitError;

fn validate_atom_count(expected: usize, actual: usize) -> Result<(), FrameError> {
    if actual != expected {
        return Err(FrameError::AtomCountMismatch { expected, actual });
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum FrameError {
    TopologyMismatch,
    InvalidIndex { index: usize },
    AtomCountMismatch { expected: usize, actual: usize },
    BondCountMismatch { expected: usize, actual: usize },
    NonFiniteVector { index: usize },
    NonFiniteTime,
    Position(PositionError),
    Property(PropertyError),
    Unit(UnitError),
}

impl fmt::Display for FrameError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("trajectory frame belongs to a different topology")
            }
            Self::InvalidIndex { index } => write!(formatter, "invalid dense frame index {index}"),
            Self::AtomCountMismatch { expected, actual } => write!(
                formatter,
                "trajectory array requires {expected} atoms, but received {actual}"
            ),
            Self::BondCountMismatch { expected, actual } => write!(
                formatter,
                "trajectory array requires {expected} bonds, but received {actual}"
            ),
            Self::NonFiniteVector { index } => {
                write!(formatter, "trajectory vector at {index} is not finite")
            }
            Self::NonFiniteTime => formatter.write_str("trajectory time must be finite"),
            Self::Position(error) => write!(formatter, "invalid frame positions: {error}"),
            Self::Property(error) => write!(formatter, "invalid frame property: {error}"),
            Self::Unit(error) => write!(formatter, "invalid frame quantity unit: {error}"),
        }
    }
}

impl std::error::Error for FrameError {}

impl From<PositionError> for FrameError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<PropertyError> for FrameError {
    fn from(error: PropertyError) -> Self {
        Self::Property(error)
    }
}

impl From<UnitError> for FrameError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

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
    MissingTime { frame: usize },
    NonMonotonicTime { frame: usize },
    UnsupportedField(&'static str),
    Frame(Box<FrameError>),
    Position(PositionError),
    Io(Box<TrajectoryIoErrorContext>),
    Codec(Box<TrajectoryCodecErrorContext>),
}

/// Failure to subset a trajectory topology or transfer frame state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TrajectorySliceError {
    Topology(TopologySubsetError),
    Position(PositionError),
    Property(PropertyError),
    Frame(Box<FrameError>),
    Trajectory(TrajectoryError),
}

impl fmt::Display for TrajectorySliceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "cannot slice trajectory: {self:?}")
    }
}

impl std::error::Error for TrajectorySliceError {}

impl From<TopologySubsetError> for TrajectorySliceError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Topology(error)
    }
}
impl From<PositionError> for TrajectorySliceError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}
impl From<PropertyError> for TrajectorySliceError {
    fn from(error: PropertyError) -> Self {
        Self::Property(error)
    }
}
impl From<FrameError> for TrajectorySliceError {
    fn from(error: FrameError) -> Self {
        Self::Frame(Box::new(error))
    }
}
impl From<TrajectoryError> for TrajectorySliceError {
    fn from(error: TrajectoryError) -> Self {
        Self::Trajectory(error)
    }
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
            Self::MissingTime { frame } => {
                write!(formatter, "trajectory frame {frame} has no time")
            }
            Self::NonMonotonicTime { frame } => {
                write!(formatter, "trajectory time decreases at frame {frame}")
            }
            Self::UnsupportedField(field) => {
                write!(formatter, "trajectory writer does not support {field}")
            }
            Self::Frame(error) => write!(formatter, "invalid trajectory frame: {error}"),
            Self::Position(error) => write!(formatter, "invalid trajectory positions: {error}"),
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

impl std::error::Error for TrajectoryError {}

impl From<FrameError> for TrajectoryError {
    fn from(error: FrameError) -> Self {
        match error {
            FrameError::TopologyMismatch => Self::TopologyMismatch,
            error => Self::Frame(Box::new(error)),
        }
    }
}

impl From<PositionError> for TrajectoryError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<ModelError> for TrajectoryError {
    fn from(_: ModelError) -> Self {
        Self::TopologyMismatch
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
