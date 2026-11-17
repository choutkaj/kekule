//! The failure of frame-table I/O: [`FrameTableError`].

use std::fmt;
use std::sync::Arc;

use crate::properties::PropertyKey;

/// What kind of failure a [`FrameTableError`] is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum FrameTableErrorKind {
    /// Reading or writing the file failed; the [`std::io::Error`] is the
    /// source.
    Io,
    /// A header or row cannot be written as given, such as a step that does
    /// not follow the last one; [`FrameTableError::column`] names the column
    /// when one is at fault.
    Invalid,
    /// A complete line of the file is not valid;
    /// [`FrameTableError::line`] is its number.
    Corrupt,
    /// The file is a frame table of another version, or has a header record
    /// this version does not know.
    Unsupported,
    /// The file's header differs from the expected one
    /// ([`append_frame_table`](super::append_frame_table)).
    Mismatch,
    /// The file's header is not complete yet, as while its writer starts;
    /// read it again later.
    Pending,
    /// Bytes already read changed, as when a run was cut and continued
    /// ([`FrameTableReader::refresh`](super::FrameTableReader::refresh)).
    /// Open the file again.
    Changed,
}

/// Failure of a frame-table operation.
///
/// Its kind classifies it, its accessors give typed payloads, and its
/// message is written from a typed reason.
#[derive(Debug, Clone)]
pub struct FrameTableError(Box<Data>);

#[derive(Debug, Clone)]
struct Data {
    kind: FrameTableErrorKind,
    reason: Reason,
    line: Option<usize>,
    source: Option<Arc<dyn std::error::Error + Send + Sync>>,
}

/// Why an operation failed, which writes the message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum Reason {
    Read,
    Write,
    ReservedStep,
    RepeatedColumn(PropertyKey),
    UnnamedUnit(PropertyKey),
    NoColumn(PropertyKey),
    NoTimes(PropertyKey),
    NoEnergies(PropertyKey),
    NoStep,
    StepOrder { step: u64, last: u64 },
    CellType(PropertyKey),
    CellUnit(PropertyKey),
    ForeignRow,
    NotFrameTable,
    Version(String),
    UnknownRecord,
    Record,
    RepeatedRecord,
    Header,
    Names,
    Encoding,
    CellCount { found: usize, expected: usize },
    Cell(PropertyKey),
    Step,
    RowOrder,
    Mismatch,
    Pending,
    Changed,
}

impl FrameTableError {
    pub(super) fn new(kind: FrameTableErrorKind, reason: Reason) -> Self {
        Self(Box::new(Data {
            kind,
            reason,
            line: None,
            source: None,
        }))
    }

    pub(super) fn invalid(reason: Reason) -> Self {
        Self::new(FrameTableErrorKind::Invalid, reason)
    }

    /// A complete line `line` (counted from 1) that is not valid.
    pub(super) fn corrupt(line: usize, reason: Reason) -> Self {
        Self::new(FrameTableErrorKind::Corrupt, reason).at_line(line)
    }

    pub(super) fn read(error: std::io::Error) -> Self {
        Self::new(FrameTableErrorKind::Io, Reason::Read).with_source(error)
    }

    pub(super) fn write(error: std::io::Error) -> Self {
        Self::new(FrameTableErrorKind::Io, Reason::Write).with_source(error)
    }

    pub(super) fn at_line(mut self, line: usize) -> Self {
        self.0.line = Some(line);
        self
    }

    pub(super) fn with_source(
        mut self,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> Self {
        self.0.source = Some(Arc::new(source));
        self
    }

    /// The kind of failure.
    pub fn kind(&self) -> FrameTableErrorKind {
        self.0.kind
    }

    /// The number of the line at fault, counted from 1, for a corrupt line
    /// or an unknown header record.
    pub fn line(&self) -> Option<usize> {
        self.0.line
    }

    /// The column at fault, when one is.
    pub fn column(&self) -> Option<&PropertyKey> {
        match &self.0.reason {
            Reason::RepeatedColumn(key)
            | Reason::UnnamedUnit(key)
            | Reason::NoColumn(key)
            | Reason::NoTimes(key)
            | Reason::NoEnergies(key)
            | Reason::CellType(key)
            | Reason::CellUnit(key)
            | Reason::Cell(key) => Some(key),
            _ => None,
        }
    }
}

impl PartialEq for FrameTableError {
    fn eq(&self, other: &Self) -> bool {
        let (left, right) = (&self.0, &other.0);
        left.kind == right.kind
            && left.reason == right.reason
            && left.line == right.line
            && left.source.as_ref().map(ToString::to_string)
                == right.source.as_ref().map(ToString::to_string)
    }
}

impl fmt::Display for FrameTableError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let line = self.0.line.unwrap_or(0);
        match &self.0.reason {
            Reason::Read => formatter.write_str("Cannot read the frame table."),
            Reason::Write => formatter.write_str("Cannot write the frame table."),
            Reason::ReservedStep => formatter.write_str(
                "The column \"step\" holds the steps of a frame table. Use another key.",
            ),
            Reason::RepeatedColumn(key) => {
                write!(formatter, "The frame table has two columns \"{key}\".")
            }
            Reason::UnnamedUnit(key) => write!(
                formatter,
                "The unit of the column \"{key}\" has no symbol. Use a named unit."
            ),
            Reason::NoColumn(key) => {
                write!(formatter, "The frame table has no column \"{key}\".")
            }
            Reason::NoTimes(key) => write!(formatter, "The column \"{key}\" holds no times."),
            Reason::NoEnergies(key) => {
                write!(formatter, "The column \"{key}\" holds no energies.")
            }
            Reason::NoStep => formatter.write_str("The row has no step."),
            Reason::StepOrder { step, last } => {
                write!(formatter, "Step {step} does not follow step {last}.")
            }
            Reason::CellType(key) => write!(
                formatter,
                "The value for the column \"{key}\" has another type."
            ),
            Reason::CellUnit(key) => write!(
                formatter,
                "The value for the column \"{key}\" has a unit of another dimension."
            ),
            Reason::ForeignRow => formatter.write_str("The row belongs to another frame table."),
            Reason::NotFrameTable => formatter.write_str("The file is not a frame table."),
            Reason::Version(version) => write!(
                formatter,
                "The frame table has version \"{version}\". Only version {} is read.",
                super::header::FRAME_TABLE_VERSION
            ),
            Reason::UnknownRecord => write!(
                formatter,
                "Line {line} of the frame table is an unknown header record."
            ),
            Reason::Record => write!(
                formatter,
                "Line {line} of the frame table is not a valid header record."
            ),
            Reason::RepeatedRecord => write!(
                formatter,
                "Line {line} of the frame table repeats a header record."
            ),
            Reason::Header => formatter.write_str("The header of the frame table is not valid."),
            Reason::Names => write!(
                formatter,
                "Line {line} of the frame table does not name its columns."
            ),
            Reason::Encoding => write!(formatter, "Line {line} of the frame table is not UTF-8."),
            Reason::CellCount { found, expected } => write!(
                formatter,
                "Line {line} of the frame table has {found} cells instead of {expected}."
            ),
            Reason::Cell(key) => write!(
                formatter,
                "Line {line} of the frame table has an invalid value in the column \"{key}\"."
            ),
            Reason::Step => write!(
                formatter,
                "Line {line} of the frame table has an invalid step."
            ),
            Reason::RowOrder => write!(
                formatter,
                "The step on line {line} of the frame table does not follow the step before it."
            ),
            Reason::Mismatch => {
                formatter.write_str("The header of the frame table differs from the expected one.")
            }
            Reason::Pending => formatter.write_str("The frame table has no complete header yet."),
            Reason::Changed => {
                formatter.write_str("The frame table changed since it was last read.")
            }
        }
    }
}

impl std::error::Error for FrameTableError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.0
            .source
            .as_deref()
            .map(|source| source as &(dyn std::error::Error + 'static))
    }
}
