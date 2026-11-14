use std::fmt;
use std::sync::Arc;

use super::LimitExceeded;

/// A position in the original input text.
///
/// Lines are 1-based and count every line terminator of the original input
/// (`\n`, `\r\n`, or `\r`). Columns, when known, are 1-based byte offsets within
/// the line.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SourceLocation {
    line: usize,
    column: Option<usize>,
}

impl SourceLocation {
    /// A location identified by its line alone.
    pub const fn at_line(line: usize) -> Self {
        Self { line, column: None }
    }

    /// A location identified by line and byte column.
    pub const fn at_column(line: usize, column: usize) -> Self {
        Self {
            line,
            column: Some(column),
        }
    }

    /// The 1-based line number.
    pub const fn line(&self) -> usize {
        self.line
    }

    /// The 1-based byte column within the line, or `None` when only the line
    /// is known.
    pub const fn column(&self) -> Option<usize> {
        self.column
    }
}

impl fmt::Display for SourceLocation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.column {
            Some(column) => write!(formatter, "line {}, column {column}", self.line),
            None => write!(formatter, "line {}", self.line),
        }
    }
}

/// What kind of problem a parse or interpretation diagnostic reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum DiagnosticKind {
    /// Reading the input failed, or it is not valid UTF-8.
    Io,
    /// The input exceeds a configured resource limit.
    Limit,
    /// The text does not follow the format grammar.
    Syntax,
    /// A well-formed field holds a value outside its domain, such as an unknown
    /// element symbol or a malformed number.
    InvalidValue,
    /// A field refers to something that does not exist, such as a bond
    /// endpoint outside the atom block.
    InvalidReference,
    /// Valid source content that Kekule does not support.
    Unsupported,
    /// Fields or records contradict each other, or a requested selection does
    /// not match the source.
    Inconsistent,
    /// The represented chemistry cannot be constructed or published.
    Chemistry,
    /// A writer requires perception the source does not have installed;
    /// perceive it and retry.
    MissingPerception,
    /// The caller's [`Cancellation`](super::Cancellation) stopped the operation; the error
    /// carries [`Cancelled`](super::Cancelled) as its source and no partial result.
    Cancelled,
}

impl fmt::Display for DiagnosticKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Io => "input or output error",
            Self::Limit => "resource limit",
            Self::Syntax => "syntax error",
            Self::InvalidValue => "invalid value",
            Self::InvalidReference => "invalid reference",
            Self::Unsupported => "unsupported content",
            Self::Inconsistent => "inconsistent content",
            Self::Chemistry => "chemistry error",
            Self::MissingPerception => "missing perception",
            Self::Cancelled => "cancelled",
        })
    }
}

/// The format a [`FormatError`] concerns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Format {
    /// SMILES and CXSMILES.
    Smiles,
    /// SMARTS queries.
    Smarts,
    /// V2000 and V3000 Molfiles.
    Molfile,
    /// SD files.
    Sdf,
    /// CIF/mmCIF, including chemical component dictionaries.
    Mmcif,
    /// PDB.
    Pdb,
    /// Kekule's native binary serialization.
    Native,
    /// Kekule's atom selection language.
    Selection,
    /// XYZ coordinate files.
    Xyz,
}

impl fmt::Display for Format {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Smiles => "SMILES",
            Self::Smarts => "SMARTS",
            Self::Molfile => "Molfile",
            Self::Sdf => "SDF",
            Self::Mmcif => "mmCIF",
            Self::Pdb => "PDB",
            Self::Native => "native",
            Self::Selection => "selection",
            Self::Xyz => "XYZ",
        })
    }
}

/// The stage of format processing at which a [`FormatError`] occurred.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum FormatStage {
    /// Reading text into a typed document, without assigning chemical meaning.
    Parse,
    /// Assigning canonical meaning to a parsed document: molecules,
    /// topologies, structures and structure sets.
    Interpret,
    /// Encoding canonical values as text, including writing to the sink.
    Write,
}

impl fmt::Display for FormatStage {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Parse => "parse",
            Self::Interpret => "interpretation",
            Self::Write => "write",
        })
    }
}

/// Typed, matchable detail of a [`FormatError`] beyond its kind.
///
/// Most failures carry no detail: their [`DiagnosticKind`], location and
/// message describe them. A variant names the failures callers can act on.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum FormatErrorDetail {
    /// A single-molecule accessor of a SMILES interpretation found `actual`
    /// connected components instead of one.
    ComponentCount {
        /// The number of components the interpretation holds.
        actual: usize,
    },
    /// The input contains no coordinate models.
    NoCoordinateModels,
    /// An explicit coordinate-model selection is empty.
    EmptyModelSelection,
    /// An explicit conformation selection is empty.
    EmptyConformationSelection,
    /// Document-level interpretation found atom-site data in more than one
    /// data block.
    MultipleAtomSiteBlocks,
    /// A coordinate model was requested more than once.
    DuplicateRequestedModel {
        /// The requested model ID.
        model_id: String,
    },
    /// A requested coordinate model is not present.
    UnknownRequestedModel {
        /// The requested model ID.
        model_id: String,
    },
    /// Interpreting this coordinate model failed; the error's kind, location
    /// and message describe the failure within it.
    Model {
        /// The model ID.
        model_id: String,
    },
    /// Interpreting the zero-based caller-supplied conformation selection
    /// `selection` failed; the error's kind, location and message describe
    /// the failure within it.
    Conformation {
        /// The zero-based selection index.
        selection: usize,
    },
    /// This coordinate model has a molecule partition, chemistry,
    /// connectivity or hierarchy inconsistent with the first one.
    InconsistentTopology {
        /// The model ID.
        model_id: String,
    },
    /// This coordinate model has an atom identity set inconsistent with the
    /// first one.
    InconsistentAtomSet {
        /// The model ID.
        model_id: String,
    },
    /// This coordinate model orders its atom identities differently from the
    /// first one.
    InconsistentDenseAtomOrder {
        /// The model ID.
        model_id: String,
    },
    /// A record title or name cannot be written: an SDF or Molfile title, or
    /// a SMILES record name.
    Title,
    /// An SDF data field cannot be written, or holds an invalid ensemble
    /// member value.
    DataField {
        /// The field name.
        name: String,
    },
    /// A failure specific to the mmCIF writer.
    MmcifWrite(crate::io::MmcifWriteIssue),
    /// A failure specific to the PDB writer.
    PdbWrite(crate::io::PdbWriteIssue),
    /// The written structure, geometry, or trajectory frame has this state,
    /// which the format cannot carry, and the writer's policy for it is
    /// [`UnrepresentablePolicy::Fail`](super::UnrepresentablePolicy::Fail).
    UnrepresentableState(GeometryField),
    /// A configured or built-in resource limit was exceeded; the kind is
    /// [`DiagnosticKind::Limit`].
    LimitExceeded(LimitExceeded),
    /// The input does not start with the native format's magic bytes; the
    /// kind is [`DiagnosticKind::Syntax`].
    NativeMagic,
    /// The input is written in a native format version this reader does not
    /// implement; the kind is [`DiagnosticKind::Unsupported`].
    NativeVersion {
        /// The version the input declares.
        found: u32,
        /// The version this reader implements.
        supported: u32,
    },
    /// The native input's checksum does not match its content; the kind is
    /// [`DiagnosticKind::Inconsistent`].
    NativeChecksum {
        /// The checksum the input stores.
        stored: u32,
        /// The checksum of the content read.
        computed: u32,
    },
    /// A native input holds a different kind of object than the projection
    /// requested; the kind is [`DiagnosticKind::Inconsistent`].
    NativeKind {
        /// The kind the projection requires.
        expected: crate::io::NativeKind,
        /// The kind the input holds.
        found: crate::io::NativeKind,
    },
    /// A native structure set record stores its geometries elsewhere, so it has
    /// none to give; the kind is [`DiagnosticKind::Inconsistent`].
    NativeExternalGeometries {
        /// The number of geometries the record describes.
        geometries: usize,
    },
}

/// Geometry, structure set, or trajectory frame state named by
/// [`FormatErrorDetail::UnrepresentableState`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GeometryField {
    /// A structure set's time order without times or steps.
    Order,
    /// A geometry's or frame's time in a time order.
    Time,
    /// A geometry's or frame's MD step in a time order.
    Step,
    /// A geometry's or frame's velocities.
    Velocities,
    /// A geometry's or frame's forces.
    Forces,
    /// A geometry's or frame's per-atom properties.
    AtomProperties,
    /// A geometry's or frame's per-bond properties.
    BondProperties,
    /// A geometry's or frame's own properties.
    Properties,
    /// A structure set's own properties.
    SetProperties,
    /// A structure's crystallographic occupancies.
    Occupancies,
    /// A structure's crystallographic B factors.
    BFactors,
    /// A geometry's step in a step order.
    OrderStep,
    /// A geometry's coordinate in a coordinate order.
    OrderCoordinate,
}

impl fmt::Display for GeometryField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Order => "geometry order",
            Self::Time => "time",
            Self::Step => "MD step",
            Self::Velocities => "velocities",
            Self::Forces => "forces",
            Self::AtomProperties => "per-atom geometry properties",
            Self::BondProperties => "per-bond geometry properties",
            Self::Properties => "geometry properties",
            Self::SetProperties => "structure set properties",
            Self::Occupancies => "occupancies",
            Self::BFactors => "B factors",
            Self::OrderStep => "order step",
            Self::OrderCoordinate => "order coordinate",
        })
    }
}

/// A failure to parse, interpret or write a format, SMARTS included.
///
/// [`kind`](Self::kind) is the stable classification callers match on;
/// [`detail`](Self::detail) adds typed data for the failures that have it.
/// The [message](Self::message) is human-readable detail, not an API
/// contract. When another error caused this one, [`source`](std::error::Error::source)
/// returns it and `Display` does not repeat it.
#[derive(Debug, Clone)]
pub struct FormatError(Box<FormatErrorData>);

#[derive(Debug, Clone)]
struct FormatErrorData {
    format: Format,
    stage: FormatStage,
    kind: DiagnosticKind,
    location: Option<SourceLocation>,
    byte_offset: Option<u64>,
    span: Option<std::ops::Range<usize>>,
    record: Option<usize>,
    message: String,
    detail: Option<FormatErrorDetail>,
    source: Option<Arc<dyn std::error::Error + Send + Sync>>,
}

impl FormatError {
    pub(crate) fn new(
        format: Format,
        stage: FormatStage,
        kind: DiagnosticKind,
        message: impl Into<String>,
    ) -> Self {
        Self(Box::new(FormatErrorData {
            format,
            stage,
            kind,
            location: None,
            byte_offset: None,
            span: None,
            record: None,
            message: message.into(),
            detail: None,
            source: None,
        }))
    }

    pub(crate) fn parse(format: Format, kind: DiagnosticKind, message: impl Into<String>) -> Self {
        Self::new(format, FormatStage::Parse, kind, message)
    }

    pub(crate) fn interpret(
        format: Format,
        kind: DiagnosticKind,
        message: impl Into<String>,
    ) -> Self {
        Self::new(format, FormatStage::Interpret, kind, message)
    }

    pub(crate) fn write(format: Format, kind: DiagnosticKind, message: impl Into<String>) -> Self {
        Self::new(format, FormatStage::Write, kind, message)
    }

    /// A failure of the output sink, keeping `error` as the source.
    pub(crate) fn sink(format: Format, error: std::io::Error) -> Self {
        Self::write(format, DiagnosticKind::Io, "").with_source(error)
    }

    #[must_use]
    pub(crate) fn at(mut self, location: SourceLocation) -> Self {
        self.0.location = Some(location);
        self
    }

    #[must_use]
    pub(crate) fn with_location(mut self, location: Option<SourceLocation>) -> Self {
        self.0.location = location;
        self
    }

    /// Moves the location to `line`, keeping its column.
    #[must_use]
    pub(crate) fn on_line(mut self, line: usize) -> Self {
        self.0.location = Some(match self.0.location.and_then(|location| location.column) {
            Some(column) => SourceLocation::at_column(line, column),
            None => SourceLocation::at_line(line),
        });
        self
    }

    /// Shifts the column by `offset` bytes, for an error located within a
    /// nested part of its record.
    #[must_use]
    pub(crate) fn shift_column(mut self, offset: usize) -> Self {
        if let Some(SourceLocation {
            line,
            column: Some(column),
        }) = self.0.location
        {
            self.0.location = Some(SourceLocation::at_column(line, column + offset));
        }
        self
    }

    /// Locates the failure at a byte offset of binary input.
    #[must_use]
    pub(crate) fn at_byte(mut self, offset: u64) -> Self {
        self.0.byte_offset = Some(offset);
        self
    }

    /// Records the byte range of the input the failure concerns.
    #[must_use]
    pub(crate) fn with_span(mut self, span: std::ops::Range<usize>) -> Self {
        self.0.span = Some(span);
        self
    }

    #[must_use]
    pub(crate) fn in_record(mut self, record: usize) -> Self {
        self.0.record = Some(record);
        self
    }

    #[must_use]
    pub(crate) fn in_format(mut self, format: Format) -> Self {
        self.0.format = format;
        self
    }

    #[must_use]
    pub(crate) fn with_detail(mut self, detail: FormatErrorDetail) -> Self {
        self.0.detail = Some(detail);
        self
    }

    /// Marks the error as exhaustion of `limit`.
    #[must_use]
    pub(crate) fn exceeding(self, limit: LimitExceeded) -> Self {
        self.with_detail(FormatErrorDetail::LimitExceeded(limit))
    }

    /// The exceeded limit, for a [`DiagnosticKind::Limit`] failure that
    /// counts one.
    pub fn limit_exceeded(&self) -> Option<LimitExceeded> {
        match self.detail() {
            Some(FormatErrorDetail::LimitExceeded(limit)) => Some(*limit),
            _ => None,
        }
    }

    #[must_use]
    pub(crate) fn with_source(
        mut self,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> Self {
        self.0.source = Some(Arc::new(source));
        self
    }

    /// The format processed when the failure occurred.
    pub fn format(&self) -> Format {
        self.0.format
    }

    /// Whether parsing, interpretation or writing failed.
    pub fn stage(&self) -> FormatStage {
        self.0.stage
    }

    /// The stable classification of the failure.
    pub fn kind(&self) -> DiagnosticKind {
        self.0.kind
    }

    /// The input location the failure concerns, when it concerns specific
    /// source text. SMILES and SMARTS locations carry the byte column within
    /// the record; other formats locate a line.
    pub fn location(&self) -> Option<SourceLocation> {
        self.0.location
    }

    /// The 0-based byte offset of binary input the failure concerns, for
    /// the native format, whose input has no lines.
    pub fn byte_offset(&self) -> Option<u64> {
        self.0.byte_offset
    }

    /// The byte range of the input the failure concerns, for inputs parsed
    /// as one expression, such as selection-language text; the location
    /// gives its start.
    pub fn span(&self) -> Option<std::ops::Range<usize>> {
        self.0.span.clone()
    }

    /// The 1-based number of the record the failure concerns, for
    /// record-oriented input and output such as SDF.
    pub fn record(&self) -> Option<usize> {
        self.0.record
    }

    /// Typed detail of the failure, when it has any.
    pub fn detail(&self) -> Option<&FormatErrorDetail> {
        self.0.detail.as_ref()
    }

    /// What went wrong, in one or two sentences, without the format, stage,
    /// location, record or source; empty when the source says it all.
    pub fn message(&self) -> &str {
        &self.0.message
    }
}

/// Writes the stage that failed, the format and the place, such as `Cannot
/// parse the SDF input in record 3 at line 12.`, followed by the message.
impl fmt::Display for FormatError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let format = self.0.format;
        match self.0.stage {
            FormatStage::Parse => write!(formatter, "Cannot parse the {format} input")?,
            FormatStage::Interpret => write!(formatter, "Cannot interpret the {format} input")?,
            FormatStage::Write => write!(formatter, "Cannot write the {format} output")?,
        }
        if let Some(record) = self.0.record {
            write!(formatter, " in record {record}")?;
        }
        if let Some(location) = self.0.location {
            write!(formatter, " at {location}")?;
        }
        if let Some(offset) = self.0.byte_offset {
            write!(formatter, " at byte {offset}")?;
        }
        formatter.write_str(".")?;
        if !self.0.message.is_empty() {
            write!(formatter, " {}", self.0.message)?;
        }
        Ok(())
    }
}

impl std::error::Error for FormatError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.0
            .source
            .as_deref()
            .map(|source| source as &(dyn std::error::Error + 'static))
    }
}

/// How interpretation treats source content it does not support.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum UnsupportedContent {
    /// Fail with a [`DiagnosticKind::Unsupported`] error.
    #[default]
    Reject,
    /// Omit the content and record a typed report entry naming it.
    Report,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_errors_display_format_stage_record_and_location_but_not_their_source() {
        let cause = std::io::Error::new(std::io::ErrorKind::BrokenPipe, "closed");
        let error = FormatError::write(Format::Sdf, DiagnosticKind::Io, "")
            .in_record(3)
            .with_source(cause);
        assert_eq!(
            error.to_string(),
            "Cannot write the SDF output in record 3."
        );
        assert_eq!(
            std::error::Error::source(&error)
                .and_then(|source| source.downcast_ref::<std::io::Error>())
                .map(std::io::Error::kind),
            Some(std::io::ErrorKind::BrokenPipe)
        );

        let located = FormatError::parse(
            Format::Smiles,
            DiagnosticKind::Syntax,
            "Ring bond 1 is not closed.",
        )
        .at(SourceLocation::at_column(1, 4))
        .on_line(7)
        .shift_column(2);
        assert_eq!(located.location(), Some(SourceLocation::at_column(7, 6)));
        assert_eq!(located.record(), None);
        assert_eq!(
            located.to_string(),
            "Cannot parse the SMILES input at line 7, column 6. Ring bond 1 is not closed."
        );
        assert!(std::error::Error::source(&located).is_none());

        let lined = FormatError::interpret(
            Format::Pdb,
            DiagnosticKind::Inconsistent,
            "Atom 3 appears twice.",
        )
        .on_line(2);
        assert_eq!(lined.location(), Some(SourceLocation::at_line(2)));
        assert_eq!(
            lined.to_string(),
            "Cannot interpret the PDB input at line 2. Atom 3 appears twice."
        );
    }

    #[test]
    fn locations_display_line_and_optional_column() {
        assert_eq!(SourceLocation::at_line(4).to_string(), "line 4");
        assert_eq!(
            SourceLocation::at_column(4, 7).to_string(),
            "line 4, column 7"
        );
    }
}
