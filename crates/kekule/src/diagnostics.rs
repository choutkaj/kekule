//! Shared diagnostics of format parsing, interpretation and writing.
//!
//! Every format reports failures as one [`FormatError`]: its [`FormatStage`],
//! a structured [`DiagnosticKind`] that callers can match on, and a
//! [`SourceLocation`] in the original input. Messages remain human-readable
//! detail, not an API contract. Writers share [`UnrepresentablePolicy`] for
//! state a format cannot represent.

use std::fmt;
use std::sync::Arc;

pub use crate::io::UnrepresentablePolicy;

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

    pub const fn line(&self) -> usize {
        self.line
    }

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
    /// topologies, models and ensembles.
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
    /// An SDF record title cannot be written.
    Title,
    /// An SDF data field cannot be written, or holds an invalid ensemble
    /// member value.
    DataField {
        /// The field name.
        name: String,
    },
    /// A failure specific to the mmCIF writer.
    MmcifWrite(crate::io::MmcifWriteIssue),
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
        Self::write(
            format,
            DiagnosticKind::Io,
            format!("{format} output failed"),
        )
        .with_source(error)
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

    /// The 1-based number of the record the failure concerns, for
    /// record-oriented input and output such as SDF.
    pub fn record(&self) -> Option<usize> {
        self.0.record
    }

    /// Typed detail of the failure, when it has any.
    pub fn detail(&self) -> Option<&FormatErrorDetail> {
        self.0.detail.as_ref()
    }

    /// Human-readable description of the failure, without its location,
    /// record or source.
    pub fn message(&self) -> &str {
        &self.0.message
    }
}

impl fmt::Display for FormatError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{} {} error", self.0.format, self.0.stage)?;
        if let Some(record) = self.0.record {
            write!(formatter, " in record {record}")?;
        }
        if let Some(location) = self.0.location {
            write!(formatter, " at {location}")?;
        }
        write!(formatter, ": {}", self.0.message)
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

/// Normalizes line terminators of `input` to `\n`, copying it once.
///
/// Every terminator (`\r\n`, `\r`, `\n`) counts as one line, so line numbers in
/// the normalized text equal line numbers in the original input.
pub(crate) fn normalize_line_endings(input: &str) -> String {
    if !input.contains('\r') {
        return input.to_owned();
    }
    let mut normalized = String::with_capacity(input.len());
    let mut chars = input.chars().peekable();
    while let Some(character) = chars.next() {
        if character == '\r' {
            if chars.peek() == Some(&'\n') {
                chars.next();
            }
            normalized.push('\n');
        } else {
            normalized.push(character);
        }
    }
    normalized
}

/// Reads at most `max_bytes` of UTF-8 text from `reader`, failing with a
/// parse-stage [`FormatError`] of `format` without a location.
pub(crate) fn read_bounded(
    reader: impl std::io::Read,
    max_bytes: usize,
    format: Format,
) -> Result<String, FormatError> {
    use std::io::Read;

    let limit = u64::try_from(max_bytes)
        .unwrap_or(u64::MAX)
        .saturating_add(1);
    let mut bytes = Vec::new();
    reader
        .take(limit)
        .read_to_end(&mut bytes)
        .map_err(|error| {
            FormatError::parse(format, DiagnosticKind::Io, "cannot read input").with_source(error)
        })?;
    if bytes.len() > max_bytes {
        return Err(FormatError::parse(
            format,
            DiagnosticKind::Limit,
            "input exceeds configured byte limit",
        ));
    }
    String::from_utf8(bytes).map_err(|error| {
        FormatError::parse(format, DiagnosticKind::Io, "input is not valid UTF-8")
            .with_source(error)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_endings_normalize_without_changing_line_numbers() {
        for (input, expected) in [
            ("a\nb", "a\nb"),
            ("a\r\nb\r\n", "a\nb\n"),
            ("a\rb\r\rc", "a\nb\n\nc"),
            ("\r\n\r", "\n\n"),
        ] {
            assert_eq!(normalize_line_endings(input), expected, "{input:?}");
        }
    }

    #[test]
    fn bounded_reads_report_limits_and_encoding_as_typed_kinds() {
        assert_eq!(
            read_bounded("abc".as_bytes(), 3, Format::Smiles).unwrap(),
            "abc"
        );
        let limit = read_bounded("abcd".as_bytes(), 3, Format::Smiles).unwrap_err();
        assert_eq!(
            (limit.format(), limit.stage(), limit.kind()),
            (Format::Smiles, FormatStage::Parse, DiagnosticKind::Limit)
        );
        let encoding = read_bounded([0xff_u8, 0xfe].as_slice(), 10, Format::Pdb).unwrap_err();
        assert_eq!(encoding.kind(), DiagnosticKind::Io);
        assert!(std::error::Error::source(&encoding)
            .and_then(|source| source.downcast_ref::<std::string::FromUtf8Error>())
            .is_some());
    }

    #[test]
    fn format_errors_display_format_stage_record_and_location_but_not_their_source() {
        let cause = std::io::Error::new(std::io::ErrorKind::BrokenPipe, "closed");
        let error = FormatError::write(Format::Sdf, DiagnosticKind::Io, "SDF output failed")
            .in_record(3)
            .with_source(cause);
        assert_eq!(
            error.to_string(),
            "SDF write error in record 3: SDF output failed"
        );
        assert_eq!(
            std::error::Error::source(&error)
                .and_then(|source| source.downcast_ref::<std::io::Error>())
                .map(std::io::Error::kind),
            Some(std::io::ErrorKind::BrokenPipe)
        );

        let located = FormatError::parse(Format::Smiles, DiagnosticKind::Syntax, "unclosed ring")
            .at(SourceLocation::at_column(1, 4))
            .on_line(7)
            .shift_column(2);
        assert_eq!(located.location(), Some(SourceLocation::at_column(7, 6)));
        assert_eq!(located.record(), None);
        assert_eq!(
            located.to_string(),
            "SMILES parse error at line 7, column 6: unclosed ring"
        );
        assert!(std::error::Error::source(&located).is_none());

        let lined = FormatError::interpret(Format::Pdb, DiagnosticKind::Inconsistent, "duplicate")
            .on_line(2);
        assert_eq!(lined.location(), Some(SourceLocation::at_line(2)));
        assert_eq!(
            lined.to_string(),
            "PDB interpretation error at line 2: duplicate"
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
