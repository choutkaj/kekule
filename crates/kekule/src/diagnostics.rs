//! Shared diagnostics of format parsing, interpretation and writing.
//!
//! Every format reports failures as one [`FormatError`]: its [`FormatStage`],
//! a structured [`DiagnosticKind`] that callers can match on, and a
//! [`SourceLocation`] in the original input. Messages remain human-readable
//! detail, not an API contract. Writers share [`UnrepresentablePolicy`] for
//! state a format cannot represent, and readers share [`UnsupportedContent`]
//! for valid content Kekule does not support.
//!
//! Bounded operations, in formats and elsewhere, report an exhausted resource
//! limit as a [`LimitExceeded`] naming the counted [`Resource`].
//!
//! # Example
//!
//! ```
//! use kekule::diagnostics::{DiagnosticKind, Format, FormatStage};
//!
//! let error = kekule::smiles::to_molecules("C1CC").unwrap_err();
//! assert_eq!(error.format(), Format::Smiles);
//! assert_eq!(error.stage(), FormatStage::Parse);
//! assert_eq!(error.kind(), DiagnosticKind::Syntax);
//! ```

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

/// A quantity a bounded algorithm or reader counts against a limit.
///
/// The error carrying a [`LimitExceeded`] names the operation, so a resource
/// says only what was counted: [`Resource::Atoms`] in a ring perception error
/// counts the atoms ring perception was asked to analyse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Resource {
    /// Bytes of input text.
    InputBytes,
    /// Bytes of one input record or line.
    RecordBytes,
    /// Bytes of one value or token.
    ValueBytes,
    /// Input records.
    Records,
    /// Data blocks of a CIF document.
    Blocks,
    /// Rows of one table, such as `_atom_site` rows.
    Rows,
    /// Values of one table or document.
    Values,
    /// Atoms, or atom storage slots.
    Atoms,
    /// Combined atom and bond storage slots.
    GraphSlots,
    /// Pending atom labels of canonical labeling.
    Labels,
    /// Bonds, or bond storage slots.
    Bonds,
    /// Hydrogen atoms an operation would add.
    Hydrogens,
    /// Molecule components or instances.
    Components,
    /// Nodes of a query expression or of all expressions of a query.
    ExpressionNodes,
    /// Nesting depth of a query expression.
    ExpressionDepth,
    /// Nesting depth of recursive queries.
    RecursionDepth,
    /// Nesting depth of branches.
    BranchDepth,
    /// Ring closures, open at once or in total.
    RingClosures,
    /// Units of work of a bounded calculation, such as graph visits.
    Work,
    /// States of a search, such as partial mappings or labelings.
    SearchStates,
    /// Candidate cycles of ring perception.
    CandidateCycles,
    /// Path expansions of ring perception.
    PathExpansions,
    /// Equivalent shortest paths of ring perception.
    ShortestPaths,
    /// Ring size, in atoms.
    RingSize,
    /// Matches reported by a search.
    Matches,
    /// Permutations of symmetry-equivalent atoms.
    Permutations,
    /// Nodes of a CIP hierarchical digraph.
    DigraphNodes,
    /// Depth of a CIP hierarchical digraph.
    DigraphDepth,
    /// Resonance contributors.
    Structures,
    /// Periodic lattice images examined.
    LatticeImages,
    /// Candidate positions, such as solvent placement sites.
    Positions,
    /// Pairwise entries of a distance matrix.
    MatrixEntries,
    /// Residues, such as protein residues of a secondary-structure assignment.
    Residues,
    /// Frames of a trajectory.
    Frames,
    /// Fields or items of one record.
    Fields,
    /// Combinations of carriers or other choices an encoder considers.
    Combinations,
    /// Candidate pairs, such as query-target atom pairs or residue pairs.
    CandidatePairs,
    /// Beta ladders of a secondary-structure assignment.
    Ladders,
}

impl fmt::Display for Resource {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::InputBytes => "input bytes",
            Self::RecordBytes => "record bytes",
            Self::ValueBytes => "value bytes",
            Self::Records => "records",
            Self::Blocks => "data blocks",
            Self::Rows => "rows",
            Self::Values => "values",
            Self::Atoms => "atoms",
            Self::GraphSlots => "graph slots",
            Self::Labels => "pending labels",
            Self::Bonds => "bonds",
            Self::Hydrogens => "hydrogens",
            Self::Components => "components",
            Self::ExpressionNodes => "expression nodes",
            Self::ExpressionDepth => "expression depth",
            Self::RecursionDepth => "recursion depth",
            Self::BranchDepth => "branch depth",
            Self::RingClosures => "ring closures",
            Self::Work => "work",
            Self::SearchStates => "search states",
            Self::CandidateCycles => "candidate cycles",
            Self::PathExpansions => "path expansions",
            Self::ShortestPaths => "equivalent shortest paths",
            Self::RingSize => "ring size",
            Self::Matches => "matches",
            Self::Permutations => "permutations",
            Self::DigraphNodes => "digraph nodes",
            Self::DigraphDepth => "digraph depth",
            Self::Structures => "resonance structures",
            Self::LatticeImages => "lattice images",
            Self::Positions => "candidate positions",
            Self::MatrixEntries => "matrix entries",
            Self::Residues => "residues",
            Self::Frames => "frames",
            Self::Fields => "fields",
            Self::Combinations => "combinations",
            Self::CandidatePairs => "candidate pairs",
            Self::Ladders => "ladders",
        })
    }
}

/// A bounded operation stopped because it would exceed a resource limit.
///
/// `observed` is the value that exceeded `limit`. An operation that stops at
/// the first excess reports that first value, so `observed` may be smaller
/// than the input would eventually have required; it is always greater than
/// `limit`. No partial result accompanies the error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LimitExceeded {
    resource: Resource,
    limit: usize,
    observed: usize,
}

impl LimitExceeded {
    /// Reports that `observed` units of `resource` exceed `limit`.
    pub const fn new(resource: Resource, limit: usize, observed: usize) -> Self {
        Self {
            resource,
            limit,
            observed,
        }
    }

    /// Reports the first unit of `resource` beyond `limit`, for an operation
    /// that counts one unit at a time.
    pub const fn first_beyond(resource: Resource, limit: usize) -> Self {
        Self::new(resource, limit, limit.saturating_add(1))
    }

    /// Returns `Err` when `observed` exceeds `limit`.
    pub(crate) const fn check(
        resource: Resource,
        limit: usize,
        observed: usize,
    ) -> Result<(), Self> {
        if observed > limit {
            Err(Self::new(resource, limit, observed))
        } else {
            Ok(())
        }
    }

    /// The bounded resource.
    pub const fn resource(&self) -> Resource {
        self.resource
    }

    /// The configured or built-in bound.
    pub const fn limit(&self) -> usize {
        self.limit
    }

    /// The value that exceeded the bound.
    pub const fn observed(&self) -> usize {
        self.observed
    }
}

impl fmt::Display for LimitExceeded {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{} limit exceeded: observed {}, limit {}",
            self.resource, self.observed, self.limit
        )
    }
}

impl std::error::Error for LimitExceeded {}

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
    /// The written structure is an ensemble member or a trajectory frame
    /// with this state, which the format cannot carry, and the writer's
    /// policy for it is [`UnrepresentablePolicy::Fail`].
    UnrepresentableState(StructureField),
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
}

/// Ensemble member or trajectory frame state named by
/// [`FormatErrorDetail::UnrepresentableState`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StructureField {
    /// An ensemble member's relative weight.
    Weight,
    /// An ensemble member's energy.
    Energy,
    /// A trajectory frame's time.
    Time,
    /// A trajectory frame's step.
    Step,
    /// A trajectory frame's velocities.
    Velocities,
    /// A trajectory frame's forces.
    Forces,
}

impl fmt::Display for StructureField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Weight => "ensemble member weight",
            Self::Energy => "ensemble member energy",
            Self::Time => "trajectory frame time",
            Self::Step => "trajectory frame step",
            Self::Velocities => "trajectory frame velocities",
            Self::Forces => "trajectory frame forces",
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
        if let Some(offset) = self.0.byte_offset {
            write!(formatter, " at byte {offset}")?;
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
        )
        .exceeding(LimitExceeded::new(
            Resource::InputBytes,
            max_bytes,
            bytes.len(),
        )));
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
    fn limits_report_resource_bound_and_first_excess() {
        assert_eq!(LimitExceeded::check(Resource::Atoms, 3, 3), Ok(()));
        let exceeded = LimitExceeded::check(Resource::Atoms, 3, 5).unwrap_err();
        assert_eq!(
            (exceeded.resource(), exceeded.limit(), exceeded.observed()),
            (Resource::Atoms, 3, 5)
        );
        assert_eq!(
            exceeded.to_string(),
            "atoms limit exceeded: observed 5, limit 3"
        );
        assert_eq!(
            LimitExceeded::first_beyond(Resource::Records, usize::MAX).observed(),
            usize::MAX
        );
        let error = FormatError::parse(Format::Sdf, DiagnosticKind::Limit, "too many records")
            .exceeding(LimitExceeded::first_beyond(Resource::Records, 2));
        assert_eq!(
            error.limit_exceeded(),
            Some(LimitExceeded::new(Resource::Records, 2, 3))
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
