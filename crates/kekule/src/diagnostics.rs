//! Shared diagnostics of format parsing, interpretation and writing.
//!
//! Every format reports failures as one [`FormatError`]: its [`FormatStage`],
//! a structured [`DiagnosticKind`] that callers can match on, and a
//! [`SourceLocation`] in the original input. Messages remain human-readable
//! detail, not an API contract. Writers share [`UnrepresentablePolicy`] for
//! state a format cannot represent, and structure writers return a
//! [`WriteReport`] of what that policy left out; readers share [`UnsupportedContent`]
//! for valid content Kekule does not support.
//!
//! Bounded operations, in formats and elsewhere, report an exhausted resource
//! limit as a [`LimitExceeded`] naming the counted [`Resource`].
//!
//! Long-running operations that accept a [`Cancellation`] stop when its check
//! turns true and report [`Cancelled`]; in a [`FormatError`] the kind is
//! [`DiagnosticKind::Cancelled`].
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

pub use crate::io::{UnrepresentablePolicy, WriteOmission, WriteReport};

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
    /// The caller's [`Cancellation`] stopped the operation; the error
    /// carries [`Cancelled`] as its source and no partial result.
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

impl Resource {
    /// The name of the limit on this resource, such as `work` in "work
    /// limit".
    const fn limit_name(self) -> &'static str {
        match self {
            Self::InputBytes => "input size",
            Self::RecordBytes => "record size",
            Self::ValueBytes => "value size",
            Self::Records => "record",
            Self::Blocks => "data block",
            Self::Rows => "row",
            Self::Values => "value",
            Self::Atoms => "atom",
            Self::GraphSlots => "graph slot",
            Self::Labels => "pending label",
            Self::Bonds => "bond",
            Self::Hydrogens => "hydrogen",
            Self::Components => "component",
            Self::ExpressionNodes => "expression node",
            Self::ExpressionDepth => "expression depth",
            Self::RecursionDepth => "recursion depth",
            Self::BranchDepth => "branch depth",
            Self::RingClosures => "ring closure",
            Self::Work => "work",
            Self::SearchStates => "search state",
            Self::CandidateCycles => "candidate cycle",
            Self::PathExpansions => "path expansion",
            Self::ShortestPaths => "shortest path",
            Self::RingSize => "ring size",
            Self::Matches => "match",
            Self::Permutations => "permutation",
            Self::DigraphNodes => "digraph node",
            Self::DigraphDepth => "digraph depth",
            Self::Structures => "resonance structure",
            Self::LatticeImages => "lattice image",
            Self::Positions => "candidate position",
            Self::MatrixEntries => "matrix entry",
            Self::Residues => "residue",
            Self::Frames => "frame",
            Self::Fields => "field",
            Self::Combinations => "combination",
            Self::CandidatePairs => "candidate pair",
            Self::Ladders => "ladder",
        }
    }

    /// The unit a count of this resource is given in.
    const fn unit(self) -> &'static str {
        match self {
            Self::InputBytes | Self::RecordBytes | Self::ValueBytes => "bytes",
            Self::Records => "records",
            Self::Blocks => "data blocks",
            Self::Rows => "rows",
            Self::Values => "values",
            Self::Atoms | Self::RingSize => "atoms",
            Self::GraphSlots => "slots",
            Self::Labels => "labels",
            Self::Bonds => "bonds",
            Self::Hydrogens => "hydrogens",
            Self::Components => "components",
            Self::ExpressionNodes | Self::DigraphNodes => "nodes",
            Self::ExpressionDepth
            | Self::RecursionDepth
            | Self::BranchDepth
            | Self::DigraphDepth => "levels",
            Self::RingClosures => "ring closures",
            Self::Work => "steps",
            Self::SearchStates => "states",
            Self::CandidateCycles => "cycles",
            Self::PathExpansions => "expansions",
            Self::ShortestPaths => "paths",
            Self::Matches => "matches",
            Self::Permutations => "permutations",
            Self::Structures => "structures",
            Self::LatticeImages => "images",
            Self::Positions => "positions",
            Self::MatrixEntries => "entries",
            Self::Residues => "residues",
            Self::Frames => "frames",
            Self::Fields => "fields",
            Self::Combinations => "combinations",
            Self::CandidatePairs => "pairs",
            Self::Ladders => "ladders",
        }
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

    /// The limit with the observed value and the bound, such as `work limit:
    /// 69,400,898 steps of 69,248,000`, for the error of an operation to name
    /// in a sentence of its own.
    ///
    /// ```
    /// use kekule::diagnostics::{LimitExceeded, Resource};
    ///
    /// let limit = LimitExceeded::new(Resource::Work, 69_248_000, 69_400_898);
    /// assert_eq!(
    ///     format!("Ring perception exceeded its {}.", limit.describe()),
    ///     "Ring perception exceeded its work limit: 69,400,898 steps of 69,248,000."
    /// );
    /// ```
    pub fn describe(&self) -> impl fmt::Display + '_ {
        DescribedLimit(self)
    }
}

struct DescribedLimit<'a>(&'a LimitExceeded);

impl fmt::Display for DescribedLimit<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let limit = self.0;
        let unit = limit.resource.unit();
        let unit = match (limit.observed, unit) {
            (1, "matches") => "match",
            (1, "entries") => "entry",
            (1, unit) => unit.strip_suffix('s').unwrap_or(unit),
            (_, unit) => unit,
        };
        write!(
            formatter,
            "{} limit: {} {unit} of {}",
            limit.resource.limit_name(),
            Count(limit.observed),
            Count(limit.limit)
        )
    }
}

/// Writes `Exceeded the work limit: 69,400,898 steps of 69,248,000.`; an
/// error that carries the limit names its operation with
/// [`LimitExceeded::describe`] instead.
impl fmt::Display for LimitExceeded {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "Exceeded the {}.", self.describe())
    }
}

/// A count written with thousands separators, such as `69,248,000`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Count(pub(crate) usize);

impl fmt::Display for Count {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let digits = self.0.to_string();
        let lead = digits.len() % 3;
        for (index, digit) in digits.char_indices() {
            if index != 0 && index % 3 == lead {
                formatter.write_str(",")?;
            }
            write!(formatter, "{digit}")?;
        }
        Ok(())
    }
}

/// A count with its noun, such as `1 atom` or `1,024 atoms`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Counted<'a> {
    count: usize,
    singular: &'a str,
    plural: &'a str,
}

/// `count` with `singular` or `plural`, such as `1 atom` or `1,024 atoms`.
pub(crate) const fn counted<'a>(count: usize, singular: &'a str, plural: &'a str) -> Counted<'a> {
    Counted {
        count,
        singular,
        plural,
    }
}

impl fmt::Display for Counted<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let noun = if self.count == 1 {
            self.singular
        } else {
            self.plural
        };
        write!(formatter, "{} {noun}", Count(self.count))
    }
}

/// Items joined as a list in a sentence, such as `"CG", "CD" and "NE"`.
pub(crate) fn listed<T: fmt::Display>(items: &[T]) -> String {
    match items {
        [] => String::new(),
        [only] => only.to_string(),
        [rest @ .., last] => {
            let rest = rest.iter().map(ToString::to_string).collect::<Vec<_>>();
            format!("{} and {last}", rest.join(", "))
        }
    }
}

/// Writes a value's display text with its first letter in upper case, for a
/// noun phrase such as an identifier (`atom 3`) that starts a sentence.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Capitalized<T>(pub(crate) T);

impl<T: fmt::Display> fmt::Display for Capitalized<T> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = self.0.to_string();
        let mut characters = text.chars();
        if let Some(first) = characters.next() {
            write!(formatter, "{}", first.to_uppercase())?;
            formatter.write_str(characters.as_str())?;
        }
        Ok(())
    }
}

/// An error followed by every error of its
/// [`source`](std::error::Error::source) chain, written as consecutive
/// sentences.
///
/// Every Kekule error displays as one or two sentences and does not repeat
/// its source, so the chain reads as one message, such as `Cannot build the
/// solvated model. Atom 12 overlaps atom 3.` Errors from outside Kekule in
/// the chain, such as I/O errors, get an upper-case first letter and a
/// closing full stop when they lack them. A cause whose text equals the
/// previous one, as a transparent wrapper writes it, appears once. The
/// alternate form, `{:#}`, writes one error per line.
///
/// ```
/// use kekule::diagnostics::ErrorChain;
///
/// let error = kekule::smiles::to_molecules("C1CC").unwrap_err();
/// let message = ErrorChain::new(&error).to_string();
/// assert!(message.starts_with("Cannot parse the SMILES input"));
/// assert!(message.ends_with('.'));
/// ```
#[derive(Clone, Copy)]
pub struct ErrorChain<'a> {
    error: &'a (dyn std::error::Error + 'static),
}

impl<'a> ErrorChain<'a> {
    /// The chain that starts at `error`.
    pub fn new(error: &'a (dyn std::error::Error + 'static)) -> Self {
        Self { error }
    }
}

impl fmt::Debug for ErrorChain<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ErrorChain")
            .field("error", &self.error)
            .finish()
    }
}

impl fmt::Display for ErrorChain<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let separator = if formatter.alternate() { "\n" } else { " " };
        let mut previous: Option<String> = None;
        let mut next = Some(self.error);
        while let Some(error) = next {
            next = error.source();
            let text = error.to_string();
            let text = text.trim();
            if text.is_empty() || previous.as_deref() == Some(text) {
                continue;
            }
            if previous.is_some() {
                formatter.write_str(separator)?;
            }
            write!(formatter, "{}", Capitalized(text))?;
            if !text.ends_with(['.', '!', '?']) {
                formatter.write_str(".")?;
            }
            previous = Some(text.to_owned());
        }
        Ok(())
    }
}

impl std::error::Error for LimitExceeded {}

/// A host's request to stop a cancellable operation partway.
///
/// It wraps a check that the operation polls between its stages and
/// periodically inside its long loops, a few times per millisecond of work at
/// most. Once the check returns `true`, the operation stops at its next poll
/// and fails with [`Cancelled`], returning no partial result. Operations
/// without a cancellation never poll. The check runs on the thread doing the
/// work, so it should be cheap, such as reading an atomic flag.
///
/// Clones share one check. Equality is identity: clones are equal, while two
/// separately constructed cancellations are not, whatever their checks do.
///
/// ```
/// use std::sync::atomic::{AtomicBool, Ordering};
/// use std::sync::Arc;
/// use kekule::diagnostics::Cancellation;
///
/// let flag = Arc::new(AtomicBool::new(false));
/// let cancellation = Cancellation::from_flag(Arc::clone(&flag));
/// assert!(!cancellation.is_cancelled());
/// flag.store(true, Ordering::Relaxed);
/// assert!(cancellation.is_cancelled());
/// ```
#[derive(Clone)]
pub struct Cancellation {
    check: Arc<dyn Fn() -> bool + Send + Sync>,
}

impl Cancellation {
    /// A cancellation requested once `check` returns `true`.
    pub fn new(check: impl Fn() -> bool + Send + Sync + 'static) -> Self {
        Self {
            check: Arc::new(check),
        }
    }

    /// A cancellation requested once `flag` is set.
    pub fn from_flag(flag: Arc<std::sync::atomic::AtomicBool>) -> Self {
        Self::new(move || flag.load(std::sync::atomic::Ordering::Relaxed))
    }

    /// Whether cancellation has been requested, running the check once.
    pub fn is_cancelled(&self) -> bool {
        (self.check)()
    }

    /// Fails with [`Cancelled`] once cancellation has been requested, running
    /// the check once; for cancellable operations outside this crate.
    ///
    /// # Errors
    ///
    /// Returns [`Cancelled`] when the check returns `true`.
    pub fn check(&self) -> Result<(), Cancelled> {
        if self.is_cancelled() {
            Err(Cancelled)
        } else {
            Ok(())
        }
    }

    /// Polls `cancellation`, which is `None` for an operation that cannot be
    /// cancelled.
    #[inline]
    pub(crate) fn poll(cancellation: Option<&Self>) -> Result<(), Cancelled> {
        match cancellation {
            Some(cancellation) if cancellation.is_cancelled() => Err(Cancelled),
            _ => Ok(()),
        }
    }

    /// Polls `cancellation` on every [`Self::POLL_STRIDE`]th step of a loop
    /// whose steps are too cheap to poll each.
    #[inline]
    pub(crate) fn poll_every(cancellation: Option<&Self>, step: usize) -> Result<(), Cancelled> {
        if cancellation.is_some() && step.is_multiple_of(Self::POLL_STRIDE) {
            Self::poll(cancellation)
        } else {
            Ok(())
        }
    }

    /// Steps between polls of a loop with cheap steps.
    pub(crate) const POLL_STRIDE: usize = 1024;
}

impl fmt::Debug for Cancellation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Cancellation")
    }
}

impl PartialEq for Cancellation {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.check, &other.check)
    }
}

impl Eq for Cancellation {}

// Kekule only calls the check and reads its answer, so a panic inside it
// cannot leave Kekule state broken, and options holding a cancellation stay
// unwind safe as they were without one.
impl std::panic::UnwindSafe for Cancellation {}
impl std::panic::RefUnwindSafe for Cancellation {}

/// An operation stopped because its [`Cancellation`] was requested.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub struct Cancelled;

impl fmt::Display for Cancelled {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("The operation was cancelled.")
    }
}

impl std::error::Error for Cancelled {}

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
    /// topologies, models and structure sets.
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
    Structure {
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
    /// The written structure, a model, a geometry, or a trajectory frame, has
    /// this state, which the format cannot carry, and the writer's policy
    /// for it is [`UnrepresentablePolicy::Fail`].
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
    /// A native structure set record stores its geometries elsewhere, so it has
    /// none to give; the kind is [`DiagnosticKind::Inconsistent`].
    NativeExternalGeometries {
        /// The number of geometries the record describes.
        geometries: usize,
    },
}

/// Structure, geometry, or trajectory frame state named by
/// [`FormatErrorDetail::UnrepresentableState`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StructureField {
    /// A geometry's relative weight in its set's distribution.
    Weight,
    /// The energy of a geometry's state.
    Energy,
    /// The state a geometry belongs to: the grouping of geometries into states.
    State,
    /// The temperature a distribution's populations describe.
    Temperature,
    /// The properties of a structure set's state.
    StateProperties,
    /// A structure set's state representative that is not the state's first
    /// member.
    Representative,
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
    /// A structure set's geometry that is outside its distribution, which a file
    /// whose every structure is a member cannot tell apart.
    Membership,
}

impl fmt::Display for StructureField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Weight => "weight",
            Self::Energy => "state energy",
            Self::State => "state",
            Self::Temperature => "temperature",
            Self::StateProperties => "state properties",
            Self::Representative => "state representative",
            Self::Order => "geometry order",
            Self::Time => "time",
            Self::Step => "MD step",
            Self::Velocities => "velocities",
            Self::Forces => "forces",
            Self::AtomProperties => "per-atom structure properties",
            Self::BondProperties => "per-bond structure properties",
            Self::Properties => "structure properties",
            Self::SetProperties => "structure set properties",
            Self::Occupancies => "occupancies",
            Self::BFactors => "B factors",
            Self::OrderStep => "order step",
            Self::OrderCoordinate => "order coordinate",
            Self::Membership => "distribution membership",
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
        .map_err(|error| FormatError::parse(format, DiagnosticKind::Io, "").with_source(error))?;
    if bytes.len() > max_bytes {
        let limit = LimitExceeded::new(Resource::InputBytes, max_bytes, bytes.len());
        return Err(FormatError::parse(
            format,
            DiagnosticKind::Limit,
            format!("The input exceeds its {}.", limit.describe()),
        )
        .exceeding(limit));
    }
    String::from_utf8(bytes).map_err(|error| {
        FormatError::parse(format, DiagnosticKind::Io, "The input is not valid UTF-8.")
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
    fn limits_report_resource_bound_and_first_excess() {
        assert_eq!(LimitExceeded::check(Resource::Atoms, 3, 3), Ok(()));
        let exceeded = LimitExceeded::check(Resource::Atoms, 3, 5).unwrap_err();
        assert_eq!(
            (exceeded.resource(), exceeded.limit(), exceeded.observed()),
            (Resource::Atoms, 3, 5)
        );
        assert_eq!(
            exceeded.to_string(),
            "Exceeded the atom limit: 5 atoms of 3."
        );
        let work = LimitExceeded::new(Resource::Work, 69_248_000, 69_400_898);
        assert_eq!(
            format!("Ring perception exceeded its {}.", work.describe()),
            "Ring perception exceeded its work limit: 69,400,898 steps of 69,248,000."
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
    fn counted_nouns_agree_with_their_count() {
        assert_eq!(counted(1, "atom", "atoms").to_string(), "1 atom");
        assert_eq!(counted(0, "atom", "atoms").to_string(), "0 atoms");
        assert_eq!(
            counted(1_024, "entry", "entries").to_string(),
            "1,024 entries"
        );
        assert_eq!(
            LimitExceeded::first_beyond(Resource::Rows, 0)
                .describe()
                .to_string(),
            "row limit: 1 row of 0"
        );
        assert_eq!(
            LimitExceeded::first_beyond(Resource::MatrixEntries, 0)
                .describe()
                .to_string(),
            "matrix entry limit: 1 entry of 0"
        );
    }

    #[test]
    fn lists_join_with_and() {
        assert_eq!(listed::<&str>(&[]), "");
        assert_eq!(listed(&["a"]), "a");
        assert_eq!(listed(&["a", "b"]), "a and b");
        assert_eq!(listed(&["a", "b", "c"]), "a, b and c");
    }

    #[test]
    fn counts_group_thousands() {
        for (count, text) in [
            (0, "0"),
            (999, "999"),
            (1_000, "1,000"),
            (12_345, "12,345"),
            (123_456, "123,456"),
            (1_234_567, "1,234,567"),
        ] {
            assert_eq!(Count(count).to_string(), text);
        }
    }

    #[test]
    fn error_chains_read_as_sentences() {
        #[derive(Debug)]
        struct Outer(std::io::Error);
        impl fmt::Display for Outer {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("Cannot build the solvated model.")
            }
        }
        impl std::error::Error for Outer {
            fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
                Some(&self.0)
            }
        }

        let error = Outer(std::io::Error::other("disk full"));
        assert_eq!(
            ErrorChain::new(&error).to_string(),
            "Cannot build the solvated model. Disk full."
        );
        assert_eq!(
            format!("{:#}", ErrorChain::new(&error)),
            "Cannot build the solvated model.
Disk full."
        );

        let sink = FormatError::sink(
            Format::Pdb,
            std::io::Error::other("Access is denied. (os error 5)"),
        );
        assert_eq!(
            ErrorChain::new(&sink).to_string(),
            "Cannot write the PDB output. Access is denied. (os error 5)."
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
