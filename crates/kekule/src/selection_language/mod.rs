//! The selection language: text expressions that evaluate to atom
//! selections. The grammar is documented in `docs/selection-language.md` and
//! on [`crate::select`].

mod ast;
mod eval;
mod grid;
mod print;
mod syntax;

use std::fmt;
use std::ops::Range;
use std::sync::Arc;

use crate::diagnostics::{DiagnosticKind, Format, FormatError};
use crate::properties::PropertyKey;
use crate::structure::{AsModelView, ModelView, Positions};
use crate::substructure::{SubstructureMatchError, SubstructureMatchOptions};
use crate::topology::{AtomSelection, Topology};
use crate::units::UnitError;

pub use ast::HostTermKind;

/// Default [`SelectParseOptions::max_input_bytes`]: 64 KiB.
pub const DEFAULT_SELECT_MAX_INPUT_BYTES: usize = 64 * 1024;
/// Default [`SelectParseOptions::max_depth`].
pub const DEFAULT_SELECT_MAX_DEPTH: usize = 256;
/// Default maximum SMARTS matches a `smarts` term visits; see
/// [`SelectEvaluateOptions::with_smarts`].
pub const DEFAULT_SELECT_MAX_SMARTS_MATCHES: usize = 1_000_000;

/// A parsed selection expression.
///
/// [`Display`](fmt::Display) prints the canonical text: parsing it gives an
/// equal expression. Equality compares structure, not source spelling, so
/// `elem C` equals `element C`.
#[derive(Debug, Clone, PartialEq)]
pub struct SelectionExpr {
    root: ast::Node,
}

impl fmt::Display for SelectionExpr {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        print::node(formatter, &self.root, 0)
    }
}

/// Displays one node canonically, for naming a failing term.
struct NodeDisplay<'a>(&'a ast::Node);

impl fmt::Display for NodeDisplay<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        print::node(formatter, self.0, 0)
    }
}

impl SelectionExpr {
    pub(crate) fn parse(text: &str, options: &SelectParseOptions) -> Result<Self, FormatError> {
        syntax::parse(text, options).map(|root| Self { root })
    }

    /// Evaluates the expression against a topology alone.
    ///
    /// # Errors
    ///
    /// [`SelectEvaluateError::RequiresStructure`] for terms that read
    /// coordinates or site values (`within`, `bfactor`, `occupancy`), and the
    /// other [`SelectEvaluateError`] cases.
    pub fn evaluate(
        &self,
        topology: impl Into<Arc<Topology>>,
        options: SelectEvaluateOptions<'_>,
    ) -> Result<AtomSelection, SelectEvaluateError> {
        let topology = topology.into();
        eval::Evaluator::new(self, &topology, None, &options).run()
    }

    /// Evaluates the expression against one structure: a model, an ensemble
    /// member, a trajectory frame, or any other [`AsModelView`]. The result
    /// binds to the structure's topology.
    ///
    /// # Errors
    ///
    /// The [`SelectEvaluateError`] cases other than
    /// [`RequiresStructure`](SelectEvaluateError::RequiresStructure).
    pub fn evaluate_structure(
        &self,
        structure: &impl AsModelView,
        options: SelectEvaluateOptions<'_>,
    ) -> Result<AtomSelection, SelectEvaluateError> {
        let view = structure.as_model_view();
        let topology = view.shared_topology();
        eval::Evaluator::new(self, &topology, Some(view), &options).run()
    }

    /// The host terms of the expression in text order, for applications that
    /// check references before evaluating.
    pub fn host_terms(&self) -> Vec<HostTerm<'_>> {
        fn visit<'a>(node: &'a ast::Node, terms: &mut Vec<HostTerm<'a>>) {
            match &node.kind {
                ast::NodeKind::Host(kind, text) => terms.push(HostTerm { kind: *kind, text }),
                ast::NodeKind::And(operands) | ast::NodeKind::Or(operands) => {
                    for operand in operands {
                        visit(operand, terms);
                    }
                }
                ast::NodeKind::Not(inner)
                | ast::NodeKind::Bonded(inner)
                | ast::NodeKind::Expand { inner, .. }
                | ast::NodeKind::Within {
                    reference: inner, ..
                } => visit(inner, terms),
                _ => {}
            }
        }
        let mut terms = Vec::new();
        visit(&self.root, &mut terms);
        terms
    }
}

/// Parse options and resource limits of the selection language.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectParseOptions {
    max_input_bytes: usize,
    max_depth: usize,
    host_keywords: Vec<String>,
}

impl Default for SelectParseOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: DEFAULT_SELECT_MAX_INPUT_BYTES,
            max_depth: DEFAULT_SELECT_MAX_DEPTH,
            host_keywords: Vec::new(),
        }
    }
}

fn option_error(message: impl Into<String>) -> FormatError {
    FormatError::parse(Format::Selection, DiagnosticKind::InvalidValue, message)
}

impl SelectParseOptions {
    /// Maximum UTF-8 byte length of the text.
    pub const fn max_input_bytes(&self) -> usize {
        self.max_input_bytes
    }

    /// Sets [`Self::max_input_bytes`].
    ///
    /// # Errors
    ///
    /// [`DiagnosticKind::InvalidValue`] when the bound is zero.
    pub fn with_max_input_bytes(mut self, max_input_bytes: usize) -> Result<Self, FormatError> {
        if max_input_bytes == 0 {
            return Err(option_error("max_input_bytes must be greater than zero"));
        }
        self.max_input_bytes = max_input_bytes;
        Ok(self)
    }

    /// Maximum nesting of parentheses and prefix operators.
    pub const fn max_depth(&self) -> usize {
        self.max_depth
    }

    /// Sets [`Self::max_depth`].
    ///
    /// # Errors
    ///
    /// [`DiagnosticKind::InvalidValue`] when the bound is zero.
    pub fn with_max_depth(mut self, max_depth: usize) -> Result<Self, FormatError> {
        if max_depth == 0 {
            return Err(option_error("max_depth must be greater than zero"));
        }
        self.max_depth = max_depth;
        Ok(self)
    }

    /// Bare keywords that parse as [`HostTermKind::Keyword`] host terms.
    pub fn host_keywords(&self) -> &[String] {
        &self.host_keywords
    }

    /// Sets [`Self::host_keywords`], such as `displayed` and `selected`.
    ///
    /// # Errors
    ///
    /// [`DiagnosticKind::InvalidValue`] when a keyword is not an ASCII
    /// identifier (a letter or `_`, then letters, digits, or `_`) or is a
    /// built-in keyword.
    pub fn with_host_keywords(
        mut self,
        keywords: impl IntoIterator<Item = impl Into<String>>,
    ) -> Result<Self, FormatError> {
        let keywords = keywords.into_iter().map(Into::into).collect::<Vec<_>>();
        for keyword in &keywords {
            let mut characters = keyword.chars();
            let identifier = characters
                .next()
                .is_some_and(|first| first.is_ascii_alphabetic() || first == '_')
                && characters.all(|next| next.is_ascii_alphanumeric() || next == '_');
            if !identifier || syntax::KEYWORDS.contains(&keyword.as_str()) {
                return Err(option_error(format!(
                    "`{keyword}` cannot be a host keyword"
                )));
            }
        }
        self.host_keywords = keywords;
        Ok(self)
    }

    fn is_host_keyword(&self, word: &str) -> bool {
        self.host_keywords.iter().any(|keyword| keyword == word)
    }
}

/// Evaluation options: the host resolver and SMARTS matching bounds.
#[derive(Clone, Copy)]
pub struct SelectEvaluateOptions<'a> {
    resolver: Option<&'a dyn SelectionResolver>,
    smarts: SubstructureMatchOptions,
}

impl fmt::Debug for SelectEvaluateOptions<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SelectEvaluateOptions")
            .field("resolver", &self.resolver.is_some())
            .field("smarts", &self.smarts)
            .finish()
    }
}

impl Default for SelectEvaluateOptions<'_> {
    fn default() -> Self {
        Self {
            resolver: None,
            smarts: SubstructureMatchOptions::default()
                .with_max_matches(DEFAULT_SELECT_MAX_SMARTS_MATCHES)
                .expect("the default SMARTS match bound is positive"),
        }
    }
}

impl<'a> SelectEvaluateOptions<'a> {
    /// The resolver of host terms, if any.
    pub fn resolver(&self) -> Option<&'a dyn SelectionResolver> {
        self.resolver
    }

    /// Sets [`Self::resolver`].
    #[must_use]
    pub fn with_resolver(mut self, resolver: &'a dyn SelectionResolver) -> Self {
        self.resolver = Some(resolver);
        self
    }

    /// Bounds of `smarts` terms. The default allows
    /// [`DEFAULT_SELECT_MAX_SMARTS_MATCHES`] matches; exceeding a bound is an
    /// error, never a partial selection.
    pub const fn smarts(&self) -> SubstructureMatchOptions {
        self.smarts
    }

    /// Sets [`Self::smarts`].
    #[must_use]
    pub const fn with_smarts(mut self, smarts: SubstructureMatchOptions) -> Self {
        self.smarts = smarts;
        self
    }
}

/// One host term of an expression: `#3` is a [`HostTermKind::Reference`]
/// with text `3`, `$site` a [`HostTermKind::Named`] with text `site`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HostTerm<'a> {
    kind: HostTermKind,
    text: &'a str,
}

impl<'a> HostTerm<'a> {
    /// The term's syntax.
    pub const fn kind(&self) -> HostTermKind {
        self.kind
    }

    /// The text after the sigil, or the keyword.
    pub const fn text(&self) -> &'a str {
        self.text
    }
}

/// Where a host term appears, and so what it may resolve to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum HostRole {
    /// A term selecting atoms of the evaluated topology.
    Atoms,
    /// The reference of `within`: atoms of the evaluated topology, or points
    /// such as another object's coordinates.
    Reference,
}

/// What a resolver is asked to resolve against.
#[derive(Debug, Clone, Copy)]
pub struct HostRequest<'a> {
    topology: &'a Arc<Topology>,
    structure: Option<ModelView<'a>>,
    role: HostRole,
}

impl<'a> HostRequest<'a> {
    /// The evaluated topology; returned atoms must share its layout.
    pub fn topology(&self) -> &'a Topology {
        self.topology
    }

    /// The evaluated topology, shared.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    /// The evaluated structure, when evaluation has one.
    pub const fn structure(&self) -> Option<ModelView<'a>> {
        self.structure
    }

    /// Where the term appears.
    pub const fn role(&self) -> HostRole {
        self.role
    }
}

/// A resolved host term.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum HostValue {
    /// Atoms of the evaluated topology.
    Atoms(AtomSelection),
    /// Points, accepted only for [`HostRole::Reference`].
    Points(Positions),
}

/// A host's failure to resolve a term.
#[derive(Debug, Clone)]
pub struct HostTermError {
    message: String,
    source: Option<Arc<dyn std::error::Error + Send + Sync>>,
}

impl HostTermError {
    /// A failure described by `message`.
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            source: None,
        }
    }

    /// Attaches the underlying cause.
    #[must_use]
    pub fn with_source(mut self, source: impl std::error::Error + Send + Sync + 'static) -> Self {
        self.source = Some(Arc::new(source));
        self
    }

    /// The host's description of the failure.
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for HostTermError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl std::error::Error for HostTermError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_deref()
            .map(|source| source as &(dyn std::error::Error + 'static))
    }
}

/// Resolves host terms: `#...` references, `$name` named selections, and
/// registered host keywords.
///
/// This is an extension point. An implementation must be deterministic for
/// one request, must return [`HostValue::Atoms`] bound to a snapshot sharing
/// the request's topology layout, and may return [`HostValue::Points`] only
/// when the request's role is [`HostRole::Reference`]. Points are Cartesian,
/// in the evaluated structure's frame; no imaging is applied.
pub trait SelectionResolver {
    /// Resolves one term.
    ///
    /// # Errors
    ///
    /// A [`HostTermError`] when the host does not know the term or cannot
    /// resolve it for this request.
    fn resolve(
        &self,
        term: &HostTerm<'_>,
        request: HostRequest<'_>,
    ) -> Result<HostValue, HostTermError>;
}

/// A failure to evaluate a selection expression.
///
/// Every case names the failing term in canonical form and its byte span in
/// the parsed text.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum SelectEvaluateError {
    /// The term reads coordinates or site values, and evaluation has only a
    /// topology.
    RequiresStructure {
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
    },
    /// The term reads perception that is not installed.
    MissingPerception {
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
    },
    /// No atom property scope holds the key.
    UnknownProperty {
        /// The key.
        key: PropertyKey,
        /// The term's span.
        span: Range<usize>,
    },
    /// More than one atom property scope holds the key.
    AmbiguousProperty {
        /// The key.
        key: PropertyKey,
        /// The term's span.
        span: Range<usize>,
    },
    /// A value does not apply to the property column's type or unit.
    PropertyType {
        /// The key.
        key: PropertyKey,
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
    },
    /// A unit conversion failed.
    Unit {
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
        /// The conversion failure.
        source: Box<UnitError>,
    },
    /// The SMARTS pattern did not parse.
    InvalidSmarts {
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
        /// The parse failure.
        source: FormatError,
    },
    /// SMARTS matching failed, for example on missing perception or an
    /// exceeded bound.
    Smarts {
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
        /// The matching failure.
        source: Box<SubstructureMatchError>,
    },
    /// The expression has a host term and the options have no resolver.
    NoResolver {
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
    },
    /// The resolver failed.
    Host {
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
        /// The resolver's failure.
        source: HostTermError,
    },
    /// The resolver returned atoms of a topology that does not share the
    /// evaluated layout.
    HostTopologyMismatch {
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
    },
    /// The resolver returned points where atoms are required.
    HostPoints {
        /// The term.
        term: String,
        /// Its span.
        span: Range<usize>,
    },
}

impl SelectEvaluateError {
    /// The byte span of the failing term in the parsed text.
    pub fn span(&self) -> Range<usize> {
        match self {
            Self::RequiresStructure { span, .. }
            | Self::MissingPerception { span, .. }
            | Self::UnknownProperty { span, .. }
            | Self::AmbiguousProperty { span, .. }
            | Self::PropertyType { span, .. }
            | Self::Unit { span, .. }
            | Self::InvalidSmarts { span, .. }
            | Self::Smarts { span, .. }
            | Self::NoResolver { span, .. }
            | Self::Host { span, .. }
            | Self::HostTopologyMismatch { span, .. }
            | Self::HostPoints { span, .. } => span.clone(),
        }
    }
}

impl fmt::Display for SelectEvaluateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RequiresStructure { term, .. } => {
                write!(
                    formatter,
                    "`{term}` needs coordinates; evaluate it on a structure"
                )
            }
            Self::MissingPerception { term, .. } => write!(
                formatter,
                "`{term}` needs installed perception; perceive the topology first"
            ),
            Self::UnknownProperty { key, .. } => {
                write!(formatter, "no atom property `{key}` is defined")
            }
            Self::AmbiguousProperty { key, .. } => write!(
                formatter,
                "atom property `{key}` is defined in more than one scope"
            ),
            Self::PropertyType { key, term, .. } => write!(
                formatter,
                "`{term}` compares a value that does not fit property `{key}`"
            ),
            Self::Unit { term, .. } => write!(formatter, "`{term}` uses an incompatible unit"),
            Self::InvalidSmarts { term, .. } => write!(formatter, "`{term}` is not valid SMARTS"),
            Self::Smarts { term, .. } => write!(formatter, "matching `{term}` failed"),
            Self::NoResolver { term, .. } => {
                write!(formatter, "`{term}` is a host term and no resolver is set")
            }
            Self::Host { term, .. } => write!(formatter, "resolving `{term}` failed"),
            Self::HostTopologyMismatch { term, .. } => write!(
                formatter,
                "`{term}` resolved to atoms of another topology layout"
            ),
            Self::HostPoints { term, .. } => write!(
                formatter,
                "`{term}` resolved to points; points are only a `within` reference"
            ),
        }
    }
}

impl std::error::Error for SelectEvaluateError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit { source, .. } => Some(source.as_ref()),
            Self::InvalidSmarts { source, .. } => Some(source),
            Self::Smarts { source, .. } => Some(source.as_ref()),
            Self::Host { source, .. } => Some(source),
            _ => None,
        }
    }
}
