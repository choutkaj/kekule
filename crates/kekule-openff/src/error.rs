//! One error type with a stable kind, readable detail, and source chain.
use std::{
    fmt,
    path::{Path, PathBuf},
    sync::Arc,
};

use kekule::topology::MoleculeDefinitionId;

/// The category of a parameterization failure, for programmatic handling.
///
/// Every category is a deliberate, explicit failure: the crate never
/// substitutes a fallback model, truncates a search, or publishes a partial
/// result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// Reading an input file failed; see [`Error::path`].
    Io,
    /// OFFXML is malformed or outside the supported SMIRNOFF subset.
    ForceField,
    /// Two force fields cannot be composed with [`crate::ForceField::append`].
    Composition,
    /// A NAGL model bundle is malformed, corrupt, or unsupported.
    Model,
    /// The supplied model is not the one the force field declares.
    ModelMismatch,
    /// The molecule is outside the supported chemical domain, for example it
    /// has implicit hydrogens, radicals, or elements the model excludes.
    UnsupportedMolecule,
    /// The force field's rules leave an atom or interaction without
    /// parameters.
    Unparameterized,
    /// Charges cannot be assigned: library charges are incomplete or do not
    /// conserve the formal charge, no charge method applies, or a model
    /// produced invalid values.
    Charges,
    /// The fixed-H InChI lookup identifier cannot be computed.
    Identity,
    /// A bounded combinatorial search exceeded its limit.
    ResourceLimit,
    /// An underlying Kekule perception, matching, or editing operation failed;
    /// see [`std::error::Error::source`].
    Chemistry,
}

/// A parameterization, loading, or composition failure.
///
/// [`Self::kind`] classifies the failure; the display text adds the specific
/// detail and context. Failures inside one molecule of a system name its
/// definition through [`Self::definition`], and file failures name the file
/// through [`Self::path`].
#[derive(Debug, Clone)]
pub struct Error {
    kind: ErrorKind,
    detail: String,
    path: Option<PathBuf>,
    definition: Option<MoleculeDefinitionId>,
    source: Option<Arc<dyn std::error::Error + Send + Sync>>,
}

impl Error {
    pub(crate) fn new(kind: ErrorKind, detail: impl fmt::Display) -> Self {
        Self {
            kind,
            detail: detail.to_string(),
            path: None,
            definition: None,
            source: None,
        }
    }

    /// Wraps an underlying error, keeping it as the source.
    pub(crate) fn wrap(
        kind: ErrorKind,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> Self {
        Self {
            detail: source.to_string(),
            source: Some(Arc::new(source)),
            ..Self::new(kind, "")
        }
    }

    pub(crate) fn chemistry(source: impl std::error::Error + Send + Sync + 'static) -> Self {
        Self::wrap(ErrorKind::Chemistry, source)
    }

    pub(crate) fn io(path: &Path, source: std::io::Error) -> Self {
        Self::wrap(ErrorKind::Io, source).at_path(path)
    }

    /// Adds the file this failure concerns, keeping a more specific path.
    pub(crate) fn at_path(mut self, path: &Path) -> Self {
        self.path.get_or_insert_with(|| path.to_owned());
        self
    }

    /// Adds the molecule definition this failure concerns.
    pub(crate) fn in_definition(mut self, definition: MoleculeDefinitionId) -> Self {
        self.definition.get_or_insert(definition);
        self
    }

    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// The specific failure, without path or definition context.
    pub fn detail(&self) -> &str {
        &self.detail
    }

    /// The file being read, for loading failures.
    pub fn path(&self) -> Option<&Path> {
        self.path.as_deref()
    }

    /// The molecule definition being parameterized, for failures inside one
    /// molecule of a system.
    pub fn definition(&self) -> Option<MoleculeDefinitionId> {
        self.definition
    }

    /// The wrapped lower-level error, such as a Kekule perception or I/O
    /// error, for downcasting.
    ///
    /// As with [`std::io::Error`], the wrapper is transparent: its display text
    /// is the wrapped error's, and [`std::error::Error::source`] continues with
    /// the wrapped error's own source, so a source chain never repeats a
    /// message.
    pub fn get_ref(&self) -> Option<&(dyn std::error::Error + Send + Sync + 'static)> {
        self.source.as_deref()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(path) = &self.path {
            write!(f, "{}: ", path.display())?;
        }
        if let Some(definition) = self.definition {
            write!(f, "molecule definition {}: ", definition.index())?;
        }
        f.write_str(&self.detail)
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source.as_deref().and_then(std::error::Error::source)
    }
}

pub(crate) type Result<T> = std::result::Result<T, Error>;
