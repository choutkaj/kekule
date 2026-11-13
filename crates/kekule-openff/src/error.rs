//! One error type with a stable kind, readable detail, and source chain.
use std::{fmt, path::Path, sync::Arc};

use kekule::diagnostics::LimitExceeded;
use kekule::molecule::AtomId;
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
    /// A lookup hit cannot be mapped onto its entry.
    Identity,
    /// A bounded combinatorial search exceeded its limit; see
    /// [`Error::limit_exceeded`].
    ResourceLimit,
    /// An underlying Kekule perception, matching, or editing operation failed;
    /// see [`std::error::Error::source`].
    Chemistry,
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Io => "reading a file failed",
            Self::ForceField => "the force field is invalid or unsupported",
            Self::Composition => "the force fields cannot be composed",
            Self::Model => "the NAGL model bundle is invalid or unsupported",
            Self::ModelMismatch => "the model is not the one the force field declares",
            Self::UnsupportedMolecule => "the molecule is outside the supported domain",
            Self::Unparameterized => "the force field leaves the molecule unparameterized",
            Self::Charges => "charges cannot be assigned",
            Self::Identity => "a lookup hit cannot be mapped onto its entry",
            Self::ResourceLimit => "a bounded search exceeded its limit",
            Self::Chemistry => "a Kekule chemistry operation failed",
        })
    }
}

/// A parameterization, loading, or composition failure.
///
/// [`Self::kind`] classifies the failure. The display text is one or two
/// sentences: the file or molecule definition that failed, such as `Cannot
/// parameterize molecule definition 2.`, then the [detail](Self::detail),
/// such as `Atoms 3 and 5 have no vdW parameter.` Failures inside one
/// molecule of a system name its definition through [`Self::definition`],
/// failures caused by particular atoms name all of them through
/// [`Self::atoms`], and file failures name the file through [`Self::path`]. A
/// failure caused by a lower-level error, such as a Kekule perception or I/O
/// error, returns it from [`source`](std::error::Error::source) and does not
/// repeat its message; `kekule::diagnostics::ErrorChain` writes both.
#[derive(Debug, Clone)]
pub struct Error {
    kind: ErrorKind,
    detail: String,
    path: Option<Box<Path>>,
    definition: Option<MoleculeDefinitionId>,
    atoms: Box<[AtomId]>,
    limit: Option<LimitExceeded>,
    source: Option<Arc<dyn std::error::Error + Send + Sync>>,
}

impl Error {
    pub(crate) fn new(kind: ErrorKind, detail: impl fmt::Display) -> Self {
        Self {
            kind,
            detail: detail.to_string(),
            path: None,
            definition: None,
            atoms: Box::default(),
            limit: None,
            source: None,
        }
    }

    /// Wraps an underlying error, keeping it as the source; the detail says
    /// what failed at this level, and the source says why.
    pub(crate) fn caused(
        kind: ErrorKind,
        detail: impl fmt::Display,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> Self {
        Self::new(kind, detail).with_source(source)
    }

    /// A [`ErrorKind::ResourceLimit`] failure with the exceeded bound.
    pub(crate) fn limit(detail: impl fmt::Display, limit: LimitExceeded) -> Self {
        Self {
            limit: Some(limit),
            ..Self::new(ErrorKind::ResourceLimit, detail)
        }
    }

    /// Adds the exceeded bound of a [`ErrorKind::ResourceLimit`] failure.
    pub(crate) fn with_limit(mut self, limit: LimitExceeded) -> Self {
        self.limit = Some(limit);
        self
    }

    /// Adds the underlying error that caused this failure.
    pub(crate) fn with_source(
        mut self,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> Self {
        self.source = Some(Arc::new(source));
        self
    }

    /// A failed Kekule operation; the detail says what failed.
    pub(crate) fn chemistry(
        detail: impl fmt::Display,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> Self {
        Self::caused(ErrorKind::Chemistry, detail, source)
    }

    /// A file that cannot be read; the path names it and the source says why.
    pub(crate) fn io(path: &Path, source: std::io::Error) -> Self {
        Self::new(ErrorKind::Io, "")
            .with_source(source)
            .at_path(path)
    }

    /// Adds the file this failure concerns, keeping a more specific path.
    pub(crate) fn at_path(mut self, path: &Path) -> Self {
        self.path.get_or_insert_with(|| path.into());
        self
    }

    /// Adds the molecule definition this failure concerns.
    pub(crate) fn in_definition(mut self, definition: MoleculeDefinitionId) -> Self {
        self.definition.get_or_insert(definition);
        self
    }

    /// Adds the atoms this failure concerns, sorted and without repeats.
    pub(crate) fn with_atoms(mut self, atoms: impl IntoIterator<Item = AtomId>) -> Self {
        let mut all = std::mem::take(&mut self.atoms).into_vec();
        all.extend(atoms);
        all.sort_unstable();
        all.dedup();
        self.atoms = all.into_boxed_slice();
        self
    }

    /// The failure category.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// The specific failure in one or two sentences, without path or
    /// definition context or the message of the
    /// [`source`](std::error::Error::source); empty when the path and the
    /// source say it all.
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

    /// Every atom the failure concerns, in ascending order: for example all
    /// atoms of elements outside a charge model's domain, every atom of the
    /// interactions a force field leaves unassigned, or the atoms with implicit
    /// hydrogens. Empty for failures that concern no particular atom, such as
    /// file, force-field, model, or composition failures.
    ///
    /// The IDs belong to the molecule of [`Self::definition`] when a system's
    /// definition failed, so atom `atom` of instance `instance` is
    /// `InstanceAtomId::new(instance, atom)`, and otherwise to the molecule
    /// passed in, for example to [`crate::ForceField::parameterize_molecule`] or
    /// [`crate::NaglModel::assign_charges`].
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    /// The exceeded bound, for [`ErrorKind::ResourceLimit`] failures whose
    /// bound is a single count.
    pub fn limit_exceeded(&self) -> Option<LimitExceeded> {
        self.limit
    }

    /// The wrapped lower-level error, such as a Kekule perception or I/O
    /// error, for downcasting; the same error
    /// [`source`](std::error::Error::source) returns.
    pub fn get_ref(&self) -> Option<&(dyn std::error::Error + Send + Sync + 'static)> {
        self.source.as_deref()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut context = true;
        if let Some(path) = &self.path {
            let path = path.display();
            match self.kind {
                ErrorKind::Io => write!(f, "Cannot read \"{path}\".")?,
                ErrorKind::ForceField => write!(f, "Cannot load the force field \"{path}\".")?,
                ErrorKind::Model => write!(f, "Cannot load the NAGL model \"{path}\".")?,
                _ => write!(f, "Cannot load \"{path}\".")?,
            }
        } else if let Some(definition) = self.definition {
            write!(
                f,
                "Cannot parameterize molecule definition {}.",
                definition.index()
            )?;
        } else {
            context = false;
        }
        if !self.detail.is_empty() {
            if context {
                f.write_str(" ")?;
            }
            f.write_str(&self.detail)?;
        }
        Ok(())
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_deref()
            .map(|source| source as &(dyn std::error::Error + 'static))
    }
}

/// Names atoms in a failure's detail, such as `atoms 3, 5 and 7`.
pub(crate) fn atom_list(atoms: &[AtomId]) -> String {
    let ids: Vec<_> = atoms.iter().map(|atom| atom.index().to_string()).collect();
    match ids.as_slice() {
        [] => "no atoms".to_owned(),
        [one] => format!("atom {one}"),
        [rest @ .., last] => format!("atoms {} and {last}", rest.join(", ")),
    }
}

/// `a` or `an` for a noun in a message.
pub(crate) fn article(noun: &str) -> &'static str {
    if noun.starts_with(['a', 'e', 'i', 'o', 'u']) {
        "an"
    } else {
        "a"
    }
}

/// A count with its noun, such as `1 atom` or `3 atoms`.
pub(crate) fn counted(count: usize, singular: &str, plural: &str) -> String {
    if count == 1 {
        format!("1 {singular}")
    } else {
        format!("{count} {plural}")
    }
}

/// `text` with an upper-case first letter, for a phrase such as
/// [`atom_list`]'s that starts a sentence.
pub(crate) fn capitalized(text: &str) -> String {
    let mut characters = text.chars();
    match characters.next() {
        Some(first) => first.to_uppercase().chain(characters).collect(),
        None => String::new(),
    }
}

pub(crate) type Result<T> = std::result::Result<T, Error>;
