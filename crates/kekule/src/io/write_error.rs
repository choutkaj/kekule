//! Failure shared by the molecule-encoding writers (SMILES and Molfile CTABs).
use std::fmt;

/// Why a SMILES or Molfile writer produced no output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MoleculeWriteErrorKind {
    /// The format, or the requested format version, cannot encode valid
    /// source content without losing or changing it.
    UnsupportedRepresentation,
    /// The source lacks perception the writer requires; perceive it and retry.
    MissingPerception,
    /// A bounded writer calculation could not complete within its work limit.
    ResourceLimit,
    /// The source is internally inconsistent or incomplete, for example a
    /// model without positions for every atom.
    InvalidModel,
    /// The output sink failed.
    Io(std::io::ErrorKind),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MoleculeWriteError {
    kind: MoleculeWriteErrorKind,
    message: String,
}

impl MoleculeWriteError {
    pub(crate) fn unsupported(message: impl Into<String>) -> Self {
        Self::with_kind(MoleculeWriteErrorKind::UnsupportedRepresentation, message)
    }

    pub(crate) fn missing_perception(message: impl Into<String>) -> Self {
        Self::with_kind(MoleculeWriteErrorKind::MissingPerception, message)
    }

    pub(crate) fn resource_limit(message: impl Into<String>) -> Self {
        Self::with_kind(MoleculeWriteErrorKind::ResourceLimit, message)
    }

    pub(crate) fn invalid_model(message: impl Into<String>) -> Self {
        Self::with_kind(MoleculeWriteErrorKind::InvalidModel, message)
    }

    pub(crate) fn io(error: std::io::Error) -> Self {
        Self::with_kind(MoleculeWriteErrorKind::Io(error.kind()), error.to_string())
    }

    fn with_kind(kind: MoleculeWriteErrorKind, message: impl Into<String>) -> Self {
        Self {
            kind,
            message: message.into(),
        }
    }

    pub const fn kind(&self) -> MoleculeWriteErrorKind {
        self.kind
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for MoleculeWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for MoleculeWriteError {}
