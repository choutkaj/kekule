//! Failures of the molecule-encoding writers (SMILES and Molfile CTABs).
use crate::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, LimitExceeded, StructureField,
};
use crate::structure::StructureKind;

/// Constructors of write-stage [`FormatError`]s for one format, following the
/// writer failure taxonomy: cannot encode, missing perception, limit,
/// inconsistent source, and sink.
#[derive(Debug, Clone, Copy)]
pub(crate) struct WriteErrors(Format);

pub(crate) const SMILES_WRITE: WriteErrors = WriteErrors(Format::Smiles);
pub(crate) const MOLFILE_WRITE: WriteErrors = WriteErrors(Format::Molfile);

impl WriteErrors {
    /// The format, or the requested format version, cannot encode valid
    /// source content without losing or changing it.
    pub(crate) fn unsupported(self, message: impl Into<String>) -> FormatError {
        FormatError::write(self.0, DiagnosticKind::Unsupported, message)
    }

    /// The source lacks perception the writer requires.
    pub(crate) fn missing_perception(self, message: impl Into<String>) -> FormatError {
        FormatError::write(self.0, DiagnosticKind::MissingPerception, message)
    }

    /// A bounded writer calculation, described by `what`, could not complete
    /// within `limit`.
    pub(crate) fn resource_limit(self, what: &str, limit: LimitExceeded) -> FormatError {
        FormatError::write(self.0, DiagnosticKind::Limit, format!("{what} {limit}"))
            .exceeding(limit)
    }

    /// The source is internally inconsistent or incomplete, for example a
    /// model without positions for every atom.
    pub(crate) fn invalid_model(self, message: impl Into<String>) -> FormatError {
        FormatError::write(self.0, DiagnosticKind::Inconsistent, message)
    }

    /// The output sink failed.
    pub(crate) fn io(self, error: std::io::Error) -> FormatError {
        FormatError::sink(self.0, error)
    }
}

/// What a structure writer does with state its format cannot represent, such
/// as a periodic cell in a Molfile. Representable state is always written.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum UnrepresentablePolicy {
    /// Fail rather than write a record that silently lacks the state.
    #[default]
    Fail,
    /// Write the record without the state.
    Omit,
}

/// Fails for ensemble member or trajectory frame state that `format` cannot
/// carry, unless the caller chose to omit it. `member` is `None` for a format
/// that writes member state itself.
pub(crate) fn check_structure_state(
    format: Format,
    kind: StructureKind<'_>,
    member: Option<UnrepresentablePolicy>,
    frame: UnrepresentablePolicy,
) -> Result<(), FormatError> {
    let field = match kind {
        StructureKind::Member { energy, .. } if member == Some(UnrepresentablePolicy::Fail) => {
            // A lone record cannot carry a weight relative to other members.
            Some(if energy.is_some() {
                StructureField::Energy
            } else {
                StructureField::Weight
            })
        }
        StructureKind::Frame {
            time,
            step,
            velocities,
            forces,
            ..
        } if frame == UnrepresentablePolicy::Fail => {
            if time.is_some() {
                Some(StructureField::Time)
            } else if step.is_some() {
                Some(StructureField::Step)
            } else if velocities.is_some() {
                Some(StructureField::Velocities)
            } else if forces.is_some() {
                Some(StructureField::Forces)
            } else {
                None
            }
        }
        _ => None,
    };
    match field {
        None => Ok(()),
        Some(field) => Err(FormatError::write(
            format,
            DiagnosticKind::Unsupported,
            format!("{format} cannot represent the {field}; choose to omit it explicitly"),
        )
        .with_detail(FormatErrorDetail::UnrepresentableState(field))),
    }
}
