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

/// What a structure writer left out of its output: content the caller chose
/// [`UnrepresentablePolicy::Omit`] for, or, for a writer without such a
/// choice, such as `kekule-traj`'s XYZ writer, everything its format cannot
/// carry.
///
/// A writer reports every omission it makes, once per kind, with the number
/// of entities it affects summed over everything one call writes. An empty
/// report means the output holds everything the source does that the format
/// can carry.
///
/// Collecting omissions into a report merges repeated kinds by summing their
/// counts and drops zero counts, as writers do.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct WriteReport {
    omissions: Vec<WriteOmission>,
}

impl WriteReport {
    /// The omissions, in the order the writer first made each kind.
    pub fn omissions(&self) -> &[WriteOmission] {
        &self.omissions
    }

    /// Whether the writer omitted nothing.
    pub fn is_empty(&self) -> bool {
        self.omissions.is_empty()
    }
}

impl FromIterator<WriteOmission> for WriteReport {
    fn from_iter<I: IntoIterator<Item = WriteOmission>>(omissions: I) -> Self {
        let mut collected = Omissions::default();
        for omission in omissions {
            collected.record(omission);
        }
        collected.into_report()
    }
}

/// One kind of content a writer omitted, with how much of it.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WriteOmission {
    /// Implicit-hydrogen counts of this many atoms that readers of the format
    /// do not recover: PDB and mmCIF readers infer counts from valence in
    /// place of fixed ones, and XYZ readers fix every count at zero.
    ImplicitHydrogens {
        /// The atoms whose counts were left out.
        atoms: usize,
    },
    /// Nonzero formal charges of this many atoms.
    FormalCharges {
        /// The atoms written without their charges.
        atoms: usize,
    },
    /// Isotopes the format cannot name, of this many atoms.
    Isotopes {
        /// The atoms written as their element alone.
        atoms: usize,
    },
    /// Radical electrons of this many atoms.
    Radicals {
        /// The atoms written without their radicals.
        atoms: usize,
    },
    /// Atom-map numbers of this many atoms.
    AtomMaps {
        /// The atoms written without their atom maps.
        atoms: usize,
    },
    /// Represented stereo of this many molecule instances.
    Stereo {
        /// The molecule instances written without their stereo.
        molecules: usize,
    },
    /// Named hierarchy atom sites of this many atoms: their atom names and
    /// their residues and chains.
    AtomSites {
        /// The atoms written without their sites.
        atoms: usize,
    },
    /// Bonds the format's connectivity records cannot state, left out.
    Bonds {
        /// The bonds left out.
        bonds: usize,
    },
    /// Residues whose bonds differ from the component dictionary that
    /// readers apply to them, so readers will see the dictionary's bonds.
    DictionaryBonds {
        /// The residues affected.
        residues: usize,
    },
    /// A periodic cell the format cannot carry, of this many written
    /// structures or blocks.
    Cell {
        /// The structures or blocks written without their cell.
        structures: usize,
    },
    /// Ensemble member or trajectory frame state the format cannot carry.
    StructureState {
        /// The state left out.
        field: StructureField,
        /// The members or frames whose state was left out.
        structures: usize,
    },
}

/// Collects a writer's omissions, merging repeated kinds by summing counts.
#[derive(Debug, Default)]
pub(crate) struct Omissions(Vec<WriteOmission>);

impl Omissions {
    /// Records `omission`, adding its count to an earlier one of its kind.
    /// Zero counts record nothing.
    pub(crate) fn record(&mut self, omission: WriteOmission) {
        use WriteOmission as O;
        let count = |omission: &WriteOmission| match *omission {
            O::ImplicitHydrogens { atoms }
            | O::FormalCharges { atoms }
            | O::Isotopes { atoms }
            | O::Radicals { atoms }
            | O::AtomMaps { atoms }
            | O::AtomSites { atoms } => atoms,
            O::Stereo { molecules } => molecules,
            O::Bonds { bonds } => bonds,
            O::DictionaryBonds { residues } => residues,
            O::Cell { structures } | O::StructureState { structures, .. } => structures,
        };
        let added = count(&omission);
        if added == 0 {
            return;
        }
        let same_kind = |existing: &WriteOmission| match (existing, &omission) {
            (O::StructureState { field: a, .. }, O::StructureState { field: b, .. }) => a == b,
            _ => std::mem::discriminant(existing) == std::mem::discriminant(&omission),
        };
        let Some(existing) = self.0.iter_mut().find(|existing| same_kind(existing)) else {
            self.0.push(omission);
            return;
        };
        match existing {
            O::ImplicitHydrogens { atoms }
            | O::FormalCharges { atoms }
            | O::Isotopes { atoms }
            | O::Radicals { atoms }
            | O::AtomMaps { atoms }
            | O::AtomSites { atoms } => *atoms += added,
            O::Stereo { molecules } => *molecules += added,
            O::Bonds { bonds } => *bonds += added,
            O::DictionaryBonds { residues } => *residues += added,
            O::Cell { structures } | O::StructureState { structures, .. } => *structures += added,
        }
    }

    /// Records the member or frame state `kind` has that the policies omit;
    /// call it after [`check_structure_state`] has passed.
    pub(crate) fn record_structure_state(
        &mut self,
        kind: StructureKind<'_>,
        member: Option<UnrepresentablePolicy>,
        frame: UnrepresentablePolicy,
    ) {
        let mut fields = Vec::new();
        match kind {
            StructureKind::Member { energy, .. } if member == Some(UnrepresentablePolicy::Omit) => {
                // A lone record cannot carry a weight relative to other
                // members, so its weight is always left out.
                fields.push(StructureField::Weight);
                if energy.is_some() {
                    fields.push(StructureField::Energy);
                }
            }
            StructureKind::Frame {
                time,
                step,
                velocities,
                forces,
                ..
            } if frame == UnrepresentablePolicy::Omit => {
                for (present, field) in [
                    (time.is_some(), StructureField::Time),
                    (step.is_some(), StructureField::Step),
                    (velocities.is_some(), StructureField::Velocities),
                    (forces.is_some(), StructureField::Forces),
                ] {
                    if present {
                        fields.push(field);
                    }
                }
            }
            _ => {}
        }
        for field in fields {
            self.record(WriteOmission::StructureState {
                field,
                structures: 1,
            });
        }
    }

    pub(crate) fn into_report(self) -> WriteReport {
        WriteReport { omissions: self.0 }
    }
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
