//! Public writer API and orchestration of planning, preparation, and emission.
use crate::core::BondOrder;
use crate::io::mmcif_interpret::{
    MmcifEnsembleInterpretation, MmcifEntityKind, MmcifInterpretationReport,
};
use crate::structure::{Ensemble, Model};
use crate::topology::{InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeInstanceId};
use std::collections::BTreeMap;
use std::fmt;
use std::io::Write;

mod emit;
mod entities;
mod prepare;
use emit::{render_model_to, write_atom_rows, write_block_end, write_block_start};
use entities::{
    ensemble_entity_plan, entity_plan_from_report, generic_entity_plan,
    normalize_entity_classifications,
};
use prepare::{prepare_model, validate_ensemble_member};

const MAX_COORDINATE_PRECISION: usize = 15;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifWriteOptions {
    pub block_name: String,
    pub coordinate_precision: usize,
}

impl Default for MmcifWriteOptions {
    fn default() -> Self {
        Self {
            block_name: "model".to_owned(),
            coordinate_precision: 3,
        }
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifWriteError {
    InvalidBlockName(String),
    CoordinatePrecisionTooLarge(usize),
    InvalidModel(String),
    InvalidHierarchy {
        message: String,
    },
    MissingEntityClassification(MoleculeInstanceId),
    DuplicateEntityClassification(MoleculeInstanceId),
    ConflictingEntityClassifications {
        molecule: MoleculeInstanceId,
        classifications: Vec<MmcifEntityKind>,
    },
    UnsupportedEntityClassification {
        molecule: MoleculeInstanceId,
        classification: String,
    },
    UnresolvedCanonicalEntityClassification {
        molecule: MoleculeInstanceId,
        classification: MoleculeClass,
    },
    ConflictingAsymEntityIds {
        asym_id: String,
        entity_ids: Vec<String>,
    },
    ConflictingAsymEntityClassifications {
        asym_id: String,
        classifications: Vec<MmcifEntityKind>,
    },
    ConflictingSourceEntityClassifications {
        entity_id: String,
        classifications: Vec<MmcifEntityKind>,
    },
    UnknownClassifiedMolecule(MoleculeInstanceId),
    DuplicateAsymId(String),
    MissingAtomSite(InstanceAtomId),
    DuplicateAtomSite(InstanceAtomId),
    InconsistentAtomSite {
        atom: InstanceAtomId,
        field: &'static str,
    },
    InvalidGroupPdb {
        atom: InstanceAtomId,
        value: String,
    },
    DuplicateAtomIdentity(InstanceAtomId),
    MissingAtomProvenance(InstanceAtomId),
    DuplicateAtomProvenance(InstanceAtomId),
    UnknownAtomProvenance(InstanceAtomId),
    UnsupportedAtomField {
        atom: InstanceAtomId,
        field: &'static str,
    },
    FormalChargeOutOfRange {
        atom: InstanceAtomId,
        charge: i8,
    },
    UnsupportedStereo(MoleculeInstanceId),
    UnsupportedBondOrder {
        bond: InstanceBondId,
        order: BondOrder,
    },
    AmbiguousConnectionSelector(InstanceAtomId),
    UnsupportedTextValue {
        field: &'static str,
    },
    EmptyEnsemble,
    ReportCountMismatch {
        expected: usize,
        actual: usize,
    },
    ClassificationCountMismatch {
        expected: usize,
        actual: usize,
    },
    IncompatibleEnsembleMember {
        member: usize,
        field: &'static str,
    },
    Io {
        kind: std::io::ErrorKind,
        message: String,
    },
}

impl fmt::Display for MmcifWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBlockName(name) => {
                write!(f, "invalid mmCIF data block name `{name}`")
            }
            Self::CoordinatePrecisionTooLarge(precision) => write!(
                f,
                "mmCIF coordinate precision {precision} exceeds the supported maximum of {MAX_COORDINATE_PRECISION}"
            ),
            Self::InvalidModel(message) => write!(f, "invalid molecular model: {message}"),
            Self::InvalidHierarchy { message } => {
                write!(f, "invalid topology hierarchy: {message}")
            }
            Self::MissingEntityClassification(molecule) => write!(
                f,
                "{molecule} has no explicit mmCIF entity classification"
            ),
            Self::DuplicateEntityClassification(molecule) => write!(
                f,
                "mmCIF entity semantics classify {molecule} more than once"
            ),
            Self::ConflictingEntityClassifications {
                molecule,
                classifications,
            } => write!(
                f,
                "{molecule} has conflicting mmCIF entity classifications {classifications:?}"
            ),
            Self::UnsupportedEntityClassification {
                molecule,
                classification,
            } => write!(
                f,
                "{molecule} has unsupported mmCIF entity classification `{classification}`"
            ),
            Self::UnresolvedCanonicalEntityClassification {
                molecule,
                classification,
            } => write!(
                f,
                "cannot derive an unambiguous mmCIF entity kind for {molecule} from canonical class {classification:?}"
            ),
            Self::ConflictingAsymEntityIds {
                asym_id,
                entity_ids,
            } => write!(
                f,
                "mmCIF structural instance `{asym_id}` has conflicting source entity IDs {entity_ids:?}"
            ),
            Self::ConflictingAsymEntityClassifications {
                asym_id,
                classifications,
            } => write!(
                f,
                "mmCIF structural instance `{asym_id}` has conflicting entity classifications {classifications:?}"
            ),
            Self::ConflictingSourceEntityClassifications {
                entity_id,
                classifications,
            } => write!(
                f,
                "source mmCIF entity `{entity_id}` has conflicting classifications {classifications:?}"
            ),
            Self::UnknownClassifiedMolecule(molecule) => write!(
                f,
                "mmCIF entity semantics reference unknown {molecule}"
            ),
            Self::DuplicateAsymId(id) => {
                write!(f, "duplicate mmCIF structural-instance ID `{id}`")
            }
            Self::MissingAtomSite(atom) => write!(f, "{atom} has no biomolecular atom site"),
            Self::DuplicateAtomSite(atom) => {
                write!(f, "{atom} appears in more than one biomolecular atom site")
            }
            Self::InconsistentAtomSite { atom, field } => {
                write!(f, "{atom} has inconsistent atom-site {field}")
            }
            Self::InvalidGroupPdb { atom, value } => write!(
                f,
                "{atom} has unsupported _atom_site.group_PDB value `{value}`"
            ),
            Self::DuplicateAtomIdentity(atom) => write!(
                f,
                "{atom} duplicates an mmCIF atom identity within one residue"
            ),
            Self::MissingAtomProvenance(atom) => {
                write!(f, "{atom} has no atom-level mmCIF source provenance")
            }
            Self::DuplicateAtomProvenance(atom) => {
                write!(f, "{atom} has duplicate atom-level mmCIF source provenance")
            }
            Self::UnknownAtomProvenance(atom) => {
                write!(f, "mmCIF source provenance references unknown {atom}")
            }
            Self::UnsupportedAtomField { atom, field } => {
                write!(f, "{atom} has unsupported atom field `{field}`")
            }
            Self::FormalChargeOutOfRange { atom, charge } => write!(
                f,
                "{atom} formal charge {charge} is outside the PDBx/mmCIF range -8..=8"
            ),
            Self::UnsupportedStereo(molecule) => write!(
                f,
                "{molecule} contains stereochemistry not represented by the foundational mmCIF writer"
            ),
            Self::UnsupportedBondOrder { bond, order } => {
                write!(f, "{bond} has unsupported mmCIF bond order {order:?}")
            }
            Self::AmbiguousConnectionSelector(atom) => write!(
                f,
                "{atom} cannot be selected unambiguously by an mmCIF struct_conn partner"
            ),
            Self::UnsupportedTextValue { field } => write!(
                f,
                "{field} contains a text value that cannot be emitted as a single mmCIF token"
            ),
            Self::EmptyEnsemble => f.write_str("cannot write an empty ensemble as mmCIF"),
            Self::ReportCountMismatch { expected, actual } => write!(
                f,
                "mmCIF writing requires {expected} interpretation reports, but received {actual}"
            ),
            Self::ClassificationCountMismatch { expected, actual } => write!(
                f,
                "mmCIF model collection writing requires {expected} classification sets, but received {actual}"
            ),
            Self::IncompatibleEnsembleMember { member, field } => write!(
                f,
                "ensemble member {member} has incompatible mmCIF {field}"
            ),
            Self::Io { message, .. } => write!(f, "mmCIF output failed: {message}"),
        }
    }
}

impl std::error::Error for MmcifWriteError {}

/// Expert mmCIF entity-kind overrides for canonical molecule instances.
///
/// Ordinary writing derives kinds from canonical topology classification.
/// Entries in this collection take precedence where exact source distinctions
/// or an expert format-specific choice is required.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MmcifEntityClassifications {
    kinds: BTreeMap<MoleculeInstanceId, MmcifEntityKind>,
}

impl MmcifEntityClassifications {
    pub const fn new() -> Self {
        Self {
            kinds: BTreeMap::new(),
        }
    }

    /// Assigns one explicit mmCIF entity kind to a molecule instance.
    pub fn insert(
        &mut self,
        molecule: MoleculeInstanceId,
        kind: MmcifEntityKind,
    ) -> Result<(), MmcifWriteError> {
        if self.kinds.contains_key(&molecule) {
            return Err(MmcifWriteError::DuplicateEntityClassification(molecule));
        }
        self.kinds.insert(molecule, kind);
        Ok(())
    }

    /// Returns the assigned kind for one molecule instance.
    pub fn get(&self, molecule: MoleculeInstanceId) -> Option<&MmcifEntityKind> {
        self.kinds.get(&molecule)
    }

    /// Iterates over explicitly classified molecule instances in ID order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &MmcifEntityKind)> {
        self.kinds.iter().map(|(&molecule, kind)| (molecule, kind))
    }

    /// Returns the number of explicitly classified molecule instances.
    pub fn len(&self) -> usize {
        self.kinds.len()
    }

    /// Returns whether no molecule instance has been classified.
    pub fn is_empty(&self) -> bool {
        self.kinds.is_empty()
    }
}

pub fn write_mmcif_model(
    model: &Model,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_model_to(&mut output, model, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_model_to(
    writer: &mut impl Write,
    model: &Model,
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    validate_options(&options)?;
    let view = model.view();
    let classifications = normalize_entity_classifications(view, std::iter::empty())?;
    let plan = generic_entity_plan(view, &classifications)?;
    let prepared = prepare_model(view, plan)?;
    render_model_to(writer, &prepared, &options)
}

pub fn write_mmcif_model_with_classifications(
    model: &Model,
    classifications: &MmcifEntityClassifications,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_model_with_classifications_to(&mut output, model, classifications, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_model_with_classifications_to(
    writer: &mut impl Write,
    model: &Model,
    classifications: &MmcifEntityClassifications,
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    validate_options(&options)?;
    let view = model.view();
    let classifications = normalize_entity_classifications(
        view,
        classifications
            .iter()
            .map(|(molecule, kind)| (molecule, vec![kind.clone()])),
    )?;
    let plan = generic_entity_plan(view, &classifications)?;
    let prepared = prepare_model(view, plan)?;
    render_model_to(writer, &prepared, &options)
}

pub fn write_mmcif_model_with_report(
    model: &Model,
    report: &MmcifInterpretationReport,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_model_with_report_to(&mut output, model, report, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_model_with_report_to(
    writer: &mut impl Write,
    model: &Model,
    report: &MmcifInterpretationReport,
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    validate_options(&options)?;
    let view = model.view();
    let plan = entity_plan_from_report(view, report)?;
    let prepared = prepare_model(view, plan)?;
    render_model_to(writer, &prepared, &options)
}

pub fn write_mmcif_models(
    models: &[Model],
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_models_to(&mut output, models, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_models_with_classifications(
    models: &[Model],
    classifications: &[MmcifEntityClassifications],
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_models_with_classifications_to(&mut output, models, classifications, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_models_with_classifications_to(
    writer: &mut impl Write,
    models: &[Model],
    classifications: &[MmcifEntityClassifications],
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    if classifications.len() != models.len() {
        return Err(MmcifWriteError::ClassificationCountMismatch {
            expected: models.len(),
            actual: classifications.len(),
        });
    }
    write_independent_models_to(writer, models, Some(classifications), None, options)
}

pub fn write_mmcif_models_with_reports(
    models: &[Model],
    reports: &[MmcifInterpretationReport],
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_models_with_reports_to(&mut output, models, reports, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_models_with_reports_to(
    writer: &mut impl Write,
    models: &[Model],
    reports: &[MmcifInterpretationReport],
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    if reports.len() != models.len() {
        return Err(MmcifWriteError::ReportCountMismatch {
            expected: models.len(),
            actual: reports.len(),
        });
    }
    write_independent_models_to(writer, models, None, Some(reports), options)
}

pub fn write_mmcif_ensemble(
    ensemble: &Ensemble,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_ensemble_to(&mut output, ensemble, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_ensemble_with_classifications(
    ensemble: &Ensemble,
    classifications: &MmcifEntityClassifications,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_ensemble_with_classifications_to(&mut output, ensemble, classifications, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_ensemble_with_classifications_to(
    writer: &mut impl Write,
    ensemble: &Ensemble,
    classifications: &MmcifEntityClassifications,
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    write_ensemble_views_to(writer, ensemble, Some(classifications), None, options)
}

pub fn write_mmcif_ensemble_with_reports(
    ensemble: &Ensemble,
    reports: &[MmcifInterpretationReport],
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_ensemble_with_reports_to(&mut output, ensemble, reports, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_ensemble_with_reports_to(
    writer: &mut impl Write,
    ensemble: &Ensemble,
    reports: &[MmcifInterpretationReport],
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    if reports.len() != ensemble.len() {
        return Err(MmcifWriteError::ReportCountMismatch {
            expected: ensemble.len(),
            actual: reports.len(),
        });
    }
    write_ensemble_views_to(writer, ensemble, None, Some(reports), options)
}

pub fn write_mmcif_ensemble_interpretation(
    interpretation: &MmcifEnsembleInterpretation,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_ensemble_interpretation_to(&mut output, interpretation, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

pub fn write_mmcif_ensemble_interpretation_to(
    writer: &mut impl Write,
    interpretation: &MmcifEnsembleInterpretation,
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    write_mmcif_ensemble_with_reports_to(
        writer,
        interpretation.ensemble(),
        interpretation.reports(),
        options,
    )
}

pub fn write_mmcif_models_to(
    writer: &mut impl Write,
    models: &[Model],
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    write_independent_models_to(writer, models, None, None, options)
}

pub fn write_mmcif_ensemble_to(
    writer: &mut impl Write,
    ensemble: &Ensemble,
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    write_ensemble_views_to(writer, ensemble, None, None, options)
}

fn write_independent_models_to(
    writer: &mut impl Write,
    models: &[Model],
    classifications: Option<&[MmcifEntityClassifications]>,
    reports: Option<&[MmcifInterpretationReport]>,
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    validate_options(&options)?;
    for (index, model) in models.iter().enumerate() {
        let block_options = MmcifWriteOptions {
            block_name: format!("{}_{}", options.block_name, index + 1),
            coordinate_precision: options.coordinate_precision,
        };
        let view = model.view();
        let plan = if let Some(reports) = reports {
            entity_plan_from_report(view, &reports[index])?
        } else {
            let normalized = normalize_entity_classifications(
                view,
                classifications
                    .map(|all| all[index].iter().map(|(id, kind)| (id, vec![kind.clone()])))
                    .into_iter()
                    .flatten(),
            )?;
            generic_entity_plan(view, &normalized)?
        };
        let prepared = prepare_model(view, plan)?;
        render_model_to(writer, &prepared, &block_options)?;
    }
    Ok(())
}

fn write_ensemble_views_to(
    writer: &mut impl Write,
    ensemble: &Ensemble,
    classifications: Option<&MmcifEntityClassifications>,
    reports: Option<&[MmcifInterpretationReport]>,
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    validate_options(&options)?;
    let mut members = ensemble.members().enumerate();
    let (_, first_member) = members.next().ok_or(MmcifWriteError::EmptyEnsemble)?;
    let first_view = first_member.as_model();
    let first_plan = ensemble_entity_plan(first_view, classifications, reports.map(|all| &all[0]))?;
    let first = prepare_model(first_view, first_plan)?;

    write_block_start(writer, &first, &options)?;
    let mut atom_serial = 1u64;
    write_atom_rows(writer, &first.atoms, 1, &mut atom_serial, &options)?;

    for (index, member) in members {
        let view = member.as_model();
        let plan = ensemble_entity_plan(view, classifications, reports.map(|all| &all[index]))?;
        let candidate = prepare_model(view, plan)?;
        validate_ensemble_member(&first, &candidate, index + 1)?;
        write_atom_rows(
            writer,
            &candidate.atoms,
            index + 1,
            &mut atom_serial,
            &options,
        )?;
    }
    write_block_end(writer, &first)
}

fn validate_options(options: &MmcifWriteOptions) -> Result<(), MmcifWriteError> {
    if options.block_name.is_empty()
        || !options
            .block_name
            .chars()
            .all(|character| character.is_ascii_alphanumeric() || "_-.".contains(character))
    {
        return Err(MmcifWriteError::InvalidBlockName(
            options.block_name.clone(),
        ));
    }
    if options.coordinate_precision > MAX_COORDINATE_PRECISION {
        return Err(MmcifWriteError::CoordinatePrecisionTooLarge(
            options.coordinate_precision,
        ));
    }
    Ok(())
}

fn one_based_serial(raw: u32) -> u64 {
    u64::from(raw) + 1
}
