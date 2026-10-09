//! Public writer API and orchestration of planning, preparation, and emission.
use crate::core::BondOrder;
use crate::io::mmcif_interpret::{
    MmcifEnsembleInterpretation, MmcifEntityKind, MmcifInterpretation, MmcifInterpretationReport,
};
use crate::structure::{AsModelView, Ensemble, Model, ModelView};
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

/// One mmCIF data block to write: a model or an ensemble, plus the entity
/// semantics used for its `_entity` and asymmetry assignments.
///
/// Convert from `&Model`, `&Ensemble`, `&MmcifInterpretation`, or
/// `&MmcifEnsembleInterpretation`; interpretations keep their source reports.
/// [`Self::model`] accepts any borrowed model view, such as an ensemble member.
#[derive(Debug, Clone, Copy)]
pub struct MmcifBlockSource<'a> {
    content: BlockContent<'a>,
    entities: EntitySemantics<'a>,
}

#[derive(Debug, Clone, Copy)]
enum BlockContent<'a> {
    Model(ModelView<'a>),
    Ensemble(&'a Ensemble),
}

#[derive(Debug, Clone, Copy)]
enum EntitySemantics<'a> {
    Canonical,
    Classifications(&'a MmcifEntityClassifications),
    Reports(&'a [MmcifInterpretationReport]),
}

impl<'a> MmcifBlockSource<'a> {
    /// One model block with entity kinds derived from canonical topology
    /// classification.
    pub fn model(model: &'a (impl AsModelView + ?Sized)) -> Self {
        Self {
            content: BlockContent::Model(model.as_model_view()),
            entities: EntitySemantics::Canonical,
        }
    }

    /// One multi-model block with one model per ensemble member.
    pub fn ensemble(ensemble: &'a Ensemble) -> Self {
        Self {
            content: BlockContent::Ensemble(ensemble),
            entities: EntitySemantics::Canonical,
        }
    }

    /// Overrides derived entity kinds for the listed molecule instances;
    /// omitted instances keep canonical classification.
    ///
    /// One mmCIF entity is assigned to each populated hierarchy chain (and one
    /// to each hierarchy-free instance), so instances touched by the same chain
    /// must have the same classification.
    #[must_use]
    pub fn with_classifications(mut self, classifications: &'a MmcifEntityClassifications) -> Self {
        self.entities = EntitySemantics::Classifications(classifications);
        self
    }

    /// Preserves interpreted mmCIF entity and asymmetry semantics: one report
    /// for a model block, one per member for an ensemble block.
    ///
    /// Atom-level provenance keeps one source entity and structural asymmetry
    /// consistent even when it spans several molecule instances; conflicting
    /// source identity is rejected. An auth-only source is normalized by
    /// copying author identifiers into the required label fields.
    #[must_use]
    pub fn with_reports(mut self, reports: &'a [MmcifInterpretationReport]) -> Self {
        self.entities = EntitySemantics::Reports(reports);
        self
    }
}

impl<'a> From<&'a Model> for MmcifBlockSource<'a> {
    fn from(model: &'a Model) -> Self {
        Self::model(model)
    }
}

impl<'a> From<&'a Ensemble> for MmcifBlockSource<'a> {
    fn from(ensemble: &'a Ensemble) -> Self {
        Self::ensemble(ensemble)
    }
}

impl<'a> From<&'a MmcifInterpretation> for MmcifBlockSource<'a> {
    fn from(interpretation: &'a MmcifInterpretation) -> Self {
        Self::model(interpretation.model())
            .with_reports(std::slice::from_ref(interpretation.report()))
    }
}

impl<'a> From<&'a MmcifEnsembleInterpretation> for MmcifBlockSource<'a> {
    fn from(interpretation: &'a MmcifEnsembleInterpretation) -> Self {
        Self::ensemble(interpretation.ensemble()).with_reports(interpretation.reports())
    }
}

/// Writes data blocks to a string. See [`write_mmcif_to`].
pub fn write_mmcif<'a, B: Into<MmcifBlockSource<'a>>>(
    blocks: impl IntoIterator<Item = B>,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = Vec::new();
    write_mmcif_to(&mut output, blocks, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

/// Writes one deterministic data block per source in input order.
///
/// A single block is named [`MmcifWriteOptions::block_name`]; several blocks
/// are named `{block_name}_1`, `{block_name}_2`, and so on. Rows follow
/// dense atom order.
pub fn write_mmcif_to<'a, B: Into<MmcifBlockSource<'a>>>(
    writer: &mut impl Write,
    blocks: impl IntoIterator<Item = B>,
    options: MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    validate_options(&options)?;
    let blocks = blocks.into_iter().map(Into::into).collect::<Vec<_>>();
    let numbered = blocks.len() > 1;
    for (index, block) in blocks.into_iter().enumerate() {
        let block_options = MmcifWriteOptions {
            block_name: if numbered {
                format!("{}_{}", options.block_name, index + 1)
            } else {
                options.block_name.clone()
            },
            coordinate_precision: options.coordinate_precision,
        };
        match block.content {
            BlockContent::Model(view) => {
                write_model_block_to(writer, view, block.entities, &block_options)?
            }
            BlockContent::Ensemble(ensemble) => {
                write_ensemble_block_to(writer, ensemble, block.entities, &block_options)?
            }
        }
    }
    Ok(())
}

fn write_model_block_to(
    writer: &mut impl Write,
    view: ModelView<'_>,
    entities: EntitySemantics<'_>,
    options: &MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    let plan = match entities {
        EntitySemantics::Reports(reports) => {
            let [report] = reports else {
                return Err(MmcifWriteError::ReportCountMismatch {
                    expected: 1,
                    actual: reports.len(),
                });
            };
            entity_plan_from_report(view, report)?
        }
        EntitySemantics::Classifications(classifications) => {
            let normalized = normalize_entity_classifications(
                view,
                classifications
                    .iter()
                    .map(|(molecule, kind)| (molecule, vec![kind.clone()])),
            )?;
            generic_entity_plan(view, &normalized)?
        }
        EntitySemantics::Canonical => {
            let normalized = normalize_entity_classifications(view, std::iter::empty())?;
            generic_entity_plan(view, &normalized)?
        }
    };
    let prepared = prepare_model(view, plan)?;
    render_model_to(writer, &prepared, options)
}

fn write_ensemble_block_to(
    writer: &mut impl Write,
    ensemble: &Ensemble,
    entities: EntitySemantics<'_>,
    options: &MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    let (classifications, reports) = match entities {
        EntitySemantics::Canonical => (None, None),
        EntitySemantics::Classifications(classifications) => (Some(classifications), None),
        EntitySemantics::Reports(reports) => {
            if reports.len() != ensemble.len() {
                return Err(MmcifWriteError::ReportCountMismatch {
                    expected: ensemble.len(),
                    actual: reports.len(),
                });
            }
            (None, Some(reports))
        }
    };
    let mut members = ensemble.iter().enumerate();
    let (_, first_member) = members.next().ok_or(MmcifWriteError::EmptyEnsemble)?;
    let first_view = first_member.as_model_view();
    let first_plan = ensemble_entity_plan(first_view, classifications, reports.map(|all| &all[0]))?;
    let first = prepare_model(first_view, first_plan)?;

    write_block_start(writer, &first, options)?;
    let mut atom_serial = 1u64;
    write_atom_rows(writer, &first.atoms, 1, &mut atom_serial, options)?;

    for (index, member) in members {
        let view = member.as_model_view();
        let plan = ensemble_entity_plan(view, classifications, reports.map(|all| &all[index]))?;
        let candidate = prepare_model(view, plan)?;
        validate_ensemble_member(&first, &candidate, index + 1)?;
        write_atom_rows(
            writer,
            &candidate.atoms,
            index + 1,
            &mut atom_serial,
            options,
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
