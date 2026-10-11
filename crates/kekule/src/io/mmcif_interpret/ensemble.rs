use std::collections::BTreeSet;
use std::sync::Arc;

use crate::structure::{Ensemble, EnsembleMember};

use super::super::{MmcifBlock, MmcifDocument};
use super::atom_site::coordinate_model_ids;
use super::types::{
    MmcifEnsembleInterpretError, MmcifEnsembleInterpretOptions, MmcifEnsembleInterpretation,
    MmcifInterpretOptions, MmcifInterpretation, MmcifInterpretationReport,
};
use super::PreparedBlock;

/// Interprets explicitly selected or all coordinate models as one
/// shared-topology non-temporal ensemble.
pub(crate) fn interpret_mmcif_ensemble(
    document: &MmcifDocument,
    options: MmcifEnsembleInterpretOptions,
) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
    interpret_mmcif_ensemble_block(atom_site_block(document)?, options)
}

fn atom_site_block(document: &MmcifDocument) -> Result<&MmcifBlock, MmcifEnsembleInterpretError> {
    let mut blocks = document
        .blocks()
        .iter()
        .filter(|block| block.has_category("_atom_site"));
    let block = blocks
        .next()
        .ok_or(MmcifEnsembleInterpretError::NoCoordinateModels)?;
    if blocks.next().is_some() {
        return Err(MmcifEnsembleInterpretError::MultipleAtomSiteBlocks);
    }
    Ok(block)
}

/// Interprets explicitly selected or all coordinate models in one block as
/// one shared-topology non-temporal ensemble.
pub(crate) fn interpret_mmcif_ensemble_block(
    block: &MmcifBlock,
    options: MmcifEnsembleInterpretOptions,
) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
    if options.model_ids.as_ref().is_some_and(Vec::is_empty) {
        return Err(MmcifEnsembleInterpretError::EmptyModelSelection);
    }
    if !block.has_category("_atom_site") {
        return Err(MmcifEnsembleInterpretError::NoCoordinateModels);
    }
    let available =
        coordinate_model_ids(block).map_err(|error| MmcifEnsembleInterpretError::Model {
            model_id: "<model inventory>".to_owned(),
            error,
        })?;
    if available.is_empty() {
        return Err(MmcifEnsembleInterpretError::NoCoordinateModels);
    }
    let selected = options.model_ids.unwrap_or_else(|| available.clone());
    let mut seen = BTreeSet::new();
    for model in &selected {
        if !seen.insert(model.clone()) {
            return Err(MmcifEnsembleInterpretError::DuplicateRequestedModel(
                model.clone(),
            ));
        }
        if !available.contains(model) {
            return Err(MmcifEnsembleInterpretError::UnknownRequestedModel(
                model.clone(),
            ));
        }
    }

    let mut prepared = PreparedBlock::new(
        block,
        options.strict_entity_metadata,
        &options.altloc_policy,
    )
    .map_err(|error| MmcifEnsembleInterpretError::Model {
        model_id: selected[0].clone(),
        error,
    })?;
    let mut selected = selected.into_iter();
    let first_id = selected.next().expect("validated nonempty model selection");
    let first = prepared.interpret_model(&first_id).map_err(|error| {
        MmcifEnsembleInterpretError::Model {
            model_id: first_id,
            error,
        }
    })?;
    let remaining = selected.map(|model_id| {
        prepared
            .interpret_model(&model_id)
            .map_err(|error| MmcifEnsembleInterpretError::Model { model_id, error })
    });
    assemble(first, remaining)
}

/// Explicit caller-selected realizations, never automatic altloc enumeration.
pub(crate) fn interpret_mmcif_conformations(
    document: &MmcifDocument,
    selections: &[MmcifInterpretOptions],
) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
    interpret_mmcif_conformations_block(atom_site_block(document)?, selections)
}

pub(crate) fn interpret_mmcif_conformations_block(
    block: &MmcifBlock,
    selections: &[MmcifInterpretOptions],
) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
    if selections.is_empty() {
        return Err(MmcifEnsembleInterpretError::EmptyConformationSelection);
    }
    let mut interpreted = selections.iter().enumerate().map(|(selection, options)| {
        super::interpret_mmcif_block(block, options.clone())
            .map_err(|error| MmcifEnsembleInterpretError::Conformation { selection, error })
    });
    let first = interpreted.next().expect("validated nonempty selections")?;
    assemble(first, interpreted)
}

fn assemble(
    first: MmcifInterpretation,
    remaining: impl Iterator<Item = Result<MmcifInterpretation, MmcifEnsembleInterpretError>>,
) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
    let shared_topology = first.model.shared_topology();
    let shared_atom_identity = provenance_identity(&first.report);
    let mut ensemble = Ensemble::new(Arc::clone(&shared_topology));
    let (first_model, first_report) = first.into_parts();
    ensemble
        .push(EnsembleMember::from_model(first_model))
        .map_err(|error| MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
    let mut reports = vec![first_report];
    for interpreted in remaining {
        let (model, report) = interpreted?.into_parts();
        let model_id = report
            .selected_model()
            .expect("interpreted model has a source ID")
            .to_owned();
        let atom_identity = provenance_identity(&report);
        if atom_identity != shared_atom_identity {
            let error = if atom_identity.sorted_atoms() != shared_atom_identity.sorted_atoms() {
                MmcifEnsembleInterpretError::InconsistentAtomSet { model_id }
            } else {
                MmcifEnsembleInterpretError::InconsistentDenseAtomOrder { model_id }
            };
            return Err(error);
        }
        if !shared_topology.same_layout(model.topology()) {
            return Err(MmcifEnsembleInterpretError::InconsistentTopology { model_id });
        }
        if shared_topology.atom_ids() != model.topology().atom_ids()
            || shared_topology.bond_ids() != model.topology().bond_ids()
        {
            return Err(MmcifEnsembleInterpretError::InconsistentDenseAtomOrder { model_id });
        }
        ensemble
            .push(EnsembleMember::from_model(model))
            .map_err(|error| MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
        reports.push(report);
    }
    Ok(MmcifEnsembleInterpretation { ensemble, reports })
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ProvenanceAtomIdentity {
    atom_name: String,
    component_id: String,
    asym_id: String,
    auth_asym_id: Option<String>,
    entity_id: Option<String>,
    label_sequence_id: Option<i32>,
    author_sequence_id: Option<String>,
    insertion_code: Option<String>,
    occurrence: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ProvenanceIdentity {
    atoms: Vec<ProvenanceAtomIdentity>,
}

impl ProvenanceIdentity {
    fn sorted_atoms(&self) -> Vec<ProvenanceAtomIdentity> {
        let mut atoms = self.atoms.clone();
        atoms.sort_unstable();
        atoms
    }
}

fn provenance_identity(report: &MmcifInterpretationReport) -> ProvenanceIdentity {
    ProvenanceIdentity {
        atoms: report
            .instances
            .iter()
            .flat_map(|instance| {
                instance.atoms.iter().map(|atom| ProvenanceAtomIdentity {
                    atom_name: atom.atom_name.clone(),
                    component_id: atom.component_id.clone(),
                    asym_id: atom.asym_id.clone(),
                    auth_asym_id: atom.auth_asym_id.clone(),
                    entity_id: atom.entity_id.clone(),
                    label_sequence_id: atom.label_sequence_id,
                    author_sequence_id: atom.author_sequence_id.clone(),
                    insertion_code: atom.insertion_code.clone(),
                    occurrence: atom.occurrence,
                })
            })
            .collect(),
    }
}
