use std::collections::BTreeSet;
use std::sync::Arc;

use crate::structure::{Ensemble, EnsembleMember};

use super::super::{MmcifBlock, MmcifDocument};
use super::atom_site::coordinate_model_ids;
use super::types::{
    MmcifEnsembleInterpretError, MmcifEnsembleInterpretOptions, MmcifEnsembleInterpretation,
    MmcifInterpretationReport,
};
use super::PreparedBlock;

/// Interprets explicitly selected or all coordinate models as one
/// shared-topology non-temporal ensemble.
pub(crate) fn interpret_mmcif_ensemble(
    document: &MmcifDocument,
    options: MmcifEnsembleInterpretOptions,
) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
    let blocks = document
        .blocks()
        .iter()
        .filter(|block| block.loop_with_tag("_atom_site.type_symbol").is_some())
        .collect::<Vec<_>>();
    if blocks.is_empty() {
        return Err(MmcifEnsembleInterpretError::NoCoordinateModels);
    }
    if blocks.len() > 1 {
        return Err(MmcifEnsembleInterpretError::MultipleAtomSiteBlocks);
    }
    interpret_mmcif_ensemble_block(blocks[0], options)
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
    if block.loop_with_tag("_atom_site.type_symbol").is_none() {
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
    let shared_topology = first.model.shared_topology();
    let shared_atom_identity = provenance_identity(&first.report);
    let mut ensemble = Ensemble::new(Arc::clone(&shared_topology));
    let (first_model, first_report) = first.to_parts();
    ensemble
        .push(EnsembleMember::from_model(first_model))
        .map_err(|error| MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
    let mut reports = vec![first_report];
    for model_id in selected {
        let (model, report) = prepared
            .interpret_model(&model_id)
            .map_err(|error| MmcifEnsembleInterpretError::Model {
                model_id: model_id.clone(),
                error,
            })?
            .to_parts();
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
