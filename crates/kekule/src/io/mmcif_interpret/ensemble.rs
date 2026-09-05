use std::collections::BTreeSet;
use std::sync::Arc;

use crate::structure::{Ensemble, EnsembleMember, Positions};

use super::super::{MmcifBlock, MmcifDocument};
use super::atom_site::coordinate_model_ids;
use super::interpret_mmcif_block;
use super::types::{
    MmcifEnsembleInterpretError, MmcifEnsembleInterpretOptions, MmcifEnsembleInterpretation,
    MmcifInterpretOptions, MmcifInterpretationReport, MmcifModelSelection,
};

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

    let mut interpreted = Vec::with_capacity(selected.len());
    for model_id in &selected {
        let interpretation = interpret_mmcif_block(
            block,
            MmcifInterpretOptions {
                strict_entity_metadata: options.strict_entity_metadata,
                altloc_policy: options.altloc_policy.clone(),
                model_selection: MmcifModelSelection::Select(model_id.clone()),
            },
        )
        .map_err(|error| MmcifEnsembleInterpretError::Model {
            model_id: model_id.clone(),
            error,
        })?;
        interpreted.push(interpretation);
    }

    let first = interpreted
        .first()
        .ok_or(MmcifEnsembleInterpretError::EmptyModelSelection)?;
    let shared_topology = first.model.shared_topology();
    let shared_atom_identity = provenance_identity(&first.report);
    let mut ensemble = Ensemble::new(Arc::clone(&shared_topology));
    let mut reports = Vec::with_capacity(interpreted.len());
    for interpretation in interpreted {
        let (model, report) = interpretation.to_parts();
        let model_id = report.selected_model().unwrap_or("<unknown>").to_owned();
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
        let positions = Positions::new(model.positions().values())
            .map_err(MmcifEnsembleInterpretError::Position)?;
        let mut member = EnsembleMember::new(positions);
        member.set_cell(model.cell().copied());
        member
            .set_properties(model.properties().clone())
            .map_err(|error| MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
        ensemble
            .push(member)
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
