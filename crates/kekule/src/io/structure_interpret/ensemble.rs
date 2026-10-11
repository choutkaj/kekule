use std::collections::BTreeSet;
use std::sync::Arc;

use crate::structure::{Ensemble, EnsembleMember};

use super::types::{
    StructureEnsembleInterpretError, StructureEnsembleInterpretOptions,
    StructureEnsembleInterpretation, StructureInterpretError, StructureInterpretOptions,
    StructureInterpretation, StructureInterpretationReport,
};
use super::{PreparedSource, StructureSource};

/// Interprets explicitly selected or all of the `available` coordinate models
/// of one source as one shared-topology non-temporal ensemble.
pub(super) fn interpret_ensemble(
    available: Vec<String>,
    options: StructureEnsembleInterpretOptions,
    source: impl FnOnce(
        &StructureEnsembleInterpretOptions,
    ) -> Result<StructureSource, StructureInterpretError>,
) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
    if options.model_ids.as_ref().is_some_and(Vec::is_empty) {
        return Err(StructureEnsembleInterpretError::EmptyModelSelection);
    }
    if available.is_empty() {
        return Err(StructureEnsembleInterpretError::NoCoordinateModels);
    }
    let selected = options
        .model_ids
        .clone()
        .unwrap_or_else(|| available.clone());
    let mut seen = BTreeSet::new();
    for model in &selected {
        if !seen.insert(model.clone()) {
            return Err(StructureEnsembleInterpretError::DuplicateRequestedModel(
                model.clone(),
            ));
        }
        if !available.contains(model) {
            return Err(StructureEnsembleInterpretError::UnknownRequestedModel(
                model.clone(),
            ));
        }
    }

    let mut prepared = source(&options)
        .and_then(|source| PreparedSource::new(source, &options.altloc_policy))
        .map_err(|error| StructureEnsembleInterpretError::Model {
            model_id: selected[0].clone(),
            error,
        })?;
    let mut selected = selected.into_iter();
    let first_id = selected.next().expect("validated nonempty model selection");
    let first = prepared.interpret_model(&first_id).map_err(|error| {
        StructureEnsembleInterpretError::Model {
            model_id: first_id,
            error,
        }
    })?;
    let remaining = selected.map(|model_id| {
        prepared
            .interpret_model(&model_id)
            .map_err(|error| StructureEnsembleInterpretError::Model { model_id, error })
    });
    assemble(first, remaining)
}

/// Explicit caller-selected alternate conformations, never automatic altloc
/// enumeration.
pub(super) fn interpret_conformations(
    selections: &[StructureInterpretOptions],
    interpret: impl Fn(
        &StructureInterpretOptions,
    ) -> Result<StructureInterpretation, StructureInterpretError>,
) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
    if selections.is_empty() {
        return Err(StructureEnsembleInterpretError::EmptyConformationSelection);
    }
    let mut interpreted = selections.iter().enumerate().map(|(selection, options)| {
        interpret(options)
            .map_err(|error| StructureEnsembleInterpretError::Conformation { selection, error })
    });
    let first = interpreted.next().expect("validated nonempty selections")?;
    assemble(first, interpreted)
}

fn assemble(
    first: StructureInterpretation,
    remaining: impl Iterator<Item = Result<StructureInterpretation, StructureEnsembleInterpretError>>,
) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
    let shared_topology = first.model.shared_topology();
    let shared_atom_identity = provenance_identity(&first.report);
    let mut ensemble = Ensemble::new(Arc::clone(&shared_topology));
    let (first_model, first_report) = first.into_parts();
    ensemble
        .push(EnsembleMember::from_model(first_model))
        .map_err(|error| StructureEnsembleInterpretError::Ensemble(Box::new(error)))?;
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
                StructureEnsembleInterpretError::InconsistentAtomSet { model_id }
            } else {
                StructureEnsembleInterpretError::InconsistentDenseAtomOrder { model_id }
            };
            return Err(error);
        }
        if !shared_topology.same_layout(model.topology()) {
            return Err(StructureEnsembleInterpretError::InconsistentTopology { model_id });
        }
        if shared_topology.atom_ids() != model.topology().atom_ids()
            || shared_topology.bond_ids() != model.topology().bond_ids()
        {
            return Err(StructureEnsembleInterpretError::InconsistentDenseAtomOrder { model_id });
        }
        ensemble
            .push(EnsembleMember::from_model(model))
            .map_err(|error| StructureEnsembleInterpretError::Ensemble(Box::new(error)))?;
        reports.push(report);
    }
    Ok(StructureEnsembleInterpretation { ensemble, reports })
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

fn provenance_identity(report: &StructureInterpretationReport) -> ProvenanceIdentity {
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
