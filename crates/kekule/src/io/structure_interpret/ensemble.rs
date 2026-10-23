//! Assembles interpreted coordinate models or explicit conformation
//! selections into one shared-topology ensemble.
use std::collections::BTreeSet;
use std::sync::Arc;

use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::structure::{Ensemble, EnsembleError, EnsembleSample, EnsembleState};

use super::types::{
    ensemble_error, StructureEnsembleInterpretOptions, StructureEnsembleInterpretation,
    StructureInterpretOptions, StructureInterpretation, StructureInterpretationReport,
};
use super::{PreparedSource, StructureSource};

/// Interprets explicitly selected or all of the `available` coordinate models
/// of one source as one shared-topology non-temporal ensemble.
pub(super) fn interpret_ensemble(
    available: Vec<String>,
    options: StructureEnsembleInterpretOptions,
    source: impl FnOnce(&StructureEnsembleInterpretOptions) -> Result<StructureSource, FormatError>,
) -> Result<StructureEnsembleInterpretation, FormatError> {
    if available.is_empty() {
        return Err(no_coordinate_models());
    }
    let selected = options
        .model_ids()
        .map(<[String]>::to_vec)
        .unwrap_or_else(|| available.clone());
    let mut seen = BTreeSet::new();
    for model in &selected {
        if !seen.insert(model.clone()) {
            return Err(ensemble_error(
                DiagnosticKind::Inconsistent,
                format!("coordinate model `{model}` was requested more than once"),
                FormatErrorDetail::DuplicateRequestedModel {
                    model_id: model.clone(),
                },
            ));
        }
        if !available.contains(model) {
            return Err(ensemble_error(
                DiagnosticKind::Inconsistent,
                format!("coordinate model `{model}` is not present"),
                FormatErrorDetail::UnknownRequestedModel {
                    model_id: model.clone(),
                },
            ));
        }
    }

    let mut prepared = source(&options)
        .and_then(|source| {
            PreparedSource::new(
                source,
                options.altloc_policy(),
                options.cancellation().cloned(),
            )
        })
        .map_err(|error| in_model(error, selected[0].clone()))?;
    let mut selected = selected.into_iter();
    let first_id = selected.next().expect("validated nonempty model selection");
    let first = prepared
        .interpret_model(&first_id)
        .map_err(|error| in_model(error, first_id))?;
    let remaining = selected.map(|model_id| {
        prepared
            .interpret_model(&model_id)
            .map_err(|error| in_model(error, model_id))
    });
    assemble(first, remaining)
}

/// Explicit caller-selected alternate conformations, never automatic altloc
/// enumeration.
pub(super) fn interpret_conformations(
    selections: &[StructureInterpretOptions],
    interpret: impl Fn(&StructureInterpretOptions) -> Result<StructureInterpretation, FormatError>,
) -> Result<StructureEnsembleInterpretation, FormatError> {
    if selections.is_empty() {
        return Err(ensemble_error(
            DiagnosticKind::Inconsistent,
            "explicit conformation selection is empty",
            FormatErrorDetail::EmptyConformationSelection,
        ));
    }
    let mut interpreted = selections.iter().enumerate().map(|(selection, options)| {
        interpret(options)
            .map_err(|error| error.with_detail(FormatErrorDetail::Conformation { selection }))
    });
    let first = interpreted.next().expect("validated nonempty selections")?;
    assemble(first, interpreted)
}

/// An explicit coordinate-model selection is empty.
pub(super) fn empty_model_selection() -> FormatError {
    ensemble_error(
        DiagnosticKind::Inconsistent,
        "explicit coordinate model selection is empty",
        FormatErrorDetail::EmptyModelSelection,
    )
}

/// The source contains no coordinate models.
pub(super) fn no_coordinate_models() -> FormatError {
    ensemble_error(
        DiagnosticKind::Inconsistent,
        "input contains no coordinate models",
        FormatErrorDetail::NoCoordinateModels,
    )
}

/// Interpreting coordinate model `model_id` failed with `error`.
fn in_model(error: FormatError, model_id: String) -> FormatError {
    error.with_detail(FormatErrorDetail::Model { model_id })
}

fn inconsistent_atom_set(model_id: String) -> FormatError {
    ensemble_error(
        DiagnosticKind::Inconsistent,
        format!("coordinate model `{model_id}` has an inconsistent atom identity set"),
        FormatErrorDetail::InconsistentAtomSet { model_id },
    )
}

fn inconsistent_dense_atom_order(model_id: String) -> FormatError {
    ensemble_error(
        DiagnosticKind::Inconsistent,
        format!("coordinate model `{model_id}` has an inconsistent dense atom identity order"),
        FormatErrorDetail::InconsistentDenseAtomOrder { model_id },
    )
}

fn assembly_error(error: EnsembleError) -> FormatError {
    FormatError::interpret(
        Format::Mmcif,
        DiagnosticKind::Chemistry,
        "cannot assemble ensemble",
    )
    .with_source(error)
}

fn assemble(
    first: StructureInterpretation,
    remaining: impl Iterator<Item = Result<StructureInterpretation, FormatError>>,
) -> Result<StructureEnsembleInterpretation, FormatError> {
    let shared_topology = first.model.shared_topology();
    let shared_atom_identity = provenance_identity(&first.report);
    let mut ensemble = Ensemble::new(Arc::clone(&shared_topology));
    let (first_model, first_report) = first.into_parts();
    // Coordinate models are draws of one distribution: one state.
    ensemble
        .push_state(
            EnsembleSample::from_model(first_model),
            EnsembleState::new(),
        )
        .map_err(assembly_error)?;
    let mut reports = vec![first_report];
    for interpreted in remaining {
        let (model, report) = interpreted?.into_parts();
        let model_id = report
            .selected_model()
            .expect("interpreted model has a source ID")
            .to_owned();
        let atom_identity = provenance_identity(&report);
        if atom_identity != shared_atom_identity {
            return Err(
                if atom_identity.sorted_atoms() != shared_atom_identity.sorted_atoms() {
                    inconsistent_atom_set(model_id)
                } else {
                    inconsistent_dense_atom_order(model_id)
                },
            );
        }
        if !shared_topology.same_layout(model.topology()) {
            return Err(ensemble_error(
                DiagnosticKind::Inconsistent,
                format!(
                    "coordinate model `{model_id}` has inconsistent molecule partition, chemistry, connectivity, or hierarchy"
                ),
                FormatErrorDetail::InconsistentTopology { model_id },
            ));
        }
        if shared_topology.atom_ids() != model.topology().atom_ids()
            || shared_topology.bond_ids() != model.topology().bond_ids()
        {
            return Err(inconsistent_dense_atom_order(model_id));
        }
        ensemble
            .push_sample(EnsembleSample::from_model(model), 0)
            .map_err(assembly_error)?;
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
