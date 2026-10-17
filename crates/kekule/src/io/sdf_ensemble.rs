//! Explicit assembly of interpreted SDF records into one ensemble.
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::structure::{Ensemble, EnsembleError, EnsembleMember};
use crate::units::{Quantity, KILOJOULE_PER_MOLE};

use super::sdf_document::SdfRecordInterpretation;
use super::sdf_write::{SDF_ENERGY_FIELD, SDF_WEIGHT_FIELD};

fn ensemble_error(error: EnsembleError) -> FormatError {
    FormatError::interpret(
        Format::Sdf,
        DiagnosticKind::Inconsistent,
        "SDF records form no ensemble",
    )
    .with_source(error)
}

/// Assembles interpreted SDF records, in order, into one ensemble on the
/// first record's topology. Every record must have the same topology layout
/// ([`Topology::same_layout`](crate::topology::Topology::same_layout)). Each
/// member takes its relative weight from the record's [`SDF_WEIGHT_FIELD`]
/// (1.0 when absent) and its energy from [`SDF_ENERGY_FIELD`] in kJ/mol (none
/// when absent), the fields the SDF writer emits for ensemble members.
pub(crate) fn sdf_ensemble<'a>(
    records: impl IntoIterator<Item = &'a SdfRecordInterpretation>,
) -> Result<Ensemble, FormatError> {
    let records = records.into_iter().collect::<Vec<_>>();
    let topology = records
        .first()
        .ok_or_else(|| ensemble_error(EnsembleError::EmptySource))?
        .model()
        .shared_topology();
    let mut members = Vec::with_capacity(records.len());
    for (index, record) in records.iter().enumerate() {
        if !record.model().topology().same_layout(&topology) {
            return Err(ensemble_error(EnsembleError::TopologyMismatch).in_record(index + 1));
        }
        let invalid = |kind, name: &str, message: String| {
            FormatError::interpret(Format::Sdf, kind, format!("data field `{name}` {message}"))
                .in_record(index + 1)
                .with_detail(FormatErrorDetail::DataField {
                    name: name.to_owned(),
                })
        };
        let number = |name: &str| -> Result<Option<f64>, FormatError> {
            let mut values = record
                .data_fields()
                .iter()
                .filter(|field| field.name() == name);
            let Some(field) = values.next() else {
                return Ok(None);
            };
            if values.next().is_some() {
                return Err(invalid(
                    DiagnosticKind::Inconsistent,
                    name,
                    "is repeated".to_owned(),
                ));
            }
            let value = field.value().trim().parse::<f64>().map_err(|_| {
                invalid(
                    DiagnosticKind::InvalidValue,
                    name,
                    format!("is not a number: `{}`", field.value()),
                )
            })?;
            Ok(Some(value))
        };
        let weight = number(SDF_WEIGHT_FIELD)?.unwrap_or(1.0);
        let mut member =
            EnsembleMember::new(record.model().positions().clone(), weight).map_err(|error| {
                invalid(
                    DiagnosticKind::InvalidValue,
                    SDF_WEIGHT_FIELD,
                    "is not a valid weight".to_owned(),
                )
                .with_source(error)
            })?;
        if let Some(energy) = number(SDF_ENERGY_FIELD)? {
            member = member
                .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
                .map_err(|error| {
                    invalid(
                        DiagnosticKind::InvalidValue,
                        SDF_ENERGY_FIELD,
                        "is not a valid energy".to_owned(),
                    )
                    .with_source(error)
                })?;
        }
        members.push(member);
    }
    Ensemble::from_items(topology, members).map_err(ensemble_error)
}
