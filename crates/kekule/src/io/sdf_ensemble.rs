//! Explicit assembly of interpreted SDF records into one ensemble.
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::structure::{Ensemble, EnsembleError, EnsemblePartition, EnsembleSample};
use crate::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE};

use super::sdf_document::SdfRecordInterpretation;
use super::sdf_write::{
    SDF_ENERGY_FIELD, SDF_REPRESENTATIVE_FIELD, SDF_STATE_FIELD, SDF_TEMPERATURE_FIELD,
    SDF_WEIGHT_FIELD,
};

fn ensemble_error(error: EnsembleError) -> FormatError {
    FormatError::interpret(
        Format::Sdf,
        DiagnosticKind::Inconsistent,
        "SDF records form no ensemble",
    )
    .with_source(error)
}

/// A data field failure of the 1-based record `record`.
fn invalid(kind: DiagnosticKind, record: usize, name: &str, message: &str) -> FormatError {
    FormatError::interpret(Format::Sdf, kind, format!("data field `{name}` {message}"))
        .in_record(record)
        .with_detail(FormatErrorDetail::DataField {
            name: name.to_owned(),
        })
}

/// The ensemble data fields of one record.
struct RecordFields {
    weight: f64,
    state: Option<usize>,
    representative: bool,
    energy: Option<f64>,
    temperature: Option<f64>,
}

impl RecordFields {
    fn read(record: &SdfRecordInterpretation, number: usize) -> Result<Self, FormatError> {
        let text = |name: &str| -> Result<Option<&str>, FormatError> {
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
                    number,
                    name,
                    "is repeated",
                ));
            }
            Ok(Some(field.value().trim()))
        };
        let real = |name: &str| -> Result<Option<f64>, FormatError> {
            text(name)?
                .map(|value| {
                    value.parse::<f64>().map_err(|_| {
                        invalid(
                            DiagnosticKind::InvalidValue,
                            number,
                            name,
                            &format!("is not a number: `{value}`"),
                        )
                    })
                })
                .transpose()
        };
        let state = text(SDF_STATE_FIELD)?
            .map(|value| {
                value.parse::<usize>().map_err(|_| {
                    invalid(
                        DiagnosticKind::InvalidValue,
                        number,
                        SDF_STATE_FIELD,
                        &format!("is not a state index: `{value}`"),
                    )
                })
            })
            .transpose()?;
        let representative = match text(SDF_REPRESENTATIVE_FIELD)? {
            None => false,
            Some("true") => true,
            Some(value) => {
                return Err(invalid(
                    DiagnosticKind::InvalidValue,
                    number,
                    SDF_REPRESENTATIVE_FIELD,
                    &format!("is not `true`: `{value}`"),
                ))
            }
        };
        Ok(Self {
            weight: real(SDF_WEIGHT_FIELD)?.unwrap_or(1.0),
            state,
            representative,
            energy: real(SDF_ENERGY_FIELD)?,
            temperature: real(SDF_TEMPERATURE_FIELD)?,
        })
    }
}

/// Assembles interpreted SDF records, in order, into one ensemble on the
/// first record's topology, record `i` becoming sample `i`. Every record must
/// have the same topology layout
/// ([`Topology::same_layout`](crate::topology::Topology::same_layout)).
///
/// Each sample takes its relative weight from the record's
/// [`SDF_WEIGHT_FIELD`] (1.0 when absent). Records without
/// [`SDF_STATE_FIELD`] are one state each, with the energy of
/// [`SDF_ENERGY_FIELD`] in kJ/mol (none when absent). Records with it, as the
/// SDF writer emits for ensemble samples, are grouped by it: states are
/// numbered in ascending order of its values, each represented by the one
/// record marked with [`SDF_REPRESENTATIVE_FIELD`], which alone may carry
/// the state's energy. [`SDF_TEMPERATURE_FIELD`], in kelvin, must be absent
/// from every record or the same on every record.
pub(crate) fn sdf_ensemble<'a>(
    records: impl IntoIterator<Item = &'a SdfRecordInterpretation>,
) -> Result<Ensemble, FormatError> {
    let records = records.into_iter().collect::<Vec<_>>();
    let topology = records
        .first()
        .ok_or_else(|| ensemble_error(EnsembleError::EmptySource))?
        .model()
        .shared_topology();
    let mut samples = Vec::with_capacity(records.len());
    let mut fields = Vec::with_capacity(records.len());
    for (index, record) in records.iter().enumerate() {
        let number = index + 1;
        if !record.model().topology().same_layout(&topology) {
            return Err(ensemble_error(EnsembleError::TopologyMismatch).in_record(number));
        }
        let read = RecordFields::read(record, number)?;
        let sample = EnsembleSample::new(record.model().positions().clone(), read.weight).map_err(
            |error| {
                invalid(
                    DiagnosticKind::InvalidValue,
                    number,
                    SDF_WEIGHT_FIELD,
                    "is not a valid weight",
                )
                .with_source(error)
            },
        )?;
        samples.push(sample);
        fields.push(read);
    }
    let partition = partition(&fields)?;
    let temperature = temperature(&fields)?;
    let mut ensemble =
        Ensemble::from_samples(topology, samples, partition).map_err(ensemble_error)?;
    for (index, read) in fields.iter().enumerate() {
        let Some(energy) = read.energy else {
            continue;
        };
        let number = index + 1;
        let state = ensemble.partition().assignments()[index];
        if ensemble.partition().representatives()[state] != index {
            return Err(invalid(
                DiagnosticKind::Inconsistent,
                number,
                SDF_ENERGY_FIELD,
                "belongs to the record of a state's representative",
            ));
        }
        ensemble
            .state_mut(state)
            .expect("every sample has a state")
            .set_energy(Some(Quantity::new(energy, KILOJOULE_PER_MOLE)))
            .map_err(|error| {
                invalid(
                    DiagnosticKind::InvalidValue,
                    number,
                    SDF_ENERGY_FIELD,
                    "is not a valid energy",
                )
                .with_source(error)
            })?;
    }
    ensemble
        .set_temperature(temperature.map(|kelvin| Quantity::new(kelvin, KELVIN)))
        .map_err(|error| {
            invalid(
                DiagnosticKind::InvalidValue,
                1,
                SDF_TEMPERATURE_FIELD,
                "is not a valid temperature",
            )
            .with_source(error)
        })?;
    Ok(ensemble)
}

/// One state per record without state fields; otherwise the states the
/// fields name, numbered in ascending order of their values.
fn partition(fields: &[RecordFields]) -> Result<EnsemblePartition, FormatError> {
    let stated = fields.iter().filter(|read| read.state.is_some()).count();
    if stated == 0 {
        return Ok(EnsemblePartition::state_per_sample(fields.len()));
    }
    if let Some(missing) = fields.iter().position(|read| read.state.is_none()) {
        return Err(invalid(
            DiagnosticKind::Inconsistent,
            missing + 1,
            SDF_STATE_FIELD,
            "is missing, but other records of the ensemble have it",
        ));
    }
    let mut values = fields
        .iter()
        .filter_map(|read| read.state)
        .collect::<Vec<_>>();
    values.sort_unstable();
    values.dedup();
    let mut assignments = Vec::with_capacity(fields.len());
    let mut representatives = vec![None; values.len()];
    for (index, read) in fields.iter().enumerate() {
        let value = read.state.expect("every record has a state");
        let state = values
            .binary_search(&value)
            .expect("every state value is listed");
        assignments.push(state);
        if read.representative && representatives[state].replace(index).is_some() {
            return Err(invalid(
                DiagnosticKind::Inconsistent,
                index + 1,
                SDF_REPRESENTATIVE_FIELD,
                &format!("marks a second representative of state {value}"),
            ));
        }
    }
    let representatives = representatives
        .into_iter()
        .enumerate()
        .map(|(state, representative)| {
            representative.ok_or_else(|| {
                let first = assignments
                    .iter()
                    .position(|assigned| *assigned == state)
                    .expect("every state has a record");
                invalid(
                    DiagnosticKind::Inconsistent,
                    first + 1,
                    SDF_REPRESENTATIVE_FIELD,
                    &format!("marks no record of state {}", values[state]),
                )
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(EnsemblePartition::new(assignments, representatives)
        .expect("every state is numbered and represented by one of its records"))
}

/// The temperature every record states, if any does.
fn temperature(fields: &[RecordFields]) -> Result<Option<f64>, FormatError> {
    let first = fields.first().and_then(|read| read.temperature);
    let bits = first.map(f64::to_bits);
    if let Some(index) = fields
        .iter()
        .position(|read| read.temperature.map(f64::to_bits) != bits)
    {
        return Err(invalid(
            DiagnosticKind::Inconsistent,
            index + 1,
            SDF_TEMPERATURE_FIELD,
            "differs from the first record's",
        ));
    }
    Ok(first)
}
