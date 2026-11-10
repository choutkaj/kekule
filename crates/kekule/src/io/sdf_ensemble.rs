//! Explicit assembly of interpreted SDF records into one ensemble or
//! sample set.
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::structure::{
    Ensemble, EnsembleError, EnsemblePartition, EnsembleSample, Order, SampleSet,
};
use crate::units::{
    Quantity, Unit, ANGSTROM, BOHR, COULOMB, DEGREE, DIMENSIONLESS, ELEMENTARY_CHARGE, FEMTOSECOND,
    JOULE, KELVIN, KILOCALORIE, KILOCALORIE_PER_MOLE, KILOJOULE, KILOJOULE_PER_MOLE, METER, MOLE,
    NANOMETER, PICOSECOND, RADIAN, SECOND, SQUARE_ANGSTROM, SQUARE_NANOMETER,
};

use super::sdf_document::SdfRecordInterpretation;
use super::sdf_write::{
    SDF_COORDINATE_FIELD, SDF_COORDINATE_UNIT_FIELD, SDF_ENERGY_FIELD, SDF_REPRESENTATIVE_FIELD,
    SDF_STATE_FIELD, SDF_STEP_FIELD, SDF_TEMPERATURE_FIELD, SDF_WEIGHT_FIELD,
};

fn ensemble_error(error: EnsembleError) -> FormatError {
    FormatError::interpret(
        Format::Sdf,
        DiagnosticKind::Inconsistent,
        "The records do not form an ensemble.",
    )
    .with_source(error)
}

/// A data field failure of the 1-based record `record`; `predicate` follows
/// the field's name in a sentence, such as `appears twice.`
fn invalid(kind: DiagnosticKind, record: usize, name: &str, predicate: &str) -> FormatError {
    FormatError::interpret(
        Format::Sdf,
        kind,
        format!("The data field {name:?} {predicate}"),
    )
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
                    "appears twice.",
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
                            &format!("holds {value:?}, not a number."),
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
                        &format!("holds {value:?}, not a state index."),
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
                    &format!("holds {value:?}, not \"true\"."),
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
                    "holds an invalid weight.",
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
                "is on a record that does not represent its state.",
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
                    "holds an invalid energy.",
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
                "holds an invalid temperature.",
            )
            .with_source(error)
        })?;
    Ok(ensemble)
}

/// Assembles interpreted SDF records as [`sdf_ensemble`] does, into an
/// unordered sample set whose distribution is that ensemble's, and orders
/// the samples by the records' [`SDF_STEP_FIELD`] or, in the unit
/// [`SDF_COORDINATE_UNIT_FIELD`] names, [`SDF_COORDINATE_FIELD`], when
/// every record carries one.
pub(crate) fn sdf_sample_set<'a>(
    records: impl IntoIterator<Item = &'a SdfRecordInterpretation>,
) -> Result<SampleSet, FormatError> {
    let records = records.into_iter().collect::<Vec<_>>();
    let ensemble = sdf_ensemble(records.iter().copied())?;
    let order = record_order(&records)?;
    let mut set = SampleSet::from_ensemble(ensemble);
    set.set_order(order).map_err(|error| {
        FormatError::interpret(
            Format::Sdf,
            DiagnosticKind::InvalidValue,
            "The records' coordinates do not form an order.",
        )
        .with_source(error)
    })?;
    Ok(set)
}

/// The units a record's [`SDF_COORDINATE_UNIT_FIELD`] can name, by symbol.
const COORDINATE_UNITS: [Unit; 21] = [
    DIMENSIONLESS,
    METER,
    NANOMETER,
    ANGSTROM,
    BOHR,
    SQUARE_ANGSTROM,
    SQUARE_NANOMETER,
    SECOND,
    PICOSECOND,
    FEMTOSECOND,
    KELVIN,
    MOLE,
    COULOMB,
    ELEMENTARY_CHARGE,
    RADIAN,
    DEGREE,
    JOULE,
    KILOJOULE,
    KILOCALORIE,
    KILOJOULE_PER_MOLE,
    KILOCALORIE_PER_MOLE,
];

/// The step or coordinate order the records' fields state: none when no
/// record has a step or a coordinate, otherwise every record's.
fn record_order(records: &[&SdfRecordInterpretation]) -> Result<Option<Order>, FormatError> {
    let text = |record: &SdfRecordInterpretation, name: &str| {
        record
            .data_fields()
            .iter()
            .find(|field| field.name() == name)
            .map(|field| field.value().trim().to_owned())
    };
    let has = |name: &str| {
        records
            .iter()
            .position(|record| text(record, name).is_some())
    };
    let required = |index: usize, name: &str| {
        text(records[index], name).ok_or_else(|| {
            invalid(
                DiagnosticKind::Inconsistent,
                index + 1,
                name,
                "is missing. Other records have it.",
            )
        })
    };
    match (has(SDF_STEP_FIELD), has(SDF_COORDINATE_FIELD)) {
        (None, None) => Ok(None),
        (Some(_), Some(placed)) => Err(invalid(
            DiagnosticKind::Inconsistent,
            placed + 1,
            SDF_COORDINATE_FIELD,
            "is on a record of a set with steps.",
        )),
        (Some(_), None) => {
            let steps = (0..records.len())
                .map(|index| {
                    let value = required(index, SDF_STEP_FIELD)?;
                    value.parse::<u64>().map_err(|_| {
                        invalid(
                            DiagnosticKind::InvalidValue,
                            index + 1,
                            SDF_STEP_FIELD,
                            &format!("holds {value:?}, not a step."),
                        )
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(Some(Order::Step(steps)))
        }
        (None, Some(_)) => {
            let symbol = required(0, SDF_COORDINATE_UNIT_FIELD)?;
            let unit = COORDINATE_UNITS
                .into_iter()
                .find(|unit| unit.to_string() == symbol)
                .ok_or_else(|| {
                    invalid(
                        DiagnosticKind::InvalidValue,
                        1,
                        SDF_COORDINATE_UNIT_FIELD,
                        &format!("holds {symbol:?}, not a known unit."),
                    )
                })?;
            let values = (0..records.len())
                .map(|index| {
                    if required(index, SDF_COORDINATE_UNIT_FIELD)? != symbol {
                        return Err(invalid(
                            DiagnosticKind::Inconsistent,
                            index + 1,
                            SDF_COORDINATE_UNIT_FIELD,
                            "differs from that of the first record.",
                        ));
                    }
                    let value = required(index, SDF_COORDINATE_FIELD)?;
                    value.parse::<f64>().map_err(|_| {
                        invalid(
                            DiagnosticKind::InvalidValue,
                            index + 1,
                            SDF_COORDINATE_FIELD,
                            &format!("holds {value:?}, not a number."),
                        )
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(Some(Order::Coordinate { unit, values }))
        }
    }
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
            "is missing. Other records of the ensemble have it.",
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
                &format!("marks a second representative of state {value}."),
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
                    &format!("marks no record of state {}.", values[state]),
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
            "differs from that of the first record.",
        ));
    }
    Ok(first)
}
