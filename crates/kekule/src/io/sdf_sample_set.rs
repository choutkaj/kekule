//! Explicit assembly of interpreted SDF records into one sample set.
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::structure::{
    Distribution, Members, Order, Partition, Sample, SampleSet, SampleSetError,
};
use crate::units::{named_unit, Quantity, KELVIN, KILOJOULE_PER_MOLE};

use super::sdf_document::SdfRecordInterpretation;
use super::sdf_write::{
    SDF_COORDINATE_FIELD, SDF_COORDINATE_UNIT_FIELD, SDF_ENERGY_FIELD, SDF_REPRESENTATIVE_FIELD,
    SDF_STATE_FIELD, SDF_STEP_FIELD, SDF_TEMPERATURE_FIELD, SDF_WEIGHT_FIELD,
};

fn set_error(error: SampleSetError) -> FormatError {
    FormatError::interpret(
        Format::Sdf,
        DiagnosticKind::Inconsistent,
        "The records do not form a sample set.",
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

/// The trimmed value of the data field `name` of the 1-based record
/// `number`, if it has one; a second field of that name fails.
fn field<'r>(
    record: &'r SdfRecordInterpretation,
    number: usize,
    name: &str,
) -> Result<Option<&'r str>, FormatError> {
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
}

/// The ensemble data fields of one record.
struct RecordFields {
    weight: Option<f64>,
    state: Option<usize>,
    representative: bool,
    energy: Option<f64>,
    temperature: Option<f64>,
}

impl RecordFields {
    fn read(record: &SdfRecordInterpretation, number: usize) -> Result<Self, FormatError> {
        let text = |name: &str| field(record, number, name);
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
            weight: real(SDF_WEIGHT_FIELD)?,
            state,
            representative,
            energy: real(SDF_ENERGY_FIELD)?,
            temperature: real(SDF_TEMPERATURE_FIELD)?,
        })
    }

    /// Whether the record has any ensemble data field.
    fn any(&self) -> bool {
        self.weight.is_some()
            || self.state.is_some()
            || self.representative
            || self.energy.is_some()
            || self.temperature.is_some()
    }
}

/// Assembles interpreted SDF records, in order, into one sample set on the
/// first record's topology, record `i` becoming sample `i`. Every record must
/// have the same topology layout
/// ([`Topology::same_layout`](crate::topology::Topology::same_layout)).
///
/// The set has a distribution when a record has a distribution data field
/// ([`SDF_WEIGHT_FIELD`], [`SDF_STATE_FIELD`], [`SDF_REPRESENTATIVE_FIELD`],
/// [`SDF_ENERGY_FIELD`], or [`SDF_TEMPERATURE_FIELD`]), and none otherwise.
/// The distribution has weights only when a record has [`SDF_WEIGHT_FIELD`],
/// each sample's weight being its record's (1.0 when absent). Records
/// without [`SDF_STATE_FIELD`] are one state each, with the energy of
/// [`SDF_ENERGY_FIELD`] in kJ/mol (none when absent). Records with it, as
/// the SDF writer emits for samples of a distribution, are grouped by it:
/// states are numbered in ascending order of its values, each represented
/// by the one record marked with [`SDF_REPRESENTATIVE_FIELD`], which alone
/// may carry the state's energy. [`SDF_TEMPERATURE_FIELD`], in kelvin, must
/// be absent from every record or the same on every record.
///
/// The samples are ordered by the records' [`SDF_STEP_FIELD`] or, in the
/// unit [`SDF_COORDINATE_UNIT_FIELD`] names, [`SDF_COORDINATE_FIELD`], when
/// every record carries one.
pub(crate) fn sdf_sample_set<'a>(
    records: impl IntoIterator<Item = &'a SdfRecordInterpretation>,
) -> Result<SampleSet, FormatError> {
    let records = records.into_iter().collect::<Vec<_>>();
    let topology = records
        .first()
        .ok_or_else(|| set_error(SampleSetError::EmptySource))?
        .model()
        .shared_topology();
    let mut samples = Vec::with_capacity(records.len());
    let mut weights = Vec::with_capacity(records.len());
    let mut fields = Vec::with_capacity(records.len());
    for (index, record) in records.iter().enumerate() {
        let number = index + 1;
        if !record.model().topology().same_layout(&topology) {
            return Err(set_error(SampleSetError::TopologyMismatch).in_record(number));
        }
        let read = RecordFields::read(record, number)?;
        let weight = read.weight.unwrap_or(1.0);
        if !(weight.is_finite() && weight > 0.0) {
            return Err(invalid(
                DiagnosticKind::InvalidValue,
                number,
                SDF_WEIGHT_FIELD,
                "holds an invalid weight.",
            ));
        }
        samples.push(Sample::new(record.model().positions().clone()));
        weights.push(weight);
        fields.push(read);
    }
    let mut set = SampleSet::from_samples(topology, samples).map_err(set_error)?;
    if fields.iter().any(RecordFields::any) {
        let mut distribution =
            Distribution::new(Members::All, partition(&fields)?).map_err(set_error)?;
        if fields.iter().any(|read| read.weight.is_some()) {
            distribution = distribution.with_weights(weights).map_err(set_error)?;
        }
        let temperature = temperature(&fields)?;
        set.set_distribution(Some(distribution))
            .map_err(set_error)?;
        let mut distribution = set
            .distribution_mut()
            .expect("the distribution was just installed");
        for (index, read) in fields.iter().enumerate() {
            let Some(energy) = read.energy else {
                continue;
            };
            let number = index + 1;
            let partition = distribution.partition();
            let state = partition.assignments()[index];
            if partition.representatives()[state] != index {
                return Err(invalid(
                    DiagnosticKind::Inconsistent,
                    number,
                    SDF_ENERGY_FIELD,
                    "is on a record that does not represent its state.",
                ));
            }
            distribution
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
        distribution
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
    }
    let order = record_order(&records)?;
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

/// The step or coordinate order the records' fields state: none when no
/// record has a step or a coordinate, otherwise every record's. The
/// coordinate unit is named by the symbol of a unit constant.
fn record_order(records: &[&SdfRecordInterpretation]) -> Result<Option<Order>, FormatError> {
    // Every record's place fields are read first, so a duplicate fails
    // wherever it is.
    let names = [
        SDF_STEP_FIELD,
        SDF_COORDINATE_FIELD,
        SDF_COORDINATE_UNIT_FIELD,
    ];
    let places = records
        .iter()
        .enumerate()
        .map(|(index, record)| {
            names
                .into_iter()
                .map(|name| Ok((name, field(record, index + 1, name)?)))
                .collect::<Result<Vec<_>, FormatError>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let text = |index: usize, name: &str| {
        places[index]
            .iter()
            .find_map(|(known, value)| (*known == name).then_some(*value))
            .flatten()
    };
    let has = |name: &str| (0..records.len()).find(|index| text(*index, name).is_some());
    let required = |index: usize, name: &str| {
        text(index, name).ok_or_else(|| {
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
            let unit = named_unit(symbol).ok_or_else(|| {
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
fn partition(fields: &[RecordFields]) -> Result<Partition, FormatError> {
    let stated = fields.iter().filter(|read| read.state.is_some()).count();
    if stated == 0 {
        return Ok(Partition::state_per_sample(fields.len()));
    }
    if let Some(missing) = fields.iter().position(|read| read.state.is_none()) {
        return Err(invalid(
            DiagnosticKind::Inconsistent,
            missing + 1,
            SDF_STATE_FIELD,
            "is missing. Other records have it.",
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
    Ok(Partition::new(assignments, representatives)
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
