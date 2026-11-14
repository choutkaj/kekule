//! Explicit assembly of interpreted SDF records into one structure set.
use crate::diagnostics::TopologyMismatch;
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::structure::{Geometry, Order, StructureSet, StructureSetError};
use crate::units::named_unit;

use super::document::SdfRecordInterpretation;
use super::write::{SDF_COORDINATE_FIELD, SDF_COORDINATE_UNIT_FIELD, SDF_STEP_FIELD};

fn set_error(error: StructureSetError) -> FormatError {
    FormatError::interpret(
        Format::Sdf,
        DiagnosticKind::Inconsistent,
        "The records do not form a structure set.",
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

/// Assembles interpreted SDF records, in order, into one structure set on the
/// first record's topology, record `i` becoming geometry `i`. Every record must
/// have the same topology layout
/// ([`Topology::same_layout`](crate::topology::Topology::same_layout)).
///
/// The geometries are ordered by the records' [`SDF_STEP_FIELD`] or, in the
/// unit [`SDF_COORDINATE_UNIT_FIELD`] names, [`SDF_COORDINATE_FIELD`], when
/// every record carries one.
pub(crate) fn sdf_structure_set<'a>(
    records: impl IntoIterator<Item = &'a SdfRecordInterpretation>,
) -> Result<StructureSet, FormatError> {
    let records = records.into_iter().collect::<Vec<_>>();
    let topology = records
        .first()
        .ok_or_else(|| set_error(StructureSetError::EmptySource))?
        .structure()
        .shared_topology();
    let mut geometries = Vec::with_capacity(records.len());
    for (index, record) in records.iter().enumerate() {
        if !record.structure().topology().same_layout(&topology) {
            return Err(set_error(StructureSetError::TopologyMismatch(
                TopologyMismatch::Structure,
            ))
            .in_record(index + 1));
        }
        geometries.push(Geometry::new(record.structure().positions().clone()));
    }
    let mut set = StructureSet::from_geometries(topology, geometries).map_err(set_error)?;
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
                            "differs from that of record 1.",
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
