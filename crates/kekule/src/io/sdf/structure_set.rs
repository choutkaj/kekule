//! Explicit assembly of interpreted SDF records into one structure set.
use crate::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, TopologyMismatch,
};
use std::collections::BTreeMap;

use crate::properties::{PropertyColumn, PropertyKey};
use crate::structure::{Geometry, Order, StructureSet, StructureSetError};
use crate::units::{named_unit, Unit, DIMENSIONLESS};

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
/// every record carries one. The other data fields become the set's columns
/// ([`record_columns`]).
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
    for (key, column) in record_columns(&records) {
        set.columns_mut().insert(key, column).map_err(set_error)?;
    }
    Ok(set)
}

/// The columns the records' data fields state, one cell per record.
///
/// A data field becomes a column when its name is a property key, it is not
/// a place field, and no record has it twice; the records keep every field.
/// A record without the field, or with an empty value, has a missing cell.
/// The column holds integers when every value is one, reals when every
/// value is a finite number, booleans when every value is `true` or
/// `false`, and the trimmed text otherwise, so a field whose type differs
/// between records reads as text. A real field `<name>` whose records, and
/// only those, carry `<name>.unit` with one symbol of a unit constant holds
/// reals in that unit, and that unit field is no column of its own.
fn record_columns(records: &[&SdfRecordInterpretation]) -> Vec<(PropertyKey, PropertyColumn)> {
    let place = [
        SDF_STEP_FIELD,
        SDF_COORDINATE_FIELD,
        SDF_COORDINATE_UNIT_FIELD,
    ];
    // A unit field names its real field's unit even when the unit field's
    // own name is too long to be a key.
    let candidate = |name: &str| {
        PropertyKey::new(name).is_ok()
            || name
                .strip_suffix(".unit")
                .is_some_and(|base| PropertyKey::new(base).is_ok())
    };
    // Each candidate name's trimmed value per record; `None` marks a name
    // that some record repeats.
    let mut fields = BTreeMap::<&str, Option<Vec<Option<&str>>>>::new();
    for (index, record) in records.iter().enumerate() {
        for field in record.data_fields() {
            let name = field.name();
            if place.contains(&name) || !candidate(name) {
                continue;
            }
            let entry = fields
                .entry(name)
                .or_insert_with(|| Some(vec![None; records.len()]));
            if let Some(values) = entry {
                if values[index].is_some() {
                    *entry = None;
                } else {
                    values[index] = Some(field.value().trim());
                }
            }
        }
    }
    let fields = fields
        .into_iter()
        .filter_map(|(name, values)| Some((name, values?)))
        .map(|(name, values)| {
            let values = values
                .into_iter()
                .map(|value| value.filter(|value| !value.is_empty()))
                .collect::<Vec<_>>();
            (name, values)
        })
        .filter(|(_, values)| values.iter().any(Option::is_some))
        .collect::<BTreeMap<_, _>>();
    let mut units = BTreeMap::<&str, Unit>::new();
    for (name, values) in &fields {
        if let Some(unit) = record_unit(&fields, name, values) {
            units.insert(name, unit);
        }
    }
    fields
        .iter()
        .filter(|(name, _)| {
            !name
                .strip_suffix(".unit")
                .is_some_and(|base| units.contains_key(base))
        })
        .filter_map(|(name, values)| {
            // A unit field of no real field is a column only under a key.
            let key = PropertyKey::new(*name).ok()?;
            let column = units.get(name).map_or_else(
                || text_column(values),
                |&unit| PropertyColumn::Real {
                    unit,
                    values: values
                        .iter()
                        .map(|value| value.and_then(|value| value.parse().ok()))
                        .collect(),
                },
            );
            Some((key, column))
        })
        .collect()
}

/// The unit of the real field `name`: the one symbol its unit field holds on
/// exactly the records that have a value of `name`.
fn record_unit(
    fields: &BTreeMap<&str, Vec<Option<&str>>>,
    name: &str,
    values: &[Option<&str>],
) -> Option<Unit> {
    if !values.iter().flatten().all(|value| real(value).is_some()) {
        return None;
    }
    let symbols = fields.get(format!("{name}.unit").as_str())?;
    if symbols
        .iter()
        .zip(values)
        .any(|(symbol, value)| symbol.is_some() != value.is_some())
    {
        return None;
    }
    let mut present = symbols.iter().flatten();
    let symbol = *present.next()?;
    if present.any(|other| *other != symbol) {
        return None;
    }
    named_unit(symbol)
}

/// A finite number.
fn real(value: &str) -> Option<f64> {
    value.parse::<f64>().ok().filter(|value| value.is_finite())
}

/// The column of a field's values, of the first type every value has:
/// integer, real, boolean, then text.
fn text_column(values: &[Option<&str>]) -> PropertyColumn {
    fn all<T>(
        values: &[Option<&str>],
        parse: impl Fn(&str) -> Option<T>,
    ) -> Option<Vec<Option<T>>> {
        values
            .iter()
            .map(|value| {
                value
                    .map(&parse)
                    .map_or(Some(None), |parsed| parsed.map(Some))
            })
            .collect()
    }
    if let Some(values) = all(values, |value| value.parse::<i64>().ok()) {
        return PropertyColumn::Int(values);
    }
    if let Some(values) = all(values, real) {
        return PropertyColumn::Real {
            unit: DIMENSIONLESS,
            values,
        };
    }
    let boolean = |value: &str| match value {
        "true" => Some(true),
        "false" => Some(false),
        _ => None,
    };
    if let Some(values) = all(values, boolean) {
        return PropertyColumn::Bool(values);
    }
    PropertyColumn::String(
        values
            .iter()
            .map(|value| value.map(str::to_owned))
            .collect(),
    )
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
