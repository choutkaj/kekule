use super::*;
use crate::molecule::AtomId;
use crate::units::{ANGSTROM, DIMENSIONLESS, KELVIN, NANOMETER};

fn key(value: &str) -> PropertyKey {
    PropertyKey::new(value).unwrap()
}

#[test]
fn detached_structure_bond_rows_bind_once_and_then_must_match() {
    let score = key("score");
    let mut properties = GeometryProperties::new(3);
    properties
        .owner_mut()
        .insert(score.clone(), PropertyValue::Int(9))
        .unwrap();
    assert_eq!(properties.atoms().len(), 3);
    assert!(properties.bonds().is_empty());
    assert!(properties.validate_bond_rows(2).is_ok());
    properties.bind_bond_rows(2);
    GeometryPropertiesMut::new(&mut properties)
        .bonds_mut()
        .set_value(
            score.clone(),
            TopologyBondIndex::new(1),
            Some(PropertyValue::Int(3)),
        )
        .unwrap();
    assert_eq!(
        properties.validate_bond_rows(1),
        Err(PropertyError::LengthMismatch(LengthMismatch::new(
            Countable::Value,
            1,
            2
        )))
    );
    let projected = properties.project(&[2, 0], &[1]).unwrap();
    assert!(projected.owner().is_empty());
    assert_eq!(
        projected
            .bonds()
            .value(&score, TopologyBondIndex::new(0))
            .unwrap(),
        Some(PropertyValue::Int(3))
    );
    assert_eq!(properties.owner().get(&score), Some(&PropertyValue::Int(9)));
}

#[test]
fn typed_tables_address_rows_and_batch_updates_are_transactional() {
    let tag = key("tag");
    let mut properties = MoleculeProperties::new(3, 0);
    let mut guard = MoleculePropertiesMut::new(&mut properties, LiveRows::All, LiveRows::All);
    guard
        .atoms_mut()
        .set_values(
            tag.clone(),
            [
                (AtomId::new(0), Some(PropertyValue::Int(1))),
                (AtomId::new(2), Some(PropertyValue::Int(3))),
            ],
        )
        .unwrap();
    let before = guard.atoms().clone();
    assert!(matches!(
        guard.atoms_mut().set_values(
            tag.clone(),
            [
                (AtomId::new(1), Some(PropertyValue::Int(2))),
                (AtomId::new(3), Some(PropertyValue::Int(4))),
            ],
        ),
        Err(PropertyError::InvalidIndex(payload)) if payload == IndexOutOfRange::new(Countable::Row, 3, 3)
    ));
    assert_eq!(guard.atoms(), &before);
    assert_eq!(
        properties.atoms().value(&tag, AtomId::new(2)).unwrap(),
        Some(PropertyValue::Int(3))
    );
    assert!(!properties.atoms().row_has_data(AtomId::new(1)).unwrap());
    assert_eq!(properties.atoms().keys().collect::<Vec<_>>(), [&tag]);
}

#[test]
fn property_keys_use_one_conservative_validation_path() {
    for valid in ["x", "_private", "partial_charge", "force-field.v1", "a-2"] {
        assert_eq!(PropertyKey::new(valid).unwrap().as_str(), valid);
        assert_eq!(valid.parse::<PropertyKey>().unwrap().as_str(), valid);
    }
    for invalid in ["", "2bad", "has space", "slash/name", "unicode_µ"] {
        assert!(PropertyKey::new(invalid).is_err(), "accepted {invalid:?}");
    }
    assert!(PropertyKey::new("x".repeat(MAX_PROPERTY_KEY_LEN + 1)).is_err());
}

#[test]
fn scalar_values_cover_the_complete_domain_and_reject_non_finite_reals() {
    let mut properties = OwnerProperties::new();
    properties
        .insert(key("bool"), PropertyValue::Bool(true))
        .unwrap();
    properties
        .insert(key("int"), PropertyValue::Int(-7))
        .unwrap();
    properties
        .insert(
            key("real"),
            PropertyValue::Real {
                value: 2.5,
                unit: DIMENSIONLESS,
            },
        )
        .unwrap();
    properties
        .insert(key("string"), PropertyValue::String("value".into()))
        .unwrap();
    assert_eq!(properties.iter().count(), 4);
    assert_eq!(
        properties
            .iter()
            .map(|(key, _)| key.as_str())
            .collect::<Vec<_>>(),
        ["bool", "int", "real", "string"]
    );
    assert!(matches!(
        properties.insert(
            key("nan"),
            PropertyValue::Real {
                value: f64::NAN,
                unit: DIMENSIONLESS,
            },
        ),
        Err(PropertyError::NonFiniteValue { index: None })
    ));
}

#[test]
fn every_column_type_supports_missing_values_and_deterministic_iteration() {
    let mut table = RawPropertyTable::new(3);
    table
        .insert(
            key("z_string"),
            PropertyColumn::String(vec![None, Some("x".into()), None]),
        )
        .unwrap();
    table
        .insert(
            key("a_bool"),
            PropertyColumn::Bool(vec![Some(true), None, Some(false)]),
        )
        .unwrap();
    table
        .insert(key("m_int"), PropertyColumn::Int(vec![None, Some(3), None]))
        .unwrap();
    table
        .insert(
            key("r_real"),
            PropertyColumn::Real {
                unit: DIMENSIONLESS,
                values: vec![Some(1.0), None, Some(2.0)],
            },
        )
        .unwrap();
    assert_eq!(
        table
            .iter()
            .map(|(key, _)| key.as_str())
            .collect::<Vec<_>>(),
        ["a_bool", "m_int", "r_real", "z_string"]
    );
    assert_eq!(table.value(&key("m_int"), 0).unwrap(), None);
    assert_eq!(
        table.value(&key("m_int"), 1).unwrap(),
        Some(PropertyValue::Int(3))
    );
}

#[test]
fn table_updates_are_transactional_unit_aware_and_normalize_all_missing() {
    let mut table = RawPropertyTable::new(2);
    let length_error = table.insert(key("bad_len"), PropertyColumn::Bool(vec![Some(true)]));
    assert!(matches!(
        length_error,
        Err(PropertyError::LengthMismatch { .. })
    ));
    assert!(!table.has_data());

    table
        .insert(
            key("distance"),
            PropertyColumn::Real {
                unit: NANOMETER,
                values: vec![Some(1.0), None],
            },
        )
        .unwrap();
    table
        .insert(
            key("distance"),
            PropertyColumn::Real {
                unit: ANGSTROM,
                values: vec![Some(5.0), Some(10.0)],
            },
        )
        .unwrap();
    let Some(PropertyValue::Real { value, unit }) = table.value(&key("distance"), 0).unwrap()
    else {
        panic!("distance should be a real value");
    };
    assert_eq!(unit, NANOMETER);
    assert!((value - 0.5).abs() < 1.0e-12);

    let before = table.clone();
    assert!(matches!(
        table.insert(
            key("distance"),
            PropertyColumn::Real {
                unit: KELVIN,
                values: vec![Some(1.0), Some(2.0)],
            },
        ),
        Err(PropertyError::Unit(_))
    ));
    assert_eq!(table, before);
    assert!(matches!(
        table.insert(
            key("finite"),
            PropertyColumn::Real {
                unit: DIMENSIONLESS,
                values: vec![Some(f64::INFINITY), None],
            },
        ),
        Err(PropertyError::NonFiniteValue { index: Some(0) })
    ));

    table.set_value(key("distance"), 0, None).unwrap();
    table.set_value(key("distance"), 1, None).unwrap();
    assert!(table.get(&key("distance")).is_none());
    table
        .insert(key("missing"), PropertyColumn::Int(vec![None, None]))
        .unwrap();
    assert!(table.get(&key("missing")).is_none());
}

#[test]
fn failed_single_cell_updates_leave_the_existing_column_unchanged() {
    let mut table = RawPropertyTable::new(3);
    table
        .insert(
            key("distance"),
            PropertyColumn::Real {
                unit: NANOMETER,
                values: vec![Some(1.0), Some(2.0), None],
            },
        )
        .unwrap();

    for invalid in [
        PropertyValue::Int(4),
        PropertyValue::Real {
            value: 3.0,
            unit: KELVIN,
        },
        PropertyValue::Real {
            value: f64::INFINITY,
            unit: NANOMETER,
        },
    ] {
        let before = table.clone();
        assert!(table.set_value(key("distance"), 1, Some(invalid)).is_err());
        assert_eq!(table, before);
    }
}

#[test]
fn compatible_single_cell_updates_convert_without_changing_column_unit() {
    let mut table = RawPropertyTable::new(2);
    table
        .insert(
            key("distance"),
            PropertyColumn::Real {
                unit: NANOMETER,
                values: vec![Some(1.0), None],
            },
        )
        .unwrap();
    table
        .set_value(
            key("distance"),
            1,
            Some(PropertyValue::Real {
                value: 5.0,
                unit: ANGSTROM,
            }),
        )
        .unwrap();

    let Some(PropertyColumn::Real { unit, values }) = table.get(&key("distance")) else {
        panic!("distance should remain a real column");
    };
    assert_eq!(*unit, NANOMETER);
    assert_eq!(values[0], Some(1.0));
    assert!((values[1].unwrap() - 0.5).abs() < 1.0e-12);
}

#[test]
fn successful_single_cell_update_keeps_the_column_allocation() {
    let mut table = RawPropertyTable::new(4);
    table
        .insert(
            key("labels"),
            PropertyColumn::String(vec![Some("a".into()), None, None, None]),
        )
        .unwrap();
    let before = match table.get(&key("labels")).unwrap() {
        PropertyColumn::String(values) => values.as_ptr(),
        _ => unreachable!(),
    };

    table
        .set_value(key("labels"), 3, Some(PropertyValue::String("d".into())))
        .unwrap();

    let after = match table.get(&key("labels")).unwrap() {
        PropertyColumn::String(values) => values.as_ptr(),
        _ => unreachable!(),
    };
    assert_eq!(after, before);
}

#[test]
fn clearing_the_final_populated_cell_removes_the_column() {
    let mut table = RawPropertyTable::new(4);
    table
        .set_value(key("flag"), 2, Some(PropertyValue::Bool(true)))
        .unwrap();
    table.set_value(key("flag"), 2, None).unwrap();
    assert_eq!(table.get(&key("flag")), None);
    assert!(!table.has_data());
}

#[test]
fn checked_projection_preserves_columns_and_missing_cells() {
    let mut table = RawPropertyTable::new(3);
    table
        .insert(
            key("flag"),
            PropertyColumn::Bool(vec![Some(true), None, Some(false)]),
        )
        .unwrap();
    let selected = table.select_indices(&[2, 1]).unwrap();
    assert_eq!(selected.len(), 2);
    assert_eq!(
        selected.value(&key("flag"), 0).unwrap(),
        Some(PropertyValue::Bool(false))
    );
    assert_eq!(selected.value(&key("flag"), 1).unwrap(), None);
    assert!(matches!(
        table.select_indices(&[3]),
        Err(PropertyError::InvalidIndex { .. })
    ));
}

#[test]
fn optional_projection_is_column_wise_and_preserves_real_units() {
    let mut table = RawPropertyTable::new(3);
    table
        .insert(
            key("distance"),
            PropertyColumn::Real {
                unit: ANGSTROM,
                values: vec![Some(1.0), None, Some(3.0)],
            },
        )
        .unwrap();
    table
        .insert(
            key("label"),
            PropertyColumn::String(vec![Some("a".into()), None, Some("c".into())]),
        )
        .unwrap();

    let selected = table
        .select_optional_indices(&[Some(2), None, Some(0), Some(1)])
        .unwrap();
    assert_eq!(
        selected.get(&key("distance")),
        Some(&PropertyColumn::Real {
            unit: ANGSTROM,
            values: vec![Some(3.0), None, Some(1.0), None],
        })
    );
    assert_eq!(
        selected.get(&key("label")),
        Some(&PropertyColumn::String(vec![
            Some("c".into()),
            None,
            Some("a".into()),
            None,
        ]))
    );
    assert!(matches!(
        table.select_optional_indices(&[Some(3)]),
        Err(PropertyError::InvalidIndex { .. })
    ));
}

#[test]
fn generic_tables_do_not_reserve_structure_semantic_names() {
    let mut table = RawPropertyTable::new(1);
    table
        .set_value(key("occupancy"), 0, Some(PropertyValue::Int(7)))
        .unwrap();
    table
        .insert(
            key("b_factor"),
            PropertyColumn::String(vec![Some("generic".into())]),
        )
        .unwrap();
    assert_eq!(
        table.value(&key("occupancy"), 0).unwrap(),
        Some(PropertyValue::Int(7))
    );
    assert_eq!(
        table.value(&key("b_factor"), 0).unwrap(),
        Some(PropertyValue::String("generic".into()))
    );
}
