use std::cell::Cell;
use std::error::Error;

use kekule::geometry::Point3;
use kekule::geometry::{PeriodicCell, Vector3};
use kekule::properties::{
    PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue, PropertyValueRef,
};
use kekule::structure::{
    EnsembleError, EnsembleSliceError, Model, ModelError, ModelSliceError, PositionError, Positions,
};
use kekule::units::{Quantity, UnitError, ANGSTROM, CANONICAL_LENGTH_UNIT, DIMENSIONLESS};

struct Alternating<'a> {
    calls: &'a Cell<usize>,
    first: &'a [Point3],
    later: &'a [Point3],
}

impl AsRef<[Point3]> for Alternating<'_> {
    fn as_ref(&self) -> &[Point3] {
        let call = self.calls.get();
        self.calls.set(call + 1);
        if call == 0 {
            self.first
        } else {
            self.later
        }
    }
}

#[test]
fn bulk_positions_copy_the_exact_slice_that_was_validated() {
    let finite = [Point3::new(1.0, 2.0, 3.0)];
    let invalid = [Point3::new(f64::NAN, 0.0, 0.0)];
    let calls = Cell::new(0);
    let mut model = Model::new(
        kekule::smiles::to_topology("C").unwrap(),
        Positions::zeros(1),
    )
    .unwrap();
    model
        .set_positions(Quantity::new(
            Alternating {
                calls: &calls,
                first: &finite,
                later: &invalid,
            },
            CANONICAL_LENGTH_UNIT,
        ))
        .unwrap();
    assert_eq!(calls.get(), 1);
    assert_eq!(model.positions().values().value(), &finite);

    calls.set(0);
    let before = model.positions().clone();
    assert!(matches!(
        model.set_positions(Quantity::new(
            Alternating {
                calls: &calls,
                first: &invalid,
                later: &finite
            },
            CANONICAL_LENGTH_UNIT
        )),
        Err(PositionError::NonFinitePosition { index: 0 })
    ));
    assert_eq!(calls.get(), 1);
    assert_eq!(model.positions(), &before);
}

#[test]
fn position_vectors_transfer_their_allocation_and_convert_units() {
    let mut values = Vec::with_capacity(8);
    values.extend([Point3::new(10.0, 20.0, -30.0), Point3::origin()]);
    let pointer = values.as_ptr();
    let capacity = values.capacity();
    let positions = Positions::from_vec(Quantity::new(values, ANGSTROM)).unwrap();
    assert_eq!(positions.values().value().as_ptr(), pointer);
    let values = positions.into_values();
    assert_eq!(values.unit(), CANONICAL_LENGTH_UNIT);
    assert_eq!(values.value().as_ptr(), pointer);
    assert_eq!(values.value().capacity(), capacity);
    assert!((values.value()[0].x - 1.0).abs() < 1.0e-14);
    assert!((values.value()[0].y - 2.0).abs() < 1.0e-14);
    assert!((values.value()[0].z + 3.0).abs() < 1.0e-14);
    assert_eq!(values.value()[1], Point3::origin());
}

#[test]
fn owned_positions_reject_invalid_units_and_nonfinite_conversions() {
    assert!(matches!(
        Positions::from_vec(Quantity::new(vec![Point3::origin()], DIMENSIONLESS)),
        Err(PositionError::Unit(_))
    ));
    assert!(matches!(
        Positions::from_vec(Quantity::new(
            vec![Point3::new(f64::INFINITY, 0.0, 0.0)],
            ANGSTROM
        )),
        Err(PositionError::NonFinitePosition { index: 0 })
    ));
    let huge_unit =
        kekule::units::Unit::new(CANONICAL_LENGTH_UNIT.dimension(), 1.0e200, None).unwrap();
    assert!(matches!(
        Positions::from_vec(Quantity::new(
            vec![Point3::new(1.0e300, 0.0, 0.0)],
            huge_unit
        )),
        Err(PositionError::NonFinitePosition { index: 0 })
    ));
}

#[test]
fn borrowed_property_reads_retain_string_storage_and_missing_cells() {
    let key = PropertyKey::new("label").unwrap();
    let text = String::from("a property that should be borrowed");
    let pointer = text.as_ptr();
    let mut table = PropertyTable::new(2);
    table
        .insert(key.clone(), PropertyColumn::String(vec![Some(text), None]))
        .unwrap();
    let value = table.value_ref(&key, 0).unwrap().unwrap();
    let PropertyValueRef::String(text) = value else {
        panic!("expected string")
    };
    assert_eq!(text.as_ptr(), pointer);
    assert_eq!(value.to_value(), table.value(&key, 0).unwrap().unwrap());
    assert_eq!(table.value_ref(&key, 1).unwrap(), None);
    assert_eq!(
        table
            .value_ref(&PropertyKey::new("missing").unwrap(), 0)
            .unwrap(),
        None
    );
    assert!(matches!(
        table.value_ref(&key, 2),
        Err(PropertyError::InvalidIndex { len: 2, index: 2 })
    ));
    assert!(matches!(
        table.value_ref(&PropertyKey::new("missing").unwrap(), 2),
        Err(PropertyError::InvalidIndex { len: 2, index: 2 })
    ));
}

#[test]
fn scalar_and_column_views_preserve_each_property_kind() {
    let cases = [
        (
            PropertyValue::Bool(true),
            PropertyColumn::Bool(vec![Some(true)]),
        ),
        (PropertyValue::Int(-3), PropertyColumn::Int(vec![Some(-3)])),
        (
            PropertyValue::real(4.5, ANGSTROM).unwrap(),
            PropertyColumn::Real {
                unit: ANGSTROM,
                values: vec![Some(4.5)],
            },
        ),
        (
            PropertyValue::String("retained".to_owned()),
            PropertyColumn::String(vec![Some("retained".to_owned())]),
        ),
    ];
    for (value, column) in cases {
        assert_eq!(value.as_ref(), column.value_ref(0).unwrap().unwrap());
        assert_eq!(value.as_ref().to_value(), value);
    }
    let scalar = PropertyValue::String("owned".into());
    let PropertyValue::String(owned) = &scalar else {
        unreachable!()
    };
    let PropertyValueRef::String(borrowed) = scalar.as_ref() else {
        unreachable!()
    };
    assert_eq!(borrowed.as_ptr(), owned.as_ptr());
}

#[test]
fn structural_error_chains_retain_the_underlying_unit_failure() {
    let position_error =
        Positions::new(Quantity::new([Point3::origin()], DIMENSIONLESS)).unwrap_err();
    assert!(position_error.source().unwrap().is::<UnitError>());
    let error = EnsembleSliceError::from(EnsembleError::from(position_error.clone()));
    let member = error.source().unwrap();
    assert!(member.is::<EnsembleError>());
    let position = member.source().unwrap();
    assert!(position.is::<PositionError>());
    assert!(position.source().unwrap().is::<UnitError>());
    let error = ModelSliceError::from(ModelError::from(position_error));
    assert!(error.source().unwrap().is::<ModelError>());
    assert!(error
        .source()
        .unwrap()
        .source()
        .unwrap()
        .source()
        .unwrap()
        .is::<UnitError>());
    let cell_error = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 1.0, 1.0), DIMENSIONLESS),
        [true; 3],
    )
    .unwrap_err();
    assert!(cell_error.source().unwrap().is::<UnitError>());
}
