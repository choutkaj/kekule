use std::cell::Cell;
use std::error::Error;

use kekule::core::MoleculeEditor;
use kekule::geometry::Point3;
use kekule::geometry::{PeriodicCell, Vector3};
use kekule::properties::{
    PropertyColumn, PropertyError, PropertyKey, PropertyValue, PropertyValueRef,
};
use kekule::smiles;
use kekule::structure::{
    Conformation, ConformationError, Model, ModelError, PositionError, Positions, RealizationError,
};
use kekule::topology::TopologyAtomIndex;
use kekule::units::{
    Quantity, ScaleValue, UnitError, ANGSTROM, CANONICAL_LENGTH_UNIT, DIMENSIONLESS, NANOMETER,
};

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
        .conformation_mut()
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
        model.conformation_mut().set_positions(Quantity::new(
            Alternating {
                calls: &calls,
                first: &invalid,
                later: &finite
            },
            CANONICAL_LENGTH_UNIT
        )),
        Err(ConformationError::Position(
            PositionError::NonFinitePosition { index: 0 }
        ))
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
    let row = TopologyAtomIndex::new;
    let mut conformation = Conformation::new(Positions::zeros(2));
    conformation
        .properties_mut()
        .atoms_mut()
        .insert(key.clone(), PropertyColumn::String(vec![Some(text), None]))
        .unwrap();
    let table = conformation.properties().atoms();
    let value = table.value_ref(&key, row(0)).unwrap().unwrap();
    let PropertyValueRef::String(text) = value else {
        panic!("expected string")
    };
    assert_eq!(text.as_ptr(), pointer);
    assert_eq!(
        value.to_value(),
        table.value(&key, row(0)).unwrap().unwrap()
    );
    assert_eq!(table.value_ref(&key, row(1)).unwrap(), None);
    assert_eq!(
        table
            .value_ref(&PropertyKey::new("missing").unwrap(), row(0))
            .unwrap(),
        None
    );
    assert!(matches!(
        table.value_ref(&key, row(2)),
        Err(PropertyError::InvalidIndex { len: 2, index: 2 })
    ));
    assert!(matches!(
        table.value_ref(&PropertyKey::new("missing").unwrap(), row(2)),
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
    let error = RealizationError::from(ConformationError::from(position_error.clone()));
    let conformation = error.source().unwrap();
    assert!(conformation.is::<ConformationError>());
    let position = conformation.source().unwrap();
    assert!(position.is::<PositionError>());
    assert!(position.source().unwrap().is::<UnitError>());
    let error = ModelError::from(position_error);
    assert!(error.source().unwrap().is::<ConformationError>());
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

#[test]
fn borrowed_quantity_conversion_retains_the_source() {
    let source = Quantity::new(vec![1.0, 2.0], NANOMETER);
    let converted = source.to_unit(ANGSTROM).unwrap();
    assert_eq!(source.value(), &[1.0, 2.0]);
    assert_eq!(source.unit(), NANOMETER);
    assert_eq!(converted.value(), &[10.0, 20.0]);
    assert_eq!(converted.unit(), ANGSTROM);
    assert_eq!(source.value_in(ANGSTROM).unwrap(), converted.into_value());
}

#[test]
fn consuming_quantity_conversions_support_values_that_cannot_be_cloned() {
    struct NonClone(f64);
    impl ScaleValue for NonClone {
        fn scaled(self, factor: f64) -> Self {
            Self(self.0 * factor)
        }
    }
    let converted = Quantity::new(NonClone(1.25), NANOMETER)
        .into_unit(ANGSTROM)
        .unwrap()
        .into_value();
    assert_eq!(converted.0, 12.5);

    let text = String::from("owned payload");
    let pointer = text.as_ptr();
    let moved = Quantity::new(text, ANGSTROM).into_value();
    assert_eq!(moved.as_ptr(), pointer);
}

#[test]
fn borrowed_document_conversion_keeps_source_and_consuming_interpretation_moves_storage() {
    let document = smiles::parse_str("CC").unwrap();
    let first = document.to_molecules().unwrap();
    let second = document.to_molecules().unwrap();
    assert_eq!(first, second);
    assert_eq!(document.source(), "CC");

    let interpretation = document.interpret().unwrap();
    let id = interpretation
        .molecule()
        .unwrap()
        .atom_ids()
        .next()
        .unwrap();
    let atom_storage = interpretation.molecule().unwrap().atom(id).unwrap() as *const _;
    let molecule = interpretation.into_molecule().unwrap();
    assert_eq!(molecule.atom(id).unwrap() as *const _, atom_storage);
    let editor = molecule.into_editor();
    assert_eq!(editor.atom(id).unwrap() as *const _, atom_storage);
}

#[test]
fn copy_view_materialization_retains_owner_and_creates_independent_geometry() {
    let molecule = smiles::to_molecules("CC").unwrap().pop().unwrap();
    let source = Model::from_molecule(
        molecule.clone(),
        &Positions::new(Quantity::new(
            vec![Point3::origin(), Point3::new(1.5, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    let view = source.as_model_view();
    let mut owned = view.to_model();
    assert_eq!(view.atom_count(), source.atom_count());
    assert!(std::sync::Arc::ptr_eq(
        &owned.shared_topology(),
        &source.shared_topology()
    ));
    assert_ne!(
        owned.positions().values().value().as_ptr(),
        source.positions().values().value().as_ptr()
    );
    let atom = source.topology().atom_ids()[0];
    owned
        .set_position(atom, Quantity::new(Point3::new(9.0, 0.0, 0.0), ANGSTROM))
        .unwrap();
    assert_eq!(view.position(atom).unwrap(), source.position(atom).unwrap());
    assert_ne!(
        owned.position(atom).unwrap(),
        source.position(atom).unwrap()
    );
}

#[test]
fn failed_publication_returns_the_editor_by_ownership_transfer() {
    let failed = MoleculeEditor::new().try_finish().unwrap_err();
    let mut editor = failed.into_editor();
    editor
        .append_molecule(&smiles::to_molecules("C").unwrap().pop().unwrap())
        .unwrap();
    assert_eq!(editor.finish().unwrap().atom_count(), 1);
}
