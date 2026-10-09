use kekule::geometry::Point3;
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{Conformation, PositionError, Positions};
use kekule::topology::TopologyAtomIndex;
use kekule::units::{Quantity, ANGSTROM, DIMENSIONLESS};

#[test]
fn public_dense_containers_are_topology_free() {
    let positions = Positions::new(Quantity::new([Point3::new(10.0, 0.0, 0.0)], ANGSTROM)).unwrap();
    assert_eq!(positions.len(), 1);
    assert!((positions.values().value()[0].x - 1.0).abs() < 1.0e-15);

    let mut conformation = Conformation::new(Positions::zeros(2));
    let mut properties = conformation.properties_mut();
    let mut atoms = properties.atoms_mut();
    let key = PropertyKey::new("score").unwrap();
    atoms
        .set_value(
            key.clone(),
            TopologyAtomIndex::new(1),
            Some(PropertyValue::Real {
                value: 0.5,
                unit: DIMENSIONLESS,
            }),
        )
        .unwrap();
    assert_eq!(atoms.value(&key, TopologyAtomIndex::new(0)).unwrap(), None);
    assert_eq!(
        atoms.value(&key, TopologyAtomIndex::new(1)).unwrap(),
        Some(PropertyValue::Real {
            value: 0.5,
            unit: DIMENSIONLESS
        })
    );
    assert!(atoms
        .insert(key, PropertyColumn::Int(vec![Some(1)]))
        .is_err());
}

#[test]
fn dense_position_projection_is_checked() {
    let positions = Positions::zeros(2);
    assert!(matches!(
        positions.select_indices(&[2]),
        Err(PositionError::InvalidIndex { index: 2 })
    ));
}
