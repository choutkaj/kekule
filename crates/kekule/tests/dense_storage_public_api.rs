use kekule::geometry::Point3;
use kekule::structure::{AtomData, BondData, PositionError, Positions};
use kekule::units::{Quantity, ANGSTROM, DIMENSIONLESS};

#[test]
fn primitive_dense_storage_uses_usize_without_topology_imports() {
    let mut positions = Positions::new(Quantity::new(
        [Point3::origin(), Point3::new(1.0, 0.0, 0.0)],
        ANGSTROM,
    ))
    .unwrap();
    assert!((positions.position_at(1).unwrap().value().x - 0.1).abs() < 1.0e-15);
    positions
        .set_position_at(0, Quantity::new(Point3::new(2.0, 0.0, 0.0), ANGSTROM))
        .unwrap();
    assert_eq!(
        positions.position_at(2),
        Err(PositionError::InvalidIndex { index: 2 })
    );

    let mut atoms = AtomData::new(2);
    atoms.set_occupancy_at(1, Some(0.5)).unwrap();
    assert_eq!(atoms.occupancy_at(1).unwrap(), Some(0.5));

    let mut bonds = BondData::new(1);
    bonds
        .set_property_value_at("score", 0, Some(Quantity::new(3.0, DIMENSIONLESS)))
        .unwrap();
    assert_eq!(
        bonds.property_value_at("score", 0).unwrap(),
        Some(Quantity::new(3.0, DIMENSIONLESS))
    );
}

#[test]
fn dimensioned_data_emptiness_and_presence_are_unambiguous() {
    let atoms = AtomData::new(7);
    assert_eq!(atoms.len(), 7);
    assert!(!atoms.is_empty());
    assert!(!atoms.has_data());
    assert!(AtomData::new(0).is_empty());

    let bonds = BondData::new(5);
    assert_eq!(bonds.len(), 5);
    assert!(!bonds.is_empty());
    assert!(!bonds.has_data());
    assert!(BondData::new(0).is_empty());
}
