use std::sync::Arc;

use super::*;
use crate::core::*;
use crate::geometry::*;
use crate::topology::*;
use crate::units::*;

fn one_atom_topology() -> Arc<Topology> {
    let mut editor = MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(definition).unwrap();
    Arc::new(builder.build().unwrap())
}

fn two_bond_instances_topology() -> Arc<Topology> {
    let mut editor = MoleculeEditor::new();
    let carbon = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let oxygen = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    editor.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    Arc::new(builder.build().unwrap())
}

fn positions(points: impl AsRef<[Point3]>) -> Positions {
    Positions::new(Quantity::new(points, ANGSTROM)).unwrap()
}

#[test]
fn positions_are_topology_free_unit_aware_and_transactional() {
    let mut values =
        Positions::new(Quantity::new([Point3::new(0.1, 0.0, 0.0)], NANOMETER)).unwrap();
    assert_eq!(values.len(), 1);
    assert_eq!(values.values().value()[0], Point3::new(1.0, 0.0, 0.0));
    let pointer = values.values().value().as_ptr();
    values
        .set_all(Quantity::new([Point3::new(2.0, 0.0, 0.0)], ANGSTROM))
        .unwrap();
    assert_eq!(values.values().value().as_ptr(), pointer);

    let before = values.clone();
    assert!(matches!(
        values.set_all(Quantity::new([Point3::new(f64::NAN, 0.0, 0.0)], ANGSTROM)),
        Err(PositionError::NonFinitePosition { .. })
    ));
    assert_eq!(values, before);
    assert!(matches!(
        Positions::new(Quantity::new([Point3::origin()], PICOSECOND)),
        Err(PositionError::Unit(UnitError::IncompatibleUnits { .. }))
    ));
}

#[test]
fn primitive_equality_compares_dense_state_not_topology_identity() {
    let first = one_atom_topology();
    let second = one_atom_topology();
    assert!(!Arc::ptr_eq(&first, &second));
    assert!(first.same_layout(&second));

    assert_eq!(
        positions([Point3::new(1.0, 2.0, 3.0)]),
        positions([Point3::new(1.0, 2.0, 3.0)])
    );
    assert_eq!(AtomData::new(1), AtomData::new(1));
    assert_eq!(BondData::new(0), BondData::new(0));
}

#[test]
fn empty_atom_and_bond_data_retain_logical_lengths() {
    let atom_data = AtomData::new(7);
    assert_eq!(atom_data.len(), 7);
    assert_eq!(atom_data.atom_count(), 7);
    assert!(atom_data.is_empty());
    assert!(atom_data.occupancies().is_none());

    let bond_data = BondData::new(5);
    assert_eq!(bond_data.len(), 5);
    assert_eq!(bond_data.bond_count(), 5);
    assert!(bond_data.is_empty());
}

#[test]
fn dense_atom_and_bond_columns_are_unit_aware_and_transactional() {
    let mut atoms = AtomData::new(2);
    atoms.set_occupancies([Some(0.5), None]).unwrap();
    atoms
        .set_b_factors(Quantity::new([None, Some(0.01)], NANOMETER.powi(2)))
        .unwrap();
    atoms
        .set_property(
            "partial_charge",
            Quantity::new([Some(-0.2), Some(0.2)], DIMENSIONLESS),
        )
        .unwrap();
    assert_eq!(
        atoms.occupancy_at(TopologyAtomIndex::new(0)).unwrap(),
        Some(0.5)
    );
    assert!(atoms
        .b_factor_at(TopologyAtomIndex::new(1))
        .unwrap()
        .unwrap()
        .is_close(&Quantity::new(1.0, SQUARE_ANGSTROM), 1.0e-12, 1.0e-12,)
        .unwrap());
    let before = atoms.clone();
    assert!(matches!(
        atoms.set_property(
            "partial_charge",
            Quantity::new([Some(f64::NAN), None], DIMENSIONLESS)
        ),
        Err(AtomDataError::NonFinitePropertyValue { .. })
    ));
    assert_eq!(atoms, before);

    let mut bonds = BondData::new(1);
    bonds
        .set_property("display_width", Quantity::new([Some(0.2)], NANOMETER))
        .unwrap();
    bonds
        .set_property_value_at(
            "display_width",
            TopologyBondIndex::new(0),
            Some(Quantity::new(3.0, ANGSTROM)),
        )
        .unwrap();
    assert!(bonds
        .property_value_at("display_width", TopologyBondIndex::new(0))
        .unwrap()
        .unwrap()
        .is_close(&Quantity::new(0.3, NANOMETER), 1.0e-12, 1.0e-12)
        .unwrap());
}

#[test]
fn model_rejects_every_dense_dimension_mismatch() {
    let topology = two_bond_instances_topology();
    assert_eq!(
        Model::new(Arc::clone(&topology), positions([Point3::origin(); 3])),
        Err(ModelError::PositionCountMismatch {
            expected: 4,
            actual: 3
        })
    );
    assert!(matches!(
        Model::with_data(
            Arc::clone(&topology),
            positions([Point3::origin(); 4]),
            None,
            AtomData::new(3),
            BondData::new(2),
        ),
        Err(ModelError::AtomDataCountMismatch {
            expected: 4,
            actual: 3
        })
    ));
    assert!(matches!(
        Model::with_data(
            Arc::clone(&topology),
            positions([Point3::origin(); 4]),
            None,
            AtomData::new(4),
            BondData::new(1),
        ),
        Err(ModelError::BondDataCountMismatch {
            expected: 2,
            actual: 1
        })
    ));
}

#[test]
fn model_and_view_resolve_semantic_atoms_into_dense_state() {
    let topology = one_atom_topology();
    let atom = topology.atom_ids()[0];
    let mut model = Model::new(
        Arc::clone(&topology),
        positions([Point3::new(1.0, 2.0, 3.0)]),
    )
    .unwrap();
    model.set_occupancy(atom, Some(0.75)).unwrap();
    model
        .set_b_factor(atom, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
        .unwrap();
    assert_eq!(model.position(atom).unwrap().value().x, 1.0);
    assert_eq!(model.occupancy(atom).unwrap(), Some(0.75));
    assert_eq!(
        model.view().b_factor(atom).unwrap(),
        Some(Quantity::new(12.5, SQUARE_ANGSTROM))
    );
    model
        .set_position(atom, Quantity::new(Point3::new(4.0, 0.0, 0.0), ANGSTROM))
        .unwrap();
    assert_eq!(model.position(atom).unwrap().value().x, 4.0);
    assert_eq!(
        model.set_atom_data(AtomData::new(2)),
        Err(ModelError::AtomDataCountMismatch {
            expected: 1,
            actual: 2
        })
    );
}

#[test]
fn arbitrary_model_views_validate_dimensions_without_copying() {
    let topology = one_atom_topology();
    let positions = positions([Point3::origin()]);
    let atom_data = AtomData::new(1);
    let bond_data = BondData::new(0);
    let view = ModelView::new(&topology, &positions, None, &atom_data, &bond_data).unwrap();
    assert_eq!(
        view.positions().values().value().as_ptr(),
        positions.values().value().as_ptr()
    );
    assert!(matches!(
        ModelView::new(&topology, &positions, None, &AtomData::new(2), &bond_data),
        Err(ModelError::AtomDataCountMismatch { .. })
    ));
}

#[test]
fn ensemble_validates_member_dimensions_and_preserves_weights_and_views() {
    let topology = one_atom_topology();
    let mut first = EnsembleMember::new(positions([Point3::origin()]), 0);
    first.set_weight(Some(1.0)).unwrap();
    first
        .atom_data_mut()
        .set_occupancy_at(TopologyAtomIndex::new(0), Some(0.5))
        .unwrap();
    let mut second = EnsembleMember::new(positions([Point3::new(2.0, 0.0, 0.0)]), 0);
    second.set_weight(Some(3.0)).unwrap();
    let mut ensemble = Ensemble::from_members(Arc::clone(&topology), [first, second]).unwrap();
    ensemble.normalize_weights().unwrap();
    assert_eq!(ensemble.member(0).unwrap().weight(), Some(0.25));
    assert_eq!(ensemble.member(1).unwrap().weight(), Some(0.75));
    assert_eq!(ensemble.views().count(), 2);
    assert!(Arc::ptr_eq(&ensemble.shared_topology(), &topology));
    assert!(matches!(
        ensemble.push(EnsembleMember::new(
            positions([Point3::origin(), Point3::origin()]),
            0
        )),
        Err(EnsembleError::PositionCountMismatch { .. })
    ));
}

#[test]
fn topology_subset_transforms_return_topology_directly_and_preserve_noops() {
    let topology = two_bond_instances_topology();
    let instances = topology.instances().map(|(id, _)| id).collect::<Vec<_>>();
    let retained =
        crate::topology::transform::retain_instances(&topology, instances.clone()).unwrap();
    assert!(Arc::ptr_eq(&retained, &topology));
    let removed_none =
        crate::topology::transform::remove_instances(&topology, std::iter::empty()).unwrap();
    assert!(Arc::ptr_eq(&removed_none, &topology));

    let subset = crate::topology::transform::retain_instances(&topology, [instances[1]]).unwrap();
    assert_eq!(subset.instance_count(), 1);
    assert_eq!(subset.atom_count(), 2);
    assert!(!Arc::ptr_eq(&subset, &topology));
}
