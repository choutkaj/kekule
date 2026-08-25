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
    assert!(!atom_data.is_empty());
    assert!(!atom_data.has_data());
    assert!(atom_data.occupancies().is_none());

    let bond_data = BondData::new(5);
    assert_eq!(bond_data.len(), 5);
    assert_eq!(bond_data.bond_count(), 5);
    assert!(!bond_data.is_empty());
    assert!(!bond_data.has_data());
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
    assert_eq!(atoms.occupancy_at(0).unwrap(), Some(0.5));
    assert!(atoms
        .b_factor_at(1)
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
        .set_property_value_at("display_width", 0, Some(Quantity::new(3.0, ANGSTROM)))
        .unwrap();
    assert!(bonds
        .property_value_at("display_width", 0)
        .unwrap()
        .unwrap()
        .is_close(&Quantity::new(0.3, NANOMETER), 1.0e-12, 1.0e-12)
        .unwrap());
}

#[test]
fn model_and_view_preserve_shared_topology_and_qualified_hierarchy_navigation() {
    let mut editor = MoleculeEditor::new();
    let carbon = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let oxygen = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    let bond = editor.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
    let molecule = editor.finish().unwrap();

    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    let first = builder.add_instance(definition).unwrap();
    let second = builder.add_instance(definition).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let first_residue = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    let second_residue = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(2), None, None)
        .unwrap();
    let carbon_site = builder
        .hierarchy_mut()
        .add_atom_site(
            first_residue,
            InstanceAtomId::new(first, carbon),
            SmcraAtomSiteMetadata::default(),
        )
        .unwrap();
    for (residue, instance, atom) in [
        (first_residue, first, oxygen),
        (second_residue, second, carbon),
        (second_residue, second, oxygen),
    ] {
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, atom),
                SmcraAtomSiteMetadata::default(),
            )
            .unwrap();
    }
    let topology = Arc::new(builder.build().unwrap());
    let mut model = Model::new(
        Arc::clone(&topology),
        positions([
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
        ]),
    )
    .unwrap();

    let first_chain = chain;
    let first_site = carbon_site;
    let first_atom = InstanceAtomId::new(first, carbon);
    let second_atom = InstanceAtomId::new(second, carbon);
    let first_bond = InstanceBondId::new(first, bond);
    let view = model.view();

    assert_eq!(model.chains().count(), 1);
    assert_eq!(model.residues().count(), 2);
    assert_eq!(model.atom_sites().count(), 4);
    assert_eq!(
        model.chains().map(InstanceChain::id).collect::<Vec<_>>(),
        view.chains().map(InstanceChain::id).collect::<Vec<_>>()
    );
    assert!(std::ptr::eq(
        model.chain(first_chain).unwrap().local(),
        view.chain(first_chain).unwrap().local()
    ));
    assert!(std::ptr::eq(
        model.residue(second_residue).unwrap().local(),
        view.residue(second_residue).unwrap().local()
    ));
    assert_eq!(model.atom_site(first_site).unwrap().atom(), first_atom);
    assert_eq!(view.atom_site(first_site).unwrap().atom(), first_atom);
    assert_eq!(
        model.atom_site_for_atom(first_atom).unwrap().unwrap().id(),
        first_site
    );
    assert_eq!(model.position(first_atom).unwrap().value().x, 1.0);
    assert_eq!(model.position(second_atom).unwrap().value().x, 3.0);
    assert_eq!(model.bond(first_bond).unwrap().order, BondOrder::Single);
    assert_eq!(view.bond(first_bond).unwrap().order, BondOrder::Single);
    assert_eq!(
        model.positions().values().value().as_ptr(),
        view.positions().values().value().as_ptr()
    );

    let cloned = model.clone();
    assert!(Arc::ptr_eq(
        &model.shared_topology(),
        &cloned.shared_topology()
    ));
    let shared = model.shared_topology();
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 11.0, 12.0), ANGSTROM),
        [true; 3],
    )
    .unwrap();
    model.set_cell(Some(cell));
    assert!(Arc::ptr_eq(&model.shared_topology(), &shared));
    assert_eq!(model.cell(), Some(&cell));
}

#[test]
fn atom_data_columns_names_units_and_updates_are_strongly_validated() {
    let mut data = AtomData::new(4);
    assert!(!data.is_empty());
    assert!(!data.has_data());
    assert!(AtomData::new(0).is_empty());

    data.set_occupancies([None; 4]).unwrap();
    data.set_b_factors(Quantity::new([None; 4], NANOMETER.powi(2)))
        .unwrap();
    assert!(!data.has_data());

    data.set_occupancies([Some(0.25), None, None, None])
        .unwrap();
    data.set_b_factors(Quantity::new(
        [None, Some(0.01), None, None],
        NANOMETER.powi(2),
    ))
    .unwrap();
    assert!(data.has_data());
    assert_eq!(data.occupancy_at(0).unwrap(), Some(0.25));
    assert!(data
        .b_factor_at(1)
        .unwrap()
        .unwrap()
        .is_close(&Quantity::new(1.0, SQUARE_ANGSTROM), 1.0e-12, 1.0e-12)
        .unwrap());
    assert_eq!(
        data.occupancies.as_ref().map(|column| column.unit),
        Some(DIMENSIONLESS)
    );
    assert_eq!(
        data.b_factors.as_ref().map(|column| column.unit),
        Some(SQUARE_ANGSTROM)
    );

    data.set_property(
        "partial_charge",
        Quantity::new([Some(-0.4), None, Some(0.2), Some(0.2)], DIMENSIONLESS),
    )
    .unwrap();
    assert_eq!(
        data.property_value_at("partial_charge", 0).unwrap(),
        Some(Quantity::new(-0.4, DIMENSIONLESS))
    );
    data.set_property_value_at("display_radius", 1, Some(Quantity::new(0.2, NANOMETER)))
        .unwrap();
    assert_eq!(
        data.properties().map(|(name, _)| name).collect::<Vec<_>>(),
        vec!["display_radius", "partial_charge"]
    );
    assert_eq!(
        data.property("display_radius").unwrap().unwrap().value(),
        &[None, Some(0.2), None, None]
    );
    assert_eq!(data.property_value_at("missing", 0).unwrap(), None);

    for reserved in ["occupancy", "Occupancy", "b_factor", "B_FACTOR"] {
        assert!(matches!(
            data.set_property(reserved, Quantity::new([None; 4], DIMENSIONLESS)),
            Err(AtomDataError::ReservedPropertyName { .. })
        ));
    }
    assert!(matches!(
        data.set_property("bad name", Quantity::new([None; 4], DIMENSIONLESS)),
        Err(AtomDataError::InvalidPropertyName { .. })
    ));
    assert_eq!(
        data.occupancy_at(4),
        Err(AtomDataError::InvalidIndex { index: 4 })
    );

    let before = data.clone();
    assert!(matches!(
        data.set_occupancies([Some(1.0)]),
        Err(AtomDataError::AtomCountMismatch { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_occupancy_at(0, Some(f64::NAN)),
        Err(AtomDataError::NonFiniteOccupancy { index: 0 })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_b_factor_at(1, Some(Quantity::new(1.0, KELVIN))),
        Err(AtomDataError::Unit(UnitError::IncompatibleUnits { .. }))
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property(
            "partial_charge",
            Quantity::new([Some(f64::INFINITY), None, None, None], DIMENSIONLESS)
        ),
        Err(AtomDataError::NonFinitePropertyValue { index: 0, .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property("partial_charge", Quantity::new([Some(1.0); 4], KELVIN)),
        Err(AtomDataError::PropertyUnit { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property("partial_charge", Quantity::new([Some(1.0)], DIMENSIONLESS)),
        Err(AtomDataError::PropertyValueCountMismatch { .. })
    ));
    assert_eq!(data, before);

    data.set_occupancy_at(0, None).unwrap();
    data.set_b_factor_at(1, None).unwrap();
    assert!(data.occupancies().is_none());
    assert!(data.b_factors().is_none());
    data.set_property_value_at("display_radius", 1, None)
        .unwrap();
    assert!(data.property("display_radius").unwrap().is_none());
    data.set_property("partial_charge", Quantity::new([None; 4], DIMENSIONLESS))
        .unwrap();
    assert!(!data.has_data());

    data.set_occupancies([Some(0.5), None, None, None]).unwrap();
    data.set_b_factors(Quantity::new(
        [Some(1.0), None, None, None],
        SQUARE_ANGSTROM,
    ))
    .unwrap();
    data.clear_occupancies();
    data.clear_b_factors();
    assert!(!data.has_data());
}

#[test]
fn bond_data_properties_are_dense_unit_aware_and_transactional() {
    let entropy_unit = KILOJOULE_PER_MOLE / KELVIN;
    let mut data = BondData::new(2);
    assert_eq!(data.len(), 2);
    assert!(!data.is_empty());
    assert!(!data.has_data());
    assert!(BondData::new(0).is_empty());

    data.set_property(
        "conformational_entropy",
        Quantity::new([Some(0.012), None], entropy_unit),
    )
    .unwrap();
    assert!(data.has_data());
    assert_eq!(
        data.property_value_at("conformational_entropy", 0).unwrap(),
        Some(Quantity::new(0.012, entropy_unit))
    );
    data.set_property_value_at("display_width", 0, Some(Quantity::new(0.2, NANOMETER)))
        .unwrap();
    data.set_property_value_at("display_width", 0, Some(Quantity::new(3.0, ANGSTROM)))
        .unwrap();
    assert!(data
        .property_value_at("display_width", 0)
        .unwrap()
        .unwrap()
        .is_close(&Quantity::new(0.3, NANOMETER), 1.0e-12, 1.0e-12)
        .unwrap());
    assert_eq!(
        data.property_value_at("display_width", 2),
        Err(BondDataError::InvalidIndex { index: 2 })
    );

    let before = data.clone();
    assert!(matches!(
        data.set_property(
            "conformational_entropy",
            Quantity::new([Some(f64::INFINITY), None], entropy_unit)
        ),
        Err(BondDataError::NonFinitePropertyValue { index: 0, .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property(
            "conformational_entropy",
            Quantity::new([Some(1.0), None], ANGSTROM)
        ),
        Err(BondDataError::PropertyUnit { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property(
            "conformational_entropy",
            Quantity::new([Some(1.0)], entropy_unit)
        ),
        Err(BondDataError::PropertyValueCountMismatch { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property("bad name", Quantity::new([None; 2], entropy_unit)),
        Err(BondDataError::InvalidPropertyName { .. })
    ));

    data.set_property_value_at("display_width", 0, None)
        .unwrap();
    assert!(data.property("display_width").unwrap().is_none());
    data.set_property(
        "conformational_entropy",
        Quantity::new([None; 2], entropy_unit),
    )
    .unwrap();
    assert!(!data.has_data());
    assert!(!data.remove_property("conformational_entropy").unwrap());
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
    let topology = two_bond_instances_topology();
    let atom = topology.atom_ids()[0];
    let mut first = EnsembleMember::new(
        positions([
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
        ]),
        topology.bond_count(),
    );
    first.set_weight(Some(1.0)).unwrap();
    first
        .atom_data_mut()
        .set_occupancy_at(0, Some(0.5))
        .unwrap();
    first
        .bond_data_mut()
        .set_property("score", Quantity::new([Some(2.0), None], DIMENSIONLESS))
        .unwrap();
    assert!(matches!(
        first.set_atom_data(AtomData::new(3)),
        Err(EnsembleError::AtomDataCountMismatch { .. })
    ));
    assert!(matches!(
        first.set_bond_data(BondData::new(1)),
        Err(EnsembleError::BondDataCountMismatch { .. })
    ));

    let mut second = EnsembleMember::new(
        positions([
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(6.0, 0.0, 0.0),
            Point3::new(7.0, 0.0, 0.0),
            Point3::new(8.0, 0.0, 0.0),
        ]),
        topology.bond_count(),
    );
    second.set_weight(Some(3.0)).unwrap();
    let mut ensemble = Ensemble::from_members(Arc::clone(&topology), [first, second]).unwrap();
    ensemble.normalize_weights().unwrap();
    assert_eq!(ensemble.member(0).unwrap().weight(), Some(0.25));
    assert_eq!(ensemble.member(1).unwrap().weight(), Some(0.75));
    assert_eq!(
        ensemble
            .views()
            .map(|view| view.positions().values().value()[0].x)
            .collect::<Vec<_>>(),
        vec![1.0, 5.0]
    );
    let first_view = ensemble.views().next().unwrap();
    assert_eq!(first_view.occupancy(atom).unwrap(), Some(0.5));
    assert_eq!(
        first_view
            .bond_data()
            .property_value_at("score", 0)
            .unwrap(),
        Some(Quantity::new(2.0, DIMENSIONLESS))
    );
    assert!(Arc::ptr_eq(&ensemble.shared_topology(), &topology));
    assert!(matches!(
        ensemble.push(EnsembleMember::new(
            positions([Point3::origin(); 3]),
            topology.bond_count()
        )),
        Err(EnsembleError::PositionCountMismatch { .. })
    ));
    assert!(matches!(
        ensemble.push(EnsembleMember::new(positions([Point3::origin(); 4]), 1)),
        Err(EnsembleError::BondDataCountMismatch { .. })
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

#[test]
fn model_and_ensemble_slices_reuse_one_subset_mapping_for_dense_state() {
    let mut editor = MoleculeEditor::new();
    let first = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let second = editor
        .add_atom(Atom::new(Element::from_symbol("N").unwrap()))
        .unwrap();
    let third = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    editor.add_bond(first, second, BondOrder::Single).unwrap();
    editor.add_bond(second, third, BondOrder::Double).unwrap();
    let molecule = editor.finish().unwrap();
    let mut builder = ModelBuilder::new();
    let instance = builder
        .add_molecule(
            &molecule,
            &positions([
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(3.0, 0.0, 0.0),
            ]),
        )
        .unwrap();
    let chain = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_chain("A", None)
        .unwrap();
    for (sequence, atom) in [(1, first), (2, second), (3, third)] {
        let residue = builder
            .topology_builder_mut()
            .hierarchy_mut()
            .add_residue(chain, "RES", Some(sequence), None, None)
            .unwrap();
        builder
            .topology_builder_mut()
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, atom),
                SmcraAtomSiteMetadata::default(),
            )
            .unwrap();
    }
    let mut model = builder.build().unwrap();
    model
        .atom_data_mut()
        .set_occupancies([Some(0.1), Some(0.2), Some(0.3)])
        .unwrap();
    model
        .bond_data_mut()
        .set_property(
            "score",
            Quantity::new([Some(10.0), Some(20.0)], DIMENSIONLESS),
        )
        .unwrap();
    let source_topology = model.shared_topology();
    let selection = AtomSelection::from_atoms(
        &source_topology,
        [
            InstanceAtomId::new(instance, first),
            InstanceAtomId::new(instance, second),
        ],
    )
    .unwrap();

    let sliced = model.slice(&selection).unwrap();
    assert_eq!(sliced.topology().instance_count(), 1);
    assert_eq!(sliced.topology().residues().count(), 2);
    assert_eq!(
        sliced.positions().values().value(),
        &[Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)]
    );
    assert_eq!(
        sliced.atom_data().occupancies(),
        Some(&[Some(0.1), Some(0.2)][..])
    );
    assert_eq!(
        sliced
            .bond_data()
            .property("score")
            .unwrap()
            .unwrap()
            .value(),
        &[Some(10.0)]
    );

    let ensemble = Ensemble::from_models(&[model.clone(), model]).unwrap();
    let sliced_ensemble = ensemble.slice(&selection).unwrap();
    assert_eq!(sliced_ensemble.len(), 2);
    assert_eq!(sliced_ensemble.topology().atom_count(), 2);
    assert_eq!(
        sliced_ensemble
            .member(1)
            .unwrap()
            .positions()
            .values()
            .value(),
        sliced.positions().values().value()
    );
}
