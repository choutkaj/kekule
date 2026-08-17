use std::sync::Arc;

use super::*;
use crate::bio::*;
use crate::core::*;
use crate::geometry::*;
use crate::small::*;
use crate::topology::*;
use crate::units::*;
fn one_atom_topology() -> Arc<Topology> {
    let mut graph = Molecule::new();
    graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .expect("atom identifier capacity");
    let molecule = SmallMolecule::from_graph(graph);
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_small_molecule_definition(&molecule).unwrap();
    builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .unwrap();
    Arc::new(builder.build().unwrap())
}

fn two_bond_instances_topology() -> (Arc<Topology>, MoleculeInstanceId, MoleculeInstanceId) {
    let mut graph = Molecule::builder();
    let carbon = graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let oxygen = graph
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    graph.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
    let molecule = SmallMolecule::from_graph(graph.build().unwrap());
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_small_molecule_definition(&molecule).unwrap();
    let first = builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .unwrap();
    let second = builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .unwrap();
    (Arc::new(builder.build().unwrap()), first, second)
}

fn positions(topology: &Arc<Topology>, x: f64) -> Positions {
    Positions::new(
        topology,
        Quantity::new(vec![Point3::new(x, 0.0, 0.0)], ANGSTROM),
    )
    .unwrap()
}

#[test]
fn positions_convert_units_reuse_storage_and_update_transactionally() {
    let topology = one_atom_topology();
    let mut positions = Positions::new(
        &topology,
        Quantity::new(vec![Point3::new(0.1, 0.0, 0.0)], NANOMETER),
    )
    .unwrap();
    assert_eq!(positions.values().value()[0], Point3::new(1.0, 0.0, 0.0));
    let pointer = positions.values().value().as_ptr();
    positions
        .set_all(
            &topology,
            Quantity::new(vec![Point3::new(2.0, 0.0, 0.0)], ANGSTROM),
        )
        .unwrap();
    assert_eq!(positions.values().value().as_ptr(), pointer);

    let before = positions.clone();
    assert!(matches!(
        positions.set_all(
            &topology,
            Quantity::new(vec![Point3::new(f64::NAN, 0.0, 0.0)], ANGSTROM)
        ),
        Err(PositionError::NonFinitePosition { .. })
    ));
    assert_eq!(positions, before);
}

#[test]
fn model_requires_exact_topology_and_views_do_not_copy_coordinates() {
    let topology = one_atom_topology();
    let independent = one_atom_topology();
    assert!(topology.same_layout(&independent));
    assert!(!Arc::ptr_eq(&topology, &independent));

    let wrong_positions = positions(&independent, 1.0);
    assert_eq!(
        Model::new(Arc::clone(&topology), wrong_positions),
        Err(ModelError::TopologyMismatch)
    );

    let mut model = Model::new(Arc::clone(&topology), positions(&topology, 1.0)).unwrap();
    assert!(model.atom_data().is_empty());
    let view = model.view();
    assert_eq!(
        view.positions().values().value().as_ptr(),
        model.positions().values().value().as_ptr()
    );
    let clone = model.clone();
    assert!(Arc::ptr_eq(
        &model.shared_topology(),
        &clone.shared_topology()
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
fn model_and_view_share_qualified_hierarchy_navigation_without_copying() {
    let mut macro_builder = MacroMolecule::builder();
    let atom = macro_builder
        .graph_mut()
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let chain = macro_builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = macro_builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    let site = macro_builder
        .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
        .unwrap();
    let macro_molecule = macro_builder.build().unwrap();

    let mut topology_builder = TopologyBuilder::new();
    let definition = topology_builder
        .add_macro_molecule_definition(&macro_molecule)
        .unwrap();
    let first = topology_builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .unwrap();
    let second = topology_builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .unwrap();
    let topology = Arc::new(topology_builder.build().unwrap());
    let positions = Positions::new(
        &topology,
        Quantity::new(
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            ANGSTROM,
        ),
    )
    .unwrap();
    let model = Model::new(Arc::clone(&topology), positions).unwrap();
    let view = model.view();
    let first_chain = InstanceChainId::new(first, chain);
    let second_residue = InstanceResidueId::new(second, residue);
    let first_site = InstanceAtomSiteId::new(first, site);
    let first_atom = InstanceAtomId::new(first, atom);

    assert_eq!(
        model.chains().map(InstanceChain::id).collect::<Vec<_>>(),
        view.chains().map(InstanceChain::id).collect::<Vec<_>>()
    );
    assert_eq!(
        model
            .residues()
            .map(InstanceResidue::id)
            .collect::<Vec<_>>(),
        view.residues().map(InstanceResidue::id).collect::<Vec<_>>()
    );
    assert_eq!(
        model
            .atom_sites()
            .map(InstanceAtomSite::id)
            .collect::<Vec<_>>(),
        view.atom_sites()
            .map(InstanceAtomSite::id)
            .collect::<Vec<_>>()
    );
    assert!(std::ptr::eq(
        model.chain(first_chain).unwrap().local(),
        view.chain(first_chain).unwrap().local()
    ));
    assert!(std::ptr::eq(
        model.residue(second_residue).unwrap().local(),
        view.residue(second_residue).unwrap().local()
    ));
    assert!(std::ptr::eq(
        model.atom_site(first_site).unwrap().local(),
        view.atom_site(first_site).unwrap().local()
    ));
    assert_eq!(
        model.atom_site_for_atom(first_atom).unwrap().unwrap().id(),
        view.atom_site_for_atom(first_atom).unwrap().unwrap().id()
    );
    let model_site = model
        .chain(first_chain)
        .unwrap()
        .residues()
        .next()
        .unwrap()
        .atom_sites()
        .next()
        .unwrap();
    let view_site = view
        .chain(first_chain)
        .unwrap()
        .residues()
        .next()
        .unwrap()
        .atom_sites()
        .next()
        .unwrap();
    assert_eq!(model_site.id(), first_site);
    assert_eq!(view_site.id(), first_site);
    assert_eq!(model_site.atom(), first_atom);
    assert_eq!(view_site.atom(), first_atom);
    assert_eq!(
        model.hierarchy(first).unwrap().unwrap().molecule(),
        view.hierarchy(first).unwrap().unwrap().molecule()
    );
    assert_eq!(model.atoms().count(), view.atoms().count());
    assert_eq!(model.instances().count(), view.instances().count());
    assert_eq!(
        model.positions().values().value().as_ptr(),
        view.positions().values().value().as_ptr()
    );
}

#[test]
fn ensembles_validate_shared_topology_atom_data_order_and_weights() {
    let topology = one_atom_topology();
    let atom = topology.atom_ids()[0];
    let mut first = EnsembleMember::new(positions(&topology, 1.0));
    first.set_weight(Some(1.0)).unwrap();
    first
        .atom_data_mut()
        .set_occupancy(&topology, atom, Some(0.5))
        .unwrap();

    let mut second = EnsembleMember::new(positions(&topology, 2.0));
    second.set_weight(Some(3.0)).unwrap();
    let mut invalid_weight = EnsembleMember::new(positions(&topology, 3.0));
    assert_eq!(
        invalid_weight.set_weight(Some(-1.0)),
        Err(EnsembleError::InvalidWeight)
    );
    let mut ensemble = Ensemble::from_members(Arc::clone(&topology), [first, second]).unwrap();
    ensemble.normalize_weights().unwrap();
    assert_eq!(ensemble.member(0).unwrap().weight(), Some(0.25));
    assert_eq!(ensemble.member(1).unwrap().weight(), Some(0.75));
    assert_eq!(
        ensemble
            .views()
            .map(|view| view.positions().values().value()[0].x)
            .collect::<Vec<_>>(),
        vec![1.0, 2.0]
    );
    assert_eq!(
        ensemble
            .member(0)
            .unwrap()
            .atom_data()
            .occupancy(&topology, atom)
            .unwrap(),
        Some(0.5)
    );

    let independent = one_atom_topology();
    assert_eq!(
        ensemble.push(EnsembleMember::new(positions(&independent, 3.0))),
        Err(EnsembleError::TopologyMismatch)
    );
}

#[test]
fn atom_data_validates_columns_supports_mutation_and_model_topology_binding() {
    let topology = one_atom_topology();
    let atom = topology.atom_ids()[0];
    let index = topology.atom_index(atom).unwrap();
    let mut data = AtomData::new(&topology);
    assert!(data.is_empty());
    assert_eq!(data.atom_count(), 1);
    assert_eq!(data.occupancy(&topology, atom).unwrap(), None);
    data.set_occupancy(&topology, atom, Some(0.75)).unwrap();
    data.set_b_factor_at(index, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
        .unwrap();
    assert_eq!(data.occupancy_at(index).unwrap(), Some(0.75));
    assert_eq!(
        data.b_factor(&topology, atom).unwrap(),
        Some(Quantity::new(12.5, SQUARE_ANGSTROM))
    );
    data.set_b_factor(
        &topology,
        atom,
        Some(Quantity::new(0.125, NANOMETER.powi(2))),
    )
    .unwrap();
    assert!(data
        .b_factor(&topology, atom)
        .unwrap()
        .unwrap()
        .is_close(&Quantity::new(12.5, SQUARE_ANGSTROM), 1.0e-12, 1.0e-12,)
        .unwrap());
    assert!(matches!(
        data.set_b_factor(&topology, atom, Some(Quantity::new(1.0, KELVIN))),
        Err(AtomDataError::Unit(UnitError::IncompatibleUnits { .. }))
    ));
    assert!(matches!(
        data.set_b_factor(
            &topology,
            atom,
            Some(Quantity::new(f64::INFINITY, SQUARE_ANGSTROM)),
        ),
        Err(AtomDataError::NonFiniteBFactor { .. })
    ));
    assert!(matches!(
        data.set_occupancies(Vec::new()),
        Err(AtomDataError::AtomCountMismatch { .. })
    ));
    data.clear_occupancies();
    data.clear_b_factors();
    assert!(data.is_empty());
    data.set_occupancy(&topology, atom, Some(0.75)).unwrap();
    data.set_b_factor_at(index, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
        .unwrap();

    let independent = one_atom_topology();
    assert_eq!(
        data.occupancy(&independent, atom),
        Err(AtomDataError::TopologyMismatch)
    );
    let mut model = Model::new(Arc::clone(&topology), positions(&topology, 1.0)).unwrap();
    assert_eq!(
        model.set_atom_data(AtomData::new(&independent)),
        Err(ModelError::TopologyMismatch)
    );
    model.set_atom_data(data).unwrap();
    assert_eq!(model.occupancy(atom).unwrap(), Some(0.75));
    assert_eq!(
        model.b_factor(atom).unwrap(),
        Some(Quantity::new(12.5, SQUARE_ANGSTROM))
    );
}

#[test]
fn atom_data_canonical_and_custom_columns_share_dense_semantics_without_namespace_overlap() {
    let (topology, _, _) = two_bond_instances_topology();
    let mut data = AtomData::new(&topology);

    assert!(data.is_empty());
    data.set_occupancies(vec![None; 4]).unwrap();
    data.set_b_factors(Quantity::new(vec![None; 4], NANOMETER.powi(2)))
        .unwrap();
    assert!(data.is_empty());

    data.set_occupancies(vec![Some(0.25), None, None, None])
        .unwrap();
    data.set_b_factors(Quantity::new(
        vec![None, Some(0.01), None, None],
        NANOMETER.powi(2),
    ))
    .unwrap();
    data.set_property(
        "analysis_score",
        Quantity::new(vec![None, None, Some(3.0), None], DIMENSIONLESS),
    )
    .unwrap();

    assert_eq!(
        data.occupancies.as_ref().map(|column| column.unit),
        Some(DIMENSIONLESS)
    );
    assert_eq!(
        data.b_factors.as_ref().map(|column| column.unit),
        Some(SQUARE_ANGSTROM)
    );
    assert_eq!(
        data.properties().map(|(name, _)| name).collect::<Vec<_>>(),
        vec!["analysis_score"]
    );
    assert_eq!(
        data.occupancies(),
        Some([Some(0.25), None, None, None].as_slice())
    );
    let b_factors = data.b_factors().unwrap();
    assert_eq!(b_factors.unit(), SQUARE_ANGSTROM);
    assert!((b_factors.value()[1].unwrap() - 1.0).abs() < 1.0e-12);
    assert!(matches!(
        data.property("Occupancy"),
        Err(AtomDataError::ReservedPropertyName { .. })
    ));
    assert!(matches!(
        data.property("B_FACTOR"),
        Err(AtomDataError::ReservedPropertyName { .. })
    ));

    let before = data.clone();
    assert!(matches!(
        data.set_occupancies(vec![Some(f64::NAN), None, None, None]),
        Err(AtomDataError::NonFiniteOccupancy { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_occupancies(vec![Some(1.0)]),
        Err(AtomDataError::AtomCountMismatch { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_b_factors(Quantity::new(
            vec![Some(f64::INFINITY), None, None, None],
            SQUARE_ANGSTROM,
        )),
        Err(AtomDataError::NonFiniteBFactor { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_b_factors(Quantity::new(vec![Some(1.0); 4], KELVIN)),
        Err(AtomDataError::Unit(UnitError::IncompatibleUnits { .. }))
    ));
    assert_eq!(data, before);

    data.set_occupancy_at(TopologyAtomIndex::new(0), None)
        .unwrap();
    data.set_b_factor_at(TopologyAtomIndex::new(1), None)
        .unwrap();
    assert!(data.occupancies().is_none());
    assert!(data.b_factors().is_none());
    assert!(!data.is_empty());
    assert!(data.remove_property("analysis_score").unwrap());
    assert!(data.is_empty());
}

#[test]
fn atom_custom_properties_are_dense_unit_aware_and_transactional() {
    let (topology, _, _) = two_bond_instances_topology();
    let first_atom = topology.atom_ids()[0];
    let first_index = topology.atom_index(first_atom).unwrap();
    let mut data = AtomData::new(&topology);

    data.set_property(
        "partial_charge",
        Quantity::new(vec![Some(-0.4), None, Some(0.2), Some(0.2)], DIMENSIONLESS),
    )
    .unwrap();
    let charges = data.property("partial_charge").unwrap().unwrap();
    assert_eq!(charges.unit(), DIMENSIONLESS);
    assert_eq!(charges.value(), &[Some(-0.4), None, Some(0.2), Some(0.2)]);
    assert_eq!(
        data.property_value(&topology, "partial_charge", first_atom)
            .unwrap(),
        Some(Quantity::new(-0.4, DIMENSIONLESS))
    );
    assert_eq!(
        data.property_value_at("missing", first_index).unwrap(),
        None
    );

    data.set_property(
        "display_radius",
        Quantity::new(vec![Some(1.0), None, None, None], ANGSTROM),
    )
    .unwrap();
    data.set_property_value_at(
        "display_radius",
        TopologyAtomIndex::new(1),
        Some(Quantity::new(0.2, NANOMETER)),
    )
    .unwrap();
    let radii = data.property("display_radius").unwrap().unwrap();
    assert_eq!(radii.unit(), ANGSTROM);
    assert_eq!(radii.value(), &[Some(1.0), Some(2.0), None, None]);

    let before = data.clone();
    assert!(matches!(
        data.set_property(
            "partial_charge",
            Quantity::new(vec![Some(1.0)], DIMENSIONLESS)
        ),
        Err(AtomDataError::PropertyValueCountMismatch { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property("partial_charge", Quantity::new(vec![Some(1.0); 4], KELVIN)),
        Err(AtomDataError::PropertyUnit { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property(
            "partial_charge",
            Quantity::new(vec![Some(f64::NAN), None, None, None], DIMENSIONLESS)
        ),
        Err(AtomDataError::NonFinitePropertyValue { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property("", Quantity::new(vec![None; 4], DIMENSIONLESS)),
        Err(AtomDataError::InvalidPropertyName { .. })
    ));
    assert!(matches!(
        data.set_property("occupancy", Quantity::new(vec![None; 4], DIMENSIONLESS)),
        Err(AtomDataError::ReservedPropertyName { .. })
    ));
    assert!(matches!(
        data.set_property("b_factor", Quantity::new(vec![None; 4], DIMENSIONLESS)),
        Err(AtomDataError::ReservedPropertyName { .. })
    ));

    data.set_property(
        "partial_charge",
        Quantity::new(vec![None; 4], DIMENSIONLESS),
    )
    .unwrap();
    assert!(data.property("partial_charge").unwrap().is_none());
    assert!(data.remove_property("display_radius").unwrap());
    assert!(!data.remove_property("display_radius").unwrap());
    assert!(data.is_empty());
}

#[test]
fn bond_data_properties_bind_exact_topology_and_support_visualization_columns() {
    let (topology, _, _) = two_bond_instances_topology();
    let independent = two_bond_instances_topology().0;
    let first_bond = topology.bond_ids()[0];
    let first_index = topology.bond_index(first_bond).unwrap();
    let entropy_unit = KILOJOULE_PER_MOLE / KELVIN;
    let mut data = BondData::new(&topology);

    assert!(data.is_empty());
    assert_eq!(data.bond_count(), 2);
    assert!(data.is_compatible(&topology));
    assert!(!data.is_compatible(&independent));
    assert!(Arc::ptr_eq(&data.shared_topology(), &topology));
    data.set_property(
        "conformational_entropy",
        Quantity::new(vec![Some(0.012), None], entropy_unit),
    )
    .unwrap();

    let (name, column) = data.properties().next().unwrap();
    assert_eq!(name, "conformational_entropy");
    assert_eq!(column.unit(), entropy_unit);
    assert_eq!(column.value(), &[Some(0.012), None]);
    assert_eq!(
        data.property_value(&topology, "conformational_entropy", first_bond)
            .unwrap(),
        Some(Quantity::new(0.012, entropy_unit))
    );
    assert_eq!(
        data.property_value_at("conformational_entropy", TopologyBondIndex::new(1))
            .unwrap(),
        None
    );
    data.set_property_value(
        &topology,
        "display_width",
        first_bond,
        Some(Quantity::new(0.2, NANOMETER)),
    )
    .unwrap();
    data.set_property_value_at(
        "display_width",
        first_index,
        Some(Quantity::new(3.0, ANGSTROM)),
    )
    .unwrap();
    assert!(data
        .property_value_at("display_width", first_index)
        .unwrap()
        .unwrap()
        .is_close(&Quantity::new(0.3, NANOMETER), 1.0e-12, 1.0e-12)
        .unwrap());
    data.set_property_value_at("display_width", first_index, None)
        .unwrap();
    assert!(data.property("display_width").unwrap().is_none());

    assert_eq!(
        data.property_value(&independent, "conformational_entropy", first_bond),
        Err(BondDataError::TopologyMismatch)
    );
    let invalid_bond =
        InstanceBondId::new(MoleculeInstanceId::new(99), crate::core::BondId::new(0));
    assert_eq!(
        data.property_value(&topology, "conformational_entropy", invalid_bond),
        Err(BondDataError::InvalidBondId(invalid_bond))
    );
    assert_eq!(
        data.property_value_at("conformational_entropy", TopologyBondIndex::new(99)),
        Err(BondDataError::InvalidBondIndex(TopologyBondIndex::new(99)))
    );

    let before = data.clone();
    assert!(matches!(
        data.set_property(
            "conformational_entropy",
            Quantity::new(vec![Some(f64::INFINITY), None], entropy_unit)
        ),
        Err(BondDataError::NonFinitePropertyValue { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property(
            "conformational_entropy",
            Quantity::new(vec![Some(1.0), None], ANGSTROM)
        ),
        Err(BondDataError::PropertyUnit { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property(
            "conformational_entropy",
            Quantity::new(vec![Some(1.0)], entropy_unit)
        ),
        Err(BondDataError::PropertyValueCountMismatch { .. })
    ));
    assert_eq!(data, before);
    assert!(matches!(
        data.set_property("bad name", Quantity::new(vec![None; 2], entropy_unit)),
        Err(BondDataError::InvalidPropertyName { .. })
    ));

    data.set_property(
        "conformational_entropy",
        Quantity::new(vec![None; 2], entropy_unit),
    )
    .unwrap();
    assert!(data.is_empty());
    assert!(data.property("conformational_entropy").unwrap().is_none());
    assert!(!data.remove_property("conformational_entropy").unwrap());
    assert_eq!(first_index, TopologyBondIndex::new(0));
}

#[test]
fn model_and_ensemble_remap_atom_and_bond_properties_in_target_dense_order() {
    let (source, first, second) = two_bond_instances_topology();
    let edit = crate::topology::transform::retain_instances(&source, [second]).unwrap();
    let target = edit.shared_topology();
    let positions = Positions::new(
        &source,
        Quantity::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(3.0, 0.0, 0.0),
            ],
            ANGSTROM,
        ),
    )
    .unwrap();
    let mut model = Model::new(Arc::clone(&source), positions).unwrap();
    model
        .atom_data_mut()
        .set_occupancies(vec![None, None, Some(0.8), None])
        .unwrap();
    model
        .atom_data_mut()
        .set_b_factors(Quantity::new(
            vec![None, None, None, Some(4.0)],
            SQUARE_ANGSTROM,
        ))
        .unwrap();
    model
        .atom_data_mut()
        .set_property(
            "partial_charge",
            Quantity::new(
                vec![Some(-0.1), Some(0.1), Some(-0.2), Some(0.2)],
                DIMENSIONLESS,
            ),
        )
        .unwrap();
    model
        .bond_data_mut()
        .set_property(
            "conformational_entropy",
            Quantity::new(vec![Some(1.0), Some(2.0)], KILOJOULE_PER_MOLE / KELVIN),
        )
        .unwrap();

    let view = model.view();
    assert!(std::ptr::eq(model.bond_data(), view.bond_data()));
    let cloned = model.clone();
    assert_eq!(cloned.atom_data(), model.atom_data());
    assert_eq!(cloned.bond_data(), model.bond_data());

    let remapped = model.remap_to(&target, edit.mapping()).unwrap();
    assert_eq!(
        remapped.atom_data().occupancies(),
        Some([Some(0.8), None].as_slice())
    );
    assert_eq!(
        remapped.atom_data().b_factors(),
        Some(Quantity::new([None, Some(4.0)].as_slice(), SQUARE_ANGSTROM,))
    );
    assert_eq!(
        remapped
            .atom_data()
            .property("partial_charge")
            .unwrap()
            .unwrap()
            .value(),
        &[Some(-0.2), Some(0.2)]
    );
    let entropy = remapped
        .bond_data()
        .property("conformational_entropy")
        .unwrap()
        .unwrap();
    assert_eq!(entropy.unit(), KILOJOULE_PER_MOLE / KELVIN);
    assert_eq!(entropy.value(), &[Some(2.0)]);
    assert_eq!(edit.mapping().removed_instances(), &[first]);
    assert_eq!(edit.mapping().removed_bonds().len(), 1);
    assert_eq!(
        edit.mapping().bond_index_pairs().next().unwrap().0.index(),
        1
    );

    let ensemble = Ensemble::from_models(&[model]).unwrap();
    assert_eq!(
        ensemble
            .member(0)
            .unwrap()
            .bond_data()
            .property("conformational_entropy")
            .unwrap()
            .unwrap()
            .value(),
        &[Some(1.0), Some(2.0)]
    );
    assert!(std::ptr::eq(
        ensemble.member(0).unwrap().bond_data(),
        ensemble.views().next().unwrap().bond_data()
    ));
    let remapped_ensemble = ensemble.remap_to(&target, edit.mapping()).unwrap();
    assert_eq!(
        remapped_ensemble
            .member(0)
            .unwrap()
            .bond_data()
            .property("conformational_entropy")
            .unwrap()
            .unwrap()
            .value(),
        &[Some(2.0)]
    );
}

#[test]
fn model_remaps_positions_atom_data_bond_data_and_cell_together() {
    let source = one_atom_topology();
    let target = one_atom_topology();
    let mapping = TopologyMapping::between_identical_layouts(&source, &target).unwrap();
    let atom = source.atom_ids()[0];
    let mut model = Model::new(Arc::clone(&source), positions(&source, 3.0)).unwrap();
    model
        .atom_data_mut()
        .set_occupancy(&source, atom, Some(0.8))
        .unwrap();
    model
        .atom_data_mut()
        .set_b_factor(&source, atom, Some(Quantity::new(21.0, SQUARE_ANGSTROM)))
        .unwrap();
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(8.0, 9.0, 10.0), ANGSTROM),
        [true; 3],
    )
    .unwrap();
    model.set_cell(Some(cell));

    let remapped = model.remap_to(&target, &mapping).unwrap();
    let target_atom = target.atom_ids()[0];
    assert_eq!(remapped.positions().values().value()[0].x, 3.0);
    assert_eq!(remapped.occupancy(target_atom).unwrap(), Some(0.8));
    assert_eq!(
        remapped.b_factor(target_atom).unwrap(),
        Some(Quantity::new(21.0, SQUARE_ANGSTROM))
    );
    assert_eq!(remapped.cell(), Some(&cell));
    let view = remapped.view();
    assert_eq!(view.position(target_atom).unwrap().value().x, 3.0);
    assert_eq!(view.occupancy(target_atom).unwrap(), Some(0.8));
    assert_eq!(view.atom_data(), remapped.atom_data());
}

#[test]
fn ensemble_from_conformers_preserves_source_order_without_copying_conformers_to_topology() {
    let mut graph = Molecule::new();
    let atom = graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .expect("atom identifier capacity");
    let mut first = Conformer::new(ANGSTROM).unwrap();
    first
        .set_position(atom, Quantity::new(Point3::new(1.0, 0.0, 0.0), ANGSTROM))
        .unwrap();
    let first = graph.add_conformer(first).unwrap();
    let mut second = Conformer::new(ANGSTROM).unwrap();
    second
        .set_position(atom, Quantity::new(Point3::new(2.0, 0.0, 0.0), ANGSTROM))
        .unwrap();
    let second = graph.add_conformer(second).unwrap();
    let molecule = SmallMolecule::from_graph(graph);

    let ensemble = Ensemble::from_small_molecule_conformers(&molecule, [second, first]).unwrap();
    assert_eq!(
        ensemble
            .views()
            .map(|view| view.positions().values().value()[0].x)
            .collect::<Vec<_>>(),
        vec![2.0, 1.0]
    );
    assert_eq!(
        ensemble
            .topology()
            .definitions()
            .next()
            .unwrap()
            .1
            .graph()
            .conformers()
            .count(),
        0
    );
}
