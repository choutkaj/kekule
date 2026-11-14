use kekule::diagnostics::{LimitExceeded, Resource};
use std::sync::Arc;

use kekule::{
    math::{PeriodicCell, Point3, Vector3},
    molecule::{Atom, AtomId, BondId, BondOrder, Element, MoleculeEditor},
    smiles,
    stereo::{
        cip::{assign_cip_descriptors, CipAssignmentOptions, CipRankingError},
        StereoDescriptor, StereoElementKind,
    },
    structure::{
        measure::{self, BondDihedral, BondDihedralError, ConnectivityCheck},
        CoordinateError, Geometry, Positions, Structure, StructureSet,
    },
    topology::{InstanceAtomId, InstanceBondId, MoleculeInstanceId},
    units::{Quantity, DEGREE, NANOMETER},
};

fn make_structure(smiles: &str, perceive: bool) -> Structure {
    let mut molecule = smiles::to_molecules(smiles).unwrap().remove(0);
    if perceive {
        molecule.perceive().unwrap();
    }
    let points = molecule
        .atom_ids()
        .enumerate()
        .map(|(i, _)| {
            let t = i as f64;
            Point3::new(t, t * t, t * t * t)
        })
        .collect::<Vec<_>>();
    Structure::from_molecule(
        molecule.clone(),
        &Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

fn bond(structure: &Structure, b: usize, c: usize) -> InstanceBondId {
    let ids = structure.topology().atom_ids();
    let instance = structure.topology().instance(ids[b].instance()).unwrap();
    instance.qualify_bond(
        instance
            .molecule()
            .bond_between(ids[b].atom(), ids[c].atom())
            .unwrap(),
    )
}

fn definition(structure: &Structure, b: usize, c: usize) -> BondDihedral {
    BondDihedral::new(
        structure.shared_topology(),
        bond(structure, b, c),
        Default::default(),
    )
    .unwrap()
}

#[test]
fn references_follow_cip_elements_isotopes_deep_ligands_and_rings() {
    for (smiles, b, c, a, d) in [
        ("FC(Cl)C(Br)I", 1, 3, 2, 5),
        ("[12CH3]C([13CH3])CC", 1, 3, 2, 4),
        ("CC(CO)CC", 1, 4, 2, 5),
        ("CC(C=O)CC", 1, 4, 2, 5),
        ("CC(C1CC1)CC", 1, 5, 2, 6),
        ("CC(c1ccccc1)CC", 1, 8, 2, 9),
        ("[H]C([2H])CC", 1, 3, 2, 4),
    ] {
        let structure = make_structure(smiles, true);
        let topology = structure.shared_topology();
        let before = topology
            .instance(structure.topology().atom_ids()[b].instance())
            .unwrap()
            .molecule()
            .perception()
            .clone();
        let selected = definition(&structure, b, c);
        let ids = structure.topology().atom_ids();
        assert_eq!(
            selected.atoms(),
            Some([ids[a], ids[b], ids[c], ids[d]]),
            "{smiles}"
        );
        assert_eq!(
            selected.measure(structure.as_structure_view()).unwrap(),
            Some(
                measure::dihedral(
                    structure.as_structure_view(),
                    ids[a],
                    ids[b],
                    ids[c],
                    ids[d],
                    ConnectivityCheck::Unrestricted
                )
                .unwrap()
            )
        );
        assert_eq!(
            topology
                .instance(ids[b].instance())
                .unwrap()
                .molecule()
                .perception(),
            &before
        );
        assert!(Arc::ptr_eq(&topology, &structure.shared_topology()));
    }
}

#[test]
fn tied_references_are_smallest_ids_and_do_not_follow_coordinates() {
    let mut structure = make_structure("CC(C)C(C)C", true);
    let selected = definition(&structure, 1, 3);
    let ids = structure.topology().atom_ids().to_vec();
    assert_eq!(selected.atoms(), Some([ids[0], ids[1], ids[3], ids[4]]));
    // The alternative methyl reference stays off axis, while the chosen one
    // becomes collinear. It must not be substituted just to obtain an angle.
    for (id, point) in [
        (ids[0], Point3::new(0.0, -1.0, 0.0)),
        (ids[1], Point3::origin()),
        (ids[3], Point3::new(0.0, 1.0, 0.0)),
        (ids[4], Point3::new(0.0, 1.0, 1.0)),
    ] {
        structure
            .set_position(id, Quantity::new(point, NANOMETER))
            .unwrap();
    }
    assert!(measure::dihedral(
        structure.as_structure_view(),
        ids[2],
        ids[1],
        ids[3],
        ids[4],
        ConnectivityCheck::Unrestricted
    )
    .is_ok());
    assert_eq!(
        selected.measure(structure.as_structure_view()).unwrap(),
        None
    );
    assert_eq!(definition(&structure, 1, 3).atoms(), selected.atoms());
    assert_eq!(
        measure::bond_dihedral(structure.as_structure_view(), selected.bond()).unwrap(),
        None
    );
    structure
        .set_position(ids[0], Quantity::new(Point3::new(1.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    let angle = selected
        .measure(structure.as_structure_view())
        .unwrap()
        .unwrap()
        .value_in(DEGREE)
        .unwrap();
    assert!((angle + 90.0).abs() < 1e-10);
}

#[test]
fn id_ties_ignore_bond_insertion_order_and_endpoint_orientation() {
    let mut editor = MoleculeEditor::new();
    let atoms = (0..6)
        .map(|_| {
            editor
                .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
                .unwrap()
        })
        .collect::<Vec<_>>();
    // Higher-ID alternatives are inserted first, and the axis is stored C-B.
    for (a, b) in [(2, 1), (5, 3), (3, 1), (4, 3), (0, 1)] {
        editor
            .add_bond(atoms[a], atoms[b], BondOrder::Single)
            .unwrap();
    }
    let mut molecule = editor.finish().unwrap();
    molecule.perceive().unwrap();
    let structure = Structure::from_molecule(molecule.clone(), &Positions::zeros(6)).unwrap();
    let ids = structure.topology().atom_ids();
    assert_eq!(
        definition(&structure, 1, 3).atoms(),
        Some([ids[0], ids[1], ids[3], ids[4]])
    );
}

#[test]
fn iterator_retains_all_bonds_and_none_is_not_a_rotatability_filter() {
    for (smiles, expected) in [
        ("CCCC", vec![false, true, false]),
        ("CC=CC", vec![false, true, false]),
        ("CC", vec![false]),
        ("C1CC1", vec![false; 3]),
    ] {
        let structure = make_structure(smiles, true);
        let iter = measure::bond_dihedrals(structure.as_structure_view());
        assert_eq!(iter.len(), structure.topology().bond_count());
        let values = iter.collect::<Vec<_>>();
        assert_eq!(
            values.iter().map(|(id, _)| *id).collect::<Vec<_>>(),
            structure.topology().bond_ids()
        );
        assert_eq!(
            values
                .into_iter()
                .map(|(_, result)| result.unwrap().is_some())
                .collect::<Vec<_>>(),
            expected,
            "{smiles}"
        );
    }
    let structure = make_structure("[He]", true);
    assert_eq!(
        measure::bond_dihedrals(structure.as_structure_view()).len(),
        0
    );
}

#[test]
fn geometry_absence_is_distinct_from_overflow_and_foreign_topologies() {
    let mut structure = make_structure("CCCC", true);
    let selected = definition(&structure, 1, 2);
    let ids = structure.topology().atom_ids().to_vec();
    structure
        .set_position(
            ids[0],
            Quantity::new(Point3::new(f64::MAX, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    structure
        .set_position(
            ids[1],
            Quantity::new(Point3::new(-f64::MAX, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    assert_eq!(
        selected.measure(structure.as_structure_view()),
        Err(BondDihedralError::Measurement(
            CoordinateError::NumericalFailure
        ))
    );
    for &id in &ids {
        structure
            .set_position(id, Quantity::new(Point3::origin(), NANOMETER))
            .unwrap();
    }
    assert_eq!(
        selected.measure(structure.as_structure_view()).unwrap(),
        None
    );
    let foreign = make_structure("CCCC", true);
    assert_eq!(
        selected.measure(foreign.as_structure_view()),
        Err(BondDihedralError::TopologyMismatch)
    );
    let absent = definition(&structure, 0, 1);
    assert_eq!(absent.atoms(), None);
    assert_eq!(
        absent.measure(foreign.as_structure_view()),
        Err(BondDihedralError::TopologyMismatch)
    );
    for invalid in [
        InstanceBondId::new(ids[0].instance(), BondId::new(99)),
        InstanceBondId::new(MoleculeInstanceId::new(99), BondId::new(0)),
    ] {
        assert!(
            matches!(BondDihedral::new(structure.shared_topology(), invalid, Default::default()), Err(BondDihedralError::InvalidBondId(id)) if id == invalid)
        );
    }
}

#[test]
fn losing_reference_ties_do_not_exhaust_ranking_bounds() {
    // The ethyl reference wins after one shell, regardless of its position
    // among the two tied methyl references.
    for (smiles, b, a) in [
        ("CC(C)(CC)CO", 1, 3),
        ("CC(CC)(C)CO", 1, 2),
        ("CCC(C)(C)CO", 2, 1),
    ] {
        let structure = make_structure(smiles, true);
        let ids = structure.topology().atom_ids();
        for options in [
            CipAssignmentOptions::default().with_max_depth(1),
            CipAssignmentOptions::default().with_max_nodes(4),
        ] {
            let selected =
                BondDihedral::new(structure.shared_topology(), bond(&structure, b, 5), options)
                    .unwrap();
            assert_eq!(selected.atoms(), Some([ids[a], ids[b], ids[5], ids[6]]));
        }
    }
    // Losing branches may also extend beyond the default depth limit.
    let chain = "C".repeat(70);
    let structure = make_structure(&format!("C({chain})({chain})(C(F)F)CO"), true);
    let ids = structure.topology().atom_ids();
    assert_eq!(
        definition(&structure, 0, 144).atoms(),
        Some([ids[141], ids[0], ids[144], ids[145]])
    );
}

#[test]
fn only_possible_maxima_consume_further_ranking_budget() {
    // The tert-butyl reference loses to both fluorinated references at shell 1.
    // Its next shell would exceed seven nodes, while the tied maxima complete
    // within that budget and must use the smaller atom ID.
    let structure = make_structure("C(C(C)(C)C)(C(F)C)(C(F)C)CO", true);
    let axis = bond(&structure, 0, 11);
    let ids = structure.topology().atom_ids();
    let selected = BondDihedral::new(
        structure.shared_topology(),
        axis,
        CipAssignmentOptions::default()
            .with_max_depth(2)
            .with_max_nodes(7),
    )
    .unwrap();
    assert_eq!(selected.atoms(), Some([ids[5], ids[0], ids[11], ids[12]]));
    // Discarding a loser must not turn an unfinished tie among winners into
    // an atom-ID tie, for either kind of resource bound.
    for (options, expected) in [
        (
            CipAssignmentOptions::default()
                .with_max_depth(1)
                .with_max_nodes(7),
            CipRankingError::LimitExceeded(LimitExceeded::first_beyond(Resource::DigraphDepth, 1)),
        ),
        (
            CipAssignmentOptions::default()
                .with_max_depth(2)
                .with_max_nodes(4),
            CipRankingError::LimitExceeded(LimitExceeded::first_beyond(Resource::DigraphNodes, 4)),
        ),
    ] {
        assert!(matches!(
            BondDihedral::new(structure.shared_topology(), axis, options),
            Err(BondDihedralError::Ranking { error, .. }) if error == expected
        ));
    }
}

#[test]
fn incomplete_ranking_is_an_error_not_an_atom_id_tie() {
    let structure = make_structure("CC(CO)CC", true);
    let axis = bond(&structure, 1, 4);
    for (options, expected) in [
        (
            CipAssignmentOptions::default()
                .with_max_depth(0)
                .with_max_nodes(100_000),
            CipRankingError::LimitExceeded(LimitExceeded::first_beyond(Resource::DigraphDepth, 0)),
        ),
        (
            CipAssignmentOptions::default()
                .with_max_depth(64)
                .with_max_nodes(1),
            CipRankingError::LimitExceeded(LimitExceeded::first_beyond(Resource::DigraphNodes, 1)),
        ),
    ] {
        assert!(
            matches!(BondDihedral::new(structure.shared_topology(), axis, options), Err(BondDihedralError::Ranking { bond, error }) if bond == axis && error == expected)
        );
    }
    let unperceived = make_structure("CC(CO)CC", false);
    assert!(matches!(
        BondDihedral::new(
            unperceived.shared_topology(),
            bond(&unperceived, 1, 4),
            Default::default()
        ),
        Err(BondDihedralError::Ranking {
            error: CipRankingError::UnknownHydrogenCount { .. },
            ..
        })
    ));
    // A proved element ordering needs no complete digraph or deep exploration.
    let structure = make_structure("FC(Cl)C(Br)I", true);
    let tiny = BondDihedral::new(
        structure.shared_topology(),
        bond(&structure, 1, 3),
        CipAssignmentOptions::default()
            .with_max_depth(0)
            .with_max_nodes(1),
    )
    .unwrap();
    assert_eq!(tiny.atoms(), definition(&structure, 1, 3).atoms());
    let unperceived_elements = make_structure("FC(Cl)C(Br)I", false);
    assert_eq!(
        definition(&unperceived_elements, 1, 3).atoms(),
        tiny.atoms()
    );
    let tied = make_structure("CC(C)C(C)C", true);
    assert!(matches!(
            BondDihedral::new(
                tied.shared_topology(),
                bond(&tied, 1, 3),
                CipAssignmentOptions::default()
    .with_max_depth(0)
    .with_max_nodes(100_000)
            ),
            Err(BondDihedralError::Ranking {
                error: CipRankingError::LimitExceeded(limit),
                ..
            }) if limit.resource() == Resource::DigraphDepth
        ));
    assert!(measure::bond_dihedrals(unperceived.as_structure_view())
        .any(|(_, result)| matches!(result, Err(BondDihedralError::Ranking { .. }))));
}

#[test]
fn references_preserve_instance_identity_and_work_on_ensemble_views() {
    let source = make_structure("CCCC", true);
    let molecule = source.topology().instances().next().unwrap().molecule();
    let mut builder = Structure::builder();
    let definition_id = builder.add_molecule_definition(molecule.clone()).unwrap();
    let first = builder
        .add_instance(definition_id, source.positions())
        .unwrap();
    let second = builder
        .add_instance(definition_id, source.positions())
        .unwrap();
    let structure = builder.build().unwrap();
    let axis = InstanceBondId::new(second, BondId::new(1));
    let selected =
        BondDihedral::new(structure.shared_topology(), axis, Default::default()).unwrap();
    assert_eq!(
        selected.atoms().unwrap(),
        [0, 1, 2, 3].map(|i| InstanceAtomId::new(second, AtomId::new(i)))
    );
    assert_ne!(first, second);
    let expected = selected.measure(structure.as_structure_view()).unwrap();
    let member = Geometry::new(structure.positions().clone()).with_cell(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    );
    let ensemble = StructureSet::from_geometries(structure.shared_topology(), [member]).unwrap();
    assert_eq!(
        selected
            .measure(ensemble.get(0).unwrap().as_structure_view())
            .unwrap(),
        expected
    );
}

#[test]
fn represented_remote_stereo_distinguishes_constitutionally_tied_ligands() {
    // Changing only the configuration of both terminal centers swaps their
    // relative CIP preference without changing any atom identifiers.
    for (left_smiles, right_smiles, c, other_reference) in [
        ("C[C@H](F)C(CC)[C@H](F)C", "C[C@@H](F)C(CC)[C@@H](F)C", 4, 6),
        // A third, lower-priority methyl reference must stay discarded when
        // the tied maxima are reranked with auxiliary stereo descriptors.
        (
            "C[C@H](F)C(C)(CC)[C@H](F)C",
            "C[C@@H](F)C(C)(CC)[C@@H](F)C",
            5,
            7,
        ),
    ] {
        let left = make_structure(left_smiles, true);
        let right = make_structure(right_smiles, true);
        let left_atoms = definition(&left, 3, c).atoms().unwrap();
        let right_atoms = definition(&right, 3, c).atoms().unwrap();
        assert!([1, other_reference].contains(&left_atoms[0].atom().raw()));
        assert!([1, other_reference].contains(&right_atoms[0].atom().raw()));
        assert_ne!(left_atoms[0].atom(), right_atoms[0].atom());
        // For this simple enantiomorphic pair Rule 5 prefers R to S. Check the
        // direction as well as sensitivity to stereo, without installing labels
        // on either source structure.
        for (structure, atoms) in [(&left, left_atoms), (&right, right_atoms)] {
            let mut labelled = structure
                .topology()
                .instances()
                .next()
                .unwrap()
                .molecule()
                .clone();
            assign_cip_descriptors(&mut labelled, Default::default()).unwrap();
            let (id, _) = labelled.stereo_elements().find(|(_, element)| matches!(&element.kind, StereoElementKind::Tetrahedral(stereo) if stereo.center == atoms[0].atom())).unwrap();
            assert_eq!(
                labelled
                    .perception()
                    .cip()
                    .and_then(|cip| cip.cip_descriptor(id)),
                Some(StereoDescriptor::R)
            );
        }
    }
}

#[test]
fn embedded_alkene_stereo_uses_sequence_priority_before_atom_ids() {
    for (smiles, expected_reference) in [(r"C/C=C/C(CC)/C=C\C", 6), (r"C/C=C\C(CC)/C=C/C", 2)] {
        let structure = make_structure(smiles, true);
        let atoms = definition(&structure, 3, 4).atoms().unwrap();
        assert_eq!(
            atoms[0],
            structure.topology().atom_ids()[expected_reference],
            "{smiles}"
        );
    }
}
