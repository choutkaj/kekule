use std::sync::Arc;

use kekule::{
    core::{
        Atom, AtomId, BondId, BondOrder, Element, MoleculeEditor, StereoDescriptor,
        StereoElementKind,
    },
    geometry::{PeriodicCell, Point3, Vector3},
    smiles,
    stereo::{assign_cip_descriptors, CipAssignmentOptions, CipRankingError},
    structure::{
        measure::{self, BondDihedral, BondDihedralError, MeasurementError},
        Ensemble, EnsembleMember, Model, Positions,
    },
    topology::{InstanceAtomId, InstanceBondId, MoleculeInstanceId},
    units::{Quantity, DEGREE, NANOMETER},
};

fn make_model(smiles: &str, perceive: bool) -> Model {
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
    Model::from_molecule(
        &molecule,
        &Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

fn bond(model: &Model, b: usize, c: usize) -> InstanceBondId {
    let ids = model.atom_ids();
    let instance = model.topology().molecule(ids[b].molecule()).unwrap();
    instance.qualify_bond(
        instance
            .molecule()
            .bond_between(ids[b].atom(), ids[c].atom())
            .unwrap()
            .unwrap(),
    )
}

fn definition(model: &Model, b: usize, c: usize) -> BondDihedral {
    BondDihedral::new(&model.shared_topology(), bond(model, b, c)).unwrap()
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
        let model = make_model(smiles, true);
        let topology = model.shared_topology();
        let before = topology
            .molecule(model.atom_ids()[b].molecule())
            .unwrap()
            .molecule()
            .perception()
            .clone();
        let selected = definition(&model, b, c);
        let ids = model.atom_ids();
        assert_eq!(
            selected.atoms(),
            Some([ids[a], ids[b], ids[c], ids[d]]),
            "{smiles}"
        );
        assert_eq!(
            selected.measure(model.view()).unwrap(),
            Some(measure::dihedral(model.view(), ids[a], ids[b], ids[c], ids[d]).unwrap())
        );
        assert_eq!(
            topology
                .molecule(ids[b].molecule())
                .unwrap()
                .molecule()
                .perception(),
            &before
        );
        assert!(Arc::ptr_eq(&topology, &model.shared_topology()));
    }
}

#[test]
fn tied_references_are_smallest_ids_and_do_not_follow_coordinates() {
    let mut model = make_model("CC(C)C(C)C", true);
    let selected = definition(&model, 1, 3);
    let ids = model.atom_ids().to_vec();
    assert_eq!(selected.atoms(), Some([ids[0], ids[1], ids[3], ids[4]]));
    // The alternative methyl reference stays off axis, while the chosen one
    // becomes collinear. It must not be substituted just to obtain an angle.
    for (id, point) in [
        (ids[0], Point3::new(0.0, -1.0, 0.0)),
        (ids[1], Point3::origin()),
        (ids[3], Point3::new(0.0, 1.0, 0.0)),
        (ids[4], Point3::new(0.0, 1.0, 1.0)),
    ] {
        model
            .set_position(id, Quantity::new(point, NANOMETER))
            .unwrap();
    }
    assert!(measure::dihedral(model.view(), ids[2], ids[1], ids[3], ids[4]).is_ok());
    assert_eq!(selected.measure(model.view()).unwrap(), None);
    assert_eq!(definition(&model, 1, 3).atoms(), selected.atoms());
    assert_eq!(
        measure::bond_dihedral(model.view(), selected.bond()).unwrap(),
        None
    );
    model
        .set_position(ids[0], Quantity::new(Point3::new(1.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    let angle = selected
        .measure(model.view())
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
    let model = Model::from_molecule(&molecule, &Positions::zeros(6)).unwrap();
    let ids = model.atom_ids();
    assert_eq!(
        definition(&model, 1, 3).atoms(),
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
        let model = make_model(smiles, true);
        let iter = measure::bond_dihedrals(model.view());
        assert_eq!(iter.len(), model.bond_count());
        let values = iter.collect::<Vec<_>>();
        assert_eq!(
            values.iter().map(|(id, _)| *id).collect::<Vec<_>>(),
            model.bond_ids()
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
    let model = make_model("[He]", true);
    assert_eq!(measure::bond_dihedrals(model.view()).len(), 0);
}

#[test]
fn geometry_absence_is_distinct_from_overflow_and_foreign_topologies() {
    let mut model = make_model("CCCC", true);
    let selected = definition(&model, 1, 2);
    let ids = model.atom_ids().to_vec();
    model
        .set_position(
            ids[0],
            Quantity::new(Point3::new(f64::MAX, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    model
        .set_position(
            ids[1],
            Quantity::new(Point3::new(-f64::MAX, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    assert_eq!(
        selected.measure(model.view()),
        Err(BondDihedralError::Measurement(
            MeasurementError::NumericalFailure
        ))
    );
    for &id in &ids {
        model
            .set_position(id, Quantity::new(Point3::origin(), NANOMETER))
            .unwrap();
    }
    assert_eq!(selected.measure(model.view()).unwrap(), None);
    let foreign = make_model("CCCC", true);
    assert_eq!(
        selected.measure(foreign.view()),
        Err(BondDihedralError::TopologyMismatch)
    );
    let absent = definition(&model, 0, 1);
    assert_eq!(absent.atoms(), None);
    assert_eq!(
        absent.measure(foreign.view()),
        Err(BondDihedralError::TopologyMismatch)
    );
    for invalid in [
        InstanceBondId::new(ids[0].molecule(), BondId::new(99)),
        InstanceBondId::new(MoleculeInstanceId::new(99), BondId::new(0)),
    ] {
        assert!(
            matches!(BondDihedral::new(&model.shared_topology(), invalid), Err(BondDihedralError::InvalidBondId(id)) if id == invalid)
        );
    }
}

#[test]
fn incomplete_ranking_is_an_error_not_an_atom_id_tie() {
    let model = make_model("CC(CO)CC", true);
    let axis = bond(&model, 1, 4);
    for (options, expected) in [
        (
            CipAssignmentOptions {
                max_depth: 0,
                max_nodes: 100_000,
            },
            CipRankingError::DepthLimitExceeded { max_depth: 0 },
        ),
        (
            CipAssignmentOptions {
                max_depth: 64,
                max_nodes: 1,
            },
            CipRankingError::ResourceLimitExceeded { max_nodes: 1 },
        ),
    ] {
        assert!(
            matches!(BondDihedral::with_options(&model.shared_topology(), axis, options), Err(BondDihedralError::Ranking { bond, error }) if bond == axis && error == expected)
        );
    }
    let unperceived = make_model("CC(CO)CC", false);
    assert!(matches!(
        BondDihedral::new(&unperceived.shared_topology(), bond(&unperceived, 1, 4)),
        Err(BondDihedralError::Ranking {
            error: CipRankingError::UnknownHydrogenCount { .. },
            ..
        })
    ));
    // A proved element ordering needs no complete digraph or deep exploration.
    let model = make_model("FC(Cl)C(Br)I", true);
    let tiny = BondDihedral::with_options(
        &model.shared_topology(),
        bond(&model, 1, 3),
        CipAssignmentOptions {
            max_depth: 0,
            max_nodes: 1,
        },
    )
    .unwrap();
    assert_eq!(tiny.atoms(), definition(&model, 1, 3).atoms());
    let unperceived_elements = make_model("FC(Cl)C(Br)I", false);
    assert_eq!(
        definition(&unperceived_elements, 1, 3).atoms(),
        tiny.atoms()
    );
    let tied = make_model("CC(C)C(C)C", true);
    assert!(matches!(
        BondDihedral::with_options(
            &tied.shared_topology(),
            bond(&tied, 1, 3),
            CipAssignmentOptions {
                max_depth: 0,
                max_nodes: 100_000
            }
        ),
        Err(BondDihedralError::Ranking {
            error: CipRankingError::DepthLimitExceeded { .. },
            ..
        })
    ));
    assert!(measure::bond_dihedrals(unperceived.view())
        .any(|(_, result)| matches!(result, Err(BondDihedralError::Ranking { .. }))));
}

#[test]
fn references_preserve_instance_identity_and_work_on_ensemble_views() {
    let source = make_model("CCCC", true);
    let molecule = source.topology().molecules().next().unwrap().molecule();
    let mut builder = Model::builder();
    let definition_id = builder.add_molecule_definition(molecule).unwrap();
    let first = builder
        .add_instance(definition_id, source.positions())
        .unwrap();
    let second = builder
        .add_instance(definition_id, source.positions())
        .unwrap();
    let model = builder.build().unwrap();
    let axis = InstanceBondId::new(second, BondId::new(1));
    let selected = BondDihedral::new(&model.shared_topology(), axis).unwrap();
    assert_eq!(
        selected.atoms().unwrap(),
        [0, 1, 2, 3].map(|i| InstanceAtomId::new(second, AtomId::new(i)))
    );
    assert_ne!(first, second);
    let expected = selected.measure(model.view()).unwrap();
    let mut ensemble = Ensemble::new(model.shared_topology());
    let mut member = EnsembleMember::new(model.positions().clone());
    member.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    ensemble.push(member).unwrap();
    assert_eq!(
        selected
            .measure(ensemble.member(0).unwrap().as_model())
            .unwrap(),
        expected
    );
}

#[test]
fn represented_remote_stereo_distinguishes_constitutionally_tied_ligands() {
    // Changing only the configuration of both terminal centers swaps their
    // relative CIP preference without changing any atom identifiers.
    let left = make_model("C[C@H](F)C(CC)[C@H](F)C", true);
    let right = make_model("C[C@@H](F)C(CC)[C@@H](F)C", true);
    let left_atoms = definition(&left, 3, 4).atoms().unwrap();
    let right_atoms = definition(&right, 3, 4).atoms().unwrap();
    assert!(matches!(left_atoms[0].atom().raw(), 1 | 6));
    assert!(matches!(right_atoms[0].atom().raw(), 1 | 6));
    assert_ne!(left_atoms[0].atom(), right_atoms[0].atom());
    // For this simple enantiomorphic pair Rule 5 prefers R to S. Check the
    // direction as well as sensitivity to stereo, without installing labels
    // on either source model.
    for (model, atoms) in [(&left, left_atoms), (&right, right_atoms)] {
        let mut labelled = model
            .topology()
            .molecules()
            .next()
            .unwrap()
            .molecule()
            .clone();
        assign_cip_descriptors(&mut labelled).unwrap();
        let (id, _) = labelled.stereo_elements().find(|(_, element)| matches!(&element.kind, StereoElementKind::Tetrahedral(stereo) if stereo.center == atoms[0].atom())).unwrap();
        assert_eq!(
            labelled.cip_descriptor(id).unwrap(),
            Some(StereoDescriptor::R)
        );
    }
}

#[test]
fn embedded_alkene_stereo_uses_sequence_priority_before_atom_ids() {
    for (smiles, expected_reference) in [(r"C/C=C/C(CC)/C=C\C", 6), (r"C/C=C\C(CC)/C=C/C", 2)] {
        let model = make_model(smiles, true);
        let atoms = definition(&model, 3, 4).atoms().unwrap();
        assert_eq!(atoms[0], model.atom_ids()[expected_reference], "{smiles}");
    }
}
