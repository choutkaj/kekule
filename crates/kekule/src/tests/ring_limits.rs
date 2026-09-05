use super::*;

#[test]
fn ring_set_reports_symmetric_cycles_for_fused_rings() {
    let (mut mol, _, _) = ring_molecule(
        &["C", "C", "C", "C", "C", "C"],
        &[
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(AtomId::new(0), a, BondOrder::Single)
        .expect("bond");
    mol.add_bond(a, b, BondOrder::Single).expect("bond");
    mol.add_bond(b, AtomId::new(3), BondOrder::Single)
        .expect("bond");

    let ring_set = rings_api::perceive_ring_set(&mut mol).expect("ring perception should succeed");

    assert_eq!(ring_set.len(), 3);
    assert!(ring_set.rings().iter().all(|ring| ring.atoms.len() >= 4));
}

#[test]
fn long_chain_ring_and_smiles_traversals_are_stack_safe() {
    let mut molecule = crate::core::MoleculeEditor::new();
    let mut previous = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    for _ in 1..20_000 {
        let atom = molecule
            .add_atom(carbon())
            .expect("atom identifier capacity");
        molecule
            .add_bond(previous, atom, BondOrder::Single)
            .expect("chain bond should be valid");
        previous = atom;
    }

    let ring_set = rings_api::perceive_ring_set(molecule.working_mut())
        .expect("long chain should perceive rings");
    assert!(ring_set.is_empty());

    let written = smiles_api::write(molecule.working()).expect("long chain should write");
    assert_eq!(written.matches('C').count(), 20_000);
}

#[test]
fn ring_input_and_work_limits_return_actionable_errors() {
    let options_and_error = [
        (
            RingPerceptionOptions {
                max_atoms: 2,
                ..RingPerceptionOptions::default()
            },
            ("atoms", 3, 2),
        ),
        (
            RingPerceptionOptions {
                max_bonds: 2,
                ..RingPerceptionOptions::default()
            },
            ("bonds", 3, 2),
        ),
        (
            RingPerceptionOptions {
                max_total_work: 5,
                ..RingPerceptionOptions::default()
            },
            ("total work", 6, 5),
        ),
        (
            RingPerceptionOptions {
                max_path_expansions: 0,
                ..RingPerceptionOptions::default()
            },
            ("path expansions", 1, 0),
        ),
        (
            RingPerceptionOptions {
                max_equivalent_shortest_paths: 0,
                ..RingPerceptionOptions::default()
            },
            ("equivalent shortest paths", 1, 0),
        ),
    ];

    for (options, expected) in options_and_error {
        let (mut molecule, _, _) = ring_molecule(
            &["C", "C", "C"],
            &[BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let error = rings_api::perceive_ring_set_with_options(&mut molecule, options)
            .expect_err("configured ring resource limit should fail");
        assert_eq!(
            error,
            RingPerceptionError::ResourceLimit {
                resource: expected.0,
                observed: expected.1,
                limit: expected.2,
            }
        );
        assert!(molecule.ring_set().is_none());
    }
}

#[test]
fn symmetric_cage_returns_named_candidate_limit_error() {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = (0..4)
        .map(|_| mol.add_atom(carbon()).expect("atom identifier capacity"))
        .collect::<Vec<_>>();
    let right = (0..4)
        .map(|_| mol.add_atom(carbon()).expect("atom identifier capacity"))
        .collect::<Vec<_>>();
    for a in &left {
        for b in &right {
            mol.add_bond(*a, *b, BondOrder::Single)
                .expect("cage bond should be valid");
        }
    }
    let error = rings_api::perceive_ring_set_with_options(
        mol.working_mut(),
        RingPerceptionOptions {
            max_candidates: 2,
            ..RingPerceptionOptions::default()
        },
    )
    .expect_err("symmetric cage should hit candidate limit");
    assert!(matches!(
        error,
        RingPerceptionError::ResourceLimit {
            resource: "candidate cycles",
            observed,
            limit,
            ..
        } if observed > limit
    ));
    assert!(mol.ring_set().is_none());
}

#[test]
fn theta_graph_with_acyclic_tail_is_deterministic() {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    for _ in 0..3 {
        let middle = mol.add_atom(carbon()).expect("atom identifier capacity");
        mol.add_bond(left, middle, BondOrder::Single)
            .expect("theta edge");
        mol.add_bond(middle, right, BondOrder::Single)
            .expect("theta edge");
    }
    let tail_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let tail_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(right, tail_a, BondOrder::Single)
        .expect("tail linker");
    mol.add_bond(tail_a, tail_b, BondOrder::Single)
        .expect("tail bond");

    let first =
        rings_api::perceive_ring_set(mol.working_mut()).expect("theta graph should perceive rings");
    let first_rings = first.rings().to_vec();
    assert_eq!(first.len(), 3);
    assert!(first.rings().iter().all(|ring| ring.atoms.len() == 4));

    let second =
        rings_api::perceive_ring_set(mol.working_mut()).expect("repeat should perceive rings");
    assert_eq!(second.rings(), first_rings);
}

#[test]
fn decorated_theta_graph_does_not_overenumerate_equal_path_cycles() {
    let mut mol = crate::core::MoleculeEditor::new();
    let atoms = (0..9)
        .map(|_| mol.add_atom(carbon()).expect("atom identifier capacity"))
        .collect::<Vec<_>>();
    for (left, right) in [
        (0, 2),
        (2, 1),
        (0, 3),
        (3, 1),
        (0, 4),
        (4, 5),
        (5, 6),
        (6, 1),
        (2, 7),
        (7, 8),
        (8, 2),
    ] {
        mol.add_bond(atoms[left], atoms[right], BondOrder::Single)
            .expect("decorated theta edge");
    }

    let ring_set = rings_api::perceive_ring_set(mol.working_mut())
        .expect("decorated theta graph should perceive rings");
    let mut atom_sets = ring_set
        .rings()
        .iter()
        .map(|ring| {
            let mut atoms = ring
                .atoms
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>();
            atoms.sort_unstable();
            atoms
        })
        .collect::<Vec<_>>();
    atom_sets.sort();

    assert_eq!(
        atom_sets,
        vec![vec![0, 1, 2, 3], vec![0, 1, 2, 4, 5, 6], vec![2, 7, 8]]
    );
}

#[test]
fn cycle_size_limit_returns_structured_error() {
    let (mut mol, _, _) = ring_molecule(
        &["C", "C", "C", "C", "C", "C", "C", "C", "C", "C"],
        &[BondOrder::Single; 10],
    );
    let error = rings_api::perceive_ring_set_with_options(
        &mut mol,
        RingPerceptionOptions {
            max_cycle_size: 5,
            ..RingPerceptionOptions::default()
        },
    )
    .expect_err("large cycle should hit cycle-size limit");
    assert!(matches!(
        error,
        RingPerceptionError::ResourceLimit {
            resource: "cycle size",
            observed: 10,
            limit: 5,
            ..
        }
    ));
}

#[test]
fn focused_ring_resource_errors_are_transactional() {
    let mut molecule = crate::core::MoleculeEditor::new();
    let atoms = (0..3)
        .map(|_| {
            molecule
                .add_atom(carbon())
                .expect("atom identifier capacity")
        })
        .collect::<Vec<_>>();
    for (left, right) in [(0, 1), (1, 2), (2, 0)] {
        molecule
            .add_bond(atoms[left], atoms[right], BondOrder::Single)
            .expect("triangle bond should be valid");
    }
    let original = molecule.clone();
    let ring_options = RingPerceptionOptions {
        max_path_expansions: 0,
        ..RingPerceptionOptions::default()
    };
    let error = rings_api::perceive_ring_set_with_options(molecule.working_mut(), ring_options)
        .expect_err("ring limit should fail focused perception");
    assert!(matches!(error, RingPerceptionError::ResourceLimit { .. }));
    assert_eq!(molecule, original);

    let mut aromatic = read_smiles("c1ccccc1").expect("benzene should parse");
    aromatic
        .canonicalize_fixture()
        .expect("benzene should normalize");
    let error = aromaticity_api::perceive_aromaticity_with_ring_options(
        &mut aromatic,
        AromaticityModel::RdkitLike,
        RingPerceptionOptions {
            max_atoms: 0,
            ..RingPerceptionOptions::default()
        },
    )
    .expect_err("aromaticity should propagate ring limit");
    assert!(matches!(error, AromaticityError::RingPerception(_)));
}

#[test]
fn standalone_aromaticity_reuses_an_installed_ring_set() {
    let mut molecule = read_smiles("C1=CC=CC=C1").expect("benzene should parse");
    valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
        .expect("valence perception");
    let installed =
        rings_api::perceive_ring_set(&mut molecule).expect("ring perception should succeed");

    aromaticity_api::perceive_aromaticity_with_ring_options(
        &mut molecule,
        AromaticityModel::RdkitLike,
        RingPerceptionOptions {
            max_atoms: 0,
            ..RingPerceptionOptions::default()
        },
    )
    .expect("installed rings should bypass the impossible ring limit");

    assert_eq!(molecule.ring_set(), Some(&installed));
    assert!(molecule.perception().has_aromaticity());
}

#[test]
fn standalone_aromaticity_computes_a_missing_ring_set() {
    let mut baseline = read_smiles("C1=CC=CC=C1").expect("benzene should parse");
    valence_api::perceive_valence(&mut baseline, ValenceModel::RdkitLike)
        .expect("valence perception");

    for membership_only in [false, true] {
        let mut molecule = baseline.clone();
        if membership_only {
            rings_api::perceive_ring_membership(&mut molecule);
            assert!(molecule.ring_membership().is_some());
            assert!(molecule.ring_set().is_none());
        }

        aromaticity_api::perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike)
            .expect("aromaticity should compute a missing ring basis");

        assert_eq!(molecule.ring_set().expect("computed ring basis").len(), 1);
        assert!(molecule.perception().has_aromaticity());
    }
}
