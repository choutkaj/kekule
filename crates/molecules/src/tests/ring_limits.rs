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
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    mol.add_bond(AtomId::new(0), a, BondOrder::Single)
        .expect("bond");
    mol.add_bond(a, b, BondOrder::Single).expect("bond");
    mol.add_bond(b, AtomId::new(3), BondOrder::Single)
        .expect("bond");

    let ring_set =
        perception_api::perceive_ring_set(&mut mol).expect("ring perception should succeed");

    assert_eq!(ring_set.len(), 3);
    assert!(ring_set.rings().iter().all(|ring| ring.atoms.len() >= 4));
}

#[test]
fn long_chain_ring_and_smiles_traversals_are_stack_safe() {
    let mut molecule = SmallMolecule::default();
    let mut previous = molecule.graph_mut().add_atom(carbon());
    for _ in 1..20_000 {
        let atom = molecule.graph_mut().add_atom(carbon());
        molecule
            .graph_mut()
            .add_bond(previous, atom, BondOrder::Single)
            .expect("chain bond should be valid");
        previous = atom;
    }

    let ring_set = perception_api::perceive_ring_set(molecule.graph_mut())
        .expect("long chain should perceive rings");
    assert!(ring_set.is_empty());
    assert_eq!(ring_set.work().atom_count, 20_000);
    assert!(ring_set.work().stack_peak >= 20_000);

    let written = smiles_api::write_with_options(&molecule, SmilesWriteOptions)
        .expect("long chain should write");
    assert_eq!(written.matches('C').count(), 20_000);
}

#[test]
fn ladder_ring_work_is_instrumented() {
    let mut mol = Molecule::new();
    let top = (0..12).map(|_| mol.add_atom(carbon())).collect::<Vec<_>>();
    let bottom = (0..12).map(|_| mol.add_atom(carbon())).collect::<Vec<_>>();
    for index in 0..11 {
        mol.add_bond(top[index], top[index + 1], BondOrder::Single)
            .expect("top rail");
        mol.add_bond(bottom[index], bottom[index + 1], BondOrder::Single)
            .expect("bottom rail");
    }
    for index in 0..12 {
        mol.add_bond(top[index], bottom[index], BondOrder::Single)
            .expect("rung");
    }

    let ring_set =
        perception_api::perceive_ring_set(&mut mol).expect("ladder should perceive rings");
    let work = ring_set.work();
    assert_eq!(ring_set.len(), 11);
    assert!(work.candidate_cycles >= ring_set.len());
    assert!(work.equivalent_shortest_paths >= work.candidate_cycles);
    assert!(work.path_expansions > 0);
    assert!(work.queue_peak > 0);
    assert!(work.stack_peak > 0);
    assert!(work.total_work >= work.atom_count + work.bond_count);
}

#[test]
fn symmetric_cage_returns_named_candidate_limit_error() {
    let mut mol = Molecule::new();
    let left = (0..4).map(|_| mol.add_atom(carbon())).collect::<Vec<_>>();
    let right = (0..4).map(|_| mol.add_atom(carbon())).collect::<Vec<_>>();
    for a in &left {
        for b in &right {
            mol.add_bond(*a, *b, BondOrder::Single)
                .expect("cage bond should be valid");
        }
    }
    let error = perception_api::perceive_ring_set_with_options(
        &mut mol,
        RingPerceptionOptions {
            max_candidates: 2,
            ..RingPerceptionOptions::default()
        },
    )
    .expect_err("symmetric cage should hit candidate limit");
    assert_eq!(error.resource, "candidate cycles");
    assert!(error.observed > error.limit);
    assert!(mol.ring_set().is_none());
}

#[test]
fn theta_graph_and_disconnected_mixture_are_deterministic() {
    let mut mol = Molecule::new();
    let left = mol.add_atom(carbon());
    let right = mol.add_atom(carbon());
    for _ in 0..3 {
        let middle = mol.add_atom(carbon());
        mol.add_bond(left, middle, BondOrder::Single)
            .expect("theta edge");
        mol.add_bond(middle, right, BondOrder::Single)
            .expect("theta edge");
    }
    let tail_a = mol.add_atom(carbon());
    let tail_b = mol.add_atom(carbon());
    mol.add_bond(tail_a, tail_b, BondOrder::Single)
        .expect("disconnected tail");

    let first =
        perception_api::perceive_ring_set(&mut mol).expect("theta graph should perceive rings");
    let first_rings = first.rings().to_vec();
    assert_eq!(first.len(), 3);
    assert!(first.rings().iter().all(|ring| ring.atoms.len() == 4));

    let second = perception_api::perceive_ring_set(&mut mol).expect("repeat should perceive rings");
    assert_eq!(second.rings(), first_rings);
}

#[test]
fn cycle_size_limit_returns_structured_error() {
    let (mut mol, _, _) = ring_molecule(
        &["C", "C", "C", "C", "C", "C", "C", "C", "C", "C"],
        &[BondOrder::Single; 10],
    );
    let error = perception_api::perceive_ring_set_with_options(
        &mut mol,
        RingPerceptionOptions {
            max_cycle_size: 5,
            ..RingPerceptionOptions::default()
        },
    )
    .expect_err("large cycle should hit cycle-size limit");
    assert_eq!(error.resource, "cycle size");
    assert_eq!(error.observed, 10);
    assert_eq!(error.limit, 5);
}

#[test]
fn ring_resource_errors_propagate_transactionally() {
    let mut molecule = SmallMolecule::default();
    let atoms = (0..3)
        .map(|_| molecule.graph_mut().add_atom(carbon()))
        .collect::<Vec<_>>();
    for (left, right) in [(0, 1), (1, 2), (2, 0)] {
        molecule
            .graph_mut()
            .add_bond(atoms[left], atoms[right], BondOrder::Single)
            .expect("triangle bond should be valid");
    }
    let original = molecule.clone();
    let ring_options = RingPerceptionOptions {
        max_path_expansions: 0,
        ..RingPerceptionOptions::default()
    };
    let error = perception_api::sanitize_with_ring_options(
        &mut molecule,
        SanitizeOptions {
            perceive_valence: false,
            perceive_rings: true,
            perceive_aromaticity: false,
        },
        ring_options,
    )
    .expect_err("ring limit should fail sanitization");
    assert!(matches!(error, SanitizeError::Rings(_)));
    assert_eq!(molecule, original);

    let mut aromatic = smiles_api::read_str_with_options("c1ccccc1", SmilesParseOptions)
        .expect("benzene should parse");
    let error = perception_api::perceive_aromaticity_with_ring_options(
        aromatic.graph_mut(),
        AromaticityModel::RdkitLike,
        RingPerceptionOptions {
            max_atoms: 0,
            ..RingPerceptionOptions::default()
        },
    )
    .expect_err("aromaticity should propagate ring limit");
    assert!(matches!(error, AromaticityError::RingPerception(_)));
}
