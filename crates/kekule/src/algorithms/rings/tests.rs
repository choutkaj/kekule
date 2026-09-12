use super::*;

#[test]
fn ring_limits_include_deleted_storage_before_allocating_scratch_arrays() {
    let mut editor = graph(4, &[(0, 1), (1, 2), (2, 0), (0, 3)]).into_editor();
    editor.delete_atom(AtomId::new(3)).unwrap();
    let mut molecule = editor.finish().unwrap();
    perceive_ring_set(&mut molecule).unwrap();
    let previous = molecule.perception().clone();
    for (options, resource, observed, limit) in [
        (
            RingPerceptionOptions {
                max_atoms: 3,
                ..Default::default()
            },
            "atoms",
            4,
            3,
        ),
        (
            RingPerceptionOptions {
                max_bonds: 3,
                ..Default::default()
            },
            "bonds",
            4,
            3,
        ),
        (
            RingPerceptionOptions {
                max_total_work: 7,
                ..Default::default()
            },
            "total work",
            8,
            7,
        ),
    ] {
        assert_eq!(
            perceive_ring_set_with_options(&mut molecule, options),
            Err(RingPerceptionError::ResourceLimit {
                resource,
                observed,
                limit
            })
        );
        assert_eq!(molecule.perception(), &previous);
    }
}

fn graph(size: usize, edges: &[(usize, usize)]) -> Molecule {
    let mut editor = MoleculeEditor::new();
    let atoms = (0..size)
        .map(|_| {
            editor
                .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
                .unwrap()
        })
        .collect::<Vec<_>>();
    for &(left, right) in edges {
        editor
            .add_bond(atoms[left], atoms[right], BondOrder::Single)
            .unwrap();
    }
    editor.finish().unwrap()
}

fn atom_sets(rings: &RingSet) -> Vec<Vec<usize>> {
    let mut sets = rings
        .rings()
        .iter()
        .map(|ring| {
            let mut atoms = ring
                .atoms
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>();
            atoms.sort();
            atoms
        })
        .collect::<Vec<_>>();
    sets.sort();
    sets
}

#[test]
fn short_ring_queries_use_the_default_ring_bond_policy() {
    for order in [BondOrder::Single, BondOrder::Zero, BondOrder::Dative] {
        let mut editor = graph(3, &[(0, 1), (1, 2)]).into_editor();
        let closure = editor
            .add_bond(AtomId::new(2), AtomId::new(0), order)
            .unwrap();
        let molecule = editor.finish().unwrap();
        let expected = order == BondOrder::Single;
        for bond in [BondId::new(0), BondId::new(1), closure] {
            assert_eq!(bond_in_ring_smaller_than(&molecule, bond, 4), expected);
            assert!(!bond_in_ring_smaller_than(&molecule, bond, 3));
        }
    }
}

#[test]
fn high_degree_graph_uses_rdkit_depth_first_fallback() {
    let edges = (0..5)
        .flat_map(|left| (left + 1..5).map(move |right| (left, right)))
        .collect::<Vec<_>>();
    let mut molecule = graph(5, &edges);
    let rings = perceive_ring_set(&mut molecule).unwrap();
    // RDKit 2026.03.3 falls back when every remaining vertex has degree four.
    assert_eq!(
        atom_sets(&rings),
        vec![
            vec![0, 1, 2],
            vec![0, 1, 2, 3],
            vec![0, 1, 2, 3, 4],
            vec![1, 2, 3],
            vec![1, 2, 3, 4],
            vec![2, 3, 4],
        ]
    );
    assert!(molecule
        .ring_membership()
        .unwrap()
        .bond_slot_flags()
        .iter()
        .all(|flag| *flag));
}

#[test]
fn fallback_preserves_coverage_when_rdkit_selected_rings_omit_a_cyclic_bond() {
    let edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (0, 9),
        (0, 3),
        (4, 9),
        (2, 7),
        (3, 8),
        (0, 6),
        (1, 8),
        (2, 9),
        (5, 7),
    ];
    let mut molecule = graph(10, &edges);
    let omitted = BondId::new(14);
    // RDKit 2026.03.3 GetSymmSSSR omits edge 0--6. The independent
    // alternative path 0--3--4--5--6 proves that this edge is cyclic.
    for (left, right) in [(0, 3), (3, 4), (4, 5), (5, 6)] {
        assert!(molecule
            .bond_between(AtomId::new(left), AtomId::new(right))
            .unwrap()
            .is_some());
    }
    let rings = perceive_ring_set(&mut molecule).unwrap();
    assert!(molecule.ring_membership().unwrap().bond_in_ring(omitted));
    assert!(rings
        .rings()
        .iter()
        .any(|ring| ring.bonds.contains(&omitted)));
    for bond in molecule.ring_membership().unwrap().ring_bond_ids() {
        assert!(rings.rings().iter().any(|ring| ring.bonds.contains(&bond)));
    }
}

#[test]
fn fallback_work_limit_preserves_installed_perception() {
    let edges = (0..5)
        .flat_map(|left| (left + 1..5).map(move |right| (left, right)))
        .collect::<Vec<_>>();
    let mut molecule = graph(5, &edges);
    perceive_ring_membership(&mut molecule);
    let before = molecule.perception().clone();
    let options = RingPerceptionOptions {
        max_total_work: 20,
        ..RingPerceptionOptions::default()
    };
    assert!(matches!(
        perceive_ring_set_with_options(&mut molecule, options),
        Err(RingPerceptionError::ResourceLimit {
            resource: "total work",
            ..
        })
    ));
    assert_eq!(molecule.perception(), &before);
}

#[test]
fn candidate_count_at_cyclomatic_number_is_not_pruned() {
    let mut molecule = graph(
        7,
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (1, 5),
            (1, 6),
            (0, 3),
            (0, 4),
            (0, 6),
            (2, 5),
        ],
    );
    let rings = perceive_ring_set(&mut molecule).unwrap();
    assert_eq!(
        atom_sets(&rings),
        vec![
            vec![0, 1, 2, 3],
            vec![0, 1, 4, 5],
            vec![0, 1, 6],
            vec![0, 3, 4],
            vec![0, 4, 5, 6],
            vec![1, 2, 5],
            vec![1, 5, 6],
            vec![2, 3, 4, 5],
        ]
    );
}

#[test]
fn pubchem_125634_ring_scaffold_keeps_the_central_six_membered_ring() {
    // Ring-only topology from externally supplied PubChem CID 125634. Removing
    // acyclic atoms and compacting IDs retains its nine RDKit-selected rings.
    let mut molecule = graph(
        28,
        &[
            (0, 2),
            (0, 4),
            (0, 6),
            (0, 12),
            (1, 3),
            (1, 5),
            (1, 7),
            (1, 13),
            (2, 3),
            (2, 10),
            (2, 14),
            (3, 11),
            (3, 15),
            (4, 5),
            (4, 8),
            (5, 9),
            (6, 9),
            (6, 11),
            (7, 8),
            (7, 10),
            (12, 16),
            (13, 17),
            (14, 18),
            (15, 19),
            (16, 18),
            (16, 20),
            (17, 19),
            (17, 21),
            (18, 22),
            (19, 23),
            (20, 24),
            (21, 25),
            (22, 26),
            (23, 27),
            (24, 26),
            (25, 27),
        ],
    );
    let rings = perceive_ring_set(&mut molecule).unwrap();
    assert_eq!(
        atom_sets(&rings),
        vec![
            vec![0, 1, 2, 3, 4, 5],
            vec![0, 2, 3, 6, 11],
            vec![0, 2, 12, 14, 16, 18],
            vec![0, 4, 5, 6, 9],
            vec![1, 2, 3, 7, 10],
            vec![1, 3, 13, 15, 17, 19],
            vec![1, 4, 5, 7, 8],
            vec![16, 18, 20, 22, 24, 26],
            vec![17, 19, 21, 23, 25, 27],
        ]
    );
}

#[test]
fn degree_three_recovery_preserves_rdkit_candidate_order() {
    let molecule = graph(
        6,
        &[
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 2),
            (1, 4),
            (4, 3),
            (2, 5),
            (5, 3),
        ],
    );
    let active = ActiveRingGraph::new(&molecule);
    let mut tracker = RingWorkTracker::new(RingPerceptionOptions::default(), 6, 8).unwrap();
    let mut rings = Vec::new();
    find_rings_from_degree_three_node(
        AtomId::new(0),
        &active,
        &mut rings,
        &mut BTreeSet::new(),
        &mut tracker,
    )
    .unwrap();
    let order = rings
        .iter()
        .map(|ring| {
            let mut atoms = ring
                .atoms
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>();
            atoms.sort();
            atoms
        })
        .collect::<Vec<_>>();
    assert_eq!(
        order,
        vec![vec![0, 1, 2], vec![0, 1, 3, 4], vec![0, 2, 3, 5]]
    );
}

#[test]
fn symmetrization_only_uses_the_original_sssr_as_replacement_witnesses() {
    let ring = |bonds: &[u32]| Ring {
        atoms: Vec::new(),
        bonds: bonds.iter().copied().map(BondId::new).collect(),
    };
    let original = vec![ring(&[0, 1, 2]), ring(&[0, 1, 3, 4]), ring(&[2, 3, 5, 6])];
    let accepted = ring(&[0, 1, 3]);
    let chained = ring(&[3, 4, 5]);
    let mut tracker = RingWorkTracker::new(RingPerceptionOptions::default(), 0, 7).unwrap();
    let result = symmetrize_ring_set(
        original.clone(),
        vec![accepted.clone(), chained],
        7,
        &mut tracker,
    )
    .unwrap();
    let mut expected = original;
    expected.push(accepted);
    assert_eq!(result, expected);
}

#[test]
fn pruning_and_symmetrization_charge_comparison_work() {
    let rings = (0..20)
        .map(|index| Ring {
            atoms: Vec::new(),
            bonds: (index * 3..index * 3 + 3).map(BondId::new).collect(),
        })
        .collect::<Vec<_>>();
    let options = RingPerceptionOptions {
        max_total_work: 100,
        ..RingPerceptionOptions::default()
    };
    let mut tracker = RingWorkTracker::new(options, 0, 0).unwrap();
    assert!(matches!(
        remove_extra_rings(rings.clone(), 60, &mut tracker),
        Err(RingPerceptionError::ResourceLimit {
            resource: "total work",
            ..
        })
    ));
    let mut tracker = RingWorkTracker::new(options, 0, 0).unwrap();
    let extra = Ring {
        atoms: Vec::new(),
        bonds: vec![BondId::new(0), BondId::new(3), BondId::new(6)],
    };
    assert!(matches!(
        symmetrize_ring_set(rings, vec![extra], 60, &mut tracker),
        Err(RingPerceptionError::ResourceLimit {
            resource: "total work",
            ..
        })
    ));
}

#[test]
fn search_and_recovery_charge_setup_before_visiting_paths() {
    let molecule = graph(3, &[(0, 1), (1, 2), (2, 0)]);
    let active = ActiveRingGraph::new(&molecule);
    let tracker_with_budget = |max_total_work| {
        RingWorkTracker::new(
            RingPerceptionOptions {
                max_total_work,
                ..RingPerceptionOptions::default()
            },
            0,
            0,
        )
        .unwrap()
    };

    let mut tracker = tracker_with_budget(8);
    assert_eq!(
        smallest_rings_bfs(AtomId::new(0), &active, &BTreeSet::new(), &mut tracker),
        Err(RingPerceptionError::ResourceLimit {
            resource: "total work",
            observed: 9,
            limit: 8,
        })
    );
    assert_eq!(tracker.path_expansions, 0);

    let mut tracker = tracker_with_budget(14);
    let mut candidates = Vec::new();
    let mut invariants = BTreeSet::new();
    let duplicate_roots = BTreeMap::from([(
        vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)],
        vec![AtomId::new(0), AtomId::new(1)],
    )]);
    assert_eq!(
        recover_duplicate_degree_two_candidates(
            &active,
            &duplicate_roots,
            &BTreeMap::new(),
            &mut candidates,
            &mut invariants,
            &mut tracker,
        ),
        Err(RingPerceptionError::ResourceLimit {
            resource: "total work",
            observed: 15,
            limit: 14,
        })
    );
    assert_eq!(tracker.path_expansions, 0);
    assert!(candidates.is_empty());
    assert!(invariants.is_empty());

    let mut tracker = tracker_with_budget(2);
    let mut discovered = RingMembership {
        atom_flags: vec![false; 3],
        bond_flags: vec![false; 3],
    };
    assert_eq!(
        recover_connecting_cycles(
            &molecule,
            &[AtomId::new(0), AtomId::new(1), AtomId::new(2)],
            &mut candidates,
            &mut invariants,
            &mut discovered,
            &mut tracker,
        ),
        Err(RingPerceptionError::ResourceLimit {
            resource: "total work",
            observed: 3,
            limit: 2,
        })
    );
    assert_eq!(tracker.path_expansions, 0);
}

#[test]
fn connecting_recovery_charges_every_bond_scan() {
    let molecule = graph(3, &[(0, 1), (1, 2), (2, 0)]);
    let atoms = vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)];
    let mut invariants = BTreeSet::from([atoms.clone()]);
    let mut candidates = Vec::new();
    let mut discovered = RingMembership {
        atom_flags: vec![true; 3],
        bond_flags: vec![false; 3],
    };
    let mut tracker = RingWorkTracker::new(RingPerceptionOptions::default(), 0, 0).unwrap();
    recover_connecting_cycles(
        &molecule,
        &atoms,
        &mut candidates,
        &mut invariants,
        &mut discovered,
        &mut tracker,
    )
    .unwrap();
    // Each existing-invariant candidate is exhausted, followed by a final scan
    // proving no candidate remains. All four scans must consume the budget.
    assert!(tracker.total_work >= 4 * molecule.bond_count() + tracker.path_expansions);
    assert!(candidates.is_empty());
}

#[test]
fn degree_two_roots_preserve_fragment_order_without_revisiting_groups() {
    let molecule = graph(
        8,
        &[
            (0, 2),
            (2, 3),
            (3, 1),
            (0, 4),
            (4, 5),
            (5, 1),
            (0, 6),
            (6, 7),
            (7, 1),
        ],
    );
    let active = ActiveRingGraph::new(&molecule);
    let fragment = [5, 4, 3, 2, 7, 6, 1, 0].map(AtomId::new);
    assert_eq!(
        pick_degree_two_nodes(&fragment, &active),
        vec![AtomId::new(5), AtomId::new(3), AtomId::new(7)]
    );
}
