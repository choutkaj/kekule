use super::*;

#[test]
fn nitrogen_stereo_requires_constrained_unconjugated_pyramidal_geometry() {
    for (source, expected) in [
        ("CN(C)CC", 0),
        ("CN1CC1C", 1),
        ("N1CC1C", 1),
        ("CN1CCC1C", 0),
        ("CN1C(=O)C1C", 0),
        ("O=CN1CC1C", 0),
        ("C=CN1CC1C", 0),
        ("N#CN1CC1C", 0),
        ("c1ccccc1N2CC2C", 0),
        ("O=PN1CC1C", 1),
        ("O=SN1CC1C", 1),
        ("CC12CCN(CC1)C2", 1),
        ("C1CN2CCC1CC2", 1),
        ("C1C2CC3CC1N2C3", 0),
        ("C1CC2CCCN2C1", 0),
        ("C[N+](CC)(CCC)CCCC", 1),
    ] {
        let mut molecule = read_smiles(source).unwrap();
        for expanded in [false, true] {
            if expanded {
                molecule.perceive().unwrap();
                crate::hydrogens::add_hydrogens(&mut molecule).unwrap();
            }
            let before = molecule.clone();
            let count = stereo_api::detect_stereo_candidates(&molecule)
                .unwrap()
                .iter()
                .filter(|candidate| {
                    matches!(candidate,
                    StereoCandidate::Tetrahedral { center, .. }
                    if molecule.atom(*center).unwrap().element.symbol() == "N")
                })
                .count();
            assert_eq!(count, expected, "{source}, expanded={expanded}");
            assert_eq!(molecule, before);
        }
    }
}

#[test]
fn nitrogen_cleanup_distinguishes_inverting_centers_from_unclassified_radicals() {
    for (source, removed, unclassified) in [
        ("C[N@](CC)CCC", 1, 0),
        ("C[N@]1CC1C", 0, 0),
        ("O=C[N@]1CC1C", 1, 0),
        ("C[N@+]1CC1C", 0, 1),
    ] {
        let molecule = read_smiles(source).unwrap();
        let mut editor = molecule.edit();
        let report = stereo_api::cleanup_stereo(&mut editor, Default::default()).unwrap();
        assert_eq!(report.removed_elements.len(), removed, "{source}");
        assert_eq!(report.unclassified_elements.len(), unclassified, "{source}");
        assert_eq!(molecule.stereo_elements().count(), 1);
        assert_eq!(
            editor.finish().unwrap().stereo_elements().count(),
            1 - removed
        );
    }
}

#[test]
fn stereo_symmetry_resolves_branches_rings_and_configuration_dependencies() {
    for (source, tetrahedra, double_bonds) in [
        ("C(CC)(CC)(F)Cl", 0, 0),
        ("C([CH2:1][CH3:2])([CH2:3][CH3:4])(F)Cl", 0, 0),
        ("C(CC)(CO)(F)Cl", 1, 0),
        ("P(=O)(OCC)(OCC)N", 0, 0),
        ("O=P1(C)CCNCC1", 0, 0),
        ("CC1CCC(C)CC1", 2, 0),
        ("C1CN2CCC1CC2", 2, 0),
        ("CC1CCCCC1", 0, 0),
        ("C([C@H](F)Cl)([C@H](F)Cl)(Br)I", 2, 0),
        ("C([C@H](F)Cl)([C@@H](F)Cl)(Br)I", 3, 0),
        ("C(C(F)Cl)(C(F)Cl)(Br)I", 3, 0),
        ("C(C(CC)(CC)F)(C(CC)(CC)F)(Br)I", 0, 0),
        ("FC=C([C@H](F)Cl)[C@H](F)Cl", 2, 0),
        ("FC=C([C@H](F)Cl)[C@@H](F)Cl", 2, 1),
        ("CC/C(CC)=C(/F)Cl", 0, 0),
        ("CC1CCC(=C(F)F)CC1", 0, 0),
    ] {
        let mut molecule = read_smiles(source).unwrap();
        for expanded in [false, true] {
            if expanded {
                molecule.perceive().unwrap();
                crate::hydrogens::add_hydrogens(&mut molecule).unwrap();
            }
            let before = molecule.clone();
            let candidates = stereo_api::detect_stereo_candidates(&molecule).unwrap();
            assert_eq!(molecule, before, "read-only analysis: {source}");
            assert_eq!(
                candidates
                    .iter()
                    .filter(|s| matches!(s, StereoCandidate::Tetrahedral { .. }))
                    .count(),
                tetrahedra,
                "{source}, expanded={expanded}"
            );
            assert_eq!(
                candidates
                    .iter()
                    .filter(|s| matches!(s, StereoCandidate::DoubleBond { .. }))
                    .count(),
                double_bonds,
                "{source}, expanded={expanded}"
            );
        }
    }
}

#[test]
fn stereo_symmetry_does_not_enumerate_unconstrained_hydrogen_permutations() {
    for tail in ["CCC", "CCCCCC", "CCCCCCCCC"] {
        let source = format!("{tail}[N+]12CC[N+](CC1)(CC2){tail}");
        let molecule = read_smiles(&source).unwrap();
        let candidates = stereo_api::detect_stereo_candidates_with_options(
            &molecule,
            stereo_api::StereoPerceptionOptions {
                max_search_states: 20_000,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(candidates.len(), 2, "{source}");
    }
}

#[test]
fn stereo_symmetry_is_invariant_to_atom_and_bond_numbering() {
    use std::collections::{BTreeMap, BTreeSet};
    for source in [
        "C(C(CC)(CC)F)(C(CC)(CC)F)(Br)I",
        "CC1CCC(C)CC1",
        "C([C@H](F)Cl)([C@@H](F)Cl)(Br)I",
        "FC=C([C@H](F)Cl)[C@H](F)Cl",
        "FC=C([C@H](F)Cl)[C@@H](F)Cl",
        "CC12CCN(CC1)C2",
        "C1C2CC3CC1N2C3",
        "O=CN1CC1C",
    ] {
        let original = read_smiles(source).unwrap();
        let expected = stereo_api::detect_stereo_candidates(&original).unwrap();
        for shift in [1, 3] {
            let mut order: Vec<_> = original.atoms().collect();
            order.reverse();
            order.rotate_left(shift);
            let mut editor = MoleculeEditor::new();
            let atoms: BTreeMap<_, _> = order
                .into_iter()
                .map(|(id, atom)| (id, editor.add_atom(atom.clone()).unwrap()))
                .collect();
            let mut order: Vec<_> = original.bonds().collect();
            order.reverse();
            let bonds: BTreeMap<_, _> = order
                .into_iter()
                .map(|(id, bond)| {
                    (
                        id,
                        editor
                            .add_bond(atoms[&bond.a()], atoms[&bond.b()], bond.order)
                            .unwrap(),
                    )
                })
                .collect();
            let carrier = |c: StereoCarrier| match c {
                StereoCarrier::Atom(id) => StereoCarrier::Atom(atoms[&id]),
                other => other,
            };
            for (_, element) in original.stereo_elements() {
                let mut element = element.clone();
                match &mut element.kind {
                    StereoElementKind::Tetrahedral(s) => {
                        s.center = atoms[&s.center];
                        s.carriers = s.carriers.iter().copied().map(carrier).collect();
                    }
                    StereoElementKind::DoubleBond(s) => {
                        s.bond = bonds[&s.bond];
                        s.left = atoms[&s.left];
                        s.right = atoms[&s.right];
                        s.left_carrier = carrier(s.left_carrier);
                        s.right_carrier = carrier(s.right_carrier);
                    }
                    StereoElementKind::Axis(_) => unreachable!(),
                }
                editor.add_stereo_element(element).unwrap();
            }
            let reordered = editor.finish().unwrap();
            let expected: BTreeSet<_> = expected
                .iter()
                .map(|candidate| match candidate {
                    StereoCandidate::Tetrahedral { center, .. } => (0, atoms[center].raw()),
                    StereoCandidate::DoubleBond { bond, .. } => (1, bonds[bond].raw()),
                })
                .collect();
            let actual: BTreeSet<_> = stereo_api::detect_stereo_candidates(&reordered)
                .unwrap()
                .iter()
                .map(|candidate| match candidate {
                    StereoCandidate::Tetrahedral { center, .. } => (0, center.raw()),
                    StereoCandidate::DoubleBond { bond, .. } => (1, bond.raw()),
                })
                .collect();
            assert_eq!(actual, expected, "{source}, shift={shift}");
        }
    }
}

#[test]
fn stereo_cleanup_is_explicit_transactional_and_preserves_surviving_groups() {
    let mut molecule = read_smiles("F[C@](CC)(CC)C[C@H](Cl)Br |&1:1,7|").unwrap();
    molecule.perceive().unwrap();
    assert_eq!(
        molecule.stereo_elements().count(),
        2,
        "default perception preserves assertions"
    );
    let before = molecule.clone();
    let elements: Vec<_> = molecule
        .stereo_elements()
        .map(|(id, element)| (id, element.clone()))
        .collect();
    let mut editor = molecule.edit();
    for options in [
        stereo_api::StereoPerceptionOptions {
            max_atoms: 1,
            ..Default::default()
        },
        stereo_api::StereoPerceptionOptions {
            max_search_states: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            stereo_api::cleanup_stereo(&mut editor, options),
            Err(stereo_api::StereoPerceptionError::ResourceLimit { .. })
        ));
        assert_eq!(editor.working(), &before);
    }
    let report = stereo_api::cleanup_stereo(&mut editor, Default::default()).unwrap();
    assert_eq!(report.removed_elements, vec![elements[0].0]);
    let cleaned = editor.finish().unwrap();
    assert_eq!(cleaned.stereo_elements().count(), 1);
    assert_eq!(
        cleaned.stereo_element(elements[1].0).unwrap(),
        &elements[1].1
    );
    let group = cleaned.stereo_groups().next().unwrap().1;
    assert_eq!(group.kind, StereoGroupKind::And);
    assert_eq!(group.members, vec![elements[1].0]);
    assert!(cleaned.perception().has_valence());
    for id in before.atom_ids() {
        assert_eq!(
            cleaned.implicit_hydrogens(id).unwrap(),
            before.implicit_hydrogens(id).unwrap()
        );
    }
    assert_eq!(molecule, before);
}

#[test]
fn stereo_cleanup_removes_assertions_on_perceived_aromatic_bonds() {
    let molecule = read_smiles("C1=CC=CC=C1").unwrap();
    let mut editor = molecule.into_editor();
    let id = editor
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: BondId::new(0),
                left: AtomId::new(0),
                right: AtomId::new(1),
                left_carrier: StereoCarrier::Atom(AtomId::new(5)),
                right_carrier: StereoCarrier::Atom(AtomId::new(2)),
                orientation: Some(DoubleBondOrientation::Together),
            },
        )))
        .unwrap();
    let mut molecule = editor.finish().unwrap();
    molecule.perceive().unwrap();
    assert!(molecule.stereo_element(id).is_ok());
    let mut editor = molecule.edit();
    assert_eq!(
        stereo_api::cleanup_stereo(&mut editor, Default::default())
            .unwrap()
            .removed_elements,
        vec![id]
    );
    assert_eq!(editor.finish().unwrap().stereo_elements().count(), 0);
    assert_eq!(molecule.stereo_elements().count(), 1);
}

#[test]
fn stereo_cleanup_preserves_and_reports_unclassified_geometry() {
    let molecule = read_smiles("C([Sb@](F)Cl)([Sb@@](F)Cl)(Br)I").unwrap();
    let before = molecule.clone();
    let mut editor = molecule.edit();
    let result = stereo_api::cleanup_stereo(&mut editor, Default::default()).unwrap();
    assert!(result.removed_elements.is_empty());
    assert_eq!(result.unclassified_elements.len(), 2);
    assert_eq!(editor.finish().unwrap(), before);
    let candidates = stereo_api::detect_stereo_candidates(&molecule).unwrap();
    assert!(candidates.iter().any(|candidate| matches!(candidate,
        StereoCandidate::Tetrahedral { center, .. } if *center == AtomId::new(0))));
}
