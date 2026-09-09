use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor};
use kekule::smiles::{self, MolWriteErrorKind};

#[test]
fn canonical_smiles_retains_charged_and_mapped_hydrogen_atoms() {
    for source in ["[H-][Na+]", "[H+]C", "[H:7]C", "[H-:7][Na+]"] {
        let original = smiles::to_molecules(source).unwrap().pop().unwrap();
        let written = smiles::write_canonical(&original).unwrap();
        let restored = smiles::to_molecules(&written).unwrap().pop().unwrap();
        assert_eq!(
            restored.formal_charge(),
            original.formal_charge(),
            "{source} -> {written}"
        );
        let asserted_hydrogens = |molecule: &Molecule| {
            molecule
                .atoms()
                .filter(|(_, atom)| atom.element.symbol() == "H")
                .map(|(_, atom)| (atom.formal_charge, atom.atom_map))
                .collect::<Vec<_>>()
        };
        assert_eq!(
            asserted_hydrogens(&restored),
            asserted_hydrogens(&original),
            "{source} -> {written}"
        );
    }
}

fn carbon_graph(
    count: usize,
    edges: &[(usize, usize)],
    permutation: &[usize],
    reverse: bool,
) -> Molecule {
    let mut editor = MoleculeEditor::new();
    let atoms = (0..count)
        .map(|_| {
            editor
                .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
                .unwrap()
        })
        .collect::<Vec<_>>();
    for index in 0..edges.len() {
        let (left, right) = edges[if reverse {
            edges.len() - 1 - index
        } else {
            index
        }];
        let (left, right) = if reverse {
            (right, left)
        } else {
            (left, right)
        };
        editor
            .add_bond(
                atoms[permutation[left]],
                atoms[permutation[right]],
                BondOrder::Single,
            )
            .unwrap();
    }
    editor.finish().unwrap()
}

#[test]
fn canonical_smiles_is_invariant_under_atom_bond_and_endpoint_permutations() {
    // Focused graph regressions, not external chemistry benchmark fixtures.
    let graphs: &[(usize, &[(usize, usize)])] = &[
        (
            6,
            &[
                (0, 3),
                (0, 4),
                (0, 5),
                (1, 3),
                (1, 4),
                (1, 5),
                (2, 3),
                (2, 4),
                (2, 5),
            ],
        ),
        (
            10,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 0),
                (0, 5),
                (1, 6),
                (2, 7),
                (3, 8),
                (4, 9),
                (5, 7),
                (7, 9),
                (9, 6),
                (6, 8),
                (8, 5),
            ],
        ),
        (
            8,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 0),
                (4, 5),
                (5, 6),
                (6, 7),
                (7, 4),
                (0, 4),
                (1, 5),
                (2, 6),
                (3, 7),
            ],
        ),
        (
            8,
            &[
                (0, 2),
                (0, 3),
                (1, 2),
                (1, 3),
                (2, 3),
                (4, 6),
                (4, 7),
                (5, 6),
                (5, 7),
                (6, 7),
                (0, 4),
                (1, 5),
            ],
        ),
    ];
    for &(count, edges) in graphs {
        let mut permutation = (0..count).collect::<Vec<_>>();
        let expected =
            smiles::write_canonical(&carbon_graph(count, edges, &permutation, false)).unwrap();
        let restored = smiles::to_molecules(&expected).unwrap().pop().unwrap();
        assert_eq!(smiles::write_canonical(&restored).unwrap(), expected);
        let mut seed = 17u64;
        for iteration in 0..24 {
            for index in 0..count {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                permutation.swap(index, (seed as usize) % count);
            }
            let molecule = carbon_graph(count, edges, &permutation, iteration % 2 == 0);
            assert_eq!(
                smiles::write_canonical(&molecule).unwrap(),
                expected,
                "permutation {permutation:?}"
            );
        }
    }
}

#[test]
fn canonical_smiles_ranks_the_emitted_hydrogen_and_isotope_projection() {
    let canonical = |source: &str| {
        let mut molecule = smiles::to_molecules(source).unwrap().pop().unwrap();
        molecule.perceive().unwrap();
        smiles::write_canonical(&molecule).unwrap()
    };
    for equivalents in [
        ["CC(C)CC", "CC(C)[CH2]C", "CC(C)[13CH2]C"],
        ["CC(=O)OC(C)CN", "CC(=[18O])OC(C)CN", "CC(=O)O[CH](C)CN"],
        ["C1CC2CCC1C2", "[13CH2]1CC2CCC1C2", "[CH2]1CC2CCC1C2"],
        ["c1ccccc1", "[H]c1ccccc1", "[2H]c1ccccc1"],
        ["N1C=CC=C1", "[NH]1C=CC=C1", "[15NH]1C=CC=C1"],
        [
            "CC1=CC(C)=CC=C1O",
            "C[13C]1=CC(C)=CC=C1O",
            "CC1=CC(C)=CC=[13C]1O",
        ],
    ] {
        let expected = canonical(equivalents[0]);
        for source in equivalents {
            let written = canonical(source);
            assert_eq!(written, expected, "{source}");
            assert_eq!(canonical(&written), written, "{source} -> {written}");
        }
    }
}

#[test]
fn canonical_isotope_free_hydrogen_projection_is_stable() {
    let mut original = smiles::to_molecules("[H]C([3H])(F)Cl")
        .unwrap()
        .pop()
        .unwrap();
    original.perceive().unwrap();
    let written = smiles::write_canonical(&original).unwrap();
    let mut restored = smiles::to_molecules(&written).unwrap().pop().unwrap();
    restored.perceive().unwrap();
    assert_eq!(restored.atom_count(), 3);
    assert_eq!(smiles::write_canonical(&restored).unwrap(), written);
    let hydrogen_count = |molecule: &Molecule| {
        molecule
            .atoms()
            .map(|(id, atom)| {
                usize::from(atom.element.symbol() == "H")
                    + usize::from(atom.hydrogens.explicit_count())
                    + usize::from(molecule.implicit_hydrogens(id).unwrap().unwrap_or(0))
            })
            .sum::<usize>()
    };
    assert_eq!(hydrogen_count(&original), 2);
    assert_eq!(hydrogen_count(&restored), hydrogen_count(&original));
}

#[test]
fn canonical_candidate_work_is_bounded_before_export() {
    let mut editor = MoleculeEditor::new();
    let mut previous = None;
    for index in 0..4096 {
        let mut atom = Atom::new(Element::from_symbol("C").unwrap());
        atom.atom_map = Some(index + 1);
        let atom = editor.add_atom(atom).unwrap();
        if let Some(previous) = previous {
            editor.add_bond(previous, atom, BondOrder::Single).unwrap();
        }
        previous = Some(atom);
    }
    let molecule = editor.finish().unwrap();
    let error = smiles::write_canonical(&molecule).unwrap_err();
    assert_eq!(error.kind(), MolWriteErrorKind::ResourceLimit);
    assert!(error.to_string().contains("candidate traversal"));
}
