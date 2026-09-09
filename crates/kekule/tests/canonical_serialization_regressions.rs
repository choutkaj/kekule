use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor};
use kekule::smiles;

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
