use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor};
use kekule::smiles::{self, MolWriteErrorKind};

#[test]
fn canonical_smiles_matches_rdkit_charge_closures_stereo_and_cx_radicals() {
    for (source, expected) in [
        ("CC([O-])=O", "CC(=O)[O-]"),
        ("[O-][N+](=O)c1ccccc1", "O=[N+]([O-])c1ccccc1"),
        (
            "CC12CCC3C(CCC4CC(O)CCC43C)C2CCC1=O",
            "CC12CCC3C(CCC4CC(O)CCC43C)C1CCC2=O",
        ),
        ("C[C@H]1CCC[C@H]1C", "C[C@@H]1CCC[C@@H]1C"),
        ("[CH3].[CH2].[CH]", "[CH2].[CH3].[CH] |^1:1,^2:0,^5:2|"),
        ("O.[Cu].[Cu]", "O.[Cu].[Cu] |^1:1,2|"),
        (
            "[CH2][C@H](F)[C@@H](C)Cl |&1:1,3|",
            "[CH2][C@H](F)[C@@H](C)Cl |^1:0,&1:1,3|",
        ),
        (
            "C[C@H](F)[C@H](Cl)[C@H](Br)[C@H](O)C |&1:1,3,o2:5,7|",
            "C[C@H](O)[C@H](Br)[C@@H](Cl)[C@H](C)F |o1:1,3,&1:5,7|",
        ),
    ] {
        let write = |text: &str| {
            let mut molecules = smiles::to_molecules(text).unwrap();
            for molecule in &mut molecules {
                molecule.perceive().unwrap();
            }
            let topology = kekule::topology::Topology::from_molecules(molecules.clone()).unwrap();
            smiles::write(
                &topology,
                smiles::SmilesWriteOptions {
                    mode: smiles::SmilesWriteMode::Canonical,
                },
            )
            .unwrap()
        };
        let actual = write(source);
        assert_eq!(actual, expected, "{source}");
        assert_eq!(write(&actual), actual, "fixed point: {source}");
    }
}

#[test]
fn cx_radical_emission_preserves_occupancy_and_rejects_explicit_spin() {
    for source in ["[CH3]", "[CH2]", "[CH]", "[C]", "[Cu]", "[O-][O]"] {
        let mut original = smiles::to_molecules(source).unwrap().pop().unwrap();
        original.perceive().unwrap();
        for text in [
            smiles::write(&original, smiles::SmilesWriteOptions::default()).unwrap(),
            smiles::write(&original, smiles::SmilesWriteOptions::isomeric()).unwrap(),
            smiles::write(&original, smiles::SmilesWriteOptions::canonical()).unwrap(),
        ] {
            let restored = smiles::to_molecules(&text).unwrap().pop().unwrap();
            let radicals = |m: &Molecule| {
                let mut states = m
                    .atoms()
                    .map(|(_, a)| {
                        (
                            a.element.atomic_number(),
                            a.radical
                                .map(|r| (r.electron_count(), r.spin_multiplicity())),
                        )
                    })
                    .collect::<Vec<_>>();
                states.sort_unstable();
                states
            };
            assert_eq!(
                radicals(&original),
                radicals(&restored),
                "{source} -> {text}"
            );
        }
    }
    for source in ["[CH2] |^3:0|", "[CH2] |^4:0|", "[CH] |^6:0|", "[CH] |^7:0|"] {
        let molecule = smiles::to_molecules(source).unwrap().pop().unwrap();
        assert!(
            smiles::write(&molecule, smiles::SmilesWriteOptions::canonical())
                .unwrap_err()
                .message()
                .contains("spin")
        );
    }
}

#[test]
fn canonical_smiles_uses_rdkit_style_roots_and_branch_continuations() {
    // Focused regressions against RDKit 2026.03.3 canonical isomeric SMILES.
    // Alternate source traversals also guard against input-order tie breaking.
    for (sources, expected) in [
        (["CC(C)SC", "CSC(C)C"], "CSC(C)C"),
        (["CC(O)=O", "OC(C)=O"], "CC(=O)O"),
        (["C1=COC=C1", "o1cccc1"], "c1ccoc1"),
        (
            ["CC(=O)OC1=CC=CC=C1C(=O)O", "O=C(O)c1ccccc1OC(C)=O"],
            "CC(=O)Oc1ccccc1C(=O)O",
        ),
        (["C1=CC=C2C=CC=CC2=C1", "c1ccc2ccccc2c1"], "c1ccc2ccccc2c1"),
        (["C(CCF)COCCCCCl", "ClCCCCOCCCCF"], "FCCCCOCCCCCl"),
    ] {
        for source in sources {
            let mut molecule = smiles::to_molecules(source).unwrap().pop().unwrap();
            molecule.perceive().unwrap();
            let written =
                smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap();
            assert_eq!(written, expected, "{source}");
            let mut restored = smiles::to_molecules(&written).unwrap().pop().unwrap();
            restored.perceive().unwrap();
            assert_eq!(
                smiles::write(&restored, smiles::SmilesWriteOptions::canonical()).unwrap(),
                written
            );
        }
    }
}

#[test]
fn canonical_smiles_retains_charged_and_mapped_hydrogen_atoms() {
    for source in ["[H-][Na+]", "[H+]C", "[H:7]C", "[H-:7][Na+]"] {
        let original = smiles::to_molecules(source).unwrap().pop().unwrap();
        let written = smiles::write(&original, smiles::SmilesWriteOptions::canonical()).unwrap();
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
        let expected = smiles::write(
            &carbon_graph(count, edges, &permutation, false),
            smiles::SmilesWriteOptions::canonical(),
        )
        .unwrap();
        let restored = smiles::to_molecules(&expected).unwrap().pop().unwrap();
        assert_eq!(
            smiles::write(&restored, smiles::SmilesWriteOptions::canonical()).unwrap(),
            expected
        );
        let mut seed = 17u64;
        for iteration in 0..24 {
            for index in 0..count {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                permutation.swap(index, (seed as usize) % count);
            }
            let molecule = carbon_graph(count, edges, &permutation, iteration % 2 == 0);
            assert_eq!(
                smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap(),
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
        smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap()
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
            if source.contains("13")
                || source.contains("18")
                || source.contains("15")
                || source.contains("2H")
            {
                assert_ne!(written, expected, "isotope must distinguish {source}");
            } else {
                assert_eq!(written, expected, "{source}");
            }
            assert_eq!(canonical(&written), written, "{source} -> {written}");
        }
    }
}

#[test]
fn canonical_hydrogen_normalization_preserves_isotope_vertices() {
    let mut original = smiles::to_molecules("[H]C([3H])(F)Cl")
        .unwrap()
        .pop()
        .unwrap();
    original.perceive().unwrap();
    let written = smiles::write(&original, smiles::SmilesWriteOptions::canonical()).unwrap();
    let mut restored = smiles::to_molecules(&written).unwrap().pop().unwrap();
    restored.perceive().unwrap();
    assert_eq!(restored.atom_count(), 4);
    assert!(restored.atoms().any(|(_, atom)| atom.isotope == Some(3)));
    assert_eq!(
        smiles::write(&restored, smiles::SmilesWriteOptions::canonical()).unwrap(),
        written
    );
    let hydrogen_count = |molecule: &Molecule| {
        molecule
            .atoms()
            .map(|(id, atom)| {
                usize::from(atom.element.symbol() == "H")
                    + usize::from(atom.hydrogens.specified_count())
                    + usize::from(molecule.inferred_hydrogens(id).unwrap().unwrap_or(0))
            })
            .sum::<usize>()
    };
    assert_eq!(hydrogen_count(&original), 2);
    assert_eq!(hydrogen_count(&restored), hydrogen_count(&original));
}

#[test]
fn canonical_input_complexity_is_bounded_before_export() {
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
    let error = smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap_err();
    assert_eq!(error.kind(), MolWriteErrorKind::ResourceLimit);
    assert!(error.to_string().contains("input complexity"));
}

#[test]
fn canonical_stereo_is_invariant_under_atom_and_bond_permutations() {
    use kekule::core::{StereoCarrier, StereoElement, StereoElementKind};
    use std::collections::BTreeMap;
    let mut seed = 41u64;
    for source in [
        "CC(=O)Oc1ccccc1C(=O)O",
        "c1ccoc1",
        "c1ccc2ccccc2c1",
        "ClCCCCOCCCCF",
        "CC([O-])=O",
        "[O-][N+](=O)c1ccccc1",
        "CC12CCC3C(CCC4CC(O)CCC43C)C2CCC1=O",
        "C[C@H]1CCC[C@H]1C",
        "[CH2][C@H](F)[C@@H](C)Cl |&1:1,3|",
        "C[C@H](F)[C@H](Cl)[C@H](Br)[C@H](O)C |&1:1,3,o2:5,7|",
        "C[C@@H](F)[C@@H](Cl)[C@@H](Br)[C@@H](O)C |&1:1,3,o2:5,7|",
        "C[C@H](F)[C@H](Cl)[C@H](Br)[C@H](O)C |&1:1,3,r|",
        "N[C@@H](C)C(=O)O",
        "O[C@H]1CC[C@@H](O)CC1",
        "C[C@H](O)[C@@H](O)C",
        "C[C@H](O)[C@H](O)C",
        "F/C=C/C=C\\Cl",
        // PubChem CIDs 444596, 445186, 445554, 445925, 446439, 446467,
        // 448343 and 448437 exposed numbering-dependent carrier choices and
        // independently flipped redundant directional marks.
        "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CN)/C)/C",
        "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(C)C=C)/C",
        "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(C)CCN)/C",
        "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=CC(=O)O)C)/C",
        r"CC1=CCCC(C1/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2C(CCC=C2C)(C)C)\C)\C)/C)/C)(C)C",
        r"C/C=C(\C)/C=C\C=C(/C)\C=C\C1=C(CCCC1(C)C)C",
        r"CC(=CCC/C(=C/CC/C(=C/C=C/C(=C/C=C\C=C(/C)\C=C\C=C(/C)\C=C\C=C(/C)\C=C\CC(C)(C)OC)/C)/C)/C)C",
        r"CC1=C(C(C[C@H](C1)O)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(\C)/C=C/C=C(\C)/C=C/[C@H]2C(=C[C@@H](CC2(C)C)O)C)/C)/C",
        "F/C=C(/F)F",
        "C/C(Cl)=C(F)/C",
        "[13CH3][C@H]([2H])O",
        "F[P@](Cl)Br",
        "[H][C@@](F)(Cl)Br",
    ] {
        let mut original = smiles::to_molecules(source).unwrap().pop().unwrap();
        original.perceive().unwrap();
        let expected = smiles::write(&original, smiles::SmilesWriteOptions::canonical()).unwrap();
        let mut restored = smiles::to_molecules(&expected).unwrap().pop().unwrap();
        restored.perceive().unwrap();
        assert_eq!(
            smiles::write(&restored, smiles::SmilesWriteOptions::canonical()).unwrap(),
            expected
        );
        for iteration in 0..24 {
            let mut order = original.atom_ids().collect::<Vec<_>>();
            for index in (1..order.len()).rev() {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                order.swap(index, seed as usize % (index + 1));
            }
            let mut editor = MoleculeEditor::new();
            if iteration % 2 != 0 {
                // Deleted atom/bond slots must not enter the canonical order
                // or CX atom indices. Exercise sparse storage as well as a
                // permutation of densely allocated live atoms.
                let carbon = Atom::new(Element::from_symbol("C").unwrap());
                let left = editor.add_atom(carbon.clone()).unwrap();
                let right = editor.add_atom(carbon).unwrap();
                editor.add_bond(left, right, BondOrder::Single).unwrap();
                editor.delete_atom(left).unwrap();
                editor.delete_atom(right).unwrap();
            }
            let atoms = order
                .into_iter()
                .map(|old| {
                    (
                        old,
                        editor
                            .add_atom(original.atom(old).unwrap().clone())
                            .unwrap(),
                    )
                })
                .collect::<BTreeMap<_, _>>();
            let bonds = original
                .bonds()
                .collect::<Vec<_>>()
                .into_iter()
                .rev()
                .map(|(old, bond)| {
                    (
                        old,
                        editor
                            .add_bond(atoms[&bond.b()], atoms[&bond.a()], bond.order)
                            .unwrap(),
                    )
                })
                .collect::<BTreeMap<_, _>>();
            let remap = |carrier: &mut StereoCarrier| {
                if let StereoCarrier::Atom(id) = carrier {
                    *id = atoms[id];
                }
            };
            let mut elements = BTreeMap::new();
            for (id, element) in original.stereo_elements() {
                let mut kind = element.kind.clone();
                match &mut kind {
                    StereoElementKind::Tetrahedral(value) => {
                        value.center = atoms[&value.center];
                        value.carriers.iter_mut().for_each(remap);
                    }
                    StereoElementKind::DoubleBond(value) => {
                        value.bond = bonds[&value.bond];
                        value.left = atoms[&value.left];
                        value.right = atoms[&value.right];
                        remap(&mut value.left_carrier);
                        remap(&mut value.right_carrier);
                    }
                    StereoElementKind::Axis(_) => unreachable!(),
                }
                elements.insert(
                    id,
                    editor.add_stereo_element(StereoElement::new(kind)).unwrap(),
                );
            }
            for (_, group) in original.stereo_groups() {
                editor
                    .add_stereo_group(kekule::core::StereoGroup {
                        kind: group.kind,
                        members: group.members.iter().map(|id| elements[id]).collect(),
                    })
                    .unwrap();
            }
            let mut molecule = editor.finish().unwrap();
            molecule.perceive().unwrap();
            assert_eq!(
                smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap(),
                expected,
                "{source}: {atoms:?}"
            );
            assert_eq!(molecule.atom_count(), original.atom_count());
        }
    }
}

#[test]
fn canonical_topology_sorts_components_without_collapsing_instances() {
    let options = smiles::SmilesWriteOptions {
        mode: smiles::SmilesWriteMode::Canonical,
    };
    let canonical = |source: &str| {
        let mut molecules = smiles::to_molecules(source).unwrap();
        for molecule in &mut molecules {
            molecule.perceive().unwrap();
        }
        let topology = kekule::topology::Topology::from_molecules(molecules.clone()).unwrap();
        smiles::write(&topology, options).unwrap()
    };
    assert_eq!(canonical("O.CC.O.[13CH4]"), canonical("[13CH4].O.O.CC"));
    assert_eq!(canonical("O.CC.O.[13CH4]").split('.').count(), 4);
}

#[test]
fn canonical_aromatic_output_does_not_depend_on_localized_double_bonds() {
    let canonical = |source: &str| {
        let mut molecule = smiles::to_molecules(source).unwrap().pop().unwrap();
        molecule.perceive().unwrap();
        smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap()
    };
    for equivalents in [
        ["CC1=CC=CC=C1O", "Cc1ccccc1O", "CC1=C(O)C=CC=C1"],
        [
            "CC1=NC2=C(C=CC=C2)C(=O)N1C",
            "Cc1nc2ccccc2c(=O)n1C",
            "Cn1c(C)nc2ccccc2c1=O",
        ],
    ] {
        let expected = canonical(equivalents[0]);
        for source in equivalents {
            assert_eq!(canonical(source), expected, "{source}");
        }
        assert_eq!(canonical(&expected), expected);
    }
}

#[test]
fn canonical_retains_isolated_and_molecular_hydrogen_without_perception() {
    for source in ["[H]", "[H][H]", "[2H][H]"] {
        let molecule = smiles::to_molecules(source).unwrap().pop().unwrap();
        let written = smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap();
        let restored = smiles::to_molecules(&written).unwrap().pop().unwrap();
        assert_eq!(restored.atom_count(), molecule.atom_count());
        assert_eq!(
            smiles::write(&restored, smiles::SmilesWriteOptions::canonical()).unwrap(),
            written
        );
    }
}

#[test]
fn smiles_uses_all_one_hundred_ring_labels_and_bounds_live_labels() {
    let complete = (0..21)
        .flat_map(|left| ((left + 1)..21).map(move |right| (left, right)))
        .collect::<Vec<_>>();
    let bounded = complete
        .iter()
        .copied()
        .filter(|&(left, right)| !(left == 0 && (11..20).contains(&right)))
        .collect::<Vec<_>>();
    let order = (0..21).collect::<Vec<_>>();
    let molecule = carbon_graph(21, &bounded, &order, false);
    let written = smiles::write(&molecule, smiles::SmilesWriteOptions::default()).unwrap();
    assert!(written.contains('0'), "{written}");
    let restored = smiles::to_molecules(&written).unwrap().pop().unwrap();
    assert_eq!(restored.bond_count(), molecule.bond_count());
    let error = smiles::write(
        &carbon_graph(21, &complete, &order, false),
        smiles::SmilesWriteOptions::default(),
    )
    .unwrap_err();
    assert_eq!(error.kind(), MolWriteErrorKind::ResourceLimit);
}

#[test]
fn canonical_labeling_prunes_equivalent_phenyl_branches() {
    let source = "COC1C(C(C(C(O1)CBr)OC(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)OC(C5=CC=CC=C5)(C6=CC=CC=C6)C7=CC=CC=C7)OC(C8=CC=CC=C8)(C9=CC=CC=C9)C1=CC=CC=C1";
    let mut molecule = smiles::to_molecules(source).unwrap().pop().unwrap();
    molecule.perceive().unwrap();
    let written = smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap();
    let mut restored = smiles::to_molecules(&written).unwrap().pop().unwrap();
    restored.perceive().unwrap();
    assert_eq!(restored.atom_count(), molecule.atom_count());
    assert_eq!(restored.bond_count(), molecule.bond_count());
    assert_eq!(
        smiles::write(&restored, smiles::SmilesWriteOptions::canonical()).unwrap(),
        written
    );
}

#[test]
fn source_order_spiro_groups_preserve_native_stereo_on_rereading() {
    // Supplied Enamine cases whose RDKit canonical identity strings depend on
    // traversal. RDKit 2026.03.3 confirms mutual chirality/enhanced-group-aware
    // matches for every source/output pair; native reconstruction must agree too.
    for source in [
        r"O[C@@H]1C[C@@]2(C[C@H](C2)NC=3N=CC=C(N3)C=4C=CC=NC4)[C@H]1O |&1:1,3,5,20,r| ID:Z9237738672",
        r"COC(=O)[C@H]1C(F)(F)[C@@]21C[C@H](C2)C3=NC(=NO3)C4=NC(C)=NN4C |&1:4,8,10,r| ID:Z5188334346",
        r"O=C(O)[C@H]1C[C@@]21C[C@H](C2)NC(=O)C3=CC=C(COC=4C=CC=CC4)O3 |&1:3,5,7,r| ID:Z3510155828",
        r"COC(=O)[C@H]1C(F)(F)[C@@]21C[C@H](C2)C3=NC(=NO3)C(C)C4=CC=CO4 |&1:4,8,10,r| ID:Z4875774646",
        r"O=C(O)[C@H]1C[C@@]21C[C@H](C2)NC(=O)C=3C=CC=C(C3)OCC4=CN5C=CC=CC5=N4 |&1:3,5,7,r| ID:Z3510152401",
    ] {
        let mut molecule = smiles::to_molecules(source).unwrap().pop().unwrap();
        molecule.perceive().unwrap();
        let expected = smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap();
        let written = smiles::write(&molecule, smiles::SmilesWriteOptions::isomeric()).unwrap();
        let mut restored = smiles::to_molecules(&written).unwrap().pop().unwrap();
        restored.perceive().unwrap();
        assert_eq!(
            smiles::write(&restored, smiles::SmilesWriteOptions::canonical()).unwrap(),
            expected,
            "{source} -> {written}"
        );
    }
}
