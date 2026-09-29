use kekule::{
    core::*,
    perception::{self, conjugation::*, resonance::*},
    smiles,
};

fn molecule(smi: &str) -> Molecule {
    let mut m = smiles::to_molecules(smi).unwrap().remove(0);
    m.perceive().unwrap();
    m
}

#[test]
fn conjugation_matches_pinned_rdkit_local_rules() {
    for (smi, expected) in [
        ("c1ccccc1", 6),
        ("C=CC=C", 3),
        ("CC(=O)[O-]", 2),
        ("C[N+](=O)[O-]", 2),
        ("Nc1ccccc1", 7),
        ("C=CC", 0),
        ("C=C=C", 2),
        ("[CH2+]C=C", 2),
        ("[CH2-]C=C", 2),
        ("CSC=C", 0),
        ("CC(=O)N", 2),
    ] {
        let m = molecule(smi);
        assert_eq!(
            m.perception().conjugation_state().unwrap().bonds().len(),
            expected,
            "{smi}"
        );
        assert!(!m.perception().has_resonance());
    }
}

#[test]
fn default_contributors_and_all_kekule_forms() {
    for (smi, expected) in [
        ("c1ccccc1", 1),
        ("C=CC=C", 1),
        ("CC(=O)[O-]", 2),
        ("C[N+](=O)[O-]", 2),
        ("CC(=O)N", 1),
        ("[CH2+]C=C", 2),
        ("[CH2-]C=C", 2),
        ("C=C=C", 1),
        ("C", 1),
    ] {
        let m = molecule(smi);
        let result = enumerate_resonance(&m, Default::default()).unwrap();
        assert_eq!(result.contributors().len(), expected, "{smi}");
        for i in 0..expected {
            let c = result.to_molecule(i).unwrap();
            assert_eq!(c.atom_count(), m.atom_count());
            assert!(!c.perception().has_conjugation());
        }
    }
    let m = molecule("c1ccccc1");
    let all = enumerate_resonance(
        &m,
        ResonanceOptions {
            flags: ResonanceFlags::KEKULE_ALL,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(all.contributors().len(), 2);
}

#[test]
fn group_preparation_is_separate_and_transactional() {
    let mut m = molecule("CC(=O)NCC(=O)[O-]");
    let original = m.clone();
    perceive_resonance(&mut m).unwrap();
    assert_eq!(m.perception().resonance_state().unwrap().groups().len(), 2);
    assert_eq!(m, original);
    let before = m.perception().clone();
    assert!(matches!(
        enumerate_resonance(
            &m,
            ResonanceOptions {
                max_total_work: 0,
                ..Default::default()
            }
        ),
        Err(ResonanceError::ResourceLimit { .. })
    ));
    assert_eq!(m.perception(), &before);
    perception::aromaticity::perceive_aromaticity(&mut m, AromaticityModel::Mdl).unwrap();
    assert!(!m.perception().has_conjugation());
    assert!(!m.perception().has_resonance());
    perceive_conjugation(&mut m, ConjugationModel::RdkitLike).unwrap();
    assert!(m.perception().has_conjugation());
}

#[test]
fn missing_prerequisites_do_not_install_empty_state() {
    let mut m = smiles::to_molecules("CC").unwrap().remove(0);
    let before = m.perception().clone();
    assert_eq!(
        perceive_conjugation(&mut m, ConjugationModel::RdkitLike),
        Err(ConjugationError::MissingAromaticity)
    );
    assert_eq!(
        perceive_resonance(&mut m),
        Err(ResonanceError::MissingConjugation)
    );
    assert_eq!(m.perception(), &before);
}

#[test]
fn all_flag_combinations_match_pinned_reference_counts() {
    // Independently evaluated with RDKit 2026.03.3. Keep all masks, including
    // redundant combinations, to exercise the implications of both ion flags.
    for (smi, counts) in [
        (
            "c1ccccc1",
            [
                1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5,
                6, 6, 6, 6,
            ],
        ),
        (
            "CC(=O)N",
            [
                1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2,
            ],
        ),
        (
            "C[N+](=O)[O-]",
            [
                2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4,
                4, 4, 4, 4,
            ],
        ),
        (
            "C=CC#N",
            [
                1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
                2, 2, 2, 2,
            ],
        ),
    ] {
        let m = molecule(smi);
        for (mask, count) in counts.into_iter().enumerate() {
            let forms = enumerate_resonance(
                &m,
                ResonanceOptions {
                    flags: ResonanceFlags::from_bits(mask as u8).unwrap(),
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(forms.contributors().len(), count, "{smi}, flags {mask}");
        }
    }
    assert_eq!(ResonanceFlags::from_bits(32), None);
}

#[test]
fn localized_contributors_keep_atom_identity_and_canonical_seed() {
    let m = molecule("CC(=O)[O-]");
    let forms = enumerate_resonance(&m, Default::default()).unwrap();
    let mut signatures: Vec<_> = forms
        .contributors()
        .iter()
        .map(|c| {
            (
                c.formal_charges().map(|(_, q)| q).collect::<Vec<_>>(),
                c.bond_orders().map(|(_, o)| o).collect::<Vec<_>>(),
            )
        })
        .collect();
    signatures.sort_by(|a, b| a.0.cmp(&b.0));
    assert_eq!(
        signatures,
        vec![
            (
                vec![0, 0, -1, 0],
                vec![BondOrder::Single, BondOrder::Single, BondOrder::Double]
            ),
            (
                vec![0, 0, 0, -1],
                vec![BondOrder::Single, BondOrder::Double, BondOrder::Single]
            )
        ]
    );
    let benzene = molecule("c1ccccc1");
    let forms = enumerate_resonance(&benzene, Default::default()).unwrap();
    assert_eq!(
        forms.contributors()[0]
            .bond_orders()
            .map(|(_, o)| o)
            .collect::<Vec<_>>(),
        vec![
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double
        ]
    );
}

#[test]
fn canonical_aromatic_seed_retains_aryne_triple_bonds() {
    let m = molecule("C1#CC=CC=C1");
    let triple = m
        .bonds()
        .find(|(_, b)| b.order == BondOrder::Triple)
        .unwrap()
        .0;
    assert_eq!(m.bond_is_aromatic(triple).unwrap(), Some(true));
    let result = enumerate_resonance(&m, Default::default()).unwrap();
    assert_eq!(result.contributors().len(), 1);
    let orders: std::collections::BTreeMap<_, _> = result.contributors()[0].bond_orders().collect();
    assert_eq!(orders[&triple], BondOrder::Triple);
    assert_eq!(
        orders.values().filter(|&&b| b == BondOrder::Double).count(),
        2
    );
    assert_eq!(
        result.contributors()[0]
            .formal_charges()
            .map(|(_, q)| q)
            .sum::<i8>(),
        0
    );
}

#[test]
fn explicit_preparation_preserves_chemistry_aromaticity_cip_and_output() {
    for model in [AromaticityModel::RdkitLike, AromaticityModel::Mdl] {
        let mut m = molecule("N[C@@H](C)c1ccccc1");
        perception::aromaticity::perceive_aromaticity(&mut m, model).unwrap();
        kekule::stereo::assign_cip_descriptors(&mut m).unwrap();
        let chemistry = m.clone();
        let before = m.perception().clone();
        let text = smiles::write_canonical(&m).unwrap();
        perceive_conjugation(&mut m, ConjugationModel::RdkitLike).unwrap();
        perceive_resonance(&mut m).unwrap();
        enumerate_resonance(&m, Default::default()).unwrap();
        assert_eq!(m, chemistry);
        assert_eq!(
            m.perception().aromaticity_state(),
            before.aromaticity_state()
        );
        assert_eq!(m.perception().stereo_state(), before.stereo_state());
        assert_eq!(smiles::write_canonical(&m).unwrap(), text);
        let state = m.perception().clone();
        m.clear_perception();
        m.install_perception(state.clone()).unwrap();
        assert_eq!(m.perception(), &state);
        perceive_resonance(&mut m).unwrap();
        assert_eq!(m.perception(), &state);
    }
}

#[test]
fn present_empty_sections_and_prerequisite_invalidation() {
    let mut m = molecule("CC");
    assert!(m.perception().has_conjugation());
    assert!(!m.perception().has_resonance());
    let a = m.atom_ids().next().unwrap();
    assert_eq!(m.perception().atom_is_conjugated(a), Some(false));
    assert_eq!(m.atom_is_conjugated(a).unwrap(), Some(false));
    assert!(m.atom_is_conjugated(AtomId::new(999)).is_err());
    let bond = m.bond_ids().next().unwrap();
    assert_eq!(m.bond_is_conjugated(bond).unwrap(), Some(false));
    assert!(m.bond_is_conjugated(BondId::new(999)).is_err());
    perceive_resonance(&mut m).unwrap();
    assert!(m
        .perception()
        .resonance_state()
        .unwrap()
        .groups()
        .is_empty());
    perception::valence::perceive_valence(&mut m, ValenceModel::RdkitLike).unwrap();
    assert!(!m.perception().has_conjugation());
    assert!(!m.perception().has_resonance());
    m.perceive().unwrap();
    perceive_resonance(&mut m).unwrap();
    let mut edit = m.edit();
    edit.atom_mut(a).unwrap().isotope = Some(13);
    let edited = edit.finish().unwrap();
    assert!(!edited.perception().has_conjugation());
    assert!(!edited.perception().has_resonance());
}

fn detached(m: &Molecule, groups: Vec<ResonanceGroup>) -> Perception {
    let c = m.perception().conjugation_state().unwrap();
    let a = m.perception().aromaticity_state().unwrap();
    Perception::builder()
        .with_aromaticity(a.model(), a.atoms().collect(), a.bonds().collect())
        .unwrap()
        .with_conjugation(c.model(), c.atoms().collect(), c.bonds().collect())
        .unwrap()
        .with_resonance_groups(groups)
        .unwrap()
        .build()
}

#[test]
fn detached_group_order_is_canonical_before_and_after_installation() {
    let mut m = molecule("CC(=O)NCC(=O)[O-]");
    perceive_resonance(&mut m).unwrap();
    let expected = m.perception().resonance_state().unwrap().clone();
    let mut groups = expected.groups().to_vec();
    assert_eq!(groups.len(), 2);
    groups.reverse();
    let state = detached(&m, groups);
    assert_eq!(state.resonance_state(), Some(&expected));
    m.install_perception(state).unwrap();
    let installed = m.perception().resonance_state().unwrap();
    assert_eq!(installed, &expected);
    for (index, group) in expected.groups().iter().enumerate() {
        for &atom in &group.atoms {
            assert_eq!(installed.atom_group(atom), Some(index));
        }
        for &bond in &group.bonds {
            assert_eq!(installed.bond_group(bond), Some(index));
        }
    }
}

#[test]
fn detached_group_validation_is_structural_and_atomic() {
    let mut m = molecule("CC(=O)NCC(=O)[O-]");
    perceive_resonance(&mut m).unwrap();
    let before = m.perception().clone();
    let groups = before.resonance_state().unwrap().groups().to_vec();
    let mut malformed = Vec::new();
    malformed.push(vec![]); // incomplete partition
    let mut empty = groups.clone();
    empty.push(ResonanceGroup {
        atoms: vec![],
        bonds: vec![],
    });
    malformed.push(empty);
    let mut duplicate = groups.clone();
    duplicate[0].atoms.push(groups[0].atoms[0]);
    malformed.push(duplicate);
    let mut duplicate = groups.clone();
    duplicate[0].bonds.push(groups[0].bonds[0]);
    malformed.push(duplicate);
    let mut ends = groups.clone();
    ends[0].atoms.pop();
    malformed.push(ends);
    let mut stale = groups.clone();
    stale[0].bonds[0] = BondId::new(999);
    malformed.push(stale);
    malformed.push(vec![ResonanceGroup {
        atoms: groups
            .iter()
            .flat_map(|g| g.atoms.iter().copied())
            .collect(),
        bonds: groups
            .iter()
            .flat_map(|g| g.bonds.iter().copied())
            .collect(),
    }]);
    for groups in malformed {
        assert!(m.install_perception(detached(&m, groups)).is_err());
        assert_eq!(m.perception(), &before);
    }
    m.install_perception(detached(&m, groups)).unwrap();
    assert!(m.perception().has_resonance());
}

#[test]
fn explicit_hydrogens_and_atom_renumbering_preserve_membership() {
    fn mapped_edges(m: &Molecule) -> std::collections::BTreeSet<[u32; 2]> {
        m.perception()
            .conjugation_state()
            .unwrap()
            .bonds()
            .map(|id| {
                let b = m.bond(id).unwrap();
                let mut ends = [
                    m.atom(b.a()).unwrap().atom_map.unwrap(),
                    m.atom(b.b()).unwrap().atom_map.unwrap(),
                ];
                ends.sort();
                ends
            })
            .collect()
    }
    let a = molecule("[CH3:1][C:2](=[O:3])[O-:4]");
    let b = molecule("[O-:4][C:2]([CH3:1])=[O:3]");
    let c = molecule("[O:4]=[C:2]([CH3:1])[O-:3]");
    assert_eq!(mapped_edges(&a), mapped_edges(&b));
    assert_eq!(mapped_edges(&a), mapped_edges(&c));
    for smi in ["CC(=O)[O-]", "CC(=O)N", "C=CC#N", "c1ccccc1"] {
        let mut explicit = molecule(smi);
        let original = explicit
            .perception()
            .conjugation_state()
            .unwrap()
            .bonds()
            .collect::<Vec<_>>();
        let count = enumerate_resonance(&explicit, Default::default())
            .unwrap()
            .contributors()
            .len();
        explicit.add_hydrogens().unwrap();
        explicit.perceive().unwrap();
        assert_eq!(
            explicit
                .perception()
                .conjugation_state()
                .unwrap()
                .bonds()
                .collect::<Vec<_>>(),
            original
        );
        assert_eq!(
            enumerate_resonance(&explicit, Default::default())
                .unwrap()
                .contributors()
                .len(),
            count
        );
    }
}

#[test]
fn saturated_and_unsupported_parts_do_not_hide_supported_groups() {
    for (smi, groups, bonds) in [
        ("C=CC(=O)C=C", 1, 5),
        ("C=CCCC=C", 0, 0),
        ("C=CC=C[SiH2]C=C", 1, 3),
        ("O=C(O)C", 1, 2),
        ("C=CC=O", 1, 3),
        ("C=CC#N", 1, 3),
        ("CC(=S)N", 1, 2),
        ("C=Cc1ccccc1", 1, 8),
    ] {
        let mut m = molecule(smi);
        perceive_resonance(&mut m).unwrap();
        assert_eq!(
            m.perception().resonance_state().unwrap().groups().len(),
            groups,
            "{smi}"
        );
        assert_eq!(
            m.perception().conjugation_state().unwrap().bonds().len(),
            bonds,
            "{smi}"
        );
    }
}

#[test]
fn tombstone_ids_survive_reconstruction_and_enumeration() {
    let mut editor = MoleculeEditor::new();
    let carbon = Atom::new(Element::from_symbol("C").unwrap());
    let dead = editor.add_atom(carbon.clone()).unwrap();
    editor.delete_atom(dead).unwrap();
    let a = editor.add_atom(carbon.clone()).unwrap();
    let b = editor.add_atom(carbon.clone()).unwrap();
    let c = editor.add_atom(carbon.clone()).unwrap();
    let d = editor.add_atom(carbon).unwrap();
    let removed = editor.add_bond(a, d, BondOrder::Single).unwrap();
    editor.delete_bond(removed).unwrap();
    editor.add_bond(a, b, BondOrder::Double).unwrap();
    editor.add_bond(b, c, BondOrder::Single).unwrap();
    editor.add_bond(c, d, BondOrder::Double).unwrap();
    let mut m = editor.finish().unwrap();
    m.perceive().unwrap();
    perceive_resonance(&mut m).unwrap();
    let state = m.perception().clone();
    m.clear_perception();
    m.install_perception(state.clone()).unwrap();
    assert_eq!(m.perception(), &state);
    let forms = enumerate_resonance(&m, Default::default()).unwrap();
    assert_eq!(
        forms.contributors()[0]
            .formal_charges()
            .map(|(id, _)| id)
            .collect::<Vec<_>>(),
        vec![a, b, c, d]
    );
    let invalid = Perception::builder()
        .with_aromaticity(AromaticityModel::RdkitLike, vec![], vec![])
        .unwrap()
        .with_conjugation(ConjugationModel::RdkitLike, vec![a, dead], vec![removed])
        .unwrap()
        .build();
    assert!(m.install_perception(invalid).is_err());
    assert_eq!(m.perception(), &state);
}

#[test]
fn structure_caps_and_unknown_hydrogens_are_explicit() {
    let m = molecule("CC(=O)[O-]");
    assert!(enumerate_resonance(
        &m,
        ResonanceOptions {
            max_structures: 0,
            max_total_work: 0,
            ..Default::default()
        }
    )
    .unwrap()
    .contributors()
    .is_empty());
    for (cap, count, reached) in [(0, 0, true), (2, 2, true), (3, 2, false)] {
        let forms = enumerate_resonance(
            &m,
            ResonanceOptions {
                max_structures: cap,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(forms.contributors().len(), count);
        assert_eq!(forms.limit_reached(), reached);
    }
    assert!(matches!(
        enumerate_resonance(
            &m,
            ResonanceOptions {
                max_structures: 1,
                ..Default::default()
            }
        ),
        Err(ResonanceError::InvalidStructureLimit)
    ));
    let mut m = smiles::to_molecules("CC").unwrap().remove(0);
    m.install_perception(
        Perception::builder()
            .with_aromaticity(AromaticityModel::RdkitLike, vec![], vec![])
            .unwrap()
            .build(),
    )
    .unwrap();
    let before = m.perception().clone();
    assert!(matches!(
        perceive_conjugation(&mut m, ConjugationModel::RdkitLike),
        Err(ConjugationError::UnknownHydrogens(_))
    ));
    assert_eq!(m.perception(), &before);
}
