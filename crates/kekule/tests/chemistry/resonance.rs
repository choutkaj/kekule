use kekule::{
    molecule::*,
    perception::{
        self, aromaticity::AromaticityModel, conjugation::*, resonance::*, valence::ValenceModel,
        Perception,
    },
    smiles,
};
use std::num::NonZeroU16;

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
            m.perception().conjugation().unwrap().bonds().len(),
            expected,
            "{smi}"
        );
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
            assert!(c.perception().conjugation().is_none());
        }
    }
    let m = molecule("c1ccccc1");
    let all = enumerate_resonance(
        &m,
        ResonanceOptions::default().with_flags(ResonanceFlags::KEKULE_ALL),
    )
    .unwrap();
    assert_eq!(all.contributors().len(), 2);
}

#[test]
fn resonance_perception_is_default_and_follows_aromaticity() {
    let mut m = molecule("CC(=O)NCC(=O)[O-]");
    let original = m.clone();
    let installed = m.perception().resonance().unwrap().clone();
    assert_eq!(installed.model(), ResonanceModel::ChargeMoves);
    assert_eq!(installed.groups().len(), 1);
    perceive_resonance(&mut m, ResonanceModel::ChargeMoves).unwrap();
    assert_eq!(m.perception().resonance(), Some(&installed));
    assert_eq!(m, original);
    let before = m.perception().clone();
    assert!(matches!(
        enumerate_resonance(&m, ResonanceOptions::default().with_max_total_work(0)),
        Err(ResonanceError::LimitExceeded(_))
    ));
    assert_eq!(m.perception(), &before);
    perception::aromaticity::perceive_aromaticity(
        &mut m,
        AromaticityModel::Mdl,
        Default::default(),
    )
    .unwrap();
    assert!(m.perception().conjugation().is_none());
    assert!(m.perception().resonance().is_none());
    perceive_conjugation(&mut m, ConjugationModel::RdkitLike).unwrap();
    assert!(m.perception().conjugation().is_some());
    assert!(m.perception().resonance().is_none());
    perceive_resonance(&mut m, ResonanceModel::ChargeMoves).unwrap();
    assert_eq!(m.perception().resonance(), Some(&installed));
    perceive_conjugation(&mut m, ConjugationModel::RdkitLike).unwrap();
    assert_eq!(m.perception().resonance(), Some(&installed));
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
        perceive_resonance(&mut m, ResonanceModel::ChargeMoves),
        Err(ResonancePerceptionError::MissingAromaticity)
    );
    assert!(matches!(
        enumerate_resonance(&m, Default::default()),
        Err(ResonanceError::MissingConjugation)
    ));
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
                ResonanceOptions::default()
                    .with_flags(ResonanceFlags::from_bits(mask as u8).unwrap()),
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
    assert_eq!(
        m.perception()
            .aromaticity()
            .map(|aromaticity| aromaticity.bond_is_aromatic(triple)),
        Some(true)
    );
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
        perception::aromaticity::perceive_aromaticity(&mut m, model, Default::default()).unwrap();
        kekule::stereo::cip::assign_cip_descriptors(&mut m, Default::default()).unwrap();
        let chemistry = m.clone();
        let before = m.perception().clone();
        let text = smiles::write(&m, smiles::SmilesWriteOptions::canonical()).unwrap();
        perceive_conjugation(&mut m, ConjugationModel::RdkitLike).unwrap();
        perceive_resonance(&mut m, ResonanceModel::ChargeMoves).unwrap();
        enumerate_resonance(&m, Default::default()).unwrap();
        assert_eq!(m, chemistry);
        assert_eq!(m.perception().aromaticity(), before.aromaticity());
        assert_eq!(m.perception().cip(), before.cip());
        assert_eq!(
            smiles::write(&m, smiles::SmilesWriteOptions::canonical()).unwrap(),
            text
        );
        let state = m.perception().clone();
        m.clear_perception();
        m.install_perception(state.clone()).unwrap();
        assert_eq!(m.perception(), &state);
        perceive_resonance(&mut m, ResonanceModel::ChargeMoves).unwrap();
        assert_eq!(m.perception(), &state);
    }
}

#[test]
fn present_empty_sections_and_prerequisite_invalidation() {
    let mut m = molecule("CC");
    assert!(m.perception().conjugation().is_some());
    assert!(m.perception().resonance().unwrap().groups().is_empty());
    let a = m.atom_ids().next().unwrap();
    let conjugation = m.perception().conjugation().unwrap();
    assert!(!conjugation.atom_is_conjugated(a));
    // The section answers membership; an ID outside the molecule is not a member.
    assert!(!conjugation.atom_is_conjugated(AtomId::new(999)));
    let bond = m.bond_ids().next().unwrap();
    assert!(!conjugation.bond_is_conjugated(bond));
    assert!(!conjugation.bond_is_conjugated(BondId::new(999)));
    let resonance = m.perception().resonance().unwrap();
    assert!(!resonance.bond_is_delocalized(bond));
    assert!(!resonance.bond_is_delocalized(BondId::new(999)));
    assert_eq!(resonance.atom_group(a), None);
    perception::valence::perceive_valence(&mut m, ValenceModel::RdkitLike, Default::default())
        .unwrap();
    assert!(m.perception().conjugation().is_none());
    assert!(m.perception().resonance().is_none());
    m.perceive().unwrap();
    let mut edit = m.to_editor();
    edit.atom_mut(a).unwrap().set_isotope(NonZeroU16::new(13));
    let edited = edit.finish().unwrap();
    assert!(edited.perception().conjugation().is_none());
    assert!(edited.perception().resonance().is_none());
}

fn detached(m: &Molecule, groups: Vec<ResonanceGroup>) -> Perception {
    let c = m.perception().conjugation().unwrap();
    let a = m.perception().aromaticity().unwrap();
    Perception::builder()
        .with_aromaticity(a.model(), a.atoms().collect(), a.bonds().collect())
        .unwrap()
        .with_conjugation(c.model(), c.atoms().collect(), c.bonds().collect())
        .unwrap()
        .with_resonance(ResonanceModel::ChargeMoves, groups)
        .unwrap()
        .build()
}

#[test]
fn detached_group_order_is_canonical_before_and_after_installation() {
    let mut m = molecule("[O-]C(=O)CCC(=O)[O-]");
    let expected = m.perception().resonance().unwrap().clone();
    let mut groups = expected.groups().to_vec();
    assert_eq!(groups.len(), 2);
    groups.reverse();
    let state = detached(&m, groups);
    assert_eq!(state.resonance(), Some(&expected));
    m.install_perception(state).unwrap();
    let installed = m.perception().resonance().unwrap();
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
    let mut m = molecule("[O-]C(=O)CCC(=O)[O-]");
    let before = m.perception().clone();
    let groups = before.resonance().unwrap().groups().to_vec();
    let mut malformed = Vec::new();
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
    let without_aromaticity = Perception::builder()
        .with_resonance(ResonanceModel::ChargeMoves, groups.clone())
        .unwrap()
        .build();
    assert!(m.install_perception(without_aromaticity).is_err());
    m.install_perception(detached(&m, vec![])).unwrap();
    assert!(m.perception().resonance().unwrap().groups().is_empty());
    m.install_perception(detached(&m, groups)).unwrap();
    assert_eq!(m.perception().resonance(), before.resonance());
    let mut phenolate = molecule("[O-]c1ccccc1");
    let ring = phenolate
        .perception()
        .aromaticity()
        .unwrap()
        .bonds()
        .next()
        .unwrap();
    let bond = phenolate.bond(ring).unwrap();
    let aromatic = vec![ResonanceGroup {
        atoms: vec![bond.a().min(bond.b()), bond.a().max(bond.b())],
        bonds: vec![ring],
    }];
    let state = detached(&phenolate, aromatic);
    assert!(phenolate.install_perception(state).is_err());
}

#[test]
fn explicit_hydrogens_and_atom_renumbering_preserve_membership() {
    fn mapped_edges(m: &Molecule) -> std::collections::BTreeSet<[u32; 2]> {
        m.perception()
            .conjugation()
            .unwrap()
            .bonds()
            .map(|id| {
                let b = m.bond(id).unwrap();
                let mut ends = [
                    m.atom(b.a()).unwrap().atom_map().unwrap().get(),
                    m.atom(b.b()).unwrap().atom_map().unwrap().get(),
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
            .conjugation()
            .unwrap()
            .bonds()
            .collect::<Vec<_>>();
        let count = enumerate_resonance(&explicit, Default::default())
            .unwrap()
            .contributors()
            .len();
        explicit.add_hydrogens(Default::default()).unwrap();
        explicit.perceive().unwrap();
        assert_eq!(
            explicit
                .perception()
                .conjugation()
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
        let m = molecule(smi);
        assert_eq!(
            enumerate_resonance(&m, Default::default())
                .unwrap()
                .groups()
                .len(),
            groups,
            "{smi}"
        );
        assert_eq!(
            m.perception().conjugation().unwrap().bonds().len(),
            bonds,
            "{smi}"
        );
    }
}

#[test]
fn edited_ids_compact_before_reconstruction_and_enumeration() {
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
    let (mut m, ids) = editor.finish_with_correspondence().unwrap();
    // Deleted draft slots do not survive publication.
    assert_eq!((ids.atom(dead), ids.bond(removed)), (None, None));
    let [a, b, c, d] = [a, b, c, d].map(|atom| ids.atom(atom).unwrap());
    assert_eq!(
        [a, b, c, d],
        [0, 1, 2, 3].map(kekule::molecule::AtomId::new)
    );
    m.perceive().unwrap();
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
        .with_conjugation(
            ConjugationModel::RdkitLike,
            vec![a, kekule::molecule::AtomId::new(4)],
            vec![kekule::molecule::BondId::new(3)],
        )
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
        ResonanceOptions::default()
            .with_max_structures(0)
            .with_max_total_work(0)
    )
    .unwrap()
    .contributors()
    .is_empty());
    for (cap, count, reached) in [(0, 0, true), (2, 2, true), (3, 2, false)] {
        let forms =
            enumerate_resonance(&m, ResonanceOptions::default().with_max_structures(cap)).unwrap();
        assert_eq!(forms.contributors().len(), count);
        assert_eq!(forms.limit_reached(), reached);
    }
    assert!(matches!(
        enumerate_resonance(&m, ResonanceOptions::default().with_max_structures(1)),
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
        Err(ConjugationError::UnknownHydrogenCount(_))
    ));
    assert_eq!(m.perception(), &before);
    // Neutral molecules need no hydrogen counts; a charge's path does.
    perceive_resonance(&mut m, ResonanceModel::ChargeMoves).unwrap();
    assert!(m.perception().resonance().unwrap().groups().is_empty());
    let mut m = smiles::to_molecules("CC(=O)[O-]").unwrap().remove(0);
    m.install_perception(
        Perception::builder()
            .with_aromaticity(AromaticityModel::RdkitLike, vec![], vec![])
            .unwrap()
            .build(),
    )
    .unwrap();
    let before = m.perception().clone();
    assert!(matches!(
        perceive_resonance(&mut m, ResonanceModel::ChargeMoves),
        Err(ResonancePerceptionError::UnknownHydrogenCount(_))
    ));
    assert_eq!(m.perception(), &before);
}

/// Delocalized bonds as sorted pairs of atom-map numbers.
fn delocalized(smi: &str) -> Vec<[u32; 2]> {
    let m = molecule(smi);
    let resonance = m.perception().resonance().unwrap();
    let map = |atom| m.atom(atom).unwrap().atom_map().map_or(0, |n| n.get());
    let mut pairs: Vec<_> = resonance
        .bonds()
        .map(|id| {
            let bond = m.bond(id).unwrap();
            let mut pair = [map(bond.a()), map(bond.b())];
            pair.sort_unstable();
            pair
        })
        .collect();
    pairs.sort_unstable();
    pairs
}

#[test]
fn moving_charges_delocalize_their_paths() {
    for (smi, groups, bonds) in [
        ("CC(=O)[O-]", 1, 2),
        ("C=CC(=O)[O-]", 1, 2),
        ("[O-]C(=O)C(=O)[O-]", 2, 4),
        ("NC(N)=[NH2+]", 1, 3),
        ("CCCNC(N)=[NH2+]", 1, 3),
        ("CC(N)=[NH2+]", 1, 2),
        ("C[N+](=O)[O-]", 1, 2),
        ("[O-][N+](=O)[O-]", 1, 3),
        ("COP(=O)([O-])[O-]", 1, 3),
        ("COP(=O)([O-])OC", 1, 2),
        ("COP(=S)([O-])[O-]", 1, 3),
        ("CS(=O)(=O)[O-]", 1, 3),
        ("CC(=O)[CH-]C(C)=O", 1, 4),
        ("CC(=O)[N-]C", 1, 2),
        ("CN=[N+]=[N-]", 1, 2),
        ("C=C[O-]", 1, 2),
        ("[CH2+]C=C", 1, 2),
        ("[CH2-]C=C", 1, 2),
        ("[O-]c1ccc(cc1)[N+](=O)[O-]", 1, 2),
    ] {
        let m = molecule(smi);
        let resonance = m.perception().resonance().unwrap();
        assert_eq!(resonance.groups().len(), groups, "{smi}");
        assert_eq!(resonance.bonds().len(), bonds, "{smi}");
        for (index, group) in resonance.groups().iter().enumerate() {
            for &bond in &group.bonds {
                assert!(resonance.bond_is_delocalized(bond));
                assert_eq!(resonance.bond_group(bond), Some(index));
            }
            for &atom in &group.atoms {
                assert_eq!(resonance.atom_group(atom), Some(index));
            }
        }
    }
}

#[test]
fn neutral_push_pull_aromatic_and_trapped_charges_stay_localized() {
    for smi in [
        "CC(=O)N",
        "CC([O-])=[NH2+]",
        "CC(=O)O",
        "COP(=O)(O)O",
        "CS(C)=O",
        "C[S+](C)[O-]",
        "C[N+](C)(C)[O-]",
        "C[NH3+]",
        "C=[N+](C)C",
        "C=CC=C",
        "C=CC=O",
        "C=Cc1ccccc1",
        "[O-]c1ccccc1",
        "Cc1c[nH+]c[nH]1",
        "c1ccccc1",
    ] {
        let m = molecule(smi);
        assert!(
            m.perception().resonance().unwrap().groups().is_empty(),
            "{smi}"
        );
    }
}

#[test]
fn ylidic_and_double_bonded_drawings_delocalize_alike() {
    for (ylidic, double) in [
        (
            "[CH3:1][O:2][P+:3]([O-:4])([O-:5])[O-:6]",
            "[CH3:1][O:2][P:3](=[O:4])([O-:5])[O-:6]",
        ),
        (
            "[CH3:1][S++:2]([O-:3])([O-:4])[O-:5]",
            "[CH3:1][S:2](=[O:3])(=[O:4])[O-:5]",
        ),
        (
            "[CH3:1][O:2][P+:3]([O-:4])([O-:5])[O:6][CH3:7]",
            "[CH3:1][O:2][P:3](=[O:4])([O-:5])[O:6][CH3:7]",
        ),
    ] {
        assert_eq!(delocalized(ylidic), delocalized(double), "{ylidic}");
        assert!(!delocalized(double).is_empty(), "{double}");
    }
    assert_eq!(
        delocalized("[CH3:1][N+:2](=[O:3])[O-:4]"),
        vec![[2, 3], [2, 4]]
    );
    assert_eq!(
        delocalized("[O-:1][C:2](=[O:3])[CH2:4][C:5](=[O:6])[NH2:7]"),
        vec![[1, 2], [2, 3]]
    );
}
