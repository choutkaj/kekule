use super::*;

#[test]
fn huckel_electron_count_does_not_narrow_large_candidate_systems() {
    let donors = vec![AromaticElectronDonorType::One; 130];

    assert_eq!(huckel_electron_count_for_donors(&donors), Some(130));
}

#[test]
fn fused_ten_electron_system_marks_its_explicit_single_fusion_bond_aromatic() {
    let source = "On2c1-c(ccc2)ccn1";
    let document = crate::smiles::parse_str(source, Default::default()).expect("parses");
    let (mut molecule, report) = document
        .interpret(Default::default())
        .expect("interprets")
        .into_parts()
        .expect("one component");
    let explicit_single_offset = source.find('-').expect("explicit single marker") + 1;
    let fusion_single = report
        .bond_mappings()
        .iter()
        .find(|mapping| mapping.source_offset() == explicit_single_offset)
        .map(|mapping| mapping.bond())
        .expect("explicit fusion single mapping");
    assert_eq!(
        molecule.bond(fusion_single).unwrap().order,
        BondOrder::Single
    );
    let valence = perceive_valence(&mut molecule, ValenceModel::RdkitLike, Default::default());
    assert!(valence.is_ok(), "{valence:#?}");
    perceive_rings(&mut molecule, RingPerceptionOptions::default()).expect("rings");

    perceive_aromaticity(
        &mut molecule,
        AromaticityModel::RdkitLike,
        Default::default(),
    )
    .expect("fused aromaticity");

    let fusion = molecule.bond(fusion_single).expect("fusion bond");
    // Represented chemistry keeps the single bond; the rings are aromatic
    // only together, and the fusion bond is aromatic with them.
    assert_eq!(fusion.order, BondOrder::Single);
    assert_eq!(
        molecule
            .perception()
            .aromaticity()
            .map(|aromaticity| aromaticity.bond_is_aromatic(fusion_single)),
        Some(true)
    );
    assert_eq!(
        molecule
            .atoms()
            .filter(|(atom_id, atom)| {
                atom.element.symbol() != "O"
                    && molecule
                        .perception()
                        .aromaticity()
                        .map(|aromaticity| aromaticity.atom_is_aromatic(*atom_id))
                        == Some(true)
            })
            .count(),
        9
    );
}

#[test]
fn canonical_localized_dye_assigns_nitrogen_hydrogens_before_aromaticity() {
    let input = "N2c1c(Nc3c2c6c(OS(=O)(=O)[O-])c7c(cccc7)c(OS(=O)(=O)[O-])c6cc3Cl)c4c(OS(=O)(=O)[O-])c5c(cccc5)c(OS(=O)(=O)[O-])c4cc1Cl";
    let mut molecule = crate::test_fixtures::read_smiles(input).expect("dye parses");
    assert!(molecule.perception().aromaticity().is_none());
    let valence = perceive_valence(&mut molecule, ValenceModel::RdkitLike, Default::default());
    assert!(valence.is_ok(), "{valence:#?}");
    let valence_nitrogens = molecule
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "N")
        .map(|(atom_id, _)| molecule.implicit_hydrogens(atom_id))
        .collect::<Vec<_>>();
    assert_eq!(valence_nitrogens, vec![Some(1), Some(1)]);
    assert!(molecule.perception().aromaticity().is_none());
    perceive_aromaticity(
        &mut molecule,
        AromaticityModel::RdkitLike,
        Default::default(),
    )
    .expect("aromaticity");

    let nitrogens = molecule
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "N")
        .map(|(atom_id, _)| {
            (
                molecule
                    .perception()
                    .aromaticity()
                    .map(|aromaticity| aromaticity.atom_is_aromatic(atom_id)),
                molecule.implicit_hydrogens(atom_id),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        nitrogens,
        vec![(Some(false), Some(1)), (Some(false), Some(1))]
    );
}

#[test]
fn neutral_carbon_radical_can_complete_an_aromatic_sextet() {
    let mut molecule =
        crate::test_fixtures::read_smiles("C1=CC(=CC=[C]1)N").expect("aminophenyl radical parses");
    crate::test_fixtures::edit(&mut molecule, |editor| {
        editor
            .atom_mut(AtomId::new(5))
            .expect("source-selected radical carbon")
            .radical = AtomRadical::new(1, Some(2));
    });
    let valence = perceive_valence(&mut molecule, ValenceModel::RdkitLike, Default::default());
    assert!(valence.is_ok(), "{valence:#?}");
    perceive_rings(&mut molecule, RingPerceptionOptions::default()).expect("ring perception");
    let radical = molecule
        .atoms()
        .find_map(|(id, atom)| atom.radical.is_some().then_some(id))
        .expect("radical carbon");
    let radical_atom = molecule.atom(radical).expect("radical atom");
    let donor = rdkit_localized_atom_donor_type(&molecule, radical, radical_atom);
    assert_eq!(donor, AromaticElectronDonorType::One);
    assert!(atom_is_rdkit_aromatic_candidate_for_donor(
        &molecule,
        radical,
        radical_atom,
        donor,
    ));

    perceive_aromaticity(
        &mut molecule,
        AromaticityModel::RdkitLike,
        Default::default(),
    )
    .expect("aromaticity perception");

    assert_eq!(
        molecule
            .atoms()
            .filter(|(atom_id, atom)| {
                atom.element.symbol() == "C"
                    && molecule
                        .perception()
                        .aromaticity()
                        .map(|aromaticity| aromaticity.atom_is_aromatic(*atom_id))
                        == Some(true)
            })
            .count(),
        6
    );
}

#[test]
fn charge_adjusted_candidate_valence_does_not_change_carbocation_electron_count() {
    for smiles in ["C1=C[C+]=CC(=C1)N", "C1=C[C+]=CC(=C1)C=O"] {
        let mut molecule = crate::test_fixtures::read_smiles(smiles).expect("carbocation parses");
        molecule.perceive().expect("carbocation perceives");
        assert!(
            molecule.atom_ids().all(|atom_id| {
                molecule
                    .perception()
                    .aromaticity()
                    .map(|aromaticity| aromaticity.atom_is_aromatic(atom_id))
                    == Some(false)
            }),
            "{smiles}"
        );
    }
}

#[test]
fn rdkit_candidate_elements_include_light_main_group_and_heavy_chalcogens() {
    // Full atom/bond masks checked against RDKit 2026.03.3 MolFromSmiles.
    for (source, aromatic) in [
        ("[Be-]1=CC=CC=C1", true),
        ("[Mg-]1=CC=CC=C1", true),
        ("[Al]1=CC=CC=C1", true),
        ("[SiH]1=CC=CC=C1", true),
        ("[SiH-]1C=CC=C1", true),
        ("[PH]1C=CC=C1", true),
        ("[Se]1C=CC=C1", true),
        ("[Te]1C=CC=C1", true),
        ("[GeH]1=CC=CC=C1", false),
        ("[AsH]1C=CC=C1", false),
    ] {
        let mut molecule = crate::test_fixtures::read_smiles(source).expect(source);
        molecule.perceive().expect(source);
        assert!(
            molecule.atom_ids().all(|atom| {
                molecule
                    .perception()
                    .aromaticity()
                    .map(|aromaticity| aromaticity.atom_is_aromatic(atom))
                    == Some(aromatic)
            }),
            "atom mask for {source}"
        );
        assert!(
            molecule.bond_ids().all(|bond| {
                molecule
                    .perception()
                    .aromaticity()
                    .map(|aromaticity| aromaticity.bond_is_aromatic(bond))
                    == Some(aromatic)
            }),
            "bond mask for {source}"
        );
    }
}

#[test]
fn exocyclic_electron_withdrawal_uses_all_neighbor_elements() {
    // Exocyclic halogens and transition metals participate in RDKit's
    // outer-electron electronegativity comparison despite not being ring
    // candidates. An equally or less electronegative neighbor does not
    // remove the carbon's electron.
    for (source, aromatic) in [
        ("C1(=[Cl+])C=CC=CC=C1", true),
        ("C1(=[I+])C=CC=CC=C1", true),
        ("C1(=[Fe])C=CC=CC=C1", true),
        ("C1(=[Se])C=CC=CC=C1", true),
        ("C1(=O)C=CC=CC=C1", true),
        ("C1(=C)C=CC=CC=C1", false),
        ("C1(=[SiH2])C=CC=CC=C1", false),
    ] {
        let mut molecule = crate::test_fixtures::read_smiles(source).expect(source);
        molecule.perceive().expect(source);
        for atom in molecule.atom_ids() {
            assert_eq!(
                molecule
                    .perception()
                    .aromaticity()
                    .map(|aromaticity| aromaticity.atom_is_aromatic(atom)),
                Some(aromatic && atom.index() != 1),
                "{source}: {atom}"
            );
        }
        for (bond_id, bond) in molecule.bonds() {
            let (left, right) = bond.endpoints();
            assert_eq!(
                molecule
                    .perception()
                    .aromaticity()
                    .map(|aromaticity| aromaticity.bond_is_aromatic(bond_id)),
                Some(aromatic && left.index() != 1 && right.index() != 1),
                "{source}: {bond_id}"
            );
        }
    }
}

#[test]
fn zero_electron_carbon_with_exocyclic_multiple_bond_is_vacant() {
    let mut molecule = crate::test_fixtures::read_smiles("C1(=O)C=CC=CC=C1").expect("tropone");
    crate::test_fixtures::edit(&mut molecule, |editor| {
        editor.atom_mut(AtomId::new(0)).unwrap().radical = AtomRadical::new(1, Some(2));
    });
    perceive_rings(&mut molecule, RingPerceptionOptions::default()).expect("rings");
    let atom = molecule.atom(AtomId::new(0)).unwrap();
    assert_eq!(
        count_rdkit_like_atom_pi_electrons(&molecule, AtomId::new(0), atom),
        Some(0)
    );
    assert_eq!(
        rdkit_localized_atom_donor_type(&molecule, AtomId::new(0), atom),
        AromaticElectronDonorType::Vacant
    );

    // Matches Kekulize(clearAromaticFlags=True), setting atom zero's
    // radical count to one, then SetAromaticity in RDKit 2026.03.3.
    perceive_aromaticity(
        &mut molecule,
        AromaticityModel::RdkitLike,
        Default::default(),
    )
    .expect("aromaticity");
    for atom in molecule.atom_ids() {
        assert_eq!(
            molecule
                .perception()
                .aromaticity()
                .map(|aromaticity| aromaticity.atom_is_aromatic(atom)),
            Some(atom.index() != 1)
        );
    }
    for (bond_id, bond) in molecule.bonds() {
        let (left, right) = bond.endpoints();
        assert_eq!(
            molecule
                .perception()
                .aromaticity()
                .map(|aromaticity| aromaticity.bond_is_aromatic(bond_id)),
            Some(left.index() != 1 && right.index() != 1)
        );
    }
}

#[test]
fn zero_bond_contributes_to_rdkit_unsaturation_degree() {
    let mut molecule = crate::test_fixtures::read_smiles("C1=C=CC=C1").expect("cyclic allene");
    let (helium, zero) = crate::test_fixtures::edit(&mut molecule, |editor| {
        let helium = editor
            .add_atom(Atom::new(Element::from_symbol("He").unwrap()))
            .unwrap();
        let zero = editor
            .add_bond_unpruned(AtomId::new(1), helium, BondOrder::Zero)
            .unwrap();
        (helium, zero)
    });
    molecule.perceive().expect("zero bond model");

    // RDKit's raw graph degree includes the zero-order bond. The resulting
    // explicit-valence-minus-degree is one, so the two double bonds do not
    // trigger the multiple-unsaturation exclusion and donate two electrons.
    for atom in molecule.atom_ids() {
        assert_eq!(
            molecule
                .perception()
                .aromaticity()
                .map(|aromaticity| aromaticity.atom_is_aromatic(atom)),
            Some(atom != helium)
        );
    }
    for bond in molecule.bond_ids() {
        assert_eq!(
            molecule
                .perception()
                .aromaticity()
                .map(|aromaticity| aromaticity.bond_is_aromatic(bond)),
            Some(bond != zero)
        );
    }
}

mod fused;
