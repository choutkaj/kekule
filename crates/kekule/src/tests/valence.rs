use super::*;

#[test]
fn valence_rejects_all_unnormalized_aromatic_bonds_before_assigning_hydrogens() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should interpret");
    let expected_issues = molecule
        .graph()
        .bond_ids()
        .map(ValenceIssue::UnsupportedBondOrder)
        .collect::<Vec<_>>();

    let error = valence_api::perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike)
        .expect_err("unnormalized aromatic bonds must be rejected");

    assert_eq!(error.issues, expected_issues);
    assert_eq!(molecule.graph().perception(), &PerceptionState::default());
    assert!(molecule
        .graph()
        .atom_ids()
        .all(|atom| molecule.graph().implicit_hydrogens(atom) == Ok(None)));

    molecule.normalize().expect("benzene should normalize");
    valence_api::perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike)
        .expect("normalized benzene valence should succeed");
    assert!(molecule
        .graph()
        .atom_ids()
        .all(|atom| molecule.graph().implicit_hydrogens(atom) == Ok(Some(1))));
}

#[test]
fn unnormalized_aromatic_valence_failure_preserves_complete_previous_perception() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should interpret");
    let atom_ids = molecule.graph().atom_ids().collect::<Vec<_>>();
    let bond_ids = molecule.graph().bond_ids().collect::<Vec<_>>();
    let expected_issues = bond_ids
        .iter()
        .copied()
        .map(ValenceIssue::UnsupportedBondOrder)
        .collect::<Vec<_>>();
    let previous = PerceptionState::builder()
        .with_valence(
            Some(ValenceModel::RdkitLike),
            atom_ids.iter().copied().map(|atom| (atom, 2)).collect(),
        )
        .expect("valid previous valence")
        .with_rings(
            RingMembership::from_slot_flags(vec![true; atom_ids.len()], vec![true; bond_ids.len()]),
            None,
        )
        .with_aromaticity(AromaticityModel::RdkitLike, atom_ids.clone(), bond_ids)
        .expect("valid previous aromaticity")
        .build();
    molecule
        .graph_mut()
        .install_perception_state(previous.clone())
        .expect("valid previous perception");
    let before = molecule.clone();

    let error = valence_api::perceive_valence_with_options(
        molecule.graph_mut(),
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    )
    .expect_err("normalization preflight cannot be made permissive");

    assert_eq!(error.issues, expected_issues);
    assert_eq!(molecule.graph().perception(), &previous);
    assert_eq!(molecule, before);
}

fn assert_aromatic_valence_pipeline(
    smiles: &str,
    expected_implicit_hydrogens: &[u8],
    expected_aromatic_atoms: usize,
) {
    let molecule = read_smiles(smiles)
        .unwrap_or_else(|error| panic!("aromatic fixture should parse: {smiles}: {error}"));
    assert_aromatic_valence_pipeline_for_molecule(
        smiles,
        molecule,
        expected_implicit_hydrogens,
        expected_aromatic_atoms,
    );
}

fn assert_aromatic_valence_pipeline_for_molecule(
    label: &str,
    mut molecule: SmallMolecule,
    expected_implicit_hydrogens: &[u8],
    expected_aromatic_atoms: usize,
) {
    let smiles = label;
    assert!(!molecule.graph().perception().has_aromaticity(), "{smiles}");

    molecule
        .normalize()
        .unwrap_or_else(|error| panic!("aromatic fixture should normalize: {smiles}: {error}"));
    assert!(molecule
        .graph()
        .bonds()
        .all(|(_, bond)| bond.order != BondOrder::Aromatic));
    assert_eq!(
        molecule.graph().perception(),
        &PerceptionState::default(),
        "{smiles}"
    );

    valence_api::perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike)
        .unwrap_or_else(|error| panic!("localized valence should succeed: {smiles}: {error}"));
    assert!(molecule.graph().perception().has_valence(), "{smiles}");
    assert!(!molecule.graph().perception().has_rings(), "{smiles}");
    assert!(!molecule.graph().perception().has_aromaticity(), "{smiles}");
    assert_eq!(
        molecule
            .graph()
            .atom_ids()
            .map(|atom| {
                molecule
                    .graph()
                    .implicit_hydrogens(atom)
                    .expect("live atom")
                    .expect("complete valence assignment")
            })
            .collect::<Vec<_>>(),
        expected_implicit_hydrogens,
        "{smiles}"
    );

    rings_api::perceive_ring_set(molecule.graph_mut())
        .unwrap_or_else(|error| panic!("ring perception should succeed: {smiles}: {error}"));
    aromaticity_api::perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
        .unwrap_or_else(|error| panic!("aromaticity perception should succeed: {smiles}: {error}"));
    assert_eq!(
        molecule
            .graph()
            .atom_ids()
            .filter(|atom| molecule.graph().atom_is_aromatic(*atom) == Ok(Some(true)))
            .count(),
        expected_aromatic_atoms,
        "{smiles}"
    );
}

#[test]
fn normalized_aromatic_systems_perceive_valence_before_rings_and_aromaticity() {
    for (smiles, implicit_hydrogens, aromatic_atoms) in [
        ("c1ccccc1", &[1, 1, 1, 1, 1, 1][..], 6),
        ("n1ccccc1", &[0, 1, 1, 1, 1, 1][..], 6),
        ("[nH]1cccc1", &[0, 1, 1, 1, 1][..], 5),
        ("c1ccoc1", &[1, 1, 1, 0, 1][..], 5),
        ("c1ccsc1", &[1, 1, 1, 0, 1][..], 5),
        ("C[n+]1ccccc1", &[3, 0, 1, 1, 1, 1, 1][..], 6),
        ("c1[n-]cnn1", &[1, 0, 1, 0, 0][..], 5),
    ] {
        assert_aromatic_valence_pipeline(smiles, implicit_hydrogens, aromatic_atoms);
    }

    let mut radical = read_smiles("c1ccccc1").expect("radical fixture syntax should interpret");
    {
        let mut radical_carbon = radical
            .graph_mut()
            .atom_mut(AtomId::new(0))
            .expect("radical carbon");
        radical_carbon.radical = Some(AtomRadical::Doublet);
        radical_carbon.no_implicit_hydrogens = true;
    }
    assert_aromatic_valence_pipeline_for_molecule(
        "explicitly represented phenyl radical",
        radical,
        &[0, 1, 1, 1, 1, 1],
        6,
    );
}

#[test]
fn normalized_pyrrole_retains_represented_hydrogen_before_valence() {
    let mut molecule = read_smiles("[nH]1cccc1").expect("pyrrole should parse");
    molecule.normalize().expect("pyrrole should normalize");
    let represented = represented_molecule_snapshot(molecule.graph());
    let represented_nitrogen = molecule
        .graph()
        .atom(AtomId::new(0))
        .expect("pyrrole nitrogen")
        .clone();

    valence_api::perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike)
        .expect("pyrrole valence should succeed without aromaticity");

    let nitrogen = molecule
        .graph()
        .atom(AtomId::new(0))
        .expect("pyrrole nitrogen");
    assert_eq!(nitrogen.explicit_hydrogens, 1);
    assert_eq!(
        nitrogen.explicit_hydrogens,
        represented_nitrogen.explicit_hydrogens
    );
    assert_eq!(
        nitrogen.no_implicit_hydrogens,
        represented_nitrogen.no_implicit_hydrogens
    );
    assert_eq!(
        molecule
            .graph()
            .implicit_hydrogens(AtomId::new(0))
            .expect("live nitrogen"),
        Some(0)
    );
    assert!(!molecule.graph().perception().has_aromaticity());
    assert_eq!(represented_molecule_snapshot(molecule.graph()), represented);

    rings_api::perceive_ring_set(molecule.graph_mut()).expect("pyrrole ring perception");
    aromaticity_api::perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
        .expect("pyrrole aromaticity perception");

    assert_eq!(represented_molecule_snapshot(molecule.graph()), represented);
    assert_eq!(
        molecule
            .graph()
            .atom_ids()
            .filter(|atom| molecule.graph().atom_is_aromatic(*atom) == Ok(Some(true)))
            .count(),
        5
    );
    let total_hydrogens = molecule
        .graph()
        .atoms()
        .map(|(atom_id, atom)| {
            usize::from(atom.explicit_hydrogens)
                + usize::from(
                    molecule
                        .graph()
                        .implicit_hydrogens(atom_id)
                        .expect("live atom")
                        .expect("complete valence assignment"),
                )
        })
        .sum::<usize>();
    assert_eq!(total_hydrogens, 5);

    let written = smiles_api::write_with_options(&molecule, SmilesWriteOptions)
        .expect("perceived pyrrole should write");
    assert!(written.contains("[nH]"), "{written}");
}

#[test]
fn valence_ignores_preinstalled_semantic_aromaticity() {
    let mut without_aromaticity = read_smiles("c1ccccc1").expect("benzene should parse");
    without_aromaticity
        .normalize()
        .expect("benzene should normalize");
    let mut with_aromaticity = without_aromaticity.clone();
    let aromatic_atoms = with_aromaticity.graph().atom_ids().collect::<Vec<_>>();
    let aromatic_bonds = with_aromaticity.graph().bond_ids().collect::<Vec<_>>();
    let previous = PerceptionState::builder()
        .with_aromaticity(AromaticityModel::RdkitLike, aromatic_atoms, aromatic_bonds)
        .expect("valid semantic aromaticity")
        .build();
    with_aromaticity
        .graph_mut()
        .install_perception_state(previous)
        .expect("valid perception state");

    valence_api::perceive_valence(without_aromaticity.graph_mut(), ValenceModel::RdkitLike)
        .expect("valence without aromaticity");
    valence_api::perceive_valence(with_aromaticity.graph_mut(), ValenceModel::RdkitLike)
        .expect("valence with preinstalled aromaticity");

    let without = without_aromaticity
        .graph()
        .atom_ids()
        .map(|atom| without_aromaticity.graph().implicit_hydrogens(atom))
        .collect::<Vec<_>>();
    let with = with_aromaticity
        .graph()
        .atom_ids()
        .map(|atom| with_aromaticity.graph().implicit_hydrogens(atom))
        .collect::<Vec<_>>();
    assert_eq!(with, without);
    assert_eq!(with, vec![Ok(Some(1)); 6]);
}

#[test]
fn fused_aromatic_valence_comes_from_localized_bond_orders() {
    let mut molecule = read_smiles("c1ccc2ccccc2c1").expect("naphthalene should parse");
    molecule.normalize().expect("naphthalene should normalize");

    valence_api::perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike)
        .expect("naphthalene valence should run first");

    let mut peripheral = 0;
    let mut fused = 0;
    for atom_id in molecule.graph().atom_ids() {
        let degree = molecule
            .graph()
            .incident_bonds(atom_id)
            .expect("live atom")
            .filter(|(_, bond)| !matches!(bond.order, BondOrder::Zero | BondOrder::Dative))
            .count();
        let implicit = molecule
            .graph()
            .implicit_hydrogens(atom_id)
            .expect("live atom")
            .expect("complete valence assignment");
        match degree {
            2 => {
                peripheral += 1;
                assert_eq!(implicit, 1);
            }
            3 => {
                fused += 1;
                assert_eq!(implicit, 0);
            }
            _ => panic!("unexpected naphthalene atom degree {degree}"),
        }
    }
    assert_eq!((peripheral, fused), (8, 2));

    rings_api::perceive_ring_set(molecule.graph_mut()).expect("naphthalene rings");
    aromaticity_api::perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
        .expect("naphthalene aromaticity");
    assert_eq!(
        molecule
            .graph()
            .atom_ids()
            .filter(|atom| molecule.graph().atom_is_aromatic(*atom) == Ok(Some(true)))
            .count(),
        10
    );
}
