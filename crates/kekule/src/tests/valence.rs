use super::*;

#[test]
fn valence_accepts_aromatic_input_localized_during_interpretation() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should interpret");
    assert_eq!(
        molecule
            .as_molecule()
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        3
    );
    assert_eq!(
        molecule
            .as_molecule()
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Single)
            .count(),
        3
    );

    valence_api::perceive_valence(molecule.as_molecule_mut(), ValenceModel::RdkitLike)
        .expect("localized benzene valence should succeed");
    assert!(molecule
        .as_molecule()
        .atom_ids()
        .all(|atom| molecule.as_molecule().implicit_hydrogens(atom) == Ok(Some(1))));
}

#[test]
fn localized_aromatic_valence_replaces_previous_valence_transactionally() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should interpret");
    let atom_ids = molecule.as_molecule().atom_ids().collect::<Vec<_>>();
    let bond_ids = molecule.as_molecule().bond_ids().collect::<Vec<_>>();
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
        .as_molecule_mut()
        .install_perception_state(previous.clone())
        .expect("valid previous perception");
    valence_api::perceive_valence_with_options(
        molecule.as_molecule_mut(),
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    )
    .expect("localized aromatic valence can be recomputed");

    assert!(molecule
        .as_molecule()
        .atom_ids()
        .all(|atom| molecule.as_molecule().implicit_hydrogens(atom) == Ok(Some(1))));
    assert!(molecule.as_molecule().perception().has_rings());
    assert!(!molecule.as_molecule().perception().has_aromaticity());
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
    assert!(
        !molecule.as_molecule().perception().has_aromaticity(),
        "{smiles}"
    );

    molecule
        .canonicalize_fixture()
        .unwrap_or_else(|error| panic!("aromatic fixture should normalize: {smiles}: {error}"));
    assert!(molecule
        .as_molecule()
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
    assert_eq!(
        molecule.as_molecule().perception(),
        &PerceptionState::default(),
        "{smiles}"
    );

    valence_api::perceive_valence(molecule.as_molecule_mut(), ValenceModel::RdkitLike)
        .unwrap_or_else(|error| panic!("localized valence should succeed: {smiles}: {error}"));
    assert!(
        molecule.as_molecule().perception().has_valence(),
        "{smiles}"
    );
    assert!(!molecule.as_molecule().perception().has_rings(), "{smiles}");
    assert!(
        !molecule.as_molecule().perception().has_aromaticity(),
        "{smiles}"
    );
    assert_eq!(
        molecule
            .as_molecule()
            .atom_ids()
            .map(|atom| {
                molecule
                    .as_molecule()
                    .implicit_hydrogens(atom)
                    .expect("live atom")
                    .expect("complete valence assignment")
            })
            .collect::<Vec<_>>(),
        expected_implicit_hydrogens,
        "{smiles}"
    );

    rings_api::perceive_ring_set(molecule.as_molecule_mut())
        .unwrap_or_else(|error| panic!("ring perception should succeed: {smiles}: {error}"));
    aromaticity_api::perceive_aromaticity(molecule.as_molecule_mut(), AromaticityModel::RdkitLike)
        .unwrap_or_else(|error| panic!("aromaticity perception should succeed: {smiles}: {error}"));
    assert_eq!(
        molecule
            .as_molecule()
            .atom_ids()
            .filter(|atom| molecule.as_molecule().atom_is_aromatic(*atom) == Ok(Some(true)))
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
            .as_molecule_mut()
            .atom_mut(AtomId::new(0))
            .expect("radical carbon");
        radical_carbon.radical = Some(AtomRadical::Doublet);
        radical_carbon.hydrogens = HydrogenDeclaration::Fixed(0);
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
    molecule
        .canonicalize_fixture()
        .expect("pyrrole should normalize");
    let represented = represented_molecule_snapshot(molecule.as_molecule());
    let represented_nitrogen = molecule
        .as_molecule()
        .atom(AtomId::new(0))
        .expect("pyrrole nitrogen")
        .clone();

    valence_api::perceive_valence(molecule.as_molecule_mut(), ValenceModel::RdkitLike)
        .expect("pyrrole valence should succeed without aromaticity");

    let nitrogen = molecule
        .as_molecule()
        .atom(AtomId::new(0))
        .expect("pyrrole nitrogen");
    assert_eq!(nitrogen.hydrogens, HydrogenDeclaration::Fixed(1));
    assert_eq!(nitrogen.hydrogens, represented_nitrogen.hydrogens);
    assert_eq!(
        molecule
            .as_molecule()
            .implicit_hydrogens(AtomId::new(0))
            .expect("live nitrogen"),
        Some(0)
    );
    assert!(!molecule.as_molecule().perception().has_aromaticity());
    assert_eq!(
        represented_molecule_snapshot(molecule.as_molecule()),
        represented
    );

    rings_api::perceive_ring_set(molecule.as_molecule_mut()).expect("pyrrole ring perception");
    aromaticity_api::perceive_aromaticity(molecule.as_molecule_mut(), AromaticityModel::RdkitLike)
        .expect("pyrrole aromaticity perception");

    assert_eq!(
        represented_molecule_snapshot(molecule.as_molecule()),
        represented
    );
    assert_eq!(
        molecule
            .as_molecule()
            .atom_ids()
            .filter(|atom| molecule.as_molecule().atom_is_aromatic(*atom) == Ok(Some(true)))
            .count(),
        5
    );
    let total_hydrogens = molecule
        .as_molecule()
        .atoms()
        .map(|(atom_id, atom)| {
            usize::from(atom.hydrogens.explicit_count())
                + usize::from(
                    molecule
                        .as_molecule()
                        .implicit_hydrogens(atom_id)
                        .expect("live atom")
                        .expect("complete valence assignment"),
                )
        })
        .sum::<usize>();
    assert_eq!(total_hydrogens, 5);

    let written = smiles_api::write(&molecule).expect("perceived pyrrole should write");
    assert!(written.contains("[nH]"), "{written}");
}

#[test]
fn valence_ignores_preinstalled_semantic_aromaticity() {
    let mut without_aromaticity = read_smiles("c1ccccc1").expect("benzene should parse");
    without_aromaticity
        .canonicalize_fixture()
        .expect("benzene should normalize");
    let mut with_aromaticity = without_aromaticity.clone();
    let aromatic_atoms = with_aromaticity
        .as_molecule()
        .atom_ids()
        .collect::<Vec<_>>();
    let aromatic_bonds = with_aromaticity
        .as_molecule()
        .bond_ids()
        .collect::<Vec<_>>();
    let previous = PerceptionState::builder()
        .with_aromaticity(AromaticityModel::RdkitLike, aromatic_atoms, aromatic_bonds)
        .expect("valid semantic aromaticity")
        .build();
    with_aromaticity
        .as_molecule_mut()
        .install_perception_state(previous)
        .expect("valid perception state");

    valence_api::perceive_valence(
        without_aromaticity.as_molecule_mut(),
        ValenceModel::RdkitLike,
    )
    .expect("valence without aromaticity");
    valence_api::perceive_valence(with_aromaticity.as_molecule_mut(), ValenceModel::RdkitLike)
        .expect("valence with preinstalled aromaticity");

    let without = without_aromaticity
        .as_molecule()
        .atom_ids()
        .map(|atom| without_aromaticity.as_molecule().implicit_hydrogens(atom))
        .collect::<Vec<_>>();
    let with = with_aromaticity
        .as_molecule()
        .atom_ids()
        .map(|atom| with_aromaticity.as_molecule().implicit_hydrogens(atom))
        .collect::<Vec<_>>();
    assert_eq!(with, without);
    assert_eq!(with, vec![Ok(Some(1)); 6]);
}

#[test]
fn fused_aromatic_valence_comes_from_localized_bond_orders() {
    let mut molecule = read_smiles("c1ccc2ccccc2c1").expect("naphthalene should parse");
    molecule
        .canonicalize_fixture()
        .expect("naphthalene should normalize");

    valence_api::perceive_valence(molecule.as_molecule_mut(), ValenceModel::RdkitLike)
        .expect("naphthalene valence should run first");

    let mut peripheral = 0;
    let mut fused = 0;
    for atom_id in molecule.as_molecule().atom_ids() {
        let degree = molecule
            .as_molecule()
            .incident_bonds(atom_id)
            .expect("live atom")
            .filter(|(_, bond)| !matches!(bond.order, BondOrder::Zero | BondOrder::Dative))
            .count();
        let implicit = molecule
            .as_molecule()
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

    rings_api::perceive_ring_set(molecule.as_molecule_mut()).expect("naphthalene rings");
    aromaticity_api::perceive_aromaticity(molecule.as_molecule_mut(), AromaticityModel::RdkitLike)
        .expect("naphthalene aromaticity");
    assert_eq!(
        molecule
            .as_molecule()
            .atom_ids()
            .filter(|atom| molecule.as_molecule().atom_is_aromatic(*atom) == Ok(Some(true)))
            .count(),
        10
    );
}
