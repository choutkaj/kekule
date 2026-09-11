use super::*;

#[test]
fn quadruple_bonds_contribute_four_to_both_endpoint_valences() {
    // RDKit 2026.03.3 UpdatePropertyCache(strict=True) assigns explicit valence
    // four and zero implicit H to both carbons connected by a QUADRUPLE bond.
    let mut editor = MoleculeEditor::new();
    let carbon = Atom::new(Element::from_symbol("C").unwrap());
    let left = editor.add_atom(carbon.clone()).unwrap();
    let right = editor.add_atom(carbon).unwrap();
    let bond = editor.add_bond(left, right, BondOrder::Quadruple).unwrap();
    let mut molecule = editor.finish().unwrap();
    let represented = represented_molecule_snapshot(&molecule);

    valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
        .expect("RDKit accepts valence four at each endpoint");

    assert_eq!(molecule.implicit_hydrogens(left), Ok(Some(0)));
    assert_eq!(molecule.implicit_hydrogens(right), Ok(Some(0)));
    assert_eq!(molecule.bond(bond).unwrap().order, BondOrder::Quadruple);
    assert_eq!(represented_molecule_snapshot(&molecule), represented);
}

fn isolated_valence_state(
    symbol: &str,
    formal_charge: i8,
    radical: AtomRadical,
    hydrogens: HydrogenDeclaration,
) -> (Molecule, AtomId) {
    let mut editor = MoleculeEditor::new();
    let mut atom = Atom::new(Element::from_symbol(symbol).unwrap());
    atom.formal_charge = formal_charge;
    atom.radical = Some(radical);
    atom.hydrogens = hydrogens;
    let id = editor.add_atom(atom).unwrap();
    (editor.finish().unwrap(), id)
}

#[test]
fn disabling_implicit_hydrogens_skips_radical_occupancy_checks() {
    // RDKit 2026.03.3 Atom::UpdatePropertyCache accepts represented CH4 with
    // an explicitly assigned radical when noImplicit is true. It rejects the
    // same occupancy when implicit-valence calculation is enabled.
    for (hydrogens, succeeds) in [
        (HydrogenDeclaration::Fixed(4), true),
        (HydrogenDeclaration::Infer { explicit: 4 }, false),
    ] {
        let (mut molecule, atom) = isolated_valence_state("C", 0, AtomRadical::Doublet, hydrogens);
        let represented = represented_molecule_snapshot(&molecule);
        let result = valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike);
        assert_eq!(result.is_ok(), succeeds);
        assert_eq!(represented_molecule_snapshot(&molecule), represented);
        if succeeds {
            assert_eq!(molecule.implicit_hydrogens(atom), Ok(Some(0)));
        } else {
            assert_eq!(
                result.unwrap_err().issues,
                [ValenceIssue::ValenceOccupancyExceeded {
                    atom,
                    explicit_valence: 4,
                    radical_electrons: 1,
                    charge_offset: 0,
                    max_allowed: 4,
                }]
            );
            assert_eq!(molecule.perception(), &Perception::default());
        }
    }
}

#[test]
fn strict_valence_rejects_excess_occupancy_even_without_bonds_or_declared_hydrogen() {
    for (symbol, charge, radical, hydrogens, charge_offset, max_allowed) in [
        (
            "O",
            0,
            AtomRadical::Quintet,
            HydrogenDeclaration::Infer { explicit: 0 },
            0,
            2,
        ),
        (
            "P",
            -8,
            AtomRadical::Singlet,
            HydrogenDeclaration::Fixed(0),
            8,
            5,
        ),
    ] {
        let (mut molecule, atom) = isolated_valence_state(symbol, charge, radical, hydrogens);
        let error = valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
            .expect_err("RDKit rejects the occupancy even when represented valence is zero");
        assert_eq!(
            error.issues,
            [ValenceIssue::ValenceOccupancyExceeded {
                atom,
                explicit_valence: 0,
                radical_electrons: usize::from(radical.unpaired_electron_count()),
                charge_offset,
                max_allowed,
            }]
        );
        assert_eq!(molecule.perception(), &Perception::default());
        valence_api::perceive_valence_with_options(
            &mut molecule,
            ValenceModel::RdkitLike,
            ValenceOptions { strict: false },
        )
        .unwrap();
        assert_eq!(molecule.implicit_hydrogens(atom), Ok(Some(0)));
    }
}

#[test]
fn original_unrestricted_and_noble_gas_valences_do_not_gain_radical_limits() {
    for (symbol, radical) in [("Li", AtomRadical::Singlet), ("He", AtomRadical::Doublet)] {
        let (mut molecule, atom) = isolated_valence_state(
            symbol,
            1,
            radical,
            HydrogenDeclaration::Infer { explicit: 1 },
        );
        valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
            .expect("RDKit preserves the original-element implicit-valence exemption");
        assert_eq!(molecule.implicit_hydrogens(atom), Ok(Some(0)));
    }
}

#[test]
fn isolated_hydrogen_charge_validation_is_specific_to_hydrogen_inference() {
    for charge in [-2, 2] {
        let (mut molecule, atom) = isolated_valence_state(
            "H",
            charge,
            AtomRadical::Singlet,
            HydrogenDeclaration::Infer { explicit: 0 },
        );
        let error =
            valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike).unwrap_err();
        assert_eq!(
            error.issues,
            [ValenceIssue::InvalidFormalCharge {
                atom,
                formal_charge: charge
            }]
        );
        valence_api::perceive_valence_with_options(
            &mut molecule,
            ValenceModel::RdkitLike,
            ValenceOptions { strict: false },
        )
        .unwrap();
        assert_eq!(molecule.implicit_hydrogens(atom), Ok(Some(0)));

        let (mut fixed, atom) = isolated_valence_state(
            "H",
            charge,
            AtomRadical::Singlet,
            HydrogenDeclaration::Fixed(0),
        );
        valence_api::perceive_valence(&mut fixed, ValenceModel::RdkitLike).unwrap();
        assert_eq!(fixed.implicit_hydrogens(atom), Ok(Some(0)));
    }
}

#[test]
fn hypervalent_anions_mapping_to_unrestricted_elements_do_not_infer_hydrogen() {
    for symbol in ["S", "Se"] {
        let (mut molecule, atom) = isolated_valence_state(
            symbol,
            -5,
            AtomRadical::Singlet,
            HydrogenDeclaration::Infer { explicit: 0 },
        );
        valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike).unwrap();
        assert_eq!(molecule.implicit_hydrogens(atom), Ok(Some(0)));
    }
}

#[test]
fn extreme_charge_adjustments_are_bounded_without_overflow() {
    for charge in [i8::MIN, i8::MAX] {
        let (mut molecule, atom) = isolated_valence_state(
            "C",
            charge,
            AtomRadical::Singlet,
            HydrogenDeclaration::Infer { explicit: 0 },
        );
        valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike).unwrap();
        assert_eq!(molecule.implicit_hydrogens(atom), Ok(Some(0)));
    }
}

#[test]
fn hydride_explicit_valence_compatibility_does_not_override_implicit_valence_rules() {
    for (hydrogens, succeeds) in [
        (HydrogenDeclaration::Fixed(2), true),
        (HydrogenDeclaration::Infer { explicit: 2 }, false),
    ] {
        let (mut molecule, atom) = isolated_valence_state("H", -1, AtomRadical::Singlet, hydrogens);
        let result = valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike);
        assert_eq!(result.is_ok(), succeeds);
        valence_api::perceive_valence_with_options(
            &mut molecule,
            ValenceModel::RdkitLike,
            ValenceOptions { strict: false },
        )
        .unwrap();
        assert_eq!(molecule.implicit_hydrogens(atom), Ok(Some(0)));
    }
}

#[test]
fn valence_accepts_aromatic_input_localized_during_interpretation() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should interpret");
    assert_eq!(
        molecule
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        3
    );
    assert_eq!(
        molecule
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Single)
            .count(),
        3
    );

    valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
        .expect("localized benzene valence should succeed");
    assert!(molecule
        .atom_ids()
        .all(|atom| molecule.implicit_hydrogens(atom) == Ok(Some(1))));
}

#[test]
fn localized_aromatic_valence_replaces_previous_valence_transactionally() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should interpret");
    let atom_ids = molecule.atom_ids().collect::<Vec<_>>();
    let bond_ids = molecule.bond_ids().collect::<Vec<_>>();
    let previous = Perception::builder()
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
        .install_perception(previous.clone())
        .expect("valid previous perception");
    valence_api::perceive_valence_with_options(
        &mut molecule,
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    )
    .expect("localized aromatic valence can be recomputed");

    assert!(molecule
        .atom_ids()
        .all(|atom| molecule.implicit_hydrogens(atom) == Ok(Some(1))));
    assert!(molecule.perception().has_rings());
    assert!(!molecule.perception().has_aromaticity());
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
    mut molecule: Molecule,
    expected_implicit_hydrogens: &[u8],
    expected_aromatic_atoms: usize,
) {
    let smiles = label;
    assert!(!molecule.perception().has_aromaticity(), "{smiles}");

    molecule
        .canonicalize_fixture()
        .unwrap_or_else(|error| panic!("aromatic fixture should normalize: {smiles}: {error}"));
    assert!(molecule
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
    assert_eq!(molecule.perception(), &Perception::default(), "{smiles}");

    valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
        .unwrap_or_else(|error| panic!("localized valence should succeed: {smiles}: {error}"));
    assert!(molecule.perception().has_valence(), "{smiles}");
    assert!(!molecule.perception().has_rings(), "{smiles}");
    assert!(!molecule.perception().has_aromaticity(), "{smiles}");
    assert_eq!(
        molecule
            .atom_ids()
            .map(|atom| {
                molecule
                    .implicit_hydrogens(atom)
                    .expect("live atom")
                    .expect("complete valence assignment")
            })
            .collect::<Vec<_>>(),
        expected_implicit_hydrogens,
        "{smiles}"
    );

    rings_api::perceive_ring_set(&mut molecule)
        .unwrap_or_else(|error| panic!("ring perception should succeed: {smiles}: {error}"));
    aromaticity_api::perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike)
        .unwrap_or_else(|error| panic!("aromaticity perception should succeed: {smiles}: {error}"));
    assert_eq!(
        molecule
            .atom_ids()
            .filter(|atom| molecule.atom_is_aromatic(*atom) == Ok(Some(true)))
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
        let mut radical_carbon = radical.atom_mut(AtomId::new(0)).expect("radical carbon");
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
    let represented = represented_molecule_snapshot(&molecule);
    let represented_nitrogen = molecule
        .atom(AtomId::new(0))
        .expect("pyrrole nitrogen")
        .clone();

    valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
        .expect("pyrrole valence should succeed without aromaticity");

    let nitrogen = molecule.atom(AtomId::new(0)).expect("pyrrole nitrogen");
    assert_eq!(nitrogen.hydrogens, HydrogenDeclaration::Fixed(1));
    assert_eq!(nitrogen.hydrogens, represented_nitrogen.hydrogens);
    assert_eq!(
        molecule
            .implicit_hydrogens(AtomId::new(0))
            .expect("live nitrogen"),
        Some(0)
    );
    assert!(!molecule.perception().has_aromaticity());
    assert_eq!(represented_molecule_snapshot(&molecule), represented);

    rings_api::perceive_ring_set(&mut molecule).expect("pyrrole ring perception");
    aromaticity_api::perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike)
        .expect("pyrrole aromaticity perception");

    assert_eq!(represented_molecule_snapshot(&molecule), represented);
    assert_eq!(
        molecule
            .atom_ids()
            .filter(|atom| molecule.atom_is_aromatic(*atom) == Ok(Some(true)))
            .count(),
        5
    );
    let total_hydrogens = molecule
        .atoms()
        .map(|(atom_id, atom)| {
            usize::from(atom.hydrogens.explicit_count())
                + usize::from(
                    molecule
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
    let aromatic_atoms = with_aromaticity.atom_ids().collect::<Vec<_>>();
    let aromatic_bonds = with_aromaticity.bond_ids().collect::<Vec<_>>();
    let previous = Perception::builder()
        .with_aromaticity(AromaticityModel::RdkitLike, aromatic_atoms, aromatic_bonds)
        .expect("valid semantic aromaticity")
        .build();
    with_aromaticity
        .install_perception(previous)
        .expect("valid perception state");

    valence_api::perceive_valence(&mut without_aromaticity, ValenceModel::RdkitLike)
        .expect("valence without aromaticity");
    valence_api::perceive_valence(&mut with_aromaticity, ValenceModel::RdkitLike)
        .expect("valence with preinstalled aromaticity");

    let without = without_aromaticity
        .atom_ids()
        .map(|atom| without_aromaticity.implicit_hydrogens(atom))
        .collect::<Vec<_>>();
    let with = with_aromaticity
        .atom_ids()
        .map(|atom| with_aromaticity.implicit_hydrogens(atom))
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

    valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
        .expect("naphthalene valence should run first");

    let mut peripheral = 0;
    let mut fused = 0;
    for atom_id in molecule.atom_ids() {
        let degree = molecule
            .incident_bonds(atom_id)
            .expect("live atom")
            .filter(|(_, bond)| !matches!(bond.order, BondOrder::Zero | BondOrder::Dative))
            .count();
        let implicit = molecule
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

    rings_api::perceive_ring_set(&mut molecule).expect("naphthalene rings");
    aromaticity_api::perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike)
        .expect("naphthalene aromaticity");
    assert_eq!(
        molecule
            .atom_ids()
            .filter(|atom| molecule.atom_is_aromatic(*atom) == Ok(Some(true)))
            .count(),
        10
    );
}
