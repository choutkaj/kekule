use super::*;

fn oxo_halide(oxo_count: usize) -> (SmallMolecule, AtomId, Vec<AtomId>, Vec<BondId>) {
    let mut molecule = Molecule::new();
    let chlorine = molecule
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let hydroxyl = molecule
        .add_atom(oxygen())
        .expect("atom identifier capacity");
    molecule
        .add_bond(chlorine, hydroxyl, BondOrder::Single)
        .expect("hydroxyl bond");

    let mut oxygens = Vec::with_capacity(oxo_count);
    let mut bonds = Vec::with_capacity(oxo_count);
    for _ in 0..oxo_count {
        let oxygen = molecule
            .add_atom(oxygen())
            .expect("atom identifier capacity");
        let bond = molecule
            .add_bond(chlorine, oxygen, BondOrder::Double)
            .expect("oxo bond");
        oxygens.push(oxygen);
        bonds.push(bond);
    }
    (
        SmallMolecule::from_molecule(molecule),
        chlorine,
        oxygens,
        bonds,
    )
}

#[test]
fn normalization_rewrites_hypervalent_oxo_halide_representation() {
    let (mut molecule, chlorine, oxygens, bonds) = oxo_halide(2);

    canonicalize_molecule_for_publication(molecule.as_molecule_mut(), &[])
        .expect("normalization should succeed");

    assert_eq!(
        molecule.as_molecule().atom(chlorine).unwrap().formal_charge,
        2
    );
    for oxygen in oxygens {
        assert_eq!(
            molecule.as_molecule().atom(oxygen).unwrap().formal_charge,
            -1
        );
    }
    for bond in bonds {
        assert_eq!(
            molecule.as_molecule().bond(bond).unwrap().order,
            BondOrder::Single
        );
    }
    assert_eq!(molecule.as_molecule().formal_charge(), 0);
}

#[test]
fn normalization_is_idempotent() {
    let (mut molecule, ..) = oxo_halide(2);

    molecule
        .canonicalize_fixture()
        .expect("first normalization should succeed");
    let once = molecule.clone();
    molecule
        .canonicalize_fixture()
        .expect("second normalization should succeed");

    assert_eq!(molecule, once);
}

#[test]
fn source_stereo_is_canonicalized_once_during_interpretation() {
    let (mut molecule, interpretation) =
        read_smiles_with_report("C/C=C\\F").expect("directional SMILES should interpret");

    assert_eq!(interpretation.created_stereo_elements().len(), 1);
    let once = molecule.clone();

    let second = molecule
        .canonicalize_fixture()
        .expect("second source-stereo normalization should succeed");
    assert!(second.created_stereo_elements.is_empty());
    assert!(second.warnings.is_empty());
    assert_eq!(molecule, once);
}

#[test]
fn direct_smiles_tetrahedral_stereo_is_preserved_without_duplication() {
    for smiles in ["F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I"] {
        let mut molecule = read_smiles(smiles).expect("tetrahedral SMILES should parse");
        assert_eq!(molecule.as_molecule().stereo_elements().count(), 1);

        let report = molecule
            .canonicalize_fixture()
            .expect("direct stereo should normalize");

        assert!(report.created_stereo_elements.is_empty());
        assert_eq!(molecule.as_molecule().stereo_elements().count(), 1);
    }
}

#[test]
fn publication_canonicalization_clears_preinstalled_perception() {
    let mut empty_state = read_smiles("C/C=C\\F").expect("directional SMILES should parse");
    let mut installed_state = empty_state.clone();
    mark_all_fresh(installed_state.as_molecule_mut());
    assert_ne!(
        installed_state.as_molecule().perception(),
        &PerceptionState::default()
    );

    let empty_report = empty_state
        .canonicalize_fixture()
        .expect("empty-state normalization should succeed");
    let installed_report = installed_state
        .canonicalize_fixture()
        .expect("installed-state normalization should succeed");

    assert_eq!(installed_report, empty_report);
    assert_eq!(
        installed_state
            .as_molecule()
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>(),
        empty_state
            .as_molecule()
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        installed_state.as_molecule().perception(),
        &PerceptionState::default()
    );
}

#[test]
fn ambiguous_directional_source_marks_return_a_structured_publication_error() {
    let mut graph = Molecule::new();
    let left = graph.add_atom(carbon()).expect("left endpoint");
    let right = graph.add_atom(carbon()).expect("right endpoint");
    let left_a = graph.add_atom(carbon()).expect("left carrier");
    let left_b = graph.add_atom(element_atom("F")).expect("left carrier");
    let right_a = graph.add_atom(element_atom("Cl")).expect("right carrier");
    let double_bond = graph
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    let left_mark_a = graph
        .add_bond(left, left_a, BondOrder::Single)
        .expect("left carrier bond");
    let left_mark_b = graph
        .add_bond(left, left_b, BondOrder::Single)
        .expect("left carrier bond");
    let right_mark = graph
        .add_bond(right, right_a, BondOrder::Single)
        .expect("right carrier bond");
    let source_stereo = [
        (left_mark_a, left),
        (left_mark_b, left),
        (right_mark, right),
    ]
    .into_iter()
    .map(|(bond, from)| SourceStereoBondMark {
        bond,
        from,
        kind: SourceStereoBondMarkKind::DirectionalUp,
    })
    .collect::<Vec<_>>();
    mark_all_fresh(&mut graph);
    let mut molecule = SmallMolecule::from_molecule(graph);
    let error = molecule
        .canonicalize_fixture_with_source_stereo(&source_stereo)
        .expect_err("same-direction marks on both left carriers are ambiguous");

    assert!(matches!(
        error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::AmbiguousDirectionalBondMarks {
                double_bond,
                endpoint: left,
                mark_count: 2,
            })
    ));
}

#[test]
fn aromatic_smiles_is_localized_during_interpretation_then_perceived() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should parse");
    assert!(!molecule.as_molecule().perception().has_aromaticity());
    assert!(molecule
        .as_molecule()
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
    assert_eq!(
        molecule
            .as_molecule()
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        3
    );
    let localized = molecule.as_molecule().clone();

    molecule
        .canonicalize_fixture()
        .expect("benzene should normalize");

    assert_eq!(
        molecule.as_molecule().perception(),
        &PerceptionState::default()
    );
    assert_eq!(molecule.as_molecule(), &localized);
    assert_eq!(
        molecule
            .as_molecule()
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        3
    );

    valence_api::perceive_valence(molecule.as_molecule_mut(), ValenceModel::RdkitLike)
        .expect("localized valence should be perceived");
    aromaticity_api::perceive_aromaticity(molecule.as_molecule_mut(), AromaticityModel::RdkitLike)
        .expect("localized aromaticity should be perceived");

    assert_eq!(
        molecule.as_molecule().perception().aromaticity_model(),
        Some(AromaticityModel::RdkitLike)
    );
    assert!(molecule.as_molecule().atom_ids().all(|atom| molecule
        .as_molecule()
        .atom_is_aromatic(atom)
        .unwrap()
        == Some(true)));
    assert!(molecule.as_molecule().bond_ids().all(|bond| molecule
        .as_molecule()
        .bond_is_aromatic(bond)
        .unwrap()
        == Some(true)));
}

#[test]
fn normalization_preserves_already_localized_aromatic_input() {
    let mut molecule = read_smiles("c1ccc2ccccc2c1").expect("naphthalene should parse");

    molecule
        .canonicalize_fixture()
        .expect("first normalization should succeed");
    let once = molecule.clone();
    molecule
        .canonicalize_fixture()
        .expect("second normalization should succeed");

    assert_eq!(molecule, once);
    assert!(molecule
        .as_molecule()
        .bonds()
        .all(|(_, bond)| { matches!(bond.order, BondOrder::Single | BondOrder::Double) }));
}

#[test]
fn invalid_aromatic_localization_fails_during_interpretation() {
    let error = read_smiles("c1cccc1")
        .expect_err("unlocalizable aromatic source must not publish a molecule");

    assert!(error
        .to_string()
        .contains("invalid imported aromatic representation"));
}

#[test]
fn aromaticity_perception_preserves_complete_primary_representation() {
    let mut molecule =
        read_smiles("C1=CC=CC=C1[C@H](F)C/C=C/Cl").expect("localized stereo fixture should parse");
    molecule
        .as_molecule_mut()
        .props_mut()
        .insert("source".to_owned(), PropValue::String("fixture".to_owned()));
    molecule
        .as_molecule_mut()
        .atom_mut(AtomId::new(0))
        .expect("first atom")
        .props
        .insert("label".to_owned(), PropValue::Int(7));
    let first_bond = molecule
        .as_molecule()
        .bond_ids()
        .next()
        .expect("first bond");
    molecule
        .as_molecule_mut()
        .bond_mut(first_bond)
        .expect("first bond")
        .props
        .insert("source".to_owned(), PropValue::Bool(true));
    let mut conformer = Conformer::new(crate::units::ANGSTROM).expect("length unit");
    for atom_id in molecule.as_molecule().atom_ids() {
        conformer
            .set_position(
                atom_id,
                crate::units::Quantity::new(
                    Point3::new(atom_id.index() as f64, 0.0, 0.0),
                    crate::units::ANGSTROM,
                ),
            )
            .expect("valid position");
    }
    molecule
        .as_molecule_mut()
        .add_conformer(conformer)
        .expect("valid conformer");
    assert!(molecule.as_molecule().stereo_elements().next().is_some());
    assert!(molecule.as_molecule().conformers().next().is_some());
    valence_api::perceive_valence(molecule.as_molecule_mut(), ValenceModel::RdkitLike)
        .expect("valence perception");

    let mut represented_before = molecule.as_molecule().clone();
    clear_test_derived_state(&mut represented_before);
    aromaticity_api::perceive_aromaticity(molecule.as_molecule_mut(), AromaticityModel::RdkitLike)
        .expect("aromaticity perception");
    let mut represented_after = molecule.as_molecule().clone();
    clear_test_derived_state(&mut represented_after);

    assert_eq!(represented_after, represented_before);
    assert!(molecule.as_molecule().perception().has_aromaticity());
}

fn clear_test_derived_state(molecule: &mut Molecule) {
    molecule
        .install_perception_state(PerceptionState::default())
        .expect("empty perception state should install");
}

#[test]
fn successful_normalization_clears_installed_perception() {
    let mut molecule = read_smiles("CCO").expect("ethanol should parse");
    mark_all_fresh(molecule.as_molecule_mut());
    assert_ne!(
        molecule.as_molecule().perception(),
        &PerceptionState::default()
    );

    molecule
        .canonicalize_fixture()
        .expect("normalization should succeed");

    assert_eq!(
        molecule.as_molecule().perception(),
        &PerceptionState::default()
    );
    assert_all_stale(molecule.as_molecule());
}

#[test]
fn failed_normalization_is_transactional() {
    let (mut molecule, chlorine, ..) = oxo_halide(128);
    mark_all_fresh(molecule.as_molecule_mut());
    let before = molecule.clone();

    let error = molecule
        .canonicalize_fixture()
        .expect_err("unrepresentable formal charge should fail");

    assert_eq!(
        error,
        NormalizationError::FormalChargeOutOfRange {
            atom: chlorine,
            charge: 128,
        }
    );
    assert_eq!(molecule, before);
}
