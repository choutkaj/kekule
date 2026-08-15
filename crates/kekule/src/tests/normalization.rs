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
        SmallMolecule::from_graph(molecule),
        chlorine,
        oxygens,
        bonds,
    )
}

#[test]
fn normalization_rewrites_hypervalent_oxo_halide_representation() {
    let (mut molecule, chlorine, oxygens, bonds) = oxo_halide(2);

    normalization_api::normalize(molecule.graph_mut()).expect("normalization should succeed");

    assert_eq!(molecule.graph().atom(chlorine).unwrap().formal_charge, 2);
    for oxygen in oxygens {
        assert_eq!(molecule.graph().atom(oxygen).unwrap().formal_charge, -1);
    }
    for bond in bonds {
        assert_eq!(
            molecule.graph().bond(bond).unwrap().order,
            BondOrder::Single
        );
    }
    assert_eq!(molecule.graph().formal_charge(), 0);
}

#[test]
fn normalization_is_idempotent() {
    let (mut molecule, ..) = oxo_halide(2);

    molecule
        .normalize()
        .expect("first normalization should succeed");
    let once = molecule.clone();
    molecule
        .normalize()
        .expect("second normalization should succeed");

    assert_eq!(molecule, once);
}

#[test]
fn aromatic_smiles_normalizes_then_perceives_semantic_aromaticity() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should parse");
    assert!(!molecule.graph().perception().has_aromaticity());
    assert!(molecule
        .graph()
        .bonds()
        .all(|(_, bond)| bond.order == BondOrder::Aromatic));

    molecule.normalize().expect("benzene should normalize");

    assert_eq!(molecule.graph().perception(), &PerceptionState::default());
    assert!(molecule
        .graph()
        .bonds()
        .all(|(_, bond)| bond.order != BondOrder::Aromatic));
    assert_eq!(
        molecule
            .graph()
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        3
    );

    valence_api::perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike)
        .expect("localized valence should be perceived");
    aromaticity_api::perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
        .expect("localized aromaticity should be perceived");

    assert_eq!(
        molecule.graph().perception().aromaticity_model(),
        Some(AromaticityModel::RdkitLike)
    );
    assert!(molecule.graph().atom_ids().all(|atom| molecule
        .graph()
        .atom_is_aromatic(atom)
        .unwrap()
        == Some(true)));
    assert!(molecule.graph().bond_ids().all(|bond| molecule
        .graph()
        .bond_is_aromatic(bond)
        .unwrap()
        == Some(true)));
}

#[test]
fn aromatic_localization_is_idempotent() {
    let mut molecule = read_smiles("c1ccc2ccccc2c1").expect("naphthalene should parse");

    molecule
        .normalize()
        .expect("first aromatic localization should succeed");
    let once = molecule.clone();
    molecule
        .normalize()
        .expect("second aromatic localization should succeed");

    assert_eq!(molecule, once);
    assert!(molecule
        .graph()
        .bonds()
        .all(|(_, bond)| bond.order != BondOrder::Aromatic));
}

#[test]
fn invalid_aromatic_localization_is_transactional() {
    let mut molecule = read_smiles("c1cccc1").expect("source syntax should parse");
    mark_all_fresh(molecule.graph_mut());
    let before = molecule.clone();

    let error = molecule
        .normalize()
        .expect_err("odd aromatic demand should fail localization");

    assert!(matches!(
        error,
        crate::normalization::NormalizationError::InvalidAromaticRepresentation(_)
    ));
    assert_eq!(molecule, before);
}

#[test]
fn aromaticity_perception_preserves_complete_primary_representation() {
    let mut molecule =
        read_smiles("C1=CC=CC=C1[C@H](F)C/C=C/Cl").expect("localized stereo fixture should parse");
    molecule
        .graph_mut()
        .props_mut()
        .insert("source".to_owned(), PropValue::String("fixture".to_owned()));
    molecule
        .graph_mut()
        .atom_mut(AtomId::new(0))
        .expect("first atom")
        .props
        .insert("label".to_owned(), PropValue::Int(7));
    let first_bond = molecule.graph().bond_ids().next().expect("first bond");
    molecule
        .graph_mut()
        .bond_mut(first_bond)
        .expect("first bond")
        .props
        .insert("source".to_owned(), PropValue::Bool(true));
    let mut conformer = Conformer::new(crate::units::ANGSTROM).expect("length unit");
    for atom_id in molecule.graph().atom_ids() {
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
        .graph_mut()
        .add_conformer(conformer)
        .expect("valid conformer");
    assert!(molecule.graph().stereo_elements().next().is_some());
    assert!(molecule.graph().stereo_bond_marks().next().is_some());
    assert!(molecule.graph().conformers().next().is_some());
    valence_api::perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike)
        .expect("valence perception");

    let mut represented_before = molecule.graph().clone();
    clear_test_derived_state(&mut represented_before);
    aromaticity_api::perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
        .expect("aromaticity perception");
    let mut represented_after = molecule.graph().clone();
    clear_test_derived_state(&mut represented_after);

    assert_eq!(represented_after, represented_before);
    assert!(molecule.graph().perception().has_aromaticity());
}

fn clear_test_derived_state(molecule: &mut Molecule) {
    molecule.perception = PerceptionState::default();
    for atom in molecule.atoms.iter_mut().flatten() {
        atom.implicit_hydrogens = None;
        atom.aromatic = false;
    }
    for bond in molecule.bonds.iter_mut().flatten() {
        bond.aromatic = false;
    }
    for element in molecule.stereo_elements.iter_mut().flatten() {
        element.descriptor = None;
    }
}

#[test]
fn successful_normalization_clears_installed_perception() {
    let mut molecule = read_smiles("CCO").expect("ethanol should parse");
    mark_all_fresh(molecule.graph_mut());
    assert_ne!(molecule.graph().perception(), &PerceptionState::default());

    molecule.normalize().expect("normalization should succeed");

    assert_eq!(molecule.graph().perception(), &PerceptionState::default());
    assert_all_stale(molecule.graph());
}

#[test]
fn failed_normalization_is_transactional() {
    let (mut molecule, chlorine, ..) = oxo_halide(128);
    mark_all_fresh(molecule.graph_mut());
    let before = molecule.clone();

    let error = molecule
        .normalize()
        .expect_err("unrepresentable formal charge should fail");

    assert_eq!(
        error,
        crate::normalization::NormalizationError::FormalChargeOutOfRange {
            atom: chlorine,
            charge: 128,
        }
    );
    assert_eq!(molecule, before);
}

#[test]
fn sanitizer_delegates_normalization_failure_transactionally() {
    let (mut molecule, chlorine, ..) = oxo_halide(128);
    let before = molecule.clone();

    let error = perception_api::sanitize_with_options(
        &mut molecule,
        SanitizeOptions {
            perceive_valence: false,
            perceive_rings: false,
            perceive_aromaticity: false,
            perceive_stereo: false,
        },
    )
    .expect_err("normalization failure should fail sanitization");

    assert_eq!(
        error,
        SanitizeError::Normalization(
            crate::normalization::NormalizationError::FormalChargeOutOfRange {
                atom: chlorine,
                charge: 128,
            }
        )
    );
    assert_eq!(molecule, before);
}
