use super::*;
use crate::properties::{PropertyKey, PropertyValue};

fn oxo_halide(oxo_count: usize) -> (MoleculeEditor, AtomId, Vec<AtomId>, Vec<BondId>) {
    let mut molecule = crate::core::MoleculeEditor::new();
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
    (molecule, chlorine, oxygens, bonds)
}

#[test]
fn normalization_rewrites_hypervalent_oxo_halide_representation() {
    let (mut molecule, chlorine, oxygens, bonds) = oxo_halide(2);

    canonicalize_molecule_for_publication(molecule.working_mut(), None, &[])
        .expect("normalization should succeed");

    assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 2);
    for oxygen in oxygens {
        assert_eq!(molecule.atom(oxygen).unwrap().formal_charge, -1);
    }
    for bond in bonds {
        assert_eq!(molecule.bond(bond).unwrap().order, BondOrder::Single);
    }
    assert_eq!(molecule.formal_charge(), 0);
}

#[test]
fn normalization_is_idempotent() {
    let (mut molecule, ..) = oxo_halide(2);

    molecule
        .working_mut()
        .canonicalize_fixture()
        .expect("first normalization should succeed");
    let once = molecule.clone();
    molecule
        .working_mut()
        .canonicalize_fixture()
        .expect("second normalization should succeed");

    assert_eq!(molecule, once);
}

#[test]
fn oxohalogen_cleanup_matches_rdkit_for_ester_and_declared_hydrogen_forms() {
    // Full represented charges and bond orders verified against RDKit 2026.03.3.
    for (source, charges, implicit_hydrogens) in [
        (
            "CC(C)(C)OCl(=O)(=O)=O",
            vec![0, 0, 0, 0, 0, 3, -1, -1, -1],
            vec![3, 0, 3, 3, 0, 0, 0, 0, 0],
        ),
        ("[IH]=O", vec![1, -1], vec![0, 0]),
        ("[O-]=I(O)", vec![-1, 1, 0], vec![0, 0, 1]),
    ] {
        let mut molecule = read_smiles(source).expect(source);
        assert_eq!(
            molecule
                .atoms()
                .map(|(_, atom)| atom.formal_charge)
                .collect::<Vec<_>>(),
            charges,
            "{source}"
        );
        assert!(molecule
            .bonds()
            .all(|(_, bond)| bond.order == BondOrder::Single));
        assert_eq!(molecule.perception(), &Perception::default());
        let represented = molecule.clone();
        molecule.perceive().expect(source);
        assert_eq!(
            molecule
                .atom_ids()
                .map(|atom| molecule.implicit_hydrogens(atom).unwrap().unwrap())
                .collect::<Vec<_>>(),
            implicit_hydrogens,
            "{source}"
        );
        assert_eq!(molecule, represented, "perception must preserve {source}");
        assert_eq!(
            molecule.edit().finish().unwrap(),
            represented,
            "idempotent {source}"
        );
    }
}

#[test]
fn oxohalogen_cleanup_requires_only_oxygen_neighbors_and_odd_modeled_valence() {
    for source in ["OI(=O)(C)C", "OI(=O)O"] {
        let mut molecule = read_smiles(source).expect(source);
        assert!(
            molecule.atoms().all(|(_, atom)| atom.formal_charge == 0),
            "{source}"
        );
        assert_eq!(
            molecule
                .bonds()
                .filter(|(_, bond)| bond.order == BondOrder::Double)
                .count(),
            1,
            "{source}"
        );
        molecule.perceive().expect(source);
    }
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
        assert_eq!(molecule.stereo_elements().count(), 1);

        let report = molecule
            .canonicalize_fixture()
            .expect("direct stereo should normalize");

        assert!(report.created_stereo_elements.is_empty());
        assert_eq!(molecule.stereo_elements().count(), 1);
    }
}

#[test]
fn publication_canonicalization_clears_preinstalled_perception() {
    let mut empty_state = read_smiles("C/C=C\\F").expect("directional SMILES should parse");
    let mut installed_state = empty_state.clone();
    mark_all_fresh(&mut installed_state);
    assert_ne!(installed_state.perception(), &Perception::default());

    let empty_report = empty_state
        .canonicalize_fixture()
        .expect("empty-state normalization should succeed");
    let installed_report = installed_state
        .canonicalize_fixture()
        .expect("installed-state normalization should succeed");

    assert_eq!(installed_report, empty_report);
    assert_eq!(
        installed_state
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>(),
        empty_state
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>()
    );
    assert_eq!(installed_state.perception(), &Perception::default());
}

#[test]
fn ambiguous_directional_source_marks_return_a_structured_publication_error() {
    let mut graph = crate::core::MoleculeEditor::new();
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
    mark_all_fresh(graph.working_mut());
    let mut molecule = graph;
    let error = molecule
        .working_mut()
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
    // This is interpreter-owned staging: its documented failure contract is
    // discard-on-error, and arbitrary perception is cleared before decoding.
    assert_eq!(molecule.perception(), &Perception::default());
}

#[test]
fn aromatic_smiles_is_localized_during_interpretation_then_perceived() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should parse");
    assert!(!molecule.perception().has_aromaticity());
    assert!(molecule
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
    assert_eq!(
        molecule
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        3
    );
    let localized = molecule.clone();

    molecule
        .canonicalize_fixture()
        .expect("benzene should normalize");

    assert_eq!(molecule.perception(), &Perception::default());
    assert_eq!(molecule, localized);
    assert_eq!(
        molecule
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        3
    );

    valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
        .expect("localized valence should be perceived");
    aromaticity_api::perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike)
        .expect("localized aromaticity should be perceived");

    assert_eq!(
        molecule.perception().aromaticity_model(),
        Some(AromaticityModel::RdkitLike)
    );
    assert!(molecule
        .atom_ids()
        .all(|atom| molecule.atom_is_aromatic(atom).unwrap() == Some(true)));
    assert!(molecule
        .bond_ids()
        .all(|bond| molecule.bond_is_aromatic(bond).unwrap() == Some(true)));
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
        .insert_property(
            PropertyKey::new("source").unwrap(),
            PropertyValue::String("fixture".to_owned()),
        )
        .unwrap();
    molecule
        .set_atom_property(
            AtomId::new(0),
            PropertyKey::new("label").unwrap(),
            Some(PropertyValue::Int(7)),
        )
        .unwrap();
    let first_bond = molecule.bond_ids().next().expect("first bond");
    molecule
        .set_bond_property(
            first_bond,
            PropertyKey::new("source").unwrap(),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    let positions = crate::structure::Positions::new(crate::units::Quantity::new(
        molecule
            .atom_ids()
            .map(|atom_id| Point3::new(atom_id.index() as f64, 0.0, 0.0))
            .collect::<Vec<_>>(),
        crate::units::ANGSTROM,
    ))
    .expect("valid positions");
    assert!(molecule.stereo_elements().next().is_some());
    assert_eq!(positions.len(), molecule.atom_count());
    valence_api::perceive_valence(&mut molecule, ValenceModel::RdkitLike)
        .expect("valence perception");

    let mut represented_before = molecule.clone();
    clear_test_derived_state(&mut represented_before);
    aromaticity_api::perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike)
        .expect("aromaticity perception");
    let mut represented_after = molecule.clone();
    clear_test_derived_state(&mut represented_after);

    assert_eq!(represented_after, represented_before);
    assert!(molecule.perception().has_aromaticity());
}

fn clear_test_derived_state(molecule: &mut Molecule) {
    molecule
        .install_perception(Perception::default())
        .expect("empty perception state should install");
}

#[test]
fn successful_normalization_clears_installed_perception() {
    let mut molecule = read_smiles("CCO").expect("ethanol should parse");
    mark_all_fresh(&mut molecule);
    assert_ne!(molecule.perception(), &Perception::default());

    molecule
        .canonicalize_fixture()
        .expect("normalization should succeed");

    assert_eq!(molecule.perception(), &Perception::default());
    assert_all_stale(&molecule);
}

#[test]
fn normalization_does_not_invent_large_oxohalogen_charges() {
    for oxo_count in [4, 128] {
        let (mut molecule, chlorine, oxygens, bonds) = oxo_halide(oxo_count);
        molecule
            .working_mut()
            .canonicalize_fixture()
            .expect("out-of-model chemistry is preserved");
        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 0);
        assert!(oxygens
            .iter()
            .all(|&oxygen| molecule.atom(oxygen).unwrap().formal_charge == 0));
        assert!(bonds
            .iter()
            .all(|&bond| molecule.bond(bond).unwrap().order == BondOrder::Double));
        assert!(
            valence_api::perceive_valence(molecule.working_mut(), ValenceModel::RdkitLike).is_err()
        );
    }
}
