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
fn source_stereo_normalization_is_idempotent_and_consumes_marks() {
    let mut molecule = read_smiles("C/C=C\\F").expect("directional SMILES should parse");

    let first = molecule
        .normalize()
        .expect("first source-stereo normalization should succeed");
    assert_eq!(first.created_stereo_elements.len(), 1);
    assert!(molecule.graph().stereo_bond_marks().next().is_none());
    let once = molecule.clone();

    let second = molecule
        .normalize()
        .expect("second source-stereo normalization should succeed");
    assert!(second.created_stereo_elements.is_empty());
    assert!(second.warnings.is_empty());
    assert_eq!(molecule, once);
}

#[test]
fn direct_smiles_tetrahedral_stereo_is_preserved_without_duplication() {
    for smiles in ["F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I"] {
        let mut molecule = read_smiles(smiles).expect("tetrahedral SMILES should parse");
        assert_eq!(molecule.graph().stereo_elements().count(), 1);

        let report = molecule
            .normalize()
            .expect("direct stereo should normalize");

        assert!(report.created_stereo_elements.is_empty());
        assert_eq!(molecule.graph().stereo_elements().count(), 1);
        assert!(molecule.graph().stereo_bond_marks().next().is_none());
    }
}

#[test]
fn installed_perception_does_not_affect_source_stereo_normalization() {
    let mut empty_state = read_smiles("C/C=C\\F").expect("directional SMILES should parse");
    let mut installed_state = empty_state.clone();
    mark_all_fresh(installed_state.graph_mut());
    assert_ne!(
        installed_state.graph().perception(),
        &PerceptionState::default()
    );

    let empty_report = empty_state
        .normalize()
        .expect("empty-state normalization should succeed");
    let installed_report = installed_state
        .normalize()
        .expect("installed-state normalization should succeed");

    assert_eq!(installed_report, empty_report);
    assert_eq!(
        installed_state
            .graph()
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>(),
        empty_state
            .graph()
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>()
    );
    assert!(installed_state.graph().stereo_bond_marks().next().is_none());
    assert_eq!(
        installed_state.graph().perception(),
        &PerceptionState::default()
    );
}

#[test]
fn ambiguous_directional_source_marks_roll_back_complete_state() {
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
    for bond in [left_mark_a, left_mark_b, right_mark] {
        graph
            .set_stereo_bond_mark(StereoBondMark {
                bond,
                kind: StereoBondMarkKind::DirectionalUp,
                source: StereoSource::Smiles,
            })
            .expect("directional mark");
    }
    mark_all_fresh(&mut graph);
    let mut molecule = SmallMolecule::from_graph(graph);
    let before = molecule.clone();

    let error = molecule
        .normalize()
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
    assert_eq!(molecule, before);
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
    molecule
        .install_perception_state(PerceptionState::default())
        .expect("empty perception state should install");
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
