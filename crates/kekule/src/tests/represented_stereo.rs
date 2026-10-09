//! Validation, candidate detection and installation of represented stereo.

use super::*;

#[test]
fn stereo_validation_reports_invalid_local_elements_without_mutating() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let a = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let b = mol
        .add_atom(element_atom("N"))
        .expect("atom identifier capacity");
    mol.add_bond(center, a, BondOrder::Single).expect("bond");
    mark_all_fresh(mol.working_mut());
    let element = insert_unchecked_stereo(
        &mut mol,
        StereoElement {
            kind: StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(a),
                    StereoCarrier::Atom(a),
                    StereoCarrier::Atom(a),
                ],
                orientation: None,
            }),
            group: None,
        },
    );
    // Inject nonadjacency through internal storage to test the diagnostic;
    // checked editing now rejects this malformed reference.
    let stored = mol.working_mut().graph.stereo_elements[element.index()]
        .as_mut()
        .expect("stored stereo element");
    let StereoElementKind::Tetrahedral(stereo) = &mut stored.kind else {
        unreachable!("test element is tetrahedral");
    };
    stereo.carriers[2] = StereoCarrier::Atom(b);
    mark_all_fresh(mol.working_mut());

    let error = stereo_api::validate_stereo(mol.working()).expect_err("invalid stored stereo");

    assert!(mol.stereo_elements().next().is_some());
    assert!(error
        .issues
        .contains(&StereoValidationIssue::InvalidTetrahedralCarrierCount {
            element,
            center,
            carrier_count: 3,
        }));
    assert!(error
        .issues
        .contains(&StereoValidationIssue::DuplicateTetrahedralCarrier {
            element,
            center,
            carrier: StereoCarrier::Atom(a),
        }));
    assert!(error
        .issues
        .contains(&StereoValidationIssue::TetrahedralCarrierNotAdjacent {
            element,
            center,
            carrier: StereoCarrier::Atom(b),
        }));
    assert!(mol
        .stereo_element(element)
        .expect("element")
        .is_explicitly_unknown());
}

#[test]
fn stereo_validation_checks_implicit_carrier_form_without_perception_state() {
    let mut tetrahedral = crate::core::MoleculeEditor::new();
    let center = tetrahedral
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let mut atom_carriers = Vec::new();
    for symbol in ["F", "Cl", "Br"] {
        let carrier = tetrahedral
            .add_atom(element_atom(symbol))
            .expect("atom identifier capacity");
        tetrahedral
            .add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
        atom_carriers.push(StereoCarrier::Atom(carrier));
    }
    let mut carriers = atom_carriers.clone();
    carriers.push(StereoCarrier::ImplicitHydrogen);
    let hydrogen_element = tetrahedral
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers,
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("tetrahedral stereo element");

    stereo_api::validate_stereo(tetrahedral.working())
        .expect("implicit hydrogen availability is chemically interpretive");

    tetrahedral
        .remove_stereo_element(hydrogen_element)
        .expect("remove hydrogen-carrier element");
    atom_carriers.push(StereoCarrier::ImplicitLonePair);
    tetrahedral
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: atom_carriers,
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("tetrahedral stereo element");
    stereo_api::validate_stereo(tetrahedral.working())
        .expect("implicit lone-pair availability is chemically interpretive");

    let mut double_bond = crate::core::MoleculeEditor::new();
    let left = double_bond
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let right = double_bond
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let bond = double_bond
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    let double_element = insert_unchecked_stereo(
        &mut double_bond,
        StereoElement::new(StereoElementKind::DoubleBond(DoubleBondStereo {
            bond,
            left,
            right,
            left_carrier: StereoCarrier::ImplicitHydrogen,
            right_carrier: StereoCarrier::ImplicitLonePair,
            orientation: Some(DoubleBondOrientation::Together),
        })),
    );

    let error = stereo_api::validate_stereo(double_bond.working())
        .expect_err("unavailable double-bond carriers should be reported");
    assert_eq!(
        error.issues,
        vec![StereoValidationIssue::UnsupportedDoubleBondCarrier {
            element: double_element,
            endpoint: right,
            carrier: StereoCarrier::ImplicitLonePair,
        }]
    );

    let mut axis = crate::core::MoleculeEditor::new();
    let axis_left = axis.add_atom(carbon()).expect("atom identifier capacity");
    let axis_right = axis.add_atom(carbon()).expect("atom identifier capacity");
    let axis_bond = axis
        .add_bond(axis_left, axis_right, BondOrder::Single)
        .expect("axis bond");
    let axis_element = insert_unchecked_stereo(
        &mut axis,
        StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis: axis_bond,
            carriers: vec![
                StereoCarrier::ImplicitHydrogen,
                StereoCarrier::ImplicitLonePair,
            ],
            orientation: Some(AxisOrientation::Clockwise),
        })),
    );

    let error = stereo_api::validate_stereo(axis.working())
        .expect_err("implicit axis carriers should be unsupported");
    assert_eq!(
        error.issues,
        vec![
            StereoValidationIssue::UnsupportedAxisCarrier {
                element: axis_element,
                axis: axis_bond,
                carrier: StereoCarrier::ImplicitHydrogen,
            },
            StereoValidationIssue::UnsupportedAxisCarrier {
                element: axis_element,
                axis: axis_bond,
                carrier: StereoCarrier::ImplicitLonePair,
            },
        ]
    );
}

#[test]
fn stereo_candidates_use_normalized_and_perceived_hydrogen_state_without_cip_assignment() {
    let mut molecule = read_smiles("CC(F)(Cl)Br").expect("smiles should parse");
    perceive(&mut molecule).expect("molecule should perceive");

    stereo_api::validate_stereo(&molecule).expect("stored stereo should be valid");
    let candidates = stereo_api::detect_stereo_candidates(&molecule).unwrap();

    assert!(candidates.iter().any(|candidate| matches!(
        candidate,
        StereoCandidate::Tetrahedral { center, carriers }
            if *center == AtomId::new(1)
                && carriers.len() == 4
                && !carriers.contains(&StereoCarrier::ImplicitHydrogen)
    )));
    assert!(molecule.stereo_elements().next().is_none());
}

#[test]
fn stereo_candidates_exclude_repeated_hydrogen_ligands_but_preserve_isotopes() {
    for (input, expected) in [
        ("C(F)Cl", 0),
        ("[H]C([H])(F)Cl", 0),
        ("[H]C(F)Cl", 0),
        ("[2H]C([2H])(F)Cl", 0),
        ("[H]C([2H])(F)Cl", 1),
        ("[2H]C([3H])(F)Cl", 1),
        ("C=CF", 0),
        ("[H]C([H])=CF", 0),
        ("[H]C([2H])=CF", 1),
    ] {
        let mut molecule = read_smiles(input).unwrap();
        perceive(&mut molecule).unwrap();
        assert_eq!(
            stereo_api::detect_stereo_candidates(&molecule)
                .unwrap()
                .len(),
            expected,
            "{input}"
        );
    }
}

#[test]
fn stereo_candidates_include_pyramidal_pnictogens_and_tetracoordinate_multiple_bonds() {
    for (input, expected, virtual_h, lone_pair) in [
        ("P(C)(CC)CCC", true, false, true),
        ("[PH](C)CC", true, true, true),
        ("[As](C)(CC)CCC", true, false, true),
        ("[AsH](C)CC", true, true, true),
        ("P(=O)(C)(CC)CCC", true, false, false),
        ("S(=O)(C)CC", true, false, true),
        ("[S+](C)(CC)CCC", true, false, true),
        ("[S](C)(CC)CCC", false, false, false),
        ("[Se](C)(CC)CCC", false, false, false),
        ("S(=O)(=O)(C)CC", false, false, false),
    ] {
        let mut molecule = read_smiles(input).unwrap();
        perceive(&mut molecule).unwrap();
        let candidates = stereo_api::detect_stereo_candidates(&molecule).unwrap();
        let carriers = candidates.iter().find_map(|candidate| match candidate {
            StereoCandidate::Tetrahedral { center, carriers } if *center == AtomId::new(0) => {
                Some(carriers)
            }
            _ => None,
        });
        assert_eq!(carriers.is_some(), expected, "{input}");
        if let Some(carriers) = carriers {
            assert_eq!(carriers.len(), 4);
            assert_eq!(
                carriers.contains(&StereoCarrier::ImplicitHydrogen),
                virtual_h,
                "{input}"
            );
            assert_eq!(
                carriers.contains(&StereoCarrier::ImplicitLonePair),
                lone_pair,
                "{input}"
            );
        }
    }
}

#[test]
fn terminal_ligand_equivalence_is_hydrogen_representation_invariant() {
    for (input, expected) in [
        ("C(C)(C)(F)Cl", 0),
        ("C(F)(F)(Cl)Br", 0),
        ("C(C)([13CH3])(F)Cl", 1),
        ("C(C)(C[2H])(F)Cl", 1),
        ("C([C@]([H])([2H])[3H])([C@@]([H])([2H])[3H])(F)Cl", 3),
        ("C([CH3:1])([CH3:2])(F)Cl", 0),
        ("CC(C)=C(F)Cl", 0),
        ("CC([13CH3])=C(F)Cl", 1),
        ("CC=N", 1),
        ("CC=[NH]", 1),
        ("CC=O", 0),
        ("CP(=O)(C)CC", 0),
    ] {
        let mut molecule = read_smiles(input).unwrap();
        perceive(&mut molecule).unwrap();
        let before = molecule.clone();
        assert_eq!(
            stereo_api::detect_stereo_candidates(&molecule)
                .unwrap()
                .len(),
            expected,
            "{input}"
        );
        assert_eq!(molecule, before, "candidate detection must be read-only");
        molecule.add_hydrogens().unwrap();
        perceive(&mut molecule).unwrap();
        assert_eq!(
            stereo_api::detect_stereo_candidates(&molecule)
                .unwrap()
                .len(),
            expected,
            "expanded {input}"
        );
    }
}

#[test]
fn ring_double_bond_eligibility_depends_on_size_not_carbon_endpoints() {
    for (input, expected) in [
        ("C1=NCCCCC1", 0),
        ("C1=NCCCCCC1", 1),
        ("C1=NNCCCCC1", 1),
        ("N1=NCCCCCC1", 1),
        ("C1=CCCCCCC1", 1),
    ] {
        let mut molecule = read_smiles(input).unwrap();
        perceive(&mut molecule).unwrap();
        let doubles = stereo_api::detect_stereo_candidates(&molecule)
            .unwrap()
            .into_iter()
            .filter(|candidate| matches!(candidate, StereoCandidate::DoubleBond { .. }))
            .count();
        assert_eq!(doubles, expected, "{input}");
    }
    let source = read_smiles(r"C1/N=C\CCCCC1").expect("eight-membered imine source geometry");
    assert_eq!(source.stereo_elements().count(), 1);
    assert!(
        read_smiles(r"C1/N=C\CCCC1").is_err(),
        "seven-membered ring remains excluded"
    );
}

#[test]
fn perception_installation_rejects_a_descriptor_for_the_wrong_stereo_geometry() {
    let (axis, positions, _) = coordinate_axis_graph(true);
    let mut editor = axis.into_editor();
    stereo_api::materialize_coordinate_stereo_with_options(
        &mut editor,
        &positions,
        CoordinateStereoOptions {
            infer_axes: true,
            ..Default::default()
        },
    )
    .unwrap();
    let axis = editor.finish().unwrap();
    for (mut molecule, accepted, rejected) in [
        (
            read_smiles("F[C@](Cl)(Br)I").unwrap(),
            StereoDescriptor::R,
            StereoDescriptor::E,
        ),
        (
            read_smiles("F/C=C/Cl").unwrap(),
            StereoDescriptor::E,
            StereoDescriptor::R,
        ),
        (axis, StereoDescriptor::M, StereoDescriptor::LowerR),
    ] {
        let element = molecule.stereo_element_ids().next().unwrap();
        let previous = Perception::builder()
            .with_cip_descriptors(vec![(element, accepted)])
            .unwrap()
            .build();
        molecule.install_perception(previous.clone()).unwrap();
        let before = molecule.clone();
        let incompatible = Perception::builder()
            .with_cip_descriptors(vec![(element, rejected)])
            .unwrap()
            .build();
        let error = molecule.install_perception(incompatible).unwrap_err();
        assert_eq!(
            error,
            PerceptionInstallError::IncompatibleStereoDescriptor {
                element,
                descriptor: rejected
            }
        );
        assert!(error.to_string().contains("incompatible"));
        assert_eq!(molecule, before);
        assert_eq!(molecule.perception(), &previous);
    }
}

#[test]
fn stereo_validation_accepts_structural_axis_elements() {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left_carrier = mol
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");
    let right_carrier = mol
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let axis = mol.add_bond(left, right, BondOrder::Single).expect("axis");
    mol.add_bond(left, left_carrier, BondOrder::Single)
        .expect("left carrier");
    mol.add_bond(right, right_carrier, BondOrder::Single)
        .expect("right carrier");
    let valid_axis = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![
                StereoCarrier::Atom(left_carrier),
                StereoCarrier::Atom(right_carrier),
            ],
            orientation: Some(AxisOrientation::CounterClockwise),
        })))
        .expect("axis element");

    stereo_api::validate_stereo(mol.working()).expect("axis should be structurally valid");

    mol.remove_stereo_element(valid_axis)
        .expect("remove valid axis");
    let invalid_axis = insert_unchecked_stereo(
        &mut mol,
        StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![StereoCarrier::Atom(left_carrier)],
            orientation: Some(AxisOrientation::CounterClockwise),
        })),
    );

    let error = stereo_api::validate_stereo(mol.working()).expect_err("axis should be invalid");

    assert_eq!(
        error.issues,
        vec![StereoValidationIssue::InvalidAxisCarrierCount {
            element: invalid_axis,
            axis,
            carrier_count: 1,
        }]
    );
}
