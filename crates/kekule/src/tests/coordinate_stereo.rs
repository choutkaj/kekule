//! Stereo inferred from 2D drawings and 3D coordinates.

use super::*;

#[test]
fn coordinate_stereo_inference_is_read_only_and_materializes_tetrahedral_stereo() {
    let (mut mol, center, carriers, _) = tetrahedral_marked_graph();
    let positions = test_positions(vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 0.0, -1.0),
    ]);
    mark_all_fresh(mol.working_mut());
    let before = mol.clone();

    let inferred = stereo_api::infer_coordinate_stereo(mol.working(), &positions)
        .expect("3D tetrahedral stereo should be inferred");
    assert_eq!(mol, before);
    assert_eq!(inferred.elements.len(), 1);
    let proposed = &inferred.elements[0];
    match &proposed.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assert_eq!(stereo.center, center);
            assert_eq!(
                stereo.carriers,
                carriers
                    .iter()
                    .copied()
                    .map(StereoCarrier::Atom)
                    .collect::<Vec<_>>()
            );
            assert_eq!(stereo.orientation, Some(TetrahedralOrientation::Clockwise));
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }

    let report = stereo_api::materialize_coordinate_stereo(&mut mol, &positions)
        .expect("3D tetrahedral stereo should materialize");
    assert_eq!(report.created_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert_eq!(element, proposed);
}

#[test]
fn coordinate_stereo_resolves_tombstoned_atom_ids_to_dense_positions() {
    let mut molecule = crate::core::MoleculeEditor::new();
    let tombstone = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let center = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let carriers = ["F", "Cl", "Br", "I"]
        .into_iter()
        .map(element_atom)
        .map(|atom| molecule.add_atom(atom).expect("atom identifier capacity"))
        .collect::<Vec<_>>();
    for carrier in &carriers {
        molecule
            .add_bond(center, *carrier, BondOrder::Single)
            .expect("tetrahedral carrier bond");
    }
    molecule.delete_atom(tombstone).expect("isolated tombstone");
    mark_all_fresh(molecule.working_mut());
    assert_eq!(center, AtomId::new(1));
    assert_eq!(molecule.atom_count(), 5);

    let positions = test_positions(vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 0.0, -1.0),
    ]);
    let inferred = stereo_api::infer_coordinate_stereo(molecule.working(), &positions)
        .expect("dense positions should resolve through non-contiguous atom IDs");

    assert_eq!(inferred.elements.len(), 1);
    assert!(matches!(
        &inferred.elements[0].kind,
        StereoElementKind::Tetrahedral(stereo)
            if stereo.center == center
                && stereo.carriers
                    == carriers
                        .iter()
                        .copied()
                        .map(StereoCarrier::Atom)
                        .collect::<Vec<_>>()
    ));
}

#[test]
fn coordinate_stereo_inference_is_read_only_and_materializes_double_bond_stereo() {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left_carrier = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let right_carrier = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let double_bond = mol.add_bond(left, right, BondOrder::Double).expect("bond");
    mol.add_bond(left, left_carrier, BondOrder::Single)
        .expect("left carrier");
    mol.add_bond(right, right_carrier, BondOrder::Single)
        .expect("right carrier");
    let positions = test_positions(vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, -1.0, 0.0),
    ]);
    let before = mol.clone();

    let inferred = stereo_api::infer_coordinate_stereo(mol.working(), &positions)
        .expect("2D double-bond stereo should be inferred");
    assert_eq!(mol, before);
    assert_eq!(inferred.elements.len(), 1);
    let proposed = &inferred.elements[0];
    match &proposed.kind {
        StereoElementKind::DoubleBond(stereo) => {
            assert_eq!(stereo.bond, double_bond);
            assert_eq!(stereo.left, left);
            assert_eq!(stereo.right, right);
            assert_eq!(stereo.left_carrier, StereoCarrier::Atom(left_carrier));
            assert_eq!(stereo.right_carrier, StereoCarrier::Atom(right_carrier));
            assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Opposite));
        }
        other => panic!("expected double-bond stereo, found {other:?}"),
    }

    let report = stereo_api::materialize_coordinate_stereo(&mut mol, &positions)
        .expect("2D double-bond stereo should materialize");
    assert_eq!(report.created_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert_eq!(element, proposed);
}

#[test]
fn coordinate_stereo_inference_assigns_axis_only_when_requested() {
    let (mol, positions, axis) = coordinate_axis_graph(true);
    let before = mol.clone();

    let default = stereo_api::infer_coordinate_stereo(&mol, &positions)
        .expect("default coordinate-stereo inference should succeed");
    assert!(default.elements.is_empty());
    let inferred = stereo_api::infer_coordinate_stereo_with_options(
        &mol,
        &positions,
        CoordinateStereoOptions {
            infer_axes: true,
            ..Default::default()
        },
    )
    .expect("3D axis stereo should be inferred");
    assert_eq!(mol, before);
    assert_eq!(inferred.elements.len(), 1);
    let element = &inferred.elements[0];
    match &element.kind {
        StereoElementKind::Axis(stereo) => {
            assert_eq!(stereo.axis, axis);
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(2)),
                    StereoCarrier::Atom(AtomId::new(4)),
                ]
            );
            assert_eq!(stereo.orientation, Some(AxisOrientation::Clockwise));
        }
        other => panic!("expected axis stereo, found {other:?}"),
    }
}

#[test]
fn coordinate_stereo_inference_skips_axis_without_3d_handedness() {
    let (mol, positions, _axis) = coordinate_axis_graph(false);

    let result = stereo_api::infer_coordinate_stereo_with_options(
        &mol,
        &positions,
        CoordinateStereoOptions {
            infer_axes: true,
            ..Default::default()
        },
    )
    .expect("flat coordinates should be a successful non-assignment");
    assert!(result.elements.is_empty());
    assert!(mol.stereo_elements().next().is_none());
}

#[test]
fn coordinate_stereo_does_not_duplicate_existing_represented_stereo() {
    let (mut mol, center, carriers, _) = tetrahedral_marked_graph();
    let positions = test_positions(vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 0.0, -1.0),
    ]);
    mol.add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
        TetrahedralStereo {
            center,
            carriers: carriers.iter().copied().map(StereoCarrier::Atom).collect(),
            orientation: Some(TetrahedralOrientation::CounterClockwise),
        },
    )))
    .expect("represented source stereo");

    let inferred = stereo_api::infer_coordinate_stereo(mol.working(), &positions)
        .expect("represented stereo should make coordinate inference a no-op");
    assert!(inferred.elements.is_empty());
    let report = stereo_api::materialize_coordinate_stereo(&mut mol, &positions)
        .expect("represented stereo should make materialization a no-op");
    assert!(report.created_elements.is_empty());
    assert_eq!(mol.stereo_elements().count(), 1);
}

#[test]
fn coordinate_stereo_materialization_is_transactional_on_invalid_representation() {
    let (mut mol, center, carriers, _) = tetrahedral_marked_graph();
    insert_unchecked_stereo(
        &mut mol,
        StereoElement::new(StereoElementKind::Tetrahedral(TetrahedralStereo {
            center,
            carriers: vec![
                StereoCarrier::Atom(carriers[0]),
                StereoCarrier::Atom(carriers[0]),
                StereoCarrier::Atom(carriers[1]),
            ],
            orientation: Some(TetrahedralOrientation::Clockwise),
        })),
    );
    let before = mol.clone();

    let positions = crate::structure::Positions::zeros(mol.atom_count());
    let error = stereo_api::materialize_coordinate_stereo(&mut mol, &positions)
        .expect_err("invalid represented stereo must reject materialization");

    assert!(matches!(error, CoordinateStereoError::InvalidStereo(_)));
    assert_eq!(mol, before);
}

#[test]
fn coordinate_stereo_infers_three_explicit_ligands_and_preserves_handedness() {
    let mut implicit = read_smiles("C(F)(Cl)Br").unwrap();
    perceive(&mut implicit).unwrap();
    let mut explicit = implicit.clone();
    explicit.add_hydrogens().unwrap();
    for scale in [1.0e-100, 1.0, 1.0e100] {
        for sign in [-1.0, 1.0] {
            let points = vec![
                Point3::origin(),
                Point3::new(scale, 0.0, 0.0),
                Point3::new(0.0, scale, 0.0),
                Point3::new(0.0, 0.0, sign * scale),
            ];
            let inferred =
                stereo_api::infer_coordinate_stereo(&implicit, &test_positions(points.clone()))
                    .unwrap();
            assert_eq!(inferred.elements.len(), 1);
            let StereoElementKind::Tetrahedral(stereo) = &inferred.elements[0].kind else {
                unreachable!()
            };
            assert_eq!(stereo.carriers[3], StereoCarrier::ImplicitHydrogen);
            let expected = if sign > 0.0 {
                TetrahedralOrientation::Clockwise
            } else {
                TetrahedralOrientation::CounterClockwise
            };
            assert_eq!(stereo.orientation, Some(expected));

            let mut explicit_points = points;
            explicit_points.push(Point3::new(-scale, -scale, -sign * scale));
            let expanded =
                stereo_api::infer_coordinate_stereo(&explicit, &test_positions(explicit_points))
                    .unwrap();
            assert_eq!(expanded.elements.len(), 1);
            let StereoElementKind::Tetrahedral(stereo) = &expanded.elements[0].kind else {
                unreachable!()
            };
            assert_eq!(stereo.orientation, Some(expected));
        }
    }
    let flat = test_positions(vec![
        Point3::origin(),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, -1.0, 0.0),
    ]);
    assert!(stereo_api::infer_coordinate_stereo(&implicit, &flat)
        .unwrap()
        .elements
        .is_empty());
}

#[test]
fn coordinate_stereo_skips_overcoordinated_double_bond_endpoints() {
    let mut molecule = read_smiles("CP(C)(C)=NC").unwrap();
    perceive(&mut molecule).unwrap();
    assert!(stereo_api::detect_stereo_candidates(&molecule)
        .unwrap()
        .is_empty());
    let positions = test_positions(vec![
        Point3::new(0.0, 1.0, 0.0),
        Point3::origin(),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, -1.0, 0.0),
    ]);
    let mut editor = molecule.edit();
    stereo_api::materialize_coordinate_stereo(&mut editor, &positions).unwrap();
    assert_eq!(editor.finish().unwrap().stereo_elements().count(), 0);
}

#[test]
fn coordinate_stereo_infers_fully_substituted_alkene_at_any_coordinate_scale() {
    let molecule = read_smiles("FC(Cl)=C(Br)I").unwrap();
    for scale in [1.0e-100, 1.0, 1.0e100] {
        let points = [
            (0.0, 1.0),
            (0.0, 0.0),
            (0.0, -1.0),
            (1.0, 0.0),
            (1.0, -1.0),
            (1.0, 1.0),
        ]
        .map(|(x, y)| Point3::new(x * scale, y * scale, 0.0));
        let result =
            stereo_api::infer_coordinate_stereo(&molecule, &test_positions(points.to_vec()))
                .unwrap();
        assert_eq!(result.elements.len(), 1);
        let StereoElementKind::DoubleBond(stereo) = &result.elements[0].kind else {
            unreachable!()
        };
        assert_eq!(stereo.left_carrier, StereoCarrier::Atom(AtomId::new(0)));
        assert_eq!(stereo.right_carrier, StereoCarrier::Atom(AtomId::new(4)));
        assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Opposite));
    }
}

#[test]
fn coordinate_stereo_infers_a_tetrahedral_sulfoxide_lone_pair() {
    for input in ["S(=O)(C)CC", "[Se](=O)(C)CC", "[S+](C)(CC)CCC"] {
        let mut molecule = read_smiles(input).unwrap();
        perceive(&mut molecule).unwrap();
        let mut points = vec![
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        points.extend(
            (3..molecule.atom_count()).map(|index| Point3::new(0.0, 0.0, (index - 2) as f64)),
        );
        let result =
            stereo_api::infer_coordinate_stereo(&molecule, &test_positions(points)).unwrap();
        let stereo = result
            .elements
            .iter()
            .find_map(|element| match &element.kind {
                StereoElementKind::Tetrahedral(stereo) if stereo.center == AtomId::new(0) => {
                    Some(stereo)
                }
                _ => None,
            })
            .expect("three-coordinate S/Se stereo");
        assert_eq!(stereo.carriers[3], StereoCarrier::ImplicitLonePair);
        assert_eq!(stereo.orientation, Some(TetrahedralOrientation::Clockwise));
    }
}

#[test]
fn coordinate_axis_inference_respects_pyramidal_lone_pair_endpoints() {
    for source in ["S(=O)(C)C1=CC=CC=C1", "[Se](=O)(C)C1=CC=CC=C1"] {
        let mut molecule = read_smiles(source).unwrap();
        perceive(&mut molecule).unwrap();
        let points = (0..molecule.atom_count())
            .map(|index| {
                let t = index as f64;
                Point3::new(t, (t * 1.7).sin(), (t * 2.3).cos())
            })
            .collect();
        let result = stereo_api::infer_coordinate_stereo_with_options(
            &molecule,
            &test_positions(points),
            CoordinateStereoOptions {
                infer_axes: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(result.elements.iter().any(|element| matches!(
            &element.kind,
            StereoElementKind::Tetrahedral(stereo) if stereo.center == AtomId::new(0)
        )));
        assert!(
            result
                .elements
                .iter()
                .all(|element| !matches!(element.kind, StereoElementKind::Axis(_))),
            "{source}: {:?}",
            result.elements
        );
    }
}

#[test]
fn coordinate_axis_inference_does_not_treat_saturated_rings_as_sp2() {
    let molecule = read_smiles("C1CCCCC1C2CCCCC2").unwrap();
    let points = (0..molecule.atom_count())
        .map(|index| {
            let t = index as f64;
            Point3::new(t, (t * 1.7).sin(), (t * 2.3).cos())
        })
        .collect();
    let result = stereo_api::infer_coordinate_stereo_with_options(
        &molecule,
        &test_positions(points),
        CoordinateStereoOptions {
            infer_axes: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(result
        .elements
        .iter()
        .all(|element| !matches!(element.kind, StereoElementKind::Axis(_))));
}

#[test]
fn coordinate_axis_handedness_is_scale_and_rotation_invariant() {
    let (molecule, positions, _) = coordinate_axis_graph(true);
    for scale in [1.0e-100, 1.0, 1.0e100] {
        let points = (0..positions.len())
            .map(|index| {
                let point = positions.position_at(index).unwrap().into_value();
                // A cyclic permutation is a proper rotation.
                Point3::new(point.y * scale, point.z * scale, point.x * scale)
            })
            .collect();
        let result = stereo_api::infer_coordinate_stereo_with_options(
            &molecule,
            &test_positions(points),
            CoordinateStereoOptions {
                infer_axes: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(result.elements.len(), 1);
        assert!(
            matches!(&result.elements[0].kind, StereoElementKind::Axis(stereo) if stereo.orientation == Some(AxisOrientation::Clockwise))
        );
    }
}
