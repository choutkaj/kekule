//! Stereo read from SMILES direction marks, Molfile wedges, 3D tetrahedra,
//! atropisomeric axes and V3000 stereo groups.

use super::*;

#[test]
fn interpretation_assembles_paired_directional_marks_into_double_bond_element() {
    let (molecule, report) =
        read_smiles_with_report("C/C=C\\F").expect("directional smiles should interpret");
    assert_eq!(report.created_stereo_elements().len(), 1);
    assert!(molecule.stereo_elements().next().is_some());
    let element = molecule
        .stereo_element(report.created_stereo_elements()[0])
        .expect("created stereo element");
    match &element.kind {
        StereoElementKind::DoubleBond(stereo) => {
            assert_eq!(stereo.bond, BondId::new(1));
            assert_eq!(stereo.left, AtomId::new(1));
            assert_eq!(stereo.right, AtomId::new(2));
            assert_eq!(stereo.left_carrier, StereoCarrier::Atom(AtomId::new(0)));
            assert_eq!(stereo.right_carrier, StereoCarrier::Atom(AtomId::new(3)));
            assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Together));
        }
        other => panic!("expected double-bond stereo, found {other:?}"),
    }
}

#[test]
fn equivalent_smiles_direction_tokens_publish_equivalent_canonical_stereo() {
    for (first, second) in [("C/C=C/C", r"C\C=C\C"), (r"C/C=C\C", r"C\C=C/C")] {
        let first = read_smiles(first).expect("first directional spelling should interpret");
        let second = read_smiles(second).expect("second directional spelling should interpret");
        let first = first
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>();
        let second = second
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>();

        assert_eq!(first, second);
        assert_eq!(first.len(), 1);
    }
}

#[test]
fn alternate_directional_source_carriers_publish_identical_double_bond_stereo() {
    let source_graph = || {
        let mut molecule = crate::core::MoleculeEditor::new();
        let left = molecule
            .add_atom(carbon())
            .expect("atom identifier capacity");
        let right = molecule
            .add_atom(carbon())
            .expect("atom identifier capacity");
        let left_reference = molecule
            .add_atom(element_atom("F"))
            .expect("atom identifier capacity");
        let left_alternative = molecule
            .add_atom(element_atom("Cl"))
            .expect("atom identifier capacity");
        let right_reference = molecule
            .add_atom(element_atom("Br"))
            .expect("atom identifier capacity");
        let right_alternative = molecule
            .add_atom(element_atom("I"))
            .expect("atom identifier capacity");
        molecule
            .add_bond(left, right, BondOrder::Double)
            .expect("double bond");
        let left_reference_bond = molecule
            .add_bond(left, left_reference, BondOrder::Single)
            .expect("left reference bond");
        let left_alternative_bond = molecule
            .add_bond(left, left_alternative, BondOrder::Single)
            .expect("left alternative bond");
        let right_reference_bond = molecule
            .add_bond(right, right_reference, BondOrder::Single)
            .expect("right reference bond");
        let right_alternative_bond = molecule
            .add_bond(right, right_alternative, BondOrder::Single)
            .expect("right alternative bond");
        (
            molecule,
            left,
            right,
            left_reference_bond,
            left_alternative_bond,
            right_reference_bond,
            right_alternative_bond,
        )
    };

    let canonicalize = |mut molecule: MoleculeEditor, marks: &[SourceStereoBondMark]| {
        let report = canonicalize_molecule_for_publication(molecule.working_mut(), None, marks)
            .expect("paired directional source marks should canonicalize");
        assert_eq!(report.created_stereo_elements.len(), 1);
        molecule
            .stereo_element(report.created_stereo_elements[0])
            .expect("created double-bond element")
            .clone()
    };

    let (molecule, left, right, left_reference, _, right_reference, _) = source_graph();
    let expected = canonicalize(
        molecule,
        &[
            SourceStereoBondMark {
                bond: left_reference,
                from: left,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
            SourceStereoBondMark {
                bond: right_reference,
                from: right,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
        ],
    );

    let (molecule, left, right, _, left_alternative, right_reference, _) = source_graph();
    let alternate_left = canonicalize(
        molecule,
        &[
            SourceStereoBondMark {
                bond: left_alternative,
                from: left,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
            SourceStereoBondMark {
                bond: right_reference,
                from: right,
                kind: SourceStereoBondMarkKind::DirectionalDown,
            },
        ],
    );

    let (molecule, left, right, _, left_alternative, _, right_alternative) = source_graph();
    let both_alternatives = canonicalize(
        molecule,
        &[
            SourceStereoBondMark {
                bond: left_alternative,
                from: left,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
            SourceStereoBondMark {
                bond: right_alternative,
                from: right,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
        ],
    );

    assert_eq!(alternate_left, expected);
    assert_eq!(both_alternatives, expected);
}

#[test]
fn smiles_ring_direction_preserves_the_textual_origin_endpoint() {
    let marked_when_opened =
        read_smiles(r"F/C=C/1CCCCC1").expect("opening ring direction should interpret");
    let marked_when_closed =
        read_smiles(r"F/C=C1CCCCC\1").expect("closing ring direction should interpret");
    let opened_stereo = marked_when_opened
        .stereo_elements()
        .map(|(_, element)| element.clone())
        .collect::<Vec<_>>();
    let closed_stereo = marked_when_closed
        .stereo_elements()
        .map(|(_, element)| element.clone())
        .collect::<Vec<_>>();

    assert_eq!(opened_stereo.len(), 1);
    assert_eq!(closed_stereo, opened_stereo);
}

#[test]
fn interpretation_enforces_small_ring_double_bond_boundary() {
    let document = smiles_api::parse_str(r"C1/C=C\CCC1").expect("source syntax should parse");
    let error = smiles_api::interpret(&document)
        .expect_err("excluded small-ring directional marks must reject interpretation");
    assert_eq!(error.offset(), 2);
    assert!(error.message().contains("UnpairedDirectionalBondMark"));

    let (cyclooctene, report) =
        read_smiles_with_report(r"C1/C=C\CCCCC1").expect("marked cyclooctene interprets");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = cyclooctene
        .stereo_element(report.created_stereo_elements()[0])
        .expect("created stereo element");
    assert!(matches!(element.kind, StereoElementKind::DoubleBond(_)));
}

#[test]
fn interpretation_assembles_molfile_wedge_into_tetrahedral_element() {
    let input = "\
wedge
kekule

  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
    1.0000    0.0000    0.0000 F   0  0  0  0  0  0
   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0
    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 I   0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
";
    let (molecule, report) =
        read_molfile_with_report(input).expect("wedge molfile should interpret");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = molecule
        .stereo_element(report.created_stereo_elements()[0])
        .expect("created stereo element");
    assert!(element.is_specified());
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assert_eq!(stereo.center, AtomId::new(0));
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(1)),
                    StereoCarrier::Atom(AtomId::new(2)),
                    StereoCarrier::Atom(AtomId::new(3)),
                    StereoCarrier::Atom(AtomId::new(4)),
                ]
            );
            assert_eq!(
                stereo.orientation,
                Some(TetrahedralOrientation::CounterClockwise)
            );
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn canonical_tetrahedral_stereo_is_identical_across_smiles_molfile_and_manual_sources() {
    let smiles = read_smiles("F[C@](Cl)(Br)I").expect("tetrahedral SMILES should interpret");
    let expected = smiles
        .stereo_elements()
        .next()
        .expect("SMILES should create canonical tetrahedral stereo")
        .1
        .clone();

    let model = tetrahedral_drawing(&smiles);
    for written in [
        molfile::write(
            &model,
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V2000,
            },
        )
        .expect("canonical stereo should project to V2000"),
        molfile::write(
            &model,
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V3000,
            },
        )
        .expect("canonical stereo should project to V3000"),
    ] {
        let interpreted = read_molfile(&written).expect("projected Molfile should interpret");
        let actual = interpreted
            .stereo_elements()
            .next()
            .expect("Molfile should recreate canonical tetrahedral stereo")
            .1;
        assert_eq!(actual, &expected);
    }

    let mut manual = crate::core::MoleculeEditor::new();
    let fluorine = manual
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let center = manual.add_atom(carbon()).expect("atom identifier capacity");
    let chlorine = manual
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let bromine = manual
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let iodine = manual
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");
    for carrier in [fluorine, chlorine, bromine, iodine] {
        manual
            .add_bond(center, carrier, BondOrder::Single)
            .expect("tetrahedral carrier bond");
    }
    let manual_id = manual
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(fluorine),
                    StereoCarrier::Atom(chlorine),
                    StereoCarrier::Atom(bromine),
                    StereoCarrier::Atom(iodine),
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("manual canonical stereo element");
    assert_eq!(manual.stereo_element(manual_id).unwrap(), &expected);
}

#[test]
fn interpretation_uses_source_declared_h_for_molfile_wedge_geometry() {
    let (molecule, report) = read_molfile_with_report(implicit_h_wedge_geometry_molblock())
        .expect("implicit-H wedge molfile should interpret");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = molecule
        .stereo_element(report.created_stereo_elements()[0])
        .expect("created stereo element");
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assert_eq!(stereo.center, AtomId::new(0));
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(1)),
                    StereoCarrier::Atom(AtomId::new(2)),
                    StereoCarrier::Atom(AtomId::new(3)),
                    StereoCarrier::ImplicitHydrogen,
                ]
            );
            assert_eq!(
                stereo.orientation,
                Some(TetrahedralOrientation::CounterClockwise)
            );
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn normalization_assembles_wedge_either_as_explicit_unknown() {
    let (mut mol, center, carriers, marked_bond) = tetrahedral_marked_graph();
    let source_stereo = [SourceStereoBondMark {
        bond: marked_bond,
        from: center,
        kind: SourceStereoBondMarkKind::WedgeEither,
    }];

    let report = canonicalize_molecule_for_publication(mol.working_mut(), None, &source_stereo)
        .expect("wedge/either should assemble as unknown stereo");
    assert_eq!(report.created_stereo_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_stereo_elements[0])
        .expect("created stereo element");
    assert!(element.is_explicitly_unknown());
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assert_eq!(stereo.center, center);
            assert_eq!(stereo.carriers[0], StereoCarrier::Atom(carriers[0]));
            assert_eq!(stereo.orientation, None);
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn alternate_tetrahedral_wedge_carriers_publish_identical_canonical_stereo() {
    struct Drawing([Point3; 5]);
    impl AtomPositionSource for Drawing {
        fn position_value(&self, atom: AtomId) -> Option<Point3> {
            self.0.get(atom.index()).copied()
        }
    }
    let canonicalize = |kind, bond| {
        let (mut molecule, center, _, _) = tetrahedral_marked_graph();
        let points = Drawing([
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
        ]);
        let report = canonicalize_molecule_for_publication(
            molecule.working_mut(),
            Some(&points),
            &[SourceStereoBondMark {
                bond,
                from: center,
                kind,
            }],
        )
        .expect("tetrahedral source wedge should canonicalize");
        assert_eq!(report.created_stereo_elements.len(), 1);
        molecule
            .stereo_element(report.created_stereo_elements[0])
            .expect("created tetrahedral element")
            .clone()
    };

    let wedge_on_first = canonicalize(SourceStereoBondMarkKind::WedgeUp, BondId::new(0));
    // In this cross-shaped drawing, opposite bonds use the same wedge direction
    // to encode the same configuration; carrier-list parity alone is insufficient.
    let wedge_on_second = canonicalize(SourceStereoBondMarkKind::WedgeUp, BondId::new(1));
    assert_eq!(wedge_on_second, wedge_on_first);

    let unknown_on_first = canonicalize(SourceStereoBondMarkKind::WedgeEither, BondId::new(0));
    let unknown_on_second = canonicalize(SourceStereoBondMarkKind::WedgeEither, BondId::new(1));
    assert_eq!(unknown_on_second, unknown_on_first);
    assert!(unknown_on_first.is_explicitly_unknown());
}

#[test]
fn source_stereo_rejects_an_origin_outside_the_marked_bond() {
    let mut molecule = crate::core::MoleculeEditor::new();
    let a = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let b = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let outside = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let bond = molecule
        .add_bond(a, b, BondOrder::Single)
        .expect("marked bond");
    let source_stereo = [SourceStereoBondMark {
        bond,
        from: outside,
        kind: SourceStereoBondMarkKind::WedgeUp,
    }];

    let error = canonicalize_molecule_for_publication(molecule.working_mut(), None, &source_stereo)
        .expect_err("the marked origin must be an endpoint of its bond");

    assert!(matches!(
        error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues == vec![SourceStereoNormalizationIssue::InvalidSourceBondMarkEndpoint {
                bond,
                from: outside,
            }]
    ));
}

#[test]
fn normalization_reports_ambiguous_tetrahedral_wedge_marks() {
    let (mut mol, center, _carriers, first_bond) = tetrahedral_marked_graph();
    let second_bond = BondId::new(1);
    let source_stereo = [
        SourceStereoBondMark {
            bond: first_bond,
            from: center,
            kind: SourceStereoBondMarkKind::WedgeUp,
        },
        SourceStereoBondMark {
            bond: second_bond,
            from: center,
            kind: SourceStereoBondMarkKind::WedgeDown,
        },
    ];

    let report = canonicalize_molecule_for_publication(mol.working_mut(), None, &source_stereo)
        .expect("ambiguous wedges should warn without failing");

    assert!(report
        .warnings
        .contains(&NormalizationWarning::AmbiguousTetrahedralWedgeMarks {
            center,
            mark_count: 2,
        }));
    assert!(report.created_stereo_elements.is_empty());
    assert!(mol.stereo_elements().next().is_none());
}

#[test]
fn a_wedge_without_coordinates_does_not_define_atom_order_parity() {
    for kind in [
        SourceStereoBondMarkKind::WedgeUp,
        SourceStereoBondMarkKind::WedgeDown,
    ] {
        let (mut molecule, center, _, bond) = tetrahedral_marked_graph();
        let report = canonicalize_molecule_for_publication(
            molecule.working_mut(),
            None,
            &[SourceStereoBondMark {
                bond,
                from: center,
                kind,
            }],
        )
        .unwrap();
        assert!(report.created_stereo_elements.is_empty());
        assert_eq!(
            report.warnings,
            vec![NormalizationWarning::AmbiguousTetrahedralWedgeMarks {
                center,
                mark_count: 1
            }]
        );
        assert!(molecule.stereo_elements().next().is_none());
    }
}

#[test]
fn invalid_source_stereo_reports_an_issue_without_publishing_a_placeholder_element() {
    let mut marked = crate::core::MoleculeEditor::new();
    let a = marked.add_atom(carbon()).expect("atom identifier capacity");
    let b = marked.add_atom(carbon()).expect("atom identifier capacity");
    let bond = marked.add_bond(a, b, BondOrder::Single).expect("bond");
    let marked_source = [SourceStereoBondMark {
        bond,
        from: a,
        kind: SourceStereoBondMarkKind::WedgeEither,
    }];

    let positions = crate::structure::Positions::zeros(marked.atom_count());
    let coordinate_result = stereo_api::infer_coordinate_stereo(marked.working(), &positions)
        .expect("coordinate inference is independent of detached source marks");
    assert!(coordinate_result.elements.is_empty());
    let marked_before = marked.clone();

    let marked_error =
        canonicalize_molecule_for_publication(marked.working_mut(), None, &marked_source)
            .expect_err("unassembled tetrahedral mark should fail");
    assert!(matches!(
        marked_error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnassembledTetrahedralBondMark {
            bond,
            kind: SourceStereoBondMarkKind::WedgeEither,
        })
    ));
    assert!(marked.stereo_elements().next().is_none());
    assert_eq!(marked, marked_before);

    let mut unsupported = crate::core::MoleculeEditor::new();
    let c = unsupported
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let d = unsupported
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let double_bond = unsupported.add_bond(c, d, BondOrder::Double).expect("bond");
    let unsupported_source = [SourceStereoBondMark {
        bond: double_bond,
        from: c,
        kind: SourceStereoBondMarkKind::DoubleBondEither,
    }];
    let unsupported_error =
        canonicalize_molecule_for_publication(unsupported.working_mut(), None, &unsupported_source)
            .expect_err("unsupported double-bond mark should fail");
    assert!(matches!(
        unsupported_error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnsupportedSourceBondMark {
            bond: double_bond,
            kind: SourceStereoBondMarkKind::DoubleBondEither,
        })
    ));

    let mut unknown = crate::core::MoleculeEditor::new();
    let left = unknown
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let right = unknown
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let left_carrier = unknown
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let right_carrier = unknown
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let unknown_bond = unknown
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    unknown
        .add_bond(left, left_carrier, BondOrder::Single)
        .expect("left carrier");
    unknown
        .add_bond(right, right_carrier, BondOrder::Single)
        .expect("right carrier");
    let unknown_source = [SourceStereoBondMark {
        bond: unknown_bond,
        from: left,
        kind: SourceStereoBondMarkKind::DoubleBondEither,
    }];

    let unknown_report =
        canonicalize_molecule_for_publication(unknown.working_mut(), None, &unknown_source)
            .expect("double-bond either should assemble as unknown stereo");
    assert_eq!(unknown_report.created_stereo_elements.len(), 1);
    let (_, element) = unknown.stereo_elements().next().expect("unknown element");
    assert!(matches!(
        &element.kind,
        StereoElementKind::DoubleBond(stereo)
            if stereo.bond == unknown_bond && stereo.orientation.is_none()
    ));

    let mut absent = crate::core::MoleculeEditor::new();
    let x = absent.add_atom(carbon()).expect("atom identifier capacity");
    let y = absent.add_atom(carbon()).expect("atom identifier capacity");
    absent.add_bond(x, y, BondOrder::Single).expect("bond");
    let absent_report = canonicalize_molecule_for_publication(absent.working_mut(), None, &[])
        .expect("unmarked molecule should normalize");
    assert!(absent_report.created_stereo_elements.is_empty());
    assert!(absent.stereo_elements().next().is_none());
}

#[test]
fn failed_source_stereo_canonicalization_reports_the_unpaired_mark() {
    let mut molecule = read_smiles("F[C@](Cl)(Br)I").expect("stereo SMILES should parse");
    perceive(&mut molecule).expect("stored stereo should prepare");
    let marked_bond = molecule.bond_ids().next().expect("single bond");
    let source_stereo = [SourceStereoBondMark {
        bond: marked_bond,
        from: molecule.bond(marked_bond).expect("marked bond").a(),
        kind: SourceStereoBondMarkKind::DirectionalUp,
    }];
    let cip =
        stereo_api::assign_cip_descriptors(&mut molecule).expect("CIP assignment should succeed");
    assert_eq!(cip.assigned.len(), 1);
    stereo_api::validate_stereo(&molecule).expect("stored-state validation should succeed");
    let error = molecule
        .canonicalize_fixture_with_source_stereo(&source_stereo)
        .expect_err("unpaired directional mark should fail canonicalization");

    assert!(matches!(
        error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnpairedDirectionalBondMark {
                bond: marked_bond,
            })
    ));
}

#[test]
fn interpretation_assembles_molfile_atropisomeric_axis() {
    let (molecule, report) = read_molfile_with_report(rdkit_rp6306_atrop_molblock())
        .expect("RDKit atropisomer fixture interprets");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = molecule
        .stereo_element(report.created_stereo_elements()[0])
        .expect("created axis element");
    match &element.kind {
        StereoElementKind::Axis(stereo) => {
            assert_eq!(stereo.axis, BondId::new(3));
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(6)),
                    StereoCarrier::Atom(AtomId::new(11)),
                ]
            );
            assert_eq!(stereo.orientation, Some(AxisOrientation::Clockwise));
        }
        other => panic!("expected axis stereo, found {other:?}"),
    }
}

#[test]
fn molfile_writers_project_tetrahedral_stereo_independent_of_bond_endpoint_storage() {
    let molecule = canonical_tetrahedral_molecule();
    let reversed = reverse_bond_endpoint_storage(&molecule);
    let expected = molecule
        .stereo_elements()
        .next()
        .expect("canonical tetrahedral element")
        .1
        .clone();

    for molecule in [&molecule, &reversed] {
        let before = molecule.clone();
        let model = tetrahedral_drawing(molecule);
        for written in [
            molfile::write(
                &model,
                molfile::MolfileWriteOptions {
                    version: molfile::MolfileWriteVersion::V2000,
                },
            )
            .expect("V2000 projects tetrahedral stereo"),
            molfile::write(
                &model,
                molfile::MolfileWriteOptions {
                    version: molfile::MolfileWriteVersion::V3000,
                },
            )
            .expect("V3000 projects tetrahedral stereo"),
        ] {
            let (reparsed, report) =
                read_molfile_with_report(&written).expect("projected tetrahedral stereo reparses");
            assert_eq!(report.created_stereo_elements().len(), 1);
            let actual = reparsed
                .stereo_element(report.created_stereo_elements()[0])
                .expect("reparsed tetrahedral element");
            assert_eq!(actual.kind, expected.kind);
        }
        assert_eq!(*molecule, before);
    }
}

#[test]
fn geometry_free_molfile_writers_reject_axis_stereo_without_mutating_the_molecule() {
    let (molecule, report) = read_molfile_with_report(rdkit_rp6306_atrop_molblock())
        .expect("RDKit atropisomer fixture interprets");
    let expected = molecule
        .stereo_element(report.created_stereo_elements()[0])
        .expect("canonical axis element");
    let axis = match &expected.kind {
        StereoElementKind::Axis(stereo) => stereo.axis,
        other => panic!("expected canonical axis stereo, found {other:?}"),
    };
    let reversed = reverse_bond_endpoint_storage_except(&molecule, &[axis]);

    for molecule in [&molecule, &reversed] {
        let before = molecule.clone();
        assert!(molfile::write(
            molecule,
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V2000
            }
        )
        .is_err());
        assert!(molfile::write(
            molecule,
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V3000
            }
        )
        .is_err());
        assert_eq!(*molecule, before);
    }
}

#[test]
fn interpretation_prefers_exocyclic_molfile_atropisomeric_axis() {
    let (molecule, report) = read_molfile_with_report(rdkit_rp6306_atrop3_molblock())
        .expect("RDKit alternate atropisomer fixture interprets");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = molecule
        .stereo_element(report.created_stereo_elements()[0])
        .expect("created axis element");
    match &element.kind {
        StereoElementKind::Axis(stereo) => {
            assert_eq!(stereo.axis, BondId::new(3));
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(6)),
                    StereoCarrier::Atom(AtomId::new(11)),
                ]
            );
            assert_eq!(stereo.orientation, Some(AxisOrientation::Clockwise));
        }
        other => panic!("expected axis stereo, found {other:?}"),
    }
}

#[test]
fn interpretation_materializes_omitted_hydrogen_alongside_atrop_stereo() {
    let (mut molecule, report) = read_molfile_with_report(rdkit_bms986142_atrop5_molblock())
        .expect("source valence defines the omitted hydrogen");
    assert_eq!(report.created_stereo_elements().len(), 2);
    assert_eq!(
        molecule.atom(AtomId::new(10)).unwrap().hydrogens,
        ImplicitHydrogens::Fixed(1)
    );
    assert!(molecule.stereo_elements().any(|(_, element)| matches!(&element.kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(8))));
    perceive(&mut molecule).unwrap();
    let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
    assert_eq!(
        assigned
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
        vec![StereoDescriptor::S, StereoDescriptor::P]
    );
}

#[test]
fn interpretation_preserves_one_ring_endpoint_atrop_and_omitted_hydrogen() {
    for fixture in [
        rdkit_zm374979_atrop1_molblock(),
        rdkit_zm374979_atrop2_molblock(),
    ] {
        let (mut molecule, report) =
            read_molfile_with_report(fixture).expect("omitted source hydrogen interprets");
        assert_eq!(report.created_stereo_elements().len(), 2);
        assert_eq!(
            molecule.atom(AtomId::new(3)).unwrap().hydrogens,
            ImplicitHydrogens::Fixed(1)
        );
        assert!(molecule.stereo_elements().any(|(_, element)| matches!(&element.kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(33))));
        perceive(&mut molecule).unwrap();
        let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
        assert_eq!(
            assigned
                .assigned
                .iter()
                .map(|assignment| assignment.descriptor)
                .collect::<Vec<_>>(),
            vec![StereoDescriptor::R, StereoDescriptor::M]
        );
    }
}

#[test]
fn interpretation_assembles_ring_internal_molfile_atrop_axis() {
    for fixture in [
        rdkit_macrocycle8_ortho_wedge_molblock(),
        rdkit_macrocycle8_ortho_hash_molblock(),
    ] {
        let (molecule, report) = read_molfile_with_report(fixture)
            .expect("RDKit macrocyclic atropisomer fixture interprets");
        assert_eq!(report.created_stereo_elements().len(), 1);
        assert!(molecule.stereo_elements().any(|(_, element)| {
            matches!(&element.kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(15))
        }));
    }
}

fn tetrahedral_drawing(molecule: &Molecule) -> Model {
    let points = molecule
        .atoms()
        .map(|(_, atom)| match atom.element.symbol() {
            "C" => Point3::origin(),
            "F" => Point3::new(1.0, 0.0, 0.0),
            "Cl" => Point3::new(-1.0, 0.0, 0.0),
            "Br" => Point3::new(0.0, 1.0, 0.0),
            "I" => Point3::new(0.0, -1.0, 0.0),
            _ => unreachable!("tetrahalomethane regression"),
        })
        .collect();
    Model::from_molecule(molecule.clone(), &test_positions(points)).unwrap()
}

fn canonical_tetrahedral_molecule() -> Molecule {
    let (mut molecule, center, carriers, _) = tetrahedral_marked_graph();
    molecule
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: carriers.into_iter().map(StereoCarrier::Atom).collect(),
                orientation: Some(TetrahedralOrientation::CounterClockwise),
            },
        )))
        .expect("canonical tetrahedral element");
    molecule.finish().expect("tetrahedral molecule publishes")
}

fn reverse_bond_endpoint_storage(molecule: &Molecule) -> Molecule {
    reverse_bond_endpoint_storage_except(molecule, &[])
}

fn reverse_bond_endpoint_storage_except(molecule: &Molecule, excluded: &[BondId]) -> Molecule {
    let mut reversed = molecule.clone();
    for (index, bond) in reversed.graph.bonds.iter_mut().enumerate() {
        if excluded.contains(&BondId::new(index as u32)) {
            continue;
        }
        let Some(bond) = bond else {
            continue;
        };
        std::mem::swap(&mut bond.a, &mut bond.b);
    }
    reversed
}

fn coordinate_tetrahedron(v3000: bool, points: &[[f64; 3]], mark: u8) -> String {
    let symbols = ["C", "F", "Cl", "Br", "I"];
    let mut source = if v3000 {
        format!("coordinate stereo\nkekule\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {} {} 0 0 0\nM  V30 BEGIN ATOM\n", points.len(), points.len() - 1)
    } else {
        format!(
            "coordinate stereo\nkekule\n\n{:3}{:3}  0  0  0  0            999 V2000\n",
            points.len(),
            points.len() - 1
        )
    };
    for (index, [x, y, z]) in points.iter().enumerate() {
        if v3000 {
            source += &format!("M  V30 {} {} {x} {y} {z} 0\n", index + 1, symbols[index]);
        } else {
            source += &format!(
                "{x:10.4}{y:10.4}{z:10.4} {:<3} 0  0  0  0  0  0\n",
                symbols[index]
            );
        }
    }
    if v3000 {
        source += "M  V30 END ATOM\nM  V30 BEGIN BOND\n";
    }
    for index in 1..points.len() {
        let mark = if index == 1 { mark } else { 0 };
        if v3000 {
            source += &format!(
                "M  V30 {index} 1 1 {}{}\n",
                index + 1,
                if mark == 0 {
                    String::new()
                } else {
                    format!(" CFG={mark}")
                }
            );
        } else {
            source += &format!("  1{:3}  1{mark:3}  0  0  0\n", index + 1);
        }
    }
    if v3000 {
        source += "M  V30 END BOND\nM  V30 END CTAB\n";
    }
    source + "M  END\n"
}

#[test]
fn molfile_unmarked_3d_tetrahedra_use_shared_geometry_without_installing_perception() {
    for v3000 in [false, true] {
        for explicit_count in [3, 4] {
            let mut descriptors = Vec::new();
            for reflected in [false, true] {
                let sign = if reflected { -1.0 } else { 1.0 };
                let points = [
                    [0.0, 0.0, 0.0],
                    [sign, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [-sign, -1.0, -1.0],
                ];
                let source = coordinate_tetrahedron(v3000, &points[..=explicit_count], 0);
                let (mut molecule, report) = read_molfile_with_report(&source).unwrap();
                assert_eq!(report.created_stereo_elements().len(), 1);
                assert!(!molecule.perception().has_valence());
                assert!(report.warnings().is_empty());
                if explicit_count == 3 {
                    assert_eq!(
                        molecule.atom(AtomId::new(0)).unwrap().hydrogens,
                        ImplicitHydrogens::Fixed(1)
                    );
                }
                perceive(&mut molecule).unwrap();
                let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
                assert_eq!(assigned.assigned.len(), 1);
                let descriptor = assigned.assigned[0].descriptor;
                // Independently checked with RDKit 2026.03.3 from these CTAB
                // coordinates, including the omitted fourth hydrogen.
                let expected = if (explicit_count == 3) != reflected {
                    StereoDescriptor::R
                } else {
                    StereoDescriptor::S
                };
                assert_eq!(descriptor, expected);
                descriptors.push(descriptor);
            }
            assert_ne!(descriptors[0], descriptors[1]);
        }
    }
}

#[test]
fn molfile_3d_tetrahedral_inference_is_scale_rotation_and_translation_invariant() {
    for scale in [1.0e-100, 1.0, 1.0e100] {
        for rotated in [false, true] {
            let points = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
            .map(|[x, y, z]| {
                let [x, y, z] = if rotated { [z, x, y] } else { [x, y, z] };
                [scale * (x + 2.0), scale * (y - 3.0), scale * (z + 4.0)]
            });
            let source = coordinate_tetrahedron(true, &points, 0);
            let (mut molecule, _) = read_molfile_with_report(&source).unwrap();
            perceive(&mut molecule).unwrap();
            let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
            assert_eq!(assigned.assigned.len(), 1);
            assert_eq!(assigned.assigned[0].descriptor, StereoDescriptor::R);
        }
    }
}

#[test]
fn molfile_3d_inference_preserves_unknown_and_rejects_degenerate_or_symmetric_sites() {
    for v3000 in [false, true] {
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let unknown = coordinate_tetrahedron(v3000, &points, if v3000 { 2 } else { 4 });
        let (molecule, _) = read_molfile_with_report(&unknown).unwrap();
        assert_eq!(molecule.stereo_elements().count(), 1);
        assert!(molecule
            .stereo_elements()
            .all(|(_, element)| !element.is_specified()));
        let source = coordinate_tetrahedron(v3000, &points, 0);
        let symmetric = source.replace("Cl", "F ");
        assert!(read_molfile_with_report(&symmetric)
            .unwrap()
            .0
            .stereo_elements()
            .next()
            .is_none());
        let flat = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [-1.0, -1.0, -2.0],
        ];
        let source = coordinate_tetrahedron(v3000, &flat, 0);
        assert!(read_molfile_with_report(&source)
            .unwrap()
            .0
            .stereo_elements()
            .next()
            .is_none());
    }
}

#[test]
fn molfile_writers_preserve_unasserted_3d_tetrahedra_as_unknown() {
    for smiles in ["C(F)(Cl)Br", "C(F)(Cl)(Br)I"] {
        let molecule = read_smiles(smiles).unwrap();
        let points = [
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(-1.0, -1.0, -1.0),
        ];
        let positions = test_positions(points[..molecule.atom_count()].to_vec());
        let model = Model::from_molecule(molecule.clone(), &positions).unwrap();
        let before = model.clone();
        for source in [
            molfile::write(
                &model,
                molfile::MolfileWriteOptions {
                    version: molfile::MolfileWriteVersion::V2000,
                },
            )
            .unwrap(),
            molfile::write(
                &model,
                molfile::MolfileWriteOptions {
                    version: molfile::MolfileWriteVersion::V3000,
                },
            )
            .unwrap(),
        ] {
            let (parsed, _) = read_molfile_with_report(&source).unwrap();
            assert_eq!(parsed.stereo_elements().count(), 1);
            assert!(parsed
                .stereo_elements()
                .all(|(_, element)| !element.is_specified()));
        }
        assert_eq!(model, before);
    }
}

fn wedged_tetrahedron(symbol: &str, charge_code: u8, mirror: bool) -> String {
    let y = if mirror { -1.0 } else { 1.0 };
    format!(
        "stereo drawing\nkekule\n\n  4  3  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 {symbol:<3} 0  {charge_code}  0  0  0  0\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0\n    0.0000{y:10.4}    0.0000 Br  0  0  0  0  0  0\n  1  2  1  1  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\nM  END\n"
    )
}

fn wedge_drawing(points: &[[f64; 2]], scale: f64, mirror: f64, v3000: bool) -> String {
    let symbols = ["C", "F", "Cl", "Br", "I"];
    let mut source = format!(
        "drawing\nkekule\n\n{:3}{:3}  0  0  0  0            999 V2000\n",
        points.len(),
        points.len() - 1
    );
    if v3000 {
        source = format!("drawing\nkekule\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {} {} 0 0 0\nM  V30 BEGIN ATOM\n", points.len(), points.len() - 1);
    }
    for (index, [x, y]) in points.iter().enumerate() {
        let (x, y) = (scale * (x + 3.0), scale * (mirror * y - 2.0));
        if v3000 {
            source += &format!("M  V30 {} {} {x} {y} 0 0\n", index + 1, symbols[index]);
        } else {
            source += &format!(
                "{x:10.4}{y:10.4}    0.0000 {:<3} 0  0  0  0  0  0\n",
                symbols[index]
            );
        }
    }
    if v3000 {
        source += "M  V30 END ATOM\nM  V30 BEGIN BOND\n";
    }
    for index in 1..points.len() {
        if v3000 {
            source += &format!(
                "M  V30 {index} 1 1 {}{}\n",
                index + 1,
                if index == 1 { " CFG=1" } else { "" }
            );
        } else {
            source += &format!(
                "  1{:3}  1{:3}  0  0  0\n",
                index + 1,
                usize::from(index == 1)
            );
        }
    }
    if v3000 {
        source += "M  V30 END BOND\nM  V30 END CTAB\n";
    }
    source + "M  END\n"
}

#[test]
fn molfile_ambiguous_wedge_geometry_never_invents_a_configuration() {
    for points in [
        // Coincident and almost coincident unmarked directions, even with
        // different bond lengths, do not establish a drawing order.
        vec![[0.0, 0.0], [1.0, 0.0], [-1.0, 1.0], [-2.0, 2.0]],
        vec![[0.0, 0.0], [1.0, 0.0], [-1.0, 1.0], [-2.0, 2.01]],
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 0.0], [0.0, 1.0]],
        vec![[0.0, 0.0]; 4],
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [-1.0, 0.0],
            [-2.0, 0.0],
            [-3.0, 0.0],
        ],
    ] {
        for v3000 in [false, true] {
            for scale in [0.1, 1.0, 10.0] {
                for mirror in [-1.0, 1.0] {
                    let source = wedge_drawing(&points, scale, mirror, v3000);
                    let (molecule, report) = read_molfile_with_report(&source).unwrap();
                    assert!(molecule.stereo_elements().next().is_none(), "{source}");
                    assert_eq!(report.warnings().len(), 1, "{source}");
                    assert!(report.created_stereo_elements().is_empty());
                    assert!(!molecule.perception().has_valence());
                }
            }
        }
    }
}

#[test]
fn molfile_specified_tetrahedral_export_requires_a_drawing() {
    let molecule = read_smiles("F[C@](Cl)(Br)I").unwrap();
    for result in [
        molfile::write(
            &molecule,
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V2000,
            },
        ),
        molfile::write(
            &molecule,
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V3000,
            },
        ),
    ] {
        assert!(
            result.is_err(),
            "zero coordinates cannot encode a wedge configuration"
        );
    }
    let model = test_model(&molecule);
    for result in [
        molfile::write(
            &model,
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V2000,
            },
        ),
        molfile::write(
            &model,
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V3000,
            },
        ),
    ] {
        assert!(
            result.is_err(),
            "a degenerate drawing cannot preserve specified stereo"
        );
    }
}

#[test]
fn molfile_valid_wedge_geometry_is_scale_and_reflection_consistent() {
    // RDKit 2026.03.3 independently checks these drawing configurations.
    for (points, original) in [
        (
            vec![[0.0, 0.0], [1.0, 0.0], [-1.0, 0.0], [0.0, 1.0]],
            StereoDescriptor::S,
        ),
        (
            vec![[0.0, 0.0], [1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]],
            StereoDescriptor::R,
        ),
        (
            vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0]],
            StereoDescriptor::R,
        ),
        (
            vec![[0.0, 0.0], [0.0, 0.0], [-1.0, 0.0], [0.0, 1.0]],
            StereoDescriptor::S,
        ),
        // Local drawings from PubChem 10524, with distinct test ligands.
        // The three unmarked bonds lie on the same side of the center;
        // their drawing order must not be read as a displaced tetrahedron.
        (
            vec![
                [6.3981, -0.4499],
                [6.3820, -1.4914],
                [7.2641, 0.0501],
                [6.1950, -0.0800],
                [5.5320, 0.0501],
            ],
            StereoDescriptor::S,
        ),
        (
            vec![
                [7.2641, 1.0501],
                [7.6891, 1.7862],
                [6.8422, 1.0410],
                [7.2641, 0.0501],
                [6.3981, 1.5501],
            ],
            StereoDescriptor::S,
        ),
    ] {
        for v3000 in [false, true] {
            for scale in [0.1, 1.0, 10.0] {
                for mirror in [-1.0, 1.0] {
                    let source = wedge_drawing(&points, scale, mirror, v3000);
                    let (mut molecule, report) = read_molfile_with_report(&source).unwrap();
                    assert!(report.warnings().is_empty(), "{source}");
                    perceive(&mut molecule).unwrap();
                    let assignment = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
                    let expected = if mirror > 0.0 {
                        original
                    } else if original == StereoDescriptor::R {
                        StereoDescriptor::S
                    } else {
                        StereoDescriptor::R
                    };
                    assert_eq!(assignment.assigned.len(), 1, "{source}");
                    assert_eq!(assignment.assigned[0].descriptor, expected, "{source}");
                }
            }
        }
    }
}

#[test]
fn crowded_wedge_drawing_respects_carrier_permutations() {
    let points = [
        [6.3981, -0.4499],
        [6.3820, -1.4914],
        [7.2641, 0.0501],
        [6.1950, -0.0800],
        [5.5320, 0.0501],
    ];
    for (order, expected) in [
        ([2, 3, 4], StereoDescriptor::S),
        ([3, 4, 2], StereoDescriptor::S),
        ([4, 2, 3], StereoDescriptor::S),
        ([2, 4, 3], StereoDescriptor::R),
        ([4, 3, 2], StereoDescriptor::R),
        ([3, 2, 4], StereoDescriptor::R),
    ] {
        let permuted = [
            points[0],
            points[1],
            points[order[0]],
            points[order[1]],
            points[order[2]],
        ];
        for v3000 in [false, true] {
            let source = wedge_drawing(&permuted, 1.0, 1.0, v3000);
            let mut molecule = read_molfile(&source).unwrap();
            perceive(&mut molecule).unwrap();
            let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
            assert_eq!(assigned.assigned.len(), 1);
            assert_eq!(assigned.assigned[0].descriptor, expected, "{source}");
        }
    }
}

#[test]
fn molfile_wedges_use_drawing_geometry_for_hydrogen_and_lone_pair_carriers() {
    // RDKit 2026.03.5 AssignCIPLabels gives S for the original drawing,
    // and R after reflection, for all four fourth-carrier chemistries.
    for (symbol, charge_code) in [("C", 0), ("S", 0), ("S", 3), ("P", 0)] {
        for (mirror, expected) in [(false, StereoDescriptor::S), (true, StereoDescriptor::R)] {
            let source = wedged_tetrahedron(symbol, charge_code, mirror);
            let mut molecule = read_molfile(&source).expect("wedged center interprets");
            perceive(&mut molecule).expect("valence perceives");
            let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
            assert_eq!(assigned.assigned.len(), 1, "{symbol} {charge_code}");
            assert_eq!(assigned.assigned[0].descriptor, expected, "{source}");
        }
    }
}

#[test]
fn molfile_model_writing_preserves_tetrahedral_drawing_orientation() {
    for (symbol, charge_code) in [("C", 0), ("S", 3), ("P", 0)] {
        for mirror in [false, true] {
            let document = molfile::parse_str(&wedged_tetrahedron(symbol, charge_code, mirror))
                .expect("valid source drawing");
            let interpreted = molfile::interpret(&document).expect("interpreted drawing");
            let original = interpreted.molecules().next().unwrap();
            for written in [
                molfile::write(
                    interpreted.model(),
                    molfile::MolfileWriteOptions {
                        version: molfile::MolfileWriteVersion::V2000,
                    },
                )
                .expect("V2000 model writes"),
                molfile::write(
                    interpreted.model(),
                    molfile::MolfileWriteOptions {
                        version: molfile::MolfileWriteVersion::V3000,
                    },
                )
                .expect("V3000 model writes"),
            ] {
                let reparsed = read_molfile(&written).expect("written drawing interprets");
                assert_eq!(
                    original
                        .stereo_elements()
                        .map(|(_, element)| &element.kind)
                        .collect::<Vec<_>>(),
                    reparsed
                        .stereo_elements()
                        .map(|(_, element)| &element.kind)
                        .collect::<Vec<_>>(),
                    "{written}"
                );
            }
        }
    }
}

#[test]
fn molfile_redundant_wedges_preserve_consistent_and_unknown_configurations() {
    let source = wedged_tetrahedron("C", 0, false)
        .replace("   -1.0000    0.0000", "   -0.5000    0.8660")
        .replace("    0.0000    1.0000", "   -0.5000   -0.8660");
    let original = read_molfile(&source).expect("single wedge interprets");
    let redundant = source.replace("  1  3  1  0", "  1  3  1  1");
    let (molecule, report) =
        read_molfile_with_report(&redundant).expect("redundant wedges interpret");
    assert!(report.warnings().is_empty());
    assert_eq!(
        molecule.stereo_elements().next().unwrap().1.kind,
        original.stereo_elements().next().unwrap().1.kind
    );

    let unknown = source.replace("  1  3  1  0", "  1  3  1  4");
    let molecule = read_molfile(&unknown).expect("wavy mark interprets");
    assert!(molecule
        .stereo_elements()
        .next()
        .unwrap()
        .1
        .is_explicitly_unknown());

    let degenerate = unknown
        .replace("    1.0000    0.0000", "    0.0000    0.0000")
        .replace("   -0.5000    0.8660", "    0.0000    0.0000")
        .replace("   -0.5000   -0.8660", "    0.0000    0.0000");
    let (molecule, report) = read_molfile_with_report(&degenerate).unwrap();
    assert!(report.warnings().is_empty());
    assert!(molecule
        .stereo_elements()
        .next()
        .unwrap()
        .1
        .is_explicitly_unknown());

    let conflicting = source.replace("  1  3  1  0", "  1  3  1  6");
    let (molecule, report) = read_molfile_with_report(&conflicting).expect("conflict is reported");
    assert_eq!(report.warnings().len(), 1);
    assert!(molecule.stereo_elements().next().is_none());
}

#[test]
fn v3000_round_trips_absolute_or_and_stereo_groups_and_promotes_auto_output() {
    let document = molfile::parse_str(&wedged_tetrahedron("C", 0, false)).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let source = molfile::write(
        interpreted.model(),
        molfile::MolfileWriteOptions {
            version: molfile::MolfileWriteVersion::V3000,
        },
    )
    .unwrap();
    for (group, expected) in [
        ("MDLV30/STEABS", StereoGroupKind::Absolute),
        ("MDLV30/STEREL1", StereoGroupKind::Or),
        ("MDLV30/STERAC1", StereoGroupKind::And),
    ] {
        let grouped = source.replace(
            "M  V30 END CTAB",
            &format!("M  V30 BEGIN COLLECTION\nM  V30 {group} ATOMS=(1 1)\nM  V30 END COLLECTION\nM  V30 END CTAB"),
        );
        let document = molfile::parse_str(&grouped).expect("collection syntax is preserved");
        let interpreted = molfile::interpret(&document).expect("group semantics are represented");
        let group = interpreted
            .molecules()
            .next()
            .unwrap()
            .stereo_groups()
            .next()
            .unwrap()
            .1;
        assert_eq!(group.kind, expected);
        assert_eq!(group.members.len(), 1);
        let mut molecule = interpreted.molecules().next().unwrap().clone();
        molecule.perceive().unwrap();
        let cx = crate::smiles::write(&molecule, crate::smiles::SmilesWriteOptions::canonical())
            .unwrap();
        let mut restored = crate::smiles::to_molecules(&cx).unwrap().pop().unwrap();
        restored.perceive().unwrap();
        assert_eq!(restored.stereo_groups().next().unwrap().1.kind, expected);
        assert_eq!(
            crate::smiles::write(&restored, crate::smiles::SmilesWriteOptions::canonical())
                .unwrap(),
            cx
        );
        assert!(interpreted.reports()[0].ignored_record_lines().is_empty());
        let written =
            molfile::write(interpreted.model(), molfile::MolfileWriteOptions::default()).unwrap();
        assert!(written.contains("V3000"));
        let document = molfile::parse_str(&written).unwrap();
        let reparsed = molfile::interpret(&document).unwrap();
        assert_eq!(
            reparsed
                .molecules()
                .next()
                .unwrap()
                .stereo_groups()
                .next()
                .unwrap()
                .1,
            group
        );
        assert!(molfile::write(
            interpreted.model(),
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V2000
            }
        )
        .unwrap_err()
        .message()
        .contains("enhanced stereo groups"));
    }
}

#[test]
fn v3000_round_trips_atropisomeric_bond_group_members() {
    let document = molfile::parse_str(rdkit_rp6306_atrop_molblock()).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let source = molfile::write(
        interpreted.model(),
        molfile::MolfileWriteOptions {
            version: molfile::MolfileWriteVersion::V3000,
        },
    )
    .unwrap();
    // RDKit represents enhanced axis membership using either endpoint atom,
    // and collapses a pair of endpoint references to one bond member.
    for (name, kind) in [
        ("STEABS", StereoGroupKind::Absolute),
        ("STERAC1", StereoGroupKind::And),
        ("STEREL1", StereoGroupKind::Or),
    ] {
        for atoms in ["1 3", "1 9", "2 3 9"] {
            let grouped = source.replace("M  V30 END CTAB", &format!("M  V30 BEGIN COLLECTION\nM  V30 MDLV30/{name} ATOMS=({atoms})\nM  V30 END COLLECTION\nM  V30 END CTAB"));
            let interpretation = molfile::parse_str(&grouped).unwrap().interpret().unwrap();
            let molecule = interpretation.molecules().next().unwrap();
            let group = molecule.stereo_groups().next().unwrap().1;
            assert_eq!(group.kind, kind);
            assert_eq!(group.members.len(), 1);
            assert!(
                matches!(&molecule.stereo_element(group.members[0]).unwrap().kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(3))
            );
            let output = molfile::write(
                interpretation.model(),
                molfile::MolfileWriteOptions {
                    version: molfile::MolfileWriteVersion::V3000,
                },
            )
            .unwrap();
            assert!(!output.contains("BONDS="));
            assert!(output.contains("ATOMS=(1 3)"));
            let reread = molfile::parse_str(&output).unwrap().interpret().unwrap();
            assert_eq!(
                reread
                    .molecules()
                    .next()
                    .unwrap()
                    .stereo_groups()
                    .next()
                    .unwrap()
                    .1,
                group
            );
        }
    }
}

#[test]
fn molfile_wedges_respect_pyramidal_lone_pair_endpoints() {
    // A focused aryl sulfoxide drawing reproduces the competing axis/tetrahedral
    // interpretations in supplied PubChem 146091, 461502 and 461520.
    for symbol in ["S", "Se"] {
        let mut descriptors = Vec::new();
        for cfg in [1, 3] {
            let source = format!(
                "pyramidal endpoint\nkekule\n\n  0  0  0     0  0            999 V3000\n\
M  V30 BEGIN CTAB\nM  V30 COUNTS 9 9 0 0 0\nM  V30 BEGIN ATOM\n\
M  V30 1 {symbol} 0 0 0 0\nM  V30 2 O 0 1 0 0\nM  V30 3 C -1 -1 0 0\n\
M  V30 4 C 1 0 0 0\nM  V30 5 C 1.5 0.866 0 0\nM  V30 6 C 2.5 0.866 0 0\n\
M  V30 7 C 3 0 0 0\nM  V30 8 C 2.5 -0.866 0 0\nM  V30 9 C 1.5 -0.866 0 0\n\
M  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 2 1 2\nM  V30 2 1 1 3 CFG={cfg}\n\
M  V30 3 1 1 4\nM  V30 4 2 4 5\nM  V30 5 1 5 6\nM  V30 6 2 6 7\n\
M  V30 7 1 7 8\nM  V30 8 2 8 9\nM  V30 9 1 9 4\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n"
            );
            let mut molecule = read_molfile(&source).unwrap();
            assert!(!molecule.perception().has_valence());
            let elements = molecule.stereo_elements().collect::<Vec<_>>();
            assert_eq!(elements.len(), 1);
            assert!(
                matches!(&elements[0].1.kind, StereoElementKind::Tetrahedral(stereo)
                if stereo.center == AtomId::new(0)
                    && stereo.carriers.contains(&StereoCarrier::ImplicitLonePair))
            );
            perceive(&mut molecule).unwrap();
            let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
            assert_eq!(assigned.assigned.len(), 1);
            descriptors.push(assigned.assigned[0].descriptor);
            let interpretation = molfile::parse_str(&source).unwrap().interpret().unwrap();
            for output in [
                molfile::write(
                    interpretation.model(),
                    molfile::MolfileWriteOptions {
                        version: molfile::MolfileWriteVersion::V2000,
                    },
                )
                .unwrap(),
                molfile::write(
                    interpretation.model(),
                    molfile::MolfileWriteOptions {
                        version: molfile::MolfileWriteVersion::V3000,
                    },
                )
                .unwrap(),
            ] {
                let mut reread = read_molfile(&output).unwrap();
                perceive(&mut reread).unwrap();
                let assigned = stereo_api::assign_cip_descriptors(&mut reread).unwrap();
                assert_eq!(assigned.assigned.len(), 1);
                assert_eq!(
                    assigned.assigned[0].descriptor,
                    *descriptors.last().unwrap()
                );
            }
        }
        assert_ne!(descriptors[0], descriptors[1]);
    }
}

#[test]
fn molfile_atropisomeric_wedges_validate_all_marks_and_preserve_unknown_stereo() {
    let source = rdkit_rp6306_atrop_molblock().replace("  9 12  1  6", "  9 12  1  0");
    let unmarked = read_molfile(&source).unwrap();
    let marked = |left, right| {
        source
            .replace("  3  7  1  0", &format!("  3  7  1  {left}"))
            .replace("  3 10  1  0", &format!("  3 10  1  {right}"))
    };
    // Opposite directions at one end specify an axis. Two wedges or two
    // hashes conflict: preserve the structure without guessing a configuration.
    for (left, right, expected) in [(1, 6, StereoDescriptor::M), (6, 1, StereoDescriptor::P)] {
        let mut molecule = read_molfile(&marked(left, right)).unwrap();
        perceive(&mut molecule).unwrap();
        let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
        assert_eq!(assigned.assigned.len(), 1);
        assert_eq!(assigned.assigned[0].descriptor, expected);
    }
    for direction in [1, 6] {
        let source = marked(direction, direction);
        let (molecule, report) = read_molfile_with_report(&source).unwrap();
        assert_eq!(molecule, unmarked);
        assert!(molecule.stereo_elements().next().is_none());
        assert!(!molecule.perception().has_valence());
        assert!(report.created_stereo_elements().is_empty());
        let [molfile::MolfileInterpretationWarning::ConflictingAtropisomericWedgeMarks {
            axis,
            source_line,
            mark_count,
        }] = report.warnings()
        else {
            panic!("missing conflict diagnostic")
        };
        assert_eq!(*axis, BondId::new(3));
        assert_eq!(*mark_count, 2);
        assert_eq!(
            *source_line,
            report
                .bond_mappings()
                .iter()
                .find(|mapping| mapping.bond() == *axis)
                .unwrap()
                .source_line()
        );
    }
    let model = molfile::parse_str(&source).unwrap().interpret().unwrap();
    let v3000 = molfile::write(
        model.model(),
        molfile::MolfileWriteOptions {
            version: molfile::MolfileWriteVersion::V3000,
        },
    )
    .unwrap();
    let unmarked_v3000 = read_molfile(&v3000).unwrap();
    for cfg in [1, 3] {
        let conflicting = v3000
            .lines()
            .map(|line| {
                if line.starts_with("M  V30 ")
                    && (line.ends_with(" 1 3 7") || line.ends_with(" 1 3 10"))
                {
                    format!("{line} CFG={cfg}\n")
                } else {
                    format!("{line}\n")
                }
            })
            .collect::<String>();
        let document = molfile::parse_str(&conflicting).unwrap();
        let interpretation = document.interpret().unwrap();
        assert_eq!(document.source(), conflicting);
        assert_eq!(interpretation.molecules().next().unwrap(), &unmarked_v3000);
        let report = &interpretation.reports()[0];
        let [molfile::MolfileInterpretationWarning::ConflictingAtropisomericWedgeMarks {
            axis,
            source_line,
            mark_count,
        }] = report.warnings()
        else {
            panic!("missing V3000 conflict diagnostic")
        };
        assert_eq!(*mark_count, 2);
        assert_eq!(
            *source_line,
            report
                .bond_mappings()
                .iter()
                .find(|mapping| mapping.bond() == *axis)
                .unwrap()
                .source_line()
        );
    }
    for (left, right) in [(4, 0), (4, 1), (4, 6), (1, 4), (6, 4), (4, 4)] {
        let mut molecule = read_molfile(&marked(left, right)).unwrap();
        let elements = molecule
            .stereo_elements()
            .map(|(_, element)| element)
            .collect::<Vec<_>>();
        assert_eq!(elements.len(), 1);
        assert!(matches!(&elements[0].kind, StereoElementKind::Axis(stereo)
            if stereo.axis == BondId::new(3) && stereo.orientation.is_none()));
        perceive(&mut molecule).unwrap();
        assert!(stereo_api::assign_cip_descriptors(&mut molecule)
            .unwrap()
            .assigned
            .is_empty());
        assert!(molfile::write(
            &molecule,
            molfile::MolfileWriteOptions {
                version: molfile::MolfileWriteVersion::V3000
            }
        )
        .unwrap_err()
        .message()
        .contains("unknown axis"));
    }
}

#[test]
fn v3000_stereo_groups_validate_members_and_preserve_source_ids_and_continuations() {
    let document = molfile::parse_str(&wedged_tetrahedron("C", 0, false)).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let source = molfile::write(
        interpreted.model(),
        molfile::MolfileWriteOptions {
            version: molfile::MolfileWriteVersion::V3000,
        },
    )
    .unwrap()
    .replace("M  V30 1 C", "M  V30 101 C")
    .replace("M  V30 1 1 1 2", "M  V30 1 1 101 2")
    .replace("M  V30 2 1 1 3", "M  V30 2 1 101 3")
    .replace("M  V30 3 1 1 4", "M  V30 3 1 101 4");
    let collection = |row: &str| {
        source.replace(
            "M  V30 END CTAB",
            &format!(
                "M  V30 BEGIN COLLECTION\nM  V30 {row}\nM  V30 END COLLECTION\nM  V30 END CTAB"
            ),
        )
    };
    let continued = collection("MDLV30/STEREL1 ATOMS=(1 -\nM  V30 101)");
    let document = molfile::parse_str(&continued).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let group = interpreted
        .molecules()
        .next()
        .unwrap()
        .stereo_groups()
        .next()
        .unwrap()
        .1;
    assert_eq!(group.kind, StereoGroupKind::Or);
    assert_eq!(group.members.len(), 1);
    assert!(interpreted.reports()[0].ignored_record_lines().is_empty());

    for row in [
        "MDLV30/STEREL0 ATOMS=(1 101)",
        "MDLV30/STEREL1 ATOMS=(2 101)",
        "MDLV30/STEREL1 ATOMS=(0)",
        "MDLV30/STEREL1 ATOMS=(1 999)",
        "MDLV30/STEREL1 ATOMS=(2 101 101)",
        "MDLV30/STEREL1 ATOMS=(1 2)",
        "MDLV30/STEREL1 ATOMS=(1 101) ATOMS=(1 101)",
        "MDLV30/STEREL1 BONDS=(1 1)",
    ] {
        let grouped = collection(row);
        if let Ok(document) = molfile::parse_str(&grouped) {
            assert!(molfile::interpret(&document).is_err(), "{row}");
        }
    }
}

#[test]
fn v3000_repeated_group_ids_preserve_one_relation() {
    let molecule = read_smiles("F[C@H](Cl)[C@H](Br)I").unwrap();
    let model = Model::from_molecule(
        molecule.clone(),
        &test_positions(vec![
            Point3::new(-1.0, 1.0, 0.0),
            Point3::origin(),
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(2.0, -1.0, 0.0),
        ]),
    )
    .unwrap();
    let source = molfile::write(
        &model,
        molfile::MolfileWriteOptions {
            version: molfile::MolfileWriteVersion::V3000,
        },
    )
    .unwrap();
    for kind in ["STEREL", "STERAC"] {
        for number in ["1", "01"] {
            let grouped = source.replace("M  V30 END CTAB", &format!("M  V30 BEGIN COLLECTION\nM  V30 MDLV30/{kind}1 ATOMS=(1 2)\nM  V30 MDLV30/{kind}{number} ATOMS=(1 4)\nM  V30 END COLLECTION\nM  V30 END CTAB"));
            let document = molfile::parse_str(&grouped).unwrap();
            let interpreted = molfile::interpret(&document).unwrap();
            let molecule = interpreted.molecules().next().unwrap();
            assert_eq!(molecule.stereo_groups().count(), 1);
            assert_eq!(molecule.stereo_groups().next().unwrap().1.members.len(), 2);
            assert!(interpreted.reports()[0].ignored_record_lines().is_empty());
        }
    }
}

#[test]
fn v3000_rejects_relative_groups_across_components_and_preserves_absolute_members() {
    let document = molfile::parse_str(&wedged_tetrahedron("C", 0, false)).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let source = molfile::write(interpreted.model(), molfile::MolfileWriteOptions { version: molfile::MolfileWriteVersion::V3000 }).unwrap()
        .replace("COUNTS 4 3", "COUNTS 8 6")
        .replace("M  V30 END ATOM", "M  V30 5 C 5 0 0 0\nM  V30 6 F 6 0 0 0\nM  V30 7 Cl 4 0 0 0\nM  V30 8 Br 5 1 0 0\nM  V30 END ATOM")
        .replace("M  V30 END BOND", "M  V30 4 1 5 6 CFG=1\nM  V30 5 1 5 7\nM  V30 6 1 5 8\nM  V30 END BOND");
    for group in ["MDLV30/STERAC1", "MDLV30/STEREL1", "MDLV30/STEABS"] {
        let source = source.replace("M  V30 END CTAB", &format!("M  V30 BEGIN COLLECTION\nM  V30 {group} ATOMS=(2 1 5)\nM  V30 END COLLECTION\nM  V30 END CTAB"));
        let document = molfile::parse_str(&source).unwrap();
        if group == "MDLV30/STEABS" {
            let interpreted = molfile::interpret(&document).unwrap();
            assert_eq!(interpreted.molecules().count(), 2);
            assert!(interpreted
                .molecules()
                .all(|molecule| molecule.stereo_groups().count() == 1));
            let output = molfile::write(
                interpreted.model(),
                molfile::MolfileWriteOptions {
                    version: molfile::MolfileWriteVersion::V3000,
                },
            )
            .unwrap();
            assert_eq!(output.matches("MDLV30/STEABS").count(), 1);
            assert!(output.contains("ATOMS=(2 1 5)"));
            assert_eq!(
                molfile::parse_str(&output)
                    .unwrap()
                    .to_molecules()
                    .unwrap()
                    .len(),
                2
            );
        } else {
            assert!(molfile::interpret(&document)
                .unwrap_err()
                .message()
                .contains("spanning disconnected molecules"));
        }
    }
}

#[test]
fn molfile_wedge_assembles_tetrahedral_p_with_a_double_bond() {
    let input = r#"tetrahedral phosphorus
  kekule

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  2  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"#;
    let molecule = read_sdf_molecules(input)
        .expect("compact phosphorus regression parses")
        .into_iter()
        .next()
        .expect("one molecule");

    assert_eq!(molecule.stereo_elements().count(), 1);
    let element = molecule
        .stereo_elements()
        .next()
        .expect("created tetrahedral element")
        .1;
    assert!(matches!(
        &element.kind,
        StereoElementKind::Tetrahedral(stereo) if stereo.center == AtomId::new(0)
    ));
}

#[test]
fn molfile_wedge_assembles_pyramidal_s_with_a_lone_pair() {
    let input = r#"pyramidal sulfur
  kekule

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  2  0  0  0  0
  1  4  1  0  0  0  0
M  END
$$$$
"#;
    let molecule = read_sdf_molecules(input)
        .expect("compact sulfur regression parses")
        .into_iter()
        .next()
        .expect("one molecule");

    assert_eq!(molecule.stereo_elements().count(), 1);
    let element = molecule
        .stereo_elements()
        .next()
        .expect("created tetrahedral element")
        .1;
    assert!(matches!(
        &element.kind,
        StereoElementKind::Tetrahedral(stereo)
            if stereo.center == AtomId::new(0)
                && stereo.carriers.contains(&StereoCarrier::ImplicitLonePair)
    ));
}
