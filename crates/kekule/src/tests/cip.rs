use super::*;

fn assign_cip(mol: &mut Molecule) -> CipAssignmentReport {
    stereo_api::assign_cip_descriptors(mol).expect("CIP assignment should succeed")
}

fn installed_cip_descriptors(mol: &Molecule) -> Vec<(StereoElementId, StereoDescriptor)> {
    mol.perception().cip_descriptors().collect()
}

fn assigned_descriptors(report: &CipAssignmentReport) -> Vec<StereoDescriptor> {
    report
        .assigned
        .iter()
        .map(|assignment| assignment.descriptor)
        .collect()
}

fn assert_cip_not_stereogenic(mol: &mut Molecule, element: StereoElementId) {
    let report = assign_cip(mol);
    assert!(report.assigned.is_empty());
    assert_eq!(
        report.skipped,
        vec![CipSkipped {
            element,
            reason: CipSkippedReason::NotStereogenic,
        }]
    );
    assert_eq!(mol.cip_descriptor(element).expect("stereo element"), None);
}

#[test]
fn successful_cip_with_no_assignments_installs_empty_stereo_section() {
    let mut mol = crate::core::MoleculeEditor::new();
    mol.add_atom(carbon()).expect("single carbon");
    assert!(!mol.perception().has_stereo());

    let report = assign_cip(mol.working_mut());

    assert!(report.assigned.is_empty());
    assert!(report.skipped.is_empty());
    assert!(mol.perception().has_stereo());
    assert!(mol
        .perception()
        .stereo_state()
        .expect("installed stereo section")
        .cip_descriptors()
        .next()
        .is_none());
}

#[test]
fn cip_assigns_tetrahedral_descriptors_from_stored_local_stereo() {
    let mut s_alanine = read_smiles("C[C@@H](C(=O)O)N").expect("alanine parses");
    perceive(&mut s_alanine).expect("alanine perceives");

    let report = assign_cip(&mut s_alanine);

    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::S,
        }]
    );
    assert_eq!(
        s_alanine
            .cip_descriptor(StereoElementId::new(0))
            .expect("stereo element"),
        Some(StereoDescriptor::S)
    );

    let mut r_alanine = read_smiles("C[C@H](C(=O)O)N").expect("alanine parses");
    perceive(&mut r_alanine).expect("alanine perceives");

    let report = assign_cip(&mut r_alanine);

    assert_eq!(report.assigned[0].descriptor, StereoDescriptor::R);
}

#[test]
fn cip_matches_rdkit_for_molfile_wedge_up_and_down() {
    for (stereo_code, expected) in [(1, StereoDescriptor::S), (6, StereoDescriptor::R)] {
        let input = format!(
            "wedge {stereo_code}
kekule

  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  {stereo_code}  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
"
        );
        let mut molecule = read_molfile(&input).expect("wedge molfile parses");
        perceive(&mut molecule).expect("wedge molfile perceives");

        let report = assign_cip(&mut molecule);

        assert_eq!(
            report.assigned,
            vec![CipAssignment {
                element: StereoElementId::new(0),
                descriptor: expected,
            }]
        );
    }
}

#[test]
fn cip_matches_rdkit_for_molfile_implicit_h_wedge_geometry() {
    let mut molecule = read_molfile(implicit_h_wedge_geometry_molblock())
        .expect("implicit-H wedge molfile parses");
    perceive(&mut molecule).expect("implicit-H wedge molfile perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::S,
        }]
    );
}

#[test]
fn cip_assigns_axis_descriptors_from_ranked_anchors() {
    for (left_reference, right_reference, orientation, expected) in [
        (
            AxisReference::Low,
            AxisReference::Low,
            AxisOrientation::CounterClockwise,
            StereoDescriptor::M,
        ),
        (
            AxisReference::High,
            AxisReference::Low,
            AxisOrientation::CounterClockwise,
            StereoDescriptor::P,
        ),
        (
            AxisReference::High,
            AxisReference::High,
            AxisOrientation::Clockwise,
            StereoDescriptor::P,
        ),
    ] {
        let (mut mol, stereo) = axis_stereo_graph(left_reference, right_reference, orientation);

        let report = assign_cip(&mut mol);

        assert_eq!(
            report.assigned,
            vec![CipAssignment {
                element: stereo,
                descriptor: expected,
            }]
        );
        assert_eq!(
            mol.cip_descriptor(stereo).expect("axis element"),
            Some(expected)
        );
    }
}

#[test]
fn cip_skips_axis_with_equivalent_endpoint_ligands_as_nonstereogenic() {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left_high = mol
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");
    let left_low = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let axis = mol.add_bond(left, right, BondOrder::Single).expect("axis");
    for carrier in [left_high, left_low] {
        mol.add_bond(left, carrier, BondOrder::Single)
            .expect("left carrier bond");
    }
    for carrier in [right_a, right_b] {
        mol.add_bond(right, carrier, BondOrder::Single)
            .expect("right carrier bond");
    }
    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![StereoCarrier::Atom(left_high), StereoCarrier::Atom(right_a)],
            orientation: Some(AxisOrientation::CounterClockwise),
        })))
        .expect("axis stereo element");

    assert_cip_not_stereogenic(mol.working_mut(), stereo);
}

#[test]
fn cip_assigns_axis_descriptor_after_coordinate_stereo_materialization() {
    let (molecule, positions, _axis) = coordinate_axis_graph(true);
    let mut mol = molecule.edit();
    let materialization_report = stereo_api::materialize_coordinate_stereo_with_options(
        &mut mol,
        &positions,
        CoordinateStereoOptions { infer_axes: true },
    )
    .expect("coordinate axis materialization");
    assert_eq!(materialization_report.created_elements.len(), 1);

    let mut mol = mol.finish().expect("materialized molecule publishes");

    let report = assign_cip(&mut mol);

    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::P,
        }]
    );
    assert_eq!(
        mol.cip_descriptor(StereoElementId::new(0))
            .expect("axis element"),
        Some(StereoDescriptor::P)
    );
}

#[test]
fn cip_assigns_pseudo_axis_descriptors_for_pseudoasymmetric_endpoint_ordering() {
    for (orientation, expected) in [
        (AxisOrientation::CounterClockwise, StereoDescriptor::LowerM),
        (AxisOrientation::Clockwise, StereoDescriptor::LowerP),
    ] {
        let (mut mol, axis_element) = pseudoasymmetric_axis_graph(orientation);

        assign_cip(&mut mol);

        assert_eq!(
            mol.cip_descriptor(axis_element).expect("axis stereo"),
            Some(expected)
        );
    }
}

#[test]
fn cip_axis_with_two_enantiomorphic_endpoint_pairs_is_absolute() {
    for (orientation, expected) in [
        (AxisOrientation::CounterClockwise, StereoDescriptor::M),
        (AxisOrientation::Clockwise, StereoDescriptor::P),
    ] {
        let mut editor = MoleculeEditor::new();
        let left = editor.add_atom(carbon()).expect("left endpoint");
        let right = editor.add_atom(carbon()).expect("right endpoint");
        let axis = editor
            .add_bond(left, right, BondOrder::Single)
            .expect("axis bond");
        let (left_r, _) = add_enantiomorphic_tetrahedral_carriers(editor.working_mut(), left);
        let (right_r, _) = add_enantiomorphic_tetrahedral_carriers(editor.working_mut(), right);
        let element = editor
            .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
                axis,
                carriers: vec![StereoCarrier::Atom(left_r), StereoCarrier::Atom(right_r)],
                orientation: Some(orientation),
            })))
            .expect("axis assertion");
        let mut molecule = editor.finish().expect("molecule");
        assign_cip(&mut molecule);
        assert_eq!(
            molecule.cip_descriptor(element).expect("axis"),
            Some(expected)
        );
    }
}

#[test]
fn cip_depth_limit_cannot_promote_isotope_rule_or_claim_nonstereogenicity() {
    for smiles in ["N[C@@H]([13CH2]CCC)CCCCO", "F[C@](Cl)(CCCC)CCCC"] {
        let mut molecule = read_smiles(smiles).expect("depth regression");
        perceive(&mut molecule).expect("perceive");
        assign_cip(&mut molecule);
        let previous = molecule.perception().clone();
        let error = stereo_api::assign_cip_descriptors_with_options(
            &mut molecule,
            CipAssignmentOptions {
                max_depth: 0,
                ..CipAssignmentOptions::default()
            },
        )
        .expect_err("a truncated constitutional tie cannot be resolved");
        assert_eq!(
            error.issues,
            vec![CipAssignmentIssue::DepthLimitExceeded {
                element: StereoElementId::new(0),
                max_depth: 0,
            }]
        );
        assert_eq!(molecule.perception(), &previous);
    }
}

#[test]
fn cip_matches_rdkit_for_molfile_atropisomeric_axis() {
    let mut molecule =
        read_molfile(rdkit_rp6306_atrop_molblock()).expect("RDKit atropisomer fixture parses");
    perceive(&mut molecule).expect("atropisomer fixture perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::P,
        }]
    );
    assert_eq!(
        molecule
            .cip_descriptor(StereoElementId::new(0))
            .expect("axis stereo element"),
        Some(StereoDescriptor::P)
    );
}

#[test]
fn cip_matches_rdkit_for_alternate_molfile_atropisomeric_wedge() {
    let mut molecule = read_molfile(rdkit_rp6306_atrop3_molblock())
        .expect("RDKit alternate atropisomer fixture parses");
    perceive(&mut molecule).expect("atropisomer fixture perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::P,
        }]
    );
}

#[test]
fn cip_axis_ranking_is_stable_across_all_carbon_aromatic_source_kekule_variants() {
    let mut molecule = read_molfile(rdkit_rp6306_atrop4_molblock())
        .expect("fully declared RDKit atropisomer fixture interprets");
    perceive(&mut molecule).expect("atropisomer fixture perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(assigned_descriptors(&report), vec![StereoDescriptor::P]);

    let mut molecule = read_molfile(rdkit_bms986142_atrop4_molblock())
        .expect("omitted tetrahedral hydrogen interprets");
    perceive(&mut molecule).expect("BMS atropisomer perceives");
    assert_eq!(
        assigned_descriptors(&assign_cip(&mut molecule)),
        vec![StereoDescriptor::S, StereoDescriptor::P]
    );
}

#[test]
fn cip_axis_ranking_preserves_heteromancude_source_kekule_guardrail() {
    let mut molecule = read_molfile(rdkit_jdq443_atrop1_molblock())
        .expect("RDKit JDQ443 atropisomer fixture parses");
    perceive(&mut molecule).expect("JDQ443 atropisomer fixture perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::M,
        }]
    );
}

#[test]
fn molfile_interpretation_preserves_atrop_fixtures_with_omitted_hydrogens() {
    for (fixture, expected) in [
        (
            rdkit_bms986142_atrop5_molblock(),
            [StereoDescriptor::S, StereoDescriptor::P],
        ),
        (
            rdkit_zm374979_atrop1_molblock(),
            [StereoDescriptor::R, StereoDescriptor::M],
        ),
        (
            rdkit_zm374979_atrop2_molblock(),
            [StereoDescriptor::R, StereoDescriptor::M],
        ),
    ] {
        let mut molecule = read_molfile(fixture).expect("source hydrogen interprets");
        perceive(&mut molecule).expect("atropisomer perceives");
        assert_eq!(assigned_descriptors(&assign_cip(&mut molecule)), expected);
    }
}

#[test]
fn cip_matches_rdkit_for_ring_internal_molfile_atropisomeric_axis() {
    for (fixture, expected) in [
        (
            rdkit_macrocycle8_ortho_wedge_molblock(),
            StereoDescriptor::M,
        ),
        (rdkit_macrocycle8_ortho_hash_molblock(), StereoDescriptor::P),
    ] {
        let mut molecule =
            read_molfile(fixture).expect("RDKit macrocycle atropisomer fixture parses");
        perceive(&mut molecule).expect("macrocycle atropisomer fixture perceives");

        let report = assign_cip(&mut molecule);

        assert_eq!(report.assigned.len(), 1);
        assert_eq!(report.assigned[0].descriptor, expected);
    }
}

#[test]
fn cip_matches_rdkit_for_pubchem_start_atom_bracket_h_tetrahedral_centers() {
    let mut molecule = read_smiles("[C@@H]([C@H](C(=O)O)O)(C(=O)O)O").expect("tartrate parses");
    perceive(&mut molecule).expect("tartrate perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![StereoDescriptor::R, StereoDescriptor::R]
    );
}

#[derive(Debug, Clone, Copy)]
enum AxisReference {
    High,
    Low,
}

fn axis_stereo_graph(
    left_reference: AxisReference,
    right_reference: AxisReference,
    orientation: AxisOrientation,
) -> (Molecule, StereoElementId) {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left_high = mol
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");
    let left_low = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right_high = mol
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let right_low = mol.add_atom(carbon()).expect("atom identifier capacity");
    let axis = mol.add_bond(left, right, BondOrder::Single).expect("axis");
    mol.add_bond(left, left_high, BondOrder::Single)
        .expect("left high");
    mol.add_bond(left, left_low, BondOrder::Single)
        .expect("left low");
    mol.add_bond(right, right_high, BondOrder::Single)
        .expect("right high");
    mol.add_bond(right, right_low, BondOrder::Single)
        .expect("right low");

    let left_carrier = match left_reference {
        AxisReference::High => left_high,
        AxisReference::Low => left_low,
    };
    let right_carrier = match right_reference {
        AxisReference::High => right_high,
        AxisReference::Low => right_low,
    };
    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![
                StereoCarrier::Atom(left_carrier),
                StereoCarrier::Atom(right_carrier),
            ],
            orientation: Some(orientation),
        })))
        .expect("axis stereo element");
    (mol.finish().expect("connected axis fixture"), stereo)
}

#[test]
fn cip_matches_rdkit_for_smiles_ring_digit_tetrahedral_order() {
    let mut molecule = read_smiles("CC(C)C[C@@H]1CN2CCC3=CC(=C(C=C3C2CC1=O)OC)O[11CH3]")
        .expect("ring chiral molecule parses");
    perceive(&mut molecule).expect("ring chiral molecule perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::R,
        }]
    );
}

#[test]
fn cip_matches_rdkit_for_branch_preserving_sugar_ligand_ranking() {
    let mut molecule =
        read_smiles("C1=C2C(=NC=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O")
            .expect("nucleotide sugar parses");
    perceive(&mut molecule).expect("nucleotide sugar perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![
            StereoDescriptor::R,
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::R,
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_fused_ring_paired_breadth_first_ranking() {
    let mut molecule =
        read_smiles("CC(=O)OC[C@]1([C@@H](CC[C@@]2(C1C[C@@H]([C@]34[C@H]2CC[C@@H](C3)C(=C)C4)OC(=O)C5=CC=C(C=C5)OC)C)OC(=O)C6=CC=C(C=C6)OC)C")
            .expect("polycycle parses");
    perceive(&mut molecule).expect("polycycle perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![
            StereoDescriptor::S,
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::S,
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::S,
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_polyene_directional_double_bonds() {
    let mut molecule =
        read_smiles("CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(C)C=C)/C").expect("polyene parses");
    perceive(&mut molecule).expect("polyene perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![
            StereoDescriptor::E,
            StereoDescriptor::E,
            StereoDescriptor::E
        ]
    );
}

#[test]
fn cip_skips_small_ring_double_bond_stereo_but_assigns_cyclooctene() {
    let mut cyclohexene = read_smiles("C1C=CCCC1").expect("cyclohexene parses");
    perceive(&mut cyclohexene).expect("cyclohexene perceives");
    assert!(cyclohexene
        .stereo_elements()
        .all(|(_, element)| !matches!(element.kind, StereoElementKind::DoubleBond(_))));

    let cip_report = assign_cip(&mut cyclohexene);
    assert!(cip_report.assigned.is_empty());

    let mut cyclooctene = read_smiles(r"C1/C=C\CCCCC1").expect("marked cyclooctene parses");
    perceive(&mut cyclooctene).expect("marked cyclooctene perceives");
    let cip_report = assign_cip(&mut cyclooctene);

    assert_eq!(assigned_descriptors(&cip_report), vec![StereoDescriptor::Z]);
}

#[test]
fn cip_labels_explicitly_asserted_small_ring_double_bond() {
    let mut mol = crate::core::MoleculeEditor::new();
    let atoms = (0..6)
        .map(|_| mol.add_atom(carbon()).expect("atom identifier capacity"))
        .collect::<Vec<_>>();
    let double_bond = mol
        .add_bond(atoms[0], atoms[1], BondOrder::Double)
        .expect("double bond");
    for (left, right) in [
        (atoms[1], atoms[2]),
        (atoms[2], atoms[3]),
        (atoms[3], atoms[4]),
        (atoms[4], atoms[5]),
        (atoms[5], atoms[0]),
    ] {
        mol.add_bond(left, right, BondOrder::Single)
            .expect("ring bond");
    }
    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: double_bond,
                left: atoms[0],
                right: atoms[1],
                left_carrier: StereoCarrier::Atom(atoms[5]),
                right_carrier: StereoCarrier::Atom(atoms[2]),
                orientation: Some(DoubleBondOrientation::Together),
            },
        )))
        .expect("double-bond stereo element");

    let report = assign_cip(mol.working_mut());
    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: stereo,
            descriptor: StereoDescriptor::Z
        }]
    );
}

#[test]
fn cip_labels_asserted_ring_bonds_independently_of_candidate_perception() {
    // RDKit 2026.03.5 AssignCIPLabels ranks explicit CIS/TRANS assertions on
    // small rings, aromatic Kekule bonds, and endocyclic imines alike.
    for smiles in [
        "C1=CCCCC1",
        "C1=CCCCCC1",
        "C1=CCCCCCC1",
        "C1=NCCCCCCC1",
        "C1=CC=CC=C1",
    ] {
        for (orientation, expected) in [
            (DoubleBondOrientation::Together, StereoDescriptor::Z),
            (DoubleBondOrientation::Opposite, StereoDescriptor::E),
        ] {
            let mut molecule = read_smiles(smiles).expect("ring molecule parses");
            perceive(&mut molecule).expect("ring molecule perceives");
            let (bond, left, right) = molecule
                .bonds()
                .find_map(|(id, bond)| {
                    (bond.order == BondOrder::Double).then_some((id, bond.a(), bond.b()))
                })
                .expect("localized double bond");
            let left_carrier = molecule
                .neighbors(left)
                .unwrap()
                .find(|atom| *atom != right)
                .unwrap();
            let right_carrier = molecule
                .neighbors(right)
                .unwrap()
                .find(|atom| *atom != left)
                .unwrap();
            let element = molecule
                .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
                    DoubleBondStereo {
                        bond,
                        left,
                        right,
                        left_carrier: StereoCarrier::Atom(left_carrier),
                        right_carrier: StereoCarrier::Atom(right_carrier),
                        orientation: Some(orientation),
                    },
                )))
                .expect("explicit bond assertion");
            let report = assign_cip(&mut molecule);
            assert_eq!(
                report.assigned,
                vec![CipAssignment {
                    element,
                    descriptor: expected
                }],
                "{smiles}"
            );
        }
    }
}

#[test]
fn cip_skips_double_bond_with_equivalent_endpoint_ligands_as_nonstereogenic() {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let double_bond = mol
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    let fluorine = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let chlorine = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let methyl_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let methyl_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    for carrier in [fluorine, chlorine] {
        mol.add_bond(left, carrier, BondOrder::Single)
            .expect("left carrier bond");
    }
    for carrier in [methyl_a, methyl_b] {
        mol.add_bond(right, carrier, BondOrder::Single)
            .expect("right carrier bond");
    }
    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: double_bond,
                left,
                right,
                left_carrier: StereoCarrier::Atom(fluorine),
                right_carrier: StereoCarrier::Atom(methyl_a),
                orientation: Some(DoubleBondOrientation::Together),
            },
        )))
        .expect("double-bond stereo element");

    assert_cip_not_stereogenic(mol.working_mut(), stereo);
}

#[test]
fn cip_skips_endocyclic_kekule_bond_stereo_after_ring_perception() {
    let mut molecule =
        read_smiles("CC\\1=C(/C/2=C/C3=C(C(=C(N3)/C=C\\4/[C@@](C(=C(N4)/C=C\\5/[C@@](C(=C(N5)/C=C1\\N2)O)(C)CC(=O)O)O)(C)CC(=O)O)C)CCC(=O)O)CCC(=O)O")
            .expect("CID 445170 parses");
    perceive(&mut molecule).expect("CID 445170 perceives");

    assign_cip(&mut molecule);

    let bond_descriptors = double_bond_descriptor_map(&molecule);
    assert_eq!(
        bond_descriptors,
        vec![
            (3, 4, StereoDescriptor::Z),
            (10, 11, StereoDescriptor::Z),
            (16, 17, StereoDescriptor::Z),
            (22, 23, StereoDescriptor::Z),
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_large_fused_ring_with_many_centers() {
    let mut molecule =
        read_smiles("CN1CC[C@@]23[C@H]4[C@H]1CC5=C2C(=C(C=C5)OC)O[C@@H]3[C@]6(C4)C(=O)C7=C8N6CCC9=C8C(=C(C=C9)OC)OC1=C7C=CC(=C1O)OC")
            .expect("fused ring parses");
    perceive(&mut molecule).expect("fused ring perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::R
        ]
    );
}

#[test]
fn cip_assigns_double_bond_descriptors_from_ranked_carriers() {
    let mut together = read_smiles("C(=C\\F)\\F").expect("alkene parses");
    perceive(&mut together).expect("alkene perceives");

    let report = assign_cip(&mut together);

    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::Z,
        }]
    );

    let mut opposite = read_smiles("C(=C/F)\\F").expect("alkene parses");
    perceive(&mut opposite).expect("alkene perceives");

    let report = assign_cip(&mut opposite);

    assert_eq!(report.assigned[0].descriptor, StereoDescriptor::E);
}

#[test]
fn cip_assigns_sequence_descriptors_for_pseudoasymmetric_double_bond_endpoints() {
    for (orientation, expected) in [
        (DoubleBondOrientation::Together, StereoDescriptor::SeqCis),
        (DoubleBondOrientation::Opposite, StereoDescriptor::SeqTrans),
    ] {
        let (mut mol, double_bond_element) = pseudoasymmetric_double_bond_graph(orientation);

        assign_cip(&mut mol);

        assert_eq!(
            mol.cip_descriptor(double_bond_element)
                .expect("double-bond stereo"),
            Some(expected)
        );
    }
}

#[test]
fn cip_uses_rule3_embedded_e_z_descriptors_to_order_ligands() {
    let mut molecule = read_smiles("Br[C@H](/C=C/F)/C=C\\F").expect("Rule 3 alkene pair parses");
    perceive(&mut molecule).expect("Rule 3 alkene pair perceives");

    assign_cip(&mut molecule);

    let atom_descriptors = tetrahedral_descriptor_map(&molecule);
    let bond_descriptors = double_bond_descriptor_map(&molecule);

    assert_eq!(atom_descriptors, vec![(1, StereoDescriptor::R)]);
    assert_eq!(
        bond_descriptors,
        vec![(2, 3, StereoDescriptor::E), (5, 6, StereoDescriptor::Z)]
    );
}

fn pseudoasymmetric_double_bond_graph(
    orientation: DoubleBondOrientation,
) -> (Molecule, StereoElementId) {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let double_bond = mol
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    let (child_r, _) = add_enantiomorphic_tetrahedral_carriers(mol.working_mut(), left);

    let chlorine = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let fluorine = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");

    for carrier in [chlorine, fluorine] {
        mol.add_bond(right, carrier, BondOrder::Single)
            .expect("right carrier bond");
    }
    let double_bond_element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: double_bond,
                left,
                right,
                left_carrier: StereoCarrier::Atom(child_r),
                right_carrier: StereoCarrier::Atom(chlorine),
                orientation: Some(orientation),
            },
        )))
        .expect("double-bond stereo element");

    (
        mol.finish().expect("connected double-bond fixture"),
        double_bond_element,
    )
}

fn pseudoasymmetric_axis_graph(orientation: AxisOrientation) -> (Molecule, StereoElementId) {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let axis = mol.add_bond(left, right, BondOrder::Single).expect("axis");
    let (child_r, _child_s) = add_enantiomorphic_tetrahedral_carriers(mol.working_mut(), left);

    let chlorine = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let fluorine = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    for carrier in [chlorine, fluorine] {
        mol.add_bond(right, carrier, BondOrder::Single)
            .expect("right carrier bond");
    }

    let axis_element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![StereoCarrier::Atom(child_r), StereoCarrier::Atom(chlorine)],
            orientation: Some(orientation),
        })))
        .expect("axis stereo element");

    (mol.finish().expect("connected axis fixture"), axis_element)
}

fn add_enantiomorphic_tetrahedral_carriers(mol: &mut Molecule, parent: AtomId) -> (AtomId, AtomId) {
    let child_r = mol.add_atom(carbon()).expect("atom identifier capacity");
    let child_r_oxygen = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let child_r_nitrogen = mol
        .add_atom(element_atom("N"))
        .expect("atom identifier capacity");

    let child_s = mol.add_atom(carbon()).expect("atom identifier capacity");
    let child_s_oxygen = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let child_s_nitrogen = mol
        .add_atom(element_atom("N"))
        .expect("atom identifier capacity");

    for carrier in [child_r, child_s] {
        mol.add_bond(parent, carrier, BondOrder::Single)
            .expect("parent carrier bond");
    }
    for (child, oxygen, nitrogen) in [
        (child_r, child_r_oxygen, child_r_nitrogen),
        (child_s, child_s_oxygen, child_s_nitrogen),
    ] {
        mol.add_bond(child, oxygen, BondOrder::Single)
            .expect("child oxygen bond");
        mol.add_bond(child, nitrogen, BondOrder::Single)
            .expect("child nitrogen bond");
    }

    mol.add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
        TetrahedralStereo {
            center: child_r,
            carriers: vec![
                StereoCarrier::Atom(parent),
                StereoCarrier::Atom(child_r_oxygen),
                StereoCarrier::Atom(child_r_nitrogen),
                StereoCarrier::ImplicitHydrogen,
            ],
            orientation: Some(TetrahedralOrientation::CounterClockwise),
        },
    )))
    .expect("R child stereo element");
    mol.add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
        TetrahedralStereo {
            center: child_s,
            carriers: vec![
                StereoCarrier::Atom(parent),
                StereoCarrier::Atom(child_s_oxygen),
                StereoCarrier::Atom(child_s_nitrogen),
                StereoCarrier::ImplicitHydrogen,
            ],
            orientation: Some(TetrahedralOrientation::Clockwise),
        },
    )))
    .expect("S child stereo element");

    mol.set_implicit_hydrogens(child_r, 1);
    mol.set_implicit_hydrogens(child_s, 1);

    (child_r, child_s)
}

#[test]
fn cip_assigns_pseudoasymmetric_lowercase_descriptor_from_enantiomorphic_ligands() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let chlorine = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let fluorine = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");

    let child_r = mol.add_atom(carbon()).expect("atom identifier capacity");
    let child_r_oxygen = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let child_r_nitrogen = mol
        .add_atom(element_atom("N"))
        .expect("atom identifier capacity");

    let child_s = mol.add_atom(carbon()).expect("atom identifier capacity");
    let child_s_oxygen = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let child_s_nitrogen = mol
        .add_atom(element_atom("N"))
        .expect("atom identifier capacity");

    for carrier in [chlorine, fluorine, child_r, child_s] {
        mol.add_bond(center, carrier, BondOrder::Single)
            .expect("parent carrier bond");
    }
    for (child, oxygen, nitrogen) in [
        (child_r, child_r_oxygen, child_r_nitrogen),
        (child_s, child_s_oxygen, child_s_nitrogen),
    ] {
        mol.add_bond(child, oxygen, BondOrder::Single)
            .expect("child oxygen bond");
        mol.add_bond(child, nitrogen, BondOrder::Single)
            .expect("child nitrogen bond");
    }

    let child_r_element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center: child_r,
                carriers: vec![
                    StereoCarrier::Atom(center),
                    StereoCarrier::Atom(child_r_oxygen),
                    StereoCarrier::Atom(child_r_nitrogen),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: Some(TetrahedralOrientation::CounterClockwise),
            },
        )))
        .expect("R child stereo element");
    let child_s_element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center: child_s,
                carriers: vec![
                    StereoCarrier::Atom(center),
                    StereoCarrier::Atom(child_s_oxygen),
                    StereoCarrier::Atom(child_s_nitrogen),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("S child stereo element");
    let parent_element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(chlorine),
                    StereoCarrier::Atom(fluorine),
                    StereoCarrier::Atom(child_r),
                    StereoCarrier::Atom(child_s),
                ],
                orientation: Some(TetrahedralOrientation::CounterClockwise),
            },
        )))
        .expect("parent pseudoasymmetric stereo element");

    mol.working_mut().set_implicit_hydrogens(center, 0);
    mol.working_mut().set_implicit_hydrogens(child_r, 1);
    mol.working_mut().set_implicit_hydrogens(child_s, 1);

    assign_cip(mol.working_mut());

    assert_eq!(
        mol.cip_descriptor(child_r_element).expect("R child stereo"),
        Some(StereoDescriptor::R)
    );
    assert_eq!(
        mol.cip_descriptor(child_s_element).expect("S child stereo"),
        Some(StereoDescriptor::S)
    );
    assert_eq!(
        mol.cip_descriptor(parent_element).expect("parent stereo"),
        Some(StereoDescriptor::LowerR)
    );
}

#[test]
fn cip_bootstraps_coupled_pseudoasymmetric_tetrahedral_centers() {
    let mut molecule = read_smiles("CC1=NC(=NN1)[C@@H]2CC[C@H](CC2)NC3CCC3CC(C)C")
        .expect("para-stereo scaffold parses");
    perceive(&mut molecule).expect("para-stereo scaffold perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(report.assigned.len(), 2);
    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![(6, StereoDescriptor::LowerR), (9, StereoDescriptor::LowerR)]
    );
}

#[test]
fn cip_matches_rdkit_for_para_stereochemistry_with_directional_double_bonds() {
    for (smiles, expected) in [
        (
            r"C\C=C/[C@@H](\C=C\O)[C@H](C)[C@H](\C=C/C)\C=C\O",
            vec![
                (3, StereoDescriptor::R),
                (7, StereoDescriptor::LowerR),
                (9, StereoDescriptor::S),
            ],
        ),
        (
            r"C\C=C/[C@@H](\C=C\C)[C@H](C)[C@H](\C=C/C)\C=C\C",
            vec![
                (3, StereoDescriptor::R),
                (7, StereoDescriptor::LowerR),
                (9, StereoDescriptor::S),
            ],
        ),
    ] {
        let mut molecule = read_smiles(smiles).expect("RDKit para-stereochemistry scaffold parses");
        perceive(&mut molecule).expect("RDKit para-stereochemistry scaffold perceives");

        assign_cip(&mut molecule);

        assert_eq!(tetrahedral_descriptor_map(&molecule), expected);
    }
}

#[test]
fn cip_matches_rdkit_for_auxiliary_stereochemistry_beyond_initial_expansion() {
    let mut molecule = read_smiles("CC[C@H](C)CCCCC[C@H]1CC[C@@H](C)CC1")
        .expect("RDKit auxiliary para-stereochemistry scaffold parses");
    perceive(&mut molecule).expect("RDKit auxiliary para-stereochemistry scaffold perceives");

    assign_cip(&mut molecule);

    assert_eq!(
        tetrahedral_descriptor_map(&molecule),
        vec![
            (2, StereoDescriptor::S),
            (9, StereoDescriptor::LowerS),
            (12, StereoDescriptor::S)
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_cyclohexane_pseudo_symmetry_examples() {
    for (smiles, expected) in [
        (
            "OC(=O)[C@H]1CC[C@@H](CC1)O[C@@H](F)Cl",
            vec![
                (3, StereoDescriptor::S),
                (6, StereoDescriptor::LowerR),
                (10, StereoDescriptor::S),
            ],
        ),
        (
            "OC(=O)[C@H]1CC[C@@H](CC1)O[CH](Cl)Cl",
            vec![(3, StereoDescriptor::LowerR), (6, StereoDescriptor::LowerR)],
        ),
    ] {
        let mut molecule = read_smiles(smiles).expect("RDKit pseudo-symmetry scaffold parses");
        perceive(&mut molecule).expect("RDKit pseudo-symmetry scaffold perceives");

        assign_cip(&mut molecule);

        assert_eq!(tetrahedral_descriptor_map(&molecule), expected);
    }
}

#[test]
fn cip_preserves_absolute_centers_next_to_pseudoasymmetric_ring_center() {
    let mut molecule = read_smiles("CCOC=1C=CC(=CC1OCC)C(C)N[C@H]2CC[C@]3(C[C@H]3C#N)CC2")
        .expect("mixed absolute and pseudoasymmetric scaffold parses");
    perceive(&mut molecule).expect("mixed scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![
            (15, StereoDescriptor::S),
            (18, StereoDescriptor::LowerS),
            (20, StereoDescriptor::R)
        ]
    );
}

#[test]
fn cip_bootstraps_coupled_pseudoasymmetric_fused_ring_centers() {
    let mut molecule = read_smiles("O=S(=O)(N[C@H]1C[C@H](C1)C2=NN=C3CCCCCN23)C4CC54CCC5")
        .expect("fused para-stereo scaffold parses");
    perceive(&mut molecule).expect("fused para-stereo scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![(4, StereoDescriptor::LowerS), (6, StereoDescriptor::LowerS)]
    );
}

#[test]
fn cip_bootstraps_coupled_pseudoasymmetric_cyclopentane_centers() {
    let mut molecule =
        read_smiles("CC=1N=CC(=CN1)C(=O)N[C@@H]2C[C@H](CNC(=O)C=3C=NC(=NC3)C(F)(F)F)C2")
            .expect("cyclopentane para-stereo scaffold parses");
    perceive(&mut molecule).expect("cyclopentane para-stereo scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![
            (10, StereoDescriptor::LowerS),
            (12, StereoDescriptor::LowerS)
        ]
    );
}

#[test]
fn cip_marks_middle_center_pseudoasymmetric_in_fused_three_center_system() {
    let mut molecule =
        read_smiles("CCC1(CCOCC1)C(=O)N2C[C@H]3[C@H](NC(=O)C4=CN(C)C(=O)C=N4)[C@H]3C2")
            .expect("three-center fused scaffold parses");
    perceive(&mut molecule).expect("three-center fused scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![
            (12, StereoDescriptor::S),
            (13, StereoDescriptor::LowerS),
            (25, StereoDescriptor::R)
        ]
    );
}

#[test]
fn cip_bootstraps_enamine_coupled_cyclobutane_pseudoasymmetric_centers() {
    let mut molecule =
        read_smiles("O=C(CCC(=O)N1CCC(=N1)C=2C=CC=CC2)N[C@@H]3C[C@H](C3)C4=CC=CC(=C4)C=5N=NNN5")
            .expect("Enamine coupled pseudoasymmetric scaffold parses");
    perceive(&mut molecule).expect("Enamine coupled scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![
            (18, StereoDescriptor::LowerR),
            (20, StereoDescriptor::LowerR)
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_enamine_quaternary_ring_center() {
    let mut molecule =
        read_smiles("C[C@]1(O)C[C@H](C1)C(=O)N2CC[C@H](CCNC(=O)C[C@@H]3CCCC[C@H]3O)C2")
            .expect("Enamine quaternary ring-center scaffold parses");
    perceive(&mut molecule).expect("Enamine quaternary scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![
            (1, StereoDescriptor::S),
            (4, StereoDescriptor::LowerS),
            (11, StereoDescriptor::S),
            (18, StereoDescriptor::S),
            (23, StereoDescriptor::R)
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_enamine_fused_three_center_pseudoasymmetry() {
    let mut molecule = read_smiles("CC1=CSC=C1C(=O)N2C[C@H]3[C@H](CNC(=O)CN4CCC(C)CC4)[C@H]3C2")
        .expect("Enamine fused three-center scaffold parses");
    perceive(&mut molecule).expect("Enamine fused three-center scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![
            (10, StereoDescriptor::S),
            (11, StereoDescriptor::LowerR),
            (24, StereoDescriptor::R)
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_enamine_fused_ring_dual_pseudoasymmetry() {
    let mut molecule =
        read_smiles("CC(C)(C)C(=O)N[C@H]1C[C@H]2C[C@H](C[C@H]2C1)NC(=O)C=3C=CC=CC3N4C=NC=N4")
            .expect("Enamine fused-ring dual pseudoasymmetric scaffold parses");
    perceive(&mut molecule).expect("Enamine fused-ring dual pseudoasymmetric scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![
            (7, StereoDescriptor::LowerR),
            (9, StereoDescriptor::S),
            (11, StereoDescriptor::LowerR),
            (13, StereoDescriptor::R)
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_enamine_spiro_fused_pseudoasymmetry() {
    let mut molecule = read_smiles("O=C(NS(=O)(=O)C=1C=NN(C1)C=2C=CC=CC2F)[C@@]34CCC[C@H]4CCC3")
        .expect("Enamine spiro-fused pseudoasymmetric scaffold parses");
    perceive(&mut molecule).expect("Enamine spiro-fused pseudoasymmetric scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![
            (18, StereoDescriptor::LowerR),
            (22, StereoDescriptor::LowerR)
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_enamine_absolute_center_in_coupled_bicycle() {
    let mut molecule = read_smiles("CN1N=CN=C1C=2C=CC(=CC2)C(=O)N3[C@H](C(=O)O)[C@@H]4CC[C@H]3CC4")
        .expect("Enamine coupled bicyclic scaffold parses");
    perceive(&mut molecule).expect("Enamine coupled bicyclic scaffold perceives");

    assign_cip(&mut molecule);

    let descriptors = tetrahedral_descriptor_map(&molecule);
    assert_eq!(
        descriptors,
        vec![
            (15, StereoDescriptor::S),
            (19, StereoDescriptor::R),
            (22, StereoDescriptor::R)
        ]
    );
}

#[test]
fn cip_applies_recursive_rule1a_before_isotope_priority() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let bromine = mol
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let mut carbon_13 = carbon();
    carbon_13.isotope = Some(13);
    let isotope_carbon = mol.add_atom(carbon_13).expect("atom identifier capacity");
    let substituted_carbon = mol.add_atom(carbon()).expect("atom identifier capacity");
    let iodine = mol
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");

    for carrier in [bromine, isotope_carbon, substituted_carbon] {
        mol.add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
    }
    mol.add_bond(substituted_carbon, iodine, BondOrder::Single)
        .expect("substituent bond");

    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(bromine),
                    StereoCarrier::Atom(isotope_carbon),
                    StereoCarrier::Atom(substituted_carbon),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("stereo element");

    mol.working_mut().set_implicit_hydrogens(center, 1);

    let report = assign_cip(mol.working_mut());

    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: stereo,
            descriptor: StereoDescriptor::R,
        }]
    );
}

#[test]
fn cip_matches_rdkit_for_pubchem_73056_recursive_rule_ordering() {
    let mut molecule =
        read_smiles("CC1=C(C(=O)O[C@@H](C1)[C@@H](C)[C@H]2CC[C@@H]3[C@@]2(CC[C@H]4[C@H]3C[C@@H]5[C@]6([C@@]4(C(=O)C=C[C@@H]6OC(=O)C)C)O5)C)COC(=O)C")
            .expect("CID 73056 parses");
    perceive(&mut molecule).expect("CID 73056 perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![
            StereoDescriptor::S,
            StereoDescriptor::S,
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::S,
            StereoDescriptor::S,
            StereoDescriptor::S,
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::R,
            StereoDescriptor::S,
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_pubchem_134556_recursive_rule_ordering() {
    let mut molecule = read_smiles("CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)[14CH2]O)O")
        .expect("CID 134556 parses");
    perceive(&mut molecule).expect("CID 134556 perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::R,
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_pubchem_246236_phosphorus_centers() {
    let mut molecule =
        read_smiles("C1COCCN1[P@]2(=NP(=N[P@@](=NP(=N2)(Cl)Cl)(N3CCOCC3)Cl)(Cl)Cl)Cl")
            .expect("CID 246236 parses");
    perceive(&mut molecule).expect("CID 246236 perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![StereoDescriptor::R, StereoDescriptor::S]
    );
}

#[test]
fn cip_matches_rdkit_for_pubchem_359164_sulfur_lone_pair() {
    let mut molecule =
        read_smiles("C1=CC=C(C=C1)N=NC2=CC3=C(C=C2)S[S@@](=O)N3").expect("CID 359164 parses");
    perceive(&mut molecule).expect("CID 359164 perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(assigned_descriptors(&report), vec![StereoDescriptor::R]);
}

#[test]
fn cip_matches_rdkit_for_pubchem_444295_with_spectators_interpreted_separately() {
    let mut molecule =
        read_smiles_component("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N.[NH2-].[NH2-].[NH2-].[NH2-].[NH2-].[OH3+].[OH3+].O.[Ac].[Ac].[Ac].[Ac].[Ac].[Ac].[Ac].[Ac].[Ac].[Ac]", 0)
            .expect("CID 444295 parses");
    perceive(&mut molecule).expect("CID 444295 perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![
            StereoDescriptor::R,
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::R,
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_pubchem_446291_with_unsupported_spectator_interpreted_separately() {
    let mut molecule =
        read_smiles_component("CCCCCCCCCCCCC(=O)CSCCNC(=O)CCNC(=O)[C@H](C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O.[Cf]", 0)
            .expect("CID 446291 parses");
    perceive(&mut molecule).expect("CID 446291 perceives");

    let report = assign_cip(&mut molecule);

    assert_eq!(
        assigned_descriptors(&report),
        vec![
            StereoDescriptor::S,
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::R,
            StereoDescriptor::R,
        ]
    );
}

#[test]
fn cip_skips_endocyclic_hetero_double_bond_stereo() {
    let mut molecule =
        read_smiles("C/C/1=C/2\\[C@@]([C@@H](/C(=C/C3=N/C(=C(\\C4=N[C@H]([C@@H]([C@@]4(C)CCC(=O)O)CC(=O)O)[C@]5([C@@]([C@@H](C1=N5)CCC(=O)O)(C)CC(=O)O)C)/C)/[C@@H](C3(C)C)CCC(=O)O)/N2)CCC(=O)O)(C)CC(=O)O")
            .expect("CID 446180 parses");
    perceive(&mut molecule).expect("CID 446180 perceives");

    assign_cip(&mut molecule);

    let bond_descriptors = double_bond_descriptor_map(&molecule);
    assert_eq!(
        bond_descriptors,
        vec![
            (1, 2, StereoDescriptor::Z),
            (5, 6, StereoDescriptor::Z),
            (9, 10, StereoDescriptor::Z),
        ]
    );
}

#[test]
fn cip_skips_equivalent_ligands_as_nonstereogenic() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let fluorine = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let chlorine = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let methyl_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let methyl_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    for carrier in [fluorine, chlorine, methyl_a, methyl_b] {
        mol.add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
    }
    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(fluorine),
                    StereoCarrier::Atom(chlorine),
                    StereoCarrier::Atom(methyl_a),
                    StereoCarrier::Atom(methyl_b),
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("stereo element");

    assert_cip_not_stereogenic(mol.working_mut(), stereo);
}

#[test]
fn cip_skips_large_complete_equivalent_ligands_as_nonstereogenic() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let fluorine = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let chlorine = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let chain_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let chain_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    for carrier in [fluorine, chlorine, chain_a, chain_b] {
        mol.add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
    }
    for chain in [chain_a, chain_b] {
        let mut previous = chain;
        for _ in 1..18 {
            let next = mol.add_atom(carbon()).expect("atom identifier capacity");
            mol.add_bond(previous, next, BondOrder::Single)
                .expect("chain bond");
            previous = next;
        }
    }
    assert!(mol.atom_count() > CipAssignmentOptions::default().max_depth);

    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(fluorine),
                    StereoCarrier::Atom(chlorine),
                    StereoCarrier::Atom(chain_a),
                    StereoCarrier::Atom(chain_b),
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("stereo element");

    assert_cip_not_stereogenic(mol.working_mut(), stereo);
}

#[test]
fn cip_skips_large_complete_equivalent_double_bond_endpoint_as_nonstereogenic() {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let double_bond = mol
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    let fluorine = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let chlorine = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let chain_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let chain_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    for carrier in [fluorine, chlorine] {
        mol.add_bond(left, carrier, BondOrder::Single)
            .expect("left carrier bond");
    }
    for carrier in [chain_a, chain_b] {
        mol.add_bond(right, carrier, BondOrder::Single)
            .expect("right carrier bond");
    }
    for chain in [chain_a, chain_b] {
        add_carbon_chain(mol.working_mut(), chain, 18);
    }
    assert!(mol.atom_count() > CipAssignmentOptions::default().max_depth);

    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: double_bond,
                left,
                right,
                left_carrier: StereoCarrier::Atom(fluorine),
                right_carrier: StereoCarrier::Atom(chain_a),
                orientation: Some(DoubleBondOrientation::Together),
            },
        )))
        .expect("double-bond stereo element");

    assert_cip_not_stereogenic(mol.working_mut(), stereo);
}

#[test]
fn cip_skips_large_complete_equivalent_axis_endpoint_as_nonstereogenic() {
    let mut mol = crate::core::MoleculeEditor::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let axis = mol.add_bond(left, right, BondOrder::Single).expect("axis");
    let iodine = mol
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");
    let bromine = mol
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let chain_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let chain_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    for carrier in [iodine, bromine] {
        mol.add_bond(left, carrier, BondOrder::Single)
            .expect("left carrier bond");
    }
    for carrier in [chain_a, chain_b] {
        mol.add_bond(right, carrier, BondOrder::Single)
            .expect("right carrier bond");
    }
    for chain in [chain_a, chain_b] {
        add_carbon_chain(mol.working_mut(), chain, 18);
    }
    assert!(mol.atom_count() > CipAssignmentOptions::default().max_depth);

    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![StereoCarrier::Atom(iodine), StereoCarrier::Atom(chain_a)],
            orientation: Some(AxisOrientation::CounterClockwise),
        })))
        .expect("axis stereo element");

    assert_cip_not_stereogenic(mol.working_mut(), stereo);
}

fn add_carbon_chain(mol: &mut Molecule, start: AtomId, length: usize) {
    let mut previous = start;
    for _ in 1..length {
        let next = mol.add_atom(carbon()).expect("atom identifier capacity");
        mol.add_bond(previous, next, BondOrder::Single)
            .expect("chain bond");
        previous = next;
    }
}

#[test]
fn cip_skips_equivalent_ring_ligands_as_nonstereogenic() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let nitrogen = mol
        .add_atom(element_atom("N"))
        .expect("atom identifier capacity");
    let hydrogen = mol
        .add_atom(element_atom("H"))
        .expect("atom identifier capacity");
    let ring_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let ring_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let ring_c = mol.add_atom(carbon()).expect("atom identifier capacity");
    let ring_d = mol.add_atom(carbon()).expect("atom identifier capacity");

    for carrier in [nitrogen, hydrogen, ring_a, ring_b] {
        mol.add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
    }
    mol.add_bond(ring_a, ring_c, BondOrder::Single)
        .expect("ring bond");
    mol.add_bond(ring_c, ring_d, BondOrder::Single)
        .expect("ring bond");
    mol.add_bond(ring_d, ring_b, BondOrder::Single)
        .expect("ring bond");

    let stereo = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(hydrogen),
                    StereoCarrier::Atom(nitrogen),
                    StereoCarrier::Atom(ring_a),
                    StereoCarrier::Atom(ring_b),
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("stereo element");

    assert_cip_not_stereogenic(mol.working_mut(), stereo);
}

#[test]
fn failed_cip_resource_limit_restores_previous_stereo_section() {
    let mut molecule = read_smiles("C[C@@H](C(=O)O)N").expect("alanine parses");
    perceive(&mut molecule).expect("alanine perceives");
    let report = assign_cip(&mut molecule);
    assert_eq!(report.assigned.len(), 1);
    let before = installed_cip_descriptors(&molecule);
    let before_state = molecule.perception().stereo_state().cloned();

    let error = stereo_api::assign_cip_descriptors_with_options(
        &mut molecule,
        CipAssignmentOptions {
            max_nodes: 1,
            ..CipAssignmentOptions::default()
        },
    )
    .expect_err("the constrained reassignment should fail");

    assert_eq!(
        error.issues,
        vec![CipAssignmentIssue::ResourceLimitExceeded {
            element: StereoElementId::new(0),
            max_nodes: 1,
        }]
    );
    assert_eq!(installed_cip_descriptors(&molecule), before);
    assert_eq!(molecule.perception().stereo_state(), before_state.as_ref());
    assert_eq!(
        molecule
            .cip_descriptor(StereoElementId::new(0))
            .expect("stereo element"),
        Some(StereoDescriptor::S)
    );

    let present_empty = Perception::builder()
        .with_cip_descriptors(Vec::new())
        .expect("empty stereo section")
        .build();
    molecule
        .install_perception(present_empty.clone())
        .expect("present-empty baseline");
    stereo_api::assign_cip_descriptors_with_options(
        &mut molecule,
        CipAssignmentOptions {
            max_nodes: 1,
            ..CipAssignmentOptions::default()
        },
    )
    .expect_err("the constrained reassignment should still fail");
    assert_eq!(molecule.perception(), &present_empty);
}

#[test]
fn failed_mixed_cip_attempt_publishes_no_partial_assignments() {
    const FIXTURE: &str = "F[C@](Cl)(Br)C[C@@H](N)C(=O)O";
    let options = CipAssignmentOptions {
        max_nodes: 1,
        ..CipAssignmentOptions::default()
    };

    let mut first_only = read_smiles(FIXTURE).expect("fixture parses");
    perceive(&mut first_only).expect("fixture perceives");
    first_only
        .remove_stereo_element(StereoElementId::new(1))
        .expect("second stereo element exists");
    let report = stereo_api::assign_cip_descriptors_with_options(&mut first_only, options)
        .expect("the first center is assignable within the resource limit");
    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::R,
        }]
    );

    let mut molecule = read_smiles(FIXTURE).expect("fixture parses");
    perceive(&mut molecule).expect("fixture perceives");
    let error = stereo_api::assign_cip_descriptors_with_options(&mut molecule, options)
        .expect_err("the later center should exceed the resource limit");
    assert_eq!(
        error.issues,
        vec![CipAssignmentIssue::ResourceLimitExceeded {
            element: StereoElementId::new(1),
            max_nodes: 1,
        }]
    );
    assert!(installed_cip_descriptors(&molecule).is_empty());
    assert!(!molecule.perception().has_stereo());
    assert!(molecule
        .stereo_elements()
        .all(|(element_id, _)| molecule.cip_descriptor(element_id) == Ok(None)));
}

#[test]
fn failed_cip_validation_preserves_previous_stereo_section() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let adjacent = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let nonadjacent = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(center, adjacent, BondOrder::Single)
        .expect("center bond");
    mol.add_bond(adjacent, nonadjacent, BondOrder::Single)
        .expect("connecting bond");
    let element = StereoElementId::new(0);
    mol.working_mut()
        .graph
        .stereo_elements
        .push(Some(StereoElement {
            kind: StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(adjacent),
                    StereoCarrier::Atom(adjacent),
                    StereoCarrier::Atom(adjacent),
                ],
                orientation: None,
            }),
            group: None,
        }));
    // Inject nonadjacency through internal storage to test defensive CIP
    // validation; checked editing now rejects this malformed reference.
    let stored = mol.working_mut().graph.stereo_elements[element.index()]
        .as_mut()
        .expect("stored stereo element");
    let StereoElementKind::Tetrahedral(stereo) = &mut stored.kind else {
        unreachable!("test element is tetrahedral");
    };
    stereo.carriers[2] = StereoCarrier::Atom(nonadjacent);
    mol.working_mut()
        .install_cip_descriptor(element, StereoDescriptor::R);
    let before = installed_cip_descriptors(mol.working());

    let error = stereo_api::assign_cip_descriptors(mol.working_mut())
        .expect_err("invalid stored stereo should reject CIP assignment");

    assert!(error
        .issues
        .iter()
        .all(|issue| matches!(issue, CipAssignmentIssue::InvalidStereo { .. })));
    assert_eq!(installed_cip_descriptors(mol.working()), before);
    assert_eq!(
        mol.cip_descriptor(element).expect("stereo element"),
        Some(StereoDescriptor::R)
    );

    let present_empty = Perception::builder()
        .with_cip_descriptors(Vec::new())
        .expect("empty stereo section")
        .build();
    mol.working_mut()
        .install_perception(present_empty.clone())
        .expect("valid empty stereo section");
    stereo_api::assign_cip_descriptors(mol.working_mut())
        .expect_err("invalid stored stereo should still reject CIP assignment");
    assert_eq!(mol.perception(), &present_empty);
    assert!(mol.perception().has_stereo());
}

#[test]
fn successful_cip_reassignment_replaces_the_complete_descriptor_set() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let carriers = ["F", "Cl", "Br", "I"]
        .into_iter()
        .map(element_atom)
        .map(|atom| mol.add_atom(atom).expect("atom identifier capacity"))
        .collect::<Vec<_>>();
    for carrier in &carriers {
        mol.add_bond(center, *carrier, BondOrder::Single)
            .expect("carrier bond");
    }
    let kind = || {
        StereoElementKind::Tetrahedral(TetrahedralStereo {
            center,
            carriers: carriers.iter().copied().map(StereoCarrier::Atom).collect(),
            orientation: Some(TetrahedralOrientation::Clockwise),
        })
    };
    let specified = mol
        .add_stereo_element(StereoElement::new(kind()))
        .expect("specified stereo element");
    let mut unknown_kind = kind();
    let StereoElementKind::Tetrahedral(stereo) = &mut unknown_kind else {
        unreachable!("test fixture is tetrahedral");
    };
    // Keep the unknown assertion on a separate center: each graph focus has
    // exactly one represented configuration, including an unknown one.
    let second_center = mol.add_atom(carbon()).unwrap();
    let mut second_carriers = Vec::new();
    for symbol in ["F", "Cl", "Br", "I"] {
        let carrier = mol.add_atom(element_atom(symbol)).unwrap();
        mol.add_bond(second_center, carrier, BondOrder::Single)
            .unwrap();
        second_carriers.push(StereoCarrier::Atom(carrier));
    }
    stereo.center = second_center;
    stereo.carriers = second_carriers;
    stereo.orientation = None;
    let unknown = mol
        .add_stereo_element(StereoElement::new(unknown_kind))
        .expect("unknown stereo element");
    mol.working_mut()
        .install_cip_descriptor(specified, StereoDescriptor::S);
    mol.working_mut()
        .install_cip_descriptor(unknown, StereoDescriptor::R);

    let report = assign_cip(mol.working_mut());

    assert_eq!(report.assigned.len(), 1);
    assert_eq!(report.assigned[0].element, specified);
    assert_eq!(
        report.skipped,
        vec![CipSkipped {
            element: unknown,
            reason: CipSkippedReason::UnknownConfiguration,
        }]
    );
    assert_eq!(
        installed_cip_descriptors(mol.working()),
        vec![(specified, report.assigned[0].descriptor)]
    );
    assert_eq!(mol.cip_descriptor(unknown).expect("unknown element"), None);
}

#[test]
fn cip_descriptors_are_cleared_by_stereo_invalidating_mutations() {
    let mut molecule = read_smiles("C[C@@H](C(=O)O)N").expect("alanine parses");
    perceive(&mut molecule).expect("alanine perceives");
    assign_cip(&mut molecule);
    assert_eq!(
        molecule
            .cip_descriptor(StereoElementId::new(0))
            .expect("stereo element"),
        Some(StereoDescriptor::S)
    );

    molecule
        .add_atom(oxygen())
        .expect("atom identifier capacity");

    assert_eq!(
        molecule
            .perception()
            .cip_descriptor(StereoElementId::new(0)),
        None
    );
}

#[test]
fn cip_rule6_descriptor_case_is_invariant_to_equivalent_smiles_neighbor_order() {
    // Published Hanson validation suite VS294-298 specifies absolute R/S.
    // These reordered encodings preserve the complete chiral graph. RDKit
    // 2026.03.5 changes only the central label's case for some encodings;
    // retain the published descriptors rather than that traversal dependence.
    use StereoDescriptor::{R, S};
    for (smiles, expected) in [
        (
            "O[C@H]1O[C@@]2(O[C@H]1O)O[C@H]([C@H](O2)O)O",
            vec![(1, S), (3, R), (5, R), (8, R), (9, S)],
        ),
        (
            "[C@]12(O[C@@H]([C@H](O)O1)O)O[C@H]([C@H](O2)O)O",
            vec![(0, R), (2, S), (3, R), (8, R), (9, S)],
        ),
        (
            "O[C@H]1O[C@@]2(O[C@H]1O)O[C@@H]([C@@H](O2)O)O",
            vec![(1, S), (3, S), (5, R), (8, S), (9, R)],
        ),
        (
            "O1[C@H]([C@H](O[C@]21O[C@H]([C@H](O2)O)O)O)O",
            vec![(1, R), (2, S), (4, S), (6, R), (7, S)],
        ),
        (
            "[C@@H]1(O)[C@@H](O[C@@]2(O[C@H](O)[C@H](O)O2)O1)O",
            vec![(0, S), (2, R), (4, S), (6, S), (8, R)],
        ),
        (
            "C1[C@@]2(OCCC1)CC[C@]3(CC2)CC[C@@]4(CC3)OCCCC4",
            vec![(1, S), (8, S), (13, S)],
        ),
        (
            "C1[C@]2(CCCCO2)CC[C@@]2(CC[C@]3(CC2)OCCCC3)C1",
            vec![(1, S), (9, S), (12, S)],
        ),
        (
            "C1[C@@]2(CC[C@]3(C1)OCCC3)CC[C@@]4(CC2)CCCO4",
            vec![(1, R), (4, S), (12, S)],
        ),
        (
            "C1[C@]2(CC[C@]3(CC[C@@]4(CC3)CCCO4)C1)OCCC2",
            vec![(1, S), (4, R), (7, S)],
        ),
        (
            "[C@]12(CC[C@]3(OCCC3)CC1)CC[C@]1(CC2)OCCC1",
            vec![(0, R), (3, S), (12, S)],
        ),
        (
            "Cl[C@H]1C[C@]2(C1)C[C@H](C2)Cl",
            vec![(1, R), (3, S), (6, R)],
        ),
        (
            "C1[C@@]2(C[C@H](Cl)C2)C[C@@H]1Cl",
            vec![(1, S), (3, R), (7, R)],
        ),
        (
            "Cl[C@H]1C[C@@]2(C[C@@H](Cl)C2)C1",
            vec![(1, R), (3, S), (5, R)],
        ),
    ] {
        let mut molecule = read_smiles(smiles).expect("validation molecule parses");
        perceive(&mut molecule).expect("validation molecule perceives");
        assign_cip(&mut molecule);
        assert_eq!(tetrahedral_descriptor_map(&molecule), expected, "{smiles}");
    }
}

#[test]
fn cip_assigns_mutually_dependent_double_bond_auxiliary_descriptors() {
    // Hanson validation suite VS188: each alkene differentiates the two
    // otherwise equivalent ring paths used to label the other alkene.
    for (smiles, expected) in [
        (
            r"C\C=C/1\C/C(/C1)=C\C",
            vec![(1, 2, StereoDescriptor::E), (4, 6, StereoDescriptor::E)],
        ),
        (
            r"C/C=C1\C/C(=C/C)C1",
            vec![(1, 2, StereoDescriptor::E), (4, 5, StereoDescriptor::E)],
        ),
    ] {
        let mut molecule = read_smiles(smiles).expect("VS188 parses");
        perceive(&mut molecule).expect("VS188 perceives");
        assign_cip(&mut molecule);
        assert_eq!(double_bond_descriptor_map(&molecule), expected, "{smiles}");
    }
}

#[test]
fn cip_uses_locally_rooted_alkene_auxiliary_labels_at_tetrahedral_centers() {
    // Hanson validation suite VS191, VS214 and VS246 exercise ordinary,
    // pseudoasymmetric and sequence-cis/trans auxiliary labels respectively.
    for (smiles, atoms, bonds) in [
        (
            r"C/C=C\C/C(=C\[C@H](C=C(C/C=C\C)C/C=C\C)O)/C/C=C/C",
            vec![(6, StereoDescriptor::S)],
            vec![
                (1, 2, StereoDescriptor::Z),
                (4, 5, StereoDescriptor::Z),
                (10, 11, StereoDescriptor::Z),
                (14, 15, StereoDescriptor::Z),
                (19, 20, StereoDescriptor::E),
            ],
        ),
        (
            r"CC[C@@H]1CC/C(/CC1)=C\2/CC[C@H](CC2)CC",
            vec![
                (2, StereoDescriptor::LowerS),
                (11, StereoDescriptor::LowerS),
            ],
            vec![(5, 8, StereoDescriptor::E)],
        ),
        (
            r"C/C=C/1\CC[C@H](C)CC1",
            vec![(5, StereoDescriptor::R)],
            vec![(1, 2, StereoDescriptor::SeqCis)],
        ),
    ] {
        let mut molecule = read_smiles(smiles).expect("validation molecule parses");
        perceive(&mut molecule).expect("validation molecule perceives");
        assign_cip(&mut molecule);
        assert_eq!(tetrahedral_descriptor_map(&molecule), atoms, "{smiles}");
        assert_eq!(double_bond_descriptor_map(&molecule), bonds, "{smiles}");
    }
}

#[test]
fn cip_progressive_expansion_labels_large_pubchem_cage_with_default_limits() {
    // External PubChem100k pack_038.smi:223. Eager complete expansion exhausted
    // 100,000 nodes even though Rule 1a resolves every center near its root.
    let mut molecule = read_smiles("C[C@H]1[C@H]2CC=C3C2(COC14[C@@H](CC(O4)(C)C)O)C(=O)C[C@H]5[C@H]3CC[C@@H]6[C@@]5(CC7=NC8=C(C[C@]9([C@H](C8)CC[C@@H]1[C@@H]9C[C@H]([C@]2(C1=C[C@H]1[C@@]2([C@@H](C2(O1)[C@@H](C[C@@](O2)(C)CO)O)C)O)C)O)C)N=C7C6)C").expect("PubChem cage parses");
    perceive(&mut molecule).expect("PubChem cage perceives");
    assign_cip(&mut molecule);
    assert_eq!(
        tetrahedral_descriptor_map(&molecule),
        vec![
            (1, StereoDescriptor::S),
            (2, StereoDescriptor::R),
            (10, StereoDescriptor::R),
            (20, StereoDescriptor::S),
            (21, StereoDescriptor::R),
            (24, StereoDescriptor::S),
            (25, StereoDescriptor::S),
            (32, StereoDescriptor::S),
            (33, StereoDescriptor::S),
            (37, StereoDescriptor::R),
            (38, StereoDescriptor::S),
            (40, StereoDescriptor::R),
            (41, StereoDescriptor::R),
            (44, StereoDescriptor::S),
            (45, StereoDescriptor::S),
            (46, StereoDescriptor::S),
            (49, StereoDescriptor::R),
            (51, StereoDescriptor::S),
        ]
    );
}

fn tetrahedral_descriptor_map(mol: &Molecule) -> Vec<(u32, StereoDescriptor)> {
    mol.stereo_elements()
        .filter_map(|(element_id, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => mol
                .cip_descriptor(element_id)
                .expect("stereo element")
                .map(|descriptor| (stereo.center.raw(), descriptor)),
            StereoElementKind::DoubleBond(_) | StereoElementKind::Axis(_) => None,
        })
        .collect()
}

fn double_bond_descriptor_map(mol: &Molecule) -> Vec<(u32, u32, StereoDescriptor)> {
    mol.stereo_elements()
        .filter_map(|(element_id, element)| match &element.kind {
            StereoElementKind::DoubleBond(stereo) => mol
                .cip_descriptor(element_id)
                .expect("stereo element")
                .map(|descriptor| (stereo.left.raw(), stereo.right.raw(), descriptor)),
            StereoElementKind::Tetrahedral(_) | StereoElementKind::Axis(_) => None,
        })
        .collect()
}
