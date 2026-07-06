use super::*;

#[test]
fn cip_assigns_tetrahedral_descriptors_from_stored_local_stereo() {
    let mut s_alanine = smiles_api::read_str("C[C@@H](C(=O)O)N").expect("alanine parses");
    perception_api::sanitize(&mut s_alanine).expect("alanine sanitizes");

    let report = stereo_api::assign_cip_descriptors(s_alanine.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::S,
        }]
    );
    assert_eq!(
        s_alanine
            .graph()
            .stereo_element(StereoElementId::new(0))
            .expect("stereo element")
            .descriptor,
        Some(StereoDescriptor::S)
    );

    let mut r_alanine = smiles_api::read_str("C[C@H](C(=O)O)N").expect("alanine parses");
    perception_api::sanitize(&mut r_alanine).expect("alanine sanitizes");

    let report = stereo_api::assign_cip_descriptors(r_alanine.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.assigned[0].descriptor, StereoDescriptor::R);
}

#[test]
fn cip_assigns_double_bond_descriptors_from_ranked_carriers() {
    let mut together = smiles_api::read_str("C(=C\\F)\\F").expect("alkene parses");
    perception_api::sanitize(&mut together).expect("alkene sanitizes");

    let report = stereo_api::assign_cip_descriptors(together.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report.assigned,
        vec![CipAssignment {
            element: StereoElementId::new(0),
            descriptor: StereoDescriptor::Z,
        }]
    );

    let mut opposite = smiles_api::read_str("C(=C/F)\\F").expect("alkene parses");
    perception_api::sanitize(&mut opposite).expect("alkene sanitizes");

    let report = stereo_api::assign_cip_descriptors(opposite.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.assigned[0].descriptor, StereoDescriptor::E);
}

#[test]
fn cip_reports_unresolved_equivalent_ligands_without_descriptor() {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon());
    let fluorine = mol.add_atom(element_atom("F"));
    let chlorine = mol.add_atom(element_atom("Cl"));
    let methyl_a = mol.add_atom(carbon());
    let methyl_b = mol.add_atom(carbon());
    for carrier in [fluorine, chlorine, methyl_a, methyl_b] {
        mol.add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
    }
    let stereo = mol
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(fluorine),
                    StereoCarrier::Atom(chlorine),
                    StereoCarrier::Atom(methyl_a),
                    StereoCarrier::Atom(methyl_b),
                ],
                orientation: TetrahedralOrientation::Clockwise,
            }),
            StereoSource::User,
        ))
        .expect("stereo element");

    let report = stereo_api::assign_cip_descriptors(&mut mol);

    assert_eq!(
        report.issues,
        vec![CipAssignmentIssue::UnresolvedPriority { element: stereo }]
    );
    assert!(report.assigned.is_empty());
    assert_eq!(
        mol.stereo_element(stereo).expect("element").descriptor,
        None
    );
}

#[test]
fn cip_respects_resource_limits_without_assigning_partial_descriptors() {
    let mut molecule = smiles_api::read_str("C[C@@H](C(=O)O)N").expect("alanine parses");
    perception_api::sanitize(&mut molecule).expect("alanine sanitizes");

    let report = stereo_api::assign_cip_descriptors_with_options(
        molecule.graph_mut(),
        CipAssignmentOptions {
            max_nodes: 1,
            ..CipAssignmentOptions::default()
        },
    );

    assert_eq!(
        report.issues,
        vec![CipAssignmentIssue::ResourceLimitExceeded {
            element: StereoElementId::new(0),
            max_nodes: 1,
        }]
    );
    assert!(report.assigned.is_empty());
    assert_eq!(
        molecule
            .graph()
            .stereo_element(StereoElementId::new(0))
            .expect("stereo element")
            .descriptor,
        None
    );
}

#[test]
fn cip_descriptors_are_cleared_by_stereo_invalidating_mutations() {
    let mut molecule = smiles_api::read_str("C[C@@H](C(=O)O)N").expect("alanine parses");
    perception_api::sanitize(&mut molecule).expect("alanine sanitizes");
    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());
    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        molecule
            .graph()
            .stereo_element(StereoElementId::new(0))
            .expect("stereo element")
            .descriptor,
        Some(StereoDescriptor::S)
    );

    molecule.graph_mut().add_atom(oxygen());

    assert_eq!(
        molecule
            .graph()
            .stereo_element(StereoElementId::new(0))
            .expect("stereo element")
            .descriptor,
        None
    );
}
