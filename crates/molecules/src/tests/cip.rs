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
fn cip_matches_rdkit_for_pubchem_start_atom_bracket_h_tetrahedral_centers() {
    let mut molecule =
        smiles_api::read_str("[C@@H]([C@H](C(=O)O)O)(C(=O)O)O").expect("tartrate parses");
    perception_api::sanitize(&mut molecule).expect("tartrate sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
        vec![StereoDescriptor::R, StereoDescriptor::R]
    );
}

#[test]
fn cip_matches_rdkit_for_smiles_ring_digit_tetrahedral_order() {
    let mut molecule = smiles_api::read_str("CC(C)C[C@@H]1CN2CCC3=CC(=C(C=C3C2CC1=O)OC)O[11CH3]")
        .expect("ring chiral molecule parses");
    perception_api::sanitize(&mut molecule).expect("ring chiral molecule sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
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
        smiles_api::read_str("C1=C2C(=NC=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O")
            .expect("nucleotide sugar parses");
    perception_api::sanitize(&mut molecule).expect("nucleotide sugar sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
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
        smiles_api::read_str("CC(=O)OC[C@]1([C@@H](CC[C@@]2(C1C[C@@H]([C@]34[C@H]2CC[C@@H](C3)C(=C)C4)OC(=O)C5=CC=C(C=C5)OC)C)OC(=O)C6=CC=C(C=C6)OC)C")
            .expect("polycycle parses");
    perception_api::sanitize(&mut molecule).expect("polycycle sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
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
        smiles_api::read_str("CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(C)C=C)/C").expect("polyene parses");
    perception_api::sanitize(&mut molecule).expect("polyene sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
        vec![
            StereoDescriptor::E,
            StereoDescriptor::E,
            StereoDescriptor::E
        ]
    );
}

#[test]
fn cip_skips_endocyclic_kekule_bond_stereo_after_ring_perception() {
    let mut molecule =
        smiles_api::read_str("CC\\1=C(/C/2=C/C3=C(C(=C(N3)/C=C\\4/[C@@](C(=C(N4)/C=C\\5/[C@@](C(=C(N5)/C=C1\\N2)O)(C)CC(=O)O)O)(C)CC(=O)O)C)CCC(=O)O)CCC(=O)O")
            .expect("CID 445170 parses");
    perception_api::sanitize(&mut molecule).expect("CID 445170 sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    let bond_descriptors = molecule
        .graph()
        .stereo_elements()
        .filter_map(|(_, element)| match &element.kind {
            StereoElementKind::DoubleBond(stereo) => element
                .descriptor
                .map(|descriptor| (stereo.left.raw(), stereo.right.raw(), descriptor)),
            StereoElementKind::Tetrahedral(_) | StereoElementKind::Axis(_) => None,
        })
        .collect::<Vec<_>>();
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
        smiles_api::read_str("CN1CC[C@@]23[C@H]4[C@H]1CC5=C2C(=C(C=C5)OC)O[C@@H]3[C@]6(C4)C(=O)C7=C8N6CCC9=C8C(=C(C=C9)OC)OC1=C7C=CC(=C1O)OC")
            .expect("fused ring parses");
    perception_api::sanitize(&mut molecule).expect("fused ring sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
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
fn cip_uses_rule3_embedded_e_z_descriptors_to_order_ligands() {
    let mut molecule =
        smiles_api::read_str("Br[C@H](/C=C/F)/C=C\\F").expect("Rule 3 alkene pair parses");
    perception_api::sanitize(&mut molecule).expect("Rule 3 alkene pair sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    let atom_descriptors = molecule
        .graph()
        .stereo_elements()
        .filter_map(|(_, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => element
                .descriptor
                .map(|descriptor| (stereo.center.raw(), descriptor)),
            StereoElementKind::DoubleBond(_) | StereoElementKind::Axis(_) => None,
        })
        .collect::<Vec<_>>();
    let bond_descriptors = molecule
        .graph()
        .stereo_elements()
        .filter_map(|(_, element)| match &element.kind {
            StereoElementKind::DoubleBond(stereo) => element
                .descriptor
                .map(|descriptor| (stereo.left.raw(), stereo.right.raw(), descriptor)),
            StereoElementKind::Tetrahedral(_) | StereoElementKind::Axis(_) => None,
        })
        .collect::<Vec<_>>();

    assert_eq!(atom_descriptors, vec![(1, StereoDescriptor::R)]);
    assert_eq!(
        bond_descriptors,
        vec![(2, 3, StereoDescriptor::E), (5, 6, StereoDescriptor::Z)]
    );
}

#[test]
fn cip_applies_recursive_rule1a_before_isotope_priority() {
    let mut mol = Molecule::new();
    let mut center_atom = carbon();
    center_atom.implicit_hydrogens = Some(1);
    let center = mol.add_atom(center_atom);
    let bromine = mol.add_atom(element_atom("Br"));
    let mut carbon_13 = carbon();
    carbon_13.isotope = Some(13);
    let isotope_carbon = mol.add_atom(carbon_13);
    let substituted_carbon = mol.add_atom(carbon());
    let iodine = mol.add_atom(element_atom("I"));

    for carrier in [bromine, isotope_carbon, substituted_carbon] {
        mol.add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
    }
    mol.add_bond(substituted_carbon, iodine, BondOrder::Single)
        .expect("substituent bond");

    let stereo = mol
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(bromine),
                    StereoCarrier::Atom(isotope_carbon),
                    StereoCarrier::Atom(substituted_carbon),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: TetrahedralOrientation::Clockwise,
            }),
            StereoSource::User,
        ))
        .expect("stereo element");

    let report = stereo_api::assign_cip_descriptors(&mut mol);

    assert!(report.is_ok(), "{:?}", report.issues);
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
        smiles_api::read_str("CC1=C(C(=O)O[C@@H](C1)[C@@H](C)[C@H]2CC[C@@H]3[C@@]2(CC[C@H]4[C@H]3C[C@@H]5[C@]6([C@@]4(C(=O)C=C[C@@H]6OC(=O)C)C)O5)C)COC(=O)C")
            .expect("CID 73056 parses");
    perception_api::sanitize(&mut molecule).expect("CID 73056 sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
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
    let mut molecule = smiles_api::read_str("CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)[14CH2]O)O")
        .expect("CID 134556 parses");
    perception_api::sanitize(&mut molecule).expect("CID 134556 sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
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
        smiles_api::read_str("C1COCCN1[P@]2(=NP(=N[P@@](=NP(=N2)(Cl)Cl)(N3CCOCC3)Cl)(Cl)Cl)Cl")
            .expect("CID 246236 parses");
    perception_api::sanitize(&mut molecule).expect("CID 246236 sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
        vec![StereoDescriptor::R, StereoDescriptor::S]
    );
}

#[test]
fn cip_matches_rdkit_for_pubchem_359164_sulfur_lone_pair() {
    let mut molecule = smiles_api::read_str("C1=CC=C(C=C1)N=NC2=CC3=C(C=C2)S[S@@](=O)N3")
        .expect("CID 359164 parses");
    perception_api::sanitize(&mut molecule).expect("CID 359164 sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
        vec![StereoDescriptor::R]
    );
}

#[test]
fn cip_matches_rdkit_for_pubchem_444295_disconnected_metal_fragments() {
    let mut molecule =
        smiles_api::read_str("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N.[NH2-].[NH2-].[NH2-].[NH2-].[NH2-].[OH3+].[OH3+].O.[Ac].[Ac].[Ac].[Ac].[Ac].[Ac].[Ac].[Ac].[Ac].[Ac]")
            .expect("CID 444295 parses");
    perception_api::sanitize(&mut molecule).expect("CID 444295 sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
        vec![
            StereoDescriptor::R,
            StereoDescriptor::R,
            StereoDescriptor::S,
            StereoDescriptor::R,
        ]
    );
}

#[test]
fn cip_matches_rdkit_for_pubchem_446291_disconnected_unsupported_spectator() {
    let mut molecule =
        smiles_api::read_str("CCCCCCCCCCCCC(=O)CSCCNC(=O)CCNC(=O)[C@H](C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O.[Cf]")
            .expect("CID 446291 parses");
    perception_api::sanitize(&mut molecule).expect("CID 446291 sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        report
            .assigned
            .iter()
            .map(|assignment| assignment.descriptor)
            .collect::<Vec<_>>(),
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
        smiles_api::read_str("C/C/1=C/2\\[C@@]([C@@H](/C(=C/C3=N/C(=C(\\C4=N[C@H]([C@@H]([C@@]4(C)CCC(=O)O)CC(=O)O)[C@]5([C@@]([C@@H](C1=N5)CCC(=O)O)(C)CC(=O)O)C)/C)/[C@@H](C3(C)C)CCC(=O)O)/N2)CCC(=O)O)(C)CC(=O)O")
            .expect("CID 446180 parses");
    perception_api::sanitize(&mut molecule).expect("CID 446180 sanitizes");

    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    let bond_descriptors = molecule
        .graph()
        .stereo_elements()
        .filter_map(|(_, element)| match &element.kind {
            StereoElementKind::DoubleBond(stereo) => element
                .descriptor
                .map(|descriptor| (stereo.left.raw(), stereo.right.raw(), descriptor)),
            StereoElementKind::Tetrahedral(_) | StereoElementKind::Axis(_) => None,
        })
        .collect::<Vec<_>>();
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
