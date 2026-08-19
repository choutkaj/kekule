use super::*;

#[test]
fn element_from_atomic_number_accepts_periodic_table_bounds() {
    assert_eq!(
        Element::from_atomic_number(1)
            .expect("hydrogen exists")
            .symbol(),
        "H"
    );
    assert_eq!(
        Element::from_atomic_number(118)
            .expect("oganesson exists")
            .symbol(),
        "Og"
    );
}

#[test]
fn element_from_atomic_number_rejects_out_of_range_values() {
    assert_eq!(Element::from_atomic_number(0), None);
    assert_eq!(Element::from_atomic_number(119), None);
}

#[test]
fn element_from_symbol_is_canonical_and_case_sensitive() {
    assert_eq!(
        Element::from_symbol("C")
            .expect("carbon exists")
            .atomic_number(),
        6
    );
    assert_eq!(
        Element::from_symbol("Cl")
            .expect("chlorine exists")
            .atomic_number(),
        17
    );
    assert_eq!(
        Element::from_symbol("Og")
            .expect("oganesson exists")
            .atomic_number(),
        118
    );
    assert_eq!(Element::from_symbol("CL"), None);
    assert_eq!(Element::from_symbol("Xx"), None);
    assert_eq!(Element::from_symbol("?"), None);
}

#[test]
fn element_symbol_and_display_are_canonical() {
    let iron = Element::from_atomic_number(26).expect("iron exists");

    assert_eq!(iron.symbol(), "Fe");
    assert_eq!(iron.to_string(), "Fe");
}

#[test]
fn element_exposes_foundational_covalent_radii() {
    let hydrogen = Element::from_symbol("H").expect("hydrogen");
    let carbon = Element::from_symbol("C").expect("carbon");
    let curium = Element::from_symbol("Cm").expect("curium");
    let oganesson = Element::from_symbol("Og").expect("oganesson");

    assert_eq!(hydrogen.covalent_radius_angstrom(), Some(0.31));
    assert_eq!(carbon.covalent_radius_angstrom(), Some(0.76));
    assert_eq!(curium.covalent_radius_angstrom(), Some(1.69));
    assert_eq!(oganesson.covalent_radius_angstrom(), None);
}

#[test]
fn atom_new_sets_chemically_general_defaults() {
    let atom = carbon();

    assert_eq!(atom.element.symbol(), "C");
    assert_eq!(atom.isotope, None);
    assert_eq!(atom.formal_charge, 0);
    assert_eq!(atom.radical, None);
    assert_eq!(atom.hydrogens, HydrogenDeclaration::Infer { explicit: 0 });
    assert_eq!(atom.atom_map, None);
    assert!(atom.props.is_empty());
}

#[test]
fn atom_payload_fields_can_be_set_and_read() {
    let mut atom = carbon();
    atom.isotope = Some(13);
    atom.formal_charge = -1;
    atom.radical = Some(AtomRadical::Doublet);
    atom.hydrogens = HydrogenDeclaration::Fixed(3);
    atom.atom_map = Some(7);
    atom.props
        .insert("label".to_owned(), PropValue::String("alpha".to_owned()));

    assert_eq!(atom.isotope, Some(13));
    assert_eq!(atom.formal_charge, -1);
    assert_eq!(atom.radical, Some(AtomRadical::Doublet));
    assert_eq!(atom.hydrogens, HydrogenDeclaration::Fixed(3));
    assert_eq!(atom.atom_map, Some(7));
    assert_eq!(
        atom.props.get("label"),
        Some(&PropValue::String("alpha".to_owned()))
    );
}

#[test]
fn hydrogen_declaration_expresses_each_canonical_policy_without_overlap() {
    for (declaration, explicit, allows_implicit) in [
        (HydrogenDeclaration::Infer { explicit: 0 }, 0, true),
        (HydrogenDeclaration::Infer { explicit: 2 }, 2, true),
        (HydrogenDeclaration::Fixed(0), 0, false),
        (HydrogenDeclaration::Fixed(3), 3, false),
    ] {
        assert_eq!(declaration.explicit_count(), explicit);
        assert_eq!(declaration.allows_implicit(), allows_implicit);
        assert_eq!(
            declaration.with_explicit_count(7),
            if allows_implicit {
                HydrogenDeclaration::Infer { explicit: 7 }
            } else {
                HydrogenDeclaration::Fixed(7)
            }
        );
    }
}

#[test]
fn radical_multiplicity_reports_unpaired_electrons() {
    assert_eq!(AtomRadical::Singlet.unpaired_electron_count(), 0);
    assert_eq!(AtomRadical::Doublet.unpaired_electron_count(), 1);
    assert_eq!(AtomRadical::Triplet.unpaired_electron_count(), 2);
    assert_eq!(AtomRadical::Quartet.unpaired_electron_count(), 3);
    assert_eq!(AtomRadical::Quintet.unpaired_electron_count(), 4);
}

#[test]
fn bond_new_sets_endpoints_and_order() {
    let a = AtomId::new(3);
    let b = AtomId::new(4);
    let single = Bond::new(a, b, BondOrder::Single);
    let double = Bond::new(a, b, BondOrder::Double);

    assert_eq!(single.a(), a);
    assert_eq!(single.b(), b);
    assert_eq!(single.endpoints(), (a, b));
    assert_eq!(single.order, BondOrder::Single);
    assert!(single.props.is_empty());
    assert_eq!(double.order, BondOrder::Double);
}

#[test]
fn bond_payload_fields_can_be_set_and_read() {
    let mut bond = Bond::new(AtomId::new(1), AtomId::new(2), BondOrder::Dative);
    bond.props
        .insert("score".to_owned(), PropValue::Float(1.25));

    assert_eq!(bond.order, BondOrder::Dative);
    assert_eq!(bond.props.get("score"), Some(&PropValue::Float(1.25)));
}

#[test]
fn stereo_elements_and_groups_live_on_molecule() {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let a = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(center, a, BondOrder::Single).expect("bond");
    mol.add_bond(center, b, BondOrder::Single).expect("bond");
    mol.add_bond(center, c, BondOrder::Single).expect("bond");
    mark_all_fresh(&mut mol);

    let element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(a),
                    StereoCarrier::Atom(b),
                    StereoCarrier::Atom(c),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("stereo element should be stored");
    assert!(!mol.perception().has_stereo());

    let stored = mol.stereo_element(element).expect("stored element");
    assert!(stored.is_specified());

    let group = mol
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: vec![element],
        })
        .expect("group should be stored");
    assert_eq!(
        mol.stereo_element(element).expect("element").group,
        Some(group)
    );
    assert_eq!(
        mol.stereo_group(group).expect("group").members,
        vec![element]
    );
}

#[test]
fn stereo_replacement_and_group_creation_preserve_graph_references() {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let a = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    let element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(a),
                    StereoCarrier::Atom(b),
                    StereoCarrier::Atom(c),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("valid stereo element");
    let before = mol.stereo_element(element).expect("element").clone();
    let mut invalid = before.clone();
    let StereoElementKind::Tetrahedral(stereo) = &mut invalid.kind else {
        unreachable!("test element is tetrahedral");
    };
    stereo.center = AtomId::new(999);

    assert!(matches!(
        mol.replace_stereo_element(element, invalid),
        Err(MoleculeError::InvalidAtomId(id)) if id == AtomId::new(999)
    ));
    assert_eq!(mol.stereo_element(element).expect("element"), &before);

    assert!(matches!(
        mol.add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: Vec::new(),
        }),
        Err(MoleculeError::InvalidStereoReference(_))
    ));
    assert!(matches!(
        mol.add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: vec![element, element],
        }),
        Err(MoleculeError::InvalidStereoReference(_))
    ));
    assert!(mol.stereo_groups().next().is_none());
    assert_eq!(mol.stereo_element(element).expect("element").group, None);
}

#[test]
fn tetrahedral_stereo_storage_canonicalizes_carrier_permutations() {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let carriers = ["F", "Cl", "Br", "I"]
        .into_iter()
        .map(element_atom)
        .map(|atom| mol.add_atom(atom).expect("atom identifier capacity"))
        .collect::<Vec<_>>();
    for carrier in &carriers {
        mol.add_bond(center, *carrier, BondOrder::Single)
            .expect("tetrahedral carrier bond");
    }
    let canonical_carriers = carriers
        .iter()
        .copied()
        .map(StereoCarrier::Atom)
        .collect::<Vec<_>>();

    let canonical = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: canonical_carriers.clone(),
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("canonical tetrahedral element");
    let permuted = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    canonical_carriers[1],
                    canonical_carriers[0],
                    canonical_carriers[2],
                    canonical_carriers[3],
                ],
                orientation: Some(TetrahedralOrientation::CounterClockwise),
            },
        )))
        .expect("equivalent permuted tetrahedral element");
    assert_eq!(
        mol.stereo_element(canonical).unwrap(),
        mol.stereo_element(permuted).unwrap()
    );

    let unknown = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: canonical_carriers.clone(),
                orientation: None,
            },
        )))
        .expect("canonical unknown tetrahedral element");
    let unknown_permuted = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    canonical_carriers[3],
                    canonical_carriers[1],
                    canonical_carriers[0],
                    canonical_carriers[2],
                ],
                orientation: None,
            },
        )))
        .expect("permuted unknown tetrahedral element");
    assert_eq!(
        mol.stereo_element(unknown).unwrap(),
        mol.stereo_element(unknown_permuted).unwrap()
    );

    let replacement = StereoElement::new(StereoElementKind::Tetrahedral(TetrahedralStereo {
        center,
        carriers: vec![
            canonical_carriers[2],
            canonical_carriers[1],
            canonical_carriers[0],
            canonical_carriers[3],
        ],
        orientation: Some(TetrahedralOrientation::CounterClockwise),
    }));
    mol.replace_stereo_element(permuted, replacement)
        .expect("replacement should use the same canonical storage boundary");
    assert_eq!(
        mol.stereo_element(canonical).unwrap(),
        mol.stereo_element(permuted).unwrap()
    );
}

#[test]
fn double_bond_stereo_storage_canonicalizes_endpoints_and_references() {
    let mut mol = Molecule::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left_reference = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let left_alternative = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let right_reference = mol
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let right_alternative = mol
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");
    let double_bond = mol
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    for (endpoint, carrier) in [
        (left, left_reference),
        (left, left_alternative),
        (right, right_reference),
        (right, right_alternative),
    ] {
        mol.add_bond(endpoint, carrier, BondOrder::Single)
            .expect("substituent bond");
    }

    let add = |mol: &mut Molecule,
               left_endpoint,
               right_endpoint,
               left_carrier,
               right_carrier,
               orientation| {
        mol.add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: double_bond,
                left: left_endpoint,
                right: right_endpoint,
                left_carrier: StereoCarrier::Atom(left_carrier),
                right_carrier: StereoCarrier::Atom(right_carrier),
                orientation,
            },
        )))
        .expect("double-bond stereo element")
    };
    let canonical = add(
        &mut mol,
        left,
        right,
        left_reference,
        right_reference,
        Some(DoubleBondOrientation::Together),
    );
    let alternate_left = add(
        &mut mol,
        left,
        right,
        left_alternative,
        right_reference,
        Some(DoubleBondOrientation::Opposite),
    );
    let reversed_and_alternate = add(
        &mut mol,
        right,
        left,
        right_alternative,
        left_alternative,
        Some(DoubleBondOrientation::Together),
    );
    for equivalent in [alternate_left, reversed_and_alternate] {
        assert_eq!(
            mol.stereo_element(canonical).unwrap(),
            mol.stereo_element(equivalent).unwrap()
        );
    }

    let unknown = add(&mut mol, left, right, left_reference, right_reference, None);
    let unknown_alternatives = add(
        &mut mol,
        right,
        left,
        right_alternative,
        left_alternative,
        None,
    );
    assert_eq!(
        mol.stereo_element(unknown).unwrap(),
        mol.stereo_element(unknown_alternatives).unwrap()
    );
}

#[test]
fn axis_stereo_storage_canonicalizes_reference_carriers() {
    let mut mol = Molecule::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left_reference = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let left_alternative = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let right_reference = mol
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let right_alternative = mol
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");
    let axis = mol
        .add_bond(left, right, BondOrder::Single)
        .expect("axis bond");
    for (endpoint, carrier) in [
        (left, left_reference),
        (left, left_alternative),
        (right, right_reference),
        (right, right_alternative),
    ] {
        mol.add_bond(endpoint, carrier, BondOrder::Single)
            .expect("axis substituent bond");
    }

    let add = |mol: &mut Molecule, carriers, orientation| {
        mol.add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers,
            orientation,
        })))
        .expect("axis stereo element")
    };
    let canonical = add(
        &mut mol,
        vec![
            StereoCarrier::Atom(left_reference),
            StereoCarrier::Atom(right_reference),
        ],
        Some(AxisOrientation::Clockwise),
    );
    let reversed = add(
        &mut mol,
        vec![
            StereoCarrier::Atom(right_reference),
            StereoCarrier::Atom(left_reference),
        ],
        Some(AxisOrientation::Clockwise),
    );
    let alternate_left = add(
        &mut mol,
        vec![
            StereoCarrier::Atom(left_alternative),
            StereoCarrier::Atom(right_reference),
        ],
        Some(AxisOrientation::CounterClockwise),
    );
    let both_alternatives = add(
        &mut mol,
        vec![
            StereoCarrier::Atom(right_alternative),
            StereoCarrier::Atom(left_alternative),
        ],
        Some(AxisOrientation::Clockwise),
    );
    for equivalent in [reversed, alternate_left, both_alternatives] {
        assert_eq!(
            mol.stereo_element(canonical).unwrap(),
            mol.stereo_element(equivalent).unwrap()
        );
    }

    let unknown = add(
        &mut mol,
        vec![
            StereoCarrier::Atom(left_reference),
            StereoCarrier::Atom(right_reference),
        ],
        None,
    );
    let unknown_alternatives = add(
        &mut mol,
        vec![
            StereoCarrier::Atom(right_alternative),
            StereoCarrier::Atom(left_alternative),
        ],
        None,
    );
    assert_eq!(
        mol.stereo_element(unknown).unwrap(),
        mol.stereo_element(unknown_alternatives).unwrap()
    );
}

#[test]
fn stereo_element_group_membership_is_transactional_and_relation_owned() {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let a = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    let element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(a),
                    StereoCarrier::Atom(b),
                    StereoCarrier::Atom(c),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("stereo element");
    let group = mol
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: vec![element],
        })
        .expect("stereo group");
    let perception = PerceptionState::builder()
        .with_cip_descriptors(vec![(element, StereoDescriptor::R)])
        .expect("unique CIP assignment")
        .build();
    mol.install_perception_state(perception.clone())
        .expect("valid perception");

    let mut pre_grouped = mol.stereo_element(element).expect("element").clone();
    let StereoElementKind::Tetrahedral(stereo) = &mut pre_grouped.kind else {
        unreachable!("test element is tetrahedral");
    };
    stereo.center = AtomId::new(999);
    let slots_before = mol.stereo_elements.clone();
    assert!(matches!(
        mol.add_stereo_element(pre_grouped),
        Err(MoleculeError::InvalidStereoReference(
            "stereo element group membership must be established through add_stereo_group"
        ))
    ));
    assert_eq!(mol.stereo_elements.len(), slots_before.len());
    assert_eq!(mol.stereo_elements, slots_before);
    assert_eq!(mol.perception(), &perception);

    let removed = mol
        .remove_stereo_element(element)
        .expect("grouped element removal");
    assert_eq!(removed.group, None);
    assert!(mol.stereo_group(group).is_err());

    let readded = mol
        .add_stereo_element(removed)
        .expect("detached element can be re-added");
    assert_eq!(readded, StereoElementId::new(1));
    let regrouped = mol
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Relative,
            members: vec![readded],
        })
        .expect("re-added element can be grouped");
    assert_eq!(
        mol.stereo_element(readded).expect("re-added element").group,
        Some(regrouped)
    );
}

#[test]
fn topology_deletions_prune_referencing_stereo_state() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let ab = mol.add_bond(a, b, BondOrder::Double).expect("double bond");
    let ac = mol.add_bond(a, c, BondOrder::Single).expect("single bond");
    let bc = mol.add_bond(b, c, BondOrder::Single).expect("single bond");

    let element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: ab,
                left: a,
                right: b,
                left_carrier: StereoCarrier::Atom(c),
                right_carrier: StereoCarrier::Atom(c),
                orientation: Some(DoubleBondOrientation::Opposite),
            },
        )))
        .expect("double-bond element");
    mol.add_stereo_group(StereoGroup {
        kind: StereoGroupKind::Relative,
        members: vec![element],
    })
    .expect("group");
    mol.delete_bond(ab).expect("delete double bond");
    assert!(mol.stereo_element(element).is_err());
    assert!(mol
        .stereo_groups()
        .all(|(_, group)| group.members.is_empty()));

    mol.delete_bond(ac).expect("delete bond");

    let atom_element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center: c,
                carriers: vec![StereoCarrier::Atom(b), StereoCarrier::ImplicitHydrogen],
                orientation: Some(TetrahedralOrientation::CounterClockwise),
            },
        )))
        .expect("atom element");
    mol.delete_atom(c).expect("delete atom");
    assert!(mol.stereo_element(atom_element).is_err());
    assert!(mol.bond(bc).is_err());
}

#[test]
fn prop_value_equality_covers_all_initial_variants() {
    assert_eq!(
        PropValue::String("value".to_owned()),
        PropValue::String("value".to_owned())
    );
    assert_eq!(PropValue::Int(42), PropValue::Int(42));
    assert_eq!(PropValue::Float(2.5), PropValue::Float(2.5));
    assert_eq!(PropValue::Bool(true), PropValue::Bool(true));
}

#[test]
fn mutable_payload_access_invalidates_fresh_perception() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");

    mark_all_fresh(&mut mol);
    mol.atom_mut(a).expect("atom exists").formal_charge = 1;
    assert_all_stale(&mol);

    mark_all_fresh(&mut mol);
    mol.bond_mut(bond).expect("bond exists").order = BondOrder::Double;
    assert_all_stale(&mol);
}

#[test]
fn perception_owned_chemistry_edits_invalidate_dependent_state() {
    let mut methane = Molecule::new();
    methane
        .add_atom(carbon())
        .expect("atom identifier capacity");
    mark_all_fresh(&mut methane);

    let report = valence_api::perceive_valence(&mut methane, ValenceModel::RdkitLike);

    assert!(report.is_ok());
    assert!(methane.perception().has_valence());
    assert!(methane.perception().has_rings());
    assert!(!methane.perception().has_aromaticity());
    assert!(!methane.perception().has_stereo());

    let (mut benzene, _, _) = ring_molecule(
        &["C", "C", "C", "C", "C", "C"],
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    mark_all_fresh(&mut benzene);

    aromaticity_api::perceive_aromaticity(&mut benzene, AromaticityModel::RdkitLike)
        .expect("benzene should be supported");

    assert!(benzene.perception().has_valence());
    assert!(benzene.perception().has_rings());
    assert!(benzene.perception().has_aromaticity());
    assert!(!benzene.perception().has_stereo());
}
