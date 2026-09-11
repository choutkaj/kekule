use super::*;
use crate::properties::{PropertyKey, PropertyValue};

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
}

#[test]
fn atom_payload_fields_can_be_set_and_read() {
    let mut atom = carbon();
    atom.isotope = Some(13);
    atom.formal_charge = -1;
    atom.radical = Some(AtomRadical::Doublet);
    atom.hydrogens = HydrogenDeclaration::Fixed(3);
    atom.atom_map = Some(7);

    assert_eq!(atom.isotope, Some(13));
    assert_eq!(atom.formal_charge, -1);
    assert_eq!(atom.radical, Some(AtomRadical::Doublet));
    assert_eq!(atom.hydrogens, HydrogenDeclaration::Fixed(3));
    assert_eq!(atom.atom_map, Some(7));
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
    assert_eq!(double.order, BondOrder::Double);
}

#[test]
fn bond_payload_fields_can_be_set_and_read() {
    let bond = Bond::new(AtomId::new(1), AtomId::new(2), BondOrder::Dative);
    assert_eq!(bond.order, BondOrder::Dative);
}

#[test]
fn stereo_elements_and_groups_live_on_molecule() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let a = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(center, a, BondOrder::Single).expect("bond");
    mol.add_bond(center, b, BondOrder::Single).expect("bond");
    mol.add_bond(center, c, BondOrder::Single).expect("bond");
    mark_all_fresh(mol.working_mut());

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
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let a = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    for carrier in [a, b, c] {
        mol.add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
    }
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
    let mut mol = crate::core::MoleculeEditor::new();
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
    let canonical = mol.remove_stereo_element(canonical).unwrap();
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
    assert_eq!(&canonical, mol.stereo_element(permuted).unwrap());

    let mut unknown_mol = mol.clone();
    unknown_mol.remove_stereo_element(permuted).unwrap();
    let unknown = unknown_mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: canonical_carriers.clone(),
                orientation: None,
            },
        )))
        .expect("canonical unknown tetrahedral element");
    let unknown = unknown_mol.remove_stereo_element(unknown).unwrap();
    let unknown_permuted = unknown_mol
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
        &unknown,
        unknown_mol.stereo_element(unknown_permuted).unwrap()
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
    assert_eq!(&canonical, mol.stereo_element(permuted).unwrap());
}

#[test]
fn double_bond_stereo_storage_canonicalizes_endpoints_and_references() {
    let mut mol = crate::core::MoleculeEditor::new();
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
        let mut variant = mol.clone();
        let id = variant
            .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
                DoubleBondStereo {
                    bond: double_bond,
                    left: left_endpoint,
                    right: right_endpoint,
                    left_carrier: StereoCarrier::Atom(left_carrier),
                    right_carrier: StereoCarrier::Atom(right_carrier),
                    orientation,
                },
            )))
            .expect("double-bond stereo element");
        variant.stereo_element(id).unwrap().clone()
    };
    let canonical = add(
        mol.working_mut(),
        left,
        right,
        left_reference,
        right_reference,
        Some(DoubleBondOrientation::Together),
    );
    let alternate_left = add(
        mol.working_mut(),
        left,
        right,
        left_alternative,
        right_reference,
        Some(DoubleBondOrientation::Opposite),
    );
    let reversed_and_alternate = add(
        mol.working_mut(),
        right,
        left,
        right_alternative,
        left_alternative,
        Some(DoubleBondOrientation::Together),
    );
    for equivalent in [alternate_left, reversed_and_alternate] {
        assert_eq!(canonical, equivalent);
    }

    let unknown = add(
        mol.working_mut(),
        left,
        right,
        left_reference,
        right_reference,
        None,
    );
    let unknown_alternatives = add(
        mol.working_mut(),
        right,
        left,
        right_alternative,
        left_alternative,
        None,
    );
    assert_eq!(unknown, unknown_alternatives);
}

#[test]
fn axis_stereo_storage_canonicalizes_reference_carriers() {
    let mut mol = crate::core::MoleculeEditor::new();
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
        let mut variant = mol.clone();
        let id = variant
            .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
                axis,
                carriers,
                orientation,
            })))
            .expect("axis stereo element");
        variant.stereo_element(id).unwrap().clone()
    };
    let canonical = add(
        mol.working_mut(),
        vec![
            StereoCarrier::Atom(left_reference),
            StereoCarrier::Atom(right_reference),
        ],
        Some(AxisOrientation::Clockwise),
    );
    let reversed = add(
        mol.working_mut(),
        vec![
            StereoCarrier::Atom(right_reference),
            StereoCarrier::Atom(left_reference),
        ],
        Some(AxisOrientation::Clockwise),
    );
    let alternate_left = add(
        mol.working_mut(),
        vec![
            StereoCarrier::Atom(left_alternative),
            StereoCarrier::Atom(right_reference),
        ],
        Some(AxisOrientation::CounterClockwise),
    );
    let both_alternatives = add(
        mol.working_mut(),
        vec![
            StereoCarrier::Atom(right_alternative),
            StereoCarrier::Atom(left_alternative),
        ],
        Some(AxisOrientation::Clockwise),
    );
    for equivalent in [reversed, alternate_left, both_alternatives] {
        assert_eq!(canonical, equivalent);
    }

    let unknown = add(
        mol.working_mut(),
        vec![
            StereoCarrier::Atom(left_reference),
            StereoCarrier::Atom(right_reference),
        ],
        None,
    );
    let unknown_alternatives = add(
        mol.working_mut(),
        vec![
            StereoCarrier::Atom(right_alternative),
            StereoCarrier::Atom(left_alternative),
        ],
        None,
    );
    assert_eq!(unknown, unknown_alternatives);
}

#[test]
fn stereo_element_group_membership_is_transactional_and_relation_owned() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let a = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    for carrier in [a, b, c] {
        mol.add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
    }
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
    let perception = Perception::builder()
        .with_cip_descriptors(vec![(element, StereoDescriptor::R)])
        .expect("unique CIP assignment")
        .build();
    mol.working_mut()
        .install_perception(perception.clone())
        .expect("valid perception");

    let mut pre_grouped = mol.stereo_element(element).expect("element").clone();
    let StereoElementKind::Tetrahedral(stereo) = &mut pre_grouped.kind else {
        unreachable!("test element is tetrahedral");
    };
    stereo.center = AtomId::new(999);
    let slots_before = mol.working().graph.stereo_elements.clone();
    assert!(matches!(
        mol.add_stereo_element(pre_grouped),
        Err(MoleculeError::InvalidStereoReference(
            "stereo element group membership must be established through add_stereo_group"
        ))
    ));
    assert_eq!(
        mol.working().graph.stereo_elements.len(),
        slots_before.len()
    );
    assert_eq!(mol.working().graph.stereo_elements, slots_before);
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
    let mut mol = crate::core::MoleculeEditor::new();
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

    let fluorine = mol.add_atom(element_atom("F")).unwrap();
    let chlorine = mol.add_atom(element_atom("Cl")).unwrap();
    mol.add_bond(c, fluorine, BondOrder::Single).unwrap();
    mol.add_bond(c, chlorine, BondOrder::Single).unwrap();
    let atom_element = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center: c,
                carriers: vec![
                    StereoCarrier::Atom(b),
                    StereoCarrier::Atom(fluorine),
                    StereoCarrier::Atom(chlorine),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: Some(TetrahedralOrientation::CounterClockwise),
            },
        )))
        .expect("atom element");
    mol.delete_atom(c).expect("delete atom");
    assert!(mol.stereo_element(atom_element).is_err());
    assert!(mol.bond(bc).is_err());
}

#[test]
fn property_value_equality_covers_all_initial_variants() {
    assert_eq!(
        PropertyValue::String("value".to_owned()),
        PropertyValue::String("value".to_owned())
    );
    assert_eq!(PropertyValue::Int(42), PropertyValue::Int(42));
    assert_eq!(
        PropertyValue::Real {
            value: 2.5,
            unit: crate::units::DIMENSIONLESS
        },
        PropertyValue::Real {
            value: 2.5,
            unit: crate::units::DIMENSIONLESS
        }
    );
    assert_eq!(PropertyValue::Bool(true), PropertyValue::Bool(true));
}

#[test]
fn checked_stereo_rejects_invalid_carrier_shapes_transactionally() {
    let mut molecule = read_smiles("F[C@](Cl)(Br)I").unwrap();
    perceive(&mut molecule).unwrap();
    stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
    let (id, element) = molecule.stereo_elements().next().unwrap();
    let element = element.clone();
    let mut editor = molecule.edit();
    for duplicate in [false, true] {
        let mut invalid = element.clone();
        let StereoElementKind::Tetrahedral(stereo) = &mut invalid.kind else {
            unreachable!();
        };
        if duplicate {
            stereo.carriers[3] = stereo.carriers[0];
        } else {
            stereo.carriers.pop();
        }
        let before = editor.clone();
        assert!(matches!(
            editor.add_stereo_element(invalid.clone()),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        assert!(matches!(
            editor.replace_stereo_element(id, invalid.clone()),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        assert_eq!(editor, before);
        assert_eq!(editor.perception(), before.perception());

        let mut unchecked = before.clone();
        unchecked.working_mut().graph.stereo_elements[id.index()] = Some(invalid);
        assert!(matches!(
            unchecked.finish(),
            Err(MoleculePublicationError::InvalidStereo(_))
        ));
    }
}

#[test]
fn axis_stereo_rejects_duplicate_and_same_endpoint_references() {
    let (molecule, _, axis) = coordinate_axis_graph(true);
    for carriers in [
        vec![StereoCarrier::Atom(AtomId::new(2))],
        vec![StereoCarrier::Atom(AtomId::new(2)); 2],
        vec![
            StereoCarrier::Atom(AtomId::new(2)),
            StereoCarrier::Atom(AtomId::new(3)),
        ],
        vec![
            StereoCarrier::ImplicitHydrogen,
            StereoCarrier::Atom(AtomId::new(4)),
        ],
    ] {
        let mut editor = molecule.edit();
        let invalid = StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers,
            orientation: Some(AxisOrientation::Clockwise),
        }));
        let before = editor.clone();
        assert!(matches!(
            editor.add_stereo_element(invalid.clone()),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        assert_eq!(editor, before);
        editor
            .working_mut()
            .graph
            .stereo_elements
            .push(Some(invalid));
        assert!(stereo_api::validate_stereo(editor.working()).is_err());
        assert!(matches!(
            editor.finish(),
            Err(MoleculePublicationError::InvalidStereo(_))
        ));
    }
}

#[test]
fn stereo_focus_is_unique_and_rejection_preserves_groups_and_perception() {
    let molecule = read_smiles("F[C@H](Cl)[C@H](F)Cl").unwrap();
    let mut editor = molecule.into_editor();
    let ids = editor
        .stereo_elements()
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    let group = editor
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Relative,
            members: ids.clone(),
        })
        .unwrap();
    perceive(editor.working_mut()).unwrap();
    stereo_api::assign_cip_descriptors(editor.working_mut()).unwrap();
    let first = editor.stereo_element(ids[0]).unwrap().clone();
    for unknown in [false, true] {
        let mut duplicate = first.clone();
        duplicate.group = None;
        if unknown {
            let StereoElementKind::Tetrahedral(stereo) = &mut duplicate.kind else {
                unreachable!()
            };
            stereo.orientation = None;
        }
        let before = editor.clone();
        assert!(matches!(
            editor.add_stereo_element(duplicate.clone()),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        duplicate.group = Some(group);
        assert!(matches!(
            editor.replace_stereo_element(ids[1], duplicate.clone()),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        assert_eq!(editor, before);
        assert_eq!(editor.perception(), before.perception());
        assert_eq!(editor.stereo_group(group).unwrap().members, ids);

        let mut unchecked = editor.clone();
        unchecked.working_mut().graph.stereo_elements[ids[1].index()] = Some(duplicate);
        let issue = StereoValidationIssue::DuplicateStereoFocus {
            element: ids[1],
            previous: ids[0],
        };
        assert!(stereo_api::validate_stereo(unchecked.working())
            .unwrap_err()
            .issues
            .contains(&issue));
        assert_eq!(
            unchecked.finish().unwrap_err(),
            MoleculePublicationError::InvalidStereo(
                StereoPublicationError::DuplicateElementFocus {
                    element: ids[1],
                    previous: ids[0]
                }
            )
        );
    }
    // Replacing the existing focus is the supported way to update an assertion.
    let mut replacement = first;
    let StereoElementKind::Tetrahedral(stereo) = &mut replacement.kind else {
        unreachable!()
    };
    stereo.orientation = None;
    editor.replace_stereo_element(ids[0], replacement).unwrap();
    assert_eq!(editor.stereo_group(group).unwrap().members, ids);
    assert!(editor
        .stereo_element(ids[0])
        .unwrap()
        .is_explicitly_unknown());
    assert!(!editor.perception().has_stereo());
}

#[test]
fn bond_stereo_focus_rejects_duplicate_or_conflicting_axis_assertions() {
    let molecule = read_smiles("F/C=C/Cl").unwrap();
    let mut editor = molecule.into_editor();
    let element = editor.stereo_elements().next().unwrap().1.clone();
    assert!(matches!(
        editor.add_stereo_element(element.clone()),
        Err(MoleculeError::InvalidStereoReference(_))
    ));
    let StereoElementKind::DoubleBond(stereo) = element.kind else {
        unreachable!()
    };
    let axis = StereoElement::new(StereoElementKind::Axis(AxisStereo {
        axis: stereo.bond,
        carriers: vec![stereo.left_carrier, stereo.right_carrier],
        orientation: Some(AxisOrientation::Clockwise),
    }));
    assert!(matches!(
        editor.add_stereo_element(axis),
        Err(MoleculeError::InvalidStereoReference(_))
    ));
    assert_eq!(editor.stereo_elements().count(), 1);
}

#[test]
fn axis_reference_canonicalization_preserves_exclusive_endpoint_adjacency() {
    let mut editor = MoleculeEditor::new();
    let atoms = (0..5)
        .map(|_| editor.add_atom(carbon()).unwrap())
        .collect::<Vec<_>>();
    let axis = editor
        .add_bond(atoms[0], atoms[1], BondOrder::Single)
        .unwrap();
    for (left, right) in [(0, 2), (1, 2), (0, 3), (1, 4)] {
        editor
            .add_bond(atoms[left], atoms[right], BondOrder::Single)
            .unwrap();
    }
    let element = StereoElement::new(StereoElementKind::Axis(AxisStereo {
        axis,
        carriers: vec![StereoCarrier::Atom(atoms[3]), StereoCarrier::Atom(atoms[4])],
        orientation: Some(AxisOrientation::Clockwise),
    }));
    let id = editor.add_stereo_element(element.clone()).unwrap();
    assert_eq!(editor.stereo_element(id).unwrap(), &element);
    stereo_api::validate_stereo(editor.working()).unwrap();
    editor.finish().unwrap();
}

#[test]
fn tetrahedral_stereo_requires_complete_explicit_neighbor_coverage() {
    let molecule = read_smiles("F[C@](Cl)(Br)I").unwrap();
    let mut editor = molecule.into_editor();
    let (id, element) = editor.stereo_elements().next().unwrap();
    let element = element.clone();
    let StereoElementKind::Tetrahedral(stereo) = &element.kind else {
        unreachable!()
    };
    let center = stereo.center;
    let extra = editor.add_atom(carbon()).unwrap();
    editor.add_bond(center, extra, BondOrder::Single).unwrap();
    let before = editor.clone();
    assert!(matches!(
        editor.replace_stereo_element(id, element.clone()),
        Err(MoleculeError::InvalidStereoReference(_))
    ));
    assert_eq!(editor, before);
    assert!(stereo_api::validate_stereo(editor.working())
        .unwrap_err()
        .issues
        .contains(&StereoValidationIssue::UnrepresentedTetrahedralNeighbor {
            element: id,
            center,
            neighbor: extra
        }));
    assert!(matches!(
        editor.clone().finish(),
        Err(MoleculePublicationError::InvalidStereo(_))
    ));
    editor.remove_stereo_element(id).unwrap();
    assert!(matches!(
        editor.add_stereo_element(element),
        Err(MoleculeError::InvalidStereoReference(_))
    ));

    // The multiple bond is one spatial ligand; a sulfoxide's fourth carrier
    // is its lone pair, and all three explicit neighbors remain represented.
    let sulfoxide = read_smiles("C[S@](=O)CC").unwrap();
    stereo_api::validate_stereo(&sulfoxide).unwrap();
    sulfoxide.into_editor().finish().unwrap();
}

#[test]
fn stereo_bond_endpoints_reject_more_than_two_explicit_substituents() {
    for order in [BondOrder::Single, BondOrder::Double] {
        let mut editor = MoleculeEditor::new();
        let atoms = ["C", "C", "F", "Cl", "Br", "I"]
            .map(|symbol| editor.add_atom(element_atom(symbol)).unwrap());
        let focus = editor.add_bond(atoms[0], atoms[1], order).unwrap();
        for neighbor in &atoms[2..5] {
            editor
                .add_bond(atoms[0], *neighbor, BondOrder::Single)
                .unwrap();
        }
        editor
            .add_bond(atoms[1], atoms[5], BondOrder::Single)
            .unwrap();
        let kind = if order == BondOrder::Double {
            StereoElementKind::DoubleBond(DoubleBondStereo {
                bond: focus,
                left: atoms[0],
                right: atoms[1],
                left_carrier: StereoCarrier::Atom(atoms[2]),
                right_carrier: StereoCarrier::Atom(atoms[5]),
                orientation: Some(DoubleBondOrientation::Together),
            })
        } else {
            StereoElementKind::Axis(AxisStereo {
                axis: focus,
                carriers: vec![StereoCarrier::Atom(atoms[2]), StereoCarrier::Atom(atoms[5])],
                orientation: Some(AxisOrientation::Clockwise),
            })
        };
        let element = StereoElement::new(kind);
        let before = editor.clone();
        assert!(matches!(
            editor.add_stereo_element(element.clone()),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        assert_eq!(editor, before);
        editor
            .working_mut()
            .graph
            .stereo_elements
            .push(Some(element));
        let issues = stereo_api::validate_stereo(editor.working())
            .unwrap_err()
            .issues;
        assert!(issues.iter().any(|issue| matches!(issue,
            StereoValidationIssue::DoubleBondEndpointOvercoordinated { endpoint, substituent_count: 3, .. }
                | StereoValidationIssue::AxisEndpointOvercoordinated { endpoint, substituent_count: 3, .. }
                if *endpoint == atoms[0]
        )));
        assert!(matches!(
            editor.finish(),
            Err(MoleculePublicationError::InvalidStereo(_))
        ));
    }
}

#[test]
fn mutable_payload_access_invalidates_fresh_perception() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");

    mark_all_fresh(mol.working_mut());
    mol.atom_mut(a).expect("atom exists").formal_charge = 1;
    assert_all_stale(mol.working());

    mark_all_fresh(mol.working_mut());
    mol.bond_mut(bond)
        .expect("bond exists")
        .set_order(BondOrder::Double);
    assert_all_stale(mol.working());
}

#[test]
fn atom_map_only_mutation_invalidates_perception_and_owner_properties() {
    let mut molecule = crate::core::MoleculeEditor::new();
    let atom = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let owner_key = PropertyKey::new("calculation_label").unwrap();
    molecule
        .insert_property(
            owner_key.clone(),
            PropertyValue::String("before mutation".to_owned()),
        )
        .unwrap();
    mark_all_fresh(molecule.working_mut());

    molecule.atom_mut(atom).unwrap().atom_map = Some(7);

    assert_all_stale(molecule.working());
    assert_eq!(molecule.properties().get(&owner_key), None);
}

#[test]
fn perception_owned_chemistry_edits_invalidate_dependent_state() {
    let mut methane = crate::core::MoleculeEditor::new();
    methane
        .add_atom(carbon())
        .expect("atom identifier capacity");
    mark_all_fresh(methane.working_mut());

    let report = valence_api::perceive_valence(methane.working_mut(), ValenceModel::RdkitLike);

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
