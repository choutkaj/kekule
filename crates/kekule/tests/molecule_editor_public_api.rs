use kekule::core::*;
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::units::{ANGSTROM, KELVIN, NANOMETER};

fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}
fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}
fn molecule(smiles: &str) -> Molecule {
    kekule::smiles::to_molecules(smiles).unwrap().pop().unwrap()
}
fn snapshot(editor: &MoleculeEditor) -> String {
    format!("{editor:#?}")
}

#[test]
fn public_editor_builds_inspects_and_edits_without_a_molecule_view() {
    let mut editor = MoleculeEditor::new();
    assert!(editor.is_empty());
    assert!(!editor.is_connected());
    let a = editor.add_atom(atom("C")).unwrap();
    let b = editor.add_atom(atom("O")).unwrap();
    assert_eq!(editor.connected_components(), vec![vec![a], vec![b]]);
    let bond = editor.add_bond(a, b, BondOrder::Single).unwrap();
    assert!(editor.is_connected());
    assert_eq!(editor.atom_ids().collect::<Vec<_>>(), vec![a, b]);
    assert_eq!(editor.bond_ids().collect::<Vec<_>>(), vec![bond]);
    assert_eq!(editor.atoms().count(), 2);
    assert_eq!(editor.bonds().count(), 1);
    assert_eq!(editor.neighbors(a).unwrap().collect::<Vec<_>>(), vec![b]);
    assert_eq!(editor.incident_bonds(b).unwrap().next().unwrap().0, bond);
    assert_eq!(editor.bond_between(a, b).unwrap(), Some(bond));
    editor.atom_mut(a).unwrap().formal_charge = 1;
    assert_eq!(editor.formal_charge(), 1);
    assert_eq!(editor.replace_atom(a, atom("N")).unwrap().formal_charge, 1);
    editor.bond_mut(bond).unwrap().set_order(BondOrder::Double);
    assert_eq!(editor.bond(bond).unwrap().order, BondOrder::Double);
    editor
        .insert_property(key("label"), PropertyValue::String("draft".into()))
        .unwrap();
    editor.validate().unwrap();
    let published = editor.finish().unwrap();
    let mut edited = published.edit();
    edited.atom_mut(a).unwrap().formal_charge = 1;
    assert_eq!(published.atom(a).unwrap().formal_charge, 0);
    assert!(edited.properties().get(&key("label")).is_none());
    let atom_pointer = published.atom(a).unwrap() as *const Atom;
    let moved = published.to_editor();
    assert_eq!(moved.atom(a).unwrap() as *const Atom, atom_pointer);
}

#[test]
fn rewiring_retains_ids_properties_and_checks_every_endpoint_before_mutation() {
    let mut editor = molecule("CCCC").to_editor();
    let ids = editor.atom_ids().collect::<Vec<_>>();
    let bonds = editor.bond_ids().collect::<Vec<_>>();
    editor
        .set_bond_property(bonds[0], key("score"), Some(PropertyValue::Int(3)))
        .unwrap();
    editor
        .insert_property(key("label"), PropertyValue::Int(8))
        .unwrap();
    let before = snapshot(&editor);
    for replacement in [
        Bond::new(ids[0], ids[0], BondOrder::Single),
        Bond::new(ids[1], ids[2], BondOrder::Single),
        Bond::new(ids[0], AtomId::new(999), BondOrder::Single),
    ] {
        assert!(editor.replace_bond(bonds[0], replacement).is_err());
        assert_eq!(snapshot(&editor), before);
    }
    editor.set_bond_endpoints(bonds[0], ids[0], ids[3]).unwrap();
    assert_eq!(editor.bond_between(ids[0], ids[1]).unwrap(), None);
    assert_eq!(editor.bond_between(ids[0], ids[3]).unwrap(), Some(bonds[0]));
    assert_eq!(
        editor.neighbors(ids[0]).unwrap().collect::<Vec<_>>(),
        vec![ids[3]]
    );
    assert_eq!(
        editor.neighbors(ids[1]).unwrap().collect::<Vec<_>>(),
        vec![ids[2]]
    );
    assert_eq!(
        editor.bond_property(bonds[0], &key("score")).unwrap(),
        Some(PropertyValue::Int(3))
    );
    assert!(editor.properties().owner_is_empty());
    editor.validate().unwrap();
    editor.finish().unwrap();
}

#[test]
fn batch_deletion_is_atomic_and_retains_surviving_identity() {
    let mut editor = molecule("CCCC").to_editor();
    let ids = editor.atom_ids().collect::<Vec<_>>();
    let bonds = editor.bond_ids().collect::<Vec<_>>();
    let before = snapshot(&editor);
    assert!(editor.delete_atoms([ids[0], AtomId::new(999)]).is_err());
    assert!(editor.delete_bonds([bonds[0], BondId::new(999)]).is_err());
    assert!(editor.retain_atoms([ids[0], AtomId::new(999)]).is_err());
    assert_eq!(snapshot(&editor), before);
    assert_eq!(editor.delete_bonds([bonds[0], bonds[0]]).unwrap().len(), 1);
    assert_eq!(editor.delete_atoms([ids[0], ids[0]]).unwrap().len(), 1);
    assert_eq!(editor.retain_atoms([ids[2], ids[3]]).unwrap().len(), 1);
    assert_eq!(editor.atom_ids().collect::<Vec<_>>(), ids[2..]);
    assert_eq!(editor.bond_ids().collect::<Vec<_>>(), vec![bonds[2]]);
    editor.validate().unwrap();
    editor.clear();
    assert!(editor.is_empty());
    assert_eq!(editor.bond_count(), 0);
    assert_eq!(editor.add_atom(atom("O")).unwrap(), AtomId::new(0));
}

#[test]
fn property_columns_use_live_order_and_batches_preserve_state_on_error() {
    let mut editor = molecule("CCC").to_editor();
    let ids = editor.atom_ids().collect::<Vec<_>>();
    let bond = editor.bond_ids().last().unwrap();
    editor.delete_atom(ids[0]).unwrap();
    editor
        .set_atom_property_column(
            key("length"),
            PropertyColumn::Real {
                unit: ANGSTROM,
                values: vec![Some(1.0), Some(2.0)],
            },
        )
        .unwrap();
    assert_eq!(
        editor.atom_properties().value(&key("length"), 0).unwrap(),
        None
    );
    assert_eq!(
        editor.atom_property(ids[1], &key("length")).unwrap(),
        Some(PropertyValue::real(1.0, ANGSTROM).unwrap())
    );
    editor
        .set_bond_property_column(key("score"), PropertyColumn::Int(vec![Some(4)]))
        .unwrap();
    assert_eq!(
        editor.bond_property(bond, &key("score")).unwrap(),
        Some(PropertyValue::Int(4))
    );
    let before = snapshot(&editor);
    assert!(editor
        .set_atom_property_column(key("length"), PropertyColumn::Int(vec![Some(1)]))
        .is_err());
    assert!(editor
        .set_atom_properties(
            key("length"),
            [
                (ids[1], Some(PropertyValue::real(2.0, ANGSTROM).unwrap())),
                (ids[2], Some(PropertyValue::real(3.0, KELVIN).unwrap()))
            ]
        )
        .is_err());
    assert!(editor
        .set_atom_properties(
            key("x"),
            [
                (ids[1], Some(PropertyValue::Int(2))),
                (ids[0], Some(PropertyValue::Int(3)))
            ]
        )
        .is_err());
    assert!(editor
        .set_bond_properties(
            key("x"),
            [
                (bond, Some(PropertyValue::Int(2))),
                (BondId::new(999), None)
            ]
        )
        .is_err());
    assert_eq!(snapshot(&editor), before);
    editor
        .set_atom_properties(
            key("length"),
            [(ids[1], Some(PropertyValue::real(0.5, NANOMETER).unwrap()))],
        )
        .unwrap();
    assert_eq!(
        editor.atom_property(ids[1], &key("length")).unwrap(),
        Some(PropertyValue::real(5.0, ANGSTROM).unwrap())
    );
    editor
        .set_bond_properties(key("score"), [(bond, Some(PropertyValue::Int(7)))])
        .unwrap();
    assert_eq!(
        editor
            .remove_bond_property_column(&key("score"))
            .unwrap()
            .value(bond.index())
            .unwrap(),
        Some(PropertyValue::Int(7))
    );
    assert!(editor.remove_atom_property_column(&key("length")).is_some());
    editor.finish().unwrap();
}

fn grouped_fragment() -> Molecule {
    let mut editor = molecule("F[C@](Cl)(Br)I").to_editor();
    let stereo = editor.stereo_element_ids().next().unwrap();
    editor
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::And,
            members: vec![stereo],
        })
        .unwrap();
    let atom = editor.atom_ids().next().unwrap();
    let bond = editor.bond_ids().next().unwrap();
    editor
        .set_atom_property(atom, key("tag"), Some(PropertyValue::Int(5)))
        .unwrap();
    editor
        .set_bond_property(bond, key("tag"), Some(PropertyValue::Int(6)))
        .unwrap();
    editor
        .insert_property(key("source"), PropertyValue::Int(9))
        .unwrap();
    editor.finish().unwrap()
}

#[test]
fn append_preserves_fragment_stereo_groups_and_entity_properties() {
    let source = grouped_fragment();
    let mut editor = molecule("C").to_editor();
    let existing = editor.atom_ids().next().unwrap();
    let map = editor.append_molecule(&source).unwrap();
    assert_eq!(editor.connected_components().len(), 2);
    assert_eq!(map.atoms().len(), source.atom_count());
    assert_eq!(map.bonds().len(), source.bond_count());
    let old_atom = source.atom_ids().next().unwrap();
    assert_eq!(
        editor
            .atom_property(map.atoms()[&old_atom], &key("tag"))
            .unwrap(),
        Some(PropertyValue::Int(5))
    );
    let old_bond = source.bond_ids().next().unwrap();
    assert_eq!(
        editor
            .bond_property(map.bonds()[&old_bond], &key("tag"))
            .unwrap(),
        Some(PropertyValue::Int(6))
    );
    let (old_id, old_element) = source.stereo_elements().next().unwrap();
    let new_element = editor
        .stereo_element(map.stereo_elements()[&old_id])
        .unwrap();
    let old_group = old_element.group.unwrap();
    let new_group = map.stereo_groups()[&old_group];
    assert_eq!(new_element.group, Some(new_group));
    assert_eq!(
        editor.stereo_group(new_group).unwrap().members,
        vec![map.stereo_elements()[&old_id]]
    );
    let StereoElementKind::Tetrahedral(old) = &old_element.kind else {
        panic!("tetrahedral fixture")
    };
    let StereoElementKind::Tetrahedral(new) = &new_element.kind else {
        panic!("tetrahedral copy")
    };
    assert_eq!(new.center, map.atoms()[&old.center]);
    assert_eq!(new.orientation, old.orientation);
    assert_eq!(
        new.carriers,
        old.carriers
            .iter()
            .map(|c| match c {
                StereoCarrier::Atom(id) => StereoCarrier::Atom(map.atoms()[id]),
                other => *other,
            })
            .collect::<Vec<_>>()
    );
    assert!(editor.properties().get(&key("source")).is_none());
    editor
        .add_bond(existing, map.atoms()[&old_atom], BondOrder::Single)
        .unwrap();
    editor.finish().unwrap();
}

#[test]
fn append_conflicts_roll_back_graph_properties_and_id_allocation() {
    let mut editor = molecule("C").to_editor();
    let id = editor.atom_ids().next().unwrap();
    editor
        .set_atom_property(id, key("tag"), Some(PropertyValue::String("text".into())))
        .unwrap();
    editor
        .insert_property(key("label"), PropertyValue::Int(4))
        .unwrap();
    let before = snapshot(&editor);
    assert!(editor.append_molecule(&grouped_fragment()).is_err());
    assert_eq!(snapshot(&editor), before);
    assert_eq!(editor.add_atom(atom("C")).unwrap(), AtomId::new(1));
}

#[test]
fn append_remaps_double_bond_and_axis_stereo_across_sparse_ids() {
    let mut source = molecule("F/C=C/F").to_editor();
    let atoms = source.atom_ids().collect::<Vec<_>>();
    let focus = source.bond_between(atoms[1], atoms[2]).unwrap().unwrap();
    let double = source.stereo_element_ids().next().unwrap();
    let axis = source
        .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis: focus,
            carriers: vec![StereoCarrier::Atom(atoms[0]), StereoCarrier::Atom(atoms[3])],
            orientation: Some(AxisOrientation::Clockwise),
        })))
        .unwrap();
    let removed = source.add_atom(atom("H")).unwrap();
    source.delete_atom(removed).unwrap();
    let source = source.finish().unwrap();
    let mut target = molecule("CC").to_editor();
    target.delete_atom(AtomId::new(0)).unwrap();
    let map = target.append_molecule(&source).unwrap();
    assert!(!map.atoms().contains_key(&removed));
    assert_eq!(map.atoms().len(), 4);
    for id in [double, axis] {
        let old = &source.stereo_element(id).unwrap().kind;
        let new = &target
            .stereo_element(map.stereo_elements()[&id])
            .unwrap()
            .kind;
        match (old, new) {
            (StereoElementKind::DoubleBond(old), StereoElementKind::DoubleBond(new)) => {
                assert_eq!(new.bond, map.bonds()[&old.bond]);
                assert_eq!(new.left, map.atoms()[&old.left]);
                assert_eq!(new.right, map.atoms()[&old.right]);
                assert_eq!(
                    new.left_carrier,
                    StereoCarrier::Atom(map.atoms()[&atoms[0]])
                );
                assert_eq!(
                    new.right_carrier,
                    StereoCarrier::Atom(map.atoms()[&atoms[3]])
                );
                assert_eq!(new.orientation, old.orientation);
            }
            (StereoElementKind::Axis(old), StereoElementKind::Axis(new)) => {
                assert_eq!(new.axis, map.bonds()[&old.axis]);
                assert_eq!(
                    new.carriers,
                    vec![
                        StereoCarrier::Atom(map.atoms()[&atoms[0]]),
                        StereoCarrier::Atom(map.atoms()[&atoms[3]])
                    ]
                );
                assert_eq!(new.orientation, old.orientation);
            }
            _ => panic!("stereo kind changed during append"),
        }
    }
    target
        .add_bond(AtomId::new(1), map.atoms()[&atoms[0]], BondOrder::Single)
        .unwrap();
    target.finish().unwrap();
}

#[test]
fn property_and_no_op_edits_preserve_perception_but_chemistry_changes_clear_it() {
    let mut source = molecule("CC");
    source.perceive().unwrap();
    let cached = source.perception().clone();
    assert_ne!(cached, Perception::default());
    let mut editor = source.to_editor();
    let id = editor.atom_ids().next().unwrap();
    let bond = editor.bond_ids().next().unwrap();
    editor
        .insert_property(key("label"), PropertyValue::Int(7))
        .unwrap();
    editor
        .set_atom_property(id, key("tag"), Some(PropertyValue::Int(2)))
        .unwrap();
    editor
        .replace_atom(id, editor.atom(id).unwrap().clone())
        .unwrap();
    editor.set_bond_order(bond, BondOrder::Single).unwrap();
    assert_eq!(editor.perception(), &cached);
    assert_eq!(
        editor.properties().get(&key("label")),
        Some(&PropertyValue::Int(7))
    );
    editor.set_bond_order(bond, BondOrder::Double).unwrap();
    assert_eq!(editor.perception(), &Perception::default());
    assert!(editor.properties().owner_is_empty());
    assert_eq!(
        editor.atom_property(id, &key("tag")).unwrap(),
        Some(PropertyValue::Int(2))
    );
}

#[test]
fn stereo_group_replacement_preserves_identity_and_rewiring_prunes_affected_stereo() {
    let mut editor = grouped_fragment().to_editor();
    let (id, group) = editor
        .stereo_groups()
        .next()
        .map(|(id, g)| (id, g.clone()))
        .unwrap();
    let before = snapshot(&editor);
    assert!(editor
        .replace_stereo_group(
            id,
            StereoGroup {
                kind: StereoGroupKind::Or,
                members: vec![]
            }
        )
        .is_err());
    assert_eq!(snapshot(&editor), before);
    editor
        .replace_stereo_group(
            id,
            StereoGroup {
                kind: StereoGroupKind::Or,
                members: group.members.clone(),
            },
        )
        .unwrap();
    assert_eq!(editor.stereo_group(id).unwrap().kind, StereoGroupKind::Or);
    assert_eq!(
        editor.stereo_element(group.members[0]).unwrap().group,
        Some(id)
    );
    let mut added = editor.stereo_element(group.members[0]).unwrap().clone();
    added.group = None;
    let added = editor.add_stereo_element(added).unwrap();
    let before = snapshot(&editor);
    assert!(editor
        .replace_stereo_group(
            id,
            StereoGroup {
                kind: StereoGroupKind::And,
                members: vec![added, added]
            }
        )
        .is_err());
    assert_eq!(snapshot(&editor), before);
    editor
        .replace_stereo_group(
            id,
            StereoGroup {
                kind: StereoGroupKind::And,
                members: vec![added],
            },
        )
        .unwrap();
    assert_eq!(editor.stereo_element(group.members[0]).unwrap().group, None);
    assert_eq!(editor.stereo_element(added).unwrap().group, Some(id));
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    let bond = editor.bond_ids().next().unwrap();
    let other = editor.add_atom(atom("C")).unwrap();
    editor.set_bond_endpoints(bond, atoms[0], other).unwrap();
    assert_eq!(editor.stereo_elements().count(), 0);
    assert_eq!(editor.stereo_groups().count(), 0);
}

#[test]
fn failed_finish_returns_exact_draft_for_repair() {
    let mut editor = MoleculeEditor::new();
    let a = editor.add_atom(atom("C")).unwrap();
    let b = editor.add_atom(atom("O")).unwrap();
    editor
        .insert_property(key("note"), PropertyValue::Int(3))
        .unwrap();
    let before = snapshot(&editor);
    assert!(editor.validate().is_err());
    assert_eq!(snapshot(&editor), before);
    let error = editor.try_finish().unwrap_err();
    assert!(matches!(
        error.error(),
        MoleculePublicationError::DisconnectedGraph(_)
    ));
    assert_eq!(snapshot(error.editor()), before);
    let mut editor = error.to_editor();
    editor.add_bond(a, b, BondOrder::Single).unwrap();
    assert_eq!(editor.try_finish().unwrap().atom_count(), 2);
}
