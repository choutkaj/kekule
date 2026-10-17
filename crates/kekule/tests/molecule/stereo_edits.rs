use std::sync::Arc;

use kekule::topology::{InstanceBondId, Topology};
use kekule::{
    molecule::{Atom, AtomId, Bond, BondId, BondOrder, Element, Molecule, MoleculeError},
    stereo::{
        AxisStereo, StereoCarrier, StereoElement, StereoElementId, StereoElementKind, StereoGroup,
        StereoGroupKind,
    },
};

fn molecule(text: &str) -> Molecule {
    kekule::smiles::to_molecules(text).unwrap().pop().unwrap()
}

#[test]
fn changed_bond_orders_prune_stereo_through_all_molecular_mutators() {
    let source = molecule("F/C=C/F");
    let bond = source
        .bonds()
        .find(|(_, bond)| bond.order == BondOrder::Double)
        .unwrap()
        .0;
    let (a, b) = source.bond(bond).unwrap().endpoints();
    for operation in 0..3 {
        let mut editor = source.to_editor();
        match operation {
            0 => editor.set_bond_order(bond, BondOrder::Single).unwrap(),
            1 => editor.bond_mut(bond).unwrap().set_order(BondOrder::Single),
            _ => {
                editor
                    .replace_bond(bond, Bond::new(b, a, BondOrder::Single))
                    .unwrap();
            }
        }
        assert_eq!(
            editor.pruned_stereo_elements(),
            [StereoElementId::new(0)],
            "operation {operation}"
        );
        editor.validate().unwrap();
        let result = editor.try_finish().unwrap();
        assert_eq!(result.bond(bond).unwrap().order, BondOrder::Single);
        assert_eq!(result.stereo_elements().count(), 0);
        let text =
            kekule::smiles::write(&result, kekule::smiles::SmilesWriteOptions::isomeric()).unwrap();
        let reparsed = molecule(&text);
        assert_eq!(reparsed.stereo_elements().count(), 0);
        assert_eq!(
            kekule::smiles::write(&result, kekule::smiles::SmilesWriteOptions::canonical())
                .unwrap(),
            "FCCF"
        );
    }
    assert_eq!(source.stereo_elements().count(), 1);
}

#[test]
fn order_changes_update_stereo_groups_without_removing_unrelated_assertions() {
    let source = molecule("F/C=C/C=C/F");
    let bonds = source
        .bonds()
        .filter_map(|(id, bond)| (bond.order == BondOrder::Double).then_some(id))
        .collect::<Vec<_>>();
    let elements = source.stereo_element_ids().collect::<Vec<_>>();
    assert_eq!(elements.len(), 2);
    let mut editor = source.into_editor();
    let group = editor
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::And,
            members: elements,
        })
        .unwrap();
    let source = editor.finish().unwrap();
    let mut editor = source.to_editor();
    editor.set_bond_order(bonds[0], BondOrder::Double).unwrap();
    assert_eq!(editor.clone().finish().unwrap(), source);
    editor.set_bond_order(bonds[0], BondOrder::Single).unwrap();
    let retained = editor.stereo_element_ids().next().unwrap();
    assert_eq!(editor.stereo_element_ids().count(), 1);
    assert_eq!(editor.stereo_group(group).unwrap().members, vec![retained]);
    assert_eq!(editor.stereo_element(retained).unwrap().group, Some(group));
    let mut editor = editor.finish().unwrap().into_editor();
    editor.set_bond_order(bonds[1], BondOrder::Single).unwrap();
    let result = editor.finish().unwrap();
    assert_eq!(result.stereo_elements().count(), 0);
    assert_eq!(result.stereo_groups().count(), 0);
}

#[test]
fn checked_stereo_insertion_and_replacement_reject_single_bond_focus() {
    let source = molecule("F/C=C/C=C/F");
    let elements = source.stereo_elements().collect::<Vec<_>>();
    let invalid = elements[0].1.clone();
    let StereoElementKind::DoubleBond(stereo) = &invalid.kind else {
        panic!("expected double-bond stereo");
    };
    let focus = stereo.bond;
    let retained_element = elements[1].0;
    let mut editor = source.to_editor();
    editor.set_bond_order(focus, BondOrder::Single).unwrap();
    let before = format!("{editor:?}");
    assert!(editor.add_stereo_element(invalid.clone()).is_err());
    assert!(editor
        .replace_stereo_element(retained_element, invalid)
        .is_err());
    assert_eq!(format!("{editor:?}"), before);
    assert_eq!(editor.finish().unwrap().stereo_elements().count(), 1);
}

#[test]
fn topology_order_changes_publish_without_stale_double_bond_stereo() {
    for replacement in [false, true] {
        let source = Arc::new(Topology::from_molecule((molecule("F/C=C/F")).clone()).unwrap());
        let mut editor = source.to_editor();
        let (id, bond) = editor
            .bonds()
            .find(|(_, bond)| bond.order == BondOrder::Double)
            .unwrap();
        if replacement {
            editor
                .replace_bond(
                    id,
                    kekule::topology::edit::EditBond::new(bond.a(), bond.b(), BondOrder::Single),
                )
                .unwrap();
        } else {
            editor.set_bond_order(id, BondOrder::Single).unwrap();
        }
        let target = editor.finish().unwrap();
        let definition = target.instances().next().unwrap();
        assert_eq!(definition.molecule().stereo_elements().count(), 0);
        assert_eq!(
            source
                .instances()
                .next()
                .unwrap()
                .molecule()
                .stereo_elements()
                .count(),
            1
        );
    }
}

#[test]
fn bond_edits_prune_exactly_the_stereo_they_invalidate() {
    // Atoms: center 0, F 1, Cl 2, ring atoms 3-5. Bonds 3 and 4 are ring bonds
    // away from the center; bond 2 joins the center to the ring.
    let source = molecule("[C@@H](F)(Cl)C1CC1");
    let element = StereoElementId::new(0);
    // A new bond to the center (atom 0) adds a neighbor its carriers do not
    // list, so the assertion is pruned, whatever unrelated edits accompany
    // it. Bonds elsewhere, even to a carrier such as F, leave it intact.
    for (anchor, deleted, pruned) in [
        (0, None, vec![element]),
        (0, Some(4), vec![element]),
        (0, Some(3), vec![element]),
        (4, None, vec![]),
        (3, Some(4), vec![]),
        (1, Some(3), vec![]),
    ] {
        let mut editor = source.to_editor();
        let oxygen = editor
            .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
            .unwrap();
        editor
            .add_bond(AtomId::new(anchor), oxygen, BondOrder::Single)
            .unwrap();
        if let Some(bond) = deleted {
            editor.delete_bond(BondId::new(bond)).unwrap();
        }
        let row = format!("bond to {anchor}, delete {deleted:?}");
        assert_eq!(editor.pruned_stereo_elements(), pruned, "{row}");
        assert_eq!(editor.take_pruned_stereo_elements(), pruned, "{row}");
        assert_eq!(editor.pruned_stereo_elements(), [], "{row}");
        let result = editor.finish().unwrap();
        assert_eq!(result.stereo_elements().count(), 1 - pruned.len(), "{row}");
    }
}

#[test]
fn ring_bond_deletion_prunes_invalid_carriers_and_preserves_unrelated_stereo_groups() {
    let mut source = molecule("F[C@]1(Cl)CCNC1C[C@H](Br)I").into_editor();
    let elements = source.stereo_element_ids().collect::<Vec<_>>();
    assert_eq!(elements.len(), 2);
    let group = source
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::And,
            members: elements.clone(),
        })
        .unwrap();
    let source = source.finish().unwrap();
    let StereoElementKind::Tetrahedral(stereo) = &source.stereo_element(elements[0]).unwrap().kind
    else {
        panic!("expected tetrahedral center");
    };
    let carrier = stereo
        .carriers
        .iter()
        .find_map(|carrier| match carrier {
            StereoCarrier::Atom(atom) if source.neighbors(*atom).unwrap().count() > 1 => {
                Some(*atom)
            }
            _ => None,
        })
        .unwrap();
    let removed = source.bond_between(stereo.center, carrier).unwrap();
    for batch in [false, true] {
        let mut editor = source.to_editor();
        if batch {
            editor.delete_bonds([removed]).unwrap();
        } else {
            editor.delete_bond(removed).unwrap();
        }
        assert!(editor.stereo_element(elements[0]).is_none());
        assert_eq!(editor.pruned_stereo_elements(), [elements[0]]);
        assert_eq!(
            editor.stereo_element(elements[1]).unwrap(),
            source.stereo_element(elements[1]).unwrap()
        );
        assert_eq!(
            editor.stereo_group(group).unwrap().members,
            vec![elements[1]]
        );
        editor.validate().unwrap();
        let result = editor.try_finish().unwrap();
        // Enhanced groups have a separate export capability; test carrier export
        // on an ungrouped copy after asserting the preserved group above.
        let mut ungrouped = result.to_editor();
        ungrouped.remove_stereo_group(group).unwrap();
        // Deleting the bond leaves the fixed-H carbon under-coordinated without
        // asserting a radical. Explicitly cap that site for the export check.
        ungrouped.atom_mut(stereo.center).unwrap().hydrogens =
            kekule::molecule::ImplicitHydrogens::Fixed(1);
        let text = kekule::smiles::write(
            &ungrouped.finish().unwrap(),
            kekule::smiles::SmilesWriteOptions::isomeric(),
        )
        .unwrap();
        assert_eq!(molecule(&text).stereo_elements().count(), 1);
    }
    let topology = Arc::new(Topology::from_molecule(source.clone()).unwrap());
    let instance = topology.instances().next().unwrap().id();
    let mut editor = topology.to_editor();
    let bond = editor
        .bond_handle(InstanceBondId::new(instance, removed))
        .unwrap();
    editor.delete_bond(bond).unwrap();
    let result = editor.finish().unwrap();
    let result = result.instances().next().unwrap();
    assert_eq!(result.molecule().stereo_elements().count(), 1);
    let mut ungrouped = result.molecule().to_editor();
    for group in result.molecule().stereo_groups().map(|(id, _)| id) {
        ungrouped.remove_stereo_group(group).unwrap();
    }
    ungrouped.atom_mut(stereo.center).unwrap().hydrogens =
        kekule::molecule::ImplicitHydrogens::Fixed(1);
    kekule::smiles::write(
        &ungrouped.finish().unwrap(),
        kekule::smiles::SmilesWriteOptions::isomeric(),
    )
    .unwrap();
    assert_eq!(source.stereo_elements().count(), 2);
}

#[test]
fn checked_stereo_boundaries_reject_nonadjacent_carriers_transactionally() {
    let mut axis = molecule("FCCF").into_editor();
    let atoms = axis.atom_ids().collect::<Vec<_>>();
    axis.add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
        axis: axis.bond_between(atoms[1], atoms[2]).unwrap(),
        carriers: vec![StereoCarrier::Atom(atoms[0]), StereoCarrier::Atom(atoms[3])],
        orientation: None,
    })))
    .unwrap();
    for source in [
        molecule("F[C@]1(Cl)CCNC1"),
        molecule("F/C=C/F"),
        axis.finish().unwrap(),
    ] {
        let (id, original) = source.stereo_elements().next().unwrap();
        let mut invalid = original.clone();
        match &mut invalid.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                stereo.carriers[0] = StereoCarrier::Atom(stereo.center)
            }
            StereoElementKind::DoubleBond(stereo) => {
                stereo.left_carrier = StereoCarrier::Atom(stereo.right)
            }
            StereoElementKind::Axis(stereo) => {
                stereo.carriers[0] = StereoCarrier::Atom(source.bond(stereo.axis).unwrap().a())
            }
            other => unreachable!("no fixture asserts {other:?}"),
        }
        let mut editor = source.to_editor();
        let before = format!("{editor:?}");
        assert!(matches!(
            editor.add_stereo_element(invalid.clone()),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        assert!(matches!(
            editor.replace_stereo_element(id, invalid),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        assert_eq!(format!("{editor:?}"), before);
        assert_eq!(editor.finish().unwrap(), source);
    }
}
