use std::sync::Arc;

use kekule::core::{
    Atom, Bond, BondOrder, Element, Molecule, MoleculeEditor, StereoElementKind, StereoGroup,
    StereoGroupKind,
};
use kekule::topology::{
    transform, AtomSelection, AtomSiteMetadata, EditAtomId, Hierarchy, InstanceAtomId,
    MoleculeClass, ResidueClass, Topology, TopologyBuilder, TopologyEditor,
};

fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}

fn molecule(text: &str) -> Molecule {
    kekule::smiles::to_molecules(text).unwrap().pop().unwrap()
}

fn single_atom(symbol: &str) -> Molecule {
    let mut editor = MoleculeEditor::new();
    editor.add_atom(atom(symbol)).unwrap();
    editor.finish().unwrap()
}

fn classes(topology: &Topology) -> Vec<MoleculeClass> {
    topology.molecules().map(|m| m.class()).collect()
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
        let mut editor = source.edit();
        match operation {
            0 => editor.set_bond_order(bond, BondOrder::Single).unwrap(),
            1 => editor.bond_mut(bond).unwrap().set_order(BondOrder::Single),
            _ => {
                editor
                    .replace_bond(bond, Bond::new(b, a, BondOrder::Single))
                    .unwrap();
            }
        }
        editor.validate().unwrap();
        let result = editor.try_finish().unwrap();
        assert_eq!(result.bond(bond).unwrap().order, BondOrder::Single);
        assert_eq!(result.stereo_elements().count(), 0);
        let text = kekule::smiles::write_isomeric(&result).unwrap();
        let reparsed = molecule(&text);
        assert_eq!(reparsed.stereo_elements().count(), 0);
        assert_eq!(kekule::smiles::write_canonical(&result).unwrap(), "FCCF");
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
    let mut editor = source.edit();
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
    let mut editor = source.edit();
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
        let source = Arc::new(Topology::from_molecule(&molecule("F/C=C/F")).unwrap());
        let mut editor = source.edit();
        let (id, bond) = editor
            .bonds()
            .find(|(_, bond)| bond.order == BondOrder::Double)
            .unwrap();
        if replacement {
            editor
                .replace_bond(
                    id,
                    kekule::topology::EditBond::new(bond.a(), bond.b(), BondOrder::Single),
                )
                .unwrap();
        } else {
            editor.set_bond_order(id, BondOrder::Single).unwrap();
        }
        let target = editor.finish().unwrap();
        let definition = target.molecules().next().unwrap();
        assert_eq!(definition.molecule().stereo_elements().count(), 0);
        assert_eq!(
            source
                .molecules()
                .next()
                .unwrap()
                .molecule()
                .stereo_elements()
                .count(),
            1
        );
    }
}

fn split_editor() -> (TopologyEditor, Vec<EditAtomId>) {
    let mut editor = Topology::from_molecule(&molecule("CCCC"))
        .unwrap()
        .into_editor();
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    let middle = editor.bond_between(atoms[1], atoms[2]).unwrap().unwrap();
    editor.delete_bond(middle).unwrap();
    editor
        .set_molecule_class(atoms[0], MoleculeClass::Other)
        .unwrap();
    editor
        .set_molecule_class(atoms[3], MoleculeClass::Ion)
        .unwrap();
    (editor, atoms)
}

#[test]
fn component_overrides_survive_chemical_changes_to_unrelated_fragments() {
    for operation in 0..5 {
        let (mut editor, atoms) = split_editor();
        let bond = editor.bond_between(atoms[2], atoms[3]).unwrap().unwrap();
        match operation {
            0 => {
                editor.replace_atom(atoms[2], atom("N")).unwrap();
            }
            1 => {
                editor.set_bond_order(bond, BondOrder::Double).unwrap();
            }
            2 => {
                editor.delete_bond(bond).unwrap();
            }
            3 => {
                editor.delete_atom(atoms[2]).unwrap();
            }
            _ => {
                editor.delete_atom(atoms[3]).unwrap();
            }
        }
        let target = editor.finish().unwrap();
        let result = classes(&target);
        assert_eq!(result[0], MoleculeClass::Other, "operation {operation}");
        assert!(result[1..]
            .iter()
            .all(|class| *class == MoleculeClass::SmallMolecule));
    }
}

#[test]
fn merging_and_rewiring_invalidate_only_touched_components_in_historical_groups() {
    for rewire in [false, true] {
        let (mut editor, atoms) = split_editor();
        let added = editor.add_atom(atom("C")).unwrap();
        editor
            .set_molecule_class(added, MoleculeClass::Water)
            .unwrap();
        if rewire {
            let bond = editor.bond_between(atoms[2], atoms[3]).unwrap().unwrap();
            editor.set_bond_endpoints(bond, atoms[2], added).unwrap();
        } else {
            editor.add_bond(atoms[3], added, BondOrder::Single).unwrap();
        }
        let target = editor.finish().unwrap();
        let result = classes(&target);
        assert_eq!(result[0], MoleculeClass::Other);
        assert_eq!(result.len(), if rewire { 3 } else { 2 });
        assert!(result[1..]
            .iter()
            .all(|class| *class == MoleculeClass::SmallMolecule));
    }

    let mut editor = Topology::from_molecule(&molecule("CCCCCC"))
        .unwrap()
        .into_editor();
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    for (left, right) in [(1, 2), (3, 4)] {
        let bond = editor
            .bond_between(atoms[left], atoms[right])
            .unwrap()
            .unwrap();
        editor.delete_bond(bond).unwrap();
    }
    for (index, class) in [
        (0, MoleculeClass::Other),
        (2, MoleculeClass::Ion),
        (4, MoleculeClass::Water),
    ] {
        editor.set_molecule_class(atoms[index], class).unwrap();
    }
    editor
        .add_bond(atoms[3], atoms[4], BondOrder::Single)
        .unwrap();
    assert_eq!(
        classes(&editor.finish().unwrap()),
        vec![MoleculeClass::Other, MoleculeClass::SmallMolecule]
    );
}

#[test]
fn failed_and_noop_edits_retain_component_override_assignments() {
    let (mut editor, atoms) = split_editor();
    let bond = editor.bond_between(atoms[2], atoms[3]).unwrap().unwrap();
    let previous = editor.atom(atoms[2]).unwrap().clone();
    editor.replace_atom(atoms[2], previous).unwrap();
    editor.set_bond_order(bond, BondOrder::Single).unwrap();
    editor.set_bond_endpoints(bond, atoms[3], atoms[2]).unwrap();
    assert!(editor
        .add_bond(atoms[2], atoms[2], BondOrder::Single)
        .is_err());
    assert!(editor
        .add_bond(atoms[2], atoms[3], BondOrder::Single)
        .is_err());
    assert!(editor.set_bond_endpoints(bond, atoms[2], atoms[2]).is_err());
    assert_eq!(
        classes(&editor.finish().unwrap()),
        vec![MoleculeClass::Other, MoleculeClass::Ion]
    );
}

#[test]
fn published_hierarchy_removal_discards_orphaned_explicit_residue_assignments() {
    let carbon = single_atom("C");
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&carbon).unwrap();
    let atom = InstanceAtomId::new(instance, carbon.atom_ids().next().unwrap());
    let chain = builder.hierarchy_mut().add_chain("OLD", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "UNL", None, None, None)
        .unwrap();
    builder
        .hierarchy_mut()
        .add_atom_site(residue, atom, AtomSiteMetadata::default())
        .unwrap();
    builder
        .set_residue_class(residue, ResidueClass::Ion)
        .unwrap();
    let source = Arc::new(builder.clone().build().unwrap());

    for recoverable in [false, true] {
        let mut builder = builder.clone().build().unwrap().into_builder();
        *builder.hierarchy_mut() = Hierarchy::new();
        builder.validate().unwrap();
        let removed = if recoverable {
            builder.try_build().unwrap()
        } else {
            builder.build().unwrap()
        };
        assert_eq!(removed.residues().count(), 0);
        let mut builder = removed.into_builder();
        let chain = builder.hierarchy_mut().add_chain("NEW", None).unwrap();
        let reused = builder
            .hierarchy_mut()
            .add_residue(chain, "UNL", None, None, None)
            .unwrap();
        assert_eq!(reused, residue);
        builder
            .hierarchy_mut()
            .add_atom_site(reused, atom, AtomSiteMetadata::default())
            .unwrap();
        let result = builder.build().unwrap();
        assert_eq!(
            result.residues().next().unwrap().class(),
            ResidueClass::Other
        );
        assert_eq!(classes(&result), vec![MoleculeClass::SmallMolecule]);
    }
    assert_eq!(source.residues().next().unwrap().class(), ResidueClass::Ion);
}

#[test]
fn builder_assignments_still_apply_to_surviving_residues_after_hierarchy_replacement() {
    let oxygen = single_atom("O");
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&oxygen).unwrap();
    let chain = builder.hierarchy_mut().add_chain("OLD", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "UNL", None, None, None)
        .unwrap();
    builder
        .set_residue_class(residue, ResidueClass::Ion)
        .unwrap();
    *builder.hierarchy_mut() = Hierarchy::new();
    let chain = builder.hierarchy_mut().add_chain("NEW", None).unwrap();
    let replaced = builder
        .hierarchy_mut()
        .add_residue(chain, "HOH", None, None, None)
        .unwrap();
    assert_eq!(replaced, residue);
    builder
        .hierarchy_mut()
        .add_atom_site(
            replaced,
            InstanceAtomId::new(instance, oxygen.atom_ids().next().unwrap()),
            AtomSiteMetadata::default(),
        )
        .unwrap();
    let result = builder.build().unwrap();
    assert_eq!(result.residues().next().unwrap().class(), ResidueClass::Ion);
    assert_eq!(
        result
            .into_builder()
            .build()
            .unwrap()
            .residues()
            .next()
            .unwrap()
            .class(),
        ResidueClass::Ion
    );
}

#[test]
fn instance_filters_reclassify_partial_residues_and_preserve_complete_assignments() {
    for explicit in [false, true] {
        let mut builder = TopologyBuilder::new();
        let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, "UNL", None, None, None)
            .unwrap();
        let mut atoms = Vec::new();
        for symbol in ["O", "H", "H"] {
            let molecule = single_atom(symbol);
            let instance = builder.add_molecule(&molecule).unwrap();
            let atom = InstanceAtomId::new(instance, molecule.atom_ids().next().unwrap());
            builder
                .hierarchy_mut()
                .add_atom_site(residue, atom, AtomSiteMetadata::default())
                .unwrap();
            atoms.push(atom);
        }
        if explicit {
            builder
                .set_residue_class(residue, ResidueClass::Water)
                .unwrap();
        }
        let whole = builder
            .hierarchy_mut()
            .add_residue(chain, "UNL", None, None, None)
            .unwrap();
        let nitrogen = single_atom("N");
        let instance = builder.add_molecule(&nitrogen).unwrap();
        let nitrogen_atom = InstanceAtomId::new(instance, nitrogen.atom_ids().next().unwrap());
        builder
            .hierarchy_mut()
            .add_atom_site(whole, nitrogen_atom, AtomSiteMetadata::default())
            .unwrap();
        builder
            .set_residue_class(whole, ResidueClass::Carbohydrate)
            .unwrap();
        let source = Arc::new(builder.build().unwrap());
        assert_eq!(
            source.residues().next().unwrap().class(),
            ResidueClass::Water
        );
        let retained =
            transform::retain_instances(&source, [atoms[0].molecule(), instance]).unwrap();
        let removed =
            transform::remove_instances(&source, [atoms[1].molecule(), atoms[2].molecule()])
                .unwrap();
        let subset = source
            .subset(&AtomSelection::from_atoms(&source, [atoms[0], nitrogen_atom]).unwrap())
            .unwrap();
        for result in [&*retained, &*removed, subset.topology()] {
            assert_eq!(
                result.residues().map(|r| r.class()).collect::<Vec<_>>(),
                vec![ResidueClass::Other, ResidueClass::Carbohydrate]
            );
        }
        let mut editor = retained.edit();
        let whole = editor.residue_handle(whole).unwrap();
        editor
            .set_residue_component_ids(whole, Some("UNK".into()), None)
            .unwrap();
        assert_eq!(
            editor.finish().unwrap().residues().nth(1).unwrap().class(),
            ResidueClass::Carbohydrate
        );
    }
}
