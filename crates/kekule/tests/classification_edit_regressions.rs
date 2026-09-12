use std::sync::Arc;

use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor};
use kekule::topology::{
    transform, AtomSelection, AtomSiteMetadata, EditAtomId, Hierarchy, InstanceAtomId,
    MoleculeClass, MoleculeDefinitionId, ResidueClass, Topology, TopologyBuilder, TopologyEditor,
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

fn oxygen_builder() -> (TopologyBuilder, MoleculeDefinitionId) {
    let mut editor = MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder
        .add_molecule_definition_owned(editor.finish().unwrap())
        .unwrap();
    (builder, definition)
}

fn add_instance(
    builder: &mut TopologyBuilder,
    definition: MoleculeDefinitionId,
    component: Option<&str>,
    residue_class: Option<ResidueClass>,
) -> kekule::topology::MoleculeInstanceId {
    let local = builder
        .definition(definition)
        .unwrap()
        .molecule()
        .atom_ids()
        .next()
        .unwrap();
    let instance = builder.add_instance(definition).unwrap();
    if let Some(component) = component {
        let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, component, None, None, None)
            .unwrap();
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, local),
                AtomSiteMetadata::default(),
            )
            .unwrap();
        if let Some(class) = residue_class {
            builder.set_residue_class(residue, class).unwrap();
        }
    }
    instance
}

#[test]
fn resumed_reused_definitions_combine_new_informative_evidence_like_fresh_builders() {
    for explicit in [false, true] {
        let (mut builder, definition) = oxygen_builder();
        add_instance(&mut builder, definition, Some("HOH"), None);
        if explicit {
            builder
                .set_molecule_class(definition, MoleculeClass::Water)
                .unwrap();
        }
        let source = builder.clone().build().unwrap();
        assert_eq!(
            source.definition(definition).unwrap().class(),
            MoleculeClass::Water
        );
        for mut builder in [source.into_builder(), builder] {
            add_instance(
                &mut builder,
                definition,
                Some("UNK"),
                Some(ResidueClass::Ion),
            );
            builder.validate().unwrap();
            let result = builder.try_build().unwrap();
            assert_eq!(
                result.definition(definition).unwrap().class(),
                if explicit {
                    MoleculeClass::Water
                } else {
                    MoleculeClass::Other
                }
            );
            assert_eq!(
                result.residues().map(|r| r.class()).collect::<Vec<_>>(),
                vec![ResidueClass::Water, ResidueClass::Ion]
            );
        }
    }
    let (mut builder, definition) = oxygen_builder();
    add_instance(&mut builder, definition, None, None);
    let mut builder = builder.build().unwrap().into_builder();
    add_instance(&mut builder, definition, Some("HOH"), None);
    assert_eq!(
        builder
            .build()
            .unwrap()
            .definition(definition)
            .unwrap()
            .class(),
        MoleculeClass::Water
    );
}

#[test]
fn uninformative_and_unrelated_appends_preserve_complete_entity_cached_classes() {
    let (mut builder, definition) = oxygen_builder();
    let retained = add_instance(&mut builder, definition, Some("HOH"), None);
    add_instance(
        &mut builder,
        definition,
        Some("UNK"),
        Some(ResidueClass::Ion),
    );
    let source = Arc::new(builder.build().unwrap());
    // The retained class deliberately remembers the complete source definition.
    // Adding an uninformative occurrence must not recompute it from just HOH.
    for component in [None, Some("UNK")] {
        let retained = transform::retain_instances(&source, [retained]).unwrap();
        assert_eq!(
            retained.definition(definition).unwrap().class(),
            MoleculeClass::Other
        );
        let mut builder = Arc::try_unwrap(retained).unwrap().into_builder();
        add_instance(&mut builder, definition, component, None);
        let independent = builder.add_molecule_definition(&molecule("CC")).unwrap();
        builder.add_instance(independent).unwrap();
        let result = builder.build().unwrap();
        assert_eq!(
            result.definition(definition).unwrap().class(),
            MoleculeClass::Other
        );
        assert_eq!(
            result.definition(independent).unwrap().class(),
            MoleculeClass::SmallMolecule
        );
    }
}
