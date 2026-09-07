use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor};
use kekule::geometry::Point3;
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{Model, ModelEditor, Positions};
use kekule::topology::{
    AtomSiteMetadata, InstanceAtomId, MoleculeClass, Topology, TopologyBuilder, TopologyEditor,
};
use kekule::units::{Quantity, ANGSTROM, NANOMETER};
use std::sync::Arc;

fn molecule(text: &str) -> Molecule {
    kekule::smiles::to_molecules(text).unwrap().pop().unwrap()
}
fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}
fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}
fn positions(xs: &[f64]) -> Positions {
    Positions::new(Quantity::new(
        xs.iter()
            .map(|&x| Point3::new(x, 0.0, 0.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap()
}
fn model(text: &str) -> Model {
    let mut builder = Model::builder();
    for molecule in kekule::smiles::to_molecules(text).unwrap() {
        let start = builder.atom_count();
        builder
            .add_molecule(
                &molecule,
                &positions(
                    &(start..start + molecule.atom_count())
                        .map(|i| i as f64)
                        .collect::<Vec<_>>(),
                ),
            )
            .unwrap();
    }
    builder.build().unwrap()
}

#[test]
fn split_then_merge_preserves_handles_geometry_and_hierarchy() {
    let molecule = molecule("CCC");
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "UNL", Some(1), None, None)
        .unwrap();
    for atom in molecule.atom_ids() {
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, atom),
                AtomSiteMetadata::default(),
            )
            .unwrap();
    }
    let source = Model::new(builder.build().unwrap(), positions(&[1.0, 2.0, 3.0])).unwrap();
    let mut editor = source.edit();
    let handles = source
        .topology()
        .atom_ids()
        .iter()
        .map(|&id| editor.atom_handle(id).unwrap())
        .collect::<Vec<_>>();
    let removed = editor.bond_handle(source.topology().bond_ids()[0]).unwrap();
    editor.delete_bond(removed).unwrap();
    editor.validate().unwrap();
    let split = editor.clone().finish_with_correspondence().unwrap();
    assert_eq!(split.model().topology().instance_count(), 2);
    assert_eq!(
        split
            .correspondence()
            .target_instances(instance)
            .unwrap()
            .len(),
        2
    );
    assert_eq!(split.model().hierarchy().chains().count(), 1);
    assert_eq!(split.model().hierarchy().residues().count(), 1);
    assert_eq!(split.model().hierarchy().atom_sites().count(), 3);
    for (&handle, &id) in handles.iter().zip(source.topology().atom_ids()) {
        assert_eq!(
            split
                .model()
                .position(split.correspondence().atom(handle).unwrap())
                .unwrap(),
            source.position(id).unwrap()
        );
    }
    editor
        .add_bond(handles[0], handles[2], BondOrder::Single)
        .unwrap();
    let merged = editor.finish_with_correspondence().unwrap();
    assert_eq!(merged.model().topology().instance_count(), 1);
    assert_eq!(
        merged
            .correspondence()
            .target_instances(instance)
            .unwrap()
            .len(),
        1
    );
    assert!(merged.correspondence().bond(removed).is_none());
    assert_eq!(source.topology().bond_count(), 2);
}

#[test]
fn cross_instance_bond_merges_but_preserves_distinct_chains() {
    let molecule = molecule("C");
    let mut builder = TopologyBuilder::new();
    for label in ["A", "B"] {
        let instance = builder.add_molecule(&molecule).unwrap();
        let chain = builder.hierarchy_mut().add_chain(label, None).unwrap();
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, "UNL", None, None, None)
            .unwrap();
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, molecule.atom_ids().next().unwrap()),
                AtomSiteMetadata::default(),
            )
            .unwrap();
    }
    builder
        .molecule_instance_properties_mut()
        .insert(key("instance"), PropertyColumn::Int(vec![Some(1), Some(2)]))
        .unwrap();
    let source = Model::new(builder.build().unwrap(), positions(&[2.0, 6.0])).unwrap();
    let mut editor = source.edit();
    let handles = editor.atom_ids().collect::<Vec<_>>();
    let bond = editor
        .add_bond(handles[0], handles[1], BondOrder::Single)
        .unwrap();
    let result = editor.finish_with_correspondence().unwrap();
    assert_eq!(result.model().topology().instance_count(), 1);
    assert_eq!(result.model().topology().bond_count(), 1);
    assert_eq!(result.model().hierarchy().chains().count(), 2);
    assert_eq!(result.model().hierarchy().residues().count(), 2);
    assert!(!result
        .model()
        .topology()
        .molecule_instance_properties()
        .has_data());
    assert!(result.correspondence().bond(bond).is_some());
    for &source_atom in source.topology().atom_ids() {
        let target = result.correspondence().target_atom(source_atom).unwrap();
        assert_eq!(
            source.position(source_atom).unwrap(),
            result.model().position(target).unwrap()
        );
        assert_eq!(
            result
                .correspondence()
                .target_instances(source_atom.molecule())
                .unwrap(),
            &[target.molecule()]
        );
    }
}

#[test]
fn editing_one_reused_occurrence_preserves_others_and_their_perception() {
    let mut water = molecule("O");
    water.perceive().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&water).unwrap();
    for _ in 0..3 {
        builder.add_instance(definition).unwrap();
    }
    let topology = Arc::new(builder.build().unwrap());
    let source_atoms = topology.atom_ids().to_vec();
    let mut editor = topology.edit();
    let changed = editor.atom_handle(source_atoms[1]).unwrap();
    editor.replace_atom(changed, atom("N")).unwrap();
    let result = editor.finish_with_correspondence().unwrap();
    assert_eq!(result.topology().definition_count(), 2);
    let untouched =
        [source_atoms[0], source_atoms[2]].map(|a| result.correspondence().target_atom(a).unwrap());
    assert_eq!(
        result
            .topology()
            .instance(untouched[0].molecule())
            .unwrap()
            .definition(),
        result
            .topology()
            .instance(untouched[1].molecule())
            .unwrap()
            .definition()
    );
    for atom in untouched {
        assert_eq!(result.topology().atom(atom).unwrap().element.symbol(), "O");
        assert_eq!(
            result
                .topology()
                .definition_for_instance(atom.molecule())
                .unwrap()
                .molecule()
                .perception(),
            water.perception()
        );
    }
    let changed = result.correspondence().atom(changed).unwrap();
    assert_eq!(
        result.topology().atom(changed).unwrap().element.symbol(),
        "N"
    );
    assert_ne!(
        result
            .topology()
            .definition_for_instance(changed.molecule())
            .unwrap()
            .class(),
        MoleculeClass::Water
    );
    assert_eq!(topology.definition_count(), 1);
    assert_eq!(
        topology.atom(source_atoms[1]).unwrap().element.symbol(),
        "O"
    );
}

#[test]
fn no_op_and_geometry_only_keep_exact_source_topology() {
    let mut source = model("CO");
    source.perceive().unwrap();
    source
        .insert_property(key("label"), PropertyValue::Int(9))
        .unwrap();
    let topology = source.shared_topology();
    let no_op = source.edit().finish().unwrap();
    assert_eq!(no_op, source);
    assert!(Arc::ptr_eq(&no_op.shared_topology(), &topology));
    let mut editor = source.edit();
    let id = editor.atom_ids().next().unwrap();
    editor
        .replace_atom(id, editor.atom(id).unwrap().clone())
        .unwrap();
    editor
        .set_position(id, Quantity::new(Point3::new(3.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    let result = editor.finish_with_correspondence().unwrap();
    assert!(Arc::ptr_eq(&result.model().shared_topology(), &topology));
    assert_eq!(
        result.model().properties().get(&key("label")),
        source.properties().get(&key("label"))
    );
    assert_eq!(
        result
            .model()
            .position(result.correspondence().atom(id).unwrap())
            .unwrap()
            .value()
            .x,
        3.0
    );
    assert_ne!(result.model().positions(), source.positions());
}

#[test]
fn deletion_prunes_hierarchy_and_projects_properties_and_coordinates() {
    let mut builder = TopologyBuilder::new();
    let molecule = molecule("C");
    for label in ["A", "B"] {
        let instance = builder.add_molecule(&molecule).unwrap();
        let chain = builder.hierarchy_mut().add_chain(label, None).unwrap();
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, "UNL", None, None, None)
            .unwrap();
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, molecule.atom_ids().next().unwrap()),
                AtomSiteMetadata::default(),
            )
            .unwrap();
    }
    builder
        .atom_properties_mut()
        .insert(key("static"), PropertyColumn::Int(vec![Some(1), Some(2)]))
        .unwrap();
    let mut source = Model::new(builder.build().unwrap(), positions(&[4.0, 9.0])).unwrap();
    let source_atoms = source.topology().atom_ids().to_vec();
    source.set_occupancy(source_atoms[1], Some(0.7)).unwrap();
    source
        .insert_property(key("energy"), PropertyValue::Int(10))
        .unwrap();
    let mut editor = source.edit();
    let removed = editor.atom_handle(source_atoms[0]).unwrap();
    editor.delete_atom(removed).unwrap();
    assert!(editor.position(removed).is_err());
    let new = editor
        .add_atom(
            atom("O"),
            Quantity::new(Point3::new(5.0, 0.0, 0.0), ANGSTROM),
        )
        .unwrap();
    editor
        .insert_atom_property_column(key("live"), PropertyColumn::Int(vec![Some(7), Some(8)]))
        .unwrap();
    let column = editor.atom_property_column(&key("live")).unwrap().unwrap();
    editor
        .insert_atom_property_column(key("live"), column)
        .unwrap();
    let result = editor.finish_with_correspondence().unwrap();
    assert_eq!(result.model().hierarchy().chains().count(), 1);
    assert_eq!(result.model().hierarchy().residues().count(), 1);
    assert_eq!(result.model().hierarchy().atom_sites().count(), 1);
    assert!(result.model().properties().owner_is_empty());
    let retained = result
        .correspondence()
        .target_atom(source_atoms[1])
        .unwrap();
    assert_eq!(
        result.model().position(retained).unwrap(),
        source.position(source_atoms[1]).unwrap()
    );
    assert_eq!(result.model().occupancy(retained).unwrap(), Some(0.7));
    assert_eq!(
        result
            .model()
            .topology()
            .atom_property(retained, &key("static"))
            .unwrap(),
        Some(PropertyValue::Int(2))
    );
    let new = result.correspondence().atom(new).unwrap();
    assert_eq!(result.model().occupancy(new).unwrap(), None);
    assert_eq!(
        result.model().atom_property(new, &key("live")).unwrap(),
        Some(PropertyValue::Int(8))
    );
    assert_eq!(
        result
            .correspondence()
            .source_atom_indices()
            .iter()
            .filter(|v| v.is_none())
            .count(),
        1
    );
}

#[test]
fn failures_are_atomic_and_empty_publication_is_recoverable() {
    let source = model("CC");
    let mut editor = source.edit();
    let handles = editor.atom_ids().collect::<Vec<_>>();
    let mut foreign = ModelEditor::new();
    let invalid = foreign
        .add_atom(atom("O"), Quantity::new(Point3::origin(), ANGSTROM))
        .unwrap();
    let before = format!("{editor:?}");
    assert!(editor.delete_atoms([handles[0], invalid]).is_err());
    assert!(editor
        .add_bond(handles[0], handles[1], BondOrder::Single)
        .is_err());
    assert!(editor
        .add_atom(
            atom("O"),
            Quantity::new(Point3::new(f64::NAN, 0.0, 0.0), ANGSTROM)
        )
        .is_err());
    assert!(editor
        .set_atom_positions([
            (handles[0], Quantity::new(Point3::origin(), ANGSTROM)),
            (invalid, Quantity::new(Point3::origin(), ANGSTROM))
        ])
        .is_err());
    assert_eq!(format!("{editor:?}"), before);
    editor.delete_atoms(handles).unwrap();
    let before = format!("{editor:?}");
    let failure = editor.try_finish().unwrap_err();
    assert_eq!(format!("{:?}", failure.editor()), before);
    let mut editor = failure.into_editor();
    editor
        .add_atom(atom("O"), Quantity::new(Point3::origin(), ANGSTROM))
        .unwrap();
    assert_eq!(editor.finish().unwrap().atom_count(), 1);
}

#[test]
fn append_preserves_sparse_local_ids_dense_order_and_definition_order() {
    let mut molecule = molecule("CCC").into_editor();
    let removed = molecule.atom_ids().next().unwrap();
    molecule.delete_atom(removed).unwrap();
    let source =
        Model::from_molecule(&molecule.finish().unwrap(), &positions(&[3.0, 7.0])).unwrap();
    let mut editor = source.edit();
    let added = editor
        .add_molecule(&self::molecule("O"), &positions(&[8.0]))
        .unwrap();
    let result = editor.finish_with_correspondence().unwrap();
    assert_eq!(
        &result.model().topology().atom_ids()[..source.atom_count()],
        source.topology().atom_ids()
    );
    assert_eq!(
        &result.model().topology().bond_ids()[..source.topology().bond_count()],
        source.topology().bond_ids()
    );
    for &id in source.topology().atom_ids() {
        assert_eq!(result.correspondence().target_atom(id), Some(id));
        assert_eq!(
            result.model().position(id).unwrap(),
            source.position(id).unwrap()
        );
    }
    assert!(result
        .correspondence()
        .atom(*added.atoms().values().next().unwrap())
        .is_some());
}

#[test]
fn merge_property_conflicts_roll_back_and_stereo_survives_unrelated_append() {
    let mut left = molecule("C");
    let id = left.atom_ids().next().unwrap();
    left.set_atom_property(id, key("tag"), Some(PropertyValue::Int(3)))
        .unwrap();
    let mut right = molecule("C");
    let id = right.atom_ids().next().unwrap();
    right
        .set_atom_property(id, key("tag"), Some(PropertyValue::String("x".into())))
        .unwrap();
    let mut editor = TopologyEditor::new();
    let left = *editor
        .add_molecule(&left)
        .unwrap()
        .atoms()
        .values()
        .next()
        .unwrap();
    let right = *editor
        .add_molecule(&right)
        .unwrap()
        .atoms()
        .values()
        .next()
        .unwrap();
    let before = format!("{editor:?}");
    assert!(editor.add_bond(left, right, BondOrder::Single).is_err());
    assert_eq!(format!("{editor:?}"), before);
    let stereo = molecule("F[C@](Cl)(Br)I");
    let mut editor = Topology::from_molecule(&stereo).unwrap().into_editor();
    editor.add_molecule(&molecule("O")).unwrap();
    let topology = editor.finish().unwrap();
    assert_eq!(
        topology.definitions().next().unwrap().1.molecule().graph(),
        stereo.graph()
    );
}

#[test]
fn unchanged_molecule_publication_keeps_perception_and_columns_round_trip() {
    let mut source = molecule("CCC");
    source.perceive().unwrap();
    assert_eq!(
        source.edit().finish().unwrap().perception(),
        source.perception()
    );
    let mut editor = source.into_editor();
    let deleted = editor.atom_ids().next().unwrap();
    editor.delete_atom(deleted).unwrap();
    editor
        .insert_atom_property_column(key("tag"), PropertyColumn::Int(vec![Some(1), Some(2)]))
        .unwrap();
    let column = editor.atom_property_column(&key("tag")).unwrap().unwrap();
    assert_eq!(column.len(), 2);
    assert_eq!(
        editor
            .insert_atom_property_column(key("tag"), column.clone())
            .unwrap(),
        Some(column.clone())
    );
    assert_eq!(
        editor.remove_atom_property_column(&key("tag")),
        Some(column)
    );
    assert!(MoleculeEditor::new().finish().is_err());
}

#[test]
fn classification_overrides_are_occurrence_and_component_local() {
    let carbon = molecule("CC");
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&carbon).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    let source = Arc::new(builder.build().unwrap());
    let mut editor = source.edit();
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    editor
        .set_molecule_class(atoms[0], MoleculeClass::Other)
        .unwrap();
    let classified = editor.clone().finish_with_correspondence().unwrap();
    let class = |result: &kekule::topology::TopologyEdit, atom| {
        result
            .topology()
            .definition_for_instance(result.correspondence().atom(atom).unwrap().molecule())
            .unwrap()
            .class()
    };
    assert_eq!(class(&classified, atoms[0]), MoleculeClass::Other);
    assert_eq!(class(&classified, atoms[2]), MoleculeClass::SmallMolecule);
    let bond = editor.bond_ids().next().unwrap();
    editor.delete_bond(bond).unwrap();
    editor
        .set_molecule_class(atoms[0], MoleculeClass::Ion)
        .unwrap();
    editor
        .set_molecule_class(atoms[1], MoleculeClass::Other)
        .unwrap();
    let split = editor.clone().finish_with_correspondence().unwrap();
    assert_eq!(class(&split, atoms[0]), MoleculeClass::Ion);
    assert_eq!(class(&split, atoms[1]), MoleculeClass::Other);
    assert_eq!(class(&split, atoms[2]), MoleculeClass::SmallMolecule);
    editor
        .add_bond(atoms[0], atoms[1], BondOrder::Single)
        .unwrap();
    let restored = editor.finish_with_correspondence().unwrap();
    assert_eq!(class(&restored, atoms[0]), MoleculeClass::SmallMolecule);
}

#[test]
fn rewiring_preserves_bond_identity_and_all_three_annotation_scopes() {
    let mut source = model("CC.O");
    let source_bond = source.bond_ids()[0];
    source
        .set_bond_property(source_bond, key("dynamic"), Some(PropertyValue::Int(10)))
        .unwrap();
    let mut editor = source.edit();
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    let bond = editor.bond_handle(source_bond).unwrap();
    editor
        .set_topology_bond_property(bond, key("static"), Some(PropertyValue::Int(20)))
        .unwrap();
    editor
        .set_definition_bond_property(bond, key("definition"), Some(PropertyValue::Int(30)))
        .unwrap();
    let before = editor.bond(bond).unwrap();
    assert!(editor.set_bond_endpoints(bond, atoms[2], atoms[2]).is_err());
    assert_eq!(editor.bond(bond).unwrap(), before);
    editor.set_bond_endpoints(bond, atoms[0], atoms[2]).unwrap();
    let result = editor.finish_with_correspondence().unwrap();
    let target = result.correspondence().bond(bond).unwrap();
    assert_eq!(
        result.correspondence().target_bond(source_bond),
        Some(target)
    );
    assert_eq!(result.model().topology().instance_count(), 2);
    let actual = result.model().bond(target).unwrap();
    let a = result.correspondence().atom(atoms[0]).unwrap();
    let b = result.correspondence().atom(atoms[2]).unwrap();
    assert_eq!(actual.endpoints(), (a.atom(), b.atom()));
    assert_eq!(
        result
            .model()
            .bond_property(target, &key("dynamic"))
            .unwrap(),
        Some(PropertyValue::Int(10))
    );
    assert_eq!(
        result
            .model()
            .topology()
            .bond_property(target, &key("static"))
            .unwrap(),
        Some(PropertyValue::Int(20))
    );
    assert_eq!(
        result
            .model()
            .topology()
            .definition_for_instance(target.molecule())
            .unwrap()
            .molecule()
            .bond_property(target.bond(), &key("definition"))
            .unwrap(),
        Some(PropertyValue::Int(30))
    );
    for &id in source.atom_ids() {
        assert_eq!(
            source.position(id).unwrap(),
            result
                .model()
                .position(result.correspondence().target_atom(id).unwrap())
                .unwrap()
        );
    }
}

#[test]
fn recovery_keeps_builder_state_and_editor_handles() {
    use std::error::Error;
    let carbon = molecule("C");
    let mut builder = Model::builder();
    let definition = builder
        .add_molecule_definition_owned(carbon.clone())
        .unwrap();
    let failed = builder.try_build().unwrap_err();
    assert!(failed.source().is_some());
    assert_eq!(failed.builder().topology_builder().definition_count(), 1);
    let mut builder = failed.into_builder();
    let instance = builder
        .add_instance(definition, &positions(&[3.0]))
        .unwrap();
    let id = InstanceAtomId::new(instance, carbon.atom_ids().next().unwrap());
    builder.set_occupancy(id, Some(0.5)).unwrap();
    builder
        .set_atom_positions([(id, Quantity::new(Point3::new(0.4, 0.0, 0.0), NANOMETER))])
        .unwrap();
    assert_eq!(
        builder
            .position(id)
            .unwrap()
            .into_unit(ANGSTROM)
            .unwrap()
            .value()
            .x,
        4.0
    );
    assert_eq!(builder.occupancy(id).unwrap(), Some(0.5));
    builder
        .set_atom_properties(key("tag"), [(id, Some(PropertyValue::Int(7)))])
        .unwrap();
    assert_eq!(
        builder.atom_property(id, &key("tag")).unwrap(),
        Some(PropertyValue::Int(7))
    );
    assert_eq!(builder.atom_ids().collect::<Vec<_>>(), vec![id]);
    builder.validate().unwrap();
    let source = builder.build().unwrap();
    let mut editor = source.into_editor();
    let handle = editor.atom_handle(id).unwrap();
    editor.delete_atom(handle).unwrap();
    let failed = editor.try_finish_with_correspondence().unwrap_err();
    assert!(failed.source().is_some());
    let mut editor = failed.into_editor();
    let replacement = editor
        .add_atom(atom("N"), Quantity::new(Point3::origin(), ANGSTROM))
        .unwrap();
    assert_ne!(handle, replacement);
    assert!(editor.atom(handle).is_err());
    let result = editor.finish_with_correspondence().unwrap();
    assert!(result.correspondence().target_atom(id).is_none());
    assert_eq!(
        result.correspondence().target_instances(instance),
        Some([].as_slice())
    );
    assert!(result.correspondence().atom(replacement).is_some());

    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&carbon).unwrap();
    let failed = builder.try_build().unwrap_err();
    assert_eq!(failed.builder().definition_count(), 1);
    let mut builder = failed.into_builder();
    builder.add_instance(definition).unwrap();
    builder.validate().unwrap();
    assert_eq!(builder.build().unwrap().atom_count(), 1);
}

#[test]
fn resumed_model_builder_retains_entity_state_and_extends_missing_rows() {
    let mut draft = molecule("NCC").into_editor();
    let removed = draft.atom_ids().next().unwrap();
    draft.delete_atom(removed).unwrap();
    let mut carbon = draft.finish().unwrap();
    carbon.perceive().unwrap();
    let mut builder = Model::builder();
    let definition = builder.add_molecule_definition(&carbon).unwrap();
    let first = builder
        .add_instance(definition, &positions(&[2.0, 3.0]))
        .unwrap();
    builder
        .add_instance(definition, &positions(&[5.0, 6.0]))
        .unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "UNL", None, None, None)
        .unwrap();
    let a = InstanceAtomId::new(first, carbon.atom_ids().next().unwrap());
    builder
        .hierarchy_mut()
        .add_atom_site(residue, a, AtomSiteMetadata::default())
        .unwrap();
    builder
        .set_atom_property(a, key("dynamic"), Some(PropertyValue::Int(4)))
        .unwrap();
    builder
        .topology_builder_mut()
        .atom_properties_mut()
        .set_value(key("static"), 0, Some(PropertyValue::Int(8)))
        .unwrap();
    builder
        .insert_property(key("old_owner"), PropertyValue::Int(1))
        .unwrap();
    builder
        .topology_builder_mut()
        .insert_property(key("old_owner"), PropertyValue::Int(2))
        .unwrap();
    let source = builder.build().unwrap();
    let unchanged = source.to_builder().build().unwrap();
    assert_eq!(unchanged.properties(), source.properties());
    assert_eq!(
        unchanged.topology().properties(),
        source.topology().properties()
    );
    let mut extended = source.clone().into_builder();
    extended
        .add_instance(definition, &positions(&[8.0, 9.0]))
        .unwrap();
    let result = extended.build().unwrap();
    assert_eq!(&result.atom_ids()[..4], source.atom_ids());
    assert_eq!(&result.bond_ids()[..2], source.bond_ids());
    assert_eq!(
        result.positions().values().value()[..4],
        source.positions().values().value()[..]
    );
    assert_eq!(result.topology().definition_count(), 1);
    assert_eq!(result.hierarchy(), source.hierarchy());
    assert_eq!(
        result
            .topology()
            .definition(definition)
            .unwrap()
            .molecule()
            .perception(),
        carbon.perception()
    );
    assert_eq!(
        result.atom_property(a, &key("dynamic")).unwrap(),
        Some(PropertyValue::Int(4))
    );
    assert_eq!(
        result.topology().atom_property(a, &key("static")).unwrap(),
        Some(PropertyValue::Int(8))
    );
    assert_eq!(
        result
            .atom_property(result.atom_ids()[4], &key("dynamic"))
            .unwrap(),
        None
    );
    assert!(result.properties().get(&key("old_owner")).is_none());
    assert!(result
        .topology()
        .properties()
        .get(&key("old_owner"))
        .is_none());
}

#[test]
fn sparse_batches_are_atomic_and_deleted_rows_do_not_poison_new_values() {
    use std::error::Error;
    let mut source = model("CCC");
    let ids = source.atom_ids().to_vec();
    let before = source.positions().clone();
    assert!(source
        .set_atom_positions([
            (ids[0], Quantity::new(Point3::new(4.0, 0.0, 0.0), ANGSTROM)),
            (
                ids[1],
                Quantity::new(Point3::new(f64::NAN, 0.0, 0.0), ANGSTROM)
            ),
        ])
        .is_err());
    assert_eq!(source.positions(), &before);
    let error = source
        .set_atom_properties(
            key("tag"),
            [
                (ids[0], Some(PropertyValue::Int(1))),
                (ids[1], Some(PropertyValue::String("bad".into()))),
            ],
        )
        .unwrap_err();
    assert!(error.source().is_some());
    assert!(source.atom_property_column(&key("tag")).is_none());
    let mut editor = source.into_editor();
    let handles = editor.atom_ids().collect::<Vec<_>>();
    editor
        .set_atom_property(handles[0], key("tag"), Some(PropertyValue::Int(5)))
        .unwrap();
    editor.delete_atom(handles[0]).unwrap();
    editor
        .set_atom_property(
            handles[1],
            key("tag"),
            Some(PropertyValue::String("valid".into())),
        )
        .unwrap();
    let snapshot = editor.positions();
    assert!(editor
        .set_atom_positions([
            (handles[1], Quantity::new(Point3::origin(), ANGSTROM)),
            (handles[0], Quantity::new(Point3::origin(), ANGSTROM)),
        ])
        .is_err());
    assert_eq!(editor.positions(), snapshot);
    let result = editor.finish_with_correspondence().unwrap();
    assert_eq!(
        result
            .model()
            .atom_property(
                result.correspondence().atom(handles[1]).unwrap(),
                &key("tag")
            )
            .unwrap(),
        Some(PropertyValue::String("valid".into()))
    );
}

#[test]
fn unchanged_setters_keep_topology_snapshot_and_perception() {
    let mut carbon = molecule("CC");
    carbon.perceive().unwrap();
    let bond = carbon.bond_ids().next().unwrap();
    let mut molecule_editor = carbon.edit();
    molecule_editor
        .set_bond_order(bond, BondOrder::Single)
        .unwrap();
    assert_eq!(
        molecule_editor.finish().unwrap().perception(),
        carbon.perception()
    );
    let source = Arc::new(Topology::from_molecule(&carbon).unwrap());
    let mut editor = source.edit();
    let a = editor.atom_ids().next().unwrap();
    let b = editor.bond_ids().next().unwrap();
    editor
        .replace_atom(a, editor.atom(a).unwrap().clone())
        .unwrap();
    editor.set_bond_order(b, BondOrder::Single).unwrap();
    editor.set_atom_property(a, key("absent"), None).unwrap();
    editor
        .set_definition_atom_property(a, key("absent"), None)
        .unwrap();
    editor
        .set_definition_bond_property(b, key("absent"), None)
        .unwrap();
    assert!(Arc::ptr_eq(&source, &editor.finish().unwrap()));
}

#[test]
fn split_publication_preserves_surviving_stereo_and_prunes_changed_centers() {
    let source = model("F[C@](Cl)(Br)CC");
    let mut editor = source.edit();
    let terminal = editor
        .atoms()
        .find(|(id, a)| a.element.symbol() == "C" && editor.neighbors(*id).unwrap().count() == 1)
        .unwrap()
        .0;
    let cut = editor.incident_bonds(terminal).unwrap().next().unwrap();
    editor.delete_bond(cut).unwrap();
    let result = editor.clone().finish().unwrap();
    assert_eq!(
        result
            .topology()
            .definitions()
            .map(|(_, d)| d.molecule().stereo_elements().count())
            .sum::<usize>(),
        1
    );
    let fluorine = editor
        .atoms()
        .find(|(_, a)| a.element.symbol() == "F")
        .unwrap()
        .0;
    editor.delete_atom(fluorine).unwrap();
    let result = editor.finish().unwrap();
    assert_eq!(
        result
            .topology()
            .definitions()
            .map(|(_, d)| d.molecule().stereo_elements().count())
            .sum::<usize>(),
        0
    );
}

#[test]
fn independently_created_or_cleared_drafts_never_accept_foreign_handles() {
    let source = model("C");
    let mut first = source.edit();
    let mut second = source.edit();
    let a = first.atom_ids().next().unwrap();
    let b = second.atom_ids().next().unwrap();
    assert_ne!(a, b);
    assert!(first.delete_atom(b).is_err());
    assert!(second.delete_atom(a).is_err());
    let mut branch = first.clone();
    let c = first
        .add_atom(atom("O"), Quantity::new(Point3::origin(), ANGSTROM))
        .unwrap();
    let d = branch
        .add_atom(atom("N"), Quantity::new(Point3::origin(), ANGSTROM))
        .unwrap();
    assert_ne!(c, d);
    assert!(branch.atom(c).is_err());
    first.clear();
    let e = first
        .add_atom(atom("C"), Quantity::new(Point3::origin(), ANGSTROM))
        .unwrap();
    assert_ne!(a, e);
    assert!(first.atom(a).is_err());
}
