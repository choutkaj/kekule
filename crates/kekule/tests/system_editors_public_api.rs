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

fn has_bond(model: &Model, a: usize, b: usize) -> bool {
    let a = model.atom_ids()[a];
    let b = model.atom_ids()[b];
    model.bonds().any(|(id, bond)| {
        let endpoints = (
            InstanceAtomId::new(id.molecule(), bond.a()),
            InstanceAtomId::new(id.molecule(), bond.b()),
        );
        endpoints == (a, b) || endpoints == (b, a)
    })
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
    let split = editor.clone().finish().unwrap();
    assert_eq!(split.topology().instance_count(), 2);
    assert_eq!(split.hierarchy().chains().count(), 1);
    assert_eq!(split.hierarchy().residues().count(), 1);
    assert_eq!(split.hierarchy().atom_sites().count(), 3);
    assert_eq!(split.positions(), source.positions());
    for (site, before) in source.hierarchy().atom_sites() {
        let after = split.hierarchy().atom_site(site).unwrap();
        assert_eq!(
            source.position(before.atom()).unwrap(),
            split.position(after.atom()).unwrap()
        );
    }
    editor
        .add_bond(handles[0], handles[2], BondOrder::Single)
        .unwrap();
    let merged = editor.finish().unwrap();
    assert_eq!(merged.topology().instance_count(), 1);
    assert!(!has_bond(&merged, 0, 1));
    assert_eq!(merged.topology().bond_count(), 2);
    assert_eq!(merged.positions(), source.positions());
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
    editor
        .add_bond(handles[0], handles[1], BondOrder::Single)
        .unwrap();
    let result = editor.finish().unwrap();
    assert_eq!(result.topology().instance_count(), 1);
    assert_eq!(result.topology().bond_count(), 1);
    assert_eq!(result.hierarchy().chains().count(), 2);
    assert_eq!(result.hierarchy().residues().count(), 2);
    assert!(!result.topology().molecule_instance_properties().has_data());
    assert_eq!(result.positions(), source.positions());
    assert!(has_bond(&result, 0, 1));
    for (site, before) in source.hierarchy().atom_sites() {
        let after = result.hierarchy().atom_site(site).unwrap();
        assert_eq!(
            source.position(before.atom()).unwrap(),
            result.position(after.atom()).unwrap()
        );
        assert_eq!(after.atom().molecule(), result.atom_ids()[0].molecule());
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
    let result = editor.finish().unwrap();
    assert_eq!(result.definition_count(), 2);
    let untouched = [result.atom_ids()[0], result.atom_ids()[2]];
    assert_eq!(
        result
            .instance(untouched[0].molecule())
            .unwrap()
            .definition(),
        result
            .instance(untouched[1].molecule())
            .unwrap()
            .definition()
    );
    for atom in untouched {
        assert_eq!(result.atom(atom).unwrap().element.symbol(), "O");
        assert_eq!(
            result
                .definition_for_instance(atom.molecule())
                .unwrap()
                .molecule()
                .perception(),
            water.perception()
        );
    }
    let changed = result.atom_ids()[1];
    assert_eq!(result.atom(changed).unwrap().element.symbol(), "N");
    assert_ne!(
        result
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
    let result = editor.finish().unwrap();
    assert!(Arc::ptr_eq(&result.shared_topology(), &topology));
    assert_eq!(
        result.properties().get(&key("label")),
        source.properties().get(&key("label"))
    );
    assert_eq!(
        result.position(source.atom_ids()[0]).unwrap().value().x,
        3.0
    );
    assert_ne!(result.positions(), source.positions());
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
    editor
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
    let result = editor.finish().unwrap();
    assert_eq!(result.hierarchy().chains().count(), 1);
    assert_eq!(result.hierarchy().residues().count(), 1);
    assert_eq!(result.hierarchy().atom_sites().count(), 1);
    assert!(result.properties().owner_is_empty());
    let retained = result.atom_ids()[0];
    assert_eq!(
        result.position(retained).unwrap(),
        source.position(source_atoms[1]).unwrap()
    );
    assert_eq!(result.occupancy(retained).unwrap(), Some(0.7));
    assert_eq!(
        result
            .topology()
            .atom_property(retained, &key("static"))
            .unwrap(),
        Some(PropertyValue::Int(2))
    );
    let new = result.atom_ids()[1];
    assert_eq!(result.occupancy(new).unwrap(), None);
    assert_eq!(
        result.atom_property(new, &key("live")).unwrap(),
        Some(PropertyValue::Int(8))
    );
    assert_eq!(result.atom_count(), 2);
    assert_eq!(result.atom(new).unwrap().element.symbol(), "O");
    assert_eq!(result.positions(), &positions(&[9.0, 5.0]));
    assert_eq!(
        result
            .topology()
            .atom_property(new, &key("static"))
            .unwrap(),
        None
    );
}

#[test]
fn finish_and_try_finish_publish_the_same_split_model_and_topology() {
    let source = model("CCC");
    let mut editor = source.edit();
    let bonds = editor.bond_ids().collect::<Vec<_>>();
    editor
        .set_bond_property(bonds[1], key("dynamic"), Some(PropertyValue::Int(17)))
        .unwrap();
    editor
        .set_topology_bond_property(bonds[1], key("static"), Some(PropertyValue::Int(23)))
        .unwrap();
    editor.delete_bond(bonds[0]).unwrap();
    let structural = editor.topology_editor().clone();
    let finished_topology: Arc<Topology> = structural.clone().finish().unwrap();
    let recovered_topology: Arc<Topology> = structural.try_finish().unwrap();
    assert!(finished_topology.same_layout(&recovered_topology));
    assert_eq!(
        finished_topology.properties(),
        recovered_topology.properties()
    );
    assert_eq!(finished_topology.instance_count(), 2);
    assert_eq!(finished_topology.bond_count(), 1);
    let finished: Model = editor.clone().finish().unwrap();
    let recovered: Model = editor.try_finish().unwrap();
    // Independent publications have different snapshot identities. Compare the
    // represented layout and stored values explicitly.
    assert!(finished.topology().same_layout(recovered.topology()));
    assert_eq!(
        finished.topology().properties(),
        recovered.topology().properties()
    );
    assert_eq!(finished.properties(), recovered.properties());
    assert_eq!(finished.positions(), recovered.positions());
    assert_eq!(finished.cell(), recovered.cell());
    assert!(finished.topology().same_layout(&finished_topology));
    assert_eq!(
        finished.topology().properties(),
        finished_topology.properties()
    );
    assert_eq!(finished.positions(), source.positions());
    let bond = finished.bond_ids()[0];
    assert_eq!(
        finished.bond_property(bond, &key("dynamic")).unwrap(),
        Some(PropertyValue::Int(17))
    );
    assert_eq!(
        finished
            .topology()
            .bond_property(bond, &key("static"))
            .unwrap(),
        Some(PropertyValue::Int(23))
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
    editor
        .add_molecule(&self::molecule("O"), &positions(&[8.0]))
        .unwrap();
    let result = editor.finish().unwrap();
    assert_eq!(
        &result.topology().atom_ids()[..source.atom_count()],
        source.topology().atom_ids()
    );
    assert_eq!(
        &result.topology().bond_ids()[..source.topology().bond_count()],
        source.topology().bond_ids()
    );
    for &id in source.topology().atom_ids() {
        assert_eq!(result.position(id).unwrap(), source.position(id).unwrap());
    }
    assert_eq!(result.atom_count(), 3);
    assert_eq!(
        result.atom(result.atom_ids()[2]).unwrap().element.symbol(),
        "O"
    );
    assert_eq!(result.positions(), &positions(&[3.0, 7.0, 8.0]));
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
    let classified = editor.clone().finish().unwrap();
    let class = |result: &Topology, row: usize| {
        result
            .definition_for_instance(result.atom_ids()[row].molecule())
            .unwrap()
            .class()
    };
    assert_eq!(class(&classified, 0), MoleculeClass::Other);
    assert_eq!(class(&classified, 2), MoleculeClass::SmallMolecule);
    let bond = editor.bond_ids().next().unwrap();
    editor.delete_bond(bond).unwrap();
    editor
        .set_molecule_class(atoms[0], MoleculeClass::Ion)
        .unwrap();
    editor
        .set_molecule_class(atoms[1], MoleculeClass::Other)
        .unwrap();
    let split = editor.clone().finish().unwrap();
    assert_eq!(class(&split, 0), MoleculeClass::Ion);
    assert_eq!(class(&split, 1), MoleculeClass::Other);
    assert_eq!(class(&split, 2), MoleculeClass::SmallMolecule);
    editor
        .add_bond(atoms[0], atoms[1], BondOrder::Single)
        .unwrap();
    let restored = editor.finish().unwrap();
    assert_eq!(class(&restored, 0), MoleculeClass::SmallMolecule);
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
    let result = editor.finish().unwrap();
    let target = result.bond_ids()[0];
    assert_eq!(result.topology().instance_count(), 2);
    let actual = result.bond(target).unwrap();
    let a = result.atom_ids()[0];
    let b = result.atom_ids()[1];
    assert_eq!(actual.endpoints(), (a.atom(), b.atom()));
    assert_eq!(
        result.bond_property(target, &key("dynamic")).unwrap(),
        Some(PropertyValue::Int(10))
    );
    assert_eq!(
        result
            .topology()
            .bond_property(target, &key("static"))
            .unwrap(),
        Some(PropertyValue::Int(20))
    );
    assert_eq!(
        result
            .topology()
            .definition_for_instance(target.molecule())
            .unwrap()
            .molecule()
            .bond_property(target.bond(), &key("definition"))
            .unwrap(),
        Some(PropertyValue::Int(30))
    );
    assert_eq!(result.bond_ids().len(), 1);
    assert_eq!(actual.order, BondOrder::Single);
    // Rewiring CC.O to C-O + C publishes the connected pair before the isolated C.
    assert_eq!(result.positions(), &positions(&[0.0, 2.0, 1.0]));
    assert_eq!(
        result
            .atoms()
            .map(|(_, a)| a.element.symbol())
            .collect::<Vec<_>>(),
        ["C", "O", "C"]
    );
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
    let failed = editor.try_finish().unwrap_err();
    assert!(failed.source().is_some());
    let mut editor = failed.into_editor();
    let replacement = editor
        .add_atom(atom("N"), Quantity::new(Point3::origin(), ANGSTROM))
        .unwrap();
    assert_ne!(handle, replacement);
    assert!(editor.atom(handle).is_err());
    let result = editor.finish().unwrap();
    assert_eq!(result.atom_count(), 1);
    assert_eq!(result.topology().instance_count(), 1);
    assert_eq!(
        result.atom(result.atom_ids()[0]).unwrap().element.symbol(),
        "N"
    );
    assert_eq!(result.positions(), &positions(&[0.0]));
    assert_eq!(result.occupancy(result.atom_ids()[0]).unwrap(), None);
    assert_eq!(
        result
            .atom_property(result.atom_ids()[0], &key("tag"))
            .unwrap(),
        None
    );

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
    let result = editor.finish().unwrap();
    assert_eq!(
        result
            .atom_property(result.atom_ids()[0], &key("tag"))
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
