use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor};
use kekule::geometry::Point3;
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{Model, ModelEditor, Positions};
use kekule::topology::{
    AtomSiteMetadata, InstanceAtomId, MoleculeClass, Topology, TopologyBuilder, TopologyEditor,
};
use kekule::units::{Quantity, ANGSTROM, BOHR, NANOMETER};
use std::sync::Arc;

fn molecule(text: &str) -> Molecule {
    kekule::smiles::to_molecules(text).unwrap().pop().unwrap()
}
fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}
/// An atom whose hydrogen count is fixed at zero, so publication needs no
/// perception.
fn bare_atom(atomic_number: u8) -> Atom {
    let mut atom = Atom::new(Element::from_atomic_number(atomic_number).unwrap());
    atom.hydrogens = kekule::core::ImplicitHydrogens::Fixed(0);
    atom
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
                molecule.clone(),
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
    let a = model.topology().atom_ids()[a];
    let b = model.topology().atom_ids()[b];
    model.topology().bonds().any(|bond| {
        let [left, right] = bond.atoms().map(|atom| atom.id());
        (left, right) == (a, b) || (left, right) == (b, a)
    })
}

#[test]
fn split_then_merge_preserves_handles_geometry_and_hierarchy() {
    let molecule = molecule("CCC");
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(molecule.clone()).unwrap();
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
    assert_eq!(split.topology().hierarchy().chains().count(), 1);
    assert_eq!(split.topology().hierarchy().residues().count(), 1);
    assert_eq!(split.topology().hierarchy().atom_sites().count(), 3);
    assert_eq!(split.positions(), source.positions());
    for (site, before) in source.topology().hierarchy().atom_sites() {
        let after = split.topology().hierarchy().atom_site(site).unwrap();
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
fn merging_fragments_of_two_split_occurrences_keeps_their_annotations() {
    // Two H2 occurrences with per-atom and per-bond definition annotations.
    let tagged = |first: i64| {
        let mut h2 = molecule("[H][H]");
        let atoms = h2.atom_ids().collect::<Vec<_>>();
        for (offset, &atom) in atoms.iter().enumerate() {
            h2.properties_mut()
                .atoms_mut()
                .set_value(
                    key("tag"),
                    atom,
                    Some(PropertyValue::Int(first + offset as i64)),
                )
                .unwrap();
        }
        let bond = h2.bond_ids().next().unwrap();
        h2.properties_mut()
            .bonds_mut()
            .set_value(key("bond_tag"), bond, Some(PropertyValue::Int(first)))
            .unwrap();
        h2
    };
    let source = Arc::new(Topology::from_molecules([tagged(1), tagged(3)]).unwrap());
    let mut editor = source.edit();
    let handles = source
        .atom_ids()
        .iter()
        .map(|&id| editor.atom_handle(id).unwrap())
        .collect::<Vec<_>>();
    // Split both occurrences, then join one fragment of each (H-H + H-H -> H + H-H + H).
    for &bond in source.bond_ids() {
        editor
            .delete_bond(editor.bond_handle(bond).unwrap())
            .unwrap();
    }
    editor
        .add_bond(handles[1], handles[2], BondOrder::Single)
        .unwrap();
    let edited = editor.finish().unwrap();

    let mut occurrences = edited
        .molecules()
        .map(|instance| {
            let molecule = instance.molecule();
            let mut tags = molecule
                .atom_ids()
                .map(|atom| {
                    molecule
                        .properties()
                        .atoms()
                        .value(&key("tag"), atom)
                        .unwrap()
                })
                .collect::<Vec<_>>();
            tags.sort_by_key(|tag| format!("{tag:?}"));
            let bond_tags = molecule
                .bond_ids()
                .map(|bond| {
                    molecule
                        .properties()
                        .bonds()
                        .value(&key("bond_tag"), bond)
                        .unwrap()
                })
                .collect::<Vec<_>>();
            (tags, bond_tags)
        })
        .collect::<Vec<_>>();
    occurrences.sort_by_key(|o| format!("{o:?}"));
    let int = |value| Some(PropertyValue::Int(value));
    assert_eq!(
        occurrences,
        vec![
            (vec![int(1)], vec![]),
            (vec![int(2), int(3)], vec![None]),
            (vec![int(4)], vec![]),
        ]
    );
    assert_eq!(edited.bond_count(), 1);
}

#[test]
fn cross_instance_bond_merges_but_preserves_distinct_chains() {
    let molecule = molecule("C");
    let mut builder = TopologyBuilder::new();
    for label in ["A", "B"] {
        let instance = builder.add_molecule(molecule.clone()).unwrap();
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
        .properties_mut()
        .molecule_instances_mut()
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
    assert_eq!(result.topology().hierarchy().chains().count(), 2);
    assert_eq!(result.topology().hierarchy().residues().count(), 2);
    assert!(!result
        .topology()
        .properties()
        .molecule_instances()
        .has_data());
    assert_eq!(result.positions(), source.positions());
    assert!(has_bond(&result, 0, 1));
    for (site, before) in source.topology().hierarchy().atom_sites() {
        let after = result.topology().hierarchy().atom_site(site).unwrap();
        assert_eq!(
            source.position(before.atom()).unwrap(),
            result.position(after.atom()).unwrap()
        );
        assert_eq!(
            after.atom().molecule(),
            result.topology().atom_ids()[0].molecule()
        );
    }
}

#[test]
fn editing_one_reused_occurrence_preserves_others_and_their_perception() {
    let mut water = molecule("O");
    water.perceive().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(water.clone()).unwrap();
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
            .molecule(untouched[0].molecule())
            .unwrap()
            .definition_id(),
        result
            .molecule(untouched[1].molecule())
            .unwrap()
            .definition_id()
    );
    for atom in untouched {
        assert_eq!(result.atom(atom).unwrap().element.symbol(), "O");
        assert_eq!(
            result
                .molecule(atom.molecule())
                .unwrap()
                .molecule()
                .perception(),
            water.perception()
        );
    }
    let changed = result.atom_ids()[1];
    assert_eq!(result.atom(changed).unwrap().element.symbol(), "N");
    assert_ne!(
        result.molecule(changed.molecule()).unwrap().class(),
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
        .conformation_mut()
        .properties_mut()
        .owner_mut()
        .insert(key("label"), PropertyValue::Int(9))
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
        result.properties().owner().get(&key("label")),
        source.properties().owner().get(&key("label"))
    );
    assert_eq!(
        result
            .position(source.topology().atom_ids()[0])
            .unwrap()
            .value()
            .x,
        3.0
    );
    assert_ne!(result.positions(), source.positions());
}

#[test]
fn deletion_prunes_hierarchy_and_projects_properties_and_coordinates() {
    let mut builder = TopologyBuilder::new();
    let molecule = molecule("C");
    for label in ["A", "B"] {
        let instance = builder.add_molecule(molecule.clone()).unwrap();
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
        .properties_mut()
        .atoms_mut()
        .insert(key("static"), PropertyColumn::Int(vec![Some(1), Some(2)]))
        .unwrap();
    let mut source = Model::new(builder.build().unwrap(), positions(&[4.0, 9.0])).unwrap();
    let source_atoms = source.topology().atom_ids().to_vec();
    let row = source.topology().atom_index(source_atoms[1]).unwrap();
    source
        .conformation_mut()
        .set_occupancy(row, Some(0.7))
        .unwrap();
    source
        .conformation_mut()
        .properties_mut()
        .owner_mut()
        .insert(key("energy"), PropertyValue::Int(10))
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
    let column = editor.atom_property_column(&key("live")).unwrap();
    editor
        .insert_atom_property_column(key("live"), column)
        .unwrap();
    let result = editor.finish().unwrap();
    assert_eq!(result.topology().hierarchy().chains().count(), 1);
    assert_eq!(result.topology().hierarchy().residues().count(), 1);
    assert_eq!(result.topology().hierarchy().atom_sites().count(), 1);
    assert!(result.properties().owner().is_empty());
    let retained = result.topology().atom_ids()[0];
    assert_eq!(
        result.position(retained).unwrap(),
        source.position(source_atoms[1]).unwrap()
    );
    assert_eq!(
        result
            .occupancy(result.topology().atom_index(retained).unwrap())
            .unwrap(),
        Some(0.7)
    );
    assert_eq!(
        result
            .topology()
            .atom(retained)
            .unwrap()
            .property(&key("static")),
        Some(PropertyValue::Int(2))
    );
    let new = result.topology().atom_ids()[1];
    assert_eq!(
        result
            .occupancy(result.topology().atom_index(new).unwrap())
            .unwrap(),
        None
    );
    assert_eq!(
        result.atom(new).unwrap().realization_property(&key("live")),
        Some(PropertyValue::Int(8))
    );
    assert_eq!(result.atom_count(), 2);
    assert_eq!(result.topology().atom(new).unwrap().element.symbol(), "O");
    assert_eq!(result.positions(), &positions(&[9.0, 5.0]));
    assert_eq!(
        result
            .topology()
            .atom(new)
            .unwrap()
            .property(&key("static")),
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
    let structural = editor.topology().clone();
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
    let bond = finished.topology().bond_ids()[0];
    assert_eq!(
        finished
            .properties()
            .bonds()
            .value(
                &key("dynamic"),
                finished.topology().bond_index(bond).unwrap()
            )
            .unwrap(),
        Some(PropertyValue::Int(17))
    );
    assert_eq!(
        finished
            .topology()
            .bond(bond)
            .unwrap()
            .property(&key("static")),
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
        Model::from_molecule(molecule.finish().unwrap().clone(), &positions(&[3.0, 7.0])).unwrap();
    let mut editor = source.edit();
    editor
        .add_molecule((self::molecule("O")).clone(), &positions(&[8.0]))
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
        result
            .topology()
            .atom(result.topology().atom_ids()[2])
            .unwrap()
            .element
            .symbol(),
        "O"
    );
    assert_eq!(result.positions(), &positions(&[3.0, 7.0, 8.0]));
}

#[test]
fn merge_property_conflicts_roll_back_and_stereo_survives_unrelated_append() {
    let mut left = molecule("C");
    let id = left.atom_ids().next().unwrap();
    left.properties_mut()
        .atoms_mut()
        .set_value(key("tag"), id, Some(PropertyValue::Int(3)))
        .unwrap();
    let mut right = molecule("C");
    let id = right.atom_ids().next().unwrap();
    right
        .properties_mut()
        .atoms_mut()
        .set_value(key("tag"), id, Some(PropertyValue::String("x".into())))
        .unwrap();
    let mut editor = TopologyEditor::new();
    let left = *editor
        .add_molecule(left.clone())
        .unwrap()
        .atoms()
        .values()
        .next()
        .unwrap();
    let right = *editor
        .add_molecule(right.clone())
        .unwrap()
        .atoms()
        .values()
        .next()
        .unwrap();
    let before = format!("{editor:?}");
    assert!(editor.add_bond(left, right, BondOrder::Single).is_err());
    assert_eq!(format!("{editor:?}"), before);
    let stereo = molecule("F[C@](Cl)(Br)I");
    let mut editor = Topology::from_molecule(stereo.clone())
        .unwrap()
        .into_editor();
    editor.add_molecule((molecule("O")).clone()).unwrap();
    let topology = editor.finish().unwrap();
    assert_eq!(
        topology.definitions().next().unwrap().molecule().graph(),
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
    // Draft rows are atom IDs; the deleted atom keeps an allocated, missing row.
    assert_eq!(
        editor.properties_mut().atoms_mut().insert(
            key("tag"),
            PropertyColumn::Int(vec![Some(0), Some(1), Some(2)])
        ),
        Err(kekule::properties::PropertyError::RemovedRow {
            index: deleted.index()
        })
    );
    editor
        .properties_mut()
        .atoms_mut()
        .insert(
            key("tag"),
            PropertyColumn::Int(vec![None, Some(1), Some(2)]),
        )
        .unwrap();
    let column = editor
        .properties()
        .atoms()
        .get(&key("tag"))
        .unwrap()
        .clone();
    assert_eq!(column.len(), 3);
    assert_eq!(
        editor
            .properties_mut()
            .atoms_mut()
            .insert(key("tag"), column.clone())
            .unwrap(),
        Some(column.clone())
    );
    assert_eq!(
        editor.properties_mut().atoms_mut().remove(&key("tag")),
        Some(column)
    );
    assert!(MoleculeEditor::new().finish().is_err());
}

#[test]
fn classification_overrides_are_occurrence_and_component_local() {
    let carbon = molecule("CC");
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(carbon.clone()).unwrap();
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
            .molecule(result.atom_ids()[row].molecule())
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
    let source_bond = source.topology().bond_ids()[0];
    let row = source.topology().bond_index(source_bond).unwrap();
    source
        .conformation_mut()
        .properties_mut()
        .bonds_mut()
        .set_value(key("dynamic"), row, Some(PropertyValue::Int(10)))
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
    let target = result.topology().bond_ids()[0];
    assert_eq!(result.topology().instance_count(), 2);
    let actual = result.topology().bond(target).unwrap();
    // Survivors keep source dense order: C0, the now-isolated C1, then O.
    let [a, isolated, b] = [0, 1, 2].map(|index| result.topology().atom_ids()[index]);
    assert_eq!(a.molecule(), b.molecule());
    assert_ne!(a.molecule(), isolated.molecule());
    assert_eq!(actual.endpoints(), (a.atom(), b.atom()));
    assert_eq!(
        result
            .properties()
            .bonds()
            .value(
                &key("dynamic"),
                result.topology().bond_index(target).unwrap()
            )
            .unwrap(),
        Some(PropertyValue::Int(10))
    );
    assert_eq!(
        result
            .topology()
            .bond(target)
            .unwrap()
            .property(&key("static")),
        Some(PropertyValue::Int(20))
    );
    assert_eq!(
        result
            .topology()
            .molecule(target.molecule())
            .unwrap()
            .molecule()
            .properties()
            .bonds()
            .value(&key("definition"), target.bond())
            .unwrap(),
        Some(PropertyValue::Int(30))
    );
    assert_eq!(result.topology().bond_ids().len(), 1);
    assert_eq!(actual.order, BondOrder::Single);
    // Rewiring CC.O to C-O + C repartitions instances without reordering atoms.
    assert_eq!(result.positions(), &positions(&[0.0, 1.0, 2.0]));
    assert_eq!(
        result
            .topology()
            .atoms()
            .map(|a| a.element.symbol())
            .collect::<Vec<_>>(),
        ["C", "C", "O"]
    );
}

#[test]
fn bond_replacement_rolls_back_failed_group_merges_and_preserves_the_handle() {
    use kekule::topology::EditBond;

    let source = model("CC.O");
    let mut editor = TopologyEditor::from_topology(source.shared_topology());
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    let bond = editor.bond_ids().next().unwrap();
    editor
        .set_bond_property(bond, key("tag"), Some(PropertyValue::Int(7)))
        .unwrap();
    let previous = editor.bond(bond).unwrap();
    let snapshot = format!("{editor:?}");
    // The oxygen group is merged before the molecular editor rejects the loop.
    assert!(editor
        .replace_bond(bond, EditBond::new(atoms[2], atoms[2], BondOrder::Double))
        .is_err());
    assert_eq!(format!("{editor:?}"), snapshot);

    let replacement = EditBond::new(atoms[0], atoms[2], BondOrder::Double);
    assert_eq!(editor.replace_bond(bond, replacement).unwrap(), previous);
    assert_eq!(editor.bond(bond).unwrap(), replacement);
    assert_eq!(
        editor.bond_property(bond, &key("tag")).unwrap(),
        Some(PropertyValue::Int(7))
    );
    let result = editor.finish().unwrap();
    assert_eq!(result.instance_count(), 2);
    assert_eq!(result.bond_count(), 1);
    assert_eq!(result.bonds().next().unwrap().order, BondOrder::Double);
}

#[test]
fn model_property_batch_recreation_preserves_unit_symbols_and_signed_zero() {
    use kekule::units::Unit;

    let mut editor = model("C").into_editor();
    let atom = editor.atom_ids().next().unwrap();
    editor
        .set_atom_property(
            atom,
            key("length"),
            Some(PropertyValue::Real {
                value: 0.0,
                unit: NANOMETER,
            }),
        )
        .unwrap();
    let alias = Unit::new(NANOMETER.dimension(), NANOMETER.scale(), Some("custom_nm")).unwrap();
    assert_eq!(alias, NANOMETER); // Unit equality deliberately ignores symbols.
    editor
        .set_atom_properties(
            key("length"),
            [
                (atom, None),
                (
                    atom,
                    Some(PropertyValue::Real {
                        value: -0.0,
                        unit: alias,
                    }),
                ),
            ],
        )
        .unwrap();
    let Some(PropertyValue::Real { value, unit }) =
        editor.atom_property(atom, &key("length")).unwrap()
    else {
        panic!("expected recreated real column");
    };
    assert!(value.is_sign_negative());
    assert_eq!(unit.symbol(), Some("custom_nm"));
    editor.finish().unwrap();
}

#[test]
fn editor_property_batches_are_transactional_and_follow_ordered_update_semantics() {
    macro_rules! check_batches {
        ($editor:expr) => {{
            let editor = &mut $editor;
            let atoms = editor.atom_ids().collect::<Vec<_>>();
            let bonds = editor.bond_ids().collect::<Vec<_>>();
            editor
                .set_atom_property(atoms[0], key("untouched"), Some(PropertyValue::Int(1)))
                .unwrap();
            editor
                .set_bond_property(bonds[0], key("untouched"), Some(PropertyValue::Int(2)))
                .unwrap();
            let untouched_atoms = editor.atom_property_column(&key("untouched"));
            let untouched_bonds = editor.bond_property_column(&key("untouched"));
            editor
                .set_atom_property(atoms[0], key("edited"), Some(PropertyValue::Int(3)))
                .unwrap();
            editor
                .set_bond_property(bonds[0], key("edited"), Some(PropertyValue::Int(3)))
                .unwrap();
            editor
                .set_atom_properties(
                    key("edited"),
                    [
                        (atoms[0], Some(PropertyValue::Int(3))),
                        (atoms[0], None),
                        (
                            atoms[1],
                            Some(PropertyValue::Real {
                                value: 10.0,
                                unit: ANGSTROM,
                            }),
                        ),
                        (
                            atoms[1],
                            Some(PropertyValue::Real {
                                value: 2.0,
                                unit: NANOMETER,
                            }),
                        ),
                    ],
                )
                .unwrap();
            assert_eq!(
                editor.atom_property(atoms[1], &key("edited")).unwrap(),
                Some(PropertyValue::Real {
                    value: 20.0,
                    unit: ANGSTROM
                })
            );
            let snapshot = editor.atom_property_column(&key("edited"));
            assert!(editor
                .set_atom_properties(
                    key("edited"),
                    [
                        (
                            atoms[0],
                            Some(PropertyValue::Real {
                                value: 1.0,
                                unit: ANGSTROM
                            })
                        ),
                        (atoms[1], Some(PropertyValue::String("bad type".into()))),
                    ]
                )
                .is_err());
            assert_eq!(editor.atom_property_column(&key("edited")), snapshot);
            editor
                .set_bond_properties(
                    key("edited"),
                    [
                        (bonds[0], Some(PropertyValue::Int(3))),
                        (bonds[0], None),
                        (bonds[1], Some(PropertyValue::String("new type".into()))),
                    ],
                )
                .unwrap();
            assert_eq!(
                editor.bond_property(bonds[1], &key("edited")).unwrap(),
                Some(PropertyValue::String("new type".into()))
            );
            let snapshot = editor.bond_property_column(&key("edited"));
            assert!(editor
                .set_bond_properties(
                    key("edited"),
                    [
                        (bonds[0], Some(PropertyValue::String("staged".into()))),
                        (bonds[1], Some(PropertyValue::Int(4))),
                    ]
                )
                .is_err());
            assert_eq!(editor.bond_property_column(&key("edited")), snapshot);
            editor
                .set_bond_properties(key("edited"), [(bonds[1], None)])
                .unwrap();
            assert!(editor.bond_property_column(&key("edited")).is_none());
            assert_eq!(
                editor.atom_property_column(&key("untouched")),
                untouched_atoms
            );
            assert_eq!(
                editor.bond_property_column(&key("untouched")),
                untouched_bonds
            );
        }};
    }

    let source = model("CCC");
    let mut topology_editor = TopologyEditor::from_topology(source.shared_topology());
    topology_editor
        .insert_property(key("owner"), PropertyValue::String("retained".into()))
        .unwrap();
    check_batches!(topology_editor);
    assert_eq!(
        topology_editor.owner_properties().get(&key("owner")),
        Some(&PropertyValue::String("retained".into()))
    );
    topology_editor.finish().unwrap();
    let mut editor = source.into_editor();
    editor
        .owner_properties_mut()
        .insert(key("owner"), PropertyValue::String("retained".into()))
        .unwrap();
    check_batches!(editor);
    assert_eq!(
        editor.owner_properties().get(&key("owner")),
        Some(&PropertyValue::String("retained".into()))
    );
    // Generic realization keys are unreserved: "occupancy" is independent of
    // the typed occupancy state.
    let atom = editor.atom_ids().next().unwrap();
    editor
        .set_atom_properties(key("occupancy"), [(atom, Some(PropertyValue::Int(1)))])
        .unwrap();
    assert_eq!(editor.occupancy(atom).unwrap(), None);
    editor.finish().unwrap();
}

#[test]
fn recovery_keeps_builder_state_and_editor_handles() {
    use std::error::Error;
    let carbon = molecule("C");
    let mut builder = Model::builder();
    let definition = builder.add_molecule_definition(carbon.clone()).unwrap();
    let failed = builder.try_build().unwrap_err();
    assert!(failed.source().is_some());
    assert_eq!(failed.builder().topology_builder().definition_count(), 1);
    let mut builder = failed.into_builder();
    let instance = builder
        .add_instance(definition, &positions(&[3.0]))
        .unwrap();
    let id = InstanceAtomId::new(instance, carbon.atom_ids().next().unwrap());
    let row = builder.atom_index(id).unwrap();
    builder
        .conformation_mut()
        .set_occupancy(row, Some(0.5))
        .unwrap();
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
    assert_eq!(
        builder
            .conformation()
            .occupancy(builder.atom_index(id).unwrap())
            .unwrap(),
        Some(0.5)
    );
    let row = builder.atom_index(id).unwrap();
    builder
        .conformation_mut()
        .properties_mut()
        .atoms_mut()
        .set_values(key("tag"), [(row, Some(PropertyValue::Int(7)))])
        .unwrap();
    assert_eq!(
        builder
            .conformation()
            .properties()
            .atoms()
            .value(&key("tag"), builder.atom_index(id).unwrap())
            .unwrap(),
        Some(PropertyValue::Int(7))
    );
    assert_eq!(builder.atom_ids(), &[id]);
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
        result
            .topology()
            .atom(result.topology().atom_ids()[0])
            .unwrap()
            .element
            .symbol(),
        "N"
    );
    assert_eq!(result.positions(), &positions(&[0.0]));
    assert_eq!(
        result
            .occupancy(
                result
                    .topology()
                    .atom_index(result.topology().atom_ids()[0])
                    .unwrap()
            )
            .unwrap(),
        None
    );
    assert_eq!(
        result
            .atom(result.topology().atom_ids()[0])
            .unwrap()
            .realization_property(&key("tag")),
        None
    );

    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(carbon.clone()).unwrap();
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
    let definition = builder.add_molecule_definition(carbon.clone()).unwrap();
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
    let row = builder.atom_index(a).unwrap();
    builder
        .conformation_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(key("dynamic"), row, Some(PropertyValue::Int(4)))
        .unwrap();
    builder
        .topology_builder_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(
            key("static"),
            kekule::topology::TopologyAtomIndex::new(0),
            Some(PropertyValue::Int(8)),
        )
        .unwrap();
    builder
        .conformation_mut()
        .properties_mut()
        .owner_mut()
        .insert(key("old_owner"), PropertyValue::Int(1))
        .unwrap();
    builder
        .topology_builder_mut()
        .properties_mut()
        .owner_mut()
        .insert(key("old_owner"), PropertyValue::Int(2))
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
    assert_eq!(
        &result.topology().atom_ids()[..4],
        source.topology().atom_ids()
    );
    assert_eq!(
        &result.topology().bond_ids()[..2],
        source.topology().bond_ids()
    );
    assert_eq!(
        result.positions().values().value()[..4],
        source.positions().values().value()[..]
    );
    assert_eq!(result.topology().definition_count(), 1);
    assert_eq!(result.topology().hierarchy(), source.topology().hierarchy());
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
        result
            .atom(a)
            .unwrap()
            .realization_property(&key("dynamic")),
        Some(PropertyValue::Int(4))
    );
    assert_eq!(
        result.topology().atom(a).unwrap().property(&key("static")),
        Some(PropertyValue::Int(8))
    );
    assert_eq!(
        result
            .atom(result.topology().atom_ids()[4])
            .unwrap()
            .realization_property(&key("dynamic")),
        None
    );
    assert!(result.properties().owner().get(&key("old_owner")).is_none());
    assert!(result
        .topology()
        .properties()
        .owner()
        .get(&key("old_owner"))
        .is_none());
}

#[test]
fn sparse_batches_are_atomic_and_deleted_rows_do_not_poison_new_values() {
    let mut source = model("CCC");
    let ids = source.topology().atom_ids().to_vec();
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
    let rows = ids
        .iter()
        .map(|id| source.topology().atom_index(*id).unwrap())
        .collect::<Vec<_>>();
    let error = source
        .conformation_mut()
        .properties_mut()
        .atoms_mut()
        .set_values(
            key("tag"),
            [
                (rows[0], Some(PropertyValue::Int(1))),
                (rows[1], Some(PropertyValue::String("bad".into()))),
            ],
        )
        .unwrap_err();
    assert_eq!(
        error,
        kekule::properties::PropertyError::TypeMismatch { key: key("tag") }
    );
    assert!(source.properties().atoms().get(&key("tag")).is_none());
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
            .atom(result.topology().atom_ids()[0])
            .unwrap()
            .realization_property(&key("tag")),
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
    let source = Arc::new(Topology::from_molecule(carbon.clone()).unwrap());
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
            .map(|d| d.molecule().stereo_elements().count())
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
            .map(|d| d.molecule().stereo_elements().count())
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

#[test]
fn publication_correspondence_tracks_reordering_deletion_and_bonds() {
    let mut editor = ModelEditor::new();
    let a = editor
        .add_atom(bare_atom(6), Quantity::new(Point3::new(1., 0., 0.), BOHR))
        .unwrap();
    let b = editor
        .add_atom(bare_atom(8), Quantity::new(Point3::new(2., 0., 0.), BOHR))
        .unwrap();
    let c = editor
        .add_atom(bare_atom(7), Quantity::new(Point3::new(3., 0., 0.), BOHR))
        .unwrap();
    let removed = editor
        .add_atom(bare_atom(1), Quantity::new(Point3::new(4., 0., 0.), BOHR))
        .unwrap();
    editor.delete_atom(removed).unwrap();
    let bond = editor.add_bond(a, c, BondOrder::Single).unwrap();
    let (model, map) = editor.finish_with_correspondence().unwrap();
    assert!(map.atom(removed).is_none());
    for (handle, x) in [(a, 1.), (b, 2.), (c, 3.)] {
        let (id, index) = map.atom(handle).unwrap();
        assert_eq!(model.topology().atom_ids()[index.index()], id);
        assert!((model.position(id).unwrap().value_in(BOHR).unwrap().x - x).abs() < 1e-12);
    }
    let (id, index) = map.bond(bond).unwrap();
    assert_eq!(model.topology().bond_ids()[index.index()], id);
    let mut next = model.edit();
    let handle = next.bond_handle(id).unwrap();
    next.delete_bond(handle).unwrap();
    let (split, second) = next.finish_with_correspondence().unwrap();
    assert!(second.bond(handle).is_none());
    assert_eq!(split.topology().instance_count(), 3);
}

#[test]
fn no_op_publication_retains_snapshot_and_foreign_handles_do_not_map() {
    let mut e = ModelEditor::new();
    e.add_atom(bare_atom(2), Quantity::new(Point3::origin(), BOHR))
        .unwrap();
    let model = e.finish().unwrap();
    let source = model.shared_topology();
    let edit = TopologyEditor::from_topology(source.clone());
    let handle = edit.atom_handle(model.topology().atom_ids()[0]).unwrap();
    let (topology, map) = edit.finish_with_correspondence().unwrap();
    assert!(Arc::ptr_eq(&source, &topology));
    assert_eq!(map.atom(handle).unwrap().0, model.topology().atom_ids()[0]);
    let mut foreign = TopologyEditor::new();
    let other = foreign.add_atom(bare_atom(2)).unwrap();
    assert!(map.atom(other).is_none());
}

#[test]
fn stereo_replacement_is_occurrence_local_and_rejects_changed_source() {
    let molecule = kekule::smiles::to_molecules("F[C@](Cl)(Br)I")
        .unwrap()
        .pop()
        .unwrap();
    assert_eq!(molecule.stereo_elements().count(), 1);
    let mut builder = Model::builder();
    builder
        .add_molecule(molecule.clone(), &Positions::zeros(5))
        .unwrap();
    builder
        .add_molecule(molecule.clone(), &Positions::zeros(5))
        .unwrap();
    let original = builder.build().unwrap();
    let instances = original
        .topology()
        .molecules()
        .map(|m| m.id())
        .collect::<Vec<_>>();
    let mut editor = original.edit();
    editor
        .replace_source_instance_stereo(instances[0], &[])
        .unwrap();
    let result = editor.finish().unwrap();
    let counts = result
        .topology()
        .molecules()
        .map(|m| m.molecule().stereo_elements().count())
        .collect::<Vec<_>>();
    assert_eq!(counts, vec![0, 1]);
    assert!(original
        .topology()
        .molecules()
        .all(|m| m.molecule().stereo_elements().count() == 1));
    let mut changed = original.edit();
    let id = changed
        .atom_handle(original.topology().atom_ids()[0])
        .unwrap();
    let mut a = changed.atom(id).unwrap().clone();
    a.isotope = Some(19);
    changed.replace_atom(id, a).unwrap();
    assert!(changed
        .replace_source_instance_stereo(instances[0], &[])
        .is_err());
}
