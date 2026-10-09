use kekule::core::{Atom, AtomId, BondOrder, Element};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyTable, PropertyValue};
use kekule::structure::{Ensemble, Model, ModelEditError, ModelEditor, Positions};
use kekule::topology::{
    AtomSiteMetadata, InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeInstanceId,
    ResidueClass, TopologyAtomIndex, TopologyBuilder,
};
use kekule::units::{Quantity, ANGSTROM, KELVIN, NANOMETER, SQUARE_ANGSTROM};
use std::sync::Arc;

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}
fn index(model: &Model, atom: InstanceAtomId) -> TopologyAtomIndex {
    model.topology().atom_index(atom).unwrap()
}
fn occupancy(model: &Model, atom: InstanceAtomId) -> Option<f64> {
    model.occupancy(index(model, atom)).unwrap()
}
fn b_factor(model: &Model, atom: InstanceAtomId) -> Option<Quantity<f64>> {
    model.b_factor(index(model, atom)).unwrap()
}
fn atom_value(model: &Model, atom: InstanceAtomId, name: &PropertyKey) -> Option<PropertyValue> {
    model
        .properties()
        .atoms()
        .value(name, index(model, atom))
        .unwrap()
}
fn set_atom_value(
    model: &mut Model,
    atom: InstanceAtomId,
    name: PropertyKey,
    value: PropertyValue,
) {
    let row = index(model, atom);
    model
        .conformation_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(name, row, Some(value))
        .unwrap();
}
fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}
fn positions(count: usize) -> Positions {
    Positions::new(Quantity::new(
        (0..count)
            .map(|i| Point3::new(i as f64 + 0.25, 2.0, -1.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap()
}
fn model(smiles: &str) -> Model {
    let molecule = kekule::smiles::to_molecules(smiles).unwrap().pop().unwrap();
    Model::from_molecule(molecule.clone(), &positions(molecule.atom_count())).unwrap()
}
fn column(len: usize, offset: i64, conflict: bool) -> PropertyColumn {
    let values = (0..len)
        .map(|i| (i % 3 != 1).then_some(offset + i as i64))
        .collect::<Vec<_>>();
    if conflict {
        PropertyColumn::String(
            values
                .into_iter()
                .map(|v| v.map(|v| v.to_string()))
                .collect(),
        )
    } else {
        PropertyColumn::Int(values)
    }
}

// Focused synthetic regression: repeated stereochemical definitions and hierarchy
// that crosses molecular boundaries, with annotations in every supported domain.
fn annotated(conflict: Option<&str>) -> Model {
    let mut molecule = kekule::smiles::to_molecules("N[C@@H](C)C(=O)O")
        .unwrap()
        .pop()
        .unwrap();
    molecule.perceive().unwrap();
    kekule::stereo::assign_cip_descriptors(&mut molecule).unwrap();
    let mut chemistry = molecule.edit();
    chemistry
        .properties_mut()
        .owner_mut()
        .insert(key("definition_note"), PropertyValue::String("keep".into()))
        .unwrap();
    for id in molecule.atom_ids() {
        chemistry
            .properties_mut()
            .atoms_mut()
            .set_value(
                key("tag"),
                id,
                Some(PropertyValue::Int(100 + id.index() as i64)),
            )
            .unwrap();
    }
    for id in molecule.bond_ids() {
        chemistry
            .properties_mut()
            .bonds_mut()
            .set_value(
                key("tag"),
                id,
                Some(PropertyValue::Int(200 + id.index() as i64)),
            )
            .unwrap();
    }
    let molecule = chemistry.finish().unwrap();
    assert!(molecule.perception().has_valence());
    assert!(molecule.perception().has_stereo());
    assert!(molecule.stereo_elements().count() > 0);
    let mut builder = TopologyBuilder::new();
    let shared = builder.add_molecule_definition(molecule.clone()).unwrap();
    let separate = builder.add_molecule_definition(molecule.clone()).unwrap();
    builder
        .set_molecule_class(shared, MoleculeClass::Other)
        .unwrap();
    builder
        .set_molecule_class(separate, MoleculeClass::SmallMolecule)
        .unwrap();
    let instances = [shared, shared, separate].map(|id| builder.add_instance(id).unwrap());
    let chain = builder
        .hierarchy_mut()
        .add_chain("A", Some("source-A".into()))
        .unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "LIG", Some(9), Some("42".into()), Some("B".into()))
        .unwrap();
    builder
        .hierarchy_mut()
        .set_residue_component_ids(residue, Some("ALA".into()), Some("ALT".into()))
        .unwrap();
    for instance in instances {
        for (id, value) in molecule.atoms() {
            builder
                .hierarchy_mut()
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, id),
                    AtomSiteMetadata {
                        type_symbol: Some(value.element.symbol().into()),
                        label_asym_id: Some("A".into()),
                        auth_asym_id: Some("source-A".into()),
                        label_atom_id: Some(format!("L{}", id.index())),
                        auth_atom_id: Some(format!("X{}", id.index())),
                    },
                )
                .unwrap();
        }
    }
    builder
        .set_residue_class(residue, ResidueClass::Other)
        .unwrap();
    let atom_count = molecule.atom_count() * 3;
    let bond_count = molecule.bond_count() * 3;
    builder
        .properties_mut()
        .molecule_instances_mut()
        .insert(key("tag"), column(3, 300, conflict == Some("instance")))
        .unwrap();
    builder
        .properties_mut()
        .atoms_mut()
        .insert(
            key("tag"),
            column(atom_count, 400, conflict == Some("topology atom")),
        )
        .unwrap();
    builder
        .properties_mut()
        .bonds_mut()
        .insert(
            key("tag"),
            column(bond_count, 500, conflict == Some("topology bond")),
        )
        .unwrap();
    builder
        .properties_mut()
        .chains_mut()
        .insert(key("tag"), column(1, 600, conflict == Some("chain")))
        .unwrap();
    builder
        .properties_mut()
        .residues_mut()
        .insert(key("tag"), column(1, 700, conflict == Some("residue")))
        .unwrap();
    builder
        .properties_mut()
        .atom_sites_mut()
        .insert(
            key("tag"),
            column(atom_count, 800, conflict == Some("atom site")),
        )
        .unwrap();
    builder
        .properties_mut()
        .owner_mut()
        .insert(
            key("system_note"),
            PropertyValue::String("original system".into()),
        )
        .unwrap();
    let mut model = Model::new(builder.build().unwrap(), positions(atom_count)).unwrap();
    let mut conformation = model.conformation_mut();
    conformation
        .properties_mut()
        .atoms_mut()
        .insert(
            key("tag"),
            column(atom_count, 900, conflict == Some("model atom")),
        )
        .unwrap();
    conformation
        .properties_mut()
        .bonds_mut()
        .insert(
            key("tag"),
            column(bond_count, 1000, conflict == Some("model bond")),
        )
        .unwrap();
    for i in (0..atom_count).step_by(2) {
        let row = TopologyAtomIndex::new(i as u32);
        conformation.set_occupancy(row, Some(0.75)).unwrap();
        conformation
            .set_b_factor(row, Some(Quantity::new(15.0 + i as f64, SQUARE_ANGSTROM)))
            .unwrap();
    }
    model
        .conformation_mut()
        .properties_mut()
        .owner_mut()
        .insert(key("energy"), PropertyValue::Int(-12))
        .unwrap();
    model
}

fn table_cell<R>(
    table: &PropertyTable<R>,
    name: &PropertyKey,
    row: usize,
) -> Option<PropertyValue> {
    table
        .get(name)
        .and_then(|column| column.value(row).unwrap())
}

fn assert_rows<R>(source: &PropertyTable<R>, target: &PropertyTable<R>, pairs: &[(usize, usize)]) {
    for name in source.keys() {
        for &(from, to) in pairs {
            assert_eq!(
                table_cell(source, name, from),
                table_cell(target, name, to),
                "{name:?}, row {from}"
            );
        }
    }
}

fn assert_same_rows<R>(source: &PropertyTable<R>, target: &PropertyTable<R>) {
    assert_rows(
        source,
        target,
        &(0..source.len()).map(|i| (i, i)).collect::<Vec<_>>(),
    );
}

// Append-only fixture checks: occurrences retain their local chemistry IDs and
// appear after the destination occurrences. This is an expected layout, not an
// editor-supplied correspondence; every row, endpoint and hierarchy link is checked.
fn assert_import(source: &Model, target: &Model, first_instance: usize) {
    let instance =
        |id: MoleculeInstanceId| MoleculeInstanceId::new((first_instance + id.index()) as u32);
    let atom = |id: InstanceAtomId| InstanceAtomId::new(instance(id.molecule()), id.atom());
    let atoms = source
        .topology()
        .atom_ids()
        .iter()
        .map(|&id| {
            let other = atom(id);
            assert_eq!(
                source.topology().atom(id).unwrap().atom(),
                target.topology().atom(other).unwrap().atom()
            );
            assert_eq!(
                source.position(id).unwrap(),
                target.position(other).unwrap()
            );
            assert_eq!(occupancy(source, id), occupancy(target, other));
            assert_eq!(b_factor(source, id), b_factor(target, other));
            let before = source.topology().hierarchy().atom_site_for_atom(id);
            let after = target.topology().hierarchy().atom_site_for_atom(other);
            assert_eq!(before.is_some(), after.is_some());
            if let (Some(before), Some(after)) = (before, after) {
                assert_eq!(before.metadata(), after.metadata());
                assert_rows(
                    source.topology().properties().atom_sites(),
                    target.topology().properties().atom_sites(),
                    &[(before.id().index(), after.id().index())],
                );
                let before = source
                    .topology()
                    .hierarchy()
                    .residue(before.residue())
                    .unwrap();
                let after = target
                    .topology()
                    .hierarchy()
                    .residue(after.residue())
                    .unwrap();
                assert_eq!(before.atom_sites().len(), after.atom_sites().len());
                assert_eq!(before.name(), after.name());
                assert_eq!(before.label_seq_id(), after.label_seq_id());
                assert_eq!(before.author_seq_id(), after.author_seq_id());
                assert_eq!(before.insertion_code(), after.insertion_code());
                assert_eq!(before.label_comp_id(), after.label_comp_id());
                assert_eq!(before.author_comp_id(), after.author_comp_id());
                assert_eq!(before.class(), after.class());
                assert_rows(
                    source.topology().properties().residues(),
                    target.topology().properties().residues(),
                    &[(before.id().index(), after.id().index())],
                );
                let before = source.topology().hierarchy().chain(before.chain()).unwrap();
                let after = target.topology().hierarchy().chain(after.chain()).unwrap();
                assert_eq!(before.residues().len(), after.residues().len());
                assert_eq!(before.label_id(), after.label_id());
                assert_eq!(before.author_id(), after.author_id());
                assert_rows(
                    source.topology().properties().chains(),
                    target.topology().properties().chains(),
                    &[(before.id().index(), after.id().index())],
                );
            }
            (
                source.topology().atom_index(id).unwrap().index(),
                target.topology().atom_index(other).unwrap().index(),
            )
        })
        .collect::<Vec<_>>();
    assert_rows(
        source.properties().atoms(),
        target.properties().atoms(),
        &atoms,
    );
    assert_rows(
        source.topology().properties().atoms(),
        target.topology().properties().atoms(),
        &atoms,
    );
    let bonds = source
        .topology()
        .bond_ids()
        .iter()
        .map(|&id| {
            let other = InstanceBondId::new(instance(id.molecule()), id.bond());
            let before = source.topology().bond(id).unwrap();
            let after = target.topology().bond(other).unwrap();
            assert_eq!(before.order, after.order);
            assert_eq!(
                atom(InstanceAtomId::new(id.molecule(), before.a())),
                InstanceAtomId::new(other.molecule(), after.a())
            );
            assert_eq!(
                atom(InstanceAtomId::new(id.molecule(), before.b())),
                InstanceAtomId::new(other.molecule(), after.b())
            );
            (
                source.topology().bond_index(id).unwrap().index(),
                target.topology().bond_index(other).unwrap().index(),
            )
        })
        .collect::<Vec<_>>();
    assert_rows(
        source.properties().bonds(),
        target.properties().bonds(),
        &bonds,
    );
    assert_rows(
        source.topology().properties().bonds(),
        target.topology().properties().bonds(),
        &bonds,
    );
    for molecule in source.topology().molecules() {
        let id = molecule.id();
        let other = instance(id);
        let before = molecule.definition();
        let after = target.topology().molecule(other).unwrap().definition();
        assert_eq!(before.class(), after.class());
        assert_eq!(before.molecule(), after.molecule());
        assert_eq!(
            before.molecule().properties(),
            after.molecule().properties()
        );
        assert_eq!(
            before.molecule().perception(),
            after.molecule().perception()
        );
        assert_rows(
            source.topology().properties().molecule_instances(),
            target.topology().properties().molecule_instances(),
            &[(id.index(), other.index())],
        );
    }
}

#[test]
fn complete_import_preserves_every_scope_and_explicit_definition_reuse() {
    let source = annotated(None);
    let before = format!("{source:?}");
    let mut editor = source.edit();
    let first = editor.append_model(&source).unwrap();
    assert_eq!(
        first.report().cleared_topology_properties,
        vec![key("system_note")]
    );
    assert_eq!(first.report().cleared_model_properties, vec![key("energy")]);
    assert_eq!(
        first.report().omitted_topology_properties,
        vec![key("system_note")]
    );
    assert_eq!(first.report().omitted_model_properties, vec![key("energy")]);
    let second = editor.append_model(&source.as_model_view()).unwrap();
    assert!(second.report().cleared_model_properties.is_empty());
    assert!(second.report().cleared_topology_properties.is_empty());
    let result = editor.finish().unwrap();
    let target = &result;
    assert_eq!(target.atom_count(), source.atom_count() * 3);
    assert_eq!(target.topology().instance_count(), 9);
    assert_eq!(target.topology().definition_count(), 6);
    assert_eq!(target.topology().chains().count(), 3); // Equal labels never merge chains.
    assert_eq!(target.topology().residues().count(), 3);
    assert_eq!(
        target.topology().atom_sites().count(),
        source.topology().atom_sites().count() * 3
    );
    assert!(target.properties().owner().is_empty());
    assert!(target.topology().properties().owner().is_empty());
    assert_import(&source, target, 3);
    assert_import(&source, target, 6);
    assert_same_rows(source.properties().atoms(), target.properties().atoms());
    assert_same_rows(source.properties().bonds(), target.properties().bonds());
    let (from, to) = (
        source.topology().properties(),
        target.topology().properties(),
    );
    assert_same_rows(from.molecule_instances(), to.molecule_instances());
    assert_same_rows(from.atoms(), to.atoms());
    assert_same_rows(from.bonds(), to.bonds());
    assert_same_rows(from.chains(), to.chains());
    assert_same_rows(from.residues(), to.residues());
    assert_same_rows(from.atom_sites(), to.atom_sites());
    for &id in source.topology().atom_ids() {
        assert_eq!(target.position(id).unwrap(), source.position(id).unwrap());
        assert_ne!(first.atom(id).unwrap(), second.atom(id).unwrap());
    }
    // Explicitly reused occurrences share a definition; equal independent ones do not.
    for first_instance in [3, 6] {
        let defs = source
            .topology()
            .molecules()
            .map(|molecule| {
                target
                    .topology()
                    .molecule(MoleculeInstanceId::new(
                        first_instance + molecule.id().index() as u32,
                    ))
                    .unwrap()
                    .definition_id()
            })
            .collect::<Vec<_>>();
        assert_eq!(defs[0], defs[1]);
        assert_ne!(defs[0], defs[2]);
    }
    assert_eq!(format!("{source:?}"), before);
    assert!(!Arc::ptr_eq(
        &target.shared_topology(),
        &source.shared_topology()
    ));
}

#[test]
fn append_composes_with_generic_atom_bond_edits_and_preserves_geometry_through_splits_and_merges() {
    let receptor = model("CC");
    let ligand = model("CO");
    let receptor_atom = receptor.topology().atom_ids()[0];
    let ligand_atoms = ligand.topology().atom_ids();
    let source_bond = ligand.topology().bond_ids()[0];
    let mut editor = receptor.edit();
    // Existing deleted slots must not shift imported coordinates or property rows.
    editor
        .delete_atom(
            editor
                .atom_handle(receptor.topology().atom_ids()[1])
                .unwrap(),
        )
        .unwrap();
    let imported = editor.append_model(&ligand).unwrap();
    let left = imported.atom(ligand_atoms[0]).unwrap();
    let right = imported.atom(ligand_atoms[1]).unwrap();
    let hydrogen = editor
        .add_atom(
            atom("H"),
            Quantity::new(Point3::new(9.0, 1.0, 2.0), ANGSTROM),
        )
        .unwrap();
    editor.add_bond(left, hydrogen, BondOrder::Single).unwrap();
    editor
        .add_bond(
            editor.atom_handle(receptor_atom).unwrap(),
            left,
            BondOrder::Single,
        )
        .unwrap();
    let linked = editor.clone().finish().unwrap();
    assert_eq!(linked.topology().instance_count(), 1);
    assert_eq!(linked.atom_count(), 4);
    assert_eq!(linked.topology().bond_ids().len(), 3);
    let oxygen = linked
        .topology()
        .atoms()
        .find(|a| a.element.symbol() == "O")
        .unwrap()
        .id();
    assert_eq!(
        linked.position(oxygen).unwrap(),
        ligand.position(ligand_atoms[1]).unwrap()
    );
    let h = linked
        .topology()
        .atoms()
        .find(|a| a.element.symbol() == "H")
        .unwrap()
        .id();
    assert_eq!(
        linked.position(h).unwrap(),
        Quantity::new(Point3::new(9.0, 1.0, 2.0), ANGSTROM)
            .into_unit(NANOMETER)
            .unwrap()
    );
    let removed = imported.bond(source_bond).unwrap();
    editor.delete_bond(removed).unwrap();
    assert!(editor.bond(removed).is_err());
    let split = editor.clone().finish().unwrap();
    assert_eq!(split.topology().instance_count(), 2);
    assert_eq!(split.atom_count(), 4);
    assert_eq!(split.topology().bond_ids().len(), 2);
    let oxygen = split
        .topology()
        .atoms()
        .find(|a| a.element.symbol() == "O")
        .unwrap()
        .id();
    let oxygen_molecule = split
        .topology()
        .molecule(oxygen.molecule())
        .unwrap()
        .molecule();
    assert_eq!(oxygen_molecule.atom_count(), 1);
    assert_eq!(
        split.position(oxygen).unwrap(),
        ligand.position(ligand_atoms[1]).unwrap()
    );
    editor.delete_atom(left).unwrap();
    editor.delete_atom(right).unwrap();
    assert!(editor.atom(left).is_err());
    assert!(editor.atom(right).is_err());
    let deleted = editor.finish().unwrap();
    assert_eq!(deleted.atom_count(), 2);
    assert_eq!(deleted.topology().instance_count(), 2);
    assert_eq!(deleted.topology().bond_ids().len(), 0);
    assert_eq!(
        deleted
            .topology()
            .atoms()
            .map(|a| a.element.symbol())
            .collect::<Vec<_>>(),
        ["C", "H"]
    );
    assert_eq!(
        deleted.position(deleted.topology().atom_ids()[0]).unwrap(),
        receptor.position(receptor_atom).unwrap()
    );
    assert_eq!(
        deleted.position(deleted.topology().atom_ids()[1]).unwrap(),
        linked.position(h).unwrap()
    );
}

#[test]
fn import_handles_reject_foreign_drafts_and_remain_invalid_after_clear_and_recovery() {
    let source = model("CO");
    let mut editor = source.edit();
    let imported = editor.append_model(&source).unwrap();
    let mut foreign = source.edit();
    foreign.append_model(&source).unwrap();
    let imported_atom = imported.atom(source.topology().atom_ids()[0]).unwrap();
    assert!(foreign.position(imported_atom).is_err());
    assert!(foreign.delete_atom(imported_atom).is_err());
    let invalid = InstanceAtomId::new(MoleculeInstanceId::new(999), AtomId::new(0));
    assert!(imported.atom(invalid).is_err());
    editor.clear();
    let failure = editor.try_finish().unwrap_err();
    let mut editor = failure.into_editor();
    editor
        .add_atom(atom("Na"), Quantity::new(Point3::origin(), NANOMETER))
        .unwrap();
    assert!(editor.position(imported_atom).is_err());
    assert!(editor.delete_atom(imported_atom).is_err());
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
        "Na"
    );
    assert_eq!(
        result
            .position(result.topology().atom_ids()[0])
            .unwrap()
            .into_value(),
        Point3::origin()
    );
}

#[test]
fn all_property_domain_conflicts_roll_back_and_can_be_retried() {
    for domain in [
        "instance",
        "topology atom",
        "topology bond",
        "chain",
        "residue",
        "atom site",
        "model atom",
        "model bond",
    ] {
        let destination = annotated(None);
        let incompatible = annotated(Some(domain));
        let mut editor = destination.edit();
        editor.append_model(&destination).unwrap();
        let before = format!("{editor:?}");
        let error = editor.append_model(&incompatible).unwrap_err();
        assert!(error.to_string().contains(domain), "{error}");
        assert!(error.to_string().contains("tag"), "{error}");
        assert_eq!(format!("{editor:?}"), before, "rollback for {domain}");
        editor.append_model(&destination).unwrap();
        let result = editor.finish().unwrap();
        assert_import(&destination, &result, 3);
        assert_import(&destination, &result, 6);
    }
}

#[test]
fn coordinates_and_property_units_are_preserved_with_missing_rows() {
    let mut destination = model("CC");
    let mut source = model("CO");
    let destination_atom = destination.topology().atom_ids()[0];
    let source_atom = source.topology().atom_ids()[0];
    set_atom_value(
        &mut destination,
        destination_atom,
        key("distance"),
        PropertyValue::real(10.0, ANGSTROM).unwrap(),
    );
    set_atom_value(
        &mut source,
        source_atom,
        key("distance"),
        PropertyValue::real(2.0, NANOMETER).unwrap(),
    );
    set_atom_value(
        &mut destination,
        destination_atom,
        key("destination_only"),
        PropertyValue::Bool(true),
    );
    set_atom_value(
        &mut source,
        source_atom,
        key("source_only"),
        PropertyValue::String("note".into()),
    );
    source
        .conformation_mut()
        .set_positions(Quantity::new(
            vec![Point3::new(2.0, 3.0, 4.0), Point3::new(-1.0, 8.0, 0.0)],
            NANOMETER,
        ))
        .unwrap();
    let mut editor = destination.edit();
    editor.append_model(&source).unwrap();
    let result = editor.finish().unwrap();
    let target = &result;
    let new_atom = target.topology().atom_ids()[2];
    assert_eq!(
        target.position(new_atom).unwrap(),
        source.position(source_atom).unwrap()
    );
    for (id, expected) in [(destination_atom, 1.0), (new_atom, 2.0)] {
        let Some(PropertyValue::Real { value, unit }) = atom_value(target, id, &key("distance"))
        else {
            panic!("real value")
        };
        assert!(
            (Quantity::new(value, unit)
                .into_unit(NANOMETER)
                .unwrap()
                .into_value()
                - expected)
                .abs()
                < 1e-12
        );
    }
    assert_eq!(atom_value(target, new_atom, &key("destination_only")), None);
    assert_eq!(
        atom_value(target, destination_atom, &key("source_only")),
        None
    );
    assert_eq!(
        atom_value(target, target.topology().atom_ids()[3], &key("distance")),
        None
    );
    // An incompatible dimension rejects even though both columns are Real.
    source
        .conformation_mut()
        .properties_mut()
        .atoms_mut()
        .remove(&key("distance"))
        .unwrap();
    set_atom_value(
        &mut source,
        source_atom,
        key("distance"),
        PropertyValue::real(300.0, KELVIN).unwrap(),
    );
    let mut editor = destination.edit();
    let before = format!("{editor:?}");
    assert!(matches!(
        editor.append_model(&source),
        Err(ModelEditError::AppendProperty {
            domain: "model atom",
            ..
        })
    ));
    assert_eq!(format!("{editor:?}"), before);
}

fn cell(size: f64, axes: [bool; 3]) -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(size, size, size), NANOMETER),
        axes,
    )
    .unwrap()
}

#[test]
fn periodic_cell_rules_are_explicit_and_transactional() {
    let mut source = model("CO");
    source
        .conformation_mut()
        .set_cell(Some(cell(4.0, [true; 3])));
    let mut empty = ModelEditor::new();
    empty.append_model(&source).unwrap();
    assert_eq!(empty.finish().unwrap().cell(), source.cell());
    let destination = model("C");
    let mut editor = destination.edit();
    let before = format!("{editor:?}");
    assert!(matches!(
        editor.append_model(&source),
        Err(ModelEditError::IncompatibleAppendCell)
    ));
    assert_eq!(format!("{editor:?}"), before);
    editor.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(40.0, 40.0, 40.0), ANGSTROM),
            [true; 3],
        )
        .unwrap(),
    ));
    editor.append_model(&source).unwrap();
    for wrong in [
        cell(5.0, [true; 3]),
        cell(4.0 + 1e-10, [true; 3]),
        cell(4.0, [true, true, false]),
    ] {
        source.conformation_mut().set_cell(Some(wrong));
        let before = format!("{editor:?}");
        assert!(matches!(
            editor.append_model(&source),
            Err(ModelEditError::IncompatibleAppendCell)
        ));
        assert_eq!(format!("{editor:?}"), before);
    }
    source.conformation_mut().set_cell(None);
    editor.append_model(&source).unwrap();
    let expected = *editor.cell().unwrap();
    assert_eq!(editor.finish().unwrap().cell(), Some(&expected));
}

#[test]
fn borrowed_ensemble_member_and_plain_finish_need_no_intermediate_model() {
    let source = annotated(None);
    let ensemble = Ensemble::from_models([source.clone()]).unwrap();
    let member = ensemble.iter().next().unwrap();
    let mut editor = ModelEditor::new();
    let imported = editor.append_model(&member).unwrap();
    for &id in source.topology().atom_ids() {
        assert_eq!(
            editor.position(imported.atom(id).unwrap()).unwrap(),
            source.position(id).unwrap()
        );
    }
    let result = editor.finish().unwrap();
    assert_eq!(result.atom_count(), source.atom_count());
    assert_eq!(
        result.topology().definition_count(),
        source.topology().definition_count()
    );
    assert_eq!(result.properties().atoms(), source.properties().atoms());
}

#[test]
fn editing_one_imported_occurrence_preserves_others_and_their_instance_annotations() {
    let source = annotated(None);
    let mut editor = ModelEditor::new();
    let imported = editor.append_model(&source).unwrap();
    let source_atom = source.topology().atom_ids()[0];
    let mut replacement = source.topology().atom(source_atom).unwrap().atom().clone();
    replacement.formal_charge = 1;
    editor
        .replace_atom(imported.atom(source_atom).unwrap(), replacement)
        .unwrap();
    let result = editor.finish().unwrap();
    assert_eq!(result.topology().definition_count(), 3);
    assert_eq!(result.positions(), source.positions());
    assert_eq!(
        result.topology().atom(source_atom).unwrap().formal_charge,
        1
    );
    for molecule in source.topology().molecules() {
        let id = molecule.id();
        let target = id; // No membership change: occurrences retain their input order.
        let before = molecule.molecule();
        let after = result.topology().molecule(target).unwrap().molecule();
        if id == source_atom.molecule() {
            assert_eq!(
                result
                    .topology()
                    .properties()
                    .molecule_instances()
                    .value(&key("tag"), target)
                    .unwrap(),
                None
            );
            assert!(!after.perception().has_valence());
        } else {
            assert_eq!(before, after);
            assert_eq!(before.perception(), after.perception());
            assert_rows(
                source.topology().properties().molecule_instances(),
                result.topology().properties().molecule_instances(),
                &[(id.index(), target.index())],
            );
        }
    }
}

#[test]
fn edited_source_ids_and_reordered_definition_occurrences_keep_dense_state_associated() {
    let molecule = kekule::smiles::to_molecules("CCO").unwrap().pop().unwrap();
    let mut edit = molecule.edit();
    edit.delete_atom(AtomId::new(0)).unwrap();
    let sparse = edit.finish().unwrap();
    // Publication renumbers the surviving atoms densely.
    assert_eq!(
        sparse.atom_ids().collect::<Vec<_>>(),
        [AtomId::new(0), AtomId::new(1)]
    );
    let sodium = kekule::smiles::to_molecules("[Na+]")
        .unwrap()
        .pop()
        .unwrap();
    let mut builder = TopologyBuilder::new();
    let organic = builder.add_molecule_definition(sparse.clone()).unwrap();
    let ion = builder.add_molecule_definition(sodium.clone()).unwrap();
    for definition in [ion, organic, organic] {
        builder.add_instance(definition).unwrap();
    }
    builder
        .properties_mut()
        .atoms_mut()
        .insert(key("tag"), column(5, 200, false))
        .unwrap();
    builder
        .properties_mut()
        .bonds_mut()
        .insert(key("tag"), column(2, 300, false))
        .unwrap();
    let mut source = Model::new(builder.build().unwrap(), positions(5)).unwrap();
    let mut conformation = source.conformation_mut();
    let mut properties = conformation.properties_mut();
    properties
        .atoms_mut()
        .insert(key("tag"), column(5, 400, false))
        .unwrap();
    properties
        .bonds_mut()
        .insert(key("tag"), column(2, 500, false))
        .unwrap();
    let mut destination = model("CCC").into_editor();
    let deleted = destination.atom_ids().next().unwrap();
    destination.delete_atom(deleted).unwrap();
    destination.append_model(&source).unwrap();
    let result = destination.finish().unwrap();
    assert_import(&source, &result, 1);
    assert_eq!(result.topology().definition_count(), 3);
}
