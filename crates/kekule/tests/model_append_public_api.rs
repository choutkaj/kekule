use kekule::core::{Atom, AtomId, BondOrder, Element};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyTable, PropertyValue};
use kekule::structure::{
    Ensemble, Model, ModelAppendCorrespondence, ModelEditError, ModelEditor, Positions,
};
use kekule::topology::{
    AtomSiteMetadata, InstanceAtomId, MoleculeClass, MoleculeInstanceId, ResidueClass,
    TopologyBuilder,
};
use kekule::units::{Quantity, ANGSTROM, KELVIN, NANOMETER, SQUARE_ANGSTROM};
use std::sync::Arc;

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
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
    Model::from_molecule(&molecule, &positions(molecule.atom_count())).unwrap()
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
        .insert_property(key("definition_note"), PropertyValue::String("keep".into()))
        .unwrap();
    for id in molecule.atom_ids() {
        chemistry
            .set_atom_property(
                id,
                key("tag"),
                Some(PropertyValue::Int(100 + id.index() as i64)),
            )
            .unwrap();
    }
    for id in molecule.bond_ids() {
        chemistry
            .set_bond_property(
                id,
                key("tag"),
                Some(PropertyValue::Int(200 + id.index() as i64)),
            )
            .unwrap();
    }
    let molecule = chemistry.finish().unwrap();
    assert!(molecule.perception().has_valence());
    assert!(molecule.perception().has_stereo());
    assert!(molecule.stereo_elements().count() > 0);
    let mut builder = TopologyBuilder::new();
    let shared = builder.add_molecule_definition(&molecule).unwrap();
    let separate = builder.add_molecule_definition(&molecule).unwrap();
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
        .molecule_instance_properties_mut()
        .insert(key("tag"), column(3, 300, conflict == Some("instance")))
        .unwrap();
    builder
        .atom_properties_mut()
        .insert(
            key("tag"),
            column(atom_count, 400, conflict == Some("topology atom")),
        )
        .unwrap();
    builder
        .bond_properties_mut()
        .insert(
            key("tag"),
            column(bond_count, 500, conflict == Some("topology bond")),
        )
        .unwrap();
    builder
        .chain_properties_mut()
        .insert(key("tag"), column(1, 600, conflict == Some("chain")))
        .unwrap();
    builder
        .residue_properties_mut()
        .insert(key("tag"), column(1, 700, conflict == Some("residue")))
        .unwrap();
    builder
        .atom_site_properties_mut()
        .insert(
            key("tag"),
            column(atom_count, 800, conflict == Some("atom site")),
        )
        .unwrap();
    builder
        .insert_property(
            key("system_note"),
            PropertyValue::String("original system".into()),
        )
        .unwrap();
    let mut model = Model::new(builder.build().unwrap(), positions(atom_count)).unwrap();
    model
        .insert_atom_property_column(
            key("tag"),
            column(atom_count, 900, conflict == Some("model atom")),
        )
        .unwrap();
    model
        .insert_bond_property_column(
            key("tag"),
            column(bond_count, 1000, conflict == Some("model bond")),
        )
        .unwrap();
    for (i, id) in model.topology().atom_ids().to_vec().into_iter().enumerate() {
        if i % 2 == 0 {
            model.set_occupancy(id, Some(0.75)).unwrap();
            model
                .set_b_factor(id, Some(Quantity::new(15.0 + i as f64, SQUARE_ANGSTROM)))
                .unwrap();
        }
    }
    model
        .insert_property(key("energy"), PropertyValue::Int(-12))
        .unwrap();
    model
}

fn assert_rows(source: &PropertyTable, target: &PropertyTable, pairs: &[(usize, usize)]) {
    for (key, _) in source.iter() {
        for &(from, to) in pairs {
            assert_eq!(
                source.value(key, from).unwrap(),
                target.value(key, to).unwrap(),
                "{key:?}, row {from}"
            );
        }
    }
}

fn assert_import(source: &Model, target: &Model, mapping: ModelAppendCorrespondence<'_>) {
    assert!(std::ptr::eq(mapping.source_topology(), source.topology()));
    assert!(std::ptr::eq(mapping.target_topology(), target.topology()));
    let atoms = source
        .topology()
        .atom_ids()
        .iter()
        .map(|&id| {
            let other = mapping.atom(id).unwrap();
            assert_eq!(source.atom(id).unwrap(), target.atom(other).unwrap());
            assert_eq!(
                source.position(id).unwrap(),
                target.position(other).unwrap()
            );
            assert_eq!(
                source.occupancy(id).unwrap(),
                target.occupancy(other).unwrap()
            );
            assert_eq!(
                source.b_factor(id).unwrap(),
                target.b_factor(other).unwrap()
            );
            (
                source.topology().atom_index(id).unwrap().index(),
                target.topology().atom_index(other).unwrap().index(),
            )
        })
        .collect::<Vec<_>>();
    assert_rows(source.atom_properties(), target.atom_properties(), &atoms);
    assert_rows(
        source.topology().atom_properties(),
        target.topology().atom_properties(),
        &atoms,
    );
    let bonds = source
        .topology()
        .bond_ids()
        .iter()
        .map(|&id| {
            let other = mapping.bond(id).unwrap();
            let before = source.bond(id).unwrap();
            let after = target.bond(other).unwrap();
            assert_eq!(before.order, after.order);
            assert_eq!(
                mapping
                    .atom(InstanceAtomId::new(id.molecule(), before.a()))
                    .unwrap(),
                InstanceAtomId::new(other.molecule(), after.a())
            );
            assert_eq!(
                mapping
                    .atom(InstanceAtomId::new(id.molecule(), before.b()))
                    .unwrap(),
                InstanceAtomId::new(other.molecule(), after.b())
            );
            (
                source.topology().bond_index(id).unwrap().index(),
                target.topology().bond_index(other).unwrap().index(),
            )
        })
        .collect::<Vec<_>>();
    assert_rows(source.bond_properties(), target.bond_properties(), &bonds);
    assert_rows(
        source.topology().bond_properties(),
        target.topology().bond_properties(),
        &bonds,
    );
    for (id, instance) in source.topology().instances() {
        let other = mapping.instances(id).unwrap();
        assert_eq!(other.len(), 1);
        let before = source.topology().definition(instance.definition()).unwrap();
        let after = target
            .topology()
            .definition(target.topology().instance(other[0]).unwrap().definition())
            .unwrap();
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
            source.topology().molecule_instance_properties(),
            target.topology().molecule_instance_properties(),
            &[(id.index(), other[0].index())],
        );
    }
    for (id, before) in source.hierarchy().chains() {
        let other = mapping.chain(id).unwrap();
        let after = target.hierarchy().chain(other).unwrap();
        assert_eq!(before.label_id(), after.label_id());
        assert_eq!(before.author_id(), after.author_id());
        assert_rows(
            source.topology().chain_properties(),
            target.topology().chain_properties(),
            &[(id.index(), other.index())],
        );
    }
    for (id, before) in source.hierarchy().residues() {
        let other = mapping.residue(id).unwrap();
        let after = target.hierarchy().residue(other).unwrap();
        assert_eq!(mapping.chain(before.chain()), Some(after.chain()));
        assert_eq!(before.name(), after.name());
        assert_eq!(before.label_seq_id(), after.label_seq_id());
        assert_eq!(before.author_seq_id(), after.author_seq_id());
        assert_eq!(before.insertion_code(), after.insertion_code());
        assert_eq!(before.label_comp_id(), after.label_comp_id());
        assert_eq!(before.author_comp_id(), after.author_comp_id());
        assert_eq!(before.class(), after.class());
        assert_rows(
            source.topology().residue_properties(),
            target.topology().residue_properties(),
            &[(id.index(), other.index())],
        );
    }
    for (id, before) in source.hierarchy().atom_sites() {
        let other = mapping.atom_site(id).unwrap();
        let after = target.hierarchy().atom_site(other).unwrap();
        assert_eq!(mapping.atom(before.atom()), Some(after.atom()));
        assert_eq!(mapping.residue(before.residue()), Some(after.residue()));
        assert_eq!(before.metadata(), after.metadata());
        assert_rows(
            source.topology().atom_site_properties(),
            target.topology().atom_site_properties(),
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
    let second = editor.append_model(source.view()).unwrap();
    assert!(second.report().cleared_model_properties.is_empty());
    assert!(second.report().cleared_topology_properties.is_empty());
    let result = editor.finish_with_correspondence().unwrap();
    let target = result.model();
    assert_eq!(target.atom_count(), source.atom_count() * 3);
    assert_eq!(target.topology().instance_count(), 9);
    assert_eq!(target.topology().definition_count(), 6);
    assert_eq!(target.chains().count(), 3); // Equal labels never merge chains.
    assert_eq!(target.residues().count(), 3);
    assert_eq!(target.atom_sites().count(), source.atom_sites().count() * 3);
    assert!(target.properties().owner_is_empty());
    assert!(target.topology().properties().owner_is_empty());
    assert_import(&source, target, first.published(&result).unwrap());
    assert_import(&source, target, second.published(&result).unwrap());
    for (before, after) in [
        (source.atom_properties(), target.atom_properties()),
        (source.bond_properties(), target.bond_properties()),
        (
            source.topology().molecule_instance_properties(),
            target.topology().molecule_instance_properties(),
        ),
        (
            source.topology().atom_properties(),
            target.topology().atom_properties(),
        ),
        (
            source.topology().bond_properties(),
            target.topology().bond_properties(),
        ),
        (
            source.topology().chain_properties(),
            target.topology().chain_properties(),
        ),
        (
            source.topology().residue_properties(),
            target.topology().residue_properties(),
        ),
        (
            source.topology().atom_site_properties(),
            target.topology().atom_site_properties(),
        ),
    ] {
        assert_rows(
            before,
            after,
            &(0..before.len()).map(|i| (i, i)).collect::<Vec<_>>(),
        );
    }
    for &id in source.topology().atom_ids() {
        assert_eq!(result.correspondence().target_atom(id), Some(id));
        assert_eq!(target.position(id).unwrap(), source.position(id).unwrap());
        assert_ne!(
            first.published(&result).unwrap().atom(id),
            second.published(&result).unwrap().atom(id)
        );
    }
    // Explicitly reused occurrences share a definition; equal independent ones do not.
    for mapping in [
        first.published(&result).unwrap(),
        second.published(&result).unwrap(),
    ] {
        let defs = source
            .topology()
            .instances()
            .map(|(id, _)| {
                target
                    .topology()
                    .instance(mapping.instances(id).unwrap()[0])
                    .unwrap()
                    .definition()
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
fn append_composes_with_generic_atom_bond_edits_and_tracks_splits_and_merges() {
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
    let linked = editor.clone().finish_with_correspondence().unwrap();
    assert_eq!(linked.model().topology().instance_count(), 1);
    let linked_mapping = imported.published(&linked).unwrap();
    assert_eq!(
        linked_mapping
            .instances(ligand_atoms[0].molecule())
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        linked
            .model()
            .position(linked_mapping.atom(ligand_atoms[1]).unwrap())
            .unwrap(),
        ligand.position(ligand_atoms[1]).unwrap()
    );
    editor
        .delete_bond(imported.bond(source_bond).unwrap())
        .unwrap();
    let split = editor.clone().finish_with_correspondence().unwrap();
    let mapping = imported.published(&split).unwrap();
    assert_eq!(
        mapping.instances(ligand_atoms[0].molecule()).unwrap().len(),
        2
    );
    assert!(mapping.bond(source_bond).is_none());
    editor.delete_atom(left).unwrap();
    editor.delete_atom(right).unwrap();
    let deleted = editor.finish_with_correspondence().unwrap();
    let mapping = imported.published(&deleted).unwrap();
    assert!(mapping.atom(ligand_atoms[0]).is_none());
    assert!(mapping.atom(ligand_atoms[1]).is_none());
    assert!(mapping
        .instances(ligand_atoms[0].molecule())
        .unwrap()
        .is_empty());
}

#[test]
fn import_mapping_rejects_foreign_publications_and_survives_clear_and_recovery() {
    let source = model("CO");
    let mut editor = source.edit();
    let imported = editor.append_model(&source).unwrap();
    let mut foreign = source.edit();
    foreign.append_model(&source).unwrap();
    let foreign = foreign.finish_with_correspondence().unwrap();
    assert!(matches!(
        imported.published(&foreign),
        Err(ModelEditError::ForeignAppend)
    ));
    let invalid = InstanceAtomId::new(MoleculeInstanceId::new(999), AtomId::new(0));
    assert!(imported.atom(invalid).is_err());
    editor.clear();
    let failure = editor.try_finish_with_correspondence().unwrap_err();
    let mut editor = failure.into_editor();
    editor
        .add_atom(atom("Na"), Quantity::new(Point3::origin(), NANOMETER))
        .unwrap();
    let result = editor.finish_with_correspondence().unwrap();
    let mapping = imported.published(&result).unwrap();
    assert!(mapping.atom(source.topology().atom_ids()[0]).is_none());
    assert!(mapping
        .instances(source.topology().atom_ids()[0].molecule())
        .unwrap()
        .is_empty());
    assert!(mapping.instances(MoleculeInstanceId::new(999)).is_err());
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
        let previous = editor.append_model(&destination).unwrap();
        let before = format!("{editor:?}");
        let error = editor.append_model(&incompatible).unwrap_err();
        assert!(error.to_string().contains(domain), "{error}");
        assert!(error.to_string().contains("tag"), "{error}");
        assert_eq!(format!("{editor:?}"), before, "rollback for {domain}");
        let imported = editor.append_model(&destination).unwrap();
        let result = editor.finish_with_correspondence().unwrap();
        assert_import(
            &destination,
            result.model(),
            previous.published(&result).unwrap(),
        );
        assert_import(
            &destination,
            result.model(),
            imported.published(&result).unwrap(),
        );
    }
}

#[test]
fn coordinates_and_property_units_are_preserved_with_missing_rows() {
    let mut destination = model("CC");
    let mut source = model("CO");
    let destination_atom = destination.topology().atom_ids()[0];
    let source_atom = source.topology().atom_ids()[0];
    destination
        .set_atom_property(
            destination_atom,
            key("distance"),
            Some(PropertyValue::real(10.0, ANGSTROM).unwrap()),
        )
        .unwrap();
    source
        .set_atom_property(
            source_atom,
            key("distance"),
            Some(PropertyValue::real(2.0, NANOMETER).unwrap()),
        )
        .unwrap();
    destination
        .set_atom_property(
            destination_atom,
            key("destination_only"),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    source
        .set_atom_property(
            source_atom,
            key("source_only"),
            Some(PropertyValue::String("note".into())),
        )
        .unwrap();
    source
        .set_positions(Quantity::new(
            vec![Point3::new(2.0, 3.0, 4.0), Point3::new(-1.0, 8.0, 0.0)],
            NANOMETER,
        ))
        .unwrap();
    let mut editor = destination.edit();
    let imported = editor.append_model(&source).unwrap();
    let result = editor.finish_with_correspondence().unwrap();
    let mapping = imported.published(&result).unwrap();
    let target = result.model();
    let new_atom = mapping.atom(source_atom).unwrap();
    assert_eq!(
        target.position(new_atom).unwrap(),
        source.position(source_atom).unwrap()
    );
    for (id, expected) in [(destination_atom, 1.0), (new_atom, 2.0)] {
        let Some(PropertyValue::Real { value, unit }) =
            target.atom_property(id, &key("distance")).unwrap()
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
    assert_eq!(
        target
            .atom_property(new_atom, &key("destination_only"))
            .unwrap(),
        None
    );
    assert_eq!(
        target
            .atom_property(destination_atom, &key("source_only"))
            .unwrap(),
        None
    );
    assert_eq!(
        target
            .atom_property(
                mapping.atom(source.topology().atom_ids()[1]).unwrap(),
                &key("distance")
            )
            .unwrap(),
        None
    );
    // An incompatible dimension rejects even though both columns are Real.
    source
        .remove_atom_property_column(&key("distance"))
        .unwrap();
    source
        .set_atom_property(
            source_atom,
            key("distance"),
            Some(PropertyValue::real(300.0, KELVIN).unwrap()),
        )
        .unwrap();
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
    source.set_cell(Some(cell(4.0, [true; 3])));
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
        source.set_cell(Some(wrong));
        let before = format!("{editor:?}");
        assert!(matches!(
            editor.append_model(&source),
            Err(ModelEditError::IncompatibleAppendCell)
        ));
        assert_eq!(format!("{editor:?}"), before);
    }
    source.set_cell(None);
    editor.append_model(&source).unwrap();
    let expected = *editor.cell().unwrap();
    assert_eq!(editor.finish().unwrap().cell(), Some(&expected));
}

#[test]
fn borrowed_ensemble_member_and_plain_finish_need_no_intermediate_model() {
    let source = annotated(None);
    let ensemble = Ensemble::from_models(std::slice::from_ref(&source)).unwrap();
    let member = ensemble.members().next().unwrap();
    let mut editor = ModelEditor::new();
    let imported = editor.append_model(member.as_model()).unwrap();
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
    assert_eq!(result.atom_properties(), source.atom_properties());
}

#[test]
fn editing_one_imported_occurrence_preserves_others_and_their_instance_annotations() {
    let source = annotated(None);
    let mut editor = ModelEditor::new();
    let imported = editor.append_model(&source).unwrap();
    let source_atom = source.topology().atom_ids()[0];
    let mut replacement = source.atom(source_atom).unwrap().clone();
    replacement.formal_charge = 1;
    editor
        .replace_atom(imported.atom(source_atom).unwrap(), replacement)
        .unwrap();
    let result = editor.finish_with_correspondence().unwrap();
    let mapping = imported.published(&result).unwrap();
    assert_eq!(result.model().topology().definition_count(), 3);
    for (id, instance) in source.topology().instances() {
        let target = mapping.instances(id).unwrap()[0];
        let before = source
            .topology()
            .definition(instance.definition())
            .unwrap()
            .molecule();
        let after = result
            .model()
            .topology()
            .definition(
                result
                    .model()
                    .topology()
                    .instance(target)
                    .unwrap()
                    .definition(),
            )
            .unwrap()
            .molecule();
        if id == source_atom.molecule() {
            assert_eq!(
                result
                    .model()
                    .topology()
                    .molecule_instance_properties()
                    .value(&key("tag"), target.index())
                    .unwrap(),
                None
            );
            assert!(!after.perception().has_valence());
        } else {
            assert_eq!(before, after);
            assert_eq!(before.perception(), after.perception());
            assert_rows(
                source.topology().molecule_instance_properties(),
                result.model().topology().molecule_instance_properties(),
                &[(id.index(), target.index())],
            );
        }
    }
}

#[test]
fn sparse_source_ids_and_reordered_definition_occurrences_keep_dense_state_associated() {
    let molecule = kekule::smiles::to_molecules("CCO").unwrap().pop().unwrap();
    let mut edit = molecule.edit();
    edit.delete_atom(AtomId::new(0)).unwrap();
    let sparse = edit.finish().unwrap();
    assert_eq!(sparse.atom_ids().next(), Some(AtomId::new(1)));
    let sodium = kekule::smiles::to_molecules("[Na+]")
        .unwrap()
        .pop()
        .unwrap();
    let mut builder = TopologyBuilder::new();
    let organic = builder.add_molecule_definition(&sparse).unwrap();
    let ion = builder.add_molecule_definition(&sodium).unwrap();
    for definition in [ion, organic, organic] {
        builder.add_instance(definition).unwrap();
    }
    builder
        .atom_properties_mut()
        .insert(key("tag"), column(5, 200, false))
        .unwrap();
    builder
        .bond_properties_mut()
        .insert(key("tag"), column(2, 300, false))
        .unwrap();
    let mut source = Model::new(builder.build().unwrap(), positions(5)).unwrap();
    source
        .insert_atom_property_column(key("tag"), column(5, 400, false))
        .unwrap();
    source
        .insert_bond_property_column(key("tag"), column(2, 500, false))
        .unwrap();
    let mut destination = model("CCC").into_editor();
    let deleted = destination.atom_ids().next().unwrap();
    destination.delete_atom(deleted).unwrap();
    let imported = destination.append_model(&source).unwrap();
    let result = destination.finish_with_correspondence().unwrap();
    assert_import(
        &source,
        result.model(),
        imported.published(&result).unwrap(),
    );
    assert_eq!(result.model().topology().definition_count(), 3);
}
