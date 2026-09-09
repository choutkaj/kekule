use std::error::Error;
use std::sync::Arc;

use kekule::core::{Atom, Element, Molecule, MoleculeEditor, MoleculePublicationError};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::topology::transform::{TopologySubsetError, TopologyTransformError};
use kekule::topology::{
    AtomSelection, AtomSiteMetadata, InstanceAtomId, MoleculeClass, ResidueClass, ResidueId,
    SelectionError, Topology, TopologyBuildError, TopologyBuilder, TopologyEditor,
};

fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}

fn oxygen() -> Molecule {
    let mut editor = MoleculeEditor::new();
    editor.add_atom(atom("O")).unwrap();
    editor.finish().unwrap()
}

fn water_builder(label: &str) -> (TopologyBuilder, ResidueId) {
    let molecule = oxygen();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "HOH", None, None, None)
        .unwrap();
    builder
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, molecule.atom_ids().next().unwrap()),
            AtomSiteMetadata::default(),
        )
        .unwrap();
    builder
        .hierarchy_mut()
        .set_residue_component_ids(residue, Some(label.into()), None)
        .unwrap();
    (builder, residue)
}

fn classes(topology: &Topology) -> (MoleculeClass, ResidueClass) {
    (
        topology.molecules().next().unwrap().class(),
        topology.residues().next().unwrap().class(),
    )
}

#[test]
fn represented_equality_ignores_adjacency_history_but_not_bond_chemistry() {
    let molecule = kekule::smiles::to_molecules("CCC").unwrap().pop().unwrap();
    let atoms = molecule.atom_ids().collect::<Vec<_>>();
    let bond = molecule.bond_ids().next().unwrap();
    let mut editor = molecule.edit();
    editor.set_bond_endpoints(bond, atoms[0], atoms[2]).unwrap();
    assert_ne!(editor.clone().finish().unwrap(), molecule);
    editor.set_bond_endpoints(bond, atoms[0], atoms[1]).unwrap();
    let rewired = editor.finish().unwrap();
    assert!(molecule.atoms().eq(rewired.atoms()));
    assert!(molecule.bonds().eq(rewired.bonds()));
    assert_ne!(
        molecule
            .incident_bonds(atoms[1])
            .unwrap()
            .map(|(id, _)| id)
            .collect::<Vec<_>>(),
        rewired
            .incident_bonds(atoms[1])
            .unwrap()
            .map(|(id, _)| id)
            .collect::<Vec<_>>()
    );
    assert_eq!(molecule, rewired);
    assert!(Topology::from_molecule(&molecule)
        .unwrap()
        .same_layout(&Topology::from_molecule(&rewired).unwrap()));
}

#[test]
fn resumed_builder_reinfers_changed_component_identity_like_fresh_construction() {
    let (builder, residue) = water_builder("HOH");
    let source = builder.build().unwrap();
    assert_eq!(
        classes(&source),
        (MoleculeClass::Water, ResidueClass::Water)
    );
    let mut builder = source.into_builder();
    builder
        .hierarchy_mut()
        .set_residue_component_ids(residue, Some("UNK".into()), None)
        .unwrap();
    let rebuilt = builder.build().unwrap();
    let fresh = water_builder("UNK").0.build().unwrap();
    assert_eq!(
        classes(&rebuilt),
        (MoleculeClass::SmallMolecule, ResidueClass::Other)
    );
    assert!(rebuilt.same_layout(&fresh));
}

#[test]
fn explicit_classes_survive_changed_component_evidence_and_rebuilding() {
    let (mut builder, residue) = water_builder("HOH");
    let definition = builder.definitions().next().unwrap().0;
    // An explicit assignment equal to the inferred value still expresses intent.
    builder
        .set_molecule_class(definition, MoleculeClass::Water)
        .unwrap();
    builder
        .set_residue_class(residue, ResidueClass::Water)
        .unwrap();
    let source = Arc::new(builder.build().unwrap());
    let inferred = water_builder("HOH").0.build().unwrap();
    assert!(source.same_layout(&inferred));
    let mut editor = source.edit();
    let residue_handle = editor.residue_handle(residue).unwrap();
    editor
        .set_residue_component_ids(residue_handle, Some("UNK".into()), None)
        .unwrap();
    assert_eq!(
        classes(&editor.finish().unwrap()),
        (MoleculeClass::Water, ResidueClass::Water)
    );
    let mut builder = Arc::try_unwrap(source).unwrap().into_builder();
    builder
        .hierarchy_mut()
        .set_residue_component_ids(residue, Some("UNK".into()), None)
        .unwrap();
    let rebuilt = builder.build().unwrap();
    assert_eq!(
        classes(&rebuilt),
        (MoleculeClass::Water, ResidueClass::Water)
    );
    assert_eq!(
        classes(&rebuilt.into_builder().build().unwrap()),
        (MoleculeClass::Water, ResidueClass::Water)
    );
}

#[test]
fn no_op_and_append_preserve_untouched_classes_and_explicit_intent() {
    let source = Arc::new(water_builder("HOH").0.build().unwrap());
    assert!(Arc::ptr_eq(&source, &source.edit().finish().unwrap()));
    let mut editor = source.edit();
    editor.add_molecule(&oxygen()).unwrap();
    let appended = editor.finish().unwrap();
    assert_eq!(classes(&appended), classes(&source));
    assert_eq!(
        appended.molecules().nth(1).unwrap().class(),
        MoleculeClass::SmallMolecule
    );
    let mut builder = Arc::try_unwrap(source).unwrap().into_builder();
    builder.add_molecule(&oxygen()).unwrap();
    let appended_builder = builder.build().unwrap();
    assert!(appended.same_layout(&appended_builder));
}

#[test]
fn hierarchy_edits_reinfer_source_molecule_classes_and_preserve_shared_snapshot() {
    let source = Arc::new(Topology::from_molecule(&oxygen()).unwrap());
    let mut editor = source.edit();
    let atom = editor.atom_ids().next().unwrap();
    let chain = editor.add_chain("A", None).unwrap();
    let residue = editor.add_residue(chain, "HOH", None, None, None).unwrap();
    editor
        .add_atom_site(residue, atom, AtomSiteMetadata::default())
        .unwrap();
    let result = editor.finish().unwrap();
    assert_eq!(
        classes(&result),
        (MoleculeClass::Water, ResidueClass::Water)
    );
    assert_eq!(
        source.molecules().next().unwrap().class(),
        MoleculeClass::SmallMolecule
    );
    assert!(source.hierarchy().is_empty());
    let mut editor = result.edit();
    let residue = editor.residues().next().unwrap().0;
    editor
        .set_residue_component_ids(residue, Some("UNK".into()), None)
        .unwrap();
    assert!(editor
        .finish()
        .unwrap()
        .same_layout(&water_builder("UNK").0.build().unwrap()));
}

#[test]
fn transformations_do_not_promote_inferred_classes_to_explicit_overrides() {
    let source = Arc::new(water_builder("HOH").0.build().unwrap());
    let subset = source.subset(&AtomSelection::all(&source)).unwrap();
    let (target, _) = subset.into_parts();
    let mut editor = target.edit();
    let residue = editor.residues().next().unwrap().0;
    editor
        .set_residue_component_ids(residue, Some("UNK".into()), None)
        .unwrap();
    assert_eq!(
        classes(&editor.finish().unwrap()),
        (MoleculeClass::SmallMolecule, ResidueClass::Other)
    );
}

#[test]
fn preserved_inferred_class_survives_no_op_rebuilding_until_its_evidence_changes() {
    let (mut builder, residue) = water_builder("HOH");
    let definition = builder.definitions().next().unwrap().0;
    let first = builder.instances().next().unwrap().0;
    let second = builder.add_instance(definition).unwrap();
    let chain = builder.hierarchy_mut().add_chain("B", None).unwrap();
    let ion = builder
        .hierarchy_mut()
        .add_residue(chain, "UNK", None, None, None)
        .unwrap();
    let local = builder
        .definition(definition)
        .unwrap()
        .molecule()
        .atom_ids()
        .next()
        .unwrap();
    builder
        .hierarchy_mut()
        .add_atom_site(
            ion,
            InstanceAtomId::new(second, local),
            AtomSiteMetadata::default(),
        )
        .unwrap();
    builder.set_residue_class(ion, ResidueClass::Ion).unwrap();
    let source = Arc::new(builder.build().unwrap());
    assert_eq!(
        source.molecules().next().unwrap().class(),
        MoleculeClass::Other
    );
    let retained = kekule::topology::transform::retain_instances(&source, [first]).unwrap();
    assert_eq!(
        classes(&retained),
        (MoleculeClass::Other, ResidueClass::Water)
    );
    let mut editor = retained.edit();
    editor
        .insert_property(PropertyKey::new("label").unwrap(), PropertyValue::Int(1))
        .unwrap();
    assert_eq!(classes(&editor.finish().unwrap()), classes(&retained));
    let mut editor = retained.edit();
    let site = editor.atom_sites().next().unwrap().0;
    let mut metadata = editor.atom_site(site).unwrap().metadata().clone();
    metadata.type_symbol = Some("O".into());
    editor.set_atom_site_metadata(site, metadata).unwrap();
    let residue_handle = editor.residues().next().unwrap().0;
    editor
        .set_residue_class(residue_handle, ResidueClass::Water)
        .unwrap();
    assert_eq!(
        editor.residue(residue_handle).unwrap().class_override(),
        Some(ResidueClass::Water)
    );
    assert_eq!(classes(&editor.finish().unwrap()), classes(&retained));
    let rebuilt = Arc::try_unwrap(retained)
        .unwrap()
        .into_builder()
        .build()
        .unwrap();
    assert_eq!(
        classes(&rebuilt),
        (MoleculeClass::Other, ResidueClass::Water)
    );
    let mut builder = rebuilt.into_builder();
    builder.add_molecule(&oxygen()).unwrap();
    let appended = builder.build().unwrap();
    assert_eq!(
        classes(&appended),
        (MoleculeClass::Other, ResidueClass::Water)
    );
    let mut builder = appended.into_builder();
    builder
        .hierarchy_mut()
        .set_residue_component_ids(residue, Some("UNK".into()), None)
        .unwrap();
    assert_eq!(
        classes(&builder.build().unwrap()),
        (MoleculeClass::SmallMolecule, ResidueClass::Other)
    );
}

#[test]
fn consuming_publication_moves_unique_source_definitions_and_modified_drafts() {
    let molecule = oxygen();
    let local = molecule.atom_ids().next().unwrap();
    let source = Topology::from_molecule(&molecule).unwrap();
    let original_pointer = source
        .molecules()
        .next()
        .unwrap()
        .molecule()
        .atom(local)
        .unwrap() as *const Atom;
    let mut editor = source.into_editor();
    editor.add_atom(atom("C")).unwrap();
    let result = editor.finish().unwrap();
    assert_eq!(
        result
            .molecules()
            .next()
            .unwrap()
            .molecule()
            .atom(local)
            .unwrap() as *const Atom,
        original_pointer
    );

    let mut editor = result.edit();
    let handle = editor.atom_ids().next().unwrap();
    editor.replace_atom(handle, atom("N")).unwrap();
    let draft_pointer = editor.atom(handle).unwrap() as *const Atom;
    editor.validate().unwrap();
    assert_eq!(editor.atom(handle).unwrap() as *const Atom, draft_pointer);
    let edited = editor.finish().unwrap();
    assert_eq!(
        edited
            .molecules()
            .next()
            .unwrap()
            .molecule()
            .atom(local)
            .unwrap() as *const Atom,
        draft_pointer
    );
    assert_eq!(
        result
            .molecules()
            .next()
            .unwrap()
            .molecule()
            .atom(local)
            .unwrap()
            .element,
        Element::from_symbol("O").unwrap()
    );
}

#[test]
fn consuming_publication_moves_added_molecules_and_preserves_owner_annotations() {
    let mut editor = TopologyEditor::new();
    editor.add_molecule(&oxygen()).unwrap();
    let atom = editor.atom_ids().next().unwrap();
    let pointer = editor.atom(atom).unwrap() as *const Atom;
    let key = PropertyKey::new("label").unwrap();
    editor
        .insert_property(key.clone(), PropertyValue::String("draft".into()))
        .unwrap();
    let topology = editor.finish().unwrap();
    assert_eq!(topology.atoms().next().unwrap().1 as *const Atom, pointer);
    assert_eq!(
        topology.properties().get(&key),
        Some(&PropertyValue::String("draft".into()))
    );
}

#[test]
fn topology_transform_errors_expose_wrapped_causes() {
    let error = TopologyTransformError::TopologyBuild(TopologyBuildError::NoMoleculeInstances);
    assert!(error.source().unwrap().is::<TopologyBuildError>());
    let error = TopologySubsetError::Selection(SelectionError::TopologyMismatch);
    assert!(error.source().unwrap().is::<SelectionError>());
    let error = TopologySubsetError::Publication(MoleculePublicationError::EmptyGraph);
    assert!(error.source().unwrap().is::<MoleculePublicationError>());
    assert!(TopologySubsetError::EmptySelection.source().is_none());
    assert!(TopologyTransformError::EmptyTargetTopology
        .source()
        .is_none());
}

#[test]
fn editor_hierarchy_membership_stays_aligned_during_bulk_assembly_and_site_edits() {
    let mut editor = TopologyEditor::new();
    let chain = editor.add_chain("A", None).unwrap();
    let mut residues = Vec::new();
    let mut sites = Vec::new();
    for _ in 0..512 {
        let atom = editor.add_atom(atom("O")).unwrap();
        let residue = editor.add_residue(chain, "HOH", None, None, None).unwrap();
        sites.push(
            editor
                .add_atom_site(residue, atom, AtomSiteMetadata::default())
                .unwrap(),
        );
        residues.push(residue);
    }
    editor.set_atom_site_residue(sites[0], residues[1]).unwrap();
    assert!(editor.residue(residues[0]).is_err());
    editor.delete_atom_site(sites[1]).unwrap();
    assert!(editor.residue(residues[1]).is_ok());
    for &site in sites.iter().skip(2).step_by(2) {
        editor.delete_atom_site(site).unwrap();
    }
    editor.validate().unwrap();
    let mut cleared = editor.clone();
    let topology = editor.finish().unwrap();
    assert_eq!(topology.atom_count(), 512);
    assert_eq!(topology.hierarchy().atom_sites().count(), 256);
    assert_eq!(topology.residues().count(), 256);
    assert_eq!(
        topology
            .molecules()
            .filter(|m| m.class() == MoleculeClass::Water)
            .count(),
        256
    );
    let mut imported = topology.edit();
    let site = imported.atom_sites().next().unwrap().0;
    imported.delete_atom_site(site).unwrap();
    assert_eq!(
        imported.finish().unwrap().hierarchy().atom_sites().count(),
        255
    );
    cleared.clear();
    cleared.add_atom(atom("C")).unwrap();
    assert!(cleared.finish().unwrap().hierarchy().is_empty());
}
