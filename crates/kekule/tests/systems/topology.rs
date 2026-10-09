//! Topology publication: consuming builders and editors, and transform errors.

use std::error::Error;

use kekule::core::{Atom, Element, Molecule, MoleculeEditor, MoleculePublicationError};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::topology::transform::{TopologySubsetError, TopologyTransformError};
use kekule::topology::{SelectionError, Topology, TopologyBuildError, TopologyEditor};

fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}

fn oxygen() -> Molecule {
    let mut editor = MoleculeEditor::new();
    editor.add_atom(atom("O")).unwrap();
    editor.finish().unwrap()
}

#[test]
fn consuming_publication_moves_unique_source_definitions_and_modified_drafts() {
    let molecule = oxygen();
    let local = molecule.atom_ids().next().unwrap();
    let source = Topology::from_molecule(molecule.clone()).unwrap();
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
    editor.add_molecule(oxygen().clone()).unwrap();
    let atom = editor.atom_ids().next().unwrap();
    let pointer = editor.atom(atom).unwrap() as *const Atom;
    let key = PropertyKey::new("label").unwrap();
    editor
        .insert_property(key.clone(), PropertyValue::String("draft".into()))
        .unwrap();
    let topology = editor.finish().unwrap();
    assert_eq!(
        topology.atoms().next().unwrap().atom() as *const Atom,
        pointer
    );
    assert_eq!(
        topology.properties().owner().get(&key),
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
