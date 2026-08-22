use std::sync::Arc;

use kekule::core::{Atom, Element, MoleculeEditor};
use kekule::topology::transform::{remove_instances, retain_instances, TopologyTransformError};
use kekule::topology::{MoleculeInstanceId, Topology, TopologyBuilder};

fn repeated_topology() -> Arc<Topology> {
    let mut editor = MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    Arc::new(builder.build().unwrap())
}

#[test]
fn retain_and_remove_return_topology_directly_in_source_order() {
    let source = repeated_topology();
    let ids = source.instances().map(|(id, _)| id).collect::<Vec<_>>();

    let retained = retain_instances(&source, [ids[2], ids[0]]).unwrap();
    assert_eq!(retained.instance_count(), 2);
    assert_eq!(retained.definition_count(), 1);
    assert_eq!(retained.atom_count(), 2);
    assert!(!Arc::ptr_eq(&source, &retained));

    let removed = remove_instances(&source, [ids[1]]).unwrap();
    assert_eq!(removed.instance_count(), 2);
    assert_eq!(removed.definition_count(), 1);
}

#[test]
fn no_op_subset_operations_preserve_the_original_arc() {
    let source = repeated_topology();
    let ids = source.instances().map(|(id, _)| id).collect::<Vec<_>>();
    let retained = retain_instances(&source, ids).unwrap();
    let removed = remove_instances(&source, std::iter::empty()).unwrap();
    assert!(Arc::ptr_eq(&retained, &source));
    assert!(Arc::ptr_eq(&removed, &source));
}

#[test]
fn subset_operations_reject_invalid_or_empty_targets() {
    let source = repeated_topology();
    assert_eq!(
        retain_instances(&source, std::iter::empty()).unwrap_err(),
        TopologyTransformError::EmptyTargetTopology
    );
    assert_eq!(
        remove_instances(&source, source.instances().map(|(id, _)| id)).unwrap_err(),
        TopologyTransformError::EmptyTargetTopology
    );
    assert_eq!(
        retain_instances(&source, [MoleculeInstanceId::new(99)]).unwrap_err(),
        TopologyTransformError::InvalidSourceInstance(MoleculeInstanceId::new(99))
    );
}
