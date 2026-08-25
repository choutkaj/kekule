use kekule::core::{
    Atom, AtomId, BondOrder, Element, Molecule, MoleculeEditor, MoleculePublicationError,
    Perception,
};
use kekule::smiles;
use kekule::topology::{InstanceAtomId, TopologyBuilder};

fn carbon() -> Atom {
    Atom::new(Element::from_symbol("C").expect("carbon exists"))
}

fn one_smiles(input: &str) -> Molecule {
    let mut molecules = smiles::to_molecules(input).expect("SMILES interprets");
    assert_eq!(molecules.len(), 1);
    molecules.pop().expect("component count was checked")
}

#[test]
fn editor_preserves_stable_local_ids_and_publishes_only_connected_graphs() {
    let mut editor = MoleculeEditor::new();
    let first = editor.add_atom(carbon()).unwrap();
    let tombstone = editor.add_atom(carbon()).unwrap();
    let last = editor.add_atom(carbon()).unwrap();
    editor
        .add_bond(first, tombstone, BondOrder::Single)
        .unwrap();
    editor.add_bond(tombstone, last, BondOrder::Single).unwrap();
    editor.delete_atom(tombstone).unwrap();
    assert!(matches!(
        editor.clone().finish(),
        Err(MoleculePublicationError::DisconnectedGraph(_))
    ));
    editor.add_bond(first, last, BondOrder::Single).unwrap();
    let molecule = editor.finish().expect("connected graph publishes");
    assert_eq!(molecule.atom_ids().collect::<Vec<_>>(), [first, last]);
    assert_eq!(last, AtomId::new(2));
}

#[test]
fn perception_is_reconstructible_and_not_part_of_represented_equality() {
    let represented = one_smiles("c1ccccc1");
    let mut perceived = represented.clone();
    perceived.perceive().expect("benzene perceives");
    assert_ne!(perceived.perception(), &Perception::default());
    assert_eq!(perceived, represented);
    let exported = perceived.perception().clone();
    perceived.clear_perception();
    perceived
        .install_perception(exported.clone())
        .expect("matching perception reinstalls");
    assert_eq!(perceived.perception(), &exported);
}

#[test]
fn hierarchy_is_owned_by_topology() {
    let mut editor = MoleculeEditor::new();
    let atom = editor.add_atom(carbon()).unwrap();
    let molecule = editor
        .finish()
        .expect("molecule publishes without hierarchy");
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    builder
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, atom),
            Default::default(),
        )
        .unwrap();
    let topology = builder.build().expect("valid hierarchy publishes");
    assert_eq!(topology.hierarchy().chains().count(), 1);
    assert_eq!(topology.hierarchy().residues().count(), 1);
    assert_eq!(topology.hierarchy().atom_sites().count(), 1);
}

#[test]
fn topology_definitions_consume_the_universal_molecule_directly() {
    let molecule = one_smiles("CCO");
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(definition).unwrap();
    let topology = builder.build().unwrap();
    let stored = topology.definition(definition).unwrap().molecule();
    assert_eq!(stored, &molecule);
    assert!(topology.hierarchy().is_empty());
}
