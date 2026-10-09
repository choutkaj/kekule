//! Represented molecular identity: what equality compares and what perception
//! may add without changing it.

use kekule::core::{Molecule, Perception};
use kekule::smiles;
use kekule::topology::Topology;

fn one_smiles(input: &str) -> Molecule {
    let mut molecules = smiles::to_molecules(input).expect("SMILES interprets");
    assert_eq!(molecules.len(), 1);
    molecules.pop().expect("component count was checked")
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
    assert!(Topology::from_molecule(molecule.clone())
        .unwrap()
        .same_layout(&Topology::from_molecule(rewired.clone()).unwrap()));
}
