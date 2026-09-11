use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor, Perception};
use kekule::perception::aromaticity::{perceive_aromaticity, AromaticityModel};

fn protonated_thiophene_with_inferred_hydrogen() -> Molecule {
    let mut editor = MoleculeEditor::new();
    let atoms = ["S", "C", "C", "C", "C"]
        .into_iter()
        .enumerate()
        .map(|(index, symbol)| {
            let mut atom = Atom::new(Element::from_symbol(symbol).unwrap());
            if index == 0 {
                atom.formal_charge = 1;
            }
            editor.add_atom(atom).unwrap()
        })
        .collect::<Vec<_>>();
    for (index, order) in [
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
    ]
    .into_iter()
    .enumerate()
    {
        editor
            .add_bond(atoms[index], atoms[(index + 1) % atoms.len()], order)
            .unwrap();
    }
    editor.finish().unwrap()
}

#[test]
fn standalone_aromaticity_uses_the_same_charged_hydrogen_rules_as_default_perception() {
    // RDKit 2026.03.3 perceives [SH+]1C=CC=C1 as five aromatic atoms and bonds.
    // This equivalent represented graph permits the sulfur hydrogen to be inferred.
    let mut standalone = protonated_thiophene_with_inferred_hydrogen();
    let represented = standalone.clone();
    let sulfur = standalone.atom_ids().next().unwrap();
    let mut complete = standalone.clone();
    complete.perceive().unwrap();
    assert_eq!(complete.implicit_hydrogens(sulfur).unwrap(), Some(1));
    assert!(complete
        .atoms()
        .all(|(atom, _)| complete.atom_is_aromatic(atom).unwrap() == Some(true)));
    assert!(complete
        .bonds()
        .all(|(bond, _)| complete.bond_is_aromatic(bond).unwrap() == Some(true)));

    perceive_aromaticity(&mut standalone, AromaticityModel::RdkitLike).unwrap();

    assert_eq!(
        standalone.perception().aromaticity_state(),
        complete.perception().aromaticity_state()
    );
    assert!(!standalone.perception().has_valence());
    assert_eq!(standalone, represented);
}

#[test]
fn standalone_aromaticity_preserves_partial_installed_hydrogen_assignments() {
    let mut molecule = protonated_thiophene_with_inferred_hydrogen();
    let sulfur = molecule.atom_ids().next().unwrap();
    let installed = Perception::builder()
        .with_valence(None, vec![(sulfur, 0)])
        .unwrap()
        .build();
    molecule.install_perception(installed.clone()).unwrap();

    perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike).unwrap();

    // An expert-installed zero takes precedence over the model's inferred one.
    assert_eq!(
        molecule.perception().valence_state(),
        installed.valence_state()
    );
    assert!(molecule
        .atoms()
        .all(|(atom, _)| molecule.atom_is_aromatic(atom).unwrap() == Some(false)));
    assert!(molecule
        .bonds()
        .all(|(bond, _)| molecule.bond_is_aromatic(bond).unwrap() == Some(false)));
}
