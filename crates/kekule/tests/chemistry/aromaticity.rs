use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor, Perception};
use kekule::perception::aromaticity::{
    perceive_aromaticity, perceive_aromaticity_with_options, AromaticityError, AromaticityModel,
    AromaticityOptions,
};
use kekule::perception::rings::RingPerceptionOptions;
use kekule::smiles;
use kekule::stereo::assign_cip_descriptors;

use crate::support::aromaticity_cases;

fn molecule(source: &str) -> Molecule {
    smiles::to_molecules(source).unwrap().pop().unwrap()
}

#[test]
fn default_perception_matches_reference_aromaticity_regressions() {
    let mut failures = Vec::new();
    for case in aromaticity_cases() {
        let mut molecule = case.molecule();
        molecule
            .perceive()
            .unwrap_or_else(|error| panic!("{}: {error}", case.label));
        let aromatic = |atom| molecule.atom_is_aromatic(atom).unwrap() == Some(true);
        let atoms = molecule
            .atom_ids()
            .filter(|&atom| aromatic(atom))
            .map(|atom| atom.index())
            .collect::<Vec<_>>();
        let mut nonaromatic_bonds = Vec::new();
        for (id, bond) in molecule.bonds() {
            let endpoints = (bond.a().index(), bond.b().index());
            let endpoints = (endpoints.0.min(endpoints.1), endpoints.0.max(endpoints.1));
            let both_aromatic = aromatic(bond.a()) && aromatic(bond.b());
            match molecule.bond_is_aromatic(id).unwrap() == Some(true) {
                true if !both_aromatic => failures.push(format!(
                    "{}: aromatic bond {endpoints:?} has a non-aromatic endpoint",
                    case.label
                )),
                false if both_aromatic => nonaromatic_bonds.push(endpoints),
                _ => {}
            }
        }
        nonaromatic_bonds.sort_unstable();
        if atoms != case.aromatic_atoms || nonaromatic_bonds != case.nonaromatic_bonds {
            failures.push(format!(
                "{}: aromatic atoms {atoms:?}, non-aromatic bonds {nonaromatic_bonds:?}; \
                 expected {:?} and {:?}",
                case.label, case.aromatic_atoms, case.nonaromatic_bonds
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

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
    assert_eq!(complete.perception().inferred_hydrogens(sulfur), Some(1));
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

#[test]
fn aromaticity_work_failures_restore_both_missing_and_installed_perception() {
    let mut perceived = molecule("F[C@H](Cl)c1ccccc1-c1ccccc1");
    perceived.perceive().unwrap();
    assign_cip_descriptors(&mut perceived).unwrap();
    for original in [molecule("F[C@H](Cl)c1ccccc1-c1ccccc1"), perceived] {
        let mut failures = 0;
        let mut successes = 0;
        // Exercise failures at successive stages, including after rings have
        // been installed and after aromaticity has replaced the old section.
        for max_total_work in (0..1_000).step_by(7) {
            let mut actual = original.clone();
            let options = AromaticityOptions {
                max_total_work,
                ..Default::default()
            };
            match perceive_aromaticity_with_options(
                &mut actual,
                AromaticityModel::RdkitLike,
                options,
            ) {
                Err(error) => {
                    assert_eq!(
                        error,
                        AromaticityError::ResourceLimit {
                            limit: max_total_work
                        }
                    );
                    assert_eq!(actual.perception(), original.perception());
                    failures += 1;
                }
                Ok(()) => {
                    assert_eq!(
                        actual
                            .atoms()
                            .filter(|(id, _)| actual.atom_is_aromatic(*id).unwrap() == Some(true))
                            .count(),
                        12
                    );
                    assert_eq!(
                        actual
                            .bonds()
                            .filter(|(id, _)| actual.bond_is_aromatic(*id).unwrap() == Some(true))
                            .count(),
                        12
                    );
                    successes += 1;
                }
            }
            assert_eq!(actual, original, "represented chemistry changed");
        }
        assert!(failures > 0 && successes > 0);
    }
}

#[test]
fn installed_rings_skip_enumeration_limits_but_not_aromaticity_limits() {
    let mut actual = molecule("c1ccccc1");
    actual.perceive().unwrap();
    let options = AromaticityOptions {
        ring_options: RingPerceptionOptions {
            max_atoms: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    perceive_aromaticity_with_options(&mut actual, AromaticityModel::RdkitLike, options).unwrap();
    let before = actual.perception().clone();
    assert!(matches!(
        perceive_aromaticity_with_options(
            &mut actual,
            AromaticityModel::RdkitLike,
            AromaticityOptions {
                max_total_work: 0,
                ..options
            }
        ),
        Err(AromaticityError::ResourceLimit { limit: 0 })
    ));
    assert_eq!(actual.perception(), &before);
}
