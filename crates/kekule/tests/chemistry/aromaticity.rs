use kekule::diagnostics::Resource;
use kekule::perception::aromaticity::{
    perceive_aromaticity, AromaticityError, AromaticityModel, AromaticityOptions,
};
use kekule::perception::rings::RingPerceptionOptions;
use kekule::perception::valence::{ValenceModel, ValenceValidation};
use kekule::smiles;
use kekule::stereo::cip::assign_cip_descriptors;
use kekule::{
    molecule::{Atom, BondOrder, Element, Molecule, MoleculeEditor},
    perception::Perception,
};

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
        let aromatic = |atom| {
            molecule
                .perception()
                .aromaticity()
                .map(|aromaticity| aromaticity.atom_is_aromatic(atom))
                == Some(true)
        };
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
            match molecule
                .perception()
                .aromaticity()
                .map(|aromaticity| aromaticity.bond_is_aromatic(id))
                == Some(true)
            {
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
                atom = atom.with_formal_charge(1);
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
    assert_eq!(
        complete
            .perception()
            .valence()
            .and_then(|valence| valence.inferred_hydrogen_count(sulfur)),
        Some(1)
    );
    assert!(complete.atoms().all(|(atom, _)| complete
        .perception()
        .aromaticity()
        .map(|aromaticity| aromaticity.atom_is_aromatic(atom))
        == Some(true)));
    assert!(complete.bonds().all(|(bond, _)| complete
        .perception()
        .aromaticity()
        .map(|aromaticity| aromaticity.bond_is_aromatic(bond))
        == Some(true)));

    perceive_aromaticity(
        &mut standalone,
        AromaticityModel::RdkitLike,
        Default::default(),
    )
    .unwrap();

    assert_eq!(
        standalone.perception().aromaticity(),
        complete.perception().aromaticity()
    );
    assert!(standalone.perception().valence().is_none());
    assert_eq!(standalone, represented);
}

#[test]
fn standalone_aromaticity_preserves_partial_installed_hydrogen_assignments() {
    let mut molecule = protonated_thiophene_with_inferred_hydrogen();
    let sulfur = molecule.atom_ids().next().unwrap();
    let installed = Perception::builder()
        .with_valence(
            ValenceModel::RdkitLike,
            ValenceValidation::Permissive,
            vec![(sulfur, 0)],
        )
        .unwrap()
        .build();
    molecule.install_perception(installed.clone()).unwrap();

    perceive_aromaticity(
        &mut molecule,
        AromaticityModel::RdkitLike,
        Default::default(),
    )
    .unwrap();

    // An expert-installed zero takes precedence over the model's inferred one.
    assert_eq!(molecule.perception().valence(), installed.valence());
    assert!(molecule.atoms().all(|(atom, _)| molecule
        .perception()
        .aromaticity()
        .map(|aromaticity| aromaticity.atom_is_aromatic(atom))
        == Some(false)));
    assert!(molecule.bonds().all(|(bond, _)| molecule
        .perception()
        .aromaticity()
        .map(|aromaticity| aromaticity.bond_is_aromatic(bond))
        == Some(false)));
}

#[test]
fn aromaticity_work_failures_restore_both_missing_and_installed_perception() {
    let mut perceived = molecule("F[C@H](Cl)c1ccccc1-c1ccccc1");
    perceived.perceive().unwrap();
    assign_cip_descriptors(&mut perceived, Default::default()).unwrap();
    for original in [molecule("F[C@H](Cl)c1ccccc1-c1ccccc1"), perceived] {
        let mut failures = 0;
        let mut successes = 0;
        // Exercise failures at successive stages, including after rings have
        // been installed and after aromaticity has replaced the old section.
        for max_total_work in (0..1_000).step_by(7) {
            let mut actual = original.clone();
            let options = AromaticityOptions::default().with_max_total_work(max_total_work);
            match perceive_aromaticity(&mut actual, AromaticityModel::RdkitLike, options) {
                Err(error) => {
                    assert!(matches!(
                        error,
                        AromaticityError::LimitExceeded(limit)
                            if limit.resource() == Resource::Work
                                && limit.limit() == max_total_work
                                && limit.observed() > max_total_work
                    ));
                    assert_eq!(actual.perception(), original.perception());
                    failures += 1;
                }
                Ok(()) => {
                    assert_eq!(
                        actual
                            .atoms()
                            .filter(|(id, _)| actual
                                .perception()
                                .aromaticity()
                                .map(|aromaticity| aromaticity.atom_is_aromatic(*id))
                                == Some(true))
                            .count(),
                        12
                    );
                    assert_eq!(
                        actual
                            .bonds()
                            .filter(|(id, _)| actual
                                .perception()
                                .aromaticity()
                                .map(|aromaticity| aromaticity.bond_is_aromatic(*id))
                                == Some(true))
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
    let options = AromaticityOptions::default()
        .with_ring_options(RingPerceptionOptions::default().with_max_atoms(0));
    perceive_aromaticity(&mut actual, AromaticityModel::RdkitLike, options).unwrap();
    let before = actual.perception().clone();
    assert!(matches!(
            perceive_aromaticity(
                &mut actual,
                AromaticityModel::RdkitLike,
                options
    .with_max_total_work(0)
            ),
            Err(AromaticityError::LimitExceeded(limit))
                if limit.resource() == Resource::Work && limit.limit() == 0
        ));
    assert_eq!(actual.perception(), &before);
}

/// Aromatic atoms and aromatic bonds as sorted source indices.
fn aromatic_masks(molecule: &Molecule) -> (Vec<usize>, Vec<(usize, usize)>) {
    let aromaticity = molecule.perception().aromaticity().unwrap();
    let atoms = molecule
        .atom_ids()
        .filter(|&atom| aromaticity.atom_is_aromatic(atom))
        .map(|atom| atom.index())
        .collect();
    let mut bonds = molecule
        .bonds()
        .filter(|&(id, _)| aromaticity.bond_is_aromatic(id))
        .map(|(_, bond)| {
            let (a, b) = (bond.a().index(), bond.b().index());
            (a.min(b), a.max(b))
        })
        .collect::<Vec<_>>();
    bonds.sort_unstable();
    (atoms, bonds)
}

/// The ring bonds of `molecule` between two of `atoms`, minus `except`.
fn ring_bonds_among(
    molecule: &Molecule,
    atoms: &[usize],
    except: &[(usize, usize)],
) -> Vec<(usize, usize)> {
    let membership = molecule.perception().rings().unwrap().membership();
    let mut bonds = molecule
        .bonds()
        .filter(|&(id, _)| membership.bond_in_ring(id))
        .map(|(_, bond)| {
            let (a, b) = (bond.a().index(), bond.b().index());
            (a.min(b), a.max(b))
        })
        .filter(|(a, b)| atoms.contains(a) && atoms.contains(b) && !except.contains(&(*a, *b)))
        .collect::<Vec<_>>();
    bonds.sort_unstable();
    bonds
}

/// Label, SMILES, aromatic atoms, and non-aromatic ring bonds between them.
type FusedCase = (
    &'static str,
    &'static str,
    &'static [usize],
    &'static [(usize, usize)],
);

#[test]
fn fused_aromatic_systems_are_aromatic_whole() {
    // Each case lists its aromatic atoms and the ring bonds between two of
    // them that stay non-aromatic; every other ring bond between aromatic
    // atoms is aromatic, and no other bond is.
    let cases: &[FusedCase] = &[
        // Ten-electron perimeters whose rings are not aromatic alone: the
        // fusion bond is aromatic with the system (RDKit leaves it out).
        (
            "azulene",
            "C1=CC=C2C=CC=C2C=C1",
            &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            &[],
        ),
        (
            "cyclopenta[c]pyran",
            "C1=COC=C2C1=CC=C2",
            &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            &[],
        ),
        (
            "pyrrolopyrimidine tautomer",
            "C1=C2C=CN=C2NC=N1",
            &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            &[],
        ),
        // Isoalloxazine: the pyrazine and pyrimidinedione rings are aromatic
        // only together, or with the benzene ring.
        (
            "lumiflavin",
            "CC1=CC2=C(C=C1C)N(C3=NC(=O)NC(=O)C3=N2)C",
            &[1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 16, 17],
            &[],
        ),
        // Rings aromatic alone, as RDKit marks them.
        (
            "naphthalene",
            "C1=CC=C2C=CC=CC2=C1",
            &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            &[],
        ),
        (
            "indole",
            "C1=CC=C2C(=C1)C=CN2",
            &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            &[],
        ),
        (
            "purine",
            "C1=C2C(=NC=N1)N=CN2",
            &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            &[],
        ),
        (
            "guanine",
            "C1=NC2=C(N1)C(=O)NC(=N2)N",
            &[0, 1, 2, 3, 4, 5, 7, 8, 9],
            &[],
        ),
        (
            "caffeine",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            &[1, 2, 3, 4, 5, 6, 8, 9, 11],
            &[],
        ),
        // A bond outside every aromatic ring stays non-aromatic, between
        // aromatic atoms too.
        (
            "biphenylene",
            "C1=CC=C2C(=C1)C3=CC=CC=C32",
            &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
            &[(3, 11), (4, 6)],
        ),
        (
            "fluorene",
            "C1=CC=C2C(=C1)CC3=CC=CC=C32",
            &[0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12],
            &[(3, 12)],
        ),
        // A fusion bond of an aromatic and a non-aromatic ring is aromatic
        // with the aromatic ring alone; the other ring's bonds are not.
        ("indane", "C1=CC=C2CCCC2=C1", &[0, 1, 2, 3, 7, 8], &[]),
        (
            "1,4-naphthoquinone",
            "C1=CC=C2C(=C1)C(=O)C=CC2=O",
            &[0, 1, 2, 3, 4, 5],
            &[],
        ),
        // An eight-electron perimeter is not aromatic.
        ("pentalene", "C1=CC2=CC=CC2=C1", &[], &[]),
    ];
    for &(label, source, atoms, except) in cases {
        let mut molecule = molecule(source);
        molecule.perceive().unwrap();
        let expected = (atoms.to_vec(), ring_bonds_among(&molecule, atoms, except));
        assert_eq!(aromatic_masks(&molecule), expected, "{label}");
    }

    // The MDL model marks only perimeters, as RDKit's does.
    let mut azulene = molecule("C1=CC=C2C=CC=C2C=C1");
    azulene.perceive().unwrap();
    perceive_aromaticity(&mut azulene, AromaticityModel::Mdl, Default::default()).unwrap();
    let (atoms, bonds) = aromatic_masks(&azulene);
    assert_eq!(atoms, (0..10).collect::<Vec<_>>());
    assert_eq!(bonds, ring_bonds_among(&azulene, &atoms, &[(3, 7)]));
}
