//! Molfile atoms whose drawing leaves their charge implicit: onium centres
//! and unbonded metals.

use kekule::molecule::Molecule;
use kekule::{molfile, smiles};

/// A V2000 record of `atoms` (symbol and charge) at spaced positions, bonded
/// by single `bonds` between one-based atom numbers.
fn v2000(atoms: &[(&str, i8)], bonds: &[(usize, usize)]) -> String {
    let mut text = format!(
        "drawn\nkekule\n\n{:>3}{:>3}  0  0  0  0            999 V2000\n",
        atoms.len(),
        bonds.len()
    );
    for (index, (symbol, charge)) in atoms.iter().enumerate() {
        let code = match charge {
            1 => 3,
            -1 => 5,
            _ => 0,
        };
        text.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {symbol:<3} 0{code:>3}  0  0  0  0  0  0  0  0  0  0\n",
            index as f64 * 1.5,
            0.0,
            0.0
        ));
    }
    for (a, b) in bonds {
        text.push_str(&format!("{a:>3}{b:>3}  1  0  0  0  0\n"));
    }
    text.push_str("M  END\n");
    text
}

fn interpret(text: &str) -> molfile::MolfileInterpretation {
    molfile::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// Each atom's symbol, formal charge, and total hydrogens after perception.
fn perceived(molecule: &Molecule) -> Vec<(String, i8, usize)> {
    let mut molecule = molecule.clone();
    molecule.perceive().unwrap();
    molecule
        .atoms()
        .map(|(id, atom)| {
            (
                atom.element().symbol().to_owned(),
                atom.formal_charge(),
                molecule.total_hydrogens(id).unwrap(),
            )
        })
        .collect()
}

#[test]
fn neutral_onium_centres_take_their_charge_not_a_hydrogen() {
    for (center, ligands) in [
        ("S", 3),
        ("Se", 3),
        ("Te", 3),
        ("P", 4),
        ("As", 4),
        ("I", 2),
    ] {
        let mut atoms = vec![(center, 0)];
        atoms.extend(std::iter::repeat_n(("C", 0), ligands));
        let bonds = (2..=ligands + 1)
            .map(|ligand| (1, ligand))
            .collect::<Vec<_>>();
        let interpretation = interpret(&v2000(&atoms, &bonds));
        let molecule = interpretation.molecules().next().unwrap();
        assert_eq!(
            perceived(molecule)[0],
            (center.to_owned(), 1, 0),
            "{center} with {ligands} carbons"
        );
        assert_eq!(
            interpretation.reports()[0].issues(),
            [molfile::MolfileInterpretIssue::ChargeFromBonds {
                atom: molecule.atom_ids().next().unwrap(),
                location: interpretation.reports()[0].atom_mappings()[0].location(),
                charge: 1,
            }],
            "{center}"
        );
    }
}

#[test]
fn stated_charges_and_default_valences_are_kept() {
    // Charged sulfonium, dimethyl sulfide, and a sulfoxide need nothing.
    for (atoms, bonds) in [
        (
            vec![("S", 1), ("C", 0), ("C", 0), ("C", 0)],
            vec![(1, 2), (1, 3), (1, 4)],
        ),
        (vec![("S", 0), ("C", 0), ("C", 0)], vec![(1, 2), (1, 3)]),
        (
            vec![("P", 0), ("C", 0), ("C", 0), ("C", 0)],
            vec![(1, 2), (1, 3), (1, 4)],
        ),
    ] {
        let expected = atoms.iter().map(|(_, charge)| *charge).collect::<Vec<_>>();
        let interpretation = interpret(&v2000(&atoms, &bonds));
        let molecule = interpretation.molecules().next().unwrap();
        let charges = perceived(molecule)
            .into_iter()
            .map(|(_, charge, _)| charge)
            .collect::<Vec<_>>();
        assert_eq!(charges, expected);
        assert!(interpretation.reports()[0].issues().is_empty());
    }
}

#[test]
fn smiles_onium_spellings_keep_their_opensmiles_hydrogen() {
    // OpenSMILES gives the organic-subset sulfur of CS(C)C one hydrogen.
    let molecule = smiles::to_molecules("CS(C)C").unwrap().remove(0);
    assert_eq!(perceived(&molecule)[1], ("S".to_owned(), 0, 1));
}

#[test]
fn molfiles_written_from_onium_hydrides_read_back_unchanged() {
    // The hydrogen of CS(C)C is not the Molfile default, so the writer states
    // the valence and the reader keeps the hydride neutral.
    let molecule = smiles::to_molecules("CS(C)C").unwrap().remove(0);
    for version in [
        molfile::MolfileWriteVersion::V2000,
        molfile::MolfileWriteVersion::V3000,
    ] {
        let text = molfile::write(
            &molecule,
            molfile::MolfileWriteOptions::default().with_version(version),
        )
        .unwrap()
        .0;
        let interpretation = interpret(&text);
        assert_eq!(
            perceived(interpretation.molecules().next().unwrap())[1],
            ("S".to_owned(), 0, 1),
            "{version:?}"
        );
        assert!(interpretation.reports()[0].issues().is_empty());
    }
}

#[test]
fn unbonded_metals_gain_no_hydrides_and_are_reported() {
    for symbol in ["Li", "Na", "K", "Mg", "Ca", "Ba", "Al", "Ga", "Sn", "Pb"] {
        let interpretation = interpret(&v2000(&[(symbol, 0)], &[]));
        let molecule = interpretation.molecules().next().unwrap();
        assert_eq!(perceived(molecule), [(symbol.to_owned(), 0, 0)], "{symbol}");
        assert_eq!(
            interpretation.reports()[0].issues(),
            [molfile::MolfileInterpretIssue::NeutralUnbondedMetal {
                atom: molecule.atom_ids().next().unwrap(),
                location: interpretation.reports()[0].atom_mappings()[0].location(),
            }],
            "{symbol}"
        );
    }
}

#[test]
fn under_bonded_and_charged_metals_gain_no_hydrides() {
    // A charged ion is complete; a metal with a bond is not reported.
    for (atoms, bonds) in [
        (vec![("Na", 1)], vec![]),
        (vec![("Mg", 0), ("C", 0)], vec![(1, 2)]),
        (vec![("Al", 0), ("C", 0)], vec![(1, 2)]),
    ] {
        let interpretation = interpret(&v2000(&atoms, &bonds));
        let molecule = interpretation.molecules().next().unwrap();
        assert_eq!(perceived(molecule)[0].2, 0, "{atoms:?}");
        assert!(interpretation.reports()[0].issues().is_empty(), "{atoms:?}");
    }
}
