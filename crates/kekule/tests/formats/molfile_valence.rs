//! Molfile atoms whose drawing leaves their charge or hydrogens implicit:
//! onium centres, their donor bonds, and metals.

use kekule::molecule::{AtomId, BondOrder, Element, Molecule, MoleculeEditor};
use kekule::{molfile, smiles};

use molfile::MolfileInterpretIssue as Issue;

/// A V2000 record of `atoms` (symbol and charge) at spaced positions, bonded
/// by single `bonds` between one-based atom numbers.
fn v2000(atoms: &[(&str, i8)], bonds: &[(usize, usize)]) -> String {
    let bonds = bonds.iter().map(|&(a, b)| (a, b, 1)).collect::<Vec<_>>();
    v2000_orders(atoms, &bonds)
}

/// [`v2000`] with a Molfile bond type for each bond.
fn v2000_orders(atoms: &[(&str, i8)], bonds: &[(usize, usize, u8)]) -> String {
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
    for (a, b, order) in bonds {
        text.push_str(&format!("{a:>3}{b:>3}{order:>3}  0  0  0  0\n"));
    }
    text.push_str("M  END\n");
    text
}

/// Atom `center` singly bonded to `carbons` carbons and then to the `extra`
/// atoms.
fn centre(center: &str, carbons: usize, extra: &[&str]) -> String {
    let mut atoms = vec![(center, 0)];
    atoms.extend(std::iter::repeat_n(("C", 0), carbons));
    atoms.extend(extra.iter().map(|symbol| (*symbol, 0)));
    let bonds = (2..=atoms.len()).map(|atom| (1, atom)).collect::<Vec<_>>();
    v2000(&atoms, &bonds)
}

fn interpret(text: &str) -> molfile::MolfileInterpretation {
    molfile::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// The published atom and source line of the `index`th source atom of a
/// one-molecule interpretation.
fn source_atom(
    interpretation: &molfile::MolfileInterpretation,
    index: usize,
) -> (AtomId, kekule::diagnostics::SourceLocation) {
    let mapping = interpretation.reports()[0].atom_mappings()[index];
    (mapping.atom(), mapping.location())
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
        let interpretation = interpret(&centre(center, ligands, &[]));
        let molecule = interpretation.molecules().next().unwrap();
        assert_eq!(
            perceived(molecule)[0],
            (center.to_owned(), 1, 0),
            "{center} with {ligands} carbons"
        );
        let (atom, location) = source_atom(&interpretation, 0);
        assert_eq!(
            interpretation.reports()[0].issues(),
            [Issue::ChargeFromBonds {
                atom,
                location,
                charge: 1,
            }],
            "{center}"
        );
    }
}

#[test]
fn single_bonds_to_terminal_chalcogens_separate_charges() {
    // DMSO, trimethylphosphine oxide, its sulfur analogue, and dimethyl
    // sulfone drawn with one single S-O bond: the chalcogen takes -1, not a
    // hydrogen.
    let sulfone = v2000_orders(
        &[("S", 0), ("C", 0), ("C", 0), ("O", 0), ("O", 0)],
        &[(1, 2, 1), (1, 3, 1), (1, 4, 2), (1, 5, 1)],
    );
    for (text, chalcogen) in [
        (centre("S", 2, &["O"]), 3),
        (centre("P", 3, &["O"]), 4),
        (centre("P", 3, &["S"]), 4),
        (sulfone, 4),
    ] {
        let interpretation = interpret(&text);
        let molecule = interpretation.molecules().next().unwrap();
        let atoms = perceived(molecule);
        assert_eq!((atoms[0].1, atoms[0].2), (1, 0), "{text}");
        assert_eq!((atoms[chalcogen].1, atoms[chalcogen].2), (-1, 0), "{text}");
        let total_charge = atoms.iter().map(|(_, charge, _)| *charge).sum::<i8>();
        assert_eq!(total_charge, 0, "{text}");
        let (atom, location) = source_atom(&interpretation, 0);
        assert_eq!(
            interpretation.reports()[0].issues(),
            [Issue::ChargesSeparated {
                atom,
                acceptor: source_atom(&interpretation, chalcogen).0,
                location,
            }],
            "{text}"
        );
    }
}

#[test]
fn donor_bonds_to_boron_separate_charges() {
    // Trimethylphosphine-borane and dimethyl sulfide-borane: P+ and S+ with
    // BH3-, not a cation with BH2.
    for (center, carbons) in [("P", 3), ("S", 2)] {
        let interpretation = interpret(&centre(center, carbons, &["B"]));
        let molecule = interpretation.molecules().next().unwrap();
        let atoms = perceived(molecule);
        assert_eq!(atoms[0], (center.to_owned(), 1, 0));
        assert_eq!(atoms[carbons + 1], ("B".to_owned(), -1, 3));
        let (atom, location) = source_atom(&interpretation, 0);
        assert_eq!(
            interpretation.reports()[0].issues(),
            [Issue::ChargesSeparated {
                atom,
                acceptor: source_atom(&interpretation, carbons + 1).0,
                location,
            }]
        );
    }
}

#[test]
fn donor_bonds_to_metals_take_no_charge_or_hydrogen() {
    // Trimethylphosphine on palladium and dimethyl sulfide on platinum.
    for (center, carbons, metal) in [("P", 3, "Pd"), ("S", 2, "Pt")] {
        let interpretation = interpret(&centre(center, carbons, &[metal]));
        let molecule = interpretation.molecules().next().unwrap();
        let atoms = perceived(molecule);
        assert_eq!(atoms[0], (center.to_owned(), 0, 0));
        assert_eq!(atoms[carbons + 1], (metal.to_owned(), 0, 0));
        let (atom, location) = source_atom(&interpretation, 0);
        assert_eq!(
            interpretation.reports()[0].issues(),
            [Issue::MetalDonor { atom, location }]
        );
    }
}

#[test]
fn stated_charges_and_default_valences_are_kept() {
    // Charged sulfonium, dimethyl sulfide, and trimethylphosphine need
    // nothing.
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

fn write(molecule: &Molecule, version: molfile::MolfileWriteVersion) -> String {
    molfile::write(
        molecule,
        molfile::MolfileWriteOptions::default().with_version(version),
    )
    .unwrap()
    .0
}

const VERSIONS: [molfile::MolfileWriteVersion; 2] = [
    molfile::MolfileWriteVersion::V2000,
    molfile::MolfileWriteVersion::V3000,
];

#[test]
fn molfiles_written_from_onium_hydrides_read_back_unchanged() {
    // The hydrogen of CS(C)C is not the Molfile default, so the writer states
    // the valence and the reader keeps the hydride neutral.
    let molecule = smiles::to_molecules("CS(C)C").unwrap().remove(0);
    for version in VERSIONS {
        let interpretation = interpret(&write(&molecule, version));
        assert_eq!(
            perceived(interpretation.molecules().next().unwrap())[1],
            ("S".to_owned(), 0, 1),
            "{version:?}"
        );
        assert!(interpretation.reports()[0].issues().is_empty());
    }
}

/// The stated valence of each atom of a written Molfile: the V2000 `vvv`
/// field or the V3000 `VAL` property, `None` when absent.
fn stated_valences(text: &str) -> Vec<Option<i32>> {
    if text.contains("V3000") {
        text.lines()
            .skip_while(|line| !line.ends_with("BEGIN ATOM"))
            .skip(1)
            .take_while(|line| !line.ends_with("END ATOM"))
            .map(|line| {
                line.split_whitespace()
                    .find_map(|field| field.strip_prefix("VAL="))
                    .map(|value| value.parse().unwrap())
            })
            .collect()
    } else {
        let count = text.lines().nth(3).unwrap()[..3].trim().parse().unwrap();
        text.lines()
            .skip(4)
            .take(count)
            .map(|line| match line[48..51].trim().parse().unwrap() {
                0 => None,
                value => Some(value),
            })
            .collect()
    }
}

#[test]
fn molfiles_state_the_hydrogens_readers_would_infer_otherwise() {
    // RDKit and MDL readers give unbonded sodium a hydride, the tin of
    // trimethylstannane drawn without its hydrogen another, and a hydride
    // joined to iron by a dative bond a hydrogen. This model infers none.
    let sodium = interpret(&v2000(&[("Na", 0)], &[]));
    let stannyl = interpret(&centre("Sn", 3, &[]));
    let mut editor = MoleculeEditor::new();
    let iron = editor
        .add_atom(kekule::molecule::Atom::new(
            Element::from_symbol("Fe").unwrap(),
        ))
        .unwrap();
    let hydride = editor
        .add_atom(kekule::molecule::Atom::new(
            Element::from_symbol("H").unwrap(),
        ))
        .unwrap();
    editor.add_bond(iron, hydride, BondOrder::Dative).unwrap();
    let iron_hydride = editor.finish().unwrap();
    for version in VERSIONS {
        let zero = match version {
            molfile::MolfileWriteVersion::V2000 => 15,
            _ => -1,
        };
        let cases = [
            (sodium.molecules().next().unwrap(), vec![Some(zero)]),
            (
                stannyl.molecules().next().unwrap(),
                vec![Some(3), None, None, None],
            ),
            (&iron_hydride, vec![None, Some(zero)]),
        ];
        for (molecule, expected) in cases {
            let text = write(molecule, version);
            assert_eq!(stated_valences(&text), expected, "{text}");
            let read = interpret(&text);
            assert_eq!(
                perceived(read.molecules().next().unwrap()),
                perceived(molecule),
                "{text}"
            );
            assert!(read.reports()[0].issues().is_empty(), "{text}");
        }
    }
}

#[test]
fn unbonded_metals_gain_no_hydrides_and_are_reported() {
    for symbol in [
        "Li", "Na", "K", "Cs", "Mg", "Ca", "Ba", "Al", "Ga", "Sn", "Pb", "Bi", "Po",
    ] {
        let interpretation = interpret(&v2000(&[(symbol, 0)], &[]));
        let molecule = interpretation.molecules().next().unwrap();
        assert_eq!(perceived(molecule), [(symbol.to_owned(), 0, 0)], "{symbol}");
        let (atom, location) = source_atom(&interpretation, 0);
        assert_eq!(
            interpretation.reports()[0].issues(),
            [Issue::MetalValenceIncomplete {
                atom,
                location,
                valence: 0,
            }],
            "{symbol}"
        );
    }
}

#[test]
fn under_bonded_metals_gain_no_hydrides_and_are_reported() {
    // Tributylstannane and DIBAL-H drawn without their hydrides, as methyl
    // analogues, a Grignard reagent without its halide, and a bismuth atom
    // with two bonds.
    for (center, carbons, valence) in [("Sn", 3, 3), ("Al", 2, 2), ("Mg", 1, 1), ("Bi", 2, 2)] {
        let interpretation = interpret(&centre(center, carbons, &[]));
        let molecule = interpretation.molecules().next().unwrap();
        assert_eq!(perceived(molecule)[0], (center.to_owned(), 0, 0));
        let (atom, location) = source_atom(&interpretation, 0);
        assert_eq!(
            interpretation.reports()[0].issues(),
            [Issue::MetalValenceIncomplete {
                atom,
                location,
                valence,
            }],
            "{center}"
        );
    }
}

#[test]
fn metals_at_allowed_valences_are_not_reported() {
    // A sodium ion, dimethylmagnesium, trimethylaluminium, a stannylene, and
    // tetramethyltin.
    for (center, charge, carbons) in [
        ("Na", 1, 0),
        ("Mg", 0, 2),
        ("Al", 0, 3),
        ("Sn", 0, 2),
        ("Sn", 0, 4),
    ] {
        let mut atoms = vec![(center, charge)];
        atoms.extend(std::iter::repeat_n(("C", 0), carbons));
        let bonds = (2..=atoms.len()).map(|atom| (1, atom)).collect::<Vec<_>>();
        let interpretation = interpret(&v2000(&atoms, &bonds));
        let molecule = interpretation.molecules().next().unwrap();
        assert_eq!(perceived(molecule)[0].2, 0, "{atoms:?}");
        assert!(interpretation.reports()[0].issues().is_empty(), "{atoms:?}");
    }
}
