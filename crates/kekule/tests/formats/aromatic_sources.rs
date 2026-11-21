//! Localization of aromatic bonds that SMILES and Molfile sources state.

use kekule::molecule::{BondOrder, Molecule};
use kekule::{molfile, smiles};

fn canonical(mut molecule: Molecule) -> String {
    molecule.perceive().unwrap();
    smiles::write(&molecule, smiles::SmilesWriteOptions::canonical()).unwrap()
}

fn smiles_molecule(source: &str) -> Molecule {
    smiles::to_molecules(source)
        .unwrap_or_else(|error| panic!("{source}: {error}"))
        .pop()
        .unwrap()
}

/// A V2000 ring of `atoms` (symbol and charge) on a regular polygon, bond
/// `i` joining atoms `i` and `i + 1` with the V2000 bond type `bonds[i]`.
fn v2000_ring(atoms: &[(&str, i8)], bonds: &[u8]) -> String {
    let count = atoms.len();
    let mut text =
        format!("ring\nkekule\n\n{count:>3}{count:>3}  0  0  0  0            999 V2000\n");
    for (index, (symbol, _)) in atoms.iter().enumerate() {
        let angle = std::f64::consts::TAU * index as f64 / count as f64;
        text.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {symbol:<3} 0  0  0  0  0  0\n",
            1.4 * angle.cos(),
            1.4 * angle.sin(),
            0.0
        ));
    }
    for (index, order) in bonds.iter().enumerate() {
        text.push_str(&format!(
            "{:>3}{:>3}{order:>3}  0  0  0  0\n",
            index + 1,
            (index + 1) % count + 1
        ));
    }
    for (index, (_, charge)) in atoms.iter().enumerate() {
        if *charge != 0 {
            text.push_str(&format!("M  CHG  1 {:>3} {charge:>3}\n", index + 1));
        }
    }
    text.push_str("M  END\n");
    text
}

fn molfile_molecule(text: &str) -> Molecule {
    molfile::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap_or_else(|error| panic!("{error}\n{text}"))
        .into_molecules()
        .pop()
        .unwrap()
}

fn double_bonds(molecule: &Molecule) -> usize {
    molecule
        .bonds()
        .filter(|(_, bond)| bond.order() == BondOrder::Double)
        .count()
}

#[test]
fn stated_double_bonds_between_aromatic_atoms_meet_their_pi_bond() {
    let benzene = canonical(smiles_molecule("c1ccccc1"));
    for source in ["c1=cc=cc=c1", "c1=ccccc1", "c1cc=ccc1", "c1:c=c:c=c:c=1"] {
        let molecule = smiles_molecule(source);
        assert_eq!(double_bonds(&molecule), 3, "{source}");
        assert_eq!(canonical(molecule), benzene, "{source}");
    }
    let pyridine = canonical(smiles_molecule("n1ccccc1"));
    for source in ["n1=cc=cc=c1", "n1=ccccc1", "n1c=cccc1"] {
        assert_eq!(canonical(smiles_molecule(source)), pyridine, "{source}");
    }
    // V2000 bond type 2 is double and type 4 aromatic.
    let carbons = [("C", 0); 6];
    for bonds in [[2, 4, 2, 4, 2, 4], [2, 4, 4, 4, 4, 4], [4, 4, 2, 4, 4, 4]] {
        let molecule = molfile_molecule(&v2000_ring(&carbons, &bonds));
        assert_eq!(double_bonds(&molecule), 3, "{bonds:?}");
        assert_eq!(canonical(molecule), benzene, "{bonds:?}");
    }
}

#[test]
fn aromatic_positive_chalcogens_take_no_hydrogen() {
    for (heteroatom, smiles_form) in [
        ("O", "[o+]1ccccc1"),
        ("S", "[s+]1ccccc1"),
        ("Se", "[se+]1ccccc1"),
    ] {
        let mut atoms = vec![(heteroatom, 1)];
        atoms.extend([("C", 0); 5]);
        let molecule = molfile_molecule(&v2000_ring(&atoms, &[4; 6]));
        assert_eq!(
            canonical(molecule),
            canonical(smiles_molecule(smiles_form)),
            "{heteroatom}"
        );
    }
    // 1,3-Dithiolium: the positive sulfur gives the ring its even demand.
    let atoms = [("S", 1), ("C", 0), ("S", 0), ("C", 0), ("C", 0)];
    let molecule = molfile_molecule(&v2000_ring(&atoms, &[4; 5]));
    assert_eq!(
        canonical(molecule),
        canonical(smiles_molecule("[s+]1cscc1"))
    );
    // Positive nitrogen with two ring bonds keeps its hydrogen.
    let mut atoms = vec![("N", 1)];
    atoms.extend([("C", 0); 5]);
    let molecule = molfile_molecule(&v2000_ring(&atoms, &[4; 6]));
    assert_eq!(
        canonical(molecule),
        canonical(smiles_molecule("[nH+]1ccccc1"))
    );
}
