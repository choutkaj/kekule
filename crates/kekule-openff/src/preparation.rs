//! Compatibility with RDKit's nitrogen/phosphorus cleanup on a private copy.
//! See THIRD_PARTY.md for the reference revision and BSD attribution.
use crate::{Error, ErrorKind, Result};
use kekule::core::{BondOrder, Molecule};

/// A perceived private copy of an explicit-hydrogen molecule, with
/// RDKit-compatible valence normalization, ready for assignment.
///
/// Atom IDs are preserved, so per-atom results computed on the copy apply
/// to the input atom for atom.
pub(crate) fn explicit(molecule: &Molecule) -> Result<Molecule> {
    let mut copy = normalize_valence(molecule)?;
    copy.perceive().map_err(Error::chemistry)?;
    if !copy.atom_ids().eq(molecule.atom_ids()) {
        return Err(Error::new(
            ErrorKind::Chemistry,
            "molecule preparation changed atom identities",
        ));
    }
    for (id, atom) in copy.atoms() {
        if copy.implicit_hydrogens(id).map_err(Error::chemistry)? != Some(0) {
            return Err(Error::new(
                ErrorKind::UnsupportedMolecule,
                format!(
                    "atom {id} has implicit hydrogens; expand hydrogens before parameterization"
                ),
            ));
        }
        if atom.radical.is_some() {
            return Err(Error::new(
                ErrorKind::UnsupportedMolecule,
                "radicals are outside the supported parameterization domain",
            ));
        }
    }
    Ok(copy)
}

pub(crate) fn normalize_valence(input: &Molecule) -> Result<Molecule> {
    let mut copy = input.clone();
    let mut nitrogen = Vec::new();
    for (id, atom) in input.atoms() {
        if atom.element.atomic_number() != 7 || atom.formal_charge != 0 {
            continue;
        }
        let valence: u32 = input
            .neighbors(id)
            .map_err(Error::chemistry)?
            .map(|next| {
                match input
                    .bond(input.bond_between(id, next).unwrap().unwrap())
                    .unwrap()
                    .order
                {
                    BondOrder::Single => 1,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    _ => 100,
                }
            })
            .sum();
        if valence == 5 {
            nitrogen.push(id);
        }
    }
    // Same two passes as RDKit: neutral N=O first, then neutral N#N.
    for (element, old, new) in [
        (8, BondOrder::Double, BondOrder::Single),
        (7, BondOrder::Triple, BondOrder::Double),
    ] {
        for &id in &nitrogen {
            let found = copy
                .neighbors(id)
                .map_err(Error::chemistry)?
                .find_map(|next| {
                    let atom = copy.atom(next).unwrap();
                    let bond = copy.bond_between(id, next).unwrap().unwrap();
                    (atom.element.atomic_number() == element
                        && atom.formal_charge == 0
                        && copy.bond(bond).unwrap().order == old)
                        .then_some((next, bond))
                });
            if let Some((next, bond)) = found {
                let mut editor = copy.edit();
                editor.atom_mut(id).map_err(Error::chemistry)?.formal_charge = 1;
                editor
                    .atom_mut(next)
                    .map_err(Error::chemistry)?
                    .formal_charge = -1;
                editor
                    .bond_mut(bond)
                    .map_err(Error::chemistry)?
                    .set_order(new);
                copy = editor.finish().map_err(Error::chemistry)?;
            }
        }
    }
    let copy = normalize_phosphorus(&copy)?;
    if copy.formal_charge() != input.formal_charge() {
        return Err(Error::new(
            ErrorKind::Chemistry,
            "input valence normalization changed total charge",
        ));
    }
    Ok(copy)
}

// A neutral, degree-three, valence-five P with both P=O and P=C/N becomes
// [P+](-[O-])=C/N. This preserves atom identities, connectivity and total charge.
// It is necessary for the fixed-H lookup key of Ash entry 5868.
pub(crate) fn normalize_phosphorus(input: &Molecule) -> Result<Molecule> {
    let mut edits = Vec::new();
    for (id, atom) in input.atoms() {
        if atom.element.atomic_number() != 15 || atom.formal_charge != 0 {
            continue;
        }
        let neighbors = input
            .neighbors(id)
            .map_err(Error::chemistry)?
            .collect::<Vec<_>>();
        if neighbors.len() != 3 {
            continue;
        }
        let mut valence = 0;
        let mut oxygen = None;
        let mut carbon_or_nitrogen = false;
        for next in neighbors {
            let bond_id = input
                .bond_between(id, next)
                .map_err(Error::chemistry)?
                .unwrap();
            let bond = input.bond(bond_id).map_err(Error::chemistry)?;
            valence += match bond.order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                _ => 100,
            };
            let neighbor = input.atom(next).map_err(Error::chemistry)?;
            if bond.order == BondOrder::Double {
                if neighbor.element.atomic_number() == 8 && neighbor.formal_charge == 0 {
                    oxygen = Some((next, bond_id));
                } else if matches!(neighbor.element.atomic_number(), 6 | 7)
                    && input.neighbors(next).map_err(Error::chemistry)?.count() >= 2
                {
                    carbon_or_nitrogen = true;
                }
            }
        }
        if valence == 5 && carbon_or_nitrogen {
            if let Some((oxygen, bond)) = oxygen {
                edits.push((id, oxygen, bond));
            }
        }
    }
    if edits.is_empty() {
        return Ok(input.clone());
    }
    let mut editor = input.edit();
    for (p, o, bond) in edits {
        editor.atom_mut(p).map_err(Error::chemistry)?.formal_charge = 1;
        editor.atom_mut(o).map_err(Error::chemistry)?.formal_charge = -1;
        editor
            .bond_mut(bond)
            .map_err(Error::chemistry)?
            .set_order(BondOrder::Single);
    }
    editor.finish().map_err(Error::chemistry)
}

#[cfg(test)]
mod tests {
    // Expected graphs are RDKit's sanitized forms, which upstream identifies.
    fn same_graph(prepared: &kekule::core::Molecule, expected: &str) -> bool {
        let mut expected = kekule::smiles::to_molecules(expected).unwrap().remove(0);
        expected.perceive().unwrap();
        crate::identity::same_graph(prepared, &expected).unwrap()
    }

    #[test]
    fn neutral_valence_five_nitrogen_matches_charge_separated_input() {
        for (s, expected) in [
            ("[H:4][N:2](=[O:1])=[O:3]", "[H][N+](=O)[O-]"),
            ("[H:4][N:2](#[N:3])[F:1]", "[H][N+](=[N-])F"),
        ] {
            let input = kekule::smiles::to_molecules(s).unwrap().remove(0);
            let prepared = crate::explicit(&input).unwrap();
            assert!(same_graph(&prepared, expected), "{s}");
            assert!(input.atoms().all(|(_, a)| a.formal_charge == 0));
        }
    }
    #[test]
    fn ash_phosphorus_entry_matches_rdkit_without_mutating_input() {
        let input = kekule::smiles::to_molecules("[H:5][N:3]=[P:2](=[O:1])[H:4]")
            .unwrap()
            .remove(0);
        let normalized = crate::explicit(&input).unwrap();
        assert!(same_graph(&normalized, "[H]N=[P+]([H])[O-]"));
        assert!(input.atoms().all(|(_, a)| a.formal_charge == 0));
        assert_eq!(
            input.atom_ids().collect::<Vec<_>>(),
            normalized.atom_ids().collect::<Vec<_>>()
        );
        assert_eq!(input.formal_charge(), normalized.formal_charge());
    }
}
