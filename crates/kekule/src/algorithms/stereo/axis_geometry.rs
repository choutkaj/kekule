use crate::core::*;

/// Conservative local trigonal geometry for the three explicit ligands needed
/// by the supported atropisomer model. This reads localized bonds, charge and
/// the caller's hydrogen count, never installed aromaticity or ring membership.
/// Source interpretation supplies declared H; coordinate inference supplies
/// total implicit H. A multiple bond is one ligand direction, while a localized
/// lone pair adds a domain (notably in sulfoxides and phosphines).
pub(crate) fn atom_is_atropisomeric_sp2_endpoint(
    mol: &Molecule,
    center: AtomId,
    hydrogens: u8,
) -> bool {
    let Ok(atom) = mol.atom(center) else {
        return false;
    };
    if hydrogens != 0 || atom.radical.is_some() {
        return false;
    }
    let Some(outer) = main_group_outer_electrons(atom) else {
        return false;
    };
    let Ok(incident) = mol.incident_bonds(center) else {
        return false;
    };
    let bonds = incident.map(|(_, bond)| bond).collect::<Vec<_>>();
    if bonds.len() != 3 {
        return false;
    }
    let Some(valence) = bonds.iter().try_fold(0_i16, |sum, bond| {
        Some(
            sum + match bond.order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                _ => return None,
            },
        )
    }) else {
        return false;
    };
    let nonbonding = outer - i16::from(atom.formal_charge) - valence;
    if nonbonding == 0 && valence <= 4 {
        return true;
    }
    // A second-row lone pair can occupy a conjugated p orbital (e.g. pyrrole
    // or amide N). Heavier pyramidal donors retain their localized lone pair.
    // Require an actual adjacent pi system; ring membership alone says nothing
    // about geometry and would accept saturated rings.
    nonbonding == 2
        && valence == 3
        && atom.element.atomic_number() <= 10
        && bonds.iter().any(|bond| {
            let neighbor = bond.other_atom(center);
            let Ok(other) = mol.atom(neighbor) else {
                return false;
            };
            other.radical.is_none()
                && main_group_outer_electrons(other)
                    .is_some_and(|outer| other.element.atomic_number() <= 10 || outer <= 4)
                && mol
                    .incident_bonds(neighbor)
                    .ok()
                    .is_some_and(|mut incident| {
                        incident.any(|(_, next)| {
                            next.other_atom(neighbor) != center
                                && matches!(next.order, BondOrder::Double | BondOrder::Triple)
                        })
                    })
        })
}

fn main_group_outer_electrons(atom: &Atom) -> Option<i16> {
    // The supported trigonal subset is p-block groups 13 through 16. Transition
    // metal coordination and expanded-valence geometry need a different model.
    Some(match atom.element.atomic_number() {
        5..=8 => i16::from(atom.element.atomic_number()) - 2,
        13..=16 => i16::from(atom.element.atomic_number()) - 10,
        31..=34 => i16::from(atom.element.atomic_number()) - 28,
        49..=52 => i16::from(atom.element.atomic_number()) - 46,
        81..=84 => i16::from(atom.element.atomic_number()) - 78,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::read_smiles;

    #[test]
    fn trigonal_axis_endpoints_use_ligand_and_lone_pair_geometry() {
        for (source, expected) in [
            ("C(=C)(F)Cl", true),
            ("C1(Cl)=CC=CC=C1", true),
            ("N1(C)C=CC=C1", true),
            ("N(C)(C)C(=O)C", true),
            ("[C+](C)(C)C", true),
            ("B(C)(C)C", true),
            ("C1(C)CCCCC1", false),
            ("N1(C)CCCCC1", false),
            ("S(=O)(C)C1=CC=CC=C1", false),
            ("[Se](=O)(C)C1=CC=CC=C1", false),
            ("[S+](C)(C)C1=CC=CC=C1", false),
            ("P(C)(C)C1=CC=CC=C1", false),
            ("N(C)(C)S(=O)C", false),
        ] {
            let mol = read_smiles(source).unwrap();
            assert_eq!(
                atom_is_atropisomeric_sp2_endpoint(&mol, AtomId::new(0), 0),
                expected,
                "{source}"
            );
            assert!(!atom_is_atropisomeric_sp2_endpoint(&mol, AtomId::new(0), 1));
        }
    }
}
