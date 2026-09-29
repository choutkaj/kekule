//! Bond conjugation, following RDKit 2026.03.3 ConjugHybrid.cpp.
// Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors.
// BSD 3-Clause; see LICENSE-RDKit in the crate root.
use super::aromaticity::{count_rdkit_like_atom_pi_electrons, rdkit_outer_electrons};
use super::valence::{explicit_valence, rdkit_default_valence};
use crate::core::*;
use std::collections::BTreeSet;
use std::fmt;

/// Missing prerequisites for standalone conjugation perception.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ConjugationError {
    MissingAromaticity,
    UnknownHydrogens(AtomId),
}
impl fmt::Display for ConjugationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingAromaticity => f.write_str("conjugation requires installed aromaticity"),
            Self::UnknownHydrogens(atom) => {
                write!(f, "conjugation requires known hydrogen counts at {atom}")
            }
        }
    }
}
impl std::error::Error for ConjugationError {}

/// Installs RDKit-like conjugation without altering represented chemistry or CIP.
///
/// Requires aromaticity and known implicit hydrogen counts. All aromatic bonds
/// are conjugated; additional flags follow RDKit's local electronic rules,
/// including its element restrictions and treatment of cumulenes and radicals.
/// Success replaces conjugation and invalidates prepared resonance groups.
/// Failure leaves all state untouched. No resonance structures are enumerated.
pub fn perceive_conjugation(
    mol: &mut Molecule,
    model: ConjugationModel,
) -> std::result::Result<(), ConjugationError> {
    let aromaticity = mol
        .perception()
        .aromaticity_state()
        .ok_or(ConjugationError::MissingAromaticity)?;
    let mut bonds: BTreeSet<_> = aromaticity.bonds().collect();
    let mut profiles = vec![(false, 0); mol.graph.atom_slot_count()];
    for (id, atom) in mol.atoms() {
        let h = mol
            .implicit_hydrogens(id)
            .expect("live atom")
            .ok_or(ConjugationError::UnknownHydrogens(id))?;
        let degree = mol.incident_bonds(id).expect("live atom").count() + h;
        let total_valence = explicit_valence(mol, id) + h;
        let outer = rdkit_outer_electrons(atom);
        let hypervalent = atom.formal_charge == 0
            && rdkit_default_valence(atom).is_some_and(|v| total_valence > usize::from(v));
        let eligible = !hypervalent
            && (atom.element.atomic_number() <= 10
                || !matches!(outer, 5 | 6)
                || outer == 6 && degree < 2)
            && count_rdkit_like_atom_pi_electrons(mol, id, atom).is_some_and(|n| n > 0);
        profiles[id.index()] = (eligible, degree);
    }
    for (id, _) in mol.atoms() {
        let (candidate, degree) = profiles[id.index()];
        if !candidate || !(2..=3).contains(&degree) {
            continue;
        }
        for (first, bond) in mol.incident_bonds(id).expect("live atom") {
            if !(mol.perception().bond_is_aromatic(first) == Some(true)
                || matches!(
                    bond.order,
                    BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
                ))
                || !profiles[bond.other_atom(id).index()].0
            {
                continue;
            }
            for (second, other) in mol.incident_bonds(id).expect("live atom") {
                let (eligible, degree) = profiles[other.other_atom(id).index()];
                if first != second && eligible && degree <= 3 {
                    bonds.extend([first, second]);
                }
            }
        }
    }
    let atoms = bonds
        .iter()
        .flat_map(|&id| {
            let b = mol.bond(id).expect("live bond");
            [b.a(), b.b()]
        })
        .collect();
    mol.perception.conjugation = Some(ConjugationPerception {
        model,
        atoms,
        bonds,
    });
    mol.perception.resonance = None;
    Ok(())
}
