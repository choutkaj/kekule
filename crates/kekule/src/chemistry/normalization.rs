use std::fmt;

use crate::core::{AtomId, BondId, BondOrder, Molecule};

/// Failure to publish Kekule's canonical represented chemistry.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NormalizationError {
    /// A meaning-preserving representation rewrite requires an unsupported charge.
    FormalChargeOutOfRange { atom: AtomId, charge: usize },
}

impl fmt::Display for NormalizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FormalChargeOutOfRange { atom, charge } => write!(
                f,
                "normalizing atom {atom} requires formal charge +{charge}, which is outside the supported range"
            ),
        }
    }
}

impl std::error::Error for NormalizationError {}

/// Deterministically normalize represented chemistry into Kekule's canonical form.
///
/// Normalization is transactional and idempotent. A successful call clears the
/// complete installed perception state, including when the primary
/// representation was already normalized.
pub fn normalize_molecule(molecule: &mut Molecule) -> Result<(), NormalizationError> {
    let mut staged = molecule.clone();
    normalize_hypervalent_oxo_halides(&mut staged)?;
    staged.invalidate_topology();
    *molecule = staged;
    Ok(())
}

fn normalize_hypervalent_oxo_halides(molecule: &mut Molecule) -> Result<(), NormalizationError> {
    let halogens = molecule
        .atoms()
        .filter_map(|(atom_id, atom)| {
            (atom.formal_charge == 0
                && matches!(atom.element.symbol(), "Cl" | "Br" | "I")
                && has_terminal_single_bond_oxygen_neighbor(molecule, atom_id))
            .then_some(atom_id)
        })
        .collect::<Vec<_>>();

    for atom_id in halogens {
        let oxo_bonds = oxo_bonds_to_neutral_oxygen(molecule, atom_id);
        if oxo_bonds.is_empty() {
            continue;
        }
        let charge = oxo_bonds.len();
        let formal_charge =
            i8::try_from(charge).map_err(|_| NormalizationError::FormalChargeOutOfRange {
                atom: atom_id,
                charge,
            })?;

        if let Some(atom) = molecule.atoms[atom_id.index()].as_mut() {
            atom.formal_charge = formal_charge;
        }
        for (oxygen_id, bond_id) in oxo_bonds {
            if let Some(atom) = molecule.atoms[oxygen_id.index()].as_mut() {
                atom.formal_charge = -1;
            }
            if let Some(bond) = molecule.bonds[bond_id.index()].as_mut() {
                bond.order = BondOrder::Single;
            }
        }
    }
    Ok(())
}

fn has_terminal_single_bond_oxygen_neighbor(molecule: &Molecule, atom_id: AtomId) -> bool {
    molecule
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(_, bond)| {
            let oxygen_id = bond.other_atom(atom_id);
            bond.order == BondOrder::Single
                && molecule
                    .atom(oxygen_id)
                    .is_ok_and(|neighbor| neighbor.element.symbol() == "O")
                && molecule.incident_bonds(oxygen_id).is_ok_and(|mut bonds| {
                    bonds.all(|(_, oxygen_bond)| {
                        let neighbor_id = oxygen_bond.other_atom(oxygen_id);
                        neighbor_id == atom_id
                            || molecule
                                .atom(neighbor_id)
                                .is_ok_and(|neighbor| neighbor.element.symbol() == "H")
                    })
                })
        })
}

fn oxo_bonds_to_neutral_oxygen(molecule: &Molecule, atom_id: AtomId) -> Vec<(AtomId, BondId)> {
    molecule
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter_map(|(bond_id, bond)| {
            if bond.order != BondOrder::Double {
                return None;
            }
            let oxygen_id = bond.other_atom(atom_id);
            let oxygen = molecule.atom(oxygen_id).ok()?;
            (oxygen.element.symbol() == "O" && oxygen.formal_charge == 0)
                .then_some((oxygen_id, bond_id))
        })
        .collect()
}
