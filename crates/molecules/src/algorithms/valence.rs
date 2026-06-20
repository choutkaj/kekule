use crate::core::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValenceModel {
    RdkitLike,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValenceIssue {
    UnsupportedElement(AtomId),
    ValenceExceeded {
        atom: AtomId,
        explicit_valence: u8,
        max_allowed: u8,
    },
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ValenceReport {
    pub issues: Vec<ValenceIssue>,
}

impl ValenceReport {
    pub fn is_ok(&self) -> bool {
        self.issues.is_empty()
    }
}

pub fn perceive_valence(mol: &mut Molecule, model: ValenceModel) -> ValenceReport {
    match model {
        ValenceModel::RdkitLike => perceive_rdkit_like_valence(mol),
    }
}

fn perceive_rdkit_like_valence(mol: &mut Molecule) -> ValenceReport {
    let mut assignments = Vec::<(AtomId, u8)>::new();
    let mut issues = Vec::new();
    for (atom_id, atom) in mol.atoms() {
        let explicit = explicit_valence(mol, atom_id).saturating_add(atom.explicit_hydrogens);
        match allowed_valences(atom) {
            Some(allowed) => {
                if let Some(max_allowed) = allowed.iter().copied().max() {
                    if explicit > max_allowed {
                        issues.push(ValenceIssue::ValenceExceeded {
                            atom: atom_id,
                            explicit_valence: explicit,
                            max_allowed,
                        });
                        assignments.push((atom_id, 0));
                    } else {
                        let target = if atom.no_implicit_hydrogens {
                            explicit
                        } else if let Some(target) = aromatic_valence_target(atom, explicit) {
                            target
                        } else if atom.radical.is_some() {
                            explicit
                        } else {
                            allowed
                                .iter()
                                .copied()
                                .find(|allowed| *allowed >= explicit)
                                .unwrap_or(explicit)
                        };
                        assignments.push((atom_id, target - explicit));
                    }
                }
            }
            None => issues.push(ValenceIssue::UnsupportedElement(atom_id)),
        }
    }
    let mut changed = false;
    for (atom_id, hydrogens) in assignments {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            changed |= atom.implicit_hydrogens != Some(hydrogens);
            atom.implicit_hydrogens = Some(hydrogens);
        }
    }
    if changed {
        mol.perception.aromaticity = invalidate(mol.perception.aromaticity);
        mol.perception.stereo = invalidate(mol.perception.stereo);
    }
    mol.perception.valence = ComputedState::Fresh;
    ValenceReport { issues }
}

fn aromatic_valence_target(atom: &Atom, explicit: u8) -> Option<u8> {
    if !atom.aromatic {
        return None;
    }
    let target = match atom.element.symbol() {
        "B" | "C" => 3,
        "N" | "O" | "S" | "Se" | "Te" => {
            if atom.explicit_hydrogens > 0 || atom.formal_charge > 0 {
                3
            } else {
                2
            }
        }
        "P" => 3,
        _ => return None,
    };
    Some(target.max(explicit))
}

pub(crate) fn explicit_valence(mol: &Molecule, atom: AtomId) -> u8 {
    mol.incident_bonds(atom)
        .ok()
        .into_iter()
        .flatten()
        .map(|(_, bond)| bond_order_valence(bond.order))
        .sum()
}

fn bond_order_valence(order: BondOrder) -> u8 {
    match order {
        BondOrder::Zero | BondOrder::Dative => 0,
        BondOrder::Single | BondOrder::Aromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
    }
}

fn allowed_valences(atom: &Atom) -> Option<&'static [u8]> {
    match (atom.element.symbol(), atom.formal_charge) {
        ("H", 0) => Some(&[1]),
        ("B", -1) => Some(&[4]),
        ("B", _) => Some(&[3]),
        ("C", 0) => Some(&[4]),
        ("C", 1 | -1) => Some(&[3]),
        ("N", 1) => Some(&[4]),
        ("N", -1) => Some(&[2]),
        ("N", 0) => Some(&[3, 5]),
        ("O", 0) => Some(&[2]),
        ("O", -1) => Some(&[1]),
        ("O", -2) => Some(&[0]),
        ("O", 1) => Some(&[1, 3]),
        ("Li" | "Na" | "K", 1) => Some(&[0]),
        ("Mg" | "Ca" | "Sr" | "Ba", 1..=3) => Some(&[0]),
        ("Nb", 5) => Some(&[0]),
        (
            "Sc" | "Ti" | "V" | "Cr" | "Mn" | "Fe" | "Co" | "Ni" | "Cu" | "Zn" | "Y" | "Zr" | "Nb"
            | "Mo" | "Tc" | "Ru" | "Rh" | "Pd" | "Ag" | "Cd" | "La" | "Ce" | "Pr" | "Nd" | "Sm"
            | "Eu" | "Gd" | "Tb" | "Dy" | "Ho" | "Er" | "Tm" | "Yb" | "Lu" | "Hf" | "Ta" | "W"
            | "Re" | "Os" | "Ir" | "Pt" | "Au" | "Hg",
            -1..=4,
        ) => Some(&[0, 1, 2, 3, 4, 5, 6]),
        ("F" | "Cl" | "Br" | "I", -1) => Some(&[0]),
        ("F", 0) => Some(&[1]),
        ("Cl" | "Br" | "I", 1) => Some(&[2, 4]),
        ("Cl" | "Br" | "I", 3) => Some(&[4]),
        ("Cl" | "Br" | "I", 0) => Some(&[1, 3, 5, 7]),
        ("P" | "As" | "Sb" | "Bi", 0) => Some(&[3, 5]),
        ("Bi", 3) => Some(&[0]),
        ("In", 3) => Some(&[0]),
        ("P", 1) => Some(&[4]),
        ("Si", 0) => Some(&[4]),
        ("Ge", 0) => Some(&[0, 4]),
        ("Sn", 0) => Some(&[0, 2, 4]),
        ("Sn", 4) => Some(&[0]),
        ("U", 2) => Some(&[4]),
        ("S" | "Se" | "Te", 0) => Some(&[2, 4, 6]),
        ("S" | "Se" | "Te", -1 | 1) => Some(&[1, 3, 5]),
        ("Tl", 0) => Some(&[3]),
        _ => None,
    }
}
