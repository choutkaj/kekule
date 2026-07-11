use crate::core::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValenceModel {
    RdkitLike,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ValenceOptions {
    pub strict: bool,
}

impl Default for ValenceOptions {
    fn default() -> Self {
        Self { strict: true }
    }
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
    perceive_valence_with_options(mol, model, ValenceOptions::default())
}

pub fn perceive_valence_with_options(
    mol: &mut Molecule,
    model: ValenceModel,
    options: ValenceOptions,
) -> ValenceReport {
    match model {
        ValenceModel::RdkitLike => perceive_rdkit_like_valence(mol, options),
    }
}

fn perceive_rdkit_like_valence(mol: &mut Molecule, options: ValenceOptions) -> ValenceReport {
    let mut assignments = Vec::<(AtomId, u8)>::new();
    let mut issues = Vec::new();
    for (atom_id, atom) in mol.atoms() {
        let explicit = explicit_valence(mol, atom_id).saturating_add(atom.explicit_hydrogens);
        if has_rdkit_unrestricted_valence(atom) {
            assignments.push((atom_id, 0));
            continue;
        }
        match allowed_valences(atom) {
            Some(allowed) => {
                if let Some(max_allowed) = allowed.iter().copied().max() {
                    if explicit > max_allowed {
                        if options.strict {
                            issues.push(ValenceIssue::ValenceExceeded {
                                atom: atom_id,
                                explicit_valence: explicit,
                                max_allowed,
                            });
                        }
                        assignments.push((atom_id, 0));
                    } else {
                        let target = if atom.no_implicit_hydrogens
                            || rdkit_suppresses_implicit_hydrogens(atom)
                        {
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
            None if explicit == 0 => assignments.push((atom_id, 0)),
            None => {
                if options.strict {
                    issues.push(ValenceIssue::UnsupportedElement(atom_id));
                }
                assignments.push((atom_id, 0));
            }
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

fn rdkit_suppresses_implicit_hydrogens(atom: &Atom) -> bool {
    matches!(
        atom.element.atomic_number(),
        3 | 4 | 11 | 12 | 19 | 20 | 37 | 38 | 55 | 56 | 87 | 88
    )
}

fn has_rdkit_unrestricted_valence(atom: &Atom) -> bool {
    rdkit_atomic_number_has_unrestricted_valence(atom.element.atomic_number())
        || (i16::from(atom.element.atomic_number()) - i16::from(atom.formal_charge))
            .try_into()
            .ok()
            .is_some_and(rdkit_atomic_number_has_unrestricted_valence)
}

fn rdkit_atomic_number_has_unrestricted_valence(atomic_number: u8) -> bool {
    matches!(atomic_number, 21..=30 | 39..=48 | 57..=81 | 89..=118)
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
        "P" => explicit,
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

pub(crate) fn allowed_valences(atom: &Atom) -> Option<&'static [u8]> {
    if atom.formal_charge != 0 {
        return rdkit_charge_adjusted_allowed_valences(atom);
    }
    match (atom.element.symbol(), atom.formal_charge) {
        ("H", 0) => Some(&[1]),
        ("H", -1 | 1) => Some(&[0]),
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
        ("Li" | "Na" | "K" | "Rb" | "Cs", 0) => Some(&[1]),
        ("Li" | "Na" | "K" | "Rb" | "Cs", 1) => Some(&[0]),
        ("Be" | "Mg" | "Ca" | "Sr" | "Ba", 1..=3) => Some(&[0]),
        ("Be" | "Mg" | "Ca" | "Sr" | "Ba", 0) => Some(&[2]),
        ("F" | "Cl" | "Br" | "I", -1) => Some(&[0]),
        ("F", 0) => Some(&[1]),
        ("Cl" | "Br" | "I", 1) => Some(&[2, 4, 6]),
        ("Cl" | "Br" | "I", 2) => Some(&[3, 5]),
        ("Cl" | "Br" | "I", 3) => Some(&[4]),
        ("Cl" | "Br", 0) => Some(&[1]),
        ("I", 0) => Some(&[1, 3, 5]),
        ("Xe", 0) => Some(&[0, 2, 4, 6]),
        ("Po", 0) => Some(&[2, 4, 6]),
        ("At", 0) => Some(&[1, 3, 5]),
        ("P" | "As" | "Sb" | "Bi", 0) => Some(&[3, 5]),
        ("P" | "As", -1) => Some(&[2, 4, 6]),
        ("P" | "As", -3) => Some(&[0]),
        ("P" | "As", 1) => Some(&[4]),
        ("As" | "Sb" | "Bi", 3 | 5) => Some(&[0]),
        ("Al" | "Ga" | "In" | "Tl", 3) => Some(&[0]),
        ("Al" | "Ga" | "In", 0) => Some(&[3]),
        ("Si", 0) => Some(&[4]),
        ("Si", -1) => Some(&[3]),
        ("Si", -2 | 4) => Some(&[0]),
        ("Ge", 0) => Some(&[0, 4]),
        ("Ge", 4) => Some(&[0]),
        ("Sn", 0) => Some(&[0, 2, 4]),
        ("Sn", -1) => Some(&[3]),
        ("Sn", 2..=4) => Some(&[0]),
        ("Pb", 0) => Some(&[0, 2, 4]),
        ("Pb", 2 | 4) => Some(&[0]),
        ("S" | "Se" | "Te", 0) => Some(&[2, 4, 6]),
        ("S" | "Se" | "Te", -2) => Some(&[0]),
        ("S" | "Se" | "Te", -1 | 1) => Some(&[1, 3, 5]),
        _ => None,
    }
}

pub(crate) fn rdkit_default_valence(atom: &Atom) -> Option<u8> {
    rdkit_default_valence_for_atomic_number(atom.element.atomic_number())
}

pub(crate) fn rdkit_charge_adjusted_default_valence(atom: &Atom) -> Option<u8> {
    let adjusted = i16::from(atom.element.atomic_number()) - i16::from(atom.formal_charge);
    rdkit_default_valence_for_atomic_number(u8::try_from(adjusted).ok()?)
}

fn rdkit_default_valence_for_atomic_number(atomic_number: u8) -> Option<u8> {
    match atomic_number {
        1 | 3 | 9 | 11 | 17 | 19 | 35 | 37 | 53 | 55 | 85 | 87 => Some(1),
        2 | 10 | 18 | 36 | 54 | 86 => Some(0),
        4 | 8 | 12 | 16 | 20 | 34 | 38 | 50 | 52 | 56 | 82 | 84 | 88 => Some(2),
        5 | 7 | 13 | 15 | 31 | 33 | 49 | 51 | 83 => Some(3),
        6 | 14 | 32 => Some(4),
        _ => None,
    }
}

fn rdkit_charge_adjusted_allowed_valences(atom: &Atom) -> Option<&'static [u8]> {
    let adjusted = i16::from(atom.element.atomic_number()) - i16::from(atom.formal_charge);
    if adjusted == 0 {
        return Some(&[0]);
    }
    let adjusted = u8::try_from(adjusted).ok()?;
    rdkit_neutral_allowed_valences(adjusted)
}

fn rdkit_neutral_allowed_valences(atomic_number: u8) -> Option<&'static [u8]> {
    match atomic_number {
        1 => Some(&[1]),
        2 | 10 | 18 | 36 | 86 => Some(&[0]),
        3 | 11 | 19 | 37 | 55 | 87 => Some(&[1]),
        4 | 12 | 20 | 38 | 56 | 88 => Some(&[2]),
        5 => Some(&[3]),
        6 => Some(&[4]),
        7 => Some(&[3]),
        8 => Some(&[2]),
        9 => Some(&[1]),
        13 | 31 | 49 => Some(&[3]),
        14 | 32 => Some(&[4]),
        15 | 33 | 51 | 83 => Some(&[3, 5]),
        16 | 34 | 52 | 84 => Some(&[2, 4, 6]),
        17 => Some(&[1]),
        35 => Some(&[1]),
        50 | 82 => Some(&[2, 4]),
        53 | 85 => Some(&[1, 3, 5]),
        54 => Some(&[0, 2, 4, 6]),
        _ => None,
    }
}
