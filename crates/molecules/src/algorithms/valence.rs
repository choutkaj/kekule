use crate::core::*;
use std::collections::BTreeMap;

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
        explicit_valence: usize,
        max_allowed: usize,
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
        let explicit = explicit_valence(mol, atom_id) + usize::from(atom.explicit_hydrogens);
        let radical_electrons = atom
            .radical
            .map_or(0, |radical| usize::from(radical.unpaired_electron_count()));

        let Some(original_rule) = rdkit_neutral_valence_rule(atom.element.atomic_number()) else {
            record_unsupported(&mut issues, atom_id, options);
            assignments.push((atom_id, 0));
            continue;
        };

        // RDKit leaves atoms whose periodic-table entry is only `-1` on that
        // unrestricted rule. All other charged atoms use the isoelectronic
        // neutral element's valence list.
        let effective_rule = if original_rule.is_only_unrestricted() {
            Some(original_rule)
        } else {
            rdkit_effective_atomic_number(atom).and_then(rdkit_neutral_valence_rule)
        };
        let Some(mut target_rule) = effective_rule else {
            record_unsupported(&mut issues, atom_id, options);
            assignments.push((atom_id, 0));
            continue;
        };

        let effective_atomic_number = rdkit_effective_atomic_number(atom);
        let hypervalent_anion = effective_atomic_number
            .is_some_and(|effective| can_be_rdkit_hypervalent_anion(atom, effective));
        let charge_offset = if hypervalent_anion {
            target_rule = original_rule;
            usize::from(atom.formal_charge.unsigned_abs())
        } else {
            0
        };
        let occupied_for_target = explicit + radical_electrons + charge_offset;

        // Explicit-valence checking in RDKit honors unrestricted sentinels in
        // either the original or effective valence list. Negatively charged
        // P/S/As/Se instead use their original hypervalent limit with the
        // charge offset applied.
        let two_coordinate_hydride = atom.element.atomic_number() == 1 && atom.formal_charge == -1;
        let explicit_limit = if two_coordinate_hydride {
            // Historical RDKit compatibility: two-coordinate hydride is
            // accepted even though it is chemically unusual.
            Some(2)
        } else if hypervalent_anion {
            target_rule
                .max_fixed()
                .map(|maximum| maximum.saturating_sub(charge_offset))
        } else if original_rule.unrestricted_above || target_rule.unrestricted_above {
            None
        } else {
            target_rule.max_fixed()
        };
        let target_limit = if two_coordinate_hydride || target_rule.unrestricted_above {
            None
        } else {
            target_rule
                .max_fixed()
                .map(|maximum| maximum.saturating_sub(radical_electrons + charge_offset))
        };
        let max_allowed = explicit_limit.into_iter().chain(target_limit).min();
        if max_allowed.is_some_and(|maximum| explicit > maximum) {
            if options.strict {
                issues.push(ValenceIssue::ValenceExceeded {
                    atom: atom_id,
                    explicit_valence: explicit,
                    max_allowed: max_allowed.expect("checked as present"),
                });
            }
            assignments.push((atom_id, 0));
            continue;
        }

        let implicit = if atom.no_implicit_hydrogens {
            0
        } else if let Some(target) =
            aromatic_valence_target(atom, mol.atom_is_aromatic(atom_id).ok().flatten(), explicit)
        {
            target.saturating_sub(explicit)
        } else if let Some(target) = target_rule
            .fixed
            .iter()
            .copied()
            .map(usize::from)
            .find(|allowed| *allowed >= occupied_for_target)
        {
            target - occupied_for_target
        } else {
            // A trailing `-1` accepts any valence above the largest fixed one
            // but never implies additional hydrogens there.
            0
        };
        assignments.push((
            atom_id,
            u8::try_from(implicit).expect("RDKit implicit valences fit in u8"),
        ));
    }
    mol.install_valence(
        model_from_options(),
        assignments.into_iter().collect::<BTreeMap<_, _>>(),
    );
    ValenceReport { issues }
}

fn model_from_options() -> ValenceModel {
    ValenceModel::RdkitLike
}

fn record_unsupported(issues: &mut Vec<ValenceIssue>, atom: AtomId, options: ValenceOptions) {
    if options.strict {
        issues.push(ValenceIssue::UnsupportedElement(atom));
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AllowedValenceRule {
    fixed: &'static [u8],
    unrestricted_above: bool,
}

impl AllowedValenceRule {
    const fn fixed(fixed: &'static [u8]) -> Self {
        Self {
            fixed,
            unrestricted_above: false,
        }
    }

    const fn with_unrestricted(fixed: &'static [u8]) -> Self {
        Self {
            fixed,
            unrestricted_above: true,
        }
    }

    fn is_only_unrestricted(self) -> bool {
        self.fixed.is_empty() && self.unrestricted_above
    }

    fn max_fixed(self) -> Option<usize> {
        self.fixed.last().copied().map(usize::from)
    }
}

fn rdkit_effective_atomic_number(atom: &Atom) -> Option<u8> {
    let effective = i16::from(atom.element.atomic_number()) - i16::from(atom.formal_charge);
    (0..=118)
        .contains(&effective)
        .then(|| u8::try_from(effective).expect("range checked"))
}

fn can_be_rdkit_hypervalent_anion(atom: &Atom, effective_atomic_number: u8) -> bool {
    match atom.element.atomic_number() {
        15 | 16 => effective_atomic_number > 16,
        33 | 34 => effective_atomic_number > 34,
        _ => false,
    }
}

fn aromatic_valence_target(atom: &Atom, aromatic: Option<bool>, explicit: usize) -> Option<usize> {
    if aromatic != Some(true) {
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

pub(crate) fn explicit_valence(mol: &Molecule, atom: AtomId) -> usize {
    mol.incident_bonds(atom)
        .ok()
        .into_iter()
        .flatten()
        .map(|(_, bond)| bond_order_valence(bond.order))
        .sum()
}

fn bond_order_valence(order: BondOrder) -> usize {
    match order {
        BondOrder::Zero | BondOrder::Dative => 0,
        BondOrder::Single | BondOrder::Aromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
    }
}

pub(crate) fn allowed_valences(atom: &Atom) -> Option<&'static [u8]> {
    let original = rdkit_neutral_valence_rule(atom.element.atomic_number())?;
    let rule = if original.is_only_unrestricted() {
        original
    } else {
        rdkit_effective_atomic_number(atom).and_then(rdkit_neutral_valence_rule)?
    };
    Some(rule.fixed)
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

fn rdkit_neutral_valence_rule(atomic_number: u8) -> Option<AllowedValenceRule> {
    match atomic_number {
        0 => Some(AllowedValenceRule::with_unrestricted(&[])),
        1 => Some(AllowedValenceRule::fixed(&[1])),
        2 | 10 | 18 | 36 | 86 => Some(AllowedValenceRule::fixed(&[0])),
        3 | 11 | 19 | 37 => Some(AllowedValenceRule::with_unrestricted(&[1])),
        4 => Some(AllowedValenceRule::fixed(&[2])),
        12 | 20 | 38 | 56 | 88 => Some(AllowedValenceRule::with_unrestricted(&[2])),
        5 | 7 | 13 | 31 | 49 => Some(AllowedValenceRule::fixed(&[3])),
        6 | 14 | 32 => Some(AllowedValenceRule::fixed(&[4])),
        8 => Some(AllowedValenceRule::fixed(&[2])),
        9 | 17 | 35 => Some(AllowedValenceRule::fixed(&[1])),
        15 | 33 | 51 | 83 => Some(AllowedValenceRule::fixed(&[3, 5])),
        16 | 34 | 52 | 84 => Some(AllowedValenceRule::fixed(&[2, 4, 6])),
        50 | 82 => Some(AllowedValenceRule::fixed(&[2, 4])),
        53 | 85 => Some(AllowedValenceRule::fixed(&[1, 3, 5])),
        54 => Some(AllowedValenceRule::fixed(&[0, 2, 4, 6])),
        55 | 87 => Some(AllowedValenceRule::fixed(&[1])),
        // Every other current RDKit periodic-table entry has only `-1`.
        21..=30 | 39..=48 | 57..=81 | 89..=118 => Some(AllowedValenceRule::with_unrestricted(&[])),
        _ => None,
    }
}
