use crate::core::*;
use std::collections::BTreeMap;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ValenceOptions {
    /// Reject valence and implicit-hydrogen states that RDKit's strict property
    /// cache calculation rejects. With `false`, excess valence implies zero
    /// additional hydrogens instead of an error.
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
    InvalidFormalCharge {
        atom: AtomId,
        formal_charge: i8,
    },
    ValenceExceeded {
        atom: AtomId,
        explicit_valence: usize,
        max_allowed: usize,
    },
    /// Represented valence plus radical electrons and/or the hypervalent-anion
    /// charge offset exceeds the unsubtracted target valence.
    ValenceOccupancyExceeded {
        atom: AtomId,
        explicit_valence: usize,
        radical_electrons: usize,
        charge_offset: usize,
        max_allowed: usize,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceError {
    pub issues: Vec<ValenceIssue>,
}

impl fmt::Display for ValenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "valence perception reported {} issue(s)",
            self.issues.len()
        )
    }
}

impl std::error::Error for ValenceError {}

/// Installs implicit-hydrogen assignments from localized bonds and represented
/// atom state using the selected model's strict valence rules.
///
/// Success replaces installed valence and clears dependent aromaticity and CIP
/// assignments while retaining ring perception. Failure preserves all previously
/// installed perception. Neither outcome changes represented graph chemistry.
pub fn perceive_valence(
    mol: &mut Molecule,
    model: ValenceModel,
) -> std::result::Result<(), ValenceError> {
    perceive_valence_with_options(mol, model, ValenceOptions::default())
}

/// Installs valence with explicit control over strict validation.
///
/// With `options.strict == false`, RDKit's permissive property-cache rules assign
/// zero additional hydrogens to excess occupancy. As with [`perceive_valence`],
/// installation is transactional and invalidates aromaticity and CIP only on
/// success; represented chemistry and ring perception are preserved.
pub fn perceive_valence_with_options(
    mol: &mut Molecule,
    model: ValenceModel,
    options: ValenceOptions,
) -> std::result::Result<(), ValenceError> {
    match model {
        ValenceModel::RdkitLike => perceive_rdkit_like_valence(mol, options),
    }
}

fn perceive_rdkit_like_valence(
    mol: &mut Molecule,
    options: ValenceOptions,
) -> std::result::Result<(), ValenceError> {
    let mut assignments = Vec::<(AtomId, u8)>::new();
    let mut issues = Vec::new();
    for (atom_id, atom) in mol.atoms() {
        let implicit = match rdkit_atom_implicit_hydrogen_count(mol, atom_id, atom, options.strict)
        {
            Ok(implicit) => implicit,
            Err(issue) => {
                if options.strict {
                    issues.push(issue);
                }
                0
            }
        };
        assignments.push((atom_id, implicit));
    }
    if !issues.is_empty() {
        return Err(ValenceError { issues });
    }
    mol.install_valence(
        ValenceModel::RdkitLike,
        assignments.into_iter().collect::<BTreeMap<_, _>>(),
    );
    Ok(())
}

/// Derive an uninstalled assignment with the same permissive behavior as
/// non-strict valence perception. Aromaticity uses this only when the atom has
/// no installed hydrogen assignment.
pub(crate) fn rdkit_implicit_hydrogen_count(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> u8 {
    rdkit_atom_implicit_hydrogen_count(mol, atom_id, atom, false).unwrap_or(0)
}

fn rdkit_atom_implicit_hydrogen_count(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    strict: bool,
) -> std::result::Result<u8, ValenceIssue> {
    let explicit = explicit_valence(mol, atom_id) + usize::from(atom.hydrogens.explicit_count());
    let radical_electrons = atom
        .radical
        .map_or(0, |radical| usize::from(radical.unpaired_electron_count()));

    let original_rule = rdkit_neutral_valence_rule(atom.element.atomic_number())
        .ok_or(ValenceIssue::UnsupportedElement(atom_id))?;

    // RDKit leaves atoms whose periodic-table entry is only `-1` on that
    // unrestricted rule. All other charged atoms use the isoelectronic
    // neutral element's valence list.
    let effective_atomic_number = if original_rule.is_only_unrestricted() {
        atom.element.atomic_number()
    } else {
        rdkit_effective_atomic_number(atom)
    };
    let effective_rule = rdkit_neutral_valence_rule(effective_atomic_number)
        .ok_or(ValenceIssue::UnsupportedElement(atom_id))?;

    let hypervalent_anion = can_be_rdkit_hypervalent_anion(atom, effective_atomic_number);
    let charge_offset = if hypervalent_anion {
        usize::from(atom.formal_charge.unsigned_abs())
    } else {
        0
    };

    // Explicit-valence checking is separate from implicit-H inference. It
    // counts represented H and bonds, but not radical electrons. An unrestricted
    // sentinel on either valence list disables this check, except that
    // hypervalent anions use the original element's charge-adjusted limit.
    let two_coordinate_hydride = atom.element.atomic_number() == 1 && atom.formal_charge == -1;
    let explicit_limit = if original_rule.unrestricted_above {
        None
    } else if two_coordinate_hydride {
        // RDKit retains historical acceptance of two-coordinate hydride.
        Some(2)
    } else if hypervalent_anion {
        original_rule.max_fixed()
    } else if effective_rule.unrestricted_above {
        None
    } else {
        effective_rule.max_fixed()
    };
    if let Some(maximum) =
        explicit_limit.filter(|maximum| strict && explicit + charge_offset > *maximum)
    {
        if charge_offset != 0 {
            return Err(ValenceIssue::ValenceOccupancyExceeded {
                atom: atom_id,
                explicit_valence: explicit,
                radical_electrons: 0,
                charge_offset,
                max_allowed: maximum,
            });
        }
        return Err(ValenceIssue::ValenceExceeded {
            atom: atom_id,
            explicit_valence: explicit,
            max_allowed: maximum,
        });
    }

    // RDKit skips the implicit-valence calculation completely when H inference
    // is disabled, including its radical occupancy check.
    if !atom.hydrogens.allows_implicit() {
        return Ok(0);
    }
    if atom.element.atomic_number() == 1 && explicit == 0 && radical_electrons == 0 {
        return match atom.formal_charge {
            0 => Ok(1),
            -1 | 1 => Ok(0),
            _ if strict => Err(ValenceIssue::InvalidFormalCharge {
                atom: atom_id,
                formal_charge: atom.formal_charge,
            }),
            _ => Ok(0),
        };
    }
    // This exit precedes the hypervalent-anion adjustment in RDKit. A charge
    // mapping to an unrestricted element never implies hydrogen.
    if effective_rule.is_only_unrestricted() {
        return Ok(0);
    }
    let target_rule = if hypervalent_anion {
        original_rule
    } else {
        effective_rule
    };
    let occupied_for_target = explicit + radical_electrons + charge_offset;
    if let Some(target) = target_rule
        .fixed
        .iter()
        .copied()
        .map(usize::from)
        .find(|allowed| *allowed >= occupied_for_target)
    {
        return Ok(
            u8::try_from(target - occupied_for_target).expect("RDKit implicit valences fit in u8")
        );
    }
    // The implicit occupancy check has a different original-element guard
    // than explicit valence: an original zero-only valence list or unrestricted
    // sentinel does not reject excess radical occupancy.
    if strict
        && !target_rule.unrestricted_above
        && !original_rule.unrestricted_above
        && original_rule.max_fixed().is_some_and(|maximum| maximum > 0)
    {
        return Err(ValenceIssue::ValenceOccupancyExceeded {
            atom: atom_id,
            explicit_valence: explicit,
            radical_electrons,
            charge_offset,
            max_allowed: target_rule.max_fixed().expect("fixed target rule"),
        });
    }
    Ok(0)
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

fn rdkit_effective_atomic_number(atom: &Atom) -> u8 {
    let effective = i16::from(atom.element.atomic_number()) - i16::from(atom.formal_charge);
    // Atom::UpdatePropertyCache clamps charge adjustment even in strict mode.
    u8::try_from(effective.clamp(0, 118)).expect("clamped to periodic-table range")
}

fn can_be_rdkit_hypervalent_anion(atom: &Atom, effective_atomic_number: u8) -> bool {
    match atom.element.atomic_number() {
        15 | 16 => effective_atomic_number > 16,
        33 | 34 => effective_atomic_number > 34,
        _ => false,
    }
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
        BondOrder::Single => 1,
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
        rdkit_neutral_valence_rule(rdkit_effective_atomic_number(atom))?
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
    rdkit_neutral_valence_rule(atomic_number)?
        .fixed
        .first()
        .copied()
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
