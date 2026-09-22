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

/// Installs inferred-hydrogen assignments from localized bonds and represented
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
    let assignments = match model {
        ValenceModel::RdkitLike => rdkit_valence_assignments(mol, options)?,
    };
    mol.install_valence(model, assignments);
    Ok(())
}

pub(crate) fn rdkit_valence_assignments(
    mol: &Molecule,
    options: ValenceOptions,
) -> std::result::Result<BTreeMap<AtomId, u8>, ValenceError> {
    let mut assignments = BTreeMap::new();
    let mut issues = Vec::new();
    for (atom_id, atom) in mol.atoms() {
        let implicit = match rdkit_atom_inferred_hydrogen_count(mol, atom_id, atom, options.strict)
        {
            Ok(implicit) => implicit,
            Err(issue) => {
                if options.strict {
                    issues.push(issue);
                }
                0
            }
        };
        assignments.insert(atom_id, implicit);
    }
    if !issues.is_empty() {
        return Err(ValenceError { issues });
    }
    Ok(assignments)
}

/// Derive an uninstalled assignment with the same permissive behavior as
/// non-strict valence perception. Aromaticity uses this only when the atom has
/// no installed hydrogen assignment.
pub(crate) fn rdkit_inferred_hydrogen_count(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> u8 {
    rdkit_atom_inferred_hydrogen_count(mol, atom_id, atom, false).unwrap_or(0)
}

/// Radical-electron count implied by a fixed-H SMILES bracket atom after
/// localization. Mirrors RDKit's octet/duet and allowed-valence conventions,
/// without assigning a spin multiplicity or changing a published graph.
/// `hydrogens` is the count that the bracket actually encodes, not a storage policy.
pub(crate) fn rdkit_bracket_radical_electrons(
    atom: &Atom,
    bond_valence: usize,
    degree: usize,
    hydrogens: u8,
) -> u8 {
    // The outer-electron entries of RDKit 2026.03.3's periodic table. These
    // belong to this valence model; they are not electronic ground states.
    const OUTER_ELECTRONS: [u8; 119] = [
        0, 1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3,
        4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1,
        2, 3, 4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2,
    ];
    let number = atom.element.atomic_number();
    let rule = rdkit_neutral_valence_rule(number)
        .expect("supported elements have RDKit-like valence rules");
    let outer = i64::from(OUTER_ELECTRONS[usize::from(number)]);
    let charge = i64::from(atom.formal_charge);
    if rule.is_only_unrestricted() {
        // With no preferred valence, bonded metal centers do not acquire
        // invented radicals. Isolated atoms use only electron-count parity.
        return if degree == 0 {
            ((outer - charge).max(0) % 2) as u8
        } else {
            0
        };
    }
    let Ok(occupied) = i32::try_from(bond_valence.saturating_add(usize::from(hydrogens))) else {
        // Arbitrarily large represented valence cannot leave an octet deficit.
        return 0;
    };
    let occupied = i64::from(occupied);
    let shell = if number <= 2 { 2 } else { 8 };
    let mut electrons = shell - outer - occupied + charge;
    if electrons < 0 {
        electrons = 0;
        if rule.fixed.len() > 1 {
            electrons = rule
                .fixed
                .iter()
                .map(|valence| i64::from(*valence) - occupied + charge)
                .find(|deficit| *deficit >= 0)
                .unwrap_or(0);
        }
    }
    let available = outer - occupied - charge;
    if available >= 0 {
        electrons = electrons.min(available);
    }
    u8::try_from(electrons).expect("octet deficits with i8 charges fit in u8")
}

fn rdkit_atom_inferred_hydrogen_count(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    strict: bool,
) -> std::result::Result<u8, ValenceIssue> {
    let explicit = explicit_valence(mol, atom_id) + usize::from(atom.hydrogens.specified_count());
    let radical_electrons = atom
        .radical
        .map_or(0, |radical| usize::from(radical.electron_count()));

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
    if !atom.hydrogens.allows_inference() {
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

/// Returns the represented valence used by this module's hydrogen inference.
///
/// Counts localized incident bond orders and the atom's declared hydrogens.
/// Inferred hydrogens are excluded. Zero-order and dative bonds contribute zero
/// at both endpoints, as specified by [`ValenceModel::RdkitLike`]. This reads
/// represented chemistry without installing or changing perception.
///
/// Returns an error if `atom` is not live in `mol`.
pub fn represented_valence(mol: &Molecule, atom: AtomId) -> Result<usize> {
    let declaration = mol.atom(atom)?.hydrogens.specified_count();
    Ok(explicit_valence(mol, atom) + usize::from(declaration))
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
