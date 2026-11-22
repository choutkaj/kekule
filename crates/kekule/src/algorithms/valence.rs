use crate::diagnostics::Capitalized;
use crate::molecule::*;
use std::collections::BTreeMap;
use std::fmt;

/// Options of [`perceive_valence`]: how invalid valence states are treated.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ValenceOptions {
    validation: ValenceValidation,
}

impl ValenceOptions {
    /// How invalid valence states are treated; strict by default.
    pub const fn validation(&self) -> ValenceValidation {
        self.validation
    }

    /// Sets [`Self::validation`].
    #[must_use]
    pub const fn with_validation(mut self, validation: ValenceValidation) -> Self {
        self.validation = validation;
        self
    }

    const fn strict(&self) -> bool {
        match self.validation {
            ValenceValidation::Strict => true,
            ValenceValidation::Permissive => false,
        }
    }
}

/// One atom's invalid valence state.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ValenceViolation {
    /// In strict mode, an unbonded hydrogen with inferred hydrogens has a
    /// formal charge other than -1, 0 or +1.
    InvalidFormalCharge {
        /// The hydrogen atom.
        atom: AtomId,
        /// Its formal charge.
        formal_charge: i8,
    },
    /// In strict mode, the represented valence exceeds the element's maximum.
    ValenceExceeded {
        /// The offending atom.
        atom: AtomId,
        /// Bond orders plus explicit hydrogens.
        explicit_valence: usize,
        /// The largest allowed explicit valence.
        max_allowed: usize,
    },
    /// Represented valence plus radical electrons and/or the hypervalent-anion
    /// charge offset exceeds the unsubtracted target valence.
    ValenceOccupancyExceeded {
        /// The offending atom.
        atom: AtomId,
        /// Bond orders plus explicit hydrogens.
        explicit_valence: usize,
        /// Radical electrons counted against the target valence.
        radical_electrons: usize,
        /// The hypervalent-anion charge offset, or zero.
        charge_offset: usize,
        /// The applicable maximum valence.
        max_allowed: usize,
    },
}

/// Valence perception found atoms whose valence state is invalid.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceError {
    /// One violation per offending atom.
    pub violations: Vec<ValenceViolation>,
}

impl ValenceViolation {
    /// The offending atom.
    pub(crate) const fn atom(&self) -> AtomId {
        match self {
            Self::InvalidFormalCharge { atom, .. }
            | Self::ValenceExceeded { atom, .. }
            | Self::ValenceOccupancyExceeded { atom, .. } => *atom,
        }
    }

    /// The same violation of `atom`, as when publication renumbers atoms.
    pub(crate) fn with_atom(mut self, atom: AtomId) -> Self {
        match &mut self {
            Self::InvalidFormalCharge { atom: old, .. }
            | Self::ValenceExceeded { atom: old, .. }
            | Self::ValenceOccupancyExceeded { atom: old, .. } => *old = atom,
        }
        self
    }

    /// Writes the violation as a sentence about `subject`, a capitalized
    /// noun phrase for the atom, such as `The atom at line 5`.
    pub(crate) fn write_about(
        &self,
        formatter: &mut fmt::Formatter<'_>,
        subject: impl fmt::Display,
    ) -> fmt::Result {
        match self {
            Self::InvalidFormalCharge { formal_charge, .. } => write!(
                formatter,
                "{subject} is an unbonded hydrogen with formal charge {formal_charge:+}."
            ),
            Self::ValenceExceeded {
                explicit_valence,
                max_allowed,
                ..
            } => write!(
                formatter,
                "{subject} has valence {explicit_valence}, more than the allowed {max_allowed}."
            ),
            Self::ValenceOccupancyExceeded {
                explicit_valence,
                radical_electrons,
                charge_offset,
                max_allowed,
                ..
            } => write!(
                formatter,
                "{subject} has valence {explicit_valence} plus {} for its radical electrons and charge, more than the allowed {max_allowed}.",
                radical_electrons + charge_offset
            ),
        }
    }
}

impl fmt::Display for ValenceViolation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_about(formatter, Capitalized(self.atom()))
    }
}

/// Writes the first violation, then the number of atoms with violations when there
/// are several.
impl fmt::Display for ValenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.violations.as_slice() {
            [] => formatter.write_str("The valence state is invalid."),
            [issue] => write!(formatter, "{issue}"),
            [issue, ..] => write!(
                formatter,
                "{issue} In all, {} atoms have an invalid valence state.",
                self.violations.len()
            ),
        }
    }
}

impl std::error::Error for ValenceError {}

/// Installs inferred-hydrogen assignments from localized bonds and represented
/// atom state using the selected model's valence rules.
///
/// The default options apply strict validation. With
/// [`ValenceValidation::Permissive`], RDKit's permissive property-cache rules
/// assign zero additional hydrogens to excess occupancy. Success replaces
/// installed valence and clears dependent aromaticity and CIP assignments
/// while retaining ring perception. Failure preserves all previously installed
/// perception. Neither outcome changes represented graph chemistry.
///
/// # Errors
///
/// Only strict validation fails: [`ValenceError`] with one [`ValenceViolation`]
/// per atom whose valence or hydrogen state is invalid.
pub fn perceive_valence(
    mol: &mut Molecule,
    model: ValenceModel,
    options: ValenceOptions,
) -> Result<(), ValenceError> {
    let assignments = match model {
        ValenceModel::RdkitLike => rdkit_valence_assignments(mol, options)?,
    };
    mol.install_valence(model, options.validation(), assignments);
    Ok(())
}

pub(crate) fn rdkit_valence_assignments(
    mol: &Molecule,
    options: ValenceOptions,
) -> Result<BTreeMap<AtomId, u8>, ValenceError> {
    let mut assignments = BTreeMap::new();
    let mut issues = Vec::new();
    for (atom_id, atom) in mol.atoms() {
        let implicit =
            match rdkit_atom_inferred_hydrogen_count(mol, atom_id, atom, options.strict()) {
                Ok(implicit) => implicit,
                Err(issue) => {
                    if options.strict() {
                        issues.push(issue);
                    }
                    0
                }
            };
        assignments.insert(atom_id, implicit);
    }
    if !issues.is_empty() {
        return Err(ValenceError { violations: issues });
    }
    Ok(assignments)
}

/// Derive an uninstalled assignment with the same permissive behavior as
/// non-strict valence perception. Aromaticity uses this only when the atom has
/// no installed hydrogen assignment, and hydrogen collapse uses it to decide
/// whether inference reproduces a collapsed count.
pub(crate) fn rdkit_inferred_hydrogen_count<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    atom: &Atom,
) -> u8 {
    rdkit_atom_inferred_hydrogen_count(mol, atom_id, atom, false).unwrap_or(0)
}

/// The violation strict valence perception would report for one atom, if any.
pub(crate) fn rdkit_valence_violation<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    atom: &Atom,
) -> Option<ValenceViolation> {
    rdkit_atom_inferred_hydrogen_count(mol, atom_id, atom, true).err()
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

fn rdkit_atom_inferred_hydrogen_count<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    atom: &Atom,
    strict: bool,
) -> Result<u8, ValenceViolation> {
    let explicit = explicit_valence(mol, atom_id) + usize::from(atom.hydrogens.represented_count());
    let bonded = mol
        .incident_bonds(atom_id)
        .is_some_and(|mut bonds| bonds.next().is_some());
    hydrogens_at_valence(atom_id, atom, explicit, bonded, strict, Hydrides::Model)
}

/// Where hydrogens are inferred: this model's rule or RDKit's.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Hydrides {
    /// No hydrides on main-group metals and no hydrogen on a bonded hydrogen
    /// atom.
    Model,
    /// RDKit's rule, which a Molfile reader also applies.
    Rdkit,
}

/// The hydrogens of `atom` at represented valence `explicit`. `bonded` says
/// whether it has bonds of any order.
fn hydrogens_at_valence(
    atom_id: AtomId,
    atom: &Atom,
    explicit: usize,
    bonded: bool,
    strict: bool,
    hydrides: Hydrides,
) -> Result<u8, ValenceViolation> {
    let radical_electrons = atom
        .radical
        .map_or(0, |radical| usize::from(radical.electron_count()));

    let original_rule = rdkit_neutral_valence_rule(atom.element.atomic_number())
        .expect("RDKit valence rules cover atomic numbers 0 to 118");

    // RDKit leaves atoms whose periodic-table entry is only `-1` on that
    // unrestricted rule. All other charged atoms use the isoelectronic
    // neutral element's valence list.
    let effective_atomic_number = if original_rule.is_only_unrestricted() {
        atom.element.atomic_number()
    } else {
        rdkit_effective_atomic_number(atom)
    };
    let effective_rule = rdkit_neutral_valence_rule(effective_atomic_number)
        .expect("RDKit valence rules cover atomic numbers 0 to 118");

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
            return Err(ValenceViolation::ValenceOccupancyExceeded {
                atom: atom_id,
                explicit_valence: explicit,
                radical_electrons: 0,
                charge_offset,
                max_allowed: maximum,
            });
        }
        return Err(ValenceViolation::ValenceExceeded {
            atom: atom_id,
            explicit_valence: explicit,
            max_allowed: maximum,
        });
    }

    // RDKit skips the implicit-valence calculation completely when H inference
    // is disabled, including its radical occupancy check.
    if !atom.hydrogens.is_inferred() {
        return Ok(0);
    }
    // Unlike RDKit, no hydride is inferred on these metals. An unbonded or
    // under-bonded one more likely lacks its charge or its bonds than carries
    // hydrogens; its valence list stays a limit.
    if hydrides == Hydrides::Model && is_main_group_metal(atom.element.atomic_number()) {
        return Ok(0);
    }
    if atom.element.atomic_number() == 1 && explicit == 0 && radical_electrons == 0 {
        return match atom.formal_charge {
            // A hydrogen atom joined only by zero-order or dative bonds, as
            // a hydride on a metal, takes no hydrogen of its own.
            0 if bonded && hydrides == Hydrides::Model => Ok(0),
            0 => Ok(1),
            -1 | 1 => Ok(0),
            _ if strict => Err(ValenceViolation::InvalidFormalCharge {
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
        return Err(ValenceViolation::ValenceOccupancyExceeded {
            atom: atom_id,
            explicit_valence: explicit,
            radical_electrons,
            charge_offset,
            max_allowed: target_rule.max_fixed().expect("fixed target rule"),
        });
    }
    Ok(0)
}

/// Whether atoms of this element are metals: alkali, alkaline-earth,
/// transition, post-transition, lanthanide, and actinide metals. Metalloids
/// are not metals.
pub(crate) fn is_metal(element: Element) -> bool {
    matches!(
        element.atomic_number(),
        3 | 4 | 11..=13 | 19..=31 | 37..=50 | 55..=84 | 87..=116
    )
}

/// Whether atoms of this element are alkali, alkaline-earth, or p-block
/// metals, on which [`ValenceModel::RdkitLike`] infers no hydrogens.
pub(crate) const fn is_main_group_metal(atomic_number: u8) -> bool {
    matches!(
        atomic_number,
        3 | 4 | 11..=13 | 19 | 20 | 31 | 37 | 38 | 49 | 50 | 55 | 56 | 81..=84 | 87 | 88
    )
}

/// The hydrogens RDKit's rule, which Molfile readers share, infers on an
/// atom with inferred hydrogens. It differs from this model on main-group
/// metals and on hydrogen atoms joined only by zero-order or dative bonds.
pub(crate) fn rdkit_reader_hydrogens<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    atom: &Atom,
) -> u8 {
    let explicit = explicit_valence(mol, atom_id) + usize::from(atom.hydrogens.represented_count());
    hydrogens_at_valence(atom_id, atom, explicit, true, false, Hydrides::Rdkit).unwrap_or(0)
}

/// The hydrogens the model infers on a neutral onium centre: an atom with
/// more bonded neighbors than its default valence, which takes hydrogens only
/// at a higher allowed valence, as the sulfur of `CS(C)C` takes one at
/// valence four. `None` for every other atom.
///
/// OpenSMILES specifies these hydrogens for its organic subset. Formats that
/// state charges, such as Molfile and mmCIF component tables, more likely
/// omit charges than draw such a hydride ([`implied_charges`]).
pub(crate) fn onium_hydrogens<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    atom: &Atom,
) -> Option<u8> {
    if atom.formal_charge != 0 || atom.radical.is_some() || !atom.hydrogens.is_inferred() {
        return None;
    }
    let default = usize::from(rdkit_default_valence(atom)?);
    let neighbors = mol
        .incident_bonds(atom_id)?
        .filter(|(_, bond)| bond_order_valence(bond.order) > 0)
        .count();
    if neighbors <= default {
        return None;
    }
    rdkit_atom_inferred_hydrogen_count(mol, atom_id, atom, true)
        .ok()
        .filter(|hydrogens| *hydrogens > 0)
}

/// The charges a neutral atom's bonds imply when its source states charges
/// but none for it ([`implied_charges`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ImpliedCharges {
    /// An onium centre ([`onium_hydrogens`]) takes +1, as the sulfur of
    /// trimethylsulfonium.
    Onium,
    /// An atom with more bonds than any neutral valence allows takes the
    /// first of +1 and -1 that makes its valence valid, as a quaternary
    /// nitrogen takes +1 and a tetrahedral boron -1.
    OverBonded(i8),
    /// The atom takes +1 and its singly bonded neighbor -1. The neighbor is
    /// a chalcogen with no other non-metal neighbor, as the oxygen of DMSO
    /// or trimethylamine N-oxide drawn with a single bond, or a boron or
    /// aluminium acceptor, as the boron of a phosphine-borane.
    Separated(AtomId),
}

/// The charges a neutral atom with inferred hydrogens takes when its source
/// states charges but none for it: an onium centre ([`onium_hydrogens`]) or
/// an atom with more bonds than any neutral valence allows. `None` for every
/// other atom, and for a metal donor ([`is_metal_donor`]). `may_charge` says
/// which neighbors may take -1.
pub(crate) fn implied_charges<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    atom: &Atom,
    may_charge: impl Fn(AtomId) -> bool,
) -> Option<ImpliedCharges> {
    if atom.formal_charge != 0
        || atom.radical.is_some()
        || !atom.hydrogens.is_inferred()
        || is_metal_donor(mol, atom_id, atom)
    {
        return None;
    }
    let over_bonded = matches!(
        rdkit_valence_violation(mol, atom_id, atom),
        Some(ValenceViolation::ValenceExceeded { .. })
    );
    if !over_bonded && onium_hydrogens(mol, atom_id, atom).is_none() {
        return None;
    }
    if valid_charge(mol, atom_id, atom, 1).is_some() {
        if let Some(partner) = charge_partner(mol, atom_id, may_charge) {
            return Some(ImpliedCharges::Separated(partner));
        }
    }
    if over_bonded {
        return [1, -1]
            .into_iter()
            .find_map(|charge| valid_charge(mol, atom_id, atom, charge))
            .map(ImpliedCharges::OverBonded);
    }
    valid_charge(mol, atom_id, atom, 1).map(|_| ImpliedCharges::Onium)
}

/// The singly bonded neighbor of `atom_id` that takes -1 with its +1
/// ([`ImpliedCharges::Separated`]): a neutral chalcogen whose only non-metal
/// neighbor is `atom_id`, or else a neutral boron or aluminium atom. Its bonds
/// to metals may be donor bonds ([`is_metal_donor`]).
fn charge_partner<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    may_charge: impl Fn(AtomId) -> bool,
) -> Option<AtomId> {
    let (mut chalcogen, mut acceptor) = (None, None);
    for (_, bond) in mol.incident_bonds(atom_id)? {
        let neighbor_id = bond.other_atom(atom_id);
        let neighbor = mol.atom(neighbor_id)?;
        if bond.order != BondOrder::Single
            || neighbor.formal_charge != 0
            || neighbor.radical.is_some()
            || !neighbor.hydrogens.is_inferred()
            || !may_charge(neighbor_id)
        {
            continue;
        }
        let mut anion = neighbor.clone();
        anion.formal_charge = -1;
        if rdkit_valence_violation(mol, neighbor_id, &anion).is_some()
            && !is_metal_donor(mol, neighbor_id, &anion)
        {
            continue;
        }
        match neighbor.element.atomic_number() {
            8 | 16 | 34 | 52 if non_metal_neighbors(mol, neighbor_id)? == 1 => {
                chalcogen = chalcogen.or(Some(neighbor_id));
            }
            5 | 13 => acceptor = acceptor.or(Some(neighbor_id)),
            _ => {}
        }
    }
    chalcogen.or(acceptor)
}

fn non_metal_neighbors<L: SlotLiveness>(mol: &Graph<L>, atom_id: AtomId) -> Option<usize> {
    let mut count = 0;
    for neighbor in mol.neighbors(atom_id)? {
        if !is_metal(mol.atom(neighbor)?.element) {
            count += 1;
        }
    }
    Some(count)
}

/// Whether a non-metal atom with inferred hydrogens donates to metals: its
/// bond valence is not one its element and charge allow, but some allowed
/// valence lies between it less its bond valence to metals and it, as for
/// the phosphorus of a phosphine on palladium, a sulfide bridging three iron
/// atoms, or the sulfur of a sulfoxide on ruthenium. Its bonds to metals are
/// then donor bonds drawn as covalent ones, and it takes no hydrogens.
pub(crate) fn is_metal_donor<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    atom: &Atom,
) -> bool {
    if is_metal(atom.element) || atom.radical.is_some() || !atom.hydrogens.is_inferred() {
        return false;
    }
    let Some(bonds) = mol.incident_bonds(atom_id) else {
        return false;
    };
    let mut metal_valence = 0;
    for (_, bond) in bonds {
        if mol
            .atom(bond.other_atom(atom_id))
            .is_some_and(|neighbor| is_metal(neighbor.element))
        {
            metal_valence += bond_order_valence(bond.order);
        }
    }
    let valence = explicit_valence(mol, atom_id);
    metal_valence > 0 && donates_to_metals(atom, valence - metal_valence, metal_valence)
}

/// Whether a non-metal donates to metals when its other bonds give it the
/// valence `other` and its bonds to metals add `metal`: `other + metal` is
/// not a valence its element and charge allow, but some allowed valence lies
/// from `other` up to, and not including, `other + metal`. This one rule
/// decides metal donors ([`is_metal_donor`]) and untyped metal bonds
/// ([`metal_bond_reading`]).
pub(crate) fn donates_to_metals(atom: &Atom, other: usize, metal: usize) -> bool {
    let Some(allowed) = allowed_valences(atom) else {
        return false;
    };
    let allowed = || allowed.iter().map(|value| usize::from(*value));
    !allowed().any(|value| value == other + metal)
        && allowed().any(|value| (other..other + metal).contains(&value))
}

/// How a reader takes a single bond between a metal and a non-metal of one
/// residue whose source states no type, as PDB `CONECT` does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MetalBondReading {
    Covalent,
    Dative,
    /// The source cannot tell: the bond is read as covalent and reported.
    Unstated,
}

/// The reading of an untyped metal bond from its non-metal `donor`.
/// `other` is the donor's valence from its other bonds, those of
/// `undecided` untyped metal bonds of the residue left out, and `hydrogens`
/// whether the residue has explicit hydrogens.
///
/// The bond is dative when the donor donates to the metal
/// ([`donates_to_metals`]), as an ammine nitrogen or a phosphine does. It is
/// covalent when the bond leaves the donor at an allowed valence, as for a
/// chloride or a phosphonate phosphorus. A residue without hydrogens, a
/// donor of several untyped metal bonds, such as a bridging chloride, and a
/// donor that neither reading satisfies, such as a carbene carbon, are
/// unstated.
pub(crate) fn metal_bond_reading(
    donor: &Atom,
    other: usize,
    undecided: usize,
    hydrogens: bool,
) -> MetalBondReading {
    if !hydrogens || undecided != 1 {
        return MetalBondReading::Unstated;
    }
    if donates_to_metals(donor, other, 1) {
        return MetalBondReading::Dative;
    }
    let covalent = allowed_valences(donor).is_some_and(|valences| {
        valences
            .iter()
            .any(|&allowed| usize::from(allowed) == other + 1)
    });
    if covalent {
        MetalBondReading::Covalent
    } else {
        MetalBondReading::Unstated
    }
}

fn valid_charge<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    atom: &Atom,
    charge: i8,
) -> Option<i8> {
    let mut charged = atom.clone();
    charged.formal_charge = charge;
    rdkit_valence_violation(mol, atom_id, &charged)
        .is_none()
        .then_some(charge)
}

/// The valence of a metal of [`is_main_group_metal`] with inferred
/// hydrogens, when its element and charge allow no such valence, as of a tin
/// atom with three bonds or an unbonded neutral sodium atom. The model infers
/// no hydrides to complete it, so the source likely omits a charge, a bond,
/// or a hydrogen. Radical electrons count toward the valence. A valence above
/// the allowed ones is a [`ValenceViolation`] instead, unless the element
/// allows any higher valence.
pub(crate) fn incomplete_metal_valence<L: SlotLiveness>(
    mol: &Graph<L>,
    atom_id: AtomId,
    atom: &Atom,
) -> Option<u8> {
    if !atom.hydrogens.is_inferred() || !is_main_group_metal(atom.element.atomic_number()) {
        return None;
    }
    let radical_electrons = atom
        .radical
        .map_or(0, |radical| usize::from(radical.electron_count()));
    let valence = explicit_valence(mol, atom_id) + radical_electrons;
    let allowed = allowed_valences(atom)?;
    let below_maximum = allowed
        .last()
        .is_some_and(|max| valence < usize::from(*max));
    let listed = allowed.iter().any(|value| usize::from(*value) == valence);
    (below_maximum && !listed).then(|| u8::try_from(valence).expect("below an allowed valence"))
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
    effective_atomic_number(atom.element.atomic_number(), atom.formal_charge)
}

fn effective_atomic_number(atomic_number: u8, charge: i8) -> u8 {
    let effective = i16::from(atomic_number) - i16::from(charge);
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
/// # Errors
///
/// [`MoleculeError::InvalidAtomId`] when `atom` is not live in `mol`.
pub fn represented_valence(mol: &Molecule, atom: AtomId) -> Result<usize, MoleculeError> {
    let declaration = mol.require_atom(atom)?.hydrogens.represented_count();
    Ok(explicit_valence(mol, atom) + usize::from(declaration))
}

pub(crate) fn explicit_valence<L: SlotLiveness>(mol: &Graph<L>, atom: AtomId) -> usize {
    mol.incident_bonds(atom)
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
    allowed_valences_with_charge(atom.element.atomic_number(), atom.formal_charge)
}

/// [`allowed_valences`] of an atom of `atomic_number` carrying `charge`.
pub(crate) fn allowed_valences_with_charge(atomic_number: u8, charge: i8) -> Option<&'static [u8]> {
    let original = rdkit_neutral_valence_rule(atomic_number)?;
    let rule = if original.is_only_unrestricted() {
        original
    } else {
        rdkit_neutral_valence_rule(effective_atomic_number(atomic_number, charge))?
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

#[cfg(test)]
mod tests;
