//! Resonance delocalization perception and explicit RDKit-like contributor
//! enumeration over connected conjugated groups.
//!
//! The enumeration search follows RDKit 2026.03.3 (Paolo Tosco, 2015), distributed under
//! the BSD 3-Clause license reproduced in LICENSE-RDKit in this crate.
use crate::diagnostics::Countable;
use crate::diagnostics::IndexOutOfRange;
use crate::diagnostics::{Cancellation, Cancelled, LimitExceeded};
use crate::molecule::{
    AtomId, BondId, BondOrder, ImplicitHydrogens, Molecule, MoleculeError,
    MoleculePublicationError, ResonanceGroup, ResonanceModel, ResonancePerception,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

mod delocalization;
mod kekule;
mod search;

/// RDKit resonance contributor policies. Combine constants with `|`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ResonanceFlags(u8);
impl ResonanceFlags {
    /// Keep contributors with less complete octets than the most
    /// octet-complete one.
    pub const ALLOW_INCOMPLETE_OCTETS: Self = Self(1);
    /// Keep charge-separated contributors even when uncharged ones exist.
    pub const ALLOW_CHARGE_SEPARATION: Self = Self(2);
    /// Enumerate every degenerate Kekulé structure instead of one.
    pub const KEKULE_ALL: Self = Self(4);
    /// Do not restrict which atoms may carry a positive charge; implies
    /// [`Self::ALLOW_INCOMPLETE_OCTETS`] and [`Self::ALLOW_CHARGE_SEPARATION`].
    pub const UNCONSTRAINED_CATIONS: Self = Self(8);
    /// Do not restrict which atoms may carry a negative charge; implies
    /// [`Self::ALLOW_CHARGE_SEPARATION`].
    pub const UNCONSTRAINED_ANIONS: Self = Self(16);
    /// Constructs flags, rejecting unknown bits.
    pub const fn from_bits(bits: u8) -> Option<Self> {
        if bits < 32 {
            Some(Self(bits))
        } else {
            None
        }
    }
    /// The raw flag bits, as accepted by [`Self::from_bits`].
    pub const fn bits(self) -> u8 {
        self.0
    }
    /// Whether every bit of `flag` is set.
    pub const fn contains(self, flag: Self) -> bool {
        self.0 & flag.0 == flag.0
    }
    /// Adds the flags implied by [`Self::UNCONSTRAINED_CATIONS`] and
    /// [`Self::UNCONSTRAINED_ANIONS`], as RDKit does.
    /// [`enumerate_resonance`] applies this itself.
    pub const fn effective(self) -> Self {
        Self(self.0 | if self.0 & 8 != 0 { 3 } else { 0 } | if self.0 & 16 != 0 { 2 } else { 0 })
    }
}
impl std::ops::BitOr for ResonanceFlags {
    type Output = Self;
    fn bitor(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }
}

/// Bounds and contributor policy for explicit enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ResonanceOptions {
    flags: ResonanceFlags,
    max_structures: usize,
    max_total_work: usize,
}

impl ResonanceOptions {
    /// Contributor policy; empty by default.
    pub const fn flags(&self) -> ResonanceFlags {
        self.flags
    }

    /// Sets [`Self::flags`].
    #[must_use]
    pub const fn with_flags(mut self, flags: ResonanceFlags) -> Self {
        self.flags = flags;
        self
    }

    /// Maximum returned contributors, clamped to RDKit's ceiling of 1,000,000.
    /// Zero returns an empty collection. The pinned reference cannot produce a
    /// valid conjugated contributor at a limit of one; use at least two.
    pub const fn max_structures(&self) -> usize {
        self.max_structures
    }

    /// Sets [`Self::max_structures`].
    #[must_use]
    pub const fn with_max_structures(mut self, max_structures: usize) -> Self {
        self.max_structures = max_structures;
        self
    }

    /// Bounds search visits, copied state, and permutation construction/sorting.
    /// Exhaustion returns an error without a partial result.
    pub const fn max_total_work(&self) -> usize {
        self.max_total_work
    }

    /// Sets [`Self::max_total_work`].
    #[must_use]
    pub const fn with_max_total_work(mut self, max_total_work: usize) -> Self {
        self.max_total_work = max_total_work;
        self
    }
}
impl ResonanceOptions {
    /// The [`Default`] bound of [`Self::max_structures`].
    pub const DEFAULT_MAX_STRUCTURES: usize = 1_000;
    /// The ceiling every [`Self::max_structures`] is clamped to, RDKit's.
    pub const MAX_STRUCTURES: usize = 1_000_000;
    /// The [`Default`] bound of [`Self::max_total_work`].
    pub const DEFAULT_MAX_TOTAL_WORK: usize = 100_000_000;
}

impl Default for ResonanceOptions {
    fn default() -> Self {
        Self {
            flags: ResonanceFlags::default(),
            max_structures: Self::DEFAULT_MAX_STRUCTURES,
            max_total_work: Self::DEFAULT_MAX_TOTAL_WORK,
        }
    }
}

/// Why resonance perception failed.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ResonancePerceptionError {
    /// Aromaticity is not installed.
    MissingAromaticity,
    /// The implicit hydrogen count of an atom the search reads is unknown.
    UnknownHydrogens(AtomId),
    /// The path search exceeded its work bound, which grows with the graph.
    LimitExceeded(LimitExceeded),
    /// The caller's [`Cancellation`] stopped the path search. Only
    /// cancellable perception, such as
    /// [`Molecule::perceive_with`](crate::molecule::Molecule::perceive_with),
    /// reports it, as [`PerceptionError::Cancelled`](crate::perception::PerceptionError::Cancelled);
    /// [`perceive_resonance`] never does.
    Cancelled(Cancelled),
}
impl fmt::Display for ResonancePerceptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingAromaticity => {
                f.write_str("Aromaticity is not perceived. Perceive it before resonance.")
            }
            Self::UnknownHydrogens(a) => write!(
                f,
                "The hydrogen count of {a} is unknown. Perceive valence first."
            ),
            Self::LimitExceeded(limit) => {
                write!(f, "Resonance perception exceeded its {}.", limit.describe())
            }
            Self::Cancelled(_) => f.write_str("Resonance perception stopped."),
        }
    }
}
impl std::error::Error for ResonancePerceptionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cancelled(cancelled) => Some(cancelled),
            _ => None,
        }
    }
}

/// Why resonance enumeration failed.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ResonanceError {
    /// Conjugation is not installed.
    MissingConjugation,
    /// The atom's implicit hydrogen count is unknown.
    UnknownHydrogens(AtomId),
    /// A resonance group contains a bond the search does not support.
    UnsupportedBond(BondId),
    /// The search reached an impossible electron count.
    InvalidElectronCount,
    /// [`ResonanceOptions::max_structures`] is one; RDKit-compatible
    /// enumeration needs zero or at least two.
    InvalidStructureLimit,
    /// The enumeration exceeded [`ResonanceOptions::max_total_work`].
    LimitExceeded(LimitExceeded),
    /// No contributor has this index.
    ContributorOutOfRange(IndexOutOfRange),
    /// Editing the contributor molecule failed; see the source.
    Molecule(MoleculeError),
    /// Publishing the contributor molecule failed; see the source.
    Publication(MoleculePublicationError),
}
impl fmt::Display for ResonanceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingConjugation => {
                f.write_str("Conjugation is not perceived. Perceive it before resonance.")
            }
            Self::UnknownHydrogens(a) => write!(
                f,
                "The hydrogen count of {a} is unknown. Perceive valence first."
            ),
            Self::UnsupportedBond(b) => {
                write!(f, "Resonance enumeration does not support {b}.")
            }
            Self::InvalidElectronCount => {
                f.write_str("Resonance enumeration reached an impossible electron count.")
            }
            Self::InvalidStructureLimit => {
                f.write_str("The resonance structure limit is 1. Use 0 or at least 2.")
            }
            Self::LimitExceeded(limit) => {
                write!(
                    f,
                    "Resonance enumeration exceeded its {}.",
                    limit.describe()
                )
            }
            Self::ContributorOutOfRange(error) => error.fmt(f),
            Self::Molecule(_) => f.write_str("Cannot edit the resonance contributor."),
            Self::Publication(_) => f.write_str("Cannot publish the resonance contributor."),
        }
    }
}
impl std::error::Error for ResonanceError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Molecule(error) => Some(error),
            Self::Publication(error) => Some(error),
            _ => None,
        }
    }
}

/// Installs resonance delocalization: which bonds resonance delocalizes, in
/// connected groups.
///
/// [`ResonanceModel::ChargeMoves`] moves each mobile formal charge of ±1 along
/// simple paths of alternately raised and lowered bond orders to an atom that
/// can carry the same kind of charge, and marks every bond of such a path.
/// Aromatic bonds are left to aromaticity, so groups never contain them.
/// Neutral push-pull systems such as amides stay localized. See
/// `docs/resonance.md` for the rules. No contributor is enumerated. Only
/// charged atoms start a search, so a neutral molecule costs one scan of its
/// formal charges.
///
/// Success replaces resonance; failure leaves all state intact.
///
/// # Errors
///
/// - [`ResonancePerceptionError::MissingAromaticity`] when aromaticity is not
///   installed.
/// - [`ResonancePerceptionError::UnknownHydrogens`] for the first atom the
///   search reads whose implicit hydrogen count is unknown.
/// - [`ResonancePerceptionError::LimitExceeded`] when the path search exceeds
///   its work bound.
///
/// ```
/// use kekule::{perception::resonance::perceive_resonance, smiles};
/// use kekule::perception::resonance::ResonanceModel;
/// let mut acetate = smiles::to_molecules("CC(=O)[O-]")?.remove(0);
/// acetate.perceive()?;
/// let resonance = acetate.perception().resonance().expect("default perception");
/// assert_eq!(resonance.model(), ResonanceModel::ChargeMoves);
/// assert_eq!(resonance.groups().len(), 1);
/// assert_eq!(resonance.bonds().len(), 2);
/// perceive_resonance(&mut acetate, ResonanceModel::ChargeMoves)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn perceive_resonance(
    mol: &mut Molecule,
    model: ResonanceModel,
) -> Result<(), ResonancePerceptionError> {
    perceive_resonance_cancellable(mol, model, None)
}

/// [`perceive_resonance`], polling `cancellation` as the path search works.
pub(crate) fn perceive_resonance_cancellable(
    mol: &mut Molecule,
    model: ResonanceModel,
    cancellation: Option<&Cancellation>,
) -> Result<(), ResonancePerceptionError> {
    let groups = match model {
        ResonanceModel::ChargeMoves => delocalization::delocalized_groups(mol, cancellation)?,
    };
    mol.perception.resonance = Some(ResonancePerception::new(model, groups));
    Ok(())
}

fn conjugated_groups(mol: &Molecule) -> Result<Vec<ResonanceGroup>, ResonanceError> {
    let state = mol
        .perception()
        .conjugation()
        .ok_or(ResonanceError::MissingConjugation)?;
    let mut remaining: BTreeSet<_> = state.bonds().collect();
    let mut groups = Vec::new();
    while let Some(&root) = remaining.first() {
        let mut atoms = BTreeSet::new();
        let mut bonds = BTreeSet::new();
        let mut stack = vec![root];
        while let Some(id) = stack.pop() {
            if !remaining.remove(&id) {
                continue;
            }
            bonds.insert(id);
            let bond = mol.bond(id).expect("installed live bond");
            for atom in [bond.a(), bond.b()] {
                if atoms.insert(atom) {
                    stack.extend(
                        mol.incident_bonds(atom)
                            .expect("live atom")
                            .map(|(id, _)| id)
                            .filter(|id| remaining.contains(id)),
                    );
                }
            }
        }
        groups.push(ResonanceGroup {
            atoms: atoms.into_iter().collect(),
            bonds: bonds.into_iter().collect(),
        });
    }
    Ok(groups)
}

/// One complete indexed contributor. Atom and bond identity is retained.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResonanceContributor {
    charges: BTreeMap<AtomId, i8>,
    orders: BTreeMap<BondId, BondOrder>,
}
impl ResonanceContributor {
    /// The formal charge of every atom of the source, in ascending atom ID
    /// order.
    pub fn formal_charges(&self) -> impl ExactSizeIterator<Item = (AtomId, i8)> + '_ {
        self.charges.iter().map(|(&id, &q)| (id, q))
    }
    /// The localized order of every bond of the source, in ascending bond ID
    /// order.
    pub fn bond_orders(&self) -> impl ExactSizeIterator<Item = (BondId, BondOrder)> + '_ {
        self.orders.iter().map(|(&id, &o)| (id, o))
    }
}

/// Enumerated contributors bound by borrowing the exact source molecule.
/// The source's represented chemistry, perception, properties and stereo are untouched.
#[derive(Debug)]
pub struct ResonanceStructures<'a> {
    source: &'a Molecule,
    groups: Vec<ResonanceGroup>,
    contributors: Vec<ResonanceContributor>,
    options: ResonanceOptions,
    limit_reached: bool,
}
impl ResonanceStructures<'_> {
    /// The connected conjugated groups the contributors vary over, which
    /// differ from installed resonance delocalization.
    pub fn groups(&self) -> &[ResonanceGroup] {
        &self.groups
    }
    /// The enumerated contributors, in deterministic order.
    pub fn contributors(&self) -> &[ResonanceContributor] {
        &self.contributors
    }
    /// The options used, with [`ResonanceFlags::effective`] flags and the
    /// structure limit clamped to [`ResonanceOptions::MAX_STRUCTURES`].
    pub const fn options(&self) -> ResonanceOptions {
        self.options
    }
    /// The configured structure limit was reached; completeness is not claimed.
    pub const fn limit_reached(&self) -> bool {
        self.limit_reached
    }
    /// Materializes one contributor, fixes hydrogen counts and clears derived state.
    /// Bond edits prune invalid represented stereo through the ordinary editor.
    ///
    /// # Errors
    ///
    /// [`ResonanceError::ContributorOutOfRange`] for an index beyond
    /// [`Self::contributors`], [`ResonanceError::InvalidElectronCount`] when
    /// an implicit hydrogen count does not fit a fixed count, and
    /// [`ResonanceError::Molecule`] or [`ResonanceError::Publication`] when
    /// editing or publishing the contributor fails.
    pub fn to_molecule(&self, index: usize) -> Result<Molecule, ResonanceError> {
        let contributor =
            self.contributors
                .get(index)
                .ok_or(ResonanceError::ContributorOutOfRange(IndexOutOfRange::new(
                    Countable::ResonanceContributor,
                    index,
                    self.contributors.len(),
                )))?;
        let mut editor = self.source.to_editor();
        for (id, _) in self.source.atoms() {
            let h = self
                .source
                .implicit_hydrogens(id)
                .ok_or(ResonanceError::UnknownHydrogens(id))?;
            let mut atom = editor
                .require_atom_mut(id)
                .map_err(ResonanceError::Molecule)?;
            atom.formal_charge = contributor.charges[&id];
            atom.hydrogens = ImplicitHydrogens::Fixed(
                u8::try_from(h).map_err(|_| ResonanceError::InvalidElectronCount)?,
            );
        }
        for (&id, &order) in &contributor.orders {
            editor
                .set_bond_order(id, order)
                .map_err(ResonanceError::Molecule)?;
        }
        editor.finish().map_err(ResonanceError::Publication)
    }
}

/// Enumerates RDKit-like resonance contributors on an immutable source.
///
/// Requires installed conjugation; contributors vary over its connected
/// groups, prepared here. All five RDKit flags are supported. Graph-indexed structures preserve localized
/// bond orders, even on aromatic atoms. This operation is never run by default
/// perception. Resource exhaustion returns no partial result.
/// At a contributor cutoff, exact ties use deterministic indexed ordering;
/// the selected subset can differ from RDKit's platform-dependent tie ordering.
///
/// # Errors
///
/// - [`ResonanceError::MissingConjugation`] when conjugation is not installed.
/// - [`ResonanceError::InvalidStructureLimit`] when the structure limit is one
///   and the molecule has a conjugated group.
/// - [`ResonanceError::UnknownHydrogens`] for the first atom whose implicit
///   hydrogen count is unknown.
/// - [`ResonanceError::LimitExceeded`] when the search exceeds
///   [`ResonanceOptions::max_total_work`].
/// - [`ResonanceError::UnsupportedBond`] or
///   [`ResonanceError::InvalidElectronCount`] when a group holds a bond or
///   electron count the search cannot represent.
///
/// ```
/// use kekule::{smiles, perception::resonance::{enumerate_resonance, ResonanceOptions}};
/// let mut acetate = smiles::to_molecules("CC(=O)[O-]")?.remove(0);
/// acetate.perceive()?;
/// let forms = enumerate_resonance(&acetate, ResonanceOptions::default())?;
/// assert_eq!(forms.contributors().len(), 2);
/// let localized = forms.to_molecule(0)?;
/// assert!(localized.perception().conjugation().is_none());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn enumerate_resonance(
    mol: &Molecule,
    mut options: ResonanceOptions,
) -> Result<ResonanceStructures<'_>, ResonanceError> {
    options.flags = options.flags.effective();
    options.max_structures = options.max_structures.min(ResonanceOptions::MAX_STRUCTURES);
    let groups = conjugated_groups(mol)?;
    if options.max_structures == 1 && !groups.is_empty() {
        return Err(ResonanceError::InvalidStructureLimit);
    }
    for id in mol.atom_ids() {
        if mol.implicit_hydrogens(id).is_none() {
            return Err(ResonanceError::UnknownHydrogens(id));
        }
    }
    let contributors = search::enumerate(mol, &groups, options)?;
    let limit_reached = contributors.len() == options.max_structures;
    Ok(ResonanceStructures {
        source: mol,
        groups,
        contributors,
        options,
        limit_reached,
    })
}
