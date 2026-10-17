//! RDKit-like connected conjugated groups and explicit resonance enumeration.
//!
//! The search follows RDKit 2026.03.3 (Paolo Tosco, 2015), distributed under
//! the BSD 3-Clause license reproduced in LICENSE-RDKit in this crate.
use crate::diagnostics::LimitExceeded;
use crate::molecule::{
    AtomId, BondId, BondOrder, ImplicitHydrogens, Molecule, MoleculeError,
    MoleculePublicationError, ResonanceGroup, ResonancePerception,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

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

/// Why resonance perception or enumeration failed.
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
    ContributorOutOfRange {
        /// The requested index.
        index: usize,
        /// The number of contributors.
        count: usize,
    },
    /// Editing the contributor molecule failed; see the source.
    Molecule(MoleculeError),
    /// Publishing the contributor molecule failed; see the source.
    Publication(MoleculePublicationError),
}
impl fmt::Display for ResonanceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingConjugation => f.write_str("resonance requires installed conjugation"),
            Self::UnknownHydrogens(a) => write!(f, "resonance requires known hydrogens at {a}"),
            Self::UnsupportedBond(b) => write!(f, "resonance group contains unsupported {b}"),
            Self::InvalidElectronCount => f.write_str("invalid electron count in resonance search"),
            Self::InvalidStructureLimit => f.write_str("RDKit-compatible conjugated enumeration requires a structure limit of zero or at least two"),
            Self::LimitExceeded(limit) => write!(f, "resonance {limit}"),
            Self::ContributorOutOfRange { index, count } => write!(
                f,
                "resonance contributor {index} is out of range for {count} contributors"
            ),
            Self::Molecule(_) => f.write_str("cannot edit the resonance contributor"),
            Self::Publication(_) => f.write_str("cannot publish the resonance contributor"),
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

/// Prepares connected conjugated groups, without enumerating contributors.
/// Success installs only the group partition; failure leaves all state intact.
///
/// # Errors
///
/// [`ResonanceError::MissingConjugation`] when conjugation perception is not
/// installed.
pub fn perceive_resonance(mol: &mut Molecule) -> Result<(), ResonanceError> {
    let groups = conjugated_groups(mol)?;
    mol.perception.resonance = Some(ResonancePerception { groups });
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
    /// The connected conjugated groups the contributors vary over.
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
                .ok_or(ResonanceError::ContributorOutOfRange {
                    index,
                    count: self.contributors.len(),
                })?;
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
/// Requires installed conjugation; groups are reused or prepared locally. All
/// five RDKit flags are supported. Graph-indexed structures preserve localized
/// bond orders, even on aromatic atoms. This operation is never run by default
/// perception. Resource exhaustion returns no partial result.
/// At a contributor cutoff, exact ties use deterministic indexed ordering;
/// the selected subset can differ from RDKit's platform-dependent tie ordering.
///
/// # Errors
///
/// - [`ResonanceError::MissingConjugation`] when neither resonance groups nor
///   conjugation are installed.
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
    let groups = match mol.perception().resonance() {
        Some(state) => state.groups().to_vec(),
        None => conjugated_groups(mol)?,
    };
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
