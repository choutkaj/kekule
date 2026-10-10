//! RDKit-like connected conjugated groups and explicit resonance enumeration.
//!
//! The search follows RDKit 2026.03.3 (Paolo Tosco, 2015), distributed under
//! the BSD 3-Clause license reproduced in LICENSE-RDKit in this crate.
use crate::core::{
    AtomId, BondId, BondOrder, ImplicitHydrogens, Molecule, ResonanceGroup, ResonancePerception,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

mod kekule;
mod search;

/// RDKit resonance contributor policies. Combine constants with `|`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ResonanceFlags(u8);
impl ResonanceFlags {
    pub const ALLOW_INCOMPLETE_OCTETS: Self = Self(1);
    pub const ALLOW_CHARGE_SEPARATION: Self = Self(2);
    pub const KEKULE_ALL: Self = Self(4);
    pub const UNCONSTRAINED_CATIONS: Self = Self(8);
    pub const UNCONSTRAINED_ANIONS: Self = Self(16);
    /// Constructs flags, rejecting unknown bits.
    pub const fn from_bits(bits: u8) -> Option<Self> {
        if bits < 32 {
            Some(Self(bits))
        } else {
            None
        }
    }
    pub const fn bits(self) -> u8 {
        self.0
    }
    pub const fn contains(self, flag: Self) -> bool {
        self.0 & flag.0 == flag.0
    }
    /// Applies RDKit's implications for unconstrained ions.
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
    pub flags: ResonanceFlags,
    /// Maximum returned contributors, clamped to RDKit's ceiling of 1,000,000.
    /// Zero returns an empty collection. The pinned reference cannot produce a
    /// valid conjugated contributor at a limit of one; use at least two.
    pub max_structures: usize,
    /// Bounds search visits, copied state, and permutation construction/sorting.
    /// Exhaustion returns an error without a partial result.
    pub max_total_work: usize,
}
impl Default for ResonanceOptions {
    fn default() -> Self {
        Self {
            flags: ResonanceFlags::default(),
            max_structures: 1000,
            max_total_work: 100_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ResonanceError {
    MissingConjugation,
    UnknownHydrogens(AtomId),
    UnsupportedBond(BondId),
    InvalidElectronCount,
    InvalidStructureLimit,
    ResourceLimit { limit: usize },
    Materialization(String),
}
impl fmt::Display for ResonanceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingConjugation => f.write_str("resonance requires installed conjugation"),
            Self::UnknownHydrogens(a) => write!(f, "resonance requires known hydrogens at {a}"),
            Self::UnsupportedBond(b) => write!(f, "resonance group contains unsupported bond {b}"),
            Self::InvalidElectronCount => f.write_str("invalid electron count in resonance search"),
            Self::InvalidStructureLimit => f.write_str("RDKit-compatible conjugated enumeration requires a structure limit of zero or at least two"),
            Self::ResourceLimit { limit } => write!(f, "resonance work limit exceeded ({limit})"),
            Self::Materialization(reason) => {
                write!(f, "cannot materialize resonance contributor: {reason}")
            }
        }
    }
}
impl std::error::Error for ResonanceError {}

/// Prepares connected conjugated groups, without enumerating contributors.
/// Success installs only the group partition; failure leaves all state intact.
pub fn perceive_resonance(mol: &mut Molecule) -> Result<(), ResonanceError> {
    let groups = conjugated_groups(mol)?;
    mol.perception.resonance = Some(ResonancePerception { groups });
    Ok(())
}

fn conjugated_groups(mol: &Molecule) -> Result<Vec<ResonanceGroup>, ResonanceError> {
    let state = mol
        .perception()
        .conjugation_state()
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
    pub fn formal_charges(&self) -> impl ExactSizeIterator<Item = (AtomId, i8)> + '_ {
        self.charges.iter().map(|(&id, &q)| (id, q))
    }
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
    pub fn groups(&self) -> &[ResonanceGroup] {
        &self.groups
    }
    pub fn contributors(&self) -> &[ResonanceContributor] {
        &self.contributors
    }
    pub const fn options(&self) -> ResonanceOptions {
        self.options
    }
    /// The configured structure limit was reached; completeness is not claimed.
    pub const fn limit_reached(&self) -> bool {
        self.limit_reached
    }
    /// Materializes one contributor, fixes hydrogen counts and clears derived state.
    /// Bond edits prune invalid represented stereo through the ordinary editor.
    pub fn to_molecule(&self, index: usize) -> Result<Molecule, ResonanceError> {
        let contributor = self.contributors.get(index).ok_or_else(|| {
            ResonanceError::Materialization("contributor index out of bounds".into())
        })?;
        let mut editor = self.source.edit();
        for (id, _) in self.source.atoms() {
            let h = self
                .source
                .implicit_hydrogens(id)
                .expect("live atom")
                .ok_or(ResonanceError::UnknownHydrogens(id))?;
            let mut atom = editor
                .atom_mut(id)
                .map_err(|e| ResonanceError::Materialization(e.to_string()))?;
            atom.formal_charge = contributor.charges[&id];
            atom.hydrogens = ImplicitHydrogens::Fixed(
                u8::try_from(h).map_err(|_| ResonanceError::InvalidElectronCount)?,
            );
        }
        for (&id, &order) in &contributor.orders {
            editor
                .set_bond_order(id, order)
                .map_err(|e| ResonanceError::Materialization(e.to_string()))?;
        }
        editor
            .finish()
            .map_err(|e| ResonanceError::Materialization(e.to_string()))
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
/// ```
/// use kekule::{smiles, perception::resonance::{enumerate_resonance, ResonanceOptions}};
/// let mut acetate = smiles::to_molecules("CC(=O)[O-]")?.remove(0);
/// acetate.perceive()?;
/// let forms = enumerate_resonance(&acetate, ResonanceOptions::default())?;
/// assert_eq!(forms.contributors().len(), 2);
/// let localized = forms.to_molecule(0)?;
/// assert!(!localized.perception().has_conjugation());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn enumerate_resonance(
    mol: &Molecule,
    mut options: ResonanceOptions,
) -> Result<ResonanceStructures<'_>, ResonanceError> {
    options.flags = options.flags.effective();
    options.max_structures = options.max_structures.min(1_000_000);
    let groups = match mol.perception().resonance_state() {
        Some(state) => state.groups().to_vec(),
        None => conjugated_groups(mol)?,
    };
    if options.max_structures == 1 && !groups.is_empty() {
        return Err(ResonanceError::InvalidStructureLimit);
    }
    for id in mol.atom_ids() {
        if mol.implicit_hydrogens(id).expect("live atom").is_none() {
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
