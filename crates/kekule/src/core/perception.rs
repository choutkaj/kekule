use super::{AtomId, BondId};

/// The valence model used to produce installed valence state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValenceModel {
    RdkitLike,
}

/// The aromaticity model used to produce installed aromaticity state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    RdkitLike,
}

/// Cycle membership over the stable atom and bond slots of a molecule.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingMembership {
    pub(crate) atom_flags: Vec<bool>,
    pub(crate) bond_flags: Vec<bool>,
}

impl RingMembership {
    /// Constructs detached ring membership over complete stable atom and bond slots.
    ///
    /// Slot lengths and live references are checked when the containing
    /// [`super::PerceptionState`] is installed on a molecule.
    pub fn from_slot_flags(atom_flags: Vec<bool>, bond_flags: Vec<bool>) -> Self {
        Self {
            atom_flags,
            bond_flags,
        }
    }

    /// Returns the complete stable atom-slot flags, including tombstones.
    pub fn atom_slot_flags(&self) -> &[bool] {
        &self.atom_flags
    }

    /// Returns the complete stable bond-slot flags, including tombstones.
    pub fn bond_slot_flags(&self) -> &[bool] {
        &self.bond_flags
    }

    pub fn atom_in_ring(&self, atom: AtomId) -> bool {
        self.atom_flags.get(atom.index()).copied().unwrap_or(false)
    }

    pub fn bond_in_ring(&self, bond: BondId) -> bool {
        self.bond_flags.get(bond.index()).copied().unwrap_or(false)
    }

    pub fn ring_atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        (0..=u32::MAX)
            .zip(self.atom_flags.iter())
            .filter_map(|(raw, in_ring)| in_ring.then_some(AtomId::new(raw)))
    }

    pub fn ring_bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        (0..=u32::MAX)
            .zip(self.bond_flags.iter())
            .filter_map(|(raw, in_ring)| in_ring.then_some(BondId::new(raw)))
    }
}

/// One ring in an installed deterministic ring basis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
}

/// A deterministic ring basis installed in molecule perception state.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingSet {
    pub(crate) rings: Vec<Ring>,
}

impl RingSet {
    /// Constructs a detached deterministic ring basis.
    ///
    /// Ring references and graph coherence are checked when the containing
    /// [`super::PerceptionState`] is installed on a molecule.
    pub fn from_rings(rings: Vec<Ring>) -> Self {
        Self { rings }
    }

    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    pub fn len(&self) -> usize {
        self.rings.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rings.is_empty()
    }
}
