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
    pub fn atom_in_ring(&self, atom: AtomId) -> bool {
        self.atom_flags.get(atom.index()).copied().unwrap_or(false)
    }

    pub fn bond_in_ring(&self, bond: BondId) -> bool {
        self.bond_flags.get(bond.index()).copied().unwrap_or(false)
    }

    pub fn ring_atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atom_flags
            .iter()
            .enumerate()
            .filter_map(|(index, in_ring)| in_ring.then_some(AtomId::new(index as u32)))
    }

    pub fn ring_bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bond_flags
            .iter()
            .enumerate()
            .filter_map(|(index, in_ring)| in_ring.then_some(BondId::new(index as u32)))
    }
}

/// One ring in an installed deterministic ring basis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
}

/// Resource accounting for an installed ring-basis calculation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RingWork {
    pub atom_count: usize,
    pub bond_count: usize,
    pub candidate_cycles: usize,
    pub equivalent_shortest_paths: usize,
    pub path_expansions: usize,
    pub queue_peak: usize,
    pub stack_peak: usize,
    pub total_work: usize,
}

/// A deterministic ring basis installed in molecule perception state.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingSet {
    pub(crate) rings: Vec<Ring>,
    pub(crate) work: RingWork,
}

impl RingSet {
    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    pub fn len(&self) -> usize {
        self.rings.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rings.is_empty()
    }

    pub fn work(&self) -> RingWork {
        self.work
    }
}
