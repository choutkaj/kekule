use super::{Atom, Bond, BondId, StereoElement, StereoGroup};

/// Authoritative represented chemistry for one molecule.
///
/// `Graph` owns local atom and bond slots, adjacency, and represented
/// stereochemistry. Structural mutation is kept crate-private; editing drafts
/// may hold deleted slots, which `MoleculeEditor::finish` removes by
/// renumbering every ID space densely.
#[derive(Debug, Clone, Default)]
pub struct Graph {
    pub(crate) atoms: Vec<Option<Atom>>,
    pub(crate) bonds: Vec<Option<Bond>>,
    pub(crate) adjacency: Vec<Vec<BondId>>,
    pub(crate) stereo_elements: Vec<Option<StereoElement>>,
    pub(crate) stereo_groups: Vec<Option<StereoGroup>>,
}

impl PartialEq for Graph {
    fn eq(&self, other: &Self) -> bool {
        // Adjacency is an index over represented bonds. Rewiring can change its
        // traversal order without changing any asserted atom, bond, or stereo.
        self.atoms == other.atoms
            && self.bonds == other.bonds
            && self.stereo_elements == other.stereo_elements
            && self.stereo_groups == other.stereo_groups
    }
}

impl Graph {
    /// Number of live atoms. Deleted draft slots are not counted.
    pub fn atom_count(&self) -> usize {
        self.atoms.iter().flatten().count()
    }

    /// Number of live bonds. Deleted draft slots are not counted.
    pub fn bond_count(&self) -> usize {
        self.bonds.iter().flatten().count()
    }

    pub(crate) fn atom_slot_count(&self) -> usize {
        self.atoms.len()
    }

    pub(crate) fn bond_slot_count(&self) -> usize {
        self.bonds.len()
    }
}
