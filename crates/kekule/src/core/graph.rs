use super::{Atom, Bond, BondId, StereoElement, StereoGroup};

/// Authoritative represented chemistry for one molecule.
///
/// `Graph` owns stable local atom and bond slots, adjacency, represented
/// stereochemistry. Structural mutation is kept
/// crate-private and is published only through `MoleculeEditor::finish`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Graph {
    pub(crate) atoms: Vec<Option<Atom>>,
    pub(crate) bonds: Vec<Option<Bond>>,
    pub(crate) adjacency: Vec<Vec<BondId>>,
    pub(crate) stereo_elements: Vec<Option<StereoElement>>,
    pub(crate) stereo_groups: Vec<Option<StereoGroup>>,
}

impl Graph {
    /// Number of live atoms. Stable slot tombstones are not counted.
    pub fn atom_count(&self) -> usize {
        self.atoms.iter().flatten().count()
    }

    /// Number of live bonds. Stable slot tombstones are not counted.
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
