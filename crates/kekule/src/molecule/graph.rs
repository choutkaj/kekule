use super::owner::MoleculeError;
use super::{
    Atom, AtomId, Bond, BondId, StereoElement, StereoElementId, StereoGroup, StereoGroupId,
};

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

    /// Returns the sum of the asserted formal charges on all live atoms.
    ///
    /// This aggregate does not require perception.
    pub fn formal_charge(&self) -> i64 {
        self.atoms()
            .map(|(_, atom)| i64::from(atom.formal_charge))
            .sum()
    }

    /// One live atom, or `None` when `id` names no atom of this graph.
    pub fn atom(&self, id: AtomId) -> Option<&Atom> {
        self.atoms.get(id.index()).and_then(Option::as_ref)
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        (0..=u32::MAX)
            .zip(self.atoms.iter())
            .filter_map(|(raw, atom)| atom.as_ref().map(|atom| (AtomId::new(raw), atom)))
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms().map(|(id, _)| id)
    }

    /// One live bond, or `None` when `id` names no bond of this graph.
    pub fn bond(&self, id: BondId) -> Option<&Bond> {
        self.bonds.get(id.index()).and_then(Option::as_ref)
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)> {
        (0..=u32::MAX)
            .zip(self.bonds.iter())
            .filter_map(|(raw, bond)| bond.as_ref().map(|bond| (BondId::new(raw), bond)))
    }

    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bonds().map(|(id, _)| id)
    }

    /// Neighbours of one atom in adjacency order, or `None` when `id` names no
    /// atom of this graph.
    pub fn neighbors(&self, id: AtomId) -> Option<impl Iterator<Item = AtomId> + '_> {
        Some(
            self.incident_bonds(id)?
                .map(move |(_, bond)| bond.other_atom(id)),
        )
    }

    /// Returns graph components for validation and graph algorithms.
    ///
    /// A completed nonempty public molecule has exactly one component. The
    /// general result shape also supports empty values and private builder,
    /// editor, and format-interpretation staging.
    pub(crate) fn connected_components(&self) -> Vec<Vec<AtomId>> {
        let mut seen = vec![false; self.atoms.len()];
        let mut components = Vec::new();
        for start in self.atom_ids() {
            if seen[start.index()] {
                continue;
            }
            seen[start.index()] = true;
            let mut stack = vec![start];
            let mut component = Vec::new();
            while let Some(atom) = stack.pop() {
                component.push(atom);
                let mut neighbors = self
                    .neighbors(atom)
                    .expect("component traversal visits live atoms")
                    .filter(|neighbor| !seen[neighbor.index()])
                    .collect::<Vec<_>>();
                neighbors.sort_unstable_by(|left, right| right.cmp(left));
                for neighbor in neighbors {
                    if !seen[neighbor.index()] {
                        seen[neighbor.index()] = true;
                        stack.push(neighbor);
                    }
                }
            }
            component.sort_unstable();
            components.push(component);
        }
        components
    }

    /// Bonds incident to one atom in adjacency order, or `None` when `id`
    /// names no atom of this graph.
    pub fn incident_bonds(&self, id: AtomId) -> Option<impl Iterator<Item = (BondId, &Bond)> + '_> {
        self.atom(id)?;
        Some(self.adjacency[id.index()].iter().map(|bond_id| {
            let bond = self
                .bond(*bond_id)
                .expect("adjacency references only live bonds");
            (*bond_id, bond)
        }))
    }

    /// The bond joining two atoms, or `None` when they are not bonded or
    /// either names no atom of this graph.
    pub fn bond_between(&self, a: AtomId, b: AtomId) -> Option<BondId> {
        self.atom(b)?;
        self.incident_bonds(a)?
            .find(|(_, bond)| bond.connects(a, b))
            .map(|(id, _)| id)
    }

    /// One stereo element, or `None` when `id` names no element of this graph.
    pub fn stereo_element(&self, id: StereoElementId) -> Option<&StereoElement> {
        self.stereo_elements
            .get(id.index())
            .and_then(Option::as_ref)
    }

    pub fn stereo_elements(&self) -> impl Iterator<Item = (StereoElementId, &StereoElement)> {
        (0..=u32::MAX)
            .zip(self.stereo_elements.iter())
            .filter_map(|(raw, element)| {
                element
                    .as_ref()
                    .map(|element| (StereoElementId::new(raw), element))
            })
    }

    pub fn stereo_element_ids(&self) -> impl Iterator<Item = StereoElementId> + '_ {
        self.stereo_elements().map(|(id, _)| id)
    }

    /// One stereo group, or `None` when `id` names no group of this graph.
    pub fn stereo_group(&self, id: StereoGroupId) -> Option<&StereoGroup> {
        self.stereo_groups.get(id.index()).and_then(Option::as_ref)
    }

    /// [`Self::atom`] for staging code that reports a missing atom.
    pub(crate) fn require_atom(&self, id: AtomId) -> Result<&Atom, MoleculeError> {
        self.atom(id).ok_or(MoleculeError::InvalidAtomId(id))
    }

    /// [`Self::bond`] for staging code that reports a missing bond.
    pub(crate) fn require_bond(&self, id: BondId) -> Result<&Bond, MoleculeError> {
        self.bond(id).ok_or(MoleculeError::InvalidBondId(id))
    }

    /// [`Self::neighbors`] for staging code that reports a missing atom.
    pub(crate) fn require_neighbors(
        &self,
        id: AtomId,
    ) -> Result<impl Iterator<Item = AtomId> + '_, MoleculeError> {
        self.neighbors(id).ok_or(MoleculeError::InvalidAtomId(id))
    }

    /// [`Self::incident_bonds`] for staging code that reports a missing atom.
    pub(crate) fn require_incident_bonds(
        &self,
        id: AtomId,
    ) -> Result<impl Iterator<Item = (BondId, &Bond)> + '_, MoleculeError> {
        self.incident_bonds(id)
            .ok_or(MoleculeError::InvalidAtomId(id))
    }

    /// [`Self::bond_between`] for staging code that reports a missing atom.
    pub(crate) fn require_bond_between(
        &self,
        a: AtomId,
        b: AtomId,
    ) -> Result<Option<BondId>, MoleculeError> {
        self.require_atom(a)?;
        self.require_atom(b)?;
        Ok(self.bond_between(a, b))
    }

    /// [`Self::stereo_element`] for staging code that reports a missing element.
    pub(crate) fn require_stereo_element(
        &self,
        id: StereoElementId,
    ) -> Result<&StereoElement, MoleculeError> {
        self.stereo_element(id)
            .ok_or(MoleculeError::InvalidStereoElementId(id))
    }

    /// [`Self::stereo_group`] for staging code that reports a missing group.
    pub(crate) fn require_stereo_group(
        &self,
        id: StereoGroupId,
    ) -> Result<&StereoGroup, MoleculeError> {
        self.stereo_group(id)
            .ok_or(MoleculeError::InvalidStereoGroupId(id))
    }

    pub fn stereo_groups(&self) -> impl Iterator<Item = (StereoGroupId, &StereoGroup)> {
        (0..=u32::MAX)
            .zip(self.stereo_groups.iter())
            .filter_map(|(raw, group)| group.as_ref().map(|group| (StereoGroupId::new(raw), group)))
    }
}
