use std::sync::Arc;

use crate::core::Bond;
use crate::topology::{
    InstanceAtomId, InstanceBondId, MoleculeDefinitionId, MoleculeInstanceId, Topology,
    TopologyBondIndex,
};

use super::{
    combine_indices, insert_index, remove_index, toggle_index, AtomSelection, SelectionError,
};

/// Endpoint rule for converting an atom selection into a bond selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondSelectionMode {
    /// Both endpoints are selected (the induced bond set).
    Internal,
    /// At least one endpoint is selected, including internal bonds.
    Incident,
    /// Exactly one endpoint is selected (the cut boundary).
    Boundary,
}

/// A topology-bound, sorted, unique dense bond selection.
///
/// Bonds can be selected independently of atoms. [`Self::to_atoms`] explicitly
/// selects their endpoints; [`AtomSelection::to_bonds`] explicitly selects bonds
/// by an endpoint rule. Neither operation edits chemistry.
///
/// Like [`AtomSelection`], equality and set operations require the exact shared
/// topology snapshot, including for empty sets. Bare IDs and dense indices are
/// interpreted in that snapshot and do not carry provenance of their own.
///
/// # Editing selected bonds
///
/// Resolve selected source IDs to draft handles before deleting them. Publication
/// can split or merge molecule instances and returns a new topology; the original
/// selection continues to refer to its original snapshot.
///
/// ```
/// use std::sync::Arc;
/// use kekule::{smiles, topology::BondSelection};
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let topology = Arc::new(smiles::to_topology("CCC")?);
/// let selected = BondSelection::from_bonds(&topology, [topology.bond_ids()[0]])?;
/// selected.ensure_compatible(&topology)?;
/// let mut editor = topology.edit();
/// let handles = selected.bond_ids()
///     .map(|id| editor.bond_handle(id))
///     .collect::<Result<Vec<_>, _>>()?;
/// editor.delete_bonds(handles)?;
/// let edited = editor.finish()?;
/// assert_eq!(edited.instance_count(), 2);
/// assert_eq!(edited.bond_count(), 1);
/// assert_eq!(topology.bond_count(), 2);
/// assert!(selected.ensure_compatible(&edited).is_err());
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct BondSelection {
    topology: Arc<Topology>,
    indices: Vec<TopologyBondIndex>,
}

impl PartialEq for BondSelection {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.indices == other.indices
    }
}

impl Eq for BondSelection {}

impl BondSelection {
    /// An empty selection bound to this exact snapshot.
    pub fn empty(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            indices: Vec::new(),
        }
    }

    /// Selects every bond in authoritative dense order. A topology containing
    /// only isolated atoms has an empty bond selection.
    pub fn all(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            indices: (0..topology.bond_count())
                .map(|index| TopologyBondIndex::new(index as u32))
                .collect(),
        }
    }

    /// Validates IDs, removes duplicates, and orders bonds by dense index.
    pub fn from_bonds(
        topology: &Arc<Topology>,
        bonds: impl IntoIterator<Item = InstanceBondId>,
    ) -> Result<Self, SelectionError> {
        let indices = bonds
            .into_iter()
            .map(|id| {
                topology
                    .bond_index(id)
                    .ok_or(SelectionError::InvalidBondId(id))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_indices(topology, indices)
    }

    /// Validates indices, removes duplicates, and orders bonds by dense index.
    pub fn from_indices(
        topology: &Arc<Topology>,
        indices: impl IntoIterator<Item = TopologyBondIndex>,
    ) -> Result<Self, SelectionError> {
        let mut indices = indices
            .into_iter()
            .map(|index| {
                topology
                    .bond_id(index)
                    .ok_or(SelectionError::InvalidBondIndex(index))?;
                Ok(index)
            })
            .collect::<Result<Vec<_>, SelectionError>>()?;
        indices.sort_unstable();
        indices.dedup();
        Ok(Self {
            topology: Arc::clone(topology),
            indices,
        })
    }

    /// Selects bonds in the requested complete molecule instances.
    pub fn for_instances(
        topology: &Arc<Topology>,
        instances: impl IntoIterator<Item = MoleculeInstanceId>,
    ) -> Result<Self, SelectionError> {
        Ok(
            AtomSelection::for_instances(topology, instances)?
                .to_bonds(BondSelectionMode::Internal),
        )
    }

    /// Selects bonds in all instances of the requested reusable definitions.
    pub fn for_definitions(
        topology: &Arc<Topology>,
        definitions: impl IntoIterator<Item = MoleculeDefinitionId>,
    ) -> Result<Self, SelectionError> {
        Ok(AtomSelection::for_definitions(topology, definitions)?
            .to_bonds(BondSelectionMode::Internal))
    }

    /// Selects bonds satisfying a predicate, in dense order. The ID permits
    /// property or perception lookup; perception is never computed implicitly.
    pub fn from_predicate(
        topology: &Arc<Topology>,
        mut predicate: impl FnMut(InstanceBondId, &Bond) -> bool,
    ) -> Self {
        Self::from_bonds(
            topology,
            topology
                .bonds()
                .filter(|(id, bond)| predicate(*id, bond))
                .map(|(id, _)| id),
        )
        .expect("predicate selects validated topology bonds")
    }

    /// Keeps selected bonds satisfying a predicate, preserving dense order.
    pub fn filter(&self, mut predicate: impl FnMut(InstanceBondId, &Bond) -> bool) -> Self {
        Self::from_bonds(
            &self.topology,
            self.bond_ids()
                .filter(|id| predicate(*id, self.topology.bond(*id).expect("validated bond"))),
        )
        .expect("filter selects validated topology bonds")
    }

    pub fn len(&self) -> usize {
        self.indices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    /// Tests membership in this snapshot. Invalid IDs return false.
    pub fn contains(&self, bond: InstanceBondId) -> bool {
        self.topology
            .bond_index(bond)
            .is_some_and(|index| self.contains_index(index))
    }

    /// Tests dense-index membership. Out-of-range indices return false.
    pub fn contains_index(&self, index: TopologyBondIndex) -> bool {
        self.indices.binary_search(&index).is_ok()
    }

    /// Selects a bond; returns whether membership changed. Invalid IDs leave
    /// the selection unchanged.
    pub fn insert(&mut self, bond: InstanceBondId) -> Result<bool, SelectionError> {
        let index = self
            .topology
            .bond_index(bond)
            .ok_or(SelectionError::InvalidBondId(bond))?;
        Ok(insert_index(&mut self.indices, index))
    }

    /// Deselects a bond; returns whether membership changed. Invalid IDs leave
    /// the selection unchanged.
    pub fn remove(&mut self, bond: InstanceBondId) -> Result<bool, SelectionError> {
        let index = self
            .topology
            .bond_index(bond)
            .ok_or(SelectionError::InvalidBondId(bond))?;
        Ok(remove_index(&mut self.indices, index))
    }

    /// Toggles a bond and returns its new selected state. Invalid IDs leave
    /// the selection unchanged.
    pub fn toggle(&mut self, bond: InstanceBondId) -> Result<bool, SelectionError> {
        let index = self
            .topology
            .bond_index(bond)
            .ok_or(SelectionError::InvalidBondId(bond))?;
        Ok(toggle_index(&mut self.indices, index))
    }

    /// Removes all members while retaining the topology binding.
    pub fn clear(&mut self) {
        self.indices.clear();
    }

    /// Combines selections from the exact same topology snapshot.
    pub fn union(&self, other: &Self) -> Result<Self, SelectionError> {
        self.combine(other, true, true, true)
    }

    /// Keeps bonds present in both selections.
    pub fn intersection(&self, other: &Self) -> Result<Self, SelectionError> {
        self.combine(other, false, true, false)
    }

    /// Keeps selected bonds absent from `other`.
    pub fn difference(&self, other: &Self) -> Result<Self, SelectionError> {
        self.combine(other, true, false, false)
    }

    /// Keeps bonds present in exactly one selection, for group toggling.
    pub fn symmetric_difference(&self, other: &Self) -> Result<Self, SelectionError> {
        self.combine(other, true, false, true)
    }

    /// Tests whether every selected bond is selected in `other`.
    /// Snapshot compatibility is checked even for empty selections.
    pub fn is_subset(&self, other: &Self) -> Result<bool, SelectionError> {
        self.ensure_compatible(&other.topology)?;
        Ok(self
            .indices
            .iter()
            .all(|index| other.contains_index(*index)))
    }

    /// Tests whether the selections share no bonds, after checking snapshots.
    pub fn is_disjoint(&self, other: &Self) -> Result<bool, SelectionError> {
        self.ensure_compatible(&other.topology)?;
        Ok(self
            .indices
            .iter()
            .all(|index| !other.contains_index(*index)))
    }

    /// Inverts membership within this exact topology.
    pub fn complement(&self) -> Self {
        Self::all(&self.topology)
            .difference(self)
            .expect("same snapshot")
    }

    fn combine(
        &self,
        other: &Self,
        left_only: bool,
        both: bool,
        right_only: bool,
    ) -> Result<Self, SelectionError> {
        self.ensure_compatible(&other.topology)?;
        Ok(Self {
            topology: Arc::clone(&self.topology),
            indices: combine_indices(&self.indices, &other.indices, left_only, both, right_only),
        })
    }

    /// Selects both endpoints of every selected bond, deduplicated in atom
    /// topology order. Converting back with `Internal` can add other bonds
    /// between those atoms, so this is not generally a lossless round trip.
    pub fn to_atoms(&self) -> AtomSelection {
        AtomSelection::from_atoms(
            &self.topology,
            self.bond_ids().flat_map(|id| {
                let bond = self.topology.bond(id).expect("validated bond");
                [
                    InstanceAtomId::new(id.molecule(), bond.a()),
                    InstanceAtomId::new(id.molecule(), bond.b()),
                ]
            }),
        )
        .expect("validated endpoints")
    }

    /// Iterates semantic bond IDs in topology order without allocating.
    pub fn bond_ids(&self) -> impl ExactSizeIterator<Item = InstanceBondId> + '_ {
        self.indices.iter().map(|index| {
            self.topology
                .bond_id(*index)
                .expect("validated selection index")
        })
    }

    pub fn indices(&self) -> &[TopologyBondIndex] {
        &self.indices
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// Clones the handle to the exact snapshot without copying topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Rejects even layout-equal snapshots unless the shared owner is identical.
    pub fn ensure_compatible(&self, topology: &Arc<Topology>) -> Result<(), SelectionError> {
        if !Arc::ptr_eq(&self.topology, topology) {
            return Err(SelectionError::TopologyMismatch);
        }
        Ok(())
    }
}
