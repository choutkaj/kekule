//! Graph storage for published molecules ([`Graph`]) and editor drafts
//! ([`DraftGraph`]).
use std::fmt;
use std::iter::{DoubleEndedIterator, Enumerate};

use super::{
    Atom, AtomId, Bond, BondId, MoleculeError, StereoElement, StereoElementId, StereoGroup,
    StereoGroupId,
};
use sealed::Value;

mod sealed {
    pub trait Sealed {}

    /// What one graph slot stores.
    pub trait Value: Clone + PartialEq + std::fmt::Debug {}

    impl<T: Clone + PartialEq + std::fmt::Debug> Value for T {}
}

/// Whether the identifier spaces of a [`Graph`] may hold deleted entities.
///
/// A [`Dense`] graph belongs to a published [`Molecule`](super::Molecule) and
/// numbers every identifier space `0..n`. A [`Sparse`] graph, the
/// [`DraftGraph`] of a [`MoleculeEditor`](super::MoleculeEditor), keeps the
/// slots of deleted atoms, bonds, and stereo so draft IDs stay stable until
/// publication renumbers them. Reads skip deleted slots under either policy.
/// The trait is sealed.
pub trait SlotLiveness: sealed::Sealed + 'static {
    #[doc(hidden)]
    type Slot<T: Value>: Value;
    #[doc(hidden)]
    fn value<T: Value>(slot: &Self::Slot<T>) -> Option<&T>;
}

/// Every slot is live: the graph of a published [`Molecule`](super::Molecule).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dense {}

/// Deleted slots remain until publication: the graph of a
/// [`MoleculeEditor`](super::MoleculeEditor) draft.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sparse {}

impl sealed::Sealed for Dense {}
impl sealed::Sealed for Sparse {}

impl SlotLiveness for Dense {
    type Slot<T: Value> = T;

    fn value<T: Value>(slot: &T) -> Option<&T> {
        Some(slot)
    }
}

impl SlotLiveness for Sparse {
    type Slot<T: Value> = Option<T>;

    fn value<T: Value>(slot: &Option<T>) -> Option<&T> {
        slot.as_ref()
    }
}

/// One identifier space: its slots and how many of them are live.
pub(crate) struct Slots<T: Value, L: SlotLiveness> {
    slots: Vec<L::Slot<T>>,
    live: usize,
}

impl<T: Value, L: SlotLiveness> Clone for Slots<T, L> {
    fn clone(&self) -> Self {
        Self {
            slots: self.slots.clone(),
            live: self.live,
        }
    }
}

impl<T: Value, L: SlotLiveness> fmt::Debug for Slots<T, L> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.slots.fmt(formatter)
    }
}

impl<T: Value, L: SlotLiveness> PartialEq for Slots<T, L> {
    fn eq(&self, other: &Self) -> bool {
        self.slots == other.slots
    }
}

impl<T: Value, L: SlotLiveness> Default for Slots<T, L> {
    fn default() -> Self {
        Self {
            slots: Vec::new(),
            live: 0,
        }
    }
}

impl<T: Value, L: SlotLiveness> Slots<T, L> {
    pub(crate) fn len(&self) -> usize {
        self.live
    }

    pub(crate) fn slot_count(&self) -> usize {
        self.slots.len()
    }

    pub(crate) fn get(&self, index: usize) -> Option<&T> {
        self.slots.get(index).and_then(L::value)
    }

    pub(crate) fn iter(&self) -> LiveSlots<'_, T, L> {
        LiveSlots {
            slots: self.slots.iter().enumerate(),
            remaining: self.live,
        }
    }

    fn to_sparse(&self) -> Slots<T, Sparse> {
        Slots {
            slots: self
                .slots
                .iter()
                .map(|slot| L::value(slot).cloned())
                .collect(),
            live: self.live,
        }
    }
}

impl<T: Value> Slots<T, Dense> {
    pub(crate) fn from_values(values: Vec<T>) -> Self {
        Self {
            live: values.len(),
            slots: values,
        }
    }

    fn get_mut(&mut self, index: usize) -> Option<&mut T> {
        self.slots.get_mut(index)
    }
}

impl<T: Value> Slots<T, Sparse> {
    /// Stores `value` in a new slot and returns its index.
    fn push(&mut self, value: T) -> usize {
        self.slots.push(Some(value));
        self.live += 1;
        self.slots.len() - 1
    }

    fn take(&mut self, index: usize) -> Option<T> {
        let value = self.slots.get_mut(index)?.take()?;
        self.live -= 1;
        Some(value)
    }

    fn get_mut(&mut self, index: usize) -> Option<&mut T> {
        self.slots.get_mut(index).and_then(Option::as_mut)
    }

    pub(crate) fn has_deleted(&self) -> bool {
        self.live != self.slots.len()
    }

    /// The dense index each slot receives at publication; deleted slots get
    /// `None`.
    fn dense_indices<Id>(&self, id: impl Fn(u32) -> Id) -> Vec<Option<Id>> {
        let mut next = 0;
        self.slots
            .iter()
            .map(|slot| {
                slot.as_ref().map(|_| {
                    next += 1;
                    id(next - 1)
                })
            })
            .collect()
    }

    /// Collects in place when `Option<T>` is no larger than `T`, so publishing a
    /// draft without deletions keeps its storage.
    fn into_dense(self, mut map: impl FnMut(T) -> T) -> Slots<T, Dense> {
        Slots::from_values(
            self.slots
                .into_iter()
                .filter_map(|slot| slot.map(&mut map))
                .collect(),
        )
    }
}

impl<T: Value> From<Slots<T, Dense>> for Slots<T, Sparse> {
    fn from(slots: Slots<T, Dense>) -> Self {
        Self {
            live: slots.live,
            slots: slots.slots.into_iter().map(Some).collect(),
        }
    }
}

impl<T: Value> crate::properties::RowLiveness for Slots<T, Sparse> {
    fn is_live(&self, row: usize) -> bool {
        self.slots.get(row).is_none_or(Option::is_some)
    }
}

/// The live slots of one identifier space, with their indices.
pub(crate) struct LiveSlots<'a, T: Value + 'a, L: SlotLiveness> {
    slots: Enumerate<std::slice::Iter<'a, L::Slot<T>>>,
    remaining: usize,
}

impl<'a, T: Value + 'a, L: SlotLiveness> Iterator for LiveSlots<'a, T, L> {
    type Item = (usize, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        let found = self
            .slots
            .by_ref()
            .find_map(|(index, slot)| L::value(slot).map(|value| (index, value)))?;
        self.remaining -= 1;
        Some(found)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<T: Value, L: SlotLiveness> DoubleEndedIterator for LiveSlots<'_, T, L> {
    fn next_back(&mut self) -> Option<Self::Item> {
        loop {
            let (index, slot) = self.slots.next_back()?;
            if let Some(value) = L::value(slot) {
                self.remaining -= 1;
                return Some((index, value));
            }
        }
    }
}

impl<T: Value, L: SlotLiveness> ExactSizeIterator for LiveSlots<'_, T, L> {}

/// The identifier of slot `index`; insertion checked that it fits.
fn id_at<Id>(index: usize, id: impl FnOnce(u32) -> Id) -> Id {
    id(u32::try_from(index).expect("slot indices fit their identifier width"))
}

/// Authoritative represented chemistry: atoms, bonds, adjacency, and
/// represented stereochemistry.
///
/// A published [`Molecule`](super::Molecule) owns a dense `Graph`; a
/// [`MoleculeEditor`](super::MoleculeEditor) owns a [`DraftGraph`] whose
/// deleted slots publication removes by renumbering every identifier space
/// densely (see [`SlotLiveness`]). Both answer the same reads. Counts are
/// constant time and iteration is exact-size, in ascending ID order.
pub struct Graph<L: SlotLiveness = Dense> {
    pub(crate) atoms: Slots<Atom, L>,
    pub(crate) bonds: Slots<Bond, L>,
    pub(crate) adjacency: Vec<Vec<BondId>>,
    pub(crate) stereo_elements: Slots<StereoElement, L>,
    pub(crate) stereo_groups: Slots<StereoGroup, L>,
}

/// The graph of a [`MoleculeEditor`](super::MoleculeEditor) draft, which may
/// hold deleted slots and be empty or disconnected.
pub type DraftGraph = Graph<Sparse>;

impl<L: SlotLiveness> Clone for Graph<L> {
    fn clone(&self) -> Self {
        Self {
            atoms: self.atoms.clone(),
            bonds: self.bonds.clone(),
            adjacency: self.adjacency.clone(),
            stereo_elements: self.stereo_elements.clone(),
            stereo_groups: self.stereo_groups.clone(),
        }
    }
}

impl<L: SlotLiveness> fmt::Debug for Graph<L> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Graph")
            .field("atoms", &self.atoms)
            .field("bonds", &self.bonds)
            .field("adjacency", &self.adjacency)
            .field("stereo_elements", &self.stereo_elements)
            .field("stereo_groups", &self.stereo_groups)
            .finish()
    }
}

impl<L: SlotLiveness> Default for Graph<L> {
    fn default() -> Self {
        Self {
            atoms: Slots::default(),
            bonds: Slots::default(),
            adjacency: Vec::new(),
            stereo_elements: Slots::default(),
            stereo_groups: Slots::default(),
        }
    }
}

impl<L: SlotLiveness> PartialEq for Graph<L> {
    fn eq(&self, other: &Self) -> bool {
        // Adjacency is an index over represented bonds. Rewiring can change its
        // traversal order without changing any asserted atom, bond, or stereo.
        self.atoms == other.atoms
            && self.bonds == other.bonds
            && self.stereo_elements == other.stereo_elements
            && self.stereo_groups == other.stereo_groups
    }
}

impl<L: SlotLiveness> Graph<L> {
    /// Number of atoms.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Number of bonds.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Atom slots, deleted ones included; per-atom scratch arrays have this
    /// length.
    pub(crate) fn atom_slot_count(&self) -> usize {
        self.atoms.slot_count()
    }

    /// Bond slots, deleted ones included.
    pub(crate) fn bond_slot_count(&self) -> usize {
        self.bonds.slot_count()
    }

    /// The sum of the asserted formal charges; it needs no perception.
    pub fn formal_charge(&self) -> i64 {
        self.atoms()
            .map(|(_, atom)| i64::from(atom.formal_charge))
            .sum()
    }

    /// One atom, or `None` when `id` names no atom of this graph.
    pub fn atom(&self, id: AtomId) -> Option<&Atom> {
        self.atoms.get(id.index())
    }

    /// Atoms with their IDs.
    pub fn atoms(
        &self,
    ) -> impl ExactSizeIterator<Item = (AtomId, &Atom)> + DoubleEndedIterator + '_ {
        self.atoms
            .iter()
            .map(|(index, atom)| (id_at(index, AtomId::new), atom))
    }

    /// Atom IDs.
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = AtomId> + DoubleEndedIterator + '_ {
        self.atoms().map(|(id, _)| id)
    }

    /// One bond, or `None` when `id` names no bond of this graph.
    pub fn bond(&self, id: BondId) -> Option<&Bond> {
        self.bonds.get(id.index())
    }

    /// Bonds with their IDs.
    pub fn bonds(
        &self,
    ) -> impl ExactSizeIterator<Item = (BondId, &Bond)> + DoubleEndedIterator + '_ {
        self.bonds
            .iter()
            .map(|(index, bond)| (id_at(index, BondId::new), bond))
    }

    /// Bond IDs.
    pub fn bond_ids(&self) -> impl ExactSizeIterator<Item = BondId> + DoubleEndedIterator + '_ {
        self.bonds().map(|(id, _)| id)
    }

    /// Neighbours of one atom in adjacency order, or `None` when `id` names no
    /// atom of this graph.
    pub fn neighbors(&self, id: AtomId) -> Option<impl ExactSizeIterator<Item = AtomId> + '_> {
        Some(
            self.incident_bonds(id)?
                .map(move |(_, bond)| bond.other_atom(id)),
        )
    }

    /// Bonds incident to one atom in adjacency order, or `None` when `id`
    /// names no atom of this graph.
    pub fn incident_bonds(
        &self,
        id: AtomId,
    ) -> Option<impl ExactSizeIterator<Item = (BondId, &Bond)> + '_> {
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
        self.stereo_elements.get(id.index())
    }

    /// Represented stereo elements with their IDs.
    pub fn stereo_elements(
        &self,
    ) -> impl ExactSizeIterator<Item = (StereoElementId, &StereoElement)> + DoubleEndedIterator + '_
    {
        self.stereo_elements
            .iter()
            .map(|(index, element)| (id_at(index, StereoElementId::new), element))
    }

    /// Represented stereo element IDs.
    pub fn stereo_element_ids(
        &self,
    ) -> impl ExactSizeIterator<Item = StereoElementId> + DoubleEndedIterator + '_ {
        self.stereo_elements().map(|(id, _)| id)
    }

    /// One stereo group, or `None` when `id` names no group of this graph.
    pub fn stereo_group(&self, id: StereoGroupId) -> Option<&StereoGroup> {
        self.stereo_groups.get(id.index())
    }

    /// Enhanced stereo groups with their IDs.
    pub fn stereo_groups(
        &self,
    ) -> impl ExactSizeIterator<Item = (StereoGroupId, &StereoGroup)> + DoubleEndedIterator + '_
    {
        self.stereo_groups
            .iter()
            .map(|(index, group)| (id_at(index, StereoGroupId::new), group))
    }

    /// Connected components in ascending order of their smallest atom, each
    /// sorted; isolated atoms are components of their own.
    ///
    /// A published molecule has exactly one. Drafts and format interpretation
    /// staging may have none or several.
    pub(crate) fn connected_components(&self) -> Vec<Vec<AtomId>> {
        let mut seen = vec![false; self.atom_slot_count()];
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
    ) -> Result<impl ExactSizeIterator<Item = AtomId> + '_, MoleculeError> {
        self.neighbors(id).ok_or(MoleculeError::InvalidAtomId(id))
    }

    /// [`Self::incident_bonds`] for staging code that reports a missing atom.
    pub(crate) fn require_incident_bonds(
        &self,
        id: AtomId,
    ) -> Result<impl ExactSizeIterator<Item = (BondId, &Bond)> + '_, MoleculeError> {
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
}

impl<L: SlotLiveness> Graph<L> {
    /// A dense copy as a [`Molecule`](super::Molecule) without perception,
    /// properties or publication's checks, for crate algorithms that perceive
    /// a draft. It may be empty or disconnected and must not escape the crate.
    pub(crate) fn to_molecule(&self) -> (super::Molecule, SlotCompaction) {
        let draft = DraftGraph {
            atoms: self.atoms.to_sparse(),
            bonds: self.bonds.to_sparse(),
            adjacency: self.adjacency.clone(),
            stereo_elements: self.stereo_elements.to_sparse(),
            stereo_groups: self.stereo_groups.to_sparse(),
        };
        let (graph, compaction) = draft.into_dense();
        let properties =
            crate::properties::MoleculeProperties::new(graph.atom_count(), graph.bond_count());
        let molecule = super::Molecule {
            graph,
            perception: super::Perception::default(),
            properties,
        };
        (molecule, compaction)
    }
}

impl Graph {
    /// A published atom's payload, for crate-private projection copies that
    /// re-derive their own perception and never escape the crate.
    pub(crate) fn projected_atom_mut(&mut self, id: AtomId) -> Option<&mut Atom> {
        self.atoms.get_mut(id.index())
    }

    /// The same chemistry as a draft whose IDs equal the published ones.
    pub(crate) fn into_draft(self) -> DraftGraph {
        Graph {
            atoms: self.atoms.into(),
            bonds: self.bonds.into(),
            adjacency: self.adjacency,
            stereo_elements: self.stereo_elements.into(),
            stereo_groups: self.stereo_groups.into(),
        }
    }
}

/// Storage primitives for [`MoleculeEditor`](super::MoleculeEditor), which
/// checks identifier capacity and keeps stereo, perception, and properties
/// consistent around them.
impl DraftGraph {
    pub(crate) fn push_atom(&mut self, atom: Atom) -> AtomId {
        self.adjacency.push(Vec::new());
        id_at(self.atoms.push(atom), AtomId::new)
    }

    /// Removes an atom whose bonds are already deleted.
    pub(crate) fn take_atom(&mut self, id: AtomId) -> Option<Atom> {
        let atom = self.atoms.take(id.index())?;
        debug_assert!(self.adjacency[id.index()].is_empty());
        Some(atom)
    }

    pub(crate) fn atom_mut(&mut self, id: AtomId) -> Option<&mut Atom> {
        self.atoms.get_mut(id.index())
    }

    /// Stores a bond between live, distinct, unbonded atoms.
    pub(crate) fn push_bond(&mut self, bond: Bond) -> BondId {
        let (a, b) = bond.endpoints();
        let id = id_at(self.bonds.push(bond), BondId::new);
        self.adjacency[a.index()].push(id);
        self.adjacency[b.index()].push(id);
        id
    }

    pub(crate) fn take_bond(&mut self, id: BondId) -> Option<Bond> {
        let bond = self.bonds.take(id.index())?;
        for atom in [bond.a, bond.b] {
            self.adjacency[atom.index()].retain(|incident| *incident != id);
        }
        Some(bond)
    }

    /// Mutable bond state; changing endpoints here would break adjacency, so
    /// rewire through [`Self::rewire_bond`].
    pub(crate) fn bond_mut(&mut self, id: BondId) -> Option<&mut Bond> {
        self.bonds.get_mut(id.index())
    }

    /// Replaces a live bond with one between live, distinct, unbonded atoms.
    pub(crate) fn rewire_bond(&mut self, id: BondId, replacement: Bond) {
        let previous = self.bonds.get(id.index()).expect("live bond").clone();
        for atom in [previous.a, previous.b] {
            self.adjacency[atom.index()].retain(|incident| *incident != id);
        }
        self.adjacency[replacement.a.index()].push(id);
        self.adjacency[replacement.b.index()].push(id);
        *self.bonds.get_mut(id.index()).expect("live bond") = replacement;
    }

    /// The ID the next stereo element receives, which issues name.
    pub(crate) fn next_stereo_element_id(&self) -> StereoElementId {
        StereoElementId::new(u32::try_from(self.stereo_elements.slot_count()).unwrap_or(u32::MAX))
    }

    /// The ID the next stereo group receives, which issues name.
    pub(crate) fn next_stereo_group_id(&self) -> StereoGroupId {
        StereoGroupId::new(u32::try_from(self.stereo_groups.slot_count()).unwrap_or(u32::MAX))
    }

    pub(crate) fn push_stereo_element(&mut self, element: StereoElement) -> StereoElementId {
        id_at(self.stereo_elements.push(element), StereoElementId::new)
    }

    pub(crate) fn stereo_element_mut(&mut self, id: StereoElementId) -> Option<&mut StereoElement> {
        self.stereo_elements.get_mut(id.index())
    }

    /// Removes an element and drops it from every group listing it, deleting
    /// a group it leaves empty.
    pub(crate) fn take_stereo_element(&mut self, id: StereoElementId) -> Option<StereoElement> {
        let element = self.stereo_elements.take(id.index())?;
        for index in 0..self.stereo_groups.slot_count() {
            let emptied = self.stereo_groups.get_mut(index).is_some_and(|group| {
                group.members.retain(|member| *member != id);
                group.members.is_empty()
            });
            if emptied {
                self.stereo_groups.take(index);
            }
        }
        Some(element)
    }

    pub(crate) fn push_stereo_group(&mut self, group: StereoGroup) -> StereoGroupId {
        id_at(self.stereo_groups.push(group), StereoGroupId::new)
    }

    pub(crate) fn stereo_group_mut(&mut self, id: StereoGroupId) -> Option<&mut StereoGroup> {
        self.stereo_groups.get_mut(id.index())
    }

    pub(crate) fn take_stereo_group(&mut self, id: StereoGroupId) -> Option<StereoGroup> {
        self.stereo_groups.take(id.index())
    }

    /// Removes every stereo element and group, restarting their IDs.
    pub(crate) fn clear_stereo(&mut self) {
        self.stereo_elements = Slots::default();
        self.stereo_groups = Slots::default();
    }

    /// Whether any identifier space holds a deleted slot.
    pub(crate) fn has_deleted(&self) -> bool {
        self.atoms.has_deleted()
            || self.bonds.has_deleted()
            || self.stereo_elements.has_deleted()
            || self.stereo_groups.has_deleted()
    }

    /// Renumbers every identifier space densely, preserving relative order, so
    /// canonical stereo carrier order survives.
    pub(crate) fn into_dense(self) -> (Graph, SlotCompaction) {
        let compaction = SlotCompaction {
            atoms: self.atoms.dense_indices(AtomId::new),
            bonds: self.bonds.dense_indices(BondId::new),
            stereo_elements: self.stereo_elements.dense_indices(StereoElementId::new),
            stereo_groups: self.stereo_groups.dense_indices(StereoGroupId::new),
        };
        if !self.has_deleted() {
            let graph = Graph {
                atoms: self.atoms.into_dense(|atom| atom),
                bonds: self.bonds.into_dense(|bond| bond),
                adjacency: self.adjacency,
                stereo_elements: self.stereo_elements.into_dense(|element| element),
                stereo_groups: self.stereo_groups.into_dense(|group| group),
            };
            return (graph, compaction);
        }
        let atom = |id: AtomId| compaction.atoms[id.index()].expect("live atom reference");
        let bond = |id: BondId| compaction.bonds[id.index()].expect("live bond reference");
        let adjacency = self
            .adjacency
            .iter()
            .enumerate()
            .filter(|(slot, _)| self.atoms.get(*slot).is_some())
            .map(|(_, incident)| incident.iter().map(|&id| bond(id)).collect())
            .collect();
        let graph = Graph {
            atoms: self.atoms.into_dense(|atom| atom),
            bonds: self
                .bonds
                .into_dense(|value| Bond::new(atom(value.a), atom(value.b), value.order)),
            adjacency,
            stereo_elements: self.stereo_elements.into_dense(|mut element| {
                element.map_references(atom, bond);
                element.group = element.group.map(|id| {
                    compaction.stereo_groups[id.index()].expect("live stereo group reference")
                });
                element
            }),
            stereo_groups: self.stereo_groups.into_dense(|mut group| {
                for member in &mut group.members {
                    *member = compaction.stereo_elements[member.index()]
                        .expect("live stereo group member");
                }
                group
            }),
        };
        (graph, compaction)
    }
}

/// Dense renumbering of one draft's identifier spaces at publication.
///
/// Each vector is indexed by a draft slot; deleted slots map to `None`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct SlotCompaction {
    pub(crate) atoms: Vec<Option<AtomId>>,
    pub(crate) bonds: Vec<Option<BondId>>,
    pub(crate) stereo_elements: Vec<Option<StereoElementId>>,
    pub(crate) stereo_groups: Vec<Option<StereoGroupId>>,
}

impl SlotCompaction {
    /// The draft slots that survive in each space, in published order.
    pub(crate) fn live_atom_slots(&self) -> Vec<usize> {
        live_slots(&self.atoms)
    }

    pub(crate) fn live_bond_slots(&self) -> Vec<usize> {
        live_slots(&self.bonds)
    }

    /// The published ID of live draft atom `draft`.
    pub(crate) fn atom(&self, draft: AtomId) -> AtomId {
        self.atoms[draft.index()].expect("live draft atom")
    }

    /// The published ID of live draft bond `draft`.
    pub(crate) fn bond(&self, draft: BondId) -> BondId {
        self.bonds[draft.index()].expect("live draft bond")
    }

    /// The draft ID of published stereo element `published`.
    pub(crate) fn draft_stereo_element(&self, published: StereoElementId) -> StereoElementId {
        id_at(
            live_slots(&self.stereo_elements)[published.index()],
            StereoElementId::new,
        )
    }

    /// Rewrites `element`'s atom and bond references from published IDs back
    /// to draft IDs.
    pub(crate) fn restore_draft_references(&self, element: &mut StereoElement) {
        let atoms = self.live_atom_slots();
        let bonds = self.live_bond_slots();
        element.map_references(
            |id| id_at(atoms[id.index()], AtomId::new),
            |id| id_at(bonds[id.index()], BondId::new),
        );
    }
}

fn live_slots<Id>(slots: &[Option<Id>]) -> Vec<usize> {
    slots
        .iter()
        .enumerate()
        .filter_map(|(slot, id)| id.as_ref().map(|_| slot))
        .collect()
}
