use std::collections::BTreeSet;
use std::fmt;
use std::ops::Range;
use std::sync::Arc;

use crate::molecule::Element;
use crate::substructure::{QueryMatch, TopologyQueryMatch};

use super::{
    AtomSiteId, AtomSiteView, AtomView, ChainId, ChainView, InstanceAtomId, InstanceBondId,
    MoleculeClass, MoleculeDefinitionId, MoleculeInstanceId, ResidueClass, ResidueId, ResidueView,
    Topology, TopologyAtomIndex, TopologyBondIndex,
};

mod bonds;
pub use bonds::{BondSelection, BondSelectionMode};

/// A topology-bound, sorted, unique dense atom selection.
///
/// Empty selections retain their topology. Equality and set operations use
/// layout identity ([`Topology::shares_layout`]), not chemical or layout equality,
/// so selections made before and after perception combine and compare equal. IDs
/// supplied to constructors and membership methods are interpreted in the bound
/// layout: bare IDs do not carry provenance. Selections are static sets, not
/// stored queries.
/// See [`BondSelection`] for independently selected bonds.
///
/// # Selection sources
///
/// | Source | Constructors |
/// | --- | --- |
/// | Picking | [`Self::from_atoms`], [`Self::from_indices`] |
/// | Dense order, such as a host's indices or bit mask | [`Self::from_sorted_indices`], [`Self::from_mask`] |
/// | All/none | [`Self::all`], [`Self::empty`] |
/// | Molecule occurrences or reused definitions | [`Self::for_instances`], [`Self::for_definitions`] |
/// | Chemistry or custom properties | [`Self::for_elements`], [`Self::for_molecule_classes`], [`Self::from_predicate`] |
/// | Hierarchy | [`Self::for_chains`], [`Self::for_chain_label`], [`Self::for_residues`], [`Self::for_residue_classes`], [`Self::for_atom_sites`] |
/// | Exact names in separate namespaces | [`Self::for_label_atom_names`], [`Self::for_author_atom_names`] |
/// | Substructure matches | [`Self::from_query_matches`], [`Self::from_topology_query_matches`] |
/// | Cartesian distance in one structure | [`crate::structure::measure::within`] |
///
/// Hierarchy constructors return atom sets, not independently selected hierarchy
/// nodes. Implicit hydrogens have no atom IDs and cannot be selected separately.
/// Predicates can inspect charge, bond environment, properties or perceived
/// chemistry through the topology, without growing a separate query language.
/// [`Self::filter`] restricts an existing set.
///
/// # Interactive selection
///
/// ```
/// use std::sync::Arc;
/// use kekule::{smiles, topology::{AtomSelection, BondSelectionMode}};
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let topology = smiles::to_topology("CCO")?;
/// let mut selected = AtomSelection::empty(topology.clone());
/// selected.toggle(topology.atom_ids()[0])?; // e.g. a picked atom
/// let neighborhood = selected.expand_bonded(1);
/// selected = selected.union(&neighborhood)?;
/// let highlighted_bonds = selected.to_bonds(BondSelectionMode::Internal);
/// assert_eq!(selected.len(), 2);
/// assert_eq!(highlighted_bonds.len(), 1);
/// let endpoints = highlighted_bonds.to_atoms();
/// assert_eq!(selected, endpoints);
/// # Ok(())
/// # }
/// ```
///
/// Use union/difference/symmetric difference for add/remove/toggle group actions,
/// and [`Self::complement`] to invert. Selection order is dense topology order,
/// not click or query order. Keep an active pick separately when order matters.
/// Membership search is logarithmic in selection size; single-member mutations
/// may shift a linear number of indices. Set algebra merges sorted sets linearly.
///
/// Structural edits publish a different layout: keep selections with their
/// source and rebuild explicitly for a new one. Owning perception and
/// geometry-only changes keep the layout, so selections remain usable. [`Topology::subset`] and model,
/// ensemble, or trajectory slicing consume atom selections as induced subsets;
/// bond selection does not change their structural-subset semantics.
#[derive(Debug, Clone)]
pub struct AtomSelection {
    topology: Arc<Topology>,
    indices: Vec<TopologyAtomIndex>,
}

impl PartialEq for AtomSelection {
    fn eq(&self, other: &Self) -> bool {
        self.topology.shares_layout(&other.topology) && self.indices == other.indices
    }
}

impl Eq for AtomSelection {}

impl AtomSelection {
    /// An empty selection bound to this exact snapshot.
    pub fn empty(topology: impl Into<Arc<Topology>>) -> Self {
        let topology = &topology.into();
        Self {
            topology: Arc::clone(topology),
            indices: Vec::new(),
        }
    }

    /// The number of selected atoms.
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    /// Whether no atom is selected.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    /// Tests membership in this snapshot. Invalid IDs return false.
    pub fn contains(&self, atom: InstanceAtomId) -> bool {
        self.topology
            .atom_index(atom)
            .is_some_and(|index| self.contains_index(index))
    }

    /// Tests dense-index membership. Out-of-range indices return false.
    pub fn contains_index(&self, index: TopologyAtomIndex) -> bool {
        self.indices.binary_search(&index).is_ok()
    }

    /// Selects an atom; returns whether membership changed. Invalid IDs leave
    /// the selection unchanged.
    pub fn insert(&mut self, atom: InstanceAtomId) -> Result<bool, SelectionError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(SelectionError::InvalidAtomId(atom))?;
        Ok(insert_index(&mut self.indices, index))
    }

    /// Deselects an atom; returns whether membership changed. Invalid IDs leave
    /// the selection unchanged.
    pub fn remove(&mut self, atom: InstanceAtomId) -> Result<bool, SelectionError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(SelectionError::InvalidAtomId(atom))?;
        Ok(remove_index(&mut self.indices, index))
    }

    /// Toggles an atom and returns its new selected state. Invalid IDs leave
    /// the selection unchanged.
    pub fn toggle(&mut self, atom: InstanceAtomId) -> Result<bool, SelectionError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(SelectionError::InvalidAtomId(atom))?;
        Ok(toggle_index(&mut self.indices, index))
    }

    /// Removes all members while retaining the topology binding.
    pub fn clear(&mut self) {
        self.indices.clear();
    }

    /// Selects all atoms for which the predicate returns true, in dense order.
    /// The [`AtomView`] exposes chemistry, hierarchy, properties, and stored
    /// perception. No perception is computed implicitly.
    pub fn from_predicate(
        topology: impl Into<Arc<Topology>>,
        mut predicate: impl FnMut(AtomView<'_>) -> bool,
    ) -> Self {
        let topology = &topology.into();
        Self::from_atoms(
            Arc::clone(topology),
            topology
                .atoms()
                .filter(|atom| predicate(*atom))
                .map(AtomView::id),
        )
        .expect("predicate selects validated topology atoms")
    }

    /// Keeps selected atoms satisfying a predicate, preserving dense order.
    pub fn filter(&self, mut predicate: impl FnMut(AtomView<'_>) -> bool) -> Self {
        Self::from_atoms(
            Arc::clone(&self.topology),
            self.indices
                .iter()
                .map(|index| self.topology.atom_at(*index).expect("validated atom"))
                .filter(|atom| predicate(*atom))
                .map(AtomView::id),
        )
        .expect("filter selects validated topology atoms")
    }

    /// Inverts membership within this exact topology.
    pub fn complement(&self) -> Self {
        Self::all(Arc::clone(&self.topology))
            .difference(self)
            .expect("same snapshot")
    }

    /// Keeps atoms present in exactly one selection. Useful for toggling a
    /// picked group; even empty operands must share one layout.
    pub fn symmetric_difference(&self, other: &Self) -> Result<Self, SelectionError> {
        self.combine(other, true, false, true)
    }

    /// Tests whether every selected atom is selected in `other`.
    /// Snapshot compatibility is checked even for empty selections.
    pub fn is_subset(&self, other: &Self) -> Result<bool, SelectionError> {
        self.ensure_compatible(&other.topology)?;
        Ok(self
            .indices
            .iter()
            .all(|index| other.contains_index(*index)))
    }

    /// Tests whether the selections share no atoms, after checking snapshots.
    pub fn is_disjoint(&self, other: &Self) -> Result<bool, SelectionError> {
        self.ensure_compatible(&other.topology)?;
        Ok(self
            .indices
            .iter()
            .all(|index| !other.contains_index(*index)))
    }

    /// Combines two selections from the exact same topology snapshot.
    pub fn union(&self, other: &Self) -> Result<Self, SelectionError> {
        self.combine(other, true, true, true)
    }

    /// Keeps atoms present in both selections, in topology order.
    pub fn intersection(&self, other: &Self) -> Result<Self, SelectionError> {
        self.combine(other, false, true, false)
    }

    /// Keeps atoms in this selection that are absent from `other`.
    pub fn difference(&self, other: &Self) -> Result<Self, SelectionError> {
        self.combine(other, true, false, false)
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

    /// Iterates semantic atom IDs in topology order without allocating.
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = InstanceAtomId> + '_ {
        self.indices.iter().map(|index| {
            self.topology
                .atom_id(*index)
                .expect("validated selection index")
        })
    }

    /// Adds every atom in each touched residue. Selected atoms without a
    /// hierarchy assignment remain selected. This never removes an atom.
    pub fn expand_to_residues(&self) -> Self {
        let residues = self
            .atom_ids()
            .filter_map(|atom| {
                self.topology
                    .atom(atom)
                    .expect("validated atom")
                    .residue()
                    .map(|residue| residue.id())
            })
            .collect::<BTreeSet<_>>();
        let added = self
            .topology
            .atom_sites()
            .filter(|site| residues.contains(&site.residue().id()))
            .map(|site| site.atom().id());
        Self::from_atoms(Arc::clone(&self.topology), self.atom_ids().chain(added))
            .expect("residue expansion uses validated topology atoms")
    }

    /// Adds all atoms in each touched molecule instance (connected component).
    pub fn expand_to_instances(&self) -> Self {
        Self::for_instances(
            Arc::clone(&self.topology),
            self.atom_ids().map(InstanceAtomId::instance),
        )
        .expect("selected atoms have valid instances")
    }

    /// Adds all atoms in each touched chain, retaining atoms without hierarchy.
    /// Chain membership is independent of covalent connectedness.
    pub fn expand_to_chains(&self) -> Self {
        let chains = self.atom_ids().filter_map(|atom| {
            self.topology
                .atom(atom)
                .expect("validated atom")
                .chain()
                .map(|chain| chain.id())
        });
        self.union(&Self::for_chains(Arc::clone(&self.topology), chains).expect("valid chains"))
            .expect("same snapshot")
    }

    /// Adds atoms at most `steps` asserted bonds away from any selected atom.
    /// Zero steps returns the original set. Traversal never crosses instances
    /// or infers spatial contacts; selected atoms are always retained.
    pub fn expand_bonded(&self, steps: usize) -> Self {
        let mut visited = vec![false; self.topology.atom_count()];
        let mut frontier = self.atom_ids().collect::<Vec<_>>();
        for index in &self.indices {
            visited[index.index()] = true;
        }
        for _ in 0..steps {
            if frontier.is_empty() {
                break;
            }
            let mut next = Vec::new();
            for atom in frontier {
                for neighbor in self
                    .topology
                    .atom(atom)
                    .expect("validated atom")
                    .neighbors()
                {
                    let index = neighbor.index().index();
                    if !visited[index] {
                        visited[index] = true;
                        next.push(neighbor.id());
                    }
                }
            }
            frontier = next;
        }
        Self {
            topology: Arc::clone(&self.topology),
            indices: visited
                .into_iter()
                .enumerate()
                .filter(|(_, selected)| *selected)
                .map(|(index, _)| TopologyAtomIndex::new(index as u32))
                .collect(),
        }
    }

    /// Selects bonds according to explicit endpoint membership. This does not
    /// modify atom membership or perform a structural subset operation.
    pub fn to_bonds(&self, mode: BondSelectionMode) -> BondSelection {
        BondSelection::from_predicate(Arc::clone(&self.topology), |bond| {
            let [a, b] = bond.atoms();
            let a = self.indices.binary_search(&a.index()).is_ok();
            let b = self.indices.binary_search(&b.index()).is_ok();
            match mode {
                BondSelectionMode::Internal => a && b,
                BondSelectionMode::Incident => a || b,
                BondSelectionMode::Boundary => a != b,
            }
        })
    }

    /// Selects sites with one of these exact label atom names. Missing names
    /// do not match. Compose with residue-class selections for biological roles.
    pub fn for_label_atom_names<'a>(
        topology: impl Into<Arc<Topology>>,
        names: impl IntoIterator<Item = &'a str>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let names = names.into_iter().collect::<BTreeSet<_>>();
        Self::from_atoms(
            Arc::clone(topology),
            topology
                .atom_sites()
                .filter(|site| {
                    site.metadata()
                        .label_atom_id
                        .as_deref()
                        .is_some_and(|name| names.contains(name))
                })
                .map(|site| site.atom().id()),
        )
    }

    /// Selects sites with one of these exact author atom names; never falls
    /// back to label names. Missing names do not match.
    pub fn for_author_atom_names<'a>(
        topology: impl Into<Arc<Topology>>,
        names: impl IntoIterator<Item = &'a str>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let names = names.into_iter().collect::<BTreeSet<_>>();
        Self::from_atoms(
            Arc::clone(topology),
            topology
                .atom_sites()
                .filter(|site| {
                    site.metadata()
                        .auth_atom_id
                        .as_deref()
                        .is_some_and(|name| names.contains(name))
                })
                .map(|site| site.atom().id()),
        )
    }

    /// Selects every atom in authoritative dense order, sharing this exact topology.
    ///
    /// This is infallible because the topology has already validated its layout.
    pub fn all(topology: impl Into<Arc<Topology>>) -> Self {
        let topology = &topology.into();
        Self {
            topology: Arc::clone(topology),
            indices: (0..topology.atom_count())
                .map(|index| TopologyAtomIndex::new(index as u32))
                .collect(),
        }
    }

    /// Selects the given qualified atoms; order and duplicates are ignored.
    ///
    /// # Errors
    ///
    /// [`SelectionError::InvalidAtomId`] when an atom is not in `topology`.
    pub fn from_atoms(
        topology: impl Into<Arc<Topology>>,
        atoms: impl IntoIterator<Item = InstanceAtomId>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let mut indices = atoms
            .into_iter()
            .map(|atom| {
                topology
                    .atom_index(atom)
                    .ok_or(SelectionError::InvalidAtomId(atom))
            })
            .collect::<Result<Vec<_>, _>>()?;
        indices.sort_unstable();
        indices.dedup();
        Ok(Self {
            topology: Arc::clone(topology),
            indices,
        })
    }

    /// Selects the atoms at the given dense indices; order and duplicates are
    /// ignored.
    ///
    /// Use [`Self::from_sorted_indices`] when the indices are already sorted
    /// and unique, or [`Self::from_mask`] for a bit mask.
    ///
    /// # Errors
    ///
    /// [`SelectionError::InvalidAtomIndex`] when an index is not below the
    /// topology's atom count.
    pub fn from_indices(
        topology: impl Into<Arc<Topology>>,
        indices: impl IntoIterator<Item = TopologyAtomIndex>,
    ) -> Result<Self, SelectionError> {
        let topology = topology.into();
        let count = topology.atom_count();
        let mut indices = indices
            .into_iter()
            .map(|index| {
                if index.index() < count {
                    Ok(index)
                } else {
                    Err(SelectionError::InvalidAtomIndex(index))
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        indices.sort_unstable();
        indices.dedup();
        Ok(Self { topology, indices })
    }

    /// Selects the atoms at dense indices that are already strictly
    /// increasing, keeping the vector as the selection's storage.
    ///
    /// The indices are checked in one pass, without sorting or translating
    /// them to atom IDs, so a host that walks dense order (for example a
    /// selection language resolver) pays time linear in the selection and no
    /// further allocation.
    ///
    /// # Errors
    ///
    /// [`SelectionError::UnsortedAtomIndex`] at the first index that is not
    /// greater than the one before it, and otherwise
    /// [`SelectionError::InvalidAtomIndex`] at the first index that is not
    /// below the topology's atom count.
    pub fn from_sorted_indices(
        topology: impl Into<Arc<Topology>>,
        indices: Vec<TopologyAtomIndex>,
    ) -> Result<Self, SelectionError> {
        let topology = topology.into();
        // A branch-free pass vectorizes; the failing index is found only on error.
        let increasing = indices
            .windows(2)
            .fold(true, |increasing, pair| increasing & (pair[0] < pair[1]));
        if !increasing {
            let pair = indices
                .windows(2)
                .find(|pair| pair[0] >= pair[1])
                .expect("an index is out of order");
            return Err(SelectionError::UnsortedAtomIndex(pair[1]));
        }
        let count = topology.atom_count();
        if indices.last().is_some_and(|last| last.index() >= count) {
            let first = indices.partition_point(|index| index.index() < count);
            return Err(SelectionError::InvalidAtomIndex(indices[first]));
        }
        Ok(Self { topology, indices })
    }

    /// Selects the atoms whose bits are set in a mask over dense atom order.
    ///
    /// Bit `i % 64` of word `i / 64` (least significant bit first) selects the
    /// atom at dense index `i`. The mask has exactly
    /// `atom_count.div_ceil(64)` words, and the bits past the last atom in
    /// its final word are clear. [`Self::to_mask`] produces this layout.
    ///
    /// # Errors
    ///
    /// [`SelectionError::AtomMaskLength`] when the mask has a different
    /// number of words, or [`SelectionError::InvalidAtomIndex`] at the first
    /// set bit past the last atom.
    pub fn from_mask(
        topology: impl Into<Arc<Topology>>,
        mask: &[u64],
    ) -> Result<Self, SelectionError> {
        let topology = topology.into();
        let count = topology.atom_count();
        let expected = count.div_ceil(64);
        if mask.len() != expected {
            return Err(SelectionError::AtomMaskLength {
                expected,
                found: mask.len(),
            });
        }
        let tail = count % 64;
        if let Some(&last) = mask.last().filter(|_| tail != 0) {
            let padding = last >> tail;
            if padding != 0 {
                let index = count + padding.trailing_zeros() as usize;
                let index = u32::try_from(index).unwrap_or(u32::MAX);
                return Err(SelectionError::InvalidAtomIndex(TopologyAtomIndex::new(
                    index,
                )));
            }
        }
        let selected = mask.iter().map(|word| word.count_ones() as usize).sum();
        let mut indices = Vec::with_capacity(selected);
        for (position, &word) in mask.iter().enumerate() {
            // Every set bit is a dense index, so its word's base fits in u32.
            let base = (position * 64) as u32;
            if word == u64::MAX {
                indices.extend((base..base + 64).map(TopologyAtomIndex::new));
                continue;
            }
            let mut word = word;
            while word != 0 {
                indices.push(TopologyAtomIndex::new(base + word.trailing_zeros()));
                word &= word - 1;
            }
        }
        Ok(Self { topology, indices })
    }

    /// The selection as a bit mask over dense atom order, in the layout
    /// [`Self::from_mask`] accepts.
    pub fn to_mask(&self) -> Vec<u64> {
        let mut mask = vec![0; self.topology.atom_count().div_ceil(64)];
        // Bits gather in a register and are stored once per word.
        let (mut position, mut word) = (0, 0u64);
        for index in &self.indices {
            let index = index.index();
            if index / 64 != position {
                mask[position] = word;
                (position, word) = (index / 64, 0);
            }
            word |= 1 << (index % 64);
        }
        if word != 0 {
            mask[position] = word;
        }
        mask
    }

    /// The selection as maximal runs of consecutive dense indices, in order.
    pub fn ranges(&self) -> impl Iterator<Item = Range<usize>> + '_ {
        let mut indices = self.indices.iter().map(|index| index.index()).peekable();
        std::iter::from_fn(move || {
            let start = indices.next()?;
            let mut end = start + 1;
            while indices.next_if_eq(&end).is_some() {
                end += 1;
            }
            Some(start..end)
        })
    }

    /// Selects every atom of the given molecule instances.
    ///
    /// # Errors
    ///
    /// [`SelectionError::InvalidMoleculeInstanceId`] when an instance is not in
    /// `topology`.
    pub fn for_instances(
        topology: impl Into<Arc<Topology>>,
        instances: impl IntoIterator<Item = MoleculeInstanceId>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let instances = instances.into_iter().collect::<BTreeSet<_>>();
        for instance in &instances {
            topology
                .instance(*instance)
                .ok_or(SelectionError::InvalidMoleculeInstanceId(*instance))?;
        }
        Self::from_atoms(
            Arc::clone(topology),
            topology
                .atom_ids()
                .iter()
                .copied()
                .filter(|atom| instances.contains(&atom.instance())),
        )
    }

    /// Selects every atom of every instance of the given definitions.
    ///
    /// # Errors
    ///
    /// [`SelectionError::InvalidMoleculeDefinitionId`] when a definition is not
    /// in `topology`.
    pub fn for_definitions(
        topology: impl Into<Arc<Topology>>,
        definitions: impl IntoIterator<Item = MoleculeDefinitionId>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let definitions = definitions.into_iter().collect::<BTreeSet<_>>();
        for definition in &definitions {
            topology
                .definition(*definition)
                .ok_or(SelectionError::InvalidMoleculeDefinitionId(*definition))?;
        }
        let instances = topology
            .instances()
            .filter(|molecule| definitions.contains(&molecule.definition_id()))
            .map(|molecule| molecule.id());
        Self::for_instances(Arc::clone(topology), instances)
    }

    /// Selects all atoms in molecule instances having one of `classes`.
    ///
    /// # Errors
    ///
    /// Currently never fails; every selected atom comes from `topology` itself.
    pub fn for_molecule_classes(
        topology: impl Into<Arc<Topology>>,
        classes: impl IntoIterator<Item = MoleculeClass>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let classes = classes.into_iter().collect::<BTreeSet<_>>();
        Self::for_instances(
            Arc::clone(topology),
            topology
                .instances()
                .filter(|molecule| classes.contains(&molecule.class()))
                .map(|molecule| molecule.id()),
        )
    }

    /// Selects every atom whose element is one of `elements`.
    ///
    /// # Errors
    ///
    /// Currently never fails; every selected atom comes from `topology` itself.
    pub fn for_elements(
        topology: impl Into<Arc<Topology>>,
        elements: impl IntoIterator<Item = Element>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let elements = elements.into_iter().collect::<BTreeSet<_>>();
        Self::from_atoms(
            Arc::clone(topology),
            topology
                .atoms()
                .filter(|atom| elements.contains(&atom.element))
                .map(AtomView::id),
        )
    }

    /// Selects the atoms placed by the given hierarchy atom sites.
    ///
    /// # Errors
    ///
    /// [`SelectionError::InvalidAtomSiteId`] when a site is not in `topology`.
    pub fn for_atom_sites(
        topology: impl Into<Arc<Topology>>,
        atom_sites: impl IntoIterator<Item = AtomSiteId>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let atom_sites = atom_sites.into_iter().collect::<BTreeSet<_>>();
        let atoms = atom_sites
            .into_iter()
            .map(|site| {
                topology
                    .atom_site(site)
                    .map(|site| site.atom().id())
                    .ok_or(SelectionError::InvalidAtomSiteId(site))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_atoms(Arc::clone(topology), atoms)
    }

    /// Selects the atoms placed in the given residues.
    ///
    /// # Errors
    ///
    /// [`SelectionError::InvalidResidueId`] when a residue is not in `topology`.
    pub fn for_residues(
        topology: impl Into<Arc<Topology>>,
        residues: impl IntoIterator<Item = ResidueId>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let residues = residues.into_iter().collect::<BTreeSet<_>>();
        for residue in &residues {
            topology
                .residue(*residue)
                .ok_or(SelectionError::InvalidResidueId(*residue))?;
        }
        Self::for_atom_sites(
            Arc::clone(topology),
            topology
                .atom_sites()
                .filter(|site| residues.contains(&site.residue().id()))
                .map(AtomSiteView::id),
        )
    }

    /// Selects all atoms in residues having one of `classes`.
    ///
    /// # Errors
    ///
    /// Currently never fails; every selected atom comes from `topology` itself.
    pub fn for_residue_classes(
        topology: impl Into<Arc<Topology>>,
        classes: impl IntoIterator<Item = ResidueClass>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let classes = classes.into_iter().collect::<BTreeSet<_>>();
        Self::for_residues(
            Arc::clone(topology),
            topology
                .residues()
                .filter(|residue| classes.contains(&residue.class()))
                .map(ResidueView::id),
        )
    }

    /// Selects the atoms placed in the residues of the given chains.
    ///
    /// # Errors
    ///
    /// [`SelectionError::InvalidChainId`] when a chain is not in `topology`.
    pub fn for_chains(
        topology: impl Into<Arc<Topology>>,
        chains: impl IntoIterator<Item = ChainId>,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        let chains = chains.into_iter().collect::<BTreeSet<_>>();
        for chain in &chains {
            topology
                .chain(*chain)
                .ok_or(SelectionError::InvalidChainId(*chain))?;
        }
        Self::for_residues(
            Arc::clone(topology),
            topology
                .residues()
                .filter(|residue| chains.contains(&residue.chain().id()))
                .map(ResidueView::id),
        )
    }

    /// Selects the atoms of every chain whose label identifier equals `label`.
    ///
    /// An unmatched label selects nothing.
    ///
    /// # Errors
    ///
    /// Currently never fails; every selected atom comes from `topology` itself.
    pub fn for_chain_label(
        topology: impl Into<Arc<Topology>>,
        label: &str,
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        Self::for_chains(
            Arc::clone(topology),
            topology
                .chains()
                .filter(|chain| chain.label_id() == label)
                .map(ChainView::id),
        )
    }

    /// Unions molecule-local matches qualified with the supplied instance.
    /// The caller must match against that instance's molecule: `QueryMatch`
    /// carries local IDs, so target provenance cannot be checked here.
    /// Use [`Self::from_topology_query_matches`] for snapshot-checked matches.
    ///
    /// # Errors
    ///
    /// [`SelectionError::InvalidMoleculeInstanceId`] when `instance` is not in
    /// `topology`, or [`SelectionError::InvalidAtomId`] when a matched atom is not
    /// in that instance.
    pub fn from_query_matches(
        topology: impl Into<Arc<Topology>>,
        instance: MoleculeInstanceId,
        matches: &[QueryMatch],
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        topology
            .instance(instance)
            .ok_or(SelectionError::InvalidMoleculeInstanceId(instance))?;
        Self::from_atoms(
            Arc::clone(topology),
            matches.iter().flat_map(|query_match| {
                query_match
                    .atoms()
                    .iter()
                    .copied()
                    .map(move |atom| InstanceAtomId::new(instance, atom))
            }),
        )
    }

    /// Unions all atoms in topology query matches after validating every
    /// match's exact snapshot. An empty match list produces a bound empty set.
    /// Query ordering, tags, and individual match boundaries are intentionally
    /// discarded; retain the original matches when those associations matter.
    ///
    /// # Errors
    ///
    /// [`SelectionError::TopologyMismatch`] when a match was made against a
    /// topology that does not share this layout.
    pub fn from_topology_query_matches(
        topology: impl Into<Arc<Topology>>,
        matches: &[TopologyQueryMatch],
    ) -> Result<Self, SelectionError> {
        let topology = &topology.into();
        for matched in matches {
            if !topology.shares_layout(matched.topology()) {
                return Err(SelectionError::TopologyMismatch);
            }
        }
        Self::from_atoms(
            Arc::clone(topology),
            matches
                .iter()
                .flat_map(|matched| matched.atoms().iter().copied()),
        )
    }

    /// Checks that `topology` shares this selection's layout
    /// ([`Topology::shares_layout`]), so its dense indices address the same atoms.
    ///
    /// # Errors
    ///
    /// [`SelectionError::TopologyMismatch`] when the layouts differ.
    pub fn ensure_compatible(&self, topology: &Topology) -> Result<(), SelectionError> {
        if !self.topology.shares_layout(topology) {
            return Err(SelectionError::TopologyMismatch);
        }
        Ok(())
    }

    /// The topology snapshot this selection is bound to.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// Clones the handle to the exact snapshot, without copying topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// The selected dense atom indices, sorted and unique.
    pub fn indices(&self) -> &[TopologyAtomIndex] {
        &self.indices
    }

    /// The selected qualified atom IDs in dense order, resolved in `topology`.
    ///
    /// # Errors
    ///
    /// [`SelectionError::TopologyMismatch`] when `topology` does not share this
    /// selection's layout.
    pub fn semantic_ids(&self, topology: &Topology) -> Result<Vec<InstanceAtomId>, SelectionError> {
        self.ensure_compatible(topology)?;
        Ok(self
            .indices
            .iter()
            .map(|index| {
                topology
                    .atom_id(*index)
                    .expect("selection contains validated dense indices")
            })
            .collect())
    }
}

/// An atom or bond selection references something outside its topology.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SelectionError {
    /// The selection belongs to a different topology layout.
    TopologyMismatch,
    /// A molecule definition is not in the topology.
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    /// A molecule instance is not in the topology.
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    /// An atom is not in the topology.
    InvalidAtomId(InstanceAtomId),
    /// A chain is not in the topology's hierarchy.
    InvalidChainId(ChainId),
    /// A residue is not in the topology's hierarchy.
    InvalidResidueId(ResidueId),
    /// An atom site is not in the topology's hierarchy.
    InvalidAtomSiteId(AtomSiteId),
    /// A dense atom index is outside the topology.
    InvalidAtomIndex(TopologyAtomIndex),
    /// A dense atom index is not greater than the one before it, in input
    /// that must be strictly increasing.
    UnsortedAtomIndex(TopologyAtomIndex),
    /// An atom mask does not have one bit per atom, rounded up to whole words.
    AtomMaskLength {
        /// The number of words the topology's atoms need.
        expected: usize,
        /// The number of words supplied.
        found: usize,
    },
    /// A bond is not in the topology.
    InvalidBondId(InstanceBondId),
    /// A dense bond index is outside the topology.
    InvalidBondIndex(TopologyBondIndex),
}

impl fmt::Display for SelectionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("selection belongs to a different topology")
            }
            Self::InvalidMoleculeDefinitionId(id) => {
                write!(formatter, "invalid selected molecule definition: {id}")
            }
            Self::InvalidMoleculeInstanceId(id) => {
                write!(formatter, "invalid selected molecule instance: {id}")
            }
            Self::InvalidAtomId(id) => write!(formatter, "invalid selected atom: {id}"),
            Self::InvalidChainId(id) => write!(formatter, "invalid selected chain: {id}"),
            Self::InvalidResidueId(id) => write!(formatter, "invalid selected residue: {id}"),
            Self::InvalidAtomSiteId(id) => write!(formatter, "invalid selected atom site: {id}"),
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid selected {index}"),
            Self::UnsortedAtomIndex(index) => {
                write!(formatter, "selected {index} is out of order or repeated")
            }
            Self::AtomMaskLength { expected, found } => write!(
                formatter,
                "atom mask has {found} words; the topology needs {expected}"
            ),
            Self::InvalidBondId(id) => write!(formatter, "invalid selected bond: {id}"),
            Self::InvalidBondIndex(index) => write!(formatter, "invalid selected {index}"),
        }
    }
}

impl std::error::Error for SelectionError {}

fn insert_index<T: Ord>(indices: &mut Vec<T>, index: T) -> bool {
    match indices.binary_search(&index) {
        Ok(_) => false,
        Err(position) => {
            indices.insert(position, index);
            true
        }
    }
}

fn remove_index<T: Ord>(indices: &mut Vec<T>, index: T) -> bool {
    match indices.binary_search(&index) {
        Ok(position) => {
            indices.remove(position);
            true
        }
        Err(_) => false,
    }
}

fn toggle_index<T: Ord>(indices: &mut Vec<T>, index: T) -> bool {
    match indices.binary_search(&index) {
        Ok(position) => {
            indices.remove(position);
            false
        }
        Err(position) => {
            indices.insert(position, index);
            true
        }
    }
}

fn combine_indices<T: Ord + Copy>(
    left: &[T],
    right: &[T],
    left_only: bool,
    both: bool,
    right_only: bool,
) -> Vec<T> {
    let mut left = left.iter().copied().peekable();
    let mut right = right.iter().copied().peekable();
    let mut indices = Vec::new();
    while let (Some(&a), Some(&b)) = (left.peek(), right.peek()) {
        match a.cmp(&b) {
            std::cmp::Ordering::Less => {
                left.next();
                if left_only {
                    indices.push(a);
                }
            }
            std::cmp::Ordering::Equal => {
                left.next();
                right.next();
                if both {
                    indices.push(a);
                }
            }
            std::cmp::Ordering::Greater => {
                right.next();
                if right_only {
                    indices.push(b);
                }
            }
        }
    }
    if left_only {
        indices.extend(left);
    }
    if right_only {
        indices.extend(right);
    }
    indices
}
