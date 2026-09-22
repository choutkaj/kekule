use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use crate::core::{Atom, Element};
use crate::substructure::{QueryMatch, TopologyQueryMatch};

use super::{
    AtomSiteId, AtomSiteView, ChainId, ChainView, InstanceAtomId, InstanceBondId, MoleculeClass,
    MoleculeDefinitionId, MoleculeInstanceId, ResidueClass, ResidueId, ResidueView, Topology,
    TopologyAtomIndex, TopologyBondIndex,
};

mod bonds;
pub use bonds::{BondSelection, BondSelectionMode};

/// A topology-bound, sorted, unique dense atom selection.
///
/// Empty selections retain their topology. Equality and set operations use exact
/// shared snapshot identity, not chemical or layout equality. IDs supplied to
/// constructors and membership methods are interpreted in this snapshot: bare
/// IDs do not carry provenance. Selections are static sets, not stored queries.
/// See [`BondSelection`] for independently selected bonds.
///
/// # Selection sources
///
/// | Source | Constructors |
/// | --- | --- |
/// | Picking or dense masks | [`Self::from_atoms`], [`Self::from_indices`] |
/// | All/none | [`Self::all`], [`Self::empty`] |
/// | Molecule occurrences or reused definitions | [`Self::for_instances`], [`Self::for_definitions`] |
/// | Chemistry or custom properties | [`Self::for_elements`], [`Self::for_molecule_classes`], [`Self::from_predicate`] |
/// | Hierarchy | [`Self::for_chains`], [`Self::for_chain_label`], [`Self::for_residues`], [`Self::for_residue_classes`], [`Self::for_atom_sites`] |
/// | Exact names in separate namespaces | [`Self::for_label_atom_names`], [`Self::for_author_atom_names`] |
/// | Substructure matches | [`Self::from_query_matches`], [`Self::from_topology_query_matches`] |
/// | Cartesian distance in one realization | [`crate::structure::measure::within`] |
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
/// let topology = Arc::new(smiles::to_topology("CCO")?);
/// let mut selected = AtomSelection::empty(&topology);
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
/// Structural edits and owning perception publish a different snapshot: keep
/// selections with their source snapshot and rebuild explicitly for a new one.
/// Geometry-only changes retain the binding. [`Topology::subset`] and model,
/// ensemble, or trajectory slicing consume atom selections as induced subsets;
/// bond selection does not change their structural-subset semantics.
#[derive(Debug, Clone)]
pub struct AtomSelection {
    topology: Arc<Topology>,
    indices: Vec<TopologyAtomIndex>,
}

impl PartialEq for AtomSelection {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.indices == other.indices
    }
}

impl Eq for AtomSelection {}

impl AtomSelection {
    /// An empty selection bound to this exact snapshot.
    pub fn empty(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            indices: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.indices.len()
    }

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
    /// The semantic ID allows inspecting hierarchy, properties, or perception
    /// through the supplied topology. No perception is computed implicitly.
    pub fn from_predicate(
        topology: &Arc<Topology>,
        mut predicate: impl FnMut(InstanceAtomId, &Atom) -> bool,
    ) -> Self {
        Self::from_atoms(
            topology,
            topology
                .atoms()
                .filter(|(id, atom)| predicate(*id, atom))
                .map(|(id, _)| id),
        )
        .expect("predicate selects validated topology atoms")
    }

    /// Keeps selected atoms satisfying a predicate, preserving dense order.
    pub fn filter(&self, mut predicate: impl FnMut(InstanceAtomId, &Atom) -> bool) -> Self {
        Self::from_atoms(
            &self.topology,
            self.atom_ids()
                .filter(|id| predicate(*id, self.topology.atom(*id).expect("validated atom"))),
        )
        .expect("filter selects validated topology atoms")
    }

    /// Inverts membership within this exact topology.
    pub fn complement(&self) -> Self {
        Self::all(&self.topology)
            .difference(self)
            .expect("same snapshot")
    }

    /// Keeps atoms present in exactly one selection. Useful for toggling a
    /// picked group; even empty operands must share the exact snapshot.
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
                    .hierarchy()
                    .atom_site_for_atom(atom)
                    .map(|site| site.residue())
            })
            .collect::<BTreeSet<_>>();
        let added = self
            .topology
            .atom_sites()
            .filter(|site| residues.contains(&site.residue().id()))
            .map(|site| site.atom());
        Self::from_atoms(&self.topology, self.atom_ids().chain(added))
            .expect("residue expansion uses validated topology atoms")
    }

    /// Adds all atoms in each touched molecule instance (connected component).
    pub fn expand_to_instances(&self) -> Self {
        Self::for_instances(
            &self.topology,
            self.atom_ids().map(InstanceAtomId::molecule),
        )
        .expect("selected atoms have valid instances")
    }

    /// Adds all atoms in each touched chain, retaining atoms without hierarchy.
    /// Chain membership is independent of covalent connectedness.
    pub fn expand_to_chains(&self) -> Self {
        let chains = self.atom_ids().filter_map(|atom| {
            self.topology
                .atom_site_for_atom(atom)
                .expect("validated atom")
                .map(|site| site.residue().chain().id())
        });
        self.union(&Self::for_chains(&self.topology, chains).expect("valid chains"))
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
                for neighbor in self.topology.neighbors(atom).expect("validated atom") {
                    let index = self
                        .topology
                        .atom_index(neighbor)
                        .expect("validated neighbor")
                        .index();
                    if !visited[index] {
                        visited[index] = true;
                        next.push(neighbor);
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
        BondSelection::from_predicate(&self.topology, |id, bond| {
            let a = self.contains(InstanceAtomId::new(id.molecule(), bond.a()));
            let b = self.contains(InstanceAtomId::new(id.molecule(), bond.b()));
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
        topology: &Arc<Topology>,
        names: impl IntoIterator<Item = &'a str>,
    ) -> Result<Self, SelectionError> {
        let names = names.into_iter().collect::<BTreeSet<_>>();
        Self::from_atoms(
            topology,
            topology
                .atom_sites()
                .filter(|site| {
                    site.metadata()
                        .label_atom_id
                        .as_deref()
                        .is_some_and(|name| names.contains(name))
                })
                .map(AtomSiteView::atom),
        )
    }

    /// Selects sites with one of these exact author atom names; never falls
    /// back to label names. Missing names do not match.
    pub fn for_author_atom_names<'a>(
        topology: &Arc<Topology>,
        names: impl IntoIterator<Item = &'a str>,
    ) -> Result<Self, SelectionError> {
        let names = names.into_iter().collect::<BTreeSet<_>>();
        Self::from_atoms(
            topology,
            topology
                .atom_sites()
                .filter(|site| {
                    site.metadata()
                        .auth_atom_id
                        .as_deref()
                        .is_some_and(|name| names.contains(name))
                })
                .map(AtomSiteView::atom),
        )
    }

    /// Selects every atom in authoritative dense order, sharing this exact topology.
    ///
    /// This is infallible because the topology has already validated its layout.
    pub fn all(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            indices: (0..topology.atom_count())
                .map(|index| TopologyAtomIndex::new(index as u32))
                .collect(),
        }
    }

    pub fn from_atoms(
        topology: &Arc<Topology>,
        atoms: impl IntoIterator<Item = InstanceAtomId>,
    ) -> Result<Self, SelectionError> {
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

    pub fn from_indices(
        topology: &Arc<Topology>,
        indices: impl IntoIterator<Item = TopologyAtomIndex>,
    ) -> Result<Self, SelectionError> {
        let atoms = indices
            .into_iter()
            .map(|index| {
                topology
                    .atom_id(index)
                    .ok_or(SelectionError::InvalidAtomIndex(index))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_atoms(topology, atoms)
    }

    pub fn for_instances(
        topology: &Arc<Topology>,
        instances: impl IntoIterator<Item = MoleculeInstanceId>,
    ) -> Result<Self, SelectionError> {
        let instances = instances.into_iter().collect::<BTreeSet<_>>();
        for instance in &instances {
            topology
                .instance(*instance)
                .map_err(|_| SelectionError::InvalidMoleculeInstanceId(*instance))?;
        }
        Self::from_atoms(
            topology,
            topology
                .atom_ids()
                .iter()
                .copied()
                .filter(|atom| instances.contains(&atom.molecule())),
        )
    }

    pub fn for_definitions(
        topology: &Arc<Topology>,
        definitions: impl IntoIterator<Item = MoleculeDefinitionId>,
    ) -> Result<Self, SelectionError> {
        let definitions = definitions.into_iter().collect::<BTreeSet<_>>();
        for definition in &definitions {
            topology
                .definition(*definition)
                .map_err(|_| SelectionError::InvalidMoleculeDefinitionId(*definition))?;
        }
        let instances = topology
            .instances()
            .filter(|(_, instance)| definitions.contains(&instance.definition()))
            .map(|(id, _)| id);
        Self::for_instances(topology, instances)
    }

    /// Selects all atoms in molecule instances having one of `classes`.
    pub fn for_molecule_classes(
        topology: &Arc<Topology>,
        classes: impl IntoIterator<Item = MoleculeClass>,
    ) -> Result<Self, SelectionError> {
        let classes = classes.into_iter().collect::<BTreeSet<_>>();
        Self::for_instances(
            topology,
            topology
                .molecules()
                .filter(|molecule| classes.contains(&molecule.class()))
                .map(|molecule| molecule.id()),
        )
    }

    pub fn for_elements(
        topology: &Arc<Topology>,
        elements: impl IntoIterator<Item = Element>,
    ) -> Result<Self, SelectionError> {
        let elements = elements.into_iter().collect::<BTreeSet<_>>();
        Self::from_atoms(
            topology,
            topology
                .atoms()
                .filter(|(_, atom)| elements.contains(&atom.element))
                .map(|(id, _)| id),
        )
    }

    pub fn for_atom_sites(
        topology: &Arc<Topology>,
        atom_sites: impl IntoIterator<Item = AtomSiteId>,
    ) -> Result<Self, SelectionError> {
        let atom_sites = atom_sites.into_iter().collect::<BTreeSet<_>>();
        let atoms = atom_sites
            .into_iter()
            .map(|site| {
                topology
                    .atom_for_site(site)
                    .map_err(|_| SelectionError::InvalidAtomSiteId(site))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_atoms(topology, atoms)
    }

    pub fn for_residues(
        topology: &Arc<Topology>,
        residues: impl IntoIterator<Item = ResidueId>,
    ) -> Result<Self, SelectionError> {
        let residues = residues.into_iter().collect::<BTreeSet<_>>();
        for residue in &residues {
            topology
                .residue(*residue)
                .map_err(|_| SelectionError::InvalidResidueId(*residue))?;
        }
        Self::for_atom_sites(
            topology,
            topology
                .atom_sites()
                .filter(|site| residues.contains(&site.residue().id()))
                .map(AtomSiteView::id),
        )
    }

    /// Selects all atoms in residues having one of `classes`.
    pub fn for_residue_classes(
        topology: &Arc<Topology>,
        classes: impl IntoIterator<Item = ResidueClass>,
    ) -> Result<Self, SelectionError> {
        let classes = classes.into_iter().collect::<BTreeSet<_>>();
        Self::for_residues(
            topology,
            topology
                .residues()
                .filter(|residue| classes.contains(&residue.class()))
                .map(ResidueView::id),
        )
    }

    pub fn for_chains(
        topology: &Arc<Topology>,
        chains: impl IntoIterator<Item = ChainId>,
    ) -> Result<Self, SelectionError> {
        let chains = chains.into_iter().collect::<BTreeSet<_>>();
        for chain in &chains {
            topology
                .chain(*chain)
                .map_err(|_| SelectionError::InvalidChainId(*chain))?;
        }
        Self::for_residues(
            topology,
            topology
                .residues()
                .filter(|residue| chains.contains(&residue.chain().id()))
                .map(ResidueView::id),
        )
    }

    pub fn for_chain_label(topology: &Arc<Topology>, label: &str) -> Result<Self, SelectionError> {
        Self::for_chains(
            topology,
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
    pub fn from_query_matches(
        topology: &Arc<Topology>,
        instance: MoleculeInstanceId,
        matches: &[QueryMatch],
    ) -> Result<Self, SelectionError> {
        topology
            .instance(instance)
            .map_err(|_| SelectionError::InvalidMoleculeInstanceId(instance))?;
        Self::from_atoms(
            topology,
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
    pub fn from_topology_query_matches(
        topology: &Arc<Topology>,
        matches: &[TopologyQueryMatch],
    ) -> Result<Self, SelectionError> {
        for matched in matches {
            if !Arc::ptr_eq(topology, matched.topology()) {
                return Err(SelectionError::TopologyMismatch);
            }
        }
        Self::from_atoms(
            topology,
            matches
                .iter()
                .flat_map(|matched| matched.atoms().iter().copied()),
        )
    }

    pub fn ensure_compatible(&self, topology: &Arc<Topology>) -> Result<(), SelectionError> {
        if !Arc::ptr_eq(&self.topology, topology) {
            return Err(SelectionError::TopologyMismatch);
        }
        Ok(())
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// Clones the handle to the exact snapshot, without copying topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn indices(&self) -> &[TopologyAtomIndex] {
        &self.indices
    }

    pub fn semantic_ids(
        &self,
        topology: &Arc<Topology>,
    ) -> Result<Vec<InstanceAtomId>, SelectionError> {
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

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SelectionError {
    TopologyMismatch,
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    InvalidAtomId(InstanceAtomId),
    InvalidChainId(ChainId),
    InvalidResidueId(ResidueId),
    InvalidAtomSiteId(AtomSiteId),
    InvalidAtomIndex(TopologyAtomIndex),
    InvalidBondId(InstanceBondId),
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
