use crate::diagnostics::{CapacityExceeded, Countable};
use std::fmt;
use std::sync::Arc;

use super::dense_map::QualifiedMap;
use super::shared_rows::{SharedRow, SharedRows};

use super::{InstanceAtomId, ResidueClass};

fixed_u32_id!(
    /// Stable topology-global semantic identifier of one hierarchy chain;
    /// [`Topology::chain`](super::Topology::chain) resolves it.
    ///
    /// Not a dense coordinate index. Chains are numbered in creation order, and
    /// the ID also selects the chain's row in topology chain properties.
    ChainId,
    "chain"
);
fixed_u32_id!(
    /// Stable topology-global semantic identifier of one hierarchy residue;
    /// [`Topology::residue`](super::Topology::residue) resolves it.
    ///
    /// Not a dense coordinate index. Residues are numbered in creation order,
    /// and the ID also selects the residue's row in topology residue properties.
    ResidueId,
    "residue"
);
fixed_u32_id!(
    /// Stable topology-global semantic identifier of one hierarchy atom site;
    /// [`Topology::atom_site`](super::Topology::atom_site) resolves it.
    ///
    /// Not a dense atom index: sites are numbered in creation order,
    /// independently of [`TopologyAtomIndex`](super::TopologyAtomIndex).
    AtomSiteId,
    "atom site"
);

/// Topology-owned storage of coordinate-independent residue, chain, and
/// atom-site organization.
///
/// A hierarchy owns no atoms or bonds. Atom sites refer to topology-global
/// [`InstanceAtomId`] values and may span molecule-instance boundaries. It is
/// valid for a topology to have an empty hierarchy. Reads go through the
/// topology views ([`super::ChainView`], [`super::ResidueView`],
/// [`super::AtomSiteView`]); staging goes through [`HierarchyBuilder`].
///
/// Residues and atom sites are stored by column. Residue labels and site
/// metadata are keys into row tables shared by every snapshot and draft that
/// keeps them, so republishing a large hierarchy copies integers.
#[derive(Debug, Clone, Default)]
pub(crate) struct Hierarchy {
    chain_label_ids: Vec<String>,
    chain_author_ids: Vec<Option<String>>,
    chain_residues: Vec<Vec<ResidueId>>,
    residue_chains: Vec<ChainId>,
    residue_classes: Vec<ResidueClass>,
    /// The explicitly assigned class of each residue, which publication keeps.
    residue_explicit: Vec<Option<ResidueClass>>,
    residue_labels: Vec<u32>,
    labels: SharedRows<ResidueLabels>,
    residue_sites: SiteLists,
    site_residues: Vec<ResidueId>,
    site_atoms: Vec<InstanceAtomId>,
    site_metadata: Vec<u32>,
    metadata: SharedRows<AtomSiteMetadata>,
    /// The site of each placed atom, built while staging so that additions
    /// reject a second placement. A published hierarchy drops it; its
    /// topology's index answers instead.
    atom_lookup: Option<QualifiedMap<InstanceAtomId, AtomSiteId>>,
}

/// The atom sites of each residue, in residue order.
#[derive(Debug, Clone)]
enum SiteLists {
    /// One list per residue, while staging.
    Staged(Vec<Vec<AtomSiteId>>),
    /// Residue `r` owns `items[starts[r]..starts[r + 1]]`.
    Packed {
        starts: Vec<u32>,
        items: Vec<AtomSiteId>,
    },
}

impl Default for SiteLists {
    fn default() -> Self {
        Self::Staged(Vec::new())
    }
}

impl SiteLists {
    fn len(&self) -> usize {
        match self {
            Self::Staged(lists) => lists.len(),
            Self::Packed { starts, .. } => starts.len() - 1,
        }
    }

    fn get(&self, residue: usize) -> &[AtomSiteId] {
        match self {
            Self::Staged(lists) => &lists[residue],
            Self::Packed { starts, items } => {
                &items[starts[residue] as usize..starts[residue + 1] as usize]
            }
        }
    }

    fn staged(&mut self) -> &mut Vec<Vec<AtomSiteId>> {
        if let Self::Packed { .. } = self {
            let lists = (0..self.len()).map(|r| self.get(r).to_vec()).collect();
            *self = Self::Staged(lists);
        }
        match self {
            Self::Staged(lists) => lists,
            Self::Packed { .. } => unreachable!("just staged"),
        }
    }

    fn pack(&mut self) {
        if let Self::Staged(lists) = self {
            let mut starts = Vec::with_capacity(lists.len() + 1);
            let mut items = Vec::with_capacity(lists.iter().map(Vec::len).sum());
            starts.push(0);
            for list in lists.iter() {
                items.extend_from_slice(list);
                starts.push(items.len() as u32);
            }
            *self = Self::Packed { starts, items };
        }
    }
}

/// The columns of a hierarchy, for constructing one directly.
pub(crate) struct HierarchyColumns {
    pub(crate) chain_label_ids: Vec<String>,
    pub(crate) chain_author_ids: Vec<Option<String>>,
    pub(crate) chain_residues: Vec<Vec<ResidueId>>,
    pub(crate) residue_chains: Vec<ChainId>,
    pub(crate) residue_classes: Vec<ResidueClass>,
    pub(crate) residue_explicit: Vec<Option<ResidueClass>>,
    pub(crate) residue_labels: Vec<u32>,
    pub(crate) labels: SharedRows<ResidueLabels>,
    /// Residue `r`'s sites are `residue_site_items[residue_site_starts[r]..][..]`.
    pub(crate) residue_site_starts: Vec<u32>,
    pub(crate) residue_site_items: Vec<AtomSiteId>,
    pub(crate) site_residues: Vec<ResidueId>,
    pub(crate) site_atoms: Vec<InstanceAtomId>,
    pub(crate) site_metadata: Vec<u32>,
    pub(crate) metadata: SharedRows<AtomSiteMetadata>,
}

impl Hierarchy {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// A published hierarchy from consistent columns.
    pub(crate) fn from_columns(columns: HierarchyColumns) -> Self {
        let HierarchyColumns {
            chain_label_ids,
            chain_author_ids,
            chain_residues,
            residue_chains,
            residue_classes,
            residue_explicit,
            residue_labels,
            mut labels,
            residue_site_starts,
            residue_site_items,
            site_residues,
            site_atoms,
            mut site_metadata,
            mut metadata,
        } = columns;
        let mut residue_labels = residue_labels;
        labels.compact(&mut residue_labels);
        metadata.compact(&mut site_metadata);
        Self {
            chain_label_ids,
            chain_author_ids,
            chain_residues,
            residue_chains,
            residue_classes,
            residue_explicit,
            residue_labels,
            labels,
            residue_sites: SiteLists::Packed {
                starts: residue_site_starts,
                items: residue_site_items,
            },
            site_residues,
            site_atoms,
            site_metadata,
            metadata,
            atom_lookup: None,
        }
    }

    /// Drops staging state once the hierarchy is published.
    pub(crate) fn freeze(&mut self) {
        self.residue_sites.pack();
        self.atom_lookup = None;
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.chain_label_ids.is_empty()
            && self.residue_chains.is_empty()
            && self.site_atoms.is_empty()
    }

    /// The residue label rows, shared with snapshots that keep them.
    pub(crate) fn label_rows(&self) -> &SharedRows<ResidueLabels> {
        &self.labels
    }

    /// The site metadata rows, shared with snapshots that keep them.
    pub(crate) fn metadata_rows(&self) -> &SharedRows<AtomSiteMetadata> {
        &self.metadata
    }

    /// Each site's atom, by site ID.
    pub(crate) fn site_atom_ids(&self) -> &[InstanceAtomId] {
        &self.site_atoms
    }

    /// Each site's residue, by site ID.
    pub(crate) fn site_residue_ids(&self) -> &[ResidueId] {
        &self.site_residues
    }

    /// Each site's metadata key, by site ID.
    pub(crate) fn site_metadata_keys(&self) -> &[u32] {
        &self.site_metadata
    }

    /// Adds a chain or returns a structured capacity error without mutation.
    pub(crate) fn add_chain(
        &mut self,
        label_id: impl Into<String>,
        author_id: Option<String>,
    ) -> Result<ChainId, HierarchyError> {
        let label_id = label_id.into();
        let id = checked_hierarchy_id(self.chain_count(), Countable::Chain, ChainId::new)?;
        self.chain_label_ids.push(label_id);
        self.chain_author_ids.push(author_id);
        self.chain_residues.push(Vec::new());
        Ok(id)
    }

    /// Adds a residue to an existing chain transactionally.
    pub(crate) fn add_residue(
        &mut self,
        chain: ChainId,
        name: impl Into<String>,
        label_seq_id: Option<i32>,
        author_seq_id: Option<String>,
        insertion_code: Option<String>,
    ) -> Result<ResidueId, HierarchyError> {
        self.require_chain(chain)?;
        let id = checked_hierarchy_id(
            self.residue_chains.len(),
            Countable::Residue,
            ResidueId::new,
        )?;
        let name = Arc::<str>::from(name.into());
        let key = self.labels.push(ResidueLabels {
            name: Arc::clone(&name),
            label_comp_id: Some(name),
            author_comp_id: None,
            label_seq_id,
            author_seq_id: author_seq_id.map(Arc::from),
            insertion_code: insertion_code.map(Arc::from),
        });
        self.residue_chains.push(chain);
        self.residue_classes.push(ResidueClass::Other);
        self.residue_explicit.push(None);
        self.residue_labels.push(key);
        self.residue_sites.staged().push(Vec::new());
        self.chain_residues[chain.index()].push(id);
        Ok(id)
    }

    /// Adds one atom site after validating residue, atom placement, and capacity.
    pub(crate) fn add_atom_site(
        &mut self,
        residue: ResidueId,
        atom: InstanceAtomId,
        metadata: AtomSiteMetadata,
    ) -> Result<AtomSiteId, HierarchyError> {
        self.require_residue(residue)?;
        if self.lookup().contains_key(&atom) {
            return Err(HierarchyError::DuplicateAtomPlacement(atom));
        }
        let id = checked_hierarchy_id(self.site_atoms.len(), Countable::AtomSite, AtomSiteId::new)?;
        let key = self.metadata.push(metadata);
        self.push_site(id, residue, atom, key);
        Ok(id)
    }

    /// [`Self::add_atom_site`] with metadata shared with another snapshot.
    pub(crate) fn add_shared_atom_site(
        &mut self,
        residue: ResidueId,
        atom: InstanceAtomId,
        metadata: &SharedRow<AtomSiteMetadata>,
    ) -> Result<AtomSiteId, HierarchyError> {
        self.require_residue(residue)?;
        if self.lookup().contains_key(&atom) {
            return Err(HierarchyError::DuplicateAtomPlacement(atom));
        }
        let id = checked_hierarchy_id(self.site_atoms.len(), Countable::AtomSite, AtomSiteId::new)?;
        let key = self.metadata.key(metadata);
        self.push_site(id, residue, atom, key);
        Ok(id)
    }

    fn push_site(&mut self, id: AtomSiteId, residue: ResidueId, atom: InstanceAtomId, key: u32) {
        self.site_residues.push(residue);
        self.site_atoms.push(atom);
        self.site_metadata.push(key);
        self.residue_sites.staged()[residue.index()].push(id);
        self.lookup().insert(atom, id);
    }

    /// Builds the lookup from atoms to sites, so [`Self::atom_site_for_atom`]
    /// answers in constant time.
    pub(crate) fn index_atoms(&mut self) {
        self.lookup();
    }

    /// Takes the lookup from atoms to sites, building it first.
    pub(crate) fn take_atom_lookup(&mut self) -> QualifiedMap<InstanceAtomId, AtomSiteId> {
        self.lookup();
        self.atom_lookup.take().expect("just built")
    }

    /// The staging lookup from atoms to sites, built from the sites first.
    fn lookup(&mut self) -> &mut QualifiedMap<InstanceAtomId, AtomSiteId> {
        let sites = &self.site_atoms;
        self.atom_lookup.get_or_insert_with(|| {
            sites
                .iter()
                .enumerate()
                .map(|(index, &atom)| (atom, AtomSiteId::new(index as u32)))
                .collect()
        })
    }

    /// Reconstructs stored hierarchy from its node lists.
    ///
    /// Node IDs are list positions. Each residue must be listed by exactly one
    /// chain and each site by exactly one residue; those lists give every
    /// residue its chain and every site its residue, in stored child order.
    /// Atom placements must be unique. Atoms are checked against a topology
    /// only at publication.
    pub(crate) fn from_raw_parts(
        chains: Vec<RawChain>,
        residues: Vec<RawResidue>,
        sites: Vec<RawAtomSite>,
    ) -> Result<Self, RawHierarchyError> {
        let id = |index: usize, kind| {
            crate::molecule::checked_raw_id(index)
                .map_err(|_| RawHierarchyError::CapacityExceeded(CapacityExceeded::new(kind)))
        };
        id(chains.len().saturating_sub(1), Countable::Chain)?;
        id(residues.len().saturating_sub(1), Countable::Residue)?;
        id(sites.len().saturating_sub(1), Countable::AtomSite)?;
        let mut residue_chain = vec![None; residues.len()];
        for (index, chain) in chains.iter().enumerate() {
            let chain_id = ChainId::new(index as u32);
            for &residue in &chain.residues {
                let slot = residue_chain
                    .get_mut(residue.index())
                    .ok_or(RawHierarchyError::InvalidResidue(residue))?;
                if slot.replace(chain_id).is_some() {
                    return Err(RawHierarchyError::ResidueListedTwice(residue));
                }
            }
        }
        let mut site_residue = vec![None; sites.len()];
        for (index, residue) in residues.iter().enumerate() {
            let residue_id = ResidueId::new(index as u32);
            for &site in &residue.atom_sites {
                let slot = site_residue
                    .get_mut(site.index())
                    .ok_or(RawHierarchyError::InvalidAtomSite(site))?;
                if slot.replace(residue_id).is_some() {
                    return Err(RawHierarchyError::AtomSiteListedTwice(site));
                }
            }
        }
        let mut hierarchy = Self::new();
        for chain in chains {
            hierarchy.chain_label_ids.push(chain.label_id);
            hierarchy.chain_author_ids.push(chain.author_id);
            hierarchy.chain_residues.push(chain.residues);
        }
        let mut lists = Vec::with_capacity(residues.len());
        let mut labels = Vec::with_capacity(residues.len());
        for (index, (residue, chain)) in residues.into_iter().zip(residue_chain).enumerate() {
            let id = ResidueId::new(index as u32);
            hierarchy
                .residue_chains
                .push(chain.ok_or(RawHierarchyError::ResidueUnlisted(id))?);
            hierarchy.residue_classes.push(residue.class);
            hierarchy
                .residue_explicit
                .push(residue.explicit.then_some(residue.class));
            hierarchy.residue_labels.push(index as u32);
            labels.push(ResidueLabels {
                name: residue.name.into(),
                label_comp_id: residue.label_comp_id.map(Arc::from),
                author_comp_id: residue.author_comp_id.map(Arc::from),
                label_seq_id: residue.label_seq_id,
                author_seq_id: residue.author_seq_id.map(Arc::from),
                insertion_code: residue.insertion_code.map(Arc::from),
            });
            lists.push(residue.atom_sites);
        }
        hierarchy.labels = labels.into_iter().collect();
        hierarchy.residue_sites = SiteLists::Staged(lists);
        let mut lookup = Vec::with_capacity(sites.len());
        let mut metadata = Vec::with_capacity(sites.len());
        for (index, (site, residue)) in sites.into_iter().zip(site_residue).enumerate() {
            let id = AtomSiteId::new(index as u32);
            let residue = residue.ok_or(RawHierarchyError::AtomSiteUnlisted(id))?;
            lookup.push((site.atom, id));
            hierarchy.site_residues.push(residue);
            hierarchy.site_atoms.push(site.atom);
            hierarchy.site_metadata.push(index as u32);
            metadata.push(site.metadata);
        }
        hierarchy.metadata = metadata.into_iter().collect();
        // Sorted entries build the lookup in bulk rather than one insertion
        // at a time.
        lookup.sort_unstable();
        if let Some(pair) = lookup.windows(2).find(|pair| pair[0].0 == pair[1].0) {
            return Err(RawHierarchyError::DuplicateAtomPlacement(pair[0].0));
        }
        hierarchy.atom_lookup = Some(lookup.into_iter().collect());
        Ok(hierarchy)
    }

    /// Mutable residue list of one chain, for corrupting hierarchies in tests.
    #[cfg(test)]
    pub(crate) fn chain_residues_mut(&mut self, chain: ChainId) -> &mut Vec<ResidueId> {
        &mut self.chain_residues[chain.index()]
    }

    /// Mutable chain of one residue, for corrupting hierarchies in tests.
    #[cfg(test)]
    pub(crate) fn residue_chain_mut(&mut self, residue: ResidueId) -> &mut ChainId {
        &mut self.residue_chains[residue.index()]
    }

    /// Mutable site list of one residue, for corrupting hierarchies in tests.
    #[cfg(test)]
    pub(crate) fn residue_sites_mut(&mut self, residue: ResidueId) -> &mut Vec<AtomSiteId> {
        &mut self.residue_sites.staged()[residue.index()]
    }

    pub(crate) fn chain_count(&self) -> usize {
        self.chain_label_ids.len()
    }

    pub(crate) fn residue_count(&self) -> usize {
        self.residue_chains.len()
    }

    pub(crate) fn atom_site_count(&self) -> usize {
        self.site_atoms.len()
    }

    /// Every chain ID in hierarchy order.
    pub(crate) fn chain_ids(&self) -> impl ExactSizeIterator<Item = ChainId> {
        (0..self.chain_count() as u32).map(ChainId::new)
    }

    /// Every residue ID in hierarchy order.
    pub(crate) fn residue_ids(&self) -> impl ExactSizeIterator<Item = ResidueId> {
        (0..self.residue_count() as u32).map(ResidueId::new)
    }

    /// Every atom-site ID in hierarchy order.
    pub(crate) fn atom_site_ids(&self) -> impl ExactSizeIterator<Item = AtomSiteId> {
        (0..self.atom_site_count() as u32).map(AtomSiteId::new)
    }

    pub(crate) fn has_chain(&self, id: ChainId) -> bool {
        id.index() < self.chain_count()
    }

    pub(crate) fn has_residue(&self, id: ResidueId) -> bool {
        id.index() < self.residue_count()
    }

    pub(crate) fn has_atom_site(&self, id: AtomSiteId) -> bool {
        id.index() < self.atom_site_count()
    }

    /// [`Self::has_chain`] for staging code that reports a missing chain.
    pub(crate) fn require_chain(&self, id: ChainId) -> Result<(), HierarchyError> {
        self.has_chain(id)
            .then_some(())
            .ok_or(HierarchyError::InvalidChainId(id))
    }

    /// [`Self::has_residue`] for staging code that reports a missing residue.
    pub(crate) fn require_residue(&self, id: ResidueId) -> Result<(), HierarchyError> {
        self.has_residue(id)
            .then_some(())
            .ok_or(HierarchyError::InvalidResidueId(id))
    }

    /// [`Self::has_atom_site`] for staging code that reports a missing site.
    pub(crate) fn require_atom_site(&self, id: AtomSiteId) -> Result<(), HierarchyError> {
        self.has_atom_site(id)
            .then_some(())
            .ok_or(HierarchyError::InvalidAtomSiteId(id))
    }

    // Column reads take an ID of this hierarchy and panic on any other.

    pub(crate) fn chain_label_id(&self, chain: ChainId) -> &str {
        &self.chain_label_ids[chain.index()]
    }

    pub(crate) fn chain_author_id(&self, chain: ChainId) -> Option<&str> {
        self.chain_author_ids[chain.index()].as_deref()
    }

    /// The chain's residues in chain order.
    pub(crate) fn chain_residues(&self, chain: ChainId) -> &[ResidueId] {
        &self.chain_residues[chain.index()]
    }

    pub(crate) fn residue_chain(&self, residue: ResidueId) -> ChainId {
        self.residue_chains[residue.index()]
    }

    pub(crate) fn residue_class(&self, residue: ResidueId) -> ResidueClass {
        self.residue_classes[residue.index()]
    }

    /// The explicitly assigned class, which publication keeps.
    pub(crate) fn explicit_residue_class(&self, residue: ResidueId) -> Option<ResidueClass> {
        self.residue_explicit[residue.index()]
    }

    pub(crate) fn residue_labels(&self, residue: ResidueId) -> &ResidueLabels {
        self.labels.get(self.residue_label_key(residue))
    }

    /// The key of the residue's labels in [`Self::label_rows`].
    pub(crate) fn residue_label_key(&self, residue: ResidueId) -> u32 {
        self.residue_labels[residue.index()]
    }

    /// The residue's labels as a row shared with this hierarchy.
    pub(crate) fn shared_residue_labels(&self, residue: ResidueId) -> SharedRow<ResidueLabels> {
        self.labels.row(self.residue_label_key(residue))
    }

    /// The residue's atom sites in residue order.
    pub(crate) fn residue_sites(&self, residue: ResidueId) -> &[AtomSiteId] {
        self.residue_sites.get(residue.index())
    }

    pub(crate) fn site_residue(&self, site: AtomSiteId) -> ResidueId {
        self.site_residues[site.index()]
    }

    pub(crate) fn site_atom(&self, site: AtomSiteId) -> InstanceAtomId {
        self.site_atoms[site.index()]
    }

    pub(crate) fn site_metadata(&self, site: AtomSiteId) -> &AtomSiteMetadata {
        self.metadata.get(self.site_metadata[site.index()])
    }

    /// The site's metadata as a row shared with this hierarchy.
    pub(crate) fn shared_site_metadata(&self, site: AtomSiteId) -> SharedRow<AtomSiteMetadata> {
        self.metadata.row(self.site_metadata[site.index()])
    }

    /// The site placing `atom`, from the staging lookup, or by scanning the
    /// sites of a published hierarchy, whose topology answers in constant
    /// time instead.
    pub(crate) fn atom_site_for_atom(&self, atom: InstanceAtomId) -> Option<AtomSiteId> {
        match &self.atom_lookup {
            Some(lookup) => lookup.get(&atom).copied(),
            None => {
                let index = self.site_atoms.iter().position(|&site| site == atom)?;
                Some(AtomSiteId::new(index as u32))
            }
        }
    }

    /// Every placed atom with its site, in qualified atom order.
    pub(crate) fn atom_lookup_entries(&self) -> Vec<(InstanceAtomId, AtomSiteId)> {
        let mut entries = self
            .site_atoms
            .iter()
            .enumerate()
            .map(|(index, &atom)| (atom, AtomSiteId::new(index as u32)))
            .collect::<Vec<_>>();
        entries.sort_unstable();
        entries
    }

    /// Restores distinct label and author component identifiers on one residue.
    pub(crate) fn set_residue_component_ids(
        &mut self,
        residue: ResidueId,
        label_comp_id: Option<String>,
        author_comp_id: Option<String>,
    ) -> Result<(), HierarchyError> {
        self.require_residue(residue)?;
        let labels = ResidueLabels {
            label_comp_id: label_comp_id.map(Arc::from),
            author_comp_id: author_comp_id.map(Arc::from),
            ..self.residue_labels(residue).clone()
        };
        self.residue_labels[residue.index()] = self.labels.push(labels);
        Ok(())
    }

    /// Records or clears an explicit class on one residue node.
    pub(crate) fn set_explicit_residue_class(
        &mut self,
        residue: ResidueId,
        class: Option<ResidueClass>,
    ) -> Result<(), HierarchyError> {
        self.require_residue(residue)?;
        self.residue_explicit[residue.index()] = class;
        Ok(())
    }

    /// Records one residue's published class.
    pub(crate) fn set_residue_class(&mut self, residue: ResidueId, class: ResidueClass) {
        self.residue_classes[residue.index()] = class;
    }

    /// Equal chains, residues, and atom sites, ignoring explicit-class
    /// intent, which is editing policy rather than layout.
    pub(crate) fn same_layout(&self, other: &Self) -> bool {
        self.equal(other, false)
    }

    fn equal(&self, other: &Self, explicit: bool) -> bool {
        self.chain_label_ids == other.chain_label_ids
            && self.chain_author_ids == other.chain_author_ids
            && self.chain_residues == other.chain_residues
            && self.residue_chains == other.residue_chains
            && self.residue_classes == other.residue_classes
            && (!explicit || self.residue_explicit == other.residue_explicit)
            && self.site_residues == other.site_residues
            && self.site_atoms == other.site_atoms
            && self.residue_ids().all(|residue| {
                self.residue_labels(residue) == other.residue_labels(residue)
                    && self.residue_sites(residue) == other.residue_sites(residue)
            })
            && self
                .atom_site_ids()
                .all(|site| self.site_metadata(site) == other.site_metadata(site))
    }
}

impl PartialEq for Hierarchy {
    fn eq(&self, other: &Self) -> bool {
        self.equal(other, true)
    }
}

/// Staged chains, residues, and atom sites of a topology under construction.
///
/// [`TopologyBuilder::hierarchy_mut`](super::TopologyBuilder::hierarchy_mut)
/// and [`StructureBuilder::hierarchy_mut`](crate::structure::StructureBuilder::hierarchy_mut)
/// return the one staged hierarchy. Additions validate their own references
/// and capacity; atom-site atoms are checked against the final instance atoms
/// when the topology is built. A published topology is read through its
/// views ([`Topology::chains`](super::Topology::chains), ...).
///
/// A residue's explicit class is staged on the residue itself, so replacing
/// the staged hierarchy with a new one also discards its explicit classes.
///
/// ```
/// use kekule::{molecule::AtomId, smiles};
/// use kekule::topology::{AtomSiteMetadata, InstanceAtomId, ResidueClass, TopologyBuilder};
/// let mut builder = TopologyBuilder::new();
/// let water = builder.add_molecule(smiles::to_molecules("O")?.remove(0))?;
/// let hierarchy = builder.hierarchy_mut();
/// let chain = hierarchy.add_chain("W", None)?;
/// let residue = hierarchy.add_residue(chain, "HOH", None, Some("1".into()), None)?;
/// let oxygen = InstanceAtomId::new(water, AtomId::new(0));
/// hierarchy.add_atom_site(residue, oxygen, AtomSiteMetadata::default())?;
/// hierarchy.set_residue_class(residue, ResidueClass::Other)?;
/// let topology = builder.build()?;
/// let residue = topology.residue(residue).unwrap();
/// assert_eq!(residue.class(), ResidueClass::Other);
/// assert!(residue.class_assignment().is_explicit());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct HierarchyBuilder {
    hierarchy: Hierarchy,
}

impl HierarchyBuilder {
    /// An empty staged hierarchy.
    pub fn new() -> Self {
        Self::default()
    }

    pub(crate) const fn from_storage(hierarchy: Hierarchy) -> Self {
        Self { hierarchy }
    }

    pub(crate) const fn storage(&self) -> &Hierarchy {
        &self.hierarchy
    }

    pub(crate) fn storage_mut(&mut self) -> &mut Hierarchy {
        &mut self.hierarchy
    }

    pub(crate) fn into_storage(self) -> Hierarchy {
        self.hierarchy
    }

    /// Whether no chain, residue, or atom site is staged.
    pub fn is_empty(&self) -> bool {
        self.hierarchy.is_empty()
    }

    /// The number of staged chains.
    pub fn chain_count(&self) -> usize {
        self.hierarchy.chain_count()
    }

    /// The number of staged residues.
    pub fn residue_count(&self) -> usize {
        self.hierarchy.residue_count()
    }

    /// The number of staged atom sites.
    pub fn atom_site_count(&self) -> usize {
        self.hierarchy.atom_site_count()
    }

    /// Adds a chain with a label identifier and an optional author identifier.
    pub fn add_chain(
        &mut self,
        label_id: impl Into<String>,
        author_id: Option<String>,
    ) -> Result<ChainId, HierarchyError> {
        self.hierarchy.add_chain(label_id, author_id)
    }

    /// Adds a residue to an existing chain. Its label component identifier is
    /// its name until [`Self::set_residue_component_ids`] changes it.
    pub fn add_residue(
        &mut self,
        chain: ChainId,
        name: impl Into<String>,
        label_seq_id: Option<i32>,
        author_seq_id: Option<String>,
        insertion_code: Option<String>,
    ) -> Result<ResidueId, HierarchyError> {
        self.hierarchy
            .add_residue(chain, name, label_seq_id, author_seq_id, insertion_code)
    }

    /// Adds one atom site to an existing residue; an atom has at most one
    /// site.
    pub fn add_atom_site(
        &mut self,
        residue: ResidueId,
        atom: InstanceAtomId,
        metadata: AtomSiteMetadata,
    ) -> Result<AtomSiteId, HierarchyError> {
        self.hierarchy.add_atom_site(residue, atom, metadata)
    }

    /// Sets distinct label and author component identifiers on one residue.
    pub fn set_residue_component_ids(
        &mut self,
        residue: ResidueId,
        label_comp_id: Option<String>,
        author_comp_id: Option<String>,
    ) -> Result<(), HierarchyError> {
        self.hierarchy
            .set_residue_component_ids(residue, label_comp_id, author_comp_id)
    }

    /// Assigns one residue's class explicitly. Publication keeps it instead
    /// of inferring one, even when the residue's evidence changes.
    pub fn set_residue_class(
        &mut self,
        residue: ResidueId,
        class: ResidueClass,
    ) -> Result<(), HierarchyError> {
        self.hierarchy
            .set_explicit_residue_class(residue, Some(class))
    }

    /// Clears one residue's explicit class, so publication infers it afresh.
    pub fn clear_residue_class(&mut self, residue: ResidueId) -> Result<(), HierarchyError> {
        self.hierarchy.set_explicit_residue_class(residue, None)
    }

    /// One residue's explicit class, or `None` when publication will infer
    /// it or `residue` is not staged.
    pub fn explicit_residue_class(&self, residue: ResidueId) -> Option<ResidueClass> {
        self.hierarchy
            .has_residue(residue)
            .then(|| self.hierarchy.explicit_residue_class(residue))?
    }
}

/// One stored chain for [`Hierarchy::from_raw_parts`].
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RawChain {
    pub(crate) label_id: String,
    pub(crate) author_id: Option<String>,
    pub(crate) residues: Vec<ResidueId>,
}

/// One stored residue for [`Hierarchy::from_raw_parts`]; `class` is its
/// published class and `explicit` whether it was assigned explicitly.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RawResidue {
    pub(crate) name: String,
    pub(crate) label_comp_id: Option<String>,
    pub(crate) author_comp_id: Option<String>,
    pub(crate) label_seq_id: Option<i32>,
    pub(crate) author_seq_id: Option<String>,
    pub(crate) insertion_code: Option<String>,
    pub(crate) class: ResidueClass,
    pub(crate) explicit: bool,
    pub(crate) atom_sites: Vec<AtomSiteId>,
}

/// One stored atom site for [`Hierarchy::from_raw_parts`].
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RawAtomSite {
    pub(crate) atom: InstanceAtomId,
    pub(crate) metadata: AtomSiteMetadata,
}

/// Why stored hierarchy node lists do not form a hierarchy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RawHierarchyError {
    CapacityExceeded(CapacityExceeded),
    InvalidResidue(ResidueId),
    InvalidAtomSite(AtomSiteId),
    ResidueListedTwice(ResidueId),
    AtomSiteListedTwice(AtomSiteId),
    ResidueUnlisted(ResidueId),
    AtomSiteUnlisted(AtomSiteId),
    DuplicateAtomPlacement(InstanceAtomId),
}

impl fmt::Display for RawHierarchyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use crate::diagnostics::Capitalized;
        match self {
            Self::CapacityExceeded(error) => error.fmt(f),
            Self::InvalidResidue(id) => write!(f, "A chain lists {id}, which does not exist."),
            Self::InvalidAtomSite(id) => {
                write!(f, "A residue lists {id}, which does not exist.")
            }
            Self::ResidueListedTwice(id) => {
                write!(f, "{} is listed by more than one chain.", Capitalized(id))
            }
            Self::AtomSiteListedTwice(id) => {
                write!(f, "{} is listed by more than one residue.", Capitalized(id))
            }
            Self::ResidueUnlisted(id) => write!(f, "{} belongs to no chain.", Capitalized(id)),
            Self::AtomSiteUnlisted(id) => write!(f, "{} belongs to no residue.", Capitalized(id)),
            Self::DuplicateAtomPlacement(atom) => {
                write!(f, "{} has more than one atom site.", Capitalized(atom))
            }
        }
    }
}

impl std::error::Error for RawHierarchyError {}

/// The source labels of one residue, shared between snapshots.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ResidueLabels {
    pub(crate) name: Arc<str>,
    pub(crate) label_comp_id: Option<Arc<str>>,
    pub(crate) author_comp_id: Option<Arc<str>>,
    pub(crate) label_seq_id: Option<i32>,
    pub(crate) author_seq_id: Option<Arc<str>>,
    pub(crate) insertion_code: Option<Arc<str>>,
}

impl ResidueLabels {
    pub(crate) fn name(&self) -> &str {
        &self.name
    }

    pub(crate) fn label_comp_id(&self) -> Option<&str> {
        self.label_comp_id.as_deref()
    }

    pub(crate) fn author_comp_id(&self) -> Option<&str> {
        self.author_comp_id.as_deref()
    }

    pub(crate) fn label_seq_id(&self) -> Option<i32> {
        self.label_seq_id
    }

    pub(crate) fn author_seq_id(&self) -> Option<&str> {
        self.author_seq_id.as_deref()
    }

    pub(crate) fn insertion_code(&self) -> Option<&str> {
        self.insertion_code.as_deref()
    }
}

/// Optional source labels of one atom site, named after mmCIF `atom_site` items.
///
/// Every field defaults to `None` and is set with the `with_*` methods.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AtomSiteMetadata {
    pub(crate) type_symbol: Option<String>,
    pub(crate) label_asym_id: Option<String>,
    pub(crate) auth_asym_id: Option<String>,
    pub(crate) label_atom_id: Option<String>,
    pub(crate) auth_atom_id: Option<String>,
}

impl AtomSiteMetadata {
    /// The element type symbol (`type_symbol`), or `None` when absent.
    pub fn type_symbol(&self) -> Option<&str> {
        self.type_symbol.as_deref()
    }

    /// Sets [`Self::type_symbol`].
    #[must_use]
    pub fn with_type_symbol(mut self, type_symbol: Option<String>) -> Self {
        self.type_symbol = type_symbol;
        self
    }

    /// The label chain identifier (`label_asym_id`), or `None` when absent.
    pub fn label_asym_id(&self) -> Option<&str> {
        self.label_asym_id.as_deref()
    }

    /// Sets [`Self::label_asym_id`].
    #[must_use]
    pub fn with_label_asym_id(mut self, label_asym_id: Option<String>) -> Self {
        self.label_asym_id = label_asym_id;
        self
    }

    /// The author chain identifier (`auth_asym_id`), or `None` when absent.
    pub fn auth_asym_id(&self) -> Option<&str> {
        self.auth_asym_id.as_deref()
    }

    /// Sets [`Self::auth_asym_id`].
    #[must_use]
    pub fn with_auth_asym_id(mut self, auth_asym_id: Option<String>) -> Self {
        self.auth_asym_id = auth_asym_id;
        self
    }

    /// The label atom name (`label_atom_id`), or `None` when absent.
    pub fn label_atom_id(&self) -> Option<&str> {
        self.label_atom_id.as_deref()
    }

    /// Sets [`Self::label_atom_id`].
    #[must_use]
    pub fn with_label_atom_id(mut self, label_atom_id: Option<String>) -> Self {
        self.label_atom_id = label_atom_id;
        self
    }

    /// The author atom name (`auth_atom_id`), or `None` when absent.
    pub fn auth_atom_id(&self) -> Option<&str> {
        self.auth_atom_id.as_deref()
    }

    /// Sets [`Self::auth_atom_id`].
    #[must_use]
    pub fn with_auth_atom_id(mut self, auth_atom_id: Option<String>) -> Self {
        self.auth_atom_id = auth_atom_id;
        self
    }
}

/// Fixed-width identifier spaces of a topology hierarchy.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HierarchyIdKind {
    /// Chain identifiers.
    Chain,
    /// Residue identifiers.
    Residue,
    /// Atom-site identifiers.
    AtomSite,
}

impl fmt::Display for HierarchyIdKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Chain => "chain",
            Self::Residue => "residue",
            Self::AtomSite => "atom site",
        })
    }
}

/// Why a hierarchy operation or lookup failed.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HierarchyError {
    /// The chain does not exist.
    InvalidChainId(ChainId),
    /// The residue does not exist.
    InvalidResidueId(ResidueId),
    /// The atom site does not exist.
    InvalidAtomSiteId(AtomSiteId),
    /// The atom already has an atom site.
    DuplicateAtomPlacement(InstanceAtomId),
    /// A new hierarchy node cannot be given a fixed-width ID.
    CapacityExceeded(CapacityExceeded),
    /// A lookup by name or label found no match.
    NotFound {
        /// The kind of node looked up.
        kind: HierarchyIdKind,
        /// The name or label looked up.
        query: String,
    },
    /// A lookup by name or label found more than one match.
    Ambiguous {
        /// The kind of node looked up.
        kind: HierarchyIdKind,
        /// The name or label looked up.
        query: String,
        /// The number of matches.
        matches: usize,
    },
}

impl fmt::Display for HierarchyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use crate::diagnostics::Capitalized;
        match self {
            Self::InvalidChainId(id) => write!(f, "{} is not in the hierarchy.", Capitalized(id)),
            Self::InvalidResidueId(id) => {
                write!(f, "{} is not in the hierarchy.", Capitalized(id))
            }
            Self::InvalidAtomSiteId(id) => {
                write!(f, "{} is not in the hierarchy.", Capitalized(id))
            }
            Self::DuplicateAtomPlacement(id) => {
                write!(f, "{} already has an atom site.", Capitalized(id))
            }
            Self::CapacityExceeded(error) => error.fmt(f),
            Self::NotFound { kind, query } => write!(f, "No {kind} matches {query}."),
            Self::Ambiguous {
                kind,
                query,
                matches,
            } => write!(f, "{matches} {kind}s match {query}."),
        }
    }
}

impl std::error::Error for HierarchyError {}

fn checked_hierarchy_id<T>(
    length: usize,
    kind: Countable,
    construct: impl FnOnce(u32) -> T,
) -> Result<T, HierarchyError> {
    crate::molecule::checked_raw_id(length)
        .map(construct)
        .map_err(|_| HierarchyError::CapacityExceeded(CapacityExceeded::new(kind)))
}
