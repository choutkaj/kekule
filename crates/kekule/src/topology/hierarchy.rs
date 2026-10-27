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
    chains: Vec<Chain>,
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
    pub(crate) chains: Vec<Chain>,
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
            chains,
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
            chains,
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
        self.chains.is_empty() && self.residue_chains.is_empty() && self.site_atoms.is_empty()
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
    ) -> std::result::Result<ChainId, HierarchyError> {
        let label_id = label_id.into();
        let id = checked_hierarchy_id(self.chains.len(), HierarchyIdKind::Chain, ChainId::new)?;
        self.chains.push(Chain {
            id,
            label_id,
            author_id,
            residues: Vec::new(),
        });
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
    ) -> std::result::Result<ResidueId, HierarchyError> {
        self.require_chain(chain)?;
        let id = checked_hierarchy_id(
            self.residue_chains.len(),
            HierarchyIdKind::Residue,
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
        self.chains[chain.index()].residues.push(id);
        Ok(id)
    }

    /// Adds one atom site after validating residue, atom placement, and capacity.
    pub(crate) fn add_atom_site(
        &mut self,
        residue: ResidueId,
        atom: InstanceAtomId,
        metadata: AtomSiteMetadata,
    ) -> std::result::Result<AtomSiteId, HierarchyError> {
        self.require_residue(residue)?;
        if self.lookup().contains_key(&atom) {
            return Err(HierarchyError::DuplicateAtomPlacement(atom));
        }
        let id = checked_hierarchy_id(
            self.site_atoms.len(),
            HierarchyIdKind::AtomSite,
            AtomSiteId::new,
        )?;
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
    ) -> std::result::Result<AtomSiteId, HierarchyError> {
        self.require_residue(residue)?;
        if self.lookup().contains_key(&atom) {
            return Err(HierarchyError::DuplicateAtomPlacement(atom));
        }
        let id = checked_hierarchy_id(
            self.site_atoms.len(),
            HierarchyIdKind::AtomSite,
            AtomSiteId::new,
        )?;
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
    ) -> std::result::Result<Self, RawHierarchyError> {
        let id = |index: usize, kind| {
            crate::molecule::checked_raw_id(index)
                .map_err(|_| RawHierarchyError::IdentifierCapacityExceeded(kind))
        };
        id(chains.len().saturating_sub(1), HierarchyIdKind::Chain)?;
        id(residues.len().saturating_sub(1), HierarchyIdKind::Residue)?;
        id(sites.len().saturating_sub(1), HierarchyIdKind::AtomSite)?;
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
        hierarchy.chains = chains
            .into_iter()
            .enumerate()
            .map(|(index, chain)| Chain {
                id: ChainId::new(index as u32),
                label_id: chain.label_id,
                author_id: chain.author_id,
                residues: chain.residues,
            })
            .collect();
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
        &mut self.chains[chain.index()].residues
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

    /// One chain, or `None` when `id` names no chain of this hierarchy.
    pub(crate) fn chain(&self, id: ChainId) -> Option<&Chain> {
        self.chains.get(id.index())
    }

    /// One residue, or `None` when `id` names no residue of this hierarchy.
    pub(crate) fn residue(&self, id: ResidueId) -> Option<Residue<'_>> {
        (id.index() < self.residue_chains.len()).then_some(Residue {
            hierarchy: self,
            id,
        })
    }

    /// One atom site, or `None` when `id` names no site of this hierarchy.
    pub(crate) fn atom_site(&self, id: AtomSiteId) -> Option<AtomSite<'_>> {
        (id.index() < self.site_atoms.len()).then_some(AtomSite {
            hierarchy: self,
            id,
        })
    }

    /// [`Self::chain`] for staging code that reports a missing chain.
    pub(crate) fn require_chain(&self, id: ChainId) -> std::result::Result<&Chain, HierarchyError> {
        self.chain(id).ok_or(HierarchyError::InvalidChainId(id))
    }

    /// [`Self::residue`] for staging code that reports a missing residue.
    pub(crate) fn require_residue(
        &self,
        id: ResidueId,
    ) -> std::result::Result<Residue<'_>, HierarchyError> {
        self.residue(id).ok_or(HierarchyError::InvalidResidueId(id))
    }

    /// [`Self::atom_site`] for staging code that reports a missing site.
    pub(crate) fn require_atom_site(
        &self,
        id: AtomSiteId,
    ) -> std::result::Result<AtomSite<'_>, HierarchyError> {
        self.atom_site(id)
            .ok_or(HierarchyError::InvalidAtomSiteId(id))
    }

    /// The site placing `atom`, from the staging lookup, or by scanning the
    /// sites of a published hierarchy, whose topology answers in constant
    /// time instead.
    pub(crate) fn atom_site_for_atom(&self, atom: InstanceAtomId) -> Option<AtomSite<'_>> {
        let id = match &self.atom_lookup {
            Some(lookup) => *lookup.get(&atom)?,
            None => AtomSiteId::new(self.site_atoms.iter().position(|&site| site == atom)? as u32),
        };
        self.atom_site(id)
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
    ) -> std::result::Result<(), HierarchyError> {
        let labels = ResidueLabels {
            label_comp_id: label_comp_id.map(Arc::from),
            author_comp_id: author_comp_id.map(Arc::from),
            ..self.require_residue(residue)?.labels().clone()
        };
        self.residue_labels[residue.index()] = self.labels.push(labels);
        Ok(())
    }

    /// Records or clears an explicit class on one residue node.
    pub(crate) fn set_explicit_residue_class(
        &mut self,
        residue: ResidueId,
        class: Option<ResidueClass>,
    ) -> std::result::Result<(), HierarchyError> {
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
        self.chains == other.chains
            && self.residue_chains == other.residue_chains
            && self.residue_classes == other.residue_classes
            && (!explicit || self.residue_explicit == other.residue_explicit)
            && self.site_residues == other.site_residues
            && self.site_atoms == other.site_atoms
            && self
                .residues()
                .zip(other.residues())
                .all(|((_, a), (_, b))| {
                    a.labels() == b.labels() && a.atom_sites() == b.atom_sites()
                })
            && self
                .atom_sites()
                .zip(other.atom_sites())
                .all(|((_, a), (_, b))| a.metadata() == b.metadata())
    }

    pub(crate) fn chains(&self) -> impl ExactSizeIterator<Item = (ChainId, &Chain)> {
        self.chains.iter().map(|chain| (chain.id, chain))
    }

    pub(crate) fn residues(&self) -> impl ExactSizeIterator<Item = (ResidueId, Residue<'_>)> {
        (0..self.residue_chains.len()).map(move |index| {
            let id = ResidueId::new(index as u32);
            (
                id,
                Residue {
                    hierarchy: self,
                    id,
                },
            )
        })
    }

    pub(crate) fn atom_sites(&self) -> impl ExactSizeIterator<Item = (AtomSiteId, AtomSite<'_>)> {
        (0..self.site_atoms.len()).map(move |index| {
            let id = AtomSiteId::new(index as u32);
            (
                id,
                AtomSite {
                    hierarchy: self,
                    id,
                },
            )
        })
    }

    /// The number of residues.
    pub(crate) fn residue_count(&self) -> usize {
        self.residue_chains.len()
    }

    /// The number of atom sites.
    pub(crate) fn atom_site_count(&self) -> usize {
        self.site_atoms.len()
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
/// and [`ModelBuilder::hierarchy_mut`](crate::structure::ModelBuilder::hierarchy_mut)
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
        self.hierarchy.chains.len()
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
        self.hierarchy.residue(residue)?.explicit_class()
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
    IdentifierCapacityExceeded(HierarchyIdKind),
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
        match self {
            Self::IdentifierCapacityExceeded(kind) => {
                write!(f, "{kind} identifier capacity exceeded")
            }
            Self::InvalidResidue(id) => write!(f, "a chain lists missing {id}"),
            Self::InvalidAtomSite(id) => write!(f, "a residue lists missing {id}"),
            Self::ResidueListedTwice(id) => write!(f, "{id} is listed by more than one chain"),
            Self::AtomSiteListedTwice(id) => {
                write!(f, "{id} is listed by more than one residue")
            }
            Self::ResidueUnlisted(id) => write!(f, "{id} belongs to no chain"),
            Self::AtomSiteUnlisted(id) => write!(f, "{id} belongs to no residue"),
            Self::DuplicateAtomPlacement(atom) => write!(f, "{atom} has more than one atom site"),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Chain {
    pub(crate) id: ChainId,
    pub(crate) label_id: String,
    pub(crate) author_id: Option<String>,
    pub(crate) residues: Vec<ResidueId>,
}

impl Chain {
    pub(crate) fn label_id(&self) -> &str {
        &self.label_id
    }

    pub(crate) fn author_id(&self) -> Option<&str> {
        self.author_id.as_deref()
    }

    pub(crate) fn residues(&self) -> &[ResidueId] {
        &self.residues
    }
}

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

/// One hierarchy residue: its chain, class, source identifiers, and atom sites.
#[derive(Clone, Copy)]
pub(crate) struct Residue<'a> {
    hierarchy: &'a Hierarchy,
    id: ResidueId,
}

impl fmt::Debug for Residue<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Residue")
            .field("id", &self.id)
            .field("chain", &self.chain())
            .field("class", &self.class())
            .field("labels", self.labels())
            .field("atom_sites", &self.atom_sites())
            .finish()
    }
}

impl<'a> Residue<'a> {
    pub(crate) const fn id(self) -> ResidueId {
        self.id
    }

    pub(crate) fn chain(self) -> ChainId {
        self.hierarchy.residue_chains[self.id.index()]
    }

    pub(crate) fn class(self) -> ResidueClass {
        self.hierarchy.residue_classes[self.id.index()]
    }

    /// The explicitly assigned class, which publication keeps.
    pub(crate) fn explicit_class(self) -> Option<ResidueClass> {
        self.hierarchy.residue_explicit[self.id.index()]
    }

    /// The residue's source labels.
    pub(crate) fn labels(self) -> &'a ResidueLabels {
        self.hierarchy.labels.get(self.label_key())
    }

    /// The residue's labels as a row shared with this hierarchy.
    pub(crate) fn shared_labels(self) -> SharedRow<ResidueLabels> {
        self.hierarchy.labels.row(self.label_key())
    }

    /// The key of the residue's labels in the hierarchy's label rows.
    pub(crate) fn label_key(self) -> u32 {
        self.hierarchy.residue_labels[self.id.index()]
    }

    pub(crate) fn name(self) -> &'a str {
        &self.labels().name
    }

    pub(crate) fn label_comp_id(self) -> Option<&'a str> {
        self.labels().label_comp_id.as_deref()
    }

    pub(crate) fn author_comp_id(self) -> Option<&'a str> {
        self.labels().author_comp_id.as_deref()
    }

    pub(crate) fn label_seq_id(self) -> Option<i32> {
        self.labels().label_seq_id
    }

    pub(crate) fn author_seq_id(self) -> Option<&'a str> {
        self.labels().author_seq_id.as_deref()
    }

    pub(crate) fn insertion_code(self) -> Option<&'a str> {
        self.labels().insertion_code.as_deref()
    }

    pub(crate) fn atom_sites(self) -> &'a [AtomSiteId] {
        self.hierarchy.residue_sites.get(self.id.index())
    }
}

/// One hierarchy atom site: the placement of one topology atom in a residue.
#[derive(Clone, Copy)]
pub(crate) struct AtomSite<'a> {
    hierarchy: &'a Hierarchy,
    id: AtomSiteId,
}

impl fmt::Debug for AtomSite<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("AtomSite")
            .field("id", &self.id)
            .field("residue", &self.residue())
            .field("atom", &self.atom())
            .field("metadata", self.metadata())
            .finish()
    }
}

impl<'a> AtomSite<'a> {
    pub(crate) const fn id(self) -> AtomSiteId {
        self.id
    }

    pub(crate) fn residue(self) -> ResidueId {
        self.hierarchy.site_residues[self.id.index()]
    }

    pub(crate) fn atom(self) -> InstanceAtomId {
        self.hierarchy.site_atoms[self.id.index()]
    }

    pub(crate) fn metadata(self) -> &'a AtomSiteMetadata {
        self.hierarchy.metadata.get(self.metadata_key())
    }

    /// The key of the site's metadata in the hierarchy's metadata rows.
    pub(crate) fn metadata_key(self) -> u32 {
        self.hierarchy.site_metadata[self.id.index()]
    }

    /// The site's metadata as a row shared with this hierarchy.
    pub(crate) fn shared_metadata(self) -> SharedRow<AtomSiteMetadata> {
        self.hierarchy.metadata.row(self.metadata_key())
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
            Self::AtomSite => "atom-site",
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
    /// A new hierarchy node cannot be represented by the fixed-width ID for `kind`.
    IdentifierCapacityExceeded(HierarchyIdKind),
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
        match self {
            Self::InvalidChainId(id) => write!(f, "invalid chain id: {}", id.raw()),
            Self::InvalidResidueId(id) => write!(f, "invalid residue id: {}", id.raw()),
            Self::InvalidAtomSiteId(id) => write!(f, "invalid atom-site id: {}", id.raw()),
            Self::DuplicateAtomPlacement(id) => write!(f, "duplicate hierarchy placement for {id}"),
            Self::IdentifierCapacityExceeded(kind) => {
                write!(f, "{kind} identifier capacity exceeded")
            }
            Self::NotFound { kind, query } => write!(f, "no {kind} matches {query}"),
            Self::Ambiguous {
                kind,
                query,
                matches,
            } => write!(f, "{matches} {kind} matches for {query}"),
        }
    }
}

impl std::error::Error for HierarchyError {}

fn checked_hierarchy_id<T>(
    length: usize,
    kind: HierarchyIdKind,
    construct: impl FnOnce(u32) -> T,
) -> std::result::Result<T, HierarchyError> {
    crate::molecule::checked_raw_id(length)
        .map(construct)
        .map_err(|_| HierarchyError::IdentifierCapacityExceeded(kind))
}
