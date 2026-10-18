use std::collections::BTreeMap;
use std::fmt;

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
#[derive(Debug, Clone, PartialEq, Default)]
pub(crate) struct Hierarchy {
    chains: Vec<Chain>,
    pub(crate) residues: Vec<Residue>,
    atom_sites: Vec<AtomSite>,
    atom_lookup: BTreeMap<InstanceAtomId, AtomSiteId>,
}

impl Hierarchy {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.chains.is_empty() && self.residues.is_empty() && self.atom_sites.is_empty()
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
        let name = name.into();
        let id = checked_hierarchy_id(
            self.residues.len(),
            HierarchyIdKind::Residue,
            ResidueId::new,
        )?;
        self.residues.push(Residue {
            id,
            chain,
            class: ResidueClass::Other,
            explicit_class: None,
            name: name.clone(),
            label_comp_id: Some(name),
            author_comp_id: None,
            label_seq_id,
            author_seq_id,
            insertion_code,
            atom_sites: Vec::new(),
        });
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
        if self.atom_lookup.contains_key(&atom) {
            return Err(HierarchyError::DuplicateAtomPlacement(atom));
        }
        let id = checked_hierarchy_id(
            self.atom_sites.len(),
            HierarchyIdKind::AtomSite,
            AtomSiteId::new,
        )?;
        self.atom_sites.push(AtomSite {
            id,
            residue,
            atom,
            metadata,
        });
        self.residues[residue.index()].atom_sites.push(id);
        self.atom_lookup.insert(atom, id);
        Ok(id)
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
        hierarchy.residues = residues
            .into_iter()
            .zip(residue_chain)
            .enumerate()
            .map(|(index, (residue, chain))| {
                let id = ResidueId::new(index as u32);
                Ok(Residue {
                    id,
                    chain: chain.ok_or(RawHierarchyError::ResidueUnlisted(id))?,
                    class: residue.class,
                    explicit_class: residue.explicit.then_some(residue.class),
                    name: residue.name,
                    label_comp_id: residue.label_comp_id,
                    author_comp_id: residue.author_comp_id,
                    label_seq_id: residue.label_seq_id,
                    author_seq_id: residue.author_seq_id,
                    insertion_code: residue.insertion_code,
                    atom_sites: residue.atom_sites,
                })
            })
            .collect::<std::result::Result<_, _>>()?;
        let mut atom_sites = Vec::with_capacity(sites.len());
        let mut lookup = Vec::with_capacity(sites.len());
        for (index, (site, residue)) in sites.into_iter().zip(site_residue).enumerate() {
            let id = AtomSiteId::new(index as u32);
            let residue = residue.ok_or(RawHierarchyError::AtomSiteUnlisted(id))?;
            lookup.push((site.atom, id));
            atom_sites.push(AtomSite {
                id,
                residue,
                atom: site.atom,
                metadata: site.metadata,
            });
        }
        // Sorted entries build the lookup in bulk rather than one insertion
        // at a time.
        lookup.sort_unstable();
        if let Some(pair) = lookup.windows(2).find(|pair| pair[0].0 == pair[1].0) {
            return Err(RawHierarchyError::DuplicateAtomPlacement(pair[0].0));
        }
        hierarchy.atom_lookup = lookup.into_iter().collect();
        hierarchy.atom_sites = atom_sites;
        Ok(hierarchy)
    }

    /// Mutable residue list of one chain, for corrupting hierarchies in tests.
    #[cfg(test)]
    pub(crate) fn chain_residues_mut(&mut self, chain: ChainId) -> &mut Vec<ResidueId> {
        &mut self.chains[chain.index()].residues
    }

    /// One chain, or `None` when `id` names no chain of this hierarchy.
    pub(crate) fn chain(&self, id: ChainId) -> Option<&Chain> {
        self.chains.get(id.index())
    }

    /// One residue, or `None` when `id` names no residue of this hierarchy.
    pub(crate) fn residue(&self, id: ResidueId) -> Option<&Residue> {
        self.residues.get(id.index())
    }

    /// One atom site, or `None` when `id` names no site of this hierarchy.
    pub(crate) fn atom_site(&self, id: AtomSiteId) -> Option<&AtomSite> {
        self.atom_sites.get(id.index())
    }

    /// [`Self::chain`] for staging code that reports a missing chain.
    pub(crate) fn require_chain(&self, id: ChainId) -> std::result::Result<&Chain, HierarchyError> {
        self.chain(id).ok_or(HierarchyError::InvalidChainId(id))
    }

    /// [`Self::residue`] for staging code that reports a missing residue.
    pub(crate) fn require_residue(
        &self,
        id: ResidueId,
    ) -> std::result::Result<&Residue, HierarchyError> {
        self.residue(id).ok_or(HierarchyError::InvalidResidueId(id))
    }

    /// [`Self::atom_site`] for staging code that reports a missing site.
    pub(crate) fn require_atom_site(
        &self,
        id: AtomSiteId,
    ) -> std::result::Result<&AtomSite, HierarchyError> {
        self.atom_site(id)
            .ok_or(HierarchyError::InvalidAtomSiteId(id))
    }

    fn residue_mut(&mut self, id: ResidueId) -> std::result::Result<&mut Residue, HierarchyError> {
        self.residues
            .get_mut(id.index())
            .ok_or(HierarchyError::InvalidResidueId(id))
    }

    pub(crate) fn atom_site_for_atom(&self, atom: InstanceAtomId) -> Option<&AtomSite> {
        self.atom_lookup
            .get(&atom)
            .and_then(|id| self.atom_sites.get(id.index()))
    }

    pub(crate) fn atom_lookup_entries(
        &self,
    ) -> impl Iterator<Item = (InstanceAtomId, AtomSiteId)> + '_ {
        self.atom_lookup.iter().map(|(atom, site)| (*atom, *site))
    }

    /// Restores distinct label and author component identifiers on one residue.
    pub(crate) fn set_residue_component_ids(
        &mut self,
        residue: ResidueId,
        label_comp_id: Option<String>,
        author_comp_id: Option<String>,
    ) -> std::result::Result<(), HierarchyError> {
        let residue = self.residue_mut(residue)?;
        residue.label_comp_id = label_comp_id;
        residue.author_comp_id = author_comp_id;
        Ok(())
    }

    /// Records or clears an explicit class on one residue node.
    pub(crate) fn set_explicit_residue_class(
        &mut self,
        residue: ResidueId,
        class: Option<ResidueClass>,
    ) -> std::result::Result<(), HierarchyError> {
        self.residue_mut(residue)?.explicit_class = class;
        Ok(())
    }

    /// Equal chains, residues, and atom sites, ignoring explicit-class
    /// intent, which is editing policy rather than layout.
    pub(crate) fn same_layout(&self, other: &Self) -> bool {
        self.chains == other.chains
            && self.atom_sites == other.atom_sites
            && self.residues.len() == other.residues.len()
            && self.residues.iter().zip(&other.residues).all(|(a, b)| {
                Residue {
                    explicit_class: b.explicit_class,
                    ..a.clone()
                } == *b
            })
    }

    pub(crate) fn chains(&self) -> impl ExactSizeIterator<Item = (ChainId, &Chain)> {
        self.chains.iter().map(|chain| (chain.id, chain))
    }

    pub(crate) fn residues(&self) -> impl ExactSizeIterator<Item = (ResidueId, &Residue)> {
        self.residues.iter().map(|residue| (residue.id, residue))
    }

    pub(crate) fn atom_sites(&self) -> impl ExactSizeIterator<Item = (AtomSiteId, &AtomSite)> {
        self.atom_sites.iter().map(|site| (site.id, site))
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
        self.hierarchy.residues.len()
    }

    /// The number of staged atom sites.
    pub fn atom_site_count(&self) -> usize {
        self.hierarchy.atom_sites.len()
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
        self.hierarchy.residue(residue)?.explicit_class
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

/// One hierarchy residue: its chain, class, source identifiers, and atom sites.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Residue {
    pub(crate) id: ResidueId,
    pub(crate) chain: ChainId,
    pub(crate) class: ResidueClass,
    /// The explicitly assigned class, which publication keeps.
    pub(crate) explicit_class: Option<ResidueClass>,
    pub(crate) name: String,
    pub(crate) label_comp_id: Option<String>,
    pub(crate) author_comp_id: Option<String>,
    pub(crate) label_seq_id: Option<i32>,
    pub(crate) author_seq_id: Option<String>,
    pub(crate) insertion_code: Option<String>,
    pub(crate) atom_sites: Vec<AtomSiteId>,
}

impl Residue {
    pub(crate) const fn id(&self) -> ResidueId {
        self.id
    }

    pub(crate) const fn chain(&self) -> ChainId {
        self.chain
    }

    pub(crate) const fn class(&self) -> ResidueClass {
        self.class
    }

    pub(crate) fn name(&self) -> &str {
        &self.name
    }

    pub(crate) fn label_comp_id(&self) -> Option<&str> {
        self.label_comp_id.as_deref()
    }

    pub(crate) fn author_comp_id(&self) -> Option<&str> {
        self.author_comp_id.as_deref()
    }

    pub(crate) const fn label_seq_id(&self) -> Option<i32> {
        self.label_seq_id
    }

    pub(crate) fn author_seq_id(&self) -> Option<&str> {
        self.author_seq_id.as_deref()
    }

    pub(crate) fn insertion_code(&self) -> Option<&str> {
        self.insertion_code.as_deref()
    }

    pub(crate) fn atom_sites(&self) -> &[AtomSiteId] {
        &self.atom_sites
    }
}

/// One hierarchy atom site: the placement of one topology atom in a residue.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AtomSite {
    pub(crate) id: AtomSiteId,
    pub(crate) residue: ResidueId,
    pub(crate) atom: InstanceAtomId,
    pub(crate) metadata: AtomSiteMetadata,
}

impl AtomSite {
    pub(crate) const fn id(&self) -> AtomSiteId {
        self.id
    }

    pub(crate) const fn residue(&self) -> ResidueId {
        self.residue
    }

    pub(crate) const fn atom(&self) -> InstanceAtomId {
        self.atom
    }

    pub(crate) fn metadata(&self) -> &AtomSiteMetadata {
        &self.metadata
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
