use std::collections::BTreeMap;
use std::fmt;

use super::InstanceAtomId;
use crate::core::PropMap;

fixed_u32_id!(ChainId, "chain");
fixed_u32_id!(ResidueId, "residue");
fixed_u32_id!(AtomSiteId, "atom-site");

/// Coordinate-independent residue, chain, polymer, and atom-site organization.
///
/// A hierarchy owns no atoms or bonds. Atom sites refer to topology-global
/// [`InstanceAtomId`] values and may span molecule-instance boundaries.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Hierarchy {
    chains: Vec<Chain>,
    pub(crate) residues: Vec<Residue>,
    atom_sites: Vec<AtomSite>,
    atom_lookup: BTreeMap<InstanceAtomId, AtomSiteId>,
    props: PropMap,
}

impl Hierarchy {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.chains.is_empty() && self.residues.is_empty() && self.atom_sites.is_empty()
    }

    /// Adds a chain or returns a structured capacity error without mutation.
    pub fn add_chain(
        &mut self,
        label_id: impl Into<String>,
        author_id: Option<String>,
    ) -> std::result::Result<ChainId, HierarchyError> {
        self.add_chain_at_slot(label_id.into(), author_id, self.chains.len())
    }

    fn add_chain_at_slot(
        &mut self,
        label_id: String,
        author_id: Option<String>,
        slot: usize,
    ) -> std::result::Result<ChainId, HierarchyError> {
        let id = checked_hierarchy_id(slot, HierarchyIdKind::Chain, ChainId::new)?;
        debug_assert_eq!(slot, self.chains.len());
        self.chains.push(Chain {
            id,
            label_id,
            author_id,
            residues: Vec::new(),
            props: PropMap::new(),
        });
        Ok(id)
    }

    /// Adds a residue to an existing chain transactionally.
    pub fn add_residue(
        &mut self,
        chain: ChainId,
        name: impl Into<String>,
        label_seq_id: Option<i32>,
        author_seq_id: Option<String>,
        insertion_code: Option<String>,
    ) -> std::result::Result<ResidueId, HierarchyError> {
        self.chain(chain)?;
        let name = name.into();
        let id = checked_hierarchy_id(
            self.residues.len(),
            HierarchyIdKind::Residue,
            ResidueId::new,
        )?;
        self.residues.push(Residue {
            id,
            chain,
            name: name.clone(),
            label_comp_id: Some(name),
            author_comp_id: None,
            label_seq_id,
            author_seq_id,
            insertion_code,
            atom_sites: Vec::new(),
            props: PropMap::new(),
        });
        self.chains[chain.index()].residues.push(id);
        Ok(id)
    }

    /// Adds one atom site after validating residue, atom placement, and capacity.
    pub fn add_atom_site(
        &mut self,
        residue: ResidueId,
        atom: InstanceAtomId,
        metadata: AtomSiteMetadata,
    ) -> std::result::Result<AtomSiteId, HierarchyError> {
        self.residue(residue)?;
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
            props: PropMap::new(),
        });
        self.residues[residue.index()].atom_sites.push(id);
        self.atom_lookup.insert(atom, id);
        Ok(id)
    }

    pub fn chain(&self, id: ChainId) -> std::result::Result<&Chain, HierarchyError> {
        self.chains
            .get(id.index())
            .ok_or(HierarchyError::InvalidChainId(id))
    }

    pub fn residue(&self, id: ResidueId) -> std::result::Result<&Residue, HierarchyError> {
        self.residues
            .get(id.index())
            .ok_or(HierarchyError::InvalidResidueId(id))
    }

    pub fn atom_site(&self, id: AtomSiteId) -> std::result::Result<&AtomSite, HierarchyError> {
        self.atom_sites
            .get(id.index())
            .ok_or(HierarchyError::InvalidAtomSiteId(id))
    }

    pub fn atom_site_for_atom(&self, atom: InstanceAtomId) -> Option<&AtomSite> {
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
    pub fn set_residue_component_ids(
        &mut self,
        residue: ResidueId,
        label_comp_id: Option<String>,
        author_comp_id: Option<String>,
    ) -> std::result::Result<(), HierarchyError> {
        let residue = self
            .residues
            .get_mut(residue.index())
            .ok_or(HierarchyError::InvalidResidueId(residue))?;
        residue.label_comp_id = label_comp_id;
        residue.author_comp_id = author_comp_id;
        Ok(())
    }

    /// Returns checked mutable access to one chain's property map.
    pub fn chain_props_mut(
        &mut self,
        chain: ChainId,
    ) -> std::result::Result<&mut PropMap, HierarchyError> {
        self.chains
            .get_mut(chain.index())
            .map(|chain| &mut chain.props)
            .ok_or(HierarchyError::InvalidChainId(chain))
    }

    /// Returns checked mutable access to one residue's property map.
    pub fn residue_props_mut(
        &mut self,
        residue: ResidueId,
    ) -> std::result::Result<&mut PropMap, HierarchyError> {
        self.residues
            .get_mut(residue.index())
            .map(|residue| &mut residue.props)
            .ok_or(HierarchyError::InvalidResidueId(residue))
    }

    /// Returns checked mutable access to one atom site's property map.
    pub fn atom_site_props_mut(
        &mut self,
        site: AtomSiteId,
    ) -> std::result::Result<&mut PropMap, HierarchyError> {
        self.atom_sites
            .get_mut(site.index())
            .map(|site| &mut site.props)
            .ok_or(HierarchyError::InvalidAtomSiteId(site))
    }

    pub fn chains(&self) -> impl Iterator<Item = (ChainId, &Chain)> {
        self.chains.iter().map(|chain| (chain.id, chain))
    }

    pub fn residues(&self) -> impl Iterator<Item = (ResidueId, &Residue)> {
        self.residues.iter().map(|residue| (residue.id, residue))
    }

    pub fn atom_sites(&self) -> impl Iterator<Item = (AtomSiteId, &AtomSite)> {
        self.atom_sites.iter().map(|site| (site.id, site))
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Chain {
    pub(crate) id: ChainId,
    pub(crate) label_id: String,
    pub(crate) author_id: Option<String>,
    pub(crate) residues: Vec<ResidueId>,
    pub(crate) props: PropMap,
}

impl Chain {
    pub const fn id(&self) -> ChainId {
        self.id
    }

    pub fn label_id(&self) -> &str {
        &self.label_id
    }

    pub fn author_id(&self) -> Option<&str> {
        self.author_id.as_deref()
    }

    pub fn residues(&self) -> &[ResidueId] {
        &self.residues
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub(crate) id: ResidueId,
    pub(crate) chain: ChainId,
    pub(crate) name: String,
    pub(crate) label_comp_id: Option<String>,
    pub(crate) author_comp_id: Option<String>,
    pub(crate) label_seq_id: Option<i32>,
    pub(crate) author_seq_id: Option<String>,
    pub(crate) insertion_code: Option<String>,
    pub(crate) atom_sites: Vec<AtomSiteId>,
    pub(crate) props: PropMap,
}

impl Residue {
    pub const fn id(&self) -> ResidueId {
        self.id
    }

    pub const fn chain(&self) -> ChainId {
        self.chain
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn label_comp_id(&self) -> Option<&str> {
        self.label_comp_id.as_deref()
    }

    pub fn author_comp_id(&self) -> Option<&str> {
        self.author_comp_id.as_deref()
    }

    pub const fn label_seq_id(&self) -> Option<i32> {
        self.label_seq_id
    }

    pub fn author_seq_id(&self) -> Option<&str> {
        self.author_seq_id.as_deref()
    }

    pub fn insertion_code(&self) -> Option<&str> {
        self.insertion_code.as_deref()
    }

    pub fn atom_sites(&self) -> &[AtomSiteId] {
        &self.atom_sites
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomSite {
    pub(crate) id: AtomSiteId,
    pub(crate) residue: ResidueId,
    pub(crate) atom: InstanceAtomId,
    pub(crate) metadata: AtomSiteMetadata,
    pub(crate) props: PropMap,
}

impl AtomSite {
    pub const fn id(&self) -> AtomSiteId {
        self.id
    }

    pub const fn residue(&self) -> ResidueId {
        self.residue
    }

    pub const fn atom(&self) -> InstanceAtomId {
        self.atom
    }

    pub fn metadata(&self) -> &AtomSiteMetadata {
        &self.metadata
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct AtomSiteMetadata {
    pub type_symbol: Option<String>,
    pub label_asym_id: Option<String>,
    pub auth_asym_id: Option<String>,
    pub label_atom_id: Option<String>,
    pub auth_atom_id: Option<String>,
}

/// Fixed-width identifier spaces owned by [`Hierarchy`].
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

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HierarchyError {
    InvalidChainId(ChainId),
    InvalidResidueId(ResidueId),
    InvalidAtomSiteId(AtomSiteId),
    DuplicateAtomPlacement(InstanceAtomId),
    /// A new hierarchy node cannot be represented by the fixed-width ID for `kind`.
    IdentifierCapacityExceeded(HierarchyIdKind),
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
        }
    }
}

impl std::error::Error for HierarchyError {}

fn checked_hierarchy_id<T>(
    length: usize,
    kind: HierarchyIdKind,
    construct: impl FnOnce(u32) -> T,
) -> std::result::Result<T, HierarchyError> {
    crate::core::checked_raw_id(length)
        .map(construct)
        .map_err(|_| HierarchyError::IdentifierCapacityExceeded(kind))
}
