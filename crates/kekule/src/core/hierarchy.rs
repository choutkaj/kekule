use std::collections::BTreeMap;
use std::fmt;

use super::{AtomId, PropMap};

fixed_u32_id!(SmcraChainId);
fixed_u32_id!(SmcraResidueId);
fixed_u32_id!(SmcraAtomSiteId);

/// Coordinate-independent residue, chain, polymer, and atom-site organization.
///
/// A hierarchy owns no atoms or bonds. Atom sites refer to stable `AtomId`
/// values in the graph of the containing molecule.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Hierarchy {
    chains: Vec<SmcraChain>,
    pub(crate) residues: Vec<SmcraResidue>,
    atom_sites: Vec<Option<SmcraAtomSite>>,
    atom_lookup: BTreeMap<AtomId, SmcraAtomSiteId>,
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
    ) -> std::result::Result<SmcraChainId, HierarchyError> {
        self.add_chain_at_slot(label_id.into(), author_id, self.chains.len())
    }

    fn add_chain_at_slot(
        &mut self,
        label_id: String,
        author_id: Option<String>,
        slot: usize,
    ) -> std::result::Result<SmcraChainId, HierarchyError> {
        let id = checked_hierarchy_id(slot, SmcraIdKind::Chain, SmcraChainId::new)?;
        debug_assert_eq!(slot, self.chains.len());
        self.chains.push(SmcraChain {
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
        chain: SmcraChainId,
        name: impl Into<String>,
        label_seq_id: Option<i32>,
        author_seq_id: Option<String>,
        insertion_code: Option<String>,
    ) -> std::result::Result<SmcraResidueId, HierarchyError> {
        self.chain(chain)?;
        let name = name.into();
        let id = checked_hierarchy_id(
            self.residues.len(),
            SmcraIdKind::Residue,
            SmcraResidueId::new,
        )?;
        self.residues.push(SmcraResidue {
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
        residue: SmcraResidueId,
        atom: AtomId,
        metadata: SmcraAtomSiteMetadata,
    ) -> std::result::Result<SmcraAtomSiteId, HierarchyError> {
        self.residue(residue)?;
        if self.atom_lookup.contains_key(&atom) {
            return Err(HierarchyError::DuplicateAtomPlacement(atom));
        }
        let id = checked_hierarchy_id(
            self.atom_sites.len(),
            SmcraIdKind::AtomSite,
            SmcraAtomSiteId::new,
        )?;
        self.atom_sites.push(Some(SmcraAtomSite {
            id,
            residue,
            atom,
            metadata,
            props: PropMap::new(),
        }));
        self.residues[residue.index()].atom_sites.push(id);
        self.atom_lookup.insert(atom, id);
        Ok(id)
    }

    pub fn chain(&self, id: SmcraChainId) -> std::result::Result<&SmcraChain, HierarchyError> {
        self.chains
            .get(id.index())
            .ok_or(HierarchyError::InvalidChainId(id))
    }

    pub fn residue(
        &self,
        id: SmcraResidueId,
    ) -> std::result::Result<&SmcraResidue, HierarchyError> {
        self.residues
            .get(id.index())
            .ok_or(HierarchyError::InvalidResidueId(id))
    }

    pub fn atom_site(
        &self,
        id: SmcraAtomSiteId,
    ) -> std::result::Result<&SmcraAtomSite, HierarchyError> {
        self.atom_sites
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(HierarchyError::InvalidAtomSiteId(id))
    }

    pub fn atom_site_for_atom(&self, atom: AtomId) -> Option<&SmcraAtomSite> {
        self.atom_lookup
            .get(&atom)
            .and_then(|id| self.atom_sites.get(id.index()).and_then(Option::as_ref))
    }

    pub(crate) fn atom_lookup_entries(
        &self,
    ) -> impl Iterator<Item = (AtomId, SmcraAtomSiteId)> + '_ {
        self.atom_lookup.iter().map(|(atom, site)| (*atom, *site))
    }

    /// Removes one atom site without renumbering any hierarchy identifier.
    pub(crate) fn remove_atom_site(
        &mut self,
        id: SmcraAtomSiteId,
    ) -> std::result::Result<SmcraAtomSite, HierarchyError> {
        let site = self
            .atom_sites
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(HierarchyError::InvalidAtomSiteId(id))?;
        let residue = self
            .residues
            .get(site.residue.index())
            .ok_or(HierarchyError::InvalidResidueId(site.residue))?;
        let Some(residue_position) = residue.atom_sites.iter().position(|site| *site == id) else {
            return Err(HierarchyError::InconsistentAtomSite(id));
        };
        if self.atom_lookup.get(&site.atom) != Some(&id) {
            return Err(HierarchyError::InconsistentAtomSite(id));
        }

        let atom = site.atom;
        let removed = self.atom_sites[id.index()]
            .take()
            .expect("validated atom-site must remain live");
        self.residues[removed.residue.index()]
            .atom_sites
            .remove(residue_position);
        self.atom_lookup.remove(&atom);
        Ok(removed)
    }

    /// Restores distinct label and author component identifiers on one residue.
    pub fn set_residue_component_ids(
        &mut self,
        residue: SmcraResidueId,
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
        chain: SmcraChainId,
    ) -> std::result::Result<&mut PropMap, HierarchyError> {
        self.chains
            .get_mut(chain.index())
            .map(|chain| &mut chain.props)
            .ok_or(HierarchyError::InvalidChainId(chain))
    }

    /// Returns checked mutable access to one residue's property map.
    pub fn residue_props_mut(
        &mut self,
        residue: SmcraResidueId,
    ) -> std::result::Result<&mut PropMap, HierarchyError> {
        self.residues
            .get_mut(residue.index())
            .map(|residue| &mut residue.props)
            .ok_or(HierarchyError::InvalidResidueId(residue))
    }

    /// Returns checked mutable access to one atom site's property map.
    pub fn atom_site_props_mut(
        &mut self,
        site: SmcraAtomSiteId,
    ) -> std::result::Result<&mut PropMap, HierarchyError> {
        self.atom_sites
            .get_mut(site.index())
            .and_then(Option::as_mut)
            .map(|site| &mut site.props)
            .ok_or(HierarchyError::InvalidAtomSiteId(site))
    }

    pub fn chains(&self) -> impl Iterator<Item = (SmcraChainId, &SmcraChain)> {
        self.chains.iter().map(|chain| (chain.id, chain))
    }

    pub fn residues(&self) -> impl Iterator<Item = (SmcraResidueId, &SmcraResidue)> {
        self.residues.iter().map(|residue| (residue.id, residue))
    }

    pub fn atom_sites(&self) -> impl Iterator<Item = (SmcraAtomSiteId, &SmcraAtomSite)> {
        self.atom_sites.iter().flatten().map(|site| (site.id, site))
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SmcraChain {
    pub(crate) id: SmcraChainId,
    pub(crate) label_id: String,
    pub(crate) author_id: Option<String>,
    pub(crate) residues: Vec<SmcraResidueId>,
    pub(crate) props: PropMap,
}

impl SmcraChain {
    pub const fn id(&self) -> SmcraChainId {
        self.id
    }

    pub fn label_id(&self) -> &str {
        &self.label_id
    }

    pub fn author_id(&self) -> Option<&str> {
        self.author_id.as_deref()
    }

    pub fn residues(&self) -> &[SmcraResidueId] {
        &self.residues
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SmcraResidue {
    pub(crate) id: SmcraResidueId,
    pub(crate) chain: SmcraChainId,
    pub(crate) name: String,
    pub(crate) label_comp_id: Option<String>,
    pub(crate) author_comp_id: Option<String>,
    pub(crate) label_seq_id: Option<i32>,
    pub(crate) author_seq_id: Option<String>,
    pub(crate) insertion_code: Option<String>,
    pub(crate) atom_sites: Vec<SmcraAtomSiteId>,
    pub(crate) props: PropMap,
}

impl SmcraResidue {
    pub const fn id(&self) -> SmcraResidueId {
        self.id
    }

    pub const fn chain(&self) -> SmcraChainId {
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

    pub fn atom_sites(&self) -> &[SmcraAtomSiteId] {
        &self.atom_sites
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SmcraAtomSite {
    pub(crate) id: SmcraAtomSiteId,
    pub(crate) residue: SmcraResidueId,
    pub(crate) atom: AtomId,
    pub(crate) metadata: SmcraAtomSiteMetadata,
    pub(crate) props: PropMap,
}

impl SmcraAtomSite {
    pub const fn id(&self) -> SmcraAtomSiteId {
        self.id
    }

    pub const fn residue(&self) -> SmcraResidueId {
        self.residue
    }

    pub const fn atom(&self) -> AtomId {
        self.atom
    }

    pub fn metadata(&self) -> &SmcraAtomSiteMetadata {
        &self.metadata
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SmcraAtomSiteMetadata {
    pub type_symbol: Option<String>,
    pub label_asym_id: Option<String>,
    pub auth_asym_id: Option<String>,
    pub label_atom_id: Option<String>,
    pub auth_atom_id: Option<String>,
}

/// Fixed-width identifier spaces owned by [`Hierarchy`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SmcraIdKind {
    /// Chain identifiers.
    Chain,
    /// Residue identifiers.
    Residue,
    /// Atom-site identifiers.
    AtomSite,
}

impl fmt::Display for SmcraIdKind {
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
    InvalidChainId(SmcraChainId),
    InvalidResidueId(SmcraResidueId),
    InvalidAtomSiteId(SmcraAtomSiteId),
    InvalidAtomId(AtomId),
    DuplicateAtomPlacement(AtomId),
    /// An atom-site slot disagrees with its residue or atom lookup.
    InconsistentAtomSite(SmcraAtomSiteId),
    /// A new hierarchy node cannot be represented by the fixed-width ID for `kind`.
    IdentifierCapacityExceeded(SmcraIdKind),
}

impl fmt::Display for HierarchyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidChainId(id) => write!(f, "invalid chain id: {}", id.raw()),
            Self::InvalidResidueId(id) => write!(f, "invalid residue id: {}", id.raw()),
            Self::InvalidAtomSiteId(id) => write!(f, "invalid atom-site id: {}", id.raw()),
            Self::InvalidAtomId(id) => write!(f, "invalid hierarchy atom id: {id}"),
            Self::DuplicateAtomPlacement(id) => write!(f, "duplicate hierarchy placement for {id}"),
            Self::InconsistentAtomSite(id) => {
                write!(f, "inconsistent hierarchy atom-site: {}", id.raw())
            }
            Self::IdentifierCapacityExceeded(kind) => {
                write!(f, "{kind} identifier capacity exceeded")
            }
        }
    }
}

impl std::error::Error for HierarchyError {}

fn checked_hierarchy_id<T>(
    length: usize,
    kind: SmcraIdKind,
    construct: impl FnOnce(u32) -> T,
) -> std::result::Result<T, HierarchyError> {
    crate::core::checked_raw_id(length)
        .map(construct)
        .map_err(|_| HierarchyError::IdentifierCapacityExceeded(kind))
}
