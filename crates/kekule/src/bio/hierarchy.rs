use std::collections::BTreeMap;
use std::fmt;

use crate::chemistry::{perceive_molecule, PerceptionError};
use crate::core::*;

/// One connected macromolecular graph with its coordinated SMCRA sidecar.
///
/// ```compile_fail
/// use kekule::bio::MacroMolecule;
///
/// let mut molecule = MacroMolecule::new();
/// let _ = molecule.as_molecule_mut();
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MacroMolecule {
    molecule: Molecule,
    hierarchy: SmcraHierarchy,
}

impl MacroMolecule {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn builder() -> MacroMoleculeBuilder {
        MacroMoleculeBuilder::new()
    }

    pub fn from_parts(
        mut molecule: Molecule,
        hierarchy: SmcraHierarchy,
    ) -> std::result::Result<Self, MacroValidateError> {
        canonicalize_represented_chemistry(&mut molecule).map_err(|error| {
            MacroValidateError::CanonicalRepresentation {
                atom: error.atom,
                charge: error.charge,
            }
        })?;
        let macro_molecule = Self {
            molecule,
            hierarchy,
        };
        macro_molecule.validate()?;
        Ok(macro_molecule)
    }

    /// Builds private interpretation staging state before component partitioning.
    pub(crate) fn from_parts_unchecked_connectedness(
        molecule: Molecule,
        hierarchy: SmcraHierarchy,
    ) -> std::result::Result<Self, MacroValidateError> {
        let macro_molecule = Self {
            molecule,
            hierarchy,
        };
        validate_macro_molecule_contents(&macro_molecule, MacroValidateOptions::default())?;
        Ok(macro_molecule)
    }

    pub fn as_molecule(&self) -> &Molecule {
        &self.molecule
    }

    pub(crate) fn molecule_mut_unchecked_connectedness(&mut self) -> &mut Molecule {
        &mut self.molecule
    }

    pub fn to_parts(self) -> (Molecule, SmcraHierarchy) {
        (self.molecule, self.hierarchy)
    }

    pub fn hierarchy(&self) -> &SmcraHierarchy {
        &self.hierarchy
    }

    pub fn edit(&mut self) -> MacroMoleculeEditor<'_> {
        MacroMoleculeEditor {
            molecule: self.molecule.clone(),
            hierarchy: self.hierarchy.clone(),
            target: self,
        }
    }

    pub(crate) fn without_conformers(mut self) -> Self {
        self.molecule = self.molecule.without_conformers();
        self
    }

    pub(crate) fn clone_without_conformers(&self) -> Self {
        Self {
            molecule: self.molecule.clone_without_conformers(),
            hierarchy: self.hierarchy.clone(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.molecule.is_empty()
    }

    pub fn atom_count(&self) -> usize {
        self.molecule.atom_count()
    }

    pub fn bond_count(&self) -> usize {
        self.molecule.bond_count()
    }

    pub fn formal_charge(&self) -> i64 {
        self.molecule.formal_charge()
    }

    pub fn atom(&self, id: AtomId) -> crate::core::Result<&Atom> {
        self.molecule.atom(id)
    }

    pub fn bond(&self, id: BondId) -> crate::core::Result<&Bond> {
        self.molecule.bond(id)
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.molecule.atoms()
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)> {
        self.molecule.bonds()
    }

    pub fn perceive(&mut self) -> std::result::Result<(), PerceptionError> {
        perceive_molecule(&mut self.molecule)
    }

    pub fn chains(&self) -> impl Iterator<Item = (SmcraChainId, &SmcraChain)> {
        self.hierarchy.chains()
    }

    pub fn residues(&self) -> impl Iterator<Item = (SmcraResidueId, &SmcraResidue)> {
        self.hierarchy.residues()
    }

    pub fn atom_sites(&self) -> impl Iterator<Item = (SmcraAtomSiteId, &SmcraAtomSite)> {
        self.hierarchy.atom_sites()
    }

    pub fn atom_site_for_atom(&self, atom: AtomId) -> Option<&SmcraAtomSite> {
        self.hierarchy.atom_site_for_atom(atom)
    }

    pub fn validate(&self) -> std::result::Result<MacroValidateReport, MacroValidateError> {
        self.validate_with_options(MacroValidateOptions::default())
    }

    pub fn validate_with_options(
        &self,
        options: MacroValidateOptions,
    ) -> std::result::Result<MacroValidateReport, MacroValidateError> {
        validate_macro_molecule(self, options)
    }
}

impl AsRef<Molecule> for MacroMolecule {
    fn as_ref(&self) -> &Molecule {
        self.as_molecule()
    }
}

/// Mutable staging for one final connected [`MacroMolecule`].
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MacroMoleculeBuilder {
    molecule: Molecule,
    hierarchy: SmcraHierarchy,
}

impl MacroMoleculeBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn as_molecule(&self) -> &Molecule {
        &self.molecule
    }

    pub fn as_molecule_mut(&mut self) -> &mut Molecule {
        &mut self.molecule
    }

    pub fn hierarchy(&self) -> &SmcraHierarchy {
        &self.hierarchy
    }

    pub fn hierarchy_mut(&mut self) -> &mut SmcraHierarchy {
        &mut self.hierarchy
    }

    pub fn add_atom_site(
        &mut self,
        residue: SmcraResidueId,
        atom: AtomId,
        metadata: SmcraAtomSiteMetadata,
    ) -> std::result::Result<SmcraAtomSiteId, SmcraHierarchyError> {
        self.molecule
            .atom(atom)
            .map_err(|_| SmcraHierarchyError::InvalidAtomId(atom))?;
        self.hierarchy.add_atom_site(residue, atom, metadata)
    }

    pub fn build(self) -> std::result::Result<MacroMolecule, MacroValidateError> {
        MacroMolecule::from_parts(self.molecule, self.hierarchy)
    }
}

/// Transactional coordinated editing for a connected [`MacroMolecule`].
pub struct MacroMoleculeEditor<'a> {
    target: &'a mut MacroMolecule,
    molecule: Molecule,
    hierarchy: SmcraHierarchy,
}

impl MacroMoleculeEditor<'_> {
    pub fn as_molecule(&self) -> &Molecule {
        &self.molecule
    }

    pub fn as_molecule_mut(&mut self) -> &mut Molecule {
        &mut self.molecule
    }

    /// Returns mutable represented atom state in this checked working copy.
    pub fn atom_mut(&mut self, atom: AtomId) -> crate::core::Result<AtomMut<'_>> {
        self.molecule.atom_mut(atom)
    }

    /// Returns mutable represented bond state in this checked working copy.
    pub fn bond_mut(&mut self, bond: BondId) -> crate::core::Result<BondMut<'_>> {
        self.molecule.bond_mut(bond)
    }

    pub fn hierarchy(&self) -> &SmcraHierarchy {
        &self.hierarchy
    }

    pub fn hierarchy_mut(&mut self) -> &mut SmcraHierarchy {
        &mut self.hierarchy
    }

    pub fn add_atom_site(
        &mut self,
        residue: SmcraResidueId,
        atom: AtomId,
        metadata: SmcraAtomSiteMetadata,
    ) -> std::result::Result<SmcraAtomSiteId, SmcraHierarchyError> {
        self.molecule
            .atom(atom)
            .map_err(|_| SmcraHierarchyError::InvalidAtomId(atom))?;
        self.hierarchy.add_atom_site(residue, atom, metadata)
    }

    pub fn commit(self) -> std::result::Result<(), MacroValidateError> {
        let candidate = MacroMolecule::from_parts(self.molecule, self.hierarchy)?;
        *self.target = candidate;
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MacroValidateOptions {
    pub validate_coordinates: bool,
}

impl Default for MacroValidateOptions {
    fn default() -> Self {
        Self {
            validate_coordinates: true,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MacroValidateReport {
    pub chains_checked: usize,
    pub residues_checked: usize,
    pub atom_sites_checked: usize,
    pub conformers_checked: usize,
    pub coordinates_checked: usize,
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MacroValidateError {
    DisconnectedGraph(MoleculeConnectivityError),
    CanonicalRepresentation {
        atom: AtomId,
        charge: usize,
    },
    InvalidResidueChain {
        residue: SmcraResidueId,
        chain: SmcraChainId,
    },
    InvalidResidueAtomSite {
        residue: SmcraResidueId,
        site: SmcraAtomSiteId,
    },
    InvalidAtomSiteResidue {
        site: SmcraAtomSiteId,
        residue: SmcraResidueId,
    },
    InvalidAtomSiteAtom {
        site: SmcraAtomSiteId,
        atom: AtomId,
    },
    MissingAtomSiteForAtom {
        atom: AtomId,
    },
    InvalidConformerAtom {
        conformer: ConformerId,
        atom: AtomId,
    },
}

impl fmt::Display for MacroValidateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DisconnectedGraph(error) => write!(f, "invalid macromolecule graph: {error}"),
            Self::CanonicalRepresentation { atom, charge } => write!(
                f,
                "publishing macromolecule atom {atom} requires formal charge +{charge}, which is outside the supported range"
            ),
            Self::InvalidResidueChain { residue, chain } => write!(
                f,
                "residue {} references invalid chain {}",
                residue.raw(),
                chain.raw()
            ),
            Self::InvalidResidueAtomSite { residue, site } => write!(
                f,
                "residue {} references invalid atom-site {}",
                residue.raw(),
                site.raw()
            ),
            Self::InvalidAtomSiteResidue { site, residue } => write!(
                f,
                "atom-site {} references invalid residue {}",
                site.raw(),
                residue.raw()
            ),
            Self::InvalidAtomSiteAtom { site, atom } => {
                write!(f, "atom-site {} references invalid atom {atom}", site.raw())
            }
            Self::MissingAtomSiteForAtom { atom } => {
                write!(f, "macro-molecule atom {atom} has no hierarchy atom-site")
            }
            Self::InvalidConformerAtom { conformer, atom } => write!(
                f,
                "conformer {} stores coordinates for invalid atom {atom}",
                conformer.raw()
            ),
        }
    }
}

impl std::error::Error for MacroValidateError {}

fn validate_macro_molecule(
    molecule: &MacroMolecule,
    options: MacroValidateOptions,
) -> std::result::Result<MacroValidateReport, MacroValidateError> {
    molecule
        .molecule
        .validate_connected()
        .map_err(MacroValidateError::DisconnectedGraph)?;
    validate_macro_molecule_contents(molecule, options)
}

fn validate_macro_molecule_contents(
    molecule: &MacroMolecule,
    options: MacroValidateOptions,
) -> std::result::Result<MacroValidateReport, MacroValidateError> {
    let mut report = MacroValidateReport {
        chains_checked: 0,
        residues_checked: 0,
        atom_sites_checked: 0,
        conformers_checked: 0,
        coordinates_checked: 0,
    };

    for _ in molecule.hierarchy.chains() {
        report.chains_checked += 1;
    }
    for (residue_id, residue) in molecule.hierarchy.residues() {
        molecule.hierarchy.chain(residue.chain).map_err(|_| {
            MacroValidateError::InvalidResidueChain {
                residue: residue_id,
                chain: residue.chain,
            }
        })?;
        for site in &residue.atom_sites {
            molecule.hierarchy.atom_site(*site).map_err(|_| {
                MacroValidateError::InvalidResidueAtomSite {
                    residue: residue_id,
                    site: *site,
                }
            })?;
        }
        report.residues_checked += 1;
    }
    for (site_id, site) in molecule.hierarchy.atom_sites() {
        molecule.hierarchy.residue(site.residue).map_err(|_| {
            MacroValidateError::InvalidAtomSiteResidue {
                site: site_id,
                residue: site.residue,
            }
        })?;
        molecule
            .molecule
            .atom(site.atom)
            .map_err(|_| MacroValidateError::InvalidAtomSiteAtom {
                site: site_id,
                atom: site.atom,
            })?;
        report.atom_sites_checked += 1;
    }
    for atom in molecule.molecule.atom_ids() {
        if molecule.hierarchy.atom_site_for_atom(atom).is_none() {
            return Err(MacroValidateError::MissingAtomSiteForAtom { atom });
        }
    }
    if options.validate_coordinates {
        for (conformer_id, conformer) in molecule.molecule.conformers() {
            report.conformers_checked += 1;
            for (atom, point) in conformer.positions() {
                let point = point.value();
                molecule.molecule.atom(atom).map_err(|_| {
                    MacroValidateError::InvalidConformerAtom {
                        conformer: conformer_id,
                        atom,
                    }
                })?;
                if point.x.is_finite() && point.y.is_finite() && point.z.is_finite() {
                    report.coordinates_checked += 1;
                } else {
                    return Err(MacroValidateError::InvalidConformerAtom {
                        conformer: conformer_id,
                        atom,
                    });
                }
            }
        }
    }
    Ok(report)
}

fixed_u32_id!(SmcraChainId);
fixed_u32_id!(SmcraResidueId);
fixed_u32_id!(SmcraAtomSiteId);

#[derive(Debug, Clone, PartialEq, Default)]
pub struct SmcraHierarchy {
    chains: Vec<SmcraChain>,
    pub(crate) residues: Vec<SmcraResidue>,
    atom_sites: Vec<SmcraAtomSite>,
    atom_lookup: BTreeMap<AtomId, SmcraAtomSiteId>,
    props: PropMap,
}

impl SmcraHierarchy {
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a chain or returns a structured capacity error without mutation.
    pub fn add_chain(
        &mut self,
        label_id: impl Into<String>,
        author_id: Option<String>,
    ) -> std::result::Result<SmcraChainId, SmcraHierarchyError> {
        self.add_chain_at_slot(label_id.into(), author_id, self.chains.len())
    }

    fn add_chain_at_slot(
        &mut self,
        label_id: String,
        author_id: Option<String>,
        slot: usize,
    ) -> std::result::Result<SmcraChainId, SmcraHierarchyError> {
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
    ) -> std::result::Result<SmcraResidueId, SmcraHierarchyError> {
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
    ) -> std::result::Result<SmcraAtomSiteId, SmcraHierarchyError> {
        self.residue(residue)?;
        if self.atom_lookup.contains_key(&atom) {
            return Err(SmcraHierarchyError::DuplicateAtomPlacement(atom));
        }
        let id = checked_hierarchy_id(
            self.atom_sites.len(),
            SmcraIdKind::AtomSite,
            SmcraAtomSiteId::new,
        )?;
        self.atom_sites.push(SmcraAtomSite {
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

    pub fn chain(&self, id: SmcraChainId) -> std::result::Result<&SmcraChain, SmcraHierarchyError> {
        self.chains
            .get(id.index())
            .ok_or(SmcraHierarchyError::InvalidChainId(id))
    }

    pub fn residue(
        &self,
        id: SmcraResidueId,
    ) -> std::result::Result<&SmcraResidue, SmcraHierarchyError> {
        self.residues
            .get(id.index())
            .ok_or(SmcraHierarchyError::InvalidResidueId(id))
    }

    pub fn atom_site(
        &self,
        id: SmcraAtomSiteId,
    ) -> std::result::Result<&SmcraAtomSite, SmcraHierarchyError> {
        self.atom_sites
            .get(id.index())
            .ok_or(SmcraHierarchyError::InvalidAtomSiteId(id))
    }

    pub fn atom_site_for_atom(&self, atom: AtomId) -> Option<&SmcraAtomSite> {
        self.atom_lookup
            .get(&atom)
            .and_then(|id| self.atom_sites.get(id.index()))
    }

    /// Restores distinct label and author component identifiers on one residue.
    pub fn set_residue_component_ids(
        &mut self,
        residue: SmcraResidueId,
        label_comp_id: Option<String>,
        author_comp_id: Option<String>,
    ) -> std::result::Result<(), SmcraHierarchyError> {
        let residue = self
            .residues
            .get_mut(residue.index())
            .ok_or(SmcraHierarchyError::InvalidResidueId(residue))?;
        residue.label_comp_id = label_comp_id;
        residue.author_comp_id = author_comp_id;
        Ok(())
    }

    /// Returns checked mutable access to one chain's property map.
    pub fn chain_props_mut(
        &mut self,
        chain: SmcraChainId,
    ) -> std::result::Result<&mut PropMap, SmcraHierarchyError> {
        self.chains
            .get_mut(chain.index())
            .map(|chain| &mut chain.props)
            .ok_or(SmcraHierarchyError::InvalidChainId(chain))
    }

    /// Returns checked mutable access to one residue's property map.
    pub fn residue_props_mut(
        &mut self,
        residue: SmcraResidueId,
    ) -> std::result::Result<&mut PropMap, SmcraHierarchyError> {
        self.residues
            .get_mut(residue.index())
            .map(|residue| &mut residue.props)
            .ok_or(SmcraHierarchyError::InvalidResidueId(residue))
    }

    /// Returns checked mutable access to one atom site's property map.
    pub fn atom_site_props_mut(
        &mut self,
        site: SmcraAtomSiteId,
    ) -> std::result::Result<&mut PropMap, SmcraHierarchyError> {
        self.atom_sites
            .get_mut(site.index())
            .map(|site| &mut site.props)
            .ok_or(SmcraHierarchyError::InvalidAtomSiteId(site))
    }

    pub fn chains(&self) -> impl Iterator<Item = (SmcraChainId, &SmcraChain)> {
        self.chains.iter().map(|chain| (chain.id, chain))
    }

    pub fn residues(&self) -> impl Iterator<Item = (SmcraResidueId, &SmcraResidue)> {
        self.residues.iter().map(|residue| (residue.id, residue))
    }

    pub fn atom_sites(&self) -> impl Iterator<Item = (SmcraAtomSiteId, &SmcraAtomSite)> {
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

/// Fixed-width identifier spaces owned by [`SmcraHierarchy`].
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
pub enum SmcraHierarchyError {
    InvalidChainId(SmcraChainId),
    InvalidResidueId(SmcraResidueId),
    InvalidAtomSiteId(SmcraAtomSiteId),
    InvalidAtomId(AtomId),
    DuplicateAtomPlacement(AtomId),
    /// A new hierarchy node cannot be represented by the fixed-width ID for `kind`.
    IdentifierCapacityExceeded(SmcraIdKind),
}

impl fmt::Display for SmcraHierarchyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidChainId(id) => write!(f, "invalid chain id: {}", id.raw()),
            Self::InvalidResidueId(id) => write!(f, "invalid residue id: {}", id.raw()),
            Self::InvalidAtomSiteId(id) => write!(f, "invalid atom-site id: {}", id.raw()),
            Self::InvalidAtomId(id) => write!(f, "invalid hierarchy atom id: {id}"),
            Self::DuplicateAtomPlacement(id) => write!(f, "duplicate hierarchy placement for {id}"),
            Self::IdentifierCapacityExceeded(kind) => {
                write!(f, "{kind} identifier capacity exceeded")
            }
        }
    }
}

impl std::error::Error for SmcraHierarchyError {}

fn checked_hierarchy_id<T>(
    length: usize,
    kind: SmcraIdKind,
    construct: impl FnOnce(u32) -> T,
) -> std::result::Result<T, SmcraHierarchyError> {
    crate::core::checked_raw_id(length)
        .map(construct)
        .map_err(|_| SmcraHierarchyError::IdentifierCapacityExceeded(kind))
}

#[cfg(test)]
mod coordinate_independence_tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::structure::{ModelBuildError, ModelBuilder};
    use crate::topology::{MoleculeInstanceMetadata, TopologyBuilder};
    use crate::units::{Quantity, ANGSTROM};

    fn macro_with_unused_invalid_conformers(
        invalid_count: usize,
    ) -> (MacroMolecule, AtomId, ConformerId, ConformerId) {
        let mut graph = Molecule::new();
        let atom = graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");

        let mut valid = Conformer::new(ANGSTROM).unwrap();
        valid
            .set_position(atom, Quantity::new(Point3::new(1.0, 2.0, 3.0), ANGSTROM))
            .unwrap();
        let valid = graph.add_conformer(valid).unwrap();

        let mut first_invalid = None;
        for _ in 0..invalid_count {
            let mut invalid = Conformer::new(ANGSTROM).unwrap();
            invalid
                .set_position(
                    atom,
                    Quantity::new(Point3::new(f64::NAN, 0.0, 0.0), ANGSTROM),
                )
                .unwrap();
            let conformer = graph.add_conformer(invalid).unwrap();
            first_invalid.get_or_insert(conformer);
        }

        let mut hierarchy = SmcraHierarchy::new();
        let chain = hierarchy.add_chain("A", None).unwrap();
        let residue = hierarchy
            .add_residue(chain, "GLY", Some(1), None, None)
            .unwrap();
        hierarchy
            .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
            .unwrap();
        // Keep public checked assembly strict; this module-private fixture
        // exercises downstream behavior for an invalid unselected conformer.
        (
            MacroMolecule {
                molecule: graph,
                hierarchy,
            },
            atom,
            valid,
            first_invalid.expect("fixture requires an invalid conformer"),
        )
    }

    #[test]
    fn topology_and_model_construction_validate_only_selected_coordinate_state() {
        const UNUSED_CONFORMERS: usize = 128;
        let (molecule, atom, valid, invalid) =
            macro_with_unused_invalid_conformers(UNUSED_CONFORMERS);

        let static_report = molecule
            .validate_with_options(MacroValidateOptions {
                validate_coordinates: false,
            })
            .unwrap();
        assert_eq!(static_report.conformers_checked, 0);
        assert_eq!(static_report.coordinates_checked, 0);

        let mut topology_builder = TopologyBuilder::new();
        topology_builder
            .add_macro_molecule_instance(&molecule, MoleculeInstanceMetadata::default())
            .unwrap();
        let topology = topology_builder.build().unwrap();
        assert_eq!(topology.atom_count(), 1);
        assert_eq!(
            topology
                .definitions()
                .next()
                .unwrap()
                .1
                .graph()
                .conformers()
                .count(),
            0
        );

        let mut valid_builder = ModelBuilder::new();
        valid_builder.add_macro_molecule(&molecule, valid).unwrap();
        let model = valid_builder.build().unwrap();
        assert_eq!(
            model.positions().values().value()[0],
            Point3::new(1.0, 2.0, 3.0)
        );

        let mut invalid_builder = ModelBuilder::new();
        assert_eq!(
            invalid_builder.add_macro_molecule(&molecule, invalid),
            Err(ModelBuildError::NonFinitePosition { atom })
        );
        assert_eq!(
            molecule.validate(),
            Err(MacroValidateError::InvalidConformerAtom {
                conformer: invalid,
                atom,
            })
        );

        assert_eq!(
            molecule.as_molecule().conformers().count(),
            UNUSED_CONFORMERS + 1
        );
        assert!(molecule
            .as_molecule()
            .conformer(invalid)
            .unwrap()
            .position(atom)
            .unwrap()
            .value()
            .x
            .is_nan());
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn hierarchy_identifier_boundaries_are_structured_and_transactional() {
        let max_slot = usize::try_from(u64::from(u32::MAX)).expect("64-bit usize");
        let first_unsupported_slot =
            usize::try_from(u64::from(u32::MAX) + 1).expect("64-bit usize");
        assert_eq!(
            checked_hierarchy_id(max_slot, SmcraIdKind::Chain, SmcraChainId::new),
            Ok(SmcraChainId::new(u32::MAX))
        );
        assert_eq!(
            checked_hierarchy_id(
                first_unsupported_slot,
                SmcraIdKind::Residue,
                SmcraResidueId::new,
            ),
            Err(SmcraHierarchyError::IdentifierCapacityExceeded(
                SmcraIdKind::Residue
            ))
        );
        assert_eq!(
            checked_hierarchy_id(
                first_unsupported_slot,
                SmcraIdKind::AtomSite,
                SmcraAtomSiteId::new,
            ),
            Err(SmcraHierarchyError::IdentifierCapacityExceeded(
                SmcraIdKind::AtomSite
            ))
        );

        let mut hierarchy = SmcraHierarchy::new();
        let before = hierarchy.clone();
        assert_eq!(
            hierarchy.add_chain_at_slot("A".to_owned(), None, first_unsupported_slot,),
            Err(SmcraHierarchyError::IdentifierCapacityExceeded(
                SmcraIdKind::Chain
            ))
        );
        assert_eq!(hierarchy, before);
    }
}
