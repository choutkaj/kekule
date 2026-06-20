use std::collections::BTreeMap;
use std::fmt;

use crate::core::*;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MacroMolecule {
    pub mol: Molecule,
    pub hierarchy: BioHierarchy,
}

impl MacroMolecule {
    pub fn add_atom_site(
        &mut self,
        residue: ResidueId,
        atom: AtomId,
        metadata: AtomSiteMetadata,
    ) -> std::result::Result<AtomSiteId, BioHierarchyError> {
        self.mol
            .atom(atom)
            .map_err(|_| BioHierarchyError::InvalidAtomId(atom))?;
        self.hierarchy.add_atom_site(residue, atom, metadata)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ModelId(u32);

impl ModelId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ChainId(u32);

impl ChainId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ResidueId(u32);

impl ResidueId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomSiteId(u32);

impl AtomSiteId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioHierarchy {
    models: Vec<Model>,
    chains: Vec<Chain>,
    pub(crate) residues: Vec<Residue>,
    atom_sites: Vec<AtomSite>,
    atom_lookup: BTreeMap<AtomId, AtomSiteId>,
    pub props: PropMap,
}

impl BioHierarchy {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_model(&mut self, model_id: impl Into<String>) -> ModelId {
        let id = ModelId::new(self.models.len() as u32);
        self.models.push(Model {
            id,
            model_id: model_id.into(),
            chains: Vec::new(),
            props: PropMap::new(),
        });
        id
    }

    pub fn add_chain(
        &mut self,
        model: ModelId,
        label_id: impl Into<String>,
        author_id: Option<String>,
    ) -> std::result::Result<ChainId, BioHierarchyError> {
        self.model(model)?;
        let id = ChainId::new(self.chains.len() as u32);
        self.chains.push(Chain {
            id,
            model,
            label_id: label_id.into(),
            author_id,
            residues: Vec::new(),
            props: PropMap::new(),
        });
        self.models[model.index()].chains.push(id);
        Ok(id)
    }

    pub fn add_residue(
        &mut self,
        chain: ChainId,
        name: impl Into<String>,
        label_seq_id: Option<i32>,
        author_seq_id: Option<String>,
        insertion_code: Option<String>,
    ) -> std::result::Result<ResidueId, BioHierarchyError> {
        self.chain(chain)?;
        let name = name.into();
        let id = ResidueId::new(self.residues.len() as u32);
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

    pub fn add_atom_site(
        &mut self,
        residue: ResidueId,
        atom: AtomId,
        metadata: AtomSiteMetadata,
    ) -> std::result::Result<AtomSiteId, BioHierarchyError> {
        self.residue(residue)?;
        if self.atom_lookup.contains_key(&atom) {
            return Err(BioHierarchyError::DuplicateAtomPlacement(atom));
        }
        let id = AtomSiteId::new(self.atom_sites.len() as u32);
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

    pub fn model(&self, id: ModelId) -> std::result::Result<&Model, BioHierarchyError> {
        self.models
            .get(id.index())
            .ok_or(BioHierarchyError::InvalidModelId(id))
    }

    pub fn chain(&self, id: ChainId) -> std::result::Result<&Chain, BioHierarchyError> {
        self.chains
            .get(id.index())
            .ok_or(BioHierarchyError::InvalidChainId(id))
    }

    pub fn residue(&self, id: ResidueId) -> std::result::Result<&Residue, BioHierarchyError> {
        self.residues
            .get(id.index())
            .ok_or(BioHierarchyError::InvalidResidueId(id))
    }

    pub fn atom_site(&self, id: AtomSiteId) -> std::result::Result<&AtomSite, BioHierarchyError> {
        self.atom_sites
            .get(id.index())
            .ok_or(BioHierarchyError::InvalidAtomSiteId(id))
    }

    pub fn atom_site_for_atom(&self, atom: AtomId) -> Option<&AtomSite> {
        self.atom_lookup
            .get(&atom)
            .and_then(|id| self.atom_sites.get(id.index()))
    }

    pub fn models(&self) -> impl Iterator<Item = (ModelId, &Model)> {
        self.models.iter().map(|model| (model.id, model))
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
}

#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    pub id: ModelId,
    pub model_id: String,
    pub chains: Vec<ChainId>,
    pub props: PropMap,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Chain {
    pub id: ChainId,
    pub model: ModelId,
    pub label_id: String,
    pub author_id: Option<String>,
    pub residues: Vec<ResidueId>,
    pub props: PropMap,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub id: ResidueId,
    pub chain: ChainId,
    pub name: String,
    pub label_comp_id: Option<String>,
    pub author_comp_id: Option<String>,
    pub label_seq_id: Option<i32>,
    pub author_seq_id: Option<String>,
    pub insertion_code: Option<String>,
    pub atom_sites: Vec<AtomSiteId>,
    pub props: PropMap,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomSite {
    pub id: AtomSiteId,
    pub residue: ResidueId,
    pub atom: AtomId,
    pub metadata: AtomSiteMetadata,
    pub props: PropMap,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct AtomSiteMetadata {
    pub group_pdb: Option<String>,
    pub atom_site_id: Option<String>,
    pub type_symbol: Option<String>,
    pub label_asym_id: Option<String>,
    pub auth_asym_id: Option<String>,
    pub label_atom_id: Option<String>,
    pub auth_atom_id: Option<String>,
    pub label_alt_id: Option<String>,
    pub occupancy: Option<f64>,
    pub occupancy_raw: Option<String>,
    pub b_factor: Option<f64>,
    pub b_factor_raw: Option<String>,
    pub cartn_x_raw: Option<String>,
    pub cartn_y_raw: Option<String>,
    pub cartn_z_raw: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BioHierarchyError {
    InvalidModelId(ModelId),
    InvalidChainId(ChainId),
    InvalidResidueId(ResidueId),
    InvalidAtomSiteId(AtomSiteId),
    InvalidAtomId(AtomId),
    DuplicateAtomPlacement(AtomId),
}

impl fmt::Display for BioHierarchyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidModelId(id) => write!(f, "invalid model id: {}", id.raw()),
            Self::InvalidChainId(id) => write!(f, "invalid chain id: {}", id.raw()),
            Self::InvalidResidueId(id) => write!(f, "invalid residue id: {}", id.raw()),
            Self::InvalidAtomSiteId(id) => write!(f, "invalid atom-site id: {}", id.raw()),
            Self::InvalidAtomId(id) => write!(f, "invalid hierarchy atom id: {id}"),
            Self::DuplicateAtomPlacement(id) => write!(f, "duplicate hierarchy placement for {id}"),
        }
    }
}

impl std::error::Error for BioHierarchyError {}
