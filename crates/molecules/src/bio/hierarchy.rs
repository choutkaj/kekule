use std::collections::BTreeMap;
use std::fmt;

use crate::core::*;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MacroMolecule {
    graph: Molecule,
    hierarchy: SmcraHierarchy,
}

impl MacroMolecule {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn builder() -> MacroMoleculeBuilder {
        MacroMoleculeBuilder::new()
    }

    pub fn try_from_parts(
        graph: Molecule,
        hierarchy: SmcraHierarchy,
    ) -> std::result::Result<Self, MacroValidateError> {
        let molecule = Self { graph, hierarchy };
        molecule.validate()?;
        Ok(molecule)
    }

    pub fn graph(&self) -> &Molecule {
        &self.graph
    }

    pub fn hierarchy(&self) -> &SmcraHierarchy {
        &self.hierarchy
    }

    pub fn edit(&mut self) -> MacroMoleculeEditor<'_> {
        MacroMoleculeEditor {
            graph: self.graph.clone(),
            hierarchy: self.hierarchy.clone(),
            target: self,
        }
    }

    pub(crate) fn without_conformers(mut self) -> Self {
        self.graph = self.graph.without_conformers();
        self
    }

    pub fn models(&self) -> impl Iterator<Item = (SmcraModelId, &SmcraModel)> {
        self.hierarchy.models()
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

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MacroMoleculeBuilder {
    graph: Molecule,
    hierarchy: SmcraHierarchy,
}

impl MacroMoleculeBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn graph(&self) -> &Molecule {
        &self.graph
    }

    pub fn graph_mut(&mut self) -> &mut Molecule {
        &mut self.graph
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
        self.graph
            .atom(atom)
            .map_err(|_| SmcraHierarchyError::InvalidAtomId(atom))?;
        self.hierarchy.add_atom_site(residue, atom, metadata)
    }

    pub fn build(self) -> std::result::Result<MacroMolecule, MacroValidateError> {
        MacroMolecule::try_from_parts(self.graph, self.hierarchy)
    }
}

pub struct MacroMoleculeEditor<'a> {
    target: &'a mut MacroMolecule,
    graph: Molecule,
    hierarchy: SmcraHierarchy,
}

impl MacroMoleculeEditor<'_> {
    pub fn graph(&self) -> &Molecule {
        &self.graph
    }

    pub fn graph_mut(&mut self) -> &mut Molecule {
        &mut self.graph
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
        self.graph
            .atom(atom)
            .map_err(|_| SmcraHierarchyError::InvalidAtomId(atom))?;
        self.hierarchy.add_atom_site(residue, atom, metadata)
    }

    pub fn commit(self) -> std::result::Result<(), MacroValidateError> {
        let candidate = MacroMolecule::try_from_parts(self.graph, self.hierarchy)?;
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
    pub models_checked: usize,
    pub chains_checked: usize,
    pub residues_checked: usize,
    pub atom_sites_checked: usize,
    pub conformers_checked: usize,
    pub coordinates_checked: usize,
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MacroValidateError {
    InvalidChainModel {
        chain: SmcraChainId,
        model: SmcraModelId,
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
    InvalidAtomSiteOccupancy {
        site: SmcraAtomSiteId,
    },
    InvalidAtomSiteBFactor {
        site: SmcraAtomSiteId,
    },
    InvalidConformerAtom {
        conformer: ConformerId,
        atom: AtomId,
    },
}

impl fmt::Display for MacroValidateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidChainModel { chain, model } => {
                write!(
                    f,
                    "chain {} references invalid model {}",
                    chain.raw(),
                    model.raw()
                )
            }
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
            Self::InvalidAtomSiteOccupancy { site } => {
                write!(f, "atom-site {} has non-finite occupancy", site.raw())
            }
            Self::InvalidAtomSiteBFactor { site } => {
                write!(f, "atom-site {} has non-finite B-factor", site.raw())
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
    let mut report = MacroValidateReport {
        models_checked: molecule.hierarchy.models().count(),
        chains_checked: 0,
        residues_checked: 0,
        atom_sites_checked: 0,
        conformers_checked: 0,
        coordinates_checked: 0,
    };

    for (chain_id, chain) in molecule.hierarchy.chains() {
        molecule.hierarchy.model(chain.model).map_err(|_| {
            MacroValidateError::InvalidChainModel {
                chain: chain_id,
                model: chain.model,
            }
        })?;
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
            .graph
            .atom(site.atom)
            .map_err(|_| MacroValidateError::InvalidAtomSiteAtom {
                site: site_id,
                atom: site.atom,
            })?;
        if site
            .metadata
            .occupancy
            .is_some_and(|value| !value.is_finite())
        {
            return Err(MacroValidateError::InvalidAtomSiteOccupancy { site: site_id });
        }
        if site
            .metadata
            .b_factor
            .is_some_and(|value| !value.is_finite())
        {
            return Err(MacroValidateError::InvalidAtomSiteBFactor { site: site_id });
        }
        report.atom_sites_checked += 1;
    }
    for atom in molecule.graph.atom_ids() {
        if molecule.hierarchy.atom_site_for_atom(atom).is_none() {
            return Err(MacroValidateError::MissingAtomSiteForAtom { atom });
        }
    }
    if options.validate_coordinates {
        for (conformer_id, conformer) in molecule.graph.conformers() {
            report.conformers_checked += 1;
            for (atom, point) in conformer.positions() {
                let point = point.value();
                molecule.graph.atom(atom).map_err(|_| {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SmcraModelId(u32);

impl SmcraModelId {
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
pub struct SmcraChainId(u32);

impl SmcraChainId {
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
pub struct SmcraResidueId(u32);

impl SmcraResidueId {
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
pub struct SmcraAtomSiteId(u32);

impl SmcraAtomSiteId {
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
pub struct SmcraHierarchy {
    models: Vec<SmcraModel>,
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

    pub fn add_model(&mut self, model_id: impl Into<String>) -> SmcraModelId {
        let id = SmcraModelId::new(self.models.len() as u32);
        self.models.push(SmcraModel {
            id,
            model_id: model_id.into(),
            chains: Vec::new(),
            props: PropMap::new(),
        });
        id
    }

    pub fn add_chain(
        &mut self,
        model: SmcraModelId,
        label_id: impl Into<String>,
        author_id: Option<String>,
    ) -> std::result::Result<SmcraChainId, SmcraHierarchyError> {
        self.model(model)?;
        let id = SmcraChainId::new(self.chains.len() as u32);
        self.chains.push(SmcraChain {
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
        chain: SmcraChainId,
        name: impl Into<String>,
        label_seq_id: Option<i32>,
        author_seq_id: Option<String>,
        insertion_code: Option<String>,
    ) -> std::result::Result<SmcraResidueId, SmcraHierarchyError> {
        self.chain(chain)?;
        let name = name.into();
        let id = SmcraResidueId::new(self.residues.len() as u32);
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
        let id = SmcraAtomSiteId::new(self.atom_sites.len() as u32);
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

    pub fn model(&self, id: SmcraModelId) -> std::result::Result<&SmcraModel, SmcraHierarchyError> {
        self.models
            .get(id.index())
            .ok_or(SmcraHierarchyError::InvalidModelId(id))
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

    pub fn models(&self) -> impl Iterator<Item = (SmcraModelId, &SmcraModel)> {
        self.models.iter().map(|model| (model.id, model))
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
pub struct SmcraModel {
    pub(crate) id: SmcraModelId,
    pub(crate) model_id: String,
    pub(crate) chains: Vec<SmcraChainId>,
    pub(crate) props: PropMap,
}

impl SmcraModel {
    pub const fn id(&self) -> SmcraModelId {
        self.id
    }

    pub fn model_id(&self) -> &str {
        &self.model_id
    }

    pub fn chains(&self) -> &[SmcraChainId] {
        &self.chains
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SmcraChain {
    pub(crate) id: SmcraChainId,
    pub(crate) model: SmcraModelId,
    pub(crate) label_id: String,
    pub(crate) author_id: Option<String>,
    pub(crate) residues: Vec<SmcraResidueId>,
    pub(crate) props: PropMap,
}

impl SmcraChain {
    pub const fn id(&self) -> SmcraChainId {
        self.id
    }

    pub const fn model(&self) -> SmcraModelId {
        self.model
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

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmcraHierarchyError {
    InvalidModelId(SmcraModelId),
    InvalidChainId(SmcraChainId),
    InvalidResidueId(SmcraResidueId),
    InvalidAtomSiteId(SmcraAtomSiteId),
    InvalidAtomId(AtomId),
    DuplicateAtomPlacement(AtomId),
}

impl fmt::Display for SmcraHierarchyError {
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

impl std::error::Error for SmcraHierarchyError {}
