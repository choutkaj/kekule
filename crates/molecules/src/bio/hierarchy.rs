use std::collections::BTreeMap;
use std::fmt;

use crate::core::*;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MacroMolecule {
    graph: Molecule,
    hierarchy: BioHierarchy,
}

impl MacroMolecule {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_parts(graph: Molecule, hierarchy: BioHierarchy) -> Self {
        Self { graph, hierarchy }
    }

    pub fn graph(&self) -> &Molecule {
        &self.graph
    }

    pub fn graph_mut(&mut self) -> &mut Molecule {
        self.graph.invalidate_topology();
        &mut self.graph
    }

    pub fn hierarchy(&self) -> &BioHierarchy {
        &self.hierarchy
    }

    pub fn hierarchy_mut(&mut self) -> &mut BioHierarchy {
        &mut self.hierarchy
    }

    pub(crate) fn without_conformers(mut self) -> Self {
        self.graph = self.graph.without_conformers();
        self
    }

    pub fn models(&self) -> impl Iterator<Item = (ModelId, &Model)> {
        self.hierarchy.models()
    }

    pub fn chains(&self) -> impl Iterator<Item = (ChainId, &Chain)> {
        self.hierarchy.chains()
    }

    pub fn residues(&self) -> impl Iterator<Item = (ResidueId, &Residue)> {
        self.hierarchy.residues()
    }

    pub fn atom_sites(&self) -> impl Iterator<Item = (AtomSiteId, &AtomSite)> {
        self.hierarchy.atom_sites()
    }

    pub fn atom_site_for_atom(&self, atom: AtomId) -> Option<&AtomSite> {
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

    pub fn sanitize(&mut self) -> std::result::Result<MacroSanitizeReport, MacroSanitizeError> {
        self.sanitize_with_options(MacroSanitizeOptions::default())
    }

    pub fn sanitize_with_options(
        &mut self,
        options: MacroSanitizeOptions,
    ) -> std::result::Result<MacroSanitizeReport, MacroSanitizeError> {
        if !matches!(options.altloc_policy, AltLocPolicy::PreserveAll) {
            return Err(MacroSanitizeError::UnsupportedOption(
                "alternate-location selection is not implemented",
            ));
        }
        if options.normalize_elements
            || options.normalize_atom_site_metadata
            || options.recognize_standard_residues
        {
            return Err(MacroSanitizeError::UnsupportedOption(
                "element normalization, atom-site metadata normalization, or residue recognition is not implemented",
            ));
        }
        if options.assign_template_bonds
            || options.assign_polymer_bonds
            || options.detect_disulfides
            || !matches!(options.ligand_policy, LigandSanitizePolicy::LeaveRaw)
        {
            return Err(MacroSanitizeError::UnsupportedOption(
                "bond, disulfide, or ligand sanitization is not implemented",
            ));
        }
        let validation = if options.validate_first || options.validate_coordinates {
            Some(
                self.validate_with_options(MacroValidateOptions {
                    validate_coordinates: options.validate_coordinates,
                })
                .map_err(MacroSanitizeError::Validate)?,
            )
        } else {
            None
        };
        Ok(MacroSanitizeReport {
            validation,
            normalized_atom_sites: 0,
            recognized_residues: 0,
            assigned_bonds: 0,
        })
    }

    pub fn add_atom_site(
        &mut self,
        residue: ResidueId,
        atom: AtomId,
        metadata: AtomSiteMetadata,
    ) -> std::result::Result<AtomSiteId, BioHierarchyError> {
        self.graph
            .atom(atom)
            .map_err(|_| BioHierarchyError::InvalidAtomId(atom))?;
        self.hierarchy.add_atom_site(residue, atom, metadata)
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MacroValidateError {
    InvalidChainModel {
        chain: ChainId,
        model: ModelId,
    },
    InvalidResidueChain {
        residue: ResidueId,
        chain: ChainId,
    },
    InvalidResidueAtomSite {
        residue: ResidueId,
        site: AtomSiteId,
    },
    InvalidAtomSiteResidue {
        site: AtomSiteId,
        residue: ResidueId,
    },
    InvalidAtomSiteAtom {
        site: AtomSiteId,
        atom: AtomId,
    },
    InvalidAtomSiteOccupancy {
        site: AtomSiteId,
    },
    InvalidAtomSiteBFactor {
        site: AtomSiteId,
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MacroSanitizeOptions {
    pub validate_first: bool,
    pub normalize_elements: bool,
    pub normalize_atom_site_metadata: bool,
    pub validate_coordinates: bool,
    pub recognize_standard_residues: bool,
    pub assign_template_bonds: bool,
    pub assign_polymer_bonds: bool,
    pub detect_disulfides: bool,
    pub altloc_policy: AltLocPolicy,
    pub ligand_policy: LigandSanitizePolicy,
}

impl Default for MacroSanitizeOptions {
    fn default() -> Self {
        Self {
            validate_first: true,
            normalize_elements: false,
            normalize_atom_site_metadata: false,
            validate_coordinates: true,
            recognize_standard_residues: false,
            assign_template_bonds: false,
            assign_polymer_bonds: false,
            detect_disulfides: false,
            altloc_policy: AltLocPolicy::PreserveAll,
            ligand_policy: LigandSanitizePolicy::LeaveRaw,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AltLocPolicy {
    PreserveAll,
    SelectHighestOccupancy,
    SelectLabel(String),
    ErrorOnAltLoc,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LigandSanitizePolicy {
    LeaveRaw,
    SanitizeNonPolymerComponents,
    SanitizeAllDisconnectedComponents,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MacroSanitizeReport {
    pub validation: Option<MacroValidateReport>,
    pub normalized_atom_sites: usize,
    pub recognized_residues: usize,
    pub assigned_bonds: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MacroSanitizeError {
    Validate(MacroValidateError),
    UnsupportedOption(&'static str),
}

impl fmt::Display for MacroSanitizeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Validate(error) => write!(f, "{error}"),
            Self::UnsupportedOption(message) => write!(f, "{message}"),
        }
    }
}

impl std::error::Error for MacroSanitizeError {}

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
    if options.validate_coordinates {
        for (conformer_id, conformer) in molecule.graph.conformers() {
            report.conformers_checked += 1;
            for (atom, point) in conformer.positions() {
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
