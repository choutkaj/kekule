use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::bio::{Hierarchy, SmcraAtomSite, SmcraResidueId};
use crate::core::{AtomId, BondOrder, HydrogenDeclaration};
use crate::geometry::Point3;
use crate::io::mmcif_interpret::{MmcifEntityKind, MmcifInterpretationReport};
use crate::structure::Model;
use crate::topology::{
    InstanceAtomId, InstanceBondId, MoleculeDefinition, MoleculeInstance, MoleculeInstanceId,
};
use crate::units::ANGSTROM;

const MAX_COORDINATE_PRECISION: usize = 15;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifWriteOptions {
    pub data_block_name: String,
    pub coordinate_precision: usize,
}

impl Default for MmcifWriteOptions {
    fn default() -> Self {
        Self {
            data_block_name: "model".to_owned(),
            coordinate_precision: 3,
        }
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifWriteError {
    InvalidDataBlockName(String),
    CoordinatePrecisionTooLarge(usize),
    InvalidModel(String),
    InvalidHierarchy {
        message: String,
    },
    MissingEntityClassification(MoleculeInstanceId),
    DuplicateEntityClassification(MoleculeInstanceId),
    ConflictingEntityClassifications {
        molecule: MoleculeInstanceId,
        classifications: Vec<MmcifEntityKind>,
    },
    UnsupportedEntityClassification {
        molecule: MoleculeInstanceId,
        classification: String,
    },
    ConflictingAsymEntityIds {
        asym_id: String,
        entity_ids: Vec<String>,
    },
    ConflictingAsymEntityClassifications {
        asym_id: String,
        classifications: Vec<MmcifEntityKind>,
    },
    ConflictingSourceEntityClassifications {
        entity_id: String,
        classifications: Vec<MmcifEntityKind>,
    },
    UnknownClassifiedMolecule(MoleculeInstanceId),
    DuplicateAsymId(String),
    MissingAtomSite(InstanceAtomId),
    DuplicateAtomSite(InstanceAtomId),
    InconsistentAtomSite {
        atom: InstanceAtomId,
        field: &'static str,
    },
    InvalidGroupPdb {
        atom: InstanceAtomId,
        value: String,
    },
    DuplicateAtomIdentity(InstanceAtomId),
    MissingAtomProvenance(InstanceAtomId),
    DuplicateAtomProvenance(InstanceAtomId),
    UnknownAtomProvenance(InstanceAtomId),
    UnsupportedAtomField {
        atom: InstanceAtomId,
        field: &'static str,
    },
    FormalChargeOutOfRange {
        atom: InstanceAtomId,
        charge: i8,
    },
    UnsupportedStereo(MoleculeInstanceId),
    UnsupportedBondOrder {
        bond: InstanceBondId,
        order: BondOrder,
    },
    AmbiguousConnectionSelector(InstanceAtomId),
    UnsupportedTextValue {
        field: &'static str,
    },
}

impl fmt::Display for MmcifWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDataBlockName(name) => {
                write!(f, "invalid mmCIF data block name `{name}`")
            }
            Self::CoordinatePrecisionTooLarge(precision) => write!(
                f,
                "mmCIF coordinate precision {precision} exceeds the supported maximum of {MAX_COORDINATE_PRECISION}"
            ),
            Self::InvalidModel(message) => write!(f, "invalid molecular model: {message}"),
            Self::InvalidHierarchy { message } => {
                write!(f, "invalid topology hierarchy: {message}")
            }
            Self::MissingEntityClassification(molecule) => write!(
                f,
                "{molecule} has no explicit mmCIF entity classification"
            ),
            Self::DuplicateEntityClassification(molecule) => write!(
                f,
                "mmCIF entity semantics classify {molecule} more than once"
            ),
            Self::ConflictingEntityClassifications {
                molecule,
                classifications,
            } => write!(
                f,
                "{molecule} has conflicting mmCIF entity classifications {classifications:?}"
            ),
            Self::UnsupportedEntityClassification {
                molecule,
                classification,
            } => write!(
                f,
                "{molecule} has unsupported mmCIF entity classification `{classification}`"
            ),
            Self::ConflictingAsymEntityIds {
                asym_id,
                entity_ids,
            } => write!(
                f,
                "mmCIF structural instance `{asym_id}` has conflicting source entity IDs {entity_ids:?}"
            ),
            Self::ConflictingAsymEntityClassifications {
                asym_id,
                classifications,
            } => write!(
                f,
                "mmCIF structural instance `{asym_id}` has conflicting entity classifications {classifications:?}"
            ),
            Self::ConflictingSourceEntityClassifications {
                entity_id,
                classifications,
            } => write!(
                f,
                "source mmCIF entity `{entity_id}` has conflicting classifications {classifications:?}"
            ),
            Self::UnknownClassifiedMolecule(molecule) => write!(
                f,
                "mmCIF entity semantics reference unknown {molecule}"
            ),
            Self::DuplicateAsymId(id) => {
                write!(f, "duplicate mmCIF structural-instance ID `{id}`")
            }
            Self::MissingAtomSite(atom) => write!(f, "{atom} has no biomolecular atom site"),
            Self::DuplicateAtomSite(atom) => {
                write!(f, "{atom} appears in more than one biomolecular atom site")
            }
            Self::InconsistentAtomSite { atom, field } => {
                write!(f, "{atom} has inconsistent atom-site {field}")
            }
            Self::InvalidGroupPdb { atom, value } => write!(
                f,
                "{atom} has unsupported _atom_site.group_PDB value `{value}`"
            ),
            Self::DuplicateAtomIdentity(atom) => write!(
                f,
                "{atom} duplicates an mmCIF atom identity within one residue"
            ),
            Self::MissingAtomProvenance(atom) => {
                write!(f, "{atom} has no atom-level mmCIF source provenance")
            }
            Self::DuplicateAtomProvenance(atom) => {
                write!(f, "{atom} has duplicate atom-level mmCIF source provenance")
            }
            Self::UnknownAtomProvenance(atom) => {
                write!(f, "mmCIF source provenance references unknown {atom}")
            }
            Self::UnsupportedAtomField { atom, field } => {
                write!(f, "{atom} has unsupported atom field `{field}`")
            }
            Self::FormalChargeOutOfRange { atom, charge } => write!(
                f,
                "{atom} formal charge {charge} is outside the PDBx/mmCIF range -8..=8"
            ),
            Self::UnsupportedStereo(molecule) => write!(
                f,
                "{molecule} contains stereochemistry not represented by the foundational mmCIF writer"
            ),
            Self::UnsupportedBondOrder { bond, order } => {
                write!(f, "{bond} has unsupported mmCIF bond order {order:?}")
            }
            Self::AmbiguousConnectionSelector(atom) => write!(
                f,
                "{atom} cannot be selected unambiguously by an mmCIF struct_conn partner"
            ),
            Self::UnsupportedTextValue { field } => write!(
                f,
                "{field} contains a text value that cannot be emitted as a single mmCIF token"
            ),
        }
    }
}

impl std::error::Error for MmcifWriteError {}

/// Explicit mmCIF entity semantics for generic molecule instances.
///
/// Generic [`Model`] and topology state intentionally do not classify
/// molecules as polymers, branched entities, non-polymers, or water.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MmcifEntityClassifications {
    kinds: BTreeMap<MoleculeInstanceId, MmcifEntityKind>,
}

impl MmcifEntityClassifications {
    pub const fn new() -> Self {
        Self {
            kinds: BTreeMap::new(),
        }
    }

    /// Assigns one explicit mmCIF entity kind to a molecule instance.
    pub fn insert(
        &mut self,
        molecule: MoleculeInstanceId,
        kind: MmcifEntityKind,
    ) -> Result<(), MmcifWriteError> {
        if self.kinds.contains_key(&molecule) {
            return Err(MmcifWriteError::DuplicateEntityClassification(molecule));
        }
        self.kinds.insert(molecule, kind);
        Ok(())
    }

    /// Returns the assigned kind for one molecule instance.
    pub fn get(&self, molecule: MoleculeInstanceId) -> Option<&MmcifEntityKind> {
        self.kinds.get(&molecule)
    }

    /// Iterates over explicitly classified molecule instances in ID order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &MmcifEntityKind)> {
        self.kinds.iter().map(|(&molecule, kind)| (molecule, kind))
    }

    /// Returns the number of explicitly classified molecule instances.
    pub fn len(&self) -> usize {
        self.kinds.len()
    }

    /// Returns whether no molecule instance has been classified.
    pub fn is_empty(&self) -> bool {
        self.kinds.is_empty()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EntityKind {
    Polymer,
    Branched,
    NonPolymer,
    Water,
}

impl EntityKind {
    const fn as_mmcif(self) -> &'static str {
        match self {
            Self::Polymer => "polymer",
            Self::Branched => "branched",
            Self::NonPolymer => "non-polymer",
            Self::Water => "water",
        }
    }

    const fn default_group_pdb(self) -> &'static str {
        match self {
            Self::Polymer => "ATOM",
            Self::Branched | Self::NonPolymer | Self::Water => "HETATM",
        }
    }

    fn as_source(self) -> MmcifEntityKind {
        match self {
            Self::Polymer => MmcifEntityKind::Polymer,
            Self::Branched => MmcifEntityKind::Branched,
            Self::NonPolymer => MmcifEntityKind::NonPolymer,
            Self::Water => MmcifEntityKind::Water,
        }
    }
}

#[derive(Debug, Clone)]
struct EntityRow {
    id: String,
    kind: EntityKind,
}

#[derive(Debug, Clone)]
struct AsymRow {
    id: String,
    entity_id: String,
}

#[derive(Debug, Clone)]
struct AtomEntityAssignment {
    entity_id: String,
    asym_id: String,
    kind: EntityKind,
}

#[derive(Debug)]
struct EntityPlan {
    entities: Vec<EntityRow>,
    asyms: Vec<AsymRow>,
    atoms: BTreeMap<InstanceAtomId, AtomEntityAssignment>,
}

#[derive(Debug, Clone)]
struct AtomRow {
    atom: InstanceAtomId,
    residue: Option<SmcraResidueId>,
    entity_id: String,
    asym_id: String,
    group_pdb: String,
    type_symbol: String,
    label_atom_id: String,
    label_alt_id: Option<String>,
    label_comp_id: String,
    label_seq_id: Option<i32>,
    insertion_code: Option<String>,
    position: Point3,
    occupancy: Option<f64>,
    b_factor: Option<f64>,
    formal_charge: i8,
    auth_seq_id: Option<String>,
    auth_comp_id: String,
    auth_asym_id: String,
    auth_atom_id: String,
}

#[derive(Debug, Clone)]
struct ConnectionRow {
    left: InstanceAtomId,
    right: InstanceAtomId,
    order: BondOrder,
}

#[derive(Debug)]
struct PreparedModel {
    entities: Vec<EntityRow>,
    asyms: Vec<AsymRow>,
    atoms: Vec<AtomRow>,
    connections: Vec<ConnectionRow>,
}

pub fn write_mmcif_model(
    model: &Model,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    validate_options(&options)?;
    let classifications = normalize_entity_classifications(model, std::iter::empty())?;
    let plan = generic_entity_plan(model, &classifications)?;
    let prepared = prepare_model(model, plan)?;
    render_model(&prepared, &options)
}

pub fn write_mmcif_model_with_classifications(
    model: &Model,
    classifications: &MmcifEntityClassifications,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    validate_options(&options)?;
    let classifications = normalize_entity_classifications(
        model,
        classifications
            .iter()
            .map(|(molecule, kind)| (molecule, vec![kind.clone()])),
    )?;
    let plan = generic_entity_plan(model, &classifications)?;
    let prepared = prepare_model(model, plan)?;
    render_model(&prepared, &options)
}

pub fn write_mmcif_model_with_report(
    model: &Model,
    report: &MmcifInterpretationReport,
    options: MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    validate_options(&options)?;
    let plan = if report
        .instances()
        .iter()
        .all(|instance| instance.atoms().is_empty())
    {
        let classifications = normalize_entity_classifications(
            model,
            report
                .instances()
                .iter()
                .map(|instance| (instance.molecule(), instance.entity_kinds().to_vec())),
        )?;
        generic_entity_plan(model, &classifications)?
    } else {
        report_entity_plan(model, report)?
    };
    let prepared = prepare_model(model, plan)?;
    render_model(&prepared, &options)
}

fn validate_options(options: &MmcifWriteOptions) -> Result<(), MmcifWriteError> {
    if options.data_block_name.is_empty()
        || !options
            .data_block_name
            .chars()
            .all(|character| character.is_ascii_alphanumeric() || "_-.".contains(character))
    {
        return Err(MmcifWriteError::InvalidDataBlockName(
            options.data_block_name.clone(),
        ));
    }
    if options.coordinate_precision > MAX_COORDINATE_PRECISION {
        return Err(MmcifWriteError::CoordinatePrecisionTooLarge(
            options.coordinate_precision,
        ));
    }
    Ok(())
}

fn hierarchy_asym_ids(model: &Model) -> Result<BTreeSet<String>, MmcifWriteError> {
    let mut reserved_asym_ids = BTreeSet::new();
    for (_, chain) in model.topology().hierarchy().chains() {
        if chain.label_id().is_empty() {
            return Err(MmcifWriteError::InvalidModel(
                "topology hierarchy chain label ID is empty".to_owned(),
            ));
        }
        if !reserved_asym_ids.insert(chain.label_id().to_owned()) {
            return Err(MmcifWriteError::DuplicateAsymId(
                chain.label_id().to_owned(),
            ));
        }
    }
    Ok(reserved_asym_ids)
}

fn generic_entity_plan(
    model: &Model,
    classifications: &BTreeMap<MoleculeInstanceId, EntityKind>,
) -> Result<EntityPlan, MmcifWriteError> {
    let hierarchy = model.topology().hierarchy();
    let mut reserved_asym_ids = hierarchy_asym_ids(model)?;
    let mut entities = Vec::new();
    let mut asyms = Vec::new();
    let mut atoms = BTreeMap::new();

    for (_, chain) in hierarchy.chains() {
        let sites = chain
            .residues()
            .iter()
            .flat_map(|residue| {
                hierarchy
                    .residue(*residue)
                    .into_iter()
                    .flat_map(|residue| residue.atom_sites().iter().copied())
            })
            .filter_map(|site| hierarchy.atom_site(site).ok())
            .collect::<Vec<_>>();
        if sites.is_empty() {
            continue;
        }
        let source_kinds = sites
            .iter()
            .filter_map(|site| classifications.get(&site.atom().molecule()).copied())
            .map(EntityKind::as_source)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let kind = match source_kinds.as_slice() {
            [kind] => entity_kind_from_source(sites[0].atom().molecule(), kind)?,
            _ => {
                return Err(MmcifWriteError::ConflictingAsymEntityClassifications {
                    asym_id: chain.label_id().to_owned(),
                    classifications: source_kinds,
                });
            }
        };
        let entity_id = (entities.len() + 1).to_string();
        entities.push(EntityRow {
            id: entity_id.clone(),
            kind,
        });
        asyms.push(AsymRow {
            id: chain.label_id().to_owned(),
            entity_id: entity_id.clone(),
        });
        for site in sites {
            let atom = site.atom();
            if atoms
                .insert(
                    atom,
                    AtomEntityAssignment {
                        entity_id: entity_id.clone(),
                        asym_id: chain.label_id().to_owned(),
                        kind,
                    },
                )
                .is_some()
            {
                return Err(MmcifWriteError::DuplicateAtomSite(atom));
            }
        }
    }

    for (id, molecule) in model.topology().instances() {
        let definition = model
            .topology()
            .definition_for_instance(id)
            .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?;
        let qualified_atoms = definition
            .molecule()
            .atoms()
            .map(|(atom, _)| molecule.qualify_atom(atom))
            .collect::<Vec<_>>();
        let assigned = qualified_atoms
            .iter()
            .filter(|atom| atoms.contains_key(atom))
            .count();
        if assigned != 0 {
            if assigned != qualified_atoms.len() {
                let missing = qualified_atoms
                    .into_iter()
                    .find(|atom| !atoms.contains_key(atom))
                    .expect("partially assigned molecule has a missing atom site");
                return Err(MmcifWriteError::MissingAtomSite(missing));
            }
            continue;
        }
        let base = format!("M{}", one_based_serial(id.raw()));
        let mut candidate = base.clone();
        let mut suffix = 2usize;
        while reserved_asym_ids.contains(&candidate) {
            candidate = format!("{base}_{suffix}");
            suffix += 1;
        }
        reserved_asym_ids.insert(candidate.clone());
        let kind = entity_kind(molecule, classifications)?;
        let entity_id = (entities.len() + 1).to_string();
        entities.push(EntityRow {
            id: entity_id.clone(),
            kind,
        });
        asyms.push(AsymRow {
            id: candidate.clone(),
            entity_id: entity_id.clone(),
        });
        for atom in qualified_atoms {
            atoms.insert(
                atom,
                AtomEntityAssignment {
                    entity_id: entity_id.clone(),
                    asym_id: candidate.clone(),
                    kind,
                },
            );
        }
    }

    Ok(EntityPlan {
        entities,
        asyms,
        atoms,
    })
}

fn report_entity_plan(
    model: &Model,
    report: &MmcifInterpretationReport,
) -> Result<EntityPlan, MmcifWriteError> {
    let _ = hierarchy_asym_ids(model)?;
    let mut seen_instances = BTreeSet::new();
    let mut provenance = BTreeMap::new();
    let mut reserved_entity_ids = BTreeSet::new();
    for instance in report.instances() {
        if model.topology().instance(instance.molecule()).is_err() {
            return Err(MmcifWriteError::UnknownClassifiedMolecule(
                instance.molecule(),
            ));
        }
        if !seen_instances.insert(instance.molecule()) {
            return Err(MmcifWriteError::DuplicateEntityClassification(
                instance.molecule(),
            ));
        }
        for atom in instance.atoms() {
            if atom.atom().molecule() != instance.molecule()
                || model.topology().atom(atom.atom()).is_err()
            {
                return Err(MmcifWriteError::UnknownAtomProvenance(atom.atom()));
            }
            if let Some(entity_id) = atom.entity_id() {
                reserved_entity_ids.insert(entity_id.to_owned());
            }
            if provenance.insert(atom.atom(), atom).is_some() {
                return Err(MmcifWriteError::DuplicateAtomProvenance(atom.atom()));
            }
        }
    }

    let hierarchy = model.topology().hierarchy();
    let mut entities = Vec::new();
    let mut entity_kinds = BTreeMap::<String, EntityKind>::new();
    let mut asyms = Vec::new();
    let mut atoms = BTreeMap::new();
    let mut generated_serial = 1usize;
    for (_, chain) in hierarchy.chains() {
        let sites = chain
            .residues()
            .iter()
            .flat_map(|residue| {
                hierarchy
                    .residue(*residue)
                    .into_iter()
                    .flat_map(|residue| residue.atom_sites().iter().copied())
            })
            .filter_map(|site| hierarchy.atom_site(site).ok())
            .collect::<Vec<_>>();
        if sites.is_empty() {
            continue;
        }
        let source_atoms = sites
            .iter()
            .map(|site| {
                provenance
                    .get(&site.atom())
                    .copied()
                    .ok_or(MmcifWriteError::MissingAtomProvenance(site.atom()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        for atom in &source_atoms {
            if atom.asym_id() != chain.label_id() {
                return Err(MmcifWriteError::InconsistentAtomSite {
                    atom: atom.atom(),
                    field: "source asym_id",
                });
            }
        }
        let source_entity_ids = source_atoms
            .iter()
            .filter_map(|atom| atom.entity_id().map(str::to_owned))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        if source_entity_ids.len() > 1 {
            return Err(MmcifWriteError::ConflictingAsymEntityIds {
                asym_id: chain.label_id().to_owned(),
                entity_ids: source_entity_ids,
            });
        }
        let source_kinds = source_atoms
            .iter()
            .map(|atom| atom.entity_kind.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let kind = match source_kinds.as_slice() {
            [kind] => entity_kind_from_source(source_atoms[0].atom().molecule(), kind)?,
            _ => {
                return Err(MmcifWriteError::ConflictingAsymEntityClassifications {
                    asym_id: chain.label_id().to_owned(),
                    classifications: source_kinds,
                });
            }
        };
        let entity_id = if let Some(entity_id) = source_entity_ids.first() {
            entity_id.clone()
        } else {
            loop {
                let candidate = format!("K{generated_serial}");
                generated_serial += 1;
                if reserved_entity_ids.insert(candidate.clone()) {
                    break candidate;
                }
            }
        };
        if let Some(existing) = entity_kinds.get(&entity_id) {
            if *existing != kind {
                return Err(MmcifWriteError::ConflictingSourceEntityClassifications {
                    entity_id,
                    classifications: vec![existing.as_source(), kind.as_source()],
                });
            }
        } else {
            entity_kinds.insert(entity_id.clone(), kind);
            entities.push(EntityRow {
                id: entity_id.clone(),
                kind,
            });
        }
        asyms.push(AsymRow {
            id: chain.label_id().to_owned(),
            entity_id: entity_id.clone(),
        });
        for site in sites {
            let atom = site.atom();
            if atoms
                .insert(
                    atom,
                    AtomEntityAssignment {
                        entity_id: entity_id.clone(),
                        asym_id: chain.label_id().to_owned(),
                        kind,
                    },
                )
                .is_some()
            {
                return Err(MmcifWriteError::DuplicateAtomSite(atom));
            }
        }
    }

    for atom in model.topology().atom_ids() {
        if !atoms.contains_key(atom) {
            return Err(MmcifWriteError::MissingAtomProvenance(*atom));
        }
    }
    for atom in provenance.keys() {
        if !atoms.contains_key(atom) {
            return Err(MmcifWriteError::UnknownAtomProvenance(*atom));
        }
    }

    Ok(EntityPlan {
        entities,
        asyms,
        atoms,
    })
}

fn prepare_model(model: &Model, plan: EntityPlan) -> Result<PreparedModel, MmcifWriteError> {
    let EntityPlan {
        entities,
        asyms,
        atoms: assignments,
    } = plan;

    let mut atoms = Vec::new();
    for (id, molecule) in model.topology().instances() {
        let definition = model
            .topology()
            .definition_for_instance(id)
            .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?;
        validate_graph_chemistry(molecule, definition)?;
        if model
            .topology()
            .molecule(id)
            .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?
            .atom_sites()
            .next()
            .is_some()
        {
            collect_macro_rows(
                model,
                molecule,
                definition,
                model.topology().hierarchy(),
                &assignments,
                &mut atoms,
            )?;
        } else {
            collect_small_rows(model, molecule, definition, &assignments, &mut atoms)?;
        }
    }

    let hierarchy_order = model
        .topology()
        .hierarchy()
        .chains()
        .flat_map(|(_, chain)| chain.residues().iter().copied())
        .filter_map(|residue| model.topology().hierarchy().residue(residue).ok())
        .flat_map(|residue| residue.atom_sites().iter().copied())
        .filter_map(|site| model.topology().hierarchy().atom_site(site).ok())
        .enumerate()
        .map(|(index, site)| (site.atom(), index))
        .collect::<BTreeMap<_, _>>();
    atoms.sort_by_key(|row| {
        hierarchy_order
            .get(&row.atom)
            .copied()
            .unwrap_or(usize::MAX)
    });

    let atom_indexes = atoms
        .iter()
        .enumerate()
        .map(|(index, row)| (row.atom, index))
        .collect::<BTreeMap<_, _>>();
    validate_atom_identities(&atoms)?;
    let mut connections = Vec::new();
    for (bond_id, bond) in model.topology().bonds() {
        let order = supported_bond_order(bond_id, bond.order)?;
        let molecule = model
            .topology()
            .instance(bond_id.molecule())
            .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?;
        let left = molecule.qualify_atom(bond.a());
        let right = molecule.qualify_atom(bond.b());
        validate_connection_selector(left, &atoms, &atom_indexes)?;
        validate_connection_selector(right, &atoms, &atom_indexes)?;
        connections.push(ConnectionRow { left, right, order });
    }
    Ok(PreparedModel {
        entities,
        asyms,
        atoms,
        connections,
    })
}

fn entity_kind(
    molecule: &MoleculeInstance,
    classifications: &BTreeMap<MoleculeInstanceId, EntityKind>,
) -> Result<EntityKind, MmcifWriteError> {
    classifications
        .get(&molecule.id())
        .copied()
        .ok_or(MmcifWriteError::MissingEntityClassification(molecule.id()))
}

fn entity_kind_from_source(
    molecule: MoleculeInstanceId,
    kind: &MmcifEntityKind,
) -> Result<EntityKind, MmcifWriteError> {
    match kind {
        MmcifEntityKind::Polymer => Ok(EntityKind::Polymer),
        MmcifEntityKind::Branched => Ok(EntityKind::Branched),
        MmcifEntityKind::NonPolymer => Ok(EntityKind::NonPolymer),
        MmcifEntityKind::Water => Ok(EntityKind::Water),
        MmcifEntityKind::Other(classification) => {
            Err(MmcifWriteError::UnsupportedEntityClassification {
                molecule,
                classification: classification.clone(),
            })
        }
    }
}

fn normalize_entity_classifications(
    model: &Model,
    entries: impl IntoIterator<Item = (MoleculeInstanceId, Vec<MmcifEntityKind>)>,
) -> Result<BTreeMap<MoleculeInstanceId, EntityKind>, MmcifWriteError> {
    let mut classifications = BTreeMap::new();
    for (molecule, source_kinds) in entries {
        if model.topology().instance(molecule).is_err() {
            return Err(MmcifWriteError::UnknownClassifiedMolecule(molecule));
        }
        let kinds = source_kinds
            .into_iter()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let kind = match kinds.as_slice() {
            [kind] => entity_kind_from_source(molecule, kind)?,
            [] => return Err(MmcifWriteError::MissingEntityClassification(molecule)),
            _ => {
                return Err(MmcifWriteError::ConflictingEntityClassifications {
                    molecule,
                    classifications: kinds,
                });
            }
        };
        if classifications.insert(molecule, kind).is_some() {
            return Err(MmcifWriteError::DuplicateEntityClassification(molecule));
        }
    }
    for (molecule, _) in model.topology().instances() {
        if !classifications.contains_key(&molecule) {
            return Err(MmcifWriteError::MissingEntityClassification(molecule));
        }
    }
    Ok(classifications)
}

fn validate_graph_chemistry(
    molecule: &MoleculeInstance,
    definition: &MoleculeDefinition,
) -> Result<(), MmcifWriteError> {
    let molecule_definition = definition.molecule();
    if molecule_definition.stereo_elements().next().is_some()
        || molecule_definition.stereo_groups().next().is_some()
    {
        return Err(MmcifWriteError::UnsupportedStereo(molecule.id()));
    }
    for (atom_id, atom) in molecule_definition.atoms() {
        let atom_id = molecule.qualify_atom(atom_id);
        if atom.isotope.is_some() {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "isotope",
            });
        }
        if atom.radical.is_some() {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "radical",
            });
        }
        if atom.hydrogens != HydrogenDeclaration::default() {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "hydrogens",
            });
        }
        if atom.atom_map.is_some() {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "atom_map",
            });
        }
        if !(-8..=8).contains(&atom.formal_charge) {
            return Err(MmcifWriteError::FormalChargeOutOfRange {
                atom: atom_id,
                charge: atom.formal_charge,
            });
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn collect_macro_rows(
    model: &Model,
    molecule: &MoleculeInstance,
    definition: &MoleculeDefinition,
    hierarchy: &Hierarchy,
    assignments: &BTreeMap<InstanceAtomId, AtomEntityAssignment>,
    rows: &mut Vec<AtomRow>,
) -> Result<(), MmcifWriteError> {
    let mut sites = BTreeMap::<InstanceAtomId, &SmcraAtomSite>::new();
    for (_, site) in hierarchy.atom_sites() {
        if site.atom.molecule() == molecule.id() && sites.insert(site.atom, site).is_some() {
            return Err(MmcifWriteError::DuplicateAtomSite(site.atom));
        }
    }
    for (atom_id, atom) in definition.molecule().atoms() {
        let qualified = molecule.qualify_atom(atom_id);
        let site = sites
            .get(&qualified)
            .copied()
            .ok_or(MmcifWriteError::MissingAtomSite(qualified))?;
        let residue = hierarchy
            .residue(site.residue)
            .map_err(|error| invalid_hierarchy(error.to_string()))?;
        let chain = hierarchy
            .chain(residue.chain)
            .map_err(|error| invalid_hierarchy(error.to_string()))?;
        let assignment = assignments
            .get(&qualified)
            .ok_or(MmcifWriteError::MissingAtomSite(qualified))?;
        if assignment.asym_id != chain.label_id {
            return Err(MmcifWriteError::InconsistentAtomSite {
                atom: qualified,
                field: "entity/asymmetry assignment",
            });
        }
        if site
            .metadata
            .label_asym_id
            .as_deref()
            .is_some_and(|value| value != chain.label_id)
        {
            return Err(MmcifWriteError::InconsistentAtomSite {
                atom: qualified,
                field: "label_asym_id",
            });
        }
        if site
            .metadata
            .type_symbol
            .as_deref()
            .is_some_and(|value| !value.eq_ignore_ascii_case(atom.element.symbol()))
        {
            return Err(MmcifWriteError::InconsistentAtomSite {
                atom: qualified,
                field: "type_symbol",
            });
        }
        let group_pdb = normalized_group_pdb(qualified, None, assignment.kind.default_group_pdb())?;
        let label_atom_id = site
            .metadata
            .label_atom_id
            .as_ref()
            .or(site.metadata.auth_atom_id.as_ref())
            .cloned()
            .unwrap_or_else(|| generated_atom_name(atom.element.symbol(), atom_id));
        let label_comp_id = residue
            .label_comp_id
            .clone()
            .unwrap_or_else(|| residue.name.clone());
        rows.push(AtomRow {
            atom: qualified,
            residue: Some(site.residue()),
            entity_id: assignment.entity_id.clone(),
            asym_id: chain.label_id.clone(),
            group_pdb,
            type_symbol: atom.element.symbol().to_owned(),
            label_atom_id: label_atom_id.clone(),
            label_alt_id: None,
            label_comp_id: label_comp_id.clone(),
            label_seq_id: residue.label_seq_id,
            insertion_code: residue.insertion_code.clone(),
            position: model
                .position(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?
                .value_in(ANGSTROM)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            occupancy: model
                .occupancy(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            b_factor: model
                .b_factor(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?
                .map(|value| value.value_in(crate::units::SQUARE_ANGSTROM))
                .transpose()
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            formal_charge: atom.formal_charge,
            auth_seq_id: residue
                .author_seq_id
                .clone()
                .or_else(|| residue.label_seq_id.map(|value| value.to_string())),
            auth_comp_id: residue.author_comp_id.clone().unwrap_or(label_comp_id),
            auth_asym_id: site
                .metadata
                .auth_asym_id
                .clone()
                .or_else(|| chain.author_id.clone())
                .unwrap_or_else(|| chain.label_id.clone()),
            auth_atom_id: site.metadata.auth_atom_id.clone().unwrap_or(label_atom_id),
        });
    }
    Ok(())
}

fn collect_small_rows(
    model: &Model,
    molecule: &MoleculeInstance,
    definition: &MoleculeDefinition,
    assignments: &BTreeMap<InstanceAtomId, AtomEntityAssignment>,
    rows: &mut Vec<AtomRow>,
) -> Result<(), MmcifWriteError> {
    for (atom_id, atom) in definition.molecule().atoms() {
        let qualified = molecule.qualify_atom(atom_id);
        let assignment = assignments
            .get(&qualified)
            .ok_or(MmcifWriteError::MissingAtomSite(qualified))?;
        let component_id = if assignment.kind == EntityKind::Water {
            "HOH"
        } else {
            "MOL"
        };
        let atom_name = generated_atom_name(atom.element.symbol(), atom_id);
        rows.push(AtomRow {
            atom: qualified,
            residue: None,
            entity_id: assignment.entity_id.clone(),
            asym_id: assignment.asym_id.clone(),
            group_pdb: normalized_group_pdb(qualified, None, assignment.kind.default_group_pdb())?,
            type_symbol: atom.element.symbol().to_owned(),
            label_atom_id: atom_name.clone(),
            label_alt_id: None,
            label_comp_id: component_id.to_owned(),
            label_seq_id: None,
            insertion_code: None,
            position: model
                .position(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?
                .value_in(ANGSTROM)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            occupancy: model
                .occupancy(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            b_factor: model
                .b_factor(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?
                .map(|value| value.value_in(crate::units::SQUARE_ANGSTROM))
                .transpose()
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            formal_charge: atom.formal_charge,
            auth_seq_id: None,
            auth_comp_id: component_id.to_owned(),
            auth_asym_id: assignment.asym_id.clone(),
            auth_atom_id: atom_name,
        });
    }
    Ok(())
}

fn normalized_group_pdb(
    atom: InstanceAtomId,
    value: Option<&str>,
    default: &str,
) -> Result<String, MmcifWriteError> {
    let value = value.unwrap_or(default);
    if value.eq_ignore_ascii_case("ATOM") {
        Ok("ATOM".to_owned())
    } else if value.eq_ignore_ascii_case("HETATM") {
        Ok("HETATM".to_owned())
    } else {
        Err(MmcifWriteError::InvalidGroupPdb {
            atom,
            value: value.to_owned(),
        })
    }
}

fn generated_atom_name(symbol: &str, atom: AtomId) -> String {
    format!("{symbol}{}", one_based_serial(atom.raw()))
}

fn one_based_serial(raw: u32) -> u64 {
    u64::from(raw) + 1
}

fn supported_bond_order(
    bond: InstanceBondId,
    order: BondOrder,
) -> Result<BondOrder, MmcifWriteError> {
    match order {
        BondOrder::Single | BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple => {
            Ok(order)
        }
        BondOrder::Zero | BondOrder::Dative => {
            Err(MmcifWriteError::UnsupportedBondOrder { bond, order })
        }
    }
}

fn validate_connection_selector(
    atom: InstanceAtomId,
    rows: &[AtomRow],
    indexes: &BTreeMap<InstanceAtomId, usize>,
) -> Result<(), MmcifWriteError> {
    let row =
        rows.get(*indexes.get(&atom).ok_or_else(|| {
            MmcifWriteError::InvalidModel(format!("missing atom row for {atom}"))
        })?)
        .expect("atom row index is valid");
    let matches = rows
        .iter()
        .filter(|candidate| {
            candidate.asym_id == row.asym_id
                && candidate.label_atom_id == row.label_atom_id
                && row
                    .label_seq_id
                    .is_none_or(|sequence| candidate.label_seq_id == Some(sequence))
        })
        .count();
    if matches != 1 {
        return Err(MmcifWriteError::AmbiguousConnectionSelector(atom));
    }
    Ok(())
}

fn validate_atom_identities(rows: &[AtomRow]) -> Result<(), MmcifWriteError> {
    let mut identities = BTreeSet::new();
    for row in rows {
        let residue = if let Some(sequence) = row.label_seq_id {
            format!(
                "label:{sequence}:{}",
                row.insertion_code.as_deref().unwrap_or("")
            )
        } else if let Some(sequence) = &row.auth_seq_id {
            format!(
                "auth:{sequence}:{}",
                row.insertion_code.as_deref().unwrap_or("")
            )
        } else {
            format!("unsequenced:{:?}", row.residue)
        };
        if !identities.insert((row.asym_id.clone(), residue, row.label_atom_id.clone())) {
            return Err(MmcifWriteError::DuplicateAtomIdentity(row.atom));
        }
    }
    Ok(())
}

fn invalid_hierarchy(message: impl Into<String>) -> MmcifWriteError {
    MmcifWriteError::InvalidHierarchy {
        message: message.into(),
    }
}

fn render_model(
    model: &PreparedModel,
    options: &MmcifWriteOptions,
) -> Result<String, MmcifWriteError> {
    let mut output = String::with_capacity(model.atoms.len().saturating_mul(160));
    output.push_str("data_");
    output.push_str(&options.data_block_name);
    output.push_str("\n#\n");

    write_loop_header(&mut output, &["_entity.id", "_entity.type"]);
    for entity in &model.entities {
        write_row(
            &mut output,
            vec![
                cif_value(&entity.id, "_entity.id")?,
                entity.kind.as_mmcif().to_owned(),
            ],
        );
    }
    output.push_str("#\n");

    write_loop_header(&mut output, &["_struct_asym.id", "_struct_asym.entity_id"]);
    for asym in &model.asyms {
        write_row(
            &mut output,
            vec![
                cif_value(&asym.id, "_struct_asym.id")?,
                cif_value(&asym.entity_id, "_struct_asym.entity_id")?,
            ],
        );
    }
    output.push_str("#\n");

    const ATOM_TAGS: &[&str] = &[
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_formal_charge",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    ];
    write_loop_header(&mut output, ATOM_TAGS);
    for (serial, atom) in (1u64..).zip(model.atoms.iter()) {
        write_row(
            &mut output,
            vec![
                atom.group_pdb.clone(),
                serial.to_string(),
                cif_value(&atom.type_symbol, "_atom_site.type_symbol")?,
                cif_value(&atom.label_atom_id, "_atom_site.label_atom_id")?,
                optional_cif_value(atom.label_alt_id.as_deref(), "_atom_site.label_alt_id")?,
                cif_value(&atom.label_comp_id, "_atom_site.label_comp_id")?,
                cif_value(&atom.asym_id, "_atom_site.label_asym_id")?,
                cif_value(&atom.entity_id, "_atom_site.label_entity_id")?,
                optional_display(atom.label_seq_id),
                optional_cif_value(
                    atom.insertion_code.as_deref(),
                    "_atom_site.pdbx_PDB_ins_code",
                )?,
                format_coordinate(atom.position.x, options.coordinate_precision),
                format_coordinate(atom.position.y, options.coordinate_precision),
                format_coordinate(atom.position.z, options.coordinate_precision),
                optional_float(atom.occupancy),
                optional_float(atom.b_factor),
                atom.formal_charge.to_string(),
                optional_cif_value(atom.auth_seq_id.as_deref(), "_atom_site.auth_seq_id")?,
                cif_value(&atom.auth_comp_id, "_atom_site.auth_comp_id")?,
                cif_value(&atom.auth_asym_id, "_atom_site.auth_asym_id")?,
                cif_value(&atom.auth_atom_id, "_atom_site.auth_atom_id")?,
                "1".to_owned(),
            ],
        );
    }
    output.push_str("#\n");

    if !model.connections.is_empty() {
        let indexes = model
            .atoms
            .iter()
            .enumerate()
            .map(|(index, row)| (row.atom, index))
            .collect::<BTreeMap<_, _>>();
        const CONNECTION_TAGS: &[&str] = &[
            "_struct_conn.id",
            "_struct_conn.conn_type_id",
            "_struct_conn.ptnr1_label_asym_id",
            "_struct_conn.ptnr1_label_comp_id",
            "_struct_conn.ptnr1_label_seq_id",
            "_struct_conn.ptnr1_label_atom_id",
            "_struct_conn.ptnr2_label_asym_id",
            "_struct_conn.ptnr2_label_comp_id",
            "_struct_conn.ptnr2_label_seq_id",
            "_struct_conn.ptnr2_label_atom_id",
            "_struct_conn.pdbx_value_order",
        ];
        write_loop_header(&mut output, CONNECTION_TAGS);
        for (serial, connection) in (1u64..).zip(model.connections.iter()) {
            let left = &model.atoms[indexes[&connection.left]];
            let right = &model.atoms[indexes[&connection.right]];
            write_row(
                &mut output,
                vec![
                    format!("conn{serial}"),
                    "covale".to_owned(),
                    cif_value(&left.asym_id, "_struct_conn.ptnr1_label_asym_id")?,
                    cif_value(&left.label_comp_id, "_struct_conn.ptnr1_label_comp_id")?,
                    optional_display(left.label_seq_id),
                    cif_value(&left.label_atom_id, "_struct_conn.ptnr1_label_atom_id")?,
                    cif_value(&right.asym_id, "_struct_conn.ptnr2_label_asym_id")?,
                    cif_value(&right.label_comp_id, "_struct_conn.ptnr2_label_comp_id")?,
                    optional_display(right.label_seq_id),
                    cif_value(&right.label_atom_id, "_struct_conn.ptnr2_label_atom_id")?,
                    bond_order_code(connection.order).to_owned(),
                ],
            );
        }
        output.push_str("#\n");
    }
    Ok(output)
}

fn write_loop_header(output: &mut String, tags: &[&str]) {
    output.push_str("loop_\n");
    for tag in tags {
        output.push_str(tag);
        output.push('\n');
    }
}

fn write_row(output: &mut String, values: Vec<String>) {
    output.push_str(&values.join(" "));
    output.push('\n');
}

fn format_coordinate(value: f64, precision: usize) -> String {
    format!("{value:.precision$}")
}

fn optional_float(value: Option<f64>) -> String {
    value.map_or_else(|| ".".to_owned(), |value| value.to_string())
}

fn optional_display(value: Option<i32>) -> String {
    value.map_or_else(|| ".".to_owned(), |value| value.to_string())
}

fn optional_cif_value(value: Option<&str>, field: &'static str) -> Result<String, MmcifWriteError> {
    value.map_or_else(|| Ok(".".to_owned()), |value| cif_value(value, field))
}

fn cif_value(value: &str, field: &'static str) -> Result<String, MmcifWriteError> {
    if value.is_empty() || value.contains(['\n', '\r']) {
        return Err(MmcifWriteError::UnsupportedTextValue { field });
    }
    let lower = value.to_ascii_lowercase();
    let is_control = lower == "loop_"
        || lower == "stop_"
        || lower == "global_"
        || lower.starts_with("data_")
        || lower.starts_with("save_")
        || value.starts_with('_');
    let bare = !value.is_empty()
        && value != "."
        && value != "?"
        && !is_control
        && !value
            .chars()
            .any(|character| character.is_ascii_whitespace() || character == '#')
        && !value.starts_with(';')
        && !value.contains(['\'', '"']);
    if bare {
        return Ok(value.to_owned());
    }
    if !value.contains('\'') {
        return Ok(format!("'{value}'"));
    }
    if !value.contains('"') {
        return Ok(format!("\"{value}\""));
    }
    Err(MmcifWriteError::UnsupportedTextValue { field })
}

fn bond_order_code(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "sing",
        BondOrder::Double => "doub",
        BondOrder::Triple => "trip",
        BondOrder::Quadruple => "quad",
        BondOrder::Zero | BondOrder::Dative => {
            unreachable!("unsupported bond order was rejected")
        }
    }
}

#[cfg(test)]
mod capacity_tests {
    use super::one_based_serial;

    #[test]
    fn one_based_serials_widen_before_incrementing() {
        assert_eq!(one_based_serial(0), 1);
        assert_eq!(one_based_serial(u32::MAX), u64::from(u32::MAX) + 1);
    }
}
