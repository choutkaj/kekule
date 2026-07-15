use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;

use crate::bio::{MacroValidateOptions, SmcraAtomSite, SmcraHierarchy};
use crate::core::{AtomId, BondOrder, Point3};
use crate::modeling::{
    InstanceAtomId, InstanceBondId, Model, MoleculeInstance, MoleculeInstanceId, MoleculeRole,
};

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

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifWriteError {
    InvalidDataBlockName(String),
    CoordinatePrecisionTooLarge(usize),
    InvalidModel(String),
    InvalidHierarchy {
        molecule: MoleculeInstanceId,
        message: String,
    },
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
    ConflictingEntityRoles(MoleculeInstanceId),
    EntityRolePayloadMismatch(MoleculeInstanceId),
    UnsupportedMoleculeRole {
        molecule: MoleculeInstanceId,
        role: MoleculeRole,
    },
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
    UnrepresentableInstanceBoundary(MoleculeInstanceId),
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
            Self::InvalidHierarchy { molecule, message } => {
                write!(f, "invalid hierarchy for {molecule}: {message}")
            }
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
            Self::ConflictingEntityRoles(molecule) => write!(
                f,
                "{molecule} has multiple mmCIF entity-kind roles that cannot be represented losslessly"
            ),
            Self::EntityRolePayloadMismatch(molecule) => write!(
                f,
                "{molecule} has an mmCIF entity-kind role inconsistent with its Small/Macro payload"
            ),
            Self::UnsupportedMoleculeRole { molecule, role } => write!(
                f,
                "{molecule} has model role {role:?}, which the foundational mmCIF writer cannot encode losslessly"
            ),
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
            Self::UnrepresentableInstanceBoundary(molecule) => write!(
                f,
                "{molecule} spans disconnected structural-instance IDs and cannot be represented losslessly"
            ),
            Self::UnsupportedTextValue { field } => write!(
                f,
                "{field} contains a text value that cannot be emitted as a single mmCIF token"
            ),
        }
    }
}

impl std::error::Error for MmcifWriteError {}

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
            Self::Polymer | Self::Branched => "ATOM",
            Self::NonPolymer | Self::Water => "HETATM",
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
struct AtomRow {
    atom: InstanceAtomId,
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
    bond: InstanceBondId,
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
    let prepared = prepare_model(model)?;
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

fn prepare_model(model: &Model) -> Result<PreparedModel, MmcifWriteError> {
    let mut reserved_asym_ids = BTreeSet::new();
    for (_, molecule) in model.topology().molecules() {
        if let Some(hierarchy) = molecule.hierarchy() {
            for (_, chain) in hierarchy.chains() {
                if chain.label_id.is_empty() {
                    return Err(invalid_hierarchy(molecule, "chain label ID is empty"));
                }
                if !reserved_asym_ids.insert(chain.label_id.clone()) {
                    return Err(MmcifWriteError::DuplicateAsymId(chain.label_id.clone()));
                }
            }
        }
    }

    let mut small_asym_ids = BTreeMap::new();
    for (id, molecule) in model.topology().molecules() {
        if molecule.small_molecule().is_none() {
            continue;
        }
        let base = format!("M{}", id.raw() + 1);
        let mut candidate = base.clone();
        let mut suffix = 2usize;
        while reserved_asym_ids.contains(&candidate) {
            candidate = format!("{base}_{suffix}");
            suffix += 1;
        }
        reserved_asym_ids.insert(candidate.clone());
        small_asym_ids.insert(id, candidate);
    }

    let mut entities = Vec::new();
    let mut asyms = Vec::new();
    let mut atoms = Vec::new();
    let mut asym_seen = BTreeSet::new();
    for (id, molecule) in model.topology().molecules() {
        validate_graph_chemistry(molecule)?;
        let entity_id = (id.raw() + 1).to_string();
        let kind = entity_kind(molecule)?;
        entities.push(EntityRow {
            id: entity_id.clone(),
            kind,
        });
        if let Some(hierarchy) = molecule.hierarchy() {
            collect_macro_rows(
                model,
                molecule,
                hierarchy,
                &entity_id,
                kind,
                &mut asyms,
                &mut asym_seen,
                &mut atoms,
            )?;
        } else {
            let asym_id = small_asym_ids
                .get(&id)
                .expect("small-molecule asym ID was allocated")
                .clone();
            asyms.push(AsymRow {
                id: asym_id.clone(),
                entity_id: entity_id.clone(),
            });
            asym_seen.insert(asym_id.clone());
            collect_small_rows(model, molecule, &entity_id, &asym_id, kind, &mut atoms)?;
        }
    }

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
            .molecule(bond_id.molecule())
            .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?;
        let left = molecule.qualify_atom(bond.a());
        let right = molecule.qualify_atom(bond.b());
        validate_connection_selector(left, &atoms, &atom_indexes)?;
        validate_connection_selector(right, &atoms, &atom_indexes)?;
        connections.push(ConnectionRow {
            bond: bond_id,
            left,
            right,
            order,
        });
    }
    validate_instance_boundaries(model, &atoms, &connections)?;

    Ok(PreparedModel {
        entities,
        asyms,
        atoms,
        connections,
    })
}

fn entity_kind(molecule: &MoleculeInstance) -> Result<EntityKind, MmcifWriteError> {
    for role in [MoleculeRole::Ligand, MoleculeRole::Cofactor] {
        if molecule.has_role(role) {
            return Err(MmcifWriteError::UnsupportedMoleculeRole {
                molecule: molecule.id(),
                role,
            });
        }
    }
    let inferred_ion = molecule.graph().atom_count() == 1
        && molecule
            .graph()
            .atoms()
            .next()
            .is_some_and(|(_, atom)| atom.formal_charge != 0);
    if molecule.has_role(MoleculeRole::Ion) != inferred_ion {
        return Err(MmcifWriteError::UnsupportedMoleculeRole {
            molecule: molecule.id(),
            role: MoleculeRole::Ion,
        });
    }
    let primary = [
        MoleculeRole::Polymer,
        MoleculeRole::Branched,
        MoleculeRole::NonPolymer,
        MoleculeRole::Solvent,
    ]
    .into_iter()
    .filter(|role| molecule.has_role(*role))
    .collect::<Vec<_>>();
    if primary.len() > 1 {
        return Err(MmcifWriteError::ConflictingEntityRoles(molecule.id()));
    }
    let kind = match primary.first().copied() {
        Some(MoleculeRole::Polymer) => EntityKind::Polymer,
        Some(MoleculeRole::Branched) => EntityKind::Branched,
        Some(MoleculeRole::NonPolymer) => EntityKind::NonPolymer,
        Some(MoleculeRole::Solvent) => EntityKind::Water,
        Some(_) => unreachable!("primary roles are exhaustive"),
        None if molecule.macro_molecule().is_some() => EntityKind::Polymer,
        None => EntityKind::NonPolymer,
    };
    let macro_kind = matches!(kind, EntityKind::Polymer | EntityKind::Branched);
    if macro_kind != molecule.macro_molecule().is_some() {
        return Err(MmcifWriteError::EntityRolePayloadMismatch(molecule.id()));
    }
    Ok(kind)
}

fn validate_graph_chemistry(molecule: &MoleculeInstance) -> Result<(), MmcifWriteError> {
    let graph = molecule.graph();
    if graph.stereo_elements().next().is_some()
        || graph.stereo_groups().next().is_some()
        || graph.stereo_bond_marks().next().is_some()
    {
        return Err(MmcifWriteError::UnsupportedStereo(molecule.id()));
    }
    for (atom_id, atom) in graph.atoms() {
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
        if atom.explicit_hydrogens != 0 {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "explicit_hydrogens",
            });
        }
        if atom.no_implicit_hydrogens {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "no_implicit_hydrogens",
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
    hierarchy: &SmcraHierarchy,
    entity_id: &str,
    kind: EntityKind,
    asyms: &mut Vec<AsymRow>,
    asym_seen: &mut BTreeSet<String>,
    rows: &mut Vec<AtomRow>,
) -> Result<(), MmcifWriteError> {
    molecule
        .macro_molecule()
        .expect("hierarchy implies macro molecule")
        .validate_with_options(MacroValidateOptions {
            validate_coordinates: false,
        })
        .map_err(|error| invalid_hierarchy(molecule, error.to_string()))?;
    if hierarchy.models().count() != 1 {
        return Err(invalid_hierarchy(
            molecule,
            "foundational mmCIF writing requires exactly one hierarchy model",
        ));
    }

    let mut sites = BTreeMap::<AtomId, &SmcraAtomSite>::new();
    for (_, site) in hierarchy.atom_sites() {
        if sites.insert(site.atom, site).is_some() {
            return Err(MmcifWriteError::DuplicateAtomSite(
                molecule.qualify_atom(site.atom),
            ));
        }
    }
    for (_, chain) in hierarchy.chains() {
        if asym_seen.insert(chain.label_id.clone()) {
            asyms.push(AsymRow {
                id: chain.label_id.clone(),
                entity_id: entity_id.to_owned(),
            });
        }
    }

    for (atom_id, atom) in molecule.graph().atoms() {
        let qualified = molecule.qualify_atom(atom_id);
        let site = sites
            .get(&atom_id)
            .copied()
            .ok_or(MmcifWriteError::MissingAtomSite(qualified))?;
        let residue = hierarchy
            .residue(site.residue)
            .map_err(|error| invalid_hierarchy(molecule, error.to_string()))?;
        if residue.label_seq_id.is_none() && residue.author_seq_id.is_none() {
            return Err(invalid_hierarchy(
                molecule,
                format!(
                    "residue {} has neither a label nor author sequence identifier",
                    site.residue.raw()
                ),
            ));
        }
        let chain = hierarchy
            .chain(residue.chain)
            .map_err(|error| invalid_hierarchy(molecule, error.to_string()))?;
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
        let group_pdb = normalized_group_pdb(
            qualified,
            site.metadata.group_pdb.as_deref(),
            kind.default_group_pdb(),
        )?;
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
            entity_id: entity_id.to_owned(),
            asym_id: chain.label_id.clone(),
            group_pdb,
            type_symbol: atom.element.symbol().to_owned(),
            label_atom_id: label_atom_id.clone(),
            label_alt_id: site.metadata.label_alt_id.clone(),
            label_comp_id: label_comp_id.clone(),
            label_seq_id: residue.label_seq_id,
            insertion_code: residue.insertion_code.clone(),
            position: model
                .position(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            occupancy: site.metadata.occupancy,
            b_factor: site.metadata.b_factor,
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
    entity_id: &str,
    asym_id: &str,
    kind: EntityKind,
    rows: &mut Vec<AtomRow>,
) -> Result<(), MmcifWriteError> {
    let component_id = if kind == EntityKind::Water {
        "HOH"
    } else {
        "MOL"
    };
    for (atom_id, atom) in molecule.graph().atoms() {
        let qualified = molecule.qualify_atom(atom_id);
        let atom_name = generated_atom_name(atom.element.symbol(), atom_id);
        rows.push(AtomRow {
            atom: qualified,
            entity_id: entity_id.to_owned(),
            asym_id: asym_id.to_owned(),
            group_pdb: kind.default_group_pdb().to_owned(),
            type_symbol: atom.element.symbol().to_owned(),
            label_atom_id: atom_name.clone(),
            label_alt_id: None,
            label_comp_id: component_id.to_owned(),
            label_seq_id: None,
            insertion_code: None,
            position: model
                .position(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            occupancy: None,
            b_factor: None,
            formal_charge: atom.formal_charge,
            auth_seq_id: None,
            auth_comp_id: component_id.to_owned(),
            auth_asym_id: asym_id.to_owned(),
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
    format!("{symbol}{}", atom.raw() + 1)
}

fn supported_bond_order(
    bond: InstanceBondId,
    order: BondOrder,
) -> Result<BondOrder, MmcifWriteError> {
    match order {
        BondOrder::Single | BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple => {
            Ok(order)
        }
        BondOrder::Zero | BondOrder::Aromatic | BondOrder::Dative => {
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
            "unsequenced".to_owned()
        };
        if !identities.insert((row.asym_id.clone(), residue, row.label_atom_id.clone())) {
            return Err(MmcifWriteError::DuplicateAtomIdentity(row.atom));
        }
    }
    Ok(())
}

fn validate_instance_boundaries(
    model: &Model,
    atoms: &[AtomRow],
    connections: &[ConnectionRow],
) -> Result<(), MmcifWriteError> {
    let atom_rows = atoms
        .iter()
        .map(|row| (row.atom, row))
        .collect::<BTreeMap<_, _>>();
    for (molecule_id, molecule) in model.topology().molecules() {
        let asym_ids = atoms
            .iter()
            .filter(|row| row.atom.molecule() == molecule_id)
            .map(|row| row.asym_id.clone())
            .collect::<BTreeSet<_>>();
        if asym_ids.len() <= 1 {
            continue;
        }
        let mut adjacency = BTreeMap::<String, BTreeSet<String>>::new();
        for asym_id in &asym_ids {
            adjacency.entry(asym_id.clone()).or_default();
        }
        for connection in connections
            .iter()
            .filter(|connection| connection.bond.molecule() == molecule_id)
        {
            let left = &atom_rows[&connection.left].asym_id;
            let right = &atom_rows[&connection.right].asym_id;
            if left != right {
                adjacency
                    .entry(left.clone())
                    .or_default()
                    .insert(right.clone());
                adjacency
                    .entry(right.clone())
                    .or_default()
                    .insert(left.clone());
            }
        }
        let start = asym_ids.iter().next().expect("multiple asym IDs").clone();
        let mut queue = VecDeque::from([start]);
        let mut visited = BTreeSet::new();
        while let Some(asym_id) = queue.pop_front() {
            if !visited.insert(asym_id.clone()) {
                continue;
            }
            queue.extend(adjacency[&asym_id].iter().cloned());
        }
        if visited != asym_ids {
            return Err(MmcifWriteError::UnrepresentableInstanceBoundary(
                molecule.id(),
            ));
        }
    }
    Ok(())
}

fn invalid_hierarchy(molecule: &MoleculeInstance, message: impl Into<String>) -> MmcifWriteError {
    MmcifWriteError::InvalidHierarchy {
        molecule: molecule.id(),
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
    for (index, atom) in model.atoms.iter().enumerate() {
        write_row(
            &mut output,
            vec![
                atom.group_pdb.clone(),
                (index + 1).to_string(),
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
        for (index, connection) in model.connections.iter().enumerate() {
            let left = &model.atoms[indexes[&connection.left]];
            let right = &model.atoms[indexes[&connection.right]];
            write_row(
                &mut output,
                vec![
                    format!("conn{}", index + 1),
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
        BondOrder::Zero | BondOrder::Aromatic | BondOrder::Dative => {
            unreachable!("unsupported bond order was rejected")
        }
    }
}
