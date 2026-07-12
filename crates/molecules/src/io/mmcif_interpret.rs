use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::bio::{AtomSiteMetadata, BioHierarchy, MacroMolecule, MolecularContents};
use crate::core::{Atom, AtomId, BondOrder, Conformer, Element, Molecule, Point3, PropValue};
use crate::small::SmallMolecule;

use super::{MmcifDataBlock, MmcifDocument, MmcifLoopTable, MmcifValue};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifAltLocPolicy {
    HighestOccupancy,
    SelectLabel(String),
    ErrorOnAlternateLocations,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifInterpretOptions {
    pub strict_entity_metadata: bool,
    pub altloc_policy: MmcifAltLocPolicy,
}

impl Default for MmcifInterpretOptions {
    fn default() -> Self {
        Self {
            strict_entity_metadata: false,
            altloc_policy: MmcifAltLocPolicy::HighestOccupancy,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum MmcifEntityKind {
    Polymer,
    Branched,
    NonPolymer,
    Water,
    Other(String),
}

impl MmcifEntityKind {
    fn from_mmcif(value: &str) -> Self {
        match value.to_ascii_lowercase().as_str() {
            "polymer" => Self::Polymer,
            "branched" => Self::Branched,
            "non-polymer" => Self::NonPolymer,
            "water" => Self::Water,
            _ => Self::Other(value.to_owned()),
        }
    }

    fn as_str(&self) -> &str {
        match self {
            Self::Polymer => "polymer",
            Self::Branched => "branched",
            Self::NonPolymer => "non-polymer",
            Self::Water => "water",
            Self::Other(value) => value,
        }
    }

    fn is_macro(&self) -> bool {
        matches!(self, Self::Polymer | Self::Branched)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifInterpretIssue {
    EntityTypeInferred {
        asym_id: String,
        kind: MmcifEntityKind,
    },
    ConnectionIgnored {
        connection_type: String,
    },
    ConnectionUnresolved {
        connection_type: String,
    },
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MmcifInterpretationReport {
    pub data_block: String,
    pub entity_definitions: usize,
    pub coordinate_models: usize,
    pub macromolecules: usize,
    pub small_molecules: usize,
    pub solvent_molecules: usize,
    pub applied_connections: usize,
    pub template_bonds_pending: usize,
    pub issues: Vec<MmcifInterpretIssue>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MmcifInterpretation {
    contents: MolecularContents,
    report: MmcifInterpretationReport,
}

impl MmcifInterpretation {
    pub fn contents(&self) -> &MolecularContents {
        &self.contents
    }

    pub fn report(&self) -> &MmcifInterpretationReport {
        &self.report
    }

    pub fn into_parts(self) -> (MolecularContents, MmcifInterpretationReport) {
        (self.contents, self.report)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifInterpretError {
    pub line: Option<usize>,
    pub message: String,
}

impl MmcifInterpretError {
    fn new(line: Option<usize>, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }
}

impl fmt::Display for MmcifInterpretError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.line {
            Some(line) => write!(
                f,
                "mmCIF interpretation error at line {line}: {}",
                self.message
            ),
            None => write!(f, "mmCIF interpretation error: {}", self.message),
        }
    }
}

impl std::error::Error for MmcifInterpretError {}

pub fn interpret_mmcif(
    document: &MmcifDocument,
    options: MmcifInterpretOptions,
) -> Result<MmcifInterpretation, MmcifInterpretError> {
    let blocks = document
        .blocks()
        .iter()
        .filter(|block| block.loop_with_tag("_atom_site.type_symbol").is_some())
        .collect::<Vec<_>>();
    if blocks.is_empty() {
        return Err(MmcifInterpretError::new(
            None,
            "document has no atom-site loop",
        ));
    }
    if blocks.len() > 1 {
        return Err(MmcifInterpretError::new(
            None,
            "document has atom-site data in more than one data block",
        ));
    }
    interpret_block(blocks[0], options)
}

fn interpret_block(
    block: &MmcifDataBlock,
    options: MmcifInterpretOptions,
) -> Result<MmcifInterpretation, MmcifInterpretError> {
    let entities = read_entity_types(block)?;
    let asym_entities = read_asym_entities(block)?;
    let atom_table = block
        .loop_with_tag("_atom_site.type_symbol")
        .expect("selected block has atom-site data");
    if atom_table.row_count() == 0 {
        return Err(MmcifInterpretError::new(
            None,
            "atom-site loop contains no rows",
        ));
    }
    let mut report = MmcifInterpretationReport {
        data_block: block.name().to_owned(),
        entity_definitions: entities.len(),
        ..MmcifInterpretationReport::default()
    };
    let rows = read_atom_rows(atom_table, &entities, &asym_entities, &options, &mut report)?;
    let selected = select_alt_locations(rows, &options.altloc_policy)?;
    let mut union = InstanceUnion::new(selected.iter().map(|row| row.instance_key.clone()));
    let connections = read_connections(block, &selected, &mut union, &mut report)?;
    let groups = group_rows(selected, &mut union);
    let model_ids = groups
        .values()
        .flat_map(|group| group.rows.iter().map(|row| row.model_id.clone()))
        .collect::<BTreeSet<_>>();
    report.coordinate_models = model_ids.len();

    let mut contents = MolecularContents::new();
    for (_, group) in groups {
        let built = build_molecule(group, &connections, &mut report)?;
        match built {
            BuiltMolecule::Macro(molecule) => contents.push_macro(molecule),
            BuiltMolecule::Small(molecule) => contents.push_small(molecule),
            BuiltMolecule::Water(molecule) => contents.solvent_mut().push(molecule),
        }
    }
    report.macromolecules = contents.macromolecules().count();
    report.small_molecules = contents.small_molecules().count();
    report.solvent_molecules = contents.solvent().len();
    Ok(MmcifInterpretation { contents, report })
}

fn read_entity_types(
    block: &MmcifDataBlock,
) -> Result<BTreeMap<String, MmcifEntityKind>, MmcifInterpretError> {
    let mut entities = BTreeMap::new();
    if let Some(table) = block.loop_with_tag("_entity.id") {
        for row in 0..table.row_count() {
            let id = required(table, row, "_entity.id")?;
            let kind = required(table, row, "_entity.type")?;
            if entities
                .insert(id.to_owned(), MmcifEntityKind::from_mmcif(kind))
                .is_some()
            {
                return Err(row_error(table, row, format!("duplicate entity `{id}`")));
            }
        }
    } else if let (Some(id), Some(kind)) = (
        block.item("_entity.id").and_then(MmcifValue::optional_text),
        block
            .item("_entity.type")
            .and_then(MmcifValue::optional_text),
    ) {
        entities.insert(id.to_owned(), MmcifEntityKind::from_mmcif(kind));
    }
    Ok(entities)
}

fn read_asym_entities(
    block: &MmcifDataBlock,
) -> Result<BTreeMap<String, String>, MmcifInterpretError> {
    let mut instances = BTreeMap::new();
    if let Some(table) = block.loop_with_tag("_struct_asym.id") {
        for row in 0..table.row_count() {
            let id = required(table, row, "_struct_asym.id")?;
            let entity = required(table, row, "_struct_asym.entity_id")?;
            if instances.insert(id.to_owned(), entity.to_owned()).is_some() {
                return Err(row_error(
                    table,
                    row,
                    format!("duplicate structural instance `{id}`"),
                ));
            }
        }
    } else if let (Some(id), Some(entity)) = (
        block
            .item("_struct_asym.id")
            .and_then(MmcifValue::optional_text),
        block
            .item("_struct_asym.entity_id")
            .and_then(MmcifValue::optional_text),
    ) {
        instances.insert(id.to_owned(), entity.to_owned());
    }
    Ok(instances)
}

#[derive(Debug, Clone)]
struct AtomRow {
    line: usize,
    row_index: usize,
    model_id: String,
    entity_id: Option<String>,
    kind: MmcifEntityKind,
    instance_key: String,
    asym_id: String,
    auth_asym_id: Option<String>,
    residue_key: String,
    label_seq_id: Option<i32>,
    auth_seq_id: Option<String>,
    insertion_code: Option<String>,
    comp_id: String,
    auth_comp_id: Option<String>,
    atom_name: String,
    auth_atom_name: Option<String>,
    atom_site_id: Option<String>,
    group_pdb: Option<String>,
    alt_id: Option<String>,
    occupancy: Option<f64>,
    occupancy_raw: Option<String>,
    b_factor: Option<f64>,
    b_factor_raw: Option<String>,
    point: Option<Point3>,
    point_raw: [Option<String>; 3],
    element: Element,
    formal_charge: i8,
}

impl AtomRow {
    fn atom_key(&self) -> String {
        format!("{}|{}|{}", self.asym_id, self.residue_key, self.atom_name)
    }
}

#[derive(Debug, Default)]
struct OccurrenceState {
    occurrence: usize,
    seen: BTreeMap<String, BTreeSet<Option<String>>>,
}

fn read_atom_rows(
    table: &MmcifLoopTable,
    entities: &BTreeMap<String, MmcifEntityKind>,
    asym_entities: &BTreeMap<String, String>,
    options: &MmcifInterpretOptions,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<AtomRow>, MmcifInterpretError> {
    let mut rows = Vec::with_capacity(table.row_count());
    let mut occurrences = BTreeMap::<(String, String, String), OccurrenceState>::new();
    let mut inferred = BTreeSet::new();
    for row in 0..table.row_count() {
        let type_symbol = required(table, row, "_atom_site.type_symbol")?;
        let type_value = table
            .value(row, "_atom_site.type_symbol")
            .expect("required");
        let element = Element::from_symbol(&canonical_mmcif_element_symbol(type_symbol))
            .ok_or_else(|| {
                MmcifInterpretError::new(
                    Some(type_value.line()),
                    format!("unknown atom-site element `{type_symbol}`"),
                )
            })?;
        let asym_id = optional(table, row, "_atom_site.label_asym_id")
            .or_else(|| optional(table, row, "_atom_site.auth_asym_id"))
            .ok_or_else(|| row_error(table, row, "missing atom-site chain identifier"))?
            .to_owned();
        let auth_asym_id = optional(table, row, "_atom_site.auth_asym_id").map(str::to_owned);
        let comp_id = optional(table, row, "_atom_site.label_comp_id")
            .or_else(|| optional(table, row, "_atom_site.auth_comp_id"))
            .ok_or_else(|| row_error(table, row, "missing atom-site component identifier"))?
            .to_owned();
        let atom_name = optional(table, row, "_atom_site.label_atom_id")
            .or_else(|| optional(table, row, "_atom_site.auth_atom_id"))
            .ok_or_else(|| row_error(table, row, "missing atom-site atom identifier"))?
            .to_owned();
        let model_id = optional(table, row, "_atom_site.pdbx_PDB_model_num")
            .unwrap_or("1")
            .to_owned();
        let entity_id = optional(table, row, "_atom_site.label_entity_id")
            .map(str::to_owned)
            .or_else(|| asym_entities.get(&asym_id).cloned());
        let group_pdb = optional(table, row, "_atom_site.group_PDB").map(str::to_owned);
        let kind = entity_id
            .as_ref()
            .and_then(|entity| entities.get(entity))
            .cloned()
            .unwrap_or_else(|| infer_entity_kind(group_pdb.as_deref(), &comp_id));
        if entity_id
            .as_ref()
            .and_then(|entity| entities.get(entity))
            .is_none()
        {
            if options.strict_entity_metadata {
                return Err(row_error(
                    table,
                    row,
                    format!("missing entity type for structural instance `{asym_id}`"),
                ));
            }
            if inferred.insert(asym_id.clone()) {
                report.issues.push(MmcifInterpretIssue::EntityTypeInferred {
                    asym_id: asym_id.clone(),
                    kind: kind.clone(),
                });
            }
        }
        let label_seq_id = optional_i32(table, row, "_atom_site.label_seq_id")?;
        let auth_seq_id = optional(table, row, "_atom_site.auth_seq_id").map(str::to_owned);
        let insertion_code =
            optional(table, row, "_atom_site.pdbx_PDB_ins_code").map(str::to_owned);
        let alt_id = optional(table, row, "_atom_site.label_alt_id").map(str::to_owned);
        let residue_key = if let Some(sequence) = label_seq_id {
            format!(
                "label:{sequence}:{}",
                insertion_code.as_deref().unwrap_or("")
            )
        } else if let Some(sequence) = &auth_seq_id {
            format!(
                "auth:{sequence}:{}",
                insertion_code.as_deref().unwrap_or("")
            )
        } else {
            let state = occurrences
                .entry((model_id.clone(), asym_id.clone(), comp_id.clone()))
                .or_default();
            let prior = state.seen.get(&atom_name);
            let repeats = prior.is_some_and(|labels| {
                alt_id.is_none() || labels.contains(&None) || labels.contains(&alt_id)
            });
            if repeats {
                state.occurrence += 1;
                state.seen.clear();
            }
            state
                .seen
                .entry(atom_name.clone())
                .or_default()
                .insert(alt_id.clone());
            format!("occurrence:{}", state.occurrence)
        };
        let instance_key = if kind.is_macro() {
            format!("macro:{asym_id}")
        } else {
            format!("small:{asym_id}:{residue_key}")
        };
        let formal_charge =
            optional_i8(table, row, "_atom_site.pdbx_formal_charge")?.unwrap_or_default();
        let occupancy_raw = optional(table, row, "_atom_site.occupancy").map(str::to_owned);
        let b_factor_raw = optional(table, row, "_atom_site.B_iso_or_equiv").map(str::to_owned);
        let x_raw = optional(table, row, "_atom_site.Cartn_x").map(str::to_owned);
        let y_raw = optional(table, row, "_atom_site.Cartn_y").map(str::to_owned);
        let z_raw = optional(table, row, "_atom_site.Cartn_z").map(str::to_owned);
        let point = optional_point(table, row)?;
        rows.push(AtomRow {
            line: type_value.line(),
            row_index: row,
            model_id,
            entity_id,
            kind,
            instance_key,
            asym_id,
            auth_asym_id,
            residue_key,
            label_seq_id,
            auth_seq_id,
            insertion_code,
            comp_id,
            auth_comp_id: optional(table, row, "_atom_site.auth_comp_id").map(str::to_owned),
            atom_name,
            auth_atom_name: optional(table, row, "_atom_site.auth_atom_id").map(str::to_owned),
            atom_site_id: optional(table, row, "_atom_site.id").map(str::to_owned),
            group_pdb,
            alt_id,
            occupancy: optional_f64(table, row, "_atom_site.occupancy")?,
            occupancy_raw,
            b_factor: optional_f64(table, row, "_atom_site.B_iso_or_equiv")?,
            b_factor_raw,
            point,
            point_raw: [x_raw, y_raw, z_raw],
            element,
            formal_charge,
        });
    }
    Ok(rows)
}

fn infer_entity_kind(group_pdb: Option<&str>, comp_id: &str) -> MmcifEntityKind {
    if ["HOH", "WAT", "DOD"]
        .iter()
        .any(|water| comp_id.eq_ignore_ascii_case(water))
    {
        MmcifEntityKind::Water
    } else if group_pdb.is_some_and(|group| group.eq_ignore_ascii_case("ATOM")) {
        MmcifEntityKind::Polymer
    } else {
        MmcifEntityKind::NonPolymer
    }
}

fn canonical_mmcif_element_symbol(symbol: &str) -> String {
    let mut chars = symbol.chars();
    let Some(first) = chars.next() else {
        return String::new();
    };
    let mut canonical = first.to_ascii_uppercase().to_string();
    canonical.extend(chars.flat_map(char::to_lowercase));
    canonical
}

fn select_alt_locations(
    rows: Vec<AtomRow>,
    policy: &MmcifAltLocPolicy,
) -> Result<Vec<AtomRow>, MmcifInterpretError> {
    let mut grouped = BTreeMap::<(String, String, String), Vec<AtomRow>>::new();
    for row in rows {
        grouped
            .entry((
                row.instance_key.clone(),
                row.atom_key(),
                row.model_id.clone(),
            ))
            .or_default()
            .push(row);
    }
    let mut selected = Vec::new();
    for (_, mut candidates) in grouped {
        candidates.sort_by_key(|row| row.row_index);
        let mut identities = BTreeSet::new();
        if let Some(duplicate) = candidates
            .iter()
            .find(|row| !identities.insert(row.alt_id.clone()))
        {
            return Err(MmcifInterpretError::new(
                Some(duplicate.line),
                format!(
                    "atom `{}` has duplicate records for one alternate location",
                    duplicate.atom_name
                ),
            ));
        }
        let labels = candidates
            .iter()
            .filter_map(|row| row.alt_id.clone())
            .collect::<BTreeSet<_>>();
        if candidates.len() > 1
            && !labels.is_empty()
            && matches!(policy, MmcifAltLocPolicy::ErrorOnAlternateLocations)
        {
            return Err(MmcifInterpretError::new(
                Some(candidates[0].line),
                format!("atom `{}` has alternate locations", candidates[0].atom_name),
            ));
        }
        let chosen = match policy {
            MmcifAltLocPolicy::HighestOccupancy => candidates.into_iter().max_by(|left, right| {
                left.occupancy
                    .unwrap_or(0.0)
                    .partial_cmp(&right.occupancy.unwrap_or(0.0))
                    .unwrap_or(Ordering::Equal)
                    .then_with(|| right.alt_id.cmp(&left.alt_id))
            }),
            MmcifAltLocPolicy::SelectLabel(label) => candidates
                .iter()
                .find(|row| row.alt_id.as_deref() == Some(label.as_str()))
                .cloned()
                .or_else(|| candidates.iter().find(|row| row.alt_id.is_none()).cloned()),
            MmcifAltLocPolicy::ErrorOnAlternateLocations => candidates.into_iter().next(),
        };
        let Some(chosen) = chosen else {
            return Err(MmcifInterpretError::new(
                None,
                "requested alternate-location label is unavailable",
            ));
        };
        selected.push(chosen);
    }
    selected.sort_by_key(|row| row.row_index);
    Ok(selected)
}

#[derive(Debug, Clone)]
struct DeclaredConnection {
    left_atom: String,
    right_atom: String,
}

fn read_connections(
    block: &MmcifDataBlock,
    rows: &[AtomRow],
    union: &mut InstanceUnion,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<DeclaredConnection>, MmcifInterpretError> {
    let Some(table) = block.loop_with_tag("_struct_conn.conn_type_id") else {
        return Ok(Vec::new());
    };
    let mut connections = Vec::new();
    for row in 0..table.row_count() {
        let kind = required(table, row, "_struct_conn.conn_type_id")?.to_owned();
        if !is_covalent_connection(&kind) {
            report.issues.push(MmcifInterpretIssue::ConnectionIgnored {
                connection_type: kind,
            });
            continue;
        }
        let left = connection_partner(table, row, 1, rows)?;
        let right = connection_partner(table, row, 2, rows)?;
        let (Some(left), Some(right)) = (left, right) else {
            report
                .issues
                .push(MmcifInterpretIssue::ConnectionUnresolved {
                    connection_type: kind,
                });
            continue;
        };
        union.union(&left.instance_key, &right.instance_key);
        connections.push(DeclaredConnection {
            left_atom: left.atom_key(),
            right_atom: right.atom_key(),
        });
        report.applied_connections += 1;
    }
    Ok(connections)
}

fn is_covalent_connection(kind: &str) -> bool {
    let kind = kind.to_ascii_lowercase();
    kind.starts_with("covale") || kind == "disulf" || kind == "modres"
}

fn connection_partner<'a>(
    table: &MmcifLoopTable,
    row: usize,
    partner: usize,
    rows: &'a [AtomRow],
) -> Result<Option<&'a AtomRow>, MmcifInterpretError> {
    let asym_tag = format!("_struct_conn.ptnr{partner}_label_asym_id");
    let atom_tag = format!("_struct_conn.ptnr{partner}_label_atom_id");
    let seq_tag = format!("_struct_conn.ptnr{partner}_label_seq_id");
    let asym = optional(table, row, &asym_tag);
    let atom = optional(table, row, &atom_tag);
    let seq = optional(table, row, &seq_tag);
    let (Some(asym), Some(atom)) = (asym, atom) else {
        return Ok(None);
    };
    Ok(rows.iter().find(|candidate| {
        candidate.asym_id == asym
            && candidate.atom_name == atom
            && seq.is_none_or(|seq| {
                candidate
                    .label_seq_id
                    .map(|value| value.to_string())
                    .as_deref()
                    == Some(seq)
            })
    }))
}

#[derive(Debug)]
struct MoleculeGroup {
    rows: Vec<AtomRow>,
    kinds: BTreeSet<MmcifEntityKind>,
    instance_keys: BTreeSet<String>,
}

fn group_rows(rows: Vec<AtomRow>, union: &mut InstanceUnion) -> BTreeMap<String, MoleculeGroup> {
    let mut groups = BTreeMap::new();
    for row in rows {
        let root = union.find(&row.instance_key);
        let group = groups.entry(root).or_insert_with(|| MoleculeGroup {
            rows: Vec::new(),
            kinds: BTreeSet::new(),
            instance_keys: BTreeSet::new(),
        });
        group.kinds.insert(row.kind.clone());
        group.instance_keys.insert(row.instance_key.clone());
        group.rows.push(row);
    }
    groups
}

enum BuiltMolecule {
    Small(SmallMolecule),
    Macro(MacroMolecule),
    Water(SmallMolecule),
}

fn build_molecule(
    group: MoleculeGroup,
    connections: &[DeclaredConnection],
    report: &mut MmcifInterpretationReport,
) -> Result<BuiltMolecule, MmcifInterpretError> {
    let is_macro = group.kinds.iter().any(MmcifEntityKind::is_macro);
    let is_water = group
        .kinds
        .iter()
        .all(|kind| *kind == MmcifEntityKind::Water);
    let mut graph = Molecule::new();
    let mut atoms = BTreeMap::new();
    let mut representative = Vec::<(String, AtomRow)>::new();
    let mut seen_atoms = BTreeMap::<String, usize>::new();
    for row in &group.rows {
        let key = row.atom_key();
        if let Some(&index) = seen_atoms.get(&key) {
            let prior = &representative[index].1;
            if prior.element != row.element
                || prior.formal_charge != row.formal_charge
                || prior.comp_id != row.comp_id
                || prior.entity_id != row.entity_id
            {
                return Err(MmcifInterpretError::new(
                    Some(row.line),
                    format!(
                        "atom `{}` has inconsistent topology payload across coordinate models",
                        row.atom_name
                    ),
                ));
            }
        } else {
            seen_atoms.insert(key.clone(), representative.len());
            representative.push((key, row.clone()));
        }
    }
    for (key, row) in &representative {
        let mut atom = Atom::new(row.element);
        atom.formal_charge = row.formal_charge;
        atom.props.insert(
            "mmcif.atom_id".into(),
            PropValue::String(row.atom_name.clone()),
        );
        atom.props.insert(
            "mmcif.comp_id".into(),
            PropValue::String(row.comp_id.clone()),
        );
        atom.props.insert(
            "mmcif.asym_id".into(),
            PropValue::String(row.asym_id.clone()),
        );
        atoms.insert(key.clone(), graph.add_atom(atom));
    }
    let mut seen_models = BTreeSet::new();
    let model_ids = group
        .rows
        .iter()
        .filter(|row| seen_models.insert(row.model_id.clone()))
        .map(|row| row.model_id.clone())
        .collect::<Vec<_>>();
    for model_id in model_ids {
        let mut conformer = Conformer::new();
        conformer
            .props_mut()
            .insert("mmcif.model_id".into(), PropValue::String(model_id.clone()));
        let mut has_positions = false;
        for row in group.rows.iter().filter(|row| row.model_id == model_id) {
            if let Some(point) = row.point {
                conformer.set_position(atoms[&row.atom_key()], point);
                has_positions = true;
            }
        }
        if has_positions {
            graph.add_conformer(conformer);
        }
    }
    for connection in connections {
        let Some(&left) = atoms.get(&connection.left_atom) else {
            continue;
        };
        let Some(&right) = atoms.get(&connection.right_atom) else {
            continue;
        };
        if graph
            .bond_between(left, right)
            .map_err(graph_error)?
            .is_none()
        {
            graph
                .add_bond(left, right, BondOrder::Single)
                .map_err(graph_error)?;
        }
    }
    let asym_ids = representative
        .iter()
        .map(|(_, row)| row)
        .map(|row| row.asym_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>()
        .join(",");
    let entity_ids = representative
        .iter()
        .map(|(_, row)| row)
        .filter_map(|row| row.entity_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>()
        .join(",");
    graph
        .props_mut()
        .insert("mmcif.asym_ids".into(), PropValue::String(asym_ids));
    if !entity_ids.is_empty() {
        graph
            .props_mut()
            .insert("mmcif.entity_ids".into(), PropValue::String(entity_ids));
    }
    graph.props_mut().insert(
        "mmcif.entity_kinds".into(),
        PropValue::String(
            group
                .kinds
                .iter()
                .map(MmcifEntityKind::as_str)
                .collect::<Vec<_>>()
                .join(","),
        ),
    );
    if graph.atom_count() > 1 {
        report.template_bonds_pending += 1;
    }

    if is_macro {
        let hierarchy = build_hierarchy(&graph, &representative, &atoms)?;
        Ok(BuiltMolecule::Macro(MacroMolecule::from_parts(
            graph, hierarchy,
        )))
    } else {
        let molecule = SmallMolecule::from_graph(graph);
        if is_water {
            Ok(BuiltMolecule::Water(molecule))
        } else {
            Ok(BuiltMolecule::Small(molecule))
        }
    }
}

fn build_hierarchy(
    graph: &Molecule,
    representative: &[(String, AtomRow)],
    atoms: &BTreeMap<String, AtomId>,
) -> Result<BioHierarchy, MmcifInterpretError> {
    let mut hierarchy = BioHierarchy::new();
    let model = hierarchy.add_model("structure");
    let mut chains = BTreeMap::new();
    let mut residues = BTreeMap::new();
    for (key, row) in representative {
        let chain = if let Some(chain) = chains.get(&row.asym_id) {
            *chain
        } else {
            let chain = hierarchy
                .add_chain(model, row.asym_id.clone(), row.auth_asym_id.clone())
                .map_err(hierarchy_error)?;
            chains.insert(row.asym_id.clone(), chain);
            chain
        };
        let residue_key = (row.asym_id.clone(), row.residue_key.clone());
        let residue = if let Some(residue) = residues.get(&residue_key) {
            *residue
        } else {
            let residue = hierarchy
                .add_residue(
                    chain,
                    row.comp_id.clone(),
                    row.label_seq_id,
                    row.auth_seq_id.clone(),
                    row.insertion_code.clone(),
                )
                .map_err(hierarchy_error)?;
            let record = &mut hierarchy.residues[residue.index()];
            record.label_comp_id = Some(row.comp_id.clone());
            record.author_comp_id = row.auth_comp_id.clone();
            residues.insert(residue_key, residue);
            residue
        };
        let atom = atoms[key];
        graph.atom(atom).map_err(graph_error)?;
        hierarchy
            .add_atom_site(
                residue,
                atom,
                AtomSiteMetadata {
                    group_pdb: row.group_pdb.clone(),
                    atom_site_id: row.atom_site_id.clone(),
                    type_symbol: Some(row.element.symbol().to_owned()),
                    label_asym_id: Some(row.asym_id.clone()),
                    auth_asym_id: row.auth_asym_id.clone(),
                    label_atom_id: Some(row.atom_name.clone()),
                    auth_atom_id: row.auth_atom_name.clone(),
                    label_alt_id: row.alt_id.clone(),
                    occupancy: row.occupancy,
                    occupancy_raw: row.occupancy_raw.clone(),
                    b_factor: row.b_factor,
                    b_factor_raw: row.b_factor_raw.clone(),
                    cartn_x_raw: row.point_raw[0].clone(),
                    cartn_y_raw: row.point_raw[1].clone(),
                    cartn_z_raw: row.point_raw[2].clone(),
                },
            )
            .map_err(hierarchy_error)?;
    }
    Ok(hierarchy)
}

fn graph_error(error: impl fmt::Display) -> MmcifInterpretError {
    MmcifInterpretError::new(None, error.to_string())
}

fn hierarchy_error(error: impl fmt::Display) -> MmcifInterpretError {
    MmcifInterpretError::new(None, error.to_string())
}

#[derive(Debug)]
struct InstanceUnion {
    parent: BTreeMap<String, String>,
}

impl InstanceUnion {
    fn new(keys: impl IntoIterator<Item = String>) -> Self {
        let parent = keys.into_iter().map(|key| (key.clone(), key)).collect();
        Self { parent }
    }

    fn find(&mut self, key: &str) -> String {
        let mut current = key.to_owned();
        let mut path = Vec::new();
        loop {
            let parent = self
                .parent
                .get(&current)
                .cloned()
                .unwrap_or_else(|| current.clone());
            if parent == current {
                break;
            }
            path.push(current);
            current = parent;
        }
        self.parent
            .entry(current.clone())
            .or_insert_with(|| current.clone());
        for node in path {
            self.parent.insert(node, current.clone());
        }
        current
    }

    fn union(&mut self, left: &str, right: &str) {
        let left_root = self.find(left);
        let right_root = self.find(right);
        if left_root != right_root {
            let (root, child) = if left_root < right_root {
                (left_root, right_root)
            } else {
                (right_root, left_root)
            };
            self.parent.insert(child, root);
        }
    }
}

fn required<'a>(
    table: &'a MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<&'a str, MmcifInterpretError> {
    let value = table
        .value(row, tag)
        .ok_or_else(|| row_error(table, row, format!("missing required {tag}")))?;
    value.optional_text().ok_or_else(|| {
        MmcifInterpretError::new(Some(value.line()), format!("missing required {tag}"))
    })
}

fn optional<'a>(table: &'a MmcifLoopTable, row: usize, tag: &str) -> Option<&'a str> {
    table.value(row, tag).and_then(MmcifValue::optional_text)
}

fn optional_f64(
    table: &MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<Option<f64>, MmcifInterpretError> {
    optional(table, row, tag)
        .map(|value| {
            let parsed = value
                .parse::<f64>()
                .map_err(|_| row_error(table, row, format!("invalid float {tag}")))?;
            if !parsed.is_finite() {
                return Err(row_error(table, row, format!("non-finite float {tag}")));
            }
            Ok(parsed)
        })
        .transpose()
}

fn optional_i32(
    table: &MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<Option<i32>, MmcifInterpretError> {
    optional(table, row, tag)
        .map(|value| {
            value
                .parse::<i32>()
                .map_err(|_| row_error(table, row, format!("invalid integer {tag}")))
        })
        .transpose()
}

fn optional_i8(
    table: &MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<Option<i8>, MmcifInterpretError> {
    optional(table, row, tag)
        .map(|value| {
            value
                .parse::<i8>()
                .map_err(|_| row_error(table, row, format!("invalid integer {tag}")))
        })
        .transpose()
}

fn optional_point(
    table: &MmcifLoopTable,
    row: usize,
) -> Result<Option<Point3>, MmcifInterpretError> {
    let x = optional_f64(table, row, "_atom_site.Cartn_x")?;
    let y = optional_f64(table, row, "_atom_site.Cartn_y")?;
    let z = optional_f64(table, row, "_atom_site.Cartn_z")?;
    match (x, y, z) {
        (Some(x), Some(y), Some(z)) => Ok(Some(Point3::new(x, y, z))),
        (None, None, None) => Ok(None),
        _ => Err(row_error(
            table,
            row,
            "partial atom-site coordinate triplet",
        )),
    }
}

fn row_error(
    table: &MmcifLoopTable,
    row: usize,
    message: impl Into<String>,
) -> MmcifInterpretError {
    MmcifInterpretError::new(
        table
            .row(row)
            .and_then(|row| row.first())
            .map(MmcifValue::line),
        message,
    )
}
