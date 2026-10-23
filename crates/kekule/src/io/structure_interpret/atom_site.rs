//! Reads mmCIF `_atom_site`, `_entity`, and `_struct_asym` rows into typed
//! atom rows, the source-neutral input of interpretation.
use super::types::{poll_every, structure_error};
use crate::diagnostics::{Cancellation, DiagnosticKind, FormatError};
use std::collections::{BTreeMap, BTreeSet};

use crate::geometry::Point3;
use crate::molecule::Element;

use super::super::element_symbol::parse_structure_element_symbol;
use super::super::mmcif_category::MmcifCategory;
use super::super::{MmcifBlock, MmcifValue};
use super::types::{EntityKind, StructureInterpretIssue, StructureInterpretationReport};

/// A tag whose column is resolved once per table, for reads in row loops.
#[derive(Debug, Clone, Copy)]
pub(super) struct Field {
    tag: &'static str,
    column: Option<usize>,
}

impl Field {
    pub(super) fn resolve(table: &MmcifCategory<'_>, tag: &'static str) -> Self {
        Self {
            tag,
            column: table.column(tag),
        }
    }
}

/// How a helper finds one cell: by tag, or by a pre-resolved [`Field`].
pub(super) trait FieldKey: Copy {
    fn tag(&self) -> &str;
    fn lookup<'a>(&self, table: &MmcifCategory<'a>, row: usize) -> Option<MmcifValue<'a>>;
}

impl FieldKey for &str {
    fn tag(&self) -> &str {
        self
    }

    fn lookup<'a>(&self, table: &MmcifCategory<'a>, row: usize) -> Option<MmcifValue<'a>> {
        table.value(row, self)
    }
}

impl FieldKey for &String {
    fn tag(&self) -> &str {
        self
    }

    fn lookup<'a>(&self, table: &MmcifCategory<'a>, row: usize) -> Option<MmcifValue<'a>> {
        table.value(row, self)
    }
}

impl FieldKey for Field {
    fn tag(&self) -> &str {
        self.tag
    }

    fn lookup<'a>(&self, table: &MmcifCategory<'a>, row: usize) -> Option<MmcifValue<'a>> {
        table.value_at(row, self.column?)
    }
}

/// Every `_atom_site` column read per row, resolved once.
struct AtomSiteFields {
    type_symbol: Field,
    label_asym_id: Field,
    auth_asym_id: Field,
    label_comp_id: Field,
    auth_comp_id: Field,
    label_atom_id: Field,
    auth_atom_id: Field,
    model_num: Field,
    label_entity_id: Field,
    group_pdb: Field,
    label_seq_id: Field,
    auth_seq_id: Field,
    insertion_code: Field,
    label_alt_id: Field,
    formal_charge: Field,
    id: Field,
    occupancy: Field,
    b_factor: Field,
    cartn: [Field; 3],
}

impl AtomSiteFields {
    fn resolve(table: &MmcifCategory<'_>) -> Self {
        let field = |tag| Field::resolve(table, tag);
        Self {
            type_symbol: field("_atom_site.type_symbol"),
            label_asym_id: field("_atom_site.label_asym_id"),
            auth_asym_id: field("_atom_site.auth_asym_id"),
            label_comp_id: field("_atom_site.label_comp_id"),
            auth_comp_id: field("_atom_site.auth_comp_id"),
            label_atom_id: field("_atom_site.label_atom_id"),
            auth_atom_id: field("_atom_site.auth_atom_id"),
            model_num: field("_atom_site.pdbx_PDB_model_num"),
            label_entity_id: field("_atom_site.label_entity_id"),
            group_pdb: field("_atom_site.group_PDB"),
            label_seq_id: field("_atom_site.label_seq_id"),
            auth_seq_id: field("_atom_site.auth_seq_id"),
            insertion_code: field("_atom_site.pdbx_PDB_ins_code"),
            label_alt_id: field("_atom_site.label_alt_id"),
            formal_charge: field("_atom_site.pdbx_formal_charge"),
            id: field("_atom_site.id"),
            occupancy: field("_atom_site.occupancy"),
            b_factor: field("_atom_site.B_iso_or_equiv"),
            cartn: [
                field("_atom_site.Cartn_x"),
                field("_atom_site.Cartn_y"),
                field("_atom_site.Cartn_z"),
            ],
        }
    }
}

pub(super) fn coordinate_model_ids(block: &MmcifBlock) -> Result<Vec<String>, FormatError> {
    let table = block.category("_atom_site")?.ok_or_else(|| {
        structure_error(
            DiagnosticKind::Inconsistent,
            None,
            "data block has no atom-site category",
        )
    })?;
    let table = &table;
    let model_num = Field::resolve(table, "_atom_site.pdbx_PDB_model_num");
    let mut seen = BTreeSet::new();
    let mut models = Vec::new();
    for row in 0..table.row_count() {
        let model = optional(table, row, model_num).unwrap_or("1").to_owned();
        if seen.insert(model.clone()) {
            models.push(model);
        }
    }
    Ok(models)
}

pub(super) fn read_entity_types(
    block: &MmcifBlock,
) -> Result<BTreeMap<String, EntityKind>, FormatError> {
    let mut entities = BTreeMap::new();
    if let Some(table) = block.category("_entity")? {
        let table = &table;
        for row in 0..table.row_count() {
            let id = required(table, row, "_entity.id")?;
            let kind = required(table, row, "_entity.type")?;
            if entities
                .insert(id.to_owned(), EntityKind::from_mmcif(kind))
                .is_some()
            {
                return Err(row_error(
                    table,
                    row,
                    DiagnosticKind::Inconsistent,
                    format!("duplicate entity `{id}`"),
                ));
            }
        }
    }
    Ok(entities)
}

pub(super) fn read_asym_entities(
    block: &MmcifBlock,
) -> Result<BTreeMap<String, String>, FormatError> {
    let mut instances = BTreeMap::new();
    if let Some(table) = block.category("_struct_asym")? {
        let table = &table;
        for row in 0..table.row_count() {
            let id = required(table, row, "_struct_asym.id")?;
            let entity = required(table, row, "_struct_asym.entity_id")?;
            if instances.insert(id.to_owned(), entity.to_owned()).is_some() {
                return Err(row_error(
                    table,
                    row,
                    DiagnosticKind::Inconsistent,
                    format!("duplicate structural instance `{id}`"),
                ));
            }
        }
    }
    Ok(instances)
}

#[derive(Debug, Clone)]
pub(super) struct AtomRow {
    pub(super) line: usize,
    pub(super) row_index: usize,
    pub(super) site_order: usize,
    pub(super) model_id: String,
    pub(super) entity_id: Option<String>,
    pub(super) kind: EntityKind,
    pub(super) instance_key: String,
    pub(super) label_asym_id: Option<String>,
    pub(super) asym_id: String,
    pub(super) auth_asym_id: Option<String>,
    pub(super) residue_key: String,
    pub(super) label_seq_id: Option<i32>,
    /// The residue's position along its polymer chain: `label_seq_id` for
    /// mmCIF, the residue's order within its chain for PDB.
    pub(super) polymer_position: Option<i32>,
    pub(super) auth_seq_id: Option<String>,
    pub(super) insertion_code: Option<String>,
    pub(super) occurrence: Option<usize>,
    pub(super) label_comp_id: Option<String>,
    pub(super) comp_id: String,
    pub(super) auth_comp_id: Option<String>,
    pub(super) label_atom_name: Option<String>,
    pub(super) atom_name: String,
    pub(super) auth_atom_name: Option<String>,
    pub(super) atom_site_id: Option<String>,
    pub(super) alt_id: Option<String>,
    pub(super) occupancy: Option<f64>,
    pub(super) b_factor: Option<f64>,
    pub(super) point: Option<Point3>,
    pub(super) element: Element,
    /// The mass number a `D` or `T` symbol states.
    pub(super) isotope: Option<u16>,
    /// The stated `pdbx_formal_charge`, if any.
    pub(super) formal_charge: Option<i8>,
}

impl AtomRow {
    pub(super) fn atom_key(&self) -> String {
        format!("{}|{}|{}", self.asym_id, self.residue_key, self.atom_name)
    }
}

#[derive(Debug, Default)]
struct OccurrenceState {
    occurrence: usize,
    seen: BTreeMap<String, BTreeSet<Option<String>>>,
}

pub(super) fn read_atom_rows(
    table: &MmcifCategory<'_>,
    entities: &BTreeMap<String, EntityKind>,
    asym_entities: &BTreeMap<String, String>,
    strict_entity_metadata: bool,
    report: &mut StructureInterpretationReport,
    cancellation: Option<&Cancellation>,
) -> Result<Vec<AtomRow>, FormatError> {
    let mut rows = Vec::with_capacity(table.row_count());
    let mut occurrences = BTreeMap::<(String, String, String), OccurrenceState>::new();
    let mut inferred = BTreeSet::new();
    let fields = AtomSiteFields::resolve(table);
    for row in 0..table.row_count() {
        poll_every(cancellation, row)?;
        let type_symbol = required(table, row, fields.type_symbol)?;
        let type_value = fields.type_symbol.lookup(table, row).expect("required");
        let (element, isotope) = parse_structure_element_symbol(type_symbol).ok_or_else(|| {
            structure_error(
                DiagnosticKind::InvalidValue,
                Some(type_value.line()),
                format!("unknown atom-site element `{type_symbol}`"),
            )
        })?;
        let label_asym_id = optional(table, row, fields.label_asym_id).map(str::to_owned);
        let asym_id = label_asym_id
            .as_deref()
            .or_else(|| optional(table, row, fields.auth_asym_id))
            .ok_or_else(|| {
                row_error(
                    table,
                    row,
                    DiagnosticKind::InvalidValue,
                    "missing atom-site chain identifier",
                )
            })?
            .to_owned();
        let auth_asym_id = optional(table, row, fields.auth_asym_id).map(str::to_owned);
        let label_comp_id = optional(table, row, fields.label_comp_id).map(str::to_owned);
        let comp_id = label_comp_id
            .as_deref()
            .or_else(|| optional(table, row, fields.auth_comp_id))
            .ok_or_else(|| {
                row_error(
                    table,
                    row,
                    DiagnosticKind::InvalidValue,
                    "missing atom-site component identifier",
                )
            })?
            .to_owned();
        let label_atom_name = optional(table, row, fields.label_atom_id).map(str::to_owned);
        let atom_name = label_atom_name
            .as_deref()
            .or_else(|| optional(table, row, fields.auth_atom_id))
            .ok_or_else(|| {
                row_error(
                    table,
                    row,
                    DiagnosticKind::InvalidValue,
                    "missing atom-site atom identifier",
                )
            })?
            .to_owned();
        let model_id = optional(table, row, fields.model_num)
            .unwrap_or("1")
            .to_owned();
        let entity_id = optional(table, row, fields.label_entity_id)
            .map(str::to_owned)
            .or_else(|| asym_entities.get(&asym_id).cloned());
        let group_pdb = optional(table, row, fields.group_pdb).map(str::to_owned);
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
            if strict_entity_metadata {
                return Err(row_error(
                    table,
                    row,
                    DiagnosticKind::InvalidValue,
                    format!("missing entity type for structural instance `{asym_id}`"),
                ));
            }
            if inferred.insert(asym_id.clone()) {
                report
                    .issues
                    .push(StructureInterpretIssue::EntityTypeInferred {
                        asym_id: asym_id.clone(),
                        kind: kind.clone(),
                    });
            }
        }
        let label_seq_id = optional_i32(table, row, fields.label_seq_id)?;
        let auth_seq_id = optional(table, row, fields.auth_seq_id).map(str::to_owned);
        let insertion_code = optional(table, row, fields.insertion_code).map(str::to_owned);
        let alt_id = optional(table, row, fields.label_alt_id).map(str::to_owned);
        let (residue_key, occurrence) = if let Some(sequence) = label_seq_id {
            (
                format!(
                    "label:{sequence}:{}",
                    insertion_code.as_deref().unwrap_or("")
                ),
                None,
            )
        } else if let Some(sequence) = &auth_seq_id {
            (
                format!(
                    "auth:{sequence}:{}",
                    insertion_code.as_deref().unwrap_or("")
                ),
                None,
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
            (
                format!("occurrence:{comp_id}:{}", state.occurrence),
                Some(state.occurrence),
            )
        };
        let instance_key = if kind.is_macro() {
            format!("macro:{asym_id}")
        } else {
            format!("small:{asym_id}:{residue_key}")
        };
        let formal_charge = optional_i8(table, row, fields.formal_charge)?;
        let point = optional_point(table, row, fields.cartn)?;
        rows.push(AtomRow {
            line: type_value.line(),
            row_index: row,
            site_order: row,
            model_id,
            entity_id,
            kind,
            instance_key,
            label_asym_id,
            asym_id,
            auth_asym_id,
            residue_key,
            label_seq_id,
            polymer_position: label_seq_id,
            auth_seq_id,
            insertion_code,
            occurrence,
            label_comp_id,
            comp_id,
            auth_comp_id: optional(table, row, fields.auth_comp_id).map(str::to_owned),
            label_atom_name,
            atom_name,
            auth_atom_name: optional(table, row, fields.auth_atom_id).map(str::to_owned),
            atom_site_id: optional(table, row, fields.id).map(str::to_owned),
            alt_id,
            occupancy: optional_f64(table, row, fields.occupancy)?,
            b_factor: optional_f64(table, row, fields.b_factor)?,
            point,
            element,
            isotope,
            formal_charge,
        });
    }
    Ok(rows)
}

pub(super) fn is_water(comp_id: &str) -> bool {
    ["HOH", "WAT", "DOD"]
        .iter()
        .any(|water| comp_id.eq_ignore_ascii_case(water))
}

fn infer_entity_kind(group_pdb: Option<&str>, comp_id: &str) -> EntityKind {
    if is_water(comp_id) {
        EntityKind::Water
    } else if group_pdb.is_some_and(|group| group.eq_ignore_ascii_case("ATOM")) {
        EntityKind::Polymer
    } else {
        EntityKind::NonPolymer
    }
}

pub(super) fn required<'a>(
    table: &'a MmcifCategory<'_>,
    row: usize,
    key: impl FieldKey,
) -> Result<&'a str, FormatError> {
    let tag = key.tag();
    let value = key.lookup(table, row).ok_or_else(|| {
        row_error(
            table,
            row,
            DiagnosticKind::InvalidValue,
            format!("missing required {tag}"),
        )
    })?;
    value.optional_text().ok_or_else(|| {
        structure_error(
            DiagnosticKind::InvalidValue,
            Some(value.line()),
            format!("missing required {tag}"),
        )
    })
}

pub(super) fn optional<'a>(
    table: &'a MmcifCategory<'_>,
    row: usize,
    key: impl FieldKey,
) -> Option<&'a str> {
    key.lookup(table, row)
        .and_then(|value| value.optional_text())
}

fn optional_f64(
    table: &MmcifCategory<'_>,
    row: usize,
    key: impl FieldKey,
) -> Result<Option<f64>, FormatError> {
    let tag = key.tag();
    optional(table, row, key)
        .map(|value| {
            let parsed = value.parse::<f64>().map_err(|_| {
                row_error(
                    table,
                    row,
                    DiagnosticKind::InvalidValue,
                    format!("invalid float {tag}"),
                )
            })?;
            if !parsed.is_finite() {
                return Err(row_error(
                    table,
                    row,
                    DiagnosticKind::InvalidValue,
                    format!("non-finite float {tag}"),
                ));
            }
            Ok(parsed)
        })
        .transpose()
}

pub(super) fn optional_i32(
    table: &MmcifCategory<'_>,
    row: usize,
    key: impl FieldKey,
) -> Result<Option<i32>, FormatError> {
    let tag = key.tag();
    optional(table, row, key)
        .map(|value| {
            value.parse::<i32>().map_err(|_| {
                row_error(
                    table,
                    row,
                    DiagnosticKind::InvalidValue,
                    format!("invalid integer {tag}"),
                )
            })
        })
        .transpose()
}

fn optional_i8(
    table: &MmcifCategory<'_>,
    row: usize,
    key: impl FieldKey,
) -> Result<Option<i8>, FormatError> {
    let tag = key.tag();
    optional(table, row, key)
        .map(|value| {
            value.parse::<i8>().map_err(|_| {
                row_error(
                    table,
                    row,
                    DiagnosticKind::InvalidValue,
                    format!("invalid integer {tag}"),
                )
            })
        })
        .transpose()
}

fn optional_point(
    table: &MmcifCategory<'_>,
    row: usize,
    [x, y, z]: [Field; 3],
) -> Result<Option<Point3>, FormatError> {
    let x = optional_f64(table, row, x)?;
    let y = optional_f64(table, row, y)?;
    let z = optional_f64(table, row, z)?;
    match (x, y, z) {
        (Some(x), Some(y), Some(z)) => Ok(Some(Point3::new(x, y, z))),
        (None, None, None) => Ok(None),
        _ => Err(row_error(
            table,
            row,
            DiagnosticKind::InvalidValue,
            "partial atom-site coordinate triplet",
        )),
    }
}

pub(super) fn row_error(
    table: &MmcifCategory<'_>,
    row: usize,
    kind: DiagnosticKind,
    message: impl Into<String>,
) -> FormatError {
    structure_error(kind, table.row_line(row), message)
}
