use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use crate::core::Element;
use crate::geometry::Point3;

use super::super::{MmcifDataBlock, MmcifLoopTable, MmcifValue};
use super::types::{
    MmcifAltLocPolicy, MmcifEntityKind, MmcifInterpretError, MmcifInterpretIssue,
    MmcifInterpretOptions, MmcifInterpretationReport, MmcifModelSelection,
};

pub(super) fn coordinate_model_ids(
    block: &MmcifDataBlock,
) -> Result<Vec<String>, MmcifInterpretError> {
    let table = block
        .loop_with_tag("_atom_site.type_symbol")
        .ok_or_else(|| MmcifInterpretError::new(None, "data block has no atom-site loop"))?;
    let mut seen = BTreeSet::new();
    let mut models = Vec::new();
    for row in 0..table.row_count() {
        let model = optional(table, row, "_atom_site.pdbx_PDB_model_num")
            .unwrap_or("1")
            .to_owned();
        if seen.insert(model.clone()) {
            models.push(model);
        }
    }
    Ok(models)
}

pub(super) fn read_entity_types(
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

pub(super) fn read_asym_entities(
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
pub(super) struct AtomRow {
    pub(super) line: usize,
    row_index: usize,
    pub(super) model_id: String,
    pub(super) entity_id: Option<String>,
    pub(super) kind: MmcifEntityKind,
    pub(super) instance_key: String,
    pub(super) label_asym_id: Option<String>,
    pub(super) asym_id: String,
    pub(super) auth_asym_id: Option<String>,
    pub(super) residue_key: String,
    pub(super) label_seq_id: Option<i32>,
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
    pub(super) formal_charge: i8,
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
        let label_asym_id = optional(table, row, "_atom_site.label_asym_id").map(str::to_owned);
        let asym_id = label_asym_id
            .as_deref()
            .or_else(|| optional(table, row, "_atom_site.auth_asym_id"))
            .ok_or_else(|| row_error(table, row, "missing atom-site chain identifier"))?
            .to_owned();
        let auth_asym_id = optional(table, row, "_atom_site.auth_asym_id").map(str::to_owned);
        let label_comp_id = optional(table, row, "_atom_site.label_comp_id").map(str::to_owned);
        let comp_id = label_comp_id
            .as_deref()
            .or_else(|| optional(table, row, "_atom_site.auth_comp_id"))
            .ok_or_else(|| row_error(table, row, "missing atom-site component identifier"))?
            .to_owned();
        let label_atom_name = optional(table, row, "_atom_site.label_atom_id").map(str::to_owned);
        let atom_name = label_atom_name
            .as_deref()
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
                format!("occurrence:{}", state.occurrence),
                Some(state.occurrence),
            )
        };
        let instance_key = if kind.is_macro() {
            format!("macro:{asym_id}")
        } else {
            format!("small:{asym_id}:{residue_key}")
        };
        let formal_charge =
            optional_i8(table, row, "_atom_site.pdbx_formal_charge")?.unwrap_or_default();
        let point = optional_point(table, row)?;
        rows.push(AtomRow {
            line: type_value.line(),
            row_index: row,
            model_id,
            entity_id,
            kind,
            instance_key,
            label_asym_id,
            asym_id,
            auth_asym_id,
            residue_key,
            label_seq_id,
            auth_seq_id,
            insertion_code,
            occurrence,
            label_comp_id,
            comp_id,
            auth_comp_id: optional(table, row, "_atom_site.auth_comp_id").map(str::to_owned),
            label_atom_name,
            atom_name,
            auth_atom_name: optional(table, row, "_atom_site.auth_atom_id").map(str::to_owned),
            atom_site_id: optional(table, row, "_atom_site.id").map(str::to_owned),
            alt_id,
            occupancy: optional_f64(table, row, "_atom_site.occupancy")?,
            b_factor: optional_f64(table, row, "_atom_site.B_iso_or_equiv")?,
            point,
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

pub(super) fn select_alt_locations(
    rows: &[AtomRow],
    policy: &MmcifAltLocPolicy,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<AtomRow>, MmcifInterpretError> {
    let mut grouped = BTreeMap::<(String, String, String), Vec<&AtomRow>>::new();
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
            MmcifAltLocPolicy::HighestOccupancy => candidates
                .iter()
                .max_by(|left, right| {
                    left.occupancy
                        .unwrap_or(0.0)
                        .partial_cmp(&right.occupancy.unwrap_or(0.0))
                        .unwrap_or(Ordering::Equal)
                        .then_with(|| right.alt_id.cmp(&left.alt_id))
                })
                .map(|row| (**row).clone()),
            MmcifAltLocPolicy::SelectLabel(label) => candidates
                .iter()
                .find(|row| row.alt_id.as_deref() == Some(label.as_str()))
                .map(|row| (**row).clone())
                .or_else(|| {
                    candidates
                        .iter()
                        .find(|row| row.alt_id.is_none())
                        .map(|row| (**row).clone())
                }),
            MmcifAltLocPolicy::ErrorOnAlternateLocations => {
                candidates.first().map(|row| (**row).clone())
            }
        };
        let Some(chosen) = chosen else {
            return Err(MmcifInterpretError::new(
                None,
                "requested alternate-location label is unavailable",
            ));
        };
        for omitted in candidates
            .iter()
            .filter(|candidate| candidate.row_index != chosen.row_index)
        {
            report
                .issues
                .push(MmcifInterpretIssue::AlternateLocationOmitted {
                    atom_name: omitted.atom_name.clone(),
                    alt_id: omitted.alt_id.clone(),
                });
        }
        selected.push(chosen);
    }
    selected.sort_by_key(|row| row.row_index);
    Ok(selected)
}

pub(super) fn select_coordinate_model(
    rows: Vec<AtomRow>,
    selection: &MmcifModelSelection,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<AtomRow>, MmcifInterpretError> {
    let mut model_ids = Vec::new();
    let mut seen = BTreeSet::new();
    for row in &rows {
        if seen.insert(row.model_id.clone()) {
            model_ids.push(row.model_id.clone());
        }
    }
    report.coordinate_models = model_ids.len();
    let selected = match selection {
        MmcifModelSelection::RequireSingle if model_ids.len() == 1 => model_ids[0].clone(),
        MmcifModelSelection::RequireSingle => {
            return Err(MmcifInterpretError::new(
                None,
                format!(
                    "coordinate data contains {} models; select one explicitly",
                    model_ids.len()
                ),
            ));
        }
        MmcifModelSelection::Select(id) if seen.contains(id) => id.clone(),
        MmcifModelSelection::Select(id) => {
            return Err(MmcifInterpretError::new(
                None,
                format!("coordinate model `{id}` is unavailable"),
            ));
        }
        MmcifModelSelection::First => model_ids
            .first()
            .cloned()
            .ok_or_else(|| MmcifInterpretError::new(None, "coordinate data contains no models"))?,
    };
    report.selected_model = Some(selected.clone());
    for ignored in model_ids.iter().filter(|id| **id != selected) {
        let atom_site_rows = rows.iter().filter(|row| row.model_id == *ignored).count();
        report.ignored_coordinate_models.push(ignored.clone());
        report
            .issues
            .push(MmcifInterpretIssue::CoordinateModelIgnored {
                model_id: ignored.clone(),
                atom_site_rows,
            });
    }
    let selected_rows = rows
        .into_iter()
        .filter(|row| row.model_id == selected)
        .collect::<Vec<_>>();
    if let Some(row) = selected_rows.iter().find(|row| row.point.is_none()) {
        return Err(MmcifInterpretError::new(
            Some(row.line),
            format!(
                "selected coordinate model `{selected}` has no complete position for atom `{}`",
                row.atom_name
            ),
        ));
    }
    Ok(selected_rows)
}

pub(super) fn required<'a>(
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

pub(super) fn optional<'a>(table: &'a MmcifLoopTable, row: usize, tag: &str) -> Option<&'a str> {
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

pub(super) fn optional_i32(
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

pub(super) fn row_error(
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
