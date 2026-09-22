use std::collections::{BTreeMap, BTreeSet};

use super::super::atom_site::{read_asym_entities, read_atom_rows, read_entity_types, AtomRow};
use super::super::{MmcifInterpretError, MmcifInterpretIssue, MmcifInterpretationReport};
use super::{
    MmcifAltLocDecision, MmcifAltLocPolicy, MmcifAltLocPreference, MmcifAltLocResidue,
    MmcifAltLocSelection, MmcifAltLocSelectionReason, MmcifResidueId,
};
use crate::io::MmcifBlock;

struct Residue<'a> {
    id: MmcifResidueId,
    rows: Vec<&'a AtomRow>,
    labels: BTreeSet<String>,
}

fn residues(rows: &[AtomRow]) -> Result<Vec<Residue<'_>>, MmcifInterpretError> {
    let mut grouped = BTreeMap::<MmcifResidueId, Vec<&AtomRow>>::new();
    for row in rows {
        grouped
            .entry(MmcifResidueId::of(row))
            .or_default()
            .push(row);
    }
    grouped
        .into_iter()
        .map(|(id, rows)| {
            let mut identities = BTreeSet::new();
            for row in &rows {
                if !identities.insert((&row.comp_id, &row.atom_name, &row.alt_id)) {
                    return Err(MmcifInterpretError::new(
                        Some(row.line),
                        format!(
                            "atom `{}` has duplicate records for one alternate location in {id:?}",
                            row.atom_name
                        ),
                    ));
                }
            }
            let labels = rows.iter().filter_map(|row| row.alt_id.clone()).collect();
            Ok(Residue { id, rows, labels })
        })
        .collect()
}

pub(crate) fn inventory(
    block: &MmcifBlock,
) -> Result<Vec<MmcifAltLocResidue>, MmcifInterpretError> {
    let entities = read_entity_types(block)?;
    let asym = read_asym_entities(block)?;
    let table = block
        .category("_atom_site")?
        .ok_or_else(|| MmcifInterpretError::new(None, "block has no atom-site category"))?;
    let rows = read_atom_rows(
        &table,
        &entities,
        &asym,
        false,
        &mut MmcifInterpretationReport::default(),
    )?;
    super::source::conformations(block, &rows)?;
    Ok(residues(&rows)?
        .into_iter()
        .filter(|residue| !residue.labels.is_empty())
        .map(|residue| MmcifAltLocResidue {
            residue: residue.id,
            labels: residue.labels.into_iter().collect(),
        })
        .collect())
}

fn configuration(policy: &MmcifAltLocPolicy) -> MmcifAltLocSelection {
    let default = match policy {
        MmcifAltLocPolicy::HighestOccupancy => MmcifAltLocPreference::HighestOccupancy,
        MmcifAltLocPolicy::SelectLabel(label) => MmcifAltLocPreference::SelectLabel(label.clone()),
        MmcifAltLocPolicy::PreferLabel(label) => MmcifAltLocPreference::PreferLabel(label.clone()),
        MmcifAltLocPolicy::ErrorOnAlternateLocations => {
            MmcifAltLocPreference::ErrorOnAlternateLocations
        }
        MmcifAltLocPolicy::Configured(selection) => return selection.clone(),
    };
    MmcifAltLocSelection {
        default,
        ..Default::default()
    }
}

fn selection_error(residue: &Residue<'_>, message: impl std::fmt::Display) -> MmcifInterpretError {
    MmcifInterpretError::new(
        Some(residue.rows[0].line),
        format!("{message} at {:?}", residue.id),
    )
}

/// A source label can cover part of a residue. Unlabelled rows are shared, but
/// atoms belonging to a different chemical component are never borrowed.
fn candidate_rows<'a>(
    residue: &Residue<'a>,
    labels: &BTreeSet<String>,
) -> Result<Vec<&'a AtomRow>, MmcifInterpretError> {
    let labelled = residue
        .rows
        .iter()
        .copied()
        .filter(|row| {
            row.alt_id
                .as_ref()
                .is_some_and(|label| labels.contains(label))
        })
        .collect::<Vec<_>>();
    if labelled.is_empty() && !residue.labels.is_empty() {
        return Ok(Vec::new());
    }
    let components = if labelled.is_empty() {
        &residue.rows
    } else {
        &labelled
    }
    .iter()
    .map(|row| row.comp_id.clone())
    .collect::<BTreeSet<_>>();
    if components.len() != 1 {
        return Err(selection_error(
            residue,
            "canonical residue identity has conflicting chemical components in alternate selection",
        ));
    }
    let mut selected = labelled;
    selected.extend(
        residue
            .rows
            .iter()
            .copied()
            .filter(|row| row.alt_id.is_none() && components.contains(row.comp_id.as_str())),
    );
    let mut atoms = BTreeSet::new();
    for row in &selected {
        if !atoms.insert(row.atom_name.as_str()) {
            return Err(selection_error(
                residue,
                format!(
                    "alternate selection contains overlapping sites for atom `{}`",
                    row.atom_name
                ),
            ));
        }
    }
    Ok(selected)
}

fn missing_atoms(residue: &Residue<'_>, selected: &[&AtomRow]) -> Vec<String> {
    let Some(first) = selected.first() else {
        return Vec::new();
    };
    let present = selected
        .iter()
        .map(|row| row.atom_name.as_str())
        .collect::<BTreeSet<_>>();
    residue
        .rows
        .iter()
        .filter(|row| row.comp_id == first.comp_id && !present.contains(row.atom_name.as_str()))
        .map(|row| row.atom_name.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn score(rows: &[&AtomRow]) -> f64 {
    let mut labelled = rows
        .iter()
        .filter(|row| row.alt_id.is_some())
        .collect::<Vec<_>>();
    if labelled.is_empty() {
        return 0.0;
    }
    labelled.sort_by(|a, b| a.atom_name.cmp(&b.atom_name));
    // Divide before summing to keep even large finite source values bounded.
    labelled
        .iter()
        .map(|row| row.occupancy.unwrap_or(0.0) / labelled.len() as f64)
        .sum()
}

fn labels_of(rows: &[&AtomRow]) -> BTreeSet<String> {
    rows.iter().filter_map(|row| row.alt_id.clone()).collect()
}

fn root(parents: &[usize], mut index: usize) -> usize {
    while parents[index] != index {
        index = parents[index];
    }
    index
}

fn groups(
    residues: &[Residue<'_>],
    selection: &MmcifAltLocSelection,
    source_declared: bool,
) -> Result<Vec<Vec<usize>>, MmcifInterpretError> {
    let indices = residues
        .iter()
        .enumerate()
        .map(|(index, residue)| (&residue.id, index))
        .collect::<BTreeMap<_, _>>();
    let mut parents = (0..residues.len()).collect::<Vec<_>>();
    for id in selection
        .overrides
        .keys()
        .chain(selection.groups.iter().flatten())
    {
        if !indices.contains_key(id) || residues[indices[id]].labels.is_empty() {
            return Err(MmcifInterpretError::new(
                None,
                format!(
                    "alternate-location selection refers to unknown or unlabelled residue {id:?}"
                ),
            ));
        }
    }
    for group in &selection.groups {
        if group.len() < 2 || group.iter().collect::<BTreeSet<_>>().len() != group.len() {
            return Err(MmcifInterpretError::new(
                None,
                "alternate-location group requires at least two distinct residues",
            ));
        }
        if group.iter().any(|id| id.model_id != group[0].model_id) {
            return Err(MmcifInterpretError::new(
                None,
                "alternate-location group crosses source coordinate models",
            ));
        }
        let first = indices[&group[0]];
        for id in &group[1..] {
            let a = root(&parents, first);
            let b = root(&parents, indices[id]);
            parents[b] = a;
        }
    }
    if source_declared {
        let mut models = BTreeMap::new();
        for (index, residue) in residues
            .iter()
            .enumerate()
            .filter(|(_, residue)| !residue.labels.is_empty())
        {
            let first = *models.entry(&residue.id.model_id).or_insert(index);
            let a = root(&parents, first);
            let b = root(&parents, index);
            parents[b] = a;
        }
    }
    let mut groups = BTreeMap::<usize, Vec<usize>>::new();
    for (index, residue) in residues.iter().enumerate() {
        if !residue.labels.is_empty() {
            groups.entry(root(&parents, index)).or_default().push(index);
        }
    }
    Ok(groups.into_values().collect())
}

struct Candidate<'a> {
    id: String,
    rows: Vec<Vec<&'a AtomRow>>,
    score: f64,
    preferred: usize,
}

pub(in crate::io::mmcif_interpret) fn select_alt_locations(
    block: &MmcifBlock,
    rows: &[AtomRow],
    policy: &MmcifAltLocPolicy,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<AtomRow>, MmcifInterpretError> {
    let residues = residues(rows)?;
    let selection = configuration(policy);
    let conformations = super::source::conformations(block, rows)?;
    if let Some(id) = &selection.source_conformation {
        if !conformations.contains_key(id) {
            return Err(MmcifInterpretError::new(
                None,
                format!("source alternate conformation `{id}` is unavailable"),
            ));
        }
    }
    let grouped = groups(&residues, &selection, !conformations.is_empty())?;
    let mut selected = Vec::new();
    for residue in &residues {
        if residue.labels.is_empty() {
            selected.extend(candidate_rows(residue, &BTreeSet::new())?);
        } else {
            if matches!(
                selection.default,
                MmcifAltLocPreference::ErrorOnAlternateLocations
            ) && !selection.overrides.contains_key(&residue.id)
            {
                return Err(selection_error(residue, "residue has alternate locations"));
            }
            let missing = residue
                .rows
                .iter()
                .filter(|row| row.alt_id.is_some() && row.occupancy.is_none())
                .map(|row| row.line)
                .collect::<Vec<_>>();
            if !missing.is_empty() {
                report
                    .issues
                    .push(MmcifInterpretIssue::AlternateLocationOccupancyMissing {
                        residue: residue.id.clone(),
                        source_lines: missing,
                    });
            }
        }
    }
    for group in grouped {
        let alternatives = if conformations.is_empty() {
            residues[group[0]]
                .labels
                .iter()
                .map(|label| (label.clone(), BTreeSet::from([label.clone()])))
                .collect()
        } else {
            conformations.clone()
        };
        let chosen = choose(
            &residues,
            &group,
            &selection,
            &alternatives,
            !conformations.is_empty(),
            report,
        )?;
        for (index, retained) in group.iter().zip(chosen.rows) {
            let residue = &residues[*index];
            let labels = labels_of(&retained);
            let kept = retained
                .iter()
                .map(|row| row.row_index)
                .collect::<BTreeSet<_>>();
            let omitted = residue
                .rows
                .iter()
                .filter(|row| !kept.contains(&row.row_index))
                .collect::<Vec<_>>();
            for row in &omitted {
                report
                    .issues
                    .push(MmcifInterpretIssue::AlternateLocationOmitted {
                        atom_name: row.atom_name.clone(),
                        alt_id: row.alt_id.clone(),
                        residue: residue.id.clone(),
                        source_line: row.line,
                    });
            }
            let reason = if !conformations.is_empty() {
                MmcifAltLocSelectionReason::SourceConformation {
                    id: chosen.id.clone(),
                }
            } else if selection.overrides.contains_key(&residue.id)
                || matches!(selection.default, MmcifAltLocPreference::SelectLabel(_))
            {
                MmcifAltLocSelectionReason::ExactSelection
            } else if let MmcifAltLocPreference::PreferLabel(label) = &selection.default {
                if labels.contains(label) {
                    MmcifAltLocSelectionReason::PreferredLabel
                } else {
                    MmcifAltLocSelectionReason::PreferredLabelUnavailable
                }
            } else {
                MmcifAltLocSelectionReason::HighestOccupancy
            };
            report.alternate_locations.push(MmcifAltLocDecision {
                residue: residue.id.clone(),
                available_labels: residue.labels.iter().cloned().collect(),
                selected_labels: labels.into_iter().collect(),
                reason,
                selected_source_lines: retained.iter().map(|row| row.line).collect(),
                omitted_source_lines: omitted.iter().map(|row| row.line).collect(),
            });
            selected.extend(retained);
        }
    }
    // Dense atom order follows each site's first source occurrence, independently
    // of which alternative's row supplied its coordinates.
    let mut order = BTreeMap::new();
    for row in rows {
        order
            .entry((MmcifResidueId::of(row), &row.comp_id, &row.atom_name))
            .or_insert(row.row_index);
    }
    let mut selected = selected
        .into_iter()
        .map(|row| {
            let mut chosen = row.clone();
            chosen.site_order = order[&(MmcifResidueId::of(row), &row.comp_id, &row.atom_name)];
            chosen
        })
        .collect::<Vec<_>>();
    selected.sort_by_key(|row| row.site_order);
    Ok(selected)
}

fn choose<'a>(
    residues: &[Residue<'a>],
    group: &[usize],
    selection: &MmcifAltLocSelection,
    alternatives: &BTreeMap<String, BTreeSet<String>>,
    source_declared: bool,
    report: &mut MmcifInterpretationReport,
) -> Result<Candidate<'a>, MmcifInterpretError> {
    let mut candidates = Vec::new();
    let mut last_error = None;
    for (id, labels) in alternatives {
        if selection
            .source_conformation
            .as_ref()
            .is_some_and(|wanted| wanted != id)
        {
            continue;
        }
        let mut candidate = Candidate {
            id: id.clone(),
            rows: Vec::new(),
            score: 0.0,
            preferred: 0,
        };
        let mut valid = true;
        for index in group {
            let residue = &residues[*index];
            let retained = match candidate_rows(residue, labels) {
                Ok(rows) => rows,
                Err(error) => {
                    if source_declared {
                        report.issues.push(
                            MmcifInterpretIssue::AlternateLocationCombinationRejected {
                                residue: residue.id.clone(),
                                alternative: id.clone(),
                                reason: error.message().to_owned(),
                            },
                        );
                    }
                    last_error = Some(error);
                    valid = false;
                    break;
                }
            };
            let missing = missing_atoms(residue, &retained);
            let incomplete = !missing.is_empty();
            if incomplete {
                report
                    .issues
                    .push(MmcifInterpretIssue::AlternateLocationIncomplete {
                        residue: residue.id.clone(),
                        alternative: id.clone(),
                        missing_atoms: missing,
                    });
            }
            if retained.is_empty() || incomplete {
                if source_declared && retained.is_empty() {
                    report
                        .issues
                        .push(MmcifInterpretIssue::AlternateLocationCombinationRejected {
                            residue: residue.id.clone(),
                            alternative: id.clone(),
                            reason: "combination contains no alternate label for this residue"
                                .into(),
                        });
                }
                valid = false;
                break;
            }
            let selected_labels = labels_of(&retained);
            let exact = selection
                .overrides
                .get(&residue.id)
                .or(match &selection.default {
                    MmcifAltLocPreference::SelectLabel(label) => Some(label),
                    _ => None,
                });
            if let Some(label) = exact {
                if selected_labels.len() != 1 || !selected_labels.contains(label) {
                    valid = false;
                    break;
                }
            } else {
                match &selection.default {
                    MmcifAltLocPreference::ErrorOnAlternateLocations => {
                        return Err(selection_error(residue, "residue has alternate locations"))
                    }
                    MmcifAltLocPreference::PreferLabel(label)
                        if selected_labels.contains(label) =>
                    {
                        candidate.preferred += 1
                    }
                    _ => {}
                }
            }
            candidate.score += score(&retained) / group.len() as f64;
            candidate.rows.push(retained);
        }
        // Source combinations must also satisfy user-declared same-label groups.
        if valid && source_declared {
            for linked in &selection.groups {
                if linked[0].model_id != residues[group[0]].id.model_id {
                    continue;
                }
                let mut common = None;
                for residue_id in linked {
                    let offset = group
                        .iter()
                        .position(|index| &residues[*index].id == residue_id)
                        .expect("validated source group");
                    let labels = labels_of(&candidate.rows[offset]);
                    if labels.len() != 1 || common.as_ref().is_some_and(|value| value != &labels) {
                        valid = false;
                        break;
                    }
                    common = Some(labels);
                }
            }
        }
        if valid {
            candidates.push(candidate);
        }
    }
    if candidates.is_empty() {
        return Err(last_error.unwrap_or_else(|| selection_error(&residues[group[0]], format!(
            "no complete alternate configuration satisfies the requested labels, residue groups, and source constraints (default: {:?}, overrides: {:?}, available: {})",
            selection.default, selection.overrides, alternatives.keys().cloned().collect::<Vec<_>>().join(", ")
        ))));
    }
    candidates.sort_by(|a, b| {
        b.preferred
            .cmp(&a.preferred)
            .then_with(|| b.score.total_cmp(&a.score))
            .then_with(|| a.id.cmp(&b.id))
    });
    let best = &candidates[0];
    let tied = candidates
        .iter()
        .filter(|candidate| candidate.preferred == best.preferred && candidate.score == best.score)
        .map(|candidate| candidate.id.clone())
        .collect::<Vec<_>>();
    if tied.len() > 1 {
        report
            .issues
            .push(MmcifInterpretIssue::AlternateLocationOccupancyTie {
                residues: group
                    .iter()
                    .map(|index| residues[*index].id.clone())
                    .collect(),
                alternatives: tied,
            });
    }
    Ok(candidates.remove(0))
}
