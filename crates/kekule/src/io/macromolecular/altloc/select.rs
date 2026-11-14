//! Alternate-location selection and inventory over typed atom rows.
use super::super::types::{poll, poll_every, structure_error};
use crate::diagnostics::{Cancellation, Capitalized, DiagnosticKind, FormatError};
use std::collections::{BTreeMap, BTreeSet};

use super::super::atom_site::AtomRow;
use super::super::{StructureInterpretIssue, StructureInterpretationReport};
use super::{
    AltLocDecision, AltLocPolicy, AltLocPreference, AltLocResidue, AltLocSelection,
    AltLocSelectionReason, SourceResidueId,
};

struct Residue<'a> {
    id: SourceResidueId,
    rows: Vec<&'a AtomRow>,
    labels: BTreeSet<String>,
}

fn residues<'a>(
    rows: &'a [AtomRow],
    cancellation: Option<&Cancellation>,
) -> Result<Vec<Residue<'a>>, FormatError> {
    // A residue's rows usually follow each other, so its identity is built
    // once per run of rows rather than once per row.
    let mut grouped = BTreeMap::<SourceResidueId, Vec<&AtomRow>>::new();
    let mut run: Option<(SourceResidueId, Vec<&AtomRow>)> = None;
    for (step, row) in rows.iter().enumerate() {
        poll_every(cancellation, step)?;
        match &mut run {
            Some((id, rows)) if id.matches(row) => rows.push(row),
            _ => {
                if let Some((id, rows)) = run.take() {
                    grouped.entry(id).or_default().extend(rows);
                }
                run = Some((SourceResidueId::of(row), vec![row]));
            }
        }
    }
    if let Some((id, rows)) = run {
        grouped.entry(id).or_default().extend(rows);
    }
    grouped
        .into_iter()
        .enumerate()
        .map(|(step, (id, rows))| {
            poll_every(cancellation, step)?;
            let mut identities = BTreeSet::new();
            for row in &rows {
                if !identities.insert((&row.comp_id, &row.atom_name, &row.alt_id)) {
                    return Err(structure_error(
                        DiagnosticKind::Inconsistent,
                        Some(row.line),
                        format!(
                            "Atom {:?} has two records for one alternate location in source conformation {id:?}.",
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

/// The residues of a source that have alternate locations.
pub(in crate::io::macromolecular) fn inventory(
    rows: &[AtomRow],
) -> Result<Vec<AltLocResidue>, FormatError> {
    Ok(residues(rows, None)?
        .into_iter()
        .filter(|residue| !residue.labels.is_empty())
        .map(|residue| AltLocResidue {
            residue: residue.id,
            labels: residue.labels.into_iter().collect(),
        })
        .collect())
}

fn configuration(policy: &AltLocPolicy) -> AltLocSelection {
    let default = match policy {
        AltLocPolicy::HighestOccupancy => AltLocPreference::HighestOccupancy,
        AltLocPolicy::SelectLabel(label) => AltLocPreference::SelectLabel(label.clone()),
        AltLocPolicy::PreferLabel(label) => AltLocPreference::PreferLabel(label.clone()),
        AltLocPolicy::ErrorOnAlternateLocations => AltLocPreference::ErrorOnAlternateLocations,
        AltLocPolicy::Configured(selection) => return selection.clone(),
    };
    AltLocSelection::default().with_default_preference(default)
}

fn selection_error(
    residue: &Residue<'_>,
    kind: DiagnosticKind,
    message: impl std::fmt::Display,
) -> FormatError {
    structure_error(
        kind,
        Some(residue.rows[0].line),
        format!("{} {message}", Capitalized(&residue.id)),
    )
}

/// A source label can cover part of a residue. Unlabelled rows are shared, but
/// atoms belonging to a different chemical component are never borrowed.
fn candidate_rows<'a>(
    residue: &Residue<'a>,
    labels: &BTreeSet<String>,
) -> Result<Vec<&'a AtomRow>, FormatError> {
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
    .map(|row| row.comp_id.as_str())
    .collect::<BTreeSet<_>>();
    if components.len() != 1 {
        return Err(selection_error(
            residue,
            DiagnosticKind::Inconsistent,
            "has different components in its selected alternate locations.",
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
                DiagnosticKind::Inconsistent,
                format!("has two selected sites for atom {:?}.", row.atom_name),
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
    selection: &AltLocSelection,
    source_declared: bool,
) -> Result<Vec<Vec<usize>>, FormatError> {
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
            return Err(structure_error(
                DiagnosticKind::Inconsistent,
                None,
                format!(
                    "The alternate-location selection names {id}, which has no alternate locations."
                ),
            ));
        }
    }
    for group in &selection.groups {
        if group.len() < 2 || group.iter().collect::<BTreeSet<_>>().len() != group.len() {
            return Err(structure_error(
                DiagnosticKind::Inconsistent,
                None,
                "An alternate-location group names fewer than two distinct residues.",
            ));
        }
        if group.iter().any(|id| id.model_id != group[0].model_id) {
            return Err(structure_error(
                DiagnosticKind::Inconsistent,
                None,
                "An alternate-location group spans more than one coordinate model.",
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

/// The selected rows as their positions in `rows` (each row's `row_index`)
/// with their sites' dense order, sorted by that order.
pub(in crate::io::macromolecular) fn select_alt_locations(
    rows: &[AtomRow],
    conformations: &BTreeMap<String, BTreeSet<String>>,
    policy: &AltLocPolicy,
    report: &mut StructureInterpretationReport,
    cancellation: Option<&Cancellation>,
) -> Result<Vec<(usize, usize)>, FormatError> {
    debug_assert!(rows
        .iter()
        .enumerate()
        .all(|(position, row)| row.row_index == position));
    let residues = residues(rows, cancellation)?;
    let selection = configuration(policy);
    if let Some(id) = &selection.source_conformation {
        if !conformations.contains_key(id) {
            return Err(structure_error(
                DiagnosticKind::Inconsistent,
                None,
                format!("Source conformation {id:?} does not exist."),
            ));
        }
    }
    let grouped = groups(&residues, &selection, !conformations.is_empty())?;
    poll(cancellation)?;
    let mut selected = Vec::new();
    for (step, residue) in residues.iter().enumerate() {
        poll_every(cancellation, step)?;
        if residue.labels.is_empty() {
            selected.extend(candidate_rows(residue, &BTreeSet::new())?);
        } else {
            if matches!(
                selection.default,
                AltLocPreference::ErrorOnAlternateLocations
            ) && !selection.overrides.contains_key(&residue.id)
            {
                return Err(selection_error(
                    residue,
                    DiagnosticKind::Unsupported,
                    "has alternate locations.",
                ));
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
                    .push(StructureInterpretIssue::AlternateLocationOccupancyMissing {
                        residue: residue.id.clone(),
                        source_lines: missing,
                    });
            }
        }
    }
    for group in grouped {
        poll(cancellation)?;
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
                    .push(StructureInterpretIssue::AlternateLocationOmitted {
                        atom_name: row.atom_name.clone(),
                        alt_id: row.alt_id.clone(),
                        residue: residue.id.clone(),
                        source_line: row.line,
                    });
            }
            let reason = if !conformations.is_empty() {
                AltLocSelectionReason::SourceConformation {
                    id: chosen.id.clone(),
                }
            } else if selection.overrides.contains_key(&residue.id)
                || matches!(selection.default, AltLocPreference::SelectLabel(_))
            {
                AltLocSelectionReason::ExactSelection
            } else if let AltLocPreference::PreferLabel(label) = &selection.default {
                if labels.contains(label) {
                    AltLocSelectionReason::PreferredLabel
                } else {
                    AltLocSelectionReason::PreferredLabelUnavailable
                }
            } else {
                AltLocSelectionReason::HighestOccupancy
            };
            report.alternate_locations.push(AltLocDecision {
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
    // of which alternative's row supplied its coordinates. Rows without
    // alternatives are their sites' only rows.
    let mut site_orders = BTreeMap::new();
    for residue in residues.iter().filter(|residue| !residue.labels.is_empty()) {
        poll(cancellation)?;
        let mut first = BTreeMap::<(&str, &str), usize>::new();
        for row in &residue.rows {
            let order = *first
                .entry((row.comp_id.as_str(), row.atom_name.as_str()))
                .or_insert(row.row_index);
            site_orders.insert(row.row_index, order);
        }
    }
    let mut selected = selected
        .into_iter()
        .enumerate()
        .map(|(step, row)| {
            poll_every(cancellation, step)?;
            let order = site_orders
                .get(&row.row_index)
                .copied()
                .unwrap_or(row.row_index);
            Ok((row.row_index, order))
        })
        .collect::<Result<Vec<_>, FormatError>>()?;
    poll(cancellation)?;
    selected.sort_by_key(|&(_, order)| order);
    Ok(selected)
}

fn choose<'a>(
    residues: &[Residue<'a>],
    group: &[usize],
    selection: &AltLocSelection,
    alternatives: &BTreeMap<String, BTreeSet<String>>,
    source_declared: bool,
    report: &mut StructureInterpretationReport,
) -> Result<Candidate<'a>, FormatError> {
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
                            StructureInterpretIssue::AlternateLocationCombinationRejected {
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
                    .push(StructureInterpretIssue::AlternateLocationIncomplete {
                        residue: residue.id.clone(),
                        alternative: id.clone(),
                        missing_atoms: missing,
                    });
            }
            if retained.is_empty() || incomplete {
                if source_declared && retained.is_empty() {
                    report.issues.push(
                        StructureInterpretIssue::AlternateLocationCombinationRejected {
                            residue: residue.id.clone(),
                            alternative: id.clone(),
                            reason: "The combination has no alternate location of this residue."
                                .into(),
                        },
                    );
                }
                valid = false;
                break;
            }
            let selected_labels = labels_of(&retained);
            let exact = selection
                .overrides
                .get(&residue.id)
                .or(match &selection.default {
                    AltLocPreference::SelectLabel(label) => Some(label),
                    _ => None,
                });
            if let Some(label) = exact {
                if selected_labels.len() != 1 || !selected_labels.contains(label) {
                    valid = false;
                    break;
                }
            } else {
                match &selection.default {
                    AltLocPreference::ErrorOnAlternateLocations => {
                        return Err(selection_error(
                            residue,
                            DiagnosticKind::Unsupported,
                            "has alternate locations.",
                        ))
                    }
                    AltLocPreference::PreferLabel(label) if selected_labels.contains(label) => {
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
        return Err(last_error.unwrap_or_else(|| selection_error(&residues[group[0]], DiagnosticKind::Inconsistent, format!(
            "has no complete alternate-location combination that meets the selection. The source conformations are {}.",
            alternatives.keys().map(|id| format!("{id:?}")).collect::<Vec<_>>().join(", ")
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
            .push(StructureInterpretIssue::AlternateLocationOccupancyTie {
                residues: group
                    .iter()
                    .map(|index| residues[*index].id.clone())
                    .collect(),
                alternatives: tied,
            });
    }
    Ok(candidates.remove(0))
}
