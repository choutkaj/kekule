use std::collections::BTreeMap;

use crate::core::BondOrder;

use super::super::{MmcifDataBlock, MmcifLoopTable, MmcifValue};
use super::atom_site::{optional, optional_i32, required, row_error, AtomRow};
use super::types::{
    MmcifConnectionResolutionReason, MmcifInterpretError, MmcifInterpretIssue,
    MmcifInterpretationReport,
};

#[derive(Debug, Clone)]
pub(super) struct DeclaredConnection {
    pub(super) left_atom: String,
    pub(super) right_atom: String,
    pub(super) order: BondOrder,
}

pub(super) fn read_connections(
    block: &MmcifDataBlock,
    rows: &[AtomRow],
    all_rows: &[AtomRow],
    selected_model: &str,
    union: &mut InstanceUnion,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<DeclaredConnection>, MmcifInterpretError> {
    let Some(table) = block.loop_with_tag("_struct_conn.conn_type_id") else {
        return Ok(Vec::new());
    };
    let mut connections = Vec::new();
    for row in 0..table.row_count() {
        let kind = required(table, row, "_struct_conn.conn_type_id")?.to_owned();
        let connection_id = optional(table, row, "_struct_conn.id").map(str::to_owned);
        let source_line = table
            .row(row)
            .and_then(|values| values.first())
            .map(MmcifValue::line);
        if !is_covalent_connection(&kind) {
            report.issues.push(MmcifInterpretIssue::ConnectionIgnored {
                connection_type: kind,
            });
            continue;
        }
        let order = connection_bond_order(table, row)?;
        let left = connection_partner(table, row, 1, rows, all_rows, selected_model)?;
        let right = connection_partner(table, row, 2, rows, all_rows, selected_model)?;
        let left = report_connection_partner_resolution(
            left,
            connection_id.as_deref(),
            &kind,
            1,
            source_line,
            report,
        );
        let right = report_connection_partner_resolution(
            right,
            connection_id.as_deref(),
            &kind,
            2,
            source_line,
            report,
        );
        let (Some(left), Some(right)) = (left, right) else {
            continue;
        };
        union.union(&left.instance_key, &right.instance_key);
        connections.push(DeclaredConnection {
            left_atom: left.atom_key(),
            right_atom: right.atom_key(),
            order,
        });
        report.applied_connections += 1;
    }
    Ok(connections)
}

fn report_connection_partner_resolution<'a>(
    resolution: ConnectionPartnerResolution<'a>,
    connection_id: Option<&str>,
    connection_type: &str,
    partner: u8,
    source_line: Option<usize>,
    report: &mut MmcifInterpretationReport,
) -> Option<&'a AtomRow> {
    match resolution {
        ConnectionPartnerResolution::Resolved(atom) => Some(atom),
        ConnectionPartnerResolution::Unresolved(reason) => {
            report
                .issues
                .push(MmcifInterpretIssue::ConnectionUnresolved {
                    connection_id: connection_id.map(str::to_owned),
                    connection_type: connection_type.to_owned(),
                    partner,
                    source_line,
                    reason,
                });
            None
        }
        ConnectionPartnerResolution::Ambiguous { candidates, reason } => {
            report
                .issues
                .push(MmcifInterpretIssue::ConnectionAmbiguous {
                    connection_id: connection_id.map(str::to_owned),
                    connection_type: connection_type.to_owned(),
                    partner,
                    source_line,
                    candidates,
                    reason,
                });
            None
        }
    }
}

fn connection_bond_order(
    table: &MmcifLoopTable,
    row: usize,
) -> Result<BondOrder, MmcifInterpretError> {
    let Some(order) = optional(table, row, "_struct_conn.pdbx_value_order") else {
        return Ok(BondOrder::Single);
    };
    match order.to_ascii_lowercase().as_str() {
        "sing" => Ok(BondOrder::Single),
        "doub" => Ok(BondOrder::Double),
        "trip" => Ok(BondOrder::Triple),
        "quad" => Ok(BondOrder::Quadruple),
        _ => Err(row_error(
            table,
            row,
            format!("unsupported struct_conn bond order `{order}`"),
        )),
    }
}

fn is_covalent_connection(kind: &str) -> bool {
    let kind = kind.to_ascii_lowercase();
    kind.starts_with("covale") || kind == "disulf" || kind == "modres"
}

#[derive(Debug, Default)]
struct LabelConnectionSelector {
    asym_id: Option<String>,
    component_id: Option<String>,
    sequence_id: Option<i32>,
    atom_id: Option<String>,
    alternate_location: Option<String>,
}

impl LabelConnectionSelector {
    fn is_empty(&self) -> bool {
        self.asym_id.is_none()
            && self.component_id.is_none()
            && self.sequence_id.is_none()
            && self.atom_id.is_none()
            && self.alternate_location.is_none()
    }

    fn matches(&self, candidate: &AtomRow) -> bool {
        self.asym_id
            .as_deref()
            .is_none_or(|expected| candidate.label_asym_id.as_deref() == Some(expected))
            && self
                .component_id
                .as_deref()
                .is_none_or(|expected| candidate.label_comp_id.as_deref() == Some(expected))
            && self
                .sequence_id
                .is_none_or(|expected| candidate.label_seq_id == Some(expected))
            && self
                .atom_id
                .as_deref()
                .is_none_or(|expected| candidate.label_atom_name.as_deref() == Some(expected))
            && self
                .alternate_location
                .as_deref()
                .is_none_or(|expected| candidate.alt_id.as_deref() == Some(expected))
    }
}

#[derive(Debug, Default)]
struct AuthorConnectionSelector {
    asym_id: Option<String>,
    component_id: Option<String>,
    sequence_id: Option<String>,
    atom_id: Option<String>,
    insertion_code: Option<String>,
    alternate_location: Option<String>,
}

impl AuthorConnectionSelector {
    fn is_empty(&self) -> bool {
        self.asym_id.is_none()
            && self.component_id.is_none()
            && self.sequence_id.is_none()
            && self.atom_id.is_none()
            && self.insertion_code.is_none()
            && self.alternate_location.is_none()
    }

    fn matches(&self, candidate: &AtomRow) -> bool {
        self.asym_id
            .as_deref()
            .is_none_or(|expected| candidate.auth_asym_id.as_deref() == Some(expected))
            && self
                .component_id
                .as_deref()
                .is_none_or(|expected| candidate.auth_comp_id.as_deref() == Some(expected))
            && self
                .sequence_id
                .as_deref()
                .is_none_or(|expected| candidate.auth_seq_id.as_deref() == Some(expected))
            && self
                .atom_id
                .as_deref()
                .is_none_or(|expected| candidate.auth_atom_name.as_deref() == Some(expected))
            && self
                .insertion_code
                .as_deref()
                .is_none_or(|expected| candidate.insertion_code.as_deref() == Some(expected))
            && self
                .alternate_location
                .as_deref()
                .is_none_or(|expected| candidate.alt_id.as_deref() == Some(expected))
    }
}

#[derive(Debug, Default)]
struct ConnectionPartnerSelector {
    label: LabelConnectionSelector,
    author: AuthorConnectionSelector,
    conflict: Option<MmcifConnectionResolutionReason>,
}

impl ConnectionPartnerSelector {
    fn is_empty(&self) -> bool {
        self.label.is_empty() && self.author.is_empty()
    }

    fn matches(&self, candidate: &AtomRow) -> bool {
        self.label.matches(candidate) && self.author.matches(candidate)
    }

    fn explicit_alternate_location(&self) -> Option<&str> {
        self.label
            .alternate_location
            .as_deref()
            .or(self.author.alternate_location.as_deref())
    }
}

#[derive(Debug)]
enum ConnectionPartnerResolution<'a> {
    Resolved(&'a AtomRow),
    Unresolved(MmcifConnectionResolutionReason),
    Ambiguous {
        candidates: usize,
        reason: MmcifConnectionResolutionReason,
    },
}

fn connection_partner<'a>(
    table: &MmcifLoopTable,
    row: usize,
    partner: u8,
    selected_rows: &'a [AtomRow],
    all_rows: &[AtomRow],
    selected_model: &str,
) -> Result<ConnectionPartnerResolution<'a>, MmcifInterpretError> {
    let symmetry_tag = format!("_struct_conn.ptnr{partner}_symmetry");
    if let Some(symmetry) = optional(table, row, &symmetry_tag) {
        if symmetry != "1_555" {
            return Ok(ConnectionPartnerResolution::Unresolved(
                MmcifConnectionResolutionReason::UnsupportedSymmetry {
                    symmetry: symmetry.to_owned(),
                },
            ));
        }
    }

    let label_alt_tag = format!("_struct_conn.ptnr{partner}_label_alt_id");
    let pdbx_label_alt_tag = format!("_struct_conn.pdbx_ptnr{partner}_label_alt_id");
    let label_alt = optional(table, row, &label_alt_tag);
    let pdbx_label_alt = optional(table, row, &pdbx_label_alt_tag);
    let (alternate_location, conflict) = match (label_alt, pdbx_label_alt) {
        (Some(left), Some(right)) if left != right => (
            Some(left.to_owned()),
            Some(MmcifConnectionResolutionReason::ConflictingSelectorValues {
                selector: "label alternate-location aliases",
            }),
        ),
        (Some(value), _) | (_, Some(value)) => (Some(value.to_owned()), None),
        (None, None) => (None, None),
    };

    let mut selector = ConnectionPartnerSelector {
        label: LabelConnectionSelector {
            asym_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_label_asym_id"),
            )
            .map(str::to_owned),
            component_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_label_comp_id"),
            )
            .map(str::to_owned),
            sequence_id: optional_i32(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_label_seq_id"),
            )?,
            atom_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_label_atom_id"),
            )
            .map(str::to_owned),
            alternate_location,
        },
        author: AuthorConnectionSelector {
            asym_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_auth_asym_id"),
            )
            .map(str::to_owned),
            component_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_auth_comp_id"),
            )
            .map(str::to_owned),
            sequence_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_auth_seq_id"),
            )
            .map(str::to_owned),
            atom_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_auth_atom_id"),
            )
            .map(str::to_owned),
            insertion_code: optional(
                table,
                row,
                &format!("_struct_conn.pdbx_ptnr{partner}_PDB_ins_code"),
            )
            .map(str::to_owned),
            alternate_location: optional(
                table,
                row,
                &format!("_struct_conn.pdbx_ptnr{partner}_auth_alt_id"),
            )
            .map(str::to_owned),
        },
        conflict,
    };
    if let (Some(label), Some(author)) = (
        selector.label.alternate_location.as_deref(),
        selector.author.alternate_location.as_deref(),
    ) {
        if label != author {
            selector.conflict = Some(MmcifConnectionResolutionReason::ConflictingSelectorValues {
                selector: "label and author alternate locations",
            });
        }
    }
    if let Some(reason) = selector.conflict.take() {
        return Ok(ConnectionPartnerResolution::Unresolved(reason));
    }
    if selector.is_empty() {
        return Ok(ConnectionPartnerResolution::Unresolved(
            MmcifConnectionResolutionReason::MissingSelector,
        ));
    }

    let candidates = selected_rows
        .iter()
        .filter(|candidate| selector.matches(candidate))
        .collect::<Vec<_>>();
    match candidates.as_slice() {
        [candidate] => Ok(ConnectionPartnerResolution::Resolved(candidate)),
        [] => {
            if let Some(alternate_location) = selector.explicit_alternate_location() {
                if all_rows.iter().any(|candidate| {
                    candidate.model_id == selected_model && selector.matches(candidate)
                }) {
                    return Ok(ConnectionPartnerResolution::Unresolved(
                        MmcifConnectionResolutionReason::AlternateLocationOmitted {
                            alternate_location: alternate_location.to_owned(),
                        },
                    ));
                }
            }
            if !selector.label.is_empty()
                && !selector.author.is_empty()
                && selected_rows
                    .iter()
                    .any(|candidate| selector.label.matches(candidate))
                && selected_rows
                    .iter()
                    .any(|candidate| selector.author.matches(candidate))
            {
                return Ok(ConnectionPartnerResolution::Unresolved(
                    MmcifConnectionResolutionReason::ConflictingLabelAndAuthorSelectors,
                ));
            }
            Ok(ConnectionPartnerResolution::Unresolved(
                MmcifConnectionResolutionReason::NoMatchingAtom,
            ))
        }
        candidates => Ok(ConnectionPartnerResolution::Ambiguous {
            candidates: candidates.len(),
            reason: MmcifConnectionResolutionReason::MultipleMatchingAtoms,
        }),
    }
}

#[derive(Debug)]
pub(super) struct InstanceUnion {
    parent: BTreeMap<String, String>,
}

impl InstanceUnion {
    pub(super) fn new(keys: impl IntoIterator<Item = String>) -> Self {
        let parent = keys.into_iter().map(|key| (key.clone(), key)).collect();
        Self { parent }
    }

    pub(super) fn find(&mut self, key: &str) -> String {
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
