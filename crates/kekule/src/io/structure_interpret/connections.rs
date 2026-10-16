//! Connection records are typed once per source and resolve against the atoms
//! of each selected coordinate model.
use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::molecule::{BondOrder, Element};

use super::atom_site::AtomRow;
use super::types::{
    ConnectionResolutionReason, StructureInterpretIssue, StructureInterpretationReport,
};

/// One connection record of a source, before a coordinate model is chosen.
#[derive(Debug, Clone)]
pub(super) enum ConnectionRecord {
    /// A record of a non-covalent type, reported and never bonded.
    Ignored {
        connection_type: String,
    },
    Bond(Box<DeclaredBond>),
}

/// A covalent connection record and the selectors of its two partners.
#[derive(Debug, Clone)]
pub(super) struct DeclaredBond {
    pub(super) connection_id: Option<String>,
    pub(super) connection_type: String,
    pub(super) line: Option<usize>,
    pub(super) order: DeclaredOrder,
    pub(super) partners: [PartnerSpec; 2],
    /// Whether the record lists connectivity only, as PDB `CONECT` does. Such
    /// a bond within one residue supplements the residue's template rather
    /// than replacing it, and a bond an earlier record declares is skipped.
    pub(super) connectivity_only: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum DeclaredOrder {
    /// The order the record states.
    Stated(BondOrder),
    /// A record without a bond type: a bond between a metal and a non-metal
    /// is coordination, and dative; any other bond has this order.
    ByElements(BondOrder),
}

/// The selector of one partner, or why the record cannot identify it.
pub(super) type PartnerSpec = Result<ConnectionPartnerSelector, ConnectionResolutionReason>;

/// A connection resolved to two atoms of the selected model.
#[derive(Debug, Clone)]
pub(super) struct DeclaredConnection {
    pub(super) left_atom: String,
    pub(super) right_atom: String,
    /// The instance of the left atom; both atoms end up in one molecule group.
    pub(super) instance_key: String,
    pub(super) order: BondOrder,
    pub(super) connection_id: Option<String>,
    pub(super) connection_type: String,
    /// A connectivity-only bond within one residue; see
    /// [`DeclaredBond::connectivity_only`].
    pub(super) supplements_template: bool,
}

/// Resolves every record against the selected rows of one coordinate model,
/// uniting the instances that resolved bonds join.
pub(super) fn resolve_connections(
    records: &[ConnectionRecord],
    rows: &[AtomRow],
    all_rows: &[AtomRow],
    selected_model: &str,
    union: &mut InstanceUnion,
    report: &mut StructureInterpretationReport,
) -> Vec<DeclaredConnection> {
    let sites = records
        .iter()
        .any(|record| {
            matches!(
                record,
                ConnectionRecord::Bond(bond)
                    if bond.partners.iter().flatten().any(|selector| selector.site_id.is_some())
            )
        })
        .then(|| SiteIndex::new(rows, all_rows, selected_model));
    let mut connections = Vec::new();
    let mut declared = BTreeSet::new();
    for record in records {
        let bond = match record {
            ConnectionRecord::Ignored { connection_type } => {
                report
                    .issues
                    .push(StructureInterpretIssue::ConnectionIgnored {
                        connection_type: connection_type.clone(),
                    });
                continue;
            }
            ConnectionRecord::Bond(bond) => bond,
        };
        let mut resolve = |partner: u8| {
            let resolution = match &bond.partners[usize::from(partner - 1)] {
                Ok(selector) => {
                    resolve_partner(selector, rows, all_rows, selected_model, sites.as_ref())
                }
                Err(reason) => ConnectionPartnerResolution::Unresolved(reason.clone()),
            };
            report_connection_partner_resolution(resolution, bond, partner, report)
        };
        let left = resolve(1);
        let right = resolve(2);
        let (Some(left), Some(right)) = (left, right) else {
            continue;
        };
        let (left_atom, right_atom) = (left.atom_key(), right.atom_key());
        let pair = if left_atom <= right_atom {
            (left_atom.clone(), right_atom.clone())
        } else {
            (right_atom.clone(), left_atom.clone())
        };
        if !declared.insert(pair) && bond.connectivity_only {
            continue;
        }
        let order = match bond.order {
            DeclaredOrder::Stated(order) => order,
            DeclaredOrder::ByElements(_) if is_metal(left.element) != is_metal(right.element) => {
                BondOrder::Dative
            }
            DeclaredOrder::ByElements(order) => order,
        };
        let supplements_template = bond.connectivity_only
            && left.asym_id == right.asym_id
            && left.residue_key == right.residue_key
            && left.comp_id == right.comp_id;
        if !supplements_template {
            union.union(&left.instance_key, &right.instance_key);
        }
        connections.push(DeclaredConnection {
            left_atom,
            right_atom,
            instance_key: left.instance_key.clone(),
            order,
            connection_id: bond.connection_id.clone(),
            connection_type: bond.connection_type.clone(),
            supplements_template,
        });
        report.applied_connections += 1;
    }
    connections
}

/// Alkali, alkaline-earth, transition, post-transition, lanthanide, and
/// actinide metals; metalloids are not metals.
fn is_metal(element: Element) -> bool {
    matches!(
        element.atomic_number(),
        3 | 4 | 11..=13 | 19..=31 | 37..=50 | 55..=84 | 87..=116
    )
}

fn report_connection_partner_resolution<'a>(
    resolution: ConnectionPartnerResolution<'a>,
    bond: &DeclaredBond,
    partner: u8,
    report: &mut StructureInterpretationReport,
) -> Option<&'a AtomRow> {
    match resolution {
        ConnectionPartnerResolution::Resolved(atom) => Some(atom),
        ConnectionPartnerResolution::Unresolved(reason) => {
            report
                .issues
                .push(StructureInterpretIssue::ConnectionUnresolved {
                    connection_id: bond.connection_id.clone(),
                    connection_type: bond.connection_type.clone(),
                    partner,
                    source_line: bond.line,
                    reason,
                });
            None
        }
        ConnectionPartnerResolution::Ambiguous { candidates, reason } => {
            report
                .issues
                .push(StructureInterpretIssue::ConnectionAmbiguous {
                    connection_id: bond.connection_id.clone(),
                    connection_type: bond.connection_type.clone(),
                    partner,
                    source_line: bond.line,
                    candidates,
                    reason,
                });
            None
        }
    }
}

#[derive(Debug, Clone, Default)]
pub(super) struct LabelConnectionSelector {
    pub(super) asym_id: Option<String>,
    pub(super) component_id: Option<String>,
    pub(super) sequence_id: Option<i32>,
    pub(super) atom_id: Option<String>,
    pub(super) alternate_location: Option<String>,
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

#[derive(Debug, Clone, Default)]
pub(super) struct AuthorConnectionSelector {
    pub(super) asym_id: Option<String>,
    pub(super) component_id: Option<String>,
    pub(super) sequence_id: Option<String>,
    pub(super) atom_id: Option<String>,
    pub(super) insertion_code: Option<String>,
    pub(super) alternate_location: Option<String>,
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

/// Selects one atom by label and author identity, or by its source record
/// identifier (a PDB atom serial number).
#[derive(Debug, Clone, Default)]
pub(super) struct ConnectionPartnerSelector {
    pub(super) label: LabelConnectionSelector,
    pub(super) author: AuthorConnectionSelector,
    pub(super) site_id: Option<String>,
}

impl ConnectionPartnerSelector {
    pub(super) fn is_empty(&self) -> bool {
        self.label.is_empty() && self.author.is_empty() && self.site_id.is_none()
    }

    fn matches(&self, candidate: &AtomRow) -> bool {
        self.label.matches(candidate)
            && self.author.matches(candidate)
            && self
                .site_id
                .as_deref()
                .is_none_or(|expected| candidate.atom_site_id.as_deref() == Some(expected))
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
    Unresolved(ConnectionResolutionReason),
    Ambiguous {
        candidates: usize,
        reason: ConnectionResolutionReason,
    },
}

/// Rows of the selected model by record identifier, for selectors that name
/// one record, such as PDB `CONECT` serial numbers.
struct SiteIndex<'a> {
    selected: HashMap<&'a str, Vec<&'a AtomRow>>,
    model: HashMap<&'a str, Vec<&'a AtomRow>>,
}

impl<'a> SiteIndex<'a> {
    fn new(selected_rows: &'a [AtomRow], all_rows: &'a [AtomRow], selected_model: &str) -> Self {
        let index = |rows: &mut dyn Iterator<Item = &'a AtomRow>| {
            let mut index = HashMap::<&str, Vec<&AtomRow>>::new();
            for row in rows {
                if let Some(site) = row.atom_site_id.as_deref() {
                    index.entry(site).or_default().push(row);
                }
            }
            index
        };
        Self {
            selected: index(&mut selected_rows.iter()),
            model: index(&mut all_rows.iter().filter(|row| row.model_id == selected_model)),
        }
    }
}

fn resolve_partner<'a>(
    selector: &ConnectionPartnerSelector,
    selected_rows: &'a [AtomRow],
    all_rows: &[AtomRow],
    selected_model: &str,
    sites: Option<&SiteIndex<'a>>,
) -> ConnectionPartnerResolution<'a> {
    let candidates = match (selector.site_id.as_deref(), sites) {
        (Some(site), Some(sites)) => sites
            .selected
            .get(site)
            .into_iter()
            .flatten()
            .copied()
            .filter(|candidate| selector.matches(candidate))
            .collect::<Vec<_>>(),
        _ => selected_rows
            .iter()
            .filter(|candidate| selector.matches(candidate))
            .collect(),
    };
    match candidates.as_slice() {
        [candidate] => ConnectionPartnerResolution::Resolved(candidate),
        [] => {
            let omitted = |alternate_location: &str| {
                ConnectionPartnerResolution::Unresolved(
                    ConnectionResolutionReason::AlternateLocationOmitted {
                        alternate_location: alternate_location.to_owned(),
                    },
                )
            };
            if let Some(alternate_location) = selector.explicit_alternate_location() {
                if all_rows.iter().any(|candidate| {
                    candidate.model_id == selected_model && selector.matches(candidate)
                }) {
                    return omitted(alternate_location);
                }
            }
            // A record identifier names one source row, and so its alternate
            // location, without stating it.
            if let (Some(site), Some(sites)) = (selector.site_id.as_deref(), sites) {
                if let Some(alternate_location) = sites
                    .model
                    .get(site)
                    .into_iter()
                    .flatten()
                    .filter(|candidate| selector.matches(candidate))
                    .find_map(|candidate| candidate.alt_id.as_deref())
                {
                    return omitted(alternate_location);
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
                return ConnectionPartnerResolution::Unresolved(
                    ConnectionResolutionReason::ConflictingLabelAndAuthorSelectors,
                );
            }
            ConnectionPartnerResolution::Unresolved(ConnectionResolutionReason::NoMatchingAtom)
        }
        candidates => ConnectionPartnerResolution::Ambiguous {
            candidates: candidates.len(),
            reason: ConnectionResolutionReason::MultipleMatchingAtoms,
        },
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
