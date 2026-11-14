//! Connection records are typed once per source and resolve against the atoms
//! of each selected coordinate model.
use std::cell::OnceCell;
use std::collections::{BTreeSet, HashMap};
use std::ops::Range;

use crate::diagnostics::{Cancellation, FormatError};
use crate::molecule::{BondOrder, Element};

use super::atom_site::AtomRow;
use super::types::{
    poll, ConnectionResolutionReason, StructureInterpretIssue, StructureInterpretationReport,
};

/// One connection record of a source, before a coordinate model is chosen.
#[derive(Debug, Clone)]
pub(super) enum ConnectionRecord {
    /// The records of one non-covalent type, such as the hydrogen bonds a
    /// file lists: never bonded, and reported once with their count.
    Ignored {
        connection_type: String,
        count: usize,
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
    /// The partners' source rows, as [`AtomRow::row_index`].
    pub(super) rows: [usize; 2],
    /// The instance of the left atom; both atoms end up in one molecule group.
    pub(super) instance_key: String,
    pub(super) order: BondOrder,
    pub(super) connection_id: Option<String>,
    pub(super) connection_type: String,
    /// A connectivity-only bond within one residue; see
    /// [`DeclaredBond::connectivity_only`].
    pub(super) supplements_template: bool,
}

/// Metal coordination a record states between atoms of different residues:
/// source connection metadata, never a bond, so it does not join molecules.
#[derive(Debug, Clone)]
pub(super) struct DeclaredCoordination {
    /// The atom keys of the partners, in record order.
    pub(super) atoms: [String; 2],
    pub(super) connection_id: Option<String>,
    pub(super) connection_type: String,
    pub(super) line: Option<usize>,
}

/// The resolved records of one coordinate model.
#[derive(Debug, Default)]
pub(super) struct ResolvedConnections {
    pub(super) bonds: Vec<DeclaredConnection>,
    pub(super) coordination: Vec<DeclaredCoordination>,
}

/// Resolves every record against the selected rows of one coordinate model,
/// uniting the instances that resolved bonds join.
pub(super) fn resolve_connections(
    records: &[ConnectionRecord],
    rows: &[&AtomRow],
    all_rows: &[AtomRow],
    selected_model: &str,
    union: &mut InstanceUnion,
    report: &mut StructureInterpretationReport,
    cancellation: Option<&Cancellation>,
) -> Result<ResolvedConnections, FormatError> {
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
    let names = NameIndex::new(rows);
    let mut resolved = ResolvedConnections::default();
    let mut declared = BTreeSet::new();
    for record in records {
        // Resolving a partner can scan every row of the model.
        poll(cancellation)?;
        let bond = match record {
            ConnectionRecord::Ignored {
                connection_type,
                count,
            } => {
                report
                    .issues
                    .push(StructureInterpretIssue::ConnectionsIgnored {
                        connection_type: connection_type.clone(),
                        count: *count,
                    });
                continue;
            }
            ConnectionRecord::Bond(bond) => bond,
        };
        let mut resolve = |partner: u8| {
            let resolution = match &bond.partners[usize::from(partner - 1)] {
                Ok(selector) => resolve_partner(
                    selector,
                    rows,
                    all_rows,
                    selected_model,
                    sites.as_ref(),
                    &names,
                ),
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
        let same_residue = left.asym_id == right.asym_id
            && left.residue_key == right.residue_key
            && left.comp_id == right.comp_id;
        report.applied_connections += 1;
        // Coordination between residues is recorded, not bonded: covalent
        // connectedness alone decides molecule boundaries.
        if order == BondOrder::Dative && !same_residue {
            resolved.coordination.push(DeclaredCoordination {
                atoms: [left_atom, right_atom],
                connection_id: bond.connection_id.clone(),
                connection_type: bond.connection_type.clone(),
                line: bond.line,
            });
            continue;
        }
        let supplements_template = bond.connectivity_only && same_residue;
        if !supplements_template {
            union.union(&left.instance_key, &right.instance_key);
        }
        resolved.bonds.push(DeclaredConnection {
            rows: [left.row_index, right.row_index],
            instance_key: left.instance_key.clone(),
            order,
            connection_id: bond.connection_id.clone(),
            connection_type: bond.connection_type.clone(),
            supplements_template,
        });
    }
    Ok(resolved)
}

/// Alkali, alkaline-earth, transition, post-transition, lanthanide, and
/// actinide metals; metalloids are not metals.
pub(crate) fn is_metal(element: Element) -> bool {
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
    fn new(selected_rows: &[&'a AtomRow], all_rows: &'a [AtomRow], selected_model: &str) -> Self {
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
            selected: index(&mut selected_rows.iter().copied()),
            model: index(&mut all_rows.iter().filter(|row| row.model_id == selected_model)),
        }
    }
}

/// Rows by asym and atom name.
type Named<'a> = HashMap<(&'a str, &'a str), Vec<&'a AtomRow>>;

/// Runs of consecutive rows, as ranges of a row list, by key.
type Runs<K> = HashMap<K, Vec<Range<usize>>>;

/// Rows of the selected model by residue and by asym and atom name, under
/// label and under author identity, so that resolving a selector reads only
/// the rows that can match it. Each lookup is built on first use.
struct NameIndex<'a> {
    rows: Vec<&'a AtomRow>,
    label: OnceCell<Named<'a>>,
    author: OnceCell<Named<'a>>,
    /// Runs of rows, as ranges of `rows`, by label asym and sequence number.
    label_residues: OnceCell<Runs<(&'a str, i32)>>,
    /// Runs of rows by author asym and sequence number.
    author_residues: OnceCell<Runs<(&'a str, &'a str)>>,
}

impl<'a> NameIndex<'a> {
    fn new(rows: &[&'a AtomRow]) -> Self {
        Self {
            rows: rows.to_vec(),
            label: OnceCell::new(),
            author: OnceCell::new(),
            label_residues: OnceCell::new(),
            author_residues: OnceCell::new(),
        }
    }

    /// The rows, in source order, among which every row matching `selector`
    /// is, or `None` for a selector stating no asym with a residue or an
    /// atom name. A residue, stated under either identity, narrows the most.
    fn candidates(&self, selector: &ConnectionPartnerSelector) -> Option<Vec<&'a AtomRow>> {
        let (label, author) = (&selector.label, &selector.author);
        if label.sequence_id.is_none() && author.sequence_id.is_some() {
            self.author_rows(author).or_else(|| self.label_rows(label))
        } else {
            self.label_rows(label).or_else(|| self.author_rows(author))
        }
    }

    /// The rows, in source order, among which every row matching `label` is:
    /// its residue's rows, or its atom name's rows within its asym. `None`
    /// for a selector stating neither.
    fn label_rows(&self, label: &LabelConnectionSelector) -> Option<Vec<&'a AtomRow>> {
        let asym = label.asym_id.as_deref()?;
        match (label.sequence_id, label.atom_id.as_deref()) {
            (Some(sequence), _) => {
                let residues = self.label_residues.get_or_init(|| {
                    runs(&self.rows, |row| {
                        row.label_asym_id.as_deref().zip(row.label_seq_id)
                    })
                });
                Some(self.residue_rows(residues.get(&(asym, sequence))))
            }
            (None, Some(atom)) => {
                let names = self.label.get_or_init(|| {
                    names(&self.rows, |row| {
                        row.label_asym_id
                            .as_deref()
                            .zip(row.label_atom_name.as_deref())
                    })
                });
                Some(names.get(&(asym, atom)).cloned().unwrap_or_default())
            }
            (None, None) => None,
        }
    }

    /// [`Self::label_rows`] for an author selector.
    fn author_rows(&self, author: &AuthorConnectionSelector) -> Option<Vec<&'a AtomRow>> {
        let asym = author.asym_id.as_deref()?;
        match (author.sequence_id.as_deref(), author.atom_id.as_deref()) {
            (Some(sequence), _) => {
                let residues = self.author_residues.get_or_init(|| {
                    runs(&self.rows, |row| {
                        row.auth_asym_id.as_deref().zip(row.auth_seq_id.as_deref())
                    })
                });
                Some(self.residue_rows(residues.get(&(asym, sequence))))
            }
            (None, Some(atom)) => {
                let names = self.author.get_or_init(|| {
                    names(&self.rows, |row| {
                        row.auth_asym_id
                            .as_deref()
                            .zip(row.auth_atom_name.as_deref())
                    })
                });
                Some(names.get(&(asym, atom)).cloned().unwrap_or_default())
            }
            (None, None) => None,
        }
    }

    /// The rows of a residue's runs, in source order.
    fn residue_rows(&self, runs: Option<&Vec<Range<usize>>>) -> Vec<&'a AtomRow> {
        runs.into_iter()
            .flatten()
            .flat_map(|run| &self.rows[run.clone()])
            .copied()
            .collect()
    }
}

/// The rows with each key, in source order.
fn names<'a>(
    rows: &[&'a AtomRow],
    key: impl Fn(&'a AtomRow) -> Option<(&'a str, &'a str)>,
) -> Named<'a> {
    let mut names = HashMap::<_, Vec<_>>::new();
    for &row in rows {
        if let Some(key) = key(row) {
            names.entry(key).or_default().push(row);
        }
    }
    names
}

/// The runs of consecutive rows with one key, as ranges of `rows`, by key.
fn runs<'a, K: std::hash::Hash + Eq + Copy>(
    rows: &[&'a AtomRow],
    key: impl Fn(&'a AtomRow) -> Option<K>,
) -> Runs<K> {
    let mut runs = Runs::<K>::new();
    let mut current: Option<(K, usize)> = None;
    for (position, &row) in rows.iter().enumerate() {
        let row_key = key(row);
        if let Some((key, start)) = current {
            if row_key == Some(key) {
                continue;
            }
            runs.entry(key).or_default().push(start..position);
        }
        current = row_key.map(|key| (key, position));
    }
    if let Some((key, start)) = current {
        runs.entry(key).or_default().push(start..rows.len());
    }
    runs
}

fn resolve_partner<'a>(
    selector: &ConnectionPartnerSelector,
    selected_rows: &[&'a AtomRow],
    all_rows: &[AtomRow],
    selected_model: &str,
    sites: Option<&SiteIndex<'a>>,
    names: &NameIndex<'a>,
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
        _ => match names.candidates(selector) {
            Some(rows) => rows
                .into_iter()
                .filter(|candidate| selector.matches(candidate))
                .collect(),
            None => selected_rows
                .iter()
                .copied()
                .filter(|candidate| selector.matches(candidate))
                .collect(),
        },
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
            // Whether any selected row matches, reading only the rows an
            // index narrows the search to where it can.
            let any = |narrowed: Option<Vec<&AtomRow>>, matches: &dyn Fn(&AtomRow) -> bool| {
                match narrowed {
                    Some(rows) => rows.into_iter().any(matches),
                    None => selected_rows.iter().any(|candidate| matches(candidate)),
                }
            };
            if !selector.label.is_empty()
                && !selector.author.is_empty()
                && any(names.label_rows(&selector.label), &|candidate| {
                    selector.label.matches(candidate)
                })
                && any(names.author_rows(&selector.author), &|candidate| {
                    selector.author.matches(candidate)
                })
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
/// Instance keys united into molecules; each set's root is its smallest
/// key.
pub(super) struct InstanceUnion {
    indices: HashMap<String, usize>,
    keys: Vec<String>,
    parents: Vec<usize>,
}

impl InstanceUnion {
    pub(super) fn new<'a>(keys: impl IntoIterator<Item = &'a str>) -> Self {
        let mut union = Self {
            indices: HashMap::new(),
            keys: Vec::new(),
            parents: Vec::new(),
        };
        for key in keys {
            union.index(key);
        }
        union
    }

    fn index(&mut self, key: &str) -> usize {
        if let Some(&index) = self.indices.get(key) {
            return index;
        }
        let index = self.keys.len();
        self.indices.insert(key.to_owned(), index);
        self.keys.push(key.to_owned());
        self.parents.push(index);
        index
    }

    fn root(&mut self, key: &str) -> usize {
        let mut current = self.index(key);
        let start = current;
        while self.parents[current] != current {
            current = self.parents[current];
        }
        let root = current;
        let mut current = start;
        while self.parents[current] != root {
            let next = self.parents[current];
            self.parents[current] = root;
            current = next;
        }
        root
    }

    /// The root key of `key`'s set.
    pub(super) fn find(&mut self, key: &str) -> &str {
        let root = self.root(key);
        &self.keys[root]
    }

    fn union(&mut self, left: &str, right: &str) {
        let left_root = self.root(left);
        let right_root = self.root(right);
        if left_root != right_root {
            let (root, child) = if self.keys[left_root] < self.keys[right_root] {
                (left_root, right_root)
            } else {
                (right_root, left_root)
            };
            self.parents[child] = root;
        }
    }
}
