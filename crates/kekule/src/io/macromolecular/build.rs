//! Builds molecules, provenance, and the topology hierarchy from the selected
//! atom rows of one coordinate model.
use super::types::{poll, poll_every, structure_error};
use crate::diagnostics::{Cancellation, DiagnosticKind, FormatError};
use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::io::canonicalize_molecule_for_publication;
use crate::math::Point3;
use crate::molecule::{Atom, AtomId, Element, Molecule, MoleculeEditor};
use crate::structure::Positions;
use crate::topology::{AtomSiteMetadata, Hierarchy};
use crate::topology::{InstanceAtomId, MoleculeInstanceId};
use crate::units::{Quantity, ANGSTROM};

use super::super::element_symbol::structure_element_symbol;
use super::super::staged_coordinates::StagedCoordinates;
use super::atom_site::AtomRow;
use super::connections::{DeclaredConnection, InstanceUnion};
use super::connectivity::{record_evidence, BondEvidenceMap, ListedBond};
use super::types::{
    AtomProvenance, BondEvidence, BondProvenance, EntityKind, InstanceProvenance,
    StructureInterpretIssue, StructureInterpretationReport,
};

#[derive(Debug)]
pub(super) struct MoleculeGroup {
    /// The instance union root shared by every row.
    root: String,
    /// The group's rows, as indices into the selected rows.
    rows: Vec<usize>,
    kinds: BTreeSet<EntityKind>,
}

impl MoleculeGroup {
    pub(super) fn root(&self) -> &str {
        &self.root
    }
}

pub(super) fn group_rows(
    rows: &[AtomRow],
    union: &mut InstanceUnion,
    polymer_asym_order: &BTreeMap<String, usize>,
    cancellation: Option<&Cancellation>,
) -> Result<Vec<MoleculeGroup>, FormatError> {
    let mut group_indices = BTreeMap::<String, usize>::new();
    let mut groups = Vec::new();
    // Rows of one instance usually follow each other, so its group is looked
    // up once per run of them.
    let mut last: Option<(String, usize)> = None;
    for (step, row) in rows.iter().enumerate() {
        poll_every(cancellation, step)?;
        let index = match &last {
            Some((key, index)) if *key == row.instance_key => *index,
            _ => {
                let root = union.find(&row.instance_key);
                let index = match group_indices.get(root) {
                    Some(&index) => index,
                    None => {
                        group_indices.insert(root.to_owned(), groups.len());
                        groups.push(MoleculeGroup {
                            root: root.to_owned(),
                            rows: Vec::new(),
                            kinds: BTreeSet::new(),
                        });
                        groups.len() - 1
                    }
                };
                last = Some((row.instance_key.clone(), index));
                index
            }
        };
        let group = &mut groups[index];
        if !group.kinds.contains(&row.kind) {
            group.kinds.insert(row.kind.clone());
        }
        group.rows.push(step);
    }
    for group in &mut groups {
        poll(cancellation)?;
        group.rows.sort_by_cached_key(|&row| {
            polymer_asym_order
                .get(&rows[row].asym_id)
                .copied()
                .unwrap_or(usize::MAX)
        });
    }
    groups.sort_by_cached_key(|group| {
        group
            .rows
            .iter()
            .filter_map(|&row| polymer_asym_order.get(&rows[row].asym_id).copied())
            .min()
            .unwrap_or(usize::MAX)
    });
    Ok(groups)
}

pub(super) struct BuiltMolecule {
    pub(super) editor: MoleculeEditor,
    pub(super) coordinates: StagedCoordinates,
    pub(super) provenance: BuiltMoleculeProvenance,
    evidence: BondEvidenceMap,
    /// Bonds connectivity-only records list within one residue, for the
    /// connectivity resolver.
    listed: Vec<ListedBond>,
}

pub(super) struct PublishedMolecule {
    pub(super) molecule: Molecule,
    pub(super) positions: Positions,
    pub(super) provenance: BuiltMoleculeProvenance,
}

#[derive(Clone)]
pub(super) struct BuiltMoleculeProvenance {
    coordinate_model_id: String,
    asym_ids: Vec<String>,
    entity_ids: Vec<String>,
    entity_kinds: Vec<EntityKind>,
    atoms: Vec<BuiltAtomProvenance>,
    bonds: Vec<([AtomId; 2], BondEvidence)>,
}

#[derive(Clone)]
struct BuiltAtomProvenance {
    atom: AtomId,
    element: Element,
    point: Point3,
    stated_charge: Option<i8>,
    type_symbol: String,
    source_line: usize,
    source_order: usize,
    atom_site_id: Option<String>,
    label_atom_name: Option<String>,
    atom_name: String,
    auth_atom_name: Option<String>,
    label_component_id: Option<String>,
    component_id: String,
    auth_component_id: Option<String>,
    label_asym_id: Option<String>,
    asym_id: String,
    auth_asym_id: Option<String>,
    entity_id: Option<String>,
    entity_kind: EntityKind,
    residue_key: String,
    label_sequence_id: Option<i32>,
    polymer_position: Option<i32>,
    author_sequence_id: Option<String>,
    insertion_code: Option<String>,
    occurrence: Option<usize>,
    selected_alternate_location: Option<String>,
    occupancy: Option<f64>,
    b_factor: Option<f64>,
}

pub(super) type QualifiedAtomProperties = Vec<(InstanceAtomId, Option<f64>, Option<f64>)>;

impl BuiltMoleculeProvenance {
    /// Each atom's component and atom name, by local atom ID: what a shared
    /// definition's instances must agree on besides chemistry.
    pub(super) fn atom_labels(&self, atom_count: usize) -> Vec<(String, String)> {
        let mut labels = vec![(String::new(), String::new()); atom_count];
        for atom in &self.atoms {
            labels[atom.atom.index()] = (atom.component_id.clone(), atom.atom_name.clone());
        }
        labels
    }

    pub(super) fn qualify(
        self,
        molecule: MoleculeInstanceId,
    ) -> (InstanceProvenance, QualifiedAtomProperties) {
        let mut atom_data = Vec::with_capacity(self.atoms.len());
        let provenance = InstanceProvenance {
            molecule,
            coordinate_model_id: self.coordinate_model_id,
            asym_ids: self.asym_ids,
            entity_ids: self.entity_ids,
            entity_kinds: self.entity_kinds,
            bonds: self
                .bonds
                .into_iter()
                .map(|([a, b], evidence)| BondProvenance {
                    atoms: [
                        InstanceAtomId::new(molecule, a),
                        InstanceAtomId::new(molecule, b),
                    ],
                    evidence,
                })
                .collect(),
            atoms: self
                .atoms
                .into_iter()
                .map(|atom| {
                    let qualified = InstanceAtomId::new(molecule, atom.atom);
                    atom_data.push((qualified, atom.occupancy, atom.b_factor));
                    AtomProvenance {
                        atom: qualified,
                        type_symbol: atom.type_symbol,
                        source_line: atom.source_line,
                        source_order: atom.source_order,
                        atom_site_id: atom.atom_site_id,
                        label_atom_name: atom.label_atom_name,
                        atom_name: atom.atom_name,
                        auth_atom_name: atom.auth_atom_name,
                        label_component_id: atom.label_component_id,
                        component_id: atom.component_id,
                        auth_component_id: atom.auth_component_id,
                        label_asym_id: atom.label_asym_id,
                        asym_id: atom.asym_id,
                        auth_asym_id: atom.auth_asym_id,
                        entity_id: atom.entity_id,
                        entity_kind: atom.entity_kind,
                        residue_key: atom.residue_key,
                        label_sequence_id: atom.label_sequence_id,
                        author_sequence_id: atom.author_sequence_id,
                        insertion_code: atom.insertion_code,
                        occurrence: atom.occurrence,
                        selected_alternate_location: atom.selected_alternate_location,
                    }
                })
                .collect(),
        };
        (provenance, atom_data)
    }
}

impl BuiltMolecule {
    pub(super) fn complete_connectivity(
        &mut self,
        catalog: &super::connectivity::ConnectivityCatalog,
        report: &mut StructureInterpretationReport,
        cancellation: Option<&Cancellation>,
    ) -> Result<(), FormatError> {
        let Self {
            editor,
            provenance,
            evidence,
            listed,
            ..
        } = self;
        let atoms = provenance
            .atoms
            .iter()
            .map(|atom| super::connectivity::StagedAtomProvenance {
                atom: atom.atom,
                atom_name: &atom.atom_name,
                component_id: &atom.component_id,
                asym_id: &atom.asym_id,
                entity_id: atom.entity_id.as_deref(),
                label_sequence_id: atom.label_sequence_id,
                polymer_position: atom.polymer_position,
                author_sequence_id: atom.author_sequence_id.as_deref(),
                insertion_code: atom.insertion_code.as_deref(),
                occurrence: atom.occurrence,
                element: atom.element,
                point: atom.point,
                line: atom.source_line,
                stated_charge: atom.stated_charge,
            });
        super::connectivity::complete_editor_connectivity(
            catalog,
            editor,
            atoms,
            listed,
            &provenance.coordinate_model_id,
            evidence,
            &mut report.issues,
            cancellation,
        )?;
        poll(cancellation)?;
        // Close contacts that no evidence bonds, counted once every bond is known.
        let sites = provenance
            .atoms
            .iter()
            .map(|atom| (atom.atom, atom.element, atom.point, atom.source_line))
            .collect::<Vec<_>>();
        let candidates = count_bond_candidates(editor.working(), &sites, cancellation)?;
        if candidates > 0 {
            report.connectivity_candidates += candidates;
            report
                .issues
                .push(StructureInterpretIssue::ConnectivityCandidatesInferred {
                    atom_count: editor.working().atom_count(),
                    candidate_count: candidates,
                });
        }
        Ok(())
    }

    pub(super) fn publish_components(
        self,
        cancellation: Option<&Cancellation>,
    ) -> Result<Vec<PublishedMolecule>, FormatError> {
        let Self {
            editor: source,
            coordinates: source_coordinates,
            provenance: source_provenance,
            evidence,
            ..
        } = self;
        let working = source.working();
        let components = working.connected_components();
        if components.is_empty() {
            return Err(structure_error(
                DiagnosticKind::Chemistry,
                None,
                "A molecule group has no atoms.",
            ));
        }
        // Each source atom's component and its ID there.
        let slots = working.graph.atom_slot_count();
        let mut component_of = vec![usize::MAX; slots];
        let mut local = vec![AtomId::new(0); slots];
        let mut editors = Vec::with_capacity(components.len());
        for (index, component) in components.iter().enumerate() {
            poll(cancellation)?;
            let mut editor = crate::molecule::MoleculeEditor::new();
            for (step, source_atom) in component.iter().copied().enumerate() {
                poll_every(cancellation, step)?;
                let atom = working
                    .require_atom(source_atom)
                    .map_err(graph_error)?
                    .clone();
                component_of[source_atom.index()] = index;
                local[source_atom.index()] = editor.add_atom(atom).map_err(graph_error)?;
            }
            editors.push(editor);
        }
        // Every bond lies within the component of its endpoints; each list
        // keeps the source bond order.
        let mut component_bonds = vec![Vec::new(); components.len()];
        for (step, (id, bond)) in working.bonds().enumerate() {
            poll_every(cancellation, step)?;
            component_bonds[component_of[bond.a().index()]].push((id, bond));
        }
        let BuiltMoleculeProvenance {
            coordinate_model_id,
            asym_ids,
            entity_ids,
            entity_kinds,
            atoms: provenance_atoms,
            ..
        } = source_provenance;
        let mut component_atoms = vec![Vec::new(); components.len()];
        for (step, mut atom) in provenance_atoms.into_iter().enumerate() {
            poll_every(cancellation, step)?;
            let Some(&component) = component_of.get(atom.atom.index()) else {
                continue;
            };
            if component == usize::MAX {
                continue;
            }
            atom.atom = local[atom.atom.index()];
            component_atoms[component].push(atom);
        }
        let mut published = Vec::with_capacity(components.len());
        for (((component, mut editor), bonds_of), atoms_of) in components
            .iter()
            .zip(editors)
            .zip(component_bonds)
            .zip(component_atoms)
        {
            poll(cancellation)?;
            let mut bonds = Vec::with_capacity(bonds_of.len());
            for (step, (_, source_bond)) in bonds_of.iter().enumerate() {
                poll_every(cancellation, step)?;
                let (left, right) = source_bond.endpoints();
                let atoms = [local[left.index()], local[right.index()]];
                editor
                    .add_bond(atoms[0], atoms[1], source_bond.order)
                    .map_err(graph_error)?;
                let evidence = evidence
                    .get(&(left.min(right), left.max(right)))
                    .cloned()
                    .ok_or_else(|| {
                        structure_error(
                            DiagnosticKind::Chemistry,
                            None,
                            "A bond has no source record.",
                        )
                    })?;
                bonds.push((atoms, evidence));
            }
            let bond_sources = bonds_of
                .iter()
                .map(|(id, _)| id.index())
                .collect::<Vec<_>>();
            let mut properties =
                crate::properties::MoleculeProperties::new(component.len(), bond_sources.len());
            *properties.atoms_mut() = working
                .properties()
                .atoms()
                .raw()
                .select_indices(&component.iter().map(|id| id.index()).collect::<Vec<_>>())
                .map_err(graph_error)?;
            *properties.bonds_mut() = working
                .properties()
                .bonds()
                .raw()
                .select_indices(&bond_sources)
                .map_err(graph_error)?;
            editor.working_mut().properties = properties;
            let mut coordinates =
                StagedCoordinates::with_atom_capacity(component.len(), source_coordinates.unit())
                    .map_err(graph_error)?;
            for source_atom in component {
                if let Some(point) = source_coordinates.position(*source_atom) {
                    coordinates
                        .set_position(local[source_atom.index()], point)
                        .map_err(graph_error)?;
                }
            }
            poll(cancellation)?;
            canonicalize_molecule_for_publication(editor.working_mut(), Some(&coordinates), &[])
                .map_err(graph_error)?;
            poll(cancellation)?;
            let molecule = editor.finish().map_err(graph_error)?;
            poll(cancellation)?;
            let positions = coordinates.to_positions(&molecule).map_err(graph_error)?;
            poll(cancellation)?;
            published.push(PublishedMolecule {
                molecule,
                positions,
                provenance: BuiltMoleculeProvenance {
                    coordinate_model_id: coordinate_model_id.clone(),
                    asym_ids: asym_ids.clone(),
                    entity_ids: entity_ids.clone(),
                    entity_kinds: entity_kinds.clone(),
                    atoms: atoms_of,
                    bonds,
                },
            });
        }
        Ok(published)
    }
}

/// Builds the molecule of `group`, whose rows' text moves out of `selected`.
pub(super) fn build_molecule(
    group: MoleculeGroup,
    selected: &mut [AtomRow],
    connections: &[&DeclaredConnection],
    cancellation: Option<&Cancellation>,
) -> Result<BuiltMolecule, FormatError> {
    let MoleculeGroup {
        rows: indices,
        kinds,
        ..
    } = group;
    let rows = RowsOf {
        rows: selected,
        indices: &indices,
    };
    // One atom per atom key (its asym, residue, and name), represented by its
    // first row; `row_atoms` maps every row to its atom, and `endpoints` the
    // connections' source rows.
    let mut endpoints = connections
        .iter()
        .flat_map(|connection| connection.rows)
        .map(|row| (row, None))
        .collect::<HashMap<usize, Option<usize>>>();
    let mut keys = HashMap::<(&str, &str, &str), usize>::with_capacity(rows.len());
    let mut representatives = Vec::<usize>::with_capacity(rows.len());
    let mut row_atoms = Vec::with_capacity(rows.len());
    for (step, row) in rows.iter().enumerate() {
        poll_every(cancellation, step)?;
        let key = (
            row.asym_id.as_str(),
            row.residue_key.as_str(),
            row.atom_name.as_str(),
        );
        if let Some(&index) = keys.get(&key) {
            let prior = &rows[representatives[index]];
            if prior.element != row.element
                || prior.isotope != row.isotope
                || prior.formal_charge != row.formal_charge
                || prior.comp_id != row.comp_id
                || prior.entity_id != row.entity_id
            {
                return Err(structure_error(
                    DiagnosticKind::Inconsistent,
                    Some(row.line),
                    format!(
                        "Atom {:?} differs in its chemistry or identity between coordinate models.",
                        row.atom_name
                    ),
                ));
            }
            row_atoms.push(index);
        } else {
            keys.insert(key, representatives.len());
            row_atoms.push(representatives.len());
            representatives.push(step);
        }
        if let Some(atom) = endpoints.get_mut(&row.row_index) {
            *atom = row_atoms.last().copied();
        }
    }
    drop(keys);
    let mut editor = crate::molecule::MoleculeEditor::new();
    let mut atoms = Vec::with_capacity(representatives.len());
    for (step, &index) in representatives.iter().enumerate() {
        poll_every(cancellation, step)?;
        let row = &rows[index];
        let mut atom = Atom::new(row.element);
        atom.isotope = row.isotope.and_then(std::num::NonZeroU16::new);
        atom.formal_charge = row.formal_charge.unwrap_or_default();
        atoms.push(editor.add_atom(atom).map_err(graph_error)?);
    }
    let model_id = rows
        .first()
        .map(|row| row.model_id.clone())
        .ok_or_else(|| {
            structure_error(
                DiagnosticKind::Chemistry,
                None,
                "A molecule group has no atoms.",
            )
        })?;
    let mut coordinates =
        StagedCoordinates::with_atom_capacity(atoms.len(), ANGSTROM).map_err(graph_error)?;
    for (step, (row, &atom)) in rows.iter().zip(&row_atoms).enumerate() {
        poll_every(cancellation, step)?;
        let point = row.point.ok_or_else(|| {
            structure_error(
                DiagnosticKind::InvalidValue,
                Some(row.line),
                format!(
                    "Atom {:?} has no position in the selected model.",
                    row.atom_name
                ),
            )
        })?;
        coordinates
            .set_position(atoms[atom], Quantity::new(point, ANGSTROM))
            .map_err(graph_error)?;
    }
    let mut evidence = BondEvidenceMap::new();
    let mut listed = Vec::new();
    for (step, connection) in connections.iter().enumerate() {
        poll_every(cancellation, step)?;
        let [Some(left), Some(right)] = connection.rows.map(|row| endpoints[&row]) else {
            continue;
        };
        let (left, right) = (atoms[left], atoms[right]);
        let source = BondEvidence::Connection {
            connection_id: connection.connection_id.clone(),
            connection_type: connection.connection_type.clone(),
        };
        if connection.supplements_template {
            listed.push(ListedBond {
                atoms: [left, right],
                order: connection.order,
                evidence: source,
            });
            continue;
        }
        record_evidence(&mut evidence, left, right, source);
        if editor
            .working()
            .require_bond_between(left, right)
            .map_err(graph_error)?
            .is_none()
        {
            editor
                .add_bond(left, right, connection.order)
                .map_err(graph_error)?;
        } else {
            let existing = editor
                .working()
                .require_bond_between(left, right)
                .map_err(graph_error)?
                .expect("existing bond was found");
            if editor
                .working()
                .require_bond(existing)
                .map_err(graph_error)?
                .order
                != connection.order
            {
                return Err(structure_error(
                    DiagnosticKind::Inconsistent,
                    None,
                    "Two connection records give one bond different orders.",
                ));
            }
        }
    }
    poll(cancellation)?;
    let asym_ids = representatives
        .iter()
        .map(|&index| rows[index].asym_id.as_str())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let entity_ids = representatives
        .iter()
        .filter_map(|&index| rows[index].entity_id.as_deref())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let entity_kinds = kinds.into_iter().collect::<Vec<_>>();
    poll(cancellation)?;
    // Representatives are distinct rows, so their text moves into provenance.
    let atom_provenance = representatives
        .iter()
        .zip(&atoms)
        .enumerate()
        .map(|(step, (&index, &atom))| {
            poll_every(cancellation, step)?;
            let row = &mut rows.rows[rows.indices[index]];
            Ok(BuiltAtomProvenance {
                atom,
                element: row.element,
                point: row
                    .point
                    .expect("selected mmCIF atom rows have complete positions"),
                stated_charge: row.formal_charge,
                type_symbol: structure_element_symbol(row.element, row.isotope)
                    .expect("read element symbols name their isotope")
                    .to_owned(),
                source_line: row.line,
                source_order: row.site_order,
                atom_site_id: row.atom_site_id.take(),
                label_atom_name: row.label_atom_name.take(),
                atom_name: std::mem::take(&mut row.atom_name),
                auth_atom_name: row.auth_atom_name.take(),
                label_component_id: row.label_comp_id.take(),
                component_id: std::mem::take(&mut row.comp_id),
                auth_component_id: row.auth_comp_id.take(),
                label_asym_id: row.label_asym_id.take(),
                asym_id: std::mem::take(&mut row.asym_id),
                auth_asym_id: row.auth_asym_id.take(),
                entity_id: row.entity_id.take(),
                entity_kind: row.kind.clone(),
                residue_key: std::mem::take(&mut row.residue_key),
                label_sequence_id: row.label_seq_id,
                polymer_position: row.polymer_position,
                author_sequence_id: row.auth_seq_id.take(),
                insertion_code: row.insertion_code.take(),
                occurrence: row.occurrence,
                selected_alternate_location: row.alt_id.take(),
                occupancy: row.occupancy,
                b_factor: row.b_factor,
            })
        })
        .collect::<Result<_, FormatError>>()?;
    let provenance = BuiltMoleculeProvenance {
        coordinate_model_id: model_id,
        asym_ids,
        entity_ids,
        entity_kinds,
        atoms: atom_provenance,
        bonds: Vec::new(),
    };
    Ok(BuiltMolecule {
        editor,
        coordinates,
        provenance,
        evidence,
        listed,
    })
}

/// A group's rows within the selected rows, indexable in group order.
struct RowsOf<'a> {
    rows: &'a mut [AtomRow],
    indices: &'a [usize],
}

impl RowsOf<'_> {
    fn first(&self) -> Option<&AtomRow> {
        self.indices.first().map(|&index| &self.rows[index])
    }

    fn len(&self) -> usize {
        self.indices.len()
    }

    fn iter(&self) -> impl Iterator<Item = &AtomRow> {
        self.indices.iter().map(|&index| &self.rows[index])
    }
}

impl std::ops::Index<usize> for RowsOf<'_> {
    type Output = AtomRow;

    fn index(&self, index: usize) -> &AtomRow {
        &self.rows[self.indices[index]]
    }
}

const COVALENT_BOND_CELL_ANGSTROM: f64 = 2.1;
const COVALENT_BOND_TOLERANCE_ANGSTROM: f64 = 0.45;
const MIN_COVALENT_BOND_DISTANCE_SQUARED: f64 = 0.16;
const FALLBACK_COVALENT_RADIUS_ANGSTROM: f64 = 0.77;

/// Counts atom pairs within covalent bonding distance that are not bonded.
fn count_bond_candidates(
    graph: &Molecule,
    sites: &[(AtomId, Element, Point3, usize)],
    cancellation: Option<&Cancellation>,
) -> Result<usize, FormatError> {
    // No cutoff exceeds the cell edge, so a pair within bonding distance lies
    // in one cell or in two adjacent ones.
    let cells = sites
        .iter()
        .enumerate()
        .map(|(step, &(_, _, point, line))| {
            poll_every(cancellation, step)?;
            covalent_bond_cell(point, line)
        })
        .collect::<Result<Vec<_>, FormatError>>()?;
    let mut order = (0..sites.len()).collect::<Vec<_>>();
    order.sort_unstable_by_key(|&index| (cells[index], index));
    // Occupied cells in order, each with its range of `order`.
    let mut occupied = Vec::<([i64; 3], usize, usize)>::new();
    for (position, &index) in order.iter().enumerate() {
        match occupied.last_mut() {
            Some((cell, _, end)) if *cell == cells[index] => *end = position + 1,
            _ => occupied.push((cells[index], position, position + 1)),
        }
    }
    let mut candidates = 0usize;
    let mut count = |left: usize, right: usize| -> Result<(), FormatError> {
        let (left, left_element, left_point, _) = sites[left];
        let (right, right_element, right_point, _) = sites[right];
        let distance_squared = point_distance_squared(left_point, right_point);
        if distance_squared <= MIN_COVALENT_BOND_DISTANCE_SQUARED {
            return Ok(());
        }
        let left_radius = left_element
            .covalent_radius_angstrom()
            .unwrap_or(FALLBACK_COVALENT_RADIUS_ANGSTROM);
        let right_radius = right_element
            .covalent_radius_angstrom()
            .unwrap_or(FALLBACK_COVALENT_RADIUS_ANGSTROM);
        let cutoff = (left_radius + right_radius + COVALENT_BOND_TOLERANCE_ANGSTROM)
            .min(COVALENT_BOND_CELL_ANGSTROM);
        if distance_squared <= cutoff * cutoff
            && graph
                .require_bond_between(left, right)
                .map_err(graph_error)?
                .is_none()
        {
            candidates += 1;
        }
        Ok(())
    };
    // Each pair of cells is visited once: a cell with itself and with the
    // 13 neighbours that follow it.
    let forward = (-1_i64..=1)
        .flat_map(|x| (-1_i64..=1).flat_map(move |y| (-1_i64..=1).map(move |z| [x, y, z])))
        .filter(|offset| *offset > [0, 0, 0])
        .collect::<Vec<_>>();
    // Pair scans poll once per atom whose partners they visit, since one
    // dense cell can hold many atoms.
    let mut scanned = 0usize;
    for &(cell, start, end) in &occupied {
        let members = &order[start..end];
        for (offset, &left) in members.iter().enumerate() {
            poll_every(cancellation, scanned)?;
            scanned += 1;
            for &right in &members[offset + 1..] {
                count(left, right)?;
            }
        }
        for offset in &forward {
            // Cells beyond the supported range hold no atoms.
            let (Some(x), Some(y), Some(z)) = (
                cell[0].checked_add(offset[0]),
                cell[1].checked_add(offset[1]),
                cell[2].checked_add(offset[2]),
            ) else {
                continue;
            };
            let Ok(neighbor) = occupied.binary_search_by_key(&[x, y, z], |&(cell, _, _)| cell)
            else {
                continue;
            };
            let (_, neighbor_start, neighbor_end) = occupied[neighbor];
            for &left in members {
                poll_every(cancellation, scanned)?;
                scanned += 1;
                for &right in &order[neighbor_start..neighbor_end] {
                    count(left, right)?;
                }
            }
        }
    }
    Ok(candidates)
}

fn covalent_bond_cell(point: Point3, line: usize) -> Result<[i64; 3], FormatError> {
    Ok([
        covalent_bond_cell_axis(point.x, "x", line)?,
        covalent_bond_cell_axis(point.y, "y", line)?,
        covalent_bond_cell_axis(point.z, "z", line)?,
    ])
}

fn covalent_bond_cell_axis(coordinate: f64, axis: &str, line: usize) -> Result<i64, FormatError> {
    let cell = (coordinate / COVALENT_BOND_CELL_ANGSTROM).floor();
    if !cell.is_finite() || cell <= i64::MIN as f64 || cell >= i64::MAX as f64 {
        return Err(structure_error(DiagnosticKind::Limit,
            Some(line),
            format!(
                "The _atom_site.Cartn_{axis} coordinate {coordinate} Å is too large to search for bonds."
            ),
        ));
    }
    Ok(cell as i64)
}

fn point_distance_squared(left: Point3, right: Point3) -> f64 {
    let dx = left.x - right.x;
    let dy = left.y - right.y;
    let dz = left.z - right.z;
    dx * dx + dy * dy + dz * dz
}

/// Builds one topology-global hierarchy from source atom identity after every
/// connected molecule component has received its final instance-qualified IDs.
pub(super) fn build_topology_hierarchy(
    instances: &[InstanceProvenance],
    polymer_asym_order: &BTreeMap<String, usize>,
    cancellation: Option<&Cancellation>,
) -> Result<Hierarchy, FormatError> {
    #[derive(Clone)]
    struct ResidueMetadata {
        component_id: String,
        label_component_id: Option<String>,
        author_component_id: Option<String>,
        label_sequence_id: Option<i32>,
        author_sequence_id: Option<String>,
        insertion_code: Option<String>,
        occurrence: Option<usize>,
    }

    let mut hierarchy = Hierarchy::new();
    let mut atoms = instances
        .iter()
        .flat_map(|instance| instance.atoms.iter())
        .collect::<Vec<_>>();
    poll(cancellation)?;
    atoms.sort_by_key(|atom| atom.source_order);
    poll(cancellation)?;
    // Chains' author identities and residues' metadata, each with its index
    // by key; atoms of one residue mostly follow each other.
    let mut chain_authors = Vec::<Option<String>>::new();
    let mut chain_indices = BTreeMap::<&str, usize>::new();
    let mut residue_metadata = Vec::<ResidueMetadata>::new();
    let mut residue_indices = BTreeMap::<(&str, &str), usize>::new();
    // Each atom's residue, and each residue's chain, as indices.
    let mut atom_residues = Vec::with_capacity(atoms.len());
    let mut residue_chains = Vec::new();
    for (step, atom) in atoms.iter().enumerate() {
        poll_every(cancellation, step)?;
        let key = (atom.asym_id.as_str(), atom.residue_key.as_str());
        let known = match atom_residues.last() {
            Some(&last)
                if {
                    let previous = atoms[step - 1];
                    previous.asym_id == atom.asym_id && previous.residue_key == atom.residue_key
                } =>
            {
                Some(last)
            }
            _ => residue_indices.get(&key).copied(),
        };
        let chain = match known {
            Some(residue) => residue_chains[residue],
            None => *chain_indices.entry(key.0).or_insert_with(|| {
                chain_authors.push(None);
                chain_authors.len() - 1
            }),
        };
        merge_optional_source_value(
            &mut chain_authors[chain],
            atom.auth_asym_id.as_ref(),
            "auth_asym_id",
            &atom.asym_id,
            atom.source_line,
        )?;
        if let Some(residue) = known {
            atom_residues.push(residue);
            let metadata = &mut residue_metadata[residue];
            if metadata.component_id != atom.component_id
                || metadata.label_sequence_id != atom.label_sequence_id
                || metadata.insertion_code != atom.insertion_code
                || metadata.occurrence != atom.occurrence
            {
                return Err(inconsistent_residue_metadata(
                    atom,
                    "component, sequence number or insertion code",
                ));
            }
            merge_optional_source_value(
                &mut metadata.label_component_id,
                atom.label_component_id.as_ref(),
                "label_comp_id",
                &atom.residue_key,
                atom.source_line,
            )?;
            merge_optional_source_value(
                &mut metadata.author_component_id,
                atom.auth_component_id.as_ref(),
                "auth_comp_id",
                &atom.residue_key,
                atom.source_line,
            )?;
            merge_optional_source_value(
                &mut metadata.author_sequence_id,
                atom.author_sequence_id.as_ref(),
                "auth_seq_id",
                &atom.residue_key,
                atom.source_line,
            )?;
        } else {
            residue_indices.insert(key, residue_metadata.len());
            atom_residues.push(residue_metadata.len());
            residue_chains.push(chain);
            residue_metadata.push(ResidueMetadata {
                component_id: atom.component_id.clone(),
                label_component_id: atom.label_component_id.clone(),
                author_component_id: atom.auth_component_id.clone(),
                label_sequence_id: atom.label_sequence_id,
                author_sequence_id: atom.author_sequence_id.clone(),
                insertion_code: atom.insertion_code.clone(),
                occurrence: atom.occurrence,
            });
        }
    }
    drop(residue_indices);
    let mut chains = vec![None; chain_authors.len()];
    let mut residues = vec![None; residue_metadata.len()];
    let mut ordered_polymer_chains = polymer_asym_order.iter().collect::<Vec<_>>();
    ordered_polymer_chains.sort_by_key(|(_, order)| **order);
    for (asym_id, _) in ordered_polymer_chains {
        poll(cancellation)?;
        // Every chain with atoms has an author entry.
        let Some(&index) = chain_indices.get(asym_id.as_str()) else {
            continue;
        };
        let chain = hierarchy
            .add_chain(asym_id.clone(), chain_authors[index].clone())
            .map_err(hierarchy_error)?;
        chains[index] = Some(chain);
    }
    for (step, (atom, &index)) in atoms.into_iter().zip(&atom_residues).enumerate() {
        poll_every(cancellation, step)?;
        let residue = if let Some(residue) = residues[index] {
            residue
        } else {
            let chain = match chains[residue_chains[index]] {
                Some(chain) => chain,
                None => {
                    let chain = hierarchy
                        .add_chain(
                            atom.asym_id.clone(),
                            chain_authors[residue_chains[index]].clone(),
                        )
                        .map_err(hierarchy_error)?;
                    chains[residue_chains[index]] = Some(chain);
                    chain
                }
            };
            let metadata = &residue_metadata[index];
            let residue = hierarchy
                .add_residue(
                    chain,
                    metadata.component_id.clone(),
                    metadata.label_sequence_id,
                    metadata.author_sequence_id.clone(),
                    metadata.insertion_code.clone(),
                )
                .map_err(hierarchy_error)?;
            hierarchy
                .set_residue_component_ids(
                    residue,
                    metadata.label_component_id.clone(),
                    metadata.author_component_id.clone(),
                )
                .map_err(hierarchy_error)?;
            residues[index] = Some(residue);
            residue
        };
        hierarchy
            .add_atom_site(
                residue,
                atom.atom,
                AtomSiteMetadata::default()
                    .with_type_symbol(Some(atom.type_symbol.clone()))
                    .with_label_asym_id(atom.label_asym_id.clone())
                    .with_auth_asym_id(atom.auth_asym_id.clone())
                    .with_label_atom_id(atom.label_atom_name.clone())
                    .with_auth_atom_id(atom.auth_atom_name.clone()),
            )
            .map_err(hierarchy_error)?;
    }
    Ok(hierarchy)
}

fn merge_optional_source_value(
    current: &mut Option<String>,
    incoming: Option<&String>,
    field: &'static str,
    identity: &str,
    source_line: usize,
) -> Result<(), FormatError> {
    let Some(incoming) = incoming else {
        return Ok(());
    };
    if let Some(current) = current {
        if current != incoming {
            return Err(structure_error(DiagnosticKind::Inconsistent,
                Some(source_line),
                format!(
                    "The atom sites of {identity:?} have two _atom_site.{field} values: {current:?} and {incoming:?}."
                ),
            ));
        }
    } else {
        *current = Some(incoming.clone());
    }
    Ok(())
}

fn inconsistent_residue_metadata(atom: &AtomProvenance, field: &'static str) -> FormatError {
    structure_error(
        DiagnosticKind::Inconsistent,
        Some(atom.source_line),
        format!(
            "Residue {:?} of structural instance {:?} has different {field} values.",
            atom.residue_key, atom.asym_id
        ),
    )
}

/// Molecule construction failed with `source` as the cause.
pub(super) fn graph_error(source: impl std::error::Error + Send + Sync + 'static) -> FormatError {
    structure_error(
        DiagnosticKind::Chemistry,
        None,
        "Cannot build the molecule.",
    )
    .with_source(source)
}

fn hierarchy_error(source: impl std::error::Error + Send + Sync + 'static) -> FormatError {
    structure_error(
        DiagnosticKind::Inconsistent,
        None,
        "Cannot build the chains and residues.",
    )
    .with_source(source)
}
