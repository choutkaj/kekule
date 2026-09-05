use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::chemistry::canonicalize_molecule_for_publication;
use crate::core::{Atom, AtomId, Molecule, MoleculeEditor};
use crate::geometry::Point3;
use crate::structure::Positions;
use crate::topology::{AtomSiteMetadata, Hierarchy};
use crate::topology::{InstanceAtomId, MoleculeInstanceId};
use crate::units::{Quantity, ANGSTROM};

use super::super::staged_coordinates::StagedCoordinates;
use super::super::MmcifBlock;
use super::atom_site::{optional, AtomRow};
use super::struct_conn::{DeclaredConnection, InstanceUnion};
use super::types::{
    MmcifAtomProvenance, MmcifEntityKind, MmcifInstanceProvenance, MmcifInterpretError,
    MmcifInterpretIssue, MmcifInterpretationReport,
};

#[derive(Debug)]
pub(super) struct MoleculeGroup {
    rows: Vec<AtomRow>,
    kinds: BTreeSet<MmcifEntityKind>,
    instance_keys: BTreeSet<String>,
}

pub(super) fn polymer_asym_order(block: &MmcifBlock) -> BTreeMap<String, usize> {
    let mut order = BTreeMap::new();
    let Some(table) = block.loop_with_tag("_pdbx_poly_seq_scheme.asym_id") else {
        return order;
    };
    for row in 0..table.row_count() {
        if let Some(asym_id) = optional(table, row, "_pdbx_poly_seq_scheme.asym_id") {
            let next = order.len();
            order.entry(asym_id.to_owned()).or_insert(next);
        }
    }
    order
}

pub(super) fn group_rows(
    rows: Vec<AtomRow>,
    union: &mut InstanceUnion,
    polymer_asym_order: &BTreeMap<String, usize>,
) -> Vec<MoleculeGroup> {
    let mut group_indices = BTreeMap::new();
    let mut groups = Vec::new();
    for row in rows {
        let root = union.find(&row.instance_key);
        let index = *group_indices.entry(root).or_insert_with(|| {
            groups.push(MoleculeGroup {
                rows: Vec::new(),
                kinds: BTreeSet::new(),
                instance_keys: BTreeSet::new(),
            });
            groups.len() - 1
        });
        let group = &mut groups[index];
        group.kinds.insert(row.kind.clone());
        group.instance_keys.insert(row.instance_key.clone());
        group.rows.push(row);
    }
    for group in &mut groups {
        group.rows.sort_by_key(|row| {
            polymer_asym_order
                .get(&row.asym_id)
                .copied()
                .unwrap_or(usize::MAX)
        });
    }
    groups.sort_by_key(|group| {
        group
            .rows
            .iter()
            .filter_map(|row| polymer_asym_order.get(&row.asym_id).copied())
            .min()
            .unwrap_or(usize::MAX)
    });
    groups
}

pub(super) struct BuiltMolecule {
    pub(super) editor: MoleculeEditor,
    pub(super) coordinates: StagedCoordinates,
    pub(super) provenance: BuiltMoleculeProvenance,
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
    entity_kinds: Vec<MmcifEntityKind>,
    atoms: Vec<BuiltAtomProvenance>,
}

#[derive(Clone)]
struct BuiltAtomProvenance {
    atom: AtomId,
    type_symbol: String,
    source_line: usize,
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
    entity_kind: MmcifEntityKind,
    residue_key: String,
    label_sequence_id: Option<i32>,
    author_sequence_id: Option<String>,
    insertion_code: Option<String>,
    occurrence: Option<usize>,
    selected_alternate_location: Option<String>,
    occupancy: Option<f64>,
    b_factor: Option<f64>,
}

pub(super) type QualifiedAtomProperties = Vec<(InstanceAtomId, Option<f64>, Option<f64>)>;

impl BuiltMoleculeProvenance {
    pub(super) fn qualify(
        self,
        molecule: MoleculeInstanceId,
    ) -> (MmcifInstanceProvenance, QualifiedAtomProperties) {
        let mut atom_data = Vec::with_capacity(self.atoms.len());
        let provenance = MmcifInstanceProvenance {
            molecule,
            coordinate_model_id: self.coordinate_model_id,
            asym_ids: self.asym_ids,
            entity_ids: self.entity_ids,
            entity_kinds: self.entity_kinds,
            atoms: self
                .atoms
                .into_iter()
                .map(|atom| {
                    let qualified = InstanceAtomId::new(molecule, atom.atom);
                    atom_data.push((qualified, atom.occupancy, atom.b_factor));
                    MmcifAtomProvenance {
                        atom: qualified,
                        type_symbol: atom.type_symbol,
                        source_line: atom.source_line,
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
        catalog: &super::super::mmcif_connectivity::ConnectivityCatalog,
    ) -> Result<(), MmcifInterpretError> {
        let Self {
            editor, provenance, ..
        } = self;
        let atoms = provenance.atoms.iter().map(|atom| {
            super::super::mmcif_connectivity::StagedAtomProvenance {
                atom: atom.atom,
                atom_name: &atom.atom_name,
                component_id: &atom.component_id,
                asym_id: &atom.asym_id,
                entity_id: atom.entity_id.as_deref(),
                label_sequence_id: atom.label_sequence_id,
                author_sequence_id: atom.author_sequence_id.as_deref(),
                insertion_code: atom.insertion_code.as_deref(),
                occurrence: atom.occurrence,
            }
        });
        super::super::mmcif_connectivity::complete_editor_connectivity(catalog, editor, atoms)
    }

    pub(super) fn publish_components(self) -> Result<Vec<PublishedMolecule>, MmcifInterpretError> {
        let components = self.editor.working().connected_components();
        if components.is_empty() {
            return Err(graph_error("mmCIF molecule group has no atoms"));
        }
        let mut published = Vec::with_capacity(components.len());
        for component in components {
            let selected = component.iter().copied().collect::<BTreeSet<_>>();
            let mut editor = crate::core::MoleculeEditor::new();
            let mut atom_map = BTreeMap::new();
            for source_atom in component.iter().copied() {
                let atom = self
                    .editor
                    .working()
                    .atom(source_atom)
                    .map_err(graph_error)?
                    .clone();
                let target_atom = editor.add_atom(atom).map_err(graph_error)?;
                atom_map.insert(source_atom, target_atom);
            }
            for (_, source_bond) in self.editor.working().bonds() {
                let (left, right) = source_bond.endpoints();
                if !selected.contains(&left) || !selected.contains(&right) {
                    continue;
                }
                let target = editor
                    .add_bond(atom_map[&left], atom_map[&right], source_bond.order)
                    .map_err(graph_error)?;
                let _ = target;
            }
            let bond_sources = self
                .editor
                .working()
                .bonds()
                .filter(|(_, bond)| selected.contains(&bond.a()) && selected.contains(&bond.b()))
                .map(|(id, _)| id.index())
                .collect::<Vec<_>>();
            let mut properties =
                crate::properties::Properties::molecule(component.len(), bond_sources.len());
            *properties.atoms_mut() = self
                .editor
                .working()
                .properties()
                .atoms()
                .select_indices(&component.iter().map(|id| id.index()).collect::<Vec<_>>())
                .map_err(graph_error)?;
            *properties.bonds_mut() = self
                .editor
                .working()
                .properties()
                .bonds()
                .select_indices(&bond_sources)
                .map_err(graph_error)?;
            *editor.working_mut().properties_mut() = properties;
            let mut coordinates =
                StagedCoordinates::with_atom_capacity(component.len(), self.coordinates.unit())
                    .map_err(graph_error)?;
            for (source_atom, target_atom) in &atom_map {
                if let Some(point) = self.coordinates.position(*source_atom) {
                    coordinates
                        .set_position(*target_atom, point)
                        .map_err(graph_error)?;
                }
            }
            canonicalize_molecule_for_publication(editor.working_mut(), Some(&coordinates), &[])
                .map_err(graph_error)?;
            let molecule = editor.finish().map_err(graph_error)?;
            let positions = coordinates.to_positions(&molecule).map_err(graph_error)?;
            let mut provenance = self.provenance.clone();
            provenance
                .atoms
                .retain(|atom| selected.contains(&atom.atom));
            for atom in &mut provenance.atoms {
                atom.atom = atom_map[&atom.atom];
            }
            published.push(PublishedMolecule {
                molecule,
                positions,
                provenance,
            });
        }
        Ok(published)
    }
}

pub(super) fn build_molecule(
    group: MoleculeGroup,
    connections: &[DeclaredConnection],
    report: &mut MmcifInterpretationReport,
) -> Result<BuiltMolecule, MmcifInterpretError> {
    let mut editor = crate::core::MoleculeEditor::new();
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
        atoms.insert(key.clone(), editor.add_atom(atom).map_err(graph_error)?);
    }
    let model_id = group
        .rows
        .first()
        .map(|row| row.model_id.clone())
        .ok_or_else(|| MmcifInterpretError::new(None, "empty molecule group"))?;
    let mut coordinates =
        StagedCoordinates::with_atom_capacity(atoms.len(), ANGSTROM).map_err(graph_error)?;
    for row in &group.rows {
        let point = row.point.ok_or_else(|| {
            MmcifInterpretError::new(
                Some(row.line),
                format!(
                    "missing selected-model position for atom `{}`",
                    row.atom_name
                ),
            )
        })?;
        coordinates
            .set_position(atoms[&row.atom_key()], Quantity::new(point, ANGSTROM))
            .map_err(graph_error)?;
    }
    for connection in connections {
        let Some(&left) = atoms.get(&connection.left_atom) else {
            continue;
        };
        let Some(&right) = atoms.get(&connection.right_atom) else {
            continue;
        };
        if editor
            .working()
            .bond_between(left, right)
            .map_err(graph_error)?
            .is_none()
        {
            editor
                .add_bond(left, right, connection.order)
                .map_err(graph_error)?;
        } else {
            let existing = editor
                .working()
                .bond_between(left, right)
                .map_err(graph_error)?
                .expect("existing bond was found");
            if editor.working().bond(existing).map_err(graph_error)?.order != connection.order {
                return Err(MmcifInterpretError::new(
                    None,
                    "duplicate struct_conn records assign conflicting bond orders",
                ));
            }
        }
    }
    let connectivity_candidates =
        infer_covalent_bonds(editor.working_mut(), &representative, &atoms)?;
    if connectivity_candidates > 0 {
        report.connectivity_candidates += connectivity_candidates;
        report
            .issues
            .push(MmcifInterpretIssue::ConnectivityCandidatesInferred {
                atom_count: editor.working().atom_count(),
                candidate_count: connectivity_candidates,
            });
    }
    let asym_ids = representative
        .iter()
        .map(|(_, row)| row)
        .map(|row| row.asym_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let entity_ids = representative
        .iter()
        .map(|(_, row)| row)
        .filter_map(|row| row.entity_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let entity_kinds = group.kinds.iter().cloned().collect::<Vec<_>>();
    let atom_provenance = representative
        .iter()
        .map(|(key, row)| BuiltAtomProvenance {
            atom: atoms[key],
            type_symbol: row.element.symbol().to_owned(),
            source_line: row.line,
            atom_site_id: row.atom_site_id.clone(),
            label_atom_name: row.label_atom_name.clone(),
            atom_name: row.atom_name.clone(),
            auth_atom_name: row.auth_atom_name.clone(),
            label_component_id: row.label_comp_id.clone(),
            component_id: row.comp_id.clone(),
            auth_component_id: row.auth_comp_id.clone(),
            label_asym_id: row.label_asym_id.clone(),
            asym_id: row.asym_id.clone(),
            auth_asym_id: row.auth_asym_id.clone(),
            entity_id: row.entity_id.clone(),
            entity_kind: row.kind.clone(),
            residue_key: row.residue_key.clone(),
            label_sequence_id: row.label_seq_id,
            author_sequence_id: row.auth_seq_id.clone(),
            insertion_code: row.insertion_code.clone(),
            occurrence: row.occurrence,
            selected_alternate_location: row.alt_id.clone(),
            occupancy: row.occupancy,
            b_factor: row.b_factor,
        })
        .collect();
    let provenance = BuiltMoleculeProvenance {
        coordinate_model_id: model_id,
        asym_ids,
        entity_ids,
        entity_kinds,
        atoms: atom_provenance,
    };
    Ok(BuiltMolecule {
        editor,
        coordinates,
        provenance,
    })
}

const COVALENT_BOND_CELL_ANGSTROM: f64 = 2.1;
const COVALENT_BOND_TOLERANCE_ANGSTROM: f64 = 0.45;
const MIN_COVALENT_BOND_DISTANCE_SQUARED: f64 = 0.16;
const FALLBACK_COVALENT_RADIUS_ANGSTROM: f64 = 0.77;

fn infer_covalent_bonds(
    graph: &Molecule,
    representative: &[(String, AtomRow)],
    atoms: &BTreeMap<String, AtomId>,
) -> Result<usize, MmcifInterpretError> {
    let mut cells = BTreeMap::<[i64; 3], Vec<usize>>::new();
    let mut candidates = 0usize;

    for (right_index, (right_key, right_row)) in representative.iter().enumerate() {
        let right_point = right_row
            .point
            .expect("selected mmCIF atom rows have complete positions");
        let right_cell = covalent_bond_cell(right_point, right_row.line)?;
        for offset_x in -1_i64..=1 {
            for offset_y in -1_i64..=1 {
                for offset_z in -1_i64..=1 {
                    let neighbor = covalent_bond_neighbor(
                        right_cell,
                        [offset_x, offset_y, offset_z],
                        right_row.line,
                    )?;
                    let Some(left_indexes) = cells.get(&neighbor) else {
                        continue;
                    };
                    for &left_index in left_indexes {
                        let (left_key, left_row) = &representative[left_index];
                        let left_point = left_row
                            .point
                            .expect("selected mmCIF atom rows have complete positions");
                        let distance_squared = point_distance_squared(left_point, right_point);
                        if distance_squared <= MIN_COVALENT_BOND_DISTANCE_SQUARED {
                            continue;
                        }
                        let left_radius = left_row
                            .element
                            .covalent_radius_angstrom()
                            .unwrap_or(FALLBACK_COVALENT_RADIUS_ANGSTROM);
                        let right_radius = right_row
                            .element
                            .covalent_radius_angstrom()
                            .unwrap_or(FALLBACK_COVALENT_RADIUS_ANGSTROM);
                        let cutoff =
                            (left_radius + right_radius + COVALENT_BOND_TOLERANCE_ANGSTROM)
                                .min(COVALENT_BOND_CELL_ANGSTROM);
                        if distance_squared > cutoff * cutoff {
                            continue;
                        }

                        let left = atoms[left_key];
                        let right = atoms[right_key];
                        if graph
                            .bond_between(left, right)
                            .map_err(graph_error)?
                            .is_none()
                        {
                            candidates += 1;
                        }
                    }
                }
            }
        }
        cells.entry(right_cell).or_default().push(right_index);
    }

    Ok(candidates)
}

fn covalent_bond_cell(point: Point3, line: usize) -> Result<[i64; 3], MmcifInterpretError> {
    Ok([
        covalent_bond_cell_axis(point.x, "x", line)?,
        covalent_bond_cell_axis(point.y, "y", line)?,
        covalent_bond_cell_axis(point.z, "z", line)?,
    ])
}

fn covalent_bond_cell_axis(
    coordinate: f64,
    axis: &str,
    line: usize,
) -> Result<i64, MmcifInterpretError> {
    let cell = (coordinate / COVALENT_BOND_CELL_ANGSTROM).floor();
    if !cell.is_finite() || cell <= i64::MIN as f64 || cell >= i64::MAX as f64 {
        return Err(MmcifInterpretError::new(
            Some(line),
            format!(
                "_atom_site.Cartn_{axis} coordinate is outside the supported covalent-connectivity diagnostic cell range"
            ),
        ));
    }
    Ok(cell as i64)
}

fn covalent_bond_neighbor(
    cell: [i64; 3],
    offset: [i64; 3],
    line: usize,
) -> Result<[i64; 3], MmcifInterpretError> {
    let checked = |axis: usize| {
        cell[axis].checked_add(offset[axis]).ok_or_else(|| {
            MmcifInterpretError::new(
                Some(line),
                "atom-site coordinate exceeds the covalent-connectivity diagnostic neighbor range",
            )
        })
    };
    Ok([checked(0)?, checked(1)?, checked(2)?])
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
    instances: &[MmcifInstanceProvenance],
    polymer_asym_order: &BTreeMap<String, usize>,
) -> Result<Hierarchy, MmcifInterpretError> {
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
    let mut chains = BTreeMap::new();
    let mut residues = BTreeMap::new();
    let mut atoms = instances
        .iter()
        .flat_map(|instance| instance.atoms.iter())
        .collect::<Vec<_>>();
    atoms.sort_by_key(|atom| atom.source_line);
    let mut chain_authors = BTreeMap::<String, Option<String>>::new();
    let mut residue_metadata = BTreeMap::<(String, String), ResidueMetadata>::new();
    for atom in &atoms {
        merge_optional_source_value(
            chain_authors.entry(atom.asym_id.clone()).or_default(),
            atom.auth_asym_id.as_ref(),
            "auth_asym_id",
            &atom.asym_id,
            atom.source_line,
        )?;
        let key = (atom.asym_id.clone(), atom.residue_key.clone());
        if let Some(metadata) = residue_metadata.get_mut(&key) {
            if metadata.component_id != atom.component_id
                || metadata.label_sequence_id != atom.label_sequence_id
                || metadata.insertion_code != atom.insertion_code
                || metadata.occurrence != atom.occurrence
            {
                return Err(inconsistent_residue_metadata(
                    atom,
                    "canonical residue identity",
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
            residue_metadata.insert(
                key,
                ResidueMetadata {
                    component_id: atom.component_id.clone(),
                    label_component_id: atom.label_component_id.clone(),
                    author_component_id: atom.auth_component_id.clone(),
                    label_sequence_id: atom.label_sequence_id,
                    author_sequence_id: atom.author_sequence_id.clone(),
                    insertion_code: atom.insertion_code.clone(),
                    occurrence: atom.occurrence,
                },
            );
        }
    }
    let mut ordered_polymer_chains = polymer_asym_order.iter().collect::<Vec<_>>();
    ordered_polymer_chains.sort_by_key(|(_, order)| **order);
    for (asym_id, _) in ordered_polymer_chains {
        let Some(_) = atoms.iter().find(|atom| &atom.asym_id == asym_id) else {
            continue;
        };
        let chain = hierarchy
            .add_chain(asym_id.clone(), chain_authors[asym_id].clone())
            .map_err(hierarchy_error)?;
        chains.insert(asym_id.clone(), chain);
    }
    for atom in atoms {
        let chain = if let Some(chain) = chains.get(&atom.asym_id) {
            *chain
        } else {
            let chain = hierarchy
                .add_chain(atom.asym_id.clone(), chain_authors[&atom.asym_id].clone())
                .map_err(hierarchy_error)?;
            chains.insert(atom.asym_id.clone(), chain);
            chain
        };
        let residue_key = (atom.asym_id.clone(), atom.residue_key.clone());
        let residue = if let Some(residue) = residues.get(&residue_key) {
            *residue
        } else {
            let metadata = &residue_metadata[&residue_key];
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
            residues.insert(residue_key, residue);
            residue
        };
        hierarchy
            .add_atom_site(
                residue,
                atom.atom,
                AtomSiteMetadata {
                    type_symbol: Some(atom.type_symbol.clone()),
                    label_asym_id: atom.label_asym_id.clone(),
                    auth_asym_id: atom.auth_asym_id.clone(),
                    label_atom_id: atom.label_atom_name.clone(),
                    auth_atom_id: atom.auth_atom_name.clone(),
                },
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
) -> Result<(), MmcifInterpretError> {
    let Some(incoming) = incoming else {
        return Ok(());
    };
    if let Some(current) = current {
        if current != incoming {
            return Err(MmcifInterpretError::new(
                Some(source_line),
                format!(
                    "canonical hierarchy identity `{identity}` has conflicting _atom_site.{field} values `{current}` and `{incoming}`"
                ),
            ));
        }
    } else {
        *current = Some(incoming.clone());
    }
    Ok(())
}

fn inconsistent_residue_metadata(
    atom: &MmcifAtomProvenance,
    field: &'static str,
) -> MmcifInterpretError {
    MmcifInterpretError::new(
        Some(atom.source_line),
        format!(
            "canonical residue `{}` in asymmetry `{}` has inconsistent {field}",
            atom.residue_key, atom.asym_id
        ),
    )
}

pub(super) fn graph_error(error: impl fmt::Display) -> MmcifInterpretError {
    MmcifInterpretError::new(None, error.to_string())
}

fn hierarchy_error(error: impl fmt::Display) -> MmcifInterpretError {
    MmcifInterpretError::new(None, error.to_string())
}
