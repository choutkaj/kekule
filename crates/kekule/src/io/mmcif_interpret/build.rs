use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::bio::{MacroMolecule, SmcraAtomSiteMetadata, SmcraHierarchy};
use crate::core::{Atom, AtomId, Conformer, ConformerId, Molecule};
use crate::geometry::Point3;
use crate::small::model::SmallMolecule;
use crate::topology::{InstanceAtomId, MoleculeInstanceId, MoleculeInstanceMetadata, MoleculeRole};
use crate::units::{Quantity, ANGSTROM};

use super::super::MmcifDataBlock;
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

pub(super) fn polymer_asym_order(block: &MmcifDataBlock) -> BTreeMap<String, usize> {
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

pub(super) enum BuiltMolecule {
    Small {
        molecule: SmallMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
        provenance: BuiltMoleculeProvenance,
    },
    Macro {
        molecule: MacroMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
        provenance: BuiltMoleculeProvenance,
    },
}

pub(super) struct BuiltMoleculeProvenance {
    coordinate_model_id: String,
    asym_ids: Vec<String>,
    entity_ids: Vec<String>,
    entity_kinds: Vec<MmcifEntityKind>,
    atoms: Vec<BuiltAtomProvenance>,
}

struct BuiltAtomProvenance {
    atom: AtomId,
    source_line: usize,
    atom_site_id: Option<String>,
    atom_name: String,
    component_id: String,
    asym_id: String,
    auth_asym_id: Option<String>,
    entity_id: Option<String>,
    label_sequence_id: Option<i32>,
    author_sequence_id: Option<String>,
    insertion_code: Option<String>,
    occurrence: Option<usize>,
    selected_alternate_location: Option<String>,
    occupancy: Option<f64>,
    b_factor: Option<f64>,
}

pub(super) type QualifiedAtomData = Vec<(InstanceAtomId, Option<f64>, Option<f64>)>;

impl BuiltMoleculeProvenance {
    pub(super) fn qualify(
        self,
        molecule: MoleculeInstanceId,
    ) -> (MmcifInstanceProvenance, QualifiedAtomData) {
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
                        source_line: atom.source_line,
                        atom_site_id: atom.atom_site_id,
                        atom_name: atom.atom_name,
                        component_id: atom.component_id,
                        asym_id: atom.asym_id,
                        auth_asym_id: atom.auth_asym_id,
                        entity_id: atom.entity_id,
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

pub(super) fn build_molecule(
    group: MoleculeGroup,
    connections: &[DeclaredConnection],
    report: &mut MmcifInterpretationReport,
) -> Result<BuiltMolecule, MmcifInterpretError> {
    let is_macro = group.kinds.iter().any(MmcifEntityKind::is_macro);
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
        atoms.insert(key.clone(), graph.add_atom(atom).map_err(graph_error)?);
    }
    let model_id = group
        .rows
        .first()
        .map(|row| row.model_id.clone())
        .ok_or_else(|| MmcifInterpretError::new(None, "empty molecule group"))?;
    let mut conformer = Conformer::new(ANGSTROM).expect("angstrom is a length unit");
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
        conformer
            .set_position(atoms[&row.atom_key()], Quantity::new(point, ANGSTROM))
            .expect("matching coordinate units");
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
                .add_bond(left, right, connection.order)
                .map_err(graph_error)?;
        } else {
            let existing = graph
                .bond_between(left, right)
                .map_err(graph_error)?
                .expect("existing bond was found");
            if graph.bond(existing).map_err(graph_error)?.order != connection.order {
                return Err(MmcifInterpretError::new(
                    None,
                    "duplicate struct_conn records assign conflicting bond orders",
                ));
            }
        }
    }
    let connectivity_candidates = infer_covalent_bonds(&graph, &representative, &atoms)?;
    if connectivity_candidates > 0 {
        report.connectivity_candidates += connectivity_candidates;
        report
            .issues
            .push(MmcifInterpretIssue::ConnectivityCandidatesInferred {
                atom_count: graph.atom_count(),
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
            source_line: row.line,
            atom_site_id: row.atom_site_id.clone(),
            atom_name: row.atom_name.clone(),
            component_id: row.comp_id.clone(),
            asym_id: row.asym_id.clone(),
            auth_asym_id: row.auth_asym_id.clone(),
            entity_id: row.entity_id.clone(),
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
    if graph.atom_count() > 1 {
        report.template_bonds_pending += 1;
    }

    let mut metadata = MoleculeInstanceMetadata::default();
    for kind in &group.kinds {
        match kind {
            MmcifEntityKind::Polymer => {
                metadata.insert_role(MoleculeRole::Polymer);
            }
            MmcifEntityKind::Branched => {
                metadata.insert_role(MoleculeRole::Branched);
            }
            MmcifEntityKind::NonPolymer => {
                metadata.insert_role(MoleculeRole::NonPolymer);
            }
            MmcifEntityKind::Water => {
                metadata.insert_role(MoleculeRole::Solvent);
            }
            MmcifEntityKind::Other(_) => {}
        }
    }
    if graph.atom_count() == 1
        && graph
            .atoms()
            .next()
            .is_some_and(|(_, atom)| atom.formal_charge != 0)
    {
        metadata.insert_role(MoleculeRole::Ion);
    }
    let conformer = graph
        .add_conformer(conformer)
        .expect("interpreted coordinates reference live atoms");
    if is_macro {
        let hierarchy = build_hierarchy(&graph, &representative, &atoms)?;
        Ok(BuiltMolecule::Macro {
            molecule: MacroMolecule::from_parts_unchecked_connectedness(graph, hierarchy)
                .map_err(graph_error)?,
            conformer,
            metadata,
            provenance,
        })
    } else {
        let molecule = SmallMolecule::from_molecule_unchecked_connectedness(graph);
        Ok(BuiltMolecule::Small {
            molecule,
            conformer,
            metadata,
            provenance,
        })
    }
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

fn build_hierarchy(
    graph: &Molecule,
    representative: &[(String, AtomRow)],
    atoms: &BTreeMap<String, AtomId>,
) -> Result<SmcraHierarchy, MmcifInterpretError> {
    let mut hierarchy = SmcraHierarchy::new();
    let mut chains = BTreeMap::new();
    let mut residues = BTreeMap::new();
    for (key, row) in representative {
        let chain = if let Some(chain) = chains.get(&row.asym_id) {
            *chain
        } else {
            let chain = hierarchy
                .add_chain(row.asym_id.clone(), row.auth_asym_id.clone())
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
                SmcraAtomSiteMetadata {
                    type_symbol: Some(row.element.symbol().to_owned()),
                    label_asym_id: Some(row.asym_id.clone()),
                    auth_asym_id: row.auth_asym_id.clone(),
                    label_atom_id: Some(row.atom_name.clone()),
                    auth_atom_id: row.auth_atom_name.clone(),
                },
            )
            .map_err(hierarchy_error)?;
    }
    Ok(hierarchy)
}

pub(super) fn graph_error(error: impl fmt::Display) -> MmcifInterpretError {
    MmcifInterpretError::new(None, error.to_string())
}

fn hierarchy_error(error: impl fmt::Display) -> MmcifInterpretError {
    MmcifInterpretError::new(None, error.to_string())
}
