//! Interpretation of macromolecular structure sources, mmCIF data blocks and
//! PDB entries, through one pipeline: atom sites, alternate locations,
//! coordinate models, connection records, and connectivity resolution.
mod altloc;
mod assemblies;
mod atom_site;
mod build;
mod connections;
pub(crate) use connections::is_metal;
mod ensemble;
mod mmcif;
mod pdb;
mod sequences;
mod types;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use crate::diagnostics::{Cancellation, DiagnosticKind, FormatError};
use crate::geometry::PeriodicCell;
use crate::structure::ModelBuilder;
use crate::topology::HierarchyBuilder;
use crate::units::{SQUARE_ANGSTROM, SQUARE_NANOMETER};

use super::definition_reuse::DefinitionReuse;
use super::structure_connectivity::ConnectivityCatalog;
use altloc::select_alt_locations;
use atom_site::{atom_key, AtomRow};
use build::{build_molecule, build_topology_hierarchy, graph_error, group_rows};
use connections::{
    resolve_connections, ConnectionRecord, DeclaredConnection, InstanceUnion, ResolvedConnections,
};

pub use altloc::{
    AltLocDecision, AltLocPolicy, AltLocPreference, AltLocResidue, AltLocSelection,
    AltLocSelectionReason, SourceResidueId, SourceResiduePosition,
};
pub use assemblies::{Assembly, AssemblyError, AssemblyGenerator, AssemblyOperation};
pub(crate) use mmcif::{
    alternate_location_inventory, interpret_mmcif, interpret_mmcif_block,
    interpret_mmcif_conformations, interpret_mmcif_conformations_block, interpret_mmcif_ensemble,
    interpret_mmcif_ensemble_block,
};
pub(crate) use pdb::{
    interpret_pdb, interpret_pdb_conformations, interpret_pdb_ensemble, pdb_alternate_locations,
};
pub use sequences::{
    EntityMonomer, PolymerEntity, PolymerEntityInstance, PolymerType, SequencePosition,
    SequenceSource,
};
pub(in crate::io) use types::{poll, poll_every, structure_error};
pub use types::{
    AssemblyOmissionReason, AtomProvenance, BondEvidence, BondProvenance, CellOmissionReason,
    ConflictingSelector, ConnectionResolutionReason, ConnectivityOptions, CoordinationProvenance,
    EntityKind, EntityMetadataPolicy, InstanceProvenance, ModelSelection,
    StructureEnsembleInterpretOptions, StructureEnsembleInterpretation, StructureInterpretIssue,
    StructureInterpretOptions, StructureInterpretation, StructureInterpretationReport,
};

/// Everything interpretation reads from one source, an mmCIF data block or a
/// PDB entry, typed once before any coordinate model or alternative is chosen.
struct StructureSource {
    /// The data block name or PDB entry code, for the report.
    name: String,
    entity_definitions: usize,
    rows: Vec<AtomRow>,
    /// Source-declared alternate-location combinations.
    conformations: BTreeMap<String, BTreeSet<String>>,
    /// Polymer chains in declared order.
    polymer_order: BTreeMap<String, usize>,
    connections: Vec<ConnectionRecord>,
    catalog: ConnectivityCatalog,
    /// The periodic cell of every coordinate model, if the source has one.
    cell: Option<PeriodicCell>,
    /// Issues found while reading the source.
    issues: Vec<StructureInterpretIssue>,
    /// The polymer entities the source states, not yet placed in a model.
    polymer_entities: Vec<PolymerEntity>,
    /// The biological assemblies the source states.
    assemblies: Vec<Assembly>,
}

fn interpret_source(
    source: StructureSource,
    options: &StructureInterpretOptions,
) -> Result<StructureInterpretation, FormatError> {
    let mut prepared = PreparedSource::new(
        source,
        options.altloc_policy(),
        options.cancellation().cloned(),
    )?;
    let model_id = match &options.model_selection() {
        ModelSelection::RequireSingle if prepared.model_sizes.len() == 1 => {
            prepared.model_sizes[0].0.clone()
        }
        ModelSelection::RequireSingle => {
            return Err(structure_error(
                DiagnosticKind::Inconsistent,
                None,
                format!(
                    "The input has {} coordinate models. Select one.",
                    prepared.model_sizes.len(),
                ),
            ))
        }
        ModelSelection::First => prepared.model_sizes[0].0.clone(),
        ModelSelection::Select(id) => id.clone(),
    };
    prepared.interpret_model(&model_id)
}

#[derive(Default)]
struct ModelRows {
    all: Vec<AtomRow>,
    /// Selected rows as positions in `all` with their sites' dense order, in
    /// that order.
    selected: Vec<(usize, usize)>,
}

/// Shared interpretation of one source; each model's rows are consumed once.
struct PreparedSource {
    models: BTreeMap<String, ModelRows>,
    model_sizes: Vec<(String, usize)>,
    report: StructureInterpretationReport,
    polymer_order: BTreeMap<String, usize>,
    connections: Vec<ConnectionRecord>,
    catalog: ConnectivityCatalog,
    cell: Option<PeriodicCell>,
    polymer_entities: Vec<PolymerEntity>,
    cancellation: Option<Cancellation>,
}

impl PreparedSource {
    fn new(
        source: StructureSource,
        altloc_policy: &AltLocPolicy,
        cancellation: Option<Cancellation>,
    ) -> Result<Self, FormatError> {
        let StructureSource {
            name,
            entity_definitions,
            rows,
            conformations,
            polymer_order,
            connections,
            catalog,
            cell,
            issues,
            polymer_entities,
            assemblies,
        } = source;
        let mut report = StructureInterpretationReport {
            block_name: name,
            entity_definitions,
            issues,
            assemblies,
            ..StructureInterpretationReport::default()
        };
        // Alternate-location validation still covers the complete source,
        // including models omitted from the eventual projection.
        let selected = select_alt_locations(
            &rows,
            &conformations,
            altloc_policy,
            &mut report,
            cancellation.as_ref(),
        )?;
        // Models in source order, each source row's model, and each model's
        // row count.
        let mut model_ids = Vec::<String>::new();
        let mut model_sizes = Vec::<usize>::new();
        let mut indices = BTreeMap::<&str, usize>::new();
        let mut row_models = Vec::with_capacity(rows.len());
        for (step, row) in rows.iter().enumerate() {
            poll_every(cancellation.as_ref(), step)?;
            let index = match row_models.last() {
                Some(&index) if model_ids[index] == row.model_id => index,
                _ => *indices.entry(&row.model_id).or_insert_with(|| {
                    model_ids.push(row.model_id.clone());
                    model_sizes.push(0);
                    model_ids.len() - 1
                }),
            };
            model_sizes[index] += 1;
            row_models.push(index);
        }
        drop(indices);
        let mut model_rows = model_sizes
            .iter()
            .map(|&size| ModelRows {
                all: Vec::with_capacity(size),
                selected: Vec::new(),
            })
            .collect::<Vec<_>>();
        // Each source row's position in its model's rows.
        let mut positions = Vec::with_capacity(rows.len());
        for (step, (row, &model)) in rows.into_iter().zip(&row_models).enumerate() {
            poll_every(cancellation.as_ref(), step)?;
            let rows = &mut model_rows[model].all;
            positions.push(rows.len());
            rows.push(row);
        }
        for (step, (row, order)) in selected.into_iter().enumerate() {
            poll_every(cancellation.as_ref(), step)?;
            model_rows[row_models[row]]
                .selected
                .push((positions[row], order));
        }
        let model_sizes = model_ids
            .iter()
            .cloned()
            .zip(model_sizes)
            .collect::<Vec<_>>();
        let models = model_ids
            .into_iter()
            .zip(model_rows)
            .collect::<BTreeMap<_, _>>();
        report.coordinate_models = model_sizes.len();
        Ok(Self {
            models,
            model_sizes,
            report,
            polymer_order,
            connections,
            catalog,
            cell,
            polymer_entities,
            cancellation,
        })
    }

    fn interpret_model(&mut self, model_id: &str) -> Result<StructureInterpretation, FormatError> {
        poll(self.cancellation.as_ref())?;
        let rows = self.models.remove(model_id).ok_or_else(|| {
            structure_error(
                DiagnosticKind::Inconsistent,
                None,
                format!("Coordinate model {model_id:?} does not exist."),
            )
        })?;
        if let Some(row) = rows
            .selected
            .iter()
            .map(|&(position, _)| &rows.all[position])
            .find(|row| row.point.is_none())
        {
            return Err(structure_error(
                DiagnosticKind::InvalidValue,
                Some(row.line),
                format!(
                    "Atom {:?} has no complete position in coordinate model {model_id:?}.",
                    row.atom_name,
                ),
            ));
        }
        let mut report = self.report.clone();
        report.selected_model = Some(model_id.to_owned());
        report
            .alternate_locations
            .retain(|choice| choice.residue.model_id == model_id);
        report.issues.retain(|issue| match issue {
            StructureInterpretIssue::AlternateLocationOmitted { residue, .. }
            | StructureInterpretIssue::AlternateLocationIncomplete { residue, .. }
            | StructureInterpretIssue::AlternateLocationCombinationRejected { residue, .. }
            | StructureInterpretIssue::AlternateLocationOccupancyMissing { residue, .. } => {
                residue.model_id == model_id
            }
            StructureInterpretIssue::AlternateLocationOccupancyTie { residues, .. } => {
                residues.iter().any(|residue| residue.model_id == model_id)
            }
            _ => true,
        });
        for (ignored, count) in &self.model_sizes {
            if ignored != model_id {
                report.ignored_coordinate_models.push(ignored.clone());
                report
                    .issues
                    .push(StructureInterpretIssue::CoordinateModelIgnored {
                        model_id: ignored.clone(),
                        atom_site_rows: *count,
                    });
            }
        }
        publish_model(
            &self.connections,
            rows,
            &self.polymer_order,
            &self.catalog,
            self.cell,
            report,
            self.cancellation.as_ref(),
        )
        .map(|mut interpretation| {
            interpretation.report.polymer_entities = sequences::place_polymer_entities(
                &self.polymer_entities,
                interpretation.model.topology(),
            );
            interpretation
        })
    }
}

fn publish_model(
    connections: &[ConnectionRecord],
    rows: ModelRows,
    polymer_asym_order: &BTreeMap<String, usize>,
    catalog: &ConnectivityCatalog,
    cell: Option<PeriodicCell>,
    mut report: StructureInterpretationReport,
    cancellation: Option<&Cancellation>,
) -> Result<StructureInterpretation, FormatError> {
    let ModelRows {
        all: rows,
        selected,
    } = rows;
    let selected_model = report
        .selected_model
        .as_deref()
        .expect("selected model is recorded")
        .to_owned();
    let selected_rows = selected
        .iter()
        .map(|&(position, _)| &rows[position])
        .collect::<Vec<_>>();
    let mut union = InstanceUnion::new(selected_rows.iter().map(|row| row.instance_key.as_str()));
    let ResolvedConnections {
        bonds: connections,
        coordination,
    } = resolve_connections(
        connections,
        &selected_rows,
        &rows,
        &selected_model,
        &mut union,
        &mut report,
        cancellation,
    )?;
    drop(selected_rows);
    // Selected rows move out of the model's rows, which are not needed again.
    let mut rows = rows.into_iter().map(Some).collect::<Vec<_>>();
    let mut selected = selected
        .into_iter()
        .map(|(position, order)| {
            let mut row = rows[position].take().expect("a row is selected once");
            row.site_order = order;
            row
        })
        .collect::<Vec<_>>();
    drop(rows);
    let groups = group_rows(&selected, &mut union, polymer_asym_order, cancellation)?;
    let mut group_connections = BTreeMap::<String, Vec<&DeclaredConnection>>::new();
    for (step, connection) in connections.iter().enumerate() {
        poll_every(cancellation, step)?;
        group_connections
            .entry(union.find(&connection.instance_key).to_owned())
            .or_default()
            .push(connection);
    }
    let mut builder = ModelBuilder::new();
    let mut definitions = DefinitionReuse::new();
    let mut qualified_atom_data = Vec::new();
    for group in groups {
        poll(cancellation)?;
        let declared = group_connections.remove(group.root()).unwrap_or_default();
        let mut built = build_molecule(group, &mut selected, &declared, cancellation)?;
        built.complete_connectivity(catalog, &mut report, cancellation)?;
        for built in built.publish_components(cancellation)? {
            let labels = built.provenance.atom_labels(built.molecule.atom_count());
            let definition = definitions
                .definition(&mut builder, built.molecule, labels)
                .map_err(graph_error)?;
            let id = builder
                .add_instance(definition, &built.positions)
                .map_err(graph_error)?;
            poll(cancellation)?;
            let (provenance, atom_data) = built.provenance.qualify(id);
            report.instances.push(provenance);
            qualified_atom_data.extend(atom_data);
        }
    }
    poll(cancellation)?;
    if !coordination.is_empty() {
        // An atom's key ends in its name, so only atoms whose name ends as a
        // partner's key does can be partners.
        let suffixes = coordination
            .iter()
            .flat_map(|record| &record.atoms)
            .filter_map(|key| key.rsplit('|').next())
            .collect::<HashSet<_>>();
        let atoms = report
            .instances
            .iter()
            .flat_map(|instance| instance.atoms.iter())
            .filter(|atom| {
                atom.atom_name
                    .rsplit('|')
                    .next()
                    .is_some_and(|suffix| suffixes.contains(suffix))
            })
            .map(|atom| {
                (
                    atom_key(&atom.asym_id, &atom.residue_key, &atom.atom_name),
                    atom.atom,
                )
            })
            .collect::<HashMap<_, _>>();
        report.coordination = coordination
            .into_iter()
            .map(|record| CoordinationProvenance {
                atoms: record.atoms.map(|key| {
                    *atoms
                        .get(&key)
                        .expect("every selected row is an interpreted atom")
                }),
                connection_id: record.connection_id,
                connection_type: record.connection_type,
                source_line: record.line,
            })
            .collect();
    }
    poll(cancellation)?;
    // Dense atom order follows the selected atom-site rows, so coordinate
    // files in this row order address the same atoms.
    let mut rows = report
        .instances
        .iter()
        .flat_map(|instance| instance.atoms.iter())
        .map(|atom| (atom.source_order, atom.atom))
        .collect::<Vec<_>>();
    rows.sort_unstable_by_key(|&(row, _)| row);
    builder
        .set_atom_order(rows.into_iter().map(|(_, atom)| atom))
        .map_err(graph_error)?;
    poll(cancellation)?;
    *builder.hierarchy_mut() = HierarchyBuilder::from_storage(build_topology_hierarchy(
        &report.instances,
        polymer_asym_order,
        cancellation,
    )?);
    poll(cancellation)?;
    let mut model = builder.build().map_err(graph_error)?;
    poll(cancellation)?;
    model.set_cell(cell);
    let mut occupancies = vec![None; model.atom_count()];
    let mut b_factors = vec![None; model.atom_count()];
    let b_factor_scale = SQUARE_ANGSTROM
        .conversion_factor_to(SQUARE_NANOMETER)
        .map_err(graph_error)?;
    for (atom, occupancy, b_factor) in qualified_atom_data {
        let index = model.topology().atom_index(atom).ok_or_else(|| {
            structure_error(
                DiagnosticKind::Chemistry,
                None,
                format!("The topology has no {atom}."),
            )
        })?;
        occupancies[index.index()] = occupancy;
        b_factors[index.index()] = b_factor.map(|value| value * b_factor_scale);
    }
    model
        .set_occupancies(Some(occupancies))
        .map_err(graph_error)?;
    model
        .set_b_factors(Some(crate::units::Quantity::new(
            b_factors,
            SQUARE_NANOMETER,
        )))
        .map_err(graph_error)?;
    report.macromolecules = report
        .instances
        .iter()
        .filter(|instance| instance.entity_kinds().iter().any(EntityKind::is_macro))
        .count();
    report.small_molecules = report.instances.len() - report.macromolecules;
    report.solvent_molecules = report
        .instances
        .iter()
        .filter(|instance| instance.entity_kinds().contains(&EntityKind::Water))
        .count();
    Ok(StructureInterpretation { model, report })
}
