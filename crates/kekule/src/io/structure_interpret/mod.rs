//! Interpretation of macromolecular structure sources, mmCIF data blocks and
//! PDB entries, through one pipeline: atom sites, alternate locations,
//! coordinate models, connection records, and connectivity resolution.
mod altloc;
mod atom_site;
mod build;
mod connections;
mod ensemble;
mod mmcif;
mod types;

use std::collections::{BTreeMap, BTreeSet};

use crate::diagnostics::DiagnosticKind;
use crate::structure::ModelBuilder;
use crate::units::{SQUARE_ANGSTROM, SQUARE_NANOMETER};

use super::structure_connectivity::ConnectivityCatalog;
use altloc::select_alt_locations;
use atom_site::AtomRow;
use build::{build_molecule, build_topology_hierarchy, graph_error, group_rows};
use connections::{resolve_connections, ConnectionRecord, DeclaredConnection, InstanceUnion};

pub use altloc::{
    AltLocDecision, AltLocPolicy, AltLocPreference, AltLocResidue, AltLocSelection,
    AltLocSelectionReason, SourceResidueId, SourceResiduePosition,
};
pub(crate) use mmcif::{
    alternate_location_inventory, interpret_mmcif, interpret_mmcif_block,
    interpret_mmcif_conformations, interpret_mmcif_conformations_block, interpret_mmcif_ensemble,
    interpret_mmcif_ensemble_block,
};
pub use types::{
    AtomProvenance, BondEvidence, BondProvenance, ConnectionResolutionReason, ConnectivityOptions,
    EntityKind, InstanceProvenance, ModelSelection, StructureEnsembleInterpretError,
    StructureEnsembleInterpretOptions, StructureEnsembleInterpretation, StructureInterpretError,
    StructureInterpretIssue, StructureInterpretOptions, StructureInterpretation,
    StructureInterpretationReport,
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
    /// Issues found while reading the source.
    issues: Vec<StructureInterpretIssue>,
}

fn interpret_source(
    source: StructureSource,
    options: &StructureInterpretOptions,
) -> Result<StructureInterpretation, StructureInterpretError> {
    let mut prepared = PreparedSource::new(source, &options.altloc_policy)?;
    let model_id = match &options.model_selection {
        ModelSelection::RequireSingle if prepared.model_sizes.len() == 1 => {
            prepared.model_sizes[0].0.clone()
        }
        ModelSelection::RequireSingle => {
            return Err(StructureInterpretError::new(
                DiagnosticKind::Inconsistent,
                None,
                format!(
                    "coordinate data contains {} models; select one explicitly",
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
    selected: Vec<AtomRow>,
}

/// Shared interpretation of one source; each model's rows are consumed once.
struct PreparedSource {
    models: BTreeMap<String, ModelRows>,
    model_sizes: Vec<(String, usize)>,
    report: StructureInterpretationReport,
    polymer_order: BTreeMap<String, usize>,
    connections: Vec<ConnectionRecord>,
    catalog: ConnectivityCatalog,
}

impl PreparedSource {
    fn new(
        source: StructureSource,
        altloc_policy: &AltLocPolicy,
    ) -> Result<Self, StructureInterpretError> {
        let StructureSource {
            name,
            entity_definitions,
            rows,
            conformations,
            polymer_order,
            connections,
            catalog,
            issues,
        } = source;
        let mut report = StructureInterpretationReport {
            block_name: name,
            entity_definitions,
            issues,
            ..StructureInterpretationReport::default()
        };
        // Alternate-location validation still covers the complete source,
        // including models omitted from the eventual projection.
        let selected = select_alt_locations(&rows, &conformations, altloc_policy, &mut report)?;
        let mut models = BTreeMap::<String, ModelRows>::new();
        let mut model_ids = Vec::new();
        for row in rows {
            let model = models.entry(row.model_id.clone()).or_default();
            if model.all.is_empty() {
                model_ids.push(row.model_id.clone());
            }
            model.all.push(row);
        }
        for row in selected {
            let model = models
                .get_mut(&row.model_id)
                .expect("selected row has a source model");
            model.selected.push(row);
        }
        let model_sizes = model_ids
            .into_iter()
            .map(|id| {
                let count = models[&id].all.len();
                (id, count)
            })
            .collect::<Vec<_>>();
        report.coordinate_models = model_sizes.len();
        Ok(Self {
            models,
            model_sizes,
            report,
            polymer_order,
            connections,
            catalog,
        })
    }

    fn interpret_model(
        &mut self,
        model_id: &str,
    ) -> Result<StructureInterpretation, StructureInterpretError> {
        let rows = self.models.remove(model_id).ok_or_else(|| {
            StructureInterpretError::new(
                DiagnosticKind::Inconsistent,
                None,
                format!("coordinate model `{model_id}` is unavailable"),
            )
        })?;
        if let Some(row) = rows.selected.iter().find(|row| row.point.is_none()) {
            return Err(StructureInterpretError::new(
                DiagnosticKind::InvalidValue,
                Some(row.line),
                format!(
                    "selected coordinate model `{model_id}` has no complete position for atom `{}`",
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
            report,
        )
    }
}

fn publish_model(
    connections: &[ConnectionRecord],
    rows: ModelRows,
    polymer_asym_order: &BTreeMap<String, usize>,
    catalog: &ConnectivityCatalog,
    mut report: StructureInterpretationReport,
) -> Result<StructureInterpretation, StructureInterpretError> {
    let ModelRows {
        all: rows,
        selected,
    } = rows;
    let selected_model = report
        .selected_model
        .as_deref()
        .expect("selected model is recorded")
        .to_owned();
    let mut union = InstanceUnion::new(selected.iter().map(|row| row.instance_key.clone()));
    let connections = resolve_connections(
        connections,
        &selected,
        &rows,
        &selected_model,
        &mut union,
        &mut report,
    );
    let groups = group_rows(selected, &mut union, polymer_asym_order);
    let mut group_connections = BTreeMap::<String, Vec<&DeclaredConnection>>::new();
    for connection in &connections {
        group_connections
            .entry(union.find(&connection.instance_key))
            .or_default()
            .push(connection);
    }
    let mut builder = ModelBuilder::new();
    let mut qualified_atom_data = Vec::new();
    for group in groups {
        let declared = group_connections.remove(group.root()).unwrap_or_default();
        let mut built = build_molecule(group, &declared)?;
        built.complete_connectivity(catalog, &mut report)?;
        for built in built.publish_components()? {
            let id = builder
                .add_molecule(built.molecule, &built.positions)
                .map_err(graph_error)?;
            let (provenance, atom_data) = built.provenance.qualify(id);
            report.instances.push(provenance);
            qualified_atom_data.extend(atom_data);
        }
    }
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
    *builder.topology_builder_mut().hierarchy_mut() =
        build_topology_hierarchy(&report.instances, polymer_asym_order)?;
    let mut model = builder.build().map_err(graph_error)?;
    let mut occupancies = vec![None; model.atom_count()];
    let mut b_factors = vec![None; model.atom_count()];
    let b_factor_scale = SQUARE_ANGSTROM
        .conversion_factor_to(SQUARE_NANOMETER)
        .map_err(graph_error)?;
    for (atom, occupancy, b_factor) in qualified_atom_data {
        let index = model
            .topology()
            .atom_index(atom)
            .ok_or_else(|| graph_error(format_args!("invalid qualified mmCIF atom: {atom}")))?;
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
