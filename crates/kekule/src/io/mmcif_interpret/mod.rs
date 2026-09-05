mod atom_site;
mod build;
mod ensemble;
mod struct_conn;
mod types;

use crate::structure::ModelBuilder;
use crate::units::{SQUARE_ANGSTROM, SQUARE_NANOMETER};

use super::{MmcifBlock, MmcifDocument};
use atom_site::{read_asym_entities, read_atom_rows, read_entity_types, select_alt_locations};
use build::{
    build_molecule, build_topology_hierarchy, graph_error, group_rows, polymer_asym_order,
};
use struct_conn::{read_connections, InstanceUnion};

pub(crate) use ensemble::{interpret_mmcif_ensemble, interpret_mmcif_ensemble_block};
pub use types::{
    MmcifAltLocPolicy, MmcifAtomProvenance, MmcifConnectionResolutionReason,
    MmcifEnsembleInterpretError, MmcifEnsembleInterpretOptions, MmcifEnsembleInterpretation,
    MmcifEntityKind, MmcifInstanceProvenance, MmcifInterpretError, MmcifInterpretIssue,
    MmcifInterpretOptions, MmcifInterpretation, MmcifInterpretationReport, MmcifModelSelection,
};

pub(crate) fn interpret_mmcif(
    document: &MmcifDocument,
    options: MmcifInterpretOptions,
) -> Result<MmcifInterpretation, MmcifInterpretError> {
    let blocks = document
        .blocks()
        .iter()
        .filter(|block| block.loop_with_tag("_atom_site.type_symbol").is_some())
        .collect::<Vec<_>>();
    if blocks.is_empty() {
        return Err(MmcifInterpretError::new(
            None,
            "document has no atom-site loop",
        ));
    }
    if blocks.len() > 1 {
        return Err(MmcifInterpretError::new(
            None,
            "document has atom-site data in more than one data block",
        ));
    }
    interpret_mmcif_block(blocks[0], options)
}

pub(crate) fn interpret_mmcif_block(
    block: &MmcifBlock,
    options: MmcifInterpretOptions,
) -> Result<MmcifInterpretation, MmcifInterpretError> {
    let mut prepared = PreparedBlock::new(
        block,
        options.strict_entity_metadata,
        &options.altloc_policy,
    )?;
    let model_id = match &options.model_selection {
        MmcifModelSelection::RequireSingle if prepared.model_sizes.len() == 1 => {
            prepared.model_sizes[0].0.clone()
        }
        MmcifModelSelection::RequireSingle => {
            return Err(MmcifInterpretError::new(
                None,
                format!(
                    "coordinate data contains {} models; select one explicitly",
                    prepared.model_sizes.len(),
                ),
            ))
        }
        MmcifModelSelection::First => prepared.model_sizes[0].0.clone(),
        MmcifModelSelection::Select(id) => id.clone(),
    };
    prepared.interpret_model(&model_id)
}

#[derive(Default)]
struct ModelRows {
    all: Vec<atom_site::AtomRow>,
    selected: Vec<atom_site::AtomRow>,
}

/// Shared source interpretation for a block; each model's rows are consumed once.
struct PreparedBlock<'a> {
    block: &'a MmcifBlock,
    models: std::collections::BTreeMap<String, ModelRows>,
    model_sizes: Vec<(String, usize)>,
    report: MmcifInterpretationReport,
    polymer_order: std::collections::BTreeMap<String, usize>,
    connectivity: super::mmcif_connectivity::ConnectivityCatalog,
}

impl<'a> PreparedBlock<'a> {
    fn new(
        block: &'a MmcifBlock,
        strict_entity_metadata: bool,
        altloc_policy: &MmcifAltLocPolicy,
    ) -> Result<Self, MmcifInterpretError> {
        let entities = read_entity_types(block)?;
        let asym_entities = read_asym_entities(block)?;
        let atom_table = block
            .loop_with_tag("_atom_site.type_symbol")
            .ok_or_else(|| MmcifInterpretError::new(None, "block has no atom-site loop"))?;
        if atom_table.row_count() == 0 {
            return Err(MmcifInterpretError::new(
                None,
                "atom-site loop contains no rows",
            ));
        }
        let mut report = MmcifInterpretationReport {
            block_name: block.name().to_owned(),
            entity_definitions: entities.len(),
            ..MmcifInterpretationReport::default()
        };
        let rows = read_atom_rows(
            atom_table,
            &entities,
            &asym_entities,
            strict_entity_metadata,
            &mut report,
        )?;
        // Alternate-location validation still covers the complete source block,
        // including models omitted from the eventual projection.
        let selected = select_alt_locations(&rows, altloc_policy, &mut report)?;
        let mut models = std::collections::BTreeMap::<String, ModelRows>::new();
        for row in rows {
            let model = models.entry(row.model_id.clone()).or_default();
            model.all.push(row);
        }
        let mut model_ids = Vec::new();
        for row in selected {
            let model = models
                .get_mut(&row.model_id)
                .expect("selected row has a source model");
            if model.selected.is_empty() {
                model_ids.push(row.model_id.clone());
            }
            model.selected.push(row);
        }
        let model_sizes = model_ids
            .into_iter()
            .map(|id| {
                let count = models[&id].selected.len();
                (id, count)
            })
            .collect::<Vec<_>>();
        report.coordinate_models = model_sizes.len();
        Ok(Self {
            block,
            models,
            model_sizes,
            report,
            polymer_order: polymer_asym_order(block),
            connectivity: super::mmcif_connectivity::ConnectivityCatalog::from_block(block)?,
        })
    }

    fn interpret_model(
        &mut self,
        model_id: &str,
    ) -> Result<MmcifInterpretation, MmcifInterpretError> {
        let rows = self.models.remove(model_id).ok_or_else(|| {
            MmcifInterpretError::new(
                None,
                format!("coordinate model `{model_id}` is unavailable"),
            )
        })?;
        if let Some(row) = rows.selected.iter().find(|row| row.point.is_none()) {
            return Err(MmcifInterpretError::new(
                Some(row.line),
                format!(
                    "selected coordinate model `{model_id}` has no complete position for atom `{}`",
                    row.atom_name,
                ),
            ));
        }
        let mut report = self.report.clone();
        report.selected_model = Some(model_id.to_owned());
        for (ignored, count) in &self.model_sizes {
            if ignored != model_id {
                report.ignored_coordinate_models.push(ignored.clone());
                report
                    .issues
                    .push(MmcifInterpretIssue::CoordinateModelIgnored {
                        model_id: ignored.clone(),
                        atom_site_rows: *count,
                    });
            }
        }
        publish_model(
            self.block,
            rows,
            &self.polymer_order,
            &self.connectivity,
            report,
        )
    }
}

fn publish_model(
    block: &MmcifBlock,
    rows: ModelRows,
    polymer_asym_order: &std::collections::BTreeMap<String, usize>,
    connectivity: &super::mmcif_connectivity::ConnectivityCatalog,
    mut report: MmcifInterpretationReport,
) -> Result<MmcifInterpretation, MmcifInterpretError> {
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
    let connections = read_connections(
        block,
        &selected,
        &rows,
        &selected_model,
        &mut union,
        &mut report,
    )?;
    let groups = group_rows(selected, &mut union, polymer_asym_order);
    let mut builder = ModelBuilder::new();
    let mut qualified_atom_data = Vec::new();
    for group in groups {
        let mut built = build_molecule(group, &connections, &mut report)?;
        built.complete_connectivity(connectivity)?;
        for built in built.publish_components()? {
            let id = builder
                .add_molecule(&built.molecule, &built.positions)
                .map_err(graph_error)?;
            let (provenance, atom_data) = built.provenance.qualify(id);
            report.instances.push(provenance);
            qualified_atom_data.extend(atom_data);
        }
    }
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
        .install_canonical_atom_properties(occupancies, b_factors)
        .map_err(graph_error)?;
    report.macromolecules = report
        .instances
        .iter()
        .filter(|instance| {
            instance
                .entity_kinds()
                .iter()
                .any(MmcifEntityKind::is_macro)
        })
        .count();
    report.small_molecules = report.instances.len() - report.macromolecules;
    report.solvent_molecules = report
        .instances
        .iter()
        .filter(|instance| instance.entity_kinds().contains(&MmcifEntityKind::Water))
        .count();
    Ok(MmcifInterpretation { model, report })
}
