mod atom_site;
mod build;
mod ensemble;
mod struct_conn;
mod types;

use crate::structure::{AtomData, ModelBuilder};
use crate::topology::MoleculeInstanceId;
use crate::units::{Quantity, SQUARE_ANGSTROM};

use super::{MmcifDataBlock, MmcifDocument};
use atom_site::{
    read_asym_entities, read_atom_rows, read_entity_types, select_alt_locations,
    select_coordinate_model,
};
use build::{build_molecule, graph_error, group_rows, polymer_asym_order};
use struct_conn::{read_connections, InstanceUnion};

pub(crate) use ensemble::interpret_mmcif_ensemble;
pub(crate) use types::MmcifInterpretation;
pub use types::{
    MmcifAltLocPolicy, MmcifAtomProvenance, MmcifConnectionResolutionReason,
    MmcifEnsembleInterpretError, MmcifEnsembleInterpretOptions, MmcifEntityKind,
    MmcifInstanceProvenance, MmcifInterpretError, MmcifInterpretIssue, MmcifInterpretOptions,
    MmcifInterpretationReport, MmcifModelSelection,
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
    interpret_block(blocks[0], options)
}

fn interpret_block(
    block: &MmcifDataBlock,
    options: MmcifInterpretOptions,
) -> Result<MmcifInterpretation, MmcifInterpretError> {
    let entities = read_entity_types(block)?;
    let asym_entities = read_asym_entities(block)?;
    let atom_table = block
        .loop_with_tag("_atom_site.type_symbol")
        .expect("selected block has atom-site data");
    if atom_table.row_count() == 0 {
        return Err(MmcifInterpretError::new(
            None,
            "atom-site loop contains no rows",
        ));
    }
    let mut report = MmcifInterpretationReport {
        data_block: block.name().to_owned(),
        entity_definitions: entities.len(),
        ..MmcifInterpretationReport::default()
    };
    let rows = read_atom_rows(atom_table, &entities, &asym_entities, &options, &mut report)?;
    let selected = select_alt_locations(&rows, &options.altloc_policy, &mut report)?;
    let selected = select_coordinate_model(selected, &options.model_selection, &mut report)?;
    let selected_model = report
        .selected_model
        .clone()
        .ok_or_else(|| MmcifInterpretError::new(None, "coordinate model selection was lost"))?;
    let mut union = InstanceUnion::new(selected.iter().map(|row| row.instance_key.clone()));
    let connections = read_connections(
        block,
        &selected,
        &rows,
        &selected_model,
        &mut union,
        &mut report,
    )?;
    let polymer_asym_order = polymer_asym_order(block);
    let groups = group_rows(selected, &mut union, &polymer_asym_order);
    let mut builder = ModelBuilder::new();
    let mut qualified_atom_data = Vec::new();
    for group in groups {
        let mut built = build_molecule(group, &connections, &mut report)?;
        let (staged_provenance, _) = built.provenance.clone().qualify(MoleculeInstanceId::new(0));
        super::mmcif_connectivity::complete_editor_connectivity(
            block,
            &mut built.editor,
            &staged_provenance,
        )?;
        for built in built.publish_components()? {
            let id = builder
                .add_molecule(&built.molecule, &built.conformer)
                .map_err(graph_error)?;
            let (provenance, atom_data) = built.provenance.qualify(id);
            report.instances.push(provenance);
            qualified_atom_data.extend(atom_data);
        }
    }
    let mut model = builder.build().map_err(graph_error)?;
    let mut occupancies = vec![None; model.topology().atom_count()];
    let mut b_factors = vec![None; model.topology().atom_count()];
    for (atom, occupancy, b_factor) in qualified_atom_data {
        let index = model
            .topology()
            .atom_index(atom)
            .expect("interpreted atom has a dense topology index");
        occupancies[index.index()] = occupancy;
        b_factors[index.index()] = b_factor;
    }
    let topology = model.shared_topology();
    let mut atom_data = AtomData::new(&topology);
    atom_data
        .set_occupancies(occupancies)
        .map_err(graph_error)?;
    atom_data
        .set_b_factors(Quantity::new(b_factors, SQUARE_ANGSTROM))
        .map_err(graph_error)?;
    model.set_atom_data(atom_data).map_err(graph_error)?;
    report.macromolecules = model
        .topology()
        .instances()
        .filter(|(id, _)| {
            model
                .topology()
                .definition_for_instance(*id)
                .is_ok_and(|definition| definition.hierarchy().is_some())
        })
        .count();
    report.small_molecules = model
        .topology()
        .instances()
        .filter(|(id, _)| {
            model
                .topology()
                .definition_for_instance(*id)
                .is_ok_and(|definition| definition.hierarchy().is_none())
        })
        .count();
    report.solvent_molecules = report
        .instances
        .iter()
        .filter(|instance| instance.entity_kinds().contains(&MmcifEntityKind::Water))
        .count();
    Ok(MmcifInterpretation { model, report })
}
