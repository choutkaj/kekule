//! Reads an mmCIF data block into the source-neutral interpretation input.
use std::collections::BTreeMap;

use crate::core::BondOrder;
use crate::diagnostics::DiagnosticKind;

use super::super::mmcif_category::MmcifCategory;
use super::super::structure_connectivity::ConnectivityCatalog;
use super::super::{MmcifBlock, MmcifDocument};
use super::altloc::{inventory, source::conformations};
use super::atom_site::{
    coordinate_model_ids, optional, optional_i32, read_asym_entities, read_atom_rows,
    read_entity_types, required, row_error,
};
use super::connections::{
    AuthorConnectionSelector, ConnectionPartnerSelector, ConnectionRecord, DeclaredBond,
    DeclaredOrder, LabelConnectionSelector, PartnerSpec,
};
use super::ensemble::{interpret_conformations, interpret_ensemble};
use super::types::{
    ConnectionResolutionReason, ConnectivityOptions, StructureEnsembleInterpretError,
    StructureEnsembleInterpretOptions, StructureEnsembleInterpretation, StructureInterpretError,
    StructureInterpretOptions, StructureInterpretation, StructureInterpretationReport,
};
use super::{interpret_source, AltLocResidue, StructureSource};

pub(crate) fn interpret_mmcif(
    document: &MmcifDocument,
    options: StructureInterpretOptions,
) -> Result<StructureInterpretation, StructureInterpretError> {
    let blocks = document
        .blocks()
        .iter()
        .filter(|block| block.has_category("_atom_site"))
        .collect::<Vec<_>>();
    if blocks.is_empty() {
        return Err(StructureInterpretError::new(
            DiagnosticKind::Inconsistent,
            None,
            "document has no atom-site category",
        ));
    }
    if blocks.len() > 1 {
        return Err(StructureInterpretError::new(
            DiagnosticKind::Inconsistent,
            None,
            "document has atom-site data in more than one data block",
        ));
    }
    interpret_mmcif_block(blocks[0], options)
}

pub(crate) fn interpret_mmcif_block(
    block: &MmcifBlock,
    options: StructureInterpretOptions,
) -> Result<StructureInterpretation, StructureInterpretError> {
    let source = mmcif_source(block, options.strict_entity_metadata, &options.connectivity)?;
    interpret_source(source, &options)
}

/// Interprets explicitly selected or all coordinate models as one
/// shared-topology non-temporal ensemble.
pub(crate) fn interpret_mmcif_ensemble(
    document: &MmcifDocument,
    options: StructureEnsembleInterpretOptions,
) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
    interpret_mmcif_ensemble_block(atom_site_block(document)?, options)
}

/// Interprets explicitly selected or all coordinate models in one block as
/// one shared-topology non-temporal ensemble.
pub(crate) fn interpret_mmcif_ensemble_block(
    block: &MmcifBlock,
    options: StructureEnsembleInterpretOptions,
) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
    if options.model_ids.as_ref().is_some_and(Vec::is_empty) {
        return Err(StructureEnsembleInterpretError::EmptyModelSelection);
    }
    if !block.has_category("_atom_site") {
        return Err(StructureEnsembleInterpretError::NoCoordinateModels);
    }
    let available =
        coordinate_model_ids(block).map_err(|error| StructureEnsembleInterpretError::Model {
            model_id: "<model inventory>".to_owned(),
            error,
        })?;
    interpret_ensemble(available, options, |options| {
        mmcif_source(block, options.strict_entity_metadata, &options.connectivity)
    })
}

/// Explicit caller-selected alternate conformations, never automatic altloc
/// enumeration.
pub(crate) fn interpret_mmcif_conformations(
    document: &MmcifDocument,
    selections: &[StructureInterpretOptions],
) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
    interpret_mmcif_conformations_block(atom_site_block(document)?, selections)
}

pub(crate) fn interpret_mmcif_conformations_block(
    block: &MmcifBlock,
    selections: &[StructureInterpretOptions],
) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
    interpret_conformations(selections, |options| {
        interpret_mmcif_block(block, options.clone())
    })
}

pub(crate) fn alternate_location_inventory(
    block: &MmcifBlock,
) -> Result<Vec<AltLocResidue>, StructureInterpretError> {
    let entities = read_entity_types(block)?;
    let asym = read_asym_entities(block)?;
    let table = block.category("_atom_site")?.ok_or_else(|| {
        StructureInterpretError::new(
            DiagnosticKind::Inconsistent,
            None,
            "block has no atom-site category",
        )
    })?;
    let rows = read_atom_rows(
        &table,
        &entities,
        &asym,
        false,
        &mut StructureInterpretationReport::default(),
    )?;
    conformations(block, &rows)?;
    inventory(&rows)
}

fn atom_site_block(
    document: &MmcifDocument,
) -> Result<&MmcifBlock, StructureEnsembleInterpretError> {
    let mut blocks = document
        .blocks()
        .iter()
        .filter(|block| block.has_category("_atom_site"));
    let block = blocks
        .next()
        .ok_or(StructureEnsembleInterpretError::NoCoordinateModels)?;
    if blocks.next().is_some() {
        return Err(StructureEnsembleInterpretError::MultipleAtomSiteBlocks);
    }
    Ok(block)
}

fn mmcif_source(
    block: &MmcifBlock,
    strict_entity_metadata: bool,
    connectivity: &ConnectivityOptions,
) -> Result<StructureSource, StructureInterpretError> {
    let entities = read_entity_types(block)?;
    let asym_entities = read_asym_entities(block)?;
    let atom_table = block.category("_atom_site")?.ok_or_else(|| {
        StructureInterpretError::new(
            DiagnosticKind::Inconsistent,
            None,
            "block has no atom-site category",
        )
    })?;
    if atom_table.row_count() == 0 {
        return Err(StructureInterpretError::new(
            DiagnosticKind::Inconsistent,
            None,
            "atom-site category contains no rows",
        ));
    }
    let mut report = StructureInterpretationReport::default();
    let rows = read_atom_rows(
        &atom_table,
        &entities,
        &asym_entities,
        strict_entity_metadata,
        &mut report,
    )?;
    Ok(StructureSource {
        name: block.name().to_owned(),
        entity_definitions: entities.len(),
        conformations: conformations(block, &rows)?,
        rows,
        polymer_order: polymer_asym_order(block)?,
        connections: read_struct_conn(block)?,
        catalog: ConnectivityCatalog::from_block(block, connectivity)?,
        issues: report.issues,
    })
}

pub(super) fn polymer_asym_order(
    block: &MmcifBlock,
) -> Result<BTreeMap<String, usize>, StructureInterpretError> {
    let mut order = BTreeMap::new();
    let Some(table) = block.category("_pdbx_poly_seq_scheme")? else {
        return Ok(order);
    };
    for row in 0..table.row_count() {
        let asym_id = required(&table, row, "_pdbx_poly_seq_scheme.asym_id")?;
        let next = order.len();
        order.entry(asym_id.to_owned()).or_insert(next);
    }
    Ok(order)
}

fn read_struct_conn(block: &MmcifBlock) -> Result<Vec<ConnectionRecord>, StructureInterpretError> {
    let Some(table) = block.category("_struct_conn")? else {
        return Ok(Vec::new());
    };
    let table = &table;
    let mut records = Vec::new();
    for row in 0..table.row_count() {
        let kind = required(table, row, "_struct_conn.conn_type_id")?.to_owned();
        if !is_covalent_connection(&kind) {
            records.push(ConnectionRecord::Ignored {
                connection_type: kind,
            });
            continue;
        }
        // Metal coordination is covalent; the dative order keeps it from
        // counting toward the valence of the ligand atom.
        let order = if is_metal_coordination(&kind) {
            BondOrder::Dative
        } else {
            connection_bond_order(table, row)?
        };
        let partners = [
            connection_partner(table, row, 1)?,
            connection_partner(table, row, 2)?,
        ];
        records.push(ConnectionRecord::Bond(Box::new(DeclaredBond {
            connection_id: optional(table, row, "_struct_conn.id").map(str::to_owned),
            connection_type: kind,
            line: table.row_line(row),
            order: DeclaredOrder::Stated(order),
            partners,
            connectivity_only: false,
        })));
    }
    Ok(records)
}

fn connection_bond_order(
    table: &MmcifCategory<'_>,
    row: usize,
) -> Result<BondOrder, StructureInterpretError> {
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
            DiagnosticKind::Unsupported,
            format!("unsupported struct_conn bond order `{order}`"),
        )),
    }
}

fn is_covalent_connection(kind: &str) -> bool {
    let kind = kind.to_ascii_lowercase();
    kind.starts_with("covale") || kind == "disulf" || kind == "modres" || kind == "metalc"
}

fn is_metal_coordination(kind: &str) -> bool {
    kind.eq_ignore_ascii_case("metalc")
}

fn connection_partner(
    table: &MmcifCategory<'_>,
    row: usize,
    partner: u8,
) -> Result<PartnerSpec, StructureInterpretError> {
    let symmetry_tag = format!("_struct_conn.ptnr{partner}_symmetry");
    if let Some(symmetry) = optional(table, row, &symmetry_tag) {
        if symmetry != "1_555" {
            return Ok(Err(ConnectionResolutionReason::UnsupportedSymmetry {
                symmetry: symmetry.to_owned(),
            }));
        }
    }

    let label_alt_tag = format!("_struct_conn.ptnr{partner}_label_alt_id");
    let pdbx_label_alt_tag = format!("_struct_conn.pdbx_ptnr{partner}_label_alt_id");
    let label_alt = optional(table, row, &label_alt_tag);
    let pdbx_label_alt = optional(table, row, &pdbx_label_alt_tag);
    let (alternate_location, mut conflict) = match (label_alt, pdbx_label_alt) {
        (Some(left), Some(right)) if left != right => (
            Some(left.to_owned()),
            Some(ConnectionResolutionReason::ConflictingSelectorValues {
                selector: "label alternate-location aliases",
            }),
        ),
        (Some(value), _) | (_, Some(value)) => (Some(value.to_owned()), None),
        (None, None) => (None, None),
    };
    let text =
        |field: &str| optional(table, row, &format!("_struct_conn.{field}")).map(str::to_owned);

    let selector = ConnectionPartnerSelector {
        label: LabelConnectionSelector {
            asym_id: text(&format!("ptnr{partner}_label_asym_id")),
            component_id: text(&format!("ptnr{partner}_label_comp_id")),
            sequence_id: optional_i32(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_label_seq_id"),
            )?,
            atom_id: text(&format!("ptnr{partner}_label_atom_id")),
            alternate_location,
        },
        author: AuthorConnectionSelector {
            asym_id: text(&format!("ptnr{partner}_auth_asym_id")),
            component_id: text(&format!("ptnr{partner}_auth_comp_id")),
            sequence_id: text(&format!("ptnr{partner}_auth_seq_id")),
            atom_id: text(&format!("ptnr{partner}_auth_atom_id")),
            insertion_code: text(&format!("pdbx_ptnr{partner}_PDB_ins_code")),
            alternate_location: text(&format!("pdbx_ptnr{partner}_auth_alt_id")),
        },
        site_id: None,
    };
    if let (Some(label), Some(author)) = (
        selector.label.alternate_location.as_deref(),
        selector.author.alternate_location.as_deref(),
    ) {
        if label != author {
            conflict = Some(ConnectionResolutionReason::ConflictingSelectorValues {
                selector: "label and author alternate locations",
            });
        }
    }
    if let Some(reason) = conflict {
        return Ok(Err(reason));
    }
    if selector.is_empty() {
        return Ok(Err(ConnectionResolutionReason::MissingSelector));
    }
    Ok(Ok(selector))
}

#[cfg(test)]
mod syntax_tests {
    #[test]
    fn singleton_polymer_scheme_retains_declared_chain_order() {
        let looped =
            "data_x\nloop_\n_pdbx_poly_seq_scheme.asym_id\n_pdbx_poly_seq_scheme.seq_id\nZ 1\n";
        for input in [
            looped.to_owned(),
            crate::tests::mmcif_syntax::singleton_loops_as_scalars(looped),
        ] {
            let doc = crate::mmcif::parse_str(&input).unwrap();
            let order = super::polymer_asym_order(&doc.blocks()[0]).unwrap();
            assert_eq!(order.len(), 1);
            assert_eq!(order["Z"], 0);
        }
    }
}
