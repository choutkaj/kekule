//! Source-declared combinations constrain selection; they do not create members.
use std::collections::{BTreeMap, BTreeSet};

use super::super::atom_site::{required, row_error, AtomRow};
use super::super::MmcifInterpretError;
use crate::io::MmcifBlock;

pub(super) fn conformations(
    block: &MmcifBlock,
    rows: &[AtomRow],
) -> Result<BTreeMap<String, BTreeSet<String>>, MmcifInterpretError> {
    let mut conformations = BTreeMap::<String, BTreeSet<String>>::new();
    let Some(table) = block.category("_atom_sites_alt_gen")? else {
        return Ok(conformations);
    };
    let available = rows
        .iter()
        .filter_map(|row| row.alt_id.as_deref())
        .collect::<BTreeSet<_>>();
    for row in 0..table.row_count() {
        let id = required(&table, row, "_atom_sites_alt_gen.ens_id")?;
        let value = table
            .value(row, "_atom_sites_alt_gen.alt_id")
            .ok_or_else(|| row_error(&table, row, "missing _atom_sites_alt_gen.alt_id"))?;
        // The dictionary examples include '.' for shared sites. Shared rows are
        // retained for every combination regardless of whether '.' is listed.
        let label = value.text();
        if label == "?" || (label != "." && !available.contains(label)) {
            return Err(row_error(
                &table,
                row,
                format!(
                    "source conformation `{id}` references unavailable alternate label `{label}`"
                ),
            ));
        }
        if !conformations
            .entry(id.to_owned())
            .or_default()
            .insert(label.to_owned())
        {
            return Err(row_error(
                &table,
                row,
                format!("duplicate alternate label `{label}` in source conformation `{id}`"),
            ));
        }
    }
    Ok(conformations)
}
