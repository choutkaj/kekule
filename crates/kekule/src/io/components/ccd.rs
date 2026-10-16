//! Reads Chemical Component Dictionary entries from CIF data blocks.

use std::collections::BTreeSet;

use crate::diagnostics::{DiagnosticKind, SourceLocation};

use super::super::element_symbol::parse_structure_element_symbol;
use super::super::mmcif_category::MmcifCategory;
use super::super::{MmcifBlock, MmcifValue};
use super::{Component, ComponentAtom, ComponentBond, ComponentBondOrder, ComponentError};

pub(super) fn component_from_block(block: &MmcifBlock) -> Result<Component, ComponentError> {
    let scalar = |tag: &str| block.item(tag).and_then(|value| value.optional_text());
    let id = scalar("_chem_comp.id").unwrap_or(block.name()).to_owned();
    let error = |kind, value: Option<MmcifValue<'_>>, message: String| {
        ComponentError::new(
            kind,
            value.map(|value| value.location()),
            format!("component `{id}`: {message}"),
        )
    };

    let atoms_table = category(block, "_chem_comp_atom", &id)?.ok_or_else(|| {
        ComponentError::new(
            DiagnosticKind::Syntax,
            None,
            format!("component `{id}` has no _chem_comp_atom category"),
        )
    })?;
    let mut atoms = Vec::with_capacity(atoms_table.row_count());
    let mut names = BTreeSet::new();
    for row in 0..atoms_table.row_count() {
        let field = |tag: &str| atoms_table.value(row, tag);
        let text = |tag: &str| field(tag).and_then(|value| value.optional_text());
        let required = |tag: &str| {
            text(tag).ok_or_else(|| {
                ComponentError::new(
                    DiagnosticKind::Syntax,
                    atoms_table.row_line(row).map(SourceLocation::at_line),
                    format!("component `{id}`: atom row lacks {tag}"),
                )
            })
        };
        check_component(&id, field("_chem_comp_atom.comp_id"))?;
        let name = required("_chem_comp_atom.atom_id")?;
        if !names.insert(name.to_owned()) {
            return Err(error(
                DiagnosticKind::Inconsistent,
                field("_chem_comp_atom.atom_id"),
                format!("duplicate atom `{name}`"),
            ));
        }
        let symbol = required("_chem_comp_atom.type_symbol")?;
        let (element, isotope) = parse_structure_element_symbol(symbol).ok_or_else(|| {
            error(
                DiagnosticKind::InvalidValue,
                field("_chem_comp_atom.type_symbol"),
                format!("unknown element `{symbol}`"),
            )
        })?;
        let formal_charge = text("_chem_comp_atom.charge")
            .map(|charge| {
                charge.parse::<i8>().map_err(|_| {
                    error(
                        DiagnosticKind::InvalidValue,
                        field("_chem_comp_atom.charge"),
                        format!("invalid charge `{charge}`"),
                    )
                })
            })
            .transpose()?;
        let alternate_name = text("_chem_comp_atom.alt_atom_id")
            .filter(|alternate| *alternate != name)
            .map(str::to_owned);
        atoms.push(ComponentAtom {
            name: name.to_owned(),
            alternate_name,
            element,
            isotope,
            formal_charge,
            aromatic: flag(text("_chem_comp_atom.pdbx_aromatic_flag")),
            leaving: flag(text("_chem_comp_atom.pdbx_leaving_atom_flag")),
        });
    }

    let mut bonds = Vec::new();
    let mut pairs = BTreeSet::new();
    if let Some(table) = category(block, "_chem_comp_bond", &id)? {
        for row in 0..table.row_count() {
            let field = |tag: &str| table.value(row, tag);
            let text = |tag: &str| field(tag).and_then(|value| value.optional_text());
            check_component(&id, field("_chem_comp_bond.comp_id"))?;
            let mut endpoints = [0usize; 2];
            for (endpoint, tag) in endpoints
                .iter_mut()
                .zip(["_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2"])
            {
                let name = text(tag).ok_or_else(|| {
                    ComponentError::new(
                        DiagnosticKind::Syntax,
                        table.row_line(row).map(SourceLocation::at_line),
                        format!("component `{id}`: bond row lacks {tag}"),
                    )
                })?;
                *endpoint = atoms
                    .iter()
                    .position(|atom| atom.name == name)
                    .ok_or_else(|| {
                        error(
                            DiagnosticKind::InvalidReference,
                            field(tag),
                            format!("bond names unknown atom `{name}`"),
                        )
                    })?;
            }
            let [a, b] = endpoints;
            if a == b || !pairs.insert((a.min(b), a.max(b))) {
                return Err(error(
                    DiagnosticKind::Inconsistent,
                    field("_chem_comp_bond.atom_id_1"),
                    format!(
                        "self or duplicate bond `{}`-`{}`",
                        atoms[a].name, atoms[b].name
                    ),
                ));
            }
            let order = match text("_chem_comp_bond.value_order")
                .map(str::to_ascii_uppercase)
                .as_deref()
            {
                None | Some("SING") => ComponentBondOrder::Single,
                Some("DOUB") => ComponentBondOrder::Double,
                Some("TRIP") => ComponentBondOrder::Triple,
                Some("QUAD") => ComponentBondOrder::Quadruple,
                Some("AROM") => ComponentBondOrder::Aromatic,
                Some("DELO") => ComponentBondOrder::Delocalized,
                Some("PI") => ComponentBondOrder::Pi,
                Some("POLY") => ComponentBondOrder::Polymeric,
                Some(other) => {
                    return Err(error(
                        DiagnosticKind::InvalidValue,
                        field("_chem_comp_bond.value_order"),
                        format!("unknown bond order `{other}`"),
                    ))
                }
            };
            bonds.push(ComponentBond {
                atoms: endpoints,
                order,
                aromatic: flag(text("_chem_comp_bond.pdbx_aromatic_flag")),
                metal_coordination: flag(text("_chem_comp_bond.pdbx_metal_coordination_flag")),
            });
        }
    }

    Ok(Component {
        name: scalar("_chem_comp.name").unwrap_or_default().to_owned(),
        kind: scalar("_chem_comp.type").unwrap_or_default().to_owned(),
        id,
        variant: None,
        atoms,
        bonds,
    })
}

fn category<'a>(
    block: &'a MmcifBlock,
    name: &str,
    id: &str,
) -> Result<Option<MmcifCategory<'a>>, ComponentError> {
    block.category(name).map_err(|error| {
        ComponentError::new(
            error.kind(),
            error.location(),
            format!("component `{id}`: {}", error.message()),
        )
    })
}

/// Rows of a multi-component block must belong to the block's component.
fn check_component(id: &str, value: Option<MmcifValue<'_>>) -> Result<(), ComponentError> {
    match value.and_then(|value| value.optional_text().map(|text| (value, text))) {
        Some((value, text)) if !text.eq_ignore_ascii_case(id) => Err(ComponentError::new(
            DiagnosticKind::Inconsistent,
            Some(value.location()),
            format!("component `{id}` contains a row of component `{text}`"),
        )),
        _ => Ok(()),
    }
}

fn flag(value: Option<&str>) -> bool {
    value.is_some_and(|value| value.eq_ignore_ascii_case("Y"))
}
