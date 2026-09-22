//! Emit prepared rows with format-specific quoting and numeric formatting.
use super::prepare::{AtomRow, PreparedModel};
use super::{MmcifWriteError, MmcifWriteOptions};
use crate::core::BondOrder;
use std::collections::BTreeMap;
use std::io::Write;

fn mmcif_io(error: std::io::Error) -> MmcifWriteError {
    MmcifWriteError::Io {
        kind: error.kind(),
        message: error.to_string(),
    }
}

pub(super) fn render_model_to(
    writer: &mut impl Write,
    model: &PreparedModel,
    options: &MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    write_block_start(writer, model, options)?;
    let mut atom_serial = 1u64;
    write_atom_rows(writer, &model.atoms, 1, &mut atom_serial, options)?;
    write_block_end(writer, model)
}

pub(super) fn write_block_start(
    writer: &mut impl Write,
    model: &PreparedModel,
    options: &MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    writer
        .write_all(format!("data_{}\n#\n", options.block_name).as_bytes())
        .map_err(mmcif_io)?;

    write_loop_header(writer, &["_entity.id", "_entity.type"])?;
    for entity in &model.entities {
        write_row(
            writer,
            vec![
                cif_value(&entity.id, "_entity.id")?,
                entity.kind.as_mmcif().to_owned(),
            ],
        )?;
    }
    writer.write_all(b"#\n").map_err(mmcif_io)?;

    write_loop_header(writer, &["_struct_asym.id", "_struct_asym.entity_id"])?;
    for asym in &model.asyms {
        write_row(
            writer,
            vec![
                cif_value(&asym.id, "_struct_asym.id")?,
                cif_value(&asym.entity_id, "_struct_asym.entity_id")?,
            ],
        )?;
    }
    writer.write_all(b"#\n").map_err(mmcif_io)?;

    const ATOM_TAGS: &[&str] = &[
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_formal_charge",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    ];
    write_loop_header(writer, ATOM_TAGS)
}

pub(super) fn write_atom_rows(
    writer: &mut impl Write,
    atoms: &[AtomRow],
    model_number: usize,
    serial: &mut u64,
    options: &MmcifWriteOptions,
) -> Result<(), MmcifWriteError> {
    for atom in atoms {
        write_row(
            writer,
            vec![
                atom.group_pdb.clone(),
                serial.to_string(),
                cif_value(&atom.type_symbol, "_atom_site.type_symbol")?,
                cif_value(&atom.label_atom_id, "_atom_site.label_atom_id")?,
                optional_cif_value(atom.label_alt_id.as_deref(), "_atom_site.label_alt_id")?,
                cif_value(&atom.label_comp_id, "_atom_site.label_comp_id")?,
                cif_value(&atom.asym_id, "_atom_site.label_asym_id")?,
                cif_value(&atom.entity_id, "_atom_site.label_entity_id")?,
                optional_display(atom.label_seq_id),
                optional_cif_value(
                    atom.insertion_code.as_deref(),
                    "_atom_site.pdbx_PDB_ins_code",
                )?,
                format_coordinate(atom.position.x, options.coordinate_precision),
                format_coordinate(atom.position.y, options.coordinate_precision),
                format_coordinate(atom.position.z, options.coordinate_precision),
                optional_float(atom.occupancy),
                optional_float(atom.b_factor),
                atom.formal_charge.to_string(),
                optional_cif_value(atom.auth_seq_id.as_deref(), "_atom_site.auth_seq_id")?,
                cif_value(&atom.auth_comp_id, "_atom_site.auth_comp_id")?,
                cif_value(&atom.auth_asym_id, "_atom_site.auth_asym_id")?,
                cif_value(&atom.auth_atom_id, "_atom_site.auth_atom_id")?,
                model_number.to_string(),
            ],
        )?;
        *serial += 1;
    }
    Ok(())
}

pub(super) fn write_block_end(
    writer: &mut impl Write,
    model: &PreparedModel,
) -> Result<(), MmcifWriteError> {
    writer.write_all(b"#\n").map_err(mmcif_io)?;
    if !model.connections.is_empty() {
        let indexes = model
            .atoms
            .iter()
            .enumerate()
            .map(|(index, row)| (row.atom, index))
            .collect::<BTreeMap<_, _>>();
        const CONNECTION_TAGS: &[&str] = &[
            "_struct_conn.id",
            "_struct_conn.conn_type_id",
            "_struct_conn.ptnr1_label_asym_id",
            "_struct_conn.ptnr1_label_comp_id",
            "_struct_conn.ptnr1_label_seq_id",
            "_struct_conn.ptnr1_label_atom_id",
            "_struct_conn.ptnr2_label_asym_id",
            "_struct_conn.ptnr2_label_comp_id",
            "_struct_conn.ptnr2_label_seq_id",
            "_struct_conn.ptnr2_label_atom_id",
            "_struct_conn.pdbx_value_order",
        ];
        write_loop_header(writer, CONNECTION_TAGS)?;
        for (serial, connection) in (1u64..).zip(model.connections.iter()) {
            let left = &model.atoms[indexes[&connection.left]];
            let right = &model.atoms[indexes[&connection.right]];
            write_row(
                writer,
                vec![
                    format!("conn{serial}"),
                    "covale".to_owned(),
                    cif_value(&left.asym_id, "_struct_conn.ptnr1_label_asym_id")?,
                    cif_value(&left.label_comp_id, "_struct_conn.ptnr1_label_comp_id")?,
                    optional_display(left.label_seq_id),
                    cif_value(&left.label_atom_id, "_struct_conn.ptnr1_label_atom_id")?,
                    cif_value(&right.asym_id, "_struct_conn.ptnr2_label_asym_id")?,
                    cif_value(&right.label_comp_id, "_struct_conn.ptnr2_label_comp_id")?,
                    optional_display(right.label_seq_id),
                    cif_value(&right.label_atom_id, "_struct_conn.ptnr2_label_atom_id")?,
                    bond_order_code(connection.order).to_owned(),
                ],
            )?;
        }
        writer.write_all(b"#\n").map_err(mmcif_io)?;
    }
    Ok(())
}

fn write_loop_header(writer: &mut impl Write, tags: &[&str]) -> Result<(), MmcifWriteError> {
    writer.write_all(b"loop_\n").map_err(mmcif_io)?;
    for tag in tags {
        writer.write_all(tag.as_bytes()).map_err(mmcif_io)?;
        writer.write_all(b"\n").map_err(mmcif_io)?;
    }
    Ok(())
}

fn write_row(writer: &mut impl Write, values: Vec<String>) -> Result<(), MmcifWriteError> {
    for (index, value) in values.into_iter().enumerate() {
        if index != 0 {
            writer.write_all(b" ").map_err(mmcif_io)?;
        }
        writer.write_all(value.as_bytes()).map_err(mmcif_io)?;
    }
    writer.write_all(b"\n").map_err(mmcif_io)
}

fn format_coordinate(value: f64, precision: usize) -> String {
    format!("{value:.precision$}")
}

fn optional_float(value: Option<f64>) -> String {
    value.map_or_else(|| ".".to_owned(), |value| value.to_string())
}

fn optional_display(value: Option<i32>) -> String {
    value.map_or_else(|| ".".to_owned(), |value| value.to_string())
}

fn optional_cif_value(value: Option<&str>, field: &'static str) -> Result<String, MmcifWriteError> {
    value.map_or_else(|| Ok(".".to_owned()), |value| cif_value(value, field))
}

fn cif_value(value: &str, field: &'static str) -> Result<String, MmcifWriteError> {
    if value.is_empty() || value.contains(['\n', '\r']) {
        return Err(MmcifWriteError::UnsupportedTextValue { field });
    }
    let lower = value.to_ascii_lowercase();
    let is_control = lower == "loop_"
        || lower == "stop_"
        || lower == "global_"
        || lower.starts_with("data_")
        || lower.starts_with("save_")
        || value.starts_with('_');
    let bare = !value.is_empty()
        && value != "."
        && value != "?"
        && !is_control
        && !value
            .chars()
            .any(|character| character.is_ascii_whitespace() || character == '#')
        && !value.starts_with(';')
        && !value.contains(['\'', '"']);
    if bare {
        return Ok(value.to_owned());
    }
    if !value.contains('\'') {
        return Ok(format!("'{value}'"));
    }
    if !value.contains('"') {
        return Ok(format!("\"{value}\""));
    }
    Err(MmcifWriteError::UnsupportedTextValue { field })
}

fn bond_order_code(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "sing",
        BondOrder::Double => "doub",
        BondOrder::Triple => "trip",
        BondOrder::Quadruple => "quad",
        BondOrder::Zero | BondOrder::Dative => {
            unreachable!("unsupported bond order was rejected")
        }
    }
}
