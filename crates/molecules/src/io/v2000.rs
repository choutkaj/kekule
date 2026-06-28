use std::collections::BTreeMap;
use std::fmt;

use crate::core::*;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MolParseOptions;

const V2000_MAX_ATOMS: usize = 999;
const V2000_MAX_BONDS: usize = 999;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfParseOptions {
    pub allow_missing_final_delimiter: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    pub title: String,
    pub molecule: SmallMolecule,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfParseError {
    pub record: usize,
    pub line: usize,
    pub message: String,
}

impl SdfParseError {
    pub(crate) fn new(record: usize, line: usize, message: impl Into<String>) -> Self {
        Self {
            record,
            line,
            message: message.into(),
        }
    }
}

impl fmt::Display for SdfParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SDF parse error in record {} at line {}: {}",
            self.record, self.line, self.message
        )
    }
}

impl std::error::Error for SdfParseError {}

pub fn read_mol_v2000_str(input: &str) -> std::result::Result<SmallMolecule, SdfParseError> {
    let normalized = input.replace("\r\n", "\n").replace('\r', "\n");
    let lines = normalized.lines().collect::<Vec<_>>();
    parse_mol_v2000_lines(1, 1, &lines)
}

pub fn read_sdf_v2000_str(
    input: &str,
    options: SdfParseOptions,
) -> std::result::Result<Vec<SmallMolecule>, SdfParseError> {
    read_sdf_v2000_records(input, options)
        .map(|records| records.into_iter().map(|record| record.molecule).collect())
}

pub fn read_sdf_v2000_records(
    input: &str,
    options: SdfParseOptions,
) -> std::result::Result<Vec<SdfRecord>, SdfParseError> {
    let normalized = input.replace("\r\n", "\n").replace('\r', "\n");
    let mut records = Vec::new();
    let mut current = Vec::new();
    let mut start_line = 1usize;
    let mut saw_delimiter = false;

    for (offset, line) in normalized.lines().enumerate() {
        let line_number = offset + 1;
        if line.trim() == "$$$$" {
            saw_delimiter = true;
            if current.iter().any(|line: &&str| !line.trim().is_empty()) {
                records.push(parse_sdf_record(records.len() + 1, start_line, &current)?);
            }
            current.clear();
            start_line = line_number + 1;
        } else {
            current.push(line);
        }
    }

    if current.iter().any(|line| !line.trim().is_empty()) {
        if saw_delimiter || options.allow_missing_final_delimiter {
            records.push(parse_sdf_record(records.len() + 1, start_line, &current)?);
        } else {
            let final_offset = current.len().saturating_sub(1);
            return Err(SdfParseError::new(
                records.len() + 1,
                checked_line_number(records.len() + 1, start_line, final_offset)?,
                "missing final $$$$ record delimiter",
            ));
        }
    }

    Ok(records)
}

fn parse_sdf_record(
    record: usize,
    start_line: usize,
    lines: &[&str],
) -> std::result::Result<SdfRecord, SdfParseError> {
    let title = lines.first().copied().unwrap_or_default().to_owned();
    let end_index = lines
        .iter()
        .position(|line| line.trim() == "M  END")
        .ok_or_else(|| SdfParseError::new(record, start_line, "missing M  END"))?;
    let mut molecule = parse_mol_v2000_lines(record, start_line, &lines[..=end_index])?;
    let data_start = end_index
        .checked_add(1)
        .ok_or_else(|| SdfParseError::new(record, start_line, "SDF data offset overflow"))?;
    parse_sdf_data_fields(
        record,
        checked_line_number(record, start_line, data_start)?,
        &mut molecule.mol,
        &lines[data_start..],
    )?;
    Ok(SdfRecord { title, molecule })
}

fn parse_mol_v2000_lines(
    record: usize,
    start_line: usize,
    lines: &[&str],
) -> std::result::Result<SmallMolecule, SdfParseError> {
    if lines.len() < 4 {
        return Err(SdfParseError::new(
            record,
            start_line,
            "record must contain three header lines and a counts line",
        ));
    }
    let title = lines[0].to_owned();
    let counts = lines[3];
    let counts_line = checked_line_number(record, start_line, 3)?;
    if counts.contains("V3000") {
        return Err(SdfParseError::new(
            record,
            counts_line,
            "V3000 records are not supported by the V2000 parser",
        ));
    }
    if !counts.contains("V2000") {
        return Err(SdfParseError::new(
            record,
            counts_line,
            "counts line must declare V2000",
        ));
    }
    let (atom_count, bond_count) = parse_counts_line(counts)
        .ok_or_else(|| SdfParseError::new(record, counts_line, "invalid V2000 counts line"))?;
    if atom_count > V2000_MAX_ATOMS || bond_count > V2000_MAX_BONDS {
        return Err(SdfParseError::new(
            record,
            counts_line,
            "V2000 counts exceed the supported 999 atom or bond limit",
        ));
    }

    let mut mol = Molecule::new();
    mol.props_mut()
        .insert("sdf.title".to_owned(), PropValue::String(title.clone()));
    mol.props_mut().insert(
        "sdf.program".to_owned(),
        PropValue::String(lines[1].to_owned()),
    );
    mol.props_mut().insert(
        "sdf.comment".to_owned(),
        PropValue::String(lines[2].to_owned()),
    );

    let atom_start = 4usize;
    let bond_start = atom_start.checked_add(atom_count).ok_or_else(|| {
        SdfParseError::new(record, start_line, "V2000 atom block offset overflow")
    })?;
    let property_start = bond_start.checked_add(bond_count).ok_or_else(|| {
        SdfParseError::new(record, start_line, "V2000 bond block offset overflow")
    })?;
    if lines.len() < property_start {
        return Err(SdfParseError::new(
            record,
            checked_line_number(record, start_line, lines.len())?,
            "record ended before declared atom and bond blocks",
        ));
    }

    let mut atom_ids = Vec::with_capacity(atom_count);
    let mut conformer = Conformer::with_atom_capacity(atom_count);
    for atom_index in 0..atom_count {
        let block_index = atom_start
            .checked_add(atom_index)
            .ok_or_else(|| SdfParseError::new(record, start_line, "V2000 atom index overflow"))?;
        let line_number = checked_line_number(record, start_line, block_index)?;
        let atom_line = lines
            .get(block_index)
            .copied()
            .ok_or_else(|| SdfParseError::new(record, line_number, "truncated V2000 atom block"))?;
        let symbol = atom_symbol_from_v2000_line(atom_line)
            .ok_or_else(|| SdfParseError::new(record, line_number, "invalid atom line"))?;
        let element = Element::from_symbol(symbol).ok_or_else(|| {
            SdfParseError::new(
                record,
                line_number,
                format!("unknown element symbol `{symbol}`"),
            )
        })?;
        let mut atom = Atom::new(element);
        apply_atom_v2000_fields(&mut atom, atom_line);
        let point = atom_coordinates_from_v2000_line(atom_line)
            .ok_or_else(|| SdfParseError::new(record, line_number, "invalid atom coordinates"))?;
        let atom_id = mol.add_atom(atom);
        conformer.set_position(atom_id, point);
        atom_ids.push(atom_id);
    }

    for bond_index in 0..bond_count {
        let block_index = bond_start
            .checked_add(bond_index)
            .ok_or_else(|| SdfParseError::new(record, start_line, "V2000 bond index overflow"))?;
        let line_number = checked_line_number(record, start_line, block_index)?;
        let bond_line = lines
            .get(block_index)
            .copied()
            .ok_or_else(|| SdfParseError::new(record, line_number, "truncated V2000 bond block"))?;
        let (a, b, order, stereo) = parse_v2000_bond_line(bond_line)
            .ok_or_else(|| SdfParseError::new(record, line_number, "invalid bond line"))?;
        let a_index = a.checked_sub(1).ok_or_else(|| {
            SdfParseError::new(record, line_number, "bond endpoint must be one-based")
        })?;
        let b_index = b.checked_sub(1).ok_or_else(|| {
            SdfParseError::new(record, line_number, "bond endpoint must be one-based")
        })?;
        let a = atom_ids.get(a_index).copied().ok_or_else(|| {
            SdfParseError::new(record, line_number, "bond endpoint outside atom block")
        })?;
        let b = atom_ids.get(b_index).copied().ok_or_else(|| {
            SdfParseError::new(record, line_number, "bond endpoint outside atom block")
        })?;
        let bond_id = mol.add_bond(a, b, order).map_err(|error| {
            SdfParseError::new(record, line_number, format!("invalid graph bond: {error}"))
        })?;
        if let Some(stereo) = stereo {
            mol.bond_mut(bond_id)
                .expect("newly added bond should be mutable")
                .stereo = Some(stereo);
        }
    }

    let property_line = checked_line_number(record, start_line, property_start)?;
    let relative_end = lines[property_start..]
        .iter()
        .position(|line| line.trim() == "M  END")
        .ok_or_else(|| SdfParseError::new(record, property_line, "missing M  END"))?;
    let end_index = property_start.checked_add(relative_end).ok_or_else(|| {
        SdfParseError::new(
            record,
            property_line,
            "V2000 property block offset overflow",
        )
    })?;
    parse_m_records(
        record,
        property_line,
        &mut mol,
        &atom_ids,
        &lines[property_start..end_index],
    )?;
    if conformer.positions().next().is_some() {
        mol.add_conformer(conformer);
    }

    Ok(SmallMolecule { mol })
}

fn parse_counts_line(line: &str) -> Option<(usize, usize)> {
    if !line.is_ascii() {
        return None;
    }
    if let (Some(atom_field), Some(bond_field)) = (ascii_field(line, 0, 3), ascii_field(line, 3, 6))
    {
        if let (Ok(atoms), Ok(bonds)) = (atom_field.trim().parse(), bond_field.trim().parse()) {
            return Some((atoms, bonds));
        }
    }
    let fields = line.split_whitespace().collect::<Vec<_>>();
    Some((fields.first()?.parse().ok()?, fields.get(1)?.parse().ok()?))
}

fn atom_symbol_from_v2000_line(line: &str) -> Option<&str> {
    if !line.is_ascii() {
        return None;
    }
    ascii_field(line, 31, 34)
        .map(str::trim)
        .filter(|symbol| !symbol.is_empty())
        .or_else(|| line.split_whitespace().nth(3))
}

fn atom_coordinates_from_v2000_line(line: &str) -> Option<Point3> {
    if !line.is_ascii() {
        return None;
    }
    if let (Some(x), Some(y), Some(z)) = (
        ascii_field(line, 0, 10),
        ascii_field(line, 10, 20),
        ascii_field(line, 20, 30),
    ) {
        if let (Ok(x), Ok(y), Ok(z)) = (
            x.trim().parse::<f64>(),
            y.trim().parse::<f64>(),
            z.trim().parse::<f64>(),
        ) {
            if x.is_finite() && y.is_finite() && z.is_finite() {
                return Some(Point3::new(x, y, z));
            }
        }
    }
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let point = Point3::new(
        fields.first()?.parse().ok()?,
        fields.get(1)?.parse().ok()?,
        fields.get(2)?.parse().ok()?,
    );
    (point.x.is_finite() && point.y.is_finite() && point.z.is_finite()).then_some(point)
}

fn apply_atom_v2000_fields(atom: &mut Atom, line: &str) {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if let Some(charge_code) = fields.get(5).and_then(|value| value.parse::<i8>().ok()) {
        atom.formal_charge = match charge_code {
            1 => 3,
            2 => 2,
            3 => 1,
            5 => -1,
            6 => -2,
            7 => -3,
            _ => 0,
        };
    }
    if let Some(atom_map) = fields
        .get(13)
        .or_else(|| fields.get(12))
        .and_then(|value| value.parse::<u32>().ok())
    {
        if atom_map != 0 {
            atom.atom_map = Some(atom_map);
        }
    }
}

fn parse_v2000_bond_line(line: &str) -> Option<(usize, usize, BondOrder, Option<BondStereo>)> {
    if !line.is_ascii() {
        return None;
    }
    let (a, b, order_code, stereo_code) = if let (Some(a), Some(b), Some(order), Some(stereo)) = (
        ascii_field(line, 0, 3),
        ascii_field(line, 3, 6),
        ascii_field(line, 6, 9),
        ascii_field(line, 9, 12),
    ) {
        (
            a.trim().parse().ok()?,
            b.trim().parse().ok()?,
            order.trim().parse().ok()?,
            stereo.trim().parse::<u8>().ok(),
        )
    } else {
        let mut fields = line.split_whitespace();
        (
            fields.next()?.parse().ok()?,
            fields.next()?.parse().ok()?,
            fields.next()?.parse().ok()?,
            fields.next().and_then(|value| value.parse::<u8>().ok()),
        )
    };
    let order = match order_code {
        0 => BondOrder::Zero,
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        9 => BondOrder::Dative,
        _ => return None,
    };
    let stereo_code = stereo_code.unwrap_or(0);
    let stereo = match (order, stereo_code) {
        (_, 0) => None,
        (BondOrder::Single, 1) => Some(BondStereo::Up),
        (BondOrder::Single, 4) => Some(BondStereo::Any),
        (BondOrder::Single, 6) => Some(BondStereo::Down),
        (BondOrder::Double, 3) => Some(BondStereo::Any),
        _ => return None,
    };
    Some((a, b, order, stereo))
}

fn ascii_field(line: &str, start: usize, end: usize) -> Option<&str> {
    std::str::from_utf8(line.as_bytes().get(start..end)?).ok()
}

fn checked_line_number(
    record: usize,
    start_line: usize,
    offset: usize,
) -> std::result::Result<usize, SdfParseError> {
    start_line
        .checked_add(offset)
        .ok_or_else(|| SdfParseError::new(record, start_line, "line number overflow"))
}

fn parse_sdf_data_fields(
    record: usize,
    start_line: usize,
    mol: &mut Molecule,
    lines: &[&str],
) -> std::result::Result<(), SdfParseError> {
    let mut index = 0;
    while index < lines.len() {
        let line = lines[index];
        if !line.trim_start().starts_with('>') {
            index += 1;
            continue;
        }
        let field_name = sdf_field_name(line).ok_or_else(|| {
            SdfParseError::new(record, start_line + index, "invalid SDF data field header")
        })?;
        index += 1;
        let mut values = Vec::new();
        while index < lines.len() && !lines[index].trim_start().starts_with('>') {
            if lines[index].is_empty() {
                index += 1;
                break;
            }
            values.push(lines[index]);
            index += 1;
        }
        mol.props_mut().insert(
            format!("sdf.field.{field_name}"),
            PropValue::String(values.join("\n")),
        );
    }
    Ok(())
}

fn parse_m_records(
    record: usize,
    start_line: usize,
    mol: &mut Molecule,
    atom_ids: &[AtomId],
    lines: &[&str],
) -> std::result::Result<(), SdfParseError> {
    for (offset, line) in lines.iter().enumerate() {
        let fields = line.split_whitespace().collect::<Vec<_>>();
        match fields.as_slice() {
            ["M", "CHG", count, rest @ ..] => {
                parse_atom_value_pairs(
                    record,
                    start_line + offset,
                    count,
                    rest,
                    mol,
                    atom_ids,
                    |atom, value| {
                        atom.formal_charge =
                            i8::try_from(value).map_err(|_| "formal charge is outside i8 range")?;
                        Ok(())
                    },
                )?;
            }
            ["M", "ISO", count, rest @ ..] => {
                parse_atom_value_pairs(
                    record,
                    start_line + offset,
                    count,
                    rest,
                    mol,
                    atom_ids,
                    |atom, value| {
                        atom.isotope = if value > 0 {
                            Some(u16::try_from(value).map_err(|_| "isotope is outside u16 range")?)
                        } else {
                            None
                        };
                        Ok(())
                    },
                )?;
            }
            ["M", "RAD", count, rest @ ..] => {
                parse_atom_value_pairs(
                    record,
                    start_line + offset,
                    count,
                    rest,
                    mol,
                    atom_ids,
                    |atom, value| {
                        atom.radical = Some(match value {
                            1 => AtomRadical::Singlet,
                            2 => AtomRadical::Doublet,
                            3 => AtomRadical::Triplet,
                            _ => return Err("unsupported M  RAD code"),
                        });
                        Ok(())
                    },
                )?;
            }
            _ => {}
        }
    }
    Ok(())
}

fn parse_atom_value_pairs<F>(
    record: usize,
    line: usize,
    count: &str,
    rest: &[&str],
    mol: &mut Molecule,
    atom_ids: &[AtomId],
    mut apply: F,
) -> std::result::Result<(), SdfParseError>
where
    F: FnMut(&mut Atom, i32) -> std::result::Result<(), &'static str>,
{
    let count = count
        .parse::<usize>()
        .map_err(|_| SdfParseError::new(record, line, "invalid M record count"))?;
    let pair_fields = count
        .checked_mul(2)
        .ok_or_else(|| SdfParseError::new(record, line, "M record pair count overflow"))?;
    if rest.len() != pair_fields {
        return Err(SdfParseError::new(
            record,
            line,
            "M record pair count does not match its fields",
        ));
    }
    for pair in rest.chunks(2).take(count) {
        let atom_index = pair[0]
            .parse::<usize>()
            .map_err(|_| SdfParseError::new(record, line, "invalid M record atom index"))?;
        let value = pair[1]
            .parse::<i32>()
            .map_err(|_| SdfParseError::new(record, line, "invalid M record value"))?;
        let atom_offset = atom_index.checked_sub(1).ok_or_else(|| {
            SdfParseError::new(record, line, "M record atom index must be one-based")
        })?;
        let atom_id = atom_ids
            .get(atom_offset)
            .copied()
            .ok_or_else(|| SdfParseError::new(record, line, "M record atom outside atom block"))?;
        let mut atom = mol
            .atom_mut(atom_id)
            .map_err(|error| SdfParseError::new(record, line, error.to_string()))?;
        apply(&mut atom, value).map_err(|message| SdfParseError::new(record, line, message))?;
    }
    Ok(())
}

fn sdf_field_name(line: &str) -> Option<String> {
    let start = line.find('<')?;
    let end = line[start + 1..].find('>')? + start + 1;
    let name = line[start + 1..end].trim();
    if name.is_empty() {
        None
    } else {
        Some(name.to_owned())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolWriteError {
    pub message: String,
}

impl MolWriteError {
    pub(crate) fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl fmt::Display for MolWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for MolWriteError {}

pub fn write_mol_v2000(molecule: &SmallMolecule) -> std::result::Result<String, MolWriteError> {
    let mol = &molecule.mol;
    if mol.atom_count() > 999 || mol.bond_count() > 999 {
        return Err(MolWriteError::new(
            "V2000 writer supports at most 999 atoms and 999 bonds",
        ));
    }
    let atoms = mol.atom_ids().collect::<Vec<_>>();
    let bonds = mol.bond_ids().collect::<Vec<_>>();
    let mut atom_index = BTreeMap::new();
    for (index, atom_id) in atoms.iter().enumerate() {
        atom_index.insert(*atom_id, index + 1);
    }

    let title = prop_string(mol, "sdf.title").unwrap_or_default();
    let program = prop_string(mol, "sdf.program").unwrap_or_else(|| "molecules".to_owned());
    let comment = prop_string(mol, "sdf.comment").unwrap_or_default();
    let conformer = mol.first_conformer().map(|(_, conformer)| conformer);
    let mut out = String::new();
    out.push_str(&format!("{title}\n{program}\n{comment}\n"));
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0            999 V2000\n",
        atoms.len(),
        bonds.len()
    ));

    for atom_id in &atoms {
        let atom = mol
            .atom(*atom_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let point = conformer
            .and_then(|conformer| conformer.position(*atom_id))
            .unwrap_or_default();
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3}{:>2}{:>3}  0  0  0  0  0  0  0{:>3}  0  0\n",
            point.x,
            point.y,
            point.z,
            atom.element.symbol(),
            0,
            v2000_charge_code(atom.formal_charge),
            atom.atom_map.unwrap_or(0)
        ));
    }

    for bond_id in &bonds {
        let bond = mol
            .bond(*bond_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let a = atom_index
            .get(&bond.a())
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from atom table"))?;
        let b = atom_index
            .get(&bond.b())
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from atom table"))?;
        let order_code = v2000_bond_code(bond.order)?;
        let stereo_code = v2000_bond_stereo_code(bond.order, bond.stereo)?;
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}{:>3}  0  0  0\n",
            a, b, order_code, stereo_code
        ));
    }

    push_m_record(
        &mut out,
        "CHG",
        atoms
            .iter()
            .filter_map(|id| {
                let atom = mol.atom(*id).ok()?;
                (atom.formal_charge != 0)
                    .then_some((*atom_index.get(id)? as i32, atom.formal_charge as i32))
            })
            .collect(),
    );
    push_m_record(
        &mut out,
        "ISO",
        atoms
            .iter()
            .filter_map(|id| {
                let atom = mol.atom(*id).ok()?;
                atom.isotope.map(|isotope| {
                    (
                        *atom_index.get(id).expect("atom indexed") as i32,
                        isotope as i32,
                    )
                })
            })
            .collect(),
    );
    push_m_record(
        &mut out,
        "RAD",
        atoms
            .iter()
            .filter_map(|id| {
                let atom = mol.atom(*id).ok()?;
                let index = *atom_index.get(id)?;
                atom.radical
                    .map(|radical| (index as i32, v2000_radical_code(radical)))
            })
            .collect(),
    );
    out.push_str("M  END\n");
    Ok(out)
}

pub fn write_sdf_v2000(molecules: &[SmallMolecule]) -> std::result::Result<String, MolWriteError> {
    let mut out = String::new();
    for molecule in molecules {
        out.push_str(&write_mol_v2000(molecule)?);
        for (key, value) in molecule.mol.props() {
            if let (Some(name), PropValue::String(text)) = (key.strip_prefix("sdf.field."), value) {
                out.push_str(&format!(">  <{name}>\n{text}\n\n"));
            }
        }
        out.push_str("$$$$\n");
    }
    Ok(out)
}

fn prop_string(mol: &Molecule, key: &str) -> Option<String> {
    match mol.props().get(key) {
        Some(PropValue::String(value)) => Some(value.clone()),
        _ => None,
    }
}

fn v2000_charge_code(charge: i8) -> i8 {
    match charge {
        3 => 1,
        2 => 2,
        1 => 3,
        -1 => 5,
        -2 => 6,
        -3 => 7,
        _ => 0,
    }
}

fn v2000_bond_code(order: BondOrder) -> std::result::Result<u8, MolWriteError> {
    match order {
        BondOrder::Zero => Ok(0),
        BondOrder::Single => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        BondOrder::Dative => Ok(9),
        BondOrder::Quadruple => Err(MolWriteError::new(
            "V2000 writer does not support quadruple bonds",
        )),
    }
}

fn v2000_bond_stereo_code(
    order: BondOrder,
    stereo: Option<BondStereo>,
) -> std::result::Result<u8, MolWriteError> {
    match (order, stereo) {
        (_, None | Some(BondStereo::Unspecified)) => Ok(0),
        (BondOrder::Single, Some(BondStereo::Up)) => Ok(1),
        (BondOrder::Single, Some(BondStereo::Any)) => Ok(4),
        (BondOrder::Single, Some(BondStereo::Down)) => Ok(6),
        (BondOrder::Double, Some(BondStereo::Any)) => Ok(3),
        (_, Some(BondStereo::E | BondStereo::Z)) => Err(MolWriteError::new(
            "V2000 cannot encode perceived E/Z bond stereo directly",
        )),
        _ => Err(MolWriteError::new(
            "V2000 bond stereo is incompatible with the bond order",
        )),
    }
}

fn v2000_radical_code(radical: AtomRadical) -> i32 {
    match radical {
        AtomRadical::Singlet => 1,
        AtomRadical::Doublet => 2,
        AtomRadical::Triplet => 3,
    }
}

fn push_m_record(out: &mut String, code: &str, pairs: Vec<(i32, i32)>) {
    for chunk in pairs.chunks(8) {
        if chunk.is_empty() {
            continue;
        }
        out.push_str(&format!("M  {code} {:>2}", chunk.len()));
        for (atom, value) in chunk {
            out.push_str(&format!("{atom:>4}{value:>4}"));
        }
        out.push('\n');
    }
}
