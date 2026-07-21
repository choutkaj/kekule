use std::collections::BTreeMap;
use std::fmt;

use crate::algorithms::explicit_valence;
use crate::core::*;
use crate::io::preserve_molfile_tetrahedral_hydrogens;
use crate::small::model::SmallMolecule;
use crate::units::{Quantity, ANGSTROM};

use super::sdf_document::{SdfDataField, SdfRecord};

pub(super) const V2000_MAX_ATOMS: usize = 999;
pub(super) const V2000_MAX_BONDS: usize = 999;

#[derive(Debug, Clone, PartialEq)]
pub(super) struct V2000Syntax {
    pub(super) atoms: Vec<V2000AtomSyntax>,
    pub(super) bonds: Vec<V2000BondSyntax>,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct V2000AtomSyntax {
    pub(super) atom: Atom,
    pub(super) point: Point3,
    pub(super) line: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct V2000BondSyntax {
    pub(super) a: usize,
    pub(super) b: usize,
    pub(super) order: BondOrder,
    pub(super) stereo: Option<StereoBondMarkKind>,
    pub(super) line: usize,
}

/// Resource and record-boundary policy for SDF parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SdfParseOptions {
    /// Accept a nonempty final record without its terminating `$$$$` line.
    pub allow_missing_final_delimiter: bool,
    /// Maximum UTF-8 byte length of the complete input document.
    pub max_input_bytes: usize,
    /// Maximum number of nonempty records.
    pub max_records: usize,
    /// Maximum normalized byte length of one record, excluding its delimiter.
    pub max_record_bytes: usize,
}

impl Default for SdfParseOptions {
    fn default() -> Self {
        Self {
            allow_missing_final_delimiter: false,
            max_input_bytes: 256 * 1024 * 1024,
            max_records: 1_000_000,
            max_record_bytes: 64 * 1024 * 1024,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfParseError {
    pub(crate) record: usize,
    pub(crate) line: usize,
    pub(crate) message: String,
}

impl SdfParseError {
    pub(crate) fn new(record: usize, line: usize, message: impl Into<String>) -> Self {
        Self {
            record,
            line,
            message: message.into(),
        }
    }

    pub const fn record(&self) -> usize {
        self.record
    }

    pub const fn line(&self) -> usize {
        self.line
    }

    pub fn message(&self) -> &str {
        &self.message
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

pub(super) fn parse_v2000_syntax(
    record: usize,
    start_line: usize,
    lines: &[&str],
) -> std::result::Result<V2000Syntax, SdfParseError> {
    if lines.len() < 4 {
        return Err(SdfParseError::new(
            record,
            start_line,
            "record must contain three header lines and a counts line",
        ));
    }
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

    let mut atoms = Vec::with_capacity(atom_count);
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
        apply_atom_v2000_fields(record, line_number, &mut atom, atom_line)?;
        let point = atom_coordinates_from_v2000_line(atom_line)
            .ok_or_else(|| SdfParseError::new(record, line_number, "invalid atom coordinates"))?;
        atoms.push(V2000AtomSyntax {
            atom,
            point,
            line: line_number,
        });
    }

    let mut bonds = Vec::with_capacity(bond_count);
    let mut endpoints = std::collections::BTreeSet::new();
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
        if a_index >= atoms.len() || b_index >= atoms.len() {
            return Err(SdfParseError::new(
                record,
                line_number,
                "bond endpoint outside atom block",
            ));
        }
        if a_index == b_index {
            return Err(SdfParseError::new(
                record,
                line_number,
                "bond endpoints must be distinct",
            ));
        }
        let ordered = if a_index < b_index {
            (a_index, b_index)
        } else {
            (b_index, a_index)
        };
        if !endpoints.insert(ordered) {
            return Err(SdfParseError::new(
                record,
                line_number,
                "duplicate bond endpoints",
            ));
        }
        bonds.push(V2000BondSyntax {
            a: a_index,
            b: b_index,
            order,
            stereo,
            line: line_number,
        });
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
        &mut atoms,
        &lines[property_start..end_index],
    )?;

    Ok(V2000Syntax { atoms, bonds })
}

pub(super) fn interpret_v2000_syntax(
    syntax: &V2000Syntax,
) -> std::result::Result<SmallMolecule, SdfParseError> {
    let mut mol = Molecule::new();
    let mut atom_ids = Vec::with_capacity(syntax.atoms.len());
    let mut conformer = Conformer::with_atom_capacity(syntax.atoms.len(), ANGSTROM)
        .expect("angstrom is a length unit");
    for record in &syntax.atoms {
        let atom_id = mol.add_atom(record.atom.clone());
        conformer
            .set_position(atom_id, Quantity::new(record.point, ANGSTROM))
            .expect("matching coordinate units");
        atom_ids.push(atom_id);
    }
    for bond in &syntax.bonds {
        let a = atom_ids.get(bond.a).copied().ok_or_else(|| {
            SdfParseError::new(1, bond.line, "bond endpoint outside parsed atom records")
        })?;
        let b = atom_ids.get(bond.b).copied().ok_or_else(|| {
            SdfParseError::new(1, bond.line, "bond endpoint outside parsed atom records")
        })?;
        let bond_id = mol.add_bond(a, b, bond.order).map_err(|error| {
            SdfParseError::new(1, bond.line, format!("invalid graph bond: {error}"))
        })?;
        if let Some(kind) = bond.stereo {
            mol.set_stereo_bond_mark(StereoBondMark {
                bond: bond_id,
                kind,
                source: StereoSource::MolfileV2000,
            })
            .expect("newly added bond should accept a stereo mark");
        }
    }

    preserve_molfile_tetrahedral_hydrogens(&mut mol);
    if conformer.positions().next().is_some() {
        mol.add_conformer(conformer)
            .expect("parsed coordinates reference live atoms");
    }

    Ok(SmallMolecule::from_graph(mol))
}

pub(super) fn parse_counts_line(line: &str) -> Option<(usize, usize)> {
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

fn apply_atom_v2000_fields(
    record: usize,
    line_number: usize,
    atom: &mut Atom,
    line: &str,
) -> std::result::Result<(), SdfParseError> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if let Some(value) = fields.get(5) {
        let charge_code = value.parse::<u8>().map_err(|_| {
            SdfParseError::new(record, line_number, "invalid V2000 atom charge code")
        })?;
        atom.formal_charge = match charge_code {
            0 | 4 => 0,
            1 => 3,
            2 => 2,
            3 => 1,
            5 => -1,
            6 => -2,
            7 => -3,
            _ => {
                return Err(SdfParseError::new(
                    record,
                    line_number,
                    "unsupported V2000 atom charge code",
                ))
            }
        };
        if charge_code == 4 {
            atom.radical = Some(AtomRadical::Doublet);
        }
    }
    if fields
        .get(9)
        .and_then(|value| value.parse::<u8>().ok())
        .is_some_and(|valence| (1..=15).contains(&valence))
    {
        atom.no_implicit_hydrogens = true;
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
    Ok(())
}

fn parse_v2000_bond_line(
    line: &str,
) -> Option<(usize, usize, BondOrder, Option<StereoBondMarkKind>)> {
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
        (BondOrder::Single, 1) => Some(StereoBondMarkKind::WedgeUp),
        (BondOrder::Single, 4) => Some(StereoBondMarkKind::WedgeEither),
        (BondOrder::Single, 6) => Some(StereoBondMarkKind::WedgeDown),
        (BondOrder::Double, 3) => Some(StereoBondMarkKind::DoubleBondEither),
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

fn parse_m_records(
    record: usize,
    start_line: usize,
    atoms: &mut [V2000AtomSyntax],
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
                    atoms,
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
                    atoms,
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
                    atoms,
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
    atoms: &mut [V2000AtomSyntax],
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
        let atom = atoms
            .get_mut(atom_offset)
            .map(|record| &mut record.atom)
            .ok_or_else(|| SdfParseError::new(record, line, "M record atom outside atom block"))?;
        apply(atom, value).map_err(|message| SdfParseError::new(record, line, message))?;
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolWriteError {
    pub(crate) message: String,
}

impl MolWriteError {
    pub(crate) fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for MolWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for MolWriteError {}

pub fn write_mol_v2000(molecule: &SmallMolecule) -> std::result::Result<String, MolWriteError> {
    let mol = molecule.graph();
    if mol.stereo_elements().next().is_some() {
        return Err(MolWriteError::new(
            "V2000 writer does not support stereo elements",
        ));
    }
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

    let title = "";
    let program = "molecular";
    let comment = "";
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
            .map(|point| point.value_in(ANGSTROM).expect("conformer length unit"))
            .unwrap_or_default();
        let valence_code = v2000_valence_code(mol, *atom_id, atom)?;
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3}{:>2}{:>3}  0  0  0{:>3}  0  0  0{:>3}  0  0\n",
            point.x,
            point.y,
            point.z,
            atom.element.symbol(),
            0,
            v2000_charge_code(atom.formal_charge),
            valence_code,
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
        let stereo_code = v2000_bond_stereo_code(bond.order, mol.stereo_bond_mark(*bond_id))?;
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
    let mut radical_records = Vec::new();
    for id in &atoms {
        let atom = mol
            .atom(*id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let Some(radical) = atom.radical else {
            continue;
        };
        let index = *atom_index
            .get(id)
            .ok_or_else(|| MolWriteError::new("atom missing from V2000 atom table"))?;
        radical_records.push((index as i32, v2000_radical_code(radical)?));
    }
    push_m_record(&mut out, "RAD", radical_records);
    out.push_str("M  END\n");
    Ok(out)
}

pub fn write_sdf_v2000(records: &[SdfRecord]) -> std::result::Result<String, MolWriteError> {
    let mut out = String::new();
    for record in records {
        validate_sdf_title(record.title())?;
        for field in record.data_fields() {
            validate_sdf_data_field(field)?;
        }
        let written = write_mol_v2000(record.molecule())?;
        let mut lines = written.lines();
        let _generated_title = lines.next();
        out.push_str(record.title());
        out.push('\n');
        for line in lines {
            out.push_str(line);
            out.push('\n');
        }
        for field in record.data_fields() {
            out.push_str(&format!(">  <{}>\n{}\n\n", field.name(), field.value()));
        }
        out.push_str("$$$$\n");
    }
    Ok(out)
}

fn validate_sdf_title(title: &str) -> std::result::Result<(), MolWriteError> {
    if title.contains(['\r', '\n']) {
        return Err(MolWriteError::new(
            "SDF record titles cannot contain line breaks",
        ));
    }
    Ok(())
}

fn validate_sdf_data_field(field: &SdfDataField) -> std::result::Result<(), MolWriteError> {
    let name = field.name();
    if name.is_empty() || name.trim() != name || name.contains(['<', '>', '\r', '\n']) {
        return Err(MolWriteError::new(
            "SDF data field names must be nonempty, trimmed, and exclude angle brackets or line breaks",
        ));
    }
    let value = field.value();
    if value.contains('\r') {
        return Err(MolWriteError::new(
            "SDF data field values cannot contain carriage returns",
        ));
    }
    if !value.is_empty() && value.split('\n').any(str::is_empty) {
        return Err(MolWriteError::new(
            "SDF data field values cannot contain blank lines",
        ));
    }
    if value.lines().any(|line| line.trim() == "$$$$") {
        return Err(MolWriteError::new(
            "SDF data field values cannot contain a record delimiter line",
        ));
    }
    Ok(())
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

fn v2000_valence_code(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
) -> std::result::Result<u8, MolWriteError> {
    if !atom.no_implicit_hydrogens {
        return Ok(0);
    }
    let valence = explicit_valence(mol, atom_id) + usize::from(atom.explicit_hydrogens);
    match valence {
        0 => Ok(15),
        1..=14 => Ok(u8::try_from(valence).expect("range checked")),
        _ => Err(MolWriteError::new(format!(
            "V2000 cannot encode explicit valence {valence} for atom {}",
            atom_id.index()
        ))),
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
    stereo: Option<&StereoBondMark>,
) -> std::result::Result<u8, MolWriteError> {
    match (order, stereo.map(|mark| mark.kind)) {
        (_, None) => Ok(0),
        (BondOrder::Single, Some(StereoBondMarkKind::WedgeUp)) => Ok(1),
        (BondOrder::Single, Some(StereoBondMarkKind::WedgeEither)) => Ok(4),
        (BondOrder::Single, Some(StereoBondMarkKind::WedgeDown)) => Ok(6),
        (BondOrder::Double, Some(StereoBondMarkKind::DoubleBondEither)) => Ok(3),
        _ => Err(MolWriteError::new(
            "V2000 bond stereo is incompatible with the bond order",
        )),
    }
}

fn v2000_radical_code(radical: AtomRadical) -> std::result::Result<i32, MolWriteError> {
    match radical {
        AtomRadical::Singlet => Ok(1),
        AtomRadical::Doublet => Ok(2),
        AtomRadical::Triplet => Ok(3),
        AtomRadical::Quartet | AtomRadical::Quintet => Err(MolWriteError::new(
            "V2000 writer cannot encode radical multiplicity above triplet",
        )),
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
