use std::collections::BTreeSet;
use std::fmt;

use crate::algorithms::explicit_valence;
use crate::chemistry::{
    localize_source_aromatic_bonds, SourceStereoBondMark, SourceStereoBondMarkKind,
};
use crate::core::*;
use crate::geometry::Point3;
use crate::structure::ModelView;
use crate::units::{Quantity, ANGSTROM};

use super::molfile_write::MolfileRecord;
use super::sdf_document::{SdfDataField, SdfRecordInterpretation};
use super::staged_coordinates::StagedCoordinates;
use super::structure_documents::{
    apply_molfile_declared_valence, checked_line_number, interpret_molfile_atom_fields,
    MolfileVersion,
};

pub(super) const V2000_MAX_ATOMS: usize = 999;
pub(super) const V2000_MAX_BONDS: usize = 999;

#[derive(Debug, Clone, PartialEq)]
pub(super) struct V2000Syntax {
    pub(super) atoms: Vec<V2000AtomSyntax>,
    pub(super) bonds: Vec<V2000BondSyntax>,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct V2000AtomSyntax {
    pub(super) symbol: String,
    pub(super) isotope: Option<i32>,
    pub(super) formal_charge: i32,
    pub(super) radical: Option<V2000RadicalSyntax>,
    pub(super) hydrogen_count: Option<u8>,
    pub(super) valence: Option<u8>,
    pub(super) atom_map: Option<u32>,
    pub(super) coordinates: [f64; 3],
    pub(super) line: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum V2000RadicalSyntax {
    Singlet,
    Doublet,
    Triplet,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct V2000BondSyntax {
    pub(super) a: usize,
    pub(super) b: usize,
    pub(super) order_code: u8,
    pub(super) stereo_code: u8,
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
        let mut atom = V2000AtomSyntax {
            symbol: symbol.to_owned(),
            isotope: None,
            formal_charge: 0,
            radical: None,
            hydrogen_count: None,
            valence: None,
            atom_map: None,
            coordinates: atom_coordinates_from_v2000_line(atom_line).ok_or_else(|| {
                SdfParseError::new(record, line_number, "invalid atom coordinates")
            })?,
            line: line_number,
        };
        apply_atom_v2000_fields(record, line_number, &mut atom, atom_line)?;
        atoms.push(atom);
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
        let (a, b, order_code, stereo_code) = parse_v2000_bond_line(bond_line)
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
            order_code,
            stereo_code,
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
) -> std::result::Result<
    (MoleculeEditor, StagedCoordinates, Vec<SourceStereoBondMark>),
    SdfParseError,
> {
    let mut editor = crate::core::MoleculeEditor::new();
    let mut atom_ids = Vec::with_capacity(syntax.atoms.len());
    let mut coordinates = StagedCoordinates::with_atom_capacity(syntax.atoms.len(), ANGSTROM)
        .map_err(|error| SdfParseError::new(1, 1, error.to_string()))?;
    for record in &syntax.atoms {
        let atom = interpret_v2000_atom(record)?;
        let atom_id = editor.add_atom(atom).map_err(|error| {
            SdfParseError::new(1, record.line, format!("invalid graph atom: {error}"))
        })?;
        coordinates
            .set_position(
                atom_id,
                Quantity::new(
                    Point3::new(
                        record.coordinates[0],
                        record.coordinates[1],
                        record.coordinates[2],
                    ),
                    ANGSTROM,
                ),
            )
            .map_err(|error| SdfParseError::new(1, record.line, error.to_string()))?;
        atom_ids.push(atom_id);
    }
    let mut source_aromatic_bonds = BTreeSet::new();
    let mut source_stereo = Vec::new();
    let mut first_aromatic_line = None;
    for bond in &syntax.bonds {
        let a = atom_ids.get(bond.a).copied().ok_or_else(|| {
            SdfParseError::new(1, bond.line, "bond endpoint outside parsed atom records")
        })?;
        let b = atom_ids.get(bond.b).copied().ok_or_else(|| {
            SdfParseError::new(1, bond.line, "bond endpoint outside parsed atom records")
        })?;
        let (order, source_aromatic) = interpret_v2000_bond_order(bond.order_code, bond.line)?;
        let bond_id = editor.add_bond(a, b, order).map_err(|error| {
            SdfParseError::new(1, bond.line, format!("invalid graph bond: {error}"))
        })?;
        if source_aromatic {
            source_aromatic_bonds.insert(bond_id);
            first_aromatic_line.get_or_insert(bond.line);
        }
        if let Some(kind) =
            interpret_v2000_bond_stereo(order, source_aromatic, bond.stereo_code, bond.line)?
        {
            source_stereo.push(SourceStereoBondMark {
                bond: bond_id,
                from: a,
                kind,
            });
        }
    }

    apply_v2000_declared_hydrogens(
        editor.working_mut(),
        syntax,
        &atom_ids,
        &source_aromatic_bonds,
    )?;
    localize_source_aromatic_bonds(editor.working_mut(), &source_aromatic_bonds).map_err(
        |error| {
            SdfParseError::new(
                1,
                first_aromatic_line.unwrap_or(1),
                format!("aromatic bond localization failed: {error}"),
            )
        },
    )?;
    Ok((editor, coordinates, source_stereo))
}

fn interpret_v2000_atom(record: &V2000AtomSyntax) -> std::result::Result<Atom, SdfParseError> {
    let radical = record.radical.map(|radical| match radical {
        V2000RadicalSyntax::Singlet => AtomRadical::Singlet,
        V2000RadicalSyntax::Doublet => AtomRadical::Doublet,
        V2000RadicalSyntax::Triplet => AtomRadical::Triplet,
    });
    let explicit = record.hydrogen_count.unwrap_or(0);
    let hydrogens = if record.hydrogen_count.is_some() || record.valence.is_some() {
        HydrogenDeclaration::Fixed(explicit)
    } else {
        HydrogenDeclaration::Infer { explicit }
    };
    interpret_molfile_atom_fields(
        &record.symbol,
        record.formal_charge,
        record.isotope,
        radical,
        hydrogens,
        record.atom_map,
        record.line,
    )
}

fn interpret_v2000_bond_order(
    code: u8,
    line: usize,
) -> std::result::Result<(BondOrder, bool), SdfParseError> {
    match code {
        0 => Ok((BondOrder::Zero, false)),
        1 => Ok((BondOrder::Single, false)),
        2 => Ok((BondOrder::Double, false)),
        3 => Ok((BondOrder::Triple, false)),
        4 => Ok((BondOrder::Single, true)),
        9 => Ok((BondOrder::Dative, false)),
        _ => Err(SdfParseError::new(
            1,
            line,
            format!("unsupported V2000 bond order code {code}"),
        )),
    }
}

fn interpret_v2000_bond_stereo(
    order: BondOrder,
    source_aromatic: bool,
    code: u8,
    line: usize,
) -> std::result::Result<Option<SourceStereoBondMarkKind>, SdfParseError> {
    if source_aromatic && code != 0 {
        return Err(SdfParseError::new(
            1,
            line,
            format!("V2000 stereo code {code} is unsupported for an aromatic source bond"),
        ));
    }
    match (order, code) {
        (_, 0) => Ok(None),
        (BondOrder::Single, 1) => Ok(Some(SourceStereoBondMarkKind::WedgeUp)),
        (BondOrder::Single, 4) => Ok(Some(SourceStereoBondMarkKind::WedgeEither)),
        (BondOrder::Single, 6) => Ok(Some(SourceStereoBondMarkKind::WedgeDown)),
        (BondOrder::Double, 3) => Ok(Some(SourceStereoBondMarkKind::DoubleBondEither)),
        _ => Err(SdfParseError::new(
            1,
            line,
            format!("V2000 stereo code {code} is unsupported for bond order {order:?}"),
        )),
    }
}

fn apply_v2000_declared_hydrogens(
    molecule: &mut Molecule,
    syntax: &V2000Syntax,
    atom_ids: &[AtomId],
    source_aromatic_bonds: &BTreeSet<BondId>,
) -> std::result::Result<(), SdfParseError> {
    for (record, atom_id) in syntax.atoms.iter().zip(atom_ids.iter().copied()) {
        let Some(declared_valence) = record.valence else {
            continue;
        };
        apply_molfile_declared_valence(
            molecule,
            atom_id,
            declared_valence,
            record.hydrogen_count,
            source_aromatic_bonds,
            MolfileVersion::V2000,
            record.line,
        )?;
    }
    Ok(())
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

fn atom_coordinates_from_v2000_line(line: &str) -> Option<[f64; 3]> {
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
                return Some([x, y, z]);
            }
        }
    }
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let point = [
        fields.first()?.parse().ok()?,
        fields.get(1)?.parse().ok()?,
        fields.get(2)?.parse().ok()?,
    ];
    point
        .iter()
        .all(|value: &f64| value.is_finite())
        .then_some(point)
}

fn apply_atom_v2000_fields(
    record: usize,
    line_number: usize,
    atom: &mut V2000AtomSyntax,
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
            atom.radical = Some(V2000RadicalSyntax::Doublet);
        }
    }
    if let Some(code) = fields.get(7).and_then(|value| value.parse::<u8>().ok()) {
        atom.hydrogen_count = match code {
            0 => None,
            1..=5 => Some(code - 1),
            _ => {
                return Err(SdfParseError::new(
                    record,
                    line_number,
                    "unsupported V2000 atom hydrogen-count code",
                ))
            }
        };
    }
    if let Some(code) = fields.get(9).and_then(|value| value.parse::<u8>().ok()) {
        atom.valence = match code {
            0 => None,
            1..=14 => Some(code),
            15 => Some(0),
            _ => {
                return Err(SdfParseError::new(
                    record,
                    line_number,
                    "unsupported V2000 atom valence code",
                ))
            }
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
    Ok(())
}

fn parse_v2000_bond_line(line: &str) -> Option<(usize, usize, u8, u8)> {
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
    Some((a, b, order_code, stereo_code.unwrap_or(0)))
}

fn ascii_field(line: &str, start: usize, end: usize) -> Option<&str> {
    std::str::from_utf8(line.as_bytes().get(start..end)?).ok()
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
                        atom.formal_charge = value;
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
                        atom.isotope = (value > 0).then_some(value);
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
                            1 => V2000RadicalSyntax::Singlet,
                            2 => V2000RadicalSyntax::Doublet,
                            3 => V2000RadicalSyntax::Triplet,
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
    F: FnMut(&mut V2000AtomSyntax, i32) -> std::result::Result<(), &'static str>,
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
            .ok_or_else(|| SdfParseError::new(record, line, "M record atom outside atom block"))?;
        apply(atom, value).map_err(|message| SdfParseError::new(record, line, message))?;
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MolWriteErrorKind {
    UnsupportedRepresentation,
    InvalidModel,
    InvalidMetadata,
    Io(std::io::ErrorKind),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolWriteError {
    pub(crate) kind: MolWriteErrorKind,
    pub(crate) message: String,
}

impl MolWriteError {
    pub(crate) fn new(message: impl Into<String>) -> Self {
        Self {
            kind: MolWriteErrorKind::UnsupportedRepresentation,
            message: message.into(),
        }
    }

    pub(crate) fn invalid_model(message: impl Into<String>) -> Self {
        Self {
            kind: MolWriteErrorKind::InvalidModel,
            message: message.into(),
        }
    }

    pub(crate) fn invalid_metadata(message: impl Into<String>) -> Self {
        Self {
            kind: MolWriteErrorKind::InvalidMetadata,
            message: message.into(),
        }
    }

    pub(crate) fn io(error: std::io::Error) -> Self {
        Self {
            kind: MolWriteErrorKind::Io(error.kind()),
            message: error.to_string(),
        }
    }

    pub const fn kind(&self) -> MolWriteErrorKind {
        self.kind
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

pub fn write_mol_v2000(molecule: &Molecule) -> std::result::Result<String, MolWriteError> {
    let record = MolfileRecord::molecule(molecule)?;
    render_mol_v2000(&record, "")
}

pub(crate) fn write_model_v2000(
    model: ModelView<'_>,
) -> std::result::Result<String, MolWriteError> {
    let record = MolfileRecord::model(model)?;
    render_mol_v2000(&record, "")
}

pub(super) fn render_mol_v2000(
    record: &MolfileRecord<'_>,
    title: &str,
) -> std::result::Result<String, MolWriteError> {
    if record.atoms.len() > V2000_MAX_ATOMS || record.bonds.len() > V2000_MAX_BONDS {
        return Err(MolWriteError::new(
            "V2000 writer supports at most 999 atoms and 999 bonds",
        ));
    }

    let program = "kekule";
    let comment = "";
    let mut out = String::new();
    out.push_str(&format!("{title}\n{program}\n{comment}\n"));
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0            999 V2000\n",
        record.atoms.len(),
        record.bonds.len()
    ));

    for record_atom in &record.atoms {
        let atom = record_atom.atom;
        let point = record_atom.position;
        if atom.atom_map.is_some_and(|atom_map| atom_map > 999) {
            return Err(MolWriteError::new(format!(
                "V2000 cannot encode atom-map number {} above 999",
                atom.atom_map.expect("atom map was checked")
            )));
        }
        if [point.x, point.y, point.z]
            .into_iter()
            .any(|coordinate| format!("{coordinate:.4}").len() > 10)
        {
            return Err(MolWriteError::new(
                "V2000 cannot encode a coordinate in its fixed-width atom field",
            ));
        }
        let valence_code =
            v2000_valence_code(record_atom.molecule, record_atom.id, record_atom.atom)?;
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

    for bond in &record.bonds {
        let order_code = v2000_bond_code(bond.order)?;
        let stereo_code = v2000_bond_stereo_code(bond.order, bond.stereo)?;
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}{:>3}  0  0  0\n",
            bond.from, bond.to, order_code, stereo_code
        ));
    }

    push_m_record(
        &mut out,
        "CHG",
        record
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(index, record_atom)| {
                (record_atom.atom.formal_charge != 0).then_some((
                    (index + 1) as u64,
                    i32::from(record_atom.atom.formal_charge),
                ))
            })
            .collect(),
    );
    push_m_record(
        &mut out,
        "ISO",
        record
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(index, record_atom)| {
                record_atom
                    .atom
                    .isotope
                    .map(|isotope| ((index + 1) as u64, i32::from(isotope)))
            })
            .collect(),
    );
    let mut radical_records = Vec::new();
    for (index, record_atom) in record.atoms.iter().enumerate() {
        let Some(radical) = record_atom.atom.radical else {
            continue;
        };
        radical_records.push(((index + 1) as u64, v2000_radical_code(radical)?));
    }
    push_m_record(&mut out, "RAD", radical_records);
    out.push_str("M  END\n");
    Ok(out)
}

pub fn write_sdf_v2000(
    records: &[SdfRecordInterpretation],
) -> std::result::Result<String, MolWriteError> {
    let mut out = String::new();
    for record in records {
        validate_sdf_title(record.title())?;
        for field in record.data_fields() {
            validate_sdf_data_field(field)?;
        }
        let structural = MolfileRecord::model(record.model().view())?;
        out.push_str(&render_mol_v2000(&structural, record.title())?);
        for field in record.data_fields() {
            out.push_str(&format!(">  <{}>\n{}\n\n", field.name(), field.value()));
        }
        out.push_str("$$$$\n");
    }
    Ok(out)
}

pub(super) fn validate_sdf_title(title: &str) -> std::result::Result<(), MolWriteError> {
    if title.contains(['\r', '\n']) {
        return Err(MolWriteError::invalid_metadata(
            "SDF record titles cannot contain line breaks",
        ));
    }
    Ok(())
}

pub(super) fn validate_sdf_data_field(
    field: &SdfDataField,
) -> std::result::Result<(), MolWriteError> {
    let name = field.name();
    if name.is_empty() || name.trim() != name || name.contains(['<', '>', '\r', '\n']) {
        return Err(MolWriteError::invalid_metadata(
            "SDF data field names must be nonempty, trimmed, and exclude angle brackets or line breaks",
        ));
    }
    let value = field.value();
    if value.contains('\r') {
        return Err(MolWriteError::invalid_metadata(
            "SDF data field values cannot contain carriage returns",
        ));
    }
    if !value.is_empty() && value.split('\n').any(str::is_empty) {
        return Err(MolWriteError::invalid_metadata(
            "SDF data field values cannot contain blank lines",
        ));
    }
    if value.lines().any(|line| line.trim() == "$$$$") {
        return Err(MolWriteError::invalid_metadata(
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
    let explicit_hydrogens = match atom.hydrogens {
        HydrogenDeclaration::Infer { explicit: 0 } => return Ok(0),
        HydrogenDeclaration::Infer { .. } => {
            return Err(MolWriteError::new(format!(
                "V2000 cannot encode represented hydrogens while leaving implicit-H inference enabled for atom {}",
                atom_id.index()
            )));
        }
        HydrogenDeclaration::Fixed(explicit) => explicit,
    };
    let valence = explicit_valence(mol, atom_id) + usize::from(explicit_hydrogens);
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
        BondOrder::Dative => Ok(9),
        BondOrder::Quadruple => Err(MolWriteError::new(
            "V2000 writer does not support quadruple bonds",
        )),
    }
}

fn v2000_bond_stereo_code(
    order: BondOrder,
    stereo: Option<SourceStereoBondMarkKind>,
) -> std::result::Result<u8, MolWriteError> {
    match (order, stereo) {
        (_, None) => Ok(0),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeUp)) => Ok(1),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeEither)) => Ok(4),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeDown)) => Ok(6),
        (BondOrder::Double, Some(SourceStereoBondMarkKind::DoubleBondEither)) => Ok(3),
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

fn push_m_record(out: &mut String, code: &str, pairs: Vec<(u64, i32)>) {
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
