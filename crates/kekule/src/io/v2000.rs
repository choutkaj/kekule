use std::collections::BTreeSet;
use std::fmt;

use crate::algorithms::explicit_valence;
use crate::chemistry::SourceStereoBondMarkKind;
use crate::core::*;
use crate::diagnostics::DiagnosticKind;
use crate::geometry::Point3;

use super::molfile_document::{
    MolfileAtom, MolfileBond, MolfileBondOrder, MolfileBondStereo, MolfileCtab, MolfileLines,
    MolfileParseError, MolfileRadical, MolfileUnsupported, MolfileUnsupportedContent,
};
use super::molfile_write::MolfileRecord;
use super::sdf_document::SdfDataField;

pub(super) const V2000_MAX_ATOMS: usize = 999;
pub(super) const V2000_MAX_BONDS: usize = 999;

/// Parses the V2000 connection table after the header and returns the index
/// of its `M  END` line.
pub(super) fn parse_v2000_ctab(
    lines: MolfileLines<'_>,
) -> std::result::Result<(MolfileCtab, usize), MolfileParseError> {
    let counts = lines.text(3).expect("the caller checked the counts line");
    let (atom_count, bond_count) = parse_counts_line(counts)
        .ok_or_else(|| lines.error(3, DiagnosticKind::Syntax, "invalid V2000 counts line"))?;
    if atom_count > V2000_MAX_ATOMS || bond_count > V2000_MAX_BONDS {
        return Err(lines.error(
            3,
            DiagnosticKind::InvalidValue,
            "V2000 counts exceed the supported 999 atom or bond limit",
        ));
    }
    let bond_start = 4 + atom_count;
    let property_start = bond_start + bond_count;
    if lines.len() < property_start {
        return Err(lines.error(
            lines.len() - 1,
            DiagnosticKind::Syntax,
            "Molfile ended before declared atom and bond blocks",
        ));
    }

    let mut ctab = MolfileCtab {
        chiral_flag: v2000_chiral_flag(counts),
        ..MolfileCtab::default()
    };
    for index in 4..bond_start {
        let line = lines.text(index).expect("block length checked");
        ctab.atoms.push(parse_v2000_atom(lines, index, line)?);
    }
    let mut endpoints = BTreeSet::new();
    for index in bond_start..property_start {
        let line = lines.text(index).expect("block length checked");
        let bond = parse_v2000_bond(lines, index, line, ctab.atoms.len())?;
        let [a, b] = bond.atoms;
        if !endpoints.insert((a.min(b), a.max(b))) {
            return Err(lines.error(
                index,
                DiagnosticKind::Inconsistent,
                "duplicate bond endpoints",
            ));
        }
        ctab.bonds.push(bond);
    }
    let end = (property_start..lines.len())
        .find(|&index| {
            lines
                .text(index)
                .is_some_and(|line| line.trim() == "M  END")
        })
        .ok_or_else(|| {
            lines.error(
                lines.len().saturating_sub(1),
                DiagnosticKind::Syntax,
                "missing M  END",
            )
        })?;
    parse_v2000_properties(lines, property_start..end, &mut ctab)?;
    Ok((ctab, end))
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

fn v2000_chiral_flag(line: &str) -> bool {
    match ascii_field(line, 12, 15) {
        Some(field) if line.is_ascii() => field.trim() == "1",
        _ => line.split_whitespace().nth(4) == Some("1"),
    }
}

fn parse_v2000_atom(
    lines: MolfileLines<'_>,
    index: usize,
    line: &str,
) -> std::result::Result<MolfileAtom, MolfileParseError> {
    let symbol = atom_symbol_from_v2000_line(line)
        .ok_or_else(|| lines.error(index, DiagnosticKind::Syntax, "invalid V2000 atom line"))?;
    let [x, y, z] = atom_coordinates_from_v2000_line(line).ok_or_else(|| {
        lines.error(
            index,
            DiagnosticKind::Syntax,
            "invalid V2000 atom coordinates",
        )
    })?;
    let mut atom = MolfileAtom::new(symbol, Point3::new(x, y, z), lines.location(index));
    apply_atom_v2000_fields(lines, index, &mut atom, line)?;
    Ok(atom)
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
    if let Some(point) = fixed_atom_coordinates(line) {
        return Some(point);
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

fn fixed_atom_coordinates(line: &str) -> Option<[f64; 3]> {
    if let (Some(x), Some(y), Some(z)) = (
        ascii_field(line, 0, 10),
        ascii_field(line, 10, 20),
        // Some supplied CTfiles use the separator column for the last
        // decimal digit of a large negative z coordinate. Read it rather
        // than silently truncating an otherwise valid number.
        ascii_field(line, 20, 31),
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
    None
}

fn apply_atom_v2000_fields(
    lines: MolfileLines<'_>,
    index: usize,
    atom: &mut MolfileAtom,
    line: &str,
) -> std::result::Result<(), MolfileParseError> {
    // Coordinate fields can touch, for example a negative y immediately
    // followed by a full-width negative z. Tokenize only the atom fields so
    // that this cannot shift charge, hydrogen, valence or atom-map columns.
    let fields = if fixed_atom_coordinates(line).is_some() {
        line[31..].split_whitespace().collect::<Vec<_>>()
    } else {
        line.split_whitespace().skip(3).collect::<Vec<_>>()
    };
    if let Some(value) = fields.get(2) {
        let charge_code = value.parse::<u8>().map_err(|_| {
            lines.error(
                index,
                DiagnosticKind::Syntax,
                "invalid V2000 atom charge code",
            )
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
                return Err(lines.error(
                    index,
                    DiagnosticKind::InvalidValue,
                    format!("invalid V2000 atom charge code {charge_code}"),
                ))
            }
        };
        if charge_code == 4 {
            atom.radical = Some(MolfileRadical::Doublet);
        }
    }
    if let Some(code) = fields.get(4).and_then(|value| value.parse::<u8>().ok()) {
        atom.hydrogen_count = match code {
            0 => None,
            1..=5 => Some(code - 1),
            _ => {
                return Err(lines.error(
                    index,
                    DiagnosticKind::InvalidValue,
                    format!("invalid V2000 atom hydrogen-count code {code}"),
                ))
            }
        };
    }
    if let Some(code) = fields.get(6).and_then(|value| value.parse::<u8>().ok()) {
        atom.valence = match code {
            0 => None,
            1..=14 => Some(code),
            15 => Some(0),
            _ => {
                return Err(lines.error(
                    index,
                    DiagnosticKind::InvalidValue,
                    format!("invalid V2000 atom valence code {code}"),
                ))
            }
        };
    }
    if let Some(atom_map) = fields
        .get(10)
        .or_else(|| fields.get(9))
        .and_then(|value| value.parse::<u32>().ok())
    {
        if atom_map != 0 {
            atom.atom_map = Some(atom_map);
        }
    }
    Ok(())
}

fn parse_v2000_bond(
    lines: MolfileLines<'_>,
    index: usize,
    line: &str,
    atom_count: usize,
) -> std::result::Result<MolfileBond, MolfileParseError> {
    let (a, b, order_code, stereo_code) = parse_v2000_bond_line(line)
        .ok_or_else(|| lines.error(index, DiagnosticKind::Syntax, "invalid V2000 bond line"))?;
    let endpoint = |serial: usize| {
        serial
            .checked_sub(1)
            .filter(|&atom| atom < atom_count)
            .ok_or_else(|| {
                lines.error(
                    index,
                    DiagnosticKind::InvalidReference,
                    "bond endpoint outside atom block",
                )
            })
    };
    let atoms = [endpoint(a)?, endpoint(b)?];
    if atoms[0] == atoms[1] {
        return Err(lines.error(
            index,
            DiagnosticKind::InvalidValue,
            "bond endpoints must be distinct",
        ));
    }
    let order = match order_code {
        0 => MolfileBondOrder::Zero,
        1 => MolfileBondOrder::Single,
        2 => MolfileBondOrder::Double,
        3 => MolfileBondOrder::Triple,
        4 => MolfileBondOrder::Aromatic,
        9 => MolfileBondOrder::Dative,
        5..=8 | 10 => {
            return Err(lines.error(
                index,
                DiagnosticKind::Unsupported,
                format!("unsupported V2000 bond order code {order_code}"),
            ))
        }
        _ => {
            return Err(lines.error(
                index,
                DiagnosticKind::InvalidValue,
                format!("invalid V2000 bond order code {order_code}"),
            ))
        }
    };
    let stereo = match (order, stereo_code) {
        (_, 0) => None,
        (MolfileBondOrder::Aromatic, _) => {
            return Err(lines.error(
                index,
                DiagnosticKind::Unsupported,
                format!("V2000 stereo code {stereo_code} is unsupported for an aromatic bond"),
            ))
        }
        (MolfileBondOrder::Single, 1) => Some(MolfileBondStereo::Up),
        (MolfileBondOrder::Single, 4) => Some(MolfileBondStereo::Either),
        (MolfileBondOrder::Single, 6) => Some(MolfileBondStereo::Down),
        (MolfileBondOrder::Double, 3) => Some(MolfileBondStereo::Either),
        _ => {
            return Err(lines.error(
                index,
                DiagnosticKind::InvalidValue,
                format!("V2000 stereo code {stereo_code} is invalid for bond order {order:?}"),
            ))
        }
    };
    Ok(MolfileBond {
        atoms,
        order,
        stereo,
        location: lines.location(index),
    })
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

/// Reads the properties block: `M  CHG`, `M  ISO`, and `M  RAD` update atoms;
/// every other record is classified as unsupported content.
fn parse_v2000_properties(
    lines: MolfileLines<'_>,
    range: std::ops::Range<usize>,
    ctab: &mut MolfileCtab,
) -> std::result::Result<(), MolfileParseError> {
    let mut replaced_legacy_charge_and_radical = false;
    let mut index = range.start;
    while index < range.end {
        let start = index;
        let line = lines
            .text(index)
            .expect("properties lie inside the Molfile");
        index += 1;
        if line.trim().is_empty() {
            continue;
        }
        let Some(record) = line.strip_prefix("M  ") else {
            let keyword = line
                .split_whitespace()
                .next()
                .expect("nonblank line has a field");
            // An atom alias or group abbreviation continues on the next line.
            if matches!(keyword, "A" | "G") && index < range.end {
                index += 1;
            }
            ctab.unsupported.push(MolfileUnsupported {
                content: MolfileUnsupportedContent::Record(keyword.to_owned()),
                location: lines.location(start),
            });
            continue;
        };
        let fields = record.split_whitespace().collect::<Vec<_>>();
        let code = fields.first().copied().unwrap_or_default();
        if matches!(code, "CHG" | "RAD") && !replaced_legacy_charge_and_radical {
            // CTfile properties replace the entire legacy charge/radical block.
            // Clear once: later property lines accumulate independent entries.
            for atom in &mut ctab.atoms {
                atom.formal_charge = 0;
                atom.radical = None;
            }
            replaced_legacy_charge_and_radical = true;
        }
        let pairs = &fields[fields.len().min(1)..];
        match code {
            "CHG" => {
                parse_atom_value_pairs(lines, start, pairs, &mut ctab.atoms, |atom, value| {
                    atom.formal_charge = value;
                    Ok(())
                })?
            }
            "ISO" => {
                parse_atom_value_pairs(lines, start, pairs, &mut ctab.atoms, |atom, value| {
                    atom.isotope = match value {
                        0 => None,
                        value if value > 0 => Some(
                            u16::try_from(value).map_err(|_| "M  ISO mass number exceeds u16")?,
                        ),
                        _ => return Err("M  ISO mass number must not be negative"),
                    };
                    Ok(())
                })?
            }
            "RAD" => {
                parse_atom_value_pairs(lines, start, pairs, &mut ctab.atoms, |atom, value| {
                    atom.radical = match value {
                        0 => None,
                        1 => Some(MolfileRadical::Singlet),
                        2 => Some(MolfileRadical::Doublet),
                        3 => Some(MolfileRadical::Triplet),
                        _ => return Err("invalid M  RAD code"),
                    };
                    Ok(())
                })?
            }
            _ => ctab.unsupported.push(MolfileUnsupported {
                content: MolfileUnsupportedContent::PropertyRecord(code.to_owned()),
                location: lines.location(start),
            }),
        }
    }
    Ok(())
}

fn parse_atom_value_pairs<F>(
    lines: MolfileLines<'_>,
    index: usize,
    fields: &[&str],
    atoms: &mut [MolfileAtom],
    mut apply: F,
) -> std::result::Result<(), MolfileParseError>
where
    F: FnMut(&mut MolfileAtom, i32) -> std::result::Result<(), &'static str>,
{
    let syntax = |message: &str| lines.error(index, DiagnosticKind::Syntax, message);
    let (count, pairs) = fields
        .split_first()
        .ok_or_else(|| syntax("invalid M record count"))?;
    let count = count
        .parse::<usize>()
        .map_err(|_| syntax("invalid M record count"))?;
    if count.checked_mul(2) != Some(pairs.len()) {
        return Err(syntax("M record pair count does not match its fields"));
    }
    for pair in pairs.chunks(2) {
        let serial = pair[0]
            .parse::<usize>()
            .map_err(|_| syntax("invalid M record atom index"))?;
        let value = pair[1]
            .parse::<i32>()
            .map_err(|_| syntax("invalid M record value"))?;
        let atom = serial
            .checked_sub(1)
            .and_then(|atom| atoms.get_mut(atom))
            .ok_or_else(|| {
                lines.error(
                    index,
                    DiagnosticKind::InvalidReference,
                    "M record atom outside atom block",
                )
            })?;
        apply(atom, value)
            .map_err(|message| lines.error(index, DiagnosticKind::InvalidValue, message))?;
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MolWriteErrorKind {
    UnsupportedRepresentation,
    /// A bounded writer calculation could not complete within its work limit.
    ResourceLimit,
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
    pub(crate) fn resource_limit(message: impl Into<String>) -> Self {
        Self {
            kind: MolWriteErrorKind::ResourceLimit,
            message: message.into(),
        }
    }

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

pub(super) fn render_mol_v2000(
    record: &MolfileRecord<'_>,
    title: &str,
) -> std::result::Result<String, MolWriteError> {
    if !record.stereo_groups.is_empty() {
        return Err(MolWriteError::new(
            "V2000 cannot encode enhanced stereo groups; use V3000",
        ));
    }
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

pub(super) fn validate_sdf_title(title: &str) -> std::result::Result<(), MolWriteError> {
    if title.contains(['\r', '\n']) {
        return Err(MolWriteError::invalid_metadata(
            "SDF record titles cannot contain line breaks",
        ));
    }
    if title.trim() == "$$$$" {
        return Err(MolWriteError::invalid_metadata(
            "SDF record titles cannot be a record delimiter",
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
    let fixed_hydrogens = match atom.hydrogens {
        ImplicitHydrogens::Inferred => return Ok(0),
        ImplicitHydrogens::Fixed(explicit) => explicit,
    };
    let valence = explicit_valence(mol, atom_id) + usize::from(fixed_hydrogens);
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
        (BondOrder::Double, Some(SourceStereoBondMarkKind::MolfileDoubleBondGeometry)) => Ok(0),
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
    match (radical.electron_count(), radical.spin_multiplicity()) {
        (2, Some(1)) => Ok(1),
        (1, Some(2)) => Ok(2),
        (2, Some(3)) => Ok(3),
        _ => Err(MolWriteError::new(
            "V2000 writer requires an encodable radical electron count and explicit spin multiplicity",
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
