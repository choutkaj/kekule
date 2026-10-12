use std::collections::BTreeSet;

use crate::diagnostics::DiagnosticKind;
use crate::geometry::Point3;

use super::molfile_document::{
    MolfileAtom, MolfileBond, MolfileBondOrder, MolfileBondStereo, MolfileCtab, MolfileLines,
    MolfileParseError, MolfileRadical, MolfileUnsupported, MolfileUnsupportedContent,
};

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
