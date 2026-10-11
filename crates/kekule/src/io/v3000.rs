use std::collections::{BTreeMap, BTreeSet};

use crate::algorithms::explicit_valence;
use crate::chemistry::SourceStereoBondMarkKind;
use crate::core::*;
use crate::diagnostics::DiagnosticKind;
use crate::geometry::Point3;
use crate::io::{MoleculeWriteError, MolfileParseOptions};

use super::molfile_document::{
    MolfileAtom, MolfileAtomParity, MolfileBond, MolfileBondOrder, MolfileBondStereo, MolfileCtab,
    MolfileLines, MolfileParseError, MolfileRadical, MolfileStereoGroup, MolfileStereoGroupKind,
    MolfileUnsupported, MolfileUnsupportedContent,
};
use super::molfile_write::MolfileRecord;

pub(super) fn render_mol_v3000(
    record: &MolfileRecord<'_>,
    title: &str,
) -> std::result::Result<String, MoleculeWriteError> {
    let program = "kekule";
    let comment = "";

    let mut out = String::new();
    out.push_str(&format!("{title}\n{program}\n{comment}\n"));
    out.push_str("  0  0  0  0  0  0            999 V3000\n");
    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {} 0 0 0\n",
        record.atoms.len(),
        record.bonds.len()
    ));
    out.push_str("M  V30 BEGIN ATOM\n");
    for (index, record_atom) in (1u64..).zip(record.atoms.iter()) {
        let atom = record_atom.atom;
        let point = record_atom.position;
        out.push_str(&format!(
            "M  V30 {index} {} {} {} {} {}",
            atom.element.symbol(),
            point.x,
            point.y,
            point.z,
            atom.atom_map.unwrap_or(0)
        ));
        if atom.formal_charge != 0 {
            out.push_str(&format!(" CHG={}", atom.formal_charge));
        }
        if let Some(isotope) = atom.isotope {
            out.push_str(&format!(" MASS={isotope}"));
        }
        if let Some(radical) = atom.radical {
            out.push_str(&format!(" RAD={}", v3000_radical_code(radical)?));
        }
        match atom.hydrogens {
            ImplicitHydrogens::Inferred => {}
            ImplicitHydrogens::Fixed(explicit) => {
                // HCOUNT is a query constraint in CTfile, not a molecular
                // hydrogen declaration. VAL preserves the fixed total valence
                // without turning the emitted atom into a query atom.
                let valence =
                    explicit_valence(record_atom.molecule, record_atom.id) + usize::from(explicit);
                // VAL=0 means unspecified, whereas -1 explicitly means zero.
                if valence == 0 {
                    out.push_str(" VAL=-1");
                } else {
                    out.push_str(&format!(" VAL={valence}"));
                }
            }
        }
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
    out.push_str("M  V30 BEGIN BOND\n");
    for (serial, bond) in (1u64..).zip(record.bonds.iter()) {
        let order_code = v3000_bond_code(bond.order)?;
        out.push_str(&format!(
            "M  V30 {serial} {order_code} {} {}",
            bond.from, bond.to
        ));
        if let Some(cfg) = v3000_bond_cfg(bond.order, bond.stereo)? {
            out.push_str(&format!(" CFG={cfg}"));
        }
        out.push('\n');
    }
    out.push_str("M  V30 END BOND\n");
    if !record.stereo_groups.is_empty() {
        out.push_str("M  V30 BEGIN COLLECTION\n");
        let mut absolute_written = false;
        for (index, group) in record.stereo_groups.iter().enumerate() {
            if group.kind == StereoGroupKind::Absolute && absolute_written {
                continue;
            }
            let atoms = if group.kind == StereoGroupKind::Absolute {
                absolute_written = true;
                record
                    .stereo_groups
                    .iter()
                    .filter(|group| group.kind == StereoGroupKind::Absolute)
                    .flat_map(|group| group.atoms.iter().copied())
                    .collect::<Vec<_>>()
            } else {
                group.atoms.clone()
            };
            let kind = match group.kind {
                StereoGroupKind::Absolute => "MDLV30/STEABS".to_owned(),
                StereoGroupKind::Or => format!("MDLV30/STEREL{}", index + 1),
                StereoGroupKind::And => format!("MDLV30/STERAC{}", index + 1),
                _ => {
                    return Err(MoleculeWriteError::new(
                        "V3000 cannot encode this stereo group kind",
                    ))
                }
            };
            out.push_str(&format!("M  V30 {kind}"));
            out.push_str(&format!(" ATOMS=({}", atoms.len()));
            for member in atoms {
                out.push_str(&format!(" {member}"));
            }
            out.push(')');
            out.push('\n');
        }
        out.push_str("M  V30 END COLLECTION\n");
    }
    out.push_str("M  V30 END CTAB\n");
    out.push_str("M  END\n");
    Ok(out)
}

/// Parses the V3000 connection table after the header and returns the index
/// of its `M  END` line.
pub(super) fn parse_v3000_ctab(
    lines: MolfileLines<'_>,
    options: MolfileParseOptions,
) -> std::result::Result<(MolfileCtab, usize), MolfileParseError> {
    let (records, end) = collect_v3000_lines(lines, options)?;
    let Some(first) = records.first() else {
        return Err(lines.error(end, DiagnosticKind::Syntax, "missing V3000 BEGIN CTAB"));
    };
    if first.body != "BEGIN CTAB" {
        return Err(lines.error(
            first.index,
            DiagnosticKind::Syntax,
            "V3000 records must be inside the CTAB",
        ));
    }
    let Some(close) = records.iter().position(|record| record.body == "END CTAB") else {
        return Err(lines.error(end, DiagnosticKind::Syntax, "missing V3000 END CTAB"));
    };
    if let Some(outside) = records.get(close + 1) {
        return Err(lines.error(
            outside.index,
            DiagnosticKind::Syntax,
            "V3000 records must be inside the CTAB",
        ));
    }

    let body = &records[1..close];
    let mut ctab = MolfileCtab::default();
    let mut counts = None::<(V3000Counts, usize)>;
    let mut atom_positions = BTreeMap::<usize, usize>::new();
    let mut atoms_read = false;
    let mut bonds_read = false;
    let mut collection_read = false;
    let mut position = 0;
    while position < body.len() {
        let record = &body[position];
        let syntax = |message: &str| lines.error(record.index, DiagnosticKind::Syntax, message);
        position += 1;
        if let Some(name) = record.body.strip_prefix("BEGIN ") {
            let name = name.trim();
            let close = format!("END {name}");
            let length = body[position..]
                .iter()
                .position(|row| row.body == close)
                .ok_or_else(|| syntax(&format!("missing V3000 END {name}")))?;
            let rows = &body[position..position + length];
            position += length + 1;
            match name {
                "ATOM" => {
                    if atoms_read {
                        return Err(syntax("V3000 must contain exactly one ATOM section"));
                    }
                    let Some((counts, counts_index)) = counts else {
                        return Err(syntax(
                            "V3000 COUNTS line must occur before the ATOM section",
                        ));
                    };
                    if rows.len() != counts.atoms {
                        return Err(lines.error(
                            counts_index,
                            DiagnosticKind::Inconsistent,
                            "V3000 COUNTS do not match the ATOM section size",
                        ));
                    }
                    for row in rows {
                        let (serial, atom) =
                            parse_v3000_atom(lines, row, ctab.atoms.len(), &mut ctab.unsupported)?;
                        if atom_positions.insert(serial, ctab.atoms.len()).is_some() {
                            return Err(lines.error(
                                row.index,
                                DiagnosticKind::Inconsistent,
                                "duplicate V3000 atom index",
                            ));
                        }
                        ctab.atoms.push(atom);
                    }
                    atoms_read = true;
                }
                "BOND" => {
                    let (Some((counts, counts_index)), true) = (counts, atoms_read) else {
                        return Err(syntax("V3000 BOND section must follow the ATOM section"));
                    };
                    if bonds_read {
                        return Err(syntax("V3000 must contain exactly one BOND section"));
                    }
                    if collection_read {
                        return Err(syntax(
                            "V3000 BOND section must precede the COLLECTION section",
                        ));
                    }
                    if rows.len() != counts.bonds {
                        return Err(lines.error(
                            counts_index,
                            DiagnosticKind::Inconsistent,
                            "V3000 COUNTS do not match the BOND section size",
                        ));
                    }
                    parse_v3000_bonds(lines, rows, &atom_positions, &mut ctab)?;
                    bonds_read = true;
                }
                "COLLECTION" => {
                    if !bonds_read && !bond_block_optional(counts, atoms_read) {
                        return Err(syntax(
                            "V3000 COLLECTION section must follow the BOND section",
                        ));
                    }
                    if collection_read {
                        return Err(syntax("V3000 permits one COLLECTION section"));
                    }
                    parse_v3000_collection(lines, rows, &atom_positions, &mut ctab)?;
                    collection_read = true;
                }
                "CTAB" => return Err(syntax("V3000 must contain exactly one CTAB")),
                _ => ctab.unsupported.push(MolfileUnsupported {
                    content: MolfileUnsupportedContent::Block(name.to_owned()),
                    location: lines.location(record.index),
                }),
            }
            continue;
        }
        let Some(keyword) = record.body.split_whitespace().next() else {
            continue;
        };
        if keyword == "COUNTS" {
            if counts.is_some() {
                return Err(syntax("V3000 CTAB must contain exactly one COUNTS line"));
            }
            let parsed = parse_v3000_counts(&record.body)
                .ok_or_else(|| syntax("invalid V3000 COUNTS line"))?;
            if parsed.atoms > options.max_v3000_atoms {
                return Err(lines.error(
                    record.index,
                    DiagnosticKind::Limit,
                    "V3000 atom count exceeds configured limit",
                ));
            }
            if parsed.bonds > options.max_v3000_bonds {
                return Err(lines.error(
                    record.index,
                    DiagnosticKind::Limit,
                    "V3000 bond count exceeds configured limit",
                ));
            }
            ctab.chiral_flag = parsed.chiral;
            counts = Some((parsed, record.index));
        } else if keyword == "END" {
            return Err(syntax(&format!("unexpected V3000 {}", record.body)));
        } else if is_v3000_stereo_group(keyword) {
            return Err(syntax("V3000 stereo group is outside COLLECTION"));
        } else {
            ctab.unsupported.push(MolfileUnsupported {
                content: MolfileUnsupportedContent::Record(keyword.to_owned()),
                location: lines.location(record.index),
            });
        }
    }
    if !atoms_read {
        return Err(lines.error(end, DiagnosticKind::Syntax, "missing V3000 ATOM section"));
    }
    if !bonds_read && !bond_block_optional(counts, atoms_read) {
        return Err(lines.error(end, DiagnosticKind::Syntax, "missing V3000 BOND section"));
    }
    Ok((ctab, end))
}

/// Writers such as RDKit omit the BOND section of a CTAB without bonds.
fn bond_block_optional(counts: Option<(V3000Counts, usize)>, atoms_read: bool) -> bool {
    atoms_read && counts.is_some_and(|(counts, _)| counts.bonds == 0)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct V3000Counts {
    atoms: usize,
    bonds: usize,
    chiral: bool,
}

/// One `M  V30` logical line after joining continuations, located at its
/// first physical line.
#[derive(Debug, Clone, PartialEq)]
struct V3000Line {
    index: usize,
    body: String,
}

fn collect_v3000_lines(
    lines: MolfileLines<'_>,
    options: MolfileParseOptions,
) -> std::result::Result<(Vec<V3000Line>, usize), MolfileParseError> {
    let limit = |index| {
        lines.error(
            index,
            DiagnosticKind::Limit,
            "V3000 logical line exceeds configured byte limit",
        )
    };
    let mut records = Vec::new();
    let mut index = 4usize;
    while let Some(line) = lines.text(index) {
        let start = index;
        if line.trim() == "M  END" {
            return Ok((records, index));
        }
        let body = v3000_body(line)
            .ok_or_else(|| lines.error(index, DiagnosticKind::Syntax, "expected M  V30 record"))?;
        if body.len() > options.max_v3000_logical_line_bytes {
            return Err(limit(start));
        }
        let mut body = body.to_owned();
        while body.ends_with('-') {
            body.pop();
            index += 1;
            let continuation_line = lines.text(index).ok_or_else(|| {
                lines.error(
                    start,
                    DiagnosticKind::Syntax,
                    "unterminated V3000 continuation",
                )
            })?;
            let continuation = v3000_body(continuation_line)
                .ok_or_else(|| {
                    lines.error(index, DiagnosticKind::Syntax, "invalid V3000 continuation")
                })?
                .trim_start();
            if body.len().saturating_add(continuation.len()) > options.max_v3000_logical_line_bytes
            {
                return Err(limit(start));
            }
            body.push_str(continuation);
        }
        records.push(V3000Line { index: start, body });
        index += 1;
    }
    Err(lines.error(
        lines.len().saturating_sub(1),
        DiagnosticKind::Syntax,
        "missing M  END",
    ))
}

fn v3000_body(line: &str) -> Option<&str> {
    let trimmed = line.strip_prefix("M  V30 ")?;
    Some(trimmed.trim())
}

/// Splits a V3000 record at whitespace outside parentheses and double quotes.
fn v3000_fields(body: &str) -> Option<Vec<&str>> {
    let mut fields = Vec::new();
    let mut depth = 0usize;
    let mut quoted = false;
    let mut start = None;
    for (offset, character) in body.char_indices() {
        match character {
            '"' => quoted = !quoted,
            '(' if !quoted => depth += 1,
            ')' if !quoted => depth = depth.checked_sub(1)?,
            _ if character.is_whitespace() && !quoted && depth == 0 => {
                if let Some(start) = start.take() {
                    fields.push(&body[start..offset]);
                }
                continue;
            }
            _ => {}
        }
        start.get_or_insert(offset);
    }
    if quoted || depth != 0 {
        return None;
    }
    fields.extend(start.map(|start| &body[start..]));
    Some(fields)
}

fn parse_v3000_counts(line: &str) -> Option<V3000Counts> {
    let mut fields = line.split_whitespace();
    (fields.next()? == "COUNTS").then_some(())?;
    let atoms = fields.next()?.parse().ok()?;
    let bonds = fields.next()?.parse().ok()?;
    let _sgroups = fields.next()?.parse::<usize>().ok()?;
    let _three_dimensional_constraints = fields.next()?.parse::<usize>().ok()?;
    let chiral = match fields.next()? {
        "0" => false,
        "1" => true,
        _ => return None,
    };
    if let Some(regno) = fields.next() {
        regno
            .strip_prefix("REGNO=")?
            .parse::<u64>()
            .ok()
            .map(|_| ())?;
    }
    fields.next().is_none().then_some(V3000Counts {
        atoms,
        bonds,
        chiral,
    })
}

fn parse_v3000_atom(
    lines: MolfileLines<'_>,
    row: &V3000Line,
    position: usize,
    unsupported: &mut Vec<MolfileUnsupported>,
) -> std::result::Result<(usize, MolfileAtom), MolfileParseError> {
    let syntax = || lines.error(row.index, DiagnosticKind::Syntax, "invalid V3000 atom line");
    let invalid = |name: &str| {
        lines.error(
            row.index,
            DiagnosticKind::InvalidValue,
            format!("invalid V3000 {name} value"),
        )
    };
    let fields = v3000_fields(&row.body).ok_or_else(syntax)?;
    let [serial, symbol, x, y, z, atom_map, options @ ..] = fields.as_slice() else {
        return Err(syntax());
    };
    let serial = serial.parse::<usize>().map_err(|_| syntax())?;
    if serial == 0 {
        return Err(lines.error(
            row.index,
            DiagnosticKind::InvalidValue,
            "V3000 atom indices must be positive",
        ));
    }
    let coordinate = |value: &str| parse_finite_f64(value).ok_or_else(syntax);
    let position_value = Point3::new(coordinate(x)?, coordinate(y)?, coordinate(z)?);
    let mut atom = MolfileAtom::new(symbol, position_value, lines.location(row.index));
    let atom_map = atom_map.parse::<u32>().map_err(|_| syntax())?;
    atom.atom_map = (atom_map != 0).then_some(atom_map);
    let mut seen = BTreeSet::new();
    for field in options {
        let (key, value) = split_v3000_option(field).ok_or_else(syntax)?;
        if !seen.insert(key) {
            return Err(lines.error(
                row.index,
                DiagnosticKind::Syntax,
                format!("duplicate V3000 atom option `{key}`"),
            ));
        }
        match key {
            "CFG" => {
                atom.parity = match value {
                    "0" => None,
                    "1" => Some(MolfileAtomParity::Odd),
                    "2" => Some(MolfileAtomParity::Even),
                    "3" => Some(MolfileAtomParity::Either),
                    _ => return Err(invalid("CFG")),
                };
            }
            "CHG" => atom.formal_charge = value.parse().map_err(|_| invalid("CHG"))?,
            "MASS" => {
                atom.isotope = match value.parse::<i32>().map_err(|_| invalid("MASS"))? {
                    0 => None,
                    mass => Some(u16::try_from(mass).map_err(|_| invalid("MASS"))?),
                };
            }
            "RAD" => {
                atom.radical = match value {
                    "0" => None,
                    "1" => Some(MolfileRadical::Singlet),
                    "2" => Some(MolfileRadical::Doublet),
                    "3" => Some(MolfileRadical::Triplet),
                    _ => return Err(invalid("RAD")),
                };
            }
            "HCOUNT" => {
                atom.hydrogen_count = parse_v3000_count(value).ok_or_else(|| invalid("HCOUNT"))?;
            }
            "VAL" => atom.valence = parse_v3000_count(value).ok_or_else(|| invalid("VAL"))?,
            _ => unsupported.push(MolfileUnsupported {
                content: MolfileUnsupportedContent::AtomOption {
                    atom: position,
                    key: key.to_owned(),
                },
                location: lines.location(row.index),
            }),
        }
    }
    Ok((serial, atom))
}

/// Decodes `HCOUNT` and `VAL`: 0 is unspecified and -1 is an explicit zero.
fn parse_v3000_count(value: &str) -> Option<Option<u8>> {
    match value.parse::<i32>().ok()? {
        -1 => Some(Some(0)),
        0 => Some(None),
        count => u8::try_from(count).ok().map(Some),
    }
}

fn parse_v3000_bonds(
    lines: MolfileLines<'_>,
    rows: &[V3000Line],
    atom_positions: &BTreeMap<usize, usize>,
    ctab: &mut MolfileCtab,
) -> std::result::Result<(), MolfileParseError> {
    let mut serials = BTreeSet::new();
    let mut endpoints = BTreeSet::new();
    for row in rows {
        let (serial, bond) = parse_v3000_bond(
            lines,
            row,
            ctab.bonds.len(),
            atom_positions,
            &mut ctab.unsupported,
        )?;
        let inconsistent =
            |message: &str| lines.error(row.index, DiagnosticKind::Inconsistent, message);
        if !serials.insert(serial) {
            return Err(inconsistent("duplicate V3000 bond index"));
        }
        let [a, b] = bond.atoms;
        if !endpoints.insert((a.min(b), a.max(b))) {
            return Err(inconsistent("duplicate bond endpoints"));
        }
        ctab.bonds.push(bond);
    }
    Ok(())
}

fn parse_v3000_bond(
    lines: MolfileLines<'_>,
    row: &V3000Line,
    position: usize,
    atom_positions: &BTreeMap<usize, usize>,
    unsupported: &mut Vec<MolfileUnsupported>,
) -> std::result::Result<(usize, MolfileBond), MolfileParseError> {
    let syntax = || lines.error(row.index, DiagnosticKind::Syntax, "invalid V3000 bond line");
    let fields = v3000_fields(&row.body).ok_or_else(syntax)?;
    let [serial, order, a, b, options @ ..] = fields.as_slice() else {
        return Err(syntax());
    };
    let serial = serial.parse::<usize>().map_err(|_| syntax())?;
    if serial == 0 {
        return Err(lines.error(
            row.index,
            DiagnosticKind::InvalidValue,
            "V3000 bond indices must be positive",
        ));
    }
    let order_code = order.parse::<u8>().map_err(|_| syntax())?;
    let order = match order_code {
        0 => MolfileBondOrder::Zero,
        1 => MolfileBondOrder::Single,
        2 => MolfileBondOrder::Double,
        3 => MolfileBondOrder::Triple,
        4 => MolfileBondOrder::Aromatic,
        9 => MolfileBondOrder::Dative,
        5..=8 | 10 => {
            return Err(lines.error(
                row.index,
                DiagnosticKind::Unsupported,
                format!("unsupported V3000 bond order code {order_code}"),
            ))
        }
        _ => {
            return Err(lines.error(
                row.index,
                DiagnosticKind::InvalidValue,
                format!("invalid V3000 bond order code {order_code}"),
            ))
        }
    };
    let endpoint = |serial: &str| {
        let serial = serial.parse::<usize>().map_err(|_| syntax())?;
        atom_positions.get(&serial).copied().ok_or_else(|| {
            lines.error(
                row.index,
                DiagnosticKind::InvalidReference,
                "bond endpoint outside atom block",
            )
        })
    };
    let atoms = [endpoint(a)?, endpoint(b)?];
    if atoms[0] == atoms[1] {
        return Err(lines.error(
            row.index,
            DiagnosticKind::InvalidValue,
            "bond endpoints must be distinct",
        ));
    }
    let mut stereo = None;
    let mut seen = BTreeSet::new();
    for field in options {
        let (key, value) = split_v3000_option(field).ok_or_else(syntax)?;
        if !seen.insert(key) {
            return Err(lines.error(
                row.index,
                DiagnosticKind::Syntax,
                format!("duplicate V3000 bond option `{key}`"),
            ));
        }
        if key != "CFG" {
            unsupported.push(MolfileUnsupported {
                content: MolfileUnsupportedContent::BondOption {
                    bond: position,
                    key: key.to_owned(),
                },
                location: lines.location(row.index),
            });
            continue;
        }
        let cfg = value.parse::<u8>().map_err(|_| syntax())?;
        stereo = match (order, cfg) {
            (_, 0) => None,
            (MolfileBondOrder::Aromatic, _) => {
                return Err(lines.error(
                    row.index,
                    DiagnosticKind::Unsupported,
                    format!("V3000 bond CFG value {cfg} is unsupported for an aromatic bond"),
                ))
            }
            (MolfileBondOrder::Single, 1) => Some(MolfileBondStereo::Up),
            (MolfileBondOrder::Single, 2) => Some(MolfileBondStereo::Either),
            (MolfileBondOrder::Single, 3) => Some(MolfileBondStereo::Down),
            (MolfileBondOrder::Double, 2) => Some(MolfileBondStereo::Either),
            _ => {
                return Err(lines.error(
                    row.index,
                    DiagnosticKind::InvalidValue,
                    format!("V3000 bond CFG value {cfg} is invalid for bond order {order:?}"),
                ))
            }
        };
    }
    Ok((
        serial,
        MolfileBond {
            atoms,
            order,
            stereo,
            location: lines.location(row.index),
        },
    ))
}

fn is_v3000_stereo_group(name: &str) -> bool {
    name == "MDLV30/STEABS"
        || name.starts_with("MDLV30/STEREL")
        || name.starts_with("MDLV30/STERAC")
}

/// Reads enhanced stereo groups; every other collection is unsupported content.
fn parse_v3000_collection(
    lines: MolfileLines<'_>,
    rows: &[V3000Line],
    atom_positions: &BTreeMap<usize, usize>,
    ctab: &mut MolfileCtab,
) -> std::result::Result<(), MolfileParseError> {
    for row in rows {
        let error = |kind, message: &str| lines.error(row.index, kind, message);
        let syntax = |message: &str| error(DiagnosticKind::Syntax, message);
        let fields =
            v3000_fields(&row.body).ok_or_else(|| syntax("invalid V3000 collection line"))?;
        let Some((&name, lists)) = fields.split_first() else {
            continue;
        };
        let number = |digits: &str| match digits.parse::<usize>() {
            Ok(number) if number > 0 => Ok(number),
            _ => Err(error(
                DiagnosticKind::InvalidValue,
                "invalid V3000 stereo group number",
            )),
        };
        let kind = if name == "MDLV30/STEABS" {
            MolfileStereoGroupKind::Absolute
        } else if let Some(digits) = name.strip_prefix("MDLV30/STEREL") {
            MolfileStereoGroupKind::Or(number(digits)?)
        } else if let Some(digits) = name.strip_prefix("MDLV30/STERAC") {
            MolfileStereoGroupKind::And(number(digits)?)
        } else {
            ctab.unsupported.push(MolfileUnsupported {
                content: MolfileUnsupportedContent::Collection(name.to_owned()),
                location: lines.location(row.index),
            });
            continue;
        };
        let mut atoms = None;
        for field in lists {
            let (key, value) = split_v3000_option(field)
                .ok_or_else(|| syntax("invalid V3000 stereo group member list"))?;
            if key != "ATOMS" {
                return Err(error(
                    DiagnosticKind::Unsupported,
                    &format!("unsupported V3000 stereo group member list `{key}`"),
                ));
            }
            if atoms.is_some() {
                return Err(syntax("duplicate V3000 stereo group member list"));
            }
            let list = value
                .strip_prefix('(')
                .and_then(|list| list.strip_suffix(')'))
                .ok_or_else(|| syntax("invalid V3000 stereo group member list"))?;
            let mut values = list.split_whitespace();
            let count = values
                .next()
                .and_then(|count| count.parse::<usize>().ok())
                .ok_or_else(|| syntax("invalid V3000 stereo group member count"))?;
            let members = values
                .map(|serial| {
                    let serial = serial
                        .parse::<usize>()
                        .map_err(|_| syntax("invalid V3000 stereo group member"))?;
                    atom_positions.get(&serial).copied().ok_or_else(|| {
                        error(
                            DiagnosticKind::InvalidReference,
                            "V3000 stereo group member is outside the atom block",
                        )
                    })
                })
                .collect::<std::result::Result<Vec<_>, _>>()?;
            if members.len() != count {
                return Err(error(
                    DiagnosticKind::Inconsistent,
                    "V3000 stereo group member count mismatch",
                ));
            }
            atoms = Some(members);
        }
        let atoms = atoms.filter(|atoms| !atoms.is_empty()).ok_or_else(|| {
            error(
                DiagnosticKind::InvalidValue,
                "V3000 stereo group must have members",
            )
        })?;
        if let Some(group) = ctab
            .stereo_groups
            .iter_mut()
            .find(|group| group.kind == kind)
        {
            group.atoms.extend(atoms);
        } else {
            ctab.stereo_groups.push(MolfileStereoGroup {
                kind,
                atoms,
                location: lines.location(row.index),
            });
        }
    }
    Ok(())
}

fn v3000_bond_code(order: BondOrder) -> std::result::Result<u8, MoleculeWriteError> {
    match order {
        BondOrder::Zero => Ok(0),
        BondOrder::Single => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Dative => Ok(9),
        BondOrder::Quadruple => Err(MoleculeWriteError::new(
            "V3000 writer does not support quadruple bonds",
        )),
    }
}

fn v3000_bond_cfg(
    order: BondOrder,
    stereo: Option<SourceStereoBondMarkKind>,
) -> std::result::Result<Option<u8>, MoleculeWriteError> {
    match (order, stereo) {
        (_, None) => Ok(None),
        (BondOrder::Double, Some(SourceStereoBondMarkKind::MolfileDoubleBondGeometry)) => Ok(None),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeUp)) => Ok(Some(1)),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeEither)) => Ok(Some(2)),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeDown)) => Ok(Some(3)),
        (BondOrder::Double, Some(SourceStereoBondMarkKind::DoubleBondEither)) => Ok(Some(2)),
        _ => Err(MoleculeWriteError::new(
            "V3000 bond CFG is incompatible with the bond order",
        )),
    }
}

fn v3000_radical_code(radical: AtomRadical) -> std::result::Result<u8, MoleculeWriteError> {
    match (radical.electron_count(), radical.spin_multiplicity()) {
        (2, Some(1)) => Ok(1),
        (1, Some(2)) => Ok(2),
        (2, Some(3)) => Ok(3),
        _ => Err(MoleculeWriteError::new(
            "V3000 writer requires an encodable radical electron count and explicit spin multiplicity",
        )),
    }
}

fn split_v3000_option(field: &str) -> Option<(&str, &str)> {
    let (key, value) = field.split_once('=')?;
    (!key.is_empty() && !value.is_empty()).then_some((key, value))
}

fn parse_finite_f64(value: &str) -> Option<f64> {
    let parsed: f64 = value.parse().ok()?;
    parsed.is_finite().then_some(parsed)
}
