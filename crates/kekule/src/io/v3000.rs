//! V3000 (`M  V30`) connection-table parsing.

use std::collections::{BTreeMap, BTreeSet};

use crate::diagnostics::{DiagnosticKind, LimitExceeded, Resource};
use crate::geometry::Point3;
use crate::io::MolfileParseOptions;

use super::molfile_document::{
    MolfileAtom, MolfileAtomParity, MolfileBond, MolfileBondOrder, MolfileBondStereo, MolfileCtab,
    MolfileLines, MolfileRadical, MolfileStereoGroup, MolfileStereoGroupKind, MolfileUnsupported,
    MolfileUnsupportedContent,
};
use crate::diagnostics::FormatError;

/// Parses the V3000 connection table after the header and returns the index
/// of its `M  END` line.
pub(super) fn parse_v3000_ctab(
    lines: MolfileLines<'_>,
    options: MolfileParseOptions,
) -> std::result::Result<(MolfileCtab, usize), FormatError> {
    let (records, end) = collect_v3000_lines(lines, options)?;
    let Some(first) = records.first() else {
        return Err(lines.error(
            end,
            DiagnosticKind::Syntax,
            "The \"BEGIN CTAB\" line is missing.",
        ));
    };
    if first.body != "BEGIN CTAB" {
        return Err(lines.error(
            first.index,
            DiagnosticKind::Syntax,
            "The V3000 record lies outside the CTAB.",
        ));
    }
    let Some(close) = records.iter().position(|record| record.body == "END CTAB") else {
        return Err(lines.error(
            end,
            DiagnosticKind::Syntax,
            "The \"END CTAB\" line is missing.",
        ));
    };
    if let Some(outside) = records.get(close + 1) {
        return Err(lines.error(
            outside.index,
            DiagnosticKind::Syntax,
            "The V3000 record lies outside the CTAB.",
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
                .ok_or_else(|| syntax(&format!("The \"END {name}\" line is missing.")))?;
            let rows = &body[position..position + length];
            position += length + 1;
            match name {
                "ATOM" => {
                    if atoms_read {
                        return Err(syntax("The CTAB has a second ATOM section."));
                    }
                    let Some((counts, counts_index)) = counts else {
                        return Err(syntax("The ATOM section comes before the COUNTS line."));
                    };
                    if rows.len() != counts.atoms {
                        return Err(lines.error(
                            counts_index,
                            DiagnosticKind::Inconsistent,
                            format!(
                                "The COUNTS line declares {} atoms, but the ATOM section has {}.",
                                counts.atoms,
                                rows.len()
                            ),
                        ));
                    }
                    for row in rows {
                        let (serial, atom) =
                            parse_v3000_atom(lines, row, ctab.atoms.len(), &mut ctab.unsupported)?;
                        if atom_positions.insert(serial, ctab.atoms.len()).is_some() {
                            return Err(lines.error(
                                row.index,
                                DiagnosticKind::Inconsistent,
                                format!("Atom index {serial} appears twice."),
                            ));
                        }
                        ctab.atoms.push(atom);
                    }
                    atoms_read = true;
                }
                "BOND" => {
                    let (Some((counts, counts_index)), true) = (counts, atoms_read) else {
                        return Err(syntax("The BOND section comes before the ATOM section."));
                    };
                    if bonds_read {
                        return Err(syntax("The CTAB has a second BOND section."));
                    }
                    if collection_read {
                        return Err(syntax(
                            "The BOND section comes after the COLLECTION section.",
                        ));
                    }
                    if rows.len() != counts.bonds {
                        return Err(lines.error(
                            counts_index,
                            DiagnosticKind::Inconsistent,
                            format!(
                                "The COUNTS line declares {} bonds, but the BOND section has {}.",
                                counts.bonds,
                                rows.len()
                            ),
                        ));
                    }
                    parse_v3000_bonds(lines, rows, &atom_positions, &mut ctab)?;
                    bonds_read = true;
                }
                "COLLECTION" => {
                    if !bonds_read && !bond_block_optional(counts, atoms_read) {
                        return Err(syntax(
                            "The COLLECTION section comes before the BOND section.",
                        ));
                    }
                    if collection_read {
                        return Err(syntax("The CTAB has a second COLLECTION section."));
                    }
                    parse_v3000_collection(lines, rows, &atom_positions, &mut ctab)?;
                    collection_read = true;
                }
                "CTAB" => return Err(syntax("The Molfile has a second CTAB.")),
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
                return Err(syntax("The CTAB has a second COUNTS line."));
            }
            let parsed = parse_v3000_counts(&record.body)
                .ok_or_else(|| syntax("The COUNTS line is malformed."))?;
            for limit in [
                LimitExceeded::check(Resource::Atoms, options.max_v3000_atoms(), parsed.atoms),
                LimitExceeded::check(Resource::Bonds, options.max_v3000_bonds(), parsed.bonds),
            ] {
                if let Err(limit) = limit {
                    return Err(lines
                        .error(
                            record.index,
                            DiagnosticKind::Limit,
                            format!("The COUNTS line exceeds the {}.", limit.describe()),
                        )
                        .exceeding(limit));
                }
            }
            ctab.chiral_flag = parsed.chiral;
            counts = Some((parsed, record.index));
        } else if keyword == "END" {
            return Err(syntax(&format!(
                "The line {:?} closes no section.",
                record.body
            )));
        } else if is_v3000_stereo_group(keyword) {
            return Err(syntax(
                "The stereo group lies outside the COLLECTION section.",
            ));
        } else {
            ctab.unsupported.push(MolfileUnsupported {
                content: MolfileUnsupportedContent::Record(keyword.to_owned()),
                location: lines.location(record.index),
            });
        }
    }
    if !atoms_read {
        return Err(lines.error(end, DiagnosticKind::Syntax, "The ATOM section is missing."));
    }
    if !bonds_read && !bond_block_optional(counts, atoms_read) {
        return Err(lines.error(end, DiagnosticKind::Syntax, "The BOND section is missing."));
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
) -> std::result::Result<(Vec<V3000Line>, usize), FormatError> {
    let limit = |index| {
        let limit = LimitExceeded::first_beyond(
            Resource::RecordBytes,
            options.max_v3000_logical_line_bytes(),
        );
        lines
            .error(
                index,
                DiagnosticKind::Limit,
                format!("The V3000 line exceeds its {}.", limit.describe()),
            )
            .exceeding(limit)
    };
    let mut records = Vec::new();
    let mut index = 4usize;
    while let Some(line) = lines.text(index) {
        let start = index;
        if line.trim() == "M  END" {
            return Ok((records, index));
        }
        let body = v3000_body(line).ok_or_else(|| {
            lines.error(
                index,
                DiagnosticKind::Syntax,
                "The line does not start with \"M  V30\".",
            )
        })?;
        if body.len() > options.max_v3000_logical_line_bytes() {
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
                    "The continued V3000 line has no next line.",
                )
            })?;
            let continuation = v3000_body(continuation_line)
                .ok_or_else(|| {
                    lines.error(
                        index,
                        DiagnosticKind::Syntax,
                        "The continuation line does not start with \"M  V30\".",
                    )
                })?
                .trim_start();
            if body.len().saturating_add(continuation.len())
                > options.max_v3000_logical_line_bytes()
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
        "The \"M  END\" line is missing.",
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
) -> std::result::Result<(usize, MolfileAtom), FormatError> {
    let syntax = || {
        lines.error(
            row.index,
            DiagnosticKind::Syntax,
            format!("The line of atom {} is malformed.", position + 1),
        )
    };
    let invalid = |name: &str| {
        lines.error(
            row.index,
            DiagnosticKind::InvalidValue,
            format!("The {name} value of atom {} is not valid.", position + 1),
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
            "The atom index is 0. Atom indices start at 1.",
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
                format!("Atom {serial} has the option {key:?} twice."),
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
            // HCOUNT is a query constraint ("at least n"), not a molecular
            // hydrogen declaration; only a nonzero value states one.
            "HCOUNT" => {
                if parse_v3000_count(value)
                    .ok_or_else(|| invalid("HCOUNT"))?
                    .is_some()
                {
                    unsupported.push(MolfileUnsupported {
                        content: MolfileUnsupportedContent::AtomOption {
                            atom: position,
                            key: key.to_owned(),
                        },
                        location: lines.location(row.index),
                    });
                }
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
) -> std::result::Result<(), FormatError> {
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
            |message: String| lines.error(row.index, DiagnosticKind::Inconsistent, message);
        if !serials.insert(serial) {
            return Err(inconsistent(format!("Bond index {serial} appears twice.")));
        }
        let [a, b] = bond.atoms;
        if !endpoints.insert((a.min(b), a.max(b))) {
            return Err(inconsistent(format!(
                "Bond {serial} bonds two atoms that another bond already joins."
            )));
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
) -> std::result::Result<(usize, MolfileBond), FormatError> {
    let syntax = || {
        lines.error(
            row.index,
            DiagnosticKind::Syntax,
            format!("The line of bond {} is malformed.", position + 1),
        )
    };
    let fields = v3000_fields(&row.body).ok_or_else(syntax)?;
    let [serial, order, a, b, options @ ..] = fields.as_slice() else {
        return Err(syntax());
    };
    let serial = serial.parse::<usize>().map_err(|_| syntax())?;
    if serial == 0 {
        return Err(lines.error(
            row.index,
            DiagnosticKind::InvalidValue,
            "The bond index is 0. Bond indices start at 1.",
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
                format!("The bond order code {order_code} of bond {serial} is not supported."),
            ))
        }
        _ => {
            return Err(lines.error(
                row.index,
                DiagnosticKind::InvalidValue,
                format!("The bond order code {order_code} of bond {serial} is not valid."),
            ))
        }
    };
    let endpoint = |atom: &str| {
        let atom = atom.parse::<usize>().map_err(|_| syntax())?;
        atom_positions.get(&atom).copied().ok_or_else(|| {
            lines.error(
                row.index,
                DiagnosticKind::InvalidReference,
                format!("Bond {serial} refers to atom {atom}, which the ATOM section lacks."),
            )
        })
    };
    let atoms = [endpoint(a)?, endpoint(b)?];
    if atoms[0] == atoms[1] {
        return Err(lines.error(
            row.index,
            DiagnosticKind::InvalidValue,
            format!("Bond {serial} joins atom {a} to itself."),
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
                format!("Bond {serial} has the option {key:?} twice."),
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
                    format!("The CFG value {cfg} of aromatic bond {serial} is not supported."),
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
                    format!(
                        "The CFG value {cfg} of bond {serial} does not apply to its bond order."
                    ),
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
) -> std::result::Result<(), FormatError> {
    for row in rows {
        let error = |kind, message: &str| lines.error(row.index, kind, message);
        let syntax = |message: &str| error(DiagnosticKind::Syntax, message);
        let fields =
            v3000_fields(&row.body).ok_or_else(|| syntax("The collection line is malformed."))?;
        let Some((&name, lists)) = fields.split_first() else {
            continue;
        };
        let number = |digits: &str| match digits.parse::<usize>() {
            Ok(number) if number > 0 => Ok(number),
            _ => Err(error(
                DiagnosticKind::InvalidValue,
                &format!("The stereo group {name:?} has no positive number."),
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
            let (key, value) = split_v3000_option(field).ok_or_else(|| {
                syntax(&format!(
                    "The atom list of stereo group {name:?} is malformed."
                ))
            })?;
            if key != "ATOMS" {
                return Err(error(
                    DiagnosticKind::Unsupported,
                    &format!("The stereo group list {key:?} is not supported."),
                ));
            }
            if atoms.is_some() {
                return Err(syntax(&format!(
                    "The stereo group {name:?} has two atom lists."
                )));
            }
            let list = value
                .strip_prefix('(')
                .and_then(|list| list.strip_suffix(')'))
                .ok_or_else(|| {
                    syntax(&format!(
                        "The atom list of stereo group {name:?} is malformed."
                    ))
                })?;
            let mut values = list.split_whitespace();
            let count = values
                .next()
                .and_then(|count| count.parse::<usize>().ok())
                .ok_or_else(|| {
                    syntax(&format!(
                        "The atom list of stereo group {name:?} has no count."
                    ))
                })?;
            let members = values
                .map(|serial| {
                    let atom = serial.parse::<usize>().map_err(|_| {
                        syntax(&format!(
                            "The stereo group {name:?} lists {serial:?}, which is not an atom index."
                        ))
                    })?;
                    atom_positions.get(&atom).copied().ok_or_else(|| {
                        error(
                            DiagnosticKind::InvalidReference,
                            &format!(
                                "The stereo group {name:?} lists atom {atom}, which the ATOM section lacks."
                            ),
                        )
                    })
                })
                .collect::<std::result::Result<Vec<_>, _>>()?;
            if members.len() != count {
                return Err(error(
                    DiagnosticKind::Inconsistent,
                    &format!(
                        "The stereo group {name:?} declares {} but lists {}.",
                        crate::diagnostics::counted(count, "atom", "atoms"),
                        members.len()
                    ),
                ));
            }
            atoms = Some(members);
        }
        let atoms = atoms.filter(|atoms| !atoms.is_empty()).ok_or_else(|| {
            error(
                DiagnosticKind::InvalidValue,
                &format!("The stereo group {name:?} has no atoms."),
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

fn split_v3000_option(field: &str) -> Option<(&str, &str)> {
    let (key, value) = field.split_once('=')?;
    (!key.is_empty() && !value.is_empty()).then_some((key, value))
}

fn parse_finite_f64(value: &str) -> Option<f64> {
    let parsed: f64 = value.parse().ok()?;
    parsed.is_finite().then_some(parsed)
}
