use std::collections::BTreeMap;

use crate::core::*;
use crate::io::SdfParseError;

pub fn read_mol_v3000_str(input: &str) -> std::result::Result<SmallMolecule, SdfParseError> {
    let normalized = input.replace("\r\n", "\n").replace('\r', "\n");
    let lines = normalized.lines().collect::<Vec<_>>();
    parse_mol_v3000_lines(1, 1, &lines)
}

fn parse_mol_v3000_lines(
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
    let counts_line = checked_line_number(record, start_line, 3)?;
    if !lines[3].contains("V3000") {
        return Err(SdfParseError::new(
            record,
            counts_line,
            "counts line must declare V3000",
        ));
    }

    let v30_lines = collect_v3000_lines(record, start_line, lines)?;
    let ctab = v3000_section(record, &v30_lines, "CTAB", 0)?;
    let counts_index = find_v3000_record(&v30_lines, "COUNTS", ctab.start + 1, ctab.end)
        .ok_or_else(|| SdfParseError::new(record, counts_line, "missing V3000 COUNTS line"))?;
    let counts = parse_v3000_counts(&v30_lines[counts_index].body).ok_or_else(|| {
        SdfParseError::new(
            record,
            v30_lines[counts_index].line,
            "invalid V3000 COUNTS line",
        )
    })?;

    let atom_section = v3000_section(record, &v30_lines, "ATOM", ctab.start + 1)?;
    let bond_section = v3000_section(record, &v30_lines, "BOND", atom_section.end + 1)?;
    if atom_section.end > ctab.end || bond_section.end > ctab.end {
        return Err(SdfParseError::new(
            record,
            v30_lines[ctab.end].line,
            "V3000 ATOM/BOND section escapes CTAB",
        ));
    }

    let atom_rows = &v30_lines[atom_section.start + 1..atom_section.end];
    let bond_rows = &v30_lines[bond_section.start + 1..bond_section.end];
    if atom_rows.len() != counts.atoms || bond_rows.len() != counts.bonds {
        return Err(SdfParseError::new(
            record,
            v30_lines[counts_index].line,
            "V3000 COUNTS do not match ATOM/BOND section sizes",
        ));
    }

    let mut mol = Molecule::new();
    mol.props_mut()
        .insert("sdf.title".to_owned(), PropValue::String(title));
    mol.props_mut().insert(
        "sdf.program".to_owned(),
        PropValue::String(lines[1].to_owned()),
    );
    mol.props_mut().insert(
        "sdf.comment".to_owned(),
        PropValue::String(lines[2].to_owned()),
    );

    let mut atom_ids = BTreeMap::<usize, AtomId>::new();
    let mut conformer = Conformer::with_atom_capacity(atom_rows.len());
    for row in atom_rows {
        let parsed = parse_v3000_atom(&row.body)
            .ok_or_else(|| SdfParseError::new(record, row.line, "invalid V3000 atom line"))?;
        if atom_ids.contains_key(&parsed.index) {
            return Err(SdfParseError::new(
                record,
                row.line,
                "duplicate V3000 atom index",
            ));
        }
        let element = Element::from_symbol(parsed.symbol).ok_or_else(|| {
            SdfParseError::new(
                record,
                row.line,
                format!("unknown element symbol `{}`", parsed.symbol),
            )
        })?;
        let mut atom = Atom::new(element);
        atom.atom_map = (parsed.atom_map != 0).then_some(parsed.atom_map);
        apply_v3000_atom_options(record, row.line, &mut atom, &parsed.options)?;
        let atom_id = mol.add_atom(atom);
        conformer.set_position(atom_id, parsed.point);
        atom_ids.insert(parsed.index, atom_id);
    }

    for row in bond_rows {
        let parsed = parse_v3000_bond(&row.body)
            .ok_or_else(|| SdfParseError::new(record, row.line, "invalid V3000 bond line"))?;
        let a = *atom_ids.get(&parsed.a).ok_or_else(|| {
            SdfParseError::new(record, row.line, "bond endpoint outside atom block")
        })?;
        let b = *atom_ids.get(&parsed.b).ok_or_else(|| {
            SdfParseError::new(record, row.line, "bond endpoint outside atom block")
        })?;
        let bond_id = mol
            .add_bond(a, b, parsed.order)
            .map_err(|error| SdfParseError::new(record, row.line, error.to_string()))?;
        if let Some(stereo) = parsed.stereo {
            mol.bond_mut(bond_id)
                .expect("newly added bond should be mutable")
                .stereo = Some(stereo);
        }
    }

    if conformer.positions().next().is_some() {
        mol.add_conformer(conformer);
    }
    Ok(SmallMolecule { mol })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct V3000Counts {
    atoms: usize,
    bonds: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct V3000Section {
    start: usize,
    end: usize,
}

#[derive(Debug, Clone, PartialEq)]
struct V3000Line {
    line: usize,
    body: String,
}

#[derive(Debug, Clone, PartialEq)]
struct V3000Atom<'a> {
    index: usize,
    symbol: &'a str,
    point: Point3,
    atom_map: u32,
    options: Vec<(&'a str, &'a str)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct V3000Bond {
    order: BondOrder,
    a: usize,
    b: usize,
    stereo: Option<BondStereo>,
}

fn collect_v3000_lines(
    record: usize,
    start_line: usize,
    lines: &[&str],
) -> std::result::Result<Vec<V3000Line>, SdfParseError> {
    let mut records = Vec::new();
    let mut index = 4usize;
    while index < lines.len() {
        let line_number = checked_line_number(record, start_line, index)?;
        let line = lines[index];
        if line.trim() == "M  END" {
            return Ok(records);
        }
        let mut body = v3000_body(line)
            .ok_or_else(|| SdfParseError::new(record, line_number, "expected M  V30 record"))?;
        while body.ends_with('-') {
            body.pop();
            index = index.checked_add(1).ok_or_else(|| {
                SdfParseError::new(record, line_number, "V3000 continuation overflow")
            })?;
            let continuation_line = lines.get(index).copied().ok_or_else(|| {
                SdfParseError::new(record, line_number, "unterminated V3000 continuation")
            })?;
            let continuation = v3000_body(continuation_line).ok_or_else(|| {
                SdfParseError::new(record, line_number, "invalid V3000 continuation")
            })?;
            body.push_str(continuation.trim_start());
        }
        records.push(V3000Line {
            line: line_number,
            body,
        });
        index += 1;
    }
    Err(SdfParseError::new(record, start_line, "missing M  END"))
}

fn v3000_body(line: &str) -> Option<String> {
    let trimmed = line.strip_prefix("M  V30 ")?;
    Some(trimmed.trim().to_owned())
}

fn v3000_section(
    record: usize,
    lines: &[V3000Line],
    name: &str,
    search_start: usize,
) -> std::result::Result<V3000Section, SdfParseError> {
    let begin = format!("BEGIN {name}");
    let end = format!("END {name}");
    let start = lines
        .iter()
        .enumerate()
        .skip(search_start)
        .find_map(|(index, line)| (line.body == begin).then_some(index))
        .ok_or_else(|| SdfParseError::new(record, 1, format!("missing V3000 BEGIN {name}")))?;
    let end = lines
        .iter()
        .enumerate()
        .skip(start + 1)
        .find_map(|(index, line)| (line.body == end).then_some(index))
        .ok_or_else(|| {
            SdfParseError::new(
                record,
                lines[start].line,
                format!("missing V3000 END {name}"),
            )
        })?;
    Ok(V3000Section { start, end })
}

fn find_v3000_record(
    lines: &[V3000Line],
    keyword: &str,
    start: usize,
    end: usize,
) -> Option<usize> {
    lines
        .iter()
        .enumerate()
        .take(end)
        .skip(start)
        .find_map(|(index, line)| line.body.starts_with(keyword).then_some(index))
}

fn parse_v3000_counts(line: &str) -> Option<V3000Counts> {
    let mut fields = line.split_whitespace();
    (fields.next()? == "COUNTS").then_some(())?;
    Some(V3000Counts {
        atoms: fields.next()?.parse().ok()?,
        bonds: fields.next()?.parse().ok()?,
    })
}

fn parse_v3000_atom(line: &str) -> Option<V3000Atom<'_>> {
    let mut fields = line.split_whitespace();
    let index = fields.next()?.parse().ok()?;
    let symbol = fields.next()?;
    let point = Point3::new(
        parse_finite_f64(fields.next()?)?,
        parse_finite_f64(fields.next()?)?,
        parse_finite_f64(fields.next()?)?,
    );
    let atom_map = fields.next()?.parse().ok()?;
    let options = fields.filter_map(split_v3000_option).collect();
    Some(V3000Atom {
        index,
        symbol,
        point,
        atom_map,
        options,
    })
}

fn apply_v3000_atom_options(
    record: usize,
    line: usize,
    atom: &mut Atom,
    options: &[(&str, &str)],
) -> std::result::Result<(), SdfParseError> {
    for (key, value) in options {
        match *key {
            "CHG" => {
                atom.formal_charge = value
                    .parse()
                    .map_err(|_| SdfParseError::new(record, line, "invalid V3000 CHG value"))?;
            }
            "MASS" => {
                let isotope = value
                    .parse::<u16>()
                    .map_err(|_| SdfParseError::new(record, line, "invalid V3000 MASS value"))?;
                atom.isotope = (isotope != 0).then_some(isotope);
            }
            "RAD" => {
                atom.radical = Some(match *value {
                    "1" => AtomRadical::Singlet,
                    "2" => AtomRadical::Doublet,
                    "3" => AtomRadical::Triplet,
                    _ => {
                        return Err(SdfParseError::new(
                            record,
                            line,
                            "unsupported V3000 RAD code",
                        ))
                    }
                });
            }
            "CFG" => {
                return Err(SdfParseError::new(
                    record,
                    line,
                    "V3000 atom stereochemistry is not supported",
                ));
            }
            _ => {}
        }
    }
    Ok(())
}

fn parse_v3000_bond(line: &str) -> Option<V3000Bond> {
    let mut fields = line.split_whitespace();
    let _index = fields.next()?.parse::<usize>().ok()?;
    let order_code = fields.next()?.parse().ok()?;
    let a = fields.next()?.parse().ok()?;
    let b = fields.next()?.parse().ok()?;
    let order = v3000_bond_order(order_code)?;
    let mut stereo = None;
    for (key, value) in fields.filter_map(split_v3000_option) {
        if key == "CFG" {
            stereo = Some(v3000_bond_stereo(order, value)?);
        }
    }
    Some(V3000Bond {
        order,
        a,
        b,
        stereo,
    })
}

fn v3000_bond_order(code: u8) -> Option<BondOrder> {
    match code {
        0 => Some(BondOrder::Zero),
        1 => Some(BondOrder::Single),
        2 => Some(BondOrder::Double),
        3 => Some(BondOrder::Triple),
        4 => Some(BondOrder::Aromatic),
        9 => Some(BondOrder::Dative),
        _ => None,
    }
}

fn v3000_bond_stereo(order: BondOrder, value: &str) -> Option<BondStereo> {
    match (order, value) {
        (_, "0") => Some(BondStereo::Unspecified),
        (BondOrder::Single, "1") => Some(BondStereo::Up),
        (BondOrder::Single, "2") => Some(BondStereo::Any),
        (BondOrder::Single, "3") => Some(BondStereo::Down),
        (BondOrder::Double, "2") => Some(BondStereo::Any),
        _ => None,
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

fn checked_line_number(
    record: usize,
    start_line: usize,
    offset: usize,
) -> std::result::Result<usize, SdfParseError> {
    start_line
        .checked_add(offset)
        .ok_or_else(|| SdfParseError::new(record, start_line, "line number overflow"))
}
