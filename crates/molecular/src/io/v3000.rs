use std::collections::BTreeMap;

use crate::core::*;
use crate::io::{
    preserve_molfile_tetrahedral_hydrogens, MolWriteError, MolfileParseOptions, SdfParseError,
};
use crate::small::model::SmallMolecule;
use crate::units::{Quantity, ANGSTROM};

#[derive(Debug, Clone, PartialEq)]
pub(super) struct V3000Syntax {
    pub(super) atoms: Vec<V3000AtomSyntax>,
    pub(super) bonds: Vec<V3000BondSyntax>,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct V3000AtomSyntax {
    pub(super) index: usize,
    pub(super) atom: Atom,
    pub(super) point: Point3,
    pub(super) line: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct V3000BondSyntax {
    pub(super) a: usize,
    pub(super) b: usize,
    pub(super) order: BondOrder,
    pub(super) stereo: Option<StereoBondMarkKind>,
    pub(super) line: usize,
}

pub fn write_mol_v3000(molecule: &SmallMolecule) -> std::result::Result<String, MolWriteError> {
    let mol = molecule.graph();
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
    out.push_str("  0  0  0  0  0  0            999 V3000\n");
    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {} 0 0 0\n",
        atoms.len(),
        bonds.len()
    ));
    out.push_str("M  V30 BEGIN ATOM\n");
    if mol.stereo_elements().next().is_some() {
        return Err(MolWriteError::new(
            "V3000 writer does not support stereo elements",
        ));
    }
    for atom_id in &atoms {
        let atom = mol
            .atom(*atom_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let point = conformer
            .and_then(|conformer| conformer.position(*atom_id))
            .map(|point| point.value_in(ANGSTROM).expect("conformer length unit"))
            .unwrap_or_default();
        let index = atom_index
            .get(atom_id)
            .ok_or_else(|| MolWriteError::new("atom missing from V3000 atom table"))?;
        out.push_str(&format!(
            "M  V30 {index} {} {:.4} {:.4} {:.4} {}",
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
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
    out.push_str("M  V30 BEGIN BOND\n");
    for (index, bond_id) in bonds.iter().enumerate() {
        let bond = mol
            .bond(*bond_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let a = atom_index
            .get(&bond.a())
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from V3000 atom table"))?;
        let b = atom_index
            .get(&bond.b())
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from V3000 atom table"))?;
        let order_code = v3000_bond_code(bond.order)?;
        out.push_str(&format!("M  V30 {} {order_code} {a} {b}", index + 1));
        if let Some(cfg) = v3000_bond_cfg(bond.order, mol.stereo_bond_mark(*bond_id))? {
            out.push_str(&format!(" CFG={cfg}"));
        }
        out.push('\n');
    }
    out.push_str("M  V30 END BOND\n");
    out.push_str("M  V30 END CTAB\n");
    out.push_str("M  END\n");
    Ok(out)
}

pub(super) fn parse_v3000_syntax(
    record: usize,
    start_line: usize,
    lines: &[&str],
    options: MolfileParseOptions,
) -> std::result::Result<V3000Syntax, SdfParseError> {
    if lines.len() < 4 {
        return Err(SdfParseError::new(
            record,
            start_line,
            "record must contain three header lines and a counts line",
        ));
    }
    let counts_line = checked_line_number(record, start_line, 3)?;
    if !lines[3].contains("V3000") {
        return Err(SdfParseError::new(
            record,
            counts_line,
            "counts line must declare V3000",
        ));
    }

    let v30_lines = collect_v3000_lines(record, start_line, lines, options)?;
    for control in [
        "BEGIN CTAB",
        "END CTAB",
        "BEGIN ATOM",
        "END ATOM",
        "BEGIN BOND",
        "END BOND",
    ] {
        if v30_lines.iter().filter(|line| line.body == control).count() != 1 {
            return Err(SdfParseError::new(
                record,
                counts_line,
                format!("V3000 must contain exactly one `{control}` control record"),
            ));
        }
    }
    let ctab = v3000_section(record, &v30_lines, "CTAB", 0)?;
    if ctab.start != 0 || ctab.end + 1 != v30_lines.len() {
        return Err(SdfParseError::new(
            record,
            v30_lines
                .get(ctab.start)
                .map_or(counts_line, |line| line.line),
            "V3000 CTAB must contain every V30 record",
        ));
    }
    let atom_section = v3000_section(record, &v30_lines, "ATOM", ctab.start + 1)?;
    let counts_indexes = v30_lines
        .iter()
        .enumerate()
        .filter_map(|(index, line)| {
            (line.body.split_whitespace().next() == Some("COUNTS")).then_some(index)
        })
        .collect::<Vec<_>>();
    let [counts_index] = counts_indexes.as_slice() else {
        return Err(SdfParseError::new(
            record,
            counts_line,
            "V3000 CTAB must contain exactly one COUNTS line before the ATOM section",
        ));
    };
    let counts_index = *counts_index;
    if counts_index <= ctab.start || counts_index >= atom_section.start {
        return Err(SdfParseError::new(
            record,
            v30_lines[counts_index].line,
            "V3000 COUNTS line must occur before the ATOM section inside the CTAB",
        ));
    }
    let counts = parse_v3000_counts(&v30_lines[counts_index].body).ok_or_else(|| {
        SdfParseError::new(
            record,
            v30_lines[counts_index].line,
            "invalid V3000 COUNTS line",
        )
    })?;
    if counts.atoms > options.max_v3000_atoms {
        return Err(SdfParseError::new(
            record,
            v30_lines[counts_index].line,
            "V3000 atom count exceeds configured limit",
        ));
    }
    if counts.bonds > options.max_v3000_bonds {
        return Err(SdfParseError::new(
            record,
            v30_lines[counts_index].line,
            "V3000 bond count exceeds configured limit",
        ));
    }

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

    let mut atoms = Vec::with_capacity(atom_rows.len());
    let mut atom_indices = BTreeMap::<usize, usize>::new();
    for row in atom_rows {
        let parsed = parse_v3000_atom(&row.body)
            .ok_or_else(|| SdfParseError::new(record, row.line, "invalid V3000 atom line"))?;
        if parsed.index == 0 {
            return Err(SdfParseError::new(
                record,
                row.line,
                "V3000 atom indices must be positive",
            ));
        }
        if atom_indices.contains_key(&parsed.index) {
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
        atom_indices.insert(parsed.index, atoms.len());
        atoms.push(V3000AtomSyntax {
            index: parsed.index,
            atom,
            point: parsed.point,
            line: row.line,
        });
    }

    let mut bonds = Vec::with_capacity(bond_rows.len());
    let mut bond_indices = std::collections::BTreeSet::new();
    let mut endpoints = std::collections::BTreeSet::new();
    for row in bond_rows {
        let parsed = parse_v3000_bond(record, row.line, &row.body)?;
        if parsed.index == 0 {
            return Err(SdfParseError::new(
                record,
                row.line,
                "V3000 bond indices must be positive",
            ));
        }
        if !bond_indices.insert(parsed.index) {
            return Err(SdfParseError::new(
                record,
                row.line,
                "duplicate V3000 bond index",
            ));
        }
        atom_indices.get(&parsed.a).ok_or_else(|| {
            SdfParseError::new(record, row.line, "bond endpoint outside atom block")
        })?;
        atom_indices.get(&parsed.b).ok_or_else(|| {
            SdfParseError::new(record, row.line, "bond endpoint outside atom block")
        })?;
        if parsed.a == parsed.b {
            return Err(SdfParseError::new(
                record,
                row.line,
                "bond endpoints must be distinct",
            ));
        }
        let ordered = if parsed.a < parsed.b {
            (parsed.a, parsed.b)
        } else {
            (parsed.b, parsed.a)
        };
        if !endpoints.insert(ordered) {
            return Err(SdfParseError::new(
                record,
                row.line,
                "duplicate bond endpoints",
            ));
        }
        bonds.push(V3000BondSyntax {
            a: parsed.a,
            b: parsed.b,
            order: parsed.order,
            stereo: parsed.stereo,
            line: row.line,
        });
    }

    Ok(V3000Syntax { atoms, bonds })
}

pub(super) fn interpret_v3000_syntax(
    syntax: &V3000Syntax,
) -> std::result::Result<SmallMolecule, SdfParseError> {
    let mut mol = Molecule::new();
    let mut atom_ids = BTreeMap::<usize, AtomId>::new();
    let mut conformer = Conformer::with_atom_capacity(syntax.atoms.len(), ANGSTROM)
        .expect("angstrom is a length unit");
    for record in &syntax.atoms {
        let atom_id = mol.add_atom(record.atom.clone());
        conformer
            .set_position(atom_id, Quantity::new(record.point, ANGSTROM))
            .expect("matching coordinate units");
        atom_ids.insert(record.index, atom_id);
    }
    for record in &syntax.bonds {
        let a = *atom_ids.get(&record.a).ok_or_else(|| {
            SdfParseError::new(1, record.line, "bond endpoint outside parsed atom records")
        })?;
        let b = *atom_ids.get(&record.b).ok_or_else(|| {
            SdfParseError::new(1, record.line, "bond endpoint outside parsed atom records")
        })?;
        let bond_id = mol
            .add_bond(a, b, record.order)
            .map_err(|error| SdfParseError::new(1, record.line, error.to_string()))?;
        if let Some(kind) = record.stereo {
            mol.set_stereo_bond_mark(StereoBondMark {
                bond: bond_id,
                kind,
                source: StereoSource::MolfileV3000,
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
    index: usize,
    order: BondOrder,
    a: usize,
    b: usize,
    stereo: Option<StereoBondMarkKind>,
}

fn collect_v3000_lines(
    record: usize,
    start_line: usize,
    lines: &[&str],
    options: MolfileParseOptions,
) -> std::result::Result<Vec<V3000Line>, SdfParseError> {
    let mut records = Vec::new();
    let mut index = 4usize;
    while index < lines.len() {
        let line_number = checked_line_number(record, start_line, index)?;
        let line = lines[index];
        if line.trim() == "M  END" {
            return Ok(records);
        }
        let body = v3000_body(line)
            .ok_or_else(|| SdfParseError::new(record, line_number, "expected M  V30 record"))?;
        if body.len() > options.max_v3000_logical_line_bytes {
            return Err(SdfParseError::new(
                record,
                line_number,
                "V3000 logical line exceeds configured byte limit",
            ));
        }
        let mut body = body.to_owned();
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
            let continuation = continuation.trim_start();
            let next_len = body.len().checked_add(continuation.len()).ok_or_else(|| {
                SdfParseError::new(record, line_number, "V3000 continuation length overflow")
            })?;
            if next_len > options.max_v3000_logical_line_bytes {
                return Err(SdfParseError::new(
                    record,
                    line_number,
                    "V3000 logical line exceeds configured byte limit",
                ));
            }
            body.push_str(continuation);
        }
        records.push(V3000Line {
            line: line_number,
            body,
        });
        index += 1;
    }
    Err(SdfParseError::new(record, start_line, "missing M  END"))
}

fn v3000_body(line: &str) -> Option<&str> {
    let trimmed = line.strip_prefix("M  V30 ")?;
    Some(trimmed.trim())
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

fn parse_v3000_counts(line: &str) -> Option<V3000Counts> {
    let mut fields = line.split_whitespace();
    (fields.next()? == "COUNTS").then_some(())?;
    let counts = V3000Counts {
        atoms: fields.next()?.parse().ok()?,
        bonds: fields.next()?.parse().ok()?,
    };
    let _sgroups = fields.next()?.parse::<usize>().ok()?;
    let _three_dimensional_constraints = fields.next()?.parse::<usize>().ok()?;
    let chiral = fields.next()?.parse::<u8>().ok()?;
    if chiral > 1 {
        return None;
    }
    if let Some(regno) = fields.next() {
        regno
            .strip_prefix("REGNO=")?
            .parse::<u64>()
            .ok()
            .map(|_| ())?;
    }
    fields.next().is_none().then_some(counts)
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
    let options = fields.map(split_v3000_option).collect::<Option<Vec<_>>>()?;
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
    let mut seen = std::collections::BTreeSet::new();
    for (key, value) in options {
        if !seen.insert(*key) {
            return Err(SdfParseError::new(
                record,
                line,
                format!("duplicate V3000 atom option `{key}`"),
            ));
        }
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
            _ => {
                return Err(SdfParseError::new(
                    record,
                    line,
                    format!("unsupported V3000 atom option `{key}`"),
                ))
            }
        }
    }
    Ok(())
}

fn parse_v3000_bond(
    record: usize,
    line_number: usize,
    line: &str,
) -> std::result::Result<V3000Bond, SdfParseError> {
    let invalid = || SdfParseError::new(record, line_number, "invalid V3000 bond line");
    let mut fields = line.split_whitespace();
    let index = fields
        .next()
        .ok_or_else(invalid)?
        .parse::<usize>()
        .map_err(|_| invalid())?;
    let order_code = fields
        .next()
        .ok_or_else(invalid)?
        .parse()
        .map_err(|_| invalid())?;
    let a = fields
        .next()
        .ok_or_else(invalid)?
        .parse()
        .map_err(|_| invalid())?;
    let b = fields
        .next()
        .ok_or_else(invalid)?
        .parse()
        .map_err(|_| invalid())?;
    let order = v3000_bond_order(order_code).ok_or_else(invalid)?;
    let mut stereo = None;
    let mut seen = std::collections::BTreeSet::new();
    for field in fields {
        let (key, value) = split_v3000_option(field).ok_or_else(invalid)?;
        if !seen.insert(key) {
            return Err(SdfParseError::new(
                record,
                line_number,
                format!("duplicate V3000 bond option `{key}`"),
            ));
        }
        if key != "CFG" {
            return Err(SdfParseError::new(
                record,
                line_number,
                format!("unsupported V3000 bond option `{key}`"),
            ));
        }
        stereo = v3000_bond_stereo(order, value).ok_or_else(|| {
            SdfParseError::new(record, line_number, "unsupported V3000 bond CFG value")
        })?;
    }
    Ok(V3000Bond {
        index,
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

fn v3000_bond_stereo(order: BondOrder, value: &str) -> Option<Option<StereoBondMarkKind>> {
    match (order, value) {
        (_, "0") => Some(None),
        (BondOrder::Single, "1") => Some(Some(StereoBondMarkKind::WedgeUp)),
        (BondOrder::Single, "2") => Some(Some(StereoBondMarkKind::WedgeEither)),
        (BondOrder::Single, "3") => Some(Some(StereoBondMarkKind::WedgeDown)),
        (BondOrder::Double, "2") => Some(Some(StereoBondMarkKind::DoubleBondEither)),
        _ => None,
    }
}

fn v3000_bond_code(order: BondOrder) -> std::result::Result<u8, MolWriteError> {
    match order {
        BondOrder::Zero => Ok(0),
        BondOrder::Single => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        BondOrder::Dative => Ok(9),
        BondOrder::Quadruple => Err(MolWriteError::new(
            "V3000 writer does not support quadruple bonds",
        )),
    }
}

fn v3000_bond_cfg(
    order: BondOrder,
    stereo: Option<&StereoBondMark>,
) -> std::result::Result<Option<u8>, MolWriteError> {
    match (order, stereo.map(|mark| mark.kind)) {
        (_, None) => Ok(None),
        (BondOrder::Single, Some(StereoBondMarkKind::WedgeUp)) => Ok(Some(1)),
        (BondOrder::Single, Some(StereoBondMarkKind::WedgeEither)) => Ok(Some(2)),
        (BondOrder::Single, Some(StereoBondMarkKind::WedgeDown)) => Ok(Some(3)),
        (BondOrder::Double, Some(StereoBondMarkKind::DoubleBondEither)) => Ok(Some(2)),
        _ => Err(MolWriteError::new(
            "V3000 bond CFG is incompatible with the bond order",
        )),
    }
}

fn v3000_radical_code(radical: AtomRadical) -> std::result::Result<u8, MolWriteError> {
    match radical {
        AtomRadical::Singlet => Ok(1),
        AtomRadical::Doublet => Ok(2),
        AtomRadical::Triplet => Ok(3),
        AtomRadical::Quartet | AtomRadical::Quintet => Err(MolWriteError::new(
            "V3000 writer cannot encode radical multiplicity above triplet",
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

fn checked_line_number(
    record: usize,
    start_line: usize,
    offset: usize,
) -> std::result::Result<usize, SdfParseError> {
    start_line
        .checked_add(offset)
        .ok_or_else(|| SdfParseError::new(record, start_line, "line number overflow"))
}
