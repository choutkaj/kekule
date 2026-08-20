use std::collections::{BTreeMap, BTreeSet};

use crate::algorithms::explicit_valence;
use crate::chemistry::{
    localize_source_aromatic_bonds, project_molfile_stereo_bond_marks, SourceStereoBondMark,
    SourceStereoBondMarkKind,
};
use crate::core::*;
use crate::geometry::Point3;
use crate::io::{MolWriteError, MolfileParseOptions, MolfileVersion, SdfParseError};
use crate::small::model::SmallMolecule;
use crate::units::{Quantity, ANGSTROM};

use super::structure_documents::{
    apply_molfile_declared_valence, checked_line_number, interpret_molfile_atom_fields,
};

#[derive(Debug, Clone, PartialEq)]
pub(super) struct V3000Syntax {
    pub(super) atoms: Vec<V3000AtomSyntax>,
    pub(super) bonds: Vec<V3000BondSyntax>,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct V3000AtomSyntax {
    pub(super) index: usize,
    pub(super) symbol: String,
    pub(super) formal_charge: i32,
    pub(super) isotope: Option<i32>,
    pub(super) radical: Option<V3000RadicalSyntax>,
    pub(super) hydrogen_count: Option<i32>,
    pub(super) valence: Option<i32>,
    pub(super) atom_map: Option<u32>,
    pub(super) unsupported_options: Vec<String>,
    pub(super) coordinates: [f64; 3],
    pub(super) line: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum V3000RadicalSyntax {
    Singlet,
    Doublet,
    Triplet,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct V3000BondSyntax {
    pub(super) a: usize,
    pub(super) b: usize,
    pub(super) order_code: u8,
    pub(super) stereo_code: Option<u8>,
    pub(super) unsupported_options: Vec<String>,
    pub(super) line: usize,
}

pub fn write_mol_v3000(molecule: &SmallMolecule) -> std::result::Result<String, MolWriteError> {
    let mol = molecule.graph();
    let atoms = mol.atom_ids().collect::<Vec<_>>();
    let bonds = mol.bond_ids().collect::<Vec<_>>();
    let mut atom_index = BTreeMap::new();
    for (serial, atom_id) in (1u64..).zip(atoms.iter()) {
        atom_index.insert(*atom_id, serial);
    }

    let title = "";
    let program = "kekule";
    let comment = "";
    let conformer = mol.first_conformer().map(|(_, conformer)| conformer);
    let projected_stereo = project_molfile_stereo_bond_marks(mol).map_err(MolWriteError::new)?;

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
        match atom.hydrogens {
            HydrogenDeclaration::Infer { explicit: 0 } => {}
            HydrogenDeclaration::Infer { .. } => {
                return Err(MolWriteError::new(format!(
                    "V3000 cannot encode represented hydrogens while leaving implicit-H inference enabled for atom {}",
                    atom_id.index()
                )));
            }
            HydrogenDeclaration::Fixed(explicit) => {
                if explicit > 0 {
                    out.push_str(&format!(" HCOUNT={explicit}"));
                } else {
                    out.push_str(&format!(" VAL={}", explicit_valence(mol, *atom_id)));
                }
            }
        }
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
    out.push_str("M  V30 BEGIN BOND\n");
    for (serial, bond_id) in (1u64..).zip(bonds.iter()) {
        let bond = mol
            .bond(*bond_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let projection = projected_stereo.get(bond_id).copied();
        let (from, to) = projection
            .map(|projection| (projection.from, bond.other_atom(projection.from)))
            .unwrap_or_else(|| bond.endpoints());
        let a = atom_index
            .get(&from)
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from V3000 atom table"))?;
        let b = atom_index
            .get(&to)
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from V3000 atom table"))?;
        let order_code = v3000_bond_code(bond.order)?;
        out.push_str(&format!("M  V30 {serial} {order_code} {a} {b}"));
        let stereo = projection.map(|projection| projection.kind);
        if let Some(cfg) = v3000_bond_cfg(bond.order, stereo)? {
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
        let mut atom = V3000AtomSyntax {
            index: parsed.index,
            symbol: parsed.symbol.to_owned(),
            formal_charge: 0,
            isotope: None,
            radical: None,
            hydrogen_count: None,
            valence: None,
            atom_map: (parsed.atom_map != 0).then_some(parsed.atom_map),
            unsupported_options: Vec::new(),
            coordinates: parsed.coordinates,
            line: row.line,
        };
        apply_v3000_atom_options(record, row.line, &mut atom, &parsed.options)?;
        atom_indices.insert(parsed.index, atoms.len());
        atoms.push(atom);
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
            order_code: parsed.order_code,
            stereo_code: parsed.stereo_code,
            unsupported_options: parsed.unsupported_options,
            line: row.line,
        });
    }

    Ok(V3000Syntax { atoms, bonds })
}

pub(super) fn interpret_v3000_syntax(
    syntax: &V3000Syntax,
) -> std::result::Result<(SmallMolecule, Vec<SourceStereoBondMark>), SdfParseError> {
    let mut mol = Molecule::new();
    let mut atom_ids = BTreeMap::<usize, AtomId>::new();
    let mut conformer = Conformer::with_atom_capacity(syntax.atoms.len(), ANGSTROM)
        .expect("angstrom is a length unit");
    for record in &syntax.atoms {
        let atom = interpret_v3000_atom(record)?;
        let atom_id = mol.add_atom(atom).map_err(|error| {
            SdfParseError::new(1, record.line, format!("invalid graph atom: {error}"))
        })?;
        conformer
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
            .expect("matching coordinate units");
        atom_ids.insert(record.index, atom_id);
    }
    let mut source_aromatic_bonds = BTreeSet::new();
    let mut source_stereo = Vec::new();
    let mut first_aromatic_line = None;
    for record in &syntax.bonds {
        let a = *atom_ids.get(&record.a).ok_or_else(|| {
            SdfParseError::new(1, record.line, "bond endpoint outside parsed atom records")
        })?;
        let b = *atom_ids.get(&record.b).ok_or_else(|| {
            SdfParseError::new(1, record.line, "bond endpoint outside parsed atom records")
        })?;
        if let Some(option) = record.unsupported_options.first() {
            return Err(SdfParseError::new(
                1,
                record.line,
                format!("unsupported V3000 bond option `{option}`"),
            ));
        }
        let (order, source_aromatic) = interpret_v3000_bond_order(record.order_code, record.line)?;
        let bond_id = mol
            .add_bond(a, b, order)
            .map_err(|error| SdfParseError::new(1, record.line, error.to_string()))?;
        if source_aromatic {
            source_aromatic_bonds.insert(bond_id);
            first_aromatic_line.get_or_insert(record.line);
        }
        if let Some(kind) =
            interpret_v3000_bond_stereo(order, source_aromatic, record.stereo_code, record.line)?
        {
            source_stereo.push(SourceStereoBondMark {
                bond: bond_id,
                from: a,
                kind,
            });
        }
    }

    apply_v3000_declared_hydrogens(&mut mol, syntax, &atom_ids, &source_aromatic_bonds)?;
    localize_source_aromatic_bonds(&mut mol, &source_aromatic_bonds).map_err(|error| {
        SdfParseError::new(
            1,
            first_aromatic_line.unwrap_or(1),
            format!("aromatic bond localization failed: {error}"),
        )
    })?;

    if conformer.positions().next().is_some() {
        mol.add_conformer(conformer)
            .expect("parsed coordinates reference live atoms");
    }
    Ok((
        SmallMolecule::from_graph_unchecked_connectedness(mol),
        source_stereo,
    ))
}

fn interpret_v3000_atom(record: &V3000AtomSyntax) -> std::result::Result<Atom, SdfParseError> {
    if let Some(option) = record.unsupported_options.first() {
        return Err(SdfParseError::new(
            1,
            record.line,
            format!("unsupported V3000 atom option `{option}`"),
        ));
    }
    let radical = record.radical.map(|radical| match radical {
        V3000RadicalSyntax::Singlet => AtomRadical::Singlet,
        V3000RadicalSyntax::Doublet => AtomRadical::Doublet,
        V3000RadicalSyntax::Triplet => AtomRadical::Triplet,
    });
    let explicit = interpret_v3000_count_declaration(record.hydrogen_count, "HCOUNT", record.line)?
        .unwrap_or(0);
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

fn interpret_v3000_bond_order(
    code: u8,
    line: usize,
) -> std::result::Result<(BondOrder, bool), SdfParseError> {
    match code {
        4 => Ok((BondOrder::Single, true)),
        _ => v3000_bond_order(code)
            .map(|order| (order, false))
            .ok_or_else(|| {
                SdfParseError::new(1, line, format!("unsupported V3000 bond order code {code}"))
            }),
    }
}

fn interpret_v3000_bond_stereo(
    order: BondOrder,
    source_aromatic: bool,
    code: Option<u8>,
    line: usize,
) -> std::result::Result<Option<SourceStereoBondMarkKind>, SdfParseError> {
    let Some(code) = code else {
        return Ok(None);
    };
    if source_aromatic {
        return Err(SdfParseError::new(
            1,
            line,
            format!("V3000 bond CFG value {code} is unsupported for an aromatic source bond"),
        ));
    }
    v3000_bond_stereo(order, &code.to_string()).ok_or_else(|| {
        SdfParseError::new(1, line, format!("unsupported V3000 bond CFG value {code}"))
    })
}

fn apply_v3000_declared_hydrogens(
    molecule: &mut Molecule,
    syntax: &V3000Syntax,
    atom_ids: &BTreeMap<usize, AtomId>,
    source_aromatic_bonds: &BTreeSet<BondId>,
) -> std::result::Result<(), SdfParseError> {
    for record in &syntax.atoms {
        let Some(declared_valence) =
            interpret_v3000_count_declaration(record.valence, "VAL", record.line)?
        else {
            continue;
        };
        let atom_id = *atom_ids.get(&record.index).ok_or_else(|| {
            SdfParseError::new(1, record.line, "V3000 atom index was not interpreted")
        })?;
        let declared_hydrogens = record.hydrogen_count.map(|_| {
            molecule
                .atom(atom_id)
                .expect("interpreted V3000 atom remains live")
                .hydrogens
                .explicit_count()
        });
        apply_molfile_declared_valence(
            molecule,
            atom_id,
            declared_valence,
            declared_hydrogens,
            source_aromatic_bonds,
            MolfileVersion::V3000,
            record.line,
        )?;
    }
    Ok(())
}

fn interpret_v3000_count_declaration(
    value: Option<i32>,
    name: &str,
    line: usize,
) -> std::result::Result<Option<u8>, SdfParseError> {
    value
        .map(|value| {
            let value = if value == -1 { 0 } else { value };
            u8::try_from(value).map_err(|_| {
                SdfParseError::new(
                    1,
                    line,
                    format!("V3000 {name} declaration is outside the represented count range"),
                )
            })
        })
        .transpose()
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
    coordinates: [f64; 3],
    atom_map: u32,
    options: Vec<(&'a str, &'a str)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct V3000Bond {
    index: usize,
    order_code: u8,
    a: usize,
    b: usize,
    stereo_code: Option<u8>,
    unsupported_options: Vec<String>,
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
    let coordinates = [
        parse_finite_f64(fields.next()?)?,
        parse_finite_f64(fields.next()?)?,
        parse_finite_f64(fields.next()?)?,
    ];
    let atom_map = fields.next()?.parse().ok()?;
    let options = fields.map(split_v3000_option).collect::<Option<Vec<_>>>()?;
    Some(V3000Atom {
        index,
        symbol,
        coordinates,
        atom_map,
        options,
    })
}

fn apply_v3000_atom_options(
    record: usize,
    line: usize,
    atom: &mut V3000AtomSyntax,
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
                    .parse::<i32>()
                    .map_err(|_| SdfParseError::new(record, line, "invalid V3000 MASS value"))?;
                atom.isotope = (isotope != 0).then_some(isotope);
            }
            "RAD" => {
                atom.radical = Some(match *value {
                    "1" => V3000RadicalSyntax::Singlet,
                    "2" => V3000RadicalSyntax::Doublet,
                    "3" => V3000RadicalSyntax::Triplet,
                    _ => {
                        return Err(SdfParseError::new(
                            record,
                            line,
                            "unsupported V3000 RAD code",
                        ))
                    }
                });
            }
            "HCOUNT" => {
                atom.hydrogen_count = parse_v3000_count_declaration(record, line, value, "HCOUNT")?;
            }
            "VAL" => {
                atom.valence = parse_v3000_count_declaration(record, line, value, "VAL")?;
            }
            _ => atom.unsupported_options.push((*key).to_owned()),
        }
    }
    Ok(())
}

fn parse_v3000_count_declaration(
    record: usize,
    line: usize,
    value: &str,
    name: &str,
) -> std::result::Result<Option<i32>, SdfParseError> {
    let value = value
        .parse::<i32>()
        .map_err(|_| SdfParseError::new(record, line, format!("invalid V3000 {name} value")))?;
    if value < -1 {
        return Err(SdfParseError::new(
            record,
            line,
            format!("invalid V3000 {name} value"),
        ));
    }
    Ok((value != 0).then_some(value))
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
    let mut stereo_code = None;
    let mut unsupported_options = Vec::new();
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
        if key == "CFG" {
            stereo_code = Some(value.parse::<u8>().map_err(|_| invalid())?);
        } else {
            unsupported_options.push(key.to_owned());
        }
    }
    Ok(V3000Bond {
        index,
        order_code,
        a,
        b,
        stereo_code,
        unsupported_options,
    })
}

fn v3000_bond_order(code: u8) -> Option<BondOrder> {
    match code {
        0 => Some(BondOrder::Zero),
        1 => Some(BondOrder::Single),
        2 => Some(BondOrder::Double),
        3 => Some(BondOrder::Triple),
        9 => Some(BondOrder::Dative),
        _ => None,
    }
}

fn v3000_bond_stereo(order: BondOrder, value: &str) -> Option<Option<SourceStereoBondMarkKind>> {
    match (order, value) {
        (_, "0") => Some(None),
        (BondOrder::Single, "1") => Some(Some(SourceStereoBondMarkKind::WedgeUp)),
        (BondOrder::Single, "2") => Some(Some(SourceStereoBondMarkKind::WedgeEither)),
        (BondOrder::Single, "3") => Some(Some(SourceStereoBondMarkKind::WedgeDown)),
        (BondOrder::Double, "2") => Some(Some(SourceStereoBondMarkKind::DoubleBondEither)),
        _ => None,
    }
}

fn v3000_bond_code(order: BondOrder) -> std::result::Result<u8, MolWriteError> {
    match order {
        BondOrder::Zero => Ok(0),
        BondOrder::Single => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Dative => Ok(9),
        BondOrder::Quadruple => Err(MolWriteError::new(
            "V3000 writer does not support quadruple bonds",
        )),
    }
}

fn v3000_bond_cfg(
    order: BondOrder,
    stereo: Option<SourceStereoBondMarkKind>,
) -> std::result::Result<Option<u8>, MolWriteError> {
    match (order, stereo) {
        (_, None) => Ok(None),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeUp)) => Ok(Some(1)),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeEither)) => Ok(Some(2)),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeDown)) => Ok(Some(3)),
        (BondOrder::Double, Some(SourceStereoBondMarkKind::DoubleBondEither)) => Ok(Some(2)),
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
