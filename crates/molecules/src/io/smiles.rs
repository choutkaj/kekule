use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::algorithms::{canonical_atom_ranking, ordered_atom_pair, CanonicalAtomRanking};
use crate::core::*;
use crate::io::MolWriteError;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SmilesParseOptions;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SmilesWriteOptions;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CanonicalSmilesWriteOptions;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesParseError {
    pub offset: usize,
    pub message: String,
}

impl SmilesParseError {
    fn new(offset: usize, message: impl Into<String>) -> Self {
        Self {
            offset,
            message: message.into(),
        }
    }
}

impl fmt::Display for SmilesParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SMILES parse error at {}: {}", self.offset, self.message)
    }
}

impl std::error::Error for SmilesParseError {}

pub fn read_smiles_str(
    input: &str,
    _options: SmilesParseOptions,
) -> std::result::Result<SmallMolecule, SmilesParseError> {
    let chars = input.char_indices().collect::<Vec<_>>();
    let mut mol = Molecule::new();
    let mut current: Option<AtomId> = None;
    let mut stack = Vec::<AtomId>::new();
    let mut pending_bond = None::<(BondOrder, usize)>;
    let mut rings = BTreeMap::<usize, (AtomId, Option<BondOrder>, usize)>::new();
    let mut component = 0usize;
    let mut previous = SmilesTokenKind::Start;
    let mut cursor = 0;
    while cursor < chars.len() {
        let (offset, ch) = chars[cursor];
        match ch {
            '(' => {
                if !matches!(
                    previous,
                    SmilesTokenKind::Atom | SmilesTokenKind::Ring | SmilesTokenKind::BranchClose
                ) || pending_bond.is_some()
                {
                    return Err(SmilesParseError::new(offset, "invalid branch start"));
                }
                let atom =
                    current.ok_or_else(|| SmilesParseError::new(offset, "branch without atom"))?;
                stack.push(atom);
                previous = SmilesTokenKind::BranchOpen;
                cursor += 1;
            }
            ')' => {
                if matches!(
                    previous,
                    SmilesTokenKind::Start
                        | SmilesTokenKind::BranchOpen
                        | SmilesTokenKind::Bond
                        | SmilesTokenKind::Dot
                ) {
                    return Err(SmilesParseError::new(offset, "empty or incomplete branch"));
                }
                current = Some(
                    stack
                        .pop()
                        .ok_or_else(|| SmilesParseError::new(offset, "unmatched branch close"))?,
                );
                previous = SmilesTokenKind::BranchClose;
                cursor += 1;
            }
            '.' => {
                if current.is_none()
                    || pending_bond.is_some()
                    || !stack.is_empty()
                    || matches!(
                        previous,
                        SmilesTokenKind::Start
                            | SmilesTokenKind::BranchOpen
                            | SmilesTokenKind::Bond
                            | SmilesTokenKind::Dot
                    )
                {
                    return Err(SmilesParseError::new(offset, "invalid component separator"));
                }
                current = None;
                component = component
                    .checked_add(1)
                    .ok_or_else(|| SmilesParseError::new(offset, "component counter overflow"))?;
                previous = SmilesTokenKind::Dot;
                cursor += 1;
            }
            '-' | '=' | '#' | ':' => {
                if current.is_none()
                    || pending_bond.is_some()
                    || !matches!(
                        previous,
                        SmilesTokenKind::Atom
                            | SmilesTokenKind::Ring
                            | SmilesTokenKind::BranchClose
                            | SmilesTokenKind::BranchOpen
                    )
                {
                    return Err(SmilesParseError::new(offset, "bond without left endpoint"));
                }
                let order = match ch {
                    '-' => BondOrder::Single,
                    '=' => BondOrder::Double,
                    '#' => BondOrder::Triple,
                    ':' => BondOrder::Aromatic,
                    _ => unreachable!(),
                };
                pending_bond = Some((order, offset));
                previous = SmilesTokenKind::Bond;
                cursor += 1;
            }
            '0'..='9' | '%' => {
                let atom = current
                    .ok_or_else(|| SmilesParseError::new(offset, "ring closure without atom"))?;
                let (label, next_cursor) = parse_smiles_ring_label(&chars, cursor)?;
                let close_order = pending_bond.take().map(|(order, _)| order);
                if let Some((other, open_order, open_component)) = rings.remove(&label) {
                    if open_component != component {
                        return Err(SmilesParseError::new(
                            offset,
                            "ring closure crosses a component separator",
                        ));
                    }
                    if open_order.is_some() && close_order.is_some() && open_order != close_order {
                        return Err(SmilesParseError::new(
                            offset,
                            "conflicting ring bond symbols",
                        ));
                    }
                    let order = match close_order.or(open_order) {
                        Some(order) => order,
                        None => default_smiles_bond_order(&mol, other, atom, offset)?,
                    };
                    mol.add_bond(other, atom, order)
                        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
                } else {
                    rings.insert(label, (atom, close_order, component));
                }
                previous = SmilesTokenKind::Ring;
                cursor = next_cursor;
            }
            '[' => {
                let (atom, next_cursor) = parse_bracket_atom(&chars, cursor)?;
                let atom_id = mol.add_atom(atom);
                if let Some(previous) = current {
                    let order = match pending_bond.take().map(|(order, _)| order) {
                        Some(order) => order,
                        None => default_smiles_bond_order(&mol, previous, atom_id, offset)?,
                    };
                    mol.add_bond(previous, atom_id, order)
                        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
                } else if pending_bond.is_some() {
                    return Err(SmilesParseError::new(offset, "bond without left endpoint"));
                }
                current = Some(atom_id);
                previous = SmilesTokenKind::Atom;
                cursor = next_cursor;
            }
            '@' | '/' | '\\' | '*' => {
                return Err(SmilesParseError::new(
                    offset,
                    "unsupported stereochemistry or query syntax",
                ));
            }
            _ => {
                let (atom, next_cursor) = parse_organic_atom(&chars, cursor)?;
                let atom_id = mol.add_atom(atom);
                if let Some(previous) = current {
                    let order = match pending_bond.take().map(|(order, _)| order) {
                        Some(order) => order,
                        None => default_smiles_bond_order(&mol, previous, atom_id, offset)?,
                    };
                    mol.add_bond(previous, atom_id, order)
                        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
                } else if pending_bond.is_some() {
                    return Err(SmilesParseError::new(offset, "bond without left endpoint"));
                }
                current = Some(atom_id);
                previous = SmilesTokenKind::Atom;
                cursor = next_cursor;
            }
        }
    }
    if !stack.is_empty() {
        return Err(SmilesParseError::new(input.len(), "unclosed branch"));
    }
    if !rings.is_empty() {
        return Err(SmilesParseError::new(input.len(), "unclosed ring closure"));
    }
    if let Some((_, offset)) = pending_bond {
        return Err(SmilesParseError::new(offset, "bond without right endpoint"));
    }
    if matches!(previous, SmilesTokenKind::Dot | SmilesTokenKind::BranchOpen) {
        return Err(SmilesParseError::new(input.len(), "incomplete SMILES"));
    }
    Ok(SmallMolecule { mol })
}

fn default_smiles_bond_order(
    mol: &Molecule,
    left: AtomId,
    right: AtomId,
    offset: usize,
) -> std::result::Result<BondOrder, SmilesParseError> {
    let left = mol
        .atom(left)
        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
    let right = mol
        .atom(right)
        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
    if left.aromatic && right.aromatic {
        Ok(BondOrder::Aromatic)
    } else {
        Ok(BondOrder::Single)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SmilesTokenKind {
    Start,
    Atom,
    BranchOpen,
    BranchClose,
    Bond,
    Ring,
    Dot,
}

fn parse_smiles_ring_label(
    chars: &[(usize, char)],
    cursor: usize,
) -> std::result::Result<(usize, usize), SmilesParseError> {
    let (offset, ch) = chars[cursor];
    if ch != '%' {
        return Ok(((ch as u8 - b'0') as usize, cursor + 1));
    }
    let first = chars
        .get(cursor + 1)
        .filter(|(_, ch)| ch.is_ascii_digit())
        .ok_or_else(|| SmilesParseError::new(offset, "malformed percent ring label"))?;
    let second = chars
        .get(cursor + 2)
        .filter(|(_, ch)| ch.is_ascii_digit())
        .ok_or_else(|| SmilesParseError::new(offset, "malformed percent ring label"))?;
    let label = first.1.to_digit(10).unwrap_or(0) as usize * 10
        + second.1.to_digit(10).unwrap_or(0) as usize;
    if label < 10 {
        return Err(SmilesParseError::new(
            offset,
            "percent ring labels must be between 10 and 99",
        ));
    }
    Ok((label, cursor + 3))
}

fn parse_organic_atom(
    chars: &[(usize, char)],
    cursor: usize,
) -> std::result::Result<(Atom, usize), SmilesParseError> {
    let (offset, ch) = chars[cursor];
    let mut symbol = ch.to_string();
    let mut aromatic = false;
    let mut next = cursor + 1;
    let following = chars.get(cursor + 1).map(|(_, c)| *c);
    if (ch == 'C' && following == Some('l')) || (ch == 'B' && following == Some('r')) {
        symbol.push(chars[cursor + 1].1);
        next += 1;
    } else if matches!(ch, 'b' | 'c' | 'n' | 'o' | 'p' | 's') {
        symbol = ch.to_ascii_uppercase().to_string();
        aromatic = true;
    } else if !matches!(ch, 'B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I') {
        return Err(SmilesParseError::new(
            offset,
            format!("unsupported organic-subset atom `{ch}`"),
        ));
    }
    let element = Element::from_symbol(&symbol)
        .ok_or_else(|| SmilesParseError::new(offset, format!("unsupported atom `{ch}`")))?;
    let mut atom = Atom::new(element);
    atom.aromatic = aromatic;
    Ok((atom, next))
}

fn parse_bracket_atom(
    chars: &[(usize, char)],
    cursor: usize,
) -> std::result::Result<(Atom, usize), SmilesParseError> {
    let start = chars[cursor].0;
    let mut end = cursor + 1;
    while end < chars.len() && chars[end].1 != ']' {
        end += 1;
    }
    if end == chars.len() {
        return Err(SmilesParseError::new(start, "unclosed bracket atom"));
    }
    let text = chars[cursor + 1..end]
        .iter()
        .map(|(_, c)| *c)
        .collect::<String>();
    if text.is_empty() {
        return Err(SmilesParseError::new(start, "empty bracket atom"));
    }
    if !text.is_ascii() {
        return Err(SmilesParseError::new(
            start,
            "bracket atom must use ASCII syntax",
        ));
    }
    let bytes = text.as_bytes();
    let mut index = 0;
    let isotope_end = ascii_digits_end(bytes, index);
    let isotope = if isotope_end > index {
        let value = text[index..isotope_end]
            .parse::<u16>()
            .map_err(|_| SmilesParseError::new(start + 1 + index, "invalid isotope"))?;
        if value == 0 {
            return Err(SmilesParseError::new(
                start + 1 + index,
                "isotope must be positive",
            ));
        }
        index = isotope_end;
        Some(value)
    } else {
        None
    };
    let symbol_start = index;
    let first = *bytes
        .get(index)
        .ok_or_else(|| SmilesParseError::new(start, "bracket atom missing element"))?;
    let aromatic = first.is_ascii_lowercase();
    let canonical_symbol = if aromatic {
        let symbol = match first {
            b'b' => "B",
            b'c' => "C",
            b'n' => "N",
            b'o' => "O",
            b'p' => "P",
            b's' => "S",
            _ => {
                return Err(SmilesParseError::new(
                    start + 1 + index,
                    "unsupported aromatic bracket element",
                ))
            }
        };
        index += 1;
        symbol.to_owned()
    } else if first.is_ascii_uppercase() {
        index += 1;
        if bytes.get(index).is_some_and(u8::is_ascii_lowercase) {
            index += 1;
        }
        text[symbol_start..index].to_owned()
    } else {
        return Err(SmilesParseError::new(
            start + 1 + index,
            "bracket atom missing element",
        ));
    };
    let element = Element::from_symbol(&canonical_symbol).ok_or_else(|| {
        SmilesParseError::new(start + 1 + symbol_start, "unsupported bracket element")
    })?;
    let mut atom = Atom::new(element);
    atom.aromatic = aromatic;
    atom.isotope = isotope;
    atom.no_implicit_hydrogens = true;
    let mut saw_hydrogen = false;
    let mut saw_charge = false;
    let mut saw_map = false;
    while index < text.len() {
        match bytes[index] {
            b'H' if !saw_hydrogen && !saw_charge && !saw_map => {
                saw_hydrogen = true;
                index += 1;
                let digit_end = ascii_digits_end(bytes, index);
                atom.explicit_hydrogens = if digit_end == index {
                    1
                } else {
                    let value = text[index..digit_end].parse::<u8>().map_err(|_| {
                        SmilesParseError::new(start + 1 + index, "invalid hydrogen count")
                    })?;
                    if value == 0 {
                        return Err(SmilesParseError::new(
                            start + 1 + index,
                            "hydrogen count must be positive",
                        ));
                    }
                    index = digit_end;
                    value
                };
            }
            b'+' | b'-' if !saw_charge && !saw_map => {
                saw_charge = true;
                let sign_byte = bytes[index];
                let sign = if sign_byte == b'+' { 1i16 } else { -1i16 };
                index += 1;
                let mut magnitude = 1u16;
                while bytes.get(index) == Some(&sign_byte) {
                    magnitude = magnitude.checked_add(1).ok_or_else(|| {
                        SmilesParseError::new(start + 1 + index, "charge overflow")
                    })?;
                    index += 1;
                }
                let digit_end = ascii_digits_end(bytes, index);
                if digit_end > index {
                    if magnitude != 1 {
                        return Err(SmilesParseError::new(
                            start + 1 + index,
                            "charge cannot mix repeated signs and digits",
                        ));
                    }
                    magnitude = text[index..digit_end]
                        .parse::<u16>()
                        .map_err(|_| SmilesParseError::new(start + 1 + index, "invalid charge"))?;
                    if magnitude == 0 {
                        return Err(SmilesParseError::new(
                            start + 1 + index,
                            "charge magnitude must be positive",
                        ));
                    }
                    index = digit_end;
                }
                let charge =
                    sign.checked_mul(i16::try_from(magnitude).map_err(|_| {
                        SmilesParseError::new(start + 1 + index, "charge overflow")
                    })?)
                    .ok_or_else(|| SmilesParseError::new(start + 1 + index, "charge overflow"))?;
                atom.formal_charge = i8::try_from(charge).map_err(|_| {
                    SmilesParseError::new(start + 1 + index, "charge is outside i8 range")
                })?;
            }
            b':' if !saw_map => {
                saw_map = true;
                index += 1;
                let digit_end = ascii_digits_end(bytes, index);
                if digit_end == index {
                    return Err(SmilesParseError::new(
                        start + 1 + index,
                        "atom map requires digits",
                    ));
                }
                let map = text[index..digit_end]
                    .parse::<u32>()
                    .map_err(|_| SmilesParseError::new(start + 1 + index, "invalid atom map"))?;
                if map == 0 {
                    return Err(SmilesParseError::new(
                        start + 1 + index,
                        "atom map must be positive",
                    ));
                }
                atom.atom_map = Some(map);
                index = digit_end;
            }
            b'@' | b'/' | b'\\' | b'*' => {
                return Err(SmilesParseError::new(
                    start + 1 + index,
                    "unsupported stereochemistry or query syntax",
                ));
            }
            _ => {
                return Err(SmilesParseError::new(
                    start + 1 + index,
                    "unsupported bracket atom syntax",
                ));
            }
        }
    }
    Ok((atom, end + 1))
}

fn ascii_digits_end(bytes: &[u8], mut index: usize) -> usize {
    while bytes.get(index).is_some_and(u8::is_ascii_digit) {
        index += 1;
    }
    index
}

pub fn write_smiles(
    molecule: &SmallMolecule,
    _options: SmilesWriteOptions,
) -> std::result::Result<String, MolWriteError> {
    let mol = &molecule.mol;
    let plan = plan_smiles_write(mol)?;
    let mut parts = Vec::new();
    for start in &plan.roots {
        parts.push(write_smiles_component(mol, *start, None, &plan)?);
    }
    Ok(parts.join("."))
}

pub fn write_canonical_smiles(
    molecule: &SmallMolecule,
    _options: CanonicalSmilesWriteOptions,
) -> std::result::Result<String, MolWriteError> {
    let mol = &molecule.mol;
    validate_smiles_writeable(mol)?;
    let ranking = canonical_atom_ranking(mol);
    let mut components = Vec::new();
    for component in smiles_connected_components(mol)? {
        let mut candidates = component
            .iter()
            .map(|root| write_canonical_smiles_component(mol, *root, &ranking))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        candidates.sort();
        if let Some(candidate) = candidates.into_iter().next() {
            components.push(candidate);
        }
    }
    components.sort();
    Ok(components.join("."))
}

#[derive(Debug, Clone)]
struct SmilesWritePlan {
    roots: Vec<AtomId>,
    tree_bonds: BTreeSet<BondId>,
    closures: BTreeMap<AtomId, Vec<SmilesRingClosure>>,
    subtree_sizes: BTreeMap<AtomId, usize>,
}

#[derive(Debug, Clone, Copy)]
struct SmilesRingClosure {
    number: usize,
    order: BondOrder,
    other: AtomId,
}

fn plan_smiles_write(mol: &Molecule) -> std::result::Result<SmilesWritePlan, MolWriteError> {
    validate_smiles_writeable(mol)?;
    let mut roots = Vec::new();
    let mut visited = BTreeSet::<AtomId>::new();
    let mut tree_bonds = BTreeSet::<BondId>::new();
    let mut ring_bonds = BTreeMap::<BondId, (AtomId, AtomId, BondOrder)>::new();

    for start in mol.atom_ids() {
        if visited.contains(&start) {
            continue;
        }
        roots.push(start);
        collect_smiles_tree(
            mol,
            start,
            None,
            &mut visited,
            &mut tree_bonds,
            &mut ring_bonds,
        )?;
    }

    let mut ring_bonds = ring_bonds
        .into_iter()
        .map(|(bond_id, (a, b, order))| {
            let (first, second) = ordered_atom_pair(a, b);
            (bond_id, first, second, order)
        })
        .collect::<Vec<_>>();
    ring_bonds.sort_by_key(|(bond_id, first, second, _)| (*first, *second, *bond_id));
    if ring_bonds.len() > 99 {
        return Err(MolWriteError::new(
            "SMILES writer supports at most 99 simultaneous ring closures",
        ));
    }

    let mut closures = BTreeMap::<AtomId, Vec<SmilesRingClosure>>::new();
    for (index, (_, first, second, order)) in ring_bonds.into_iter().enumerate() {
        let number = index + 1;
        closures.entry(first).or_default().push(SmilesRingClosure {
            number,
            order,
            other: second,
        });
        closures.entry(second).or_default().push(SmilesRingClosure {
            number,
            order,
            other: first,
        });
    }

    let mut subtree_sizes = BTreeMap::new();
    for root in &roots {
        compute_smiles_subtree_sizes(mol, *root, None, &tree_bonds, &mut subtree_sizes)?;
    }

    Ok(SmilesWritePlan {
        roots,
        tree_bonds,
        closures,
        subtree_sizes,
    })
}

fn validate_smiles_writeable(mol: &Molecule) -> std::result::Result<(), MolWriteError> {
    for (_, atom) in mol.atoms() {
        if atom.chiral.is_some() {
            return Err(MolWriteError::new(
                "SMILES writer cannot encode atom stereochemistry",
            ));
        }
        if atom.radical.is_some() {
            return Err(MolWriteError::new(
                "SMILES writer cannot encode atom radicals",
            ));
        }
    }
    for (_, bond) in mol.bonds() {
        match bond.order {
            BondOrder::Single | BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic => {}
            BondOrder::Zero | BondOrder::Dative | BondOrder::Quadruple => {
                return Err(MolWriteError::new(
                    "SMILES writer cannot encode zero, dative, or quadruple bonds",
                ));
            }
        }
        if !matches!(bond.stereo, None | Some(BondStereo::Unspecified)) {
            return Err(MolWriteError::new(
                "SMILES writer cannot encode bond stereochemistry",
            ));
        }
    }
    Ok(())
}

fn smiles_connected_components(
    mol: &Molecule,
) -> std::result::Result<Vec<Vec<AtomId>>, MolWriteError> {
    let mut components = Vec::new();
    let mut visited = BTreeSet::new();
    for start in mol.atom_ids() {
        if !visited.insert(start) {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![start];
        while let Some(atom) = stack.pop() {
            component.push(atom);
            for (_, _, neighbor) in smiles_incident_bonds(mol, atom)? {
                if visited.insert(neighbor) {
                    stack.push(neighbor);
                }
            }
        }
        component.sort();
        components.push(component);
    }
    Ok(components)
}

fn collect_smiles_tree(
    mol: &Molecule,
    atom_id: AtomId,
    parent_bond: Option<BondId>,
    visited: &mut BTreeSet<AtomId>,
    tree_bonds: &mut BTreeSet<BondId>,
    ring_bonds: &mut BTreeMap<BondId, (AtomId, AtomId, BondOrder)>,
) -> std::result::Result<(), MolWriteError> {
    struct Frame {
        parent_bond: Option<BondId>,
        incident: Vec<(BondId, BondOrder, AtomId)>,
        next_edge: usize,
    }

    visited.insert(atom_id);
    let mut stack = vec![Frame {
        parent_bond,
        incident: smiles_incident_bonds(mol, atom_id)?,
        next_edge: 0,
    }];
    while let Some(frame) = stack.last_mut() {
        if frame.next_edge >= frame.incident.len() {
            stack.pop();
            continue;
        }
        let (bond_id, order, neighbor) = frame.incident[frame.next_edge];
        frame.next_edge += 1;
        if Some(bond_id) == frame.parent_bond {
            continue;
        }
        if visited.contains(&neighbor) {
            if !tree_bonds.contains(&bond_id) {
                let bond = mol
                    .bond(bond_id)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                ring_bonds
                    .entry(bond_id)
                    .or_insert((bond.a(), bond.b(), order));
            }
            continue;
        }
        tree_bonds.insert(bond_id);
        visited.insert(neighbor);
        stack.push(Frame {
            parent_bond: Some(bond_id),
            incident: smiles_incident_bonds(mol, neighbor)?,
            next_edge: 0,
        });
    }
    Ok(())
}

fn write_canonical_smiles_component(
    mol: &Molecule,
    root: AtomId,
    ranking: &CanonicalAtomRanking,
) -> std::result::Result<String, MolWriteError> {
    let plan = plan_canonical_smiles_component(mol, root, ranking)?;
    write_canonical_smiles_component_with_plan(mol, root, &plan, ranking)
}

fn plan_canonical_smiles_component(
    mol: &Molecule,
    root: AtomId,
    ranking: &CanonicalAtomRanking,
) -> std::result::Result<SmilesWritePlan, MolWriteError> {
    struct Frame {
        parent_bond: Option<BondId>,
        incident: Vec<(BondId, BondOrder, AtomId)>,
        next_edge: usize,
    }

    let mut visited = BTreeSet::<AtomId>::new();
    let mut tree_bonds = BTreeSet::<BondId>::new();
    let mut ring_bonds = BTreeMap::<BondId, (AtomId, AtomId, BondOrder)>::new();
    visited.insert(root);
    let mut stack = vec![Frame {
        parent_bond: None,
        incident: canonical_smiles_incident_bonds(mol, root, ranking)?,
        next_edge: 0,
    }];
    while let Some(frame) = stack.last_mut() {
        if frame.next_edge >= frame.incident.len() {
            stack.pop();
            continue;
        }
        let (bond_id, order, neighbor) = frame.incident[frame.next_edge];
        frame.next_edge += 1;
        if Some(bond_id) == frame.parent_bond {
            continue;
        }
        if visited.contains(&neighbor) {
            if !tree_bonds.contains(&bond_id) {
                let bond = mol
                    .bond(bond_id)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                ring_bonds
                    .entry(bond_id)
                    .or_insert((bond.a(), bond.b(), order));
            }
            continue;
        }
        tree_bonds.insert(bond_id);
        visited.insert(neighbor);
        stack.push(Frame {
            parent_bond: Some(bond_id),
            incident: canonical_smiles_incident_bonds(mol, neighbor, ranking)?,
            next_edge: 0,
        });
    }

    let mut ring_bonds = ring_bonds
        .into_iter()
        .map(|(bond_id, (a, b, order))| {
            let (first, second) = ordered_atom_pair(a, b);
            (bond_id, first, second, order)
        })
        .collect::<Vec<_>>();
    ring_bonds.sort_by_key(|(bond_id, first, second, order)| {
        (
            canonical_rank(ranking, *first),
            canonical_rank(ranking, *second),
            bond_order_code(*order),
            *first,
            *second,
            *bond_id,
        )
    });
    if ring_bonds.len() > 99 {
        return Err(MolWriteError::new(
            "SMILES writer supports at most 99 simultaneous ring closures",
        ));
    }

    let mut closures = BTreeMap::<AtomId, Vec<SmilesRingClosure>>::new();
    for (index, (_, first, second, order)) in ring_bonds.into_iter().enumerate() {
        let number = index + 1;
        closures.entry(first).or_default().push(SmilesRingClosure {
            number,
            order,
            other: second,
        });
        closures.entry(second).or_default().push(SmilesRingClosure {
            number,
            order,
            other: first,
        });
    }
    for (atom, closures) in &mut closures {
        closures.sort_by_key(|closure| {
            (
                canonical_rank(ranking, closure.other),
                bond_order_code(closure.order),
                closure.other,
                *atom,
            )
        });
    }

    Ok(SmilesWritePlan {
        roots: vec![root],
        tree_bonds,
        closures,
        subtree_sizes: BTreeMap::new(),
    })
}

fn write_canonical_smiles_component_with_plan(
    mol: &Molecule,
    root: AtomId,
    plan: &SmilesWritePlan,
    ranking: &CanonicalAtomRanking,
) -> std::result::Result<String, MolWriteError> {
    enum Action {
        Node {
            atom: AtomId,
            parent: Option<AtomId>,
        },
        Bond {
            order: BondOrder,
            left: AtomId,
            right: AtomId,
        },
        OpenBranch,
        CloseBranch,
    }

    let mut out = String::new();
    let mut actions = vec![Action::Node {
        atom: root,
        parent: None,
    }];
    while let Some(action) = actions.pop() {
        match action {
            Action::OpenBranch => out.push('('),
            Action::CloseBranch => out.push(')'),
            Action::Bond { order, left, right } => {
                out.push_str(smiles_bond_between(mol, order, left, right)?);
            }
            Action::Node { atom, parent } => {
                let atom_record = mol
                    .atom(atom)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                out.push_str(&smiles_atom(atom_record));
                if let Some(closures) = plan.closures.get(&atom) {
                    for closure in closures {
                        out.push_str(smiles_bond_between(
                            mol,
                            closure.order,
                            atom,
                            closure.other,
                        )?);
                        out.push_str(&smiles_ring_number(closure.number));
                    }
                }

                let mut children = canonical_smiles_incident_bonds(mol, atom, ranking)?
                    .into_iter()
                    .filter(|(bond_id, _, neighbor)| {
                        plan.tree_bonds.contains(bond_id) && Some(*neighbor) != parent
                    })
                    .collect::<Vec<_>>();
                children.sort_by_key(|(bond_id, order, child)| {
                    (
                        canonical_rank(ranking, *child),
                        smiles_atom(
                            mol.atom(*child)
                                .expect("canonical tree child should remain live"),
                        ),
                        bond_order_code(*order),
                        *child,
                        *bond_id,
                    )
                });
                for (_, order, child) in children.into_iter().rev() {
                    actions.push(Action::CloseBranch);
                    actions.push(Action::Node {
                        atom: child,
                        parent: Some(atom),
                    });
                    actions.push(Action::Bond {
                        order,
                        left: atom,
                        right: child,
                    });
                    actions.push(Action::OpenBranch);
                }
            }
        }
    }
    Ok(out)
}

fn compute_smiles_subtree_sizes(
    mol: &Molecule,
    atom_id: AtomId,
    parent: Option<AtomId>,
    tree_bonds: &BTreeSet<BondId>,
    subtree_sizes: &mut BTreeMap<AtomId, usize>,
) -> std::result::Result<usize, MolWriteError> {
    let mut stack = vec![(atom_id, parent, false)];
    while let Some((current, parent, expanded)) = stack.pop() {
        if expanded {
            let mut size = 1usize;
            for (bond_id, _, neighbor) in smiles_incident_bonds(mol, current)? {
                if tree_bonds.contains(&bond_id) && Some(neighbor) != parent {
                    size = size
                        .saturating_add(subtree_sizes.get(&neighbor).copied().unwrap_or_default());
                }
            }
            subtree_sizes.insert(current, size);
            continue;
        }
        stack.push((current, parent, true));
        let mut children = smiles_incident_bonds(mol, current)?
            .into_iter()
            .filter(|(bond_id, _, neighbor)| {
                tree_bonds.contains(bond_id) && Some(*neighbor) != parent
            })
            .map(|(_, _, neighbor)| neighbor)
            .collect::<Vec<_>>();
        children.sort();
        for child in children.into_iter().rev() {
            stack.push((child, Some(current), false));
        }
    }
    Ok(subtree_sizes.get(&atom_id).copied().unwrap_or_default())
}

fn write_smiles_component(
    mol: &Molecule,
    atom_id: AtomId,
    parent: Option<AtomId>,
    plan: &SmilesWritePlan,
) -> std::result::Result<String, MolWriteError> {
    enum Action {
        Node {
            atom: AtomId,
            parent: Option<AtomId>,
        },
        Bond {
            order: BondOrder,
            left: AtomId,
            right: AtomId,
        },
        OpenBranch,
        CloseBranch,
    }

    let mut out = String::new();
    let mut actions = vec![Action::Node {
        atom: atom_id,
        parent,
    }];
    while let Some(action) = actions.pop() {
        match action {
            Action::OpenBranch => out.push('('),
            Action::CloseBranch => out.push(')'),
            Action::Bond { order, left, right } => {
                out.push_str(smiles_bond_between(mol, order, left, right)?);
            }
            Action::Node { atom, parent } => {
                let atom_record = mol
                    .atom(atom)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                out.push_str(&smiles_atom(atom_record));
                if let Some(closures) = plan.closures.get(&atom) {
                    for closure in closures {
                        out.push_str(smiles_bond_between(
                            mol,
                            closure.order,
                            atom,
                            closure.other,
                        )?);
                        out.push_str(&smiles_ring_number(closure.number));
                    }
                }

                let mut children = smiles_incident_bonds(mol, atom)?
                    .into_iter()
                    .filter(|(bond_id, _, neighbor)| {
                        plan.tree_bonds.contains(bond_id) && Some(*neighbor) != parent
                    })
                    .collect::<Vec<_>>();
                children.sort_by_key(|(bond_id, _, child)| (*child, *bond_id));
                let main_child_index = children
                    .iter()
                    .enumerate()
                    .max_by_key(|(_, child_entry)| {
                        let child = child_entry.2;
                        (plan.subtree_sizes.get(&child).copied().unwrap_or(0), child)
                    })
                    .map(|(index, _)| index);

                if let Some(index) = main_child_index {
                    let (_, order, child) = children[index];
                    actions.push(Action::Node {
                        atom: child,
                        parent: Some(atom),
                    });
                    actions.push(Action::Bond {
                        order,
                        left: atom,
                        right: child,
                    });
                }
                for (index, (_, order, child)) in children.into_iter().enumerate().rev() {
                    if Some(index) == main_child_index {
                        continue;
                    }
                    actions.push(Action::CloseBranch);
                    actions.push(Action::Node {
                        atom: child,
                        parent: Some(atom),
                    });
                    actions.push(Action::Bond {
                        order,
                        left: atom,
                        right: child,
                    });
                    actions.push(Action::OpenBranch);
                }
            }
        }
    }
    Ok(out)
}

fn smiles_incident_bonds(
    mol: &Molecule,
    atom_id: AtomId,
) -> std::result::Result<Vec<(BondId, BondOrder, AtomId)>, MolWriteError> {
    let mut incident = Vec::new();
    for (bond_id, bond) in mol
        .incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        incident.push((bond_id, bond.order, bond.other_atom(atom_id)));
    }
    incident.sort_by_key(|(bond_id, _, atom)| (*atom, *bond_id));
    Ok(incident)
}

fn canonical_smiles_incident_bonds(
    mol: &Molecule,
    atom_id: AtomId,
    ranking: &CanonicalAtomRanking,
) -> std::result::Result<Vec<(BondId, BondOrder, AtomId)>, MolWriteError> {
    let mut incident = smiles_incident_bonds(mol, atom_id)?;
    incident.sort_by_key(|(bond_id, order, atom)| {
        (
            canonical_rank(ranking, *atom),
            smiles_atom(
                mol.atom(*atom)
                    .expect("incident atom from live bond should remain live"),
            ),
            bond_order_code(*order),
            *atom,
            *bond_id,
        )
    });
    Ok(incident)
}

fn canonical_rank(ranking: &CanonicalAtomRanking, atom: AtomId) -> u32 {
    ranking
        .rank_of(atom)
        .expect("canonical ranking should cover every live atom")
}

fn bond_order_code(order: BondOrder) -> u8 {
    match order {
        BondOrder::Zero => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 5,
        BondOrder::Dative => 6,
    }
}

fn smiles_ring_number(number: usize) -> String {
    if number < 10 {
        number.to_string()
    } else {
        format!("%{number}")
    }
}

fn smiles_bond(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Aromatic => ":",
        BondOrder::Zero | BondOrder::Dative | BondOrder::Quadruple => "-",
    }
}

fn smiles_bond_between(
    mol: &Molecule,
    order: BondOrder,
    left: AtomId,
    right: AtomId,
) -> std::result::Result<&'static str, MolWriteError> {
    if matches!(order, BondOrder::Single | BondOrder::Aromatic) {
        let left = mol
            .atom(left)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let right = mol
            .atom(right)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        if left.aromatic && right.aromatic {
            return Ok(if order == BondOrder::Single { "-" } else { "" });
        }
    }
    Ok(smiles_bond(order))
}

fn smiles_atom(atom: &Atom) -> String {
    let explicit_hydrogens = smiles_atom_explicit_hydrogens(atom);
    let organic = atom.isotope.is_none()
        && atom.formal_charge == 0
        && explicit_hydrogens == 0
        && !atom.no_implicit_hydrogens
        && atom.atom_map.is_none()
        && matches!(
            atom.element.symbol(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        );
    if organic {
        if atom.aromatic {
            atom.element.symbol().to_ascii_lowercase()
        } else {
            atom.element.symbol().to_owned()
        }
    } else {
        let mut out = String::from("[");
        if let Some(isotope) = atom.isotope {
            out.push_str(&isotope.to_string());
        }
        if atom.aromatic {
            out.push_str(&atom.element.symbol().to_ascii_lowercase());
        } else {
            out.push_str(atom.element.symbol());
        }
        if explicit_hydrogens > 0 {
            out.push('H');
            if explicit_hydrogens > 1 {
                out.push_str(&explicit_hydrogens.to_string());
            }
        }
        if atom.formal_charge > 0 {
            out.push('+');
            if atom.formal_charge > 1 {
                out.push_str(&atom.formal_charge.to_string());
            }
        } else if atom.formal_charge < 0 {
            out.push('-');
            if atom.formal_charge < -1 {
                out.push_str(&(-atom.formal_charge).to_string());
            }
        }
        if let Some(map) = atom.atom_map {
            out.push(':');
            out.push_str(&map.to_string());
        }
        out.push(']');
        out
    }
}

fn smiles_atom_explicit_hydrogens(atom: &Atom) -> u8 {
    if atom.element.symbol() == "N"
        && atom.aromatic
        && atom.explicit_hydrogens == 0
        && atom.implicit_hydrogens == Some(1)
    {
        1
    } else {
        atom.explicit_hydrogens
    }
}
