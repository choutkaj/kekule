use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::algorithms::{canonical_atom_ranking, ordered_atom_pair, CanonicalAtomRanking};
use crate::core::*;
use crate::io::MolWriteError;
use crate::small::SmallMolecule;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct SmilesParseOptions;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct SmilesWriteOptions;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct IsomericSmilesWriteOptions;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
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
    let mut pending_bond = None::<(BondOrder, Option<StereoBondMarkKind>, usize)>;
    let mut rings = BTreeMap::<
        usize,
        (
            AtomId,
            Option<(BondOrder, Option<StereoBondMarkKind>)>,
            usize,
        ),
    >::new();
    let mut pending_tetrahedral = Vec::<PendingTetrahedral>::new();
    let mut tetrahedral_carriers = BTreeMap::<AtomId, Vec<PendingStereoCarrier>>::new();
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
            '-' | '=' | '#' | ':' | '/' | '\\' => {
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
                    '/' | '\\' => BondOrder::Single,
                    _ => unreachable!(),
                };
                let stereo = match ch {
                    '/' => Some(StereoBondMarkKind::DirectionalUp),
                    '\\' => Some(StereoBondMarkKind::DirectionalDown),
                    _ => None,
                };
                pending_bond = Some((order, stereo, offset));
                previous = SmilesTokenKind::Bond;
                cursor += 1;
            }
            '0'..='9' | '%' => {
                let atom = current
                    .ok_or_else(|| SmilesParseError::new(offset, "ring closure without atom"))?;
                let (label, next_cursor) = parse_smiles_ring_label(&chars, cursor)?;
                let close_bond = pending_bond
                    .take()
                    .map(|(order, stereo, _)| (order, stereo));
                if let Some((other, open_bond, open_component)) = rings.remove(&label) {
                    if open_component != component {
                        return Err(SmilesParseError::new(
                            offset,
                            "ring closure crosses a component separator",
                        ));
                    }
                    if open_bond.is_some() && close_bond.is_some() && open_bond != close_bond {
                        return Err(SmilesParseError::new(
                            offset,
                            "conflicting ring bond symbols",
                        ));
                    }
                    let (order, stereo) = match close_bond.or(open_bond) {
                        Some((order, stereo)) => (order, stereo),
                        None => (default_smiles_bond_order(&mol, other, atom, offset)?, None),
                    };
                    add_smiles_bond(&mut mol, other, atom, order, stereo, offset)?;
                    resolve_tetrahedral_ring_carrier(
                        &mut tetrahedral_carriers,
                        other,
                        label,
                        component,
                        atom,
                    );
                    push_tetrahedral_carrier(
                        &mut tetrahedral_carriers,
                        atom,
                        StereoCarrier::Atom(other),
                    );
                } else {
                    rings.insert(label, (atom, close_bond, component));
                    push_tetrahedral_ring_carrier(
                        &mut tetrahedral_carriers,
                        atom,
                        label,
                        component,
                    );
                }
                previous = SmilesTokenKind::Ring;
                cursor = next_cursor;
            }
            '[' => {
                let (atom, chirality, next_cursor) = parse_bracket_atom(&chars, cursor)?;
                let explicit_hydrogens = atom.explicit_hydrogens;
                let atom_id = mol.add_atom(atom);
                if let Some(orientation) = chirality {
                    pending_tetrahedral.push(PendingTetrahedral {
                        center: atom_id,
                        orientation,
                    });
                    tetrahedral_carriers.insert(
                        atom_id,
                        initial_tetrahedral_carriers(current, explicit_hydrogens),
                    );
                }
                if let Some(previous) = current {
                    let (order, stereo) = match pending_bond
                        .take()
                        .map(|(order, stereo, _)| (order, stereo))
                    {
                        Some((order, stereo)) => (order, stereo),
                        None => (
                            default_smiles_bond_order(&mol, previous, atom_id, offset)?,
                            None,
                        ),
                    };
                    add_smiles_bond(&mut mol, previous, atom_id, order, stereo, offset)?;
                    push_tetrahedral_carrier(
                        &mut tetrahedral_carriers,
                        previous,
                        StereoCarrier::Atom(atom_id),
                    );
                } else if pending_bond.is_some() {
                    return Err(SmilesParseError::new(offset, "bond without left endpoint"));
                }
                current = Some(atom_id);
                previous = SmilesTokenKind::Atom;
                cursor = next_cursor;
            }
            '@' | '*' => {
                return Err(SmilesParseError::new(
                    offset,
                    "unsupported stereochemistry or query syntax",
                ));
            }
            _ => {
                let (atom, next_cursor) = parse_organic_atom(&chars, cursor)?;
                let atom_id = mol.add_atom(atom);
                if let Some(previous) = current {
                    let (order, stereo) = match pending_bond
                        .take()
                        .map(|(order, stereo, _)| (order, stereo))
                    {
                        Some((order, stereo)) => (order, stereo),
                        None => (
                            default_smiles_bond_order(&mol, previous, atom_id, offset)?,
                            None,
                        ),
                    };
                    add_smiles_bond(&mut mol, previous, atom_id, order, stereo, offset)?;
                    push_tetrahedral_carrier(
                        &mut tetrahedral_carriers,
                        previous,
                        StereoCarrier::Atom(atom_id),
                    );
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
    if let Some((_, _, offset)) = pending_bond {
        return Err(SmilesParseError::new(offset, "bond without right endpoint"));
    }
    if matches!(previous, SmilesTokenKind::Dot | SmilesTokenKind::BranchOpen) {
        return Err(SmilesParseError::new(input.len(), "incomplete SMILES"));
    }
    add_smiles_tetrahedral_elements(
        &mut mol,
        pending_tetrahedral,
        tetrahedral_carriers,
        input.len(),
    )?;
    Ok(SmallMolecule::from_graph(mol))
}

fn add_smiles_bond(
    mol: &mut Molecule,
    left: AtomId,
    right: AtomId,
    order: BondOrder,
    stereo: Option<StereoBondMarkKind>,
    offset: usize,
) -> std::result::Result<(), SmilesParseError> {
    let bond_id = mol
        .add_bond(left, right, order)
        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
    if let Some(kind) = stereo {
        mol.set_stereo_bond_mark(StereoBondMark {
            bond: bond_id,
            kind,
            source: StereoSource::Smiles,
        })
        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
    }
    Ok(())
}

fn add_smiles_tetrahedral_elements(
    mol: &mut Molecule,
    centers: Vec<PendingTetrahedral>,
    mut carriers_by_center: BTreeMap<AtomId, Vec<PendingStereoCarrier>>,
    offset: usize,
) -> std::result::Result<(), SmilesParseError> {
    for pending in centers {
        let carriers = resolve_smiles_tetrahedral_carriers(
            mol,
            pending.center,
            carriers_by_center
                .remove(&pending.center)
                .unwrap_or_default(),
            offset,
        )?;
        mol.add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center: pending.center,
                carriers,
                orientation: pending.orientation,
            }),
            StereoSource::Smiles,
        ))
        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
    }
    Ok(())
}

fn initial_tetrahedral_carriers(
    previous: Option<AtomId>,
    explicit_hydrogens: u8,
) -> Vec<PendingStereoCarrier> {
    let mut carriers = Vec::new();
    if let Some(previous) = previous {
        carriers.push(PendingStereoCarrier::Resolved(StereoCarrier::Atom(
            previous,
        )));
    }
    for _ in 0..explicit_hydrogens {
        carriers.push(PendingStereoCarrier::Resolved(
            StereoCarrier::ImplicitHydrogen,
        ));
    }
    carriers
}

fn push_tetrahedral_carrier(
    carriers_by_center: &mut BTreeMap<AtomId, Vec<PendingStereoCarrier>>,
    center: AtomId,
    carrier: StereoCarrier,
) {
    if let Some(carriers) = carriers_by_center.get_mut(&center) {
        carriers.push(PendingStereoCarrier::Resolved(carrier));
    }
}

fn push_tetrahedral_ring_carrier(
    carriers_by_center: &mut BTreeMap<AtomId, Vec<PendingStereoCarrier>>,
    center: AtomId,
    label: usize,
    component: usize,
) {
    if let Some(carriers) = carriers_by_center.get_mut(&center) {
        carriers.push(PendingStereoCarrier::Ring { label, component });
    }
}

fn resolve_tetrahedral_ring_carrier(
    carriers_by_center: &mut BTreeMap<AtomId, Vec<PendingStereoCarrier>>,
    center: AtomId,
    label: usize,
    component: usize,
    carrier: AtomId,
) {
    let Some(carriers) = carriers_by_center.get_mut(&center) else {
        return;
    };
    if let Some(pending) = carriers.iter_mut().find(|pending| {
        matches!(
            pending,
            PendingStereoCarrier::Ring {
                label: pending_label,
                component: pending_component,
            } if *pending_label == label && *pending_component == component
        )
    }) {
        *pending = PendingStereoCarrier::Resolved(StereoCarrier::Atom(carrier));
    }
}

fn resolve_smiles_tetrahedral_carriers(
    mol: &Molecule,
    center: AtomId,
    carriers: Vec<PendingStereoCarrier>,
    offset: usize,
) -> std::result::Result<Vec<StereoCarrier>, SmilesParseError> {
    let mut carriers = carriers
        .into_iter()
        .map(|carrier| match carrier {
            PendingStereoCarrier::Resolved(carrier) => Ok(carrier),
            PendingStereoCarrier::Ring { .. } => Err(SmilesParseError::new(
                offset,
                "unresolved tetrahedral ring carrier",
            )),
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    if carriers.len() == 3 && smiles_tetrahedral_center_can_have_lone_pair(mol, center) {
        carriers.push(StereoCarrier::ImplicitLonePair);
    }
    Ok(carriers)
}

fn smiles_tetrahedral_center_can_have_lone_pair(mol: &Molecule, center: AtomId) -> bool {
    mol.atom(center)
        .map(|atom| {
            matches!(
                atom.element.symbol(),
                "N" | "P" | "As" | "Sb" | "O" | "S" | "Se" | "Te"
            ) && atom.explicit_hydrogens == 0
                && atom.implicit_hydrogens.unwrap_or(0) == 0
        })
        .unwrap_or(false)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PendingTetrahedral {
    center: AtomId,
    orientation: TetrahedralOrientation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PendingStereoCarrier {
    Resolved(StereoCarrier),
    Ring { label: usize, component: usize },
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
) -> std::result::Result<(Atom, Option<TetrahedralOrientation>, usize), SmilesParseError> {
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
        let (symbol, symbol_len) =
            parse_aromatic_bracket_element(bytes, index).ok_or_else(|| {
                SmilesParseError::new(start + 1 + index, "unsupported aromatic bracket element")
            })?;
        index += symbol_len;
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
    let mut saw_chirality = false;
    let mut chirality = None;
    let mut saw_hydrogen = false;
    let mut saw_charge = false;
    let mut saw_map = false;
    while index < text.len() {
        match bytes[index] {
            b'@' if !saw_chirality && !saw_hydrogen && !saw_charge && !saw_map => {
                saw_chirality = true;
                index += 1;
                chirality = if bytes.get(index) == Some(&b'@') {
                    index += 1;
                    Some(TetrahedralOrientation::CounterClockwise)
                } else {
                    Some(TetrahedralOrientation::Clockwise)
                };
            }
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
            b'/' | b'\\' | b'*' => {
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
    Ok((atom, chirality, end + 1))
}

fn parse_aromatic_bracket_element(bytes: &[u8], index: usize) -> Option<(&'static str, usize)> {
    match bytes.get(index)? {
        b'b' => Some(("B", 1)),
        b'c' => Some(("C", 1)),
        b'n' => Some(("N", 1)),
        b'o' => Some(("O", 1)),
        b'p' => Some(("P", 1)),
        b's' if bytes.get(index + 1) == Some(&b'e') => Some(("Se", 2)),
        b's' => Some(("S", 1)),
        b't' if bytes.get(index + 1) == Some(&b'e') => Some(("Te", 2)),
        _ => None,
    }
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
    let mol = molecule.graph();
    let plan = plan_smiles_write(mol, StereoWriteMode::Reject)?;
    let mut parts = Vec::new();
    for start in &plan.roots {
        parts.push(write_smiles_component(mol, *start, None, &plan, None)?);
    }
    Ok(parts.join("."))
}

pub fn write_isomeric_smiles(
    molecule: &SmallMolecule,
    _options: IsomericSmilesWriteOptions,
) -> std::result::Result<String, MolWriteError> {
    let mol = molecule.graph();
    let plan = plan_smiles_write(mol, StereoWriteMode::Encode)?;
    let stereo = SmilesStereoWriteContext::new(mol)?;
    let mut parts = Vec::new();
    for start in &plan.roots {
        parts.push(write_smiles_component(
            mol,
            *start,
            None,
            &plan,
            Some(&stereo),
        )?);
    }
    Ok(parts.join("."))
}

pub fn write_canonical_smiles(
    molecule: &SmallMolecule,
    _options: CanonicalSmilesWriteOptions,
) -> std::result::Result<String, MolWriteError> {
    let mol = molecule.graph();
    validate_smiles_writeable(mol, StereoWriteMode::Ignore)?;
    let ranking = canonical_atom_ranking(mol);
    let mut components = Vec::new();
    for component in smiles_connected_components(mol)? {
        let atom_style = canonical_component_atom_style(mol, &component)?;
        let mut candidates = Vec::new();
        for preference in [
            CanonicalBondTraversal::HighOrderFirst,
            CanonicalBondTraversal::LowOrderFirst,
        ] {
            candidates.extend(
                component
                    .iter()
                    .map(|root| {
                        write_canonical_smiles_component(
                            mol, *root, &ranking, preference, atom_style,
                        )
                    })
                    .collect::<std::result::Result<Vec<_>, _>>()?,
            );
        }
        candidates.sort_by_key(|candidate| canonical_smiles_candidate_key(candidate));
        candidates.dedup();
        if let Some(candidate) = candidates.into_iter().next() {
            components.push(candidate);
        }
    }
    components.sort();
    Ok(components.join("."))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CanonicalBondTraversal {
    HighOrderFirst,
    LowOrderFirst,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CanonicalAtomStyle {
    Aromatic,
    StoredKekule,
}

impl CanonicalBondTraversal {
    fn order_key(self, order: BondOrder) -> u8 {
        match self {
            Self::HighOrderFirst => reverse_bond_order_code(order),
            Self::LowOrderFirst => bond_order_code(order),
        }
    }
}

fn canonical_component_atom_style(
    mol: &Molecule,
    atom_ids: &[AtomId],
) -> std::result::Result<CanonicalAtomStyle, MolWriteError> {
    if canonical_component_has_aromatic_shorthand_sensitive_atom(mol, atom_ids)?
        && canonical_component_has_stored_kekule_orders(mol, atom_ids)?
    {
        Ok(CanonicalAtomStyle::StoredKekule)
    } else {
        Ok(CanonicalAtomStyle::Aromatic)
    }
}

fn canonical_component_has_aromatic_shorthand_sensitive_atom(
    mol: &Molecule,
    atom_ids: &[AtomId],
) -> std::result::Result<bool, MolWriteError> {
    let atom_set = atom_ids.iter().copied().collect::<BTreeSet<_>>();
    let component_has_aromatic_atom = atom_ids.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| atom.aromatic)
            .unwrap_or(false)
    });
    if !component_has_aromatic_atom {
        return Ok(false);
    }
    for atom_id in atom_ids {
        let atom = mol
            .atom(*atom_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        if atom.aromatic && atom.formal_charge != 0 && matches!(atom.element.symbol(), "B" | "C") {
            return Ok(true);
        }
        if atom.aromatic && atom_has_exocyclic_hetero_multiple_bond(mol, *atom_id, &atom_set)? {
            return Ok(true);
        }
        if atom.aromatic || atom.formal_charge != 0 {
            continue;
        }
        let mut aromatic_neighbors = 0usize;
        let mut pi_framework_neighbors = 0usize;
        let mut multiple_bond_to_non_aromatic_neighbor = false;
        for (_, bond) in mol
            .incident_bonds(*atom_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?
        {
            let neighbor_id = bond.other_atom(*atom_id);
            let neighbor = mol
                .atom(neighbor_id)
                .map_err(|error| MolWriteError::new(error.to_string()))?;
            if atom_set.contains(&neighbor_id) && neighbor.aromatic {
                aromatic_neighbors += 1;
            }
            if atom_set.contains(&neighbor_id)
                && matches!(neighbor.element.symbol(), "B" | "C" | "N" | "P" | "S")
            {
                pi_framework_neighbors += 1;
            }
            if matches!(bond.order, BondOrder::Double | BondOrder::Triple) && !neighbor.aromatic {
                multiple_bond_to_non_aromatic_neighbor = true;
            }
        }
        if (aromatic_neighbors > 0 || pi_framework_neighbors >= 3)
            && pi_framework_neighbors >= 2
            && multiple_bond_to_non_aromatic_neighbor
        {
            return Ok(true);
        }
    }
    Ok(false)
}

fn atom_has_exocyclic_hetero_multiple_bond(
    mol: &Molecule,
    atom_id: AtomId,
    atom_set: &BTreeSet<AtomId>,
) -> std::result::Result<bool, MolWriteError> {
    for (_, bond) in mol
        .incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        if !matches!(bond.order, BondOrder::Double | BondOrder::Triple) {
            continue;
        }
        let neighbor_id = bond.other_atom(atom_id);
        let neighbor = mol
            .atom(neighbor_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        if !atom_set.contains(&neighbor_id)
            || !neighbor.aromatic && !matches!(neighbor.element.symbol(), "B" | "C")
        {
            return Ok(true);
        }
    }
    Ok(false)
}

fn canonical_component_has_stored_kekule_orders(
    mol: &Molecule,
    atom_ids: &[AtomId],
) -> std::result::Result<bool, MolWriteError> {
    let atom_set = atom_ids.iter().copied().collect::<BTreeSet<_>>();
    for atom_id in atom_ids {
        for (_, bond) in mol
            .incident_bonds(*atom_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?
        {
            let neighbor_id = bond.other_atom(*atom_id);
            if *atom_id < neighbor_id
                && atom_set.contains(&neighbor_id)
                && bond.aromatic
                && matches!(bond.order, BondOrder::Aromatic)
            {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn canonical_smiles_candidate_key(candidate: &str) -> (usize, usize, usize, String) {
    (
        candidate.matches('(').count(),
        explicit_ring_bond_marker_count(candidate),
        leading_ring_label_count(candidate),
        candidate.to_owned(),
    )
}

fn leading_ring_label_count(candidate: &str) -> usize {
    let bytes = candidate.as_bytes();
    let mut index = smiles_atom_token_end(candidate);
    let mut count = 0usize;
    while let Some(byte) = bytes.get(index) {
        if byte.is_ascii_digit() {
            count += 1;
            index += 1;
        } else if *byte == b'%' && bytes.get(index + 1).is_some_and(u8::is_ascii_digit) {
            count += 1;
            index += 3;
        } else {
            break;
        }
    }
    count
}

fn smiles_atom_token_end(candidate: &str) -> usize {
    let bytes = candidate.as_bytes();
    if bytes.first() == Some(&b'[') {
        return bytes
            .iter()
            .position(|byte| *byte == b']')
            .map(|index| index + 1)
            .unwrap_or(candidate.len());
    }
    if matches!(bytes.first(), Some(b'B' | b'C')) && matches!(bytes.get(1), Some(b'l' | b'r')) {
        2
    } else {
        bytes.first().map(|_| 1).unwrap_or(0)
    }
}

fn explicit_ring_bond_marker_count(candidate: &str) -> usize {
    let bytes = candidate.as_bytes();
    bytes
        .windows(2)
        .filter(|pair| matches!(pair[0], b'-' | b'=' | b'#' | b':') && pair[1].is_ascii_digit())
        .count()
        + bytes
            .windows(2)
            .filter(|pair| matches!(pair[0], b'-' | b'=' | b'#' | b':') && pair[1] == b'%')
            .count()
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
    bond: BondId,
    number: usize,
    order: BondOrder,
    other: AtomId,
}

fn plan_smiles_write(
    mol: &Molecule,
    stereo: StereoWriteMode,
) -> std::result::Result<SmilesWritePlan, MolWriteError> {
    validate_smiles_writeable(mol, stereo)?;
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
    for (index, (bond_id, first, second, order)) in ring_bonds.into_iter().enumerate() {
        let number = index + 1;
        closures.entry(first).or_default().push(SmilesRingClosure {
            bond: bond_id,
            number,
            order,
            other: second,
        });
        closures.entry(second).or_default().push(SmilesRingClosure {
            bond: bond_id,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StereoWriteMode {
    Reject,
    Ignore,
    Encode,
}

fn validate_smiles_writeable(
    mol: &Molecule,
    stereo: StereoWriteMode,
) -> std::result::Result<(), MolWriteError> {
    match stereo {
        StereoWriteMode::Reject if mol.stereo_elements().next().is_some() => {
            return Err(MolWriteError::new(
                "SMILES writer cannot encode atom stereochemistry",
            ));
        }
        StereoWriteMode::Reject if mol.stereo_bond_marks().next().is_some() => {
            return Err(MolWriteError::new(
                "SMILES writer cannot encode bond stereochemistry",
            ));
        }
        StereoWriteMode::Encode => validate_isomeric_smiles_stereo(mol)?,
        StereoWriteMode::Reject | StereoWriteMode::Ignore => {}
    }
    for (_, atom) in mol.atoms() {
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
    }
    Ok(())
}

fn validate_isomeric_smiles_stereo(mol: &Molecule) -> std::result::Result<(), MolWriteError> {
    if mol.stereo_groups().next().is_some() {
        return Err(MolWriteError::new(
            "isomeric SMILES writer cannot encode enhanced stereo groups",
        ));
    }
    for (_, element) in mol.stereo_elements() {
        if element.specifiedness != StereoSpecifiedness::Specified {
            return Err(MolWriteError::new(
                "isomeric SMILES writer cannot encode unspecified or unknown stereo",
            ));
        }
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                if stereo.carriers.len() != 4 {
                    return Err(MolWriteError::new(
                        "isomeric SMILES writer cannot encode invalid tetrahedral stereo",
                    ));
                }
                let hydrogen_count = stereo
                    .carriers
                    .iter()
                    .filter(|carrier| matches!(carrier, StereoCarrier::ImplicitHydrogen))
                    .count();
                if hydrogen_count > 1 {
                    return Err(MolWriteError::new(
                        "isomeric SMILES writer cannot encode tetrahedral stereo with repeated implicit hydrogens",
                    ));
                }
            }
            StereoElementKind::DoubleBond(stereo) => {
                validate_isomeric_double_bond_endpoint(
                    mol,
                    stereo.left,
                    stereo.right,
                    stereo.bond,
                    stereo.left_carrier,
                )?;
                validate_isomeric_double_bond_endpoint(
                    mol,
                    stereo.right,
                    stereo.left,
                    stereo.bond,
                    stereo.right_carrier,
                )?;
            }
            StereoElementKind::Axis(_) => {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode axial stereochemistry yet",
                ));
            }
        }
    }
    Ok(())
}

fn validate_isomeric_double_bond_endpoint(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
    carrier: StereoCarrier,
) -> std::result::Result<(), MolWriteError> {
    match carrier {
        StereoCarrier::Atom(atom) => {
            let bond = mol
                .bond_between(endpoint, atom)
                .map_err(|error| MolWriteError::new(error.to_string()))?
                .ok_or_else(|| MolWriteError::new("double-bond stereo carrier is not bonded"))?;
            let order = mol
                .bond(bond)
                .map_err(|error| MolWriteError::new(error.to_string()))?
                .order;
            if order != BondOrder::Single || atom == other_endpoint {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode invalid double-bond stereo carrier",
                ));
            }
        }
        StereoCarrier::ImplicitHydrogen => {
            let atom = mol
                .atom(endpoint)
                .map_err(|error| MolWriteError::new(error.to_string()))?;
            let hydrogens = atom
                .explicit_hydrogens
                .saturating_add(atom.implicit_hydrogens.unwrap_or(0));
            if hydrogens == 0 {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode unavailable implicit double-bond hydrogen carrier",
                ));
            }
            if implicit_double_bond_printable_carrier_bond(
                mol,
                endpoint,
                other_endpoint,
                focus_bond,
            )?
            .is_none()
            {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode implicit double-bond carrier without a unique explicit substituent bond",
                ));
            }
        }
        StereoCarrier::ImplicitLonePair => {
            return Err(MolWriteError::new(
                "isomeric SMILES writer cannot encode lone-pair double-bond carrier",
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
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<String, MolWriteError> {
    let plan = plan_canonical_smiles_component(mol, root, ranking, preference, atom_style)?;
    write_canonical_smiles_component_with_plan(mol, root, &plan, ranking, preference, atom_style)
}

fn plan_canonical_smiles_component(
    mol: &Molecule,
    root: AtomId,
    ranking: &CanonicalAtomRanking,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
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
        incident: canonical_smiles_incident_bonds(mol, root, ranking, preference, atom_style)?,
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
            incident: canonical_smiles_incident_bonds(
                mol, neighbor, ranking, preference, atom_style,
            )?,
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
    for (index, (bond_id, first, second, order)) in ring_bonds.into_iter().enumerate() {
        let number = index + 1;
        closures.entry(first).or_default().push(SmilesRingClosure {
            bond: bond_id,
            number,
            order,
            other: second,
        });
        closures.entry(second).or_default().push(SmilesRingClosure {
            bond: bond_id,
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
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
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
                out.push_str(&canonical_smiles_atom(mol, atom, atom_record, atom_style)?);
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

                let mut children =
                    canonical_smiles_incident_bonds(mol, atom, ranking, preference, atom_style)?
                        .into_iter()
                        .filter(|(bond_id, _, neighbor)| {
                            plan.tree_bonds.contains(bond_id) && Some(*neighbor) != parent
                        })
                        .collect::<Vec<_>>();
                children.sort_by_key(|(bond_id, order, child)| {
                    (
                        !canonical_smiles_aromatic_continuation(mol, atom, *child, *order),
                        canonical_rank(ranking, *child),
                        canonical_smiles_atom_for_sort(mol, *child, atom_style),
                        preference.order_key(*order),
                        *child,
                        *bond_id,
                    )
                });
                let main_child = children.first().copied();
                if let Some((_, order, child)) = main_child {
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
                    if index == 0 {
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

fn canonical_smiles_aromatic_continuation(
    mol: &Molecule,
    left: AtomId,
    right: AtomId,
    order: BondOrder,
) -> bool {
    matches!(order, BondOrder::Aromatic)
        && mol.atom(left).is_ok_and(|atom| atom.aromatic)
        && mol.atom(right).is_ok_and(|atom| atom.aromatic)
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

#[derive(Debug, Clone)]
struct SmilesStereoWriteContext {
    tetrahedral: BTreeMap<AtomId, TetrahedralSmilesState>,
    directional: BTreeMap<BondId, Vec<DirectionalSmilesConstraint>>,
}

#[derive(Debug, Clone)]
struct TetrahedralSmilesState {
    carriers: Vec<StereoCarrier>,
    orientation: TetrahedralOrientation,
}

#[derive(Debug, Clone, Copy)]
struct ChiralAtomWriteState {
    orientation: TetrahedralOrientation,
    force_hydrogen: bool,
}

#[derive(Debug, Clone, Copy)]
struct DirectionalSmilesConstraint {
    endpoint: AtomId,
    direction_at_endpoint: StereoBondMarkKind,
}

impl SmilesStereoWriteContext {
    fn new(mol: &Molecule) -> std::result::Result<Self, MolWriteError> {
        let mut tetrahedral = BTreeMap::new();
        let mut directional = BTreeMap::<BondId, Vec<DirectionalSmilesConstraint>>::new();
        for (_, element) in mol.stereo_elements() {
            match &element.kind {
                StereoElementKind::Tetrahedral(stereo) => {
                    if tetrahedral
                        .insert(
                            stereo.center,
                            TetrahedralSmilesState {
                                carriers: stereo.carriers.clone(),
                                orientation: stereo.orientation,
                            },
                        )
                        .is_some()
                    {
                        return Err(MolWriteError::new(
                            "isomeric SMILES writer cannot encode multiple tetrahedral elements on one atom",
                        ));
                    }
                }
                StereoElementKind::DoubleBond(stereo) => {
                    add_double_bond_directional_constraints(mol, stereo, &mut directional)?;
                }
                StereoElementKind::Axis(_) => {}
            }
        }
        validate_isomeric_source_marks(mol, &directional)?;
        Ok(Self {
            tetrahedral,
            directional,
        })
    }

    fn atom_chirality(
        &self,
        atom: AtomId,
        parent: Option<AtomId>,
        closures: Option<&[SmilesRingClosure]>,
        children: &[(BondId, BondOrder, AtomId)],
        main_child_index: Option<usize>,
    ) -> Option<std::result::Result<ChiralAtomWriteState, MolWriteError>> {
        let stereo = self.tetrahedral.get(&atom)?;
        Some(tetrahedral_chirality_for_smiles_order(
            stereo,
            parent,
            closures,
            children,
            main_child_index,
        ))
    }

    fn directional_bond(
        &self,
        bond: BondId,
        left: AtomId,
        right: AtomId,
    ) -> std::result::Result<Option<StereoBondMarkKind>, MolWriteError> {
        let Some(constraints) = self.directional.get(&bond) else {
            return Ok(None);
        };
        let mut concrete = None;
        for constraint in constraints {
            let mark = directional_mark_for_emitted_bond(
                constraint.direction_at_endpoint,
                constraint.endpoint,
                left,
                right,
            )?;
            if let Some(previous) = concrete {
                if previous != mark {
                    return Err(MolWriteError::new(
                        "isomeric SMILES writer cannot encode conflicting double-bond stereo constraints",
                    ));
                }
            } else {
                concrete = Some(mark);
            }
        }
        Ok(concrete)
    }
}

fn validate_isomeric_source_marks(
    mol: &Molecule,
    directional: &BTreeMap<BondId, Vec<DirectionalSmilesConstraint>>,
) -> std::result::Result<(), MolWriteError> {
    for mark in mol.stereo_bond_marks() {
        if !matches!(
            mark.kind,
            StereoBondMarkKind::DirectionalUp | StereoBondMarkKind::DirectionalDown
        ) {
            return Err(MolWriteError::new(
                "isomeric SMILES writer cannot encode non-directional source bond marks",
            ));
        }
        if !directional.contains_key(&mark.bond) {
            return Err(MolWriteError::new(
                "isomeric SMILES writer requires perceived double-bond stereo elements for source bond marks",
            ));
        }
    }
    Ok(())
}

fn add_double_bond_directional_constraints(
    mol: &Molecule,
    stereo: &DoubleBondStereo,
    directional: &mut BTreeMap<BondId, Vec<DirectionalSmilesConstraint>>,
) -> std::result::Result<(), MolWriteError> {
    let left_carrier_bond = double_bond_printable_carrier_bond(
        mol,
        stereo.left,
        stereo.right,
        stereo.bond,
        stereo.left_carrier,
    )?;
    let right_carrier_bond = double_bond_printable_carrier_bond(
        mol,
        stereo.right,
        stereo.left,
        stereo.bond,
        stereo.right_carrier,
    )?;
    let left_direction = StereoBondMarkKind::DirectionalUp;
    let right_direction = match stereo.orientation {
        DoubleBondOrientation::Together => left_direction,
        DoubleBondOrientation::Opposite => invert_directional_mark(left_direction),
    };
    directional
        .entry(left_carrier_bond.bond)
        .or_default()
        .push(DirectionalSmilesConstraint {
            endpoint: stereo.left,
            direction_at_endpoint: if left_carrier_bond.invert_direction {
                invert_directional_mark(left_direction)
            } else {
                left_direction
            },
        });
    directional
        .entry(right_carrier_bond.bond)
        .or_default()
        .push(DirectionalSmilesConstraint {
            endpoint: stereo.right,
            direction_at_endpoint: if right_carrier_bond.invert_direction {
                invert_directional_mark(right_direction)
            } else {
                right_direction
            },
        });
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct DoubleBondPrintableCarrierBond {
    bond: BondId,
    invert_direction: bool,
}

fn double_bond_printable_carrier_bond(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
    carrier: StereoCarrier,
) -> std::result::Result<DoubleBondPrintableCarrierBond, MolWriteError> {
    match carrier {
        StereoCarrier::Atom(atom) => {
            let bond = mol
                .bond_between(endpoint, atom)
                .map_err(|error| MolWriteError::new(error.to_string()))?
                .ok_or_else(|| MolWriteError::new("double-bond stereo carrier is not bonded"))?;
            Ok(DoubleBondPrintableCarrierBond {
                bond,
                invert_direction: false,
            })
        }
        StereoCarrier::ImplicitHydrogen => {
            let Some(bond) = implicit_double_bond_printable_carrier_bond(
                mol,
                endpoint,
                other_endpoint,
                focus_bond,
            )?
            else {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode implicit double-bond carrier without a unique explicit substituent bond",
                ));
            };
            Ok(DoubleBondPrintableCarrierBond {
                bond,
                invert_direction: true,
            })
        }
        StereoCarrier::ImplicitLonePair => Err(MolWriteError::new(
            "isomeric SMILES writer cannot encode lone-pair double-bond carrier",
        )),
    }
}

fn implicit_double_bond_printable_carrier_bond(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
) -> std::result::Result<Option<BondId>, MolWriteError> {
    let mut candidates = Vec::new();
    for (bond_id, bond) in mol
        .incident_bonds(endpoint)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        if bond_id == focus_bond || bond.order != BondOrder::Single {
            continue;
        }
        let other = bond.other_atom(endpoint);
        if other != other_endpoint {
            candidates.push(bond_id);
        }
    }
    match candidates.as_slice() {
        [bond] => Ok(Some(*bond)),
        [] => Ok(None),
        _ => Err(MolWriteError::new(
            "isomeric SMILES writer cannot encode implicit double-bond carrier with multiple explicit substituent bonds",
        )),
    }
}

fn directional_mark_for_emitted_bond(
    direction_at_endpoint: StereoBondMarkKind,
    endpoint: AtomId,
    left: AtomId,
    right: AtomId,
) -> std::result::Result<StereoBondMarkKind, MolWriteError> {
    if endpoint == left {
        Ok(direction_at_endpoint)
    } else if endpoint == right {
        Ok(invert_directional_mark(direction_at_endpoint))
    } else {
        Err(MolWriteError::new(
            "double-bond stereo endpoint is not on emitted directional bond",
        ))
    }
}

fn invert_directional_mark(kind: StereoBondMarkKind) -> StereoBondMarkKind {
    match kind {
        StereoBondMarkKind::DirectionalUp => StereoBondMarkKind::DirectionalDown,
        StereoBondMarkKind::DirectionalDown => StereoBondMarkKind::DirectionalUp,
        _ => kind,
    }
}

fn tetrahedral_chirality_for_smiles_order(
    stereo: &TetrahedralSmilesState,
    parent: Option<AtomId>,
    closures: Option<&[SmilesRingClosure]>,
    children: &[(BondId, BondOrder, AtomId)],
    main_child_index: Option<usize>,
) -> std::result::Result<ChiralAtomWriteState, MolWriteError> {
    let force_hydrogen = stereo
        .carriers
        .iter()
        .any(|carrier| matches!(carrier, StereoCarrier::ImplicitHydrogen));
    let mut emitted = Vec::with_capacity(stereo.carriers.len());
    if let Some(parent) = parent {
        emitted.push(StereoCarrier::Atom(parent));
    }
    if force_hydrogen {
        emitted.push(StereoCarrier::ImplicitHydrogen);
    }
    if let Some(closures) = closures {
        emitted.extend(
            closures
                .iter()
                .map(|closure| StereoCarrier::Atom(closure.other)),
        );
    }
    emitted.extend(
        children
            .iter()
            .enumerate()
            .filter(|(index, _)| Some(*index) != main_child_index)
            .map(|(_, (_, _, child))| StereoCarrier::Atom(*child)),
    );
    if let Some(index) = main_child_index {
        emitted.push(StereoCarrier::Atom(children[index].2));
    }
    if stereo
        .carriers
        .iter()
        .any(|carrier| matches!(carrier, StereoCarrier::ImplicitLonePair))
    {
        emitted.push(StereoCarrier::ImplicitLonePair);
    }
    if emitted != stereo.carriers {
        let Some(odd) = carrier_permutation_is_odd(&stereo.carriers, &emitted) else {
            return Err(MolWriteError::new(
                "isomeric SMILES writer cannot encode tetrahedral carrier order",
            ));
        };
        Ok(ChiralAtomWriteState {
            orientation: if odd {
                flip_tetrahedral_orientation(stereo.orientation)
            } else {
                stereo.orientation
            },
            force_hydrogen,
        })
    } else {
        Ok(ChiralAtomWriteState {
            orientation: stereo.orientation,
            force_hydrogen,
        })
    }
}

fn carrier_permutation_is_odd(from: &[StereoCarrier], to: &[StereoCarrier]) -> Option<bool> {
    if from.len() != to.len() {
        return None;
    }
    let mut positions = Vec::with_capacity(to.len());
    let mut used = vec![false; to.len()];
    for carrier in from {
        let position = to
            .iter()
            .enumerate()
            .find(|(index, candidate)| !used[*index] && *candidate == carrier)
            .map(|(index, _)| index)?;
        used[position] = true;
        positions.push(position);
    }
    let mut odd = false;
    for left in 0..positions.len() {
        for right in (left + 1)..positions.len() {
            if positions[left] > positions[right] {
                odd = !odd;
            }
        }
    }
    Some(odd)
}

fn flip_tetrahedral_orientation(orientation: TetrahedralOrientation) -> TetrahedralOrientation {
    match orientation {
        TetrahedralOrientation::Clockwise => TetrahedralOrientation::CounterClockwise,
        TetrahedralOrientation::CounterClockwise => TetrahedralOrientation::Clockwise,
    }
}

fn write_smiles_component(
    mol: &Molecule,
    atom_id: AtomId,
    parent: Option<AtomId>,
    plan: &SmilesWritePlan,
    stereo: Option<&SmilesStereoWriteContext>,
) -> std::result::Result<String, MolWriteError> {
    enum Action {
        Node {
            atom: AtomId,
            parent: Option<AtomId>,
        },
        Bond {
            bond: BondId,
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
            Action::Bond {
                bond,
                order,
                left,
                right,
            } => {
                let directional = stereo
                    .map(|context| context.directional_bond(bond, left, right))
                    .transpose()?
                    .flatten();
                out.push_str(smiles_bond_between_with_direction(
                    mol,
                    order,
                    left,
                    right,
                    directional,
                )?);
            }
            Action::Node { atom, parent } => {
                let atom_record = mol
                    .atom(atom)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                let closures = plan.closures.get(&atom).map(Vec::as_slice);
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
                let chirality = stereo
                    .and_then(|context| {
                        context.atom_chirality(atom, parent, closures, &children, main_child_index)
                    })
                    .transpose()?;
                out.push_str(&smiles_atom_with_chirality(
                    atom_record,
                    chirality.map(|state| state.orientation),
                    chirality.is_some_and(|state| state.force_hydrogen),
                ));
                if let Some(closures) = closures {
                    for closure in closures {
                        let directional = stereo
                            .map(|context| {
                                context.directional_bond(closure.bond, atom, closure.other)
                            })
                            .transpose()?
                            .flatten();
                        out.push_str(smiles_bond_between_with_direction(
                            mol,
                            closure.order,
                            atom,
                            closure.other,
                            directional,
                        )?);
                        out.push_str(&smiles_ring_number(closure.number));
                    }
                }

                if let Some(index) = main_child_index {
                    let (bond, order, child) = children[index];
                    actions.push(Action::Node {
                        atom: child,
                        parent: Some(atom),
                    });
                    actions.push(Action::Bond {
                        bond,
                        order,
                        left: atom,
                        right: child,
                    });
                }
                for (index, (bond, order, child)) in children.into_iter().enumerate().rev() {
                    if Some(index) == main_child_index {
                        continue;
                    }
                    actions.push(Action::CloseBranch);
                    actions.push(Action::Node {
                        atom: child,
                        parent: Some(atom),
                    });
                    actions.push(Action::Bond {
                        bond,
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
    smiles_incident_bonds_for_style(mol, atom_id, CanonicalAtomStyle::Aromatic)
}

fn smiles_incident_bonds_for_style(
    mol: &Molecule,
    atom_id: AtomId,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<Vec<(BondId, BondOrder, AtomId)>, MolWriteError> {
    let mut incident = Vec::new();
    for (bond_id, bond) in mol
        .incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        let order = match atom_style {
            CanonicalAtomStyle::Aromatic if bond.aromatic => BondOrder::Aromatic,
            CanonicalAtomStyle::Aromatic | CanonicalAtomStyle::StoredKekule => bond.order,
        };
        incident.push((bond_id, order, bond.other_atom(atom_id)));
    }
    incident.sort_by_key(|(bond_id, _, atom)| (*atom, *bond_id));
    Ok(incident)
}

fn canonical_smiles_incident_bonds(
    mol: &Molecule,
    atom_id: AtomId,
    ranking: &CanonicalAtomRanking,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<Vec<(BondId, BondOrder, AtomId)>, MolWriteError> {
    let mut incident = smiles_incident_bonds_for_style(mol, atom_id, atom_style)?;
    incident.sort_by_key(|(bond_id, order, atom)| {
        (
            canonical_rank(ranking, *atom),
            canonical_smiles_atom_for_sort(mol, *atom, atom_style),
            preference.order_key(*order),
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

fn reverse_bond_order_code(order: BondOrder) -> u8 {
    u8::MAX - bond_order_code(order)
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

fn smiles_bond_between_with_direction(
    mol: &Molecule,
    order: BondOrder,
    left: AtomId,
    right: AtomId,
    directional: Option<StereoBondMarkKind>,
) -> std::result::Result<&'static str, MolWriteError> {
    if let Some(directional) = directional {
        if order != BondOrder::Single {
            return Err(MolWriteError::new(
                "isomeric SMILES writer cannot place directional stereo on a non-single bond",
            ));
        }
        return match directional {
            StereoBondMarkKind::DirectionalUp => Ok("/"),
            StereoBondMarkKind::DirectionalDown => Ok("\\"),
            _ => Err(MolWriteError::new(
                "isomeric SMILES writer received a non-directional stereo mark",
            )),
        };
    }
    smiles_bond_between(mol, order, left, right)
}

fn smiles_atom(atom: &Atom) -> String {
    smiles_atom_with_chirality(atom, None, false)
}

fn smiles_atom_with_chirality(
    atom: &Atom,
    chirality: Option<TetrahedralOrientation>,
    force_hydrogen: bool,
) -> String {
    let explicit_hydrogens = if force_hydrogen {
        smiles_atom_explicit_hydrogens(atom).max(1)
    } else {
        smiles_atom_explicit_hydrogens(atom)
    };
    let organic = atom.isotope.is_none()
        && atom.formal_charge == 0
        && explicit_hydrogens == 0
        && !atom.no_implicit_hydrogens
        && atom.atom_map.is_none()
        && chirality.is_none()
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
        if let Some(chirality) = chirality {
            out.push('@');
            if chirality == TetrahedralOrientation::CounterClockwise {
                out.push('@');
            }
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

fn canonical_smiles_atom(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<String, MolWriteError> {
    if matches!(atom_style, CanonicalAtomStyle::StoredKekule) && atom.aromatic {
        let mut normalized = atom.clone();
        normalized.isotope = None;
        normalized.aromatic = false;
        if !matches!(atom.element.symbol(), "B" | "C") && atom.implicit_hydrogens.unwrap_or(0) > 0 {
            normalized.explicit_hydrogens = atom
                .explicit_hydrogens
                .saturating_add(atom.implicit_hydrogens.unwrap_or(0));
            normalized.implicit_hydrogens = Some(0);
            normalized.no_implicit_hydrogens = true;
        }
        return Ok(smiles_atom(&normalized));
    }
    if canonical_smiles_should_bracket_metal_bound_hydrogens(mol, atom_id, atom)? {
        let mut normalized = atom.clone();
        normalized.isotope = None;
        normalized.explicit_hydrogens = atom.implicit_hydrogens.unwrap_or(0);
        normalized.implicit_hydrogens = Some(0);
        normalized.no_implicit_hydrogens = true;
        return Ok(smiles_atom(&normalized));
    }
    if canonical_smiles_should_bracket_metal_bound_zero_hydrogens(mol, atom_id, atom)? {
        let mut normalized = atom.clone();
        normalized.isotope = None;
        normalized.implicit_hydrogens = Some(0);
        normalized.no_implicit_hydrogens = true;
        return Ok(smiles_atom(&normalized));
    }
    if canonical_smiles_can_use_organic_form(mol, atom_id, atom)? {
        let mut normalized = atom.clone();
        normalized.isotope = None;
        normalized.explicit_hydrogens = 0;
        normalized.no_implicit_hydrogens = false;
        return Ok(smiles_atom(&normalized));
    }
    let mut normalized = atom.clone();
    normalized.isotope = None;
    Ok(smiles_atom(&normalized))
}

fn canonical_smiles_should_bracket_metal_bound_hydrogens(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
) -> std::result::Result<bool, MolWriteError> {
    Ok(atom.formal_charge == 0
        && atom.radical.is_none()
        && atom.atom_map.is_none()
        && !atom.aromatic
        && !atom.no_implicit_hydrogens
        && atom.explicit_hydrogens == 0
        && atom.implicit_hydrogens.unwrap_or(0) > 0
        && matches!(atom.element.symbol(), "B" | "C" | "N" | "O" | "P" | "S")
        && atom_has_metal_neighbor(mol, atom_id)?)
}

fn canonical_smiles_should_bracket_metal_bound_zero_hydrogens(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
) -> std::result::Result<bool, MolWriteError> {
    Ok(atom.formal_charge == 0
        && atom.radical.is_none()
        && atom.atom_map.is_none()
        && atom.explicit_hydrogens == 0
        && atom.implicit_hydrogens == Some(0)
        && matches!(
            atom.element.symbol(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        )
        && atom_has_metal_neighbor(mol, atom_id)?)
}

fn canonical_smiles_atom_for_sort(
    mol: &Molecule,
    atom_id: AtomId,
    atom_style: CanonicalAtomStyle,
) -> String {
    let atom = mol
        .atom(atom_id)
        .expect("canonical atom sort should only use live atoms");
    canonical_smiles_atom(mol, atom_id, atom, atom_style)
        .expect("canonical atom sort should be encodable")
}

fn canonical_smiles_can_use_organic_form(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
) -> std::result::Result<bool, MolWriteError> {
    if atom.formal_charge != 0
        || atom.radical.is_some()
        || atom.atom_map.is_some()
        || (atom.aromatic && atom.explicit_hydrogens > 0)
    {
        return Ok(false);
    }
    if !matches!(
        atom.element.symbol(),
        "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
    ) {
        return Ok(false);
    }
    let Some(target) = canonical_organic_valence_target(atom) else {
        return Ok(false);
    };
    if (atom.no_implicit_hydrogens || atom.implicit_hydrogens == Some(0))
        && atom_has_metal_neighbor(mol, atom_id)?
    {
        return Ok(false);
    }
    let bond_valence = smiles_bond_valence_sum(mol, atom_id)?;
    Ok(bond_valence.saturating_add(atom.explicit_hydrogens) == target)
}

fn atom_has_metal_neighbor(
    mol: &Molecule,
    atom_id: AtomId,
) -> std::result::Result<bool, MolWriteError> {
    for (_, bond) in mol
        .incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        let neighbor_id = bond.other_atom(atom_id);
        let neighbor = mol
            .atom(neighbor_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        if is_smiles_metal_like(neighbor.element.symbol()) {
            return Ok(true);
        }
    }
    Ok(false)
}

fn is_smiles_metal_like(symbol: &str) -> bool {
    matches!(
        symbol,
        "Li" | "Na"
            | "K"
            | "Rb"
            | "Cs"
            | "Fr"
            | "Be"
            | "Mg"
            | "Ca"
            | "Sr"
            | "Ba"
            | "Ra"
            | "Al"
            | "Ge"
            | "Ga"
            | "In"
            | "Tl"
            | "Sn"
            | "Pb"
            | "Bi"
            | "Po"
            | "Sc"
            | "Ti"
            | "V"
            | "Cr"
            | "Mn"
            | "Fe"
            | "Co"
            | "Ni"
            | "Cu"
            | "Zn"
            | "Y"
            | "Zr"
            | "Nb"
            | "Mo"
            | "Tc"
            | "Ru"
            | "Rh"
            | "Pd"
            | "Ag"
            | "Cd"
            | "La"
            | "Ce"
            | "Pr"
            | "Nd"
            | "Sm"
            | "Eu"
            | "Gd"
            | "Tb"
            | "Dy"
            | "Ho"
            | "Er"
            | "Tm"
            | "Yb"
            | "Lu"
            | "Ac"
            | "Th"
            | "Pa"
            | "U"
            | "Np"
            | "Pu"
            | "Am"
            | "Cm"
            | "Bk"
            | "Cf"
            | "Es"
            | "Fm"
            | "Md"
            | "No"
            | "Lr"
            | "Hf"
            | "Ta"
            | "W"
            | "Re"
            | "Os"
            | "Ir"
            | "Pt"
            | "Au"
            | "Hg"
    )
}

fn canonical_organic_valence_target(atom: &Atom) -> Option<u8> {
    match (atom.element.symbol(), atom.aromatic) {
        ("B", false) => Some(3),
        ("C", false) => Some(4),
        ("N", false) | ("P", false) => Some(3),
        ("O", false) | ("S", false) => Some(2),
        ("F" | "Cl" | "Br" | "I", false) => Some(1),
        ("B" | "C", true) => Some(3),
        ("N" | "O" | "S" | "P", true) => Some(3),
        _ => None,
    }
}

fn smiles_bond_valence_sum(
    mol: &Molecule,
    atom_id: AtomId,
) -> std::result::Result<u8, MolWriteError> {
    mol.incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
        .map(|(_, bond)| {
            Ok(match bond.order {
                BondOrder::Zero | BondOrder::Dative => 0,
                BondOrder::Single | BondOrder::Aromatic => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Quadruple => 4,
            })
        })
        .try_fold(0u8, |sum, value: std::result::Result<u8, MolWriteError>| {
            Ok(sum.saturating_add(value?))
        })
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
