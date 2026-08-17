use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::Range;

#[derive(Debug, Clone, PartialEq)]
pub struct SmilesDocument {
    source: String,
    tokens: Vec<SmilesDocumentToken>,
    components: Vec<Range<usize>>,
    pub(super) program: SmilesProgram,
}

impl SmilesDocument {
    pub fn source(&self) -> &str {
        &self.source
    }

    pub fn tokens(&self) -> &[SmilesDocumentToken] {
        &self.tokens
    }

    pub fn component_token_ranges(&self) -> &[Range<usize>] {
        &self.components
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct SmilesProgram {
    pub(super) atoms: Vec<SmilesProgramAtom>,
    pub(super) bonds: Vec<SmilesProgramBond>,
    pub(super) imported_aromatic_atoms: BTreeSet<usize>,
    pub(super) tetrahedral: Vec<PendingTetrahedral>,
    pub(super) tetrahedral_carriers: BTreeMap<usize, Vec<PendingStereoCarrier>>,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct SmilesProgramAtom {
    pub(super) syntax: SmilesAtomSyntax,
    pub(super) span: Range<usize>,
    pub(super) component: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SmilesAtomSyntax {
    pub(super) symbol: String,
    pub(super) isotope: Option<u16>,
    pub(super) explicit_hydrogens: u8,
    pub(super) formal_charge: i8,
    pub(super) atom_map: Option<u32>,
    pub(super) aromatic: bool,
    pub(super) bracketed: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SmilesProgramBond {
    pub(super) left: usize,
    pub(super) right: usize,
    pub(super) token: SmilesBondToken,
    pub(super) direction: Option<SmilesDirectionToken>,
    pub(super) offset: usize,
    pub(super) component: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SmilesBondToken {
    Single,
    Double,
    Triple,
    Aromatic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SmilesDirectionToken {
    Up,
    Down,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SmilesChiralityToken {
    At,
    AtAt,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesDocumentToken {
    kind: SmilesDocumentTokenKind,
    span: Range<usize>,
}

impl SmilesDocumentToken {
    pub const fn kind(&self) -> SmilesDocumentTokenKind {
        self.kind
    }

    pub fn span(&self) -> Range<usize> {
        self.span.clone()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SmilesDocumentTokenKind {
    Atom,
    Bond,
    BranchOpen,
    BranchClose,
    Ring,
    ComponentSeparator,
    Unsupported,
}

/// Resource limits for parsing one SMILES record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmilesParseOptions {
    /// Maximum UTF-8 byte length of the complete SMILES record.
    pub max_input_bytes: usize,
    /// Maximum parsed atom count.
    pub max_atoms: usize,
    /// Maximum parsed bond count.
    pub max_bonds: usize,
}

impl Default for SmilesParseOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 * 1024 * 1024,
            max_atoms: 1_000_000,
            max_bonds: 2_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesParseError {
    offset: usize,
    message: String,
}

impl SmilesParseError {
    fn new(offset: usize, message: impl Into<String>) -> Self {
        Self {
            offset,
            message: message.into(),
        }
    }

    pub const fn offset(&self) -> usize {
        self.offset
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for SmilesParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SMILES parse error at {}: {}", self.offset, self.message)
    }
}

impl std::error::Error for SmilesParseError {}

pub fn parse_smiles_document(input: &str) -> std::result::Result<SmilesDocument, SmilesParseError> {
    parse_smiles_document_with_options(input, SmilesParseOptions::default())
}

pub fn parse_smiles_document_with_options(
    input: &str,
    options: SmilesParseOptions,
) -> std::result::Result<SmilesDocument, SmilesParseError> {
    if options.max_atoms > u32::MAX as usize || options.max_bonds > u32::MAX as usize {
        return Err(SmilesParseError::new(
            0,
            "configured atom and bond limits must fit the public identifier range",
        ));
    }
    if input.len() > options.max_input_bytes {
        return Err(SmilesParseError::new(
            0,
            "input exceeds configured byte limit",
        ));
    }
    if input.is_empty() {
        return Err(SmilesParseError::new(0, "empty SMILES document"));
    }
    let chars = input.char_indices().collect::<Vec<_>>();
    let mut tokens = Vec::new();
    let mut components = Vec::new();
    let mut component_start = 0usize;
    let mut cursor = 0usize;
    let mut branch_depth = 0usize;
    while cursor < chars.len() {
        let (start, ch) = chars[cursor];
        let (kind, next) = match ch {
            '[' => {
                let mut end = cursor + 1;
                while end < chars.len() && chars[end].1 != ']' {
                    end += 1;
                }
                if end == chars.len() {
                    return Err(SmilesParseError::new(start, "unclosed bracket atom"));
                }
                (SmilesDocumentTokenKind::Atom, end + 1)
            }
            '(' => {
                branch_depth += 1;
                (SmilesDocumentTokenKind::BranchOpen, cursor + 1)
            }
            ')' => {
                if branch_depth == 0 {
                    return Err(SmilesParseError::new(start, "unmatched branch close"));
                }
                branch_depth -= 1;
                (SmilesDocumentTokenKind::BranchClose, cursor + 1)
            }
            '.' => {
                if component_start == tokens.len() {
                    return Err(SmilesParseError::new(start, "empty SMILES component"));
                }
                components.push(component_start..tokens.len());
                component_start = tokens.len() + 1;
                (SmilesDocumentTokenKind::ComponentSeparator, cursor + 1)
            }
            '-' | '=' | '#' | ':' | '/' | '\\' => (SmilesDocumentTokenKind::Bond, cursor + 1),
            '0'..='9' => (SmilesDocumentTokenKind::Ring, cursor + 1),
            '%' => {
                let mut next = cursor + 1;
                while next < chars.len() && chars[next].1.is_ascii_digit() {
                    next += 1;
                }
                if next == cursor + 1 {
                    return Err(SmilesParseError::new(start, "invalid ring label"));
                }
                (SmilesDocumentTokenKind::Ring, next)
            }
            '*' | '@' => (SmilesDocumentTokenKind::Unsupported, cursor + 1),
            ch if ch.is_ascii_alphabetic() => {
                let next = if matches!(ch, 'C' | 'B')
                    && chars
                        .get(cursor + 1)
                        .is_some_and(|(_, next)| matches!((ch, *next), ('C', 'l') | ('B', 'r')))
                {
                    cursor + 2
                } else {
                    cursor + 1
                };
                (SmilesDocumentTokenKind::Atom, next)
            }
            _ => {
                return Err(SmilesParseError::new(
                    start,
                    format!("invalid SMILES syntax character `{ch}`"),
                ));
            }
        };
        let end = chars
            .get(next)
            .map(|(offset, _)| *offset)
            .unwrap_or(input.len());
        tokens.push(SmilesDocumentToken {
            kind,
            span: start..end,
        });
        cursor = next;
    }
    if branch_depth != 0 {
        return Err(SmilesParseError::new(input.len(), "unclosed branch"));
    }
    if component_start == tokens.len() {
        return Err(SmilesParseError::new(
            input.len(),
            "empty trailing SMILES component",
        ));
    }
    components.push(component_start..tokens.len());
    let program = parse_smiles_program(input, &chars, options)?;
    Ok(SmilesDocument {
        source: input.to_owned(),
        tokens,
        components,
        program,
    })
}

fn parse_smiles_program(
    input: &str,
    chars: &[(usize, char)],
    options: SmilesParseOptions,
) -> std::result::Result<SmilesProgram, SmilesParseError> {
    let mut atoms = Vec::<SmilesProgramAtom>::new();
    let mut bonds = Vec::<SmilesProgramBond>::new();
    let mut imported_aromatic_atoms = BTreeSet::new();
    let mut current: Option<usize> = None;
    let mut stack = Vec::<usize>::new();
    let mut pending_bond = None::<(SmilesBondToken, Option<SmilesDirectionToken>, usize)>;
    let mut rings = BTreeMap::<
        usize,
        (
            usize,
            Option<(SmilesBondToken, Option<SmilesDirectionToken>)>,
            usize,
        ),
    >::new();
    let mut pending_tetrahedral = Vec::<PendingTetrahedral>::new();
    let mut tetrahedral_carriers = BTreeMap::<usize, Vec<PendingStereoCarrier>>::new();
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
                let token = match ch {
                    '-' => SmilesBondToken::Single,
                    '=' => SmilesBondToken::Double,
                    '#' => SmilesBondToken::Triple,
                    ':' => SmilesBondToken::Aromatic,
                    '/' | '\\' => SmilesBondToken::Single,
                    _ => unreachable!(),
                };
                let direction = match ch {
                    '/' => Some(SmilesDirectionToken::Up),
                    '\\' => Some(SmilesDirectionToken::Down),
                    _ => None,
                };
                pending_bond = Some((token, direction, offset));
                previous = SmilesTokenKind::Bond;
                cursor += 1;
            }
            '0'..='9' | '%' => {
                let atom = current
                    .ok_or_else(|| SmilesParseError::new(offset, "ring closure without atom"))?;
                let (label, next_cursor) = parse_smiles_ring_label(chars, cursor)?;
                let close_bond = pending_bond
                    .take()
                    .map(|(token, direction, _)| (token, direction));
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
                    let (token, direction) = match close_bond.or(open_bond) {
                        Some((token, direction)) => (token, direction),
                        None => (
                            default_smiles_bond_order(&imported_aromatic_atoms, other, atom),
                            None,
                        ),
                    };
                    add_smiles_program_bond(
                        &mut bonds,
                        (other, atom),
                        token,
                        direction,
                        offset,
                        component,
                        options.max_bonds,
                    )?;
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
                        SmilesStereoCarrier::Atom(other),
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
                let (atom, aromatic, chirality, next_cursor) = parse_bracket_atom(chars, cursor)?;
                let explicit_hydrogens = atom.explicit_hydrogens;
                let atom_id = next_smiles_atom_index(atoms.len(), offset, options.max_atoms)?;
                let end = chars
                    .get(next_cursor)
                    .map(|(offset, _)| *offset)
                    .unwrap_or(input.len());
                atoms.push(SmilesProgramAtom {
                    syntax: atom,
                    span: offset..end,
                    component,
                });
                if aromatic {
                    imported_aromatic_atoms.insert(atom_id);
                }
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
                    let (token, direction) = match pending_bond
                        .take()
                        .map(|(token, direction, _)| (token, direction))
                    {
                        Some((token, direction)) => (token, direction),
                        None => (
                            default_smiles_bond_order(&imported_aromatic_atoms, previous, atom_id),
                            None,
                        ),
                    };
                    add_smiles_program_bond(
                        &mut bonds,
                        (previous, atom_id),
                        token,
                        direction,
                        offset,
                        component,
                        options.max_bonds,
                    )?;
                    push_tetrahedral_carrier(
                        &mut tetrahedral_carriers,
                        previous,
                        SmilesStereoCarrier::Atom(atom_id),
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
                let (atom, aromatic, next_cursor) = parse_organic_atom(chars, cursor)?;
                let atom_id = next_smiles_atom_index(atoms.len(), offset, options.max_atoms)?;
                let end = chars
                    .get(next_cursor)
                    .map(|(offset, _)| *offset)
                    .unwrap_or(input.len());
                atoms.push(SmilesProgramAtom {
                    syntax: atom,
                    span: offset..end,
                    component,
                });
                if aromatic {
                    imported_aromatic_atoms.insert(atom_id);
                }
                if let Some(previous) = current {
                    let (token, direction) = match pending_bond
                        .take()
                        .map(|(token, direction, _)| (token, direction))
                    {
                        Some((token, direction)) => (token, direction),
                        None => (
                            default_smiles_bond_order(&imported_aromatic_atoms, previous, atom_id),
                            None,
                        ),
                    };
                    add_smiles_program_bond(
                        &mut bonds,
                        (previous, atom_id),
                        token,
                        direction,
                        offset,
                        component,
                        options.max_bonds,
                    )?;
                    push_tetrahedral_carrier(
                        &mut tetrahedral_carriers,
                        previous,
                        SmilesStereoCarrier::Atom(atom_id),
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
    Ok(SmilesProgram {
        atoms,
        bonds,
        imported_aromatic_atoms,
        tetrahedral: pending_tetrahedral,
        tetrahedral_carriers,
    })
}

fn next_smiles_atom_index(
    atom_count: usize,
    offset: usize,
    max_atoms: usize,
) -> std::result::Result<usize, SmilesParseError> {
    if atom_count >= max_atoms {
        return Err(SmilesParseError::new(
            offset,
            "SMILES atom count exceeds configured limit",
        ));
    }
    if atom_count > u32::MAX as usize {
        return Err(SmilesParseError::new(
            offset,
            "SMILES atom index capacity exceeded",
        ));
    }
    Ok(atom_count)
}

fn add_smiles_program_bond(
    bonds: &mut Vec<SmilesProgramBond>,
    endpoints: (usize, usize),
    token: SmilesBondToken,
    direction: Option<SmilesDirectionToken>,
    offset: usize,
    component: usize,
    max_bonds: usize,
) -> std::result::Result<(), SmilesParseError> {
    let (left, right) = endpoints;
    if bonds.len() >= max_bonds {
        return Err(SmilesParseError::new(
            offset,
            "SMILES bond count exceeds configured limit",
        ));
    }
    if left == right {
        return Err(SmilesParseError::new(offset, "self bond"));
    }
    let endpoints = if left < right {
        (left, right)
    } else {
        (right, left)
    };
    if bonds.iter().any(|bond| {
        let bond_endpoints = if bond.left < bond.right {
            (bond.left, bond.right)
        } else {
            (bond.right, bond.left)
        };
        bond_endpoints == endpoints
    }) {
        return Err(SmilesParseError::new(offset, "duplicate bond"));
    }
    bonds.push(SmilesProgramBond {
        left,
        right,
        token,
        direction,
        offset,
        component,
    });
    Ok(())
}

fn initial_tetrahedral_carriers(
    previous: Option<usize>,
    explicit_hydrogens: u8,
) -> Vec<PendingStereoCarrier> {
    let mut carriers = Vec::new();
    if let Some(previous) = previous {
        carriers.push(PendingStereoCarrier::Resolved(SmilesStereoCarrier::Atom(
            previous,
        )));
    }
    for _ in 0..explicit_hydrogens {
        carriers.push(PendingStereoCarrier::Resolved(
            SmilesStereoCarrier::ImplicitHydrogen,
        ));
    }
    carriers
}

fn push_tetrahedral_carrier(
    carriers_by_center: &mut BTreeMap<usize, Vec<PendingStereoCarrier>>,
    center: usize,
    carrier: SmilesStereoCarrier,
) {
    if let Some(carriers) = carriers_by_center.get_mut(&center) {
        carriers.push(PendingStereoCarrier::Resolved(carrier));
    }
}

fn push_tetrahedral_ring_carrier(
    carriers_by_center: &mut BTreeMap<usize, Vec<PendingStereoCarrier>>,
    center: usize,
    label: usize,
    component: usize,
) {
    if let Some(carriers) = carriers_by_center.get_mut(&center) {
        carriers.push(PendingStereoCarrier::Ring { label, component });
    }
}

fn resolve_tetrahedral_ring_carrier(
    carriers_by_center: &mut BTreeMap<usize, Vec<PendingStereoCarrier>>,
    center: usize,
    label: usize,
    component: usize,
    carrier: usize,
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
        *pending = PendingStereoCarrier::Resolved(SmilesStereoCarrier::Atom(carrier));
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct PendingTetrahedral {
    pub(super) center: usize,
    pub(super) orientation: SmilesChiralityToken,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum PendingStereoCarrier {
    Resolved(SmilesStereoCarrier),
    Ring { label: usize, component: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SmilesStereoCarrier {
    Atom(usize),
    ImplicitHydrogen,
}

fn default_smiles_bond_order(
    imported_aromatic_atoms: &BTreeSet<usize>,
    left: usize,
    right: usize,
) -> SmilesBondToken {
    if imported_aromatic_atoms.contains(&left) && imported_aromatic_atoms.contains(&right) {
        SmilesBondToken::Aromatic
    } else {
        SmilesBondToken::Single
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
) -> std::result::Result<(SmilesAtomSyntax, bool, usize), SmilesParseError> {
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
    Ok((
        SmilesAtomSyntax {
            symbol,
            isotope: None,
            explicit_hydrogens: 0,
            formal_charge: 0,
            atom_map: None,
            aromatic,
            bracketed: false,
        },
        aromatic,
        next,
    ))
}

fn parse_bracket_atom(
    chars: &[(usize, char)],
    cursor: usize,
) -> std::result::Result<
    (SmilesAtomSyntax, bool, Option<SmilesChiralityToken>, usize),
    SmilesParseError,
> {
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
    let mut explicit_hydrogens = 0;
    let mut formal_charge = 0;
    let mut atom_map = None;
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
                    Some(SmilesChiralityToken::AtAt)
                } else {
                    Some(SmilesChiralityToken::At)
                };
            }
            b'H' if !saw_hydrogen && !saw_charge && !saw_map => {
                saw_hydrogen = true;
                index += 1;
                let digit_end = ascii_digits_end(bytes, index);
                explicit_hydrogens = if digit_end == index {
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
                formal_charge = i8::try_from(charge).map_err(|_| {
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
                atom_map = Some(map);
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
    Ok((
        SmilesAtomSyntax {
            symbol: canonical_symbol,
            isotope,
            explicit_hydrogens,
            formal_charge,
            atom_map,
            aromatic,
            bracketed: true,
        },
        aromatic,
        chirality,
        end + 1,
    ))
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
