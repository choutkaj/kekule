use std::cell::Cell;
use std::collections::BTreeMap;
use std::fmt;
use std::ops::Range;
use std::rc::Rc;

use crate::core::{BondOrder, Element};

use super::{
    AtomExpression, AtomPredicate, BondExpression, BondPredicate, QueryAtomId,
    QueryExpressionError, QueryGraph, QueryGraphBuilder, QueryGraphError,
};

mod bond;
mod stereo;
use stereo::{AtomStereo, BondSyntax, StereoSyntax};

/// Resource bounds for deterministic SMARTS parsing.
///
/// The defaults accept ordinary small-molecule queries while bounding input,
/// graph, branch, ring-closure, and expression complexity before allocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmartsParseOptions {
    /// Maximum UTF-8 byte length; the bounded grammar itself is ASCII.
    pub max_input_bytes: usize,
    /// Maximum query atoms.
    pub max_atoms: usize,
    /// Maximum query bonds, including ring closures.
    pub max_bonds: usize,
    /// Maximum nested branch depth.
    pub max_branch_depth: usize,
    /// Maximum simultaneously open and total completed ring closures.
    pub max_ring_closures: usize,
    /// Maximum nodes in each atom expression after normalization.
    pub max_expression_nodes: usize,
    /// Maximum depth of each atom expression after normalization.
    pub max_expression_depth: usize,
    /// Maximum nested recursive SMARTS depth (hard ceiling 32).
    pub max_recursive_depth: usize,
    /// Aggregate atom and bond expression nodes, including recursive queries.
    pub max_total_expression_nodes: usize,
}

impl Default for SmartsParseOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: 16_384,
            max_atoms: 256,
            max_bonds: 512,
            max_branch_depth: 64,
            max_ring_closures: 128,
            max_expression_nodes: 512,
            max_expression_depth: 32,
            max_recursive_depth: 32,
            max_total_expression_nodes: 16_384,
        }
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SmartsParseErrorKind {
    Empty,
    InvalidSyntax,
    Unsupported,
    ResourceLimit,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmartsParseError {
    kind: SmartsParseErrorKind,
    span: Range<usize>,
    message: String,
}

impl SmartsParseError {
    pub const fn kind(&self) -> SmartsParseErrorKind {
        self.kind
    }

    pub fn span(&self) -> Range<usize> {
        self.span.clone()
    }

    pub fn message(&self) -> &str {
        &self.message
    }

    fn new(kind: SmartsParseErrorKind, span: Range<usize>, message: impl Into<String>) -> Self {
        Self {
            kind,
            span,
            message: message.into(),
        }
    }

    fn syntax(span: Range<usize>, message: impl Into<String>) -> Self {
        Self::new(SmartsParseErrorKind::InvalidSyntax, span, message)
    }

    fn unsupported(span: Range<usize>, message: impl Into<String>) -> Self {
        Self::new(SmartsParseErrorKind::Unsupported, span, message)
    }

    fn limit(span: Range<usize>, resource: &'static str, observed: usize, limit: usize) -> Self {
        Self::new(
            SmartsParseErrorKind::ResourceLimit,
            span,
            format!("SMARTS {resource} limit exceeded: observed {observed}, limit {limit}"),
        )
    }
}

impl fmt::Display for SmartsParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SMARTS parse error at {}..{}: {}",
            self.span.start, self.span.end, self.message
        )
    }
}

impl std::error::Error for SmartsParseError {}

/// Parses one SMARTS pattern into a syntax-independent [`QueryGraph`].
///
/// This uses [`SmartsParseOptions::default`]. Parsing does not match a target or
/// run chemical perception.
/// `X` counts all connections, including declared and inferred hydrogens, and
/// defaults to one. `xN` counts incident cyclic bonds; bare `x` means any ring
/// membership. Matching these requires installed valence and ring perception,
/// respectively. Neither predicate selects or counts a ring basis.
/// Single and double bond types exclude perceived aromatic bonds. Aromatic
/// bond types require aromatic membership on a localized single or double bond;
/// higher orders retain their own type. Matching single, double, or aromatic
/// SMARTS bond types requires installed aromaticity. A programmatic
/// [`BondPredicate::Order`] alone only inspects represented order.
///
/// Tetrahedral `@`, `@@`, `@TH1`, and `@TH2` and paired `/`/`\\` bond
/// directions retain their Boolean context and compare represented carrier order
/// under a complete atom mapping, without CIP assignment. An isolated direction
/// adds no stereo relationship. Unspecified and non-tetrahedral stereo syntax
/// remain unsupported. Boolean stereo is evaluated literally; this intentionally
/// differs from RDKit's handling of some negations and alternatives.
/// An optional separated CX extension supports tetrahedral `a`, `oN`, `&N`
/// groups and the `r` flag, with zero-based query atom indices. Other CX fields
/// are rejected. Enhanced groups compare correlated configurations during matching.
/// See the crate's SMARTS capability matrix for the complete dialect contract.
pub fn parse_smarts(input: &str) -> Result<QueryGraph, SmartsParseError> {
    parse_smarts_with_options(input, SmartsParseOptions::default())
}

/// Parses one SMARTS pattern with explicit resource limits.
pub fn parse_smarts_with_options(
    input: &str,
    options: SmartsParseOptions,
) -> Result<QueryGraph, SmartsParseError> {
    validate_options(options)?;
    if input.is_empty() {
        return Err(SmartsParseError::new(
            SmartsParseErrorKind::Empty,
            0..0,
            "empty SMARTS query",
        ));
    }
    if input.len() > options.max_input_bytes {
        return Err(SmartsParseError::limit(
            0..input.len(),
            "input bytes",
            input.len(),
            options.max_input_bytes,
        ));
    }
    if let Some((offset, ch)) = input.char_indices().find(|(_, ch)| !ch.is_ascii()) {
        return Err(SmartsParseError::unsupported(
            offset..offset + ch.len_utf8(),
            "non-ASCII query syntax is outside the bounded SMARTS subset",
        ));
    }
    if let Some(start) = input.find('|') {
        if start == 0
            || start + 1 >= input.len()
            || !input.as_bytes()[start - 1].is_ascii_whitespace()
            || !input.ends_with('|')
        {
            return Err(SmartsParseError::syntax(
                start..input.len(),
                "CXSMARTS requires a separated, terminated extension",
            ));
        }
        let query = Parser::new(
            input[..start].trim_end(),
            options,
            Rc::new(ParseBudget::default()),
            0,
        )
        .parse()?;
        return install_cx_groups(
            query,
            &input[start + 1..input.len() - 1],
            start..input.len(),
        );
    }
    Parser::new(input, options, Rc::new(ParseBudget::default()), 0).parse()
}

fn install_cx_groups(
    query: QueryGraph,
    fields: &str,
    span: Range<usize>,
) -> Result<QueryGraph, SmartsParseError> {
    if fields.trim().is_empty() {
        return Ok(query);
    }
    use crate::core::StereoGroupKind as K;
    let fail = |message| SmartsParseError::syntax(span.clone(), message);
    let mut groups = BTreeMap::<(u8, usize), (K, Vec<QueryAtomId>)>::new();
    let mut current = None;
    let mut relative = false;
    for part in fields.split(',').map(str::trim) {
        if part == "r" {
            if relative {
                return Err(fail("duplicate CXSMARTS relative flag"));
            }
            relative = true;
            current = None;
            continue;
        }
        let atom = if let Some((header, atom)) = part.split_once(':') {
            let (tag, number, kind) = if header == "a" {
                (b'a', 0, K::Absolute)
            } else if header.starts_with('o') || header.starts_with('&') {
                if header.len() == 1 || !header[1..].bytes().all(|b| b.is_ascii_digit()) {
                    return Err(fail("invalid CXSMARTS stereo group number"));
                }
                let number = header[1..]
                    .parse::<usize>()
                    .map_err(|_| fail("invalid CXSMARTS stereo group number"))?;
                let tag = header.as_bytes()[0];
                (tag, number, if tag == b'o' { K::Or } else { K::And })
            } else {
                return Err(SmartsParseError::unsupported(
                    span.clone(),
                    "unsupported CXSMARTS field",
                ));
            };
            current = Some((tag, number));
            groups
                .entry((tag, number))
                .or_insert_with(|| (kind, Vec::new()));
            atom
        } else {
            part
        };
        let key = current.ok_or_else(|| fail("expected CXSMARTS stereo group field"))?;
        if atom.is_empty() || !atom.bytes().all(|b| b.is_ascii_digit()) {
            return Err(fail("invalid CXSMARTS member index"));
        }
        let index = atom
            .parse::<u32>()
            .map_err(|_| fail("CXSMARTS member index out of range"))?;
        groups
            .get_mut(&key)
            .unwrap()
            .1
            .push(QueryAtomId::new(index));
    }
    if relative {
        let grouped: std::collections::BTreeSet<_> = groups
            .values()
            .flat_map(|(_, m)| m.iter().copied())
            .collect();
        let members = query
            .atom_ids()
            .filter(|id| !grouped.contains(id) && query.atom(*id).unwrap().stereo_frame().is_some())
            .collect::<Vec<_>>();
        if !members.is_empty() {
            groups.insert((b'r', 0), (K::Relative, members));
        }
    }
    let mut builder = query.to_builder();
    for (_, (kind, members)) in groups {
        builder
            .add_stereo_group(super::QueryStereoGroup { kind, members })
            .map_err(|e| SmartsParseError::syntax(span.clone(), e.to_string()))?;
    }
    builder
        .build()
        .map_err(|e| SmartsParseError::syntax(span.clone(), e.to_string()))
}

fn validate_options(options: SmartsParseOptions) -> Result<(), SmartsParseError> {
    for (name, value) in [
        ("max_input_bytes", options.max_input_bytes),
        ("max_atoms", options.max_atoms),
        ("max_bonds", options.max_bonds),
        ("max_branch_depth", options.max_branch_depth),
        ("max_ring_closures", options.max_ring_closures),
        ("max_expression_nodes", options.max_expression_nodes),
        ("max_expression_depth", options.max_expression_depth),
        ("max_recursive_depth", options.max_recursive_depth),
        (
            "max_total_expression_nodes",
            options.max_total_expression_nodes,
        ),
    ] {
        if value == 0 {
            return Err(SmartsParseError::new(
                SmartsParseErrorKind::ResourceLimit,
                0..0,
                format!("SMARTS option {name} must be greater than zero"),
            ));
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PreviousToken {
    Start,
    Atom,
    Bond,
    Ring,
    BranchOpen,
    BranchClose,
    Dot,
}

impl PreviousToken {
    fn can_end_atom(self) -> bool {
        matches!(self, Self::Atom | Self::Ring | Self::BranchClose)
    }
}

struct RingOpen {
    atom: QueryAtomId,
    bond: Option<BondSyntax>,
    span: Range<usize>,
}

#[derive(Default)]
struct ParseBudget {
    atoms: Cell<usize>,
    bonds: Cell<usize>,
    nodes: Cell<usize>,
}

fn charge(
    counter: &Cell<usize>,
    amount: usize,
    limit: usize,
    name: &'static str,
    span: Range<usize>,
) -> Result<(), SmartsParseError> {
    let observed = counter.get().saturating_add(amount);
    if observed > limit {
        return Err(SmartsParseError::limit(span, name, observed, limit));
    }
    counter.set(observed);
    Ok(())
}

struct Parser<'a> {
    input: &'a str,
    bytes: &'a [u8],
    options: SmartsParseOptions,
    cursor: usize,
    builder: QueryGraphBuilder,
    budget: Rc<ParseBudget>,
    recursive_depth: usize,
    current: Option<QueryAtomId>,
    pending_bond: Option<BondSyntax>,
    stereo: StereoSyntax,
    branches: Vec<(QueryAtomId, usize)>,
    rings: BTreeMap<u16, RingOpen>,
    ring_closure_count: usize,
    previous: PreviousToken,
}

impl<'a> Parser<'a> {
    fn new(
        input: &'a str,
        options: SmartsParseOptions,
        budget: Rc<ParseBudget>,
        recursive_depth: usize,
    ) -> Self {
        Self {
            input,
            bytes: input.as_bytes(),
            options,
            cursor: 0,
            builder: QueryGraphBuilder::new(),
            budget,
            recursive_depth,
            current: None,
            pending_bond: None,
            stereo: StereoSyntax::default(),
            branches: Vec::new(),
            rings: BTreeMap::new(),
            ring_closure_count: 0,
            previous: PreviousToken::Start,
        }
    }

    fn parse(mut self) -> Result<QueryGraph, SmartsParseError> {
        while self.cursor < self.bytes.len() {
            match self.bytes[self.cursor] {
                b'(' => self.open_branch()?,
                b')' => self.close_branch()?,
                b'.' => self.component_separator()?,
                byte if is_bond_start(byte) => self.read_bond()?,
                b'0'..=b'9' | b'%' => self.read_ring()?,
                _ => self.read_atom()?,
            }
        }

        if matches!(
            self.previous,
            PreviousToken::Start
                | PreviousToken::Bond
                | PreviousToken::BranchOpen
                | PreviousToken::Dot
        ) {
            return Err(SmartsParseError::syntax(
                self.input.len()..self.input.len(),
                "incomplete SMARTS query",
            ));
        }
        if let Some((_, offset)) = self.branches.last() {
            return Err(SmartsParseError::syntax(
                *offset..(*offset + 1),
                "unclosed branch",
            ));
        }
        if let Some(open) = self.rings.values().next() {
            return Err(SmartsParseError::syntax(
                open.span.clone(),
                "unclosed ring label",
            ));
        }
        self.stereo.install(&mut self.builder)?;
        self.builder.build().map_err(|error| {
            let kind = if matches!(error, QueryGraphError::ResourceLimit { .. }) {
                SmartsParseErrorKind::ResourceLimit
            } else {
                SmartsParseErrorKind::InvalidSyntax
            };
            SmartsParseError::new(
                kind,
                0..self.input.len(),
                format!("invalid query graph: {error}"),
            )
        })
    }

    fn open_branch(&mut self) -> Result<(), SmartsParseError> {
        let offset = self.cursor;
        if matches!(self.previous, PreviousToken::Start | PreviousToken::Dot) {
            return Err(SmartsParseError::unsupported(
                offset..offset + 1,
                "component-level SMARTS grouping is not supported",
            ));
        }
        if !self.previous.can_end_atom() || self.pending_bond.is_some() {
            return Err(SmartsParseError::syntax(
                offset..offset + 1,
                "branch must follow an atom",
            ));
        }
        let parent = self
            .current
            .expect("an atom-ending token has a current atom");
        let observed = self.branches.len().saturating_add(1);
        if observed > self.options.max_branch_depth {
            return Err(SmartsParseError::limit(
                offset..offset + 1,
                "branch depth",
                observed,
                self.options.max_branch_depth,
            ));
        }
        self.branches.push((parent, offset));
        self.previous = PreviousToken::BranchOpen;
        self.cursor += 1;
        Ok(())
    }

    fn close_branch(&mut self) -> Result<(), SmartsParseError> {
        let offset = self.cursor;
        if !self.previous.can_end_atom() || self.pending_bond.is_some() {
            return Err(SmartsParseError::syntax(
                offset..offset + 1,
                "empty or incomplete branch",
            ));
        }
        let (parent, _) = self.branches.pop().ok_or_else(|| {
            SmartsParseError::syntax(offset..offset + 1, "unmatched branch close")
        })?;
        self.current = Some(parent);
        self.previous = PreviousToken::BranchClose;
        self.cursor += 1;
        Ok(())
    }

    fn component_separator(&mut self) -> Result<(), SmartsParseError> {
        let offset = self.cursor;
        if !self.previous.can_end_atom()
            || self.pending_bond.is_some()
            || !self.branches.is_empty()
            || !self.rings.is_empty()
        {
            return Err(SmartsParseError::syntax(
                offset..offset + 1,
                "component separator must follow a complete top-level component",
            ));
        }
        self.current = None;
        self.previous = PreviousToken::Dot;
        self.cursor += 1;
        Ok(())
    }

    fn read_bond(&mut self) -> Result<(), SmartsParseError> {
        let start = self.cursor;
        if !(self.previous.can_end_atom() || self.previous == PreviousToken::BranchOpen)
            || self.pending_bond.is_some()
        {
            return Err(SmartsParseError::syntax(
                start..start + 1,
                "bond expression must follow an atom or ring label",
            ));
        }
        let (expression, direction, end) = bond::parse(self.input, self.cursor, self.options)?;
        self.cursor = end;
        charge(
            &self.budget.nodes,
            expression.node_count(),
            self.options.max_total_expression_nodes,
            "total expression nodes",
            start..end,
        )?;
        self.pending_bond = Some(BondSyntax {
            expression,
            direction,
            span: start..self.cursor,
        });
        self.previous = PreviousToken::Bond;
        Ok(())
    }

    fn read_ring(&mut self) -> Result<(), SmartsParseError> {
        let start = self.cursor;
        if !self.previous.can_end_atom() && self.previous != PreviousToken::Bond {
            return Err(SmartsParseError::syntax(
                start..start + 1,
                "ring label must follow an atom",
            ));
        }
        let current = self
            .current
            .expect("an atom-ending token has a current atom");
        let label = self.parse_ring_label()?;
        let span = start..self.cursor;
        if let Some(open) = self.rings.remove(&label) {
            if open.atom == current {
                return Err(SmartsParseError::syntax(
                    span,
                    "ring label cannot create a self-bond",
                ));
            }
            let closing = self.pending_bond.take();
            let bond = BondSyntax::ring(open.bond, closing, span.clone())?;
            self.add_bond(open.atom, current, bond, span.clone())?;
            self.stereo.close_ring(open.atom, current, label);
            self.ring_closure_count = self.ring_closure_count.saturating_add(1);
            if self.ring_closure_count > self.options.max_ring_closures {
                return Err(SmartsParseError::limit(
                    span,
                    "ring closures",
                    self.ring_closure_count,
                    self.options.max_ring_closures,
                ));
            }
        } else {
            if self.rings.len() >= self.options.max_ring_closures {
                return Err(SmartsParseError::limit(
                    span,
                    "open ring labels",
                    self.rings.len().saturating_add(1),
                    self.options.max_ring_closures,
                ));
            }
            self.stereo.open_ring(current, label);
            self.rings.insert(
                label,
                RingOpen {
                    atom: current,
                    bond: self.pending_bond.take(),
                    span,
                },
            );
        }
        self.previous = PreviousToken::Ring;
        Ok(())
    }

    fn parse_ring_label(&mut self) -> Result<u16, SmartsParseError> {
        let start = self.cursor;
        if self.bytes[self.cursor] == b'%' {
            if self.cursor + 2 >= self.bytes.len()
                || !self.bytes[self.cursor + 1].is_ascii_digit()
                || !self.bytes[self.cursor + 2].is_ascii_digit()
            {
                return Err(SmartsParseError::syntax(
                    start..(start + 1).min(self.bytes.len()),
                    "percent ring labels require exactly two digits",
                ));
            }
            let label = u16::from(self.bytes[self.cursor + 1] - b'0') * 10
                + u16::from(self.bytes[self.cursor + 2] - b'0');
            self.cursor += 3;
            Ok(label)
        } else {
            let label = u16::from(self.bytes[self.cursor] - b'0');
            self.cursor += 1;
            Ok(label)
        }
    }

    fn read_atom(&mut self) -> Result<(), SmartsParseError> {
        let start = self.cursor;
        let (expression, stereo, tag) = if self.bytes[self.cursor] == b'[' {
            let mut level = 1usize;
            let mut close = self.cursor + 1;
            while close < self.bytes.len() {
                match self.bytes[close] {
                    b'[' => level += 1,
                    b']' => {
                        level -= 1;
                        if level == 0 {
                            break;
                        }
                    }
                    _ => {}
                }
                close += 1;
            }
            if close == self.bytes.len() {
                return Err(SmartsParseError::syntax(
                    start..close,
                    "unclosed bracket atom",
                ));
            }
            if close == self.cursor + 1 {
                return Err(SmartsParseError::syntax(
                    start..close + 1,
                    "empty bracket atom expression",
                ));
            }
            let expression = BracketParser::new(
                &self.input[self.cursor + 1..close],
                self.cursor + 1,
                self.options,
                self.budget.clone(),
                self.recursive_depth,
            )
            .parse()?;
            self.cursor = close + 1;
            expression
        } else {
            (self.parse_simple_atom()?, None, None)
        };

        charge(
            &self.budget.atoms,
            1,
            self.options.max_atoms,
            "atoms",
            start..self.cursor,
        )?;
        charge(
            &self.budget.nodes,
            expression.node_count(),
            self.options.max_total_expression_nodes,
            "total expression nodes",
            start..self.cursor,
        )?;
        let atom = self.builder.add_atom(expression).map_err(|error| {
            SmartsParseError::syntax(start..self.cursor, format!("invalid query atom: {error}"))
        })?;
        self.builder
            .set_atom_tag(atom, tag)
            .expect("new query atom");
        self.stereo.add_atom(atom, self.current, stereo);
        if let Some(previous_atom) = self.current {
            let bond = self
                .pending_bond
                .take()
                .map(Ok)
                .unwrap_or_else(|| BondSyntax::default_at(start..self.cursor))?;
            self.add_bond(previous_atom, atom, bond, start..self.cursor)?;
        } else if self.pending_bond.is_some() {
            return Err(SmartsParseError::syntax(
                start..self.cursor,
                "bond expression has no preceding atom",
            ));
        }
        self.current = Some(atom);
        self.previous = PreviousToken::Atom;
        Ok(())
    }

    fn parse_simple_atom(&mut self) -> Result<AtomExpression, SmartsParseError> {
        let start = self.cursor;
        match self.bytes[self.cursor] {
            b'*' => {
                self.cursor += 1;
                Ok(AtomExpression::always())
            }
            b'A' => {
                self.cursor += 1;
                Ok(AtomExpression::predicate(AtomPredicate::Aromatic(false)))
            }
            b'a' => {
                self.cursor += 1;
                Ok(AtomExpression::predicate(AtomPredicate::Aromatic(true)))
            }
            byte if byte.is_ascii_uppercase() && byte != b'X' => {
                // Only Cl and Br are two-letter unbracketed elements. In
                // `Cn`, for example, n starts a separate aromatic atom.
                let length = if matches!(self.bytes.get(start..start + 2), Some(b"Cl" | b"Br")) {
                    2
                } else {
                    1
                };
                let end = start + length;
                let symbol = &self.input[start..end];
                if !matches!(
                    symbol,
                    "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
                ) {
                    return Err(SmartsParseError::unsupported(
                        start..end,
                        "elements outside the organic subset must be bracketed",
                    ));
                }
                let element = Element::from_symbol(symbol).expect("organic subset element");
                self.cursor = end;
                atom_element_expression(element, false, start..end)
            }
            b'b' | b'c' | b'n' | b'o' | b'p' | b's' => {
                let symbol = match self.bytes[self.cursor] {
                    b'b' => "B",
                    b'c' => "C",
                    b'n' => "N",
                    b'o' => "O",
                    b'p' => "P",
                    b's' => "S",
                    _ => unreachable!(),
                };
                self.cursor += 1;
                atom_element_expression(
                    Element::from_symbol(symbol).expect("known aromatic element"),
                    true,
                    start..self.cursor,
                )
            }
            b'/' | b'\\' | b'@' => {
                self.cursor += 1;
                Err(SmartsParseError::unsupported(
                    start..self.cursor,
                    "stereochemical SMARTS is not supported",
                ))
            }
            _ => {
                self.cursor += 1;
                Err(SmartsParseError::syntax(
                    start..self.cursor,
                    format!(
                        "unexpected SMARTS character `{}`",
                        self.bytes[start] as char
                    ),
                ))
            }
        }
    }

    fn add_bond(
        &mut self,
        a: QueryAtomId,
        b: QueryAtomId,
        syntax: BondSyntax,
        span: Range<usize>,
    ) -> Result<(), SmartsParseError> {
        charge(
            &self.budget.bonds,
            1,
            self.options.max_bonds,
            "bonds",
            span.clone(),
        )?;
        let is_double = syntax.expression
            == non_aromatic_bond_expression(BondOrder::Double)
                .map_err(|error| expression_error(span.clone(), error))?;
        let id = self
            .builder
            .add_bond(a, b, syntax.expression)
            .map_err(|error| {
                SmartsParseError::syntax(span, format!("invalid query bond: {error}"))
            })?;
        self.stereo
            .add_bond(id, a, b, is_double, syntax.direction, syntax.span);
        Ok(())
    }
}

fn is_bond_start(byte: u8) -> bool {
    matches!(
        byte,
        b'-' | b'=' | b'#' | b'$' | b':' | b'~' | b'@' | b'!' | b'/' | b'\\'
    )
}

fn default_bond_expression() -> Result<BondExpression, SmartsParseError> {
    BondExpression::any([
        non_aromatic_bond_expression(BondOrder::Single)
            .map_err(|error| expression_error(0..0, error))?,
        aromatic_bond_expression().map_err(|error| expression_error(0..0, error))?,
    ])
    .map_err(|error| expression_error(0..0, error))
}

fn non_aromatic_bond_expression(order: BondOrder) -> Result<BondExpression, QueryExpressionError> {
    BondExpression::all([
        BondExpression::predicate(BondPredicate::Order(order)),
        BondExpression::predicate(BondPredicate::Aromatic(false)),
    ])
}

fn aromatic_bond_expression() -> Result<BondExpression, QueryExpressionError> {
    // SMARTS aromatic bond type is distinct from a localized single/double
    // order. Higher orders retain their type even when aromaticity flags them.
    BondExpression::all([
        BondExpression::predicate(BondPredicate::Aromatic(true)),
        BondExpression::any([
            BondExpression::predicate(BondPredicate::Order(BondOrder::Single)),
            BondExpression::predicate(BondPredicate::Order(BondOrder::Double)),
        ])?,
    ])
}

fn atom_element_expression(
    element: Element,
    aromatic: bool,
    span: Range<usize>,
) -> Result<AtomExpression, SmartsParseError> {
    AtomExpression::all([
        AtomExpression::predicate(AtomPredicate::Element(element)),
        AtomExpression::predicate(AtomPredicate::Aromatic(aromatic)),
    ])
    .map_err(|error| expression_error(span, error))
}

fn expression_error(span: Range<usize>, error: QueryExpressionError) -> SmartsParseError {
    SmartsParseError::new(SmartsParseErrorKind::ResourceLimit, span, error.to_string())
}

fn parse_element_symbol(input: &str, start: usize) -> Option<(Element, usize)> {
    let bytes = input.as_bytes();
    if !bytes.get(start)?.is_ascii_uppercase() {
        return None;
    }
    if bytes.get(start + 1).is_some_and(u8::is_ascii_lowercase) {
        let candidate = &input[start..start + 2];
        if let Some(element) = Element::from_symbol(candidate) {
            return Some((element, start + 2));
        }
    }
    let candidate = &input[start..start + 1];
    Element::from_symbol(candidate).map(|element| (element, start + 1))
}

struct BracketParser<'a> {
    source: &'a str,
    bytes: &'a [u8],
    base: usize,
    cursor: usize,
    options: SmartsParseOptions,
    elemental_hydrogen: bool,
    stereo: Option<AtomStereo>,
    tag: Option<u32>,
    budget: Rc<ParseBudget>,
    recursive_depth: usize,
}

impl<'a> BracketParser<'a> {
    fn new(
        source: &'a str,
        base: usize,
        options: SmartsParseOptions,
        budget: Rc<ParseBudget>,
        recursive_depth: usize,
    ) -> Self {
        Self {
            source,
            bytes: source.as_bytes(),
            base,
            cursor: 0,
            options,
            elemental_hydrogen: is_elemental_hydrogen_expression(
                source.split(':').next().unwrap_or(source),
            ),
            stereo: None,
            tag: None,
            budget,
            recursive_depth,
        }
    }

    fn parse(
        mut self,
    ) -> Result<(AtomExpression, Option<AtomStereo>, Option<u32>), SmartsParseError> {
        let expression = self.parse_low_and()?;
        if self.peek() == Some(b':') {
            self.cursor += 1;
            self.tag = Some(
                u32::try_from(self.read_number("atom tag")?)
                    .map_err(|_| self.syntax_here("atom tag exceeds u32"))?,
            );
        }
        if self.cursor != self.bytes.len() {
            return Err(self.syntax_here("unexpected atom-expression syntax"));
        }
        if expression.node_count() > self.options.max_expression_nodes {
            return Err(SmartsParseError::limit(
                self.base..self.base + self.bytes.len(),
                "expression nodes",
                expression.node_count(),
                self.options.max_expression_nodes,
            ));
        }
        if expression.depth() > self.options.max_expression_depth {
            return Err(SmartsParseError::limit(
                self.base..self.base + self.bytes.len(),
                "expression depth",
                expression.depth(),
                self.options.max_expression_depth,
            ));
        }
        Ok((expression, self.stereo, self.tag))
    }

    fn parse_low_and(&mut self) -> Result<AtomExpression, SmartsParseError> {
        let mut terms = vec![self.parse_or()?];
        while self.peek() == Some(b';') {
            self.cursor += 1;
            if self.at_expression_end() {
                return Err(self.syntax_here("missing expression after `;`"));
            }
            terms.push(self.parse_or()?);
        }
        self.compose_all(terms)
    }

    fn parse_or(&mut self) -> Result<AtomExpression, SmartsParseError> {
        let mut terms = vec![self.parse_high_and()?];
        while self.peek() == Some(b',') {
            self.cursor += 1;
            if self.at_expression_end() {
                return Err(self.syntax_here("missing expression after `,`"));
            }
            terms.push(self.parse_high_and()?);
        }
        self.compose_any(terms)
    }

    fn parse_high_and(&mut self) -> Result<AtomExpression, SmartsParseError> {
        let mut terms = vec![self.parse_unary()?];
        loop {
            if self.peek() == Some(b'&') {
                self.cursor += 1;
                if self.at_expression_end() {
                    return Err(self.syntax_here("missing expression after `&`"));
                }
                terms.push(self.parse_unary()?);
            } else if self.peek().is_some_and(is_primitive_start) {
                terms.push(self.parse_unary()?);
            } else {
                break;
            }
        }
        self.compose_all(terms)
    }

    fn parse_unary(&mut self) -> Result<AtomExpression, SmartsParseError> {
        let start = self.cursor;
        let mut negations = 0usize;
        while self.peek() == Some(b'!') {
            negations = negations.saturating_add(1);
            self.cursor += 1;
        }
        if self.at_expression_end() {
            return Err(self.syntax_here("negation must precede an atom primitive"));
        }
        let mut expression = self.parse_primitive()?;
        for _ in 0..negations {
            expression = expression
                .negate()
                .map_err(|error| expression_error(self.absolute(start..self.cursor), error))?;
        }
        Ok(expression)
    }

    fn parse_primitive(&mut self) -> Result<AtomExpression, SmartsParseError> {
        let start = self.cursor;
        let byte = self
            .peek()
            .ok_or_else(|| self.syntax_here("missing atom primitive"))?;
        if byte.is_ascii_uppercase()
            && self
                .bytes
                .get(self.cursor + 1)
                .is_some_and(u8::is_ascii_lowercase)
        {
            let absolute_start = self.base + self.cursor;
            if let Some((element, absolute_end)) =
                parse_element_symbol_at_bytes(self.source, self.base, self.cursor)
            {
                self.cursor = absolute_end - self.base;
                return atom_element_expression(element, false, absolute_start..absolute_end);
            }
        }
        if let Some(symbol) = match self.bytes.get(self.cursor..self.cursor + 2) {
            Some(b"as") => Some("As"),
            Some(b"se") => Some("Se"),
            _ => None,
        } {
            self.cursor += 2;
            return atom_element_expression(
                Element::from_symbol(symbol).expect("known aromatic element"),
                true,
                self.absolute(start..self.cursor),
            );
        }
        match byte {
            b'*' => {
                self.cursor += 1;
                Ok(AtomExpression::always())
            }
            b'#' => {
                self.cursor += 1;
                let atomic_number = self.read_number("atomic number")?;
                let atomic_number = u8::try_from(atomic_number)
                    .ok()
                    .and_then(Element::from_atomic_number)
                    .ok_or_else(|| {
                        SmartsParseError::syntax(
                            self.absolute(start..self.cursor),
                            "atomic number must be in 1..=118",
                        )
                    })?;
                Ok(AtomExpression::predicate(AtomPredicate::Element(
                    atomic_number,
                )))
            }
            b'0'..=b'9' => {
                let isotope = self.read_number("isotope")?;
                let isotope = u16::try_from(isotope).map_err(|_| {
                    SmartsParseError::syntax(
                        self.absolute(start..self.cursor),
                        "isotope exceeds u16 range",
                    )
                })?;
                Ok(AtomExpression::predicate(AtomPredicate::Isotope(isotope)))
            }
            b'A' => {
                self.cursor += 1;
                Ok(AtomExpression::predicate(AtomPredicate::Aromatic(false)))
            }
            b'a' => {
                self.cursor += 1;
                Ok(AtomExpression::predicate(AtomPredicate::Aromatic(true)))
            }
            b'D' => {
                self.cursor += 1;
                let degree = self.read_optional_u8(1, "degree")?;
                Ok(AtomExpression::predicate(AtomPredicate::Degree(degree)))
            }
            b'X' => {
                self.cursor += 1;
                let count = self.read_optional_u8(1, "total connectivity")?;
                Ok(AtomExpression::predicate(AtomPredicate::TotalConnectivity(
                    count,
                )))
            }
            b'x' => {
                self.cursor += 1;
                if self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    let count = self.read_optional_u8(0, "ring-bond count")?;
                    Ok(AtomExpression::predicate(AtomPredicate::RingBondCount(
                        count,
                    )))
                } else {
                    Ok(AtomExpression::predicate(AtomPredicate::RingMembership(
                        true,
                    )))
                }
            }
            b'H' if self.elemental_hydrogen => {
                self.cursor += 1;
                Ok(AtomExpression::predicate(AtomPredicate::Element(
                    Element::from_symbol("H").expect("hydrogen element"),
                )))
            }
            b'H' => {
                self.cursor += 1;
                let hydrogens = self.read_optional_u8(1, "hydrogen count")?;
                Ok(AtomExpression::predicate(AtomPredicate::TotalHydrogens(
                    hydrogens,
                )))
            }
            b'h' => {
                self.cursor += 1;
                let count = if self.peek().is_some_and(|b| b.is_ascii_digit()) {
                    Some(self.read_optional_u8(0, "implicit hydrogens")?)
                } else {
                    None
                };
                Ok(AtomExpression::predicate(AtomPredicate::ImplicitHydrogens(
                    count,
                )))
            }
            b'v' => {
                self.cursor += 1;
                Ok(AtomExpression::predicate(AtomPredicate::TotalValence(
                    self.read_optional_u8(1, "valence")?,
                )))
            }
            b'R' | b'r' => {
                self.cursor += 1;
                if !self.peek().is_some_and(|b| b.is_ascii_digit()) {
                    return Ok(AtomExpression::predicate(AtomPredicate::RingMembership(
                        true,
                    )));
                }
                let count = self.read_optional_u8(0, "ring predicate")?;
                Ok(AtomExpression::predicate(if count == 0 {
                    AtomPredicate::RingMembership(false)
                } else if byte == b'R' {
                    AtomPredicate::RingCount(count)
                } else {
                    AtomPredicate::SmallestRingSize(count)
                }))
            }
            b'+' | b'-' => self.parse_charge(),
            b'b' | b'c' | b'n' | b'o' | b'p' | b's' => {
                let symbol = match byte {
                    b'b' => "B",
                    b'c' => "C",
                    b'n' => "N",
                    b'o' => "O",
                    b'p' => "P",
                    b's' => "S",
                    _ => unreachable!(),
                };
                self.cursor += 1;
                atom_element_expression(
                    Element::from_symbol(symbol).expect("known aromatic element"),
                    true,
                    self.absolute(start..self.cursor),
                )
            }
            byte if byte.is_ascii_uppercase() && byte != b'X' => {
                let absolute_start = self.base + self.cursor;
                let (element, absolute_end) =
                    parse_element_symbol_at_bytes(self.source, self.base, self.cursor).ok_or_else(
                        || {
                            SmartsParseError::syntax(
                                self.absolute(start..start + 1),
                                "unknown element symbol",
                            )
                        },
                    )?;
                self.cursor = absolute_end - self.base;
                atom_element_expression(element, false, absolute_start..absolute_end)
            }
            b'$' => {
                self.cursor += 1;
                if self.peek() != Some(b'(') {
                    return Err(self.syntax_here("recursive query requires $(...)"));
                }
                let depth = self.recursive_depth + 1;
                let limit = self.options.max_recursive_depth.min(32);
                if depth > limit {
                    return Err(SmartsParseError::limit(
                        self.absolute(start..self.cursor),
                        "recursive depth",
                        depth,
                        limit,
                    ));
                }
                self.cursor += 1;
                let begin = self.cursor;
                let mut nesting = 1usize;
                while let Some(b) = self.peek() {
                    match b {
                        b'(' => nesting += 1,
                        b')' => {
                            nesting -= 1;
                            if nesting == 0 {
                                break;
                            }
                        }
                        _ => {}
                    }
                    self.cursor += 1;
                }
                if self.peek() != Some(b')') {
                    return Err(self.syntax_here("unclosed recursive query"));
                }
                let query = Parser::new(
                    &self.source[begin..self.cursor],
                    self.options,
                    self.budget.clone(),
                    depth,
                )
                .parse()
                .map_err(|mut e| {
                    e.span = self.absolute(begin + e.span.start..begin + e.span.end);
                    e
                })?;
                self.cursor += 1;
                Ok(AtomExpression::predicate(AtomPredicate::Recursive(
                    Box::new(query),
                )))
            }
            b'@' => {
                let (mut stereo, end) = AtomStereo::parse(self.source, self.cursor, self.base)?;
                let same = stereo.orientation == crate::core::TetrahedralOrientation::Clockwise;
                stereo.orientation = crate::core::TetrahedralOrientation::Clockwise;
                if let Some(existing) = &self.stereo {
                    if existing.inline_hydrogen != stereo.inline_hydrogen {
                        return Err(self.syntax_here("inconsistent inline hydrogen stereo frames"));
                    }
                } else {
                    self.stereo = Some(stereo);
                }
                self.cursor = end;
                Ok(AtomExpression::predicate(AtomPredicate::Tetrahedral(same)))
            }
            b'^' => {
                self.cursor += 1;
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.cursor += 1;
                }
                Err(SmartsParseError::unsupported(
                    self.absolute(start..self.cursor),
                    "valence and hybridization primitives are outside the bounded SMARTS subset",
                ))
            }
            _ => {
                self.cursor += 1;
                Err(SmartsParseError::syntax(
                    self.absolute(start..self.cursor),
                    format!("unexpected atom primitive `{}`", byte as char),
                ))
            }
        }
    }

    fn parse_charge(&mut self) -> Result<AtomExpression, SmartsParseError> {
        let start = self.cursor;
        let sign = self.bytes[self.cursor];
        self.cursor += 1;
        let magnitude = if self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
            self.read_number("formal charge")?
        } else {
            let mut magnitude = 1usize;
            while self.peek() == Some(sign) {
                magnitude = magnitude.saturating_add(1);
                self.cursor += 1;
            }
            magnitude
        };
        let magnitude = i8::try_from(magnitude).map_err(|_| {
            SmartsParseError::syntax(
                self.absolute(start..self.cursor),
                "formal charge exceeds i8 range",
            )
        })?;
        let charge = if sign == b'+' { magnitude } else { -magnitude };
        Ok(AtomExpression::predicate(AtomPredicate::FormalCharge(
            charge,
        )))
    }

    fn read_optional_u8(
        &mut self,
        default: u8,
        description: &'static str,
    ) -> Result<u8, SmartsParseError> {
        let start = self.cursor;
        match self.read_optional_number() {
            None => Ok(default),
            Some(value) => u8::try_from(value).map_err(|_| {
                SmartsParseError::syntax(
                    self.absolute(start..self.cursor),
                    format!("{description} exceeds u8 range"),
                )
            }),
        }
    }

    fn read_number(&mut self, description: &'static str) -> Result<usize, SmartsParseError> {
        let start = self.cursor;
        self.read_optional_number().ok_or_else(|| {
            SmartsParseError::syntax(
                self.absolute(start..self.cursor),
                format!("missing {description}"),
            )
        })
    }

    fn read_optional_number(&mut self) -> Option<usize> {
        let start = self.cursor;
        let mut value = 0usize;
        while let Some(byte) = self.peek().filter(u8::is_ascii_digit) {
            value = value
                .saturating_mul(10)
                .saturating_add(usize::from(byte - b'0'));
            self.cursor += 1;
        }
        (self.cursor > start).then_some(value)
    }

    fn compose_all(
        &self,
        expressions: Vec<AtomExpression>,
    ) -> Result<AtomExpression, SmartsParseError> {
        AtomExpression::all(expressions)
            .map_err(|error| expression_error(self.base..self.base + self.bytes.len(), error))
    }

    fn compose_any(
        &self,
        expressions: Vec<AtomExpression>,
    ) -> Result<AtomExpression, SmartsParseError> {
        AtomExpression::any(expressions)
            .map_err(|error| expression_error(self.base..self.base + self.bytes.len(), error))
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.cursor).copied()
    }

    fn at_expression_end(&self) -> bool {
        self.cursor >= self.bytes.len() || matches!(self.peek(), Some(b',' | b';' | b':'))
    }

    fn absolute(&self, span: Range<usize>) -> Range<usize> {
        self.base + span.start..self.base + span.end
    }

    fn syntax_here(&self, message: impl Into<String>) -> SmartsParseError {
        let end = (self.cursor + 1).min(self.bytes.len());
        SmartsParseError::syntax(self.absolute(self.cursor..end), message)
    }
}

fn is_primitive_start(byte: u8) -> bool {
    !matches!(byte, b',' | b';' | b'&' | b':')
}

fn parse_element_symbol_at_bytes(
    source: &str,
    base: usize,
    cursor: usize,
) -> Option<(Element, usize)> {
    let local = &source[cursor..];
    let (element, end) = parse_element_symbol(local, 0)?;
    Some((element, base + cursor + end))
}

fn is_elemental_hydrogen_expression(source: &str) -> bool {
    let bytes = source.as_bytes();
    let mut cursor = 0usize;
    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        cursor += 1;
    }
    if bytes.get(cursor) != Some(&b'H') {
        return false;
    }
    cursor += 1;
    if cursor == bytes.len() {
        return true;
    }
    let Some(sign @ (b'+' | b'-')) = bytes.get(cursor).copied() else {
        return false;
    };
    cursor += 1;
    if bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
        while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
            cursor += 1;
        }
    } else {
        while bytes.get(cursor) == Some(&sign) {
            cursor += 1;
        }
    }
    cursor == bytes.len()
}
