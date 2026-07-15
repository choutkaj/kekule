use std::collections::BTreeMap;
use std::fmt;
use std::ops::Range;

use crate::core::{BondOrder, Element};

use super::{
    AtomExpression, AtomPredicate, BondExpression, BondPredicate, QueryAtomId,
    QueryExpressionError, QueryGraph, QueryGraphBuilder,
};

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
        }
    }
}

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

pub fn parse_smarts(input: &str) -> Result<QueryGraph, SmartsParseError> {
    parse_smarts_with_options(input, SmartsParseOptions::default())
}

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
    Parser::new(input, options).parse()
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
    bond: Option<(BondExpression, Range<usize>)>,
    span: Range<usize>,
}

struct Parser<'a> {
    input: &'a str,
    bytes: &'a [u8],
    options: SmartsParseOptions,
    cursor: usize,
    builder: QueryGraphBuilder,
    atom_count: usize,
    bond_count: usize,
    current: Option<QueryAtomId>,
    pending_bond: Option<(BondExpression, Range<usize>)>,
    branches: Vec<(QueryAtomId, usize)>,
    rings: BTreeMap<u16, RingOpen>,
    ring_closure_count: usize,
    previous: PreviousToken,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str, options: SmartsParseOptions) -> Self {
        Self {
            input,
            bytes: input.as_bytes(),
            options,
            cursor: 0,
            builder: QueryGraphBuilder::new(),
            atom_count: 0,
            bond_count: 0,
            current: None,
            pending_bond: None,
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
        self.builder.build().map_err(|error| {
            SmartsParseError::syntax(0..self.input.len(), format!("invalid query graph: {error}"))
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
        let expression = match self.bytes[self.cursor] {
            b'-' => {
                self.cursor += 1;
                BondExpression::predicate(BondPredicate::Order(BondOrder::Single))
            }
            b'=' => {
                self.cursor += 1;
                BondExpression::predicate(BondPredicate::Order(BondOrder::Double))
            }
            b'#' => {
                self.cursor += 1;
                BondExpression::predicate(BondPredicate::Order(BondOrder::Triple))
            }
            b'$' => {
                self.cursor += 1;
                BondExpression::predicate(BondPredicate::Order(BondOrder::Quadruple))
            }
            b':' => {
                self.cursor += 1;
                BondExpression::predicate(BondPredicate::Aromatic(true))
            }
            b'~' => {
                self.cursor += 1;
                BondExpression::always()
            }
            b'@' => {
                self.cursor += 1;
                BondExpression::predicate(BondPredicate::RingMembership(true))
            }
            b'!' if self.bytes.get(self.cursor + 1) == Some(&b'@') => {
                self.cursor += 2;
                BondExpression::predicate(BondPredicate::RingMembership(false))
            }
            b'/' | b'\\' => {
                self.cursor += 1;
                return Err(SmartsParseError::unsupported(
                    start..self.cursor,
                    "directional and stereochemical bond queries are not supported",
                ));
            }
            _ => {
                self.cursor += 1;
                return Err(SmartsParseError::unsupported(
                    start..self.cursor,
                    "bounded SMARTS supports only -, =, #, $, :, ~, @, and !@ bond expressions",
                ));
            }
        };
        self.pending_bond = Some((expression, start..self.cursor));
        self.previous = PreviousToken::Bond;
        Ok(())
    }

    fn read_ring(&mut self) -> Result<(), SmartsParseError> {
        let start = self.cursor;
        if !self.previous.can_end_atom() {
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
            let expression = match (open.bond, closing) {
                (None, None) => default_bond_expression()?,
                (Some((expression, _)), None) | (None, Some((expression, _))) => expression,
                (Some((left, _)), Some((right, _))) if left == right => left,
                (Some((_, left_span)), Some((_, right_span))) => {
                    return Err(SmartsParseError::syntax(
                        left_span.start..right_span.end,
                        "conflicting bond expressions on ring closure",
                    ));
                }
            };
            self.add_bond(open.atom, current, expression, span.clone())?;
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
            self.rings.insert(
                label,
                RingOpen {
                    atom: current,
                    bond: self.pending_bond.take(),
                    span: span.clone(),
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
        let expression = if self.bytes[self.cursor] == b'[' {
            let close = self.bytes[self.cursor + 1..]
                .iter()
                .position(|byte| *byte == b']')
                .map(|relative| self.cursor + 1 + relative)
                .ok_or_else(|| {
                    SmartsParseError::syntax(start..self.input.len(), "unclosed bracket atom")
                })?;
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
            )
            .parse()?;
            self.cursor = close + 1;
            expression
        } else {
            self.parse_simple_atom()?
        };

        self.atom_count = self.atom_count.saturating_add(1);
        if self.atom_count > self.options.max_atoms {
            return Err(SmartsParseError::limit(
                start..self.cursor,
                "atoms",
                self.atom_count,
                self.options.max_atoms,
            ));
        }
        let atom = self.builder.add_atom(expression).map_err(|error| {
            SmartsParseError::syntax(start..self.cursor, format!("invalid query atom: {error}"))
        })?;
        if let Some(previous_atom) = self.current {
            let bond = self
                .pending_bond
                .take()
                .map(|(expression, _)| expression)
                .map(Ok)
                .unwrap_or_else(default_bond_expression)?;
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
                let (element, end) =
                    parse_element_symbol(self.input, self.cursor).ok_or_else(|| {
                        SmartsParseError::syntax(
                            start..start + 1,
                            "unbracketed atom is not in the organic subset",
                        )
                    })?;
                if !matches!(
                    element.symbol(),
                    "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
                ) {
                    return Err(SmartsParseError::unsupported(
                        start..end,
                        "elements outside the organic subset must be bracketed",
                    ));
                }
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
        expression: BondExpression,
        span: Range<usize>,
    ) -> Result<(), SmartsParseError> {
        self.bond_count = self.bond_count.saturating_add(1);
        if self.bond_count > self.options.max_bonds {
            return Err(SmartsParseError::limit(
                span,
                "bonds",
                self.bond_count,
                self.options.max_bonds,
            ));
        }
        self.builder.add_bond(a, b, expression).map_err(|error| {
            SmartsParseError::syntax(span, format!("invalid query bond: {error}"))
        })?;
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
        BondExpression::predicate(BondPredicate::Order(BondOrder::Single)),
        BondExpression::predicate(BondPredicate::Aromatic(true)),
    ])
    .map_err(|error| expression_error(0..0, error))
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
}

impl<'a> BracketParser<'a> {
    fn new(source: &'a str, base: usize, options: SmartsParseOptions) -> Self {
        Self {
            source,
            bytes: source.as_bytes(),
            base,
            cursor: 0,
            options,
            elemental_hydrogen: is_elemental_hydrogen_expression(source),
        }
    }

    fn parse(mut self) -> Result<AtomExpression, SmartsParseError> {
        let expression = self.parse_low_and()?;
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
        Ok(expression)
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
                if isotope == 0 {
                    return Err(SmartsParseError::syntax(
                        self.absolute(start..self.cursor),
                        "isotope must be greater than zero",
                    ));
                }
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
            b'R' => {
                self.cursor += 1;
                let value = self.read_optional_number();
                match value {
                    None => Ok(AtomExpression::predicate(AtomPredicate::RingMembership(
                        true,
                    ))),
                    Some(0) => Ok(AtomExpression::predicate(AtomPredicate::RingMembership(
                        false,
                    ))),
                    Some(_) => Err(SmartsParseError::unsupported(
                        self.absolute(start..self.cursor),
                        "exact ring-membership counts are outside the bounded SMARTS subset",
                    )),
                }
            }
            b'r' => {
                self.cursor += 1;
                if self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.read_optional_number();
                    return Err(SmartsParseError::unsupported(
                        self.absolute(start..self.cursor),
                        "ring-size predicates are outside the bounded SMARTS subset",
                    ));
                }
                Ok(AtomExpression::predicate(AtomPredicate::RingMembership(
                    true,
                )))
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
            b'$' | b'(' | b')' => {
                self.cursor += 1;
                Err(SmartsParseError::unsupported(
                    self.absolute(start..self.cursor),
                    "recursive SMARTS and atom-expression grouping are not supported",
                ))
            }
            b'@' => {
                self.cursor += 1;
                Err(SmartsParseError::unsupported(
                    self.absolute(start..self.cursor),
                    "stereochemical atom queries are not supported",
                ))
            }
            b':' => {
                self.cursor += 1;
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.cursor += 1;
                }
                Err(SmartsParseError::unsupported(
                    self.absolute(start..self.cursor),
                    "atom maps are not part of substructure predicate semantics",
                ))
            }
            b'X' | b'x' | b'v' | b'^' => {
                self.cursor += 1;
                while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                    self.cursor += 1;
                }
                Err(SmartsParseError::unsupported(
                    self.absolute(start..self.cursor),
                    "connectivity, ring-bond-count, valence, and hybridization primitives are outside the bounded SMARTS subset",
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
        self.cursor >= self.bytes.len() || matches!(self.peek(), Some(b',' | b';'))
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
    !matches!(byte, b',' | b';' | b'&')
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
