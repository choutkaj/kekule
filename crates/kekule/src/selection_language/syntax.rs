//! Lexing and recursive-descent parsing of selection text.

use std::ops::Range;

use super::ast::{
    CompareOp, DistanceUnit, Field, FieldType, Flag, HostTermKind, Level, Literal, Node, NodeKind,
    Span, ValueMatch,
};
use super::{SelectParseOptions, SelectionExpr};
use crate::diagnostics::{
    DiagnosticKind, Format, FormatError, LimitExceeded, Resource, SourceLocation,
};
use crate::molecule::Element;
use crate::properties::PropertyKey;
use crate::topology::{MoleculeClass, ResidueClass};
use crate::units::{self, Dimension, Unit};

/// Every reserved keyword; an unquoted value equal to one ends a value list.
pub(crate) const KEYWORDS: &[&str] = &[
    "and",
    "or",
    "not",
    "to",
    "of",
    "as",
    "within",
    "byres",
    "bychain",
    "bymolecule",
    "same",
    "residue",
    "chain",
    "molecule",
    "bonded",
    "all",
    "none",
    "smarts",
    "property",
    "aromatic",
    "ring",
    "hydrogen",
    "heavy",
    "polar_hydrogen",
    "protein",
    "nucleic",
    "water",
    "ion",
    "ligand",
    "backbone",
    "sidechain",
    "index",
    "element",
    "elem",
    "name",
    "label_name",
    "charge",
    "isotope",
    "bfactor",
    "occupancy",
    "resname",
    "resid",
    "label_resid",
    "resclass",
    "label_chain",
    "molclass",
    "instance",
    "definition",
];

#[derive(Debug, Clone, PartialEq)]
enum TokenKind {
    Open,
    Close,
    Op(CompareOp),
    Word { text: String, quoted: bool },
}

#[derive(Debug, Clone)]
struct Token {
    kind: TokenKind,
    span: Range<usize>,
}

pub(crate) fn error(
    input: &str,
    kind: DiagnosticKind,
    span: Range<usize>,
    message: impl Into<String>,
) -> FormatError {
    let start = span.start.min(input.len());
    let before = &input.as_bytes()[..start];
    let line = 1 + before.iter().filter(|&&byte| byte == b'\n').count();
    let line_start = before
        .iter()
        .rposition(|&byte| byte == b'\n')
        .map_or(0, |index| index + 1);
    FormatError::parse(Format::Selection, kind, message)
        .at(SourceLocation::at_column(line, start - line_start + 1))
        .with_span(span)
}

fn syntax(input: &str, span: Range<usize>, message: impl Into<String>) -> FormatError {
    error(input, DiagnosticKind::Syntax, span, message)
}

fn invalid(input: &str, span: Range<usize>, message: impl Into<String>) -> FormatError {
    error(input, DiagnosticKind::InvalidValue, span, message)
}

fn is_word_byte(character: char) -> bool {
    !character.is_whitespace() && !matches!(character, '(' | ')' | '"' | '<' | '>' | '=' | '!')
}

fn lex(input: &str) -> Result<Vec<Token>, FormatError> {
    let mut tokens = Vec::new();
    let mut chars = input.char_indices().peekable();
    while let Some(&(start, character)) = chars.peek() {
        if character.is_whitespace() {
            chars.next();
            continue;
        }
        match character {
            '(' | ')' => {
                chars.next();
                tokens.push(Token {
                    kind: if character == '(' {
                        TokenKind::Open
                    } else {
                        TokenKind::Close
                    },
                    span: start..start + 1,
                });
            }
            '<' | '>' | '=' | '!' => {
                chars.next();
                let equals = chars.peek().is_some_and(|&(_, next)| next == '=');
                if equals {
                    chars.next();
                }
                let op = match (character, equals) {
                    ('<', false) => CompareOp::Less,
                    ('<', true) => CompareOp::LessEqual,
                    ('>', false) => CompareOp::Greater,
                    ('>', true) => CompareOp::GreaterEqual,
                    ('=', _) => CompareOp::Equal,
                    ('!', true) => CompareOp::NotEqual,
                    _ => {
                        return Err(syntax(
                            input,
                            start..start + 1,
                            "`!` must be followed by `=`; write `not` for negation",
                        ))
                    }
                };
                tokens.push(Token {
                    kind: TokenKind::Op(op),
                    span: start..start + if equals { 2 } else { 1 },
                });
            }
            '"' | '\'' => {
                chars.next();
                let mut text = String::new();
                let mut end = None;
                while let Some((index, next)) = chars.next() {
                    if next == character {
                        end = Some(index + next.len_utf8());
                        break;
                    }
                    if next == '\\' {
                        match chars.next() {
                            Some((_, escaped)) => text.push(escaped),
                            None => break,
                        }
                    } else {
                        text.push(next);
                    }
                }
                let Some(end) = end else {
                    return Err(syntax(
                        input,
                        start..input.len(),
                        "unterminated quoted value",
                    ));
                };
                tokens.push(Token {
                    kind: TokenKind::Word { text, quoted: true },
                    span: start..end,
                });
            }
            _ => {
                let mut end = start;
                while let Some(&(index, next)) = chars.peek() {
                    if !is_word_byte(next) {
                        break;
                    }
                    end = index + next.len_utf8();
                    chars.next();
                }
                tokens.push(Token {
                    kind: TokenKind::Word {
                        text: input[start..end].to_owned(),
                        quoted: false,
                    },
                    span: start..end,
                });
            }
        }
    }
    Ok(tokens)
}

pub(crate) fn parse(input: &str, options: &SelectParseOptions) -> Result<Node, FormatError> {
    if input.len() > options.max_input_bytes() {
        return Err(error(
            input,
            DiagnosticKind::Limit,
            0..0,
            "selection text exceeds the configured byte limit",
        )
        .exceeding(LimitExceeded::new(
            Resource::InputBytes,
            options.max_input_bytes(),
            input.len(),
        )));
    }
    let tokens = lex(input)?;
    let mut parser = Parser {
        input,
        tokens,
        position: 0,
        depth: 0,
        options,
    };
    if parser.tokens.is_empty() {
        return Err(syntax(input, 0..input.len(), "empty selection"));
    }
    let node = parser.or()?;
    if let Some(token) = parser.tokens.get(parser.position) {
        let message = match &token.kind {
            TokenKind::Close => "unmatched `)`".to_owned(),
            _ => "expected `and`, `or`, or the end; terms are never joined implicitly".to_owned(),
        };
        return Err(syntax(input, token.span.clone(), message));
    }
    Ok(node)
}

struct Parser<'a> {
    input: &'a str,
    tokens: Vec<Token>,
    position: usize,
    depth: usize,
    options: &'a SelectParseOptions,
}

/// A recognized unit spelling and its canonical spelling.
fn unit(word: &str) -> Option<(Unit, &'static str)> {
    Some(match word {
        "A" | "Å" | "angstrom" => (units::ANGSTROM, "A"),
        "nm" => (units::NANOMETER, "nm"),
        "pm" => (
            Unit::new(Dimension::LENGTH, 1.0e-12, Some("pm")).ok()?,
            "pm",
        ),
        "bohr" => (units::BOHR, "bohr"),
        "Da" => (units::DALTON, "Da"),
        "e" => (units::ELEMENTARY_CHARGE, "e"),
        "K" => (units::KELVIN, "K"),
        "ps" => (units::PICOSECOND, "ps"),
        "fs" => (units::FEMTOSECOND, "fs"),
        "kJ/mol" => (units::KILOJOULE_PER_MOLE, "kJ/mol"),
        "kcal/mol" => (units::KILOCALORIE_PER_MOLE, "kcal/mol"),
        "deg" => (units::DEGREE, "deg"),
        "rad" => (units::RADIAN, "rad"),
        _ => return None,
    })
}

fn distance_unit(word: &str) -> Option<DistanceUnit> {
    Some(match word {
        "A" | "Å" | "angstrom" => DistanceUnit::Angstrom,
        "nm" => DistanceUnit::Nanometer,
        "pm" => DistanceUnit::Picometer,
        "bohr" => DistanceUnit::Bohr,
        _ => return None,
    })
}

fn flag(word: &str) -> Option<Flag> {
    Some(match word {
        "aromatic" => Flag::Aromatic,
        "ring" => Flag::Ring,
        "hydrogen" => Flag::Hydrogen,
        "heavy" => Flag::Heavy,
        "polar_hydrogen" => Flag::PolarHydrogen,
        "protein" => Flag::Protein,
        "nucleic" => Flag::Nucleic,
        "water" => Flag::Water,
        "ion" => Flag::Ion,
        "ligand" => Flag::Ligand,
        "backbone" => Flag::Backbone,
        "sidechain" => Flag::Sidechain,
        _ => return None,
    })
}

fn field(word: &str) -> Option<Field> {
    Some(match word {
        "index" => Field::Index,
        "element" | "elem" => Field::Element,
        "name" => Field::Name,
        "label_name" => Field::LabelName,
        "charge" => Field::Charge,
        "isotope" => Field::Isotope,
        "bfactor" => Field::BFactor,
        "occupancy" => Field::Occupancy,
        "resname" => Field::ResName,
        "resid" => Field::ResId,
        "label_resid" => Field::LabelResId,
        "resclass" => Field::ResClass,
        "chain" => Field::Chain,
        "label_chain" => Field::LabelChain,
        "molclass" => Field::MolClass,
        "instance" => Field::Instance,
        "definition" => Field::Definition,
        _ => return None,
    })
}

pub(crate) fn residue_class(word: &str) -> Option<ResidueClass> {
    Some(match word {
        "amino_acid" => ResidueClass::AminoAcid,
        "dna" => ResidueClass::DnaNucleotide,
        "rna" => ResidueClass::RnaNucleotide,
        "carbohydrate" => ResidueClass::Carbohydrate,
        "water" => ResidueClass::Water,
        "ion" => ResidueClass::Ion,
        "other" => ResidueClass::Other,
        _ => return None,
    })
}

pub(crate) fn molecule_class(word: &str) -> Option<MoleculeClass> {
    Some(match word {
        "protein" => MoleculeClass::Protein,
        "dna" => MoleculeClass::Dna,
        "rna" => MoleculeClass::Rna,
        "carbohydrate" => MoleculeClass::Carbohydrate,
        "water" => MoleculeClass::Water,
        "ion" => MoleculeClass::Ion,
        "small_molecule" => MoleculeClass::SmallMolecule,
        "other" => MoleculeClass::Other,
        _ => return None,
    })
}

fn element(word: &str) -> Option<Element> {
    let mut characters = word.chars();
    let first = characters.next()?.to_ascii_uppercase();
    let normalized = std::iter::once(first)
        .chain(characters.map(|character| character.to_ascii_lowercase()))
        .collect::<String>();
    Element::from_symbol(&normalized)
}

/// A residue number with an optional insertion code, such as `52` or `52A`.
fn residue_number(word: &str) -> Option<(i64, Option<String>)> {
    let digits = word
        .char_indices()
        .find(|&(index, character)| {
            !(character.is_ascii_digit() || (index == 0 && character == '-'))
        })
        .map_or(word.len(), |(index, _)| index);
    let number = word[..digits].parse().ok()?;
    let code = &word[digits..];
    if !code
        .chars()
        .all(|character| character.is_ascii_alphabetic())
    {
        return None;
    }
    Some((number, (!code.is_empty()).then(|| code.to_owned())))
}

impl Parser<'_> {
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.position)
    }

    fn end_span(&self) -> Range<usize> {
        self.input.len()..self.input.len()
    }

    fn keyword_at(&self, position: usize) -> Option<&str> {
        match self.tokens.get(position) {
            Some(Token {
                kind:
                    TokenKind::Word {
                        text,
                        quoted: false,
                    },
                ..
            }) => Some(text.as_str()),
            _ => None,
        }
    }

    fn peek_keyword(&self) -> Option<&str> {
        self.keyword_at(self.position)
    }

    fn expect_keyword(&mut self, keyword: &str) -> Result<Range<usize>, FormatError> {
        match self.tokens.get(self.position) {
            Some(token) if matches!(&token.kind, TokenKind::Word { text, quoted: false } if text == keyword) =>
            {
                self.position += 1;
                Ok(token.span.clone())
            }
            Some(token) => Err(syntax(
                self.input,
                token.span.clone(),
                format!("expected `{keyword}`"),
            )),
            None => Err(syntax(
                self.input,
                self.end_span(),
                format!("expected `{keyword}` before the end"),
            )),
        }
    }

    fn enter(&mut self, span: Range<usize>) -> Result<(), FormatError> {
        self.depth += 1;
        if self.depth > self.options.max_depth() {
            return Err(error(
                self.input,
                DiagnosticKind::Limit,
                span,
                "selection nesting exceeds the configured depth limit",
            )
            .exceeding(LimitExceeded::new(
                Resource::ExpressionDepth,
                self.options.max_depth(),
                self.depth,
            )));
        }
        Ok(())
    }

    fn or(&mut self) -> Result<Node, FormatError> {
        let first = self.and()?;
        if self.peek_keyword() != Some("or") {
            return Ok(first);
        }
        let start = first.span.0.start;
        let mut operands = Vec::new();
        push_flat(&mut operands, first, true);
        while self.peek_keyword() == Some("or") {
            self.position += 1;
            let operand = self.and()?;
            push_flat(&mut operands, operand, true);
        }
        let end = operands.last().map_or(start, |node| node.span.0.end);
        Ok(Node {
            kind: NodeKind::Or(operands),
            span: Span(start..end),
        })
    }

    fn and(&mut self) -> Result<Node, FormatError> {
        let first = self.unary()?;
        if self.peek_keyword() != Some("and") {
            return Ok(first);
        }
        let start = first.span.0.start;
        let mut operands = Vec::new();
        push_flat(&mut operands, first, false);
        while self.peek_keyword() == Some("and") {
            self.position += 1;
            let operand = self.unary()?;
            push_flat(&mut operands, operand, false);
        }
        let end = operands.last().map_or(start, |node| node.span.0.end);
        Ok(Node {
            kind: NodeKind::And(operands),
            span: Span(start..end),
        })
    }

    fn unary(&mut self) -> Result<Node, FormatError> {
        let Some(token) = self.peek().cloned() else {
            return Err(syntax(
                self.input,
                self.end_span(),
                "expected a term before the end",
            ));
        };
        self.enter(token.span.clone())?;
        let start = token.span.start;
        let node = match self.peek_keyword() {
            Some("not") => {
                self.position += 1;
                let inner = self.unary()?;
                let end = inner.span.0.end;
                Node {
                    kind: NodeKind::Not(Box::new(inner)),
                    span: Span(start..end),
                }
            }
            Some("within") => {
                self.position += 1;
                let (distance, unit) = self.distance()?;
                self.expect_keyword("of")?;
                let reference = self.unary()?;
                let end = reference.span.0.end;
                Node {
                    kind: NodeKind::Within {
                        distance,
                        unit,
                        reference: Box::new(SelectionExpr { root: reference }),
                    },
                    span: Span(start..end),
                }
            }
            Some(word @ ("byres" | "bychain" | "bymolecule")) => {
                let level = match word {
                    "byres" => Level::Residue,
                    "bychain" => Level::Chain,
                    _ => Level::Molecule,
                };
                self.position += 1;
                self.expand(start, level)?
            }
            Some("same") => {
                self.position += 1;
                let level = match self.peek_keyword() {
                    Some("residue") => Level::Residue,
                    Some("chain") => Level::Chain,
                    Some("molecule") => Level::Molecule,
                    _ => {
                        let span = self
                            .peek()
                            .map_or(self.end_span(), |token| token.span.clone());
                        return Err(syntax(
                            self.input,
                            span,
                            "expected `residue`, `chain`, or `molecule` after `same`",
                        ));
                    }
                };
                self.position += 1;
                self.expect_keyword("as")?;
                self.expand(start, level)?
            }
            Some("bonded") => {
                self.position += 1;
                self.expect_keyword("to")?;
                let inner = self.unary()?;
                let end = inner.span.0.end;
                Node {
                    kind: NodeKind::Bonded(Box::new(inner)),
                    span: Span(start..end),
                }
            }
            _ => self.primary()?,
        };
        self.depth -= 1;
        Ok(node)
    }

    fn expand(&mut self, start: usize, level: Level) -> Result<Node, FormatError> {
        let inner = self.unary()?;
        let end = inner.span.0.end;
        Ok(Node {
            kind: NodeKind::Expand {
                level,
                inner: Box::new(inner),
            },
            span: Span(start..end),
        })
    }

    fn distance(&mut self) -> Result<(f64, DistanceUnit), FormatError> {
        let Some(token) = self.peek().cloned() else {
            return Err(syntax(
                self.input,
                self.end_span(),
                "expected a distance after `within`",
            ));
        };
        let TokenKind::Word {
            text,
            quoted: false,
        } = &token.kind
        else {
            return Err(syntax(
                self.input,
                token.span,
                "expected a distance after `within`",
            ));
        };
        let distance = text
            .parse::<f64>()
            .ok()
            .filter(|distance| distance.is_finite() && *distance >= 0.0)
            .ok_or_else(|| {
                invalid(
                    self.input,
                    token.span.clone(),
                    "a distance is a finite non-negative number, separated from its unit",
                )
            })?;
        self.position += 1;
        let unit = match self.peek_keyword().and_then(distance_unit) {
            Some(unit) => {
                self.position += 1;
                unit
            }
            None => DistanceUnit::Angstrom,
        };
        Ok((distance, unit))
    }

    /// Whether the next token can be a value rather than the end of a list.
    fn at_value(&self) -> bool {
        match self.peek() {
            Some(Token {
                kind: TokenKind::Word { quoted: true, .. },
                ..
            }) => true,
            Some(Token {
                kind:
                    TokenKind::Word {
                        text,
                        quoted: false,
                    },
                ..
            }) => {
                !KEYWORDS.contains(&text.as_str())
                    && !text.starts_with(['#', '$'])
                    && !self.options.is_host_keyword(text)
            }
            _ => false,
        }
    }

    /// Whether the next word is a class name after a class keyword: class
    /// names such as `water` are also keywords, and are values here.
    fn at_class_name(&self, field: &Field) -> bool {
        self.peek_keyword()
            .is_some_and(|word| match field.value_type() {
                FieldType::ResidueClass => residue_class(word).is_some(),
                FieldType::MoleculeClass => molecule_class(word).is_some(),
                _ => false,
            })
    }

    fn primary(&mut self) -> Result<Node, FormatError> {
        let token = self.peek().cloned().ok_or_else(|| {
            syntax(
                self.input,
                self.end_span(),
                "expected a term before the end",
            )
        })?;
        let start = token.span.start;
        match &token.kind {
            TokenKind::Open => {
                self.position += 1;
                let mut inner = self.or()?;
                match self.peek() {
                    Some(Token {
                        kind: TokenKind::Close,
                        span,
                    }) => {
                        inner.span = Span(start..span.end);
                        self.position += 1;
                        Ok(inner)
                    }
                    Some(other) => Err(syntax(
                        self.input,
                        other.span.clone(),
                        "expected `)`, `and`, or `or`",
                    )),
                    None => Err(syntax(self.input, start..start + 1, "unclosed `(`")),
                }
            }
            TokenKind::Close => Err(syntax(self.input, token.span, "expected a term before `)`")),
            TokenKind::Op(_) => Err(syntax(
                self.input,
                token.span,
                "a comparison needs a keyword before it, such as `bfactor > 50`",
            )),
            TokenKind::Word { text, quoted: true } => Err(syntax(
                self.input,
                token.span.clone(),
                format!("expected a keyword, not the value {text:?}"),
            )),
            TokenKind::Word {
                text,
                quoted: false,
            } => {
                let text = text.clone();
                self.position += 1;
                let leaf = |kind| Node {
                    kind,
                    span: Span(token.span.clone()),
                };
                if text == "all" {
                    return Ok(leaf(NodeKind::All));
                }
                if text == "none" {
                    return Ok(leaf(NodeKind::None));
                }
                if let Some(flag) = flag(&text) {
                    return Ok(leaf(NodeKind::Flag(flag)));
                }
                if let Some(reference) = text.strip_prefix('#') {
                    if reference.is_empty() {
                        return Err(syntax(self.input, token.span, "`#` needs a reference"));
                    }
                    return Ok(leaf(NodeKind::Host(
                        HostTermKind::Reference,
                        reference.to_owned(),
                    )));
                }
                if let Some(name) = text.strip_prefix('$') {
                    if name.is_empty() {
                        return Err(syntax(self.input, token.span, "`$` needs a name"));
                    }
                    return Ok(leaf(NodeKind::Host(HostTermKind::Named, name.to_owned())));
                }
                if self.options.is_host_keyword(&text) {
                    return Ok(leaf(NodeKind::Host(HostTermKind::Keyword, text)));
                }
                match text.as_str() {
                    "smarts" => return self.smarts(start),
                    "property" => return self.property(start),
                    _ => {}
                }
                let Some(field) = field(&text) else {
                    return Err(syntax(
                        self.input,
                        token.span,
                        format!("unknown keyword `{text}`"),
                    ));
                };
                self.field(start, field)
            }
        }
    }

    /// The next word, quoted or not, as a raw value.
    fn word(&mut self, what: &str) -> Result<(String, bool, Range<usize>), FormatError> {
        match self.peek().cloned() {
            Some(Token {
                kind: TokenKind::Word { text, quoted },
                span,
            }) => {
                self.position += 1;
                Ok((text, quoted, span))
            }
            Some(token) => Err(syntax(self.input, token.span, format!("expected {what}"))),
            None => Err(syntax(
                self.input,
                self.end_span(),
                format!("expected {what} before the end"),
            )),
        }
    }

    fn smarts(&mut self, start: usize) -> Result<Node, FormatError> {
        let (text, _, span) = self.word("a SMARTS pattern")?;
        if let Err(error) = crate::smarts::parse_str(&text, Default::default()) {
            return Err(invalid(
                self.input,
                span.clone(),
                format!("invalid SMARTS pattern: {}", error.message()),
            )
            .with_source(error));
        }
        Ok(Node {
            kind: NodeKind::Smarts(text),
            span: Span(start..span.end),
        })
    }

    fn property(&mut self, start: usize) -> Result<Node, FormatError> {
        let (text, _, span) = self.word("a property key")?;
        let key = PropertyKey::new(text).map_err(|error| {
            invalid(self.input, span.clone(), "invalid property key").with_source(error)
        })?;
        let field = Field::Property(key.clone());
        if matches!(self.peek().map(|token| &token.kind), Some(TokenKind::Op(_))) || self.at_value()
        {
            return self.field(start, field);
        }
        Ok(Node {
            kind: NodeKind::PropertyTrue(key),
            span: Span(start..span.end),
        })
    }

    fn field(&mut self, start: usize, field: Field) -> Result<Node, FormatError> {
        if let Some(Token {
            kind: TokenKind::Op(op),
            span,
        }) = self.peek().cloned()
        {
            self.position += 1;
            if !field.is_ordered() && !matches!(op, CompareOp::Equal | CompareOp::NotEqual) {
                return Err(syntax(
                    self.input,
                    span,
                    "ordering comparisons apply only to numeric keywords",
                ));
            }
            let (text, quoted, value_span) = self.word("a value after the comparison")?;
            let value = self.literal(&field, &text, quoted, value_span.clone())?;
            return Ok(Node {
                kind: NodeKind::Compare { field, op, value },
                span: Span(start..value_span.end.max(self.previous_end())),
            });
        }
        let mut values = Vec::new();
        while self.at_value() || self.at_class_name(&field) {
            let (text, quoted, span) = self.word("a value")?;
            if field.is_ordered() && !quoted && self.peek_keyword() != Some("to") {
                if let Some(split) = text
                    .char_indices()
                    .skip(1)
                    .find(|&(_, c)| c == '-')
                    .map(|(i, _)| i)
                {
                    let (low, high) = (&text[..split], &text[split + 1..]);
                    let low_span = span.start..span.start + split;
                    let high_span = span.start + split + 1..span.end;
                    // A property range needs numbers on both sides; `a-b` is
                    // otherwise one text value.
                    let numeric = |text: &str| text.parse::<f64>().is_ok_and(f64::is_finite);
                    let splits = field.value_type() != FieldType::Property
                        || (numeric(low) && numeric(high));
                    if splits {
                        if let (Ok(low), Ok(high)) = (
                            self.literal(&field, low, false, low_span),
                            self.literal(&field, high, false, high_span),
                        ) {
                            values.push(range(low, high));
                            continue;
                        }
                    }
                }
            }
            // The value, with any unit after it, then possibly `to` and the
            // range's end.
            let low = self.literal(&field, &text, quoted, span.clone())?;
            if self.peek_keyword() == Some("to") {
                if !field.is_ordered() {
                    return Err(syntax(
                        self.input,
                        self.peek().map_or(span, |token| token.span.clone()),
                        "ranges apply only to numeric keywords",
                    ));
                }
                self.position += 1;
                let (high_text, high_quoted, high_span) = self.word("the end of a range")?;
                let high = self.literal(&field, &high_text, high_quoted, high_span)?;
                values.push(range(low, high));
                continue;
            }
            values.push(ValueMatch::Exact(low));
        }
        if values.is_empty() {
            let span = self
                .peek()
                .map_or(self.end_span(), |token| token.span.clone());
            return Err(syntax(self.input, span, "expected a value"));
        }
        Ok(Node {
            kind: NodeKind::Match { field, values },
            span: Span(start..self.previous_end()),
        })
    }

    fn previous_end(&self) -> usize {
        self.position
            .checked_sub(1)
            .and_then(|index| self.tokens.get(index))
            .map_or(0, |token| token.span.end)
    }

    fn literal(
        &mut self,
        field: &Field,
        text: &str,
        quoted: bool,
        span: Range<usize>,
    ) -> Result<Literal, FormatError> {
        let fail =
            |what: &str| invalid(self.input, span.clone(), format!("`{text}` is not {what}"));
        Ok(match field.value_type() {
            FieldType::Integer => Literal::Integer(text.parse().map_err(|_| fail("an integer"))?),
            FieldType::Real => Literal::Real {
                value: text
                    .parse::<f64>()
                    .ok()
                    .filter(|value| value.is_finite())
                    .ok_or_else(|| fail("a finite number"))?,
                unit: None,
            },
            FieldType::Text => Literal::Text(text.to_owned()),
            FieldType::Element => {
                Literal::Element(element(text).ok_or_else(|| fail("an element symbol"))?)
            }
            FieldType::ResId => {
                let (number, insertion_code) =
                    residue_number(text).ok_or_else(|| fail("a residue number"))?;
                Literal::ResId {
                    number,
                    insertion_code,
                }
            }
            FieldType::ResidueClass => {
                Literal::ResidueClass(residue_class(text).ok_or_else(|| fail("a residue class"))?)
            }
            FieldType::MoleculeClass => Literal::MoleculeClass(
                molecule_class(text).ok_or_else(|| fail("a molecule class"))?,
            ),
            FieldType::Property => {
                let number = if quoted {
                    None
                } else {
                    text.parse::<f64>().ok().filter(|value| value.is_finite())
                };
                let unit = match (number, self.peek_keyword().and_then(unit)) {
                    (Some(_), Some(unit)) => {
                        self.position += 1;
                        Some(unit)
                    }
                    _ => None,
                };
                Literal::Property {
                    text: text.to_owned(),
                    number,
                    unit,
                }
            }
        })
    }
}

/// A range whose property bounds share a unit written on only one of them,
/// as in `property q 1-2 e`.
fn range(mut low: Literal, mut high: Literal) -> ValueMatch {
    if let (
        Literal::Property {
            unit: low_unit,
            number: Some(_),
            ..
        },
        Literal::Property {
            unit: high_unit,
            number: Some(_),
            ..
        },
    ) = (&mut low, &mut high)
    {
        match (&*low_unit, &*high_unit) {
            (None, Some(unit)) => *low_unit = Some(*unit),
            (Some(unit), None) => *high_unit = Some(*unit),
            _ => {}
        }
    }
    ValueMatch::Range(low, high)
}

/// Appends `node`, splicing in the operands of a nested operator of the
/// same kind, so equal expressions have one form.
fn push_flat(operands: &mut Vec<Node>, node: Node, or: bool) {
    match node.kind {
        NodeKind::Or(inner) if or => operands.extend(inner),
        NodeKind::And(inner) if !or => operands.extend(inner),
        kind => operands.push(Node {
            kind,
            span: node.span,
        }),
    }
}
