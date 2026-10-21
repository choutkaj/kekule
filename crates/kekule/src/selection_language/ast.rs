//! The parsed form of a selection expression.

use std::ops::Range;

use super::SelectionExpr;
use crate::molecule::Element;
use crate::properties::PropertyKey;
use crate::topology::{MoleculeClass, ResidueClass};
use crate::units::Unit;

/// A byte range of the parsed text. Spans locate evaluation errors; they do
/// not take part in expression equality.
#[derive(Debug, Clone, Default)]
pub(crate) struct Span(pub(crate) Range<usize>);

impl PartialEq for Span {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

/// One node of a selection expression.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Node {
    pub(crate) kind: NodeKind,
    pub(crate) span: Span,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum NodeKind {
    All,
    None,
    /// Two or more operands; nested conjunctions are flattened.
    And(Vec<Node>),
    /// Two or more operands; nested disjunctions are flattened.
    Or(Vec<Node>),
    Not(Box<Node>),
    Within {
        distance: f64,
        unit: DistanceUnit,
        /// An expression of its own, so a resolver can evaluate it on other
        /// structures.
        reference: Box<SelectionExpr>,
    },
    Expand {
        level: Level,
        inner: Box<Node>,
    },
    Bonded(Box<Node>),
    Flag(Flag),
    /// A field matching any of its values.
    Match {
        field: Field,
        values: Vec<ValueMatch>,
    },
    Compare {
        field: Field,
        op: CompareOp,
        value: Literal,
    },
    /// `property key` alone: a boolean property that is true.
    PropertyTrue(PropertyKey),
    Smarts(String),
    Host(HostTermKind, String),
}

/// The unit of a `within` distance; `Angstrom` is the default and prints
/// without a unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DistanceUnit {
    Angstrom,
    Nanometer,
    Picometer,
    Bohr,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Level {
    Residue,
    Chain,
    Molecule,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Flag {
    Aromatic,
    Ring,
    Hydrogen,
    Heavy,
    PolarHydrogen,
    Protein,
    Nucleic,
    Water,
    Ion,
    Ligand,
    Backbone,
    Sidechain,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Field {
    Index,
    Element,
    Name,
    LabelName,
    Charge,
    Isotope,
    BFactor,
    Occupancy,
    Property(PropertyKey),
    ResName,
    ResId,
    LabelResId,
    ResClass,
    Chain,
    LabelChain,
    MolClass,
    Instance,
    Definition,
}

/// What kind of value a field takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FieldType {
    Integer,
    Real,
    Text,
    Element,
    ResId,
    ResidueClass,
    MoleculeClass,
    /// A property: interpreted per column type at evaluation.
    Property,
}

impl Field {
    pub(crate) fn value_type(&self) -> FieldType {
        match self {
            Self::Index | Self::Charge | Self::Isotope | Self::Instance | Self::Definition => {
                FieldType::Integer
            }
            Self::LabelResId => FieldType::Integer,
            Self::BFactor | Self::Occupancy => FieldType::Real,
            Self::Name | Self::LabelName | Self::ResName | Self::Chain | Self::LabelChain => {
                FieldType::Text
            }
            Self::Element => FieldType::Element,
            Self::ResId => FieldType::ResId,
            Self::ResClass => FieldType::ResidueClass,
            Self::MolClass => FieldType::MoleculeClass,
            Self::Property(_) => FieldType::Property,
        }
    }

    /// Whether ordering comparisons and ranges apply.
    pub(crate) fn is_ordered(&self) -> bool {
        matches!(
            self.value_type(),
            FieldType::Integer | FieldType::Real | FieldType::ResId | FieldType::Property
        )
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum ValueMatch {
    Exact(Literal),
    /// Inclusive bounds.
    Range(Literal, Literal),
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Literal {
    Integer(i64),
    /// A number with an optional explicit unit and its spelling, for
    /// printing.
    Real {
        value: f64,
        unit: Option<(Unit, &'static str)>,
    },
    Text(String),
    Element(Element),
    ResId {
        number: i64,
        insertion_code: Option<String>,
    },
    ResidueClass(ResidueClass),
    MoleculeClass(MoleculeClass),
    /// A property value whose type the column decides: the word as written,
    /// its number when it is an unquoted number, and a unit.
    Property {
        text: String,
        number: Option<f64>,
        unit: Option<(Unit, &'static str)>,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CompareOp {
    Less,
    LessEqual,
    Greater,
    GreaterEqual,
    Equal,
    NotEqual,
}

impl CompareOp {
    pub(crate) const fn symbol(self) -> &'static str {
        match self {
            Self::Less => "<",
            Self::LessEqual => "<=",
            Self::Greater => ">",
            Self::GreaterEqual => ">=",
            Self::Equal => "==",
            Self::NotEqual => "!=",
        }
    }

    pub(crate) fn holds(self, ordering: std::cmp::Ordering) -> bool {
        use std::cmp::Ordering::{Equal, Greater, Less};
        match self {
            Self::Less => ordering == Less,
            Self::LessEqual => ordering != Greater,
            Self::Greater => ordering == Greater,
            Self::GreaterEqual => ordering != Less,
            Self::Equal => ordering == Equal,
            Self::NotEqual => ordering != Equal,
        }
    }
}

/// The syntax of a host term.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum HostTermKind {
    /// `#...`, such as `#3` or `#2-5`: the text after `#`.
    Reference,
    /// `$name`: the text after `$`.
    Named,
    /// A bare keyword registered with
    /// [`SelectParseOptions::with_host_keywords`](super::SelectParseOptions::with_host_keywords).
    Keyword,
}
