use std::collections::BTreeMap;

use super::*;

pub type PropMap = BTreeMap<String, PropValue>;

#[derive(Debug, Clone, PartialEq)]
pub enum PropValue {
    String(String),
    Int(i64),
    Float(f64),
    Bool(bool),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub isotope: Option<u16>,
    pub formal_charge: i8,
    pub radical: Option<AtomRadical>,
    pub hydrogens: HydrogenDeclaration,
    pub atom_map: Option<u32>,
    pub props: PropMap,
}

impl Atom {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            isotope: None,
            formal_charge: 0,
            radical: None,
            hydrogens: HydrogenDeclaration::default(),
            atom_map: None,
            props: PropMap::new(),
        }
    }
}

/// The complete non-graph hydrogen statement represented on an atom.
///
/// Graph hydrogen atoms are separate atoms and are not counted here. Hydrogens
/// inferred by valence perception are stored in [`PerceptionState`] rather
/// than this declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HydrogenDeclaration {
    /// The represented count is present and valence perception may infer
    /// additional implicit hydrogens.
    Infer { explicit: u8 },
    /// Exactly this many non-graph hydrogens are represented.
    Fixed(u8),
}

impl HydrogenDeclaration {
    pub const fn explicit_count(self) -> u8 {
        match self {
            Self::Infer { explicit } | Self::Fixed(explicit) => explicit,
        }
    }

    pub const fn allows_implicit(self) -> bool {
        matches!(self, Self::Infer { .. })
    }

    /// Returns the same inference policy with a different represented count.
    pub const fn with_explicit_count(self, explicit: u8) -> Self {
        match self {
            Self::Infer { .. } => Self::Infer { explicit },
            Self::Fixed(_) => Self::Fixed(explicit),
        }
    }
}

impl Default for HydrogenDeclaration {
    fn default() -> Self {
        Self::Infer { explicit: 0 }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomRadical {
    Singlet,
    Doublet,
    Triplet,
    Quartet,
    Quintet,
}

impl AtomRadical {
    pub const fn unpaired_electron_count(self) -> u8 {
        match self {
            Self::Singlet => 0,
            Self::Doublet => 1,
            Self::Triplet => 2,
            Self::Quartet => 3,
            Self::Quintet => 4,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub(crate) a: AtomId,
    pub(crate) b: AtomId,
    pub order: BondOrder,
    pub props: PropMap,
}

impl Bond {
    pub fn new(a: AtomId, b: AtomId, order: BondOrder) -> Self {
        Self {
            a,
            b,
            order,
            props: PropMap::new(),
        }
    }

    pub const fn a(&self) -> AtomId {
        self.a
    }

    pub const fn b(&self) -> AtomId {
        self.b
    }

    pub const fn endpoints(&self) -> (AtomId, AtomId) {
        (self.a, self.b)
    }
}

/// A canonical localized bond order stored in a [`Molecule`](super::Molecule).
///
/// Aromaticity is perceived state, not a represented bond order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Zero,
    Single,
    Double,
    Triple,
    Quadruple,
    Dative,
}
