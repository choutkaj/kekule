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
    pub explicit_hydrogens: u8,
    pub implicit_hydrogens: Option<u8>,
    pub no_implicit_hydrogens: bool,
    pub aromatic: bool,
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
            explicit_hydrogens: 0,
            implicit_hydrogens: None,
            no_implicit_hydrogens: false,
            aromatic: false,
            atom_map: None,
            props: PropMap::new(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomRadical {
    Singlet,
    Doublet,
    Triplet,
}

impl AtomRadical {
    pub const fn unpaired_electron_count(self) -> u8 {
        match self {
            Self::Singlet => 0,
            Self::Doublet => 1,
            Self::Triplet => 2,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub(crate) a: AtomId,
    pub(crate) b: AtomId,
    pub order: BondOrder,
    pub aromatic: bool,
    pub props: PropMap,
}

impl Bond {
    pub fn new(a: AtomId, b: AtomId, order: BondOrder) -> Self {
        Self {
            a,
            b,
            order,
            aromatic: matches!(order, BondOrder::Aromatic),
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Zero,
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Dative,
}
