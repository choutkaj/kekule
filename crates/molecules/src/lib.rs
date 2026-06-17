#![forbid(unsafe_code)]

use std::collections::BTreeMap;
use std::fmt;

pub mod prelude {
    pub use crate::{
        Atom, AtomId, AtomStereo, BioHierarchy, Bond, BondId, BondOrder, BondStereo, ComputedState,
        Element, MacroMolecule, Molecule, MoleculeError, PropMap, PropValue, Result, SmallMolecule,
    };
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomId(u32);

impl AtomId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl fmt::Display for AtomId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "a{}", self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondId(u32);

impl BondId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl fmt::Display for BondId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "b{}", self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Element {
    atomic_number: u8,
}

impl Element {
    pub fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        if (1..=118).contains(&atomic_number) {
            Some(Self { atomic_number })
        } else {
            None
        }
    }

    pub fn from_symbol(symbol: &str) -> Option<Self> {
        let atomic_number = match symbol {
            "H" => 1,
            "B" => 5,
            "C" => 6,
            "N" => 7,
            "O" => 8,
            "F" => 9,
            "P" => 15,
            "S" => 16,
            "Cl" => 17,
            "Br" => 35,
            "I" => 53,
            _ => return None,
        };
        Some(Self { atomic_number })
    }

    pub const fn atomic_number(self) -> u8 {
        self.atomic_number
    }

    pub fn symbol(self) -> &'static str {
        match self.atomic_number {
            1 => "H",
            5 => "B",
            6 => "C",
            7 => "N",
            8 => "O",
            9 => "F",
            15 => "P",
            16 => "S",
            17 => "Cl",
            35 => "Br",
            53 => "I",
            _ => "?",
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

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
    pub radical_electrons: u8,
    pub explicit_hydrogens: u8,
    pub implicit_hydrogens: Option<u8>,
    pub aromatic: bool,
    pub chiral: Option<AtomStereo>,
    pub atom_map: Option<u32>,
    pub props: PropMap,
}

impl Atom {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            isotope: None,
            formal_charge: 0,
            radical_electrons: 0,
            explicit_hydrogens: 0,
            implicit_hydrogens: None,
            aromatic: false,
            chiral: None,
            atom_map: None,
            props: PropMap::new(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomStereo {
    TetrahedralClockwise,
    TetrahedralCounterClockwise,
    Unspecified,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub a: AtomId,
    pub b: AtomId,
    pub order: BondOrder,
    pub aromatic: bool,
    pub stereo: Option<BondStereo>,
    pub props: PropMap,
}

impl Bond {
    pub fn new(a: AtomId, b: AtomId, order: BondOrder) -> Self {
        Self {
            a,
            b,
            order,
            aromatic: matches!(order, BondOrder::Aromatic),
            stereo: None,
            props: PropMap::new(),
        }
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondStereo {
    E,
    Z,
    Up,
    Down,
    Unspecified,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub enum ComputedState {
    #[default]
    Absent,
    Stale,
    Fresh,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PerceptionState {
    pub valence: ComputedState,
    pub rings: ComputedState,
    pub aromaticity: ComputedState,
    pub stereo: ComputedState,
}

impl PerceptionState {
    pub fn invalidate_all(&mut self) {
        self.valence = invalidate(self.valence);
        self.rings = invalidate(self.rings);
        self.aromaticity = invalidate(self.aromaticity);
        self.stereo = invalidate(self.stereo);
    }
}

fn invalidate(state: ComputedState) -> ComputedState {
    match state {
        ComputedState::Fresh => ComputedState::Stale,
        ComputedState::Stale | ComputedState::Absent => state,
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Molecule {
    atoms: Vec<Option<Atom>>,
    bonds: Vec<Option<Bond>>,
    props: PropMap,
    perception: PerceptionState,
}

impl Molecule {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.iter().flatten().count()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.iter().flatten().count()
    }

    pub fn add_atom(&mut self, atom: Atom) -> AtomId {
        let id = AtomId::new(self.atoms.len() as u32);
        self.atoms.push(Some(atom));
        self.invalidate_topology();
        id
    }

    pub fn atom(&self, id: AtomId) -> Result<&Atom> {
        self.atoms
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidAtomId(id))
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms
            .iter()
            .enumerate()
            .filter_map(|(index, atom)| atom.as_ref().map(|atom| (AtomId::new(index as u32), atom)))
    }

    pub fn add_bond(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<BondId> {
        self.atom(a)?;
        self.atom(b)?;
        if a == b {
            return Err(MoleculeError::SelfBond(a));
        }
        if self.has_bond_between(a, b) {
            return Err(MoleculeError::DuplicateBond { a, b });
        }
        let id = BondId::new(self.bonds.len() as u32);
        self.bonds.push(Some(Bond::new(a, b, order)));
        self.invalidate_topology();
        Ok(id)
    }

    pub fn bond(&self, id: BondId) -> Result<&Bond> {
        self.bonds
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidBondId(id))
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)> {
        self.bonds
            .iter()
            .enumerate()
            .filter_map(|(index, bond)| bond.as_ref().map(|bond| (BondId::new(index as u32), bond)))
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    pub fn perception(&self) -> &PerceptionState {
        &self.perception
    }

    pub fn perception_mut(&mut self) -> &mut PerceptionState {
        &mut self.perception
    }

    pub fn invalidate_topology(&mut self) {
        self.perception.invalidate_all();
    }

    fn has_bond_between(&self, a: AtomId, b: AtomId) -> bool {
        self.bonds
            .iter()
            .flatten()
            .any(|bond| (bond.a == a && bond.b == b) || (bond.a == b && bond.b == a))
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SmallMolecule {
    pub mol: Molecule,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct BioHierarchy {
    pub props: PropMap,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MacroMolecule {
    pub mol: Molecule,
    pub hierarchy: BioHierarchy,
}

pub type Result<T> = std::result::Result<T, MoleculeError>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoleculeError {
    InvalidAtomId(AtomId),
    InvalidBondId(BondId),
    SelfBond(AtomId),
    DuplicateBond { a: AtomId, b: AtomId },
    UnsupportedFeature(&'static str),
}

impl fmt::Display for MoleculeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(id) => write!(f, "invalid atom id: {id}"),
            Self::InvalidBondId(id) => write!(f, "invalid bond id: {id}"),
            Self::SelfBond(id) => write!(f, "cannot create a bond from atom {id} to itself"),
            Self::DuplicateBond { a, b } => write!(f, "duplicate bond between {a} and {b}"),
            Self::UnsupportedFeature(name) => write!(f, "unsupported feature: {name}"),
        }
    }
}

impl std::error::Error for MoleculeError {}

#[cfg(test)]
mod tests {
    use super::*;

    fn carbon() -> Atom {
        Atom::new(Element::from_symbol("C").expect("carbon should be available"))
    }

    #[test]
    fn add_atoms_and_bond() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");

        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(mol.bond(bond).expect("bond should exist").a, a);
    }

    #[test]
    fn duplicate_bond_is_rejected() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        mol.add_bond(a, b, BondOrder::Single)
            .expect("first bond should be valid");

        let err = mol
            .add_bond(a, b, BondOrder::Double)
            .expect_err("duplicate should fail");
        assert_eq!(err, MoleculeError::DuplicateBond { a, b });
    }

    #[test]
    fn fresh_perception_becomes_stale_after_topology_mutation() {
        let mut mol = Molecule::new();
        mol.perception_mut().rings = ComputedState::Fresh;
        mol.add_atom(carbon());
        assert_eq!(mol.perception().rings, ComputedState::Stale);
    }
}
