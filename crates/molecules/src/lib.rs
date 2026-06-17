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
    a: AtomId,
    b: AtomId,
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
    adjacency: Vec<Vec<BondId>>,
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
        self.adjacency.push(Vec::new());
        self.invalidate_topology();
        id
    }

    pub fn delete_atom(&mut self, id: AtomId) -> Result<Atom> {
        self.atom(id)?;
        let incident = self.adjacency[id.index()].clone();
        for bond_id in incident {
            if self
                .bonds
                .get(bond_id.index())
                .and_then(Option::as_ref)
                .is_some()
            {
                self.delete_bond(bond_id)?;
            }
        }
        self.adjacency[id.index()].clear();
        let atom = self.atoms[id.index()]
            .take()
            .ok_or(MoleculeError::InvalidAtomId(id))?;
        self.invalidate_topology();
        Ok(atom)
    }

    pub fn atom(&self, id: AtomId) -> Result<&Atom> {
        self.atoms
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidAtomId(id))
    }

    pub fn atom_mut(&mut self, id: AtomId) -> Result<&mut Atom> {
        self.atom(id)?;
        self.invalidate_topology();
        self.atoms[id.index()]
            .as_mut()
            .ok_or(MoleculeError::InvalidAtomId(id))
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms
            .iter()
            .enumerate()
            .filter_map(|(index, atom)| atom.as_ref().map(|atom| (AtomId::new(index as u32), atom)))
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms().map(|(id, _)| id)
    }

    pub fn add_bond(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<BondId> {
        self.atom(a)?;
        self.atom(b)?;
        if a == b {
            return Err(MoleculeError::SelfBond(a));
        }
        if self.bond_between(a, b)?.is_some() {
            return Err(MoleculeError::DuplicateBond { a, b });
        }
        let id = BondId::new(self.bonds.len() as u32);
        self.bonds.push(Some(Bond::new(a, b, order)));
        self.adjacency[a.index()].push(id);
        self.adjacency[b.index()].push(id);
        self.invalidate_topology();
        Ok(id)
    }

    pub fn delete_bond(&mut self, id: BondId) -> Result<Bond> {
        let bond = self
            .bonds
            .get_mut(id.index())
            .and_then(Option::take)
            .ok_or(MoleculeError::InvalidBondId(id))?;
        self.remove_incident_bond(bond.a, id);
        self.remove_incident_bond(bond.b, id);
        self.invalidate_topology();
        Ok(bond)
    }

    pub fn bond(&self, id: BondId) -> Result<&Bond> {
        self.bonds
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidBondId(id))
    }

    pub fn bond_mut(&mut self, id: BondId) -> Result<&mut Bond> {
        self.bond(id)?;
        self.invalidate_topology();
        self.bonds[id.index()]
            .as_mut()
            .ok_or(MoleculeError::InvalidBondId(id))
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)> {
        self.bonds
            .iter()
            .enumerate()
            .filter_map(|(index, bond)| bond.as_ref().map(|bond| (BondId::new(index as u32), bond)))
    }

    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bonds().map(|(id, _)| id)
    }

    pub fn neighbors(&self, id: AtomId) -> Result<impl Iterator<Item = AtomId> + '_> {
        self.atom(id)?;
        Ok(self.adjacency[id.index()]
            .iter()
            .filter_map(|bond_id| self.bond(*bond_id).ok())
            .map(move |bond| bond.other_atom(id)))
    }

    pub fn incident_bonds(&self, id: AtomId) -> Result<impl Iterator<Item = (BondId, &Bond)> + '_> {
        self.atom(id)?;
        Ok(self.adjacency[id.index()]
            .iter()
            .filter_map(|bond_id| self.bond(*bond_id).ok().map(|bond| (*bond_id, bond))))
    }

    pub fn bond_between(&self, a: AtomId, b: AtomId) -> Result<Option<BondId>> {
        self.atom(a)?;
        self.atom(b)?;
        Ok(self.adjacency[a.index()].iter().copied().find(|bond_id| {
            self.bond(*bond_id)
                .map(|bond| bond.connects(a, b))
                .unwrap_or(false)
        }))
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

    fn remove_incident_bond(&mut self, atom: AtomId, bond: BondId) {
        if let Some(incident) = self.adjacency.get_mut(atom.index()) {
            incident.retain(|id| *id != bond);
        }
    }
}

impl Bond {
    fn connects(&self, a: AtomId, b: AtomId) -> bool {
        (self.a == a && self.b == b) || (self.a == b && self.b == a)
    }

    fn other_atom(&self, atom: AtomId) -> AtomId {
        if self.a == atom {
            self.b
        } else {
            self.a
        }
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

    fn oxygen() -> Atom {
        Atom::new(Element::from_symbol("O").expect("oxygen should be available"))
    }

    fn sorted_atom_ids(ids: impl IntoIterator<Item = AtomId>) -> Vec<AtomId> {
        let mut ids = ids.into_iter().collect::<Vec<_>>();
        ids.sort();
        ids
    }

    fn sorted_bond_ids(ids: impl IntoIterator<Item = BondId>) -> Vec<BondId> {
        let mut ids = ids.into_iter().collect::<Vec<_>>();
        ids.sort();
        ids
    }

    fn mark_all_fresh(mol: &mut Molecule) {
        mol.perception_mut().valence = ComputedState::Fresh;
        mol.perception_mut().rings = ComputedState::Fresh;
        mol.perception_mut().aromaticity = ComputedState::Fresh;
        mol.perception_mut().stereo = ComputedState::Fresh;
    }

    fn assert_all_stale(mol: &Molecule) {
        assert_eq!(mol.perception().valence, ComputedState::Stale);
        assert_eq!(mol.perception().rings, ComputedState::Stale);
        assert_eq!(mol.perception().aromaticity, ComputedState::Stale);
        assert_eq!(mol.perception().stereo, ComputedState::Stale);
    }

    #[test]
    fn empty_molecule_has_no_atoms_or_bonds() {
        let mol = Molecule::new();

        assert_eq!(mol.atom_count(), 0);
        assert_eq!(mol.bond_count(), 0);
        assert!(mol.atoms().next().is_none());
        assert!(mol.bonds().next().is_none());
    }

    #[test]
    fn atom_insertion_assigns_stable_typed_ids() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());

        assert_eq!(a.raw(), 0);
        assert_eq!(b.raw(), 1);
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(
            mol.atom(a).expect("first atom exists").element.symbol(),
            "C"
        );
        assert_eq!(
            mol.atom(b).expect("second atom exists").element.symbol(),
            "O"
        );
        assert_eq!(mol.atom_ids().collect::<Vec<_>>(), vec![a, b]);
    }

    #[test]
    fn bond_insertion_assigns_stable_typed_ids() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");

        assert_eq!(bond.raw(), 0);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(
            mol.bond(bond).expect("bond should exist").endpoints(),
            (a, b)
        );
        assert_eq!(mol.bond_ids().collect::<Vec<_>>(), vec![bond]);
    }

    #[test]
    fn invalid_atom_ids_are_rejected() {
        let mut mol = Molecule::new();
        let atom = mol.add_atom(carbon());

        assert_eq!(
            mol.atom(AtomId::new(99))
                .expect_err("missing atom should fail"),
            MoleculeError::InvalidAtomId(AtomId::new(99))
        );
        mol.delete_atom(atom).expect("atom should delete");
        assert_eq!(
            mol.atom(atom).expect_err("deleted atom should fail"),
            MoleculeError::InvalidAtomId(atom)
        );
        assert_eq!(
            mol.add_bond(atom, AtomId::new(99), BondOrder::Single)
                .expect_err("deleted endpoint should fail"),
            MoleculeError::InvalidAtomId(atom)
        );
    }

    #[test]
    fn invalid_bond_ids_are_rejected() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");

        assert_eq!(
            mol.bond(BondId::new(99))
                .expect_err("missing bond should fail"),
            MoleculeError::InvalidBondId(BondId::new(99))
        );
        mol.delete_bond(bond).expect("bond should delete");
        assert_eq!(
            mol.bond(bond).expect_err("deleted bond should fail"),
            MoleculeError::InvalidBondId(bond)
        );
        assert_eq!(
            mol.delete_bond(bond)
                .expect_err("deleting bond twice should fail"),
            MoleculeError::InvalidBondId(bond)
        );
    }

    #[test]
    fn self_bonds_are_rejected() {
        let mut mol = Molecule::new();
        let atom = mol.add_atom(carbon());

        let err = mol
            .add_bond(atom, atom, BondOrder::Single)
            .expect_err("self-bond should fail");
        assert_eq!(err, MoleculeError::SelfBond(atom));
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

        let reverse_err = mol
            .add_bond(b, a, BondOrder::Double)
            .expect_err("reverse duplicate should fail");
        assert_eq!(reverse_err, MoleculeError::DuplicateBond { a: b, b: a });
    }

    #[test]
    fn neighbor_iteration_reports_live_adjacent_atoms() {
        let mut mol = Molecule::new();
        let center = mol.add_atom(carbon());
        let left = mol.add_atom(carbon());
        let right = mol.add_atom(oxygen());
        let isolated = mol.add_atom(carbon());
        mol.add_bond(center, left, BondOrder::Single)
            .expect("left bond should be valid");
        mol.add_bond(center, right, BondOrder::Double)
            .expect("right bond should be valid");

        assert_eq!(
            sorted_atom_ids(mol.neighbors(center).expect("center exists")),
            vec![left, right]
        );
        assert_eq!(
            mol.neighbors(isolated)
                .expect("isolated atom exists")
                .collect::<Vec<_>>(),
            Vec::<AtomId>::new()
        );
        match mol.neighbors(AtomId::new(99)) {
            Ok(_) => panic!("missing atom should fail"),
            Err(err) => assert_eq!(err, MoleculeError::InvalidAtomId(AtomId::new(99))),
        };
    }

    #[test]
    fn incident_bond_iteration_reports_live_bonds() {
        let mut mol = Molecule::new();
        let center = mol.add_atom(carbon());
        let left = mol.add_atom(carbon());
        let right = mol.add_atom(oxygen());
        let left_bond = mol
            .add_bond(center, left, BondOrder::Single)
            .expect("left bond should be valid");
        let right_bond = mol
            .add_bond(center, right, BondOrder::Double)
            .expect("right bond should be valid");

        assert_eq!(
            sorted_bond_ids(
                mol.incident_bonds(center)
                    .expect("center exists")
                    .map(|(id, _)| id)
            ),
            vec![left_bond, right_bond]
        );

        mol.delete_bond(left_bond).expect("left bond should delete");
        assert_eq!(
            mol.incident_bonds(center)
                .expect("center still exists")
                .map(|(id, _)| id)
                .collect::<Vec<_>>(),
            vec![right_bond]
        );
        match mol.incident_bonds(AtomId::new(99)) {
            Ok(_) => panic!("missing atom should fail"),
            Err(err) => assert_eq!(err, MoleculeError::InvalidAtomId(AtomId::new(99))),
        };
    }

    #[test]
    fn bond_between_finds_live_undirected_bonds() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());
        let c = mol.add_atom(carbon());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");

        assert_eq!(mol.bond_between(a, b).expect("atoms exist"), Some(bond));
        assert_eq!(mol.bond_between(b, a).expect("atoms exist"), Some(bond));
        assert_eq!(mol.bond_between(a, c).expect("atoms exist"), None);
    }

    #[test]
    fn bond_deletion_preserves_remaining_ids_and_counts() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let c = mol.add_atom(oxygen());
        let first = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("first bond should be valid");
        let second = mol
            .add_bond(b, c, BondOrder::Double)
            .expect("second bond should be valid");

        let removed = mol.delete_bond(first).expect("first bond should delete");

        assert_eq!(removed.a(), a);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(mol.bond(first), Err(MoleculeError::InvalidBondId(first)));
        assert_eq!(
            mol.bond(second).expect("second bond remains").order,
            BondOrder::Double
        );
        assert_eq!(mol.bond_ids().collect::<Vec<_>>(), vec![second]);
        assert_eq!(
            mol.neighbors(b)
                .expect("middle atom exists")
                .collect::<Vec<_>>(),
            vec![c]
        );
    }

    #[test]
    fn atom_deletion_removes_incident_bonds_and_preserves_remaining_ids() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let c = mol.add_atom(oxygen());
        let first = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("first bond should be valid");
        let second = mol
            .add_bond(b, c, BondOrder::Double)
            .expect("second bond should be valid");

        let removed = mol.delete_atom(b).expect("middle atom should delete");

        assert_eq!(removed.element.symbol(), "C");
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 0);
        assert_eq!(mol.atom(b), Err(MoleculeError::InvalidAtomId(b)));
        assert_eq!(
            mol.atom(a).expect("first atom remains").element.symbol(),
            "C"
        );
        assert_eq!(
            mol.atom(c).expect("third atom remains").element.symbol(),
            "O"
        );
        assert_eq!(mol.bond(first), Err(MoleculeError::InvalidBondId(first)));
        assert_eq!(mol.bond(second), Err(MoleculeError::InvalidBondId(second)));
        assert_eq!(mol.atom_ids().collect::<Vec<_>>(), vec![a, c]);
        assert_eq!(
            mol.neighbors(a)
                .expect("first atom exists")
                .collect::<Vec<_>>(),
            Vec::<AtomId>::new()
        );
    }

    #[test]
    fn adding_after_deletion_allocates_new_ids() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let first_bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");
        mol.delete_bond(first_bond).expect("bond should delete");
        mol.delete_atom(a).expect("atom should delete");

        let c = mol.add_atom(oxygen());
        let second_bond = mol
            .add_bond(b, c, BondOrder::Double)
            .expect("new bond should be valid");

        assert_eq!(c.raw(), 2);
        assert_eq!(second_bond.raw(), 1);
        assert_eq!(mol.atom_ids().collect::<Vec<_>>(), vec![b, c]);
        assert_eq!(mol.bond_ids().collect::<Vec<_>>(), vec![second_bond]);
    }

    #[test]
    fn every_topology_mutation_invalidates_fresh_perception() {
        let mut mol = Molecule::new();
        mark_all_fresh(&mut mol);
        let a = mol.add_atom(carbon());
        assert_all_stale(&mol);

        mark_all_fresh(&mut mol);
        let b = mol.add_atom(oxygen());
        assert_all_stale(&mol);

        mark_all_fresh(&mut mol);
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");
        assert_all_stale(&mol);

        mark_all_fresh(&mut mol);
        mol.delete_bond(bond).expect("bond should delete");
        assert_all_stale(&mol);

        mark_all_fresh(&mut mol);
        mol.delete_atom(a).expect("atom should delete");
        assert_all_stale(&mol);
    }

    #[test]
    fn absent_perception_remains_absent_after_topology_mutation() {
        let mut mol = Molecule::new();

        mol.add_atom(carbon());

        assert_eq!(mol.perception().valence, ComputedState::Absent);
        assert_eq!(mol.perception().rings, ComputedState::Absent);
        assert_eq!(mol.perception().aromaticity, ComputedState::Absent);
        assert_eq!(mol.perception().stereo, ComputedState::Absent);
    }

    #[test]
    fn property_maps_can_be_mutated_without_topology_changes() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");
        mol.props_mut().insert(
            "name".to_owned(),
            PropValue::String("carbon monoxide".to_owned()),
        );
        mol.atom_mut(a)
            .expect("atom exists")
            .props
            .insert("role".to_owned(), PropValue::String("donor".to_owned()));
        mol.bond_mut(bond)
            .expect("bond exists")
            .props
            .insert("source".to_owned(), PropValue::Bool(true));

        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(
            mol.props().get("name"),
            Some(&PropValue::String("carbon monoxide".to_owned()))
        );
        assert_eq!(
            mol.atom(a).expect("atom exists").props.get("role"),
            Some(&PropValue::String("donor".to_owned()))
        );
        assert_eq!(
            mol.bond(bond).expect("bond exists").props.get("source"),
            Some(&PropValue::Bool(true))
        );
    }

    #[test]
    fn wrappers_share_the_core_molecule_graph() {
        let mut small = SmallMolecule::default();
        let a = small.mol.add_atom(carbon());
        let b = small.mol.add_atom(oxygen());
        small
            .mol
            .add_bond(a, b, BondOrder::Single)
            .expect("small molecule graph should accept bonds");

        let mut macro_mol = MacroMolecule::default();
        let c = macro_mol.mol.add_atom(carbon());

        assert_eq!(small.mol.atom_count(), 2);
        assert_eq!(small.mol.bond_count(), 1);
        assert_eq!(
            macro_mol
                .mol
                .atom(c)
                .expect("macro atom exists")
                .element
                .symbol(),
            "C"
        );
    }
}
