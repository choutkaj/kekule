use std::fmt;
use std::ops::{Deref, DerefMut};

use crate::algorithms::{RingMembership, RingSet};

use super::*;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub(crate) enum ComputedState {
    #[default]
    Absent,
    Stale,
    Fresh,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct PerceptionState {
    pub valence: ComputedState,
    pub rings: ComputedState,
    pub aromaticity: ComputedState,
    pub stereo: ComputedState,
}

impl PerceptionState {
    pub(crate) fn invalidate_all(&mut self) {
        self.valence = invalidate(self.valence);
        self.rings = invalidate(self.rings);
        self.aromaticity = invalidate(self.aromaticity);
        self.stereo = invalidate(self.stereo);
    }
}

pub(crate) fn invalidate(state: ComputedState) -> ComputedState {
    match state {
        ComputedState::Fresh => ComputedState::Stale,
        ComputedState::Stale | ComputedState::Absent => state,
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Molecule {
    pub(crate) atoms: Vec<Option<Atom>>,
    pub(crate) bonds: Vec<Option<Bond>>,
    pub(crate) adjacency: Vec<Vec<BondId>>,
    pub(crate) conformers: Vec<Option<Conformer>>,
    pub(crate) props: PropMap,
    pub(crate) perception: PerceptionState,
    pub(crate) ring_membership: Option<RingMembership>,
    pub(crate) ring_set: Option<RingSet>,
}

pub struct AtomMut<'a> {
    molecule: &'a mut Molecule,
    id: AtomId,
    original: AtomChemistry,
}

impl Deref for AtomMut<'_> {
    type Target = Atom;

    fn deref(&self) -> &Self::Target {
        self.molecule.atoms[self.id.index()]
            .as_ref()
            .expect("validated atom must remain live while borrowed")
    }
}

impl DerefMut for AtomMut<'_> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.molecule.atoms[self.id.index()]
            .as_mut()
            .expect("validated atom must remain live while borrowed")
    }
}

impl Drop for AtomMut<'_> {
    fn drop(&mut self) {
        if AtomChemistry::from(&**self) != self.original {
            self.molecule.invalidate_topology();
        }
    }
}

pub struct BondMut<'a> {
    molecule: &'a mut Molecule,
    id: BondId,
    original: BondChemistry,
}

impl Deref for BondMut<'_> {
    type Target = Bond;

    fn deref(&self) -> &Self::Target {
        self.molecule.bonds[self.id.index()]
            .as_ref()
            .expect("validated bond must remain live while borrowed")
    }
}

impl DerefMut for BondMut<'_> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.molecule.bonds[self.id.index()]
            .as_mut()
            .expect("validated bond must remain live while borrowed")
    }
}

impl Drop for BondMut<'_> {
    fn drop(&mut self) {
        if BondChemistry::from(&**self) != self.original {
            self.molecule.invalidate_topology();
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct AtomChemistry {
    element: Element,
    isotope: Option<u16>,
    formal_charge: i8,
    radical: Option<AtomRadical>,
    explicit_hydrogens: u8,
    implicit_hydrogens: Option<u8>,
    no_implicit_hydrogens: bool,
    aromatic: bool,
    chiral: Option<AtomStereo>,
}

impl From<&Atom> for AtomChemistry {
    fn from(atom: &Atom) -> Self {
        Self {
            element: atom.element,
            isotope: atom.isotope,
            formal_charge: atom.formal_charge,
            radical: atom.radical,
            explicit_hydrogens: atom.explicit_hydrogens,
            implicit_hydrogens: atom.implicit_hydrogens,
            no_implicit_hydrogens: atom.no_implicit_hydrogens,
            aromatic: atom.aromatic,
            chiral: atom.chiral,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct BondChemistry {
    order: BondOrder,
    aromatic: bool,
    stereo: Option<BondStereo>,
}

impl From<&Bond> for BondChemistry {
    fn from(bond: &Bond) -> Self {
        Self {
            order: bond.order,
            aromatic: bond.aromatic,
            stereo: bond.stereo,
        }
    }
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
        for conformer in self.conformers.iter_mut().flatten() {
            conformer.clear_position(id);
        }
        self.invalidate_topology();
        Ok(atom)
    }

    pub fn atom(&self, id: AtomId) -> Result<&Atom> {
        self.atoms
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidAtomId(id))
    }

    pub fn atom_mut(&mut self, id: AtomId) -> Result<AtomMut<'_>> {
        let original = AtomChemistry::from(self.atom(id)?);
        Ok(AtomMut {
            molecule: self,
            id,
            original,
        })
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

    pub fn bond_mut(&mut self, id: BondId) -> Result<BondMut<'_>> {
        let original = BondChemistry::from(self.bond(id)?);
        Ok(BondMut {
            molecule: self,
            id,
            original,
        })
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

    #[cfg(test)]
    pub(crate) fn perception(&self) -> &PerceptionState {
        &self.perception
    }

    pub fn ring_membership(&self) -> Option<&RingMembership> {
        if self.perception.rings == ComputedState::Fresh {
            self.ring_membership.as_ref()
        } else {
            None
        }
    }

    pub fn ring_set(&self) -> Option<&RingSet> {
        if self.perception.rings == ComputedState::Fresh {
            self.ring_set.as_ref()
        } else {
            None
        }
    }

    pub fn add_conformer(&mut self, mut conformer: Conformer) -> ConformerId {
        if conformer.positions.len() < self.atoms.len() {
            conformer.positions.resize(self.atoms.len(), None);
        }
        let id = ConformerId::new(self.conformers.len() as u32);
        self.conformers.push(Some(conformer));
        id
    }

    pub fn conformer(&self, id: ConformerId) -> Result<&Conformer> {
        self.conformers
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidConformerId(id))
    }

    pub fn conformer_mut(&mut self, id: ConformerId) -> Result<&mut Conformer> {
        self.conformers
            .get_mut(id.index())
            .and_then(Option::as_mut)
            .ok_or(MoleculeError::InvalidConformerId(id))
    }

    pub fn conformers(&self) -> impl Iterator<Item = (ConformerId, &Conformer)> {
        self.conformers
            .iter()
            .enumerate()
            .filter_map(|(index, conformer)| {
                conformer
                    .as_ref()
                    .map(|conformer| (ConformerId::new(index as u32), conformer))
            })
    }

    pub fn first_conformer(&self) -> Option<(ConformerId, &Conformer)> {
        self.conformers().next()
    }

    pub fn invalidate_topology(&mut self) {
        self.perception.invalidate_all();
        self.ring_membership = None;
        self.ring_set = None;
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

    pub(crate) fn other_atom(&self, atom: AtomId) -> AtomId {
        if self.a == atom {
            self.b
        } else {
            self.a
        }
    }
}

pub type Result<T> = std::result::Result<T, MoleculeError>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoleculeError {
    InvalidAtomId(AtomId),
    InvalidBondId(BondId),
    InvalidConformerId(ConformerId),
    SelfBond(AtomId),
    DuplicateBond { a: AtomId, b: AtomId },
    UnsupportedFeature(&'static str),
}

impl fmt::Display for MoleculeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(id) => write!(f, "invalid atom id: {id}"),
            Self::InvalidBondId(id) => write!(f, "invalid bond id: {id}"),
            Self::InvalidConformerId(id) => write!(f, "invalid conformer id: {id}"),
            Self::SelfBond(id) => write!(f, "cannot create a bond from atom {id} to itself"),
            Self::DuplicateBond { a, b } => write!(f, "duplicate bond between {a} and {b}"),
            Self::UnsupportedFeature(name) => write!(f, "unsupported feature: {name}"),
        }
    }
}

impl std::error::Error for MoleculeError {}
