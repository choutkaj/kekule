use std::fmt;

use crate::core::{Atom, Bond};
use crate::properties::{PropertyKey, PropertyValue, PropertyValueRef};

use super::{
    AtomSiteView, ChainView, InstanceAtomId, InstanceBondId, MoleculeInstanceView, ResidueView,
    Topology, TopologyAtomIndex, TopologyBondIndex,
};

/// One topology atom: its instance-qualified identity, dense index, and
/// chemistry, with perception, connectivity, hierarchy, and static
/// annotations of this occurrence.
///
/// Atom fields read through `Deref` (`atom.element`, `atom.formal_charge`).
#[derive(Clone, Copy)]
pub struct AtomView<'a> {
    topology: &'a Topology,
    id: InstanceAtomId,
    index: TopologyAtomIndex,
    atom: &'a Atom,
}

impl fmt::Debug for AtomView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("AtomView")
            .field("id", &self.id)
            .field("index", &self.index)
            .field("atom", self.atom)
            .finish()
    }
}

impl std::ops::Deref for AtomView<'_> {
    type Target = Atom;

    fn deref(&self) -> &Self::Target {
        self.atom
    }
}

impl<'a> AtomView<'a> {
    pub(super) fn new(topology: &'a Topology, index: TopologyAtomIndex) -> Self {
        let id = topology.atom_ids()[index.index()];
        let atom = topology
            .definition_molecule(id.molecule())
            .atom(id.atom())
            .expect("published topology atoms are live");
        Self {
            topology,
            id,
            index,
            atom,
        }
    }

    pub const fn id(self) -> InstanceAtomId {
        self.id
    }

    /// Position in the topology's dense atom order.
    pub const fn index(self) -> TopologyAtomIndex {
        self.index
    }

    /// The atom payload with the view's full lifetime.
    pub const fn atom(self) -> &'a Atom {
        self.atom
    }

    /// The molecule occurrence containing this atom.
    pub fn molecule(self) -> MoleculeInstanceView<'a> {
        MoleculeInstanceView::new(self.topology, self.id.molecule())
    }

    /// Bonded neighbors in local adjacency order.
    pub fn neighbors(self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let topology = self.topology;
        let molecule = self.id.molecule();
        topology
            .definition_molecule(molecule)
            .neighbors(self.id.atom())
            .expect("published topology atoms are live")
            .map(move |neighbor| {
                topology
                    .atom(InstanceAtomId::new(molecule, neighbor))
                    .expect("bonded neighbors are topology atoms")
            })
    }

    /// Incident bonds in local adjacency order.
    pub fn bonds(self) -> impl Iterator<Item = BondView<'a>> + 'a {
        let topology = self.topology;
        let molecule = self.id.molecule();
        topology
            .definition_molecule(molecule)
            .incident_bonds(self.id.atom())
            .expect("published topology atoms are live")
            .map(move |(bond, _)| {
                topology
                    .bond(InstanceBondId::new(molecule, bond))
                    .expect("incident bonds are topology bonds")
            })
    }

    /// Explicit graph hydrogen neighbors.
    pub fn explicit_hydrogens(self) -> usize {
        self.molecule()
            .molecule()
            .explicit_hydrogens(self.id.atom())
            .expect("published topology atoms are live")
    }

    /// Implicit hydrogens; see [`crate::core::Molecule::implicit_hydrogens`].
    pub fn implicit_hydrogens(self) -> Option<usize> {
        self.molecule()
            .molecule()
            .implicit_hydrogens(self.id.atom())
            .expect("published topology atoms are live")
    }

    /// Explicit plus implicit hydrogens; see [`crate::core::Molecule::total_hydrogens`].
    pub fn total_hydrogens(self) -> Option<usize> {
        self.molecule()
            .molecule()
            .total_hydrogens(self.id.atom())
            .expect("published topology atoms are live")
    }

    /// Only the inferred hydrogen contribution, for valence-model diagnostics.
    pub fn inferred_hydrogens(self) -> Option<u8> {
        self.molecule()
            .molecule()
            .inferred_hydrogens(self.id.atom())
            .expect("published topology atoms are live")
    }

    /// Perceived aromaticity; `None` before perception.
    pub fn is_aromatic(self) -> Option<bool> {
        self.molecule()
            .molecule()
            .atom_is_aromatic(self.id.atom())
            .expect("published topology atoms are live")
    }

    pub fn atom_site(self) -> Option<AtomSiteView<'a>> {
        self.topology
            .hierarchy()
            .atom_site_for_atom(self.id)
            .map(|site| AtomSiteView::new(self.topology, site.id()))
    }

    pub fn residue(self) -> Option<ResidueView<'a>> {
        self.atom_site().map(AtomSiteView::residue)
    }

    pub fn chain(self) -> Option<ChainView<'a>> {
        self.residue().map(ResidueView::chain)
    }

    /// One static topology annotation of this atom.
    pub fn property(self, key: &PropertyKey) -> Option<PropertyValue> {
        self.property_ref(key).map(PropertyValueRef::to_value)
    }

    /// Borrows one static annotation without copying string storage.
    pub fn property_ref(self, key: &PropertyKey) -> Option<PropertyValueRef<'a>> {
        self.topology
            .properties()
            .atoms()
            .value_ref(key, self.index)
            .expect("dense atom rows cover every atom")
    }
}

/// One topology bond: its instance-qualified identity, dense index, and
/// chemistry, with endpoint atoms, perception, and static annotations.
///
/// Bond fields read through `Deref` (`bond.order`).
#[derive(Clone, Copy)]
pub struct BondView<'a> {
    topology: &'a Topology,
    id: InstanceBondId,
    index: TopologyBondIndex,
    bond: &'a Bond,
}

impl fmt::Debug for BondView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("BondView")
            .field("id", &self.id)
            .field("index", &self.index)
            .field("bond", self.bond)
            .finish()
    }
}

impl std::ops::Deref for BondView<'_> {
    type Target = Bond;

    fn deref(&self) -> &Self::Target {
        self.bond
    }
}

impl<'a> BondView<'a> {
    pub(super) fn new(topology: &'a Topology, index: TopologyBondIndex) -> Self {
        let id = topology.bond_ids()[index.index()];
        let bond = topology
            .definition_molecule(id.molecule())
            .bond(id.bond())
            .expect("published topology bonds are live");
        Self {
            topology,
            id,
            index,
            bond,
        }
    }

    pub const fn id(self) -> InstanceBondId {
        self.id
    }

    /// Position in the topology's dense bond order.
    pub const fn index(self) -> TopologyBondIndex {
        self.index
    }

    /// The bond payload with the view's full lifetime.
    pub const fn bond(self) -> &'a Bond {
        self.bond
    }

    pub fn molecule(self) -> MoleculeInstanceView<'a> {
        MoleculeInstanceView::new(self.topology, self.id.molecule())
    }

    /// Both endpoint atoms in stored bond order.
    pub fn atoms(self) -> [AtomView<'a>; 2] {
        let (a, b) = self.bond.endpoints();
        let molecule = self.id.molecule();
        [a, b].map(|atom| {
            self.topology
                .atom(InstanceAtomId::new(molecule, atom))
                .expect("bond endpoints are topology atoms")
        })
    }

    /// Perceived aromaticity; `None` before perception.
    pub fn is_aromatic(self) -> Option<bool> {
        self.molecule()
            .molecule()
            .bond_is_aromatic(self.id.bond())
            .expect("published topology bonds are live")
    }

    /// One static topology annotation of this bond.
    pub fn property(self, key: &PropertyKey) -> Option<PropertyValue> {
        self.property_ref(key).map(PropertyValueRef::to_value)
    }

    /// Borrows one static annotation without copying string storage.
    pub fn property_ref(self, key: &PropertyKey) -> Option<PropertyValueRef<'a>> {
        self.topology
            .properties()
            .bonds()
            .value_ref(key, self.index)
            .expect("dense bond rows cover every bond")
    }
}
