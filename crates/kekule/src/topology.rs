//! Immutable coordinate-free molecular systems, qualified identities, dense
//! orderings, and compiled selections.

mod hierarchy;
pub mod transform;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::sync::Arc;

use crate::core::{Atom, AtomId, Bond, BondId, Element, Molecule, MoleculeConnectivityError};
use crate::properties::{Properties, PropertyError, PropertyKey, PropertyTable, PropertyValue};
use crate::substructure::QueryMatch;
pub use hierarchy::{
    AtomSite, AtomSiteId, AtomSiteMetadata, Chain, ChainId, Hierarchy, HierarchyError,
    HierarchyIdKind, Residue, ResidueId,
};

fixed_u32_id!(MoleculeDefinitionId, "definition");
fixed_u32_id!(MoleculeInstanceId, "molecule");
fixed_u32_id!(TopologyAtomIndex, "atom-index");
fixed_u32_id!(TopologyBondIndex, "bond-index");

/// The local atom of one explicit molecule instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InstanceAtomId {
    molecule: MoleculeInstanceId,
    atom: AtomId,
}

impl InstanceAtomId {
    pub const fn new(molecule: MoleculeInstanceId, atom: AtomId) -> Self {
        Self { molecule, atom }
    }

    pub const fn molecule(self) -> MoleculeInstanceId {
        self.molecule
    }

    pub const fn atom(self) -> AtomId {
        self.atom
    }
}

impl fmt::Display for InstanceAtomId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}:{}", self.molecule, self.atom)
    }
}

/// The local bond of one explicit molecule instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InstanceBondId {
    molecule: MoleculeInstanceId,
    bond: BondId,
}

impl InstanceBondId {
    pub const fn new(molecule: MoleculeInstanceId, bond: BondId) -> Self {
        Self { molecule, bond }
    }

    pub const fn molecule(self) -> MoleculeInstanceId {
        self.molecule
    }

    pub const fn bond(self) -> BondId {
        self.bond
    }
}

impl fmt::Display for InstanceBondId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}:{}", self.molecule, self.bond)
    }
}

/// A topology-global hierarchy chain borrowed from a [`Topology`].
#[derive(Clone, Copy)]
pub struct ChainView<'a> {
    topology: &'a Topology,
    id: ChainId,
}

impl fmt::Debug for ChainView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ChainView")
            .field("id", &self.id)
            .finish()
    }
}

impl<'a> ChainView<'a> {
    const fn new(topology: &'a Topology, id: ChainId) -> Self {
        Self { topology, id }
    }

    pub const fn id(self) -> ChainId {
        self.id
    }

    pub fn label_id(self) -> &'a str {
        self.local().label_id()
    }

    pub fn author_id(self) -> Option<&'a str> {
        self.local().author_id()
    }

    pub fn residues(self) -> impl ExactSizeIterator<Item = ResidueView<'a>> + 'a {
        let topology = self.topology;
        self.local()
            .residues()
            .iter()
            .copied()
            .map(move |residue| ResidueView::new(topology, residue))
    }

    pub fn property(self, key: &PropertyKey) -> Result<Option<PropertyValue>, PropertyError> {
        self.topology
            .properties
            .chains()
            .value(key, self.id.index())
    }

    pub(crate) fn local(self) -> &'a Chain {
        self.topology
            .hierarchy
            .chain(self.id)
            .expect("chain view references a validated topology hierarchy")
    }
}

/// A topology-global hierarchy residue borrowed from a [`Topology`].
#[derive(Clone, Copy)]
pub struct ResidueView<'a> {
    topology: &'a Topology,
    id: ResidueId,
}

impl fmt::Debug for ResidueView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ResidueView")
            .field("id", &self.id)
            .finish()
    }
}

impl<'a> ResidueView<'a> {
    const fn new(topology: &'a Topology, id: ResidueId) -> Self {
        Self { topology, id }
    }

    pub const fn id(self) -> ResidueId {
        self.id
    }

    pub fn chain(self) -> ChainView<'a> {
        ChainView::new(self.topology, self.local().chain())
    }

    pub fn name(self) -> &'a str {
        self.local().name()
    }

    pub fn label_comp_id(self) -> Option<&'a str> {
        self.local().label_comp_id()
    }

    pub fn author_comp_id(self) -> Option<&'a str> {
        self.local().author_comp_id()
    }

    pub fn label_seq_id(self) -> Option<i32> {
        self.local().label_seq_id()
    }

    pub fn author_seq_id(self) -> Option<&'a str> {
        self.local().author_seq_id()
    }

    pub fn insertion_code(self) -> Option<&'a str> {
        self.local().insertion_code()
    }

    pub fn atom_sites(self) -> impl ExactSizeIterator<Item = AtomSiteView<'a>> + 'a {
        let topology = self.topology;
        self.local()
            .atom_sites()
            .iter()
            .copied()
            .map(move |site| AtomSiteView::new(topology, site))
    }

    pub fn property(self, key: &PropertyKey) -> Result<Option<PropertyValue>, PropertyError> {
        self.topology
            .properties
            .residues()
            .value(key, self.id.index())
    }

    pub(crate) fn local(self) -> &'a Residue {
        self.topology
            .hierarchy
            .residue(self.id)
            .expect("residue view references a validated topology hierarchy")
    }
}

/// A topology-global hierarchy atom site borrowed from a [`Topology`].
#[derive(Clone, Copy)]
pub struct AtomSiteView<'a> {
    topology: &'a Topology,
    id: AtomSiteId,
}

impl fmt::Debug for AtomSiteView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("AtomSiteView")
            .field("id", &self.id)
            .finish()
    }
}

impl<'a> AtomSiteView<'a> {
    const fn new(topology: &'a Topology, id: AtomSiteId) -> Self {
        Self { topology, id }
    }

    pub const fn id(self) -> AtomSiteId {
        self.id
    }

    pub fn atom(self) -> InstanceAtomId {
        self.local().atom()
    }

    pub fn residue(self) -> ResidueView<'a> {
        ResidueView::new(self.topology, self.local().residue())
    }

    pub fn metadata(self) -> &'a AtomSiteMetadata {
        self.local().metadata()
    }

    pub fn property(self, key: &PropertyKey) -> Result<Option<PropertyValue>, PropertyError> {
        self.topology
            .properties
            .atom_sites()
            .value(key, self.id.index())
    }

    pub(crate) fn local(self) -> &'a AtomSite {
        self.topology
            .hierarchy
            .atom_site(self.id)
            .expect("atom-site view references a validated topology hierarchy")
    }
}

/// One reusable coordinate-free molecule definition.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeDefinition {
    id: MoleculeDefinitionId,
    molecule: Molecule,
}

impl MoleculeDefinition {
    pub const fn id(&self) -> MoleculeDefinitionId {
        self.id
    }

    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
}

/// One explicit occurrence of a reusable molecule definition.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeInstance {
    id: MoleculeInstanceId,
    definition: MoleculeDefinitionId,
}

impl MoleculeInstance {
    pub const fn id(&self) -> MoleculeInstanceId {
        self.id
    }

    pub const fn definition(&self) -> MoleculeDefinitionId {
        self.definition
    }

    pub const fn qualify_atom(&self, atom: AtomId) -> InstanceAtomId {
        InstanceAtomId::new(self.id, atom)
    }

    pub const fn qualify_bond(&self, bond: BondId) -> InstanceBondId {
        InstanceBondId::new(self.id, bond)
    }
}

/// One explicit molecule occurrence borrowed from a [`Topology`].
///
/// This is the instance-first system view. The underlying [`Molecule`] retains
/// definition-local identities, while atoms and bonds reached through this
/// view are qualified by this occurrence's [`MoleculeInstanceId`]. Hierarchy
/// methods are projections over the one topology-owned hierarchy.
#[derive(Clone, Copy)]
pub struct MoleculeInstanceView<'a> {
    topology: &'a Topology,
    id: MoleculeInstanceId,
}

impl fmt::Debug for MoleculeInstanceView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MoleculeInstanceView")
            .field("id", &self.id)
            .field("definition", &self.definition_id())
            .finish()
    }
}

impl<'a> MoleculeInstanceView<'a> {
    const fn new(topology: &'a Topology, id: MoleculeInstanceId) -> Self {
        Self { topology, id }
    }

    /// Returns this occurrence's topology-wide molecule identity.
    pub const fn id(self) -> MoleculeInstanceId {
        self.id
    }

    /// Returns the reusable definition referenced by this occurrence.
    pub fn definition_id(self) -> MoleculeDefinitionId {
        self.instance().definition()
    }

    /// Returns the minimal stored instance record.
    pub fn instance(self) -> &'a MoleculeInstance {
        &self.topology.instances[self.id.index()]
    }

    /// Returns the reusable definition referenced by this occurrence.
    pub fn definition(self) -> &'a MoleculeDefinition {
        &self.topology.definitions[self.definition_id().index()]
    }

    /// Returns the definition-owned molecular state for this occurrence.
    pub fn molecule(self) -> &'a Molecule {
        self.definition().molecule()
    }

    pub fn property(self, key: &PropertyKey) -> Result<Option<PropertyValue>, PropertyError> {
        self.topology
            .properties
            .molecule_instances()
            .value(key, self.id.index())
    }

    /// Iterates atoms with topology-wide, instance-qualified identities.
    pub fn atoms(self) -> impl Iterator<Item = (InstanceAtomId, &'a Atom)> + 'a {
        let instance = self.id;
        self.molecule()
            .atoms()
            .map(move |(atom, payload)| (InstanceAtomId::new(instance, atom), payload))
    }

    /// Iterates bonds with topology-wide, instance-qualified identities.
    pub fn bonds(self) -> impl Iterator<Item = (InstanceBondId, &'a Bond)> + 'a {
        let instance = self.id;
        self.molecule()
            .bonds()
            .map(move |(bond, payload)| (InstanceBondId::new(instance, bond), payload))
    }

    /// Iterates complete chains that contain at least one atom in this molecule.
    /// A returned chain is not clipped and may also contain atoms from other
    /// molecule instances.
    pub fn chains(self) -> impl Iterator<Item = ChainView<'a>> + 'a {
        let molecule = self.id;
        self.topology.chains().filter(move |chain| {
            chain
                .residues()
                .flat_map(ResidueView::atom_sites)
                .any(|site| site.atom().molecule() == molecule)
        })
    }

    /// Iterates complete residues that contain at least one atom in this molecule.
    /// A returned residue is not clipped and may contain sites from other
    /// molecule instances.
    pub fn residues(self) -> impl Iterator<Item = ResidueView<'a>> + 'a {
        let molecule = self.id;
        self.topology.residues().filter(move |residue| {
            residue
                .atom_sites()
                .any(|site| site.atom().molecule() == molecule)
        })
    }

    /// Iterates atom sites whose atoms belong to this molecule instance.
    pub fn atom_sites(self) -> impl Iterator<Item = AtomSiteView<'a>> + 'a {
        let molecule = self.id;
        self.topology
            .atom_sites()
            .filter(move |site| site.atom().molecule() == molecule)
    }

    pub const fn qualify_atom(self, atom: AtomId) -> InstanceAtomId {
        InstanceAtomId::new(self.id, atom)
    }

    pub const fn qualify_bond(self, bond: BondId) -> InstanceBondId {
        InstanceBondId::new(self.id, bond)
    }
}

/// An immutable, coordinate-free molecular system.
///
/// `Topology` owns molecule instances, their dense identity and ordering, and
/// the system [`Hierarchy`]. Shared exact ownership uses [`Arc<Topology>`].
/// Topology-changing operations construct new topology values; generic
/// topology remapping is not part of the core architecture.
#[derive(Debug)]
pub struct Topology {
    definitions: Vec<MoleculeDefinition>,
    instances: Vec<MoleculeInstance>,
    instance_atoms: Vec<InstanceAtomId>,
    instance_bonds: Vec<InstanceBondId>,
    atom_indices: BTreeMap<InstanceAtomId, TopologyAtomIndex>,
    bond_indices: BTreeMap<InstanceBondId, TopologyBondIndex>,
    hierarchy: Hierarchy,
    properties: Properties,
}

impl Topology {
    pub fn builder() -> TopologyBuilder {
        TopologyBuilder::new()
    }

    /// Returns whether two topologies have the same complete static layout.
    ///
    /// Layout equality includes chemical and hierarchy content, definition and
    /// instance partitioning, semantic identifiers,
    /// authoritative dense atom and bond order, and the corresponding index
    /// maps. Whether two values share one `Arc` allocation is deliberately
    /// excluded.
    ///
    /// This is stricter than order-independent structural equivalence. It does
    /// not perform graph isomorphism, reorder definitions or instances, or
    /// resolve repeated indistinguishable content.
    pub fn same_layout(&self, other: &Self) -> bool {
        self.definitions == other.definitions
            && self.instances == other.instances
            && self.instance_atoms == other.instance_atoms
            && self.instance_bonds == other.instance_bonds
            && self.atom_indices == other.atom_indices
            && self.bond_indices == other.bond_indices
            && self.hierarchy == other.hierarchy
    }

    pub fn definition(
        &self,
        id: MoleculeDefinitionId,
    ) -> Result<&MoleculeDefinition, TopologyError> {
        self.definitions
            .get(id.index())
            .ok_or(TopologyError::InvalidMoleculeDefinitionId(id))
    }

    pub fn definitions(
        &self,
    ) -> impl ExactSizeIterator<Item = (MoleculeDefinitionId, &MoleculeDefinition)> {
        self.definitions
            .iter()
            .map(|definition| (definition.id, definition))
    }

    pub fn definition_count(&self) -> usize {
        self.definitions.len()
    }

    pub fn instance(&self, id: MoleculeInstanceId) -> Result<&MoleculeInstance, TopologyError> {
        self.instances
            .get(id.index())
            .ok_or(TopologyError::InvalidMoleculeInstanceId(id))
    }

    pub fn instances(
        &self,
    ) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &MoleculeInstance)> {
        self.instances
            .iter()
            .map(|instance| (instance.id, instance))
    }

    pub fn instance_count(&self) -> usize {
        self.instances.len()
    }

    /// Returns one instance-qualified molecular view.
    pub fn molecule(
        &self,
        id: MoleculeInstanceId,
    ) -> Result<MoleculeInstanceView<'_>, TopologyError> {
        self.instance(id)?;
        Ok(MoleculeInstanceView::new(self, id))
    }

    /// Iterates explicit molecules in authoritative instance order.
    pub fn molecules(&self) -> impl ExactSizeIterator<Item = MoleculeInstanceView<'_>> {
        self.instances
            .iter()
            .map(|instance| MoleculeInstanceView::new(self, instance.id))
    }

    pub fn definition_for_instance(
        &self,
        instance: MoleculeInstanceId,
    ) -> Result<&MoleculeDefinition, TopologyError> {
        let instance = self.instance(instance)?;
        self.definition(instance.definition)
    }

    pub fn instances_for_definition(
        &self,
        definition: MoleculeDefinitionId,
    ) -> Result<impl Iterator<Item = &MoleculeInstance>, TopologyError> {
        self.definition(definition)?;
        Ok(self
            .instances
            .iter()
            .filter(move |instance| instance.definition == definition))
    }

    /// Returns the one authoritative system-level hierarchy.
    pub const fn hierarchy(&self) -> &Hierarchy {
        &self.hierarchy
    }

    pub const fn properties(&self) -> &Properties {
        &self.properties
    }

    pub const fn molecule_instance_properties(&self) -> &PropertyTable {
        self.properties.molecule_instances()
    }

    pub const fn atom_properties(&self) -> &PropertyTable {
        self.properties.atoms()
    }

    pub const fn bond_properties(&self) -> &PropertyTable {
        self.properties.bonds()
    }

    pub const fn chain_properties(&self) -> &PropertyTable {
        self.properties.chains()
    }

    pub const fn residue_properties(&self) -> &PropertyTable {
        self.properties.residues()
    }

    pub const fn atom_site_properties(&self) -> &PropertyTable {
        self.properties.atom_sites()
    }

    pub fn molecule_instance_property(
        &self,
        instance: MoleculeInstanceId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyError> {
        self.instance(instance)?;
        self.properties
            .molecule_instances()
            .value(key, instance.index())
            .map_err(|error| TopologyError::Property(Box::new(error)))
    }

    pub fn atom_property(
        &self,
        atom: InstanceAtomId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyError> {
        let index = self
            .atom_index(atom)
            .ok_or(TopologyError::InvalidAtomId(atom))?;
        self.properties
            .atoms()
            .value(key, index.index())
            .map_err(|error| TopologyError::Property(Box::new(error)))
    }

    pub fn bond_property(
        &self,
        bond: InstanceBondId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyError> {
        let index = self
            .bond_index(bond)
            .ok_or(TopologyError::InvalidBondId(bond))?;
        self.properties
            .bonds()
            .value(key, index.index())
            .map_err(|error| TopologyError::Property(Box::new(error)))
    }

    /// Iterates every topology-global hierarchy chain in hierarchy order.
    pub fn chains(&self) -> impl Iterator<Item = ChainView<'_>> {
        self.hierarchy
            .chains()
            .map(move |(id, _)| ChainView::new(self, id))
    }

    /// Iterates every topology-global hierarchy residue in hierarchy order.
    pub fn residues(&self) -> impl Iterator<Item = ResidueView<'_>> {
        self.hierarchy
            .residues()
            .map(move |(id, _)| ResidueView::new(self, id))
    }

    /// Iterates every topology-global hierarchy atom site in hierarchy order.
    pub fn atom_sites(&self) -> impl Iterator<Item = AtomSiteView<'_>> {
        self.hierarchy
            .atom_sites()
            .map(move |(id, _)| AtomSiteView::new(self, id))
    }

    pub fn chain(&self, id: ChainId) -> Result<ChainView<'_>, TopologyError> {
        self.hierarchy
            .chain(id)
            .map_err(|_| TopologyError::InvalidChainId(id))?;
        Ok(ChainView::new(self, id))
    }

    pub fn residue(&self, id: ResidueId) -> Result<ResidueView<'_>, TopologyError> {
        self.hierarchy
            .residue(id)
            .map_err(|_| TopologyError::InvalidResidueId(id))?;
        Ok(ResidueView::new(self, id))
    }

    pub fn atom_site(&self, id: AtomSiteId) -> Result<AtomSiteView<'_>, TopologyError> {
        self.hierarchy
            .atom_site(id)
            .map_err(|_| TopologyError::InvalidAtomSiteId(id))?;
        Ok(AtomSiteView::new(self, id))
    }

    pub fn atom_for_site(&self, site: AtomSiteId) -> Result<InstanceAtomId, TopologyError> {
        let atom = self.atom_site(site)?.atom();
        self.atom(atom)?;
        Ok(atom)
    }

    pub fn atom_site_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<AtomSiteView<'_>>, TopologyError> {
        self.atom(atom)?;
        Ok(self
            .hierarchy
            .atom_site_for_atom(atom)
            .map(|site| AtomSiteView::new(self, site.id())))
    }

    pub fn residue_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<ResidueView<'_>>, TopologyError> {
        let Some(site) = self.atom_site_for_atom(atom)? else {
            return Ok(None);
        };
        Ok(Some(site.residue()))
    }

    pub fn chain_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<ChainView<'_>>, TopologyError> {
        let Some(residue) = self.residue_for_atom(atom)? else {
            return Ok(None);
        };
        Ok(Some(residue.chain()))
    }

    pub fn residue_for_site(&self, site: AtomSiteId) -> Result<ResidueView<'_>, TopologyError> {
        Ok(self.atom_site(site)?.residue())
    }

    pub fn chain_for_residue(&self, residue: ResidueId) -> Result<ChainView<'_>, TopologyError> {
        Ok(self.residue(residue)?.chain())
    }

    pub fn atom(&self, id: InstanceAtomId) -> Result<&Atom, TopologyError> {
        self.molecule(id.molecule)?
            .molecule()
            .atom(id.atom)
            .map_err(|_| TopologyError::InvalidAtomId(id))
    }

    pub fn bond(&self, id: InstanceBondId) -> Result<&Bond, TopologyError> {
        self.molecule(id.molecule)?
            .molecule()
            .bond(id.bond)
            .map_err(|_| TopologyError::InvalidBondId(id))
    }

    pub fn atoms(&self) -> impl ExactSizeIterator<Item = (InstanceAtomId, &Atom)> {
        self.instance_atoms.iter().copied().map(|id| {
            (
                id,
                self.atom(id)
                    .expect("built topology instance atoms contain only live atoms"),
            )
        })
    }

    pub fn bonds(&self) -> impl ExactSizeIterator<Item = (InstanceBondId, &Bond)> {
        self.instance_bonds.iter().copied().map(|id| {
            (
                id,
                self.bond(id)
                    .expect("built topology instance bonds contain only live bonds"),
            )
        })
    }

    pub fn atom_count(&self) -> usize {
        self.instance_atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.instance_bonds.len()
    }

    pub fn atom_ids(&self) -> &[InstanceAtomId] {
        &self.instance_atoms
    }

    pub fn bond_ids(&self) -> &[InstanceBondId] {
        &self.instance_bonds
    }

    pub fn atom_index(&self, atom: InstanceAtomId) -> Option<TopologyAtomIndex> {
        self.atom_indices.get(&atom).copied()
    }

    pub fn atom_id(&self, index: TopologyAtomIndex) -> Option<InstanceAtomId> {
        self.instance_atoms.get(index.index()).copied()
    }

    pub fn bond_index(&self, bond: InstanceBondId) -> Option<TopologyBondIndex> {
        self.bond_indices.get(&bond).copied()
    }

    pub fn bond_id(&self, index: TopologyBondIndex) -> Option<InstanceBondId> {
        self.instance_bonds.get(index.index()).copied()
    }

    pub fn neighbors(
        &self,
        atom: InstanceAtomId,
    ) -> Result<impl Iterator<Item = InstanceAtomId> + '_, TopologyError> {
        let molecule = self.molecule(atom.molecule)?.molecule();
        molecule
            .atom(atom.atom)
            .map_err(|_| TopologyError::InvalidAtomId(atom))?;
        Ok(molecule
            .neighbors(atom.atom)
            .expect("validated atom has valid local adjacency")
            .map(move |neighbor| InstanceAtomId::new(atom.molecule, neighbor)))
    }

    pub fn incident_bonds(
        &self,
        atom: InstanceAtomId,
    ) -> Result<impl Iterator<Item = (InstanceBondId, &Bond)> + '_, TopologyError> {
        let molecule = self.molecule(atom.molecule)?.molecule();
        molecule
            .atom(atom.atom)
            .map_err(|_| TopologyError::InvalidAtomId(atom))?;
        Ok(molecule
            .incident_bonds(atom.atom)
            .expect("validated atom has valid local adjacency")
            .map(move |(bond, payload)| (InstanceBondId::new(atom.molecule, bond), payload)))
    }

    pub fn implicit_hydrogens(&self, atom: InstanceAtomId) -> Result<Option<u8>, TopologyError> {
        self.atom(atom)?;
        self.molecule(atom.molecule)?
            .molecule()
            .implicit_hydrogens(atom.atom)
            .map_err(|_| TopologyError::InvalidAtomId(atom))
    }

    pub fn atom_is_aromatic(&self, atom: InstanceAtomId) -> Result<Option<bool>, TopologyError> {
        self.atom(atom)?;
        self.molecule(atom.molecule)?
            .molecule()
            .atom_is_aromatic(atom.atom)
            .map_err(|_| TopologyError::InvalidAtomId(atom))
    }

    pub fn bond_is_aromatic(&self, bond: InstanceBondId) -> Result<Option<bool>, TopologyError> {
        self.bond(bond)?;
        self.molecule(bond.molecule)?
            .molecule()
            .bond_is_aromatic(bond.bond)
            .map_err(|_| TopologyError::InvalidBondId(bond))
    }
}

/// Linear, validate-then-commit builder for coordinate-free topology.
///
/// Molecule definitions and instances are staged before the one system-level
/// hierarchy is validated against their final instance-qualified atom IDs.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TopologyBuilder {
    definitions: Vec<MoleculeDefinition>,
    instances: Vec<MoleculeInstance>,
    hierarchy: Hierarchy,
    properties: Properties,
}

impl TopologyBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns staged system-level hierarchy state.
    pub const fn hierarchy(&self) -> &Hierarchy {
        &self.hierarchy
    }

    /// Returns mutable staged hierarchy state.
    ///
    /// References are checked transactionally by [`Self::build`]; published
    /// topologies never expose mutable hierarchy access.
    pub fn hierarchy_mut(&mut self) -> &mut Hierarchy {
        &mut self.hierarchy
    }

    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, PropertyError> {
        self.properties.insert(key, value)
    }

    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.properties.remove(key)
    }

    pub fn clear_properties(&mut self) {
        self.properties.clear_owner();
    }

    pub fn molecule_instance_properties_mut(&mut self) -> &mut PropertyTable {
        self.sync_property_dimensions();
        self.properties.molecule_instances_mut()
    }

    pub fn atom_properties_mut(&mut self) -> &mut PropertyTable {
        self.sync_property_dimensions();
        self.properties.atoms_mut()
    }

    pub fn bond_properties_mut(&mut self) -> &mut PropertyTable {
        self.sync_property_dimensions();
        self.properties.bonds_mut()
    }

    pub fn chain_properties_mut(&mut self) -> &mut PropertyTable {
        self.sync_property_dimensions();
        self.properties.chains_mut()
    }

    pub fn residue_properties_mut(&mut self) -> &mut PropertyTable {
        self.sync_property_dimensions();
        self.properties.residues_mut()
    }

    pub fn atom_site_properties_mut(&mut self) -> &mut PropertyTable {
        self.sync_property_dimensions();
        self.properties.atom_sites_mut()
    }

    pub fn reserve_definitions(&mut self, additional: usize) -> Result<(), TopologyBuildError> {
        checked_future_len(
            self.definitions.len(),
            additional,
            TopologyIdKind::MoleculeDefinition,
        )?;
        self.definitions.try_reserve(additional).map_err(|_| {
            TopologyBuildError::IdentifierCapacityExceeded(TopologyIdKind::MoleculeDefinition)
        })
    }

    pub fn reserve_instances(&mut self, additional: usize) -> Result<(), TopologyBuildError> {
        checked_future_len(
            self.instances.len(),
            additional,
            TopologyIdKind::MoleculeInstance,
        )?;
        self.instances.try_reserve(additional).map_err(|_| {
            TopologyBuildError::IdentifierCapacityExceeded(TopologyIdKind::MoleculeInstance)
        })
    }

    pub fn definition(
        &self,
        id: MoleculeDefinitionId,
    ) -> Result<&MoleculeDefinition, TopologyBuildError> {
        self.definitions
            .get(id.index())
            .ok_or(TopologyBuildError::InvalidMoleculeDefinitionId(id))
    }

    pub fn add_molecule_definition(
        &mut self,
        molecule: &Molecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        validate_graph(molecule)?;
        self.commit_definition(molecule.clone())
    }

    pub fn add_molecule_definition_owned(
        &mut self,
        molecule: Molecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        validate_graph(&molecule)?;
        self.commit_definition(molecule)
    }

    pub fn add_instance(
        &mut self,
        definition: MoleculeDefinitionId,
    ) -> Result<MoleculeInstanceId, TopologyBuildError> {
        self.definition(definition)?;
        self.reserve_instances(1)?;
        let id = checked_id::<MoleculeInstanceId>(
            self.instances.len(),
            TopologyIdKind::MoleculeInstance,
        )?;
        self.instances.push(MoleculeInstance { id, definition });
        Ok(id)
    }

    /// Adds one fresh definition and one instance in a single operation.
    pub fn add_molecule(
        &mut self,
        molecule: &Molecule,
    ) -> Result<MoleculeInstanceId, TopologyBuildError> {
        validate_graph(molecule)?;
        self.commit_definition_and_instance(molecule.clone())
            .map(|(_, instance)| instance)
    }

    pub fn build(mut self) -> Result<Topology, TopologyBuildError> {
        if self.instances.is_empty() {
            return Err(TopologyBuildError::NoMoleculeInstances);
        }
        for instance in &self.instances {
            if self.definitions.get(instance.definition.index()).is_none() {
                return Err(TopologyBuildError::InvalidMoleculeDefinitionId(
                    instance.definition,
                ));
            }
        }
        let mut referenced_definitions = vec![false; self.definitions.len()];
        for instance in &self.instances {
            referenced_definitions[instance.definition.index()] = true;
        }
        if let Some(index) = referenced_definitions
            .iter()
            .position(|referenced| !referenced)
        {
            return Err(TopologyBuildError::UnusedMoleculeDefinition(
                self.definitions[index].id,
            ));
        }
        for definition in &self.definitions {
            validate_graph(definition.molecule())?;
        }

        let atom_count = self.instances.iter().try_fold(0usize, |count, instance| {
            count
                .checked_add(
                    self.definitions[instance.definition.index()]
                        .molecule()
                        .atom_count(),
                )
                .ok_or(TopologyBuildError::IdentifierCapacityExceeded(
                    TopologyIdKind::Atom,
                ))
        })?;
        checked_future_len(0, atom_count, TopologyIdKind::Atom)?;
        let bond_count = self.instances.iter().try_fold(0usize, |count, instance| {
            count
                .checked_add(
                    self.definitions[instance.definition.index()]
                        .molecule()
                        .bond_count(),
                )
                .ok_or(TopologyBuildError::IdentifierCapacityExceeded(
                    TopologyIdKind::Bond,
                ))
        })?;
        checked_future_len(0, bond_count, TopologyIdKind::Bond)?;

        let mut instance_atoms = Vec::new();
        let mut instance_bonds = Vec::new();
        let mut atom_indices = BTreeMap::new();
        let mut bond_indices = BTreeMap::new();
        instance_atoms
            .try_reserve_exact(atom_count)
            .map_err(|_| TopologyBuildError::IdentifierCapacityExceeded(TopologyIdKind::Atom))?;
        instance_bonds
            .try_reserve_exact(bond_count)
            .map_err(|_| TopologyBuildError::IdentifierCapacityExceeded(TopologyIdKind::Bond))?;

        for instance in &self.instances {
            let molecule = self.definitions[instance.definition.index()].molecule();
            for atom in molecule.atom_ids() {
                let qualified = instance.qualify_atom(atom);
                let index =
                    checked_id::<TopologyAtomIndex>(instance_atoms.len(), TopologyIdKind::Atom)?;
                atom_indices.insert(qualified, index);
                instance_atoms.push(qualified);
            }
            for bond in molecule.bond_ids() {
                let qualified = instance.qualify_bond(bond);
                let index =
                    checked_id::<TopologyBondIndex>(instance_bonds.len(), TopologyIdKind::Bond)?;
                bond_indices.insert(qualified, index);
                instance_bonds.push(qualified);
            }
        }

        validate_hierarchy(&self.hierarchy, &atom_indices)
            .map_err(TopologyBuildError::InvalidHierarchy)?;

        self.properties.resize_domains(
            self.instances.len(),
            atom_count,
            bond_count,
            self.hierarchy.chains().count(),
            self.hierarchy.residues().count(),
            self.hierarchy.atom_sites().count(),
        );

        Ok(Topology {
            definitions: self.definitions,
            instances: self.instances,
            instance_atoms,
            instance_bonds,
            atom_indices,
            bond_indices,
            hierarchy: self.hierarchy,
            properties: self.properties,
        })
    }

    fn sync_property_dimensions(&mut self) {
        let atom_count = self
            .instances
            .iter()
            .map(|instance| {
                self.definitions[instance.definition.index()]
                    .molecule()
                    .atom_count()
            })
            .sum();
        let bond_count = self
            .instances
            .iter()
            .map(|instance| {
                self.definitions[instance.definition.index()]
                    .molecule()
                    .bond_count()
            })
            .sum();
        self.properties.resize_domains(
            self.instances.len(),
            atom_count,
            bond_count,
            self.hierarchy.chains().count(),
            self.hierarchy.residues().count(),
            self.hierarchy.atom_sites().count(),
        );
    }

    fn commit_definition(
        &mut self,
        molecule: Molecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        self.reserve_definitions(1)?;
        let id = checked_id::<MoleculeDefinitionId>(
            self.definitions.len(),
            TopologyIdKind::MoleculeDefinition,
        )?;
        self.definitions.push(MoleculeDefinition { id, molecule });
        Ok(id)
    }

    fn commit_definition_and_instance(
        &mut self,
        molecule: Molecule,
    ) -> Result<(MoleculeDefinitionId, MoleculeInstanceId), TopologyBuildError> {
        self.reserve_definitions(1)?;
        self.reserve_instances(1)?;
        let definition = checked_id::<MoleculeDefinitionId>(
            self.definitions.len(),
            TopologyIdKind::MoleculeDefinition,
        )?;
        let instance = checked_id::<MoleculeInstanceId>(
            self.instances.len(),
            TopologyIdKind::MoleculeInstance,
        )?;
        self.definitions.push(MoleculeDefinition {
            id: definition,
            molecule,
        });
        self.instances.push(MoleculeInstance {
            id: instance,
            definition,
        });
        Ok((definition, instance))
    }
}

trait FromRawId {
    fn from_raw(raw: u32) -> Self;
}

impl FromRawId for MoleculeDefinitionId {
    fn from_raw(raw: u32) -> Self {
        Self::new(raw)
    }
}

impl FromRawId for MoleculeInstanceId {
    fn from_raw(raw: u32) -> Self {
        Self::new(raw)
    }
}

impl FromRawId for TopologyAtomIndex {
    fn from_raw(raw: u32) -> Self {
        Self::new(raw)
    }
}

impl FromRawId for TopologyBondIndex {
    fn from_raw(raw: u32) -> Self {
        Self::new(raw)
    }
}

fn checked_id<T: FromRawId>(length: usize, kind: TopologyIdKind) -> Result<T, TopologyBuildError> {
    let raw = crate::core::checked_raw_id(length)
        .map_err(|_| TopologyBuildError::IdentifierCapacityExceeded(kind))?;
    Ok(T::from_raw(raw))
}

fn checked_future_len(
    current: usize,
    additional: usize,
    kind: TopologyIdKind,
) -> Result<(), TopologyBuildError> {
    crate::core::checked_fixed_id_collection_len(current, additional)
        .map_err(|_| TopologyBuildError::IdentifierCapacityExceeded(kind))
}

fn validate_graph(graph: &Molecule) -> Result<(), TopologyBuildError> {
    validate_nonempty_graph(graph)?;
    graph
        .validate_connected()
        .map_err(TopologyBuildError::DisconnectedMoleculeDefinition)
}

fn validate_nonempty_graph(graph: &Molecule) -> Result<(), TopologyBuildError> {
    if graph.atom_count() == 0 {
        return Err(TopologyBuildError::EmptyMoleculeDefinition);
    }
    Ok(())
}

/// A hierarchy reference or reverse lookup that cannot be published in a topology.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologyHierarchyError {
    InvalidChainIdentifier {
        slot: usize,
        id: ChainId,
    },
    InvalidResidueIdentifier {
        slot: usize,
        id: ResidueId,
    },
    InvalidAtomSiteIdentifier {
        slot: usize,
        id: AtomSiteId,
    },
    InvalidChainResidue {
        chain: ChainId,
        residue: ResidueId,
    },
    InvalidResidueChain {
        residue: ResidueId,
        chain: ChainId,
    },
    InconsistentChainResidue {
        chain: ChainId,
        residue: ResidueId,
    },
    InvalidResidueAtomSite {
        residue: ResidueId,
        site: AtomSiteId,
    },
    InvalidAtomSiteResidue {
        site: AtomSiteId,
        residue: ResidueId,
    },
    InconsistentResidueAtomSite {
        residue: ResidueId,
        site: AtomSiteId,
    },
    InvalidAtomSiteAtom {
        site: AtomSiteId,
        atom: InstanceAtomId,
    },
    InconsistentAtomLookup {
        site: AtomSiteId,
        atom: InstanceAtomId,
    },
}

impl fmt::Display for TopologyHierarchyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidChainIdentifier { slot, id } => write!(
                formatter,
                "hierarchy chain slot {slot} stores non-matching identifier {id}"
            ),
            Self::InvalidResidueIdentifier { slot, id } => write!(
                formatter,
                "hierarchy residue slot {slot} stores non-matching identifier {id}"
            ),
            Self::InvalidAtomSiteIdentifier { slot, id } => write!(
                formatter,
                "hierarchy atom-site slot {slot} stores non-matching identifier {id}"
            ),
            Self::InvalidChainResidue { chain, residue } => write!(
                formatter,
                "hierarchy chain {chain} references missing residue {residue}"
            ),
            Self::InvalidResidueChain { residue, chain } => write!(
                formatter,
                "hierarchy residue {residue} references missing chain {chain}"
            ),
            Self::InconsistentChainResidue { chain, residue } => write!(
                formatter,
                "hierarchy chain {chain} and residue {residue} do not reference each other"
            ),
            Self::InvalidResidueAtomSite { residue, site } => write!(
                formatter,
                "hierarchy residue {residue} references missing atom site {site}"
            ),
            Self::InvalidAtomSiteResidue { site, residue } => write!(
                formatter,
                "hierarchy atom site {site} references missing residue {residue}"
            ),
            Self::InconsistentResidueAtomSite { residue, site } => write!(
                formatter,
                "hierarchy residue {residue} and atom site {site} do not reference each other"
            ),
            Self::InvalidAtomSiteAtom { site, atom } => write!(
                formatter,
                "hierarchy atom site {site} references missing topology atom {atom}"
            ),
            Self::InconsistentAtomLookup { site, atom } => write!(
                formatter,
                "hierarchy atom lookup for {atom} is inconsistent with atom site {site}"
            ),
        }
    }
}

impl std::error::Error for TopologyHierarchyError {}

fn validate_hierarchy(
    hierarchy: &Hierarchy,
    atom_indices: &BTreeMap<InstanceAtomId, TopologyAtomIndex>,
) -> Result<(), TopologyHierarchyError> {
    for (slot, (chain_id, chain)) in hierarchy.chains().enumerate() {
        if chain_id.index() != slot {
            return Err(TopologyHierarchyError::InvalidChainIdentifier { slot, id: chain_id });
        }
        for residue_id in chain.residues() {
            let residue = hierarchy.residue(*residue_id).map_err(|_| {
                TopologyHierarchyError::InvalidChainResidue {
                    chain: chain_id,
                    residue: *residue_id,
                }
            })?;
            if residue.chain() != chain_id
                || chain
                    .residues()
                    .iter()
                    .filter(|candidate| **candidate == *residue_id)
                    .count()
                    != 1
            {
                return Err(TopologyHierarchyError::InconsistentChainResidue {
                    chain: chain_id,
                    residue: *residue_id,
                });
            }
        }
    }
    for (slot, (residue_id, residue)) in hierarchy.residues().enumerate() {
        if residue_id.index() != slot {
            return Err(TopologyHierarchyError::InvalidResidueIdentifier {
                slot,
                id: residue_id,
            });
        }
        let chain = hierarchy.chain(residue.chain()).map_err(|_| {
            TopologyHierarchyError::InvalidResidueChain {
                residue: residue_id,
                chain: residue.chain(),
            }
        })?;
        if chain
            .residues()
            .iter()
            .filter(|candidate| **candidate == residue_id)
            .count()
            != 1
        {
            return Err(TopologyHierarchyError::InconsistentChainResidue {
                chain: residue.chain(),
                residue: residue_id,
            });
        }
        for site_id in residue.atom_sites() {
            let site = hierarchy.atom_site(*site_id).map_err(|_| {
                TopologyHierarchyError::InvalidResidueAtomSite {
                    residue: residue_id,
                    site: *site_id,
                }
            })?;
            if site.residue() != residue_id
                || residue
                    .atom_sites()
                    .iter()
                    .filter(|candidate| **candidate == *site_id)
                    .count()
                    != 1
            {
                return Err(TopologyHierarchyError::InconsistentResidueAtomSite {
                    residue: residue_id,
                    site: *site_id,
                });
            }
        }
    }
    for (slot, (site_id, site)) in hierarchy.atom_sites().enumerate() {
        if site_id.index() != slot {
            return Err(TopologyHierarchyError::InvalidAtomSiteIdentifier { slot, id: site_id });
        }
        let residue = hierarchy.residue(site.residue()).map_err(|_| {
            TopologyHierarchyError::InvalidAtomSiteResidue {
                site: site_id,
                residue: site.residue(),
            }
        })?;
        if residue
            .atom_sites()
            .iter()
            .filter(|candidate| **candidate == site_id)
            .count()
            != 1
        {
            return Err(TopologyHierarchyError::InconsistentResidueAtomSite {
                residue: site.residue(),
                site: site_id,
            });
        }
        if !atom_indices.contains_key(&site.atom()) {
            return Err(TopologyHierarchyError::InvalidAtomSiteAtom {
                site: site_id,
                atom: site.atom(),
            });
        }
        if hierarchy
            .atom_site_for_atom(site.atom())
            .is_none_or(|mapped| mapped.id() != site_id)
        {
            return Err(TopologyHierarchyError::InconsistentAtomLookup {
                site: site_id,
                atom: site.atom(),
            });
        }
    }
    for (atom, site_id) in hierarchy.atom_lookup_entries() {
        if !hierarchy
            .atom_site(site_id)
            .is_ok_and(|site| site.atom() == atom)
        {
            return Err(TopologyHierarchyError::InconsistentAtomLookup {
                site: site_id,
                atom,
            });
        }
    }
    Ok(())
}

/// Fixed-width identifier spaces owned by [`Topology`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopologyIdKind {
    /// Reusable molecule definitions.
    MoleculeDefinition,
    /// Explicit molecule instances.
    MoleculeInstance,
    /// Authoritative dense atom indices.
    Atom,
    /// Authoritative dense bond indices.
    Bond,
}

impl fmt::Display for TopologyIdKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::MoleculeDefinition => "molecule definition",
            Self::MoleculeInstance => "molecule instance",
            Self::Atom => "topology atom",
            Self::Bond => "topology bond",
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologyBuildError {
    NoMoleculeInstances,
    EmptyMoleculeDefinition,
    DisconnectedMoleculeDefinition(MoleculeConnectivityError),
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    /// A staged definition was not instantiated and cannot be published.
    UnusedMoleculeDefinition(MoleculeDefinitionId),
    /// The staged system hierarchy is inconsistent with itself or the topology.
    InvalidHierarchy(TopologyHierarchyError),
    /// A topology collection exceeded the fixed-width identifier space for `kind`.
    IdentifierCapacityExceeded(TopologyIdKind),
}

impl fmt::Display for TopologyBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoMoleculeInstances => {
                formatter.write_str("topology must contain at least one molecule instance")
            }
            Self::EmptyMoleculeDefinition => {
                formatter.write_str("molecule definition must contain at least one atom")
            }
            Self::DisconnectedMoleculeDefinition(error) => {
                write!(formatter, "molecule definition is disconnected: {error}")
            }
            Self::InvalidMoleculeDefinitionId(id) => {
                write!(formatter, "invalid molecule definition: {id}")
            }
            Self::UnusedMoleculeDefinition(id) => {
                write!(
                    formatter,
                    "molecule definition {id} is not referenced by any instance"
                )
            }
            Self::InvalidHierarchy(error) => {
                write!(formatter, "invalid topology hierarchy: {error}")
            }
            Self::IdentifierCapacityExceeded(kind) => {
                write!(formatter, "{kind} identifier capacity exceeded")
            }
        }
    }
}

impl std::error::Error for TopologyBuildError {}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologyError {
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    InvalidAtomId(InstanceAtomId),
    InvalidBondId(InstanceBondId),
    InvalidChainId(ChainId),
    InvalidResidueId(ResidueId),
    InvalidAtomSiteId(AtomSiteId),
    InvalidAtomIndex(TopologyAtomIndex),
    InvalidBondIndex(TopologyBondIndex),
    Property(Box<PropertyError>),
}

impl fmt::Display for TopologyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMoleculeDefinitionId(id) => {
                write!(formatter, "invalid molecule definition: {id}")
            }
            Self::InvalidMoleculeInstanceId(id) => {
                write!(formatter, "invalid molecule instance: {id}")
            }
            Self::InvalidAtomId(id) => write!(formatter, "invalid topology atom: {id}"),
            Self::InvalidBondId(id) => write!(formatter, "invalid topology bond: {id}"),
            Self::InvalidChainId(id) => write!(formatter, "invalid topology chain: {id}"),
            Self::InvalidResidueId(id) => write!(formatter, "invalid topology residue: {id}"),
            Self::InvalidAtomSiteId(id) => write!(formatter, "invalid topology atom site: {id}"),
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid {index}"),
            Self::InvalidBondIndex(index) => write!(formatter, "invalid {index}"),
            Self::Property(error) => write!(formatter, "invalid topology property: {error}"),
        }
    }
}

impl std::error::Error for TopologyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Property(error) => Some(error),
            _ => None,
        }
    }
}

impl Eq for TopologyError {}

/// A topology-bound, sorted, unique dense atom selection.
#[derive(Debug, Clone)]
pub struct AtomSelection {
    topology: Arc<Topology>,
    indices: Vec<TopologyAtomIndex>,
}

impl PartialEq for AtomSelection {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.indices == other.indices
    }
}

impl Eq for AtomSelection {}

impl AtomSelection {
    pub fn from_atoms(
        topology: &Arc<Topology>,
        atoms: impl IntoIterator<Item = InstanceAtomId>,
    ) -> Result<Self, SelectionError> {
        let mut indices = atoms
            .into_iter()
            .map(|atom| {
                topology
                    .atom_index(atom)
                    .ok_or(SelectionError::InvalidAtomId(atom))
            })
            .collect::<Result<Vec<_>, _>>()?;
        indices.sort_unstable();
        indices.dedup();
        Ok(Self {
            topology: Arc::clone(topology),
            indices,
        })
    }

    pub fn from_indices(
        topology: &Arc<Topology>,
        indices: impl IntoIterator<Item = TopologyAtomIndex>,
    ) -> Result<Self, SelectionError> {
        let atoms = indices
            .into_iter()
            .map(|index| {
                topology
                    .atom_id(index)
                    .ok_or(SelectionError::InvalidAtomIndex(index))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_atoms(topology, atoms)
    }

    pub fn for_instances(
        topology: &Arc<Topology>,
        instances: impl IntoIterator<Item = MoleculeInstanceId>,
    ) -> Result<Self, SelectionError> {
        let instances = instances.into_iter().collect::<BTreeSet<_>>();
        for instance in &instances {
            topology
                .instance(*instance)
                .map_err(|_| SelectionError::InvalidMoleculeInstanceId(*instance))?;
        }
        Self::from_atoms(
            topology,
            topology
                .atom_ids()
                .iter()
                .copied()
                .filter(|atom| instances.contains(&atom.molecule())),
        )
    }

    pub fn for_definitions(
        topology: &Arc<Topology>,
        definitions: impl IntoIterator<Item = MoleculeDefinitionId>,
    ) -> Result<Self, SelectionError> {
        let definitions = definitions.into_iter().collect::<BTreeSet<_>>();
        for definition in &definitions {
            topology
                .definition(*definition)
                .map_err(|_| SelectionError::InvalidMoleculeDefinitionId(*definition))?;
        }
        let instances = topology
            .instances()
            .filter(|(_, instance)| definitions.contains(&instance.definition()))
            .map(|(id, _)| id);
        Self::for_instances(topology, instances)
    }

    pub fn for_elements(
        topology: &Arc<Topology>,
        elements: impl IntoIterator<Item = Element>,
    ) -> Result<Self, SelectionError> {
        let elements = elements.into_iter().collect::<BTreeSet<_>>();
        Self::from_atoms(
            topology,
            topology
                .atoms()
                .filter(|(_, atom)| elements.contains(&atom.element))
                .map(|(id, _)| id),
        )
    }

    pub fn for_atom_sites(
        topology: &Arc<Topology>,
        atom_sites: impl IntoIterator<Item = AtomSiteId>,
    ) -> Result<Self, SelectionError> {
        let atom_sites = atom_sites.into_iter().collect::<BTreeSet<_>>();
        let atoms = atom_sites
            .into_iter()
            .map(|site| {
                topology
                    .atom_for_site(site)
                    .map_err(|_| SelectionError::InvalidAtomSiteId(site))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_atoms(topology, atoms)
    }

    pub fn for_residues(
        topology: &Arc<Topology>,
        residues: impl IntoIterator<Item = ResidueId>,
    ) -> Result<Self, SelectionError> {
        let residues = residues.into_iter().collect::<BTreeSet<_>>();
        for residue in &residues {
            topology
                .residue(*residue)
                .map_err(|_| SelectionError::InvalidResidueId(*residue))?;
        }
        Self::for_atom_sites(
            topology,
            topology
                .atom_sites()
                .filter(|site| residues.contains(&site.residue().id()))
                .map(AtomSiteView::id),
        )
    }

    pub fn for_chains(
        topology: &Arc<Topology>,
        chains: impl IntoIterator<Item = ChainId>,
    ) -> Result<Self, SelectionError> {
        let chains = chains.into_iter().collect::<BTreeSet<_>>();
        for chain in &chains {
            topology
                .chain(*chain)
                .map_err(|_| SelectionError::InvalidChainId(*chain))?;
        }
        Self::for_residues(
            topology,
            topology
                .residues()
                .filter(|residue| chains.contains(&residue.chain().id()))
                .map(ResidueView::id),
        )
    }

    pub fn for_chain_label(topology: &Arc<Topology>, label: &str) -> Result<Self, SelectionError> {
        Self::for_chains(
            topology,
            topology
                .chains()
                .filter(|chain| chain.label_id() == label)
                .map(ChainView::id),
        )
    }

    pub fn from_query_matches(
        topology: &Arc<Topology>,
        instance: MoleculeInstanceId,
        matches: &[QueryMatch],
    ) -> Result<Self, SelectionError> {
        topology
            .instance(instance)
            .map_err(|_| SelectionError::InvalidMoleculeInstanceId(instance))?;
        Self::from_atoms(
            topology,
            matches.iter().flat_map(|query_match| {
                query_match
                    .atoms()
                    .iter()
                    .copied()
                    .map(move |atom| InstanceAtomId::new(instance, atom))
            }),
        )
    }

    pub fn ensure_compatible(&self, topology: &Arc<Topology>) -> Result<(), SelectionError> {
        if !Arc::ptr_eq(&self.topology, topology) {
            return Err(SelectionError::TopologyMismatch);
        }
        Ok(())
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub(crate) fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn indices(&self) -> &[TopologyAtomIndex] {
        &self.indices
    }

    pub fn semantic_ids(
        &self,
        topology: &Arc<Topology>,
    ) -> Result<Vec<InstanceAtomId>, SelectionError> {
        self.ensure_compatible(topology)?;
        Ok(self
            .indices
            .iter()
            .map(|index| {
                topology
                    .atom_id(*index)
                    .expect("selection contains validated dense indices")
            })
            .collect())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SelectionError {
    TopologyMismatch,
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    InvalidAtomId(InstanceAtomId),
    InvalidChainId(ChainId),
    InvalidResidueId(ResidueId),
    InvalidAtomSiteId(AtomSiteId),
    InvalidAtomIndex(TopologyAtomIndex),
}

impl fmt::Display for SelectionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("atom selection belongs to a different topology")
            }
            Self::InvalidMoleculeDefinitionId(id) => {
                write!(formatter, "invalid selected molecule definition: {id}")
            }
            Self::InvalidMoleculeInstanceId(id) => {
                write!(formatter, "invalid selected molecule instance: {id}")
            }
            Self::InvalidAtomId(id) => write!(formatter, "invalid selected atom: {id}"),
            Self::InvalidChainId(id) => write!(formatter, "invalid selected chain: {id}"),
            Self::InvalidResidueId(id) => write!(formatter, "invalid selected residue: {id}"),
            Self::InvalidAtomSiteId(id) => write!(formatter, "invalid selected atom site: {id}"),
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid selected {index}"),
        }
    }
}

impl std::error::Error for SelectionError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::BondOrder;
    use crate::properties::{PropertyColumn, PropertyKey, PropertyValue};
    use crate::query;
    use crate::substructure;
    use crate::topology::AtomSiteMetadata;

    fn perceived_molecule(smiles: &str) -> Molecule {
        let mut molecule = crate::tests::read_smiles(smiles).expect("SMILES should parse");
        molecule.perceive().expect("molecule should perceive");
        molecule
    }

    #[test]
    fn hierarchy_errors_have_diagnostic_display_messages() {
        let error = TopologyHierarchyError::InconsistentResidueAtomSite {
            residue: ResidueId::new(2),
            site: AtomSiteId::new(7),
        };
        let message = error.to_string();
        assert!(message.contains("residue2"));
        assert!(message.contains("atom-site7"));
        assert!(message.contains("do not reference each other"));
    }

    fn tombstoned_molecule() -> (Molecule, AtomId, AtomId, BondId) {
        let mut graph = crate::core::MoleculeEditor::new();
        let carbon = graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");
        let tombstone = graph
            .add_atom(Atom::new(Element::from_symbol("H").unwrap()))
            .expect("atom identifier capacity");
        let oxygen = graph
            .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
            .expect("atom identifier capacity");
        graph.delete_atom(tombstone).unwrap();
        let deleted_bond = graph.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
        graph.delete_bond(deleted_bond).unwrap();
        let bond = graph.add_bond(carbon, oxygen, BondOrder::Double).unwrap();
        (graph.finish().unwrap(), carbon, oxygen, bond)
    }

    fn topology_with_reused_definition() -> (Arc<Topology>, AtomId, AtomId, BondId) {
        let (molecule, carbon, oxygen, bond) = tombstoned_molecule();
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(&molecule).unwrap();
        builder.add_instance(definition).unwrap();
        builder.add_instance(definition).unwrap();
        (Arc::new(builder.build().unwrap()), carbon, oxygen, bond)
    }

    fn topology_with_two_distinct_definitions(reverse: bool) -> Arc<Topology> {
        let carbon = perceived_molecule("C");
        let carbon_oxygen = perceived_molecule("CO");
        let molecules = if reverse {
            [&carbon_oxygen, &carbon]
        } else {
            [&carbon, &carbon_oxygen]
        };
        let mut builder = TopologyBuilder::new();
        for molecule in molecules {
            let definition = builder.add_molecule_definition(molecule).unwrap();
            builder.add_instance(definition).unwrap();
        }
        Arc::new(builder.build().unwrap())
    }

    #[test]
    fn topology_reuses_definitions_and_preserves_qualified_dense_order() {
        let (topology, carbon, oxygen, bond) = topology_with_reused_definition();
        assert_eq!(topology.definition_count(), 1);
        assert_eq!(topology.instance_count(), 2);
        assert_eq!(topology.atom_count(), 4);
        assert_eq!(topology.bond_count(), 2);

        let first = MoleculeInstanceId::new(0);
        let second = MoleculeInstanceId::new(1);
        assert!(std::ptr::eq(
            topology.definition_for_instance(first).unwrap(),
            topology.definition_for_instance(second).unwrap()
        ));
        assert_eq!(
            topology.atom_ids(),
            &[
                InstanceAtomId::new(first, carbon),
                InstanceAtomId::new(first, oxygen),
                InstanceAtomId::new(second, carbon),
                InstanceAtomId::new(second, oxygen),
            ]
        );
        assert_eq!(
            topology.bond_ids(),
            &[
                InstanceBondId::new(first, bond),
                InstanceBondId::new(second, bond),
            ]
        );
        for (raw, atom) in topology.atom_ids().iter().copied().enumerate() {
            let index = topology.atom_index(atom).unwrap();
            assert_eq!(index.index(), raw);
            assert_eq!(topology.atom_id(index), Some(atom));
        }
        for (raw, bond) in topology.bond_ids().iter().copied().enumerate() {
            let index = topology.bond_index(bond).unwrap();
            assert_eq!(index.index(), raw);
            assert_eq!(topology.bond_id(index), Some(bond));
        }
    }

    #[test]
    fn molecule_views_are_instance_qualified_and_share_definition_state() {
        let (topology, carbon, oxygen, bond) = topology_with_reused_definition();
        let molecules = topology.molecules().collect::<Vec<_>>();
        assert_eq!(molecules.len(), 2);
        assert_eq!(molecules[0].id(), MoleculeInstanceId::new(0));
        assert_eq!(molecules[1].id(), MoleculeInstanceId::new(1));
        assert_eq!(molecules[0].definition_id(), molecules[1].definition_id());
        assert!(std::ptr::eq(
            molecules[0].molecule(),
            molecules[1].molecule()
        ));
        assert_eq!(
            molecules[0].atoms().map(|(id, _)| id).collect::<Vec<_>>(),
            vec![
                InstanceAtomId::new(molecules[0].id(), carbon),
                InstanceAtomId::new(molecules[0].id(), oxygen),
            ]
        );
        assert_eq!(
            molecules[1].bonds().map(|(id, _)| id).collect::<Vec<_>>(),
            vec![InstanceBondId::new(molecules[1].id(), bond)]
        );
        assert_eq!(
            topology.molecule(molecules[1].id()).unwrap().id(),
            molecules[1].id()
        );
    }

    #[test]
    fn builder_rejects_empty_topologies_and_unused_definitions() {
        assert!(matches!(
            TopologyBuilder::new().build(),
            Err(TopologyBuildError::NoMoleculeInstances)
        ));

        let molecule = perceived_molecule("O");
        let mut builder = TopologyBuilder::new();
        let used = builder.add_molecule_definition(&molecule).unwrap();
        builder.add_instance(used).unwrap();
        let unused = builder.add_molecule_definition(&molecule).unwrap();
        assert!(matches!(
            builder.build(),
            Err(TopologyBuildError::UnusedMoleculeDefinition(id)) if id == unused
        ));
    }

    #[test]
    fn builder_add_molecule_is_the_concise_single_instance_path() {
        let molecule = perceived_molecule("CO");
        let mut builder = TopologyBuilder::new();
        let instance = builder.add_molecule(&molecule).unwrap();
        let topology = builder.build().unwrap();
        assert_eq!(instance, MoleculeInstanceId::new(0));
        assert_eq!(topology.definition_count(), 1);
        assert_eq!(topology.instance_count(), 1);
        assert_eq!(topology.molecule(instance).unwrap().molecule(), &molecule);
    }

    #[test]
    fn topology_directly_owns_its_layout_collections() {
        let (topology, ..) = topology_with_reused_definition();
        assert_eq!(topology.definitions.len(), 1);
        assert_eq!(topology.instances.len(), 2);
        assert_eq!(topology.instance_atoms.len(), 4);
        assert_eq!(topology.instance_bonds.len(), 2);
        assert_eq!(topology.atom_indices.len(), 4);
        assert_eq!(topology.bond_indices.len(), 2);
        for &atom in &topology.instance_atoms {
            assert_eq!(
                topology.instance_atoms[topology.atom_indices[&atom].index()],
                atom
            );
        }
        for &bond in &topology.instance_bonds {
            assert_eq!(
                topology.instance_bonds[topology.bond_indices[&bond].index()],
                bond
            );
        }

        let debug = format!("{topology:?}");
        assert!(debug.contains("instance_atoms"));
        assert!(debug.contains("instance_bonds"));
    }

    #[test]
    fn shared_allocation_is_exact_while_layout_can_match() {
        let (topology, ..) = topology_with_reused_definition();
        let clone = Arc::clone(&topology);
        let (independent, ..) = topology_with_reused_definition();

        assert!(Arc::ptr_eq(&topology, &clone));
        assert!(topology.same_layout(&clone));
        assert!(!Arc::ptr_eq(&topology, &independent));
        assert!(topology.same_layout(&independent));
        assert_eq!(clone.atom_ids(), topology.atom_ids());
    }

    #[test]
    fn layout_equality_does_not_reorder_definitions_instances_or_dense_state() {
        let forward = topology_with_two_distinct_definitions(false);
        let reverse = topology_with_two_distinct_definitions(true);

        assert!(!Arc::ptr_eq(&forward, &reverse));
        assert!(!forward.same_layout(&reverse));
    }

    #[test]
    fn topology_properties_cover_every_domain_and_do_not_change_layout_identity() {
        let (molecule, carbon, oxygen, _) = tombstoned_molecule();
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(&molecule).unwrap();
        let instance = builder.add_instance(definition).unwrap();
        let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, "LIG", Some(1), None, None)
            .unwrap();
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, carbon),
                AtomSiteMetadata {
                    label_atom_id: Some("C1".into()),
                    ..AtomSiteMetadata::default()
                },
            )
            .unwrap();

        let owner_key = PropertyKey::new("source").unwrap();
        let value_key = PropertyKey::new("tag").unwrap();
        builder
            .insert_property(owner_key.clone(), PropertyValue::String("test".into()))
            .unwrap();
        fn insert_tag(table: &mut crate::properties::PropertyTable, key: &PropertyKey) {
            table
                .insert(key.clone(), PropertyColumn::Int(vec![Some(1); table.len()]))
                .unwrap();
        }
        insert_tag(builder.molecule_instance_properties_mut(), &value_key);
        insert_tag(builder.atom_properties_mut(), &value_key);
        insert_tag(builder.bond_properties_mut(), &value_key);
        insert_tag(builder.chain_properties_mut(), &value_key);
        insert_tag(builder.residue_properties_mut(), &value_key);
        insert_tag(builder.atom_site_properties_mut(), &value_key);
        let enriched = builder.build().unwrap();

        assert_eq!(
            enriched.properties().get(&owner_key),
            Some(&PropertyValue::String("test".into()))
        );
        assert_eq!(enriched.molecule_instance_properties().len(), 1);
        assert_eq!(enriched.atom_properties().len(), 2);
        assert_eq!(enriched.bond_properties().len(), 1);
        assert_eq!(enriched.chain_properties().len(), 1);
        assert_eq!(enriched.residue_properties().len(), 1);
        assert_eq!(enriched.atom_site_properties().len(), 1);
        assert_eq!(
            enriched
                .molecule_instance_property(instance, &value_key)
                .unwrap(),
            Some(PropertyValue::Int(1))
        );
        assert_eq!(
            enriched
                .molecule(instance)
                .unwrap()
                .property(&value_key)
                .unwrap(),
            Some(PropertyValue::Int(1))
        );
        assert_eq!(
            enriched
                .atom_property(InstanceAtomId::new(instance, oxygen), &value_key)
                .unwrap(),
            Some(PropertyValue::Int(1))
        );

        let mut plain_builder = TopologyBuilder::new();
        let definition = plain_builder.add_molecule_definition(&molecule).unwrap();
        let instance = plain_builder.add_instance(definition).unwrap();
        let chain = plain_builder.hierarchy_mut().add_chain("A", None).unwrap();
        let residue = plain_builder
            .hierarchy_mut()
            .add_residue(chain, "LIG", Some(1), None, None)
            .unwrap();
        plain_builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, carbon),
                AtomSiteMetadata {
                    label_atom_id: Some("C1".into()),
                    ..AtomSiteMetadata::default()
                },
            )
            .unwrap();
        let plain = plain_builder.build().unwrap();
        assert!(enriched.same_layout(&plain));
    }

    #[test]
    fn builder_is_transactional_and_does_not_intern_equal_definitions() {
        let (molecule, ..) = tombstoned_molecule();
        let mut builder = TopologyBuilder::new();
        let first = builder.add_molecule_definition(&molecule).unwrap();
        assert_eq!(
            builder.add_instance(MoleculeDefinitionId::new(99)),
            Err(TopologyBuildError::InvalidMoleculeDefinitionId(
                MoleculeDefinitionId::new(99)
            ))
        );
        assert!(builder.instances.is_empty());
        builder.add_instance(first).unwrap();
        let second = builder.add_molecule_definition(&molecule).unwrap();
        builder.add_instance(second).unwrap();
        let topology = Arc::new(builder.build().unwrap());
        assert_eq!(topology.definition_count(), 2);
        assert_eq!(topology.instance_count(), 2);
        assert_eq!(
            checked_future_len(usize::MAX, 1, TopologyIdKind::Atom),
            Err(TopologyBuildError::IdentifierCapacityExceeded(
                TopologyIdKind::Atom
            ))
        );
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn every_topology_identifier_space_checks_its_boundary() {
        let max_slot = usize::try_from(u64::from(u32::MAX)).expect("64-bit usize");
        let first_unsupported_slot =
            usize::try_from(u64::from(u32::MAX) + 1).expect("64-bit usize");

        assert_eq!(
            checked_id::<MoleculeDefinitionId>(max_slot, TopologyIdKind::MoleculeDefinition,),
            Ok(MoleculeDefinitionId::new(u32::MAX))
        );
        assert_eq!(
            checked_id::<MoleculeDefinitionId>(
                first_unsupported_slot,
                TopologyIdKind::MoleculeDefinition,
            ),
            Err(TopologyBuildError::IdentifierCapacityExceeded(
                TopologyIdKind::MoleculeDefinition
            ))
        );
        assert_eq!(
            checked_id::<MoleculeInstanceId>(
                first_unsupported_slot,
                TopologyIdKind::MoleculeInstance,
            ),
            Err(TopologyBuildError::IdentifierCapacityExceeded(
                TopologyIdKind::MoleculeInstance
            ))
        );
        assert_eq!(
            checked_id::<TopologyAtomIndex>(first_unsupported_slot, TopologyIdKind::Atom,),
            Err(TopologyBuildError::IdentifierCapacityExceeded(
                TopologyIdKind::Atom
            ))
        );
        assert_eq!(
            checked_id::<TopologyBondIndex>(first_unsupported_slot, TopologyIdKind::Bond,),
            Err(TopologyBuildError::IdentifierCapacityExceeded(
                TopologyIdKind::Bond
            ))
        );
        assert_eq!(
            checked_future_len(first_unsupported_slot, 1, TopologyIdKind::Atom),
            Err(TopologyBuildError::IdentifierCapacityExceeded(
                TopologyIdKind::Atom
            ))
        );
    }

    #[test]
    fn selections_distinguish_instances_elements_and_queries() {
        let ethane = perceived_molecule("CC");
        let water = perceived_molecule("O");
        let mut builder = TopologyBuilder::new();
        let ethane_definition = builder.add_molecule_definition(&ethane).unwrap();
        let water_definition = builder.add_molecule_definition(&water).unwrap();
        let ethane_instance = builder.add_instance(ethane_definition).unwrap();
        let water_instance = builder.add_instance(water_definition).unwrap();
        let topology = Arc::new(builder.build().unwrap());

        let selected =
            AtomSelection::for_instances(&topology, [ethane_instance, water_instance]).unwrap();
        assert_eq!(selected.indices().len(), 3);
        let ethane_selection = AtomSelection::for_instances(&topology, [ethane_instance]).unwrap();
        let water_selection = AtomSelection::for_instances(&topology, [water_instance]).unwrap();
        assert_eq!(ethane_selection.indices().len(), 2);
        assert_eq!(water_selection.indices().len(), 1);
        let oxygen =
            AtomSelection::for_elements(&topology, [Element::from_symbol("O").unwrap()]).unwrap();
        assert_eq!(oxygen.indices().len(), 1);

        let query = query::parse_smarts("O").unwrap();
        let matches = substructure::find_substructure_matches(&water, &query).unwrap();
        let from_query =
            AtomSelection::from_query_matches(&topology, water_instance, &matches).unwrap();
        assert_eq!(
            from_query.semantic_ids(&topology).unwrap(),
            oxygen.semantic_ids(&topology).unwrap()
        );

        let mut independent_builder = TopologyBuilder::new();
        let definition = independent_builder
            .add_molecule_definition(&ethane)
            .unwrap();
        independent_builder.add_instance(definition).unwrap();
        let independent = Arc::new(independent_builder.build().unwrap());
        assert_eq!(
            selected.ensure_compatible(&independent),
            Err(SelectionError::TopologyMismatch)
        );
    }

    #[test]
    fn topology_global_hierarchy_crosses_molecule_boundaries() {
        let mut macro_builder = crate::core::MoleculeEditor::new();
        let atom = macro_builder
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");
        let molecule = macro_builder.finish().unwrap();

        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(&molecule).unwrap();
        let first = builder.add_instance(definition).unwrap();
        let small = perceived_molecule("O");
        let small_definition = builder.add_molecule_definition(&small).unwrap();
        let small_instance = builder.add_instance(small_definition).unwrap();
        let second = builder.add_instance(definition).unwrap();
        let chain = builder
            .hierarchy_mut()
            .add_chain("A", Some("AUTH".into()))
            .unwrap();
        let first_residue = builder
            .hierarchy_mut()
            .add_residue(chain, "GLY", Some(1), Some("10".into()), None)
            .unwrap();
        let second_residue = builder
            .hierarchy_mut()
            .add_residue(chain, "GLY", Some(2), Some("11".into()), None)
            .unwrap();
        let first_atom = InstanceAtomId::new(first, atom);
        let second_atom = InstanceAtomId::new(second, atom);
        let first_site = builder
            .hierarchy_mut()
            .add_atom_site(first_residue, first_atom, AtomSiteMetadata::default())
            .unwrap();
        let second_site = builder
            .hierarchy_mut()
            .add_atom_site(second_residue, second_atom, AtomSiteMetadata::default())
            .unwrap();
        let topology = Arc::new(builder.build().unwrap());
        assert_eq!(topology.definition_count(), 2);
        assert_eq!(topology.hierarchy().chains().count(), 1);
        let first_molecule = topology.molecule(first).unwrap();
        assert_eq!(first_molecule.molecule(), &molecule);
        assert_eq!(
            first_molecule
                .chains()
                .map(ChainView::id)
                .collect::<Vec<_>>(),
            vec![chain]
        );
        assert!(topology
            .molecule(small_instance)
            .unwrap()
            .chains()
            .next()
            .is_none());

        assert_eq!(
            topology.chains().map(ChainView::id).collect::<Vec<_>>(),
            vec![chain]
        );
        assert_eq!(
            topology.residues().map(ResidueView::id).collect::<Vec<_>>(),
            vec![first_residue, second_residue]
        );
        assert_eq!(
            topology
                .atom_sites()
                .map(AtomSiteView::id)
                .collect::<Vec<_>>(),
            vec![first_site, second_site]
        );

        let chain_view = topology.chain(chain).unwrap();
        assert_eq!(chain_view.residues().count(), 2);

        assert_eq!(topology.atom_for_site(first_site).unwrap(), first_atom);
        assert_eq!(
            topology
                .atom_site_for_atom(first_atom)
                .unwrap()
                .unwrap()
                .id(),
            first_site
        );
        assert_eq!(
            topology.residue_for_atom(first_atom).unwrap().unwrap().id(),
            first_residue
        );
        assert_eq!(
            topology.chain_for_atom(first_atom).unwrap().unwrap().id(),
            chain
        );
        assert_eq!(
            topology.residue_for_site(first_site).unwrap().id(),
            first_residue
        );
        assert_eq!(
            topology.chain_for_residue(first_residue).unwrap().id(),
            chain
        );

        let small_atom = topology
            .atom_ids()
            .iter()
            .copied()
            .find(|id| id.molecule() == small_instance)
            .unwrap();
        assert!(topology.atom_site_for_atom(small_atom).unwrap().is_none());
        assert!(topology.residue_for_atom(small_atom).unwrap().is_none());
        assert!(topology.chain_for_atom(small_atom).unwrap().is_none());

        assert_eq!(
            AtomSelection::for_chains(&topology, [chain])
                .unwrap()
                .semantic_ids(&topology)
                .unwrap(),
            vec![first_atom, second_atom]
        );
        assert_eq!(
            AtomSelection::for_residues(&topology, [second_residue])
                .unwrap()
                .semantic_ids(&topology)
                .unwrap(),
            vec![second_atom]
        );
        assert_eq!(
            AtomSelection::for_atom_sites(&topology, [first_site])
                .unwrap()
                .semantic_ids(&topology)
                .unwrap(),
            vec![first_atom]
        );
        assert_eq!(
            AtomSelection::for_chain_label(&topology, "A")
                .unwrap()
                .semantic_ids(&topology)
                .unwrap(),
            vec![first_atom, second_atom]
        );

        let chain_selection = AtomSelection::for_chains(&topology, [chain]).unwrap();
        let subset = topology.subset(&chain_selection).unwrap();
        assert_eq!(subset.topology().instance_count(), 2);
        assert_eq!(subset.topology().chains().count(), 1);
        assert_eq!(subset.topology().residues().count(), 2);
        assert_eq!(subset.topology().atom_sites().count(), 2);
        assert_eq!(
            subset.correspondence().source_atom_indices(),
            [TopologyAtomIndex::new(0), TopologyAtomIndex::new(2)]
        );
        assert!(subset.correspondence().target_atom(small_atom).is_none());
    }

    #[test]
    fn induced_subset_splits_molecules_and_filters_hierarchy_deterministically() {
        let mut editor = crate::core::MoleculeEditor::new();
        let first = editor
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        let middle = editor
            .add_atom(Atom::new(Element::from_symbol("N").unwrap()))
            .unwrap();
        let tombstone = editor
            .add_atom(Atom::new(Element::from_symbol("H").unwrap()))
            .unwrap();
        editor.delete_atom(tombstone).unwrap();
        let last = editor
            .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
            .unwrap();
        editor.add_bond(first, middle, BondOrder::Single).unwrap();
        editor.add_bond(middle, last, BondOrder::Double).unwrap();
        let molecule = editor.finish().unwrap();
        assert_eq!(last.raw(), 3);

        let mut builder = TopologyBuilder::new();
        let instance = builder.add_molecule(&molecule).unwrap();
        let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
        for (sequence, atom) in [(1, first), (2, middle), (3, last)] {
            let residue = builder
                .hierarchy_mut()
                .add_residue(chain, "RES", Some(sequence), None, None)
                .unwrap();
            builder
                .hierarchy_mut()
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, atom),
                    AtomSiteMetadata::default(),
                )
                .unwrap();
        }
        let owner_key = PropertyKey::new("owner_note").unwrap();
        let value_key = PropertyKey::new("source_index").unwrap();
        builder
            .insert_property(owner_key.clone(), PropertyValue::Bool(true))
            .unwrap();
        builder
            .molecule_instance_properties_mut()
            .insert(value_key.clone(), PropertyColumn::Int(vec![Some(10)]))
            .unwrap();
        builder
            .atom_properties_mut()
            .insert(
                value_key.clone(),
                PropertyColumn::Int(vec![Some(1), Some(2), Some(3)]),
            )
            .unwrap();
        builder
            .bond_properties_mut()
            .insert(
                value_key.clone(),
                PropertyColumn::Int(vec![Some(4), Some(5)]),
            )
            .unwrap();
        builder
            .chain_properties_mut()
            .insert(value_key.clone(), PropertyColumn::Int(vec![Some(6)]))
            .unwrap();
        builder
            .residue_properties_mut()
            .insert(
                value_key.clone(),
                PropertyColumn::Int(vec![Some(7), Some(8), Some(9)]),
            )
            .unwrap();
        builder
            .atom_site_properties_mut()
            .insert(
                value_key.clone(),
                PropertyColumn::Int(vec![Some(10), Some(11), Some(12)]),
            )
            .unwrap();
        let source = Arc::new(builder.build().unwrap());
        let whole_selection = AtomSelection::from_atoms(
            &source,
            [first, middle, last].map(|atom| InstanceAtomId::new(instance, atom)),
        )
        .unwrap();
        let whole = source.subset(&whole_selection).unwrap();
        assert_eq!(
            whole
                .topology()
                .molecule_instance_properties()
                .get(&value_key),
            Some(&PropertyColumn::Int(vec![Some(10)]))
        );

        let partial_selection = AtomSelection::from_atoms(
            &source,
            [first, middle].map(|atom| InstanceAtomId::new(instance, atom)),
        )
        .unwrap();
        let partial = source.subset(&partial_selection).unwrap();
        assert!(!partial.topology().molecule_instance_properties().has_data());

        let selection = AtomSelection::from_atoms(
            &source,
            [
                InstanceAtomId::new(instance, first),
                InstanceAtomId::new(instance, last),
            ],
        )
        .unwrap();

        let subset = source.subset(&selection).unwrap();
        assert_eq!(subset.topology().instance_count(), 2);
        assert_eq!(subset.topology().atom_count(), 2);
        assert_eq!(subset.topology().bond_count(), 0);
        assert_eq!(subset.topology().chains().count(), 1);
        assert_eq!(subset.topology().residues().count(), 2);
        assert_eq!(subset.topology().atom_sites().count(), 2);
        assert_eq!(
            subset
                .topology()
                .residues()
                .map(|residue| residue.label_seq_id())
                .collect::<Vec<_>>(),
            [Some(1), Some(3)]
        );
        assert_eq!(
            subset.correspondence().source_atom_indices(),
            [TopologyAtomIndex::new(0), TopologyAtomIndex::new(2)]
        );
        assert!(subset
            .correspondence()
            .target_atom(InstanceAtomId::new(instance, middle))
            .is_none());
        let target_last = subset
            .correspondence()
            .target_atom(InstanceAtomId::new(instance, last))
            .expect("selected tombstone-separated atom is mapped");
        assert_eq!(subset.topology().atom_index(target_last).unwrap().raw(), 1);
        let projected = subset.topology();
        assert!(projected.properties().get(&owner_key).is_none());
        assert!(!projected.molecule_instance_properties().has_data());
        assert_eq!(
            projected.atom_properties().get(&value_key).unwrap(),
            &PropertyColumn::Int(vec![Some(1), Some(3)])
        );
        assert!(!projected.bond_properties().has_data());
        assert_eq!(
            projected.chain_properties().get(&value_key).unwrap(),
            &PropertyColumn::Int(vec![Some(6)])
        );
        assert_eq!(
            projected.residue_properties().get(&value_key).unwrap(),
            &PropertyColumn::Int(vec![Some(7), Some(9)])
        );
        assert_eq!(
            projected.atom_site_properties().get(&value_key).unwrap(),
            &PropertyColumn::Int(vec![Some(10), Some(12)])
        );
    }

    #[test]
    fn retaining_one_whole_instance_projects_exactly_one_instance_property_row() {
        let mut editor = crate::core::MoleculeEditor::new();
        editor
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        let molecule = editor.finish().unwrap();
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(&molecule).unwrap();
        builder.add_instance(definition).unwrap();
        let second = builder.add_instance(definition).unwrap();
        let key = PropertyKey::new("instance_score").unwrap();
        builder
            .molecule_instance_properties_mut()
            .insert(key.clone(), PropertyColumn::Int(vec![Some(10), Some(20)]))
            .unwrap();
        let source = Arc::new(builder.build().unwrap());

        let retained = transform::retain_instances(&source, [second]).unwrap();
        assert_eq!(retained.instance_count(), 1);
        assert_eq!(
            retained.molecule_instance_properties().get(&key),
            Some(&PropertyColumn::Int(vec![Some(20)]))
        );
    }
}
