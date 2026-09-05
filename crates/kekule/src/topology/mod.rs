//! Immutable coordinate-free molecular systems and biological hierarchy.
//!
//! A [`Topology`] turns connected molecular definitions into explicit
//! occurrences in one system. Repeated definitions can be reused—for example,
//! a solvent box can contain many instances of one water definition—while each
//! instance still receives distinct [`InstanceAtomId`] and [`InstanceBondId`]
//! identities.
//!
//! The system-wide [`Hierarchy`] organizes those atoms into chains, residues,
//! and atom sites. Hierarchy is independent of covalent connectedness: one
//! chain may span several molecule instances, and one connected molecule may
//! contribute to several chains.
//!
//! Topologies are immutable after publication. Use [`TopologyBuilder`] for
//! assembly, [`Topology::into_builder`] for append-oriented transformation, and
//! [`transform`] or selections for checked structural subsets. Explicit
//! [`Topology::perceived`] derives chemistry in a new snapshot with the same layout.

mod builder;
mod classification;
mod hierarchy;
mod perception;
mod selection;
pub mod transform;

#[cfg(test)]
mod tests;

use std::collections::BTreeMap;
use std::fmt;

use crate::core::{Atom, AtomId, Bond, BondId, Molecule};
use crate::properties::{Properties, PropertyError, PropertyKey, PropertyTable, PropertyValue};
pub use builder::{TopologyBuildError, TopologyBuilder, TopologyHierarchyError, TopologyIdKind};
pub use hierarchy::{
    AtomSite, AtomSiteId, AtomSiteMetadata, Chain, ChainId, Hierarchy, HierarchyError,
    HierarchyIdKind, Residue, ResidueId,
};
pub use perception::TopologyPerceptionError;
pub use selection::{AtomSelection, SelectionError};

fixed_u32_id!(MoleculeDefinitionId, "definition");
fixed_u32_id!(MoleculeInstanceId, "molecule");
fixed_u32_id!(TopologyAtomIndex, "atom-index");
fixed_u32_id!(TopologyBondIndex, "bond-index");

/// Broad intrinsic classification of one reusable topology molecule definition.
///
/// The class is inferred when a [`Topology`] is published, may use hierarchy
/// context, and is shared by every instance of the definition. It is distinct
/// from contextual roles such as ligand or receptor and from format-specific
/// categories such as mmCIF entity kinds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MoleculeClass {
    Protein,
    Dna,
    Rna,
    Carbohydrate,
    Water,
    Ion,
    SmallMolecule,
    Other,
}

/// Broad canonical classification of one topology-owned hierarchy residue.
///
/// Residue classes are inferred during topology publication and can be
/// explicitly overridden through [`TopologyBuilder`]. They describe component
/// identity rather than contextual roles or source-format entity semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ResidueClass {
    AminoAcid,
    DnaNucleotide,
    RnaNucleotide,
    Carbohydrate,
    Water,
    Ion,
    Other,
}

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

    /// Returns this topology-owned residue's canonical class.
    pub fn class(self) -> ResidueClass {
        self.local().class()
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
    class: MoleculeClass,
}

impl MoleculeDefinition {
    pub const fn id(&self) -> MoleculeDefinitionId {
        self.id
    }

    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }

    /// Returns the canonical class shared by every instance of this definition.
    pub const fn class(&self) -> MoleculeClass {
        self.class
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

    /// Returns the canonical class of this occurrence's reusable definition.
    pub fn class(self) -> MoleculeClass {
        self.definition().class()
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
/// `Topology` owns reusable molecule definitions, explicit instances of those
/// definitions, topology-global dense atom and bond order, one system
/// [`Hierarchy`], and topology-scoped properties. It does not own coordinates.
/// Attach positions through [`crate::structure::Model`] or another realization
/// type.
///
/// Local [`crate::core::AtomId`] and [`crate::core::BondId`] values identify
/// entities inside one definition. At system scope they are qualified with a
/// [`MoleculeInstanceId`] as [`InstanceAtomId`] and [`InstanceBondId`]. This
/// distinction matters whenever a definition is reused.
///
/// Exact shared ownership conventionally uses [`std::sync::Arc<Topology>`].
/// Some coordinate-dependent APIs require the same allocation, not merely an
/// independently constructed topology with equal contents. [`Self::same_layout`]
/// checks complete static layout equality when shared identity is not required.
/// Topology-changing operations return new published values.
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

    /// Consumes this published topology and stages a topology transformation.
    ///
    /// Existing definitions, instances, semantic identifiers, dense order,
    /// hierarchy, and static properties are retained. Appending through the
    /// returned builder assigns new identifiers after the retained identity
    /// spaces, and [`TopologyBuilder::build`] reconstructs derived lookups.
    pub fn into_builder(self) -> TopologyBuilder {
        TopologyBuilder::from_topology(self)
    }

    /// Builds a topology containing one explicit occurrence of `molecule`.
    ///
    /// The molecule is installed as its own definition and instance. No
    /// hierarchy is fabricated and no chemical perception is run.
    pub fn from_molecule(molecule: &Molecule) -> Result<Self, TopologyBuildError> {
        Self::from_molecules(std::slice::from_ref(molecule))
    }

    /// Builds a topology containing one explicit occurrence per input molecule.
    ///
    /// Input order becomes authoritative instance order. Each input is
    /// installed as a fresh definition; definition reuse and interning remain
    /// explicit [`TopologyBuilder`] policies. Empty input fails with
    /// [`TopologyBuildError::NoMoleculeInstances`]. No hierarchy is fabricated
    /// and no chemical perception is run.
    pub fn from_molecules(molecules: &[Molecule]) -> Result<Self, TopologyBuildError> {
        let mut builder = TopologyBuilder::new();
        builder.reserve_definitions(molecules.len())?;
        builder.reserve_instances(molecules.len())?;
        for molecule in molecules {
            builder.add_molecule(molecule)?;
        }
        builder.build()
    }

    /// Returns whether two topologies have the same complete static layout.
    ///
    /// Layout equality includes chemical and hierarchy content, definition and
    /// instance partitioning, semantic identifiers,
    /// authoritative dense atom and bond order, and the corresponding index
    /// maps. Whether two values share one `Arc` allocation is deliberately
    /// excluded, as are installed perception and generic properties. In
    /// particular, [`Self::perceived`] preserves layout equality without
    /// preserving the source snapshot's shared allocation identity.
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
    pub fn molecules(
        &self,
    ) -> impl ExactSizeIterator<Item = MoleculeInstanceView<'_>> + DoubleEndedIterator {
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

    /// Returns the canonical class shared by the instance's definition.
    pub fn molecule_class(
        &self,
        instance: MoleculeInstanceId,
    ) -> Result<MoleculeClass, TopologyError> {
        Ok(self.definition_for_instance(instance)?.class())
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
