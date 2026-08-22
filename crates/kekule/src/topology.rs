//! Immutable coordinate-free molecular systems, qualified identities, dense
//! orderings, and compiled selections.

pub mod transform;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::sync::Arc;

use crate::bio::{
    Hierarchy, SmcraAtomSite, SmcraAtomSiteId, SmcraChain, SmcraChainId, SmcraResidue,
    SmcraResidueId,
};
use crate::core::{
    Atom, AtomId, Bond, BondId, Element, Molecule, MoleculeConnectivityError, PropMap,
};
use crate::substructure::QueryMatch;

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

/// One definition-local SMCRA chain qualified by its molecule instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InstanceChainId {
    molecule: MoleculeInstanceId,
    chain: SmcraChainId,
}

impl InstanceChainId {
    pub const fn new(molecule: MoleculeInstanceId, chain: SmcraChainId) -> Self {
        Self { molecule, chain }
    }

    pub const fn molecule(self) -> MoleculeInstanceId {
        self.molecule
    }

    pub const fn chain(self) -> SmcraChainId {
        self.chain
    }
}

impl fmt::Display for InstanceChainId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}:chain{}", self.molecule, self.chain.raw())
    }
}

/// One definition-local SMCRA residue qualified by its molecule instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InstanceResidueId {
    molecule: MoleculeInstanceId,
    residue: SmcraResidueId,
}

impl InstanceResidueId {
    pub const fn new(molecule: MoleculeInstanceId, residue: SmcraResidueId) -> Self {
        Self { molecule, residue }
    }

    pub const fn molecule(self) -> MoleculeInstanceId {
        self.molecule
    }

    pub const fn residue(self) -> SmcraResidueId {
        self.residue
    }
}

impl fmt::Display for InstanceResidueId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}:residue{}", self.molecule, self.residue.raw())
    }
}

/// One definition-local SMCRA atom site qualified by its molecule instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InstanceAtomSiteId {
    molecule: MoleculeInstanceId,
    atom_site: SmcraAtomSiteId,
}

impl InstanceAtomSiteId {
    pub const fn new(molecule: MoleculeInstanceId, atom_site: SmcraAtomSiteId) -> Self {
        Self {
            molecule,
            atom_site,
        }
    }

    pub const fn molecule(self) -> MoleculeInstanceId {
        self.molecule
    }

    pub const fn atom_site(self) -> SmcraAtomSiteId {
        self.atom_site
    }
}

impl fmt::Display for InstanceAtomSiteId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}:atom-site{}",
            self.molecule,
            self.atom_site.raw()
        )
    }
}

/// A definition-owned SMCRA chain borrowed through one molecule instance.
///
/// Hierarchical navigation remains instance-qualified. Use [`Self::local`]
/// only when definition-local identity is explicitly required.
#[derive(Clone, Copy)]
pub struct InstanceChain<'a> {
    topology: &'a Topology,
    id: InstanceChainId,
}

impl fmt::Debug for InstanceChain<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("InstanceChain")
            .field("id", &self.id)
            .finish()
    }
}

impl<'a> InstanceChain<'a> {
    const fn new(topology: &'a Topology, id: InstanceChainId) -> Self {
        Self { topology, id }
    }

    pub const fn id(self) -> InstanceChainId {
        self.id
    }

    pub fn label_id(self) -> &'a str {
        self.local().label_id()
    }

    pub fn author_id(self) -> Option<&'a str> {
        self.local().author_id()
    }

    pub fn residues(self) -> impl ExactSizeIterator<Item = InstanceResidue<'a>> + 'a {
        let topology = self.topology;
        let molecule = self.id.molecule();
        self.local().residues().iter().copied().map(move |residue| {
            InstanceResidue::new(topology, InstanceResidueId::new(molecule, residue))
        })
    }

    pub fn props(self) -> &'a PropMap {
        self.local().props()
    }

    /// Returns the underlying definition-local node.
    pub fn local(self) -> &'a SmcraChain {
        self.topology
            .local_chain(self.id)
            .expect("instance chain view references a validated local chain")
    }
}

/// A definition-owned SMCRA residue borrowed through one molecule instance.
///
/// Hierarchical navigation remains instance-qualified. Use [`Self::local`]
/// only when definition-local identity is explicitly required.
#[derive(Clone, Copy)]
pub struct InstanceResidue<'a> {
    topology: &'a Topology,
    id: InstanceResidueId,
}

impl fmt::Debug for InstanceResidue<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("InstanceResidue")
            .field("id", &self.id)
            .finish()
    }
}

impl<'a> InstanceResidue<'a> {
    const fn new(topology: &'a Topology, id: InstanceResidueId) -> Self {
        Self { topology, id }
    }

    pub const fn id(self) -> InstanceResidueId {
        self.id
    }

    pub fn chain(self) -> InstanceChain<'a> {
        InstanceChain::new(
            self.topology,
            InstanceChainId::new(self.id.molecule(), self.local().chain()),
        )
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

    pub fn atom_sites(self) -> impl ExactSizeIterator<Item = InstanceAtomSite<'a>> + 'a {
        let topology = self.topology;
        let molecule = self.id.molecule();
        self.local().atom_sites().iter().copied().map(move |site| {
            InstanceAtomSite::new(topology, InstanceAtomSiteId::new(molecule, site))
        })
    }

    pub fn props(self) -> &'a PropMap {
        self.local().props()
    }

    /// Returns the underlying definition-local node.
    pub fn local(self) -> &'a SmcraResidue {
        self.topology
            .local_residue(self.id)
            .expect("instance residue view references a validated local residue")
    }
}

/// A definition-owned SMCRA atom site borrowed through one molecule instance.
///
/// Hierarchical navigation remains instance-qualified. Use [`Self::local`]
/// only when definition-local identity is explicitly required.
#[derive(Clone, Copy)]
pub struct InstanceAtomSite<'a> {
    topology: &'a Topology,
    id: InstanceAtomSiteId,
}

impl fmt::Debug for InstanceAtomSite<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("InstanceAtomSite")
            .field("id", &self.id)
            .finish()
    }
}

impl<'a> InstanceAtomSite<'a> {
    const fn new(topology: &'a Topology, id: InstanceAtomSiteId) -> Self {
        Self { topology, id }
    }

    pub const fn id(self) -> InstanceAtomSiteId {
        self.id
    }

    pub fn atom(self) -> InstanceAtomId {
        InstanceAtomId::new(self.id.molecule(), self.local().atom())
    }

    pub fn residue(self) -> InstanceResidue<'a> {
        InstanceResidue::new(
            self.topology,
            InstanceResidueId::new(self.id.molecule(), self.local().residue()),
        )
    }

    pub fn metadata(self) -> &'a crate::bio::SmcraAtomSiteMetadata {
        self.local().metadata()
    }

    pub fn props(self) -> &'a PropMap {
        self.local().props()
    }

    /// Returns the underlying definition-local node.
    pub fn local(self) -> &'a SmcraAtomSite {
        self.topology
            .local_atom_site(self.id)
            .expect("instance atom-site view references a validated local atom site")
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

    pub fn hierarchy(&self) -> Option<&Hierarchy> {
        (!self.molecule.hierarchy().is_empty()).then(|| self.molecule.hierarchy())
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

    pub const fn qualify_chain(&self, chain: SmcraChainId) -> InstanceChainId {
        InstanceChainId::new(self.id, chain)
    }

    pub const fn qualify_residue(&self, residue: SmcraResidueId) -> InstanceResidueId {
        InstanceResidueId::new(self.id, residue)
    }

    pub const fn qualify_atom_site(&self, atom_site: SmcraAtomSiteId) -> InstanceAtomSiteId {
        InstanceAtomSiteId::new(self.id, atom_site)
    }
}

/// One explicit molecule occurrence borrowed from a [`Topology`].
///
/// This is the instance-first system view. The underlying [`Molecule`] retains
/// definition-local identities, while atoms, bonds, and hierarchy reached
/// through this view are qualified by this occurrence's
/// [`MoleculeInstanceId`].
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

    /// Returns this occurrence's qualified hierarchy when one is present.
    pub fn hierarchy(self) -> Option<InstanceHierarchy<'a>> {
        self.definition().hierarchy().map(|_| InstanceHierarchy {
            molecule: self.id,
            topology: self.topology,
        })
    }

    pub const fn qualify_atom(self, atom: AtomId) -> InstanceAtomId {
        InstanceAtomId::new(self.id, atom)
    }

    pub const fn qualify_bond(self, bond: BondId) -> InstanceBondId {
        InstanceBondId::new(self.id, bond)
    }

    pub const fn qualify_chain(self, chain: SmcraChainId) -> InstanceChainId {
        InstanceChainId::new(self.id, chain)
    }

    pub const fn qualify_residue(self, residue: SmcraResidueId) -> InstanceResidueId {
        InstanceResidueId::new(self.id, residue)
    }

    pub const fn qualify_atom_site(self, atom_site: SmcraAtomSiteId) -> InstanceAtomSiteId {
        InstanceAtomSiteId::new(self.id, atom_site)
    }
}

/// A hierarchy borrowed through one qualified molecule instance.
#[derive(Debug, Clone, Copy)]
pub struct InstanceHierarchy<'a> {
    molecule: MoleculeInstanceId,
    topology: &'a Topology,
}

impl<'a> InstanceHierarchy<'a> {
    pub const fn molecule(&self) -> MoleculeInstanceId {
        self.molecule
    }

    pub fn chains(self) -> impl Iterator<Item = InstanceChain<'a>> + 'a {
        let molecule = self.molecule;
        let topology = self.topology;
        self.local_hierarchy().chains().map(move |(chain, _)| {
            InstanceChain::new(topology, InstanceChainId::new(molecule, chain))
        })
    }

    pub fn residues(self) -> impl Iterator<Item = InstanceResidue<'a>> + 'a {
        let molecule = self.molecule;
        let topology = self.topology;
        self.local_hierarchy().residues().map(move |(residue, _)| {
            InstanceResidue::new(topology, InstanceResidueId::new(molecule, residue))
        })
    }

    pub fn atom_sites(self) -> impl Iterator<Item = InstanceAtomSite<'a>> + 'a {
        let molecule = self.molecule;
        let topology = self.topology;
        self.local_hierarchy().atom_sites().map(move |(site, _)| {
            InstanceAtomSite::new(topology, InstanceAtomSiteId::new(molecule, site))
        })
    }

    pub fn chain(self, id: InstanceChainId) -> Result<InstanceChain<'a>, TopologyError> {
        self.ensure_molecule(id.molecule())?;
        self.topology.chain(id)
    }

    pub fn residue(self, id: InstanceResidueId) -> Result<InstanceResidue<'a>, TopologyError> {
        self.ensure_molecule(id.molecule())?;
        self.topology.residue(id)
    }

    pub fn atom_site(self, id: InstanceAtomSiteId) -> Result<InstanceAtomSite<'a>, TopologyError> {
        self.ensure_molecule(id.molecule())?;
        self.topology.atom_site(id)
    }

    pub fn atom_for_site(self, site: InstanceAtomSiteId) -> Result<InstanceAtomId, TopologyError> {
        self.ensure_molecule(site.molecule())?;
        self.topology.atom_for_site(site)
    }

    pub fn atom_site_for_atom(
        self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceAtomSite<'a>>, TopologyError> {
        self.ensure_molecule(atom.molecule())?;
        self.topology.atom_site_for_atom(atom)
    }

    pub fn residue_for_atom(
        self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceResidue<'a>>, TopologyError> {
        self.ensure_molecule(atom.molecule())?;
        self.topology.residue_for_atom(atom)
    }

    pub fn chain_for_atom(
        self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceChain<'a>>, TopologyError> {
        self.ensure_molecule(atom.molecule())?;
        self.topology.chain_for_atom(atom)
    }

    pub fn residue_for_site(
        self,
        site: InstanceAtomSiteId,
    ) -> Result<InstanceResidue<'a>, TopologyError> {
        self.ensure_molecule(site.molecule())?;
        self.topology.residue_for_site(site)
    }

    pub fn chain_for_residue(
        self,
        residue: InstanceResidueId,
    ) -> Result<InstanceChain<'a>, TopologyError> {
        self.ensure_molecule(residue.molecule())?;
        self.topology.chain_for_residue(residue)
    }

    fn local_hierarchy(self) -> &'a Hierarchy {
        self.topology
            .definition_for_instance(self.molecule)
            .expect("qualified hierarchy has a live molecule instance")
            .hierarchy()
            .expect("qualified hierarchy belongs to a macromolecule definition")
    }

    fn ensure_molecule(self, actual: MoleculeInstanceId) -> Result<(), TopologyError> {
        if actual != self.molecule {
            return Err(TopologyError::HierarchyInstanceMismatch {
                expected: self.molecule,
                actual,
            });
        }
        Ok(())
    }
}

/// An immutable, coordinate-free molecular system.
///
/// `Topology` directly owns its definitions, instances, and authoritative
/// dense layouts. It intentionally does not implement [`Clone`]; shared exact
/// topology ownership uses [`Arc<Topology>`]. [`Topology::same_layout`]
/// compares complete static layout without changing allocation identity.
///
/// Generic topology mapping is deliberately not part of the public model:
///
/// ```compile_fail
/// use kekule::topology::TopologyMapping;
/// ```
///
/// ```compile_fail
/// use kekule::topology::TopologyRemapError;
/// ```
///
/// Selection remapping tied to the generic mapping layer is also unavailable:
///
/// ```compile_fail
/// # use std::sync::Arc;
/// # use kekule::topology::{AtomSelection, Topology};
/// # let selection: AtomSelection = todo!();
/// # let topology: Arc<Topology> = todo!();
/// let _ = selection.remap_to(&topology);
/// ```
///
/// The removed identity-handle API is deliberately unavailable:
///
/// ```compile_fail
/// fn removed_identity_api(topology: &kekule::topology::Topology) {
///     use kekule::topology::TopologyIdentity;
///     let _ = topology.identity();
/// }
/// ```
///
/// Raw topologies are also deliberately not cloneable:
///
/// ```compile_fail
/// fn raw_topology_is_not_cloneable(topology: kekule::topology::Topology) {
///     let _ = topology.clone();
/// }
/// ```
///
/// Generic instance metadata and roles are deliberately not core concepts:
///
/// ```compile_fail
/// use kekule::topology::{MoleculeInstanceMetadata, MoleculeRole};
/// ```
///
/// Molecules are connected by publication, so per-instance component queries
/// are deliberately unavailable:
///
/// ```compile_fail
/// fn removed_component_api(
///     topology: &kekule::topology::Topology,
///     instance: kekule::topology::MoleculeInstanceId,
/// ) {
///     let _ = topology.connected_components(instance);
/// }
/// ```
#[derive(Debug)]
pub struct Topology {
    definitions: Vec<MoleculeDefinition>,
    instances: Vec<MoleculeInstance>,
    instance_atoms: Vec<InstanceAtomId>,
    instance_bonds: Vec<InstanceBondId>,
    atom_indices: BTreeMap<InstanceAtomId, TopologyAtomIndex>,
    bond_indices: BTreeMap<InstanceBondId, TopologyBondIndex>,
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

    pub fn hierarchy(
        &self,
        instance: MoleculeInstanceId,
    ) -> Result<Option<InstanceHierarchy<'_>>, TopologyError> {
        Ok(self
            .definition_for_instance(instance)?
            .hierarchy()
            .map(|_| InstanceHierarchy {
                molecule: instance,
                topology: self,
            }))
    }

    /// Iterates every qualified SMCRA chain in instance then hierarchy order.
    pub fn chains(&self) -> impl Iterator<Item = InstanceChain<'_>> {
        self.instances.iter().flat_map(move |instance| {
            let molecule = instance.id;
            self.definitions[instance.definition.index()]
                .hierarchy()
                .into_iter()
                .flat_map(move |hierarchy| {
                    hierarchy.chains().map(move |(chain, _)| {
                        InstanceChain::new(self, InstanceChainId::new(molecule, chain))
                    })
                })
        })
    }

    /// Iterates every qualified SMCRA residue in instance then hierarchy order.
    pub fn residues(&self) -> impl Iterator<Item = InstanceResidue<'_>> {
        self.instances.iter().flat_map(move |instance| {
            let molecule = instance.id;
            self.definitions[instance.definition.index()]
                .hierarchy()
                .into_iter()
                .flat_map(move |hierarchy| {
                    hierarchy.residues().map(move |(residue, _)| {
                        InstanceResidue::new(self, InstanceResidueId::new(molecule, residue))
                    })
                })
        })
    }

    /// Iterates every qualified SMCRA atom site in instance then hierarchy order.
    pub fn atom_sites(&self) -> impl Iterator<Item = InstanceAtomSite<'_>> {
        self.instances.iter().flat_map(move |instance| {
            let molecule = instance.id;
            self.definitions[instance.definition.index()]
                .hierarchy()
                .into_iter()
                .flat_map(move |hierarchy| {
                    hierarchy.atom_sites().map(move |(site, _)| {
                        InstanceAtomSite::new(self, InstanceAtomSiteId::new(molecule, site))
                    })
                })
        })
    }

    pub fn chain(&self, id: InstanceChainId) -> Result<InstanceChain<'_>, TopologyError> {
        self.local_chain(id)?;
        Ok(InstanceChain::new(self, id))
    }

    fn local_chain(&self, id: InstanceChainId) -> Result<&SmcraChain, TopologyError> {
        let hierarchy = self
            .definition_for_instance(id.molecule())?
            .hierarchy()
            .ok_or(TopologyError::InvalidChainId(id))?;
        hierarchy
            .chain(id.chain())
            .map_err(|_| TopologyError::InvalidChainId(id))
    }

    pub fn residue(&self, id: InstanceResidueId) -> Result<InstanceResidue<'_>, TopologyError> {
        self.local_residue(id)?;
        Ok(InstanceResidue::new(self, id))
    }

    fn local_residue(&self, id: InstanceResidueId) -> Result<&SmcraResidue, TopologyError> {
        let hierarchy = self
            .definition_for_instance(id.molecule())?
            .hierarchy()
            .ok_or(TopologyError::InvalidResidueId(id))?;
        hierarchy
            .residue(id.residue())
            .map_err(|_| TopologyError::InvalidResidueId(id))
    }

    pub fn atom_site(&self, id: InstanceAtomSiteId) -> Result<InstanceAtomSite<'_>, TopologyError> {
        self.local_atom_site(id)?;
        Ok(InstanceAtomSite::new(self, id))
    }

    fn local_atom_site(&self, id: InstanceAtomSiteId) -> Result<&SmcraAtomSite, TopologyError> {
        let hierarchy = self
            .definition_for_instance(id.molecule())?
            .hierarchy()
            .ok_or(TopologyError::InvalidAtomSiteId(id))?;
        hierarchy
            .atom_site(id.atom_site())
            .map_err(|_| TopologyError::InvalidAtomSiteId(id))
    }

    pub fn atom_for_site(&self, site: InstanceAtomSiteId) -> Result<InstanceAtomId, TopologyError> {
        let atom = self.atom_site(site)?.atom();
        self.atom(atom)?;
        Ok(atom)
    }

    pub fn atom_site_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceAtomSite<'_>>, TopologyError> {
        self.atom(atom)?;
        let Some(hierarchy) = self.definition_for_instance(atom.molecule())?.hierarchy() else {
            return Ok(None);
        };
        Ok(hierarchy.atom_site_for_atom(atom.atom()).map(|site| {
            InstanceAtomSite::new(self, InstanceAtomSiteId::new(atom.molecule(), site.id()))
        }))
    }

    pub fn residue_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceResidue<'_>>, TopologyError> {
        let Some(site) = self.atom_site_for_atom(atom)? else {
            return Ok(None);
        };
        Ok(Some(site.residue()))
    }

    pub fn chain_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceChain<'_>>, TopologyError> {
        let Some(residue) = self.residue_for_atom(atom)? else {
            return Ok(None);
        };
        Ok(Some(residue.chain()))
    }

    pub fn residue_for_site(
        &self,
        site: InstanceAtomSiteId,
    ) -> Result<InstanceResidue<'_>, TopologyError> {
        Ok(self.atom_site(site)?.residue())
    }

    pub fn chain_for_residue(
        &self,
        residue: InstanceResidueId,
    ) -> Result<InstanceChain<'_>, TopologyError> {
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
/// Macro-molecule insertion validates static graph/hierarchy consistency with
/// coordinate validation disabled. Source conformers are not scanned or
/// cloned into topology definitions.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TopologyBuilder {
    definitions: Vec<MoleculeDefinition>,
    instances: Vec<MoleculeInstance>,
}

impl TopologyBuilder {
    pub fn new() -> Self {
        Self::default()
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

    pub fn build(self) -> Result<Topology, TopologyBuildError> {
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

        Ok(Topology {
            definitions: self.definitions,
            instances: self.instances,
            instance_atoms,
            instance_bonds,
            atom_indices,
            bond_indices,
        })
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
            Self::IdentifierCapacityExceeded(kind) => {
                write!(formatter, "{kind} identifier capacity exceeded")
            }
        }
    }
}

impl std::error::Error for TopologyBuildError {}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum TopologyError {
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    InvalidAtomId(InstanceAtomId),
    InvalidBondId(InstanceBondId),
    InvalidChainId(InstanceChainId),
    InvalidResidueId(InstanceResidueId),
    InvalidAtomSiteId(InstanceAtomSiteId),
    HierarchyInstanceMismatch {
        expected: MoleculeInstanceId,
        actual: MoleculeInstanceId,
    },
    InvalidAtomIndex(TopologyAtomIndex),
    InvalidBondIndex(TopologyBondIndex),
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
            Self::HierarchyInstanceMismatch { expected, actual } => write!(
                formatter,
                "hierarchy for {expected} cannot resolve an identifier qualified by {actual}"
            ),
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid {index}"),
            Self::InvalidBondIndex(index) => write!(formatter, "invalid {index}"),
        }
    }
}

impl std::error::Error for TopologyError {}

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
        atom_sites: impl IntoIterator<Item = InstanceAtomSiteId>,
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
        residues: impl IntoIterator<Item = InstanceResidueId>,
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
                .map(InstanceAtomSite::id),
        )
    }

    pub fn for_chains(
        topology: &Arc<Topology>,
        chains: impl IntoIterator<Item = InstanceChainId>,
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
                .map(InstanceResidue::id),
        )
    }

    pub fn for_chain_label(topology: &Arc<Topology>, label: &str) -> Result<Self, SelectionError> {
        Self::for_chains(
            topology,
            topology
                .chains()
                .filter(|chain| chain.label_id() == label)
                .map(InstanceChain::id),
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
    InvalidChainId(InstanceChainId),
    InvalidResidueId(InstanceResidueId),
    InvalidAtomSiteId(InstanceAtomSiteId),
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
    use crate::bio::SmcraAtomSiteMetadata;
    use crate::core::BondOrder;
    use crate::query;
    use crate::substructure;

    fn perceived_molecule(smiles: &str) -> Molecule {
        let mut molecule = crate::tests::read_smiles(smiles).expect("SMILES should parse");
        molecule.perceive().expect("molecule should perceive");
        molecule
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
    fn reused_macro_hierarchy_navigation_and_selections_are_instance_qualified() {
        let mut macro_builder = crate::core::MoleculeEditor::new();
        let atom = macro_builder
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");
        let chain = macro_builder.hierarchy_mut().add_chain("A", None).unwrap();
        let residue = macro_builder
            .hierarchy_mut()
            .add_residue(chain, "GLY", Some(1), None, None)
            .unwrap();
        let site = macro_builder
            .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
            .unwrap();
        let molecule = macro_builder.finish().unwrap();

        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(&molecule).unwrap();
        let first = builder.add_instance(definition).unwrap();
        let small = perceived_molecule("O");
        let small_definition = builder.add_molecule_definition(&small).unwrap();
        let small_instance = builder.add_instance(small_definition).unwrap();
        let second = builder.add_instance(definition).unwrap();
        let topology = Arc::new(builder.build().unwrap());
        assert_eq!(topology.definition_count(), 2);
        assert!(topology
            .definition(definition)
            .unwrap()
            .hierarchy()
            .is_some());
        assert!(topology
            .definition(small_definition)
            .unwrap()
            .hierarchy()
            .is_none());
        let first_molecule = topology.molecule(first).unwrap();
        assert_eq!(first_molecule.molecule(), &molecule);
        assert_eq!(
            first_molecule
                .hierarchy()
                .unwrap()
                .chains()
                .map(InstanceChain::id)
                .collect::<Vec<_>>(),
            vec![InstanceChainId::new(first, chain)]
        );
        assert!(topology
            .molecule(small_instance)
            .unwrap()
            .hierarchy()
            .is_none());

        let first_chain = InstanceChainId::new(first, chain);
        let second_chain = InstanceChainId::new(second, chain);
        let first_residue = InstanceResidueId::new(first, residue);
        let second_residue = InstanceResidueId::new(second, residue);
        let first_site = InstanceAtomSiteId::new(first, site);
        let second_site = InstanceAtomSiteId::new(second, site);
        let first_atom = InstanceAtomId::new(first, atom);
        let second_atom = InstanceAtomId::new(second, atom);

        assert_ne!(first_chain, second_chain);
        assert_ne!(first_residue, second_residue);
        assert_ne!(first_site, second_site);
        assert_eq!(first_chain.chain(), second_chain.chain());
        assert_eq!(first_residue.residue(), second_residue.residue());
        assert_eq!(first_site.atom_site(), second_site.atom_site());
        assert_eq!(first_residue.to_string(), "molecule0:residue0");

        assert_eq!(
            topology.chains().map(InstanceChain::id).collect::<Vec<_>>(),
            vec![first_chain, second_chain]
        );
        assert_eq!(
            topology
                .residues()
                .map(InstanceResidue::id)
                .collect::<Vec<_>>(),
            vec![first_residue, second_residue]
        );
        assert_eq!(
            topology
                .atom_sites()
                .map(InstanceAtomSite::id)
                .collect::<Vec<_>>(),
            vec![first_site, second_site]
        );

        for (chain_id, residue_id, site_id, atom_id) in [
            (first_chain, first_residue, first_site, first_atom),
            (second_chain, second_residue, second_site, second_atom),
        ] {
            let chain_view = topology.chain(chain_id).unwrap();
            assert_eq!(chain_view.id(), chain_id);
            assert_eq!(chain_view.label_id(), "A");
            assert_eq!(chain_view.local().id(), chain);
            let residues = chain_view.residues().collect::<Vec<_>>();
            assert_eq!(residues.len(), 1);
            assert_eq!(residues[0].id(), residue_id);
            assert_eq!(residues[0].chain().id(), chain_id);
            assert_eq!(residues[0].name(), "GLY");
            let atom_sites = residues[0].atom_sites().collect::<Vec<_>>();
            assert_eq!(atom_sites.len(), 1);
            assert_eq!(atom_sites[0].id(), site_id);
            assert_eq!(atom_sites[0].residue().id(), residue_id);
            assert_eq!(atom_sites[0].atom(), atom_id);
        }
        assert!(std::ptr::eq(
            topology.chain(first_chain).unwrap().local(),
            topology.chain(second_chain).unwrap().local()
        ));

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
            first_chain
        );
        assert_eq!(
            topology.residue_for_site(first_site).unwrap().id(),
            first_residue
        );
        assert_eq!(
            topology.chain_for_residue(first_residue).unwrap().id(),
            first_chain
        );

        let scoped = topology.hierarchy(first).unwrap().unwrap();
        assert_eq!(
            scoped.chains().map(InstanceChain::id).collect::<Vec<_>>(),
            vec![first_chain]
        );
        assert_eq!(
            scoped.atom_site_for_atom(first_atom).unwrap().unwrap().id(),
            first_site
        );
        assert!(matches!(
            scoped.chain(second_chain),
            Err(TopologyError::HierarchyInstanceMismatch {
                expected,
                actual,
            }) if expected == first && actual == second
        ));

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
            AtomSelection::for_chains(&topology, [first_chain])
                .unwrap()
                .semantic_ids(&topology)
                .unwrap(),
            vec![first_atom]
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
    }
}
