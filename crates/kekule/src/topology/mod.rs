//! Coordinate-free molecular systems and biological hierarchy.
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
//! Topology structure is immutable after publication. Use [`TopologyBuilder`]
//! for assembly, [`Topology::into_builder`] for append-oriented transformation,
//! [`edit::TopologyEditor`] for structural edits, and [`Topology::subset`]
//! with selections for checked structural subsets. Explicit
//! [`Topology::perceived`] derives chemistry in a new snapshot with the same
//! layout, and [`Topology::properties_mut`] edits annotations without changing
//! it.

mod builder;
mod classification;
mod components;
mod dense_order;
mod editor;
mod entity_views;
mod hierarchy;
mod lookup;
mod perception;
mod selection;
mod transform;

#[cfg(test)]
mod tests;

use std::collections::BTreeMap;
use std::fmt;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use crate::molecule::{AtomId, BondId, Molecule};
use crate::properties::{
    PropertyKey, PropertyValue, PropertyValueRef, TopologyProperties, TopologyPropertiesMut,
};
pub use builder::{TopologyBuildError, TopologyBuilder, TopologyIdKind};
use dense_order::DenseOrder;
pub(crate) use editor::AppendMapping;
pub use entity_views::{AtomView, BondView};
pub use hierarchy::{
    AtomSite, AtomSiteId, AtomSiteMetadata, Chain, ChainId, Hierarchy, HierarchyError,
    HierarchyIdKind, Residue, ResidueId,
};
pub use perception::TopologyPerceptionError;
pub use selection::{AtomSelection, BondSelection, BondSelectionMode, SelectionError};
pub use transform::{TopologySubset, TopologySubsetCorrespondence, TopologySubsetError};

/// Transactional structural editing of a published topology.
///
/// A [`TopologyEditor`](edit::TopologyEditor) resolves source identifiers to stable draft-only
/// handles, applies atom, bond, molecule, and hierarchy edits, and publishes a
/// new topology through `finish()` together with an
/// [`EditCorrespondence`](edit::EditCorrespondence) from source to target
/// identities.
pub mod edit {
    pub use super::editor::{
        EditAtomId, EditAtomSite, EditAtomSiteId, EditBond, EditBondId, EditChain, EditChainId,
        EditCorrespondence, EditMolecule, EditResidue, EditResidueId, TopologyEditError,
        TopologyEditor,
    };
}

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

    /// One static topology annotation of this node.
    pub fn property(self, key: &PropertyKey) -> Option<PropertyValue> {
        self.property_ref(key).map(PropertyValueRef::to_value)
    }

    /// Borrows one static annotation without copying string storage.
    pub fn property_ref(self, key: &PropertyKey) -> Option<PropertyValueRef<'a>> {
        self.topology
            .properties
            .chains()
            .value_ref(key, self.id)
            .expect("dense hierarchy rows cover every node")
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

    /// One static topology annotation of this node.
    pub fn property(self, key: &PropertyKey) -> Option<PropertyValue> {
        self.property_ref(key).map(PropertyValueRef::to_value)
    }

    /// Borrows one static annotation without copying string storage.
    pub fn property_ref(self, key: &PropertyKey) -> Option<PropertyValueRef<'a>> {
        self.topology
            .properties
            .residues()
            .value_ref(key, self.id)
            .expect("dense hierarchy rows cover every node")
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

    /// The topology atom this site describes.
    pub fn atom(self) -> AtomView<'a> {
        self.topology
            .atom(self.local().atom())
            .expect("atom sites name topology atoms")
    }

    pub fn residue(self) -> ResidueView<'a> {
        ResidueView::new(self.topology, self.local().residue())
    }

    pub fn metadata(self) -> &'a AtomSiteMetadata {
        self.local().metadata()
    }

    /// One static topology annotation of this node.
    pub fn property(self, key: &PropertyKey) -> Option<PropertyValue> {
        self.property_ref(key).map(PropertyValueRef::to_value)
    }

    /// Borrows one static annotation without copying string storage.
    pub fn property_ref(self, key: &PropertyKey) -> Option<PropertyValueRef<'a>> {
        self.topology
            .properties
            .atom_sites()
            .value_ref(key, self.id)
            .expect("dense hierarchy rows cover every node")
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
pub(crate) struct MoleculeInstance {
    id: MoleculeInstanceId,
    definition: MoleculeDefinitionId,
}

impl MoleculeInstance {
    pub(crate) const fn id(&self) -> MoleculeInstanceId {
        self.id
    }

    pub(crate) const fn definition(&self) -> MoleculeDefinitionId {
        self.definition
    }

    pub(crate) const fn qualify_atom(&self, atom: AtomId) -> InstanceAtomId {
        InstanceAtomId::new(self.id, atom)
    }
}

/// One explicit molecule occurrence borrowed from a [`Topology`].
///
/// This is the instance-first system view. The underlying [`Molecule`] retains
/// definition-local identities, while atoms and bonds reached through this
/// view are qualified by this occurrence's [`MoleculeInstanceId`]. Hierarchy
/// methods are projections over the one topology-owned hierarchy.
#[derive(Clone, Copy)]
pub struct InstanceView<'a> {
    topology: &'a Topology,
    id: MoleculeInstanceId,
}

impl fmt::Debug for InstanceView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("InstanceView")
            .field("id", &self.id)
            .field("definition", &self.definition_id())
            .finish()
    }
}

impl<'a> InstanceView<'a> {
    const fn new(topology: &'a Topology, id: MoleculeInstanceId) -> Self {
        Self { topology, id }
    }

    /// Returns this occurrence's topology-wide molecule identity.
    pub const fn id(self) -> MoleculeInstanceId {
        self.id
    }

    /// Returns the reusable definition referenced by this occurrence.
    pub fn definition_id(self) -> MoleculeDefinitionId {
        self.topology.instances[self.id.index()].definition()
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

    /// One static topology annotation of this occurrence.
    pub fn property(self, key: &PropertyKey) -> Option<PropertyValue> {
        self.property_ref(key).map(PropertyValueRef::to_value)
    }

    /// Borrows one static annotation without copying string storage.
    pub fn property_ref(self, key: &PropertyKey) -> Option<PropertyValueRef<'a>> {
        self.topology
            .properties
            .molecule_instances()
            .value_ref(key, self.id)
            .expect("dense instance rows cover every occurrence")
    }

    /// This occurrence's atoms in local atom-ID order.
    pub fn atoms(self) -> impl Iterator<Item = AtomView<'a>> + 'a {
        let topology = self.topology;
        let instance = self.id;
        self.molecule().atom_ids().map(move |atom| {
            topology
                .atom(InstanceAtomId::new(instance, atom))
                .expect("instance atoms are topology atoms")
        })
    }

    /// This occurrence's bonds in local bond-ID order.
    pub fn bonds(self) -> impl Iterator<Item = BondView<'a>> + 'a {
        let topology = self.topology;
        let instance = self.id;
        self.molecule().bond_ids().map(move |bond| {
            topology
                .bond(InstanceBondId::new(instance, bond))
                .expect("instance bonds are topology bonds")
        })
    }

    /// Complete chains that contain at least one atom of this occurrence, in
    /// hierarchy order. A chain is not clipped and may contain other atoms.
    pub fn chains(self) -> impl Iterator<Item = ChainView<'a>> + 'a {
        let topology = self.topology;
        let chains = self
            .residue_ids()
            .into_iter()
            .map(|residue| {
                topology
                    .hierarchy
                    .residue(residue)
                    .expect("site residue")
                    .chain()
            })
            .collect::<std::collections::BTreeSet<_>>();
        chains
            .into_iter()
            .map(move |chain| ChainView::new(topology, chain))
    }

    /// Complete residues that contain at least one atom of this occurrence, in
    /// hierarchy order. A residue is not clipped and may contain other atoms.
    pub fn residues(self) -> impl Iterator<Item = ResidueView<'a>> + 'a {
        let topology = self.topology;
        self.residue_ids()
            .into_iter()
            .map(move |residue| ResidueView::new(topology, residue))
    }

    /// Atom sites of this occurrence's atoms, in hierarchy order.
    ///
    /// Costs one indexed lookup per atom of this occurrence, independent of
    /// the size of the rest of the hierarchy.
    pub fn atom_sites(self) -> impl Iterator<Item = AtomSiteView<'a>> + 'a {
        let topology = self.topology;
        self.site_ids()
            .into_iter()
            .map(move |site| AtomSiteView::new(topology, site))
    }

    fn site_ids(self) -> std::collections::BTreeSet<AtomSiteId> {
        let hierarchy = &self.topology.hierarchy;
        let instance = self.id;
        self.molecule()
            .atom_ids()
            .filter_map(|atom| hierarchy.atom_site_for_atom(InstanceAtomId::new(instance, atom)))
            .map(|site| site.id())
            .collect()
    }

    fn residue_ids(self) -> std::collections::BTreeSet<ResidueId> {
        let hierarchy = &self.topology.hierarchy;
        self.site_ids()
            .into_iter()
            .map(|site| hierarchy.atom_site(site).expect("indexed site").residue())
            .collect()
    }

    pub const fn qualify_atom(self, atom: AtomId) -> InstanceAtomId {
        InstanceAtomId::new(self.id, atom)
    }

    pub const fn qualify_bond(self, bond: BondId) -> InstanceBondId {
        InstanceBondId::new(self.id, bond)
    }
}

/// A coordinate-free molecular system with immutable structure.
///
/// `Topology` owns reusable molecule definitions, explicit instances of those
/// definitions, topology-global dense atom and bond order, one system
/// [`Hierarchy`], and topology-scoped properties. It does not own coordinates.
/// Attach positions through [`crate::structure::Model`] or another structure
/// type.
///
/// Local [`crate::molecule::AtomId`] and [`crate::molecule::BondId`] values identify
/// entities inside one definition. At system scope they are qualified with a
/// [`MoleculeInstanceId`] as [`InstanceAtomId`] and [`InstanceBondId`]. This
/// distinction matters whenever a definition is reused.
///
/// Dense atom order is stored explicitly and chosen at publication; it need not
/// keep an instance's atoms contiguous. Format interpretations keep source
/// atom-row order, so coordinate-file index `i` is dense atom `i`. Builders
/// default to instance order and accept [`TopologyBuilder::set_atom_order`];
/// editors and subsets keep surviving atoms in their source order. Dense bond
/// order is instance order, then local bond ID. Both orders map to and from
/// qualified IDs in constant time.
///
/// Shared ownership conventionally uses [`std::sync::Arc<Topology>`].
///
/// Every publication mints a new layout identity. [`Self::perceived`] (and
/// model, ensemble, or trajectory perception) installs new chemistry in a new
/// snapshot that keeps the same layout identity, because perception never
/// changes atoms, bonds, IDs, dense order, or hierarchy. Clones and annotation
/// edits through [`Self::properties_mut`] keep it as well. Consumers that only
/// address atoms and bonds by index, such as selections, structures, frame
/// buffers, readers, alignment, and potentials, accept any snapshot sharing
/// their layout ([`Self::shares_layout`]). Consumers that read perception, such
/// as prepared substructure targets, bind the exact snapshot. Independently
/// published topologies never share a layout, even when [`Self::same_layout`]
/// reports equal contents. Structure-changing operations return new values.
#[derive(Debug, Clone)]
pub struct Topology {
    definitions: Vec<MoleculeDefinition>,
    instances: Vec<MoleculeInstance>,
    // Static state a perception snapshot shares instead of copying.
    atom_order: Arc<DenseOrder<InstanceAtomId>>,
    bond_order: Arc<DenseOrder<InstanceBondId>>,
    hierarchy: Arc<Hierarchy>,
    // Retain explicit assignment intent separately from inferred class values.
    // It affects future editing policy, not current topology layout equality.
    molecule_class_overrides: BTreeMap<MoleculeDefinitionId, MoleculeClass>,
    residue_class_overrides: BTreeMap<ResidueId, ResidueClass>,
    // Shared with perception snapshots until edited.
    properties: Arc<TopologyProperties>,
    layout: LayoutId,
}

/// The identity of one published topology layout.
///
/// Minted once per publication and kept by every snapshot that preserves
/// atoms, bonds, IDs, dense order, and hierarchy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct LayoutId(u64);

impl LayoutId {
    fn mint() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        Self(NEXT.fetch_add(1, Ordering::Relaxed))
    }
}

impl Topology {
    pub fn builder() -> TopologyBuilder {
        TopologyBuilder::new()
    }

    /// Consumes this published topology and stages a topology transformation.
    ///
    /// Existing definitions, instances, semantic identifiers, dense order,
    /// hierarchy, and entity properties are retained. Appending clears inherited
    /// owner annotations because the system has changed. Appending through the
    /// returned builder assigns new identifiers after the retained identity
    /// spaces, and [`TopologyBuilder::build`] reconstructs derived lookups.
    /// Explicit class assignments survive rebuilding. Automatically inferred
    /// classes are reevaluated when hierarchy evidence for retained atoms or
    /// residues changes; unrelated append-only extension preserves them.
    pub fn into_builder(self) -> TopologyBuilder {
        TopologyBuilder::from_topology(self)
    }

    /// Builds a topology containing one explicit occurrence of `molecule`.
    ///
    /// The molecule is installed as its own definition and instance. No
    /// hierarchy is fabricated and no chemical perception is run.
    pub fn from_molecule(molecule: Molecule) -> Result<Self, TopologyBuildError> {
        Self::from_molecules([molecule])
    }

    /// Builds a topology containing one explicit occurrence per input molecule.
    ///
    /// Input order becomes authoritative instance order. Each input is
    /// installed as a fresh definition; definition reuse and interning remain
    /// explicit [`TopologyBuilder`] policies. Empty input fails with
    /// [`TopologyBuildError::NoMoleculeInstances`]. No hierarchy is fabricated
    /// and no chemical perception is run.
    pub fn from_molecules(
        molecules: impl IntoIterator<Item = Molecule>,
    ) -> Result<Self, TopologyBuildError> {
        let molecules = molecules.into_iter();
        let mut builder = TopologyBuilder::new();
        builder.reserve_definitions(molecules.size_hint().0)?;
        builder.reserve_instances(molecules.size_hint().0)?;
        for molecule in molecules {
            builder.add_molecule(molecule)?;
        }
        builder.build()
    }

    /// Returns whether two topologies have the same complete static layout.
    ///
    /// Layout equality includes chemical and hierarchy content, definition and
    /// instance partitioning, semantic identifiers, and authoritative dense
    /// atom and bond order. Layout identity is deliberately excluded, as are
    /// installed perception and generic properties: independently published
    /// equal topologies have the same layout without sharing it.
    ///
    /// This is stricter than order-independent structural equivalence. It does
    /// not perform graph isomorphism, reorder definitions or instances, or
    /// resolve repeated indistinguishable content.
    pub fn same_layout(&self, other: &Self) -> bool {
        self.shares_layout(other)
            || (self.definitions == other.definitions
                && self.instances == other.instances
                && self.atom_order == other.atom_order
                && self.bond_order == other.bond_order
                && self.hierarchy == other.hierarchy)
    }

    /// Returns whether both values are snapshots of one published layout.
    ///
    /// This holds for one topology, its clones, and snapshots derived from it
    /// by perception or property edits, which preserve every atom, bond, ID,
    /// dense index, and the hierarchy. Index-based consumers accept any
    /// snapshot sharing their layout. Independently published topologies never
    /// share a layout, even with equal contents; compare those with
    /// [`Self::same_layout`].
    ///
    /// ```
    /// use std::sync::Arc;
    /// use kekule::{smiles, topology::{AtomSelection, Topology}};
    ///
    /// let source = smiles::to_topology("c1ccccc1")?;
    /// let selection = AtomSelection::all(&source);
    /// let perceived = Arc::new(source.perceived()?);
    /// assert!(perceived.shares_layout(&source));
    /// // Selections stay usable with the perceived snapshot.
    /// assert!(selection.ensure_compatible(&perceived).is_ok());
    ///
    /// let rebuilt = smiles::to_topology("c1ccccc1")?;
    /// assert!(rebuilt.same_layout(&source) && !rebuilt.shares_layout(&source));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn shares_layout(&self, other: &Self) -> bool {
        self.layout == other.layout
    }

    pub fn definition(&self, id: MoleculeDefinitionId) -> Option<&MoleculeDefinition> {
        self.definitions.get(id.index())
    }

    /// Reusable definitions in definition-ID order.
    pub fn definitions(&self) -> impl ExactSizeIterator<Item = &MoleculeDefinition> {
        self.definitions.iter()
    }

    pub fn definition_count(&self) -> usize {
        self.definitions.len()
    }

    pub fn instance_count(&self) -> usize {
        self.instances.len()
    }

    /// One instance-qualified molecule occurrence.
    pub fn molecule(&self, id: MoleculeInstanceId) -> Option<InstanceView<'_>> {
        (id.index() < self.instances.len()).then(|| InstanceView::new(self, id))
    }

    /// Molecule occurrences in authoritative instance order.
    pub fn molecules(
        &self,
    ) -> impl ExactSizeIterator<Item = InstanceView<'_>> + DoubleEndedIterator {
        self.instances
            .iter()
            .map(|instance| InstanceView::new(self, instance.id))
    }

    /// Occurrences of one definition in instance order.
    pub fn instances_of(
        &self,
        definition: MoleculeDefinitionId,
    ) -> impl Iterator<Item = InstanceView<'_>> {
        self.molecules()
            .filter(move |molecule| molecule.definition_id() == definition)
    }

    /// Returns the one authoritative system-level hierarchy.
    pub fn hierarchy(&self) -> &Hierarchy {
        &self.hierarchy
    }

    /// Static annotations: owner values plus dense rows for molecule
    /// instances, atoms, bonds, chains, residues, and atom sites.
    pub fn properties(&self) -> &TopologyProperties {
        &self.properties
    }

    /// Length-preserving mutable access to static annotations.
    ///
    /// Annotations do not define atoms, bonds, IDs, dense order, or hierarchy,
    /// so editing them keeps this topology's layout identity, as perception
    /// does. Shared topologies are edited through [`std::sync::Arc::make_mut`]
    /// or the owning model or collection, which replace a shared value with an
    /// edited snapshot.
    ///
    /// ```
    /// use std::sync::Arc;
    /// use kekule::{properties::{PropertyKey, PropertyValue}, smiles, topology::AtomSelection};
    ///
    /// let mut topology = smiles::to_topology("CCO")?;
    /// let selection = AtomSelection::all(&topology);
    /// let key = PropertyKey::new("source")?;
    /// Arc::make_mut(&mut topology)
    ///     .properties_mut()
    ///     .owner_mut()
    ///     .insert(key.clone(), PropertyValue::String("example".into()))?;
    ///
    /// // The selection keeps its unedited snapshot of the same layout.
    /// assert_eq!(selection.topology().properties().owner().get(&key), None);
    /// assert_eq!(
    ///     topology.properties().owner().get(&key),
    ///     Some(&PropertyValue::String("example".into()))
    /// );
    /// assert!(selection.ensure_compatible(&topology).is_ok());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn properties_mut(&mut self) -> TopologyPropertiesMut<'_> {
        TopologyPropertiesMut::new(Arc::make_mut(&mut self.properties))
    }

    /// Every hierarchy chain in hierarchy order.
    pub fn chains(&self) -> impl Iterator<Item = ChainView<'_>> {
        self.hierarchy
            .chains()
            .map(move |(id, _)| ChainView::new(self, id))
    }

    /// Every hierarchy residue in hierarchy order.
    pub fn residues(&self) -> impl Iterator<Item = ResidueView<'_>> {
        self.hierarchy
            .residues()
            .map(move |(id, _)| ResidueView::new(self, id))
    }

    /// Every hierarchy atom site in hierarchy order.
    pub fn atom_sites(&self) -> impl Iterator<Item = AtomSiteView<'_>> {
        self.hierarchy
            .atom_sites()
            .map(move |(id, _)| AtomSiteView::new(self, id))
    }

    pub fn chain(&self, id: ChainId) -> Option<ChainView<'_>> {
        self.hierarchy.chain(id).map(|_| ChainView::new(self, id))
    }

    pub fn residue(&self, id: ResidueId) -> Option<ResidueView<'_>> {
        self.hierarchy
            .residue(id)
            .map(|_| ResidueView::new(self, id))
    }

    pub fn atom_site(&self, id: AtomSiteId) -> Option<AtomSiteView<'_>> {
        self.hierarchy
            .atom_site(id)
            .map(|_| AtomSiteView::new(self, id))
    }

    /// One atom by qualified ID.
    pub fn atom(&self, id: InstanceAtomId) -> Option<AtomView<'_>> {
        self.atom_index(id).map(|index| AtomView::new(self, index))
    }

    /// One atom by dense index.
    pub fn atom_at(&self, index: TopologyAtomIndex) -> Option<AtomView<'_>> {
        (index.index() < self.atom_count()).then(|| AtomView::new(self, index))
    }

    /// One bond by qualified ID.
    pub fn bond(&self, id: InstanceBondId) -> Option<BondView<'_>> {
        self.bond_index(id).map(|index| BondView::new(self, index))
    }

    /// One bond by dense index.
    pub fn bond_at(&self, index: TopologyBondIndex) -> Option<BondView<'_>> {
        (index.index() < self.bond_count()).then(|| BondView::new(self, index))
    }

    /// Atoms in authoritative dense order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomView<'_>> + DoubleEndedIterator {
        (0..self.atom_count())
            .map(|index| AtomView::new(self, TopologyAtomIndex::new(index as u32)))
    }

    /// Bonds in authoritative dense order.
    pub fn bonds(&self) -> impl ExactSizeIterator<Item = BondView<'_>> + DoubleEndedIterator {
        (0..self.bond_count())
            .map(|index| BondView::new(self, TopologyBondIndex::new(index as u32)))
    }

    pub fn atom_count(&self) -> usize {
        self.atom_order.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bond_order.len()
    }

    /// Qualified atom IDs in authoritative dense order.
    pub fn atom_ids(&self) -> &[InstanceAtomId] {
        self.atom_order.order()
    }

    /// Qualified bond IDs in authoritative dense order.
    pub fn bond_ids(&self) -> &[InstanceBondId] {
        self.bond_order.order()
    }

    pub fn atom_index(&self, atom: InstanceAtomId) -> Option<TopologyAtomIndex> {
        self.atom_order
            .index(atom)
            .map(|index| TopologyAtomIndex::new(index as u32))
    }

    pub fn atom_id(&self, index: TopologyAtomIndex) -> Option<InstanceAtomId> {
        self.atom_order.id(index.index())
    }

    pub fn bond_index(&self, bond: InstanceBondId) -> Option<TopologyBondIndex> {
        self.bond_order
            .index(bond)
            .map(|index| TopologyBondIndex::new(index as u32))
    }

    pub fn bond_id(&self, index: TopologyBondIndex) -> Option<InstanceBondId> {
        self.bond_order.id(index.index())
    }

    /// Stored instance records in instance order.
    pub(crate) fn instances(
        &self,
    ) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &MoleculeInstance)> {
        self.instances
            .iter()
            .map(|instance| (instance.id, instance))
    }

    pub(crate) fn instance(&self, id: MoleculeInstanceId) -> Option<&MoleculeInstance> {
        self.instances.get(id.index())
    }

    /// The definition molecule of one valid occurrence.
    pub(crate) fn definition_molecule(&self, instance: MoleculeInstanceId) -> &Molecule {
        let definition = self.instances[instance.index()].definition;
        self.definitions[definition.index()].molecule()
    }
}
