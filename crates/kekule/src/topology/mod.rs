//! Coordinate-free molecular systems and biological hierarchy.
//!
//! A [`Topology`] turns connected molecular definitions
//! ([`MoleculeDefinition`]) into explicit occurrences in one system. Repeated
//! definitions can be reused—for example, a solvent box can contain many
//! instances of one water definition—while each instance still receives
//! distinct [`InstanceAtomId`] and [`InstanceBondId`] identities.
//!
//! Two kinds of number address a topology's atoms and bonds. Semantic IDs
//! ([`MoleculeDefinitionId`], [`MoleculeInstanceId`], [`InstanceAtomId`],
//! [`InstanceBondId`], [`ChainId`], [`ResidueId`], [`AtomSiteId`]) name
//! entities. Dense indices ([`TopologyAtomIndex`], [`TopologyBondIndex`]) are
//! the coordinate and property row order shared with structures. The topology
//! translates between them in constant time; they are never interchangeable.
//!
//! The system-wide hierarchy organizes those atoms into chains, residues,
//! and atom sites, read through [`ChainView`], [`ResidueView`], and
//! [`AtomSiteView`] and staged through [`HierarchyBuilder`]. Hierarchy is
//! independent of covalent connectedness: one chain may span several molecule
//! instances, and one connected molecule may contribute to several chains.
//! Each definition has one [`MoleculeClass`] and each residue one
//! [`ResidueClass`], inferred at publication or assigned explicitly
//! ([`ClassAssignment`]).
//!
//! Read access is through borrowed views: [`InstanceView`], [`AtomView`],
//! [`BondView`], [`ChainView`], [`ResidueView`], and [`AtomSiteView`].
//!
//! Topology structure is immutable after publication. Use [`TopologyBuilder`]
//! for assembly, [`Topology::into_builder`] for append-oriented transformation,
//! the [`edit`] module's [`TopologyEditor`](edit::TopologyEditor) for
//! structural edits, and [`Topology::subset`] with an [`AtomSelection`] or
//! [`BondSelection`] for checked structural subsets. Explicit
//! [`Topology::perceived`] derives chemistry in a new snapshot with the same
//! layout, and [`Topology::properties_mut`] edits annotations without changing
//! it.
//!
//! ```
//! use kekule::{smiles, topology::{MoleculeClass, Topology}};
//!
//! // Two separate components become two instances of two definitions.
//! let topology = Topology::from_molecules(smiles::to_molecules("CCO.[OH2]")?)?;
//! assert_eq!(topology.instance_count(), 2);
//! assert_eq!(topology.atom_count(), 4);
//! let water = topology.instances().nth(1).unwrap();
//! assert_eq!(water.class(), MoleculeClass::Water);
//! let atom = topology.atom_ids()[3];
//! assert_eq!(atom.instance(), water.id());
//! assert_eq!(topology.atom_index(atom).map(|index| index.index()), Some(3));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod builder;
mod classification;
mod components;
mod dense_map;
mod dense_order;
mod editor;
mod entity_views;
mod hierarchy;
mod hierarchy_index;
mod lookup;
mod perception;
mod reconstruct;
mod selection;
mod shared_rows;
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
pub(crate) use editor::{AppendMapping, AttachedInOccurrence, DraftHydrogenError};
pub use entity_views::{AtomView, BondView};
pub(crate) use hierarchy::{
    AtomSite, Chain, Hierarchy, RawAtomSite, RawChain, RawHierarchyError, RawResidue, Residue,
    ResidueLabels,
};
pub use hierarchy::{
    AtomSiteId, AtomSiteMetadata, ChainId, HierarchyBuilder, HierarchyError, HierarchyIdKind,
    ResidueId,
};
use hierarchy_index::HierarchyIndex;
pub(crate) use hierarchy_index::SiteNames;
pub use perception::TopologyPerceptionError;
pub(crate) use reconstruct::{RawDefinition, TopologyParts};
pub use selection::{AtomSelection, BondSelection, BondSelectionMode, SelectionError};
pub use transform::{TopologySubset, TopologySubsetCorrespondence, TopologySubsetError};

/// Transactional structural editing of a published topology.
///
/// A [`TopologyEditor`](edit::TopologyEditor) resolves source identifiers to
/// stable draft-only handles ([`EditAtomId`](edit::EditAtomId),
/// [`EditBondId`](edit::EditBondId), [`EditChainId`](edit::EditChainId),
/// [`EditResidueId`](edit::EditResidueId), and
/// [`EditAtomSiteId`](edit::EditAtomSiteId)), applies atom, bond, molecule,
/// hierarchy, class, and property edits, and publishes a new topology through
/// [`finish`](edit::TopologyEditor::finish), or through
/// [`finish_with_correspondence`](edit::TopologyEditor::finish_with_correspondence)
/// together with an [`EditCorrespondence`](edit::EditCorrespondence) from
/// handles to published identities. Failures are
/// [`TopologyEditError`](edit::TopologyEditError) values.
///
/// ```
/// use kekule::{molecule::BondOrder, smiles, topology::{Topology, edit::TopologyEditor}};
///
/// // Bonding two separate occurrences merges them into one molecule.
/// let topology = Topology::from_molecules(smiles::to_molecules("C.O")?)?;
/// let [carbon, oxygen] = [topology.atom_ids()[0], topology.atom_ids()[1]];
/// let mut editor = TopologyEditor::from_topology(topology);
/// let a = editor.atom_handle(carbon).expect("source carbon");
/// let b = editor.atom_handle(oxygen).expect("source oxygen");
/// let bond = editor.add_bond(a, b, BondOrder::Single)?;
/// let (edited, correspondence) = editor.finish_with_correspondence()?;
/// assert_eq!(edited.instance_count(), 1);
/// assert_eq!(edited.bond_count(), 1);
/// let (_, index) = correspondence.bond(bond).expect("published bond");
/// assert_eq!(index.index(), 0);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod edit {
    pub use super::editor::{
        DenseCorrespondence, EditAtomId, EditAtomSite, EditAtomSiteId, EditBond, EditBondId,
        EditChain, EditChainId, EditCorrespondence, EditMolecule, EditProperties,
        EditPropertiesMut, EditPropertyTable, EditPropertyTableMut, EditResidue, EditResidueId,
        TopologyEditError, TopologyEditor,
    };
}

fixed_u32_id!(
    /// Stable semantic identifier of one reusable [`MoleculeDefinition`] in a
    /// [`Topology`].
    ///
    /// Not a dense atom or bond index. Published definition IDs are contiguous
    /// from zero in definition order.
    MoleculeDefinitionId,
    "definition"
);
fixed_u32_id!(
    /// Stable semantic identifier of one molecule instance (occurrence) in a
    /// [`Topology`].
    ///
    /// Not a dense atom or bond index. Published instance IDs are contiguous
    /// from zero in authoritative instance order, and qualify local IDs as
    /// [`InstanceAtomId`] and [`InstanceBondId`].
    MoleculeInstanceId,
    "instance"
);
fixed_u32_id!(
    /// Dense position of one atom in a [`Topology`]'s authoritative atom order.
    ///
    /// The coordinate and property row index shared with structures, not a
    /// semantic ID; translate with [`Topology::atom_id`] and
    /// [`Topology::atom_index`]. Meaningful only for topologies sharing one
    /// layout ([`Topology::shares_layout`]).
    TopologyAtomIndex,
    "atom index"
);
fixed_u32_id!(
    /// Dense position of one bond in a [`Topology`]'s authoritative bond order.
    ///
    /// The property row index shared with structures, not a semantic ID;
    /// translate with [`Topology::bond_id`] and [`Topology::bond_index`].
    /// Meaningful only for topologies sharing one layout
    /// ([`Topology::shares_layout`]).
    TopologyBondIndex,
    "bond index"
);

/// Broad intrinsic classification of one reusable topology molecule definition.
///
/// The class is inferred when a [`Topology`] is published, may use hierarchy
/// context, and is shared by every instance of the definition. It is distinct
/// from contextual roles such as ligand or receptor and from format-specific
/// categories such as mmCIF entity kinds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum MoleculeClass {
    /// A polypeptide.
    Protein,
    /// A DNA polymer.
    Dna,
    /// An RNA polymer.
    Rna,
    /// A saccharide or polysaccharide.
    Carbohydrate,
    /// A water molecule.
    Water,
    /// A free ion.
    Ion,
    /// Any other recognized small molecule, such as a ligand or solvent.
    SmallMolecule,
    /// Unclassified, or classified with conflicting strong evidence.
    Other,
}

/// Broad canonical classification of one topology-owned hierarchy residue.
///
/// Residue classes are inferred during topology publication and can be
/// explicitly assigned through [`HierarchyBuilder::set_residue_class`]. They describe component
/// identity rather than contextual roles or source-format entity semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum ResidueClass {
    /// An amino-acid residue.
    AminoAcid,
    /// A deoxyribonucleotide residue.
    DnaNucleotide,
    /// A ribonucleotide residue.
    RnaNucleotide,
    /// A saccharide residue.
    Carbohydrate,
    /// A water residue.
    Water,
    /// An ion residue.
    Ion,
    /// Any other or unrecognized residue.
    Other,
}

/// A published class and whether it was assigned explicitly or inferred.
///
/// Explicit assignments keep their class through later publication whatever
/// the evidence; inferred classes are reevaluated when their evidence
/// changes. Builders and editors clear an explicit assignment
/// (`clear_molecule_class`, `clear_residue_class`) to return to inference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ClassAssignment<C> {
    /// Assigned explicitly by a builder or editor.
    Explicit(C),
    /// Inferred at publication.
    Inferred(C),
}

impl<C: Copy> ClassAssignment<C> {
    pub(crate) const fn new(class: C, explicit: bool) -> Self {
        if explicit {
            Self::Explicit(class)
        } else {
            Self::Inferred(class)
        }
    }

    /// The class, however it was assigned.
    pub const fn class(self) -> C {
        match self {
            Self::Explicit(class) | Self::Inferred(class) => class,
        }
    }

    /// Whether the class was assigned explicitly rather than inferred.
    pub const fn is_explicit(self) -> bool {
        matches!(self, Self::Explicit(_))
    }
}

/// The local atom of one explicit molecule instance.
///
/// A stable semantic system atom ID: a definition-local [`AtomId`] qualified by
/// the [`MoleculeInstanceId`] containing it. Not a dense index; see
/// [`Topology::atom_index`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InstanceAtomId {
    molecule: MoleculeInstanceId,
    atom: AtomId,
}

impl InstanceAtomId {
    /// Qualifies local `atom` with instance `molecule`. Nothing is validated.
    pub const fn new(molecule: MoleculeInstanceId, atom: AtomId) -> Self {
        Self { molecule, atom }
    }

    /// The molecule occurrence containing the atom.
    pub const fn instance(self) -> MoleculeInstanceId {
        self.molecule
    }

    /// The definition-local atom ID.
    pub const fn atom(self) -> AtomId {
        self.atom
    }
    /// The stable encoding of this ID for persistence: the instance's raw
    /// value in the high 32 bits and the local atom's raw value in the low 32
    /// bits. [`Self::from_bits`] restores it, and the encoding does not change
    /// between Kekule versions.
    ///
    /// ```
    /// use kekule::molecule::AtomId;
    /// use kekule::topology::{InstanceAtomId, MoleculeInstanceId};
    ///
    /// let id = InstanceAtomId::new(MoleculeInstanceId::new(2), AtomId::new(7));
    /// assert_eq!(id.to_bits(), (2 << 32) | 7);
    /// assert_eq!(InstanceAtomId::from_bits(id.to_bits()), id);
    /// ```
    pub const fn to_bits(self) -> u64 {
        ((self.molecule.raw() as u64) << 32) | self.atom.raw() as u64
    }

    /// Restores an ID from [`Self::to_bits`]. Nothing is validated.
    pub const fn from_bits(bits: u64) -> Self {
        Self::new(
            MoleculeInstanceId::new((bits >> 32) as u32),
            AtomId::new(bits as u32),
        )
    }
}

impl fmt::Display for InstanceAtomId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}:{}", self.molecule, self.atom)
    }
}

/// The local bond of one explicit molecule instance.
///
/// A stable semantic system bond ID: a definition-local [`BondId`] qualified by
/// the [`MoleculeInstanceId`] containing it. Not a dense index; see
/// [`Topology::bond_index`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InstanceBondId {
    molecule: MoleculeInstanceId,
    bond: BondId,
}

impl InstanceBondId {
    /// Qualifies local `bond` with instance `molecule`. Nothing is validated.
    pub const fn new(molecule: MoleculeInstanceId, bond: BondId) -> Self {
        Self { molecule, bond }
    }

    /// The molecule occurrence containing the bond.
    pub const fn instance(self) -> MoleculeInstanceId {
        self.molecule
    }

    /// The definition-local bond ID.
    pub const fn bond(self) -> BondId {
        self.bond
    }
    /// The stable encoding of this ID for persistence: the instance's raw
    /// value in the high 32 bits and the local bond's raw value in the low 32
    /// bits. [`Self::from_bits`] restores it, and the encoding does not change
    /// between Kekule versions.
    ///
    /// ```
    /// use kekule::molecule::BondId;
    /// use kekule::topology::{InstanceBondId, MoleculeInstanceId};
    ///
    /// let id = InstanceBondId::new(MoleculeInstanceId::new(2), BondId::new(7));
    /// assert_eq!(id.to_bits(), (2 << 32) | 7);
    /// assert_eq!(InstanceBondId::from_bits(id.to_bits()), id);
    /// ```
    pub const fn to_bits(self) -> u64 {
        ((self.molecule.raw() as u64) << 32) | self.bond.raw() as u64
    }

    /// Restores an ID from [`Self::to_bits`]. Nothing is validated.
    pub const fn from_bits(bits: u64) -> Self {
        Self::new(
            MoleculeInstanceId::new((bits >> 32) as u32),
            BondId::new(bits as u32),
        )
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

    /// This chain's identifier.
    pub const fn id(self) -> ChainId {
        self.id
    }

    /// The label chain identifier (mmCIF `label_asym_id`).
    pub fn label_id(self) -> &'a str {
        self.local().label_id()
    }

    /// The author chain identifier, or `None` when the source gave none.
    pub fn author_id(self) -> Option<&'a str> {
        self.local().author_id()
    }

    /// This chain's residues in chain order.
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

    /// This residue's identifier.
    pub const fn id(self) -> ResidueId {
        self.id
    }

    /// The chain containing this residue.
    pub fn chain(self) -> ChainView<'a> {
        ChainView::new(self.topology, self.local().chain())
    }

    /// The residue (component) name, such as `ALA` or `HOH`.
    pub fn name(self) -> &'a str {
        self.local().name()
    }

    /// Returns this topology-owned residue's canonical class.
    pub fn class(self) -> ResidueClass {
        self.local().class()
    }

    /// This residue's class and whether it was assigned explicitly.
    pub fn class_assignment(self) -> ClassAssignment<ResidueClass> {
        let residue = self.local();
        ClassAssignment::new(residue.class(), residue.explicit_class().is_some())
    }

    /// The label component identifier, or `None` when absent.
    pub fn label_comp_id(self) -> Option<&'a str> {
        self.local().label_comp_id()
    }

    /// The author component identifier, or `None` when absent.
    pub fn author_comp_id(self) -> Option<&'a str> {
        self.local().author_comp_id()
    }

    /// The label sequence number, or `None` when absent.
    pub fn label_seq_id(self) -> Option<i32> {
        self.local().label_seq_id()
    }

    /// The author sequence identifier as written by the source, or `None`.
    pub fn author_seq_id(self) -> Option<&'a str> {
        self.local().author_seq_id()
    }

    /// The insertion code, or `None` when the residue has none.
    pub fn insertion_code(self) -> Option<&'a str> {
        self.local().insertion_code()
    }

    /// This residue's atom sites in residue order.
    pub fn atom_sites(self) -> impl ExactSizeIterator<Item = AtomSiteView<'a>> + 'a {
        let topology = self.topology;
        self.local()
            .atom_sites()
            .iter()
            .map(move |&site| AtomSiteView::new(topology, site))
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

    pub(crate) fn local(self) -> Residue<'a> {
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

    /// This atom site's identifier.
    pub const fn id(self) -> AtomSiteId {
        self.id
    }

    /// The dense index of this site's atom. Constant time, and cheaper than
    /// [`Self::atom`] because it reads no atom.
    pub fn atom_index(self) -> TopologyAtomIndex {
        self.topology.hierarchy_index.site_atom(self.id)
    }

    /// The topology atom this site describes.
    pub fn atom(self) -> AtomView<'a> {
        AtomView::new(self.topology, self.atom_index())
    }

    /// The residue containing this site.
    pub fn residue(self) -> ResidueView<'a> {
        ResidueView::new(self.topology, self.local().residue())
    }

    /// Source atom-site labels.
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

    pub(crate) fn local(self) -> AtomSite<'a> {
        self.topology
            .hierarchy
            .atom_site(self.id)
            .expect("atom-site view references a validated topology hierarchy")
    }
}

/// One reusable coordinate-free molecule definition.
///
/// Snapshots share an unchanged definition's molecule instead of copying it:
/// publishing an edit copies only the definitions it changes, and perception
/// shares the definitions it does not perceive again.
#[derive(Debug, Clone)]
pub struct MoleculeDefinition {
    id: MoleculeDefinitionId,
    molecule: Arc<Molecule>,
    class: MoleculeClass,
    // Whether the installed perception is exactly what `Topology::perceived`
    // installs on this graph, so perceiving again would change nothing.
    // Published molecules never change, so the flag cannot go stale.
    default_perception: bool,
}

impl PartialEq for MoleculeDefinition {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id && self.class == other.class && self.molecule == other.molecule
    }
}

impl MoleculeDefinition {
    /// This definition's identifier.
    pub const fn id(&self) -> MoleculeDefinitionId {
        self.id
    }

    /// The definition's molecule, with definition-local atom and bond IDs.
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

    /// The class of this occurrence's definition and whether it was assigned
    /// explicitly.
    pub fn class_assignment(self) -> ClassAssignment<MoleculeClass> {
        let explicit = self
            .topology
            .molecule_class_overrides
            .contains_key(&self.definition_id());
        ClassAssignment::new(self.class(), explicit)
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
    ///
    /// Costs what it returns: the chains of every occurrence are indexed on
    /// first use, once per layout.
    pub fn chains(self) -> impl ExactSizeIterator<Item = ChainView<'a>> + 'a {
        let topology = self.topology;
        topology
            .hierarchy_index
            .instance_chains(&topology.hierarchy, self.id.index())
            .iter()
            .map(move |&chain| ChainView::new(topology, chain))
    }

    /// Complete residues that contain at least one atom of this occurrence, in
    /// hierarchy order. A residue is not clipped and may contain other atoms.
    ///
    /// Costs what it returns: the residues of every occurrence are indexed on
    /// first use, once per layout.
    pub fn residues(self) -> impl ExactSizeIterator<Item = ResidueView<'a>> + 'a {
        let topology = self.topology;
        topology
            .hierarchy_index
            .instance_residues(&topology.hierarchy, self.id.index())
            .iter()
            .map(move |&residue| ResidueView::new(topology, residue))
    }

    /// Atom sites of this occurrence's atoms, in hierarchy order.
    ///
    /// Costs what it returns: the sites of every occurrence are indexed on
    /// first use, once per layout.
    pub fn atom_sites(self) -> impl ExactSizeIterator<Item = AtomSiteView<'a>> + 'a {
        let topology = self.topology;
        topology
            .hierarchy_index
            .instance_sites(&topology.hierarchy, self.id.index())
            .iter()
            .map(move |&site| AtomSiteView::new(topology, site))
    }

    /// Qualifies a local atom ID with this occurrence. Nothing is validated.
    pub const fn qualify_atom(self, atom: AtomId) -> InstanceAtomId {
        InstanceAtomId::new(self.id, atom)
    }

    /// Qualifies a local bond ID with this occurrence. Nothing is validated.
    pub const fn qualify_bond(self, bond: BondId) -> InstanceBondId {
        InstanceBondId::new(self.id, bond)
    }
}

/// A coordinate-free molecular system with immutable structure.
///
/// `Topology` owns reusable molecule definitions, explicit instances of those
/// definitions, topology-global dense atom and bond order, one system
/// hierarchy of chains, residues, and atom sites, and topology-scoped
/// properties. It does not own coordinates.
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
    // Residue nodes record their own explicit classes.
    hierarchy: Arc<Hierarchy>,
    // Derived from `hierarchy` and `atom_order` at publication.
    hierarchy_index: Arc<HierarchyIndex>,
    // Retain explicit assignment intent separately from inferred class values.
    // It affects future editing policy, not current topology layout equality.
    molecule_class_overrides: BTreeMap<MoleculeDefinitionId, MoleculeClass>,
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
    /// An empty [`TopologyBuilder`].
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
    ///
    /// # Errors
    ///
    /// A [`TopologyBuildError`] when publication fails, as for
    /// [`Self::from_molecules`].
    pub fn from_molecule(molecule: Molecule) -> Result<Self, TopologyBuildError> {
        Self::from_molecules([molecule])
    }

    /// Builds a topology containing one explicit occurrence per input molecule.
    ///
    /// Input order becomes authoritative instance order. Each input is
    /// installed as a fresh definition; definition reuse and interning remain
    /// explicit [`TopologyBuilder`] policies. No hierarchy is fabricated and no
    /// chemical perception is run.
    ///
    /// # Errors
    ///
    /// [`TopologyBuildError::NoMoleculeInstances`] for empty input,
    /// [`TopologyBuildError::IdentifierCapacityExceeded`] when an ID space
    /// overflows, or another [`TopologyBuildError`] when publication fails.
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
                && self.hierarchy.same_layout(&other.hierarchy))
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
    /// let selection = AtomSelection::all(source.clone());
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

    /// One reusable definition, or `None` when `id` names no definition here.
    pub fn definition(&self, id: MoleculeDefinitionId) -> Option<&MoleculeDefinition> {
        self.definitions.get(id.index())
    }

    /// Reusable definitions in definition-ID order.
    pub fn definitions(&self) -> impl ExactSizeIterator<Item = &MoleculeDefinition> {
        self.definitions.iter()
    }

    /// The number of reusable definitions.
    pub fn definition_count(&self) -> usize {
        self.definitions.len()
    }

    /// The number of molecule instances (occurrences).
    pub fn instance_count(&self) -> usize {
        self.instances.len()
    }

    /// One instance-qualified molecule occurrence, or `None` when `id` names no
    /// instance here.
    pub fn instance(&self, id: MoleculeInstanceId) -> Option<InstanceView<'_>> {
        (id.index() < self.instances.len()).then(|| InstanceView::new(self, id))
    }

    /// Molecule occurrences in authoritative instance order.
    pub fn instances(
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
        self.instances()
            .filter(move |molecule| molecule.definition_id() == definition)
    }

    /// The one authoritative system-level hierarchy storage.
    pub(crate) fn hierarchy(&self) -> &Hierarchy {
        &self.hierarchy
    }

    /// The site placing a qualified atom, or `None` when it has none or is
    /// not in this topology. Constant time.
    pub(crate) fn atom_site_of(&self, atom: InstanceAtomId) -> Option<AtomSite<'_>> {
        let index = self.atom_index(atom)?;
        let site = self.hierarchy_index.atom_site(index)?;
        self.hierarchy.atom_site(site)
    }

    /// The dense atom of every hierarchy atom site, by site ID.
    pub(crate) fn site_atoms(&self) -> &[TopologyAtomIndex] {
        self.hierarchy_index.site_atoms()
    }

    /// The atom names of the hierarchy's sites, numbered on first use and
    /// shared by every snapshot of this layout.
    pub(crate) fn site_names(&self) -> &SiteNames {
        self.hierarchy_index.names(&self.hierarchy)
    }

    /// Whether `definition`'s class was assigned explicitly.
    pub(crate) fn definition_class_is_explicit(&self, definition: MoleculeDefinitionId) -> bool {
        self.molecule_class_overrides.contains_key(&definition)
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
    /// let selection = AtomSelection::all(topology.clone());
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
    pub fn chains(&self) -> impl ExactSizeIterator<Item = ChainView<'_>> {
        self.hierarchy
            .chains()
            .map(move |(id, _)| ChainView::new(self, id))
    }

    /// Every hierarchy residue in hierarchy order.
    pub fn residues(&self) -> impl ExactSizeIterator<Item = ResidueView<'_>> {
        (0..self.hierarchy.residues().len())
            .map(move |index| ResidueView::new(self, ResidueId::new(index as u32)))
    }

    /// Every hierarchy atom site in hierarchy order.
    pub fn atom_sites(&self) -> impl ExactSizeIterator<Item = AtomSiteView<'_>> {
        (0..self.hierarchy.atom_sites().len())
            .map(move |index| AtomSiteView::new(self, AtomSiteId::new(index as u32)))
    }

    /// One hierarchy chain, or `None` when `id` names no chain here.
    pub fn chain(&self, id: ChainId) -> Option<ChainView<'_>> {
        self.hierarchy.chain(id).map(|_| ChainView::new(self, id))
    }

    /// One hierarchy residue, or `None` when `id` names no residue here.
    pub fn residue(&self, id: ResidueId) -> Option<ResidueView<'_>> {
        self.hierarchy
            .residue(id)
            .map(|_| ResidueView::new(self, id))
    }

    /// One hierarchy atom site, or `None` when `id` names no site here.
    pub fn atom_site(&self, id: AtomSiteId) -> Option<AtomSiteView<'_>> {
        self.hierarchy
            .atom_site(id)
            .map(|_| AtomSiteView::new(self, id))
    }

    /// The atom site placing the atom at a dense index, or `None` when the
    /// atom has no site or `index` is not below [`Self::atom_count`]. Constant
    /// time; `atom_at(index)?.atom_site()` also reads the atom.
    pub fn atom_site_at(&self, index: TopologyAtomIndex) -> Option<AtomSiteView<'_>> {
        (index.index() < self.atom_count())
            .then(|| self.hierarchy_index.atom_site(index))
            .flatten()
            .map(|site| AtomSiteView::new(self, site))
    }

    /// One atom by qualified ID, or `None` when the atom is not in this topology.
    pub fn atom(&self, id: InstanceAtomId) -> Option<AtomView<'_>> {
        self.atom_index(id).map(|index| AtomView::new(self, index))
    }

    /// One atom by dense index, or `None` when `index` is not below
    /// [`Self::atom_count`].
    pub fn atom_at(&self, index: TopologyAtomIndex) -> Option<AtomView<'_>> {
        (index.index() < self.atom_count()).then(|| AtomView::new(self, index))
    }

    /// One bond by qualified ID, or `None` when the bond is not in this topology.
    pub fn bond(&self, id: InstanceBondId) -> Option<BondView<'_>> {
        self.bond_index(id).map(|index| BondView::new(self, index))
    }

    /// One bond by dense index, or `None` when `index` is not below
    /// [`Self::bond_count`].
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

    /// The number of atoms across all instances.
    pub fn atom_count(&self) -> usize {
        self.atom_order.len()
    }

    /// The number of bonds across all instances.
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

    /// The dense index of a qualified atom, or `None` when it is not in this
    /// topology. Constant time.
    pub fn atom_index(&self, atom: InstanceAtomId) -> Option<TopologyAtomIndex> {
        self.atom_order
            .index(atom)
            .map(|index| TopologyAtomIndex::new(index as u32))
    }

    /// The qualified atom at a dense index, or `None` when `index` is not below
    /// [`Self::atom_count`].
    pub fn atom_id(&self, index: TopologyAtomIndex) -> Option<InstanceAtomId> {
        self.atom_order.id(index.index())
    }

    /// The dense index of a qualified bond, or `None` when it is not in this
    /// topology. Constant time.
    pub fn bond_index(&self, bond: InstanceBondId) -> Option<TopologyBondIndex> {
        self.bond_order
            .index(bond)
            .map(|index| TopologyBondIndex::new(index as u32))
    }

    /// The qualified bond at a dense index, or `None` when `index` is not below
    /// [`Self::bond_count`].
    pub fn bond_id(&self, index: TopologyBondIndex) -> Option<InstanceBondId> {
        self.bond_order.id(index.index())
    }

    /// Stored instance records in instance order.
    pub(crate) fn instance_records(
        &self,
    ) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &MoleculeInstance)> {
        self.instances
            .iter()
            .map(|instance| (instance.id, instance))
    }

    pub(crate) fn instance_record(&self, id: MoleculeInstanceId) -> Option<&MoleculeInstance> {
        self.instances.get(id.index())
    }

    /// The definition molecule of one valid occurrence.
    pub(crate) fn definition_molecule(&self, instance: MoleculeInstanceId) -> &Molecule {
        let definition = self.instances[instance.index()].definition;
        self.definitions[definition.index()].molecule()
    }
}
