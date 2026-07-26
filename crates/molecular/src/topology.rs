//! Immutable coordinate-free molecular systems, qualified identities, dense
//! orderings, mappings, and compiled selections.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::bio::{
    MacroMolecule, MacroValidateError, MacroValidateOptions, SmcraAtomSite, SmcraAtomSiteId,
    SmcraHierarchy,
};
use crate::core::{Atom, AtomId, Bond, BondId, Element, Molecule, PropMap};
use crate::small::SmallMolecule;
use crate::substructure::QueryMatch;

macro_rules! fixed_id {
    ($name:ident, $display:literal) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
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

        impl fmt::Display for $name {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(formatter, concat!($display, "{}"), self.0)
            }
        }
    };
}

fixed_id!(MoleculeDefinitionId, "definition");
fixed_id!(MoleculeInstanceId, "molecule");
fixed_id!(TopologyAtomIndex, "atom-index");
fixed_id!(TopologyBondIndex, "bond-index");

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

/// Opaque exact identity shared by clones of one topology.
#[derive(Clone)]
pub struct TopologyIdentity(Arc<IdentityToken>);

#[derive(Debug)]
struct IdentityToken;

impl fmt::Debug for TopologyIdentity {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("TopologyIdentity(..)")
    }
}

impl PartialEq for TopologyIdentity {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl Eq for TopologyIdentity {}

impl Hash for TopologyIdentity {
    fn hash<H: Hasher>(&self, state: &mut H) {
        Arc::as_ptr(&self.0).hash(state);
    }
}

/// Coordinate-free payload stored once per topology definition.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum MoleculeDefinitionPayload {
    Small(SmallMolecule),
    Macro(MacroMolecule),
}

impl MoleculeDefinitionPayload {
    pub fn graph(&self) -> &Molecule {
        match self {
            Self::Small(molecule) => molecule.graph(),
            Self::Macro(molecule) => molecule.graph(),
        }
    }

    pub fn small_molecule(&self) -> Option<&SmallMolecule> {
        match self {
            Self::Small(molecule) => Some(molecule),
            Self::Macro(_) => None,
        }
    }

    pub fn macro_molecule(&self) -> Option<&MacroMolecule> {
        match self {
            Self::Macro(molecule) => Some(molecule),
            Self::Small(_) => None,
        }
    }

    pub fn hierarchy(&self) -> Option<&SmcraHierarchy> {
        self.macro_molecule().map(MacroMolecule::hierarchy)
    }
}

/// One reusable coordinate-free molecule definition.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeDefinition {
    id: MoleculeDefinitionId,
    payload: MoleculeDefinitionPayload,
}

impl MoleculeDefinition {
    pub const fn id(&self) -> MoleculeDefinitionId {
        self.id
    }

    pub fn payload(&self) -> &MoleculeDefinitionPayload {
        &self.payload
    }

    pub fn graph(&self) -> &Molecule {
        self.payload.graph()
    }

    pub fn small_molecule(&self) -> Option<&SmallMolecule> {
        self.payload.small_molecule()
    }

    pub fn macro_molecule(&self) -> Option<&MacroMolecule> {
        self.payload.macro_molecule()
    }

    pub fn hierarchy(&self) -> Option<&SmcraHierarchy> {
        self.payload.hierarchy()
    }
}

/// Conservative semantic roles attached to one molecule instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum MoleculeRole {
    Polymer,
    Branched,
    NonPolymer,
    Solvent,
    Ion,
    Ligand,
    Cofactor,
}

/// Static metadata unique to one occurrence of a molecule definition.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MoleculeInstanceMetadata {
    roles: BTreeSet<MoleculeRole>,
    props: PropMap,
}

impl MoleculeInstanceMetadata {
    pub fn roles(&self) -> &BTreeSet<MoleculeRole> {
        &self.roles
    }

    pub fn has_role(&self, role: MoleculeRole) -> bool {
        self.roles.contains(&role)
    }

    pub fn insert_role(&mut self, role: MoleculeRole) -> bool {
        self.roles.insert(role)
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }
}

/// One explicit occurrence of a reusable molecule definition.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeInstance {
    id: MoleculeInstanceId,
    definition: MoleculeDefinitionId,
    metadata: MoleculeInstanceMetadata,
}

impl MoleculeInstance {
    pub const fn id(&self) -> MoleculeInstanceId {
        self.id
    }

    pub const fn definition(&self) -> MoleculeDefinitionId {
        self.definition
    }

    pub fn metadata(&self) -> &MoleculeInstanceMetadata {
        &self.metadata
    }

    pub fn roles(&self) -> &BTreeSet<MoleculeRole> {
        self.metadata.roles()
    }

    pub fn has_role(&self, role: MoleculeRole) -> bool {
        self.metadata.has_role(role)
    }

    pub fn props(&self) -> &PropMap {
        self.metadata.props()
    }

    pub const fn qualify_atom(&self, atom: AtomId) -> InstanceAtomId {
        InstanceAtomId::new(self.id, atom)
    }

    pub const fn qualify_bond(&self, bond: BondId) -> InstanceBondId {
        InstanceBondId::new(self.id, bond)
    }
}

/// A hierarchy borrowed through one qualified molecule instance.
#[derive(Debug, Clone, Copy)]
pub struct InstanceSmcraHierarchy<'a> {
    molecule: MoleculeInstanceId,
    hierarchy: &'a SmcraHierarchy,
}

impl InstanceSmcraHierarchy<'_> {
    pub const fn molecule(&self) -> MoleculeInstanceId {
        self.molecule
    }

    pub fn hierarchy(&self) -> &SmcraHierarchy {
        self.hierarchy
    }

    pub fn atom_for_site(&self, site: SmcraAtomSiteId) -> Result<InstanceAtomId, TopologyError> {
        let site = self
            .hierarchy
            .atom_site(site)
            .map_err(|_| TopologyError::InvalidAtomSiteId(site))?;
        Ok(InstanceAtomId::new(self.molecule, site.atom()))
    }

    pub fn atom_site_for_atom(&self, atom: InstanceAtomId) -> Option<&SmcraAtomSite> {
        (atom.molecule == self.molecule)
            .then(|| self.hierarchy.atom_site_for_atom(atom.atom))
            .flatten()
    }
}

#[derive(Debug, PartialEq)]
struct TopologyData {
    definitions: Vec<MoleculeDefinition>,
    instances: Vec<MoleculeInstance>,
    atom_order: Vec<InstanceAtomId>,
    bond_order: Vec<InstanceBondId>,
    atom_indices: BTreeMap<InstanceAtomId, TopologyAtomIndex>,
    bond_indices: BTreeMap<InstanceBondId, TopologyBondIndex>,
}

/// An immutable, coordinate-free molecular system.
///
/// Cloning this handle is constant-time and retains exact identity.
/// [`Topology::same_layout`] compares the complete static layout, including
/// semantic identifiers and dense order, but does not make independently
/// constructed topologies compatible with topology-bound state.
#[derive(Clone)]
pub struct Topology {
    data: Arc<TopologyData>,
    identity: TopologyIdentity,
}

impl fmt::Debug for Topology {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Topology")
            .field("definitions", &self.data.definitions)
            .field("instances", &self.data.instances)
            .field("atom_order", &self.data.atom_order)
            .field("bond_order", &self.data.bond_order)
            .finish()
    }
}

impl Topology {
    pub fn builder() -> TopologyBuilder {
        TopologyBuilder::new()
    }

    pub fn identity(&self) -> TopologyIdentity {
        self.identity.clone()
    }

    pub fn same_identity(&self, other: &Self) -> bool {
        self.identity == other.identity
    }

    /// Returns whether two topologies have the same complete static layout.
    ///
    /// Layout equality includes chemical and hierarchy content, definition and
    /// instance partitioning, instance metadata, semantic identifiers,
    /// authoritative dense atom and bond order, and the corresponding index
    /// maps. Exact topology identity is deliberately excluded.
    ///
    /// This is stricter than order-independent structural equivalence. It does
    /// not perform graph isomorphism, reorder definitions or instances, or
    /// resolve repeated indistinguishable content.
    pub fn same_layout(&self, other: &Self) -> bool {
        self.data.as_ref() == other.data.as_ref()
    }

    pub fn definition(
        &self,
        id: MoleculeDefinitionId,
    ) -> Result<&MoleculeDefinition, TopologyError> {
        self.data
            .definitions
            .get(id.index())
            .ok_or(TopologyError::InvalidMoleculeDefinitionId(id))
    }

    pub fn definitions(
        &self,
    ) -> impl ExactSizeIterator<Item = (MoleculeDefinitionId, &MoleculeDefinition)> {
        self.data
            .definitions
            .iter()
            .map(|definition| (definition.id, definition))
    }

    pub fn definition_count(&self) -> usize {
        self.data.definitions.len()
    }

    pub fn instance(&self, id: MoleculeInstanceId) -> Result<&MoleculeInstance, TopologyError> {
        self.data
            .instances
            .get(id.index())
            .ok_or(TopologyError::InvalidMoleculeInstanceId(id))
    }

    pub fn instances(
        &self,
    ) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &MoleculeInstance)> {
        self.data
            .instances
            .iter()
            .map(|instance| (instance.id, instance))
    }

    pub fn instance_count(&self) -> usize {
        self.data.instances.len()
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
            .data
            .instances
            .iter()
            .filter(move |instance| instance.definition == definition))
    }

    pub fn graph_for_instance(
        &self,
        instance: MoleculeInstanceId,
    ) -> Result<&Molecule, TopologyError> {
        Ok(self.definition_for_instance(instance)?.graph())
    }

    pub fn hierarchy(
        &self,
        instance: MoleculeInstanceId,
    ) -> Result<Option<InstanceSmcraHierarchy<'_>>, TopologyError> {
        Ok(self
            .definition_for_instance(instance)?
            .hierarchy()
            .map(|hierarchy| InstanceSmcraHierarchy {
                molecule: instance,
                hierarchy,
            }))
    }

    pub fn atom(&self, id: InstanceAtomId) -> Result<&Atom, TopologyError> {
        self.graph_for_instance(id.molecule)?
            .atom(id.atom)
            .map_err(|_| TopologyError::InvalidAtomId(id))
    }

    pub fn bond(&self, id: InstanceBondId) -> Result<&Bond, TopologyError> {
        self.graph_for_instance(id.molecule)?
            .bond(id.bond)
            .map_err(|_| TopologyError::InvalidBondId(id))
    }

    pub fn atoms(&self) -> impl ExactSizeIterator<Item = (InstanceAtomId, &Atom)> {
        self.data.atom_order.iter().copied().map(|id| {
            (
                id,
                self.atom(id)
                    .expect("built topology atom order contains only live atoms"),
            )
        })
    }

    pub fn bonds(&self) -> impl ExactSizeIterator<Item = (InstanceBondId, &Bond)> {
        self.data.bond_order.iter().copied().map(|id| {
            (
                id,
                self.bond(id)
                    .expect("built topology bond order contains only live bonds"),
            )
        })
    }

    pub fn atom_count(&self) -> usize {
        self.data.atom_order.len()
    }

    pub fn bond_count(&self) -> usize {
        self.data.bond_order.len()
    }

    pub fn atom_ids(&self) -> &[InstanceAtomId] {
        &self.data.atom_order
    }

    pub fn bond_ids(&self) -> &[InstanceBondId] {
        &self.data.bond_order
    }

    pub fn atom_index(&self, atom: InstanceAtomId) -> Option<TopologyAtomIndex> {
        self.data.atom_indices.get(&atom).copied()
    }

    pub fn atom_id(&self, index: TopologyAtomIndex) -> Option<InstanceAtomId> {
        self.data.atom_order.get(index.index()).copied()
    }

    pub fn bond_index(&self, bond: InstanceBondId) -> Option<TopologyBondIndex> {
        self.data.bond_indices.get(&bond).copied()
    }

    pub fn bond_id(&self, index: TopologyBondIndex) -> Option<InstanceBondId> {
        self.data.bond_order.get(index.index()).copied()
    }

    pub fn neighbors(
        &self,
        atom: InstanceAtomId,
    ) -> Result<impl Iterator<Item = InstanceAtomId> + '_, TopologyError> {
        let graph = self.graph_for_instance(atom.molecule)?;
        graph
            .atom(atom.atom)
            .map_err(|_| TopologyError::InvalidAtomId(atom))?;
        Ok(graph
            .neighbors(atom.atom)
            .expect("validated atom has valid local adjacency")
            .map(move |neighbor| InstanceAtomId::new(atom.molecule, neighbor)))
    }

    pub fn incident_bonds(
        &self,
        atom: InstanceAtomId,
    ) -> Result<impl Iterator<Item = (InstanceBondId, &Bond)> + '_, TopologyError> {
        let graph = self.graph_for_instance(atom.molecule)?;
        graph
            .atom(atom.atom)
            .map_err(|_| TopologyError::InvalidAtomId(atom))?;
        Ok(graph
            .incident_bonds(atom.atom)
            .expect("validated atom has valid local adjacency")
            .map(move |(bond, payload)| (InstanceBondId::new(atom.molecule, bond), payload)))
    }

    pub fn connected_components(
        &self,
        instance: MoleculeInstanceId,
    ) -> Result<Vec<Vec<InstanceAtomId>>, TopologyError> {
        Ok(self
            .graph_for_instance(instance)?
            .connected_components()
            .into_iter()
            .map(|component| {
                component
                    .into_iter()
                    .map(|atom| InstanceAtomId::new(instance, atom))
                    .collect()
            })
            .collect())
    }

    pub fn implicit_hydrogens(&self, atom: InstanceAtomId) -> Result<Option<u8>, TopologyError> {
        self.atom(atom)?;
        self.graph_for_instance(atom.molecule)?
            .implicit_hydrogens(atom.atom)
            .map_err(|_| TopologyError::InvalidAtomId(atom))
    }

    pub fn atom_is_aromatic(&self, atom: InstanceAtomId) -> Result<Option<bool>, TopologyError> {
        self.atom(atom)?;
        self.graph_for_instance(atom.molecule)?
            .atom_is_aromatic(atom.atom)
            .map_err(|_| TopologyError::InvalidAtomId(atom))
    }

    pub fn bond_is_aromatic(&self, bond: InstanceBondId) -> Result<Option<bool>, TopologyError> {
        self.bond(bond)?;
        self.graph_for_instance(bond.molecule)?
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

    pub fn add_small_molecule_definition(
        &mut self,
        molecule: &SmallMolecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        validate_graph(molecule.graph())?;
        let payload = MoleculeDefinitionPayload::Small(molecule.clone_without_conformers());
        self.commit_definition(payload)
    }

    pub fn add_small_molecule_definition_owned(
        &mut self,
        molecule: SmallMolecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        validate_graph(molecule.graph())?;
        self.commit_definition(MoleculeDefinitionPayload::Small(
            molecule.without_conformers(),
        ))
    }

    pub fn add_macro_molecule_definition(
        &mut self,
        molecule: &MacroMolecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        validate_macro(molecule)?;
        let payload = MoleculeDefinitionPayload::Macro(molecule.clone_without_conformers());
        self.commit_definition(payload)
    }

    pub fn add_macro_molecule_definition_owned(
        &mut self,
        molecule: MacroMolecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        validate_macro(&molecule)?;
        self.commit_definition(MoleculeDefinitionPayload::Macro(
            molecule.without_conformers(),
        ))
    }

    pub fn add_instance(
        &mut self,
        definition: MoleculeDefinitionId,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, TopologyBuildError> {
        self.definition(definition)?;
        self.reserve_instances(1)?;
        let id = checked_id::<MoleculeInstanceId>(
            self.instances.len(),
            TopologyIdKind::MoleculeInstance,
        )?;
        self.instances.push(MoleculeInstance {
            id,
            definition,
            metadata,
        });
        Ok(id)
    }

    pub fn add_small_molecule_instance(
        &mut self,
        molecule: &SmallMolecule,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<(MoleculeDefinitionId, MoleculeInstanceId), TopologyBuildError> {
        validate_graph(molecule.graph())?;
        let payload = MoleculeDefinitionPayload::Small(molecule.clone_without_conformers());
        self.commit_definition_and_instance(payload, metadata)
    }

    pub fn add_macro_molecule_instance(
        &mut self,
        molecule: &MacroMolecule,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<(MoleculeDefinitionId, MoleculeInstanceId), TopologyBuildError> {
        validate_macro(molecule)?;
        let payload = MoleculeDefinitionPayload::Macro(molecule.clone_without_conformers());
        self.commit_definition_and_instance(payload, metadata)
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

        let atom_count = self.instances.iter().try_fold(0usize, |count, instance| {
            count
                .checked_add(
                    self.definitions[instance.definition.index()]
                        .graph()
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
                        .graph()
                        .bond_count(),
                )
                .ok_or(TopologyBuildError::IdentifierCapacityExceeded(
                    TopologyIdKind::Bond,
                ))
        })?;
        checked_future_len(0, bond_count, TopologyIdKind::Bond)?;

        let mut atom_order = Vec::new();
        let mut bond_order = Vec::new();
        let mut atom_indices = BTreeMap::new();
        let mut bond_indices = BTreeMap::new();
        atom_order
            .try_reserve_exact(atom_count)
            .map_err(|_| TopologyBuildError::IdentifierCapacityExceeded(TopologyIdKind::Atom))?;
        bond_order
            .try_reserve_exact(bond_count)
            .map_err(|_| TopologyBuildError::IdentifierCapacityExceeded(TopologyIdKind::Bond))?;

        for instance in &self.instances {
            let graph = self.definitions[instance.definition.index()].graph();
            for atom in graph.atom_ids() {
                let qualified = instance.qualify_atom(atom);
                let index =
                    checked_id::<TopologyAtomIndex>(atom_order.len(), TopologyIdKind::Atom)?;
                atom_indices.insert(qualified, index);
                atom_order.push(qualified);
            }
            for bond in graph.bond_ids() {
                let qualified = instance.qualify_bond(bond);
                let index =
                    checked_id::<TopologyBondIndex>(bond_order.len(), TopologyIdKind::Bond)?;
                bond_indices.insert(qualified, index);
                bond_order.push(qualified);
            }
        }

        Ok(Topology {
            data: Arc::new(TopologyData {
                definitions: self.definitions,
                instances: self.instances,
                atom_order,
                bond_order,
                atom_indices,
                bond_indices,
            }),
            identity: TopologyIdentity(Arc::new(IdentityToken)),
        })
    }

    fn commit_definition(
        &mut self,
        payload: MoleculeDefinitionPayload,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        self.reserve_definitions(1)?;
        let id = checked_id::<MoleculeDefinitionId>(
            self.definitions.len(),
            TopologyIdKind::MoleculeDefinition,
        )?;
        self.definitions.push(MoleculeDefinition { id, payload });
        Ok(id)
    }

    fn commit_definition_and_instance(
        &mut self,
        payload: MoleculeDefinitionPayload,
        metadata: MoleculeInstanceMetadata,
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
            payload,
        });
        self.instances.push(MoleculeInstance {
            id: instance,
            definition,
            metadata,
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
    if graph.atom_count() == 0 {
        return Err(TopologyBuildError::EmptyMoleculeDefinition);
    }
    Ok(())
}

fn validate_macro(molecule: &MacroMolecule) -> Result<(), TopologyBuildError> {
    validate_graph(molecule.graph())?;
    molecule
        .validate_with_options(MacroValidateOptions {
            validate_coordinates: false,
        })
        .map_err(TopologyBuildError::InvalidMacroMolecule)?;
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
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    InvalidMacroMolecule(MacroValidateError),
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
            Self::InvalidMoleculeDefinitionId(id) => {
                write!(formatter, "invalid molecule definition: {id}")
            }
            Self::InvalidMacroMolecule(error) => {
                write!(formatter, "invalid macromolecule definition: {error}")
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
    InvalidAtomIndex(TopologyAtomIndex),
    InvalidBondIndex(TopologyBondIndex),
    InvalidAtomSiteId(SmcraAtomSiteId),
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
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid {index}"),
            Self::InvalidBondIndex(index) => write!(formatter, "invalid {index}"),
            Self::InvalidAtomSiteId(id) => {
                write!(formatter, "invalid hierarchy atom site: {}", id.raw())
            }
        }
    }
}

impl std::error::Error for TopologyError {}

/// A topology-bound, sorted, unique dense atom selection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomSelection {
    topology: TopologyIdentity,
    indices: Vec<TopologyAtomIndex>,
}

impl AtomSelection {
    pub fn from_atoms(
        topology: &Topology,
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
            topology: topology.identity(),
            indices,
        })
    }

    pub fn from_indices(
        topology: &Topology,
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
        topology: &Topology,
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
        topology: &Topology,
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

    pub fn for_roles(
        topology: &Topology,
        roles: impl IntoIterator<Item = MoleculeRole>,
    ) -> Result<Self, SelectionError> {
        let roles = roles.into_iter().collect::<BTreeSet<_>>();
        let instances = topology
            .instances()
            .filter(|(_, instance)| instance.roles().iter().any(|role| roles.contains(role)))
            .map(|(id, _)| id);
        Self::for_instances(topology, instances)
    }

    pub fn for_elements(
        topology: &Topology,
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

    pub fn for_chain_label(topology: &Topology, label: &str) -> Result<Self, SelectionError> {
        let mut atoms = Vec::new();
        for (instance_id, _) in topology.instances() {
            let Some(qualified) = topology
                .hierarchy(instance_id)
                .map_err(|_| SelectionError::InvalidMoleculeInstanceId(instance_id))?
            else {
                continue;
            };
            let hierarchy = qualified.hierarchy();
            for (_, chain) in hierarchy
                .chains()
                .filter(|(_, chain)| chain.label_id() == label)
            {
                for residue_id in chain.residues() {
                    let residue = hierarchy
                        .residue(*residue_id)
                        .expect("validated hierarchy residue");
                    for site_id in residue.atom_sites() {
                        atoms.push(
                            qualified
                                .atom_for_site(*site_id)
                                .expect("validated hierarchy atom site"),
                        );
                    }
                }
            }
        }
        Self::from_atoms(topology, atoms)
    }

    pub fn connected_component(
        topology: &Topology,
        instance: MoleculeInstanceId,
        component: usize,
    ) -> Result<Self, SelectionError> {
        let components = topology
            .connected_components(instance)
            .map_err(|_| SelectionError::InvalidMoleculeInstanceId(instance))?;
        let atoms = components
            .get(component)
            .ok_or(SelectionError::InvalidConnectedComponent {
                instance,
                component,
            })?;
        Self::from_atoms(topology, atoms.iter().copied())
    }

    pub fn from_query_matches(
        topology: &Topology,
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

    pub fn ensure_compatible(&self, topology: &Topology) -> Result<(), SelectionError> {
        if self.topology != topology.identity {
            return Err(SelectionError::TopologyIdentityMismatch);
        }
        Ok(())
    }

    pub fn indices(&self) -> &[TopologyAtomIndex] {
        &self.indices
    }

    pub fn semantic_ids(&self, topology: &Topology) -> Result<Vec<InstanceAtomId>, SelectionError> {
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
    TopologyIdentityMismatch,
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    InvalidAtomId(InstanceAtomId),
    InvalidAtomIndex(TopologyAtomIndex),
    InvalidConnectedComponent {
        instance: MoleculeInstanceId,
        component: usize,
    },
}

impl fmt::Display for SelectionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyIdentityMismatch => {
                formatter.write_str("atom selection belongs to a different topology")
            }
            Self::InvalidMoleculeDefinitionId(id) => {
                write!(formatter, "invalid selected molecule definition: {id}")
            }
            Self::InvalidMoleculeInstanceId(id) => {
                write!(formatter, "invalid selected molecule instance: {id}")
            }
            Self::InvalidAtomId(id) => write!(formatter, "invalid selected atom: {id}"),
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid selected {index}"),
            Self::InvalidConnectedComponent {
                instance,
                component,
            } => write!(
                formatter,
                "molecule instance {instance} has no connected component {component}"
            ),
        }
    }
}

impl std::error::Error for SelectionError {}

/// Explicit old-to-new lineage for a topology-changing operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologyMapping {
    old: TopologyIdentity,
    new: TopologyIdentity,
    definitions: BTreeMap<MoleculeDefinitionId, MoleculeDefinitionId>,
    instances: BTreeMap<MoleculeInstanceId, MoleculeInstanceId>,
    atoms: BTreeMap<InstanceAtomId, InstanceAtomId>,
    bonds: BTreeMap<InstanceBondId, InstanceBondId>,
    atom_indices: BTreeMap<TopologyAtomIndex, TopologyAtomIndex>,
    bond_indices: BTreeMap<TopologyBondIndex, TopologyBondIndex>,
    removed_definitions: Vec<MoleculeDefinitionId>,
    added_definitions: Vec<MoleculeDefinitionId>,
    removed_instances: Vec<MoleculeInstanceId>,
    added_instances: Vec<MoleculeInstanceId>,
    removed_atoms: Vec<InstanceAtomId>,
    added_atoms: Vec<InstanceAtomId>,
    removed_bonds: Vec<InstanceBondId>,
    added_bonds: Vec<InstanceBondId>,
}

impl TopologyMapping {
    /// Constructs identity mappings between two topologies with identical
    /// static layouts.
    ///
    /// Independently constructed topologies may satisfy this precondition, but
    /// all semantic identifiers and authoritative dense orderings must already
    /// be identical. This method does not solve a graph isomorphism or infer a
    /// mapping between reordered layouts.
    pub fn between_identical_layouts(
        old: &Topology,
        new: &Topology,
    ) -> Result<Self, TopologyMappingError> {
        if !old.same_layout(new) {
            return Err(TopologyMappingError::NotSameLayout);
        }
        Self::from_pairs(
            old,
            new,
            old.definitions().map(|(id, _)| (id, id)),
            old.instances().map(|(id, _)| (id, id)),
            old.atom_ids().iter().copied().map(|id| (id, id)),
            old.bond_ids().iter().copied().map(|id| (id, id)),
        )
    }

    /// Constructs a checked explicit lineage mapping.
    ///
    /// Each map must be injective and reference live source and target values.
    /// Definition mappings must agree with mapped instances, instance mappings
    /// must agree with mapped atoms and bonds, and mapped bond endpoints must
    /// agree with the atom mapping. Added and removed values are derived as the
    /// unmapped complements in authoritative topology order.
    pub fn from_pairs(
        old: &Topology,
        new: &Topology,
        definitions: impl IntoIterator<Item = (MoleculeDefinitionId, MoleculeDefinitionId)>,
        instances: impl IntoIterator<Item = (MoleculeInstanceId, MoleculeInstanceId)>,
        atoms: impl IntoIterator<Item = (InstanceAtomId, InstanceAtomId)>,
        bonds: impl IntoIterator<Item = (InstanceBondId, InstanceBondId)>,
    ) -> Result<Self, TopologyMappingError> {
        let definitions = collect_mapping(
            definitions,
            |id| old.definition(id).is_ok(),
            |id| new.definition(id).is_ok(),
            MappingKind::Definition,
        )?;
        let instances = collect_mapping(
            instances,
            |id| old.instance(id).is_ok(),
            |id| new.instance(id).is_ok(),
            MappingKind::Instance,
        )?;
        let atoms = collect_mapping(
            atoms,
            |id| old.atom(id).is_ok(),
            |id| new.atom(id).is_ok(),
            MappingKind::Atom,
        )?;
        let bonds = collect_mapping(
            bonds,
            |id| old.bond(id).is_ok(),
            |id| new.bond(id).is_ok(),
            MappingKind::Bond,
        )?;

        validate_mapping_relationships(old, new, &definitions, &instances, &atoms, &bonds)?;

        let mut atom_indices = BTreeMap::new();
        for (source, target) in &atoms {
            let source_index = old
                .atom_index(*source)
                .ok_or(TopologyMappingError::InvalidSource(MappingKind::Atom))?;
            let target_index = new
                .atom_index(*target)
                .ok_or(TopologyMappingError::InvalidTarget(MappingKind::Atom))?;
            atom_indices.insert(source_index, target_index);
        }
        let mut bond_indices = BTreeMap::new();
        for (source, target) in &bonds {
            let source_index = old
                .bond_index(*source)
                .ok_or(TopologyMappingError::InvalidSource(MappingKind::Bond))?;
            let target_index = new
                .bond_index(*target)
                .ok_or(TopologyMappingError::InvalidTarget(MappingKind::Bond))?;
            bond_indices.insert(source_index, target_index);
        }
        let mapped_new_definitions = definitions.values().copied().collect::<BTreeSet<_>>();
        let mapped_new_instances = instances.values().copied().collect::<BTreeSet<_>>();
        let mapped_new_atoms = atoms.values().copied().collect::<BTreeSet<_>>();
        let mapped_new_bonds = bonds.values().copied().collect::<BTreeSet<_>>();
        let removed_definitions = old
            .definitions()
            .map(|(definition, _)| definition)
            .filter(|definition| !definitions.contains_key(definition))
            .collect();
        let added_definitions = new
            .definitions()
            .map(|(definition, _)| definition)
            .filter(|definition| !mapped_new_definitions.contains(definition))
            .collect();
        let removed_instances = old
            .instances()
            .map(|(instance, _)| instance)
            .filter(|instance| !instances.contains_key(instance))
            .collect();
        let added_instances = new
            .instances()
            .map(|(instance, _)| instance)
            .filter(|instance| !mapped_new_instances.contains(instance))
            .collect();

        Ok(Self {
            old: old.identity(),
            new: new.identity(),
            definitions,
            instances,
            atom_indices,
            bond_indices,
            removed_definitions,
            added_definitions,
            removed_instances,
            added_instances,
            removed_atoms: old
                .atom_ids()
                .iter()
                .copied()
                .filter(|atom| !atoms.contains_key(atom))
                .collect(),
            added_atoms: new
                .atom_ids()
                .iter()
                .copied()
                .filter(|atom| !mapped_new_atoms.contains(atom))
                .collect(),
            removed_bonds: old
                .bond_ids()
                .iter()
                .copied()
                .filter(|bond| !bonds.contains_key(bond))
                .collect(),
            added_bonds: new
                .bond_ids()
                .iter()
                .copied()
                .filter(|bond| !mapped_new_bonds.contains(bond))
                .collect(),
            atoms,
            bonds,
        })
    }

    pub fn old_identity(&self) -> &TopologyIdentity {
        &self.old
    }

    pub fn new_identity(&self) -> &TopologyIdentity {
        &self.new
    }

    pub fn map_definition(&self, id: MoleculeDefinitionId) -> Option<MoleculeDefinitionId> {
        self.definitions.get(&id).copied()
    }

    pub fn map_instance(&self, id: MoleculeInstanceId) -> Option<MoleculeInstanceId> {
        self.instances.get(&id).copied()
    }

    pub fn map_atom(&self, id: InstanceAtomId) -> Option<InstanceAtomId> {
        self.atoms.get(&id).copied()
    }

    pub fn map_bond(&self, id: InstanceBondId) -> Option<InstanceBondId> {
        self.bonds.get(&id).copied()
    }

    pub fn map_atom_index(&self, index: TopologyAtomIndex) -> Option<TopologyAtomIndex> {
        self.atom_indices.get(&index).copied()
    }

    pub fn map_bond_index(&self, index: TopologyBondIndex) -> Option<TopologyBondIndex> {
        self.bond_indices.get(&index).copied()
    }

    pub fn removed_definitions(&self) -> &[MoleculeDefinitionId] {
        &self.removed_definitions
    }

    pub fn added_definitions(&self) -> &[MoleculeDefinitionId] {
        &self.added_definitions
    }

    pub fn removed_instances(&self) -> &[MoleculeInstanceId] {
        &self.removed_instances
    }

    pub fn added_instances(&self) -> &[MoleculeInstanceId] {
        &self.added_instances
    }

    pub fn removed_atoms(&self) -> &[InstanceAtomId] {
        &self.removed_atoms
    }

    pub fn added_atoms(&self) -> &[InstanceAtomId] {
        &self.added_atoms
    }

    pub fn removed_bonds(&self) -> &[InstanceBondId] {
        &self.removed_bonds
    }

    pub fn added_bonds(&self) -> &[InstanceBondId] {
        &self.added_bonds
    }
}

fn validate_mapping_relationships(
    old: &Topology,
    new: &Topology,
    definitions: &BTreeMap<MoleculeDefinitionId, MoleculeDefinitionId>,
    instances: &BTreeMap<MoleculeInstanceId, MoleculeInstanceId>,
    atoms: &BTreeMap<InstanceAtomId, InstanceAtomId>,
    bonds: &BTreeMap<InstanceBondId, InstanceBondId>,
) -> Result<(), TopologyMappingError> {
    for (source, target) in instances {
        let source_definition = old
            .instance(*source)
            .map_err(|_| TopologyMappingError::InvalidSource(MappingKind::Instance))?
            .definition();
        let target_definition = new
            .instance(*target)
            .map_err(|_| TopologyMappingError::InvalidTarget(MappingKind::Instance))?
            .definition();
        let Some(expected_target_definition) = definitions.get(&source_definition).copied() else {
            return Err(TopologyMappingError::MissingDefinitionMappingForInstance {
                source: *source,
                definition: source_definition,
            });
        };
        if target_definition != expected_target_definition {
            return Err(TopologyMappingError::InconsistentInstanceDefinition {
                source: *source,
                target: *target,
                expected_target_definition,
                actual_target_definition: target_definition,
            });
        }
    }

    for (source, target) in atoms {
        let Some(expected_target_instance) = instances.get(&source.molecule()).copied() else {
            return Err(TopologyMappingError::MissingInstanceMappingForAtom { source: *source });
        };
        if target.molecule() != expected_target_instance {
            return Err(TopologyMappingError::InconsistentAtomInstance {
                source: *source,
                target: *target,
                expected_target_instance,
            });
        }
    }

    for (source, target) in bonds {
        let Some(expected_target_instance) = instances.get(&source.molecule()).copied() else {
            return Err(TopologyMappingError::MissingInstanceMappingForBond { source: *source });
        };
        if target.molecule() != expected_target_instance {
            return Err(TopologyMappingError::InconsistentBondInstance {
                source: *source,
                target: *target,
                expected_target_instance,
            });
        }

        let source_bond = old
            .bond(*source)
            .map_err(|_| TopologyMappingError::InvalidSource(MappingKind::Bond))?;
        let target_bond = new
            .bond(*target)
            .map_err(|_| TopologyMappingError::InvalidTarget(MappingKind::Bond))?;
        let (source_a, source_b) = source_bond.endpoints();
        let source_a = InstanceAtomId::new(source.molecule(), source_a);
        let source_b = InstanceAtomId::new(source.molecule(), source_b);
        let Some(mapped_a) = atoms.get(&source_a).copied() else {
            return Err(TopologyMappingError::MissingAtomMappingForBondEndpoint {
                source_bond: *source,
                source_atom: source_a,
            });
        };
        let Some(mapped_b) = atoms.get(&source_b).copied() else {
            return Err(TopologyMappingError::MissingAtomMappingForBondEndpoint {
                source_bond: *source,
                source_atom: source_b,
            });
        };
        let (target_a, target_b) = target_bond.endpoints();
        let target_a = InstanceAtomId::new(target.molecule(), target_a);
        let target_b = InstanceAtomId::new(target.molecule(), target_b);
        if !unordered_pair_eq((mapped_a, mapped_b), (target_a, target_b)) {
            return Err(TopologyMappingError::InconsistentBondEndpoints {
                source: *source,
                target: *target,
            });
        }
    }
    Ok(())
}

fn unordered_pair_eq<T: PartialEq>((left_a, left_b): (T, T), (right_a, right_b): (T, T)) -> bool {
    (left_a == right_a && left_b == right_b) || (left_a == right_b && left_b == right_a)
}

fn collect_mapping<K>(
    pairs: impl IntoIterator<Item = (K, K)>,
    source_exists: impl Fn(K) -> bool,
    target_exists: impl Fn(K) -> bool,
    kind: MappingKind,
) -> Result<BTreeMap<K, K>, TopologyMappingError>
where
    K: Copy + Ord,
{
    let mut mapping = BTreeMap::new();
    let mut targets = BTreeSet::new();
    for (source, target) in pairs {
        if !source_exists(source) {
            return Err(TopologyMappingError::InvalidSource(kind));
        }
        if !target_exists(target) {
            return Err(TopologyMappingError::InvalidTarget(kind));
        }
        if mapping.insert(source, target).is_some() {
            return Err(TopologyMappingError::DuplicateSource(kind));
        }
        if !targets.insert(target) {
            return Err(TopologyMappingError::DuplicateTarget(kind));
        }
    }
    Ok(mapping)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MappingKind {
    Definition,
    Instance,
    Atom,
    Bond,
}

impl fmt::Display for MappingKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Definition => formatter.write_str("definition"),
            Self::Instance => formatter.write_str("instance"),
            Self::Atom => formatter.write_str("atom"),
            Self::Bond => formatter.write_str("bond"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
/// Failure to construct a self-consistent explicit topology mapping or edit
/// result.
pub enum TopologyMappingError {
    /// Identity mapping was requested for layouts that differ.
    NotSameLayout,
    /// A mapped source identifier does not exist in the source topology.
    InvalidSource(MappingKind),
    /// A mapped target identifier does not exist in the target topology.
    InvalidTarget(MappingKind),
    /// A source identifier occurs more than once.
    DuplicateSource(MappingKind),
    /// More than one source identifier maps to the same target identifier.
    DuplicateTarget(MappingKind),
    /// A mapped instance's source definition has no definition mapping.
    MissingDefinitionMappingForInstance {
        source: MoleculeInstanceId,
        definition: MoleculeDefinitionId,
    },
    /// A mapped instance's target definition disagrees with the definition
    /// mapping.
    InconsistentInstanceDefinition {
        source: MoleculeInstanceId,
        target: MoleculeInstanceId,
        expected_target_definition: MoleculeDefinitionId,
        actual_target_definition: MoleculeDefinitionId,
    },
    /// A mapped atom's source molecule instance has no instance mapping.
    MissingInstanceMappingForAtom { source: InstanceAtomId },
    /// A mapped atom belongs to the wrong target molecule instance.
    InconsistentAtomInstance {
        source: InstanceAtomId,
        target: InstanceAtomId,
        expected_target_instance: MoleculeInstanceId,
    },
    /// A mapped bond's source molecule instance has no instance mapping.
    MissingInstanceMappingForBond { source: InstanceBondId },
    /// A mapped bond belongs to the wrong target molecule instance.
    InconsistentBondInstance {
        source: InstanceBondId,
        target: InstanceBondId,
        expected_target_instance: MoleculeInstanceId,
    },
    /// A mapped bond has a source endpoint missing from the atom mapping.
    MissingAtomMappingForBondEndpoint {
        source_bond: InstanceBondId,
        source_atom: InstanceAtomId,
    },
    /// A mapped target bond does not connect the mapped source endpoints.
    InconsistentBondEndpoints {
        source: InstanceBondId,
        target: InstanceBondId,
    },
    /// A topology edit result was paired with a topology other than the
    /// mapping target.
    TargetTopologyMismatch,
}

impl fmt::Display for TopologyMappingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotSameLayout => formatter.write_str("topologies do not have the same layout"),
            Self::InvalidSource(kind) => write!(formatter, "invalid source {kind} mapping"),
            Self::InvalidTarget(kind) => write!(formatter, "invalid target {kind} mapping"),
            Self::DuplicateSource(kind) => write!(formatter, "duplicate source {kind} mapping"),
            Self::DuplicateTarget(kind) => write!(formatter, "duplicate target {kind} mapping"),
            Self::MissingDefinitionMappingForInstance { source, definition } => write!(
                formatter,
                "mapped source instance {source} has no mapping for its definition {definition}"
            ),
            Self::InconsistentInstanceDefinition {
                source,
                target,
                expected_target_definition,
                actual_target_definition,
            } => write!(
                formatter,
                "mapped instance {source} -> {target} requires target definition \
                 {expected_target_definition}, but the target instance uses \
                 {actual_target_definition}"
            ),
            Self::MissingInstanceMappingForAtom { source } => write!(
                formatter,
                "mapped source atom {source} has no mapping for its molecule instance"
            ),
            Self::InconsistentAtomInstance {
                source,
                target,
                expected_target_instance,
            } => write!(
                formatter,
                "mapped atom {source} -> {target} must target molecule instance \
                 {expected_target_instance}"
            ),
            Self::MissingInstanceMappingForBond { source } => write!(
                formatter,
                "mapped source bond {source} has no mapping for its molecule instance"
            ),
            Self::InconsistentBondInstance {
                source,
                target,
                expected_target_instance,
            } => write!(
                formatter,
                "mapped bond {source} -> {target} must target molecule instance \
                 {expected_target_instance}"
            ),
            Self::MissingAtomMappingForBondEndpoint {
                source_bond,
                source_atom,
            } => write!(
                formatter,
                "mapped source bond {source_bond} has an unmapped endpoint {source_atom}"
            ),
            Self::InconsistentBondEndpoints { source, target } => write!(
                formatter,
                "mapped bond {source} -> {target} does not agree with the atom mapping"
            ),
            Self::TargetTopologyMismatch => {
                formatter.write_str("topology edit result does not match the mapping target")
            }
        }
    }
}

impl std::error::Error for TopologyMappingError {}

/// Result of an immutable topology-changing operation.
#[derive(Debug, Clone)]
pub struct TopologyEditResult {
    topology: Topology,
    mapping: TopologyMapping,
}

impl TopologyEditResult {
    /// Constructs a topology edit result after checking target identity.
    pub fn new(topology: Topology, mapping: TopologyMapping) -> Result<Self, TopologyMappingError> {
        if &topology.identity() != mapping.new_identity() {
            return Err(TopologyMappingError::TargetTopologyMismatch);
        }
        Ok(Self { topology, mapping })
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn mapping(&self) -> &TopologyMapping {
        &self.mapping
    }

    pub fn into_parts(self) -> (Topology, TopologyMapping) {
        (self.topology, self.mapping)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bio::SmcraAtomSiteMetadata;
    use crate::core::{BondOrder, Conformer};
    use crate::geometry::Point3;
    use crate::query;
    use crate::substructure;
    use crate::units::{Quantity, ANGSTROM};

    fn tombstoned_molecule() -> (SmallMolecule, AtomId, AtomId, BondId) {
        let mut graph = Molecule::new();
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
        (SmallMolecule::from_graph(graph), carbon, oxygen, bond)
    }

    fn topology_with_reused_definition() -> (Topology, AtomId, AtomId, BondId) {
        let (molecule, carbon, oxygen, bond) = tombstoned_molecule();
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&molecule).unwrap();
        let mut first = MoleculeInstanceMetadata::default();
        first.insert_role(MoleculeRole::Ligand);
        builder.add_instance(definition, first).unwrap();
        let mut second = MoleculeInstanceMetadata::default();
        second.insert_role(MoleculeRole::Solvent);
        builder.add_instance(definition, second).unwrap();
        (builder.build().unwrap(), carbon, oxygen, bond)
    }

    fn topology_from_smiles(
        smiles: &str,
    ) -> (
        Topology,
        MoleculeDefinitionId,
        MoleculeInstanceId,
        Vec<AtomId>,
        Vec<BondId>,
    ) {
        let molecule = SmallMolecule::from_smiles_sanitized(smiles).unwrap();
        let atoms = molecule.graph().atom_ids().collect();
        let bonds = molecule.graph().bond_ids().collect();
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&molecule).unwrap();
        let instance = builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        (builder.build().unwrap(), definition, instance, atoms, bonds)
    }

    fn topology_with_two_distinct_definitions(reverse: bool) -> Topology {
        let carbon = SmallMolecule::from_smiles_sanitized("C").unwrap();
        let carbon_oxygen = SmallMolecule::from_smiles_sanitized("CO").unwrap();
        let molecules = if reverse {
            [&carbon_oxygen, &carbon]
        } else {
            [&carbon, &carbon_oxygen]
        };
        let mut builder = TopologyBuilder::new();
        for molecule in molecules {
            let definition = builder.add_small_molecule_definition(molecule).unwrap();
            builder
                .add_instance(definition, MoleculeInstanceMetadata::default())
                .unwrap();
        }
        builder.build().unwrap()
    }

    fn topology_with_repeated_identical_definitions() -> Topology {
        let water = SmallMolecule::from_smiles_sanitized("O").unwrap();
        let mut builder = TopologyBuilder::new();
        for _ in 0..2 {
            let definition = builder.add_small_molecule_definition(&water).unwrap();
            for _ in 0..2 {
                builder
                    .add_instance(definition, MoleculeInstanceMetadata::default())
                    .unwrap();
            }
        }
        builder.build().unwrap()
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
    fn topology_identity_is_exact_while_layout_can_match() {
        let (topology, ..) = topology_with_reused_definition();
        let clone = topology.clone();
        let (independent, ..) = topology_with_reused_definition();

        assert!(topology.same_identity(&clone));
        assert!(topology.same_layout(&clone));
        assert!(!topology.same_identity(&independent));
        assert!(topology.same_layout(&independent));
        assert_eq!(clone.atom_ids(), topology.atom_ids());
    }

    #[test]
    fn layout_equality_does_not_reorder_definitions_instances_or_dense_state() {
        let forward = topology_with_two_distinct_definitions(false);
        let reverse = topology_with_two_distinct_definitions(true);

        assert!(!forward.same_identity(&reverse));
        assert!(!forward.same_layout(&reverse));
        assert_eq!(
            TopologyMapping::between_identical_layouts(&forward, &reverse),
            Err(TopologyMappingError::NotSameLayout)
        );
    }

    #[test]
    fn identical_layout_mapping_preserves_repeated_definitions_and_instances() {
        let source = topology_with_repeated_identical_definitions();
        let target = topology_with_repeated_identical_definitions();

        assert!(source.same_layout(&target));
        assert!(!source.same_identity(&target));
        let mapping = TopologyMapping::between_identical_layouts(&source, &target).unwrap();
        for (definition, _) in source.definitions() {
            assert_eq!(mapping.map_definition(definition), Some(definition));
        }
        for (instance, _) in source.instances() {
            assert_eq!(mapping.map_instance(instance), Some(instance));
        }
        assert!(mapping.removed_definitions().is_empty());
        assert!(mapping.added_definitions().is_empty());
        assert!(mapping.removed_instances().is_empty());
        assert!(mapping.added_instances().is_empty());
    }

    #[test]
    fn builder_is_transactional_and_does_not_intern_equal_definitions() {
        let (molecule, ..) = tombstoned_molecule();
        let mut conformer = Conformer::new(ANGSTROM).unwrap();
        for atom in molecule.graph().atom_ids() {
            conformer
                .set_position(atom, Quantity::new(Point3::origin(), ANGSTROM))
                .unwrap();
        }
        let mut molecule = molecule;
        molecule.graph_mut().add_conformer(conformer).unwrap();

        let mut builder = TopologyBuilder::new();
        let first = builder.add_small_molecule_definition(&molecule).unwrap();
        assert_eq!(
            builder.definitions[first.index()]
                .graph()
                .conformers()
                .count(),
            0
        );
        assert_eq!(molecule.graph().conformers().count(), 1);
        assert_eq!(
            builder.add_instance(
                MoleculeDefinitionId::new(99),
                MoleculeInstanceMetadata::default()
            ),
            Err(TopologyBuildError::InvalidMoleculeDefinitionId(
                MoleculeDefinitionId::new(99)
            ))
        );
        assert!(builder.instances.is_empty());
        builder
            .add_instance(first, MoleculeInstanceMetadata::default())
            .unwrap();
        let second = builder.add_small_molecule_definition(&molecule).unwrap();
        builder
            .add_instance(second, MoleculeInstanceMetadata::default())
            .unwrap();
        let topology = builder.build().unwrap();
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
    fn selections_distinguish_instances_components_elements_and_queries() {
        let molecule = SmallMolecule::from_smiles_sanitized("CC.O").unwrap();
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&molecule).unwrap();
        let mut metadata = MoleculeInstanceMetadata::default();
        metadata.insert_role(MoleculeRole::Ligand);
        let instance = builder.add_instance(definition, metadata).unwrap();
        let topology = builder.build().unwrap();

        let ligand = AtomSelection::for_roles(&topology, [MoleculeRole::Ligand]).unwrap();
        assert_eq!(ligand.indices().len(), 3);
        let first_component = AtomSelection::connected_component(&topology, instance, 0).unwrap();
        let second_component = AtomSelection::connected_component(&topology, instance, 1).unwrap();
        assert_eq!(first_component.indices().len(), 2);
        assert_eq!(second_component.indices().len(), 1);
        let oxygen =
            AtomSelection::for_elements(&topology, [Element::from_symbol("O").unwrap()]).unwrap();
        assert_eq!(oxygen.indices().len(), 1);

        let query = query::parse_smarts("O").unwrap();
        let matches = substructure::find_substructure_matches(molecule.graph(), &query).unwrap();
        let from_query = AtomSelection::from_query_matches(&topology, instance, &matches).unwrap();
        assert_eq!(
            from_query.semantic_ids(&topology).unwrap(),
            oxygen.semantic_ids(&topology).unwrap()
        );

        let mut independent_builder = TopologyBuilder::new();
        let definition = independent_builder
            .add_small_molecule_definition(&molecule)
            .unwrap();
        independent_builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let independent = independent_builder.build().unwrap();
        assert_eq!(
            ligand.ensure_compatible(&independent),
            Err(SelectionError::TopologyIdentityMismatch)
        );
    }

    #[test]
    fn hierarchy_label_selection_uses_instance_qualified_atoms() {
        let mut macro_builder = MacroMolecule::builder();
        let atom = macro_builder
            .graph_mut()
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");
        let chain = macro_builder.hierarchy_mut().add_chain("A", None).unwrap();
        let residue = macro_builder
            .hierarchy_mut()
            .add_residue(chain, "GLY", Some(1), None, None)
            .unwrap();
        macro_builder
            .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
            .unwrap();
        let molecule = macro_builder.build().unwrap();

        let mut builder = TopologyBuilder::new();
        let definition = builder.add_macro_molecule_definition(&molecule).unwrap();
        let first = builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let small = SmallMolecule::from_smiles_sanitized("O").unwrap();
        let small_definition = builder.add_small_molecule_definition(&small).unwrap();
        builder
            .add_instance(small_definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let second = builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let topology = builder.build().unwrap();
        assert_eq!(topology.definition_count(), 2);
        assert!(topology
            .definition(definition)
            .unwrap()
            .macro_molecule()
            .is_some());
        assert!(topology
            .definition(small_definition)
            .unwrap()
            .small_molecule()
            .is_some());
        assert_eq!(
            AtomSelection::for_chain_label(&topology, "A")
                .unwrap()
                .semantic_ids(&topology)
                .unwrap(),
            vec![
                InstanceAtomId::new(first, atom),
                InstanceAtomId::new(second, atom)
            ]
        );
    }

    #[test]
    fn topology_mapping_reports_retained_added_and_removed_state() {
        let (old_molecule, carbon, oxygen, bond) = tombstoned_molecule();
        let mut old_builder = TopologyBuilder::new();
        let old_definition = old_builder
            .add_small_molecule_definition(&old_molecule)
            .unwrap();
        let old_instance = old_builder
            .add_instance(old_definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let old_extra_molecule = SmallMolecule::from_smiles_sanitized("O").unwrap();
        let old_extra_definition = old_builder
            .add_small_molecule_definition(&old_extra_molecule)
            .unwrap();
        let old_extra_instance = old_builder
            .add_instance(old_extra_definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let old = old_builder.build().unwrap();

        let mut new_graph = Molecule::new();
        let new_carbon = new_graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");
        let new_molecule = SmallMolecule::from_graph(new_graph);
        let mut new_builder = TopologyBuilder::new();
        let new_definition = new_builder
            .add_small_molecule_definition(&new_molecule)
            .unwrap();
        let new_instance = new_builder
            .add_instance(new_definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let new_extra_molecule = SmallMolecule::from_smiles_sanitized("N").unwrap();
        let new_extra_definition = new_builder
            .add_small_molecule_definition(&new_extra_molecule)
            .unwrap();
        let new_extra_instance = new_builder
            .add_instance(new_extra_definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let new = new_builder.build().unwrap();

        let old_carbon = InstanceAtomId::new(old_instance, carbon);
        let mapping = TopologyMapping::from_pairs(
            &old,
            &new,
            [(old_definition, new_definition)],
            [(old_instance, new_instance)],
            [(old_carbon, InstanceAtomId::new(new_instance, new_carbon))],
            [],
        )
        .unwrap();
        assert_eq!(
            mapping.map_atom(old_carbon),
            Some(InstanceAtomId::new(new_instance, new_carbon))
        );
        assert_eq!(
            mapping.map_atom_index(old.atom_index(old_carbon).unwrap()),
            new.atom_index(InstanceAtomId::new(new_instance, new_carbon))
        );
        assert_eq!(
            mapping.removed_atoms(),
            &[
                InstanceAtomId::new(old_instance, oxygen),
                *old.atom_ids()
                    .iter()
                    .find(|atom| atom.molecule() == old_extra_instance)
                    .unwrap(),
            ]
        );
        assert_eq!(
            mapping.removed_bonds(),
            &[InstanceBondId::new(old_instance, bond)]
        );
        assert_eq!(mapping.removed_definitions(), &[old_extra_definition]);
        assert_eq!(mapping.added_definitions(), &[new_extra_definition]);
        assert_eq!(mapping.removed_instances(), &[old_extra_instance]);
        assert_eq!(mapping.added_instances(), &[new_extra_instance]);
        assert_eq!(
            mapping.added_atoms(),
            &[*new
                .atom_ids()
                .iter()
                .find(|atom| atom.molecule() == new_extra_instance)
                .unwrap()]
        );
        assert!(mapping.added_bonds().is_empty());

        let equivalent = topology_with_reused_definition().0;
        let independently_equivalent = topology_with_reused_definition().0;
        let round_trip =
            TopologyMapping::between_identical_layouts(&equivalent, &independently_equivalent)
                .unwrap();
        for atom in equivalent.atom_ids() {
            assert_eq!(round_trip.map_atom(*atom), Some(*atom));
        }
    }

    #[test]
    fn topology_mapping_rejects_duplicate_and_cross_instance_atom_targets() {
        let (source, carbon, _, bond) = topology_with_reused_definition();
        let (target, ..) = topology_with_reused_definition();
        let definition = MoleculeDefinitionId::new(0);
        let first = MoleculeInstanceId::new(0);
        let second = MoleculeInstanceId::new(1);
        let first_carbon = InstanceAtomId::new(first, carbon);
        let second_carbon = InstanceAtomId::new(second, carbon);

        assert_eq!(
            TopologyMapping::from_pairs(
                &source,
                &target,
                [(definition, definition)],
                [(first, first), (second, second)],
                [(first_carbon, first_carbon), (second_carbon, first_carbon),],
                [],
            ),
            Err(TopologyMappingError::DuplicateTarget(MappingKind::Atom))
        );
        assert_eq!(
            TopologyMapping::from_pairs(
                &source,
                &target,
                [(definition, definition)],
                [(first, first), (second, second)],
                [(first_carbon, second_carbon)],
                [],
            ),
            Err(TopologyMappingError::InconsistentAtomInstance {
                source: first_carbon,
                target: second_carbon,
                expected_target_instance: first,
            })
        );
        let first_bond = InstanceBondId::new(first, bond);
        let second_bond = InstanceBondId::new(second, bond);
        assert_eq!(
            TopologyMapping::from_pairs(
                &source,
                &target,
                [(definition, definition)],
                [(first, first), (second, second)],
                [],
                [(first_bond, second_bond)],
            ),
            Err(TopologyMappingError::InconsistentBondInstance {
                source: first_bond,
                target: second_bond,
                expected_target_instance: first,
            })
        );
    }

    #[test]
    fn topology_mapping_rejects_inconsistent_definition_and_bond_relationships() {
        let source_layout = topology_with_two_distinct_definitions(false);
        let target_layout = topology_with_two_distinct_definitions(false);
        assert_eq!(
            TopologyMapping::from_pairs(
                &source_layout,
                &target_layout,
                [(MoleculeDefinitionId::new(0), MoleculeDefinitionId::new(1),)],
                [(MoleculeInstanceId::new(0), MoleculeInstanceId::new(0))],
                [],
                [],
            ),
            Err(TopologyMappingError::InconsistentInstanceDefinition {
                source: MoleculeInstanceId::new(0),
                target: MoleculeInstanceId::new(0),
                expected_target_definition: MoleculeDefinitionId::new(1),
                actual_target_definition: MoleculeDefinitionId::new(0),
            })
        );

        let (source, source_definition, source_instance, source_atoms, source_bonds) =
            topology_from_smiles("CCC");
        let (target, target_definition, target_instance, target_atoms, target_bonds) =
            topology_from_smiles("CCC");
        let source_bond = source_bonds[0];
        let (source_a, source_b) = source
            .definition(source_definition)
            .unwrap()
            .graph()
            .bond(source_bond)
            .unwrap()
            .endpoints();
        let source_other = source_atoms
            .iter()
            .copied()
            .find(|atom| *atom != source_a && *atom != source_b)
            .unwrap();
        let target_other = target_atoms
            .iter()
            .copied()
            .find(|atom| *atom != source_a && *atom != source_b)
            .unwrap();
        let atom_pairs = [
            (
                InstanceAtomId::new(source_instance, source_a),
                InstanceAtomId::new(target_instance, source_a),
            ),
            (
                InstanceAtomId::new(source_instance, source_b),
                InstanceAtomId::new(target_instance, target_other),
            ),
            (
                InstanceAtomId::new(source_instance, source_other),
                InstanceAtomId::new(target_instance, source_b),
            ),
        ];
        let source_bond = InstanceBondId::new(source_instance, source_bond);
        let target_bond = InstanceBondId::new(target_instance, target_bonds[0]);
        assert_eq!(
            TopologyMapping::from_pairs(
                &source,
                &target,
                [(source_definition, target_definition)],
                [(source_instance, target_instance)],
                atom_pairs,
                [(source_bond, target_bond)],
            ),
            Err(TopologyMappingError::InconsistentBondEndpoints {
                source: source_bond,
                target: target_bond,
            })
        );
    }

    #[test]
    fn topology_edit_result_checks_mapping_target_identity() {
        let source = topology_with_repeated_identical_definitions();
        let target = topology_with_repeated_identical_definitions();
        let wrong_target = topology_with_repeated_identical_definitions();
        let mapping = TopologyMapping::between_identical_layouts(&source, &target).unwrap();

        assert_eq!(
            TopologyEditResult::new(wrong_target, mapping).unwrap_err(),
            TopologyMappingError::TargetTopologyMismatch
        );
        let mapping = TopologyMapping::between_identical_layouts(&source, &target).unwrap();
        let result = TopologyEditResult::new(target.clone(), mapping).unwrap();
        assert!(result.topology().same_identity(&target));
    }
}
