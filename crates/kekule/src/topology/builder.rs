use std::collections::BTreeMap;
use std::fmt;

use crate::core::Molecule;
use crate::properties::{Properties, PropertyError, PropertyKey, PropertyTableMut, PropertyValue};

use super::{
    AtomSiteId, ChainId, Hierarchy, InstanceAtomId, MoleculeClass, MoleculeDefinition,
    MoleculeDefinitionId, MoleculeInstance, MoleculeInstanceId, ResidueClass, ResidueId, Topology,
    TopologyAtomIndex, TopologyBondIndex,
};

/// Linear, validate-then-commit builder for coordinate-free topology.
///
/// Add connected [`Molecule`] values as reusable definitions, then add one or
/// more explicit instances of each definition. The builder stages the one
/// system-level hierarchy and validates all hierarchy references against final
/// instance-qualified atom IDs when [`Self::build`] publishes the topology.
///
/// Use [`Self::add_molecule`] when definition reuse is unimportant. Use
/// [`Self::add_molecule_definition`] followed by [`Self::add_instance`] when
/// several occurrences should share one definition.
///
/// # Example
///
/// ```
/// use kekule::{smiles, topology::TopologyBuilder};
///
/// let water = smiles::to_molecules("O")?.pop().unwrap();
/// let mut builder = TopologyBuilder::new();
/// let definition = builder.add_molecule_definition(&water)?;
/// builder.add_instance(definition)?;
/// builder.add_instance(definition)?;
/// let topology = builder.build()?;
///
/// assert_eq!(topology.definition_count(), 1);
/// assert_eq!(topology.instance_count(), 2);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TopologyBuilder {
    definitions: Vec<MoleculeDefinition>,
    pub(super) instances: Vec<MoleculeInstance>,
    hierarchy: Hierarchy,
    properties: Properties,
    molecule_class_overrides: BTreeMap<MoleculeDefinitionId, MoleculeClass>,
    residue_class_overrides: BTreeMap<ResidueId, ResidueClass>,
}

impl TopologyBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub(super) fn from_topology(topology: Topology) -> Self {
        let Topology {
            definitions,
            instances,
            hierarchy,
            properties,
            ..
        } = topology;
        let molecule_class_overrides = definitions
            .iter()
            .map(|definition| (definition.id(), definition.class()))
            .collect();
        let residue_class_overrides = hierarchy
            .residues()
            .map(|(id, residue)| (id, residue.class()))
            .collect();
        Self {
            definitions,
            instances,
            hierarchy,
            properties,
            molecule_class_overrides,
            residue_class_overrides,
        }
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

    pub fn molecule_instance_properties_mut(&mut self) -> PropertyTableMut<'_> {
        self.sync_property_dimensions();
        PropertyTableMut::new(self.properties.molecule_instances_mut())
    }

    pub fn atom_properties_mut(&mut self) -> PropertyTableMut<'_> {
        self.sync_property_dimensions();
        PropertyTableMut::new(self.properties.atoms_mut())
    }

    pub fn bond_properties_mut(&mut self) -> PropertyTableMut<'_> {
        self.sync_property_dimensions();
        PropertyTableMut::new(self.properties.bonds_mut())
    }

    pub fn chain_properties_mut(&mut self) -> PropertyTableMut<'_> {
        self.sync_property_dimensions();
        PropertyTableMut::new(self.properties.chains_mut())
    }

    pub fn residue_properties_mut(&mut self) -> PropertyTableMut<'_> {
        self.sync_property_dimensions();
        PropertyTableMut::new(self.properties.residues_mut())
    }

    pub fn atom_site_properties_mut(&mut self) -> PropertyTableMut<'_> {
        self.sync_property_dimensions();
        PropertyTableMut::new(self.properties.atom_sites_mut())
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

    /// Overrides automatic classification for one staged molecule definition.
    ///
    /// Explicit assignments have precedence over inference performed by
    /// [`Self::build`]. The class remains definition-scoped and is shared by
    /// every instance of the definition.
    pub fn set_molecule_class(
        &mut self,
        definition: MoleculeDefinitionId,
        class: MoleculeClass,
    ) -> Result<(), TopologyBuildError> {
        self.definition(definition)?;
        self.molecule_class_overrides.insert(definition, class);
        Ok(())
    }

    /// Overrides automatic classification for one staged hierarchy residue.
    pub fn set_residue_class(
        &mut self,
        residue: ResidueId,
        class: ResidueClass,
    ) -> Result<(), TopologyBuildError> {
        self.hierarchy
            .residue(residue)
            .map_err(|_| TopologyBuildError::InvalidResidueId(residue))?;
        self.residue_class_overrides.insert(residue, class);
        Ok(())
    }

    pub fn add_molecule_definition(
        &mut self,
        molecule: &Molecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        self.commit_definition(molecule.clone())
    }

    pub fn add_molecule_definition_owned(
        &mut self,
        molecule: Molecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
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

        super::classification::finalize(
            &mut self.definitions,
            &self.instances,
            &mut self.hierarchy,
            &self.molecule_class_overrides,
            &self.residue_class_overrides,
        );

        self.properties.resize_domains(
            self.instances.len(),
            atom_count,
            bond_count,
            self.hierarchy.chains().count(),
            self.hierarchy.residues().count(),
            self.hierarchy.atom_sites().count(),
        );
        self.properties
            .validate_topology_dimensions([
                self.instances.len(),
                atom_count,
                bond_count,
                self.hierarchy.chains().count(),
                self.hierarchy.residues().count(),
                self.hierarchy.atom_sites().count(),
            ])
            .map_err(|error| TopologyBuildError::Property(Box::new(error)))?;

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
        self.definitions.push(MoleculeDefinition {
            id,
            molecule,
            class: MoleculeClass::SmallMolecule,
        });
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
            class: MoleculeClass::SmallMolecule,
        });
        self.instances.push(MoleculeInstance {
            id: instance,
            definition,
        });
        Ok((definition, instance))
    }
}

pub(super) trait FromRawId {
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

pub(super) fn checked_id<T: FromRawId>(
    length: usize,
    kind: TopologyIdKind,
) -> Result<T, TopologyBuildError> {
    let raw = crate::core::checked_raw_id(length)
        .map_err(|_| TopologyBuildError::IdentifierCapacityExceeded(kind))?;
    Ok(T::from_raw(raw))
}

pub(super) fn checked_future_len(
    current: usize,
    additional: usize,
    kind: TopologyIdKind,
) -> Result<(), TopologyBuildError> {
    crate::core::checked_fixed_id_collection_len(current, additional)
        .map_err(|_| TopologyBuildError::IdentifierCapacityExceeded(kind))
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
    Property(Box<PropertyError>),
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    InvalidResidueId(ResidueId),
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
            Self::Property(error) => write!(formatter, "invalid topology properties: {error}"),
            Self::InvalidMoleculeDefinitionId(id) => {
                write!(formatter, "invalid molecule definition: {id}")
            }
            Self::InvalidResidueId(id) => write!(formatter, "invalid topology residue: {id}"),
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

impl std::error::Error for TopologyBuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Property(error) => Some(error),
            Self::InvalidHierarchy(error) => Some(error),
            _ => None,
        }
    }
}
