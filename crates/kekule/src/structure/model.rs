use std::fmt;
use std::sync::Arc;

use crate::core::{Atom, Bond, Molecule};
use crate::geometry::{PeriodicCell, Point3};
use crate::properties::{
    Properties, PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue,
};
use crate::topology::transform::TopologySubsetError;
use crate::topology::{
    AtomSelection, AtomSiteId, AtomSiteView, ChainId, ChainView, Hierarchy, InstanceAtomId,
    InstanceBondId, MoleculeClass, MoleculeDefinitionId, MoleculeInstance, MoleculeInstanceId,
    ResidueClass, ResidueId, ResidueView, Topology, TopologyAtomIndex, TopologyBuildError,
    TopologyBuilder, TopologyError, TopologyPerceptionError,
};
use crate::units::Quantity;

use super::{PositionError, Positions};

/// One concrete realization of one immutable topology.
///
/// A model owns a shared [`Topology`], one position for every topology atom, an
/// optional [`PeriodicCell`], and realization-scoped properties. The topology
/// remains coordinate-free and may be shared with other models, ensemble
/// members, trajectory frames, selections, and prepared modelling systems.
///
/// Use [`Self::from_molecule`] for a single connected molecule,
/// [`Self::builder`] to assemble several molecules and their positions together,
/// or [`Self::new`] when an existing topology and correctly ordered dense
/// positions are already available.
#[derive(Debug, Clone)]
pub struct Model {
    pub(super) topology: Arc<Topology>,
    pub(super) positions: Positions,
    pub(super) cell: Option<PeriodicCell>,
    pub(super) properties: Properties,
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.positions == other.positions
            && self.cell == other.cell
            && self.properties == other.properties
    }
}

impl Model {
    /// Creates a non-periodic model with empty atom and bond property tables.
    pub fn new(
        topology: impl Into<Arc<Topology>>,
        positions: Positions,
    ) -> Result<Self, ModelError> {
        let topology = topology.into();
        validate_positions_len(&topology, &positions)?;
        let properties = Properties::realization(topology.atom_count(), topology.bond_count());
        Ok(Self {
            topology,
            positions,
            cell: None,
            properties,
        })
    }

    /// Creates a model from complete geometry and realization properties.
    pub fn with_properties(
        topology: impl Into<Arc<Topology>>,
        positions: Positions,
        cell: Option<PeriodicCell>,
        properties: Properties,
    ) -> Result<Self, ModelError> {
        let topology = topology.into();
        validate_dimensions(&topology, &positions, &properties)?;
        properties.validate_realization_canonical_properties()?;
        Ok(Self {
            topology,
            positions,
            cell,
            properties,
        })
    }

    pub fn builder() -> ModelBuilder {
        ModelBuilder::new()
    }

    /// Builds a single-molecule model from dense positions in molecule atom order.
    ///
    /// Position construction performs unit conversion and finite-value
    /// validation before this topology-owning operation.
    pub fn from_molecule(
        molecule: &Molecule,
        positions: &Positions,
    ) -> Result<Self, ModelBuildError> {
        let mut builder = ModelBuilder::new();
        builder.add_molecule(molecule, positions)?;
        builder.build()
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Installs default perception through a new topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`], perceiving each definition once.
    /// Positions, the periodic cell, and all properties are retained without
    /// copying realization state. Failure leaves the entire model unchanged.
    /// Other owners keep their original topology. Existing selections and
    /// prepared calculations remain bound to that original snapshot.
    ///
    /// ```
    /// use kekule::{smiles, structure::{Model, Positions}};
    ///
    /// let topology = smiles::to_topology("c1ccccc1")?;
    /// let mut model = Model::new(topology, Positions::zeros(6))?;
    /// model.perceive()?;
    /// assert!(model.topology().molecules().all(|m| m.molecule().perception().has_rings()));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// Constructs an induced structural slice and transfers all dense model state.
    pub fn slice(&self, selection: &AtomSelection) -> Result<Self, ModelSliceError> {
        let subset = self.topology.subset(selection)?;
        let correspondence = subset.correspondence();
        let atom_indices = correspondence
            .source_atom_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let bond_indices = correspondence
            .source_bond_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let positions = self
            .positions
            .select_indices(&atom_indices)
            .map_err(ModelError::from)?;
        let properties = self
            .properties
            .project_realization(&atom_indices, &bond_indices)
            .map_err(ModelError::Property)?;
        Ok(Self::with_properties(
            subset.shared_topology(),
            positions,
            self.cell,
            properties,
        )?)
    }

    pub fn atom(&self, atom: InstanceAtomId) -> Result<&Atom, TopologyError> {
        self.topology.atom(atom)
    }

    pub fn bond(&self, bond: InstanceBondId) -> Result<&Bond, TopologyError> {
        self.topology.bond(bond)
    }

    pub fn atoms(&self) -> impl ExactSizeIterator<Item = (InstanceAtomId, &Atom)> {
        self.topology.atoms()
    }

    pub fn bonds(&self) -> impl ExactSizeIterator<Item = (InstanceBondId, &Bond)> {
        self.topology.bonds()
    }

    pub fn instances(
        &self,
    ) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &MoleculeInstance)> {
        self.topology.instances()
    }

    pub fn hierarchy(&self) -> &Hierarchy {
        self.topology.hierarchy()
    }

    pub fn chains(&self) -> impl Iterator<Item = ChainView<'_>> {
        self.topology.chains()
    }

    pub fn residues(&self) -> impl Iterator<Item = ResidueView<'_>> {
        self.topology.residues()
    }

    pub fn atom_sites(&self) -> impl Iterator<Item = AtomSiteView<'_>> {
        self.topology.atom_sites()
    }

    pub fn chain(&self, chain: ChainId) -> Result<ChainView<'_>, TopologyError> {
        self.topology.chain(chain)
    }

    pub fn residue(&self, residue: ResidueId) -> Result<ResidueView<'_>, TopologyError> {
        self.topology.residue(residue)
    }

    pub fn atom_site(&self, atom_site: AtomSiteId) -> Result<AtomSiteView<'_>, TopologyError> {
        self.topology.atom_site(atom_site)
    }

    pub fn atom_for_site(&self, atom_site: AtomSiteId) -> Result<InstanceAtomId, TopologyError> {
        self.topology.atom_for_site(atom_site)
    }

    pub fn atom_site_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<AtomSiteView<'_>>, TopologyError> {
        self.topology.atom_site_for_atom(atom)
    }

    pub fn residue_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<ResidueView<'_>>, TopologyError> {
        self.topology.residue_for_atom(atom)
    }

    pub fn chain_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<ChainView<'_>>, TopologyError> {
        self.topology.chain_for_atom(atom)
    }

    pub fn residue_for_site(
        &self,
        atom_site: AtomSiteId,
    ) -> Result<ResidueView<'_>, TopologyError> {
        self.topology.residue_for_site(atom_site)
    }

    pub fn chain_for_residue(&self, residue: ResidueId) -> Result<ChainView<'_>, TopologyError> {
        self.topology.chain_for_residue(residue)
    }

    pub const fn positions(&self) -> &Positions {
        &self.positions
    }

    pub fn position(&self, atom: InstanceAtomId) -> Result<Quantity<Point3>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.positions.position_at(index.index())?)
    }

    pub fn position_at(&self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position_at(index.index())
    }

    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.positions.set_position_at(index.index(), position)?)
    }

    pub fn set_positions<T>(&mut self, positions: Quantity<T>) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        self.positions.set_all(positions)
    }

    pub const fn cell(&self) -> Option<&PeriodicCell> {
        self.cell.as_ref()
    }

    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.cell = cell;
    }

    pub const fn properties(&self) -> &Properties {
        &self.properties
    }

    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, ModelError> {
        Ok(self.properties.insert(key, value)?)
    }

    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.properties.remove(key)
    }

    pub fn clear_properties(&mut self) {
        self.properties.clear_owner();
    }

    pub const fn atom_properties(&self) -> &PropertyTable {
        self.properties.realization_atom_properties()
    }

    pub const fn bond_properties(&self) -> &PropertyTable {
        self.properties.realization_bond_properties()
    }

    pub fn set_properties(&mut self, properties: Properties) -> Result<(), ModelError> {
        if properties.atoms().len() != self.topology.atom_count() {
            return Err(ModelError::AtomPropertyCountMismatch {
                expected: self.topology.atom_count(),
                actual: properties.atoms().len(),
            });
        }
        if properties.bonds().len() != self.topology.bond_count() {
            return Err(ModelError::BondPropertyCountMismatch {
                expected: self.topology.bond_count(),
                actual: properties.bonds().len(),
            });
        }
        properties.validate_realization_canonical_properties()?;
        self.properties = properties;
        Ok(())
    }

    pub fn occupancy(&self, atom: InstanceAtomId) -> Result<Option<f64>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.properties.occupancy_at(index.index())?)
    }

    pub fn b_factor(&self, atom: InstanceAtomId) -> Result<Option<Quantity<f64>>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.properties.b_factor_at(index.index())?)
    }

    pub fn set_occupancy(
        &mut self,
        atom: InstanceAtomId,
        value: Option<f64>,
    ) -> Result<(), ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        self.properties.set_occupancy_at(index.index(), value)?;
        Ok(())
    }

    pub fn set_b_factor(
        &mut self,
        atom: InstanceAtomId,
        value: Option<Quantity<f64>>,
    ) -> Result<(), ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        self.properties.set_b_factor_at(index.index(), value)?;
        Ok(())
    }

    pub(crate) fn install_canonical_atom_properties(
        &mut self,
        occupancies: Vec<Option<f64>>,
        b_factors: Vec<Option<f64>>,
    ) -> Result<(), ModelError> {
        Ok(self
            .properties
            .install_canonical_realization_atom_columns(occupancies, b_factors)?)
    }

    pub fn atom_property(
        &self,
        atom: InstanceAtomId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.atom_properties().value(key, index.index())?)
    }

    pub fn set_atom_property(
        &mut self,
        atom: InstanceAtomId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self
            .properties
            .set_realization_atom_value(key, index.index(), value)?)
    }

    pub fn insert_atom_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, ModelError> {
        Ok(self
            .properties
            .insert_realization_atom_column(key, column)?)
    }

    pub fn remove_atom_property_column(
        &mut self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, ModelError> {
        Ok(self.properties.remove_realization_atom_column(key)?)
    }

    pub fn bond_property(
        &self,
        bond: InstanceBondId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, ModelError> {
        let index = self
            .topology
            .bond_index(bond)
            .ok_or(ModelError::InvalidBondId(bond))?;
        Ok(self.bond_properties().value(key, index.index())?)
    }

    pub fn set_bond_property(
        &mut self,
        bond: InstanceBondId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelError> {
        let index = self
            .topology
            .bond_index(bond)
            .ok_or(ModelError::InvalidBondId(bond))?;
        Ok(self
            .properties
            .set_realization_bond_value(key, index.index(), value)?)
    }

    pub fn insert_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, ModelError> {
        Ok(self
            .properties
            .insert_realization_bond_column(key, column)?)
    }

    pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.properties.remove_realization_bond_column(key)
    }

    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    pub fn view(&self) -> ModelView<'_> {
        ModelView {
            topology: &self.topology,
            positions: &self.positions,
            cell: self.cell.as_ref(),
            properties: &self.properties,
        }
    }
}

fn validate_positions_len(topology: &Topology, positions: &Positions) -> Result<(), ModelError> {
    if positions.len() != topology.atom_count() {
        return Err(ModelError::PositionCountMismatch {
            expected: topology.atom_count(),
            actual: positions.len(),
        });
    }
    Ok(())
}

fn validate_dimensions(
    topology: &Topology,
    positions: &Positions,
    properties: &Properties,
) -> Result<(), ModelError> {
    validate_positions_len(topology, positions)?;
    if properties.atoms().len() != topology.atom_count() {
        return Err(ModelError::AtomPropertyCountMismatch {
            expected: topology.atom_count(),
            actual: properties.atoms().len(),
        });
    }
    if properties.bonds().len() != topology.bond_count() {
        return Err(ModelError::BondPropertyCountMismatch {
            expected: topology.bond_count(),
            actual: properties.bonds().len(),
        });
    }
    Ok(())
}

/// Borrowed topology, positions, cell, and realization properties for structural
/// kernels.
///
/// `ModelView` is the common coordinate-dependent input contract. It can borrow
/// an owned [`Model`], an ensemble member, or a trajectory frame without
/// allocating a temporary model. Views preserve exact topology identity.
#[derive(Debug, Clone, Copy)]
pub struct ModelView<'a> {
    topology: &'a Arc<Topology>,
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
    properties: &'a Properties,
}

impl<'a> ModelView<'a> {
    pub fn new(
        topology: &'a Arc<Topology>,
        positions: &'a Positions,
        cell: Option<&'a PeriodicCell>,
        properties: &'a Properties,
    ) -> Result<Self, ModelError> {
        validate_dimensions(topology, positions, properties)?;
        properties.validate_realization_canonical_properties()?;
        Ok(Self {
            topology,
            positions,
            cell,
            properties,
        })
    }

    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub(crate) const fn topology_arc(self) -> &'a Arc<Topology> {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    /// Materializes this validated borrowed realization as an owned model.
    pub fn to_model(self) -> Model {
        Model {
            topology: Arc::clone(self.topology),
            positions: self.positions.clone(),
            cell: self.cell.cloned(),
            properties: self.properties.clone(),
        }
    }

    pub fn atom(self, atom: InstanceAtomId) -> Result<&'a Atom, TopologyError> {
        self.topology.atom(atom)
    }

    pub fn bond(self, bond: InstanceBondId) -> Result<&'a Bond, TopologyError> {
        self.topology.bond(bond)
    }

    pub fn atoms(self) -> impl ExactSizeIterator<Item = (InstanceAtomId, &'a Atom)> + 'a {
        self.topology.atoms()
    }

    pub fn bonds(self) -> impl ExactSizeIterator<Item = (InstanceBondId, &'a Bond)> + 'a {
        self.topology.bonds()
    }

    pub fn instances(
        self,
    ) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &'a MoleculeInstance)> + 'a {
        self.topology.instances()
    }

    pub fn hierarchy(self) -> &'a Hierarchy {
        self.topology.hierarchy()
    }

    pub fn chains(self) -> impl Iterator<Item = ChainView<'a>> + 'a {
        self.topology.chains()
    }

    pub fn residues(self) -> impl Iterator<Item = ResidueView<'a>> + 'a {
        self.topology.residues()
    }

    pub fn atom_sites(self) -> impl Iterator<Item = AtomSiteView<'a>> + 'a {
        self.topology.atom_sites()
    }

    pub fn chain(self, chain: ChainId) -> Result<ChainView<'a>, TopologyError> {
        self.topology.chain(chain)
    }

    pub fn residue(self, residue: ResidueId) -> Result<ResidueView<'a>, TopologyError> {
        self.topology.residue(residue)
    }

    pub fn atom_site(self, atom_site: AtomSiteId) -> Result<AtomSiteView<'a>, TopologyError> {
        self.topology.atom_site(atom_site)
    }

    pub fn atom_for_site(self, atom_site: AtomSiteId) -> Result<InstanceAtomId, TopologyError> {
        self.topology.atom_for_site(atom_site)
    }

    pub fn atom_site_for_atom(
        self,
        atom: InstanceAtomId,
    ) -> Result<Option<AtomSiteView<'a>>, TopologyError> {
        self.topology.atom_site_for_atom(atom)
    }

    pub fn residue_for_atom(
        self,
        atom: InstanceAtomId,
    ) -> Result<Option<ResidueView<'a>>, TopologyError> {
        self.topology.residue_for_atom(atom)
    }

    pub fn chain_for_atom(
        self,
        atom: InstanceAtomId,
    ) -> Result<Option<ChainView<'a>>, TopologyError> {
        self.topology.chain_for_atom(atom)
    }

    pub fn residue_for_site(self, atom_site: AtomSiteId) -> Result<ResidueView<'a>, TopologyError> {
        self.topology.residue_for_site(atom_site)
    }

    pub fn chain_for_residue(self, residue: ResidueId) -> Result<ChainView<'a>, TopologyError> {
        self.topology.chain_for_residue(residue)
    }

    pub const fn positions(self) -> &'a Positions {
        self.positions
    }

    pub fn position(self, atom: InstanceAtomId) -> Result<Quantity<Point3>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.positions.position_at(index.index())?)
    }

    pub fn position_at(self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position_at(index.index())
    }

    pub const fn cell(self) -> Option<&'a PeriodicCell> {
        self.cell
    }

    pub const fn properties(self) -> &'a Properties {
        self.properties
    }

    pub const fn atom_properties(self) -> &'a PropertyTable {
        self.properties.realization_atom_properties()
    }

    pub const fn bond_properties(self) -> &'a PropertyTable {
        self.properties.realization_bond_properties()
    }

    pub fn atom_property(
        self,
        atom: InstanceAtomId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.atom_properties().value(key, index.index())?)
    }

    pub fn bond_property(
        self,
        bond: InstanceBondId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, ModelError> {
        let index = self
            .topology
            .bond_index(bond)
            .ok_or(ModelError::InvalidBondId(bond))?;
        Ok(self.bond_properties().value(key, index.index())?)
    }

    pub fn occupancy(self, atom: InstanceAtomId) -> Result<Option<f64>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.properties.occupancy_at(index.index())?)
    }

    pub fn b_factor(self, atom: InstanceAtomId) -> Result<Option<Quantity<f64>>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.properties.b_factor_at(index.index())?)
    }

    pub fn atom_count(self) -> usize {
        self.topology.atom_count()
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelError {
    PositionCountMismatch { expected: usize, actual: usize },
    AtomPropertyCountMismatch { expected: usize, actual: usize },
    BondPropertyCountMismatch { expected: usize, actual: usize },
    InvalidAtomId(InstanceAtomId),
    InvalidBondId(InstanceBondId),
    Position(PositionError),
    Property(PropertyError),
}

impl fmt::Display for ModelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PositionCountMismatch { expected, actual } => write!(
                formatter,
                "topology requires {expected} positions, but received {actual}"
            ),
            Self::AtomPropertyCountMismatch { expected, actual } => write!(
                formatter,
                "topology requires atom properties of length {expected}, but received {actual}"
            ),
            Self::BondPropertyCountMismatch { expected, actual } => write!(
                formatter,
                "topology requires bond properties of length {expected}, but received {actual}"
            ),
            Self::InvalidAtomId(atom) => write!(formatter, "invalid topology atom: {atom}"),
            Self::InvalidBondId(bond) => write!(formatter, "invalid topology bond: {bond}"),
            Self::Position(error) => write!(formatter, "invalid position data: {error}"),
            Self::Property(error) => write!(formatter, "invalid model property: {error}"),
        }
    }
}

impl std::error::Error for ModelError {}

/// Failure to subset topology or transfer one model's dense state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelSliceError {
    Topology(TopologySubsetError),
    Model(Box<ModelError>),
}

impl fmt::Display for ModelSliceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Topology(error) => write!(formatter, "cannot subset model topology: {error}"),
            Self::Model(error) => write!(formatter, "cannot transfer model state: {error}"),
        }
    }
}

impl std::error::Error for ModelSliceError {}

impl From<TopologySubsetError> for ModelSliceError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Topology(error)
    }
}

impl From<ModelError> for ModelSliceError {
    fn from(error: ModelError) -> Self {
        Self::Model(Box::new(error))
    }
}

impl From<PositionError> for ModelError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<PropertyError> for ModelError {
    fn from(error: PropertyError) -> Self {
        Self::Property(error)
    }
}

/// Convenience builder that assembles topology and one complete model.
///
/// Each call to [`Self::add_molecule`] adds one explicit molecule occurrence and
/// its dense positions. Geometry is supplied explicitly and remains outside the
/// molecule. [`Self::build`] publishes both the coordinate-free topology and the
/// matching model transactionally.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ModelBuilder {
    topology: TopologyBuilder,
    positions: Vec<Point3>,
}

impl ModelBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn topology_builder(&self) -> &TopologyBuilder {
        &self.topology
    }

    /// Returns mutable coordinate-free staging state, including hierarchy.
    ///
    /// This does not update the staged positions. Add molecule instances through
    /// [`Self::add_molecule`] or [`Self::add_instance`] to stage their coordinates
    /// together. [`Self::build`] returns an error if the final topology's atom
    /// count does not match the staged positions.
    pub fn topology_builder_mut(&mut self) -> &mut TopologyBuilder {
        &mut self.topology
    }

    pub fn add_molecule_definition(
        &mut self,
        molecule: &Molecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self.topology.add_molecule_definition(molecule)?)
    }

    /// Overrides automatic classification for one staged molecule definition.
    pub fn set_molecule_class(
        &mut self,
        definition: MoleculeDefinitionId,
        class: MoleculeClass,
    ) -> Result<(), ModelBuildError> {
        Ok(self.topology.set_molecule_class(definition, class)?)
    }

    /// Overrides automatic classification for one staged hierarchy residue.
    pub fn set_residue_class(
        &mut self,
        residue: ResidueId,
        class: ResidueClass,
    ) -> Result<(), ModelBuildError> {
        Ok(self.topology.set_residue_class(residue, class)?)
    }

    /// Adds an instance with dense positions in its definition's atom order.
    pub fn add_instance(
        &mut self,
        definition: MoleculeDefinitionId,
        positions: &Positions,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        let expected = self
            .topology
            .definition(definition)?
            .molecule()
            .atom_count();
        validate_position_count(expected, positions.len())?;
        self.positions
            .try_reserve(positions.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let instance = self.topology.add_instance(definition)?;
        self.positions.extend_from_slice(positions.values().value());
        Ok(instance)
    }

    /// Adds a molecule and dense positions in that molecule's atom order.
    pub fn add_molecule(
        &mut self,
        molecule: &Molecule,
        positions: &Positions,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        if molecule.atom_count() == 0 {
            return Err(ModelBuildError::Topology(
                TopologyBuildError::EmptyMoleculeDefinition,
            ));
        }
        validate_position_count(molecule.atom_count(), positions.len())?;
        self.positions
            .try_reserve(positions.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let instance = self.topology.add_molecule(molecule)?;
        self.positions.extend_from_slice(positions.values().value());
        Ok(instance)
    }

    /// Validates and publishes the staged topology and model.
    ///
    /// Returns an error if topology validation fails or the final topology's
    /// atom count differs from the number of staged positions, including after
    /// edits through [`Self::topology_builder_mut`].
    pub fn build(self) -> Result<Model, ModelBuildError> {
        let topology = Arc::new(self.topology.build()?);
        let positions = Positions::from_canonical_values(self.positions);
        Model::new(topology, positions).map_err(ModelBuildError::from)
    }
}

fn validate_position_count(expected: usize, actual: usize) -> Result<(), ModelBuildError> {
    if actual != expected {
        return Err(ModelBuildError::InstancePositionCountMismatch { expected, actual });
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelBuildError {
    InstancePositionCountMismatch {
        expected: usize,
        actual: usize,
    },
    CapacityOverflow,
    Topology(TopologyBuildError),
    Hierarchy(crate::topology::HierarchyError),
    /// Final model validation failed after topology publication.
    Model(Box<ModelError>),
}

impl fmt::Display for ModelBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InstancePositionCountMismatch { expected, actual } => write!(
                formatter,
                "definition instance requires {expected} positions, but received {actual}"
            ),
            Self::CapacityOverflow => {
                formatter.write_str("model construction exceeds coordinate capacity")
            }
            Self::Topology(error) => write!(formatter, "cannot build topology: {error}"),
            Self::Hierarchy(error) => write!(formatter, "cannot build hierarchy: {error}"),
            Self::Model(error) => write!(formatter, "cannot build model: {error}"),
        }
    }
}

impl std::error::Error for ModelBuildError {}

impl From<ModelError> for ModelBuildError {
    fn from(error: ModelError) -> Self {
        Self::Model(Box::new(error))
    }
}

impl From<TopologyBuildError> for ModelBuildError {
    fn from(error: TopologyBuildError) -> Self {
        Self::Topology(error)
    }
}

impl From<crate::topology::HierarchyError> for ModelBuildError {
    fn from(error: crate::topology::HierarchyError) -> Self {
        Self::Hierarchy(error)
    }
}
