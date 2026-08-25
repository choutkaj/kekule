use std::fmt;
use std::sync::Arc;

use crate::core::{Atom, Bond, Molecule};
use crate::geometry::{PeriodicCell, Point3};
use crate::topology::transform::TopologySubsetError;
use crate::topology::{
    AtomSelection, AtomSiteId, AtomSiteView, ChainId, ChainView, Hierarchy, InstanceAtomId,
    InstanceBondId, MoleculeDefinitionId, MoleculeInstance, MoleculeInstanceId, ResidueId,
    ResidueView, Topology, TopologyAtomIndex, TopologyBuildError, TopologyBuilder, TopologyError,
};
use crate::units::Quantity;

use super::{AtomData, AtomDataError, BondData, BondDataError, PositionError, Positions};

/// One concrete realization of one immutable topology.
#[derive(Debug, Clone)]
pub struct Model {
    pub(super) topology: Arc<Topology>,
    pub(super) positions: Positions,
    pub(super) cell: Option<PeriodicCell>,
    pub(super) atom_data: AtomData,
    pub(super) bond_data: BondData,
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.positions == other.positions
            && self.cell == other.cell
            && self.atom_data == other.atom_data
            && self.bond_data == other.bond_data
    }
}

impl Model {
    /// Creates a non-periodic model with empty atom and bond data.
    pub fn new(topology: Arc<Topology>, positions: Positions) -> Result<Self, ModelError> {
        validate_positions_len(&topology, &positions)?;
        let atom_data = AtomData::new(topology.atom_count());
        let bond_data = BondData::new(topology.bond_count());
        Ok(Self {
            topology,
            positions,
            cell: None,
            atom_data,
            bond_data,
        })
    }

    /// Creates a model from complete dimensioned state.
    pub fn with_atom_data(
        topology: Arc<Topology>,
        positions: Positions,
        cell: Option<PeriodicCell>,
        atom_data: AtomData,
    ) -> Result<Self, ModelError> {
        let bond_data = BondData::new(topology.bond_count());
        Self::with_data(topology, positions, cell, atom_data, bond_data)
    }

    /// Creates a model from complete dimensioned atom and bond state.
    pub fn with_data(
        topology: Arc<Topology>,
        positions: Positions,
        cell: Option<PeriodicCell>,
        atom_data: AtomData,
        bond_data: BondData,
    ) -> Result<Self, ModelError> {
        validate_dimensions(&topology, &positions, &atom_data, &bond_data)?;
        Ok(Self {
            topology,
            positions,
            cell,
            atom_data,
            bond_data,
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
        let atom_data = self
            .atom_data
            .select_indices(&atom_indices)
            .map_err(ModelError::from)?;
        let bond_data = self
            .bond_data
            .select_indices(&bond_indices)
            .map_err(ModelError::from)?;
        Ok(Self::with_data(
            subset.shared_topology(),
            positions,
            self.cell,
            atom_data,
            bond_data,
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

    pub fn atom_data(&self) -> &AtomData {
        &self.atom_data
    }

    pub fn atom_data_mut(&mut self) -> &mut AtomData {
        &mut self.atom_data
    }

    pub fn set_atom_data(&mut self, atom_data: AtomData) -> Result<(), ModelError> {
        if atom_data.len() != self.topology.atom_count() {
            return Err(ModelError::AtomDataCountMismatch {
                expected: self.topology.atom_count(),
                actual: atom_data.len(),
            });
        }
        self.atom_data = atom_data;
        Ok(())
    }

    pub fn bond_data(&self) -> &BondData {
        &self.bond_data
    }

    pub fn bond_data_mut(&mut self) -> &mut BondData {
        &mut self.bond_data
    }

    pub fn set_bond_data(&mut self, bond_data: BondData) -> Result<(), ModelError> {
        if bond_data.len() != self.topology.bond_count() {
            return Err(ModelError::BondDataCountMismatch {
                expected: self.topology.bond_count(),
                actual: bond_data.len(),
            });
        }
        self.bond_data = bond_data;
        Ok(())
    }

    pub fn occupancy(&self, atom: InstanceAtomId) -> Result<Option<f64>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.atom_data.occupancy_at(index.index())?)
    }

    pub fn b_factor(&self, atom: InstanceAtomId) -> Result<Option<Quantity<f64>>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.atom_data.b_factor_at(index.index())?)
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
        Ok(self.atom_data.set_occupancy_at(index.index(), value)?)
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
        Ok(self.atom_data.set_b_factor_at(index.index(), value)?)
    }

    pub fn atom_property_value(
        &self,
        name: &str,
        atom: InstanceAtomId,
    ) -> Result<Option<Quantity<f64>>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.atom_data.property_value_at(name, index.index())?)
    }

    pub fn set_atom_property_value(
        &mut self,
        name: &str,
        atom: InstanceAtomId,
        value: Option<Quantity<f64>>,
    ) -> Result<(), ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self
            .atom_data
            .set_property_value_at(name, index.index(), value)?)
    }

    pub fn bond_property_value(
        &self,
        name: &str,
        bond: InstanceBondId,
    ) -> Result<Option<Quantity<f64>>, ModelError> {
        let index = self
            .topology
            .bond_index(bond)
            .ok_or(ModelError::InvalidBondId(bond))?;
        Ok(self.bond_data.property_value_at(name, index.index())?)
    }

    pub fn set_bond_property_value(
        &mut self,
        name: &str,
        bond: InstanceBondId,
        value: Option<Quantity<f64>>,
    ) -> Result<(), ModelError> {
        let index = self
            .topology
            .bond_index(bond)
            .ok_or(ModelError::InvalidBondId(bond))?;
        Ok(self
            .bond_data
            .set_property_value_at(name, index.index(), value)?)
    }

    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    pub fn view(&self) -> ModelView<'_> {
        ModelView {
            topology: &self.topology,
            positions: &self.positions,
            cell: self.cell.as_ref(),
            atom_data: &self.atom_data,
            bond_data: &self.bond_data,
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
    atom_data: &AtomData,
    bond_data: &BondData,
) -> Result<(), ModelError> {
    validate_positions_len(topology, positions)?;
    if atom_data.len() != topology.atom_count() {
        return Err(ModelError::AtomDataCountMismatch {
            expected: topology.atom_count(),
            actual: atom_data.len(),
        });
    }
    if bond_data.len() != topology.bond_count() {
        return Err(ModelError::BondDataCountMismatch {
            expected: topology.bond_count(),
            actual: bond_data.len(),
        });
    }
    Ok(())
}

/// Borrowed topology, positions, cell, atom data, and bond data for structural
/// kernels.
#[derive(Debug, Clone, Copy)]
pub struct ModelView<'a> {
    topology: &'a Arc<Topology>,
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
    atom_data: &'a AtomData,
    bond_data: &'a BondData,
}

impl<'a> ModelView<'a> {
    pub fn new(
        topology: &'a Arc<Topology>,
        positions: &'a Positions,
        cell: Option<&'a PeriodicCell>,
        atom_data: &'a AtomData,
        bond_data: &'a BondData,
    ) -> Result<Self, ModelError> {
        validate_dimensions(topology, positions, atom_data, bond_data)?;
        Ok(Self {
            topology,
            positions,
            cell,
            atom_data,
            bond_data,
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

    pub const fn atom_data(self) -> &'a AtomData {
        self.atom_data
    }

    pub const fn bond_data(self) -> &'a BondData {
        self.bond_data
    }

    pub fn occupancy(self, atom: InstanceAtomId) -> Result<Option<f64>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.atom_data.occupancy_at(index.index())?)
    }

    pub fn b_factor(self, atom: InstanceAtomId) -> Result<Option<Quantity<f64>>, ModelError> {
        let index = self
            .topology
            .atom_index(atom)
            .ok_or(ModelError::InvalidAtomId(atom))?;
        Ok(self.atom_data.b_factor_at(index.index())?)
    }

    pub fn atom_count(self) -> usize {
        self.topology.atom_count()
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelError {
    PositionCountMismatch { expected: usize, actual: usize },
    AtomDataCountMismatch { expected: usize, actual: usize },
    BondDataCountMismatch { expected: usize, actual: usize },
    InvalidAtomId(InstanceAtomId),
    InvalidBondId(InstanceBondId),
    Position(PositionError),
    AtomData(AtomDataError),
    BondData(BondDataError),
}

impl fmt::Display for ModelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PositionCountMismatch { expected, actual } => write!(
                formatter,
                "topology requires {expected} positions, but received {actual}"
            ),
            Self::AtomDataCountMismatch { expected, actual } => write!(
                formatter,
                "topology requires atom data of length {expected}, but received {actual}"
            ),
            Self::BondDataCountMismatch { expected, actual } => write!(
                formatter,
                "topology requires bond data of length {expected}, but received {actual}"
            ),
            Self::InvalidAtomId(atom) => write!(formatter, "invalid topology atom: {atom}"),
            Self::InvalidBondId(bond) => write!(formatter, "invalid topology bond: {bond}"),
            Self::Position(error) => write!(formatter, "invalid position data: {error}"),
            Self::AtomData(error) => write!(formatter, "invalid atom data: {error}"),
            Self::BondData(error) => write!(formatter, "invalid bond data: {error}"),
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

impl From<AtomDataError> for ModelError {
    fn from(error: AtomDataError) -> Self {
        Self::AtomData(error)
    }
}

impl From<BondDataError> for ModelError {
    fn from(error: BondDataError) -> Self {
        Self::BondData(error)
    }
}

/// Convenience builder that assembles topology and one complete model.
///
/// Geometry is supplied explicitly and remains outside the molecule.
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
    pub fn topology_builder_mut(&mut self) -> &mut TopologyBuilder {
        &mut self.topology
    }

    pub fn add_molecule_definition(
        &mut self,
        molecule: &Molecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self.topology.add_molecule_definition(molecule)?)
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

    pub fn build(self) -> Result<Model, ModelBuildError> {
        let topology = Arc::new(self.topology.build()?);
        let positions = Positions::from_canonical_values(self.positions);
        Ok(Model::new(topology, positions).expect("builder creates dimensionally valid state"))
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
    InstancePositionCountMismatch { expected: usize, actual: usize },
    CapacityOverflow,
    Topology(TopologyBuildError),
    Hierarchy(crate::topology::HierarchyError),
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
        }
    }
}

impl std::error::Error for ModelBuildError {}

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
