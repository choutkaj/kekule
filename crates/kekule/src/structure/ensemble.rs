use std::fmt;
use std::sync::Arc;

use crate::core::Molecule;
use crate::geometry::PeriodicCell;
use crate::properties::{
    Properties, PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue,
};
use crate::topology::transform::TopologySubsetError;
use crate::topology::{AtomSelection, Topology, TopologyBuildError, TopologyBuilder};
use crate::units::Quantity;

use super::{Model, ModelView, PositionError, Positions};

/// One finite non-temporal ensemble member.
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleMember {
    positions: Positions,
    cell: Option<PeriodicCell>,
    properties: Properties,
    weight: Option<f64>,
}

impl EnsembleMember {
    pub fn new(positions: Positions, bond_count: usize) -> Self {
        let properties = Properties::realization(positions.len(), bond_count);
        Self {
            positions,
            cell: None,
            properties,
            weight: None,
        }
    }

    pub fn positions(&self) -> &Positions {
        &self.positions
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
    ) -> Result<Option<PropertyValue>, EnsembleError> {
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

    pub fn atom_property_value(
        &self,
        key: &PropertyKey,
        index: usize,
    ) -> Result<Option<PropertyValue>, EnsembleError> {
        Ok(self.atom_properties().value(key, index)?)
    }

    pub fn set_atom_property_value(
        &mut self,
        key: PropertyKey,
        index: usize,
        value: Option<PropertyValue>,
    ) -> Result<(), EnsembleError> {
        Ok(self
            .properties
            .set_realization_atom_value(key, index, value)?)
    }

    pub fn insert_atom_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, EnsembleError> {
        Ok(self
            .properties
            .insert_realization_atom_column(key, column)?)
    }

    pub fn remove_atom_property_column(
        &mut self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, EnsembleError> {
        Ok(self.properties.remove_realization_atom_column(key)?)
    }

    pub fn bond_property_value(
        &self,
        key: &PropertyKey,
        index: usize,
    ) -> Result<Option<PropertyValue>, EnsembleError> {
        Ok(self.bond_properties().value(key, index)?)
    }

    pub fn set_bond_property_value(
        &mut self,
        key: PropertyKey,
        index: usize,
        value: Option<PropertyValue>,
    ) -> Result<(), EnsembleError> {
        Ok(self
            .properties
            .set_realization_bond_value(key, index, value)?)
    }

    pub fn insert_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, EnsembleError> {
        Ok(self
            .properties
            .insert_realization_bond_column(key, column)?)
    }

    pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.properties.remove_realization_bond_column(key)
    }

    pub fn occupancy_at(&self, index: usize) -> Result<Option<f64>, EnsembleError> {
        Ok(self.properties.occupancy_at(index)?)
    }

    pub fn set_occupancy_at(
        &mut self,
        index: usize,
        value: Option<f64>,
    ) -> Result<(), EnsembleError> {
        Ok(self.properties.set_occupancy_at(index, value)?)
    }

    pub fn b_factor_at(&self, index: usize) -> Result<Option<Quantity<f64>>, EnsembleError> {
        Ok(self.properties.b_factor_at(index)?)
    }

    pub fn set_b_factor_at(
        &mut self,
        index: usize,
        value: Option<Quantity<f64>>,
    ) -> Result<(), EnsembleError> {
        Ok(self.properties.set_b_factor_at(index, value)?)
    }

    pub fn set_properties(&mut self, properties: Properties) -> Result<(), EnsembleError> {
        if properties.realization_atom_properties().len() != self.positions.len() {
            return Err(EnsembleError::AtomPropertyCountMismatch {
                expected: self.positions.len(),
                actual: properties.realization_atom_properties().len(),
            });
        }
        if properties.realization_bond_properties().len() != self.bond_properties().len() {
            return Err(EnsembleError::BondPropertyCountMismatch {
                expected: self.bond_properties().len(),
                actual: properties.realization_bond_properties().len(),
            });
        }
        properties.validate_realization_canonical_properties()?;
        self.properties = properties;
        Ok(())
    }

    pub const fn weight(&self) -> Option<f64> {
        self.weight
    }

    pub fn set_weight(&mut self, weight: Option<f64>) -> Result<(), EnsembleError> {
        if weight.is_some_and(|weight| !weight.is_finite() || weight < 0.0) {
            return Err(EnsembleError::InvalidWeight);
        }
        self.weight = weight;
        Ok(())
    }

    pub fn view<'a>(&'a self, topology: &'a Arc<Topology>) -> Result<ModelView<'a>, EnsembleError> {
        ModelView::new(
            topology,
            &self.positions,
            self.cell.as_ref(),
            &self.properties,
        )
        .map_err(|error| EnsembleError::Model(Box::new(error)))
    }
}

/// A finite stable-order collection of non-temporal models.
#[derive(Debug, Clone)]
pub struct Ensemble {
    topology: Arc<Topology>,
    properties: Properties,
    members: Vec<EnsembleMember>,
}

impl Ensemble {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            topology,
            properties: Properties::new(),
            members: Vec::new(),
        }
    }

    pub fn from_members(
        topology: Arc<Topology>,
        members: impl IntoIterator<Item = EnsembleMember>,
    ) -> Result<Self, EnsembleError> {
        let mut ensemble = Self::new(topology);
        for member in members {
            ensemble.push(member)?;
        }
        Ok(ensemble)
    }

    pub fn from_models(models: &[Model]) -> Result<Self, EnsembleError> {
        let first = models.first().ok_or(EnsembleError::EmptySource)?;
        let mut ensemble = Self::new(Arc::clone(&first.topology));
        for model in models {
            if !Arc::ptr_eq(&first.topology, &model.topology) {
                return Err(EnsembleError::TopologyMismatch);
            }
            let member = EnsembleMember {
                positions: model.positions.clone(),
                cell: model.cell,
                properties: model.properties.clone(),
                weight: None,
            };
            ensemble.members.push(member);
        }
        Ok(ensemble)
    }

    /// Builds an ensemble from dense member positions in molecule atom order.
    pub fn from_molecule_positions(
        molecule: &Molecule,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<Self, EnsembleError> {
        let mut topology_builder = TopologyBuilder::new();
        let definition = topology_builder.add_molecule_definition(molecule)?;
        topology_builder.add_instance(definition)?;
        let topology = Arc::new(topology_builder.build()?);
        let mut ensemble = Self::new(Arc::clone(&topology));
        for positions in positions {
            ensemble.push(EnsembleMember::new(positions, topology.bond_count()))?;
        }
        Ok(ensemble)
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub const fn properties(&self) -> &Properties {
        &self.properties
    }

    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, EnsembleError> {
        Ok(self.properties.insert(key, value)?)
    }

    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.properties.remove(key)
    }

    pub fn clear_properties(&mut self) {
        self.properties.clear_owner();
    }

    /// Constructs one topology subset and applies its dense mapping to every member.
    pub fn slice(&self, selection: &AtomSelection) -> Result<Self, EnsembleSliceError> {
        let subset = self.topology.subset(selection)?;
        let atom_indices = subset
            .correspondence()
            .source_atom_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let bond_indices = subset
            .correspondence()
            .source_bond_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let mut target = Self::new(subset.shared_topology());
        for member in &self.members {
            target.push(EnsembleMember {
                positions: member.positions.select_indices(&atom_indices)?,
                cell: member.cell,
                properties: member
                    .properties
                    .project_realization(&atom_indices, &bond_indices)?,
                weight: member.weight,
            })?;
        }
        Ok(target)
    }

    pub fn members(&self) -> impl ExactSizeIterator<Item = &EnsembleMember> {
        self.members.iter()
    }

    pub fn member(&self, index: usize) -> Option<&EnsembleMember> {
        self.members.get(index)
    }

    pub fn member_mut(&mut self, index: usize) -> Option<&mut EnsembleMember> {
        self.members.get_mut(index)
    }

    pub fn len(&self) -> usize {
        self.members.len()
    }

    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    pub fn push(&mut self, member: EnsembleMember) -> Result<(), EnsembleError> {
        if member.positions.len() != self.topology.atom_count() {
            return Err(EnsembleError::PositionCountMismatch {
                expected: self.topology.atom_count(),
                actual: member.positions.len(),
            });
        }
        if member.properties.atoms().len() != self.topology.atom_count() {
            return Err(EnsembleError::AtomPropertyCountMismatch {
                expected: self.topology.atom_count(),
                actual: member.properties.atoms().len(),
            });
        }
        if member.properties.bonds().len() != self.topology.bond_count() {
            return Err(EnsembleError::BondPropertyCountMismatch {
                expected: self.topology.bond_count(),
                actual: member.properties.bonds().len(),
            });
        }
        member
            .properties
            .validate_realization_canonical_properties()?;
        self.members.push(member);
        Ok(())
    }

    pub fn views(&self) -> impl ExactSizeIterator<Item = ModelView<'_>> {
        self.members.iter().map(|member| {
            member
                .view(&self.topology)
                .expect("ensemble validates shared topology allocation")
        })
    }

    pub fn normalize_weights(&mut self) -> Result<(), EnsembleError> {
        if self.members.is_empty() {
            return Err(EnsembleError::EmptySource);
        }
        let mut total = 0.0;
        for (index, member) in self.members.iter().enumerate() {
            let weight = member
                .weight
                .ok_or(EnsembleError::MissingWeight { member: index })?;
            total += weight;
        }
        if !total.is_finite() || total <= 0.0 {
            return Err(EnsembleError::ZeroTotalWeight);
        }
        for member in &mut self.members {
            member.weight = member.weight.map(|weight| weight / total);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum EnsembleError {
    EmptySource,
    TopologyMismatch,
    PositionCountMismatch { expected: usize, actual: usize },
    AtomPropertyCountMismatch { expected: usize, actual: usize },
    BondPropertyCountMismatch { expected: usize, actual: usize },
    InvalidWeight,
    MissingWeight { member: usize },
    ZeroTotalWeight,
    TopologyBuild(TopologyBuildError),
    Position(PositionError),
    Model(Box<super::ModelError>),
    Property(Box<PropertyError>),
}

impl fmt::Display for EnsembleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySource => formatter.write_str("ensemble source is empty"),
            Self::TopologyMismatch => {
                formatter.write_str("ensemble member belongs to a different topology")
            }
            Self::PositionCountMismatch { expected, actual } => write!(
                formatter,
                "ensemble topology requires {expected} positions, but received {actual}"
            ),
            Self::AtomPropertyCountMismatch { expected, actual } => write!(
                formatter,
                "ensemble topology requires atom properties of length {expected}, but received {actual}"
            ),
            Self::BondPropertyCountMismatch { expected, actual } => write!(
                formatter,
                "ensemble topology requires bond properties of length {expected}, but received {actual}"
            ),
            Self::InvalidWeight => {
                formatter.write_str("ensemble weight must be finite and non-negative")
            }
            Self::MissingWeight { member } => {
                write!(formatter, "ensemble member {member} has no weight")
            }
            Self::ZeroTotalWeight => {
                formatter.write_str("ensemble weights must have a positive finite total")
            }
            Self::TopologyBuild(error) => write!(formatter, "cannot build topology: {error}"),
            Self::Position(error) => write!(formatter, "cannot build member positions: {error}"),
            Self::Model(error) => write!(formatter, "invalid ensemble member state: {error}"),
            Self::Property(error) => write!(formatter, "invalid ensemble member property: {error}"),
        }
    }
}

impl std::error::Error for EnsembleError {}

impl From<TopologyBuildError> for EnsembleError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}

impl From<PositionError> for EnsembleError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<PropertyError> for EnsembleError {
    fn from(error: PropertyError) -> Self {
        Self::Property(Box::new(error))
    }
}

/// Failure to subset an ensemble topology or transfer member state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum EnsembleSliceError {
    Topology(TopologySubsetError),
    Position(PositionError),
    Property(PropertyError),
    Ensemble(EnsembleError),
}

impl fmt::Display for EnsembleSliceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "cannot slice ensemble: {self:?}")
    }
}

impl std::error::Error for EnsembleSliceError {}

impl From<TopologySubsetError> for EnsembleSliceError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Topology(error)
    }
}
impl From<PositionError> for EnsembleSliceError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}
impl From<PropertyError> for EnsembleSliceError {
    fn from(error: PropertyError) -> Self {
        Self::Property(error)
    }
}
impl From<EnsembleError> for EnsembleSliceError {
    fn from(error: EnsembleError) -> Self {
        Self::Ensemble(error)
    }
}
