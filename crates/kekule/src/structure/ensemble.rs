use std::fmt;
use std::sync::Arc;

use crate::core::{Molecule, PropMap};
use crate::geometry::PeriodicCell;
use crate::topology::{Topology, TopologyBuildError, TopologyBuilder};

use super::{AtomData, BondData, Model, ModelView, PositionError, Positions};

/// One finite non-temporal ensemble member.
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleMember {
    positions: Positions,
    cell: Option<PeriodicCell>,
    atom_data: AtomData,
    bond_data: BondData,
    weight: Option<f64>,
    props: PropMap,
}

impl EnsembleMember {
    pub fn new(positions: Positions, bond_count: usize) -> Self {
        let atom_data = AtomData::new(positions.len());
        let bond_data = BondData::new(bond_count);
        Self {
            positions,
            cell: None,
            atom_data,
            bond_data,
            weight: None,
            props: PropMap::new(),
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

    pub fn atom_data(&self) -> &AtomData {
        &self.atom_data
    }

    pub fn atom_data_mut(&mut self) -> &mut AtomData {
        &mut self.atom_data
    }

    pub fn set_atom_data(&mut self, atom_data: AtomData) -> Result<(), EnsembleError> {
        if atom_data.len() != self.positions.len() {
            return Err(EnsembleError::AtomDataCountMismatch {
                expected: self.positions.len(),
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

    pub fn set_bond_data(&mut self, bond_data: BondData) -> Result<(), EnsembleError> {
        if bond_data.len() != self.bond_data.len() {
            return Err(EnsembleError::BondDataCountMismatch {
                expected: self.bond_data.len(),
                actual: bond_data.len(),
            });
        }
        self.bond_data = bond_data;
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

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    pub fn view<'a>(&'a self, topology: &'a Arc<Topology>) -> Result<ModelView<'a>, EnsembleError> {
        ModelView::new(
            topology,
            &self.positions,
            self.cell.as_ref(),
            &self.atom_data,
            &self.bond_data,
        )
        .map_err(|error| EnsembleError::Model(Box::new(error)))
    }
}

/// A finite stable-order collection of non-temporal models.
#[derive(Debug, Clone)]
pub struct Ensemble {
    topology: Arc<Topology>,
    members: Vec<EnsembleMember>,
}

impl Ensemble {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            topology,
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
                atom_data: model.atom_data.clone(),
                bond_data: model.bond_data.clone(),
                weight: None,
                props: PropMap::new(),
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
        if member.atom_data.len() != self.topology.atom_count() {
            return Err(EnsembleError::AtomDataCountMismatch {
                expected: self.topology.atom_count(),
                actual: member.atom_data.len(),
            });
        }
        if member.bond_data.len() != self.topology.bond_count() {
            return Err(EnsembleError::BondDataCountMismatch {
                expected: self.topology.bond_count(),
                actual: member.bond_data.len(),
            });
        }
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
    AtomDataCountMismatch { expected: usize, actual: usize },
    BondDataCountMismatch { expected: usize, actual: usize },
    InvalidWeight,
    MissingWeight { member: usize },
    ZeroTotalWeight,
    TopologyBuild(TopologyBuildError),
    Position(PositionError),
    Model(Box<super::ModelError>),
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
            Self::AtomDataCountMismatch { expected, actual } => write!(
                formatter,
                "ensemble topology requires atom data of length {expected}, but received {actual}"
            ),
            Self::BondDataCountMismatch { expected, actual } => write!(
                formatter,
                "ensemble topology requires bond data of length {expected}, but received {actual}"
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
