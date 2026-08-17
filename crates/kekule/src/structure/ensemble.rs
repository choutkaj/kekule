use std::fmt;
use std::sync::Arc;

use crate::core::{ConformerId, PropMap};
use crate::geometry::PeriodicCell;
use crate::small::SmallMolecule;
use crate::topology::{
    MoleculeInstanceMetadata, Topology, TopologyBuildError, TopologyBuilder, TopologyMapping,
};
use crate::units::{Quantity, MODEL_LENGTH_UNIT};

use super::model::stage_conformer_positions;
use super::{
    AtomData, BondData, Model, ModelBuildError, ModelView, PositionError, Positions,
    TopologyRemapError,
};

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
    pub fn new(positions: Positions) -> Self {
        let atom_data = AtomData::new(positions.topology_arc());
        let bond_data = BondData::new(positions.topology_arc());
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
        if !Arc::ptr_eq(atom_data.topology_arc(), self.positions.topology_arc()) {
            return Err(EnsembleError::TopologyMismatch);
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
        if !Arc::ptr_eq(bond_data.topology_arc(), self.positions.topology_arc()) {
            return Err(EnsembleError::TopologyMismatch);
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
        .map_err(|_| EnsembleError::TopologyMismatch)
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

    pub fn from_small_molecule_conformers(
        molecule: &SmallMolecule,
        conformers: impl IntoIterator<Item = ConformerId>,
    ) -> Result<Self, EnsembleError> {
        let mut topology_builder = TopologyBuilder::new();
        let definition = topology_builder.add_small_molecule_definition(molecule)?;
        topology_builder.add_instance(definition, MoleculeInstanceMetadata::default())?;
        let topology = Arc::new(topology_builder.build()?);
        let mut ensemble = Self::new(Arc::clone(&topology));
        for conformer in conformers {
            let positions = stage_conformer_positions(molecule.graph(), conformer)
                .map_err(|error| EnsembleError::ModelBuild(Box::new(error)))?;
            let positions = Positions::new(&topology, Quantity::new(positions, MODEL_LENGTH_UNIT))?;
            ensemble.push(EnsembleMember::new(positions))?;
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
        if !member.positions.is_compatible(&self.topology)
            || !member.atom_data.is_compatible(&self.topology)
            || !member.bond_data.is_compatible(&self.topology)
        {
            return Err(EnsembleError::TopologyMismatch);
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

    /// Remaps every member to one exact target topology without renormalizing
    /// weights or changing member order.
    pub fn remap_to(
        &self,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TopologyRemapError> {
        if !mapping.is_source(&self.topology) {
            return Err(TopologyRemapError::MappingSourceMismatch);
        }
        if !mapping.is_target(target) {
            return Err(TopologyRemapError::MappingTargetMismatch);
        }
        let mut members = Vec::with_capacity(self.members.len());
        for (member_index, member) in self.members.iter().enumerate() {
            let remap_member = || -> Result<EnsembleMember, TopologyRemapError> {
                let positions = member.positions.remap_to(&self.topology, target, mapping)?;
                let atom_data = member.atom_data.remap_to(&self.topology, target, mapping)?;
                let bond_data = member.bond_data.remap_to(&self.topology, target, mapping)?;
                Ok(EnsembleMember {
                    positions,
                    cell: member.cell,
                    atom_data,
                    bond_data,
                    weight: member.weight,
                    props: member.props.clone(),
                })
            };
            members.push(remap_member().map_err(|error| TopologyRemapError::Member {
                member: member_index,
                error: Box::new(error),
            })?);
        }
        Ok(Self {
            topology: Arc::clone(target),
            members,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum EnsembleError {
    EmptySource,
    TopologyMismatch,
    InvalidWeight,
    MissingWeight { member: usize },
    ZeroTotalWeight,
    TopologyBuild(TopologyBuildError),
    ModelBuild(Box<ModelBuildError>),
    Position(PositionError),
}

impl fmt::Display for EnsembleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySource => formatter.write_str("ensemble source is empty"),
            Self::TopologyMismatch => {
                formatter.write_str("ensemble member belongs to a different topology")
            }
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
            Self::ModelBuild(error) => write!(formatter, "cannot build ensemble member: {error}"),
            Self::Position(error) => write!(formatter, "cannot build member positions: {error}"),
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
