use std::fmt;
use std::sync::Arc;

use crate::core::Molecule;
use crate::geometry::PeriodicCell;
use crate::properties::{
    Properties, PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue,
};
use crate::topology::transform::TopologySubsetError;
use crate::topology::{
    AtomSelection, Topology, TopologyBuildError, TopologyBuilder, TopologyPerceptionError,
};
use crate::units::Quantity;

use super::{Model, ModelView, PositionError, Positions};

/// One finite non-temporal ensemble member.
///
/// The member stores realization state but no topology. It becomes
/// topology-bound only when inserted into an [`Ensemble`] or viewed as an
/// [`EnsembleMemberView`].
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleMember {
    positions: Positions,
    cell: Option<PeriodicCell>,
    properties: Properties,
    weight: Option<f64>,
}

impl EnsembleMember {
    /// Constructs a detached member. Empty property domains are sized on insertion.
    pub fn new(positions: Positions) -> Self {
        let properties = Properties::realization(positions.len(), 0);
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

    pub fn atom_property(
        &self,
        index: usize,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, EnsembleError> {
        Ok(self.atom_properties().value(key, index)?)
    }

    pub fn set_atom_property(
        &mut self,
        index: usize,
        key: PropertyKey,
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

    pub fn bond_property(
        &self,
        index: usize,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, EnsembleError> {
        Ok(self.bond_properties().value(key, index)?)
    }

    pub fn set_bond_property(
        &mut self,
        index: usize,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), EnsembleError> {
        Ok(self
            .properties
            .set_realization_bond_value(key, index, value)?)
    }

    /// Supplies a complete bond column; the first retained column establishes
    /// the detached domain size. The ensemble checks that size at insertion.
    pub fn insert_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, EnsembleError> {
        if !self.bond_properties().has_data() {
            let mut properties = self.properties.clone();
            properties.normalize_realization_dimensions(self.positions.len(), column.len())?;
            let previous = properties.insert_realization_bond_column(key, column)?;
            self.properties = properties;
            return Ok(previous);
        }
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

    /// Replaces detached properties. Empty atom domains adopt the position count;
    /// populated atom columns are checked now and bond columns at insertion.
    pub fn set_properties(&mut self, mut properties: Properties) -> Result<(), EnsembleError> {
        if properties.realization_atom_properties().has_data()
            && properties.realization_atom_properties().len() != self.positions.len()
        {
            return Err(EnsembleError::AtomPropertyCountMismatch {
                expected: self.positions.len(),
                actual: properties.realization_atom_properties().len(),
            });
        }
        properties.normalize_realization_dimensions(
            self.positions.len(),
            properties.realization_bond_properties().len(),
        )?;
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

    /// Borrows already dimensioned state. Prefer [`Ensemble::member`] to obtain
    /// a view after the collection has established empty property domains.
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

/// Restricted mutable access to a member stored in an [`Ensemble`].
///
/// Read access uses [`std::ops::Deref`]. Only dimension-preserving mutation is
/// exposed: there is deliberately no `DerefMut` or mutable payload accessor.
/// Edits take effect immediately; dropping or forgetting the editor does not
/// perform validation. Use [`Ensemble::replace_member`] to replace the payload.
#[derive(Debug)]
pub struct EnsembleMemberMut<'a> {
    member: &'a mut EnsembleMember,
}

impl std::ops::Deref for EnsembleMemberMut<'_> {
    type Target = EnsembleMember;

    fn deref(&self) -> &Self::Target {
        self.member
    }
}

impl EnsembleMemberMut<'_> {
    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.member.set_cell(cell);
    }

    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, EnsembleError> {
        self.member.insert_property(key, value)
    }

    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.member.remove_property(key)
    }

    pub fn clear_properties(&mut self) {
        self.member.clear_properties();
    }

    pub fn set_atom_property(
        &mut self,
        index: usize,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), EnsembleError> {
        self.member.set_atom_property(index, key, value)
    }

    pub fn insert_atom_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, EnsembleError> {
        self.member.insert_atom_property_column(key, column)
    }

    pub fn remove_atom_property_column(
        &mut self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, EnsembleError> {
        self.member.remove_atom_property_column(key)
    }

    pub fn set_bond_property(
        &mut self,
        index: usize,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), EnsembleError> {
        self.member.set_bond_property(index, key, value)
    }

    pub fn insert_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, EnsembleError> {
        Ok(self
            .member
            .properties
            .insert_realization_bond_column(key, column)?)
    }

    pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.member.remove_bond_property_column(key)
    }

    pub fn set_occupancy_at(
        &mut self,
        index: usize,
        value: Option<f64>,
    ) -> Result<(), EnsembleError> {
        self.member.set_occupancy_at(index, value)
    }

    pub fn set_b_factor_at(
        &mut self,
        index: usize,
        value: Option<Quantity<f64>>,
    ) -> Result<(), EnsembleError> {
        self.member.set_b_factor_at(index, value)
    }

    pub fn set_properties(&mut self, mut properties: Properties) -> Result<(), EnsembleError> {
        let atom_count = self.member.positions.len();
        let bond_count = self.member.bond_properties().len();
        if properties.realization_atom_properties().len() != atom_count {
            return Err(EnsembleError::AtomPropertyCountMismatch {
                expected: atom_count,
                actual: properties.realization_atom_properties().len(),
            });
        }
        if properties.realization_bond_properties().len() != bond_count {
            return Err(EnsembleError::BondPropertyCountMismatch {
                expected: bond_count,
                actual: properties.realization_bond_properties().len(),
            });
        }
        properties.normalize_realization_dimensions(atom_count, bond_count)?;
        self.member.properties = properties;
        Ok(())
    }

    pub fn set_weight(&mut self, weight: Option<f64>) -> Result<(), EnsembleError> {
        self.member.set_weight(weight)
    }
}

/// Borrowed topology-bound view of one ensemble member.
#[derive(Debug, Clone, Copy)]
pub struct EnsembleMemberView<'a> {
    topology: &'a Arc<Topology>,
    member: &'a EnsembleMember,
}

impl<'a> EnsembleMemberView<'a> {
    fn new(topology: &'a Arc<Topology>, member: &'a EnsembleMember) -> Self {
        Self { topology, member }
    }

    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    pub fn positions(self) -> &'a Positions {
        self.member.positions()
    }

    pub fn cell(self) -> Option<&'a PeriodicCell> {
        self.member.cell()
    }

    pub fn properties(self) -> &'a Properties {
        self.member.properties()
    }

    pub fn atom_properties(self) -> &'a PropertyTable {
        self.member.atom_properties()
    }

    pub fn bond_properties(self) -> &'a PropertyTable {
        self.member.bond_properties()
    }

    pub fn atom_property(
        self,
        index: usize,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, EnsembleError> {
        self.member.atom_property(index, key)
    }

    pub fn bond_property(
        self,
        index: usize,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, EnsembleError> {
        self.member.bond_property(index, key)
    }

    pub fn occupancy_at(self, index: usize) -> Result<Option<f64>, EnsembleError> {
        self.member.occupancy_at(index)
    }

    pub fn b_factor_at(self, index: usize) -> Result<Option<Quantity<f64>>, EnsembleError> {
        self.member.b_factor_at(index)
    }

    pub fn weight(self) -> Option<f64> {
        self.member.weight()
    }

    /// Projects this member into zero-copy borrowed model semantics.
    pub fn as_model(self) -> ModelView<'a> {
        self.member
            .view(self.topology)
            .expect("ensemble member view has validated topology")
    }

    /// Materializes this member as an owned model, excluding ensemble weight.
    pub fn to_model(self) -> Model {
        self.as_model().to_model()
    }
}

/// A finite stable-order collection of non-temporal realizations.
///
/// Every member is validated against one shared [`Topology`]. Use an ensemble
/// for conformers, alternate experimental models, or another non-temporal
/// finite set of realizations. Use the `kekule-traj` companion crate when frame
/// order has temporal meaning.
#[derive(Debug, Clone)]
pub struct Ensemble {
    topology: Arc<Topology>,
    properties: Properties,
    members: Vec<EnsembleMember>,
}

impl Ensemble {
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self {
            topology: topology.into(),
            properties: Properties::new(),
            members: Vec::new(),
        }
    }

    pub fn from_members(
        topology: impl Into<Arc<Topology>>,
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

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once for the collection, independent
    /// of member count. All member positions, cells, weights, and properties,
    /// and collection properties are retained without copying members. Failure
    /// leaves the entire ensemble unchanged. Other owners and existing bound
    /// selections or prepared calculations retain their original snapshot.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
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

    /// Iterates topology-bound member views in stable collection order.
    pub fn members(&self) -> impl ExactSizeIterator<Item = EnsembleMemberView<'_>> {
        self.members
            .iter()
            .map(|member| EnsembleMemberView::new(&self.topology, member))
    }

    /// Returns one topology-bound member view by stable collection index.
    pub fn member(&self, index: usize) -> Option<EnsembleMemberView<'_>> {
        self.members
            .get(index)
            .map(|member| EnsembleMemberView::new(&self.topology, member))
    }

    /// Returns an editor that preserves the member's validated dimensions.
    ///
    /// Bind the editor with `let mut member = ensemble.member_mut(index).unwrap()`
    /// when making several edits. Whole-member replacement must go through
    /// [`Self::replace_member`].
    ///
    /// ```compile_fail,E0594
    /// use kekule::structure::{Ensemble, EnsembleMember, Positions};
    /// fn overwrite(ensemble: &mut Ensemble) {
    ///     *ensemble.member_mut(0).unwrap() = EnsembleMember::new(Positions::zeros(1));
    /// }
    /// ```
    pub fn member_mut(&mut self, index: usize) -> Option<EnsembleMemberMut<'_>> {
        self.members
            .get_mut(index)
            .map(|member| EnsembleMemberMut { member })
    }

    /// Replaces a member after validating it against this ensemble's topology.
    ///
    /// Returns the previous member. An invalid index or incompatible replacement
    /// returns an error and leaves the ensemble unchanged. The index is checked
    /// first. Positions and property columns must follow this topology's dense
    /// atom and bond order, just as for [`Self::push`].
    pub fn replace_member(
        &mut self,
        index: usize,
        mut member: EnsembleMember,
    ) -> Result<EnsembleMember, EnsembleError> {
        if index >= self.members.len() {
            return Err(EnsembleError::MemberIndexOutOfBounds {
                index,
                len: self.members.len(),
            });
        }
        self.prepare_member(&mut member)?;
        Ok(std::mem::replace(&mut self.members[index], member))
    }

    pub fn len(&self) -> usize {
        self.members.len()
    }

    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    pub fn push(&mut self, mut member: EnsembleMember) -> Result<(), EnsembleError> {
        self.prepare_member(&mut member)?;
        self.members.push(member);
        Ok(())
    }

    fn prepare_member(&self, member: &mut EnsembleMember) -> Result<(), EnsembleError> {
        if member.positions.len() != self.topology.atom_count() {
            return Err(EnsembleError::PositionCountMismatch {
                expected: self.topology.atom_count(),
                actual: member.positions.len(),
            });
        }
        if member.properties.atoms().has_data()
            && member.properties.atoms().len() != self.topology.atom_count()
        {
            return Err(EnsembleError::AtomPropertyCountMismatch {
                expected: self.topology.atom_count(),
                actual: member.properties.atoms().len(),
            });
        }
        if member.properties.bonds().has_data()
            && member.properties.bonds().len() != self.topology.bond_count()
        {
            return Err(EnsembleError::BondPropertyCountMismatch {
                expected: self.topology.bond_count(),
                actual: member.properties.bonds().len(),
            });
        }
        member.properties.normalize_realization_dimensions(
            self.topology.atom_count(),
            self.topology.bond_count(),
        )?;
        Ok(())
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
    MemberIndexOutOfBounds { index: usize, len: usize },
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
            Self::MemberIndexOutOfBounds { index, len } => {
                write!(formatter, "ensemble member index {index} is out of bounds for {len} members")
            }
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
