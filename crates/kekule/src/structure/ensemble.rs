use std::sync::Arc;

use crate::core::Molecule;
use crate::geometry::{PeriodicCell, RigidTransform};
use crate::properties::{OwnerProperties, RealizationProperties, TopologyPropertiesMut};
use crate::topology::{AtomSelection, Topology, TopologyBuilder, TopologyPerceptionError};

use super::scalars::AtomScalars;
use super::state::{
    is_bound, realization_mutators, realization_readers, site_mutators, site_readers, subset_rows,
    transform_cell, transform_positions, validate_binding,
};
use super::{AsModelView, ConformationError, Model, ModelView, Positions, RealizationError};

/// A weighted, unordered sample of realizations of one shared topology, such
/// as conformers, Monte Carlo samples, or alternate experimental models.
///
/// Every [`EnsembleMember`] carries a finite positive statistical weight.
/// Weights are relative: only their ratios are meaningful, so selecting,
/// duplicating, or removing members needs no renormalization. Member order is
/// stable but carries no physical meaning; use [`super::Trajectory`] when order
/// is temporal. Statistics over an ensemble must use its weights.
///
/// Members are validated against the shared topology on insertion, so they
/// can be read as [`super::ModelView`]s without further checks. Superposition
/// and RMSD are inherent methods; see [`crate::alignment`].
///
/// ```
/// use kekule::structure::{Ensemble, EnsembleMember, Positions};
/// let topology = kekule::smiles::to_topology("O")?;
/// let atoms = topology.atom_count();
/// let mut ensemble = Ensemble::new(topology);
/// ensemble.push(EnsembleMember::new(Positions::zeros(atoms), 3.0)?)?;
/// ensemble.push(EnsembleMember::new(Positions::zeros(atoms), 1.0)?)?;
/// ensemble.normalize_weights()?;
/// assert_eq!(ensemble.get(0).unwrap().weight(), 0.75);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// Member editors preserve dimensions and cannot overwrite a whole member:
///
/// ```compile_fail,E0594
/// use kekule::structure::{Ensemble, EnsembleMember, Positions};
/// fn overwrite(ensemble: &mut Ensemble) {
///     *ensemble.get_mut(0).unwrap() = EnsembleMember::new(Positions::zeros(1), 1.0).unwrap();
/// }
/// ```
#[derive(Debug, Clone)]
pub struct Ensemble {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    members: Vec<EnsembleMember>,
}

impl Ensemble {
    /// An empty ensemble of one shared topology.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            members: Vec::new(),
        }
    }

    /// Validates and collects each member in order, binding its bond property
    /// rows. Fails on the first invalid member.
    pub fn from_items(
        topology: impl Into<Arc<Topology>>,
        members: impl IntoIterator<Item = EnsembleMember>,
    ) -> Result<Self, RealizationError> {
        let mut ensemble = Self::new(topology);
        for member in members {
            ensemble.push(member)?;
        }
        Ok(ensemble)
    }

    /// Collects models that share one topology layout as equally weighted
    /// members (weight `1.0`), keeping the first model's snapshot. Positions,
    /// cells, occupancies, B-factors, and properties move without copying.
    pub fn from_models(models: impl IntoIterator<Item = Model>) -> Result<Self, RealizationError> {
        let mut models = models.into_iter();
        let first = models.next().ok_or(RealizationError::EmptySource)?;
        let mut ensemble = Self::new(first.shared_topology());
        ensemble.members.push(EnsembleMember::from_model(first));
        for model in models {
            if !model.topology().shares_layout(&ensemble.topology) {
                return Err(RealizationError::TopologyMismatch);
            }
            ensemble.members.push(EnsembleMember::from_model(model));
        }
        Ok(ensemble)
    }

    /// Builds equally weighted conformers of one molecule from dense
    /// positions in molecule atom order.
    pub fn from_molecule_positions(
        molecule: Molecule,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<Self, RealizationError> {
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(molecule)?;
        builder.add_instance(definition)?;
        Self::from_items(
            builder.build()?,
            positions.into_iter().map(EnsembleMember::unit),
        )
    }

    /// Members already bound to `topology`, for projections from another
    /// collection of the same topology.
    pub(crate) fn from_bound_members(
        topology: Arc<Topology>,
        properties: OwnerProperties,
        members: Vec<EnsembleMember>,
    ) -> Self {
        debug_assert!(members.iter().all(|member| is_bound(
            member.atom_count(),
            &member.properties,
            &topology
        )));
        Self {
            topology,
            properties,
            members,
        }
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Ensemble-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }

    pub fn len(&self) -> usize {
        self.members.len()
    }

    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    /// One topology-bound member by stable index.
    pub fn get(&self, index: usize) -> Option<EnsembleMemberView<'_>> {
        self.members.get(index).map(|member| EnsembleMemberView {
            topology: &self.topology,
            member,
        })
    }

    /// A dimension-preserving editor for one member. Whole-member replacement
    /// goes through [`Self::replace`].
    pub fn get_mut(&mut self, index: usize) -> Option<EnsembleMemberMut<'_>> {
        self.members
            .get_mut(index)
            .map(|member| EnsembleMemberMut { member })
    }

    /// Topology-bound members in stable order.
    pub fn iter(&self) -> EnsembleIter<'_> {
        EnsembleIter {
            topology: &self.topology,
            members: self.members.iter(),
        }
    }

    /// Validates and appends one member, binding its bond property rows.
    pub fn push(&mut self, mut member: EnsembleMember) -> Result<(), RealizationError> {
        self.bind(&mut member)?;
        self.members.push(member);
        Ok(())
    }

    /// Validates and replaces one member, returning the previous one. The
    /// index is checked first; failure leaves the ensemble unchanged.
    pub fn replace(
        &mut self,
        index: usize,
        mut member: EnsembleMember,
    ) -> Result<EnsembleMember, RealizationError> {
        self.check_index(index)?;
        self.bind(&mut member)?;
        Ok(std::mem::replace(&mut self.members[index], member))
    }

    /// Removes and returns one member, shifting later ones down.
    pub fn remove(&mut self, index: usize) -> Result<EnsembleMember, RealizationError> {
        self.check_index(index)?;
        Ok(self.members.remove(index))
    }

    /// Replaces the positions of every member transactionally, in order.
    /// Other state is kept.
    pub fn replace_positions(
        &mut self,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<(), RealizationError> {
        let positions = positions.into_iter().collect::<Vec<_>>();
        if positions.len() != self.members.len() {
            return Err(RealizationError::ItemCountMismatch {
                expected: self.members.len(),
                actual: positions.len(),
            });
        }
        for values in &positions {
            if values.len() != self.topology.atom_count() {
                return Err(ConformationError::AtomCountMismatch {
                    expected: self.topology.atom_count(),
                    actual: values.len(),
                }
                .into());
            }
        }
        for (member, positions) in self.members.iter_mut().zip(positions) {
            member.positions = positions;
        }
        Ok(())
    }

    /// Copies members in the requested order, sharing the topology. Ranges
    /// and strides can be written as `(start..end).step_by(stride)`;
    /// duplicates and reordering are allowed. Ensemble properties are kept.
    pub fn select(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, RealizationError> {
        let members = indices
            .into_iter()
            .map(|index| {
                self.check_index(index)?;
                Ok(self.members[index].clone())
            })
            .collect::<Result<Vec<_>, RealizationError>>()?;
        Ok(Self {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            members,
        })
    }

    /// Publishes one induced topology subset and projects every member onto
    /// it. Ensemble and per-member owner properties are dropped.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, RealizationError> {
        let subset = self.topology.subset(selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let members = self
            .members
            .iter()
            .map(|member| member.project(&atoms, &bonds))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology: subset.shared_topology(),
            properties: OwnerProperties::new(),
            members,
        })
    }

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once and keeps every member and
    /// the ensemble properties without copying. Failure leaves the ensemble
    /// unchanged. The new snapshot shares the original layout, so selections,
    /// buffers, and prepared potentials stay usable.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// Length-preserving mutable access to the shared topology's static
    /// annotations.
    ///
    /// A topology shared with other owners is first replaced by a snapshot of
    /// it ([`Arc::make_mut`]), so those owners are unaffected. Either way the
    /// layout is kept, so every member, selection, buffer, and prepared
    /// potential remains usable. See [`Topology::properties_mut`].
    pub fn topology_properties_mut(&mut self) -> TopologyPropertiesMut<'_> {
        Arc::make_mut(&mut self.topology).properties_mut()
    }

    /// Consumes the ensemble, transferring topology, ensemble properties, and
    /// every member without copying.
    pub fn into_parts(self) -> (Arc<Topology>, OwnerProperties, Vec<EnsembleMember>) {
        (self.topology, self.properties, self.members)
    }

    pub fn into_items(self) -> Vec<EnsembleMember> {
        self.members
    }

    /// Scales every weight so the weights sum to one, keeping their ratios.
    ///
    /// Fails on an empty ensemble, which has no probabilities, or when a
    /// scaled weight would underflow to zero. Failure changes nothing.
    pub fn normalize_weights(&mut self) -> Result<(), RealizationError> {
        if self.is_empty() {
            return Err(RealizationError::EmptySource);
        }
        // Scale by the largest weight first so the total cannot overflow.
        let maximum = self
            .members
            .iter()
            .map(|member| member.weight)
            .fold(0.0_f64, f64::max);
        let total = self
            .members
            .iter()
            .map(|member| member.weight / maximum)
            .sum::<f64>();
        let normalized = self
            .members
            .iter()
            .map(|member| checked_weight(member.weight / maximum / total))
            .collect::<Result<Vec<_>, _>>()?;
        for (member, weight) in self.members.iter_mut().zip(normalized) {
            member.weight = weight;
        }
        Ok(())
    }

    /// Applies one rigid transform per member transactionally.
    pub(crate) fn transform_members(
        &mut self,
        transforms: &[RigidTransform],
    ) -> Result<(), ConformationError> {
        debug_assert_eq!(transforms.len(), self.members.len());
        let members = self
            .members
            .iter()
            .zip(transforms)
            .map(|(member, transform)| member.transformed(*transform))
            .collect::<Result<Vec<_>, _>>()?;
        self.members = members;
        Ok(())
    }

    fn check_index(&self, index: usize) -> Result<(), RealizationError> {
        if index >= self.members.len() {
            return Err(RealizationError::IndexOutOfRange {
                index,
                len: self.members.len(),
            });
        }
        Ok(())
    }

    fn bind(&self, member: &mut EnsembleMember) -> Result<(), RealizationError> {
        validate_binding(member.atom_count(), &member.properties, &self.topology)?;
        member.properties.bind_bond_rows(self.topology.bond_count());
        Ok(())
    }
}

impl<'a> IntoIterator for &'a Ensemble {
    type Item = EnsembleMemberView<'a>;
    type IntoIter = EnsembleIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// One ensemble member: positions in dense topology atom order, an optional
/// periodic cell, optional crystallographic occupancies and B-factors,
/// realization properties, and a statistical weight.
///
/// The member carries no topology; it is bound when inserted into an
/// [`Ensemble`].
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleMember {
    positions: Positions,
    cell: Option<PeriodicCell>,
    occupancies: AtomScalars,
    b_factors: AtomScalars,
    properties: RealizationProperties,
    weight: f64,
}

impl EnsembleMember {
    /// Creates a non-periodic member without occupancies, B-factors, or
    /// properties, with a finite positive relative weight.
    pub fn new(positions: Positions, weight: f64) -> Result<Self, ConformationError> {
        let mut member = Self::unit(positions);
        member.weight = checked_weight(weight)?;
        Ok(member)
    }

    /// A member of an equally weighted collection.
    fn unit(positions: Positions) -> Self {
        let properties = RealizationProperties::new(positions.len());
        Self {
            positions,
            cell: None,
            occupancies: AtomScalars::default(),
            b_factors: AtomScalars::default(),
            properties,
            weight: 1.0,
        }
    }

    /// An equally weighted member holding a model's realization state; the
    /// model's topology is dropped.
    pub(crate) fn from_model(model: Model) -> Self {
        Self {
            positions: model.positions,
            cell: model.cell,
            occupancies: model.occupancies,
            b_factors: model.b_factors,
            properties: model.properties,
            weight: 1.0,
        }
    }

    /// An equally weighted member without occupancies or B-factors, for
    /// frames reinterpreted as an ensemble.
    pub(crate) fn from_frame_state(
        positions: Positions,
        cell: Option<PeriodicCell>,
        properties: RealizationProperties,
    ) -> Self {
        Self {
            positions,
            cell,
            occupancies: AtomScalars::default(),
            b_factors: AtomScalars::default(),
            properties,
            weight: 1.0,
        }
    }

    /// Returns this member with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// The relative statistical weight; always finite and positive.
    pub const fn weight(&self) -> f64 {
        self.weight
    }

    /// Sets a finite positive relative weight.
    pub fn set_weight(&mut self, weight: f64) -> Result<(), ConformationError> {
        self.weight = checked_weight(weight)?;
        Ok(())
    }

    /// Applies a rigid transform transactionally: positions move and the
    /// cell rotates.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), ConformationError> {
        *self = self.transformed(*transform)?;
        Ok(())
    }

    fn model_view<'a>(&'a self, topology: &'a Arc<Topology>) -> ModelView<'a> {
        ModelView::bound(
            topology,
            &self.positions,
            self.cell.as_ref(),
            &self.occupancies,
            &self.b_factors,
            &self.properties,
        )
    }

    /// Projects every dense row through a subset correspondence.
    fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, ConformationError> {
        Ok(Self {
            positions: self.positions.select_indices(atoms)?,
            cell: self.cell,
            occupancies: self.occupancies.select(atoms),
            b_factors: self.b_factors.select(atoms),
            properties: self.properties.project(atoms, bonds)?,
            weight: self.weight,
        })
    }

    /// Positions move and the cell rotates.
    fn transformed(&self, transform: RigidTransform) -> Result<Self, ConformationError> {
        Ok(Self {
            cell: transform_cell(self.cell, transform)?,
            positions: transform_positions(&self.positions, transform)?,
            occupancies: self.occupancies.clone(),
            b_factors: self.b_factors.clone(),
            properties: self.properties.clone(),
            weight: self.weight,
        })
    }

    fn storage_ref(&self) -> &Self {
        self
    }

    fn storage(&mut self) -> &mut Self {
        self
    }

    realization_readers!();
    site_readers!();
    realization_mutators!();
    site_mutators!();
}

fn checked_weight(weight: f64) -> Result<f64, ConformationError> {
    if weight.is_finite() && weight > 0.0 {
        Ok(weight)
    } else {
        Err(ConformationError::InvalidWeight)
    }
}

/// One topology-bound [`EnsembleMember`] borrowed from an [`Ensemble`].
///
/// Reads member state through `Deref` and converts to a [`ModelView`]
/// through [`AsModelView`] without copying.
#[derive(Debug, Clone, Copy)]
pub struct EnsembleMemberView<'a> {
    topology: &'a Arc<Topology>,
    member: &'a EnsembleMember,
}

impl std::ops::Deref for EnsembleMemberView<'_> {
    type Target = EnsembleMember;

    fn deref(&self) -> &Self::Target {
        self.member
    }
}

impl<'a> EnsembleMemberView<'a> {
    /// Binds a member that already has every row `topology` requires.
    pub fn new(
        topology: &'a Arc<Topology>,
        member: &'a EnsembleMember,
    ) -> Result<Self, ConformationError> {
        validate_binding(member.atom_count(), &member.properties, topology)?;
        if !is_bound(member.atom_count(), &member.properties, topology) {
            return Err(ConformationError::BondCountMismatch {
                expected: topology.bond_count(),
                actual: member.properties.bonds().len(),
            });
        }
        Ok(Self { topology, member })
    }

    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    /// The borrowed member with the view's full lifetime.
    pub fn payload(self) -> &'a EnsembleMember {
        self.member
    }

    /// Borrows the member as the common kernel input, with the view's full
    /// lifetime.
    pub fn as_model_view(self) -> ModelView<'a> {
        self.member.model_view(self.topology)
    }

    /// Materializes the member as an owned model, dropping its weight.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

impl AsModelView for EnsembleMemberView<'_> {
    fn as_model_view(&self) -> ModelView<'_> {
        EnsembleMemberView::as_model_view(*self)
    }
}

/// Dimension-preserving mutable access to one member of an [`Ensemble`].
///
/// Reads use `Deref`. There is deliberately no `DerefMut`; replace a whole
/// member through [`Ensemble::replace`].
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
    /// Sets a finite positive relative weight.
    pub fn set_weight(&mut self, weight: f64) -> Result<(), ConformationError> {
        self.member.set_weight(weight)
    }

    /// Applies a rigid transform transactionally: positions move and the
    /// cell rotates.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), ConformationError> {
        self.member.apply_rigid_transform(transform)
    }

    fn storage(&mut self) -> &mut EnsembleMember {
        self.member
    }

    realization_mutators!();
    site_mutators!();
}

/// Iterator over the topology-bound members of an [`Ensemble`].
#[derive(Debug, Clone)]
pub struct EnsembleIter<'a> {
    topology: &'a Arc<Topology>,
    members: std::slice::Iter<'a, EnsembleMember>,
}

impl<'a> Iterator for EnsembleIter<'a> {
    type Item = EnsembleMemberView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let topology = self.topology;
        self.members
            .next()
            .map(|member| EnsembleMemberView { topology, member })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.members.size_hint()
    }
}

impl ExactSizeIterator for EnsembleIter<'_> {}
