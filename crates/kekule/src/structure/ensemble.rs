use std::sync::Arc;

use crate::core::Molecule;
use crate::geometry::{PeriodicCell, RigidTransform};
use crate::properties::RealizationProperties;
use crate::topology::Topology;

use super::realizations::sealed::Payload;
use super::realizations::{realization_collection, single_molecule_topology, RealizationStore};
use super::scalars::AtomScalars;
use super::state::{
    realization_mutators, realization_readers, site_mutators, site_readers, transform_cell,
    transform_positions,
};
use super::{
    ConformationError, Model, ModelView, Positions, Realization, RealizationError, RealizationMut,
    RealizationView,
};

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
/// Item editors preserve dimensions and cannot overwrite a whole member:
///
/// ```compile_fail,E0594
/// use kekule::structure::{Ensemble, EnsembleMember, Positions};
/// fn overwrite(ensemble: &mut Ensemble) {
///     *ensemble.get_mut(0).unwrap() = EnsembleMember::new(Positions::zeros(1), 1.0).unwrap();
/// }
/// ```
#[derive(Debug, Clone)]
pub struct Ensemble {
    store: RealizationStore<EnsembleMember>,
}

/// One topology-bound [`EnsembleMember`] borrowed from an [`Ensemble`].
pub type EnsembleMemberView<'a> = RealizationView<'a, EnsembleMember>;

/// Dimension-preserving editor for one [`EnsembleMember`] of an [`Ensemble`].
pub type EnsembleMemberMut<'a> = RealizationMut<'a, EnsembleMember>;

realization_collection!(
    Ensemble,
    EnsembleMember,
    EnsembleMemberView<'_>,
    EnsembleMemberMut<'_>,
    "member"
);

impl Ensemble {
    /// Collects models that share one topology layout as equally weighted
    /// members (weight `1.0`), keeping the first model's snapshot. Positions,
    /// cells, occupancies, B-factors, and properties move without copying.
    pub fn from_models(models: impl IntoIterator<Item = Model>) -> Result<Self, RealizationError> {
        Ok(Self {
            store: RealizationStore::from_models(models, EnsembleMember::from_model)?,
        })
    }

    /// Builds equally weighted conformers of one molecule from dense
    /// positions in molecule atom order.
    pub fn from_molecule_positions(
        molecule: Molecule,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<Self, RealizationError> {
        Self::from_items(
            single_molecule_topology(molecule)?,
            positions.into_iter().map(EnsembleMember::unit),
        )
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
            .iter()
            .map(|member| member.weight())
            .fold(0.0_f64, f64::max);
        let total = self
            .iter()
            .map(|member| member.weight() / maximum)
            .sum::<f64>();
        let normalized = self
            .iter()
            .map(|member| checked_weight(member.weight() / maximum / total))
            .collect::<Result<Vec<_>, _>>()?;
        for (member, weight) in self.store.items_mut().iter_mut().zip(normalized) {
            member.weight = weight;
        }
        Ok(())
    }

    pub(crate) fn store(&self) -> &RealizationStore<EnsembleMember> {
        &self.store
    }

    pub(crate) fn store_mut(&mut self) -> &mut RealizationStore<EnsembleMember> {
        &mut self.store
    }

    pub(crate) const fn from_store(store: RealizationStore<EnsembleMember>) -> Self {
        Self { store }
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
    pub(crate) fn unit(positions: Positions) -> Self {
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
}

impl EnsembleMember {
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

impl Payload for EnsembleMember {
    fn row_count(&self) -> usize {
        self.positions.len()
    }

    fn realization_properties(&self) -> &RealizationProperties {
        &self.properties
    }

    fn positions_storage_mut(&mut self) -> &mut Positions {
        &mut self.positions
    }

    fn properties_storage_mut(&mut self) -> &mut RealizationProperties {
        &mut self.properties
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
}

impl Realization for EnsembleMember {}

impl RealizationMut<'_, EnsembleMember> {
    /// Sets a finite positive relative weight.
    pub fn set_weight(&mut self, weight: f64) -> Result<(), ConformationError> {
        self.item.set_weight(weight)
    }

    fn storage(&mut self) -> &mut EnsembleMember {
        self.item
    }

    realization_mutators!();
    site_mutators!();
}
