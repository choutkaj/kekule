use crate::core::Molecule;
use crate::geometry::RigidTransform;

use super::realizations::sealed::Payload;
use super::realizations::{realization_collection, single_molecule_topology, RealizationStore};
use super::{
    Conformation, ConformationError, ConformationMut, Model, Positions, Realization,
    RealizationError, RealizationMut, RealizationView,
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
    /// members (weight `1.0`), keeping the first model's snapshot. Model
    /// conformations move without copying.
    pub fn from_models(models: impl IntoIterator<Item = Model>) -> Result<Self, RealizationError> {
        Ok(Self {
            store: RealizationStore::from_models(models, EnsembleMember::unit)?,
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
            positions
                .into_iter()
                .map(|positions| EnsembleMember::unit(positions.into())),
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

/// One ensemble member: a [`Conformation`] and its statistical weight.
///
/// The member carries no topology; it is bound when inserted into an
/// [`Ensemble`]. Conformation state reads through `Deref`.
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleMember {
    conformation: Conformation,
    weight: f64,
}

impl std::ops::Deref for EnsembleMember {
    type Target = Conformation;

    fn deref(&self) -> &Self::Target {
        &self.conformation
    }
}

impl EnsembleMember {
    /// Creates a member with a finite positive relative weight; accepts
    /// [`Conformation`] or [`super::Positions`].
    pub fn new(
        conformation: impl Into<Conformation>,
        weight: f64,
    ) -> Result<Self, ConformationError> {
        Ok(Self {
            conformation: conformation.into(),
            weight: checked_weight(weight)?,
        })
    }

    /// A member of an equally weighted collection.
    pub(crate) fn unit(conformation: Conformation) -> Self {
        Self {
            conformation,
            weight: 1.0,
        }
    }

    pub fn conformation(&self) -> &Conformation {
        &self.conformation
    }

    pub fn conformation_mut(&mut self) -> ConformationMut<'_> {
        ConformationMut::new(&mut self.conformation)
    }

    pub fn into_conformation(self) -> Conformation {
        self.conformation
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

fn checked_weight(weight: f64) -> Result<f64, ConformationError> {
    if weight.is_finite() && weight > 0.0 {
        Ok(weight)
    } else {
        Err(ConformationError::InvalidWeight)
    }
}

impl Payload for EnsembleMember {
    fn conformation_storage(&mut self) -> &mut Conformation {
        &mut self.conformation
    }

    fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, ConformationError> {
        Ok(Self {
            conformation: self.conformation.project(atoms, bonds)?,
            weight: self.weight,
        })
    }

    fn transformed(&self, transform: RigidTransform) -> Result<Self, ConformationError> {
        Ok(Self {
            conformation: self.conformation.transformed(transform)?,
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
}
