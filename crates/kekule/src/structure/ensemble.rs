use crate::geometry::RigidTransform;

use super::realizations::sealed::Payload;
use super::{
    Conformation, ConformationError, ConformationMut, Realization, RealizationError,
    RealizationMut, RealizationView, Realizations,
};

/// A finite stable-order collection of non-temporal realizations, such as
/// conformers or alternate experimental models, of one shared topology.
///
/// Use [`super::Trajectory`] when item order has temporal meaning.
pub type Ensemble = Realizations<EnsembleMember>;

/// One topology-bound [`EnsembleMember`] borrowed from an [`Ensemble`].
pub type EnsembleMemberView<'a> = RealizationView<'a, EnsembleMember>;

/// Dimension-preserving editor for one [`EnsembleMember`] of an [`Ensemble`].
pub type EnsembleMemberMut<'a> = RealizationMut<'a, EnsembleMember>;

/// One non-temporal ensemble member: a [`Conformation`] plus an optional
/// statistical weight.
///
/// The member carries no topology; it is bound when inserted into an
/// [`Ensemble`]. Conformation state reads through `Deref`.
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleMember {
    conformation: Conformation,
    weight: Option<f64>,
}

impl From<Conformation> for EnsembleMember {
    fn from(conformation: Conformation) -> Self {
        Self {
            conformation,
            weight: None,
        }
    }
}

impl std::ops::Deref for EnsembleMember {
    type Target = Conformation;

    fn deref(&self) -> &Self::Target {
        &self.conformation
    }
}

impl EnsembleMember {
    /// Creates an unweighted member; accepts [`Conformation`] or
    /// [`super::Positions`].
    pub fn new(conformation: impl Into<Conformation>) -> Self {
        Self::from(conformation.into())
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

    pub const fn weight(&self) -> Option<f64> {
        self.weight
    }

    /// Sets a finite non-negative weight, or clears it.
    pub fn set_weight(&mut self, weight: Option<f64>) -> Result<(), ConformationError> {
        if weight.is_some_and(|weight| !weight.is_finite() || weight < 0.0) {
            return Err(ConformationError::InvalidWeight);
        }
        self.weight = weight;
        Ok(())
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

impl EnsembleMember {
    /// Applies a rigid transform transactionally: positions move; cells,
    /// velocities, and forces rotate.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), ConformationError> {
        *self = self.transformed(*transform)?;
        Ok(())
    }
}

impl RealizationMut<'_, EnsembleMember> {
    /// Sets a finite non-negative weight, or clears it.
    pub fn set_weight(&mut self, weight: Option<f64>) -> Result<(), ConformationError> {
        self.item.set_weight(weight)
    }
}

impl Ensemble {
    /// Scales every weight so they sum to one. Every member needs a weight.
    pub fn normalize_weights(&mut self) -> Result<(), RealizationError> {
        if self.is_empty() {
            return Err(RealizationError::EmptySource);
        }
        let mut total = 0.0;
        for (index, member) in self.iter().enumerate() {
            total += member
                .weight()
                .ok_or(RealizationError::MissingWeight { member: index })?;
        }
        if !total.is_finite() || total <= 0.0 {
            return Err(RealizationError::ZeroTotalWeight);
        }
        for index in 0..self.len() {
            let mut member = self.get_mut(index).expect("index is in range");
            let weight = member.weight().map(|weight| weight / total);
            member.set_weight(weight)?;
        }
        Ok(())
    }
}
