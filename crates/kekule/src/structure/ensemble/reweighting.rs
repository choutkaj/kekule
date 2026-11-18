//! Reduced biases of a structure set's geometries, derived from its
//! sampling record against a target temperature.

use crate::properties::{PropertyColumn, PropertyKey, PropertyTable};
use crate::structure::{
    GeometryIndex, Order, Requirement, Sampling, StructureSetError, ThermodynamicState,
};
use crate::units::{Quantity, CANONICAL_ENERGY_UNIT};

use super::distribution::REDUCED_BIAS_LIMIT;

/// One energy column the sampling record names, read in
/// [`CANONICAL_ENERGY_UNIT`]: `None` values when the column is absent.
struct EnergyColumn<'a> {
    key: &'a PropertyKey,
    values: Option<&'a [Option<f64>]>,
    factor: f64,
}

impl EnergyColumn<'_> {
    /// The value of geometry `geometry` in kJ/mol.
    fn value(&self, geometry: usize) -> Result<f64, StructureSetError> {
        self.values
            .and_then(|values| values.get(geometry).copied().flatten())
            .map(|value| value * self.factor)
            .ok_or_else(|| StructureSetError::MissingColumnValue {
                key: self.key.clone(),
                geometry,
            })
    }
}

/// The reduced bias of each geometry of a set at a target temperature,
/// from the set's sampling record:
/// `b = β_s (U + ΣΔV) − β_t U = β_s ΣΔV + (β_s − β_t) U`, with the sampled
/// state's `β_s`, the target's `β_t`, the unbiased potential energy `U`, and
/// the sum of the bias columns `ΣΔV`. At the sampled temperature `b = β_s ΣΔV`,
/// and the potential energy is not read. Pressure is not reweighted: the
/// target keeps the sampled pressure, and no `pV` term enters.
pub(crate) struct ReducedBias<'a> {
    target: ThermodynamicState,
    from_step: u64,
    steps: Option<&'a [u64]>,
    bias: Vec<EnergyColumn<'a>>,
    potential: Option<EnergyColumn<'a>>,
    sampled_beta: f64,
    beta_difference: f64,
}

impl<'a> ReducedBias<'a> {
    /// The reduced biases of a set with `sampling`, `columns`, and `order`
    /// at `temperature`, or at the sampled temperature for `None`.
    ///
    /// # Errors
    ///
    /// `Missing(Sampling)` without a sampling record, `Missing(PotentialEnergy)`
    /// for another temperature when the record names no potential-energy
    /// column, `Missing(Steps)` when the record's first step is not 0 and the
    /// set has no steps, the errors of
    /// [`ThermodynamicState::with_temperature`], and
    /// [`StructureSetError::SamplingColumn`] for a named column that holds no
    /// real energies.
    pub(crate) fn new(
        sampling: Option<&'a Sampling>,
        columns: Option<&'a PropertyTable<GeometryIndex>>,
        order: Option<&'a Order>,
        temperature: Option<Quantity<f64>>,
    ) -> Result<Self, StructureSetError> {
        let sampling = sampling.ok_or(Requirement::Sampling)?;
        let sampled = *sampling.state();
        let target = match temperature {
            None => sampled,
            Some(temperature) => sampled.with_temperature(temperature)?,
        };
        let sampled_beta = 1.0 / *sampled.thermal_energy().value();
        let target_beta = 1.0 / *target.thermal_energy().value();
        let column = |key: &'a PropertyKey| -> Result<EnergyColumn<'a>, StructureSetError> {
            match columns.and_then(|columns| columns.get(key)) {
                None => Ok(EnergyColumn {
                    key,
                    values: None,
                    factor: 1.0,
                }),
                Some(PropertyColumn::Real { unit, values }) => Ok(EnergyColumn {
                    key,
                    values: Some(values),
                    factor: unit.conversion_factor_to(CANONICAL_ENERGY_UNIT)?,
                }),
                Some(_) => Err(StructureSetError::SamplingColumn { key: key.clone() }),
            }
        };
        let potential = if target.temperature().value() == sampled.temperature().value() {
            None
        } else {
            let key = sampling
                .potential_energy()
                .ok_or(Requirement::PotentialEnergy)?;
            Some(column(key)?)
        };
        let bias = sampling
            .bias()
            .iter()
            .map(column)
            .collect::<Result<Vec<_>, _>>()?;
        let steps = if sampling.from_step() == 0 {
            None
        } else {
            let steps = match order {
                Some(Order::Time(axis)) => axis.steps(),
                Some(Order::Step(steps)) => Some(steps.as_slice()),
                _ => None,
            };
            Some(steps.ok_or(Requirement::Steps)?)
        };
        Ok(Self {
            target,
            from_step: sampling.from_step(),
            steps,
            bias,
            potential,
            sampled_beta,
            beta_difference: sampled_beta - target_beta,
        })
    }

    /// The target: the sampled state at the target temperature.
    pub(crate) const fn target(&self) -> ThermodynamicState {
        self.target
    }

    /// The reduced bias of geometry `geometry`, an index of the set.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::BeforeSampling`] for a geometry whose step lies
    /// before the record's first step,
    /// [`StructureSetError::MissingColumnValue`] for a missing bias or
    /// potential-energy cell, and [`StructureSetError::InvalidReducedBias`]
    /// for a reduced bias that is not finite or beyond
    /// [`REDUCED_BIAS_LIMIT`].
    pub(crate) fn reduced(&self, geometry: usize) -> Result<f64, StructureSetError> {
        if let Some(steps) = self.steps {
            if steps
                .get(geometry)
                .is_some_and(|&step| step < self.from_step)
            {
                return Err(StructureSetError::BeforeSampling { geometry });
            }
        }
        let mut bias = 0.0;
        for column in &self.bias {
            bias += column.value(geometry)?;
        }
        let mut reduced = self.sampled_beta * bias;
        if let Some(potential) = &self.potential {
            reduced += self.beta_difference * potential.value(geometry)?;
        }
        if !reduced.is_finite() || reduced.abs() > REDUCED_BIAS_LIMIT {
            return Err(StructureSetError::InvalidReducedBias { geometry });
        }
        Ok(reduced)
    }
}
