//! Weights of a biased run from its sampling record:
//! [`Weights::from_bias`].

use crate::properties::{PropertyColumn, PropertyKey};
use crate::structure::{Order, Requirement, StructureSetRead};
use crate::units::CANONICAL_ENERGY_UNIT;

use super::{WeightSource, Weights, WeightsError};

/// Why [`Weights::from_bias`] cannot weigh a geometry, so it is no sample of
/// the sampled state.
///
/// The reasons are checked in this order, and a geometry has the first that
/// applies.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Unweighable {
    /// Its step lies before the sampling record's first step.
    BeforeSampling,
    /// It has no value in a bias column.
    MissingBias {
        /// The first bias column, in the record's order, without its value.
        key: PropertyKey,
    },
}

impl Unweighable {
    /// How far a geometry with this reason got through the checks.
    const fn stage(&self) -> u8 {
        match self {
            Self::BeforeSampling => 0,
            Self::MissingBias { .. } => 1,
        }
    }
}

/// A geometry [`Weights::from_bias`] leaves out, and why.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeftOut {
    geometry: usize,
    reason: Unweighable,
}

impl LeftOut {
    /// The geometry's index in the set.
    pub const fn geometry(&self) -> usize {
        self.geometry
    }

    /// Why it cannot be weighed.
    pub const fn reason(&self) -> &Unweighable {
        &self.reason
    }
}

/// The weights of a biased run ([`Weights::from_bias`]), with the
/// geometries left out.
#[derive(Debug, Clone, PartialEq)]
pub struct BiasWeights {
    weights: Weights,
    left_out: Vec<LeftOut>,
}

impl BiasWeights {
    /// The weights: every weighable geometry a sample.
    pub fn weights(&self) -> &Weights {
        &self.weights
    }

    /// The geometries left out, by increasing index, each with its reason.
    pub fn left_out(&self) -> &[LeftOut] {
        &self.left_out
    }

    /// The weights, without the geometries left out.
    pub fn into_weights(self) -> Weights {
        self.weights
    }
}

impl Weights {
    /// The weights of a biased run from the sampling record of `set`: each
    /// geometry's log weight is `β ΣΔV`, the sum of its bias columns over
    /// the sampled state's thermal energy, its exact log importance weight.
    /// The weights describe the sampled state ([`WeightSource::Bias`]).
    /// Without bias columns, every geometry weighs the same.
    ///
    /// A geometry is left out when its step lies before the record's first
    /// step or when it lacks a bias value; the first reason that applies is
    /// its [`Unweighable`]. The potential-energy column is not read.
    ///
    /// ```
    /// use kekule::properties::{PropertyColumn, PropertyKey};
    /// use kekule::statistics::{Unweighable, Weights};
    /// use kekule::structure::{
    ///     Geometry, Order, Positions, Sampling, StructureSet, ThermodynamicState,
    /// };
    /// use kekule::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE};
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let atoms = topology.atom_count();
    /// let geometries = (0..3).map(|_| Geometry::new(Positions::zeros(atoms)));
    /// let mut set = StructureSet::from_geometries(topology, geometries)?;
    /// set.set_order(Some(Order::Step(vec![0, 500, 1_000])))?;
    /// let boost = PropertyKey::new("boost")?;
    /// set.columns_mut().insert(
    ///     boost.clone(),
    ///     PropertyColumn::Real { unit: KILOJOULE_PER_MOLE, values: vec![Some(1.0), Some(2.5), None] },
    /// )?;
    /// let state = ThermodynamicState::canonical(Quantity::new(300.0, KELVIN))?;
    /// set.set_sampling(Some(Sampling::new(state).with_from_step(500).with_bias([boost.clone()])?))?;
    ///
    /// let bias = Weights::from_bias(&set)?;
    /// assert_eq!(bias.weights().samples().collect::<Vec<_>>(), [1]);
    /// let rt = *state.thermal_energy().value();
    /// assert!((bias.weights().log_weight(1).unwrap() - 2.5 / rt).abs() < 1e-12);
    /// assert_eq!(bias.left_out()[0].reason(), &Unweighable::BeforeSampling);
    /// assert_eq!(bias.left_out()[1].reason(), &Unweighable::MissingBias { key: boost });
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// `Missing(Sampling)` for a set without a sampling record,
    /// `Missing(Steps)` when the record's first step is not 0 and the set has
    /// no steps, [`WeightsError::NonEnergyColumn`] for a bias column that
    /// holds no real energies, [`WeightsError::InvalidBias`] for the first
    /// geometry whose bias is not finite, and
    /// [`WeightsError::NoWeighableSample`] when the set has geometries and
    /// none can be weighed, with the reason of the first geometry that got
    /// furthest through the checks.
    pub fn from_bias(set: &(impl StructureSetRead + ?Sized)) -> Result<BiasWeights, WeightsError> {
        let sampling = set
            .sampling()
            .ok_or(WeightsError::Missing(Requirement::Sampling))?;
        let state = *sampling.state();
        let beta = 1.0 / *state.thermal_energy().value();
        let columns = set.columns();
        let bias = sampling
            .bias()
            .iter()
            .map(|key| match columns.and_then(|columns| columns.get(key)) {
                None => Ok((key, None, 1.0)),
                Some(PropertyColumn::Real { unit, values }) => unit
                    .conversion_factor_to(CANONICAL_ENERGY_UNIT)
                    .map(|factor| (key, Some(values.as_slice()), factor))
                    .map_err(|_| WeightsError::NonEnergyColumn { key: key.clone() }),
                Some(_) => Err(WeightsError::NonEnergyColumn { key: key.clone() }),
            })
            .collect::<Result<Vec<_>, _>>()?;
        let steps = if sampling.from_step() == 0 {
            None
        } else {
            let steps = match set.order() {
                Some(Order::Time(axis)) => axis.steps(),
                Some(Order::Step(steps)) => Some(steps.as_slice()),
                _ => None,
            };
            Some(steps.ok_or(WeightsError::Missing(Requirement::Steps))?)
        };
        let mut log = Vec::with_capacity(set.len());
        let mut left_out = Vec::new();
        for geometry in 0..set.len() {
            let before = steps
                .and_then(|steps| steps.get(geometry))
                .is_some_and(|&step| step < sampling.from_step());
            let mut reason = before.then_some(Unweighable::BeforeSampling);
            let mut sum = 0.0;
            if reason.is_none() {
                for (key, values, factor) in &bias {
                    match values.and_then(|values| values.get(geometry).copied().flatten()) {
                        Some(value) => sum += value * factor,
                        None => {
                            reason = Some(Unweighable::MissingBias {
                                key: (*key).clone(),
                            });
                            break;
                        }
                    }
                }
            }
            match reason {
                Some(reason) => {
                    log.push(None);
                    left_out.push(LeftOut { geometry, reason });
                }
                None => {
                    let log_weight = beta * sum;
                    if !log_weight.is_finite() {
                        return Err(WeightsError::InvalidBias { geometry });
                    }
                    log.push(Some(log_weight));
                }
            }
        }
        if left_out.len() == set.len() {
            let furthest = left_out
                .iter()
                .min_by_key(|left_out| std::cmp::Reverse(left_out.reason.stage()));
            if let Some(left_out) = furthest {
                return Err(WeightsError::NoWeighableSample(left_out.reason.clone()));
            }
        }
        Ok(BiasWeights {
            weights: Self::from_parts_unchecked(log, WeightSource::Bias { state }),
            left_out,
        })
    }
}
