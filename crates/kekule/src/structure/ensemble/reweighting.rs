//! Reduced biases of a structure set's geometries, derived from its
//! sampling record against a target temperature, and which geometries can
//! be weighed by them.

use crate::diagnostics::{Countable, IndexOutOfRange};
use crate::properties::{PropertyColumn, PropertyKey, PropertyTable};
use crate::structure::{
    GeometryIndex, Order, Requirement, Sampling, StructureSetError, ThermodynamicState,
};
use crate::units::{Quantity, CANONICAL_ENERGY_UNIT};

use super::distribution::REDUCED_BIAS_LIMIT;
use super::Samples;

/// Why reweighting cannot weigh a geometry, so it is no sample of the
/// sampled state ([`Distribution::weighable`](super::Distribution::weighable)).
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
    /// It has no value in the potential-energy column, which reweighting to
    /// another temperature reads.
    MissingPotentialEnergy,
}

impl Unweighable {
    /// How far a geometry with this reason got through the checks.
    const fn stage(&self) -> u8 {
        match self {
            Self::BeforeSampling => 0,
            Self::MissingBias { .. } => 1,
            Self::MissingPotentialEnergy => 2,
        }
    }
}

/// A proposed sample that reweighting leaves out, and why.
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

/// Which proposed samples a reweighting can weigh
/// ([`Distribution::weighable`](super::Distribution::weighable)): the
/// weighable geometries and the left-out ones, each by increasing index.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeighableSamples {
    samples: Vec<usize>,
    left_out: Vec<LeftOut>,
}

impl WeighableSamples {
    /// The geometries that can be weighed, by increasing index.
    pub fn samples(&self) -> &[usize] {
        &self.samples
    }

    /// The geometries left out, by increasing index, each with its reason.
    pub fn left_out(&self) -> &[LeftOut] {
        &self.left_out
    }
}

/// What [`Ensemble::reweight`](super::Ensemble::reweight) left out: the
/// samples it could not weigh, and where each state went.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reweighting {
    left_out: Vec<LeftOut>,
    states: Vec<Option<usize>>,
}

impl Reweighting {
    /// The samples left out, by increasing geometry index, each with its
    /// reason.
    pub fn left_out(&self) -> &[LeftOut] {
        &self.left_out
    }

    /// The new index of each state before reweighting, in the old order:
    /// `None` for a state whose samples were all left out. The kept states
    /// keep their order.
    pub fn states(&self) -> &[Option<usize>] {
        &self.states
    }
}

/// One energy column the sampling record names, read in
/// [`CANONICAL_ENERGY_UNIT`]: `None` values when the column is absent.
struct EnergyColumn<'a> {
    key: &'a PropertyKey,
    values: Option<&'a [Option<f64>]>,
    factor: f64,
}

impl EnergyColumn<'_> {
    /// The value of geometry `geometry` in kJ/mol, if it has one.
    fn value(&self, geometry: usize) -> Option<f64> {
        self.values
            .and_then(|values| values.get(geometry).copied().flatten())
            .map(|value| value * self.factor)
    }
}

/// The geometries a [`ReducedBias`] weighs among those it was given.
pub(crate) struct Weighed {
    /// The positions of the weighed geometries among those given.
    pub(crate) positions: Vec<usize>,
    /// The weighed geometries.
    pub(crate) samples: Vec<usize>,
    /// Their reduced biases.
    pub(crate) reduced: Vec<f64>,
    /// The geometries left out.
    pub(crate) left_out: Vec<LeftOut>,
}

impl Weighed {
    /// What the query reports.
    pub(crate) fn into_weighable(self) -> WeighableSamples {
        WeighableSamples {
            samples: self.samples,
            left_out: self.left_out,
        }
    }

    /// Fails when geometries were given and none can be weighed, with the
    /// reason of the first geometry that got furthest through the checks.
    pub(crate) fn require_sample(&self) -> Result<(), StructureSetError> {
        if !self.samples.is_empty() {
            return Ok(());
        }
        let furthest = self
            .left_out
            .iter()
            .min_by_key(|left_out| std::cmp::Reverse(left_out.reason.stage()));
        match furthest {
            None => Ok(()),
            Some(left_out) => Err(StructureSetError::NoWeighableSample(
                left_out.reason.clone(),
            )),
        }
    }
}

/// The reduced bias of each geometry of a set at a target temperature,
/// from the set's sampling record:
/// `b = β_s (U + ΣΔV) − β_t U = β_s ΣΔV + (β_s − β_t) U`, with the sampled
/// state's `β_s`, the target's `β_t`, the unbiased potential energy `U`, and
/// the sum of the bias columns `ΣΔV`. At the sampled temperature `b = β_s ΣΔV`,
/// and the potential energy is not read. Pressure is not reweighted: the
/// target keeps the sampled pressure, and the `pV` term is neglected, which
/// is negligible near ambient pressure.
///
/// It is the one rule of which geometries can be weighed: every reweighting
/// and [`Distribution::weighable`](super::Distribution::weighable) go
/// through [`Self::weigh`].
pub(crate) struct ReducedBias<'a> {
    target: ThermodynamicState,
    len: usize,
    from_step: u64,
    steps: Option<&'a [u64]>,
    bias: Vec<EnergyColumn<'a>>,
    potential: Option<EnergyColumn<'a>>,
    sampled_beta: f64,
    beta_difference: f64,
}

impl<'a> ReducedBias<'a> {
    /// The reduced biases of a set of `len` geometries with `sampling`,
    /// `columns`, and `order` at `temperature`, or at the sampled
    /// temperature for `None`.
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
        len: usize,
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
                    factor: unit
                        .conversion_factor_to(CANONICAL_ENERGY_UNIT)
                        .map_err(|_| StructureSetError::SamplingColumn { key: key.clone() })?,
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
            len,
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

    /// The reduced bias of geometry `geometry`, an index of the set, or why
    /// it cannot be weighed.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::InvalidReducedBias`] for a reduced bias that is
    /// not finite or beyond [`REDUCED_BIAS_LIMIT`].
    fn reduced(&self, geometry: usize) -> Result<Result<f64, Unweighable>, StructureSetError> {
        if let Some(steps) = self.steps {
            if steps
                .get(geometry)
                .is_some_and(|&step| step < self.from_step)
            {
                return Ok(Err(Unweighable::BeforeSampling));
            }
        }
        let mut bias = 0.0;
        for column in &self.bias {
            match column.value(geometry) {
                Some(value) => bias += value,
                None => {
                    return Ok(Err(Unweighable::MissingBias {
                        key: column.key.clone(),
                    }))
                }
            }
        }
        let mut reduced = self.sampled_beta * bias;
        if let Some(potential) = &self.potential {
            match potential.value(geometry) {
                Some(value) => reduced += self.beta_difference * value,
                None => return Ok(Err(Unweighable::MissingPotentialEnergy)),
            }
        }
        if !reduced.is_finite() || reduced.abs() > REDUCED_BIAS_LIMIT {
            return Err(StructureSetError::InvalidReducedBias { geometry });
        }
        Ok(Ok(reduced))
    }

    /// Weighs `geometries`, indices of the set, and leaves out those that
    /// cannot be weighed.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::InvalidReducedBias`] for the first geometry
    /// whose reduced bias is not finite or beyond [`REDUCED_BIAS_LIMIT`].
    pub(crate) fn weigh(
        &self,
        geometries: impl IntoIterator<Item = usize>,
    ) -> Result<Weighed, StructureSetError> {
        let mut weighed = Weighed {
            positions: Vec::new(),
            samples: Vec::new(),
            reduced: Vec::new(),
            left_out: Vec::new(),
        };
        for (position, geometry) in geometries.into_iter().enumerate() {
            match self.reduced(geometry)? {
                Ok(reduced) => {
                    weighed.positions.push(position);
                    weighed.samples.push(geometry);
                    weighed.reduced.push(reduced);
                }
                Err(reason) => weighed.left_out.push(LeftOut { geometry, reason }),
            }
        }
        Ok(weighed)
    }

    /// Weighs the geometries `samples` lists.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::UnsortedSamples`] and
    /// [`StructureSetError::IndexOutOfRange`] for a subset that is not
    /// increasing or lists a geometry outside the set, and the errors of
    /// [`Self::weigh`].
    pub(crate) fn weigh_samples(&self, samples: &Samples) -> Result<Weighed, StructureSetError> {
        match samples {
            Samples::All => self.weigh(0..self.len),
            Samples::Subset(indices) => {
                if indices.windows(2).any(|pair| pair[0] >= pair[1]) {
                    return Err(StructureSetError::UnsortedSamples);
                }
                if let Some(&last) = indices.last().filter(|&&last| last >= self.len) {
                    return Err(StructureSetError::IndexOutOfRange(IndexOutOfRange::new(
                        Countable::Geometry,
                        last,
                        self.len,
                    )));
                }
                self.weigh(indices.iter().copied())
            }
        }
    }

    /// The error of [`Distribution::reweighted`](super::Distribution::reweighted)
    /// for a sample it cannot weigh.
    pub(crate) fn error(&self, left_out: &LeftOut) -> StructureSetError {
        let geometry = left_out.geometry;
        match &left_out.reason {
            Unweighable::BeforeSampling => StructureSetError::BeforeSampling { geometry },
            Unweighable::MissingBias { key } => StructureSetError::MissingColumnValue {
                key: key.clone(),
                geometry,
            },
            Unweighable::MissingPotentialEnergy => StructureSetError::MissingColumnValue {
                key: self
                    .potential
                    .as_ref()
                    .expect("only a potential-energy column that is read can lack a value")
                    .key
                    .clone(),
                geometry,
            },
        }
    }
}

/// A [`Reweighting`] report.
pub(crate) fn reweighting(left_out: Vec<LeftOut>, states: Vec<Option<usize>>) -> Reweighting {
    Reweighting { left_out, states }
}
