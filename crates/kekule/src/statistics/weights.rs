//! The weights of a structure set's geometries: [`Weights`].

use crate::diagnostics::{Countable, IndexOutOfRange};
use crate::properties::{PropertyColumn, PropertyKey};
use crate::structure::{StructureSetError, StructureSetRead, ThermodynamicState};
use crate::units::{Quantity, CANONICAL_ENERGY_UNIT};

use super::log_sum_exp;

mod bias;
mod error;

pub use bias::{BiasWeights, LeftOut, Unweighable};
pub use error::WeightsError;

/// Where the log weights of [`Weights`] come from.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum WeightSource {
    /// Every sample weighs the same ([`Weights::equal`],
    /// [`Weights::equal_over`]).
    Equal,
    /// Log weights the caller gave ([`Weights::given`]), with the
    /// temperature they describe, if stated.
    Given {
        /// The temperature, in
        /// [`CANONICAL_TEMPERATURE_UNIT`](crate::units::CANONICAL_TEMPERATURE_UNIT).
        temperature: Option<Quantity<f64>>,
    },
    /// Boltzmann factors of an energy column at a temperature
    /// ([`Weights::boltzmann`]).
    Boltzmann {
        /// The temperature, in
        /// [`CANONICAL_TEMPERATURE_UNIT`](crate::units::CANONICAL_TEMPERATURE_UNIT).
        temperature: Quantity<f64>,
        /// The set column of the energies.
        energy: PropertyKey,
    },
    /// The bias of a biased run, from the set's sampling record
    /// ([`Weights::from_bias`]).
    Bias {
        /// The sampled state, which the weights describe.
        state: ThermodynamicState,
    },
}

/// One optional natural-logarithm weight per geometry of a structure set,
/// with its source: which geometries are samples and what each weighs.
///
/// A weight is inferred, never recorded, so it is a value beside its set and
/// never a field of it; one set often carries several readings. The weights
/// are indexed by geometry: a geometry without a log weight is not a sample.
/// A host that selects geometries of the set selects them in the weights too
/// ([`Self::select`]).
///
/// Only ratios of weights are meaningful, so they are natural logarithms,
/// which any finite value is: weights from a biased run span hundreds of
/// orders of magnitude, beyond what a linear weight can hold. Every
/// population is exact, the log-sum-exp of its samples' log weights, so the
/// populations of a group's parts add up to the group's. Probabilities are
/// derived on request and never stored.
///
/// ```
/// use kekule::statistics::Weights;
///
/// // Geometries 1 and 2 of three, the first weighing three times the second.
/// let weights = Weights::given(vec![None, Some(3.0_f64.ln()), Some(0.0)], None)?;
/// assert_eq!(weights.sample_count(), 2);
/// let probabilities = weights.probabilities()?;
/// assert_eq!(probabilities[0], 0.0);
/// assert!((probabilities[1] - 0.75).abs() < 1e-15);
/// assert_eq!(weights.ranked(), [1, 2]);
/// # Ok::<(), kekule::statistics::WeightsError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Weights {
    log: Vec<Option<f64>>,
    source: WeightSource,
}

impl Weights {
    /// Every geometry of a set of `len` as a sample, each weighing the same:
    /// the frames of an equilibrium trajectory or the models of an NMR
    /// bundle.
    pub fn equal(len: usize) -> Self {
        Self {
            log: vec![Some(0.0); len],
            source: WeightSource::Equal,
        }
    }

    /// The listed geometries of a set of `len` as samples, each weighing the
    /// same; the others are not samples. A geometry listed twice is one
    /// sample.
    ///
    /// # Errors
    ///
    /// [`WeightsError::IndexOutOfRange`] for the first geometry past `len`.
    pub fn equal_over(
        len: usize,
        geometries: impl IntoIterator<Item = usize>,
    ) -> Result<Self, WeightsError> {
        let mut log = vec![None; len];
        for geometry in geometries {
            *log.get_mut(geometry)
                .ok_or_else(|| out_of_range(geometry, len))? = Some(0.0);
        }
        Ok(Self {
            log,
            source: WeightSource::Equal,
        })
    }

    /// One natural-logarithm weight per geometry, `None` for a geometry that
    /// is not a sample, at `temperature` when one is stated, such as the
    /// weights of a reweighting done elsewhere.
    ///
    /// # Errors
    ///
    /// [`WeightsError::NonFiniteLogWeight`] for the first log weight that is
    /// not finite, [`WeightsError::Unit`] for a temperature that is not a
    /// temperature, and [`WeightsError::InvalidTemperature`] for one that is
    /// not finite and positive.
    pub fn given(
        log_weights: Vec<Option<f64>>,
        temperature: Option<Quantity<f64>>,
    ) -> Result<Self, WeightsError> {
        if let Some(geometry) = log_weights
            .iter()
            .position(|log_weight| log_weight.is_some_and(|value| !value.is_finite()))
        {
            return Err(WeightsError::NonFiniteLogWeight { geometry });
        }
        let temperature = temperature
            .map(|temperature| canonical(temperature).map(|state| state.temperature()))
            .transpose()?;
        Ok(Self {
            log: log_weights,
            source: WeightSource::Given { temperature },
        })
    }

    /// Boltzmann weights of the real energy column `energy` of `set` at
    /// `temperature`: `ln w = −(E − E_min) / RT`, relative to the lowest
    /// energy, so no weight underflows in log space. A geometry whose cell
    /// is missing is not a sample.
    ///
    /// The energies are of one scale, such as the minimized energies of
    /// conformers; samples drawn from a simulation already follow their
    /// distribution and weigh the same.
    ///
    /// ```
    /// use kekule::properties::{PropertyColumn, PropertyKey};
    /// use kekule::statistics::Weights;
    /// use kekule::structure::{Geometry, Positions, StructureSet};
    /// use kekule::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE};
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let atoms = topology.atom_count();
    /// let geometries = (0..3).map(|_| Geometry::new(Positions::zeros(atoms)));
    /// let mut set = StructureSet::from_geometries(topology, geometries)?;
    /// let energy = PropertyKey::new("energy")?;
    /// set.columns_mut().insert(
    ///     energy.clone(),
    ///     PropertyColumn::Real { unit: KILOJOULE_PER_MOLE, values: vec![Some(0.0), None, Some(5.0)] },
    /// )?;
    /// let weights = Weights::boltzmann(&set, &energy, Quantity::new(300.0, KELVIN))?;
    /// assert_eq!(weights.log_weight(0), Some(0.0));
    /// assert_eq!(weights.log_weight(1), None);
    /// let rt = 8.314_462_618_153_24e-3 * 300.0;
    /// assert!((weights.log_weight(2).unwrap() + 5.0 / rt).abs() < 1e-12);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`WeightsError::MissingColumn`] when the set has no column `energy`,
    /// [`WeightsError::NonEnergyColumn`] when it holds no real energies, the
    /// temperature errors of [`Self::given`], and
    /// [`WeightsError::NonFiniteLogWeight`] for the first energy so far
    /// above the lowest at this temperature that its log weight is not
    /// finite.
    pub fn boltzmann(
        set: &(impl StructureSetRead + ?Sized),
        energy: &PropertyKey,
        temperature: Quantity<f64>,
    ) -> Result<Self, WeightsError> {
        let state = canonical(temperature)?;
        let thermal_energy = *state.thermal_energy().value();
        let column = set
            .columns()
            .and_then(|columns| columns.get(energy))
            .ok_or_else(|| WeightsError::MissingColumn {
                key: energy.clone(),
            })?;
        let PropertyColumn::Real { unit, values } = column else {
            return Err(WeightsError::NonEnergyColumn {
                key: energy.clone(),
            });
        };
        let factor = unit
            .conversion_factor_to(CANONICAL_ENERGY_UNIT)
            .map_err(|_| WeightsError::NonEnergyColumn {
                key: energy.clone(),
            })?;
        let energies = values
            .iter()
            .map(|value| value.map(|value| value * factor))
            .collect::<Vec<_>>();
        let minimum = energies
            .iter()
            .flatten()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let log = energies
            .iter()
            .map(|energy| energy.map(|energy| -(energy - minimum) / thermal_energy))
            .collect::<Vec<_>>();
        if let Some(geometry) = log
            .iter()
            .position(|log_weight| log_weight.is_some_and(|value| !value.is_finite()))
        {
            return Err(WeightsError::NonFiniteLogWeight { geometry });
        }
        Ok(Self {
            log,
            source: WeightSource::Boltzmann {
                temperature: state.temperature(),
                energy: energy.clone(),
            },
        })
    }

    /// These weights with only the listed geometries left as samples, each
    /// keeping its log weight, and the source kept: the frames after
    /// equilibration, or those a filter keeps.
    ///
    /// # Errors
    ///
    /// [`WeightsError::IndexOutOfRange`] for the first geometry outside the
    /// weights.
    pub fn restricted(
        &self,
        geometries: impl IntoIterator<Item = usize>,
    ) -> Result<Self, WeightsError> {
        let mut log = vec![None; self.len()];
        for geometry in geometries {
            let geometry = self.check(geometry)?;
            log[geometry] = self.log[geometry];
        }
        Ok(Self {
            log,
            source: self.source.clone(),
        })
    }

    /// The weights of a set of the geometries `indices`, in that order, as
    /// [`StructureSet::select`](crate::structure::StructureSet::select) copies
    /// them: each keeps its log weight or stays no sample, and the source is
    /// kept. A geometry listed twice is two samples.
    ///
    /// # Errors
    ///
    /// [`WeightsError::IndexOutOfRange`] for the first index outside the
    /// weights.
    pub fn select(&self, indices: impl IntoIterator<Item = usize>) -> Result<Self, WeightsError> {
        let log = indices
            .into_iter()
            .map(|geometry| self.check(geometry).map(|geometry| self.log[geometry]))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            log,
            source: self.source.clone(),
        })
    }

    /// These weights in a set of `len` geometries whose first ones are those
    /// they weigh, such as a file that grew after they were made: the new
    /// geometries are not samples.
    ///
    /// # Errors
    ///
    /// [`WeightsError::CoverageMismatch`] when the weights cover more than
    /// `len` geometries.
    pub fn extended_to(mut self, len: usize) -> Result<Self, WeightsError> {
        if len < self.len() {
            return Err(WeightsError::CoverageMismatch {
                covered: self.len(),
                geometries: len,
            });
        }
        self.log.resize(len, None);
        Ok(self)
    }

    /// The number of geometries the weights cover.
    pub fn len(&self) -> usize {
        self.log.len()
    }

    /// Whether the weights cover no geometry.
    pub fn is_empty(&self) -> bool {
        self.log.is_empty()
    }

    /// The number of samples.
    pub fn sample_count(&self) -> usize {
        self.log.iter().flatten().count()
    }

    /// The samples, by increasing geometry index.
    pub fn samples(&self) -> impl Iterator<Item = usize> + '_ {
        self.log
            .iter()
            .enumerate()
            .filter(|(_, log_weight)| log_weight.is_some())
            .map(|(geometry, _)| geometry)
    }

    /// The natural-logarithm weight of geometry `geometry`, or `None` when it
    /// is not a sample or lies outside the weights.
    pub fn log_weight(&self, geometry: usize) -> Option<f64> {
        self.log.get(geometry).copied().flatten()
    }

    /// Where the log weights come from.
    pub fn source(&self) -> &WeightSource {
        &self.source
    }

    /// The temperature the weights describe, in
    /// [`CANONICAL_TEMPERATURE_UNIT`](crate::units::CANONICAL_TEMPERATURE_UNIT):
    /// a Boltzmann weighting's, the sampled state's of a biased run, or the
    /// one given; `None` for equal weights and given weights without one.
    pub fn temperature(&self) -> Option<Quantity<f64>> {
        match &self.source {
            WeightSource::Equal => None,
            WeightSource::Given { temperature } => *temperature,
            WeightSource::Boltzmann { temperature, .. } => Some(*temperature),
            WeightSource::Bias { state } => Some(state.temperature()),
        }
    }

    /// The probability of each geometry, in geometry order: its weight
    /// divided by the total weight, and `0.0` for a geometry that is not a
    /// sample. A probability below the smallest positive `f64` is `0.0`.
    ///
    /// # Errors
    ///
    /// [`WeightsError::NoSamples`] without samples.
    pub fn probabilities(&self) -> Result<Vec<f64>, WeightsError> {
        let maximum = self.maximum().ok_or(WeightsError::NoSamples)?;
        let terms = self
            .log
            .iter()
            .map(|log_weight| log_weight.map_or(0.0, |value| (value - maximum).exp()))
            .collect::<Vec<_>>();
        let total = terms.iter().sum::<f64>();
        Ok(terms.into_iter().map(|term| term / total).collect())
    }

    /// Kish's effective sample size, `(Σw)² / Σw²` over the weights
    /// `w = exp(ℓ)` of the samples: the number of equally weighted samples
    /// that would give the same statistical precision. It is the sample
    /// count for equal weights and approaches one when a single sample
    /// dominates, so it shows how far a reweighting degenerates.
    ///
    /// # Errors
    ///
    /// [`WeightsError::NoSamples`] without samples.
    pub fn effective_sample_size(&self) -> Result<f64, WeightsError> {
        if self.maximum().is_none() {
            return Err(WeightsError::NoSamples);
        }
        Ok(kish(self.log.iter().flatten().copied()))
    }

    /// The natural logarithm of the probability of the listed geometries
    /// together: the log-sum-exp of their log weights relative to the total.
    /// Geometries that are not samples add nothing, and a geometry listed
    /// twice counts twice.
    ///
    /// `None` when a geometry lies outside the weights or none is a sample.
    pub fn log_population(&self, geometries: &[usize]) -> Option<f64> {
        let shift = self.maximum()?;
        let mut values = Vec::with_capacity(geometries.len());
        for &geometry in geometries {
            if let Some(log_weight) = *self.log.get(geometry)? {
                values.push(log_weight - shift);
            }
        }
        if values.is_empty() {
            return None;
        }
        Some(log_sum_exp(values.into_iter()) - self.shifted_total(shift))
    }

    /// The samples by descending weight; equal weights keep geometry order.
    pub fn ranked(&self) -> Vec<usize> {
        let mut samples = self.samples().collect::<Vec<_>>();
        let key = |geometry: usize| self.log[geometry].unwrap_or(f64::NEG_INFINITY);
        samples.sort_by(|&a, &b| key(b).total_cmp(&key(a)));
        samples
    }

    /// The fewest samples, by descending weight, whose weights together
    /// reach `fraction` of the total, in [`Self::ranked`] order. `fraction`
    /// lies in `(0, 1]`; `1` takes every sample, including those too light
    /// to change a floating-point total.
    ///
    /// # Errors
    ///
    /// [`WeightsError::InvalidFraction`] for a fraction outside `(0, 1]` and
    /// [`WeightsError::NoSamples`] without samples.
    pub fn covering(&self, fraction: f64) -> Result<Vec<usize>, WeightsError> {
        check_fraction(fraction)?;
        let probabilities = self.probabilities()?;
        let mut ranked = self.ranked();
        let shares = ranked
            .iter()
            .map(|&geometry| probabilities[geometry])
            .collect::<Vec<_>>();
        ranked.truncate(covered_count(fraction, &shares));
        Ok(ranked)
    }

    /// Fails unless the weights cover `len` geometries.
    pub(crate) fn check_len(&self, len: usize) -> Result<(), WeightsError> {
        if self.len() != len {
            return Err(WeightsError::CoverageMismatch {
                covered: self.len(),
                geometries: len,
            });
        }
        Ok(())
    }

    /// The log weight of every geometry, `None` off the samples.
    pub(crate) fn log_weights(&self) -> &[Option<f64>] {
        &self.log
    }

    /// Weights whose log weights the caller has checked to be finite.
    pub(crate) fn from_parts_unchecked(log: Vec<Option<f64>>, source: WeightSource) -> Self {
        debug_assert!(log.iter().flatten().all(|value| value.is_finite()));
        Self { log, source }
    }

    /// The largest log weight, or `None` without samples.
    pub(crate) fn maximum(&self) -> Option<f64> {
        self.log.iter().flatten().copied().reduce(f64::max)
    }

    /// `ln Σ exp(ℓ − shift)` over every sample.
    fn shifted_total(&self, shift: f64) -> f64 {
        log_sum_exp(self.log.iter().flatten().map(move |value| value - shift))
    }

    fn check(&self, geometry: usize) -> Result<usize, WeightsError> {
        if geometry >= self.len() {
            return Err(out_of_range(geometry, self.len()));
        }
        Ok(geometry)
    }
}

/// Fails unless `fraction` lies in `(0, 1]`.
pub(crate) fn check_fraction(fraction: f64) -> Result<(), WeightsError> {
    if !(fraction > 0.0 && fraction <= 1.0) {
        return Err(WeightsError::InvalidFraction);
    }
    Ok(())
}

/// How many of some shares, in descending order, reach `fraction` of their
/// total: all of them for `1`.
pub(crate) fn covered_count(fraction: f64, shares: &[f64]) -> usize {
    if fraction == 1.0 {
        return shares.len();
    }
    // Sum in ranked order so the running total ends exactly at the total.
    let target = fraction * shares.iter().sum::<f64>();
    let mut covered = 0.0;
    shares
        .iter()
        .position(|share| {
            covered += share;
            covered >= target
        })
        .map_or(shares.len(), |position| position + 1)
}

/// The canonical state at `temperature`, with its errors as weights' errors.
pub(crate) fn canonical(temperature: Quantity<f64>) -> Result<ThermodynamicState, WeightsError> {
    ThermodynamicState::canonical(temperature).map_err(|error| match error {
        StructureSetError::Unit(error) => WeightsError::Unit(error),
        _ => WeightsError::InvalidTemperature,
    })
}

fn out_of_range(geometry: usize, len: usize) -> WeightsError {
    WeightsError::IndexOutOfRange(IndexOutOfRange::new(Countable::Geometry, geometry, len))
}

/// Kish's effective sample size `(Σw)² / Σw²` of the weights `w = exp(ℓ)`
/// for some finite log weights, at least one; exactly their count when they
/// are equal.
pub(crate) fn kish(log_weights: impl Iterator<Item = f64> + Clone) -> f64 {
    // Shifted by the largest log weight, the ratio is unchanged and the
    // largest term is one, so neither sum overflows or cancels.
    let maximum = log_weights.clone().fold(f64::NEG_INFINITY, f64::max);
    let (total, squares) = log_weights.fold((0.0, 0.0), |(total, squares), log_weight| {
        let weight = (log_weight - maximum).exp();
        (total + weight, squares + weight * weight)
    });
    total * (total / squares)
}
