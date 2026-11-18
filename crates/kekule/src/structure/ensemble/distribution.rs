//! The distribution of an [`Ensemble`](super::Ensemble): which geometries
//! are samples, their weights, and their grouping into states.

use crate::diagnostics::{Countable, IndexOutOfRange, LengthMismatch};
use crate::properties::OwnerProperties;
use crate::statistics::{log_sum_exp, Estimator};
use crate::structure::{
    Partition, Requirement, StructureSetError, StructureSetRead, ThermodynamicState,
};
use crate::units::{Quantity, CANONICAL_ENERGY_UNIT};
use std::sync::OnceLock;

use super::grouping::Grouping;
use super::reweighting::ReducedBias;

/// Which geometries of a set a [`Distribution`] covers: its samples.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Samples {
    /// Every geometry of the set.
    All,
    /// The listed geometries, by strictly increasing index; the others are
    /// not samples.
    Subset(Vec<usize>),
}

/// The data a [`Distribution`] stores for one state: an optional energy and
/// owner properties.
///
/// Which samples a state groups, and which of them represents it, belong to
/// the distribution's [`Partition`]; [`StateView`](super::StateView) reads
/// them together with the state's probability.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StateData {
    energy: Option<Quantity<f64>>,
    properties: OwnerProperties,
}

impl StateData {
    /// A state without an energy or properties.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns this state with a finite energy.
    pub fn with_energy(mut self, energy: Quantity<f64>) -> Result<Self, StructureSetError> {
        self.set_energy(Some(energy))?;
        Ok(self)
    }

    /// The state's energy in [`CANONICAL_ENERGY_UNIT`], if known.
    ///
    /// Energies of one distribution share one scale, such as the minimized
    /// energies of conformers or computed free energies at one level of
    /// theory, never a snapshot's instantaneous energy; only differences
    /// between states are meaningful. Record the scale in the set
    /// properties.
    pub const fn energy(&self) -> Option<Quantity<f64>> {
        self.energy
    }

    /// Sets or clears a finite energy, converted to the canonical energy
    /// unit. No weight changes.
    pub fn set_energy(&mut self, energy: Option<Quantity<f64>>) -> Result<(), StructureSetError> {
        self.energy = energy
            .map(|energy| {
                let energy = energy.into_unit(CANONICAL_ENERGY_UNIT)?;
                if !energy.value().is_finite() {
                    return Err(StructureSetError::NonFiniteEnergy);
                }
                Ok(energy)
            })
            .transpose()?;
        Ok(())
    }

    /// State-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    /// Mutable state-level annotations.
    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }
}

/// How a [`Distribution`] weighs its samples
/// ([`Distribution::weight_model`]).
///
/// Log weights and reduced biases are natural logarithms, one finite value
/// per sample, in sample order. Only their differences are meaningful.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum WeightModel<'a> {
    /// Every sample weighs the same.
    Equal,
    /// Log weights the caller gave.
    LogWeights(&'a [f64]),
    /// Log weights assigned as Boltzmann populations of the state energies
    /// at the target ([`DistributionMut::assign_boltzmann_populations`]).
    /// Later edits to samples, states, or energies keep the tag; setting
    /// weights, or a target at another temperature, replaces it.
    Boltzmann(&'a [f64]),
    /// Weights derived from a structure set's sampling record against the
    /// target ([`Distribution::reweighted`]): each sample's reduced bias
    /// `b_i`, its exact log importance weight, with the estimator that sums
    /// groups of them into populations.
    Reweighted {
        /// The estimator of group populations.
        estimator: Estimator,
        /// The reduced bias of each sample.
        reduced: &'a [f64],
    },
}

/// The largest magnitude of a reduced bias. Within it, every group
/// population is finite under every [`Estimator`]: the third cumulant of
/// such values stays below `f64::MAX`.
pub(crate) const REDUCED_BIAS_LIMIT: f64 = 1e100;

/// The stored weights of a [`Distribution`]; [`WeightModel`] reads them.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Weights {
    Equal,
    Logarithmic(Vec<f64>),
    /// Logarithmic weights assigned as Boltzmann populations; the
    /// distribution has a target.
    Boltzmann(Vec<f64>),
    /// The reduced bias of each sample, each at most
    /// [`REDUCED_BIAS_LIMIT`] in magnitude, summed into group populations by
    /// `estimator`; the distribution has a target.
    Reweighted {
        estimator: Estimator,
        reduced: Vec<f64>,
    },
}

impl Weights {
    /// The exact log weights in sample order: given, Boltzmann, or reduced
    /// biases; `None` when every sample weighs the same.
    fn exact(&self) -> Option<&[f64]> {
        match self {
            Self::Equal => None,
            Self::Logarithmic(values)
            | Self::Boltzmann(values)
            | Self::Reweighted {
                reduced: values, ..
            } => Some(values),
        }
    }

    /// The exact log weights in sample order for removing or reordering
    /// samples, or `None` when every sample weighs the same.
    fn exact_mut(&mut self) -> Option<&mut Vec<f64>> {
        match self {
            Self::Equal => None,
            Self::Logarithmic(values)
            | Self::Boltzmann(values)
            | Self::Reweighted {
                reduced: values, ..
            } => Some(values),
        }
    }

    /// The estimator of group populations: exponential unless reweighted.
    fn estimator(&self) -> Estimator {
        match self {
            Self::Reweighted { estimator, .. } => *estimator,
            _ => Estimator::Exponential,
        }
    }

    /// The log weights of `count` samples, with equal weights made `0.0`
    /// each; the tag is kept. Reweighted weights give their reduced biases,
    /// which only a common shift keeps meaningful.
    fn materialized(&mut self, count: usize) -> &mut Vec<f64> {
        if matches!(self, Self::Equal) {
            *self = Self::Logarithmic(vec![0.0; count]);
        }
        match self {
            Self::Equal => unreachable!("equal weights were just replaced"),
            Self::Logarithmic(values)
            | Self::Boltzmann(values)
            | Self::Reweighted {
                reduced: values, ..
            } => values,
        }
    }

    /// The same model over the samples at `positions`, in that order.
    fn selected(&self, positions: impl Iterator<Item = usize>) -> Self {
        let pick = |values: &[f64]| positions.map(|position| values[position]).collect();
        match self {
            Self::Equal => Self::Equal,
            Self::Logarithmic(values) => Self::Logarithmic(pick(values)),
            Self::Boltzmann(values) => Self::Boltzmann(pick(values)),
            Self::Reweighted { estimator, reduced } => Self::Reweighted {
                estimator: *estimator,
                reduced: pick(reduced),
            },
        }
    }
}

/// The state populations of a distribution, from its weights and partition.
///
/// Estimators are shift-equivariant: the estimate over `ℓ_i − c` is the
/// estimate over `ℓ_i` less `c`. Every value is taken relative to the
/// largest exact log weight, so normalizations keep their digits under a
/// large common offset.
#[derive(Debug, Clone)]
struct Populations {
    /// The largest exact log weight, `0.0` for equal weights or no samples.
    shift: f64,
    /// Each state's estimate over its shifted exact log weights.
    estimates: Vec<f64>,
    /// Each state's exact `ln Σ exp` over its shifted exact log weights.
    sums: Vec<f64>,
    /// `ln Σ exp` over the estimates.
    total: f64,
}

/// [`Populations`] computed on first use. Every change to the weights or
/// the partition clears it. It carries no information of its own, so it
/// never makes two distributions differ.
#[derive(Clone, Default)]
struct PopulationCache(OnceLock<Populations>);

impl PartialEq for PopulationCache {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

impl std::fmt::Debug for PopulationCache {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("PopulationCache")
    }
}

/// Geometries `samples` of a set of `len`, as [`Samples::All`] when they
/// are every geometry.
fn normalized(samples: Vec<usize>, len: usize) -> Samples {
    if samples.len() == len {
        Samples::All
    } else {
        Samples::Subset(samples)
    }
}

/// A distribution over some or all geometries of a structure set: weighted
/// samples grouped into states.
///
/// It has two levels, as a distribution of conformations does.
///
/// - A *sample* is a geometry the distribution covers ([`Samples`]), with a
///   *weight*: the share of probability it carries. Only ratios of weights
///   are meaningful, so they are stored as natural logarithms, which any
///   finite value is: weights from bias reweighting or importance sampling
///   span hundreds of orders of magnitude, beyond what a linear weight can
///   hold. The weights are equal, given, Boltzmann populations of the
///   states' energies, or reweighted from a biased run ([`WeightModel`]).
/// - A *state* is a distinct conformation: a group of samples, one of them
///   its representative. Its *population* is derived and never stored. Its
///   data ([`StateData`]) may carry
///   an *energy* on one scale for the whole distribution, such as a
///   conformer's minimized energy or a computed free energy, never a
///   snapshot's instantaneous energy, which stays a property of its
///   geometry. [`StateView`](super::StateView) reads a state whole.
///
/// Every sample is in exactly one state, and every state has at least one
/// sample, its representative ([`Partition`]). The distribution may record
/// its *target*, the [`ThermodynamicState`] its populations describe.
///
/// The population of any group of samples, a state, a histogram bin, or a
/// map cell, is the [`Estimator`] over the exact log weights of its samples
/// ([`Self::estimator`]): the exponential sum `ln Σ exp(ℓ_i)`, unless the
/// weights are reweighted with another. A sample's share within its group
/// is always exact, `exp(ℓ_i) / Σ exp(ℓ_j)` over the group. A sample's
/// probability is its state's probability times its share within the
/// state, so the probabilities of a state's samples sum to the state's.
/// State probabilities are the state populations normalized to sum to one
/// ([`Self::log_population`] gives any group's on that scale), and
/// [`Self::group_log_populations`] normalizes bins among themselves. Under
/// a cumulant estimator the populations of a group's parts need not add up
/// to the group's, which is inherent to the method.
///
/// Two index spaces meet here. Methods that take or return a `geometry` use
/// geometry indices of the set. [`Self::weight_model`], [`Self::partition`],
/// and [`Self::probabilities`] are in sample order, by sample position: the
/// `i`th sample is the `i`th geometry [`Samples`] lists.
/// [`Self::sample_probabilities`] joins the two.
///
/// A conformer ensemble is one state per conformer; an NMR bundle or the
/// frames of a trajectory are one state of many samples, and clustering
/// makes several.
///
/// ```
/// use kekule::structure::{Distribution, Partition, Samples};
/// // Frames 1 and 2 of a three-frame trajectory, as one state, the first
/// // weighing three times the second.
/// let distribution = Distribution::new(Samples::Subset(vec![1, 2]), Partition::one_state(2))?
///     .with_log_weights(vec![3.0_f64.ln(), 0.0])?;
/// assert_eq!(distribution.log_weight(0), None);
/// let probabilities = distribution.probabilities()?;
/// assert!((probabilities[0] - 0.75).abs() < 1e-15);
/// assert!((probabilities[1] - 0.25).abs() < 1e-15);
/// assert_eq!(distribution.ranked_by_weight(), [1, 2]);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Distribution {
    samples: Samples,
    weights: Weights,
    partition: Partition,
    states: Vec<StateData>,
    target: Option<ThermodynamicState>,
    populations: PopulationCache,
}

impl Distribution {
    /// A distribution of `samples` grouped by `partition`, with equal
    /// weights, without state energies, state properties, or a target.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::UnsortedSamples`] when a subset's indices do not
    /// increase, and [`StructureSetError::PartitionMismatch`] when the
    /// partition covers another number of samples than the subset lists. A
    /// distribution of [`Samples::All`] is checked against its set when an
    /// [`Ensemble`](super::Ensemble) is built.
    pub fn new(samples: Samples, partition: Partition) -> Result<Self, StructureSetError> {
        if let Samples::Subset(indices) = &samples {
            if indices.windows(2).any(|pair| pair[0] >= pair[1]) {
                return Err(StructureSetError::UnsortedSamples);
            }
            if indices.len() != partition.sample_count() {
                return Err(StructureSetError::PartitionMismatch(LengthMismatch::new(
                    Countable::Sample,
                    indices.len(),
                    partition.sample_count(),
                )));
            }
        }
        Ok(Self {
            samples,
            weights: Weights::Equal,
            states: vec![StateData::new(); partition.state_count()],
            partition,
            target: None,
            populations: PopulationCache::default(),
        })
    }

    /// A distribution of `samples` of `set` grouped by `partition`, each
    /// weighted by its reduced bias from the set's sampling record, with the
    /// estimator `estimator` of group populations. Its target is the sampled
    /// state at `temperature`, or the sampled state itself for `None`.
    ///
    /// A sample's reduced bias is `b = β_s (U + ΣΔV) − β_t U`, its exact log
    /// importance weight, with `β = 1/RT` at the sampled temperature `s` and
    /// the target temperature `t`, its potential energy `U`, and the sum
    /// `ΣΔV` of its bias columns. At the sampled temperature
    /// `b = β_s ΣΔV`, and the potential energy is not read. Reweighting
    /// changes the temperature only: the target keeps the sampled pressure,
    /// and the `pV` term is neglected, which is negligible near ambient
    /// pressure. Without bias columns, samples at the sampled temperature
    /// weigh the same.
    ///
    /// ```
    /// use kekule::properties::{PropertyColumn, PropertyKey};
    /// use kekule::statistics::Estimator;
    /// use kekule::structure::{
    ///     Distribution, Geometry, Partition, Positions, Samples, Sampling, StructureSet,
    ///     ThermodynamicState,
    /// };
    /// use kekule::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE};
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let atoms = topology.atom_count();
    /// let geometries = (0..2).map(|_| Geometry::new(Positions::zeros(atoms)));
    /// let mut set = StructureSet::from_geometries(topology, geometries)?;
    /// let boost = PropertyKey::new("boost")?;
    /// set.columns_mut().insert(
    ///     boost.clone(),
    ///     PropertyColumn::Real {
    ///         unit: KILOJOULE_PER_MOLE,
    ///         values: vec![Some(0.0), Some(2.5)],
    ///     },
    /// )?;
    /// let state = ThermodynamicState::canonical(Quantity::new(300.0, KELVIN))?;
    /// set.set_sampling(Some(Sampling::new(state).with_bias([boost])?))?;
    /// let distribution = Distribution::reweighted(
    ///     &set,
    ///     Samples::All,
    ///     Partition::one_state(2),
    ///     None,
    ///     Estimator::Exponential,
    /// )?;
    /// // The boosted frame weighs exp(2.5 kJ/mol / RT) times the other.
    /// let ratio = distribution.log_weight(1).unwrap() - distribution.log_weight(0).unwrap();
    /// assert!((ratio - 2.5 / *state.thermal_energy().value()).abs() < 1e-12);
    /// assert_eq!(distribution.target(), Some(state));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// The errors of [`Self::new`] and of
    /// [`Ensemble::new`](super::Ensemble::new) against `set`;
    /// `Missing(Sampling)` for a set without a sampling record;
    /// `Missing(PotentialEnergy)` at another temperature when the record
    /// names no potential-energy column; `Missing(Steps)` when the record's
    /// first step is not 0 and the set has no steps; the errors of
    /// [`ThermodynamicState::with_temperature`] for `temperature`;
    /// [`StructureSetError::SamplingColumn`] for a named column that holds no
    /// real energies; and, for the first sample that cannot be one,
    /// [`StructureSetError::BeforeSampling`] when its step lies before the
    /// record's first step, [`StructureSetError::MissingColumnValue`] when a
    /// cell it reads is missing, and [`StructureSetError::InvalidReducedBias`]
    /// when its reduced bias is not finite or exceeds `1e100` in magnitude.
    pub fn reweighted(
        set: &(impl StructureSetRead + ?Sized),
        samples: Samples,
        partition: Partition,
        temperature: Option<Quantity<f64>>,
        estimator: Estimator,
    ) -> Result<Self, StructureSetError> {
        let bias = ReducedBias::new(set.sampling(), set.columns(), set.order(), temperature)?;
        Self::reweighted_by(&bias, set.len(), samples, partition, estimator)
    }

    /// A distribution of `samples` of a set of `len` geometries grouped by
    /// `partition`, weighted by the reduced biases `bias` gives, as
    /// [`Self::reweighted`] makes it.
    pub(crate) fn reweighted_by(
        bias: &ReducedBias<'_>,
        len: usize,
        samples: Samples,
        partition: Partition,
        estimator: Estimator,
    ) -> Result<Self, StructureSetError> {
        let mut distribution = Self::new(samples, partition)?;
        distribution.bind(len)?;
        let reduced = (0..distribution.sample_count())
            .map(|position| bias.reduced(distribution.sample(position)))
            .collect::<Result<Vec<_>, _>>()?;
        distribution.set_weights(Weights::Reweighted { estimator, reduced })?;
        distribution.target = Some(bias.target());
        Ok(distribution)
    }

    /// A distribution of `samples` of a set of `len` geometries as one state,
    /// weighted by their reduced biases `reduced` (from `bias`, checked)
    /// toward its target.
    pub(crate) fn reweighted_one_state(
        bias: &ReducedBias<'_>,
        len: usize,
        samples: Vec<usize>,
        reduced: Vec<f64>,
        estimator: Estimator,
    ) -> Result<Self, StructureSetError> {
        let partition = Partition::one_state(samples.len());
        let mut distribution = Self::new(Samples::Subset(samples), partition)?;
        distribution.bind(len)?;
        distribution.set_weights(Weights::Reweighted { estimator, reduced })?;
        distribution.target = Some(bias.target());
        Ok(distribution)
    }

    /// Returns this distribution with one finite natural-logarithm weight
    /// per sample, in sample order; see [`DistributionMut::set_log_weights`].
    pub fn with_log_weights(mut self, log_weights: Vec<f64>) -> Result<Self, StructureSetError> {
        self.set_weights(Weights::Logarithmic(log_weights))?;
        Ok(self)
    }

    /// Returns this distribution with stored weights, checked as
    /// [`Self::with_log_weights`] checks them.
    pub(crate) fn with_weights(mut self, weights: Weights) -> Result<Self, StructureSetError> {
        self.set_weights(weights)?;
        Ok(self)
    }

    /// Returns this distribution with one [`StateData`] (energy and
    /// properties) per state of its partition, in state order.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::StateCountMismatch`] for another number of states.
    pub fn with_states(mut self, states: Vec<StateData>) -> Result<Self, StructureSetError> {
        if states.len() != self.states.len() {
            return Err(StructureSetError::StateCountMismatch(LengthMismatch::new(
                Countable::State,
                self.states.len(),
                states.len(),
            )));
        }
        self.states = states;
        Ok(self)
    }

    /// Returns this distribution with a target; see
    /// [`DistributionMut::set_target`].
    ///
    /// # Errors
    ///
    /// [`StructureSetError::ReweightedTarget`] for reweighted samples.
    pub fn with_target(mut self, target: ThermodynamicState) -> Result<Self, StructureSetError> {
        self.set_target(Some(target))?;
        Ok(self)
    }

    /// Which geometries are samples.
    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    /// Number of samples.
    pub fn sample_count(&self) -> usize {
        self.partition.sample_count()
    }

    /// The position of geometry `geometry` among the samples, which indexes
    /// the log weights and the partition, or `None` when it is not a sample.
    pub fn sample_position(&self, geometry: usize) -> Option<usize> {
        match &self.samples {
            Samples::All => (geometry < self.sample_count()).then_some(geometry),
            Samples::Subset(indices) => indices.binary_search(&geometry).ok(),
        }
    }

    /// The geometry index of the sample at `position`.
    pub(crate) fn sample(&self, position: usize) -> usize {
        match &self.samples {
            Samples::All => position,
            Samples::Subset(indices) => indices[position],
        }
    }

    /// How the samples are weighed, with the log weights in sample order.
    pub fn weight_model(&self) -> WeightModel<'_> {
        match &self.weights {
            Weights::Equal => WeightModel::Equal,
            Weights::Logarithmic(log_weights) => WeightModel::LogWeights(log_weights),
            Weights::Boltzmann(log_weights) => WeightModel::Boltzmann(log_weights),
            Weights::Reweighted { estimator, reduced } => WeightModel::Reweighted {
                estimator: *estimator,
                reduced,
            },
        }
    }

    /// The estimator of group populations: the reweighted weights' own, and
    /// [`Estimator::Exponential`], the exact sum, for every other model.
    pub fn estimator(&self) -> Estimator {
        self.weights.estimator()
    }

    /// The natural-logarithm weight of geometry `geometry`, which its
    /// probability follows: its probability times a factor shared by every
    /// sample. It is its own log weight, `0.0` when every sample weighs the
    /// same, or `None` when it is not a sample. Under reweighting with an
    /// estimator other than the exponential sum, it is the natural logarithm
    /// of its probability: its state's log probability plus its exact log
    /// share within the state.
    pub fn log_weight(&self, geometry: usize) -> Option<f64> {
        self.sample_position(geometry)
            .map(|position| self.sample_log_weight(position))
    }

    /// The grouping of the samples into states, by sample position.
    pub fn partition(&self) -> &Partition {
        &self.partition
    }

    /// Number of states.
    pub fn state_count(&self) -> usize {
        self.states.len()
    }

    /// The data of every state, its energy and properties, in state
    /// order, as [`Self::with_states`] takes them.
    pub fn states(&self) -> &[StateData] {
        &self.states
    }

    /// The state of geometry `geometry`, or `None` when it is not a sample.
    pub fn state_of(&self, geometry: usize) -> Option<usize> {
        self.sample_position(geometry)
            .map(|position| self.partition.assignments()[position])
    }

    /// The thermodynamic state the populations describe, if recorded.
    /// Changing weights does not clear it; the caller states what its
    /// populations describe.
    pub const fn target(&self) -> Option<ThermodynamicState> {
        self.target
    }

    /// The target's temperature, in
    /// [`CANONICAL_TEMPERATURE_UNIT`](crate::units::CANONICAL_TEMPERATURE_UNIT),
    /// if a target is recorded.
    pub fn temperature(&self) -> Option<Quantity<f64>> {
        self.target.map(|target| target.temperature())
    }

    /// This distribution in a set of `len` geometries whose first ones are
    /// those it describes, such as a file that grew after the distribution
    /// was saved. A distribution of every geometry becomes one of the
    /// geometries it covered, [`Samples::Subset`] of `0..n` for `n` samples,
    /// so the geometries after them are not samples; a subset keeps its
    /// samples. A subset listing all `len` geometries becomes
    /// [`Samples::All`]. The weights, the partition, the states, and the
    /// target are kept bit for bit. A subset does not record how many
    /// geometries it was made for, so a `len` shorter than that is not
    /// detected.
    ///
    /// ```
    /// use kekule::structure::{Distribution, Partition, Samples};
    ///
    /// let saved = Distribution::new(Samples::All, Partition::one_state(3))?
    ///     .with_log_weights(vec![0.0, 1.0, 2.0])?;
    /// let grown = saved.clone().extended_to(5)?;
    /// assert_eq!(grown.samples(), &Samples::Subset(vec![0, 1, 2]));
    /// assert_eq!(grown.weight_model(), saved.weight_model());
    /// assert_eq!(grown.log_weight(4), None);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`StructureSetError::CoverageMismatch`] when a distribution of every
    /// geometry covers more than `len`, and
    /// [`StructureSetError::IndexOutOfRange`] when a subset lists a geometry
    /// past them.
    pub fn extended_to(mut self, len: usize) -> Result<Self, StructureSetError> {
        self.extend(len)?;
        Ok(self)
    }

    /// Re-expresses this distribution in place for a set of `len`
    /// geometries, as [`Self::extended_to`] does; failure changes nothing.
    pub(crate) fn extend(&mut self, len: usize) -> Result<(), StructureSetError> {
        if self.samples == Samples::All {
            let covered = self.sample_count();
            if covered > len {
                return Err(StructureSetError::CoverageMismatch {
                    covered,
                    geometries: len,
                });
            }
            self.samples = normalized((0..covered).collect(), len);
        }
        self.bind(len)
    }

    /// The probability of each sample, in sample order: its state's
    /// probability times its exact share within the state, which for every
    /// model but reweighting with a cumulant estimator is its weight divided
    /// by the total weight. A probability below the smallest positive `f64`
    /// is `0.0`.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptyDistribution`] without samples.
    pub fn probabilities(&self) -> Result<Vec<f64>, StructureSetError> {
        if self.sample_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        Ok(normalized_exponentials(
            self.sample_log_weights().into_iter(),
        ))
    }

    /// Each sample's geometry index with its probability, in sample order,
    /// as [`Self::probabilities`] gives them: what weighted statistics over
    /// the distribution observe.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptyDistribution`] without samples.
    pub fn sample_probabilities(&self) -> Result<Vec<(usize, f64)>, StructureSetError> {
        Ok(self
            .probabilities()?
            .into_iter()
            .enumerate()
            .map(|(position, probability)| (self.sample(position), probability))
            .collect())
    }

    /// The probability of each state, in state order: its population
    /// divided by the sum of the state populations.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptyDistribution`] without samples.
    pub fn state_probabilities(&self) -> Result<Vec<f64>, StructureSetError> {
        if self.sample_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        Ok(normalized_exponentials(
            self.populations().estimates.iter().copied(),
        ))
    }

    /// The natural logarithm of the probability of the samples at
    /// `positions`, sample positions as [`Self::partition`] indexes them:
    /// the estimator over their exact log weights, relative to the sum of
    /// the state populations. For a state's samples it is the logarithm of
    /// the state's probability. A position listed twice counts twice.
    ///
    /// `None` for no positions or a position that is not a sample's.
    pub fn log_population(&self, positions: &[usize]) -> Option<f64> {
        let populations = self.populations();
        let estimate = self.group_estimate(positions, populations.shift)?;
        Some(estimate - populations.total)
    }

    /// The natural logarithm of each group's probability, in group order,
    /// for groups of samples such as the bins of a histogram: the estimator
    /// over each group's exact log weights, normalized so that the
    /// probabilities of the groups sum to one. `group_of` gives each
    /// sample's group, by sample position, or `None` for a sample in no
    /// group; a group without samples has no population.
    ///
    /// Groups that are the states give the state probabilities.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::GroupingMismatch`] when `group_of` covers another
    /// number of samples than there are, and
    /// [`StructureSetError::GroupOutOfRange`] for the first group of
    /// `group_of` that is not below `groups`.
    pub fn group_log_populations(
        &self,
        group_of: &[Option<u32>],
        groups: usize,
    ) -> Result<Vec<Option<f64>>, StructureSetError> {
        if group_of.len() != self.sample_count() {
            return Err(StructureSetError::GroupingMismatch(LengthMismatch::new(
                Countable::Sample,
                self.sample_count(),
                group_of.len(),
            )));
        }
        let grouping = Grouping::new(group_of, groups).map_err(|group| {
            StructureSetError::GroupOutOfRange(IndexOutOfRange::new(
                Countable::Group,
                group,
                groups,
            ))
        })?;
        Ok(self.grouped_log_populations(&grouping))
    }

    /// The natural logarithm of each group's probability, as
    /// [`Self::group_log_populations`] gives it, for a grouping of this
    /// distribution's samples.
    pub(crate) fn grouped_log_populations(&self, grouping: &Grouping) -> Vec<Option<f64>> {
        let shift = self.populations().shift;
        let estimates = (0..grouping.len())
            .map(|group| self.group_estimate(grouping.members(group), shift))
            .collect::<Vec<_>>();
        let total = log_sum_exp(estimates.iter().flatten().copied());
        estimates
            .into_iter()
            .map(|estimate| estimate.map(|estimate| estimate - total))
            .collect()
    }

    /// Kish's effective sample size, `(Σw)² / Σw²`, over the exact weights
    /// `w = exp(ℓ)`: the number of equally weighted samples that would give
    /// the same statistical precision. For reweighted samples the exact
    /// weights are the exponentials of their reduced biases, whatever the
    /// estimator, so it measures how far the reweighting degenerates.
    ///
    /// It is the sample count for equal weights and approaches one when a
    /// single sample dominates.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptyDistribution`] without samples.
    pub fn effective_sample_size(&self) -> Result<f64, StructureSetError> {
        if self.sample_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        Ok(match self.weights.exact() {
            None => self.sample_count() as f64,
            Some(log_weights) => kish(log_weights.iter().copied()),
        })
    }

    /// Kish's effective sample size over the exact weights of the samples at
    /// `positions`, as [`Self::effective_sample_size`] gives it for all of
    /// them: the effective count of a group such as a histogram bin. A
    /// position listed twice counts twice.
    ///
    /// `None` for no positions or a position that is not a sample's.
    pub fn group_effective_sample_size(&self, positions: &[usize]) -> Option<f64> {
        if positions.is_empty() || positions.iter().any(|&p| p >= self.sample_count()) {
            return None;
        }
        Some(match self.weights.exact() {
            None => positions.len() as f64,
            Some(log_weights) => kish(positions.iter().map(|&position| log_weights[position])),
        })
    }

    /// Sample geometry indices, not sample positions, by descending
    /// probability; equal probabilities keep geometry order.
    pub fn ranked_by_weight(&self) -> Vec<usize> {
        let log_weights = self.sample_log_weights();
        let mut positions = (0..self.sample_count()).collect::<Vec<_>>();
        positions.sort_by(|&a, &b| log_weights[b].total_cmp(&log_weights[a]));
        positions
            .into_iter()
            .map(|position| self.sample(position))
            .collect()
    }

    /// State indices by descending population; equal populations keep
    /// state order.
    pub fn ranked_by_population(&self) -> Vec<usize> {
        let populations = &self.populations().estimates;
        let mut indices = (0..populations.len()).collect::<Vec<_>>();
        indices.sort_by(|&a, &b| populations[b].total_cmp(&populations[a]));
        indices
    }

    /// The fewest states, by descending population, whose populations
    /// together reach `fraction` of the total, in
    /// [`Self::ranked_by_population`] order.
    ///
    /// `fraction` must lie in `(0, 1]`; `1` selects every state, including
    /// those too light to change a floating-point total. Pass the result to
    /// [`Ensemble::select_states`](super::Ensemble::select_states) to keep,
    /// for example, 95% of a Boltzmann population.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::InvalidFraction`] for a fraction outside
    /// `(0, 1]` and [`StructureSetError::EmptyDistribution`] without samples.
    pub fn covering_population(&self, fraction: f64) -> Result<Vec<usize>, StructureSetError> {
        if !(fraction > 0.0 && fraction <= 1.0) {
            return Err(StructureSetError::InvalidFraction);
        }
        if self.sample_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        let ranked = self.ranked_by_population();
        if fraction == 1.0 {
            return Ok(ranked);
        }
        let shares = self.state_probabilities()?;
        // Sum in ranked order so the running total ends exactly at the total.
        let total = ranked.iter().map(|&index| shares[index]).sum::<f64>();
        let target = fraction * total;
        let mut covered = 0.0;
        let count = ranked
            .iter()
            .position(|&index| {
                covered += shares[index];
                covered >= target
            })
            .map_or(ranked.len(), |position| position + 1);
        Ok(ranked[..count].to_vec())
    }

    /// State indices by ascending energy; equal energies keep state order.
    ///
    /// # Errors
    ///
    /// `Missing(StateEnergies)` if any state has no energy.
    pub fn ranked_by_energy(&self) -> Result<Vec<usize>, StructureSetError> {
        let energies = self.state_energies()?;
        let mut indices = (0..energies.len()).collect::<Vec<_>>();
        indices.sort_by(|&a, &b| energies[a].total_cmp(&energies[b]));
        Ok(indices)
    }

    /// Indices, in state order, of the states whose energy lies at most
    /// `window` above the lowest state energy.
    ///
    /// `window` is a finite, nonnegative energy, such as the 6 kcal/mol
    /// window common in conformer searches.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::InvalidEnergyWindow`] and
    /// [`StructureSetError::Unit`] for a window that is not a finite,
    /// nonnegative energy, and `Missing(StateEnergies)` if any state has no
    /// energy.
    pub fn within_energy_window(
        &self,
        window: Quantity<f64>,
    ) -> Result<Vec<usize>, StructureSetError> {
        let window = window.value_in(CANONICAL_ENERGY_UNIT)?;
        if !(window.is_finite() && window >= 0.0) {
            return Err(StructureSetError::InvalidEnergyWindow);
        }
        let energies = self.state_energies()?;
        let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
        Ok(energies
            .iter()
            .enumerate()
            .filter(|(_, energy)| *energy - minimum <= window)
            .map(|(index, _)| index)
            .collect())
    }

    /// Records or clears the target. Boltzmann populations lose their tag
    /// unless the temperature stays the same; reweighted weights refuse.
    fn set_target(&mut self, target: Option<ThermodynamicState>) -> Result<(), StructureSetError> {
        if matches!(self.weights, Weights::Reweighted { .. }) {
            return Err(StructureSetError::ReweightedTarget);
        }
        let same_temperature = match (self.target, target) {
            (Some(old), Some(new)) => old.temperature().value() == new.temperature().value(),
            _ => false,
        };
        if !same_temperature {
            self.untag();
        }
        self.target = target;
        Ok(())
    }

    /// Makes Boltzmann populations plain log weights, and reweighted
    /// weights the log weights [`Self::log_weight`] gives, so every
    /// probability is kept.
    fn untag(&mut self) {
        self.populations = PopulationCache::default();
        match &mut self.weights {
            Weights::Boltzmann(log_weights) => {
                self.weights = Weights::Logarithmic(std::mem::take(log_weights));
            }
            Weights::Reweighted { .. } => {
                self.weights = Weights::Logarithmic(self.sample_log_weights());
            }
            Weights::Equal | Weights::Logarithmic(_) => {}
        }
    }

    /// Replaces the weights after checking that their values, if any, are
    /// one finite value per sample, and reduced biases within
    /// [`REDUCED_BIAS_LIMIT`].
    pub(crate) fn set_weights(&mut self, weights: Weights) -> Result<(), StructureSetError> {
        if let Some(values) = weights.exact() {
            if values.len() != self.sample_count() {
                return Err(StructureSetError::LogWeightCountMismatch(
                    LengthMismatch::new(Countable::LogWeight, self.sample_count(), values.len()),
                ));
            }
            let reweighted = matches!(weights, Weights::Reweighted { .. });
            if let Some(position) = values.iter().position(|value| {
                !value.is_finite() || (reweighted && value.abs() > REDUCED_BIAS_LIMIT)
            }) {
                let geometry = self.sample(position);
                return Err(if reweighted {
                    StructureSetError::InvalidReducedBias { geometry }
                } else {
                    StructureSetError::NonFiniteLogWeight { geometry }
                });
            }
        }
        self.weights = weights;
        self.populations = PopulationCache::default();
        Ok(())
    }

    /// Every state energy in [`CANONICAL_ENERGY_UNIT`].
    pub(crate) fn state_energies(&self) -> Result<Vec<f64>, StructureSetError> {
        self.states
            .iter()
            .map(|state| {
                state
                    .energy()
                    .map(Quantity::into_value)
                    .ok_or(StructureSetError::Missing(Requirement::StateEnergies))
            })
            .collect()
    }

    /// The exact log weight of the sample at `position`: its own, its
    /// reduced bias, or `0.0` when every sample weighs the same.
    fn exact_log_weight(&self, position: usize) -> f64 {
        self.weights.exact().map_or(0.0, |values| values[position])
    }

    /// Every sample's exact log weight in sample order: given, Boltzmann,
    /// or reduced biases, and `0.0` each when every sample weighs the same.
    /// Unlike [`Self::sample_log_weights`], it does not depend on the
    /// partition.
    pub(crate) fn exact_log_weights(&self) -> Vec<f64> {
        self.weights
            .exact()
            .map_or_else(|| vec![0.0; self.sample_count()], <[f64]>::to_vec)
    }

    /// The estimator's `ln Σ exp(ℓ_i − shift)` over the exact log weights of
    /// the samples at `positions`, or `None` for none or a position past the
    /// last sample.
    fn group_estimate(&self, positions: &[usize], shift: f64) -> Option<f64> {
        if positions
            .iter()
            .any(|&position| position >= self.sample_count())
        {
            return None;
        }
        let values = positions
            .iter()
            .map(|&position| self.exact_log_weight(position) - shift)
            .collect::<Vec<_>>();
        self.weights.estimator().log_sum(&values)
    }

    /// The state populations, computed once per change of the weights or
    /// the partition.
    fn populations(&self) -> &Populations {
        self.populations.0.get_or_init(|| {
            let shift = match self.weights.exact() {
                Some(values) if !values.is_empty() => {
                    values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                }
                _ => 0.0,
            };
            let estimator = self.weights.estimator();
            let (estimates, sums) = self
                .partition
                .samples_by_state()
                .iter()
                .map(|positions| {
                    let values = positions
                        .iter()
                        .map(|&position| self.exact_log_weight(position) - shift)
                        .collect::<Vec<_>>();
                    let sum = log_sum_exp(values.iter().copied());
                    let estimate = match estimator {
                        Estimator::Exponential => sum,
                        estimator => estimator.log_sum(&values).expect(
                            "a state has samples, and bounded reduced biases give finite estimates",
                        ),
                    };
                    (estimate, sum)
                })
                .unzip::<_, _, Vec<_>, Vec<_>>();
            let total = log_sum_exp(estimates.iter().copied());
            Populations {
                shift,
                estimates,
                sums,
                total,
            }
        })
    }

    /// The log weight of the sample at `position`, as [`Self::log_weight`]
    /// gives it.
    fn sample_log_weight(&self, position: usize) -> f64 {
        let exact = self.exact_log_weight(position);
        if !self.reweighted_by_moments() {
            return exact;
        }
        let state = self.partition.assignments()[position];
        let populations = self.populations();
        (exact - populations.shift - populations.sums[state])
            + (populations.estimates[state] - populations.total)
    }

    /// Every sample's log weight in sample order, as [`Self::log_weight`]
    /// gives them: `0.0` each when every sample weighs the same.
    pub(crate) fn sample_log_weights(&self) -> Vec<f64> {
        let mut log_weights = self.exact_log_weights();
        if self.reweighted_by_moments() {
            let populations = self.populations();
            for (log_weight, &state) in log_weights.iter_mut().zip(self.partition.assignments()) {
                *log_weight = (*log_weight - populations.shift - populations.sums[state])
                    + (populations.estimates[state] - populations.total);
            }
        }
        log_weights
    }

    /// Whether the weights are reweighted with an estimator other than the
    /// exact sum, so that log weights differ from reduced biases.
    fn reweighted_by_moments(&self) -> bool {
        matches!(
            self.weights,
            Weights::Reweighted { estimator, .. } if estimator != Estimator::Exponential
        )
    }

    /// The natural logarithm of each state's population, in state order, on
    /// the scale of the exact log weights.
    fn state_log_populations(&self) -> Vec<f64> {
        self.partition
            .samples_by_state()
            .iter()
            .map(|positions| {
                self.group_estimate(positions, 0.0)
                    .expect("a state has samples, and bounded reduced biases give finite estimates")
            })
            .collect()
    }

    /// The natural logarithm of one state's population, which must exist,
    /// on the scale of the exact log weights.
    pub(crate) fn state_log_population(&self, state: usize) -> f64 {
        self.group_estimate(&self.partition.samples(state).collect::<Vec<_>>(), 0.0)
            .expect("a state has samples, and bounded reduced biases give finite estimates")
    }

    /// The data of one state, which must exist.
    pub(crate) fn state(&self, state: usize) -> &StateData {
        &self.states[state]
    }

    /// Checks this distribution against a set of `len` geometries, and lists
    /// every geometry as [`Samples::All`].
    pub(crate) fn bind(&mut self, len: usize) -> Result<(), StructureSetError> {
        match &self.samples {
            Samples::All => {
                if self.sample_count() != len {
                    return Err(StructureSetError::CoverageMismatch {
                        covered: self.sample_count(),
                        geometries: len,
                    });
                }
            }
            Samples::Subset(indices) => {
                if let Some(&last) = indices.last() {
                    if last >= len {
                        return Err(StructureSetError::IndexOutOfRange(IndexOutOfRange::new(
                            Countable::Geometry,
                            last,
                            len,
                        )));
                    }
                }
                if indices.len() == len {
                    self.samples = Samples::All;
                }
            }
        }
        Ok(())
    }

    /// Checks a new sample's log weight before [`Self::push_sample`] or
    /// [`Self::push_state`].
    pub(crate) fn check_new_log_weight(
        &self,
        log_weight: Option<f64>,
        geometry: usize,
    ) -> Result<(), StructureSetError> {
        match (log_weight, &self.weights) {
            (_, Weights::Reweighted { .. }) => {
                Err(StructureSetError::ReweightedSample { geometry })
            }
            (None, Weights::Logarithmic(_) | Weights::Boltzmann(_)) => {
                Err(StructureSetError::LogWeightRequired)
            }
            (Some(log_weight), _) if !log_weight.is_finite() => {
                Err(StructureSetError::NonFiniteLogWeight { geometry })
            }
            _ => Ok(()),
        }
    }

    /// Checks that `state` exists.
    pub(crate) fn check_state(&self, state: usize) -> Result<(), StructureSetError> {
        if state >= self.states.len() {
            return Err(StructureSetError::StateOutOfRange(IndexOutOfRange::new(
                Countable::State,
                state,
                self.states.len(),
            )));
        }
        Ok(())
    }

    /// Adds geometry `geometry`, the new last geometry of the set, to an
    /// existing state, after [`Self::check_state`] and
    /// [`Self::check_new_log_weight`].
    pub(crate) fn push_sample(&mut self, geometry: usize, state: usize, log_weight: Option<f64>) {
        self.push_log_weight(log_weight);
        self.push_position(geometry);
        let (mut assignments, representatives) = self.take_partition();
        assignments.push(state);
        self.partition = Partition::from_parts_unchecked(assignments, representatives);
    }

    /// Adds geometry `geometry`, the new last geometry of the set, as the
    /// representative of a new state, after [`Self::check_new_log_weight`];
    /// returns the state's index.
    pub(crate) fn push_state(
        &mut self,
        geometry: usize,
        state: StateData,
        log_weight: Option<f64>,
    ) -> usize {
        let index = self.states.len();
        let position = self.sample_count();
        self.push_log_weight(log_weight);
        self.push_position(geometry);
        let (mut assignments, mut representatives) = self.take_partition();
        assignments.push(index);
        representatives.push(position);
        self.partition = Partition::from_parts_unchecked(assignments, representatives);
        self.states.push(state);
        index
    }

    /// Adds a new last sample's log weight, after
    /// [`Self::check_new_log_weight`].
    fn push_log_weight(&mut self, log_weight: Option<f64>) {
        self.populations = PopulationCache::default();
        match log_weight {
            Some(log_weight) => self
                .weights
                .materialized(self.partition.sample_count())
                .push(log_weight),
            None => debug_assert!(
                self.weights == Weights::Equal,
                "new log weights are checked first"
            ),
        }
    }

    fn push_position(&mut self, geometry: usize) {
        if let Samples::Subset(indices) = &mut self.samples {
            indices.push(geometry);
        }
    }

    /// Checks that geometry `geometry` can be removed: it is not the
    /// representative of a state with other samples.
    pub(crate) fn check_remove(&self, geometry: usize) -> Result<(), StructureSetError> {
        let Some(position) = self.sample_position(geometry) else {
            return Ok(());
        };
        let state = self.partition.assignments()[position];
        let representative = self.partition.representatives()[state] == position;
        let alone = self.partition.samples(state).nth(1).is_none();
        if representative && !alone {
            return Err(StructureSetError::RepresentativeRemoval { geometry, state });
        }
        Ok(())
    }

    /// Removes geometry `geometry` from a set that keeps `len` geometries
    /// after [`Self::check_remove`], shifting later samples down. A sample
    /// alone in its state removes the state.
    pub(crate) fn remove(&mut self, geometry: usize, len: usize) {
        if let Some(position) = self.sample_position(geometry) {
            self.populations = PopulationCache::default();
            let state = self.partition.assignments()[position];
            let alone = self.partition.samples(state).nth(1).is_none();
            let (mut assignments, mut representatives) = self.take_partition();
            assignments.remove(position);
            if alone {
                representatives.remove(state);
                self.states.remove(state);
                for assigned in &mut assignments {
                    if *assigned > state {
                        *assigned -= 1;
                    }
                }
            }
            for represented in &mut representatives {
                if *represented > position {
                    *represented -= 1;
                }
            }
            self.partition = Partition::from_parts_unchecked(assignments, representatives);
            if let Some(log_weights) = self.weights.exact_mut() {
                log_weights.remove(position);
            }
        }
        if let Samples::Subset(indices) = &self.samples {
            let samples = indices
                .iter()
                .filter(|index| **index != geometry)
                .map(|index| if *index > geometry { index - 1 } else { *index })
                .collect();
            self.samples = normalized(samples, len);
        }
    }

    /// The distribution of a set of the geometries `indices` (checked), in
    /// order: the selected samples, with their weights. Every state with
    /// a selected sample becomes a state of the copy, in state order,
    /// represented by the first copy of its representative or, when that is
    /// not selected, by its first selected sample. A state whose samples are
    /// all selected keeps its energy and properties; one selected only in
    /// part has neither.
    pub(crate) fn select(&self, indices: &[usize]) -> Self {
        // (geometry index in the copy, sample position here)
        let picked = indices
            .iter()
            .enumerate()
            .filter_map(|(copy, geometry)| {
                self.sample_position(*geometry)
                    .map(|position| (copy, position))
            })
            .collect::<Vec<_>>();
        let old = self.partition.assignments();
        let mut selected = vec![0; self.states.len()];
        let mut chosen = vec![false; self.sample_count()];
        for (_, position) in &picked {
            if !std::mem::replace(&mut chosen[*position], true) {
                selected[old[*position]] += 1;
            }
        }
        let sizes = self.partition.sizes();
        let mut renumbered = vec![usize::MAX; self.states.len()];
        let mut states = Vec::new();
        for (state, count) in selected.iter().enumerate() {
            if *count > 0 {
                renumbered[state] = states.len();
                states.push(if *count == sizes[state] {
                    self.states[state].clone()
                } else {
                    StateData::new()
                });
            }
        }
        let mut assignments = Vec::with_capacity(picked.len());
        let mut representatives = vec![None; states.len()];
        let mut first = vec![None; states.len()];
        for (sample, (_, position)) in picked.iter().enumerate() {
            let state = renumbered[old[*position]];
            assignments.push(state);
            first[state].get_or_insert(sample);
            if *position == self.partition.representatives()[old[*position]] {
                representatives[state].get_or_insert(sample);
            }
        }
        let representatives = representatives
            .into_iter()
            .zip(first)
            .map(|(representative, first)| {
                representative
                    .or(first)
                    .expect("every renumbered state has a selected sample")
            })
            .collect();
        Self {
            samples: normalized(
                picked.iter().map(|(copy, _)| *copy).collect(),
                indices.len(),
            ),
            weights: self
                .weights
                .selected(picked.iter().map(|(_, position)| *position)),
            partition: Partition::from_parts_unchecked(assignments, representatives),
            states,
            target: self.target,
            populations: PopulationCache::default(),
        }
    }

    /// The samples of whole states `indices` (checked), in the requested
    /// order and each state's samples in geometry order, as their geometry
    /// indices, with the distribution of a set of just those geometries.
    pub(crate) fn select_states(&self, indices: &[usize]) -> (Vec<usize>, Self) {
        let mut geometries = Vec::new();
        let mut positions = Vec::new();
        let mut assignments = Vec::new();
        let mut representatives = Vec::with_capacity(indices.len());
        let mut states = Vec::with_capacity(indices.len());
        let by_state = self.partition.samples_by_state();
        for (state, old) in indices.iter().enumerate() {
            let representative = self.partition.representatives()[*old];
            for &position in &by_state[*old] {
                if position == representative {
                    representatives.push(geometries.len());
                }
                geometries.push(self.sample(position));
                positions.push(position);
                assignments.push(state);
            }
            states.push(self.states[*old].clone());
        }
        let distribution = Self {
            samples: Samples::All,
            weights: self.weights.selected(positions.into_iter()),
            partition: Partition::from_parts_unchecked(assignments, representatives),
            states,
            target: self.target,
            populations: PopulationCache::default(),
        };
        (geometries, distribution)
    }

    /// Each state's representative, as a geometry index, with the
    /// distribution of a set of just those geometries: one state each,
    /// weighted by its state's population.
    pub(crate) fn collapse(&self) -> (Vec<usize>, Self) {
        let representatives = self
            .partition
            .representatives()
            .iter()
            .map(|position| self.sample(*position))
            .collect::<Vec<_>>();
        let distribution = Self {
            samples: Samples::All,
            weights: match self.weights {
                Weights::Boltzmann(_) => Weights::Boltzmann(self.state_log_populations()),
                _ => Weights::Logarithmic(self.state_log_populations()),
            },
            partition: Partition::state_per_sample(representatives.len()),
            states: self.states.clone(),
            target: self.target,
            populations: PopulationCache::default(),
        };
        (representatives, distribution)
    }

    /// This distribution with the states' properties dropped, for a set
    /// projected onto an atom subset.
    pub(crate) fn without_state_properties(&self) -> Self {
        let mut distribution = self.clone();
        for state in &mut distribution.states {
            *state.properties_mut() = OwnerProperties::new();
        }
        distribution
    }

    /// Moves the partition out as its assignments and representatives, to
    /// be rebuilt and put back.
    fn take_partition(&mut self) -> (Vec<usize>, Vec<usize>) {
        std::mem::replace(&mut self.partition, Partition::state_per_sample(0)).into_parts()
    }
}

/// Mutable access to the distribution of an [`Ensemble`](super::Ensemble)
/// ([`Ensemble::distribution_mut`](super::Ensemble::distribution_mut)).
///
/// Reads use `Deref`. There is deliberately no `DerefMut`: the guard changes
/// the weights, the grouping into states, the states' energies and
/// properties, and the target, but never which geometries are samples, so
/// the distribution always fits its set.
///
/// ```compile_fail,E0594
/// use kekule::structure::{Distribution, Ensemble, Partition, Samples};
/// fn replace(ensemble: &mut Ensemble) {
///     let foreign = Distribution::new(Samples::Subset(vec![7]), Partition::one_state(1));
///     *ensemble.distribution_mut() = foreign.unwrap();
/// }
/// ```
#[derive(Debug)]
pub struct DistributionMut<'a> {
    distribution: &'a mut Distribution,
}

impl<'a> DistributionMut<'a> {
    pub(crate) fn new(distribution: &'a mut Distribution) -> Self {
        Self { distribution }
    }
}

impl std::ops::Deref for DistributionMut<'_> {
    type Target = Distribution;

    fn deref(&self) -> &Self::Target {
        self.distribution
    }
}

impl DistributionMut<'_> {
    /// Makes every sample weigh the same.
    pub fn set_equal_weights(&mut self) {
        self.distribution.weights = Weights::Equal;
        self.distribution.populations = PopulationCache::default();
    }

    /// Replaces every log weight, in sample order. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::LogWeightCountMismatch`] for another number of
    /// log weights than samples, and [`StructureSetError::NonFiniteLogWeight`]
    /// for the first one that is not finite.
    pub fn set_log_weights(&mut self, log_weights: Vec<f64>) -> Result<(), StructureSetError> {
        self.distribution
            .set_weights(Weights::Logarithmic(log_weights))
    }

    /// Sets the log weight of one sample, a finite value. With equal
    /// weights, the other samples keep log weight `0.0`; Boltzmann
    /// populations lose their tag, and reweighted samples take the log
    /// weights [`Distribution::log_weight`] gives them first.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::NotSample`] when `geometry` is not a sample, and
    /// [`StructureSetError::NonFiniteLogWeight`] for a log weight that is not
    /// finite.
    pub fn set_log_weight(
        &mut self,
        geometry: usize,
        log_weight: f64,
    ) -> Result<(), StructureSetError> {
        let distribution = &mut *self.distribution;
        let position = distribution
            .sample_position(geometry)
            .ok_or(StructureSetError::NotSample { geometry })?;
        if !log_weight.is_finite() {
            return Err(StructureSetError::NonFiniteLogWeight { geometry });
        }
        let count = distribution.sample_count();
        distribution.untag();
        distribution.weights.materialized(count)[position] = log_weight;
        distribution.populations = PopulationCache::default();
        Ok(())
    }

    /// Shifts every log weight so the weights sum to one, so each log weight
    /// is the natural logarithm of its sample's probability. Ratios, and so
    /// every probability and population share, are kept, and so are
    /// Boltzmann populations; equal weights become log weights. Reweighted
    /// samples shift their reduced biases together and keep their
    /// estimator, which no common shift changes.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptyDistribution`] without samples, and
    /// [`StructureSetError::InvalidReducedBias`] for the first reduced bias
    /// the shift takes beyond `1e100` in magnitude; failure changes nothing.
    pub fn normalize_log_weights(&mut self) -> Result<(), StructureSetError> {
        let distribution = &mut *self.distribution;
        if distribution.sample_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        let count = distribution.sample_count();
        let total = log_sum_exp(distribution.weights.materialized(count).iter().copied());
        if matches!(distribution.weights, Weights::Reweighted { .. }) {
            if let Some(position) = distribution
                .exact_log_weights()
                .iter()
                .position(|reduced| (reduced - total).abs() > REDUCED_BIAS_LIMIT)
            {
                return Err(StructureSetError::InvalidReducedBias {
                    geometry: distribution.sample(position),
                });
            }
        }
        for log_weight in distribution.weights.materialized(count) {
            *log_weight -= total;
        }
        distribution.populations = PopulationCache::default();
        Ok(())
    }

    /// Records or clears the thermodynamic state the populations describe.
    /// No weight changes. Boltzmann populations lose their tag unless the
    /// temperature stays the same; a pressure alone does not change them.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::ReweightedTarget`] for reweighted samples, whose
    /// weights belong to their target: reweight them to another temperature
    /// with [`Ensemble::reweight`](super::Ensemble::reweight). Failure
    /// changes nothing.
    pub fn set_target(
        &mut self,
        target: Option<ThermodynamicState>,
    ) -> Result<(), StructureSetError> {
        self.distribution.set_target(target)
    }

    /// Replaces the estimator of reweighted samples' group populations,
    /// keeping their reduced biases and the target.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::NotReweighted`] for any other weight model,
    /// whose populations are exact sums. Failure changes nothing.
    pub fn set_estimator(&mut self, estimator: Estimator) -> Result<(), StructureSetError> {
        match &mut self.distribution.weights {
            Weights::Reweighted {
                estimator: current, ..
            } => {
                *current = estimator;
                self.distribution.populations = PopulationCache::default();
                Ok(())
            }
            _ => Err(StructureSetError::NotReweighted),
        }
    }

    /// Sets each state's population in proportion to the Boltzmann factor
    /// of its energy relative to the lowest one, `exp(-(E - E_min) / RT)` at
    /// the temperature of `target`, and records `target`.
    ///
    /// The log weights of each state's samples shift together, so their
    /// ratios within the state are kept, and the total weight is kept. A
    /// distribution of one state therefore keeps its probabilities, with
    /// log weights `0.0` each when it had equal weights. In log space no
    /// population underflows, however far a state lies above the lowest.
    /// Samples drawn from a simulation already follow their distribution;
    /// group them into a state without an energy, which this refuses,
    /// rather than weighting them a second time.
    ///
    /// # Errors
    ///
    /// `Missing(StateEnergies)` when a state has no energy, and
    /// [`StructureSetError::UnrepresentablePopulation`] when a state lies so
    /// far above the lowest at this temperature that its log weights are not
    /// finite. Failure changes nothing.
    pub fn assign_boltzmann_populations(
        &mut self,
        target: ThermodynamicState,
    ) -> Result<(), StructureSetError> {
        let thermal_energy = *target.thermal_energy().value();
        let distribution = &mut *self.distribution;
        let energies = distribution.state_energies()?;
        let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
        let targets = energies
            .iter()
            .map(|energy| -(energy - minimum) / thermal_energy)
            .collect::<Vec<_>>();
        let target_total = log_sum_exp(targets.iter().copied());
        let populations = distribution.populations();
        // Exact shares within each state, and the exact total weight.
        let total = populations.shift + log_sum_exp(populations.sums.iter().copied());
        let log_weights = distribution
            .exact_log_weights()
            .into_iter()
            .zip(distribution.partition.assignments())
            .map(|(exact, state)| {
                (exact - populations.shift - populations.sums[*state])
                    + (targets[*state] - target_total)
                    + total
            })
            .collect::<Vec<_>>();
        if let Some(position) = log_weights
            .iter()
            .position(|log_weight| !log_weight.is_finite())
        {
            return Err(StructureSetError::UnrepresentablePopulation {
                state: distribution.partition.assignments()[position],
            });
        }
        distribution.weights = Weights::Boltzmann(log_weights);
        distribution.populations = PopulationCache::default();
        distribution.target = Some(target);
        Ok(())
    }

    /// The energy and properties of one state. Its samples and
    /// representative change only through [`Self::regroup`].
    pub fn state_mut(&mut self, state: usize) -> Option<&mut StateData> {
        self.distribution.states.get_mut(state)
    }

    /// Groups the samples into the states of `partition`, by sample
    /// position, changing no weight and no sample.
    ///
    /// The samples and their weights stay the same; only the grouping
    /// changes, and each new state's population is the estimator over its
    /// samples ([`Distribution::estimator`]), so under a cumulant estimator
    /// a merged state's population need not be the sum of its parts'. A new state made of whole old states
    /// has their free energy at the target's temperature `T`,
    /// `F = -RT ln Σ exp(-E_i / RT)` over the old states' energies `E_i`.
    /// With Boltzmann populations at `T`
    /// ([`Self::assign_boltzmann_populations`]) the populations of two such
    /// states then satisfy `F_a - F_b = -RT ln(p_a / p_b)`, and merging in
    /// steps gives the energy of merging at once. A state that takes one
    /// old state whole keeps its energy bit for bit, so an unchanged
    /// grouping keeps everything. A new state has no energy when any of its
    /// old states has none, when it takes only part of an old state that is
    /// split apart, or when it merges several old states and the
    /// distribution has no target.
    ///
    /// The free energy belongs to the temperature it was merged at:
    /// [`Self::set_target`] and [`Self::assign_boltzmann_populations`] at
    /// another temperature keep it, and assign populations from it as from
    /// any state energy. For free energies at another temperature, reweight
    /// before merging.
    ///
    /// A new state keeps the properties of its representative's old state
    /// when that state stays whole, and has none otherwise. The weights and
    /// the target are kept.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::PartitionMismatch`] when `partition` covers
    /// another number of samples than there are; failure changes nothing.
    pub fn regroup(&mut self, partition: Partition) -> Result<(), StructureSetError> {
        if partition.sample_count() != self.distribution.sample_count() {
            return Err(StructureSetError::PartitionMismatch(LengthMismatch::new(
                Countable::Sample,
                self.distribution.sample_count(),
                partition.sample_count(),
            )));
        }
        let old = self.distribution.partition.assignments();
        let new = partition.assignments();
        // The one new state every sample of an old state went to, if it
        // stayed whole.
        let mut target = vec![None; self.distribution.states.len()];
        let mut split = vec![false; self.distribution.states.len()];
        for (old, new) in old.iter().zip(new) {
            match target[*old] {
                None => target[*old] = Some(*new),
                Some(existing) if existing != *new => split[*old] = true,
                Some(_) => {}
            }
        }
        // The whole old states of each new state, and whether it takes part
        // of a split one.
        let mut members = vec![Vec::new(); partition.representatives().len()];
        let mut partial = vec![false; members.len()];
        for (state, target) in target.iter().enumerate() {
            if let Some(target) = target {
                if !split[state] {
                    members[*target].push(state);
                }
            }
        }
        for (old, new) in old.iter().zip(new) {
            partial[*new] |= split[*old];
        }
        let thermal_energy = self
            .distribution
            .target
            .map(|target| *target.thermal_energy().value());
        let energies = members
            .iter()
            .zip(&partial)
            .map(|(members, partial)| {
                if *partial {
                    return None;
                }
                match members.as_slice() {
                    [state] => self.distribution.states[*state].energy(),
                    members => {
                        let energies = members
                            .iter()
                            .map(|state| {
                                self.distribution.states[*state]
                                    .energy()
                                    .map(|energy| *energy.value())
                            })
                            .collect::<Option<Vec<_>>>()?;
                        free_energy(&energies, thermal_energy?)
                            .map(|energy| Quantity::new(energy, CANONICAL_ENERGY_UNIT))
                    }
                }
            })
            .collect::<Vec<_>>();
        let mut previous = std::mem::take(&mut self.distribution.states);
        self.distribution.states = partition
            .representatives()
            .iter()
            .zip(energies)
            .map(|(representative, energy)| {
                let state = old[*representative];
                let properties = if split[state] {
                    OwnerProperties::new()
                } else {
                    std::mem::take(&mut previous[state].properties)
                };
                StateData { energy, properties }
            })
            .collect();
        self.distribution.partition = partition;
        self.distribution.populations = PopulationCache::default();
        Ok(())
    }
}

/// `-RT ln Σ exp(-E / RT)` over some energies in
/// [`CANONICAL_ENERGY_UNIT`], relative to their lowest so that no term
/// overflows, when it is finite.
fn free_energy(energies: &[f64], thermal_energy: f64) -> Option<f64> {
    let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
    let free_energy = minimum
        - thermal_energy
            * log_sum_exp(
                energies
                    .iter()
                    .map(|energy| -(energy - minimum) / thermal_energy),
            );
    free_energy.is_finite().then_some(free_energy)
}

/// Kish's effective sample size `(Σw)² / Σw²` of the weights `w = exp(ℓ)`
/// for some finite log weights, at least one.
fn kish(log_weights: impl Iterator<Item = f64> + Clone) -> f64 {
    // Shifted by the largest log weight, the ratio is unchanged and the
    // largest term is one, so neither sum overflows or cancels.
    let maximum = log_weights.clone().fold(f64::NEG_INFINITY, f64::max);
    let shifted = log_weights.map(move |log_weight| log_weight - maximum);
    let total = log_sum_exp(shifted.clone());
    let squares = log_sum_exp(shifted.map(|log_weight| 2.0 * log_weight));
    (2.0 * total - squares).exp()
}

/// `exp(x) / Σ exp(x)` for each of some finite values, shifted by their
/// largest so that the largest term is exactly one and none overflows.
fn normalized_exponentials(values: impl Iterator<Item = f64> + Clone) -> Vec<f64> {
    let maximum = values.clone().fold(f64::NEG_INFINITY, f64::max);
    let terms = values
        .map(|value| (value - maximum).exp())
        .collect::<Vec<_>>();
    let total = terms.iter().sum::<f64>();
    terms.into_iter().map(|term| term / total).collect()
}
