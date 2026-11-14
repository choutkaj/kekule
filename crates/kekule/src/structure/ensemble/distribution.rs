//! The distribution of an [`Ensemble`](super::Ensemble): which geometries
//! are samples, their log weights, and their grouping into states.

use crate::properties::OwnerProperties;
use crate::structure::{GeometryError, Partition, Requirement, StructureSetError};
use crate::units::{Quantity, CANONICAL_ENERGY_UNIT, CANONICAL_TEMPERATURE_UNIT};

/// Which geometries of a set a [`Distribution`] covers: its samples.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Members {
    /// Every geometry of the set.
    All,
    /// The listed geometries, by strictly increasing index; the others are
    /// not samples.
    Subset(Vec<usize>),
}

/// One state of a [`Distribution`]: an optional energy and owner properties.
///
/// Which samples a state groups, and which of them represents it, belong to
/// the distribution's [`Partition`]; [`StateView`](super::StateView) reads
/// them together with the state's probability.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct State {
    energy: Option<Quantity<f64>>,
    properties: OwnerProperties,
}

impl State {
    /// A state without an energy or properties.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns this state with a finite energy.
    pub fn with_energy(mut self, energy: Quantity<f64>) -> Result<Self, GeometryError> {
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
    pub fn set_energy(&mut self, energy: Option<Quantity<f64>>) -> Result<(), GeometryError> {
        self.energy = energy
            .map(|energy| {
                let energy = energy.into_unit(CANONICAL_ENERGY_UNIT)?;
                if !energy.value().is_finite() {
                    return Err(GeometryError::NonFiniteEnergy);
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

/// Members covering geometries `members` of a set of `len`, as
/// [`Members::All`] when they are every geometry.
fn normalized(members: Vec<usize>, len: usize) -> Members {
    if members.len() == len {
        Members::All
    } else {
        Members::Subset(members)
    }
}

/// A distribution over some or all geometries of a structure set: weighted
/// samples grouped into states.
///
/// It has two levels, as a distribution of conformations does.
///
/// - A *sample* is a geometry the distribution covers ([`Members`]), with a
///   *weight*: the share of probability it carries. Only ratios of weights
///   are meaningful, so they are stored as natural logarithms, which any
///   finite value is: weights from bias reweighting or importance sampling
///   span hundreds of orders of magnitude, beyond what a linear weight can
///   hold. Without weights, every sample weighs the same.
/// - A *state* ([`State`]) is a distinct conformation: a group of samples,
///   one of them its representative. Its *population* is the total weight
///   of its samples, derived and never stored. A state may carry an
///   *energy* on one scale for the whole distribution, such as a
///   conformer's minimized energy or a computed free energy, never a
///   snapshot's instantaneous energy, which stays a property of its
///   geometry.
///
/// Every sample is in exactly one state, and every state has at least one
/// sample, its representative ([`Partition`]). The distribution may record
/// the *temperature* its populations describe.
///
/// Two index spaces meet here. Methods that take or return a `geometry` use
/// geometry indices of the set. [`Self::log_weights`], [`Self::partition`],
/// and [`Self::probabilities`] are in sample order, by sample position: the
/// `i`th sample is the `i`th geometry [`Members`] lists.
/// [`Self::sample_probabilities`] joins the two.
///
/// A conformer ensemble is one state per conformer; an NMR bundle or the
/// frames of a trajectory are one state of many samples, and clustering
/// makes several.
///
/// ```
/// use kekule::structure::{Distribution, Members, Partition};
/// // Frames 1 and 2 of a three-frame trajectory, as one state, the first
/// // weighing three times the second.
/// let distribution = Distribution::new(Members::Subset(vec![1, 2]), Partition::one_state(2))?
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
    members: Members,
    log_weights: Option<Vec<f64>>,
    partition: Partition,
    states: Vec<State>,
    temperature: Option<Quantity<f64>>,
}

impl Distribution {
    /// A distribution of `members` grouped by `partition`, without weights,
    /// state energies, state properties, or a temperature.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::UnsortedMembers`] when a subset's indices do not
    /// increase, and [`StructureSetError::PartitionMismatch`] when the
    /// partition covers another number of samples than the subset lists. A
    /// distribution of [`Members::All`] is checked against its set when an
    /// [`Ensemble`](super::Ensemble) is built.
    pub fn new(members: Members, partition: Partition) -> Result<Self, StructureSetError> {
        if let Members::Subset(indices) = &members {
            if indices.windows(2).any(|pair| pair[0] >= pair[1]) {
                return Err(StructureSetError::UnsortedMembers);
            }
            if indices.len() != partition.sample_count() {
                return Err(StructureSetError::PartitionMismatch {
                    expected: indices.len(),
                    actual: partition.sample_count(),
                });
            }
        }
        Ok(Self {
            members,
            log_weights: None,
            states: vec![State::new(); partition.state_count()],
            partition,
            temperature: None,
        })
    }

    /// Returns this distribution with one finite natural-logarithm weight
    /// per sample, in sample order; see [`DistributionMut::set_log_weights`].
    pub fn with_log_weights(mut self, log_weights: Vec<f64>) -> Result<Self, StructureSetError> {
        self.set_log_weights(Some(log_weights))?;
        Ok(self)
    }

    /// Returns this distribution with one state payload (energy and
    /// properties) per state of its partition, in state order.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::StateCountMismatch`] for another number of states.
    pub fn with_states(mut self, states: Vec<State>) -> Result<Self, StructureSetError> {
        if states.len() != self.states.len() {
            return Err(StructureSetError::StateCountMismatch {
                expected: self.states.len(),
                actual: states.len(),
            });
        }
        self.states = states;
        Ok(self)
    }

    /// Returns this distribution with a temperature; see
    /// [`DistributionMut::set_temperature`].
    pub fn with_temperature(
        mut self,
        temperature: Quantity<f64>,
    ) -> Result<Self, StructureSetError> {
        self.set_temperature(Some(temperature))?;
        Ok(self)
    }

    /// Which geometries are samples.
    pub fn members(&self) -> &Members {
        &self.members
    }

    /// Number of samples.
    pub fn sample_count(&self) -> usize {
        self.partition.sample_count()
    }

    /// The position of geometry `geometry` among the samples, which indexes
    /// the log weights and the partition, or `None` when it is not a sample.
    pub fn sample_position(&self, geometry: usize) -> Option<usize> {
        match &self.members {
            Members::All => (geometry < self.sample_count()).then_some(geometry),
            Members::Subset(indices) => indices.binary_search(&geometry).ok(),
        }
    }

    /// The geometry index of the sample at `position`.
    pub(crate) fn sample(&self, position: usize) -> usize {
        match &self.members {
            Members::All => position,
            Members::Subset(indices) => indices[position],
        }
    }

    /// The log weights in sample order, or `None` when every sample weighs
    /// the same.
    pub fn log_weights(&self) -> Option<&[f64]> {
        self.log_weights.as_deref()
    }

    /// The natural-logarithm weight of geometry `geometry`: its own, `0.0`
    /// when every sample weighs the same, or `None` when it is not a sample.
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

    /// The payload of every state, its energy and properties, in state
    /// order, as [`Self::with_states`] takes them.
    pub fn states(&self) -> &[State] {
        &self.states
    }

    /// The state of geometry `geometry`, or `None` when it is not a sample.
    pub fn state_of(&self, geometry: usize) -> Option<usize> {
        self.sample_position(geometry)
            .map(|position| self.partition.assignments()[position])
    }

    /// The temperature the populations describe, in
    /// [`CANONICAL_TEMPERATURE_UNIT`], if recorded. Changing weights does
    /// not clear it; the caller states what its populations describe.
    pub const fn temperature(&self) -> Option<Quantity<f64>> {
        self.temperature
    }

    /// This distribution in a set of `len` geometries whose first ones are
    /// those it describes, such as a file that grew after the distribution
    /// was saved. A distribution of every geometry becomes one of the
    /// geometries it covered, [`Members::Subset`] of `0..n` for `n` samples,
    /// so the geometries after them are not samples; a subset keeps its
    /// samples. A subset listing all `len` geometries becomes
    /// [`Members::All`]. The log weights, the partition, the states, and the
    /// temperature are kept bit for bit. A subset does not record how many
    /// geometries it was made for, so a `len` shorter than that is not
    /// detected.
    ///
    /// ```
    /// use kekule::structure::{Distribution, Members, Partition};
    ///
    /// let saved = Distribution::new(Members::All, Partition::one_state(3))?
    ///     .with_log_weights(vec![0.0, 1.0, 2.0])?;
    /// let grown = saved.clone().extended_to(5)?;
    /// assert_eq!(grown.members(), &Members::Subset(vec![0, 1, 2]));
    /// assert_eq!(grown.log_weights(), saved.log_weights());
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
        if self.members == Members::All {
            let covered = self.sample_count();
            if covered > len {
                return Err(StructureSetError::CoverageMismatch {
                    covered,
                    geometries: len,
                });
            }
            self.members = normalized((0..covered).collect(), len);
        }
        self.bind(len)
    }

    /// The probability of each sample, in sample order: its weight divided
    /// by the total weight. A probability below the smallest positive
    /// `f64` is `0.0`.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptyDistribution`] without samples.
    pub fn probabilities(&self) -> Result<Vec<f64>, StructureSetError> {
        if self.sample_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        Ok(normalized_exponentials(
            (0..self.sample_count()).map(|position| self.sample_log_weight(position)),
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

    /// The population share of each state, in state order: the total weight
    /// of its samples divided by the total weight.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptyDistribution`] without samples.
    pub fn state_probabilities(&self) -> Result<Vec<f64>, StructureSetError> {
        if self.sample_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        Ok(normalized_exponentials(
            self.state_log_populations().into_iter(),
        ))
    }

    /// Kish's effective sample size, `(Σw)² / Σw²`: the number of equally
    /// weighted samples that would give the same statistical precision.
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
        let Some(log_weights) = &self.log_weights else {
            return Ok(self.sample_count() as f64);
        };
        // Shifted by the largest log weight, the ratio is unchanged and the
        // largest term is one, so neither sum overflows or cancels.
        let maximum = log_weights
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let shifted = log_weights.iter().map(|log_weight| log_weight - maximum);
        let total = log_sum_exp(shifted.clone());
        let squares = log_sum_exp(shifted.map(|log_weight| 2.0 * log_weight));
        Ok((2.0 * total - squares).exp())
    }

    /// Sample geometry indices, not sample positions, by descending weight;
    /// equal weights keep geometry order.
    pub fn ranked_by_weight(&self) -> Vec<usize> {
        let mut positions = (0..self.sample_count()).collect::<Vec<_>>();
        positions.sort_by(|&a, &b| {
            self.sample_log_weight(b)
                .total_cmp(&self.sample_log_weight(a))
        });
        positions
            .into_iter()
            .map(|position| self.sample(position))
            .collect()
    }

    /// State indices by descending population; equal populations keep
    /// state order.
    pub fn ranked_by_population(&self) -> Vec<usize> {
        let populations = self.state_log_populations();
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

    fn set_temperature(
        &mut self,
        temperature: Option<Quantity<f64>>,
    ) -> Result<(), StructureSetError> {
        self.temperature = temperature
            .map(|temperature| {
                let kelvin = checked_temperature(temperature)?;
                Ok::<_, StructureSetError>(Quantity::new(kelvin, CANONICAL_TEMPERATURE_UNIT))
            })
            .transpose()?;
        Ok(())
    }

    fn set_log_weights(&mut self, log_weights: Option<Vec<f64>>) -> Result<(), StructureSetError> {
        if let Some(log_weights) = &log_weights {
            if log_weights.len() != self.sample_count() {
                return Err(StructureSetError::LogWeightCountMismatch {
                    expected: self.sample_count(),
                    actual: log_weights.len(),
                });
            }
            if let Some(position) = log_weights
                .iter()
                .position(|log_weight| !log_weight.is_finite())
            {
                return Err(StructureSetError::NonFiniteLogWeight {
                    geometry: self.sample(position),
                });
            }
        }
        self.log_weights = log_weights;
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

    /// The log weight of the sample at `position`.
    fn sample_log_weight(&self, position: usize) -> f64 {
        self.log_weights
            .as_ref()
            .map_or(0.0, |log_weights| log_weights[position])
    }

    /// The natural logarithm of each state's population, in state order, on
    /// the scale of the log weights.
    fn state_log_populations(&self) -> Vec<f64> {
        let mut members = vec![Vec::new(); self.states.len()];
        for (position, state) in self.partition.assignments().iter().enumerate() {
            members[*state].push(self.sample_log_weight(position));
        }
        members
            .into_iter()
            .map(|log_weights| log_sum_exp(log_weights.into_iter()))
            .collect()
    }

    /// The natural logarithm of one state's population, which must exist.
    pub(crate) fn state_log_population(&self, state: usize) -> f64 {
        log_sum_exp(
            self.partition
                .samples(state)
                .map(|position| self.sample_log_weight(position)),
        )
    }

    /// The payload of one state, which must exist.
    pub(crate) fn state(&self, state: usize) -> &State {
        &self.states[state]
    }

    /// Checks this distribution against a set of `len` geometries, and lists
    /// every geometry as [`Members::All`].
    pub(crate) fn bind(&mut self, len: usize) -> Result<(), StructureSetError> {
        match &self.members {
            Members::All => {
                if self.sample_count() != len {
                    return Err(StructureSetError::CoverageMismatch {
                        covered: self.sample_count(),
                        geometries: len,
                    });
                }
            }
            Members::Subset(indices) => {
                if let Some(&last) = indices.last() {
                    if last >= len {
                        return Err(StructureSetError::IndexOutOfRange { index: last, len });
                    }
                }
                if indices.len() == len {
                    self.members = Members::All;
                }
            }
        }
        Ok(())
    }

    /// Checks a new sample's log weight before [`Self::push_member`] or
    /// [`Self::push_state`].
    pub(crate) fn check_new_log_weight(
        &self,
        log_weight: Option<f64>,
        geometry: usize,
    ) -> Result<(), StructureSetError> {
        match (log_weight, &self.log_weights) {
            (None, Some(_)) => Err(StructureSetError::LogWeightRequired),
            (Some(log_weight), _) if !log_weight.is_finite() => {
                Err(StructureSetError::NonFiniteLogWeight { geometry })
            }
            _ => Ok(()),
        }
    }

    /// Checks that `state` exists.
    pub(crate) fn check_state(&self, state: usize) -> Result<(), StructureSetError> {
        if state >= self.states.len() {
            return Err(StructureSetError::StateOutOfRange {
                state,
                states: self.states.len(),
            });
        }
        Ok(())
    }

    /// Adds geometry `geometry`, the new last geometry of the set, to an
    /// existing state, after [`Self::check_state`] and
    /// [`Self::check_new_log_weight`].
    pub(crate) fn push_member(&mut self, geometry: usize, state: usize, log_weight: Option<f64>) {
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
        state: State,
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

    fn push_log_weight(&mut self, log_weight: Option<f64>) {
        match (log_weight, &mut self.log_weights) {
            (Some(log_weight), Some(log_weights)) => log_weights.push(log_weight),
            (Some(log_weight), None) => {
                let mut log_weights = vec![0.0; self.partition.sample_count()];
                log_weights.push(log_weight);
                self.log_weights = Some(log_weights);
            }
            (None, None) => {}
            (None, Some(_)) => unreachable!("new log weights are checked first"),
        }
    }

    fn push_position(&mut self, geometry: usize) {
        if let Members::Subset(indices) = &mut self.members {
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
            if let Some(log_weights) = &mut self.log_weights {
                log_weights.remove(position);
            }
        }
        if let Members::Subset(indices) = &self.members {
            let members = indices
                .iter()
                .filter(|index| **index != geometry)
                .map(|index| if *index > geometry { index - 1 } else { *index })
                .collect();
            self.members = normalized(members, len);
        }
    }

    /// The distribution of a set of the geometries `indices` (checked), in
    /// order: the selected samples, with their log weights. Every state with
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
                    State::new()
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
            members: normalized(
                picked.iter().map(|(copy, _)| *copy).collect(),
                indices.len(),
            ),
            log_weights: self.log_weights.as_ref().map(|log_weights| {
                picked
                    .iter()
                    .map(|(_, position)| log_weights[*position])
                    .collect()
            }),
            partition: Partition::from_parts_unchecked(assignments, representatives),
            states,
            temperature: self.temperature,
        }
    }

    /// The samples of whole states `indices` (checked), in the requested
    /// order and each state's samples in geometry order, as their geometry
    /// indices, with the distribution of a set of just those geometries.
    pub(crate) fn select_states(&self, indices: &[usize]) -> (Vec<usize>, Self) {
        let mut geometries = Vec::new();
        let mut log_weights = self.log_weights.as_ref().map(|_| Vec::new());
        let mut assignments = Vec::new();
        let mut representatives = Vec::with_capacity(indices.len());
        let mut states = Vec::with_capacity(indices.len());
        for (state, old) in indices.iter().enumerate() {
            let representative = self.partition.representatives()[*old];
            for position in self.partition.samples(*old) {
                if position == representative {
                    representatives.push(geometries.len());
                }
                geometries.push(self.sample(position));
                if let Some(log_weights) = &mut log_weights {
                    log_weights.push(self.sample_log_weight(position));
                }
                assignments.push(state);
            }
            states.push(self.states[*old].clone());
        }
        let distribution = Self {
            members: Members::All,
            log_weights,
            partition: Partition::from_parts_unchecked(assignments, representatives),
            states,
            temperature: self.temperature,
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
            members: Members::All,
            log_weights: Some(self.state_log_populations()),
            partition: Partition::state_per_sample(representatives.len()),
            states: self.states.clone(),
            temperature: self.temperature,
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
/// log weights, the grouping into states, the states' energies and
/// properties, and the temperature, but never which geometries are samples,
/// so the distribution always fits its set.
///
/// ```compile_fail,E0594
/// use kekule::structure::{Distribution, Ensemble, Members, Partition};
/// fn replace(ensemble: &mut Ensemble) {
///     let foreign = Distribution::new(Members::Subset(vec![7]), Partition::one_state(1));
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
    /// Replaces every log weight, in sample order, or makes every sample
    /// weigh the same. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::LogWeightCountMismatch`] for another number of
    /// log weights than samples, and [`StructureSetError::NonFiniteLogWeight`]
    /// for the first one that is not finite.
    pub fn set_log_weights(
        &mut self,
        log_weights: Option<Vec<f64>>,
    ) -> Result<(), StructureSetError> {
        self.distribution.set_log_weights(log_weights)
    }

    /// Sets the log weight of one sample, a finite value. In a distribution
    /// without weights, the other samples keep log weight `0.0`.
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
        distribution
            .log_weights
            .get_or_insert_with(|| vec![0.0; count])[position] = log_weight;
        Ok(())
    }

    /// Shifts every log weight so the weights sum to one, so each log weight
    /// is the natural logarithm of its sample's probability. Ratios, and so
    /// every probability and population share, are kept.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptyDistribution`] without samples; failure
    /// changes nothing.
    pub fn normalize_log_weights(&mut self) -> Result<(), StructureSetError> {
        let distribution = &mut *self.distribution;
        if distribution.sample_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        let log_weights = (0..distribution.sample_count())
            .map(|position| distribution.sample_log_weight(position))
            .collect::<Vec<_>>();
        let total = log_sum_exp(log_weights.iter().copied());
        distribution.log_weights = Some(
            log_weights
                .into_iter()
                .map(|log_weight| log_weight - total)
                .collect(),
        );
        Ok(())
    }

    /// Records or clears the temperature the populations describe, which
    /// must be finite and positive. No weight changes.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::InvalidTemperature`] and
    /// [`StructureSetError::Unit`]; failure changes nothing.
    pub fn set_temperature(
        &mut self,
        temperature: Option<Quantity<f64>>,
    ) -> Result<(), StructureSetError> {
        self.distribution.set_temperature(temperature)
    }

    /// Sets each state's population in proportion to the Boltzmann factor
    /// of its energy relative to the lowest one, `exp(-(E - E_min) / RT)`,
    /// and records `temperature`.
    ///
    /// The log weights of each state's samples shift together, so their
    /// ratios within the state are kept, and the total weight is kept. A
    /// distribution of one state is therefore unchanged, apart from gaining
    /// log weights (`0.0` each) when it had none. In log space no
    /// population underflows, however far a state lies above the lowest.
    /// Samples drawn from a simulation already follow their distribution;
    /// group them into a state without an energy, which this refuses,
    /// rather than weighting them a second time.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::InvalidTemperature`] and
    /// [`StructureSetError::Unit`] for a temperature that is not finite and
    /// positive, `Missing(StateEnergies)` when a state has no energy, and
    /// [`StructureSetError::UnrepresentablePopulation`] when a state lies so
    /// far above the lowest at this temperature that its log weights are not
    /// finite. Failure changes nothing.
    pub fn assign_boltzmann_populations(
        &mut self,
        temperature: Quantity<f64>,
    ) -> Result<(), StructureSetError> {
        let kelvin = checked_temperature(temperature)?;
        let thermal_energy = thermal_energy(kelvin).ok_or(StructureSetError::InvalidTemperature)?;
        let distribution = &mut *self.distribution;
        let energies = distribution.state_energies()?;
        let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
        let targets = energies
            .iter()
            .map(|energy| -(energy - minimum) / thermal_energy)
            .collect::<Vec<_>>();
        let target_total = log_sum_exp(targets.iter().copied());
        let populations = distribution.state_log_populations();
        let total = log_sum_exp(populations.iter().copied());
        let log_weights = distribution
            .partition
            .assignments()
            .iter()
            .enumerate()
            .map(|(position, state)| {
                (distribution.sample_log_weight(position) - populations[*state])
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
        distribution.log_weights = Some(log_weights);
        distribution.temperature = Some(Quantity::new(kelvin, CANONICAL_TEMPERATURE_UNIT));
        Ok(())
    }

    /// The energy and properties of one state. Its samples and
    /// representative change only through [`Self::regroup`].
    pub fn state_mut(&mut self, state: usize) -> Option<&mut State> {
        self.distribution.states.get_mut(state)
    }

    /// Groups the samples into the states of `partition`, by sample
    /// position, changing no weight and no sample.
    ///
    /// The distribution stays the same; only its description changes, and
    /// populations follow the samples. A new state that contains every
    /// sample of an old state keeps the energy and properties of its
    /// representative's old state, so merging whole states keeps the
    /// representative's energy, and an unchanged grouping keeps everything.
    /// A state whose old state is split apart has neither. The temperature
    /// is kept.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::PartitionMismatch`] when `partition` covers
    /// another number of samples than there are; failure changes nothing.
    pub fn regroup(&mut self, partition: Partition) -> Result<(), StructureSetError> {
        if partition.sample_count() != self.distribution.sample_count() {
            return Err(StructureSetError::PartitionMismatch {
                expected: self.distribution.sample_count(),
                actual: partition.sample_count(),
            });
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
        let mut previous = std::mem::take(&mut self.distribution.states);
        self.distribution.states = partition
            .representatives()
            .iter()
            .map(|representative| {
                let state = old[*representative];
                if split[state] {
                    State::new()
                } else {
                    std::mem::take(&mut previous[state])
                }
            })
            .collect();
        self.distribution.partition = partition;
        Ok(())
    }
}

/// The molar gas constant `N_A k_B` in kJ/(mol K), exact since the 2019 SI.
const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24e-3;

/// `R T` for a temperature in kelvin, when it is finite and positive.
fn thermal_energy(kelvin: f64) -> Option<f64> {
    let thermal_energy = MOLAR_GAS_CONSTANT * kelvin;
    (thermal_energy.is_finite() && thermal_energy > 0.0).then_some(thermal_energy)
}

/// `ln Σ exp(x)` over finite values, shifted by their largest so that no
/// term overflows; negative infinity for none.
pub(crate) fn log_sum_exp(values: impl Iterator<Item = f64> + Clone) -> f64 {
    let maximum = values.clone().fold(f64::NEG_INFINITY, f64::max);
    if maximum == f64::NEG_INFINITY {
        return maximum;
    }
    maximum
        + values
            .map(|value| (value - maximum).exp())
            .sum::<f64>()
            .ln()
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

/// A finite, positive temperature in [`CANONICAL_TEMPERATURE_UNIT`].
fn checked_temperature(temperature: Quantity<f64>) -> Result<f64, StructureSetError> {
    let temperature = temperature.value_in(CANONICAL_TEMPERATURE_UNIT)?;
    if temperature.is_finite() && temperature > 0.0 {
        Ok(temperature)
    } else {
        Err(StructureSetError::InvalidTemperature)
    }
}
