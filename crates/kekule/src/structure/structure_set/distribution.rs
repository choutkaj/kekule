//! The distribution of a [`super::StructureSet`]: which samples are members,
//! their weights, and their grouping into states.

use crate::properties::OwnerProperties;
use crate::structure::{GeometryError, Partition};
use crate::units::{Quantity, CANONICAL_ENERGY_UNIT, CANONICAL_TEMPERATURE_UNIT};

use super::{Membership, Requirement, StructureSetError};

/// Which samples of a set a [`Distribution`] covers.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Members {
    /// Every sample of the set.
    All,
    /// The listed samples, by strictly increasing index; the others have no
    /// membership.
    Subset(Vec<usize>),
}

/// One state of a [`Distribution`]: an optional energy and owner properties.
///
/// Which samples a state groups, and which of them represents it, belong to
/// the distribution's [`Partition`]; [`StateView`](super::StateView) reads
/// them together with the state's population.
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

/// Members covering samples `members` of a set of `len`, as [`Members::All`]
/// when they are every sample.
fn normalized(members: Vec<usize>, len: usize) -> Members {
    if members.len() == len {
        Members::All
    } else {
        Members::Subset(members)
    }
}

/// A distribution over some or all samples of a set: weighted members
/// grouped into states.
///
/// It has two levels, as a distribution of conformations does.
///
/// - A *member* is a sample with a finite, positive, relative *weight*: the
///   share of probability it carries. Only ratios of weights are
///   meaningful; without weights, every member weighs the same.
/// - A *state* ([`State`]) is a distinct conformation: a group of
///   members, one of them its representative. Its *population* is the total
///   weight of its members, derived and never stored. A state may carry an
///   *energy* on one scale for the whole distribution, such as a
///   conformer's minimized energy or a computed free energy. It is not a
///   sample's instantaneous energy, which stays a property of the sample.
///
/// Every member is in exactly one state, and every state has at least one
/// member, its representative among them ([`Partition`]). The
/// distribution may record the *temperature* its populations describe.
///
/// Two index spaces meet here. Methods that take or return a `sample` use
/// sample indices of the set. [`Self::weights`], [`Self::partition`], and
/// [`Self::probabilities`] are in member order, by member position: the
/// `i`th member is the `i`th sample [`Members`] lists.
/// [`Self::weighted_members`] joins the two, giving each member's sample
/// index with its weight.
///
/// A conformer set is one state per conformer; an NMR bundle or the frames
/// of a trajectory are one state of many members, and clustering makes
/// several. A distribution installed in a set is edited through the
/// [`DistributionMut`] guard of [`super::StructureSet::distribution_mut`],
/// which never changes which samples are members; the set's own editing,
/// selecting, and removing do, and remap the members.
///
/// ```
/// use kekule::structure::{Distribution, Partition, Members};
/// // Frames 1 and 2 of a three-frame trajectory, as one state.
/// let distribution = Distribution::new(
///     Members::Subset(vec![1, 2]),
///     Partition::one_state(2),
/// )?
/// .with_weights(vec![3.0, 1.0])?;
/// assert_eq!(distribution.weight(1), Some(3.0));
/// assert_eq!(distribution.weight(0), None);
/// assert_eq!(distribution.probabilities()?, [0.75, 0.25]);
/// assert_eq!(
///     distribution.weighted_members().collect::<Vec<_>>(),
///     [(1, 3.0), (2, 1.0)],
/// );
/// assert_eq!(distribution.ranked_by_weight(), [1, 2]);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Distribution {
    pub(super) members: Members,
    pub(super) weights: Option<Vec<f64>>,
    pub(super) partition: Partition,
    pub(super) states: Vec<State>,
    pub(super) temperature: Option<Quantity<f64>>,
}

impl Distribution {
    /// A distribution of `members` grouped by `partition`, without weights,
    /// state energies, state properties, or a temperature.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::UnsortedMembers`] when a subset's indices do not
    /// increase, and [`StructureSetError::PartitionMismatch`] when the partition
    /// covers another number of samples than the subset lists. A
    /// distribution of [`Members::All`] is checked against its set when it
    /// is installed.
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
            weights: None,
            states: vec![State::new(); partition.state_count()],
            partition,
            temperature: None,
        })
    }

    /// Returns this distribution with one finite, positive weight per
    /// member, in member order; see [`DistributionMut::set_weights`].
    pub fn with_weights(mut self, weights: Vec<f64>) -> Result<Self, StructureSetError> {
        self.set_weights(Some(weights))?;
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

    /// Which samples are members.
    pub fn members(&self) -> &Members {
        &self.members
    }

    /// Number of members.
    pub fn member_count(&self) -> usize {
        self.partition.sample_count()
    }

    /// The position of sample `sample` among the members, which indexes the
    /// weights and the partition, or `None` when it is not a member.
    pub fn member_position(&self, sample: usize) -> Option<usize> {
        match &self.members {
            Members::All => (sample < self.member_count()).then_some(sample),
            Members::Subset(indices) => indices.binary_search(&sample).ok(),
        }
    }

    /// Each member's sample index with its relative weight, in member
    /// order, so the `i`th item belongs to the `i`th weight, probability,
    /// and partition entry. Every weight is `1.0` when the members weigh the
    /// same. This is what weighted statistics over the distribution see;
    /// non-members are left out.
    pub fn weighted_members(&self) -> impl ExactSizeIterator<Item = (usize, f64)> + '_ {
        (0..self.member_count())
            .map(move |position| (self.member(position), self.member_weight(position)))
    }

    /// The sample index of the member at `position`.
    pub(crate) fn member(&self, position: usize) -> usize {
        match &self.members {
            Members::All => position,
            Members::Subset(indices) => indices[position],
        }
    }

    /// The weights in member order, or `None` when every member weighs the
    /// same.
    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// The relative weight of sample `sample`: its own, `1.0` when every
    /// member weighs the same, or `None` when it is not a member.
    pub fn weight(&self, sample: usize) -> Option<f64> {
        self.member_position(sample)
            .map(|position| self.member_weight(position))
    }

    /// The grouping of the members into states, by member position.
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

    /// This distribution in a set of `len` samples whose first samples are
    /// the ones it describes, such as a file that grew after the
    /// distribution was saved. A distribution of every sample becomes one
    /// of the samples it covered, [`Members::Subset`] of `0..n` for `n`
    /// members, so the samples after them are outside it; a subset keeps
    /// its members. A subset listing all `len` samples becomes
    /// [`Members::All`]. The weights, the partition, the states with their
    /// energies and properties, and the temperature are kept bit for bit.
    /// [`StructureSet::refresh`](super::StructureSet::refresh) gives a set's
    /// distribution so when the file grows. A subset does not record how
    /// many samples it was made for, so a `len` shorter than that is not
    /// detected; [`StructureSetDescription::extended_to`](crate::native::StructureSetDescription::extended_to)
    /// checks it against the saved sample count.
    ///
    /// ```
    /// use kekule::structure::{Distribution, Members, Partition};
    ///
    /// let saved = Distribution::new(Members::All, Partition::one_state(3))?
    ///     .with_weights(vec![1.0, 2.0, 3.0])?;
    /// let grown = saved.clone().extended_to(5)?;
    /// assert_eq!(grown.members(), &Members::Subset(vec![0, 1, 2]));
    /// assert_eq!(grown.weights(), saved.weights());
    /// assert_eq!(grown.weight(4), None);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`StructureSetError::CoverageMismatch`] when a distribution of every
    /// sample covers more than `len`, and
    /// [`StructureSetError::IndexOutOfRange`] when a subset lists a sample past
    /// them.
    pub fn extended_to(mut self, len: usize) -> Result<Self, StructureSetError> {
        self.extend(len)?;
        Ok(self)
    }

    /// Re-expresses this distribution in place for a set of `len` samples,
    /// as [`Self::extended_to`] does; failure changes nothing.
    pub(crate) fn extend(&mut self, len: usize) -> Result<(), StructureSetError> {
        if self.members == Members::All {
            let covered = self.member_count();
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

    /// The state of sample `sample`, or `None` when it is not a member.
    pub fn state_of(&self, sample: usize) -> Option<usize> {
        self.member_position(sample)
            .map(|position| self.partition.assignments()[position])
    }

    /// The temperature the populations describe, in
    /// [`CANONICAL_TEMPERATURE_UNIT`], if recorded. Changing weights does
    /// not clear it; the caller states what its populations describe.
    pub const fn temperature(&self) -> Option<Quantity<f64>> {
        self.temperature
    }

    /// The probability of each member, in member order: its weight divided
    /// by the total weight.
    ///
    /// Fails without members or when a probability would underflow to zero.
    pub fn probabilities(&self) -> Result<Vec<f64>, StructureSetError> {
        if self.member_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        let scaled = self.scaled_weights();
        let total = scaled.iter().sum::<f64>();
        scaled
            .into_iter()
            .enumerate()
            .map(|(position, weight)| self.checked_weight(weight / total, position))
            .collect()
    }

    /// The population share of each state, in state order: the total weight
    /// of its members divided by the total weight. The shares sum to one.
    ///
    /// Fails without members or when a share would underflow to zero.
    pub fn state_probabilities(&self) -> Result<Vec<f64>, StructureSetError> {
        if self.member_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        let scaled = self.scaled_state_populations();
        let total = scaled.iter().sum::<f64>();
        scaled
            .into_iter()
            .enumerate()
            .map(|(state, population)| {
                let representative = self.member(self.partition.representatives()[state]);
                checked_weight(population / total).ok_or(StructureSetError::InvalidWeight {
                    geometry: representative,
                })
            })
            .collect()
    }

    /// Kish's effective sample size over members, `(Σw)² / Σw²`: the number
    /// of equally weighted members that would give the same statistical
    /// precision.
    ///
    /// It is the member count for equal weights and approaches one when a
    /// single member dominates. Fails without members.
    pub fn effective_sample_size(&self) -> Result<f64, StructureSetError> {
        if self.member_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        let scaled = self.scaled_weights();
        let total = scaled.iter().sum::<f64>();
        let squares = scaled.iter().map(|weight| weight * weight).sum::<f64>();
        Ok(total * total / squares)
    }

    /// Member sample indices, not member positions, by descending weight;
    /// equal weights keep sample order.
    ///
    /// Keep the ten most probable members with
    /// `set.select(distribution.ranked_by_weight().into_iter().take(10))`.
    pub fn ranked_by_weight(&self) -> Vec<usize> {
        let mut positions = (0..self.member_count()).collect::<Vec<_>>();
        positions.sort_by(|&a, &b| self.member_weight(b).total_cmp(&self.member_weight(a)));
        positions
            .into_iter()
            .map(|position| self.member(position))
            .collect()
    }

    /// State indices by descending population; equal populations keep
    /// state order.
    ///
    /// Keep the ten most populated states with
    /// `set.select_states(distribution.ranked_by_population().into_iter().take(10))`.
    pub fn ranked_by_population(&self) -> Vec<usize> {
        let populations = self.scaled_state_populations();
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
    /// [`super::StructureSet::select_states`] to keep, for example, 95% of a
    /// Boltzmann population. Fails without members.
    pub fn covering_population(&self, fraction: f64) -> Result<Vec<usize>, StructureSetError> {
        if !(fraction > 0.0 && fraction <= 1.0) {
            return Err(StructureSetError::InvalidFraction);
        }
        if self.member_count() == 0 {
            return Err(StructureSetError::EmptyDistribution);
        }
        let ranked = self.ranked_by_population();
        if fraction == 1.0 {
            return Ok(ranked);
        }
        let scaled = self.scaled_state_populations();
        // Sum in ranked order so the running total ends exactly at the total.
        let total = ranked.iter().map(|&index| scaled[index]).sum::<f64>();
        let target = fraction * total;
        let mut covered = 0.0;
        let count = ranked
            .iter()
            .position(|&index| {
                covered += scaled[index];
                covered >= target
            })
            .map_or(ranked.len(), |position| position + 1);
        Ok(ranked[..count].to_vec())
    }

    /// State indices by ascending energy; equal energies keep state order.
    ///
    /// Fails with `Missing(StateEnergies)` if any state has no energy.
    pub fn ranked_by_energy(&self) -> Result<Vec<usize>, StructureSetError> {
        let energies = self.state_energies()?;
        let mut indices = (0..energies.len()).collect::<Vec<_>>();
        indices.sort_by(|&a, &b| {
            energies[a]
                .partial_cmp(&energies[b])
                .expect("state energies are finite")
        });
        Ok(indices)
    }

    /// Indices, in state order, of the states whose energy lies at most
    /// `window` above the lowest state energy.
    ///
    /// `window` is a finite, nonnegative energy, such as the 6 kcal/mol
    /// window common in conformer searches. Fails with
    /// `Missing(StateEnergies)` if any state has no energy.
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

    fn set_weights(&mut self, weights: Option<Vec<f64>>) -> Result<(), StructureSetError> {
        if let Some(weights) = &weights {
            if weights.len() != self.member_count() {
                return Err(StructureSetError::WeightCountMismatch {
                    expected: self.member_count(),
                    actual: weights.len(),
                });
            }
            for (position, weight) in weights.iter().enumerate() {
                self.checked_weight(*weight, position)?;
            }
        }
        self.weights = weights;
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

    /// The weight of the member at `position`.
    pub(super) fn member_weight(&self, position: usize) -> f64 {
        self.weights
            .as_ref()
            .map_or(1.0, |weights| weights[position])
    }

    /// Weights divided by the largest one, so sums cannot overflow.
    pub(crate) fn scaled_weights(&self) -> Vec<f64> {
        let maximum = self.largest_weight();
        (0..self.member_count())
            .map(|position| self.member_weight(position) / maximum)
            .collect()
    }

    /// The summed [`Self::scaled_weights`] of each state.
    fn scaled_state_populations(&self) -> Vec<f64> {
        let mut populations = vec![0.0; self.states.len()];
        for (weight, state) in self
            .scaled_weights()
            .into_iter()
            .zip(self.partition.assignments())
        {
            populations[*state] += weight;
        }
        populations
    }

    /// The population of one state, which must exist, relative to the
    /// largest member weight, so that it cannot overflow.
    pub(crate) fn population(&self, state: usize) -> f64 {
        let maximum = self.largest_weight();
        self.partition
            .samples(state)
            .map(|position| self.member_weight(position) / maximum)
            .sum()
    }

    /// The largest member weight, `1.0` when the members weigh the same or
    /// there are none.
    fn largest_weight(&self) -> f64 {
        self.weights
            .as_ref()
            .and_then(|weights| weights.iter().copied().reduce(f64::max))
            .unwrap_or(1.0)
    }

    /// The payload of one state, which must exist.
    pub(crate) fn state(&self, state: usize) -> &State {
        &self.states[state]
    }

    /// What a view of sample `sample` reports, or `None` for a non-member.
    pub(crate) fn membership(&self, sample: usize) -> Option<Membership> {
        let position = self.member_position(sample)?;
        let state = self.partition.assignments()[position];
        Some(Membership::new(
            self.member_weight(position),
            state,
            self.partition.representatives()[state] == position,
            self.states[state].energy(),
            self.temperature,
        ))
    }

    /// `weight` when it is finite and positive; the error names the member
    /// at `position`.
    fn checked_weight(&self, weight: f64, position: usize) -> Result<f64, StructureSetError> {
        checked_weight(weight).ok_or(StructureSetError::InvalidWeight {
            geometry: self.member(position),
        })
    }

    /// Checks this distribution against a set of `len` samples, and lists
    /// every sample as [`Members::All`].
    pub(crate) fn bind(&mut self, len: usize) -> Result<(), StructureSetError> {
        match &self.members {
            Members::All => {
                if self.member_count() != len {
                    return Err(StructureSetError::CoverageMismatch {
                        covered: self.member_count(),
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

    /// Checks a new member's weight before [`Self::push_member`] or
    /// [`Self::push_state`].
    pub(crate) fn check_new_weight(
        &self,
        weight: Option<f64>,
        sample: usize,
    ) -> Result<(), StructureSetError> {
        match (weight, &self.weights) {
            (None, Some(_)) => Err(StructureSetError::WeightRequired),
            (Some(weight), _) => checked_weight(weight)
                .map(|_| ())
                .ok_or(StructureSetError::InvalidWeight { geometry: sample }),
            (None, None) => Ok(()),
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

    /// Adds sample `sample`, the new last sample of the set, to an existing
    /// state, after [`Self::check_state`] and [`Self::check_new_weight`].
    pub(crate) fn push_member(&mut self, sample: usize, state: usize, weight: Option<f64>) {
        self.push_weight(weight);
        self.push_position(sample);
        let (mut assignments, representatives) = self.take_partition();
        assignments.push(state);
        self.partition = Partition::from_parts_unchecked(assignments, representatives);
    }

    /// Adds sample `sample`, the new last sample of the set, as the
    /// representative of a new state, after [`Self::check_new_weight`];
    /// returns the state's index.
    pub(crate) fn push_state(&mut self, sample: usize, state: State, weight: Option<f64>) -> usize {
        let index = self.states.len();
        let position = self.member_count();
        self.push_weight(weight);
        self.push_position(sample);
        let (mut assignments, mut representatives) = self.take_partition();
        assignments.push(index);
        representatives.push(position);
        self.partition = Partition::from_parts_unchecked(assignments, representatives);
        self.states.push(state);
        index
    }

    fn push_weight(&mut self, weight: Option<f64>) {
        match (weight, &mut self.weights) {
            (Some(weight), Some(weights)) => weights.push(weight),
            (Some(weight), None) => {
                let mut weights = vec![1.0; self.partition.sample_count()];
                weights.push(weight);
                self.weights = Some(weights);
            }
            (None, None) => {}
            (None, Some(_)) => unreachable!("new weights are checked first"),
        }
    }

    fn push_position(&mut self, sample: usize) {
        if let Members::Subset(indices) = &mut self.members {
            indices.push(sample);
        }
    }

    /// Checks that sample `sample` can be removed: it is not the
    /// representative of a state with other members.
    pub(crate) fn check_remove(&self, sample: usize) -> Result<(), StructureSetError> {
        let Some(position) = self.member_position(sample) else {
            return Ok(());
        };
        let state = self.partition.assignments()[position];
        let representative = self.partition.representatives()[state] == position;
        let alone = self.partition.samples(state).nth(1).is_none();
        if representative && !alone {
            return Err(StructureSetError::RepresentativeRemoval {
                geometry: sample,
                state,
            });
        }
        Ok(())
    }

    /// Removes sample `sample` from a set that keeps `len` samples after
    /// [`Self::check_remove`], shifting later members down. A member alone
    /// in its state removes the state.
    pub(crate) fn remove(&mut self, sample: usize, len: usize) {
        if let Some(position) = self.member_position(sample) {
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
            if let Some(weights) = &mut self.weights {
                weights.remove(position);
            }
        }
        if let Members::Subset(indices) = &self.members {
            let members = indices
                .iter()
                .filter(|index| **index != sample)
                .map(|index| if *index > sample { index - 1 } else { *index })
                .collect();
            self.members = normalized(members, len);
        }
    }

    /// The distribution of a set of the samples `indices` (checked), in
    /// order: the selected members, with their weights. Every state with a
    /// selected member becomes a state of the copy, in state order,
    /// represented by the first copy of its representative or, when that is
    /// not selected, by its first selected member. A state whose members are
    /// all selected keeps its energy and properties; one selected only in
    /// part has neither.
    pub(crate) fn select(&self, indices: &[usize]) -> Self {
        // (sample index in the copy, member position here)
        let picked = indices
            .iter()
            .enumerate()
            .filter_map(|(copy, sample)| {
                self.member_position(*sample)
                    .map(|position| (copy, position))
            })
            .collect::<Vec<_>>();
        let old = self.partition.assignments();
        let mut selected = vec![0; self.states.len()];
        let mut chosen = vec![false; self.member_count()];
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
        for (member, (_, position)) in picked.iter().enumerate() {
            let state = renumbered[old[*position]];
            assignments.push(state);
            first[state].get_or_insert(member);
            if *position == self.partition.representatives()[old[*position]] {
                representatives[state].get_or_insert(member);
            }
        }
        let representatives = representatives
            .into_iter()
            .zip(first)
            .map(|(representative, first)| {
                representative
                    .or(first)
                    .expect("every renumbered state has a selected member")
            })
            .collect();
        Self {
            members: normalized(
                picked.iter().map(|(copy, _)| *copy).collect(),
                indices.len(),
            ),
            weights: self.weights.as_ref().map(|weights| {
                picked
                    .iter()
                    .map(|(_, position)| weights[*position])
                    .collect()
            }),
            partition: Partition::from_parts_unchecked(assignments, representatives),
            states,
            temperature: self.temperature,
        }
    }

    /// The members of whole states `indices` (checked), in the requested
    /// order and each state's members in sample order, as their sample
    /// indices, with the distribution of a set of just those samples.
    pub(crate) fn select_states(&self, indices: &[usize]) -> (Vec<usize>, Self) {
        let mut samples = Vec::new();
        let mut weights = self.weights.as_ref().map(|_| Vec::new());
        let mut assignments = Vec::new();
        let mut representatives = Vec::with_capacity(indices.len());
        let mut states = Vec::with_capacity(indices.len());
        for (state, old) in indices.iter().enumerate() {
            let representative = self.partition.representatives()[*old];
            for position in self.partition.samples(*old) {
                if position == representative {
                    representatives.push(samples.len());
                }
                samples.push(self.member(position));
                if let Some(weights) = &mut weights {
                    weights.push(self.member_weight(position));
                }
                assignments.push(state);
            }
            states.push(self.states[*old].clone());
        }
        let distribution = Self {
            members: Members::All,
            weights,
            partition: Partition::from_parts_unchecked(assignments, representatives),
            states,
            temperature: self.temperature,
        };
        (samples, distribution)
    }

    /// Each state's representative, as a sample index, with the
    /// distribution of a set of just those samples: one state each, weighted
    /// by its population relative to the largest member weight.
    pub(crate) fn collapse(&self) -> Result<(Vec<usize>, Self), StructureSetError> {
        let representatives = self
            .partition
            .representatives()
            .iter()
            .map(|position| self.member(*position))
            .collect::<Vec<_>>();
        let weights = self
            .scaled_state_populations()
            .into_iter()
            .zip(&representatives)
            .map(|(population, sample)| {
                checked_weight(population)
                    .ok_or(StructureSetError::InvalidWeight { geometry: *sample })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let distribution = Self {
            members: Members::All,
            weights: Some(weights),
            partition: Partition::state_per_sample(representatives.len()),
            states: self.states.clone(),
            temperature: self.temperature,
        };
        Ok((representatives, distribution))
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

/// Mutable access to the distribution of a [`super::StructureSet`]
/// ([`super::StructureSet::distribution_mut`]).
///
/// Reads use `Deref`. There is deliberately no `DerefMut`: the guard changes
/// weights, the grouping into states, the states' energies and properties,
/// and the temperature, but never which samples are members, so the
/// distribution always fits its set. Replace or clear the whole
/// distribution with [`super::StructureSet::set_distribution`], which checks
/// it against the set.
///
/// ```compile_fail,E0594
/// use kekule::structure::{Distribution, Partition, Members, StructureSet};
/// fn replace(set: &mut StructureSet) {
///     let foreign = Distribution::new(Members::Subset(vec![7]), Partition::one_state(1));
///     *set.distribution_mut().unwrap() = foreign.unwrap();
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
    /// Replaces every weight, in member order, or makes every member weigh
    /// the same. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::WeightCountMismatch`] for another number of weights
    /// than members, and [`StructureSetError::InvalidWeight`] for the first
    /// weight that is not finite and positive.
    pub fn set_weights(&mut self, weights: Option<Vec<f64>>) -> Result<(), StructureSetError> {
        self.distribution.set_weights(weights)
    }

    /// Sets the weight of one member, a finite positive value. In a
    /// distribution without weights, the other members keep weight `1.0`.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::NotMember`] when `sample` is not a member, and
    /// [`StructureSetError::InvalidWeight`] for a weight that is not finite and
    /// positive.
    pub fn set_weight(&mut self, sample: usize, weight: f64) -> Result<(), StructureSetError> {
        let distribution = &mut *self.distribution;
        let position = distribution
            .member_position(sample)
            .ok_or(StructureSetError::NotMember { geometry: sample })?;
        let weight = distribution.checked_weight(weight, position)?;
        let count = distribution.member_count();
        distribution.weights.get_or_insert_with(|| vec![1.0; count])[position] = weight;
        Ok(())
    }

    /// Scales every weight so the weights sum to one, keeping their ratios
    /// and therefore every probability and population share.
    ///
    /// Fails without members, or when a scaled weight would underflow to
    /// zero. Failure changes nothing.
    pub fn normalize_weights(&mut self) -> Result<(), StructureSetError> {
        self.distribution.weights = Some(self.distribution.probabilities()?);
        Ok(())
    }

    /// Records or clears the temperature the populations describe, which
    /// must be finite and positive. No weight changes.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::InvalidTemperature`] and [`StructureSetError::Unit`];
    /// failure changes nothing.
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
    /// The weights of each state's members are scaled together, so their
    /// ratios within the state are kept, and the total weight is kept. A
    /// distribution of one state is therefore unchanged, apart from gaining
    /// weights (`1.0` each) when it had none. Samples drawn from a simulation
    /// already follow their distribution; group them into a state without an
    /// energy, which this refuses, rather than weighting them a second time.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::InvalidTemperature`] and [`StructureSetError::Unit`]
    /// for a temperature that is not finite and positive,
    /// `Missing(StateEnergies)` when a state has no energy, and
    /// [`StructureSetError::InvalidWeight`] when a weight would underflow to
    /// zero (keep only the states within an energy window first, with
    /// [`Distribution::within_energy_window`] and
    /// [`super::StructureSet::select_states`]). Failure changes nothing.
    pub fn assign_boltzmann_populations(
        &mut self,
        temperature: Quantity<f64>,
    ) -> Result<(), StructureSetError> {
        let kelvin = checked_temperature(temperature)?;
        let thermal_energy = thermal_energy(kelvin).ok_or(StructureSetError::InvalidTemperature)?;
        let distribution = &mut *self.distribution;
        let energies = distribution.state_energies()?;
        let weights = (0..distribution.member_count())
            .map(|position| distribution.member_weight(position))
            .collect::<Vec<_>>();
        let weights = boltzmann_weights(
            thermal_energy,
            &energies,
            &weights,
            distribution.partition.assignments(),
        )
        .map_err(|position| StructureSetError::InvalidWeight {
            geometry: distribution.member(position),
        })?;
        distribution.weights = Some(weights);
        distribution.temperature = Some(Quantity::new(kelvin, CANONICAL_TEMPERATURE_UNIT));
        Ok(())
    }

    /// The energy and properties of one state. Its members and
    /// representative change only through [`Self::regroup`].
    pub fn state_mut(&mut self, state: usize) -> Option<&mut State> {
        self.distribution.states.get_mut(state)
    }

    /// Groups the members into the states of `partition`, by member
    /// position, changing no weight and no member.
    ///
    /// The distribution stays the same; only its description changes, and
    /// populations follow the members. A new state that contains every
    /// member of an old state keeps the energy and properties of its
    /// representative's old state, so merging whole states keeps the
    /// representative's energy, and an unchanged grouping keeps everything.
    /// A state whose old state is split apart has neither. The temperature
    /// is kept.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::PartitionMismatch`] when `partition` covers another
    /// number of samples than there are members; failure changes nothing.
    pub fn regroup(&mut self, partition: Partition) -> Result<(), StructureSetError> {
        if partition.sample_count() != self.distribution.member_count() {
            return Err(StructureSetError::PartitionMismatch {
                expected: self.distribution.member_count(),
                actual: partition.sample_count(),
            });
        }
        let old = self.distribution.partition.assignments();
        let new = partition.assignments();
        // The one new state every member of an old state went to, if it
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

/// New relative weights giving each state the population share
/// `exp(-(E - E_min) / RT)` normalized over the states, at `thermal_energy`
/// `RT`, from state `energies` in [`CANONICAL_ENERGY_UNIT`] and member
/// `weights` grouped into states by `assignments`.
///
/// Each state's members scale by one factor, so their ratios within the
/// state are kept, and the total weight is kept. Fails with the position of
/// the first member whose weight would not be finite and positive.
fn boltzmann_weights(
    thermal_energy: f64,
    energies: &[f64],
    weights: &[f64],
    assignments: &[usize],
) -> Result<Vec<f64>, usize> {
    let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
    let targets = energies
        .iter()
        .map(|energy| (-(energy - minimum) / thermal_energy).exp())
        .collect::<Vec<_>>();
    let target_total = targets.iter().sum::<f64>();
    let maximum = weights.iter().copied().fold(0.0_f64, f64::max);
    let scaled = weights
        .iter()
        .map(|weight| weight / maximum)
        .collect::<Vec<_>>();
    let total = scaled.iter().sum::<f64>();
    let mut populations = vec![0.0; energies.len()];
    for (weight, state) in scaled.iter().zip(assignments) {
        populations[*state] += weight;
    }
    // Each state's share of the total weight becomes its target share: its
    // members scale by that share times the total over its population. A
    // state far lighter than the heaviest member has a population that
    // underflows when scaled by the heaviest weight, so its members are
    // normalized within the state instead, keeping every intermediate value
    // representable.
    let mut state_maxima = vec![0.0_f64; energies.len()];
    for (weight, state) in weights.iter().zip(assignments) {
        state_maxima[*state] = state_maxima[*state].max(*weight);
    }
    let mut within = vec![0.0; energies.len()];
    for (weight, state) in weights.iter().zip(assignments) {
        within[*state] += weight / state_maxima[*state];
    }
    let shares = targets
        .iter()
        .map(|target| target / target_total)
        .collect::<Vec<_>>();
    weights
        .iter()
        .zip(assignments)
        .enumerate()
        .map(|(position, (weight, state))| {
            let factor = shares[*state] * (total / populations[*state]);
            let weight = if factor.is_finite() {
                weight * factor
            } else {
                (weight / state_maxima[*state])
                    * (shares[*state] * total / within[*state])
                    * maximum
            };
            checked_weight(weight).ok_or(position)
        })
        .collect()
}

/// `weight` when it is finite and positive.
fn checked_weight(weight: f64) -> Option<f64> {
    (weight.is_finite() && weight > 0.0).then_some(weight)
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
