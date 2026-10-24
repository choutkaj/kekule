//! Weighted structural statistics: per-atom RMS fluctuation and contact
//! occupancy over weighted structures.
//!
//! An observation is one structure with an explicit, finite, positive weight.
//! Weights are relative: scaling every weight by one factor changes no
//! statistic. [`RmsfAccumulator`] and [`ContactOccupancyAccumulator`] observe
//! any borrowed structure ([`AsModelView`]: a model, an ensemble sample, a
//! trajectory frame, or a streamed frame buffer through its `as_model_view()`)
//! one at a time, with memory bounded by the selected atoms or pairs, not the
//! observation count.
//!
//! Ensembles and trajectories reduce themselves with the same accumulators.
//! [`Ensemble::rmsf`] and [`Ensemble::contact_occupancy`] observe each sample
//! with its probability, its weight divided by the total. [`Trajectory::rmsf`]
//! and [`Trajectory::contact_occupancy`] observe every frame with weight one,
//! taking the frames as equally weighted equilibrium samples; to reweight
//! frames, for example after biased sampling, observe them through an
//! accumulator with explicit weights.
//!
//! Coordinates are used as stored, ignoring cells. Align and perform any
//! periodic preprocessing explicitly before observation. Lengths are reported
//! in the canonical length unit (nm).
//!
//! # Examples
//!
//! ```
//! use kekule::geometry::Point3;
//! use kekule::structure::{Ensemble, EnsembleSample, EnsembleState, Positions};
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, NANOMETER};
//!
//! let topology = kekule::smiles::to_topology("CC")?;
//! let mut ensemble = Ensemble::new(topology.clone());
//! // The first atom sits at x = 0 with weight 3 and at x = 1 with weight 1.
//! for (x, weight) in [(0.0, 3.0), (1.0, 1.0)] {
//!     let points = [Point3::new(x, 0.0, 0.0), Point3::new(0.15, 0.0, 0.0)];
//!     let positions = Positions::new(Quantity::new(points, NANOMETER))?;
//!     ensemble.push_state(EnsembleSample::new(positions, weight)?, EnsembleState::new())?;
//! }
//!
//! // Weighted mean 0.25, weighted variance 0.75 × 0.25² + 0.25 × 0.75² = 0.1875.
//! let fluctuation = ensemble.rmsf(&AtomSelection::all(topology.clone()))?;
//! assert!((fluctuation.values().value()[0] - 0.1875_f64.sqrt()).abs() < 1e-12);
//! assert_eq!(fluctuation.values().value()[1], 0.0);
//!
//! // The pair is within 0.2 nm only in the first sample, carrying 3/4 of the weight.
//! let atoms = topology.atom_ids();
//! let cutoff = Quantity::new(0.2, NANOMETER);
//! let contacts = ensemble.contact_occupancy([(atoms[0], atoms[1])], cutoff)?;
//! assert_eq!(contacts.contacts().next(), Some(((atoms[0], atoms[1]), 0.75)));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use crate::geometry::Vector3;
use crate::structure::{AsModelView, Ensemble, Trajectory};
use crate::topology::{AtomSelection, InstanceAtomId, SelectionError, Topology, TopologyAtomIndex};
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};

/// Invalid statistics inputs or a failed observation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum StatisticsError {
    /// The atom selection is invalid; see the source.
    Selection(SelectionError),
    /// A quantity could not be converted to the canonical unit; see the
    /// source.
    Unit(UnitError),
    /// The selection names no atom or pair.
    EmptySelection,
    /// The statistic was finished before any observation.
    NoObservations,
    /// The contact cutoff is not finite and nonnegative.
    InvalidCutoff,
    /// A contact pair names the same atom twice.
    SelfPair {
        /// The index of the pair in the input.
        pair: usize,
    },
    /// A contact pair repeats an earlier unordered pair.
    DuplicatePair {
        /// The index of the repeated pair in the input.
        pair: usize,
    },
    /// An observed structure belongs to a different topology snapshot.
    TopologyMismatch {
        /// The caller's index of the observation.
        index: usize,
    },
    /// An observation weight is not finite and positive.
    InvalidWeight {
        /// The caller's index of the observation.
        index: usize,
    },
    /// The total observed weight left the finite numerical range.
    WeightOverflow {
        /// The caller's index of the observation.
        index: usize,
    },
    /// An accumulated value left the finite numerical range.
    NumericalFailure {
        /// The caller's index of the observation.
        index: usize,
        /// The atom whose value overflowed.
        atom: InstanceAtomId,
    },
}

impl fmt::Display for StatisticsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Selection(_) => f.write_str("invalid statistics selection"),
            Self::Unit(_) => f.write_str("invalid statistics unit"),
            Self::EmptySelection => f.write_str("statistic requires at least one atom or pair"),
            Self::NoObservations => f.write_str("statistic requires at least one observation"),
            Self::InvalidCutoff => f.write_str("contact cutoff must be finite and nonnegative"),
            Self::SelfPair { pair } => write!(f, "contact pair {pair} repeats the same atom"),
            Self::DuplicatePair { pair } => write!(
                f,
                "contact pair {pair} duplicates an earlier unordered pair"
            ),
            Self::TopologyMismatch { index } => write!(
                f,
                "observation {index} belongs to a different topology snapshot"
            ),
            Self::InvalidWeight { index } => {
                write!(f, "observation {index} weight must be finite and positive")
            }
            Self::WeightOverflow { index } => write!(
                f,
                "observation {index} total weight exceeds finite numerical range"
            ),
            Self::NumericalFailure { index, atom } => write!(
                f,
                "observation {index} statistic exceeds finite numerical range at {atom}"
            ),
        }
    }
}

impl std::error::Error for StatisticsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Selection(error) => Some(error),
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

/// The total weight after observing `weight`, which must be finite and
/// positive, with a total that stays finite.
fn next_total(index: usize, total: f64, weight: f64) -> Result<f64, StatisticsError> {
    if !weight.is_finite() || weight <= 0.0 {
        return Err(StatisticsError::InvalidWeight { index });
    }
    let next = total + weight;
    if !next.is_finite() {
        return Err(StatisticsError::WeightOverflow { index });
    }
    Ok(next)
}

#[derive(Debug, Clone, Copy, Default)]
struct Moments {
    mean: Vector3,
    squared_deviations: f64,
}

/// Online, weighted, per-atom RMS fluctuation about each atom's weighted mean
/// Cartesian position.
///
/// For observations `k` with weights `wₖ` and positions `xₖ` of an atom, the
/// weighted mean is `μ = Σ wₖ xₖ / W` and the RMSF is
/// `√(Σ wₖ |xₖ − μ|² / W)`, where `W = Σ wₖ`. Dividing by `W` treats the
/// weights as probabilities, without a small-sample correction; with equal
/// weights it is the population RMSF, division by the observation count. A
/// single observation yields zero.
///
/// The input selection is fixed and must be nonempty. Coordinates are used as
/// stored, ignoring cells. Align and perform any periodic preprocessing
/// explicitly before observation. Moments are updated with West's weighted
/// form of Welford's algorithm, so a large common translation does not cancel
/// the fluctuation.
///
/// State uses O(selected atoms) memory, independent of the observation count.
/// Each observation is transactional. Observation indices are diagnostic
/// labels; every successful call counts with its weight.
///
/// Reweighting trajectory frames, for example by bias factors from enhanced
/// sampling, observes each frame with its weight:
///
/// ```
/// use kekule::statistics::{RmsfAccumulator, RmsfResult};
/// use kekule::structure::Trajectory;
/// use kekule::topology::AtomSelection;
/// # fn reduce(trajectory: &Trajectory, atoms: &AtomSelection, weights: &[f64])
/// # -> Result<RmsfResult, Box<dyn std::error::Error>> {
/// let mut rmsf = RmsfAccumulator::new(atoms)?;
/// for (index, (frame, weight)) in trajectory.iter().zip(weights).enumerate() {
///     rmsf.observe(index, &frame, *weight)?;
/// }
/// let result = rmsf.finish()?;
/// for (atom, fluctuation) in result.atoms() {
///     let residue = result.selection().topology().atom(atom).and_then(|atom| atom.residue());
///     // Atom identity, optional residue, and the unit-bearing value stay associated.
///     # let _ = (residue, fluctuation);
/// }
/// # Ok(result)
/// # }
/// ```
///
/// A streamed `kekule-traj` frame buffer is observed through
/// `buffer.as_model_view()`.
#[derive(Debug)]
pub struct RmsfAccumulator {
    selection: AtomSelection,
    moments: Vec<Moments>,
    scratch: Vec<Moments>,
    total_weight: f64,
    observations: usize,
}

impl RmsfAccumulator {
    /// An empty accumulator over `selection`, which fixes the atoms and the
    /// topology layout of every observed structure.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::EmptySelection`] when `selection` has no
    /// atom.
    pub fn new(selection: &AtomSelection) -> Result<Self, StatisticsError> {
        let len = selection.indices().len();
        if len == 0 {
            return Err(StatisticsError::EmptySelection);
        }
        Ok(Self {
            selection: selection.clone(),
            moments: vec![Moments::default(); len],
            scratch: vec![Moments::default(); len],
            total_weight: 0.0,
            observations: 0,
        })
    }

    /// Observes one structure with a relative `weight`; failure preserves all
    /// accumulated statistics. `index` is the caller's label for the
    /// observation, used only in error diagnostics.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::InvalidWeight`] when `weight` is not finite
    /// and positive, [`StatisticsError::WeightOverflow`] when the total weight
    /// would leave the finite range, [`StatisticsError::TopologyMismatch`]
    /// when the structure's topology is incompatible with the selection, and
    /// [`StatisticsError::NumericalFailure`] when a moment leaves the finite
    /// range.
    pub fn observe(
        &mut self,
        index: usize,
        structure: &(impl AsModelView + ?Sized),
        weight: f64,
    ) -> Result<(), StatisticsError> {
        let structure = structure.as_model_view();
        let total = next_total(index, self.total_weight, weight)?;
        self.selection
            .ensure_compatible(&structure.shared_topology())
            .map_err(|_| StatisticsError::TopologyMismatch { index })?;
        let positions = structure.positions().values();
        // The new observation's share of the total, and the factor of its
        // squared deviation from the previous mean, `w W / (W + w)`, written
        // so that neither product overflows.
        let share = weight / total;
        let spread = weight * (self.total_weight / total);
        for ((atom, previous), pending) in self
            .selection
            .indices()
            .iter()
            .zip(&self.moments)
            .zip(&mut self.scratch)
        {
            let p = positions.value()[atom.index()];
            let p = Vector3::new(p.x, p.y, p.z);
            *pending = if self.observations == 0 {
                Moments {
                    mean: p,
                    squared_deviations: 0.0,
                }
            } else {
                // West's weighted update, with the squared increment written
                // as a nonnegative quantity to avoid subtracting nearby
                // moments.
                let delta = p - previous.mean;
                let length = delta.x.hypot(delta.y).hypot(delta.z);
                let increment = length * spread.sqrt();
                Moments {
                    mean: previous.mean + delta * share,
                    squared_deviations: previous.squared_deviations + increment * increment,
                }
            };
            if !pending.mean.is_finite() || !pending.squared_deviations.is_finite() {
                return Err(StatisticsError::NumericalFailure {
                    index,
                    atom: self
                        .selection
                        .topology()
                        .atom_id(*atom)
                        .expect("validated selection"),
                });
            }
        }
        std::mem::swap(&mut self.moments, &mut self.scratch);
        self.total_weight = total;
        self.observations += 1;
        Ok(())
    }

    /// The number of structures observed so far.
    pub fn observation_count(&self) -> usize {
        self.observations
    }

    /// The total weight observed so far.
    pub fn total_weight(&self) -> f64 {
        self.total_weight
    }

    /// Finishes with values in selection order, retaining source atom
    /// identity.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::NoObservations`] when nothing was observed.
    pub fn finish(self) -> Result<RmsfResult, StatisticsError> {
        if self.observations == 0 {
            return Err(StatisticsError::NoObservations);
        }
        let values = self
            .moments
            .into_iter()
            .map(|m| (m.squared_deviations / self.total_weight).sqrt())
            .collect();
        Ok(RmsfResult {
            selection: self.selection,
            values: Quantity::new(values, CANONICAL_LENGTH_UNIT),
            total_weight: self.total_weight,
            observations: self.observations,
        })
    }
}

/// Per-atom weighted RMSF with its original topology-bound selection,
/// observation count, and total weight. Selecting one CA per residue gives CA
/// RMSF, not an implicit residue average.
#[derive(Debug, Clone)]
pub struct RmsfResult {
    selection: AtomSelection,
    values: Quantity<Vec<f64>>,
    total_weight: f64,
    observations: usize,
}

impl RmsfResult {
    /// The selection the values follow, bound to its topology.
    pub fn selection(&self) -> &AtomSelection {
        &self.selection
    }
    /// One RMSF per selected atom in selection order, in the canonical length
    /// unit (nm).
    pub fn values(&self) -> &Quantity<Vec<f64>> {
        &self.values
    }
    /// The number of structures reduced.
    pub fn observation_count(&self) -> usize {
        self.observations
    }
    /// The total weight of the structures reduced.
    pub fn total_weight(&self) -> f64 {
        self.total_weight
    }

    /// Resolves each value to its source atom. Residues can be obtained through
    /// `result.selection().topology().atom(atom).and_then(|atom| atom.residue())`.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = (InstanceAtomId, Quantity<f64>)> + '_ {
        self.selection.atom_ids().zip(
            self.values
                .value()
                .iter()
                .map(|value| Quantity::new(*value, CANONICAL_LENGTH_UNIT)),
        )
    }
}

/// Weighted occupancy of contacts between specified atom pairs: the fraction
/// of the total observed weight in which a pair is at Cartesian distance
/// **<= cutoff**. This is a geometric contact, not a hydrogen-bond definition.
///
/// Pairs retain input order and must be nonempty, distinct, and non-self;
/// reversed duplicates are rejected. Each structure is measured anew. Cells
/// are ignored; apply periodic preprocessing first. With equal weights of one,
/// the weight in contact is the number of observations in contact. State uses
/// O(pairs) memory; failed observations preserve every sum. Cutoff
/// comparisons use canonical floating-point values without a tolerance.
#[derive(Debug)]
pub struct ContactOccupancyAccumulator {
    topology: Arc<Topology>,
    pairs: Vec<(InstanceAtomId, InstanceAtomId)>,
    indices: Vec<(TopologyAtomIndex, TopologyAtomIndex)>,
    contact_weights: Vec<f64>,
    scratch: Vec<bool>,
    cutoff: Quantity<f64>,
    total_weight: f64,
    observations: usize,
}

impl ContactOccupancyAccumulator {
    /// An empty accumulator for `pairs` of `topology` atoms in contact at
    /// distances up to `cutoff`, a length.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::Unit`] when `cutoff` is not a length,
    /// [`StatisticsError::InvalidCutoff`] when it is not finite and
    /// nonnegative, [`StatisticsError::EmptySelection`] for no pairs,
    /// [`StatisticsError::Selection`] for an atom not in `topology`,
    /// [`StatisticsError::SelfPair`] for a pair naming one atom twice, and
    /// [`StatisticsError::DuplicatePair`] for a repeated unordered pair.
    pub fn new(
        topology: &Arc<Topology>,
        pairs: impl IntoIterator<Item = (InstanceAtomId, InstanceAtomId)>,
        cutoff: Quantity<f64>,
    ) -> Result<Self, StatisticsError> {
        let cutoff = cutoff
            .into_unit(CANONICAL_LENGTH_UNIT)
            .map_err(StatisticsError::Unit)?;
        if !cutoff.value().is_finite() || *cutoff.value() < 0.0 {
            return Err(StatisticsError::InvalidCutoff);
        }
        let pairs = pairs.into_iter().collect::<Vec<_>>();
        if pairs.is_empty() {
            return Err(StatisticsError::EmptySelection);
        }
        let mut seen = BTreeSet::new();
        let mut indices = Vec::with_capacity(pairs.len());
        for (pair, &(a, b)) in pairs.iter().enumerate() {
            let index = |atom| {
                topology.atom_index(atom).ok_or(StatisticsError::Selection(
                    SelectionError::InvalidAtomId(atom),
                ))
            };
            let (a_index, b_index) = (index(a)?, index(b)?);
            if a == b {
                return Err(StatisticsError::SelfPair { pair });
            }
            if !seen.insert((a.min(b), a.max(b))) {
                return Err(StatisticsError::DuplicatePair { pair });
            }
            indices.push((a_index, b_index));
        }
        let count = pairs.len();
        Ok(Self {
            topology: Arc::clone(topology),
            pairs,
            indices,
            contact_weights: vec![0.0; count],
            scratch: vec![false; count],
            cutoff,
            total_weight: 0.0,
            observations: 0,
        })
    }

    /// Observes one structure with a relative `weight`; failure preserves
    /// every sum. `index` is the caller's label for the observation, used
    /// only in error diagnostics.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::InvalidWeight`] when `weight` is not finite
    /// and positive, [`StatisticsError::WeightOverflow`] when the total weight
    /// would leave the finite range, [`StatisticsError::TopologyMismatch`]
    /// when the structure does not share the accumulator's layout, and
    /// [`StatisticsError::NumericalFailure`] when a distance is not finite.
    pub fn observe(
        &mut self,
        index: usize,
        structure: &(impl AsModelView + ?Sized),
        weight: f64,
    ) -> Result<(), StatisticsError> {
        let structure = structure.as_model_view();
        let total = next_total(index, self.total_weight, weight)?;
        if !self.topology.shares_layout(structure.topology()) {
            return Err(StatisticsError::TopologyMismatch { index });
        }
        let positions = structure.positions().values();
        for ((&(a, b), pair), hit) in self.indices.iter().zip(&self.pairs).zip(&mut self.scratch) {
            let delta = positions.value()[a.index()] - positions.value()[b.index()];
            let distance = delta.x.hypot(delta.y).hypot(delta.z);
            if !distance.is_finite() {
                return Err(StatisticsError::NumericalFailure {
                    index,
                    atom: pair.0,
                });
            }
            *hit = distance <= *self.cutoff.value();
        }
        // A pair's weight in contact sums a subsequence of the same weights
        // as the total, so with monotone rounding it never exceeds it.
        for (sum, hit) in self.contact_weights.iter_mut().zip(&self.scratch) {
            if *hit {
                *sum += weight;
            }
        }
        self.total_weight = total;
        self.observations += 1;
        Ok(())
    }

    /// The number of structures observed so far.
    pub fn observation_count(&self) -> usize {
        self.observations
    }

    /// The total weight observed so far.
    pub fn total_weight(&self) -> f64 {
        self.total_weight
    }

    /// Finishes with sums in input pair order.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::NoObservations`] when nothing was observed.
    pub fn finish(self) -> Result<ContactOccupancyResult, StatisticsError> {
        if self.observations == 0 {
            return Err(StatisticsError::NoObservations);
        }
        Ok(ContactOccupancyResult {
            topology: self.topology,
            pairs: self.pairs,
            contact_weights: self.contact_weights,
            cutoff: self.cutoff,
            total_weight: self.total_weight,
            observations: self.observations,
        })
    }
}

/// Contact weights and occupancies associated with the original atom pairs.
#[derive(Debug, Clone)]
pub struct ContactOccupancyResult {
    topology: Arc<Topology>,
    pairs: Vec<(InstanceAtomId, InstanceAtomId)>,
    contact_weights: Vec<f64>,
    cutoff: Quantity<f64>,
    total_weight: f64,
    observations: usize,
}

impl ContactOccupancyResult {
    /// The topology the pairs belong to.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }
    /// The atom pairs in input order.
    pub fn pairs(&self) -> &[(InstanceAtomId, InstanceAtomId)] {
        &self.pairs
    }
    /// The total weight of the structures in which each pair was in contact,
    /// in pair order. With weights of one it is the number of structures.
    pub fn contact_weights(&self) -> &[f64] {
        &self.contact_weights
    }
    /// The number of structures reduced.
    pub fn observation_count(&self) -> usize {
        self.observations
    }
    /// The total weight of the structures reduced.
    pub fn total_weight(&self) -> f64 {
        self.total_weight
    }
    /// The contact cutoff in the canonical length unit (nm).
    pub fn cutoff(&self) -> Quantity<f64> {
        self.cutoff
    }

    /// Returns (source atom pair, occupancy): the fraction of the total weight
    /// in which the pair was in contact, between zero and one.
    pub fn contacts(
        &self,
    ) -> impl ExactSizeIterator<Item = ((InstanceAtomId, InstanceAtomId), f64)> + '_ {
        self.pairs.iter().copied().zip(
            self.contact_weights
                .iter()
                .map(|weight| weight / self.total_weight),
        )
    }
}

impl Ensemble {
    /// Per-atom weighted RMSF about each atom's weighted mean stored position
    /// over every sample, as defined by [`RmsfAccumulator`], each sample
    /// observed with its probability. Align and perform periodic
    /// preprocessing beforehand.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::Selection`] when `selection` is incompatible
    /// with the ensemble's topology, [`StatisticsError::NoObservations`] for
    /// an empty ensemble, [`StatisticsError::InvalidWeight`] at a sample whose
    /// probability underflows to zero, and the errors of
    /// [`RmsfAccumulator::new`] and [`RmsfAccumulator::observe`].
    pub fn rmsf(&self, selection: &AtomSelection) -> Result<RmsfResult, StatisticsError> {
        selection
            .ensure_compatible(&self.shared_topology())
            .map_err(StatisticsError::Selection)?;
        let mut accumulator = RmsfAccumulator::new(selection)?;
        for (index, (sample, probability)) in self.iter().zip(self.observed_weights()).enumerate() {
            accumulator.observe(index, &sample, probability)?;
        }
        accumulator.finish()
    }

    /// Weighted Cartesian contact occupancies over every sample, as defined
    /// by [`ContactOccupancyAccumulator`], each sample observed with its
    /// probability.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::NoObservations`] for an empty ensemble,
    /// [`StatisticsError::InvalidWeight`] at a sample whose probability
    /// underflows to zero, and the errors of
    /// [`ContactOccupancyAccumulator::new`] and
    /// [`ContactOccupancyAccumulator::observe`].
    pub fn contact_occupancy(
        &self,
        pairs: impl IntoIterator<Item = (InstanceAtomId, InstanceAtomId)>,
        cutoff: Quantity<f64>,
    ) -> Result<ContactOccupancyResult, StatisticsError> {
        let mut accumulator =
            ContactOccupancyAccumulator::new(&self.shared_topology(), pairs, cutoff)?;
        for (index, (sample, probability)) in self.iter().zip(self.observed_weights()).enumerate() {
            accumulator.observe(index, &sample, probability)?;
        }
        accumulator.finish()
    }

    /// Each sample's probability, computed from weights scaled by the largest
    /// so that the total cannot overflow.
    fn observed_weights(&self) -> Vec<f64> {
        let scaled = self.scaled_weights();
        let total = scaled.iter().sum::<f64>();
        scaled.into_iter().map(|weight| weight / total).collect()
    }
}

impl Trajectory {
    /// Per-atom population RMSF about mean stored positions over every frame,
    /// as defined by [`RmsfAccumulator`] with every frame observed with weight
    /// one: the frames are taken as equally weighted equilibrium samples.
    /// Align and perform periodic preprocessing beforehand. To weight frames
    /// differently, observe them through an [`RmsfAccumulator`].
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::Selection`] when `selection` is incompatible
    /// with the trajectory's topology, [`StatisticsError::NoObservations`] for
    /// an empty trajectory, and the errors of [`RmsfAccumulator::new`] and
    /// [`RmsfAccumulator::observe`].
    pub fn rmsf(&self, selection: &AtomSelection) -> Result<RmsfResult, StatisticsError> {
        selection
            .ensure_compatible(&self.shared_topology())
            .map_err(StatisticsError::Selection)?;
        let mut accumulator = RmsfAccumulator::new(selection)?;
        for (index, frame) in self.iter().enumerate() {
            accumulator.observe(index, &frame, 1.0)?;
        }
        accumulator.finish()
    }

    /// Cartesian contact occupancies over every frame, as defined by
    /// [`ContactOccupancyAccumulator`] with every frame observed with weight
    /// one, so each contact weight is a frame count. To weight frames
    /// differently, observe them through a [`ContactOccupancyAccumulator`].
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::NoObservations`] for an empty trajectory and
    /// the errors of [`ContactOccupancyAccumulator::new`] and
    /// [`ContactOccupancyAccumulator::observe`].
    pub fn contact_occupancy(
        &self,
        pairs: impl IntoIterator<Item = (InstanceAtomId, InstanceAtomId)>,
        cutoff: Quantity<f64>,
    ) -> Result<ContactOccupancyResult, StatisticsError> {
        let mut accumulator =
            ContactOccupancyAccumulator::new(&self.shared_topology(), pairs, cutoff)?;
        for (index, frame) in self.iter().enumerate() {
            accumulator.observe(index, &frame, 1.0)?;
        }
        accumulator.finish()
    }
}
