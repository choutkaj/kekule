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
//! Accumulators merge, so a host can observe in parallel: it splits the
//! observations into ranges, observes each range into its own accumulator,
//! and merges the partial accumulators in range order. For a given split the
//! result is the same whichever threads observed the ranges, and it equals
//! observing everything in one accumulator up to rounding.
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

use crate::geometry::{Point3, Vector3};
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
    /// Two accumulators measure different things and cannot be merged.
    MergeMismatch(MergeMismatch),
    /// The merged total weight would leave the finite numerical range.
    MergeWeightOverflow,
    /// A merged value would leave the finite numerical range.
    MergeNumericalFailure {
        /// The atom whose value overflowed.
        atom: InstanceAtomId,
    },
}

/// What differs between two accumulators that cannot be merged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MergeMismatch {
    /// They observe structures of different topology layouts.
    Topology,
    /// They select different atoms.
    Selection,
    /// They measure different atom pairs, or the same pairs in another order.
    Pairs,
    /// They count contacts at different cutoffs.
    Cutoff,
}

impl fmt::Display for MergeMismatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Topology => "topology layouts",
            Self::Selection => "atom selections",
            Self::Pairs => "atom pairs",
            Self::Cutoff => "contact cutoffs",
        })
    }
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
            Self::MergeMismatch(mismatch) => {
                write!(f, "cannot merge accumulators with different {mismatch}")
            }
            Self::MergeWeightOverflow => {
                f.write_str("merged total weight exceeds finite numerical range")
            }
            Self::MergeNumericalFailure { atom } => {
                write!(
                    f,
                    "merged statistic exceeds finite numerical range at {atom}"
                )
            }
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

/// The total weight after merging a total of `other`, which must stay finite.
fn merged_total(total: f64, other: f64) -> Result<f64, StatisticsError> {
    let next = total + other;
    if !next.is_finite() {
        return Err(StatisticsError::MergeWeightOverflow);
    }
    Ok(next)
}

/// The share `w / (W + w)` of a weight `w` added to a total `W`, and the
/// factor `W w / (W + w)` of its squared deviation from the previous mean,
/// written so that neither product overflows.
fn shares(previous: f64, weight: f64, total: f64) -> (f64, f64) {
    (weight / total, weight * (previous / total))
}

#[derive(Debug, Clone, Copy, Default)]
struct Moments {
    mean: Vector3,
    squared_deviations: f64,
}

impl Moments {
    /// These moments combined with those of other observations, of mean
    /// `mean` and squared deviations `squared_deviations`, where `share` and
    /// `spread` are the others' [`shares`]. This is Chan et al.'s pairwise
    /// update; for one observation, without squared deviations, it is West's
    /// weighted update.
    fn combine(self, mean: Vector3, squared_deviations: f64, share: f64, spread: f64) -> Self {
        let delta = mean - self.mean;
        Self {
            mean: self.mean + delta * share,
            squared_deviations: self.squared_deviations
                + squared_deviations
                + squared_increment(delta, spread),
        }
    }

    fn is_finite(&self) -> bool {
        self.mean.is_finite() && self.squared_deviations.is_finite()
    }
}

/// `|delta|² spread`, a nonnegative quantity, so that combining moments never
/// subtracts nearby values.
fn squared_increment(delta: Vector3, spread: f64) -> f64 {
    let squared = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
    if squared.is_normal() {
        squared * spread
    } else if delta.x == 0.0 && delta.y == 0.0 && delta.z == 0.0 {
        0.0
    } else {
        // The square underflowed or overflowed; scale the length instead.
        let increment = delta.x.hypot(delta.y).hypot(delta.z) * spread.sqrt();
        increment * increment
    }
}

/// The relative margin around the squared cutoff within which a contact is
/// decided by the distance rather than its square: far more than the few
/// rounding errors by which the two can disagree.
const SQUARED_MARGIN: f64 = 1e-12;

/// Replaces each slot's moments by `combine(slot, moments)`, for a total
/// weight of `total` once every coordinate observed is at most `extent` in
/// magnitude, or `None` when a new coordinate is not finite. Returns the first
/// slot whose result is not finite and changes nothing then.
///
/// Means stay within the extent, squared deviations from them are at most
/// 12 extent², and the squared deviations sum to at most 12 extent² times the
/// total weight, up to rounding. When that is far from overflowing, no result
/// can fail, so the moments are updated in place without pending copies.
fn update_moments(
    moments: &mut Vec<Moments>,
    scratch: &mut Vec<Moments>,
    extent: Option<f64>,
    total: f64,
    combine: impl Fn(usize, Moments) -> Moments,
) -> Result<(), usize> {
    match extent {
        Some(extent) if extent * extent * total <= f64::MAX / 48.0 => {
            for (slot, moments) in moments.iter_mut().enumerate() {
                *moments = combine(slot, *moments);
            }
        }
        _ => {
            scratch.resize(moments.len(), Moments::default());
            if let Some(slot) = combine_into(moments, scratch, combine) {
                return Err(slot);
            }
            std::mem::swap(moments, scratch);
        }
    }
    Ok(())
}

/// At least the magnitude of every coordinate of `atoms` and at least
/// `extent`, or `None` when a coordinate is not finite.
///
/// Coordinates within `extent` keep it, after one comparison per coordinate
/// that does not wait on the one before. Otherwise the new extent is twice
/// the largest magnitude, so that coordinates moving a little further keep it
/// too.
fn observed_extent(extent: f64, atoms: &[TopologyAtomIndex], positions: &[Point3]) -> Option<f64> {
    let within = atoms.iter().fold(true, |within, atom| {
        let p = positions[atom.index()];
        within & (p.x.abs() <= extent) & (p.y.abs() <= extent) & (p.z.abs() <= extent)
    });
    if within {
        return Some(extent);
    }
    let mut largest = extent;
    for atom in atoms {
        let p = positions[atom.index()];
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            return None;
        }
        largest = largest.max(p.x.abs()).max(p.y.abs()).max(p.z.abs());
    }
    Some(2.0 * largest)
}

/// Writes `combine(slot, moments[slot])` to `pending[slot]` for every slot
/// and returns the first slot whose result is not finite.
fn combine_into(
    moments: &[Moments],
    pending: &mut [Moments],
    combine: impl Fn(usize, Moments) -> Moments,
) -> Option<usize> {
    for (slot, (previous, pending)) in moments.iter().zip(pending.iter_mut()).enumerate() {
        *pending = combine(slot, *previous);
        if !pending.is_finite() {
            return Some(slot);
        }
    }
    None
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
/// Accumulators over one selection [merge](Self::merge), so ranges of
/// observations can be observed in parallel into partial accumulators:
///
/// ```
/// use std::num::NonZeroUsize;
/// use std::sync::Mutex;
/// use kekule::parallel::{Executor, ScopedThreads};
/// use kekule::statistics::{RmsfAccumulator, RmsfResult};
/// use kekule::structure::Trajectory;
/// use kekule::topology::AtomSelection;
/// # fn reduce(trajectory: &Trajectory, atoms: &AtomSelection)
/// # -> Result<RmsfResult, Box<dyn std::error::Error>> {
/// // Four ranges of frames, whatever the thread count.
/// let (frames, ranges) = (trajectory.len(), 4);
/// let partials = (0..ranges).map(|_| Mutex::new(None)).collect::<Vec<_>>();
/// let executor = ScopedThreads::new(NonZeroUsize::new(2).unwrap());
/// executor.run(ranges, &|range| {
///     let partial = RmsfAccumulator::new(atoms).and_then(|mut rmsf| {
///         for index in range * frames / ranges..(range + 1) * frames / ranges {
///             rmsf.observe(index, &trajectory.get(index).unwrap(), 1.0)?;
///         }
///         Ok(rmsf)
///     });
///     *partials[range].lock().unwrap() = Some(partial);
/// });
/// // Merged in range order, the result depends on the split alone.
/// let mut rmsf = RmsfAccumulator::new(atoms)?;
/// for partial in partials {
///     rmsf.merge(&partial.into_inner().unwrap().unwrap()?)?;
/// }
/// # Ok(rmsf.finish()?)
/// # }
/// ```
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
#[derive(Debug, Clone)]
pub struct RmsfAccumulator {
    selection: AtomSelection,
    moments: Vec<Moments>,
    /// Pending moments of an update that might fail, allocated on first use.
    scratch: Vec<Moments>,
    /// At least the magnitude of every coordinate observed so far.
    extent: f64,
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
            scratch: Vec::new(),
            extent: 0.0,
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
        let positions = *positions.value();
        let atoms = self.selection.indices();
        let (share, spread) = shares(self.total_weight, weight, total);
        let first = self.observations == 0;
        let extent = observed_extent(self.extent, atoms, positions);
        let updated = update_moments(
            &mut self.moments,
            &mut self.scratch,
            extent,
            total,
            |slot, previous| {
                let p = positions[atoms[slot].index()];
                let p = Vector3::new(p.x, p.y, p.z);
                if first {
                    Moments {
                        mean: p,
                        squared_deviations: 0.0,
                    }
                } else {
                    previous.combine(p, 0.0, share, spread)
                }
            },
        );
        if let Err(slot) = updated {
            return Err(StatisticsError::NumericalFailure {
                index,
                atom: self.atom_id(slot),
            });
        }
        self.extent = extent.unwrap_or(f64::INFINITY);
        self.total_weight = total;
        self.observations += 1;
        Ok(())
    }

    /// Merges the observations of `other`, an accumulator over the same
    /// selection, as if they were observed after this accumulator's own;
    /// failure changes nothing, and `other` is left as it was.
    ///
    /// The weighted moments merge with Chan et al.'s pairwise update. The
    /// result depends on the two accumulators alone, so observations split
    /// into ranges, observed into one accumulator each and merged in range
    /// order give the same result whichever threads observed them. It equals
    /// observing every structure in one accumulator up to rounding. Merging
    /// an accumulator without observations changes nothing.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::MergeMismatch`] when `other` selects atoms
    /// of another topology layout ([`MergeMismatch::Topology`]) or other atoms
    /// ([`MergeMismatch::Selection`]),
    /// [`StatisticsError::MergeWeightOverflow`] when the total weight would
    /// leave the finite range, and [`StatisticsError::MergeNumericalFailure`]
    /// when a moment would.
    pub fn merge(&mut self, other: &Self) -> Result<(), StatisticsError> {
        if !self
            .selection
            .topology()
            .shares_layout(other.selection.topology())
        {
            return Err(StatisticsError::MergeMismatch(MergeMismatch::Topology));
        }
        if self.selection != other.selection {
            return Err(StatisticsError::MergeMismatch(MergeMismatch::Selection));
        }
        let total = merged_total(self.total_weight, other.total_weight)?;
        if other.observations == 0 {
            return Ok(());
        }
        if self.observations == 0 {
            self.moments.copy_from_slice(&other.moments);
            self.extent = other.extent;
        } else {
            let (share, spread) = shares(self.total_weight, other.total_weight, total);
            let extent = self.extent.max(other.extent);
            let updated = update_moments(
                &mut self.moments,
                &mut self.scratch,
                Some(extent),
                total,
                |slot, previous| {
                    let theirs = other.moments[slot];
                    previous.combine(theirs.mean, theirs.squared_deviations, share, spread)
                },
            );
            if let Err(slot) = updated {
                return Err(StatisticsError::MergeNumericalFailure {
                    atom: self.atom_id(slot),
                });
            }
            self.extent = extent;
        }
        self.total_weight = total;
        self.observations += other.observations;
        Ok(())
    }

    fn atom_id(&self, slot: usize) -> InstanceAtomId {
        self.selection
            .topology()
            .atom_id(self.selection.indices()[slot])
            .expect("validated selection")
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
/// Accumulators of the same pairs and cutoff [merge](Self::merge), so ranges
/// of observations can be observed in parallel, as shown for
/// [`RmsfAccumulator`].
#[derive(Debug, Clone)]
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
        let positions = *positions.value();
        let cutoff = *self.cutoff.value();
        // Squared distances decide pairs clearly inside or outside the
        // cutoff; a pair near it, or whose square is not a normal number, is
        // decided by the distance itself.
        let squared_cutoff = cutoff * cutoff;
        let (inside, outside) = if squared_cutoff.is_normal() {
            (
                squared_cutoff * (1.0 - SQUARED_MARGIN),
                squared_cutoff * (1.0 + SQUARED_MARGIN),
            )
        } else {
            (f64::NAN, f64::NAN)
        };
        for ((&(a, b), pair), hit) in self.indices.iter().zip(&self.pairs).zip(&mut self.scratch) {
            let delta = positions[a.index()] - positions[b.index()];
            let squared = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
            *hit = if squared.is_normal() && squared < inside {
                true
            } else if squared.is_normal() && squared > outside {
                false
            } else {
                let distance = delta.x.hypot(delta.y).hypot(delta.z);
                if !distance.is_finite() {
                    return Err(StatisticsError::NumericalFailure {
                        index,
                        atom: pair.0,
                    });
                }
                distance <= cutoff
            };
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

    /// Merges the observations of `other`, an accumulator of the same pairs
    /// at the same cutoff, as if they were observed after this accumulator's
    /// own; failure changes nothing, and `other` is left as it was.
    ///
    /// Weights in contact and total weights add, so observations split into
    /// ranges, observed into one accumulator each and merged in range order
    /// give the same result whichever threads observed them. It equals
    /// observing every structure in one accumulator up to rounding, and
    /// exactly when the weights are integers totalling less than 2⁵³, such as
    /// weights of one.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::MergeMismatch`] when `other` measures atoms
    /// of another topology layout ([`MergeMismatch::Topology`]), other pairs
    /// or the same pairs in another order ([`MergeMismatch::Pairs`]), or at
    /// another cutoff ([`MergeMismatch::Cutoff`]), and
    /// [`StatisticsError::MergeWeightOverflow`] when the total weight would
    /// leave the finite range.
    pub fn merge(&mut self, other: &Self) -> Result<(), StatisticsError> {
        if !self.topology.shares_layout(&other.topology) {
            return Err(StatisticsError::MergeMismatch(MergeMismatch::Topology));
        }
        if self.pairs != other.pairs {
            return Err(StatisticsError::MergeMismatch(MergeMismatch::Pairs));
        }
        if self.cutoff.value() != other.cutoff.value() {
            return Err(StatisticsError::MergeMismatch(MergeMismatch::Cutoff));
        }
        let total = merged_total(self.total_weight, other.total_weight)?;
        // Each sum stays at most its total, as for observations, so a finite
        // total keeps every sum finite.
        for (sum, theirs) in self.contact_weights.iter_mut().zip(&other.contact_weights) {
            *sum += theirs;
        }
        self.total_weight = total;
        self.observations += other.observations;
        Ok(())
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
