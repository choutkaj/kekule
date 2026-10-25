//! Weighted structural statistics: per-atom RMS fluctuation and contact
//! occupancy of atom pairs or atom groups over weighted structures.
//!
//! An observation is one structure with an explicit, finite, positive weight.
//! Weights are relative: scaling every weight by one factor changes no
//! statistic. [`RmsfAccumulator`], [`ContactOccupancyAccumulator`], and
//! [`GroupContactOccupancyAccumulator`] observe any borrowed structure
//! ([`AsModelView`]: a model, an ensemble sample, a trajectory frame, or a
//! streamed frame buffer through its `as_model_view()`) one at a time, with
//! memory bounded by the selected atoms, pairs, or groups, not the
//! observation count.
//!
//! Ensembles and trajectories reduce themselves with the same accumulators.
//! [`Ensemble::rmsf`], [`Ensemble::contact_occupancy`], and
//! [`Ensemble::group_contact_occupancy`] observe each sample with its
//! probability, its weight divided by the total. [`Trajectory::rmsf`],
//! [`Trajectory::contact_occupancy`], and
//! [`Trajectory::group_contact_occupancy`] observe every frame with weight one,
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
    /// A contact group has no atom.
    EmptyGroup {
        /// The index of the group in the input.
        group: usize,
    },
    /// A pair names a group past the last.
    UnknownGroup {
        /// The index of the pair in the input.
        pair: usize,
    },
    /// A pair names two groups that share an atom.
    OverlappingGroups {
        /// The index of the pair in the input.
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
            Self::EmptyGroup { group } => write!(f, "contact group {group} has no atom"),
            Self::UnknownGroup { pair } => {
                write!(f, "contact pair {pair} names a group past the last")
            }
            Self::OverlappingGroups { pair } => {
                write!(f, "contact pair {pair} names groups that share an atom")
            }
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

/// Group pairs whose atom pair counts are at most this are measured pair by
/// pair; larger ones search a grid of the larger group.
const DIRECT_PAIRS: usize = 256;

/// Weighted occupancy of contacts between specified pairs of atom groups: the
/// fraction of the total observed weight in which a pair of groups is in
/// contact, meaning that some atom of one group is at Cartesian distance
/// **<= cutoff** from some atom of the other. A residue is in contact with a
/// ligand in an observation when any of its atoms is near any ligand atom,
/// which the atom pairs' occupancies cannot give, since their contacts overlap
/// across observations.
///
/// Groups are atom selections of the topology, such as a pocket's residues and
/// a ligand, or the residues on both sides of an interface. Each must be
/// nonempty; they may overlap, and groups no pair names are ignored. Pairs name
/// groups by their input index, retain input order, and must be nonempty,
/// distinct, and between disjoint groups; reversed duplicates are rejected.
/// A pair of groups of one atom each has the occupancy that
/// [`ContactOccupancyAccumulator`] gives the atom pair.
///
/// Each structure is measured anew, reading only the grouped atoms: a pair is
/// rejected when the groups' bounding boxes lie farther apart than the cutoff,
/// and a pair of large groups searches a grid of the larger group, built at
/// most once per observation, so the cost does not grow with the structure's
/// other atoms or with every atom pair of the groups. Cells are ignored; apply
/// periodic preprocessing first. With equal weights of one, the weight in
/// contact is the number of observations in contact. State uses O(pairs)
/// sums and O(grouped atoms) scratch, independent of the observation count;
/// failed observations preserve every sum. Distances are compared with the
/// cutoff as [`ContactOccupancyAccumulator`] compares them, without a
/// tolerance.
///
/// ```
/// use kekule::geometry::Point3;
/// use kekule::statistics::GroupContactOccupancyAccumulator;
/// use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
/// use kekule::topology::AtomSelection;
/// use kekule::units::{Quantity, NANOMETER};
///
/// let topology = kekule::smiles::to_topology("CCCC")?;
/// let atoms = topology.atom_ids();
/// let group = |members: &[usize]| {
///     AtomSelection::from_atoms(topology.clone(), members.iter().map(|&atom| atoms[atom]))
/// };
/// // Atoms 0 and 1 against atoms 2 and 3. In the first frame atom 1 touches
/// // atom 2; in the second atom 0 touches atom 3; in the third none touches.
/// let frames = [[0.0, 1.0, 1.1, 3.0], [0.0, 1.0, 2.0, 0.1], [0.0, 1.0, 2.0, 3.0]];
/// let trajectory = Trajectory::from_items(
///     topology.clone(),
///     frames.iter().map(|xs| {
///         let points = xs.map(|x| Point3::new(x, 0.0, 0.0));
///         TrajectoryFrame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap())
///     }),
/// )?;
/// let cutoff = Quantity::new(0.2, NANOMETER);
///
/// let mut contacts = GroupContactOccupancyAccumulator::new(
///     &topology,
///     [group(&[0, 1])?, group(&[2, 3])?],
///     [(0, 1)],
///     cutoff,
/// )?;
/// for (index, frame) in trajectory.iter().enumerate() {
///     contacts.observe(index, &frame, 1.0)?;
/// }
/// let result = contacts.finish()?;
/// assert_eq!(result.contact_weights(), [2.0]);
/// assert_eq!(result.contacts().next(), Some(((0, 1), 2.0 / 3.0)));
///
/// // Each atom pair is in contact a third of the time; whether in the same
/// // frames or not, which decides the groups' occupancy, they do not say.
/// let pairs = trajectory.contact_occupancy([(atoms[1], atoms[2]), (atoms[0], atoms[3])], cutoff)?;
/// assert!(pairs.contacts().all(|(_, occupancy)| occupancy == 1.0 / 3.0));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct GroupContactOccupancyAccumulator {
    topology: Arc<Topology>,
    groups: Vec<AtomSelection>,
    pairs: Vec<(usize, usize)>,
    contact_weights: Vec<f64>,
    cutoff: Quantity<f64>,
    total_weight: f64,
    observations: usize,
    scratch: GroupScratch,
}

/// Per-observation working memory, sized by the groups that pairs name.
#[derive(Debug)]
struct GroupScratch {
    /// Whether some pair names each group.
    paired: Vec<bool>,
    /// The coordinates of each group's atoms in selection order, group after
    /// group: group `g` is `points[starts[g]..starts[g + 1]]`.
    points: Vec<Point3>,
    starts: Vec<usize>,
    /// Each group's bounding box in this observation.
    bounds: Vec<(Point3, Point3)>,
    /// Each group's grid, valid in this observation when `built` says so.
    grids: Vec<GroupGrid>,
    built: Vec<bool>,
    hits: Vec<bool>,
}

type Cell = (i64, i64, i64);

/// One group's points binned into cubic cells at least as wide as the cutoff,
/// sorted by cell, so the points within the cutoff of a point lie in the 27
/// cells around its own.
#[derive(Debug, Default)]
struct GroupGrid {
    low: Point3,
    scale: f64,
    /// The highest cell on each axis.
    high: Cell,
    entries: Vec<(Cell, Point3)>,
}

impl GroupGrid {
    /// Rebuilds the grid over `points`, whose bounding box is `bounds`.
    fn build(&mut self, points: &[Point3], (low, high): (Point3, Point3), cutoff: f64) {
        // Cells slightly wider than the cutoff keep every contact within the
        // adjacent cells despite rounding, and at least a millionth of the
        // box, so cell coordinates stay small enough to round alike. Any
        // positive edge finds coincident points.
        let extent = (high.x - low.x).max(high.y - low.y).max(high.z - low.z);
        let edge = (cutoff * (1.0 + 1e-6)).max(extent * 1e-6);
        let edge = if edge > 0.0 { edge } else { 1.0 };
        self.low = low;
        self.scale = edge.recip();
        let mut entries = std::mem::take(&mut self.entries);
        entries.clear();
        entries.extend(points.iter().map(|point| (self.cell(point), *point)));
        self.entries = entries;
        self.entries.sort_unstable_by_key(|(cell, _)| *cell);
        self.high = self.cell(&high);
    }

    /// The cell of `point`; one monotone computation places every point, so
    /// points nearer than an edge along an axis lie in adjacent cells.
    fn cell(&self, point: &Point3) -> Cell {
        let axis = |value: f64, origin: f64| ((value - origin) * self.scale).floor() as i64;
        (
            axis(point.x, self.low.x),
            axis(point.y, self.low.y),
            axis(point.z, self.low.z),
        )
    }

    /// Whether any of `candidates` is within `cutoff` of a grid point.
    fn any_within(&self, candidates: &[Point3], cutoff: f64) -> bool {
        candidates.iter().any(|candidate| {
            let (x, y, z) = self.cell(candidate);
            if x < -1
                || y < -1
                || z < -1
                || x > self.high.0.saturating_add(1)
                || y > self.high.1.saturating_add(1)
                || z > self.high.2.saturating_add(1)
            {
                return false;
            }
            // Cells sort x-major, so the three cells along z around each
            // (x, y) column are one run.
            (x - 1..=x + 1).any(|x| {
                (y - 1..=y + 1).any(|y| {
                    let first = self
                        .entries
                        .partition_point(|(cell, _)| *cell < (x, y, z - 1));
                    let end = self
                        .entries
                        .partition_point(|(cell, _)| *cell <= (x, y, z + 1));
                    self.entries[first..end]
                        .iter()
                        .any(|(_, point)| within(candidate, point, cutoff))
                })
            })
        })
    }
}

/// The contact test of [`ContactOccupancyAccumulator`].
fn within(a: &Point3, b: &Point3, cutoff: f64) -> bool {
    let delta = *a - *b;
    delta.x.hypot(delta.y).hypot(delta.z) <= cutoff
}

/// The bounding box of nonempty `points`.
fn bounds(points: &[Point3]) -> (Point3, Point3) {
    points
        .iter()
        .fold((points[0], points[0]), |(low, high), point| {
            (
                Point3::new(low.x.min(point.x), low.y.min(point.y), low.z.min(point.z)),
                Point3::new(
                    high.x.max(point.x),
                    high.y.max(point.y),
                    high.z.max(point.z),
                ),
            )
        })
}

/// The length of the diagonal of the bounding box of two bounding boxes.
fn diagonal((low, high): (Point3, Point3), (other_low, other_high): (Point3, Point3)) -> f64 {
    let side = |low: f64, high: f64, other_low: f64, other_high: f64| {
        high.max(other_high) - low.min(other_low)
    };
    side(low.x, high.x, other_low.x, other_high.x)
        .hypot(side(low.y, high.y, other_low.y, other_high.y))
        .hypot(side(low.z, high.z, other_low.z, other_high.z))
}

/// Whether two bounding boxes lie farther apart than `cutoff` along an axis,
/// so that no two of their points are within it: a coordinate difference
/// rounds to at least the boxes' gap, and a distance to at least each
/// coordinate difference.
fn apart(
    (low, high): (Point3, Point3),
    (other_low, other_high): (Point3, Point3),
    cutoff: f64,
) -> bool {
    let gap = |low: f64, high: f64, other_low: f64, other_high: f64| {
        (other_low - high).max(low - other_high)
    };
    gap(low.x, high.x, other_low.x, other_high.x) > cutoff
        || gap(low.y, high.y, other_low.y, other_high.y) > cutoff
        || gap(low.z, high.z, other_low.z, other_high.z) > cutoff
}

impl GroupContactOccupancyAccumulator {
    /// An empty accumulator for `pairs` of `groups`, atom selections of
    /// `topology` named by their index in `groups`, in contact at distances up
    /// to `cutoff`, a length.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::Unit`] when `cutoff` is not a length,
    /// [`StatisticsError::InvalidCutoff`] when it is not finite and
    /// nonnegative, [`StatisticsError::Selection`] for a group whose topology
    /// does not share `topology`'s layout, [`StatisticsError::EmptyGroup`] for
    /// a group with no atom, [`StatisticsError::EmptySelection`] for no pairs,
    /// [`StatisticsError::UnknownGroup`] for a pair naming a group past the
    /// last, [`StatisticsError::SelfPair`] for a pair naming one group twice,
    /// [`StatisticsError::DuplicatePair`] for a repeated unordered pair, and
    /// [`StatisticsError::OverlappingGroups`] for a pair of groups sharing an
    /// atom.
    pub fn new(
        topology: &Arc<Topology>,
        groups: impl IntoIterator<Item = AtomSelection>,
        pairs: impl IntoIterator<Item = (usize, usize)>,
        cutoff: Quantity<f64>,
    ) -> Result<Self, StatisticsError> {
        let cutoff = cutoff
            .into_unit(CANONICAL_LENGTH_UNIT)
            .map_err(StatisticsError::Unit)?;
        if !cutoff.value().is_finite() || *cutoff.value() < 0.0 {
            return Err(StatisticsError::InvalidCutoff);
        }
        let groups = groups.into_iter().collect::<Vec<_>>();
        for (group, selection) in groups.iter().enumerate() {
            selection
                .ensure_compatible(topology)
                .map_err(StatisticsError::Selection)?;
            if selection.is_empty() {
                return Err(StatisticsError::EmptyGroup { group });
            }
        }
        let pairs = pairs.into_iter().collect::<Vec<_>>();
        if pairs.is_empty() {
            return Err(StatisticsError::EmptySelection);
        }
        let mut seen = BTreeSet::new();
        let mut paired = vec![false; groups.len()];
        for (pair, &(a, b)) in pairs.iter().enumerate() {
            if a >= groups.len() || b >= groups.len() {
                return Err(StatisticsError::UnknownGroup { pair });
            }
            if a == b {
                return Err(StatisticsError::SelfPair { pair });
            }
            if !seen.insert((a.min(b), a.max(b))) {
                return Err(StatisticsError::DuplicatePair { pair });
            }
            if !groups[a]
                .is_disjoint(&groups[b])
                .map_err(StatisticsError::Selection)?
            {
                return Err(StatisticsError::OverlappingGroups { pair });
            }
            paired[a] = true;
            paired[b] = true;
        }
        let mut starts = Vec::with_capacity(groups.len() + 1);
        let mut start = 0;
        starts.push(start);
        for (selection, paired) in groups.iter().zip(&paired) {
            if *paired {
                start += selection.len();
            }
            starts.push(start);
        }
        let count = pairs.len();
        let scratch = GroupScratch {
            paired,
            points: Vec::with_capacity(start),
            starts,
            bounds: vec![(Point3::default(), Point3::default()); groups.len()],
            grids: std::iter::repeat_with(GroupGrid::default)
                .take(groups.len())
                .collect(),
            built: vec![false; groups.len()],
            hits: vec![false; count],
        };
        Ok(Self {
            topology: Arc::clone(topology),
            groups,
            pairs,
            contact_weights: vec![0.0; count],
            cutoff,
            total_weight: 0.0,
            observations: 0,
            scratch,
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
    /// [`StatisticsError::NumericalFailure`], naming the first atom of the
    /// pair's first group, when distances between a pair's groups may leave
    /// the finite range: when the diagonal of their joint bounding box is not
    /// finite. For groups of one atom each, that diagonal is their distance.
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
        let positions = positions.value();
        let cutoff = *self.cutoff.value();
        let scratch = &mut self.scratch;
        scratch.points.clear();
        for (group, selection) in self.groups.iter().enumerate() {
            if !scratch.paired[group] {
                continue;
            }
            scratch.points.extend(
                selection
                    .indices()
                    .iter()
                    .map(|atom| positions[atom.index()]),
            );
            let points = &scratch.points[scratch.starts[group]..];
            scratch.bounds[group] = bounds(points);
            scratch.built[group] = false;
        }
        for (&(a, b), hit) in self.pairs.iter().zip(&mut scratch.hits) {
            // Every distance between the groups is at most the diagonal of
            // their joint bounding box, as rounded.
            if !diagonal(scratch.bounds[a], scratch.bounds[b]).is_finite() {
                return Err(StatisticsError::NumericalFailure {
                    index,
                    atom: self
                        .topology
                        .atom_id(self.groups[a].indices()[0])
                        .expect("validated selection"),
                });
            }
            *hit = if apart(scratch.bounds[a], scratch.bounds[b], cutoff) {
                false
            } else {
                let group = |g: usize| &scratch.points[scratch.starts[g]..scratch.starts[g + 1]];
                // The smaller group probes the larger.
                let (small, large) = if group(a).len() <= group(b).len() {
                    (a, b)
                } else {
                    (b, a)
                };
                if group(small).len().saturating_mul(group(large).len()) <= DIRECT_PAIRS {
                    group(small).iter().any(|point| {
                        group(large)
                            .iter()
                            .any(|other| within(point, other, cutoff))
                    })
                } else {
                    if !scratch.built[large] {
                        let points =
                            &scratch.points[scratch.starts[large]..scratch.starts[large + 1]];
                        scratch.grids[large].build(points, scratch.bounds[large], cutoff);
                        scratch.built[large] = true;
                    }
                    let points = &scratch.points[scratch.starts[small]..scratch.starts[small + 1]];
                    scratch.grids[large].any_within(points, cutoff)
                }
            };
        }
        // A pair's weight in contact sums a subsequence of the same weights
        // as the total, so with monotone rounding it never exceeds it.
        for (sum, hit) in self.contact_weights.iter_mut().zip(&scratch.hits) {
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
    pub fn finish(self) -> Result<GroupContactOccupancyResult, StatisticsError> {
        if self.observations == 0 {
            return Err(StatisticsError::NoObservations);
        }
        Ok(GroupContactOccupancyResult {
            topology: self.topology,
            groups: self.groups,
            pairs: self.pairs,
            contact_weights: self.contact_weights,
            cutoff: self.cutoff,
            total_weight: self.total_weight,
            observations: self.observations,
        })
    }
}

/// Contact weights and occupancies associated with the original pairs of atom
/// groups.
#[derive(Debug, Clone)]
pub struct GroupContactOccupancyResult {
    topology: Arc<Topology>,
    groups: Vec<AtomSelection>,
    pairs: Vec<(usize, usize)>,
    contact_weights: Vec<f64>,
    cutoff: Quantity<f64>,
    total_weight: f64,
    observations: usize,
}

impl GroupContactOccupancyResult {
    /// The topology the groups belong to.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }
    /// The atom groups in input order, which pairs name by index.
    pub fn groups(&self) -> &[AtomSelection] {
        &self.groups
    }
    /// The pairs of group indices in input order.
    pub fn pairs(&self) -> &[(usize, usize)] {
        &self.pairs
    }
    /// The total weight of the structures in which each pair of groups was
    /// in contact, in pair order. With weights of one it is the number of
    /// structures.
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

    /// Returns (pair of group indices, occupancy): the fraction of the total
    /// weight in which the groups were in contact, between zero and one.
    pub fn contacts(&self) -> impl ExactSizeIterator<Item = ((usize, usize), f64)> + '_ {
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

    /// Weighted Cartesian contact occupancies of pairs of atom groups over
    /// every sample, as defined by [`GroupContactOccupancyAccumulator`], each
    /// sample observed with its probability.
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::NoObservations`] for an empty ensemble,
    /// [`StatisticsError::InvalidWeight`] at a sample whose probability
    /// underflows to zero, and the errors of
    /// [`GroupContactOccupancyAccumulator::new`] and
    /// [`GroupContactOccupancyAccumulator::observe`].
    pub fn group_contact_occupancy(
        &self,
        groups: impl IntoIterator<Item = AtomSelection>,
        pairs: impl IntoIterator<Item = (usize, usize)>,
        cutoff: Quantity<f64>,
    ) -> Result<GroupContactOccupancyResult, StatisticsError> {
        let mut accumulator =
            GroupContactOccupancyAccumulator::new(&self.shared_topology(), groups, pairs, cutoff)?;
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

    /// Cartesian contact occupancies of pairs of atom groups over every
    /// frame, as defined by [`GroupContactOccupancyAccumulator`] with every
    /// frame observed with weight one, so each contact weight is a frame
    /// count. To weight frames differently, observe them through a
    /// [`GroupContactOccupancyAccumulator`].
    ///
    /// # Errors
    ///
    /// Returns [`StatisticsError::NoObservations`] for an empty trajectory and
    /// the errors of [`GroupContactOccupancyAccumulator::new`] and
    /// [`GroupContactOccupancyAccumulator::observe`].
    pub fn group_contact_occupancy(
        &self,
        groups: impl IntoIterator<Item = AtomSelection>,
        pairs: impl IntoIterator<Item = (usize, usize)>,
        cutoff: Quantity<f64>,
    ) -> Result<GroupContactOccupancyResult, StatisticsError> {
        let mut accumulator =
            GroupContactOccupancyAccumulator::new(&self.shared_topology(), groups, pairs, cutoff)?;
        for (index, frame) in self.iter().enumerate() {
            accumulator.observe(index, &frame, 1.0)?;
        }
        accumulator.finish()
    }
}
