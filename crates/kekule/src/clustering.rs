//! Partitions of the members of a [`SampleSet`]'s distribution, or of an
//! [`Ensemble`]'s samples, into states of similar structures.
//!
//! A sample set and an ensemble cluster through one implementation, so a
//! set converted from an ensemble ([`SampleSet::from_ensemble`]) clusters
//! exactly as the ensemble does. A set clusters the members of its
//! distribution, leaving the other samples out, or every sample when it has
//! no distribution. Its results are indexed by member position, the `i`th
//! member being the `i`th sample
//! [`Distribution::members`](crate::structure::Distribution::members) lists
//! (sample `i` without a distribution): the [`DistanceMatrix`] of
//! [`SampleSet::pairwise_fitted_rmsd`] and each partition, whose
//! assignments and representatives are member positions. A partition
//! installs with
//! [`DistributionMut::regroup`](crate::structure::DistributionMut::regroup),
//! which changes no weight, or as a new distribution of every sample through
//! [`SampleSet::set_distribution`]. Errors name sample indices.
//!
//! The examples below use an ensemble; the same calls work on a sample set.
//!
//! Clustering never changes the ensemble. It returns an
//! [`EnsemblePartition`] that assigns every sample to one state and names one
//! representative sample per state. [`Ensemble::regroup`] applies it without
//! changing any weight, so populations follow the samples, and
//! [`Ensemble::collapse`] then keeps one sample per state, carrying the
//! state's population.
//!
//! [`Ensemble::cluster_hierarchical`] cuts the [`Dendrogram`] of a
//! [`DistanceMatrix`], such as [`Ensemble::pairwise_fitted_rmsd`], at a height or
//! into a number of clusters, and represents each cluster by its weighted
//! medoid. Complete linkage cut at height `h` bounds every cluster's
//! diameter by `h`.
//!
//! [`Ensemble::pairwise_fitted_rmsd_with`] and
//! [`Ensemble::cluster_by_radius_with`] take a [`ClusteringRun`]: a
//! [`parallel::Executor`](crate::parallel::Executor) that spreads the fits
//! over threads, a [`Cancellation`](crate::diagnostics::Cancellation), and a
//! [`ClusteringProgress`] callback. Fits run in contiguous ranges whose
//! results combine in a fixed order, so every executor and thread count
//! returns exactly the single-threaded result, including which error is
//! reported.
//!
//! [`Ensemble::cluster_k_medoids`] finds a given number of clusters directly,
//! choosing medoid samples that minimize the weight-summed distance from
//! every sample to its nearest medoid, so heavy samples attract medoids.
//!
//! [`Ensemble::cluster_by_radius`] is leader clustering on fitted RMSD. It
//! needs no matrix, so it scales to larger ensembles. With
//! a [`SymmetricSelection`](crate::alignment::SymmetricSelection) and a small
//! radius it merges rotamers, copies of one conformer that differ only by
//! exchanged equivalent atoms. Weight the conformers by energy first; the
//! merged state keeps its representative's energy and gains the population
//! of every rotamer, its degeneracy:
//!
//! ```
//! use std::error::Error;
//!
//! use kekule::alignment::SymmetricSelection;
//! use kekule::geometry::Point3;
//! use kekule::structure::{Ensemble, EnsembleSample, EnsembleState, Positions};
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, ANGSTROM, KELVIN, KILOJOULE_PER_MOLE};
//!
//! // Acetate, then a rotamer with its oxygens exchanged, then a distinct
//! // conformer 5 kJ/mol higher, each its own state.
//! let topology = kekule::smiles::to_topology("CC(=O)[O-]")?;
//! let mut ensemble = Ensemble::new(topology.clone());
//! for (o2, o3, energy) in [(1.1, -1.1, 0.0), (-1.1, 1.1, 0.0), (1.1, -0.3, 5.0)] {
//!     let points = [
//!         Point3::new(0.0, 0.0, 0.0),
//!         Point3::new(1.5, 0.0, 0.0),
//!         Point3::new(2.2, o2, 0.5),
//!         Point3::new(2.2, o3, 0.5),
//!     ];
//!     let positions = Positions::new(Quantity::new(points, ANGSTROM))?;
//!     let state = EnsembleState::new().with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))?;
//!     ensemble.push_state(EnsembleSample::new(positions, 1.0)?, state)?;
//! }
//! ensemble.assign_boltzmann_populations(Quantity::new(298.15, KELVIN))?;
//!
//! let atoms = SymmetricSelection::new(&AtomSelection::all(topology.clone()))?;
//! let radius = Quantity::new(0.125, ANGSTROM);
//! let rotamers = ensemble.cluster_by_radius(radius, &atoms, Default::default())?;
//! assert_eq!(rotamers.assignments(), [0, 0, 1]);
//!
//! ensemble.regroup(rotamers)?;
//! assert_eq!(ensemble.state_count(), 2);
//! assert!(ensemble.state(0).unwrap().energy().is_some());
//! let conformers = ensemble.collapse()?;
//! assert_eq!(conformers.len(), 2);
//! assert_eq!(conformers.get(0).unwrap().weight(), ensemble.state(0).unwrap().population());
//! # Ok::<(), Box<dyn Error>>(())
//! ```

use std::fmt;

mod hierarchical;
mod medoids;
mod run;
pub use hierarchical::{Cut, Dendrogram, DistanceMatrix, Linkage, Merge};
pub use run::{ClusteringProgress, ClusteringRun};

use crate::alignment::{
    best_kabsch, check_fit_count, AlignmentError, AlignmentOptions, FitAtoms, NormalizedWeights,
};
use crate::diagnostics::Cancelled;
use crate::parallel::map_ranges;
use crate::structure::{
    Ensemble, EnsembleError, EnsemblePartition, Members, ModelView, Requirement, SampleSet,
    SampleSetError, StructureError,
};
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};
use run::{FirstFailure, Stop, Tracker};

/// The order in which leader clustering visits members.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum ClusterOrder {
    /// Descending sample weight, so the heaviest sample of a cluster leads
    /// it.
    #[default]
    Weight,
    /// Ascending energy of each sample's state, so a sample of the
    /// lowest-energy state leads its cluster; samples of one state keep
    /// sample order.
    Energy,
}

/// Options for [`Ensemble::cluster_by_radius`].
#[derive(Debug, Clone, Copy, Default)]
pub struct ClusteringOptions<'a> {
    alignment: AlignmentOptions<'a>,
    order: ClusterOrder,
}

impl<'a> ClusteringOptions<'a> {
    /// Fit weighting and periodic policy of every pairwise fit.
    pub const fn alignment(&self) -> AlignmentOptions<'a> {
        self.alignment
    }

    /// Sets [`Self::alignment`].
    #[must_use]
    pub const fn with_alignment(mut self, alignment: AlignmentOptions<'a>) -> Self {
        self.alignment = alignment;
        self
    }

    /// The order in which members are visited.
    pub const fn order(&self) -> ClusterOrder {
        self.order
    }

    /// Sets [`Self::order`].
    #[must_use]
    pub const fn with_order(mut self, order: ClusterOrder) -> Self {
        self.order = order;
        self
    }
}

/// The structures a clustering groups, in member order: an ensemble's
/// samples, or the members of a sample set's distribution (every sample of a
/// set without one). Every clustering of either type runs on it, so both
/// share one implementation.
struct Clustered<'a> {
    views: Vec<ModelView<'a>>,
    /// The sample index of each member, which errors name.
    samples: &'a Members,
    /// Each member's relative weight.
    weights: Vec<f64>,
    /// Each member's state energy in the canonical unit, or why there is
    /// none.
    energies: Result<Vec<f64>, ClusteringError>,
}

/// The members of a set without a distribution, or of an ensemble: every
/// sample.
static EVERY_SAMPLE: Members = Members::All;

impl<'a> Clustered<'a> {
    fn of_ensemble(ensemble: &'a Ensemble) -> Self {
        let energies = ensemble
            .state_energies()
            .map(|energies| member_energies(&energies, ensemble.partition()))
            .map_err(ClusteringError::from);
        Self {
            views: ensemble
                .iter()
                .map(|sample| sample.as_model_view())
                .collect(),
            samples: &EVERY_SAMPLE,
            weights: ensemble.iter().map(|sample| sample.weight()).collect(),
            energies,
        }
    }

    fn of_set(set: &'a SampleSet) -> Self {
        let Some(distribution) = set.distribution() else {
            return Self {
                views: set.iter().map(|sample| sample.as_model_view()).collect(),
                samples: &EVERY_SAMPLE,
                weights: vec![1.0; set.len()],
                energies: Err(SampleSetError::Missing(Requirement::Distribution).into()),
            };
        };
        let (views, weights) = distribution
            .weighted_members()
            .map(|(sample, weight)| {
                let view = set.get(sample).expect("members are samples of the set");
                (view.as_model_view(), weight)
            })
            .unzip();
        let energies = distribution
            .state_energies()
            .map(|energies| member_energies(&energies, distribution.partition()))
            .map_err(ClusteringError::from);
        Self {
            views,
            samples: distribution.members(),
            weights,
            energies,
        }
    }

    /// The number of members.
    fn len(&self) -> usize {
        self.views.len()
    }

    /// Weights divided by the largest one, so sums cannot overflow.
    fn scaled_weights(&self) -> Vec<f64> {
        let maximum = self.weights.iter().copied().fold(0.0_f64, f64::max);
        self.weights.iter().map(|weight| weight / maximum).collect()
    }

    /// Member positions in visiting order.
    fn visiting_order(&self, order: ClusterOrder) -> Result<Vec<usize>, ClusteringError> {
        let mut positions = (0..self.len()).collect::<Vec<_>>();
        match order {
            ClusterOrder::Weight => {
                positions.sort_by(|&a, &b| self.weights[b].total_cmp(&self.weights[a]))
            }
            ClusterOrder::Energy => {
                let energies = self.energies.as_ref().map_err(Clone::clone)?;
                positions.sort_by(|&a, &b| energies[a].total_cmp(&energies[b]));
            }
        }
        Ok(positions)
    }

    fn cluster_by_radius(
        &self,
        radius: Quantity<f64>,
        atoms: FitAtoms<'_>,
        options: ClusteringOptions<'_>,
        run: ClusteringRun<'_>,
    ) -> Result<EnsemblePartition, ClusteringError> {
        let atoms = undirected(atoms)?;
        let radius = radius.value_in(CANONICAL_LENGTH_UNIT)?;
        if !(radius.is_finite() && radius >= 0.0) {
            return Err(ClusteringError::InvalidRadius);
        }
        let order = self.visiting_order(options.order)?;
        let (views, samples) = (&self.views, self.samples);
        let Some(first) = views.first() else {
            Tracker::start(run, 0, 0)?;
            return Ok(EnsemblePartition::state_per_sample(0));
        };
        let sets = atoms.pair_sets(*first, *first)?;
        check_fit_count(atoms.pair_count())?;
        let weights = NormalizedWeights::new(options.alignment.weighting(), atoms.pair_count())?;
        let policy = options.alignment.periodic_policy();
        let tracker = Tracker::start(run, views.len(), atoms.pair_count())?;

        // The nearest of `representatives`, numbered from `first_cluster`,
        // within `radius` of `member`; the earlier cluster on ties.
        let search = |member: usize, representatives: &[usize], first_cluster: usize| {
            let mut nearest: Option<(usize, f64)> = None;
            for (offset, representative) in representatives.iter().enumerate() {
                if offset > 0 && offset % tracker.fits_per_checkpoint() == 0 {
                    tracker.poll()?;
                }
                let rmsd = best_kabsch(
                    views[member],
                    views[*representative],
                    &sets,
                    weights,
                    policy,
                )
                .map_err(|source| AlignmentError::Item {
                    index: sample_index(samples, member),
                    source: Box::new(source),
                })?
                .rmsd()
                .into_value();
                if rmsd <= radius && nearest.is_none_or(|(_, best)| rmsd < best) {
                    nearest = Some((first_cluster + offset, rmsd));
                }
            }
            Ok::<_, ClusteringError>(nearest)
        };

        let parallelism = tracker.executor().parallelism().get();
        let mut assignments = vec![0; views.len()];
        let mut representatives: Vec<usize> = Vec::new();
        // Members placed and fits made since the last checkpoint.
        let (mut placed, mut fits) = (0, 0);
        let mut pending = order.as_slice();
        while !pending.is_empty() {
            // A batch of members is fitted at once against the `known`
            // representatives founded before it: a work item per thread, and
            // no more members than representatives, so that fitting members
            // to those founded within the batch stays the smaller share.
            let known = representatives.len();
            let members_per_item = tracker.fits_per_work_item().div_ceil(known.max(1));
            let size = (parallelism * members_per_item)
                .min(known.max(parallelism))
                .min(pending.len());
            let (batch, rest) = pending.split_at(size);
            pending = rest;
            let failure = FirstFailure::new();
            let ranges = map_ranges(
                tracker.executor(),
                batch.len(),
                members_per_item,
                |positions| {
                    let mut found = Vec::with_capacity(positions.len());
                    for position in positions {
                        if failure.precedes(position) {
                            return (found, Some(Stop::Superseded));
                        }
                        match search(batch[position], &representatives[..known], 0) {
                            Ok(nearest) => found.push(nearest),
                            Err(error) => {
                                failure.record(position);
                                return (found, Some(Stop::Error(error)));
                            }
                        }
                    }
                    (found, None)
                },
            );
            // Placing the batch in visiting order, each member also fitted to
            // the representatives founded before it within the batch, whose
            // clusters come after the known ones, makes every choice and
            // error the serial one.
            let mut members = batch.iter();
            for (found, stop) in ranges {
                for nearest in found {
                    let member = *members.next().expect("one result per member");
                    let founded = search(member, &representatives[known..], known)?;
                    fits += known.max(1) + (representatives.len() - known);
                    let nearest = match (nearest, founded) {
                        (Some(old), Some(new)) if new.1 < old.1 => Some(new),
                        (Some(old), _) => Some(old),
                        (None, new) => new,
                    };
                    assignments[member] = match nearest {
                        Some((cluster, _)) => cluster,
                        None => {
                            representatives.push(member);
                            representatives.len() - 1
                        }
                    };
                    placed += 1;
                    if fits >= tracker.fits_per_checkpoint() {
                        tracker.checkpoint(placed)?;
                        (placed, fits) = (0, 0);
                    }
                }
                if let Some(stop) = stop {
                    return Err(stop.into_error());
                }
            }
        }
        tracker.finish(placed);
        Ok(EnsemblePartition::from_parts_unchecked(
            assignments,
            representatives,
        ))
    }
}

/// Each member's state energy, from every state's energy and the members'
/// grouping into states.
fn member_energies(energies: &[f64], partition: &EnsemblePartition) -> Vec<f64> {
    partition
        .assignments()
        .iter()
        .map(|state| energies[*state])
        .collect()
}

/// The sample index of the member at `position` of `members`.
fn sample_index(members: &Members, position: usize) -> usize {
    match members {
        Members::All => position,
        Members::Subset(indices) => indices[position],
    }
}

impl SampleSet {
    /// Leader clustering on fitted RMSD within `radius`, a finite,
    /// nonnegative length, of the members of the set's distribution, or of
    /// every sample of a set without one.
    ///
    /// The result groups exactly those samples, by member position: entry
    /// `i` of its assignments, and each representative, is the `i`th member
    /// ([`Distribution::members`](crate::structure::Distribution::members)),
    /// or sample `i` without a distribution. Install it with
    /// [`DistributionMut::regroup`](crate::structure::DistributionMut::regroup),
    /// which keeps every weight, or, for a set without a distribution, as
    /// `Distribution::new(Members::All, partition)` through
    /// [`Self::set_distribution`]. Samples outside the distribution are not
    /// clustered.
    ///
    /// Visiting order, joining, and numbering are those of
    /// [`Ensemble::cluster_by_radius`]: by default members are visited by
    /// descending weight (in member order without a distribution, every
    /// member weighing the same), and [`ClusterOrder::Energy`] visits them by
    /// ascending energy of their states. Fits run on the calling thread;
    /// [`Self::cluster_by_radius_with`] spreads them over threads, reports
    /// progress, and can be cancelled.
    ///
    /// # Errors
    ///
    /// As [`Ensemble::cluster_by_radius`], with a failed fit's
    /// [`AlignmentError::Item`] naming the sample index of the member being
    /// placed, and [`ClusteringError::SampleSet`] wrapping
    /// `Missing(Distribution)` or `Missing(StateEnergies)` when ordering by
    /// [`ClusterOrder::Energy`] a set without a distribution or with a state
    /// without an energy.
    pub fn cluster_by_radius<'a>(
        &self,
        radius: Quantity<f64>,
        atoms: impl Into<FitAtoms<'a>>,
        options: ClusteringOptions<'_>,
    ) -> Result<EnsemblePartition, ClusteringError> {
        self.cluster_by_radius_with(radius, atoms, options, ClusteringRun::default())
    }

    /// [`Self::cluster_by_radius`], with its fits run, followed and stopped
    /// as `run` says, as [`Ensemble::cluster_by_radius_with`] describes.
    ///
    /// # Errors
    ///
    /// As [`Self::cluster_by_radius`], and [`ClusteringError::Cancelled`]
    /// once `run`'s cancellation is requested.
    ///
    /// # Panics
    ///
    /// When a work item panics, which
    /// [`Executor::run`](crate::parallel::Executor::run) passes on, or when
    /// the executor skips or repeats a work item.
    pub fn cluster_by_radius_with<'a>(
        &self,
        radius: Quantity<f64>,
        atoms: impl Into<FitAtoms<'a>>,
        options: ClusteringOptions<'_>,
        run: ClusteringRun<'_>,
    ) -> Result<EnsemblePartition, ClusteringError> {
        Clustered::of_set(self).cluster_by_radius(radius, atoms.into(), options, run)
    }
}

impl Ensemble {
    /// Leader clustering on fitted RMSD within `radius`, a finite,
    /// nonnegative length.
    ///
    /// Members are visited in the [`ClusterOrder`] of `options`, by default
    /// descending weight; ties keep member order. Each member joins the
    /// cluster whose representative is nearest by RMSD after a
    /// [`fit`](crate::alignment::fit) over `atoms` with the fit weights of
    /// `options`, by default uniform, if
    /// that RMSD is at most `radius`; ties keep the earlier cluster.
    /// Otherwise it founds a new cluster and represents it. Every member is
    /// therefore within `radius` of its representative, and representatives
    /// are farther than `radius` from one another. Clusters are numbered in
    /// order of founding.
    ///
    /// `atoms` is a selection or a
    /// [`SymmetricSelection`](crate::alignment::SymmetricSelection) of this
    /// ensemble's layout, so each distance is the same whichever member
    /// moves; an [`AtomCorrespondence`](crate::alignment::AtomCorrespondence)
    /// is rejected. Each member is fitted to at most every representative, so
    /// the cost grows with members times clusters (times symmetry
    /// permutations). An empty ensemble gives an empty partition. Fits run
    /// on the calling thread; [`Self::cluster_by_radius_with`] spreads them
    /// over threads, reports progress, and can be cancelled.
    ///
    /// # Errors
    ///
    /// - [`ClusteringError::Unit`] when `radius` is not a length, and
    ///   [`ClusteringError::InvalidRadius`] when it is negative or not finite.
    /// - [`ClusteringError::DirectedCorrespondence`] when `atoms` is an
    ///   [`AtomCorrespondence`](crate::alignment::AtomCorrespondence).
    /// - [`ClusteringError::Alignment`] when `atoms` does not belong to this
    ///   ensemble's layout, has fewer than three pairs, or the explicit
    ///   weights of `options` are invalid.
    ///   A failed fit is wrapped in [`AlignmentError::Item`] naming the
    ///   member being placed.
    /// - [`ClusteringError::Ensemble`] wrapping
    ///   [`EnsembleError::MissingEnergy`] when ordering by
    ///   [`ClusterOrder::Energy`] and a state has no energy.
    pub fn cluster_by_radius<'a>(
        &self,
        radius: Quantity<f64>,
        atoms: impl Into<FitAtoms<'a>>,
        options: ClusteringOptions<'_>,
    ) -> Result<EnsemblePartition, ClusteringError> {
        self.cluster_by_radius_with(radius, atoms, options, ClusteringRun::default())
    }

    /// [`Self::cluster_by_radius`], with its fits run, followed and stopped
    /// as `run` says.
    ///
    /// Each member joins or founds a cluster depending on those placed
    /// before it, so members are still placed one at a time in visiting
    /// order. Their fits run ahead in batches: every member of a batch is
    /// fitted, in ranges of members on `run`'s executor, to the
    /// representatives founded before the batch, and then, as it is placed,
    /// to the few founded within it. The nearest is chosen in cluster order
    /// and the first failed fit in visiting order is reported, so the
    /// partition or error is the one [`Self::cluster_by_radius`] returns.
    /// Progress counts members placed, of the member count.
    ///
    /// # Errors
    ///
    /// As [`Self::cluster_by_radius`], and [`ClusteringError::Cancelled`]
    /// once `run`'s cancellation is requested.
    ///
    /// # Panics
    ///
    /// When a work item panics, which
    /// [`Executor::run`](crate::parallel::Executor::run) passes on, or when
    /// the executor skips or repeats a work item.
    pub fn cluster_by_radius_with<'a>(
        &self,
        radius: Quantity<f64>,
        atoms: impl Into<FitAtoms<'a>>,
        options: ClusteringOptions<'_>,
        run: ClusteringRun<'_>,
    ) -> Result<EnsemblePartition, ClusteringError> {
        Clustered::of_ensemble(self).cluster_by_radius(radius, atoms.into(), options, run)
    }
}

/// Clustering distances must not depend on which member moves. Fitting
/// member `a` onto `b` pairs `a`'s moving atoms with `b`'s reference atoms,
/// so a correspondence pairing different atoms of one layout gives a
/// different RMSD in each direction, and results would follow member order.
/// Selections pair atoms with themselves and symmetry permutations form a
/// group, so both give the same RMSD either way.
fn undirected(atoms: FitAtoms<'_>) -> Result<FitAtoms<'_>, ClusteringError> {
    match atoms {
        FitAtoms::Correspondence(_) => Err(ClusteringError::DirectedCorrespondence),
        atoms => Ok(atoms),
    }
}

/// Failure of a sample set or ensemble clustering operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ClusteringError {
    /// A clustering radius is not finite and nonnegative.
    InvalidRadius,
    /// A distance matrix covers a different member count than there are
    /// members to cluster.
    MatrixMismatch {
        /// The number of members to cluster.
        expected: usize,
        /// The matrix's member count.
        actual: usize,
    },
    /// Condensed distances have the wrong length for their member count.
    CondensedLength {
        /// The stated member count.
        members: usize,
        /// The condensed length for that member count.
        expected: usize,
        /// The number of distances supplied.
        actual: usize,
    },
    /// A distance is not finite and nonnegative; `index` is condensed.
    InvalidDistance {
        /// The condensed index of the invalid distance.
        index: usize,
    },
    /// The pairwise distances of this many members are not addressable.
    MatrixTooLarge {
        /// The requested member count.
        members: usize,
    },
    /// Fitting one pair of members failed.
    PairFit {
        /// The sample index of the first member of the pair, the one fitted.
        first: usize,
        /// The sample index of the second member of the pair.
        second: usize,
        /// The fit failure.
        source: Box<AlignmentError>,
    },
    /// A dendrogram cut height is not finite and nonnegative.
    InvalidHeight,
    /// A requested cluster count is zero or exceeds the member count.
    InvalidClusterCount {
        /// The requested cluster count.
        requested: usize,
        /// The number of members.
        members: usize,
    },
    /// Distances between members were requested over an
    /// [`AtomCorrespondence`](crate::alignment::AtomCorrespondence), whose
    /// fits depend on which member moves.
    DirectedCorrespondence,
    /// Fitting members or validating the fit atoms failed; see the source.
    Alignment(Box<AlignmentError>),
    /// Reading or building the ensemble failed; see the source.
    Ensemble(Box<EnsembleError>),
    /// The sample set lacks what the clustering needs, such as an energy for
    /// every state; see the source.
    SampleSet(Box<SampleSetError>),
    /// A quantity could not be converted to the required unit; see the source.
    Unit(Box<UnitError>),
    /// The run's [`Cancellation`](crate::diagnostics::Cancellation) stopped
    /// it; no partial result is returned.
    Cancelled(Cancelled),
}

impl fmt::Display for ClusteringError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRadius => {
                formatter.write_str("The clustering radius is negative or not finite.")
            }
            Self::MatrixMismatch { expected, actual } => write!(
                formatter,
                "The distance matrix covers {}, but there are {expected} to cluster.",
                crate::diagnostics::counted(*actual, "member", "members")
            ),
            Self::CondensedLength {
                members,
                expected,
                actual,
            } => write!(
                formatter,
                "{members} members need {expected} condensed distances. {actual} were given."
            ),
            Self::InvalidDistance { index } => write!(
                formatter,
                "Condensed distance {index} is negative or not finite."
            ),
            Self::MatrixTooLarge { members } => write!(
                formatter,
                "The pairwise distances of {members} members are too many to store."
            ),
            Self::PairFit { first, second, .. } => {
                write!(formatter, "Cannot fit samples {first} and {second}.")
            }
            Self::InvalidHeight => {
                formatter.write_str("The dendrogram cut height is negative or not finite.")
            }
            Self::InvalidClusterCount { requested, members } => write!(
                formatter,
                "Cannot cut {members} members into {requested} clusters."
            ),
            Self::DirectedCorrespondence => formatter.write_str(
                "Distances between members need an atom selection, not an atom correspondence.",
            ),
            Self::Alignment(_) => formatter.write_str("Cannot fit the members."),
            Self::Ensemble(_) => formatter.write_str("Cannot read the ensemble."),
            Self::SampleSet(_) => formatter.write_str("Cannot cluster the sample set."),
            Self::Unit(_) => formatter.write_str("A clustering quantity has the wrong unit."),
            Self::Cancelled(_) => formatter.write_str("Clustering stopped."),
        }
    }
}

impl std::error::Error for ClusteringError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Alignment(error) | Self::PairFit { source: error, .. } => Some(error.as_ref()),
            Self::Ensemble(error) => Some(error.as_ref()),
            Self::SampleSet(error) => Some(error.as_ref()),
            Self::Unit(error) => Some(error.as_ref()),
            Self::Cancelled(cancelled) => Some(cancelled),
            _ => None,
        }
    }
}

impl From<AlignmentError> for ClusteringError {
    fn from(error: AlignmentError) -> Self {
        Self::Alignment(Box::new(error))
    }
}

impl From<EnsembleError> for ClusteringError {
    fn from(error: EnsembleError) -> Self {
        Self::Ensemble(Box::new(error))
    }
}

impl From<SampleSetError> for ClusteringError {
    fn from(error: SampleSetError) -> Self {
        Self::SampleSet(Box::new(error))
    }
}

impl From<StructureError> for ClusteringError {
    fn from(error: StructureError) -> Self {
        Self::Ensemble(Box::new(EnsembleError::Structure(error)))
    }
}

impl From<UnitError> for ClusteringError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
