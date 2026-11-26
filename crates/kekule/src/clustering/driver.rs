//! Leader clustering of structure sets, the members a clustering groups, and
//! errors.

use std::fmt;

use super::run::Tracker;
use super::Clustering;
use crate::alignment::{
    best_kabsch, check_fit_count, AlignmentError, AlignmentOptions, Coordinates, FitAtoms,
    NormalizedWeights, Pairs,
};
use crate::diagnostics::{Cancelled, IndexOutOfRange};
use crate::fixed_u32_id::FixedId;
use crate::math::Point3;
use crate::parallel::{map_ranges, FirstFailure, Run, Stop};
use crate::statistics::{Weights, WeightsError};
use crate::structure::{
    visit_geometry, GeometryStorage, StructureSet, StructureSetRead, StructureSetReadError,
};
use crate::topology::TopologyAtomIndex;
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};

/// The geometries a clustering groups, in member order: the samples of the
/// weights by increasing index, or every geometry without weights.
pub(super) struct Clustered {
    /// The number of geometries of the set, which the clustering covers.
    len: usize,
    /// The geometry index of each member.
    pub(super) members: Vec<usize>,
    /// Each member's natural-logarithm weight.
    log_weights: Vec<f64>,
}

impl Clustered {
    /// The members of a set of `len` geometries under `weights`.
    pub(super) fn new(len: usize, weights: Option<&Weights>) -> Result<Self, ClusteringError> {
        let Some(weights) = weights else {
            return Ok(Self {
                len,
                members: (0..len).collect(),
                log_weights: vec![0.0; len],
            });
        };
        weights.check_len(len)?;
        let members = weights.samples().collect::<Vec<_>>();
        let log_weights = members
            .iter()
            .map(|geometry| {
                weights
                    .log_weight(*geometry)
                    .expect("a sample has a log weight")
            })
            .collect();
        Ok(Self {
            len,
            members,
            log_weights,
        })
    }

    /// The number of members.
    pub(super) fn len(&self) -> usize {
        self.members.len()
    }

    /// Weights divided by the largest one, so the largest is one and sums
    /// cannot overflow.
    pub(super) fn scaled_weights(&self) -> Vec<f64> {
        let maximum = self
            .log_weights
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        self.log_weights
            .iter()
            .map(|log_weight| (log_weight - maximum).exp())
            .collect()
    }

    /// The member heaviest of `positions`, the first of equal weights.
    pub(super) fn heaviest(&self, positions: &[usize]) -> usize {
        let mut best = positions[0];
        for &position in &positions[1..] {
            if self.log_weights[position] > self.log_weights[best] {
                best = position;
            }
        }
        best
    }

    /// Checks that `distances` covers the members.
    pub(super) fn check_matrix(&self, member_count: usize) -> Result<(), ClusteringError> {
        if member_count != self.len() {
            return Err(ClusteringError::MatrixMismatch {
                expected: self.len(),
                actual: member_count,
            });
        }
        Ok(())
    }

    /// The clustering of the set from each member's cluster and each
    /// cluster's representative member, both by member position; clusters
    /// are numbered densely from zero and each contains its representative.
    pub(super) fn clustering(
        &self,
        assignments: &[usize],
        representatives: &[usize],
    ) -> Clustering {
        let mut cluster_of = vec![None; self.len];
        for (position, cluster) in assignments.iter().enumerate() {
            cluster_of[self.members[position]] = Some(*cluster);
        }
        let representatives = representatives
            .iter()
            .map(|position| self.members[*position])
            .collect();
        Clustering::from_parts_unchecked(cluster_of, representatives)
    }

    fn cluster_by_radius(
        &self,
        reader: &dyn StructureSetRead,
        radius: Quantity<f64>,
        atoms: FitAtoms<'_>,
        options: AlignmentOptions<'_>,
        run: Run<'_>,
    ) -> Result<Clustering, ClusteringError> {
        let atoms = undirected(atoms)?;
        let radius = radius.value_in(CANONICAL_LENGTH_UNIT)?;
        if !(radius.is_finite() && radius >= 0.0) {
            return Err(ClusteringError::InvalidRadius);
        }
        // Descending weight; ties keep member order.
        let mut order = (0..self.len()).collect::<Vec<_>>();
        order.sort_by(|&a, &b| self.log_weights[b].total_cmp(&self.log_weights[a]));
        let Some(&first) = self.members.first() else {
            Tracker::start(run, 0, 0)?;
            return Ok(self.clustering(&[], &[]));
        };
        let fitted = FittedAtoms::new(reader, first, atoms)?;
        let weights = NormalizedWeights::new(options.weighting(), atoms.pair_count())?;
        let policy = options.periodic_policy();
        let tracker = Tracker::start(run, self.len(), atoms.pair_count())?;
        let loaded = fitted.load(reader, &self.members, &tracker)?;
        let coordinates = FittedAtoms::coordinates(&loaded);
        let sets = fitted.sets.as_slice();

        // The nearest of `representatives`, numbered from `first_cluster`,
        // within `radius` of `member`; the earlier cluster on ties.
        let search = |member: usize, representatives: &[usize], first_cluster: usize| {
            let mut nearest: Option<(usize, f64)> = None;
            for (offset, representative) in representatives.iter().enumerate() {
                if offset > 0 && offset % tracker.fits_per_checkpoint() == 0 {
                    tracker.poll()?;
                }
                let rmsd = best_kabsch(
                    coordinates[member],
                    coordinates[*representative],
                    sets,
                    weights,
                    policy,
                )
                .map_err(|source| AlignmentError::Item {
                    index: self.members[member],
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
        let mut assignments = vec![0; self.len()];
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
        Ok(self.clustering(&assignments, &representatives))
    }
}

/// Members loaded by one work item at least: reading a geometry costs
/// about as much as many fits.
const MEMBERS_PER_WORK_ITEM: usize = 16;

/// The atoms fits between members read, and the pair sets over them.
///
/// Only the fitted atoms of each member are kept: every member is lent or
/// read once and its fitted positions copied out, so memory grows with the
/// members times the fitted atoms rather than with whole geometries.
pub(super) struct FittedAtoms {
    /// The fitted atoms, sorted.
    fitted: Vec<usize>,
    /// The candidate pair sets, indexing [`Self::fitted`].
    pub(super) sets: Vec<Pairs>,
}

/// One member's fitted positions and whether it is periodic.
pub(super) type Loaded = (Vec<Point3>, bool);

impl FittedAtoms {
    /// Checks `atoms` against `geometry` of `reader` and compacts its pair
    /// sets onto the fitted atoms.
    pub(super) fn new(
        reader: &dyn StructureSetRead,
        geometry: usize,
        atoms: FitAtoms<'_>,
    ) -> Result<Self, ClusteringError> {
        let sets = visit_geometry(reader, geometry, ClusteringError::Read, |first| {
            Ok(atoms.pair_sets(first, first)?)
        })?;
        check_fit_count(atoms.pair_count())?;
        let mut fitted = sets
            .iter()
            .flatten()
            .flat_map(|(moving, reference)| [moving.index(), reference.index()])
            .collect::<Vec<_>>();
        fitted.sort_unstable();
        fitted.dedup();
        let compact = |atom: TopologyAtomIndex| {
            let position = fitted
                .binary_search(&atom.index())
                .expect("every paired atom is fitted");
            TopologyAtomIndex::from_index(position)
        };
        let sets = sets
            .iter()
            .map(|pairs| {
                pairs
                    .iter()
                    .map(|(moving, reference)| (compact(*moving), compact(*reference)))
                    .collect()
            })
            .collect::<Vec<Pairs>>();
        Ok(Self { fitted, sets })
    }

    /// The fitted positions of the geometries `members` of `reader`, in
    /// that order, loaded in ranges on the tracker's executor; cancellation
    /// is polled before each one, and the first failing read in order is
    /// reported.
    pub(super) fn load(
        &self,
        reader: &dyn StructureSetRead,
        members: &[usize],
        tracker: &Tracker<'_>,
    ) -> Result<Vec<Loaded>, ClusteringError> {
        let loaded = map_ranges(
            tracker.executor(),
            members.len(),
            MEMBERS_PER_WORK_ITEM,
            |positions| {
                positions
                    .map(|position| {
                        tracker.poll()?;
                        visit_geometry(reader, members[position], ClusteringError::Read, |view| {
                            let points = *view.positions().values().value();
                            let points = self
                                .fitted
                                .iter()
                                .map(|atom| points[*atom])
                                .collect::<Vec<Point3>>();
                            Ok((points, view.cell().is_some()))
                        })
                    })
                    .collect::<Result<Vec<_>, ClusteringError>>()
            },
        );
        let mut all = Vec::with_capacity(members.len());
        for range in loaded {
            all.extend(range?);
        }
        Ok(all)
    }

    /// The coordinates the kernels read from loaded members.
    pub(super) fn coordinates(loaded: &[Loaded]) -> Vec<Coordinates<'_>> {
        loaded
            .iter()
            .map(|(points, periodic)| Coordinates {
                points,
                periodic: *periodic,
            })
            .collect()
    }
}

impl<S: GeometryStorage> StructureSet<S>
where
    StructureSet<S>: StructureSetRead,
{
    /// Leader clustering on fitted RMSD within `radius`, a finite,
    /// nonnegative length.
    ///
    /// With `weights`, it clusters their samples, visited by descending
    /// weight, so the heaviest sample of a cluster leads it; geometries that
    /// are not samples are in no cluster. Without, it clusters every
    /// geometry in index order. Ties keep index order.
    ///
    /// Each member joins the cluster whose representative is nearest by
    /// RMSD after a [`fit`](crate::alignment::fit) over `atoms` with the fit
    /// weights of `options`, by default uniform, if that RMSD is at most
    /// `radius`; ties keep the earlier cluster. Otherwise it founds a new
    /// cluster and represents it. Every member is therefore within `radius`
    /// of its representative, and representatives are farther than `radius`
    /// from one another. Clusters are numbered in order of founding.
    ///
    /// `atoms` is a selection or a
    /// [`SymmetricSelection`](crate::alignment::SymmetricSelection) of this
    /// set's layout, so each distance is the same whichever member moves; an
    /// [`AtomCorrespondence`](crate::alignment::AtomCorrespondence) is
    /// rejected. Each member is fitted to at most every representative, so
    /// the cost grows with members times clusters (times symmetry
    /// permutations). Each member is lent or read once, and only its fitted
    /// atoms are kept, so a set over a file clusters too. No member gives
    /// no cluster.
    ///
    /// # Errors
    ///
    /// - [`ClusteringError::Weights`] when `weights` cover another number of
    ///   geometries.
    /// - [`ClusteringError::Unit`] when `radius` is not a length, and
    ///   [`ClusteringError::InvalidRadius`] when it is negative or not finite.
    /// - [`ClusteringError::DirectedCorrespondence`] when `atoms` is an
    ///   [`AtomCorrespondence`](crate::alignment::AtomCorrespondence).
    /// - [`ClusteringError::Alignment`] when `atoms` does not belong to this
    ///   set's layout, has fewer than three pairs, or the explicit weights of
    ///   `options` are invalid. A failed fit is wrapped in
    ///   [`AlignmentError::Item`] naming the geometry being placed.
    /// - [`ClusteringError::Read`] when a member cannot be read.
    pub fn cluster_by_radius<'a>(
        &self,
        radius: Quantity<f64>,
        atoms: impl Into<FitAtoms<'a>>,
        weights: Option<&Weights>,
        options: AlignmentOptions<'_>,
    ) -> Result<Clustering, ClusteringError> {
        self.cluster_by_radius_with(radius, atoms, weights, options, Run::default())
    }

    /// [`Self::cluster_by_radius`] under a [`Run`].
    ///
    /// Each member joins or founds a cluster depending on those placed
    /// before it, so members are still placed one at a time in visiting
    /// order. Their fits run ahead in batches: every member of a batch is
    /// fitted, in ranges of members on `run`'s executor, to the
    /// representatives founded before the batch, and then, as it is placed,
    /// to the few founded within it. The nearest is chosen in cluster order
    /// and the first failed fit in visiting order is reported, so the
    /// clustering or error is the one [`Self::cluster_by_radius`] returns.
    /// Progress counts members placed, of the member count.
    ///
    /// # Errors
    ///
    /// As [`Self::cluster_by_radius`], and [`ClusteringError::Cancelled`]
    /// once `run`'s cancellation is requested.
    pub fn cluster_by_radius_with<'a>(
        &self,
        radius: Quantity<f64>,
        atoms: impl Into<FitAtoms<'a>>,
        weights: Option<&Weights>,
        options: AlignmentOptions<'_>,
        run: Run<'_>,
    ) -> Result<Clustering, ClusteringError> {
        Clustered::new(self.len(), weights)?.cluster_by_radius(
            self,
            radius,
            atoms.into(),
            options,
            run,
        )
    }
}

/// Clustering distances must not depend on which member moves. Fitting
/// member `a` onto `b` pairs `a`'s moving atoms with `b`'s reference atoms,
/// so a correspondence pairing different atoms of one layout gives a
/// different RMSD in each direction, and results would follow member order.
/// Selections pair atoms with themselves and symmetry permutations form a
/// group, so both give the same RMSD either way.
pub(super) fn undirected(atoms: FitAtoms<'_>) -> Result<FitAtoms<'_>, ClusteringError> {
    match atoms {
        FitAtoms::Correspondence(_) => Err(ClusteringError::DirectedCorrespondence),
        atoms => Ok(atoms),
    }
}

/// Failure of a clustering operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ClusteringError {
    /// A clustering radius is not finite and nonnegative.
    InvalidRadius,
    /// A distance matrix or dendrogram covers a different member count than
    /// there are members to cluster.
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
        /// The geometry index of the first member of the pair, the one fitted.
        first: usize,
        /// The geometry index of the second member of the pair.
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
    /// A geometry or cluster index is outside the clustering.
    IndexOutOfRange(IndexOutOfRange),
    /// A geometry of a [`Clustering`] is in a cluster without a
    /// representative.
    ClusterOutOfRange {
        /// The geometry's index.
        geometry: usize,
        /// Its cluster.
        cluster: usize,
        /// The number of clusters, one per representative.
        clusters: usize,
    },
    /// A representative of a [`Clustering`] is not in its cluster.
    RepresentativeNotMember {
        /// The cluster.
        cluster: usize,
        /// The representative geometry.
        geometry: usize,
    },
    /// The weights do not fit, or cannot give what was asked; see the
    /// cause.
    Weights(WeightsError),
    /// Fitting members or validating the fit atoms failed; see the source.
    Alignment(Box<AlignmentError>),
    /// A quantity could not be converted to the required unit; see the source.
    Unit(Box<UnitError>),
    /// The run's [`Cancellation`](crate::diagnostics::Cancellation) stopped
    /// it; no partial result is returned.
    Cancelled(Cancelled),
    /// A member cannot be read; see the source.
    Read(StructureSetReadError),
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
                write!(formatter, "Cannot fit geometries {first} and {second}.")
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
            Self::IndexOutOfRange(error) => fmt::Display::fmt(error, formatter),
            Self::ClusterOutOfRange {
                geometry, cluster, ..
            } => write!(
                formatter,
                "Geometry {geometry} is in cluster {cluster}, which has no representative."
            ),
            Self::RepresentativeNotMember { cluster, geometry } => write!(
                formatter,
                "Geometry {geometry} represents cluster {cluster} but is not in it."
            ),
            Self::Weights(error) => fmt::Display::fmt(error, formatter),
            Self::Alignment(_) => formatter.write_str("Cannot fit the members."),
            Self::Unit(_) => formatter.write_str("A clustering quantity has the wrong unit."),
            Self::Cancelled(_) => formatter.write_str("Clustering stopped."),
            Self::Read(error) => fmt::Display::fmt(error, formatter),
        }
    }
}

impl std::error::Error for ClusteringError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Alignment(error) | Self::PairFit { source: error, .. } => Some(error.as_ref()),
            Self::Unit(error) => Some(error.as_ref()),
            Self::Cancelled(cancelled) => Some(cancelled),
            Self::Weights(error) => std::error::Error::source(error),
            Self::Read(error) => std::error::Error::source(error),
            _ => None,
        }
    }
}

impl From<AlignmentError> for ClusteringError {
    fn from(error: AlignmentError) -> Self {
        Self::Alignment(Box::new(error))
    }
}

impl From<UnitError> for ClusteringError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}

impl From<WeightsError> for ClusteringError {
    fn from(error: WeightsError) -> Self {
        Self::Weights(error)
    }
}
