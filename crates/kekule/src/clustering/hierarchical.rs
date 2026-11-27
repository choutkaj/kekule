//! Pairwise distance matrices and agglomerative hierarchical clustering.

use crate::alignment::{
    best_kabsch, AlignmentOptions, Coordinates, FitAtoms, NormalizedWeights, Pairs, PeriodicPolicy,
};
use crate::parallel::{map_ranges, FirstFailure, Run, Stop};
use crate::statistics::Weights;
use crate::structure::{GeometryStorage, StructureSet, StructureSetRead};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};

use super::driver::{undirected, Clustered, FittedAtoms};
use super::run::Tracker;
use super::{Clustering, ClusteringError};

/// Symmetric pairwise distances between the members of one clustering,
/// stored once per unordered pair in [`CANONICAL_LENGTH_UNIT`].
///
/// Members are numbered by member position, from zero. The condensed order
/// lists the pairs `(i, j)` with `i < j` row by row:
/// `(0, 1), (0, 2), ..., (0, n - 1), (1, 2), ...`, the order SciPy's
/// `squareform` uses. Build one from fitted RMSDs with
/// [`StructureSet::pairwise_fitted_rmsd`] or [`pairwise_fitted_rmsd`], or from
/// any other distances with
/// [`Self::from_condensed`]. It stores `n (n - 1) / 2` values, so memory
/// grows quadratically with the member count.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceMatrix {
    members: usize,
    values: Vec<f64>,
}

impl DistanceMatrix {
    /// Validates condensed distances between `members` members: exactly
    /// `members (members - 1) / 2` finite, nonnegative lengths.
    ///
    /// # Errors
    ///
    /// - [`ClusteringError::MatrixTooLarge`] when the pair count overflows.
    /// - [`ClusteringError::Unit`] when `distances` are not lengths.
    /// - [`ClusteringError::CondensedLength`] for the wrong number of
    ///   distances, and [`ClusteringError::InvalidDistance`] naming the first
    ///   negative or non-finite one.
    pub fn from_condensed(
        members: usize,
        distances: Quantity<Vec<f64>>,
    ) -> Result<Self, ClusteringError> {
        let expected = condensed_len(members)?;
        let values = distances.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
        if values.len() != expected {
            return Err(ClusteringError::CondensedLength {
                members,
                expected,
                actual: values.len(),
            });
        }
        if let Some(index) = values
            .iter()
            .position(|value| !(value.is_finite() && *value >= 0.0))
        {
            return Err(ClusteringError::InvalidDistance { index });
        }
        Ok(Self { members, values })
    }

    /// The number of members.
    pub fn member_count(&self) -> usize {
        self.members
    }

    /// Whether the matrix covers no members.
    pub fn is_empty(&self) -> bool {
        self.members == 0
    }

    /// The distance between two members in either order; zero for a member
    /// and itself, and `None` outside the matrix.
    pub fn get(&self, first: usize, second: usize) -> Option<Quantity<f64>> {
        if first >= self.members || second >= self.members {
            return None;
        }
        Some(Quantity::new(self.at(first, second), CANONICAL_LENGTH_UNIT))
    }

    /// Every distance in condensed order.
    pub fn condensed(&self) -> Quantity<&[f64]> {
        Quantity::new(self.values.as_slice(), CANONICAL_LENGTH_UNIT)
    }

    /// Agglomerative clustering of the members under `linkage`.
    ///
    /// Repeatedly merges the two closest clusters until one remains, using
    /// the nearest-neighbor chain algorithm in `O(n^2)` time on a copy of
    /// the distances. Every supported linkage is monotone, so merges come
    /// out in nondecreasing height. Equal distances are resolved
    /// deterministically: a nearest-neighbor search keeps the previous chain
    /// member, then the lowest index, and equal heights keep merge order.
    pub fn dendrogram(&self, linkage: Linkage) -> Dendrogram {
        let joins = nearest_neighbor_chain(self, linkage);
        Dendrogram::from_joins(self.members, joins)
    }

    /// The distance between two members in canonical units.
    pub(super) fn at(&self, first: usize, second: usize) -> f64 {
        match first.cmp(&second) {
            std::cmp::Ordering::Equal => 0.0,
            std::cmp::Ordering::Less => self.values[condensed_index(self.members, first, second)],
            std::cmp::Ordering::Greater => {
                self.values[condensed_index(self.members, second, first)]
            }
        }
    }
}

/// Distances rescaled for summing weighted objectives over the members.
///
/// Every distance is multiplied by one power of two, the largest that keeps
/// `4 n` times the largest distance finite: one unless distances come within
/// that factor of `f64::MAX`. Objectives sum at most `3 n` distances
/// weighted by at most one, so they stay finite with room for rounding.
/// Scaling by a power of two is exact barring underflow, and rounding
/// commutes with it, so scaled objectives compare and tie exactly as
/// unscaled ones would without overflow. Scaled values are not lengths in
/// any unit and only ever compare with one another.
#[derive(Debug, Clone, Copy)]
pub(super) struct ObjectiveDistances<'a> {
    matrix: &'a DistanceMatrix,
    scale: f64,
}

impl<'a> ObjectiveDistances<'a> {
    pub(super) fn new(matrix: &'a DistanceMatrix) -> Self {
        let largest = matrix.values.iter().copied().fold(0.0, f64::max);
        let terms = 4.0 * matrix.members as f64;
        let mut scale = 1.0;
        while largest * scale * terms > f64::MAX {
            scale *= 0.5;
        }
        Self { matrix, scale }
    }

    /// The number of members.
    pub(super) fn len(&self) -> usize {
        self.matrix.members
    }

    /// The scaled distance between two members in either order.
    pub(super) fn at(&self, first: usize, second: usize) -> f64 {
        self.matrix.at(first, second) * self.scale
    }
}

/// How the distance between two clusters follows from member distances.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Linkage {
    /// The closest pair of members; clusters chain through close neighbors.
    Single,
    /// The farthest pair of members; a cut at height `h` leaves clusters
    /// whose members are all within `h` of one another.
    Complete,
    /// The mean over all member pairs (UPGMA), counting each member once
    /// regardless of its weight.
    Average,
}

/// Where to cut a [`Dendrogram`] into flat clusters.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Cut {
    /// Apply every merge at most this finite, nonnegative height.
    Height(Quantity<f64>),
    /// Apply merges in order until this many clusters remain, from one to
    /// the member count (zero only for no members).
    Count(usize),
}

/// The merge tree of an agglomerative clustering.
///
/// Members are clusters `0..member_count`; merge `k` creates cluster
/// `member_count + k`, the numbering SciPy's linkage matrices use. Merges
/// are in nondecreasing height.
#[derive(Debug, Clone, PartialEq)]
pub struct Dendrogram {
    members: usize,
    merges: Vec<Merge>,
}

/// One merge of two clusters in a [`Dendrogram`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Merge {
    clusters: [usize; 2],
    height: f64,
    size: usize,
    /// One member of each merged cluster.
    members: [usize; 2],
}

impl Merge {
    /// The two merged clusters, lower number first.
    pub const fn clusters(&self) -> [usize; 2] {
        self.clusters
    }

    /// The linkage distance between them in [`CANONICAL_LENGTH_UNIT`].
    pub fn height(&self) -> Quantity<f64> {
        Quantity::new(self.height, CANONICAL_LENGTH_UNIT)
    }

    /// The number of members in the merged cluster.
    pub const fn size(&self) -> usize {
        self.size
    }
}

impl Dendrogram {
    /// The number of clustered members, the leaves of the tree.
    pub fn member_count(&self) -> usize {
        self.members
    }

    /// Every merge, in nondecreasing height.
    pub fn merges(&self) -> &[Merge] {
        &self.merges
    }

    /// The flat cluster of every member after `cut`. Clusters are numbered
    /// in order of their lowest member.
    ///
    /// # Errors
    ///
    /// For [`Cut::Height`], [`ClusteringError::Unit`] when the height is not a
    /// length and [`ClusteringError::InvalidHeight`] when it is negative or not
    /// finite; for [`Cut::Count`], [`ClusteringError::InvalidClusterCount`]
    /// when the count is outside one to the member count (zero only for no
    /// members).
    pub fn cut(&self, cut: Cut) -> Result<Vec<usize>, ClusteringError> {
        let applied = match cut {
            Cut::Height(height) => {
                let height = height.value_in(CANONICAL_LENGTH_UNIT)?;
                if !(height.is_finite() && height >= 0.0) {
                    return Err(ClusteringError::InvalidHeight);
                }
                self.merges
                    .iter()
                    .take_while(|merge| merge.height <= height)
                    .count()
            }
            Cut::Count(count) => {
                let valid = if self.members == 0 {
                    count == 0
                } else {
                    (1..=self.members).contains(&count)
                };
                if !valid {
                    return Err(ClusteringError::InvalidClusterCount {
                        requested: count,
                        members: self.members,
                    });
                }
                self.members.saturating_sub(count)
            }
        };
        let mut sets = DisjointSets::new(self.members);
        for merge in &self.merges[..applied] {
            sets.union(merge.members[0], merge.members[1]);
        }
        let mut labels = vec![usize::MAX; self.members];
        let mut assignments = Vec::with_capacity(self.members);
        let mut clusters = 0;
        for member in 0..self.members {
            let root = sets.find(member);
            if labels[root] == usize::MAX {
                labels[root] = clusters;
                clusters += 1;
            }
            assignments.push(labels[root]);
        }
        Ok(assignments)
    }

    /// The clustering of a set after `cut`, each cluster represented by its
    /// heaviest member, the lowest of equal weights. The members are the
    /// samples of `weights` by increasing index, and the clustering covers
    /// the weights' geometries; without weights, the members are the
    /// geometries `0..member_count`. It needs no distances, so a dendrogram
    /// kept from [`DistanceMatrix::dendrogram`] cuts again cheaply.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::cut`], and
    /// [`ClusteringError::MatrixMismatch`] when the dendrogram has another
    /// number of members than `weights` has samples.
    pub fn clustering(
        &self,
        cut: Cut,
        weights: Option<&Weights>,
    ) -> Result<Clustering, ClusteringError> {
        let clustered = Clustered::new(weights.map_or(self.members, Weights::len), weights)?;
        clustered.check_matrix(self.members)?;
        let assignments = self.cut(cut)?;
        let clusters = assignments.iter().max().map_or(0, |last| last + 1);
        let mut members = vec![Vec::new(); clusters];
        for (member, cluster) in assignments.iter().enumerate() {
            members[*cluster].push(member);
        }
        let representatives = members
            .iter()
            .map(|members| clustered.heaviest(members))
            .collect::<Vec<_>>();
        Ok(clustered.clustering(&assignments, &representatives))
    }

    /// Numbers merges found in any order, given as joined members and
    /// heights, after a stable sort by height.
    fn from_joins(members: usize, mut joins: Vec<(usize, usize, f64)>) -> Self {
        joins.sort_by(|a, b| a.2.total_cmp(&b.2));
        let mut sets = DisjointSets::new(members);
        let mut cluster = (0..members).collect::<Vec<_>>();
        let mut size = vec![1; members];
        let merges = joins
            .into_iter()
            .enumerate()
            .map(|(index, (first, second, height))| {
                let (a, b) = (sets.find(first), sets.find(second));
                let mut clusters = [cluster[a], cluster[b]];
                clusters.sort_unstable();
                let root = sets.union(a, b);
                cluster[root] = members + index;
                size[root] = size[a] + size[b];
                Merge {
                    clusters,
                    height,
                    size: size[root],
                    members: [first, second],
                }
            })
            .collect();
        Self { members, merges }
    }
}

impl Clustered {
    fn pairwise_fitted_rmsd(
        &self,
        reader: &dyn StructureSetRead,
        atoms: FitAtoms<'_>,
        options: AlignmentOptions<'_>,
        run: Run<'_>,
    ) -> Result<DistanceMatrix, ClusteringError> {
        let atoms = undirected(atoms)?;
        let count = self.len();
        let expected = condensed_len(count)?;
        let Some(&first) = self.members.first() else {
            Tracker::start(run, 0, 0)?;
            return Ok(DistanceMatrix {
                members: 0,
                values: Vec::new(),
            });
        };
        let fitted = FittedAtoms::new(reader, first, atoms)?;
        let weights = NormalizedWeights::new(options.weighting(), atoms.pair_count())?;
        let tracker = Tracker::start(run, expected, atoms.pair_count())?;
        let loaded = fitted.load(reader, &self.members, &tracker)?;
        let fits = PairFits {
            sets: &fitted.sets,
            weights,
            policy: options.periodic_policy(),
        };
        let coordinates = FittedAtoms::coordinates(&loaded);
        let values = fits.all_pairs(&coordinates, &self.members, &tracker)?;
        Ok(DistanceMatrix {
            members: count,
            values,
        })
    }

    fn cluster_hierarchical(
        &self,
        distances: &DistanceMatrix,
        linkage: Linkage,
        cut: Cut,
    ) -> Result<Clustering, ClusteringError> {
        self.check_matrix(distances.member_count())?;
        let assignments = distances.dendrogram(linkage).cut(cut)?;
        let representatives = weighted_medoids(&self.scaled_weights(), distances, &assignments);
        Ok(self.clustering(&assignments, &representatives))
    }
}

/// The fits between members: the pair sets, which index the members'
/// coordinates, their weights, and the periodic policy.
struct PairFits<'s, 'w> {
    sets: &'s [Pairs],
    weights: NormalizedWeights<'w>,
    policy: PeriodicPolicy,
}

impl PairFits<'_, '_> {
    /// The fitted RMSD of every unordered pair of `members`, in condensed
    /// order, with failures naming the geometries `geometries` lists; the
    /// one loop behind every pairwise matrix.
    fn all_pairs(
        &self,
        members: &[Coordinates<'_>],
        geometries: &[usize],
        tracker: &Tracker<'_>,
    ) -> Result<Vec<f64>, ClusteringError> {
        let count = members.len();
        let expected = condensed_len(count)?;
        let failure = FirstFailure::new();
        let ranges = map_ranges(
            tracker.executor(),
            expected,
            tracker.fits_per_work_item(),
            |pairs| {
                let (mut moving, mut reference) = condensed_pair(count, pairs.start);
                let mut values = Vec::with_capacity(pairs.len());
                let mut fits = 0;
                for pair in pairs {
                    if failure.precedes(pair) {
                        return Err(Stop::Superseded);
                    }
                    let fit = best_kabsch(
                        members[moving],
                        members[reference],
                        self.sets,
                        self.weights,
                        self.policy,
                    )
                    .map_err(|source| {
                        failure.record(pair);
                        ClusteringError::PairFit {
                            first: geometries[moving],
                            second: geometries[reference],
                            source: Box::new(source),
                        }
                    })?;
                    values.push(fit.rmsd().into_value());
                    fits += 1;
                    if fits == tracker.fits_per_checkpoint() {
                        tracker.checkpoint(fits)?;
                        fits = 0;
                    }
                    reference += 1;
                    if reference == count {
                        moving += 1;
                        reference = moving + 1;
                    }
                }
                tracker.finish(fits);
                Ok(values)
            },
        );
        let mut values = Vec::with_capacity(expected);
        for range in ranges {
            values.extend(range.map_err(Stop::into_error)?);
        }
        Ok(values)
    }
}

/// Fitted RMSD between every pair of members of `reader`, as
/// [`StructureSet::pairwise_fitted_rmsd`] computes it for a set, with the
/// same matrix and errors: the samples of `weights`, or every geometry
/// without weights. It fits any reader, such as a
/// [`Processed`](crate::traj::Processed) one.
///
/// ```
/// use kekule::clustering::pairwise_fitted_rmsd;
/// use kekule::math::Point3;
/// use kekule::structure::{Positions, Geometry, StructureSet};
/// use kekule::topology::AtomSelection;
/// use kekule::units::{Quantity, ANGSTROM};
///
/// let topology = kekule::smiles::to_topology("CCO")?;
/// let mut set = StructureSet::new(topology.clone());
/// for stretch in [1.0, 1.1, 1.2] {
///     let points = [
///         Point3::new(0.0, 0.0, 0.0),
///         Point3::new(1.5 * stretch, 0.0, 0.0),
///         Point3::new(2.0, 1.2, 0.0),
///     ];
///     set.push(Geometry::new(Positions::new(Quantity::new(points, ANGSTROM))?), None)?;
/// }
/// let atoms = AtomSelection::all(topology);
/// let matrix = pairwise_fitted_rmsd(&set, &atoms, None, Default::default())?;
/// assert_eq!(matrix, set.pairwise_fitted_rmsd(&atoms, None, Default::default())?);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
///
/// As [`StructureSet::pairwise_fitted_rmsd`].
pub fn pairwise_fitted_rmsd<'a>(
    reader: &dyn StructureSetRead,
    atoms: impl Into<FitAtoms<'a>>,
    weights: Option<&Weights>,
    options: AlignmentOptions<'_>,
) -> Result<DistanceMatrix, ClusteringError> {
    pairwise_fitted_rmsd_with(reader, atoms, weights, options, Run::default())
}

/// [`pairwise_fitted_rmsd`] under a [`Run`], as
/// [`StructureSet::pairwise_fitted_rmsd_with`] runs.
///
/// # Errors
///
/// As [`pairwise_fitted_rmsd`], and [`ClusteringError::Cancelled`] once
/// `run`'s cancellation is requested.
pub fn pairwise_fitted_rmsd_with<'a>(
    reader: &dyn StructureSetRead,
    atoms: impl Into<FitAtoms<'a>>,
    weights: Option<&Weights>,
    options: AlignmentOptions<'_>,
    run: Run<'_>,
) -> Result<DistanceMatrix, ClusteringError> {
    Clustered::new(reader.len(), weights)?.pairwise_fitted_rmsd(reader, atoms.into(), options, run)
}

/// The weighted medoid of every cluster of `assignments`, whose clusters
/// are numbered densely from zero, with members weighted by `weights`.
fn weighted_medoids(
    weights: &[f64],
    distances: &DistanceMatrix,
    assignments: &[usize],
) -> Vec<usize> {
    let distances = ObjectiveDistances::new(distances);
    let clusters = assignments.iter().max().map_or(0, |last| last + 1);
    let mut members = vec![Vec::new(); clusters];
    for (member, cluster) in assignments.iter().enumerate() {
        members[*cluster].push(member);
    }
    members
        .iter()
        .map(|cluster| {
            let mut best = (cluster[0], f64::INFINITY);
            for candidate in cluster {
                let cost = cluster
                    .iter()
                    .map(|other| weights[*other] * distances.at(*candidate, *other))
                    .sum::<f64>();
                if cost < best.1 {
                    best = (*candidate, cost);
                }
            }
            best.0
        })
        .collect()
}

impl<S: GeometryStorage> StructureSet<S>
where
    StructureSet<S>: StructureSetRead,
{
    /// Fitted RMSD between every pair of members, after a
    /// [`fit`](crate::alignment::fit) over `atoms` with the weights and
    /// periodic policy of `options`: the samples of `weights` by increasing
    /// index, or every geometry without weights.
    ///
    /// The matrix is indexed by member position: its member `i` is the
    /// `i`th sample, or geometry `i`. It fits [`Self::cluster_hierarchical`]
    /// and [`Self::cluster_k_medoids`] with the same weights.
    ///
    /// `atoms` is a selection or a
    /// [`SymmetricSelection`](crate::alignment::SymmetricSelection) of this
    /// set's layout, so each distance is the same whichever member moves; an
    /// [`AtomCorrespondence`](crate::alignment::AtomCorrespondence) is
    /// rejected. Each unordered pair is fitted once, lower member moving, so
    /// the cost is `n (n - 1) / 2` fits (times symmetry permutations). Each
    /// member is lent or read once and only its fitted atoms are kept, so
    /// memory grows with the members times the fitted atoms, and a set over
    /// a file fits too.
    ///
    /// ```
    /// use kekule::clustering::{Cut, Linkage};
    /// use kekule::math::Point3;
    /// use kekule::statistics::Weights;
    /// use kekule::structure::{Geometry, Positions, StructureSet};
    /// use kekule::topology::AtomSelection;
    /// use kekule::units::{Quantity, ANGSTROM};
    ///
    /// let topology = kekule::smiles::to_topology("CCO")?;
    /// let mut set = StructureSet::new(topology.clone());
    /// for stretch in [1.0, 1.0, 1.01, 1.3, 1.31] {
    ///     let points = [
    ///         Point3::new(0.0, 0.0, 0.0),
    ///         Point3::new(1.5 * stretch, 0.0, 0.0),
    ///         Point3::new(2.0, 1.2, 0.0),
    ///     ];
    ///     set.push(Geometry::new(Positions::new(Quantity::new(points, ANGSTROM))?), None)?;
    /// }
    /// // Geometries 1 to 4 are the samples; geometry 0 is left out.
    /// let weights = Weights::equal_over(5, [1, 2, 3, 4])?;
    ///
    /// let atoms = AtomSelection::all(topology);
    /// let matrix = set.pairwise_fitted_rmsd(&atoms, Some(&weights), Default::default())?;
    /// assert_eq!(matrix.member_count(), 4);
    /// let cut = Cut::Height(Quantity::new(0.1, ANGSTROM));
    /// let clustering = set.cluster_hierarchical(&matrix, Some(&weights), Linkage::Complete, cut)?;
    /// // By geometry: geometries 1 and 2, then 3 and 4.
    /// assert_eq!(clustering.cluster_of(2), Some(0));
    /// assert_eq!(clustering.cluster_of(3), Some(1));
    /// assert_eq!(clustering.cluster_of(0), None);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// - [`ClusteringError::Weights`] when `weights` cover another number of
    ///   geometries.
    /// - [`ClusteringError::DirectedCorrespondence`] when `atoms` is an
    ///   [`AtomCorrespondence`](crate::alignment::AtomCorrespondence).
    /// - [`ClusteringError::Alignment`] when `atoms` does not belong to this
    ///   set's layout, has fewer than three pairs, or the explicit weights of
    ///   `options` are invalid.
    /// - [`ClusteringError::PairFit`] naming both geometries when a fit
    ///   fails.
    /// - [`ClusteringError::MatrixTooLarge`] when the pair count overflows.
    /// - [`ClusteringError::Read`] when a member cannot be read; the first
    ///   failing member in order is reported.
    pub fn pairwise_fitted_rmsd<'a>(
        &self,
        atoms: impl Into<FitAtoms<'a>>,
        weights: Option<&Weights>,
        options: AlignmentOptions<'_>,
    ) -> Result<DistanceMatrix, ClusteringError> {
        self.pairwise_fitted_rmsd_with(atoms, weights, options, Run::default())
    }

    /// [`Self::pairwise_fitted_rmsd`] under a [`Run`].
    ///
    /// Members load in ranges on `run`'s executor, and the pairs split, in
    /// condensed order, into contiguous ranges whose distances join in
    /// order. Each distance is one fit whatever thread makes it, and a
    /// failed fit reports the first failing pair in condensed order, so the
    /// matrix and any error are the ones [`Self::pairwise_fitted_rmsd`]
    /// returns. Progress counts pairs fitted, of `n (n - 1) / 2`;
    /// cancellation is also polled before each member is loaded.
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use std::sync::atomic::{AtomicU64, Ordering};
    /// use kekule::math::Point3;
    /// use kekule::parallel::{Progress, Run, ScopedThreads};
    /// use kekule::structure::{Positions, Geometry, StructureSet};
    /// use kekule::topology::AtomSelection;
    /// use kekule::units::{Quantity, ANGSTROM};
    ///
    /// let topology = kekule::smiles::to_topology("CCO")?;
    /// let mut set = StructureSet::new(topology.clone());
    /// for stretch in [1.0, 1.1, 1.2, 1.3] {
    ///     let points = [
    ///         Point3::new(0.0, 0.0, 0.0),
    ///         Point3::new(1.5 * stretch, 0.0, 0.0),
    ///         Point3::new(2.0, 1.2, 0.0),
    ///     ];
    ///     set.push(Geometry::new(Positions::new(Quantity::new(points, ANGSTROM))?), None)?;
    /// }
    /// let atoms = AtomSelection::all(topology);
    ///
    /// let threads = ScopedThreads::new(NonZeroUsize::new(2).unwrap());
    /// let fitted = AtomicU64::new(0);
    /// let progress = |progress: Progress| {
    ///     fitted.store(progress.completed(), Ordering::Relaxed);
    /// };
    /// let run = Run::default()
    ///     .with_executor(&threads)
    ///     .with_progress(&progress);
    /// let matrix = set.pairwise_fitted_rmsd_with(&atoms, None, Default::default(), run)?;
    /// assert_eq!(matrix, set.pairwise_fitted_rmsd(&atoms, None, Default::default())?);
    /// assert_eq!(fitted.into_inner(), 6);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// As [`Self::pairwise_fitted_rmsd`], and [`ClusteringError::Cancelled`]
    /// once `run`'s cancellation is requested.
    pub fn pairwise_fitted_rmsd_with<'a>(
        &self,
        atoms: impl Into<FitAtoms<'a>>,
        weights: Option<&Weights>,
        options: AlignmentOptions<'_>,
        run: Run<'_>,
    ) -> Result<DistanceMatrix, ClusteringError> {
        pairwise_fitted_rmsd_with(self, atoms, weights, options, run)
    }
}

impl<S: GeometryStorage> StructureSet<S> {
    /// Hierarchical clustering of the members, on `distances` indexed by
    /// member position ([`Self::pairwise_fitted_rmsd`] with the same
    /// weights), cut into flat clusters, each represented by its weighted
    /// medoid: the samples of `weights`, or every geometry, each weighing
    /// the same, without weights.
    ///
    /// The dendrogram depends on distances only; see
    /// [`DistanceMatrix::dendrogram`]. The representative of a cluster is
    /// the member minimizing the weight-summed distance to the cluster's
    /// members, `sum_j w_j d(i, j)`, with exact weights, ties to the lowest
    /// index. Clusters are numbered in order of their lowest member.
    /// Geometries that are not samples are in no cluster.
    ///
    /// # Errors
    ///
    /// [`ClusteringError::Weights`] when `weights` cover another number of
    /// geometries, [`ClusteringError::MatrixMismatch`] when the matrix
    /// covers a different member count, and the errors of
    /// [`Dendrogram::cut`] for `cut`.
    pub fn cluster_hierarchical(
        &self,
        distances: &DistanceMatrix,
        weights: Option<&Weights>,
        linkage: Linkage,
        cut: Cut,
    ) -> Result<Clustering, ClusteringError> {
        Clustered::new(self.len(), weights)?.cluster_hierarchical(distances, linkage, cut)
    }
}

/// The number of unordered pairs of `members` members.
fn condensed_len(members: usize) -> Result<usize, ClusteringError> {
    members
        .checked_mul(members.saturating_sub(1))
        .map(|pairs| pairs / 2)
        .ok_or(ClusteringError::MatrixTooLarge { members })
}

/// The condensed position of the pair `(first, second)` with `first < second`.
fn condensed_index(members: usize, first: usize, second: usize) -> usize {
    first * members - first * (first + 1) / 2 + (second - first - 1)
}

/// The pair `(first, second)` at condensed position `index`, the inverse of
/// [`condensed_index`]; `(members - 1, members)` just past the last pair.
fn condensed_pair(members: usize, mut index: usize) -> (usize, usize) {
    let mut first = 0;
    while first + 1 < members && index >= members - 1 - first {
        index -= members - 1 - first;
        first += 1;
    }
    (first, first + 1 + index)
}

/// Joins of the nearest-neighbor chain algorithm (Müllner 2011, sec. 3.2),
/// as `(member, member, height)` in the order they are found.
///
/// Slot `s` always holds the cluster containing member `s`; a merge keeps
/// the lower slot and retires the higher one.
fn nearest_neighbor_chain(
    distances: &DistanceMatrix,
    linkage: Linkage,
) -> Vec<(usize, usize, f64)> {
    let members = distances.members;
    let mut working = distances.values.clone();
    let at = |working: &[f64], first: usize, second: usize| {
        let (low, high) = if first < second {
            (first, second)
        } else {
            (second, first)
        };
        working[condensed_index(members, low, high)]
    };
    let mut active = vec![true; members];
    let mut size = vec![1_usize; members];
    let mut joins = Vec::with_capacity(members.saturating_sub(1));
    let mut chain: Vec<usize> = Vec::new();
    for _ in 1..members {
        if chain.is_empty() {
            chain.push(
                active
                    .iter()
                    .position(|alive| *alive)
                    .expect("two clusters remain"),
            );
        }
        let (a, b) = loop {
            let current = *chain.last().expect("nonempty chain");
            let previous = chain.len().checked_sub(2).map(|index| chain[index]);
            // Keeping the previous chain member on ties guarantees progress.
            let mut nearest = previous.map(|slot| (slot, at(&working, current, slot)));
            for (slot, alive) in active.iter().enumerate() {
                if !*alive || slot == current {
                    continue;
                }
                let distance = at(&working, current, slot);
                if nearest.is_none_or(|(_, best)| distance < best) {
                    nearest = Some((slot, distance));
                }
            }
            let (next, _) = nearest.expect("another active cluster");
            if Some(next) == previous {
                chain.truncate(chain.len() - 2);
                break (current, next);
            }
            chain.push(next);
        };
        let (low, high) = if a < b { (a, b) } else { (b, a) };
        joins.push((low, high, at(&working, low, high)));
        for (slot, alive) in active.iter().enumerate() {
            if !*alive || slot == low || slot == high {
                continue;
            }
            let (to_low, to_high) = (at(&working, low, slot), at(&working, high, slot));
            let merged = match linkage {
                Linkage::Single => to_low.min(to_high),
                Linkage::Complete => to_low.max(to_high),
                Linkage::Average => {
                    // The size-weighted mean, stepped from one input towards
                    // the other so no intermediate exceeds the larger input;
                    // summing scaled inputs first overflows near `f64::MAX`.
                    // Clamping keeps rounding from leaving the inputs' range.
                    let (low_size, high_size) = (size[low] as f64, size[high] as f64);
                    let fraction = high_size / (low_size + high_size);
                    (to_low + (to_high - to_low) * fraction)
                        .clamp(to_low.min(to_high), to_low.max(to_high))
                }
            };
            let (first, second) = if low < slot { (low, slot) } else { (slot, low) };
            working[condensed_index(members, first, second)] = merged;
        }
        size[low] += size[high];
        active[high] = false;
    }
    joins
}

/// Union-find over members with path halving.
struct DisjointSets {
    parent: Vec<usize>,
}

impl DisjointSets {
    fn new(members: usize) -> Self {
        Self {
            parent: (0..members).collect(),
        }
    }

    fn find(&mut self, mut member: usize) -> usize {
        while self.parent[member] != member {
            self.parent[member] = self.parent[self.parent[member]];
            member = self.parent[member];
        }
        member
    }

    /// Joins two sets and returns the root of the union, the lower root.
    fn union(&mut self, first: usize, second: usize) -> usize {
        let (a, b) = (self.find(first), self.find(second));
        let (root, child) = if a < b { (a, b) } else { (b, a) };
        self.parent[child] = root;
        root
    }
}

#[cfg(test)]
mod tests;
