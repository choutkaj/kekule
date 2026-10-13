//! Partitions of an [`Ensemble`] into clusters of similar members.
//!
//! Clustering never changes the ensemble. It returns an
//! [`EnsemblePartition`], a separate result that assigns every member to one
//! cluster and names one representative member per cluster.
//! [`Ensemble::collapse`] then keeps one member per cluster, carrying the
//! summed weight of its cluster.
//!
//! [`Ensemble::cluster_hierarchical`] cuts the [`Dendrogram`] of a
//! [`DistanceMatrix`], such as [`Ensemble::pairwise_rmsd`], at a height or
//! into a number of clusters, and represents each cluster by its weighted
//! medoid. Complete linkage cut at height `h` bounds every cluster's
//! diameter by `h`.
//!
//! [`Ensemble::cluster_by_radius`] is leader clustering on fitted RMSD. It
//! needs no matrix, so it scales to larger ensembles. With
//! a [`SymmetricSelection`](crate::alignment::SymmetricSelection) and a small
//! radius it merges rotamers, copies of one conformer that differ only by
//! exchanged equivalent atoms. Weight members by energy before collapsing, so
//! each conformer's weight includes its degeneracy:
//!
//! ```
//! use std::error::Error;
//!
//! use kekule::alignment::SymmetricSelection;
//! use kekule::geometry::Point3;
//! use kekule::structure::{Ensemble, EnsembleMember, Positions};
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, ANGSTROM, KELVIN, KILOJOULE_PER_MOLE};
//!
//! // Acetate, then a rotamer with its oxygens exchanged, then a distinct
//! // conformer 5 kJ/mol higher.
//! let topology = kekule::smiles::to_topology("CC(=O)[O-]")?;
//! let member = |o2: f64, o3: f64, energy: f64| -> Result<EnsembleMember, Box<dyn Error>> {
//!     let points = [
//!         Point3::new(0.0, 0.0, 0.0),
//!         Point3::new(1.5, 0.0, 0.0),
//!         Point3::new(2.2, o2, 0.5),
//!         Point3::new(2.2, o3, 0.5),
//!     ];
//!     let positions = Positions::new(Quantity::new(points, ANGSTROM))?;
//!     Ok(EnsembleMember::new(positions, 1.0)?
//!         .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))?)
//! };
//! let mut ensemble = Ensemble::from_items(
//!     topology.clone(),
//!     [member(1.1, -1.1, 0.0)?, member(-1.1, 1.1, 0.0)?, member(1.1, -0.3, 5.0)?],
//! )?;
//! ensemble.assign_boltzmann_weights(Quantity::new(298.15, KELVIN))?;
//!
//! let atoms = SymmetricSelection::new(&AtomSelection::all(&topology))?;
//! let rotamers = ensemble.cluster_by_radius(Quantity::new(0.125, ANGSTROM), &atoms)?;
//! assert_eq!(rotamers.assignments(), [0, 0, 1]);
//!
//! let conformers = ensemble.collapse(&rotamers)?;
//! assert_eq!(conformers.len(), 2);
//! assert_eq!(conformers.get(0).unwrap().weight(), 2.0);
//! # Ok::<(), Box<dyn Error>>(())
//! ```

use std::fmt;

mod hierarchical;
pub use hierarchical::{Cut, Dendrogram, DistanceMatrix, Linkage, Merge};

use crate::alignment::{
    best_kabsch, check_fit_count, AlignmentError, AlignmentOptions, FitAtoms, NormalizedWeights,
};
use crate::structure::{Ensemble, EnsembleError, ModelView, StructureError};
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};

/// An assignment of every member of one ensemble to exactly one cluster,
/// with one representative member per cluster.
///
/// Clusters are numbered from zero; a cluster always contains its
/// representative, so none is empty. The partition records member indices
/// only: it applies to any ensemble with the same member count, and
/// [`Ensemble::collapse`] checks that count.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsemblePartition {
    assignments: Vec<usize>,
    representatives: Vec<usize>,
}

impl EnsemblePartition {
    /// Validates a partition built elsewhere: `assignments[m]` is the cluster
    /// of member `m`, and `representatives[c]` is a member assigned to
    /// cluster `c`.
    pub fn new(
        assignments: Vec<usize>,
        representatives: Vec<usize>,
    ) -> Result<Self, PartitionError> {
        for (member, cluster) in assignments.iter().enumerate() {
            if *cluster >= representatives.len() {
                return Err(PartitionError::ClusterOutOfRange {
                    member,
                    cluster: *cluster,
                    clusters: representatives.len(),
                });
            }
        }
        for (cluster, member) in representatives.iter().enumerate() {
            match assignments.get(*member) {
                None => {
                    return Err(PartitionError::RepresentativeOutOfRange {
                        cluster,
                        member: *member,
                        members: assignments.len(),
                    })
                }
                Some(assigned) if *assigned != cluster => {
                    return Err(PartitionError::RepresentativeNotInCluster {
                        cluster,
                        member: *member,
                    })
                }
                Some(_) => {}
            }
        }
        Ok(Self {
            assignments,
            representatives,
        })
    }

    /// The number of partitioned members.
    pub fn member_count(&self) -> usize {
        self.assignments.len()
    }

    pub fn cluster_count(&self) -> usize {
        self.representatives.len()
    }

    /// The cluster of every member, in member order.
    pub fn assignments(&self) -> &[usize] {
        &self.assignments
    }

    /// The representative member of every cluster, in cluster order.
    pub fn representatives(&self) -> &[usize] {
        &self.representatives
    }

    /// The members of one cluster in member order; empty for a cluster that
    /// does not exist.
    pub fn members(&self, cluster: usize) -> impl Iterator<Item = usize> + '_ {
        self.assignments
            .iter()
            .enumerate()
            .filter(move |(_, assigned)| **assigned == cluster)
            .map(|(member, _)| member)
    }
}

/// The order in which leader clustering visits members.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum ClusterOrder {
    /// Descending weight, so the most populated member of a cluster leads it.
    #[default]
    Weight,
    /// Ascending energy, so the lowest-energy member of a cluster leads it.
    Energy,
}

/// Options for [`Ensemble::cluster_by_radius_with_options`].
#[derive(Debug, Clone, Copy, Default)]
pub struct ClusteringOptions<'a> {
    /// Fit weighting and periodic policy of every pairwise fit.
    pub alignment: AlignmentOptions<'a>,
    /// The order in which members are visited.
    pub order: ClusterOrder,
}

impl Ensemble {
    /// Leader clustering on fitted RMSD within `radius`.
    ///
    /// See [`Self::cluster_by_radius_with_options`]; members are visited by
    /// descending weight with uniform fit weights.
    pub fn cluster_by_radius<'a>(
        &self,
        radius: Quantity<f64>,
        atoms: impl Into<FitAtoms<'a>>,
    ) -> Result<EnsemblePartition, ClusteringError> {
        self.cluster_by_radius_with_options(radius, atoms, ClusteringOptions::default())
    }

    /// Leader clustering on fitted RMSD within `radius`, a finite,
    /// nonnegative length.
    ///
    /// Members are visited in [`ClusterOrder`]; ties keep member order. Each
    /// member joins the cluster whose representative is nearest by RMSD
    /// after a [`kabsch`](crate::alignment::kabsch) fit over `atoms`, if
    /// that RMSD is at most `radius`; ties keep the earlier cluster.
    /// Otherwise it founds a new cluster and represents it. Every member is
    /// therefore within `radius` of its representative, and representatives
    /// are farther than `radius` from one another. Clusters are numbered in
    /// order of founding.
    ///
    /// `atoms` pairs members of this ensemble's layout, typically a
    /// selection or a [`SymmetricSelection`](crate::alignment::SymmetricSelection).
    /// Each member is fitted to at most every representative, so the cost
    /// grows with members times clusters (times symmetry permutations).
    /// Fit failures name the member being placed.
    pub fn cluster_by_radius_with_options<'a>(
        &self,
        radius: Quantity<f64>,
        atoms: impl Into<FitAtoms<'a>>,
        options: ClusteringOptions<'_>,
    ) -> Result<EnsemblePartition, ClusteringError> {
        let radius = radius.value_in(CANONICAL_LENGTH_UNIT)?;
        if !(radius.is_finite() && radius >= 0.0) {
            return Err(ClusteringError::InvalidRadius);
        }
        let order = match options.order {
            ClusterOrder::Weight => self.ranked_by_weight(),
            ClusterOrder::Energy => self.ranked_by_energy()?,
        };
        let views = self
            .iter()
            .map(|member| member.as_model_view())
            .collect::<Vec<ModelView<'_>>>();
        let Some(first) = views.first() else {
            return Ok(EnsemblePartition {
                assignments: Vec::new(),
                representatives: Vec::new(),
            });
        };
        let atoms = atoms.into();
        let sets = atoms.pair_sets(*first, *first)?;
        check_fit_count(atoms.len())?;
        let weights = NormalizedWeights::new(options.alignment.weighting, atoms.len())?;
        let policy = options.alignment.periodic_policy;

        let mut assignments = vec![0; views.len()];
        let mut representatives: Vec<usize> = Vec::new();
        for member in order {
            let mut nearest: Option<(usize, f64)> = None;
            for (cluster, representative) in representatives.iter().enumerate() {
                let rmsd = best_kabsch(
                    views[member],
                    views[*representative],
                    &sets,
                    weights,
                    policy,
                )
                .map_err(|source| AlignmentError::Item {
                    index: member,
                    source: Box::new(source),
                })?
                .rmsd()
                .into_value();
                if rmsd <= radius && nearest.is_none_or(|(_, best)| rmsd < best) {
                    nearest = Some((cluster, rmsd));
                }
            }
            assignments[member] = match nearest {
                Some((cluster, _)) => cluster,
                None => {
                    representatives.push(member);
                    representatives.len() - 1
                }
            };
        }
        Ok(EnsemblePartition {
            assignments,
            representatives,
        })
    }

    /// One member per cluster of `partition`, in cluster order: a copy of
    /// the representative carrying the summed weight of its cluster.
    ///
    /// The representative keeps its own positions, energy, and other state,
    /// and the topology and ensemble properties are shared. Summed weights
    /// stay relative, so collapsing Boltzmann-weighted rotamers gives each
    /// conformer its degeneracy-weighted population. Fails when the
    /// partition has a different member count or a summed weight overflows.
    pub fn collapse(&self, partition: &EnsemblePartition) -> Result<Self, ClusteringError> {
        if partition.member_count() != self.len() {
            return Err(ClusteringError::PartitionMismatch {
                expected: self.len(),
                actual: partition.member_count(),
            });
        }
        let mut totals = vec![0.0; partition.cluster_count()];
        for (member, cluster) in self.iter().zip(partition.assignments()) {
            totals[*cluster] += member.weight();
        }
        let mut collapsed = self.select(partition.representatives().iter().copied())?;
        for (cluster, total) in totals.into_iter().enumerate() {
            collapsed
                .get_mut(cluster)
                .expect("one member per cluster")
                .set_weight(total)
                .map_err(EnsembleError::from)?;
        }
        Ok(collapsed)
    }
}

/// An invalid [`EnsemblePartition`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PartitionError {
    /// A member is assigned to a cluster without a representative.
    ClusterOutOfRange {
        member: usize,
        cluster: usize,
        clusters: usize,
    },
    /// A representative is not a partitioned member.
    RepresentativeOutOfRange {
        cluster: usize,
        member: usize,
        members: usize,
    },
    /// A representative is assigned to another cluster.
    RepresentativeNotInCluster { cluster: usize, member: usize },
}

impl fmt::Display for PartitionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ClusterOutOfRange {
                member,
                cluster,
                clusters,
            } => write!(
                formatter,
                "member {member} is assigned to cluster {cluster}, but there are {clusters} clusters"
            ),
            Self::RepresentativeOutOfRange {
                cluster,
                member,
                members,
            } => write!(
                formatter,
                "cluster {cluster} is represented by member {member}, but there are {members} members"
            ),
            Self::RepresentativeNotInCluster { cluster, member } => write!(
                formatter,
                "member {member} represents cluster {cluster} but is assigned to another"
            ),
        }
    }
}

impl std::error::Error for PartitionError {}

/// Failure of an ensemble clustering operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ClusteringError {
    /// A clustering radius is not finite and nonnegative.
    InvalidRadius,
    /// A partition has a different member count than the ensemble.
    PartitionMismatch {
        expected: usize,
        actual: usize,
    },
    /// A distance matrix covers a different member count than the ensemble.
    MatrixMismatch {
        expected: usize,
        actual: usize,
    },
    /// Condensed distances have the wrong length for their member count.
    CondensedLength {
        members: usize,
        expected: usize,
        actual: usize,
    },
    /// A distance is not finite and nonnegative; `index` is condensed.
    InvalidDistance {
        index: usize,
    },
    /// The pairwise distances of this many members are not addressable.
    MatrixTooLarge {
        members: usize,
    },
    /// Fitting one pair of members failed.
    PairFit {
        first: usize,
        second: usize,
        source: Box<AlignmentError>,
    },
    /// A dendrogram cut height is not finite and nonnegative.
    InvalidHeight,
    /// A requested cluster count is zero or exceeds the member count.
    InvalidClusterCount {
        requested: usize,
        members: usize,
    },
    Alignment(Box<AlignmentError>),
    Ensemble(Box<EnsembleError>),
    Unit(Box<UnitError>),
}

impl fmt::Display for ClusteringError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRadius => {
                formatter.write_str("clustering radius must be finite and nonnegative")
            }
            Self::PartitionMismatch { expected, actual } => write!(
                formatter,
                "partition covers {actual} members, but the ensemble has {expected}"
            ),
            Self::MatrixMismatch { expected, actual } => write!(
                formatter,
                "distance matrix covers {actual} members, but the ensemble has {expected}"
            ),
            Self::CondensedLength {
                members,
                expected,
                actual,
            } => write!(
                formatter,
                "{members} members require {expected} condensed distances, but received {actual}"
            ),
            Self::InvalidDistance { index } => write!(
                formatter,
                "condensed distance {index} must be finite and nonnegative"
            ),
            Self::MatrixTooLarge { members } => write!(
                formatter,
                "pairwise distances of {members} members are not addressable"
            ),
            Self::PairFit {
                first,
                second,
                source,
            } => write!(
                formatter,
                "cannot fit members {first} and {second}: {source}"
            ),
            Self::InvalidHeight => {
                formatter.write_str("dendrogram cut height must be finite and nonnegative")
            }
            Self::InvalidClusterCount { requested, members } => write!(
                formatter,
                "cannot cut {members} members into {requested} clusters"
            ),
            Self::Alignment(error) => write!(formatter, "cannot fit members: {error}"),
            Self::Ensemble(error) => error.fmt(formatter),
            Self::Unit(error) => write!(formatter, "invalid quantity unit: {error}"),
        }
    }
}

impl std::error::Error for ClusteringError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Alignment(error) | Self::PairFit { source: error, .. } => Some(error.as_ref()),
            Self::Ensemble(error) => Some(error.as_ref()),
            Self::Unit(error) => Some(error.as_ref()),
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
