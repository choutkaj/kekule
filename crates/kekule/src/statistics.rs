//! Weights of a structure set's geometries and weighted structural
//! statistics: [`Weights`], per-atom RMS fluctuation and contact occupancy
//! of atom pairs or atom groups over weighted structures, and free-energy
//! and mean landscapes.
//!
//! # Weights
//!
//! [`Weights`] give each geometry of a set an optional natural-logarithm
//! weight: which geometries are samples and what each weighs. They are a
//! value beside the set, indexed by geometry, never a field of it, so one
//! set carries several readings. They come from a trajectory's equal frames
//! ([`Weights::equal`]), an energy column at a temperature
//! ([`Weights::boltzmann`]), the sampling record of a biased run
//! ([`Weights::from_bias`]), or a host's own log weights
//! ([`Weights::given`]). Every population is exact, the log-sum-exp of its
//! samples' log weights, so the populations of a group's parts add up to
//! the group's. A [`Clustering`](crate::clustering::Clustering) reads its
//! clusters' populations from them.
//!
//! # Accumulators
//!
//! An observation is one structure with an explicit, finite, positive weight.
//! Weights are relative: scaling every weight by one factor changes no
//! statistic. [`RmsfAccumulator`], [`ContactOccupancyAccumulator`], and
//! [`GroupContactOccupancyAccumulator`] observe any borrowed structure
//! ([`AsStructureView`]: a structure, a geometry of a set, or a codec frame such as a streamed
//! frame buffer's `frame_view()`) one at a time, with memory bounded by the
//! selected atoms, pairs, or groups, not the observation count.
//!
//! Accumulators merge, so a host can observe in parallel: it splits the
//! observations into ranges, observes each range into its own accumulator,
//! and merges the partial accumulators in range order. For a given split the
//! result is the same whichever threads observed the ranges, and it equals
//! observing everything in one accumulator up to rounding.
//!
//! Structure sets reduce themselves with the same accumulators.
//! [`StructureSet::rmsf`], [`StructureSet::contact_occupancy`], and
//! [`StructureSet::group_contact_occupancy`] observe the samples of
//! weights, each with its probability, and leave the other geometries out.
//! Without weights they observe every geometry with weight one, taking a
//! trajectory's frames as equally weighted equilibrium samples.
//!
//! An accumulator also observes the geometries of any
//! [`StructureSetRead`], such as a set read from a
//! file: [`SampleObservations`] lists the geometries and their weights, and
//! `observe_samples` reads a range of them in order, as the statistics of a
//! set do.
//!
//! Coordinates are used as stored, ignoring cells. Align and perform any
//! periodic preprocessing explicitly before observation. Lengths are reported
//! in the canonical length unit (nm).
//!
//! # Examples
//!
//! ```
//! use kekule::math::Point3;
//! use kekule::statistics::Weights;
//! use kekule::structure::{Positions, Geometry, StructureSet};
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, NANOMETER};
//!
//! let topology = kekule::smiles::to_topology("CC")?;
//! let mut set = StructureSet::new(topology.clone());
//! // The first atom sits at x = 0 with weight 3 and at x = 1 with weight 1.
//! for x in [0.0, 1.0] {
//!     let points = [Point3::new(x, 0.0, 0.0), Point3::new(0.15, 0.0, 0.0)];
//!     set.push(Geometry::new(Positions::new(Quantity::new(points, NANOMETER))?), None)?;
//! }
//! let weights = Weights::given(vec![Some(3.0_f64.ln()), Some(0.0)], None)?;
//!
//! // Weighted mean 0.25, weighted variance 0.75 × 0.25² + 0.25 × 0.75² = 0.1875.
//! let fluctuation = set.rmsf(&AtomSelection::all(topology.clone()), Some(&weights))?;
//! assert!((fluctuation.values().value()[0] - 0.1875_f64.sqrt()).abs() < 1e-12);
//! assert_eq!(fluctuation.values().value()[1], 0.0);
//!
//! // The pair is within 0.2 nm only in the first sample, carrying 3/4 of the weight.
//! let atoms = topology.atom_ids();
//! let cutoff = Quantity::new(0.2, NANOMETER);
//! let contacts = set.contact_occupancy([(atoms[0], atoms[1])], cutoff, Some(&weights))?;
//! // ln 3 is not exact in binary, so the share is 3/4 up to rounding.
//! let ((first, second), occupancy) = contacts.contacts().next().unwrap();
//! assert_eq!((first, second), (atoms[0], atoms[1]));
//! assert!((occupancy - 0.75).abs() < 1e-15);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Energy surfaces
//!
//! A [`Landscape`] bins the samples of weights on one or two [`GridAxis`]es
//! and gives each cell a free energy or a mean value. Its cell populations
//! are exact unless the options set another [`Estimator`], and its means
//! weigh each sample by its exact share in the cell.
//!
//! An [`Estimator`] estimates the population `Σ exp(b_i)` of a cell from
//! dimensionless reduced values `b_i`: the log importance weights of its
//! samples, for example a bias potential over the thermal energy.
//! [`Estimator::log_sum`] returns its logarithm, so populations spanning
//! hundreds of orders of magnitude stay finite. The cumulant expansion
//! lowers the variance of a sparse cell's estimate at the cost of bias.

mod accumulators;
mod estimators;
mod landscape;
mod weights;

pub(crate) use accumulators::{merged_total, next_total, observe_samples_of};
pub use accumulators::{
    ContactOccupancyAccumulator, ContactOccupancyResult, GroupContactOccupancyAccumulator,
    GroupContactOccupancyResult, MergeMismatch, RmsfAccumulator, RmsfResult, SampleObservations,
    StatisticsError,
};
pub(crate) use estimators::log_sum_exp;
pub use estimators::{CumulantOrder, Estimator};
pub use landscape::{
    Binning, Cell, GridAxis, Landscape, LandscapeError, LandscapeOptions, SurfaceQuantity,
};
pub(crate) use weights::{canonical, check_fraction, covered_count, kish};
pub use weights::{BiasWeights, LeftOut, Unweighable, WeightSource, Weights, WeightsError};

#[cfg(doc)]
use crate::structure::{AsStructureView, StructureSet, StructureSetRead};
