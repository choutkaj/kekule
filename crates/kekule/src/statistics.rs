//! Weighted structural statistics: per-atom RMS fluctuation and contact
//! occupancy of atom pairs or atom groups over weighted structures.
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
//! Ensembles and structure sets reduce themselves with the same
//! accumulators. [`Ensemble::rmsf`], [`Ensemble::contact_occupancy`], and
//! [`Ensemble::group_contact_occupancy`] observe the samples of the
//! ensemble's distribution, each with its probability, and leave the other
//! geometries out. The same methods of a [`StructureSet`] observe every
//! geometry with weight one, taking a trajectory's frames as equally
//! weighted equilibrium samples; to reweight frames, for example after
//! biased sampling, read the set as an [`Ensemble::weighted`] or observe the
//! frames through an accumulator with explicit weights.
//!
//! An accumulator also observes the geometries of any
//! [`StructureSetRead`], such as a set read from a
//! file: [`SampleObservations`] lists the geometries and their weights, and
//! `observe_samples` reads a range of them in order, as the statistics of a
//! set or an ensemble do.
//!
//! Coordinates are used as stored, ignoring cells. Align and perform any
//! periodic preprocessing explicitly before observation. Lengths are reported
//! in the canonical length unit (nm).
//!
//! # Examples
//!
//! ```
//! use kekule::math::Point3;
//! use kekule::structure::{Ensemble, Positions, Geometry, StructureSet};
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
//! let ensemble = Ensemble::weighted(set, vec![3.0_f64.ln(), 0.0])?;
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
//! // ln 3 is not exact in binary, so the share is 3/4 up to rounding.
//! let ((first, second), occupancy) = contacts.contacts().next().unwrap();
//! assert_eq!((first, second), (atoms[0], atoms[1]));
//! assert!((occupancy - 0.75).abs() < 1e-15);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod accumulators;

pub use accumulators::{
    ContactOccupancyAccumulator, ContactOccupancyResult, GroupContactOccupancyAccumulator,
    GroupContactOccupancyResult, MergeMismatch, RmsfAccumulator, RmsfResult, SampleObservations,
    StatisticsError,
};

#[cfg(doc)]
use crate::structure::{AsStructureView, Ensemble, StructureSet, StructureSetRead};
