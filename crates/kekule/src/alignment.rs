//! Rigid alignment and RMSD through selections or explicit atom correspondence.
//!
//! Atoms are paired by [`FitAtoms`]: an [`AtomSelection`] pairs each selected
//! atom with itself (both sides share one topology layout), an
//! [`AtomCorrespondence`] pairs atoms of independent topologies (given
//! explicitly, by dense order, or by hierarchy labels with
//! [`AtomCorrespondence::from_labels`]), and a
//! [`SymmetricSelection`] pairs a selection under every permutation of
//! symmetry-equivalent atoms, keeping the smallest RMSD. Every operation
//! accepts any of them through `impl Into<FitAtoms>`.
//!
//! - [`fit`] fits one view onto another with the Kabsch algorithm and returns
//!   the moving-to-reference transform and post-fit RMSD.
//! - [`rmsd`] measures coordinates exactly as stored, without fitting.
//! - [`aligned_rmsd`] fits on one set of pairs and measures another.
//! - `superpose`, `rmsd`, and `aligned_rmsd` on a
//!   [`StructureSet`](crate::structure::StructureSet::superpose) apply the same
//!   kernels to every geometry against a [`Reference`]: a geometry of the set
//!   or any external view.
//! - [`geometry_rmsd`] and [`geometry_aligned_rmsd`] measure listed geometries of
//!   any [`StructureSetRead`](crate::structure::StructureSetRead), such as a set read
//!   from a file, one at a time.
//!
//! Each takes its options, [`AlignmentOptions`] or [`AlignedRmsdOptions`], as
//! a parameter. The returned transform maps moving coordinates into reference
//! coordinates:
//!
//! ```text
//! aligned = result.transform().transform_point(canonical_moving)
//! ```
//!
//! Alignment never images periodic coordinates. The default fit uses stored
//! Cartesian coordinates, including when a cell is present; use
//! [`PeriodicPolicy::RejectPeriodic`] to reject periodic inputs.
//!
//! # Examples
//!
//! ```
//! use kekule::alignment::fit;
//! use kekule::math::Point3;
//! use kekule::smiles;
//! use kekule::structure::{Structure, Positions};
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, ANGSTROM};
//!
//! let topology = smiles::to_topology("CCC")?;
//! let moving_points = [
//!     Point3::new(0.0, 0.0, 0.0),
//!     Point3::new(1.0, 0.0, 0.0),
//!     Point3::new(0.0, 1.0, 0.0),
//! ];
//! let reference_points = moving_points.map(|point| Point3::new(
//!     point.x + 2.0,
//!     point.y - 1.0,
//!     point.z + 0.5,
//! ));
//! let moving = Structure::new(
//!     topology.clone(),
//!     Positions::new(Quantity::new(moving_points, ANGSTROM))?,
//! )?;
//! let reference = Structure::new(
//!     topology.clone(),
//!     Positions::new(Quantity::new(reference_points, ANGSTROM))?,
//! )?;
//! let selection = AtomSelection::all(topology.clone());
//!
//! let result = fit(
//!     moving.as_structure_view(),
//!     reference.as_structure_view(),
//!     &selection,
//!     Default::default(),
//! )?;
//! let canonical_moving = moving.positions().get(1).expect("atom 1").into_value();
//! let canonical_reference = reference.positions().get(1).expect("atom 1").into_value();
//! let aligned = result.transform().transform_point(canonical_moving);
//! assert!((aligned.x - canonical_reference.x).abs() < 1.0e-12);
//! assert!(result.rmsd().into_value() < 1.0e-12);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod collection;
mod correspondence;
mod rigid;
mod symmetry;

#[cfg(test)]
mod correspondence_label_tests;
#[cfg(test)]
mod correspondence_tests;
#[cfg(test)]
mod symmetry_tests;

pub use collection::{
    geometry_aligned_rmsd, geometry_aligned_rmsd_with, geometry_rmsd, geometry_rmsd_with,
    AlignedRmsdOptions, Reference, SuperpositionReport,
};
pub use correspondence::{
    AmbiguousLabel, AtomCorrespondence, AtomCorrespondenceError, AtomLabel, CorrespondenceSide,
    LabelCorrespondence, LabelScheme,
};
pub use rigid::{
    aligned_rmsd, fit, rmsd, AlignmentError, AlignmentGeometry, AlignmentOptions, FitAtoms,
    PeriodicPolicy, RigidAlignment, Weighting,
};
pub(crate) use rigid::{
    best_kabsch, best_measure, check_fit_count, check_measurement_count, Coordinates,
    NormalizedWeights, Pairs,
};
pub use symmetry::{SymmetricSelection, SymmetryError, SymmetryLimits};

#[cfg(doc)]
use crate::topology::AtomSelection;
