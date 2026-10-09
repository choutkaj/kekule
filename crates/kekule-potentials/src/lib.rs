#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]

//! Energies, gradients, and geometry optimization for topology-bound Kekule
//! models.
//!
//! The crate separates three layers:
//!
//! - [`Potential`] is the evaluation contract. It promises an energy and a
//!   Cartesian gradient for any [`kekule::structure::ModelView`] sharing one
//!   topology layout, including snapshots perceived after preparation, and
//!   assumes nothing about how the energy is built.
//! - Backends prepare a potential from explicit parameters. The default
//!   `openff` feature provides `openff::OpenFfPotential`, which evaluates a
//!   `kekule_openff::ParameterizedTopology`. Classical backends share private
//!   functional-form kernels that follow the singular-geometry policy below.
//! - [`minimize()`] optimizes coordinates under any potential.
//!
//! Preparation is explicit: nothing here parameterizes, perceives, adds
//! hydrogens, or otherwise changes chemistry. All values use Kekule's canonical
//! units: nm, kJ/mol, kJ/mol/nm, radians, and elementary charges.
//!
//! # Singular-geometry policy
//!
//! Gradients are exact derivatives of the evaluated energy and are never capped,
//! so near-singular coordinates produce large but finite gradients. Coordinates
//! are rejected with [`EvaluationError::InvalidGeometry`] only when a requested
//! quantity is mathematically undefined: nonbonded atoms closer than `1e-12` nm;
//! a zero-length angle arm (or one whose squared length underflows); a torsion with an axis shorter than `1e-12` nm or
//! with outer atoms whose perpendicular distances from the axis multiply to less
//! than `1e-18` nm^2; and, when a gradient is requested, a bond shorter than
//! `1e-12` nm or an exactly linear angle away from its equilibrium value.
//! [`Potential::energy`] accepts the last two, whose energies are defined.
//!
//! # Typical workflow
//!
//! ```no_run
//! # #[cfg(feature = "openff")]
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! use kekule::structure::Model;
//! use kekule_openff::{ForceField, NaglModel};
//! use kekule_potentials::{minimize, openff::OpenFfPotential, MinimizeOptions, Potential};
//!
//! # fn explicit_hydrogen_model() -> Model { unimplemented!() }
//! let model = explicit_hydrogen_model();
//! let nagl = NaglModel::ash()?;
//! let parameters = ForceField::rosemary()?.parameterize(model.shared_topology(), &nagl)?;
//! let potential = OpenFfPotential::new(&parameters)?;
//!
//! let energy = potential.energy(model.as_model_view())?;
//! println!("single point: {}", energy.total().into_value());
//!
//! let result = minimize(&potential, model.as_model_view(), &MinimizeOptions::default())?;
//! let minimized = result.to_model(model.as_model_view())?;
//! # let _ = minimized;
//! # Ok(())
//! # }
//! # #[cfg(not(feature = "openff"))]
//! # fn main() {}
//! ```

mod minimize;
mod mm;
mod potential;

#[cfg(feature = "openff")]
pub mod openff;

pub use minimize::{
    minimize, minimize_with_observer, Minimization, MinimizationError, MinimizationStatus,
    MinimizationStep, MinimizeOptions,
};
pub use potential::{
    ComponentKind, Energy, EnergyComponent, Evaluation, EvaluationError, Potential,
    SingularGeometry,
};
