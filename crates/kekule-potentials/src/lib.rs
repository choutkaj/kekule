#![forbid(unsafe_code)]

//! Concrete prepared potentials for topology-bound Kekule structural views.
//!
//! The lightweight evaluation contract lives in
//! [`kekule::modeling::potential`]. This companion crate owns force-field
//! preparation and parameters so the foundational `kekule` crate remains
//! independent of any particular potential model.
//!
//! The default `dreiding` feature provides [`dreiding::DreidingPotential`].
//! Preparation is explicit and bound to one exact shared topology. It does not
//! parse, perceive, add hydrogens, or otherwise alter input chemistry.
//!
//! # Typical workflow
//!
//! ```no_run
//! use kekule::{modeling::potential::Potential, structure::Model};
//! use kekule_potentials::dreiding::{
//!     DreidingPotential, DreidingPrepareOptions,
//! };
//!
//! # fn prepared_model() -> Model { unimplemented!() }
//! let model = prepared_model();
//! let topology = model.shared_topology();
//! let mut potential = DreidingPotential::prepare(
//!     &topology,
//!     model.view(),
//!     DreidingPrepareOptions::default(),
//! )?;
//! let evaluation = potential.evaluate(model.view())?;
//! assert!(evaluation.energy().is_finite());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#[cfg(feature = "dreiding")]
pub mod dreiding;
