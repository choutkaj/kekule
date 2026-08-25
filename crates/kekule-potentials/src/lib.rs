#![forbid(unsafe_code)]

//! Concrete potential implementations for topology-bound Kekule structural
//! views.
//!
//! The dependency-light evaluation contract remains in
//! [`kekule::modeling::potential`]. Implementations in this companion crate are
//! organized by preparation model so additional potentials can be added
//! without widening the foundational `kekule` crate.
//!
//! Import implementations through their focused module:
//!
//! ```
//! use kekule_potentials::dreiding::{DreidingPotential, DreidingPrepareOptions};
//! ```

#[cfg(feature = "dreiding")]
pub mod dreiding;
