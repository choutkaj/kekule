#![forbid(unsafe_code)]

//! DREIDING force-field preparation and evaluation for topology-bound
//! structural views.
//!
//! This adapter keeps automatic force-field preparation outside the lightweight
//! `kekule` core crate. Preparation is explicit: it never sanitizes input,
//! adds hydrogens, changes topology, or updates charges during evaluation.
//!
//! # Periodic-cell capability
//!
//! The current adapter is explicitly nonperiodic. Preparation rejects a
//! reference view with a periodic cell using
//! [`DreidingPrepareError::UnsupportedPeriodicCell`], and evaluation rejects
//! periodic views using
//! [`kekule::modeling::potential::PotentialError::UnsupportedPeriodicCell`].
//! No cell is silently ignored, and no orthorhombic-only minimum-image shortcut
//! is applied.
//!
mod error;
mod evaluate;
mod geometry;
mod prepare;

pub use error::DreidingPrepareError;
pub use prepare::{DreidingPotential, DreidingPrepareOptions, QeqGrouping};

#[cfg(test)]
mod tests;
