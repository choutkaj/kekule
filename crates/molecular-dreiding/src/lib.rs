#![forbid(unsafe_code)]

//! DREIDING force-field preparation and evaluation for [`Model`].
//!
//! This adapter keeps automatic force-field preparation outside the lightweight
//! `molecular` core crate. Preparation is explicit: it never sanitizes input,
//! adds hydrogens, changes topology, or updates charges during evaluation.
//!
//! [`Model`]: molecular::modeling::Model

mod error;
mod evaluate;
mod geometry;
mod prepare;

pub use error::DreidingPrepareError;
pub use prepare::DreidingPotential;

#[cfg(test)]
mod tests;
