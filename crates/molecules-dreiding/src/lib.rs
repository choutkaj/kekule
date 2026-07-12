#![forbid(unsafe_code)]

//! DREIDING force-field preparation and evaluation for [`MolecularModel`].
//!
//! This adapter keeps automatic force-field preparation outside the lightweight
//! `molecules` core crate. Preparation is explicit: it never sanitizes input,
//! adds hydrogens, changes topology, or updates charges during evaluation.
//!
//! [`MolecularModel`]: molecules::modeling::MolecularModel

mod error;
mod evaluate;
mod geometry;
mod prepare;

pub use error::DreidingPrepareError;
pub use prepare::DreidingPotential;

#[cfg(test)]
mod tests;
