//! Coordinate-dependent modelling contracts and minimization.
//!
//! [`potential::Potential`] is the lightweight evaluation interface consumed by
//! [`minimize`]. Concrete force-field preparation belongs in companion crates
//! such as `kekule-potentials`; the foundational crate does not select,
//! parameterize, or mutate chemistry implicitly.

mod minimize;
pub mod potential;

pub use minimize::*;

#[cfg(test)]
mod tests;
