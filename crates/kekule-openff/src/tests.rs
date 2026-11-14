//! Tests that compare model internals, which are not public, with reference
//! data: lookup keys, input features, and network charges.

#[cfg(feature = "ash")]
mod bundles;
#[cfg(feature = "ash")]
mod charges;
mod models;
