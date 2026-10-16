//! Periodic box construction and force-field-independent TIP3P solvation.
//!
//! These operations create starting structures, not parameterized or equilibrated
//! simulation systems. Existing coordinates, chemistry, and entity annotations are
//! preserved. Solvation publishes one new topology layout; old selections and
//! parameter assignments remain bound to the original one. Existing owner properties
//! are cleared under the ordinary append contract and listed in the optional report.
//!
//! [`add_periodic_box`] sets a fully periodic cell on a
//! [`Model`](crate::structure::Model), and [`solvate`] fills that cell with
//! water and ions.

mod solvent;

pub use solvent::{
    add_periodic_box, solvate, BoxShape, NegativeIon, PeriodicBoxError, PeriodicBoxOptions,
    PositiveIon, SolvationChargeBasis, SolvationError, SolvationRadii, SolvationReport,
    SolventOptions, WaterModel,
};
