//! Explicit small-molecule hydrogen topology transforms.
//!
//! Explicit hydrogens are graph atoms; an atom's implicit hydrogens are either
//! a fixed count or wholly inferred by valence perception. Addition converts
//! implicit H to explicit atoms; removal suppresses eligible explicit atoms
//! without losing chemical information. These operations do not change
//! protonation states.
//!
//! Addition requires current inference for inferred counts, unless
//! [`AddHydrogensScope::FixedOnly`] is selected. Fixed counts need no perception. Removal needs current inference
//! on affected parents whose counts are inferred, and keeps a parent inferred
//! only when inference reproduces its collapsed count; otherwise it fixes it.
//! Successful topology changes invalidate perception state; recompute it before
//! reading inferred counts.
//!
//! ```
//! use kekule::smiles;
//!
//! let mut ethanol = smiles::to_molecules("CCO")?.remove(0);
//! ethanol.perceive()?;
//! let added = ethanol.add_hydrogens(Default::default())?;
//! assert_eq!(added.added().len(), 6);
//! assert_eq!(ethanol.atom_count(), 9);
//!
//! ethanol.perceive()?;
//! let removed = ethanol.remove_hydrogens()?;
//! assert_eq!(removed.removed().len(), 6);
//! assert_eq!(ethanol.atom_count(), 3);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::algorithms::{
    AddHydrogensOptions, AddHydrogensReport, AddHydrogensScope, AddedHydrogen, AddedHydrogenOrigin,
    HydrogenCountAdjustment, HydrogenTransformError, RemoveHydrogensReport, RemovedHydrogen,
    RetainedHydrogen, RetainedHydrogenReason,
};
