//! Chemical component templates: residue and ligand definitions in the model of
//! the wwPDB Chemical Component Dictionary.
//!
//! Structural formats name each residue's component but rarely carry its bonds.
//! [`ComponentLibrary::bundled`](crate::components::ComponentLibrary::bundled) holds
//! the standard amino acids and nucleotides, water, common ions, cofactors, and
//! crystallization additives; read further definitions, up to the complete
//! dictionary, with
//! [`ComponentLibrary::from_mmcif`](crate::components::ComponentLibrary::from_mmcif).
//! Structures that molecular-mechanics tools write use force-field residue and
//! atom names; [`ComponentLibrary::amber`](crate::components::ComponentLibrary::amber)
//! and [`ComponentLibrary::charmm`](crate::components::ComponentLibrary::charmm)
//! define them.
//!
//! ```
//! use kekule::components::ComponentLibrary;
//!
//! let alanine = ComponentLibrary::bundled().get("ALA").expect("bundled");
//! assert_eq!(alanine.name(), "ALANINE");
//! assert_eq!(alanine.atoms().len(), 13);
//! assert!(alanine.atom("CA").is_some());
//! assert!(ComponentLibrary::bundled().get("XYZ").is_none());
//! ```

pub use crate::io::{
    Component, ComponentAtom, ComponentBond, ComponentBondOrder, ComponentLibrary, ComponentLinking,
};
