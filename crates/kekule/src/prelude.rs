//! Common foundational types for small examples and interactive use.
//!
//! The prelude is intentionally small: the molecule, its editor and parts, the
//! topology and structure that hold it, and
//! [`Quantity`]. It is the only module that repeats
//! items exported elsewhere. Format APIs, identifiers, and algorithms remain
//! available through their focused modules.
//!
//! ```
//! use kekule::prelude::*;
//! use kekule::{math::Point3, units::ANGSTROM};
//!
//! let mut editor = MoleculeEditor::new();
//! let carbon = editor.add_atom(Atom::new(Element::from_symbol("C").unwrap()))?;
//! let oxygen = Atom::new(Element::O).with_hydrogens(ImplicitHydrogens::Fixed(1));
//! let oxygen = editor.add_atom(oxygen)?;
//! editor.add_bond(carbon, oxygen, BondOrder::Single)?;
//! let molecule: Molecule = editor.finish()?;
//! let (_, bond): (_, &Bond) = molecule.bonds().next().unwrap();
//! assert_eq!(bond.order, BondOrder::Single);
//!
//! let topology = Topology::from_molecule(molecule.clone())?;
//! let positions = Positions::new(Quantity::new(
//!     vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.43, 0.0, 0.0)],
//!     ANGSTROM,
//! ))?;
//! let structure = Structure::from_molecule(molecule, &positions)?;
//! assert_eq!(structure.atom_count(), topology.atom_count());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::molecule::{
    Atom, Bond, BondOrder, Element, ImplicitHydrogens, Molecule, MoleculeEditor,
};
pub use crate::structure::{Positions, Structure};
pub use crate::topology::Topology;
pub use crate::units::Quantity;
