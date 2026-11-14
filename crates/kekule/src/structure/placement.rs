//! Builder geometry: coordinates for the atoms an editor adds.
//!
//! Kekule never invents geometry implicitly; placement is an operation a host
//! asks for, on a [`StructureEditor`] draft, as part of the edit it publishes:
//!
//! - [`StructureEditor::add_hydrogens`] makes the implicit hydrogens of chosen
//!   atoms explicit at ideal positions, and [`StructureEditor::remove_hydrogens`]
//!   collapses their explicit hydrogens back into implicit counts, each
//!   reporting per atom what changed ([`PlacedHydrogensReport`],
//!   [`RemovedHydrogensReport`]).
//! - [`StructureEditor::atom_domains`] perceives an atom's [`AtomDomains`] and
//!   [`StructureEditor::free_directions`] gives where its next substituents go.
//! - [`StructureEditor::attach_atom`] bonds a new atom to an existing one, in place
//!   of an explicit hydrogen or at the free position nearest a direction
//!   ([`AttachmentSite`]); [`StructureEditor::place_atom`] starts a new molecule
//!   at a point.
//! - [`StructureEditor::attach_fragment`] attaches a named [`Fragment`] of the
//!   library ([`fragments`], [`fragment`]) the same ways, turned to the
//!   clearest of a few staggered torsions about the new bond, and
//!   [`StructureEditor::place_fragment`] places one as a new capped molecule at a
//!   point.
//! - [`ideal_bond_length`] gives the length placement uses for an element pair
//!   and bond order.
//!
//! # Ideal geometry
//!
//! An atom's shape is its number of electron domains: its neighbors, its
//! implicit hydrogens, and the lone pairs that its element's valence
//! electrons, formal charge, radical electrons, and bond orders leave
//! ([`AtomShape`]). A neutral three-coordinate nitrogen with single bonds to a
//! boron, carbon, nitrogen, or oxygen that has a multiple bond (an amide,
//! aniline, or pyrrole nitrogen) is conjugated and planar. Elements outside
//! the s and p blocks have no lone pairs here, so their shape follows their
//! neighbors. Up to six domains are supported.
//!
//! New sites keep the ideal angles of the shape toward the atom's existing
//! neighbors. With one neighbor, the first new site is anti to that
//! neighbor's main substituent, a heavy atom preferring the one with the most
//! heavy neighbors, so chains grow anti and substituents stay staggered;
//! trigonal atoms stay in the plane of that substituent, as across a double
//! bond. Lone pairs of a trigonal bipyramid are equatorial; those of an
//! octahedron are trans to each other or to a neighbor. Placement never moves
//! an existing atom; relaxing the result is a separate step.
//!
//! # Names
//!
//! An atom a placement adds joins the residue of the atom it bonds to, named
//! by its element's uppercase symbol and the next free number in that
//! residue, as chemical components name theirs (`C21`, `H31`): one more than
//! the largest number following that symbol in a name of the residue. An
//! atom without a residue gives its new neighbors none. A new molecule gets a
//! chain and a residue `UNL` of its own; a fragment placed alone names its
//! ring atoms by locant (pyridine `N1`, `C2`, …; naphthalene `C4A`) and its
//! hydrogens after the atoms carrying them (`H2`; `H21`, `H22`), as
//! [`Fragment`] describes. [`StructureEditor::next_atom_name`] and
//! [`StructureEditor::rename_atom`] apply the same rule to existing atoms. Names
//! are label and author atom names alike, and the host may rename.
//!
//! ```
//! use kekule::{smiles, structure::{Structure, Positions}};
//! use kekule::hydrogens::AddHydrogensOptions;
//! use kekule::units::{Quantity, NANOMETER};
//! use kekule::math::Point3;
//!
//! // Ethanol's heavy atoms, then its hydrogens placed in the same draft.
//! let topology = smiles::to_topology("CCO")?;
//! let positions = Positions::new(Quantity::new(
//!     [Point3::origin(), Point3::new(0.153, 0.0, 0.0), Point3::new(0.2, 0.135, 0.0)],
//!     NANOMETER,
//! ))?;
//! let mut editor = Structure::new(topology, positions)?.into_editor();
//! let atoms = editor.topology().atom_ids().collect::<Vec<_>>();
//! let report = editor.add_hydrogens(atoms, AddHydrogensOptions::default())?;
//! assert_eq!(report.added_count(), 6);
//! assert_eq!(editor.finish()?.atom_count(), 9);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod domains;
mod fragments;
mod lengths;
mod names;
mod reports;
mod sites;

pub use domains::{ideal_bond_length, AtomDomains, AtomShape};
pub(crate) use domains::{perceive_domains, DomainInput};
pub use fragments::{fragment, fragments, Fragment, FragmentKind};
pub(crate) use names::{fragment_names, name_prefix, next_free_number, NameAllocator};
pub use reports::{
    AttachedAtom, AttachedFragment, AttachmentSite, PlacedFragment, PlacedHydrogens,
    PlacedHydrogensReport, PlacementError, RemovedHydrogens, RemovedHydrogensReport,
};
pub(crate) use sites::{free_sites, unit, Rotation, SiteQuery};

#[cfg(doc)]
use crate::structure::StructureEditor;
