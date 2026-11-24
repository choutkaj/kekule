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
//! - [`StructureEditor::settle_atoms`] gives atoms that an edit changed their
//!   builder hydrogen count ([`builder_hydrogen_count`]) at ideal positions and
//!   the angle of their shape between two neighbors, and
//!   [`StructureEditor::settle_bonds`] gives bonds their ideal length, each
//!   moving the smaller side ([`SettledAtomsReport`], [`SettledBondsReport`]).
//! - [`StructureEditor::atom_domains`] perceives an atom's [`AtomDomains`] and
//!   [`StructureEditor::free_directions`] gives where a new substituent bonded
//!   with a given order goes.
//! - [`StructureEditor::attach_atom`] bonds a new atom to an existing one, in place
//!   of an explicit hydrogen or at the free position nearest a direction
//!   ([`AttachmentSite`]); [`StructureEditor::place_atom`] starts a new molecule
//!   at a point.
//! - [`StructureEditor::attach_fragment`] attaches a named [`Fragment`] of the
//!   library ([`fragments`], [`fragment`]) the same ways, turned to the
//!   clearest of a few staggered torsions about the new bond, and
//!   [`StructureEditor::place_fragment`] places one as a new capped molecule at a
//!   point.
//! - [`StructureEditor::fuse_fragment`] fuses a ring fragment onto a bond
//!   through its [fusion bond](Fragment::fusion_bond), pairing the double
//!   bonds again and reporting a [`FusedFragment`].
//! - [`ideal_bond_length`] gives the length placement uses for an element pair
//!   and bond order.
//!
//! # Builder valences
//!
//! Reading a file infers hydrogens with the valence model, which counts the
//! hydrogens an atom has. A builder instead gives each atom the hydrogens its
//! element and charge take beside its other bonds, from a table of the
//! states it makes ([`builder_charges`], [`builder_valences`]):
//!
//! | Element | Charge: valences |
//! | --- | --- |
//! | H | 0: 1 |
//! | B | −1: 4; 0: 3 |
//! | C | −1: 3; 0: 4; +1: 3 |
//! | N | −1: 2; 0: 3; +1: 4 |
//! | O | −1: 1; 0: 2; +1: 3 |
//! | F, Cl, Br, I | −1: 0; 0: 1 |
//! | Si | 0: 4 |
//! | P | 0: 3, 5; +1: 4 |
//! | S, Se | −1: 1; 0: 2, 4, 6; +1: 3 |
//!
//! An atom takes the lowest valence at or above the bond orders to its other
//! neighbors plus its radical electrons, and the difference in hydrogens, so
//! a sulfide's sulfur takes none, a sulfonium's one, and a sulfoxide's none.
//! Every state is one the valence model accepts, so a settled structure
//! publishes. Elements outside the table follow the valence model.
//!
//! # Ideal geometry
//!
//! An atom's shape is its number of electron domains: its neighbors, its
//! implicit hydrogens, and the lone pairs that its element's valence
//! electrons, formal charge, radical electrons, and bond orders leave
//! ([`AtomShape`]). Two radical electrons take one domain, so a carbene or
//! nitrene is bent; a single one takes none, so a methyl radical is planar.
//! A neutral three-coordinate nitrogen with single bonds to a
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
//! octahedron are trans to each other or to a neighbor. Two opposite
//! neighbors of a trigonal or tetrahedral atom leave the plane of its new
//! sites to a deterministic perpendicular.
//!
//! Placement never moves an existing atom. Settling moves only what it
//! reports: an atom's hydrogens, and the smaller side of a bond or of an
//! atom's two neighbors, never a side that lies in a ring with them.
//! Relaxing the result is a separate step.
//!
//! # Names
//!
//! An atom a placement adds joins the residue of the atom it bonds to, named
//! by its element's uppercase symbol and the next free number in that
//! residue, as chemical components name theirs (`C21`, `H31`): one more than
//! the largest number following that symbol in a name of the residue. Names
//! that the bundled template of the residue's component places elsewhere are
//! skipped, as are the terminal hydrogen names readers bond by convention
//! (`H1` to `H3` of a peptide's `N`), so that a reader does not bond the atom
//! where its name says. A hydrogen takes the first free name the template
//! gives a hydrogen on its heavy neighbor instead (`HA` on a `CA`, `HB2` and
//! `HB3` on a `CB`), or a free terminal name that readers bond to that atom.
//! An atom without a residue gives its new neighbors none. A new molecule gets a
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
mod fusion;
mod lengths;
mod names;
mod reports;
mod sites;
mod valences;

pub use crate::topology::editor::RemovedHydrogens;
pub use domains::{ideal_bond_length, AtomDomains, AtomShape};
pub(crate) use domains::{perceive_domains, DomainInput};
pub use fragments::{fragment, fragments, Fragment, FragmentKind};
pub(crate) use fusion::{
    fusion_dropped, fusion_orders, single_heavy_neighbor, FusionEnd, FusionOrders,
};
pub(crate) use names::{fragment_names, name_prefix, residue_names, NameAllocator, ResidueNames};
pub use reports::{
    AttachedAtom, AttachedFragment, AttachmentSite, FusedFragment, PlacedFragment, PlacedHydrogens,
    PlacedHydrogensReport, PlacementError, RemovedHydrogensReport, SettledAtom, SettledAtomsReport,
    SettledBond, SettledBondsReport,
};
pub(crate) use sites::{free_sites, perpendicular, unit, Rotation, SiteQuery};
pub(crate) use valences::builder_max_valence;
pub use valences::{
    builder_charges, builder_hydrogen_count, builder_valences, BuilderValenceError,
};

#[cfg(doc)]
use crate::structure::StructureEditor;
