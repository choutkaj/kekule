//! Geometry-bearing structures of coordinate-free topologies.
//!
//! - Dense per-atom arrays: [`Positions`] (canonical length unit, nm),
//!   [`Velocities`] (nm/ps) and [`Forces`] (kJ/(mol·nm)). They carry no atom
//!   identity and are validated against a topology by their owner.
//! - [`Structure`]: one shared, immutable [`Topology`](crate::topology::Topology)
//!   plus one complete structure of it: positions, an optional periodic cell,
//!   optional occupancies and B-factors, and structure-scoped properties.
//! - Publication: [`StructureBuilder`] assembles a new topology and its positions
//!   together; [`StructureEditor`] edits a detached draft of a structure whose
//!   positions follow atoms through structural edits, and
//!   [`StructureEditor::append_structure`] imports another structure
//!   ([`StructureAppend`], [`StructureAppendReport`]). [`StructureEditor::topology_mut`]
//!   ([`StructureEditorTopologyMut`]) and [`StructureEditor::properties_mut`]
//!   ([`StructureEditorPropertiesMut`]) reach annotations without geometry.
//! - [`StructureSet`]: many structures of one shared topology as [`Geometry`]s,
//!   each shared and copied on write, with an optional [`Order`] that places
//!   every geometry ([`TimeAxis`], steps, or coordinate values; [`Place`]) and
//!   an optional [`Distribution`] over all geometries or some of them
//!   ([`Members`]) grouped into [`State`]s by a [`Partition`]. A trajectory,
//!   a conformer ensemble, and a scan are structure sets. [`GeometryView`],
//!   [`GeometryMut`], and [`StateView`] read and edit them, and
//!   [`DistributionMut`] edits the distribution; the set's operations fail
//!   with a [`StructureSetError`]. Its storage ([`GeometryStorage`]) is a type
//!   parameter: [`InMemory`] by default, or a file's
//!   ([`crate::traj::io::IndexedFile`]), whose geometries are read on demand.
//!   Sets of either storage, and a host's own readers, read through
//!   [`StructureSetRead`], which fails with [`StructureSetReadError`];
//!   [`StructureSet::from_reader`] loads chosen geometries of any reader into
//!   memory. File codecs and streaming live in [`crate::traj`].
//! - Borrowed access: [`StructureView`] is the common input of coordinate-dependent
//!   kernels. Structures, geometries ([`GeometryView`]), and codec frames
//!   ([`crate::traj::FrameView`]) all borrow as one without copying
//!   ([`AsStructureView`]); [`StructureAtomView`] pairs one topology atom with its
//!   structure state. [`StructureView::kind`] says whether a view borrows a structure
//!   or a geometry ([`StructureKind`]); a geometry's [`SetContext`] reports
//!   its place, [`Membership`], velocities, and forces, so that state is
//!   never silently dropped.
//! - [`measure`]: checked distances, angles, dihedrals, bond dihedrals and
//!   distance-based selections.
//! - [`placement`]: builder geometry on a [`StructureEditor`] draft: hydrogens
//!   completed at ideal positions or collapsed back into counts, an atom's
//!   ideal shape and free directions, atoms attached to existing ones or
//!   placed as new molecules, and names for the atoms placement adds.
//! - Geometry edits: [`DistanceEdit`], [`AngleEdit`] and [`DihedralEdit`] (see
//!   below), and rigid motions with [`Structure::apply_rigid_transform`] and
//!   [`Structure::transform_atoms`].
//!
//! Owners accept semantic atom and bond IDs (and, in `*_at` methods, dense
//! indexes), translate IDs to dense indexes, and check every array's length.
//! Each owner has its own error type ([`StructureError`], [`StructureBuildError`],
//! [`StructureEditError`], [`StructureSetError`], [`PartitionError`],
//! [`GeometryEditError`]); [`GeometryError`] and [`DenseArrayError`] report
//! invalid structure state shared by all of them.
//!
//! ```
//! use kekule::{math::Point3, smiles, units::{Quantity, NANOMETER}};
//! use kekule::structure::{Structure, Positions};
//!
//! let positions = Positions::new(Quantity::new(
//!     [Point3::origin(), Point3::new(0.125, 0.0, 0.0)], NANOMETER,
//! ))?;
//! let structure = Structure::new(smiles::to_topology("CO")?, positions)?;
//! assert_eq!(structure.atom_count(), 2);
//! let oxygen = structure.topology().atom_ids()[1];
//! let position = structure.position(oxygen).unwrap().into_unit(NANOMETER)?;
//! assert_eq!(position.into_value(), Point3::new(0.125, 0.0, 0.0));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Geometry editing and scans
//!
//! [`Structure::set_distance`], [`Structure::set_angle`], and [`Structure::set_dihedral`]
//! move connected fragments. Prepare a [`DihedralEdit`] (or [`DistanceEdit`] /
//! [`AngleEdit`]) once for repeated absolute targets on the same topology:
//! ```
//! use kekule::{smiles, math::Point3, structure::{Structure, Positions, DihedralEdit},
//!     units::{Quantity, ANGSTROM, DEGREE}};
//! let topology = smiles::to_topology("CCCC")?;
//! let positions = Positions::new(Quantity::new([
//!     Point3::new(0.0, 1.0, 0.0), Point3::origin(),
//!     Point3::new(1.5, 0.0, 0.0), Point3::new(1.5, 1.0, 1.0),
//! ], ANGSTROM))?;
//! let mut structure = Structure::new(topology, positions)?;
//! let [a, b, c, d] = structure.topology().atom_ids().try_into().unwrap();
//! structure.set_distance(b, c, Quantity::new(1.6, ANGSTROM))?;
//! let edit = DihedralEdit::new(structure.shared_topology(), [a, b, c, d])?;
//! for degrees in [-180.0, -60.0, 60.0, 180.0] {
//!     let mut geometry = structure.clone();
//!     edit.apply(&mut geometry, Quantity::new(degrees, DEGREE))?;
//!     // Evaluate or collect geometry; structure remains the scan baseline.
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Explicit selections can intentionally deform rings or boundary bonds:
//! ```
//! use kekule::{structure::{AngleEdit, Structure}, topology::AtomSelection,
//!     units::{Quantity, DEGREE}};
//! # fn example(structure: &mut Structure) -> Result<(), Box<dyn std::error::Error>> {
//! let [a, b, c] = structure.topology().atom_ids()[..3].try_into().unwrap();
//! let topology = structure.shared_topology();
//! let moving = AtomSelection::from_atoms(topology.clone(), [c])?;
//! let edit = AngleEdit::with_moving_atoms(topology.clone(), [a, b, c], &moving)?;
//! edit.apply(structure, Quantity::new(120.0, DEGREE))?;
//! // Only C moves; bonds to other unselected neighbors can deform.
//! # Ok(())
//! # }
//! ```
//!
//! Geometry edits are atomic, use stored Cartesian coordinates, and preserve
//! topology, properties, and cell. Apply periodic preprocessing explicitly.
//! [`Structure::translate`], [`Structure::rotate`], and [`Structure::transform_atoms`] also
//! support rigid selection manipulation and keep the cell;
//! [`Structure::apply_rigid_transform`] moves the whole structure and rotates its cell.
//! No relaxation is performed.

mod dynamics;
mod geometry_edit;
pub mod measure;
mod owner;
mod partition;
pub mod placement;
mod positions;
pub(crate) mod scalars;
mod state;
mod structure_editor;
mod structure_set;
mod structure_state;

pub use dynamics::{Forces, Velocities};
pub use geometry_edit::{AngleEdit, DihedralEdit, DistanceEdit, GeometryEditError};
pub(crate) use owner::sealed::Sealed as StructureViewSealed;
pub use owner::{
    AsStructureView, Structure, StructureAtomView, StructureBuildError, StructureBuilder,
    StructureError, StructureKind, StructureTopologyMut, StructureView,
};
pub use partition::{Partition, PartitionError};
pub use positions::{DenseArrayError, Positions};
pub(crate) use state::subset_rows;
pub use state::GeometryError;
pub use structure_editor::{
    StructureAppend, StructureAppendReport, StructureEditError, StructureEditor,
    StructureEditorProperties, StructureEditorPropertiesMut, StructureEditorTopologyMut,
};
pub(crate) use structure_set::{
    check_bound, sealed as geometry_storage, visit_geometry, DecodedGeometry,
};
pub use structure_set::{
    Distribution, DistributionMut, Geometry, GeometryIter, GeometryMut, GeometryStorage,
    GeometryView, InMemory, Members, Membership, Order, Place, Requirement, SetContext, State,
    StateIter, StateView, StructureSet, StructureSetError, StructureSetRead, StructureSetReadError,
    TimeAxis,
};

#[cfg(test)]
mod structure_tests;
#[cfg(test)]
mod tests;
