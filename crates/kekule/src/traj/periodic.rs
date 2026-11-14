//! Periodic coordinate reconstruction, imaging, imaging onto a reference
//! frame, and temporal unwrapping.
//!
//! These operations change positions only, in place. The shared topology, cells,
//! velocities, forces, order, distribution, and all properties are retained.
//! Set functions publish only after every geometry succeeds; clone the set first
//! to keep the original. Every processed geometry must have a periodic cell.
//!
//! Molecules and bonds come from the authoritative topology. No bonds are guessed
//! from distances or hierarchy. Orthorhombic, triclinic, rotated, and partially
//! periodic cells are supported. Molecular reconstruction uses shortest Cartesian
//! bond images; temporal unwrapping uses continuity in fractional cell coordinates.
//! Molecular image ties are resolved deterministically. The exact nearest-image
//! search permits up to one million candidates per displacement and fails explicitly
//! for cells requiring more; it never substitutes an approximate image silently.
//!
//! Set functions operate on an in-memory [`StructureSet`], whose geometries they
//! call frames; [`MoleculeImager`] and [`TrajectoryUnwrapper`] apply the same
//! operations to streamed [`FrameBuffer`](crate::traj::FrameBuffer)s.
//!
//! # Examples
//!
//! An ethane C–C bond split across a 1 nm periodic boundary is made whole:
//!
//! ```
//! use std::sync::Arc;
//!
//! use kekule::math::{PeriodicCell, Point3, Vector3};
//! use kekule::structure::{Positions, Geometry, StructureSet};
//! use kekule::units::{Quantity, NANOMETER};
//! use kekule::{smiles, topology::Topology};
//! use kekule::traj::periodic;
//!
//! let molecule = smiles::to_molecules("CC")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//! let lengths = Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER);
//! let cell = PeriodicCell::orthorhombic(lengths, [true; 3])?;
//! let points = vec![Point3::new(0.125, 0.0, 0.0), Point3::new(0.875, 0.0, 0.0)];
//! let positions = Positions::from_vec(Quantity::new(points, NANOMETER))?;
//! let frame = Geometry::new(positions).with_cell(cell);
//! let mut trajectory = StructureSet::from_geometries(topology, [frame])?;
//!
//! periodic::make_molecules_whole(&mut trajectory)?;
//! let frame = trajectory.get(0).unwrap();
//! let moved = frame.positions().get(1).unwrap().value_in(NANOMETER)?;
//! assert_eq!(moved, Point3::new(-0.125, 0.0, 0.0));
//! let anchor = frame.positions().get(0).unwrap().value_in(NANOMETER)?;
//! assert_eq!(anchor, Point3::new(0.125, 0.0, 0.0));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Reading a file, then reconstructing, unwrapping and imaging:
//!
//! ```no_run
//! use kekule::{mmcif, topology::AtomSelection};
//! use kekule::traj::{io::read_structure_set, periodic};
//!
//! let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?, Default::default())?;
//! let topology = document.interpret(Default::default())?.into_topology();
//! let mut trajectory = read_structure_set("trajectory.xtc", topology.clone(), Default::default())?;
//! periodic::make_molecules_whole(&mut trajectory)?;
//!
//! // A continuous path through time, retaining the initial molecular images.
//! let mut continuous = trajectory.clone();
//! periodic::unwrap_trajectory(&mut continuous)?;
//!
//! // A separate per-frame view with molecules centered around selected anchors.
//! let anchors = AtomSelection::all(topology.clone()); // Or select the solute's atoms.
//! periodic::image_molecules(&mut trajectory, &anchors)?;
//! trajectory.superpose(0, &anchors, Default::default())?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! For a large trajectory, retain the bond plan and temporal state across reads:
//!
//! ```no_run
//! use kekule::mmcif;
//! use kekule::traj::{
//!     io::open_trajectory, periodic::{MoleculeImager, TrajectoryUnwrapper},
//!     TrajectoryReader,
//! };
//! let topology = mmcif::parse_str(&std::fs::read_to_string("system.cif")?, Default::default())?
//!     .interpret(Default::default())?.into_topology();
//! let mut reader = open_trajectory("trajectory.xtc", topology.clone(), Default::default())?;
//! let mut frame = reader.frame_buffer();
//! let imager = MoleculeImager::new(topology.clone());
//! let mut unwrapper = TrajectoryUnwrapper::new(topology);
//! let mut index = 0;
//! while reader.read_next(&mut frame)? {
//!     imager.make_whole(index, &mut frame)?;
//!     unwrapper.unwrap(index, &mut frame)?;
//!     // Analyze or write this frame here. Downsample only after unwrapping.
//!     index += 1;
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod imaging;
mod stream;

pub use imaging::{
    image_molecules, image_onto, make_molecules_whole, unwrap_trajectory, PeriodicError,
};
pub use stream::{MoleculeImager, TrajectoryUnwrapper};

#[cfg(doc)]
use crate::structure::StructureSet;
