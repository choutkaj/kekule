//! Checked measurements of stored Cartesian coordinates through [`StructureView`].
//!
//! Cells are ignored: these operations never wrap, image, unwrap, or align
//! coordinates. Apply the desired preprocessing before measuring. Angles use
//! radians and distances use the canonical length unit; results support ordinary
//! [`Quantity`] conversions. No structural ownership is materialized.
//!
//! [`distance`], [`angle`], and [`dihedral`] take a [`ConnectivityCheck`]:
//! [`ConnectivityCheck::Unrestricted`] accepts arbitrary atoms, including atoms
//! in different molecule instances, and
//! [`ConnectivityCheck::ConsecutiveBonds`] requires direct bonds between
//! consecutive atoms. Connectivity comes from the topology,
//! never from coordinate proximity; geometric validity is checked separately.
//! [`bond_dihedral`] selects a deterministic quartet for one bond; prepare a
//! [`BondDihedral`] once to reuse its topology-bound references across frames.
//!
//! ```
//! use kekule::{math::Point3, smiles, units::{Quantity, NANOMETER}};
//! use kekule::structure::{measure::{self, ConnectivityCheck}, Structure, Positions};
//!
//! let positions = Positions::new(Quantity::new(
//!     [Point3::origin(), Point3::new(0.125, 0.0, 0.0)], NANOMETER,
//! ))?;
//! let structure = Structure::new(smiles::to_topology("CO")?, positions)?;
//! let [c, o] = structure.topology().atom_ids().try_into().unwrap();
//! let length = measure::distance(
//!     structure.as_structure_view(), c, o, ConnectivityCheck::ConsecutiveBonds,
//! )?;
//! assert_eq!(length.into_unit(NANOMETER)?.into_value(), 0.125);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! A whole-residue pocket around one source-identified ligand:
//! ```
//! use kekule::{
//!     molecule::Element,
//!     structure::{measure, Structure},
//!     topology::AtomSelection,
//!     units::{Quantity, ANGSTROM},
//! };
//! # fn pocket(structure: &Structure) -> Result<Structure, Box<dyn std::error::Error>> {
//! let topology = structure.shared_topology();
//! let ligand = topology.residue_by_author("A", "215", None)?;
//! let reference = AtomSelection::for_residues(topology.clone(), [ligand.id()])?;
//! let hydrogen = AtomSelection::for_elements(
//!     topology.clone(), [Element::from_symbol("H").unwrap()],
//! )?;
//! let candidates = AtomSelection::all(topology.clone()).difference(&hydrogen)?;
//! let atoms = measure::within(
//!     structure.as_structure_view(), &candidates, &reference, Quantity::new(5.0, ANGSTROM),
//! )?.expand_to_residues();
//! let pocket = structure.subset(&atoms)?;
//! # Ok(pocket)
//! # }
//! ```

mod bond_dihedrals;
mod measurements;

pub use bond_dihedrals::{bond_dihedral, bond_dihedrals, BondDihedral, BondDihedralError};
pub use measurements::{angle, dihedral, distance, within, ConnectivityCheck, MeasurementError};
pub(super) use measurements::{angle_points, dihedral_points, norm, normalized};

#[cfg(doc)]
use crate::structure::StructureView;
#[cfg(doc)]
use crate::units::Quantity;
