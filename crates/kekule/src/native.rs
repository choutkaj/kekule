//! Kekule's native binary serialization: exact, versioned persistence.
//!
//! Scientific formats exchange chemistry with other tools; they are not exact
//! persistence. The native format stores one [`Molecule`](crate::molecule::Molecule),
//! [`Topology`](crate::topology::Topology), [`Structure`](crate::structure::Structure),
//! [`StructureSet`](crate::structure::StructureSet), or
//! [`Ensemble`](crate::structure::Ensemble) losslessly: represented
//! chemistry, definitions and instances, semantic IDs, dense order,
//! hierarchy, classifications with their explicit or inferred assignment,
//! properties at every scope with units, coordinates as exact bits, and every
//! structure, geometry, order, and ensemble distribution field. A structure or structure
//! set may carry its structure interpretation
//! reports: write a
//! [`StructureInterpretation`](crate::macromolecular::StructureInterpretation)
//! or [`StructureSetInterpretation`](crate::macromolecular::StructureSetInterpretation),
//! or pair a structure or set with reports held apart through
//! [`NativeSource::interpreted_structure`]
//! and [`NativeSource::interpreted_structure_set`],
//! and read them back with
//! [`NativeObject::into_interpretation`]
//! or [`NativeObject::into_structure_set_interpretation`].
//! A [`StructureSetDescription`] stores a structure
//! set whose geometries are kept elsewhere: its topology, properties, order,
//! and the number of its geometries; an
//! [`EnsembleDescription`] adds a distribution.
//! Installed perception is stored with its provenance unless
//! [`NativeWriteOptions::with_perception`]
//! omits it.
//!
//! [`native::write_to`](crate::native::write_to) streams to a [`std::io::Write`] and [`native::read`](crate::native::read) reads exactly one
//! object from a [`std::io::Read`]; neither buffers the caller's I/O, so wrap
//! unbuffered files in [`std::io::BufWriter`] or [`std::io::BufReader`].
//! Decoding validates untrusted input through the owners' reconstruction rules
//! under the resource limits of [`NativeReadOptions`] and fails with a
//! [`FormatError`](crate::diagnostics::FormatError), never a panic. A decoded
//! topology keeps every semantic ID but is a new publication with a fresh
//! layout identity. The byte layout and versioning policy are specified in
//! `docs/native-format.md`.
//!
//! ```
//! use kekule::native::{self, NativeReadOptions, NativeWriteOptions};
//! use kekule::smiles;
//!
//! let mut molecule = smiles::to_molecules("c1ccccc1O")?.remove(0);
//! molecule.perceive()?;
//! let bytes = native::write(&molecule, NativeWriteOptions::default())?;
//! let decoded = native::read(bytes.as_slice(), NativeReadOptions::default())?.into_molecule()?;
//! assert_eq!(decoded, molecule);
//! assert_eq!(decoded.perception(), molecule.perception());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::io::native::entry_points::{read, write, write_to};
pub use crate::io::{
    EnsembleDescription, NativeKind, NativeObject, NativeReadOptions, NativeSource,
    NativeWriteOptions, StructureSetDescription, DEFAULT_NATIVE_MAX_ATOMS,
    DEFAULT_NATIVE_MAX_BONDS, DEFAULT_NATIVE_MAX_INPUT_BYTES, DEFAULT_NATIVE_MAX_ITEMS,
    DEFAULT_NATIVE_MAX_STRING_BYTES, NATIVE_FORMAT_VERSION, NATIVE_MAGIC,
};
