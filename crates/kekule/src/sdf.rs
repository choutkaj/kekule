//! Record-oriented SDF parsing, interpretation, and writing.
//!
//! An SDF input is a sequence of independent records. [`sdf::records`](crate::sdf::records) streams
//! them one at a time, so a malformed record fails alone; [`sdf::parse_str`](crate::sdf::parse_str)
//! collects every record into a [`sdf::SdfDocument`](crate::sdf::SdfDocument). Interpret each record
//! explicitly with [`sdf::SdfRecord::interpret`](crate::sdf::SdfRecord::interpret); sibling records are not
//! merged into one structure or reinterpreted as a structure set.
//!
//! ```
//! use kekule::math::Point3;
//! use kekule::structure::{Structure, Positions};
//! use kekule::units::{Quantity, ANGSTROM};
//! use kekule::{sdf, smiles};
//!
//! let water = smiles::to_molecules("O")?.remove(0);
//! let positions = Positions::new(Quantity::new(vec![Point3::new(0.0, 0.0, 0.0)], ANGSTROM))?;
//! let structure = Structure::from_molecule(water, &positions)?;
//!
//! let (text, report) = sdf::write([&structure, &structure], Default::default())?;
//! assert!(report.is_empty());
//! let document = sdf::parse_str(&text, Default::default())?;
//! assert_eq!(document.records().len(), 2);
//! let first = document.records()[0].interpret(Default::default())?;
//! assert_eq!(first.structure().atom_count(), 1);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::io::sdf::entry_points::{
    parse_reader, parse_str, records, structure_set, write, write_to,
};
pub use crate::io::{
    SdfDataField, SdfDocument, SdfParseOptions, SdfRecord, SdfRecordInterpretation,
    SdfRecordInterpretationReport, SdfRecordSource, SdfRecords, SdfWriteOptions,
    SDF_COORDINATE_FIELD as COORDINATE_FIELD, SDF_COORDINATE_UNIT_FIELD as COORDINATE_UNIT_FIELD,
    SDF_STEP_FIELD as STEP_FIELD,
};
