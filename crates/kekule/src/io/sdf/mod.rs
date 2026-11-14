//! SDF: record parsing and interpretation, structure sets from records, and
//! writing.

mod document;
pub(crate) mod entry_points;
mod structure_set;
mod write;

pub(crate) use document::{parse_sdf_document, parse_sdf_reader};
pub use document::{
    SdfDataField, SdfDocument, SdfParseOptions, SdfRecord, SdfRecordInterpretation,
    SdfRecordInterpretationReport, SdfRecords,
};
pub(crate) use structure_set::sdf_structure_set;
pub(crate) use write::{write_sdf, write_sdf_to};
pub use write::{
    SdfRecordSource, SdfWriteOptions, SDF_COORDINATE_FIELD, SDF_COORDINATE_UNIT_FIELD,
    SDF_STEP_FIELD,
};
