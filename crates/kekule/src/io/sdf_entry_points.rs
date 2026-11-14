//! The public SDF entry points that `crate::sdf` re-exports.

use crate::diagnostics::FormatError;
use crate::io::{
    SdfDocument, SdfParseOptions, SdfRecordInterpretation, SdfRecordSource, SdfRecords,
    SdfWriteOptions,
};

/// Assembles interpreted records, in order, into one structure set whose
/// geometries share the first record's topology, record `i` becoming geometry
/// `i`. Sibling records never form a set implicitly, and the set states
/// only what the records state.
///
/// A record's data fields stay with its interpretation; none says how
/// much a geometry weighs or which state it belongs to, which an
/// [`Ensemble`](crate::structure::Ensemble) of the set states.
///
/// When the records carry [`STEP_FIELD`](crate::sdf::STEP_FIELD), or
/// [`COORDINATE_FIELD`](crate::sdf::COORDINATE_FIELD) with
/// [`COORDINATE_UNIT_FIELD`](crate::sdf::COORDINATE_UNIT_FIELD), as
/// [`write()`] emits for geometries of a step
/// or coordinate order, the geometries are ordered by them. The unit is the
/// symbol of a unit constant of [`crate::units`].
///
/// ```
/// use kekule::math::Point3;
/// use kekule::sdf;
/// use kekule::structure::{Structure, Order, Positions};
/// use kekule::units::{Quantity, ANGSTROM};
///
/// let water = kekule::smiles::to_molecules("O")?.remove(0);
/// let positions = Positions::new(Quantity::new(vec![Point3::new(0.0, 0.0, 0.0)], ANGSTROM))?;
/// let structure = Structure::from_molecule(water, &positions)?;
/// let (text, _) = sdf::write([&structure, &structure], Default::default())?;
/// let records = sdf::parse_str(&text, Default::default())?
///     .records()
///     .iter()
///     .map(|record| record.interpret(Default::default()))
///     .collect::<Result<Vec<_>, _>>()?;
/// let set = sdf::structure_set(&records)?;
/// assert_eq!(set.len(), 2);
/// assert_eq!(set.order(), None);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
///
/// Returns a [`FormatError`] of kind
/// [`DiagnosticKind::Inconsistent`](crate::diagnostics::DiagnosticKind::Inconsistent)
/// when the records do not form a structure set (none, or differing
/// topologies), a field is repeated, some records carry a step or
/// coordinate and others do not, steps and coordinates are mixed, or
/// coordinate units differ, and
/// [`DiagnosticKind::InvalidValue`](crate::diagnostics::DiagnosticKind::InvalidValue)
/// when a field is not a valid step, coordinate, or unit.
pub fn structure_set<'a>(
    records: impl IntoIterator<Item = &'a SdfRecordInterpretation>,
) -> Result<crate::structure::StructureSet, FormatError> {
    crate::io::sdf_structure_set(records)
}

/// Streams the records of an SDF input from `reader`, each bounded by
/// `options`.
///
/// Each item is one record or the failure of one record; reading resumes at
/// the next `$$$$` delimiter after a malformed, oversized, or non-UTF-8
/// record. Memory is bounded by one record.
pub fn records<R: std::io::Read>(reader: R, options: SdfParseOptions) -> SdfRecords<R> {
    SdfRecords::new(reader, options)
}

/// Parses every record of an SDF input, failing on the first malformed
/// record. Use [`records`] to skip past malformed records instead.
///
/// # Errors
///
/// Returns the [`FormatError`] of the first record that fails: kind
/// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
/// when the input or a record exceeds `options`, and the errors of
/// [`molfile::parse_str`](crate::molfile::parse_str) for a malformed
/// connection table or data field.
pub fn parse_str(input: &str, options: SdfParseOptions) -> Result<SdfDocument, FormatError> {
    crate::io::parse_sdf_document(input, options)
}

/// Reads and parses every record of an SDF input from `reader`.
///
/// # Errors
///
/// As [`parse_str`]. The configured byte limit is enforced while reading;
/// a read failure or non-UTF-8 input fails with
/// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn parse_reader(
    reader: impl std::io::Read,
    options: SdfParseOptions,
) -> Result<SdfDocument, FormatError> {
    crate::io::parse_sdf_reader(reader, options)
}

/// Writes independent records in input order, with the
/// [`WriteReport`](crate::diagnostics::WriteReport) of what they left
/// out under [`UnrepresentablePolicy::Omit`](crate::diagnostics::UnrepresentablePolicy::Omit):
/// structures' cells, geometry state, occupancies, B factors, and
/// properties ([`SdfWriteOptions::cell`],
/// [`SdfWriteOptions::geometry_metadata`], [`SdfWriteOptions::properties`]),
/// counted over all records.
///
/// Records convert from structures, geometries of a structure set, and
/// trajectory frames (`sdf::write(&set, options)` writes one record per
/// geometry, which [`structure_set`] reads back with its step or
/// coordinate), or
/// interpreted records, which keep their titles and data fields.
///
/// # Errors
///
/// The errors of [`molfile::write`](crate::molfile::write) for any record,
/// and
/// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
/// or
/// [`DiagnosticKind::Inconsistent`](crate::diagnostics::DiagnosticKind::Inconsistent)
/// for a title or data field the format cannot hold.
pub fn write<'a, R: Into<SdfRecordSource<'a>>>(
    records: impl IntoIterator<Item = R>,
    options: SdfWriteOptions,
) -> Result<(String, crate::diagnostics::WriteReport), FormatError> {
    crate::io::write_sdf(records, options)
}

/// Streams the records [`write()`] would return, each as soon as it is
/// rendered, and returns their report; wrap an unbuffered sink such as
/// a `File` in a [`BufWriter`](std::io::BufWriter).
///
/// # Errors
///
/// The errors of [`write()`]; records before the failing one have already
/// reached `writer`. A failing `writer` gives
/// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn write_to<'a, R: Into<SdfRecordSource<'a>>>(
    writer: &mut impl std::io::Write,
    records: impl IntoIterator<Item = R>,
    options: SdfWriteOptions,
) -> Result<crate::diagnostics::WriteReport, FormatError> {
    crate::io::write_sdf_to(writer, records, options)
}
