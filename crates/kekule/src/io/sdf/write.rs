//! SDF writing: independent Molfile records, each with a title and data fields.
use std::io::Write;

use crate::structure::{
    AsStructureView, GeometryView, Place, Structure, StructureKind, StructureView,
};
use crate::traj::FrameView;
use crate::units::named_unit;

use super::document::{SdfDataField, SdfRecordInterpretation};
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail, GeometryField};
use crate::io::molfile::write::{check_cell, render_structure};
use crate::io::write_policy::{check_geometry_state, Omissions, StatePolicies};
use crate::io::{MolfileWriteVersion, UnrepresentablePolicy, WriteOmission, WriteReport};

/// Options of SDF writing; the default writes
/// [`MolfileWriteVersion::Auto`] records and fails on content an SDF record
/// cannot carry.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfWriteOptions {
    version: MolfileWriteVersion,
    cell: UnrepresentablePolicy,
    geometry_metadata: UnrepresentablePolicy,
    properties: UnrepresentablePolicy,
}

impl SdfWriteOptions {
    /// Connection-table version of each record's CTAB; `Auto` chooses per
    /// record.
    pub const fn version(&self) -> MolfileWriteVersion {
        self.version
    }

    /// Sets [`Self::version`].
    #[must_use]
    pub const fn with_version(mut self, version: MolfileWriteVersion) -> Self {
        self.version = version;
        self
    }

    /// A structure's periodic cell, which an SDF record cannot carry: by
    /// default a record whose structure has a cell fails to write, and under
    /// [`UnrepresentablePolicy::Omit`] the report counts the records written
    /// without one ([`WriteOmission::Cell`]).
    pub const fn cell(&self) -> UnrepresentablePolicy {
        self.cell
    }

    /// Sets [`Self::cell`].
    #[must_use]
    pub const fn with_cell(mut self, cell: UnrepresentablePolicy) -> Self {
        self.cell = cell;
        self
    }

    /// The state that an SDF record cannot carry of a geometry, whether of a
    /// structure set or a trajectory frame: its time and MD step in a time
    /// order, its velocities and forces, and a coordinate unit that is not a
    /// unit constant; and the occupancies and
    /// B factors of any structure. By default writing a structure with any of
    /// them fails. Under [`UnrepresentablePolicy::Omit`] the report counts
    /// the records written without each ([`WriteOmission::GeometryState`],
    /// [`WriteOmission::Occupancies`], [`WriteOmission::BFactors`]). A
    /// geometry's step or coordinate is written as data fields.
    pub const fn geometry_metadata(&self) -> UnrepresentablePolicy {
        self.geometry_metadata
    }

    /// Sets [`Self::geometry_metadata`].
    #[must_use]
    pub const fn with_geometry_metadata(
        mut self,
        geometry_metadata: UnrepresentablePolicy,
    ) -> Self {
        self.geometry_metadata = geometry_metadata;
        self
    }

    /// A structure's per-atom, per-bond, and own properties and the
    /// properties of a geometry's set, which an SDF record cannot carry: by
    /// default a record of a structure with any fails to write, and under
    /// [`UnrepresentablePolicy::Omit`] the report counts the records written
    /// without each ([`WriteOmission::GeometryState`]). A caller writes the
    /// values a record should carry as data fields.
    pub const fn properties(&self) -> UnrepresentablePolicy {
        self.properties
    }

    /// Sets [`Self::properties`].
    #[must_use]
    pub const fn with_properties(mut self, properties: UnrepresentablePolicy) -> Self {
        self.properties = properties;
        self
    }
}

/// The data field a geometry's step in a step order is written to.
pub const SDF_STEP_FIELD: &str = "kekule.step";
/// The data field a geometry's coordinate in a coordinate order is written
/// to, in the order's unit, which [`SDF_COORDINATE_UNIT_FIELD`] names.
pub const SDF_COORDINATE_FIELD: &str = "kekule.coordinate";
/// The data field naming the unit of [`SDF_COORDINATE_FIELD`], as the
/// unit displays itself: its symbol, such as `A` or `deg`, `1` for a
/// dimensionless coordinate, or its scale and dimension.
pub const SDF_COORDINATE_UNIT_FIELD: &str = "kekule.coordinate.unit";

/// The data fields a geometry's step or coordinate place uses.
const PLACE_FIELDS: [&str; 3] = [
    SDF_STEP_FIELD,
    SDF_COORDINATE_FIELD,
    SDF_COORDINATE_UNIT_FIELD,
];

/// One SDF record: a structure view with an optional title and data fields.
///
/// Convert from `&Structure`, a borrowed geometry of a structure set, a trajectory
/// frame, or `&SdfRecordInterpretation` (which keeps its title and data
/// fields); `sdf::write(&set, options)` writes a record per geometry. A
/// geometry's step place follows its data fields as [`SDF_STEP_FIELD`], and
/// a coordinate place as [`SDF_COORDINATE_FIELD`] with
/// [`SDF_COORDINATE_UNIT_FIELD`]. A coordinate
/// unit must be one of the unit constants of [`crate::units`], which the
/// reader knows by symbol; another fails to write unless
/// [`SdfWriteOptions::geometry_metadata`] omits the coordinate. A geometry's
/// time, MD step, velocities, and forces, and a structure's
/// occupancies and B factors, fail to write unless
/// [`SdfWriteOptions::geometry_metadata`] omits them; structure and set
/// properties fail unless [`SdfWriteOptions::properties`] omits them. Add
/// the values a record should carry, such as a geometry's energy, as data
/// fields.
#[derive(Debug, Clone, Copy)]
pub struct SdfRecordSource<'a> {
    structure: StructureView<'a>,
    title: &'a str,
    data_fields: &'a [SdfDataField],
}

impl<'a> SdfRecordSource<'a> {
    /// A record of any structure view, with an empty title and no data fields.
    pub fn structure(structure: &'a (impl AsStructureView + ?Sized)) -> Self {
        Self::view(structure.as_structure_view())
    }

    const fn view(structure: StructureView<'a>) -> Self {
        Self {
            structure,
            title: "",
            data_fields: &[],
        }
    }

    /// Sets the Molfile title line of the record. Writing fails with
    /// `Unsupported` for a title containing a line break or equal to `$$$$`.
    #[must_use]
    pub fn with_title(mut self, title: &'a str) -> Self {
        self.title = title;
        self
    }

    /// Sets the data fields, written in order; a geometry's place fields
    /// still follow them. Writing fails with `Unsupported`
    /// for a field SDF cannot frame (an empty, untrimmed, or bracketed name,
    /// or a value with a carriage return, blank line, or `$$$$` line) and
    /// with `Inconsistent` for a field named like the geometry metadata.
    #[must_use]
    pub fn with_data_fields(mut self, data_fields: &'a [SdfDataField]) -> Self {
        self.data_fields = data_fields;
        self
    }
}

impl<'a> From<&'a Structure> for SdfRecordSource<'a> {
    fn from(structure: &'a Structure) -> Self {
        Self::structure(structure)
    }
}

impl<'a> From<FrameView<'a>> for SdfRecordSource<'a> {
    fn from(frame: FrameView<'a>) -> Self {
        Self::view(frame.as_structure_view())
    }
}

impl<'a, R: ?Sized> From<GeometryView<'a, R>> for SdfRecordSource<'a> {
    fn from(geometry: GeometryView<'a, R>) -> Self {
        Self::view(geometry.as_structure_view())
    }
}

impl<'a> From<&'a SdfRecordInterpretation> for SdfRecordSource<'a> {
    fn from(record: &'a SdfRecordInterpretation) -> Self {
        Self::structure(record.structure())
            .with_title(record.title())
            .with_data_fields(record.data_fields())
    }
}

/// Writes records to a string. See [`write_sdf_to`].
pub(crate) fn write_sdf<'a, R: Into<SdfRecordSource<'a>>>(
    records: impl IntoIterator<Item = R>,
    options: SdfWriteOptions,
) -> Result<(String, WriteReport), FormatError> {
    let mut output = Vec::new();
    let report = write_sdf_to(&mut output, records, options)?;
    let text = String::from_utf8(output).expect("SDF writer emits UTF-8");
    Ok((text, report))
}

/// Writes independent records in input order and reports what they leave
/// out. A structure set writes one record per geometry (`sdf::write(&set,
/// options)`).
pub(crate) fn write_sdf_to<'a, R: Into<SdfRecordSource<'a>>>(
    writer: &mut impl Write,
    records: impl IntoIterator<Item = R>,
    options: SdfWriteOptions,
) -> Result<WriteReport, FormatError> {
    let mut omissions = Omissions::default();
    for (number, record) in (1..).zip(records) {
        write_record(writer, record.into(), options, &mut omissions)
            .map_err(|error| error.in_record(number))?;
    }
    Ok(omissions.into_report())
}

fn write_record(
    writer: &mut impl Write,
    record: SdfRecordSource<'_>,
    options: SdfWriteOptions,
    omissions: &mut Omissions,
) -> Result<(), FormatError> {
    validate_title(record.title)?;
    let policies = StatePolicies {
        writes_place: true,
        ..StatePolicies::positions_only(options.geometry_metadata, options.properties)
    };
    check_geometry_state(Format::Sdf, record.structure, policies)?;
    // A geometry's step or coordinate place, which travels as data fields.
    let place = match record.structure.kind() {
        StructureKind::InSet(context) => context
            .place()
            .filter(|place| matches!(place, Place::Step(_) | Place::Coordinate(_))),
        StructureKind::Standalone => None,
    };
    // Whether a geometry's own data uses the name.
    let reserved = |name: &str| place.is_some() && PLACE_FIELDS.contains(&name);
    // A coordinate unit travels as the symbol of a unit constant, which the
    // reader turns back into that unit.
    let written_place = match place {
        Some(Place::Coordinate(coordinate))
            if named_unit(&coordinate.unit().to_string()) != Some(coordinate.unit()) =>
        {
            if options.geometry_metadata == UnrepresentablePolicy::Fail {
                return Err(FormatError::write(
                    Format::Sdf,
                    DiagnosticKind::Unsupported,
                    format!(
                        "SDF cannot represent the coordinate unit {}.",
                        coordinate.unit()
                    ),
                )
                .with_detail(FormatErrorDetail::UnrepresentableState(
                    GeometryField::OrderCoordinate,
                )));
            }
            omissions.record(WriteOmission::GeometryState {
                field: GeometryField::OrderCoordinate,
                geometries: 1,
            });
            None
        }
        place => place,
    };
    let own_fields = written_place.map_or_else(Vec::new, place_fields);
    for field in record.data_fields {
        validate_data_field(field)?;
        // A placed geometry owns these names, so a caller's field never
        // adds a second place.
        if reserved(field.name()) {
            return Err(FormatError::write(
                Format::Sdf,
                DiagnosticKind::Inconsistent,
                format!(
                    "The data field {:?} is reserved for the data of the geometry.",
                    field.name()
                ),
            )
            .with_detail(FormatErrorDetail::DataField {
                name: field.name().to_owned(),
            }));
        }
    }
    let molfile = |error: FormatError| error.in_format(Format::Sdf);
    check_cell(record.structure, options.cell).map_err(molfile)?;
    let ctab =
        render_structure(record.structure, record.title, options.version).map_err(molfile)?;
    if record.structure.cell().is_some() {
        omissions.record(WriteOmission::Cell { structures: 1 });
    }
    omissions.record_geometry_state(record.structure, policies);
    let sink = |error| FormatError::sink(Format::Sdf, error);
    writer.write_all(ctab.as_bytes()).map_err(sink)?;
    for field in record.data_fields.iter().chain(&own_fields) {
        writeln!(writer, ">  <{}>\n{}\n", field.name(), field.value()).map_err(sink)?;
    }
    writer.write_all(b"$$$$\n").map_err(sink)
}

/// A geometry's step, or its coordinate and the coordinate's unit, as data
/// fields, the number in the shortest decimal text that reads back to the
/// same value.
fn place_fields(place: Place) -> Vec<SdfDataField> {
    match place {
        Place::Step(step) => vec![SdfDataField::new(SDF_STEP_FIELD, step.to_string())],
        Place::Coordinate(coordinate) => vec![
            SdfDataField::new(SDF_COORDINATE_FIELD, coordinate.value().to_string()),
            SdfDataField::new(SDF_COORDINATE_UNIT_FIELD, coordinate.unit().to_string()),
        ],
        _ => Vec::new(),
    }
}

fn validate_title(title: &str) -> Result<(), FormatError> {
    let invalid = |message: &str| {
        Err(
            FormatError::write(Format::Sdf, DiagnosticKind::Unsupported, message)
                .with_detail(FormatErrorDetail::Title),
        )
    };
    if title.contains(['\r', '\n']) {
        return invalid("The title contains a line break.");
    }
    if title.trim() == "$$$$" {
        return invalid("The title is \"$$$$\".");
    }
    Ok(())
}

fn validate_data_field(field: &SdfDataField) -> Result<(), FormatError> {
    let invalid = |message: String| {
        Err(
            FormatError::write(Format::Sdf, DiagnosticKind::Unsupported, message).with_detail(
                FormatErrorDetail::DataField {
                    name: field.name().to_owned(),
                },
            ),
        )
    };
    let name = field.name();
    if name.is_empty() || name.trim() != name || name.contains(['<', '>', '\r', '\n']) {
        return invalid(format!(
            "The data field name {name:?} is empty, padded, or holds angle brackets or line breaks."
        ));
    }
    let value = field.value();
    if value.contains('\r') {
        return invalid(format!(
            "The value of data field {name:?} contains a carriage return."
        ));
    }
    if !value.is_empty() && value.split('\n').any(str::is_empty) {
        return invalid(format!(
            "The value of data field {name:?} contains a blank line."
        ));
    }
    if value.lines().any(|line| line.trim() == "$$$$") {
        return invalid(format!(
            "The value of data field {name:?} contains a \"$$$$\" line."
        ));
    }
    Ok(())
}
