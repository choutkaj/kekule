//! SDF writing: independent Molfile records, each with a title and data fields.
use std::io::Write;

use crate::structure::{
    AsModelView, Membership, Model, ModelView, Place, SampleView, StructureKind,
};
use crate::traj::FrameView;
use crate::units::{named_unit, KELVIN, KILOJOULE_PER_MOLE};

use super::molfile_write::{check_cell, render_model};
use super::sdf_document::{SdfDataField, SdfRecordInterpretation};
use super::write_error::{check_structure_state, StatePolicies};
use super::{MolfileWriteVersion, UnrepresentablePolicy};
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail, StructureField};

/// Options of SDF writing; the default writes
/// [`MolfileWriteVersion::Auto`] records and fails on a periodic cell.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfWriteOptions {
    version: MolfileWriteVersion,
    cell: UnrepresentablePolicy,
    sample_metadata: UnrepresentablePolicy,
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

    /// A model's periodic cell, which an SDF record cannot carry: by default
    /// a record whose model has a cell fails to write.
    pub const fn cell(&self) -> UnrepresentablePolicy {
        self.cell
    }

    /// Sets [`Self::cell`].
    #[must_use]
    pub const fn with_cell(mut self, cell: UnrepresentablePolicy) -> Self {
        self.cell = cell;
        self
    }

    /// The state of a sample that an SDF record cannot carry, whether of a
    /// sample set or a trajectory frame: its time and MD step in a time
    /// order, its velocities and forces, a coordinate unit that is not a unit
    /// constant, and its being outside the set's distribution, which a
    /// record cannot tell from membership. By default writing a sample with
    /// any of them fails. Its weight, state, state energy, and temperature,
    /// and its step or coordinate, are otherwise written as data fields.
    pub const fn sample_metadata(&self) -> UnrepresentablePolicy {
        self.sample_metadata
    }

    /// Sets [`Self::sample_metadata`].
    #[must_use]
    pub const fn with_sample_metadata(mut self, sample_metadata: UnrepresentablePolicy) -> Self {
        self.sample_metadata = sample_metadata;
        self
    }
}

/// The data field a sample's relative weight is written to and
/// [`sdf::sample_set`](crate::sdf::sample_set) reads it from.
pub const SDF_WEIGHT_FIELD: &str = "kekule.weight";
/// The data field the index of a sample's state is written to.
pub const SDF_STATE_FIELD: &str = "kekule.state";
/// The data field, `true`, that marks the record of a state's representative
/// sample.
pub const SDF_REPRESENTATIVE_FIELD: &str = "kekule.representative";
/// The data field a state's energy is written to, in kJ/mol, on the record of
/// its representative sample.
pub const SDF_ENERGY_FIELD: &str = "kekule.energy.kJ_per_mol";
/// The data field the temperature a distribution's populations describe is
/// written to, in kelvin, on the record of every member.
pub const SDF_TEMPERATURE_FIELD: &str = "kekule.temperature.K";

/// The data field a sample's step in a step order is written to.
pub const SDF_STEP_FIELD: &str = "kekule.step";
/// The data field a sample's coordinate in a coordinate order is written
/// to, in the order's unit, which [`SDF_COORDINATE_UNIT_FIELD`] names.
pub const SDF_COORDINATE_FIELD: &str = "kekule.coordinate";
/// The data field naming the unit of [`SDF_COORDINATE_FIELD`], as the
/// unit displays itself: its symbol, such as `A` or `deg`, `1` for a
/// dimensionless coordinate, or its scale and dimension.
pub const SDF_COORDINATE_UNIT_FIELD: &str = "kekule.coordinate.unit";

/// The data fields a sample's membership uses.
const SAMPLE_FIELDS: [&str; 5] = [
    SDF_WEIGHT_FIELD,
    SDF_STATE_FIELD,
    SDF_REPRESENTATIVE_FIELD,
    SDF_ENERGY_FIELD,
    SDF_TEMPERATURE_FIELD,
];

/// The data fields a sample's step or coordinate place uses.
const PLACE_FIELDS: [&str; 3] = [
    SDF_STEP_FIELD,
    SDF_COORDINATE_FIELD,
    SDF_COORDINATE_UNIT_FIELD,
];

/// One SDF record: a model view with an optional title and data fields.
///
/// Convert from `&Model`, a borrowed sample of a sample set, a trajectory
/// frame, or `&SdfRecordInterpretation` (which keeps its title and data
/// fields); `sdf::write(&set, options)` writes a record per sample. A
/// sample's membership in a distribution follows its data fields: its
/// weight, its state, and the distribution's temperature as
/// [`SDF_WEIGHT_FIELD`], [`SDF_STATE_FIELD`], and [`SDF_TEMPERATURE_FIELD`];
/// a representative's record is marked with [`SDF_REPRESENTATIVE_FIELD`]
/// and carries its state's energy as [`SDF_ENERGY_FIELD`], whichever way
/// its view is passed. A step place follows as [`SDF_STEP_FIELD`], and a
/// coordinate place as [`SDF_COORDINATE_FIELD`] with
/// [`SDF_COORDINATE_UNIT_FIELD`]; a sample with both a membership and a
/// place writes both groups. A coordinate unit must be one of the unit
/// constants of [`crate::units`], which the reader knows by symbol; another
/// fails to write unless [`SdfWriteOptions::sample_metadata`] omits the
/// coordinate. State, set, and structure properties are not written; add
/// the values a record should carry, such as a sample's energy, as data
/// fields. A sample's time, MD step, velocities, and forces fail to write
/// unless [`SdfWriteOptions::sample_metadata`] omits them.
#[derive(Debug, Clone, Copy)]
pub struct SdfRecordSource<'a> {
    model: ModelView<'a>,
    title: &'a str,
    data_fields: &'a [SdfDataField],
}

impl<'a> SdfRecordSource<'a> {
    /// A record of any model view, with an empty title and no data fields.
    pub fn model(model: &'a (impl AsModelView + ?Sized)) -> Self {
        Self::view(model.as_model_view())
    }

    const fn view(model: ModelView<'a>) -> Self {
        Self {
            model,
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

    /// Sets the data fields, written in order; a sample's membership and
    /// place fields still follow them. Writing fails with `Unsupported`
    /// for a field SDF cannot frame (an empty, untrimmed, or bracketed name,
    /// or a value with a carriage return, blank line, or `$$$$` line) and
    /// with `Inconsistent` for a field named like the sample metadata.
    #[must_use]
    pub fn with_data_fields(mut self, data_fields: &'a [SdfDataField]) -> Self {
        self.data_fields = data_fields;
        self
    }
}

impl<'a> From<&'a Model> for SdfRecordSource<'a> {
    fn from(model: &'a Model) -> Self {
        Self::model(model)
    }
}

impl<'a> From<FrameView<'a>> for SdfRecordSource<'a> {
    fn from(frame: FrameView<'a>) -> Self {
        Self::view(frame.as_model_view())
    }
}

impl<'a, R: ?Sized> From<SampleView<'a, R>> for SdfRecordSource<'a> {
    fn from(sample: SampleView<'a, R>) -> Self {
        Self::view(sample.as_model_view())
    }
}

impl<'a> From<&'a SdfRecordInterpretation> for SdfRecordSource<'a> {
    fn from(record: &'a SdfRecordInterpretation) -> Self {
        Self::model(record.model())
            .with_title(record.title())
            .with_data_fields(record.data_fields())
    }
}

/// Writes records to a string. See [`write_sdf_to`].
pub(crate) fn write_sdf<'a, R: Into<SdfRecordSource<'a>>>(
    records: impl IntoIterator<Item = R>,
    options: SdfWriteOptions,
) -> Result<String, FormatError> {
    let mut output = Vec::new();
    write_sdf_to(&mut output, records, options)?;
    Ok(String::from_utf8(output).expect("SDF writer emits UTF-8"))
}

/// Writes independent records in input order. A sample set writes one
/// record per sample (`sdf::write(&set, options)`).
pub(crate) fn write_sdf_to<'a, R: Into<SdfRecordSource<'a>>>(
    writer: &mut impl Write,
    records: impl IntoIterator<Item = R>,
    options: SdfWriteOptions,
) -> Result<(), FormatError> {
    for (number, record) in (1..).zip(records) {
        write_record(writer, record.into(), options).map_err(|error| error.in_record(number))?;
    }
    Ok(())
}

fn write_record(
    writer: &mut impl Write,
    record: SdfRecordSource<'_>,
    options: SdfWriteOptions,
) -> Result<(), FormatError> {
    validate_title(record.title)?;
    check_structure_state(
        Format::Sdf,
        record.model.kind(),
        StatePolicies {
            policy: options.sample_metadata,
            writes_membership: true,
            writes_place: true,
        },
    )?;
    // A sample's membership and its step or coordinate place, which travel
    // as data fields.
    let (membership, place) = match record.model.kind() {
        StructureKind::Sample(context) => (
            context.membership(),
            context
                .place()
                .filter(|place| matches!(place, Place::Step(_) | Place::Coordinate(_))),
        ),
        _ => (None, None),
    };
    // Whether a sample's own data uses the name.
    let reserved = |name: &str| {
        (membership.is_some() && SAMPLE_FIELDS.contains(&name))
            || (place.is_some() && PLACE_FIELDS.contains(&name))
    };
    // A coordinate unit travels as the symbol of a unit constant, which the
    // reader turns back into that unit.
    let written_place = match place {
        Some(Place::Coordinate(coordinate))
            if named_unit(&coordinate.unit().to_string()) != Some(coordinate.unit()) =>
        {
            if options.sample_metadata == UnrepresentablePolicy::Fail {
                return Err(FormatError::write(
                    Format::Sdf,
                    DiagnosticKind::Unsupported,
                    format!(
                        "SDF cannot represent the coordinate unit {}.",
                        coordinate.unit()
                    ),
                )
                .with_detail(FormatErrorDetail::UnrepresentableState(
                    StructureField::OrderCoordinate,
                )));
            }
            None
        }
        place => place,
    };
    let mut own_fields = membership.map_or_else(Vec::new, sample_fields);
    own_fields.extend(written_place.map_or_else(Vec::new, place_fields));
    for field in record.data_fields {
        validate_data_field(field)?;
        // Every sample owns these names, whether or not it writes them, so a
        // caller's field never adds a second representative, energy, or
        // place.
        if reserved(field.name()) {
            return Err(FormatError::write(
                Format::Sdf,
                DiagnosticKind::Inconsistent,
                format!(
                    "The data field {:?} is reserved for the data of the sample.",
                    field.name()
                ),
            )
            .with_detail(FormatErrorDetail::DataField {
                name: field.name().to_owned(),
            }));
        }
    }
    let molfile = |error: FormatError| error.in_format(Format::Sdf);
    check_cell(record.model, options.cell).map_err(molfile)?;
    let ctab = render_model(record.model, record.title, options.version).map_err(molfile)?;
    let sink = |error| FormatError::sink(Format::Sdf, error);
    writer.write_all(ctab.as_bytes()).map_err(sink)?;
    for field in record.data_fields.iter().chain(&own_fields) {
        writeln!(writer, ">  <{}>\n{}\n", field.name(), field.value()).map_err(sink)?;
    }
    writer.write_all(b"$$$$\n").map_err(sink)
}

/// A sample's step, or its coordinate and the coordinate's unit, as data
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

/// A sample's weight and state, its state's energy on a representative's
/// record, and the distribution's temperature as data fields,
/// numbers in the shortest decimal text that reads back to the same value.
fn sample_fields(membership: Membership) -> Vec<SdfDataField> {
    let mut fields = vec![
        SdfDataField::new(SDF_WEIGHT_FIELD, membership.weight().to_string()),
        SdfDataField::new(SDF_STATE_FIELD, membership.state().to_string()),
    ];
    if membership.is_representative() {
        fields.push(SdfDataField::new(SDF_REPRESENTATIVE_FIELD, "true"));
        if let Some(energy) = membership.state_energy() {
            let kj_per_mol = energy
                .value_in(KILOJOULE_PER_MOLE)
                .expect("state energies have an energy-per-amount dimension");
            fields.push(SdfDataField::new(SDF_ENERGY_FIELD, kj_per_mol.to_string()));
        }
    }
    if let Some(temperature) = membership.temperature() {
        let kelvin = temperature
            .value_in(KELVIN)
            .expect("distribution temperatures have a temperature dimension");
        fields.push(SdfDataField::new(SDF_TEMPERATURE_FIELD, kelvin.to_string()));
    }
    fields
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
