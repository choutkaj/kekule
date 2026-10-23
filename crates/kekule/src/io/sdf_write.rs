//! SDF writing: independent Molfile records, each with a title and data fields.
use std::io::Write;

use crate::structure::{
    AsModelView, EnsembleSampleView, Model, ModelView, StructureKind, TrajectoryFrameView,
};
use crate::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE};

use super::molfile_write::{check_cell, render_model};
use super::sdf_document::{SdfDataField, SdfRecordInterpretation};
use super::write_error::check_structure_state;
use super::{MolfileWriteVersion, UnrepresentablePolicy};
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};

/// Options of SDF writing; the default writes
/// [`MolfileWriteVersion::Auto`] records and fails on a periodic cell.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfWriteOptions {
    version: MolfileWriteVersion,
    cell: UnrepresentablePolicy,
    frame_metadata: UnrepresentablePolicy,
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

    /// A trajectory frame's time, step, velocities, and forces, which an SDF
    /// record cannot carry: by default writing a frame with any of them
    /// fails. An ensemble sample's weight, state, state energy, and
    /// temperature are always written as data fields.
    pub const fn frame_metadata(&self) -> UnrepresentablePolicy {
        self.frame_metadata
    }

    /// Sets [`Self::frame_metadata`].
    #[must_use]
    pub const fn with_frame_metadata(mut self, frame_metadata: UnrepresentablePolicy) -> Self {
        self.frame_metadata = frame_metadata;
        self
    }
}

/// The data field an ensemble sample's relative weight is written to and
/// [`sdf::ensemble`](crate::sdf::ensemble) reads it from.
pub const SDF_WEIGHT_FIELD: &str = "kekule.weight";
/// The data field the index of an ensemble sample's state is written to.
pub const SDF_STATE_FIELD: &str = "kekule.state";
/// The data field, `true`, that marks the record of a state's representative
/// sample.
pub const SDF_REPRESENTATIVE_FIELD: &str = "kekule.representative";
/// The data field a state's energy is written to, in kJ/mol, on the record of
/// its representative sample.
pub const SDF_ENERGY_FIELD: &str = "kekule.energy.kJ_per_mol";
/// The data field the temperature an ensemble's populations describe is
/// written to, in kelvin, on every record of the ensemble.
pub const SDF_TEMPERATURE_FIELD: &str = "kekule.temperature.K";

/// The data fields an ensemble sample's own metadata uses.
const SAMPLE_FIELDS: [&str; 5] = [
    SDF_WEIGHT_FIELD,
    SDF_STATE_FIELD,
    SDF_REPRESENTATIVE_FIELD,
    SDF_ENERGY_FIELD,
    SDF_TEMPERATURE_FIELD,
];

/// One SDF record: a model view with an optional title and data fields.
///
/// Convert from `&Model`, a borrowed ensemble sample or trajectory frame, or
/// `&SdfRecordInterpretation` (which keeps its title and data fields). An
/// ensemble sample's weight, its state, and the ensemble's temperature follow
/// its data fields as [`SDF_WEIGHT_FIELD`], [`SDF_STATE_FIELD`], and
/// [`SDF_TEMPERATURE_FIELD`]; a representative's record is marked with
/// [`SDF_REPRESENTATIVE_FIELD`] and carries its state's energy as
/// [`SDF_ENERGY_FIELD`], whichever way its view is passed. State and
/// ensemble properties are not written. A trajectory frame's time, step,
/// velocities, and forces fail to write unless
/// [`SdfWriteOptions::frame_metadata`] omits them.
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

    /// Sets the data fields, written in order; an ensemble sample's fields
    /// still follow them. Writing fails with `Unsupported`
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

impl<'a> From<EnsembleSampleView<'a>> for SdfRecordSource<'a> {
    fn from(sample: EnsembleSampleView<'a>) -> Self {
        Self::view(sample.as_model_view())
    }
}

impl<'a> From<TrajectoryFrameView<'a>> for SdfRecordSource<'a> {
    fn from(frame: TrajectoryFrameView<'a>) -> Self {
        Self::view(frame.as_model_view())
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

/// Writes independent records in input order. An ensemble or trajectory
/// writes one record per sample or frame (`sdf::write(&ensemble, options)`).
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
        None,
        options.frame_metadata,
    )?;
    let sample = matches!(record.model.kind(), StructureKind::Sample { .. });
    let sample_fields = match record.model.kind() {
        StructureKind::Sample {
            weight,
            state,
            representative,
            energy,
            temperature,
            ..
        } => sample_fields(weight, state, representative, energy, temperature),
        _ => Vec::new(),
    };
    for field in record.data_fields {
        validate_data_field(field)?;
        // Every sample owns these names, whether or not it writes them, so
        // a caller's field never adds a second representative or energy.
        if sample && SAMPLE_FIELDS.contains(&field.name()) {
            return Err(FormatError::write(
                Format::Sdf,
                DiagnosticKind::Inconsistent,
                "this field carries the ensemble sample's own metadata",
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
    for field in record.data_fields.iter().chain(&sample_fields) {
        writeln!(writer, ">  <{}>\n{}\n", field.name(), field.value()).map_err(sink)?;
    }
    writer.write_all(b"$$$$\n").map_err(sink)
}

/// An ensemble sample's weight and state, its state's energy on a
/// representative's record, and the ensemble's temperature as data fields,
/// numbers in the shortest decimal text that reads back to the same value.
fn sample_fields(
    weight: f64,
    state: usize,
    representative: bool,
    energy: Option<&Quantity<f64>>,
    temperature: Option<&Quantity<f64>>,
) -> Vec<SdfDataField> {
    let mut fields = vec![
        SdfDataField::new(SDF_WEIGHT_FIELD, weight.to_string()),
        SdfDataField::new(SDF_STATE_FIELD, state.to_string()),
    ];
    if representative {
        fields.push(SdfDataField::new(SDF_REPRESENTATIVE_FIELD, "true"));
        if let Some(energy) = energy {
            let kj_per_mol = energy
                .value_in(KILOJOULE_PER_MOLE)
                .expect("state energies have an energy-per-amount dimension");
            fields.push(SdfDataField::new(SDF_ENERGY_FIELD, kj_per_mol.to_string()));
        }
    }
    if let Some(temperature) = temperature {
        let kelvin = temperature
            .value_in(KELVIN)
            .expect("ensemble temperatures have a temperature dimension");
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
        return invalid("SDF record titles cannot contain line breaks");
    }
    if title.trim() == "$$$$" {
        return invalid("SDF record titles cannot be a record delimiter");
    }
    Ok(())
}

fn validate_data_field(field: &SdfDataField) -> Result<(), FormatError> {
    let invalid = |message: &str| {
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
        return invalid(
            "SDF data field names must be nonempty, trimmed, and exclude angle brackets or line breaks",
        );
    }
    let value = field.value();
    if value.contains('\r') {
        return invalid("SDF data field values cannot contain carriage returns");
    }
    if !value.is_empty() && value.split('\n').any(str::is_empty) {
        return invalid("SDF data field values cannot contain blank lines");
    }
    if value.lines().any(|line| line.trim() == "$$$$") {
        return invalid("SDF data field values cannot contain a record delimiter line");
    }
    Ok(())
}
