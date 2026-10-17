//! SDF writing: independent Molfile records, each with a title and data fields.
use std::io::Write;

use crate::structure::{AsModelView, EnsembleMemberView, Model, ModelView, TrajectoryFrameView};
use crate::units::{Quantity, KILOJOULE_PER_MOLE};

use super::molfile_write::{check_cell, render_model};
use super::sdf_document::{SdfDataField, SdfRecordInterpretation};
use super::{MolfileWriteVersion, UnrepresentablePolicy};
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfWriteOptions {
    version: MolfileWriteVersion,
    cell: UnrepresentablePolicy,
}

impl SdfWriteOptions {
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
}

/// The data field an ensemble member's relative weight is written to and
/// [`sdf::ensemble`](crate::sdf::ensemble) reads it from.
pub const SDF_WEIGHT_FIELD: &str = "kekule.weight";
/// The data field an ensemble member's energy is written to, in kJ/mol.
pub const SDF_ENERGY_FIELD: &str = "kekule.energy.kJ_per_mol";

/// One SDF record: a model view with an optional title and data fields.
///
/// Convert from `&Model`, a borrowed ensemble member or trajectory frame, or
/// `&SdfRecordInterpretation` (which keeps its title and data fields). An
/// ensemble member's weight and energy follow its data fields as
/// [`SDF_WEIGHT_FIELD`] and [`SDF_ENERGY_FIELD`].
#[derive(Debug, Clone, Copy)]
pub struct SdfRecordSource<'a> {
    model: ModelView<'a>,
    title: &'a str,
    data_fields: &'a [SdfDataField],
    member: Option<(f64, Option<Quantity<f64>>)>,
}

impl<'a> SdfRecordSource<'a> {
    pub fn model(model: &'a (impl AsModelView + ?Sized)) -> Self {
        Self {
            model: model.as_model_view(),
            title: "",
            data_fields: &[],
            member: None,
        }
    }

    #[must_use]
    pub fn with_title(mut self, title: &'a str) -> Self {
        self.title = title;
        self
    }

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

impl<'a> From<EnsembleMemberView<'a>> for SdfRecordSource<'a> {
    fn from(member: EnsembleMemberView<'a>) -> Self {
        Self {
            model: member.as_model_view(),
            title: "",
            data_fields: &[],
            member: Some((member.weight(), member.energy())),
        }
    }
}

impl<'a> From<TrajectoryFrameView<'a>> for SdfRecordSource<'a> {
    fn from(frame: TrajectoryFrameView<'a>) -> Self {
        Self {
            model: frame.as_model_view(),
            title: "",
            data_fields: &[],
            member: None,
        }
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
/// writes one record per item (`sdf::write(&ensemble, options)`).
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
    let member_fields = record.member.map_or_else(Vec::new, member_fields);
    for field in record.data_fields {
        validate_data_field(field)?;
        if member_fields
            .iter()
            .any(|member| member.name() == field.name())
        {
            return Err(FormatError::write(
                Format::Sdf,
                DiagnosticKind::Inconsistent,
                "this field carries the ensemble member's own metadata",
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
    for field in record.data_fields.iter().chain(&member_fields) {
        writeln!(writer, ">  <{}>\n{}\n", field.name(), field.value()).map_err(sink)?;
    }
    writer.write_all(b"$$$$\n").map_err(sink)
}

/// An ensemble member's weight and energy as data fields, in the shortest
/// decimal text that reads back to the same value.
fn member_fields((weight, energy): (f64, Option<Quantity<f64>>)) -> Vec<SdfDataField> {
    let mut fields = vec![SdfDataField::new(SDF_WEIGHT_FIELD, weight.to_string())];
    if let Some(energy) = energy {
        let kj_per_mol = energy
            .value_in(KILOJOULE_PER_MOLE)
            .expect("member energies have an energy-per-amount dimension");
        fields.push(SdfDataField::new(SDF_ENERGY_FIELD, kj_per_mol.to_string()));
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
