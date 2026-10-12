//! SDF writing: independent Molfile records, each with a title and data fields.
use std::fmt;
use std::io::Write;

use crate::structure::{AsModelView, EnsembleMemberView, Model, ModelView, TrajectoryFrameView};

use super::molfile_write::render_model;
use super::sdf_document::{SdfDataField, SdfRecordInterpretation};
use super::{MoleculeWriteError, MolfileWriteVersion};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfWriteOptions {
    pub version: MolfileWriteVersion,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SdfWriteError {
    Molfile(MoleculeWriteError),
    InvalidTitle(String),
    InvalidDataField {
        name: String,
        message: String,
    },
    Io {
        kind: std::io::ErrorKind,
        message: String,
    },
}

impl fmt::Display for SdfWriteError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Molfile(error) => write!(formatter, "{error}"),
            Self::InvalidTitle(message) => write!(formatter, "invalid SDF title: {message}"),
            Self::InvalidDataField { name, message } => {
                write!(formatter, "invalid SDF data field `{name}`: {message}")
            }
            Self::Io { message, .. } => write!(formatter, "SDF output failed: {message}"),
        }
    }
}

impl std::error::Error for SdfWriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Molfile(error) => Some(error),
            Self::InvalidTitle(_) | Self::InvalidDataField { .. } | Self::Io { .. } => None,
        }
    }
}

impl From<MoleculeWriteError> for SdfWriteError {
    fn from(error: MoleculeWriteError) -> Self {
        Self::Molfile(error)
    }
}

/// One SDF record: a model view with an optional title and data fields.
///
/// Convert from `&Model`, a borrowed ensemble member or trajectory frame, or
/// `&SdfRecordInterpretation` (which keeps its title and data fields).
#[derive(Debug, Clone, Copy)]
pub struct SdfRecordSource<'a> {
    model: ModelView<'a>,
    title: &'a str,
    data_fields: &'a [SdfDataField],
}

impl<'a> SdfRecordSource<'a> {
    pub fn model(model: &'a (impl AsModelView + ?Sized)) -> Self {
        Self {
            model: model.as_model_view(),
            title: "",
            data_fields: &[],
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
        }
    }
}

impl<'a> From<TrajectoryFrameView<'a>> for SdfRecordSource<'a> {
    fn from(frame: TrajectoryFrameView<'a>) -> Self {
        Self {
            model: frame.as_model_view(),
            title: "",
            data_fields: &[],
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
pub fn write_sdf<'a, R: Into<SdfRecordSource<'a>>>(
    records: impl IntoIterator<Item = R>,
    options: SdfWriteOptions,
) -> Result<String, SdfWriteError> {
    let mut output = Vec::new();
    write_sdf_to(&mut output, records, options)?;
    Ok(String::from_utf8(output).expect("SDF writer emits UTF-8"))
}

/// Writes independent records in input order. An ensemble or trajectory
/// writes one record per item (`sdf::write(&ensemble, options)`).
pub fn write_sdf_to<'a, R: Into<SdfRecordSource<'a>>>(
    writer: &mut impl Write,
    records: impl IntoIterator<Item = R>,
    options: SdfWriteOptions,
) -> Result<(), SdfWriteError> {
    for record in records {
        let record = record.into();
        validate_title(record.title)?;
        for field in record.data_fields {
            validate_data_field(field)?;
        }
        let ctab = render_model(record.model, record.title, options.version)?;
        writer.write_all(ctab.as_bytes()).map_err(sdf_io)?;
        for field in record.data_fields {
            writeln!(writer, ">  <{}>\n{}\n", field.name(), field.value()).map_err(sdf_io)?;
        }
        writer.write_all(b"$$$$\n").map_err(sdf_io)?;
    }
    Ok(())
}

fn validate_title(title: &str) -> Result<(), SdfWriteError> {
    let invalid = |message: &str| Err(SdfWriteError::InvalidTitle(message.to_owned()));
    if title.contains(['\r', '\n']) {
        return invalid("SDF record titles cannot contain line breaks");
    }
    if title.trim() == "$$$$" {
        return invalid("SDF record titles cannot be a record delimiter");
    }
    Ok(())
}

fn validate_data_field(field: &SdfDataField) -> Result<(), SdfWriteError> {
    let invalid = |message: &str| {
        Err(SdfWriteError::InvalidDataField {
            name: field.name().to_owned(),
            message: message.to_owned(),
        })
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

fn sdf_io(error: std::io::Error) -> SdfWriteError {
    SdfWriteError::Io {
        kind: error.kind(),
        message: error.to_string(),
    }
}
