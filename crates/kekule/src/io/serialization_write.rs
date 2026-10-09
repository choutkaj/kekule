use std::fmt;
use std::io::Write;

use crate::core::Molecule;
use crate::structure::{AsModelView, Model, ModelView, Realization, RealizationView};

use super::molfile_write::MolfileRecord;
use super::sdf_document::{SdfDataField, SdfRecordInterpretation};
use super::v2000::{render_mol_v2000, validate_sdf_data_field, validate_sdf_title};
use super::v3000::render_mol_v3000;
use super::{MolWriteError, MolfileVersion};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum MolfileWriteVersion {
    #[default]
    Auto,
    V2000,
    V3000,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MolfileWriteOptions {
    pub version: MolfileWriteVersion,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfWriteOptions {
    pub version: MolfileWriteVersion,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SdfWriteError {
    Molfile(MolWriteError),
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

impl From<MolWriteError> for SdfWriteError {
    fn from(error: MolWriteError) -> Self {
        Self::Molfile(error)
    }
}

fn render_molfile_record(
    record: &MolfileRecord<'_>,
    title: &str,
    version: MolfileWriteVersion,
) -> Result<String, MolWriteError> {
    match version {
        MolfileWriteVersion::V2000 => render_mol_v2000(record, title),
        MolfileWriteVersion::V3000 => render_mol_v3000(record, title),
        MolfileWriteVersion::Auto => {
            render_mol_v2000(record, title).or_else(|_| render_mol_v3000(record, title))
        }
    }
}

fn render_molfile_model(
    model: ModelView<'_>,
    title: &str,
    version: MolfileWriteVersion,
) -> Result<String, MolWriteError> {
    match version {
        MolfileWriteVersion::V2000 => {
            render_mol_v2000(&MolfileRecord::model(model, MolfileVersion::V2000)?, title)
        }
        MolfileWriteVersion::V3000 => {
            render_mol_v3000(&MolfileRecord::model(model, MolfileVersion::V3000)?, title)
        }
        // A fallback must rebuild its stereo projection from the V3000
        // coordinates, not reuse a rounded V2000 drawing.
        MolfileWriteVersion::Auto => render_molfile_model(model, title, MolfileWriteVersion::V2000)
            .or_else(|_| render_molfile_model(model, title, MolfileWriteVersion::V3000)),
    }
}

/// What one Molfile CTAB describes: a coordinate-free molecule written with
/// zero coordinates, or a geometry-bearing model.
///
/// Convert from `&Molecule`, `&Model`, or any borrowed model view through
/// [`Self::model`]. Specified stereo requires coordinates.
#[derive(Debug, Clone, Copy)]
pub enum MolfileSource<'a> {
    Molecule(&'a Molecule),
    Model(ModelView<'a>),
}

impl<'a> MolfileSource<'a> {
    pub fn model(model: &'a (impl AsModelView + ?Sized)) -> Self {
        Self::Model(model.as_model_view())
    }
}

impl<'a> From<&'a Molecule> for MolfileSource<'a> {
    fn from(molecule: &'a Molecule) -> Self {
        Self::Molecule(molecule)
    }
}

impl<'a> From<&'a Model> for MolfileSource<'a> {
    fn from(model: &'a Model) -> Self {
        Self::model(model)
    }
}

/// Writes one CTAB. V2000 rounds coordinates to four decimal places in
/// angstroms; V3000 keeps round-trip decimal text. Stereo is projected
/// against the coordinates emitted by the chosen version.
pub fn write_molfile<'a>(
    source: impl Into<MolfileSource<'a>>,
    options: MolfileWriteOptions,
) -> Result<String, MolWriteError> {
    match source.into() {
        MolfileSource::Molecule(molecule) => {
            render_molfile_record(&MolfileRecord::molecule(molecule)?, "", options.version)
        }
        MolfileSource::Model(model) => render_molfile_model(model, "", options.version),
    }
}

pub fn write_molfile_to<'a>(
    writer: &mut impl Write,
    source: impl Into<MolfileSource<'a>>,
    options: MolfileWriteOptions,
) -> Result<(), MolWriteError> {
    writer
        .write_all(write_molfile(source, options)?.as_bytes())
        .map_err(MolWriteError::io)
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

impl<'a, P: Realization> From<RealizationView<'a, P>> for SdfRecordSource<'a> {
    fn from(item: RealizationView<'a, P>) -> Self {
        Self {
            model: item.as_model_view(),
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
            validate_sdf_data_field(field).map_err(|error| SdfWriteError::InvalidDataField {
                name: field.name().to_owned(),
                message: error.to_string(),
            })?;
        }
        let ctab = render_molfile_model(record.model, record.title, options.version)?;
        writer.write_all(ctab.as_bytes()).map_err(sdf_io)?;
        for field in record.data_fields {
            writeln!(writer, ">  <{}>\n{}\n", field.name(), field.value()).map_err(sdf_io)?;
        }
        writer.write_all(b"$$$$\n").map_err(sdf_io)?;
    }
    Ok(())
}

fn validate_title(title: &str) -> Result<(), SdfWriteError> {
    validate_sdf_title(title).map_err(|error| SdfWriteError::InvalidTitle(error.to_string()))
}

fn sdf_io(error: std::io::Error) -> SdfWriteError {
    SdfWriteError::Io {
        kind: error.kind(),
        message: error.to_string(),
    }
}
