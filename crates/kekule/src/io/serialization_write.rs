use std::fmt;
use std::io::Write;

use crate::core::Molecule;
use crate::structure::{Ensemble, Model, ModelView};

use super::molfile_write::MolfileRecord;
use super::sdf_document::SdfRecordInterpretation;
use super::v2000::{render_mol_v2000, validate_sdf_data_field, validate_sdf_title};
use super::v3000::render_mol_v3000;
use super::MolWriteError;

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

pub fn write_molfile_molecule(
    molecule: &Molecule,
    options: MolfileWriteOptions,
) -> Result<String, MolWriteError> {
    render_molfile_record(&MolfileRecord::molecule(molecule)?, "", options.version)
}

pub fn write_molfile_model(
    model: ModelView<'_>,
    options: MolfileWriteOptions,
) -> Result<String, MolWriteError> {
    render_molfile_record(&MolfileRecord::model(model)?, "", options.version)
}

pub fn write_molfile_model_to(
    writer: &mut impl Write,
    model: ModelView<'_>,
    options: MolfileWriteOptions,
) -> Result<(), MolWriteError> {
    writer
        .write_all(write_molfile_model(model, options)?.as_bytes())
        .map_err(MolWriteError::io)
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

pub fn write_sdf_model(model: &Model, options: SdfWriteOptions) -> Result<String, SdfWriteError> {
    write_sdf_views(std::iter::once(model.view()), options)
}

pub fn write_sdf_models(
    models: &[Model],
    options: SdfWriteOptions,
) -> Result<String, SdfWriteError> {
    write_sdf_views(models.iter().map(Model::view), options)
}

pub fn write_sdf_ensemble(
    ensemble: &Ensemble,
    options: SdfWriteOptions,
) -> Result<String, SdfWriteError> {
    write_sdf_views(ensemble.members().map(|member| member.as_model()), options)
}

pub fn write_sdf_model_to(
    writer: &mut impl Write,
    model: &Model,
    options: SdfWriteOptions,
) -> Result<(), SdfWriteError> {
    write_sdf_views_to(writer, std::iter::once(model.view()), options)
}

pub fn write_sdf_models_to(
    writer: &mut impl Write,
    models: &[Model],
    options: SdfWriteOptions,
) -> Result<(), SdfWriteError> {
    write_sdf_views_to(writer, models.iter().map(Model::view), options)
}

pub fn write_sdf_ensemble_to(
    writer: &mut impl Write,
    ensemble: &Ensemble,
    options: SdfWriteOptions,
) -> Result<(), SdfWriteError> {
    write_sdf_views_to(
        writer,
        ensemble.members().map(|member| member.as_model()),
        options,
    )
}

pub fn write_sdf_records(
    records: &[SdfRecordInterpretation],
    options: SdfWriteOptions,
) -> Result<String, SdfWriteError> {
    let mut output = Vec::new();
    write_sdf_records_to(&mut output, records, options)?;
    Ok(String::from_utf8(output).expect("SDF writer emits UTF-8"))
}

pub fn write_sdf_records_to(
    writer: &mut impl Write,
    records: &[SdfRecordInterpretation],
    options: SdfWriteOptions,
) -> Result<(), SdfWriteError> {
    for record in records {
        validate_title(record.title())?;
        for field in record.data_fields() {
            validate_sdf_data_field(field).map_err(|error| SdfWriteError::InvalidDataField {
                name: field.name().to_owned(),
                message: error.to_string(),
            })?;
        }
        write_sdf_record_to(
            writer,
            record.model().view(),
            record.title(),
            record
                .data_fields()
                .iter()
                .map(|field| (field.name(), field.value())),
            options,
        )?;
    }
    Ok(())
}

fn write_sdf_views<'a>(
    views: impl IntoIterator<Item = ModelView<'a>>,
    options: SdfWriteOptions,
) -> Result<String, SdfWriteError> {
    let mut output = Vec::new();
    write_sdf_views_to(&mut output, views, options)?;
    Ok(String::from_utf8(output).expect("SDF writer emits UTF-8"))
}

fn write_sdf_views_to<'a>(
    writer: &mut impl Write,
    views: impl IntoIterator<Item = ModelView<'a>>,
    options: SdfWriteOptions,
) -> Result<(), SdfWriteError> {
    for view in views {
        write_sdf_record_to(writer, view, "", std::iter::empty(), options)?;
    }
    Ok(())
}

fn write_sdf_record_to<'a>(
    writer: &mut impl Write,
    model: ModelView<'a>,
    title: &str,
    fields: impl IntoIterator<Item = (&'a str, &'a str)>,
    options: SdfWriteOptions,
) -> Result<(), SdfWriteError> {
    let structural = MolfileRecord::model(model)?;
    let ctab = render_molfile_record(&structural, title, options.version)?;
    writer.write_all(ctab.as_bytes()).map_err(sdf_io)?;
    for (name, value) in fields {
        writeln!(writer, ">  <{name}>\n{value}\n").map_err(sdf_io)?;
    }
    writer.write_all(b"$$$$\n").map_err(sdf_io)
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
