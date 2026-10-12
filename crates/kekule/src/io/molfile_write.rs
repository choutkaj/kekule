//! Molfile CTAB writing: prepare one record, then render it as V2000 or V3000.
use std::io::Write;

use crate::core::Molecule;
use crate::structure::{AsModelView, Model, ModelView};

use super::{MoleculeWriteError, MolfileVersion};

mod record;
mod v2000;
mod v3000;

use record::MolfileRecord;

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
) -> Result<String, MoleculeWriteError> {
    match source.into() {
        MolfileSource::Molecule(molecule) => {
            render_molecule(&MolfileRecord::molecule(molecule)?, "", options.version)
        }
        MolfileSource::Model(model) => render_model(model, "", options.version),
    }
}

pub fn write_molfile_to<'a>(
    writer: &mut impl Write,
    source: impl Into<MolfileSource<'a>>,
    options: MolfileWriteOptions,
) -> Result<(), MoleculeWriteError> {
    writer
        .write_all(write_molfile(source, options)?.as_bytes())
        .map_err(MoleculeWriteError::io)
}

fn render_molecule(
    record: &MolfileRecord<'_>,
    title: &str,
    version: MolfileWriteVersion,
) -> Result<String, MoleculeWriteError> {
    match version {
        MolfileWriteVersion::V2000 => v2000::render(record, title),
        MolfileWriteVersion::V3000 => v3000::render(record, title),
        MolfileWriteVersion::Auto => {
            v2000::render(record, title).or_else(|_| v3000::render(record, title))
        }
    }
}

/// Renders one model CTAB with an SDF or Molfile title.
pub(super) fn render_model(
    model: ModelView<'_>,
    title: &str,
    version: MolfileWriteVersion,
) -> Result<String, MoleculeWriteError> {
    match version {
        MolfileWriteVersion::V2000 => {
            v2000::render(&MolfileRecord::model(model, MolfileVersion::V2000)?, title)
        }
        MolfileWriteVersion::V3000 => {
            v3000::render(&MolfileRecord::model(model, MolfileVersion::V3000)?, title)
        }
        // A fallback must rebuild its stereo projection from the V3000
        // coordinates, not reuse a rounded V2000 drawing.
        MolfileWriteVersion::Auto => render_model(model, title, MolfileWriteVersion::V2000)
            .or_else(|_| render_model(model, title, MolfileWriteVersion::V3000)),
    }
}
