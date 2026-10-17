//! Molfile CTAB writing: prepare one record, then render it as V2000 or V3000.
use std::io::Write;

use crate::molecule::Molecule;
use crate::structure::{AsModelView, Model, ModelView};

use super::write_error::MOLFILE_WRITE;
use super::{MolfileVersion, UnrepresentablePolicy};
use crate::diagnostics::{DiagnosticKind, FormatError};

mod record;
mod v2000;
mod v3000;

use record::MolfileRecord;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum MolfileWriteVersion {
    /// V2000 when it can represent the content, otherwise V3000. Only an
    /// [`UnsupportedRepresentation`](DiagnosticKind::Unsupported)
    /// failure promotes; any other failure is returned as is.
    #[default]
    Auto,
    /// Fixed-width V2000, failing on content it cannot represent.
    V2000,
    /// Free-format V3000, with records continued past 80 columns.
    V3000,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MolfileWriteOptions {
    version: MolfileWriteVersion,
    cell: UnrepresentablePolicy,
}

impl MolfileWriteOptions {
    pub const fn version(&self) -> MolfileWriteVersion {
        self.version
    }

    /// Sets [`Self::version`].
    #[must_use]
    pub const fn with_version(mut self, version: MolfileWriteVersion) -> Self {
        self.version = version;
        self
    }

    /// A model's periodic cell, which a Molfile cannot carry: by default a
    /// model with a cell fails to write.
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

/// What one Molfile CTAB describes: a coordinate-free molecule written with
/// zero coordinates, or a geometry-bearing model.
///
/// Convert from `&Molecule`, `&Model`, or any borrowed model view through
/// [`Self::model`]. Specified stereo requires coordinates.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
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
pub(crate) fn write_molfile<'a>(
    source: impl Into<MolfileSource<'a>>,
    options: MolfileWriteOptions,
) -> Result<String, FormatError> {
    match source.into() {
        MolfileSource::Molecule(molecule) => {
            render_molecule(&MolfileRecord::molecule(molecule)?, "", options.version)
        }
        MolfileSource::Model(model) => {
            check_cell(model, options.cell)?;
            render_model(model, "", options.version)
        }
    }
}

pub(crate) fn write_molfile_to<'a>(
    writer: &mut impl Write,
    source: impl Into<MolfileSource<'a>>,
    options: MolfileWriteOptions,
) -> Result<(), FormatError> {
    writer
        .write_all(write_molfile(source, options)?.as_bytes())
        .map_err(|error| MOLFILE_WRITE.io(error))
}

fn render_molecule(
    record: &MolfileRecord<'_>,
    title: &str,
    version: MolfileWriteVersion,
) -> Result<String, FormatError> {
    match version {
        MolfileWriteVersion::V2000 => v2000::render(record, title),
        MolfileWriteVersion::V3000 => v3000::render(record, title),
        MolfileWriteVersion::Auto => match v2000::render(record, title) {
            Err(error) if is_representational(&error) => v3000::render(record, title),
            rendered => rendered,
        },
    }
}

/// Fails for a model with a periodic cell, which no CTAB can carry, unless
/// the caller chose to omit it.
pub(super) fn check_cell(
    model: ModelView<'_>,
    policy: UnrepresentablePolicy,
) -> Result<(), FormatError> {
    if model.cell().is_some() && policy == UnrepresentablePolicy::Fail {
        return Err(MOLFILE_WRITE.unsupported(
            "a Molfile cannot represent a periodic cell; set the cell option to omit it",
        ));
    }
    Ok(())
}

/// Renders one model CTAB with an SDF or Molfile title.
pub(super) fn render_model(
    model: ModelView<'_>,
    title: &str,
    version: MolfileWriteVersion,
) -> Result<String, FormatError> {
    match version {
        MolfileWriteVersion::V2000 => {
            v2000::render(&MolfileRecord::model(model, MolfileVersion::V2000)?, title)
        }
        MolfileWriteVersion::V3000 => {
            v3000::render(&MolfileRecord::model(model, MolfileVersion::V3000)?, title)
        }
        MolfileWriteVersion::Auto => {
            // Size and enhanced stereo are known before any stereo projection,
            // so content V2000 cannot hold skips the V2000 attempt entirely.
            let topology = model.topology();
            let fits = v2000::check_capacity(
                topology.atom_count(),
                topology.bond_count(),
                topology
                    .instances()
                    .any(|occurrence| occurrence.molecule().stereo_groups().next().is_some()),
            );
            if fits.is_ok() {
                match render_model(model, title, MolfileWriteVersion::V2000) {
                    Err(error) if is_representational(&error) => {}
                    rendered => return rendered,
                }
            }
            // The fallback rebuilds its stereo projection from the V3000
            // coordinates rather than reusing a rounded V2000 drawing.
            render_model(model, title, MolfileWriteVersion::V3000)
        }
    }
}

/// Only a format limit promotes Auto to V3000; any other failure would
/// recur in V3000 and is reported as V2000 found it.
fn is_representational(error: &FormatError) -> bool {
    error.kind() == DiagnosticKind::Unsupported
}
