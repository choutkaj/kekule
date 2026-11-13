//! Molfile CTAB writing: prepare one record, then render it as V2000 or V3000.
use std::io::Write;

use crate::molecule::Molecule;
use crate::structure::{AsStructureView, Structure, StructureView};

use super::write_error::{check_structure_state, Omissions, StatePolicies, MOLFILE_WRITE};
use super::{MolfileVersion, UnrepresentablePolicy, WriteOmission, WriteReport};
use crate::diagnostics::{DiagnosticKind, Format, FormatError};

mod record;
mod v2000;
mod v3000;

use record::MolfileRecord;

/// Connection-table version a Molfile or SDF writer emits.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum MolfileWriteVersion {
    /// V2000 when it can represent the content, otherwise V3000. Only a
    /// [`DiagnosticKind::Unsupported`] failure promotes; any other failure is
    /// returned as is.
    #[default]
    Auto,
    /// Fixed-width V2000, failing on content it cannot represent.
    V2000,
    /// Free-format V3000, with records continued past 80 columns.
    V3000,
}

/// Options of Molfile writing; the default writes
/// [`MolfileWriteVersion::Auto`] and fails on content a Molfile cannot
/// carry.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MolfileWriteOptions {
    version: MolfileWriteVersion,
    cell: UnrepresentablePolicy,
    set_metadata: UnrepresentablePolicy,
    properties: UnrepresentablePolicy,
}

impl MolfileWriteOptions {
    /// Connection-table version to emit.
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
    /// model with a cell fails to write, and under
    /// [`UnrepresentablePolicy::Omit`] the report names it
    /// ([`WriteOmission::Cell`]).
    pub const fn cell(&self) -> UnrepresentablePolicy {
        self.cell
    }

    /// Sets [`Self::cell`].
    #[must_use]
    pub const fn with_cell(mut self, cell: UnrepresentablePolicy) -> Self {
        self.cell = cell;
        self
    }

    /// The state that a Molfile cannot carry of a geometry, whether of a
    /// structure set or a trajectory frame: its membership in the set's
    /// distribution, its state's properties, its place in the set's order (a
    /// frame's time and MD step among them), and its velocities and forces;
    /// and the occupancies and B factors of any model. By default writing a
    /// structure with any of them fails, and under
    /// [`UnrepresentablePolicy::Omit`] the report names each
    /// ([`WriteOmission::StructureState`], [`WriteOmission::Occupancies`],
    /// [`WriteOmission::BFactors`]).
    pub const fn set_metadata(&self) -> UnrepresentablePolicy {
        self.set_metadata
    }

    /// Sets [`Self::set_metadata`].
    #[must_use]
    pub const fn with_set_metadata(mut self, set_metadata: UnrepresentablePolicy) -> Self {
        self.set_metadata = set_metadata;
        self
    }

    /// A structure's per-atom, per-bond, and own properties and the
    /// properties of a geometry's set, which a Molfile cannot carry: by
    /// default a structure with any fails to write, and under
    /// [`UnrepresentablePolicy::Omit`] the report names each
    /// ([`WriteOmission::StructureState`]).
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

/// What one Molfile CTAB describes: a coordinate-free molecule written with
/// zero coordinates, or a geometry-bearing model, with the title written on
/// its first header line.
///
/// Convert from `&Molecule` or `&Structure`, or build one with [`Self::molecule`]
/// or [`Self::structure`], which accepts any borrowed structure view. Specified stereo
/// requires coordinates.
///
/// ```
/// use kekule::molfile::{self, MolfileSource};
/// let ethanol = kekule::smiles::to_molecules("CCO")?.remove(0);
/// let (text, _) = molfile::write(MolfileSource::molecule(&ethanol).with_title("ethanol"), Default::default())?;
/// assert_eq!(molfile::parse_str(&text, Default::default())?.header().title(), "ethanol");
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct MolfileSource<'a> {
    content: MolfileContent<'a>,
    title: &'a str,
}

#[derive(Debug, Clone, Copy)]
enum MolfileContent<'a> {
    Molecule(&'a Molecule),
    /// A model snapshot, written with its coordinates in angstroms.
    Structure(StructureView<'a>),
}

impl<'a> MolfileSource<'a> {
    /// A coordinate-free molecule, written with zero coordinates.
    pub const fn molecule(molecule: &'a Molecule) -> Self {
        Self {
            content: MolfileContent::Molecule(molecule),
            title: "",
        }
    }

    /// A geometry-bearing model, geometry of a structure set, or trajectory frame.
    pub fn structure(structure: &'a (impl AsStructureView + ?Sized)) -> Self {
        Self {
            content: MolfileContent::Structure(structure.as_structure_view()),
            title: "",
        }
    }

    /// Writes `title` on the first header line; empty by default. Writing
    /// fails for a title spanning lines.
    #[must_use]
    pub const fn with_title(mut self, title: &'a str) -> Self {
        self.title = title;
        self
    }
}

impl<'a> From<&'a Molecule> for MolfileSource<'a> {
    fn from(molecule: &'a Molecule) -> Self {
        Self::molecule(molecule)
    }
}

impl<'a> From<&'a Structure> for MolfileSource<'a> {
    fn from(model: &'a Structure) -> Self {
        Self::structure(model)
    }
}

/// Writes one CTAB and reports what it leaves out. V2000 rounds
/// coordinates to four decimal places in angstroms; V3000 keeps round-trip
/// decimal text. Stereo is projected against the coordinates emitted by the
/// chosen version.
pub(crate) fn write_molfile<'a>(
    source: impl Into<MolfileSource<'a>>,
    options: MolfileWriteOptions,
) -> Result<(String, WriteReport), FormatError> {
    let source = source.into();
    if source.title.contains(['\r', '\n']) {
        return Err(MOLFILE_WRITE
            .unsupported("The title contains a line break.")
            .with_detail(crate::diagnostics::FormatErrorDetail::Title));
    }
    match source.content {
        MolfileContent::Molecule(molecule) => {
            let text = render_molecule(
                &MolfileRecord::molecule(molecule)?,
                source.title,
                options.version,
            )?;
            Ok((text, WriteReport::default()))
        }
        MolfileContent::Structure(model) => {
            let policies = StatePolicies::positions_only(options.set_metadata, options.properties);
            check_cell(model, options.cell)?;
            check_structure_state(Format::Molfile, model, policies)?;
            let text = render_structure(model, source.title, options.version)?;
            let mut omissions = Omissions::default();
            if model.cell().is_some() {
                omissions.record(WriteOmission::Cell { structures: 1 });
            }
            omissions.record_structure_state(model, policies);
            Ok((text, omissions.into_report()))
        }
    }
}

pub(crate) fn write_molfile_to<'a>(
    writer: &mut impl Write,
    source: impl Into<MolfileSource<'a>>,
    options: MolfileWriteOptions,
) -> Result<WriteReport, FormatError> {
    let (text, report) = write_molfile(source, options)?;
    writer
        .write_all(text.as_bytes())
        .map_err(|error| MOLFILE_WRITE.io(error))?;
    Ok(report)
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
    model: StructureView<'_>,
    policy: UnrepresentablePolicy,
) -> Result<(), FormatError> {
    if model.cell().is_some() && policy == UnrepresentablePolicy::Fail {
        return Err(MOLFILE_WRITE.unsupported("A Molfile cannot hold a periodic cell."));
    }
    Ok(())
}

/// Renders one model CTAB with an SDF or Molfile title.
pub(super) fn render_structure(
    model: StructureView<'_>,
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
                match render_structure(model, title, MolfileWriteVersion::V2000) {
                    Err(error) if is_representational(&error) => {}
                    rendered => return rendered,
                }
            }
            // The fallback rebuilds its stereo projection from the V3000
            // coordinates rather than reusing a rounded V2000 drawing.
            render_structure(model, title, MolfileWriteVersion::V3000)
        }
    }
}

/// Only a format limit promotes Auto to V3000; any other failure would
/// recur in V3000 and is reported as V2000 found it.
fn is_representational(error: &FormatError) -> bool {
    error.kind() == DiagnosticKind::Unsupported
}
