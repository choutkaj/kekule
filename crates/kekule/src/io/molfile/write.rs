//! Molfile CTAB writing: prepare one record, then render it as V2000 or V3000.
use std::io::Write;

use crate::molecule::Molecule;
use crate::structure::{AsStructureView, Structure, StructureView};

use crate::diagnostics::{DiagnosticKind, Format, FormatError};
use crate::io::write_policy::{check_geometry_state, Carries, StatePolicies, MOLFILE_WRITE};
use crate::io::{MolfileVersion, UnrepresentablePolicy, WriteReport};

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
    geometry_metadata: UnrepresentablePolicy,
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

    /// A structure's periodic cell, which a Molfile cannot carry: by default a
    /// structure with a cell fails to write, and under
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

    /// The state that a Molfile cannot carry: a geometry's place in its set's
    /// order (a frame's time and MD step among them), and the velocities,
    /// forces, occupancies, and B factors of a geometry or structure. By default writing a
    /// structure with any of them fails, and under
    /// [`UnrepresentablePolicy::Omit`] the report names each
    /// ([`WriteOmission::GeometryState`], [`WriteOmission::Occupancies`],
    /// [`WriteOmission::BFactors`]).
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
    /// properties of a geometry's set, which a Molfile cannot carry: by
    /// default a structure with any fails to write, and under
    /// [`UnrepresentablePolicy::Omit`] the report names each
    /// ([`WriteOmission::GeometryState`]).
    pub const fn properties(&self) -> UnrepresentablePolicy {
        self.properties
    }

    /// Sets [`Self::properties`].
    #[must_use]
    pub const fn with_properties(mut self, properties: UnrepresentablePolicy) -> Self {
        self.properties = properties;
        self
    }

    /// The policies of a record whose format carries `carries`.
    pub(in crate::io) const fn state_policies(&self, carries: Carries) -> StatePolicies {
        StatePolicies {
            carries,
            geometry_metadata: self.geometry_metadata,
            properties: self.properties,
            cell: self.cell,
        }
    }
}

/// What one Molfile CTAB describes: a coordinate-free molecule written with
/// zero coordinates, or a structure, with the title written on its first
/// header line.
///
/// Convert from `&Molecule` or `&Structure`, or build one with [`Self::molecule`]
/// or [`Self::structure`], which accepts any borrowed structure view. Specified stereo
/// requires coordinates. The writer never invents stereo assertions from an
/// unasserted drawing and never hides a conflict between coordinates and
/// assertions.
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
    /// A structure snapshot, written with its coordinates in angstroms.
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

    /// A structure, geometry of a structure set, or trajectory frame.
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
    fn from(structure: &'a Structure) -> Self {
        Self::structure(structure)
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
            Ok((text, WriteReport::empty()))
        }
        MolfileContent::Structure(structure) => {
            let omitted = check_geometry_state(
                Format::Molfile,
                structure,
                options.state_policies(Carries::POSITIONS),
            )?;
            let text = render_structure(structure, source.title, options.version)?;
            Ok((text, omitted.into_iter().collect()))
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

/// Renders one structure CTAB with an SDF or Molfile title.
pub(in crate::io) fn render_structure(
    structure: StructureView<'_>,
    title: &str,
    version: MolfileWriteVersion,
) -> Result<String, FormatError> {
    match version {
        MolfileWriteVersion::V2000 => v2000::render(
            &MolfileRecord::structure(structure, MolfileVersion::V2000)?,
            title,
        ),
        MolfileWriteVersion::V3000 => v3000::render(
            &MolfileRecord::structure(structure, MolfileVersion::V3000)?,
            title,
        ),
        MolfileWriteVersion::Auto => {
            // Size and enhanced stereo are known before any stereo projection,
            // so content V2000 cannot hold skips the V2000 attempt entirely.
            let topology = structure.topology();
            let fits = v2000::check_capacity(
                topology.atom_count(),
                topology.bond_count(),
                topology
                    .instances()
                    .any(|occurrence| occurrence.molecule().stereo_groups().next().is_some()),
            );
            if fits.is_ok() {
                match render_structure(structure, title, MolfileWriteVersion::V2000) {
                    Err(error) if is_representational(&error) => {}
                    rendered => return rendered,
                }
            }
            // The fallback rebuilds its stereo projection from the V3000
            // coordinates rather than reusing a rounded V2000 drawing.
            render_structure(structure, title, MolfileWriteVersion::V3000)
        }
    }
}

/// Only a format limit promotes Auto to V3000; any other failure would
/// recur in V3000 and is reported as V2000 found it.
fn is_representational(error: &FormatError) -> bool {
    error.kind() == DiagnosticKind::Unsupported
}
