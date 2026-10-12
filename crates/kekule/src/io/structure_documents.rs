use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::chemistry::{
    canonicalize_molecule_for_publication, localize_source_aromatic_bonds, NormalizationError,
    NormalizationWarning, SourceStereoBondMark, SourceStereoBondMarkKind,
};
use crate::core::{
    Atom, AtomId, AtomRadical, BondId, BondOrder, Element, ImplicitHydrogens, Molecule,
    MoleculeEditor, StereoElement, StereoElementId, StereoElementKind, StereoGroup,
    StereoGroupKind, TetrahedralOrientation, TetrahedralStereo,
};
use crate::diagnostics::{DiagnosticKind, SourceLocation, UnsupportedContent};
use crate::geometry::Point3;
use crate::structure::{Model, ModelBuildError, ModelBuilder, Positions};
use crate::topology::Topology;
use crate::units::{Quantity, ANGSTROM};

use super::molfile_document::{
    MolfileAtom, MolfileAtomParity, MolfileBondOrder, MolfileBondStereo, MolfileDocument,
    MolfileRadical, MolfileStereoGroup, MolfileStereoGroupKind, MolfileUnsupported, MolfileVersion,
};
use super::staged_coordinates::StagedCoordinates;

/// Policy for Molfile and SDF record interpretation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MolfileInterpretOptions {
    /// What to do with [`MolfileDocument::unsupported`] content. It fails
    /// interpretation by default; [`UnsupportedContent::Report`] omits it and
    /// lists it in [`MolfileInterpretation::omitted`].
    pub unsupported: UnsupportedContent,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileInterpretError {
    kind: DiagnosticKind,
    location: Option<SourceLocation>,
    message: String,
}

impl MolfileInterpretError {
    pub(super) fn new(
        kind: DiagnosticKind,
        location: Option<SourceLocation>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            kind,
            location,
            message: message.into(),
        }
    }

    fn at(kind: DiagnosticKind, location: SourceLocation, message: impl Into<String>) -> Self {
        Self::new(kind, Some(location), message)
    }

    fn chemistry(error: impl fmt::Display) -> Self {
        Self::new(DiagnosticKind::Chemistry, None, error.to_string())
    }

    pub const fn kind(&self) -> DiagnosticKind {
        self.kind
    }

    /// The source line the failure concerns, when one does.
    pub const fn location(&self) -> Option<SourceLocation> {
        self.location
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for MolfileInterpretError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.location {
            Some(location) => write!(
                f,
                "Molfile interpretation error at {location}: {}",
                self.message
            ),
            None => write!(f, "Molfile interpretation error: {}", self.message),
        }
    }
}

impl std::error::Error for MolfileInterpretError {}

fn molfile_atom(atom: &MolfileAtom) -> Result<Atom, MolfileInterpretError> {
    let invalid = |message: String| {
        MolfileInterpretError::at(DiagnosticKind::InvalidValue, atom.location, message)
    };
    let element = Element::from_symbol(&atom.symbol)
        .ok_or_else(|| invalid(format!("unsupported element symbol `{}`", atom.symbol)))?;
    let mut value = Atom::new(element);
    value.formal_charge = i8::try_from(atom.formal_charge)
        .map_err(|_| invalid("formal charge is outside i8 range".to_owned()))?;
    value.isotope = atom.isotope;
    value.radical = atom.radical.and_then(|radical| match radical {
        MolfileRadical::Singlet => AtomRadical::new(2, Some(1)),
        MolfileRadical::Doublet => AtomRadical::new(1, Some(2)),
        MolfileRadical::Triplet => AtomRadical::new(2, Some(3)),
    });
    value.hydrogens = if atom.hydrogen_count.is_some() || atom.valence.is_some() {
        ImplicitHydrogens::Fixed(atom.hydrogen_count.unwrap_or(0))
    } else {
        ImplicitHydrogens::Inferred
    };
    value.atom_map = atom.atom_map;
    Ok(value)
}

fn apply_molfile_declared_valence(
    molecule: &mut Molecule,
    atom: AtomId,
    declaration: &MolfileAtom,
    declared_valence: u8,
    source_aromatic_bonds: &BTreeSet<BondId>,
    version: MolfileVersion,
) -> Result<(), MolfileInterpretError> {
    let version_name = version.name();
    let location = declaration.location;
    let unsupported_order = || {
        MolfileInterpretError::at(
            DiagnosticKind::Unsupported,
            location,
            format!(
                "declared {version_name} valence cannot determine hydrogens for this bond order"
            ),
        )
    };
    let represented_valence = molecule
        .incident_bonds(atom)
        .map_err(|error| {
            MolfileInterpretError::at(DiagnosticKind::Chemistry, location, error.to_string())
        })?
        .try_fold(0usize, |total, (bond, value)| {
            if source_aromatic_bonds.contains(&bond) {
                return Err(unsupported_order());
            }
            let contribution = match value.order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Zero | BondOrder::Dative => 0,
                BondOrder::Quadruple => return Err(unsupported_order()),
            };
            Ok(total + contribution)
        })?;
    let inferred_hydrogens = usize::from(declared_valence)
        .checked_sub(represented_valence)
        .ok_or_else(|| {
            MolfileInterpretError::at(
                DiagnosticKind::Inconsistent,
                location,
                format!("declared {version_name} valence is below represented bond valence"),
            )
        })?;
    if let Some(declared_hydrogens) = declaration.hydrogen_count {
        if inferred_hydrogens != usize::from(declared_hydrogens) {
            let message = match version {
                MolfileVersion::V2000 => "V2000 hydrogen-count and valence declarations conflict",
                MolfileVersion::V3000 => "V3000 HCOUNT and VAL declarations conflict",
            };
            return Err(MolfileInterpretError::at(
                DiagnosticKind::Inconsistent,
                location,
                message,
            ));
        }
        return Ok(());
    }

    let explicit = u8::try_from(inferred_hydrogens).map_err(|_| {
        MolfileInterpretError::at(
            DiagnosticKind::InvalidValue,
            location,
            format!("declared {version_name} hydrogen count exceeds u8"),
        )
    })?;
    molecule
        .atom_mut(atom)
        .expect("interpreted Molfile atom remains live")
        .hydrogens = ImplicitHydrogens::Fixed(explicit);
    Ok(())
}

/// The represented graph before component partitioning and publication.
struct StagedMolfile {
    editor: MoleculeEditor,
    geometry: StagedCoordinates,
    source_stereo: Vec<SourceStereoBondMark>,
}

/// Builds the represented graph from decoded atoms and bonds in source order,
/// so staging IDs index the document's atom and bond blocks.
fn stage_molfile(document: &MolfileDocument) -> Result<StagedMolfile, MolfileInterpretError> {
    let mut editor = MoleculeEditor::new();
    let mut geometry = StagedCoordinates::with_atom_capacity(document.atoms().len(), ANGSTROM)
        .map_err(MolfileInterpretError::chemistry)?;
    let mut atom_ids = Vec::with_capacity(document.atoms().len());
    for atom in document.atoms() {
        let chemistry = |error: String| {
            MolfileInterpretError::at(DiagnosticKind::Chemistry, atom.location, error)
        };
        let id = editor
            .add_atom(molfile_atom(atom)?)
            .map_err(|error| chemistry(format!("invalid graph atom: {error}")))?;
        geometry
            .set_position(id, Quantity::new(atom.position, ANGSTROM))
            .map_err(|error| chemistry(error.to_string()))?;
        atom_ids.push(id);
    }
    let mut source_aromatic_bonds = BTreeSet::new();
    let mut source_stereo = Vec::new();
    let mut first_aromatic = None;
    for bond in document.bonds() {
        let [a, b] = bond.atoms.map(|atom| atom_ids[atom]);
        let order = match bond.order {
            MolfileBondOrder::Zero => BondOrder::Zero,
            MolfileBondOrder::Single | MolfileBondOrder::Aromatic => BondOrder::Single,
            MolfileBondOrder::Double => BondOrder::Double,
            MolfileBondOrder::Triple => BondOrder::Triple,
            MolfileBondOrder::Dative => BondOrder::Dative,
        };
        let id = editor.add_bond(a, b, order).map_err(|error| {
            MolfileInterpretError::at(
                DiagnosticKind::Chemistry,
                bond.location,
                format!("invalid graph bond: {error}"),
            )
        })?;
        if bond.order == MolfileBondOrder::Aromatic {
            source_aromatic_bonds.insert(id);
            first_aromatic.get_or_insert(bond.location);
        }
        if let Some(stereo) = bond.stereo {
            let kind = match (bond.order, stereo) {
                (MolfileBondOrder::Double, _) => SourceStereoBondMarkKind::DoubleBondEither,
                (_, MolfileBondStereo::Up) => SourceStereoBondMarkKind::WedgeUp,
                (_, MolfileBondStereo::Down) => SourceStereoBondMarkKind::WedgeDown,
                (_, MolfileBondStereo::Either) => SourceStereoBondMarkKind::WedgeEither,
            };
            source_stereo.push(SourceStereoBondMark {
                bond: id,
                from: a,
                kind,
            });
        }
    }
    for (atom, &id) in document.atoms().iter().zip(&atom_ids) {
        if let Some(valence) = atom.valence {
            apply_molfile_declared_valence(
                editor.working_mut(),
                id,
                atom,
                valence,
                &source_aromatic_bonds,
                document.version(),
            )?;
        }
    }
    localize_source_aromatic_bonds(editor.working_mut(), &source_aromatic_bonds).map_err(
        |error| {
            MolfileInterpretError::new(
                DiagnosticKind::Chemistry,
                first_aromatic,
                format!("aromatic bond localization failed: {error}"),
            )
        },
    )?;
    Ok(StagedMolfile {
        editor,
        geometry,
        source_stereo,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolfileAtomMapping {
    atom: AtomId,
    location: SourceLocation,
}

impl MolfileAtomMapping {
    pub const fn atom(&self) -> AtomId {
        self.atom
    }

    /// The source atom line.
    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolfileBondMapping {
    bond: BondId,
    location: SourceLocation,
}

impl MolfileBondMapping {
    pub const fn bond(&self) -> BondId {
        self.bond
    }

    /// The source bond line.
    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// Source correspondence and diagnostics of one interpreted component.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MolfileInterpretationReport {
    atom_mappings: Vec<MolfileAtomMapping>,
    bond_mappings: Vec<MolfileBondMapping>,
    created_stereo_elements: Vec<StereoElementId>,
    warnings: Vec<MolfileInterpretationWarning>,
}

impl MolfileInterpretationReport {
    pub fn atom_mappings(&self) -> &[MolfileAtomMapping] {
        &self.atom_mappings
    }

    pub fn bond_mappings(&self) -> &[MolfileBondMapping] {
        &self.bond_mappings
    }

    /// Canonical stereo elements decoded from source bond marks.
    pub fn created_stereo_elements(&self) -> &[StereoElementId] {
        &self.created_stereo_elements
    }

    pub fn warnings(&self) -> &[MolfileInterpretationWarning] {
        &self.warnings
    }
}

/// Nonfatal source-representation diagnostic produced during interpretation.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum MolfileInterpretationWarning {
    /// Incompatible wedge marks cannot specify an axis configuration. The
    /// source document retains the marks; the interpreted graph omits that
    /// assertion and preserves its connectivity and unrelated stereo.
    ConflictingAtropisomericWedgeMarks {
        axis: BondId,
        location: SourceLocation,
        mark_count: usize,
    },
    /// The represented graph is retained, but its valence is outside the
    /// coordinate-stereo model. No new configuration is inferred. Resource
    /// exhaustion and structural failures remain interpretation errors.
    CoordinateStereoValenceUnsupported {
        error: crate::algorithms::ValenceError,
    },
    /// Wedge geometry does not establish one configuration. A single mark may
    /// be ambiguous; no represented stereo is created from an invalid drawing.
    AmbiguousTetrahedralWedgeMarks {
        center: AtomId,
        location: SourceLocation,
        mark_count: usize,
    },
    /// A V3000 atom `CFG` parity contradicts the configuration of the atom's
    /// 3D coordinates. CTfile parity is derived from coordinates, so, as for
    /// wedges in a conformation, the coordinates decide and the parity is not
    /// asserted.
    AtomParityContradictsCoordinates {
        center: AtomId,
        location: SourceLocation,
    },
}

/// One final canonical model and source reports in molecule-instance order.
///
/// Each connected component is one instance. Dense atom order is the source
/// atom-block order, so components may interleave in the dense order.
#[derive(Debug, Clone, PartialEq)]
pub struct MolfileInterpretation {
    model: Model,
    reports: Vec<MolfileInterpretationReport>,
    omitted: Vec<MolfileUnsupported>,
}

impl MolfileInterpretation {
    pub fn model(&self) -> &Model {
        &self.model
    }

    pub fn topology(&self) -> &Topology {
        self.model.topology()
    }

    /// One component report per molecule instance, in topology instance order.
    /// Atom, bond, and stereo IDs in each report are local to that molecule.
    pub fn reports(&self) -> &[MolfileInterpretationReport] {
        &self.reports
    }

    /// Unsupported source content omitted under
    /// [`UnsupportedContent::Report`]; empty otherwise.
    pub fn omitted(&self) -> &[MolfileUnsupported] {
        &self.omitted
    }

    pub fn molecules(&self) -> impl ExactSizeIterator<Item = &Molecule> + DoubleEndedIterator {
        self.topology()
            .molecules()
            .map(|instance| instance.molecule())
    }

    pub fn into_molecules(self) -> Vec<Molecule> {
        self.molecules().cloned().collect()
    }

    /// Retains the exact topology allocation, including its synthetic hierarchy.
    pub fn into_topology(self) -> std::sync::Arc<Topology> {
        self.model.shared_topology()
    }

    pub fn into_model(self) -> Model {
        self.model
    }

    pub fn into_parts(
        self,
    ) -> (
        Model,
        Vec<MolfileInterpretationReport>,
        Vec<MolfileUnsupported>,
    ) {
        (self.model, self.reports, self.omitted)
    }
}

fn publish_molfile_components(
    components: Vec<(Molecule, Positions, MolfileInterpretationReport)>,
    omitted: Vec<MolfileUnsupported>,
) -> Result<MolfileInterpretation, ModelBuildError> {
    let mut builder = ModelBuilder::new();
    let mut reports = Vec::with_capacity(components.len());
    let mut lines = Vec::new();
    let chain = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_chain("A", None)?;
    for (index, (molecule, positions, report)) in components.into_iter().enumerate() {
        let atom_ids = molecule.atom_ids().collect::<Vec<_>>();
        let instance = builder.add_molecule(molecule, &positions)?;
        let sequence = i32::try_from(index + 1).map_err(|_| ModelBuildError::CapacityOverflow)?;
        let hierarchy = builder.topology_builder_mut().hierarchy_mut();
        let residue = hierarchy.add_residue(
            chain,
            "UNL",
            Some(sequence),
            Some(sequence.to_string()),
            None,
        )?;
        for atom in atom_ids {
            hierarchy.add_atom_site(
                residue,
                crate::topology::InstanceAtomId::new(instance, atom),
                crate::topology::AtomSiteMetadata::default(),
            )?;
        }
        lines.extend(report.atom_mappings.iter().map(|mapping| {
            (
                mapping.location,
                crate::topology::InstanceAtomId::new(instance, mapping.atom),
            )
        }));
        reports.push(report);
    }
    // Dense atom order follows the atom block, even when components interleave.
    lines.sort_unstable_by_key(|&(location, _)| location);
    builder.set_atom_order(lines.into_iter().map(|(_, atom)| atom))?;
    Ok(MolfileInterpretation {
        model: builder.build()?,
        reports,
        omitted,
    })
}

/// Interprets source assertions and tetrahedral configurations carried by 3D
/// coordinates. Coordinate inference uses the shared chemical eligibility and
/// symmetry model; explicit unknown configurations remain unknown. Derived
/// perception is temporary and is not installed on the published molecules.
/// Unsupported valence retains the graph with a diagnostic instead of imposing
/// strict perception on document import. Resource failures remain errors.
pub fn interpret_molfile_document(
    document: &MolfileDocument,
    options: MolfileInterpretOptions,
) -> Result<MolfileInterpretation, MolfileInterpretError> {
    let omitted = match options.unsupported {
        UnsupportedContent::Reject => {
            if let Some(unsupported) = document.unsupported().first() {
                return Err(MolfileInterpretError::at(
                    DiagnosticKind::Unsupported,
                    unsupported.location(),
                    format!("unsupported {}", unsupported.content()),
                ));
            }
            Vec::new()
        }
        UnsupportedContent::Report => document.unsupported().to_vec(),
    };
    let StagedMolfile {
        editor: mut staging,
        geometry,
        mut source_stereo,
    } = stage_molfile(document)?;
    let atom_lines = document
        .atoms()
        .iter()
        .map(MolfileAtom::location)
        .collect::<Vec<_>>();
    let bond_lines = document
        .bonds()
        .iter()
        .map(|bond| bond.location())
        .collect::<Vec<_>>();
    let atom_cfg = staging
        .atom_ids()
        .zip(document.atoms())
        .filter_map(|(id, atom)| atom.parity.map(|parity| (id, parity, atom.location)))
        .collect::<Vec<_>>();
    let tetrahedral_centers = source_stereo
        .iter()
        .filter(|mark| {
            matches!(
                mark.kind,
                SourceStereoBondMarkKind::WedgeUp
                    | SourceStereoBondMarkKind::WedgeDown
                    | SourceStereoBondMarkKind::WedgeEither
            )
        })
        .map(|mark| mark.from)
        .chain(atom_cfg.iter().map(|(atom, _, _)| *atom));
    materialize_molfile_stereo_hydrogens(staging.working_mut(), tetrahedral_centers);
    source_stereo.extend(crate::chemistry::molfile_double_bond_geometry_marks(
        staging.working(),
        &source_stereo,
    ));
    let mut components = Vec::new();
    for mut raw in partition_molfile_staging(staging, &geometry, &source_stereo)? {
        let component_location = raw
            .old_atoms
            .first()
            .and_then(|atom| atom_lines.get(atom.index()))
            .copied();
        let mut editor = raw.editor;
        let mut contradicted_parities = Vec::new();
        let mut publication_report = canonicalize_molecule_for_publication(
            editor.working_mut(),
            Some(&raw.geometry),
            &raw.source_stereo,
        )
        .map_err(|error| {
            MolfileInterpretError::new(
                DiagnosticKind::Chemistry,
                canonicalization_error_location(&error, &atom_lines, &bond_lines),
                format!("could not publish canonical molecule: {error}"),
            )
        })?;
        publication_report
            .created_stereo_elements
            .extend(install_v3000_atom_cfg(
                editor.working_mut(),
                &raw.atom_map,
                &atom_cfg,
                &raw.geometry,
                &mut contradicted_parities,
            )?);
        let mut coordinate_warnings = Vec::new();
        let inferred = molfile_coordinate_tetrahedra(
            editor.working(),
            &raw.geometry,
            &mut coordinate_warnings,
        )?;
        let centers = inferred.iter().filter_map(|element| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => Some(stereo.center),
            _ => None,
        });
        materialize_molfile_stereo_hydrogens(editor.working_mut(), centers);
        for element in inferred {
            publication_report.created_stereo_elements.push(
                editor.add_stereo_element(element).map_err(|error| {
                    MolfileInterpretError::new(
                        DiagnosticKind::Chemistry,
                        component_location,
                        format!("could not materialize Molfile coordinate stereo: {error}"),
                    )
                })?,
            );
        }
        install_molfile_stereo_groups(
            editor.working_mut(),
            &raw.atom_map,
            document.stereo_groups(),
        )?;
        let (molecule, ids) = editor.finish_with_correspondence().map_err(|error| {
            MolfileInterpretError::new(
                DiagnosticKind::Chemistry,
                component_location,
                error.to_string(),
            )
        })?;
        // Publication renumbers IDs densely and may prune canonicalized stereo.
        for target in raw.atom_map.values_mut() {
            *target = ids
                .atom(*target)
                .expect("interpreted atoms survive publication");
        }
        for target in raw.bond_map.values_mut() {
            *target = ids
                .bond(*target)
                .expect("interpreted bonds survive publication");
        }
        publication_report.created_stereo_elements = publication_report
            .created_stereo_elements
            .into_iter()
            .filter_map(|element| ids.stereo_element(element))
            .collect();
        let positions = raw.geometry.to_positions(&molecule).map_err(|error| {
            MolfileInterpretError::new(
                DiagnosticKind::Chemistry,
                component_location,
                format!("could not retain published coordinates: {error}"),
            )
        })?;
        let mut warnings = publication_report
            .warnings
            .into_iter()
            .map(|warning| match warning {
                NormalizationWarning::ConflictingAtropisomericWedgeMarks { axis, mark_count } => {
                    let axis = ids.bond(axis).expect("warned bonds survive publication");
                    let old = raw
                        .bond_map
                        .iter()
                        .find_map(|(old, new)| (*new == axis).then_some(*old))
                        .expect("published bonds map to source bonds");
                    MolfileInterpretationWarning::ConflictingAtropisomericWedgeMarks {
                        axis,
                        location: bond_lines[old.index()],
                        mark_count,
                    }
                }
                NormalizationWarning::AmbiguousTetrahedralWedgeMarks { center, mark_count } => {
                    let center = ids.atom(center).expect("warned atoms survive publication");
                    let old = raw
                        .atom_map
                        .iter()
                        .find_map(|(old, new)| (*new == center).then_some(*old))
                        .expect("published atoms map to source atoms");
                    MolfileInterpretationWarning::AmbiguousTetrahedralWedgeMarks {
                        center,
                        location: atom_lines[old.index()],
                        mark_count,
                    }
                }
            })
            .collect::<Vec<_>>();
        warnings.extend(coordinate_warnings);
        warnings.extend(contradicted_parities.into_iter().map(|(center, location)| {
            MolfileInterpretationWarning::AtomParityContradictsCoordinates {
                center: ids.atom(center).expect("warned atoms survive publication"),
                location,
            }
        }));
        let atom_mappings = raw
            .atom_map
            .iter()
            .map(|(old, atom)| MolfileAtomMapping {
                atom: *atom,
                location: atom_lines[old.index()],
            })
            .collect();
        let bond_mappings = raw
            .bond_map
            .iter()
            .map(|(old, bond)| MolfileBondMapping {
                bond: *bond,
                location: bond_lines[old.index()],
            })
            .collect();
        components.push((
            molecule,
            positions,
            MolfileInterpretationReport {
                atom_mappings,
                bond_mappings,
                created_stereo_elements: publication_report.created_stereo_elements,
                warnings,
            },
        ));
    }
    publish_molfile_components(components, omitted).map_err(|error| {
        MolfileInterpretError::new(
            DiagnosticKind::Chemistry,
            None,
            format!("could not build Molfile model: {error}"),
        )
    })
}

/// Installs V3000 atom parities. A specified parity that contradicts the
/// atom's 3D coordinates is skipped and recorded in `contradicted`, leaving the
/// center to coordinate inference; flat drawings never contradict a parity.
fn install_v3000_atom_cfg(
    molecule: &mut Molecule,
    atom_map: &BTreeMap<AtomId, AtomId>,
    declarations: &[(AtomId, MolfileAtomParity, SourceLocation)],
    geometry: &dyn crate::chemistry::AtomPositionSource,
    contradicted: &mut Vec<(AtomId, SourceLocation)>,
) -> Result<Vec<StereoElementId>, MolfileInterpretError> {
    let mut created = Vec::new();
    for (source, parity, location) in declarations {
        let Some(&center) = atom_map.get(source) else {
            continue;
        };
        let error = |kind, message: String| MolfileInterpretError::at(kind, *location, message);
        // Component remapping preserves atom-block order. CTfile parity uses
        // that order, except hydrogen (explicit or omitted) is numbered last.
        let carriers =
            crate::chemistry::source_tetrahedral_carriers(molecule, center).ok_or_else(|| {
                error(
                    DiagnosticKind::Chemistry,
                    "V3000 atom CFG requires four supported tetrahedral carriers".into(),
                )
            })?;
        let mut orientation = match parity {
            MolfileAtomParity::Odd => Some(TetrahedralOrientation::CounterClockwise),
            MolfileAtomParity::Even => Some(TetrahedralOrientation::Clockwise),
            MolfileAtomParity::Either => None,
        };
        let is_hydrogen = |carrier: &crate::core::StereoCarrier| match carrier {
            crate::core::StereoCarrier::Atom(atom) => molecule
                .atom(*atom)
                .is_ok_and(|atom| atom.element.symbol() == "H"),
            crate::core::StereoCarrier::ImplicitHydrogen => true,
            crate::core::StereoCarrier::ImplicitLonePair => false,
        };
        let hydrogen_swaps = carriers
            .iter()
            .enumerate()
            .filter(|(_, carrier)| is_hydrogen(carrier))
            .map(|(index, _)| {
                carriers[index + 1..]
                    .iter()
                    .filter(|carrier| !is_hydrogen(carrier))
                    .count()
            })
            .sum::<usize>();
        if hydrogen_swaps % 2 != 0 {
            orientation = orientation.map(TetrahedralOrientation::inverted);
        }
        if orientation.is_some_and(|asserted| {
            crate::algorithms::tetrahedral_orientation_at(geometry, center, &carriers)
                .is_some_and(|observed| observed != asserted)
        }) {
            contradicted.push((center, *location));
            continue;
        }
        let existing = molecule
            .stereo_elements()
            .find_map(|(id, element)| match &element.kind {
                StereoElementKind::Tetrahedral(stereo) if stereo.center == center => {
                    Some((id, stereo.clone()))
                }
                _ => None,
            });
        if let Some((id, mut stereo)) = existing {
            if let (Some(existing), Some(asserted)) = (stereo.orientation, orientation) {
                if existing != asserted {
                    return Err(error(
                        DiagnosticKind::Inconsistent,
                        "V3000 atom CFG conflicts with bond wedge configuration".into(),
                    ));
                }
            } else if orientation.is_none() && stereo.orientation.is_some() {
                stereo.orientation = None;
                molecule
                    .replace_stereo_element(
                        id,
                        StereoElement::new(StereoElementKind::Tetrahedral(stereo)),
                    )
                    .map_err(|cause| error(DiagnosticKind::Chemistry, cause.to_string()))?;
            }
        } else {
            let id = molecule
                .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
                    TetrahedralStereo {
                        center,
                        carriers,
                        orientation,
                    },
                )))
                .map_err(|cause| error(DiagnosticKind::Chemistry, cause.to_string()))?;
            created.push(id);
        }
    }
    Ok(created)
}

fn install_molfile_stereo_groups(
    molecule: &mut Molecule,
    atom_map: &BTreeMap<AtomId, AtomId>,
    groups: &[MolfileStereoGroup],
) -> Result<(), MolfileInterpretError> {
    for group in groups {
        let error = |kind, message: &str| MolfileInterpretError::at(kind, group.location, message);
        let kind = match group.kind {
            MolfileStereoGroupKind::Absolute => StereoGroupKind::Absolute,
            MolfileStereoGroupKind::Or(_) => StereoGroupKind::Or,
            MolfileStereoGroupKind::And(_) => StereoGroupKind::And,
        };
        let atoms = group
            .atoms
            .iter()
            .filter_map(|source| {
                atom_map
                    .iter()
                    .find_map(|(old, new)| (old.index() == *source).then_some(*new))
            })
            .collect::<Vec<_>>();
        let present = atoms.len();
        if present == 0 {
            continue;
        }
        if present != group.atoms.len() && kind != StereoGroupKind::Absolute {
            return Err(error(
                DiagnosticKind::Unsupported,
                "relative V3000 stereo groups spanning disconnected molecules cannot be represented",
            ));
        }
        let mut seen_atoms = BTreeSet::new();
        let mut members = BTreeSet::new();
        for atom in atoms {
            if !seen_atoms.insert(atom) {
                return Err(error(
                    DiagnosticKind::Inconsistent,
                    "duplicate V3000 stereo group atom",
                ));
            }
            let candidates = molfile_stereo_group_members_at_atom(molecule, atom);
            let [member] = candidates.as_slice() else {
                return Err(error(
                    DiagnosticKind::Chemistry,
                    "V3000 stereo group atom must identify one tetrahedral center or atropisomeric axis",
                ));
            };
            members.insert(*member);
        }
        molecule
            .add_stereo_group(StereoGroup {
                kind,
                members: members.into_iter().collect(),
            })
            .map_err(|cause| {
                error(
                    DiagnosticKind::Chemistry,
                    &format!("invalid V3000 stereo group: {cause}"),
                )
            })?;
    }
    Ok(())
}

pub(super) fn molfile_stereo_group_members_at_atom(
    molecule: &Molecule,
    atom: AtomId,
) -> Vec<StereoElementId> {
    molecule
        .stereo_elements()
        .filter_map(|(id, element)| {
            let matches = match &element.kind {
                StereoElementKind::Tetrahedral(stereo) => stereo.center == atom,
                StereoElementKind::Axis(stereo) => molecule
                    .bond(stereo.axis)
                    .is_ok_and(|bond| bond.a() == atom || bond.b() == atom),
                StereoElementKind::DoubleBond(_) => false,
            };
            matches.then_some(id)
        })
        .collect()
}

pub(super) fn materialize_molfile_stereo_hydrogens(
    molecule: &mut Molecule,
    centers: impl IntoIterator<Item = AtomId>,
) {
    for center in centers {
        let Ok(atom) = molecule.atom(center) else {
            continue;
        };
        // Molfile's omitted fourth substituent is an implicit hydrogen when
        // its atom valence permits one. Resolve that format convention into a
        // fixed carrier before the represented-only normalization kernel runs.
        // No installed perception or unmarked coordinate inference is involved.
        if !atom.hydrogens.is_inferred()
            || molecule.incident_bonds(center).ok().map(Iterator::count) != Some(3)
            || crate::algorithms::rdkit_inferred_hydrogen_count(molecule, center, atom) != 1
        {
            continue;
        }
        molecule
            .atom_mut(center)
            .expect("source stereo center exists")
            .hydrogens = ImplicitHydrogens::Fixed(1);
    }
}

/// Whether Molfile coordinates form a flat XY drawing rather than a
/// conformation. Import and the writer's header dimension code share it.
pub(super) fn molfile_coordinates_are_flat(points: &[Point3]) -> bool {
    points
        .first()
        .is_none_or(|first| points.iter().all(|point| point.z == first.z))
}

/// Uses the same detached inference for import and for lossless writer planning.
/// A flat XY drawing carries tetrahedral stereo through annotations only.
pub(super) fn molfile_coordinate_tetrahedra(
    molecule: &Molecule,
    geometry: &dyn crate::chemistry::AtomPositionSource,
    warnings: &mut Vec<MolfileInterpretationWarning>,
) -> Result<Vec<StereoElement>, MolfileInterpretError> {
    let points = molecule
        .atom_ids()
        .map(|atom| {
            geometry.position_value(atom).ok_or_else(|| {
                MolfileInterpretError::chemistry(format!("missing source coordinates for {atom}"))
            })
        })
        .collect::<Result<Vec<Point3>, _>>()?;
    if molfile_coordinates_are_flat(&points) {
        return Ok(Vec::new());
    }
    let positions = Positions::from_vec(Quantity::new(points, ANGSTROM))
        .map_err(MolfileInterpretError::chemistry)?;
    let inferred = match crate::algorithms::infer_coordinate_stereo(molecule, &positions) {
        Ok(inferred) => inferred,
        Err(crate::algorithms::CoordinateStereoError::Perception(
            crate::algorithms::StereoPerceptionError::Perception(
                crate::chemistry::PerceptionError::Valence(error),
            ),
        )) => {
            warnings
                .push(MolfileInterpretationWarning::CoordinateStereoValenceUnsupported { error });
            return Ok(Vec::new());
        }
        Err(error) => {
            return Err(MolfileInterpretError::chemistry(format!(
                "could not infer Molfile coordinate stereo: {error}"
            )))
        }
    };
    Ok(inferred
        .elements
        .into_iter()
        .filter(|element| matches!(element.kind, StereoElementKind::Tetrahedral(_)))
        .collect())
}

struct RawMolfileComponent {
    editor: MoleculeEditor,
    geometry: StagedCoordinates,
    source_stereo: Vec<SourceStereoBondMark>,
    old_atoms: Vec<AtomId>,
    atom_map: BTreeMap<AtomId, AtomId>,
    bond_map: BTreeMap<BondId, BondId>,
}

fn partition_molfile_staging(
    staging: MoleculeEditor,
    geometry: &StagedCoordinates,
    source_stereo: &[SourceStereoBondMark],
) -> Result<Vec<RawMolfileComponent>, MolfileInterpretError> {
    let graph = staging.working();
    let mut result = Vec::new();
    for old_atoms in graph.connected_components() {
        let atom_set = old_atoms.iter().copied().collect::<BTreeSet<_>>();
        let mut editor = crate::core::MoleculeEditor::new();
        let mut atom_map = BTreeMap::new();
        let mut component_geometry =
            StagedCoordinates::with_atom_capacity(old_atoms.len(), geometry.unit())
                .map_err(MolfileInterpretError::chemistry)?;
        for old in &old_atoms {
            let atom = graph.atom(*old).map_err(MolfileInterpretError::chemistry)?;
            let new = editor
                .add_atom(atom.clone())
                .map_err(MolfileInterpretError::chemistry)?;
            if let Some(point) = geometry.position(*old) {
                component_geometry
                    .set_position(new, point)
                    .map_err(MolfileInterpretError::chemistry)?;
            }
            atom_map.insert(*old, new);
        }
        let mut bond_map = BTreeMap::new();
        for (old, bond) in graph.bonds() {
            if !atom_set.contains(&bond.a()) || !atom_set.contains(&bond.b()) {
                continue;
            }
            let new = editor
                .add_bond(atom_map[&bond.a()], atom_map[&bond.b()], bond.order)
                .map_err(MolfileInterpretError::chemistry)?;
            bond_map.insert(old, new);
        }
        let mut properties =
            crate::properties::MoleculeProperties::new(atom_map.len(), bond_map.len());
        *properties.atoms_mut() = graph
            .properties()
            .atoms()
            .raw()
            .select_indices(&old_atoms.iter().map(|id| id.index()).collect::<Vec<_>>())
            .map_err(MolfileInterpretError::chemistry)?;
        *properties.bonds_mut() = graph
            .properties()
            .bonds()
            .raw()
            .select_indices(&bond_map.keys().map(|id| id.index()).collect::<Vec<_>>())
            .map_err(MolfileInterpretError::chemistry)?;
        editor.working_mut().properties = properties;
        let remapped_stereo = source_stereo
            .iter()
            .filter_map(|mark| {
                Some(SourceStereoBondMark {
                    bond: *bond_map.get(&mark.bond)?,
                    from: *atom_map.get(&mark.from)?,
                    kind: mark.kind,
                })
            })
            .collect();
        result.push(RawMolfileComponent {
            editor,
            geometry: component_geometry,
            source_stereo: remapped_stereo,
            old_atoms,
            atom_map,
            bond_map,
        });
    }
    Ok(result)
}

fn canonicalization_error_location(
    error: &NormalizationError,
    atom_lines: &[SourceLocation],
    bond_lines: &[SourceLocation],
) -> Option<SourceLocation> {
    error
        .bond_location_hint()
        .and_then(|bond| bond_lines.get(bond.index()).copied())
        .or_else(|| {
            error
                .atom_location_hint()
                .and_then(|atom| atom_lines.get(atom.index()).copied())
        })
}

#[cfg(test)]
mod stereo_group_tests {
    use super::*;
    use crate::core::{AxisOrientation, AxisStereo, StereoCarrier};

    #[test]
    fn enhanced_group_rejects_an_endpoint_shared_by_two_axes() {
        // This represented graph tests source-member identity, independently
        // of the chemical stability or coordinate perception of either axis.
        let mut editor = MoleculeEditor::new();
        let atoms = (0..5)
            .map(|_| {
                editor
                    .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let bonds = atoms
            .windows(2)
            .map(|pair| {
                editor
                    .add_bond(pair[0], pair[1], BondOrder::Single)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let mut axes = Vec::new();
        for (bond, left, right) in [
            (bonds[1], atoms[0], atoms[3]),
            (bonds[2], atoms[1], atoms[4]),
        ] {
            axes.push(
                editor
                    .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
                        axis: bond,
                        carriers: vec![StereoCarrier::Atom(left), StereoCarrier::Atom(right)],
                        orientation: Some(AxisOrientation::Clockwise),
                    })))
                    .unwrap(),
            );
        }
        let mut molecule = editor.finish().unwrap();
        assert_eq!(
            molfile_stereo_group_members_at_atom(&molecule, atoms[2]),
            axes
        );
        assert_eq!(
            molfile_stereo_group_members_at_atom(&molecule, atoms[1]),
            [axes[0]]
        );
        assert_eq!(
            molfile_stereo_group_members_at_atom(&molecule, atoms[3]),
            [axes[1]]
        );
        let atom_map = atoms.iter().copied().map(|atom| (atom, atom)).collect();
        let mut group = MolfileStereoGroup {
            kind: MolfileStereoGroupKind::And(1),
            atoms: vec![atoms[2].index()],
            location: SourceLocation::at_line(73),
        };
        let error =
            install_molfile_stereo_groups(&mut molecule, &atom_map, &[group.clone()]).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Chemistry);
        assert_eq!(error.location(), Some(SourceLocation::at_line(73)));
        assert_eq!(
            error.message(),
            "V3000 stereo group atom must identify one tetrahedral center or atropisomeric axis"
        );
        assert_eq!(molecule.stereo_groups().count(), 0);

        // The other endpoint of each axis identifies that member uniquely.
        group.atoms = vec![atoms[1].index(), atoms[3].index()];
        install_molfile_stereo_groups(&mut molecule, &atom_map, &[group]).unwrap();
        assert_eq!(molecule.stereo_groups().next().unwrap().1.members, axes);
    }
}
