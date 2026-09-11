use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::chemistry::{
    canonicalize_molecule_for_publication, NormalizationError, NormalizationWarning,
};
use crate::core::{
    Atom, AtomId, AtomRadical, BondId, BondOrder, Element, HydrogenDeclaration, Molecule,
    MoleculeEditor, StereoElement, StereoElementId, StereoElementKind, StereoGroup,
    StereoGroupKind, TetrahedralOrientation, TetrahedralStereo,
};
use crate::structure::{Model, ModelBuildError, ModelBuilder, Positions};
use crate::topology::Topology;

use super::staged_coordinates::StagedCoordinates;
use super::v2000::{interpret_v2000_syntax, parse_counts_line, parse_v2000_syntax, V2000Syntax};
use super::v3000::{interpret_v3000_syntax, parse_v3000_syntax, V3000Syntax};
use super::SdfParseError;

pub(super) fn checked_line_number(
    record: usize,
    start_line: usize,
    offset: usize,
) -> Result<usize, SdfParseError> {
    start_line
        .checked_add(offset)
        .ok_or_else(|| SdfParseError::new(record, start_line, "line number overflow"))
}

pub(super) fn interpret_molfile_atom_fields(
    symbol: &str,
    formal_charge: i32,
    isotope: Option<i32>,
    radical: Option<AtomRadical>,
    hydrogens: HydrogenDeclaration,
    atom_map: Option<u32>,
    line: usize,
) -> Result<Atom, SdfParseError> {
    let element = Element::from_symbol(symbol).ok_or_else(|| {
        SdfParseError::new(1, line, format!("unsupported element symbol `{symbol}`"))
    })?;
    let mut atom = Atom::new(element);
    atom.formal_charge = i8::try_from(formal_charge)
        .map_err(|_| SdfParseError::new(1, line, "formal charge is outside i8 range"))?;
    atom.isotope = match isotope {
        Some(value) if value > 0 => Some(
            u16::try_from(value)
                .map_err(|_| SdfParseError::new(1, line, "isotope is outside u16 range"))?,
        ),
        _ => None,
    };
    atom.radical = radical;
    atom.hydrogens = hydrogens;
    atom.atom_map = atom_map;
    Ok(atom)
}

pub(super) fn apply_molfile_declared_valence(
    molecule: &mut Molecule,
    atom: AtomId,
    declared_valence: u8,
    declared_hydrogens: Option<u8>,
    source_aromatic_bonds: &BTreeSet<BondId>,
    version: MolfileVersion,
    line: usize,
) -> Result<(), SdfParseError> {
    let version_name = match version {
        MolfileVersion::V2000 => "V2000",
        MolfileVersion::V3000 => "V3000",
    };
    let unsupported_order = || {
        SdfParseError::new(
            1,
            line,
            format!(
                "declared {version_name} valence cannot determine hydrogens for this bond order"
            ),
        )
    };
    let represented_valence = molecule
        .incident_bonds(atom)
        .map_err(|error| SdfParseError::new(1, line, error.to_string()))?
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
            SdfParseError::new(
                1,
                line,
                format!("declared {version_name} valence is below represented bond valence"),
            )
        })?;
    if let Some(declared_hydrogens) = declared_hydrogens {
        if inferred_hydrogens != usize::from(declared_hydrogens) {
            let message = match version {
                MolfileVersion::V2000 => "V2000 hydrogen-count and valence declarations conflict",
                MolfileVersion::V3000 => "V3000 HCOUNT and VAL declarations conflict",
            };
            return Err(SdfParseError::new(1, line, message));
        }
        return Ok(());
    }

    let explicit = u8::try_from(inferred_hydrogens).map_err(|_| {
        SdfParseError::new(
            1,
            line,
            format!("declared {version_name} hydrogen count exceeds u8"),
        )
    })?;
    molecule
        .atom_mut(atom)
        .expect("interpreted Molfile atom remains live")
        .hydrogens = HydrogenDeclaration::Fixed(explicit);
    Ok(())
}

/// Resource limits shared by version-autodetected Molfile parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolfileParseOptions {
    /// Maximum UTF-8 byte length of the complete Molfile document.
    pub max_input_bytes: usize,
    /// Maximum declared V3000 atom count.
    pub max_v3000_atoms: usize,
    /// Maximum declared V3000 bond count.
    pub max_v3000_bonds: usize,
    /// Maximum byte length after joining one continued `M  V30` logical line.
    pub max_v3000_logical_line_bytes: usize,
}

impl Default for MolfileParseOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 * 1024 * 1024,
            max_v3000_atoms: 1_000_000,
            max_v3000_bonds: 2_000_000,
            max_v3000_logical_line_bytes: 1024 * 1024,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MolfileVersion {
    V2000,
    V3000,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileHeader {
    title: String,
    program: String,
    comment: String,
}

impl MolfileHeader {
    pub fn title(&self) -> &str {
        &self.title
    }

    pub fn program(&self) -> &str {
        &self.program
    }

    pub fn comment(&self) -> &str {
        &self.comment
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileLine {
    number: usize,
    text: String,
}

impl MolfileLine {
    pub const fn number(&self) -> usize {
        self.number
    }

    pub fn text(&self) -> &str {
        &self.text
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MolfileDocument {
    source: String,
    version: MolfileVersion,
    header: MolfileHeader,
    atom_records: Vec<MolfileLine>,
    bond_records: Vec<MolfileLine>,
    property_records: Vec<MolfileLine>,
    unsupported_records: Vec<MolfileLine>,
    syntax: MolfileSyntax,
}

#[derive(Debug, Clone, PartialEq)]
enum MolfileSyntax {
    V2000(V2000Syntax),
    V3000(V3000Syntax),
}

impl MolfileDocument {
    pub fn source(&self) -> &str {
        &self.source
    }

    pub const fn version(&self) -> MolfileVersion {
        self.version
    }

    pub fn header(&self) -> &MolfileHeader {
        &self.header
    }

    pub fn atom_records(&self) -> &[MolfileLine] {
        &self.atom_records
    }

    pub fn bond_records(&self) -> &[MolfileLine] {
        &self.bond_records
    }

    pub fn property_records(&self) -> &[MolfileLine] {
        &self.property_records
    }

    pub fn unsupported_records(&self) -> &[MolfileLine] {
        &self.unsupported_records
    }

    /// Interprets this document once into connected chemistry with matching
    /// component geometry and source reports.
    pub fn interpret(&self) -> Result<MolfileInterpretation, MolfileInterpretError> {
        interpret_molfile_document(self)
    }

    /// Interprets this source document as connected canonical molecules.
    ///
    /// Coordinates may participate in source-stereo normalization, but are not
    /// retained in the returned molecules. No chemical perception is run.
    pub fn to_molecules(&self) -> Result<Vec<Molecule>, MolfileInterpretError> {
        Ok(self.interpret()?.into_molecules())
    }

    /// Interprets this source document and projects its complete static model
    /// layout, including the deterministic synthetic hierarchy.
    pub fn to_topology(&self) -> Result<std::sync::Arc<Topology>, MolfileInterpretError> {
        Ok(self.interpret()?.into_topology())
    }

    /// Interprets this source document as one geometry-bearing model.
    ///
    /// Each disconnected source component becomes one molecule instance in
    /// the model topology. No chemical perception is run.
    pub fn to_model(&self) -> Result<Model, MolfileInterpretError> {
        Ok(self.interpret()?.into_model())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileParseError {
    pub(crate) line: usize,
    pub(crate) message: String,
}

impl MolfileParseError {
    fn new(line: usize, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }

    pub const fn line(&self) -> usize {
        self.line
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for MolfileParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Molfile parse error at line {}: {}",
            self.line, self.message
        )
    }
}

impl std::error::Error for MolfileParseError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileInterpretError {
    pub(crate) line: usize,
    pub(crate) message: String,
}

impl MolfileInterpretError {
    pub const fn line(&self) -> usize {
        self.line
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for MolfileInterpretError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Molfile interpretation error at line {}: {}",
            self.line, self.message
        )
    }
}

impl std::error::Error for MolfileInterpretError {}

impl From<SdfParseError> for MolfileInterpretError {
    fn from(error: SdfParseError) -> Self {
        Self {
            line: error.line,
            message: error.message,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolfileAtomMapping {
    atom: crate::core::AtomId,
    source_line: usize,
}

impl MolfileAtomMapping {
    pub const fn atom(&self) -> crate::core::AtomId {
        self.atom
    }

    pub const fn source_line(&self) -> usize {
        self.source_line
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolfileBondMapping {
    bond: crate::core::BondId,
    source_line: usize,
}

impl MolfileBondMapping {
    pub const fn bond(&self) -> crate::core::BondId {
        self.bond
    }

    pub const fn source_line(&self) -> usize {
        self.source_line
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MolfileInterpretationReport {
    atom_mappings: Vec<MolfileAtomMapping>,
    bond_mappings: Vec<MolfileBondMapping>,
    ignored_record_lines: Vec<usize>,
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

    pub fn ignored_record_lines(&self) -> &[usize] {
        &self.ignored_record_lines
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
pub enum MolfileInterpretationWarning {
    AmbiguousTetrahedralWedgeMarks {
        center: AtomId,
        source_line: usize,
        mark_count: usize,
    },
}

/// One final canonical model and source reports in molecule-instance order.
#[derive(Debug, Clone, PartialEq)]
pub struct MolfileInterpretation {
    model: Model,
    reports: Vec<MolfileInterpretationReport>,
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

    pub fn into_parts(self) -> (Model, Vec<MolfileInterpretationReport>) {
        (self.model, self.reports)
    }
}

fn publish_molfile_components(
    components: Vec<(Molecule, Positions, MolfileInterpretationReport)>,
) -> Result<MolfileInterpretation, ModelBuildError> {
    let mut builder = ModelBuilder::new();
    let mut reports = Vec::with_capacity(components.len());
    let chain = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_chain("A", None)?;
    for (index, (molecule, positions, report)) in components.into_iter().enumerate() {
        let instance = builder.add_molecule(&molecule, &positions)?;
        let sequence = i32::try_from(index + 1).map_err(|_| ModelBuildError::CapacityOverflow)?;
        let hierarchy = builder.topology_builder_mut().hierarchy_mut();
        let residue = hierarchy.add_residue(
            chain,
            "UNL",
            Some(sequence),
            Some(sequence.to_string()),
            None,
        )?;
        for atom in molecule.atom_ids() {
            hierarchy.add_atom_site(
                residue,
                crate::topology::InstanceAtomId::new(instance, atom),
                crate::topology::AtomSiteMetadata::default(),
            )?;
        }
        reports.push(report);
    }
    Ok(MolfileInterpretation {
        model: builder.build()?,
        reports,
    })
}
pub fn parse_molfile_document(input: &str) -> Result<MolfileDocument, MolfileParseError> {
    parse_molfile_document_with_options(input, MolfileParseOptions::default())
}

pub fn parse_molfile_document_with_options(
    input: &str,
    options: MolfileParseOptions,
) -> Result<MolfileDocument, MolfileParseError> {
    if input.len() > options.max_input_bytes {
        return Err(MolfileParseError::new(
            1,
            "input exceeds configured byte limit",
        ));
    }
    let source = input.replace("\r\n", "\n").replace('\r', "\n");
    let lines = source.lines().collect::<Vec<_>>();
    if lines.len() < 4 {
        return Err(MolfileParseError::new(
            1,
            "document must contain three header lines and a counts line",
        ));
    }
    let version = if lines[3].contains("V3000") {
        MolfileVersion::V3000
    } else if lines[3].contains("V2000") {
        MolfileVersion::V2000
    } else {
        return Err(MolfileParseError::new(
            4,
            "counts line must declare V2000 or V3000",
        ));
    };
    let end = lines
        .iter()
        .enumerate()
        .skip(4)
        .find_map(|(index, line)| (line.trim() == "M  END").then_some(index))
        .ok_or_else(|| MolfileParseError::new(lines.len(), "missing M  END"))?;
    let header = MolfileHeader {
        title: lines[0].to_owned(),
        program: lines[1].to_owned(),
        comment: lines[2].to_owned(),
    };
    let mut atom_records = Vec::new();
    let mut bond_records = Vec::new();
    let mut property_records = Vec::new();
    let mut unsupported_records = Vec::new();
    match version {
        MolfileVersion::V2000 => {
            let (atom_count, bond_count) = parse_counts_line(lines[3])
                .ok_or_else(|| MolfileParseError::new(4, "invalid V2000 counts line"))?;
            let atom_end = 4usize
                .checked_add(atom_count)
                .ok_or_else(|| MolfileParseError::new(4, "atom count overflow"))?;
            let bond_end = atom_end
                .checked_add(bond_count)
                .ok_or_else(|| MolfileParseError::new(4, "bond count overflow"))?;
            if bond_end > end {
                return Err(MolfileParseError::new(
                    4,
                    "counts exceed records before M  END",
                ));
            }
            atom_records.extend(lines[4..atom_end].iter().enumerate().map(|(offset, line)| {
                MolfileLine {
                    number: offset + 5,
                    text: (*line).to_owned(),
                }
            }));
            bond_records.extend(lines[atom_end..bond_end].iter().enumerate().map(
                |(offset, line)| MolfileLine {
                    number: atom_end + offset + 1,
                    text: (*line).to_owned(),
                },
            ));
            for (offset, line) in lines[bond_end..end].iter().enumerate() {
                let record = MolfileLine {
                    number: bond_end + offset + 1,
                    text: (*line).to_owned(),
                };
                if line.starts_with("M  ") {
                    property_records.push(record);
                } else if !line.trim().is_empty() {
                    unsupported_records.push(record);
                }
            }
        }
        MolfileVersion::V3000 => {
            let mut section = None::<&str>;
            for (offset, line) in lines[4..end].iter().enumerate() {
                let number = offset + 5;
                let body = line.strip_prefix("M  V30 ").map(str::trim);
                match body {
                    Some("BEGIN CTAB" | "END CTAB") => {}
                    Some(body) if body.starts_with("COUNTS ") => {}
                    Some("BEGIN ATOM") => section = Some("ATOM"),
                    Some("END ATOM") => section = None,
                    Some("BEGIN BOND") => section = Some("BOND"),
                    Some("END BOND") => section = None,
                    Some(_) => {
                        let record = MolfileLine {
                            number,
                            text: (*line).to_owned(),
                        };
                        match section {
                            Some("ATOM") => atom_records.push(record),
                            Some("BOND") => bond_records.push(record),
                            _ => property_records.push(record),
                        }
                    }
                    None if !line.trim().is_empty() => unsupported_records.push(MolfileLine {
                        number,
                        text: (*line).to_owned(),
                    }),
                    None => {}
                }
            }
        }
    }
    unsupported_records.extend(
        lines[end + 1..]
            .iter()
            .enumerate()
            .filter(|(_, line)| !line.trim().is_empty())
            .map(|(offset, line)| MolfileLine {
                number: end + offset + 2,
                text: (*line).to_owned(),
            }),
    );
    let syntax = match version {
        MolfileVersion::V2000 => MolfileSyntax::V2000(
            parse_v2000_syntax(1, 1, &lines[..=end])
                .map_err(|error| MolfileParseError::new(error.line, error.message))?,
        ),
        MolfileVersion::V3000 => MolfileSyntax::V3000(
            parse_v3000_syntax(1, 1, &lines[..=end], options)
                .map_err(|error| MolfileParseError::new(error.line, error.message))?,
        ),
    };
    Ok(MolfileDocument {
        source,
        version,
        header,
        atom_records,
        bond_records,
        property_records,
        unsupported_records,
        syntax,
    })
}

pub fn interpret_molfile_document(
    document: &MolfileDocument,
) -> Result<MolfileInterpretation, MolfileInterpretError> {
    let ((mut staging, geometry, mut source_stereo), atom_lines, bond_lines) =
        match &document.syntax {
            MolfileSyntax::V2000(syntax) => (
                interpret_v2000_syntax(syntax)?,
                syntax
                    .atoms
                    .iter()
                    .map(|record| record.line)
                    .collect::<Vec<_>>(),
                syntax
                    .bonds
                    .iter()
                    .map(|record| record.line)
                    .collect::<Vec<_>>(),
            ),
            MolfileSyntax::V3000(syntax) => (
                interpret_v3000_syntax(syntax)?,
                syntax
                    .atoms
                    .iter()
                    .map(|record| record.line)
                    .collect::<Vec<_>>(),
                syntax
                    .bonds
                    .iter()
                    .map(|record| record.line)
                    .collect::<Vec<_>>(),
            ),
        };
    let stereo_groups = match &document.syntax {
        MolfileSyntax::V2000(_) => &[][..],
        MolfileSyntax::V3000(syntax) => syntax.stereo_groups.as_slice(),
    };
    let atom_cfg = match &document.syntax {
        MolfileSyntax::V2000(_) => Vec::new(),
        MolfileSyntax::V3000(syntax) => staging
            .atom_ids()
            .zip(&syntax.atoms)
            .filter_map(|(id, atom)| atom.stereo_cfg.map(|cfg| (id, cfg, atom.line)))
            .collect(),
    };
    let tetrahedral_centers = source_stereo
        .iter()
        .filter(|mark| {
            matches!(
                mark.kind,
                crate::chemistry::SourceStereoBondMarkKind::WedgeUp
                    | crate::chemistry::SourceStereoBondMarkKind::WedgeDown
                    | crate::chemistry::SourceStereoBondMarkKind::WedgeEither
            )
        })
        .map(|mark| mark.from)
        .chain(atom_cfg.iter().map(|(atom, _, _)| *atom));
    materialize_molfile_stereo_hydrogens(staging.working_mut(), tetrahedral_centers);
    source_stereo.extend(crate::chemistry::molfile_double_bond_geometry_marks(
        staging.working(),
        &source_stereo,
    ));
    let ignored_record_lines: Vec<usize> = document
        .property_records
        .iter()
        .filter(|record| {
            if stereo_groups
                .iter()
                .any(|group| group.source_lines.contains(&record.number))
                || (!stereo_groups.is_empty()
                    && matches!(
                        record.text.trim(),
                        "M  V30 BEGIN COLLECTION" | "M  V30 END COLLECTION"
                    ))
            {
                return false;
            }
            let mut fields = record.text.split_whitespace();
            !matches!(
                (fields.next(), fields.next()),
                (Some("M"), Some("CHG" | "ISO" | "RAD"))
            )
        })
        .chain(document.unsupported_records.iter())
        .map(|record| record.number)
        .collect();
    let mut components = Vec::new();
    for raw in partition_molfile_staging(staging, &geometry, &source_stereo)? {
        let mut editor = raw.editor;
        let mut publication_report = canonicalize_molecule_for_publication(
            editor.working_mut(),
            Some(&raw.geometry),
            &raw.source_stereo,
        )
        .map_err(|error| MolfileInterpretError {
            line: canonicalization_error_line(&error, &atom_lines, &bond_lines),
            message: format!("could not publish canonical molecule: {error}"),
        })?;
        publication_report
            .created_stereo_elements
            .extend(install_v3000_atom_cfg(
                editor.working_mut(),
                &raw.atom_map,
                &atom_cfg,
            )?);
        install_molfile_stereo_groups(editor.working_mut(), &raw.atom_map, stereo_groups)?;
        let molecule = editor.finish().map_err(|error| MolfileInterpretError {
            line: raw
                .old_atoms
                .first()
                .and_then(|atom| atom_lines.get(atom.index()))
                .copied()
                .unwrap_or(1),
            message: error.to_string(),
        })?;
        let positions =
            raw.geometry
                .to_positions(&molecule)
                .map_err(|error| MolfileInterpretError {
                    line: raw
                        .old_atoms
                        .first()
                        .and_then(|atom| atom_lines.get(atom.index()))
                        .copied()
                        .unwrap_or(1),
                    message: format!("could not retain published coordinates: {error}"),
                })?;
        let warnings = publication_report
            .warnings
            .into_iter()
            .map(|warning| match warning {
                NormalizationWarning::AmbiguousTetrahedralWedgeMarks { center, mark_count } => {
                    let source_line = raw
                        .atom_map
                        .iter()
                        .find_map(|(old, new)| (*new == center).then(|| atom_lines[old.index()]))
                        .unwrap_or(1);
                    MolfileInterpretationWarning::AmbiguousTetrahedralWedgeMarks {
                        center,
                        source_line,
                        mark_count,
                    }
                }
            })
            .collect();
        let atom_mappings = raw
            .atom_map
            .iter()
            .map(|(old, atom)| MolfileAtomMapping {
                atom: *atom,
                source_line: atom_lines[old.index()],
            })
            .collect();
        let bond_mappings = raw
            .bond_map
            .iter()
            .map(|(old, bond)| MolfileBondMapping {
                bond: *bond,
                source_line: bond_lines[old.index()],
            })
            .collect();
        components.push((
            molecule,
            positions,
            MolfileInterpretationReport {
                atom_mappings,
                bond_mappings,
                ignored_record_lines: ignored_record_lines.clone(),
                created_stereo_elements: publication_report.created_stereo_elements,
                warnings,
            },
        ));
    }
    publish_molfile_components(components).map_err(|error| MolfileInterpretError {
        line: 4,
        message: format!("could not build Molfile model: {error}"),
    })
}

fn install_v3000_atom_cfg(
    molecule: &mut Molecule,
    atom_map: &BTreeMap<AtomId, AtomId>,
    declarations: &[(AtomId, u8, usize)],
) -> Result<Vec<StereoElementId>, MolfileInterpretError> {
    let mut created = Vec::new();
    for (source, cfg, line) in declarations {
        let Some(&center) = atom_map.get(source) else {
            continue;
        };
        let error = |message: String| MolfileInterpretError {
            line: *line,
            message,
        };
        // Component remapping preserves atom-block order. CTfile parity uses
        // that order, except hydrogen (explicit or omitted) is numbered last.
        let carriers =
            crate::chemistry::source_tetrahedral_carriers(molecule, center).ok_or_else(|| {
                error("V3000 atom CFG requires four supported tetrahedral carriers".into())
            })?;
        let mut orientation = match cfg {
            1 => Some(TetrahedralOrientation::CounterClockwise),
            2 => Some(TetrahedralOrientation::Clockwise),
            3 => None,
            _ => unreachable!("parser validates atom CFG"),
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
                    .map_err(|cause| error(cause.to_string()))?;
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
                .map_err(|cause| error(cause.to_string()))?;
            created.push(id);
        }
    }
    Ok(created)
}

fn install_molfile_stereo_groups(
    molecule: &mut Molecule,
    atom_map: &BTreeMap<AtomId, AtomId>,
    groups: &[super::v3000::V3000StereoGroupSyntax],
) -> Result<(), MolfileInterpretError> {
    for group in groups {
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
        if present != group.atoms.len() && group.kind != StereoGroupKind::Absolute {
            return Err(MolfileInterpretError { line: group.line, message: "relative V3000 stereo groups spanning disconnected molecules cannot be represented".to_owned() });
        }
        let mut seen_atoms = BTreeSet::new();
        let mut members = BTreeSet::new();
        for atom in atoms {
            if !seen_atoms.insert(atom) {
                return Err(MolfileInterpretError {
                    line: group.line,
                    message: "duplicate V3000 stereo group atom".to_owned(),
                });
            }
            let candidates = molfile_stereo_group_members_at_atom(molecule, atom);
            let [member] = candidates.as_slice() else {
                return Err(MolfileInterpretError {
                line: group.line,
                message: "V3000 stereo group atom must identify one tetrahedral center or atropisomeric axis".to_owned(),
                });
            };
            members.insert(*member);
        }
        molecule
            .add_stereo_group(StereoGroup {
                kind: group.kind,
                members: members.into_iter().collect(),
            })
            .map_err(|error| MolfileInterpretError {
                line: group.line,
                message: format!("invalid V3000 stereo group: {error}"),
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

fn materialize_molfile_stereo_hydrogens(
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
        if !atom.hydrogens.allows_implicit()
            || molecule.incident_bonds(center).ok().map(Iterator::count) != Some(3)
            || crate::algorithms::rdkit_implicit_hydrogen_count(molecule, center, atom) != 1
        {
            continue;
        }
        molecule
            .atom_mut(center)
            .expect("source stereo center exists")
            .hydrogens = HydrogenDeclaration::Fixed(1);
    }
}

struct RawMolfileComponent {
    editor: MoleculeEditor,
    geometry: StagedCoordinates,
    source_stereo: Vec<crate::chemistry::SourceStereoBondMark>,
    old_atoms: Vec<AtomId>,
    atom_map: BTreeMap<AtomId, AtomId>,
    bond_map: BTreeMap<BondId, BondId>,
}

fn partition_molfile_staging(
    staging: MoleculeEditor,
    geometry: &StagedCoordinates,
    source_stereo: &[crate::chemistry::SourceStereoBondMark],
) -> Result<Vec<RawMolfileComponent>, MolfileInterpretError> {
    let graph = staging.working();
    let mut result = Vec::new();
    for old_atoms in graph.connected_components() {
        let atom_set = old_atoms.iter().copied().collect::<BTreeSet<_>>();
        let mut editor = crate::core::MoleculeEditor::new();
        let mut atom_map = BTreeMap::new();
        let mut component_geometry =
            StagedCoordinates::with_atom_capacity(old_atoms.len(), geometry.unit()).map_err(
                |error| MolfileInterpretError {
                    line: 1,
                    message: error.to_string(),
                },
            )?;
        for old in &old_atoms {
            let atom = graph.atom(*old).map_err(|error| MolfileInterpretError {
                line: 1,
                message: error.to_string(),
            })?;
            let new = editor
                .add_atom(atom.clone())
                .map_err(|error| MolfileInterpretError {
                    line: 1,
                    message: error.to_string(),
                })?;
            if let Some(point) = geometry.position(*old) {
                component_geometry
                    .set_position(new, point)
                    .map_err(|error| MolfileInterpretError {
                        line: 1,
                        message: error.to_string(),
                    })?;
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
                .map_err(|error| MolfileInterpretError {
                    line: 1,
                    message: error.to_string(),
                })?;
            bond_map.insert(old, new);
        }
        let mut properties =
            crate::properties::Properties::molecule(atom_map.len(), bond_map.len());
        *properties.atoms_mut() = graph
            .properties()
            .atoms()
            .select_indices(&old_atoms.iter().map(|id| id.index()).collect::<Vec<_>>())
            .map_err(|error| MolfileInterpretError {
                line: 1,
                message: error.to_string(),
            })?;
        *properties.bonds_mut() = graph
            .properties()
            .bonds()
            .select_indices(&bond_map.keys().map(|id| id.index()).collect::<Vec<_>>())
            .map_err(|error| MolfileInterpretError {
                line: 1,
                message: error.to_string(),
            })?;
        *editor.working_mut().properties_mut() = properties;
        let remapped_stereo = source_stereo
            .iter()
            .filter_map(|mark| {
                Some(crate::chemistry::SourceStereoBondMark {
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

fn canonicalization_error_line(
    error: &NormalizationError,
    atom_lines: &[usize],
    bond_lines: &[usize],
) -> usize {
    error
        .bond_location_hint()
        .and_then(|bond| bond_lines.get(bond.index()).copied())
        .or_else(|| {
            error
                .atom_location_hint()
                .and_then(|atom| atom_lines.get(atom.index()).copied())
        })
        .or_else(|| bond_lines.first().copied())
        .or_else(|| atom_lines.first().copied())
        .unwrap_or(1)
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
        let mut group = super::super::v3000::V3000StereoGroupSyntax {
            kind: StereoGroupKind::And,
            atoms: vec![atoms[2].index()],
            line: 73,
            source_lines: vec![73],
        };
        let error =
            install_molfile_stereo_groups(&mut molecule, &atom_map, &[group.clone()]).unwrap_err();
        assert_eq!(error.line(), 73);
        assert!(error
            .message()
            .contains("must identify one tetrahedral center or atropisomeric axis"));
        assert_eq!(molecule.stereo_groups().count(), 0);

        // The other endpoint of each axis identifies that member uniquely.
        group.atoms = vec![atoms[1].index(), atoms[3].index()];
        install_molfile_stereo_groups(&mut molecule, &atom_map, &[group]).unwrap();
        assert_eq!(molecule.stereo_groups().next().unwrap().1.members, axes);
    }
}
