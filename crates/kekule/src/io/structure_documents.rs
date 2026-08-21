use std::collections::BTreeSet;
use std::fmt;

use crate::chemistry::{
    canonicalize_molecule_for_publication, NormalizationError, NormalizationWarning,
};
use crate::core::{
    Atom, AtomId, AtomRadical, BondId, BondOrder, Element, HydrogenDeclaration, Molecule,
    StereoElementId,
};
use crate::small::model::SmallMolecule;

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

#[derive(Debug, Clone, PartialEq)]
pub struct MolfileInterpretation {
    molecule: SmallMolecule,
    report: MolfileInterpretationReport,
}

impl MolfileInterpretation {
    pub fn molecule(&self) -> &SmallMolecule {
        &self.molecule
    }

    pub fn report(&self) -> &MolfileInterpretationReport {
        &self.report
    }

    pub fn to_molecule(self) -> SmallMolecule {
        self.molecule
    }

    pub fn to_parts(self) -> (SmallMolecule, MolfileInterpretationReport) {
        (self.molecule, self.report)
    }
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
        .position(|line| line.trim() == "M  END")
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
    let ((mut molecule, source_stereo), atom_lines, bond_lines) = match &document.syntax {
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
    let publication_report =
        canonicalize_molecule_for_publication(molecule.as_molecule_mut(), &source_stereo).map_err(
            |error| MolfileInterpretError {
                line: canonicalization_error_line(&error, &atom_lines, &bond_lines),
                message: format!("could not publish canonical molecule: {error}"),
            },
        )?;
    let warnings = publication_report
        .warnings
        .into_iter()
        .map(|warning| match warning {
            NormalizationWarning::AmbiguousTetrahedralWedgeMarks { center, mark_count } => {
                MolfileInterpretationWarning::AmbiguousTetrahedralWedgeMarks {
                    center,
                    source_line: atom_lines.get(center.index()).copied().unwrap_or(1),
                    mark_count,
                }
            }
        })
        .collect();
    let atom_mappings = atom_lines
        .into_iter()
        .zip(molecule.as_molecule().atom_ids())
        .map(|(source_line, atom)| MolfileAtomMapping { atom, source_line })
        .collect();
    let bond_mappings = bond_lines
        .into_iter()
        .zip(molecule.as_molecule().bond_ids())
        .map(|(source_line, bond)| MolfileBondMapping { bond, source_line })
        .collect();
    let ignored_record_lines = document
        .property_records
        .iter()
        .filter(|record| {
            let mut fields = record.text.split_whitespace();
            !matches!(
                (fields.next(), fields.next()),
                (Some("M"), Some("CHG" | "ISO" | "RAD"))
            )
        })
        .chain(document.unsupported_records.iter())
        .map(|record| record.number)
        .collect();
    Ok(MolfileInterpretation {
        molecule,
        report: MolfileInterpretationReport {
            atom_mappings,
            bond_mappings,
            ignored_record_lines,
            created_stereo_elements: publication_report.created_stereo_elements,
            warnings,
        },
    })
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
