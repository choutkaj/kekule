use std::fmt;

use crate::small::model::SmallMolecule;

use super::v2000::{interpret_v2000_syntax, parse_counts_line, parse_v2000_syntax, V2000Syntax};
use super::v3000::{interpret_v3000_syntax, parse_v3000_syntax, V3000Syntax};
use super::SdfParseError;

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

    pub fn into_molecule(self) -> SmallMolecule {
        self.molecule
    }

    pub fn into_parts(self) -> (SmallMolecule, MolfileInterpretationReport) {
        (self.molecule, self.report)
    }
}

pub fn parse_molfile_document(input: &str) -> Result<MolfileDocument, MolfileParseError> {
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
    let syntax = match version {
        MolfileVersion::V2000 => MolfileSyntax::V2000(
            parse_v2000_syntax(1, 1, &lines[..=end])
                .map_err(|error| MolfileParseError::new(error.line, error.message))?,
        ),
        MolfileVersion::V3000 => MolfileSyntax::V3000(
            parse_v3000_syntax(1, 1, &lines[..=end])
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
    let (molecule, atom_lines, bond_lines) = match &document.syntax {
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
    let atom_mappings = atom_lines
        .into_iter()
        .enumerate()
        .map(|(index, source_line)| MolfileAtomMapping {
            atom: crate::core::AtomId::new(index as u32),
            source_line,
        })
        .collect();
    let bond_mappings = bond_lines
        .into_iter()
        .enumerate()
        .map(|(index, source_line)| MolfileBondMapping {
            bond: crate::core::BondId::new(index as u32),
            source_line,
        })
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
        },
    })
}
