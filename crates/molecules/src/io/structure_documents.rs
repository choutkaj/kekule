use std::fmt;

use crate::small::SmallMolecule;

use super::{read_mol_v2000_str, read_mol_v3000_str, SdfParseError};

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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileDocument {
    source: String,
    version: MolfileVersion,
    header: MolfileHeader,
    atom_records: Vec<MolfileLine>,
    bond_records: Vec<MolfileLine>,
    property_records: Vec<MolfileLine>,
    unsupported_records: Vec<MolfileLine>,
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
    pub line: usize,
    pub message: String,
}

impl MolfileParseError {
    fn new(line: usize, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
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
    pub line: usize,
    pub message: String,
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
            let fields = lines[3].split_whitespace().collect::<Vec<_>>();
            let atom_count = fields
                .first()
                .and_then(|value| value.parse::<usize>().ok())
                .ok_or_else(|| MolfileParseError::new(4, "invalid V2000 atom count"))?;
            let bond_count = fields
                .get(1)
                .and_then(|value| value.parse::<usize>().ok())
                .ok_or_else(|| MolfileParseError::new(4, "invalid V2000 bond count"))?;
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
    Ok(MolfileDocument {
        source,
        version,
        header,
        atom_records,
        bond_records,
        property_records,
        unsupported_records,
    })
}

pub fn interpret_molfile_document(
    document: &MolfileDocument,
) -> Result<SmallMolecule, MolfileInterpretError> {
    match document.version {
        MolfileVersion::V2000 => read_mol_v2000_str(&document.source),
        MolfileVersion::V3000 => read_mol_v3000_str(&document.source),
    }
    .map_err(Into::into)
}
