use std::fmt;

use crate::small::SmallMolecule;

use super::{
    interpret_molfile_document, parse_molfile_document, MolfileDocument, SdfParseError,
    SdfParseOptions,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfDataField {
    name: String,
    value: String,
    line: usize,
}

impl SdfDataField {
    pub fn new(name: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            value: value.into(),
            line: 0,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn value(&self) -> &str {
        &self.value
    }

    pub const fn line(&self) -> usize {
        self.line
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfRecordDocument {
    molfile: MolfileDocument,
    data_fields: Vec<SdfDataField>,
}

impl SdfRecordDocument {
    pub fn molfile(&self) -> &MolfileDocument {
        &self.molfile
    }

    pub fn data_fields(&self) -> &[SdfDataField] {
        &self.data_fields
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfDocument {
    records: Vec<SdfRecordDocument>,
}

impl SdfDocument {
    pub fn records(&self) -> &[SdfRecordDocument] {
        &self.records
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    title: String,
    molecule: SmallMolecule,
    data_fields: Vec<SdfDataField>,
}

impl SdfRecord {
    pub fn new(
        title: impl Into<String>,
        molecule: SmallMolecule,
        data_fields: Vec<SdfDataField>,
    ) -> Self {
        Self {
            title: title.into(),
            molecule,
            data_fields,
        }
    }

    pub fn title(&self) -> &str {
        &self.title
    }

    pub fn molecule(&self) -> &SmallMolecule {
        &self.molecule
    }

    pub fn into_molecule(self) -> SmallMolecule {
        self.molecule
    }

    pub fn data_fields(&self) -> &[SdfDataField] {
        &self.data_fields
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfInterpretError {
    pub record: usize,
    pub line: usize,
    pub message: String,
}

impl fmt::Display for SdfInterpretError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SDF interpretation error in record {} at line {}: {}",
            self.record, self.line, self.message
        )
    }
}

impl std::error::Error for SdfInterpretError {}

pub fn parse_sdf_document(
    input: &str,
    options: SdfParseOptions,
) -> Result<SdfDocument, SdfParseError> {
    let normalized = input.replace("\r\n", "\n").replace('\r', "\n");
    let mut records = Vec::new();
    let mut current = Vec::<(usize, &str)>::new();
    let mut saw_delimiter = false;
    for (offset, line) in normalized.lines().enumerate() {
        let line_number = offset + 1;
        if line.trim() == "$$$$" {
            saw_delimiter = true;
            if current.iter().any(|(_, line)| !line.trim().is_empty()) {
                records.push(parse_record_document(records.len() + 1, &current)?);
            }
            current.clear();
        } else {
            current.push((line_number, line));
        }
    }
    if current.iter().any(|(_, line)| !line.trim().is_empty()) {
        if saw_delimiter || options.allow_missing_final_delimiter {
            records.push(parse_record_document(records.len() + 1, &current)?);
        } else {
            return Err(SdfParseError::new(
                records.len() + 1,
                current.last().map(|(line, _)| *line).unwrap_or(1),
                "missing final $$$$ record delimiter",
            ));
        }
    }
    Ok(SdfDocument { records })
}

fn parse_record_document(
    record: usize,
    lines: &[(usize, &str)],
) -> Result<SdfRecordDocument, SdfParseError> {
    let end = lines
        .iter()
        .position(|(_, line)| line.trim() == "M  END")
        .ok_or_else(|| {
            SdfParseError::new(
                record,
                lines.first().map(|(line, _)| *line).unwrap_or(1),
                "missing M  END",
            )
        })?;
    let molfile_source = lines[..=end]
        .iter()
        .map(|(_, line)| *line)
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let molfile = parse_molfile_document(&molfile_source)
        .map_err(|error| SdfParseError::new(record, lines[0].0 + error.line - 1, error.message))?;
    let mut data_fields = Vec::new();
    let mut index = end + 1;
    while index < lines.len() {
        let (line_number, line) = lines[index];
        if !line.trim_start().starts_with('>') {
            index += 1;
            continue;
        }
        let name = sdf_field_name(line).ok_or_else(|| {
            SdfParseError::new(record, line_number, "invalid SDF data field header")
        })?;
        index += 1;
        let mut values = Vec::new();
        while index < lines.len() && !lines[index].1.trim_start().starts_with('>') {
            if lines[index].1.is_empty() {
                index += 1;
                break;
            }
            values.push(lines[index].1);
            index += 1;
        }
        data_fields.push(SdfDataField {
            name,
            value: values.join("\n"),
            line: line_number,
        });
    }
    Ok(SdfRecordDocument {
        molfile,
        data_fields,
    })
}

fn sdf_field_name(line: &str) -> Option<String> {
    let start = line.find('<')? + 1;
    let end = line[start..].find('>')? + start;
    let name = line[start..end].trim();
    (!name.is_empty()).then(|| name.to_owned())
}

pub fn interpret_sdf_document(document: &SdfDocument) -> Result<Vec<SdfRecord>, SdfInterpretError> {
    document
        .records
        .iter()
        .enumerate()
        .map(|(index, record)| {
            let molecule =
                interpret_molfile_document(&record.molfile).map_err(|error| SdfInterpretError {
                    record: index + 1,
                    line: error.line,
                    message: error.message,
                })?;
            Ok(SdfRecord::new(
                record.molfile.header().title(),
                molecule,
                record.data_fields.clone(),
            ))
        })
        .collect()
}
