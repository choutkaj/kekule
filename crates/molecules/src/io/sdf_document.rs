use std::fmt;

use crate::small::model::SmallMolecule;

use super::{
    interpret_molfile_document, parse_molfile_document_with_options, MolfileDocument,
    MolfileInterpretationReport, SdfParseError, SdfParseOptions,
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

#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecordDocument {
    molfile: MolfileDocument,
    data_fields: Vec<SdfDataField>,
    source_start_line: usize,
}

impl SdfRecordDocument {
    pub fn molfile(&self) -> &MolfileDocument {
        &self.molfile
    }

    pub fn data_fields(&self) -> &[SdfDataField] {
        &self.data_fields
    }

    pub const fn source_start_line(&self) -> usize {
        self.source_start_line
    }
}

#[derive(Debug, Clone, PartialEq)]
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
    pub(crate) record: usize,
    pub(crate) line: usize,
    pub(crate) message: String,
}

impl SdfInterpretError {
    pub const fn record(&self) -> usize {
        self.record
    }

    pub const fn line(&self) -> usize {
        self.line
    }

    pub fn message(&self) -> &str {
        &self.message
    }
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfRecordInterpretationReport {
    record: usize,
    source_start_line: usize,
    molfile: MolfileInterpretationReport,
}

impl SdfRecordInterpretationReport {
    pub const fn record(&self) -> usize {
        self.record
    }

    pub const fn source_start_line(&self) -> usize {
        self.source_start_line
    }

    pub fn molfile(&self) -> &MolfileInterpretationReport {
        &self.molfile
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SdfInterpretationReport {
    records: Vec<SdfRecordInterpretationReport>,
}

impl SdfInterpretationReport {
    pub fn records(&self) -> &[SdfRecordInterpretationReport] {
        &self.records
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfInterpretation {
    records: Vec<SdfRecord>,
    report: SdfInterpretationReport,
}

impl SdfInterpretation {
    pub fn records(&self) -> &[SdfRecord] {
        &self.records
    }

    pub fn report(&self) -> &SdfInterpretationReport {
        &self.report
    }

    pub fn into_records(self) -> Vec<SdfRecord> {
        self.records
    }

    pub fn into_parts(self) -> (Vec<SdfRecord>, SdfInterpretationReport) {
        (self.records, self.report)
    }
}

pub fn parse_sdf_document(
    input: &str,
    options: SdfParseOptions,
) -> Result<SdfDocument, SdfParseError> {
    if input.len() > options.max_input_bytes {
        return Err(SdfParseError::new(
            1,
            1,
            "input exceeds configured byte limit",
        ));
    }
    let normalized = input.replace("\r\n", "\n").replace('\r', "\n");
    let mut records = Vec::new();
    let mut current = Vec::<(usize, &str)>::new();
    let mut current_bytes = 0usize;
    for (offset, line) in normalized.lines().enumerate() {
        let line_number = offset + 1;
        if line.trim() == "$$$$" {
            if current.iter().any(|(_, line)| !line.trim().is_empty()) {
                push_record_document(&mut records, &current, options)?;
            }
            current.clear();
            current_bytes = 0;
        } else {
            current_bytes = current_bytes
                .checked_add(line.len())
                .and_then(|bytes| bytes.checked_add(1))
                .ok_or_else(|| {
                    SdfParseError::new(
                        records.len() + 1,
                        line_number,
                        "SDF record byte count overflow",
                    )
                })?;
            if current_bytes > options.max_record_bytes {
                return Err(SdfParseError::new(
                    records.len() + 1,
                    line_number,
                    "SDF record exceeds configured byte limit",
                ));
            }
            current.push((line_number, line));
        }
    }
    if current.iter().any(|(_, line)| !line.trim().is_empty()) {
        if options.allow_missing_final_delimiter {
            push_record_document(&mut records, &current, options)?;
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

fn push_record_document(
    records: &mut Vec<SdfRecordDocument>,
    lines: &[(usize, &str)],
    options: SdfParseOptions,
) -> Result<(), SdfParseError> {
    if records.len() >= options.max_records {
        return Err(SdfParseError::new(
            records.len() + 1,
            lines.first().map(|(line, _)| *line).unwrap_or(1),
            "SDF record count exceeds configured limit",
        ));
    }
    records.push(parse_record_document(records.len() + 1, lines, options)?);
    Ok(())
}

fn parse_record_document(
    record: usize,
    lines: &[(usize, &str)],
    options: SdfParseOptions,
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
    let molfile = parse_molfile_document_with_options(
        &molfile_source,
        super::MolfileParseOptions {
            max_input_bytes: options.max_record_bytes,
            ..super::MolfileParseOptions::default()
        },
    )
    .map_err(|error| SdfParseError::new(record, lines[0].0 + error.line - 1, error.message))?;
    let mut data_fields = Vec::new();
    let mut index = end + 1;
    while index < lines.len() {
        let (line_number, line) = lines[index];
        if line.trim().is_empty() {
            index += 1;
            continue;
        }
        if !line.trim_start().starts_with('>') {
            return Err(SdfParseError::new(
                record,
                line_number,
                "unexpected content outside an SDF data field",
            ));
        }
        let name = sdf_field_name(line).ok_or_else(|| {
            SdfParseError::new(record, line_number, "invalid SDF data field header")
        })?;
        index += 1;
        let mut values = Vec::new();
        let mut terminated = false;
        while index < lines.len() {
            if lines[index].1.is_empty() {
                index += 1;
                terminated = true;
                break;
            }
            values.push(lines[index].1);
            index += 1;
        }
        if !terminated {
            return Err(SdfParseError::new(
                record,
                lines.last().map(|(line, _)| *line).unwrap_or(line_number),
                "SDF data field is missing its terminating blank line",
            ));
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
        source_start_line: lines.first().map(|(line, _)| *line).unwrap_or(1),
    })
}

fn sdf_field_name(line: &str) -> Option<String> {
    let start = line.find('<')? + 1;
    let end = line[start..].find('>')? + start;
    let name = line[start..end].trim();
    (!name.is_empty()).then(|| name.to_owned())
}

pub fn interpret_sdf_document(
    document: &SdfDocument,
) -> Result<SdfInterpretation, SdfInterpretError> {
    let mut records = Vec::with_capacity(document.records.len());
    let mut reports = Vec::with_capacity(document.records.len());
    for (index, record) in document.records.iter().enumerate() {
        let interpretation =
            interpret_molfile_document(&record.molfile).map_err(|error| SdfInterpretError {
                record: index + 1,
                line: record.source_start_line + error.line.saturating_sub(1),
                message: error.message,
            })?;
        let (molecule, molfile) = interpretation.into_parts();
        records.push(SdfRecord::new(
            record.molfile.header().title(),
            molecule,
            record.data_fields.clone(),
        ));
        reports.push(SdfRecordInterpretationReport {
            record: index + 1,
            source_start_line: record.source_start_line,
            molfile,
        });
    }
    Ok(SdfInterpretation {
        records,
        report: SdfInterpretationReport { records: reports },
    })
}
