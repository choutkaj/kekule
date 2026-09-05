use std::fmt;
use std::sync::Arc;

use crate::core::Molecule;
use crate::structure::Model;
use crate::topology::Topology;

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

/// One independently interpretable parsed SDF source record.
#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    molfile: MolfileDocument,
    data_fields: Vec<SdfDataField>,
    source_record_number: usize,
    source_start_line: usize,
}

impl SdfRecord {
    pub fn title(&self) -> &str {
        self.molfile.header().title()
    }

    pub fn molfile(&self) -> &MolfileDocument {
        &self.molfile
    }

    pub fn data_fields(&self) -> &[SdfDataField] {
        &self.data_fields
    }

    pub const fn source_record_number(&self) -> usize {
        self.source_record_number
    }

    pub const fn source_start_line(&self) -> usize {
        self.source_start_line
    }

    /// Interprets this independently meaningful record once into its richest
    /// canonical one-realization state.
    pub fn interpret(&self) -> Result<SdfRecordInterpretation, SdfInterpretError> {
        interpret_sdf_record(self)
    }

    /// Interprets this record as connected canonical molecules.
    ///
    /// Source coordinates may assist stereo normalization, but are not
    /// retained. No chemical perception is run.
    pub fn to_molecules(&self) -> Result<Vec<Molecule>, SdfInterpretError> {
        Ok(self.interpret()?.to_molecules())
    }

    /// Projects this record's complete static model layout, including the
    /// deterministic synthetic Molfile hierarchy.
    pub fn to_topology(&self) -> Result<Arc<Topology>, SdfInterpretError> {
        Ok(self.interpret()?.to_topology())
    }

    /// Interprets this record as one model with one instance per component.
    ///
    /// SDF data fields remain source metadata and are not copied into the
    /// canonical model. No chemical perception is run. Call [`Model::perceive`]
    /// on the returned model to explicitly install the default perception profile.
    pub fn to_model(&self) -> Result<Model, SdfInterpretError> {
        Ok(self.interpret()?.to_model())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfDocument {
    records: Vec<SdfRecord>,
}

impl SdfDocument {
    pub fn records(&self) -> &[SdfRecord] {
        &self.records
    }

    /// Interprets every record independently while preserving record boundaries.
    pub fn interpret(&self) -> Result<SdfInterpretation, SdfInterpretError> {
        interpret_sdf_document(self)
    }
}

/// Rich canonical interpretation of one SDF record.
///
/// The owned model is the single authoritative topology/geometry state. Its
/// molecule and topology projections never reinterpret the source.
#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecordInterpretation {
    title: String,
    model: Model,
    data_fields: Vec<SdfDataField>,
    report: SdfRecordInterpretationReport,
}

impl SdfRecordInterpretation {
    /// Creates a writable SDF record from an already canonical model.
    ///
    /// Source correspondence is empty because this value did not originate
    /// from SDF parsing. Parsed records receive a populated report through the
    /// authoritative interpreter.
    pub fn new(title: impl Into<String>, model: Model, data_fields: Vec<SdfDataField>) -> Self {
        Self {
            title: title.into(),
            model,
            data_fields,
            report: SdfRecordInterpretationReport::default(),
        }
    }

    fn with_report(
        title: impl Into<String>,
        model: Model,
        data_fields: Vec<SdfDataField>,
        report: SdfRecordInterpretationReport,
    ) -> Self {
        Self {
            title: title.into(),
            model,
            data_fields,
            report,
        }
    }

    pub fn title(&self) -> &str {
        &self.title
    }

    /// Borrows molecule definitions once per occurrence in authoritative
    /// topology instance order.
    pub fn molecules(&self) -> impl ExactSizeIterator<Item = &Molecule> + DoubleEndedIterator {
        self.model
            .topology()
            .molecules()
            .map(|occurrence| occurrence.molecule())
    }

    pub fn to_molecules(self) -> Vec<Molecule> {
        self.model
            .topology()
            .molecules()
            .map(|occurrence| occurrence.molecule().clone())
            .collect()
    }

    pub fn model(&self) -> &Model {
        &self.model
    }

    pub fn topology(&self) -> &Topology {
        self.model.topology()
    }

    pub fn to_model(self) -> Model {
        self.model
    }

    /// Consumes the format wrapper and retains shared ownership of its exact
    /// model topology.
    pub fn to_topology(self) -> Arc<Topology> {
        self.model.shared_topology()
    }

    pub fn data_fields(&self) -> &[SdfDataField] {
        &self.data_fields
    }

    pub fn report(&self) -> &SdfRecordInterpretationReport {
        &self.report
    }

    pub fn to_parts(
        self,
    ) -> (
        String,
        Model,
        Vec<SdfDataField>,
        SdfRecordInterpretationReport,
    ) {
        (self.title, self.model, self.data_fields, self.report)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfInterpretError {
    pub(crate) record: usize,
    pub(crate) line: usize,
    pub(crate) message: String,
    pub(crate) kind: SdfInterpretErrorKind,
}

/// Structured stage that caused SDF record interpretation to fail.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SdfInterpretErrorKind {
    Molfile(super::MolfileInterpretError),
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

    /// Returns the structured underlying interpretation failure.
    pub fn kind(&self) -> &SdfInterpretErrorKind {
        &self.kind
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

impl std::error::Error for SdfInterpretError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.kind {
            SdfInterpretErrorKind::Molfile(error) => Some(error),
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SdfRecordInterpretationReport {
    record: usize,
    source_start_line: usize,
    molfile_components: Vec<MolfileInterpretationReport>,
}

impl SdfRecordInterpretationReport {
    pub const fn record(&self) -> usize {
        self.record
    }

    pub const fn source_start_line(&self) -> usize {
        self.source_start_line
    }

    pub fn molfile_components(&self) -> &[MolfileInterpretationReport] {
        &self.molfile_components
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfInterpretation {
    records: Vec<SdfRecordInterpretation>,
}

impl SdfInterpretation {
    pub fn records(&self) -> &[SdfRecordInterpretation] {
        &self.records
    }

    /// Borrows each record's report in source order, without duplicate ownership.
    pub fn reports(
        &self,
    ) -> impl ExactSizeIterator<Item = &SdfRecordInterpretationReport> + DoubleEndedIterator {
        self.records.iter().map(SdfRecordInterpretation::report)
    }

    pub fn to_records(self) -> Vec<SdfRecordInterpretation> {
        self.records
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
                push_record_document(&mut records, &current, options, true)?;
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
            push_record_document(&mut records, &current, options, false)?;
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
    records: &mut Vec<SdfRecord>,
    lines: &[(usize, &str)],
    options: SdfParseOptions,
    ended_by_delimiter: bool,
) -> Result<(), SdfParseError> {
    if records.len() >= options.max_records {
        return Err(SdfParseError::new(
            records.len() + 1,
            lines.first().map(|(line, _)| *line).unwrap_or(1),
            "SDF record count exceeds configured limit",
        ));
    }
    records.push(parse_record_document(
        records.len() + 1,
        lines,
        options,
        ended_by_delimiter,
    )?);
    Ok(())
}

fn parse_record_document(
    record: usize,
    lines: &[(usize, &str)],
    options: SdfParseOptions,
    ended_by_delimiter: bool,
) -> Result<SdfRecord, SdfParseError> {
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
        if !terminated && !ended_by_delimiter {
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
    Ok(SdfRecord {
        molfile,
        data_fields,
        source_record_number: record,
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
    for record in &document.records {
        let interpretation = record.interpret()?;
        records.push(interpretation);
    }
    Ok(SdfInterpretation { records })
}

fn interpret_sdf_record(record: &SdfRecord) -> Result<SdfRecordInterpretation, SdfInterpretError> {
    let (model, molfile_components) = interpret_sdf_record_molfile(record)?.to_parts();
    let report = SdfRecordInterpretationReport {
        record: record.source_record_number,
        source_start_line: record.source_start_line,
        molfile_components,
    };
    Ok(SdfRecordInterpretation::with_report(
        record.title(),
        model,
        record.data_fields.clone(),
        report,
    ))
}

fn interpret_sdf_record_molfile(
    record: &SdfRecord,
) -> Result<super::MolfileInterpretation, SdfInterpretError> {
    interpret_molfile_document(&record.molfile).map_err(|error| SdfInterpretError {
        record: record.source_record_number,
        line: record.source_start_line + error.line.saturating_sub(1),
        message: error.message.clone(),
        kind: SdfInterpretErrorKind::Molfile(error),
    })
}
