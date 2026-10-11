use std::fmt;
use std::io::{BufReader, Read};
use std::sync::Arc;

use crate::core::Molecule;
use crate::diagnostics::{read_bounded, DiagnosticKind, SourceLocation};
use crate::structure::Model;
use crate::topology::Topology;

use super::lines::read_line;
use super::molfile_document::{
    parse_molfile_lines, MolfileDocument, MolfileLines, MolfileParseError, MolfileParseOptions,
    MolfileUnsupported,
};
use super::structure_documents::{
    MolfileInterpretError, MolfileInterpretOptions, MolfileInterpretationReport,
};

/// Resource and record-boundary policy for SDF parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SdfParseOptions {
    /// Accept a nonempty final record without its terminating `$$$$` line.
    pub allow_missing_final_delimiter: bool,
    /// Maximum UTF-8 byte length of a complete input parsed as one document.
    /// Streaming record readers hold one record at a time and ignore it.
    pub max_input_bytes: usize,
    /// Maximum number of nonempty records.
    pub max_records: usize,
    /// Maximum normalized byte length of one record, excluding its delimiter.
    pub max_record_bytes: usize,
}

impl Default for SdfParseOptions {
    fn default() -> Self {
        Self {
            allow_missing_final_delimiter: false,
            max_input_bytes: 256 * 1024 * 1024,
            max_records: 1_000_000,
            max_record_bytes: 64 * 1024 * 1024,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfParseError {
    record: usize,
    kind: DiagnosticKind,
    location: SourceLocation,
    message: String,
}

impl SdfParseError {
    fn new(record: usize, kind: DiagnosticKind, line: usize, message: impl Into<String>) -> Self {
        Self {
            record,
            kind,
            location: SourceLocation::at_line(line),
            message: message.into(),
        }
    }

    fn molfile(record: usize, error: MolfileParseError) -> Self {
        Self {
            record,
            kind: error.kind(),
            location: error.location(),
            message: error.message().to_owned(),
        }
    }

    /// The 1-based number of the failing record among nonempty records.
    pub const fn record(&self) -> usize {
        self.record
    }

    pub const fn kind(&self) -> DiagnosticKind {
        self.kind
    }

    pub const fn location(&self) -> SourceLocation {
        self.location
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for SdfParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SDF parse error in record {} at {}: {}",
            self.record, self.location, self.message
        )
    }
}

impl std::error::Error for SdfParseError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfDataField {
    name: String,
    value: String,
    location: Option<SourceLocation>,
}

impl SdfDataField {
    pub fn new(name: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            value: value.into(),
            location: None,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn value(&self) -> &str {
        &self.value
    }

    /// The field header line of a parsed field.
    pub const fn location(&self) -> Option<SourceLocation> {
        self.location
    }
}

/// One independently interpretable parsed SDF source record.
#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    molfile: MolfileDocument,
    data_fields: Vec<SdfDataField>,
    number: usize,
    location: SourceLocation,
}

impl SdfRecord {
    pub fn title(&self) -> &str {
        self.molfile.header().title()
    }

    /// The record's Molfile, located in the SDF input.
    pub fn molfile(&self) -> &MolfileDocument {
        &self.molfile
    }

    pub fn data_fields(&self) -> &[SdfDataField] {
        &self.data_fields
    }

    /// The 1-based number of this record among nonempty records.
    pub const fn number(&self) -> usize {
        self.number
    }

    /// The first line of the record.
    pub const fn location(&self) -> SourceLocation {
        self.location
    }

    /// Interprets this independently meaningful record once into its richest
    /// canonical one-structure state, rejecting unsupported content.
    pub fn interpret(&self) -> Result<SdfRecordInterpretation, SdfInterpretError> {
        self.interpret_with_options(MolfileInterpretOptions::default())
    }

    /// Interprets this record with an explicit unsupported-content policy.
    pub fn interpret_with_options(
        &self,
        options: MolfileInterpretOptions,
    ) -> Result<SdfRecordInterpretation, SdfInterpretError> {
        let (model, molfile_components, omitted) = self
            .molfile
            .interpret_with_options(options)
            .map_err(|error| SdfInterpretError {
                record: self.number,
                error,
            })?
            .into_parts();
        let report = SdfRecordInterpretationReport {
            record: self.number,
            location: Some(self.location),
            molfile_components,
            omitted,
        };
        Ok(SdfRecordInterpretation {
            title: self.title().to_owned(),
            model,
            data_fields: self.data_fields.clone(),
            report,
        })
    }

    /// Interprets this record as connected canonical molecules.
    ///
    /// Source coordinates may assist stereo normalization, but are not
    /// retained. No chemical perception is run.
    pub fn to_molecules(&self) -> Result<Vec<Molecule>, SdfInterpretError> {
        Ok(self.interpret()?.into_molecules())
    }

    /// Projects this record's complete static model layout, including the
    /// deterministic synthetic Molfile hierarchy.
    pub fn to_topology(&self) -> Result<Arc<Topology>, SdfInterpretError> {
        Ok(self.interpret()?.into_topology())
    }

    /// Interprets this record as one model with one instance per component.
    ///
    /// SDF data fields remain source metadata and are not copied into the
    /// canonical model. No chemical perception is run. Call [`Model::perceive`]
    /// on the returned model to explicitly install the default perception profile.
    pub fn to_model(&self) -> Result<Model, SdfInterpretError> {
        Ok(self.interpret()?.into_model())
    }
}

/// Every record of one SDF input, parsed. Interpret records individually;
/// a record's interpretation never depends on its siblings.
#[derive(Debug, Clone, PartialEq)]
pub struct SdfDocument {
    records: Vec<SdfRecord>,
}

impl SdfDocument {
    pub fn records(&self) -> &[SdfRecord] {
        &self.records
    }

    pub fn into_records(self) -> Vec<SdfRecord> {
        self.records
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

    pub fn into_molecules(self) -> Vec<Molecule> {
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

    pub fn into_model(self) -> Model {
        self.model
    }

    /// Consumes the format wrapper and retains shared ownership of its exact
    /// model topology.
    pub fn into_topology(self) -> Arc<Topology> {
        self.model.shared_topology()
    }

    pub fn data_fields(&self) -> &[SdfDataField] {
        &self.data_fields
    }

    pub fn report(&self) -> &SdfRecordInterpretationReport {
        &self.report
    }

    pub fn into_parts(
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

/// Failure to interpret one SDF record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfInterpretError {
    record: usize,
    error: MolfileInterpretError,
}

impl SdfInterpretError {
    /// The 1-based number of the failing record.
    pub const fn record(&self) -> usize {
        self.record
    }

    pub const fn kind(&self) -> DiagnosticKind {
        self.error.kind()
    }

    /// The SDF input line the failure concerns, when one does.
    pub const fn location(&self) -> Option<SourceLocation> {
        self.error.location()
    }

    pub fn message(&self) -> &str {
        self.error.message()
    }
}

impl fmt::Display for SdfInterpretError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SDF interpretation error in record {}", self.record)?;
        if let Some(location) = self.location() {
            write!(f, " at {location}")?;
        }
        write!(f, ": {}", self.message())
    }
}

impl std::error::Error for SdfInterpretError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.error)
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SdfRecordInterpretationReport {
    record: usize,
    location: Option<SourceLocation>,
    molfile_components: Vec<MolfileInterpretationReport>,
    omitted: Vec<MolfileUnsupported>,
}

impl SdfRecordInterpretationReport {
    /// The 1-based source record number; 0 for records not produced by parsing.
    pub const fn record(&self) -> usize {
        self.record
    }

    /// The first line of the source record.
    pub const fn location(&self) -> Option<SourceLocation> {
        self.location
    }

    pub fn molfile_components(&self) -> &[MolfileInterpretationReport] {
        &self.molfile_components
    }

    /// Unsupported Molfile content omitted under
    /// [`UnsupportedContent::Report`](crate::diagnostics::UnsupportedContent::Report).
    pub fn omitted(&self) -> &[MolfileUnsupported] {
        &self.omitted
    }
}

pub fn parse_sdf_document(
    input: &str,
    options: SdfParseOptions,
) -> Result<SdfDocument, SdfParseError> {
    if input.len() > options.max_input_bytes {
        return Err(SdfParseError::new(
            1,
            DiagnosticKind::Limit,
            1,
            "input exceeds configured byte limit",
        ));
    }
    let records = SdfRecords::new(input.as_bytes(), options).collect::<Result<_, _>>()?;
    Ok(SdfDocument { records })
}

/// Reads and parses a complete SDF input, enforcing the byte limit while reading.
pub fn parse_sdf_reader(
    reader: impl Read,
    options: SdfParseOptions,
) -> Result<SdfDocument, SdfParseError> {
    let input = read_bounded(reader, options.max_input_bytes)
        .map_err(|(kind, message)| SdfParseError::new(1, kind, 1, message))?;
    parse_sdf_document(&input, options)
}

/// Streams SDF records one at a time.
///
/// Records are split at `$$$$` lines without reading their content, so a
/// malformed, oversized, or non-UTF-8 record yields one error and reading
/// resumes at the next record. Memory is bounded by one record. A failure to
/// read the input, or exceeding the record-count limit, ends the stream.
#[derive(Debug)]
pub struct SdfRecords<R> {
    reader: BufReader<R>,
    options: SdfParseOptions,
    line: Vec<u8>,
    lines_read: usize,
    records: usize,
    done: bool,
}

impl<R: Read> SdfRecords<R> {
    pub fn new(reader: R, options: SdfParseOptions) -> Self {
        Self {
            reader: BufReader::new(reader),
            options,
            line: Vec::new(),
            lines_read: 0,
            records: 0,
            done: false,
        }
    }

    fn record_error(
        &self,
        kind: DiagnosticKind,
        line: usize,
        message: impl Into<String>,
    ) -> SdfParseError {
        SdfParseError::new(self.records + 1, kind, line, message)
    }

    /// Numbers a completed record and enforces the record-count limit.
    fn finish(&mut self, first_line: usize) -> Result<usize, SdfParseError> {
        if self.records >= self.options.max_records {
            self.done = true;
            return Err(self.record_error(
                DiagnosticKind::Limit,
                first_line,
                "SDF record count exceeds configured limit",
            ));
        }
        self.records += 1;
        Ok(self.records)
    }
}

impl<R: Read> Iterator for SdfRecords<R> {
    type Item = Result<SdfRecord, SdfParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let mut text = String::new();
        let mut lines = Vec::<std::ops::Range<usize>>::new();
        let mut first_line = self.lines_read + 1;
        let mut nonblank = false;
        let mut failure = None::<SdfParseError>;
        let mut bytes = 0usize;
        loop {
            let truncated = match read_line(
                &mut self.reader,
                &mut self.line,
                self.options.max_record_bytes,
            ) {
                Ok(Some(truncated)) => truncated,
                Ok(None) => {
                    self.done = true;
                    if !nonblank && failure.is_none() {
                        return None;
                    }
                    let record = self.finish(first_line);
                    return Some(match (record, failure) {
                        (Err(error), _) | (Ok(_), Some(error)) => Err(error),
                        (Ok(record), None) if self.options.allow_missing_final_delimiter => {
                            parse_sdf_record(record, &text, &lines, first_line, false)
                        }
                        (Ok(record), None) => Err(SdfParseError::new(
                            record,
                            DiagnosticKind::Syntax,
                            self.lines_read,
                            "missing final $$$$ record delimiter",
                        )),
                    });
                }
                Err(error) => {
                    self.done = true;
                    return Some(Err(self.record_error(
                        DiagnosticKind::Io,
                        self.lines_read + 1,
                        format!("cannot read input: {error}"),
                    )));
                }
            };
            self.lines_read += 1;
            let number = self.lines_read;
            if !truncated && self.line.trim_ascii() == b"$$$$" {
                if !nonblank && failure.is_none() {
                    // Blank lines between delimiters are not a record.
                    text.clear();
                    lines.clear();
                    first_line = number + 1;
                    bytes = 0;
                    continue;
                }
                return Some(match (self.finish(first_line), failure) {
                    (Err(error), _) | (Ok(_), Some(error)) => Err(error),
                    (Ok(record), None) => parse_sdf_record(record, &text, &lines, first_line, true),
                });
            }
            if failure.is_some() {
                continue;
            }
            bytes = bytes.saturating_add(self.line.len() + 1);
            if truncated || bytes > self.options.max_record_bytes {
                failure = Some(self.record_error(
                    DiagnosticKind::Limit,
                    number,
                    "SDF record exceeds configured byte limit",
                ));
                continue;
            }
            match std::str::from_utf8(&self.line) {
                Ok(line) => {
                    nonblank |= !line.trim().is_empty();
                    let start = text.len();
                    text.push_str(line);
                    lines.push(start..text.len());
                }
                Err(error) => {
                    failure = Some(self.record_error(
                        DiagnosticKind::Io,
                        number,
                        format!("SDF record is not valid UTF-8: {error}"),
                    ));
                }
            }
        }
    }
}

impl<R: Read> std::iter::FusedIterator for SdfRecords<R> {}

fn parse_sdf_record(
    record: usize,
    text: &str,
    ranges: &[std::ops::Range<usize>],
    first_line: usize,
    ended_by_delimiter: bool,
) -> Result<SdfRecord, SdfParseError> {
    let lines = ranges
        .iter()
        .map(|range| &text[range.clone()])
        .collect::<Vec<_>>();
    let molfile_lines = MolfileLines::new(&lines, first_line);
    let (molfile, end) = parse_molfile_lines(molfile_lines, MolfileParseOptions::default())
        .map_err(|error| SdfParseError::molfile(record, error))?;
    let syntax = |index: usize, message: &str| {
        SdfParseError::new(
            record,
            DiagnosticKind::Syntax,
            molfile_lines.number(index),
            message,
        )
    };
    let mut data_fields = Vec::new();
    let mut index = end + 1;
    while index < lines.len() {
        let header = index;
        let line = lines[index];
        if line.trim().is_empty() {
            index += 1;
            continue;
        }
        if !line.trim_start().starts_with('>') {
            return Err(syntax(
                index,
                "unexpected content outside an SDF data field",
            ));
        }
        let name =
            sdf_field_name(line).ok_or_else(|| syntax(index, "invalid SDF data field header"))?;
        index += 1;
        let mut values = Vec::new();
        let mut terminated = false;
        while index < lines.len() {
            if lines[index].is_empty() {
                index += 1;
                terminated = true;
                break;
            }
            values.push(lines[index]);
            index += 1;
        }
        if !terminated && !ended_by_delimiter {
            return Err(syntax(
                lines.len() - 1,
                "SDF data field is missing its terminating blank line",
            ));
        }
        data_fields.push(SdfDataField {
            name,
            value: values.join("\n"),
            location: Some(molfile_lines.location(header)),
        });
    }
    Ok(SdfRecord {
        molfile,
        data_fields,
        number: record,
        location: SourceLocation::at_line(first_line),
    })
}

fn sdf_field_name(line: &str) -> Option<String> {
    let start = line.find('<')? + 1;
    let end = line[start..].find('>')? + start;
    let name = line[start..end].trim();
    (!name.is_empty()).then(|| name.to_owned())
}
