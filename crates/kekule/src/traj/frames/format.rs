//! The text of a frame table: writing and parsing the header and the rows.
//!
//! `docs/frame-table.md` specifies the layout.

use std::fmt::Write as _;

use crate::properties::{PropertyColumn, PropertyKey, PropertyValue};
use crate::structure::{Sampling, ThermodynamicState};
use crate::units::{
    named_unit, Quantity, Unit, CANONICAL_PRESSURE_UNIT, CANONICAL_TEMPERATURE_UNIT,
};

use super::error::{FrameTableError, FrameTableErrorKind, Reason};
use super::header::{
    unit_symbol, ColumnSpec, ColumnType, FrameTableHeader, FRAME_TABLE_VERSION, STEP,
};

/// The first field of a frame table's first line.
const MAGIC: &str = "#kekule-frames";

/// The text of a missing text cell.
const MISSING_TEXT: &str = "\\N";

/// The text of `header`: its records and the line naming the columns.
pub(super) fn header_text(header: &FrameTableHeader) -> String {
    let mut text = format!("{MAGIC}\t{FRAME_TABLE_VERSION}\n#column\t{STEP}\tint\n");
    for column in header.columns() {
        let _ = write!(text, "#column\t{}\t", column.key());
        match column.kind() {
            ColumnType::Bool => text.push_str("bool"),
            ColumnType::Int => text.push_str("int"),
            ColumnType::Real(unit) => {
                let symbol = unit_symbol(unit).expect("a header's real columns have named units");
                let _ = write!(text, "real\t{symbol}");
            }
            ColumnType::Text => text.push_str("text"),
        }
        text.push('\n');
    }
    if let Some(time) = header.time() {
        let _ = writeln!(text, "#time\t{time}");
    }
    if let Some(sampling) = header.sampling() {
        let state = sampling.state();
        text.push_str("#state\ttemperature\t");
        real(&mut text, *state.temperature().value());
        text.push_str("\tK\n");
        if let Some(pressure) = state.pressure() {
            text.push_str("#state\tpressure\t");
            real(&mut text, *pressure.value());
            text.push_str("\tbar\n");
        }
        let _ = writeln!(text, "#sampling\tfrom_step\t{}", sampling.from_step());
        if let Some(column) = sampling.potential_energy() {
            let _ = writeln!(text, "#sampling\tpotential_energy\t{column}");
        }
        if let Some(method) = sampling.method() {
            text.push_str("#sampling\tmethod\t");
            escape(&mut text, method);
            text.push('\n');
        }
        for column in sampling.bias() {
            let _ = writeln!(text, "#bias\t{column}");
        }
    }
    text.push_str(STEP);
    for column in header.columns() {
        let _ = write!(text, "\t{}", column.key());
    }
    text.push('\n');
    text
}

/// Appends the line of one row to `line`.
pub(super) fn row_text(
    line: &mut String,
    step: u64,
    columns: &[ColumnSpec],
    cells: &[Option<PropertyValue>],
) {
    let _ = write!(line, "{step}");
    for (column, cell) in columns.iter().zip(cells) {
        line.push('\t');
        match cell {
            None if column.kind() == ColumnType::Text => line.push_str(MISSING_TEXT),
            None => {}
            Some(PropertyValue::Bool(value)) => {
                line.push_str(if *value { "true" } else { "false" })
            }
            Some(PropertyValue::Int(value)) => {
                let _ = write!(line, "{value}");
            }
            Some(PropertyValue::Real { value, .. }) => real(line, *value),
            Some(PropertyValue::String(value)) => escape(line, value),
        }
    }
    line.push('\n');
}

/// Writes `value` in Rust's shortest round-trip form: plain for zero and
/// magnitudes from 1e-4 to 1e16, in scientific notation otherwise.
fn real(text: &mut String, value: f64) {
    let magnitude = value.abs();
    if value == 0.0 || (1e-4..1e16).contains(&magnitude) {
        let _ = write!(text, "{value}");
    } else {
        let _ = write!(text, "{value:e}");
    }
}

/// Writes `value` with tabs, line breaks, and backslashes escaped.
fn escape(text: &mut String, value: &str) {
    for character in value.chars() {
        match character {
            '\t' => text.push_str("\\t"),
            '\n' => text.push_str("\\n"),
            '\r' => text.push_str("\\r"),
            '\\' => text.push_str("\\\\"),
            other => text.push(other),
        }
    }
}

/// The text an escaped field holds, or `None` for an unknown escape.
fn unescape(field: &str) -> Option<String> {
    let mut text = String::with_capacity(field.len());
    let mut characters = field.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            text.push(character);
            continue;
        }
        text.push(match characters.next()? {
            't' => '\t',
            'n' => '\n',
            'r' => '\r',
            '\\' => '\\',
            _ => return None,
        });
    }
    Some(text)
}

/// A real number of the cell grammar: finite.
fn parse_real(field: &str) -> Option<f64> {
    field.parse::<f64>().ok().filter(|value| value.is_finite())
}

/// The header of a frame table and where it ends.
#[derive(Debug)]
pub(super) struct ParsedHeader {
    pub(super) header: FrameTableHeader,
    /// The byte offset after the line naming the columns.
    pub(super) end: usize,
    /// The number of lines it spans.
    pub(super) lines: usize,
}

/// The fields a header's records state, before they are checked together.
#[derive(Default)]
struct Records {
    columns: Vec<ColumnSpec>,
    time: Option<PropertyKey>,
    temperature: Option<f64>,
    pressure: Option<f64>,
    from_step: Option<u64>,
    potential_energy: Option<PropertyKey>,
    method: Option<String>,
    bias: Vec<PropertyKey>,
    sampled: bool,
}

/// Parses the header at the start of `bytes`.
///
/// # Errors
///
/// `Pending` when the line naming the columns is not complete,
/// `Unsupported` for another version or an unknown record, and `Corrupt`
/// for a complete line that is not valid.
pub(super) fn parse_header(bytes: &[u8]) -> Result<ParsedHeader, FrameTableError> {
    let pending = || FrameTableError::new(FrameTableErrorKind::Pending, Reason::Pending);
    let mut offset = 0;
    let mut number = 0;
    let mut records = Records::default();
    loop {
        let Some(length) = bytes[offset..].iter().position(|&byte| byte == b'\n') else {
            // A first line that cannot become the magic and its tab is no
            // frame table; trailing NUL bytes are blocks lost at power loss.
            if number == 0 {
                let end = bytes
                    .iter()
                    .rposition(|&byte| byte != 0)
                    .map_or(0, |last| last + 1);
                let start = &bytes[..end];
                let lead = format!("{MAGIC}\t");
                if !(lead.as_bytes().starts_with(start) || start.starts_with(lead.as_bytes())) {
                    return Err(FrameTableError::corrupt(1, Reason::NotFrameTable));
                }
            }
            return Err(pending());
        };
        let line = std::str::from_utf8(&bytes[offset..offset + length])
            .map_err(|_| FrameTableError::corrupt(number + 1, Reason::Encoding))?;
        let line = line.strip_suffix('\r').unwrap_or(line);
        offset += length + 1;
        number += 1;
        if number == 1 {
            check_magic(line)?;
            continue;
        }
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            records.record(line, number)?;
            continue;
        }
        let header = records.header(number)?;
        let expected = std::iter::once(STEP)
            .chain(header.columns().iter().map(|column| column.key().as_str()));
        if !expected.eq(line.split('\t')) {
            return Err(FrameTableError::corrupt(number, Reason::Names));
        }
        return Ok(ParsedHeader {
            header,
            end: offset,
            lines: number,
        });
    }
}

fn check_magic(line: &str) -> Result<(), FrameTableError> {
    let mut fields = line.split('\t');
    if fields.next() != Some(MAGIC) {
        return Err(FrameTableError::corrupt(1, Reason::NotFrameTable));
    }
    match (fields.next(), fields.next()) {
        (Some(version), None) if version == FRAME_TABLE_VERSION.to_string() => Ok(()),
        (Some(version), None) => Err(FrameTableError::new(
            FrameTableErrorKind::Unsupported,
            Reason::Version(version.to_owned()),
        )),
        _ => Err(FrameTableError::corrupt(1, Reason::NotFrameTable)),
    }
}

impl Records {
    /// Takes one header record, line `number`.
    fn record(&mut self, line: &str, number: usize) -> Result<(), FrameTableError> {
        let fields = line.split('\t').collect::<Vec<_>>();
        let invalid = || FrameTableError::corrupt(number, Reason::Record);
        let repeated = || FrameTableError::corrupt(number, Reason::RepeatedRecord);
        let key = |text: &str| PropertyKey::new(text).map_err(|_| invalid());
        match fields.as_slice() {
            ["#column", name, kind, rest @ ..] => {
                let kind = match (*kind, rest) {
                    ("bool", []) => ColumnType::Bool,
                    ("int", []) => ColumnType::Int,
                    ("text", []) => ColumnType::Text,
                    ("real", [symbol]) => ColumnType::Real(named_unit(symbol).ok_or_else(invalid)?),
                    _ => return Err(invalid()),
                };
                // The first column record is the step column's.
                if self.columns.is_empty() && (*name != STEP || kind != ColumnType::Int) {
                    return Err(invalid());
                }
                self.columns.push(ColumnSpec::new(key(name)?, kind));
            }
            ["#time", name] => {
                if self.time.replace(key(name)?).is_some() {
                    return Err(repeated());
                }
            }
            ["#state", "temperature", value, unit] => {
                let value =
                    state_value(value, unit, CANONICAL_TEMPERATURE_UNIT).ok_or_else(invalid)?;
                self.sampled = true;
                if self.temperature.replace(value).is_some() {
                    return Err(repeated());
                }
            }
            ["#state", "pressure", value, unit] => {
                let value =
                    state_value(value, unit, CANONICAL_PRESSURE_UNIT).ok_or_else(invalid)?;
                self.sampled = true;
                if self.pressure.replace(value).is_some() {
                    return Err(repeated());
                }
            }
            ["#sampling", "from_step", value] => {
                let value = value.parse().map_err(|_| invalid())?;
                self.sampled = true;
                if self.from_step.replace(value).is_some() {
                    return Err(repeated());
                }
            }
            ["#sampling", "potential_energy", name] => {
                self.sampled = true;
                if self.potential_energy.replace(key(name)?).is_some() {
                    return Err(repeated());
                }
            }
            ["#sampling", "method", value] => {
                let value = unescape(value).ok_or_else(invalid)?;
                self.sampled = true;
                if self.method.replace(value).is_some() {
                    return Err(repeated());
                }
            }
            ["#bias", name] => {
                self.sampled = true;
                self.bias.push(key(name)?);
            }
            ["#column" | "#time" | "#bias", ..]
            | ["#state", "temperature" | "pressure", ..]
            | ["#sampling", "from_step" | "potential_energy" | "method", ..] => {
                return Err(invalid())
            }
            _ => {
                return Err(FrameTableError::new(
                    FrameTableErrorKind::Unsupported,
                    Reason::UnknownRecord,
                )
                .at_line(number))
            }
        }
        Ok(())
    }

    /// The header the records state; the names line is line `number`.
    fn header(&mut self, number: usize) -> Result<FrameTableHeader, FrameTableError> {
        let header_error = |error: FrameTableError| {
            FrameTableError::corrupt(number, Reason::Header).with_source(error)
        };
        if self.columns.is_empty() {
            return Err(FrameTableError::corrupt(number, Reason::Names));
        }
        let columns = self.columns.split_off(1);
        let mut header = FrameTableHeader::new(columns).map_err(header_error)?;
        if let Some(time) = self.time.take() {
            header = header.with_time(time).map_err(header_error)?;
        }
        if self.sampled {
            let sampling = self.sampling().map_err(|error| {
                FrameTableError::corrupt(number, Reason::Header).with_source(error)
            })?;
            let Some(sampling) = sampling else {
                return Err(FrameTableError::corrupt(number, Reason::Header));
            };
            header = header.with_sampling(sampling).map_err(header_error)?;
        }
        Ok(header)
    }

    /// The sampling record, or `None` when a required record is absent.
    fn sampling(&mut self) -> Result<Option<Sampling>, crate::structure::StructureSetError> {
        let (Some(temperature), Some(from_step)) = (self.temperature, self.from_step) else {
            return Ok(None);
        };
        let temperature = Quantity::new(temperature, CANONICAL_TEMPERATURE_UNIT);
        let state = match self.pressure {
            None => ThermodynamicState::canonical(temperature)?,
            Some(pressure) => ThermodynamicState::isothermal_isobaric(
                temperature,
                Quantity::new(pressure, CANONICAL_PRESSURE_UNIT),
            )?,
        };
        let mut sampling = Sampling::new(state).with_from_step(from_step);
        if let Some(method) = self.method.take() {
            sampling = sampling.with_method(method);
        }
        if let Some(column) = self.potential_energy.take() {
            sampling = sampling.with_potential_energy(column)?;
        }
        sampling.with_bias(std::mem::take(&mut self.bias)).map(Some)
    }
}

/// A state value in `unit`'s named unit, converted to `canonical`.
fn state_value(value: &str, unit: &str, canonical: Unit) -> Option<f64> {
    let unit = named_unit(unit).filter(|unit| unit.is_compatible(canonical))?;
    Quantity::new(parse_real(value)?, unit)
        .value_in(canonical)
        .ok()
}

/// The rows of a frame table: a step and a column of cells per header
/// column.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Rows {
    pub(super) steps: Vec<u64>,
    pub(super) columns: Vec<PropertyColumn>,
}

impl Rows {
    pub(super) fn new(header: &FrameTableHeader) -> Self {
        let columns = header
            .columns()
            .iter()
            .map(|column| match column.kind() {
                ColumnType::Bool => PropertyColumn::Bool(Vec::new()),
                ColumnType::Int => PropertyColumn::Int(Vec::new()),
                ColumnType::Real(unit) => PropertyColumn::Real {
                    unit,
                    values: Vec::new(),
                },
                ColumnType::Text => PropertyColumn::String(Vec::new()),
            })
            .collect();
        Self {
            steps: Vec::new(),
            columns,
        }
    }

    /// Parses the complete lines of `bytes`, which start at line `first`,
    /// and appends their rows; on error nothing is appended.
    pub(super) fn extend(
        &mut self,
        header: &FrameTableHeader,
        bytes: &[u8],
        first: usize,
    ) -> Result<usize, FrameTableError> {
        let mut staged = Self::new(header);
        let mut last = self.steps.last().copied();
        for (offset, line) in bytes.split_inclusive(|&byte| byte == b'\n').enumerate() {
            let number = first + offset;
            let line = line
                .strip_suffix(b"\n")
                .expect("only complete lines are parsed");
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            let line = std::str::from_utf8(line)
                .map_err(|_| FrameTableError::corrupt(number, Reason::Encoding))?;
            let step = staged.push(header, line, number)?;
            if last.is_some_and(|last| step <= last) {
                return Err(FrameTableError::corrupt(number, Reason::RowOrder));
            }
            last = Some(step);
        }
        let rows = staged.steps.len();
        self.steps.append(&mut staged.steps);
        for (column, staged) in self.columns.iter_mut().zip(staged.columns) {
            match (column, staged) {
                (PropertyColumn::Bool(values), PropertyColumn::Bool(mut more)) => {
                    values.append(&mut more);
                }
                (PropertyColumn::Int(values), PropertyColumn::Int(mut more)) => {
                    values.append(&mut more);
                }
                (
                    PropertyColumn::Real { values, .. },
                    PropertyColumn::Real {
                        values: mut more, ..
                    },
                ) => values.append(&mut more),
                (PropertyColumn::String(values), PropertyColumn::String(mut more)) => {
                    values.append(&mut more);
                }
                _ => unreachable!("staged rows have the columns' types"),
            }
        }
        Ok(rows)
    }

    /// Parses one row and appends it, returning its step.
    fn push(
        &mut self,
        header: &FrameTableHeader,
        line: &str,
        number: usize,
    ) -> Result<u64, FrameTableError> {
        let fields = line.split('\t').collect::<Vec<_>>();
        let expected = header.columns().len() + 1;
        if fields.len() != expected {
            return Err(FrameTableError::corrupt(
                number,
                Reason::CellCount {
                    found: fields.len(),
                    expected,
                },
            ));
        }
        let step = fields[0]
            .parse::<u64>()
            .map_err(|_| FrameTableError::corrupt(number, Reason::Step))?;
        // Check every cell before appending any.
        for (spec, field) in header.columns().iter().zip(&fields[1..]) {
            let valid = match spec.kind() {
                ColumnType::Bool => matches!(*field, "" | "true" | "false"),
                ColumnType::Int => field.is_empty() || field.parse::<i64>().is_ok(),
                ColumnType::Real(_) => field.is_empty() || parse_real(field).is_some(),
                ColumnType::Text => *field == MISSING_TEXT || unescape(field).is_some(),
            };
            if !valid {
                return Err(FrameTableError::corrupt(
                    number,
                    Reason::Cell(spec.key().clone()),
                ));
            }
        }
        self.steps.push(step);
        for (column, field) in self.columns.iter_mut().zip(&fields[1..]) {
            let field = *field;
            match column {
                PropertyColumn::Bool(values) => {
                    values.push((!field.is_empty()).then(|| field == "true"));
                }
                PropertyColumn::Int(values) => values.push(field.parse().ok()),
                PropertyColumn::Real { values, .. } => values.push(parse_real(field)),
                PropertyColumn::String(values) => {
                    values.push((field != MISSING_TEXT).then(|| unescape(field).expect("checked")));
                }
            }
        }
        Ok(step)
    }
}
