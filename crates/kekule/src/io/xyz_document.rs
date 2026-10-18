//! XYZ files as typed frames.
//!
//! An XYZ file is a sequence of frames. Each frame is an atom-count line, a
//! free-text comment line, and one line per atom: an element token, three
//! Cartesian coordinates, and optionally further whitespace-separated columns,
//! as extended XYZ writes them. Parsing types the count, element token, and
//! coordinates and keeps the comment and any further columns verbatim;
//! interpretation decides what they mean.
use crate::diagnostics::{
    read_bounded, DiagnosticKind, Format, FormatError, LimitExceeded, Resource, SourceLocation,
};
use crate::geometry::Point3;

use super::xyz_interpret::{interpret_xyz, interpret_xyz_ensemble};
use super::{
    XyzEnsembleInterpretOptions, XyzEnsembleInterpretation, XyzInterpretOptions, XyzInterpretation,
};

/// Syntax and resource limits for parsing an XYZ file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct XyzParseOptions {
    max_input_bytes: usize,
    max_atom_records: usize,
}

impl XyzParseOptions {
    /// The [`Default`] bound of [`Self::max_input_bytes`].
    pub const DEFAULT_MAX_INPUT_BYTES: usize = 1024 * 1024 * 1024;
    /// The [`Default`] bound of [`Self::max_atom_records`].
    pub const DEFAULT_MAX_ATOM_RECORDS: usize = 5_000_000;

    /// Maximum UTF-8 byte length of the input.
    pub const fn max_input_bytes(&self) -> usize {
        self.max_input_bytes
    }

    /// Sets [`Self::max_input_bytes`].
    #[must_use]
    pub const fn with_max_input_bytes(mut self, max_input_bytes: usize) -> Self {
        self.max_input_bytes = max_input_bytes;
        self
    }

    /// Maximum number of atom lines over all frames.
    pub const fn max_atom_records(&self) -> usize {
        self.max_atom_records
    }

    /// Sets [`Self::max_atom_records`].
    #[must_use]
    pub const fn with_max_atom_records(mut self, max_atom_records: usize) -> Self {
        self.max_atom_records = max_atom_records;
        self
    }
}

impl Default for XyzParseOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: Self::DEFAULT_MAX_INPUT_BYTES,
            max_atom_records: Self::DEFAULT_MAX_ATOM_RECORDS,
        }
    }
}

/// One atom line of an XYZ frame.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzAtom {
    element: String,
    position: Point3,
    extra: Option<String>,
    line: usize,
}

impl XyzAtom {
    /// The element token as written: a symbol, matched case-insensitively at
    /// interpretation, or an atomic number.
    pub fn element(&self) -> &str {
        &self.element
    }

    /// The Cartesian coordinates as written, in the file's length unit.
    pub const fn position(&self) -> Point3 {
        self.position
    }

    /// The columns after the coordinates, verbatim from the first to the
    /// last non-whitespace character, or `None` when the line has exactly
    /// four columns.
    pub fn extra(&self) -> Option<&str> {
        self.extra.as_deref()
    }

    /// The source location of the atom line.
    pub const fn location(&self) -> SourceLocation {
        SourceLocation::at_line(self.line)
    }
}

/// One XYZ frame: its comment line and its atom lines in source order.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzFrame {
    comment: String,
    atoms: Vec<XyzAtom>,
    line: usize,
}

impl XyzFrame {
    /// The comment line verbatim, without its line terminator.
    pub fn comment(&self) -> &str {
        &self.comment
    }

    /// The atom lines in source order.
    pub fn atoms(&self) -> &[XyzAtom] {
        &self.atoms
    }

    /// The source location of the frame's atom-count line.
    pub const fn location(&self) -> SourceLocation {
        SourceLocation::at_line(self.line)
    }
}

/// A parsed XYZ file: its frames in source order.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzDocument {
    frames: Vec<XyzFrame>,
}

impl XyzDocument {
    /// The frames in source order.
    pub fn frames(&self) -> &[XyzFrame] {
        &self.frames
    }

    /// Interprets one frame as a [`Model`](crate::structure::Model) with
    /// connectivity inferred from interatomic distances.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of format
    /// [`Format::Xyz`]: `InvalidValue` for an unknown element, `Unsupported`
    /// for further columns unless `options` reports them, `Inconsistent` when
    /// the frame selection does not match the document, `Limit` when an atom
    /// has more contacts than `options` allows, and `Chemistry` when the
    /// structure cannot be published.
    pub fn interpret(
        &self,
        options: XyzInterpretOptions,
    ) -> Result<XyzInterpretation, FormatError> {
        interpret_xyz(self, &options)
    }

    /// Interprets several frames over the same atoms as one equally weighted
    /// [`Ensemble`](crate::structure::Ensemble). Each frame infers its own
    /// connectivity, and every frame must infer the same bonds as the first.
    ///
    /// # Errors
    ///
    /// As [`Self::interpret`]. A selected frame whose atom count or element
    /// sequence differs from the first selected frame fails with
    /// [`FormatErrorDetail::InconsistentAtomSet`](crate::diagnostics::FormatErrorDetail::InconsistentAtomSet),
    /// and one whose distances infer different bonds fails with
    /// [`FormatErrorDetail::InconsistentTopology`](crate::diagnostics::FormatErrorDetail::InconsistentTopology);
    /// both name the frame by its zero-based index as the model ID.
    pub fn interpret_ensemble(
        &self,
        options: XyzEnsembleInterpretOptions,
    ) -> Result<XyzEnsembleInterpretation, FormatError> {
        interpret_xyz_ensemble(self, &options)
    }
}

fn xyz_parse_error(kind: DiagnosticKind, line: usize, message: impl Into<String>) -> FormatError {
    FormatError::parse(Format::Xyz, kind, message).at(SourceLocation::at_line(line))
}

pub(crate) fn parse_xyz_reader(
    reader: impl std::io::Read,
    options: XyzParseOptions,
) -> Result<XyzDocument, FormatError> {
    let input = read_bounded(reader, options.max_input_bytes, Format::Xyz)
        .map_err(|error| error.at(SourceLocation::at_line(1)))?;
    parse_xyz_str(&input, options)
}

pub(crate) fn parse_xyz_str(
    input: &str,
    options: XyzParseOptions,
) -> Result<XyzDocument, FormatError> {
    if input.len() > options.max_input_bytes {
        return Err(xyz_parse_error(
            DiagnosticKind::Limit,
            1,
            "input exceeds configured byte limit",
        )
        .exceeding(LimitExceeded::new(
            Resource::InputBytes,
            options.max_input_bytes,
            input.len(),
        )));
    }
    // Every terminator (`\r\n`, `\r`, `\n`) ends one line; a final terminator
    // does not start another.
    let body = input.strip_suffix('\n').unwrap_or(input);
    let body = body.strip_suffix('\r').unwrap_or(body);
    let mut lines = (!input.is_empty())
        .then(|| {
            body.split('\n')
                .flat_map(|chunk| chunk.strip_suffix('\r').unwrap_or(chunk).split('\r'))
        })
        .into_iter()
        .flatten()
        .enumerate()
        .map(|(index, line)| (index + 1, line))
        .peekable();
    let mut frames = Vec::new();
    let mut atom_records = 0usize;
    while let Some((count_line, count_text)) = lines.next() {
        if count_text.trim().is_empty() {
            // Blank lines may only trail the last frame.
            if let Some((line, _)) = lines.find(|(_, text)| !text.trim().is_empty()) {
                return Err(xyz_parse_error(
                    DiagnosticKind::Syntax,
                    line,
                    "content follows a blank line after the last frame",
                ));
            }
            break;
        }
        let count = count_text.trim().parse::<usize>().map_err(|_| {
            xyz_parse_error(
                DiagnosticKind::Syntax,
                count_line,
                "the atom-count line is not one unsigned integer",
            )
        })?;
        if count == 0 {
            return Err(xyz_parse_error(
                DiagnosticKind::InvalidValue,
                count_line,
                "an XYZ frame cannot contain zero atoms",
            ));
        }
        let total = atom_records.saturating_add(count);
        if total > options.max_atom_records {
            return Err(xyz_parse_error(
                DiagnosticKind::Limit,
                count_line,
                "atom lines exceed the configured limit",
            )
            .exceeding(LimitExceeded::new(
                Resource::Atoms,
                options.max_atom_records,
                total,
            )));
        }
        atom_records = total;
        let Some((_, comment)) = lines.next() else {
            return Err(xyz_parse_error(
                DiagnosticKind::Syntax,
                count_line,
                "the frame ends before its comment line",
            ));
        };
        let mut atoms = Vec::with_capacity(count);
        for atom_index in 0..count {
            let Some((line, text)) = lines.next() else {
                return Err(xyz_parse_error(
                    DiagnosticKind::Syntax,
                    count_line,
                    format!("the frame declares {count} atoms but ends after {atom_index}"),
                ));
            };
            atoms.push(parse_atom(line, text)?);
        }
        frames.push(XyzFrame {
            comment: comment.to_owned(),
            atoms,
            line: count_line,
        });
    }
    Ok(XyzDocument { frames })
}

fn parse_atom(line: usize, text: &str) -> Result<XyzAtom, FormatError> {
    let mut rest = text.trim_start();
    let mut field = || {
        let end = rest.find(char::is_whitespace).unwrap_or(rest.len());
        let (value, tail) = rest.split_at(end);
        rest = tail.trim_start();
        (!value.is_empty()).then_some(value)
    };
    let (Some(element), Some(x), Some(y), Some(z)) = (field(), field(), field(), field()) else {
        return Err(xyz_parse_error(
            DiagnosticKind::Syntax,
            line,
            "an atom line needs an element and three coordinates",
        ));
    };
    let coordinate = |value: &str| {
        value
            .parse::<f64>()
            .ok()
            .filter(|value| value.is_finite())
            .ok_or_else(|| {
                xyz_parse_error(
                    DiagnosticKind::InvalidValue,
                    line,
                    format!("`{value}` is not a finite coordinate"),
                )
            })
    };
    let position = Point3::new(coordinate(x)?, coordinate(y)?, coordinate(z)?);
    let extra = rest.trim_end();
    Ok(XyzAtom {
        element: element.to_owned(),
        position,
        extra: (!extra.is_empty()).then(|| extra.to_owned()),
        line,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frames_keep_comments_extra_columns_and_lines() {
        let text = "2\r\n water? \r\nO 0 0 0\r\nH 0.96 0 0 0.41 charge\r\n";
        let document = parse_xyz_str(text, XyzParseOptions::default()).unwrap();
        let [frame] = document.frames() else {
            panic!("one frame")
        };
        assert_eq!(frame.comment(), " water? ");
        assert_eq!(frame.atoms()[1].element(), "H");
        assert_eq!(frame.atoms()[1].position(), Point3::new(0.96, 0.0, 0.0));
        assert_eq!(frame.atoms()[1].extra(), Some("0.41 charge"));
        assert_eq!(frame.atoms()[0].extra(), None);
        assert_eq!(frame.atoms()[1].location().line(), 4);
    }

    #[test]
    fn malformed_frames_fail_with_typed_kinds_and_lines() {
        for (text, kind, line) in [
            ("x\n\n", DiagnosticKind::Syntax, 1),
            ("0\n\n", DiagnosticKind::InvalidValue, 1),
            ("2\n\nHe 0 0 0\n", DiagnosticKind::Syntax, 1),
            ("1\n\nHe 0 0\n", DiagnosticKind::Syntax, 3),
            ("1\n\nHe 0 nan 0\n", DiagnosticKind::InvalidValue, 3),
            ("1\n\nHe 0 0 0\n\n1\n", DiagnosticKind::Syntax, 5),
            ("1\n", DiagnosticKind::Syntax, 1),
        ] {
            let error = parse_xyz_str(text, XyzParseOptions::default()).unwrap_err();
            assert_eq!(error.kind(), kind, "{text:?}");
            assert_eq!(error.location().map(|at| at.line()), Some(line), "{text:?}");
        }
        assert!(parse_xyz_str("", XyzParseOptions::default())
            .unwrap()
            .frames()
            .is_empty());
        assert_eq!(
            parse_xyz_str("1\n\nHe 0 0 0\n\n \n", XyzParseOptions::default())
                .unwrap()
                .frames()
                .len(),
            1
        );
    }

    #[test]
    fn limits_bound_input_and_atom_lines() {
        let options = XyzParseOptions::default().with_max_atom_records(1);
        let error = parse_xyz_str("2\n\nHe 0 0 0\nHe 3 0 0\n", options).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Limit);
        assert_eq!(error.limit_exceeded().unwrap().resource(), Resource::Atoms);
        // A huge declared count fails on the limit before allocating.
        let error =
            parse_xyz_str("18446744073709551615\n\n", XyzParseOptions::default()).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Limit);
        let options = XyzParseOptions::default().with_max_input_bytes(3);
        let error = parse_xyz_str("1\n\nHe 0 0 0\n", options).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Limit);
    }
}
