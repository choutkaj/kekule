//! XYZ files as typed frames.
//!
//! An XYZ file is a sequence of frames. Each frame is an atom-count line, a
//! free-text comment line, and one line per atom: an element token, three
//! Cartesian coordinates, and optionally further whitespace-separated columns,
//! as extended XYZ writes them. Parsing types the count, element token, and
//! coordinates and keeps the comment and any further columns verbatim;
//! interpretation decides what they mean. [`XyzFrames`] parses one frame at
//! a time, so a reader of a long file can stop after the frames it needs.
//!
//! A host that holds atoms as elements and positions, such as one reading a
//! format that lists a molecule beside other data, builds the same document
//! with [`XyzAtom::new`], [`XyzFrame::new`], and [`XyzDocument::new`] and
//! interprets it as it would a parsed file.
use std::io::{BufRead, BufReader, Read};

use crate::diagnostics::{
    read_bounded, DiagnosticKind, Format, FormatError, LimitExceeded, Resource, SourceLocation,
};
use crate::io::lines::read_line;
use crate::math::Point3;
use crate::molecule::Element;

use super::atom_line::{parse_atom_line, AtomLineError};
use super::interpret::{interpret_xyz, interpret_xyz_structure_set};
use super::{
    XyzInterpretOptions, XyzInterpretation, XyzStructureSetInterpretOptions,
    XyzStructureSetInterpretation,
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
    parsed: bool,
}

impl XyzAtom {
    /// An atom of `element` at `position`, in the length unit interpretation
    /// is given, as a host that holds atoms rather than XYZ text builds it.
    ///
    /// Its element token is the element's symbol and it has no further
    /// columns. Its [location](Self::location) is its 1-based position in
    /// the frame [`XyzFrame::new`] places it in, as line 1 for the first
    /// atom, and line 1 until then.
    ///
    /// The position is not checked here: interpretation rejects a
    /// coordinate that is not finite, as parsing does.
    pub fn new(element: Element, position: Point3) -> Self {
        Self {
            element: element.symbol().to_owned(),
            position,
            extra: None,
            line: 1,
            parsed: false,
        }
    }

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

    /// The source location of the atom line: for a parsed atom, its line in
    /// the input; for an atom from [`Self::new`], its 1-based position in its
    /// frame.
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
    parsed: bool,
}

impl XyzFrame {
    /// A frame of `atoms` in order with `comment` as its comment line, as a
    /// host that holds atoms rather than XYZ text builds it.
    ///
    /// The comment is kept verbatim, as a parsed one is, even if it is
    /// empty or spans several lines. Each atom from [`XyzAtom::new`] takes
    /// its 1-based position in `atoms` as its location; a parsed atom keeps
    /// its line. The frame's [location](Self::location) is its 1-based
    /// position in the document [`XyzDocument::new`] places it in, and line 1
    /// until then.
    ///
    /// A frame may be built without atoms, but interpretation rejects it, as
    /// parsing rejects a count of zero.
    pub fn new(comment: impl Into<String>, atoms: Vec<XyzAtom>) -> Self {
        let mut atoms = atoms;
        for (index, atom) in atoms.iter_mut().enumerate() {
            if !atom.parsed {
                atom.line = index + 1;
            }
        }
        Self {
            comment: comment.into(),
            atoms,
            line: 1,
            parsed: false,
        }
    }

    /// The comment line verbatim, without its line terminator.
    pub fn comment(&self) -> &str {
        &self.comment
    }

    /// The atom lines in source order.
    pub fn atoms(&self) -> &[XyzAtom] {
        &self.atoms
    }

    /// The source location of the frame's atom-count line: for a parsed
    /// frame, its line in the input; for a frame from [`Self::new`], its
    /// 1-based position in its document.
    pub const fn location(&self) -> SourceLocation {
        SourceLocation::at_line(self.line)
    }
}

/// A parsed XYZ file, or a document a host built: its frames in source
/// order.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzDocument {
    frames: Vec<XyzFrame>,
}

impl XyzDocument {
    /// A document of `frames` in order, as a host that holds atoms rather
    /// than XYZ text builds it, to interpret with [`Self::interpret`] or
    /// [`Self::interpret_structure_set`] as a parsed file.
    ///
    /// Each frame from [`XyzFrame::new`] takes its 1-based position in
    /// `frames` as its location; a parsed frame keeps its line. A document
    /// without frames is valid, as empty input is, and interpretation
    /// rejects it.
    ///
    /// ```
    /// use kekule::math::Point3;
    /// use kekule::molecule::Element;
    /// use kekule::xyz::{XyzAtom, XyzDocument, XyzFrame};
    ///
    /// let atoms = vec![
    ///     XyzAtom::new(Element::O, Point3::new(0.0, 0.0, 0.0)),
    ///     XyzAtom::new(Element::H, Point3::new(0.957, 0.0, 0.0)),
    ///     XyzAtom::new(Element::H, Point3::new(-0.240, 0.927, 0.0)),
    /// ];
    /// let document = XyzDocument::new(vec![XyzFrame::new("water", atoms)]);
    /// let interpretation = document.interpret(Default::default())?;
    /// assert_eq!(interpretation.report().inferred_bonds().len(), 2);
    /// assert_eq!(document.frames()[0].atoms()[2].location().line(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(frames: Vec<XyzFrame>) -> Self {
        let mut frames = frames;
        for (index, frame) in frames.iter_mut().enumerate() {
            if !frame.parsed {
                frame.line = index + 1;
            }
        }
        Self { frames }
    }

    /// The frames in source order.
    pub fn frames(&self) -> &[XyzFrame] {
        &self.frames
    }

    /// Interprets one frame as a [`Structure`](crate::structure::Structure) with
    /// connectivity inferred from interatomic distances.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of format
    /// [`Format::Xyz`]: `InvalidValue` for an unknown element, a frame
    /// without atoms, or a coordinate that is not finite or lies beyond
    /// 1e15 Å, `Unsupported`
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

    /// Interprets several frames over the same atoms as an unordered
    /// [`StructureSet`](crate::structure::StructureSet), one geometry per
    /// selected frame with one report each: an XYZ file states no weights or
    /// order. Each frame infers its own
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
    pub fn interpret_structure_set(
        &self,
        options: XyzStructureSetInterpretOptions,
    ) -> Result<XyzStructureSetInterpretation, FormatError> {
        interpret_xyz_structure_set(self, &options)
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
        return Err(input_limit(1, options.max_input_bytes, input.len()));
    }
    // Every terminator (`\r\n`, `\r`, `\n`) ends one line; a final terminator
    // does not start another.
    let body = input.strip_suffix('\n').unwrap_or(input);
    let body = body.strip_suffix('\r').unwrap_or(body);
    let lines = (!input.is_empty())
        .then(|| {
            body.split('\n')
                .flat_map(|chunk| chunk.strip_suffix('\r').unwrap_or(chunk).split('\r'))
        })
        .into_iter()
        .flatten()
        .enumerate()
        .map(|(index, line)| (index + 1, line));
    let mut lines = TextLines(lines, std::marker::PhantomData);
    let mut frames = Vec::new();
    let mut atom_records = 0usize;
    while let Some(frame) = next_frame(&mut lines, &options, &mut atom_records)? {
        frames.push(frame);
    }
    Ok(XyzDocument { frames })
}

/// Streams the frames of an XYZ input one at a time.
///
/// Each frame's atom-count line says where the next frame begins, so memory
/// is bounded by one frame, and the first malformed frame ends the stream:
/// no count line follows it to resume from. Blank lines may follow the last
/// frame; any other content after them is a syntax error. The limits of
/// [`XyzParseOptions`] apply to the input read so far.
#[derive(Debug)]
pub struct XyzFrames<R> {
    lines: ReadLines<R>,
    options: XyzParseOptions,
    atom_records: usize,
    done: bool,
}

impl<R: Read> XyzFrames<R> {
    /// Streams the frames of `reader` under the limits of `options`.
    pub fn new(reader: R, options: XyzParseOptions) -> Self {
        Self {
            lines: ReadLines {
                reader: Counted {
                    inner: BufReader::new(reader),
                    consumed: 0,
                },
                line: Vec::new(),
                number: 0,
                max_bytes: options.max_input_bytes,
            },
            options,
            atom_records: 0,
            done: false,
        }
    }
}

impl<R: Read> Iterator for XyzFrames<R> {
    type Item = Result<XyzFrame, FormatError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let frame = next_frame(&mut self.lines, &self.options, &mut self.atom_records);
        self.done = !matches!(frame, Ok(Some(_)));
        frame.transpose()
    }
}

impl<R: Read> std::iter::FusedIterator for XyzFrames<R> {}

/// The lines of an XYZ input, numbered from 1, without their terminators.
trait Lines {
    fn next_line(&mut self) -> Result<Option<(usize, &str)>, FormatError>;
}

/// The lines of an input held whole.
struct TextLines<'a, I>(I, std::marker::PhantomData<&'a str>);

impl<'a, I: Iterator<Item = (usize, &'a str)>> Lines for TextLines<'a, I> {
    fn next_line(&mut self) -> Result<Option<(usize, &str)>, FormatError> {
        Ok(self.0.next())
    }
}

/// The lines of an input read as they are needed.
#[derive(Debug)]
struct ReadLines<R> {
    reader: Counted<BufReader<R>>,
    line: Vec<u8>,
    number: usize,
    max_bytes: usize,
}

impl<R: Read> Lines for ReadLines<R> {
    fn next_line(&mut self) -> Result<Option<(usize, &str)>, FormatError> {
        // A line is read no further than the bytes left and one more, which
        // tells that the input exceeds them, so an endless line ends too.
        let room = self.max_bytes.saturating_sub(self.reader.consumed);
        let bound = u64::try_from(room).map_or(u64::MAX, |room| room.saturating_add(1));
        let mut bounded = (&mut self.reader).take(bound);
        let read = read_line(&mut bounded, &mut self.line, room).map_err(|error| {
            xyz_parse_error(DiagnosticKind::Io, self.number + 1, "").with_source(error)
        })?;
        if read.is_none() {
            return Ok(None);
        }
        self.number += 1;
        if self.reader.consumed > self.max_bytes {
            return Err(input_limit(
                self.number,
                self.max_bytes,
                self.reader.consumed,
            ));
        }
        let text = std::str::from_utf8(&self.line).map_err(|error| {
            xyz_parse_error(
                DiagnosticKind::Io,
                self.number,
                "The input is not valid UTF-8.",
            )
            .with_source(error)
        })?;
        Ok(Some((self.number, text)))
    }
}

/// A buffered reader that counts the bytes consumed from it, terminators
/// included.
#[derive(Debug)]
struct Counted<B> {
    inner: B,
    consumed: usize,
}

impl<B: BufRead> Read for Counted<B> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.consumed = self.consumed.saturating_add(read);
        Ok(read)
    }
}

impl<B: BufRead> BufRead for Counted<B> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amount: usize) {
        self.inner.consume(amount);
        self.consumed = self.consumed.saturating_add(amount);
    }
}

fn input_limit(line: usize, max_bytes: usize, observed: usize) -> FormatError {
    let limit = LimitExceeded::new(Resource::InputBytes, max_bytes, observed);
    xyz_parse_error(
        DiagnosticKind::Limit,
        line,
        format!("The input exceeds its {}.", limit.describe()),
    )
    .exceeding(limit)
}

/// The next frame of `lines`, `atom_records` counting the atom lines of the
/// frames before; `None` after the last.
fn next_frame(
    lines: &mut impl Lines,
    options: &XyzParseOptions,
    atom_records: &mut usize,
) -> Result<Option<XyzFrame>, FormatError> {
    let (count_line, blank, count) = match lines.next_line()? {
        None => return Ok(None),
        Some((line, text)) => (line, text.trim().is_empty(), text.trim().parse::<usize>()),
    };
    if blank {
        // Blank lines may only trail the last frame.
        while let Some((line, text)) = lines.next_line()? {
            if !text.trim().is_empty() {
                return Err(xyz_parse_error(
                    DiagnosticKind::Syntax,
                    line,
                    "Content follows a blank line after the last frame.",
                ));
            }
        }
        return Ok(None);
    }
    let count = count.map_err(|_| {
        xyz_parse_error(
            DiagnosticKind::Syntax,
            count_line,
            "The atom count line does not hold one whole number.",
        )
    })?;
    if count == 0 {
        return Err(xyz_parse_error(
            DiagnosticKind::InvalidValue,
            count_line,
            "The frame has zero atoms.",
        ));
    }
    let total = atom_records.saturating_add(count);
    if total > options.max_atom_records {
        let limit = LimitExceeded::new(Resource::Atoms, options.max_atom_records, total);
        return Err(xyz_parse_error(
            DiagnosticKind::Limit,
            count_line,
            format!("The input exceeds its {}.", limit.describe()),
        )
        .exceeding(limit));
    }
    *atom_records = total;
    let Some((_, comment)) = lines.next_line()? else {
        return Err(xyz_parse_error(
            DiagnosticKind::Syntax,
            count_line,
            "The frame ends before its comment line.",
        ));
    };
    let comment = comment.to_owned();
    let mut atoms = Vec::with_capacity(count);
    for atom_index in 0..count {
        let Some((line, text)) = lines.next_line()? else {
            return Err(xyz_parse_error(
                DiagnosticKind::Syntax,
                count_line,
                format!(
                    "The frame declares {} but ends after {atom_index}.",
                    crate::diagnostics::counted(count, "atom", "atoms")
                ),
            ));
        };
        atoms.push(parse_atom(line, text)?);
    }
    Ok(Some(XyzFrame {
        comment,
        atoms,
        line: count_line,
        parsed: true,
    }))
}

fn parse_atom(line: usize, text: &str) -> Result<XyzAtom, FormatError> {
    let atom = parse_atom_line(text).map_err(|error| match error {
        AtomLineError::MissingField => xyz_parse_error(
            DiagnosticKind::Syntax,
            line,
            "The atom line lacks an element or one of three coordinates.",
        ),
        AtomLineError::NotNumber(value) | AtomLineError::NotFinite(value) => xyz_parse_error(
            DiagnosticKind::InvalidValue,
            line,
            format!("The coordinate {value:?} is not a finite number."),
        ),
    })?;
    Ok(XyzAtom {
        element: atom.element.to_owned(),
        position: atom.position,
        extra: atom.extra.map(str::to_owned),
        line,
        parsed: true,
    })
}

#[cfg(test)]
mod tests;
