//! Typed Molfile documents: what a V2000 or V3000 connection table states,
//! decoded from its fixed-width fields, codes, and options, before chemistry.

use std::fmt;

use crate::diagnostics::{
    normalize_line_endings, read_bounded, DiagnosticKind, Format, FormatError, LimitExceeded,
    Resource, SourceLocation,
};
use crate::geometry::Point3;

use super::structure_documents::{
    interpret_molfile_document, MolfileInterpretOptions, MolfileInterpretation,
};
use super::v2000::parse_v2000_ctab;
use super::v3000::parse_v3000_ctab;

/// Resource limits shared by version-autodetected Molfile parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolfileParseOptions {
    max_input_bytes: usize,
    max_v3000_atoms: usize,
    max_v3000_bonds: usize,
    max_v3000_logical_line_bytes: usize,
}

impl MolfileParseOptions {
    /// Maximum UTF-8 byte length of the complete Molfile document.
    pub const fn max_input_bytes(&self) -> usize {
        self.max_input_bytes
    }

    /// Sets [`Self::max_input_bytes`].
    #[must_use]
    pub const fn with_max_input_bytes(mut self, max_input_bytes: usize) -> Self {
        self.max_input_bytes = max_input_bytes;
        self
    }

    /// Maximum declared V3000 atom count.
    pub const fn max_v3000_atoms(&self) -> usize {
        self.max_v3000_atoms
    }

    /// Sets [`Self::max_v3000_atoms`].
    #[must_use]
    pub const fn with_max_v3000_atoms(mut self, max_v3000_atoms: usize) -> Self {
        self.max_v3000_atoms = max_v3000_atoms;
        self
    }

    /// Maximum declared V3000 bond count.
    pub const fn max_v3000_bonds(&self) -> usize {
        self.max_v3000_bonds
    }

    /// Sets [`Self::max_v3000_bonds`].
    #[must_use]
    pub const fn with_max_v3000_bonds(mut self, max_v3000_bonds: usize) -> Self {
        self.max_v3000_bonds = max_v3000_bonds;
        self
    }

    /// Maximum byte length after joining one continued `M  V30` logical line.
    pub const fn max_v3000_logical_line_bytes(&self) -> usize {
        self.max_v3000_logical_line_bytes
    }

    /// Sets [`Self::max_v3000_logical_line_bytes`].
    #[must_use]
    pub const fn with_max_v3000_logical_line_bytes(
        mut self,
        max_v3000_logical_line_bytes: usize,
    ) -> Self {
        self.max_v3000_logical_line_bytes = max_v3000_logical_line_bytes;
        self
    }
}

impl MolfileParseOptions {
    /// The [`Default`] bound of [`Self::max_input_bytes`].
    pub const DEFAULT_MAX_INPUT_BYTES: usize = 64 * 1024 * 1024;
    /// The [`Default`] bound of [`Self::max_v3000_atoms`].
    pub const DEFAULT_MAX_V3000_ATOMS: usize = 1_000_000;
    /// The [`Default`] bound of [`Self::max_v3000_bonds`].
    pub const DEFAULT_MAX_V3000_BONDS: usize = 2_000_000;
    /// The [`Default`] bound of [`Self::max_v3000_logical_line_bytes`].
    pub const DEFAULT_MAX_V3000_LOGICAL_LINE_BYTES: usize = 1024 * 1024;
}

impl Default for MolfileParseOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: Self::DEFAULT_MAX_INPUT_BYTES,
            max_v3000_atoms: Self::DEFAULT_MAX_V3000_ATOMS,
            max_v3000_bonds: Self::DEFAULT_MAX_V3000_BONDS,
            max_v3000_logical_line_bytes: Self::DEFAULT_MAX_V3000_LOGICAL_LINE_BYTES,
        }
    }
}

/// Connection-table version declared on the counts line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MolfileVersion {
    /// Fixed-width V2000 connection table.
    V2000,
    /// Extended `M  V30` connection table.
    V3000,
}

impl MolfileVersion {
    pub(super) const fn name(self) -> &'static str {
        match self {
            Self::V2000 => "V2000",
            Self::V3000 => "V3000",
        }
    }
}

/// The three header lines of a Molfile, kept verbatim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileHeader {
    pub(super) title: String,
    pub(super) program: String,
    pub(super) comment: String,
}

impl MolfileHeader {
    /// Header line 1, the molecule name; may be empty.
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Header line 2, the program/user/date line; may be empty.
    pub fn program(&self) -> &str {
        &self.program
    }

    /// Header line 3, the comment line; may be empty.
    pub fn comment(&self) -> &str {
        &self.comment
    }
}

/// Radical state declared by a V2000 `M  RAD` record or charge code 4, or a
/// V3000 `RAD` option.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MolfileRadical {
    /// Radical code 1: a singlet divalent radical (carbene).
    Singlet,
    /// Radical code 2, or V2000 charge code 4: a monovalent doublet radical.
    Doublet,
    /// Radical code 3: a triplet divalent radical.
    Triplet,
}

/// Atom parity declared by a V3000 atom `CFG` option.
///
/// V2000 atom-block parity is ignored when read, as CTfile specifies, so V2000
/// atoms never carry one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MolfileAtomParity {
    /// `CFG=1`.
    Odd,
    /// `CFG=2`.
    Even,
    /// `CFG=3`: the configuration is explicitly unknown.
    Either,
}

/// One atom of the connection table with its decoded fields and properties.
#[derive(Debug, Clone, PartialEq)]
pub struct MolfileAtom {
    pub(super) symbol: String,
    pub(super) position: Point3,
    pub(super) formal_charge: i32,
    pub(super) isotope: Option<u16>,
    pub(super) mass_difference: Option<i8>,
    pub(super) radical: Option<MolfileRadical>,
    pub(super) valence: Option<u8>,
    pub(super) atom_map: Option<u32>,
    pub(super) parity: Option<MolfileAtomParity>,
    pub(super) location: SourceLocation,
}

impl MolfileAtom {
    pub(super) fn new(symbol: &str, position: Point3, location: SourceLocation) -> Self {
        Self {
            symbol: symbol.to_owned(),
            position,
            formal_charge: 0,
            isotope: None,
            mass_difference: None,
            radical: None,
            valence: None,
            atom_map: None,
            parity: None,
            location,
        }
    }

    /// The atom symbol as written.
    pub fn symbol(&self) -> &str {
        &self.symbol
    }

    /// Coordinates in angstroms.
    pub const fn position(&self) -> Point3 {
        self.position
    }

    /// Formal charge in elementary charges (e), from the V2000 charge code,
    /// an `M  CHG` record, or V3000 `CHG`; zero when none is declared.
    pub const fn formal_charge(&self) -> i32 {
        self.formal_charge
    }

    /// Declared mass number.
    pub const fn isotope(&self) -> Option<u16> {
        self.isotope
    }

    /// A nonzero legacy V2000 mass difference from the element's nominal
    /// mass, which interpretation resolves to a mass number. An `M  ISO`
    /// record replaces every atom-block mass difference.
    pub const fn mass_difference(&self) -> Option<i8> {
        self.mass_difference
    }

    /// Declared radical state; `None` when none is declared.
    pub const fn radical(&self) -> Option<MolfileRadical> {
        self.radical
    }

    /// Declared total valence: the V2000 valence field or V3000 `VAL`.
    /// `Some(0)` is an explicit zero valence.
    pub const fn valence(&self) -> Option<u8> {
        self.valence
    }

    /// Atom-atom mapping number; `None` when the field is absent or zero.
    pub const fn atom_map(&self) -> Option<u32> {
        self.atom_map
    }

    /// Declared V3000 `CFG` parity; always `None` for V2000 atoms.
    pub const fn parity(&self) -> Option<MolfileAtomParity> {
        self.parity
    }

    /// The atom line; for a continued V3000 line, its first physical line.
    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// Decoded bond type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MolfileBondOrder {
    /// Zero-order bond (type 0).
    Zero,
    /// Single bond (type 1).
    Single,
    /// Double bond (type 2).
    Double,
    /// Triple bond (type 3).
    Triple,
    /// Aromatic bond (type 4).
    Aromatic,
    /// Coordination bond (type 9).
    Dative,
}

/// Decoded bond stereo mark, drawn from the bond's first atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MolfileBondStereo {
    /// Wedge toward the viewer on a single bond.
    Up,
    /// Hashed wedge away from the viewer on a single bond.
    Down,
    /// Wavy single bond, or crossed double bond: the configuration is unknown.
    Either,
}

/// One bond of the connection table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileBond {
    pub(super) atoms: [usize; 2],
    pub(super) order: MolfileBondOrder,
    pub(super) stereo: Option<MolfileBondStereo>,
    pub(super) location: SourceLocation,
}

impl MolfileBond {
    /// Endpoints as indices into [`MolfileDocument::atoms`], in source order.
    pub const fn atoms(&self) -> [usize; 2] {
        self.atoms
    }

    /// Decoded bond type.
    pub const fn order(&self) -> MolfileBondOrder {
        self.order
    }

    /// Stereo mark; `None` when the bond carries none.
    pub const fn stereo(&self) -> Option<MolfileBondStereo> {
        self.stereo
    }

    /// The bond line; for a continued V3000 line, its first physical line.
    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// Kind and number of a V3000 enhanced stereo group.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MolfileStereoGroupKind {
    /// `MDLV30/STEABS`.
    Absolute,
    /// `MDLV30/STERELn`, carrying the positive group number `n`.
    Or(usize),
    /// `MDLV30/STERACn`, carrying the positive group number `n`.
    And(usize),
}

/// One V3000 enhanced stereo group. Collection lines naming the same group
/// are merged.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolfileStereoGroup {
    pub(super) kind: MolfileStereoGroupKind,
    pub(super) atoms: Vec<usize>,
    pub(super) location: SourceLocation,
}

impl MolfileStereoGroup {
    /// Kind and number of the group.
    pub const fn kind(&self) -> MolfileStereoGroupKind {
        self.kind
    }

    /// Members as indices into [`MolfileDocument::atoms`].
    pub fn atoms(&self) -> &[usize] {
        &self.atoms
    }

    /// The first collection line of the group.
    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// Valid source content that Molfile interpretation does not support.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MolfileUnsupportedContent {
    /// A V2000 property record other than `M  CHG`, `M  ISO`, and `M  RAD`,
    /// identified by its code, such as `STY`.
    PropertyRecord(String),
    /// Any other record identified by its leading keyword: a V2000 atom alias
    /// (`A`), atom value (`V`), or group abbreviation (`G`), or a V3000 CTAB
    /// record such as `LINKNODE`.
    Record(String),
    /// A V3000 atom option other than `CFG`, `CHG`, `MASS`, `RAD`, and `VAL`,
    /// including a nonzero `HCOUNT`, which is a query constraint. `atom`
    /// indexes [`MolfileDocument::atoms`].
    AtomOption {
        /// Index into [`MolfileDocument::atoms`].
        atom: usize,
        /// The option keyword as written, such as `HCOUNT`.
        key: String,
    },
    /// A V3000 bond option other than `CFG`. `bond` indexes
    /// [`MolfileDocument::bonds`].
    BondOption {
        /// Index into [`MolfileDocument::bonds`].
        bond: usize,
        /// The option keyword as written, such as `TOPO`.
        key: String,
    },
    /// A nonzero V2000 atom-block query or reaction field: hydrogen count,
    /// stereo care, H0 designator, inversion/retention or exact change. `atom` indexes
    /// [`MolfileDocument::atoms`].
    AtomField {
        /// Index into [`MolfileDocument::atoms`].
        atom: usize,
        /// Human-readable name of the field.
        field: String,
    },
    /// A nonzero V2000 bond-block ring/chain topology or reacting-center
    /// field. `bond` indexes [`MolfileDocument::bonds`].
    BondField {
        /// Index into [`MolfileDocument::bonds`].
        bond: usize,
        /// Human-readable name of the field.
        field: String,
    },
    /// A V3000 collection other than enhanced stereo groups, such as
    /// `MDLV30/HILITE`.
    Collection(String),
    /// A V3000 block such as `SGROUP` or `OBJ3D`.
    Block(String),
    /// Nonblank text after `M  END` of a standalone Molfile.
    TrailingContent,
}

impl fmt::Display for MolfileUnsupportedContent {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PropertyRecord(code) => write!(formatter, "property record `M  {code}`"),
            Self::Record(keyword) => write!(formatter, "record `{keyword}`"),
            Self::AtomOption { key, .. } => write!(formatter, "V3000 atom option `{key}`"),
            Self::BondOption { key, .. } => write!(formatter, "V3000 bond option `{key}`"),
            Self::AtomField { field, .. } => write!(formatter, "V2000 atom {field} field"),
            Self::BondField { field, .. } => write!(formatter, "V2000 bond {field} field"),
            Self::Collection(name) => write!(formatter, "V3000 collection `{name}`"),
            Self::Block(name) => write!(formatter, "V3000 block `{name}`"),
            Self::TrailingContent => formatter.write_str("content after `M  END`"),
        }
    }
}

/// Unsupported content and where it occurs.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MolfileUnsupported {
    pub(super) content: MolfileUnsupportedContent,
    pub(super) location: SourceLocation,
}

impl MolfileUnsupported {
    /// Pairs unsupported content with the source location it occurs at.
    pub const fn new(content: MolfileUnsupportedContent, location: SourceLocation) -> Self {
        Self { content, location }
    }

    /// What is unsupported.
    pub fn content(&self) -> &MolfileUnsupportedContent {
        &self.content
    }

    /// The line where the content occurs.
    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// The decoded connection table of one Molfile.
#[derive(Debug, Clone, Default, PartialEq)]
pub(super) struct MolfileCtab {
    pub(super) chiral_flag: bool,
    pub(super) atoms: Vec<MolfileAtom>,
    pub(super) bonds: Vec<MolfileBond>,
    pub(super) stereo_groups: Vec<MolfileStereoGroup>,
    pub(super) unsupported: Vec<MolfileUnsupported>,
}

/// One parsed V2000 or V3000 Molfile: what the source states, before chemistry.
///
/// Locations refer to the original input, so a Molfile inside an SDF record
/// reports lines of the SDF file.
#[derive(Debug, Clone, PartialEq)]
pub struct MolfileDocument {
    pub(super) version: MolfileVersion,
    pub(super) header: MolfileHeader,
    pub(super) ctab: MolfileCtab,
}

impl MolfileDocument {
    /// Version declared on the counts line.
    pub const fn version(&self) -> MolfileVersion {
        self.version
    }

    /// The three header lines.
    pub fn header(&self) -> &MolfileHeader {
        &self.header
    }

    /// The counts-line chiral flag. It is retained as source metadata;
    /// interpretation does not translate it into stereo groups.
    pub const fn chiral_flag(&self) -> bool {
        self.ctab.chiral_flag
    }

    /// Atoms in atom-block order.
    pub fn atoms(&self) -> &[MolfileAtom] {
        &self.ctab.atoms
    }

    /// Bonds in bond-block order.
    pub fn bonds(&self) -> &[MolfileBond] {
        &self.ctab.bonds
    }

    /// V3000 enhanced stereo groups in order of first appearance; always
    /// empty for V2000.
    pub fn stereo_groups(&self) -> &[MolfileStereoGroup] {
        &self.ctab.stereo_groups
    }

    /// Content that interpretation rejects by default, in source order.
    pub fn unsupported(&self) -> &[MolfileUnsupported] {
        &self.ctab.unsupported
    }

    /// Interprets this document once into connected chemistry with matching
    /// component geometry and source reports. Unsupported content fails or is
    /// omitted and reported, as `options` selects.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] located at the offending atom or record, with kind
    /// [`DiagnosticKind::Unsupported`] for content in
    /// [`Self::unsupported`] under the default rejecting policy or valences
    /// Kekule cannot resolve, [`DiagnosticKind::InvalidValue`] for an unknown
    /// element symbol or an out-of-range charge or mass,
    /// [`DiagnosticKind::Inconsistent`] for a declared valence below the
    /// bonded valence, and [`DiagnosticKind::Chemistry`] when the represented
    /// chemistry cannot be published.
    pub fn interpret(
        &self,
        options: MolfileInterpretOptions,
    ) -> Result<MolfileInterpretation, FormatError> {
        interpret_molfile_document(self, options)
    }
}

/// A Molfile parse failure on `line` of the input.
pub(super) fn molfile_parse_error(
    kind: DiagnosticKind,
    line: usize,
    message: impl Into<String>,
) -> FormatError {
    FormatError::parse(Format::Molfile, kind, message).at(SourceLocation::at_line(line))
}

/// The lines of one Molfile, numbered as in the original input.
#[derive(Debug, Clone, Copy)]
pub(super) struct MolfileLines<'a> {
    lines: &'a [&'a str],
    first_line: usize,
}

impl<'a> MolfileLines<'a> {
    pub(super) const fn new(lines: &'a [&'a str], first_line: usize) -> Self {
        Self { lines, first_line }
    }

    pub(super) const fn len(&self) -> usize {
        self.lines.len()
    }

    pub(super) fn text(&self, index: usize) -> Option<&'a str> {
        self.lines.get(index).copied()
    }

    /// The original line number of `index`, which may be one past the end.
    pub(super) const fn number(&self, index: usize) -> usize {
        self.first_line.saturating_add(index)
    }

    pub(super) const fn location(&self, index: usize) -> SourceLocation {
        SourceLocation::at_line(self.number(index))
    }

    pub(super) fn error(
        &self,
        index: usize,
        kind: DiagnosticKind,
        message: impl Into<String>,
    ) -> FormatError {
        molfile_parse_error(kind, self.number(index), message)
    }
}

pub(crate) fn parse_molfile_document(
    input: &str,
    options: MolfileParseOptions,
) -> Result<MolfileDocument, FormatError> {
    if input.len() > options.max_input_bytes {
        return Err(molfile_parse_error(
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
    let source = normalize_line_endings(input);
    let lines = source.lines().collect::<Vec<_>>();
    let lines = MolfileLines::new(&lines, 1);
    let (mut document, end) = parse_molfile_lines(lines, options)?;
    document.ctab.unsupported.extend(
        (end + 1..lines.len())
            .filter(|&index| {
                lines
                    .text(index)
                    .is_some_and(|line| !line.trim().is_empty())
            })
            .map(|index| MolfileUnsupported {
                content: MolfileUnsupportedContent::TrailingContent,
                location: lines.location(index),
            }),
    );
    Ok(document)
}

/// Reads and parses one Molfile, enforcing the byte limit while reading.
pub(crate) fn parse_molfile_reader(
    reader: impl std::io::Read,
    options: MolfileParseOptions,
) -> Result<MolfileDocument, FormatError> {
    let input = read_bounded(reader, options.max_input_bytes, Format::Molfile)
        .map_err(|error| error.at(SourceLocation::at_line(1)))?;
    parse_molfile_document(&input, options)
}

/// Parses the Molfile at the start of `lines` and returns the index of its
/// `M  END` line. Content after `M  END` is the caller's.
pub(super) fn parse_molfile_lines(
    lines: MolfileLines<'_>,
    options: MolfileParseOptions,
) -> Result<(MolfileDocument, usize), FormatError> {
    let (Some(title), Some(program), Some(comment), Some(counts)) =
        (lines.text(0), lines.text(1), lines.text(2), lines.text(3))
    else {
        return Err(lines.error(
            0,
            DiagnosticKind::Syntax,
            "Molfile must contain three header lines and a counts line",
        ));
    };
    let version = if counts.contains("V3000") {
        MolfileVersion::V3000
    } else if counts.contains("V2000") {
        MolfileVersion::V2000
    } else {
        return Err(lines.error(
            3,
            DiagnosticKind::Syntax,
            "counts line must declare V2000 or V3000",
        ));
    };
    let (ctab, end) = match version {
        MolfileVersion::V2000 => parse_v2000_ctab(lines)?,
        MolfileVersion::V3000 => parse_v3000_ctab(lines, options)?,
    };
    let header = MolfileHeader {
        title: title.to_owned(),
        program: program.to_owned(),
        comment: comment.to_owned(),
    };
    Ok((
        MolfileDocument {
            version,
            header,
            ctab,
        },
        end,
    ))
}
