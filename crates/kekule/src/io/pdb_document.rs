//! PDB-format coordinate entries as typed records.
//!
//! A PDB file holds one entry. Its coordinate (`ATOM`, `HETATM`), chain
//! termination (`TER`), model (`MODEL`, `ENDMDL`), and connection (`CONECT`,
//! `LINK`, `SSBOND`) records are typed by their fixed columns; every other
//! record is kept verbatim. Numbers wider than their columns follow the
//! hybrid-36 convention.
use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use crate::core::Molecule;
use crate::diagnostics::{read_bounded, DiagnosticKind, SourceLocation};
use crate::geometry::Point3;
use crate::structure::Model;
use crate::topology::Topology;

use super::structure_interpret::{
    interpret_pdb, interpret_pdb_conformations, interpret_pdb_ensemble, pdb_alternate_locations,
};
use super::{
    AltLocResidue, StructureEnsembleInterpretError, StructureEnsembleInterpretOptions,
    StructureEnsembleInterpretation, StructureInterpretError, StructureInterpretOptions,
    StructureInterpretation,
};

/// Syntax and resource limits for parsing a PDB entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PdbParseOptions {
    /// Maximum UTF-8 byte length of the input.
    pub max_input_bytes: usize,
    /// Maximum number of `ATOM` and `HETATM` records.
    pub max_atom_records: usize,
}

impl Default for PdbParseOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: 1024 * 1024 * 1024,
            max_atom_records: 5_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbParseError {
    kind: DiagnosticKind,
    location: SourceLocation,
    message: String,
}

impl PdbParseError {
    fn new(kind: DiagnosticKind, location: SourceLocation, message: impl Into<String>) -> Self {
        Self {
            kind,
            location,
            message: message.into(),
        }
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

impl fmt::Display for PdbParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PDB parse error at {}: {}", self.location, self.message)
    }
}

impl std::error::Error for PdbParseError {}

/// Whether a coordinate record is an `ATOM` or a `HETATM` record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PdbAtomKind {
    Atom,
    Hetatm,
}

/// One `ATOM` or `HETATM` record.
#[derive(Debug, Clone, PartialEq)]
pub struct PdbAtom {
    kind: PdbAtomKind,
    serial: Option<u32>,
    /// Columns 13-16 as written; their alignment encodes the element.
    name_field: String,
    alternate_location: Option<char>,
    residue_name: String,
    chain_id: Option<char>,
    residue_sequence: i32,
    insertion_code: Option<char>,
    position: Point3,
    occupancy: Option<f64>,
    temperature_factor: Option<f64>,
    segment_id: Option<String>,
    element: Option<String>,
    formal_charge: Option<i8>,
    model: Option<u32>,
    location: SourceLocation,
}

impl PdbAtom {
    pub const fn kind(&self) -> PdbAtomKind {
        self.kind
    }

    /// The atom serial number, or `None` when blank or overflowed (`*****`).
    pub const fn serial(&self) -> Option<u32> {
        self.serial
    }

    /// The atom name, without the alignment blanks of its columns.
    pub fn name(&self) -> &str {
        self.name_field.trim()
    }

    /// Columns 13-16 as written. By convention two-letter element symbols
    /// start in column 13 and one-letter symbols in column 14.
    pub fn name_field(&self) -> &str {
        &self.name_field
    }

    pub const fn alternate_location(&self) -> Option<char> {
        self.alternate_location
    }

    pub fn residue_name(&self) -> &str {
        &self.residue_name
    }

    pub const fn chain_id(&self) -> Option<char> {
        self.chain_id
    }

    pub const fn residue_sequence(&self) -> i32 {
        self.residue_sequence
    }

    pub const fn insertion_code(&self) -> Option<char> {
        self.insertion_code
    }

    /// Orthogonal coordinates in ångström.
    pub const fn position(&self) -> Point3 {
        self.position
    }

    pub const fn occupancy(&self) -> Option<f64> {
        self.occupancy
    }

    /// The isotropic temperature factor, in square ångström.
    pub const fn temperature_factor(&self) -> Option<f64> {
        self.temperature_factor
    }

    pub fn segment_id(&self) -> Option<&str> {
        self.segment_id.as_deref()
    }

    /// The element symbol as written, if any.
    pub fn element(&self) -> Option<&str> {
        self.element.as_deref()
    }

    pub const fn formal_charge(&self) -> Option<i8> {
        self.formal_charge
    }

    /// The serial number of the enclosing `MODEL`, if the entry has models.
    pub const fn model(&self) -> Option<u32> {
        self.model
    }

    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// A `TER` record. It ends the polymer chain of the coordinate record
/// before it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbTer {
    serial: Option<u32>,
    residue_name: Option<String>,
    chain_id: Option<char>,
    residue_sequence: Option<i32>,
    insertion_code: Option<char>,
    preceding_atoms: usize,
    model: Option<u32>,
    location: SourceLocation,
}

impl PdbTer {
    pub const fn serial(&self) -> Option<u32> {
        self.serial
    }

    pub fn residue_name(&self) -> Option<&str> {
        self.residue_name.as_deref()
    }

    pub const fn chain_id(&self) -> Option<char> {
        self.chain_id
    }

    pub const fn residue_sequence(&self) -> Option<i32> {
        self.residue_sequence
    }

    pub const fn insertion_code(&self) -> Option<char> {
        self.insertion_code
    }

    /// The number of coordinate records before this one; the last of them is
    /// the record whose chain this one terminates.
    pub const fn preceding_atoms(&self) -> usize {
        self.preceding_atoms
    }

    pub const fn model(&self) -> Option<u32> {
        self.model
    }

    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// A `MODEL` record, which opens one coordinate model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PdbModel {
    serial: u32,
    location: SourceLocation,
}

impl PdbModel {
    pub const fn serial(&self) -> u32 {
        self.serial
    }

    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// A `CONECT` record: the atoms bonded to one atom, by serial number.
///
/// The record states connectivity without bond orders. By a widespread
/// convention, a partner listed more than once has a multiple bond.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbConect {
    serial: u32,
    bonded: Vec<u32>,
    location: SourceLocation,
}

impl PdbConect {
    pub const fn serial(&self) -> u32 {
        self.serial
    }

    /// Serial numbers of the bonded atoms, in column order.
    pub fn bonded(&self) -> &[u32] {
        &self.bonded
    }

    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// One atom or residue a `LINK` or `SSBOND` record names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbPartner {
    atom_name: Option<String>,
    alternate_location: Option<char>,
    residue_name: String,
    chain_id: Option<char>,
    residue_sequence: i32,
    insertion_code: Option<char>,
    symmetry: Option<String>,
}

impl PdbPartner {
    /// The atom name; `SSBOND` records name residues only.
    pub fn atom_name(&self) -> Option<&str> {
        self.atom_name.as_deref()
    }

    pub const fn alternate_location(&self) -> Option<char> {
        self.alternate_location
    }

    pub fn residue_name(&self) -> &str {
        &self.residue_name
    }

    pub const fn chain_id(&self) -> Option<char> {
        self.chain_id
    }

    pub const fn residue_sequence(&self) -> i32 {
        self.residue_sequence
    }

    pub const fn insertion_code(&self) -> Option<char> {
        self.insertion_code
    }

    /// The crystallographic symmetry operator, such as `1555` for identity.
    pub fn symmetry(&self) -> Option<&str> {
        self.symmetry.as_deref()
    }
}

/// A `LINK` record: a covalent or coordination bond between residues.
#[derive(Debug, Clone, PartialEq)]
pub struct PdbLink {
    partners: [PdbPartner; 2],
    distance: Option<f64>,
    location: SourceLocation,
}

impl PdbLink {
    pub const fn partners(&self) -> &[PdbPartner; 2] {
        &self.partners
    }

    /// The stated link distance in ångström.
    pub const fn distance(&self) -> Option<f64> {
        self.distance
    }

    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// An `SSBOND` record: a disulfide bond between two cysteine residues.
#[derive(Debug, Clone, PartialEq)]
pub struct PdbSsbond {
    serial: Option<u32>,
    partners: [PdbPartner; 2],
    distance: Option<f64>,
    location: SourceLocation,
}

impl PdbSsbond {
    pub const fn serial(&self) -> Option<u32> {
        self.serial
    }

    pub const fn partners(&self) -> &[PdbPartner; 2] {
        &self.partners
    }

    /// The stated disulfide distance in ångström.
    pub const fn distance(&self) -> Option<f64> {
        self.distance
    }

    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// A record kept verbatim, such as `HEADER`, `REMARK`, or `SEQRES`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbRecord {
    name: String,
    text: String,
    location: SourceLocation,
}

impl PdbRecord {
    /// The record name: columns 1-6 without trailing blanks.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The whole line.
    pub fn text(&self) -> &str {
        &self.text
    }

    pub const fn location(&self) -> SourceLocation {
        self.location
    }
}

/// One parsed PDB entry: its records, typed where interpretation reads them.
#[derive(Debug, Clone, PartialEq)]
pub struct PdbDocument {
    id_code: Option<String>,
    models: Vec<PdbModel>,
    atoms: Vec<PdbAtom>,
    ter_records: Vec<PdbTer>,
    conect_records: Vec<PdbConect>,
    link_records: Vec<PdbLink>,
    ssbond_records: Vec<PdbSsbond>,
    other_records: Vec<PdbRecord>,
}

impl PdbDocument {
    /// The entry's ID code, from columns 63-66 of its `HEADER` record.
    pub fn id_code(&self) -> Option<&str> {
        self.id_code.as_deref()
    }

    /// `MODEL` records in source order; empty for a single-model entry.
    pub fn models(&self) -> &[PdbModel] {
        &self.models
    }

    /// `ATOM` and `HETATM` records in source order.
    pub fn atoms(&self) -> &[PdbAtom] {
        &self.atoms
    }

    pub fn ter_records(&self) -> &[PdbTer] {
        &self.ter_records
    }

    pub fn conect_records(&self) -> &[PdbConect] {
        &self.conect_records
    }

    pub fn link_records(&self) -> &[PdbLink] {
        &self.link_records
    }

    pub fn ssbond_records(&self) -> &[PdbSsbond] {
        &self.ssbond_records
    }

    /// Every other record, verbatim and in source order.
    pub fn other_records(&self) -> &[PdbRecord] {
        &self.other_records
    }

    /// Interprets one selected coordinate model using default policy.
    pub fn interpret(&self) -> Result<StructureInterpretation, StructureInterpretError> {
        self.interpret_with_options(StructureInterpretOptions::default())
    }

    /// Interprets one selected coordinate model using explicit policy.
    pub fn interpret_with_options(
        &self,
        options: StructureInterpretOptions,
    ) -> Result<StructureInterpretation, StructureInterpretError> {
        interpret_pdb(self, options)
    }

    pub fn to_molecules(&self) -> Result<Vec<Molecule>, StructureInterpretError> {
        Ok(self.interpret()?.into_molecules())
    }

    pub fn to_molecules_with_options(
        &self,
        options: StructureInterpretOptions,
    ) -> Result<Vec<Molecule>, StructureInterpretError> {
        Ok(self.interpret_with_options(options)?.into_molecules())
    }

    pub fn to_topology(&self) -> Result<Arc<Topology>, StructureInterpretError> {
        Ok(self.interpret()?.into_topology())
    }

    pub fn to_topology_with_options(
        &self,
        options: StructureInterpretOptions,
    ) -> Result<Arc<Topology>, StructureInterpretError> {
        Ok(self.interpret_with_options(options)?.into_topology())
    }

    pub fn to_model(&self) -> Result<Model, StructureInterpretError> {
        Ok(self.interpret()?.into_model())
    }

    pub fn to_model_with_options(
        &self,
        options: StructureInterpretOptions,
    ) -> Result<Model, StructureInterpretError> {
        Ok(self.interpret_with_options(options)?.into_model())
    }

    /// Interprets all coordinate models using default policy.
    pub fn interpret_ensemble(
        &self,
    ) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
        self.interpret_ensemble_with_options(StructureEnsembleInterpretOptions::default())
    }

    /// Interprets selected coordinate models using explicit policy.
    pub fn interpret_ensemble_with_options(
        &self,
        options: StructureEnsembleInterpretOptions,
    ) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
        interpret_pdb_ensemble(self, options)
    }

    /// Interprets explicit model and alternate-location selections, in the
    /// supplied order, as an unweighted ensemble. Atom identities, chemistry,
    /// and topology layout must agree between selections.
    pub fn interpret_conformations(
        &self,
        selections: &[StructureInterpretOptions],
    ) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
        interpret_pdb_conformations(self, selections)
    }

    /// Residues with alternate locations, by the identities that
    /// [`AltLocSelection`](super::AltLocSelection) overrides use.
    pub fn alternate_locations(&self) -> Result<Vec<AltLocResidue>, StructureInterpretError> {
        pdb_alternate_locations(self)
    }
}

pub fn parse_pdb_reader(
    reader: impl std::io::Read,
    options: PdbParseOptions,
) -> Result<PdbDocument, PdbParseError> {
    let input = read_bounded(reader, options.max_input_bytes)
        .map_err(|(kind, message)| PdbParseError::new(kind, SourceLocation::at_line(1), message))?;
    parse_pdb_str(&input, options)
}

pub fn parse_pdb_str(input: &str, options: PdbParseOptions) -> Result<PdbDocument, PdbParseError> {
    if input.len() > options.max_input_bytes {
        return Err(PdbParseError::new(
            DiagnosticKind::Limit,
            SourceLocation::at_line(1),
            "input exceeds configured byte limit",
        ));
    }
    let mut parser = Parser::default();
    // Every terminator (`\r\n`, `\r`, `\n`) ends one line.
    let lines = input
        .split('\n')
        .flat_map(|chunk| chunk.strip_suffix('\r').unwrap_or(chunk).split('\r'));
    for (index, line) in lines.enumerate() {
        parser.line(line, index + 1, options)?;
    }
    Ok(parser.finish())
}

#[derive(Default)]
struct Parser {
    document: Option<PdbDocument>,
    model: Option<u32>,
    model_serials: BTreeSet<u32>,
    /// The first coordinate record outside a model, if any.
    unmodelled_atom: Option<usize>,
    ended: bool,
}

impl Parser {
    fn document(&mut self) -> &mut PdbDocument {
        self.document.get_or_insert_with(|| PdbDocument {
            id_code: None,
            models: Vec::new(),
            atoms: Vec::new(),
            ter_records: Vec::new(),
            conect_records: Vec::new(),
            link_records: Vec::new(),
            ssbond_records: Vec::new(),
            other_records: Vec::new(),
        })
    }

    fn finish(mut self) -> PdbDocument {
        self.document();
        self.document.expect("document exists")
    }

    fn line(
        &mut self,
        line: &str,
        number: usize,
        options: PdbParseOptions,
    ) -> Result<(), PdbParseError> {
        if line.trim().is_empty() {
            return Ok(());
        }
        let at_line = SourceLocation::at_line(number);
        if self.ended {
            return Err(PdbParseError::new(
                DiagnosticKind::Unsupported,
                at_line,
                "content after the END record; split concatenated entries before parsing",
            ));
        }
        let end = (0..=line.len().min(6))
            .rev()
            .find(|&index| line.is_char_boundary(index))
            .unwrap_or(0);
        let name = line[..end].trim_end();
        let typed = matches!(
            name,
            "ATOM" | "HETATM" | "TER" | "MODEL" | "ENDMDL" | "CONECT" | "LINK" | "SSBOND"
        );
        if typed && !line.is_ascii() {
            return Err(PdbParseError::new(
                DiagnosticKind::Syntax,
                at_line,
                format!("{name} record contains non-ASCII text"),
            ));
        }
        let fields = Fields { line, number };
        match name {
            "ATOM" | "HETATM" => {
                if self.model.is_none() && !self.model_serials.is_empty() {
                    return Err(PdbParseError::new(
                        DiagnosticKind::Syntax,
                        at_line,
                        "coordinate record outside a MODEL block",
                    ));
                }
                if self.model.is_none() {
                    self.unmodelled_atom.get_or_insert(number);
                }
                let kind = if name == "ATOM" {
                    PdbAtomKind::Atom
                } else {
                    PdbAtomKind::Hetatm
                };
                let atom = fields.atom(kind, self.model)?;
                let document = self.document();
                if document.atoms.len() >= options.max_atom_records {
                    return Err(PdbParseError::new(
                        DiagnosticKind::Limit,
                        at_line,
                        "coordinate records exceed configured limit",
                    ));
                }
                document.atoms.push(atom);
            }
            "TER" => {
                let model = self.model;
                let document = self.document();
                let ter = fields.ter(document.atoms.len(), model)?;
                document.ter_records.push(ter);
            }
            "MODEL" => {
                if self.model.is_some() {
                    return Err(PdbParseError::new(
                        DiagnosticKind::Syntax,
                        at_line,
                        "MODEL record inside an open MODEL block",
                    ));
                }
                if let Some(line) = self.unmodelled_atom {
                    return Err(PdbParseError::new(
                        DiagnosticKind::Syntax,
                        SourceLocation::at_line(line),
                        "coordinate record outside a MODEL block",
                    ));
                }
                let serial = fields.model_serial()?;
                if !self.model_serials.insert(serial) {
                    return Err(PdbParseError::new(
                        DiagnosticKind::Inconsistent,
                        at_line,
                        format!("duplicate MODEL serial number {serial}"),
                    ));
                }
                self.model = Some(serial);
                self.document().models.push(PdbModel {
                    serial,
                    location: at_line,
                });
            }
            "ENDMDL" => {
                if self.model.take().is_none() {
                    return Err(PdbParseError::new(
                        DiagnosticKind::Syntax,
                        at_line,
                        "ENDMDL record without an open MODEL block",
                    ));
                }
            }
            "CONECT" => {
                let record = fields.conect()?;
                self.document().conect_records.push(record);
            }
            "LINK" => {
                let record = fields.link()?;
                self.document().link_records.push(record);
            }
            "SSBOND" => {
                let record = fields.ssbond()?;
                self.document().ssbond_records.push(record);
            }
            "END" => self.ended = true,
            _ => {
                let document = self.document();
                if name == "HEADER" && document.id_code.is_none() {
                    let code = fields.text(63, 66).trim();
                    if !code.is_empty() {
                        document.id_code = Some(code.to_owned());
                    }
                }
                document.other_records.push(PdbRecord {
                    name: name.to_owned(),
                    text: line.to_owned(),
                    location: at_line,
                });
            }
        }
        Ok(())
    }
}

/// The fixed columns of one ASCII record line, numbered from 1 as the format
/// specification numbers them.
struct Fields<'a> {
    line: &'a str,
    number: usize,
}

impl<'a> Fields<'a> {
    /// Columns `first..=last`, truncated where the line ends.
    fn text(&self, first: usize, last: usize) -> &'a str {
        let end = last.min(self.line.len());
        if first > end {
            return "";
        }
        self.line.get(first - 1..end).unwrap_or("")
    }

    fn character(&self, column: usize) -> Option<char> {
        self.text(column, column)
            .chars()
            .next()
            .filter(|character| *character != ' ')
    }

    fn error(&self, column: usize, message: impl Into<String>) -> PdbParseError {
        PdbParseError::new(
            DiagnosticKind::InvalidValue,
            SourceLocation::at_column(self.number, column),
            message,
        )
    }

    fn optional_text(&self, first: usize, last: usize) -> Option<String> {
        let text = self.text(first, last).trim();
        (!text.is_empty()).then(|| text.to_owned())
    }

    fn integer(&self, first: usize, last: usize, what: &str) -> Result<Option<i64>, PdbParseError> {
        hybrid36(self.text(first, last), last - first + 1)
            .map_err(|()| self.error(first, format!("invalid {what}")))
    }

    fn serial(&self, first: usize, last: usize, what: &str) -> Result<Option<u32>, PdbParseError> {
        // Some writers fill an overflowing serial number with asterisks.
        if self
            .text(first, last)
            .trim()
            .bytes()
            .all(|byte| byte == b'*')
        {
            return Ok(None);
        }
        self.integer(first, last, what)?
            .map(|value| {
                u32::try_from(value).map_err(|_| self.error(first, format!("invalid {what}")))
            })
            .transpose()
    }

    fn sequence(&self, first: usize, last: usize) -> Result<Option<i32>, PdbParseError> {
        self.integer(first, last, "residue sequence number")?
            .map(|value| {
                i32::try_from(value)
                    .map_err(|_| self.error(first, "residue sequence number out of range"))
            })
            .transpose()
    }

    fn required_sequence(&self, first: usize, last: usize) -> Result<i32, PdbParseError> {
        self.sequence(first, last)?
            .ok_or_else(|| self.error(first, "missing residue sequence number"))
    }

    fn real(&self, first: usize, last: usize, what: &str) -> Result<Option<f64>, PdbParseError> {
        let text = self.text(first, last).trim();
        if text.is_empty() {
            return Ok(None);
        }
        match text.parse::<f64>() {
            Ok(value) if value.is_finite() => Ok(Some(value)),
            _ => Err(self.error(first, format!("invalid {what}"))),
        }
    }

    fn coordinate(&self, first: usize, axis: &str) -> Result<f64, PdbParseError> {
        self.real(first, first + 7, &format!("{axis} coordinate"))?
            .ok_or_else(|| self.error(first, format!("missing {axis} coordinate")))
    }

    fn atom(&self, kind: PdbAtomKind, model: Option<u32>) -> Result<PdbAtom, PdbParseError> {
        let name_field = format!("{:<4}", self.text(13, 16));
        if name_field.trim().is_empty() {
            return Err(self.error(13, "missing atom name"));
        }
        Ok(PdbAtom {
            kind,
            serial: self.serial(7, 11, "atom serial number")?,
            name_field,
            alternate_location: self.character(17),
            residue_name: self.text(18, 21).trim().to_owned(),
            chain_id: self.character(22),
            residue_sequence: self.required_sequence(23, 26)?,
            insertion_code: self.character(27),
            position: Point3::new(
                self.coordinate(31, "x")?,
                self.coordinate(39, "y")?,
                self.coordinate(47, "z")?,
            ),
            occupancy: self.real(55, 60, "occupancy")?,
            temperature_factor: self.real(61, 66, "temperature factor")?,
            segment_id: self.optional_text(73, 76),
            element: self.optional_text(77, 78),
            formal_charge: self.charge()?,
            model,
            location: SourceLocation::at_line(self.number),
        })
    }

    /// A charge such as `2+` or `1-`; the sign may also lead.
    fn charge(&self) -> Result<Option<i8>, PdbParseError> {
        let text = self.text(79, 80).trim();
        if text.is_empty() {
            return Ok(None);
        }
        let (digits, sign) = match text.as_bytes() {
            [digit, sign @ (b'+' | b'-')] | [sign @ (b'+' | b'-'), digit]
                if digit.is_ascii_digit() =>
            {
                (digit - b'0', *sign)
            }
            _ => return Err(self.error(79, format!("invalid charge `{text}`"))),
        };
        let magnitude = i8::try_from(digits).expect("one digit");
        Ok(Some(if sign == b'-' { -magnitude } else { magnitude }))
    }

    fn ter(&self, preceding_atoms: usize, model: Option<u32>) -> Result<PdbTer, PdbParseError> {
        Ok(PdbTer {
            serial: self.serial(7, 11, "TER serial number")?,
            residue_name: self.optional_text(18, 21),
            chain_id: self.character(22),
            residue_sequence: self.sequence(23, 26)?,
            insertion_code: self.character(27),
            preceding_atoms,
            model,
            location: SourceLocation::at_line(self.number),
        })
    }

    fn model_serial(&self) -> Result<u32, PdbParseError> {
        // The serial belongs in columns 11-14; writers vary its alignment.
        let text = self.text(7, self.line.len()).trim();
        text.parse::<u32>()
            .map_err(|_| self.error(7, format!("invalid MODEL serial number `{text}`")))
    }

    fn conect(&self) -> Result<PdbConect, PdbParseError> {
        let serial = self
            .integer(7, 11, "CONECT serial number")?
            .ok_or_else(|| self.error(7, "missing CONECT serial number"))?;
        let serial =
            u32::try_from(serial).map_err(|_| self.error(7, "invalid CONECT serial number"))?;
        let mut bonded = Vec::new();
        // Columns 32-61 held hydrogen bonds and salt bridges in older versions
        // of the format; they are not bonds and are not read.
        for first in [12, 17, 22, 27] {
            if let Some(partner) = self.integer(first, first + 4, "CONECT serial number")? {
                bonded.push(
                    u32::try_from(partner)
                        .map_err(|_| self.error(first, "invalid CONECT serial number"))?,
                );
            }
        }
        Ok(PdbConect {
            serial,
            bonded,
            location: SourceLocation::at_line(self.number),
        })
    }

    fn partner(&self, columns: PartnerColumns) -> Result<PdbPartner, PdbParseError> {
        let residue_name = self
            .text(columns.residue_name, columns.residue_name + 3)
            .trim()
            .to_owned();
        if residue_name.is_empty() {
            return Err(self.error(columns.residue_name, "missing residue name"));
        }
        Ok(PdbPartner {
            atom_name: columns
                .atom_name
                .and_then(|first| self.optional_text(first, first + 3)),
            alternate_location: columns
                .atom_name
                .and_then(|first| self.character(first + 4)),
            residue_name,
            chain_id: self.character(columns.chain_id),
            residue_sequence: self
                .required_sequence(columns.residue_sequence, columns.residue_sequence + 3)?,
            insertion_code: self.character(columns.residue_sequence + 4),
            symmetry: self.optional_text(columns.symmetry, columns.symmetry + 5),
        })
    }

    fn link(&self) -> Result<PdbLink, PdbParseError> {
        let partners = [
            self.partner(PartnerColumns {
                atom_name: Some(13),
                residue_name: 18,
                chain_id: 22,
                residue_sequence: 23,
                symmetry: 60,
            })?,
            self.partner(PartnerColumns {
                atom_name: Some(43),
                residue_name: 48,
                chain_id: 52,
                residue_sequence: 53,
                symmetry: 67,
            })?,
        ];
        if partners.iter().any(|partner| partner.atom_name.is_none()) {
            return Err(self.error(13, "LINK record without an atom name"));
        }
        Ok(PdbLink {
            partners,
            distance: self.real(74, 78, "link distance")?,
            location: SourceLocation::at_line(self.number),
        })
    }

    fn ssbond(&self) -> Result<PdbSsbond, PdbParseError> {
        Ok(PdbSsbond {
            serial: self.serial(8, 10, "SSBOND serial number")?,
            partners: [
                self.partner(PartnerColumns {
                    atom_name: None,
                    residue_name: 12,
                    chain_id: 16,
                    residue_sequence: 18,
                    symmetry: 60,
                })?,
                self.partner(PartnerColumns {
                    atom_name: None,
                    residue_name: 26,
                    chain_id: 30,
                    residue_sequence: 32,
                    symmetry: 67,
                })?,
            ],
            distance: self.real(74, 78, "disulfide distance")?,
            location: SourceLocation::at_line(self.number),
        })
    }
}

/// First columns of one partner of a `LINK` or `SSBOND` record. The atom name
/// is followed by the alternate location, the sequence number by the
/// insertion code.
struct PartnerColumns {
    atom_name: Option<usize>,
    residue_name: usize,
    chain_id: usize,
    residue_sequence: usize,
    symmetry: usize,
}

/// Decodes a right-justified decimal or a hybrid-36 number of `width`
/// columns; blank fields are `None`.
///
/// Hybrid-36 continues the decimal range with base-36 numbers that begin
/// with an upper-case letter, then with ones that begin with a lower-case
/// letter, so a five-column field counts to 87,440,031.
fn hybrid36(field: &str, width: usize) -> Result<Option<i64>, ()> {
    let trimmed = field.trim();
    let Some(first) = trimmed.bytes().next() else {
        return Ok(None);
    };
    if first.is_ascii_digit() || first == b'-' || first == b'+' {
        return trimmed.parse::<i64>().map(Some).map_err(|_| ());
    }
    if field.len() != width || trimmed.len() != width {
        return Err(());
    }
    let upper = first.is_ascii_uppercase();
    if !upper && !first.is_ascii_lowercase() {
        return Err(());
    }
    let mut value = 0i64;
    for byte in field.bytes() {
        let digit = match byte {
            b'0'..=b'9' => byte - b'0',
            b'A'..=b'Z' if upper => byte - b'A' + 10,
            b'a'..=b'z' if !upper => byte - b'a' + 10,
            _ => return Err(()),
        };
        value = value * 36 + i64::from(digit);
    }
    let exponent = u32::try_from(width).map_err(|_| ())?;
    let span = 36i64.pow(exponent - 1);
    let decimal = 10i64.pow(exponent);
    Ok(Some(if upper {
        value - 10 * span + decimal
    } else {
        value + 16 * span + decimal
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hybrid36_continues_the_decimal_range() {
        for (field, width, value) in [
            ("99999", 5, 99_999),
            ("A0000", 5, 100_000),
            ("A0001", 5, 100_001),
            ("ZZZZZ", 5, 100_000 + 26 * 36i64.pow(4) - 1),
            ("a0000", 5, 100_000 + 26 * 36i64.pow(4)),
            ("zzzzz", 5, 87_440_031),
            ("9999", 4, 9_999),
            ("A000", 4, 10_000),
            ("  -1", 4, -1),
        ] {
            assert_eq!(hybrid36(field, width), Ok(Some(value)), "{field}");
        }
        assert_eq!(hybrid36("    ", 4), Ok(None));
        for invalid in ["A0a0", " A00", "1A00", "!000"] {
            assert_eq!(hybrid36(invalid, 4), Err(()), "{invalid}");
        }
    }
}
