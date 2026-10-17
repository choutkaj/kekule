//! CXSMILES extension fields: typed while parsing, then interpreted against
//! record-global SMILES atom indices.
use std::collections::{BTreeMap, BTreeSet};
use std::ops::Range;

use crate::diagnostics::{DiagnosticKind, UnsupportedContent};
use crate::molecule::{AtomId, AtomRadical, MoleculeEditor};
use crate::stereo::{StereoElementKind, StereoGroup, StereoGroupKind};

use super::interpret::interpret_error;
use super::parse::{parse_error, SmilesDocument};
use crate::diagnostics::FormatError;

/// Kind and number of a CXSMILES enhanced stereo group.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SmilesCxStereoGroupKind {
    /// `a:`.
    Absolute,
    /// `on:`.
    Or(usize),
    /// `&n:`.
    And(usize),
}

/// What one CXSMILES field states. Atom indices are zero-based positions in
/// the base SMILES, as written.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SmilesCxFieldKind {
    /// `^n:` radical declaration, decoded from its code into an unpaired
    /// electron count and, when the code states one, a spin multiplicity.
    Radical {
        electrons: u8,
        multiplicity: Option<u8>,
        atoms: Vec<usize>,
    },
    /// Enhanced stereo group membership.
    StereoGroup {
        kind: SmilesCxStereoGroupKind,
        atoms: Vec<usize>,
    },
    /// `r`: tetrahedral centers outside stereo groups share one relative
    /// configuration.
    RelativeConfiguration,
    /// A field Kekule does not interpret, as written.
    Unsupported(String),
}

/// One field of a CXSMILES extension and its byte span in the record.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SmilesCxField {
    kind: SmilesCxFieldKind,
    span: Range<usize>,
}

impl SmilesCxField {
    pub fn kind(&self) -> &SmilesCxFieldKind {
        &self.kind
    }

    pub fn span(&self) -> Range<usize> {
        self.span.clone()
    }
}

/// Parses the fields of the extension at `span`, including its vertical bars.
pub(super) fn parse_cx_fields(
    source: &str,
    span: Range<usize>,
) -> Result<Vec<SmilesCxField>, FormatError> {
    let mut reader = CxReader {
        text: &source[span.start + 1..span.end - 1],
        offset: span.start + 1,
        cursor: 0,
    };
    let mut fields = Vec::new();
    reader.skip_space();
    while reader.peek().is_some() {
        let start = reader.position();
        let kind = reader.field()?;
        fields.push(SmilesCxField {
            kind,
            span: start..reader.position(),
        });
        reader.skip_space();
        if reader.peek().is_some() {
            reader.expect(b',')?;
            reader.skip_space();
            if reader.peek().is_none() {
                return Err(syntax(
                    reader.position(),
                    "trailing CXSMILES field separator",
                ));
            }
        }
    }
    Ok(fields)
}

fn syntax(offset: usize, message: &str) -> FormatError {
    parse_error(DiagnosticKind::Syntax, offset, message)
}

struct CxReader<'a> {
    text: &'a str,
    offset: usize,
    cursor: usize,
}

impl CxReader<'_> {
    fn position(&self) -> usize {
        self.offset + self.cursor
    }

    fn peek(&self) -> Option<u8> {
        self.text.as_bytes().get(self.cursor).copied()
    }

    fn skip_space(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t')) {
            self.cursor += 1;
        }
    }

    fn expect(&mut self, expected: u8) -> Result<(), FormatError> {
        self.skip_space();
        if self.peek() != Some(expected) {
            return Err(syntax(self.position(), "invalid CXSMILES field syntax"));
        }
        self.cursor += 1;
        Ok(())
    }

    /// Reads one field. Radicals, stereo groups, and the relative flag are
    /// typed; any other field is kept as written.
    fn field(&mut self) -> Result<SmilesCxFieldKind, FormatError> {
        let start = self.cursor;
        let rest = &self.text[start..];
        let tag_len = rest.find([':', ',', ' ', '\t']).unwrap_or(rest.len());
        let tag = &rest[..tag_len];
        let has_atoms = rest[tag_len..].starts_with(':');
        if tag.is_empty() {
            return Err(syntax(self.position(), "empty CXSMILES field"));
        }
        if tag == "r" && !has_atoms {
            self.cursor += 1;
            return Ok(SmilesCxFieldKind::RelativeConfiguration);
        }
        let kind = match (tag.as_bytes()[0], has_atoms) {
            (b'a', true) if tag == "a" => None,
            (b'o' | b'&' | b'^', true) => Some(self.tag_number(&tag[1..])?),
            _ => {
                self.skip_unsupported();
                let text = self.text[start..self.cursor].trim_end();
                self.cursor = start + text.len();
                return Ok(SmilesCxFieldKind::Unsupported(text.to_owned()));
            }
        };
        let code_offset = self.position() + 1;
        self.cursor += tag_len + 1;
        let atoms = self.atoms()?;
        Ok(match (tag.as_bytes()[0], kind) {
            (b'^', Some(code)) => {
                let (electrons, multiplicity) = match code {
                    1 => (1, None),
                    2 => (2, None),
                    3 => (2, Some(1)),
                    4 => (2, Some(3)),
                    5 => (3, None),
                    6 => (3, Some(2)),
                    7 => (3, Some(4)),
                    _ => {
                        return Err(parse_error(
                            DiagnosticKind::InvalidValue,
                            code_offset,
                            "invalid CXSMILES radical code",
                        ))
                    }
                };
                SmilesCxFieldKind::Radical {
                    electrons,
                    multiplicity,
                    atoms,
                }
            }
            (b'o', Some(number)) => SmilesCxFieldKind::StereoGroup {
                kind: SmilesCxStereoGroupKind::Or(number),
                atoms,
            },
            (b'&', Some(number)) => SmilesCxFieldKind::StereoGroup {
                kind: SmilesCxStereoGroupKind::And(number),
                atoms,
            },
            _ => SmilesCxFieldKind::StereoGroup {
                kind: SmilesCxStereoGroupKind::Absolute,
                atoms,
            },
        })
    }

    /// The number after a group or radical tag character.
    fn tag_number(&self, digits: &str) -> Result<usize, FormatError> {
        let offset = self.position() + 1;
        if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
            return Err(syntax(offset, "expected unsigned CXSMILES index"));
        }
        digits.parse().map_err(|_| {
            parse_error(
                DiagnosticKind::InvalidValue,
                offset,
                "CXSMILES index exceeds supported range",
            )
        })
    }

    /// Skips an unsupported field. Commas inside parentheses or `$` labels
    /// do not end it, nor does a comma that continues an index list.
    fn skip_unsupported(&mut self) {
        let mut depth = 0usize;
        let mut label = false;
        while let Some(byte) = self.peek() {
            match byte {
                b'$' => label = !label,
                b'(' if !label => depth += 1,
                b')' if !label => depth = depth.saturating_sub(1),
                b',' if !label && depth == 0 => {
                    let next = self.text.as_bytes()[self.cursor + 1..]
                        .iter()
                        .find(|byte| !matches!(byte, b' ' | b'\t'));
                    if !next.is_some_and(u8::is_ascii_digit) {
                        break;
                    }
                }
                _ => {}
            }
            self.cursor += 1;
        }
    }

    fn number(&mut self) -> Result<usize, FormatError> {
        self.skip_space();
        let start = self.position();
        if !self.peek().is_some_and(|ch| ch.is_ascii_digit()) {
            return Err(syntax(start, "expected unsigned CXSMILES index"));
        }
        let mut value = 0usize;
        while let Some(ch @ b'0'..=b'9') = self.peek() {
            value = value
                .checked_mul(10)
                .and_then(|value| value.checked_add(usize::from(ch - b'0')))
                .ok_or_else(|| {
                    parse_error(
                        DiagnosticKind::InvalidValue,
                        start,
                        "CXSMILES index exceeds supported range",
                    )
                })?;
            self.cursor += 1;
        }
        Ok(value)
    }

    fn atoms(&mut self) -> Result<Vec<usize>, FormatError> {
        let mut atoms = Vec::new();
        loop {
            atoms.push(self.number()?);
            self.skip_space();
            if self.peek() != Some(b',') {
                break;
            }
            let separator = self.cursor;
            self.cursor += 1;
            self.skip_space();
            if !self.peek().is_some_and(|ch| ch.is_ascii_digit()) {
                self.cursor = separator;
                break;
            }
        }
        Ok(atoms)
    }
}

#[derive(Default)]
pub(super) struct CxFeatures {
    pub(super) radicals: BTreeMap<usize, AtomRadical>,
    pub(super) groups: Vec<CxStereoGroup>,
    pub(super) omitted: Vec<SmilesCxField>,
}

pub(super) struct CxStereoGroup {
    pub(super) kind: StereoGroupKind,
    pub(super) atoms: Vec<usize>,
    pub(super) offset: usize,
}

impl CxFeatures {
    pub(super) fn interpret(
        document: &SmilesDocument,
        unsupported: UnsupportedContent,
    ) -> Result<Self, FormatError> {
        let atom_count = document.program.atoms.len();
        let mut features = Self::default();
        let mut group_indices = BTreeMap::new();
        let mut grouped_atoms = BTreeSet::new();
        let mut relative = None;
        for field in document.cx_fields() {
            let offset = field.span.start;
            let inconsistent =
                |message: &str| error(DiagnosticKind::Inconsistent, offset, message.to_owned());
            let atoms = match &field.kind {
                SmilesCxFieldKind::Radical { atoms, .. }
                | SmilesCxFieldKind::StereoGroup { atoms, .. } => atoms.as_slice(),
                _ => &[],
            };
            if atoms.iter().any(|&atom| atom >= atom_count) {
                return Err(error(
                    DiagnosticKind::InvalidReference,
                    offset,
                    "CXSMILES atom index is outside the base SMILES".to_owned(),
                ));
            }
            match &field.kind {
                SmilesCxFieldKind::Radical {
                    electrons,
                    multiplicity,
                    atoms,
                } => {
                    let radical = AtomRadical::new(*electrons, *multiplicity)
                        .expect("CX radical codes specify valid electron/spin states");
                    for &atom in atoms {
                        if features.radicals.insert(atom, radical).is_some() {
                            return Err(inconsistent("duplicate CXSMILES radical atom"));
                        }
                    }
                }
                SmilesCxFieldKind::StereoGroup { kind, atoms } => {
                    for &atom in atoms {
                        if !grouped_atoms.insert(atom) {
                            return Err(inconsistent("duplicate CXSMILES stereo group atom"));
                        }
                    }
                    let index = *group_indices.entry(*kind).or_insert_with(|| {
                        features.groups.push(CxStereoGroup {
                            kind: match kind {
                                SmilesCxStereoGroupKind::Absolute => StereoGroupKind::Absolute,
                                SmilesCxStereoGroupKind::Or(_) => StereoGroupKind::Or,
                                SmilesCxStereoGroupKind::And(_) => StereoGroupKind::And,
                            },
                            atoms: Vec::new(),
                            offset,
                        });
                        features.groups.len() - 1
                    });
                    features.groups[index].atoms.extend(atoms);
                }
                SmilesCxFieldKind::RelativeConfiguration => {
                    if relative.replace(offset).is_some() {
                        return Err(inconsistent(
                            "duplicate CXSMILES relative configuration flag",
                        ));
                    }
                }
                SmilesCxFieldKind::Unsupported(text) => match unsupported {
                    UnsupportedContent::Report => features.omitted.push(field.clone()),
                    UnsupportedContent::Reject => {
                        return Err(error(
                            DiagnosticKind::Unsupported,
                            offset,
                            format!("unsupported CXSMILES field `{text}`"),
                        ))
                    }
                },
            }
        }
        if let Some(offset) = relative {
            // Enhanced groups take precedence over the legacy chiral flag.
            // Preserve the common relative configuration of the remaining
            // tetrahedral assertions, without asserting a pure or racemic sample.
            let atoms = document
                .program
                .tetrahedral
                .iter()
                .map(|stereo| stereo.center)
                .filter(|atom| !grouped_atoms.contains(atom))
                .collect::<Vec<_>>();
            if !atoms.is_empty() {
                features.groups.push(CxStereoGroup {
                    kind: StereoGroupKind::Relative,
                    atoms,
                    offset,
                });
            }
        }
        for group in &features.groups {
            let component = document.program.atoms[group.atoms[0]].component;
            if group.kind != StereoGroupKind::Absolute
                && group
                    .atoms
                    .iter()
                    .any(|atom| document.program.atoms[*atom].component != component)
            {
                return Err(error(
                    DiagnosticKind::Unsupported,
                    group.offset,
                    "CXSMILES stereo relationships spanning disconnected molecules cannot be represented"
                        .to_owned(),
                ));
            }
        }
        Ok(features)
    }
}

fn error(kind: DiagnosticKind, offset: usize, message: String) -> FormatError {
    interpret_error(kind, offset, message)
}

pub(super) fn install_stereo_groups(
    editor: &mut MoleculeEditor,
    source_to_atom: &BTreeMap<usize, AtomId>,
    groups: &[CxStereoGroup],
) -> Result<(), FormatError> {
    if groups.is_empty() {
        return Ok(());
    }
    let centers = editor
        .stereo_elements()
        .filter_map(|(id, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => Some((stereo.center, id)),
            _ => None,
        })
        .collect::<BTreeMap<_, _>>();
    for group in groups {
        let members = group
            .atoms
            .iter()
            .map(|source| {
                centers
                    .get(&source_to_atom[source])
                    .copied()
                    .ok_or_else(|| {
                        error(
                            DiagnosticKind::Chemistry,
                            group.offset,
                            "CXSMILES stereo group atom must identify a represented tetrahedral center"
                                .to_owned(),
                        )
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        editor
            .add_stereo_group(StereoGroup {
                kind: group.kind,
                members,
            })
            .map_err(|cause| {
                error(
                    DiagnosticKind::Chemistry,
                    group.offset,
                    format!("invalid CXSMILES stereo group: {cause}"),
                )
            })?;
    }
    Ok(())
}
