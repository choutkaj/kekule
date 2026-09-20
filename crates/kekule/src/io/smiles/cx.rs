//! Chemical CX fields, resolved against record-global SMILES atom indices.
use std::collections::{BTreeMap, BTreeSet};

use crate::core::{
    AtomId, AtomRadical, MoleculeEditor, StereoElementKind, StereoGroup, StereoGroupKind,
};

use super::interpret::SmilesInterpretError;
use super::parse::SmilesDocument;

#[derive(Default)]
pub(super) struct CxFeatures {
    pub(super) radicals: BTreeMap<usize, AtomRadical>,
    pub(super) groups: Vec<CxStereoGroup>,
}

pub(super) struct CxStereoGroup {
    pub(super) kind: StereoGroupKind,
    pub(super) atoms: Vec<usize>,
    pub(super) offset: usize,
}

impl CxFeatures {
    pub(super) fn interpret(document: &SmilesDocument) -> Result<Self, SmilesInterpretError> {
        let Some(span) = document.cx_extension_span() else {
            return Ok(Self::default());
        };
        let mut reader = CxReader {
            text: &document.source()[span.start + 1..span.end - 1],
            offset: span.start + 1,
            cursor: 0,
        };
        let mut features = Self::default();
        let mut group_indices = BTreeMap::new();
        let mut grouped_atoms = BTreeSet::new();
        let mut relative = None;
        reader.skip_space();
        while let Some(tag) = reader.peek() {
            let offset = reader.position();
            reader.cursor += 1;
            match tag {
                b'^' => {
                    let code = reader.number()?;
                    let (electrons, spin) = match code {
                        1 => (1, None),
                        2 => (2, None),
                        3 => (2, Some(1)),
                        4 => (2, Some(3)),
                        5 => (3, None),
                        6 => (3, Some(2)),
                        7 => (3, Some(4)),
                        _ => return Err(error(offset, "invalid CXSMILES radical code")),
                    };
                    reader.expect(b':')?;
                    for atom in reader.atoms(document)? {
                        if features
                            .radicals
                            .insert(
                                atom,
                                AtomRadical::new(electrons, spin)
                                    .expect("CX radical codes specify valid electron/spin states"),
                            )
                            .is_some()
                        {
                            return Err(error(offset, "duplicate CXSMILES radical atom"));
                        }
                    }
                }
                b'a' | b'o' | b'&' => {
                    let number = if tag == b'a' { 0 } else { reader.number()? };
                    reader.expect(b':')?;
                    let atoms = reader.atoms(document)?;
                    for atom in &atoms {
                        if !grouped_atoms.insert(*atom) {
                            return Err(error(offset, "duplicate CXSMILES stereo group atom"));
                        }
                    }
                    let kind = match tag {
                        b'a' => StereoGroupKind::Absolute,
                        b'o' => StereoGroupKind::Or,
                        _ => StereoGroupKind::And,
                    };
                    let index = *group_indices.entry((tag, number)).or_insert_with(|| {
                        features.groups.push(CxStereoGroup {
                            kind,
                            atoms: Vec::new(),
                            offset,
                        });
                        features.groups.len() - 1
                    });
                    features.groups[index].atoms.extend(atoms);
                }
                b'r' if matches!(reader.peek(), None | Some(b',' | b' ' | b'\t')) => {
                    if relative.replace(offset).is_some() {
                        return Err(error(
                            offset,
                            "duplicate CXSMILES relative configuration flag",
                        ));
                    }
                }
                _ => return Err(error(offset, "unsupported CXSMILES chemical field")),
            }
            reader.skip_space();
            if reader.peek().is_some() {
                reader.expect(b',')?;
                reader.skip_space();
                if reader.peek().is_none() {
                    return Err(error(
                        reader.position(),
                        "trailing CXSMILES field separator",
                    ));
                }
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
                return Err(error(group.offset,
                    "CXSMILES stereo relationships spanning disconnected molecules cannot be represented"));
            }
        }
        Ok(features)
    }
}

fn error(offset: usize, message: &str) -> SmilesInterpretError {
    SmilesInterpretError {
        offset,
        message: message.to_owned(),
    }
}

pub(super) fn install_stereo_groups(
    editor: &mut MoleculeEditor,
    source_to_atom: &BTreeMap<usize, AtomId>,
    groups: &[CxStereoGroup],
) -> Result<(), SmilesInterpretError> {
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
                        error(group.offset,
                "CXSMILES stereo group atom must identify a represented tetrahedral center")
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
                    group.offset,
                    &format!("invalid CXSMILES stereo group: {cause}"),
                )
            })?;
    }
    Ok(())
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

    fn expect(&mut self, expected: u8) -> Result<(), SmilesInterpretError> {
        self.skip_space();
        if self.peek() != Some(expected) {
            return Err(error(self.position(), "invalid CXSMILES field syntax"));
        }
        self.cursor += 1;
        Ok(())
    }

    fn number(&mut self) -> Result<usize, SmilesInterpretError> {
        self.skip_space();
        let start = self.position();
        if !self.peek().is_some_and(|ch| ch.is_ascii_digit()) {
            return Err(error(start, "expected unsigned CXSMILES index"));
        }
        let mut value = 0usize;
        while let Some(ch @ b'0'..=b'9') = self.peek() {
            value = value
                .checked_mul(10)
                .and_then(|value| value.checked_add(usize::from(ch - b'0')))
                .ok_or_else(|| error(start, "CXSMILES index exceeds supported range"))?;
            self.cursor += 1;
        }
        Ok(value)
    }

    fn atoms(&mut self, document: &SmilesDocument) -> Result<Vec<usize>, SmilesInterpretError> {
        let mut atoms = Vec::new();
        loop {
            self.skip_space();
            let offset = self.position();
            let atom = self.number()?;
            if atom >= document.program.atoms.len() {
                return Err(error(
                    offset,
                    "CXSMILES atom index is outside the base SMILES",
                ));
            }
            atoms.push(atom);
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
