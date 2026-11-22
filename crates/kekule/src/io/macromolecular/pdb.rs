//! Reads a PDB entry into the source-neutral interpretation input.
//!
//! PDB records carry author identities only. A residue is a run of
//! consecutive coordinate records with one chain, sequence number, and
//! insertion code, which a `TER` record also ends; a number that recurs after
//! other residues or a `TER` record names a new residue, as when numbering
//! wraps around or a segment restarts it. Without entity records, a chain's
//! polymer is read from its `TER` records, and its residues' order gives
//! polymer adjacency. A residue of `HETATM` records in a run a `TER` record
//! ends is polymer when its template is a polymer one or it has a linking
//! atom (`N`, `C`, `P`, `O3'`, or `O3*`). Without `TER`, it is polymer when
//! its template is a polymer one or one of its linking atoms lies within
//! bonding distance of one of a neighbour's, as for a capping group or a
//! 5'-terminal nucleotide triphosphate, so a ligand listed between two
//! residues does not separate them. Such residues bonded to the first or last
//! `ATOM` residue extend the polymer. The polymer residues of a chain whose
//! `SEQRES` sequence places them take their places in it as label sequence
//! numbers. Only such a sequence says which residues a polymer lacks.
use std::collections::{BTreeMap, BTreeSet, HashMap};

use super::types::{poll, poll_every, structure_error};
use crate::diagnostics::{Cancellation, DiagnosticKind, FormatError};
use crate::molecule::{BondOrder, Element};

use super::super::element_symbol::parse_structure_element_symbol;
use super::super::{PdbAtom, PdbAtomKind, PdbConect, PdbDocument, PdbPartner};
use super::altloc::inventory;
use super::assemblies::read_pdb_assemblies;
use super::atom_site::AtomRow;
use super::connections::{
    AuthorConnectionSelector, ConnectionPartnerSelector, ConnectionRecord, DeclaredBond,
    DeclaredOrder, PartnerSpec,
};
use super::connectivity::{within_bonding_distance, ConnectivityCatalog};
use super::sequences::sequence_entities;
use super::structure_set::{interpret_conformations, interpret_models};
use super::types::{
    ConnectionResolutionReason, ConnectivityOptions, EntityKind, NumberField,
    StructureInterpretIssue, StructureInterpretOptions, StructureInterpretation,
    StructureSetInterpretOptions, StructureSetInterpretation,
};
use super::{interpret_source, stated_cell, AltLocResidue, StructureSource};
use crate::components::{ComponentLibrary, ComponentLinking};
use crate::math::Point3;
use crate::topology::{is_ion, is_water};

pub(crate) fn interpret_pdb(
    document: &PdbDocument,
    options: StructureInterpretOptions,
) -> Result<StructureInterpretation, FormatError> {
    let source = pdb_source(
        document,
        options.requires_entity_metadata(),
        options.connectivity(),
        options.cancellation(),
    )?;
    interpret_source(source, &options)
}

/// Interprets explicitly selected or all coordinate models as one
/// shared-topology unordered structure set.
pub(crate) fn interpret_pdb_structure_set(
    document: &PdbDocument,
    options: StructureSetInterpretOptions,
) -> Result<StructureSetInterpretation, FormatError> {
    let mut seen = BTreeSet::new();
    let available = document
        .atoms()
        .iter()
        .map(PdbAtom::model)
        .filter(|model| seen.insert(*model))
        .map(model_id)
        .collect();
    interpret_models(available, options, |options| {
        pdb_source(
            document,
            options.requires_entity_metadata(),
            options.connectivity(),
            options.cancellation(),
        )
    })
}

/// Explicit caller-selected alternate conformations, never automatic altloc
/// enumeration.
pub(crate) fn interpret_pdb_conformations(
    document: &PdbDocument,
    selections: &[StructureInterpretOptions],
) -> Result<StructureSetInterpretation, FormatError> {
    interpret_conformations(selections, |options| {
        interpret_pdb(document, options.clone())
    })
}

pub(crate) fn pdb_alternate_locations(
    document: &PdbDocument,
) -> Result<Vec<AltLocResidue>, FormatError> {
    inventory(&read_rows(document, None, None)?.0)
}

fn pdb_source(
    document: &PdbDocument,
    strict_entity_metadata: bool,
    connectivity: &ConnectivityOptions,
    cancellation: Option<&Cancellation>,
) -> Result<StructureSource, FormatError> {
    if strict_entity_metadata {
        return Err(structure_error(
            DiagnosticKind::Inconsistent,
            None,
            "A PDB entry declares no entity types. Turn off strict entity metadata.",
        ));
    }
    let (mut rows, mut issues) = read_rows(
        document,
        connectivity.components().map(AsRef::as_ref),
        cancellation,
    )?;
    let sequences = align_sequences(document, &mut rows);
    let chain_lengths = sequences
        .iter()
        .filter_map(|(chain, monomers)| Some((chain.clone(), i32::try_from(monomers.len()).ok()?)))
        .collect();
    let polymer_entities =
        sequence_entities(&sequences, connectivity.components().map(AsRef::as_ref));
    poll(cancellation)?;
    let mut polymer_order = BTreeMap::new();
    for row in rows.iter().filter(|row| row.polymer_position.is_some()) {
        let next = polymer_order.len();
        polymer_order.entry(row.asym_id.clone()).or_insert(next);
    }
    let cell = match document.cryst1() {
        None => None,
        Some(record) => {
            let [a, b, c] = record.lengths();
            let [alpha, beta, gamma] = record.angles();
            let space_groups = record.space_group().into_iter().collect::<Vec<_>>();
            match stated_cell([a, b, c, alpha, beta, gamma], &space_groups) {
                Ok(Ok(cell)) => Some(cell),
                Ok(Err(reason)) => {
                    issues.push(StructureInterpretIssue::CellOmitted {
                        reason,
                        location: record.location(),
                    });
                    None
                }
                Err(error) => {
                    return Err(structure_error(
                        DiagnosticKind::InvalidValue,
                        Some(record.location().line()),
                        "The CRYST1 record is not a valid periodic cell.",
                    )
                    .with_source(error))
                }
            }
        }
    };
    poll(cancellation)?;
    let connections = connection_records(document)?;
    let assemblies = read_pdb_assemblies(document.biomolecules(), &mut issues);
    poll(cancellation)?;
    Ok(StructureSource {
        name: document.id_code().unwrap_or_default().to_owned(),
        entity_definitions: 0,
        rows,
        conformations: BTreeMap::new(),
        polymer_order,
        connections,
        catalog: ConnectivityCatalog::new(connectivity).with_chain_lengths(chain_lengths),
        cell,
        issues,
        polymer_entities,
        assemblies,
    })
}

/// Numbers the polymer residues of each chain whose `SEQRES` sequence
/// places them by their places in it, which become their label sequence
/// numbers. Returns, by chain, each such sequence's residue names. Other
/// chains keep their record-order positions and have no label sequence
/// numbers: without a sequence, nothing says which residues their polymer
/// lacks.
fn align_sequences(document: &PdbDocument, rows: &mut [AtomRow]) -> BTreeMap<String, Vec<String>> {
    let (sequences, places) = sequence_places(document, rows);
    for row in rows.iter_mut() {
        if let (Some(position), Some(places)) = (row.polymer_position, places.get(&row.asym_id)) {
            row.polymer_position = places.get(&position).copied().or(Some(position));
            row.label_seq_id = row.polymer_position;
        }
    }
    sequences
}

/// The sequences of [`align_sequences`] and, by chain, each
/// record-order position's place in the sequence.
type SequencePlaces = (
    BTreeMap<String, Vec<String>>,
    BTreeMap<String, BTreeMap<i32, i32>>,
);

/// The modelled residues of each chain form segments, split by `TER`
/// records and by jumps in residue numbers, and each segment matches a run
/// of the sequence's residue names (see [`segment_starts`]). A chain whose
/// segments match no such runs keeps its record-order positions.
fn sequence_places(document: &PdbDocument, rows: &[AtomRow]) -> SequencePlaces {
    let mut sequences = BTreeMap::<String, Vec<&str>>::new();
    for record in document.other_records() {
        if record.name() != "SEQRES" {
            continue;
        }
        let text = record.text();
        let Some(chain) = text.get(11..12) else {
            continue;
        };
        sequences
            .entry(chain.to_owned())
            .or_default()
            .extend(text.get(19..).unwrap_or("").split_whitespace());
    }
    if sequences.is_empty() {
        return Default::default();
    }
    let first_model = rows.first().map(|row| row.model_id.as_str());
    // Each chain's modelled polymer residues in order: name, position, and
    // number.
    let mut chains = BTreeMap::<&str, Vec<(&str, i32, Option<i32>)>>::new();
    let mut last = None;
    for row in rows {
        let Some(position) = row.polymer_position else {
            continue;
        };
        if Some(row.model_id.as_str()) != first_model {
            continue;
        }
        let residue = (row.asym_id.as_str(), position);
        if last == Some(residue) {
            continue;
        }
        last = Some(residue);
        let number = row
            .auth_seq_id
            .as_deref()
            .and_then(|text| text.parse().ok());
        chains
            .entry(&row.asym_id)
            .or_default()
            .push((&row.comp_id, position, number));
    }
    let mut aligned = BTreeMap::new();
    let mut places = BTreeMap::new();
    for (chain, residues) in chains {
        let Some(sequence) = sequences.get(chain) else {
            continue;
        };
        let mut segments = Vec::<Segment<'_>>::new();
        let mut positions = Vec::<Vec<i32>>::new();
        for (index, &(name, position, number)) in residues.iter().enumerate() {
            let continues = index > 0 && {
                let (_, previous_position, previous_number) = residues[index - 1];
                position == previous_position + 1
                    && match (previous_number, number) {
                        (Some(previous), Some(number)) => matches!(number - previous, 0 | 1),
                        _ => true,
                    }
            };
            match (segments.last_mut(), positions.last_mut()) {
                (Some(segment), Some(held)) if continues => {
                    segment.names.push(name);
                    held.push(position);
                }
                _ => {
                    segments.push(Segment {
                        names: vec![name],
                        number,
                    });
                    positions.push(vec![position]);
                }
            }
        }
        let Some(starts) = segment_starts(sequence, &segments) else {
            continue;
        };
        aligned.insert(
            chain.to_owned(),
            sequence.iter().map(|&name| name.to_owned()).collect(),
        );
        let mut placed = BTreeMap::new();
        for (held, &start) in positions.iter().zip(&starts) {
            for (offset, &position) in held.iter().enumerate() {
                let Ok(place) = i32::try_from(start + offset + 1) else {
                    continue;
                };
                placed.insert(position, place);
            }
        }
        places.insert(chain.to_owned(), placed);
    }
    (aligned, places)
}

/// A run of modelled residues that follow each other in a chain: their
/// names, and the residue number of the first.
pub(crate) struct Segment<'a> {
    pub(crate) names: Vec<&'a str>,
    pub(crate) number: Option<i32>,
}

/// Where each segment starts in `sequence`, by zero-based index, or `None`
/// when they do not all fit in order. Segments follow the residue numbers
/// when one offset between number and place fits every segment, the
/// smallest such offset. Otherwise the first segment takes its earliest
/// place, a lone residue its latest before the second segment, later ones
/// the earliest after the one before, and the last its latest.
pub(crate) fn segment_starts(sequence: &[&str], segments: &[Segment<'_>]) -> Option<Vec<usize>> {
    let fits = |segment: &Segment<'_>, start: usize| {
        sequence.get(start..start + segment.names.len()) == Some(&segment.names[..])
    };
    let first = segments.first()?;
    if let Some(numbers) = segments
        .iter()
        .map(|segment| segment.number)
        .collect::<Option<Vec<_>>>()
    {
        for start in (0..sequence.len()).filter(|&start| fits(first, start)) {
            let offset = i64::try_from(start).ok()? - i64::from(numbers[0]);
            let mut end = 0;
            let starts = segments
                .iter()
                .zip(&numbers)
                .map(|(segment, &number)| {
                    let start = usize::try_from(i64::from(number) + offset).ok()?;
                    (start >= end && fits(segment, start)).then(|| {
                        end = start + segment.names.len();
                        start
                    })
                })
                .collect::<Option<Vec<_>>>();
            if starts.is_some() {
                return starts;
            }
        }
    }
    let mut starts = Vec::with_capacity(segments.len());
    let mut from = 0;
    for segment in segments {
        let start = (from..=sequence.len()).find(|&start| fits(segment, start))?;
        starts.push(start);
        from = start + segment.names.len();
    }
    let latest = |index: usize, limit: usize| {
        (starts[index]..=limit.saturating_sub(segments[index].names.len()))
            .rev()
            .find(|&start| fits(&segments[index], start))
    };
    if segments.len() > 1 {
        let last = segments.len() - 1;
        let lone = (first.names.len() == 1)
            .then(|| latest(0, starts[1]))
            .flatten();
        let end = latest(last, sequence.len());
        if let Some(start) = lone {
            starts[0] = start;
        }
        if let Some(start) = end {
            starts[last] = start;
        }
    }
    Some(starts)
}

fn model_id(model: Option<u32>) -> String {
    model.map_or_else(|| "1".to_owned(), |serial| serial.to_string())
}

/// The chain identifier, or the segment identifier when the chain is blank,
/// as in files that name chains by segment.
fn asym_id(atom: &PdbAtom) -> String {
    match (atom.chain_id(), atom.segment_id()) {
        (Some(chain), _) => chain.to_string(),
        (None, Some(segment)) => segment.to_owned(),
        (None, None) => String::new(),
    }
}

/// One run of coordinate records forming a residue.
struct Residue {
    asym_id: String,
    /// Earlier residues of the model with the same chain, number, and
    /// insertion code.
    repeat: usize,
    atom_record: bool,
    water: bool,
    /// A water or an ion, which is never polymer.
    solvent: bool,
    /// Whether its template is a polymer one.
    polymer_template: bool,
    /// Its atoms that a polymer link joins to the residue before it (`N`,
    /// `P`) and to the one after it (`C`, `O3'`).
    heads: Vec<(Element, Point3)>,
    tails: Vec<(Element, Point3)>,
    polymer_position: Option<i32>,
}

impl Residue {
    /// Whether it has an atom a polymer link joins.
    fn has_linking_atom(&self) -> bool {
        !self.heads.is_empty() || !self.tails.is_empty()
    }

    /// Whether a linking atom of this residue lies within bonding distance
    /// of one of `next`'s.
    fn links_to(&self, next: &Self) -> bool {
        self.tails.iter().any(|&tail| {
            next.heads
                .iter()
                .any(|&head| within_bonding_distance(tail, head))
        })
    }
}

/// The residues of each chain since its last `TER` record, within one model.
#[derive(Default)]
struct ChainRuns {
    open: BTreeMap<String, Vec<usize>>,
    next_position: BTreeMap<String, i32>,
}

impl ChainRuns {
    /// Ends a chain's run. A run a `TER` record ends is polymer; otherwise its
    /// polymer spans its `ATOM` records and the residues bonded to their
    /// ends. Water and ions are never polymer, nor a single residue, nor a
    /// residue of `HETATM` records without a polymer template and, in a run
    /// a `TER` record ends, without linking atoms, or otherwise without a
    /// link to a neighbour.
    fn close(&mut self, asym_id: &str, terminated: bool, residues: &mut [Residue]) {
        let Some(run) = self.open.remove(asym_id) else {
            return;
        };
        let span = if terminated {
            &run[..]
        } else {
            let first = run
                .iter()
                .position(|&residue| residues[residue].atom_record);
            let last = run
                .iter()
                .rposition(|&residue| residues[residue].atom_record);
            match (first, last) {
                (Some(mut first), Some(mut last)) => {
                    while first > 0 && residues[run[first - 1]].links_to(&residues[run[first]]) {
                        first -= 1;
                    }
                    while last + 1 < run.len()
                        && residues[run[last]].links_to(&residues[run[last + 1]])
                    {
                        last += 1;
                    }
                    &run[first..=last]
                }
                _ => &[],
            }
        };
        let polymer = span
            .iter()
            .enumerate()
            .filter(|&(index, &residue)| {
                let residue = &residues[residue];
                let linked = if terminated {
                    residue.has_linking_atom()
                } else {
                    (index > 0 && residues[span[index - 1]].links_to(residue))
                        || span
                            .get(index + 1)
                            .is_some_and(|&next| residue.links_to(&residues[next]))
                };
                !residue.solvent && (residue.atom_record || residue.polymer_template || linked)
            })
            .map(|(_, &residue)| residue)
            .collect::<Vec<_>>();
        if polymer.len() < 2 {
            return;
        }
        let next = self.next_position.entry(asym_id.to_owned()).or_insert(1);
        for residue in polymer {
            residues[residue].polymer_position = Some(*next);
            *next += 1;
        }
        // A gap keeps separately terminated segments of one chain apart.
        *next += 1;
    }

    fn close_all(&mut self, residues: &mut [Residue]) {
        let chains = self.open.keys().cloned().collect::<Vec<_>>();
        for chain in chains {
            self.close(&chain, false, residues);
        }
    }
}

/// Atom rows of every model, and the issues of their records: elements
/// taken from templates or atom names, charges without a sign, and numbers
/// that read as hexadecimal or go back.
fn read_rows(
    document: &PdbDocument,
    components: Option<&ComponentLibrary>,
    cancellation: Option<&Cancellation>,
) -> Result<(Vec<AtomRow>, Vec<StructureInterpretIssue>), FormatError> {
    let atoms = document.atoms();
    if atoms.is_empty() {
        return Err(structure_error(
            DiagnosticKind::Inconsistent,
            None,
            "The entry has no coordinate records.",
        ));
    }
    let terminated = document
        .ter_records()
        .iter()
        .filter_map(|ter| {
            let index = ter.preceding_atoms().checked_sub(1)?;
            (atoms[index].model() == ter.model()).then_some(index)
        })
        .collect::<BTreeSet<_>>();

    let mut residues = Vec::<Residue>::new();
    let mut residue_of = Vec::with_capacity(atoms.len());
    let mut model = None;
    let mut runs = ChainRuns::default();
    let mut current = None;
    let mut repeats = HashMap::<(String, i32, Option<char>), usize>::new();
    for (index, atom) in atoms.iter().enumerate() {
        poll_every(cancellation, index)?;
        if model != Some(atom.model()) {
            runs.close_all(&mut residues);
            runs = ChainRuns::default();
            current = None;
            repeats.clear();
            model = Some(atom.model());
        }
        let identity = (
            asym_id(atom),
            atom.residue_sequence(),
            atom.insertion_code(),
        );
        let residue = match &current {
            Some((current_identity, residue)) if *current_identity == identity => *residue,
            _ => {
                let repeat = repeats.entry(identity.clone()).or_insert(0);
                let residue = residues.len();
                residues.push(Residue {
                    asym_id: identity.0.clone(),
                    repeat: *repeat,
                    atom_record: false,
                    water: is_water(atom.residue_name()),
                    solvent: is_water(atom.residue_name()) || is_ion(atom.residue_name()),
                    polymer_template: polymer_template(atom.residue_name(), components),
                    heads: Vec::new(),
                    tails: Vec::new(),
                    polymer_position: None,
                });
                *repeat += 1;
                runs.open
                    .entry(identity.0.clone())
                    .or_default()
                    .push(residue);
                current = Some((identity, residue));
                residue
            }
        };
        residues[residue].atom_record |= atom.kind() == PdbAtomKind::Atom;
        if let Some((head, atomic_number)) = linking_role(atom.name()) {
            let element = Element::from_atomic_number(atomic_number).expect("linking element");
            let residue = &mut residues[residue];
            let ends = if head {
                &mut residue.heads
            } else {
                &mut residue.tails
            };
            ends.push((element, atom.position()));
        }
        residue_of.push(residue);
        if terminated.contains(&index) {
            let chain = residues[residue].asym_id.clone();
            runs.close(&chain, true, &mut residues);
            // `TER` also ends the residue: records after it start a new one,
            // even when they repeat its chain, number, and insertion code.
            current = None;
        }
    }
    runs.close_all(&mut residues);

    let conventional = follows_name_alignment(atoms);
    let (mut from_templates, mut from_names, mut unsigned) = (Vec::new(), Vec::new(), Vec::new());
    let mut rows = Vec::with_capacity(atoms.len());
    for (index, atom) in atoms.iter().enumerate() {
        poll_every(cancellation, index)?;
        let residue = &residues[residue_of[index]];
        let line = atom.location().line();
        if atom.charge_is_unsigned() {
            unsigned.push(line);
        }
        let (element, isotope) = match atom.element() {
            Some(symbol) => parse_structure_element_symbol(symbol).ok_or_else(|| {
                structure_error(
                    DiagnosticKind::InvalidValue,
                    Some(line),
                    format!("The element symbol {symbol:?} is unknown."),
                )
            })?,
            None => {
                if let Some(element) =
                    ComponentLibrary::template_element(components, atom.residue_name(), atom.name())
                {
                    from_templates.push(line);
                    (element, None)
                } else {
                    from_names.push(line);
                    (name_element(atom, conventional)?, None)
                }
            }
        };
        let sequence = atom.residue_sequence().to_string();
        let insertion_code = atom.insertion_code().map(String::from);
        let insertion = insertion_code.as_deref().unwrap_or("");
        let residue_key = match residue.repeat {
            0 => format!("auth:{sequence}:{insertion}"),
            repeat => format!("auth:{sequence}:{insertion}:{repeat}"),
        };
        let kind = if residue.polymer_position.is_some() {
            EntityKind::Polymer
        } else if residue.water {
            EntityKind::Water
        } else {
            EntityKind::NonPolymer
        };
        let instance_key = if kind.is_macro() {
            format!("macro:{}", residue.asym_id)
        } else {
            format!("small:{}:{residue_key}", residue.asym_id)
        };
        rows.push(AtomRow {
            line,
            row_index: index,
            site_order: index,
            model_id: model_id(atom.model()),
            entity_id: None,
            kind,
            instance_key,
            label_asym_id: None,
            asym_id: residue.asym_id.clone(),
            auth_asym_id: atom.chain_id().map(String::from),
            residue_key,
            label_seq_id: None,
            polymer_position: residue.polymer_position,
            auth_seq_id: Some(sequence),
            insertion_code,
            occurrence: (residue.repeat > 0).then_some(residue.repeat),
            label_comp_id: None,
            comp_id: atom.residue_name().to_owned(),
            auth_comp_id: Some(atom.residue_name().to_owned()),
            label_atom_name: None,
            atom_name: atom.name().to_owned(),
            auth_atom_name: Some(atom.name().to_owned()),
            atom_site_id: atom.serial().map(|serial| serial.to_string()),
            alt_id: atom.alternate_location().map(String::from),
            occupancy: atom.occupancy(),
            b_factor: atom.temperature_factor(),
            point: Some(atom.position()),
            element,
            isotope,
            formal_charge: atom.formal_charge(),
        });
    }
    let mut issues = Vec::new();
    if !from_templates.is_empty() {
        issues.push(StructureInterpretIssue::ElementsFromTemplates {
            source_lines: from_templates,
        });
    }
    if !from_names.is_empty() {
        issues.push(StructureInterpretIssue::ElementsInferredFromAtomNames {
            source_lines: from_names,
        });
    }
    if !unsigned.is_empty() {
        issues.push(StructureInterpretIssue::UnsignedCharges {
            source_lines: unsigned,
        });
    }
    issues.extend(numbering_issues(document));
    Ok((rows, issues))
}

/// Whether the template of a residue named `name` is a polymer one: the
/// caller's library first, then the bundled one.
fn polymer_template(name: &str, components: Option<&ComponentLibrary>) -> bool {
    components
        .and_then(|library| library.get(name))
        .or_else(|| ComponentLibrary::bundled().get(name))
        .is_some_and(|template| {
            matches!(
                template.linking(),
                ComponentLinking::Peptide | ComponentLinking::Nucleotide
            )
        })
}

/// Whether an atom of this name is a polymer link's head (`N`, `P`), which
/// joins the residue before, or its tail (`C`, `O3'`), which joins the one
/// after, with its atomic number.
fn linking_role(name: &str) -> Option<(bool, u8)> {
    match name {
        "N" => Some((true, 7)),
        "P" => Some((true, 15)),
        "C" => Some((false, 6)),
        "O3'" | "O3*" => Some((false, 8)),
        _ => None,
    }
}

/// Reports where numbering turns hexadecimal and where it goes back: atom
/// serial numbers that decrease within a model, and residue numbers that
/// fall back after 9,999.
fn numbering_issues(document: &PdbDocument) -> Vec<StructureInterpretIssue> {
    let mut issues = document
        .hexadecimal_numbers()
        .iter()
        .map(
            |(field, location)| StructureInterpretIssue::HexadecimalNumbers {
                field: *field,
                source_line: location.line(),
            },
        )
        .collect::<Vec<_>>();
    let mut restarts = [None; 2];
    let mut previous: Option<&PdbAtom> = None;
    for atom in document.atoms() {
        let Some(last) = previous
            .replace(atom)
            .filter(|last| last.model() == atom.model())
        else {
            continue;
        };
        let line = atom.location().line();
        if matches!((last.serial(), atom.serial()), (Some(last), Some(serial)) if serial < last) {
            restarts[0].get_or_insert(line);
        }
        if last.residue_sequence() == 9999 && atom.residue_sequence() < 9999 {
            restarts[1].get_or_insert(line);
        }
    }
    for (field, line) in [NumberField::AtomSerial, NumberField::ResidueNumber]
        .into_iter()
        .zip(restarts)
    {
        if let Some(source_line) = line {
            issues.push(StructureInterpretIssue::NumbersRestart { field, source_line });
        }
    }
    issues
}

/// What an atom name gives as its element.
#[derive(Debug, PartialEq, Eq)]
enum NameElement {
    Element(Element),
    /// A one-letter and a two-letter reading that alignment cannot tell
    /// apart.
    Ambiguous([Element; 2]),
    None,
}

/// The element of an atom no template names, from its name; see
/// [`element_from_atom_name`].
fn name_element(atom: &PdbAtom, conventional: bool) -> Result<Element, FormatError> {
    let line = Some(atom.location().line());
    match element_from_atom_name(atom.name_field(), atom.residue_name(), conventional) {
        NameElement::Element(element) => Ok(element),
        NameElement::Ambiguous([one, two]) => Err(structure_error(
            DiagnosticKind::InvalidValue,
            line,
            format!(
                "Atom {:?} has no element symbol, and its name reads as {} or {}. Add element \
                 symbols in columns 77-78.",
                atom.name(),
                one.symbol(),
                two.symbol()
            ),
        )),
        NameElement::None => Err(structure_error(
            DiagnosticKind::InvalidValue,
            line,
            format!(
                "Atom {:?} has no element symbol, and its name gives none.",
                atom.name()
            ),
        )),
    }
}

/// Whether the entry aligns atom names by the PDB element convention: among
/// atoms without element symbols, some name shorter than four characters
/// starts in column 13 with a two-letter element symbol, as `FE` does, and
/// none starts there without one, as `C1` would. Writers that pad every
/// short name into column 14 show no such name.
fn follows_name_alignment(atoms: &[PdbAtom]) -> bool {
    let mut evidence = false;
    for atom in atoms.iter().filter(|atom| atom.element().is_none()) {
        let starts_in_13 = atom
            .name_field()
            .as_bytes()
            .first()
            .is_some_and(u8::is_ascii_alphabetic);
        if starts_in_13 && atom.name().len() < 4 {
            if two_letter_symbol(atom.name()).is_none() {
                return false;
            }
            evidence = true;
        }
    }
    evidence
}

/// An element symbol without a mass number.
fn element_symbol(text: &str) -> Option<Element> {
    match parse_structure_element_symbol(text) {
        Some((element, None)) => Some(element),
        _ => None,
    }
}

/// The element whose symbol is the first two characters of `text`, when
/// both are letters.
fn two_letter_symbol(text: &str) -> Option<Element> {
    text.get(..2)
        .filter(|two| two.bytes().all(|byte| byte.is_ascii_alphabetic()))
        .and_then(element_symbol)
}

/// The element a PDB atom name implies by convention, for an atom no
/// template names. A one-atom residue named for its element, such as `CA`,
/// is that element. A name starting with a digit, as in `1HG1`, has a
/// one-letter symbol in column 14. A four-character name starting with `H`
/// is a hydrogen. In an entry that follows the
/// [alignment convention](follows_name_alignment), a name starting in
/// column 14 has a one-letter symbol and a shorter name starting in column
/// 13 a two-letter one. Otherwise, and for other four-character names, a
/// name whose first letter and first two letters are both element symbols
/// is ambiguous, as ` CL1` is in an entry that pads every name, unless its
/// first letter is `H`: such a name, as ` HG1`, is a hydrogen.
fn element_from_atom_name(name_field: &str, residue_name: &str, conventional: bool) -> NameElement {
    let name = name_field.trim();
    if name.eq_ignore_ascii_case(residue_name) {
        if let Some(element) = element_symbol(name) {
            return NameElement::Element(element);
        }
    }
    let Some(start) = name_field.find(|character: char| character.is_ascii_alphabetic()) else {
        return NameElement::None;
    };
    let one = element_symbol(&name_field[start..=start]);
    let two = two_letter_symbol(&name_field[start..]);
    let four = start == 0 && name.len() == 4;
    let digit_led = name_field[..start].bytes().any(|byte| byte != b' ');
    let ambiguous = |one: Element, two: Element| {
        if one.atomic_number() == 1 {
            NameElement::Element(one)
        } else {
            NameElement::Ambiguous([one, two])
        }
    };
    match (one, two) {
        // A digit in column 13 leaves column 14 a one-letter symbol.
        (Some(one), _) if digit_led => NameElement::Element(one),
        (Some(one), Some(two)) if start > 0 && !conventional => ambiguous(one, two),
        (Some(one), _) if start > 0 => NameElement::Element(one),
        _ if four && name_field.as_bytes()[0].eq_ignore_ascii_case(&b'H') => {
            NameElement::Element(Element::from_atomic_number(1).expect("hydrogen"))
        }
        (Some(one), Some(two)) if four || !conventional => ambiguous(one, two),
        (_, Some(two)) => NameElement::Element(two),
        (Some(one), None) => NameElement::Element(one),
        (None, None) => NameElement::None,
    }
}

fn connection_records(document: &PdbDocument) -> Result<Vec<ConnectionRecord>, FormatError> {
    let mut records = Vec::new();
    for record in document.ssbond_records() {
        let [left, right] = record.partners();
        records.push(ConnectionRecord::Bond(Box::new(DeclaredBond {
            connection_id: record.serial().map(|serial| serial.to_string()),
            connection_type: "SSBOND".to_owned(),
            line: Some(record.location().line()),
            order: DeclaredOrder::Stated(BondOrder::Single),
            partners: [partner(left, Some("SG")), partner(right, Some("SG"))],
            connectivity_only: false,
            model: None,
        })));
    }
    // LINK records state no bond type; a link between a metal and a
    // non-metal of different residues is coordination.
    for record in document.link_records() {
        let [left, right] = record.partners();
        records.push(ConnectionRecord::Bond(Box::new(DeclaredBond {
            connection_id: None,
            connection_type: "LINK".to_owned(),
            line: Some(record.location().line()),
            order: DeclaredOrder::ByElements(None),
            partners: [partner(left, None), partner(right, None)],
            connectivity_only: false,
            model: None,
        })));
    }
    // Each bond is usually listed from both atoms. Listing a partner twice or
    // three times states a double or triple bond, by convention. Files that
    // repeat the records in every model state each model's bonds; others
    // state them once for all models.
    let models = document
        .conect_records()
        .iter()
        .map(PdbConect::model)
        .collect::<BTreeSet<_>>();
    let per_model = models.len() > 1;
    let mut listed = BTreeMap::<(Option<u32>, (u32, u32)), ([u8; 2], usize)>::new();
    for record in document.conect_records() {
        let serial = record.serial();
        let line = record.location().line();
        for &bonded in record.bonded() {
            if bonded == serial {
                return Err(structure_error(
                    DiagnosticKind::InvalidValue,
                    Some(line),
                    format!("CONECT record bonds atom {serial} to itself"),
                ));
            }
            let model = record.model().filter(|_| per_model);
            let (counts, _) = listed
                .entry((model, (serial.min(bonded), serial.max(bonded))))
                .or_insert(([0, 0], line));
            let direction = usize::from(serial > bonded);
            counts[direction] = counts[direction].saturating_add(1);
        }
    }
    // Atoms the records bond, directly or through others, that list no pair
    // twice, as wwPDB files do, state no orders. Groups are counted apart, so
    // that a molecule a program appends with orders leaves the others alone.
    let mut groups = ConectGroups::of(listed.keys().map(|&(_, pair)| pair));
    let mut stating = BTreeSet::new();
    for (&(_, (left, _)), (counts, _)) in &listed {
        if counts.iter().any(|&count| count > 1) {
            stating.insert(groups.root(left));
        }
    }
    let mut listed = listed.into_iter().collect::<Vec<_>>();
    listed.sort_by_key(|&((model, pair), (_, line))| (line, model, pair));
    for ((model, (left, right)), (counts, line)) in listed {
        let states_orders = stating.contains(&groups.root(left));
        let order = match counts[0].max(counts[1]) {
            _ if !states_orders => None,
            1 => Some(BondOrder::Single),
            2 => Some(BondOrder::Double),
            3 => Some(BondOrder::Triple),
            count => {
                return Err(structure_error(
                    DiagnosticKind::InvalidValue,
                    Some(line),
                    format!(
                        "The CONECT records list the bond between atoms {left} and {right} \
                         {count} times."
                    ),
                ))
            }
        };
        records.push(ConnectionRecord::Bond(Box::new(DeclaredBond {
            connection_id: None,
            connection_type: "CONECT".to_owned(),
            line: Some(line),
            order: DeclaredOrder::ByElements(order),
            partners: [serial_partner(left), serial_partner(right)],
            connectivity_only: true,
            model: model.map(|model| model_id(Some(model))),
        })));
    }
    Ok(records)
}

/// The groups of atoms `CONECT` records bond to each other, directly or
/// through others, by atom serial.
struct ConectGroups {
    parent: BTreeMap<u32, u32>,
}

impl ConectGroups {
    fn of(pairs: impl IntoIterator<Item = (u32, u32)>) -> Self {
        let mut groups = Self {
            parent: BTreeMap::new(),
        };
        for (left, right) in pairs {
            let (left, right) = (groups.root(left), groups.root(right));
            if left != right {
                groups.parent.insert(left.max(right), left.min(right));
            }
        }
        groups
    }

    /// The serial that stands for the group of `serial`.
    fn root(&mut self, serial: u32) -> u32 {
        let mut root = serial;
        while let Some(&parent) = self.parent.get(&root) {
            root = parent;
        }
        // Later lookups go straight to the root.
        let mut current = serial;
        while let Some(parent) = self.parent.insert(current, root) {
            if parent == root {
                break;
            }
            current = parent;
        }
        if current == root {
            self.parent.remove(&root);
        }
        root
    }
}

fn partner(partner: &PdbPartner, atom_name: Option<&str>) -> PartnerSpec {
    if let Some(symmetry) = partner.symmetry() {
        if symmetry != "1555" && symmetry != "1_555" {
            return Err(ConnectionResolutionReason::UnsupportedSymmetry {
                symmetry: symmetry.to_owned(),
            });
        }
    }
    Ok(ConnectionPartnerSelector {
        author: AuthorConnectionSelector {
            asym_id: partner.chain_id().map(String::from),
            component_id: Some(partner.residue_name().to_owned()),
            sequence_id: Some(partner.residue_sequence().to_string()),
            atom_id: atom_name.or(partner.atom_name()).map(str::to_owned),
            insertion_code: partner.insertion_code().map(String::from),
            alternate_location: partner.alternate_location().map(String::from),
        },
        ..ConnectionPartnerSelector::default()
    })
}

fn serial_partner(serial: u32) -> PartnerSpec {
    Ok(ConnectionPartnerSelector {
        site_id: Some(serial.to_string()),
        ..ConnectionPartnerSelector::default()
    })
}

#[cfg(test)]
mod tests;
