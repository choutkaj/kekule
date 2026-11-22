//! Reads a PDB entry into the source-neutral interpretation input.
//!
//! PDB records carry author identities only. A residue is a run of
//! consecutive coordinate records with one chain, sequence number, and
//! insertion code, which a `TER` record also ends; a number that recurs after
//! other residues or a `TER` record names a new residue, as when numbering
//! wraps around or a segment restarts it. Without entity records, a chain's
//! polymer is read from its `TER` records, and its residues' order gives
//! polymer adjacency.
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
use super::connectivity::ConnectivityCatalog;
use super::structure_set::{interpret_conformations, interpret_models};
use super::types::{
    ConnectionResolutionReason, ConnectivityOptions, EntityKind, NumberField,
    StructureInterpretIssue, StructureInterpretOptions, StructureInterpretation,
    StructureSetInterpretOptions, StructureSetInterpretation,
};
use super::{interpret_source, stated_cell, AltLocResidue, StructureSource};
use crate::components::ComponentLibrary;
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
    let (rows, mut issues) = read_rows(
        document,
        connectivity.components().map(AsRef::as_ref),
        cancellation,
    )?;
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
        catalog: ConnectivityCatalog::new(connectivity),
        cell,
        issues,
        polymer_entities: Vec::new(),
        assemblies,
    })
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
    polymer_position: Option<i32>,
}

/// The residues of each chain since its last `TER` record, within one model.
#[derive(Default)]
struct ChainRuns {
    open: BTreeMap<String, Vec<usize>>,
    next_position: BTreeMap<String, i32>,
}

impl ChainRuns {
    /// Ends a chain's run. A run a `TER` record ends is polymer; otherwise its
    /// polymer spans its `ATOM` records. Water and ions are never polymer,
    /// and neither is a single residue.
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
                (Some(first), Some(last)) => &run[first..=last],
                _ => &[],
            }
        };
        let polymer = span
            .iter()
            .copied()
            .filter(|&residue| !residues[residue].solvent)
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
            order: DeclaredOrder::ByElements(BondOrder::Single),
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
    let mut listed = listed.into_iter().collect::<Vec<_>>();
    listed.sort_by_key(|&((model, pair), (_, line))| (line, model, pair));
    for ((model, (left, right)), (counts, line)) in listed {
        let order = match counts[0].max(counts[1]) {
            1 => BondOrder::Single,
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
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
