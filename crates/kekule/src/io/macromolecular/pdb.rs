//! Reads a PDB entry into the source-neutral interpretation input.
//!
//! PDB records carry author identities only. A residue is a run of
//! consecutive coordinate records with one chain, sequence number, and
//! insertion code, which a `TER` record also ends; a number that recurs after
//! other residues or a `TER` record names a new residue, as when numbering
//! wraps around or a segment restarts it. Without entity records, a chain's
//! polymer is read from its `TER` records, and its residues' order gives
//! polymer adjacency. A residue of `HETATM` records in it is polymer only
//! when its template is a polymer one or one of its linking atoms (`N`, `C`,
//! `P`, `O3'`, or `O3*`) lies within bonding distance of one of a
//! neighbour's, as for a capping group or a 5'-terminal nucleotide
//! triphosphate. A ligand listed between two residues then does not separate
//! them. Such residues bonded to the first or last `ATOM` residue of a chain
//! without `TER` records extend its polymer.
use std::collections::{BTreeMap, BTreeSet, HashMap};

use super::types::{poll, poll_every, structure_error};
use crate::diagnostics::{Cancellation, DiagnosticKind, FormatError};
use crate::molecule::{BondOrder, Element};

use super::super::element_symbol::parse_structure_element_symbol;
use super::super::{PdbAtom, PdbAtomKind, PdbDocument, PdbPartner};
use super::altloc::inventory;
use super::assemblies::read_pdb_assemblies;
use super::atom_site::{is_water, AtomRow};
use super::connections::{
    AuthorConnectionSelector, ConnectionPartnerSelector, ConnectionRecord, DeclaredBond,
    DeclaredOrder, PartnerSpec,
};
use super::connectivity::{within_bonding_distance, ConnectivityCatalog};
use super::structure_set::{interpret_conformations, interpret_models};
use super::types::{
    ConnectionResolutionReason, ConnectivityOptions, EntityKind, StructureInterpretIssue,
    StructureInterpretOptions, StructureInterpretation, StructureSetInterpretOptions,
    StructureSetInterpretation,
};
use super::{interpret_source, stated_cell, AltLocResidue, StructureSource};
use crate::components::{ComponentLibrary, ComponentLinking};
use crate::math::Point3;

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
    let (rows, inferred) = read_rows(
        document,
        connectivity.components().map(|library| library.as_ref()),
        cancellation,
    )?;
    poll(cancellation)?;
    let mut polymer_order = BTreeMap::new();
    for row in rows.iter().filter(|row| row.polymer_position.is_some()) {
        let next = polymer_order.len();
        polymer_order.entry(row.asym_id.clone()).or_insert(next);
    }
    let mut issues = Vec::new();
    if !inferred.is_empty() {
        issues.push(StructureInterpretIssue::ElementsInferredFromAtomNames {
            source_lines: inferred,
        });
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
    /// Whether its template is a polymer one.
    polymer_template: bool,
    /// Its atoms that a polymer link joins to the residue before it (`N`,
    /// `P`) and to the one after it (`C`, `O3'`).
    heads: Vec<(Element, Point3)>,
    tails: Vec<(Element, Point3)>,
    polymer_position: Option<i32>,
}

impl Residue {
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
    /// ends. Water is never polymer, nor is a residue of `HETATM` records
    /// with neither a polymer template nor a link to a neighbour, nor a
    /// single residue.
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
                let linked = (index > 0 && residues[span[index - 1]].links_to(residue))
                    || span
                        .get(index + 1)
                        .is_some_and(|&next| residue.links_to(&residues[next]));
                !residue.water && (residue.atom_record || residue.polymer_template || linked)
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

/// Atom rows of every model, and the lines of records whose element was
/// inferred from the atom name.
fn read_rows(
    document: &PdbDocument,
    components: Option<&ComponentLibrary>,
    cancellation: Option<&Cancellation>,
) -> Result<(Vec<AtomRow>, Vec<usize>), FormatError> {
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

    let mut inferred = Vec::new();
    let mut rows = Vec::with_capacity(atoms.len());
    for (index, atom) in atoms.iter().enumerate() {
        poll_every(cancellation, index)?;
        let residue = &residues[residue_of[index]];
        let line = atom.location().line();
        let (element, isotope) = match atom.element() {
            Some(symbol) => parse_structure_element_symbol(symbol).ok_or_else(|| {
                structure_error(
                    DiagnosticKind::InvalidValue,
                    Some(line),
                    format!("The element symbol {symbol:?} is unknown."),
                )
            })?,
            None => {
                inferred.push(line);
                let element = element_from_atom_name(atom.name_field(), atom.residue_name())
                    .ok_or_else(|| {
                        structure_error(
                            DiagnosticKind::InvalidValue,
                            Some(line),
                            format!(
                                "Atom {:?} has no element symbol, and its name gives none.",
                                atom.name()
                            ),
                        )
                    })?;
                (element, None)
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
    Ok((rows, inferred))
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

/// The element a PDB atom name implies by convention: two-letter symbols
/// start in column 13, one-letter symbols in column 14, and four-character
/// hydrogen names fill columns 13-16. A one-atom residue named for its
/// element, such as `ZN`, is that element.
fn element_from_atom_name(name_field: &str, residue_name: &str) -> Option<Element> {
    let symbol = |text: &str| match parse_structure_element_symbol(text) {
        Some((element, None)) => Some(element),
        _ => None,
    };
    let name = name_field.trim();
    if name.eq_ignore_ascii_case(residue_name) {
        if let Some(element) = symbol(name) {
            return Some(element);
        }
    }
    let bytes = name_field.as_bytes();
    let first = *bytes.first()?;
    if first.is_ascii_alphabetic() {
        let hydrogen_name = first.eq_ignore_ascii_case(&b'H') && name.len() == 4;
        if !hydrogen_name && bytes.get(1).is_some_and(u8::is_ascii_alphabetic) {
            if let Some(element) = symbol(&name_field[..2]) {
                return Some(element);
            }
        }
        return symbol(&name_field[..1]);
    }
    let index = name_field[1..].find(|character: char| character.is_ascii_alphabetic())? + 1;
    symbol(&name_field[index..=index])
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
        })));
    }
    // LINK records state no bond type; links to a metal are coordination.
    for record in document.link_records() {
        let [left, right] = record.partners();
        records.push(ConnectionRecord::Bond(Box::new(DeclaredBond {
            connection_id: None,
            connection_type: "LINK".to_owned(),
            line: Some(record.location().line()),
            order: DeclaredOrder::ByElements(None),
            partners: [partner(left, None), partner(right, None)],
            connectivity_only: false,
        })));
    }
    // Each bond is usually listed from both atoms. Listing a partner twice or
    // three times states a double or triple bond, by convention, as Maestro
    // and Kekule write.
    let mut listed = BTreeMap::<(u32, u32), ([u8; 2], usize)>::new();
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
            let (counts, _) = listed
                .entry((serial.min(bonded), serial.max(bonded)))
                .or_insert(([0, 0], line));
            let direction = usize::from(serial > bonded);
            counts[direction] = counts[direction].saturating_add(1);
        }
    }
    // Atoms the records bond, directly or through others, that list no pair
    // twice, as wwPDB files do, state no orders. Groups are counted apart, so
    // that a molecule a program appends with orders leaves the others alone.
    let mut groups = ConectGroups::of(listed.keys().copied());
    let mut stating = BTreeSet::new();
    for (&(left, _), (counts, _)) in &listed {
        if counts.iter().any(|&count| count > 1) {
            stating.insert(groups.root(left));
        }
    }
    let mut listed = listed.into_iter().collect::<Vec<_>>();
    listed.sort_by_key(|&(pair, (_, line))| (line, pair));
    for ((left, right), (counts, line)) in listed {
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
                    "The CONECT records list the bond between atoms {left} and {right} {count} times."
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
