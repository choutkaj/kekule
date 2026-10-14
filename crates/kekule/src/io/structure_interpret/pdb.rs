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

use crate::core::{BondOrder, Element};
use crate::diagnostics::DiagnosticKind;

use super::super::structure_connectivity::ConnectivityCatalog;
use super::super::{PdbAtom, PdbAtomKind, PdbDocument, PdbPartner};
use super::altloc::inventory;
use super::atom_site::{canonical_mmcif_element_symbol, is_water, AtomRow};
use super::connections::{
    AuthorConnectionSelector, ConnectionPartnerSelector, ConnectionRecord, DeclaredBond,
    DeclaredOrder, PartnerSpec,
};
use super::ensemble::{interpret_conformations, interpret_ensemble};
use super::types::{
    ConnectionResolutionReason, ConnectivityOptions, EntityKind, StructureEnsembleInterpretError,
    StructureEnsembleInterpretOptions, StructureEnsembleInterpretation, StructureInterpretError,
    StructureInterpretIssue, StructureInterpretOptions, StructureInterpretation,
};
use super::{interpret_source, AltLocResidue, StructureSource};

pub(crate) fn interpret_pdb(
    document: &PdbDocument,
    options: StructureInterpretOptions,
) -> Result<StructureInterpretation, StructureInterpretError> {
    let source = pdb_source(
        document,
        options.strict_entity_metadata,
        &options.connectivity,
    )?;
    interpret_source(source, &options)
}

/// Interprets explicitly selected or all coordinate models as one
/// shared-topology non-temporal ensemble.
pub(crate) fn interpret_pdb_ensemble(
    document: &PdbDocument,
    options: StructureEnsembleInterpretOptions,
) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
    let mut seen = BTreeSet::new();
    let available = document
        .atoms()
        .iter()
        .map(PdbAtom::model)
        .filter(|model| seen.insert(*model))
        .map(model_id)
        .collect();
    interpret_ensemble(available, options, |options| {
        pdb_source(
            document,
            options.strict_entity_metadata,
            &options.connectivity,
        )
    })
}

/// Explicit caller-selected alternate conformations, never automatic altloc
/// enumeration.
pub(crate) fn interpret_pdb_conformations(
    document: &PdbDocument,
    selections: &[StructureInterpretOptions],
) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
    interpret_conformations(selections, |options| {
        interpret_pdb(document, options.clone())
    })
}

pub(crate) fn pdb_alternate_locations(
    document: &PdbDocument,
) -> Result<Vec<AltLocResidue>, StructureInterpretError> {
    inventory(&read_rows(document)?.0)
}

fn pdb_source(
    document: &PdbDocument,
    strict_entity_metadata: bool,
    connectivity: &ConnectivityOptions,
) -> Result<StructureSource, StructureInterpretError> {
    if strict_entity_metadata {
        return Err(StructureInterpretError::new(
            DiagnosticKind::Inconsistent,
            None,
            "PDB entries declare no entity types; strict entity metadata cannot be satisfied",
        ));
    }
    let (rows, inferred) = read_rows(document)?;
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
    Ok(StructureSource {
        name: document.id_code().unwrap_or_default().to_owned(),
        entity_definitions: 0,
        rows,
        conformations: BTreeMap::new(),
        polymer_order,
        connections: connection_records(document)?,
        catalog: ConnectivityCatalog::new(connectivity),
        issues,
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
    /// polymer spans its `ATOM` records. Water is never polymer, and neither
    /// is a single residue.
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
            .filter(|&residue| !residues[residue].water)
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
) -> Result<(Vec<AtomRow>, Vec<usize>), StructureInterpretError> {
    let atoms = document.atoms();
    if atoms.is_empty() {
        return Err(StructureInterpretError::new(
            DiagnosticKind::Inconsistent,
            None,
            "PDB entry has no coordinate records",
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

    let mut inferred = Vec::new();
    let mut rows = Vec::with_capacity(atoms.len());
    for (index, atom) in atoms.iter().enumerate() {
        let residue = &residues[residue_of[index]];
        let line = atom.location().line();
        let element = match atom.element() {
            Some(symbol) => Element::from_symbol(&canonical_mmcif_element_symbol(symbol))
                .ok_or_else(|| {
                    StructureInterpretError::new(
                        DiagnosticKind::InvalidValue,
                        Some(line),
                        format!("unknown element `{symbol}`"),
                    )
                })?,
            None => {
                inferred.push(line);
                element_from_atom_name(atom.name_field(), atom.residue_name()).ok_or_else(|| {
                    StructureInterpretError::new(
                        DiagnosticKind::InvalidValue,
                        Some(line),
                        format!(
                            "atom `{}` has no element symbol, and its name implies none",
                            atom.name()
                        ),
                    )
                })?
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
            formal_charge: atom.formal_charge(),
        });
    }
    Ok((rows, inferred))
}

/// The element a PDB atom name implies by convention: two-letter symbols
/// start in column 13, one-letter symbols in column 14, and four-character
/// hydrogen names fill columns 13-16. A one-atom residue named for its
/// element, such as `ZN`, is that element.
fn element_from_atom_name(name_field: &str, residue_name: &str) -> Option<Element> {
    let symbol = |text: &str| Element::from_symbol(&canonical_mmcif_element_symbol(text));
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

fn connection_records(
    document: &PdbDocument,
) -> Result<Vec<ConnectionRecord>, StructureInterpretError> {
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
            order: DeclaredOrder::ByElements(BondOrder::Single),
            partners: [partner(left, None), partner(right, None)],
            connectivity_only: false,
        })));
    }
    // Each bond is usually listed from both atoms. Listing a partner twice or
    // three times states a double or triple bond, by convention.
    let mut listed = BTreeMap::<(u32, u32), ([u8; 2], usize)>::new();
    for record in document.conect_records() {
        let serial = record.serial();
        let line = record.location().line();
        for &bonded in record.bonded() {
            if bonded == serial {
                return Err(StructureInterpretError::new(
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
    let mut listed = listed.into_iter().collect::<Vec<_>>();
    listed.sort_by_key(|&(pair, (_, line))| (line, pair));
    for ((left, right), (counts, line)) in listed {
        let order = match counts[0].max(counts[1]) {
            1 => BondOrder::Single,
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            count => {
                return Err(StructureInterpretError::new(
                    DiagnosticKind::InvalidValue,
                    Some(line),
                    format!(
                    "CONECT records list the bond between atoms {left} and {right} {count} times"
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
mod tests {
    use super::*;

    #[test]
    fn atom_name_alignment_implies_the_element() {
        for (name, residue, symbol) in [
            (" CA ", "ALA", "C"),
            ("CA  ", "CA", "Ca"),
            ("FE  ", "HEM", "Fe"),
            (" N  ", "ALA", "N"),
            ("1HB ", "ALA", "H"),
            ("HG21", "THR", "H"),
            ("HG  ", "HG", "Hg"),
            (" HG ", "CYS", "H"),
            ("ZN  ", "ZN", "Zn"),
            ("CL1 ", "LIG", "Cl"),
            ("OXT ", "ALA", "O"),
            ("SE  ", "MSE", "Se"),
        ] {
            assert_eq!(
                element_from_atom_name(name, residue).map(Element::symbol),
                Some(symbol),
                "{name:?} in {residue}"
            );
        }
        assert_eq!(element_from_atom_name(" 1  ", "LIG"), None);
    }
}
