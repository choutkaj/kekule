//! Encoding of structure interpretation reports.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::Write;

use super::{Encoder, Result};
use crate::geometry::Vector3;
use crate::io::{
    AltLocDecision, AltLocSelectionReason, Assembly, AssemblyOmissionReason, AtomProvenance,
    BondEvidence, CellOmissionReason, ConflictingSelector, ConnectionResolutionReason, EntityKind,
    PolymerEntity, PolymerType, SequenceSource, SourceResidueId, SourceResiduePosition,
    StructureInterpretIssue, StructureInterpretationReport,
};
use crate::sequence::CodeOrigin;
use crate::topology::InstanceAtomId;

impl<W: Write> Encoder<'_, W> {
    /// Unsigned values in the narrowest of 1, 2, 4, or 8 bytes that holds
    /// the largest.
    pub(super) fn packed<I: Iterator<Item = u64>>(&mut self, values: impl Fn() -> I) -> Result {
        let width = match values().max().unwrap_or(0) {
            0..=0xFF => 1,
            0x100..=0xFFFF => 2,
            0x1_0000..=0xFFFF_FFFF => 4,
            _ => 8,
        };
        self.u8(width)?;
        let width = usize::from(width);
        let mut chunk = Vec::with_capacity(PACKED_CHUNK + 8);
        for value in values() {
            chunk.extend_from_slice(&value.to_le_bytes()[..width]);
            if chunk.len() >= PACKED_CHUNK {
                self.bytes(&chunk)?;
                chunk.clear();
            }
        }
        self.bytes(&chunk)
    }

    fn usize(&mut self, value: usize) -> Result {
        self.u64(value as u64)
    }

    fn opt_usize(&mut self, value: Option<usize>) -> Result {
        match value {
            None => self.u8(0),
            Some(value) => {
                self.u8(1)?;
                self.usize(value)
            }
        }
    }

    /// A counted list of unsigned values.
    fn usizes(&mut self, values: &[usize]) -> Result {
        self.count(values.len())?;
        self.packed(|| values.iter().map(|&value| value as u64))
    }

    fn str_list(&mut self, values: &[String]) -> Result {
        self.count(values.len())?;
        for value in values {
            self.str(value)?;
        }
        Ok(())
    }

    /// A string column through a table of its distinct values, in order of
    /// first appearance, and each row's index into it.
    fn dictionary<'s>(&mut self, values: impl Iterator<Item = Option<&'s str>>) -> Result {
        let mut column = Dictionary::default();
        for value in values {
            column.push(value);
        }
        self.dictionary_column(&column)
    }

    fn dictionary_column(&mut self, column: &Dictionary<'_>) -> Result {
        self.count(column.table.len())?;
        self.strings(|| column.table.iter().copied())?;
        self.packed(|| column.indices.iter().map(|&index| u64::from(index)))
    }

    /// Atom references of one instance each, as local atom IDs when every
    /// reference is to its instance and as full IDs otherwise.
    fn atom_refs(&mut self, refs: &[(InstanceAtomId, u32)]) -> Result {
        if refs
            .iter()
            .all(|(atom, instance)| atom.instance().raw() == *instance)
        {
            self.u8(0)?;
            self.packed(|| refs.iter().map(|(atom, _)| u64::from(atom.atom().raw())))
        } else {
            self.u8(1)?;
            self.packed(|| refs.iter().map(|(atom, _)| atom.to_bits()))
        }
    }

    /// A column of entity kinds: every kind's code, then the values of the
    /// other kinds.
    fn entity_kinds(&mut self, kinds: &[&EntityKind]) -> Result {
        self.packed(|| kinds.iter().map(|kind| u64::from(entity_kind_code(kind))))?;
        self.dictionary(kinds.iter().filter_map(|kind| match kind {
            EntityKind::Other(value) => Some(Some(value.as_str())),
            _ => None,
        }))
    }

    fn entity_kind(&mut self, kind: &EntityKind) -> Result {
        self.u8(entity_kind_code(kind))?;
        if let EntityKind::Other(value) = kind {
            self.str(value)?;
        }
        Ok(())
    }

    pub(super) fn report(&mut self, report: &StructureInterpretationReport) -> Result {
        self.str(&report.block_name)?;
        self.opt_str(report.selected_model.as_deref())?;
        for value in [
            report.entity_definitions,
            report.coordinate_models,
            report.macromolecules,
            report.small_molecules,
            report.solvent_molecules,
            report.applied_connections,
            report.connectivity_candidates,
        ] {
            self.usize(value)?;
        }
        self.str_list(&report.ignored_coordinate_models)?;
        self.instance_provenance(report)?;
        self.count(report.coordination.len())?;
        for record in &report.coordination {
            self.u64(record.atoms[0].to_bits())?;
            self.u64(record.atoms[1].to_bits())?;
            self.opt_str(record.connection_id.as_deref())?;
            self.str(&record.connection_type)?;
            self.opt_usize(record.source_line)?;
        }
        self.count(report.alternate_locations.len())?;
        for decision in &report.alternate_locations {
            self.alternate_location(decision)?;
        }
        self.count(report.issues.len())?;
        for issue in &report.issues {
            self.issue(issue)?;
        }
        self.count(report.polymer_entities.len())?;
        for entity in &report.polymer_entities {
            self.polymer_entity(entity)?;
        }
        self.count(report.assemblies.len())?;
        for assembly in &report.assemblies {
            self.assembly(assembly)?;
        }
        Ok(())
    }

    fn assembly(&mut self, assembly: &Assembly) -> Result {
        self.str(&assembly.id)?;
        self.opt_str(assembly.details.as_deref())?;
        self.opt_str(assembly.oligomeric_details.as_deref())?;
        self.opt_usize(assembly.oligomeric_count)?;
        self.count(assembly.generators.len())?;
        for generator in &assembly.generators {
            self.str_list(&generator.chains)?;
            self.count(generator.operations.len())?;
            for operation in &generator.operations {
                self.str_list(&operation.operators)?;
                let [a, b, c] = operation.rotation.columns();
                for vector in [a, b, c, operation.translation] {
                    let Vector3 { x, y, z } = vector;
                    for value in [x, y, z] {
                        self.f64(value)?;
                    }
                }
            }
        }
        Ok(())
    }

    fn instance_provenance(&mut self, report: &StructureInterpretationReport) -> Result {
        let instances = &report.instances;
        self.count(instances.len())?;
        self.packed(|| {
            instances
                .iter()
                .map(|instance| u64::from(instance.molecule.raw()))
        })?;
        self.dictionary(
            instances
                .iter()
                .map(|instance| Some(instance.coordinate_model_id.as_str())),
        )?;
        for lists in [
            instances
                .iter()
                .map(|instance| &instance.asym_ids)
                .collect::<Vec<_>>(),
            instances
                .iter()
                .map(|instance| &instance.entity_ids)
                .collect(),
        ] {
            self.packed(|| lists.iter().map(|list| list.len() as u64))?;
            self.dictionary(
                lists
                    .iter()
                    .flat_map(|list| list.iter().map(|value| Some(value.as_str()))),
            )?;
        }
        self.packed(|| {
            instances
                .iter()
                .map(|instance| instance.entity_kinds.len() as u64)
        })?;
        let kinds = instances
            .iter()
            .flat_map(|instance| &instance.entity_kinds)
            .collect::<Vec<_>>();
        self.entity_kinds(&kinds)?;
        self.packed(|| instances.iter().map(|instance| instance.atoms.len() as u64))?;
        self.packed(|| instances.iter().map(|instance| instance.bonds.len() as u64))?;

        // One pass over the atoms gathers every column, since the strings of
        // one atom lie together in memory.
        let atom_count = instances.iter().map(|instance| instance.atoms.len()).sum();
        let mut refs = Vec::with_capacity(atom_count);
        let mut lines = Vec::with_capacity(atom_count);
        let mut orders = Vec::with_capacity(atom_count);
        let mut kinds = Vec::with_capacity(atom_count);
        let mut sequences = Vec::with_capacity(atom_count);
        let mut occurrences = Vec::with_capacity(atom_count);
        let mut texts: [Dictionary<'_>; ATOM_TEXTS] = Default::default();
        for instance in instances {
            for atom in &instance.atoms {
                refs.push((atom.atom, instance.molecule.raw()));
                lines.push(atom.source_line as u64);
                orders.push(atom.source_order as u64);
                kinds.push(&atom.entity_kind);
                sequences.push(atom.label_sequence_id);
                occurrences.push(atom.occurrence);
                for (column, value) in texts.iter_mut().zip(atom_texts(atom)) {
                    column.push(value);
                }
            }
        }
        self.atom_refs(&refs)?;
        self.packed(|| lines.iter().copied())?;
        self.packed(|| orders.iter().copied())?;
        self.strings(|| {
            instances
                .iter()
                .flat_map(|instance| &instance.atoms)
                .map(|atom| atom.atom_site_id.as_deref())
        })?;
        for column in &texts {
            self.dictionary_column(column)?;
        }
        self.entity_kinds(&kinds)?;
        self.bits(sequences.iter().map(Option::is_some))?;
        for &sequence in sequences.iter().flatten() {
            self.i32(sequence)?;
        }
        self.bits(occurrences.iter().map(Option::is_some))?;
        self.packed(|| occurrences.iter().flatten().map(|&value| value as u64))?;

        let bonds = instances
            .iter()
            .flat_map(|instance| {
                instance
                    .bonds
                    .iter()
                    .map(move |bond| (bond, instance.molecule.raw()))
            })
            .collect::<Vec<_>>();
        let refs = bonds
            .iter()
            .flat_map(|(bond, instance)| bond.atoms.map(|atom| (atom, *instance)))
            .collect::<Vec<_>>();
        self.atom_refs(&refs)?;
        let evidence = bonds
            .iter()
            .map(|(bond, _)| &bond.evidence)
            .collect::<Vec<_>>();
        self.packed(|| {
            evidence
                .iter()
                .map(|evidence| u64::from(evidence_code(evidence)))
        })?;
        let templates = || {
            evidence.iter().filter_map(|evidence| match evidence {
                BondEvidence::Template {
                    component_id,
                    variant,
                } => Some((component_id.as_str(), variant.as_deref())),
                _ => None,
            })
        };
        self.dictionary(templates().map(|(component, _)| Some(component)))?;
        self.dictionary(templates().map(|(_, variant)| variant))?;
        let connections = || {
            evidence.iter().filter_map(|evidence| match evidence {
                BondEvidence::Connection {
                    connection_id,
                    connection_type,
                } => Some((connection_id.as_deref(), connection_type.as_str())),
                _ => None,
            })
        };
        self.dictionary(connections().map(|(id, _)| id))?;
        self.dictionary(connections().map(|(_, kind)| Some(kind)))
    }

    fn residue(&mut self, residue: &SourceResidueId) -> Result {
        self.str(&residue.model_id)?;
        self.str(&residue.asym_id)?;
        match &residue.position {
            SourceResiduePosition::Label(sequence) => {
                self.u8(0)?;
                self.i32(*sequence)?;
            }
            SourceResiduePosition::Author(sequence) => {
                self.u8(1)?;
                self.str(sequence)?;
            }
            SourceResiduePosition::RepeatedAuthor { sequence, repeat } => {
                self.u8(2)?;
                self.str(sequence)?;
                self.usize(*repeat)?;
            }
            SourceResiduePosition::Occurrence {
                component_id,
                index,
            } => {
                self.u8(3)?;
                self.str(component_id)?;
                self.usize(*index)?;
            }
        }
        self.opt_str(residue.insertion_code.as_deref())
    }

    fn alternate_location(&mut self, decision: &AltLocDecision) -> Result {
        self.residue(&decision.residue)?;
        self.str_list(&decision.available_labels)?;
        self.str_list(&decision.selected_labels)?;
        match &decision.reason {
            AltLocSelectionReason::HighestOccupancy => self.u8(0)?,
            AltLocSelectionReason::ExactSelection => self.u8(1)?,
            AltLocSelectionReason::PreferredLabel => self.u8(2)?,
            AltLocSelectionReason::PreferredLabelUnavailable => self.u8(3)?,
            AltLocSelectionReason::SourceConformation { id } => {
                self.u8(4)?;
                self.str(id)?;
            }
        }
        self.usizes(&decision.selected_source_lines)?;
        self.usizes(&decision.omitted_source_lines)
    }

    fn charged_atoms(&mut self, atoms: &[(String, i8)]) -> Result {
        self.count(atoms.len())?;
        for (name, charge) in atoms {
            self.str(name)?;
            self.bytes(&charge.to_le_bytes())?;
        }
        Ok(())
    }

    fn connection_reason(&mut self, reason: &ConnectionResolutionReason) -> Result {
        match reason {
            ConnectionResolutionReason::UnsupportedSymmetry { symmetry } => {
                self.u8(0)?;
                self.str(symmetry)
            }
            ConnectionResolutionReason::MissingSelector => self.u8(1),
            ConnectionResolutionReason::ConflictingSelectorValues { selector } => {
                self.u8(2)?;
                self.u8(match selector {
                    ConflictingSelector::LabelAlternateLocationAliases => 0,
                    ConflictingSelector::LabelAndAuthorAlternateLocations => 1,
                })
            }
            ConnectionResolutionReason::ConflictingLabelAndAuthorSelectors => self.u8(3),
            ConnectionResolutionReason::AlternateLocationOmitted { alternate_location } => {
                self.u8(4)?;
                self.str(alternate_location)
            }
            ConnectionResolutionReason::NoMatchingAtom => self.u8(5),
            ConnectionResolutionReason::MultipleMatchingAtoms => self.u8(6),
        }
    }

    fn issue(&mut self, issue: &StructureInterpretIssue) -> Result {
        use StructureInterpretIssue as Issue;
        match issue {
            Issue::EntityTypeInferred { asym_id, kind } => {
                self.u8(0)?;
                self.str(asym_id)?;
                self.entity_kind(kind)
            }
            Issue::ElementsInferredFromAtomNames { source_lines } => {
                self.u8(1)?;
                self.usizes(source_lines)
            }
            Issue::ConnectionsIgnored {
                connection_type,
                count,
            } => {
                self.u8(2)?;
                self.str(connection_type)?;
                self.usize(*count)
            }
            Issue::ConnectionUnresolved {
                connection_id,
                connection_type,
                partner,
                source_line,
                reason,
            } => {
                self.u8(3)?;
                self.opt_str(connection_id.as_deref())?;
                self.str(connection_type)?;
                self.u8(*partner)?;
                self.opt_usize(*source_line)?;
                self.connection_reason(reason)
            }
            Issue::ConnectionAmbiguous {
                connection_id,
                connection_type,
                partner,
                source_line,
                candidates,
                reason,
            } => {
                self.u8(4)?;
                self.opt_str(connection_id.as_deref())?;
                self.str(connection_type)?;
                self.u8(*partner)?;
                self.opt_usize(*source_line)?;
                self.usize(*candidates)?;
                self.connection_reason(reason)
            }
            Issue::CoordinateModelIgnored {
                model_id,
                atom_site_rows,
            } => {
                self.u8(5)?;
                self.str(model_id)?;
                self.usize(*atom_site_rows)
            }
            Issue::AlternateLocationOmitted {
                atom_name,
                alt_id,
                residue,
                source_line,
            } => {
                self.u8(6)?;
                self.str(atom_name)?;
                self.opt_str(alt_id.as_deref())?;
                self.residue(residue)?;
                self.usize(*source_line)
            }
            Issue::AlternateLocationIncomplete {
                residue,
                alternative,
                missing_atoms,
            } => {
                self.u8(7)?;
                self.residue(residue)?;
                self.str(alternative)?;
                self.str_list(missing_atoms)
            }
            Issue::AlternateLocationOccupancyMissing {
                residue,
                source_lines,
            } => {
                self.u8(8)?;
                self.residue(residue)?;
                self.usizes(source_lines)
            }
            Issue::AlternateLocationCombinationRejected {
                residue,
                alternative,
                reason,
            } => {
                self.u8(9)?;
                self.residue(residue)?;
                self.str(alternative)?;
                self.str(reason)
            }
            Issue::AlternateLocationOccupancyTie {
                residues,
                alternatives,
            } => {
                self.u8(10)?;
                self.count(residues.len())?;
                for residue in residues {
                    self.residue(residue)?;
                }
                self.str_list(alternatives)
            }
            Issue::ConnectivityCandidatesInferred {
                atom_count,
                candidate_count,
            } => {
                self.u8(11)?;
                self.usize(*atom_count)?;
                self.usize(*candidate_count)
            }
            Issue::ResidueBondsInferred {
                residue,
                component_id,
                bonds,
            } => {
                self.u8(12)?;
                self.residue(residue)?;
                self.str(component_id)?;
                self.usize(*bonds)
            }
            Issue::HydrogenAttachedByDistance {
                residue,
                atom_name,
                heavy_atom_name,
            } => {
                self.u8(13)?;
                self.residue(residue)?;
                self.str(atom_name)?;
                self.str(heavy_atom_name)
            }
            Issue::ChainBreak { residue, next } => {
                self.u8(14)?;
                self.residue(residue)?;
                self.residue(next)
            }
            Issue::ChargesFromBonds {
                residue,
                component_id,
                charged_atoms,
            } => {
                self.u8(15)?;
                self.residue(residue)?;
                self.str(component_id)?;
                self.charged_atoms(charged_atoms)
            }
            Issue::ProtonationFromHydrogens {
                residue,
                component_id,
                variant,
                charged_atoms,
            } => {
                self.u8(16)?;
                self.residue(residue)?;
                self.str(component_id)?;
                self.opt_str(variant.as_deref())?;
                self.charged_atoms(charged_atoms)
            }
            Issue::ResidueFragmented { residue, fragments } => {
                self.u8(17)?;
                self.residue(residue)?;
                self.usize(*fragments)
            }
            Issue::CellOmitted { reason, location } => {
                self.u8(18)?;
                match reason {
                    CellOmissionReason::Incomplete => self.u8(0)?,
                    CellOmissionReason::SpaceGroup { name } => {
                        self.u8(1)?;
                        self.str(name)?;
                    }
                    CellOmissionReason::Placeholder => self.u8(2)?,
                }
                self.usize(location.line())?;
                self.opt_usize(location.column())
            }
            Issue::AssemblyOmitted {
                assembly,
                reason,
                line,
            } => {
                self.u8(19)?;
                self.opt_str(assembly.as_deref())?;
                match reason {
                    AssemblyOmissionReason::Unidentified => self.u8(0)?,
                    AssemblyOmissionReason::NoGenerators => self.u8(1)?,
                    AssemblyOmissionReason::Expression { expression } => {
                        self.u8(2)?;
                        self.str(expression)?;
                    }
                    AssemblyOmissionReason::UnknownOperator { operator } => {
                        self.u8(3)?;
                        self.str(operator)?;
                    }
                    AssemblyOmissionReason::UnreadableOperator { operator } => {
                        self.u8(4)?;
                        self.str(operator)?;
                    }
                    AssemblyOmissionReason::NotRotation { operator } => {
                        self.u8(5)?;
                        self.str(operator)?;
                    }
                }
                self.opt_usize(*line)
            }
        }
    }

    fn polymer_entity(&mut self, entity: &PolymerEntity) -> Result {
        self.str(&entity.id)?;
        match &entity.polymer_type {
            None => self.u8(0)?,
            Some(polymer_type) => {
                self.u8(1)?;
                self.u8(polymer_type_code(polymer_type))?;
                if let PolymerType::Other(value) = polymer_type {
                    self.str(value)?;
                }
            }
        }
        self.u8(sequence_source_code(entity.monomers_source))?;
        self.str(&entity.sequence)?;
        self.u8(sequence_source_code(entity.sequence_source))?;
        self.count(entity.monomers.len())?;
        for monomer in &entity.monomers {
            self.i32(monomer.number)?;
            self.str(&monomer.component_id)?;
            self.str_list(&monomer.alternatives)?;
            self.str(&monomer.code.code)?;
            self.u8(match monomer.code.origin {
                CodeOrigin::Component => 0,
                CodeOrigin::Parent => 1,
                CodeOrigin::Sequence => 2,
                CodeOrigin::Unknown => 3,
            })?;
            self.count(monomer.code.parents.len())?;
            for parent in &monomer.code.parents {
                self.str(parent)?;
            }
        }
        self.count(entity.instances.len())?;
        for instance in &entity.instances {
            self.str(&instance.asym_id)?;
            self.opt_str(instance.author_chain_id.as_deref())?;
            self.opt_u32(instance.chain.map(|chain| chain.raw()))?;
            self.count(instance.positions.len())?;
            for position in &instance.positions {
                self.opt_u32(position.residue.map(|residue| residue.raw()))?;
                self.opt_str(position.author_seq_id.as_deref())?;
                self.opt_str(position.insertion_code.as_deref())?;
            }
        }
        Ok(())
    }

    fn opt_u32(&mut self, value: Option<u32>) -> Result {
        match value {
            None => self.u8(0),
            Some(value) => {
                self.u8(1)?;
                self.u32(value)
            }
        }
    }
}

/// A string column being gathered: its table of distinct values, in order
/// of first appearance, and each row's index into it.
#[derive(Default)]
struct Dictionary<'s> {
    table: Vec<Option<&'s str>>,
    positions: HashMap<Option<&'s str>, u32, BuildHasherDefault<WordHasher>>,
    previous: Option<(Option<&'s str>, u32)>,
    indices: Vec<u32>,
}

impl<'s> Dictionary<'s> {
    fn push(&mut self, value: Option<&'s str>) {
        // Rows of one residue or instance repeat their value, so a repeat of
        // the previous row's value skips the lookup.
        let index = match self.previous {
            Some((last, index)) if last == value => index,
            _ => {
                let table = &mut self.table;
                let index = *self.positions.entry(value).or_insert_with(|| {
                    table.push(value);
                    // A column has at most one value per row, and rows are
                    // counted in `u32`.
                    u32::try_from(table.len() - 1).unwrap_or(u32::MAX)
                });
                self.previous = Some((value, index));
                index
            }
        };
        self.indices.push(index);
    }
}

/// The dictionary columns of atom provenance.
const ATOM_TEXTS: usize = 15;

/// The dictionary columns of one atom, in their stored order.
fn atom_texts(atom: &AtomProvenance) -> [Option<&str>; ATOM_TEXTS] {
    [
        Some(&atom.type_symbol),
        atom.label_atom_name.as_deref(),
        Some(&atom.atom_name),
        atom.auth_atom_name.as_deref(),
        atom.label_component_id.as_deref(),
        Some(&atom.component_id),
        atom.auth_component_id.as_deref(),
        atom.label_asym_id.as_deref(),
        Some(&atom.asym_id),
        atom.auth_asym_id.as_deref(),
        atom.entity_id.as_deref(),
        Some(&atom.residue_key),
        atom.author_sequence_id.as_deref(),
        atom.insertion_code.as_deref(),
        atom.selected_alternate_location.as_deref(),
    ]
}

/// Bytes of packed values staged before they reach the encoder.
const PACKED_CHUNK: usize = 8 * 1024;

/// A fast hash of short strings, which the encoder's own data, not input,
/// feeds.
#[derive(Default)]
struct WordHasher(u64);

impl Hasher for WordHasher {
    fn write(&mut self, bytes: &[u8]) {
        let (words, rest) = bytes.as_chunks::<8>();
        for word in words {
            self.add(u64::from_le_bytes(*word));
        }
        let mut last = [0; 8];
        last[..rest.len()].copy_from_slice(rest);
        self.add(u64::from_le_bytes(last) ^ ((rest.len() as u64) << 59));
    }

    fn finish(&self) -> u64 {
        self.0
    }
}

impl WordHasher {
    fn add(&mut self, word: u64) {
        self.0 = (self.0.rotate_left(5) ^ word).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
}

fn entity_kind_code(kind: &EntityKind) -> u8 {
    match kind {
        EntityKind::Polymer => 0,
        EntityKind::Branched => 1,
        EntityKind::NonPolymer => 2,
        EntityKind::Water => 3,
        EntityKind::Other(_) => 4,
    }
}

fn evidence_code(evidence: &BondEvidence) -> u8 {
    match evidence {
        BondEvidence::ComponentTable => 0,
        BondEvidence::Template { .. } => 1,
        BondEvidence::PolymerLink => 2,
        BondEvidence::BranchLink => 3,
        BondEvidence::Connection { .. } => 4,
        BondEvidence::Distance => 5,
        BondEvidence::HydrogenDistance => 6,
    }
}

fn polymer_type_code(polymer_type: &PolymerType) -> u8 {
    match polymer_type {
        PolymerType::PolypeptideL => 0,
        PolymerType::PolypeptideD => 1,
        PolymerType::Polydeoxyribonucleotide => 2,
        PolymerType::Polyribonucleotide => 3,
        PolymerType::NucleotideHybrid => 4,
        PolymerType::PeptideNucleicAcid => 5,
        PolymerType::CyclicPseudoPeptide => 6,
        PolymerType::PolysaccharideD => 7,
        PolymerType::PolysaccharideL => 8,
        PolymerType::Other(_) => 9,
    }
}

const fn sequence_source_code(source: SequenceSource) -> u8 {
    match source {
        SequenceSource::Stated => 0,
        SequenceSource::Derived => 1,
    }
}
