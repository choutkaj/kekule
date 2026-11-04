//! Decoding of structure interpretation reports.

use std::io::Read;

use super::{Decoder, Result};
use crate::diagnostics::{Resource, SourceLocation};
use crate::io::{
    AltLocDecision, AltLocSelectionReason, AtomProvenance, BondEvidence, BondProvenance,
    CellOmissionReason, ConflictingSelector, ConnectionResolutionReason, CoordinationProvenance,
    EntityKind, EntityMonomer, InstanceProvenance, PolymerEntity, PolymerEntityInstance,
    PolymerType, SequencePosition, SequenceSource, SourceResidueId, SourceResiduePosition,
    StructureInterpretIssue, StructureInterpretationReport,
};
use crate::molecule::AtomId;
use crate::sequence::{CodeOrigin, ResidueCode};
use crate::topology::{ChainId, InstanceAtomId, MoleculeInstanceId, ResidueId};

impl<R: Read> Decoder<'_, R> {
    /// `count` unsigned values of the width the input states.
    pub(super) fn packed(&mut self, count: usize) -> Result<Vec<u64>> {
        let width = match self.u8()? {
            width @ (1 | 2 | 4 | 8) => usize::from(width),
            other => {
                return Err(self.invalid(format!(
                    "The integer width {other} is not 1, 2, 4, or 8 bytes."
                )))
            }
        };
        let mut values = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        self.chunks(count.saturating_mul(width), |bytes| {
            values.extend(bytes.chunks_exact(width).map(|raw| {
                let mut value = [0; 8];
                value[..width].copy_from_slice(raw);
                u64::from_le_bytes(value)
            }));
            Ok(())
        })?;
        Ok(values)
    }

    /// `count` packed values, each at most `max`.
    fn packed_codes(&mut self, count: usize, max: u64, what: &str) -> Result<Vec<u64>> {
        let values = self.packed(count)?;
        if let Some(value) = values.iter().find(|&&value| value > max) {
            return Err(self.invalid(format!("The {what} {value} is unknown.")));
        }
        Ok(values)
    }

    fn to_usize(&self, value: u64) -> Result<usize> {
        usize::try_from(value)
            .map_err(|_| self.invalid(format!("The value {value} is too large for this platform.")))
    }

    fn to_u32(&self, value: u64) -> Result<u32> {
        u32::try_from(value).map_err(|_| self.invalid(format!("The ID {value} is out of range.")))
    }

    fn usize(&mut self) -> Result<usize> {
        let value = self.u64()?;
        self.to_usize(value)
    }

    fn opt_usize(&mut self) -> Result<Option<usize>> {
        Ok(if self.tag("optional value")? {
            Some(self.usize()?)
        } else {
            None
        })
    }

    fn opt_u32(&mut self) -> Result<Option<u32>> {
        Ok(if self.tag("optional value")? {
            Some(self.u32()?)
        } else {
            None
        })
    }

    fn items(&mut self) -> Result<usize> {
        self.count(self.options.max_items(), Resource::Values)
    }

    fn usizes(&mut self) -> Result<Vec<usize>> {
        let count = self.items()?;
        self.packed(count)?
            .into_iter()
            .map(|value| self.to_usize(value))
            .collect()
    }

    fn str_list(&mut self) -> Result<Vec<String>> {
        let count = self.items()?;
        let mut values = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        for _ in 0..count {
            values.push(self.str()?);
        }
        Ok(values)
    }

    fn i32_array(&mut self, count: usize) -> Result<Vec<i32>> {
        let mut values = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        self.chunks(count.saturating_mul(4), |bytes| {
            values.extend(
                bytes
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .map(|raw| i32::from_le_bytes(*raw)),
            );
            Ok(())
        })?;
        Ok(values)
    }

    /// A string column of `rows` entries stored as a table of distinct
    /// values and one index per row; absent values are rejected when
    /// `required`.
    fn dictionary_column(&mut self, rows: usize, required: bool) -> Result<Dictionary> {
        let count = self.count(rows, Resource::Values)?;
        let table = self.strings(count, required)?;
        let indices = self.packed(rows)?;
        if let Some(index) = indices.iter().find(|&&index| index >= count as u64) {
            return Err(self.invalid(format!("The string index {index} is out of range.")));
        }
        Ok(Dictionary {
            table,
            // Every index is below the table length, a `u32` count.
            indices: indices.into_iter().map(|index| index as u32).collect(),
        })
    }

    fn dictionary(&mut self, rows: usize, required: bool) -> Result<Vec<Option<String>>> {
        let column = self.dictionary_column(rows, required)?;
        Ok((0..rows).map(|row| column.get(row)).collect())
    }

    fn required_dictionary(&mut self, rows: usize) -> Result<Vec<String>> {
        let column = self.dictionary_column(rows, true)?;
        Ok((0..rows).map(|row| column.required(row)).collect())
    }

    /// A column of atom references, `per` for each of the `lengths` items of
    /// each instance in `instances`. A reference stored as a local atom ID is
    /// to the instance of its item.
    fn atom_refs(
        &mut self,
        instances: &[MoleculeInstanceId],
        lengths: &[usize],
        per: usize,
    ) -> Result<Vec<InstanceAtomId>> {
        let count = lengths.iter().sum::<usize>() * per;
        let local = !self.tag("atom reference mode")?;
        // The values are read first, so that storage grows only as input
        // bytes arrive.
        let values = self.packed(count)?;
        if local {
            let instances = instances
                .iter()
                .zip(lengths)
                .flat_map(|(&instance, &length)| std::iter::repeat_n(instance, length * per));
            values
                .into_iter()
                .zip(instances)
                .map(|(value, instance)| {
                    Ok(InstanceAtomId::new(
                        instance,
                        AtomId::new(self.to_u32(value)?),
                    ))
                })
                .collect()
        } else {
            Ok(values.into_iter().map(InstanceAtomId::from_bits).collect())
        }
    }

    fn entity_kinds(&mut self, count: usize) -> Result<Vec<EntityKind>> {
        let codes = self.packed_codes(count, 4, "entity kind")?;
        let others = codes.iter().filter(|&&code| code == 4).count();
        let mut others = self.required_dictionary(others)?.into_iter();
        Ok(codes
            .into_iter()
            .map(|code| match code {
                0 => EntityKind::Polymer,
                1 => EntityKind::Branched,
                2 => EntityKind::NonPolymer,
                3 => EntityKind::Water,
                _ => EntityKind::Other(others.next().unwrap_or_default()),
            })
            .collect())
    }

    fn entity_kind(&mut self) -> Result<EntityKind> {
        Ok(match self.u8()? {
            0 => EntityKind::Polymer,
            1 => EntityKind::Branched,
            2 => EntityKind::NonPolymer,
            3 => EntityKind::Water,
            4 => EntityKind::Other(self.str()?),
            other => return Err(self.invalid(format!("The entity kind {other} is unknown."))),
        })
    }

    pub(super) fn report(&mut self) -> Result<StructureInterpretationReport> {
        let block_name = self.str()?;
        let selected_model = self.opt_str()?;
        let mut counts = [0usize; 7];
        for count in &mut counts {
            *count = self.usize()?;
        }
        let [entity_definitions, coordinate_models, macromolecules, small_molecules, solvent_molecules, applied_connections, connectivity_candidates] =
            counts;
        let ignored_coordinate_models = self.str_list()?;
        let instances = self.instance_provenance()?;
        let count = self.items()?;
        let mut coordination = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        for _ in 0..count {
            let atoms = [
                InstanceAtomId::from_bits(self.u64()?),
                InstanceAtomId::from_bits(self.u64()?),
            ];
            coordination.push(CoordinationProvenance {
                atoms,
                connection_id: self.opt_str()?,
                connection_type: self.str()?,
                source_line: self.opt_usize()?,
            });
        }
        let count = self.items()?;
        let mut alternate_locations = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        for _ in 0..count {
            alternate_locations.push(self.alternate_location()?);
        }
        let count = self.items()?;
        let mut issues = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        for _ in 0..count {
            issues.push(self.issue()?);
        }
        let count = self.items()?;
        let mut polymer_entities = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        for _ in 0..count {
            polymer_entities.push(self.polymer_entity()?);
        }
        Ok(StructureInterpretationReport {
            block_name,
            entity_definitions,
            coordinate_models,
            selected_model,
            ignored_coordinate_models,
            macromolecules,
            small_molecules,
            solvent_molecules,
            applied_connections,
            connectivity_candidates,
            instances,
            coordination,
            alternate_locations,
            issues,
            polymer_entities,
        })
    }

    /// Splits `values` into consecutive lists of the given lengths.
    fn split<T>(values: impl IntoIterator<Item = T>, lengths: &[usize]) -> Vec<Vec<T>> {
        let mut values = values.into_iter();
        lengths
            .iter()
            .map(|&length| values.by_ref().take(length).collect())
            .collect()
    }

    /// `count` lengths whose sum is at most `limit`.
    fn lengths(
        &mut self,
        count: usize,
        limit: usize,
        resource: Resource,
    ) -> Result<(Vec<usize>, usize)> {
        let lengths = self
            .packed(count)?
            .into_iter()
            .map(|value| self.to_usize(value))
            .collect::<Result<Vec<_>>>()?;
        let total = lengths
            .iter()
            .fold(0usize, |total, &length| total.saturating_add(length));
        if total > limit {
            return Err(self.limit(resource, limit, total));
        }
        Ok((lengths, total))
    }

    fn instance_provenance(&mut self) -> Result<Vec<InstanceProvenance>> {
        let max_items = self.options.max_items();
        let count = self.count(max_items, Resource::Components)?;
        let molecules = self
            .packed(count)?
            .into_iter()
            .map(|value| Ok(MoleculeInstanceId::new(self.to_u32(value)?)))
            .collect::<Result<Vec<_>>>()?;
        let coordinate_models = self.required_dictionary(count)?;
        let (lengths, total) = self.lengths(count, max_items, Resource::Values)?;
        let asym_ids = Self::split(self.required_dictionary(total)?, &lengths);
        let (lengths, total) = self.lengths(count, max_items, Resource::Values)?;
        let entity_ids = Self::split(self.required_dictionary(total)?, &lengths);
        let (lengths, total) = self.lengths(count, max_items, Resource::Values)?;
        let entity_kinds = Self::split(self.entity_kinds(total)?, &lengths);
        let (atom_lengths, atom_total) =
            self.lengths(count, self.options.max_atoms(), Resource::Atoms)?;
        let (bond_lengths, bond_total) =
            self.lengths(count, self.options.max_bonds(), Resource::Bonds)?;

        let refs = self.atom_refs(&molecules, &atom_lengths, 1)?;
        let source_lines = self.packed(atom_total)?;
        let source_orders = self.packed(atom_total)?;
        let atom_site_ids = self.strings(atom_total, false)?;
        // The dictionary columns in the order of the encoder's `atom_texts`.
        let mut texts = Vec::with_capacity(ATOM_TEXTS_REQUIRED.len());
        for required in ATOM_TEXTS_REQUIRED {
            texts.push(self.dictionary_column(atom_total, required)?);
        }
        let mut entity_kinds_of_atoms = self.entity_kinds(atom_total)?.into_iter();
        let has_label_sequence = self.bits(atom_total)?;
        let present = has_label_sequence
            .iter()
            .filter(|&&present| present)
            .count();
        let mut label_sequences = self.i32_array(present)?.into_iter();
        let has_occurrence = self.bits(atom_total)?;
        let present = has_occurrence.iter().filter(|&&present| present).count();
        let mut occurrences = self.packed(present)?.into_iter();

        let mut atom_site_ids = atom_site_ids.into_iter();
        let mut atoms = Vec::with_capacity(count);
        let mut row = 0;
        for &length in &atom_lengths {
            let mut instance_atoms = Vec::with_capacity(length);
            for _ in 0..length {
                let text = |column: usize| texts[column].get(row);
                let required = |column: usize| texts[column].required(row);
                let occurrence = match has_occurrence[row] {
                    true => Some(self.to_usize(occurrences.next().unwrap_or_default())?),
                    false => None,
                };
                instance_atoms.push(AtomProvenance {
                    atom: refs[row],
                    type_symbol: required(0),
                    source_line: self.to_usize(source_lines[row])?,
                    source_order: self.to_usize(source_orders[row])?,
                    atom_site_id: atom_site_ids.next().flatten(),
                    label_atom_name: text(1),
                    atom_name: required(2),
                    auth_atom_name: text(3),
                    label_component_id: text(4),
                    component_id: required(5),
                    auth_component_id: text(6),
                    label_asym_id: text(7),
                    asym_id: required(8),
                    auth_asym_id: text(9),
                    entity_id: text(10),
                    entity_kind: entity_kinds_of_atoms.next().unwrap_or(EntityKind::Water),
                    residue_key: required(11),
                    label_sequence_id: has_label_sequence[row]
                        .then(|| label_sequences.next().unwrap_or_default()),
                    author_sequence_id: text(12),
                    insertion_code: text(13),
                    occurrence,
                    selected_alternate_location: text(14),
                });
                row += 1;
            }
            atoms.push(instance_atoms);
        }

        let refs = self.atom_refs(&molecules, &bond_lengths, 2)?;
        let codes = self.packed_codes(bond_total, 6, "bond evidence")?;
        let templates = codes.iter().filter(|&&code| code == 1).count();
        let mut template_components = self.required_dictionary(templates)?.into_iter();
        let mut template_variants = self.dictionary(templates, false)?.into_iter();
        let connections = codes.iter().filter(|&&code| code == 4).count();
        let mut connection_ids = self.dictionary(connections, false)?.into_iter();
        let mut connection_types = self.required_dictionary(connections)?.into_iter();
        let bonds = refs
            .as_chunks::<2>()
            .0
            .iter()
            .zip(codes)
            .map(|(&atoms, code)| BondProvenance {
                atoms,
                evidence: match code {
                    0 => BondEvidence::ComponentTable,
                    1 => BondEvidence::Template {
                        component_id: template_components.next().unwrap_or_default(),
                        variant: template_variants.next().flatten(),
                    },
                    2 => BondEvidence::PolymerLink,
                    3 => BondEvidence::BranchLink,
                    4 => BondEvidence::Connection {
                        connection_id: connection_ids.next().flatten(),
                        connection_type: connection_types.next().unwrap_or_default(),
                    },
                    5 => BondEvidence::Distance,
                    _ => BondEvidence::HydrogenDistance,
                },
            })
            .collect::<Vec<_>>();
        let bonds = Self::split(bonds, &bond_lengths);
        Ok(molecules
            .into_iter()
            .zip(coordinate_models)
            .zip(asym_ids.into_iter().zip(entity_ids).zip(entity_kinds))
            .zip(atoms.into_iter().zip(bonds))
            .map(
                |(
                    ((molecule, coordinate_model_id), ((asym_ids, entity_ids), entity_kinds)),
                    (atoms, bonds),
                )| {
                    InstanceProvenance {
                        molecule,
                        coordinate_model_id,
                        asym_ids,
                        entity_ids,
                        entity_kinds,
                        atoms,
                        bonds,
                    }
                },
            )
            .collect())
    }

    fn residue(&mut self) -> Result<SourceResidueId> {
        let model_id = self.str()?;
        let asym_id = self.str()?;
        let position = match self.u8()? {
            0 => SourceResiduePosition::Label(self.i32()?),
            1 => SourceResiduePosition::Author(self.str()?),
            2 => SourceResiduePosition::RepeatedAuthor {
                sequence: self.str()?,
                repeat: self.usize()?,
            },
            3 => SourceResiduePosition::Occurrence {
                component_id: self.str()?,
                index: self.usize()?,
            },
            other => {
                return Err(self.invalid(format!("The residue position kind {other} is unknown.")))
            }
        };
        Ok(SourceResidueId {
            model_id,
            asym_id,
            position,
            insertion_code: self.opt_str()?,
        })
    }

    fn alternate_location(&mut self) -> Result<AltLocDecision> {
        let residue = self.residue()?;
        let available_labels = self.str_list()?;
        let selected_labels = self.str_list()?;
        let reason = match self.u8()? {
            0 => AltLocSelectionReason::HighestOccupancy,
            1 => AltLocSelectionReason::ExactSelection,
            2 => AltLocSelectionReason::PreferredLabel,
            3 => AltLocSelectionReason::PreferredLabelUnavailable,
            4 => AltLocSelectionReason::SourceConformation { id: self.str()? },
            other => {
                return Err(self.invalid(format!(
                    "The alternate-location selection reason {other} is unknown."
                )))
            }
        };
        Ok(AltLocDecision {
            residue,
            available_labels,
            selected_labels,
            reason,
            selected_source_lines: self.usizes()?,
            omitted_source_lines: self.usizes()?,
        })
    }

    fn charged_atoms(&mut self) -> Result<Vec<(String, i8)>> {
        let count = self.items()?;
        let mut atoms = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        for _ in 0..count {
            let name = self.str()?;
            atoms.push((name, i8::from_le_bytes(self.array()?)));
        }
        Ok(atoms)
    }

    fn connection_reason(&mut self) -> Result<ConnectionResolutionReason> {
        Ok(match self.u8()? {
            0 => ConnectionResolutionReason::UnsupportedSymmetry {
                symmetry: self.str()?,
            },
            1 => ConnectionResolutionReason::MissingSelector,
            2 => ConnectionResolutionReason::ConflictingSelectorValues {
                selector: match self.u8()? {
                    0 => ConflictingSelector::LabelAlternateLocationAliases,
                    1 => ConflictingSelector::LabelAndAuthorAlternateLocations,
                    other => {
                        return Err(
                            self.invalid(format!("The connection selector {other} is unknown."))
                        )
                    }
                },
            },
            3 => ConnectionResolutionReason::ConflictingLabelAndAuthorSelectors,
            4 => ConnectionResolutionReason::AlternateLocationOmitted {
                alternate_location: self.str()?,
            },
            5 => ConnectionResolutionReason::NoMatchingAtom,
            6 => ConnectionResolutionReason::MultipleMatchingAtoms,
            other => {
                return Err(self.invalid(format!(
                    "The connection resolution reason {other} is unknown."
                )))
            }
        })
    }

    fn issue(&mut self) -> Result<StructureInterpretIssue> {
        use StructureInterpretIssue as Issue;
        Ok(match self.u8()? {
            0 => Issue::EntityTypeInferred {
                asym_id: self.str()?,
                kind: self.entity_kind()?,
            },
            1 => Issue::ElementsInferredFromAtomNames {
                source_lines: self.usizes()?,
            },
            2 => Issue::ConnectionsIgnored {
                connection_type: self.str()?,
                count: self.usize()?,
            },
            3 => Issue::ConnectionUnresolved {
                connection_id: self.opt_str()?,
                connection_type: self.str()?,
                partner: self.u8()?,
                source_line: self.opt_usize()?,
                reason: self.connection_reason()?,
            },
            4 => Issue::ConnectionAmbiguous {
                connection_id: self.opt_str()?,
                connection_type: self.str()?,
                partner: self.u8()?,
                source_line: self.opt_usize()?,
                candidates: self.usize()?,
                reason: self.connection_reason()?,
            },
            5 => Issue::CoordinateModelIgnored {
                model_id: self.str()?,
                atom_site_rows: self.usize()?,
            },
            6 => Issue::AlternateLocationOmitted {
                atom_name: self.str()?,
                alt_id: self.opt_str()?,
                residue: self.residue()?,
                source_line: self.usize()?,
            },
            7 => Issue::AlternateLocationIncomplete {
                residue: self.residue()?,
                alternative: self.str()?,
                missing_atoms: self.str_list()?,
            },
            8 => Issue::AlternateLocationOccupancyMissing {
                residue: self.residue()?,
                source_lines: self.usizes()?,
            },
            9 => Issue::AlternateLocationCombinationRejected {
                residue: self.residue()?,
                alternative: self.str()?,
                reason: self.str()?,
            },
            10 => {
                let count = self.items()?;
                let mut residues = Vec::with_capacity(count.min(super::CHUNK_BYTES));
                for _ in 0..count {
                    residues.push(self.residue()?);
                }
                Issue::AlternateLocationOccupancyTie {
                    residues,
                    alternatives: self.str_list()?,
                }
            }
            11 => Issue::ConnectivityCandidatesInferred {
                atom_count: self.usize()?,
                candidate_count: self.usize()?,
            },
            12 => Issue::ResidueBondsInferred {
                residue: self.residue()?,
                component_id: self.str()?,
                bonds: self.usize()?,
            },
            13 => Issue::HydrogenAttachedByDistance {
                residue: self.residue()?,
                atom_name: self.str()?,
                heavy_atom_name: self.str()?,
            },
            14 => Issue::ChainBreak {
                residue: self.residue()?,
                next: self.residue()?,
            },
            15 => Issue::ChargesFromBonds {
                residue: self.residue()?,
                component_id: self.str()?,
                charged_atoms: self.charged_atoms()?,
            },
            16 => Issue::ProtonationFromHydrogens {
                residue: self.residue()?,
                component_id: self.str()?,
                variant: self.opt_str()?,
                charged_atoms: self.charged_atoms()?,
            },
            17 => Issue::ResidueFragmented {
                residue: self.residue()?,
                fragments: self.usize()?,
            },
            18 => {
                let reason = match self.u8()? {
                    0 => CellOmissionReason::Incomplete,
                    1 => CellOmissionReason::SpaceGroup { name: self.str()? },
                    2 => CellOmissionReason::Placeholder,
                    other => {
                        return Err(
                            self.invalid(format!("The cell omission reason {other} is unknown."))
                        )
                    }
                };
                let line = self.usize()?;
                let location = match self.opt_usize()? {
                    Some(column) => SourceLocation::at_column(line, column),
                    None => SourceLocation::at_line(line),
                };
                Issue::CellOmitted { reason, location }
            }
            other => {
                return Err(
                    self.invalid(format!("The interpretation issue kind {other} is unknown."))
                )
            }
        })
    }

    fn sequence_source(&mut self) -> Result<SequenceSource> {
        match self.u8()? {
            0 => Ok(SequenceSource::Stated),
            1 => Ok(SequenceSource::Derived),
            other => Err(self.invalid(format!("The sequence source {other} is unknown."))),
        }
    }

    fn polymer_entity(&mut self) -> Result<PolymerEntity> {
        let id = self.str()?;
        let polymer_type = if self.tag("polymer type")? {
            Some(match self.u8()? {
                0 => PolymerType::PolypeptideL,
                1 => PolymerType::PolypeptideD,
                2 => PolymerType::Polydeoxyribonucleotide,
                3 => PolymerType::Polyribonucleotide,
                4 => PolymerType::NucleotideHybrid,
                5 => PolymerType::PeptideNucleicAcid,
                6 => PolymerType::CyclicPseudoPeptide,
                7 => PolymerType::PolysaccharideD,
                8 => PolymerType::PolysaccharideL,
                9 => PolymerType::Other(self.str()?),
                other => return Err(self.invalid(format!("The polymer type {other} is unknown."))),
            })
        } else {
            None
        };
        let monomers_source = self.sequence_source()?;
        let sequence = self.str()?;
        let sequence_source = self.sequence_source()?;
        let count = self.items()?;
        let mut monomers = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        for _ in 0..count {
            let number = self.i32()?;
            let component_id = self.str()?;
            let alternatives = self.str_list()?;
            let code = self.str()?;
            let origin = match self.u8()? {
                0 => CodeOrigin::Component,
                1 => CodeOrigin::Parent,
                2 => CodeOrigin::Sequence,
                3 => CodeOrigin::Unknown,
                other => return Err(self.invalid(format!("The code origin {other} is unknown."))),
            };
            let parents = self.str_list()?;
            monomers.push(EntityMonomer {
                number,
                component_id,
                alternatives,
                code: ResidueCode {
                    code: code.into(),
                    origin,
                    parents: parents.into_iter().map(String::into_boxed_str).collect(),
                },
            });
        }
        let count = self.items()?;
        let mut instances = Vec::with_capacity(count.min(super::CHUNK_BYTES));
        for _ in 0..count {
            let asym_id = self.str()?;
            let author_chain_id = self.opt_str()?;
            let chain = self.opt_u32()?.map(ChainId::new);
            let count = self.items()?;
            let mut positions = Vec::with_capacity(count.min(super::CHUNK_BYTES));
            for _ in 0..count {
                positions.push(SequencePosition {
                    residue: self.opt_u32()?.map(ResidueId::new),
                    author_seq_id: self.opt_str()?,
                    insertion_code: self.opt_str()?,
                });
            }
            instances.push(PolymerEntityInstance {
                asym_id,
                author_chain_id,
                chain,
                positions,
            });
        }
        Ok(PolymerEntity {
            id,
            polymer_type,
            monomers,
            monomers_source,
            sequence,
            sequence_source,
            instances,
        })
    }
}

/// Whether each dictionary column of atom provenance is required, in the
/// order of the encoder's `atom_texts`.
const ATOM_TEXTS_REQUIRED: [bool; 15] = [
    true, false, true, false, false, true, false, false, true, false, false, true, false, false,
    false,
];

/// A decoded string column: its distinct values and each row's index into
/// them, every index checked.
struct Dictionary {
    table: Vec<Option<String>>,
    indices: Vec<u32>,
}

impl Dictionary {
    fn get(&self, row: usize) -> Option<String> {
        self.table[self.indices[row] as usize].clone()
    }

    /// The value of a row of a required column, whose table has no absent
    /// values.
    fn required(&self, row: usize) -> String {
        self.get(row).unwrap_or_default()
    }
}
