//! Native round trips of interpretation reports that hold every variant of
//! every report field, and rejection of corrupt report sections.

use crate::diagnostics::{DiagnosticKind, Format, SourceLocation};
use crate::io::{
    AltLocDecision, AltLocSelectionReason, Assembly, AssemblyGenerator, AssemblyOmissionReason,
    AssemblyOperation, AtomProvenance, BondEvidence, BondProvenance, CellOmissionReason,
    ConflictingSelector, ConnectionResolutionReason, CoordinationProvenance, EntityKind,
    EntityMonomer, InstanceProvenance, PolymerEntity, PolymerEntityInstance, PolymerType,
    SequencePosition, SequenceSource, SourceResidueId, SourceResiduePosition,
    StructureInterpretIssue, StructureInterpretationReport, StructureSetInterpretation,
};
use crate::math::{Matrix3, Vector3};
use crate::molecule::AtomId;
use crate::native::{self, NativeReadOptions, NativeSource, NativeWriteOptions};
use crate::sequence::{CodeOrigin, ResidueCode};
use crate::structure::{Structure, StructureSet};
use crate::topology::{ChainId, InstanceAtomId, MoleculeInstanceId, ResidueId};

const MODELS: &str = "data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
O O HOH A 1 0.0 0.0 0.0 1
O O HOH B 2 3.0 0.0 0.0 1
O O HOH A 1 0.1 0.0 0.0 2
O O HOH B 2 3.1 0.0 0.0 2
";

fn atom(instance: u32, atom: u32) -> InstanceAtomId {
    InstanceAtomId::new(MoleculeInstanceId::new(instance), AtomId::new(atom))
}

fn residue(position: SourceResiduePosition, insertion_code: Option<&str>) -> SourceResidueId {
    SourceResidueId::from_parts("1", "A", position, insertion_code.map(str::to_owned))
}

fn label(sequence: i32) -> SourceResidueId {
    residue(SourceResiduePosition::Label(sequence), None)
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_owned()).collect()
}

fn atom_provenance(id: InstanceAtomId, order: usize, full: bool) -> AtomProvenance {
    let text = |value: &str| full.then(|| value.to_owned());
    AtomProvenance {
        atom: id,
        type_symbol: "O".into(),
        source_line: 1_000_000 + order,
        source_order: order,
        atom_site_id: text(&format!("{}", order + 1)),
        label_atom_name: text("O"),
        atom_name: "O".into(),
        auth_atom_name: text("OW"),
        label_component_id: text("HOH"),
        component_id: "HOH".into(),
        auth_component_id: text("WAT"),
        label_asym_id: text("A"),
        asym_id: "A".into(),
        auth_asym_id: text("Z"),
        entity_id: text("1"),
        entity_kind: if full {
            EntityKind::Water
        } else {
            EntityKind::Other("macrolide".into())
        },
        residue_key: format!("{order}"),
        label_sequence_id: full.then_some(-7),
        author_sequence_id: text("12"),
        insertion_code: text("B"),
        occurrence: (!full).then_some(u32::MAX as usize + 3),
        selected_alternate_location: text("A"),
    }
}

/// A report with every issue, bond evidence, selection reason, residue
/// position, polymer type, and code origin. Its atoms all refer to their own
/// instance unless `foreign`.
pub(super) fn report(foreign: bool) -> StructureInterpretationReport {
    use StructureInterpretIssue as Issue;
    let other = if foreign { atom(9, 4) } else { atom(1, 1) };
    let instances = vec![
        InstanceProvenance {
            molecule: MoleculeInstanceId::new(0),
            coordinate_model_id: "MODELX".into(),
            asym_ids: strings(&["A"]),
            entity_ids: strings(&["1", "2"]),
            entity_kinds: vec![EntityKind::Water, EntityKind::Other("odd".into())],
            atoms: vec![atom_provenance(atom(0, 0), 0, true)],
            bonds: Vec::new(),
        },
        InstanceProvenance {
            molecule: MoleculeInstanceId::new(1),
            coordinate_model_id: "MODELX".into(),
            asym_ids: strings(&["B", "C"]),
            entity_ids: Vec::new(),
            entity_kinds: vec![EntityKind::Polymer, EntityKind::Branched],
            atoms: vec![
                atom_provenance(atom(1, 0), 1, false),
                atom_provenance(other, 2, true),
            ],
            bonds: [
                BondEvidence::ComponentTable,
                BondEvidence::Template {
                    component_id: "HIS".into(),
                    variant: Some("HIS_LSN".into()),
                },
                BondEvidence::Template {
                    component_id: "ALA".into(),
                    variant: None,
                },
                BondEvidence::PolymerLink,
                BondEvidence::BranchLink,
                BondEvidence::Connection {
                    connection_id: Some("covale1".into()),
                    connection_type: "covale".into(),
                },
                BondEvidence::Connection {
                    connection_id: None,
                    connection_type: "LINK".into(),
                },
                BondEvidence::Distance,
                BondEvidence::HydrogenDistance,
            ]
            .into_iter()
            .map(|evidence| BondProvenance {
                atoms: [atom(1, 0), other],
                evidence,
            })
            .collect(),
        },
    ];
    let residues = [
        label(-3),
        residue(SourceResiduePosition::Author("12".into()), Some("A")),
        residue(
            SourceResiduePosition::RepeatedAuthor {
                sequence: "12".into(),
                repeat: 2,
            },
            None,
        ),
        residue(
            SourceResiduePosition::Occurrence {
                component_id: "HOH".into(),
                index: 4,
            },
            None,
        ),
    ];
    let reasons = [
        ConnectionResolutionReason::UnsupportedSymmetry {
            symmetry: "2_555".into(),
        },
        ConnectionResolutionReason::MissingSelector,
        ConnectionResolutionReason::ConflictingSelectorValues {
            selector: ConflictingSelector::LabelAlternateLocationAliases,
        },
        ConnectionResolutionReason::ConflictingSelectorValues {
            selector: ConflictingSelector::LabelAndAuthorAlternateLocations,
        },
        ConnectionResolutionReason::ConflictingLabelAndAuthorSelectors,
        ConnectionResolutionReason::AlternateLocationOmitted {
            alternate_location: "B".into(),
        },
        ConnectionResolutionReason::NoMatchingAtom,
        ConnectionResolutionReason::MultipleMatchingAtoms,
    ];
    let mut issues = vec![
        Issue::EntityTypeInferred {
            asym_id: "B".into(),
            kind: EntityKind::NonPolymer,
        },
        Issue::EntityTypeInferred {
            asym_id: "B".into(),
            kind: EntityKind::Other("odd".into()),
        },
        Issue::ElementsInferredFromAtomNames {
            source_lines: vec![12, 13, 70_000],
        },
        Issue::ConnectionsIgnored {
            connection_type: "hydrog".into(),
            count: 1_204,
        },
        Issue::CoordinateModelIgnored {
            model_id: "QQQQ".into(),
            atom_site_rows: 1_632,
        },
        Issue::AlternateLocationOmitted {
            atom_name: "CG".into(),
            alt_id: Some("B".into()),
            residue: residues[1].clone(),
            source_line: 300,
        },
        Issue::AlternateLocationOmitted {
            atom_name: "CG".into(),
            alt_id: None,
            residue: residues[2].clone(),
            source_line: 301,
        },
        Issue::AlternateLocationIncomplete {
            residue: residues[3].clone(),
            alternative: "B".into(),
            missing_atoms: strings(&["CD", "NE"]),
        },
        Issue::AlternateLocationOccupancyMissing {
            residue: label(12),
            source_lines: vec![300],
        },
        Issue::AlternateLocationCombinationRejected {
            residue: label(12),
            alternative: "1".into(),
            reason: "The combination has no alternate location of this residue.".into(),
        },
        Issue::AlternateLocationOccupancyTie {
            residues: residues.to_vec(),
            alternatives: strings(&["A", "B"]),
        },
        Issue::ConnectivityCandidatesInferred {
            atom_count: 2_410,
            candidate_count: 3,
        },
        Issue::ResidueBondsInferred {
            residue: label(401),
            component_id: "LIG".into(),
            bonds: 24,
        },
        Issue::HydrogenAttachedByDistance {
            residue: label(401),
            atom_name: "H12".into(),
            heavy_atom_name: "C5".into(),
        },
        Issue::ChainBreak {
            residue: label(12),
            next: label(13),
        },
        Issue::ChargesFromBonds {
            residue: label(5),
            component_id: "XYZ".into(),
            charged_atoms: vec![("N1".into(), 1), ("B2".into(), -1)],
        },
        Issue::ProtonationFromHydrogens {
            residue: label(57),
            component_id: "HIS".into(),
            variant: Some("HIS_LSN".into()),
            charged_atoms: Vec::new(),
        },
        Issue::ProtonationFromHydrogens {
            residue: label(15),
            component_id: "LYS".into(),
            variant: None,
            charged_atoms: vec![("NZ".into(), 0)],
        },
        Issue::ResidueFragmented {
            residue: label(7),
            fragments: 2,
        },
        Issue::CellOmitted {
            reason: CellOmissionReason::Incomplete,
            location: SourceLocation::at_line(8),
        },
        Issue::CellOmitted {
            reason: CellOmissionReason::SpaceGroup {
                name: "P 21 21 21".into(),
            },
            location: SourceLocation::at_column(8, 3),
        },
        Issue::CellOmitted {
            reason: CellOmissionReason::Placeholder,
            location: SourceLocation::at_line(9),
        },
        Issue::AssemblyOmitted {
            assembly: None,
            reason: AssemblyOmissionReason::Unidentified,
            line: Some(30),
        },
        Issue::AssemblyOmitted {
            assembly: Some("1".into()),
            reason: AssemblyOmissionReason::NoGenerators,
            line: None,
        },
        Issue::AssemblyOmitted {
            assembly: Some("2".into()),
            reason: AssemblyOmissionReason::Expression {
                expression: "(1-".into(),
            },
            line: Some(31),
        },
        Issue::AssemblyOmitted {
            assembly: Some("3".into()),
            reason: AssemblyOmissionReason::UnknownOperator {
                operator: "61".into(),
            },
            line: Some(32),
        },
        Issue::AssemblyOmitted {
            assembly: Some("4".into()),
            reason: AssemblyOmissionReason::UnreadableOperator {
                operator: "P".into(),
            },
            line: Some(33),
        },
        Issue::AssemblyOmitted {
            assembly: Some("5".into()),
            reason: AssemblyOmissionReason::NotRotation {
                operator: "X0".into(),
            },
            line: Some(34),
        },
        Issue::AssemblyOmitted {
            assembly: Some("6".into()),
            reason: AssemblyOmissionReason::TooManyOperations,
            line: None,
        },
    ];
    for (index, reason) in reasons.into_iter().enumerate() {
        issues.push(if index % 2 == 0 {
            Issue::ConnectionUnresolved {
                connection_id: Some(format!("metalc{index}")),
                connection_type: "metalc".into(),
                partner: 2,
                source_line: Some(40 + index),
                reason,
            }
        } else {
            Issue::ConnectionAmbiguous {
                connection_id: None,
                connection_type: "LINK".into(),
                partner: 1,
                source_line: None,
                candidates: 2,
                reason,
            }
        });
    }
    let alternate_locations = [
        AltLocSelectionReason::HighestOccupancy,
        AltLocSelectionReason::ExactSelection,
        AltLocSelectionReason::PreferredLabel,
        AltLocSelectionReason::PreferredLabelUnavailable,
        AltLocSelectionReason::SourceConformation { id: "2".into() },
    ]
    .into_iter()
    .zip(residues.iter().cycle())
    .map(|(reason, residue)| AltLocDecision {
        residue: residue.clone(),
        available_labels: strings(&["A", "B"]),
        selected_labels: strings(&["A"]),
        reason,
        selected_source_lines: vec![10, 11],
        omitted_source_lines: vec![12, 1 << 40],
    })
    .collect();
    let code = |code: &str, origin, parents: &[&str]| ResidueCode {
        code: code.into(),
        origin,
        parents: parents.iter().map(|parent| (*parent).into()).collect(),
    };
    let monomers = vec![
        EntityMonomer {
            number: 1,
            component_id: "MSE".into(),
            alternatives: Vec::new(),
            code: code("M", CodeOrigin::Parent, &["MET"]),
        },
        EntityMonomer {
            number: 2,
            component_id: "GLY".into(),
            alternatives: strings(&["ALA"]),
            code: code("G", CodeOrigin::Component, &[]),
        },
        EntityMonomer {
            number: 3,
            component_id: "CRO".into(),
            alternatives: Vec::new(),
            code: code("TYG", CodeOrigin::Sequence, &["THR", "TYR", "GLY"]),
        },
        EntityMonomer {
            number: 4,
            component_id: "UNK".into(),
            alternatives: Vec::new(),
            code: code("X", CodeOrigin::Unknown, &[]),
        },
    ];
    let positions = vec![
        SequencePosition {
            residue: Some(ResidueId::new(0)),
            author_seq_id: Some("1".into()),
            insertion_code: None,
        },
        SequencePosition::default(),
        SequencePosition {
            residue: None,
            author_seq_id: Some("3".into()),
            insertion_code: Some("A".into()),
        },
        SequencePosition::default(),
    ];
    let polymer_types = [
        Some(PolymerType::PolypeptideL),
        Some(PolymerType::PolypeptideD),
        Some(PolymerType::Polydeoxyribonucleotide),
        Some(PolymerType::Polyribonucleotide),
        Some(PolymerType::NucleotideHybrid),
        Some(PolymerType::PeptideNucleicAcid),
        Some(PolymerType::CyclicPseudoPeptide),
        Some(PolymerType::PolysaccharideD),
        Some(PolymerType::PolysaccharideL),
        Some(PolymerType::Other("other".into())),
        None,
    ];
    let polymer_entities = polymer_types
        .into_iter()
        .enumerate()
        .map(|(index, polymer_type)| PolymerEntity {
            id: format!("{}", index + 1),
            polymer_type,
            monomers: monomers.clone(),
            monomers_source: SequenceSource::Stated,
            sequence: "MGTYGX".into(),
            sequence_source: if index % 2 == 0 {
                SequenceSource::Derived
            } else {
                SequenceSource::Stated
            },
            instances: vec![
                PolymerEntityInstance {
                    asym_id: "A".into(),
                    author_chain_id: Some("Z".into()),
                    chain: Some(ChainId::new(0)),
                    positions: positions.clone(),
                },
                PolymerEntityInstance {
                    asym_id: "B".into(),
                    author_chain_id: None,
                    chain: None,
                    positions: vec![SequencePosition::default(); monomers.len()],
                },
            ],
        })
        .collect();
    StructureInterpretationReport {
        block_name: "demo".into(),
        entity_definitions: 2,
        coordinate_models: 3,
        selected_model: Some("1".into()),
        ignored_coordinate_models: strings(&["2", "3"]),
        macromolecules: 4,
        small_molecules: 5,
        solvent_molecules: 6,
        applied_connections: 7,
        connectivity_candidates: 8,
        instances,
        coordination: vec![
            CoordinationProvenance {
                atoms: [atom(0, 0), atom(1, 0)],
                connection_id: Some("metalc1".into()),
                connection_type: "metalc".into(),
                source_line: Some(77),
            },
            CoordinationProvenance {
                atoms: [atom(1, 0), atom(0, 0)],
                connection_id: None,
                connection_type: "LINK".into(),
                source_line: None,
            },
        ],
        alternate_locations,
        issues,
        polymer_entities,
        assemblies: vec![
            Assembly {
                id: "1".into(),
                details: Some("complete icosahedral assembly".into()),
                oligomeric_details: Some("60-meric".into()),
                oligomeric_count: Some(60),
                generators: vec![
                    AssemblyGenerator {
                        chains: strings(&["A", "B"]),
                        operations: vec![
                            AssemblyOperation {
                                operators: strings(&["1"]),
                                rotation: Matrix3::identity(),
                                translation: Vector3::zero(),
                            },
                            AssemblyOperation {
                                operators: strings(&["2", "61"]),
                                rotation: Matrix3::from_columns(
                                    Vector3::new(0.0, 1.0, 0.0),
                                    Vector3::new(-1.0, 0.0, 0.0),
                                    Vector3::new(0.0, 0.0, 1.0),
                                ),
                                translation: Vector3::new(-0.1, 1.0 / 3.0, f64::MIN_POSITIVE),
                            },
                        ],
                    },
                    AssemblyGenerator {
                        chains: strings(&["C"]),
                        operations: Vec::new(),
                    },
                ],
            },
            Assembly {
                id: "2".into(),
                details: None,
                oligomeric_details: None,
                oligomeric_count: None,
                generators: Vec::new(),
            },
        ],
    }
}

/// The first coordinate model of [`MODELS`] and the set of both models,
/// whose two geometries carry stored equal weights.
pub(super) fn parts() -> (Structure, StructureSet) {
    let document = crate::mmcif::parse_str(MODELS, Default::default()).unwrap();
    let structure = document
        .interpret(
            crate::macromolecular::StructureInterpretOptions::default()
                .with_model_selection(crate::macromolecular::ModelSelection::First),
        )
        .unwrap()
        .into_structure();
    let mut set = document
        .interpret_structure_set(Default::default())
        .unwrap()
        .into_structure_set();
    set.distribution_mut()
        .unwrap()
        .set_weights(Some(vec![1.0, 1.0]))
        .unwrap();
    (structure, set)
}

pub(super) fn read(
    bytes: &[u8],
) -> Result<crate::native::NativeObject, crate::diagnostics::FormatError> {
    native::read(bytes, NativeReadOptions::default())
}

#[test]
fn every_report_field_round_trips() {
    let (structure, set) = parts();
    assert_eq!(super::debug_fingerprint(&set), 0xbff9_29d2_6f24_4c64);
    let report = report(false);
    let bytes = native::write(
        NativeSource::interpreted_structure(&structure, &report),
        NativeWriteOptions::default(),
    )
    .unwrap();
    let decoded = read(&bytes).unwrap().into_interpretation().unwrap();
    assert_eq!(decoded.report(), &report);
    assert_eq!(decoded.structure().positions(), structure.positions());
    let again = native::write(&decoded, NativeWriteOptions::default()).unwrap();
    assert!(
        again == bytes,
        "re-encoding a decoded report changed its bytes"
    );

    let reports = vec![report, self::report(true)];
    let bytes = native::write(
        NativeSource::interpreted_structure_set(&set, &reports),
        NativeWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(super::debug_fingerprint(&bytes), 0xde2b_74f7_fb64_0b01);
    let decoded = read(&bytes)
        .unwrap()
        .into_structure_set_interpretation()
        .unwrap();
    assert_eq!(super::debug_fingerprint(&decoded), 0x63e9_9b07_512b_83a0);
    assert_eq!(decoded.reports(), reports.as_slice());
    assert_eq!(decoded.structure_set().len(), 2);
    let again = native::write(&decoded, NativeWriteOptions::default()).unwrap();
    assert!(
        again == bytes,
        "re-encoding decoded reports changed their bytes"
    );
}

/// CRC-32 (IEEE), the native trailer's checksum.
fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = u32::MAX;
    for &byte in bytes {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = if crc & 1 == 1 {
                (crc >> 1) ^ 0xEDB8_8320
            } else {
                crc >> 1
            };
        }
    }
    !crc
}

/// `bytes` with its checksum recomputed, so that only the decoder's
/// validation can reject a corruption.
pub(super) fn resealed(mut bytes: Vec<u8>) -> Vec<u8> {
    let end = bytes.len() - 4;
    let crc = crc32(&bytes[..end]);
    bytes[end..].copy_from_slice(&crc.to_le_bytes());
    bytes
}

pub(super) fn find(bytes: &[u8], needle: &[u8]) -> usize {
    bytes
        .windows(needle.len())
        .position(|window| window == needle)
        .expect("the needle is in the encoding")
}

#[test]
fn decoded_assembly_operations_are_rigid() {
    let (structure, _) = parts();
    let mut report = report(false);
    // A mirror, as no reader would have read it.
    report.assemblies[0].generators[0].operations[1].rotation = Matrix3::from_columns(
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    );
    let bytes = native::write(
        NativeSource::interpreted_structure(&structure, &report),
        NativeWriteOptions::default(),
    )
    .unwrap();
    let error = read(&bytes).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::InvalidValue, "{error}");
    assert!(
        error
            .to_string()
            .contains("Assembly operation \"2-61\" is not a rotation and a translation."),
        "{error}"
    );
}

#[test]
fn corrupt_report_sections_fail_cleanly() {
    let (structure, _) = parts();
    let report = report(true);
    let plain = native::write(&structure, NativeWriteOptions::default()).unwrap();
    let bytes = native::write(
        NativeSource::interpreted_structure(&structure, &report),
        NativeWriteOptions::default(),
    )
    .unwrap();
    // The report follows the structure, before the trailer. Only the kind
    // byte of the header differs.
    let start = plain.len() - 8;
    assert_eq!((plain[12], bytes[12]), (3, 7));
    assert_eq!(bytes[13..start], plain[13..start]);

    for len in start..bytes.len() {
        let error = read(&bytes[..len]).unwrap_err();
        assert_eq!(error.format(), Format::Native);
    }
    for index in start..bytes.len() - 8 {
        for flip in [0x01, 0xFF] {
            let mut corrupted = bytes.clone();
            corrupted[index] ^= flip;
            // Any outcome but a panic: a flipped byte inside a string still
            // decodes.
            let _ = read(&resealed(corrupted));
        }
    }

    let cases = [
        // The block name's length marks an absent string.
        (
            start,
            u32::MAX.to_le_bytes().to_vec(),
            "A required string is absent.",
        ),
        // The issue kind of `CoordinateModelIgnored`.
        (
            find(&bytes, b"QQQQ") - 5,
            vec![99],
            "The interpretation issue kind 99 is unknown.",
        ),
        // The width of the coordinate models' indices.
        (
            find(&bytes, b"MODELX") + 6,
            vec![3],
            "The integer width 3 is not 1, 2, 4, or 8 bytes.",
        ),
        // The coordinate model index of the first instance.
        (
            find(&bytes, b"MODELX") + 7,
            vec![5],
            "The string index 5 is out of range.",
        ),
    ];
    for (offset, patch, message) in cases {
        let mut corrupted = bytes.clone();
        corrupted[offset..offset + patch.len()].copy_from_slice(&patch);
        let error = read(&resealed(corrupted)).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::InvalidValue, "{error}");
        assert!(error.to_string().contains(message), "{error}");
    }
}

#[test]
fn interpreted_structure_sets_need_one_report_per_geometry() {
    let (_, set) = parts();
    let reports = vec![report(false)];
    let error = native::write(
        NativeSource::interpreted_structure_set(&set, &reports),
        NativeWriteOptions::default(),
    )
    .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.to_string(),
        "Cannot write the native output. The structure set has 2 geometries and 1 interpretation \
         report. Give one report per geometry."
    );
    let interpretation = StructureSetInterpretation {
        set,
        reports: Vec::new(),
    };
    assert!(native::write(&interpretation, NativeWriteOptions::default()).is_err());
}
