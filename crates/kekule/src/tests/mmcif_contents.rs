use super::{deterministic_text_mutations, test_positions};
use crate::core::{Atom, AtomId, BondOrder, Element};
use crate::geometry::Point3;
use crate::mmcif::{
    self, MmcifAltLocPolicy, MmcifEntityClassifications, MmcifEntityKind, MmcifInstanceProvenance,
    MmcifInterpretOptions, MmcifInterpretationReport, MmcifModelSelection, MmcifParseOptions,
    MmcifWriteError, MmcifWriteOptions,
};
use crate::properties::{PropertyColumn, PropertyKey, PropertyValue};
use crate::structure::{Model, ModelBuildError, ModelBuilder, Positions};
use crate::topology::AtomSiteMetadata;
use crate::topology::{InstanceAtomId, MoleculeInstanceId, TopologyBuildError};
use crate::units::{Quantity, DIMENSIONLESS, NANOMETER, SQUARE_NANOMETER};

const MIXED: &str = r#"
data_mixed
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 water
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
L 2
W 3
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 N N GLY A 1 1 0.0 0.0 0.0 1
ATOM 2 C CA GLY A 1 1 1.0 0.0 0.0 1
HETATM 3 C C1 LIG L 2 . 2.0 0.0 0.0 1
HETATM 4 O O HOH W 3 . 3.0 0.0 0.0 1
loop_
_audit_author.name
_audit_author.pdbx_ordinal
'Example Author' 1
"#;

const SPLIT_SOURCE_CHAIN: &str = r#"
data_split_source_chain
loop_
_entity.id
_entity.type
7 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 7
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA GLY A 7 1 0.0 0.0 0.0
ATOM 2 C CA ALA A 7 2 10.0 0.0 0.0
"#;

const MULTI_MODEL: &str = r#"
data_multi
loop_
_entity.id
_entity.type
1 non-polymer
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 LIG A 1 0.0 0.0 0.0 0.40 10.0 1
HETATM 2 O O1 LIG A 1 1.2 0.0 0.0 0.50 11.0 1
HETATM 3 C C1 LIG A 1 5.0 0.0 0.0 0.80 20.0 2
HETATM 4 O O1 LIG A 1 6.2 0.0 0.0 0.90 21.0 2
"#;

const MULTI_BLOCK: &str = r#"
data_FIRST
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
C C1 ONE A 0.0 0.0 0.0

data_SECOND
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
N N1 TWO B 5.0 0.0 0.0
O O1 THREE C 10.0 0.0 0.0
"#;

const REPEATED_RESIDUE_ENSEMBLE: &str = r#"
data_repeated_residues
loop_
_entity.id
_entity.type
1 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA GLY A X 1 7 100 . . 0.0 0.0 0.0 1
ATOM 2 C CA GLY A X 1 7 100 A . 2.0 0.0 0.0 1
ATOM 3 C CA GLY A X 1 8 101 . . 4.0 0.0 0.0 1
ATOM 4 C CA GLY A X 1 7 100 . . 10.0 0.0 0.0 2
ATOM 5 C CA GLY A X 1 7 100 A . 12.0 0.0 0.0 2
ATOM 6 C CA GLY A X 1 8 101 . . 14.0 0.0 0.0 2
"#;

const REPEATED_NONPOLYMER_ENSEMBLE: &str = r#"
data_repeated_nonpolymer
loop_
_entity.id
_entity.type
1 non-polymer
loop_
_struct_asym.id
_struct_asym.entity_id
L 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 LIG L X 1 0.0 0.0 0.0 1
HETATM 2 O O1 LIG L X 1 1.2 0.0 0.0 1
HETATM 3 C C1 LIG L X 1 5.0 0.0 0.0 1
HETATM 4 O O1 LIG L X 1 6.2 0.0 0.0 1
HETATM 5 C C1 LIG L X 1 10.0 0.0 0.0 2
HETATM 6 O O1 LIG L X 1 11.2 0.0 0.0 2
HETATM 7 C C1 LIG L X 1 15.0 0.0 0.0 2
HETATM 8 O O1 LIG L X 1 16.2 0.0 0.0 2
"#;

const EXTREME_COORDINATE: &str = r#"
data_extreme
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG A 1.936908127739503e19 0.0 0.0
"#;

const CONNECTION_ATOMS: &str = r#"
data_connection_resolution
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
L 2
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.occupancy
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 N N GLY A N GLY X 1 1 10 . . 1.0 0.0 0.0 0.0 1
ATOM 2 C CA GLY A CA GLY X 1 1 10 . . 1.0 2.0 0.0 0.0 1
ATOM 3 N N ALA A N ALA X 1 2 11 . . 1.0 4.0 0.0 0.0 1
ATOM 4 C CA ALA A CA ALA X 1 2 11 . . 1.0 6.0 0.0 0.0 1
ATOM 5 S SG CYS A SG CYS X 1 3 20 A . 1.0 8.0 0.0 0.0 1
ATOM 6 S SG CYS A SG CYS X 1 4 20 B . 1.0 10.0 0.0 0.0 1
ATOM 7 O OG SER A OG SER X 1 5 30 . A 0.8 12.0 0.0 0.0 1
ATOM 8 O OG SER A OG SER X 1 5 30 . B 0.2 14.0 0.0 0.0 1
HETATM 9 C C1 LIG L C7 LIG Y 2 . 50 . . 1.0 20.0 0.0 0.0 1
HETATM 10 O O1 LIG L O7 LIG Y 2 . 50 . . 1.0 22.0 0.0 0.0 1
"#;

const AMBIGUOUS_CONNECTION_FUZZ_SEED: &str = r#"
data_ambiguous_connection
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA GLY A 1 0.0 0.0 0.0
ATOM 2 C CA ALA A 2 2.0 0.0 0.0
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
ambiguous covale A CA A CA
"#;

const AUTH_ONLY_CONNECTION_FUZZ_SEED: &str = r#"
data_auth_connection
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG L C7 LIG X 50 0.0 0.0 0.0
HETATM 2 O O1 LIG L O7 LIG X 50 2.0 0.0 0.0
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_auth_atom_id
auth covale X LIG 50 C7 X LIG 50 O7
"#;

fn parse(input: &str) -> mmcif::MmcifDocument {
    mmcif::parse_str(input).expect("mmCIF parses")
}

fn connection_input(atom_sites: &str, tags: &str, values: &str) -> String {
    format!("{atom_sites}\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n{tags}\n{values}\n")
}

fn connection_atom(
    interpretation: &mmcif::MmcifInterpretation,
    atom_name: &str,
    label_sequence_id: Option<i32>,
    insertion_code: Option<&str>,
) -> crate::topology::InstanceAtomId {
    interpretation
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .find(|atom| {
            atom.atom_name() == atom_name
                && atom.label_sequence_id() == label_sequence_id
                && atom.insertion_code() == insertion_code
        })
        .unwrap_or_else(|| {
            panic!(
                "missing provenance atom {atom_name} at label sequence {label_sequence_id:?} insertion {insertion_code:?}"
            )
        })
        .atom()
}

fn assert_declared_bond(
    interpretation: &mmcif::MmcifInterpretation,
    left: crate::topology::InstanceAtomId,
    right: crate::topology::InstanceAtomId,
    order: BondOrder,
) {
    assert_eq!(left.molecule(), right.molecule());
    let definition = interpretation
        .model()
        .topology()
        .definition_for_instance(left.molecule())
        .expect("connection molecule definition");
    let bond = definition
        .molecule()
        .bond_between(left.atom(), right.atom())
        .expect("valid connection atoms")
        .expect("declared connection bond");
    assert_eq!(
        definition
            .molecule()
            .bond(bond)
            .expect("declared bond")
            .order,
        order
    );
}

#[test]
fn deterministic_mmcif_parser_and_interpreters_fuzz_smoke_are_panic_free() {
    for seed in [
        MIXED,
        MULTI_MODEL,
        EXTREME_COORDINATE,
        AMBIGUOUS_CONNECTION_FUZZ_SEED,
        AUTH_ONLY_CONNECTION_FUZZ_SEED,
    ] {
        for input in deterministic_text_mutations(seed) {
            std::panic::catch_unwind(|| {
                let Ok(document) = mmcif::parse_str_with_options(
                    &input,
                    MmcifParseOptions {
                        max_input_bytes: 64 * 1024,
                        max_tokens: 16 * 1024,
                        max_token_bytes: 16 * 1024,
                        max_atom_site_rows: 4 * 1024,
                    },
                ) else {
                    return;
                };
                let _ = mmcif::interpret(
                    &document,
                    MmcifInterpretOptions {
                        model_selection: MmcifModelSelection::First,
                        ..MmcifInterpretOptions::default()
                    },
                );
                let _ = mmcif::interpret_ensemble(
                    &document,
                    mmcif::MmcifEnsembleInterpretOptions::default(),
                );
                let _ = mmcif::interpret_ensemble(
                    &document,
                    mmcif::MmcifEnsembleInterpretOptions {
                        model_ids: Some(Vec::new()),
                        ..mmcif::MmcifEnsembleInterpretOptions::default()
                    },
                );
            })
            .expect("mmCIF parser or interpreter smoke mutation panicked");
        }
    }
}

#[test]
fn mmcif_parse_preserves_unknown_categories_without_chemistry() {
    let document = parse(MIXED);
    let block = &document.blocks()[0];
    assert!(block.loop_with_tag("_audit_author.name").is_some());
    assert_eq!(
        block
            .loop_with_tag("_atom_site.type_symbol")
            .unwrap()
            .row_count(),
        4
    );
}

#[test]
fn mmcif_parse_preserves_multiple_named_blocks() {
    let document = parse(MULTI_BLOCK);

    assert_eq!(document.blocks().len(), 2);
    assert_eq!(document.blocks()[0].name(), "FIRST");
    assert_eq!(document.blocks()[1].name(), "SECOND");
    assert_eq!(
        document.block("FIRST").map(|block| block.name()),
        Some("FIRST")
    );
    assert_eq!(
        document.block("second").map(|block| block.name()),
        Some("SECOND")
    );
}

#[test]
fn sibling_blocks_are_interpreted_independently() {
    let document = parse(MULTI_BLOCK);
    let first = mmcif::interpret_block(
        document.block("FIRST").expect("first block"),
        MmcifInterpretOptions::default(),
    )
    .expect("first block interprets");
    let second = mmcif::interpret_block(
        document.block("SECOND").expect("second block"),
        MmcifInterpretOptions::default(),
    )
    .expect("second block interprets");

    assert_eq!(first.report().block_name(), "FIRST");
    assert_eq!(second.report().block_name(), "SECOND");
    assert_eq!(first.model().atom_count(), 1);
    assert_eq!(second.model().atom_count(), 2);
    assert_eq!(first.model().topology().instance_count(), 1);
    assert_eq!(second.model().topology().instance_count(), 2);
    assert_eq!(first.model().positions().values().value()[0].x, 0.0);
    assert!((second.model().positions().values().value()[0].x - 0.5).abs() < 1.0e-15);
}

#[test]
fn block_interpreters_reject_a_block_without_atom_site_data() {
    let document = parse("data_METADATA\n_entry.id metadata-only\n");
    let block = document.block("METADATA").expect("metadata block");

    let model_error = mmcif::interpret_block(block, MmcifInterpretOptions::default())
        .expect_err("metadata block has no model");
    assert!(model_error.message().contains("no atom-site loop"));
    assert!(matches!(
        mmcif::interpret_ensemble_block(block, mmcif::MmcifEnsembleInterpretOptions::default()),
        Err(mmcif::MmcifEnsembleInterpretError::NoCoordinateModels)
    ));
}

#[test]
fn document_model_interpretation_delegates_for_exactly_one_atom_site_block() {
    let document = parse(MIXED);
    let from_document =
        mmcif::interpret(&document, MmcifInterpretOptions::default()).expect("document interprets");
    let from_block =
        mmcif::interpret_block(&document.blocks()[0], MmcifInterpretOptions::default())
            .expect("block interprets");

    assert_eq!(from_document.report(), from_block.report());
    assert!(from_document.topology().same_layout(from_block.topology()));
    assert_eq!(
        from_document.model().positions(),
        from_block.model().positions()
    );
    assert_eq!(
        from_document.model().properties(),
        from_block.model().properties()
    );
}

#[test]
fn mmcif_parse_preserves_hashes_inside_bare_values() {
    let document = parse("data_hash\nloop_\n_example.id\n_example.label\n1 sample-d2o#1\n#\n");
    let table = document.blocks()[0]
        .loop_with_tag("_example.label")
        .expect("example loop");
    assert_eq!(table.row_count(), 1);
    assert_eq!(
        table.value(0, "_example.label").map(|value| value.text()),
        Some("sample-d2o#1")
    );
}

#[test]
fn interpretation_preserves_first_source_occurrence_for_model_instances() {
    let input = r#"
data_order
loop_
_entity.id
_entity.type
1 polymer
2 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
Z 1
A 2
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.seq_id
Z 1
A 1
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
covale Z C1 1 A C1 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C C1 GLY A 2 1 0.0 0.0 0.0 1
ATOM 2 C C1 GLY Z 1 1 1.0 0.0 0.0 1
"#;
    let result = mmcif::interpret(&parse(input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.model().topology().instance_count(), 1);
    let hierarchy = result.model().topology().hierarchy();
    let chain_order = hierarchy
        .chains()
        .map(|(_, chain)| chain.label_id())
        .collect::<Vec<_>>();
    assert_eq!(chain_order, ["Z", "A"]);
    assert_eq!(hierarchy.atom_sites().count(), 2);
    assert!(hierarchy
        .atom_sites()
        .all(|(_, site)| site.atom().molecule() == MoleculeInstanceId::new(0)));

    let written = mmcif::write_with_report(
        result.model(),
        result.report(),
        MmcifWriteOptions::default(),
    )
    .expect("one connected molecule may retain two source asymmetries");
    let round_trip = mmcif::interpret(&parse(&written), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(round_trip.topology().instance_count(), 1);
    assert_eq!(
        round_trip
            .topology()
            .hierarchy()
            .chains()
            .map(|(_, chain)| chain.label_id())
            .collect::<Vec<_>>(),
        ["Z", "A"]
    );
    assert_eq!(round_trip.topology().bond_count(), 1);
}

#[test]
fn interpretation_builds_connected_typed_instances_and_complete_positions() {
    let interpreted = mmcif::interpret(&parse(MIXED), MmcifInterpretOptions::default()).unwrap();
    let model = interpreted.model();
    assert_eq!(model.topology().instance_count(), 4);
    assert_eq!(model.atom_count(), 4);
    assert_eq!(model.positions().len(), 4);
    assert!(model
        .positions()
        .values()
        .value()
        .iter()
        .all(|point| point.x.is_finite()));
    assert!(!model.topology().hierarchy().is_empty());
    let hierarchy = model.topology().hierarchy();
    assert_eq!(
        hierarchy
            .chains()
            .map(|(_, chain)| chain.label_id())
            .collect::<Vec<_>>(),
        ["A", "L", "W"]
    );
    let chain_a = hierarchy
        .chains()
        .find(|(_, chain)| chain.label_id() == "A")
        .unwrap()
        .1;
    let mut chain_a_instances = chain_a
        .residues()
        .iter()
        .flat_map(|residue| hierarchy.residue(*residue).unwrap().atom_sites())
        .map(|site| hierarchy.atom_site(*site).unwrap().atom().molecule())
        .collect::<Vec<_>>();
    chain_a_instances.sort_unstable();
    chain_a_instances.dedup();
    assert_eq!(chain_a_instances.len(), 2);
    for (label, component) in [("L", "LIG"), ("W", "HOH")] {
        let chain = hierarchy
            .chains()
            .find(|(_, chain)| chain.label_id() == label)
            .unwrap()
            .1;
        assert_eq!(chain.residues().len(), 1);
        assert_eq!(
            hierarchy.residue(chain.residues()[0]).unwrap().name(),
            component
        );
    }
    assert!(
        model
            .topology()
            .definitions()
            .nth(2)
            .unwrap()
            .1
            .molecule()
            .atom_count()
            > 0
    );
    assert_eq!(interpreted.report().selected_model.as_deref(), Some("1"));
    assert_eq!(interpreted.report().instances.len(), 4);
    assert_eq!(
        interpreted
            .report()
            .instances
            .iter()
            .map(|instance| instance.atoms.len())
            .sum::<usize>(),
        4
    );
    assert!(interpreted.report().instances[0]
        .entity_kinds()
        .contains(&MmcifEntityKind::Polymer));
    assert!(interpreted.report().instances[1]
        .entity_kinds()
        .contains(&MmcifEntityKind::Polymer));
    assert!(interpreted.report().instances[2]
        .entity_kinds()
        .contains(&MmcifEntityKind::NonPolymer));
    assert!(interpreted.report().instances[3]
        .entity_kinds()
        .contains(&MmcifEntityKind::Water));
    for (_, definition) in model.topology().definitions() {
        assert!(definition.molecule().properties().is_empty());
        assert!(definition
            .molecule()
            .properties()
            .atoms()
            .iter()
            .all(|(key, _)| !key.as_str().starts_with("mmcif.")));
    }
    let first_provenance = &interpreted.report().instances[0];
    assert_eq!(first_provenance.coordinate_model_id, "1");
    assert_eq!(first_provenance.asym_ids, vec!["A"]);
    assert_eq!(first_provenance.entity_ids, vec!["1"]);
    assert_eq!(first_provenance.atoms[0].atom_name, "N");
    assert_eq!(model.topology().bonds().count(), 0);
    assert_eq!(interpreted.report().connectivity_candidates(), 1);
    assert!(interpreted.report().issues().iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::ConnectivityCandidatesInferred {
            candidate_count: 1,
            ..
        }
    )));
}

#[test]
fn interpretation_preserves_distinct_label_and_author_hierarchy_identity() {
    let input = r#"
data_identity
loop_
_entity.id
_entity.type
1 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA CAX GLY GLC A X 1 7 100 B 0.0 0.0 0.0
"#;
    let result = mmcif::interpret(&parse(input), MmcifInterpretOptions::default()).unwrap();
    let hierarchy = result.model().topology().hierarchy();
    let (_, chain) = hierarchy.chains().next().unwrap();
    assert_eq!(chain.label_id(), "A");
    assert_eq!(chain.author_id(), Some("X"));
    let residue = hierarchy.residue(chain.residues()[0]).unwrap();
    assert_eq!(residue.label_comp_id(), Some("GLY"));
    assert_eq!(residue.author_comp_id(), Some("GLC"));
    assert_eq!(residue.label_seq_id(), Some(7));
    assert_eq!(residue.author_seq_id(), Some("100"));
    assert_eq!(residue.insertion_code(), Some("B"));
    let site = hierarchy.atom_site(residue.atom_sites()[0]).unwrap();
    assert_eq!(site.metadata().label_atom_id.as_deref(), Some("CA"));
    assert_eq!(site.metadata().auth_atom_id.as_deref(), Some("CAX"));
    assert_eq!(site.metadata().label_asym_id.as_deref(), Some("A"));
    assert_eq!(site.metadata().auth_asym_id.as_deref(), Some("X"));
}

#[test]
fn auth_only_hierarchy_identity_is_preserved_then_deterministically_label_normalized() {
    let input = r#"
data_auth_only
loop_
_entity.id
_entity.type
1 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
X 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C . CAX . GLC . X 1 . 100 0.0 0.0 0.0
"#;
    let interpreted = mmcif::interpret(&parse(input), MmcifInterpretOptions::default()).unwrap();
    let atom = &interpreted.report().instances()[0].atoms()[0];
    assert_eq!(atom.label_asym_id(), None);
    assert_eq!(atom.auth_asym_id(), Some("X"));
    assert_eq!(atom.label_atom_name(), None);
    assert_eq!(atom.auth_atom_name(), Some("CAX"));
    assert_eq!(atom.label_component_id(), None);
    assert_eq!(atom.auth_component_id(), Some("GLC"));
    let hierarchy = interpreted.topology().hierarchy();
    let (_, chain) = hierarchy.chains().next().unwrap();
    assert_eq!(chain.label_id(), "X");
    assert_eq!(chain.author_id(), Some("X"));
    let residue = hierarchy.residue(chain.residues()[0]).unwrap();
    assert_eq!(residue.label_comp_id(), None);
    assert_eq!(residue.author_comp_id(), Some("GLC"));
    let site = hierarchy.atom_site(residue.atom_sites()[0]).unwrap();
    assert_eq!(site.metadata().label_asym_id, None);
    assert_eq!(site.metadata().label_atom_id, None);

    let written = mmcif::write_with_report(
        interpreted.model(),
        interpreted.report(),
        MmcifWriteOptions::default(),
    )
    .unwrap();
    let normalized = mmcif::interpret(&parse(&written), MmcifInterpretOptions::default()).unwrap();
    let atom = &normalized.report().instances()[0].atoms()[0];
    assert_eq!(atom.label_asym_id(), Some("X"));
    assert_eq!(atom.label_atom_name(), Some("CAX"));
    assert_eq!(atom.label_component_id(), Some("GLC"));
    assert_eq!(atom.auth_asym_id(), Some("X"));
    assert_eq!(atom.auth_atom_name(), Some("CAX"));
    assert_eq!(atom.auth_component_id(), Some("GLC"));
}

#[test]
fn canonical_label_residue_rejects_conflicting_author_aliases() {
    let input = r#"
data_conflicting_residue_alias
loop_
_entity.id
_entity.type
1 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY GLC A X 1 7 100 0.0 0.0 0.0
ATOM 2 C CA GLY GLC A X 1 7 101 10.0 0.0 0.0
"#;
    let error = mmcif::interpret(&parse(input), MmcifInterpretOptions::default()).unwrap_err();
    assert!(error
        .message()
        .contains("conflicting _atom_site.auth_seq_id"));

    let conflicting_label_component = input.replace(
        "ATOM 2 C CA GLY GLC A X 1 7 101",
        "ATOM 2 C CA ALA GLC A X 1 7 100",
    );
    let error = mmcif::interpret(
        &parse(&conflicting_label_component),
        MmcifInterpretOptions::default(),
    )
    .unwrap_err();
    assert!(error.message().contains("canonical residue identity"));

    let conflicting_author_component = input.replace(
        "ATOM 2 C CA GLY GLC A X 1 7 101",
        "ATOM 2 C CA GLY ALC A X 1 7 100",
    );
    let error = mmcif::interpret(
        &parse(&conflicting_author_component),
        MmcifInterpretOptions::default(),
    )
    .unwrap_err();
    assert!(error
        .message()
        .contains("conflicting _atom_site.auth_comp_id"));

    let conflicting_chain_alias = input.replace(
        "ATOM 2 C CA GLY GLC A X 1 7 101",
        "ATOM 2 C CA GLY GLC A Y 1 7 100",
    );
    let error = mmcif::interpret(
        &parse(&conflicting_chain_alias),
        MmcifInterpretOptions::default(),
    )
    .unwrap_err();
    assert!(error
        .message()
        .contains("conflicting _atom_site.auth_asym_id"));
}

#[test]
fn canonical_label_residue_uses_one_compatible_author_alias() {
    let input = r#"
data_compatible_residue_alias
loop_
_entity.id
_entity.type
1 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY GLC A X 1 7 . 0.0 0.0 0.0
ATOM 2 C CA GLY GLC A X 1 7 100 10.0 0.0 0.0
"#;
    let interpreted = mmcif::interpret(&parse(input), MmcifInterpretOptions::default()).unwrap();
    let hierarchy = interpreted.topology().hierarchy();
    let (_, chain) = hierarchy.chains().next().unwrap();
    assert_eq!(chain.residues().len(), 1);
    let residue = hierarchy.residue(chain.residues()[0]).unwrap();
    assert_eq!(residue.label_seq_id(), Some(7));
    assert_eq!(residue.author_seq_id(), Some("100"));
}

#[test]
fn unsequenced_component_identity_distinguishes_source_residues() {
    let input = r#"
data_unsequenced_components
loop_
_entity.id
_entity.type
1 non-polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG A 1 0.0 0.0 0.0
HETATM 2 C C1 ION A 1 10.0 0.0 0.0
"#;
    let interpreted = mmcif::interpret(&parse(input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(interpreted.topology().instance_count(), 2);
    let hierarchy = interpreted.topology().hierarchy();
    let (_, chain) = hierarchy.chains().next().unwrap();
    assert_eq!(chain.residues().len(), 2);
    assert_eq!(
        chain
            .residues()
            .iter()
            .map(|residue| hierarchy.residue(*residue).unwrap().name())
            .collect::<Vec<_>>(),
        ["LIG", "ION"]
    );
}

#[test]
fn disconnected_waters_share_one_source_hierarchy_chain() {
    let input = r#"
data_waters
loop_
_entity.id
_entity.type
1 water
loop_
_struct_asym.id
_struct_asym.entity_id
W 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 O O HOH W 1 0.0 0.0 0.0
HETATM 2 O O HOH W 1 5.0 0.0 0.0
"#;
    let result = mmcif::interpret(&parse(input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.model().topology().instance_count(), 2);
    let hierarchy = result.model().topology().hierarchy();
    let (_, chain) = hierarchy.chains().next().unwrap();
    assert_eq!(hierarchy.chains().count(), 1);
    assert_eq!(chain.label_id(), "W");
    assert_eq!(chain.residues().len(), 2);
    assert_eq!(hierarchy.atom_sites().count(), 2);

    let written = mmcif::write_with_report(
        result.model(),
        result.report(),
        MmcifWriteOptions::default(),
    )
    .unwrap();
    let round_trip = mmcif::interpret(&parse(&written), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(round_trip.topology().instance_count(), 2);
    assert_eq!(round_trip.topology().hierarchy().chains().count(), 1);
    assert!(round_trip
        .report()
        .instances()
        .iter()
        .all(|instance| instance.entity_kinds() == [MmcifEntityKind::Water]));
}

#[test]
fn mmcif_connectivity_candidates_do_not_create_bonds() {
    let input = MIXED.replace(
        "ATOM 2 C CA GLY A 1 1 1.0 0.0 0.0 1",
        "ATOM 2 C CA GLY A 1 1 1.45 0.0 0.0 1\nATOM 5 C C GLY A 1 1 2.90 0.0 0.0 1",
    );
    let interpreted = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(interpreted.model().topology().instance_count(), 5);
    assert_eq!(interpreted.model().topology().atom_count(), 5);
    assert_eq!(interpreted.model().topology().bonds().count(), 0);
    assert!(interpreted
        .model()
        .topology()
        .definitions()
        .all(|(_, definition)| definition.molecule().validate_connected().is_ok()));
    assert_eq!(interpreted.report().connectivity_candidates(), 2);
}

#[test]
fn connectivity_diagnostic_cell_boundaries_are_finite_and_panic_free() {
    const CELL_ANGSTROM: f64 = 2.1;
    let near_positive = ((i64::MAX as f64) - 4_096.0) * CELL_ANGSTROM;
    let near_negative = ((i64::MIN as f64) + 4_096.0) * CELL_ANGSTROM;
    let beyond_positive = (i64::MAX as f64) * CELL_ANGSTROM;
    let beyond_negative = (i64::MIN as f64) * CELL_ANGSTROM;

    let with_coordinate = |coordinate: f64| {
        MIXED.replace(
            "ATOM 1 N N GLY A 1 1 0.0 0.0 0.0 1",
            &format!("ATOM 1 N N GLY A 1 1 {coordinate:e} 0.0 0.0 1"),
        )
    };

    for coordinate in [near_positive, near_negative] {
        let document = parse(&with_coordinate(coordinate));
        let result = std::panic::catch_unwind(|| {
            mmcif::interpret(&document, MmcifInterpretOptions::default())
        });
        assert!(
            matches!(result, Ok(Ok(_))),
            "supported finite coordinate {coordinate:e} must interpret without panicking"
        );
    }

    for coordinate in [beyond_positive, beyond_negative] {
        let document = parse(&with_coordinate(coordinate));
        let result = std::panic::catch_unwind(|| {
            mmcif::interpret(&document, MmcifInterpretOptions::default())
        });
        let error = result
            .expect("out-of-range finite coordinate must not panic")
            .expect_err("out-of-range finite coordinate must be rejected");
        assert!(error.line().is_some());
        assert!(error
            .message()
            .contains("supported covalent-connectivity diagnostic cell range"));
    }
}

#[test]
fn multiple_coordinate_models_require_explicit_selection() {
    let input = MIXED.replace(
        "HETATM 4 O O HOH W 3 . 3.0 0.0 0.0 1",
        "HETATM 4 O O HOH W 3 . 3.0 0.0 0.0 1\nHETATM 5 O O HOH W 3 . 8.0 0.0 0.0 2",
    );
    let document = parse(&input);
    let error = mmcif::interpret(&document, MmcifInterpretOptions::default()).unwrap_err();
    assert!(error.message.contains("select one explicitly"));

    let selected = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::Select("2".into()),
            ..MmcifInterpretOptions::default()
        },
    )
    .unwrap();
    assert_eq!(selected.report().selected_model.as_deref(), Some("2"));
    assert_eq!(selected.model().atom_count(), 1);
    assert!((selected.model().positions().values().value()[0].x - 0.8).abs() < 1.0e-15);
    assert_eq!(selected.report().ignored_coordinate_models, vec!["1"]);

    let first = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..MmcifInterpretOptions::default()
        },
    )
    .unwrap();
    assert_eq!(first.report().selected_model.as_deref(), Some("1"));
}

#[test]
fn ensemble_interpretation_rejects_empty_explicit_model_selection_without_panicking() {
    let document = parse(MULTI_MODEL);
    let result = std::panic::catch_unwind(|| {
        mmcif::interpret_ensemble(
            &document,
            mmcif::MmcifEnsembleInterpretOptions {
                model_ids: Some(Vec::new()),
                ..mmcif::MmcifEnsembleInterpretOptions::default()
            },
        )
    });
    assert!(matches!(
        result,
        Ok(Err(mmcif::MmcifEnsembleInterpretError::EmptyModelSelection))
    ));
}

#[test]
fn document_interpreters_reject_multiple_atom_site_blocks() {
    let document = parse(MULTI_BLOCK);
    let single_error = mmcif::interpret(&document, MmcifInterpretOptions::default())
        .expect_err("ambiguous blocks");
    assert!(single_error
        .message()
        .contains("atom-site data in more than one data block"));
    assert!(matches!(
        mmcif::interpret_ensemble(&document, mmcif::MmcifEnsembleInterpretOptions::default(),),
        Err(mmcif::MmcifEnsembleInterpretError::MultipleAtomSiteBlocks)
    ));
}

#[test]
fn ensemble_block_interpretation_matches_exactly_one_document_helper() {
    let document = parse(MULTI_MODEL);
    let from_document =
        mmcif::interpret_ensemble(&document, mmcif::MmcifEnsembleInterpretOptions::default())
            .expect("document ensemble interprets");
    let from_block = mmcif::interpret_ensemble_block(
        &document.blocks()[0],
        mmcif::MmcifEnsembleInterpretOptions::default(),
    )
    .expect("block ensemble interprets");
    let method_default = document.blocks()[0]
        .interpret_ensemble()
        .expect("default block method interprets");
    let method_explicit = document.blocks()[0]
        .interpret_ensemble_with_options(mmcif::MmcifEnsembleInterpretOptions::default())
        .expect("explicit block method interprets");
    let document_method = document
        .interpret_ensemble()
        .expect("document method delegates to its structural block");

    assert_eq!(from_document.reports(), from_block.reports());
    assert_eq!(method_default.reports(), method_explicit.reports());
    assert_eq!(document_method.reports(), method_default.reports());
    assert!(method_default
        .topology()
        .same_layout(method_default.ensemble().topology()));
    assert!(method_default
        .clone()
        .to_ensemble()
        .topology()
        .same_layout(method_default.topology()));
    assert!(from_document
        .ensemble()
        .topology()
        .same_layout(from_block.ensemble().topology()));
    assert_eq!(from_document.ensemble().len(), from_block.ensemble().len());
    for (document_member, block_member) in from_document
        .ensemble()
        .members()
        .zip(from_block.ensemble().members())
    {
        assert_eq!(document_member.positions(), block_member.positions());
        assert_eq!(document_member.properties(), block_member.properties());
    }
}

#[test]
fn multimodel_interpretation_builds_shared_topology_with_distinct_properties() {
    let interpreted = mmcif::interpret_ensemble(
        &parse(MULTI_MODEL),
        mmcif::MmcifEnsembleInterpretOptions::default(),
    )
    .expect("consistent coordinate models form an ensemble");
    let ensemble = interpreted.ensemble();
    assert_eq!(ensemble.len(), 2);
    assert_eq!(interpreted.reports().len(), 2);
    let shared_topology = ensemble.shared_topology();
    assert!(ensemble.members().all(|member| {
        member.positions().len() == shared_topology.atom_count()
            && member.atom_properties().len() == shared_topology.atom_count()
            && member.bond_properties().len() == shared_topology.bond_count()
    }));
    let first_positions = ensemble
        .members()
        .map(|member| member.positions().values().value()[0].x)
        .collect::<Vec<_>>();
    assert!((first_positions[0] - 0.0).abs() < 1.0e-15);
    assert!((first_positions[1] - 0.5).abs() < 1.0e-15);

    assert_eq!(interpreted.reports()[0].selected_model(), Some("1"));
    assert_eq!(interpreted.reports()[1].selected_model(), Some("2"));
    let member_properties = ensemble.members().collect::<Vec<_>>();
    let occupancy = PropertyKey::new("occupancy").unwrap();
    let b_factor = PropertyKey::new("b_factor").unwrap();
    for (member, expected_occupancies, expected_b_factors) in [
        (member_properties[0], [0.4, 0.5], [0.1, 0.11]),
        (member_properties[1], [0.8, 0.9], [0.2, 0.21]),
    ] {
        for index in 0..2 {
            assert_eq!(
                member.atom_properties().value(&occupancy, index).unwrap(),
                Some(PropertyValue::Real {
                    value: expected_occupancies[index],
                    unit: DIMENSIONLESS
                })
            );
            let Some(PropertyValue::Real { value, unit }) =
                member.atom_properties().value(&b_factor, index).unwrap()
            else {
                panic!("B-factor")
            };
            assert_eq!(unit, SQUARE_NANOMETER);
            assert!((value - expected_b_factors[index]).abs() < 1.0e-15);
        }
    }
}

#[test]
fn multimodel_interpretation_rejects_inconsistent_atom_sets() {
    let inconsistent = MULTI_MODEL.replace("HETATM 4 O O1 LIG A 1 6.2 0.0 0.0 0.90 21.0 2\n", "");
    assert!(matches!(
        mmcif::interpret_ensemble(
            &parse(&inconsistent),
            mmcif::MmcifEnsembleInterpretOptions::default(),
        ),
        Err(mmcif::MmcifEnsembleInterpretError::InconsistentAtomSet { model_id })
            if model_id == "2"
    ));
}

#[test]
fn ensemble_identity_distinguishes_repeated_residues_and_insertion_codes() {
    let interpreted = mmcif::interpret_ensemble(
        &parse(REPEATED_RESIDUE_ENSEMBLE),
        mmcif::MmcifEnsembleInterpretOptions::default(),
    )
    .expect("stable repeated-residue identities form an ensemble");

    assert_eq!(interpreted.ensemble().len(), 2);
    assert_eq!(interpreted.ensemble().topology().atom_count(), 3);
    let atoms = interpreted.reports()[0]
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .collect::<Vec<_>>();
    assert_eq!(atoms.len(), 3);
    assert_eq!(atoms[0].label_sequence_id(), Some(7));
    assert_eq!(atoms[0].author_sequence_id(), Some("100"));
    assert_eq!(atoms[0].insertion_code(), None);
    assert_eq!(atoms[0].auth_asym_id(), Some("X"));
    assert_eq!(atoms[0].occurrence(), None);
    assert_eq!(atoms[1].label_sequence_id(), Some(7));
    assert_eq!(atoms[1].insertion_code(), Some("A"));
    assert_eq!(atoms[2].label_sequence_id(), Some(8));
}

#[test]
fn ensemble_identity_preserves_repeated_nonpolymer_occurrences() {
    let interpreted = mmcif::interpret_ensemble(
        &parse(REPEATED_NONPOLYMER_ENSEMBLE),
        mmcif::MmcifEnsembleInterpretOptions::default(),
    )
    .expect("stable occurrence discriminators form an ensemble");

    assert_eq!(interpreted.ensemble().len(), 2);
    assert_eq!(interpreted.ensemble().topology().instance_count(), 4);
    for report in interpreted.reports() {
        let occurrences = report
            .instances()
            .iter()
            .flat_map(|instance| instance.atoms())
            .map(|atom| atom.occurrence())
            .collect::<Vec<_>>();
        assert_eq!(occurrences, [Some(0), Some(0), Some(1), Some(1)]);
    }
}

#[test]
fn ensemble_identity_classifies_reordered_rows_as_dense_order_mismatch() {
    let reordered = r#"
data_reordered_instances
loop_
_entity.id
_entity.type
1 non-polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 LIG A XA 1 0.0 0.0 0.0 1
HETATM 2 C C1 LIG B XB 1 5.0 0.0 0.0 1
HETATM 3 C C1 LIG B XB 1 15.0 0.0 0.0 2
HETATM 4 C C1 LIG A XA 1 10.0 0.0 0.0 2
"#;
    assert!(matches!(
        mmcif::interpret_ensemble(
            &parse(reordered),
            mmcif::MmcifEnsembleInterpretOptions::default(),
        ),
        Err(mmcif::MmcifEnsembleInterpretError::InconsistentDenseAtomOrder { model_id })
            if model_id == "2"
    ));
}

#[test]
fn ensemble_identity_detects_true_repeated_atom_set_mismatch() {
    let mismatched = REPEATED_RESIDUE_ENSEMBLE.replace(
        "ATOM 6 C CA GLY A X 1 8 101 . . 14.0 0.0 0.0 2",
        "ATOM 6 C CA GLY A X 1 9 101 . . 14.0 0.0 0.0 2",
    );
    assert!(matches!(
        mmcif::interpret_ensemble(
            &parse(&mismatched),
            mmcif::MmcifEnsembleInterpretOptions::default(),
        ),
        Err(mmcif::MmcifEnsembleInterpretError::InconsistentAtomSet { model_id })
            if model_id == "2"
    ));
}

#[test]
fn ensemble_identity_excludes_selected_alternate_location_but_reports_preserve_it() {
    let input = r#"
data_altloc_models
loop_
_entity.id
_entity.type
1 polymer
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.label_alt_id
_atom_site.occupancy
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA GLY A 1 1 A 0.8 0.0 0.0 0.0 1
ATOM 2 C CA GLY A 1 1 B 0.2 1.0 0.0 0.0 1
ATOM 3 C CA GLY A 1 1 A 0.1 10.0 0.0 0.0 2
ATOM 4 C CA GLY A 1 1 B 0.9 11.0 0.0 0.0 2
"#;
    let interpreted = mmcif::interpret_ensemble(
        &parse(input),
        mmcif::MmcifEnsembleInterpretOptions::default(),
    )
    .expect("selected altloc is provenance rather than atom identity");
    let ensemble = interpreted.ensemble();
    assert_eq!(ensemble.len(), 2);
    assert_eq!(ensemble.topology().atom_count(), 1);
    let selected_positions = ensemble
        .members()
        .map(|member| member.positions().values().value()[0].x)
        .collect::<Vec<_>>();
    assert!((selected_positions[0] - 0.0).abs() < 1.0e-15);
    assert!((selected_positions[1] - 1.1).abs() < 1.0e-15);
    let topology = ensemble.shared_topology();
    let atom = topology.atom_ids()[0];
    assert_eq!(
        interpreted.reports()[0].instances()[0].atoms()[0].atom(),
        interpreted.reports()[1].instances()[0].atoms()[0].atom()
    );
    assert_eq!(
        ensemble
            .member(0)
            .unwrap()
            .as_model()
            .occupancy(atom)
            .unwrap(),
        Some(0.8)
    );
    assert_eq!(
        ensemble
            .member(1)
            .unwrap()
            .as_model()
            .occupancy(atom)
            .unwrap(),
        Some(0.9)
    );
    let selected_altlocs = interpreted
        .reports()
        .iter()
        .map(|report| {
            report.instances()[0].atoms()[0]
                .selected_alternate_location()
                .unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(selected_altlocs, ["A", "B"]);
}

#[test]
fn alternate_location_policy_is_explicit_and_reported() {
    let input = MIXED
        .replace(
            "_atom_site.Cartn_x",
            "_atom_site.label_alt_id\n_atom_site.occupancy\n_atom_site.Cartn_x",
        )
        .replace(
            "ATOM 1 N N GLY A 1 1 0.0 0.0 0.0 1",
            "ATOM 1 N N GLY A 1 1 A 0.4 0.0 0.0 0.0 1\nATOM 5 N N GLY A 1 1 B 0.6 5.0 0.0 0.0 1",
        )
        .replace(" 1.0 0.0 0.0 1", " . 1.0 1.0 0.0 0.0 1")
        .replace(" 2.0 0.0 0.0 1", " . 1.0 2.0 0.0 0.0 1")
        .replace(" 3.0 0.0 0.0 1", " . 1.0 3.0 0.0 0.0 1");
    let document = parse(&input);
    let result = mmcif::interpret(&document, MmcifInterpretOptions::default()).unwrap();
    assert!((result.model().positions().values().value()[0].x - 0.5).abs() < 1.0e-15);
    assert!(result.report().issues.iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::AlternateLocationOmitted { alt_id: Some(id), .. } if id == "A"
    )));
    assert!(mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            altloc_policy: MmcifAltLocPolicy::ErrorOnAlternateLocations,
            ..MmcifInterpretOptions::default()
        }
    )
    .is_err());
}

#[test]
fn selected_model_requires_complete_positions() {
    let input = MIXED.replace("3.0 0.0 0.0 1", ". . . 1");
    let error = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap_err();
    assert!(error.message.contains("complete position"));
}

#[test]
fn struct_conn_label_sequence_distinguishes_repeated_atom_names() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_value_order"#,
        "label-sequence covale A GLY 1 CA A ALA 2 CA doub",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections(), 1);
    assert_declared_bond(
        &result,
        connection_atom(&result, "CA", Some(1), None),
        connection_atom(&result, "CA", Some(2), None),
        BondOrder::Double,
    );
}

#[test]
fn struct_conn_author_insertion_code_distinguishes_repeated_author_sequence() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_atom_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_auth_atom_id
_struct_conn.pdbx_ptnr2_PDB_ins_code"#,
        "insertion covale X CYS 20 SG A X CYS 20 SG B",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections(), 1);
    assert_declared_bond(
        &result,
        connection_atom(&result, "SG", Some(3), Some("A")),
        connection_atom(&result, "SG", Some(4), Some("B")),
        BondOrder::Single,
    );
}

#[test]
fn struct_conn_author_fields_resolve_when_label_sequence_is_absent() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_auth_atom_id"#,
        "author-only covale Y LIG 50 C7 Y LIG 50 O7",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections(), 1);
    assert_declared_bond(
        &result,
        connection_atom(&result, "C1", None, None),
        connection_atom(&result, "O1", None, None),
        BondOrder::Single,
    );
}

#[test]
fn struct_conn_consistent_label_and_author_selectors_resolve_same_atom() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_auth_atom_id"#,
        "consistent covale A GLY 1 N X GLY 10 N A GLY 1 CA X GLY 10 CA",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections(), 1);
    assert_declared_bond(
        &result,
        connection_atom(&result, "N", Some(1), None),
        connection_atom(&result, "CA", Some(1), None),
        BondOrder::Single,
    );
}

#[test]
fn struct_conn_conflicting_label_and_author_selectors_are_unresolved() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id"#,
        "conflict covale A GLY 1 N X ALA 11 CA A GLY 1 CA",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections(), 0);
    assert!(result.report().issues().iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::ConnectionUnresolved {
            connection_id: Some(id),
            connection_type,
            partner: 1,
            source_line: Some(_),
            reason:
                mmcif::MmcifConnectionResolutionReason::ConflictingLabelAndAuthorSelectors,
        } if id == "conflict" && connection_type == "covale"
    )));
}

#[test]
fn struct_conn_under_specified_selector_is_reported_ambiguous() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id"#,
        "ambiguous covale A CA A GLY 1 N",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections(), 0);
    assert!(result.report().issues().iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::ConnectionAmbiguous {
            connection_id: Some(id),
            connection_type,
            partner: 1,
            source_line: Some(_),
            candidates: 2,
            reason: mmcif::MmcifConnectionResolutionReason::MultipleMatchingAtoms,
        } if id == "ambiguous" && connection_type == "covale"
    )));
}

#[test]
fn struct_conn_zero_candidate_selector_is_reported_unresolved() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id"#,
        "missing covale Z CA A GLY 1 N",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections(), 0);
    assert!(result.report().issues().iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::ConnectionUnresolved {
            connection_id: Some(id),
            partner: 1,
            reason: mmcif::MmcifConnectionResolutionReason::NoMatchingAtom,
            ..
        } if id == "missing"
    )));
}

#[test]
fn struct_conn_named_omitted_altloc_does_not_bind_retained_altloc() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_alt_id
_struct_conn.pdbx_ptnr1_label_alt_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_auth_atom_id
_struct_conn.pdbx_ptnr1_auth_alt_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id"#,
        "omitted-alt covale A SER 5 OG B B X SER 30 OG B A GLY 1 N",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections(), 0);
    assert!(result.report().issues().iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::ConnectionUnresolved {
            connection_id: Some(id),
            partner: 1,
            reason:
                mmcif::MmcifConnectionResolutionReason::AlternateLocationOmitted {
                    alternate_location
                },
            ..
        } if id == "omitted-alt" && alternate_location == "B"
    )));
}

#[test]
fn struct_conn_resolution_is_independent_of_atom_site_row_order() {
    let tags = r#"_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id"#;
    let values = "reordered covale A GLY 1 CA A ALA 2 CA";
    let reordered_atoms = CONNECTION_ATOMS.replace(
        "ATOM 1 N N GLY A N GLY X 1 1 10 . . 1.0 0.0 0.0 0.0 1\nATOM 2 C CA GLY A CA GLY X 1 1 10 . . 1.0 2.0 0.0 0.0 1\nATOM 3 N N ALA A N ALA X 1 2 11 . . 1.0 4.0 0.0 0.0 1\nATOM 4 C CA ALA A CA ALA X 1 2 11 . . 1.0 6.0 0.0 0.0 1",
        "ATOM 4 C CA ALA A CA ALA X 1 2 11 . . 1.0 6.0 0.0 0.0 1\nATOM 3 N N ALA A N ALA X 1 2 11 . . 1.0 4.0 0.0 0.0 1\nATOM 2 C CA GLY A CA GLY X 1 1 10 . . 1.0 2.0 0.0 0.0 1\nATOM 1 N N GLY A N GLY X 1 1 10 . . 1.0 0.0 0.0 0.0 1",
    );
    let original = mmcif::interpret(
        &parse(&connection_input(CONNECTION_ATOMS, tags, values)),
        MmcifInterpretOptions::default(),
    )
    .unwrap();
    let reordered = mmcif::interpret(
        &parse(&connection_input(&reordered_atoms, tags, values)),
        MmcifInterpretOptions::default(),
    )
    .unwrap();
    for interpretation in [&original, &reordered] {
        assert_eq!(interpretation.report().applied_connections(), 1);
        assert_declared_bond(
            interpretation,
            connection_atom(interpretation, "CA", Some(1), None),
            connection_atom(interpretation, "CA", Some(2), None),
            BondOrder::Single,
        );
    }
}

#[test]
fn ordinary_unambiguous_struct_conn_retains_current_behavior() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_value_order"#,
        "ordinary covale A CYS 3 SG A GLY 1 N sing",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections(), 1);
    assert_declared_bond(
        &result,
        connection_atom(&result, "SG", Some(3), Some("A")),
        connection_atom(&result, "N", Some(1), None),
        BondOrder::Single,
    );
}

#[test]
fn missing_struct_conn_value_order_defaults_to_single_bond() {
    let input = connection_input(
        CONNECTION_ATOMS,
        r#"_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id"#,
        "default-order covale A GLY 1 N A GLY 1 CA",
    );
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_declared_bond(
        &result,
        connection_atom(&result, "N", Some(1), None),
        connection_atom(&result, "CA", Some(1), None),
        BondOrder::Single,
    );
}

#[test]
fn declared_covalent_links_merge_entities_but_noncovalent_links_do_not() {
    let connections = r#"
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
covale A CA 1 L C1 .
hydrog A N 1 W O .
"#;
    let input = format!("{MIXED}\n{connections}");
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.model().topology().instance_count(), 3);
    let merged = result
        .report()
        .instances()
        .iter()
        .find(|instance| {
            instance.entity_kinds().contains(&MmcifEntityKind::Polymer)
                && instance
                    .entity_kinds()
                    .contains(&MmcifEntityKind::NonPolymer)
        })
        .expect("covalently linked entities should share one instance");
    let merged_id = merged.molecule();
    assert!(result
        .model()
        .topology()
        .molecule(merged_id)
        .unwrap()
        .atom_sites()
        .next()
        .is_some());
    assert_eq!(result.report().applied_connections, 1);

    let written = mmcif::write_with_report(
        result.model(),
        result.report(),
        MmcifWriteOptions::default(),
    )
    .expect("one connected molecule may span distinct source entities and asymmetries");
    let round_trip = mmcif::interpret(&parse(&written), MmcifInterpretOptions::default()).unwrap();
    assert!(round_trip.report().instances().iter().any(|instance| {
        instance.entity_kinds().contains(&MmcifEntityKind::Polymer)
            && instance
                .entity_kinds()
                .contains(&MmcifEntityKind::NonPolymer)
    }));
    assert_eq!(round_trip.report().applied_connections(), 1);
}

#[test]
fn symmetry_mate_connections_are_reported_unresolved() {
    let connection = r#"
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_symmetry
symmetry-link disulf A N 1 1_555 A N 1 15_545
"#;
    let input = format!("{MIXED}\n{connection}");
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections, 0);
    assert!(result.report().issues.iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::ConnectionUnresolved {
            connection_id: Some(id),
            connection_type,
            partner: 2,
            source_line: Some(_),
            reason: mmcif::MmcifConnectionResolutionReason::UnsupportedSymmetry {
                symmetry
            },
        } if id == "symmetry-link" && connection_type == "disulf" && symmetry == "15_545"
    )));
}

#[test]
fn struct_conn_bond_order_is_interpreted_and_rejected_when_unknown() {
    let connection = r#"
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.pdbx_value_order
covale A N 1 A CA 1 doub
"#;
    let input = format!("{MIXED}\n{connection}");
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    let first = result.model().topology().definitions().next().unwrap().1;
    assert_eq!(
        first
            .molecule()
            .bonds()
            .next()
            .expect("declared bond")
            .1
            .order,
        BondOrder::Double
    );

    let error = mmcif::interpret(
        &parse(&input.replace("doub", "arom")),
        MmcifInterpretOptions::default(),
    )
    .unwrap_err();
    assert!(error.line().is_some());
    assert!(error
        .message
        .contains("unsupported struct_conn bond order `arom`"));
}

#[test]
fn mmcif_writer_round_trips_supported_model_content() {
    let connection = r#"
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.pdbx_value_order
covale A N 1 A CA 1 doub
"#;
    let (mut original, report) = mmcif::interpret(
        &parse(&format!("{MIXED}\n{connection}")),
        MmcifInterpretOptions::default(),
    )
    .unwrap()
    .to_parts();
    let topology = original.shared_topology();
    let first_atom = topology.atom_ids()[0];
    original.set_occupancy(first_atom, Some(0.625)).unwrap();
    original
        .set_b_factor(first_atom, Some(Quantity::new(0.125, NANOMETER.powi(2))))
        .unwrap();
    original
        .insert_atom_property_column(
            PropertyKey::new("analysis_score").unwrap(),
            PropertyColumn::Real {
                unit: DIMENSIONLESS,
                values: vec![Some(3.0), None, None, None],
            },
        )
        .unwrap();
    let written = mmcif::write_with_report(
        &original,
        &report,
        MmcifWriteOptions {
            block_name: "round_trip".to_owned(),
            coordinate_precision: 4,
        },
    )
    .expect("supported model should write");
    assert!(written.starts_with("data_round_trip\n"));
    assert!(written.contains("_struct_conn.pdbx_value_order"));
    assert!(written.contains("doub"));

    let document = parse(&written);
    let atom_sites = document.blocks()[0]
        .loop_with_tag("_atom_site.type_symbol")
        .expect("writer emits atom-site loop");
    assert_eq!(atom_sites.row_count(), 4);
    let round_trip = mmcif::interpret(&document, MmcifInterpretOptions::default()).unwrap();
    assert_eq!(round_trip.model().topology().instance_count(), 3);
    assert_eq!(
        round_trip.model().positions().values(),
        original.positions().values()
    );
    for atom in original.topology().atom_ids() {
        assert_eq!(
            round_trip.model().occupancy(*atom).unwrap(),
            original.occupancy(*atom).unwrap()
        );
        assert_eq!(
            round_trip.model().b_factor(*atom).unwrap(),
            original.b_factor(*atom).unwrap()
        );
    }
    assert!(round_trip
        .model()
        .properties()
        .atoms()
        .get(&PropertyKey::new("analysis_score").unwrap())
        .is_none());
    let (first_id, _) = round_trip.model().topology().instances().next().unwrap();
    let first = round_trip
        .model()
        .topology()
        .definition_for_instance(first_id)
        .unwrap();
    assert!(round_trip.report().instances()[0]
        .entity_kinds()
        .contains(&MmcifEntityKind::Polymer));
    assert_eq!(
        first
            .molecule()
            .bonds()
            .next()
            .expect("round-trip bond")
            .1
            .order,
        BondOrder::Double
    );
}

#[test]
fn mmcif_writer_keeps_one_source_asym_entity_across_disconnected_molecules() {
    let interpreted =
        mmcif::interpret(&parse(SPLIT_SOURCE_CHAIN), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(interpreted.topology().instance_count(), 2);
    assert_eq!(interpreted.topology().hierarchy().chains().count(), 1);
    assert!(interpreted
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .all(|atom| {
            atom.entity_id() == Some("7") && atom.entity_kind() == &MmcifEntityKind::Polymer
        }));

    let written = mmcif::write_with_report(
        interpreted.model(),
        interpreted.report(),
        MmcifWriteOptions::default(),
    )
    .expect("a source asymmetry may span disconnected molecule instances");
    let document = parse(&written);
    let block = &document.blocks()[0];
    let struct_asym = block
        .loop_with_tag("_struct_asym.id")
        .expect("writer emits structural asymmetries");
    assert_eq!(struct_asym.row_count(), 1);
    let asym_entity = struct_asym
        .value(0, "_struct_asym.entity_id")
        .unwrap()
        .optional_text()
        .unwrap();
    let atom_sites = block
        .loop_with_tag("_atom_site.label_entity_id")
        .expect("writer emits atom sites");
    assert!((0..atom_sites.row_count()).all(|row| {
        atom_sites
            .value(row, "_atom_site.label_entity_id")
            .and_then(|value| value.optional_text())
            == Some(asym_entity)
    }));

    let round_trip = mmcif::interpret(&document, MmcifInterpretOptions::default()).unwrap();
    assert_eq!(round_trip.topology().instance_count(), 2);
    assert_eq!(round_trip.topology().hierarchy().chains().count(), 1);
    assert_eq!(round_trip.topology().bond_count(), 0);
    assert!(round_trip
        .report()
        .instances()
        .iter()
        .all(|instance| { instance.entity_kinds() == [MmcifEntityKind::Polymer] }));

    let mut generic = MmcifEntityClassifications::new();
    for (molecule, _) in interpreted.topology().instances() {
        generic.insert(molecule, MmcifEntityKind::Polymer).unwrap();
    }
    let generic_written = mmcif::write_with_classifications(
        interpreted.model(),
        &generic,
        MmcifWriteOptions::default(),
    )
    .expect("consistent generic classifications represent a cross-molecule chain");
    let generic_round_trip =
        mmcif::interpret(&parse(&generic_written), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(generic_round_trip.topology().instance_count(), 2);
    assert_eq!(
        generic_round_trip.topology().hierarchy().chains().count(),
        1
    );

    let mut conflicting = MmcifEntityClassifications::new();
    let molecules = interpreted
        .topology()
        .instances()
        .map(|(molecule, _)| molecule)
        .collect::<Vec<_>>();
    conflicting
        .insert(molecules[0], MmcifEntityKind::Polymer)
        .unwrap();
    conflicting
        .insert(molecules[1], MmcifEntityKind::NonPolymer)
        .unwrap();
    assert!(matches!(
        mmcif::write_with_classifications(
            interpreted.model(),
            &conflicting,
            MmcifWriteOptions::default(),
        ),
        Err(MmcifWriteError::ConflictingAsymEntityClassifications {
            asym_id,
            ..
        }) if asym_id == "A"
    ));

    let mut contradictory_report = interpreted.report().clone();
    contradictory_report.instances[1].atoms[0].entity_id = Some("8".to_owned());
    assert!(matches!(
        mmcif::write_with_report(
            interpreted.model(),
            &contradictory_report,
            MmcifWriteOptions::default(),
        ),
        Err(MmcifWriteError::ConflictingAsymEntityIds { asym_id, .. })
            if asym_id == "A"
    ));
}

#[test]
fn mmcif_writer_rejects_unsupported_chemistry_and_topology_rejects_invalid_hierarchy() {
    let dative = small_model_with_bond(BondOrder::Dative);
    let classifications = classifications_for(&dative, MmcifEntityKind::NonPolymer);
    assert!(matches!(
        mmcif::write_with_classifications(&dative, &classifications, MmcifWriteOptions::default()),
        Err(MmcifWriteError::UnsupportedBondOrder {
            order: BondOrder::Dative,
            ..
        })
    ));

    let mut graph = crate::core::MoleculeEditor::new();
    let atom = graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .expect("atom identifier capacity");
    let molecule = graph.finish().unwrap();
    let mut builder = ModelBuilder::new();
    let instance = builder
        .add_molecule(&molecule, &Positions::zeros(1))
        .unwrap();
    let chain = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_chain("A", None)
        .unwrap();
    let residue = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, AtomId::new(99)),
            AtomSiteMetadata::default(),
        )
        .unwrap();
    assert!(matches!(
        builder.build(),
        Err(ModelBuildError::Topology(
            TopologyBuildError::InvalidHierarchy(_)
        ))
    ));
    let _ = atom;
}

#[test]
fn mmcif_writer_preserves_supported_bond_orders() {
    for order in [
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Quadruple,
    ] {
        let model = small_model_with_bond(order);
        let classifications = classifications_for(&model, MmcifEntityKind::NonPolymer);
        let written = mmcif::write_with_classifications(
            &model,
            &classifications,
            MmcifWriteOptions::default(),
        )
        .unwrap();
        let interpreted =
            mmcif::interpret(&parse(&written), MmcifInterpretOptions::default()).unwrap();
        let round_trip = interpreted
            .model()
            .topology()
            .definitions()
            .next()
            .unwrap()
            .1
            .molecule()
            .bonds()
            .next()
            .unwrap()
            .1
            .order;
        assert_eq!(round_trip, order);
    }
}

#[test]
fn mmcif_writer_rejects_ambiguous_atom_identity() {
    let carbon = Element::from_symbol("C").unwrap();
    let mut graph = crate::core::MoleculeEditor::new();
    let left = graph
        .add_atom(Atom::new(carbon))
        .expect("atom identifier capacity");
    let right = graph
        .add_atom(Atom::new(carbon))
        .expect("atom identifier capacity");
    graph
        .add_bond(left, right, BondOrder::Single)
        .expect("connected duplicate-identity fixture");
    let positions = test_positions(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]);
    let macro_molecule = graph.finish().unwrap();
    let mut builder = ModelBuilder::new();
    let instance = builder.add_molecule(&macro_molecule, &positions).unwrap();
    let chain = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_chain("A", None)
        .unwrap();
    let residue = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    for atom in [left, right] {
        builder
            .topology_builder_mut()
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, atom),
                AtomSiteMetadata {
                    label_atom_id: Some("CA".to_owned()),
                    ..AtomSiteMetadata::default()
                },
            )
            .unwrap();
    }
    let model = builder.build().unwrap();
    let report = report_with_entity_kinds(&model, &[MmcifEntityKind::NonPolymer]);
    assert!(matches!(
        mmcif::write_with_report(&model, &report, MmcifWriteOptions::default()),
        Err(MmcifWriteError::DuplicateAtomIdentity(_))
    ));
}

#[test]
fn mmcif_writer_requires_explicit_classification_for_hierarchy() {
    let model = hierarchical_single_atom_model("LIG", "C1", "C");
    let molecule = model.topology().instances().next().unwrap().0;
    assert_eq!(
        mmcif::write(&model, MmcifWriteOptions::default()),
        Err(MmcifWriteError::MissingEntityClassification(molecule))
    );

    let classifications = classifications_for(&model, MmcifEntityKind::NonPolymer);
    let written =
        mmcif::write_with_classifications(&model, &classifications, MmcifWriteOptions::default())
            .unwrap();
    assert!(written.contains("1 non-polymer"));
    assert!(!written.contains("1 polymer"));
}

#[test]
fn mmcif_writer_requires_explicit_classification_without_hierarchy() {
    let model = small_single_atom_model("C");
    let molecule = model.topology().instances().next().unwrap().0;
    assert_eq!(
        mmcif::write(&model, MmcifWriteOptions::default()),
        Err(MmcifWriteError::MissingEntityClassification(molecule))
    );
}

#[test]
fn mmcif_writer_does_not_infer_water_from_neutral_oxygen() {
    let model = small_single_atom_model("O");
    let molecule = model.topology().instances().next().unwrap().0;
    assert_eq!(
        mmcif::write(&model, MmcifWriteOptions::default()),
        Err(MmcifWriteError::MissingEntityClassification(molecule))
    );

    let classifications = classifications_for(&model, MmcifEntityKind::NonPolymer);
    let written =
        mmcif::write_with_classifications(&model, &classifications, MmcifWriteOptions::default())
            .unwrap();
    assert!(written.contains("1 non-polymer"));
    assert!(!written.contains("1 water"));

    let classifications = classifications_for(&model, MmcifEntityKind::Water);
    let written =
        mmcif::write_with_classifications(&model, &classifications, MmcifWriteOptions::default())
            .unwrap();
    assert!(written.contains("1 water"));
}

#[test]
fn mmcif_writer_with_report_preserves_every_supported_source_kind() {
    assert_report_kind(
        &hierarchical_single_atom_model("GLY", "CA", "C"),
        MmcifEntityKind::Polymer,
        "polymer",
    );
    assert_report_kind(
        &hierarchical_single_atom_model("NAG", "C1", "C"),
        MmcifEntityKind::Branched,
        "branched",
    );
    assert_report_kind(
        &small_single_atom_model("C"),
        MmcifEntityKind::NonPolymer,
        "non-polymer",
    );
    assert_report_kind(
        &small_single_atom_model("O"),
        MmcifEntityKind::Water,
        "water",
    );
}

#[test]
fn mmcif_writer_uses_explicit_polymer_and_branched_kinds() {
    let model = hierarchical_single_atom_model("NAG", "C1", "C");
    let classifications = classifications_for(&model, MmcifEntityKind::Branched);
    let written =
        mmcif::write_with_classifications(&model, &classifications, MmcifWriteOptions::default())
            .unwrap();
    assert!(written.contains("1 branched"));

    let classifications = classifications_for(&model, MmcifEntityKind::Polymer);
    let written =
        mmcif::write_with_classifications(&model, &classifications, MmcifWriteOptions::default())
            .unwrap();
    assert!(written.contains("1 polymer"));
}

#[test]
fn mmcif_writer_rejects_conflicting_missing_duplicate_and_unknown_classifications() {
    let model = hierarchical_single_atom_model("NAG", "C1", "C");
    let molecule = model.topology().instances().next().unwrap().0;
    let report = report_with_entity_kinds(
        &model,
        &[MmcifEntityKind::Polymer, MmcifEntityKind::Branched],
    );
    assert!(matches!(
        mmcif::write_with_report(&model, &report, MmcifWriteOptions::default()),
        Err(MmcifWriteError::ConflictingEntityClassifications {
            molecule: conflicted,
            ..
        }) if conflicted == molecule
    ));

    let report = report_with_entity_kinds(&model, &[]);
    assert_eq!(
        mmcif::write_with_report(&model, &report, MmcifWriteOptions::default()),
        Err(MmcifWriteError::MissingEntityClassification(molecule))
    );

    let report = report_with_entity_kinds(
        &model,
        &[MmcifEntityKind::Other("unsupported-kind".to_owned())],
    );
    assert_eq!(
        mmcif::write_with_report(&model, &report, MmcifWriteOptions::default()),
        Err(MmcifWriteError::UnsupportedEntityClassification {
            molecule,
            classification: "unsupported-kind".to_owned(),
        })
    );

    let mut report = report_with_entity_kinds(&model, &[MmcifEntityKind::Branched]);
    report.instances.push(report.instances[0].clone());
    assert_eq!(
        mmcif::write_with_report(&model, &report, MmcifWriteOptions::default()),
        Err(MmcifWriteError::DuplicateEntityClassification(molecule))
    );

    let mut duplicate = MmcifEntityClassifications::new();
    duplicate
        .insert(molecule, MmcifEntityKind::Branched)
        .unwrap();
    assert_eq!(
        duplicate.insert(molecule, MmcifEntityKind::Polymer),
        Err(MmcifWriteError::DuplicateEntityClassification(molecule))
    );

    let unknown = MoleculeInstanceId::new(u32::MAX);
    let mut classifications = MmcifEntityClassifications::new();
    classifications
        .insert(unknown, MmcifEntityKind::NonPolymer)
        .unwrap();
    assert_eq!(
        mmcif::write_with_classifications(&model, &classifications, MmcifWriteOptions::default()),
        Err(MmcifWriteError::UnknownClassifiedMolecule(unknown))
    );
}

#[test]
fn mmcif_writer_requires_one_classification_for_every_instance() {
    let interpreted = mmcif::interpret(&parse(MIXED), MmcifInterpretOptions::default()).unwrap();
    let model = interpreted.model();
    let mut instances = model.topology().instances().map(|(id, _)| id);
    let classified = instances.next().unwrap();
    let missing = instances.next().unwrap();
    let mut classifications = MmcifEntityClassifications::new();
    classifications
        .insert(classified, MmcifEntityKind::Polymer)
        .unwrap();
    assert_eq!(
        mmcif::write_with_classifications(model, &classifications, MmcifWriteOptions::default()),
        Err(MmcifWriteError::MissingEntityClassification(missing))
    );
}

fn classifications_for(model: &Model, kind: MmcifEntityKind) -> MmcifEntityClassifications {
    let mut classifications = MmcifEntityClassifications::new();
    for (molecule, _) in model.topology().instances() {
        classifications.insert(molecule, kind.clone()).unwrap();
    }
    classifications
}

fn assert_report_kind(model: &Model, kind: MmcifEntityKind, expected: &str) {
    let report = report_with_entity_kinds(model, &[kind]);
    let written = mmcif::write_with_report(model, &report, MmcifWriteOptions::default()).unwrap();
    assert!(written.contains(&format!("1 {expected}")));
}

fn report_with_entity_kinds(model: &Model, kinds: &[MmcifEntityKind]) -> MmcifInterpretationReport {
    MmcifInterpretationReport {
        instances: model
            .topology()
            .instances()
            .map(|(molecule, _)| MmcifInstanceProvenance {
                molecule,
                coordinate_model_id: "1".to_owned(),
                asym_ids: Vec::new(),
                entity_ids: Vec::new(),
                entity_kinds: kinds.to_vec(),
                atoms: Vec::new(),
            })
            .collect(),
        ..MmcifInterpretationReport::default()
    }
}

fn small_single_atom_model(element: &str) -> Model {
    let mut graph = crate::core::MoleculeEditor::new();
    graph
        .add_atom(Atom::new(Element::from_symbol(element).unwrap()))
        .unwrap();
    let positions = Positions::zeros(1);
    let molecule = graph.finish().unwrap();
    let mut builder = ModelBuilder::new();
    builder.add_molecule(&molecule, &positions).unwrap();
    builder.build().unwrap()
}

fn hierarchical_single_atom_model(component: &str, atom_name: &str, element: &str) -> Model {
    let mut graph = crate::core::MoleculeEditor::new();
    let atom = graph
        .add_atom(Atom::new(Element::from_symbol(element).unwrap()))
        .unwrap();
    let positions = Positions::zeros(1);
    let molecule = graph.finish().unwrap();
    let mut builder = ModelBuilder::new();
    let instance = builder.add_molecule(&molecule, &positions).unwrap();
    let chain = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_chain("A", None)
        .unwrap();
    let residue = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_residue(chain, component, Some(1), None, None)
        .unwrap();
    builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, atom),
            AtomSiteMetadata {
                label_atom_id: Some(atom_name.to_owned()),
                ..AtomSiteMetadata::default()
            },
        )
        .unwrap();
    builder.build().unwrap()
}

fn small_model_with_bond(order: BondOrder) -> Model {
    let mut graph = crate::core::MoleculeEditor::new();
    let left = graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .expect("atom identifier capacity");
    let right = graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .expect("atom identifier capacity");
    graph.add_bond(left, right, order).unwrap();
    let positions = test_positions(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]);
    let molecule = graph.finish().unwrap();
    let mut builder = ModelBuilder::new();
    builder.add_molecule(&molecule, &positions).unwrap();
    builder.build().unwrap()
}
