use crate::core::PropValue;
use crate::mmcif::{
    self, MmcifAltLocPolicy, MmcifEntry, MmcifInterpretIssue, MmcifInterpretOptions,
    MmcifParseOptions,
};

use super::deterministic_text_mutations;

fn interpret(source: &str) -> mmcif::MmcifInterpretation {
    let document = mmcif::parse_str(source, MmcifParseOptions::default()).expect("parse mmCIF");
    mmcif::interpret(&document, MmcifInterpretOptions::default()).expect("interpret mmCIF")
}

#[test]
fn mmcif_document_preserves_unknown_content_and_multiple_blocks() {
    let source = r#"
data_first
_entry.id demo
_custom.control 'loop_'
_custom.quoted_missing '.'
_custom.missing .
loop_
_custom.id
_custom.value
1 '_looks_like_a_tag'
2
;
multiline
value
;
stop_

data_second
_other.value 42
"#;
    let document = mmcif::parse_str(source, MmcifParseOptions::default()).expect("document");

    assert_eq!(document.blocks().len(), 2);
    let first = document
        .block("FIRST")
        .expect("case-insensitive block lookup");
    assert_eq!(first.item("_entry.id").expect("entry").text(), "demo");
    assert_eq!(
        first.item("_CUSTOM.CONTROL").expect("custom item").text(),
        "loop_"
    );
    assert_eq!(
        first
            .item("_custom.quoted_missing")
            .expect("quoted dot")
            .optional_text(),
        Some(".")
    );
    assert!(first
        .item("_custom.missing")
        .expect("missing dot")
        .is_missing());
    let table = first.loop_with_tag("_custom.value").expect("custom loop");
    assert_eq!(table.row_count(), 2);
    assert_eq!(
        table.value(0, "_custom.value").expect("row zero").text(),
        "_looks_like_a_tag"
    );
    assert_eq!(
        table.value(1, "_custom.value").expect("row one").text(),
        "multiline\nvalue"
    );
    assert!(matches!(first.entries()[4], MmcifEntry::Loop(_)));
    assert_eq!(
        document
            .block("second")
            .and_then(|block| block.item("_other.value"))
            .expect("second block item")
            .text(),
        "42"
    );
}

#[test]
fn mmcif_document_rejects_malformed_structure() {
    for source in [
        "_entry.id before_data",
        "data_x\n_entry.id",
        "data_x\nloop_\n_a.id\n_a.value\n1",
        "data_x\n_entry.id 1\n_ENTRY.ID 2",
        "data_x\nloop_\n_a.id\n_A.ID\n1 2",
        "data_x\n_entry.id 1\nDATA_X\n_entry.id 2",
    ] {
        assert!(
            mmcif::parse_str(source, MmcifParseOptions::default()).is_err(),
            "source should fail: {source:?}"
        );
    }

    let limited = MmcifParseOptions {
        max_atom_site_rows: 1,
        ..MmcifParseOptions::default()
    };
    assert!(mmcif::parse_str(
        "DATA_X\nLOOP_\n_atom_site.type_symbol\nC\nO\nSTOP_",
        limited
    )
    .is_err());
}

#[test]
fn mmcif_document_resource_limits_cover_the_canonical_reader() {
    let input = "data_x\nloop_\n_atom_site.type_symbol\nC\n";

    for options in [
        MmcifParseOptions {
            max_input_bytes: input.len() - 1,
            ..MmcifParseOptions::default()
        },
        MmcifParseOptions {
            max_tokens: 2,
            ..MmcifParseOptions::default()
        },
        MmcifParseOptions {
            max_token_bytes: 2,
            ..MmcifParseOptions::default()
        },
        MmcifParseOptions {
            max_atom_site_rows: 0,
            ..MmcifParseOptions::default()
        },
    ] {
        assert!(mmcif::parse_str(input, options).is_err());
    }
}

#[test]
fn deterministic_mmcif_document_mutations_are_panic_free() {
    let seed = "data_tiny\nloop_\n_entity.id\n_entity.type\n1 non-polymer\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\nC C1 LIG A 1\n";
    for input in deterministic_text_mutations(seed) {
        std::panic::catch_unwind(|| {
            if let Ok(document) = mmcif::parse_str(&input, MmcifParseOptions::default()) {
                let _ = mmcif::interpret(&document, MmcifInterpretOptions::default());
            }
        })
        .expect("mmCIF parse-then-interpret mutation panicked");
    }
}

const MIXED_CONTENTS: &str = r#"
data_mixed
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 non-polymer
4 water
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
C 3
W 4
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_model_num
_atom_site.pdbx_formal_charge
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM N N ALA A 1 1 10 1 ? 0 0 0
ATOM C CA ALA A 1 1 10 1 ? 1 0 0
HETATM C C1 LIG B 2 . 501 1 ? 2 0 0
HETATM O O1 LIG B 2 . 501 1 ? 3 0 0
HETATM Mg MG MG C 3 . 601 1 2 4 0 0
HETATM O O HOH W 4 . 701 1 ? 5 0 0
HETATM O O HOH W 4 . 702 1 ? 6 0 0
"#;

#[test]
fn interpretation_separates_macro_small_ion_and_solvent_molecules() {
    let interpretation = interpret(MIXED_CONTENTS);
    let contents = interpretation.contents();

    let macros = contents.macromolecules().collect::<Vec<_>>();
    assert_eq!(macros.len(), 1);
    assert_eq!(macros[0].graph().atom_count(), 2);
    assert_eq!(macros[0].hierarchy().chains().count(), 1);

    let small = contents.small_molecules().collect::<Vec<_>>();
    assert_eq!(small.len(), 2);
    let ligand = small
        .iter()
        .find(|molecule| molecule.atom_count() == 2)
        .expect("ligand");
    assert_eq!(ligand.bond_count(), 0);
    let ion = small
        .iter()
        .find(|molecule| molecule.atom_count() == 1)
        .expect("ion");
    assert_eq!(ion.atoms().next().expect("ion atom").1.formal_charge, 2);

    assert_eq!(contents.solvent().len(), 2);
    assert!(contents
        .solvent()
        .molecules()
        .all(|water| water.atom_count() == 1));

    let report = interpretation.report();
    assert_eq!(report.macromolecules, 1);
    assert_eq!(report.small_molecules, 2);
    assert_eq!(report.solvent_molecules, 2);
    assert_eq!(report.coordinate_models, 1);
    assert_eq!(report.template_bonds_pending, 2);

    let (owned_contents, owned_report) = interpretation.clone().into_parts();
    let (owned_small, owned_macro, owned_solvent) = owned_contents.into_parts();
    assert_eq!(owned_small.len(), 2);
    assert_eq!(owned_macro.len(), 1);
    assert_eq!(owned_solvent.into_molecules().len(), 2);
    assert_eq!(owned_report, *report);
}

#[test]
fn interpretation_preserves_macromolecular_hierarchy_metadata() {
    let interpretation = interpret(
        r#"
data_metadata
loop_
_entity.id
_entity.type
1 polymer
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA CAY . GLY GLY A X 1 10 42 A 0.50 12.25 1.25 2.50 3.75 7
ATOM 2 O O O . GLY GLY A X 1 10 42 A 1.00 10.00 4.25 5.50 6.75 7
"#,
    );
    let molecule = interpretation
        .contents()
        .macromolecules()
        .next()
        .expect("macromolecule");

    let (_, chain) = molecule.hierarchy().chains().next().expect("chain");
    assert_eq!(chain.label_id, "A");
    assert_eq!(chain.author_id.as_deref(), Some("X"));
    let (_, residue) = molecule.hierarchy().residues().next().expect("residue");
    assert_eq!(residue.label_comp_id.as_deref(), Some("GLY"));
    assert_eq!(residue.author_comp_id.as_deref(), Some("GLY"));
    assert_eq!(residue.label_seq_id, Some(10));
    assert_eq!(residue.author_seq_id.as_deref(), Some("42"));
    assert_eq!(residue.insertion_code.as_deref(), Some("A"));
    let (_, site) = molecule.hierarchy().atom_sites().next().expect("atom site");
    assert_eq!(site.metadata.label_atom_id.as_deref(), Some("CA"));
    assert_eq!(site.metadata.auth_atom_id.as_deref(), Some("CAY"));
    assert_eq!(site.metadata.occupancy, Some(0.5));
    assert_eq!(site.metadata.b_factor, Some(12.25));
    let (_, conformer) = molecule.graph().first_conformer().expect("conformer");
    assert_eq!(
        conformer.props().get("mmcif.model_id"),
        Some(&PropValue::String("7".into()))
    );
    assert_eq!(
        conformer.position(site.atom),
        Some(crate::core::Point3::new(1.25, 2.5, 3.75))
    );
}

#[test]
fn coordinate_models_become_conformers_without_duplicate_atoms() {
    let interpretation = interpret(
        r#"
data_models
loop_
_entity.id
_entity.type
1 non-polymer
loop_
_struct_asym.id
_struct_asym.entity_id
L 1
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_model_num
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
C C1 LIG L 1 5 1 0 0 0
O O1 LIG L 1 5 1 1 0 0
C C1 LIG L 1 5 2 0 1 0
O O1 LIG L 1 5 2 1 1 0
"#,
    );
    let molecule = interpretation
        .contents()
        .small_molecules()
        .next()
        .expect("molecule");

    assert_eq!(molecule.atom_count(), 2);
    let conformers = molecule.graph().conformers().collect::<Vec<_>>();
    assert_eq!(conformers.len(), 2);
    assert_eq!(
        conformers[0].1.props().get("mmcif.model_id"),
        Some(&PropValue::String("1".into()))
    );
    assert_eq!(
        conformers[1].1.props().get("mmcif.model_id"),
        Some(&PropValue::String("2".into()))
    );
    assert_eq!(interpretation.report().coordinate_models, 2);
}

const ALT_LOCS: &str = r#"
data_alt
loop_
_entity.id
_entity.type
1 non-polymer
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.auth_seq_id
_atom_site.label_alt_id
_atom_site.occupancy
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
C C1 LIG L 1 1 A 0.4 1 0 0
C C1 LIG L 1 1 B 0.6 2 0 0
"#;

#[test]
fn alternate_location_policy_is_explicit_and_deterministic() {
    let document = mmcif::parse_str(ALT_LOCS, MmcifParseOptions::default()).expect("document");
    let highest = mmcif::interpret(&document, MmcifInterpretOptions::default()).expect("highest");
    assert_eq!(first_x(highest.contents()), 2.0);

    let selected = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            altloc_policy: MmcifAltLocPolicy::SelectLabel("A".into()),
            ..MmcifInterpretOptions::default()
        },
    )
    .expect("select A");
    assert_eq!(first_x(selected.contents()), 1.0);

    let error = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            altloc_policy: MmcifAltLocPolicy::ErrorOnAlternateLocations,
            ..MmcifInterpretOptions::default()
        },
    )
    .expect_err("alternate locations should be rejected");
    assert!(error.message.contains("alternate locations"));

    let unavailable = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            altloc_policy: MmcifAltLocPolicy::SelectLabel("C".into()),
            ..MmcifInterpretOptions::default()
        },
    )
    .expect_err("missing alternate-location label");
    assert!(unavailable.message.contains("unavailable"));
}

#[test]
fn duplicate_atom_records_are_not_silently_treated_as_alternates() {
    let source = ALT_LOCS.replace("B 0.6", "A 0.6");
    let document = mmcif::parse_str(&source, MmcifParseOptions::default()).expect("document");
    let error = mmcif::interpret(&document, MmcifInterpretOptions::default())
        .expect_err("duplicate alternate label");
    assert!(error.message.contains("duplicate records"));
}

fn first_x(contents: &crate::bio::MolecularContents) -> f64 {
    let molecule = contents.small_molecules().next().expect("small molecule");
    let atom = molecule.graph().atom_ids().next().expect("atom");
    molecule
        .graph()
        .first_conformer()
        .expect("conformer")
        .1
        .position(atom)
        .expect("position")
        .x
}

#[test]
fn declared_covalent_connections_merge_molecular_instances() {
    let interpretation = interpret(
        r#"
data_linked
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM S SG CYS A 1 1 1 0 0 0
HETATM C C1 LIG B 2 . 9 1 0 0
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
covale A 1 SG B . C1
metalc A 1 SG B . C1
"#,
    );

    assert_eq!(interpretation.contents().macromolecules().count(), 1);
    assert_eq!(interpretation.contents().small_molecules().count(), 0);
    let molecule = interpretation
        .contents()
        .macromolecules()
        .next()
        .expect("linked macromolecule");
    assert_eq!(molecule.graph().atom_count(), 2);
    assert_eq!(molecule.graph().bond_count(), 1);
    assert_eq!(molecule.hierarchy().chains().count(), 2);
    assert_eq!(interpretation.report().applied_connections, 1);
    assert!(interpretation.report().issues.iter().any(|issue| matches!(
        issue,
        MmcifInterpretIssue::ConnectionIgnored { connection_type } if connection_type == "metalc"
    )));
}

#[test]
fn missing_entity_metadata_is_reported_or_rejected_in_strict_mode() {
    let source = r#"
data_inferred
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM Na NA NA I 1 0 0 0
"#;
    let document = mmcif::parse_str(source, MmcifParseOptions::default()).expect("document");
    let inferred = mmcif::interpret(&document, MmcifInterpretOptions::default()).expect("infer");
    assert!(matches!(
        inferred.report().issues.as_slice(),
        [MmcifInterpretIssue::EntityTypeInferred { asym_id, .. }] if asym_id == "I"
    ));

    let error = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            strict_entity_metadata: true,
            ..MmcifInterpretOptions::default()
        },
    )
    .expect_err("strict metadata");
    assert!(error.message.contains("missing entity type"));
}

#[test]
fn interpretation_rejects_ambiguous_documents_and_inconsistent_models() {
    let no_atoms = mmcif::parse_str("data_empty\n_entry.id x", MmcifParseOptions::default())
        .expect("format document");
    assert!(mmcif::interpret(&no_atoms, MmcifInterpretOptions::default()).is_err());

    let empty_atoms = mmcif::parse_str(
        "data_empty\nloop_\n_atom_site.type_symbol",
        MmcifParseOptions::default(),
    )
    .expect("empty atom table");
    assert!(mmcif::interpret(&empty_atoms, MmcifInterpretOptions::default()).is_err());

    let two_blocks = mmcif::parse_str(
        "data_a\nloop_\n_atom_site.type_symbol\nC\ndata_b\nloop_\n_atom_site.type_symbol\nO",
        MmcifParseOptions::default(),
    )
    .expect("two coordinate blocks");
    assert!(mmcif::interpret(&two_blocks, MmcifInterpretOptions::default()).is_err());

    let inconsistent = r#"
data_models
loop_
_entity.id
_entity.type
1 non-polymer
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_model_num
_atom_site.pdbx_formal_charge
C C1 LIG L 1 1 1 0
C C1 LIG L 1 1 2 1
"#;
    let document = mmcif::parse_str(inconsistent, MmcifParseOptions::default()).expect("document");
    let error = mmcif::interpret(&document, MmcifInterpretOptions::default())
        .expect_err("inconsistent charge");
    assert!(error.message.contains("inconsistent topology payload"));
}

#[test]
fn interpretation_rejects_invalid_atom_site_payloads() {
    for (source, expected) in [
        (
            "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\nC C1 LIG A 1 2",
            "partial atom-site coordinate",
        ),
        (
            "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\nXx C1 LIG A",
            "unknown atom-site element",
        ),
        (
            "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\nC C1 LIG A 999999999999999999",
            "invalid integer",
        ),
        (
            "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nC C1 LIG A 1e999 0 0",
            "non-finite float",
        ),
    ] {
        let document = mmcif::parse_str(source, MmcifParseOptions::default()).expect("document");
        let error = mmcif::interpret(&document, MmcifInterpretOptions::default())
            .expect_err("invalid atom-site payload");
        assert!(error.message.contains(expected), "{error}");
    }
}
