use crate::mmcif::{
    self, MmcifAltLocPolicy, MmcifInterpretOptions, MmcifModelSelection, MmcifParseOptions,
};
use crate::modeling::MoleculeRole;

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

fn parse(input: &str) -> mmcif::MmcifDocument {
    mmcif::parse_str(input, MmcifParseOptions::default()).expect("mmCIF parses")
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
fn interpretation_builds_distinct_typed_instances_and_complete_positions() {
    let interpreted = mmcif::interpret(&parse(MIXED), MmcifInterpretOptions::default()).unwrap();
    let model = interpreted.model();
    assert_eq!(model.topology().molecule_count(), 3);
    assert_eq!(model.atom_count(), 4);
    assert_eq!(model.positions().len(), 4);
    assert!(model.positions().iter().all(|point| point.x.is_finite()));
    let instances = model
        .topology()
        .molecules()
        .map(|(_, molecule)| molecule)
        .collect::<Vec<_>>();
    assert!(instances[0].macro_molecule().is_some());
    assert!(instances[0].has_role(MoleculeRole::Polymer));
    assert!(instances[1].small_molecule().is_some());
    assert!(instances[1].has_role(MoleculeRole::NonPolymer));
    assert!(instances[2].has_role(MoleculeRole::Solvent));
    assert_eq!(interpreted.report().selected_model.as_deref(), Some("1"));
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
    assert_eq!(selected.model().positions()[0].x, 8.0);
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
    assert_eq!(result.model().positions()[0].x, 5.0);
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
    assert_eq!(result.model().topology().molecule_count(), 2);
    let first = result.model().topology().molecules().next().unwrap().1;
    assert!(first.macro_molecule().is_some());
    assert!(first.has_role(MoleculeRole::Polymer));
    assert!(first.has_role(MoleculeRole::NonPolymer));
    assert_eq!(result.report().applied_connections, 1);
}
