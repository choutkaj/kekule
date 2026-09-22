use std::collections::BTreeMap;

use crate::mmcif::*;

fn document(rows: &str) -> MmcifDocument {
    parse_str(&format!(
        "data_alternates\nloop_\n\
_atom_site.id\n\
_atom_site.type_symbol\n\
_atom_site.label_atom_id\n\
_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n\
_atom_site.auth_seq_id\n\
_atom_site.pdbx_PDB_ins_code\n\
_atom_site.label_alt_id\n\
_atom_site.occupancy\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
_atom_site.pdbx_PDB_model_num\n{rows}\n"
    ))
    .unwrap()
}

// Synthetic coordinate examples are focused unit regressions, not benchmarks.
const CROSSED: &str = "\
1 C C1 LIG X 42 142 . A 0.9 0 0 0 1
2 C C1 LIG X 42 142 . B 0.1 10 0 0 1
3 O O1 LIG X 42 142 . A 0.2 1.2 0 0 1
4 O O1 LIG X 42 142 . B 0.8 11.2 0 0 1
5 N N1 LIG X 42 142 . . 1.0 2.4 0 0 1";

fn selected(
    document: &MmcifDocument,
    policy: MmcifAltLocPolicy,
) -> Result<MmcifInterpretation, MmcifInterpretError> {
    document.interpret_with_options(MmcifInterpretOptions {
        altloc_policy: policy,
        ..Default::default()
    })
}

fn choices(result: &MmcifInterpretation) -> Vec<Vec<String>> {
    result
        .report()
        .alternate_locations()
        .iter()
        .map(|choice| choice.selected_labels.clone())
        .collect()
}

#[test]
fn mmcif_altloc_default_keeps_residue_coherent_and_shared_atoms() {
    let source = document(CROSSED);
    let original = source.clone();
    let result = source.interpret().unwrap();
    assert_eq!(source, original);
    assert_eq!(
        source.blocks()[0]
            .loop_with_tag("_atom_site.id")
            .unwrap()
            .row_count(),
        5
    );
    assert_eq!(choices(&result), [vec!["A"]]);
    let atoms = result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .collect::<Vec<_>>();
    assert_eq!(
        atoms
            .iter()
            .map(|atom| atom.selected_alternate_location())
            .collect::<Vec<_>>(),
        [Some("A"), Some("A"), None]
    );
    let values = result.model().positions().values();
    let points = values.value();
    assert!((points[1].x - 0.12).abs() < 1e-12);
    let report = &result.report().alternate_locations()[0];
    assert_eq!(report.selected_source_lines.len(), 3);
    assert_eq!(report.omitted_source_lines.len(), 2);
}

#[test]
fn mmcif_altloc_exact_and_preferred_selection_are_distinct() {
    let source = document(CROSSED);
    let exact = selected(&source, MmcifAltLocPolicy::SelectLabel("B".into())).unwrap();
    assert_eq!(choices(&exact), [vec!["B"]]);
    assert_eq!(
        exact.report().alternate_locations()[0].reason,
        MmcifAltLocSelectionReason::ExactSelection
    );
    assert!(selected(&source, MmcifAltLocPolicy::SelectLabel("Z".into())).is_err());
    let preferred = selected(&source, MmcifAltLocPolicy::PreferLabel("B".into())).unwrap();
    assert_eq!(choices(&preferred), [vec!["B"]]);
    assert_eq!(
        preferred.report().alternate_locations()[0].reason,
        MmcifAltLocSelectionReason::PreferredLabel
    );
    let fallback = selected(&source, MmcifAltLocPolicy::PreferLabel("Z".into())).unwrap();
    assert_eq!(choices(&fallback), [vec!["A"]]);
    assert_eq!(
        fallback.report().alternate_locations()[0].reason,
        MmcifAltLocSelectionReason::PreferredLabelUnavailable
    );
}

#[test]
fn mmcif_altloc_ties_and_unknown_occupancies_are_reported() {
    let source = document("1 C C1 LIG X 1 1 . B ? 1 0 0 1\n2 C C1 LIG X 1 1 . A . 0 0 0 1");
    let result = source.interpret().unwrap();
    assert_eq!(choices(&result), [vec!["A"]]);
    assert!(result.report().issues().iter().any(|issue| matches!(issue, MmcifInterpretIssue::AlternateLocationOccupancyTie { alternatives, .. } if alternatives == &["A", "B"])));
    assert!(result.report().issues().iter().any(|issue| matches!(issue, MmcifInterpretIssue::AlternateLocationOccupancyMissing { source_lines, .. } if source_lines.len() == 2)));
    assert_eq!(
        result
            .model()
            .occupancy(result.model().topology().atom_ids()[0])
            .unwrap(),
        None
    );
}

#[test]
fn mmcif_altloc_incomplete_alternatives_are_never_patched() {
    let source = document("1 C C1 LIG X 1 1 . A 0.9 0 0 0 1\n2 C C1 LIG X 1 1 . B 0.1 5 0 0 1\n3 O O1 LIG X 1 1 . B 0.1 6.2 0 0 1");
    let result = source.interpret().unwrap();
    assert_eq!(choices(&result), [vec!["B"]]);
    assert_eq!(result.model().atom_count(), 2);
    assert!(result.report().issues().iter().any(|issue| matches!(issue, MmcifInterpretIssue::AlternateLocationIncomplete { alternative, missing_atoms, .. } if alternative == "A" && missing_atoms == &["O1"])));
    assert!(selected(&source, MmcifAltLocPolicy::SelectLabel("A".into())).is_err());
    let disjoint = document("1 C C1 LIG X 1 1 . A 0.9 0 0 0 1\n2 O O1 LIG X 1 1 . B 0.1 5 0 0 1");
    assert!(disjoint.interpret().is_err());
}

fn two_residues() -> MmcifDocument {
    document("1 C C1 LIG X 1 11 . A 0.9 0 0 0 1\n2 C C1 LIG X 1 11 . B 0.1 1 0 0 1\n3 C C1 LIG X 2 12 . A 0.2 5 0 0 1\n4 C C1 LIG X 2 12 . B 0.8 6 0 0 1")
}

#[test]
fn mmcif_altloc_local_overrides_and_linked_groups() {
    let source = two_residues();
    let ids = source.blocks()[0]
        .alternate_locations()
        .unwrap()
        .into_iter()
        .map(|residue| residue.residue)
        .collect::<Vec<_>>();
    assert_eq!(
        choices(&source.interpret().unwrap()),
        [vec!["A"], vec!["B"]]
    );
    let result = selected(
        &source,
        MmcifAltLocPolicy::Configured(MmcifAltLocSelection {
            default: MmcifAltLocPreference::SelectLabel("A".into()),
            overrides: BTreeMap::from([(ids[1].clone(), "B".into())]),
            ..Default::default()
        }),
    )
    .unwrap();
    assert_eq!(choices(&result), [vec!["A"], vec!["B"]]);
    let group = MmcifAltLocSelection {
        groups: vec![ids.clone()],
        ..Default::default()
    };
    assert_eq!(
        choices(&selected(&source, MmcifAltLocPolicy::Configured(group.clone())).unwrap()),
        [vec!["A"], vec!["A"]]
    );
    let forced = MmcifAltLocSelection {
        overrides: BTreeMap::from([(ids[1].clone(), "B".into())]),
        ..group.clone()
    };
    assert_eq!(
        choices(&selected(&source, MmcifAltLocPolicy::Configured(forced)).unwrap()),
        [vec!["B"], vec!["B"]]
    );
    let conflicting = MmcifAltLocSelection {
        overrides: BTreeMap::from([(ids[0].clone(), "A".into()), (ids[1].clone(), "B".into())]),
        ..group
    };
    assert!(selected(&source, MmcifAltLocPolicy::Configured(conflicting)).is_err());
}

#[test]
fn mmcif_altloc_invalid_selectors_and_groups_fail() {
    let source = two_residues();
    let ids = source.blocks()[0]
        .alternate_locations()
        .unwrap()
        .into_iter()
        .map(|residue| residue.residue)
        .collect::<Vec<_>>();
    let mut unknown = ids[0].clone();
    unknown.insertion_code = Some("A".into());
    for selection in [
        MmcifAltLocSelection {
            overrides: BTreeMap::from([(unknown.clone(), "A".into())]),
            ..Default::default()
        },
        MmcifAltLocSelection {
            groups: vec![vec![]],
            ..Default::default()
        },
        MmcifAltLocSelection {
            groups: vec![vec![ids[0].clone()]],
            ..Default::default()
        },
        MmcifAltLocSelection {
            groups: vec![vec![ids[0].clone(), ids[0].clone()]],
            ..Default::default()
        },
        MmcifAltLocSelection {
            groups: vec![vec![ids[0].clone(), unknown]],
            ..Default::default()
        },
    ] {
        assert!(selected(&source, MmcifAltLocPolicy::Configured(selection)).is_err());
    }
}

#[test]
fn mmcif_altloc_residue_identity_preserves_namespaces_and_insertion_codes() {
    let source = document("1 C C1 LIG X 7 17 . A 1 0 0 0 1\n2 C C1 LIG X 7 17 I B 1 3 0 0 1\n3 C C1 LIG X . 17 . A 1 6 0 0 1\n4 C C1 LIG Y . . . A 1 9 0 0 1\n5 C C1 LIG Y . . . A 1 12 0 0 1");
    let inventory = source.blocks()[0].alternate_locations().unwrap();
    assert_eq!(inventory.len(), 5);
    assert_eq!(
        inventory[0].residue.position,
        MmcifResiduePosition::Label(7)
    );
    assert_eq!(inventory[1].residue.insertion_code.as_deref(), Some("I"));
    assert_eq!(
        inventory[2].residue.position,
        MmcifResiduePosition::Author("17".into())
    );
    assert_eq!(
        inventory[3].residue.position,
        MmcifResiduePosition::Occurrence {
            component_id: "LIG".into(),
            index: 0
        }
    );
    assert_eq!(
        inventory[4].residue.position,
        MmcifResiduePosition::Occurrence {
            component_id: "LIG".into(),
            index: 1
        }
    );
}

#[test]
fn mmcif_altloc_microheterogeneity_selects_one_component() {
    let source = document("1 C CA GLY X 1 1 . A 0.4 0 0 0 1\n2 C CA ALA X 1 1 . B 0.6 5 0 0 1\n3 C CB ALA X 1 1 . B 0.6 6.5 0 0 1");
    let result = source.interpret().unwrap();
    assert_eq!(choices(&result), [vec!["B"]]);
    assert_eq!(result.model().atom_count(), 2);
    assert!(result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .all(|atom| atom.component_id() == "ALA"));
    let gly = selected(&source, MmcifAltLocPolicy::SelectLabel("A".into())).unwrap();
    assert_eq!(gly.model().atom_count(), 1);
}

#[test]
fn mmcif_altloc_rejects_single_label_and_duplicate_sites() {
    let source = document("1 C C1 LIG X 1 1 . A 1 0 0 0 1");
    assert!(selected(&source, MmcifAltLocPolicy::ErrorOnAlternateLocations).is_err());
    let duplicate = document("1 C C1 LIG X 1 1 . A 1 0 0 0 1\n2 C C1 LIG X 1 1 . A 1 1 0 0 1");
    assert!(duplicate
        .interpret()
        .unwrap_err()
        .message()
        .contains("duplicate"));
    let overlapping = document("1 C C1 LIG X 1 1 . A 1 0 0 0 1\n2 C C1 LIG X 1 1 . . 1 1 0 0 1");
    assert!(overlapping
        .interpret()
        .unwrap_err()
        .message()
        .contains("overlapping"));
}

fn with_conformations(rows: &str, mapping: &str) -> MmcifDocument {
    let source = document(rows);
    // Build the same minimal syntax with an additional independent category.
    let atom_table = source.blocks()[0].loop_with_tag("_atom_site.id").unwrap();
    let mut text = String::from("data_alternates\nloop_\n");
    for tag in atom_table.tags() {
        text.push_str(tag);
        text.push('\n');
    }
    text.push_str(rows);
    text.push_str("\nloop_\n_atom_sites_alt_gen.ens_id\n_atom_sites_alt_gen.alt_id\n");
    text.push_str(mapping);
    parse_str(&text).unwrap()
}

const DECLARED: &str = "1 C C1 LIG X 1 1 . A 0.9 0 0 0 1\n2 C C1 LIG X 1 1 . B 0.1 1 0 0 1\n3 C C1 LIG X 2 2 . C 0.2 5 0 0 1\n4 C C1 LIG X 2 2 . D 0.8 6 0 0 1";

#[test]
fn mmcif_altloc_source_combinations_constrain_default_and_exact_selection() {
    let source = with_conformations(DECLARED, "first A\nfirst C\nfirst .\nsecond B\nsecond D\n");
    let original = source.clone();
    let result = source.interpret().unwrap();
    assert_eq!(choices(&result), [vec!["A"], vec!["C"]]);
    assert_eq!(
        result.report().alternate_locations()[0].reason,
        MmcifAltLocSelectionReason::SourceConformation { id: "first".into() }
    );
    let second = selected(
        &source,
        MmcifAltLocPolicy::Configured(MmcifAltLocSelection {
            source_conformation: Some("second".into()),
            ..Default::default()
        }),
    )
    .unwrap();
    assert_eq!(choices(&second), [vec!["B"], vec!["D"]]);
    assert_eq!(source, original);
    let ids = source.blocks()[0].alternate_locations().unwrap();
    let conflict = MmcifAltLocSelection {
        overrides: BTreeMap::from([
            (ids[0].residue.clone(), "A".into()),
            (ids[1].residue.clone(), "D".into()),
        ]),
        ..Default::default()
    };
    assert!(selected(&source, MmcifAltLocPolicy::Configured(conflict)).is_err());
    assert_eq!(source.interpret_ensemble().unwrap().ensemble().len(), 1);
}

#[test]
fn mmcif_altloc_source_combinations_can_cover_disjoint_parts_of_a_residue() {
    let source = with_conformations("1 C C1 LIG X 1 1 . A 0.7 0 0 0 1\n2 C C1 LIG X 1 1 . B 0.3 5 0 0 1\n3 O O1 LIG X 1 1 . C 0.7 1.2 0 0 1\n4 O O1 LIG X 1 1 . D 0.3 6.2 0 0 1", "first A\nfirst C\nsecond B\nsecond D\n");
    let result = source.interpret().unwrap();
    assert_eq!(choices(&result), [vec!["A", "C"]]);
    assert_eq!(result.model().atom_count(), 2);
}

#[test]
fn mmcif_altloc_invalid_source_combinations_fail() {
    for mapping in [
        "one Z\n",
        "one ?\n",
        "one A\none A\n",
        "one A\n",
        "one A\none B\none C\n",
    ] {
        assert!(
            with_conformations(DECLARED, mapping).interpret().is_err(),
            "{mapping}"
        );
    }
    assert!(selected(
        &document(CROSSED),
        MmcifAltLocPolicy::Configured(MmcifAltLocSelection {
            source_conformation: Some("absent".into()),
            ..Default::default()
        })
    )
    .is_err());
}

#[test]
fn mmcif_altloc_ensemble_uses_source_models_and_stable_site_order() {
    let rows = "1 C C1 LIG X 1 1 . A 0.9 0 0 0 1\n2 O O1 LIG X 1 1 . A 0.9 1.2 0 0 1\n3 O O1 LIG X 1 1 . B 0.1 6.2 0 0 1\n4 C C1 LIG X 1 1 . B 0.1 5 0 0 1\n5 C C1 LIG X 1 1 . A 0.1 10 0 0 2\n6 O O1 LIG X 1 1 . A 0.1 11.2 0 0 2\n7 O O1 LIG X 1 1 . B 0.9 16.2 0 0 2\n8 C C1 LIG X 1 1 . B 0.9 15 0 0 2";
    let source = document(rows);
    let result = source.interpret_ensemble().unwrap();
    assert_eq!(result.ensemble().len(), 2);
    assert_eq!(
        result.reports()[0].alternate_locations()[0].selected_labels,
        ["A"]
    );
    assert_eq!(
        result.reports()[1].alternate_locations()[0].selected_labels,
        ["B"]
    );
    for (index, report) in result.reports().iter().enumerate() {
        assert_eq!(report.alternate_locations().len(), 1);
        assert_eq!(
            report.alternate_locations()[0].residue.model_id,
            (index + 1).to_string()
        );
        assert!(report
            .issues()
            .iter()
            .filter_map(|issue| match issue {
                MmcifInterpretIssue::AlternateLocationOmitted { residue, .. } => Some(residue),
                _ => None,
            })
            .all(|residue| residue.model_id == (index + 1).to_string()));
    }
    let ids = source.blocks()[0]
        .alternate_locations()
        .unwrap()
        .into_iter()
        .map(|residue| residue.residue)
        .collect::<Vec<_>>();
    let result = source
        .interpret_ensemble_with_options(MmcifEnsembleInterpretOptions {
            altloc_policy: MmcifAltLocPolicy::Configured(MmcifAltLocSelection {
                overrides: BTreeMap::from([(ids[1].clone(), "A".into())]),
                ..Default::default()
            }),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        result.reports()[1].alternate_locations()[0].selected_labels,
        ["A"]
    );
    assert!(source
        .interpret_ensemble_with_options(MmcifEnsembleInterpretOptions {
            altloc_policy: MmcifAltLocPolicy::Configured(MmcifAltLocSelection {
                groups: vec![ids],
                ..Default::default()
            }),
            ..Default::default()
        })
        .is_err());
}

#[test]
fn mmcif_altloc_ensemble_rejects_changed_chemistry_and_keeps_single_plain_model() {
    let source = document("1 C CA GLY X 1 1 . A 0.9 0 0 0 1\n2 C CA ALA X 1 1 . B 0.1 5 0 0 1\n3 C CB ALA X 1 1 . B 0.1 6.5 0 0 1\n4 C CA GLY X 1 1 . A 0.1 0 0 0 2\n5 C CA ALA X 1 1 . B 0.9 5 0 0 2\n6 C CB ALA X 1 1 . B 0.9 6.5 0 0 2");
    assert!(matches!(
        source.interpret_ensemble(),
        Err(MmcifEnsembleInterpretError::InconsistentAtomSet { .. })
    ));
    let plain = document("1 C C1 LIG X 1 1 . . 1 0 0 0 1");
    assert!(plain.blocks()[0].alternate_locations().unwrap().is_empty());
    assert_eq!(plain.interpret_ensemble().unwrap().ensemble().len(), 1);
}

#[test]
fn mmcif_altloc_explicit_conformer_ensemble_is_unweighted_and_preserves_order() {
    let source = document(CROSSED);
    let options = ["B", "A", "B"].map(|label| MmcifInterpretOptions {
        altloc_policy: MmcifAltLocPolicy::SelectLabel(label.into()),
        ..Default::default()
    });
    let result = interpret_conformations(&source, &options).unwrap();
    assert_eq!(result.ensemble().len(), 3);
    assert!(result
        .ensemble()
        .members()
        .all(|member| member.weight().is_none()));
    assert_eq!(
        result
            .reports()
            .iter()
            .map(|report| report.alternate_locations()[0].selected_labels[0].as_str())
            .collect::<Vec<_>>(),
        ["B", "A", "B"]
    );
    assert!(result
        .reports()
        .iter()
        .all(|report| report.selected_model() == Some("1")));
    let block = interpret_conformations_block(&source.blocks()[0], &options).unwrap();
    assert_eq!(block.reports(), result.reports());
    assert!(block
        .ensemble()
        .topology()
        .same_layout(result.ensemble().topology()));
    assert!(matches!(
        source.interpret_conformations(&[]),
        Err(MmcifEnsembleInterpretError::EmptyConformationSelection)
    ));
    let bad = [MmcifInterpretOptions {
        altloc_policy: MmcifAltLocPolicy::SelectLabel("Z".into()),
        ..Default::default()
    }];
    assert!(matches!(
        source.interpret_conformations(&bad),
        Err(MmcifEnsembleInterpretError::Conformation { selection: 0, .. })
    ));
}

#[test]
fn mmcif_altloc_explicit_conformer_ensemble_rejects_microheterogeneity() {
    let source = document("1 C CA GLY X 1 1 . A 0.4 0 0 0 1\n2 C CA ALA X 1 1 . B 0.6 5 0 0 1\n3 C CB ALA X 1 1 . B 0.6 6.5 0 0 1");
    let options = ["A", "B"].map(|label| MmcifInterpretOptions {
        altloc_policy: MmcifAltLocPolicy::SelectLabel(label.into()),
        ..Default::default()
    });
    assert!(matches!(
        source.interpret_conformations(&options),
        Err(MmcifEnsembleInterpretError::InconsistentAtomSet { .. })
    ));
}

#[test]
fn mmcif_altloc_unselected_models_still_validate_alternate_sites() {
    let source = document("1 C C1 LIG X 1 1 . . 1 0 0 0 1\n2 C C1 LIG X 1 1 . A 1 0 0 0 2\n3 C C1 LIG X 1 1 . A 1 1 0 0 2");
    let error = source
        .interpret_with_options(MmcifInterpretOptions {
            model_selection: MmcifModelSelection::Select("1".into()),
            ..Default::default()
        })
        .unwrap_err();
    assert!(error.message().contains("duplicate"));
}

#[test]
fn mmcif_altloc_mean_occupancy_does_not_favour_larger_components() {
    let source = document("1 C CA GLY X 1 1 . A 0.6 0 0 0 1\n2 C CA ALA X 1 1 . B 0.4 5 0 0 1\n3 C CB ALA X 1 1 . B 0.4 6.5 0 0 1");
    assert_eq!(choices(&source.interpret().unwrap()), [vec!["A"]]);
}

#[test]
fn mmcif_altloc_source_combinations_respect_user_group_constraints() {
    let rows = "1 C C1 LIG X 1 1 . A 0.9 0 0 0 1\n2 C C1 LIG X 1 1 . B 0.1 1 0 0 1\n3 C C1 LIG X 2 2 . C 0.2 5 0 0 1\n4 C C1 LIG X 2 2 . D 0.8 6 0 0 1";
    let source = with_conformations(rows, "first A\nfirst C\nsecond B\nsecond D\n");
    let ids = source.blocks()[0]
        .alternate_locations()
        .unwrap()
        .into_iter()
        .map(|residue| residue.residue)
        .collect();
    assert!(selected(
        &source,
        MmcifAltLocPolicy::Configured(MmcifAltLocSelection {
            groups: vec![ids],
            ..Default::default()
        })
    )
    .is_err());
}

#[test]
fn mmcif_altloc_overlapping_groups_are_transitive() {
    let source = document("1 C C1 LIG X 1 1 . A 0.9 0 0 0 1\n2 C C1 LIG X 1 1 . B 0.1 1 0 0 1\n3 C C1 LIG X 2 2 . A 0.2 5 0 0 1\n4 C C1 LIG X 2 2 . B 0.8 6 0 0 1\n5 C C1 LIG X 3 3 . A 0.1 10 0 0 1\n6 C C1 LIG X 3 3 . B 0.9 11 0 0 1");
    let ids = source.blocks()[0]
        .alternate_locations()
        .unwrap()
        .into_iter()
        .map(|residue| residue.residue)
        .collect::<Vec<_>>();
    let result = selected(
        &source,
        MmcifAltLocPolicy::Configured(MmcifAltLocSelection {
            groups: vec![ids[..2].to_vec(), ids[1..].to_vec()],
            ..Default::default()
        }),
    )
    .unwrap();
    assert_eq!(choices(&result), [vec!["B"], vec!["B"], vec!["B"]]);
}

#[test]
fn mmcif_altloc_rejected_source_combinations_are_reported() {
    let source = with_conformations(DECLARED, "bad A\nbad B\nbad C\ngood A\ngood C\n");
    let result = source.interpret().unwrap();
    assert_eq!(choices(&result), [vec!["A"], vec!["C"]]);
    assert!(result.report().issues().iter().any(|issue| matches!(issue,
        MmcifInterpretIssue::AlternateLocationCombinationRejected { alternative, reason, .. }
        if alternative == "bad" && reason.contains("overlapping")
    )));
    let source = with_conformations(DECLARED, "bad A\ngood A\ngood C\n");
    assert!(source.interpret().unwrap().report().issues().iter().any(|issue| matches!(issue,
        MmcifInterpretIssue::AlternateLocationCombinationRejected { alternative, .. } if alternative == "bad"
    )));
}

#[test]
fn mmcif_altloc_selection_does_not_change_first_source_model() {
    let source = document("1 C CA GLY X 1 1 . A 0.1 0 0 0 7\n2 C CA GLY X 1 1 . . 1 0 0 0 2\n3 C CA ALA X 1 1 . B 0.9 0 0 0 7");
    let result = source
        .interpret_with_options(MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(result.report().selected_model(), Some("7"));
    let other = source
        .interpret_with_options(MmcifInterpretOptions {
            model_selection: MmcifModelSelection::Select("2".into()),
            ..Default::default()
        })
        .unwrap();
    assert!(other.report().issues().iter().any(|issue| matches!(issue, MmcifInterpretIssue::CoordinateModelIgnored { model_id, atom_site_rows: 2 } if model_id == "7")));
}
