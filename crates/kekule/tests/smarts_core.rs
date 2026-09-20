use kekule::{
    core::*, hydrogens, perception::aromaticity, query::*, smiles, substructure::*,
    topology::Topology,
};
use std::ops::ControlFlow;
use std::sync::Arc;

fn molecule(source: &str) -> Molecule {
    let mut m = smiles::to_molecules(source).unwrap().pop().unwrap();
    m.perceive().unwrap();
    m
}
fn options() -> SubstructureMatchOptions {
    SubstructureMatchOptions {
        max_matches: 100_000,
        uniquify: false,
        ..Default::default()
    }
}
fn count(source: &str, pattern: &str) -> usize {
    find_substructure_matches_complete(
        &molecule(source),
        &parse_smarts(pattern).unwrap(),
        options(),
    )
    .unwrap()
    .len()
}

#[test]
fn recursive_queries_are_anchored_nested_and_independent_of_outer_embedding() {
    for (s, q, n) in [
        ("CC(=O)O", "[$(C=O)]", 1),
        ("CC(=O)O", "[!$(C=O)]", 3),
        ("CCO", "[$(C[$(CO)])]", 1),
        ("CCO", "[$(CO),$(OC)]", 2),
        ("CC", "C[$(CC)]", 2),
        ("CC", "[$(CC);$(CC)]", 2),
        ("CO", "[$(O)]C", 1),
        ("CCC", "[$(C.C)]", 3),
    ] {
        assert_eq!(count(s, q), n, "{s} {q}");
    }
}

#[test]
fn bond_boolean_operators_and_precedence() {
    for (s, q, n) in [
        ("CC=O", "C-,=O", 1),
        ("CC=O", "C!#O", 1),
        ("C1CC1C", "C-;!@C", 2),
        ("C1CC1C", "C-,=;@C", 6),
        ("CC=O", "C!-,=O", 1),
        ("c1ccccc1", "c!:c", 0),
        ("C1CC1", "C1-;@C-;@C-;@1", 6),
    ] {
        assert_eq!(count(s, q), n, "{s} {q}");
    }
}

#[test]
fn hydrogen_valence_and_ring_primitives_are_distinct() {
    for (s, q, n) in [
        ("[CH4]", "[h4]", 1),
        ("[nH]1cccc1", "[h1]", 5),
        ("CCO", "[h]", 3),
        ("CCO", "[v4]", 2),
        ("CCO", "[v2]", 1),
        ("C1CC1C", "[r3]", 3),
        ("C1CC1C", "[R0]", 1),
        ("C1CCC2CCCCC2C1", "[R2]", 2),
        ("C1CCC2CCCCC2C1", "[r6]", 10),
        ("C12C3C4C1C5C2C3C45", "[R3;r4]", 8),
    ] {
        assert_eq!(count(s, q), n, "{s} {q}");
    }
    let mut m = molecule("C");
    hydrogens::add_hydrogens(&mut m).unwrap();
    m.perceive().unwrap();
    for (q, n) in [("[h0]", 5), ("[H4]", 1), ("[H:1]", 4), ("[v4]", 1)] {
        assert_eq!(
            find_substructure_matches_complete(&m, &parse_smarts(q).unwrap(), options())
                .unwrap()
                .len(),
            n,
            "{q}"
        );
    }
}

#[test]
fn tags_are_outputs_and_do_not_match_target_map_values() {
    let m = molecule("[CH3:99][OH:88]");
    let q = parse_smarts("[#6:20]-[#8:3]").unwrap();
    let tagged = TaggedQuery::new(&q).unwrap();
    assert_eq!(
        tagged.tags().iter().map(|(t, _)| *t).collect::<Vec<_>>(),
        vec![3, 20]
    );
    let result = tagged
        .find_matches(&PreparedTarget::new(&m), options(), true)
        .unwrap();
    assert_eq!(result, vec![vec![AtomId::new(1), AtomId::new(0)]]);
    for (source, error) in [
        ("[C:0]", TaggedMatchError::ZeroTag),
        ("[C:1][C:1]", TaggedMatchError::DuplicateTag(1)),
        ("CC", TaggedMatchError::MissingTags),
    ] {
        assert_eq!(
            TaggedQuery::new(&parse_smarts(source).unwrap()).unwrap_err(),
            error
        );
    }
    let q = parse_smarts("[$([C:9]):1]").unwrap();
    assert_eq!(TaggedQuery::new(&q).unwrap().tags().len(), 1);
}

#[test]
fn tagged_automorphisms_preserve_roles_but_can_deduplicate_context() {
    let m = molecule("CCC");
    let p = PreparedTarget::new(&m);
    let q = parse_smarts("[C:1][C:2]").unwrap();
    assert_eq!(
        TaggedQuery::new(&q)
            .unwrap()
            .find_matches(&p, options(), true)
            .unwrap()
            .len(),
        4
    );
    let q = parse_smarts("[C:1](C)C").unwrap();
    let tagged = TaggedQuery::new(&q).unwrap();
    assert_eq!(tagged.find_matches(&p, options(), false).unwrap().len(), 2);
    assert_eq!(tagged.find_matches(&p, options(), true).unwrap().len(), 1);
}

#[test]
fn complete_enumeration_never_silently_truncates() {
    let m = molecule(&"C".repeat(34));
    let q = parse_smarts("[#6].[#6]").unwrap();
    assert_eq!(
        find_substructure_matches_complete(&m, &q, options())
            .unwrap()
            .len(),
        34 * 33
    );
    let capped = SubstructureMatchOptions {
        max_matches: 1000,
        ..options()
    };
    assert_eq!(
        find_substructure_matches_with_options(&m, &q, capped)
            .unwrap()
            .len(),
        1000
    );
    assert!(matches!(
        find_substructure_matches_complete(&m, &q, capped),
        Err(SubstructureMatchError::ResourceLimit {
            resource: "matches",
            ..
        })
    ));
    let mut visits = 0;
    assert_eq!(
        visit_substructure_matches(&m, &q, capped, |_| {
            visits += 1;
            ControlFlow::Break(())
        })
        .unwrap(),
        MatchCompletion::Stopped
    );
    assert_eq!(visits, 1);
    let small = molecule("CC");
    assert_eq!(
        find_substructure_matches_complete(
            &small,
            &q,
            SubstructureMatchOptions {
                max_matches: 2,
                ..options()
            }
        )
        .unwrap()
        .len(),
        2
    );
}

#[test]
fn recursive_limits_errors_and_prerequisites_survive_negation() {
    let q = parse_smarts("[!$([H3])]").unwrap();
    let unperceived = smiles::to_molecules("CC").unwrap().pop().unwrap();
    assert!(matches!(
        find_substructure_match(&unperceived, &q),
        Err(SubstructureMatchError::MissingPerception(
            QueryPerception::Valence
        ))
    ));
    assert!(matches!(
        find_substructure_matches_complete(
            &molecule("CCCC"),
            &q,
            SubstructureMatchOptions {
                max_search_states: 1,
                ..options()
            }
        ),
        Err(SubstructureMatchError::ResourceLimit { .. })
    ));
    for (source, opts) in [
        (
            "[$(CC)]",
            SmartsParseOptions {
                max_atoms: 2,
                ..Default::default()
            },
        ),
        (
            "[$(C-C)]-C",
            SmartsParseOptions {
                max_bonds: 1,
                ..Default::default()
            },
        ),
        (
            "[$([$([C])])]",
            SmartsParseOptions {
                max_recursive_depth: 1,
                ..Default::default()
            },
        ),
        (
            "[$([C,N])]",
            SmartsParseOptions {
                max_total_expression_nodes: 2,
                ..Default::default()
            },
        ),
    ] {
        assert_eq!(
            parse_smarts_with_options(source, opts).unwrap_err().kind(),
            SmartsParseErrorKind::ResourceLimit
        );
    }
    for source in [
        "[$()]", "[$([C)]", "[$(C)", "[C:]", "[C:1:2]", "C-,C", "C;&C",
    ] {
        let e = parse_smarts(source).unwrap_err();
        assert!(
            e.span().start <= e.span().end && e.span().end <= source.len(),
            "{source}: {e}"
        );
    }
}

#[test]
fn topology_disconnected_queries_keep_occurrences_and_snapshot() {
    let molecules = vec![molecule("CC"), molecule("C")];
    let topology = Arc::new(Topology::from_molecules(&molecules).unwrap());
    let prepared = PreparedTopologyTarget::new(&topology);
    let q = parse_smarts("[#6:1].[#6:2]").unwrap();
    let matches = prepared.find_matches_complete(&q, options()).unwrap();
    assert_eq!(matches.len(), 6);
    assert!(matches.iter().all(|m| Arc::ptr_eq(m.topology(), &topology)));
    assert_eq!(
        matches
            .iter()
            .filter(|m| m.atoms()[0].molecule() != m.atoms()[1].molecule())
            .count(),
        4
    );
    let q = parse_smarts("[#6]-[#6]").unwrap();
    assert_eq!(
        prepared.find_matches_complete(&q, options()).unwrap().len(),
        2
    );
}

#[test]
fn boolean_stereo_keeps_atom_expression_context() {
    assert_eq!(count("N[C@](F)(Cl)Br", "N[C@,C@@](F)(Cl)Br"), 1);
    assert_eq!(count("N[C@@](F)(Cl)Br", "N[C@,C@@](F)(Cl)Br"), 1);
    assert_eq!(count("NC(F)(Cl)Br", "N[C@,C@@](F)(Cl)Br"), 0);
    assert_eq!(count("N[C@@](F)(Cl)Br", "N[C;!@](F)(Cl)Br"), 1);
    assert_eq!(count("N[C@](F)(Cl)Br", "N[C;!@](F)(Cl)Br"), 0);
    assert_eq!(count("N[C@@](F)(Cl)Br", "[$(N[C@@](F)(Cl)Br)]"), 1);
}

#[test]
fn mdl_aromaticity_is_explicit_and_does_not_rewrite_chemistry() {
    for (source, expected) in [
        ("c1ccccc1", 6),
        ("n1ccccc1", 6),
        ("c1ccoc1", 0),
        ("c1cc[nH]c1", 0),
        ("c1ccc2ccccc2c1", 10),
        ("O=C1C=CC(=O)C=C1", 0),
        ("C1=CC=C2C=CC=C2C=C1", 10),
    ] {
        let mut m = molecule(source);
        let graph = m.graph().clone();
        aromaticity::perceive_aromaticity(&mut m, AromaticityModel::Mdl).unwrap();
        assert_eq!(m.graph(), &graph);
        assert_eq!(
            m.perception().aromaticity_model(),
            Some(AromaticityModel::Mdl)
        );
        assert_eq!(
            find_substructure_matches_complete(&m, &parse_smarts("[a]").unwrap(), options())
                .unwrap()
                .len(),
            expected,
            "{source}"
        );
    }
}

#[test]
fn incompatible_ring_models_fail_and_mdl_failure_is_transactional() {
    for model in [Some(RingBasisModel::DepthFirstFallback), None] {
        let mut m = molecule("C1CC1");
        let state = Perception::builder()
            .with_rings(
                m.ring_membership().unwrap().clone(),
                Some(RingBasisState::new(model, m.ring_set().unwrap().clone())),
            )
            .build();
        m.install_perception(state).unwrap();
        let before = m.perception().clone();
        assert!(matches!(
            find_substructure_match(&m, &parse_smarts("[R1]").unwrap()),
            Err(SubstructureMatchError::IncompatiblePerception(
                QueryPerception::RingBasis
            ))
        ));
        assert!(aromaticity::perceive_aromaticity(&mut m, AromaticityModel::Mdl).is_err());
        assert_eq!(m.perception(), &before);
        assert_eq!(
            find_substructure_matches_complete(&m, &parse_smarts("[R]").unwrap(), options())
                .unwrap()
                .len(),
            3
        );
    }
}

#[test]
fn topology_reuses_connected_matches_and_shares_budgets() {
    let mut builder = kekule::topology::TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule("CC")).unwrap();
    for _ in 0..100 {
        builder.add_instance(definition).unwrap();
    }
    let topology = Arc::new(builder.build().unwrap());
    let query = parse_smarts("[#6:1]-[#6:2]").unwrap();
    let small_budget = SubstructureMatchOptions {
        max_candidate_pairs: 4,
        max_search_states: 20,
        ..options()
    };
    assert_eq!(
        find_topology_substructure_matches_complete(&topology, &query, small_budget)
            .unwrap()
            .len(),
        200
    );
    assert!(matches!(
        find_topology_substructure_matches_complete(
            &topology,
            &query,
            SubstructureMatchOptions {
                max_matches: 199,
                ..small_budget
            }
        ),
        Err(SubstructureMatchError::ResourceLimit {
            resource: "matches",
            ..
        })
    ));
}

#[test]
fn recursive_disconnected_environment_can_span_topology_instances() {
    let topology = Arc::new(Topology::from_molecules(&[molecule("C"), molecule("C")]).unwrap());
    let query = parse_smarts("[$(C.C)]").unwrap();
    assert_eq!(
        find_topology_substructure_matches_complete(&topology, &query, options())
            .unwrap()
            .len(),
        2
    );
    assert_eq!(count("C", "[$(C.C)]"), 0);
}

#[test]
fn directional_alternatives_keep_boolean_context() {
    for target in ["F/C=C/Cl", "F/C=C\\Cl", "FC=CCl"] {
        assert_eq!(count(target, "F-,/C=C/Cl"), 1);
    }
    assert_eq!(count("F/C=C/Cl", "F/,\\C=C/Cl"), 1);
    assert_eq!(count("F/C=C\\Cl", "F/,\\C=C/Cl"), 1);
    assert_eq!(count("FC=CCl", "F/,\\C=C/Cl"), 0);
}

#[test]
fn prepared_results_follow_identity_after_source_atom_reordering() {
    let q = parse_smarts("[#6:1]-[#8:2]").unwrap();
    let tagged = TaggedQuery::new(&q).unwrap();
    for source in ["[CH3:4][CH2:5][OH:6]", "[OH:6][CH2:5][CH3:4]"] {
        let m = molecule(source);
        let prepared = PreparedTarget::new(&m);
        let matches = tagged.find_matches(&prepared, options(), true).unwrap();
        assert_eq!(matches.len(), 1);
        assert_eq!(
            matches[0]
                .iter()
                .map(|id| m.atom(*id).unwrap().atom_map)
                .collect::<Vec<_>>(),
            vec![Some(5), Some(6)]
        );
        assert_eq!(
            tagged.find_matches(&prepared, options(), true).unwrap(),
            matches
        );
    }
    assert_eq!(count("CC", "[0C]"), 2);
}

#[test]
fn graph_complexity_exhaustion_remains_a_resource_error() {
    let input = "*".repeat(4097);
    let error = parse_smarts_with_options(
        &input,
        SmartsParseOptions {
            max_atoms: 5000,
            max_bonds: 5000,
            max_total_expression_nodes: 65536,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert_eq!(error.kind(), SmartsParseErrorKind::ResourceLimit);
    assert_eq!(error.span(), 0..input.len());
}
