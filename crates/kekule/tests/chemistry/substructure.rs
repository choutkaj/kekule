use kekule::diagnostics::{DiagnosticKind, LimitExceeded, Resource, SourceLocation};
use kekule::{
    molecule::*,
    perception::{
        aromaticity::{self, AromaticityModel},
        rings::{RingBasisModel, RingBasisState},
        Perception,
    },
    smarts::{self, SmartsParseOptions},
    smiles,
    substructure::*,
    topology::Topology,
};
use std::ops::ControlFlow;
use std::sync::Arc;

#[test]
fn connected_query_extensions_scale_with_local_adjacency() {
    // Visiting length-three windows must not retry every remote atom at each
    // extension. This bounded test used to exhaust its search-state budget.
    let m = smiles::to_molecules(&"C".repeat(1500)).unwrap().remove(0);
    let q = smarts::parse_str("[#6:1]~[#6:2]~[#6:3]", Default::default()).unwrap();
    let matches = PreparedTarget::new(&m)
        .all(
            &q,
            SubstructureMatchOptions::default()
                .with_max_matches(10_000)
                .unwrap()
                .with_max_search_states(30_000)
                .unwrap()
                .with_uniquify(false),
        )
        .unwrap();
    assert_eq!(matches.len(), 2 * (1500 - 2));
    let ids = m.atom_ids().collect::<Vec<_>>();
    let actual = matches
        .iter()
        .map(|m| m.atoms().to_vec())
        .collect::<std::collections::BTreeSet<_>>();
    let expected = ids
        .windows(3)
        .flat_map(|a| [a.to_vec(), a.iter().rev().copied().collect()])
        .collect();
    assert_eq!(actual, expected);
}

fn molecule(source: &str) -> Molecule {
    let mut m = smiles::to_molecules(source).unwrap().pop().unwrap();
    m.perceive().unwrap();
    m
}
fn options() -> SubstructureMatchOptions {
    SubstructureMatchOptions::default()
        .with_max_matches(100_000)
        .unwrap()
        .with_uniquify(false)
}
fn count(source: &str, pattern: &str) -> usize {
    PreparedTarget::new(&molecule(source))
        .all(
            &smarts::parse_str(pattern, Default::default()).unwrap(),
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
    m.add_hydrogens(Default::default()).unwrap();
    m.perceive().unwrap();
    for (q, n) in [("[h0]", 5), ("[H4]", 1), ("[H:1]", 4), ("[v4]", 1)] {
        assert_eq!(
            PreparedTarget::new(&m)
                .all(
                    &smarts::parse_str(q, Default::default()).unwrap(),
                    options()
                )
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
    let q = smarts::parse_str("[#6:20]-[#8:3]", Default::default()).unwrap();
    let tagged = TaggedQuery::new(&q).unwrap();
    assert_eq!(
        tagged.tags().iter().map(|(t, _)| *t).collect::<Vec<_>>(),
        vec![3, 20]
    );
    let result = tagged
        .find_matches(&PreparedTarget::new(&m), options(), TupleMode::Distinct)
        .unwrap();
    assert_eq!(result, vec![vec![AtomId::new(1), AtomId::new(0)]]);
    for (source, error) in [
        ("[C:0]", TaggedMatchError::ZeroTag),
        ("[C:1][C:1]", TaggedMatchError::DuplicateTag(1)),
        ("CC", TaggedMatchError::MissingTags),
    ] {
        assert_eq!(
            TaggedQuery::new(&smarts::parse_str(source, Default::default()).unwrap()).unwrap_err(),
            error
        );
    }
    let q = smarts::parse_str("[$([C:9]):1]", Default::default()).unwrap();
    assert_eq!(TaggedQuery::new(&q).unwrap().tags().len(), 1);
}

#[test]
fn tagged_automorphisms_preserve_roles_but_can_deduplicate_context() {
    let m = molecule("CCC");
    let p = PreparedTarget::new(&m);
    let q = smarts::parse_str("[C:1][C:2]", Default::default()).unwrap();
    assert_eq!(
        TaggedQuery::new(&q)
            .unwrap()
            .find_matches(&p, options(), TupleMode::Distinct)
            .unwrap()
            .len(),
        4
    );
    let q = smarts::parse_str("[C:1](C)C", Default::default()).unwrap();
    let tagged = TaggedQuery::new(&q).unwrap();
    assert_eq!(
        tagged
            .find_matches(&p, options(), TupleMode::All)
            .unwrap()
            .len(),
        2
    );
    assert_eq!(
        tagged
            .find_matches(&p, options(), TupleMode::Distinct)
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn complete_enumeration_never_silently_truncates() {
    let m = molecule(&"C".repeat(34));
    let q = smarts::parse_str("[#6].[#6]", Default::default()).unwrap();
    assert_eq!(
        PreparedTarget::new(&m).all(&q, options()).unwrap().len(),
        34 * 33
    );
    let capped = options().with_max_matches(1000).unwrap();
    // A match limit is a resource bound, never a silent truncation; callers
    // that want a prefix stop the visitor explicitly.
    let mut prefix = 0;
    assert_eq!(
        PreparedTarget::new(&m)
            .visit(&q, options(), |_| {
                prefix += 1;
                if prefix == 1000 {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(())
                }
            })
            .unwrap(),
        MatchCompletion::Stopped
    );
    assert_eq!(prefix, 1000);
    assert!(matches!(
        PreparedTarget::new(&m).all(&q, capped),
        Err(SubstructureMatchError::LimitExceeded { limit, .. }) if limit.resource() == Resource::Matches
    ));
    let mut visits = 0;
    assert_eq!(
        PreparedTarget::new(&m)
            .visit(&q, capped, |_| {
                visits += 1;
                ControlFlow::Break(())
            })
            .unwrap(),
        MatchCompletion::Stopped
    );
    assert_eq!(visits, 1);
    let small = molecule("CC");
    assert_eq!(
        PreparedTarget::new(&small)
            .all(&q, options().with_max_matches(2).unwrap())
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn recursive_limits_errors_and_prerequisites_survive_negation() {
    let q = smarts::parse_str("[!$([H3])]", Default::default()).unwrap();
    let unperceived = smiles::to_molecules("CC").unwrap().pop().unwrap();
    assert!(matches!(
        PreparedTarget::new(&unperceived).first(&q, Default::default()),
        Err(SubstructureMatchError::MissingPerception(
            QueryPerception::Valence
        ))
    ));
    assert!(matches!(
        PreparedTarget::new(&molecule("CCCC"))
            .all(&q, options().with_max_search_states(1).unwrap()),
        Err(SubstructureMatchError::LimitExceeded { .. })
    ));
    for (source, opts) in [
        (
            "[$(CC)]",
            SmartsParseOptions::default().with_max_atoms(2).unwrap(),
        ),
        (
            "[$(C-C)]-C",
            SmartsParseOptions::default().with_max_bonds(1).unwrap(),
        ),
        (
            "[$([$([C])])]",
            SmartsParseOptions::default()
                .with_max_recursive_depth(1)
                .unwrap(),
        ),
        (
            "[$([C,N])]",
            SmartsParseOptions::default()
                .with_max_total_expression_nodes(2)
                .unwrap(),
        ),
    ] {
        assert_eq!(
            smarts::parse_str(source, opts).unwrap_err().kind(),
            DiagnosticKind::Limit
        );
    }
    for source in [
        "[$()]", "[$([C)]", "[$(C)", "[C:]", "[C:1:2]", "C-,C", "C;&C",
    ] {
        let e = smarts::parse_str(source, Default::default()).unwrap_err();
        let location = e.location().expect("SMARTS errors are located");
        assert_eq!(location.line(), 1, "{source}: {e}");
        assert!(
            (1..=source.len() + 1).contains(&location.column().unwrap()),
            "{source}: {e}"
        );
    }
}

#[test]
fn topology_disconnected_queries_keep_occurrences_and_snapshot() {
    let molecules = vec![molecule("CC"), molecule("C")];
    let topology = Arc::new(Topology::from_molecules(molecules.clone()).unwrap());
    let prepared = PreparedTopologyTarget::new(&topology);
    let q = smarts::parse_str("[#6:1].[#6:2]", Default::default()).unwrap();
    let matches = prepared.all(&q, options()).unwrap();
    assert_eq!(matches.len(), 6);
    assert!(matches
        .iter()
        .all(|m| Arc::ptr_eq(&m.shared_topology(), &topology)));
    assert_eq!(
        matches
            .iter()
            .filter(|m| m.atoms()[0].instance() != m.atoms()[1].instance())
            .count(),
        4
    );
    let q = smarts::parse_str("[#6]-[#6]", Default::default()).unwrap();
    assert_eq!(prepared.all(&q, options()).unwrap().len(), 2);
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
        aromaticity::perceive_aromaticity(&mut m, AromaticityModel::Mdl, Default::default())
            .unwrap();
        assert_eq!(m.graph(), &graph);
        assert_eq!(
            m.perception()
                .aromaticity()
                .map(|aromaticity| aromaticity.model()),
            Some(AromaticityModel::Mdl)
        );
        assert_eq!(
            PreparedTarget::new(&m)
                .all(
                    &smarts::parse_str("[a]", Default::default()).unwrap(),
                    options()
                )
                .unwrap()
                .len(),
            expected,
            "{source}"
        );
    }
}

#[test]
fn incompatible_ring_models_fail_and_mdl_failure_is_transactional() {
    {
        let mut m = molecule("C1CC1");
        let rings = m.perception().rings().unwrap();
        let state = Perception::builder()
            .with_rings(
                rings.membership().clone(),
                Some(RingBasisState::with_depth_first_fallback(
                    RingBasisModel::FiguerasSssrLike,
                    rings.ring_set().unwrap().clone(),
                )),
            )
            .build();
        m.install_perception(state).unwrap();
        let before = m.perception().clone();
        assert!(matches!(
            PreparedTarget::new(&m).first(
                &smarts::parse_str("[R1]", Default::default()).unwrap(),
                Default::default()
            ),
            Err(SubstructureMatchError::IncompatiblePerception(
                QueryPerception::RingBasis
            ))
        ));
        assert!(aromaticity::perceive_aromaticity(
            &mut m,
            AromaticityModel::Mdl,
            Default::default()
        )
        .is_err());
        assert_eq!(m.perception(), &before);
        assert_eq!(
            PreparedTarget::new(&m)
                .all(
                    &smarts::parse_str("[R]", Default::default()).unwrap(),
                    options()
                )
                .unwrap()
                .len(),
            3
        );
    }
}

#[test]
fn topology_reuses_connected_matches_and_shares_budgets() {
    let mut builder = kekule::topology::TopologyBuilder::new();
    let definition = builder
        .add_molecule_definition((molecule("CC")).clone())
        .unwrap();
    for _ in 0..100 {
        builder.add_instance(definition).unwrap();
    }
    let topology = Arc::new(builder.build().unwrap());
    let query = smarts::parse_str("[#6:1]-[#6:2]", Default::default()).unwrap();
    let small_budget = options()
        .with_max_candidate_pairs(4)
        .unwrap()
        .with_max_search_states(20)
        .unwrap();
    assert_eq!(
        PreparedTopologyTarget::new(&topology)
            .all(&query, small_budget)
            .unwrap()
            .len(),
        200
    );
    assert!(matches!(
            PreparedTopologyTarget::new(&topology).all(
                &query,
                small_budget
    .with_max_matches(199)
    .unwrap()
            ),
            Err(SubstructureMatchError::LimitExceeded { limit, .. }) if limit.resource() == Resource::Matches
        ));
}

#[test]
fn recursive_disconnected_environment_can_span_topology_instances() {
    let topology =
        Arc::new(Topology::from_molecules([molecule("C"), molecule("C")].clone()).unwrap());
    let query = smarts::parse_str("[$(C.C)]", Default::default()).unwrap();
    assert_eq!(
        PreparedTopologyTarget::new(&topology)
            .all(&query, options())
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
    let q = smarts::parse_str("[#6:1]-[#8:2]", Default::default()).unwrap();
    let tagged = TaggedQuery::new(&q).unwrap();
    for source in ["[CH3:4][CH2:5][OH:6]", "[OH:6][CH2:5][CH3:4]"] {
        let m = molecule(source);
        let prepared = PreparedTarget::new(&m);
        let matches = tagged
            .find_matches(&prepared, options(), TupleMode::Distinct)
            .unwrap();
        assert_eq!(matches.len(), 1);
        assert_eq!(
            matches[0]
                .iter()
                .map(|id| m.atom(*id).unwrap().atom_map().map(std::num::NonZero::get))
                .collect::<Vec<_>>(),
            vec![Some(5), Some(6)]
        );
        assert_eq!(
            tagged
                .find_matches(&prepared, options(), TupleMode::Distinct)
                .unwrap(),
            matches
        );
    }
    assert_eq!(count("CC", "[0C]"), 2);
}

#[test]
fn graph_complexity_exhaustion_remains_a_resource_error() {
    let input = "*".repeat(4097);
    let error = smarts::parse_str(
        &input,
        SmartsParseOptions::default()
            .with_max_atoms(5000)
            .unwrap()
            .with_max_bonds(5000)
            .unwrap()
            .with_max_total_expression_nodes(65536)
            .unwrap(),
    )
    .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    assert_eq!(error.location(), Some(SourceLocation::at_column(1, 1)));
    // The aggregate query-atom bound trips, not the expression-node bound.
    assert_eq!(
        error.limit_exceeded(),
        Some(LimitExceeded::new(
            Resource::Atoms,
            kekule::query::MAX_QUERY_TOTAL_ATOMS,
            4097
        ))
    );
}
