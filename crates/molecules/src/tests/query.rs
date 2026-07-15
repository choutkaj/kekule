use crate::query::*;
use crate::substructure::{self, *};

use super::*;

fn sanitized(input: &str) -> SmallMolecule {
    SmallMolecule::from_smiles_sanitized(input).expect("test target should sanitize")
}

fn manual_element_query(symbol: &str) -> QueryGraph {
    let mut builder = QueryGraph::builder();
    builder
        .add_atom(AtomExpression::predicate(AtomPredicate::Element(
            Element::from_symbol(symbol).expect("known test element"),
        )))
        .expect("query atom");
    builder.build().expect("non-empty query")
}

#[test]
fn query_graph_is_distinct_from_the_concrete_molecule_kernel() {
    let mut builder = QueryGraphBuilder::with_capacity(2, 1);
    let carbon = builder
        .add_atom(AtomExpression::predicate(AtomPredicate::Element(
            Element::from_symbol("C").unwrap(),
        )))
        .unwrap();
    let hetero = builder
        .add_atom(
            AtomExpression::any([
                AtomExpression::predicate(AtomPredicate::Element(
                    Element::from_symbol("N").unwrap(),
                )),
                AtomExpression::predicate(AtomPredicate::Element(
                    Element::from_symbol("O").unwrap(),
                )),
            ])
            .unwrap(),
        )
        .unwrap();
    let bond = builder
        .add_bond(carbon, hetero, BondExpression::always())
        .unwrap();
    assert_eq!(
        builder.add_bond(carbon, hetero, BondExpression::always()),
        Err(QueryGraphError::DuplicateBond {
            a: carbon,
            b: hetero
        })
    );

    let graph = builder.build().unwrap();
    assert_eq!((graph.atom_count(), graph.bond_count()), (2, 1));
    assert_eq!(graph.bond_between(carbon, hetero).unwrap(), Some(bond));
    assert_eq!(
        graph.neighbors(carbon).unwrap().collect::<Vec<_>>(),
        vec![hetero]
    );
}

#[test]
fn query_expressions_normalize_constants_and_enforce_depth_bounds() {
    let carbon =
        AtomExpression::predicate(AtomPredicate::Element(Element::from_symbol("C").unwrap()));
    let expression = AtomExpression::all([
        AtomExpression::always(),
        AtomExpression::all([carbon.clone()]).unwrap(),
    ])
    .unwrap();
    assert_eq!(expression, carbon);

    let oxygen =
        AtomExpression::predicate(AtomPredicate::Element(Element::from_symbol("O").unwrap()));
    let mut deep = AtomExpression::all([carbon.clone(), oxygen]).unwrap();
    while deep.depth() < MAX_QUERY_EXPRESSION_DEPTH {
        deep = AtomExpression::all([deep.negate().unwrap(), carbon.clone()]).unwrap();
    }
    assert_eq!(deep.depth(), MAX_QUERY_EXPRESSION_DEPTH);
    assert!(matches!(
        deep.negate(),
        Err(QueryExpressionError::ResourceLimit {
            resource: "expression depth",
            ..
        })
    ));
}

#[test]
fn bounded_smarts_builds_branches_rings_components_and_boolean_atoms() {
    let carbonyl = parse_smarts("[#6](=[O,N])-[#7;+1]").unwrap();
    assert_eq!((carbonyl.atom_count(), carbonyl.bond_count()), (3, 2));

    let benzene = parse_smarts("c1ccccc1").unwrap();
    assert_eq!((benzene.atom_count(), benzene.bond_count()), (6, 6));

    let salt = parse_smarts("[Na+].[Cl-]").unwrap();
    assert_eq!((salt.atom_count(), salt.bond_count()), (2, 0));

    let selenium = parse_smarts("[se]1cccc1").unwrap();
    assert_eq!((selenium.atom_count(), selenium.bond_count()), (5, 5));
}

#[test]
fn smarts_logical_precedence_matches_daylight_operator_order() {
    let target = sanitized("CN");

    // High-precedence &: C OR (N AND H2). Both atoms match in methylamine.
    let high_and = parse_smarts("[C,N&H2]").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(target.graph(), &high_and)
            .unwrap()
            .len(),
        2
    );

    // Low-precedence ;: (C OR N) AND H2. Only nitrogen matches.
    let low_and = parse_smarts("[C,N;H2]").unwrap();
    let matches = substructure::find_substructure_matches(target.graph(), &low_and).unwrap();
    assert_eq!(matches.len(), 1);
    assert_eq!(
        target
            .graph()
            .atom(matches[0].atoms()[0])
            .unwrap()
            .element
            .symbol(),
        "N"
    );
}

#[test]
fn smarts_hydrogen_primitive_disambiguation_matches_rdkit() {
    let ethanol = sanitized("CCO");
    for (smarts, expected) in [("[H,D]", 2), ("[C,H]", 3), ("[!H]", 2)] {
        let query = parse_smarts(smarts).unwrap();
        assert_eq!(
            substructure::find_substructure_matches(ethanol.graph(), &query)
                .unwrap()
                .len(),
            expected,
            "{smarts}"
        );
    }

    let mut methane = sanitized("C");
    crate::hydrogens::add_hydrogens(&mut methane).unwrap();
    methane.sanitize().unwrap();
    let elemental_hydrogen = parse_smarts("[H]").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(methane.graph(), &elemental_hydrogen)
            .unwrap()
            .len(),
        4
    );
}

#[test]
fn bounded_smarts_rejects_unsupported_semantics_instead_of_approximating() {
    for (input, expected_fragment) in [
        ("[$(C=O)]", "recursive SMARTS"),
        ("[C@H]", "stereochemical atom"),
        ("C/C=C\\C", "stereochemical bond"),
        ("[C:R]", "atom maps"),
        ("[R2]", "ring-membership counts"),
        ("[r5]", "ring-size predicates"),
        ("[X4]", "connectivity"),
        ("(C.C)", "component-level"),
    ] {
        let error = parse_smarts(input).expect_err(input);
        assert_eq!(error.kind(), SmartsParseErrorKind::Unsupported, "{input}");
        assert!(
            error.message().contains(expected_fragment),
            "{input}: {error}"
        );
        assert!(error.span().start < input.len(), "{input}: {error}");
    }
}

#[test]
fn malformed_smarts_returns_structured_syntax_errors() {
    let empty = parse_smarts("").unwrap_err();
    assert_eq!(empty.kind(), SmartsParseErrorKind::Empty);
    assert_eq!(empty.span(), 0..0);

    for input in ["C(", "C1", "C..C", "C=", "C11", "C1.C1", "[C,N,]"] {
        let error = parse_smarts(input).expect_err(input);
        assert_eq!(
            error.kind(),
            SmartsParseErrorKind::InvalidSyntax,
            "{input}: {error}"
        );
    }
}

#[test]
fn smarts_parser_is_total_over_deterministic_text_mutations() {
    for seed in ["c1ccccc1", "[#6,#7;H1](=O)!@C", "[Na+].[Cl-]"] {
        for input in deterministic_text_mutations(seed) {
            if let Err(error) = parse_smarts(&input) {
                let span = error.span();
                assert!(span.start <= span.end, "{input:?}: {error}");
                assert!(span.end <= input.len(), "{input:?}: {error}");
            }
        }
    }
}

#[test]
fn smarts_limits_cover_input_topology_and_expressions() {
    let input_error = parse_smarts_with_options(
        "CC",
        SmartsParseOptions {
            max_input_bytes: 1,
            ..SmartsParseOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(input_error.kind(), SmartsParseErrorKind::ResourceLimit);

    let atom_error = parse_smarts_with_options(
        "CC",
        SmartsParseOptions {
            max_atoms: 1,
            ..SmartsParseOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(atom_error.kind(), SmartsParseErrorKind::ResourceLimit);

    let branch_error = parse_smarts_with_options(
        "C(C(C))",
        SmartsParseOptions {
            max_branch_depth: 1,
            ..SmartsParseOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(branch_error.kind(), SmartsParseErrorKind::ResourceLimit);

    let expression_error = parse_smarts_with_options(
        "[!C]",
        SmartsParseOptions {
            max_expression_depth: 2,
            ..SmartsParseOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(expression_error.kind(), SmartsParseErrorKind::ResourceLimit);
}

#[test]
fn matcher_handles_elements_bonds_hydrogens_degree_and_negation() {
    let ethanol = sanitized("CCO");
    for (smarts, expected) in [
        ("[#6]", 2),
        ("CO", 1),
        ("C=O", 0),
        ("[OH1]", 1),
        ("[C;D1]", 1),
        ("[!#6]", 1),
        ("[#6]-[#8]", 1),
    ] {
        let query = parse_smarts(smarts).unwrap();
        let matches = substructure::find_substructure_matches(ethanol.graph(), &query).unwrap();
        assert_eq!(matches.len(), expected, "{smarts}: {matches:?}");
    }
}

#[test]
fn matcher_handles_aromatic_cycles_ring_bonds_and_uniqueness() {
    let benzene = sanitized("c1ccccc1");
    let query = parse_smarts("c1ccccc1").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(benzene.graph(), &query)
            .unwrap()
            .len(),
        1
    );
    let embeddings = substructure::find_substructure_matches_with_options(
        benzene.graph(),
        &query,
        SubstructureMatchOptions {
            uniquify: false,
            ..SubstructureMatchOptions::default()
        },
    )
    .unwrap();
    assert_eq!(embeddings.len(), 12);

    let cyclohexane = sanitized("C1CCCCC1");
    let ring_atoms = parse_smarts("[#6;R]").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(cyclohexane.graph(), &ring_atoms)
            .unwrap()
            .len(),
        6
    );
    let ring_bond = parse_smarts("C@C").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(cyclohexane.graph(), &ring_bond)
            .unwrap()
            .len(),
        6
    );
}

#[test]
fn matcher_supports_disconnected_queries_and_non_induced_subgraphs() {
    let disconnected = sanitized("O.O");
    let query = parse_smarts("[#8].[#8]").unwrap();
    let matches = substructure::find_substructure_matches(disconnected.graph(), &query).unwrap();
    assert_eq!(matches.len(), 1);
    assert_ne!(matches[0].atoms()[0], matches[0].atoms()[1]);

    let cyclopropane = sanitized("C1CC1");
    let edge = parse_smarts("C-C").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(cyclopropane.graph(), &edge)
            .unwrap()
            .len(),
        3
    );
}

#[test]
fn matcher_requires_only_the_perception_used_by_the_ir() {
    let mut raw = Molecule::new();
    let first = raw.add_atom(carbon());
    let second = raw.add_atom(carbon());
    raw.add_bond(first, second, BondOrder::Single).unwrap();
    let elemental = manual_element_query("C");
    assert_eq!(
        substructure::find_substructure_matches(&raw, &elemental)
            .unwrap()
            .len(),
        2
    );

    for (smarts, perception) in [
        ("C", QueryPerception::Aromaticity),
        ("[R]", QueryPerception::RingMembership),
        ("[H1]", QueryPerception::Valence),
    ] {
        let query = parse_smarts(smarts).unwrap();
        assert_eq!(
            substructure::find_substructure_matches(&raw, &query),
            Err(SubstructureMatchError::MissingPerception(perception)),
            "{smarts}"
        );
    }
}

#[test]
fn matcher_search_and_candidate_limits_are_hard_failures() {
    let target = sanitized("CCCC");
    let query = parse_smarts("[#6]").unwrap();
    let candidate_error = substructure::find_substructure_matches_with_options(
        target.graph(),
        &query,
        SubstructureMatchOptions {
            max_candidate_pairs: 3,
            ..SubstructureMatchOptions::default()
        },
    )
    .unwrap_err();
    assert!(matches!(
        candidate_error,
        SubstructureMatchError::ResourceLimit {
            resource: "candidate pairs",
            ..
        }
    ));

    let state_error = substructure::find_substructure_matches_with_options(
        target.graph(),
        &query,
        SubstructureMatchOptions {
            max_search_states: 1,
            max_matches: 4,
            ..SubstructureMatchOptions::default()
        },
    )
    .unwrap_err();
    assert!(matches!(
        state_error,
        SubstructureMatchError::ResourceLimit {
            resource: "search states",
            ..
        }
    ));

    assert!(matches!(
        substructure::find_substructure_matches_with_options(
            target.graph(),
            &query,
            SubstructureMatchOptions {
                max_query_atoms: MAX_SUBSTRUCTURE_QUERY_ATOMS + 1,
                ..SubstructureMatchOptions::default()
            },
        ),
        Err(SubstructureMatchError::InvalidOptions(_))
    ));
}
