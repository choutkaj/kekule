use crate::query::*;
use crate::substructure::{self, *};

use super::*;

fn perceived(input: &str) -> Molecule {
    let mut molecule = read_smiles(input).expect("test target should parse");
    perceive(&mut molecule).expect("test target should perceive");
    molecule
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
fn tetrahedral_smarts_matches_local_carrier_order_and_partial_environments() {
    // These are local configurations, independent of any CIP assignment.
    for (query, target, expected) in [
        ("N[C@H](F)Cl", "N[C@H](F)Cl", true),
        ("N[C@H](F)Cl", "N[C@@H](F)Cl", false),
        ("N[C@H](F)Cl", "NC(F)Cl", false),
        ("NC(F)Cl", "N[C@@H](F)Cl", true),
        ("N[C@H](F)Cl", "N[C@](F)(Cl)Br", false),
        ("N[C@](F)Cl", "N[C@](F)(Cl)Br", true),
        ("N[C@](F)Cl", "N[C@@](F)(Cl)Br", false),
        ("N[C@](F)Cl", "N[C@H](F)Cl", true),
        ("N[C@TH1H](F)Cl", "N[C@H](F)Cl", true),
        ("N[C@TH2H](F)Cl", "N[C@H](F)Cl", false),
        ("[C@H](N)(F)Cl", "N[C@@H](F)Cl", true),
        ("[C@H](N)(F)Cl", "N[C@H](F)Cl", false),
        ("[C@;H1](N)(F)Cl", "N[C@H](F)Cl", true),
        ("[C;H1;@](N)(F)Cl", "N[C@H](F)Cl", true),
        ("[C@H0](N)(F)Cl", "N[C@](F)(Cl)Br", true),
        ("[C@H]([H])(F)(Cl)Br", "[H][C@](F)(Cl)Br", true),
        ("[C@H]([H])(F)(Cl)Br", "[H][C@@](F)(Cl)Br", false),
        ("[C@H]([H])(F)Cl", "[H][C@@](F)(Cl)Br", true),
        ("[C@]", "N[C@H](F)Cl", true),
        ("[C@]", "N[C@@H](F)Cl", true),
        ("[C@]", "NC(F)Cl", false),
        ("N[C@]F", "N[C@@H](F)Cl", true),
        ("F[C@]1(Cl)CCCO1", "F[C@]1(Cl)CCCO1", true),
        ("F[C@]1(Cl)CCCO1", "F[C@@]1(Cl)CCCO1", false),
        ("F[C@]1(Cl)CCCO1", "O1CCC[C@@]1(F)Cl", false),
        ("O1CCC[C@]1(F)Cl", "F[C@]1(Cl)CCCO1", true),
    ] {
        let target = perceived(target);
        let query_graph = parse_smarts(query).unwrap();
        assert_eq!(
            find_substructure_match(&target, &query_graph)
                .unwrap()
                .is_some(),
            expected,
            "{query}"
        );
    }
}

#[test]
fn directional_smarts_matches_selected_carriers_in_either_endpoint_order() {
    for (query, target, expected) in [
        ("F/C=C/Cl", "F/C=C/Cl", true),
        ("F/C=C/Cl", "F/C=C\\Cl", false),
        ("F/C=C/Cl", "FC=CCl", false),
        ("FC=CCl", "F/C=C\\Cl", true),
        ("F/C=C/Cl", "Cl/C=C/F", true),
        ("Cl/C=C/F", "F/C=C/Cl", true),
        ("F/C=C/Cl", "F/C(Br)=C(/Cl)I", true),
        ("Br/C=C/Cl", "F/C(Br)=C(/Cl)I", false),
        ("Br/C=C\\Cl", "F/C(Br)=C(/Cl)I", true),
        ("F/C(/Br)=C/Cl", "F/C(Br)=C(/Cl)I", true),
        ("F/C=C", "F/C=C/Cl", true),
        ("F/C=C", "F/C=C\\Cl", true),
        ("F/C=C", "FC=CCl", true),
        ("C/C", "CC", true),
        ("c/c", "c1ccccc1", true),
        ("c-c", "c1ccccc1", false),
        ("F/C=C\\1CCCNCC1", "F/C=C\\1CCCNCC1", true),
        ("F/C=C\\1CCCNCC1", "F/C=C1CCCNCC/1", true),
        ("F/C=C\\1CCCNCC1", "F/C=C1CCCNCC\\1", false),
        ("F/C=C/C=C/Cl", "F/C=C/C=C/Cl", true),
        ("F/C=C/C=C/Cl", "F/C=C/C=C\\Cl", false),
    ] {
        let target = perceived(target);
        let query_graph = parse_smarts(query).unwrap();
        assert_eq!(
            find_substructure_match(&target, &query_graph)
                .unwrap()
                .is_some(),
            expected,
            "{query}"
        );
    }
}

#[test]
fn stereo_mapping_filter_precedes_uniqueness_and_match_limit() {
    let query = parse_smarts("[*@](*)(*)(*)*").unwrap();
    for input in ["[C@](N)(F)(Cl)Br", "[C@@](N)(F)(Cl)Br"] {
        let target = perceived(input);
        let all = find_substructure_matches_with_options(
            &target,
            &query,
            SubstructureMatchOptions {
                uniquify: false,
                ..SubstructureMatchOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            all.len(),
            12,
            "half of the 24 carrier permutations satisfy parity"
        );
        assert_eq!(find_substructure_matches(&target, &query).unwrap().len(), 1);
        assert_eq!(
            find_substructure_match(&target, &query).unwrap().unwrap(),
            all[0]
        );
    }
}

#[test]
fn query_stereo_builder_checks_ownership_and_revalidates_later_topology() {
    let mut builder = QueryGraph::builder();
    let center = builder.add_atom(AtomExpression::always()).unwrap();
    let carriers: Vec<_> = (0..3)
        .map(|_| builder.add_atom(AtomExpression::always()).unwrap())
        .collect();
    for carrier in &carriers {
        builder
            .add_bond(center, *carrier, BondExpression::always())
            .unwrap();
    }
    let constraint = QueryStereoConstraint::Tetrahedral {
        center,
        carriers: carriers.clone(),
        orientation: TetrahedralOrientation::Clockwise,
    };
    let mut invalid = constraint.clone();
    if let QueryStereoConstraint::Tetrahedral { carriers, .. } = &mut invalid {
        carriers[0] = center;
    }
    assert!(matches!(
        builder.add_stereo_constraint(invalid),
        Err(QueryGraphError::InvalidStereo(_))
    ));
    let mut invalid = constraint.clone();
    if let QueryStereoConstraint::Tetrahedral { center, .. } = &mut invalid {
        *center = QueryAtomId::new(99);
    }
    assert!(matches!(
        builder.add_stereo_constraint(invalid),
        Err(QueryGraphError::InvalidAtomId(_))
    ));
    let mut permuted = builder.clone();
    let mut alternate = constraint.clone();
    if let QueryStereoConstraint::Tetrahedral {
        carriers,
        orientation,
        ..
    } = &mut alternate
    {
        carriers.swap(0, 1);
        *orientation = orientation.inverted();
    }
    permuted.add_stereo_constraint(alternate).unwrap();
    builder.add_stereo_constraint(constraint.clone()).unwrap();
    assert_eq!(builder.clone().build().unwrap(), permuted.build().unwrap());
    assert!(matches!(
        builder.add_stereo_constraint(constraint),
        Err(QueryGraphError::InvalidStereo(_))
    ));

    // Matching represented stereo needs no perception when other predicates do not.
    let query = builder.clone().build().unwrap();
    let target = read_smiles("N[C@H](F)Cl").unwrap();
    assert!(find_substructure_match(&target, &query).unwrap().is_some());
    let extra = builder.add_atom(AtomExpression::always()).unwrap();
    builder
        .add_bond(center, extra, BondExpression::always())
        .unwrap();
    assert!(matches!(
        builder.build(),
        Err(QueryGraphError::InvalidStereo(_))
    ));
}

#[test]
fn double_bond_query_builder_requires_live_adjacent_carriers() {
    let mut builder = QueryGraph::builder();
    let atoms: Vec<_> = (0..4)
        .map(|_| builder.add_atom(AtomExpression::always()).unwrap())
        .collect();
    builder
        .add_bond(atoms[0], atoms[1], BondExpression::always())
        .unwrap();
    let focus = builder
        .add_bond(atoms[1], atoms[2], BondExpression::always())
        .unwrap();
    builder
        .add_bond(atoms[2], atoms[3], BondExpression::always())
        .unwrap();
    for (bond, left_carrier, right_carrier) in [
        (QueryBondId::new(99), atoms[0], atoms[3]),
        (focus, QueryAtomId::new(99), atoms[3]),
        (focus, atoms[1], atoms[3]),
        (focus, atoms[3], atoms[0]),
        (focus, atoms[0], atoms[0]),
    ] {
        assert!(builder
            .add_stereo_constraint(QueryStereoConstraint::DoubleBond {
                bond,
                left_carrier,
                right_carrier,
                orientation: DoubleBondOrientation::Opposite,
            })
            .is_err());
    }
    builder
        .add_stereo_constraint(QueryStereoConstraint::DoubleBond {
            bond: focus,
            left_carrier: atoms[0],
            right_carrier: atoms[3],
            orientation: DoubleBondOrientation::Opposite,
        })
        .unwrap();
    let query = builder.build().unwrap();
    assert!(
        find_substructure_match(&read_smiles("F/C=C/Cl").unwrap(), &query)
            .unwrap()
            .is_some()
    );
    assert!(
        find_substructure_match(&read_smiles("F/C=C\\Cl").unwrap(), &query)
            .unwrap()
            .is_none()
    );
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
    let target = perceived("CN");

    // High-precedence &: C OR (N AND H2). Both atoms match in methylamine.
    let high_and = parse_smarts("[C,N&H2]").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(&target, &high_and)
            .unwrap()
            .len(),
        2
    );

    // Low-precedence ;: (C OR N) AND H2. Only nitrogen matches.
    let low_and = parse_smarts("[C,N;H2]").unwrap();
    let matches = substructure::find_substructure_matches(&target, &low_and).unwrap();
    assert_eq!(matches.len(), 1);
    assert_eq!(
        target.atom(matches[0].atoms()[0]).unwrap().element.symbol(),
        "N"
    );
}

#[test]
fn unbracketed_smarts_elements_do_not_consume_a_following_aromatic_atom() {
    for (input, atoms) in [
        ("Cn1cccc1", 6),
        ("Sc1ccccc1", 7),
        ("Sn1cccc1", 6),
        ("Clc1ccccc1", 7),
        ("Brc1ccccc1", 7),
    ] {
        let query = parse_smarts(input).unwrap();
        assert_eq!(query.atom_count(), atoms, "{input}");
        assert!(
            find_substructure_match(&perceived(input), &query)
                .unwrap()
                .is_some(),
            "{input}"
        );
    }
    assert_eq!(parse_smarts("[Cn]").unwrap().atom_count(), 1);
    assert_eq!(
        parse_smarts("Ca").unwrap().atom_count(),
        2,
        "a is an aromatic wildcard outside brackets"
    );
}

#[test]
fn smarts_hydrogen_primitive_disambiguation_matches_rdkit() {
    let ethanol = perceived("CCO");
    for (smarts, expected) in [("[H,D]", 2), ("[C,H]", 3), ("[!H]", 2)] {
        let query = parse_smarts(smarts).unwrap();
        assert_eq!(
            substructure::find_substructure_matches(&ethanol, &query)
                .unwrap()
                .len(),
            expected,
            "{smarts}"
        );
    }

    let mut methane = perceived("C");
    crate::hydrogens::add_hydrogens(&mut methane).unwrap();
    perceive(&mut methane).unwrap();
    let elemental_hydrogen = parse_smarts("[H]").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(&methane, &elemental_hydrogen)
            .unwrap()
            .len(),
        4
    );
}

#[test]
fn bounded_smarts_rejects_unsupported_semantics_instead_of_approximating() {
    for (input, expected_fragment) in [
        ("[C@AL1]", "stereochemical atom"),
        ("[C@?]", "stereochemical atom"),
        ("C/?C=C\\C", "unspecified directional"),
        ("[^2]", "hybridization"),
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

    for input in [
        "C(",
        "C1",
        "C..C",
        "C=",
        "C11",
        "C1.C1",
        "[C,N,]",
        "[C:R]",
        "[C@@@]",
        "[C@](F)(Cl)(Br)(I)N",
        "F/C(\\Br)=C/Cl",
        "F/C(/Cl)(/Br)=C/I",
        "C/1CCCCC/1",
    ] {
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
    for seed in [
        "c1ccccc1",
        "[#6,#7;H1](=O)!@C",
        "[Na+].[Cl-]",
        "N[C@H]1CCCO1",
        "F/C(Br)=C(/Cl)I",
    ] {
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
    let ethanol = perceived("CCO");
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
        let matches = substructure::find_substructure_matches(&ethanol, &query).unwrap();
        assert_eq!(matches.len(), expected, "{smarts}: {matches:?}");
    }
}

#[test]
fn total_connectivity_counts_graph_declared_and_inferred_hydrogens() {
    // Counts independently checked with RDKit 2026.03.3. Degree counts graph
    // neighbors; connectivity also counts nongraph H, but not radical electrons.
    for (target, smarts, expected) in [
        ("C", "[X4;D0]", 1),
        ("[CH4]", "[X4;D0]", 1),
        ("[CH3]", "[X3]", 1),
        ("[NH4+]", "[X4]", 1),
        ("CCO", "[X4]", 2),
        ("CCO", "[!X4]", 1),
        ("CCO", "[X2,X4]", 3),
        ("C#N", "[X2]", 1),
        ("C#N", "[X]", 1),
        ("[2H]C", "[X4]", 1),
        ("[2H]C", "[X1]", 1),
        ("c1cc[nH]c1", "[X3]", 5),
        ("[Xe]", "[Xe;X0]", 1),
    ] {
        let molecule = perceived(target);
        let query = parse_smarts(smarts).unwrap();
        let matches = find_substructure_matches(&molecule, &query).unwrap();
        assert_eq!(matches.len(), expected, "{target} / {smarts}");
    }
}

#[test]
fn connectivity_includes_zero_and_dative_graph_neighbors() {
    for order in [BondOrder::Zero, BondOrder::Dative] {
        let mut editor = crate::core::MoleculeEditor::new();
        let nitrogen = editor
            .add_atom(Atom::new(Element::from_symbol("N").unwrap()))
            .unwrap();
        let copper = editor
            .add_atom(Atom::new(Element::from_symbol("Cu").unwrap()))
            .unwrap();
        editor.add_bond(nitrogen, copper, order).unwrap();
        let mut molecule = editor.try_finish().unwrap();
        perceive(&mut molecule).unwrap();
        for (smarts, expected) in [("[X4]", nitrogen), ("[X]", copper)] {
            let matches =
                find_substructure_matches(&molecule, &parse_smarts(smarts).unwrap()).unwrap();
            assert_eq!(matches.len(), 1, "{order:?} / {smarts}");
            assert_eq!(matches[0].atoms(), &[expected]);
        }
        assert_eq!(
            find_substructure_matches(&molecule, &parse_smarts("[x0]").unwrap())
                .unwrap()
                .len(),
            2
        );
    }
}

#[test]
fn connectivity_and_hydrogen_queries_survive_hydrogen_materialization() {
    for target in [
        "C",
        "[CH4]",
        "[CH3]",
        "[NH4+]",
        "CCO",
        "[2H]C",
        "c1cc[nH]c1",
    ] {
        let mut molecule = perceived(target);
        let queries = ["[!#1;X3]", "[!#1;X4]", "[!#1;X2]", "[!#1;H3]", "[!#1;H1]"];
        let before: Vec<_> = queries
            .iter()
            .map(|smarts| {
                find_substructure_matches(&molecule, &parse_smarts(smarts).unwrap()).unwrap()
            })
            .collect();
        crate::hydrogens::add_hydrogens(&mut molecule).unwrap();
        perceive(&mut molecule).unwrap();
        for (smarts, expected) in queries.into_iter().zip(before) {
            let actual =
                find_substructure_matches(&molecule, &parse_smarts(smarts).unwrap()).unwrap();
            assert_eq!(actual, expected, "{target} / {smarts}");
        }
    }
}

#[test]
fn ring_bond_count_uses_cycle_membership_without_a_selected_ring_basis() {
    // Simple, fused, bridged, spiro and cage graphs distinguish cyclic degree
    // from either total degree or the number of selected rings through an atom.
    for (target, counts) in [
        ("CC", [2, 0, 0, 0, 0]),
        ("CC1CC1", [1, 0, 3, 0, 0]),
        ("c1ccc2ccccc2c1", [0, 0, 8, 2, 0]),
        ("C1CC2CCC1C2", [0, 0, 5, 2, 0]),
        ("C1CCC2(CC1)CCCC2", [0, 0, 9, 0, 1]),
        ("C12C3C4C1C5C2C3C45", [0, 0, 0, 8, 0]),
    ] {
        let mut molecule = read_smiles(target).unwrap();
        crate::perception::rings::perceive_ring_membership(&mut molecule);
        assert!(molecule.ring_set().is_none());
        for (count, expected) in counts.into_iter().enumerate() {
            let smarts = format!("[x{count}]");
            let actual =
                find_substructure_matches(&molecule, &parse_smarts(&smarts).unwrap()).unwrap();
            assert_eq!(actual.len(), expected, "{target} / {smarts}");
        }
        for (smarts, expected) in [
            ("[x]", molecule.atom_count() - counts[0]),
            ("[!x]", counts[0]),
        ] {
            let actual =
                find_substructure_matches(&molecule, &parse_smarts(smarts).unwrap()).unwrap();
            assert_eq!(actual.len(), expected, "{target} / {smarts}");
        }
    }
}

#[test]
fn connectivity_and_ring_count_overflow_are_errors() {
    for smarts in [
        "[X256]",
        "[x256]",
        "[X99999999999999999999]",
        "[x99999999999999999999]",
    ] {
        assert!(parse_smarts(smarts).is_err(), "{smarts}");
    }
    for smarts in ["[X0]", "[x0]", "[X255]", "[x255]"] {
        assert!(parse_smarts(smarts).is_ok(), "{smarts}");
    }
}

#[test]
fn matcher_handles_aromatic_cycles_ring_bonds_and_uniqueness() {
    let benzene = perceived("c1ccccc1");
    let query = parse_smarts("c1ccccc1").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(&benzene, &query)
            .unwrap()
            .len(),
        1
    );
    let embeddings = substructure::find_substructure_matches_with_options(
        &benzene,
        &query,
        SubstructureMatchOptions {
            uniquify: false,
            ..SubstructureMatchOptions::default()
        },
    )
    .unwrap();
    assert_eq!(embeddings.len(), 12);

    let cyclohexane = perceived("C1CCCCC1");
    let ring_atoms = parse_smarts("[#6;R]").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(&cyclohexane, &ring_atoms)
            .unwrap()
            .len(),
        6
    );
    let ring_bond = parse_smarts("C@C").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(&cyclohexane, &ring_bond)
            .unwrap()
            .len(),
        6
    );
}

#[test]
fn matcher_supports_disconnected_queries_and_non_induced_subgraphs() {
    let connected_target = perceived("OCO");
    let query = parse_smarts("[#8].[#8]").unwrap();
    let matches = substructure::find_substructure_matches(&connected_target, &query).unwrap();
    assert_eq!(matches.len(), 1);
    assert_ne!(matches[0].atoms()[0], matches[0].atoms()[1]);

    let cyclopropane = perceived("C1CC1");
    let edge = parse_smarts("C-C").unwrap();
    assert_eq!(
        substructure::find_substructure_matches(&cyclopropane, &edge)
            .unwrap()
            .len(),
        3
    );
}

#[test]
fn smarts_bonds_distinguish_aromatic_types_from_localized_orders() {
    // Counts include both query-to-target orientations of every matching edge.
    let queries = ["*-*", "*=*", "*#*", "*:*", "**", "*~*"];
    for (source, expected) in [
        ("c1ccoc1", [0, 0, 0, 10, 10, 10]),
        ("c1ccccc1-c2ccccc2", [2, 0, 0, 24, 26, 26]),
        ("C1=CC=CC#C1", [0, 0, 2, 10, 10, 12]),
        ("C1=CCCCC1", [10, 2, 0, 0, 10, 12]),
        ("CC#CC", [4, 0, 2, 0, 4, 6]),
    ] {
        let target = perceived(source);
        for (smarts, count) in queries.into_iter().zip(expected) {
            let query = parse_smarts(smarts).unwrap();
            let matches = find_substructure_matches_with_options(
                &target,
                &query,
                SubstructureMatchOptions {
                    uniquify: false,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(matches.len(), count, "{source}: {smarts}");
        }
    }
    let aryne = perceived("C1=CC=CC#C1");
    assert!(
        find_substructure_matches(&aryne, &parse_smarts("c1ccccc1").unwrap())
            .unwrap()
            .is_empty()
    );
    let benzene = perceived("c1ccccc1");
    for smarts in ["c1ccccc-1", "c-1ccccc1"] {
        assert!(
            find_substructure_matches(&benzene, &parse_smarts(smarts).unwrap())
                .unwrap()
                .is_empty()
        );
    }
    assert!(parse_smarts("C=1CCCCC-1").is_err());
}

#[test]
fn programmatic_bond_predicates_keep_their_represented_meaning() {
    let make_query = |predicate| {
        let mut builder = QueryGraph::builder();
        let a = builder.add_atom(AtomExpression::always()).unwrap();
        let b = builder.add_atom(AtomExpression::always()).unwrap();
        builder
            .add_bond(a, b, BondExpression::predicate(predicate))
            .unwrap();
        builder.build().unwrap()
    };
    let mut raw = read_smiles("CC").unwrap();
    let order = make_query(BondPredicate::Order(BondOrder::Single));
    assert_eq!(find_substructure_matches(&raw, &order).unwrap().len(), 1);
    assert_eq!(
        find_substructure_matches(&raw, &parse_smarts("*-*").unwrap()),
        Err(SubstructureMatchError::MissingPerception(
            QueryPerception::Aromaticity
        ))
    );
    perceive(&mut raw).unwrap();
    assert_eq!(
        find_substructure_matches(&raw, &parse_smarts("*-*").unwrap())
            .unwrap()
            .len(),
        1
    );

    let aryne = perceived("C1=CC=CC#C1");
    let membership = make_query(BondPredicate::Aromatic(true));
    assert_eq!(
        find_substructure_matches(&aryne, &membership)
            .unwrap()
            .len(),
        6
    );
    assert_eq!(
        find_substructure_matches(&aryne, &parse_smarts("*:*").unwrap())
            .unwrap()
            .len(),
        5
    );
    let benzene = perceived("c1ccccc1");
    assert_eq!(
        find_substructure_matches(&benzene, &order).unwrap().len(),
        3
    );
}

#[test]
fn matcher_requires_only_the_perception_used_by_the_ir() {
    let mut raw = crate::core::MoleculeEditor::new();
    let first = raw.add_atom(carbon()).expect("atom identifier capacity");
    let second = raw.add_atom(carbon()).expect("atom identifier capacity");
    raw.add_bond(first, second, BondOrder::Single).unwrap();
    let elemental = manual_element_query("C");
    assert_eq!(
        substructure::find_substructure_matches(raw.working(), &elemental)
            .unwrap()
            .len(),
        2
    );

    for (smarts, perception) in [
        ("C", QueryPerception::Aromaticity),
        ("[R]", QueryPerception::RingMembership),
        ("[H1]", QueryPerception::Valence),
        ("[X4]", QueryPerception::Valence),
        ("[!X4]", QueryPerception::Valence),
        ("[x2]", QueryPerception::RingMembership),
        ("[!x0]", QueryPerception::RingMembership),
    ] {
        let query = parse_smarts(smarts).unwrap();
        assert_eq!(
            substructure::find_substructure_matches(raw.working(), &query),
            Err(SubstructureMatchError::MissingPerception(perception)),
            "{smarts}"
        );
    }
}

#[test]
fn matcher_search_and_candidate_limits_are_hard_failures() {
    let target = perceived("CCCC");
    let query = parse_smarts("[#6]").unwrap();
    let candidate_error = substructure::find_substructure_matches_with_options(
        &target,
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
        &target,
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
            &target,
            &query,
            SubstructureMatchOptions {
                max_query_atoms: MAX_SUBSTRUCTURE_QUERY_ATOMS + 1,
                ..SubstructureMatchOptions::default()
            },
        ),
        Err(SubstructureMatchError::InvalidOptions(_))
    ));
}
