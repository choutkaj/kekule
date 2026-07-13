use super::*;

#[test]
fn smiles_document_preserves_spans_and_dot_boundaries_before_interpretation() {
    let input = "[Na+].[Cl-]";
    let document = smiles_api::parse_str(input).expect("document parses");
    assert_eq!(document.source(), input);
    assert_eq!(document.component_token_ranges().len(), 2);
    assert!(document.tokens().iter().all(|token| {
        let span = token.span();
        span.start <= span.end && span.end <= input.len()
    }));
    let molecule = smiles_api::interpret(&document).expect("document interprets");
    assert_eq!(molecule.atom_count(), 2);
    assert_eq!(molecule.bond_count(), 0);
    assert_eq!(molecule.graph().connected_components().len(), 2);
    assert!(!molecule.graph().perception().has_valence());
}

#[test]
fn smiles_parses_branches_rings_brackets_and_fragments_without_sanitizing() {
    let small =
        read_smiles("C(C)O.C1=CC=CC=C1.[13NH4+:7].[C@@H](N)O").expect("smiles should parse");

    assert_eq!(small.graph().atom_count(), 13);
    assert_eq!(small.graph().bond_count(), 10);
    assert!(!small.graph().perception().has_valence());
    let bracket_atom = small.graph().atom(AtomId::new(9)).expect("bracket atom");
    assert_eq!(bracket_atom.isotope, Some(13));
    assert_eq!(bracket_atom.explicit_hydrogens, 4);
    assert!(bracket_atom.no_implicit_hydrogens);
    assert_eq!(bracket_atom.formal_charge, 1);
    assert_eq!(bracket_atom.atom_map, Some(7));
    let chiral_atom = small
        .graph()
        .atom(AtomId::new(10))
        .expect("chiral bracket atom");
    assert_eq!(chiral_atom.explicit_hydrogens, 1);
    let stereo = small
        .graph()
        .stereo_elements()
        .map(|(_, element)| element)
        .collect::<Vec<_>>();
    assert_eq!(stereo.len(), 1);
    match &stereo[0].kind {
        StereoElementKind::Tetrahedral(tetrahedral) => {
            assert_eq!(tetrahedral.center, AtomId::new(10));
            assert_eq!(
                tetrahedral.orientation,
                TetrahedralOrientation::CounterClockwise
            );
            assert!(tetrahedral
                .carriers
                .contains(&StereoCarrier::ImplicitHydrogen));
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn smiles_parses_directional_bond_markers_without_sanitizing_stereo() {
    let small = read_smiles("C/C=C\\C").expect("directional bond markers should parse");

    assert_eq!(small.graph().atom_count(), 4);
    assert_eq!(small.graph().bond_count(), 3);
    assert_eq!(
        small
            .graph()
            .stereo_bond_mark(BondId::new(0))
            .expect("first directional mark")
            .kind,
        StereoBondMarkKind::DirectionalUp
    );
    assert_eq!(
        small
            .graph()
            .stereo_bond_mark(BondId::new(2))
            .expect("second directional mark")
            .kind,
        StereoBondMarkKind::DirectionalDown
    );
    let canonical = smiles_api::write_canonical_with_options(&small, CanonicalSmilesWriteOptions)
        .expect("non-isomeric canonical SMILES should ignore directional bond markers");
    let mut reparsed = read_smiles(&canonical).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .expect("canonical output should sanitize");
}

#[test]
fn metal_bound_organic_subset_halogen_keeps_rdkit_no_implicit_state() {
    let mut small = read_smiles("Br[Pt+2]Br").expect("platinum bromide salt parses");
    perception_api::sanitize_with_options(&mut small, SanitizeOptions::default())
        .expect("platinum bromide salt sanitizes");

    let bromines = small
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "Br")
        .map(|(_, atom)| {
            (
                atom.no_implicit_hydrogens,
                atom.implicit_hydrogens.unwrap_or(0),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(bromines, vec![(false, 0), (false, 0)]);

    let mut aryl_bromide = read_smiles("c1ccccc1Br").expect("aryl bromide should parse");
    perception_api::sanitize_with_options(&mut aryl_bromide, SanitizeOptions::default())
        .expect("aryl bromide should sanitize");
    let bromine = aryl_bromide
        .graph()
        .atoms()
        .find_map(|(_, atom)| (atom.element.symbol() == "Br").then_some(atom))
        .expect("bromine atom");
    assert!(!bromine.no_implicit_hydrogens);
}

#[test]
fn metal_bound_organic_subset_atoms_rely_on_valence_hydrogens() {
    let mut aryl_mercury = read_smiles("c1ccccc1[Hg]").expect("aryl mercury should parse");
    perception_api::sanitize_with_options(&mut aryl_mercury, SanitizeOptions::default())
        .expect("aryl mercury should sanitize");
    let aryl_mercury_carbon = aryl_mercury
        .graph()
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "C"
                && aryl_mercury.graph().incident_bonds(id).is_ok_and(|bonds| {
                    bonds.into_iter().any(|(_, bond)| {
                        aryl_mercury
                            .graph()
                            .atom(bond.other_atom(id))
                            .is_ok_and(|neighbor| neighbor.element.symbol() == "Hg")
                    })
                }))
            .then_some(atom)
        })
        .expect("aryl carbon bound to mercury");
    assert!(!aryl_mercury_carbon.no_implicit_hydrogens);
    assert_eq!(aryl_mercury_carbon.implicit_hydrogens, Some(0));

    let methyl_sodium = read_smiles("C[Na]").expect("methyl sodium should parse");
    let carbon = methyl_sodium
        .graph()
        .atoms()
        .find_map(|(_, atom)| (atom.element.symbol() == "C").then_some(atom))
        .expect("carbon atom");
    assert!(!carbon.no_implicit_hydrogens);
    assert_eq!(carbon.implicit_hydrogens, None);
}

#[test]
fn aromatic_chalcogen_bracket_atoms_parse_without_sanitizing() {
    let small = read_smiles("[se]1cccc1.[te]1cccc1")
        .expect("aromatic selenium and tellurium bracket atoms should parse");

    let chalcogens = small
        .graph()
        .atoms()
        .filter(|(_, atom)| matches!(atom.element.symbol(), "Se" | "Te"))
        .map(|(_, atom)| {
            (
                atom.element.symbol().to_owned(),
                atom.aromatic,
                atom.no_implicit_hydrogens,
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        chalcogens,
        vec![("Se".to_owned(), true, true), ("Te".to_owned(), true, true)]
    );
    assert!(small.graph().perception().has_aromaticity());
}

#[test]
fn malformed_smiles_returns_errors_without_panicking() {
    let cases = [
        "C(",
        "C1",
        "C%1",
        "C%a1",
        "C=",
        "=C",
        "C..C",
        "C1.C1",
        "C=1CCCCC-1",
        "[]",
        "[13]",
        "[é]",
        "[C@@@H]",
        "[C/]",
        "[*]",
        "[C+999]",
        "[C:]",
        "[Clx]",
        "[si]1ccccc1",
        "Cé",
    ];

    for input in cases {
        let parsed = std::panic::catch_unwind(|| read_smiles(input))
            .unwrap_or_else(|_| panic!("`{input}` panicked"));
        let error = parsed.expect_err("malformed SMILES should fail");
        assert!(!error.to_string().is_empty(), "message for `{input}`");
    }
}

#[test]
fn smiles_writer_round_trips_graph_shape() {
    let small = read_smiles("CC(=O)O").expect("smiles should parse");
    let text =
        smiles_api::write_with_options(&small, SmilesWriteOptions).expect("smiles should write");
    let reparsed = read_smiles(&text).expect("written smiles should parse");

    assert_eq!(reparsed.graph().atom_count(), small.graph().atom_count());
    assert_eq!(reparsed.graph().bond_count(), small.graph().bond_count());
}

#[test]
fn canonical_smiles_is_stable_across_atom_order_for_tree_roles() {
    let mut first = SmallMolecule::new();
    let first_terminal_a = first.graph_mut().add_atom(carbon());
    let first_center = first.graph_mut().add_atom(carbon());
    let first_terminal_b = first.graph_mut().add_atom(carbon());
    first
        .graph_mut()
        .add_bond(first_terminal_a, first_center, BondOrder::Single)
        .expect("bond should be valid");
    first
        .graph_mut()
        .add_bond(first_center, first_terminal_b, BondOrder::Single)
        .expect("bond should be valid");
    perception_api::sanitize_with_options(&mut first, SanitizeOptions::default())
        .expect("propane sanitizes");

    let mut second = SmallMolecule::new();
    let second_center = second.graph_mut().add_atom(carbon());
    let second_terminal_a = second.graph_mut().add_atom(carbon());
    let second_terminal_b = second.graph_mut().add_atom(carbon());
    second
        .graph_mut()
        .add_bond(second_center, second_terminal_a, BondOrder::Single)
        .expect("bond should be valid");
    second
        .graph_mut()
        .add_bond(second_center, second_terminal_b, BondOrder::Single)
        .expect("bond should be valid");
    perception_api::sanitize_with_options(&mut second, SanitizeOptions::default())
        .expect("propane sanitizes");

    let first_written =
        smiles_api::write_canonical_with_options(&first, CanonicalSmilesWriteOptions)
            .expect("canonical SMILES should write");
    let second_written =
        smiles_api::write_canonical_with_options(&second, CanonicalSmilesWriteOptions)
            .expect("canonical SMILES should write");

    assert_eq!(first_written, second_written);
    assert_eq!(first_written, "CCC");
    read_smiles(&first_written).expect("canonical output should parse");
}

#[test]
fn canonical_smiles_sorts_disconnected_components() {
    let mut first = read_smiles("O.C").expect("SMILES parses");
    let mut second = read_smiles("C.O").expect("SMILES parses");
    perception_api::sanitize_with_options(&mut first, SanitizeOptions::default())
        .expect("first sanitizes");
    perception_api::sanitize_with_options(&mut second, SanitizeOptions::default())
        .expect("second sanitizes");

    assert_eq!(
        smiles_api::write_canonical_with_options(&first, CanonicalSmilesWriteOptions)
            .expect("canonical SMILES should write"),
        smiles_api::write_canonical_with_options(&second, CanonicalSmilesWriteOptions)
            .expect("canonical SMILES should write")
    );
}

#[test]
fn canonical_smiles_ignores_stereo_for_non_isomeric_output() {
    let mut molecule = read_smiles("N[C@H](O)C").expect("chiral SMILES parses");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("chiral molecule sanitizes");

    assert!(
        smiles_api::write_with_options(&molecule, SmilesWriteOptions)
            .expect_err("ordinary writer should reject lossy atom stereo")
            .message
            .contains("atom stereochemistry")
    );

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("non-isomeric canonical SMILES should ignore atom stereo");
    let reparsed = read_smiles(&written).expect("canonical output should parse");

    assert!(!written.contains('['), "{written}");
    assert_eq!(reparsed.graph().stereo_elements().count(), 0);
    assert_eq!(reparsed.graph().atom_count(), molecule.graph().atom_count());
    assert_eq!(reparsed.graph().bond_count(), molecule.graph().bond_count());

    let isotope = read_smiles("[11CH3]OC").expect("isotope parses");
    assert_eq!(
        smiles_api::write_canonical_with_options(&isotope, CanonicalSmilesWriteOptions)
            .expect("non-isomeric canonical SMILES should ignore isotope labels"),
        "COC"
    );

    let mut aromatic_isotope = read_smiles("C1=CC=[14CH]C=C1").expect("aromatic isotope parses");
    perception_api::sanitize_with_options(&mut aromatic_isotope, SanitizeOptions::default())
        .expect("aromatic isotope sanitizes");
    assert_eq!(
        smiles_api::write_canonical_with_options(&aromatic_isotope, CanonicalSmilesWriteOptions,)
            .expect("aromatic isotope canonicalizes"),
        "c1ccccc1"
    );

    let mut explicit_hydrogens =
        read_smiles("[H]C([3H])(F)Cl").expect("explicit hydrogen isotopologue parses");
    perception_api::sanitize_with_options(&mut explicit_hydrogens, SanitizeOptions::default())
        .expect("explicit hydrogen isotopologue sanitizes");
    let written =
        smiles_api::write_canonical_with_options(&explicit_hydrogens, CanonicalSmilesWriteOptions)
            .expect("explicit hydrogen isotopologue canonicalizes");
    assert_eq!(written.matches("[H]").count(), 1, "{written}");
    let reparsed = read_smiles(&written).expect("normalized explicit hydrogen output reparses");
    assert_eq!(reparsed.graph().atom_count(), 4, "{written}");
}

#[test]
fn canonical_smiles_round_trips_supported_branch_and_ring_graphs() {
    for input in ["CC(=O)O", "C1CCCCC1", "c1ccccc1"] {
        let mut molecule =
            read_smiles(input).unwrap_or_else(|_| panic!("SMILES should parse: {input}"));
        perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
            .unwrap_or_else(|_| panic!("SMILES should sanitize: {input}"));
        let written =
            smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
                .unwrap_or_else(|_| panic!("canonical SMILES should write: {input}"));
        let reparsed = read_smiles(&written)
            .unwrap_or_else(|_| panic!("canonical output should parse: {written}"));

        assert_eq!(reparsed.graph().atom_count(), molecule.graph().atom_count());
        assert_eq!(reparsed.graph().bond_count(), molecule.graph().bond_count());
    }
}

#[test]
fn canonical_smiles_prefers_clean_simple_ring_closure() {
    let molecule = read_smiles("C1=CC=CC=C1").expect("benzene parses");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");

    assert_eq!(written, "C1=CC=CC=C1");
}

#[test]
fn canonical_smiles_converges_after_aromaticity_perception() {
    let mut aromatic = read_smiles("c1ccccc1").expect("aromatic benzene parses");
    let mut kekule = read_smiles("C1=CC=CC=C1").expect("Kekule benzene parses");
    perception_api::sanitize_with_options(&mut aromatic, SanitizeOptions::default())
        .expect("aromatic benzene sanitizes");
    perception_api::sanitize_with_options(&mut kekule, SanitizeOptions::default())
        .expect("Kekule benzene sanitizes");

    let aromatic_written =
        smiles_api::write_canonical_with_options(&aromatic, CanonicalSmilesWriteOptions)
            .expect("aromatic benzene canonicalizes");
    let kekule_written =
        smiles_api::write_canonical_with_options(&kekule, CanonicalSmilesWriteOptions)
            .expect("perceived Kekule benzene canonicalizes");

    assert_eq!(aromatic_written, kekule_written);
    assert_eq!(aromatic_written, "c1ccccc1");
}

#[test]
fn canonical_smiles_preserves_aromatic_high_order_bonds() {
    let mut molecule = read_smiles("C1=CC#CC=C1").expect("cyclohexyne parses");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("cyclohexyne sanitizes");

    let written = smiles_api::write_canonical(&molecule).expect("cyclohexyne canonicalizes");
    assert!(written.contains('#'), "{written}");
    let mut reparsed = read_smiles(&written).expect("cyclohexyne canonical output reparses");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .expect("cyclohexyne canonical output sanitizes");
    assert!(reparsed
        .graph()
        .bonds()
        .any(|(_, bond)| bond.order == BondOrder::Triple && bond.aromatic));
}

#[test]
fn canonical_smiles_implementation_avoids_sanitizer_feedback() {
    let source = include_str!("../io/smiles.rs");

    assert!(!source.contains("canonical_smiles_candidate_sanitize_rank"));
    assert!(!source.contains("canonical_smiles_semantic_signature"));
    assert!(!source.contains("KekuleWhenStored"));
}

#[test]
fn aromatic_smiles_omitted_bonds_sanitize_with_expected_hydrogens() {
    let mut benzene = read_smiles("c1ccccc1").expect("benzene should parse");
    assert!(benzene
        .graph()
        .bonds()
        .all(|(_, bond)| bond.order == BondOrder::Aromatic));
    perception_api::sanitize_with_options(&mut benzene, SanitizeOptions::default())
        .expect("benzene should sanitize");
    for (_, atom) in benzene.graph().atoms() {
        assert_eq!(atom.implicit_hydrogens, Some(1));
        assert!(atom.aromatic);
    }

    let mut pyridine = read_smiles("n1ccccc1").expect("pyridine should parse");
    perception_api::sanitize_with_options(&mut pyridine, SanitizeOptions::default())
        .expect("pyridine should sanitize");
    let nitrogen = pyridine.graph().atom(AtomId::new(0)).expect("nitrogen");
    assert_eq!(nitrogen.implicit_hydrogens, Some(0));
    for atom_id in 1..6 {
        let atom = pyridine.graph().atom(AtomId::new(atom_id)).expect("carbon");
        assert_eq!(atom.implicit_hydrogens, Some(1));
    }

    let mut pyridinium = read_smiles("[n+]1ccccc1").expect("pyridinium should parse");
    perception_api::sanitize_with_options(&mut pyridinium, SanitizeOptions::default())
        .expect("pyridinium should sanitize");
    let nitrogen = pyridinium.graph().atom(AtomId::new(0)).expect("nitrogen");
    assert!(!nitrogen.aromatic);
    assert_eq!(nitrogen.formal_charge, 1);
    assert_eq!(nitrogen.radical, Some(AtomRadical::Doublet));
    assert_eq!(nitrogen.implicit_hydrogens, Some(0));
    assert!(pyridinium.graph().bonds().all(|(_, bond)| !bond.aromatic));
    assert_eq!(
        pyridinium
            .graph()
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        3
    );

    for smiles in [
        "[nH]1cccc1",
        "c1ccoc1",
        "c1ccsc1",
        "c1ccc2ccccc2c1",
        "Cc1ccccc1",
        "c1ccccc1.CC",
        "C%10CCCCC%10",
    ] {
        let mut molecule = read_smiles(smiles)
            .unwrap_or_else(|_| panic!("supported aromatic SMILES should parse: {smiles}"));
        perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
            .unwrap_or_else(|_| panic!("supported aromatic SMILES should sanitize: {smiles}"));
        let written = smiles_api::write_with_options(&molecule, SmilesWriteOptions)
            .unwrap_or_else(|_| panic!("supported aromatic SMILES should write: {smiles}"));
        read_smiles(&written).unwrap_or_else(|_| panic!("writer output should parse: {written}"));
    }
}

#[test]
fn invalid_lowercase_aromatic_ring_returns_structured_error() {
    for smiles in ["c1cccc1", "c1ccccc1.c1cccc1"] {
        let mut molecule = read_smiles(smiles).expect("raw syntax should parse");
        let error =
            perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
                .expect_err("invalid aromatic ring should fail sanitization");
        assert!(matches!(
            error,
            SanitizeError::Aromaticity(AromaticityError::InvalidAromaticRepresentation(_))
        ));
    }
}

#[test]
fn thiocarbonyl_chalcogen_ring_sanitizes_aromatic_like_rdkit() {
    let mut molecule = read_smiles("CCN(CC)C1=NC(=S)N(C(=S)S1)C(=S)N(CC)CC")
        .expect("thiocarbonyl heterocycle should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("thiocarbonyl heterocycle should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    let aromatic_bonds = molecule
        .graph()
        .bonds()
        .filter(|(_, bond)| bond.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 6);
    assert_eq!(aromatic_bonds, 6);

    let written = smiles_api::write_with_options(&molecule, SmilesWriteOptions)
        .expect("sanitized thiocarbonyl heterocycle should write");
    let reparsed = read_smiles(&written).expect("writer output should parse");
    assert_eq!(reparsed.graph().atom_count(), molecule.graph().atom_count());
    assert_eq!(reparsed.graph().bond_count(), molecule.graph().bond_count());

    let canonical =
        smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
            .expect("sanitized thiocarbonyl heterocycle should canonicalize");
    let mut canonical_reparsed = read_smiles(&canonical).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut canonical_reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {canonical}: {error:?}"));
    assert_eq!(
        canonical_reparsed
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        6,
        "{canonical}"
    );
}

#[test]
fn fused_chalcogen_bridge_does_not_over_aromatize_hetero_bridge() {
    let mut molecule = read_smiles("CSC1=CC2=C(C=C1)SC3=CC=CC=C3N2")
        .expect("phenothiazine-like heterocycle should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("phenothiazine-like heterocycle should sanitize");

    let sulfur_bridge = molecule
        .graph()
        .atom(AtomId::new(8))
        .expect("bridge sulfur");
    let nitrogen_bridge = molecule
        .graph()
        .atom(AtomId::new(15))
        .expect("bridge nitrogen");
    assert!(!sulfur_bridge.aromatic);
    assert!(!nitrogen_bridge.aromatic);

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 12);
}

#[test]
fn bracket_carbon_suppresses_implicit_hydrogens() {
    let mut molecule = read_smiles("C1=CC=C2C(=C1)[CH]C3=CC=CC=C32")
        .expect("bracket carbon fused aromatic should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("bracket carbon fused aromatic should sanitize");

    let bracket_carbon = molecule
        .graph()
        .atom(AtomId::new(6))
        .expect("bracket carbon");
    assert!(bracket_carbon.no_implicit_hydrogens);
    assert_eq!(bracket_carbon.explicit_hydrogens, 1);
    assert_eq!(bracket_carbon.implicit_hydrogens, Some(0));
}

#[test]
fn fused_polycycle_aromatic_core_extends_to_fused_edge() {
    let mut molecule =
        read_smiles("C1CCC2=C(C1)C3=C(C=CC4=C3C5=C(C=C4)C=CC(=C25)[N+](=O)[O-])[N+](=O)[O-]")
            .expect("fused polycycle should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused polycycle should sanitize");

    for atom_id in 0..3 {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("saturated atom")
                .aromatic,
            "saturated ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("fused core atom")
                .aromatic,
            "fused core atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn rdkit_source_comment_fused_system_matches_reference_counts() {
    let mut molecule = read_smiles("O=C3C2=CC1=CC=COC1=CC2=CC=C3")
        .expect("RDKit source fused example should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("RDKit source fused example should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    let aromatic_bonds = molecule
        .graph()
        .bonds()
        .filter(|(_, bond)| bond.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 14);
    assert_eq!(aromatic_bonds, 16);
}

#[test]
fn ring_atom_with_multiple_pi_bonds_is_not_aromatic_candidate() {
    let mut molecule = read_smiles("C1=C=NC=N1").expect("multiple-pi-bond ring should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("multiple-pi-bond ring should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    let aromatic_bonds = molecule
        .graph()
        .bonds()
        .filter(|(_, bond)| bond.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 0);
    assert_eq!(aromatic_bonds, 0);
}

#[test]
fn fused_aromatic_component_preserves_explicit_single_bond() {
    let mut molecule = read_smiles(
        "[H]c1c([H])c([H])c2c3c([H])c([H])n(C([H])([H])[H])c(C([H])([H])[H])c-3nc2c1[H]",
    )
    .expect("fused aromatic system should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused aromatic system should sanitize");

    let explicit_single_between_aromatic_atoms = molecule
        .graph()
        .bonds()
        .filter(|(_, bond)| {
            matches!(bond.order, BondOrder::Single)
                && molecule
                    .graph()
                    .atom(bond.a())
                    .is_ok_and(|atom| atom.aromatic)
                && molecule
                    .graph()
                    .atom(bond.b())
                    .is_ok_and(|atom| atom.aromatic)
        })
        .collect::<Vec<_>>();
    assert_eq!(explicit_single_between_aromatic_atoms.len(), 1);
    assert!(!explicit_single_between_aromatic_atoms[0].1.aromatic);

    let written = smiles_api::write_with_options(&molecule, SmilesWriteOptions)
        .expect("fused aromatic system should write");
    assert!(written.contains('-'));
    let reparsed = read_smiles(&written).expect("writer output should parse");
    assert_eq!(
        reparsed
            .graph()
            .bonds()
            .filter(|(_, bond)| {
                matches!(bond.order, BondOrder::Single)
                    && reparsed
                        .graph()
                        .atom(bond.a())
                        .is_ok_and(|atom| atom.aromatic)
                    && reparsed
                        .graph()
                        .atom(bond.b())
                        .is_ok_and(|atom| atom.aromatic)
            })
            .count(),
        1
    );
}

#[test]
fn fused_subset_marks_perimeter_without_aromatizing_internal_shared_bond() {
    let mut molecule = read_smiles("O=C(NC1=CC=CC=C1)N1CCCC(C(=O)N2CCN(C3=C4C=CN=C4NC=N3)CC2)C1")
        .expect("fused Enamine regression should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused Enamine regression should sanitize");

    for atom_id in [20, 21, 22, 23, 24, 25, 26, 27, 28] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("fused aromatic atom")
                .aromatic,
            "fused atom {atom_id} should be aromatic"
        );
    }
    assert!(
        molecule
            .graph()
            .bond(BondId::new(30))
            .expect("accepted fused perimeter bond")
            .aromatic,
        "accepted fused-subset perimeter bond should be aromatic"
    );
    assert!(
        !molecule
            .graph()
            .bond(BondId::new(26))
            .expect("internal fused shared bond")
            .aromatic,
        "internal fused shared bond should remain aliphatic"
    );
}

#[test]
fn fused_sdf_five_electron_neighbor_shares_aromatic_internal_bond() {
    let input = r#"
     RDKit          2D

 27 29  0  0  0  0  0  0  0  0999 V2000
    4.5000   -5.1962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -3.8971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500   -3.8971    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -2.5981    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.2010   -3.3481    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7990   -1.8481    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7500    3.8971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    2.5981    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    3.8971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
  2 11  1  0
 11 12  1  0
 12 13  2  0
 12 14  2  0
 12 15  1  0
 15 16  2  0
 16 17  1  0
 17 18  2  0
 18 19  1  0
 19 20  2  0
 19 21  1  0
 21 22  2  0
 21 23  1  0
 23 24  1  0
 24 25  2  0
 24 26  1  0
 26 27  1  0
  5 10  1  0
 15 20  1  0
 26 18  1  0
M  END
$$$$
"#;
    let mut molecule = read_sdf_molecules(input)
        .expect("regression SDF parses")
        .into_iter()
        .next()
        .expect("one SDF molecule");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused SDF regression sanitizes");

    assert!(
        molecule
            .graph()
            .bond(BondId::new(17))
            .expect("shared fused bond")
            .aromatic,
        "candidate five-electron fused neighbor should not suppress the accepted shared bond"
    );
}

#[test]
fn fused_dione_partner_shares_aromatic_internal_bond() {
    let input = r#"10250
  -OEChem-06192605442D

 16 17  0     0  0  0  0  0  0999 V2000
    3.7321    1.8100    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000   -1.1900    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321   -1.1900    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660    0.3100    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.4920    0.8447    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.4920   -1.2247    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981    0.3100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981   -0.6900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321    0.8100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -0.6900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3981    0.3308    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3981   -0.7108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321   -1.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3291    0.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.9338    0.6429    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.9338   -1.0229    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  9  2  0  0  0  0
  2 10  2  0  0  0  0
  3  8  1  0  0  0  0
  3 10  1  0  0  0  0
  3 13  1  0  0  0  0
  4  9  1  0  0  0  0
  4 10  1  0  0  0  0
  4 14  1  0  0  0  0
  5  7  2  0  0  0  0
  5 11  1  0  0  0  0
  6  8  2  0  0  0  0
  6 12  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  1  0  0  0  0
 11 12  2  0  0  0  0
 11 15  1  0  0  0  0
 12 16  1  0  0  0  0
M  END
$$$$
"#;
    let mut molecule = read_sdf_molecules(input)
        .expect("fused dione SDF parses")
        .into_iter()
        .next()
        .expect("one SDF molecule");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused dione regression sanitizes");

    assert!(
        molecule
            .graph()
            .bond(BondId::new(12))
            .expect("shared dione fused bond")
            .aromatic,
        "accepted lone-pair dione fused partner should not suppress the shared bond"
    );
    assert_eq!(
        molecule
            .graph()
            .bonds()
            .filter(|(_, bond)| bond.aromatic)
            .count(),
        11
    );
}

#[test]
fn fused_chalcogen_subset_with_exocyclic_pi_links_becomes_aromatic() {
    let mut molecule = read_smiles("CC1=C2OC3=C(C)C=CC(C(=O)NC4C(=O)NC(C(C)C)C(=O)N5CCCC5C(=O)N(C)CC(=O)N(C)C(C(C)C)C(=O)OC4C)=C3N=C2C(C(=O)NC2C(=O)NC(C(C)C)C(=O)N3CCCC3C(=O)N(C)CC(=O)N(C)C(C(C)C)C(=O)OC2C)=C(N)C1=O")
    .expect("PubChem fused chalcogen regression should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("PubChem fused chalcogen regression should sanitize");

    assert!(
        molecule
            .graph()
            .atom(AtomId::new(3))
            .expect("fused bridge oxygen")
            .aromatic,
        "fused bridge oxygen should be aromatic"
    );
    for bond_id in [2, 3] {
        assert!(
            molecule
                .graph()
                .bond(BondId::new(bond_id))
                .expect("oxygen fused-subset perimeter bond")
                .aromatic,
            "oxygen perimeter bond {bond_id} should be aromatic"
        );
    }
    assert!(
        !molecule
            .graph()
            .atom(AtomId::new(10))
            .expect("carbonyl center")
            .aromatic,
        "adjacent carbonyl center should stay aliphatic"
    );
}

#[test]
fn fused_simple_aromatic_member_rings_can_share_aromatic_single_bond() {
    let mut molecule = read_smiles("CC(CCC1=CC=CC=C1)NS(=O)(=O)C1=CC2=C(N=C1)N(C)C(=O)NC2=O")
        .expect("fused pyrimidinedione regression should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused pyrimidinedione regression should sanitize");

    assert!(
        molecule
            .graph()
            .bond(BondId::new(17))
            .expect("shared simple-ring bond")
            .aromatic,
        "shared bond between simple aromatic member rings should be aromatic"
    );
}

#[test]
fn fused_quinone_cn_core_excludes_carbonyl_centers() {
    let mut molecule = read_smiles("C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C4=C(C=C3)C(=O)C5=CC=CC=C5N4")
        .expect("fused quinone should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused quinone should sanitize");

    for atom_id in [6, 10] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("carbonyl center")
                .aromatic,
            "carbonyl center {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [
        0, 1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24,
    ] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("fused aromatic atom")
                .aromatic,
            "fused aromatic atom {atom_id} should be aromatic"
        );
    }
    assert_eq!(
        molecule
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        20
    );
}

#[test]
fn canonical_fused_quinone_cn_core_round_trip_matches_aromatic_shape() {
    let mut molecule = read_smiles("C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C4=C(C=C3)C(=O)C5=CC=CC=C5N4")
        .expect("fused quinone should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused quinone should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical fused quinone should write");
    let mut reparsed = read_smiles(&written).expect("canonical fused quinone output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    assert_eq!(reparsed.graph().atom_count(), molecule.graph().atom_count());
    assert_eq!(reparsed.graph().bond_count(), molecule.graph().bond_count());
    let original_aromatic_carbonyl_centers = aromatic_carbonyl_center_count(molecule.graph());
    let reparsed_aromatic_carbonyl_centers = aromatic_carbonyl_center_count(reparsed.graph());
    assert_eq!(
        reparsed_aromatic_carbonyl_centers, original_aromatic_carbonyl_centers,
        "{written}"
    );
    let aromatic_n_h_count = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| {
            atom.element.symbol() == "N"
                && atom.aromatic
                && atom
                    .explicit_hydrogens
                    .saturating_add(atom.implicit_hydrogens.unwrap_or(0))
                    == 1
        })
        .count();
    assert_eq!(aromatic_n_h_count, 1, "{written}");
}

#[test]
fn canonical_aromatic_carbonyl_component_uses_representable_kekule_form() {
    let mut molecule =
        read_smiles("CN(C)CCOC(=O)CCNC1=CC=CC=CC1=O").expect("aromatic carbonyl ring should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("aromatic carbonyl ring should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical aromatic carbonyl should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    assert_eq!(
        local_atom_neighbor_signatures(molecule.graph()),
        local_atom_neighbor_signatures(reparsed.graph()),
        "{written}"
    );
}

#[test]
fn canonical_charged_aromatic_carbon_component_uses_representable_kekule_form() {
    let mut molecule = read_smiles("C1CCOC1.[CH-]1[C-]=[C-][C-]=[C-]1.Cl[Cr]Cl")
        .expect("cyclopentadienyl salt should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("cyclopentadienyl salt should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical cyclopentadienyl salt should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    assert_eq!(
        local_atom_neighbor_signatures_ignoring_halogen_no_implicit(molecule.graph()),
        local_atom_neighbor_signatures_ignoring_halogen_no_implicit(reparsed.graph()),
        "{written}"
    );
    assert!(written.contains("[Cl][Cr][Cl]"), "{written}");
}

fn aromatic_carbonyl_center_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(atom_id, atom)| {
            atom.element.symbol() == "C"
                && atom.aromatic
                && mol
                    .incident_bonds(*atom_id)
                    .expect("atom should be live")
                    .any(|(_, bond)| {
                        bond.order == BondOrder::Double
                            && mol
                                .atom(bond.other_atom(*atom_id))
                                .is_ok_and(|neighbor| neighbor.element.symbol() == "O")
                    })
        })
        .count()
}

#[test]
fn indole_quinone_keeps_carbonyl_ring_atoms_aliphatic() {
    let mut molecule =
        read_smiles("C1=CC=C(C=C1)C2=CC3=C(N2)C(=O)C=CC3=O").expect("indole quinone should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("indole quinone should sanitize");

    for atom_id in [11, 13, 14, 15] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("quinone ring atom")
                .aromatic,
            "quinone ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("fused aromatic atom")
                .aromatic,
            "fused aromatic atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_imine_and_pyrimidinedione_aromaticity_matches_reference_shape() {
    let mut molecule = read_smiles("CC1=NC2=CC=CC=C2C1=CC3=C(NC(=O)NC3=O)O")
        .expect("fused imine and pyrimidinedione should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused imine and pyrimidinedione should sanitize");

    for atom_id in [1, 2, 9] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("imine-ring atom")
                .aromatic,
            "imine-ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 16, 17] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("aromatic fused atom")
                .aromatic,
            "fused aromatic atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn exocyclic_iminium_sulfur_ring_remains_aromatic() {
    let mut molecule = read_smiles("CN(C1=NC(=[N+](C)C)SS1)C(=S)SC")
        .expect("exocyclic iminium sulfur ring should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("exocyclic iminium sulfur ring should sanitize");

    for atom_id in [2, 3, 4, 8, 9] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("heteroaromatic ring atom")
                .aromatic,
            "heteroaromatic ring atom {atom_id} should be aromatic"
        );
    }
    let exocyclic_iminium = molecule.graph().atom(AtomId::new(5)).expect("iminium N");
    assert!(!exocyclic_iminium.aromatic);
    assert_eq!(exocyclic_iminium.formal_charge, 1);
}

#[test]
fn fused_exocyclic_imine_sulfur_ring_remains_aromatic() {
    let mut molecule = read_smiles("CCCCCCCCCCCCCCCCS(=O)(=O)N(C(=O)OCC)N=C1N(C2=CC=CC=C2S1)C")
        .expect("fused exocyclic imine sulfur ring should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused exocyclic imine sulfur ring should sanitize");

    for atom_id in [26, 27, 28, 29, 30, 31, 32, 33, 34] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("fused benzothiazine atom")
                .aromatic,
            "fused benzothiazine atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_sulfoxide_ring_does_not_follow_benzene_aromaticity() {
    let mut molecule =
        read_smiles("CCCCCN1SC2=CC=CC=C2S1=O").expect("fused sulfoxide ring should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused sulfoxide ring should sanitize");

    for atom_id in [5, 6, 13] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("sulfoxide ring atom")
                .aromatic,
            "sulfoxide ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [7, 8, 9, 10, 11, 12] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("benzene ring atom")
                .aromatic,
            "benzene ring atom {atom_id} should stay aromatic"
        );
    }
}

#[test]
fn neutral_exocyclic_alkene_sulfur_ring_stays_aliphatic() {
    let mut molecule = read_smiles(
        "C1=CC=C2C(=C1)C=CC3=C2[N+](=C(S3)C=C4N(C5=CC=CC=C5S4)CCCS(=O)(=O)O)CCCS(=O)(=O)O",
    )
    .expect("mixed sulfur fused system should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("mixed sulfur fused system should sanitize");

    for atom_id in [14, 15, 22] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("neutral sulfur ring atom")
                .aromatic,
            "neutral sulfur ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [10, 12] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("cationic sulfur ring atom")
                .aromatic,
            "cationic sulfur ring atom {atom_id} should be aromatic"
        );
    }
    let exocyclic_alkene_ring_carbons = molecule
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "C"
                && !atom.aromatic
                && molecule.graph().incident_bonds(*id).is_ok_and(|bonds| {
                    let bonds = bonds.collect::<Vec<_>>();
                    bonds.len() == 3
                        && bonds.iter().any(|(_, bond)| {
                            matches!(bond.order, BondOrder::Double)
                                && molecule
                                    .graph()
                                    .atom(bond.other_atom(*id))
                                    .is_ok_and(|other| other.element.symbol() == "C")
                        })
                        && bonds.iter().any(|(_, bond)| {
                            molecule
                                .graph()
                                .atom(bond.other_atom(*id))
                                .is_ok_and(|other| other.element.symbol() == "N")
                        })
                        && bonds.iter().any(|(_, bond)| {
                            molecule
                                .graph()
                                .atom(bond.other_atom(*id))
                                .is_ok_and(|other| other.element.symbol() == "S")
                        })
                })
        })
        .count();
    assert_eq!(exocyclic_alkene_ring_carbons, 1);
    let aliphatic_neutral_ring_nitrogens = molecule
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "N"
                && atom.formal_charge == 0
                && !atom.aromatic
                && molecule.graph().incident_bonds(*id).is_ok_and(|bonds| {
                    bonds.into_iter().any(|(_, bond)| {
                        molecule
                            .graph()
                            .atom(bond.other_atom(*id))
                            .is_ok_and(|other| {
                                other.element.symbol() == "C"
                                    && !other.aromatic
                                    && molecule
                                        .graph()
                                        .incident_bonds(bond.other_atom(*id))
                                        .is_ok_and(|carbon_bonds| {
                                            carbon_bonds.into_iter().any(|(_, carbon_bond)| {
                                                molecule
                                                    .graph()
                                                    .atom(
                                                        carbon_bond
                                                            .other_atom(bond.other_atom(*id)),
                                                    )
                                                    .is_ok_and(|neighbor| {
                                                        neighbor.element.symbol() == "S"
                                                    })
                                            })
                                        })
                            })
                    })
                })
        })
        .count();
    assert_eq!(aliphatic_neutral_ring_nitrogens, 1);
}

#[test]
fn fused_seven_membered_ether_ring_stays_aliphatic() {
    let mut molecule = read_smiles("CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C6(C4)C(=O)C7=C8N6CCC9=C8C(=C(C=C9)OC)OC1=C7C=CC(=C1O)OC")
        .expect("fused ether polycycle should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused ether polycycle should sanitize");

    assert!(molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "O")
        .all(|(_, atom)| !atom.aromatic));
}

#[test]
fn charged_bracket_halogen_and_bismuth_salt_sanitizes() {
    let mut molecule =
        read_smiles("C1CC2CCC[N-]C2C(C1)[OH2+].C1C=CC2=CC=CC(C2=N1)[OH2+].[ClH2+].Cl.[Bi+3]")
            .expect("charged bracket salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("charged bracket salt should sanitize");

    let protonated_chlorine = molecule.graph().atom(AtomId::new(22)).expect("chlorine");
    assert_eq!(protonated_chlorine.element.symbol(), "Cl");
    assert_eq!(protonated_chlorine.formal_charge, 1);
    assert_eq!(protonated_chlorine.explicit_hydrogens, 2);
    assert_eq!(protonated_chlorine.implicit_hydrogens, Some(0));

    let bismuth = molecule.graph().atom(AtomId::new(24)).expect("bismuth");
    assert_eq!(bismuth.element.symbol(), "Bi");
    assert_eq!(bismuth.formal_charge, 3);
    assert_eq!(bismuth.implicit_hydrogens, Some(0));
}

#[test]
fn oxide_dianion_transition_metal_salt_sanitizes() {
    let mut molecule = read_smiles("[O-2].[O-2].[O-2].[Cr+3].[Fe+3]")
        .expect("oxide transition-metal salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("oxide transition-metal salt should sanitize");

    for atom_id in [0, 1, 2] {
        let oxygen = molecule.graph().atom(AtomId::new(atom_id)).expect("oxide");
        assert_eq!(oxygen.element.symbol(), "O");
        assert_eq!(oxygen.formal_charge, -2);
        assert_eq!(oxygen.implicit_hydrogens, Some(0));
    }
}

#[test]
fn hydroxide_niobium_v_salt_sanitizes() {
    let mut molecule = read_smiles("[OH-].[Nb+5]").expect("niobium hydroxide salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("niobium hydroxide salt should sanitize");

    let niobium = molecule.graph().atom(AtomId::new(1)).expect("niobium");
    assert_eq!(niobium.element.symbol(), "Nb");
    assert_eq!(niobium.formal_charge, 5);
    assert_eq!(niobium.implicit_hydrogens, Some(0));
}

#[test]
fn formate_indium_salt_sanitizes() {
    let mut molecule = read_smiles("C(=O)[O-].C(=O)[O-].C(=O)[O-].[In+3]")
        .expect("indium formate salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("indium formate salt should sanitize");

    let indium = molecule.graph().atom(AtomId::new(9)).expect("indium");
    assert_eq!(indium.element.symbol(), "In");
    assert_eq!(indium.formal_charge, 3);
    assert_eq!(indium.implicit_hydrogens, Some(0));
}

#[test]
fn periodate_cleanup_sanitizes_iodine_plus_three() {
    let mut molecule = read_smiles("[O-]I(=O)(=O)=O").expect("periodate should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("periodate should sanitize");

    let iodine = molecule.graph().atom(AtomId::new(1)).expect("iodine");
    assert_eq!(iodine.element.symbol(), "I");
    assert_eq!(iodine.formal_charge, 3);
    assert_eq!(iodine.implicit_hydrogens, Some(0));
}

#[test]
fn sodium_chlorate_sanitizes_without_aromaticity() {
    let mut molecule = read_smiles("[O-]Cl(=O)=O.[Na+]").expect("sodium chlorate should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("sodium chlorate should sanitize");
    assert!(
        molecule.graph().atoms().all(|(_, atom)| !atom.aromatic)
            && molecule.graph().bonds().all(|(_, bond)| !bond.aromatic)
    );
}

#[test]
fn oxohalogen_cleanup_distinguishes_oxyacids_from_carbon_substituents() {
    let mut iodous_acid = read_smiles("OI=O").expect("iodous acid should parse");
    perception_api::sanitize_with_options(&mut iodous_acid, SanitizeOptions::default())
        .expect("iodous acid should sanitize");
    let iodine = iodous_acid
        .graph()
        .atoms()
        .find(|(_, atom)| atom.element.symbol() == "I")
        .expect("iodine")
        .1;
    assert_eq!(iodine.formal_charge, 1);
    assert!(iodous_acid
        .graph()
        .atoms()
        .any(|(_, atom)| atom.element.symbol() == "O" && atom.formal_charge == -1));

    let mut iodyl_methane = read_smiles("CI(=O)=O").expect("iodyl methane should parse");
    perception_api::sanitize_with_options(&mut iodyl_methane, SanitizeOptions::default())
        .expect("iodyl methane should sanitize");
    let iodine = iodyl_methane
        .graph()
        .atoms()
        .find(|(_, atom)| atom.element.symbol() == "I")
        .expect("iodine")
        .1;
    assert_eq!(iodine.formal_charge, 0);
    assert!(iodyl_methane
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "O")
        .all(|(_, atom)| atom.formal_charge == 0));

    let mut cyclic_iodane_fragment =
        read_smiles("COI(=O)(N)C").expect("iodane with a bridging oxygen should parse");
    perception_api::sanitize_with_options(&mut cyclic_iodane_fragment, SanitizeOptions::default())
        .expect("neutral lambda-five iodane should sanitize");
    let iodine = cyclic_iodane_fragment
        .graph()
        .atoms()
        .find(|(_, atom)| atom.element.symbol() == "I")
        .expect("iodine")
        .1;
    assert_eq!(iodine.formal_charge, 0);
    assert!(cyclic_iodane_fragment
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "O")
        .all(|(_, atom)| atom.formal_charge == 0));
}

#[test]
fn uranyl_beta_diketonate_salt_sanitizes() {
    let mut molecule = read_smiles(
        "C1=CC=C(C=C1)C(=O)[CH-]C(=O)C2=CC=CC=C2.C1=CC=C(C=C1)C(=O)[CH-]C(=O)C2=CC=CC=C2.O=[U+2]=O",
    )
    .expect("uranyl salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("uranyl salt should sanitize");

    let uranium = molecule.graph().atom(AtomId::new(35)).expect("uranium");
    assert_eq!(uranium.element.symbol(), "U");
    assert_eq!(uranium.formal_charge, 2);
    assert_eq!(uranium.implicit_hydrogens, Some(0));
}

#[test]
fn cyclopentadienyl_anion_sanitizes_aromatic() {
    let mut molecule = read_smiles("C1CCOC1.[CH-]1[C-]=[C-][C-]=[C-]1.Cl[Cr]Cl")
        .expect("cyclopentadienyl chromium salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("cyclopentadienyl chromium salt should sanitize");

    for atom_id in 5..=9 {
        let atom = molecule
            .graph()
            .atom(AtomId::new(atom_id))
            .expect("cyclopentadienyl atom");
        assert_eq!(atom.element.symbol(), "C");
        assert_eq!(atom.formal_charge, -1);
        assert!(
            atom.aromatic,
            "anion ring atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_quinone_ring_does_not_follow_benzene_aromaticity() {
    let mut molecule =
        read_smiles("CC(C)(C)NN=C(C1C=CCS1(=O)=O)C(=O)NC2=C(C(=O)C3=CC=CC=C3C2=O)Cl")
            .expect("fused quinone sulfone should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused quinone sulfone should sanitize");

    for atom_id in [17, 18, 19, 27] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("quinone ring atom")
                .aromatic,
            "quinone ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [21, 22, 23, 24, 25, 26] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("benzene ring atom")
                .aromatic,
            "benzene ring atom {atom_id} should stay aromatic"
        );
    }
}

#[test]
fn singly_carbonylated_fused_ring_stays_aromatic() {
    let mut molecule = read_smiles("CNCCN=C1C=CC2=C3C1=C(C4=C(C=CC(=O)C4=C3NN2CCNCCO)O)O.O.Cl.Cl")
        .expect("singly carbonylated fused ring should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("singly carbonylated fused ring should sanitize");

    for atom_id in [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("aromatic fused atom")
                .aromatic,
            "fused atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn saturated_fused_ring_does_not_follow_aromatic_core() {
    let mut molecule = read_smiles("C1CCC2=NC3=CC=CC=C3C(=C2C1)[NH2+]CCSCCCl.[Cl-]")
        .expect("saturated fused ring salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("saturated fused ring salt should sanitize");

    for atom_id in [0, 1, 2, 13] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("saturated fused atom")
                .aromatic,
            "saturated fused atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("aromatic core atom")
                .aromatic,
            "aromatic core atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn canonical_saturated_fused_ring_round_trip_stays_aliphatic() {
    let mut molecule = read_smiles("C1CCC2=NC3=CC=CC=C3C(=C2C1)[NH2+]CCSCCCl.[Cl-]")
        .expect("saturated fused ring salt should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("saturated fused ring salt should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("saturated fused ring salt should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let saturated_carbons = reparsed
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "C"
                && atom.implicit_hydrogens == Some(2)
                && reparsed
                    .graph()
                    .incident_bonds(*id)
                    .is_ok_and(|bonds| bonds.count() == 2)
        })
        .filter(|(_, atom)| !atom.aromatic)
        .count();
    assert!(
        saturated_carbons >= 4,
        "canonical output should keep saturated fused carbons aliphatic: {written}"
    );
}

#[test]
fn canonical_fused_chromanone_round_trip_keeps_lactone_ring_aliphatic() {
    let mut molecule =
        read_smiles("C1C(C(=O)C2=CC=CC=C2O1)C3=CC=CC=C3").expect("fused chromanone should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused chromanone should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused chromanone should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(
        aromatic_atoms, 12,
        "canonical output should keep only the phenyl and fused benzene rings aromatic: {written}"
    );
}

#[test]
fn conjugated_fused_benzopyrone_round_trip_keeps_lactone_ring_aromatic() {
    let mut molecule = read_smiles(
        "C1=CC=C(C=C1)C2=C(C(=O)C3=CC=CC=C3O2)OC(=O)C4=CC5=C(C=C4Cl)SC6=NC=CN6S5(=O)=O",
    )
    .expect("conjugated benzopyrone should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("conjugated benzopyrone should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("conjugated benzopyrone should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(
        aromatic_atoms, 27,
        "canonical output should preserve the conjugated benzopyrone aromatic system: {written}"
    );
}

#[test]
fn fused_fluorenone_round_trip_keeps_carbonyl_bridge_aliphatic() {
    let mut molecule = read_smiles("C1=CC=C2C(=C1)C3=C(C2=O)C=C(C=C3)[N+]#N.C(=O)(C(F)(F)F)O")
        .expect("fused fluorenone salt should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused fluorenone salt should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused fluorenone salt should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(
        aromatic_atoms, 12,
        "canonical output should keep the fluorenone carbonyl bridge aliphatic: {written}"
    );
}

#[test]
fn fused_saturated_carbonyl_bridge_round_trip_stays_aliphatic() {
    let mut molecule = read_smiles("CC1(CC(C(=O)C2=CC=CC=C21)(C(C3=CC=C(C=C3)[N+](=O)[O-])O)Cl)C")
        .expect("saturated carbonyl bridge should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("saturated carbonyl bridge should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("saturated carbonyl bridge should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(
        aromatic_atoms, 12,
        "canonical output should preserve only the two aromatic rings: {written}"
    );
    assert!(reparsed.graph().atoms().any(|(id, atom)| {
        atom.element.symbol() == "C"
            && !atom.aromatic
            && atom.implicit_hydrogens == Some(0)
            && reparsed.graph().incident_bonds(id).is_ok_and(|bonds| {
                let bonds = bonds.collect::<Vec<_>>();
                bonds.len() == 3
                    && bonds.iter().any(|(_, bond)| {
                        matches!(bond.order, BondOrder::Double)
                            && reparsed
                                .graph()
                                .atom(bond.other_atom(id))
                                .is_ok_and(|other| other.element.symbol() == "O")
                    })
            })
    }));
}

#[test]
fn fused_multi_quinone_bridge_round_trip_sanitizes() {
    let mut molecule = read_smiles("C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C4=C(C=C3)C(=O)C5=CC=CC=C5N4")
        .expect("fused multi-quinone should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused multi-quinone should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused multi-quinone should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let rewritten =
        smiles_api::write_canonical_with_options(&reparsed, CanonicalSmilesWriteOptions)
            .expect("reparsed fused multi-quinone should canonicalize");
    assert!(!rewritten.is_empty());
}

#[test]
fn canonical_tellurophene_round_trip_preserves_aromatic_chalcogen() {
    let mut molecule = read_smiles("C1=C[Te]C=C1").expect("tellurophene should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("tellurophene should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("tellurophene should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 5, "{written}");
    let tellurium = reparsed
        .graph()
        .atoms()
        .find_map(|(_, atom)| (atom.element.symbol() == "Te").then_some(atom))
        .expect("tellurium atom");
    assert!(tellurium.aromatic, "{written}");
    assert!(tellurium.no_implicit_hydrogens, "{written}");
}

#[test]
fn canonical_aryl_mercury_round_trip_preserves_no_implicit_aromatic_carbon() {
    let mut molecule =
        read_smiles("C1=CC=C(C(=C1)[N+](=O)[O-])[Hg]").expect("aryl mercury should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("aryl mercury should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("aryl mercury should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let mercury_bound_carbon = reparsed
        .graph()
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "C"
                && reparsed.graph().incident_bonds(id).is_ok_and(|bonds| {
                    bonds.into_iter().any(|(_, bond)| {
                        reparsed
                            .graph()
                            .atom(bond.other_atom(id))
                            .is_ok_and(|neighbor| neighbor.element.symbol() == "Hg")
                    })
                }))
            .then_some(atom)
        })
        .expect("aromatic carbon bound to mercury");
    assert!(mercury_bound_carbon.aromatic, "{written}");
    assert!(mercury_bound_carbon.no_implicit_hydrogens, "{written}");
    assert_eq!(
        mercury_bound_carbon.implicit_hydrogens,
        Some(0),
        "{written}"
    );
}

#[test]
fn fused_sulfonamide_tertiary_amine_round_trip_keeps_ring_nitrogen_aliphatic() {
    let mut molecule = read_smiles("COC1=CC2=C(C=C1)OC(=C2)S(=O)(=O)N3CC(C4=C3C=C(C=C4)N)CCl")
        .expect("fused sulfonamide tertiary amine should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused sulfonamide tertiary amine should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused sulfonamide tertiary amine should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aliphatic_ring_nitrogens = reparsed
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "N"
                && !atom.aromatic
                && reparsed.graph().incident_bonds(*id).is_ok_and(|bonds| {
                    let bonds = bonds.collect::<Vec<_>>();
                    bonds.len() == 3
                        && bonds.iter().any(|(_, bond)| {
                            reparsed
                                .graph()
                                .atom(bond.other_atom(*id))
                                .is_ok_and(|neighbor| neighbor.element.symbol() == "S")
                        })
                })
        })
        .count();
    assert_eq!(aliphatic_ring_nitrogens, 1, "{written}");
}

#[test]
fn cationic_fused_imide_round_trip_clears_carbonyl_ring_atoms() {
    let mut molecule = read_smiles(
        "CC1=C(C(=[N+]2N1C(=O)C(C2=O)C3=CC=CC=C3)C)C4=C(N5C(=O)C(C(=O)[N+]5=C4C)C6=CC=CC=C6)C",
    )
    .expect("cationic fused imide should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("cationic fused imide should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("cationic fused imide should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 22, "{written}");
    let aliphatic_carbonyl_ring_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "C"
                && !atom.aromatic
                && atom.implicit_hydrogens == Some(0)
                && reparsed.graph().incident_bonds(*id).is_ok_and(|bonds| {
                    let bonds = bonds.collect::<Vec<_>>();
                    bonds.len() == 3
                        && bonds.iter().any(|(_, bond)| {
                            matches!(bond.order, BondOrder::Double)
                                && reparsed
                                    .graph()
                                    .atom(bond.other_atom(*id))
                                    .is_ok_and(|other| other.element.symbol() == "O")
                        })
                })
        })
        .count();
    assert_eq!(aliphatic_carbonyl_ring_atoms, 4, "{written}");
}

#[test]
fn fused_quinone_ring_round_trip_sanitizes() {
    let mut molecule =
        read_smiles("C1=CC=C(C=C1)C2=CC3=C(N2)C(=O)C=CC3=O").expect("fused quinone should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused quinone should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused quinone should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let rewritten =
        smiles_api::write_canonical_with_options(&reparsed, CanonicalSmilesWriteOptions)
            .expect("reparsed fused quinone should canonicalize");
    assert!(!rewritten.is_empty());
}

#[test]
fn thiofuran_pyrimidinedione_canonical_round_trip_sanitizes() {
    let mut molecule = read_smiles("CC1=CN(C(=O)NC1=O)[C@H]2C=C(CS2)CO")
        .expect("thiofuran pyrimidinedione should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("thiofuran pyrimidinedione should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("thiofuran pyrimidinedione should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));
}

#[test]
fn fused_thiadiazolopyrimidinone_canonical_round_trip_preserves_aromatic_nitrogen_valence() {
    let mut molecule =
        read_smiles("C1=CC=C2C(=C1)C=CC(=C2C=CC3=NN=C4N(C3=O)N=C(S4)C5=CC(=CC=C5)[N+](=O)[O-])O")
            .expect("fused thiadiazolopyrimidinone should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused thiadiazolopyrimidinone should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused thiadiazolopyrimidinone should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));
}

#[test]
fn imine_fused_benzene_with_exocyclic_pyrimidinedione_keeps_imine_ring_aliphatic() {
    let mut molecule = read_smiles("CC1=NC2=CC=CC=C2C1=CC3=C(NC(=O)NC3=O)O")
        .expect("imine fused benzene should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("imine fused benzene should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("imine fused benzene should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    let aromatic_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 12, "{written}");
    let aliphatic_imine_nitrogens = reparsed
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "N"
                && !atom.aromatic
                && reparsed.graph().incident_bonds(*id).is_ok_and(|bonds| {
                    bonds
                        .into_iter()
                        .any(|(_, bond)| matches!(bond.order, BondOrder::Double))
                })
        })
        .count();
    assert_eq!(aliphatic_imine_nitrogens, 1, "{written}");
}

#[test]
fn fused_naphthalimide_canonical_round_trip_sanitizes() {
    let mut molecule =
        read_smiles("C1=CC(=CN=C1)CN2C(=O)C3=C(C2=O)C=C(C=C3)N(C4=CC=C(C=C4)Cl)C5=CC=C(C=C5)Cl")
            .expect("fused naphthalimide should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused naphthalimide should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused naphthalimide should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    let rewritten =
        smiles_api::write_canonical_with_options(&reparsed, CanonicalSmilesWriteOptions)
            .expect("reparsed fused naphthalimide should canonicalize");
    assert!(!rewritten.is_empty());
}

#[test]
fn partially_saturated_fused_amide_enone_ring_stays_aliphatic() {
    let mut molecule =
        read_smiles("CC1=CC=C(C=C1)C2=CC3=C(CCC(=C3)C(=O)NC4=CC=C(C=C4)C[N+]5(CCCCC5)C)C=C2")
            .expect("partially saturated fused amide enone should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("partially saturated fused amide enone should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("partially saturated fused amide enone should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    let aliphatic_enone_ring_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "C"
                && !atom.aromatic
                && reparsed.graph().incident_bonds(*id).is_ok_and(|bonds| {
                    let bonds = bonds.collect::<Vec<_>>();
                    bonds
                        .iter()
                        .any(|(_, bond)| matches!(bond.order, BondOrder::Double))
                        && bonds.iter().any(|(_, bond)| {
                            reparsed
                                .graph()
                                .atom(bond.other_atom(*id))
                                .is_ok_and(|other| other.element.symbol() == "C" && !other.aromatic)
                        })
                })
        })
        .count();
    assert!(
        aliphatic_enone_ring_atoms >= 2,
        "canonical output should keep the fused enone ring aliphatic: {written}"
    );
}

#[test]
fn fused_lactam_enone_canonical_round_trip_keeps_bridge_carbon_aliphatic() {
    let mut molecule = read_smiles("CN1CC[C@@]23[C@H]4[C@H]1CC5=C2C(=C(C=C5)OC)O[C@@H]3[C@]6(C4)C(=O)C7=C8N6CCC9=C8C(=C(C=C9)OC)OC1=C7C=CC(=C1O)OC")
    .expect("fused lactam enone should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused lactam enone should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused lactam enone should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    let aliphatic_lactam_enone_carbons = reparsed
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "C"
                && !atom.aromatic
                && reparsed.graph().incident_bonds(*id).is_ok_and(|bonds| {
                    let bonds = bonds.collect::<Vec<_>>();
                    bonds.iter().any(|(_, bond)| {
                        matches!(bond.order, BondOrder::Double)
                            && reparsed
                                .graph()
                                .atom(bond.other_atom(*id))
                                .is_ok_and(|other| other.element.symbol() == "C" && !other.aromatic)
                    }) && bonds.iter().any(|(_, bond)| {
                        matches!(bond.order, BondOrder::Single)
                            && reparsed
                                .graph()
                                .atom(bond.other_atom(*id))
                                .is_ok_and(|other| other.element.symbol() == "N")
                    })
                })
        })
        .count();
    assert!(
        aliphatic_lactam_enone_carbons >= 1,
        "canonical output should keep the lactam enone bridge aliphatic: {written}"
    );
    let aromatic_oxygens = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "O" && atom.aromatic)
        .count();
    assert_eq!(aromatic_oxygens, 0, "{written}");
}

#[test]
fn spiro_saturated_fused_hydrocarbon_bridge_stays_aliphatic() {
    let mut molecule = read_smiles("CC1=C2C=CC3=C(C2=CC=C1)CCC4(C3CCCC4)C")
        .expect("spiro saturated fused hydrocarbon should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("spiro saturated fused hydrocarbon should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("spiro saturated fused hydrocarbon should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    let saturated_aromatic_carbons = reparsed
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "C"
                && atom.aromatic
                && reparsed.graph().incident_bonds(*id).is_ok_and(|bonds| {
                    let bonds = bonds.collect::<Vec<_>>();
                    bonds
                        .iter()
                        .all(|(_, bond)| matches!(bond.order, BondOrder::Single))
                        && bonds.iter().any(|(_, bond)| !bond.aromatic)
                })
        })
        .count();
    assert_eq!(saturated_aromatic_carbons, 0, "{written}");
}

#[test]
fn fused_cyclic_imine_round_trip_keeps_imine_carbon_aliphatic() {
    let mut molecule =
        read_smiles("C1CN2CC3=CC=CC=C3N=C2[C@@H]1O").expect("fused cyclic imine should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused cyclic imine should sanitize");
    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused cyclic imine should canonicalize");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    let aromatic_imine_carbons = reparsed
        .graph()
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "C"
                && atom.aromatic
                && reparsed.graph().incident_bonds(*id).is_ok_and(|bonds| {
                    let bonds = bonds.collect::<Vec<_>>();
                    bonds.iter().any(|(_, bond)| {
                        matches!(bond.order, BondOrder::Double)
                            && reparsed
                                .graph()
                                .atom(bond.other_atom(*id))
                                .is_ok_and(|other| other.element.symbol() == "N")
                    }) && bonds.iter().any(|(_, bond)| {
                        matches!(bond.order, BondOrder::Single)
                            && reparsed
                                .graph()
                                .atom(bond.other_atom(*id))
                                .is_ok_and(|other| other.element.symbol() == "N")
                    })
                })
        })
        .count();
    assert_eq!(aromatic_imine_carbons, 0, "{written}");
}

#[test]
fn partially_saturated_carbonyl_fused_rings_stay_aliphatic() {
    let mut molecule =
        read_smiles("CC1(CC2=C(C(=O)C1)OC3=C(C2C4=CC=CC=C4[N+](=O)[O-])C(=O)CC(C3)(C)C)C")
            .expect("partially saturated fused carbonyl system should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("partially saturated fused carbonyl system should sanitize");

    for atom_id in [1, 2, 3, 4, 5, 7, 9, 10, 21, 23, 24, 25] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("partially saturated ring atom")
                .aromatic,
            "partially saturated ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [12, 13, 14, 15, 16, 17] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("phenyl atom")
                .aromatic,
            "phenyl atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_lactam_bridge_ring_stays_aliphatic() {
    let mut molecule = read_smiles("CCN1C2=C(C=C(C=C2OC3=C(C1=O)C=CC=N3)C)[N+](=O)[O-]")
        .expect("fused lactam bridge should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused lactam bridge should sanitize");

    for atom_id in [2, 9, 12] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("lactam bridge atom")
                .aromatic,
            "lactam bridge atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [3, 4, 5, 6, 7, 8, 10, 11, 14, 15, 16, 17] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("aromatic fused atom")
                .aromatic,
            "aromatic fused atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_pubchem_subset_aromaticity_remains_additive() {
    let mut molecule =
        read_smiles("CN1CCN(CC1)CCC2=CC3=C4N2C=C(C(=O)C4=CC(=C3)CN5CCOCC5)C(=O)NCC6=CC=C(C=C6)Cl")
            .expect("PubChem fused subset boundary should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("PubChem fused subset boundary should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    let aromatic_bonds = molecule
        .graph()
        .bonds()
        .filter(|(_, bond)| bond.aromatic)
        .count();
    assert_eq!(
        (aromatic_atoms, aromatic_bonds),
        (18, 20),
        "accepted fused subsystems should be marked additively"
    );
}

#[test]
fn fused_imine_sulfonamide_neighbor_ring_stays_aliphatic() {
    let mut molecule = read_smiles("O=C(NC1=CC2=NC=CN=C2C=C1)C1=CC=CN2CCS(=O)(=O)N=C12")
        .expect("fused imine sulfonamide record should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused imine sulfonamide record should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    let aromatic_bonds = molecule
        .graph()
        .bonds()
        .filter(|(_, bond)| bond.aromatic)
        .count();
    assert_eq!((aromatic_atoms, aromatic_bonds), (10, 11));

    for atom_id in [13, 14, 15, 16, 17, 24] {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("sulfonamide-adjacent ring atom")
                .aromatic,
            "sulfonamide-adjacent fused atom {atom_id} should stay aliphatic"
        );
    }
}

#[test]
fn fused_imide_heterocycle_keeps_only_phenyl_rings_aromatic() {
    let mut molecule =
        read_smiles("OOC1(CC2=CC=C(O)C=C2)N=C2C(CC3=CC=CC=C3)=NC(C3=CC=C(O)C=C3)=CN2C1=O")
            .expect("fused imide heterocycle should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused imide heterocycle should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    let aromatic_bonds = molecule
        .graph()
        .bonds()
        .filter(|(_, bond)| bond.aromatic)
        .count();
    let aromatic_nitrogens = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "N" && atom.aromatic)
        .count();
    assert_eq!((aromatic_atoms, aromatic_bonds), (18, 18));
    assert_eq!(aromatic_nitrogens, 0);
}

#[test]
fn fused_four_member_diketone_ring_can_be_aromatic() {
    let mut molecule = read_smiles("C1CSC2(C3=C(C=CC(=C3)Cl)OC4=C2C(=O)C4=O)SC1")
        .expect("fused four-member diketone should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused four-member diketone should sanitize");

    for atom_id in [12, 13, 14, 16] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("four-member diketone atom")
                .aromatic,
            "four-member diketone atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn large_conjugated_macrocycle_aromatic_core_is_not_size_skipped() {
    let mut molecule = read_smiles("CN(C)CCO.C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=N8)N=C2[N-]3.[Cu+2]")
        .expect("macrocycle salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("macrocycle salt should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert!(
        aromatic_atoms > 24,
        "large fused systems should not be constrained by a pre-RDKit-size cap"
    );
    assert_eq!(aromatic_atoms, 40);

    let copper = molecule.graph().atom(AtomId::new(46)).expect("copper atom");
    assert!(!copper.aromatic);
    assert_eq!(copper.formal_charge, 2);
}

#[test]
fn tetrahydroporphyrin_marks_each_conjugated_pyrrole_ring_aromatic() {
    let mut molecule = read_smiles(
        "CC=C1C(=C2C=C3C(=CC)C(=C(N3)C=C4C(=C(C(=CC5=C(C(=C(N5)C=C1N2)C)CCC(=O)O)N4)CCC(=O)O)C)C)C",
    )
    .expect("tetrahydroporphyrin should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("tetrahydroporphyrin should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .map(|(atom_id, atom)| (atom_id.index(), atom.element.symbol()))
        .collect::<Vec<_>>();
    let rings = molecule
        .graph()
        .ring_set()
        .expect("sanitization should retain the ring set")
        .rings()
        .iter()
        .map(|ring| {
            ring.atoms
                .iter()
                .map(|atom| {
                    let payload = molecule
                        .graph()
                        .atom(*atom)
                        .expect("ring atom should exist");
                    (
                        atom.index(),
                        payload.element.symbol(),
                        payload.explicit_hydrogens,
                        payload.implicit_hydrogens,
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(
        aromatic_atoms.len(),
        20,
        "atoms={aromatic_atoms:?} rings={rings:?}"
    );
    assert_eq!(
        aromatic_atoms
            .iter()
            .filter(|(_, symbol)| *symbol == "N")
            .count(),
        4,
        "{aromatic_atoms:?}"
    );
}

#[test]
fn fused_lone_pair_five_ring_with_macrocycle_pi_links_stays_aromatic() {
    let mut molecule = read_smiles(
        "CC1=C(C2=CC3=NC(=CC4=NC(=CC5=C(C(=C(N5)C=C1N2)C=C)C)C(=C4CCC(=O)O)C)C(=C3C)CCC(=O)O)C=C",
    )
    .expect("porphyrinoid macrocycle should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("porphyrinoid macrocycle should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    let aromatic_bonds = molecule
        .graph()
        .bonds()
        .filter(|(_, bond)| bond.aromatic)
        .count();
    assert_eq!((aromatic_atoms, aromatic_bonds), (20, 22));
}

#[test]
fn fused_five_electron_support_ring_keeps_outer_perimeter_aliphatic() {
    let mut molecule = read_smiles("CC1=C(C2=CC3=NC(=CC4=C(C(=C([N-]4)C=C5C(=C(C(=N5)C=C1[N-]2)C)C=C)C)C=C)C(=C3CCC(=O)O)C)CCC(=O)O.[Fe+2]")
    .expect("anionic macrocycle salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("anionic macrocycle salt should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    let aromatic_bonds = molecule
        .graph()
        .bonds()
        .filter(|(_, bond)| bond.aromatic)
        .count();
    assert_eq!((aromatic_atoms, aromatic_bonds), (20, 22));
}

#[test]
fn neutral_aza_macrocycle_core_stays_aliphatic() {
    let mut molecule = read_smiles(
        "C1=CC=C2C(=C1)C3=NC4=NC(=NC5=NC(=NC6=NC(=NC2=N3)C7=CC=CC=C76)C8=CC=CC=C85)C9=CC=CC=C94",
    )
    .expect("neutral aza macrocycle should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("neutral aza macrocycle should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 24);

    for atom_id in 6..=21 {
        assert!(
            !molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("neutral aza macrocycle atom")
                .aromatic,
            "neutral aza macrocycle atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [0, 1, 2, 3, 4, 5, 22, 23, 24, 25, 26, 27] {
        assert!(
            molecule
                .graph()
                .atom(AtomId::new(atom_id))
                .expect("benzene atom")
                .aromatic,
            "benzene atom {atom_id} should stay aromatic"
        );
    }
}

#[test]
fn fused_azo_indole_ring_keeps_explicit_hydrogen_nitrogen_aromatic() {
    let mut molecule = read_smiles("CN1C=NN(C)C1N=NC1=C(C2=CC=CC=C2)NC2=CC=CC=C12.[Cl-]")
        .expect("fused azo indole salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused azo indole salt should sanitize");

    let aromatic_atom_ids = molecule
        .graph()
        .atoms()
        .filter_map(|(atom_id, atom)| atom.aromatic.then_some(atom_id.index()))
        .collect::<Vec<_>>();
    assert_eq!(aromatic_atom_ids.len(), 15, "{aromatic_atom_ids:?}");
}

#[test]
fn fused_tertiary_amine_ring_does_not_extend_aromatic_core() {
    let mut molecule = read_smiles("CC(C)C[C@@H]1CN2CCC3=CC(=C(C=C3C2CC1=O)OC)O[11CH3]")
        .expect("fused tertiary amine record should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused tertiary amine record should sanitize");

    assert_eq!(
        molecule
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        6
    );
    assert!(molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "N")
        .all(|(_, atom)| !atom.aromatic));

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));
    assert_eq!(
        reparsed
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        6,
        "{written}"
    );
}

#[test]
fn fused_n_hydroxy_lactam_ring_stays_aromatic() {
    let mut molecule = read_smiles("CCCCCCCC1=CC2=C(C=C1)N(C=CC2=O)O")
        .expect("fused N-hydroxy lactam should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused N-hydroxy lactam should sanitize");

    let aromatic_atoms = molecule
        .graph()
        .atoms()
        .filter_map(|(atom_id, atom)| {
            atom.aromatic
                .then_some((atom_id.index(), atom.element.symbol()))
        })
        .collect::<Vec<_>>();
    assert_eq!(aromatic_atoms.len(), 10, "{aromatic_atoms:?}");
    assert!(molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "N")
        .all(|(_, atom)| atom.aromatic));

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));
    assert_eq!(
        reparsed
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        10,
        "{written}"
    );
}

#[test]
fn n_aryl_fused_pyrrole_ring_stays_aromatic() {
    let mut molecule =
        read_smiles("CCOC(=O)C1=C(N(C2=C1C=C(C=C2)OCC(C[NH2+]CC3=CC=CC=C3)O)C4=CC=CC=C4)C.[Cl-]")
            .expect("N-aryl fused pyrrole salt should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("N-aryl fused pyrrole salt should sanitize");

    let aromatic_neutral_nitrogens = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| {
            atom.element.symbol() == "N" && atom.formal_charge == 0 && atom.aromatic
        })
        .count();
    assert_eq!(aromatic_neutral_nitrogens, 1);

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));
    let reparsed_aromatic_neutral_nitrogens = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| {
            atom.element.symbol() == "N" && atom.formal_charge == 0 && atom.aromatic
        })
        .count();
    assert_eq!(reparsed_aromatic_neutral_nitrogens, 1, "{written}");
}

#[test]
fn fused_saturated_thioether_bridge_stays_aliphatic() {
    let mut molecule = read_smiles(
        "C1=CC=C(C=C1)C2=C(C(=O)C3=CC=CC=C3O2)OC(=O)C4=CC5=C(C=C4Cl)SC6=NC=CN6S5(=O)=O",
    )
    .expect("fused thioether bridge should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("fused thioether bridge should sanitize");

    let neutral_sulfur_aromatic_count = molecule
        .graph()
        .atoms()
        .filter(|(_, atom)| {
            atom.element.symbol() == "S" && atom.formal_charge == 0 && atom.aromatic
        })
        .count();
    assert_eq!(neutral_sulfur_aromatic_count, 0);

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));
    let reparsed_neutral_sulfur_aromatic_count = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| {
            atom.element.symbol() == "S" && atom.formal_charge == 0 && atom.aromatic
        })
        .count();
    assert_eq!(reparsed_neutral_sulfur_aromatic_count, 0, "{written}");
}

#[test]
fn canonical_smiles_prefers_sanitizable_lactone_candidate() {
    let mut molecule = read_smiles("CC[C@H]1[C@H](COC1=O)CC2=CN=CN2C.C=CC(=O)O")
        .expect("lactone imidazole mixture should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("lactone imidazole mixture should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    assert_eq!(reparsed.graph().atom_count(), molecule.graph().atom_count());
    assert_eq!(reparsed.graph().bond_count(), molecule.graph().bond_count());
}

#[test]
fn saturated_fused_benzodiazepinone_lactam_round_trip_stays_aliphatic() {
    let mut molecule = read_smiles("CN(C)CCN1C(NC(=O)C2=C1C=C(C=C2)Cl)C3=CC=C(C=C3)Cl.Cl")
        .expect("benzodiazepinone should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("benzodiazepinone should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(
        aromatic_atoms, 12,
        "canonical output should keep only the two benzene rings aromatic: {written}"
    );
    let aromatic_nitrogens = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "N" && atom.aromatic)
        .count();
    assert_eq!(aromatic_nitrogens, 0, "{written}");
}

#[test]
fn aromatic_pyridinium_smiles_sanitizes() {
    let mut molecule =
        read_smiles("CCCCCC(=O)C[n+]1ccccc1").expect("aromatic pyridinium should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("aromatic pyridinium should sanitize");

    let cationic_nitrogen = molecule
        .graph()
        .atoms()
        .find(|(_, atom)| atom.element.symbol() == "N" && atom.formal_charge > 0)
        .expect("pyridinium nitrogen should exist")
        .1;
    assert!(cationic_nitrogen.aromatic);

    let mut protonated = read_smiles("Nc1ccc[nH+]c1").expect("protonated pyridinium should parse");
    perception_api::sanitize_with_options(&mut protonated, SanitizeOptions::default())
        .expect("protonated pyridinium should sanitize");
    assert!(protonated
        .graph()
        .atoms()
        .any(|(_, atom)| atom.element.symbol() == "N" && atom.formal_charge > 0 && atom.aromatic));

    let mut anionic = read_smiles("c1[n-]cnn1").expect("anionic aromatic nitrogen should parse");
    perception_api::sanitize_with_options(&mut anionic, SanitizeOptions::default())
        .expect("anionic aromatic nitrogen should sanitize");
    assert!(anionic
        .graph()
        .atoms()
        .any(|(_, atom)| atom.element.symbol() == "N" && atom.formal_charge < 0 && atom.aromatic));
}

#[test]
fn aromatic_pyrone_canonical_smiles_sanitizes() {
    let mut molecule = read_smiles("CC#CC#Cc1cccc(=O)o1").expect("aromatic pyrone should parse");

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("aromatic pyrone should sanitize");
}

#[test]
fn canonical_smiles_preserves_metal_bound_bracket_hydrogens() {
    let mut molecule = read_smiles("CC[Hg+]").expect("organomercury SMILES parses");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("organomercury SMILES sanitizes");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");

    assert!(written.contains("[CH2][Hg+]"), "{written}");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .expect("canonical output should sanitize");
    let metal_bound_carbon = reparsed
        .graph()
        .atoms()
        .find(|(atom_id, atom)| {
            atom.element.symbol() == "C"
                && reparsed
                    .graph()
                    .incident_bonds(*atom_id)
                    .expect("atom should be live")
                    .any(|(_, bond)| {
                        let neighbor_id = bond.other_atom(*atom_id);
                        reparsed
                            .graph()
                            .atom(neighbor_id)
                            .is_ok_and(|neighbor| neighbor.element.symbol() == "Hg")
                    })
        })
        .expect("canonical output should retain a carbon-mercury bond")
        .1;
    assert_eq!(metal_bound_carbon.explicit_hydrogens, 2);
    assert!(metal_bound_carbon.no_implicit_hydrogens);
    assert_eq!(metal_bound_carbon.implicit_hydrogens, Some(0));

    let mut thallium = read_smiles("C[Tl](C)C").expect("organothallium SMILES parses");
    perception_api::sanitize_with_options(&mut thallium, SanitizeOptions::default())
        .expect("organothallium SMILES sanitizes");
    let thallium_written =
        smiles_api::write_canonical_with_options(&thallium, CanonicalSmilesWriteOptions)
            .expect("organothallium canonical SMILES should write");
    assert_eq!(
        thallium_written.matches("[CH3]").count(),
        3,
        "{thallium_written}"
    );

    let mut antimony = read_smiles("C[Sb](C)C").expect("organoantimony SMILES parses");
    perception_api::sanitize_with_options(&mut antimony, SanitizeOptions::default())
        .expect("organoantimony SMILES sanitizes");
    let antimony_written = smiles_api::write_canonical(&antimony)
        .expect("organoantimony canonical SMILES should write");
    assert_eq!(
        antimony_written.matches("[CH3]").count(),
        3,
        "{antimony_written}"
    );
}

#[test]
fn canonical_smiles_materializes_hydrogen_on_bracketed_hypervalent_phosphorus() {
    let mut molecule = read_smiles("OP(=O)O").expect("phosphorous acid should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("phosphorous acid should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical phosphorous acid should write");
    assert!(written.contains("[PH]"), "{written}");

    let mut reparsed = read_smiles(&written).expect("canonical phosphorous acid should reparse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .expect("canonical phosphorous acid should resanitize");
    let phosphorus = reparsed
        .graph()
        .atoms()
        .find(|(_, atom)| atom.element.symbol() == "P")
        .expect("phosphorus should remain")
        .1;
    assert_eq!(phosphorus.explicit_hydrogens, 1);
    assert_eq!(phosphorus.implicit_hydrogens, Some(0));
    assert!(phosphorus.no_implicit_hydrogens);
}

#[test]
fn canonical_substituted_pyridinium_round_trip_sanitizes() {
    let input = "CCCCCC(=O)C[N+]1=CC=CC=C1.C1(C(=O)NC(=O)NC1=O)[N+](=O)[O-]";
    let mut molecule = read_smiles(input).expect("pyridinium regression parses");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("pyridinium regression sanitizes");
    let written = smiles_api::write_canonical(&molecule).expect("canonical SMILES writes");
    let mut reparsed = read_smiles(&written).expect("canonical SMILES reparses");
    let result = perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default());
    assert!(result.is_ok(), "{written}: {result:#?}");
}

#[test]
fn canonical_pubchem_100k_main_group_regressions_resanitize() {
    for input in [
        "C1=CC=C2C(=C1)NC3=CC=CC=C3[As]2O[As]4C5=CC=CC=C5NC6=CC=CC=C64",
        "C1=CC2=C(C=C1Cl)[I+]C3=C(O2)C=CC(=C3)Cl.OS(=O)(=O)[O-]",
        "CC(C)CC(=O)OCC1=COC=C2C1=CC=C2C=O",
        "C1=CC=C(C=C1)C2=CC(=[O+]C(=C2)C3=CC=CC=C3)C4=CC=CC=C4.[O-]Cl(=O)(=O)=O",
        "C1=CC=PC=C1",
        "C1=C(C2=COC=C(C2=C1)CO)C=O",
        "CC(=O)OCC1=COC=C2C1=CC=C2C=O",
        "CC1=CC(=CC(=[O+]1)C)SCC=C",
        "C1=CC=C2C(=C1)C3=CC=CC=C3[Si]2(C4=CC(=CC=C4)Br)F",
        "CSC1=C(C(=[S+]S1)SC)O",
        "C[Si](C)(C)C1=CC=C(S1)C2=C3C[Si]4(CC3=C(S2)C5=CC=C(S5)[Si](C)(C)C)CC6=C(SC(=C6C4)C7=CC=C(S7)[Si](C)(C)C)C8=CC=C(S8)[Si](C)(C)C",
        "C1C2=CC=CC=C2[Sn]3(C4=CC=CC=C4CN1CC5=CC=CC=C53)Br",
        "C1=CC=C2C(=C1)NC3=CC=CC=C3[AsH]2(N)Cl",
        "C[N+](C)(C)C.C[Si-]12(C3=CC=CC=C3C(O1)(C(F)(F)F)C(F)(F)F)C4=CC=CC=C4C(O2)(C(F)(F)F)C(F)(F)F",
        "C1=CC2=C3C(=C1)[I+]C4=CC=CC(=C43)[I+]2",
        "C[Si]1(CC2=CC=CC=C2C1)[Si](C)(C)C",
    ] {
        let mut molecule = read_smiles(input)
            .unwrap_or_else(|error| panic!("input should parse: {input}: {error}"));
        perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
            .unwrap_or_else(|error| panic!("input should sanitize: {input}: {error:#?}"));
        let written = smiles_api::write_canonical(&molecule)
            .unwrap_or_else(|error| panic!("canonical output should write: {input}: {error}"));
        let mut reparsed = read_smiles(&written)
            .unwrap_or_else(|error| panic!("canonical output should parse: {written}: {error}"));
        perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
            .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error:#?}"));
    }
}

#[test]
fn canonical_aryl_germanium_round_trip_preserves_no_implicit_aromatic_carbon() {
    let mut molecule =
        read_smiles("C1=CC=C(C=C1)[Ge](Cl)(Cl)Cl").expect("aryl germanium SMILES parses");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("aryl germanium SMILES sanitizes");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("aryl germanium canonical SMILES should write");
    assert!(written.contains("[c]"), "{written}");

    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .expect("canonical output should sanitize");

    let germanium_bound_carbon = reparsed
        .graph()
        .atoms()
        .find(|(atom_id, atom)| {
            atom.element.symbol() == "C"
                && atom.aromatic
                && reparsed
                    .graph()
                    .incident_bonds(*atom_id)
                    .expect("atom should be live")
                    .any(|(_, bond)| {
                        let neighbor_id = bond.other_atom(*atom_id);
                        reparsed
                            .graph()
                            .atom(neighbor_id)
                            .is_ok_and(|neighbor| neighbor.element.symbol() == "Ge")
                    })
        })
        .expect("canonical output should retain an aryl germanium bond")
        .1;
    assert!(germanium_bound_carbon.no_implicit_hydrogens, "{written}");
}

#[test]
fn canonical_aryl_tin_round_trip_preserves_no_implicit_aromatic_carbons() {
    let mut molecule =
        read_smiles("C1=CC=C(C=C1)[SnH](C2=CC=CC=C2)Cl").expect("aryl tin SMILES parses");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("aryl tin SMILES sanitizes");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("aryl tin canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .expect("canonical output should sanitize");

    let tin_bound_aromatic_carbons = reparsed
        .graph()
        .atoms()
        .filter(|(atom_id, atom)| {
            atom.element.symbol() == "C"
                && atom.aromatic
                && reparsed
                    .graph()
                    .incident_bonds(*atom_id)
                    .expect("atom should be live")
                    .any(|(_, bond)| {
                        let neighbor_id = bond.other_atom(*atom_id);
                        reparsed
                            .graph()
                            .atom(neighbor_id)
                            .is_ok_and(|neighbor| neighbor.element.symbol() == "Sn")
                    })
        })
        .map(|(_, atom)| atom)
        .collect::<Vec<_>>();
    assert_eq!(tin_bound_aromatic_carbons.len(), 2, "{written}");
    assert!(
        tin_bound_aromatic_carbons
            .iter()
            .all(|atom| atom.no_implicit_hydrogens && atom.implicit_hydrogens == Some(0)),
        "{written}"
    );
}

#[test]
fn cationic_thiadiazolium_imine_canonical_round_trip_sanitizes() {
    let mut molecule = read_smiles("CN(C1=NC(=[N+](C)C)SS1)C(=S)SC")
        .expect("cationic thiadiazolium imine should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("cationic thiadiazolium imine should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    let aromatic_ring_atoms = reparsed
        .graph()
        .atoms()
        .filter(|(_, atom)| matches!(atom.element.symbol(), "C" | "N" | "S") && atom.aromatic)
        .count();
    assert_eq!(aromatic_ring_atoms, 5, "{written}");
}

#[test]
fn canonical_multicomponent_oxygen_neighbors_match_after_round_trip() {
    let mut molecule = read_smiles("CC(CO)O.CC(C)(C)CCCCC(CC1CO1)C(=O)O.C1=CC=C2C(=C1)C(=O)OC2=O.C1=CC2=C(C=C1C(=O)O)C(=O)OC2=O.C(CCC(=O)O)CC(=O)O")
    .expect("oxygen-rich mixture should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("oxygen-rich mixture should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    assert_eq!(
        local_atom_neighbor_signatures(molecule.graph()),
        local_atom_neighbor_signatures(reparsed.graph()),
        "{written}"
    );
}

#[test]
fn canonical_pubchem_macrocycle_anionic_nitrogen_round_trip_matches_neighbors() {
    let mut molecule = read_smiles("CN(C)CCO.C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=N8)N=C2[N-]3.[Cu+2]")
    .expect("PubChem macrocycle mixture should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("PubChem macrocycle mixture should sanitize");

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    assert_eq!(
        local_atom_neighbor_signatures(molecule.graph()),
        local_atom_neighbor_signatures(reparsed.graph()),
        "{written}"
    );
}

#[test]
fn canonical_substituted_pyrrole_uses_aromatic_nitrogen_form() {
    let mut molecule = read_smiles("CCOC(=O)C1=C(C(=C(N1)C)C(=O)OC(C)(C)C)C")
        .expect("substituted pyrrole should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("substituted pyrrole should sanitize");
    let nitrogen = molecule
        .graph()
        .atoms()
        .find_map(|(_, atom)| (atom.element.symbol() == "N").then_some(atom))
        .expect("substituted pyrrole nitrogen");
    assert!(nitrogen.aromatic);
    assert_eq!(nitrogen.explicit_hydrogens, 1);
    assert_eq!(nitrogen.implicit_hydrogens, Some(0));
    assert!(!nitrogen.no_implicit_hydrogens);

    let written = smiles_api::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed = read_smiles(&written).expect("canonical output should parse");
    perception_api::sanitize_with_options(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|error| panic!("canonical output should sanitize: {written}: {error}"));

    assert!(written.contains("[nH]"), "{written}");
    assert!(!written.contains("-N-"), "{written}");
}

type TestAtomStateSignature = (u8, i8, u16, u8, u8, bool, bool);
type TestAtomNeighborSignature = (
    TestAtomStateSignature,
    Vec<(TestAtomStateSignature, u8, bool)>,
);

fn local_atom_neighbor_signatures(mol: &Molecule) -> Vec<TestAtomNeighborSignature> {
    local_atom_neighbor_signatures_with(mol, test_atom_state_signature)
}

fn local_atom_neighbor_signatures_ignoring_halogen_no_implicit(
    mol: &Molecule,
) -> Vec<TestAtomNeighborSignature> {
    local_atom_neighbor_signatures_with(mol, test_atom_state_signature_ignoring_halogen_no_implicit)
}

fn local_atom_neighbor_signatures_with(
    mol: &Molecule,
    atom_signature: fn(&Atom) -> TestAtomStateSignature,
) -> Vec<TestAtomNeighborSignature> {
    let mut atoms = mol
        .atoms()
        .map(|(id, atom)| {
            let mut neighbors = mol
                .incident_bonds(id)
                .expect("atom should be live")
                .map(|(_, bond)| {
                    let neighbor = mol
                        .atom(bond.other_atom(id))
                        .expect("neighbor should be live");
                    (
                        atom_signature(neighbor),
                        test_semantic_bond_order_code(bond),
                        bond.aromatic,
                    )
                })
                .collect::<Vec<_>>();
            neighbors.sort_unstable();
            (atom_signature(atom), neighbors)
        })
        .collect::<Vec<_>>();
    atoms.sort_unstable();
    atoms
}

fn test_atom_state_signature(atom: &Atom) -> TestAtomStateSignature {
    (
        atom.element.atomic_number(),
        atom.formal_charge,
        atom.isotope.unwrap_or_default(),
        atom.explicit_hydrogens,
        atom.implicit_hydrogens.unwrap_or_default(),
        atom.no_implicit_hydrogens,
        atom.aromatic,
    )
}

fn test_atom_state_signature_ignoring_halogen_no_implicit(atom: &Atom) -> TestAtomStateSignature {
    let mut signature = test_atom_state_signature(atom);
    if matches!(atom.element.symbol(), "F" | "Cl" | "Br" | "I") {
        signature.5 = false;
    }
    signature
}

fn test_semantic_bond_order_code(bond: &Bond) -> u8 {
    if bond.aromatic {
        test_bond_order_code(BondOrder::Aromatic)
    } else {
        test_bond_order_code(bond.order)
    }
}

fn test_bond_order_code(order: BondOrder) -> u8 {
    match order {
        BondOrder::Zero => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 5,
        BondOrder::Dative => 6,
    }
}

#[test]
fn smiles_writer_rejects_lossy_bonds_and_stereo() {
    let mut molecule = SmallMolecule::default();
    let a = molecule.graph_mut().add_atom(carbon());
    let b = molecule.graph_mut().add_atom(carbon());
    let bond = molecule
        .graph_mut()
        .add_bond(a, b, BondOrder::Dative)
        .expect("bond");
    assert!(
        smiles_api::write_with_options(&molecule, SmilesWriteOptions)
            .expect_err("dative bond should be rejected")
            .message
            .contains("cannot encode")
    );

    molecule.graph_mut().bond_mut(bond).expect("bond").order = BondOrder::Single;
    molecule
        .graph_mut()
        .set_stereo_bond_mark(StereoBondMark {
            bond,
            kind: StereoBondMarkKind::DirectionalUp,
            source: StereoSource::Smiles,
        })
        .expect("stereo mark");
    assert!(
        smiles_api::write_with_options(&molecule, SmilesWriteOptions)
            .expect_err("stereo should be rejected")
            .message
            .contains("stereochemistry")
    );

    molecule
        .graph_mut()
        .clear_stereo_bond_mark(bond)
        .expect("clear mark");
    molecule
        .graph_mut()
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center: a,
                carriers: vec![StereoCarrier::Atom(b), StereoCarrier::ImplicitHydrogen],
                orientation: TetrahedralOrientation::Clockwise,
            }),
            StereoSource::User,
        ))
        .expect("atom stereo");
    assert!(
        smiles_api::write_with_options(&molecule, SmilesWriteOptions)
            .expect_err("atom chirality should be rejected")
            .message
            .contains("atom stereochemistry")
    );

    let element = molecule
        .graph()
        .stereo_element_ids()
        .next()
        .expect("stereo element");
    molecule
        .graph_mut()
        .remove_stereo_element(element)
        .expect("remove atom stereo");
    {
        let mut atom = molecule.graph_mut().atom_mut(a).expect("atom");
        atom.radical = Some(AtomRadical::Doublet);
        atom.explicit_hydrogens = 2;
        atom.no_implicit_hydrogens = true;
    }
    let radical_smiles = smiles_api::write_with_options(&molecule, SmilesWriteOptions)
        .expect("valence-consistent radical should write");
    assert!(radical_smiles.contains("[CH2]"));
    let radical_reparsed =
        read_smiles(&radical_smiles).expect("radical writer output should parse");
    assert!(radical_reparsed
        .graph()
        .atoms()
        .any(|(_, atom)| atom.radical == Some(AtomRadical::Doublet)));

    {
        let mut atom = molecule.graph_mut().atom_mut(a).expect("atom");
        atom.radical = None;
        atom.explicit_hydrogens = 0;
        atom.no_implicit_hydrogens = true;
    }
    let written = smiles_api::write_with_options(&molecule, SmilesWriteOptions)
        .expect("no-implicit-hydrogen atom should write");
    assert!(written.contains("[C]"));
    let reparsed = read_smiles(&written).expect("writer output should parse");
    assert!(reparsed
        .graph()
        .atoms()
        .any(|(_, atom)| atom.no_implicit_hydrogens));
}

#[test]
fn bracket_atoms_infer_rdkit_radical_multiplicity_from_valence_deficit() {
    for (smiles, atom_index, expected) in [
        ("[C]", 0, AtomRadical::Quintet),
        ("[C]C", 0, AtomRadical::Quartet),
        ("C=[C]", 1, AtomRadical::Triplet),
        ("C#[C]", 1, AtomRadical::Doublet),
        ("[N]", 0, AtomRadical::Quartet),
        ("[O]", 0, AtomRadical::Triplet),
    ] {
        let molecule = read_smiles(smiles).expect("radical SMILES should parse");
        assert_eq!(
            molecule
                .graph()
                .atom(AtomId::new(atom_index))
                .expect("bracket atom")
                .radical,
            Some(expected),
            "{smiles}"
        );
    }

    let aromatic_radical = read_smiles("[c]1ccccc1").expect("aromatic carbon radical parses");
    assert_eq!(
        aromatic_radical
            .graph()
            .atom(AtomId::new(0))
            .expect("aromatic radical")
            .radical,
        Some(AtomRadical::Doublet)
    );
    let substituted_pyridinium =
        read_smiles("C[n+]1ccccc1").expect("substituted pyridinium parses");
    assert_eq!(
        substituted_pyridinium
            .graph()
            .atom(AtomId::new(1))
            .expect("pyridinium nitrogen")
            .radical,
        None
    );
}

#[test]
fn isomeric_smiles_writes_tetrahedral_elements_from_stereo_model() {
    let molecule = read_smiles("F[C@H](Cl)Br").expect("tetrahedral SMILES should parse");

    let written = smiles_api::write_isomeric(&molecule).expect("tetrahedral stereo should write");

    assert_eq!(written, "F[C@H](Cl)Br");
    let reparsed = read_smiles(&written).expect("isomeric output should parse");
    let stereo = reparsed
        .graph()
        .stereo_elements()
        .map(|(_, element)| element)
        .collect::<Vec<_>>();
    assert_eq!(stereo.len(), 1);
    match &stereo[0].kind {
        StereoElementKind::Tetrahedral(tetrahedral) => {
            assert_eq!(tetrahedral.center, AtomId::new(1));
            assert_eq!(tetrahedral.orientation, TetrahedralOrientation::Clockwise);
            assert_eq!(
                tetrahedral.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(0)),
                    StereoCarrier::ImplicitHydrogen,
                    StereoCarrier::Atom(AtomId::new(2)),
                    StereoCarrier::Atom(AtomId::new(3)),
                ]
            );
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn isomeric_smiles_flips_tetrahedral_marker_for_odd_writer_carrier_order() {
    let mut molecule = read_smiles("F[C@H](Cl)Br").expect("tetrahedral SMILES should parse");
    let element = molecule
        .graph()
        .stereo_element_ids()
        .next()
        .expect("stereo element");
    molecule
        .graph_mut()
        .remove_stereo_element(element)
        .expect("remove parsed stereo");
    molecule
        .graph_mut()
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center: AtomId::new(1),
                carriers: vec![
                    StereoCarrier::Atom(AtomId::new(0)),
                    StereoCarrier::ImplicitHydrogen,
                    StereoCarrier::Atom(AtomId::new(3)),
                    StereoCarrier::Atom(AtomId::new(2)),
                ],
                orientation: TetrahedralOrientation::Clockwise,
            }),
            StereoSource::User,
        ))
        .expect("replacement stereo element");

    let written = smiles_api::write_isomeric(&molecule).expect("tetrahedral stereo should write");

    assert_eq!(written, "F[C@@H](Cl)Br");
}

#[test]
fn isomeric_smiles_rejects_unencoded_stereo_layers() {
    let directional = read_smiles("C/C=C\\C").expect("directional bond markers should parse");
    assert!(smiles_api::write_isomeric(&directional)
        .expect_err("unperceived source bond marks should be rejected")
        .message
        .contains("source bond marks"));

    let mut unknown = read_smiles("F[C@H](Cl)Br").expect("tetrahedral SMILES should parse");
    let element = unknown
        .graph()
        .stereo_element_ids()
        .next()
        .expect("stereo element");
    unknown
        .graph_mut()
        .stereo_element_mut(element)
        .expect("stereo element")
        .specifiedness = StereoSpecifiedness::Unknown;
    assert!(smiles_api::write_isomeric(&unknown)
        .expect_err("unknown stereo should be rejected")
        .message
        .contains("unknown stereo"));
}

#[test]
fn isomeric_smiles_writes_directional_double_bond_elements() {
    for (input, expected_output, expected_orientation) in [
        ("C/C=C\\C", "C\\C=C/C", DoubleBondOrientation::Together),
        ("C/C=C/C", "C\\C=C\\C", DoubleBondOrientation::Opposite),
    ] {
        let mut molecule = read_smiles(input).expect("directional alkene should parse");
        perception_api::sanitize(&mut molecule).expect("directional alkene should sanitize");

        let written =
            smiles_api::write_isomeric(&molecule).expect("double-bond stereo should write");

        assert_eq!(written, expected_output);
        let mut reparsed = read_smiles(&written).expect("isomeric alkene output should parse");
        perception_api::sanitize(&mut reparsed).expect("isomeric alkene output should sanitize");
        let stereo = reparsed
            .graph()
            .stereo_elements()
            .filter_map(|(_, element)| match &element.kind {
                StereoElementKind::DoubleBond(stereo) => Some(stereo),
                StereoElementKind::Tetrahedral(_) | StereoElementKind::Axis(_) => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(stereo.len(), 1);
        assert_eq!(stereo[0].orientation, expected_orientation);
    }
}

#[test]
fn isomeric_smiles_writes_pubchem_conjugated_directional_polyene() {
    let mut molecule = read_smiles("CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(C)C=C)/C")
        .expect("directional polyene should parse");
    perception_api::sanitize(&mut molecule).expect("directional polyene should sanitize");

    let written = smiles_api::write_isomeric(&molecule).expect("directional polyene should write");

    let mut reparsed = read_smiles(&written).expect("isomeric polyene output should parse");
    perception_api::sanitize(&mut reparsed).expect("isomeric polyene output should sanitize");
    assert!(reparsed.graph().stereo_elements().next().is_some());
}

#[test]
fn isomeric_smiles_preserves_pubchem_fused_quaternary_center() {
    let mut molecule = read_smiles("C[C@]12CCCC(C1CCC3=CC(=C(C=C23)C(=O)OC)C(=O)OC)(C)C")
        .expect("fused quaternary center should parse");
    perception_api::sanitize(&mut molecule).expect("fused quaternary center should sanitize");
    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut());
    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.assigned[0].descriptor, StereoDescriptor::S);

    let written =
        smiles_api::write_isomeric(&molecule).expect("fused quaternary center should write");
    let mut reparsed = read_smiles(&written).expect("isomeric fused center output should parse");
    perception_api::sanitize(&mut reparsed).expect("isomeric fused center output should sanitize");
    let report = stereo_api::assign_cip_descriptors(reparsed.graph_mut());
    assert!(report.is_ok(), "{:?}", report.issues);

    assert_eq!(report.assigned[0].descriptor, StereoDescriptor::S);
}

#[test]
fn isomeric_smiles_round_trips_pubchem_anthraquinone_aromatic_shape() {
    let mut molecule = read_smiles("CC1C(C(CC(O1)O[C@H]2C[C@@](CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)(C(=O)C)O)N=C(CCSSCCC(=NC6CC(OC(C6O)C)O[C@H]7C[C@@](CC8=C7C(=C9C(=C8O)C(=O)C1=C(C9=O)C(=CC=C1)OC)O)(C(=O)C)O)N)N)O")
    .expect("anthraquinone source should parse");
    perception_api::sanitize(&mut molecule).expect("anthraquinone source should sanitize");
    let written = smiles_api::write_isomeric(&molecule).expect("anthraquinone should write");
    let mut reparsed = read_smiles(&written).expect("anthraquinone isomeric output should parse");
    perception_api::sanitize(&mut reparsed).expect("anthraquinone isomeric output should sanitize");

    assert_eq!(
        reparsed
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        molecule
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        "{written}"
    );
    assert_eq!(
        reparsed
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.no_implicit_hydrogens)
            .count(),
        molecule
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.no_implicit_hydrogens)
            .count(),
        "{written}"
    );
}

#[test]
fn isomeric_smiles_writes_implicit_carrier_double_bond_elements() {
    for (left_carrier, right_carrier, stored_orientation, expected_orientation) in [
        (
            StereoCarrier::ImplicitHydrogen,
            StereoCarrier::ImplicitHydrogen,
            DoubleBondOrientation::Together,
            DoubleBondOrientation::Together,
        ),
        (
            StereoCarrier::ImplicitHydrogen,
            StereoCarrier::Atom(AtomId::new(3)),
            DoubleBondOrientation::Together,
            DoubleBondOrientation::Opposite,
        ),
    ] {
        let mut molecule = SmallMolecule::default();
        let mut left_atom = carbon();
        left_atom.implicit_hydrogens = Some(1);
        let left = molecule.graph_mut().add_atom(left_atom);
        let mut right_atom = carbon();
        right_atom.implicit_hydrogens = Some(1);
        let right = molecule.graph_mut().add_atom(right_atom);
        let fluorine = molecule.graph_mut().add_atom(element_atom("F"));
        let chlorine = molecule.graph_mut().add_atom(element_atom("Cl"));
        molecule
            .graph_mut()
            .add_bond(left, fluorine, BondOrder::Single)
            .expect("left carrier bond");
        let double_bond = molecule
            .graph_mut()
            .add_bond(left, right, BondOrder::Double)
            .expect("double bond");
        molecule
            .graph_mut()
            .add_bond(right, chlorine, BondOrder::Single)
            .expect("right carrier bond");
        molecule
            .graph_mut()
            .add_stereo_element(StereoElement::specified(
                StereoElementKind::DoubleBond(DoubleBondStereo {
                    bond: double_bond,
                    left,
                    right,
                    left_carrier,
                    right_carrier,
                    orientation: stored_orientation,
                }),
                StereoSource::User,
            ))
            .expect("double-bond stereo");

        let written =
            smiles_api::write_isomeric(&molecule).expect("implicit-carrier stereo should write");
        assert!(
            written.contains('/') || written.contains('\\'),
            "isomeric output should contain directional marks: {written}"
        );

        let mut reparsed = read_smiles(&written).expect("isomeric alkene output should parse");
        perception_api::sanitize(&mut reparsed).expect("isomeric alkene output should sanitize");
        let stereo = reparsed
            .graph()
            .stereo_elements()
            .filter_map(|(_, element)| match &element.kind {
                StereoElementKind::DoubleBond(stereo) => Some(stereo),
                StereoElementKind::Tetrahedral(_) | StereoElementKind::Axis(_) => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(stereo.len(), 1);
        assert_eq!(stereo[0].orientation, expected_orientation);
    }
}

#[test]
fn smiles_writer_rejects_more_ring_labels_than_parser_supports() {
    let mut molecule = SmallMolecule::default();
    let atoms = (0..16)
        .map(|_| molecule.graph_mut().add_atom(carbon()))
        .collect::<Vec<_>>();
    for left in 0..atoms.len() {
        for right in (left + 1)..atoms.len() {
            molecule
                .graph_mut()
                .add_bond(atoms[left], atoms[right], BondOrder::Single)
                .expect("complete graph bond should be valid");
        }
    }

    assert!(
        smiles_api::write_with_options(&molecule, SmilesWriteOptions)
            .expect_err("more than 99 ring closures should be rejected")
            .message
            .contains("at most 99")
    );
}
