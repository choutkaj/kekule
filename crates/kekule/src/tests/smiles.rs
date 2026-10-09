use super::*;

fn assigned_descriptors(molecule: &mut Molecule) -> Vec<StereoDescriptor> {
    stereo_api::assign_cip_descriptors(molecule)
        .unwrap()
        .assigned
        .into_iter()
        .map(|assignment| assignment.descriptor)
        .collect()
}

#[test]
fn smiles_accepts_zero_counts_and_explicit_tetrahedral_classes() {
    for (source, equivalent) in [
        ("[0CH0+0]", "[C]"),
        ("[C@TH1H](F)(Cl)Br", "[C@H](F)(Cl)Br"),
        ("[C@TH2H](F)(Cl)Br", "[C@@H](F)(Cl)Br"),
    ] {
        let left = read_smiles(source).unwrap();
        let right = read_smiles(equivalent).unwrap();
        assert_eq!(
            left.atoms().map(|(_, atom)| atom).collect::<Vec<_>>(),
            right.atoms().map(|(_, atom)| atom).collect::<Vec<_>>()
        );
        assert_eq!(
            left.stereo_elements()
                .map(|(_, value)| value)
                .collect::<Vec<_>>(),
            right
                .stereo_elements()
                .map(|(_, value)| value)
                .collect::<Vec<_>>()
        );
    }
    let mapped = read_smiles("[CH4:0]").unwrap();
    assert_eq!(mapped.atoms().next().unwrap().1.atom_map, Some(0));
    assert_eq!(
        smiles_api::write(&mapped, smiles_api::SmilesWriteOptions::isomeric()).unwrap(),
        "[CH4:0]"
    );
}

#[test]
fn smiles_components_follow_connectivity_across_fragments_and_branches() {
    for (source, counts) in [
        ("C1.C1", vec![2]),
        ("C(C.C)", vec![2, 1]),
        ("C(.O)N", vec![2, 1]),
        ("C1.C2.C12", vec![3]),
    ] {
        let components = read_smiles_components(source).unwrap();
        assert_eq!(
            components
                .iter()
                .map(|molecule| molecule.atom_count())
                .collect::<Vec<_>>(),
            counts,
            "{source}"
        );
        assert_eq!(
            components
                .iter()
                .map(|molecule| molecule.bond_count())
                .sum::<usize>(),
            counts.iter().sum::<usize>() - counts.len()
        );
    }
}

#[test]
fn smiles_ring_direction_is_relative_to_each_written_endpoint() {
    let source = r"C/1=C/CCCCCC\1";
    let mut mol = read_smiles(source).unwrap();
    perceive(&mut mol).unwrap();
    let descriptors = assigned_descriptors(&mut mol);
    assert_eq!(descriptors.len(), 1);
    for output in [
        smiles_api::write(&mol, smiles_api::SmilesWriteOptions::isomeric()).unwrap(),
        smiles_api::write(&mol, smiles_api::SmilesWriteOptions::canonical()).unwrap(),
    ] {
        let mut restored = read_smiles(&output).unwrap();
        perceive(&mut restored).unwrap();
        assert_eq!(assigned_descriptors(&mut restored), descriptors, "{output}");
    }
    assert!(smiles_api::parse_str("C/1=C/CCCCCC/1").is_err());
}

#[test]
fn smiles_rejects_directional_projection_that_invents_an_alkene_assertion() {
    let mut mol = read_smiles("F/C=C/C=C/C=C/F").unwrap();
    let middle = mol.stereo_elements().find_map(|(id, element)| matches!(&element.kind, StereoElementKind::DoubleBond(value) if value.bond == BondId::new(3)).then_some(id)).unwrap();
    mol.remove_stereo_element(middle).unwrap();
    perceive(&mut mol).unwrap();
    assert_eq!(mol.stereo_elements().count(), 2);
    assert!(
        smiles_api::write(&mol, smiles_api::SmilesWriteOptions::isomeric())
            .unwrap_err()
            .message()
            .contains("stereo assertion")
    );
    assert!(
        smiles_api::write(&mol, smiles_api::SmilesWriteOptions::canonical())
            .unwrap_err()
            .message()
            .contains("stereo assertion")
    );
}

#[test]
fn smiles_quadruple_bonds_round_trip_without_loss() {
    let molecule = read_smiles("[Mo]$[Mo]").unwrap();
    for written in [
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::default()).unwrap(),
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()).unwrap(),
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()).unwrap(),
    ] {
        assert_eq!(
            read_smiles(&written)
                .unwrap()
                .bonds()
                .next()
                .unwrap()
                .1
                .order,
            BondOrder::Quadruple
        );
    }
}

#[test]
fn phosphine_smiles_preserve_bracket_hydrogen_and_lone_pair() {
    for (source, expected) in [
        ("C[P@H]C1CCCCC1", StereoDescriptor::R),
        ("C[P@@H]C1CCCCC1", StereoDescriptor::S),
        ("[P@H](C)C1CCCCC1", StereoDescriptor::R),
        ("[P@@H](C)C1CCCCC1", StereoDescriptor::S),
    ] {
        let mut molecule = read_smiles(source).expect("bracket-H phosphine interprets");
        let element = molecule.stereo_elements().next().unwrap().1;
        let StereoElementKind::Tetrahedral(stereo) = &element.kind else {
            panic!("expected tetrahedral phosphine")
        };
        assert!(stereo.carriers.contains(&StereoCarrier::ImplicitHydrogen));
        assert!(stereo.carriers.contains(&StereoCarrier::ImplicitLonePair));
        perceive(&mut molecule).unwrap();
        let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
        assert_eq!(assigned.assigned[0].descriptor, expected, "{source}");

        let written = smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric())
            .expect("phosphine writes");
        let mut reparsed = read_smiles(&written).expect("phosphine output interprets");
        perceive(&mut reparsed).unwrap();
        let assigned = stereo_api::assign_cip_descriptors(&mut reparsed).unwrap();
        assert_eq!(assigned.assigned[0].descriptor, expected, "{written}");
    }
}

fn aromatic_atom(molecule: &Molecule, atom: AtomId) -> bool {
    molecule.atom_is_aromatic(atom).expect("atom exists") == Some(true)
}

fn aromatic_bond(molecule: &Molecule, bond: BondId) -> bool {
    molecule.bond_is_aromatic(bond).expect("bond exists") == Some(true)
}

fn inferred_hydrogens(molecule: &Molecule, atom: AtomId) -> Option<u8> {
    molecule.inferred_hydrogens(atom).expect("atom exists")
}

fn aromatic_atom_count(molecule: &Molecule) -> usize {
    molecule
        .atom_ids()
        .filter(|atom| aromatic_atom(molecule, *atom))
        .count()
}

fn aromatic_bond_count(molecule: &Molecule) -> usize {
    molecule
        .bond_ids()
        .filter(|bond| aromatic_bond(molecule, *bond))
        .count()
}

#[test]
fn smiles_parse_options_bound_input_atoms_and_bonds() {
    let input_error = smiles_api::parse_str_with_options(
        "CC",
        SmilesParseOptions {
            max_input_bytes: 1,
            ..SmilesParseOptions::default()
        },
    )
    .expect_err("SMILES input byte limit should apply");
    assert!(input_error.message().contains("input"));

    let atom_error = smiles_api::parse_str_with_options(
        "CC",
        SmilesParseOptions {
            max_atoms: 1,
            ..SmilesParseOptions::default()
        },
    )
    .expect_err("SMILES atom limit should apply");
    assert!(atom_error.message().contains("atom count"));

    let bond_error = smiles_api::parse_str_with_options(
        "CC",
        SmilesParseOptions {
            max_bonds: 0,
            ..SmilesParseOptions::default()
        },
    )
    .expect_err("SMILES bond limit should apply");
    assert!(bond_error.message().contains("bond count"));
}

#[test]
fn smiles_record_preserves_name_extension_and_original_byte_spans() {
    let source = " \t[Na+].[Cl-]\t|$ion&#124;label;$,unknown:0|  sůl | sample\t ";
    let document = smiles_api::parse_str(source).unwrap();
    assert_eq!(document.source(), source);
    assert_eq!(document.base_smiles(), "[Na+].[Cl-]");
    assert_eq!(document.base_smiles_span(), 2..13);
    assert_eq!(
        document.cx_extension(),
        Some("|$ion&#124;label;$,unknown:0|")
    );
    assert_eq!(document.name(), Some("sůl | sample"));
    assert_eq!(&source[document.name_span().unwrap()], "sůl | sample");
    assert_eq!(
        &source[document.cx_extension_span().unwrap()],
        document.cx_extension().unwrap()
    );
    assert_eq!(
        document
            .tokens()
            .iter()
            .map(|token| &source[token.span()])
            .collect::<Vec<_>>(),
        ["[Na+]", ".", "[Cl-]"]
    );
    let error = document.interpret().unwrap_err();
    assert_eq!(
        error.offset(),
        document.cx_extension_span().unwrap().start + 1
    );
    assert!(error.message().contains("CXSMILES"));
    assert!(document.to_molecules().is_err());

    let projected = smiles_api::interpret_base(&document).unwrap();
    assert_eq!(projected, document.interpret_base().unwrap());
    assert_eq!(projected.name(), document.name());
    assert_eq!(projected.cx_extension(), document.cx_extension());
    assert_eq!(
        projected.omitted_cx_extension_span(),
        document.cx_extension_span()
    );
    assert_eq!(projected.components().len(), 2);
    for (index, component) in projected.components().iter().enumerate() {
        assert_eq!(component.molecule().atom_count(), 1);
        let mapping = &component.report().atom_mappings()[0];
        assert_eq!(mapping.source_index(), index);
        assert_eq!(mapping.atom(), AtomId::new(0));
        assert_eq!(&source[mapping.source_span()], ["[Na+]", "[Cl-]"][index]);
        assert!(!component.molecule().perception().has_valence());
    }
}

#[test]
fn smiles_record_metadata_does_not_change_source_atom_order_or_connectivity() {
    for (base, indices) in [
        ("C(.O)N", vec![vec![0, 2], vec![1]]),
        ("C1.C1", vec![vec![0, 1]]),
    ] {
        let source = format!("\t{base}  molécula\twith a name ");
        let document = smiles_api::parse_str(&source).unwrap();
        let interpreted = document.interpret().unwrap();
        assert_eq!(interpreted, document.interpret_base().unwrap());
        assert_eq!(interpreted.name(), Some("molécula\twith a name"));
        assert_eq!(interpreted.cx_extension(), None);
        assert_eq!(interpreted.omitted_cx_extension_span(), None);
        let bare = smiles_api::parse_str(base).unwrap().interpret().unwrap();
        for ((component, expected_indices), bare_component) in interpreted
            .components()
            .iter()
            .zip(indices)
            .zip(bare.components())
        {
            assert_eq!(component.molecule(), bare_component.molecule());
            let mappings = component.report().atom_mappings();
            assert_eq!(
                mappings
                    .iter()
                    .map(|mapping| mapping.source_index())
                    .collect::<Vec<_>>(),
                expected_indices
            );
            for (mapping, original) in mappings.iter().zip(bare_component.report().atom_mappings())
            {
                assert_eq!(
                    mapping.source_span(),
                    original.source_span().start + 1..original.source_span().end + 1
                );
            }
        }
    }
}

#[test]
fn smiles_full_interpretation_never_silently_drops_unsupported_cx_extensions() {
    for extension in ["|a:0|", "|(0,0,0)|", "|futureField:0|", "|r:0|"] {
        let source = format!("C {extension} name");
        let document = smiles_api::parse_str(&source).unwrap();
        assert!(document.interpret().is_err(), "{source}");
        assert!(smiles_api::to_molecules(&source).is_err(), "{source}");
        let projection = document.interpret_base().unwrap();
        assert_eq!(projection.cx_extension(), Some(extension));
        assert_eq!(
            projection.omitted_cx_extension_span(),
            document.cx_extension_span()
        );
    }
    for source in ["C ||", "C | \t | name"] {
        let document = smiles_api::parse_str(source).unwrap();
        assert_eq!(
            document.interpret().unwrap().omitted_cx_extension_span(),
            None
        );
        assert_eq!(
            document
                .interpret_base()
                .unwrap()
                .omitted_cx_extension_span(),
            document.cx_extension_span()
        );
    }
}

#[test]
fn cxsmiles_radicals_preserve_electron_occupancy_and_explicit_spin() {
    for (code, electrons, spin) in [
        (1, 1, None),
        (2, 2, None),
        (3, 2, Some(1)),
        (4, 2, Some(3)),
        (5, 3, None),
        (6, 3, Some(2)),
        (7, 3, Some(4)),
    ] {
        for base in ["C".to_owned(), format!("[CH{}]", 4 - electrons)] {
            let source = format!("{base} |^{code}:0| radical");
            let interpretation = smiles_api::parse_str(&source).unwrap().interpret().unwrap();
            assert_eq!(interpretation.omitted_cx_extension_span(), None);
            let mut molecule = interpretation.into_molecule().unwrap();
            let atom = molecule.atom(AtomId::new(0)).unwrap();
            assert_eq!(atom.radical, AtomRadical::new(electrons, spin), "{source}");
            assert!(!molecule.perception().has_valence());
            perceive(&mut molecule).unwrap();
            assert_eq!(
                molecule
                    .atom(AtomId::new(0))
                    .unwrap()
                    .hydrogens
                    .specified_count()
                    + molecule
                        .inferred_hydrogens(AtomId::new(0))
                        .unwrap()
                        .unwrap(),
                4 - electrons
            );
        }
    }
    let source = "C(.O)N |^1:2,^2:1|";
    let interpreted = smiles_api::parse_str(source).unwrap().interpret().unwrap();
    assert_eq!(
        interpreted.components()[0]
            .molecule()
            .atom(AtomId::new(1))
            .unwrap()
            .radical,
        AtomRadical::new(1, None)
    );
    assert_eq!(
        interpreted.components()[1]
            .molecule()
            .atom(AtomId::new(0))
            .unwrap()
            .radical,
        AtomRadical::new(2, None)
    );
    assert_eq!(
        interpreted.components()[0]
            .molecule()
            .atom(AtomId::new(0))
            .unwrap()
            .radical,
        None
    );
    let document = smiles_api::parse_str("C |^1:0|").unwrap();
    assert_eq!(
        document
            .interpret_base()
            .unwrap()
            .molecule()
            .unwrap()
            .atom(AtomId::new(0))
            .unwrap()
            .radical,
        None
    );
}

#[test]
fn cxsmiles_radicals_override_inference_without_overwriting_source_assertions() {
    for source in ["c1ccccc1 |^1:0|", "[c]1ccccc1 |^1:0|"] {
        let mut molecule = read_smiles(source).unwrap();
        perceive(&mut molecule).unwrap();
        assert_eq!(
            molecule.atom(AtomId::new(0)).unwrap().radical,
            AtomRadical::new(1, None)
        );
        assert_eq!(
            molecule.inferred_hydrogens(AtomId::new(0)).unwrap(),
            Some(0)
        );
        assert_eq!(
            molecule
                .atom(AtomId::new(0))
                .unwrap()
                .hydrogens
                .specified_count(),
            0
        );
        assert_eq!(aromatic_atom_count(&molecule), 6);
    }
    // Import and perception must not silently rewrite contradictory source
    // assertions. The RDKit-like fixed-H valence policy is not a general
    // electronic-state validator.
    let mut molecule = read_smiles("[CH4] |^1:0|").unwrap();
    assert_eq!(
        molecule.atom(AtomId::new(0)).unwrap().radical,
        AtomRadical::new(1, None)
    );
    assert_eq!(
        molecule.atom(AtomId::new(0)).unwrap().hydrogens,
        HydrogenDeclaration::Fixed(4)
    );
    molecule.perceive().unwrap();
    assert_eq!(
        molecule.atom(AtomId::new(0)).unwrap().radical,
        AtomRadical::new(1, None)
    );
}

#[test]
fn cxsmiles_enhanced_groups_preserve_members_and_legacy_relative_configuration() {
    let base = "F[C@H](Cl)[C@@H](Br)I";
    for (fields, kinds, centers) in [
        ("a:1,3", vec![StereoGroupKind::Absolute], vec![vec![1, 3]]),
        ("&17:1,3,r", vec![StereoGroupKind::And], vec![vec![1, 3]]),
        ("r,o17:1,3", vec![StereoGroupKind::Or], vec![vec![1, 3]]),
        ("&1:1,&1:3", vec![StereoGroupKind::And], vec![vec![1, 3]]),
        (
            "&1:1,o1:3",
            vec![StereoGroupKind::And, StereoGroupKind::Or],
            vec![vec![1], vec![3]],
        ),
        ("r", vec![StereoGroupKind::Relative], vec![vec![1, 3]]),
        (
            "a:1,r",
            vec![StereoGroupKind::Absolute, StereoGroupKind::Relative],
            vec![vec![1], vec![3]],
        ),
    ] {
        let source = format!("{base} |{fields}|");
        let molecule = read_smiles(&source).unwrap();
        assert!(!molecule.perception().has_valence());
        let groups = molecule.stereo_groups().collect::<Vec<_>>();
        assert_eq!(groups.len(), kinds.len(), "{source}");
        for (((group_id, group), kind), expected_centers) in
            groups.into_iter().zip(kinds).zip(centers)
        {
            assert_eq!(group.kind, kind, "{source}");
            let actual_centers = group
                .members
                .iter()
                .map(|id| {
                    let element = molecule.stereo_element(*id).unwrap();
                    assert_eq!(element.group, Some(group_id));
                    let StereoElementKind::Tetrahedral(stereo) = &element.kind else {
                        panic!("tetrahedral member")
                    };
                    stereo.center.raw()
                })
                .collect::<Vec<_>>();
            assert_eq!(actual_centers, expected_centers, "{source}");
        }
        let written =
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()).unwrap();
        let restored = read_smiles(&written).unwrap();
        assert_eq!(
            restored.stereo_groups().count(),
            molecule.stereo_groups().count()
        );
        assert!(
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::default()).is_err(),
            "ordinary mode rejects represented stereo"
        );
    }
    let molecule = read_smiles("F/C=C/F |r|").unwrap();
    assert_eq!(molecule.stereo_elements().count(), 1);
    assert_eq!(
        molecule.stereo_groups().count(),
        0,
        "relative chirality does not invert alkene geometry"
    );
}

#[test]
fn cxsmiles_stereo_references_use_source_indices_and_respect_component_ownership() {
    let source = "F(.O)[C@H:17](Cl)Br |o1:2|";
    let interpreted = smiles_api::parse_str(source).unwrap().interpret().unwrap();
    let molecule = interpreted.components()[0].molecule();
    let group = molecule.stereo_groups().next().unwrap().1;
    let StereoElementKind::Tetrahedral(stereo) =
        &molecule.stereo_element(group.members[0]).unwrap().kind
    else {
        panic!("tetrahedral member")
    };
    assert_eq!(stereo.center, AtomId::new(1));
    assert_eq!(molecule.atom(stereo.center).unwrap().atom_map, Some(17));
    assert_eq!(
        interpreted.components()[1]
            .molecule()
            .stereo_groups()
            .count(),
        0
    );

    let base = "F[C@H](Cl)Br.F[C@H](Cl)I";
    let document = smiles_api::parse_str(&format!("{base} |a:1,5|")).unwrap();
    let absolute = document.interpret().unwrap();
    assert!(absolute
        .molecules()
        .all(|molecule| molecule.stereo_groups().count() == 1));
    for field in ["&1:1,5", "o1:1,5", "r"] {
        let document = smiles_api::parse_str(&format!("{base} |{field}|")).unwrap();
        assert!(document
            .interpret()
            .unwrap_err()
            .message()
            .contains("disconnected"));
    }
    assert!(smiles_api::parse_str(&format!("{base} |&1:1,o1:5|"))
        .unwrap()
        .interpret()
        .is_ok());
}

#[test]
fn cxsmiles_rejects_ambiguous_invalid_and_unsupported_chemical_fields() {
    for fields in [
        "^0:1",
        "^8:1",
        "^1:",
        "^1:4",
        "^1:-1",
        "^1:1,1",
        "^1:1,^2:1",
        "a:1,1",
        "a:1,&2:1",
        "a:0",
        "&:1",
        "&999999999999999999999999999999:1",
        "&1:1,",
        "&1:1,,r",
        "r,r",
        "r:0",
        "&1:1,unknown:0",
        "w:0.0",
        "C:0.0",
        "(0,0,0)",
    ] {
        let source = format!("F[C@H](Cl)Br |{fields}|");
        let document = smiles_api::parse_str(&source).unwrap();
        let error = document.interpret().unwrap_err();
        let span = document.cx_extension_span().unwrap();
        assert!(span.contains(&error.offset()), "{source}: {error}");
        assert!(error.message().contains("CXSMILES"), "{source}: {error}");
        assert!(
            document.interpret_base().is_ok(),
            "explicit projection retains original fields"
        );
    }
}

#[test]
fn smiles_record_rejects_malformed_framing_and_multiple_records() {
    for (source, offset) in [
        ("C |r", 2),
        ("C |r|name", 5),
        ("C ||||", 4),
        (" \t", 2),
        ("C\nO", 1),
        ("C |$line\nbreak$|", 8),
        ("C name\rnext", 6),
    ] {
        let error = smiles_api::parse_str(source).unwrap_err();
        assert_eq!(error.offset(), offset, "{source:?}: {error}");
    }
    assert!(smiles_api::parse_str("C|r|").is_err());
    let named = smiles_api::parse_str("C a name |with literal bars|").unwrap();
    assert_eq!(named.cx_extension(), None);
    assert_eq!(named.name(), Some("a name |with literal bars|"));
}

#[test]
fn smiles_record_limits_cover_metadata_and_syntax_errors_point_into_base() {
    let source = "\tCC |unknown:999| molecule";
    let options = SmilesParseOptions {
        max_input_bytes: source.len(),
        max_atoms: 2,
        max_bonds: 1,
    };
    assert!(smiles_api::parse_str_with_options(source, options).is_ok());
    for (options, message) in [
        (
            SmilesParseOptions {
                max_input_bytes: source.len() - 1,
                ..options
            },
            "byte limit",
        ),
        (
            SmilesParseOptions {
                max_atoms: 1,
                ..options
            },
            "atom count",
        ),
        (
            SmilesParseOptions {
                max_bonds: 0,
                ..options
            },
            "bond count",
        ),
    ] {
        assert!(smiles_api::parse_str_with_options(source, options)
            .unwrap_err()
            .message()
            .contains(message));
    }
    assert_eq!(smiles_api::parse_str("\tC( title").unwrap_err().offset(), 3);
    assert_eq!(
        smiles_api::parse_str("\t[CH3 title").unwrap_err().offset(),
        1
    );
    assert_eq!(smiles_api::parse_str("\tC. title").unwrap_err().offset(), 3);
}

#[test]
fn smiles_document_preserves_spans_and_dot_boundaries_before_interpretation() {
    let input = "[Na+].[Cl-]";
    let document = smiles_api::parse_str(input).expect("document parses");
    assert_eq!(document.source(), input);
    assert_eq!(document.fragment_token_ranges().len(), 2);
    assert!(document.tokens().iter().all(|token| {
        let span = token.span();
        span.start <= span.end && span.end <= input.len()
    }));
    let interpretation = smiles_api::interpret(&document).expect("document interprets");
    let error = interpretation
        .molecule()
        .expect_err("multi-component SMILES is not one molecule");
    assert_eq!(error.actual(), 2);
    assert_eq!(interpretation.components().len(), 2);
    for component in interpretation.components() {
        assert_eq!(component.molecule().atom_count(), 1);
        assert_eq!(component.molecule().bond_count(), 0);
        component
            .molecule()
            .validate_connected()
            .expect("each interpreted component is connected");
        assert!(!component.molecule().perception().has_valence());
        assert_eq!(component.report().atom_mappings().len(), 1);
        assert!(component.report().bond_mappings().is_empty());
    }
}

#[test]
fn smiles_interprets_branches_rings_brackets_and_fragments_canonically_without_perceiving() {
    let components = read_smiles_components("C(C)O.C1=CC=CC=C1.[13NH4+:7].F[C@@H](N)O")
        .expect("smiles should parse");

    assert_eq!(components.len(), 4);
    assert_eq!(
        components
            .iter()
            .map(|molecule| molecule.atom_count())
            .sum::<usize>(),
        14
    );
    assert_eq!(
        components
            .iter()
            .map(|molecule| molecule.bond_count())
            .sum::<usize>(),
        11
    );
    for molecule in &components {
        assert_all_stale(molecule);
    }
    let bracket_atom = components[2].atom(AtomId::new(0)).expect("bracket atom");
    assert_eq!(bracket_atom.isotope, Some(13));
    assert_eq!(bracket_atom.hydrogens, HydrogenDeclaration::Fixed(4));
    assert_eq!(bracket_atom.formal_charge, 1);
    assert_eq!(bracket_atom.atom_map, Some(7));
    let chiral_atom = components[3]
        .atom(AtomId::new(1))
        .expect("chiral bracket atom");
    assert_eq!(chiral_atom.hydrogens, HydrogenDeclaration::Fixed(1));
    let stereo = components[3]
        .stereo_elements()
        .map(|(_, element)| element)
        .collect::<Vec<_>>();
    assert_eq!(stereo.len(), 1);
    match &stereo[0].kind {
        StereoElementKind::Tetrahedral(tetrahedral) => {
            assert_eq!(tetrahedral.center, AtomId::new(1));
            assert_eq!(
                tetrahedral.orientation,
                Some(TetrahedralOrientation::CounterClockwise)
            );
            assert!(tetrahedral
                .carriers
                .contains(&StereoCarrier::ImplicitHydrogen));
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn smiles_brackets_publish_exact_hydrogen_declarations() {
    for (source, declaration, perceived_implicit) in [
        ("C", HydrogenDeclaration::Infer { specified: 0 }, 4),
        ("[C]", HydrogenDeclaration::Fixed(0), 0),
        ("[CH]", HydrogenDeclaration::Fixed(1), 0),
        ("[NH4+]", HydrogenDeclaration::Fixed(4), 0),
    ] {
        let mut molecule = read_smiles(source)
            .unwrap_or_else(|error| panic!("{source} should interpret: {error}"));
        let atom = molecule.atom(AtomId::new(0)).expect("source atom");
        assert_eq!(atom.hydrogens, declaration, "{source}");
        assert!(!molecule.perception().has_valence(), "{source}");

        perceive(&mut molecule).unwrap_or_else(|error| panic!("{source} should perceive: {error}"));
        assert_eq!(
            molecule
                .atom(AtomId::new(0))
                .expect("perceived atom")
                .hydrogens,
            declaration,
            "perception must not rewrite {source}"
        );
        assert_eq!(
            molecule.inferred_hydrogens(AtomId::new(0)),
            Ok(Some(perceived_implicit)),
            "{source}"
        );
    }
}

#[test]
fn metal_bound_organic_subset_halogen_keeps_rdkit_no_implicit_state() {
    let mut small = read_smiles("Br[Pt+2]Br").expect("platinum bromide salt parses");
    perceive(&mut small).expect("platinum bromide salt perceives");

    let bromines = small
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "Br")
        .map(|(atom_id, atom)| {
            (
                !atom.hydrogens.allows_inference(),
                inferred_hydrogens(&small, atom_id).unwrap_or(0),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(bromines, vec![(false, 0), (false, 0)]);

    let mut aryl_bromide = read_smiles("c1ccccc1Br").expect("aryl bromide should parse");
    perceive(&mut aryl_bromide).expect("aryl bromide should perceive");
    let bromine = aryl_bromide
        .atoms()
        .find_map(|(_, atom)| (atom.element.symbol() == "Br").then_some(atom))
        .expect("bromine atom");
    assert!(bromine.hydrogens.allows_inference());
}

#[test]
fn metal_bound_organic_subset_atoms_rely_on_valence_hydrogens() {
    let mut aryl_mercury = read_smiles("c1ccccc1[Hg]").expect("aryl mercury should parse");
    perceive(&mut aryl_mercury).expect("aryl mercury should perceive");
    let aryl_mercury_carbon = aryl_mercury
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "C"
                && aryl_mercury.incident_bonds(id).is_ok_and(|bonds| {
                    bonds.into_iter().any(|(_, bond)| {
                        aryl_mercury
                            .atom(bond.other_atom(id))
                            .is_ok_and(|neighbor| neighbor.element.symbol() == "Hg")
                    })
                }))
            .then_some((id, atom))
        })
        .expect("aryl carbon bound to mercury");
    assert!(aryl_mercury_carbon.1.hydrogens.allows_inference());
    assert_eq!(
        inferred_hydrogens(&aryl_mercury, aryl_mercury_carbon.0),
        Some(0)
    );

    let methyl_sodium = read_smiles("C[Na]").expect("methyl sodium should parse");
    let carbon = methyl_sodium
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "C").then_some((id, atom)))
        .expect("carbon atom");
    assert!(carbon.1.hydrogens.allows_inference());
    assert_eq!(inferred_hydrogens(&methyl_sodium, carbon.0), None);
}

#[test]
fn aromatic_chalcogen_bracket_atoms_localize_without_perceiving() {
    let components = read_smiles_components("[se]1cccc1.[te]1cccc1")
        .expect("aromatic selenium and tellurium bracket atoms should parse");

    let chalcogens = components
        .iter()
        .flat_map(|molecule| molecule.atoms())
        .filter(|(_, atom)| matches!(atom.element.symbol(), "Se" | "Te"))
        .map(|(_, atom)| {
            (
                atom.element.symbol().to_owned(),
                !atom.hydrogens.allows_inference(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        chalcogens,
        vec![("Se".to_owned(), true), ("Te".to_owned(), true)]
    );
    assert!(components
        .iter()
        .all(|molecule| !molecule.perception().has_aromaticity()));
    assert!(components.iter().all(|molecule| molecule
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double))));
}

#[test]
fn malformed_smiles_grammar_is_rejected_by_the_document_parser() {
    for (input, offset) in [
        ("C(", 2),
        ("C1", 2),
        ("C%1", 1),
        ("C%a1", 1),
        ("C=", 1),
        ("=C", 0),
        ("C..C", 2),
        ("C=1CCCCC-1", 9),
        ("[]", 0),
        ("[13]", 0),
        ("[é]", 0),
        ("[C@@@H]", 4),
        ("[C/]", 2),
        ("[*]", 1),
        ("[C+999]", 6),
        ("[C:]", 3),
        ("[Clx]", 3),
        ("[si]1ccccc1", 2),
        ("Cé", 1),
    ] {
        let parsed = std::panic::catch_unwind(|| smiles_api::parse_str(input))
            .unwrap_or_else(|_| panic!("`{input}` panicked"));
        let error = parsed.expect_err("malformed SMILES must fail document parsing");
        assert_eq!(error.offset(), offset, "{input}: {error}");
    }
    assert_eq!(
        smiles_api::parse_str("C(").unwrap_err().to_string(),
        "SMILES parse error at 2: unclosed branch"
    );
}

#[test]
fn canonical_smiles_preserves_stereo_and_isotopes() {
    for source in [
        "N[C@H](O)C",
        "[11CH3]OC",
        "C1=CC=[14CH]C=C1",
        "[H]C([3H])(F)Cl",
        "[H][C@](F)(Cl)Br",
        "C/C=C\\C",
    ] {
        let mut molecule = read_smiles(source).unwrap();
        perceive(&mut molecule).unwrap();
        let before = assigned_descriptors(&mut molecule);
        let isotopes = molecule
            .atoms()
            .filter_map(|(_, atom)| atom.isotope)
            .collect::<Vec<_>>();
        let written =
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()).unwrap();
        let mut restored = read_smiles(&written).unwrap();
        perceive(&mut restored).unwrap();
        assert_eq!(
            assigned_descriptors(&mut restored),
            before,
            "{source} -> {written}"
        );
        assert_eq!(
            restored
                .atoms()
                .filter_map(|(_, atom)| atom.isotope)
                .collect::<Vec<_>>(),
            isotopes
        );
        assert_eq!(
            smiles_api::write(&restored, smiles_api::SmilesWriteOptions::canonical()).unwrap(),
            written
        );
    }
    let canonical = |source| {
        let mut mol = read_smiles(source).unwrap();
        perceive(&mut mol).unwrap();
        smiles_api::write(&mol, smiles_api::SmilesWriteOptions::canonical()).unwrap()
    };
    assert_ne!(canonical("N[C@H](O)C"), canonical("N[C@@H](O)C"));
    assert_ne!(canonical("[11CH3]OC"), canonical("COC"));
}

#[test]
fn source_order_smiles_retains_branch_order_and_omits_kekule_single_bonds() {
    // Independently checked with RDKit 2026.03.3, canonical=False.
    for source in [
        "CC(CCC)O",
        "N[C@](CCC)(F)Cl",
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "C1(O)C(O)C(O)C(OP(=O)(O)O)C(O)C1O",
        "C12(CCCCC1)CCCCC2",
    ] {
        let mut molecule = read_smiles(source).unwrap();
        perceive(&mut molecule).unwrap();
        let written =
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()).unwrap();
        assert_eq!(written, source);
        let mut round_trip = read_smiles(&written).unwrap();
        perceive(&mut round_trip).unwrap();
        assert_eq!(
            smiles_api::write(&round_trip, smiles_api::SmilesWriteOptions::canonical()).unwrap(),
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()).unwrap()
        );
    }
    let mut branched = read_smiles("CC(CCC)O").unwrap();
    perceive(&mut branched).unwrap();
    assert_eq!(
        smiles_api::write(&branched, smiles_api::SmilesWriteOptions::default()).unwrap(),
        "CC(CCC)O"
    );
    // Aromatic spelling still needs the explicit single bond between rings.
    let mut biphenyl = read_smiles("c1ccccc1-c1ccccc1").unwrap();
    perceive(&mut biphenyl).unwrap();
    assert!(
        smiles_api::write(&biphenyl, smiles_api::SmilesWriteOptions::default())
            .unwrap()
            .contains('-')
    );
    assert!(
        !smiles_api::write(&biphenyl, smiles_api::SmilesWriteOptions::isomeric())
            .unwrap()
            .contains('-')
    );
}

#[test]
fn canonical_smiles_prefers_clean_simple_ring_closure() {
    let molecule = read_smiles("C1=CC=CC=C1").expect("benzene parses");

    let written = smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical())
        .expect("canonical SMILES should write");

    assert_eq!(written, "C1=CC=CC=C1");
}

#[test]
fn canonical_smiles_converges_after_aromaticity_perception() {
    let mut aromatic = read_smiles("c1ccccc1").expect("aromatic benzene parses");
    let mut kekule = read_smiles("C1=CC=CC=C1").expect("Kekule benzene parses");
    perceive(&mut aromatic).expect("aromatic benzene perceives");
    perceive(&mut kekule).expect("Kekule benzene perceives");

    let aromatic_written =
        smiles_api::write(&aromatic, smiles_api::SmilesWriteOptions::canonical())
            .expect("aromatic benzene canonicalizes");
    let kekule_written = smiles_api::write(&kekule, smiles_api::SmilesWriteOptions::canonical())
        .expect("perceived Kekule benzene canonicalizes");

    assert_eq!(aromatic_written, kekule_written);
    assert_eq!(aromatic_written, "c1ccccc1");
}

#[test]
fn canonical_smiles_preserves_aromatic_high_order_bonds() {
    let mut molecule = read_smiles("C1=CC#CC=C1").expect("cyclohexyne parses");
    perceive(&mut molecule).expect("cyclohexyne perceives");

    let (written, reparsed) = canonical_smiles_round_trip(&molecule);
    assert!(written.contains('#'), "{written}");
    assert!(reparsed.bonds().any(
        |(bond_id, bond)| bond.order == BondOrder::Triple && aromatic_bond(&reparsed, bond_id)
    ));
}

#[test]
fn aromatic_smiles_omitted_bonds_perceive_with_expected_hydrogens() {
    let mut benzene = read_smiles("c1ccccc1").expect("benzene should parse");
    assert_eq!(
        benzene
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        3
    );
    perceive(&mut benzene).expect("benzene should perceive");
    for atom_id in benzene.atom_ids() {
        assert_eq!(inferred_hydrogens(&benzene, atom_id), Some(1));
        assert!(aromatic_atom(&benzene, atom_id));
    }

    let mut pyridine = read_smiles("n1ccccc1").expect("pyridine should parse");
    perceive(&mut pyridine).expect("pyridine should perceive");
    assert_eq!(inferred_hydrogens(&pyridine, AtomId::new(0)), Some(0));
    for atom_id in 1..6 {
        assert_eq!(inferred_hydrogens(&pyridine, AtomId::new(atom_id)), Some(1));
    }

    let mut pyridinium = read_smiles("[nH+]1ccccc1").expect("pyridinium should parse");
    perceive(&mut pyridinium).expect("pyridinium should perceive");
    let nitrogen = pyridinium.atom(AtomId::new(0)).expect("nitrogen");
    assert!(aromatic_atom(&pyridinium, AtomId::new(0)));
    assert_eq!(nitrogen.formal_charge, 1);
    assert_eq!(nitrogen.radical, None);
    assert_eq!(nitrogen.hydrogens, HydrogenDeclaration::Fixed(1));
    assert_eq!(inferred_hydrogens(&pyridinium, AtomId::new(0)), Some(0));
    assert_eq!(aromatic_bond_count(&pyridinium), pyridinium.bond_count());
    assert_eq!(
        pyridinium
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
        let mut components = read_smiles_components(smiles)
            .unwrap_or_else(|_| panic!("supported aromatic SMILES should parse: {smiles}"));
        for molecule in &mut components {
            perceive(molecule).unwrap_or_else(|_| {
                panic!("supported aromatic component should perceive: {smiles}")
            });
            let written = smiles_api::write(&*molecule, smiles_api::SmilesWriteOptions::default())
                .unwrap_or_else(|_| panic!("supported aromatic component should write: {smiles}"));
            read_smiles(&written)
                .unwrap_or_else(|_| panic!("writer output should parse: {written}"));
        }
    }
}

#[test]
fn invalid_lowercase_aromatic_ring_returns_structured_error() {
    for smiles in ["c1cccc1", "c1ccccc1.c1cccc1"] {
        let error = read_smiles_components(smiles)
            .expect_err("unlocalizable aromatic source must fail interpretation");
        assert!(error
            .to_string()
            .contains("invalid imported aromatic representation"));
    }
}

#[test]
fn smiles_interpretation_rejects_inconsistent_source_aromaticity() {
    let cases = [
        (
            "c",
            0,
            "source-aromatic atom is not part of a source-aromatic bond",
        ),
        (
            "c-C",
            0,
            "source-aromatic atom is not part of a source-aromatic bond",
        ),
        (
            "c-c",
            0,
            "source-aromatic atom is not part of a source-aromatic bond",
        ),
        (
            "c=C",
            0,
            "source-aromatic atom is not part of a source-aromatic bond",
        ),
        (
            "C:C",
            1,
            "source-aromatic bond requires source-aromatic atom syntax at both endpoints",
        ),
        (
            "c:C",
            1,
            "source-aromatic bond requires source-aromatic atom syntax at both endpoints",
        ),
        (
            "C:c",
            1,
            "source-aromatic bond requires source-aromatic atom syntax at both endpoints",
        ),
    ];

    for (source, offset, message) in cases {
        let document = smiles_api::parse_str(source)
            .unwrap_or_else(|_| panic!("source syntax should parse: {source}"));
        let error = match smiles_api::interpret(&document) {
            Ok(_) => panic!("inconsistent source aromaticity should fail: {source}"),
            Err(error) => error,
        };

        assert_eq!(error.offset(), offset, "{source}: {error}");
        assert_eq!(error.message(), message, "{source}: {error}");
    }
}

#[test]
fn smiles_source_aromaticity_validation_preserves_supported_forms() {
    for source in [
        "cc",
        "c:c",
        "c1:c:c:c:c:c:1",
        "n1ccccc1",
        "c1ccoc1",
        "c1ccc2ccccc2c1",
        "c1ccccc1-C",
        "c1ccccc1-c1ccccc1",
    ] {
        let molecule = read_smiles(source)
            .unwrap_or_else(|_| panic!("consistent source aromaticity should publish: {source}"));

        assert!(molecule
            .bonds()
            .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
        assert!(!molecule.perception().has_aromaticity());
    }
}

#[test]
fn canonical_smiles_brackets_atoms_whose_hydrogens_are_not_inferred() {
    // Atoms bonded to metals and hypervalent or heavy main-group atoms keep
    // explicit hydrogen counts. Every expected string matches RDKit 2026.09.1
    // canonical SMILES for the same input.
    for (source, expected) in [
        ("c1ccccc1[Hg]", "[Hg][c]1ccccc1"),
        (
            "C1=CC=C(C(=C1)[N+](=O)[O-])[Hg]",
            "O=[N+]([O-])c1cccc[c]1[Hg]",
        ),
        ("CC[Hg+]", "C[CH2][Hg+]"),
        ("C[Tl](C)C", "[CH3][Tl]([CH3])[CH3]"),
        ("C[Sb](C)C", "[CH3][Sb]([CH3])[CH3]"),
        ("C[Na]", "[CH3][Na]"),
        ("Br[Pt+2]Br", "[Br][Pt+2][Br]"),
        ("Cl[Cr]Cl", "[Cl][Cr][Cl]"),
        ("OP(=O)O", "O=[PH](O)O"),
        (
            "C1=CC=C(C=C1)[Ge](Cl)(Cl)Cl",
            "[Cl][Ge]([Cl])([Cl])[c]1ccccc1",
        ),
        (
            "C1=CC=C(C=C1)[SnH](C2=CC=CC=C2)Cl",
            "[Cl][SnH]([c]1ccccc1)[c]1ccccc1",
        ),
        ("C1=C[Te]C=C1", "c1cc[te]c1"),
        (
            "CCOC(=O)C1=C(C(=C(N1)C)C(=O)OC(C)(C)C)C",
            "CCOC(=O)c1[nH]c(C)c(C(=O)OC(C)(C)C)c1C",
        ),
    ] {
        let mut molecule = read_smiles(source).unwrap();
        perceive(&mut molecule).unwrap();
        let (written, _) = canonical_smiles_round_trip(&molecule);
        assert_eq!(written, expected, "{source}");
    }
}

#[test]
fn smiles_writer_rejects_lossy_bonds_and_stereo() {
    let mut molecule = crate::core::MoleculeEditor::new();
    let a = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let b = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let bond = molecule.add_bond(a, b, BondOrder::Dative).expect("bond");
    assert!(smiles_api::write(
        molecule.working(),
        smiles_api::SmilesWriteOptions::default()
    )
    .expect_err("dative bond should be rejected")
    .message
    .contains("cannot encode"));

    molecule
        .bond_mut(bond)
        .expect("bond")
        .set_order(BondOrder::Single);
    let c = molecule.add_atom(carbon()).expect("third atom");
    let d = molecule.add_atom(carbon()).expect("fourth atom");
    molecule.add_bond(a, c, BondOrder::Single).expect("bond");
    molecule.add_bond(a, d, BondOrder::Single).expect("bond");
    molecule
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center: a,
                carriers: vec![
                    StereoCarrier::Atom(b),
                    StereoCarrier::Atom(c),
                    StereoCarrier::Atom(d),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("atom stereo");
    assert!(smiles_api::write(
        molecule.working(),
        smiles_api::SmilesWriteOptions::default()
    )
    .expect_err("atom chirality should be rejected")
    .message
    .contains("stereochemistry"));

    let element = molecule
        .stereo_element_ids()
        .next()
        .expect("stereo element");
    molecule
        .remove_stereo_element(element)
        .expect("remove atom stereo");
    {
        let mut atom = molecule.atom_mut(a).expect("atom");
        atom.radical = AtomRadical::new(1, Some(2));
        atom.hydrogens = HydrogenDeclaration::Fixed(2);
    }
    assert!(smiles_api::write(
        molecule.working(),
        smiles_api::SmilesWriteOptions::default()
    )
    .expect_err("ordinary SMILES cannot preserve an explicit spin assertion")
    .message
    .contains("explicit radical spin"));

    {
        let mut atom = molecule.atom_mut(a).expect("atom");
        atom.radical = AtomRadical::new(1, None);
        atom.hydrogens = HydrogenDeclaration::Fixed(0);
    }
    let written = smiles_api::write(
        molecule.working(),
        smiles_api::SmilesWriteOptions::default(),
    )
    .expect("no-implicit-hydrogen atom should write");
    assert!(written.contains("[C]"));
    let reparsed = read_smiles(&written).expect("writer output should parse");
    assert!(reparsed
        .atoms()
        .any(|(_, atom)| !atom.hydrogens.allows_inference()));
}

#[test]
fn all_smiles_writers_round_trip_lossless_hydrogen_declarations() {
    for (source, expected) in [
        ("C", HydrogenDeclaration::Infer { specified: 0 }),
        ("[C]", HydrogenDeclaration::Fixed(0)),
        ("[CH]", HydrogenDeclaration::Fixed(1)),
        ("[NH4+]", HydrogenDeclaration::Fixed(4)),
    ] {
        let molecule = read_smiles(source).unwrap_or_else(|error| panic!("{source}: {error}"));
        for (writer, written) in [
            (
                "regular",
                smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::default()),
            ),
            (
                "canonical",
                smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()),
            ),
            (
                "isomeric",
                smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()),
            ),
        ] {
            let written = written
                .unwrap_or_else(|error| panic!("{writer} writer rejected {source}: {error}"));
            let reparsed = read_smiles(&written).unwrap_or_else(|error| {
                panic!("{writer} output for {source} did not parse ({written}): {error}")
            });
            assert_eq!(
                reparsed
                    .atom(AtomId::new(0))
                    .expect("single atom")
                    .hydrogens,
                expected,
                "{writer} writer changed {source} through {written}"
            );
        }
    }
}

#[test]
fn all_smiles_writers_require_known_total_for_declared_and_inferred_hydrogens() {
    let mut atom = carbon();
    atom.hydrogens = HydrogenDeclaration::Infer { specified: 1 };
    let mut graph = crate::core::MoleculeEditor::new();
    graph.add_atom(atom).expect("carbon");
    let molecule = graph.finish().expect("single atom molecule");

    for (writer, result) in [
        (
            "regular",
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::default()),
        ),
        (
            "canonical",
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()),
        ),
        (
            "isomeric",
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()),
        ),
    ] {
        let error = match result {
            Err(error) => error,
            Ok(written) => panic!(
                "{writer} writer must not coerce Infer {{ explicit: 1 }} to Fixed(1): {written}"
            ),
        };
        assert!(
            error.message().contains("hydrogen perception"),
            "{writer}: {error}"
        );
    }
}

#[test]
fn all_smiles_writers_preserve_total_declared_and_inferred_hydrogens() {
    for (symbol, declared, total) in [("C", 1, 4), ("C", 4, 4), ("N", 1, 3), ("O", 1, 2)] {
        let mut atom = Atom::new(Element::from_symbol(symbol).unwrap());
        atom.hydrogens = HydrogenDeclaration::Infer {
            specified: declared,
        };
        let mut editor = MoleculeEditor::new();
        let id = editor.add_atom(atom).unwrap();
        let mut molecule = editor.finish().unwrap();
        perceive(&mut molecule).unwrap();
        assert_eq!(molecule.inferred_hydrogens(id), Ok(Some(total - declared)));
        let before = molecule.clone();
        for written in [
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::default()),
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()),
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()),
        ] {
            let written = written.unwrap();
            let mut reparsed = read_smiles(&written).unwrap();
            perceive(&mut reparsed).unwrap();
            let (id, atom) = reparsed.atoms().next().unwrap();
            assert_eq!(
                atom.hydrogens.specified_count()
                    + reparsed.inferred_hydrogens(id).unwrap().unwrap(),
                total,
                "{symbol}: {written}"
            );
            assert_eq!(atom.element.symbol(), symbol);
            assert_eq!(atom.formal_charge, 0);
            assert_eq!(atom.radical, None);
            assert_eq!(molecule, before);
            assert_eq!(molecule.perception(), before.perception());
        }
    }
}

#[test]
fn bracket_atoms_interpret_radical_electrons_without_asserting_spin() {
    for (smiles, atom_index, electrons) in [
        ("[C]", 0, 4),
        ("[C]C", 0, 3),
        ("C=[C]", 1, 2),
        ("C#[C]", 1, 1),
        ("[N]", 0, 3),
        ("[O]", 0, 2),
        ("[c]1ccccc1", 0, 1),
        ("[CH4]", 0, 0),
        ("[NH4+]", 0, 0),
        ("[CH2-]", 0, 1),
        ("[Co]", 0, 1),
        ("[Co+]", 0, 0),
        ("[Co]C", 0, 0),
        ("[H]", 0, 1),
        ("[He]", 0, 0),
        ("[BH2]", 0, 1),
        ("[BH4-]", 0, 0),
        ("[F]", 0, 1),
        ("[Cl]", 0, 1),
        ("[PH4]", 0, 1),
        ("[SH3]", 0, 1),
        ("[SH]", 0, 1),
        ("[Na+]", 0, 0),
        ("[CoH]", 0, 1),
        ("[Co+127]", 0, 0),
        ("[C+127]", 0, 131),
        ("[C-128]", 0, 0),
    ] {
        let molecule = read_smiles(smiles).expect("bracket SMILES should interpret");
        assert_eq!(
            molecule
                .atom(AtomId::new(atom_index))
                .expect("bracket atom")
                .radical,
            AtomRadical::new(electrons, None),
            "{smiles}"
        );
    }

    let document = smiles_api::parse_str("[Xx]").expect("element spelling is valid syntax");
    assert!(smiles_api::interpret(&document)
        .expect_err("unsupported core element belongs to interpretation")
        .message()
        .contains("unsupported element"));
}

#[test]
fn smiles_writers_round_trip_bracket_radical_electrons_and_unspecified_spin() {
    for source in [
        "[C]", "[CH]", "[CH2]", "[CH3]", "[13CH3]", "[CH2-]", "[NH2]", "[OH]", "[Co]",
    ] {
        let molecule = read_smiles(source).unwrap();
        let before = molecule.clone();
        let original = molecule.atoms().next().unwrap();
        assert!(original.1.radical.is_some(), "{source}");
        for written in [
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::default()),
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()),
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()),
        ] {
            let written = written.unwrap();
            let reparsed = read_smiles(&written).unwrap();
            let atom = reparsed.atoms().next().unwrap();
            assert_eq!(atom, original, "{source} -> {written}");
            assert_eq!(molecule, before);
            assert_eq!(molecule.perception(), before.perception());
        }
    }
}

#[test]
fn bracket_radicals_are_consistent_across_aromatic_and_localized_notation() {
    for (aromatic, localized) in [
        ("[c]1ccccc1", "[C]1=CC=CC=C1"),
        ("[c]1ccncc1", "[C]1=CC=NC=C1"),
    ] {
        let mut canonical = Vec::new();
        for source in [aromatic, localized] {
            let mut molecule = read_smiles(source).unwrap();
            perceive(&mut molecule).unwrap();
            assert_eq!(
                molecule.atom(AtomId::new(0)).unwrap().radical,
                AtomRadical::new(1, None),
                "{source}"
            );
            let expected =
                smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()).unwrap();
            for written in [
                smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::default()).unwrap(),
                smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()).unwrap(),
                expected.clone(),
            ] {
                let mut restored = read_smiles(&written).unwrap();
                perceive(&mut restored).unwrap();
                assert_eq!(
                    smiles_api::write(&restored, smiles_api::SmilesWriteOptions::canonical())
                        .unwrap(),
                    expected,
                    "{source} -> {written}"
                );
                let radicals = restored
                    .atoms()
                    .filter_map(|(_, atom)| atom.radical)
                    .collect::<Vec<_>>();
                assert_eq!(radicals, vec![AtomRadical::new(1, None).unwrap()]);
            }
            canonical.push(expected);
        }
        assert_eq!(canonical[0], canonical[1]);
    }
}

#[test]
fn smiles_writers_reject_brackets_that_would_change_radical_occupancy() {
    for (hydrogens, radical) in [(0, None), (2, AtomRadical::new(1, None))] {
        let mut atom = carbon();
        atom.hydrogens = HydrogenDeclaration::Fixed(hydrogens);
        atom.radical = radical;
        let mut editor = MoleculeEditor::new();
        editor.add_atom(atom).unwrap();
        let molecule = editor.finish().unwrap();
        for result in [
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::default()),
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()),
            smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()),
        ] {
            assert!(result
                .unwrap_err()
                .message()
                .contains("radical electron count or spin"));
        }
    }
}

#[test]
fn isomeric_smiles_writes_tetrahedral_elements_from_stereo_model() {
    let molecule = read_smiles("F[C@H](Cl)Br").expect("tetrahedral SMILES should parse");
    assert_eq!(
        molecule
            .atom(AtomId::new(1))
            .expect("stereo center")
            .hydrogens,
        HydrogenDeclaration::Fixed(1)
    );

    let written = smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric())
        .expect("tetrahedral stereo should write");

    assert_eq!(written, "F[C@H](Cl)Br");
    let reparsed = read_smiles(&written).expect("isomeric output should parse");
    assert_eq!(
        reparsed
            .atom(AtomId::new(1))
            .expect("reparsed stereo center")
            .hydrogens,
        HydrogenDeclaration::Fixed(1)
    );
    let stereo = reparsed
        .stereo_elements()
        .map(|(_, element)| element)
        .collect::<Vec<_>>();
    assert_eq!(stereo.len(), 1);
    match &stereo[0].kind {
        StereoElementKind::Tetrahedral(tetrahedral) => {
            assert_eq!(tetrahedral.center, AtomId::new(1));
            assert_eq!(
                tetrahedral.orientation,
                Some(TetrahedralOrientation::Clockwise)
            );
            assert_eq!(
                tetrahedral.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(0)),
                    StereoCarrier::Atom(AtomId::new(2)),
                    StereoCarrier::Atom(AtomId::new(3)),
                    StereoCarrier::ImplicitHydrogen,
                ]
            );
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn isomeric_smiles_materializes_required_tetrahedral_hydrogen_without_mutating_source() {
    let mut molecule = read_smiles("FC(Cl)Br").expect("tetrahedral graph should parse");
    perceive(&mut molecule).expect("tetrahedral graph should perceive");
    let center = AtomId::new(1);
    assert_eq!(
        molecule.atom(center).expect("center").hydrogens,
        HydrogenDeclaration::Infer { specified: 0 }
    );
    assert_eq!(molecule.inferred_hydrogens(center), Ok(Some(1)));
    molecule
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(AtomId::new(0)),
                    StereoCarrier::Atom(AtomId::new(2)),
                    StereoCarrier::Atom(AtomId::new(3)),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("canonical tetrahedral element");

    let before = assigned_descriptors(&mut molecule);
    for written in [
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()).unwrap(),
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()).unwrap(),
    ] {
        let mut restored = read_smiles(&written).unwrap();
        perceive(&mut restored).unwrap();
        assert_eq!(assigned_descriptors(&mut restored), before);
        assert!(written.contains("@H"), "{written}");
    }
    assert_eq!(
        molecule.atom(center).expect("center").hydrogens,
        HydrogenDeclaration::Infer { specified: 0 }
    );
}

#[test]
fn isomeric_smiles_flips_tetrahedral_marker_for_odd_writer_carrier_order() {
    let mut molecule = read_smiles("F[C@H](Cl)Br").expect("tetrahedral SMILES should parse");
    let element = molecule
        .stereo_element_ids()
        .next()
        .expect("stereo element");
    molecule
        .remove_stereo_element(element)
        .expect("remove parsed stereo");
    molecule
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center: AtomId::new(1),
                carriers: vec![
                    StereoCarrier::Atom(AtomId::new(0)),
                    StereoCarrier::ImplicitHydrogen,
                    StereoCarrier::Atom(AtomId::new(3)),
                    StereoCarrier::Atom(AtomId::new(2)),
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("replacement stereo element");

    let written = smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric())
        .expect("tetrahedral stereo should write");

    assert_eq!(written, "F[C@@H](Cl)Br");
}

#[test]
fn isomeric_smiles_accepts_interpreted_source_stereo_and_rejects_unknown_stereo() {
    let directional = read_smiles("C/C=C\\C").expect("directional bond markers should parse");
    let written = smiles_api::write(&directional, smiles_api::SmilesWriteOptions::isomeric())
        .expect("interpreted directional stereo should write without perception");
    assert!(written.contains('/') || written.contains('\\'));

    let mut unknown = read_smiles("F[C@H](Cl)Br").expect("tetrahedral SMILES should parse");
    let element = unknown.stereo_element_ids().next().expect("stereo element");
    let mut replacement = unknown
        .stereo_element(element)
        .expect("stereo element")
        .clone();
    match &mut replacement.kind {
        StereoElementKind::Tetrahedral(stereo) => stereo.orientation = None,
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
    unknown
        .replace_stereo_element(element, replacement)
        .expect("valid replacement");
    assert!(
        smiles_api::write(&unknown, smiles_api::SmilesWriteOptions::isomeric())
            .expect_err("unknown stereo should be rejected")
            .message
            .contains("unknown stereo")
    );
}

#[test]
fn isomeric_smiles_writes_directional_double_bond_elements() {
    for (input, expected_output, expected_orientation) in [
        ("C/C=C\\C", "C/C=C\\C", DoubleBondOrientation::Together),
        ("C/C=C/C", "C/C=C/C", DoubleBondOrientation::Opposite),
    ] {
        let mut molecule = read_smiles(input).expect("directional alkene should parse");
        perceive(&mut molecule).expect("directional alkene should perceive");

        let written = smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric())
            .expect("double-bond stereo should write");

        assert_eq!(written, expected_output);
        let mut reparsed = read_smiles(&written).expect("isomeric alkene output should parse");
        perceive(&mut reparsed).expect("isomeric alkene output should perceive");
        let stereo = reparsed
            .stereo_elements()
            .filter_map(|(_, element)| match &element.kind {
                StereoElementKind::DoubleBond(stereo) => Some(stereo),
                StereoElementKind::Tetrahedral(_) | StereoElementKind::Axis(_) => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(stereo.len(), 1);
        assert_eq!(stereo[0].orientation, Some(expected_orientation));
    }
}

#[test]
fn isomeric_smiles_preserves_pubchem_fused_quaternary_center() {
    let mut molecule = read_smiles("C[C@]12CCCC(C1CCC3=CC(=C(C=C23)C(=O)OC)C(=O)OC)(C)C")
        .expect("fused quaternary center should parse");
    perceive(&mut molecule).expect("fused quaternary center should perceive");
    let report =
        stereo_api::assign_cip_descriptors(&mut molecule).expect("CIP assignment should succeed");
    assert_eq!(report.assigned[0].descriptor, StereoDescriptor::S);

    let written = smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric())
        .expect("fused quaternary center should write");
    let mut reparsed = read_smiles(&written).expect("isomeric fused center output should parse");
    perceive(&mut reparsed).expect("isomeric fused center output should perceive");
    let report =
        stereo_api::assign_cip_descriptors(&mut reparsed).expect("CIP reassignment should succeed");

    assert_eq!(report.assigned[0].descriptor, StereoDescriptor::S);
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
        let mut molecule = crate::core::MoleculeEditor::new();
        let left = molecule
            .add_atom(carbon())
            .expect("atom identifier capacity");
        let right = molecule
            .add_atom(carbon())
            .expect("atom identifier capacity");
        let fluorine = molecule
            .add_atom(element_atom("F"))
            .expect("atom identifier capacity");
        let chlorine = molecule
            .add_atom(element_atom("Cl"))
            .expect("atom identifier capacity");
        molecule
            .add_bond(left, fluorine, BondOrder::Single)
            .expect("left carrier bond");
        let double_bond = molecule
            .add_bond(left, right, BondOrder::Double)
            .expect("double bond");
        molecule
            .add_bond(right, chlorine, BondOrder::Single)
            .expect("right carrier bond");
        molecule
            .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
                DoubleBondStereo {
                    bond: double_bond,
                    left,
                    right,
                    left_carrier,
                    right_carrier,
                    orientation: Some(stored_orientation),
                },
            )))
            .expect("double-bond stereo");
        molecule.working_mut().set_inferred_hydrogens(left, 1);
        molecule.working_mut().set_inferred_hydrogens(right, 1);

        let written = smiles_api::write(
            molecule.working(),
            smiles_api::SmilesWriteOptions::isomeric(),
        )
        .expect("implicit-carrier stereo should write");
        assert!(
            written.contains('/') || written.contains('\\'),
            "isomeric output should contain directional marks: {written}"
        );

        let mut reparsed = read_smiles(&written).expect("isomeric alkene output should parse");
        perceive(&mut reparsed).expect("isomeric alkene output should perceive");
        let stereo = reparsed
            .stereo_elements()
            .filter_map(|(_, element)| match &element.kind {
                StereoElementKind::DoubleBond(stereo) => Some(stereo),
                StereoElementKind::Tetrahedral(_) | StereoElementKind::Axis(_) => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(stereo.len(), 1);
        assert_eq!(stereo[0].orientation, Some(expected_orientation));
    }
}

#[test]
fn smiles_writer_reuses_ring_labels_after_closure() {
    let mut molecule = crate::core::MoleculeEditor::new();
    let atoms = (0..16)
        .map(|_| {
            molecule
                .add_atom(carbon())
                .expect("atom identifier capacity")
        })
        .collect::<Vec<_>>();
    for left in 0..atoms.len() {
        for right in (left + 1)..atoms.len() {
            molecule
                .add_bond(atoms[left], atoms[right], BondOrder::Single)
                .expect("complete graph bond should be valid");
        }
    }

    let written = smiles_api::write(
        molecule.working(),
        smiles_api::SmilesWriteOptions::default(),
    )
    .unwrap();
    let restored = read_smiles(&written).unwrap();
    assert_eq!(restored.atom_count(), molecule.atom_count());
    assert_eq!(restored.bond_count(), molecule.bond_count());
}

#[test]
fn smiles_directional_carrier_search_preserves_partial_branched_polyene_stereo() {
    let mut molecule = read_smiles("F/C=C(C=C(C)/C=C/F)/C").unwrap();
    perceive(&mut molecule).unwrap();
    assert_eq!(molecule.stereo_elements().count(), 2);
    let expected = assigned_descriptors(&mut molecule);
    for written in [
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()).unwrap(),
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()).unwrap(),
    ] {
        let mut restored = read_smiles(&written).unwrap();
        perceive(&mut restored).unwrap();
        assert_eq!(restored.stereo_elements().count(), 2, "{written}");
        let mut actual = assigned_descriptors(&mut restored);
        let mut expected = expected.clone();
        actual.sort_by_key(|value| format!("{value:?}"));
        expected.sort_by_key(|value| format!("{value:?}"));
        assert_eq!(actual, expected, "{written}");
    }
}

#[test]
fn smiles_aromatic_arsenic_round_trip() {
    let mut molecule = read_smiles("[as]1ccccc1").unwrap();
    perceive(&mut molecule).unwrap();
    for written in [
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::isomeric()).unwrap(),
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical()).unwrap(),
    ] {
        let mut restored = read_smiles(&written).unwrap();
        perceive(&mut restored).unwrap();
        assert_eq!(restored.atom_count(), 6);
        let hydrogens = |mol: &Molecule| {
            mol.atoms()
                .map(|(id, atom)| {
                    usize::from(atom.hydrogens.specified_count())
                        + usize::from(mol.inferred_hydrogens(id).unwrap().unwrap_or(0))
                })
                .sum::<usize>()
        };
        assert_eq!(hydrogens(&restored), hydrogens(&molecule));
    }
}

#[test]
fn smiles_canonicalization_retains_supplied_nonstereogenic_assertions() {
    let mut molecule = read_smiles("F/C=C(/F)F").unwrap();
    perceive(&mut molecule).unwrap();
    assert_eq!(molecule.stereo_elements().count(), 1);
    assert!(assigned_descriptors(&mut molecule).is_empty());
    for options in [
        smiles_api::SmilesWriteOptions::isomeric(),
        smiles_api::SmilesWriteOptions::canonical(),
    ] {
        let written = smiles_api::write(&molecule, options).unwrap();
        let mut restored = read_smiles(&written).unwrap();
        perceive(&mut restored).unwrap();
        assert_eq!(restored.stereo_elements().count(), 1);
        assert!(assigned_descriptors(&mut restored).is_empty());
    }
}

#[test]
fn smiles_percent_ring_tokens_consume_exactly_two_digits() {
    let source = "C%123CCCCC%12CCC3";
    let document = smiles_api::parse_str(source).unwrap();
    let rings = document
        .tokens()
        .iter()
        .filter(|token| token.kind() == smiles_api::SmilesDocumentTokenKind::Ring)
        .map(|token| &source[token.span()])
        .collect::<Vec<_>>();
    assert_eq!(rings, ["%12", "3", "%12", "3"]);
    let molecule = document.to_molecules().unwrap().pop().unwrap();
    assert_eq!(molecule.atom_count(), 9);
    assert_eq!(molecule.bond_count(), 10);
}
