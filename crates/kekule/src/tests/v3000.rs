use super::*;
use crate::properties::{PropertyKey, PropertyValue};

fn atom_cfg_tetrahedron(cfg: u8, fourth: Option<&str>) -> String {
    let atoms = if fourth.is_some() { 5 } else { 4 };
    let extra_atom = fourth
        .map(|symbol| format!("M  V30 15 {symbol} 0 -1 0 0\n"))
        .unwrap_or_default();
    let extra_bond = fourth.map(|_| "M  V30 4 1 99 15\n").unwrap_or_default();
    // Deliberately nonmonotonic atom serials and bond rows: parity follows
    // atom-block position, not serial number or incident-bond insertion order.
    format!("parity\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {atoms} {} 0 0 0\nM  V30 BEGIN ATOM\nM  V30 99 C 0 0 0 0 CFG={cfg}\nM  V30 20 F 1 0 0 0\nM  V30 70 Cl -1 0 0 0\nM  V30 4 Br 0 1 0 0\n{extra_atom}M  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 3 1 99 4\nM  V30 1 1 99 20\nM  V30 2 1 99 70\n{extra_bond}M  V30 END BOND\nM  V30 END CTAB\nM  END\n", atoms - 1)
}

#[test]
fn v3000_atom_cfg_preserves_parity_unknown_and_atom_block_order() {
    // Reference: RDKit 2026.03.6 AssignAtomChiralTagsFromMolParity,
    // followed by rdCIPLabeler.AssignCIPLabels (the default reader only
    // retains molParity, so the explicit parity conversion is required).
    for (fourth, expected) in [
        (None, StereoDescriptor::S),
        (Some("I"), StereoDescriptor::R),
        (Some("H"), StereoDescriptor::S),
    ] {
        for cfg in 0..=3 {
            let source = atom_cfg_tetrahedron(cfg, fourth);
            let (mut molecule, report) = read_molfile_with_report(&source).unwrap();
            assert!(!molecule.perception().has_valence());
            assert_eq!(
                report.created_stereo_elements().len(),
                usize::from(cfg != 0)
            );
            if cfg == 0 {
                continue;
            }
            assert_eq!(
                molecule
                    .stereo_elements()
                    .next()
                    .unwrap()
                    .1
                    .is_explicitly_unknown(),
                cfg == 3
            );
            perceive(&mut molecule).unwrap();
            let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
            if cfg == 3 {
                assert!(assigned.assigned.is_empty());
            } else {
                let expected = if cfg == 1 {
                    expected
                } else if expected == StereoDescriptor::R {
                    StereoDescriptor::S
                } else {
                    StereoDescriptor::R
                };
                assert_eq!(assigned.assigned[0].descriptor, expected);
            }
            for output in [
                molfile::write_v2000(&molecule).unwrap(),
                molfile::write_v3000(&molecule).unwrap(),
            ] {
                let mut reread = read_molfile(&output).unwrap();
                perceive(&mut reread).unwrap();
                assert_eq!(
                    stereo_api::assign_cip_descriptors(&mut reread)
                        .unwrap()
                        .assigned,
                    assigned.assigned
                );
            }
        }
    }
    for cfg in ["4", "-1", "1 CFG=1"] {
        assert!(molfile::parse_str(
            &atom_cfg_tetrahedron(1, None).replace("CFG=1", &format!("CFG={cfg}"))
        )
        .is_err());
    }
}

#[test]
fn v3000_atom_cfg_checks_redundant_wedges_and_unknown_precedence() {
    let source = atom_cfg_tetrahedron(1, None);
    let molecule = read_molfile(&source).unwrap();
    let wedged = molfile::write_v3000(&molecule).unwrap();
    let with_cfg = |cfg| {
        wedged.replace(
            "M  V30 1 C 0.0000 0.0000 0.0000 0",
            &format!("M  V30 1 C 0.0000 0.0000 0.0000 0 CFG={cfg}"),
        )
    };
    assert_eq!(
        read_molfile(&with_cfg(1))
            .unwrap()
            .stereo_elements()
            .next()
            .unwrap()
            .1,
        molecule.stereo_elements().next().unwrap().1
    );
    assert!(read_molfile(&with_cfg(2))
        .unwrap_err()
        .to_string()
        .contains("conflicts with bond wedge"));
    assert!(read_molfile(&with_cfg(3))
        .unwrap()
        .stereo_elements()
        .next()
        .unwrap()
        .1
        .is_explicitly_unknown());
    let wavy = with_cfg(1)
        .lines()
        .map(|line| {
            if line.starts_with("M  V30 ") && line.contains(" CFG=") && !line.contains(" C ") {
                line.replace("CFG=1", "CFG=2").replace("CFG=3", "CFG=2")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    assert!(read_molfile(&wavy)
        .unwrap()
        .stereo_elements()
        .next()
        .unwrap()
        .1
        .is_explicitly_unknown());
}

#[test]
fn v3000_atom_cfg_numbers_explicit_hydrogen_last_regardless_of_atom_row() {
    // CTfile Appendix A makes hydrogen highest numbered. RDKit's explicit
    // AssignAtomChiralTagsFromMolParity helper omits this special case in
    // 2026.03.6, so this assertion follows the format specification.
    let source = atom_cfg_tetrahedron(1, Some("H"));
    let hydrogen_row = "M  V30 15 H 0 -1 0 0\n";
    for before in [
        "M  V30 20 F",
        "M  V30 70 Cl",
        "M  V30 4 Br",
        "M  V30 END ATOM",
    ] {
        let reordered = source
            .replace(hydrogen_row, "")
            .replace(before, &format!("{hydrogen_row}{before}"));
        let mut molecule = read_molfile(&reordered).unwrap();
        perceive(&mut molecule).unwrap();
        assert_eq!(
            stereo_api::assign_cip_descriptors(&mut molecule)
                .unwrap()
                .assigned[0]
                .descriptor,
            StereoDescriptor::S
        );
    }
}

#[test]
fn molfile_model_preserves_drawn_e_z_and_rejects_inconsistent_or_degenerate_output() {
    for (smiles, last_y, expected) in [
        ("F/C=C/Cl", -1.0, StereoDescriptor::E),
        ("F/C=C\\Cl", 1.0, StereoDescriptor::Z),
    ] {
        let molecule = read_smiles(smiles).unwrap();
        let points = vec![
            Point3::new(-1.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, last_y, 0.0),
        ];
        let model = Model::from_molecule(&molecule, &test_positions(points.clone())).unwrap();
        let original = model.clone();
        for output in [
            molfile::write_model_v2000(&model).unwrap(),
            molfile::write_model_v3000(&model).unwrap(),
        ] {
            let mut reread = read_molfile(&output).unwrap();
            assert_eq!(reread.stereo_elements().count(), 1);
            assert!(!reread.perception().has_valence());
            perceive(&mut reread).unwrap();
            assert_eq!(
                stereo_api::assign_cip_descriptors(&mut reread)
                    .unwrap()
                    .assigned[0]
                    .descriptor,
                expected
            );
        }
        assert_eq!(model, original);
        for output in [
            molfile::write_v2000(&molecule),
            molfile::write_v3000(&molecule),
        ] {
            assert!(output.unwrap_err().message().contains("requires a Model"));
        }
        for invalid_y in [-last_y, 0.0, 0.00000001 * last_y] {
            let mut invalid_points = points.clone();
            invalid_points[3].y = invalid_y;
            let model = Model::from_molecule(&molecule, &test_positions(invalid_points)).unwrap();
            for output in [
                molfile::write_model_v2000(&model),
                molfile::write_model_v3000(&model),
            ] {
                assert!(output
                    .unwrap_err()
                    .message()
                    .contains("emitted coordinates"));
            }
        }
    }
}

#[test]
fn molfile_drawn_double_bond_unknown_annotations_override_coordinates() {
    let molecule = read_smiles("F/C=C/Cl").unwrap();
    let model = Model::from_molecule(
        &molecule,
        &test_positions(vec![
            Point3::new(-1.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, -1.0, 0.0),
        ]),
    )
    .unwrap();
    let source = molfile::write_model_v3000(&model).unwrap();
    for marked in [
        source.replace("M  V30 2 2 2 3", "M  V30 2 2 2 3 CFG=2"),
        source.replace("M  V30 1 1 1 2", "M  V30 1 1 1 2 CFG=2"),
    ] {
        for text in [
            marked.clone(),
            marked
                .replace("1.0000", "0.0000")
                .replace("2.0000", "0.0000"),
        ] {
            let parsed = read_molfile(&text).unwrap();
            assert_eq!(parsed.stereo_elements().count(), 1);
            assert!(parsed
                .stereo_elements()
                .next()
                .unwrap()
                .1
                .is_explicitly_unknown());
        }
    }
    assert!(read_molfile(
        &source
            .replace("1.0000", "0.0000")
            .replace("2.0000", "0.0000")
    )
    .unwrap()
    .stereo_elements()
    .next()
    .is_none());
}

#[test]
fn molfile_model_does_not_promote_unasserted_alkene_geometry_to_specified_stereo() {
    let molecule = read_smiles("FC=CCl").unwrap();
    assert_eq!(molecule.stereo_elements().count(), 0);
    let model = Model::from_molecule(
        &molecule,
        &test_positions(vec![
            Point3::new(-1.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, -1.0, 0.0),
        ]),
    )
    .unwrap();
    for output in [
        molfile::write_model_v2000(&model).unwrap(),
        molfile::write_model_v3000(&model).unwrap(),
    ] {
        let mut reread = read_molfile(&output).unwrap();
        assert_eq!(reread.stereo_elements().count(), 1);
        assert!(reread
            .stereo_elements()
            .next()
            .unwrap()
            .1
            .is_explicitly_unknown());
        perceive(&mut reread).unwrap();
        assert!(stereo_api::assign_cip_descriptors(&mut reread)
            .unwrap()
            .assigned
            .is_empty());
    }
    assert_eq!(molecule.stereo_elements().count(), 0);
}

#[test]
fn molfile_model_e_z_preserves_explicit_carrier_pairs_and_explicit_hydrogen() {
    for smiles in ["F/C(Cl)=C(Br)/I", "[H]/C(F)=C(Cl)/Br"] {
        let mut molecule = read_smiles(smiles).unwrap();
        perceive(&mut molecule).unwrap();
        let expected = stereo_api::assign_cip_descriptors(&mut molecule)
            .unwrap()
            .assigned;
        let model = Model::from_molecule(
            &molecule,
            &test_positions(vec![
                Point3::new(-1.0, 1.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(-1.0, -1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 1.0, 0.0),
                Point3::new(2.0, -1.0, 0.0),
            ]),
        )
        .unwrap();
        for output in [
            molfile::write_model_v2000(&model).unwrap(),
            molfile::write_model_v3000(&model).unwrap(),
        ] {
            let mut reread = read_molfile(&output).unwrap();
            perceive(&mut reread).unwrap();
            assert_eq!(
                stereo_api::assign_cip_descriptors(&mut reread)
                    .unwrap()
                    .assigned,
                expected
            );
        }
    }
}

#[test]
fn mol_v3000_parses_raw_atoms_bonds_coordinates_and_metadata() {
    let input = "\
charged radical
kekule benchmark
metadata fixture
  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.1000 0.2000 0.3000 7 CHG=1 RAD=2
M  V30 2 C 1.4000 0.0000 0.0000 0 MASS=13
M  V30 3 O 2.5000 0.0000 0.0000 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
";

    let small = read_molfile(input).expect("V3000 should parse");
    let mol = small;

    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
    assert!(mol
        .properties()
        .get(&PropertyKey::new("sdf.title").unwrap())
        .is_none());
    let atom0 = mol.atom(AtomId::new(0)).expect("atom exists");
    let atom1 = mol.atom(AtomId::new(1)).expect("atom exists");
    let atom2 = mol.atom(AtomId::new(2)).expect("atom exists");
    assert_eq!(atom0.element.symbol(), "N");
    assert_eq!(atom0.formal_charge, 1);
    assert_eq!(atom0.radical, Some(AtomRadical::Doublet));
    assert_eq!(atom0.atom_map, Some(7));
    assert_eq!(atom1.isotope, Some(13));
    assert_eq!(atom2.formal_charge, -1);
    let bond0 = mol.bond(BondId::new(0)).expect("bond exists");
    let bond1 = mol.bond(BondId::new(1)).expect("bond exists");
    assert_eq!(bond0.order, BondOrder::Single);
    assert_eq!(bond1.order, BondOrder::Double);
    assert_eq!(mol.atom_count(), 3);
}

#[test]
fn sdf_v3000_record_interpretation_retains_model_geometry() {
    let input = "\
v3000 geometry
kekule

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.2500 -2.5000 3.7500 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
";
    let document = sdf::parse_str(input).expect("V3000 SDF parses");
    let interpretation = document.records()[0]
        .interpret()
        .expect("V3000 SDF interprets");
    assert_eq!(interpretation.molecules().count(), 1);
    assert_eq!(interpretation.report().molfile_components().len(), 1);
    let point = interpretation.model().positions().values().value()[0];
    assert!((point.x - 0.125).abs() < 1.0e-15);
    assert!((point.y + 0.25).abs() < 1.0e-15);
    assert!((point.z - 0.375).abs() < 1.0e-15);
    assert!(interpretation
        .topology()
        .same_layout(interpretation.model().topology()));
}

#[test]
fn v3000_preserves_source_declared_tetrahedral_hydrogen_carrier() {
    for (cfg, expected_specified) in [(1, true), (3, true), (2, false)] {
        let input = format!(
            "\
stereo hydrogen
kekule

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0 HCOUNT=1
M  V30 2 F 1 0 0 0
M  V30 3 Cl -1 0 0 0
M  V30 4 Br 0 1 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG={cfg}
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
"
        );

        let (parsed, report) = read_molfile_with_report(&input).expect("V3000 should interpret");

        assert_eq!(
            parsed
                .atom(AtomId::new(0))
                .expect("stereo center")
                .hydrogens,
            HydrogenDeclaration::Fixed(1)
        );
        assert!(!parsed.perception().has_valence());
        assert_eq!(report.created_stereo_elements().len(), 1);
        assert_eq!(parsed.stereo_elements().count(), 1);
        assert_eq!(
            parsed
                .stereo_elements()
                .next()
                .expect("canonical stereo element")
                .1
                .is_specified(),
            expected_specified
        );

        let written = molfile::write_v3000(&parsed).expect("canonical stereo should project");
        let (reparsed, report) =
            read_molfile_with_report(&written).expect("projected V3000 stereo should re-interpret");
        assert_eq!(report.created_stereo_elements().len(), 1);
        assert_eq!(reparsed.stereo_elements().count(), 1);
        assert_eq!(
            reparsed
                .stereo_elements()
                .next()
                .expect("reparsed canonical stereo element")
                .1
                .is_specified(),
            expected_specified
        );
    }
}

#[test]
fn v3000_either_double_bond_publishes_unknown_canonical_stereo() {
    let input = "\
unknown double bond
kekule

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 F 0 0 0 0
M  V30 2 C 1 0 0 0
M  V30 3 C 2 0 0 0
M  V30 4 Cl 3 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3 CFG=2
M  V30 3 1 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
";

    let (molecule, report) =
        read_molfile_with_report(input).expect("V3000 either bond should interpret");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = molecule
        .stereo_element(report.created_stereo_elements()[0])
        .expect("canonical double-bond stereo element");
    assert!(element.is_explicitly_unknown());
    assert!(matches!(
        &element.kind,
        StereoElementKind::DoubleBond(stereo) if stereo.orientation.is_none()
    ));

    let written = molfile::write_v3000(&molecule).expect("unknown stereo should project");
    assert!(written.contains("CFG=2"));
    let reparsed = read_molfile(&written).expect("projected unknown stereo should interpret");
    assert!(reparsed.stereo_elements().any(|(_, element)| matches!(
        &element.kind,
        StereoElementKind::DoubleBond(stereo) if stereo.orientation.is_none()
    )));
}

#[test]
fn v3000_valence_is_source_semantics_but_unsupported_chemistry_is_interpretation_owned() {
    let valence = "\
declared valence
kekule

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0 VAL=4
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let document = molfile::parse_str(valence).expect("VAL is valid V3000 syntax");
    let molecule = molfile::interpret(&document)
        .expect("VAL can be interpreted from source semantics")
        .into_molecule();
    let carbon = molecule.atom(AtomId::new(0)).expect("carbon");
    assert_eq!(carbon.hydrogens, HydrogenDeclaration::Fixed(4));
    assert!(!molecule.perception().has_valence());

    let zero_declarations = valence.replace("VAL=4", "HCOUNT=-1 VAL=-1");
    let document =
        molfile::parse_str(&zero_declarations).expect("zero-count sentinels are valid syntax");
    let molecule = molfile::interpret(&document)
        .expect("zero-count sentinels have exact source semantics")
        .into_molecule();
    let carbon = molecule.atom(AtomId::new(0)).expect("carbon");
    assert_eq!(carbon.hydrogens, HydrogenDeclaration::Fixed(0));

    let undeclared = valence.replace(" VAL=4", "");
    let document = molfile::parse_str(&undeclared).expect("undeclared atom is valid syntax");
    let molecule = molfile::interpret(&document)
        .expect("undeclared hydrogen policy interprets")
        .into_molecule();
    assert_eq!(
        molecule.atom(AtomId::new(0)).expect("carbon").hydrogens,
        HydrogenDeclaration::Infer { explicit: 0 }
    );

    let unsupported = valence.replace("1 C 0 0 0 0 VAL=4", "1 Xx 0 0 0 0");
    let document = molfile::parse_str(&unsupported).expect("unknown symbol remains valid syntax");
    assert!(molfile::interpret(&document)
        .expect_err("core element support belongs to interpretation")
        .message()
        .contains("unsupported element"));
}

#[test]
fn mol_v3000_line_continuations_and_aromatic_bonds_localize_without_perception() {
    let input = "\
benzene-ish
kekule

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 -
M  V30 0
M  V30 2 C 1.4 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
";

    let small = read_molfile(input).expect("V3000 should parse");
    let mol = small;

    assert_eq!(
        mol.bond(BondId::new(0)).expect("bond").order,
        BondOrder::Double
    );
    assert_all_stale(&mol);
}

#[test]
fn malformed_mol_v3000_returns_errors_without_panicking() {
    let cases = [
        (
            "bad counts",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS nope 0 0 0 0\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "count mismatch",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "non-finite coordinates",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 1e999 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "bad endpoint",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "unsupported atom stereo",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 CFG=1\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "unsupported bond type",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 8 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "incomplete counts",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "zero atom index",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 0 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate bond index",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 3 C 2 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 1 1 2 3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate counts",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "counts after atom section",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate atom section",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 BEGIN ATOM\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate bond section",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "record outside CTAB",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 NOTE=outside\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "malformed atom option",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 BROKEN\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate atom option",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 CHG=1 CHG=2\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "unsupported bond option",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 TOPO=1\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate bond option",
            "Bad\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 CFG=1 CFG=1\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
    ];

    for (name, input) in cases {
        let parsed = std::panic::catch_unwind(|| read_molfile(input))
            .unwrap_or_else(|_| panic!("{name} panicked"));
        let error = parsed.expect_err("malformed V3000 input should fail");
        assert!(!error.to_string().is_empty(), "message for {name}");
    }
}

#[test]
fn mol_v3000_reports_only_nonstructural_unsupported_records_as_ignored() {
    let input = "\
collection
kekule

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/HILITE ATOMS=(1 1)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
";

    let document = molfile::parse_str(input).expect("unsupported collection is loss-preserved");
    assert_eq!(document.property_records().len(), 3);
    let interpretation =
        molfile::interpret(&document).expect("unsupported collection is reported, not hidden");
    assert_eq!(
        interpretation.report().ignored_record_lines(),
        &[12, 13, 14]
    );
}

#[test]
fn mol_v3000_parse_options_bound_input_counts_and_logical_lines() {
    let input = "\
bounded
kekule

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";

    let input_error = molfile::parse_str_with_options(
        input,
        molfile::MolfileParseOptions {
            max_input_bytes: input.len() - 1,
            ..molfile::MolfileParseOptions::default()
        },
    )
    .expect_err("Molfile input limit should apply");
    assert!(input_error.message().contains("input"));

    let atom_error = molfile::parse_str_with_options(
        input,
        molfile::MolfileParseOptions {
            max_v3000_atoms: 0,
            ..molfile::MolfileParseOptions::default()
        },
    )
    .expect_err("V3000 atom limit should apply");
    assert!(atom_error.message().contains("atom count"));

    let line_error = molfile::parse_str_with_options(
        input,
        molfile::MolfileParseOptions {
            max_v3000_logical_line_bytes: 4,
            ..molfile::MolfileParseOptions::default()
        },
    )
    .expect_err("V3000 logical line limit should apply");
    assert!(line_error.message().contains("logical line"));
}

#[test]
fn mol_v3000_writer_round_trips_supported_metadata() {
    let mut molecule = crate::core::MoleculeEditor::new();
    molecule
        .insert_property(
            PropertyKey::new("sdf.title").unwrap(),
            PropertyValue::String("metadata title".to_owned()),
        )
        .unwrap();
    molecule
        .insert_property(
            PropertyKey::new("sdf.program").unwrap(),
            PropertyValue::String("metadata program".to_owned()),
        )
        .unwrap();
    molecule
        .insert_property(
            PropertyKey::new("sdf.comment").unwrap(),
            PropertyValue::String("metadata comment".to_owned()),
        )
        .unwrap();

    let mut nitrogen = Atom::new(Element::from_symbol("N").expect("N"));
    nitrogen.formal_charge = 1;
    nitrogen.radical = Some(AtomRadical::Doublet);
    nitrogen.atom_map = Some(42);
    let n = molecule
        .add_atom(nitrogen)
        .expect("atom identifier capacity");

    let mut carbon = carbon();
    carbon.isotope = Some(13);
    let c = molecule.add_atom(carbon).expect("atom identifier capacity");

    let mut oxygen = oxygen();
    oxygen.formal_charge = -1;
    let o = molecule.add_atom(oxygen).expect("atom identifier capacity");

    molecule
        .add_bond(n, c, BondOrder::Single)
        .expect("single bond");
    molecule
        .add_bond(c, o, BondOrder::Double)
        .expect("double bond");

    let written = molfile::write_v3000(molecule.working()).expect("V3000 should write");
    assert_eq!(written.lines().nth(1), Some("kekule"));
    assert!(written.contains("V3000"));
    assert!(written.contains("CHG=1"));
    assert!(written.contains("MASS=13"));
    assert!(written.contains("RAD=2"));

    let reparsed = read_molfile(&written).expect("written V3000 should parse");
    assert!(reparsed
        .properties()
        .get(&PropertyKey::new("sdf.title").unwrap())
        .is_none());
    assert_eq!(
        reparsed.atom(AtomId::new(0)).expect("atom").formal_charge,
        1
    );
    assert_eq!(
        reparsed.atom(AtomId::new(0)).expect("atom").radical,
        Some(AtomRadical::Doublet)
    );
    assert_eq!(
        reparsed.atom(AtomId::new(0)).expect("atom").atom_map,
        Some(42)
    );
    assert_eq!(
        reparsed.atom(AtomId::new(1)).expect("atom").isotope,
        Some(13)
    );
    assert_eq!(reparsed.atom_count(), 3);
}

#[test]
fn mol_v3000_writer_rejects_unsupported_stereo_and_bonds() {
    let mut molecule = crate::core::MoleculeEditor::new();
    let a = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    molecule
        .working_mut()
        .graph
        .stereo_elements
        .push(Some(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center: a,
                carriers: vec![StereoCarrier::ImplicitHydrogen],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        ))));
    assert!(molfile::write_v3000(molecule.working())
        .expect_err("invalid stereo element should be rejected")
        .message
        .contains("cannot encode"));

    let mut molecule = crate::core::MoleculeEditor::new();
    let a = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let b = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let bond = molecule.add_bond(a, b, BondOrder::Double).expect("bond");
    let left_carrier = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let right_carrier = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    molecule
        .add_bond(a, left_carrier, BondOrder::Single)
        .expect("bond");
    molecule
        .add_bond(b, right_carrier, BondOrder::Single)
        .expect("bond");
    molecule
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond,
                left: a,
                right: b,
                left_carrier: StereoCarrier::Atom(left_carrier),
                right_carrier: StereoCarrier::Atom(right_carrier),
                orientation: Some(DoubleBondOrientation::Together),
            },
        )))
        .expect("double-bond stereo");
    assert!(molfile::write_v3000(molecule.working())
        .expect_err("specified double-bond stereo should be rejected")
        .message
        .contains("specified double-bond stereo"));

    let element = molecule
        .stereo_element_ids()
        .next()
        .expect("stereo element");
    molecule
        .remove_stereo_element(element)
        .expect("remove stereo element");
    molecule
        .bond_mut(bond)
        .expect("bond")
        .set_order(BondOrder::Quadruple);
    assert!(molfile::write_v3000(molecule.working())
        .expect_err("quadruple should be rejected")
        .message
        .contains("quadruple"));
}
