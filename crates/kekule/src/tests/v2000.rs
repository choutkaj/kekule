use super::*;
use crate::properties::{PropertyKey, PropertyValue};

fn wedged_tetrahedron(symbol: &str, charge_code: u8, mirror: bool) -> String {
    let y = if mirror { -1.0 } else { 1.0 };
    format!(
        "stereo drawing\nkekule\n\n  4  3  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 {symbol:<3} 0  {charge_code}  0  0  0  0\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0\n    0.0000{y:10.4}    0.0000 Br  0  0  0  0  0  0\n  1  2  1  1  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\nM  END\n"
    )
}

#[test]
fn molfile_wedges_use_drawing_geometry_for_hydrogen_and_lone_pair_carriers() {
    // RDKit 2026.03.5 AssignCIPLabels gives S for the original drawing,
    // and R after reflection, for all four fourth-carrier chemistries.
    for (symbol, charge_code) in [("C", 0), ("S", 0), ("S", 3), ("P", 0)] {
        for (mirror, expected) in [(false, StereoDescriptor::S), (true, StereoDescriptor::R)] {
            let source = wedged_tetrahedron(symbol, charge_code, mirror);
            let mut molecule = read_molfile(&source).expect("wedged center interprets");
            perceive(&mut molecule).expect("valence perceives");
            let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
            assert_eq!(assigned.assigned.len(), 1, "{symbol} {charge_code}");
            assert_eq!(assigned.assigned[0].descriptor, expected, "{source}");
        }
    }
}

#[test]
fn molfile_model_writing_preserves_tetrahedral_drawing_orientation() {
    for (symbol, charge_code) in [("C", 0), ("S", 3), ("P", 0)] {
        for mirror in [false, true] {
            let document = molfile::parse_str(&wedged_tetrahedron(symbol, charge_code, mirror))
                .expect("valid source drawing");
            let interpreted = molfile::interpret(&document).expect("interpreted drawing");
            let original = interpreted.molecules().next().unwrap();
            for written in [
                molfile::write_model_v2000(interpreted.model()).expect("V2000 model writes"),
                molfile::write_model_v3000(interpreted.model()).expect("V3000 model writes"),
            ] {
                let reparsed = read_molfile(&written).expect("written drawing interprets");
                assert_eq!(
                    original
                        .stereo_elements()
                        .map(|(_, element)| &element.kind)
                        .collect::<Vec<_>>(),
                    reparsed
                        .stereo_elements()
                        .map(|(_, element)| &element.kind)
                        .collect::<Vec<_>>(),
                    "{written}"
                );
            }
        }
    }
}

#[test]
fn molfile_redundant_wedges_preserve_consistent_and_unknown_configurations() {
    let source = wedged_tetrahedron("C", 0, false)
        .replace("   -1.0000    0.0000", "   -0.5000    0.8660")
        .replace("    0.0000    1.0000", "   -0.5000   -0.8660");
    let original = read_molfile(&source).expect("single wedge interprets");
    let redundant = source.replace("  1  3  1  0", "  1  3  1  1");
    let (molecule, report) =
        read_molfile_with_report(&redundant).expect("redundant wedges interpret");
    assert!(report.warnings().is_empty());
    assert_eq!(
        molecule.stereo_elements().next().unwrap().1.kind,
        original.stereo_elements().next().unwrap().1.kind
    );

    let unknown = source.replace("  1  3  1  0", "  1  3  1  4");
    let molecule = read_molfile(&unknown).expect("wavy mark interprets");
    assert!(molecule
        .stereo_elements()
        .next()
        .unwrap()
        .1
        .is_explicitly_unknown());

    let degenerate = unknown
        .replace("    1.0000    0.0000", "    0.0000    0.0000")
        .replace("   -0.5000    0.8660", "    0.0000    0.0000")
        .replace("   -0.5000   -0.8660", "    0.0000    0.0000");
    let (molecule, report) = read_molfile_with_report(&degenerate).unwrap();
    assert!(report.warnings().is_empty());
    assert!(molecule
        .stereo_elements()
        .next()
        .unwrap()
        .1
        .is_explicitly_unknown());

    let conflicting = source.replace("  1  3  1  0", "  1  3  1  6");
    let (molecule, report) = read_molfile_with_report(&conflicting).expect("conflict is reported");
    assert_eq!(report.warnings().len(), 1);
    assert!(molecule.stereo_elements().next().is_none());
}

#[test]
fn v3000_round_trips_absolute_or_and_stereo_groups_and_promotes_auto_output() {
    let document = molfile::parse_str(&wedged_tetrahedron("C", 0, false)).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let source = molfile::write_model_v3000(interpreted.model()).unwrap();
    for (group, expected) in [
        ("MDLV30/STEABS", StereoGroupKind::Absolute),
        ("MDLV30/STEREL1", StereoGroupKind::Or),
        ("MDLV30/STERAC1", StereoGroupKind::And),
    ] {
        let grouped = source.replace(
            "M  V30 END CTAB",
            &format!("M  V30 BEGIN COLLECTION\nM  V30 {group} ATOMS=(1 1)\nM  V30 END COLLECTION\nM  V30 END CTAB"),
        );
        let document = molfile::parse_str(&grouped).expect("collection syntax is preserved");
        let interpreted = molfile::interpret(&document).expect("group semantics are represented");
        let group = interpreted
            .molecules()
            .next()
            .unwrap()
            .stereo_groups()
            .next()
            .unwrap()
            .1;
        assert_eq!(group.kind, expected);
        assert_eq!(group.members.len(), 1);
        assert!(interpreted.reports()[0].ignored_record_lines().is_empty());
        let written =
            molfile::write_model(interpreted.model(), molfile::MolfileWriteOptions::default())
                .unwrap();
        assert!(written.contains("V3000"));
        let document = molfile::parse_str(&written).unwrap();
        let reparsed = molfile::interpret(&document).unwrap();
        assert_eq!(
            reparsed
                .molecules()
                .next()
                .unwrap()
                .stereo_groups()
                .next()
                .unwrap()
                .1,
            group
        );
        assert!(molfile::write_model_v2000(interpreted.model())
            .unwrap_err()
            .message()
            .contains("enhanced stereo groups"));
    }
}

#[test]
fn v3000_round_trips_atropisomeric_bond_group_members() {
    let document = molfile::parse_str(rdkit_rp6306_atrop_molblock()).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let source = molfile::write_model_v3000(interpreted.model()).unwrap();
    // RDKit represents enhanced axis membership using either endpoint atom,
    // and collapses a pair of endpoint references to one bond member.
    for (name, kind) in [
        ("STEABS", StereoGroupKind::Absolute),
        ("STERAC1", StereoGroupKind::And),
        ("STEREL1", StereoGroupKind::Or),
    ] {
        for atoms in ["1 3", "1 9", "2 3 9"] {
            let grouped = source.replace("M  V30 END CTAB", &format!("M  V30 BEGIN COLLECTION\nM  V30 MDLV30/{name} ATOMS=({atoms})\nM  V30 END COLLECTION\nM  V30 END CTAB"));
            let interpretation = molfile::parse_str(&grouped).unwrap().interpret().unwrap();
            let molecule = interpretation.molecules().next().unwrap();
            let group = molecule.stereo_groups().next().unwrap().1;
            assert_eq!(group.kind, kind);
            assert_eq!(group.members.len(), 1);
            assert!(
                matches!(&molecule.stereo_element(group.members[0]).unwrap().kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(3))
            );
            let output = molfile::write_model_v3000(interpretation.model()).unwrap();
            assert!(!output.contains("BONDS="));
            assert!(output.contains("ATOMS=(1 3)"));
            let reread = molfile::parse_str(&output).unwrap().interpret().unwrap();
            assert_eq!(
                reread
                    .molecules()
                    .next()
                    .unwrap()
                    .stereo_groups()
                    .next()
                    .unwrap()
                    .1,
                group
            );
        }
    }
}
#[test]
fn molfile_atropisomeric_wedges_validate_all_marks_and_preserve_unknown_stereo() {
    let source = rdkit_rp6306_atrop_molblock().replace("  9 12  1  6", "  9 12  1  0");
    let marked = |left, right| {
        source
            .replace("  3  7  1  0", &format!("  3  7  1  {left}"))
            .replace("  3 10  1  0", &format!("  3 10  1  {right}"))
    };
    // RDKit 2026.03.5 accepts opposite directions at one end and rejects
    // two wedges or two hashes at that end.
    for (left, right, expected) in [(1, 6, StereoDescriptor::M), (6, 1, StereoDescriptor::P)] {
        let mut molecule = read_molfile(&marked(left, right)).unwrap();
        perceive(&mut molecule).unwrap();
        let assigned = stereo_api::assign_cip_descriptors(&mut molecule).unwrap();
        assert_eq!(assigned.assigned.len(), 1);
        assert_eq!(assigned.assigned[0].descriptor, expected);
    }
    for direction in [1, 6] {
        let document = molfile::parse_str(&marked(direction, direction)).unwrap();
        let error = molfile::interpret(&document).unwrap_err();
        assert!(error
            .message()
            .contains("ConflictingAtropisomericWedgeMarks"));
    }
    for (left, right) in [(4, 0), (4, 1), (4, 6), (1, 4), (6, 4), (4, 4)] {
        let mut molecule = read_molfile(&marked(left, right)).unwrap();
        let elements = molecule
            .stereo_elements()
            .map(|(_, element)| element)
            .collect::<Vec<_>>();
        assert_eq!(elements.len(), 1);
        assert!(matches!(&elements[0].kind, StereoElementKind::Axis(stereo)
            if stereo.axis == BondId::new(3) && stereo.orientation.is_none()));
        perceive(&mut molecule).unwrap();
        assert!(stereo_api::assign_cip_descriptors(&mut molecule)
            .unwrap()
            .assigned
            .is_empty());
        assert!(molfile::write_v3000(&molecule)
            .unwrap_err()
            .message()
            .contains("unknown axis"));
    }
}

#[test]
fn v3000_stereo_groups_validate_members_and_preserve_source_ids_and_continuations() {
    let document = molfile::parse_str(&wedged_tetrahedron("C", 0, false)).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let source = molfile::write_model_v3000(interpreted.model())
        .unwrap()
        .replace("M  V30 1 C", "M  V30 101 C")
        .replace("M  V30 1 1 1 2", "M  V30 1 1 101 2")
        .replace("M  V30 2 1 1 3", "M  V30 2 1 101 3")
        .replace("M  V30 3 1 1 4", "M  V30 3 1 101 4");
    let collection = |row: &str| {
        source.replace(
            "M  V30 END CTAB",
            &format!(
                "M  V30 BEGIN COLLECTION\nM  V30 {row}\nM  V30 END COLLECTION\nM  V30 END CTAB"
            ),
        )
    };
    let continued = collection("MDLV30/STEREL1 ATOMS=(1 -\nM  V30 101)");
    let document = molfile::parse_str(&continued).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let group = interpreted
        .molecules()
        .next()
        .unwrap()
        .stereo_groups()
        .next()
        .unwrap()
        .1;
    assert_eq!(group.kind, StereoGroupKind::Or);
    assert_eq!(group.members.len(), 1);
    assert!(interpreted.reports()[0].ignored_record_lines().is_empty());

    for row in [
        "MDLV30/STEREL0 ATOMS=(1 101)",
        "MDLV30/STEREL1 ATOMS=(2 101)",
        "MDLV30/STEREL1 ATOMS=(0)",
        "MDLV30/STEREL1 ATOMS=(1 999)",
        "MDLV30/STEREL1 ATOMS=(2 101 101)",
        "MDLV30/STEREL1 ATOMS=(1 2)",
        "MDLV30/STEREL1 ATOMS=(1 101) ATOMS=(1 101)",
        "MDLV30/STEREL1 BONDS=(1 1)",
    ] {
        let grouped = collection(row);
        if let Ok(document) = molfile::parse_str(&grouped) {
            assert!(molfile::interpret(&document).is_err(), "{row}");
        }
    }
}

#[test]
fn v3000_repeated_group_ids_preserve_one_relation() {
    let molecule = read_smiles("F[C@H](Cl)[C@H](Br)I").unwrap();
    let source = molfile::write_v3000(&molecule).unwrap();
    for kind in ["STEREL", "STERAC"] {
        for number in ["1", "01"] {
            let grouped = source.replace("M  V30 END CTAB", &format!("M  V30 BEGIN COLLECTION\nM  V30 MDLV30/{kind}1 ATOMS=(1 2)\nM  V30 MDLV30/{kind}{number} ATOMS=(1 4)\nM  V30 END COLLECTION\nM  V30 END CTAB"));
            let document = molfile::parse_str(&grouped).unwrap();
            let interpreted = molfile::interpret(&document).unwrap();
            let molecule = interpreted.molecules().next().unwrap();
            assert_eq!(molecule.stereo_groups().count(), 1);
            assert_eq!(molecule.stereo_groups().next().unwrap().1.members.len(), 2);
            assert!(interpreted.reports()[0].ignored_record_lines().is_empty());
        }
    }
}

#[test]
fn v3000_rejects_relative_groups_across_components_and_preserves_absolute_members() {
    let document = molfile::parse_str(&wedged_tetrahedron("C", 0, false)).unwrap();
    let interpreted = molfile::interpret(&document).unwrap();
    let source = molfile::write_model_v3000(interpreted.model()).unwrap()
        .replace("COUNTS 4 3", "COUNTS 8 6")
        .replace("M  V30 END ATOM", "M  V30 5 C 5 0 0 0\nM  V30 6 F 6 0 0 0\nM  V30 7 Cl 4 0 0 0\nM  V30 8 Br 5 1 0 0\nM  V30 END ATOM")
        .replace("M  V30 END BOND", "M  V30 4 1 5 6 CFG=1\nM  V30 5 1 5 7\nM  V30 6 1 5 8\nM  V30 END BOND");
    for group in ["MDLV30/STERAC1", "MDLV30/STEREL1", "MDLV30/STEABS"] {
        let source = source.replace("M  V30 END CTAB", &format!("M  V30 BEGIN COLLECTION\nM  V30 {group} ATOMS=(2 1 5)\nM  V30 END COLLECTION\nM  V30 END CTAB"));
        let document = molfile::parse_str(&source).unwrap();
        if group == "MDLV30/STEABS" {
            let interpreted = molfile::interpret(&document).unwrap();
            assert_eq!(interpreted.molecules().count(), 2);
            assert!(interpreted
                .molecules()
                .all(|molecule| molecule.stereo_groups().count() == 1));
            let output = molfile::write_model_v3000(interpreted.model()).unwrap();
            assert_eq!(output.matches("MDLV30/STEABS").count(), 1);
            assert!(output.contains("ATOMS=(2 1 5)"));
            assert_eq!(
                molfile::parse_str(&output)
                    .unwrap()
                    .to_molecules()
                    .unwrap()
                    .len(),
                2
            );
        } else {
            assert!(molfile::interpret(&document)
                .unwrap_err()
                .message()
                .contains("spanning disconnected molecules"));
        }
    }
}

#[test]
fn molfile_and_sdf_documents_preserve_record_metadata_before_interpretation() {
    let molfile_text = "Header title\nprogram line\ncomment line\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nX  UNSUPPORTED\nM  END\n";
    let document = molfile::parse_str(molfile_text).expect("Molfile document parses");
    assert_eq!(document.header().title(), "Header title");
    assert_eq!(document.unsupported_records().len(), 1);
    let interpretation = molfile::interpret(&document).expect("Molfile interprets");
    let molecule = interpretation.molecule();
    assert!(molecule
        .properties()
        .get(&PropertyKey::new("sdf.title").unwrap())
        .is_none());
    assert_eq!(interpretation.report().atom_mappings().len(), 1);
    assert_eq!(interpretation.report().ignored_record_lines(), &[6]);

    let sdf_text = format!("{molfile_text}>  <FIELD>\nvalue\n\n$$$$\n");
    let document = sdf::parse_str(&sdf_text).expect("SDF parses");
    assert_eq!(document.records()[0].data_fields()[0].value(), "value");
    let interpretation = sdf::interpret(&document).expect("SDF interprets");
    let records = interpretation.records();
    assert_eq!(records[0].title(), "Header title");
    assert_eq!(records[0].data_fields()[0].name(), "FIELD");
    assert!(records[0]
        .molecule()
        .properties()
        .get(&PropertyKey::new("sdf.field.FIELD").unwrap())
        .is_none());
    assert_eq!(interpretation.reports().len(), 1);
}

#[test]
fn molfile_document_reports_nonempty_content_after_m_end_as_unsupported() {
    let input = "Header title\nprogram line\ncomment line\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\ntrailing content\n";
    let document = molfile::parse_str(input).expect("Molfile document parses");
    assert_eq!(document.unsupported_records().len(), 1);
    assert_eq!(document.unsupported_records()[0].number(), 7);
    assert_eq!(document.unsupported_records()[0].text(), "trailing content");
    let interpretation = molfile::interpret(&document).expect("Molfile interprets");
    assert_eq!(interpretation.report().ignored_record_lines(), &[7]);
}

#[test]
fn molfile_and_sdf_documents_parse_adjacent_three_digit_counts() {
    let mut molfile_text =
        String::from("Large\nprogram\ncomment\n999999  0  0  0  0            999 V2000\n");
    for _ in 0..999 {
        molfile_text.push_str("    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n");
    }
    for atom in 1..999 {
        molfile_text.push_str(&format!("{atom:>3}{:>3}  1  0  0  0  0\n", atom + 1));
    }
    molfile_text.push_str("999  1  1  0  0  0  0\n");
    molfile_text.push_str("M  END\n");

    let document = molfile::parse_str(&molfile_text).expect("fixed-width counts parse");
    assert_eq!(document.atom_records().len(), 999);
    assert_eq!(document.bond_records().len(), 999);

    let sdf_text = format!("{molfile_text}$$$$\n");
    let document =
        sdf::parse_str(&sdf_text).expect("SDF delegates to fixed-width Molfile counts parsing");
    assert_eq!(document.records()[0].molfile().atom_records().len(), 999);
    assert_eq!(document.records()[0].molfile().bond_records().len(), 999);
}

#[test]
fn molfile_document_parser_validates_declared_atom_and_bond_records() {
    let invalid_atom =
        "Bad\nprogram\ncomment\n  1  0  0  0  0  0            999 V2000\natom record\nM  END\n";
    let error =
        molfile::parse_str(invalid_atom).expect_err("invalid atom syntax must fail parsing");
    assert_eq!(error.line, 5);
    assert!(error.message.contains("atom"));

    let invalid_bond = "Bad\nprogram\ncomment\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0\nbond record\nM  END\n";
    let error =
        molfile::parse_str(invalid_bond).expect_err("invalid bond syntax must fail parsing");
    assert_eq!(error.line, 7);
    assert!(error.message.contains("bond"));
}

#[test]
fn sdf_v2000_parses_single_record_atoms_bonds_and_fields() {
    let input = "\
Water
  kekule
comment
  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0
    1.0000    0.0000    0.0000 H   0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
>  <NAME>
water

$$$$
";

    let records = read_sdf_records(input).expect("record should parse");
    let mol = records[0].molecule();

    assert_eq!(records.len(), 1);
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(
        mol.atom(AtomId::new(0))
            .expect("atom exists")
            .element
            .symbol(),
        "O"
    );
    assert_eq!(
        mol.bond(BondId::new(0)).expect("bond exists").order,
        BondOrder::Single
    );
    assert_eq!(records[0].data_fields()[0].value(), "water");
}

#[test]
fn sdf_v2000_parses_multiple_records_in_order() {
    let input = "\
One
  kekule

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
$$$$
Two
  kekule

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0
M  END
$$$$
";

    let records = read_sdf_records(input).expect("records should parse");

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].title(), "One");
    assert_eq!(records[1].title(), "Two");
    assert_eq!(
        records[1]
            .molecule()
            .atom(AtomId::new(0))
            .expect("atom exists")
            .element
            .symbol(),
        "O"
    );
}

#[test]
fn sdf_v2000_can_allow_missing_final_delimiter() {
    let input = "\
Methane
  kekule

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
";

    let molecules = read_sdf_molecules_with_options(
        input,
        SdfParseOptions {
            allow_missing_final_delimiter: true,
            ..SdfParseOptions::default()
        },
    )
    .expect("record should parse");

    assert_eq!(molecules.len(), 1);
    assert_eq!(molecules[0].atom_count(), 1);
}

#[test]
fn sdf_v2000_requires_the_final_record_delimiter_by_default() {
    let complete = "\
One
  kekule

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
$$$$
";
    let unterminated = "\
Two
  kekule

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0
M  END
";
    let input = format!("{complete}{unterminated}");

    let error = sdf::parse_str(&input)
        .expect_err("a previous delimiter must not waive the final delimiter");
    assert_eq!(error.record(), 2);
    assert!(error.message().contains("missing final"));

    let document = sdf::parse_str_with_options(
        &input,
        SdfParseOptions {
            allow_missing_final_delimiter: true,
            ..SdfParseOptions::default()
        },
    )
    .expect("the explicit permissive option accepts the final record");
    assert_eq!(document.records().len(), 2);
}

#[test]
fn sdf_v2000_rejects_unstructured_post_ctab_text_and_truly_unterminated_fields() {
    let molfile = "\
One
  kekule

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
";
    let stray = format!("{molfile}orphan text\n$$$$\n");
    let error =
        sdf::parse_str(&stray).expect_err("unstructured post-CTAB content must not be discarded");
    assert!(error.message().contains("unexpected content"));

    let delimited_field = format!("{molfile}>  <FIELD>\nvalue\n$$$$\n");
    let document = sdf::parse_str(&delimited_field)
        .expect("the record delimiter unambiguously terminates the final field");
    assert_eq!(document.records()[0].data_fields()[0].value(), "value");

    let unterminated_field = format!("{molfile}>  <FIELD>\nvalue\n");
    let error = sdf::parse_str_with_options(
        &unterminated_field,
        SdfParseOptions {
            allow_missing_final_delimiter: true,
            ..SdfParseOptions::default()
        },
    )
    .expect_err("a field at bare end-of-input still requires a blank terminator");
    assert!(error.message().contains("terminating blank line"));
}

#[test]
fn sdf_v2000_parse_limits_bound_input_records_and_record_size() {
    let record = "\
One
  kekule

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
$$$$
";
    let two_records = format!("{record}{record}");

    let input_error = sdf::parse_str_with_options(
        record,
        SdfParseOptions {
            max_input_bytes: record.len() - 1,
            ..SdfParseOptions::default()
        },
    )
    .expect_err("input byte limit should apply before parsing");
    assert!(input_error.message().contains("input"));

    let record_count_error = sdf::parse_str_with_options(
        &two_records,
        SdfParseOptions {
            max_records: 1,
            ..SdfParseOptions::default()
        },
    )
    .expect_err("record count limit should reject the second record");
    assert_eq!(record_count_error.record(), 2);
    assert!(record_count_error.message().contains("record count"));

    let record_size_error = sdf::parse_str_with_options(
        record,
        SdfParseOptions {
            max_record_bytes: 1,
            ..SdfParseOptions::default()
        },
    )
    .expect_err("record byte limit should apply while scanning");
    assert_eq!(record_size_error.record(), 1);
    assert!(record_size_error.message().contains("record exceeds"));
}

#[test]
fn sdf_v2000_rejects_v3000_and_bad_endpoints() {
    let v3000 = "\
V3000
  kekule

  0  0  0  0  0  0            999 V3000
M  END
$$$$
";
    let err = read_sdf_molecules(v3000).expect_err("V3000 should fail");
    assert!(!err.to_string().is_empty());

    let bad_endpoint = "\
Bad
  kekule

  1  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$$$$
";
    let err = read_sdf_molecules(bad_endpoint).expect_err("bad endpoint should fail");
    assert!(err.to_string().contains("outside atom block"));
}

#[test]
fn v2000_malformed_structural_fields_return_errors_without_panicking() {
    let cases = [
            (
                "zero endpoint",
                "Bad\nkekule\n\n  1  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n  0  1  1  0  0  0  0\nM  END\n",
            ),
            (
                "non-ASCII counts",
                "Bad\nkekule\n\né  1  0  0  0  0            999 V2000\nM  END\n",
            ),
            (
                "non-ASCII atom",
                "Bad\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 Cé  0  0  0  0  0  0\nM  END\n",
            ),
            (
                "truncated atom",
                "Bad\nkekule\n\n  1  0  0  0  0  0            999 V2000\n0.0 C\nM  END\n",
            ),
            (
                "non-ASCII bond",
                "Bad\nkekule\n\n  1  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n  1  é  1  0\nM  END\n",
            ),
            (
                "count over format limit",
                "Bad\nkekule\n\n1000 0 V2000\nM  END\n",
            ),
            (
                "inconsistent counts",
                "Bad\nkekule\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n",
            ),
            (
                "truncated M record",
                "Bad\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  CHG  2   1   1\nM  END\n",
            ),
            (
                "zero M-record atom",
                "Bad\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  CHG  1   0   1\nM  END\n",
            ),
        ];

    for (name, input) in cases {
        let parsed = std::panic::catch_unwind(|| read_molfile(input))
            .unwrap_or_else(|_| panic!("{name} panicked"));
        let error = parsed.expect_err("malformed V2000 input should fail");
        assert!(!error.to_string().is_empty(), "message for {name}");
    }
}

#[test]
fn sdf_v2000_aromatic_source_is_localized_without_perception() {
    let input = "\
Benzene-ish
  kekule

  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0
  1  2  4  0  0  0  0
M  END
$$$$
";

    let molecules = read_sdf_molecules(input).expect("record should parse");
    let mol = &molecules[0];

    assert_all_stale(mol);
    assert_eq!(
        mol.bond(BondId::new(0)).expect("bond exists").order,
        BondOrder::Double
    );
}

#[test]
fn mol_v2000_preserves_coordinates_charges_isotopes_radicals_and_atom_maps() {
    let input = "\
charged radical
kekule benchmark
metadata fixture
  2  1  0  0  0  0            999 V2000
    0.1000    0.2000    0.3000 N   0  0  0  0  0  0  0  0  0  7  0  0
    1.4000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  CHG  1   1   1
M  ISO  1   2  13
M  RAD  1   1   2
M  END
";

    let small = read_molfile(input).expect("mol should parse");
    let atom0 = small.atom(AtomId::new(0)).expect("atom exists");
    let atom1 = small.atom(AtomId::new(1)).expect("atom exists");
    assert_eq!(atom0.formal_charge, 1);
    assert_eq!(atom0.radical, Some(AtomRadical::Doublet));
    assert_eq!(atom0.atom_map, Some(7));
    assert_eq!(atom1.isotope, Some(13));
    assert_eq!(small.atom_count(), 2);
}

#[test]
fn v2000_atom_block_charge_code_four_preserves_a_doublet_radical() {
    let input = "\
doublet
kekule

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  4  0  0  0  0  0  0  0  0  0  0
M  END
";

    let molecule = read_molfile(input).expect("atom-block doublet radical should parse");
    let atom = molecule.atom(AtomId::new(0)).expect("radical atom");
    assert_eq!(atom.formal_charge, 0);
    assert_eq!(atom.radical, Some(AtomRadical::Doublet));
}

#[test]
fn sdf_v2000_fields_round_trip_leading_greater_than_lines_and_reject_unsafe_metadata() {
    let molecule = read_smiles("C").expect("methane parses");
    let record = SdfRecordInterpretation::new(
        "safe title",
        test_model(&molecule),
        vec![SdfDataField::new("NOTES", "> leading marker\nsecond line")],
    );
    let written = sdf::write_v2000(&[record]).expect("representable field should write");
    let reparsed = read_sdf_records(&written).expect("written field should parse");
    assert_eq!(
        reparsed[0].data_fields()[0].value(),
        "> leading marker\nsecond line"
    );

    for (title, field, expected) in [
        (
            "unsafe\ntitle",
            SdfDataField::new("FIELD", "value"),
            "titles",
        ),
        ("safe", SdfDataField::new(" BAD ", "value"), "field names"),
        (
            "safe",
            SdfDataField::new("FIELD", "first\n\nthird"),
            "blank lines",
        ),
        (
            "safe",
            SdfDataField::new("FIELD", "first\n$$$$\nthird"),
            "record delimiter",
        ),
    ] {
        let record = SdfRecordInterpretation::new(title, test_model(&molecule), vec![field]);
        let error =
            sdf::write_v2000(&[record]).expect_err("unrepresentable SDF metadata must fail");
        assert!(error.message().contains(expected), "{expected}: {error}");
    }
}

#[test]
fn v2000_radical_codes_round_trip_exact_multiplicity() {
    for (code, expected) in [
        (1, AtomRadical::Singlet),
        (2, AtomRadical::Doublet),
        (3, AtomRadical::Triplet),
    ] {
        let input = format!(
                "radical {code}\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  RAD  1   1   {code}\nM  END\n"
            );
        let parsed = read_molfile(&input).expect("radical record should parse");
        assert_eq!(
            parsed.atom(AtomId::new(0)).expect("atom").radical,
            Some(expected)
        );

        let written = molfile::write_v2000(&parsed).expect("radical record should write");
        assert!(
            written.contains(&format!("M  RAD  1   1   {code}")),
            "written code {code}: {written}"
        );
        let reparsed = read_molfile(&written).expect("written radical record should parse");
        assert_eq!(
            reparsed.atom(AtomId::new(0)).expect("atom").radical,
            Some(expected)
        );
    }
}

#[test]
fn v2000_bond_stereo_requires_enough_source_context_to_canonicalize() {
    for (order_code, stereo_code) in [(1, 1), (1, 4), (1, 6), (2, 3)] {
        let input = format!(
                "stereo\nkekule\n\n  2  1  0  0  0  0            999 V2000\n   -1.2500    0.0000    0.0000 C   0  0  0  0  0  0\n    1.2500    0.0000    0.0000 C   0  0  0  0  0  0\n  1  2  {order_code}  {stereo_code}  0  0  0\nM  END\n"
            );
        let document = molfile::parse_str(&input).expect("bond stereo syntax should parse");
        let error = molfile::interpret(&document)
            .expect_err("under-specified stereo must not publish a molecule");
        assert_eq!(error.line(), 7);
        assert!(error.message().contains("source-stereo canonicalization"));
    }
}

#[test]
fn v2000_materializes_omitted_tetrahedral_hydrogen_from_source_valence() {
    for symbol in ["C", "N", "S", "P"] {
        let input = format!(
            "stereo hydrogen\nkekule\n\n  4  3  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 {symbol:<3} 0  0  0  0  0  0\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0\n    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0\n  1  2  1  1  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\nM  END\n"
        );

        if symbol != "N" {
            let molecule = read_molfile(&input).expect("source fourth carrier interprets");
            let expected_carrier = if symbol == "P" {
                StereoCarrier::ImplicitLonePair
            } else {
                assert_eq!(
                    molecule.atom(AtomId::new(0)).unwrap().hydrogens,
                    HydrogenDeclaration::Fixed(1)
                );
                StereoCarrier::ImplicitHydrogen
            };
            assert!(!molecule.perception().has_valence());
            assert_eq!(molecule.stereo_elements().count(), 1);
            assert!(molecule.stereo_elements().any(|(_, element)| {
                matches!(
                    &element.kind,
                    StereoElementKind::Tetrahedral(stereo)
                        if stereo.carriers.contains(&expected_carrier)
                )
            }));
        } else {
            let error = read_molfile(&input)
                .expect_err("wedge without a declared fourth carrier must fail");
            assert!(error.to_string().contains("UnassembledTetrahedralBondMark"));
        }
    }
}

#[test]
fn v2000_source_hydrogen_and_valence_declarations_define_stereo_carriers() {
    for (declaration_fields, expected_hydrogens) in
        [("0  0  0  2  0  0", 1), ("0  0  0  0  0  4", 1)]
    {
        for (stereo_code, expected_specified) in [(1, true), (6, true), (4, false)] {
            let input = format!(
                "declared stereo hydrogen\nkekule\n\n  4  3  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   {declaration_fields}\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0\n    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0\n  1  2  1  {stereo_code}  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\nM  END\n"
            );

            let document = molfile::parse_str(&input).expect("source syntax parses");
            let interpreted = molfile::interpret(&document).expect("source declaration interprets");
            assert_eq!(interpreted.report().created_stereo_elements().len(), 1);
            let molecule = interpreted.into_molecule();
            let center = molecule.atom(AtomId::new(0)).expect("stereo center");
            assert_eq!(
                center.hydrogens,
                HydrogenDeclaration::Fixed(expected_hydrogens)
            );
            assert!(!molecule.perception().has_valence());
            assert_eq!(molecule.stereo_elements().count(), 1);
            assert_eq!(
                molecule
                    .stereo_elements()
                    .next()
                    .expect("canonical stereo element")
                    .1
                    .is_specified(),
                expected_specified
            );

            let written = molfile::write_v2000(&molecule).expect("canonical stereo should project");
            let (reparsed, report) =
                read_molfile_with_report(&written).expect("projected stereo should re-interpret");
            assert_eq!(report.created_stereo_elements().len(), 1);
            assert_eq!(reparsed.stereo_elements().count(), 1);
            assert_eq!(
                reparsed
                    .atom(AtomId::new(0))
                    .expect("reparsed center")
                    .hydrogens,
                HydrogenDeclaration::Fixed(expected_hydrogens)
            );
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

    let undeclared = "undeclared hydrogen policy\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n";
    let molecule = read_molfile(undeclared).expect("undeclared V2000 atom interprets");
    assert_eq!(
        molecule.atom(AtomId::new(0)).expect("carbon").hydrogens,
        HydrogenDeclaration::Infer { explicit: 0 }
    );
}

#[test]
fn molfile_and_sdf_parse_supported_syntax_before_chemistry_interpretation() {
    let molfile_source = "unknown element\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 Xx  0  0  0  0  0  0\nM  END\n";
    let document = molfile::parse_str(molfile_source)
        .expect("a syntactically valid atom record should parse independently");
    assert!(molfile::interpret(&document)
        .expect_err("unsupported core elements belong to interpretation")
        .message()
        .contains("unsupported element"));

    let sdf_source = format!("{molfile_source}$$$$\n");
    let document =
        sdf::parse_str(&sdf_source).expect("SDF record structure should parse independently");
    let error =
        sdf::interpret(&document).expect_err("SDF delegates chemistry interpretation to Molfile");
    assert_eq!(error.record(), 1);
    assert!(error.message().contains("unsupported element"));
}

#[test]
fn v2000_rejects_unsupported_stereo_and_bond_representations() {
    for bond_line in ["  1  2  1  3  0  0  0", "  1  2  2  4  0  0  0"] {
        let input = format!(
                "bad stereo\nkekule\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0\n{bond_line}\nM  END\n"
            );
        assert!(read_molfile(&input).is_err());
    }

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
                orientation: Some(DoubleBondOrientation::Opposite),
            },
        )))
        .expect("double-bond stereo");
    assert!(molfile::write_v2000(molecule.working())
        .expect_err("specified double-bond stereo should be rejected")
        .message
        .contains("requires a Model"));

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
    assert!(molfile::write_v2000(molecule.working())
        .expect_err("quadruple bond should be rejected")
        .message
        .contains("quadruple"));
}

#[test]
fn mol_and_sdf_v2000_writers_round_trip_metadata_and_fields() {
    let input = "\
ammonium_acetate_like
kekule benchmark
M CHG and M ISO fixture
  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6000    0.7000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6000   -0.7000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  CHG  2   1   1   4  -1
M  ISO  1   2  13
M  END
>  <fixture_id>
charged_isotope_records

$$$$
";

    let records = read_sdf_records(input).expect("sdf should parse");
    let sdf = sdf::write_v2000(&records).expect("sdf should write");
    let reparsed = read_sdf_records(&sdf).expect("written sdf parses");

    assert_eq!(reparsed.len(), 1);
    assert_eq!(
        reparsed[0]
            .molecule()
            .atom(AtomId::new(0))
            .expect("atom")
            .formal_charge,
        1
    );
    assert_eq!(
        reparsed[0].data_fields()[0].value(),
        "charged_isotope_records"
    );
}

#[test]
fn v2000_charge_codes_and_chunked_metadata_round_trip_semantically() {
    for (charge_code, expected_charge) in
        [(1, 3), (2, 2), (3, 1), (0, 0), (5, -1), (6, -2), (7, -3)]
    {
        let input = format!(
                "charge\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 N   0  {charge_code}  0  0  0  0\nM  END\n"
            );
        let parsed = read_molfile(&input).expect("charge code should parse");
        assert_eq!(
            parsed.atom(AtomId::new(0)).expect("atom").formal_charge,
            expected_charge
        );
        let written = molfile::write_v2000(&parsed).expect("charge should write");
        let reparsed = read_molfile(&written).expect("charge should reparse");
        assert_eq!(
            reparsed.atom(AtomId::new(0)).expect("atom").formal_charge,
            expected_charge
        );
    }

    let mut graph_builder = crate::core::MoleculeEditor::new();
    let mut atom_ids = Vec::new();
    for index in 0..9u32 {
        let mut atom = carbon();
        atom.formal_charge = 1;
        atom.isotope = Some(13 + index as u16);
        atom.radical = Some(AtomRadical::Doublet);
        atom.atom_map = Some(index + 1);
        let atom_id = graph_builder
            .add_atom(atom)
            .expect("atom identifier capacity");
        if let Some(previous) = atom_ids.last().copied() {
            graph_builder
                .add_bond(previous, atom_id, BondOrder::Single)
                .expect("chain bond");
        }
        atom_ids.push(atom_id);
    }
    let mut molecule = graph_builder
        .finish()
        .expect("metadata fixture should be connected");
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
    molecule
        .insert_property(
            PropertyKey::new("sdf.field.NOTES").unwrap(),
            PropertyValue::String("line one\nline two".to_owned()),
        )
        .unwrap();
    let mol_text = molfile::write_v2000(&molecule).expect("metadata molecule should write");
    assert_eq!(mol_text.lines().nth(1), Some("kekule"));
    assert_eq!(mol_text.matches("M  CHG").count(), 2);
    assert_eq!(mol_text.matches("M  ISO").count(), 2);
    assert_eq!(mol_text.matches("M  RAD").count(), 2);

    let fields = vec![SdfDataField::new("NOTES", "line one\nline two")];
    let records = vec![
        SdfRecordInterpretation::new("metadata title", test_model(&molecule), fields.clone()),
        SdfRecordInterpretation::new("metadata title", test_model(&molecule), fields),
    ];
    let sdf_text = sdf::write_v2000(&records).expect("two records should write");
    assert_eq!(sdf_text.lines().nth(1), Some("kekule"));
    let records = read_sdf_records(&sdf_text).expect("written records should parse");
    assert_eq!(records.len(), 2);
    for record in records {
        assert_eq!(record.title(), "metadata title");
        assert_eq!(record.data_fields()[0].name(), "NOTES");
        assert_eq!(record.data_fields()[0].value(), "line one\nline two");
        for index in 0..9u32 {
            let atom = record.molecule().atom(AtomId::new(index)).expect("atom");
            assert_eq!(atom.formal_charge, 1);
            assert_eq!(atom.isotope, Some(13 + index as u16));
            assert_eq!(atom.radical, Some(AtomRadical::Doublet));
            assert_eq!(atom.atom_map, Some(index + 1));
        }
    }
}
