use super::*;
use crate::properties::{PropertyKey, PropertyValue};

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
    assert_eq!(interpretation.report().records().len(), 1);
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
fn v2000_does_not_infer_tetrahedral_hydrogens_without_a_source_declaration() {
    for symbol in ["C", "N", "S"] {
        let input = format!(
            "stereo hydrogen\nkekule\n\n  4  3  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 {symbol:<3} 0  0  0  0  0  0\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0\n    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0\n  1  2  1  1  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\nM  END\n"
        );

        if symbol == "S" {
            let molecule = read_molfile(&input).expect("sulfur lone-pair stereo interprets");
            assert_eq!(molecule.stereo_elements().count(), 1);
            assert!(molecule.stereo_elements().any(|(_, element)| {
                matches!(
                    &element.kind,
                    StereoElementKind::Tetrahedral(stereo)
                        if stereo.carriers.contains(&StereoCarrier::ImplicitLonePair)
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
            let molecule = interpreted.to_molecule();
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
    molecule
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond,
                left: a,
                right: b,
                left_carrier: StereoCarrier::Atom(a),
                right_carrier: StereoCarrier::Atom(b),
                orientation: Some(DoubleBondOrientation::Opposite),
            },
        )))
        .expect("double-bond stereo");
    assert!(molfile::write_v2000(molecule.working())
        .expect_err("invalid stereo element should be rejected")
        .message
        .contains("cannot encode"));

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
