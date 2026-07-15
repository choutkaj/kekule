use super::*;

#[test]
fn molfile_and_sdf_documents_preserve_record_metadata_before_interpretation() {
    let molfile_text = "Header title\nprogram line\ncomment line\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nX  UNSUPPORTED\nM  END\n";
    let document = molfile::parse_str(molfile_text).expect("Molfile document parses");
    assert_eq!(document.header().title(), "Header title");
    assert_eq!(document.unsupported_records().len(), 1);
    let molecule = molfile::interpret(&document).expect("Molfile interprets");
    assert!(molecule.graph().props().get("sdf.title").is_none());

    let sdf_text = format!("{molfile_text}>  <FIELD>\nvalue\n\n$$$$\n");
    let document = sdf::parse_str(&sdf_text, SdfParseOptions::default()).expect("SDF parses");
    assert_eq!(document.records()[0].data_fields()[0].value(), "value");
    let records = sdf::interpret(&document).expect("SDF interprets");
    assert_eq!(records[0].title(), "Header title");
    assert_eq!(records[0].data_fields()[0].name(), "FIELD");
    assert!(records[0]
        .molecule()
        .graph()
        .props()
        .get("sdf.field.FIELD")
        .is_none());
}

#[test]
fn molfile_and_sdf_documents_parse_adjacent_three_digit_counts() {
    let mut molfile_text =
        String::from("Large\nprogram\ncomment\n999999  0  0  0  0            999 V2000\n");
    for _ in 0..999 {
        molfile_text.push_str("atom record\n");
    }
    for _ in 0..999 {
        molfile_text.push_str("bond record\n");
    }
    molfile_text.push_str("M  END\n");

    let document = molfile::parse_str(&molfile_text).expect("fixed-width counts parse");
    assert_eq!(document.atom_records().len(), 999);
    assert_eq!(document.bond_records().len(), 999);

    let sdf_text = format!("{molfile_text}$$$$\n");
    let document = sdf::parse_str(&sdf_text, SdfParseOptions::default())
        .expect("SDF delegates to fixed-width Molfile counts parsing");
    assert_eq!(document.records()[0].molfile().atom_records().len(), 999);
    assert_eq!(document.records()[0].molfile().bond_records().len(), 999);
}

#[test]
fn sdf_v2000_parses_single_record_atoms_bonds_and_fields() {
    let input = "\
Water
  molecules
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
    let mol = records[0].molecule().graph();

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
  molecules

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
$$$$
Two
  molecules

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
            .graph()
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
  molecules

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
";

    let molecules = read_sdf_molecules_with_options(
        input,
        SdfParseOptions {
            allow_missing_final_delimiter: true,
        },
    )
    .expect("record should parse");

    assert_eq!(molecules.len(), 1);
    assert_eq!(molecules[0].graph().atom_count(), 1);
}

#[test]
fn sdf_v2000_rejects_v3000_and_bad_endpoints() {
    let v3000 = "\
V3000
  molecules

  0  0  0  0  0  0            999 V3000
M  END
$$$$
";
    let err = read_sdf_molecules(v3000).expect_err("V3000 should fail");
    assert!(!err.to_string().is_empty());

    let bad_endpoint = "\
Bad
  molecules

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
                "Bad\nmolecules\n\n  1  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n  0  1  1  0  0  0  0\nM  END\n",
            ),
            (
                "non-ASCII counts",
                "Bad\nmolecules\n\né  1  0  0  0  0            999 V2000\nM  END\n",
            ),
            (
                "non-ASCII atom",
                "Bad\nmolecules\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 Cé  0  0  0  0  0  0\nM  END\n",
            ),
            (
                "truncated atom",
                "Bad\nmolecules\n\n  1  0  0  0  0  0            999 V2000\n0.0 C\nM  END\n",
            ),
            (
                "non-ASCII bond",
                "Bad\nmolecules\n\n  1  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n  1  é  1  0\nM  END\n",
            ),
            (
                "count over format limit",
                "Bad\nmolecules\n\n1000 0 V2000\nM  END\n",
            ),
            (
                "inconsistent counts",
                "Bad\nmolecules\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n",
            ),
            (
                "truncated M record",
                "Bad\nmolecules\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  CHG  2   1   1\nM  END\n",
            ),
            (
                "zero M-record atom",
                "Bad\nmolecules\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  CHG  1   0   1\nM  END\n",
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
fn sdf_v2000_parse_does_not_perceive_chemistry() {
    let input = "\
Benzene-ish
  molecules

  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0
  1  2  4  0  0  0  0
M  END
$$$$
";

    let molecules = read_sdf_molecules(input).expect("record should parse");
    let mol = &molecules[0].graph();

    assert!(!mol.perception().has_rings());
    assert!(!mol.perception().has_aromaticity());
    assert_eq!(
        mol.bond(BondId::new(0)).expect("bond exists").order,
        BondOrder::Aromatic
    );
}

#[test]
fn mol_v2000_preserves_coordinates_charges_isotopes_radicals_and_atom_maps() {
    let input = "\
charged radical
molecules validation
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
    let atom0 = small.graph().atom(AtomId::new(0)).expect("atom exists");
    let atom1 = small.graph().atom(AtomId::new(1)).expect("atom exists");
    assert_eq!(atom0.formal_charge, 1);
    assert_eq!(atom0.radical, Some(AtomRadical::Doublet));
    assert_eq!(atom0.atom_map, Some(7));
    assert_eq!(atom1.isotope, Some(13));
    let (_, conformer) = small.graph().first_conformer().expect("conformer exists");
    assert_eq!(
        conformer.position(AtomId::new(0)),
        Some(Point3::new(0.1, 0.2, 0.3))
    );
}

#[test]
fn v2000_radical_codes_round_trip_exact_multiplicity() {
    for (code, expected) in [
        (1, AtomRadical::Singlet),
        (2, AtomRadical::Doublet),
        (3, AtomRadical::Triplet),
    ] {
        let input = format!(
                "radical {code}\nmolecules\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  RAD  1   1   {code}\nM  END\n"
            );
        let parsed = read_molfile(&input).expect("radical record should parse");
        assert_eq!(
            parsed.graph().atom(AtomId::new(0)).expect("atom").radical,
            Some(expected)
        );

        let written = molfile::write_v2000(&parsed).expect("radical record should write");
        assert!(
            written.contains(&format!("M  RAD  1   1   {code}")),
            "written code {code}: {written}"
        );
        let reparsed = read_molfile(&written).expect("written radical record should parse");
        assert_eq!(
            reparsed.graph().atom(AtomId::new(0)).expect("atom").radical,
            Some(expected)
        );
    }
}

#[test]
fn v2000_supported_bond_stereo_codes_round_trip_exactly() {
    for (order_code, stereo_code, order, expected) in [
        (1, 1, BondOrder::Single, StereoBondMarkKind::WedgeUp),
        (1, 4, BondOrder::Single, StereoBondMarkKind::WedgeEither),
        (1, 6, BondOrder::Single, StereoBondMarkKind::WedgeDown),
        (
            2,
            3,
            BondOrder::Double,
            StereoBondMarkKind::DoubleBondEither,
        ),
    ] {
        let input = format!(
                "stereo\nmolecules\n\n  2  1  0  0  0  0            999 V2000\n   -1.2500    0.0000    0.0000 C   0  0  0  0  0  0\n    1.2500    0.0000    0.0000 C   0  0  0  0  0  0\n  1  2  {order_code}  {stereo_code}  0  0  0\nM  END\n"
            );
        let parsed = read_molfile(&input).expect("stereo record should parse");
        let bond = parsed.graph().bond(BondId::new(0)).expect("bond");
        assert_eq!(bond.order, order);
        assert_eq!(
            parsed
                .graph()
                .stereo_bond_mark(BondId::new(0))
                .expect("stereo mark")
                .kind,
            expected
        );

        let written = molfile::write_v2000(&parsed).expect("stereo record should write");
        let reparsed = read_molfile(&written).expect("written stereo record should parse");
        assert_eq!(
            reparsed
                .graph()
                .stereo_bond_mark(BondId::new(0))
                .expect("stereo mark")
                .kind,
            expected
        );
        assert_eq!(
            reparsed
                .graph()
                .first_conformer()
                .expect("conformer")
                .1
                .position(AtomId::new(0)),
            Some(Point3::new(-1.25, 0.0, 0.0))
        );
    }
}

#[test]
fn v2000_preserves_valence_implied_tetrahedral_hydrogen_carriers() {
    for (symbol, expected_hydrogens) in [("C", 1), ("N", 0), ("S", 1)] {
        let input = format!(
            "stereo hydrogen\nmolecules\n\n  4  3  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 {symbol:<3} 0  0  0  0  0  0\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0\n    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0\n  1  2  1  1  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\nM  END\n"
        );

        let parsed = read_molfile(&input).expect("stereo record should parse");

        assert_eq!(
            parsed
                .graph()
                .atom(AtomId::new(0))
                .expect("stereo center")
                .explicit_hydrogens,
            expected_hydrogens,
            "{symbol}"
        );
    }
}

#[test]
fn v2000_rejects_unsupported_stereo_and_bond_representations() {
    for bond_line in ["  1  2  1  3  0  0  0", "  1  2  2  4  0  0  0"] {
        let input = format!(
                "bad stereo\nmolecules\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0\n{bond_line}\nM  END\n"
            );
        assert!(read_molfile(&input).is_err());
    }

    let mut molecule = SmallMolecule::default();
    let a = molecule.graph_mut().add_atom(carbon());
    let b = molecule.graph_mut().add_atom(carbon());
    let bond = molecule
        .graph_mut()
        .add_bond(a, b, BondOrder::Double)
        .expect("bond");
    molecule
        .graph_mut()
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::DoubleBond(DoubleBondStereo {
                bond,
                left: a,
                right: b,
                left_carrier: StereoCarrier::Atom(a),
                right_carrier: StereoCarrier::Atom(b),
                orientation: DoubleBondOrientation::Opposite,
            }),
            StereoSource::User,
        ))
        .expect("double-bond stereo");
    assert!(molfile::write_v2000(&molecule)
        .expect_err("stereo elements should be rejected")
        .message
        .contains("stereo elements"));

    let element = molecule
        .graph()
        .stereo_element_ids()
        .next()
        .expect("stereo element");
    molecule
        .graph_mut()
        .remove_stereo_element(element)
        .expect("remove stereo element");
    molecule
        .graph_mut()
        .set_stereo_bond_mark(StereoBondMark {
            bond,
            kind: StereoBondMarkKind::WedgeUp,
            source: StereoSource::MolfileV2000,
        })
        .expect("mark");
    assert!(molfile::write_v2000(&molecule)
        .expect_err("double wedge should be rejected")
        .message
        .contains("incompatible"));

    molecule
        .graph_mut()
        .clear_stereo_bond_mark(bond)
        .expect("clear mark");
    molecule.graph_mut().bond_mut(bond).expect("bond").order = BondOrder::Quadruple;
    assert!(molfile::write_v2000(&molecule)
        .expect_err("quadruple bond should be rejected")
        .message
        .contains("quadruple"));
}

#[test]
fn mol_and_sdf_v2000_writers_round_trip_metadata_and_fields() {
    let input = "\
ammonium_acetate_like
molecules validation
M CHG and M ISO fixture
  4  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6000    0.7000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6000   -0.7000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
            .graph()
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
                "charge\nmolecules\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 N   0  {charge_code}  0  0  0  0\nM  END\n"
            );
        let parsed = read_molfile(&input).expect("charge code should parse");
        assert_eq!(
            parsed
                .graph()
                .atom(AtomId::new(0))
                .expect("atom")
                .formal_charge,
            expected_charge
        );
        let written = molfile::write_v2000(&parsed).expect("charge should write");
        let reparsed = read_molfile(&written).expect("charge should reparse");
        assert_eq!(
            reparsed
                .graph()
                .atom(AtomId::new(0))
                .expect("atom")
                .formal_charge,
            expected_charge
        );
    }

    let mut molecule = SmallMolecule::default();
    molecule.graph_mut().props_mut().insert(
        "sdf.title".to_owned(),
        PropValue::String("metadata title".to_owned()),
    );
    molecule.graph_mut().props_mut().insert(
        "sdf.program".to_owned(),
        PropValue::String("metadata program".to_owned()),
    );
    molecule.graph_mut().props_mut().insert(
        "sdf.comment".to_owned(),
        PropValue::String("metadata comment".to_owned()),
    );
    molecule.graph_mut().props_mut().insert(
        "sdf.field.NOTES".to_owned(),
        PropValue::String("line one\nline two".to_owned()),
    );
    let mut conformer = Conformer::new();
    for index in 0..9u32 {
        let mut atom = carbon();
        atom.formal_charge = 1;
        atom.isotope = Some(13 + index as u16);
        atom.radical = Some(AtomRadical::Doublet);
        atom.atom_map = Some(index + 1);
        let atom_id = molecule.graph_mut().add_atom(atom);
        conformer.set_position(atom_id, Point3::new(-(index as f64), index as f64, 0.0));
    }
    molecule
        .graph_mut()
        .add_conformer(conformer)
        .expect("valid conformer");

    let mol_text = molfile::write_v2000(&molecule).expect("metadata molecule should write");
    assert_eq!(mol_text.matches("M  CHG").count(), 2);
    assert_eq!(mol_text.matches("M  ISO").count(), 2);
    assert_eq!(mol_text.matches("M  RAD").count(), 2);

    let fields = vec![SdfDataField::new("NOTES", "line one\nline two")];
    let records = vec![
        SdfRecord::new("metadata title", molecule.clone(), fields.clone()),
        SdfRecord::new("metadata title", molecule, fields),
    ];
    let sdf_text = sdf::write_v2000(&records).expect("two records should write");
    let records = read_sdf_records(&sdf_text).expect("written records should parse");
    assert_eq!(records.len(), 2);
    for record in records {
        assert_eq!(record.title(), "metadata title");
        assert_eq!(record.data_fields()[0].name(), "NOTES");
        assert_eq!(record.data_fields()[0].value(), "line one\nline two");
        for index in 0..9u32 {
            let atom = record
                .molecule()
                .graph()
                .atom(AtomId::new(index))
                .expect("atom");
            assert_eq!(atom.formal_charge, 1);
            assert_eq!(atom.isotope, Some(13 + index as u16));
            assert_eq!(atom.radical, Some(AtomRadical::Doublet));
            assert_eq!(atom.atom_map, Some(index + 1));
        }
    }
}
