use super::*;

#[test]
fn mol_v3000_parses_raw_atoms_bonds_coordinates_and_metadata() {
    let input = "\
charged radical
molecular validation
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
M  V30 1 1 1 2 CFG=1
M  V30 2 2 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
";

    let small = read_molfile(input).expect("V3000 should parse");
    let mol = small.graph();

    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
    assert!(mol.props().get("sdf.title").is_none());
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
    assert_eq!(
        mol.stereo_bond_mark(BondId::new(0))
            .expect("stereo mark")
            .kind,
        StereoBondMarkKind::WedgeUp
    );
    assert_eq!(bond1.order, BondOrder::Double);
    let (_, conformer) = mol.first_conformer().expect("conformer exists");
    assert_eq!(
        conformer.position(AtomId::new(0)),
        Some(crate::units::Quantity::new(
            Point3::new(0.1, 0.2, 0.3),
            crate::units::ANGSTROM,
        ))
    );
}

#[test]
fn v3000_preserves_valence_implied_tetrahedral_hydrogen_carrier() {
    let input = "\
stereo hydrogen
molecular

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 F 1 0 0 0
M  V30 3 Cl -1 0 0 0
M  V30 4 Br 0 1 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
";

    let parsed = read_molfile(input).expect("V3000 should parse");

    assert_eq!(
        parsed
            .graph()
            .atom(AtomId::new(0))
            .expect("stereo center")
            .explicit_hydrogens,
        1
    );
}

#[test]
fn mol_v3000_line_continuations_and_aromatic_bonds_parse_without_perception() {
    let input = "\
benzene-ish
molecular

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
    let mol = small.graph();

    assert_eq!(
        mol.bond(BondId::new(0)).expect("bond").order,
        BondOrder::Aromatic
    );
    assert!(!mol.perception().has_rings());
    assert!(!mol.perception().has_aromaticity());
}

#[test]
fn malformed_mol_v3000_returns_errors_without_panicking() {
    let cases = [
        (
            "bad counts",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS nope 0 0 0 0\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "count mismatch",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "non-finite coordinates",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 1e999 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "bad endpoint",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "unsupported atom stereo",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 CFG=1\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "unsupported bond type",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 8 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "incomplete counts",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "zero atom index",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 0 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate bond index",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 3 C 2 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 1 1 2 3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate counts",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "counts after atom section",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate atom section",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 BEGIN ATOM\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate bond section",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "record outside CTAB",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 NOTE=outside\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "unsupported atom option",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 VAL=4\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "malformed atom option",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 BROKEN\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate atom option",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 CHG=1 CHG=2\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "unsupported bond option",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 TOPO=1\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "duplicate bond option",
            "Bad\nmolecular\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 CFG=1 CFG=1\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
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
molecular

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 1)
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
molecular

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

    let mut nitrogen = Atom::new(Element::from_symbol("N").expect("N"));
    nitrogen.formal_charge = 1;
    nitrogen.radical = Some(AtomRadical::Doublet);
    nitrogen.atom_map = Some(42);
    let n = molecule.graph_mut().add_atom(nitrogen);

    let mut carbon = carbon();
    carbon.isotope = Some(13);
    let c = molecule.graph_mut().add_atom(carbon);

    let mut oxygen = oxygen();
    oxygen.formal_charge = -1;
    let o = molecule.graph_mut().add_atom(oxygen);

    let wedge = molecule
        .graph_mut()
        .add_bond(n, c, BondOrder::Single)
        .expect("single bond");
    molecule
        .graph_mut()
        .set_stereo_bond_mark(StereoBondMark {
            bond: wedge,
            kind: StereoBondMarkKind::WedgeUp,
            source: StereoSource::MolfileV3000,
        })
        .expect("stereo mark");
    molecule
        .graph_mut()
        .add_bond(c, o, BondOrder::Double)
        .expect("double bond");

    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            n,
            crate::units::Quantity::new(Point3::new(0.1, 0.2, 0.3), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            c,
            crate::units::Quantity::new(Point3::new(1.4, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            o,
            crate::units::Quantity::new(Point3::new(2.5, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    molecule
        .graph_mut()
        .add_conformer(conformer)
        .expect("valid conformer");

    let written = molfile::write_v3000(&molecule).expect("V3000 should write");
    assert_eq!(written.lines().nth(1), Some("molecular"));
    assert!(written.contains("V3000"));
    assert!(written.contains("CHG=1"));
    assert!(written.contains("MASS=13"));
    assert!(written.contains("RAD=2"));
    assert!(written.contains("CFG=1"));

    let reparsed = read_molfile(&written).expect("written V3000 should parse");
    assert!(reparsed.graph().props().get("sdf.title").is_none());
    assert_eq!(
        reparsed
            .graph()
            .atom(AtomId::new(0))
            .expect("atom")
            .formal_charge,
        1
    );
    assert_eq!(
        reparsed.graph().atom(AtomId::new(0)).expect("atom").radical,
        Some(AtomRadical::Doublet)
    );
    assert_eq!(
        reparsed
            .graph()
            .atom(AtomId::new(0))
            .expect("atom")
            .atom_map,
        Some(42)
    );
    assert_eq!(
        reparsed.graph().atom(AtomId::new(1)).expect("atom").isotope,
        Some(13)
    );
    assert_eq!(
        reparsed
            .graph()
            .stereo_bond_mark(BondId::new(0))
            .expect("stereo mark")
            .kind,
        StereoBondMarkKind::WedgeUp
    );
    let (_, conformer) = reparsed.graph().first_conformer().expect("conformer");
    assert_eq!(
        conformer.position(AtomId::new(2)),
        Some(crate::units::Quantity::new(
            Point3::new(2.5, 0.0, 0.0),
            crate::units::ANGSTROM,
        ))
    );
}

#[test]
fn mol_v3000_writer_rejects_unsupported_stereo_and_bonds() {
    let mut molecule = SmallMolecule::default();
    let a = molecule.graph_mut().add_atom(carbon());
    molecule
        .graph_mut()
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center: a,
                carriers: vec![StereoCarrier::ImplicitHydrogen],
                orientation: TetrahedralOrientation::Clockwise,
            }),
            StereoSource::User,
        ))
        .expect("stereo element");
    assert!(molfile::write_v3000(&molecule)
        .expect_err("stereo elements should be rejected")
        .message
        .contains("stereo elements"));

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
                orientation: DoubleBondOrientation::Together,
            }),
            StereoSource::User,
        ))
        .expect("double-bond stereo");
    assert!(molfile::write_v3000(&molecule)
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
            source: StereoSource::MolfileV3000,
        })
        .expect("stereo mark");
    assert!(molfile::write_v3000(&molecule)
        .expect_err("double wedge should be rejected")
        .message
        .contains("incompatible"));

    molecule
        .graph_mut()
        .clear_stereo_bond_mark(bond)
        .expect("clear mark");
    molecule.graph_mut().bond_mut(bond).expect("bond").order = BondOrder::Quadruple;
    assert!(molfile::write_v3000(&molecule)
        .expect_err("quadruple should be rejected")
        .message
        .contains("quadruple"));
}
