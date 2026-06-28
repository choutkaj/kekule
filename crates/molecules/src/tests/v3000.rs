use super::*;

#[test]
fn mol_v3000_parses_raw_atoms_bonds_coordinates_and_metadata() {
    let input = "\
charged radical
molecules validation
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

    let small = read_mol_v3000_str(input).expect("V3000 should parse");
    let mol = &small.mol;

    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
    assert_eq!(
        mol.props().get("sdf.title"),
        Some(&PropValue::String("charged radical".to_owned()))
    );
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
    assert_eq!(bond0.stereo, Some(BondStereo::Up));
    assert_eq!(bond1.order, BondOrder::Double);
    let (_, conformer) = mol.first_conformer().expect("conformer exists");
    assert_eq!(
        conformer.position(AtomId::new(0)),
        Some(Point3::new(0.1, 0.2, 0.3))
    );
}

#[test]
fn mol_v3000_line_continuations_and_aromatic_bonds_parse_without_perception() {
    let input = "\
benzene-ish
molecules

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

    let small = read_mol_v3000_str(input).expect("V3000 should parse");
    let mol = &small.mol;

    assert_eq!(
        mol.bond(BondId::new(0)).expect("bond").order,
        BondOrder::Aromatic
    );
    assert_eq!(mol.perception().rings, ComputedState::Absent);
    assert_eq!(mol.perception().aromaticity, ComputedState::Absent);
}

#[test]
fn malformed_mol_v3000_returns_errors_without_panicking() {
    let cases = [
        (
            "missing V3000 declaration",
            "Bad\nmolecules\n\n  0  0  0  0  0  0            999 V2000\nM  END\n",
        ),
        (
            "bad counts",
            "Bad\nmolecules\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS nope 0 0 0 0\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "count mismatch",
            "Bad\nmolecules\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "non-finite coordinates",
            "Bad\nmolecules\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 1e999 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "bad endpoint",
            "Bad\nmolecules\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "unsupported atom stereo",
            "Bad\nmolecules\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 CFG=1\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
        (
            "unsupported bond type",
            "Bad\nmolecules\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 8 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ),
    ];

    for (name, input) in cases {
        let parsed = std::panic::catch_unwind(|| read_mol_v3000_str(input))
            .unwrap_or_else(|_| panic!("{name} panicked"));
        let error = parsed.expect_err("malformed V3000 input should fail");
        assert!(error.line >= 1, "line for {name}");
        assert!(!error.message.is_empty(), "message for {name}");
    }
}
