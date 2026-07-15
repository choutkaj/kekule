use super::*;

#[test]
fn valence_and_sanitization_are_explicit() {
    let mut small = read_smiles("CCO").expect("smiles should parse");
    assert!(!small.graph().perception().has_valence());

    let report = perception_api::sanitize_with_options(&mut small, SanitizeOptions::default())
        .expect("ethanol should sanitize");

    assert!(report.valence.expect("valence report").is_ok());
    assert!(small.graph().perception().has_valence());
    assert!(small.graph().perception().has_rings());
    assert_eq!(
        small
            .graph()
            .atom(AtomId::new(2))
            .expect("oxygen")
            .implicit_hydrogens,
        Some(1)
    );
}

#[test]
fn sanitize_options_do_not_leave_skipped_passes_fresh() {
    let mut baseline = read_smiles("C1=CC=CC=C1").expect("benzene should parse");
    perception_api::sanitize_with_options(&mut baseline, SanitizeOptions::default())
        .expect("benzene should sanitize");

    for mask in 0..16 {
        let options = SanitizeOptions {
            perceive_valence: mask & 1 != 0,
            perceive_rings: mask & 2 != 0,
            perceive_aromaticity: mask & 4 != 0,
            perceive_stereo: mask & 8 != 0,
        };
        let mut molecule = baseline.clone();
        perception_api::sanitize_with_options(&mut molecule, options)
            .unwrap_or_else(|error| panic!("options {mask:04b} should succeed: {error}"));

        assert_eq!(
            molecule.graph().perception().has_valence(),
            options.perceive_valence,
            "valence state for options {mask:04b}"
        );
        assert_eq!(
            molecule.graph().perception().has_rings(),
            options.perceive_rings,
            "ring state for options {mask:04b}"
        );
        assert_eq!(
            molecule.graph().perception().has_aromaticity(),
            options.perceive_aromaticity,
            "aromaticity state for options {mask:04b}"
        );
        assert_eq!(
            molecule.graph().ring_set().is_some(),
            options.perceive_rings,
            "ring cache exposure for options {mask:04b}"
        );
    }
}

#[test]
fn sanitization_perceives_stereo_by_default_and_can_skip_it() {
    let mut molecule = read_smiles("C/C=C\\F").expect("directional smiles should parse");

    let report = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("directional molecule should sanitize");

    assert!(report.stereo.expect("stereo report").is_ok());
    assert_eq!(molecule.graph().stereo_elements().count(), 1);

    let mut skipped = read_smiles("C/C=C\\F").expect("directional smiles should parse");
    perception_api::sanitize_with_options(
        &mut skipped,
        SanitizeOptions {
            perceive_stereo: false,
            ..SanitizeOptions::default()
        },
    )
    .expect("stereo-skipped molecule should sanitize");

    assert!(!skipped.graph().perception().has_cip_descriptors());
    assert_eq!(skipped.graph().stereo_elements().count(), 0);
}

#[test]
fn sanitization_preserves_unknown_double_bond_stereo() {
    let mut molecule = read_smiles("CC=CC").expect("alkene should parse");
    let double_bond = molecule
        .graph()
        .bonds()
        .find_map(|(bond_id, bond)| (bond.order == BondOrder::Double).then_some(bond_id))
        .expect("double bond");
    molecule
        .graph_mut()
        .set_stereo_bond_mark(StereoBondMark {
            bond: double_bond,
            kind: StereoBondMarkKind::DoubleBondEither,
            source: StereoSource::MolfileV2000,
        })
        .expect("double bond either mark");

    let report = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("unknown double-bond stereo should sanitize");

    assert!(report.stereo.expect("stereo report").is_ok());
    let (_, element) = molecule
        .graph()
        .stereo_elements()
        .next()
        .expect("unknown stereo element");
    assert_eq!(element.specifiedness, StereoSpecifiedness::Unknown);
    assert!(matches!(
        &element.kind,
        StereoElementKind::DoubleBond(stereo) if stereo.bond == double_bond
    ));
}

#[test]
fn sanitization_does_not_assign_coordinate_only_stereo() {
    let mut mol = Molecule::new();
    let left = mol.add_atom(carbon());
    let right = mol.add_atom(carbon());
    let left_carrier = mol.add_atom(element_atom("F"));
    let right_carrier = mol.add_atom(element_atom("Cl"));
    mol.add_bond(left, right, BondOrder::Double).expect("bond");
    mol.add_bond(left, left_carrier, BondOrder::Single)
        .expect("left carrier");
    mol.add_bond(right, right_carrier, BondOrder::Single)
        .expect("right carrier");
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            left,
            crate::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            right,
            crate::units::Quantity::new(Point3::new(1.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            left_carrier,
            crate::units::Quantity::new(Point3::new(0.0, 1.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            right_carrier,
            crate::units::Quantity::new(Point3::new(1.0, -1.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    mol.add_conformer(conformer).expect("valid conformer");
    let mut molecule = SmallMolecule::from_graph(mol);

    let report = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("coordinate-only molecule should sanitize");

    assert!(report.stereo.expect("stereo report").is_ok());
    assert_eq!(molecule.graph().stereo_elements().count(), 0);

    let direct_report = stereo_api::perceive_stereo(molecule.graph_mut());
    assert!(direct_report.is_ok(), "{:?}", direct_report.issues);
    assert_eq!(direct_report.created_elements.len(), 1);
    assert_eq!(molecule.graph().stereo_elements().count(), 1);
}

#[test]
fn failed_stereo_sanitization_is_transactional() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    let bond = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    mol.set_stereo_bond_mark(StereoBondMark {
        bond,
        kind: StereoBondMarkKind::WedgeEither,
        source: StereoSource::MolfileV2000,
    })
    .expect("wedge mark");
    let mut molecule = SmallMolecule::from_graph(mol);
    let before = molecule.clone();

    let error = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect_err("unassembled stereo mark should fail sanitization");

    assert!(matches!(error, SanitizeError::Stereo(_)));
    assert_eq!(molecule, before);
}

#[test]
fn sanitization_treats_conflicting_wedges_as_nonfatal_ambiguity() {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon());
    let mut marked_bonds = Vec::new();
    for symbol in ["F", "Cl", "Br", "I"] {
        let carrier = mol.add_atom(element_atom(symbol));
        marked_bonds.push(
            mol.add_bond(center, carrier, BondOrder::Single)
                .expect("carrier bond"),
        );
    }
    for (index, bond) in marked_bonds.into_iter().enumerate() {
        mol.set_stereo_bond_mark(StereoBondMark {
            bond,
            kind: if index % 2 == 0 {
                StereoBondMarkKind::WedgeUp
            } else {
                StereoBondMarkKind::WedgeDown
            },
            source: StereoSource::MolfileV2000,
        })
        .expect("wedge mark");
    }
    let mut molecule = SmallMolecule::from_graph(mol);

    let report = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("ambiguous drawing wedges should not reject valid chemistry");
    let stereo = report.stereo.expect("stereo report");

    assert!(stereo
        .issues
        .contains(&StereoPerceptionIssue::AmbiguousTetrahedralWedgeMarks {
            center,
            mark_count: 4,
        }));
    assert_eq!(stereo.issues.len(), 1);
    assert!(molecule.graph().stereo_elements().next().is_none());
}

#[test]
fn failed_valence_sanitization_is_transactional() {
    let mut mol = Molecule::new();
    let carbon = mol.add_atom(carbon());
    for _ in 0..5 {
        let hydrogen = mol.add_atom(Atom::new(Element::from_symbol("H").expect("hydrogen")));
        mol.add_bond(carbon, hydrogen, BondOrder::Single)
            .expect("bond");
    }
    rings_api::perceive_ring_set(&mut mol).expect("ring perception should succeed");
    let mut molecule = SmallMolecule::from_graph(mol);
    let before = molecule.clone();

    let error = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect_err("pentavalent carbon should fail valence");

    assert!(matches!(error, SanitizeError::Valence(_)));
    assert_eq!(molecule, before);
}

#[test]
fn failed_aromaticity_sanitization_is_transactional() {
    let mut molecule = read_smiles("c1cccc1").expect("raw invalid aromatic representation parses");
    let before = molecule.clone();

    let error = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect_err("unmatchable aromatic representation should fail");

    assert!(matches!(error, SanitizeError::Aromaticity(_)));
    assert_eq!(molecule, before);
}

#[test]
fn failed_direct_aromaticity_perception_is_transactional() {
    let mut molecule = read_smiles("c1cccc1").expect("raw invalid aromatic representation parses");
    let before = molecule.graph().clone();

    let error =
        aromaticity_api::perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
            .expect_err("unmatchable aromatic representation should fail");

    assert!(matches!(
        error,
        AromaticityError::InvalidAromaticRepresentation(_)
    ));
    assert_eq!(molecule.graph(), &before);
}

#[test]
fn successful_sanitization_is_idempotent() {
    let mut molecule = read_smiles("CCO").expect("ethanol should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("first sanitize should succeed");
    let once = molecule.clone();

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("second sanitize should succeed");

    assert_eq!(molecule, once);
}

#[test]
fn sanitize_cleanup_invalidates_preexisting_perception() {
    let mut mol = Molecule::new();
    let chlorine = mol.add_atom(Atom::new(Element::from_symbol("Cl").expect("chlorine")));
    let oxo = mol.add_atom(oxygen());
    let hydroxyl = mol.add_atom(oxygen());
    mol.add_bond(chlorine, oxo, BondOrder::Double)
        .expect("bond");
    mol.add_bond(chlorine, hydroxyl, BondOrder::Single)
        .expect("hydroxyl bond");
    mark_all_fresh(&mut mol);
    let mut molecule = SmallMolecule::from_graph(mol);

    perception_api::sanitize_with_options(
        &mut molecule,
        SanitizeOptions {
            perceive_valence: false,
            perceive_rings: false,
            perceive_aromaticity: false,
            perceive_stereo: false,
        },
    )
    .expect("cleanup-only sanitize should succeed");

    assert_all_stale(molecule.graph());
    assert_eq!(
        molecule
            .graph()
            .atom(chlorine)
            .expect("chlorine")
            .formal_charge,
        1
    );
    assert_eq!(
        molecule.graph().atom(oxo).expect("oxygen").formal_charge,
        -1
    );
    assert_eq!(
        molecule
            .graph()
            .atom(hydroxyl)
            .expect("hydroxyl oxygen")
            .formal_charge,
        0
    );
    assert_eq!(
        molecule.graph().bond(BondId::new(0)).expect("bond").order,
        BondOrder::Single
    );
}

#[test]
fn valence_reports_excess_common_valence() {
    let mut mol = Molecule::new();
    let c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("C")));
    for _ in 0..5 {
        let h = mol.add_atom(Atom::new(Element::from_symbol("H").expect("H")));
        mol.add_bond(c, h, BondOrder::Single).expect("bond");
    }

    let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

    assert_eq!(report.issues.len(), 1);
    assert!(!report.is_ok());

    let report = valence_api::perceive_valence_with_options(
        &mut mol,
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    );
    assert!(report.is_ok());
    assert_eq!(mol.atom(c).expect("carbon").implicit_hydrogens, Some(0));
}

#[test]
fn valence_counts_high_degree_atoms_without_narrowing_or_panicking() {
    let mut mol = Molecule::new();
    let carbon = mol.add_atom(element_atom("C"));
    for _ in 0..300 {
        let hydrogen = mol.add_atom(element_atom("H"));
        mol.add_bond(carbon, hydrogen, BondOrder::Single)
            .expect("bond");
    }

    let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

    assert_eq!(
        report.issues,
        vec![ValenceIssue::ValenceExceeded {
            atom: carbon,
            explicit_valence: 300,
            max_allowed: 4,
        }]
    );
    assert_eq!(
        mol.atom(carbon).expect("carbon").implicit_hydrogens,
        Some(0)
    );
}

#[test]
fn valence_uses_rdkit_periodic_table_rules_for_electropositive_atoms() {
    for (symbol, expected_implicit_hydrogens) in [
        ("Li", 1),
        ("Be", 2),
        ("Na", 1),
        ("Mg", 2),
        ("K", 1),
        ("Ca", 2),
        ("Rb", 1),
        ("Sr", 2),
        ("Cs", 1),
        ("Ba", 2),
        ("Fr", 1),
        ("Ra", 2),
    ] {
        let mut mol = Molecule::new();
        let atom_id = mol.add_atom(element_atom(symbol));

        let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

        assert!(report.is_ok(), "neutral {symbol} should be supported");
        assert_eq!(
            mol.atom(atom_id).expect("atom").implicit_hydrogens,
            Some(expected_implicit_hydrogens),
            "neutral {symbol} implicit hydrogens"
        );
    }
}

#[test]
fn valence_keeps_rdkit_hypervalent_anion_limits() {
    for (symbol, charge, accepted, rejected) in [
        ("P", -2, 3, 4),
        ("S", -1, 5, 6),
        ("As", -2, 3, 4),
        ("Se", -1, 5, 6),
    ] {
        let mut accepted_mol = Molecule::new();
        let accepted_center = accepted_mol.add_atom(charged_atom(symbol, charge));
        for _ in 0..accepted {
            let hydrogen = accepted_mol.add_atom(element_atom("H"));
            accepted_mol
                .add_bond(accepted_center, hydrogen, BondOrder::Single)
                .expect("bond");
        }
        let accepted_report =
            valence_api::perceive_valence(&mut accepted_mol, ValenceModel::RdkitLike);
        assert!(
            accepted_report.is_ok(),
            "{symbol}{charge:+} valence {accepted}"
        );

        let mut rejected_mol = Molecule::new();
        let rejected_center = rejected_mol.add_atom(charged_atom(symbol, charge));
        for _ in 0..rejected {
            let hydrogen = rejected_mol.add_atom(element_atom("H"));
            rejected_mol
                .add_bond(rejected_center, hydrogen, BondOrder::Single)
                .expect("bond");
        }
        let rejected_report =
            valence_api::perceive_valence(&mut rejected_mol, ValenceModel::RdkitLike);
        assert!(
            matches!(
                rejected_report.issues.as_slice(),
                [ValenceIssue::ValenceExceeded { atom, .. }] if *atom == rejected_center
            ),
            "{symbol}{charge:+} valence {rejected}"
        );
    }
}

#[test]
fn valence_accepts_rdkit_phosphorus_minus_one_and_hydride_compatibility_cases() {
    let mut hexafluorophosphate = Molecule::new();
    let phosphorus = hexafluorophosphate.add_atom(charged_atom("P", -1));
    for _ in 0..6 {
        let fluorine =
            hexafluorophosphate.add_atom(Atom::new(Element::from_symbol("F").expect("fluorine")));
        hexafluorophosphate
            .add_bond(phosphorus, fluorine, BondOrder::Single)
            .expect("P-F bond");
    }
    assert!(
        valence_api::perceive_valence(&mut hexafluorophosphate, ValenceModel::RdkitLike).is_ok()
    );

    let mut bridged_hydride = Molecule::new();
    let hydrogen = bridged_hydride.add_atom(charged_atom("H", -1));
    let boron_a = bridged_hydride.add_atom(Atom::new(Element::from_symbol("B").expect("boron")));
    let boron_b = bridged_hydride.add_atom(Atom::new(Element::from_symbol("B").expect("boron")));
    bridged_hydride
        .add_bond(hydrogen, boron_a, BondOrder::Single)
        .expect("first hydride bond");
    bridged_hydride
        .add_bond(hydrogen, boron_b, BondOrder::Single)
        .expect("second hydride bond");
    assert!(valence_api::perceive_valence(&mut bridged_hydride, ValenceModel::RdkitLike).is_ok());
}

#[test]
fn valence_supports_simple_pubchem_main_group_ions_and_salts() {
    for (symbol, charge, expected_implicit_hydrogens) in [
        ("H", 1, 0),
        ("H", -1, 0),
        ("Rb", 1, 0),
        ("Cs", 1, 0),
        ("Be", 2, 0),
        ("Al", 3, 0),
        ("Ga", 3, 0),
        ("Tl", 1, 0),
        ("U", 2, 0),
        ("Pb", 2, 0),
        ("S", -2, 0),
        ("Se", -2, 0),
    ] {
        let mut mol = Molecule::new();
        let atom_id = mol.add_atom(charged_atom(symbol, charge));

        let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

        assert!(report.is_ok(), "{symbol}{charge:+} should be supported");
        assert_eq!(
            mol.atom(atom_id).expect("atom").implicit_hydrogens,
            Some(expected_implicit_hydrogens),
            "{symbol}{charge:+} implicit hydrogens"
        );
    }

    for symbol in ["Ac", "Cf"] {
        let mut mol = Molecule::new();
        let atom_id = mol.add_atom(element_atom(symbol));

        let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

        assert!(
            report.is_ok(),
            "isolated unsupported spectator {symbol} should be accepted"
        );
        assert_eq!(
            mol.atom(atom_id).expect("atom").implicit_hydrogens,
            Some(0),
            "{symbol} implicit hydrogens"
        );
    }

    let mut mercury_cyanide = Molecule::new();
    let mercury = mercury_cyanide.add_atom(charged_atom("Hg", -2));
    for _ in 0..4 {
        let carbon = mercury_cyanide.add_atom(element_atom("C"));
        let nitrogen = mercury_cyanide.add_atom(element_atom("N"));
        mercury_cyanide
            .add_bond(mercury, carbon, BondOrder::Single)
            .expect("mercury-carbon bond");
        mercury_cyanide
            .add_bond(carbon, nitrogen, BondOrder::Triple)
            .expect("cyanide bond");
    }
    let report = valence_api::perceive_valence(&mut mercury_cyanide, ValenceModel::RdkitLike);
    assert!(report.is_ok(), "tetracyanomercurate should be supported");
    assert_eq!(
        mercury_cyanide
            .atom(mercury)
            .expect("mercury")
            .implicit_hydrogens,
        Some(0)
    );

    let mut covalent_aluminum = Molecule::new();
    let aluminum = covalent_aluminum.add_atom(element_atom("Al"));
    for _ in 0..3 {
        let chlorine = covalent_aluminum.add_atom(element_atom("Cl"));
        covalent_aluminum
            .add_bond(aluminum, chlorine, BondOrder::Single)
            .expect("bond");
    }

    let report = valence_api::perceive_valence(&mut covalent_aluminum, ValenceModel::RdkitLike);

    assert!(
        report.is_ok(),
        "neutral trivalent aluminum should be supported"
    );
    assert_eq!(
        covalent_aluminum
            .atom(aluminum)
            .expect("aluminum")
            .implicit_hydrogens,
        Some(0)
    );

    let mut neutral_magnesium = Molecule::new();
    let magnesium = neutral_magnesium.add_atom(element_atom("Mg"));
    for _ in 0..2 {
        let chlorine = neutral_magnesium.add_atom(element_atom("Cl"));
        neutral_magnesium
            .add_bond(magnesium, chlorine, BondOrder::Single)
            .expect("bond");
    }

    let report = valence_api::perceive_valence(&mut neutral_magnesium, ValenceModel::RdkitLike);

    assert!(
        report.is_ok(),
        "neutral divalent magnesium should be supported"
    );
    assert_eq!(
        neutral_magnesium
            .atom(magnesium)
            .expect("magnesium")
            .implicit_hydrogens,
        Some(0)
    );
}

#[test]
fn molfile_wedge_assembles_tetrahedral_p_with_a_double_bond() {
    let input = r#"tetrahedral phosphorus
  molecules

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  2  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"#;
    let mut molecule = read_sdf_molecules(input)
        .expect("compact phosphorus regression parses")
        .into_iter()
        .next()
        .expect("one molecule");

    let report = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("tetracoordinate phosphorus should sanitize");
    let stereo = report.stereo.expect("stereo report");

    assert!(stereo.issues.is_empty());
    assert_eq!(stereo.assembled_elements.len(), 1);
    assert!(matches!(
        &stereo.assembled_elements[0].kind,
        StereoElementKind::Tetrahedral(stereo) if stereo.center == AtomId::new(0)
    ));
}

#[test]
fn molfile_wedge_assembles_pyramidal_s_with_a_lone_pair() {
    let input = r#"pyramidal sulfur
  molecules

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  2  0  0  0  0
  1  4  1  0  0  0  0
M  END
$$$$
"#;
    let mut molecule = read_sdf_molecules(input)
        .expect("compact sulfur regression parses")
        .into_iter()
        .next()
        .expect("one molecule");

    let report = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("pyramidal sulfur should sanitize");
    let stereo = report.stereo.expect("stereo report");

    assert!(stereo.issues.is_empty());
    assert_eq!(stereo.assembled_elements.len(), 1);
    assert!(matches!(
        &stereo.assembled_elements[0].kind,
        StereoElementKind::Tetrahedral(stereo)
            if stereo.center == AtomId::new(0)
                && stereo.carriers.contains(&StereoCarrier::ImplicitLonePair)
    ));
}
