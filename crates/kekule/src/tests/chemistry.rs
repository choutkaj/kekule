use super::*;

#[test]
fn interpretation_and_default_perception_are_separate() {
    let mut small = read_smiles("CCO").expect("smiles should parse");
    assert!(!small.as_molecule().perception().has_valence());

    small.perceive().expect("ethanol perception");

    assert!(small.as_molecule().perception().has_valence());
    assert!(small.as_molecule().perception().has_rings());
    assert!(small.as_molecule().perception().has_aromaticity());
    assert_eq!(
        small
            .as_molecule()
            .implicit_hydrogens(AtomId::new(2))
            .expect("oxygen"),
        Some(1)
    );
}

#[test]
fn default_perception_installs_one_ring_basis_for_aromaticity() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should parse");
    molecule.perceive().expect("benzene perception");

    let ring_set = molecule
        .as_molecule()
        .ring_set()
        .expect("installed ring basis");
    assert_eq!(ring_set.len(), 1);
    assert_eq!(
        molecule.as_molecule().perception().ring_basis_model(),
        Some(RingBasisModel::FiguerasSssrLike)
    );
    assert!(molecule.as_molecule().perception().has_aromaticity());
    assert_eq!(
        molecule
            .as_molecule()
            .atoms()
            .filter(|(atom, _)| molecule.as_molecule().atom_is_aromatic(*atom) == Ok(Some(true)))
            .count(),
        6
    );
}

#[test]
fn interpretation_owns_source_stereo_before_default_perception() {
    let (mut molecule, report) =
        read_smiles_with_report("C/C=C\\F").expect("directional smiles should interpret");

    assert_eq!(report.created_stereo_elements().len(), 1);
    assert_eq!(molecule.as_molecule().stereo_elements().count(), 1);

    molecule.perceive().expect("directional perception");

    assert_eq!(molecule.as_molecule().stereo_elements().count(), 1);
    assert!(!molecule.as_molecule().perception().has_stereo());
    assert!(molecule.as_molecule().perception().stereo_state().is_none());
}

#[test]
fn normalization_preserves_unknown_double_bond_stereo() {
    let mut molecule = read_smiles("CC=CC").expect("alkene should parse");
    let double_bond = molecule
        .as_molecule()
        .bonds()
        .find_map(|(bond_id, bond)| (bond.order == BondOrder::Double).then_some(bond_id))
        .expect("double bond");
    let source_stereo = [SourceStereoBondMark {
        bond: double_bond,
        from: molecule
            .as_molecule()
            .bond(double_bond)
            .expect("double bond")
            .a(),
        kind: SourceStereoBondMarkKind::DoubleBondEither,
    }];

    let report = molecule
        .canonicalize_fixture_with_source_stereo(&source_stereo)
        .expect("unknown double-bond stereo should normalize");

    assert_eq!(report.created_stereo_elements.len(), 1);
    let (_, element) = molecule
        .as_molecule()
        .stereo_elements()
        .next()
        .expect("unknown stereo element");
    assert!(matches!(
        &element.kind,
        StereoElementKind::DoubleBond(stereo)
            if stereo.bond == double_bond && stereo.orientation.is_none()
    ));
}

#[test]
fn default_perception_does_not_assign_coordinate_only_stereo() {
    let mut mol = Molecule::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left_carrier = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let right_carrier = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
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
    let mut molecule = SmallMolecule::from_molecule(mol);

    molecule
        .perceive()
        .expect("coordinate-only molecule should perceive discrete chemistry");

    assert_eq!(molecule.as_molecule().stereo_elements().count(), 0);

    let inferred = stereo_api::infer_coordinate_stereo(molecule.as_molecule())
        .expect("direct coordinate inference should succeed");
    assert_eq!(inferred.elements.len(), 1);
    assert_eq!(molecule.as_molecule().stereo_elements().count(), 0);

    let materialized = stereo_api::materialize_coordinate_stereo(molecule.as_molecule_mut())
        .expect("explicit coordinate materialization should succeed");
    assert_eq!(materialized.created_elements.len(), 1);
    assert_eq!(molecule.as_molecule().stereo_elements().count(), 1);
}

#[test]
fn failed_source_stereo_normalization_is_transactional() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let bond = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    let source_stereo = [SourceStereoBondMark {
        bond,
        from: a,
        kind: SourceStereoBondMarkKind::WedgeEither,
    }];
    let mut molecule = SmallMolecule::from_molecule(mol);
    let before = molecule.clone();

    let error = molecule
        .canonicalize_fixture_with_source_stereo(&source_stereo)
        .expect_err("unassembled stereo mark should fail normalization");

    assert!(matches!(
        error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnassembledTetrahedralBondMark {
            bond,
            kind: SourceStereoBondMarkKind::WedgeEither,
        })
    ));
    assert_eq!(molecule, before);
}

#[test]
fn normalization_treats_conflicting_wedges_as_nonfatal_ambiguity() {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let mut marked_bonds = Vec::new();
    for symbol in ["F", "Cl", "Br", "I"] {
        let carrier = mol
            .add_atom(element_atom(symbol))
            .expect("atom identifier capacity");
        marked_bonds.push(
            mol.add_bond(center, carrier, BondOrder::Single)
                .expect("carrier bond"),
        );
    }
    let source_stereo = marked_bonds
        .into_iter()
        .enumerate()
        .map(|(index, bond)| SourceStereoBondMark {
            bond,
            from: center,
            kind: if index % 2 == 0 {
                SourceStereoBondMarkKind::WedgeUp
            } else {
                SourceStereoBondMarkKind::WedgeDown
            },
        })
        .collect::<Vec<_>>();
    let mut molecule = SmallMolecule::from_molecule(mol);

    let report = molecule
        .canonicalize_fixture_with_source_stereo(&source_stereo)
        .expect("ambiguous drawing wedges should not reject valid chemistry");
    assert!(report
        .warnings
        .contains(&NormalizationWarning::AmbiguousTetrahedralWedgeMarks {
            center,
            mark_count: 4,
        }));
    assert_eq!(report.warnings.len(), 1);
    assert!(molecule.as_molecule().stereo_elements().next().is_none());
}

#[test]
fn failed_default_valence_perception_is_transactional() {
    let mut mol = Molecule::new();
    let carbon = mol.add_atom(carbon()).expect("atom identifier capacity");
    for _ in 0..5 {
        let hydrogen = mol
            .add_atom(Atom::new(Element::from_symbol("H").expect("hydrogen")))
            .expect("atom identifier capacity");
        mol.add_bond(carbon, hydrogen, BondOrder::Single)
            .expect("bond");
    }
    rings_api::perceive_ring_set(&mut mol).expect("ring perception should succeed");
    let mut molecule = SmallMolecule::from_molecule(mol);
    let before = molecule.clone();

    let error = molecule
        .perceive()
        .expect_err("pentavalent carbon should fail valence");

    assert!(matches!(
        error,
        perception_api::PerceptionError::Valence(ValenceError { issues })
            if matches!(issues.as_slice(), [ValenceIssue::ValenceExceeded { atom, .. }] if *atom == carbon)
    ));
    assert_eq!(molecule, before);
}

#[test]
fn invalid_aromatic_source_is_rejected_before_publication() {
    let error = read_smiles("c1cccc1")
        .expect_err("unmatchable aromatic representation must fail interpretation");

    assert!(error
        .to_string()
        .contains("invalid imported aromatic representation"));
}

#[test]
fn direct_aromaticity_perception_accepts_localized_aromatic_input() {
    let mut molecule = read_smiles("c1ccccc1").expect("aromatic representation localizes");

    aromaticity_api::perceive_aromaticity(molecule.as_molecule_mut(), AromaticityModel::RdkitLike)
        .expect("localized aromatic representation should be perceived");

    assert!(molecule
        .as_molecule()
        .bond_ids()
        .all(|bond| molecule.as_molecule().bond_is_aromatic(bond) == Ok(Some(true))));
}

#[test]
fn successful_default_perception_is_idempotent() {
    let mut molecule = read_smiles("CCO").expect("ethanol should parse");
    molecule
        .canonicalize_fixture()
        .expect("ethanol normalization");
    molecule
        .perceive()
        .expect("first perception should succeed");
    let once = molecule.clone();

    molecule
        .perceive()
        .expect("second perception should succeed");

    assert_eq!(molecule, once);
}

#[test]
fn normalization_cleanup_invalidates_preexisting_perception() {
    let mut mol = Molecule::new();
    let chlorine = mol
        .add_atom(Atom::new(Element::from_symbol("Cl").expect("chlorine")))
        .expect("atom identifier capacity");
    let oxo = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let hydroxyl = mol.add_atom(oxygen()).expect("atom identifier capacity");
    mol.add_bond(chlorine, oxo, BondOrder::Double)
        .expect("bond");
    mol.add_bond(chlorine, hydroxyl, BondOrder::Single)
        .expect("hydroxyl bond");
    mark_all_fresh(&mut mol);
    let mut molecule = SmallMolecule::from_molecule(mol);

    molecule
        .canonicalize_fixture()
        .expect("representation cleanup");

    assert_all_stale(molecule.as_molecule());
    assert_eq!(
        molecule
            .as_molecule()
            .atom(chlorine)
            .expect("chlorine")
            .formal_charge,
        1
    );
    assert_eq!(
        molecule
            .as_molecule()
            .atom(oxo)
            .expect("oxygen")
            .formal_charge,
        -1
    );
    assert_eq!(
        molecule
            .as_molecule()
            .atom(hydroxyl)
            .expect("hydroxyl oxygen")
            .formal_charge,
        0
    );
    assert_eq!(
        molecule
            .as_molecule()
            .bond(BondId::new(0))
            .expect("bond")
            .order,
        BondOrder::Single
    );
}

#[test]
fn valence_reports_excess_common_valence() {
    let mut mol = Molecule::new();
    let c = mol
        .add_atom(Atom::new(Element::from_symbol("C").expect("C")))
        .expect("atom identifier capacity");
    for _ in 0..5 {
        let h = mol
            .add_atom(Atom::new(Element::from_symbol("H").expect("H")))
            .expect("atom identifier capacity");
        mol.add_bond(c, h, BondOrder::Single).expect("bond");
    }

    let error = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike)
        .expect_err("pentavalent carbon should fail strict perception");

    assert_eq!(error.issues.len(), 1);

    valence_api::perceive_valence_with_options(
        &mut mol,
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    )
    .expect("permissive valence inspection should succeed");
    assert_eq!(mol.implicit_hydrogens(c).expect("carbon"), Some(0));
}

#[test]
fn failed_strict_valence_perception_preserves_complete_previous_perception_state() {
    let mut mol = Molecule::new();
    let carbon = mol
        .add_atom(element_atom("C"))
        .expect("atom identifier capacity");
    for _ in 0..5 {
        let hydrogen = mol
            .add_atom(element_atom("H"))
            .expect("atom identifier capacity");
        mol.add_bond(carbon, hydrogen, BondOrder::Single)
            .expect("bond");
    }
    let previous = PerceptionState::builder()
        .with_valence(
            Some(ValenceModel::RdkitLike),
            mol.atom_ids()
                .map(|atom| (atom, if atom == carbon { 2 } else { 0 }))
                .collect(),
        )
        .expect("previous valence")
        .with_rings(
            RingMembership::from_slot_flags(
                vec![false; mol.atom_count()],
                vec![false; mol.bond_count()],
            ),
            None,
        )
        .with_aromaticity(AromaticityModel::RdkitLike, Vec::new(), Vec::new())
        .expect("previous aromaticity")
        .build();
    mol.install_perception_state(previous.clone())
        .expect("previous perception");

    let error = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike)
        .expect_err("pentavalent carbon should fail strict perception");

    assert!(matches!(
        error.issues.as_slice(),
        [ValenceIssue::ValenceExceeded { atom, .. }] if *atom == carbon
    ));
    assert_eq!(mol.perception(), &previous);
    assert_eq!(mol.implicit_hydrogens(carbon).unwrap(), Some(2));
}

#[test]
fn unsupported_valence_target_remains_strictly_diagnostic_and_permissively_installable() {
    let mut strict = Molecule::new();
    let carbon = strict
        .add_atom(charged_atom("C", 7))
        .expect("atom identifier capacity");

    let error = valence_api::perceive_valence(&mut strict, ValenceModel::RdkitLike)
        .expect_err("out-of-range charge adjustment should be unsupported");

    assert_eq!(error.issues, vec![ValenceIssue::UnsupportedElement(carbon)]);
    assert_eq!(strict.perception(), &PerceptionState::default());
    assert_eq!(strict.implicit_hydrogens(carbon).unwrap(), None);

    valence_api::perceive_valence_with_options(
        &mut strict,
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    )
    .expect("permissive unsupported-element inspection should install");
    assert!(strict.perception().has_valence());
    assert_eq!(strict.implicit_hydrogens(carbon).unwrap(), Some(0));
}

#[test]
fn valence_counts_high_degree_atoms_without_narrowing_or_panicking() {
    let mut mol = Molecule::new();
    let carbon = mol
        .add_atom(element_atom("C"))
        .expect("atom identifier capacity");
    for _ in 0..300 {
        let hydrogen = mol
            .add_atom(element_atom("H"))
            .expect("atom identifier capacity");
        mol.add_bond(carbon, hydrogen, BondOrder::Single)
            .expect("bond");
    }

    let error = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike)
        .expect_err("high-degree carbon should fail strict perception");

    assert_eq!(
        error.issues,
        vec![ValenceIssue::ValenceExceeded {
            atom: carbon,
            explicit_valence: 300,
            max_allowed: 4,
        }]
    );
    assert_eq!(mol.implicit_hydrogens(carbon).expect("carbon"), None);
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
        let atom_id = mol
            .add_atom(element_atom(symbol))
            .expect("atom identifier capacity");

        let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

        assert!(report.is_ok(), "neutral {symbol} should be supported");
        assert_eq!(
            mol.implicit_hydrogens(atom_id).expect("atom"),
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
        let accepted_center = accepted_mol
            .add_atom(charged_atom(symbol, charge))
            .expect("atom identifier capacity");
        for _ in 0..accepted {
            let hydrogen = accepted_mol
                .add_atom(element_atom("H"))
                .expect("atom identifier capacity");
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
        let rejected_center = rejected_mol
            .add_atom(charged_atom(symbol, charge))
            .expect("atom identifier capacity");
        for _ in 0..rejected {
            let hydrogen = rejected_mol
                .add_atom(element_atom("H"))
                .expect("atom identifier capacity");
            rejected_mol
                .add_bond(rejected_center, hydrogen, BondOrder::Single)
                .expect("bond");
        }
        let rejected_error =
            valence_api::perceive_valence(&mut rejected_mol, ValenceModel::RdkitLike)
                .expect_err("excess hypervalent anion valence should fail");
        assert!(
            matches!(
                rejected_error.issues.as_slice(),
                [ValenceIssue::ValenceExceeded { atom, .. }] if *atom == rejected_center
            ),
            "{symbol}{charge:+} valence {rejected}"
        );
    }
}

#[test]
fn valence_accepts_rdkit_phosphorus_minus_one_and_hydride_compatibility_cases() {
    let mut hexafluorophosphate = Molecule::new();
    let phosphorus = hexafluorophosphate
        .add_atom(charged_atom("P", -1))
        .expect("atom identifier capacity");
    for _ in 0..6 {
        let fluorine = hexafluorophosphate
            .add_atom(Atom::new(Element::from_symbol("F").expect("fluorine")))
            .expect("atom identifier capacity");
        hexafluorophosphate
            .add_bond(phosphorus, fluorine, BondOrder::Single)
            .expect("P-F bond");
    }
    assert!(
        valence_api::perceive_valence(&mut hexafluorophosphate, ValenceModel::RdkitLike).is_ok()
    );

    let mut bridged_hydride = Molecule::new();
    let hydrogen = bridged_hydride
        .add_atom(charged_atom("H", -1))
        .expect("atom identifier capacity");
    let boron_a = bridged_hydride
        .add_atom(Atom::new(Element::from_symbol("B").expect("boron")))
        .expect("atom identifier capacity");
    let boron_b = bridged_hydride
        .add_atom(Atom::new(Element::from_symbol("B").expect("boron")))
        .expect("atom identifier capacity");
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
        let atom_id = mol
            .add_atom(charged_atom(symbol, charge))
            .expect("atom identifier capacity");

        let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

        assert!(report.is_ok(), "{symbol}{charge:+} should be supported");
        assert_eq!(
            mol.implicit_hydrogens(atom_id).expect("atom"),
            Some(expected_implicit_hydrogens),
            "{symbol}{charge:+} implicit hydrogens"
        );
    }

    for symbol in ["Ac", "Cf"] {
        let mut mol = Molecule::new();
        let atom_id = mol
            .add_atom(element_atom(symbol))
            .expect("atom identifier capacity");

        let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

        assert!(
            report.is_ok(),
            "isolated unsupported spectator {symbol} should be accepted"
        );
        assert_eq!(
            mol.implicit_hydrogens(atom_id).expect("atom"),
            Some(0),
            "{symbol} implicit hydrogens"
        );
    }

    let mut mercury_cyanide = Molecule::new();
    let mercury = mercury_cyanide
        .add_atom(charged_atom("Hg", -2))
        .expect("atom identifier capacity");
    for _ in 0..4 {
        let carbon = mercury_cyanide
            .add_atom(element_atom("C"))
            .expect("atom identifier capacity");
        let nitrogen = mercury_cyanide
            .add_atom(element_atom("N"))
            .expect("atom identifier capacity");
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
            .implicit_hydrogens(mercury)
            .expect("mercury"),
        Some(0)
    );

    let mut covalent_aluminum = Molecule::new();
    let aluminum = covalent_aluminum
        .add_atom(element_atom("Al"))
        .expect("atom identifier capacity");
    for _ in 0..3 {
        let chlorine = covalent_aluminum
            .add_atom(element_atom("Cl"))
            .expect("atom identifier capacity");
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
            .implicit_hydrogens(aluminum)
            .expect("aluminum"),
        Some(0)
    );

    let mut neutral_magnesium = Molecule::new();
    let magnesium = neutral_magnesium
        .add_atom(element_atom("Mg"))
        .expect("atom identifier capacity");
    for _ in 0..2 {
        let chlorine = neutral_magnesium
            .add_atom(element_atom("Cl"))
            .expect("atom identifier capacity");
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
            .implicit_hydrogens(magnesium)
            .expect("magnesium"),
        Some(0)
    );
}

#[test]
fn molfile_wedge_assembles_tetrahedral_p_with_a_double_bond() {
    let input = r#"tetrahedral phosphorus
  kekule

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
    let molecule = read_sdf_molecules(input)
        .expect("compact phosphorus regression parses")
        .into_iter()
        .next()
        .expect("one molecule");

    assert_eq!(molecule.as_molecule().stereo_elements().count(), 1);
    let element = molecule
        .as_molecule()
        .stereo_elements()
        .next()
        .expect("created tetrahedral element")
        .1;
    assert!(matches!(
        &element.kind,
        StereoElementKind::Tetrahedral(stereo) if stereo.center == AtomId::new(0)
    ));
}

#[test]
fn molfile_wedge_assembles_pyramidal_s_with_a_lone_pair() {
    let input = r#"pyramidal sulfur
  kekule

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
    let molecule = read_sdf_molecules(input)
        .expect("compact sulfur regression parses")
        .into_iter()
        .next()
        .expect("one molecule");

    assert_eq!(molecule.as_molecule().stereo_elements().count(), 1);
    let element = molecule
        .as_molecule()
        .stereo_elements()
        .next()
        .expect("created tetrahedral element")
        .1;
    assert!(matches!(
        &element.kind,
        StereoElementKind::Tetrahedral(stereo)
            if stereo.center == AtomId::new(0)
                && stereo.carriers.contains(&StereoCarrier::ImplicitLonePair)
    ));
}
