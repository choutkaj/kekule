use super::*;

#[test]
fn smiles_parses_branches_rings_brackets_and_fragments_without_sanitizing() {
    let small = read_smiles_str(
        "C(C)O.C1=CC=CC=C1.[13NH4+:7].[C@@H](N)O",
        SmilesParseOptions,
    )
    .expect("smiles should parse");

    assert_eq!(small.mol.atom_count(), 13);
    assert_eq!(small.mol.bond_count(), 10);
    assert_eq!(small.mol.perception().valence, ComputedState::Absent);
    let bracket_atom = small.mol.atom(AtomId::new(9)).expect("bracket atom");
    assert_eq!(bracket_atom.isotope, Some(13));
    assert_eq!(bracket_atom.explicit_hydrogens, 4);
    assert!(bracket_atom.no_implicit_hydrogens);
    assert_eq!(bracket_atom.formal_charge, 1);
    assert_eq!(bracket_atom.atom_map, Some(7));
    let chiral_atom = small
        .mol
        .atom(AtomId::new(10))
        .expect("chiral bracket atom");
    assert_eq!(
        chiral_atom.chiral,
        Some(AtomStereo::TetrahedralCounterClockwise)
    );
    assert_eq!(chiral_atom.explicit_hydrogens, 1);
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
        "C/C=C\\C",
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
        let parsed = std::panic::catch_unwind(|| read_smiles_str(input, SmilesParseOptions))
            .unwrap_or_else(|_| panic!("`{input}` panicked"));
        let error = parsed.expect_err("malformed SMILES should fail");
        assert!(error.offset <= input.len(), "offset for `{input}`");
        assert!(!error.message.is_empty(), "message for `{input}`");
    }
}

#[test]
fn smiles_writer_round_trips_graph_shape() {
    let small = read_smiles_str("CC(=O)O", SmilesParseOptions).expect("smiles should parse");
    let text = write_smiles(&small, SmilesWriteOptions).expect("smiles should write");
    let reparsed = read_smiles_str(&text, SmilesParseOptions).expect("written smiles should parse");

    assert_eq!(reparsed.mol.atom_count(), small.mol.atom_count());
    assert_eq!(reparsed.mol.bond_count(), small.mol.bond_count());
}

#[test]
fn canonical_smiles_is_stable_across_atom_order_for_tree_roles() {
    let mut first = SmallMolecule {
        mol: Molecule::new(),
    };
    let first_terminal_a = first.mol.add_atom(carbon());
    let first_center = first.mol.add_atom(carbon());
    let first_terminal_b = first.mol.add_atom(carbon());
    first
        .mol
        .add_bond(first_terminal_a, first_center, BondOrder::Single)
        .expect("bond should be valid");
    first
        .mol
        .add_bond(first_center, first_terminal_b, BondOrder::Single)
        .expect("bond should be valid");
    sanitize_small_molecule(&mut first, SanitizeOptions::default()).expect("propane sanitizes");

    let mut second = SmallMolecule {
        mol: Molecule::new(),
    };
    let second_center = second.mol.add_atom(carbon());
    let second_terminal_a = second.mol.add_atom(carbon());
    let second_terminal_b = second.mol.add_atom(carbon());
    second
        .mol
        .add_bond(second_center, second_terminal_a, BondOrder::Single)
        .expect("bond should be valid");
    second
        .mol
        .add_bond(second_center, second_terminal_b, BondOrder::Single)
        .expect("bond should be valid");
    sanitize_small_molecule(&mut second, SanitizeOptions::default()).expect("propane sanitizes");

    let first_written = write_canonical_smiles(&first, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let second_written = write_canonical_smiles(&second, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");

    assert_eq!(first_written, second_written);
    assert_eq!(first_written, "CCC");
    read_smiles_str(&first_written, SmilesParseOptions).expect("canonical output should parse");
}

#[test]
fn canonical_smiles_sorts_disconnected_components() {
    let mut first = read_smiles_str("O.C", SmilesParseOptions).expect("SMILES parses");
    let mut second = read_smiles_str("C.O", SmilesParseOptions).expect("SMILES parses");
    sanitize_small_molecule(&mut first, SanitizeOptions::default()).expect("first sanitizes");
    sanitize_small_molecule(&mut second, SanitizeOptions::default()).expect("second sanitizes");

    assert_eq!(
        write_canonical_smiles(&first, CanonicalSmilesWriteOptions)
            .expect("canonical SMILES should write"),
        write_canonical_smiles(&second, CanonicalSmilesWriteOptions)
            .expect("canonical SMILES should write")
    );
}

#[test]
fn canonical_smiles_ignores_stereo_for_non_isomeric_output() {
    let mut molecule =
        read_smiles_str("N[C@H](O)C", SmilesParseOptions).expect("chiral SMILES parses");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("chiral molecule sanitizes");

    assert!(write_smiles(&molecule, SmilesWriteOptions)
        .expect_err("ordinary writer should reject lossy atom stereo")
        .message
        .contains("atom stereochemistry"));

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("non-isomeric canonical SMILES should ignore atom stereo");
    let reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");

    assert!(!written.contains('['), "{written}");
    assert!(reparsed.mol.atoms().all(|(_, atom)| atom.chiral.is_none()));
    assert_eq!(reparsed.mol.atom_count(), molecule.mol.atom_count());
    assert_eq!(reparsed.mol.bond_count(), molecule.mol.bond_count());

    let isotope = read_smiles_str("[11CH3]OC", SmilesParseOptions).expect("isotope parses");
    assert_eq!(
        write_canonical_smiles(&isotope, CanonicalSmilesWriteOptions)
            .expect("non-isomeric canonical SMILES should ignore isotope labels"),
        "COC"
    );
}

#[test]
fn canonical_smiles_round_trips_supported_branch_and_ring_graphs() {
    for input in ["CC(=O)O", "C1CCCCC1", "c1ccccc1"] {
        let mut molecule = read_smiles_str(input, SmilesParseOptions)
            .unwrap_or_else(|_| panic!("SMILES should parse: {input}"));
        sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
            .unwrap_or_else(|_| panic!("SMILES should sanitize: {input}"));
        let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
            .unwrap_or_else(|_| panic!("canonical SMILES should write: {input}"));
        let reparsed = read_smiles_str(&written, SmilesParseOptions)
            .unwrap_or_else(|_| panic!("canonical output should parse: {written}"));

        assert_eq!(reparsed.mol.atom_count(), molecule.mol.atom_count());
        assert_eq!(reparsed.mol.bond_count(), molecule.mol.bond_count());
    }
}

#[test]
fn canonical_smiles_prefers_clean_simple_ring_closure() {
    let molecule = read_smiles_str("C1=CC=CC=C1", SmilesParseOptions).expect("benzene parses");

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");

    assert_eq!(written, "C1=CC=CC=C1");
}

#[test]
fn aromatic_smiles_omitted_bonds_sanitize_with_expected_hydrogens() {
    let mut benzene =
        read_smiles_str("c1ccccc1", SmilesParseOptions).expect("benzene should parse");
    assert!(benzene
        .mol
        .bonds()
        .all(|(_, bond)| bond.order == BondOrder::Aromatic));
    sanitize_small_molecule(&mut benzene, SanitizeOptions::default())
        .expect("benzene should sanitize");
    for (_, atom) in benzene.mol.atoms() {
        assert_eq!(atom.implicit_hydrogens, Some(1));
        assert!(atom.aromatic);
    }

    let mut pyridine =
        read_smiles_str("n1ccccc1", SmilesParseOptions).expect("pyridine should parse");
    sanitize_small_molecule(&mut pyridine, SanitizeOptions::default())
        .expect("pyridine should sanitize");
    let nitrogen = pyridine.mol.atom(AtomId::new(0)).expect("nitrogen");
    assert_eq!(nitrogen.implicit_hydrogens, Some(0));
    for atom_id in 1..6 {
        let atom = pyridine.mol.atom(AtomId::new(atom_id)).expect("carbon");
        assert_eq!(atom.implicit_hydrogens, Some(1));
    }

    let mut pyridinium =
        read_smiles_str("[n+]1ccccc1", SmilesParseOptions).expect("pyridinium should parse");
    sanitize_small_molecule(&mut pyridinium, SanitizeOptions::default())
        .expect("pyridinium should sanitize");
    let nitrogen = pyridinium.mol.atom(AtomId::new(0)).expect("nitrogen");
    assert!(nitrogen.aromatic);
    assert_eq!(nitrogen.formal_charge, 1);
    assert_eq!(nitrogen.implicit_hydrogens, Some(0));

    for smiles in [
        "[nH]1cccc1",
        "c1ccoc1",
        "c1ccsc1",
        "c1ccc2ccccc2c1",
        "Cc1ccccc1",
        "c1ccccc1.CC",
        "C%10CCCCC%10",
    ] {
        let mut molecule = read_smiles_str(smiles, SmilesParseOptions)
            .unwrap_or_else(|_| panic!("supported aromatic SMILES should parse: {smiles}"));
        sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
            .unwrap_or_else(|_| panic!("supported aromatic SMILES should sanitize: {smiles}"));
        let written = write_smiles(&molecule, SmilesWriteOptions)
            .unwrap_or_else(|_| panic!("supported aromatic SMILES should write: {smiles}"));
        read_smiles_str(&written, SmilesParseOptions)
            .unwrap_or_else(|_| panic!("writer output should parse: {written}"));
    }
}

#[test]
fn invalid_lowercase_aromatic_ring_returns_structured_error() {
    for smiles in ["c1cccc1", "c1ccccc1.c1cccc1"] {
        let mut molecule =
            read_smiles_str(smiles, SmilesParseOptions).expect("raw syntax should parse");
        let error = sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
            .expect_err("invalid aromatic ring should fail sanitization");
        assert!(matches!(
            error,
            SanitizeError::Aromaticity(AromaticityError::InvalidAromaticRepresentation(_))
        ));
    }
}

#[test]
fn thiocarbonyl_chalcogen_ring_sanitizes_aromatic_like_rdkit() {
    let mut molecule =
        read_smiles_str("CCN(CC)C1=NC(=S)N(C(=S)S1)C(=S)N(CC)CC", SmilesParseOptions)
            .expect("thiocarbonyl heterocycle should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("thiocarbonyl heterocycle should sanitize");

    let aromatic_atoms = molecule
        .mol
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    let aromatic_bonds = molecule
        .mol
        .bonds()
        .filter(|(_, bond)| bond.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 6);
    assert_eq!(aromatic_bonds, 6);

    let written = write_smiles(&molecule, SmilesWriteOptions)
        .expect("sanitized thiocarbonyl heterocycle should write");
    let reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("writer output should parse");
    assert_eq!(reparsed.mol.atom_count(), molecule.mol.atom_count());
    assert_eq!(reparsed.mol.bond_count(), molecule.mol.bond_count());
}

#[test]
fn fused_chalcogen_bridge_does_not_over_aromatize_hetero_bridge() {
    let mut molecule = read_smiles_str("CSC1=CC2=C(C=C1)SC3=CC=CC=C3N2", SmilesParseOptions)
        .expect("phenothiazine-like heterocycle should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("phenothiazine-like heterocycle should sanitize");

    let sulfur_bridge = molecule.mol.atom(AtomId::new(8)).expect("bridge sulfur");
    let nitrogen_bridge = molecule.mol.atom(AtomId::new(15)).expect("bridge nitrogen");
    assert!(!sulfur_bridge.aromatic);
    assert!(!nitrogen_bridge.aromatic);

    let aromatic_atoms = molecule
        .mol
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 12);
}

#[test]
fn bracket_carbon_suppresses_implicit_hydrogens() {
    let mut molecule = read_smiles_str("C1=CC=C2C(=C1)[CH]C3=CC=CC=C32", SmilesParseOptions)
        .expect("bracket carbon fused aromatic should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("bracket carbon fused aromatic should sanitize");

    let bracket_carbon = molecule.mol.atom(AtomId::new(6)).expect("bracket carbon");
    assert!(bracket_carbon.no_implicit_hydrogens);
    assert_eq!(bracket_carbon.explicit_hydrogens, 1);
    assert_eq!(bracket_carbon.implicit_hydrogens, Some(0));
}

#[test]
fn fused_polycycle_aromatic_core_extends_to_fused_edge() {
    let mut molecule = read_smiles_str(
        "C1CCC2=C(C1)C3=C(C=CC4=C3C5=C(C=C4)C=CC(=C25)[N+](=O)[O-])[N+](=O)[O-]",
        SmilesParseOptions,
    )
    .expect("fused polycycle should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused polycycle should sanitize");

    for atom_id in 0..3 {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("saturated atom")
                .aromatic,
            "saturated ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("fused core atom")
                .aromatic,
            "fused core atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_aromatic_component_preserves_explicit_single_bond() {
    let mut molecule = read_smiles_str(
        "[H]c1c([H])c([H])c2c3c([H])c([H])n(C([H])([H])[H])c(C([H])([H])[H])c-3nc2c1[H]",
        SmilesParseOptions,
    )
    .expect("fused aromatic system should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused aromatic system should sanitize");

    let explicit_single_between_aromatic_atoms = molecule
        .mol
        .bonds()
        .filter(|(_, bond)| {
            matches!(bond.order, BondOrder::Single)
                && molecule.mol.atom(bond.a()).is_ok_and(|atom| atom.aromatic)
                && molecule.mol.atom(bond.b()).is_ok_and(|atom| atom.aromatic)
        })
        .collect::<Vec<_>>();
    assert_eq!(explicit_single_between_aromatic_atoms.len(), 1);
    assert!(!explicit_single_between_aromatic_atoms[0].1.aromatic);

    let written =
        write_smiles(&molecule, SmilesWriteOptions).expect("fused aromatic system should write");
    assert!(written.contains('-'));
    let reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("writer output should parse");
    assert_eq!(
        reparsed
            .mol
            .bonds()
            .filter(|(_, bond)| {
                matches!(bond.order, BondOrder::Single)
                    && reparsed.mol.atom(bond.a()).is_ok_and(|atom| atom.aromatic)
                    && reparsed.mol.atom(bond.b()).is_ok_and(|atom| atom.aromatic)
            })
            .count(),
        1
    );
}

#[test]
fn fused_quinone_cn_core_excludes_carbonyl_centers() {
    let mut molecule = read_smiles_str(
        "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C4=C(C=C3)C(=O)C5=CC=CC=C5N4",
        SmilesParseOptions,
    )
    .expect("fused quinone should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused quinone should sanitize");

    for atom_id in [6, 10] {
        assert!(
            !molecule
                .mol
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
                .mol
                .atom(AtomId::new(atom_id))
                .expect("fused aromatic atom")
                .aromatic,
            "fused aromatic atom {atom_id} should be aromatic"
        );
    }
    assert_eq!(
        molecule
            .mol
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        20
    );
}

#[test]
fn indole_quinone_keeps_carbonyl_ring_atoms_aliphatic() {
    let mut molecule = read_smiles_str("C1=CC=C(C=C1)C2=CC3=C(N2)C(=O)C=CC3=O", SmilesParseOptions)
        .expect("indole quinone should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("indole quinone should sanitize");

    for atom_id in [11, 13, 14, 15] {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("quinone ring atom")
                .aromatic,
            "quinone ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("fused aromatic atom")
                .aromatic,
            "fused aromatic atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_imine_and_pyrimidinedione_aromaticity_matches_reference_shape() {
    let mut molecule =
        read_smiles_str("CC1=NC2=CC=CC=C2C1=CC3=C(NC(=O)NC3=O)O", SmilesParseOptions)
            .expect("fused imine and pyrimidinedione should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused imine and pyrimidinedione should sanitize");

    for atom_id in [1, 2, 9] {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("imine-ring atom")
                .aromatic,
            "imine-ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 16, 17] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("aromatic fused atom")
                .aromatic,
            "fused aromatic atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn exocyclic_iminium_sulfur_ring_remains_aromatic() {
    let mut molecule = read_smiles_str("CN(C1=NC(=[N+](C)C)SS1)C(=S)SC", SmilesParseOptions)
        .expect("exocyclic iminium sulfur ring should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("exocyclic iminium sulfur ring should sanitize");

    for atom_id in [2, 3, 4, 8, 9] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("heteroaromatic ring atom")
                .aromatic,
            "heteroaromatic ring atom {atom_id} should be aromatic"
        );
    }
    let exocyclic_iminium = molecule.mol.atom(AtomId::new(5)).expect("iminium N");
    assert!(!exocyclic_iminium.aromatic);
    assert_eq!(exocyclic_iminium.formal_charge, 1);
}

#[test]
fn fused_exocyclic_imine_sulfur_ring_remains_aromatic() {
    let mut molecule = read_smiles_str(
        "CCCCCCCCCCCCCCCCS(=O)(=O)N(C(=O)OCC)N=C1N(C2=CC=CC=C2S1)C",
        SmilesParseOptions,
    )
    .expect("fused exocyclic imine sulfur ring should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused exocyclic imine sulfur ring should sanitize");

    for atom_id in [26, 27, 28, 29, 30, 31, 32, 33, 34] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("fused benzothiazine atom")
                .aromatic,
            "fused benzothiazine atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_sulfoxide_ring_does_not_follow_benzene_aromaticity() {
    let mut molecule = read_smiles_str("CCCCCN1SC2=CC=CC=C2S1=O", SmilesParseOptions)
        .expect("fused sulfoxide ring should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused sulfoxide ring should sanitize");

    for atom_id in [5, 6, 13] {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("sulfoxide ring atom")
                .aromatic,
            "sulfoxide ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [7, 8, 9, 10, 11, 12] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("benzene ring atom")
                .aromatic,
            "benzene ring atom {atom_id} should stay aromatic"
        );
    }
}

#[test]
fn neutral_exocyclic_alkene_sulfur_ring_stays_aliphatic() {
    let mut molecule = read_smiles_str(
        "C1=CC=C2C(=C1)C=CC3=C2[N+](=C(S3)C=C4N(C5=CC=CC=C5S4)CCCS(=O)(=O)O)CCCS(=O)(=O)O",
        SmilesParseOptions,
    )
    .expect("mixed sulfur fused system should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("mixed sulfur fused system should sanitize");

    for atom_id in [14, 15, 22] {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("neutral sulfur ring atom")
                .aromatic,
            "neutral sulfur ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [10, 11, 12] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("cationic sulfur ring atom")
                .aromatic,
            "cationic sulfur ring atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_seven_membered_ether_ring_stays_aliphatic() {
    let mut molecule = read_smiles_str(
            "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C6(C4)C(=O)C7=C8N6CCC9=C8C(=C(C=C9)OC)OC1=C7C=CC(=C1O)OC",
            SmilesParseOptions,
        )
        .expect("fused ether polycycle should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused ether polycycle should sanitize");

    assert!(molecule
        .mol
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "O")
        .all(|(_, atom)| !atom.aromatic));
}

#[test]
fn charged_bracket_halogen_and_bismuth_salt_sanitizes() {
    let mut molecule = read_smiles_str(
        "C1CC2CCC[N-]C2C(C1)[OH2+].C1C=CC2=CC=CC(C2=N1)[OH2+].[ClH2+].Cl.[Bi+3]",
        SmilesParseOptions,
    )
    .expect("charged bracket salt should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("charged bracket salt should sanitize");

    let protonated_chlorine = molecule.mol.atom(AtomId::new(22)).expect("chlorine");
    assert_eq!(protonated_chlorine.element.symbol(), "Cl");
    assert_eq!(protonated_chlorine.formal_charge, 1);
    assert_eq!(protonated_chlorine.explicit_hydrogens, 2);
    assert_eq!(protonated_chlorine.implicit_hydrogens, Some(0));

    let bismuth = molecule.mol.atom(AtomId::new(24)).expect("bismuth");
    assert_eq!(bismuth.element.symbol(), "Bi");
    assert_eq!(bismuth.formal_charge, 3);
    assert_eq!(bismuth.implicit_hydrogens, Some(0));
}

#[test]
fn oxide_dianion_transition_metal_salt_sanitizes() {
    let mut molecule = read_smiles_str("[O-2].[O-2].[O-2].[Cr+3].[Fe+3]", SmilesParseOptions)
        .expect("oxide transition-metal salt should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("oxide transition-metal salt should sanitize");

    for atom_id in [0, 1, 2] {
        let oxygen = molecule.mol.atom(AtomId::new(atom_id)).expect("oxide");
        assert_eq!(oxygen.element.symbol(), "O");
        assert_eq!(oxygen.formal_charge, -2);
        assert_eq!(oxygen.implicit_hydrogens, Some(0));
    }
}

#[test]
fn hydroxide_niobium_v_salt_sanitizes() {
    let mut molecule = read_smiles_str("[OH-].[Nb+5]", SmilesParseOptions)
        .expect("niobium hydroxide salt should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("niobium hydroxide salt should sanitize");

    let niobium = molecule.mol.atom(AtomId::new(1)).expect("niobium");
    assert_eq!(niobium.element.symbol(), "Nb");
    assert_eq!(niobium.formal_charge, 5);
    assert_eq!(niobium.implicit_hydrogens, Some(0));
}

#[test]
fn formate_indium_salt_sanitizes() {
    let mut molecule = read_smiles_str("C(=O)[O-].C(=O)[O-].C(=O)[O-].[In+3]", SmilesParseOptions)
        .expect("indium formate salt should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("indium formate salt should sanitize");

    let indium = molecule.mol.atom(AtomId::new(9)).expect("indium");
    assert_eq!(indium.element.symbol(), "In");
    assert_eq!(indium.formal_charge, 3);
    assert_eq!(indium.implicit_hydrogens, Some(0));
}

#[test]
fn periodate_cleanup_sanitizes_iodine_plus_three() {
    let mut molecule =
        read_smiles_str("[O-]I(=O)(=O)=O", SmilesParseOptions).expect("periodate should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("periodate should sanitize");

    let iodine = molecule.mol.atom(AtomId::new(1)).expect("iodine");
    assert_eq!(iodine.element.symbol(), "I");
    assert_eq!(iodine.formal_charge, 3);
    assert_eq!(iodine.implicit_hydrogens, Some(0));
}

#[test]
fn uranyl_beta_diketonate_salt_sanitizes() {
    let mut molecule = read_smiles_str(
        "C1=CC=C(C=C1)C(=O)[CH-]C(=O)C2=CC=CC=C2.C1=CC=C(C=C1)C(=O)[CH-]C(=O)C2=CC=CC=C2.O=[U+2]=O",
        SmilesParseOptions,
    )
    .expect("uranyl salt should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("uranyl salt should sanitize");

    let uranium = molecule.mol.atom(AtomId::new(35)).expect("uranium");
    assert_eq!(uranium.element.symbol(), "U");
    assert_eq!(uranium.formal_charge, 2);
    assert_eq!(uranium.implicit_hydrogens, Some(0));
}

#[test]
fn cyclopentadienyl_anion_sanitizes_aromatic() {
    let mut molecule = read_smiles_str(
        "C1CCOC1.[CH-]1[C-]=[C-][C-]=[C-]1.Cl[Cr]Cl",
        SmilesParseOptions,
    )
    .expect("cyclopentadienyl chromium salt should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("cyclopentadienyl chromium salt should sanitize");

    for atom_id in 5..=9 {
        let atom = molecule
            .mol
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
    let mut molecule = read_smiles_str(
        "CC(C)(C)NN=C(C1C=CCS1(=O)=O)C(=O)NC2=C(C(=O)C3=CC=CC=C3C2=O)Cl",
        SmilesParseOptions,
    )
    .expect("fused quinone sulfone should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused quinone sulfone should sanitize");

    for atom_id in [17, 18, 19, 27] {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("quinone ring atom")
                .aromatic,
            "quinone ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [21, 22, 23, 24, 25, 26] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("benzene ring atom")
                .aromatic,
            "benzene ring atom {atom_id} should stay aromatic"
        );
    }
}

#[test]
fn singly_carbonylated_fused_ring_stays_aromatic() {
    let mut molecule = read_smiles_str(
        "CNCCN=C1C=CC2=C3C1=C(C4=C(C=CC(=O)C4=C3NN2CCNCCO)O)O.O.Cl.Cl",
        SmilesParseOptions,
    )
    .expect("singly carbonylated fused ring should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("singly carbonylated fused ring should sanitize");

    for atom_id in [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("aromatic fused atom")
                .aromatic,
            "fused atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn saturated_fused_ring_does_not_follow_aromatic_core() {
    let mut molecule = read_smiles_str(
        "C1CCC2=NC3=CC=CC=C3C(=C2C1)[NH2+]CCSCCCl.[Cl-]",
        SmilesParseOptions,
    )
    .expect("saturated fused ring salt should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("saturated fused ring salt should sanitize");

    for atom_id in [0, 1, 2, 13] {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("saturated fused atom")
                .aromatic,
            "saturated fused atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("aromatic core atom")
                .aromatic,
            "aromatic core atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn canonical_saturated_fused_ring_round_trip_stays_aliphatic() {
    let mut molecule = read_smiles_str(
        "C1CCC2=NC3=CC=CC=C3C(=C2C1)[NH2+]CCSCCCl.[Cl-]",
        SmilesParseOptions,
    )
    .expect("saturated fused ring salt should parse");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("saturated fused ring salt should sanitize");

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("saturated fused ring salt should canonicalize");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let saturated_carbons = reparsed
        .mol
        .atoms()
        .filter(|(id, atom)| {
            atom.element.symbol() == "C"
                && atom.implicit_hydrogens == Some(2)
                && reparsed
                    .mol
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
    let mut molecule = read_smiles_str("C1C(C(=O)C2=CC=CC=C2O1)C3=CC=CC=C3", SmilesParseOptions)
        .expect("fused chromanone should parse");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused chromanone should sanitize");

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused chromanone should canonicalize");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .mol
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
    let mut molecule = read_smiles_str(
        "C1=CC=C(C=C1)C2=C(C(=O)C3=CC=CC=C3O2)OC(=O)C4=CC5=C(C=C4Cl)SC6=NC=CN6S5(=O)=O",
        SmilesParseOptions,
    )
    .expect("conjugated benzopyrone should parse");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("conjugated benzopyrone should sanitize");

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("conjugated benzopyrone should canonicalize");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .mol
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
    let mut molecule = read_smiles_str(
        "C1=CC=C2C(=C1)C3=C(C2=O)C=C(C=C3)[N+]#N.C(=O)(C(F)(F)F)O",
        SmilesParseOptions,
    )
    .expect("fused fluorenone salt should parse");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused fluorenone salt should sanitize");

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("fused fluorenone salt should canonicalize");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    let aromatic_atoms = reparsed
        .mol
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(
        aromatic_atoms, 12,
        "canonical output should keep the fluorenone carbonyl bridge aliphatic: {written}"
    );
}

#[test]
fn partially_saturated_carbonyl_fused_rings_stay_aliphatic() {
    let mut molecule = read_smiles_str(
        "CC1(CC2=C(C(=O)C1)OC3=C(C2C4=CC=CC=C4[N+](=O)[O-])C(=O)CC(C3)(C)C)C",
        SmilesParseOptions,
    )
    .expect("partially saturated fused carbonyl system should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("partially saturated fused carbonyl system should sanitize");

    for atom_id in [1, 2, 3, 4, 5, 7, 9, 10, 21, 23, 24, 25] {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("partially saturated ring atom")
                .aromatic,
            "partially saturated ring atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [12, 13, 14, 15, 16, 17] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("phenyl atom")
                .aromatic,
            "phenyl atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_lactam_bridge_ring_stays_aliphatic() {
    let mut molecule = read_smiles_str(
        "CCN1C2=C(C=C(C=C2OC3=C(C1=O)C=CC=N3)C)[N+](=O)[O-]",
        SmilesParseOptions,
    )
    .expect("fused lactam bridge should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused lactam bridge should sanitize");

    for atom_id in [2, 9, 12] {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("lactam bridge atom")
                .aromatic,
            "lactam bridge atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [3, 4, 5, 6, 7, 8, 10, 11, 14, 15, 16, 17] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("aromatic fused atom")
                .aromatic,
            "aromatic fused atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn fused_four_member_diketone_ring_can_be_aromatic() {
    let mut molecule = read_smiles_str(
        "C1CSC2(C3=C(C=CC(=C3)Cl)OC4=C2C(=O)C4=O)SC1",
        SmilesParseOptions,
    )
    .expect("fused four-member diketone should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused four-member diketone should sanitize");

    for atom_id in [12, 13, 14, 16] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("four-member diketone atom")
                .aromatic,
            "four-member diketone atom {atom_id} should be aromatic"
        );
    }
}

#[test]
fn large_conjugated_macrocycle_aromatic_core_is_not_size_skipped() {
    let mut molecule = read_smiles_str(
            "CN(C)CCO.C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=N8)N=C2[N-]3.[Cu+2]",
            SmilesParseOptions,
        )
        .expect("macrocycle salt should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("macrocycle salt should sanitize");

    let aromatic_atoms = molecule
        .mol
        .atoms()
        .filter(|(_, atom)| atom.aromatic)
        .count();
    assert_eq!(aromatic_atoms, 40);

    let copper = molecule.mol.atom(AtomId::new(46)).expect("copper atom");
    assert!(!copper.aromatic);
    assert_eq!(copper.formal_charge, 2);
}

#[test]
fn neutral_aza_macrocycle_core_stays_aliphatic() {
    let mut molecule = read_smiles_str(
        "C1=CC=C2C(=C1)C3=NC4=NC(=NC5=NC(=NC6=NC(=NC2=N3)C7=CC=CC=C76)C8=CC=CC=C85)C9=CC=CC=C94",
        SmilesParseOptions,
    )
    .expect("neutral aza macrocycle should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("neutral aza macrocycle should sanitize");

    for atom_id in 6..=21 {
        assert!(
            !molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("neutral aza macrocycle atom")
                .aromatic,
            "neutral aza macrocycle atom {atom_id} should stay aliphatic"
        );
    }
    for atom_id in [0, 1, 2, 3, 4, 5, 22, 23, 24, 25, 26, 27] {
        assert!(
            molecule
                .mol
                .atom(AtomId::new(atom_id))
                .expect("benzene atom")
                .aromatic,
            "benzene atom {atom_id} should stay aromatic"
        );
    }
}

#[test]
fn fused_tertiary_amine_ring_does_not_extend_aromatic_core() {
    let mut molecule = read_smiles_str(
        "CC(C)C[C@@H]1CN2CCC3=CC(=C(C=C3C2CC1=O)OC)O[11CH3]",
        SmilesParseOptions,
    )
    .expect("fused tertiary amine record should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused tertiary amine record should sanitize");

    assert_eq!(
        molecule
            .mol
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        6
    );
    assert!(molecule
        .mol
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "N")
        .all(|(_, atom)| !atom.aromatic));

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));
    assert_eq!(
        reparsed
            .mol
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        6,
        "{written}"
    );
}

#[test]
fn fused_n_hydroxy_lactam_ring_stays_aromatic() {
    let mut molecule = read_smiles_str("CCCCCCCC1=CC2=C(C=C1)N(C=CC2=O)O", SmilesParseOptions)
        .expect("fused N-hydroxy lactam should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused N-hydroxy lactam should sanitize");

    let aromatic_atoms = molecule
        .mol
        .atoms()
        .filter_map(|(atom_id, atom)| {
            atom.aromatic
                .then_some((atom_id.index(), atom.element.symbol()))
        })
        .collect::<Vec<_>>();
    assert_eq!(aromatic_atoms.len(), 10, "{aromatic_atoms:?}");
    assert!(molecule
        .mol
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "N")
        .all(|(_, atom)| atom.aromatic));

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));
    assert_eq!(
        reparsed
            .mol
            .atoms()
            .filter(|(_, atom)| atom.aromatic)
            .count(),
        10,
        "{written}"
    );
}

#[test]
fn n_aryl_fused_pyrrole_ring_stays_aromatic() {
    let mut molecule = read_smiles_str(
        "CCOC(=O)C1=C(N(C2=C1C=C(C=C2)OCC(C[NH2+]CC3=CC=CC=C3)O)C4=CC=CC=C4)C.[Cl-]",
        SmilesParseOptions,
    )
    .expect("N-aryl fused pyrrole salt should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("N-aryl fused pyrrole salt should sanitize");

    let aromatic_neutral_nitrogens = molecule
        .mol
        .atoms()
        .filter(|(_, atom)| {
            atom.element.symbol() == "N" && atom.formal_charge == 0 && atom.aromatic
        })
        .count();
    assert_eq!(aromatic_neutral_nitrogens, 1);

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));
    let reparsed_aromatic_neutral_nitrogens = reparsed
        .mol
        .atoms()
        .filter(|(_, atom)| {
            atom.element.symbol() == "N" && atom.formal_charge == 0 && atom.aromatic
        })
        .count();
    assert_eq!(reparsed_aromatic_neutral_nitrogens, 1, "{written}");
}

#[test]
fn fused_saturated_thioether_bridge_stays_aliphatic() {
    let mut molecule = read_smiles_str(
        "C1=CC=C(C=C1)C2=C(C(=O)C3=CC=CC=C3O2)OC(=O)C4=CC5=C(C=C4Cl)SC6=NC=CN6S5(=O)=O",
        SmilesParseOptions,
    )
    .expect("fused thioether bridge should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("fused thioether bridge should sanitize");

    let neutral_sulfur_aromatic_count = molecule
        .mol
        .atoms()
        .filter(|(_, atom)| {
            atom.element.symbol() == "S" && atom.formal_charge == 0 && atom.aromatic
        })
        .count();
    assert_eq!(neutral_sulfur_aromatic_count, 0);

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));
    let reparsed_neutral_sulfur_aromatic_count = reparsed
        .mol
        .atoms()
        .filter(|(_, atom)| {
            atom.element.symbol() == "S" && atom.formal_charge == 0 && atom.aromatic
        })
        .count();
    assert_eq!(reparsed_neutral_sulfur_aromatic_count, 0, "{written}");
}

#[test]
fn canonical_smiles_prefers_sanitizable_lactone_candidate() {
    let mut molecule = read_smiles_str(
        "CC[C@H]1[C@H](COC1=O)CC2=CN=CN2C.C=CC(=O)O",
        SmilesParseOptions,
    )
    .expect("lactone imidazole mixture should parse");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("lactone imidazole mixture should sanitize");

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .unwrap_or_else(|_| panic!("canonical output should sanitize: {written}"));

    assert_eq!(reparsed.mol.atom_count(), molecule.mol.atom_count());
    assert_eq!(reparsed.mol.bond_count(), molecule.mol.bond_count());
}

#[test]
fn aromatic_pyridinium_smiles_sanitizes() {
    let mut molecule = read_smiles_str("CCCCCC(=O)C[n+]1ccccc1", SmilesParseOptions)
        .expect("aromatic pyridinium should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("aromatic pyridinium should sanitize");

    let cationic_nitrogen = molecule
        .mol
        .atoms()
        .find(|(_, atom)| atom.element.symbol() == "N" && atom.formal_charge > 0)
        .expect("pyridinium nitrogen should exist")
        .1;
    assert!(cationic_nitrogen.aromatic);
}

#[test]
fn aromatic_pyrone_canonical_smiles_sanitizes() {
    let mut molecule = read_smiles_str("CC#CC#Cc1cccc(=O)o1", SmilesParseOptions)
        .expect("aromatic pyrone should parse");

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("aromatic pyrone should sanitize");
}

#[test]
fn canonical_smiles_preserves_metal_bound_bracket_hydrogens() {
    let mut molecule =
        read_smiles_str("CC[Hg+]", SmilesParseOptions).expect("organomercury SMILES parses");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("organomercury SMILES sanitizes");

    let written = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)
        .expect("canonical SMILES should write");

    assert!(written.contains("[CH2][Hg+]"), "{written}");
    let mut reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("canonical output should parse");
    sanitize_small_molecule(&mut reparsed, SanitizeOptions::default())
        .expect("canonical output should sanitize");
    let metal_bound_carbon = reparsed
        .mol
        .atoms()
        .find(|(atom_id, atom)| {
            atom.element.symbol() == "C"
                && reparsed
                    .mol
                    .incident_bonds(*atom_id)
                    .expect("atom should be live")
                    .any(|(_, bond)| {
                        let neighbor_id = bond.other_atom(*atom_id);
                        reparsed
                            .mol
                            .atom(neighbor_id)
                            .is_ok_and(|neighbor| neighbor.element.symbol() == "Hg")
                    })
        })
        .expect("canonical output should retain a carbon-mercury bond")
        .1;
    assert_eq!(metal_bound_carbon.explicit_hydrogens, 2);
    assert!(metal_bound_carbon.no_implicit_hydrogens);
    assert_eq!(metal_bound_carbon.implicit_hydrogens, Some(0));
}

#[test]
fn smiles_writer_rejects_lossy_bonds_and_stereo() {
    let mut molecule = SmallMolecule::default();
    let a = molecule.mol.add_atom(carbon());
    let b = molecule.mol.add_atom(carbon());
    let bond = molecule
        .mol
        .add_bond(a, b, BondOrder::Dative)
        .expect("bond");
    assert!(write_smiles(&molecule, SmilesWriteOptions)
        .expect_err("dative bond should be rejected")
        .message
        .contains("cannot encode"));

    molecule.mol.bond_mut(bond).expect("bond").order = BondOrder::Single;
    molecule.mol.bond_mut(bond).expect("bond").stereo = Some(BondStereo::Up);
    assert!(write_smiles(&molecule, SmilesWriteOptions)
        .expect_err("stereo should be rejected")
        .message
        .contains("stereochemistry"));

    molecule.mol.bond_mut(bond).expect("bond").stereo = None;
    molecule.mol.atom_mut(a).expect("atom").chiral = Some(AtomStereo::TetrahedralClockwise);
    assert!(write_smiles(&molecule, SmilesWriteOptions)
        .expect_err("atom chirality should be rejected")
        .message
        .contains("atom stereochemistry"));

    molecule.mol.atom_mut(a).expect("atom").chiral = None;
    molecule.mol.atom_mut(a).expect("atom").radical = Some(AtomRadical::Doublet);
    assert!(write_smiles(&molecule, SmilesWriteOptions)
        .expect_err("radical should be rejected")
        .message
        .contains("radicals"));

    molecule.mol.atom_mut(a).expect("atom").radical = None;
    molecule
        .mol
        .atom_mut(a)
        .expect("atom")
        .no_implicit_hydrogens = true;
    let written = write_smiles(&molecule, SmilesWriteOptions)
        .expect("no-implicit-hydrogen atom should write");
    assert!(written.contains("[C]"));
    let reparsed =
        read_smiles_str(&written, SmilesParseOptions).expect("writer output should parse");
    assert!(reparsed
        .mol
        .atoms()
        .any(|(_, atom)| atom.no_implicit_hydrogens));
}

#[test]
fn smiles_writer_rejects_more_ring_labels_than_parser_supports() {
    let mut molecule = SmallMolecule::default();
    let atoms = (0..16)
        .map(|_| molecule.mol.add_atom(carbon()))
        .collect::<Vec<_>>();
    for left in 0..atoms.len() {
        for right in (left + 1)..atoms.len() {
            molecule
                .mol
                .add_bond(atoms[left], atoms[right], BondOrder::Single)
                .expect("complete graph bond should be valid");
        }
    }

    assert!(write_smiles(&molecule, SmilesWriteOptions)
        .expect_err("more than 99 ring closures should be rejected")
        .message
        .contains("at most 99"));
}
