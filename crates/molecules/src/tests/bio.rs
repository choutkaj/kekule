use super::*;

#[test]
fn bio_hierarchy_adds_models_chains_residues_and_atom_sites() {
    let mut hierarchy = BioHierarchy::new();
    let model = hierarchy.add_model("1");
    let chain = hierarchy
        .add_chain(model, "A", Some("authA".to_owned()))
        .expect("chain should be valid");
    let residue = hierarchy
        .add_residue(
            chain,
            "GLY",
            Some(10),
            Some("42".to_owned()),
            Some("A".to_owned()),
        )
        .expect("residue should be valid");
    let metadata = AtomSiteMetadata {
        group_pdb: Some("ATOM".to_owned()),
        atom_site_id: Some("1".to_owned()),
        type_symbol: Some("C".to_owned()),
        label_asym_id: Some("A".to_owned()),
        auth_asym_id: Some("authA".to_owned()),
        label_atom_id: Some("CA".to_owned()),
        auth_atom_id: Some("CAY".to_owned()),
        label_alt_id: Some("B".to_owned()),
        occupancy: Some(0.5),
        occupancy_raw: Some("0.50".to_owned()),
        b_factor: Some(12.25),
        b_factor_raw: Some("12.25".to_owned()),
        cartn_x_raw: None,
        cartn_y_raw: None,
        cartn_z_raw: None,
    };
    let site = hierarchy
        .add_atom_site(residue, AtomId::new(7), metadata.clone())
        .expect("atom site should be valid");

    assert_eq!(model.raw(), 0);
    assert_eq!(chain.raw(), 0);
    assert_eq!(residue.raw(), 0);
    assert_eq!(site.raw(), 0);
    assert_eq!(
        hierarchy.model(model).expect("model exists").chains,
        vec![chain]
    );
    assert_eq!(
        hierarchy.chain(chain).expect("chain exists").residues,
        vec![residue]
    );
    assert_eq!(
        hierarchy
            .residue(residue)
            .expect("residue exists")
            .atom_sites,
        vec![site]
    );
    assert_eq!(
        hierarchy
            .atom_site_for_atom(AtomId::new(7))
            .expect("site exists")
            .metadata,
        metadata
    );
}

#[test]
fn bio_hierarchy_iteration_is_insertion_order() {
    let mut hierarchy = BioHierarchy::new();
    let first_model = hierarchy.add_model("1");
    let second_model = hierarchy.add_model("2");
    let first_chain = hierarchy.add_chain(first_model, "A", None).expect("chain");
    let second_chain = hierarchy.add_chain(second_model, "B", None).expect("chain");

    assert_eq!(
        hierarchy.models().map(|(id, _)| id).collect::<Vec<_>>(),
        vec![first_model, second_model]
    );
    assert_eq!(
        hierarchy.chains().map(|(id, _)| id).collect::<Vec<_>>(),
        vec![first_chain, second_chain]
    );
}

#[test]
fn bio_hierarchy_rejects_missing_parents_and_duplicate_atom_placement() {
    let mut hierarchy = BioHierarchy::new();
    assert_eq!(
        hierarchy
            .add_chain(ModelId::new(99), "A", None)
            .expect_err("missing model should fail"),
        BioHierarchyError::InvalidModelId(ModelId::new(99))
    );

    let model = hierarchy.add_model("1");
    let chain = hierarchy.add_chain(model, "A", None).expect("chain");
    assert_eq!(
        hierarchy
            .add_residue(ChainId::new(99), "GLY", None, None, None)
            .expect_err("missing chain should fail"),
        BioHierarchyError::InvalidChainId(ChainId::new(99))
    );
    let residue = hierarchy
        .add_residue(chain, "GLY", None, None, None)
        .expect("residue");
    let atom = AtomId::new(2);
    hierarchy
        .add_atom_site(residue, atom, AtomSiteMetadata::default())
        .expect("first placement should work");
    assert_eq!(
        hierarchy
            .add_atom_site(residue, atom, AtomSiteMetadata::default())
            .expect_err("duplicate atom placement should fail"),
        BioHierarchyError::DuplicateAtomPlacement(atom)
    );
}

#[test]
fn macro_molecule_validates_atom_site_atom_ids() {
    let mut macro_mol = MacroMolecule::default();
    let atom = macro_mol.graph_mut().add_atom(carbon());
    let model = macro_mol.hierarchy_mut().add_model("1");
    let chain = macro_mol
        .hierarchy_mut()
        .add_chain(model, "A", Some("authA".to_owned()))
        .expect("chain");
    let residue = macro_mol
        .hierarchy_mut()
        .add_residue(chain, "ALA", Some(1), Some("1".to_owned()), None)
        .expect("residue");

    macro_mol
        .add_atom_site(
            residue,
            atom,
            AtomSiteMetadata {
                group_pdb: Some("ATOM".to_owned()),
                atom_site_id: Some("1".to_owned()),
                type_symbol: Some("C".to_owned()),
                label_asym_id: Some("A".to_owned()),
                auth_asym_id: Some("authA".to_owned()),
                label_atom_id: Some("CA".to_owned()),
                auth_atom_id: Some("CA".to_owned()),
                label_alt_id: None,
                occupancy: Some(1.0),
                occupancy_raw: Some("1.0".to_owned()),
                b_factor: Some(10.0),
                b_factor_raw: Some("10.0".to_owned()),
                cartn_x_raw: None,
                cartn_y_raw: None,
                cartn_z_raw: None,
            },
        )
        .expect("valid atom should attach");
    assert_eq!(
        macro_mol
            .add_atom_site(residue, AtomId::new(99), AtomSiteMetadata::default())
            .expect_err("missing atom should fail"),
        BioHierarchyError::InvalidAtomId(AtomId::new(99))
    );
}

#[test]
fn macro_molecule_validates_and_sanitizes_separate_from_small_molecule_chemistry() {
    let mut macro_mol = MacroMolecule::default();
    let atom = macro_mol.graph_mut().add_atom(carbon());
    let mut conformer = Conformer::new();
    conformer.set_position(atom, Point3::new(1.0, 2.0, 3.0));
    macro_mol.graph_mut().add_conformer(conformer);

    let model = macro_mol.hierarchy_mut().add_model("1");
    let chain = macro_mol
        .hierarchy_mut()
        .add_chain(model, "A", None)
        .expect("chain");
    let residue = macro_mol
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), Some("1".to_owned()), None)
        .expect("residue");
    macro_mol
        .add_atom_site(
            residue,
            atom,
            AtomSiteMetadata {
                occupancy: Some(1.0),
                b_factor: Some(12.0),
                ..AtomSiteMetadata::default()
            },
        )
        .expect("atom site");

    assert_eq!(macro_mol.models().count(), 1);
    assert_eq!(macro_mol.chains().count(), 1);
    assert_eq!(macro_mol.residues().count(), 1);
    assert_eq!(macro_mol.atom_sites().count(), 1);
    assert_eq!(macro_mol.atom_site_for_atom(atom).expect("site").atom, atom);

    let report = macro_mol.validate().expect("macro molecule validates");
    assert_eq!(report.models_checked, 1);
    assert_eq!(report.chains_checked, 1);
    assert_eq!(report.residues_checked, 1);
    assert_eq!(report.atom_sites_checked, 1);
    assert_eq!(report.conformers_checked, 1);
    assert_eq!(report.coordinates_checked, 1);

    let sanitize = macro_mol.sanitize().expect("macro molecule sanitizes");
    assert_eq!(sanitize.validation, Some(report));
    assert_eq!(sanitize.normalized_atom_sites, 0);
    assert_eq!(sanitize.recognized_residues, 0);
    assert_eq!(sanitize.assigned_bonds, 0);
    assert_eq!(macro_mol.graph().bond_count(), 0);
}

#[test]
fn macro_molecule_validation_rejects_cross_layer_inconsistency() {
    let graph = Molecule::new();
    let mut hierarchy = BioHierarchy::new();
    let model = hierarchy.add_model("1");
    let chain = hierarchy.add_chain(model, "A", None).expect("chain");
    let residue = hierarchy
        .add_residue(chain, "GLY", None, None, None)
        .expect("residue");
    let site = hierarchy
        .add_atom_site(residue, AtomId::new(0), AtomSiteMetadata::default())
        .expect("hierarchy accepts graph-external atom ids");
    let macro_mol = MacroMolecule::from_parts(graph, hierarchy);

    assert_eq!(
        macro_mol.validate().expect_err("graph-external atom fails"),
        MacroValidateError::InvalidAtomSiteAtom {
            site,
            atom: AtomId::new(0)
        }
    );
}

#[test]
fn macro_molecule_sanitize_rejects_unsupported_preparation_options() {
    let mut macro_mol = MacroMolecule::default();
    let atom = macro_mol.graph_mut().add_atom(carbon());
    let model = macro_mol.hierarchy_mut().add_model("1");
    let chain = macro_mol
        .hierarchy_mut()
        .add_chain(model, "A", None)
        .expect("chain");
    let residue = macro_mol
        .hierarchy_mut()
        .add_residue(chain, "LIG", None, None, None)
        .expect("residue");
    macro_mol
        .add_atom_site(residue, atom, AtomSiteMetadata::default())
        .expect("atom site");

    let options = MacroSanitizeOptions {
        ligand_policy: LigandSanitizePolicy::SanitizeAllDisconnectedComponents,
        ..MacroSanitizeOptions::default()
    };

    assert_eq!(
        macro_mol
            .sanitize_with_options(options)
            .expect_err("unsupported ligand policy fails"),
        MacroSanitizeError::UnsupportedOption(
            "bond, disulfide, or ligand sanitization is not implemented"
        )
    );
}

#[test]
fn core_atom_does_not_store_biomolecular_labels() {
    let atom = carbon();

    assert_eq!(atom.element.symbol(), "C");
    assert!(atom.props.is_empty());
}

#[test]
fn mmcif_parse_builds_macro_molecule_hierarchy() {
    let input = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA CAY . GLY GLY A X 10 42 A 0.50 12.25 1.25 2.50 3.75 1
ATOM 2 O O O . GLY GLY A X 10 42 A 1.00 10.00 4.25 5.50 6.75 1
"#;

    let macro_mol =
        bio_api::read_mmcif_str(input, MmcifParseOptions::default()).expect("mmCIF should parse");

    assert_eq!(macro_mol.graph().atom_count(), 2);
    assert_eq!(macro_mol.graph().bond_count(), 0);
    assert_eq!(macro_mol.hierarchy().models().count(), 1);
    assert_eq!(macro_mol.hierarchy().chains().count(), 1);
    assert_eq!(macro_mol.hierarchy().residues().count(), 1);
    assert_eq!(macro_mol.hierarchy().atom_sites().count(), 2);
    let (_, chain) = macro_mol.hierarchy().chains().next().expect("chain exists");
    assert_eq!(chain.label_id, "A");
    assert_eq!(chain.author_id, Some("X".to_owned()));
    let (_, residue) = macro_mol
        .hierarchy()
        .residues()
        .next()
        .expect("residue exists");
    assert_eq!(residue.name, "GLY");
    assert_eq!(residue.label_comp_id, Some("GLY".to_owned()));
    assert_eq!(residue.author_comp_id, Some("GLY".to_owned()));
    assert_eq!(residue.label_seq_id, Some(10));
    assert_eq!(residue.author_seq_id, Some("42".to_owned()));
    assert_eq!(residue.insertion_code, Some("A".to_owned()));
    let site = macro_mol
        .hierarchy()
        .atom_site_for_atom(AtomId::new(0))
        .expect("site exists");
    assert_eq!(site.metadata.label_atom_id, Some("CA".to_owned()));
    assert_eq!(site.metadata.auth_atom_id, Some("CAY".to_owned()));
    assert_eq!(site.metadata.occupancy, Some(0.5));
    assert_eq!(site.metadata.b_factor, Some(12.25));
    let (_, conformer) = macro_mol
        .graph()
        .first_conformer()
        .expect("conformer exists");
    assert_eq!(
        conformer.position(AtomId::new(0)),
        Some(Point3::new(1.25, 2.50, 3.75))
    );
    assert_eq!(
        conformer.position(AtomId::new(1)),
        Some(Point3::new(4.25, 5.50, 6.75))
    );
}

#[test]
fn mmcif_parse_handles_missing_values_and_quotes() {
    let input = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
C "C A" ? "LIG" "AA" . 7 ? ? ? . 2
"#;

    let macro_mol =
        bio_api::read_mmcif_str(input, MmcifParseOptions::default()).expect("mmCIF should parse");
    let (_, model) = macro_mol.hierarchy().models().next().expect("model exists");
    let site = macro_mol
        .hierarchy()
        .atom_site_for_atom(AtomId::new(0))
        .expect("site exists");

    assert_eq!(model.model_id, "2");
    assert_eq!(site.metadata.label_atom_id, Some("C A".to_owned()));
    assert_eq!(site.metadata.auth_atom_id, None);
    assert_eq!(site.metadata.label_alt_id, None);
    assert_eq!(site.metadata.occupancy, None);
    assert_eq!(site.metadata.b_factor, None);
    assert!(macro_mol.graph().first_conformer().is_none());
}

#[test]
fn mmcif_auth_sequence_distinguishes_label_sequence_less_residues() {
    let input = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
O O HOH HOH A W . 10
O O HOH HOH A W . 11
"#;

    let macro_mol =
        bio_api::read_mmcif_str(input, MmcifParseOptions::default()).expect("mmCIF should parse");

    assert_eq!(macro_mol.hierarchy().residues().count(), 2);
    let seq_ids = macro_mol
        .hierarchy()
        .residues()
        .map(|(_, residue)| residue.author_seq_id.clone())
        .collect::<Vec<_>>();
    assert_eq!(seq_ids, vec![Some("10".to_owned()), Some("11".to_owned())]);
}

#[test]
fn mmcif_lenient_occurrences_group_atoms_and_keep_altlocs_together() {
    let input = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
C C1 AC1 LBL AUT A X . . ? A
C C1 AC1 LBL AUT A X . . ? B
O O1 AO1 LBL AUT A X . . ? .
C C1 AC1 LBL AUT A X . . ? .
O O1 AO1 LBL AUT A X . . ? .
"#;

    let macro_mol = bio_api::read_mmcif_str(
        input,
        MmcifParseOptions {
            strict: false,
            ..MmcifParseOptions::default()
        },
    )
    .expect("lenient ambiguous residues should parse");

    let residues = macro_mol
        .hierarchy()
        .residues()
        .map(|(_, residue)| residue)
        .collect::<Vec<_>>();
    assert_eq!(residues.len(), 2);
    assert_eq!(residues[0].atom_sites.len(), 3);
    assert_eq!(residues[1].atom_sites.len(), 2);
    assert_eq!(residues[0].label_comp_id, Some("LBL".to_owned()));
    assert_eq!(residues[0].author_comp_id, Some("AUT".to_owned()));
    let first_altlocs = residues[0]
        .atom_sites
        .iter()
        .filter_map(|site_id| {
            macro_mol
                .hierarchy()
                .atom_site(*site_id)
                .ok()
                .and_then(|site| site.metadata.label_alt_id.clone())
        })
        .collect::<Vec<_>>();
    assert_eq!(first_altlocs, vec!["A", "B"]);
}

#[test]
fn mmcif_insertion_codes_and_models_keep_distinct_residues_and_coordinates() {
    let input = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
C C1 LIG LG1 A X 7 70 A 1.0 2.0 3.0 1
O O1 LIG LG1 A X 7 70 B 4.0 5.0 6.0 1
N N1 LIG LG1 A X 7 70 A 7.0 8.0 9.0 2
"#;

    let macro_mol =
        bio_api::read_mmcif_str(input, MmcifParseOptions::default()).expect("mmCIF should parse");

    assert_eq!(macro_mol.hierarchy().models().count(), 2);
    assert_eq!(macro_mol.hierarchy().residues().count(), 3);
    assert_eq!(
        macro_mol
            .hierarchy()
            .residues()
            .map(|(_, residue)| residue.insertion_code.clone())
            .collect::<Vec<_>>(),
        vec![
            Some("A".to_owned()),
            Some("B".to_owned()),
            Some("A".to_owned())
        ]
    );
    let (_, conformer) = macro_mol
        .graph()
        .first_conformer()
        .expect("conformer exists");
    assert_eq!(
        macro_mol
            .graph()
            .atom_ids()
            .filter_map(|atom_id| conformer.position(atom_id))
            .collect::<Vec<_>>(),
        vec![
            Point3::new(1.0, 2.0, 3.0),
            Point3::new(4.0, 5.0, 6.0),
            Point3::new(7.0, 8.0, 9.0)
        ]
    );
    assert_eq!(macro_mol.graph().bond_count(), 0);
    assert_eq!(macro_mol.graph().perception().rings, ComputedState::Absent);
}

#[test]
fn mmcif_strict_mode_rejects_partial_coordinates() {
    let input = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
C C1 BEN A 1 1.0 2.0
"#;

    let err = bio_api::read_mmcif_str(input, MmcifParseOptions::default())
        .expect_err("partial coordinates should fail");
    assert!(err.message.contains("partial atom-site coordinate"));
}

#[test]
fn mmcif_parse_rejects_missing_strict_atom_id_and_unknown_element() {
    let missing_atom_id = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_comp_id
_atom_site.label_asym_id
C GLY A
"#;
    let err = bio_api::read_mmcif_str(missing_atom_id, MmcifParseOptions::default())
        .expect_err("strict mode should require label atom id");
    assert!(err.message.contains("label atom id"));

    let unknown_element = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
Xx CA GLY A
"#;
    let err = bio_api::read_mmcif_str(unknown_element, MmcifParseOptions::default())
        .expect_err("unknown element should fail");
    assert!(err.message.contains("unknown atom-site element"));
}

#[test]
fn malformed_mmcif_returns_located_errors_without_panicking() {
    let cases = [
            ("unterminated quote", "data_x\n_tag 'unterminated\n", 2),
            ("unterminated semicolon", "data_x\n;\nvalue\n", 2),
            (
                "ragged loop",
                "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\nC C1 O\n",
                5,
            ),
            (
                "duplicate atom-site tag",
                "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.type_symbol\n_atom_site.label_atom_id\nC C C1\n",
                4,
            ),
            (
                "integer overflow",
                "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\nC C1 BEN A 999999999999999999999\n",
                8,
            ),
            (
                "float overflow",
                "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nC C1 BEN A 1 1e999 0 0\n",
                11,
            ),
        ];

    for (name, input, expected_line) in cases {
        let parsed = std::panic::catch_unwind(|| {
            bio_api::read_mmcif_str(input, MmcifParseOptions::default())
        })
        .unwrap_or_else(|_| panic!("{name} panicked"));
        let error = parsed.expect_err("malformed mmCIF should fail");
        assert_eq!(error.line, expected_line, "line for {name}");
        assert!(!error.message.is_empty(), "message for {name}");
    }
}

#[test]
fn deterministic_parser_fuzz_smoke_is_panic_free() {
    let mol_seed = "Methane\n  molecules\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n";
    for input in deterministic_text_mutations(mol_seed) {
        std::panic::catch_unwind(|| {
            if let Ok(molecule) = molfile::read_v2000_str(&input) {
                if let Ok(output) = molfile::write_v2000(&molecule) {
                    let _ = molfile::read_v2000_str(&output);
                }
            }
        })
        .expect("Molfile parser smoke mutation panicked");
    }

    let sdf_seed = format!("{mol_seed}$$$$\n");
    for input in deterministic_text_mutations(&sdf_seed) {
        std::panic::catch_unwind(|| {
            if let Ok(records) = sdf::read_v2000_records(
                &input,
                SdfParseOptions {
                    allow_missing_final_delimiter: true,
                },
            ) {
                let molecules = records
                    .into_iter()
                    .map(|record| record.molecule)
                    .collect::<Vec<_>>();
                if let Ok(output) = sdf::write_v2000(&molecules) {
                    let _ = sdf::read_v2000_records(&output, SdfParseOptions::default());
                }
            }
        })
        .expect("SDF parser smoke mutation panicked");
    }

    for input in deterministic_text_mutations("CC(=O)O") {
        std::panic::catch_unwind(|| {
            if let Ok(molecule) = smiles_api::read_str_with_options(&input, SmilesParseOptions) {
                if let Ok(output) = smiles_api::write_with_options(&molecule, SmilesWriteOptions) {
                    let _ = smiles_api::read_str_with_options(&output, SmilesParseOptions);
                }
            }
        })
        .expect("SMILES parser smoke mutation panicked");
    }

    let mmcif_seed = "data_tiny\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\nC C1 LIG A 1\n";
    for input in deterministic_text_mutations(mmcif_seed) {
        std::panic::catch_unwind(|| {
            if let Ok(molecule) = bio_api::read_mmcif_str(&input, MmcifParseOptions::default()) {
                for atom in molecule.graph().atom_ids() {
                    let _ = molecule.graph().atom(atom);
                    let _ = molecule.hierarchy().atom_site_for_atom(atom);
                }
            }
        })
        .expect("mmCIF parser smoke mutation panicked");
    }
}

#[test]
fn mmcif_parse_options_enforce_documented_resource_limits() {
    let input = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\nC C1 BEN A 1\n";

    let input_error = bio_api::read_mmcif_str(
        input,
        MmcifParseOptions {
            max_input_bytes: input.len() - 1,
            ..MmcifParseOptions::default()
        },
    )
    .expect_err("input byte limit should fail");
    assert!(input_error.message.contains("byte limit"));

    let token_error = bio_api::read_mmcif_str(
        input,
        MmcifParseOptions {
            max_tokens: 2,
            ..MmcifParseOptions::default()
        },
    )
    .expect_err("token count limit should fail");
    assert!(token_error.message.contains("token count"));

    let value_error = bio_api::read_mmcif_str(
        input,
        MmcifParseOptions {
            max_token_bytes: 2,
            ..MmcifParseOptions::default()
        },
    )
    .expect_err("token byte limit should fail");
    assert!(value_error.message.contains("token limit"));

    let row_error = bio_api::read_mmcif_str(
        input,
        MmcifParseOptions {
            max_atom_site_rows: 0,
            ..MmcifParseOptions::default()
        },
    )
    .expect_err("row limit should fail");
    assert!(row_error.message.contains("row limit"));
}

#[test]
fn mmcif_parse_does_not_infer_bonds_or_perception() {
    let input = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
C C1 BEN A 1
C C2 BEN A 1
"#;

    let macro_mol =
        bio_api::read_mmcif_str(input, MmcifParseOptions::default()).expect("mmCIF should parse");

    assert_eq!(macro_mol.graph().atom_count(), 2);
    assert_eq!(macro_mol.graph().bond_count(), 0);
    assert_eq!(macro_mol.graph().perception().rings, ComputedState::Absent);
    assert_eq!(
        macro_mol.graph().perception().aromaticity,
        ComputedState::Absent
    );
}

#[test]
fn wrappers_share_the_core_molecule_graph() {
    let mut small = SmallMolecule::default();
    let a = small.graph_mut().add_atom(carbon());
    let b = small.graph_mut().add_atom(oxygen());
    small
        .graph_mut()
        .add_bond(a, b, BondOrder::Single)
        .expect("small molecule graph should accept bonds");

    let mut macro_mol = MacroMolecule::default();
    let c = macro_mol.graph_mut().add_atom(carbon());

    assert_eq!(small.graph().atom_count(), 2);
    assert_eq!(small.graph().bond_count(), 1);
    assert_eq!(
        macro_mol
            .graph()
            .atom(c)
            .expect("macro atom exists")
            .element
            .symbol(),
        "C"
    );
}
