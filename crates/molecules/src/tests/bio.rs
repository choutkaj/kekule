use super::*;

#[test]
fn smcra_hierarchy_adds_models_chains_residues_and_atom_sites() {
    let mut hierarchy = SmcraHierarchy::new();
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
    let metadata = SmcraAtomSiteMetadata {
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
fn smcra_hierarchy_iteration_is_insertion_order() {
    let mut hierarchy = SmcraHierarchy::new();
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
fn smcra_hierarchy_rejects_missing_parents_and_duplicate_atom_placement() {
    let mut hierarchy = SmcraHierarchy::new();
    assert_eq!(
        hierarchy
            .add_chain(SmcraModelId::new(99), "A", None)
            .expect_err("missing model should fail"),
        SmcraHierarchyError::InvalidModelId(SmcraModelId::new(99))
    );

    let model = hierarchy.add_model("1");
    let chain = hierarchy.add_chain(model, "A", None).expect("chain");
    assert_eq!(
        hierarchy
            .add_residue(SmcraChainId::new(99), "GLY", None, None, None)
            .expect_err("missing chain should fail"),
        SmcraHierarchyError::InvalidChainId(SmcraChainId::new(99))
    );
    let residue = hierarchy
        .add_residue(chain, "GLY", None, None, None)
        .expect("residue");
    let atom = AtomId::new(2);
    hierarchy
        .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
        .expect("first placement should work");
    assert_eq!(
        hierarchy
            .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
            .expect_err("duplicate atom placement should fail"),
        SmcraHierarchyError::DuplicateAtomPlacement(atom)
    );
}

#[test]
fn macro_molecule_validates_atom_site_atom_ids() {
    let mut builder = MacroMolecule::builder();
    let atom = builder.graph_mut().add_atom(carbon());
    let model = builder.hierarchy_mut().add_model("1");
    let chain = builder
        .hierarchy_mut()
        .add_chain(model, "A", Some("authA".to_owned()))
        .expect("chain");
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "ALA", Some(1), Some("1".to_owned()), None)
        .expect("residue");

    builder
        .add_atom_site(
            residue,
            atom,
            SmcraAtomSiteMetadata {
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
        builder
            .add_atom_site(residue, AtomId::new(99), SmcraAtomSiteMetadata::default())
            .expect_err("missing atom should fail"),
        SmcraHierarchyError::InvalidAtomId(AtomId::new(99))
    );
    builder.build().expect("checked macro molecule");
}

fn macro_molecule_with_valid_atom_site() -> (MacroMolecule, AtomId) {
    let mut builder = MacroMolecule::builder();
    let atom = builder.graph_mut().add_atom(carbon());
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            atom,
            crate::units::Quantity::new(Point3::new(1.0, 2.0, 3.0), crate::units::ANGSTROM),
        )
        .unwrap();
    builder
        .graph_mut()
        .add_conformer(conformer)
        .expect("valid conformer");

    let model = builder.hierarchy_mut().add_model("1");
    let chain = builder
        .hierarchy_mut()
        .add_chain(model, "A", None)
        .expect("chain");
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), Some("1".to_owned()), None)
        .expect("residue");
    builder
        .add_atom_site(
            residue,
            atom,
            SmcraAtomSiteMetadata {
                occupancy: Some(1.0),
                b_factor: Some(12.0),
                ..SmcraAtomSiteMetadata::default()
            },
        )
        .expect("atom site");

    (builder.build().expect("checked macro molecule"), atom)
}

#[test]
fn macro_molecule_validates_separately_from_small_molecule_chemistry() {
    let (macro_mol, atom) = macro_molecule_with_valid_atom_site();

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
    assert_eq!(macro_mol.graph().bond_count(), 0);
}

#[test]
fn macro_molecule_validation_rejects_cross_layer_inconsistency() {
    let graph = Molecule::new();
    let mut hierarchy = SmcraHierarchy::new();
    let model = hierarchy.add_model("1");
    let chain = hierarchy.add_chain(model, "A", None).expect("chain");
    let residue = hierarchy
        .add_residue(chain, "GLY", None, None, None)
        .expect("residue");
    let site = hierarchy
        .add_atom_site(residue, AtomId::new(0), SmcraAtomSiteMetadata::default())
        .expect("hierarchy accepts graph-external atom ids");
    assert_eq!(
        MacroMolecule::try_from_parts(graph, hierarchy).expect_err("graph-external atom fails"),
        MacroValidateError::InvalidAtomSiteAtom {
            site,
            atom: AtomId::new(0)
        }
    );
}

#[test]
fn macro_molecule_editor_is_transactional() {
    let (mut macro_mol, atom) = macro_molecule_with_valid_atom_site();
    let before = macro_mol.clone();
    let mut editor = macro_mol.edit();
    editor
        .graph_mut()
        .delete_atom(atom)
        .expect("staged graph atom exists");
    assert_eq!(
        editor.commit().expect_err("orphaned hierarchy must fail"),
        MacroValidateError::InvalidAtomSiteAtom {
            site: SmcraAtomSiteId::new(0),
            atom,
        }
    );
    assert_eq!(macro_mol, before);
}

#[test]
fn core_atom_does_not_store_biomolecular_labels() {
    let atom = carbon();

    assert_eq!(atom.element.symbol(), "C");
    assert!(atom.props.is_empty());
}

#[test]
fn deterministic_parser_fuzz_smoke_is_panic_free() {
    let mol_seed = "Methane\n  molecules\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n";
    for input in deterministic_text_mutations(mol_seed) {
        std::panic::catch_unwind(|| {
            if let Ok(molecule) = read_molfile(&input) {
                if let Ok(output) = molfile::write_v2000(&molecule) {
                    let _ = read_molfile(&output);
                }
            }
        })
        .expect("Molfile parser smoke mutation panicked");
    }

    let sdf_seed = format!("{mol_seed}$$$$\n");
    for input in deterministic_text_mutations(&sdf_seed) {
        std::panic::catch_unwind(|| {
            if let Ok(records) = read_sdf_records_with_options(
                &input,
                SdfParseOptions {
                    allow_missing_final_delimiter: true,
                },
            ) {
                if let Ok(output) = sdf::write_v2000(&records) {
                    let _ = read_sdf_records(&output);
                }
            }
        })
        .expect("SDF parser smoke mutation panicked");
    }

    for input in deterministic_text_mutations("CC(=O)O") {
        std::panic::catch_unwind(|| {
            if let Ok(molecule) = read_smiles(&input) {
                if let Ok(output) = smiles_api::write_with_options(&molecule, SmilesWriteOptions) {
                    let _ = read_smiles(&output);
                }
            }
        })
        .expect("SMILES parser smoke mutation panicked");
    }
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

    let mut builder = MacroMolecule::builder();
    let c = builder.graph_mut().add_atom(carbon());
    let model = builder.hierarchy_mut().add_model("1");
    let chain = builder
        .hierarchy_mut()
        .add_chain(model, "A", None)
        .expect("chain");
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .expect("residue");
    builder
        .add_atom_site(residue, c, SmcraAtomSiteMetadata::default())
        .expect("site");
    let macro_mol = builder.build().expect("checked macro molecule");

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
