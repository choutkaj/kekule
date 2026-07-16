use molecules::prelude::*;

#[test]
fn quantity_and_unit_public_api() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::units::{Dimension, Quantity, Unit, ANGSTROM, NANOMETER};

    let length = 1.0 * NANOMETER;
    assert_eq!(length.value_in(ANGSTROM)?, 10.0);

    let picometer = Unit::new(Dimension::LENGTH, 1.0e-12, Some("pm"))?;
    let coordinates = Quantity::new(vec![[100.0, 200.0, 300.0]], picometer);
    assert_eq!(coordinates.value_in(ANGSTROM)?, vec![[1.0, 2.0, 3.0]]);
    Ok(())
}

#[test]
fn small_molecule_happy_path() -> Result<(), Box<dyn std::error::Error>> {
    let mut mol = SmallMolecule::from_smiles("c1ccccc1O")?;
    mol.sanitize()?;
    assert_eq!(mol.atom_count(), 7);
    assert_eq!(mol.bond_count(), 7);
    let smiles = mol.to_canonical_smiles()?;
    assert!(!smiles.is_empty());
    Ok(())
}

#[test]
fn namespaced_small_molecule_api() -> Result<(), Box<dyn std::error::Error>> {
    let document = molecules::smiles::parse_str("CC(=O)O")?;
    let interpreted = molecules::smiles::interpret(&document)?;
    assert_eq!(interpreted.report().atom_mappings().len(), 4);
    assert_eq!(interpreted.report().bond_mappings().len(), 3);
    let mut mol = interpreted.into_molecule();
    molecules::perception::sanitize(&mut mol)?;
    let smiles = molecules::smiles::write_canonical(&mol)?;
    assert!(!smiles.is_empty());
    Ok(())
}

#[test]
fn molfile_and_sdf_interpretation_reports_are_public() -> Result<(), Box<dyn std::error::Error>> {
    let molfile = "\
Report
molecules

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
";
    let document = molecules::molfile::parse_str(molfile)?;
    let interpreted = molecules::molfile::interpret(&document)?;
    assert_eq!(interpreted.report().atom_mappings().len(), 1);
    assert_eq!(interpreted.report().atom_mappings()[0].source_line(), 5);

    let sdf = format!("{molfile}>  <SOURCE>\nrelease-test\n\n$$$$\n");
    let document = molecules::sdf::parse_str(&sdf, molecules::sdf::SdfParseOptions::default())?;
    let interpreted = molecules::sdf::interpret(&document)?;
    assert_eq!(interpreted.records().len(), 1);
    assert_eq!(interpreted.report().records()[0].record(), 1);
    assert_eq!(interpreted.report().records()[0].source_start_line(), 1);
    assert_eq!(
        interpreted.report().records()[0]
            .molfile()
            .atom_mappings()
            .len(),
        1
    );
    Ok(())
}

#[test]
fn hydrogen_normalization_public_api() -> Result<(), Box<dyn std::error::Error>> {
    let mut molecule = SmallMolecule::from_smiles_sanitized("C")?;
    let added = molecules::hydrogens::add_hydrogens(&mut molecule)?;
    assert_eq!(added.added.len(), 4);

    molecule.sanitize()?;
    let removed = molecule.remove_hydrogens()?;
    assert_eq!(removed.removed.len(), 4);
    assert_eq!(molecule.atom_count(), 1);
    Ok(())
}

#[test]
fn query_graph_smarts_and_substructure_public_api() -> Result<(), Box<dyn std::error::Error>> {
    let target = SmallMolecule::from_smiles_sanitized("CC(=O)O")?;
    let query = molecules::query::parse_smarts("[C](=O)[O;H1]")?;
    let matches = molecules::substructure::find_substructure_matches(target.graph(), &query)?;

    assert_eq!(query.atom_count(), 3);
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].atoms().len(), 3);
    Ok(())
}

#[test]
fn low_level_graph_api() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::core::*;

    let mut graph = Molecule::new();
    let carbon = graph.add_atom(Atom::new(
        Element::from_symbol("C").expect("carbon is a known element"),
    ));
    let oxygen = graph.add_atom(Atom::new(
        Element::from_symbol("O").expect("oxygen is a known element"),
    ));

    let bond = graph.add_bond(carbon, oxygen, BondOrder::Double)?;

    assert_eq!(graph.atom_count(), 2);
    assert_eq!(graph.bond_count(), 1);
    assert_eq!(graph.bond_between(carbon, oxygen)?, Some(bond));
    Ok(())
}

#[test]
fn macro_molecule_public_api() -> Result<(), Box<dyn std::error::Error>> {
    let mut builder = MacroMolecule::builder();
    let atom = builder.graph_mut().add_atom(Atom::new(
        Element::from_symbol("C").expect("carbon is a known element"),
    ));

    let model = builder.hierarchy_mut().add_model("1");
    let chain = builder.hierarchy_mut().add_chain(model, "A", None)?;
    let residue =
        builder
            .hierarchy_mut()
            .add_residue(chain, "GLY", Some(1), Some("1".to_owned()), None)?;
    builder.add_atom_site(
        residue,
        atom,
        molecules::bio::SmcraAtomSiteMetadata::default(),
    )?;
    let macro_mol = builder.build()?;

    let validate = macro_mol.validate()?;
    assert_eq!(validate.models_checked, 1);
    assert_eq!(validate.atom_sites_checked, 1);
    Ok(())
}

#[test]
fn model_and_smcra_model_names_coexist() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::bio::{SmcraHierarchy, SmcraModel};
    use molecules::modeling::Model;

    let mut hierarchy = SmcraHierarchy::new();
    let hierarchy_model_id = hierarchy.add_model("1");
    let hierarchy_model: &SmcraModel = hierarchy.model(hierarchy_model_id)?;
    assert_eq!(hierarchy_model.model_id(), "1");

    let mut graph = Molecule::new();
    let atom = graph.add_atom(Atom::new(
        Element::from_symbol("C").expect("carbon is a known element"),
    ));
    let mut conformer = Conformer::new(molecules::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            atom,
            molecules::units::Quantity::new(
                molecules::core::Point3::new(0.0, 0.0, 0.0),
                molecules::units::ANGSTROM,
            ),
        )
        .unwrap();
    let conformer = graph.add_conformer(conformer)?;
    let model = Model::from_small_molecule(&SmallMolecule::from_graph(graph), conformer)?;

    assert_eq!(model.atom_count(), 1);
    Ok(())
}

#[test]
fn small_molecule_modeling_public_api() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::modeling::potential::{HarmonicBondParameter, HarmonicBondPotential};
    use molecules::modeling::{minimize, MinimizationStatus, MinimizeOptions, Model};

    let mut graph = Molecule::new();
    let carbon = graph.add_atom(Atom::new(
        Element::from_symbol("C").expect("carbon is a known element"),
    ));
    let oxygen = graph.add_atom(Atom::new(
        Element::from_symbol("O").expect("oxygen is a known element"),
    ));
    let bond = graph.add_bond(carbon, oxygen, BondOrder::Single)?;
    let mut conformer = Conformer::new(molecules::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            carbon,
            molecules::units::Quantity::new(
                molecules::core::Point3::new(0.0, 0.0, 0.0),
                molecules::units::ANGSTROM,
            ),
        )
        .unwrap();
    conformer
        .set_position(
            oxygen,
            molecules::units::Quantity::new(
                molecules::core::Point3::new(2.0, 0.0, 0.0),
                molecules::units::ANGSTROM,
            ),
        )
        .unwrap();
    let conformer = graph.add_conformer(conformer).unwrap();
    let molecule = SmallMolecule::from_graph(graph);

    let mut builder = Model::builder();
    let instance = builder.add_small_molecule(&molecule, conformer)?;
    let model = builder.build()?;
    let cloned = model.clone();
    assert_eq!(model.definition_key(), cloned.definition_key());
    let model_bond = molecules::modeling::InstanceBondId::new(instance, bond);
    let mut potential = HarmonicBondPotential::new(
        &model,
        [HarmonicBondParameter::new(
            model_bond,
            molecules::units::Quantity::new(1.2, molecules::units::ANGSTROM),
            molecules::units::Quantity::new(100.0, molecules::units::MODEL_FORCE_CONSTANT_UNIT),
        )],
    )?;
    let result = minimize(&model, &mut potential, MinimizeOptions::default())?;

    assert_eq!(result.status, MinimizationStatus::Converged);
    assert!(result.final_energy < result.initial_energy);
    assert_eq!(model.positions()[1].x, 2.0);
    Ok(())
}

#[test]
fn production_smiles_stereo_uses_installed_perception_state(
) -> Result<(), Box<dyn std::error::Error>> {
    use molecules::perception::{self, stereo, SanitizeOptions};

    let document = molecules::smiles::parse_str(r"C(=C\F)\F")?;
    let mut molecule = molecules::smiles::interpret(&document)?.into_molecule();
    perception::sanitize_with_options(
        &mut molecule,
        SanitizeOptions {
            perceive_stereo: false,
            ..SanitizeOptions::default()
        },
    )?;

    let graph = molecule.graph();
    assert_eq!(graph.implicit_hydrogens(AtomId::new(0))?, Some(1));
    assert_eq!(graph.implicit_hydrogens(AtomId::new(1))?, Some(1));

    let report = stereo::perceive_stereo(molecule.graph_mut());
    assert!(report.is_ok(), "{:?}", report.issues);
    assert!(
        report.candidates.iter().any(|candidate| matches!(
            candidate,
            molecules::perception::stereo::StereoCandidate::DoubleBond {
                left_carriers,
                right_carriers,
                ..
            } if left_carriers.len() == 2 && right_carriers.len() == 2
        )),
        "{:?}",
        report.candidates
    );
    Ok(())
}

#[test]
fn production_atrop_cip_matches_pinned_reference() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::core::{StereoDescriptor, StereoElementId};
    use molecules::perception::{self, stereo};

    let input = include_str!(
        "../../../validation/corpora/smoke/data/rdkit_atropisomers/RP-6306_atrop4.mol"
    );
    let document = molecules::molfile::parse_str(input)?;
    let mut molecule = molecules::molfile::interpret(&document)?.into_molecule();
    perception::sanitize(&mut molecule)?;
    let report = stereo::assign_cip_descriptors(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(
        molecule.graph().cip_descriptor(StereoElementId::new(0))?,
        Some(StereoDescriptor::P)
    );
    Ok(())
}

#[test]
fn production_canonical_smiles_preserves_collapsed_hydrogen_without_perception(
) -> Result<(), Box<dyn std::error::Error>> {
    let document = molecules::smiles::parse_str("[H][C](F)(Cl)Br")?;
    let molecule = molecules::smiles::interpret(&document)?.into_molecule();
    assert!(!molecule.graph().perception().has_valence());

    let written = molecules::smiles::write_canonical(&molecule)?;
    let mut reparsed = SmallMolecule::from_smiles(&written)?;
    reparsed.sanitize()?;
    let carbon = reparsed
        .graph()
        .atoms()
        .find_map(|(atom_id, atom)| (atom.element.symbol() == "C").then_some(atom_id))
        .expect("canonical output retains carbon");
    assert_eq!(
        reparsed.graph().implicit_hydrogens(carbon)?,
        Some(1),
        "canonical output was {written}"
    );
    Ok(())
}
