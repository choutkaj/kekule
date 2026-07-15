use molecules::prelude::*;

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
    let mut mol = molecules::smiles::interpret(&document)?;
    molecules::perception::sanitize(&mut mol)?;
    let smiles = molecules::smiles::write_canonical(&mol)?;
    assert!(!smiles.is_empty());
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
    let mut macro_mol = MacroMolecule::new();
    let atom = macro_mol.graph_mut().add_atom(Atom::new(
        Element::from_symbol("C").expect("carbon is a known element"),
    ));

    let model = macro_mol.hierarchy_mut().add_model("1");
    let chain = macro_mol.hierarchy_mut().add_chain(model, "A", None)?;
    let residue =
        macro_mol
            .hierarchy_mut()
            .add_residue(chain, "GLY", Some(1), Some("1".to_owned()), None)?;
    macro_mol.add_atom_site(residue, atom, molecules::bio::AtomSiteMetadata::default())?;

    let validate = macro_mol.validate()?;
    assert_eq!(validate.models_checked, 1);
    assert_eq!(validate.atom_sites_checked, 1);

    let sanitize = macro_mol.sanitize()?;
    assert!(sanitize.validation.is_some());
    Ok(())
}

#[test]
fn small_molecule_modeling_public_api() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::modeling::potential::{HarmonicBondParameter, HarmonicBondPotential};
    use molecules::modeling::{minimize, MinimizationStatus, MinimizeOptions, MolecularModel};

    let mut graph = Molecule::new();
    let carbon = graph.add_atom(Atom::new(
        Element::from_symbol("C").expect("carbon is a known element"),
    ));
    let oxygen = graph.add_atom(Atom::new(
        Element::from_symbol("O").expect("oxygen is a known element"),
    ));
    let bond = graph.add_bond(carbon, oxygen, BondOrder::Single)?;
    let mut conformer = Conformer::new();
    conformer.set_position(carbon, molecules::core::Point3::new(0.0, 0.0, 0.0));
    conformer.set_position(oxygen, molecules::core::Point3::new(2.0, 0.0, 0.0));
    let conformer = graph.add_conformer(conformer).unwrap();
    let molecule = SmallMolecule::from_graph(graph);

    let mut builder = MolecularModel::builder();
    let instance = builder.add_small_molecule(&molecule, conformer)?;
    let model = builder.build()?;
    let cloned = model.clone();
    assert_eq!(model.definition_key(), cloned.definition_key());
    let model_bond = molecules::modeling::InstanceBondId::new(instance, bond);
    let mut potential =
        HarmonicBondPotential::new(&model, [HarmonicBondParameter::new(model_bond, 1.2, 100.0)])?;
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
    let mut molecule = molecules::smiles::interpret(&document)?;
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
    let mut molecule = molecules::molfile::interpret(&document)?;
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
    let molecule = molecules::smiles::interpret(&document)?;
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
