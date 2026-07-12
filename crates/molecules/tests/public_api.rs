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
    let conformer = graph.add_conformer(conformer);
    let molecule = SmallMolecule::from_graph(graph);

    let mut builder = MolecularModel::builder();
    let instance = builder.add_small_molecule(&molecule, conformer)?;
    let model = builder.build()?;
    let model_bond = molecules::modeling::InstanceBondId::new(instance, bond);
    let mut potential =
        HarmonicBondPotential::new(&model, [HarmonicBondParameter::new(model_bond, 1.2, 100.0)])?;
    let result = minimize(&model, &mut potential, MinimizeOptions::default())?;

    assert_eq!(result.status, MinimizationStatus::Converged);
    assert!(result.final_energy < result.initial_energy);
    assert_eq!(model.positions()[1].x, 2.0);
    Ok(())
}
