use kekule::prelude::*;

#[test]
fn quantity_and_unit_public_api() -> Result<(), Box<dyn std::error::Error>> {
    use kekule::units::{Dimension, Quantity, Unit, ANGSTROM, NANOMETER};

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
    let formal_charge: i64 = mol.graph().formal_charge();
    assert_eq!(formal_charge, 0);
    let smiles = mol.to_canonical_smiles()?;
    assert!(!smiles.is_empty());
    Ok(())
}

#[test]
fn molecular_descriptor_public_api() -> Result<(), Box<dyn std::error::Error>> {
    use kekule::descriptors::{
        average_mass, molecular_formula, monoisotopic_mass, HydrogenCountPolicy,
        MolecularDescriptorError, MolecularFormula,
    };
    use kekule::units::DALTON;

    let molecule = SmallMolecule::from_smiles_sanitized("[13CH3]CO")?;
    let formula: MolecularFormula =
        molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived)?;
    assert_eq!(formula.to_string(), "C[13C]H6O");
    assert_eq!(
        formula.isotope_count(Element::from_symbol("C").unwrap(), 13),
        1
    );
    assert_eq!(formula.formal_charge(), 0);

    let average = average_mass(&molecule, HydrogenCountPolicy::IncludePerceived)?;
    let monoisotopic = monoisotopic_mass(&molecule, HydrogenCountPolicy::IncludePerceived)?;
    assert!(average.value_in(DALTON)? > monoisotopic.value_in(DALTON)?);

    let raw = SmallMolecule::from_smiles("C")?;
    let error: MolecularDescriptorError =
        molecular_formula(&raw, HydrogenCountPolicy::IncludePerceived).unwrap_err();
    assert!(matches!(
        error,
        MolecularDescriptorError::MissingImplicitHydrogens { .. }
    ));
    Ok(())
}

#[test]
fn namespaced_small_molecule_api() -> Result<(), Box<dyn std::error::Error>> {
    let document = kekule::smiles::parse_str("CC(=O)O")?;
    let interpreted = kekule::smiles::interpret(&document)?;
    assert_eq!(interpreted.report().atom_mappings().len(), 4);
    assert_eq!(interpreted.report().bond_mappings().len(), 3);
    let mut mol = interpreted.into_molecule();
    kekule::perception::sanitize(&mut mol)?;
    let smiles = kekule::smiles::write_canonical(&mol)?;
    assert!(!smiles.is_empty());
    Ok(())
}

#[test]
fn molfile_and_sdf_interpretation_reports_are_public() -> Result<(), Box<dyn std::error::Error>> {
    let molfile = "\
Report
kekule

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
";
    let document = kekule::molfile::parse_str(molfile)?;
    let interpreted = kekule::molfile::interpret(&document)?;
    assert_eq!(interpreted.report().atom_mappings().len(), 1);
    assert_eq!(interpreted.report().atom_mappings()[0].source_line(), 5);

    let sdf = format!("{molfile}>  <SOURCE>\nrelease-test\n\n$$$$\n");
    let document = kekule::sdf::parse_str(&sdf, kekule::sdf::SdfParseOptions::default())?;
    let interpreted = kekule::sdf::interpret(&document)?;
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
fn parser_resource_options_are_public() -> Result<(), Box<dyn std::error::Error>> {
    let smiles = kekule::smiles::parse_str_with_options(
        "CC",
        kekule::smiles::SmilesParseOptions::default(),
    )?;
    assert_eq!(smiles.tokens().len(), 2);

    let molfile = "methane\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n";
    let document = kekule::molfile::parse_str_with_options(
        molfile,
        kekule::molfile::MolfileParseOptions::default(),
    )?;
    assert_eq!(document.atom_records().len(), 1);

    let sdf = format!("{molfile}$$$$\n");
    let document = kekule::sdf::parse_str(&sdf, kekule::sdf::SdfParseOptions::default())?;
    assert_eq!(document.records().len(), 1);
    Ok(())
}

#[test]
fn hydrogen_normalization_public_api() -> Result<(), Box<dyn std::error::Error>> {
    let mut molecule = SmallMolecule::from_smiles_sanitized("C")?;
    let added = kekule::hydrogens::add_hydrogens(&mut molecule)?;
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
    let query = kekule::query::parse_smarts("[C](=O)[O;H1]")?;
    let matches = kekule::substructure::find_substructure_matches(target.graph(), &query)?;

    assert_eq!(query.atom_count(), 3);
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].atoms().len(), 3);
    Ok(())
}

#[test]
fn low_level_graph_api() -> Result<(), Box<dyn std::error::Error>> {
    use kekule::core::*;

    let mut graph = Molecule::new();
    let carbon = graph.add_atom(Atom::new(
        Element::from_symbol("C").expect("carbon is a known element"),
    ))?;
    let oxygen = graph.add_atom(Atom::new(
        Element::from_symbol("O").expect("oxygen is a known element"),
    ))?;

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
    ))?;

    let chain = builder.hierarchy_mut().add_chain("A", None)?;
    let residue =
        builder
            .hierarchy_mut()
            .add_residue(chain, "GLY", Some(1), Some("1".to_owned()), None)?;
    builder.add_atom_site(residue, atom, kekule::bio::SmcraAtomSiteMetadata::default())?;
    let macro_mol = builder.build()?;

    let validate = macro_mol.validate()?;
    assert_eq!(validate.chains_checked, 1);
    assert_eq!(validate.atom_sites_checked, 1);
    Ok(())
}

#[test]
fn model_and_static_smcra_hierarchy_coexist() -> Result<(), Box<dyn std::error::Error>> {
    use kekule::bio::SmcraHierarchy;
    use kekule::structure::Model;

    let mut hierarchy = SmcraHierarchy::new();
    let chain = hierarchy.add_chain("A", None)?;
    assert_eq!(hierarchy.chain(chain)?.label_id(), "A");

    let mut graph = Molecule::new();
    let atom = graph.add_atom(Atom::new(
        Element::from_symbol("C").expect("carbon is a known element"),
    ))?;
    let mut conformer = Conformer::new(kekule::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            atom,
            kekule::units::Quantity::new(
                kekule::geometry::Point3::new(0.0, 0.0, 0.0),
                kekule::units::ANGSTROM,
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
    use kekule::geometry::{PeriodicCell, Vector3};
    use kekule::modeling::potential::{
        HarmonicBondParameter, HarmonicBondPotential, Potential, PotentialError,
    };
    use kekule::modeling::{minimize, MinimizationStatus, MinimizeOptions};
    use kekule::structure::Model;
    use kekule::topology::InstanceBondId;

    let mut graph = Molecule::new();
    let carbon = graph.add_atom(Atom::new(
        Element::from_symbol("C").expect("carbon is a known element"),
    ))?;
    let oxygen = graph.add_atom(Atom::new(
        Element::from_symbol("O").expect("oxygen is a known element"),
    ))?;
    let bond = graph.add_bond(carbon, oxygen, BondOrder::Single)?;
    let mut conformer = Conformer::new(kekule::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            carbon,
            kekule::units::Quantity::new(
                kekule::geometry::Point3::new(0.0, 0.0, 0.0),
                kekule::units::ANGSTROM,
            ),
        )
        .unwrap();
    conformer
        .set_position(
            oxygen,
            kekule::units::Quantity::new(
                kekule::geometry::Point3::new(2.0, 0.0, 0.0),
                kekule::units::ANGSTROM,
            ),
        )
        .unwrap();
    let conformer = graph.add_conformer(conformer).unwrap();
    let mut molecule = SmallMolecule::from_graph(graph);

    let mut builder = Model::builder();
    let instance = builder.add_small_molecule(&molecule, conformer)?;
    let model = builder.build()?;
    let cloned = model.clone();
    assert!(model.topology().same_identity(cloned.topology()));
    let model_bond = InstanceBondId::new(instance, bond);
    let mut potential = HarmonicBondPotential::new(
        model.topology(),
        [HarmonicBondParameter::new(
            model_bond,
            kekule::units::Quantity::new(1.2, kekule::units::ANGSTROM),
            kekule::units::Quantity::new(100.0, kekule::units::MODEL_FORCE_CONSTANT_UNIT),
        )],
    )?;
    let result = minimize(&model, &mut potential, MinimizeOptions::default())?;
    let mut periodic = model.clone();
    periodic.set_cell(Some(PeriodicCell::orthorhombic(
        kekule::units::Quantity::new(Vector3::new(10.0, 10.0, 10.0), kekule::units::ANGSTROM),
        [true; 3],
    )?));
    assert_eq!(
        potential.evaluate(periodic.view()),
        Err(PotentialError::UnsupportedPeriodicCell)
    );

    result
        .model
        .instance_to_conformer(instance, molecule.graph_mut(), conformer)?;

    assert_eq!(result.status, MinimizationStatus::Converged);
    assert!(result.final_energy < result.initial_energy);
    assert_eq!(model.positions()[1].x, 2.0);
    assert!(
        molecule
            .graph()
            .conformer(conformer)?
            .position(oxygen)
            .expect("oxygen position")
            .x
            < 2.0
    );
    Ok(())
}

#[test]
fn topology_and_ensemble_public_api() -> Result<(), Box<dyn std::error::Error>> {
    use kekule::geometry::Point3;
    use kekule::structure::{Configuration, Ensemble, EnsembleMember, Model, Positions};
    use kekule::topology::{
        AtomSelection, MoleculeInstanceMetadata, MoleculeRole, TopologyBuilder,
    };
    use kekule::units::{Quantity, ANGSTROM};

    let water = SmallMolecule::from_smiles_sanitized("O")?;
    let mut topology_builder = TopologyBuilder::new();
    let definition = topology_builder.add_small_molecule_definition(&water)?;
    let mut metadata = MoleculeInstanceMetadata::default();
    metadata.insert_role(MoleculeRole::Solvent);
    topology_builder.add_instance(definition, metadata.clone())?;
    topology_builder.add_instance(definition, metadata)?;
    let topology = topology_builder.build()?;

    let first_positions = Positions::new(
        &topology,
        Quantity::new(
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 0.0, 0.0)],
            ANGSTROM,
        ),
    )?;
    let second_positions = Positions::new(
        &topology,
        Quantity::new(
            vec![Point3::new(0.1, 0.0, 0.0), Point3::new(3.1, 0.0, 0.0)],
            ANGSTROM,
        ),
    )?;
    let first = Model::new(
        topology.clone(),
        Configuration::new(first_positions.clone()),
    )?;
    let second = Model::new(
        topology.clone(),
        Configuration::new(second_positions.clone()),
    )?;
    assert!(first.topology().same_identity(second.topology()));

    let solvent = AtomSelection::for_roles(&topology, [MoleculeRole::Solvent])?;
    assert_eq!(solvent.indices().len(), 2);
    let ensemble = Ensemble::from_members(
        topology.clone(),
        [
            EnsembleMember::new(Configuration::new(first_positions)),
            EnsembleMember::new(Configuration::new(second_positions.clone())),
        ],
    )?;
    assert_eq!(ensemble.views().count(), 2);

    Ok(())
}

#[test]
fn topology_layout_and_checked_mapping_public_api() -> Result<(), Box<dyn std::error::Error>> {
    use kekule::topology::{
        MoleculeInstanceMetadata, TopologyBuilder, TopologyEditResult, TopologyMapping,
    };

    let water = SmallMolecule::from_smiles_sanitized("O")?;
    let build = || -> Result<_, Box<dyn std::error::Error>> {
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&water)?;
        builder.add_instance(definition, MoleculeInstanceMetadata::default())?;
        Ok(builder.build()?)
    };
    let source = build()?;
    let target = build()?;

    assert!(!source.same_identity(&target));
    assert!(source.same_layout(&target));
    let mapping = TopologyMapping::between_identical_layouts(&source, &target)?;
    let result = TopologyEditResult::new(target.clone(), mapping)?;
    assert!(result.topology().same_identity(&target));
    Ok(())
}

#[test]
fn production_smiles_stereo_uses_installed_perception_state(
) -> Result<(), Box<dyn std::error::Error>> {
    use kekule::perception::{self, stereo, SanitizeOptions};

    let document = kekule::smiles::parse_str(r"C(=C\F)\F")?;
    let mut molecule = kekule::smiles::interpret(&document)?.into_molecule();
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
            kekule::perception::stereo::StereoCandidate::DoubleBond {
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
    use kekule::core::{StereoDescriptor, StereoElementId};
    use kekule::perception::{self, stereo};

    let input = include_str!(
        "../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/RP-6306_atrop4.mol"
    );
    let document = kekule::molfile::parse_str(input)?;
    let mut molecule = kekule::molfile::interpret(&document)?.into_molecule();
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
    let document = kekule::smiles::parse_str("[H][C](F)(Cl)Br")?;
    let molecule = kekule::smiles::interpret(&document)?.into_molecule();
    assert!(!molecule.graph().perception().has_valence());

    let written = kekule::smiles::write_canonical(&molecule)?;
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
