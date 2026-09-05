use std::sync::Arc;

use kekule::core::{Atom, BondOrder, Element, HydrogenDeclaration, Molecule, MoleculeEditor};
use kekule::structure::{Model, Positions};
use kekule::topology::{
    AtomSelection, AtomSiteMetadata, InstanceAtomId, MoleculeClass, ResidueClass, Topology,
    TopologyBuilder,
};

fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}

fn linear_residue_topology(
    residue_names: &[&str],
    atom_names: (&str, &str),
) -> (Topology, Vec<[kekule::core::AtomId; 2]>) {
    let mut editor = MoleculeEditor::new();
    let mut residue_atoms: Vec<[kekule::core::AtomId; 2]> = Vec::new();
    for _ in residue_names {
        let first = editor.add_atom(atom("C")).unwrap();
        let second = editor.add_atom(atom("N")).unwrap();
        editor.add_bond(first, second, BondOrder::Single).unwrap();
        if let Some(previous) = residue_atoms.last().copied() {
            editor
                .add_bond(previous[1], first, BondOrder::Single)
                .unwrap();
        }
        residue_atoms.push([first, second]);
    }
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    for (index, (&name, atoms)) in residue_names.iter().zip(&residue_atoms).enumerate() {
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, name, Some(index as i32 + 1), None, None)
            .unwrap();
        for (atom, name) in atoms.iter().copied().zip([atom_names.0, atom_names.1]) {
            builder
                .hierarchy_mut()
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, atom),
                    AtomSiteMetadata {
                        label_atom_id: Some(name.to_owned()),
                        ..AtomSiteMetadata::default()
                    },
                )
                .unwrap();
        }
    }
    (builder.build().unwrap(), residue_atoms)
}

fn one_residue_topology(
    molecule: &Molecule,
    component: &str,
) -> (Topology, kekule::topology::ResidueId) {
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(molecule).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, component, Some(1), None, None)
        .unwrap();
    for (index, atom) in molecule.atom_ids().enumerate() {
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, atom),
                AtomSiteMetadata {
                    label_atom_id: Some(format!("X{}", index + 1)),
                    ..AtomSiteMetadata::default()
                },
            )
            .unwrap();
    }
    (builder.build().unwrap(), residue)
}

#[test]
fn ordinary_connected_molecule_defaults_to_small_molecule() {
    let molecule = kekule::smiles::to_molecules("CCO").unwrap().pop().unwrap();
    let (topology, residue) = one_residue_topology(&molecule, "UNL");
    let instance = topology.molecules().next().unwrap();
    assert_eq!(instance.class(), MoleculeClass::SmallMolecule);
    assert_eq!(
        topology.molecule_class(instance.id()).unwrap(),
        MoleculeClass::SmallMolecule
    );
    assert_eq!(
        topology.residue(residue).unwrap().class(),
        ResidueClass::Other
    );
}

#[test]
fn water_and_monoatomic_ions_use_strong_local_evidence() {
    let mut water_editor = MoleculeEditor::new();
    let mut oxygen = atom("O");
    oxygen.hydrogens = HydrogenDeclaration::Fixed(2);
    water_editor.add_atom(oxygen).unwrap();
    let explicit_water = Topology::from_molecule(&water_editor.finish().unwrap()).unwrap();
    assert_eq!(
        explicit_water.molecules().next().unwrap().class(),
        MoleculeClass::Water
    );

    let mut oxygen_only = MoleculeEditor::new();
    oxygen_only.add_atom(atom("O")).unwrap();
    let (oxygen_only, residue) = one_residue_topology(&oxygen_only.finish().unwrap(), "HOH");
    assert_eq!(
        oxygen_only.residue(residue).unwrap().class(),
        ResidueClass::Water
    );
    assert_eq!(
        oxygen_only.molecules().next().unwrap().class(),
        MoleculeClass::Water
    );

    let mut ion_editor = MoleculeEditor::new();
    let mut sodium = atom("Na");
    sodium.formal_charge = 1;
    ion_editor.add_atom(sodium).unwrap();
    let ion = Topology::from_molecule(&ion_editor.finish().unwrap()).unwrap();
    assert_eq!(ion.molecules().next().unwrap().class(), MoleculeClass::Ion);
}

#[test]
fn peptide_and_modified_peptide_classify_as_protein() {
    let (peptide, _) = linear_residue_topology(&["ALA", "GLY"], ("N", "C"));
    assert_eq!(
        peptide.molecules().next().unwrap().class(),
        MoleculeClass::Protein
    );
    assert!(peptide
        .residues()
        .all(|residue| residue.class() == ResidueClass::AminoAcid));

    let (modified, _) = linear_residue_topology(&["ALA", "MSE", "GLY"], ("N", "C"));
    assert_eq!(
        modified.molecules().next().unwrap().class(),
        MoleculeClass::Protein
    );
    assert_eq!(
        modified.residues().nth(1).unwrap().class(),
        ResidueClass::Other
    );
}

#[test]
fn dna_rna_and_modified_backbones_classify_from_phosphodiester_links() {
    let (dna, _) = linear_residue_topology(&["DA", "DC"], ("P", "O3'"));
    assert_eq!(dna.molecules().next().unwrap().class(), MoleculeClass::Dna);
    let (rna, _) = linear_residue_topology(&["A", "U"], ("P", "O3'"));
    assert_eq!(rna.molecules().next().unwrap().class(), MoleculeClass::Rna);

    let (modified_dna, _) = linear_residue_topology(&["DA", "5MC", "DG"], ("P", "O3'"));
    assert_eq!(
        modified_dna.molecules().next().unwrap().class(),
        MoleculeClass::Dna
    );
    let (modified_rna, _) = linear_residue_topology(&["A", "PSU", "G"], ("P", "O3'"));
    assert_eq!(
        modified_rna.molecules().next().unwrap().class(),
        MoleculeClass::Rna
    );

    let (hybrid, _) = linear_residue_topology(&["DA", "A"], ("P", "O3'"));
    assert_eq!(
        hybrid.molecules().next().unwrap().class(),
        MoleculeClass::Other
    );
}

#[test]
fn recognized_carbohydrate_component_classifies_its_molecule() {
    let mut editor = MoleculeEditor::new();
    editor.add_atom(atom("C")).unwrap();
    let (topology, residue) = one_residue_topology(&editor.finish().unwrap(), "GLC");
    assert_eq!(
        topology.residue(residue).unwrap().class(),
        ResidueClass::Carbohydrate
    );
    assert_eq!(
        topology.molecules().next().unwrap().class(),
        MoleculeClass::Carbohydrate
    );
}

#[test]
fn explicit_definition_and_residue_overrides_win() {
    let molecule = kekule::smiles::to_molecules("CC").unwrap().pop().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    let instance = builder.add_instance(definition).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "UNL", Some(1), None, None)
        .unwrap();
    for atom in molecule.atom_ids() {
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, atom),
                AtomSiteMetadata::default(),
            )
            .unwrap();
    }
    builder
        .set_molecule_class(definition, MoleculeClass::Protein)
        .unwrap();
    builder
        .set_residue_class(residue, ResidueClass::AminoAcid)
        .unwrap();
    let topology = builder.build().unwrap();
    assert_eq!(
        topology.definition(definition).unwrap().class(),
        MoleculeClass::Protein
    );
    assert_eq!(
        topology.residue(residue).unwrap().class(),
        ResidueClass::AminoAcid
    );
}

#[test]
fn reused_definition_shares_class_and_conflicting_strong_instances_become_other() {
    let (template, _) = linear_residue_topology(&["ALA", "GLY"], ("N", "C"));
    let molecule = template.molecules().next().unwrap().molecule().clone();
    let atoms = molecule.atom_ids().collect::<Vec<_>>();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    let peptide = builder.add_instance(definition).unwrap();
    let dna = builder.add_instance(definition).unwrap();
    for (chain_name, instance, residues, names) in [
        ("P", peptide, ["ALA", "GLY"], ["N", "C"]),
        ("D", dna, ["DA", "DC"], ["P", "O3'"]),
    ] {
        let chain = builder.hierarchy_mut().add_chain(chain_name, None).unwrap();
        for index in 0..2 {
            let residue = builder
                .hierarchy_mut()
                .add_residue(chain, residues[index], Some(index as i32 + 1), None, None)
                .unwrap();
            for offset in 0..2 {
                builder
                    .hierarchy_mut()
                    .add_atom_site(
                        residue,
                        InstanceAtomId::new(instance, atoms[index * 2 + offset]),
                        AtomSiteMetadata {
                            label_atom_id: Some(names[offset].to_owned()),
                            ..AtomSiteMetadata::default()
                        },
                    )
                    .unwrap();
            }
        }
    }
    let topology = builder.build().unwrap();
    assert_eq!(
        topology.definition(definition).unwrap().class(),
        MoleculeClass::Other
    );
    assert!(topology
        .molecules()
        .all(|instance| instance.class() == MoleculeClass::Other));
}

#[test]
fn append_preserves_class_and_structural_subset_reinfers() {
    let (protein, residue_atoms) = linear_residue_topology(&["ALA", "GLY"], ("N", "C"));
    let source = Arc::new(protein);
    assert_eq!(
        source.molecules().next().unwrap().class(),
        MoleculeClass::Protein
    );
    let instance = source.molecules().next().unwrap().id();
    let first_residue = AtomSelection::from_atoms(
        &source,
        residue_atoms[0]
            .iter()
            .copied()
            .map(|atom| InstanceAtomId::new(instance, atom)),
    )
    .unwrap();
    let subset = source.subset(&first_residue).unwrap();
    assert_eq!(
        subset.topology().molecules().next().unwrap().class(),
        MoleculeClass::SmallMolecule
    );
    drop(subset);
    drop(first_residue);

    let existing_class = source.molecules().next().unwrap().class();
    let existing_residue_classes = source
        .residues()
        .map(|residue| residue.class())
        .collect::<Vec<_>>();
    let mut builder = Arc::try_unwrap(source).unwrap().into_builder();
    let ligand = kekule::smiles::to_molecules("CCO").unwrap().pop().unwrap();
    builder.add_molecule(&ligand).unwrap();
    let appended = builder.build().unwrap();
    assert_eq!(appended.molecules().next().unwrap().class(), existing_class);
    assert_eq!(
        appended
            .residues()
            .map(|residue| residue.class())
            .collect::<Vec<_>>(),
        existing_residue_classes
    );
}

#[test]
fn subsets_reuse_compact_definitions_and_preserve_only_complete_entity_classes() {
    let mut editor = MoleculeEditor::new();
    let removed = editor.add_atom(atom("C")).unwrap();
    let first = editor.add_atom(atom("C")).unwrap();
    let second = editor.add_atom(atom("C")).unwrap();
    editor.add_bond(removed, first, BondOrder::Single).unwrap();
    let bond = editor.add_bond(first, second, BondOrder::Single).unwrap();
    editor.delete_atom(removed).unwrap();
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder
        .set_molecule_class(definition, MoleculeClass::Other)
        .unwrap();
    let mut instances = Vec::new();
    for name in ["A", "B", "C"] {
        let instance = builder.add_instance(definition).unwrap();
        instances.push(instance);
        let chain = builder.hierarchy_mut().add_chain(name, None).unwrap();
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, "UNL", None, None, None)
            .unwrap();
        for atom in [first, second] {
            builder
                .hierarchy_mut()
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, atom),
                    AtomSiteMetadata::default(),
                )
                .unwrap();
        }
        builder
            .set_residue_class(residue, ResidueClass::Water)
            .unwrap();
    }
    let source = Arc::new(builder.build().unwrap());
    let selection = AtomSelection::from_atoms(
        &source,
        instances.iter().enumerate().flat_map(|(index, &instance)| {
            [first, second]
                .into_iter()
                .take(if index == 2 { 1 } else { 2 })
                .map(move |atom| InstanceAtomId::new(instance, atom))
        }),
    )
    .unwrap();
    let subset = source.subset(&selection).unwrap();
    let target = subset.topology();
    assert_eq!(target.definition_count(), 2);
    let molecules = target.molecules().collect::<Vec<_>>();
    assert_eq!(molecules[0].definition_id(), molecules[1].definition_id());
    assert_ne!(molecules[0].definition_id(), molecules[2].definition_id());
    assert_eq!(
        molecules
            .iter()
            .map(|molecule| molecule.class())
            .collect::<Vec<_>>(),
        [
            MoleculeClass::Other,
            MoleculeClass::Other,
            MoleculeClass::SmallMolecule
        ]
    );
    assert_eq!(
        target
            .residues()
            .map(|residue| residue.class())
            .collect::<Vec<_>>(),
        [
            ResidueClass::Water,
            ResidueClass::Water,
            ResidueClass::Other
        ]
    );
    for &instance in &instances {
        assert_eq!(
            subset
                .correspondence()
                .target_atom(InstanceAtomId::new(instance, first))
                .unwrap()
                .atom()
                .raw(),
            0
        );
    }
    for &instance in &instances[..2] {
        assert_eq!(
            subset
                .correspondence()
                .target_atom(InstanceAtomId::new(instance, second))
                .unwrap()
                .atom()
                .raw(),
            1
        );
        assert_eq!(
            subset
                .correspondence()
                .target_bond(kekule::topology::InstanceBondId::new(instance, bond))
                .unwrap()
                .bond()
                .raw(),
            0
        );
    }
}

#[test]
fn typed_class_selections_select_instances_and_residues() {
    let (protein, _) = linear_residue_topology(&["ALA", "GLY"], ("N", "C"));
    let topology = Arc::new(protein);
    assert_eq!(
        AtomSelection::for_molecule_classes(&topology, [MoleculeClass::Protein])
            .unwrap()
            .indices()
            .len(),
        topology.atom_count()
    );
    assert_eq!(
        AtomSelection::for_residue_classes(&topology, [ResidueClass::AminoAcid])
            .unwrap()
            .indices()
            .len(),
        topology.atom_count()
    );
}

#[test]
fn generic_mmcif_writing_maps_canonical_protein_water_and_small_molecule_classes() {
    let molecule = kekule::smiles::to_molecules("CC").unwrap().pop().unwrap();
    let model = Model::from_molecule(&molecule, &Positions::zeros(molecule.atom_count())).unwrap();
    let written = kekule::mmcif::write_model(&model, Default::default()).unwrap();
    assert!(written.contains("1 non-polymer"));

    let (protein, _) = linear_residue_topology(&["ALA", "GLY"], ("N", "C"));
    let protein = Model::new(protein, Positions::zeros(4)).unwrap();
    let written = kekule::mmcif::write_model(&protein, Default::default()).unwrap();
    assert!(written.contains("1 polymer"));

    let mut water_editor = MoleculeEditor::new();
    water_editor.add_atom(atom("O")).unwrap();
    let (water, _) = one_residue_topology(&water_editor.finish().unwrap(), "HOH");
    let water = Model::new(water, Positions::zeros(1)).unwrap();
    let written = kekule::mmcif::write_model(&water, Default::default()).unwrap();
    assert!(written.contains("1 water"));
}

#[test]
fn readme_combined_model_workflow_needs_no_mmcif_classification_sidecar() {
    let first = kekule::smiles::to_molecules("CC").unwrap().pop().unwrap();
    let second = kekule::smiles::to_molecules("CO").unwrap().pop().unwrap();
    let mut builder = Model::builder();
    let first_id = builder
        .add_molecule(&first, &Positions::zeros(first.atom_count()))
        .unwrap();
    let second_id = builder
        .add_molecule(&second, &Positions::zeros(second.atom_count()))
        .unwrap();
    let combined = builder.build().unwrap();
    assert_eq!(
        combined.topology().molecule_class(first_id).unwrap(),
        MoleculeClass::SmallMolecule
    );
    assert_eq!(
        combined.topology().molecule_class(second_id).unwrap(),
        MoleculeClass::SmallMolecule
    );
    let mut output = Vec::new();
    kekule::mmcif::write_model_to(&mut output, &combined, Default::default()).unwrap();
    assert_eq!(
        String::from_utf8(output)
            .unwrap()
            .matches(" non-polymer")
            .count(),
        2
    );
}
