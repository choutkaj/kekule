use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor};
use kekule::geometry::Point3;
use kekule::sdf;
use kekule::smiles;
use kekule::structure::{Model, Positions};
use kekule::units::{Quantity, ANGSTROM};

fn one_smiles(input: &str) -> Molecule {
    let mut molecules = smiles::to_molecules(input).expect("SMILES interprets");
    assert_eq!(molecules.len(), 1);
    molecules.pop().expect("component count was checked")
}

#[test]
fn smiles_partitions_components_in_source_order() {
    let molecules = smiles::to_molecules("CC.O.[Na+]").expect("dot SMILES interprets");
    assert_eq!(molecules.len(), 3);
    assert_eq!(molecules[0].atom_count(), 2);
    assert_eq!(molecules[1].atoms().next().unwrap().1.element.symbol(), "O");
    assert_eq!(
        molecules[2].atoms().next().unwrap().1.element.symbol(),
        "Na"
    );
}

#[test]
fn smiles_convenience_preserves_cardinality_and_matches_explicit_pipeline() {
    let ethanol = smiles::to_molecules("CCO").expect("ethanol interprets");
    assert_eq!(ethanol.len(), 1);

    let salt = smiles::to_molecules("[Na+].[Cl-]").expect("salt interprets");
    assert_eq!(salt.len(), 2);

    let document = smiles::parse_str("[Na+].[Cl-]").expect("salt parses");
    let explicit = smiles::interpret(&document)
        .expect("salt interprets explicitly")
        .to_molecules();
    assert_eq!(salt, explicit);
}

#[test]
fn component_interpretation_types_are_nameable_through_format_facades() {
    let smiles_document = smiles::parse_str("CC.O").expect("SMILES parses");
    let smiles_interpretation = smiles_document.interpret().expect("SMILES interprets");
    let smiles_components: &[smiles::SmilesComponentInterpretation] =
        smiles_interpretation.components();
    assert_eq!(smiles_components.len(), 2);
    let count_error: smiles::SmilesComponentCountError = smiles_interpretation
        .molecule()
        .expect_err("singular access rejects multiple components");
    assert_eq!(count_error.actual(), 2);

    let molfile_source = "two components\nkekule\n\n  2  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n    1.0000    0.0000    0.0000 O   0  0  0  0  0  0\nM  END\n";
    let molfile_document = kekule::molfile::parse_str(molfile_source).expect("Molfile parses");
    let molfile_interpretation = molfile_document.interpret().expect("Molfile interprets");
    let molfile_components: &[kekule::molfile::MolfileComponentInterpretation] =
        molfile_interpretation.components();
    assert_eq!(molfile_components.len(), 2);
}

#[test]
fn sdf_interpretation_error_kind_is_public_and_structured() {
    let source = "unknown element\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 Xx  0  0  0  0  0  0\nM  END\n$$$$\n";
    let document = sdf::parse_str(source).expect("SDF syntax parses");
    let error = document.records()[0]
        .interpret()
        .expect_err("unsupported element fails Molfile interpretation");
    let kind: &sdf::SdfInterpretErrorKind = error.kind();

    match kind {
        sdf::SdfInterpretErrorKind::Molfile(source) => {
            assert!(source.message().contains("unsupported element"));
        }
        sdf::SdfInterpretErrorKind::ModelBuild(_) => {
            panic!("fixture should fail before model construction");
        }
        _ => panic!("unexpected future SDF interpretation error kind"),
    }
    assert_eq!(error.record(), 1);
    assert!(error.line() >= 1);
    assert!(error.message().contains("unsupported element"));
    assert!(std::error::Error::source(&error).is_some());
}

#[test]
fn smiles_default_options_and_topology_projection_are_consistent() {
    let default_document = smiles::parse_str("CCO.[Na+]").expect("default parse");
    let explicit_document =
        smiles::parse_str_with_options("CCO.[Na+]", smiles::SmilesParseOptions::default())
            .expect("explicit default parse");
    assert_eq!(default_document, explicit_document);

    let interpretation = default_document.interpret().expect("document interprets");
    assert_eq!(interpretation.molecules().count(), 2);
    assert!(interpretation
        .molecules()
        .all(|molecule| molecule.perception() == &kekule::core::Perception::default()));

    let topology = interpretation.to_topology().expect("topology builds");
    assert_eq!(topology.instance_count(), 2);
    assert_eq!(topology.atom_count(), 4);
    assert_eq!(
        topology
            .molecules()
            .map(|occurrence| occurrence.molecule().atom_count())
            .collect::<Vec<_>>(),
        vec![3, 1]
    );
    assert!(topology.hierarchy().is_empty());

    let concise = smiles::to_topology("CCO.[Na+]").expect("concise topology builds");
    assert!(concise.same_layout(&topology));
    assert_eq!(
        smiles::to_topology("CCO")
            .expect("connected topology")
            .instance_count(),
        1
    );
}

#[test]
fn smiles_writers_remain_available_through_the_format_namespace() {
    let ethanol = one_smiles("CCO");
    assert_eq!(smiles::write(&ethanol).expect("SMILES writes"), "CCO");
    assert_eq!(
        smiles::write_canonical(&ethanol).expect("canonical SMILES writes"),
        "CCO"
    );

    let chiral = one_smiles("F[C@H](Cl)Br");
    assert_eq!(
        smiles::write_isomeric(&chiral).expect("isomeric SMILES writes"),
        "F[C@H](Cl)Br"
    );
}

#[test]
fn sdf_preserves_record_boundaries_and_component_order() {
    let input = "first\nkekule\n\n  2  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\nsecond\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n";
    let document = sdf::parse_str(input).unwrap();
    let records = sdf::interpret(&document).unwrap().to_records();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].title(), "first");
    assert_eq!(records[0].molecules().len(), 2);
    assert_eq!(records[1].title(), "second");
    assert_eq!(records[1].molecules().len(), 1);
}

#[test]
fn model_consumes_explicit_detached_geometry() {
    let molecule = one_smiles("CO");
    let positions = Positions::new(Quantity::new(
        molecule
            .atom_ids()
            .enumerate()
            .map(|(index, _)| Point3::new(index as f64, 0.0, 0.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap();
    let model = Model::from_molecule(&molecule, &positions).unwrap();
    assert_eq!(model.atom_count(), molecule.atom_count());
    assert_eq!(model.positions().len(), molecule.atom_count());
}

#[test]
fn topology_changes_require_an_editor_and_invalidate_perception() {
    let mut molecule = one_smiles("CC");
    molecule.perceive().unwrap();
    assert!(molecule.perception().has_valence());
    let mut editor = molecule.edit();
    let oxygen = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    editor
        .add_bond(kekule::core::AtomId::new(1), oxygen, BondOrder::Single)
        .unwrap();
    let edited = editor.finish().unwrap();
    assert_eq!(edited.atom_count(), 3);
    assert!(!edited.perception().has_valence());
}

#[test]
fn editor_enforces_molecule_publication_invariants() {
    assert!(matches!(
        MoleculeEditor::new().finish(),
        Err(kekule::core::MoleculePublicationError::EmptyGraph)
    ));

    let mut single_atom = MoleculeEditor::new();
    single_atom
        .add_atom(Atom::new(Element::from_symbol("He").unwrap()))
        .unwrap();
    assert_eq!(single_atom.finish().unwrap().atom_count(), 1);

    let mut disconnected = MoleculeEditor::new();
    disconnected
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    disconnected
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    assert!(matches!(
        disconnected.finish(),
        Err(kekule::core::MoleculePublicationError::DisconnectedGraph(_))
    ));
}
