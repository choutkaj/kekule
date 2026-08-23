use kekule::core::{Element, Molecule, Perception};
use kekule::geometry::Point3;
use kekule::{molfile, sdf, smiles};

const DISCONNECTED_MOLFILE: &str = "salt-like\nkekule\n\n  2  0  0  0  0  0            999 V2000\n    1.2500    2.5000    3.7500 Na  0  3  0  0  0  0  0  0  0  0  0  0\n   -4.0000    5.5000   -6.2500 Cl  0  5  0  0  0  0  0  0  0  0  0  0\nM  END\n";

fn element(molecule: &Molecule) -> Element {
    molecule
        .atoms()
        .next()
        .expect("one-atom component")
        .1
        .element
}

fn model_molecules(model: &kekule::structure::Model) -> Vec<&Molecule> {
    model
        .topology()
        .molecules()
        .map(|instance| instance.molecule())
        .collect()
}

#[test]
fn parsed_smiles_document_converts_components_in_source_order() {
    let document = smiles::parse_str("[Na+].[Cl-]").expect("SMILES parses");
    let molecules = document.to_molecules().expect("SMILES interprets");

    assert_eq!(molecules.len(), 2);
    assert_eq!(element(&molecules[0]).symbol(), "Na");
    assert_eq!(element(&molecules[1]).symbol(), "Cl");
    assert!(molecules.iter().all(|molecule| molecule.atom_count() == 1));
    assert!(molecules
        .iter()
        .all(|molecule| molecule.perception() == &Perception::default()));
}

#[test]
fn molfile_document_model_retains_published_component_geometry() {
    let document = molfile::parse_str(DISCONNECTED_MOLFILE).expect("Molfile parses");
    let interpretation = molfile::interpret(&document).expect("Molfile report interpretation");
    let expected = [Point3::new(1.25, 2.5, 3.75), Point3::new(-4.0, 5.5, -6.25)];
    for (component, expected) in interpretation.components().iter().zip(expected) {
        let mapping = component
            .report()
            .atom_mappings()
            .first()
            .expect("one source atom mapping");
        assert_eq!(
            component
                .conformer()
                .position(mapping.atom())
                .expect("mapped source coordinate")
                .to_value(),
            expected
        );
    }
    let molecules = document.to_molecules().expect("Molfile interprets");
    let model = document.to_model().expect("Molfile model builds");

    assert_eq!(molecules.len(), 2);
    assert_eq!(model.topology().instance_count(), 2);
    assert_eq!(model.atom_count(), 2);
    assert_eq!(
        model_molecules(&model),
        molecules.iter().collect::<Vec<_>>()
    );
    assert_eq!(
        model.positions().values().value(),
        &[Point3::new(1.25, 2.5, 3.75), Point3::new(-4.0, 5.5, -6.25),]
    );
    assert_eq!(element(model_molecules(&model)[0]).symbol(), "Na");
    assert_eq!(element(model_molecules(&model)[1]).symbol(), "Cl");
    assert!(model_molecules(&model)
        .iter()
        .all(|molecule| molecule.perception() == &Perception::default()));
}

#[test]
fn sdf_parsed_records_remain_independent_conversion_boundaries() {
    let input = format!(
        "{DISCONNECTED_MOLFILE}>  <ROLE>\nions\n\n$$$$\nsecond\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    9.0000    8.0000    7.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n"
    );
    let document = sdf::parse_str(&input, sdf::SdfParseOptions::default()).expect("SDF parses");
    let records: &[sdf::SdfRecord] = document.records();

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].source_record_number(), 1);
    assert_eq!(records[1].source_record_number(), 2);
    assert_eq!(records[0].data_fields()[0].value(), "ions");

    let first_molecules = records[0].to_molecules().expect("first record interprets");
    let first_model = records[0].to_model().expect("first record model builds");
    let second_molecules = records[1].to_molecules().expect("second record interprets");
    let second_model = records[1].to_model().expect("second record model builds");

    assert_eq!(first_molecules.len(), 2);
    assert_eq!(first_model.topology().instance_count(), 2);
    assert_eq!(
        model_molecules(&first_model),
        first_molecules.iter().collect::<Vec<_>>()
    );
    assert_eq!(
        first_model.positions().values().value(),
        &[Point3::new(1.25, 2.5, 3.75), Point3::new(-4.0, 5.5, -6.25),]
    );

    assert_eq!(second_molecules.len(), 1);
    assert_eq!(second_model.topology().instance_count(), 1);
    assert_eq!(
        second_model.positions().values().value(),
        &[Point3::new(9.0, 8.0, 7.0)]
    );
    assert!(first_molecules
        .iter()
        .chain(second_molecules.iter())
        .all(|molecule| molecule.perception() == &Perception::default()));
    assert!(model_molecules(&first_model)
        .into_iter()
        .chain(model_molecules(&second_model))
        .all(|molecule| molecule.perception() == &Perception::default()));
}
