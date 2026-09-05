use kekule::core::{Element, Molecule, Perception};
use kekule::geometry::Point3;
use kekule::{molfile, sdf, smiles};

const DISCONNECTED_MOLFILE: &str = "salt-like\nkekule\n\n  2  0  0  0  0  0            999 V2000\n    1.2500    2.5000    3.7500 Na  0  3  0  0  0  0  0  0  0  0  0  0\n   -4.0000    5.5000   -6.2500 Cl  0  5  0  0  0  0  0  0  0  0  0  0\nM  END\n";

const INTERLEAVED_COMPONENT_MOLFILE: &str = "interleaved\nkekule\n\n  4  1  0  0  0  0            999 V2000\n   10.0000   11.0000   12.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   20.0000   21.0000   22.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n   30.0000   31.0000   32.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   40.0000   41.0000   42.0000 Cl  0  5  0  0  0  0  0  0  0  0  0  0\n  1  3  1  0  0  0  0\nM  END\n";

#[test]
fn sdf_model_can_install_perception_without_losing_source_order_or_geometry() {
    let text = format!("{INTERLEAVED_COMPONENT_MOLFILE}$$$$\n");
    let document = sdf::parse_str(&text).unwrap();
    let mut model = document.records()[0].to_model().unwrap();
    let original = model.clone();
    let positions_ptr = model.positions().values().value().as_ptr();
    assert!(model_molecules(&model)
        .iter()
        .all(|m| m.perception() == &Perception::default()));

    model.perceive().unwrap();

    assert!(model.topology().same_layout(original.topology()));
    assert_eq!(model.topology().atom_ids(), original.topology().atom_ids());
    assert_eq!(
        model.topology().hierarchy(),
        original.topology().hierarchy()
    );
    assert_eq!(model.positions(), original.positions());
    assert_eq!(model.positions().values().value().as_ptr(), positions_ptr);
    for (perceived, source) in model_molecules(&model)
        .into_iter()
        .zip(model_molecules(&original))
    {
        let mut expected = source.clone();
        expected.perceive().unwrap();
        assert_eq!(perceived, source);
        assert_eq!(perceived.perception(), expected.perception());
        assert_eq!(source.perception(), &Perception::default());
    }
}

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

fn assert_point_close(actual: Point3, expected: Point3) {
    assert!((actual.x - expected.x).abs() < 1.0e-15);
    assert!((actual.y - expected.y).abs() < 1.0e-15);
    assert!((actual.z - expected.z).abs() < 1.0e-15);
}

fn assert_points_close(actual: &[Point3], expected: &[Point3]) {
    assert_eq!(actual.len(), expected.len());
    for (&actual, &expected) in actual.iter().zip(expected) {
        assert_point_close(actual, expected);
    }
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
    assert_eq!(
        document,
        molfile::parse_str_with_options(
            DISCONNECTED_MOLFILE,
            molfile::MolfileParseOptions::default(),
        )
        .expect("explicit default Molfile parse")
    );
    let interpretation = document.interpret().expect("Molfile report interpretation");
    let expected = [
        Point3::new(0.125, 0.25, 0.375),
        Point3::new(-0.4, 0.55, -0.625),
    ];
    for (component, expected) in interpretation.components().iter().zip(expected) {
        let mapping = component
            .report()
            .atom_mappings()
            .first()
            .expect("one source atom mapping");
        let position_index = component
            .molecule()
            .atom_ids()
            .position(|atom| atom == mapping.atom())
            .expect("mapped canonical atom");
        assert_point_close(
            component
                .positions()
                .position_at(position_index)
                .expect("mapped source coordinate")
                .to_value(),
            expected,
        );
    }
    let molecules = document.to_molecules().expect("Molfile interprets");
    let topology = document.to_topology().expect("Molfile topology builds");
    let model = document.to_model().expect("Molfile model builds");

    assert_eq!(molecules.len(), 2);
    assert_eq!(model.topology().instance_count(), 2);
    assert_eq!(model.atom_count(), 2);
    assert!(topology.same_layout(model.topology()));
    assert_eq!(
        model_molecules(&model),
        molecules.iter().collect::<Vec<_>>()
    );
    assert_points_close(
        model.positions().values().value(),
        &[
            Point3::new(0.125, 0.25, 0.375),
            Point3::new(-0.4, 0.55, -0.625),
        ],
    );
    assert_eq!(element(model_molecules(&model)[0]).symbol(), "Na");
    assert_eq!(element(model_molecules(&model)[1]).symbol(), "Cl");
    assert!(model_molecules(&model)
        .iter()
        .all(|molecule| molecule.perception() == &Perception::default()));
    let hierarchy = model.topology().hierarchy();
    let (_, chain) = hierarchy.chains().next().expect("synthetic chain");
    assert_eq!(chain.label_id(), "A");
    assert_eq!(chain.author_id(), None);
    assert_eq!(chain.residues().len(), 2);
    for (index, residue_id) in chain.residues().iter().enumerate() {
        let residue = hierarchy.residue(*residue_id).unwrap();
        assert_eq!(residue.name(), "UNL");
        assert_eq!(residue.label_seq_id(), Some((index + 1) as i32));
        assert_eq!(
            residue.author_seq_id(),
            Some((index + 1).to_string().as_str())
        );
        assert_eq!(residue.atom_sites().len(), 1);
    }
}

#[test]
fn molfile_model_remaps_interleaved_source_atoms_to_component_positions() {
    let document = molfile::parse_str(INTERLEAVED_COMPONENT_MOLFILE).expect("Molfile parses");
    let molecules = document.to_molecules().expect("Molfile interprets");
    let model = document.to_model().expect("Molfile model builds");

    assert_eq!(molecules.len(), 3);
    assert_eq!(model.topology().instance_count(), 3);
    assert_eq!(
        model_molecules(&model),
        molecules.iter().collect::<Vec<_>>()
    );
    assert_eq!(
        molecules[0]
            .atoms()
            .map(|(_, atom)| atom.element.symbol())
            .collect::<Vec<_>>(),
        ["C", "O"]
    );
    assert_eq!(element(&molecules[1]).symbol(), "Na");
    assert_eq!(element(&molecules[2]).symbol(), "Cl");
    assert_points_close(
        model.positions().values().value(),
        &[
            Point3::new(1.0, 1.1, 1.2),
            Point3::new(3.0, 3.1, 3.2),
            Point3::new(2.0, 2.1, 2.2),
            Point3::new(4.0, 4.1, 4.2),
        ],
    );
}

#[test]
fn sdf_parsed_records_remain_independent_conversion_boundaries() {
    let input = format!(
        "{DISCONNECTED_MOLFILE}>  <ROLE>\nions\n\n$$$$\nsecond\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    9.0000    8.0000    7.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n"
    );
    let document = sdf::parse_str(&input).expect("SDF parses");
    assert_eq!(
        document,
        sdf::parse_str_with_options(&input, sdf::SdfParseOptions::default())
            .expect("explicit default SDF parse")
    );
    let records: &[sdf::SdfRecord] = document.records();

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].source_record_number(), 1);
    assert_eq!(records[1].source_record_number(), 2);
    assert_eq!(records[0].data_fields()[0].value(), "ions");

    let first_interpretation = records[0].interpret().expect("rich record interpretation");
    assert_eq!(first_interpretation.title(), "salt-like");
    assert_eq!(first_interpretation.data_fields()[0].value(), "ions");
    assert_eq!(first_interpretation.report().record(), 1);
    assert_eq!(first_interpretation.report().molfile_components().len(), 2);
    assert_eq!(first_interpretation.molecules().count(), 2);
    assert!(first_interpretation
        .topology()
        .same_layout(first_interpretation.model().topology()));
    assert_points_close(
        first_interpretation.model().positions().values().value(),
        &[
            Point3::new(0.125, 0.25, 0.375),
            Point3::new(-0.4, 0.55, -0.625),
        ],
    );
    let projected_topology = first_interpretation.clone().to_topology();
    let projected_model = first_interpretation.clone().to_model();
    let projected_molecules = first_interpretation.to_molecules();
    assert!(projected_topology.same_layout(projected_model.topology()));
    assert_eq!(
        projected_topology
            .molecules()
            .map(|occurrence| occurrence.molecule())
            .collect::<Vec<_>>(),
        projected_molecules.iter().collect::<Vec<_>>()
    );
    assert_eq!(projected_topology.hierarchy().chains().count(), 1);
    assert_eq!(projected_topology.hierarchy().residues().count(), 2);

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
    assert_points_close(
        first_model.positions().values().value(),
        &[
            Point3::new(0.125, 0.25, 0.375),
            Point3::new(-0.4, 0.55, -0.625),
        ],
    );

    assert_eq!(second_molecules.len(), 1);
    assert_eq!(second_model.topology().instance_count(), 1);
    assert_points_close(
        second_model.positions().values().value(),
        &[Point3::new(0.9, 0.8, 0.7)],
    );
    assert!(first_molecules
        .iter()
        .chain(second_molecules.iter())
        .all(|molecule| molecule.perception() == &Perception::default()));
    assert!(model_molecules(&first_model)
        .into_iter()
        .chain(model_molecules(&second_model))
        .all(|molecule| molecule.perception() == &Perception::default()));

    let document_interpretation = document
        .interpret()
        .expect("document interprets per record");
    assert_eq!(document_interpretation.records().len(), 2);
    assert_eq!(document_interpretation.report().records().len(), 2);
    assert_eq!(document_interpretation.records()[0].report().record(), 1);
    assert_eq!(document_interpretation.records()[1].report().record(), 2);

    let written = sdf::write_v2000(&[document_interpretation.records()[1].clone()])
        .expect("single-component rich record writes");
    let round_trip_document = sdf::parse_str(&written).expect("written record parses");
    let round_trip = round_trip_document.records()[0]
        .interpret()
        .expect("written record interprets");
    assert_eq!(round_trip.title(), "second");
    assert_points_close(
        round_trip.model().positions().values().value(),
        &[Point3::new(0.9, 0.8, 0.7)],
    );
}
