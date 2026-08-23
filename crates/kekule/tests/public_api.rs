use kekule::core::{Atom, BondOrder, Element, Molecule, MoleculeEditor};
use kekule::geometry::Point3;
use kekule::sdf::{self, SdfParseOptions};
use kekule::structure::{Model, Positions};
use kekule::units::{Quantity, ANGSTROM};

fn one_smiles(input: &str) -> Molecule {
    let mut molecules = Molecule::from_smiles(input).expect("SMILES interprets");
    assert_eq!(molecules.len(), 1);
    molecules.pop().expect("component count was checked")
}

#[test]
fn smiles_partitions_components_in_source_order() {
    let molecules = Molecule::from_smiles("CC.O.[Na+]").expect("dot SMILES interprets");
    assert_eq!(molecules.len(), 3);
    assert_eq!(molecules[0].atom_count(), 2);
    assert_eq!(molecules[1].atoms().next().unwrap().1.element.symbol(), "O");
    assert_eq!(
        molecules[2].atoms().next().unwrap().1.element.symbol(),
        "Na"
    );
}

#[test]
fn sdf_preserves_record_boundaries_and_component_order() {
    let input = "first\nkekule\n\n  2  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\nsecond\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n";
    let document = sdf::parse_str(input, SdfParseOptions::default()).unwrap();
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
fn editor_rejects_empty_publication() {
    assert!(MoleculeEditor::new().finish().is_err());
}
