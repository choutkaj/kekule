use kekule::core::MoleculeEditor;
use kekule::geometry::Point3;
use kekule::smiles;
use kekule::structure::{Model, Positions};
use kekule::units::{Quantity, ScaleValue, ANGSTROM, NANOMETER};

#[test]
fn borrowed_quantity_conversion_retains_the_source() {
    let source = Quantity::new(vec![1.0, 2.0], NANOMETER);
    let converted = source.to_unit(ANGSTROM).unwrap();
    assert_eq!(source.value(), &[1.0, 2.0]);
    assert_eq!(source.unit(), NANOMETER);
    assert_eq!(converted.value(), &[10.0, 20.0]);
    assert_eq!(converted.unit(), ANGSTROM);
    assert_eq!(source.value_in(ANGSTROM).unwrap(), converted.into_value());
}

#[test]
fn consuming_quantity_conversions_support_values_that_cannot_be_cloned() {
    struct NonClone(f64);
    impl ScaleValue for NonClone {
        fn scaled(self, factor: f64) -> Self {
            Self(self.0 * factor)
        }
    }
    let converted = Quantity::new(NonClone(1.25), NANOMETER)
        .into_unit(ANGSTROM)
        .unwrap()
        .into_value();
    assert_eq!(converted.0, 12.5);

    let text = String::from("owned payload");
    let pointer = text.as_ptr();
    let moved = Quantity::new(text, ANGSTROM).into_value();
    assert_eq!(moved.as_ptr(), pointer);
}

#[test]
fn borrowed_document_conversion_keeps_source_and_consuming_interpretation_moves_storage() {
    let document = smiles::parse_str("CC").unwrap();
    let first = document.to_molecules().unwrap();
    let second = document.to_molecules().unwrap();
    assert_eq!(first, second);
    assert_eq!(document.source(), "CC");

    let interpretation = document.interpret().unwrap();
    let id = interpretation
        .molecule()
        .unwrap()
        .atom_ids()
        .next()
        .unwrap();
    let atom_storage = interpretation.molecule().unwrap().atom(id).unwrap() as *const _;
    let molecule = interpretation.into_molecule().unwrap();
    assert_eq!(molecule.atom(id).unwrap() as *const _, atom_storage);
    let editor = molecule.into_editor();
    assert_eq!(editor.atom(id).unwrap() as *const _, atom_storage);
}

#[test]
fn copy_view_materialization_retains_owner_and_creates_independent_geometry() {
    let molecule = smiles::to_molecules("CC").unwrap().pop().unwrap();
    let source = Model::from_molecule(
        &molecule,
        &Positions::new(Quantity::new(
            vec![Point3::origin(), Point3::new(1.5, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    let view = source.view();
    let mut owned = view.to_model();
    assert_eq!(view.atom_count(), source.atom_count());
    assert!(std::sync::Arc::ptr_eq(
        &owned.shared_topology(),
        &source.shared_topology()
    ));
    assert_ne!(
        owned.positions().values().value().as_ptr(),
        source.positions().values().value().as_ptr()
    );
    let atom = source.atom_ids()[0];
    owned
        .set_position(atom, Quantity::new(Point3::new(9.0, 0.0, 0.0), ANGSTROM))
        .unwrap();
    assert_eq!(view.position(atom).unwrap(), source.position(atom).unwrap());
    assert_ne!(
        owned.position(atom).unwrap(),
        source.position(atom).unwrap()
    );
}

#[test]
fn failed_publication_returns_the_editor_by_ownership_transfer() {
    let failed = MoleculeEditor::new().try_finish().unwrap_err();
    let mut editor = failed.into_editor();
    editor
        .append_molecule(&smiles::to_molecules("C").unwrap().pop().unwrap())
        .unwrap();
    assert_eq!(editor.finish().unwrap().atom_count(), 1);
}
