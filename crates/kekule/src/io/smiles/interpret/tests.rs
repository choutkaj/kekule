use super::super::parse::parse_smiles_document;
use super::*;

#[test]
fn dot_smiles_interprets_as_connected_molecules() {
    let document =
        parse_smiles_document("CC(=O)[O-].[Na+]", Default::default()).expect("valid salt");
    let interpretation = interpret_smiles_document(&document, SmilesInterpretOptions::default())
        .expect("interpret salt");
    assert_eq!(interpretation.components().len(), 2);
    assert!(interpretation
        .molecules()
        .all(|molecule| molecule.connected_components().len() == 1));
    assert_eq!(interpretation.components()[0].source_span(), 0..10);
    assert_eq!(interpretation.components()[1].source_span(), 11..16);
}

#[test]
fn component_mappings_retain_document_offsets() {
    let document = parse_smiles_document("C.[Na+]", Default::default()).expect("valid components");
    let interpretation = interpret_smiles_document(&document, SmilesInterpretOptions::default())
        .expect("interpret components");
    assert_eq!(
        interpretation.components()[1].report().atom_mappings()[0].source_span(),
        2..7
    );
}

#[test]
fn single_component_convenience_rejects_dot_smiles_without_panicking() {
    let document = parse_smiles_document("C.O", Default::default()).expect("valid components");
    let error = interpret_smiles_document(&document, SmilesInterpretOptions::default())
        .expect("interpret components")
        .into_molecule()
        .expect_err("single-component access must reject two components");
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::ComponentCount { actual: 2 })
    );
}
