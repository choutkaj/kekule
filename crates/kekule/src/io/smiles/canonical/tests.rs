use super::*;
use crate::diagnostics::DiagnosticKind;

#[test]
fn canonical_input_complexity_accepts_the_boundary() {
    let mut molecule = crate::smiles::to_molecules("CC").unwrap().pop().unwrap();
    molecule.perceive().unwrap();
    assert_eq!(
        write_canonical_smiles_with_limits(&molecule, 16, usize::MAX).unwrap(),
        "CC"
    );
    let error = write_canonical_smiles_with_limits(&molecule, 15, usize::MAX).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    assert_eq!(
        error.limit_exceeded(),
        Some(LimitExceeded::new(Resource::Work, 15, 16))
    );
}

#[test]
fn canonical_export_bounds_compact_graph_storage_before_scratch_allocation() {
    let mut editor = MoleculeEditor::new();
    let carbon = Atom::new(Element::from_symbol("C").unwrap());
    editor.add_atom(carbon.clone()).unwrap();
    for _ in 0..8 {
        let deleted = editor.add_atom(carbon.clone()).unwrap();
        editor.delete_atom(deleted).unwrap();
    }
    let mut molecule = editor.finish().unwrap();
    // Deleted draft slots do not survive publication into scratch sizing.
    assert_eq!(molecule.graph().atom_slot_count(), 1);
    molecule.perceive().unwrap();
    assert_eq!(
        write_canonical_smiles_with_limits(&molecule, usize::MAX, 1).unwrap(),
        "C"
    );
    let mut chain = crate::smiles::to_molecules("CCC").unwrap().pop().unwrap();
    chain.perceive().unwrap();
    let error = write_canonical_smiles_with_limits(&chain, usize::MAX, 4).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    assert_eq!(
        error.limit_exceeded(),
        Some(LimitExceeded::new(Resource::GraphSlots, 4, 5))
    );
}
