use crate::diagnostics::CapacityExceeded;
use crate::molecule::editor::*;

fn first_unsupported_slot() -> usize {
    usize::try_from(u64::from(u32::MAX) + 1).expect("64-bit usize")
}

#[test]
fn capacity_rejection_does_not_mutate_molecule_state() {
    let mut molecule = crate::molecule::MoleculeEditor::new();
    let before = molecule.clone();
    assert_eq!(
        molecule.add_atom_at_slot(
            Atom::new(Element::from_atomic_number(6).expect("carbon")),
            first_unsupported_slot(),
        ),
        Err(MoleculeError::CapacityExceeded(CapacityExceeded::new(
            Countable::Atom
        )))
    );
    assert_eq!(molecule, before);
}
