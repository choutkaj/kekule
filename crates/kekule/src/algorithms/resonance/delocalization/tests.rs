use super::*;

#[test]
fn long_paths_mark_each_bond_once_and_poll_cancellation() {
    let text = format!("[CH2-]{}", "C=C".repeat(5_000));
    let mut mol = crate::smiles::to_molecules(&text).unwrap().remove(0);
    mol.perceive().unwrap();
    let groups = delocalized_groups(&mol, None).unwrap();
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].bonds.len(), mol.bond_count());
    let cancelled = Cancellation::new(|| true);
    assert!(matches!(
        delocalized_groups(&mol, Some(&cancelled)),
        Err(ResonancePerceptionError::Cancelled(_))
    ));
}
