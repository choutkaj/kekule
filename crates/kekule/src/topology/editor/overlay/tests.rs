use super::*;
use crate::topology::edit::EditAtomId;

fn handle(position: u64) -> EditAtomId {
    EditAtomId::from_position(position)
}

#[test]
fn overlays_derive_source_entries_until_changed() {
    let block = HandleBlock { base: 10, len: 4 };
    let derive = |index: usize| index * 100;
    let mut overlay = Overlay::<EditAtomId, usize>::new(block);
    assert_eq!(overlay.len(), 4);
    assert_eq!(overlay.get(&handle(11), derive), Some(100));
    *overlay.get_mut(&handle(12), derive).unwrap() += 1;
    assert_eq!(overlay.remove(&handle(10), derive), Some(0));
    assert_eq!(overlay.remove(&handle(10), derive), None);
    overlay.insert(handle(20), 7);
    assert_eq!(
        overlay.iter(derive).collect::<Vec<_>>(),
        [
            (handle(11), 100),
            (handle(12), 201),
            (handle(13), 300),
            (handle(20), 7)
        ]
    );
    assert_eq!(overlay.keys().len(), 4);
    let runs = overlay
        .runs()
        .map(|run| match run {
            Run::Source(range) => format!("{range:?}"),
            Run::Stored(key, value) => format!("{}={value}", key.position()),
        })
        .collect::<Vec<_>>();
    assert_eq!(runs, ["1..2", "12=201", "3..4", "20=7"]);
    overlay.clear();
    assert!(overlay.is_empty());
    assert_eq!(overlay.iter(derive).count(), 0);
    assert_eq!(overlay.get(&handle(11), derive), None);
}

#[test]
fn bit_sets_fill_partial_words() {
    let mut set = BitSet::default();
    set.fill(70);
    assert_eq!(set.len(), 70);
    assert!(set.contains(69) && !set.contains(70));
    assert!(set.insert(70));
    assert!(!set.insert(3));
    assert_eq!(set.len(), 71);
}
