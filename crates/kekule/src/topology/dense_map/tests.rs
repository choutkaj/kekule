use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Key(u64);

impl DenseKey for Key {
    fn position(self) -> u64 {
        self.0
    }
    fn from_position(position: u64) -> Self {
        Self(position)
    }
}

#[test]
fn block_maps_keep_key_order_across_dense_and_sparse_keys() {
    let mut map = DenseMap::block(10, 5);
    for key in [12, 30, 3, 10, 14, 11, 99] {
        assert_eq!(map.insert(Key(key), key * 2), None);
    }
    assert_eq!(map.insert(Key(12), 0), Some(24));
    assert_eq!(map.remove(&Key(11)), Some(22));
    assert_eq!(map.remove(&Key(11)), None);
    assert_eq!(map.remove(&Key(30)), Some(60));
    let expected = [(3, 6), (10, 20), (12, 0), (14, 28), (99, 198)];
    let entries = map.iter().map(|(k, &v)| (k.0, v)).collect::<Vec<_>>();
    assert_eq!(entries, expected);
    assert_eq!(map.len(), expected.len());
    assert_eq!(map.iter().len(), expected.len());
    let reversed = map.iter().rev().map(|(k, _)| k.0).collect::<Vec<_>>();
    assert_eq!(reversed, [99, 14, 12, 10, 3]);
    assert!(map.contains_key(&Key(14)) && !map.contains_key(&Key(13)));
    assert_eq!(map[&Key(99)], 198);
    *map.get_mut(&Key(3)).unwrap() += 1;
    let owned = map
        .clone()
        .into_iter()
        .map(|(k, v)| (k.0, v))
        .collect::<Vec<_>>();
    assert_eq!(owned, [(3, 7), (10, 20), (12, 0), (14, 28), (99, 198)]);
    let mut other = DenseMap::<Key, u64>::block(10, 5);
    for (k, &v) in &map {
        other.insert(*k, v);
    }
    assert_eq!(other, map);
    map.clear();
    assert!(map.is_empty() && map.iter().next().is_none());
}

#[test]
fn growing_maps_stay_dense_over_short_gaps() {
    let mut map = DenseMap::new();
    for raw in [0, 1, 2, 5, 4] {
        map.insert(AtomId::new(raw), raw);
    }
    assert!(map.sparse.is_empty());
    map.insert(AtomId::new(10_000), 7);
    assert_eq!(map.sparse.len(), 1);
    let keys = map.iter().map(|(id, _)| id.raw()).collect::<Vec<_>>();
    assert_eq!(keys, [0, 1, 2, 4, 5, 10_000]);
}

#[test]
fn qualified_maps_iterate_in_instance_then_local_order() {
    use crate::topology::{InstanceAtomId, MoleculeInstanceId};
    let id =
        |instance, atom| InstanceAtomId::new(MoleculeInstanceId::new(instance), AtomId::new(atom));
    let map = [(id(2, 1), 'c'), (id(0, 3), 'b'), (id(0, 0), 'a')]
        .into_iter()
        .collect::<QualifiedMap<_, _>>();
    let entries = map.iter().map(|(id, &v)| (id, v)).collect::<Vec<_>>();
    assert_eq!(entries, [(id(0, 0), 'a'), (id(0, 3), 'b'), (id(2, 1), 'c')]);
    assert_eq!(map.instance(0).count(), 2);
    assert_eq!(map.instance(1).count(), 0);
    assert_eq!(map.get(&id(2, 1)), Some(&'c'));
    assert_eq!(map.get(&id(5, 0)), None);
}

#[test]
fn qualified_maps_store_far_instances_sparsely() {
    use crate::topology::{InstanceAtomId, MoleculeInstanceId};
    let id =
        |instance, atom| InstanceAtomId::new(MoleculeInstanceId::new(instance), AtomId::new(atom));
    let mut map = QualifiedMap::default();
    map.insert(id(0, 0), 'a');
    map.insert(id(u32::MAX, 0), 'z');
    map.insert(id(1, u32::MAX), 'y');
    assert!(map.instances.len() < 100);
    assert_eq!(map.get(&id(u32::MAX, 0)), Some(&'z'));
    assert_eq!(map.get(&id(1, u32::MAX)), Some(&'y'));
    let keys = map.iter().map(|(id, _)| id).collect::<Vec<_>>();
    assert_eq!(keys, [id(0, 0), id(1, u32::MAX), id(u32::MAX, 0)]);
    assert_eq!(map.instance(u32::MAX as usize).count(), 1);
}
