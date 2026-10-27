//! Ordered maps whose keys are mostly dense integers.
//!
//! A draft of a large topology holds one entry per source atom, bond and
//! hierarchy node. Their keys are a contiguous block of handles allocated
//! together, or small local and source IDs, so vectors index them directly;
//! keys outside that range fall back to a `BTreeMap`. Iteration is in key
//! order, as for a `BTreeMap`.

use std::collections::{btree_map, BTreeMap};
use std::fmt;
use std::iter::FusedIterator;
use std::ops::Index;

/// A key with a position in a dense integer space.
pub(crate) trait DenseKey: Copy + Ord {
    fn position(self) -> u64;
    fn from_position(position: u64) -> Self;
}

/// An ordered map storing keys in `[base, base + slots.len())` in a vector.
///
/// With `grow`, inserting a key shortly past the vector extends it, which
/// suits local IDs allocated in order; handles allocated after a block go to
/// the sparse map instead.
#[derive(Clone)]
pub(crate) struct DenseMap<K, V> {
    base: u64,
    slots: Vec<Option<(K, V)>>,
    dense_len: usize,
    sparse: BTreeMap<K, V>,
    grow: bool,
}

impl<K, V> Default for DenseMap<K, V> {
    fn default() -> Self {
        Self {
            base: 0,
            slots: Vec::new(),
            dense_len: 0,
            sparse: BTreeMap::new(),
            grow: true,
        }
    }
}

impl<K: DenseKey + fmt::Debug, V: fmt::Debug> fmt::Debug for DenseMap<K, V> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.debug_map().entries(self.iter()).finish()
    }
}

impl<K: DenseKey, V> DenseMap<K, V> {
    /// An empty map whose dense block starts at zero and grows with the keys.
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// An empty map reserving vector slots for the `len` consecutive keys
    /// from `base`; it stores other keys sparsely.
    pub(crate) fn block(base: u64, len: usize) -> Self {
        Self {
            base,
            slots: std::iter::repeat_with(|| None).take(len).collect(),
            dense_len: 0,
            sparse: BTreeMap::new(),
            grow: false,
        }
    }

    fn slot(&self, key: K) -> Option<usize> {
        let offset = key.position().checked_sub(self.base)?;
        usize::try_from(offset)
            .ok()
            .filter(|&offset| offset < self.slots.len())
    }

    pub(crate) fn len(&self) -> usize {
        self.dense_len + self.sparse.len()
    }

    /// An empty map with the same dense block as `other`.
    pub(crate) fn block_of<W>(other: &DenseMap<K, W>) -> Self {
        Self {
            grow: other.grow,
            ..Self::block(other.base, other.slots.len())
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub(crate) fn get(&self, key: &K) -> Option<&V> {
        match self.slot(*key) {
            Some(slot) => self.slots[slot].as_ref().map(|(_, value)| value),
            None => self.sparse.get(key),
        }
    }

    pub(crate) fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        match self.slot(*key) {
            Some(slot) => self.slots[slot].as_mut().map(|(_, value)| value),
            None => self.sparse.get_mut(key),
        }
    }

    pub(crate) fn contains_key(&self, key: &K) -> bool {
        self.get(key).is_some()
    }

    pub(crate) fn insert(&mut self, key: K, value: V) -> Option<V> {
        // Growing over a short gap keeps local IDs dense; the sparse map must
        // hold nothing the grown block would cover.
        if let Some(offset) = key.position().checked_sub(self.base) {
            let len = self.slots.len() as u64;
            if self.grow && self.sparse.is_empty() && offset >= len && offset <= 2 * len + 64 {
                self.slots.resize_with(offset as usize + 1, || None);
            }
        }
        match self.slot(key) {
            Some(slot) => {
                let previous = self.slots[slot].replace((key, value));
                if previous.is_none() {
                    self.dense_len += 1;
                }
                previous.map(|(_, value)| value)
            }
            None => self.sparse.insert(key, value),
        }
    }

    pub(crate) fn remove(&mut self, key: &K) -> Option<V> {
        match self.slot(*key) {
            Some(slot) => {
                let previous = self.slots[slot].take();
                if previous.is_some() {
                    self.dense_len -= 1;
                }
                previous.map(|(_, value)| value)
            }
            None => self.sparse.remove(key),
        }
    }

    pub(crate) fn clear(&mut self) {
        self.slots.clear();
        self.dense_len = 0;
        self.sparse.clear();
    }

    pub(crate) fn iter(&self) -> Iter<'_, K, V> {
        // Keys in the dense block are never sparse, so the sparse keys below
        // and above it complete the order.
        Iter {
            below: self.sparse.range(..K::from_position(self.base)),
            dense: self.slots.iter(),
            above: self.sparse.range(K::from_position(self.end())..),
            remaining: self.len(),
        }
    }

    fn end(&self) -> u64 {
        self.base.saturating_add(self.slots.len() as u64)
    }

    pub(crate) fn keys(&self) -> impl ExactSizeIterator<Item = &K> + DoubleEndedIterator + '_ {
        self.iter().map(|(key, _)| key)
    }

    pub(crate) fn values(&self) -> impl ExactSizeIterator<Item = &V> + DoubleEndedIterator + '_ {
        self.iter().map(|(_, value)| value)
    }
}

impl<K: DenseKey, V: PartialEq> PartialEq for DenseMap<K, V> {
    fn eq(&self, other: &Self) -> bool {
        self.len() == other.len() && self.iter().eq(other.iter())
    }
}

impl<K: DenseKey, V> Index<&K> for DenseMap<K, V> {
    type Output = V;

    fn index(&self, key: &K) -> &V {
        self.get(key).expect("key present in editor map")
    }
}

impl<'a, K: DenseKey, V> IntoIterator for &'a DenseMap<K, V> {
    type Item = (&'a K, &'a V);
    type IntoIter = Iter<'a, K, V>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<K: DenseKey, V> IntoIterator for DenseMap<K, V> {
    type Item = (K, V);
    type IntoIter = std::iter::Chain<
        std::iter::Chain<
            btree_map::IntoIter<K, V>,
            std::iter::Flatten<std::vec::IntoIter<Option<(K, V)>>>,
        >,
        btree_map::IntoIter<K, V>,
    >;

    fn into_iter(mut self) -> Self::IntoIter {
        let above = self.sparse.split_off(&K::from_position(self.end()));
        self.sparse
            .into_iter()
            .chain(self.slots.into_iter().flatten())
            .chain(above)
    }
}

impl<K: DenseKey, V> FromIterator<(K, V)> for DenseMap<K, V> {
    fn from_iter<I: IntoIterator<Item = (K, V)>>(entries: I) -> Self {
        let mut map = Self::new();
        for (key, value) in entries {
            map.insert(key, value);
        }
        map
    }
}

/// Entries of a [`DenseMap`] in key order.
pub(crate) struct Iter<'a, K, V> {
    below: btree_map::Range<'a, K, V>,
    dense: std::slice::Iter<'a, Option<(K, V)>>,
    above: btree_map::Range<'a, K, V>,
    remaining: usize,
}

impl<'a, K, V> Iterator for Iter<'a, K, V> {
    type Item = (&'a K, &'a V);

    fn next(&mut self) -> Option<Self::Item> {
        let next = self
            .below
            .next()
            .or_else(|| {
                self.dense
                    .by_ref()
                    .find_map(|slot| slot.as_ref().map(|(key, value)| (key, value)))
            })
            .or_else(|| self.above.next());
        if next.is_some() {
            self.remaining -= 1;
        }
        next
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<K, V> DoubleEndedIterator for Iter<'_, K, V> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let next = self
            .above
            .next_back()
            .or_else(|| {
                self.dense
                    .by_ref()
                    .rev()
                    .find_map(|slot| slot.as_ref().map(|(key, value)| (key, value)))
            })
            .or_else(|| self.below.next_back());
        if next.is_some() {
            self.remaining -= 1;
        }
        next
    }
}

impl<K, V> ExactSizeIterator for Iter<'_, K, V> {}
impl<K, V> FusedIterator for Iter<'_, K, V> {}

macro_rules! dense_id {
    ($($id:ty),*) => {$(
        impl DenseKey for $id {
            fn position(self) -> u64 {
                u64::from(self.raw())
            }
            fn from_position(position: u64) -> Self {
                // Out-of-range bounds only delimit ranges past every key.
                Self::new(u32::try_from(position).unwrap_or(u32::MAX))
            }
        }
    )*};
}
dense_id!(
    crate::molecule::AtomId,
    crate::molecule::BondId,
    crate::topology::ChainId,
    crate::topology::ResidueId,
    crate::topology::AtomSiteId
);

/// An instance-qualified source ID with a dense local part.
pub(crate) trait Qualified: Copy + Ord {
    type Local: DenseKey;
    fn split(self) -> (usize, Self::Local);
    fn join(instance: usize, local: Self::Local) -> Self;
}

impl Qualified for crate::topology::InstanceAtomId {
    type Local = crate::molecule::AtomId;
    fn split(self) -> (usize, Self::Local) {
        (self.instance().index(), self.atom())
    }
    fn join(instance: usize, local: Self::Local) -> Self {
        Self::new(
            crate::topology::MoleculeInstanceId::new(instance as u32),
            local,
        )
    }
}

impl Qualified for crate::topology::InstanceBondId {
    type Local = crate::molecule::BondId;
    fn split(self) -> (usize, Self::Local) {
        (self.instance().index(), self.bond())
    }
    fn join(instance: usize, local: Self::Local) -> Self {
        Self::new(
            crate::topology::MoleculeInstanceId::new(instance as u32),
            local,
        )
    }
}

/// A map keyed by instance-qualified source IDs: one dense local map per
/// instance, iterated in instance and then local order.
///
/// Like [`DenseMap`], it grows its instance vector only over short gaps, so
/// an untrusted instance ID far past the others is stored sparsely rather
/// than allocating every instance before it.
#[derive(Clone)]
pub(crate) struct QualifiedMap<K: Qualified, V> {
    instances: Vec<DenseMap<K::Local, V>>,
    // Keys of instances at or past `instances.len()`.
    sparse: BTreeMap<K, V>,
}

impl<K: Qualified, V> Default for QualifiedMap<K, V> {
    fn default() -> Self {
        Self {
            instances: Vec::new(),
            sparse: BTreeMap::new(),
        }
    }
}

impl<K: Qualified + fmt::Debug, V: fmt::Debug> fmt::Debug for QualifiedMap<K, V> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.debug_map().entries(self.iter()).finish()
    }
}

impl<K: Qualified, V> QualifiedMap<K, V> {
    pub(crate) fn get(&self, key: &K) -> Option<&V> {
        let (instance, local) = key.split();
        match self.instances.get(instance) {
            Some(map) => map.get(&local),
            None => self.sparse.get(key),
        }
    }

    pub(crate) fn contains_key(&self, key: &K) -> bool {
        self.get(key).is_some()
    }

    pub(crate) fn insert(&mut self, key: K, value: V) -> Option<V> {
        let (instance, local) = key.split();
        let len = self.instances.len();
        // The sparse map must hold nothing the grown vector would cover.
        if instance >= len && self.sparse.is_empty() && instance <= 2 * len + 64 {
            self.instances.resize_with(instance + 1, DenseMap::new);
        }
        match self.instances.get_mut(instance) {
            Some(map) => map.insert(local, value),
            None => self.sparse.insert(key, value),
        }
    }

    /// The entries of one instance in local order.
    pub(crate) fn instance(&self, instance: usize) -> impl Iterator<Item = (K, &V)> + '_ {
        let dense = self
            .instances
            .get(instance)
            .into_iter()
            .flat_map(move |map| {
                map.iter()
                    .map(move |(&local, value)| (K::join(instance, local), value))
            });
        let sparse = self
            .sparse
            .iter()
            .filter(move |(key, _)| key.split().0 == instance)
            .map(|(&key, value)| (key, value));
        dense.chain(sparse)
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = (K, &V)> + '_ {
        // Sparse keys all follow the dense instances.
        (0..self.instances.len())
            .flat_map(|instance| self.instance(instance))
            .chain(self.sparse.iter().map(|(&key, value)| (key, value)))
    }
}

impl<K: Qualified + PartialEq, V: PartialEq> PartialEq for QualifiedMap<K, V> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().eq(other.iter())
    }
}

impl<K: Qualified, V> FromIterator<(K, V)> for QualifiedMap<K, V> {
    fn from_iter<I: IntoIterator<Item = (K, V)>>(entries: I) -> Self {
        let mut map = Self::default();
        for (key, value) in entries {
            map.insert(key, value);
        }
        map
    }
}

impl<K: Qualified, V> Index<&K> for QualifiedMap<K, V> {
    type Output = V;

    fn index(&self, key: &K) -> &V {
        self.get(key).expect("source ID present in editor map")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::AtomId;

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
        let mut other = DenseMap::<Key, u64>::block_of(&map);
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
        let keys = map.keys().map(|id| id.raw()).collect::<Vec<_>>();
        assert_eq!(keys, [0, 1, 2, 4, 5, 10_000]);
    }

    #[test]
    fn qualified_maps_iterate_in_instance_then_local_order() {
        use crate::topology::{InstanceAtomId, MoleculeInstanceId};
        let id = |instance, atom| {
            InstanceAtomId::new(MoleculeInstanceId::new(instance), AtomId::new(atom))
        };
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
        let id = |instance, atom| {
            InstanceAtomId::new(MoleculeInstanceId::new(instance), AtomId::new(atom))
        };
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
}
