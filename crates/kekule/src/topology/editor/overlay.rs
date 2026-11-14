//! Draft entity tables that read a source snapshot until an edit changes them.
//!
//! A draft started from a topology names every source atom, bond and
//! hierarchy node by one handle of a consecutive block, in source order. An
//! [`Overlay`] stores nothing for an entity the draft has not changed: its
//! value is derived from the source when read. Changed and removed source
//! entities, and entities the draft added, are stored, so starting a draft
//! costs nothing per entity and an edit costs what it touches.

use std::collections::{btree_map, BTreeMap};
use std::iter::{FusedIterator, Peekable};
use std::ops::Range;

use super::super::dense_map::{DenseKey, DenseMap};

/// A set of indices below some bound, one bit each.
#[derive(Debug, Clone, Default)]
pub(super) struct BitSet {
    words: Vec<u64>,
    len: usize,
}

impl BitSet {
    pub(super) fn contains(&self, index: usize) -> bool {
        self.words
            .get(index / 64)
            .is_some_and(|word| word >> (index % 64) & 1 == 1)
    }

    /// Adds `index`, returning whether it was absent.
    pub(super) fn insert(&mut self, index: usize) -> bool {
        let word = index / 64;
        if word >= self.words.len() {
            self.words.resize(word + 1, 0);
        }
        let bit = 1 << (index % 64);
        let absent = self.words[word] & bit == 0;
        self.words[word] |= bit;
        self.len += usize::from(absent);
        absent
    }

    /// The number of indices in the set.
    pub(super) fn len(&self) -> usize {
        self.len
    }

    /// The first index at or after `index` in the set, if any.
    fn next_at_or_after(&self, index: usize) -> Option<usize> {
        let mut word = index / 64;
        let mut bits = *self.words.get(word)? >> (index % 64) << (index % 64);
        loop {
            if bits != 0 {
                return Some(word * 64 + bits.trailing_zeros() as usize);
            }
            word += 1;
            bits = *self.words.get(word)?;
        }
    }

    /// Adds every index below `bound`.
    pub(super) fn fill(&mut self, bound: usize) {
        self.words = vec![u64::MAX; bound / 64];
        if !bound.is_multiple_of(64) {
            self.words.push((1 << (bound % 64)) - 1);
        }
        self.len = bound;
    }
}

/// One block of consecutive handles, one per source entity in source order.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct HandleBlock {
    pub(super) base: u64,
    pub(super) len: usize,
}

impl HandleBlock {
    /// A fresh block of `len` handles.
    pub(super) fn reserve(len: usize) -> Self {
        Self {
            base: super::identity::fresh_block(len),
            len,
        }
    }

    /// The source index of a handle in this block.
    pub(super) fn source(self, handle: impl DenseKey) -> Option<u32> {
        let offset = handle.position().checked_sub(self.base)?;
        (offset < self.len as u64).then_some(offset as u32)
    }

    /// The handle of the source entity at `index`.
    pub(super) fn handle<K: DenseKey>(self, index: usize) -> K {
        debug_assert!(index < self.len);
        K::from_position(self.base + index as u64)
    }
}

/// Draft entries keyed by handle: one per source entity of a handle block,
/// derived from the source until changed or removed, then the entries the
/// draft added. Iteration is in handle order: source entities in source
/// order, then added entities in creation order.
///
/// Reads of a source entity take `derive`, which computes its unchanged
/// value from its source index.
#[derive(Debug, Clone)]
pub(super) struct Overlay<K: DenseKey, V> {
    block: HandleBlock,
    removed: BitSet,
    changed: BTreeMap<K, V>,
    added: DenseMap<K, V>,
}

impl<K: DenseKey, V> Default for Overlay<K, V> {
    fn default() -> Self {
        Self::new(HandleBlock::default())
    }
}

impl<K: DenseKey, V> Overlay<K, V> {
    /// Entries for every source entity of `block`, none changed.
    pub(super) fn new(block: HandleBlock) -> Self {
        Self {
            block,
            removed: BitSet::default(),
            changed: BTreeMap::new(),
            added: DenseMap::new(),
        }
    }
}

impl<K: DenseKey, V: Clone> Overlay<K, V> {
    pub(super) fn block(&self) -> HandleBlock {
        self.block
    }

    /// The number of live entries.
    pub(super) fn len(&self) -> usize {
        self.block.len - self.removed.len() + self.added.len()
    }

    pub(super) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Whether `key` names a live entry.
    pub(super) fn contains_key(&self, key: &K) -> bool {
        match self.block.source(*key) {
            Some(index) => !self.removed.contains(index as usize),
            None => self.added.contains_key(key),
        }
    }

    pub(super) fn get(&self, key: &K, derive: impl FnOnce(usize) -> V) -> Option<V> {
        match self.block.source(*key) {
            Some(index) => (!self.removed.contains(index as usize)).then(|| {
                self.changed
                    .get(key)
                    .cloned()
                    .unwrap_or_else(|| derive(index as usize))
            }),
            None => self.added.get(key).cloned(),
        }
    }

    /// The stored value of a live entry, copying a source entity's derived
    /// value into storage first.
    pub(super) fn get_mut(&mut self, key: &K, derive: impl FnOnce(usize) -> V) -> Option<&mut V> {
        match self.block.source(*key) {
            Some(index) => (!self.removed.contains(index as usize)).then(|| {
                self.changed
                    .entry(*key)
                    .or_insert_with(|| derive(index as usize))
            }),
            None => self.added.get_mut(key),
        }
    }

    /// Adds an entry under a handle outside the source block.
    pub(super) fn insert(&mut self, key: K, value: V) {
        debug_assert!(self.block.source(key).is_none());
        self.added.insert(key, value);
    }

    /// Removes a live entry, returning its value.
    pub(super) fn remove(&mut self, key: &K, derive: impl FnOnce(usize) -> V) -> Option<V> {
        match self.block.source(*key) {
            Some(index) => self.removed.insert(index as usize).then(|| {
                self.changed
                    .remove(key)
                    .unwrap_or_else(|| derive(index as usize))
            }),
            None => self.added.remove(key),
        }
    }

    /// Removes every entry. Source handles are not reused.
    pub(super) fn clear(&mut self) {
        self.removed.fill(self.block.len);
        self.changed.clear();
        self.added.clear();
    }

    /// Live entries in handle order.
    pub(super) fn iter<F: Fn(usize) -> V>(&self, derive: F) -> Iter<'_, K, V, F> {
        Iter {
            overlay: self,
            next: 0,
            changed: self.changed.iter().peekable(),
            added: self.added.iter(),
            remaining: self.len(),
            derive,
        }
    }

    /// Live entries in handle order, naming unchanged source entities by
    /// source index instead of deriving their values.
    pub(super) fn entries(&self) -> Entries<'_, K, V> {
        Entries {
            overlay: self,
            next: 0,
            changed: self.changed.iter().peekable(),
            added: self.added.iter(),
            remaining: self.len(),
        }
    }

    /// Live entries in handle order, grouping consecutive unchanged source
    /// entities into runs of source indices.
    pub(super) fn runs(&self) -> Runs<'_, K, V> {
        Runs {
            overlay: self,
            next: 0,
            changed: self.changed.iter().peekable(),
            added: self.added.iter(),
        }
    }

    /// Live handles in handle order.
    pub(super) fn keys(&self) -> Keys<'_, K, V> {
        Keys {
            overlay: self,
            next: 0,
            added: self.added.iter(),
            remaining: self.len(),
        }
    }

    /// Entries added by the draft, in handle order.
    pub(super) fn added(&self) -> &DenseMap<K, V> {
        &self.added
    }
}

impl<K: DenseKey, V> Overlay<K, V> {
    fn next_live(&self, mut index: usize) -> Option<usize> {
        // Skip whole removed words, so a cleared draft iterates quickly.
        while index < self.block.len {
            let word = self.removed.words.get(index / 64).copied().unwrap_or(0);
            if word == u64::MAX && index.is_multiple_of(64) {
                index += 64;
                continue;
            }
            if word >> (index % 64) & 1 == 0 {
                return Some(index);
            }
            index += 1;
        }
        None
    }
}

/// A stored draft row, which names its property storage slot.
pub(super) trait Slotted {
    fn slot(&self) -> usize;
}

// An unchanged source entity's slot is its source index.
impl<K: DenseKey, V: Clone + Slotted> crate::properties::HandleSlots for Overlay<K, V> {
    fn slot(&self, position: u64) -> Option<usize> {
        let key = K::from_position(position);
        match self.block.source(key) {
            Some(index) => {
                let index = index as usize;
                (!self.removed.contains(index))
                    .then(|| self.changed.get(&key).map_or(index, Slotted::slot))
            }
            None => self.added.get(&key).map(Slotted::slot),
        }
    }

    fn live_slots(&self) -> Vec<usize> {
        self.entries()
            .map(|(_, entry)| match entry {
                Entry::Source(index) => index,
                Entry::Stored(value) => value.slot(),
            })
            .collect()
    }
}

/// Live entries of an [`Overlay`] in handle order.
pub(super) struct Iter<'a, K: DenseKey, V, F> {
    overlay: &'a Overlay<K, V>,
    next: usize,
    changed: Peekable<btree_map::Iter<'a, K, V>>,
    added: super::super::dense_map::Iter<'a, K, V>,
    remaining: usize,
    derive: F,
}

impl<K: DenseKey, V: Clone, F: Fn(usize) -> V> Iterator for Iter<'_, K, V, F> {
    type Item = (K, V);

    fn next(&mut self) -> Option<(K, V)> {
        let item = match self.overlay.next_live(self.next) {
            Some(index) => {
                self.next = index + 1;
                let key = self.overlay.block.handle::<K>(index);
                // Changed entries are live source entries in key order.
                while self
                    .changed
                    .peek()
                    .is_some_and(|(&changed, _)| changed < key)
                {
                    self.changed.next();
                }
                let value = match self.changed.peek() {
                    Some(&(&changed, value)) if changed == key => {
                        self.changed.next();
                        value.clone()
                    }
                    _ => (self.derive)(index),
                };
                (key, value)
            }
            None => {
                self.next = self.overlay.block.len;
                self.added
                    .next()
                    .map(|(&key, value)| (key, value.clone()))?
            }
        };
        self.remaining -= 1;
        Some(item)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<K: DenseKey, V: Clone, F: Fn(usize) -> V> ExactSizeIterator for Iter<'_, K, V, F> {}
impl<K: DenseKey, V: Clone, F: Fn(usize) -> V> FusedIterator for Iter<'_, K, V, F> {}

/// One live entry of an [`Overlay`].
#[derive(Debug, Clone, Copy)]
pub(super) enum Entry<'a, V> {
    /// An unchanged source entity, by source index.
    Source(usize),
    /// A changed source entity or an added one.
    Stored(&'a V),
}

/// One stretch of live entries of an [`Overlay`].
#[derive(Debug, Clone)]
pub(super) enum Run<'a, K, V> {
    /// Consecutive unchanged live source entities, by source index.
    Source(Range<usize>),
    /// A changed source entity or an added one.
    Stored(K, &'a V),
}

/// Live entries of an [`Overlay`] in handle order, as runs.
pub(super) struct Runs<'a, K: DenseKey, V> {
    overlay: &'a Overlay<K, V>,
    next: usize,
    changed: Peekable<btree_map::Iter<'a, K, V>>,
    added: super::super::dense_map::Iter<'a, K, V>,
}

impl<'a, K: DenseKey, V> Iterator for Runs<'a, K, V> {
    type Item = Run<'a, K, V>;

    fn next(&mut self) -> Option<Run<'a, K, V>> {
        let overlay = self.overlay;
        let Some(start) = overlay.next_live(self.next) else {
            self.next = overlay.block.len;
            let (&key, value) = self.added.next()?;
            return Some(Run::Stored(key, value));
        };
        // Changed entries are live source entries in key order.
        let key = overlay.block.handle::<K>(start);
        while self
            .changed
            .peek()
            .is_some_and(|(&changed, _)| changed < key)
        {
            self.changed.next();
        }
        if let Some(&(&changed, value)) = self.changed.peek() {
            if changed == key {
                self.changed.next();
                self.next = start + 1;
                return Some(Run::Stored(key, value));
            }
        }
        let changed = self
            .changed
            .peek()
            .and_then(|(&changed, _)| overlay.block.source(changed))
            .map_or(overlay.block.len, |index| index as usize);
        let removed = overlay
            .removed
            .next_at_or_after(start)
            .unwrap_or(overlay.block.len);
        let end = changed.min(removed).min(overlay.block.len);
        self.next = end;
        Some(Run::Source(start..end))
    }
}

impl<K: DenseKey, V> FusedIterator for Runs<'_, K, V> {}

/// Live entries of an [`Overlay`] in handle order.
pub(super) struct Entries<'a, K: DenseKey, V> {
    overlay: &'a Overlay<K, V>,
    next: usize,
    changed: Peekable<btree_map::Iter<'a, K, V>>,
    added: super::super::dense_map::Iter<'a, K, V>,
    remaining: usize,
}

impl<'a, K: DenseKey, V> Iterator for Entries<'a, K, V> {
    type Item = (K, Entry<'a, V>);

    fn next(&mut self) -> Option<(K, Entry<'a, V>)> {
        let item = match self.overlay.next_live(self.next) {
            Some(index) => {
                self.next = index + 1;
                let key = self.overlay.block.handle::<K>(index);
                // Changed entries are live source entries in key order.
                while self
                    .changed
                    .peek()
                    .is_some_and(|(&changed, _)| changed < key)
                {
                    self.changed.next();
                }
                match self.changed.peek() {
                    Some(&(&changed, value)) if changed == key => {
                        self.changed.next();
                        (key, Entry::Stored(value))
                    }
                    _ => (key, Entry::Source(index)),
                }
            }
            None => {
                self.next = self.overlay.block.len;
                let (&key, value) = self.added.next()?;
                (key, Entry::Stored(value))
            }
        };
        self.remaining -= 1;
        Some(item)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<K: DenseKey, V> ExactSizeIterator for Entries<'_, K, V> {}
impl<K: DenseKey, V> FusedIterator for Entries<'_, K, V> {}

/// Live handles of an [`Overlay`] in handle order.
pub(super) struct Keys<'a, K: DenseKey, V> {
    overlay: &'a Overlay<K, V>,
    next: usize,
    added: super::super::dense_map::Iter<'a, K, V>,
    remaining: usize,
}

impl<K: DenseKey, V> Iterator for Keys<'_, K, V> {
    type Item = K;

    fn next(&mut self) -> Option<K> {
        let key = match self.overlay.next_live(self.next) {
            Some(index) => {
                self.next = index + 1;
                self.overlay.block.handle::<K>(index)
            }
            None => {
                self.next = self.overlay.block.len;
                *self.added.next()?.0
            }
        };
        self.remaining -= 1;
        Some(key)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<K: DenseKey, V> ExactSizeIterator for Keys<'_, K, V> {}
impl<K: DenseKey, V> FusedIterator for Keys<'_, K, V> {}

#[cfg(test)]
mod tests;
