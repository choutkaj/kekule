//! Dense topology orders with constant-time qualified-ID lookup.

use super::{InstanceAtomId, InstanceBondId, MoleculeInstanceId};

/// An instance-qualified local identifier stored in a [`DenseOrder`].
pub(super) trait QualifiedId: Copy + Eq {
    fn instance(self) -> MoleculeInstanceId;
    fn local(self) -> usize;
}

impl QualifiedId for InstanceAtomId {
    fn instance(self) -> MoleculeInstanceId {
        self.instance()
    }
    fn local(self) -> usize {
        self.atom().index()
    }
}

impl QualifiedId for InstanceBondId {
    fn instance(self) -> MoleculeInstanceId {
        self.instance()
    }
    fn local(self) -> usize {
        self.bond().index()
    }
}

/// The lookup entry of a local slot without a live ID.
pub(super) const ABSENT: u32 = u32::MAX;

/// One dense order over every live local ID of every instance.
///
/// `order[dense]` is the qualified ID at a dense index. Instance `i` owns
/// `lookup[starts[i]..starts[i + 1]]`, one entry per local slot of its
/// definition, holding the slot's dense index or [`ABSENT`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct DenseOrder<Id> {
    order: Vec<Id>,
    starts: Vec<usize>,
    lookup: Vec<u32>,
}

impl<Id> Default for DenseOrder<Id> {
    fn default() -> Self {
        Self {
            order: Vec::new(),
            starts: Vec::new(),
            lookup: Vec::new(),
        }
    }
}

/// A rejected replacement order for a [`DenseOrder`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum OrderError<Id> {
    Length { expected: usize, actual: usize },
    Unknown(Id),
    Duplicate(Id),
}

impl<Id: QualifiedId> DenseOrder<Id> {
    pub(super) fn len(&self) -> usize {
        self.order.len()
    }

    pub(super) fn order(&self) -> &[Id] {
        &self.order
    }

    pub(super) fn id(&self, dense: usize) -> Option<Id> {
        self.order.get(dense).copied()
    }

    pub(super) fn index(&self, id: Id) -> Option<usize> {
        let instance = id.instance().index();
        let start = *self.starts.get(instance)?;
        let end = self
            .starts
            .get(instance + 1)
            .copied()
            .unwrap_or(self.lookup.len());
        let slot = start.checked_add(id.local()).filter(|slot| *slot < end)?;
        let dense = self.lookup[slot];
        (dense != ABSENT).then_some(dense as usize)
    }

    /// Whether `additional` more IDs keep every dense index below the sentinel.
    pub(super) fn can_extend(&self, additional: usize) -> bool {
        self.order
            .len()
            .checked_add(additional)
            .is_some_and(|len| len <= ABSENT as usize)
    }

    /// Appends the next instance, ordering its live IDs after existing ones.
    ///
    /// `slots` is its definition's local slot count. Callers check
    /// [`Self::can_extend`] first.
    pub(super) fn push_instance(&mut self, slots: usize, ids: impl IntoIterator<Item = Id>) {
        let start = self.lookup.len();
        self.starts.push(start);
        self.lookup.resize(start + slots, ABSENT);
        for id in ids {
            debug_assert_eq!(id.instance().index() + 1, self.starts.len());
            debug_assert!(self.order.len() < ABSENT as usize);
            self.lookup[start + id.local()] = self.order.len() as u32;
            self.order.push(id);
        }
    }

    /// The order `order` over instances whose definitions have the given
    /// local slot counts, in instance order. Each ID names a distinct slot
    /// of its instance.
    pub(super) fn from_order(order: Vec<Id>, slots: impl IntoIterator<Item = usize>) -> Self {
        let mut starts = Vec::new();
        let mut total = 0;
        for count in slots {
            starts.push(total);
            total += count;
        }
        let mut lookup = vec![ABSENT; total];
        for (dense, &id) in order.iter().enumerate() {
            let slot = starts[id.instance().index()] + id.local();
            debug_assert_eq!(lookup[slot], ABSENT, "each ID is listed once");
            lookup[slot] = dense as u32;
        }
        Self {
            order,
            starts,
            lookup,
        }
    }

    /// The dense index of each local slot of a valid instance, or [`ABSENT`].
    pub(super) fn instance_lookup(&self, instance: usize) -> &[u32] {
        let start = self.starts[instance];
        let end = self
            .starts
            .get(instance + 1)
            .copied()
            .unwrap_or(self.lookup.len());
        &self.lookup[start..end]
    }

    /// Replaces the order with a permutation of the current IDs.
    ///
    /// Returns, for each new dense index, its former dense index. A rejected
    /// order leaves the dense order unchanged.
    pub(super) fn reorder(
        &mut self,
        order: impl IntoIterator<Item = Id>,
    ) -> Result<Vec<usize>, OrderError<Id>> {
        let order = order.into_iter().collect::<Vec<_>>();
        if order.len() != self.order.len() {
            return Err(OrderError::Length {
                expected: self.order.len(),
                actual: order.len(),
            });
        }
        let mut seen = vec![false; self.order.len()];
        let mut previous = Vec::with_capacity(order.len());
        for &id in &order {
            let index = self.index(id).ok_or(OrderError::Unknown(id))?;
            if std::mem::replace(&mut seen[index], true) {
                return Err(OrderError::Duplicate(id));
            }
            previous.push(index);
        }
        for (dense, &id) in order.iter().enumerate() {
            let slot = self.starts[id.instance().index()] + id.local();
            self.lookup[slot] = dense as u32;
        }
        self.order = order;
        Ok(previous)
    }
}
