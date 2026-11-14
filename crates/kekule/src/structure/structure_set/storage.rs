//! Where a [`super::StructureSet`] keeps its geometries: [`InMemory`], or a file.

use std::fmt;

use super::Geometry;

pub(crate) mod sealed {
    /// Restricts [`super::GeometryStorage`] to Kekule's storages.
    pub trait Storage {
        /// The number of geometries stored.
        fn geometry_count(&self) -> usize;
    }
}

/// Where a [`StructureSet`](super::StructureSet) keeps its geometries, its type
/// parameter: [`InMemory`], the default, or an indexed trajectory file
/// ([`IndexedFile`](crate::traj::io::IndexedFile)).
///
/// Every storage shares the topology, the set properties, the order, and
/// the operations on them. Only an in-memory set lends
/// its geometries by reference; any set is read through
/// [`StructureSetRead`](super::StructureSetRead). The trait is sealed.
pub trait GeometryStorage: sealed::Storage + fmt::Debug + Send + Sync {}

/// Geometries held in memory, the default storage of a
/// [`StructureSet`](super::StructureSet): they are lent by reference, and reading
/// them cannot fail.
#[derive(Debug, Clone, Default)]
pub struct InMemory {
    pub(super) geometries: Vec<Geometry>,
}

impl InMemory {
    pub(super) const fn new(geometries: Vec<Geometry>) -> Self {
        Self { geometries }
    }
}

impl sealed::Storage for InMemory {
    fn geometry_count(&self) -> usize {
        self.geometries.len()
    }
}

impl GeometryStorage for InMemory {}
