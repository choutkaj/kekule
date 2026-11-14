//! The crate's `u32` identifier types: the `fixed_u32_id!` macro that defines
//! them, the [`FixedId`] trait generic code builds them through, and checked
//! ID allocation.

use crate::diagnostics::{CapacityExceeded, Countable};

/// An identifier defined by `fixed_u32_id!`: its inherent `new`, `raw`, and
/// `index` for generic code, and `from_index` for the crate.
pub(crate) trait FixedId: Copy + Ord {
    fn new(raw: u32) -> Self;
    fn raw(self) -> u32;
    fn index(self) -> usize {
        self.raw() as usize
    }
    /// The ID at a position its owner already holds; allocation checked it
    /// with [`checked_id`] or [`checked_len`].
    fn from_index(index: usize) -> Self {
        Self::new(u32::try_from(index).expect("held positions fit their identifier width"))
    }
}

macro_rules! fixed_u32_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
            /// Wraps a raw value. Nothing is validated here; a lookup with a
            /// value its owner does not hold returns `None`.
            pub const fn new(raw: u32) -> Self {
                Self(raw)
            }

            /// Returns the raw value. It is this identifier's stable
            /// encoding for persistence: `new(id.raw())` restores `id`, and
            /// the encoding does not change between Kekule versions.
            pub const fn raw(self) -> u32 {
                self.0
            }

            /// Returns the raw value as a `usize` position.
            pub const fn index(self) -> usize {
                self.0 as usize
            }
        }

        impl crate::fixed_u32_id::FixedId for $name {
            fn new(raw: u32) -> Self {
                Self(raw)
            }
            fn raw(self) -> u32 {
                self.0
            }
        }
    };
    ($(#[$meta:meta])* $name:ident, $display:literal) => {
        fixed_u32_id!($(#[$meta])* $name);

        impl std::fmt::Display for $name {
            fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(formatter, concat!($display, " {}"), self.0)
            }
        }
    };
}

/// The ID at position `index`, which a collection of `index` items allocates
/// next.
pub(crate) fn checked_id<T: FixedId>(index: usize, kind: Countable) -> Result<T, CapacityExceeded> {
    u32::try_from(index)
        .map(T::new)
        .map_err(|_| CapacityExceeded::new(kind))
}

/// Whether every item of a collection of `len` items has an ID.
pub(crate) fn ids_fit(len: usize) -> bool {
    len.checked_sub(1)
        .is_none_or(|last| u32::try_from(last).is_ok())
}

/// Checks that `additional` more items after `current` still all have IDs.
pub(crate) fn checked_len(
    current: usize,
    additional: usize,
    kind: Countable,
) -> Result<(), CapacityExceeded> {
    match current.checked_add(additional) {
        Some(len) if ids_fit(len) => Ok(()),
        _ => Err(CapacityExceeded::new(kind)),
    }
}

#[cfg(all(test, target_pointer_width = "64"))]
mod tests {
    use super::*;
    use crate::molecule::{AtomId, BondId, StereoElementId, StereoGroupId};
    use crate::query::{QueryAtomId, QueryBondId};
    use crate::topology::{
        AtomSiteId, ChainId, MoleculeDefinitionId, MoleculeInstanceId, ResidueId,
        TopologyAtomIndex, TopologyBondIndex,
    };

    const MAX_SLOT: usize = u32::MAX as usize;
    const FIRST_UNSUPPORTED_SLOT: usize = MAX_SLOT + 1;

    fn check_boundary<T: FixedId + std::fmt::Debug>(kind: Countable) {
        assert_eq!(checked_id::<T>(MAX_SLOT, kind), Ok(T::new(u32::MAX)));
        assert_eq!(
            checked_id::<T>(FIRST_UNSUPPORTED_SLOT, kind),
            Err(CapacityExceeded::new(kind))
        );
    }

    #[test]
    fn every_identifier_space_checks_its_boundary() {
        check_boundary::<AtomId>(Countable::Atom);
        check_boundary::<BondId>(Countable::Bond);
        check_boundary::<StereoElementId>(Countable::StereoElement);
        check_boundary::<StereoGroupId>(Countable::StereoGroup);
        check_boundary::<QueryAtomId>(Countable::Atom);
        check_boundary::<QueryBondId>(Countable::Bond);
        check_boundary::<MoleculeDefinitionId>(Countable::MoleculeDefinition);
        check_boundary::<MoleculeInstanceId>(Countable::MoleculeInstance);
        check_boundary::<TopologyAtomIndex>(Countable::Atom);
        check_boundary::<TopologyBondIndex>(Countable::Bond);
        check_boundary::<ChainId>(Countable::Chain);
        check_boundary::<ResidueId>(Countable::Residue);
        check_boundary::<AtomSiteId>(Countable::AtomSite);
    }

    #[test]
    fn collection_lengths_stop_at_the_identifier_space() {
        assert_eq!(checked_len(0, 0, Countable::Atom), Ok(()));
        assert_eq!(checked_len(MAX_SLOT, 1, Countable::Atom), Ok(()));
        assert_eq!(
            checked_len(FIRST_UNSUPPORTED_SLOT, 1, Countable::Atom),
            Err(CapacityExceeded::new(Countable::Atom))
        );
        assert_eq!(
            checked_len(usize::MAX, 1, Countable::Atom),
            Err(CapacityExceeded::new(Countable::Atom))
        );
    }
}
