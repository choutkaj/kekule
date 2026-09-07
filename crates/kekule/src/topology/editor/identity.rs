use std::fmt;
use std::sync::atomic::{AtomicU64, Ordering};

// Cloneable drafts share existing handles, but allocate distinct new handles.
// Keeping allocation global prevents aliasing across independent or cloned drafts.
static NEXT_HANDLE: AtomicU64 = AtomicU64::new(1);
fn fresh() -> u64 {
    NEXT_HANDLE
        .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |id| id.checked_add(1))
        .expect("editing handle space exhausted")
}

macro_rules! handle {
    ($name:ident, $label:literal) => {
        #[doc = concat!("Opaque ", $label, " identity in a structural editing draft.\n\nHandles are stable through splits and merges; deleted or foreign handles are rejected.")]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u64);
        impl $name { pub(super) fn new() -> Self { Self(fresh()) } }
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, concat!($label, "{}"), self.0) }
        }
    };
}
handle!(EditAtomId, "edit-atom");
handle!(EditBondId, "edit-bond");
handle!(EditChainId, "edit-chain");
handle!(EditResidueId, "edit-residue");
handle!(EditAtomSiteId, "edit-site");
