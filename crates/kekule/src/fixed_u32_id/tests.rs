use super::*;
use crate::molecule::{AtomId, BondId, StereoElementId, StereoGroupId};
use crate::query::{QueryAtomId, QueryBondId};
use crate::topology::{
    AtomSiteId, ChainId, MoleculeDefinitionId, MoleculeInstanceId, ResidueId, TopologyAtomIndex,
    TopologyBondIndex,
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
