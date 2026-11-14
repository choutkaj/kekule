//! Whole-entity expansions over a dense atom bit mask, the one implementation
//! behind [`AtomSelection`](super::AtomSelection)'s expansions and the
//! selection language's `by` terms.

use crate::fixed_u32_id::FixedId;
use crate::topology::{Topology, TopologyAtomIndex};

/// The entity an expansion completes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Expansion {
    /// Every atom placed in a touched residue.
    Residues,
    /// Every atom placed in a residue of a touched chain.
    Chains,
    /// Every atom of a touched molecule instance.
    Instances,
}

/// Adds to `mask`, laid out as
/// [`AtomSelection::to_mask`](super::AtomSelection::to_mask) over
/// `topology`, every atom of each entity at `level` that holds a set atom.
/// Set atoms stay set, including atoms without a hierarchy site.
pub(crate) fn expand_mask(topology: &Topology, mask: &mut [u64], level: Expansion) {
    match level {
        Expansion::Instances => {
            let atoms = topology.atom_ids();
            let mut instances = vec![false; topology.instance_count()];
            for index in ones(mask) {
                instances[atoms[index].instance().index()] = true;
            }
            for instance in topology.instances() {
                if instances[instance.id().index()] {
                    for atom in instance.atoms() {
                        set(mask, atom.index().index());
                    }
                }
            }
        }
        Expansion::Residues | Expansion::Chains => {
            let mut residues = vec![false; topology.residues().len()];
            let mut chains = vec![false; topology.chains().len()];
            for index in ones(mask) {
                if let Some(site) = topology.atom_site_at(TopologyAtomIndex::from_index(index)) {
                    let residue = site.residue();
                    residues[residue.id().index()] = true;
                    chains[residue.chain().id().index()] = true;
                }
            }
            let site_atoms = topology.site_atoms();
            for residue in topology.residues() {
                let selected = match level {
                    Expansion::Residues => residues[residue.id().index()],
                    _ => chains[residue.chain().id().index()],
                };
                if selected {
                    for site in residue.atom_sites() {
                        set(mask, site_atoms[site.id().index()].index());
                    }
                }
            }
        }
    }
}

fn set(mask: &mut [u64], index: usize) {
    mask[index / 64] |= 1 << (index % 64);
}

fn ones(mask: &[u64]) -> impl Iterator<Item = usize> + '_ {
    mask.iter().enumerate().flat_map(|(position, &word)| {
        let mut word = word;
        std::iter::from_fn(move || {
            (word != 0).then(|| {
                let bit = word.trailing_zeros() as usize;
                word &= word - 1;
                position * 64 + bit
            })
        })
    })
}
