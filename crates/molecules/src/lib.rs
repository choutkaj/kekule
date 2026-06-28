#![forbid(unsafe_code)]

mod algorithms;
mod bio;
mod chemistry;
mod core;
mod io;

pub use algorithms::*;
pub use bio::*;
pub use chemistry::*;
pub use core::*;
pub use io::*;

pub mod prelude {
    pub use crate::{
        perceive_aromaticity, perceive_aromaticity_with_ring_options, perceive_ring_membership,
        perceive_ring_set, perceive_ring_set_with_options, perceive_valence, read_mmcif_str,
        read_mol_v2000_str, read_mol_v3000_str, read_sdf_v2000_str, read_smiles_str,
        sanitize_small_molecule, sanitize_small_molecule_with_ring_options, write_mol_v2000,
        write_mol_v3000, write_sdf_v2000, write_smiles, AromaticityError, AromaticityModel, Atom,
        AtomId, AtomMut, AtomRadical, AtomSite, AtomSiteId, AtomSiteMetadata, AtomStereo,
        BioHierarchy, BioHierarchyError, Bond, BondId, BondMut, BondOrder, BondStereo, Chain,
        ChainId, ComputedState, Conformer, ConformerId, Element, MacroMolecule, MmcifParseError,
        MmcifParseOptions, Model, ModelId, MolWriteError, Molecule, MoleculeError, Point3, PropMap,
        PropValue, Residue, ResidueId, Result, Ring, RingMembership, RingPerceptionError,
        RingPerceptionOptions, RingSet, RingWork, SanitizeError, SanitizeOptions, SanitizeReport,
        SdfParseError, SdfParseOptions, SdfRecord, SmallMolecule, SmilesParseError,
        SmilesParseOptions, SmilesWriteOptions, ValenceIssue, ValenceModel, ValenceReport,
    };
}

#[cfg(test)]
mod tests;
