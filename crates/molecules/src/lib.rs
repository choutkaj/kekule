#![forbid(unsafe_code)]

mod algorithms;
pub mod bio;
mod chemistry;
pub mod core;
mod io;
pub mod small;

pub mod smiles {
    pub use crate::io::{
        CanonicalSmilesWriteOptions, MolWriteError, SmilesParseError, SmilesParseOptions,
        SmilesWriteOptions,
    };
    pub use crate::small::SmallMoleculeReadError;

    use crate::small::{SanitizeOptions, SmallMolecule};

    pub fn read_str(input: &str) -> Result<SmallMolecule, SmilesParseError> {
        crate::io::read_smiles_str(input, SmilesParseOptions)
    }

    pub fn read_str_with_options(
        input: &str,
        options: SmilesParseOptions,
    ) -> Result<SmallMolecule, SmilesParseError> {
        crate::io::read_smiles_str(input, options)
    }

    pub fn read_sanitized_str(input: &str) -> Result<SmallMolecule, SmallMoleculeReadError> {
        let mut molecule = read_str(input)?;
        crate::perception::sanitize_with_options(&mut molecule, SanitizeOptions::default())?;
        Ok(molecule)
    }

    pub fn write(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_smiles(molecule, SmilesWriteOptions)
    }

    pub fn write_with_options(
        molecule: &SmallMolecule,
        options: SmilesWriteOptions,
    ) -> Result<String, MolWriteError> {
        crate::io::write_smiles(molecule, options)
    }

    pub fn write_canonical(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_canonical_smiles(molecule, CanonicalSmilesWriteOptions)
    }

    pub fn write_canonical_with_options(
        molecule: &SmallMolecule,
        options: CanonicalSmilesWriteOptions,
    ) -> Result<String, MolWriteError> {
        crate::io::write_canonical_smiles(molecule, options)
    }
}

pub mod molfile {
    pub use crate::io::{MolParseOptions, MolWriteError, SdfParseError};

    use crate::small::SmallMolecule;

    pub fn read_v2000_str(input: &str) -> Result<SmallMolecule, SdfParseError> {
        crate::io::read_mol_v2000_str(input)
    }

    pub fn write_v2000(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_mol_v2000(molecule)
    }

    pub fn read_v3000_str(input: &str) -> Result<SmallMolecule, SdfParseError> {
        crate::io::read_mol_v3000_str(input)
    }

    pub fn write_v3000(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_mol_v3000(molecule)
    }
}

pub mod sdf {
    pub use crate::io::{MolWriteError, SdfParseError, SdfParseOptions, SdfRecord};

    use crate::small::SmallMolecule;

    pub fn read_v2000_str(
        input: &str,
        options: SdfParseOptions,
    ) -> Result<Vec<SmallMolecule>, SdfParseError> {
        crate::io::read_sdf_v2000_str(input, options)
    }

    pub fn read_v2000_records(
        input: &str,
        options: SdfParseOptions,
    ) -> Result<Vec<SdfRecord>, SdfParseError> {
        crate::io::read_sdf_v2000_records(input, options)
    }

    pub fn write_v2000(molecules: &[SmallMolecule]) -> Result<String, MolWriteError> {
        crate::io::write_sdf_v2000(molecules)
    }
}

pub mod perception {
    pub use crate::algorithms::{
        perceive_aromaticity, perceive_aromaticity_with_ring_options, perceive_ring_membership,
        perceive_ring_set, perceive_ring_set_with_options, perceive_valence, AromaticityError,
        AromaticityModel, Ring, RingMembership, RingPerceptionError, RingPerceptionOptions,
        RingSet, RingWork, ValenceIssue, ValenceModel, ValenceReport,
    };
    pub use crate::chemistry::{SanitizeError, SanitizeOptions, SanitizeReport};

    use crate::small::SmallMolecule;

    pub fn sanitize(molecule: &mut SmallMolecule) -> Result<SanitizeReport, SanitizeError> {
        crate::chemistry::sanitize_small_molecule(molecule, SanitizeOptions::default())
    }

    pub fn sanitize_with_options(
        molecule: &mut SmallMolecule,
        options: SanitizeOptions,
    ) -> Result<SanitizeReport, SanitizeError> {
        crate::chemistry::sanitize_small_molecule(molecule, options)
    }

    pub fn sanitize_with_ring_options(
        molecule: &mut SmallMolecule,
        options: SanitizeOptions,
        ring_options: RingPerceptionOptions,
    ) -> Result<SanitizeReport, SanitizeError> {
        crate::chemistry::sanitize_small_molecule_with_ring_options(molecule, options, ring_options)
    }
}

pub mod canon {
    pub use crate::algorithms::CanonicalAtomRanking;

    use crate::core::Molecule;

    pub fn atom_ranking(molecule: &Molecule) -> CanonicalAtomRanking {
        crate::algorithms::canonical_atom_ranking(molecule)
    }
}

pub mod prelude {
    pub use crate::bio::{BioHierarchy, MacroMolecule};
    pub use crate::core::{Atom, AtomId, Bond, BondId, BondOrder, Conformer, Element, Molecule};
    pub use crate::small::{SanitizeOptions, SanitizeReport, SmallMolecule};
    pub use crate::smiles::{CanonicalSmilesWriteOptions, SmilesParseOptions, SmilesWriteOptions};
}

#[cfg(test)]
mod tests;
