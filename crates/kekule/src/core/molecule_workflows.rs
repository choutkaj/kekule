use std::fmt;

use crate::algorithms::{
    add_hydrogens_to_molecule, remove_hydrogens_from_molecule, AddHydrogensOptions,
    AddHydrogensReport, HydrogenTransformError, RemoveHydrogensReport,
};
use crate::chemistry::{perceive_molecule, PerceptionError};
use crate::io::{
    interpret_smiles_document, parse_smiles_document, write_canonical_smiles,
    write_isomeric_smiles, write_smiles, MolWriteError, SmilesInterpretError, SmilesParseError,
};

use super::Molecule;

impl Molecule {
    /// Parses and interprets one SMILES record into source-ordered connected components.
    pub fn from_smiles(input: &str) -> Result<Vec<Self>, MoleculeReadError> {
        let document = parse_smiles_document(input)?;
        Ok(interpret_smiles_document(&document)?.to_molecules())
    }

    /// Install the transactional default valence, ring-set, and aromaticity profile.
    ///
    /// This derives perception from the canonical represented chemistry and
    /// never rewrites atoms, bonds, or represented stereochemistry.
    ///
    /// There is no public normalization step between interpretation and
    /// perception:
    ///
    pub fn perceive(&mut self) -> Result<(), PerceptionError> {
        perceive_molecule(self)
    }

    /// Materialize stored and perceived hydrogens as graph atoms.
    pub fn add_hydrogens(&mut self) -> Result<AddHydrogensReport, HydrogenTransformError> {
        self.add_hydrogens_with_options(AddHydrogensOptions::default())
    }

    /// Materialize hydrogens under the supplied count and growth policy.
    pub fn add_hydrogens_with_options(
        &mut self,
        options: AddHydrogensOptions,
    ) -> Result<AddHydrogensReport, HydrogenTransformError> {
        add_hydrogens_to_molecule(self, options)
    }

    /// Collapse ordinary graph hydrogens and report retained protected atoms.
    pub fn remove_hydrogens(&mut self) -> Result<RemoveHydrogensReport, HydrogenTransformError> {
        remove_hydrogens_from_molecule(self)
    }

    pub fn to_smiles(&self) -> Result<String, MolWriteError> {
        write_smiles(self)
    }

    pub fn to_isomeric_smiles(&self) -> Result<String, MolWriteError> {
        write_isomeric_smiles(self)
    }

    pub fn to_canonical_smiles(&self) -> Result<String, MolWriteError> {
        write_canonical_smiles(self)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum MoleculeReadError {
    Parse(SmilesParseError),
    Interpret(SmilesInterpretError),
}

impl fmt::Display for MoleculeReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(error) => write!(f, "{error}"),
            Self::Interpret(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for MoleculeReadError {}

impl From<SmilesParseError> for MoleculeReadError {
    fn from(error: SmilesParseError) -> Self {
        Self::Parse(error)
    }
}

impl From<SmilesInterpretError> for MoleculeReadError {
    fn from(error: SmilesInterpretError) -> Self {
        Self::Interpret(error)
    }
}
