use std::fmt;

use crate::algorithms::{
    add_hydrogens_to_molecule, remove_hydrogens_from_molecule, AddHydrogensOptions,
    AddHydrogensReport, HydrogenNormalizationError, RemoveHydrogensReport,
};
use crate::chemistry::{
    sanitize_small_molecule, sanitize_small_molecule_with_ring_options, SanitizeError,
    SanitizeOptions, SanitizeReport,
};
use crate::io::{
    interpret_smiles_document, parse_smiles_document, write_canonical_smiles,
    write_isomeric_smiles, write_smiles, CanonicalSmilesWriteOptions, IsomericSmilesWriteOptions,
    MolWriteError, SmilesInterpretError, SmilesParseError, SmilesWriteOptions,
};

use super::model::SmallMolecule;

impl SmallMolecule {
    pub fn from_smiles(input: &str) -> Result<Self, SmallMoleculeReadError> {
        let document = parse_smiles_document(input)?;
        Ok(interpret_smiles_document(&document)?.into_molecule())
    }

    pub fn from_smiles_sanitized(input: &str) -> Result<Self, SmallMoleculeReadError> {
        let mut molecule = Self::from_smiles(input)?;
        molecule.sanitize()?;
        Ok(molecule)
    }

    pub fn sanitize(&mut self) -> Result<SanitizeReport, SanitizeError> {
        sanitize_small_molecule(self, SanitizeOptions::default())
    }

    pub fn sanitize_with_options(
        &mut self,
        options: SanitizeOptions,
    ) -> Result<SanitizeReport, SanitizeError> {
        sanitize_small_molecule_with_ring_options(
            self,
            options,
            crate::algorithms::RingPerceptionOptions::default(),
        )
    }

    /// Materialize stored and perceived hydrogens as graph atoms.
    pub fn add_hydrogens(&mut self) -> Result<AddHydrogensReport, HydrogenNormalizationError> {
        self.add_hydrogens_with_options(AddHydrogensOptions::default())
    }

    /// Materialize hydrogens under the supplied count and growth policy.
    pub fn add_hydrogens_with_options(
        &mut self,
        options: AddHydrogensOptions,
    ) -> Result<AddHydrogensReport, HydrogenNormalizationError> {
        add_hydrogens_to_molecule(self.graph_mut_raw(), options)
    }

    /// Collapse ordinary graph hydrogens and report retained protected atoms.
    pub fn remove_hydrogens(
        &mut self,
    ) -> Result<RemoveHydrogensReport, HydrogenNormalizationError> {
        remove_hydrogens_from_molecule(self.graph_mut_raw())
    }

    pub fn to_smiles(&self) -> Result<String, MolWriteError> {
        write_smiles(self, SmilesWriteOptions)
    }

    pub fn to_isomeric_smiles(&self) -> Result<String, MolWriteError> {
        write_isomeric_smiles(self, IsomericSmilesWriteOptions)
    }

    pub fn to_canonical_smiles(&self) -> Result<String, MolWriteError> {
        write_canonical_smiles(self, CanonicalSmilesWriteOptions)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SmallMoleculeReadError {
    Parse(SmilesParseError),
    Interpret(SmilesInterpretError),
    Sanitize(SanitizeError),
}

impl fmt::Display for SmallMoleculeReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(error) => write!(f, "{error}"),
            Self::Interpret(error) => write!(f, "{error}"),
            Self::Sanitize(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for SmallMoleculeReadError {}

impl From<SmilesParseError> for SmallMoleculeReadError {
    fn from(error: SmilesParseError) -> Self {
        Self::Parse(error)
    }
}

impl From<SmilesInterpretError> for SmallMoleculeReadError {
    fn from(error: SmilesInterpretError) -> Self {
        Self::Interpret(error)
    }
}

impl From<SanitizeError> for SmallMoleculeReadError {
    fn from(error: SanitizeError) -> Self {
        Self::Sanitize(error)
    }
}
