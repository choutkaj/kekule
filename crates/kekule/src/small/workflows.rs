use std::fmt;

use crate::algorithms::{
    add_hydrogens_to_molecule, remove_hydrogens_from_molecule, AddHydrogensOptions,
    AddHydrogensReport, HydrogenTransformError, RemoveHydrogensReport,
};
use crate::chemistry::{
    normalize_molecule, perceive_molecule, NormalizationError, NormalizationReport, PerceptionError,
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
        interpret_smiles_document(&document)?
            .into_molecule()
            .map_err(|error| SmallMoleculeReadError::ComponentCount {
                actual: error.actual(),
            })
    }

    /// Normalize represented chemistry and source stereo into canonical form.
    pub fn normalize(&mut self) -> Result<NormalizationReport, NormalizationError> {
        normalize_molecule(self.graph_mut_raw())
    }

    /// Install the transactional default valence, ring-set, and aromaticity profile.
    ///
    /// Call [`Self::normalize`] first when the molecule came from source
    /// representation that has not yet been normalized.
    pub fn perceive(&mut self) -> Result<(), PerceptionError> {
        perceive_molecule(self.graph_mut_raw())
    }

    /// Normalize represented chemistry and install the default perception profile.
    ///
    /// The complete composition is transactional. Coordinate-derived stereo
    /// and CIP assignment remain separate explicit operations.
    pub fn normalize_and_perceive(
        &mut self,
    ) -> Result<NormalizationReport, NormalizeAndPerceiveError> {
        let mut staged = self.clone();
        let report = normalize_molecule(staged.graph_mut_raw())
            .map_err(NormalizeAndPerceiveError::Normalization)?;
        perceive_molecule(staged.graph_mut_raw()).map_err(NormalizeAndPerceiveError::Perception)?;
        *self = staged;
        Ok(report)
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
        add_hydrogens_to_molecule(self.graph_mut_raw(), options)
    }

    /// Collapse ordinary graph hydrogens and report retained protected atoms.
    pub fn remove_hydrogens(&mut self) -> Result<RemoveHydrogensReport, HydrogenTransformError> {
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
    /// A `SmallMolecule` convenience requires one connected SMILES component.
    ComponentCount {
        actual: usize,
    },
}

impl fmt::Display for SmallMoleculeReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(error) => write!(f, "{error}"),
            Self::Interpret(error) => write!(f, "{error}"),
            Self::ComponentCount { actual } => write!(
                f,
                "SmallMolecule requires exactly one connected SMILES component, found {actual}"
            ),
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

/// Failure from the ordinary normalization plus default-perception workflow.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum NormalizeAndPerceiveError {
    Normalization(NormalizationError),
    Perception(PerceptionError),
}

impl fmt::Display for NormalizeAndPerceiveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Normalization(error) => write!(f, "{error}"),
            Self::Perception(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for NormalizeAndPerceiveError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Normalization(error) => Some(error),
            Self::Perception(error) => Some(error),
        }
    }
}
