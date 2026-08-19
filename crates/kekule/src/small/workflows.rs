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

    /// Install the transactional default valence, ring-set, and aromaticity profile.
    ///
    /// This derives perception from the canonical represented chemistry and
    /// never rewrites atoms, bonds, or represented stereochemistry.
    ///
    /// There is no public normalization step between interpretation and
    /// perception:
    ///
    /// ```compile_fail
    /// use kekule::small::SmallMolecule;
    /// let mut molecule = SmallMolecule::from_smiles("CC").unwrap();
    /// molecule.normalize().unwrap();
    /// ```
    ///
    /// ```compile_fail
    /// let mut molecule = kekule::core::Molecule::new();
    /// kekule::normalization::normalize(&mut molecule).unwrap();
    /// ```
    pub fn perceive(&mut self) -> Result<(), PerceptionError> {
        perceive_molecule(self.graph_mut())
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
        add_hydrogens_to_molecule(self.graph_mut(), options)
    }

    /// Collapse ordinary graph hydrogens and report retained protected atoms.
    pub fn remove_hydrogens(&mut self) -> Result<RemoveHydrogensReport, HydrogenTransformError> {
        remove_hydrogens_from_molecule(self.graph_mut())
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
