use std::fmt;

use crate::core::{Atom, AtomId, Bond, BondId, Molecule};
use crate::io::{MolWriteError, SmilesInterpretError, SmilesParseError};
use crate::{perception, smiles};

pub use crate::chemistry::{SanitizeError, SanitizeOptions, SanitizeReport};

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SmallMolecule {
    pub(crate) graph: Molecule,
    pub(crate) data: SmallMoleculeData,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct SmallMoleculeData;

impl SmallMolecule {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_graph(graph: Molecule) -> Self {
        Self {
            graph,
            data: SmallMoleculeData,
        }
    }

    pub fn into_graph(self) -> Molecule {
        self.graph
    }

    pub fn graph(&self) -> &Molecule {
        &self.graph
    }

    pub fn graph_mut(&mut self) -> &mut Molecule {
        &mut self.graph
    }

    pub(crate) fn graph_mut_raw(&mut self) -> &mut Molecule {
        &mut self.graph
    }

    pub(crate) fn without_conformers(mut self) -> Self {
        self.graph = self.graph.without_conformers();
        self
    }

    pub fn from_smiles(input: &str) -> Result<Self, SmallMoleculeReadError> {
        let document = smiles::parse_str(input)?;
        Ok(smiles::interpret(&document)?)
    }

    pub fn from_smiles_sanitized(input: &str) -> Result<Self, SmallMoleculeReadError> {
        let mut molecule = Self::from_smiles(input)?;
        molecule.sanitize()?;
        Ok(molecule)
    }

    pub fn sanitize(&mut self) -> Result<SanitizeReport, SanitizeError> {
        perception::sanitize(self)
    }

    pub fn sanitize_with_options(
        &mut self,
        options: SanitizeOptions,
    ) -> Result<SanitizeReport, SanitizeError> {
        perception::sanitize_with_options(self, options)
    }

    pub fn atom_count(&self) -> usize {
        self.graph.atom_count()
    }

    pub fn bond_count(&self) -> usize {
        self.graph.bond_count()
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.graph.atoms()
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)> {
        self.graph.bonds()
    }

    pub fn to_smiles(&self) -> Result<String, MolWriteError> {
        smiles::write(self)
    }

    pub fn to_isomeric_smiles(&self) -> Result<String, MolWriteError> {
        smiles::write_isomeric(self)
    }

    pub fn to_canonical_smiles(&self) -> Result<String, MolWriteError> {
        smiles::write_canonical(self)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
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
