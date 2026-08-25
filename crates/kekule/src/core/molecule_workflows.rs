use crate::algorithms::{
    add_hydrogens_to_molecule, remove_hydrogens_from_molecule, AddHydrogensOptions,
    AddHydrogensReport, HydrogenTransformError, RemoveHydrogensReport,
};
use crate::chemistry::{perceive_molecule, PerceptionError};

use super::Molecule;

impl Molecule {
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
}
