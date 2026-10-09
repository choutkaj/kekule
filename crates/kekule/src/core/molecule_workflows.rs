use crate::algorithms::{
    add_hydrogens_to_molecule, remove_hydrogens_from_molecule, AddHydrogensOptions,
    AddHydrogensReport, HydrogenTransformError, RemoveHydrogensReport,
};
use crate::chemistry::{perceive_molecule, PerceptionError};

use super::Molecule;

impl Molecule {
    /// Install default valence, rings, aromaticity and conjugation transactionally.
    ///
    /// This derives perception from the canonical represented chemistry and
    /// never rewrites atoms, bonds, or represented stereochemistry.
    ///
    /// Interpretation publishes canonical represented chemistry directly;
    /// there is no separate public normalization step before perception.
    pub fn perceive(&mut self) -> Result<(), PerceptionError> {
        perceive_molecule(self)
    }

    /// Convert resolved implicit hydrogens into explicit graph atoms.
    /// Fixed counts require no perception; inference-enabled counts must be
    /// perceived first. Success invalidates perception but preserves composition.
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
    /// Parents that permit inference require an installed hydrogen count; fixed
    /// declarations need no perception. A graph with nothing removable is unchanged.
    /// The result has dense IDs, so atoms after a removed hydrogen are renumbered;
    /// [`RemoveHydrogensReport::correspondence`] translates the input IDs.
    pub fn remove_hydrogens(&mut self) -> Result<RemoveHydrogensReport, HydrogenTransformError> {
        remove_hydrogens_from_molecule(self)
    }
}
