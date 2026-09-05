use std::fmt;

use crate::chemistry::PerceptionError;

use super::{MoleculeDefinitionId, Topology};

/// Default perception failed for one reusable molecule definition.
///
/// The definition ID belongs to the source topology. Every source owner and
/// its installed perception remain unchanged on failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologyPerceptionError {
    pub definition: MoleculeDefinitionId,
    pub source: PerceptionError,
}

impl fmt::Display for TopologyPerceptionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "perception failed for molecule definition {}: {}",
            self.definition, self.source
        )
    }
}

impl std::error::Error for TopologyPerceptionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

impl Topology {
    /// Returns a new snapshot with the default molecular perception profile.
    ///
    /// Each reusable definition is perceived once, in definition order, using
    /// [`crate::core::Molecule::perceive`]: RDKit-like valence and implicit
    /// hydrogens, the default deterministic ring set, and RDKit-like aromaticity.
    /// Existing perception is recomputed; dependent CIP state is cleared.
    /// Failure identifies the first failing definition and leaves `self` intact.
    ///
    /// Represented graphs, definition reuse, all semantic IDs, dense ordering,
    /// hierarchy, classifications, and properties are preserved exactly. Thus
    /// [`Self::same_layout`] remains true, but the result is a distinct snapshot.
    /// Selections and prepared calculations bound to the source allocation stay
    /// bound to that source; this operation does not transfer their bindings.
    ///
    /// ```
    /// use kekule::{smiles, topology::Topology};
    ///
    /// let molecules = smiles::to_molecules("c1ccccc1.[Na+]")?;
    /// let source = Topology::from_molecules(&molecules)?;
    /// let perceived = source.perceived()?;
    /// assert!(source.same_layout(&perceived));
    /// assert!(perceived.molecules().all(|m| m.molecule().perception().has_aromaticity()));
    /// assert!(source.molecules().all(|m| !m.molecule().perception().has_aromaticity()));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn perceived(&self) -> Result<Self, TopologyPerceptionError> {
        let definitions = self
            .definitions
            .iter()
            .map(|definition| {
                let mut perceived = definition.clone();
                perceived
                    .molecule
                    .perceive()
                    .map_err(|source| TopologyPerceptionError {
                        definition: definition.id,
                        source,
                    })?;
                Ok(perceived)
            })
            .collect::<Result<_, TopologyPerceptionError>>()?;

        // Copy the published layout directly: perception never requires
        // reclassification, reindexing, or rebuilding hierarchy correspondence.
        Ok(Self {
            definitions,
            instances: self.instances.clone(),
            instance_atoms: self.instance_atoms.clone(),
            instance_bonds: self.instance_bonds.clone(),
            atom_indices: self.atom_indices.clone(),
            bond_indices: self.bond_indices.clone(),
            hierarchy: self.hierarchy.clone(),
            properties: self.properties.clone(),
        })
    }
}
