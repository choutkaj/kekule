//! Shared definitions for identical molecules of one interpreted structure.
//!
//! A structure file lists every occurrence of a molecule, such as each water
//! of a crystal or a solvent box, as atoms of its own. Interpretation stages
//! one definition per distinct molecule and an instance per occurrence, so
//! occurrences map atom for atom and per-definition work, such as perception,
//! runs once per distinct molecule.
use std::collections::HashMap;
use std::hash::{BuildHasher, Hash, Hasher};

use crate::molecule::Molecule;
use crate::structure::{StructureBuildError, StructureBuilder};
use crate::topology::MoleculeDefinitionId;

/// Stages one definition per distinct molecule of a model under construction.
///
/// Two molecules share a definition when they are identical in local atom
/// order: equal atoms (element, isotope, formal charge, radical, implicit
/// hydrogens, atom map), equal bonds between the same local atoms with the
/// same orders, equal represented stereo, equal molecule properties, and
/// equal per-atom labels. Labels carry what the format names each atom,
/// such as its component and atom name, which classification reads from the
/// hierarchy; requiring them equal keeps one class valid for every instance
/// of a definition. The labels themselves stay hierarchy data of each
/// instance.
pub(crate) struct DefinitionReuse<L> {
    hasher: std::collections::hash_map::RandomState,
    staged: HashMap<u64, Vec<(MoleculeDefinitionId, Vec<L>)>>,
}

impl<L: Hash + Eq> DefinitionReuse<L> {
    pub(crate) fn new() -> Self {
        Self {
            hasher: std::collections::hash_map::RandomState::new(),
            staged: HashMap::new(),
        }
    }

    /// The definition for `molecule`, whose atom `i` is labeled `labels[i]`:
    /// a staged definition identical to it, or a new one.
    pub(crate) fn definition(
        &mut self,
        builder: &mut StructureBuilder,
        molecule: Molecule,
        labels: Vec<L>,
    ) -> Result<MoleculeDefinitionId, StructureBuildError> {
        debug_assert_eq!(labels.len(), molecule.atom_count());
        let key = self.fingerprint(&molecule, &labels);
        let candidates = self.staged.entry(key).or_default();
        for (definition, staged_labels) in candidates.iter() {
            let staged = builder
                .topology_builder()
                .definition(*definition)
                .expect("reused definitions are staged")
                .molecule();
            if *staged_labels == labels
                && *staged == molecule
                && staged.properties == molecule.properties
            {
                return Ok(*definition);
            }
        }
        let definition = builder.add_molecule_definition(molecule)?;
        candidates.push((definition, labels));
        Ok(definition)
    }

    /// A hash of everything [`Self::definition`] compares, except stereo
    /// and properties, which only the full comparison reads.
    fn fingerprint(&self, molecule: &Molecule, labels: &[L]) -> u64 {
        let mut hasher = self.hasher.build_hasher();
        molecule.hash_graph(&mut hasher);
        labels.hash(&mut hasher);
        hasher.finish()
    }
}

#[cfg(test)]
mod tests {
    use super::DefinitionReuse;
    use crate::molecule::Molecule;
    use crate::structure::StructureBuilder;

    fn molecule(smiles: &str) -> Molecule {
        crate::smiles::to_molecules(smiles).unwrap().remove(0)
    }

    #[test]
    fn only_identical_molecules_with_equal_labels_share_a_definition() {
        let mut builder = StructureBuilder::new();
        let mut reuse = DefinitionReuse::new();
        let mut definition = |smiles: &str, label: &str| {
            let molecule = molecule(smiles);
            let labels = vec![label.to_owned(); molecule.atom_count()];
            reuse.definition(&mut builder, molecule, labels).unwrap()
        };
        let ethanol = definition("CCO", "EOH");
        assert_eq!(definition("CCO", "EOH"), ethanol);
        // Each difference gives a definition of its own: atom order, charge,
        // isotope, stated hydrogens, atom map, bond order, and labels.
        let different = [
            ("OCC", "EOH"),
            ("CC[O-]", "EOH"),
            ("CC[18OH]", "EOH"),
            ("CC[OH]", "EOH"),
            ("[CH3:1]CO", "EOH"),
            ("C=CO", "EOH"),
            ("CCO", "ETH"),
        ]
        .map(|(smiles, label)| definition(smiles, label));
        let mut all = different.to_vec();
        all.push(ethanol);
        all.sort();
        all.dedup();
        assert_eq!(all.len(), different.len() + 1);

        // Stereo is chemistry: enantiomers do not share, equal stereo does.
        let left = definition("C[C@H](N)O", "ALO");
        let right = definition("C[C@@H](N)O", "ALO");
        assert_ne!(left, right);
        assert_eq!(definition("C[C@H](N)O", "ALO"), left);
    }
}
