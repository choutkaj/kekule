//! Atoms whose fixed hydrogen counts leave valence open.
use crate::molecule::{AtomId, ImplicitHydrogens, Molecule};

use super::rdkit_inferred_hydrogen_count;

/// An atom whose fixed implicit-hydrogen count is below what its bonds,
/// charge and declared radical electrons allow.
///
/// Such an atom has undeclared open valence: an edit or a file left it short
/// of hydrogens. Valence perception accepts it as stated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpenValence {
    atom: AtomId,
    missing: u8,
}

impl OpenValence {
    /// The atom.
    pub const fn atom(&self) -> AtomId {
        self.atom
    }

    /// How many hydrogens the atom lacks for its default valence.
    pub const fn missing(&self) -> u8 {
        self.missing
    }
}

/// The atoms of `molecule` with fixed counts that leave open valence, in atom
/// order. Inferred counts always fill an atom's valence, so only fixed counts
/// appear.
pub fn open_valences(molecule: &Molecule) -> Vec<OpenValence> {
    molecule
        .atoms()
        .filter_map(|(atom, value)| {
            let ImplicitHydrogens::Fixed(count) = value.hydrogens else {
                return None;
            };
            let mut inferred = value.clone();
            inferred.hydrogens = ImplicitHydrogens::Inferred;
            let full = rdkit_inferred_hydrogen_count(molecule.graph(), atom, &inferred);
            let missing = full.checked_sub(count).filter(|&missing| missing > 0)?;
            Some(OpenValence { atom, missing })
        })
        .collect()
}
