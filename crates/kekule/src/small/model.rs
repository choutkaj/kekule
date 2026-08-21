use crate::core::{Atom, AtomId, Bond, BondId, Molecule, MoleculeEditor, Result};

/// The small-molecule domain wrapper around one connected molecular graph.
///
/// Empty and single-atom values retain the core graph's valid boundary-case
/// semantics.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SmallMolecule {
    pub(crate) molecule: Molecule,
}

impl SmallMolecule {
    pub fn new() -> Self {
        Self::default()
    }

    /// Wraps a completed connected core graph.
    pub fn from_molecule(molecule: Molecule) -> Self {
        debug_assert!(
            molecule.is_connected(),
            "completed molecule must be connected"
        );
        Self { molecule }
    }

    /// Wraps private interpretation staging before its graph is partitioned or
    /// rejected at the public format boundary.
    pub(crate) fn from_molecule_unchecked_connectedness(molecule: Molecule) -> Self {
        Self { molecule }
    }

    pub fn to_molecule(self) -> Molecule {
        self.molecule
    }

    pub fn as_molecule(&self) -> &Molecule {
        &self.molecule
    }

    pub fn as_molecule_mut(&mut self) -> &mut Molecule {
        &mut self.molecule
    }

    pub fn edit(&mut self) -> MoleculeEditor<'_> {
        self.molecule.edit()
    }

    pub(crate) fn without_conformers(mut self) -> Self {
        self.molecule = self.molecule.without_conformers();
        self
    }

    pub(crate) fn clone_without_conformers(&self) -> Self {
        Self {
            molecule: self.molecule.clone_without_conformers(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.molecule.is_empty()
    }

    pub fn atom_count(&self) -> usize {
        self.molecule.atom_count()
    }

    pub fn bond_count(&self) -> usize {
        self.molecule.bond_count()
    }

    pub fn formal_charge(&self) -> i64 {
        self.molecule.formal_charge()
    }

    pub fn atom(&self, id: AtomId) -> Result<&Atom> {
        self.molecule.atom(id)
    }

    pub fn bond(&self, id: BondId) -> Result<&Bond> {
        self.molecule.bond(id)
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.molecule.atoms()
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)> {
        self.molecule.bonds()
    }
}

impl AsRef<Molecule> for SmallMolecule {
    fn as_ref(&self) -> &Molecule {
        self.as_molecule()
    }
}
