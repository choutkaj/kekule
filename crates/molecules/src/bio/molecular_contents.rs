use crate::small::SmallMolecule;

use super::MacroMolecule;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Solvent {
    molecules: Vec<SmallMolecule>,
}

impl Solvent {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.molecules.len()
    }

    pub fn is_empty(&self) -> bool {
        self.molecules.is_empty()
    }

    pub fn molecules(&self) -> impl Iterator<Item = &SmallMolecule> {
        self.molecules.iter()
    }

    pub fn into_molecules(self) -> Vec<SmallMolecule> {
        self.molecules
    }

    pub(crate) fn push(&mut self, molecule: SmallMolecule) {
        self.molecules.push(molecule);
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MolecularContents {
    small_molecules: Vec<SmallMolecule>,
    macromolecules: Vec<MacroMolecule>,
    solvent: Solvent,
}

impl MolecularContents {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn small_molecules(&self) -> impl Iterator<Item = &SmallMolecule> {
        self.small_molecules.iter()
    }

    pub fn macromolecules(&self) -> impl Iterator<Item = &MacroMolecule> {
        self.macromolecules.iter()
    }

    pub fn solvent(&self) -> &Solvent {
        &self.solvent
    }

    pub fn into_parts(self) -> (Vec<SmallMolecule>, Vec<MacroMolecule>, Solvent) {
        (self.small_molecules, self.macromolecules, self.solvent)
    }

    pub(crate) fn push_small(&mut self, molecule: SmallMolecule) {
        self.small_molecules.push(molecule);
    }

    pub(crate) fn push_macro(&mut self, molecule: MacroMolecule) {
        self.macromolecules.push(molecule);
    }

    pub(crate) fn solvent_mut(&mut self) -> &mut Solvent {
        &mut self.solvent
    }
}
