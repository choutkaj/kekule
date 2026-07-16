use crate::core::{Atom, AtomId, Bond, BondId, Molecule};

/// The small-molecule domain wrapper around one raw molecular graph.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SmallMolecule {
    pub(crate) graph: Molecule,
}

impl SmallMolecule {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_graph(graph: Molecule) -> Self {
        Self { graph }
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
}
