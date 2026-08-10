use std::fmt;

use crate::core::{Atom, AtomId, BondId, BondOrder, MoleculeError};

use super::{MacroMoleculeBuilder, MacroMoleculeEditor};

/// Errors from public topology-changing operations on macro-molecule builders.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MacroGraphEditError {
    /// Additional atoms must be attached as part of the same operation so a
    /// public macro-molecule construction path never exposes a disconnected graph.
    DisconnectedAtomInsertion,
    Molecule(MoleculeError),
}

impl fmt::Display for MacroGraphEditError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DisconnectedAtomInsertion => formatter.write_str(
                "additional macro-molecule atoms must be inserted with add_atom_bonded_to",
            ),
            Self::Molecule(error) => write!(formatter, "{error}"),
        }
    }
}

impl std::error::Error for MacroGraphEditError {}

impl From<MoleculeError> for MacroGraphEditError {
    fn from(error: MoleculeError) -> Self {
        Self::Molecule(error)
    }
}

impl MacroMoleculeBuilder {
    /// Adds the first atom of a macro-molecule under construction.
    ///
    /// Once the graph is nonempty, use [`Self::add_atom_bonded_to`] so the
    /// public builder remains connected by construction.
    pub fn add_atom(&mut self, atom: Atom) -> Result<AtomId, MacroGraphEditError> {
        if self.graph().atom_count() != 0 {
            return Err(MacroGraphEditError::DisconnectedAtomInsertion);
        }
        Ok(self.graph_mut().add_atom(atom)?)
    }

    /// Adds one atom and its connecting bond atomically.
    pub fn add_atom_bonded_to(
        &mut self,
        parent: AtomId,
        atom: Atom,
        order: BondOrder,
    ) -> Result<(AtomId, BondId), MacroGraphEditError> {
        self.graph().atom(parent)?;
        let before = self.graph().clone();
        let atom_id = self.graph_mut().add_atom(atom)?;
        match self.graph_mut().add_bond(parent, atom_id, order) {
            Ok(bond_id) => Ok((atom_id, bond_id)),
            Err(error) => {
                *self.graph_mut() = before;
                Err(error.into())
            }
        }
    }

    /// Adds an additional bond between existing atoms.
    pub fn add_bond(
        &mut self,
        left: AtomId,
        right: AtomId,
        order: BondOrder,
    ) -> Result<BondId, MacroGraphEditError> {
        Ok(self.graph_mut().add_bond(left, right, order)?)
    }
}

impl MacroMoleculeEditor<'_> {
    /// Adds an atom to an empty macro-molecule editor.
    ///
    /// For a nonempty graph use [`Self::add_atom_bonded_to`].
    pub fn add_atom(&mut self, atom: Atom) -> Result<AtomId, MacroGraphEditError> {
        if self.graph().atom_count() != 0 {
            return Err(MacroGraphEditError::DisconnectedAtomInsertion);
        }
        Ok(self.graph_mut().add_atom(atom)?)
    }

    /// Adds one atom and its connecting bond atomically while preserving graph
    /// connectedness throughout the public editing operation.
    pub fn add_atom_bonded_to(
        &mut self,
        parent: AtomId,
        atom: Atom,
        order: BondOrder,
    ) -> Result<(AtomId, BondId), MacroGraphEditError> {
        self.graph().atom(parent)?;
        let before = self.graph().clone();
        let atom_id = self.graph_mut().add_atom(atom)?;
        match self.graph_mut().add_bond(parent, atom_id, order) {
            Ok(bond_id) => Ok((atom_id, bond_id)),
            Err(error) => {
                *self.graph_mut() = before;
                Err(error.into())
            }
        }
    }

    /// Adds an additional bond between existing atoms.
    pub fn add_bond(
        &mut self,
        left: AtomId,
        right: AtomId,
        order: BondOrder,
    ) -> Result<BondId, MacroGraphEditError> {
        Ok(self.graph_mut().add_bond(left, right, order)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Element;

    fn carbon() -> Atom {
        Atom::new(Element::from_symbol("C").expect("carbon"))
    }

    #[test]
    fn public_macro_builder_is_connected_by_construction() {
        let mut builder = MacroMoleculeBuilder::new();
        let first = builder.add_atom(carbon()).expect("first atom");
        assert_eq!(
            builder.add_atom(carbon()),
            Err(MacroGraphEditError::DisconnectedAtomInsertion)
        );
        let (second, _) = builder
            .add_atom_bonded_to(first, carbon(), BondOrder::Single)
            .expect("bonded atom");
        assert!(builder.graph().bond_between(first, second).unwrap().is_some());
        assert!(builder.graph().is_connected());
    }
}
