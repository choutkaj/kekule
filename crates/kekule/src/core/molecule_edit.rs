use std::fmt;

use super::{Atom, AtomId, BondId, BondOrder, Molecule, Result};

/// A connectedness violation at a public [`Molecule`] boundary.
///
/// Empty and single-atom molecules are accepted. Every nontrivial molecule must
/// contain exactly one graph component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MoleculeConnectivityError {
    components: usize,
}

impl MoleculeConnectivityError {
    pub const fn components(self) -> usize {
        self.components
    }
}

impl fmt::Display for MoleculeConnectivityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "molecule must be connected, but contains {} graph components",
            self.components
        )
    }
}

impl std::error::Error for MoleculeConnectivityError {}

impl Molecule {
    /// Starts construction of one molecule.
    ///
    /// The builder may be temporarily disconnected. [`MoleculeBuilder::build`]
    /// validates the final connectedness invariant before returning a public
    /// `Molecule`.
    pub fn builder() -> MoleculeBuilder {
        MoleculeBuilder::new()
    }

    /// Starts a transactional topology edit.
    ///
    /// The working copy may be temporarily disconnected. [`MoleculeEditor::commit`]
    /// publishes it only when the final graph is connected, otherwise the
    /// original molecule is left unchanged.
    pub fn edit(&mut self) -> MoleculeEditor<'_> {
        MoleculeEditor {
            working: self.clone(),
            target: self,
        }
    }

    /// Returns whether this molecule satisfies the public connectedness invariant.
    ///
    /// Empty and single-atom molecules are treated as valid connected boundary
    /// cases so `Default`/`new` remain useful lightweight values.
    pub fn is_connected(&self) -> bool {
        self.atom_count() <= 1 || self.connected_components().len() == 1
    }

    /// Validates the public connectedness invariant.
    pub fn validate_connected(&self) -> std::result::Result<(), MoleculeConnectivityError> {
        if self.is_connected() {
            return Ok(());
        }
        Err(MoleculeConnectivityError {
            components: self.connected_components().len(),
        })
    }
}

/// Mutable construction state for one final connected [`Molecule`].
///
/// Connectivity is intentionally checked only by [`Self::build`], allowing
/// ordinary graph construction such as adding all atoms before adding bonds.
/// The staging graph is not publicly borrowable, so even a cloned or replaced
/// builder must pass through `build` before yielding a `Molecule`.
///
/// ```compile_fail
/// use kekule::core::{Atom, Element, Molecule};
///
/// let mut builder = Molecule::builder();
/// builder.add_atom(Atom::new(Element::from_atomic_number(6).unwrap())).unwrap();
/// builder.add_atom(Atom::new(Element::from_atomic_number(8).unwrap())).unwrap();
/// let leaked: Molecule = builder.molecule().clone();
/// # let _ = leaked;
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MoleculeBuilder {
    molecule: Molecule,
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_atom(&mut self, atom: Atom) -> Result<AtomId> {
        self.molecule.add_atom(atom)
    }

    pub fn delete_atom(&mut self, atom: AtomId) -> Result<Atom> {
        self.molecule.delete_atom(atom)
    }

    pub fn add_bond(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<BondId> {
        self.molecule.add_bond(a, b, order)
    }

    pub fn delete_bond(&mut self, bond: BondId) -> Result<super::Bond> {
        self.molecule.delete_bond(bond)
    }

    pub fn build(self) -> std::result::Result<Molecule, MoleculeConnectivityError> {
        self.molecule.validate_connected()?;
        Ok(self.molecule)
    }
}

/// Transactional topology edit over one existing connected [`Molecule`].
///
/// The private working graph can be temporarily disconnected, but it cannot be
/// borrowed or extracted as a completed `Molecule`; only [`Self::commit`] can
/// publish it.
///
/// ```compile_fail
/// use kekule::core::{Atom, BondOrder, Element, Molecule};
///
/// let mut builder = Molecule::builder();
/// let carbon = builder.add_atom(Atom::new(Element::from_atomic_number(6).unwrap())).unwrap();
/// let oxygen = builder.add_atom(Atom::new(Element::from_atomic_number(8).unwrap())).unwrap();
/// let bond = builder.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
/// let mut molecule = builder.build().unwrap();
/// let mut editor = molecule.edit();
/// editor.delete_bond(bond).unwrap();
/// let leaked: Molecule = std::mem::take(editor.molecule_mut());
/// # let _ = leaked;
/// ```
pub struct MoleculeEditor<'a> {
    working: Molecule,
    target: &'a mut Molecule,
}

impl MoleculeEditor<'_> {
    pub fn add_atom(&mut self, atom: Atom) -> Result<AtomId> {
        self.working.add_atom(atom)
    }

    pub fn delete_atom(&mut self, atom: AtomId) -> Result<Atom> {
        self.working.delete_atom(atom)
    }

    pub fn add_bond(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<BondId> {
        self.working.add_bond(a, b, order)
    }

    pub fn delete_bond(&mut self, bond: BondId) -> Result<super::Bond> {
        self.working.delete_bond(bond)
    }

    pub fn commit(self) -> std::result::Result<(), MoleculeConnectivityError> {
        self.working.validate_connected()?;
        *self.target = self.working;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Element;

    fn carbon() -> Atom {
        Atom::new(Element::from_atomic_number(6).expect("carbon"))
    }

    #[test]
    fn builder_rejects_disconnected_final_graph() {
        let mut builder = Molecule::builder();
        builder.add_atom(carbon()).expect("first atom");
        builder.add_atom(carbon()).expect("second atom");
        assert_eq!(
            builder.build(),
            Err(MoleculeConnectivityError { components: 2 })
        );
    }

    #[test]
    fn builder_accepts_connected_graph() {
        let mut builder = Molecule::builder();
        let a = builder.add_atom(carbon()).expect("first atom");
        let b = builder.add_atom(carbon()).expect("second atom");
        builder
            .add_bond(a, b, BondOrder::Single)
            .expect("connecting bond");
        assert!(builder.build().expect("connected molecule").is_connected());
    }

    #[test]
    fn failed_editor_commit_leaves_original_unchanged() {
        let mut builder = Molecule::builder();
        let a = builder.add_atom(carbon()).expect("first atom");
        let b = builder.add_atom(carbon()).expect("second atom");
        let bond = builder
            .add_bond(a, b, BondOrder::Single)
            .expect("connecting bond");
        let mut molecule = builder.build().expect("connected molecule");
        let before = molecule.clone();

        let mut editor = molecule.edit();
        editor.delete_bond(bond).expect("delete in working copy");
        assert_eq!(
            editor.commit(),
            Err(MoleculeConnectivityError { components: 2 })
        );
        assert_eq!(molecule, before);
    }

    #[test]
    fn editor_can_pass_through_temporarily_disconnected_state() {
        let mut builder = Molecule::builder();
        let a = builder.add_atom(carbon()).expect("first atom");
        let b = builder.add_atom(carbon()).expect("second atom");
        let c = builder.add_atom(carbon()).expect("third atom");
        let ab = builder
            .add_bond(a, b, BondOrder::Single)
            .expect("first bond");
        builder
            .add_bond(b, c, BondOrder::Single)
            .expect("second bond");
        let mut molecule = builder.build().expect("connected molecule");

        let mut editor = molecule.edit();
        editor.delete_bond(ab).expect("temporary disconnect");
        editor
            .add_bond(a, c, BondOrder::Single)
            .expect("reconnect through another edge");
        editor.commit().expect("connected final graph");
        assert!(molecule.is_connected());
    }
}
