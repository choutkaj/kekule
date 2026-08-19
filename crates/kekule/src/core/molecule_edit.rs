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

/// Failure to publish checked molecule construction or editing state.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum MoleculePublicationError {
    DisconnectedGraph(MoleculeConnectivityError),
    FormalChargeOutOfRange { atom: AtomId, charge: usize },
}

impl fmt::Display for MoleculePublicationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DisconnectedGraph(error) => write!(formatter, "{error}"),
            Self::FormalChargeOutOfRange { atom, charge } => write!(
                formatter,
                "publishing atom {atom} requires formal charge +{charge}, which is outside the supported range"
            ),
        }
    }
}

impl std::error::Error for MoleculePublicationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::DisconnectedGraph(error) => Some(error),
            Self::FormalChargeOutOfRange { .. } => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CanonicalRepresentationError {
    pub(crate) atom: AtomId,
    pub(crate) charge: usize,
}

impl Molecule {
    /// Starts construction of one molecule.
    ///
    /// The builder may be temporarily disconnected. [`MoleculeBuilder::build`]
    /// validates connectedness and canonicalizes represented chemistry before
    /// returning a public `Molecule`.
    pub fn builder() -> MoleculeBuilder {
        MoleculeBuilder::new()
    }

    /// Starts a transactional topology edit.
    ///
    /// The working copy may be temporarily disconnected. [`MoleculeEditor::commit`]
    /// publishes it only when the final graph is connected and representable in
    /// canonical form; otherwise the original molecule is left unchanged.
    /// Canonicality-sensitive atom, bond, and topology mutation is available on
    /// the editor rather than directly on a published molecule.
    ///
    /// ```compile_fail
    /// use kekule::core::{Atom, BondOrder, Element, Molecule};
    ///
    /// let mut builder = Molecule::builder();
    /// let chlorine = builder.add_atom(Atom::new(Element::from_symbol("Cl").unwrap())).unwrap();
    /// let oxygen = builder.add_atom(Atom::new(Element::from_symbol("O").unwrap())).unwrap();
    /// let bond = builder.add_bond(chlorine, oxygen, BondOrder::Single).unwrap();
    /// let mut molecule = builder.build().unwrap();
    /// molecule.bond_mut(bond).unwrap().order = BondOrder::Double;
    /// ```
    ///
    /// ```compile_fail
    /// use kekule::core::{Atom, Element, Molecule};
    ///
    /// let mut builder = Molecule::builder();
    /// let carbon = builder.add_atom(Atom::new(Element::from_symbol("C").unwrap())).unwrap();
    /// let mut molecule = builder.build().unwrap();
    /// molecule.atom_mut(carbon).unwrap().formal_charge = 1;
    /// ```
    ///
    /// ```compile_fail
    /// use kekule::core::{Atom, BondOrder, Element, Molecule};
    ///
    /// let mut builder = Molecule::builder();
    /// let first = builder.add_atom(Atom::new(Element::from_symbol("C").unwrap())).unwrap();
    /// let middle = builder.add_atom(Atom::new(Element::from_symbol("C").unwrap())).unwrap();
    /// let last = builder.add_atom(Atom::new(Element::from_symbol("C").unwrap())).unwrap();
    /// builder.add_bond(first, middle, BondOrder::Single).unwrap();
    /// builder.add_bond(middle, last, BondOrder::Single).unwrap();
    /// let mut molecule = builder.build().unwrap();
    /// molecule.add_bond(first, last, BondOrder::Single).unwrap();
    /// ```
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
/// Connectivity and canonical representation are checked only by
/// [`Self::build`], allowing ordinary graph construction such as adding all
/// atoms before adding bonds.
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

    pub fn build(self) -> std::result::Result<Molecule, MoleculePublicationError> {
        let mut molecule = self.molecule;
        molecule
            .validate_connected()
            .map_err(MoleculePublicationError::DisconnectedGraph)?;
        canonicalize_represented_chemistry(&mut molecule).map_err(|error| {
            MoleculePublicationError::FormalChargeOutOfRange {
                atom: error.atom,
                charge: error.charge,
            }
        })?;
        Ok(molecule)
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
    /// Returns mutable represented atom state in this private working copy.
    pub fn atom_mut(&mut self, atom: AtomId) -> Result<super::AtomMut<'_>> {
        self.working.atom_mut(atom)
    }

    /// Returns mutable represented bond state in this private working copy.
    pub fn bond_mut(&mut self, bond: BondId) -> Result<super::BondMut<'_>> {
        self.working.bond_mut(bond)
    }

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

    pub fn commit(self) -> std::result::Result<(), MoleculePublicationError> {
        let mut working = self.working;
        working
            .validate_connected()
            .map_err(MoleculePublicationError::DisconnectedGraph)?;
        canonicalize_represented_chemistry(&mut working).map_err(|error| {
            MoleculePublicationError::FormalChargeOutOfRange {
                atom: error.atom,
                charge: error.charge,
            }
        })?;
        *self.target = working;
        Ok(())
    }
}

pub(crate) fn canonicalize_represented_chemistry(
    molecule: &mut Molecule,
) -> std::result::Result<(), CanonicalRepresentationError> {
    let halogens = molecule
        .atoms()
        .filter_map(|(atom_id, atom)| {
            (atom.formal_charge == 0
                && matches!(atom.element.symbol(), "Cl" | "Br" | "I")
                && has_terminal_single_bond_oxygen_neighbor(molecule, atom_id))
            .then_some(atom_id)
        })
        .collect::<Vec<_>>();

    let mut rewritten = false;
    for atom_id in halogens {
        let oxo_bonds = oxo_bonds_to_neutral_oxygen(molecule, atom_id);
        if oxo_bonds.is_empty() {
            continue;
        }
        rewritten = true;
        let charge = oxo_bonds.len();
        let formal_charge = i8::try_from(charge).map_err(|_| CanonicalRepresentationError {
            atom: atom_id,
            charge,
        })?;

        if let Some(atom) = molecule.atoms[atom_id.index()].as_mut() {
            atom.formal_charge = formal_charge;
        }
        for (oxygen_id, bond_id) in oxo_bonds {
            if let Some(atom) = molecule.atoms[oxygen_id.index()].as_mut() {
                atom.formal_charge = -1;
            }
            if let Some(bond) = molecule.bonds[bond_id.index()].as_mut() {
                bond.order = BondOrder::Single;
            }
        }
    }
    if rewritten {
        molecule.invalidate_topology();
    }
    Ok(())
}

fn has_terminal_single_bond_oxygen_neighbor(molecule: &Molecule, atom_id: AtomId) -> bool {
    molecule
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(_, bond)| {
            let oxygen_id = bond.other_atom(atom_id);
            bond.order == BondOrder::Single
                && molecule
                    .atom(oxygen_id)
                    .is_ok_and(|neighbor| neighbor.element.symbol() == "O")
                && molecule.incident_bonds(oxygen_id).is_ok_and(|mut bonds| {
                    bonds.all(|(_, oxygen_bond)| {
                        let neighbor_id = oxygen_bond.other_atom(oxygen_id);
                        neighbor_id == atom_id
                            || molecule
                                .atom(neighbor_id)
                                .is_ok_and(|neighbor| neighbor.element.symbol() == "H")
                    })
                })
        })
}

fn oxo_bonds_to_neutral_oxygen(molecule: &Molecule, atom_id: AtomId) -> Vec<(AtomId, BondId)> {
    molecule
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter_map(|(bond_id, bond)| {
            if bond.order != BondOrder::Double {
                return None;
            }
            let oxygen_id = bond.other_atom(atom_id);
            let oxygen = molecule.atom(oxygen_id).ok()?;
            (oxygen.element.symbol() == "O" && oxygen.formal_charge == 0)
                .then_some((oxygen_id, bond_id))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Element;

    fn carbon() -> Atom {
        Atom::new(Element::from_atomic_number(6).expect("carbon"))
    }

    fn atom(symbol: &str) -> Atom {
        Atom::new(Element::from_symbol(symbol).expect("fixture element"))
    }

    #[test]
    fn builder_rejects_disconnected_final_graph() {
        let mut builder = Molecule::builder();
        builder.add_atom(carbon()).expect("first atom");
        builder.add_atom(carbon()).expect("second atom");
        assert_eq!(
            builder.build(),
            Err(MoleculePublicationError::DisconnectedGraph(
                MoleculeConnectivityError { components: 2 }
            ))
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
    fn builder_publishes_canonical_hypervalent_representation() {
        let mut builder = Molecule::builder();
        let chlorine = builder.add_atom(atom("Cl")).expect("chlorine");
        let oxo = builder.add_atom(atom("O")).expect("oxo oxygen");
        let hydroxyl = builder.add_atom(atom("O")).expect("hydroxyl oxygen");
        builder
            .add_bond(chlorine, oxo, BondOrder::Double)
            .expect("source-convention oxo bond");
        builder
            .add_bond(chlorine, hydroxyl, BondOrder::Single)
            .expect("hydroxyl bond");

        let molecule = builder.build().expect("canonical publication");

        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
        assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
        let oxo_bond = molecule.bond_between(chlorine, oxo).unwrap().unwrap();
        assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
        assert_eq!(
            molecule.perception(),
            &super::super::PerceptionState::default()
        );
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
            Err(MoleculePublicationError::DisconnectedGraph(
                MoleculeConnectivityError { components: 2 }
            ))
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

    #[test]
    fn editor_canonicalizes_before_transactional_publication() {
        let mut builder = Molecule::builder();
        let chlorine = builder.add_atom(atom("Cl")).expect("chlorine");
        let hydroxyl = builder.add_atom(atom("O")).expect("hydroxyl oxygen");
        builder
            .add_bond(chlorine, hydroxyl, BondOrder::Single)
            .expect("initial bond");
        let mut molecule = builder.build().expect("initial molecule");

        let mut editor = molecule.edit();
        let oxo = editor.add_atom(atom("O")).expect("oxo oxygen");
        editor
            .add_bond(chlorine, oxo, BondOrder::Double)
            .expect("source-convention oxo bond");
        editor.commit().expect("canonical edit publication");

        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
        assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
        let oxo_bond = molecule.bond_between(chlorine, oxo).unwrap().unwrap();
        assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
    }

    #[test]
    fn editor_canonicalizes_direct_bond_order_changes_before_publication() {
        let mut builder = Molecule::builder();
        let chlorine = builder.add_atom(atom("Cl")).unwrap();
        let anchor = builder.add_atom(atom("O")).unwrap();
        let oxo = builder.add_atom(atom("O")).unwrap();
        builder
            .add_bond(chlorine, anchor, BondOrder::Single)
            .unwrap();
        let oxo_bond = builder.add_bond(chlorine, oxo, BondOrder::Single).unwrap();
        let mut molecule = builder.build().unwrap();

        let mut editor = molecule.edit();
        editor.bond_mut(oxo_bond).unwrap().order = BondOrder::Double;
        editor.commit().expect("canonical edit publication");

        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
        assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
        assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
    }

    #[test]
    fn failed_editor_canonicalization_leaves_the_published_molecule_unchanged() {
        let mut builder = Molecule::builder();
        let chlorine = builder.add_atom(atom("Cl")).unwrap();
        let anchor = builder.add_atom(atom("O")).unwrap();
        builder
            .add_bond(chlorine, anchor, BondOrder::Single)
            .unwrap();
        let mut oxo_bonds = Vec::new();
        for _ in 0..128 {
            let oxygen = builder.add_atom(atom("O")).unwrap();
            oxo_bonds.push(
                builder
                    .add_bond(chlorine, oxygen, BondOrder::Single)
                    .unwrap(),
            );
        }
        let mut molecule = builder.build().unwrap();
        let before = molecule.clone();

        let mut editor = molecule.edit();
        for bond in oxo_bonds {
            editor.bond_mut(bond).unwrap().order = BondOrder::Double;
        }
        assert_eq!(
            editor.commit(),
            Err(MoleculePublicationError::FormalChargeOutOfRange {
                atom: chlorine,
                charge: 128,
            })
        );
        assert_eq!(molecule, before);
    }
}
