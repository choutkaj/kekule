//! The published [`Molecule`] and its non-structural state: perception and
//! properties.
use std::collections::{BTreeMap, BTreeSet};
use std::ops::Deref;

use crate::algorithms::{
    add_hydrogens_to_molecule, perceive_molecule, perceive_molecule_with,
    remove_hydrogens_from_molecule, AddHydrogensOptions, AddHydrogensReport,
    HydrogenTransformError, PerceptionError, RemoveHydrogensReport,
};
use crate::parallel::Run;
use crate::properties::{LiveRows, MoleculeProperties, MoleculePropertiesMut};

use super::*;

/// One published, non-empty, connected, geometry-independent molecular entity.
///
/// A molecule owns three kinds of state with different semantics:
///
/// - [`Graph`] is authoritative represented chemistry: atoms, bonds,
///   connectivity, and represented stereochemistry.
/// - [`Perception`] is reconstructible derived chemistry such as valence, rings,
///   aromaticity, and installed CIP descriptors.
/// - [`MoleculeProperties`] contains extensible geometry-independent annotations scoped
///   to this molecular definition.
///
/// Coordinates, periodic cells, velocities, and biological hierarchy do not
/// belong to a molecule. Put explicit molecule instances and hierarchy in a
/// [`crate::topology::Topology`], then attach geometry with
/// [`crate::structure::Structure`].
///
/// Construction and structural editing publish through
/// [`MoleculeEditor::finish`]. Temporary disconnectedness is allowed inside the
/// editor but never in a published `Molecule`. Published atom, bond, and stereo
/// IDs are dense: `0..atom_count()` and likewise for the other ID spaces.
/// Equality compares authoritative graph chemistry; cached perception and
/// generic properties do not change molecular identity.
#[derive(Debug, Clone)]
pub struct Molecule {
    pub(crate) graph: Graph,
    pub(crate) perception: Perception,
    pub(crate) properties: MoleculeProperties,
}

impl PartialEq for Molecule {
    fn eq(&self, other: &Self) -> bool {
        self.graph == other.graph
    }
}

impl Deref for Molecule {
    type Target = Graph;

    /// Graph reads (`atom`, `atoms`, `neighbors`, `bond_between`, stereo, ...)
    /// are available on a molecule directly.
    fn deref(&self) -> &Graph {
        &self.graph
    }
}

impl Molecule {
    /// The authoritative represented chemistry, also reachable through `Deref`.
    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    /// Hashes the graph in local atom order: atoms (element, isotope, formal
    /// charge, radical, implicit hydrogens, atom map), bonds with their
    /// endpoints and orders, and the number of stereo elements. Equal graphs
    /// hash equally, so this keys a search for identical molecules.
    pub(crate) fn hash_graph(&self, state: &mut impl std::hash::Hasher) {
        use std::hash::Hash;
        let graph = &self.graph;
        graph.atom_count().hash(state);
        for (_, atom) in graph.atoms() {
            atom.element.hash(state);
            atom.isotope.hash(state);
            atom.formal_charge.hash(state);
            atom.radical.hash(state);
            atom.hydrogens.hash(state);
            atom.atom_map.hash(state);
        }
        graph.bond_count().hash(state);
        for (_, bond) in graph.bonds() {
            bond.a.hash(state);
            bond.b.hash(state);
            bond.order.hash(state);
        }
        graph.stereo_elements().len().hash(state);
    }

    /// Definition-scoped annotations: owner values plus one row per atom and
    /// bond, indexed by [`AtomId::index`] and [`BondId::index`].
    pub fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }

    /// Length-preserving mutable access to the annotations. Property edits
    /// never invalidate perception or change molecular equality.
    pub fn properties_mut(&mut self) -> MoleculePropertiesMut<'_> {
        MoleculePropertiesMut::new(&mut self.properties, LiveRows::All, LiveRows::All)
    }

    /// The installed derived chemistry; sections not yet perceived are absent.
    pub fn perception(&self) -> &Perception {
        &self.perception
    }

    /// Validates and atomically replaces the complete installed perception state.
    ///
    /// The detached state is checked against the graph and its stereo before
    /// mutation. Failure leaves the previous state exactly unchanged.
    pub fn install_perception(&mut self, state: Perception) -> Result<(), PerceptionInstallError> {
        validate_perception(self, &state)?;
        self.perception = state;
        Ok(())
    }

    /// Replaces one atom's incident-bond traversal order, which publication
    /// does not normalize, for exact reconstruction. Returns `false` and
    /// changes nothing unless `order` is a permutation of the atom's current
    /// incident bonds.
    pub(crate) fn restore_adjacency_order(&mut self, atom: AtomId, order: Vec<BondId>) -> bool {
        let Some(current) = self.graph.adjacency.get(atom.index()) else {
            return false;
        };
        if current.len() != order.len() {
            return false;
        }
        let mut sorted_current = current.clone();
        let mut sorted_order = order.clone();
        sorted_current.sort_unstable();
        sorted_order.sort_unstable();
        if sorted_current != sorted_order {
            return false;
        }
        self.graph.adjacency[atom.index()] = order;
        true
    }

    /// Counts explicit hydrogen atoms among this atom's graph neighbors.
    /// Includes isotopic and charged hydrogen neighbors; requires no perception.
    /// A hydrogen joined by a zero-order or dative bond, which valence does not
    /// count, is not counted either.
    /// Returns `None` when `atom` names no atom of this molecule.
    pub fn explicit_hydrogens(&self, atom: AtomId) -> Option<usize> {
        explicit_hydrogen_count(&self.graph, atom)
    }

    /// Counts implicit (non-graph) hydrogens, excluding graph H.
    /// Fixed counts are available without perception. For inferred counts,
    /// `None` means the count has not been perceived or was invalidated by a
    /// chemical edit; it also means `atom` names no atom of this molecule.
    /// This method never runs perception implicitly.
    pub fn implicit_hydrogens(&self, atom: AtomId) -> Option<usize> {
        implicit_hydrogen_count(&self.graph, &self.perception, atom)
    }

    /// Counts explicit hydrogen neighbors plus all implicit hydrogens.
    /// Returns `None` when the implicit count is unresolved or `atom` names no
    /// atom of this molecule.
    pub fn total_hydrogens(&self, atom: AtomId) -> Option<usize> {
        let explicit = self.explicit_hydrogens(atom)?;
        self.implicit_hydrogens(atom)
            .map(|implicit| explicit + implicit)
    }

    /// Removes all installed derived perception without changing represented chemistry.
    pub fn clear_perception(&mut self) {
        self.perception = Perception::default();
    }

    pub(crate) fn invalidate_cip(&mut self) {
        self.perception.cip = None;
    }

    fn invalidate_aromaticity(&mut self) {
        self.perception.aromaticity = None;
        self.perception.conjugation = None;
        self.perception.resonance = None;
        self.invalidate_cip();
    }

    pub(crate) fn install_valence(
        &mut self,
        model: ValenceModel,
        validation: ValenceValidation,
        inferred_hydrogens: BTreeMap<AtomId, u8>,
    ) {
        self.perception.valence = Some(ValencePerception {
            model,
            validation,
            inferred_hydrogens,
        });
        self.invalidate_aromaticity();
    }

    /// Overrides one inferred count of a private projection copy, keeping the
    /// installed provenance. A copy without valence takes the permissive
    /// RDKit-like provenance under which projected counts are written.
    pub(crate) fn set_inferred_hydrogens(&mut self, atom: AtomId, count: u8) {
        self.perception
            .valence
            .get_or_insert_with(|| ValencePerception {
                model: ValenceModel::RdkitLike,
                validation: ValenceValidation::Permissive,
                inferred_hydrogens: BTreeMap::new(),
            })
            .inferred_hydrogens
            .insert(atom, count);
        self.invalidate_aromaticity();
    }

    pub(crate) fn install_ring_membership(&mut self, membership: RingMembership) {
        self.perception.rings = Some(RingPerception {
            membership,
            basis: None,
        });
        self.invalidate_aromaticity();
    }

    pub(crate) fn install_ring_basis(&mut self, membership: RingMembership, basis: RingBasisState) {
        self.perception.rings = Some(RingPerception {
            membership,
            basis: Some(basis),
        });
        self.invalidate_aromaticity();
    }

    pub(crate) fn begin_aromaticity(&mut self, model: AromaticityModel) {
        self.perception.conjugation = None;
        self.perception.resonance = None;
        self.perception.aromaticity = Some(AromaticityPerception {
            model,
            atoms: BTreeSet::new(),
            bonds: BTreeSet::new(),
        });
        self.perception.cip = None;
    }

    pub(crate) fn set_atom_aromatic(&mut self, atom: AtomId, aromatic: bool) {
        let Some(state) = self.perception.aromaticity.as_mut() else {
            return;
        };
        if aromatic {
            state.atoms.insert(atom);
        } else {
            state.atoms.remove(&atom);
        }
    }

    pub(crate) fn set_bond_aromatic(&mut self, bond: BondId, aromatic: bool) {
        let Some(state) = self.perception.aromaticity.as_mut() else {
            return;
        };
        if aromatic {
            state.bonds.insert(bond);
        } else {
            state.bonds.remove(&bond);
        }
    }

    pub(crate) fn install_cip_descriptor(
        &mut self,
        element: StereoElementId,
        descriptor: StereoDescriptor,
    ) {
        self.perception
            .cip
            .get_or_insert_default()
            .cip_descriptors
            .insert(element, descriptor);
    }

    pub(crate) fn replace_cip_perception(
        &mut self,
        state: Option<CipPerception>,
    ) -> Option<CipPerception> {
        std::mem::replace(&mut self.perception.cip, state)
    }

    /// Install default valence, rings, aromaticity, conjugation and resonance
    /// delocalization transactionally.
    ///
    /// This derives perception from the canonical represented chemistry and
    /// never rewrites atoms, bonds, or represented stereochemistry.
    ///
    /// Interpretation publishes canonical represented chemistry directly;
    /// there is no separate public normalization step before perception.
    pub fn perceive(&mut self) -> Result<(), PerceptionError> {
        perceive_molecule(self)
    }

    /// [`Self::perceive`] under a [`Run`], on the calling thread. Its
    /// cancellation is also polled during ring and resonance perception,
    /// which dominate on large molecules. Progress counts the stages done,
    /// of five: valence, rings, aromaticity, conjugation, and resonance.
    ///
    /// # Errors
    ///
    /// As [`Self::perceive`], and [`PerceptionError::Cancelled`] once `run`'s
    /// cancellation is requested.
    pub fn perceive_with(&mut self, run: Run<'_>) -> Result<(), PerceptionError> {
        perceive_molecule_with(self, run)
    }

    /// Convert resolved implicit hydrogens into explicit graph atoms under the
    /// count and growth policy of `options`.
    /// Fixed counts require no perception; inference-enabled counts must be
    /// perceived first. Success invalidates perception but preserves composition.
    pub fn add_hydrogens(
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

/// [`Molecule::explicit_hydrogens`] for a published or draft graph.
pub(crate) fn explicit_hydrogen_count<L: SlotLiveness>(
    graph: &Graph<L>,
    atom: AtomId,
) -> Option<usize> {
    Some(
        graph
            .incident_bonds(atom)?
            .filter(|(_, bond)| {
                // Valence counts neither bond, so neither does the count.
                !matches!(bond.order, BondOrder::Zero | BondOrder::Dative)
                    && graph
                        .atom(bond.other_atom(atom))
                        .expect("bonds join live atoms")
                        .element
                        .atomic_number()
                        == 1
            })
            .count(),
    )
}

/// [`Molecule::implicit_hydrogens`] for a published or draft graph and its
/// perception.
pub(crate) fn implicit_hydrogen_count<L: SlotLiveness>(
    graph: &Graph<L>,
    perception: &Perception,
    atom: AtomId,
) -> Option<usize> {
    graph.atom(atom)?.hydrogens.resolve(
        perception
            .valence()
            .and_then(|valence| valence.inferred_hydrogen_count(atom)),
    )
}

#[cfg(test)]
mod tests;
