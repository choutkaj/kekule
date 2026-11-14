//! Publication of a draft as a [`Molecule`], and the canonical oxohalogen
//! representation it requires.
use crate::draft::{Draft, DraftError};
use crate::fixed_u32_id::FixedId;

use super::super::*;
use super::stereo::stereo_violations;
use super::MoleculeEditor;

/// An internal graph-storage invariant violation. The public API cannot
/// produce one, so publication checks it only in debug builds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum GraphValidationError {
    AdjacencySlotCount,
    AtomPropertySlotCount { expected: usize, actual: usize },
    BondPropertySlotCount { expected: usize, actual: usize },
    TombstonedAtomHasAdjacency { atom: AtomId },
    InvalidBondEndpoint { bond: BondId },
    InvalidAdjacencyBond { atom: AtomId, bond: BondId },
    AdjacencyEndpointMismatch { atom: AtomId, bond: BondId },
    DuplicateAdjacencyEntry { atom: AtomId, bond: BondId },
    MissingAdjacencyEntry { atom: AtomId, bond: BondId },
}

impl MoleculeEditor {
    /// Publishes represented chemistry with dense IDs.
    ///
    /// Live atoms, bonds, stereo elements, and stereo groups are renumbered
    /// `0..n` in their draft order, so deleted draft slots never survive
    /// publication. Use [`Self::finish_with_correspondence`] to translate draft
    /// IDs. Unchanged chemistry retains its installed perception; graph edits
    /// invalidate affected perception when they occur. Install reconstructed
    /// perception on the resulting [`Molecule::install_perception`].
    ///
    /// ```compile_fail
    /// use kekule::{molecule::MoleculeEditor, perception::Perception};
    /// let mut editor = MoleculeEditor::new();
    /// editor.install_perception(Perception::default()).unwrap();
    /// ```
    pub fn finish(self) -> Result<Molecule, MoleculePublicationError> {
        self.publish().map(|(molecule, _)| molecule)
    }

    /// Publishes like [`Self::finish`] and maps every surviving draft ID to its
    /// published ID. Without deletions the correspondence is the identity.
    ///
    /// ```
    /// use kekule::molecule::{Atom, AtomId, BondOrder, Element, MoleculeEditor};
    ///
    /// let carbon = || Atom::new(Element::from_symbol("C").unwrap());
    /// let mut editor = MoleculeEditor::new();
    /// let first = editor.add_atom(carbon())?;
    /// let removed = editor.add_atom(carbon())?;
    /// let last = editor.add_atom(carbon())?;
    /// editor.add_bond(first, last, BondOrder::Single)?;
    /// editor.delete_atom(removed)?;
    /// let (molecule, ids) = editor.finish_with_correspondence()?;
    ///
    /// assert_eq!(ids.atom(last), Some(AtomId::new(1)));
    /// assert_eq!(ids.atom(removed), None);
    /// assert_eq!(molecule.atom_ids().collect::<Vec<_>>(), [AtomId::new(0), AtomId::new(1)]);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn finish_with_correspondence(
        self,
    ) -> Result<(Molecule, MoleculeCorrespondence), MoleculePublicationError> {
        self.publish().map(|(molecule, compaction)| {
            (
                molecule,
                MoleculeCorrespondence::from_compaction(&compaction),
            )
        })
    }

    /// The error [`Self::finish`] would return, leaving the draft unchanged.
    pub fn validate(&self) -> Result<(), MoleculePublicationError> {
        self.check()
    }

    /// Publishes, or returns the original editor for repair through
    /// [`DraftError::into_draft`].
    pub fn try_finish(self) -> Result<Molecule, DraftError<Self, MoleculePublicationError>> {
        self.attempt()
    }
}

// Validation is cheaper than publication, so recovery checks the draft
// instead of publishing a copy.
impl Draft for MoleculeEditor {
    type Published = Molecule;
    type Error = MoleculePublicationError;

    fn into_published(self) -> Result<Molecule, MoleculePublicationError> {
        self.finish()
    }

    fn check(&self) -> Result<(), MoleculePublicationError> {
        if self.graph.atom_count() == 0 {
            return Err(MoleculePublicationError::EmptyGraph);
        }
        debug_assert_eq!(
            validate_graph(self),
            Ok(()),
            "editors keep graph storage consistent"
        );
        let components = self.graph.connected_components().len();
        if components != 1 {
            return Err(MoleculePublicationError::Disconnected { components });
        }
        if let Some(issue) = stereo_violations(&self.graph).into_iter().next() {
            return Err(MoleculePublicationError::InvalidStereo(issue));
        }
        if let Some(atom) = unnormalized_oxohalogens(&self.graph).first() {
            return Err(MoleculePublicationError::UnnormalizedOxohalogen { atom: *atom });
        }
        Ok(())
    }

    fn attempt(self) -> Result<Molecule, DraftError<Self, MoleculePublicationError>> {
        match self.check() {
            Ok(()) => Ok(self.publish_validated().0),
            Err(error) => Err(DraftError::new(error, self)),
        }
    }
}

impl MoleculeEditor {
    /// Applies Kekule's canonical oxohalogen representation (the RDKit halogen
    /// cleanup): each neutral Cl, Br, or I with only oxygen neighbors and
    /// explicit valence 3, 5, or 7 becomes `[X+n]` with its `n` oxo double
    /// bonds rewritten as single bonds to `[O-]`. Returns the rewritten
    /// halogens; a rewrite clears perception and owner annotations, like any
    /// chemistry edit. Publication rejects unnormalized oxohalogens instead of
    /// rewriting them silently.
    ///
    /// ```
    /// use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor, MoleculePublicationError};
    /// let mut editor = MoleculeEditor::new();
    /// let chlorine = editor.add_atom(Atom::new(Element::from_symbol("Cl").unwrap()))?;
    /// for order in [BondOrder::Single, BondOrder::Double, BondOrder::Double, BondOrder::Double] {
    ///     let oxygen = editor.add_atom(Atom::new(Element::from_symbol("O").unwrap()))?;
    ///     editor.add_bond(chlorine, oxygen, order)?;
    /// }
    /// assert!(matches!(
    ///     editor.validate(),
    ///     Err(MoleculePublicationError::UnnormalizedOxohalogen { .. })
    /// ));
    /// assert_eq!(editor.normalize_oxohalogens(), [chlorine]);
    /// let perchlorate = editor.finish()?;
    /// assert_eq!(perchlorate.atom(chlorine).map(|atom| atom.formal_charge()), Some(3));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn normalize_oxohalogens(&mut self) -> Vec<AtomId> {
        let halogens = unnormalized_oxohalogens(&self.graph);
        for &atom_id in &halogens {
            let oxo_bonds = oxo_bonds_to_oxygen(&self.graph, atom_id);
            // Explicit valence is at most seven, hence at most three oxo bonds.
            // Out-of-model graphs are left for valence validation.
            let formal_charge = i8::try_from(oxo_bonds.len()).expect("at most three oxo bonds");
            self.graph
                .atom_mut(atom_id)
                .expect("live halogen")
                .formal_charge = formal_charge;
            for (oxygen_id, bond_id) in oxo_bonds {
                self.graph
                    .atom_mut(oxygen_id)
                    .expect("live oxygen")
                    .formal_charge = -1;
                self.graph.bond_mut(bond_id).expect("live oxo bond").order = BondOrder::Single;
                self.prune_stereo_for_bond(bond_id);
            }
        }
        if !halogens.is_empty() {
            self.chemistry_changed();
        }
        halogens
    }

    fn publish(self) -> Result<(Molecule, SlotCompaction), MoleculePublicationError> {
        self.validate()?;
        Ok(self.publish_validated())
    }

    /// Canonicalizes stored stereo, which drops perception computed before
    /// it changed, and renumbers densely.
    fn publish_validated(mut self) -> (Molecule, SlotCompaction) {
        if self.canonicalize_stored_stereo_elements() {
            self.clear_perception();
        }
        self.into_dense()
    }

    /// A dense copy of the draft without publication's checks or
    /// canonicalization, for internal algorithms that need a [`Molecule`]
    /// reading a draft. It may be empty or disconnected and must not escape
    /// the crate. IDs are renumbered as at publication; a draft without
    /// deleted slots keeps them.
    pub(crate) fn snapshot(&self) -> (Molecule, SlotCompaction) {
        self.clone().into_dense()
    }

    /// Publishes without checks or canonicalization, for tests that perceive
    /// or inject what publication would reject.
    #[cfg(test)]
    pub(crate) fn into_unchecked(self) -> Molecule {
        self.into_dense().0
    }

    /// [`Self::snapshot`] of a draft that has deleted no atoms or bonds, so
    /// the copy keeps their IDs.
    pub(crate) fn snapshot_keeping_ids(&self) -> Molecule {
        debug_assert!(!self.graph.atoms.has_deleted() && !self.graph.bonds.has_deleted());
        self.snapshot().0
    }

    /// Renumbers the draft densely. Renumbered atoms or bonds clear
    /// perception, which their deletion already invalidated; renumbered stereo
    /// elements clear only CIP state.
    fn into_dense(self) -> (Molecule, SlotCompaction) {
        let atoms_changed = self.graph.atoms.has_deleted();
        let bonds_changed = self.graph.bonds.has_deleted();
        let elements_changed = self.graph.stereo_elements.has_deleted();
        let mut perception = self.perception;
        let mut properties = self.properties;
        let (graph, compaction) = self.graph.into_dense();
        if atoms_changed {
            *properties.atoms_mut() = properties
                .atoms()
                .raw()
                .select_indices(&compaction.live_atom_slots())
                .expect("draft atom properties match atom slots");
        }
        if bonds_changed {
            *properties.bonds_mut() = properties
                .bonds()
                .raw()
                .select_indices(&compaction.live_bond_slots())
                .expect("draft bond properties match bond slots");
        }
        if atoms_changed || bonds_changed {
            perception = Perception::default();
        } else if elements_changed {
            perception.stereo = None;
        }
        let molecule = Molecule {
            graph,
            perception,
            properties,
        };
        (molecule, compaction)
    }
}

pub(super) fn validate_graph(editor: &MoleculeEditor) -> Result<(), GraphValidationError> {
    let graph = &editor.graph;
    if editor.properties.atoms().len() != graph.atom_slot_count() {
        return Err(GraphValidationError::AtomPropertySlotCount {
            expected: graph.atom_slot_count(),
            actual: editor.properties.atoms().len(),
        });
    }
    if editor.properties.bonds().len() != graph.bond_slot_count() {
        return Err(GraphValidationError::BondPropertySlotCount {
            expected: graph.bond_slot_count(),
            actual: editor.properties.bonds().len(),
        });
    }
    if graph.adjacency.len() != graph.atom_slot_count() {
        return Err(GraphValidationError::AdjacencySlotCount);
    }
    for (index, adjacency) in graph.adjacency.iter().enumerate() {
        let atom_id = AtomId::from_index(index);
        if graph.atom(atom_id).is_none() && !adjacency.is_empty() {
            return Err(GraphValidationError::TombstonedAtomHasAdjacency { atom: atom_id });
        }
        let mut seen = std::collections::BTreeSet::new();
        for &bond_id in adjacency {
            if !seen.insert(bond_id) {
                return Err(GraphValidationError::DuplicateAdjacencyEntry {
                    atom: atom_id,
                    bond: bond_id,
                });
            }
            let Some(bond) = graph.bond(bond_id) else {
                return Err(GraphValidationError::InvalidAdjacencyBond {
                    atom: atom_id,
                    bond: bond_id,
                });
            };
            if ![bond.a(), bond.b()].contains(&atom_id) {
                return Err(GraphValidationError::AdjacencyEndpointMismatch {
                    atom: atom_id,
                    bond: bond_id,
                });
            }
        }
    }
    for (bond_id, bond) in graph.bonds() {
        for atom in [bond.a(), bond.b()] {
            if graph.atom(atom).is_none() {
                return Err(GraphValidationError::InvalidBondEndpoint { bond: bond_id });
            }
            if !graph.adjacency[atom.index()].contains(&bond_id) {
                return Err(GraphValidationError::MissingAdjacencyEntry {
                    atom,
                    bond: bond_id,
                });
            }
        }
    }
    Ok(())
}

/// Rewrites every in-model neutral oxohalogen and canonicalizes stored stereo.
/// Format interpreters apply this to imported chemistry before publication.
pub(crate) fn canonicalize_represented_chemistry(editor: &mut MoleculeEditor) {
    editor.normalize_oxohalogens();
    editor.canonicalize_stored_stereo_elements();
}

/// Neutral in-model oxohalogens that still carry oxo double bonds.
fn unnormalized_oxohalogens<L: SlotLiveness>(graph: &Graph<L>) -> Vec<AtomId> {
    graph
        .atoms()
        .filter_map(|(atom_id, atom)| {
            (atom.formal_charge == 0
                && matches!(atom.element.symbol(), "Cl" | "Br" | "I")
                && is_rdkit_oxohalogen(graph, atom_id, atom)
                && !oxo_bonds_to_oxygen(graph, atom_id).is_empty())
            .then_some(atom_id)
        })
        .collect()
}

fn is_rdkit_oxohalogen<L: SlotLiveness>(graph: &Graph<L>, atom_id: AtomId, atom: &Atom) -> bool {
    // RDKit 2026.03.3 halogenCleanup: oxygen-only neighbors and represented
    // explicit valence 3, 5 or 7. Bridging ester oxygen is allowed. This reads
    // represented bond/H contributions, never installed perception.
    let explicit = crate::algorithms::explicit_valence(graph, atom_id)
        + usize::from(atom.hydrogens.represented_count());
    if !matches!(explicit, 3 | 5 | 7) {
        return false;
    }
    graph
        .incident_bonds(atom_id)
        .into_iter()
        .flatten()
        .all(|(_, bond)| {
            graph
                .atom(bond.other_atom(atom_id))
                .is_some_and(|neighbor| neighbor.element.symbol() == "O")
        })
}

fn oxo_bonds_to_oxygen<L: SlotLiveness>(
    graph: &Graph<L>,
    atom_id: AtomId,
) -> Vec<(AtomId, BondId)> {
    graph
        .incident_bonds(atom_id)
        .into_iter()
        .flatten()
        .filter_map(|(bond_id, bond)| {
            if bond.order != BondOrder::Double {
                return None;
            }
            let oxygen_id = bond.other_atom(atom_id);
            let oxygen = graph.atom(oxygen_id)?;
            (oxygen.element.symbol() == "O").then_some((oxygen_id, bond_id))
        })
        .collect()
}
