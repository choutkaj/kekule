use std::fmt;

use crate::properties::Properties;

use super::{
    Atom, AtomId, Bond, BondId, BondOrder, Graph, Molecule, Perception, Result, StereoElement,
    StereoElementId, StereoGroup, StereoGroupId,
};

/// A connectedness violation at a public [`Molecule`] boundary.
///
/// Every published molecule must contain exactly one graph component.
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
    EmptyGraph,
    DisconnectedGraph(MoleculeConnectivityError),
    InvalidGraph(GraphValidationError),
    InvalidStereo(StereoPublicationError),
    FormalChargeOutOfRange { atom: AtomId, charge: usize },
}

impl fmt::Display for MoleculePublicationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyGraph => formatter.write_str("molecule must contain at least one atom"),
            Self::DisconnectedGraph(error) => write!(formatter, "{error}"),
            Self::InvalidGraph(error) => write!(formatter, "invalid molecular graph: {error}"),
            Self::InvalidStereo(error) => write!(formatter, "invalid represented stereo: {error}"),
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
            Self::InvalidGraph(error) => Some(error),
            Self::InvalidStereo(error) => Some(error),
            Self::EmptyGraph | Self::FormalChargeOutOfRange { .. } => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum GraphValidationError {
    AdjacencySlotCount,
    TombstonedAtomHasAdjacency { atom: AtomId },
    InvalidBondEndpoint { bond: BondId },
    InvalidAdjacencyBond { atom: AtomId, bond: BondId },
    AdjacencyEndpointMismatch { atom: AtomId, bond: BondId },
    DuplicateAdjacencyEntry { atom: AtomId, bond: BondId },
    MissingAdjacencyEntry { atom: AtomId, bond: BondId },
}

impl fmt::Display for GraphValidationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AdjacencySlotCount => formatter.write_str(
                "molecular graph adjacency has a different slot count than atom storage",
            ),
            Self::TombstonedAtomHasAdjacency { atom } => {
                write!(
                    formatter,
                    "tombstoned atom atom{} still has adjacency entries",
                    atom.raw()
                )
            }
            Self::InvalidBondEndpoint { bond } => {
                write!(formatter, "bond bond{} references a missing atom", bond.raw())
            }
            Self::InvalidAdjacencyBond { atom, bond } => write!(
                formatter,
                "atom atom{} adjacency references missing bond bond{}",
                atom.raw(),
                bond.raw()
            ),
            Self::AdjacencyEndpointMismatch { atom, bond } => write!(
                formatter,
                "atom atom{} adjacency contains bond bond{}, but that bond is not incident to the atom",
                atom.raw(),
                bond.raw()
            ),
            Self::DuplicateAdjacencyEntry { atom, bond } => write!(
                formatter,
                "atom atom{} adjacency contains bond bond{} more than once",
                atom.raw(),
                bond.raw()
            ),
            Self::MissingAdjacencyEntry { atom, bond } => write!(
                formatter,
                "bond bond{} is missing from atom atom{} adjacency",
                bond.raw(),
                atom.raw()
            ),
        }
    }
}

impl std::error::Error for GraphValidationError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum StereoPublicationError {
    InvalidElementReference {
        element: StereoElementId,
    },
    InvalidElementGroup {
        element: StereoElementId,
    },
    EmptyGroup {
        group: StereoGroupId,
    },
    InvalidGroupMember {
        group: StereoGroupId,
        element: StereoElementId,
    },
    DuplicateGroupMember {
        group: StereoGroupId,
        element: StereoElementId,
    },
    InconsistentGroupMembership {
        group: StereoGroupId,
        element: StereoElementId,
    },
}

impl fmt::Display for StereoPublicationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidElementReference { element } => write!(
                formatter,
                "stereo element stereo{} references an invalid atom or bond",
                element.raw()
            ),
            Self::InvalidElementGroup { element } => write!(
                formatter,
                "stereo element stereo{} references a missing stereo group",
                element.raw()
            ),
            Self::EmptyGroup { group } => {
                write!(
                    formatter,
                    "stereo group group{} contains no members",
                    group.raw()
                )
            }
            Self::InvalidGroupMember { group, element } => write!(
                formatter,
                "stereo group group{} references missing stereo element stereo{}",
                group.raw(),
                element.raw()
            ),
            Self::DuplicateGroupMember { group, element } => write!(
                formatter,
                "stereo group group{} contains stereo element stereo{} more than once",
                group.raw(),
                element.raw()
            ),
            Self::InconsistentGroupMembership { group, element } => write!(
                formatter,
                "stereo group group{} and stereo element stereo{} disagree about membership",
                group.raw(),
                element.raw()
            ),
        }
    }
}

impl std::error::Error for StereoPublicationError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CanonicalRepresentationError {
    pub(crate) atom: AtomId,
    pub(crate) charge: usize,
}

impl Molecule {
    /// Starts a detached transaction from this published molecule.
    pub fn edit(&self) -> MoleculeEditor {
        MoleculeEditor::from_molecule(self)
    }

    #[cfg(test)]
    pub(crate) fn is_connected(&self) -> bool {
        self.validate_connected().is_ok()
    }

    /// Validates the public connectedness invariant.
    pub(crate) fn validate_connected(&self) -> std::result::Result<(), MoleculeConnectivityError> {
        let components = self.connected_components().len();
        if components == 1 {
            Ok(())
        } else {
            Err(MoleculeConnectivityError { components })
        }
    }
}
/// Transactional construction and structural editing state.
///
/// The working graph may be empty, disconnected, or chemically incomplete.
/// Only [`Self::finish`] can publish a `Molecule`.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeEditor {
    working: Molecule,
}

// Unit tests exercise chemistry algorithms against in-progress editor state.
#[cfg(test)]
impl std::ops::Deref for MoleculeEditor {
    type Target = Molecule;

    fn deref(&self) -> &Self::Target {
        &self.working
    }
}

#[cfg(test)]
impl std::ops::DerefMut for MoleculeEditor {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.working
    }
}

impl Default for MoleculeEditor {
    fn default() -> Self {
        Self {
            working: Molecule {
                graph: Graph::default(),
                perception: Perception::default(),
                properties: Properties::molecule(0, 0),
            },
        }
    }
}

impl MoleculeEditor {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_molecule(molecule: &Molecule) -> Self {
        Self {
            working: molecule.clone(),
        }
    }

    pub fn graph(&self) -> &Graph {
        &self.working.graph
    }

    pub(crate) fn working(&self) -> &Molecule {
        &self.working
    }

    pub(crate) fn working_mut(&mut self) -> &mut Molecule {
        &mut self.working
    }
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

    pub fn delete_bond(&mut self, bond: BondId) -> Result<Bond> {
        self.working.delete_bond(bond)
    }

    pub fn add_stereo_element(&mut self, element: StereoElement) -> Result<StereoElementId> {
        self.working.add_stereo_element(element)
    }

    pub fn replace_stereo_element(
        &mut self,
        id: StereoElementId,
        replacement: StereoElement,
    ) -> Result<StereoElement> {
        self.working.replace_stereo_element(id, replacement)
    }

    pub fn remove_stereo_element(&mut self, id: StereoElementId) -> Result<StereoElement> {
        self.working.remove_stereo_element(id)
    }

    pub fn add_stereo_group(&mut self, group: StereoGroup) -> Result<StereoGroupId> {
        self.working.add_stereo_group(group)
    }

    pub fn remove_stereo_group(&mut self, id: StereoGroupId) -> Result<StereoGroup> {
        self.working.remove_stereo_group(id)
    }

    pub fn append_stereo_group_tombstone(&mut self) -> Result<StereoGroupId> {
        self.working.append_stereo_group_tombstone()
    }

    pub fn finish(self) -> std::result::Result<Molecule, MoleculePublicationError> {
        publish_molecule(self.working)
    }
}

fn publish_molecule(
    mut molecule: Molecule,
) -> std::result::Result<Molecule, MoleculePublicationError> {
    if molecule.atom_count() == 0 {
        return Err(MoleculePublicationError::EmptyGraph);
    }
    validate_graph(&molecule).map_err(MoleculePublicationError::InvalidGraph)?;
    molecule
        .validate_connected()
        .map_err(MoleculePublicationError::DisconnectedGraph)?;
    validate_stereo(&molecule).map_err(MoleculePublicationError::InvalidStereo)?;
    canonicalize_represented_chemistry(&mut molecule).map_err(|error| {
        MoleculePublicationError::FormalChargeOutOfRange {
            atom: error.atom,
            charge: error.charge,
        }
    })?;
    molecule.clear_perception();
    Ok(molecule)
}

fn validate_graph(molecule: &Molecule) -> std::result::Result<(), GraphValidationError> {
    if molecule.graph.adjacency.len() != molecule.graph.atoms.len() {
        return Err(GraphValidationError::AdjacencySlotCount);
    }
    for (index, atom) in molecule.graph.atoms.iter().enumerate() {
        let atom_id = AtomId::new(index as u32);
        let adjacency = &molecule.graph.adjacency[index];
        if atom.is_none() && !adjacency.is_empty() {
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
            let Some(Some(bond)) = molecule.graph.bonds.get(bond_id.index()) else {
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
    for (index, bond) in molecule.graph.bonds.iter().enumerate() {
        let Some(bond) = bond else { continue };
        let bond_id = BondId::new(index as u32);
        for atom in [bond.a(), bond.b()] {
            if molecule.atom(atom).is_err() {
                return Err(GraphValidationError::InvalidBondEndpoint { bond: bond_id });
            }
            if !molecule.graph.adjacency[atom.index()].contains(&bond_id) {
                return Err(GraphValidationError::MissingAdjacencyEntry {
                    atom,
                    bond: bond_id,
                });
            }
        }
    }
    Ok(())
}

fn validate_stereo(molecule: &Molecule) -> std::result::Result<(), StereoPublicationError> {
    for (element_id, element) in molecule.stereo_elements() {
        molecule
            .validate_stereo_element_refs(element)
            .map_err(|_| StereoPublicationError::InvalidElementReference {
                element: element_id,
            })?;
        if let Some(group_id) = element.group {
            let group = molecule.stereo_group(group_id).map_err(|_| {
                StereoPublicationError::InvalidElementGroup {
                    element: element_id,
                }
            })?;
            if !group.members.contains(&element_id) {
                return Err(StereoPublicationError::InconsistentGroupMembership {
                    group: group_id,
                    element: element_id,
                });
            }
        }
    }
    for (group_id, group) in molecule.stereo_groups() {
        if group.members.is_empty() {
            return Err(StereoPublicationError::EmptyGroup { group: group_id });
        }
        let mut seen = std::collections::BTreeSet::new();
        for &element_id in &group.members {
            if !seen.insert(element_id) {
                return Err(StereoPublicationError::DuplicateGroupMember {
                    group: group_id,
                    element: element_id,
                });
            }
            let element = molecule.stereo_element(element_id).map_err(|_| {
                StereoPublicationError::InvalidGroupMember {
                    group: group_id,
                    element: element_id,
                }
            })?;
            if element.group != Some(group_id) {
                return Err(StereoPublicationError::InconsistentGroupMembership {
                    group: group_id,
                    element: element_id,
                });
            }
        }
    }
    Ok(())
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

        if let Some(atom) = molecule.graph.atoms[atom_id.index()].as_mut() {
            atom.formal_charge = formal_charge;
        }
        for (oxygen_id, bond_id) in oxo_bonds {
            if let Some(atom) = molecule.graph.atoms[oxygen_id.index()].as_mut() {
                atom.formal_charge = -1;
            }
            if let Some(bond) = molecule.graph.bonds[bond_id.index()].as_mut() {
                bond.order = BondOrder::Single;
            }
        }
    }
    if rewritten {
        molecule.clear_perception();
    }
    molecule.canonicalize_stored_stereo_elements();
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
    use crate::core::{
        DoubleBondOrientation, DoubleBondStereo, Element, StereoCarrier, StereoElement,
        StereoElementKind,
    };

    fn carbon() -> Atom {
        Atom::new(Element::from_atomic_number(6).expect("carbon"))
    }

    fn atom(symbol: &str) -> Atom {
        Atom::new(Element::from_symbol(symbol).expect("fixture element"))
    }

    #[test]
    fn graph_validation_errors_have_diagnostic_display_messages() {
        let error = GraphValidationError::InvalidAdjacencyBond {
            atom: AtomId::new(4),
            bond: BondId::new(9),
        };
        let message = error.to_string();
        assert_eq!(
            message,
            "atom atom4 adjacency references missing bond bond9"
        );
        assert!(!message.contains("InvalidAdjacencyBond"));
    }

    #[test]
    fn stereo_publication_errors_have_diagnostic_display_messages() {
        let error = StereoPublicationError::DuplicateGroupMember {
            group: StereoGroupId::new(2),
            element: StereoElementId::new(7),
        };
        let message = error.to_string();
        assert_eq!(
            message,
            "stereo group group2 contains stereo element stereo7 more than once"
        );
        assert!(!message.contains("DuplicateGroupMember"));
    }

    #[test]
    fn builder_rejects_disconnected_final_graph() {
        let mut builder = crate::core::MoleculeEditor::new();
        builder.add_atom(carbon()).expect("first atom");
        builder.add_atom(carbon()).expect("second atom");
        assert_eq!(
            builder.finish(),
            Err(MoleculePublicationError::DisconnectedGraph(
                MoleculeConnectivityError { components: 2 }
            ))
        );
    }

    #[test]
    fn builder_accepts_connected_graph() {
        let mut builder = crate::core::MoleculeEditor::new();
        let a = builder.add_atom(carbon()).expect("first atom");
        let b = builder.add_atom(carbon()).expect("second atom");
        builder
            .add_bond(a, b, BondOrder::Single)
            .expect("connecting bond");
        assert!(builder.finish().expect("connected molecule").is_connected());
    }

    #[test]
    fn builder_publishes_canonical_hypervalent_representation() {
        let mut builder = crate::core::MoleculeEditor::new();
        let chlorine = builder.add_atom(atom("Cl")).expect("chlorine");
        let oxo = builder.add_atom(atom("O")).expect("oxo oxygen");
        let hydroxyl = builder.add_atom(atom("O")).expect("hydroxyl oxygen");
        builder
            .add_bond(chlorine, oxo, BondOrder::Double)
            .expect("source-convention oxo bond");
        builder
            .add_bond(chlorine, hydroxyl, BondOrder::Single)
            .expect("hydroxyl bond");

        let molecule = builder.finish().expect("canonical publication");

        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
        assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
        let oxo_bond = molecule.bond_between(chlorine, oxo).unwrap().unwrap();
        assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
        assert_eq!(molecule.perception(), &super::super::Perception::default());
    }

    #[test]
    fn failed_editor_commit_leaves_original_unchanged() {
        let mut builder = crate::core::MoleculeEditor::new();
        let a = builder.add_atom(carbon()).expect("first atom");
        let b = builder.add_atom(carbon()).expect("second atom");
        let bond = builder
            .add_bond(a, b, BondOrder::Single)
            .expect("connecting bond");
        let molecule = builder.finish().expect("connected molecule");
        let before = molecule.clone();

        let mut editor = molecule.edit();
        editor.delete_bond(bond).expect("delete in working copy");
        assert_eq!(
            editor.finish(),
            Err(MoleculePublicationError::DisconnectedGraph(
                MoleculeConnectivityError { components: 2 }
            ))
        );
        assert_eq!(molecule, before);
    }

    #[test]
    fn editor_can_pass_through_temporarily_disconnected_state() {
        let mut builder = crate::core::MoleculeEditor::new();
        let a = builder.add_atom(carbon()).expect("first atom");
        let b = builder.add_atom(carbon()).expect("second atom");
        let c = builder.add_atom(carbon()).expect("third atom");
        let ab = builder
            .add_bond(a, b, BondOrder::Single)
            .expect("first bond");
        builder
            .add_bond(b, c, BondOrder::Single)
            .expect("second bond");
        let molecule = builder.finish().expect("connected molecule");

        let mut editor = molecule.edit();
        editor.delete_bond(ab).expect("temporary disconnect");
        editor
            .add_bond(a, c, BondOrder::Single)
            .expect("reconnect through another edge");
        let molecule = editor.finish().expect("connected final graph");
        assert!(molecule.is_connected());
    }

    #[test]
    fn editor_canonicalizes_before_transactional_publication() {
        let mut builder = crate::core::MoleculeEditor::new();
        let chlorine = builder.add_atom(atom("Cl")).expect("chlorine");
        let hydroxyl = builder.add_atom(atom("O")).expect("hydroxyl oxygen");
        builder
            .add_bond(chlorine, hydroxyl, BondOrder::Single)
            .expect("initial bond");
        let molecule = builder.finish().expect("initial molecule");

        let mut editor = molecule.edit();
        let oxo = editor.add_atom(atom("O")).expect("oxo oxygen");
        editor
            .add_bond(chlorine, oxo, BondOrder::Double)
            .expect("source-convention oxo bond");
        let molecule = editor.finish().expect("canonical edit publication");

        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
        assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
        let oxo_bond = molecule.bond_between(chlorine, oxo).unwrap().unwrap();
        assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
    }

    #[test]
    fn editor_canonicalizes_direct_bond_order_changes_before_publication() {
        let mut builder = crate::core::MoleculeEditor::new();
        let chlorine = builder.add_atom(atom("Cl")).unwrap();
        let anchor = builder.add_atom(atom("O")).unwrap();
        let oxo = builder.add_atom(atom("O")).unwrap();
        builder
            .add_bond(chlorine, anchor, BondOrder::Single)
            .unwrap();
        let oxo_bond = builder.add_bond(chlorine, oxo, BondOrder::Single).unwrap();
        let molecule = builder.finish().unwrap();

        let mut editor = molecule.edit();
        editor.bond_mut(oxo_bond).unwrap().order = BondOrder::Double;
        let molecule = editor.finish().expect("canonical edit publication");

        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
        assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
        assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
    }

    #[test]
    fn editor_recanonicalizes_stereo_references_after_topology_changes() {
        let mut builder = crate::core::MoleculeEditor::new();
        let left = builder.add_atom(carbon()).unwrap();
        let right = builder.add_atom(carbon()).unwrap();
        let hydrogen = builder.add_atom(atom("H")).unwrap();
        let chlorine = builder.add_atom(atom("Cl")).unwrap();
        let bromine = builder.add_atom(atom("Br")).unwrap();
        let double_bond = builder.add_bond(left, right, BondOrder::Double).unwrap();
        builder.add_bond(left, hydrogen, BondOrder::Single).unwrap();
        builder
            .add_bond(right, chlorine, BondOrder::Single)
            .unwrap();
        builder.add_bond(right, bromine, BondOrder::Single).unwrap();
        let molecule = builder.finish().unwrap();
        let mut editor = molecule.edit();
        let element = editor
            .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
                DoubleBondStereo {
                    bond: double_bond,
                    left,
                    right,
                    left_carrier: StereoCarrier::Atom(hydrogen),
                    right_carrier: StereoCarrier::Atom(chlorine),
                    orientation: Some(DoubleBondOrientation::Together),
                },
            )))
            .unwrap();
        let molecule = editor.finish().unwrap();

        let mut editor = molecule.edit();
        let fluorine = editor.add_atom(atom("F")).unwrap();
        editor.add_bond(left, fluorine, BondOrder::Single).unwrap();
        let molecule = editor.finish().expect("canonical edit publication");

        let StereoElementKind::DoubleBond(stereo) = &molecule.stereo_element(element).unwrap().kind
        else {
            panic!("expected double-bond stereo");
        };
        assert_eq!(stereo.left_carrier, StereoCarrier::Atom(fluorine));
        assert_eq!(stereo.right_carrier, StereoCarrier::Atom(chlorine));
        assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Opposite));
    }

    #[test]
    fn failed_editor_canonicalization_leaves_the_published_molecule_unchanged() {
        let mut builder = crate::core::MoleculeEditor::new();
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
        let molecule = builder.finish().unwrap();
        let before = molecule.clone();

        let mut editor = molecule.edit();
        for bond in oxo_bonds {
            editor.bond_mut(bond).unwrap().order = BondOrder::Double;
        }
        assert_eq!(
            editor.finish(),
            Err(MoleculePublicationError::FormalChargeOutOfRange {
                atom: chlorine,
                charge: 128,
            })
        );
        assert_eq!(molecule, before);
    }
}
