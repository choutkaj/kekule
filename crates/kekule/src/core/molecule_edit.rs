use std::fmt;

use crate::properties::{Properties, PropertyKey, PropertyTable, PropertyValue};

use super::{
    Atom, AtomId, Bond, BondId, BondOrder, Graph, Molecule, Perception, Result, RingMembership,
    RingSet, StereoDescriptor, StereoElement, StereoElementId, StereoGroup, StereoGroupId,
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
    AtomPropertySlotCount { expected: usize, actual: usize },
    BondPropertySlotCount { expected: usize, actual: usize },
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
            Self::AtomPropertySlotCount { expected, actual } => write!(formatter,
                "atom property storage requires {expected} slots, but has {actual}"),
            Self::BondPropertySlotCount { expected, actual } => write!(formatter,
                "bond property storage requires {expected} slots, but has {actual}"),
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

    /// Moves this molecule into an editor without cloning its graph or properties.
    pub fn to_editor(self) -> MoleculeEditor {
        MoleculeEditor { working: self }
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
/// [`Self::finish`] validates graph references, connectedness, stereochemistry,
/// and canonical represented chemistry before publishing a [`Molecule`]. A
/// failed finish returns structured errors and publishes nothing.
/// Use [`Self::try_finish`] when a failed draft must remain available for repair.
///
/// | Task | Operations |
/// | --- | --- |
/// | Start | [`Self::new`], [`Molecule::edit`], [`Molecule::to_editor`] |
/// | Inspect | [`Self::atoms`], [`Self::bonds`], [`Self::neighbors`], [`Self::connected_components`] |
/// | Change graph | [`Self::replace_atom`], [`Self::replace_bond`], [`Self::delete_atoms`], [`Self::retain_atoms`] |
/// | Combine fragments | [`Self::append_molecule`] with returned ID mappings |
/// | Annotate | [`Self::set_atom_property`], [`Self::set_atom_properties`], [`Self::set_atom_property_column`] and bond counterparts |
/// | Edit stereo | [`Self::replace_stereo_element`], [`Self::replace_stereo_group`] |
///
/// Structural changes invalidate perception and clear owner properties while
/// retaining annotations on surviving atom/bond identities. Generic property
/// edits leave perception intact. Atom and bond IDs remain stable until
/// [`Self::clear`]; iteration skips deleted slots.
///
/// An unfinished editor cannot be passed to algorithms requiring a published molecule:
/// ```compile_fail
/// use kekule::core::{Molecule, MoleculeEditor};
/// fn published(editor: &MoleculeEditor) -> &Molecule { editor }
/// ```
///
/// # Example
///
/// ```
/// use kekule::core::{Atom, BondOrder, Element, MoleculeEditor};
///
/// let mut editor = MoleculeEditor::new();
/// let carbon = editor.add_atom(Atom::new(Element::from_symbol("C").unwrap()))?;
/// let oxygen = editor.add_atom(Atom::new(Element::from_symbol("O").unwrap()))?;
/// editor.add_bond(carbon, oxygen, BondOrder::Single)?;
/// let molecule = editor.finish()?;
///
/// assert_eq!(molecule.atom_count(), 2);
/// assert_eq!(molecule.bond_count(), 1);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeEditor {
    pub(super) working: Molecule,
}

/// A failed recoverable finish, retaining the original, unmodified editing state.
#[derive(Debug)]
pub struct MoleculeFinishError {
    error: MoleculePublicationError,
    editor: Box<MoleculeEditor>,
}

impl MoleculeFinishError {
    pub const fn error(&self) -> &MoleculePublicationError {
        &self.error
    }
    pub fn editor(&self) -> &MoleculeEditor {
        &self.editor
    }
    pub fn to_editor(self) -> MoleculeEditor {
        *self.editor
    }
}

impl fmt::Display for MoleculeFinishError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.error.fmt(f)
    }
}

impl std::error::Error for MoleculeFinishError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.error)
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

    pub fn atom_count(&self) -> usize {
        self.working.atom_count()
    }

    pub fn bond_count(&self) -> usize {
        self.working.bond_count()
    }

    pub fn formal_charge(&self) -> i64 {
        self.working.formal_charge()
    }

    pub fn atom(&self, id: AtomId) -> Result<&Atom> {
        self.working.atom(id)
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.working.atoms()
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.working.atom_ids()
    }

    pub fn bond(&self, id: BondId) -> Result<&Bond> {
        self.working.bond(id)
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)> {
        self.working.bonds()
    }

    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        self.working.bond_ids()
    }

    pub fn neighbors(&self, id: AtomId) -> Result<impl Iterator<Item = AtomId> + '_> {
        self.working.neighbors(id)
    }

    pub fn incident_bonds(&self, id: AtomId) -> Result<impl Iterator<Item = (BondId, &Bond)> + '_> {
        self.working.incident_bonds(id)
    }

    pub fn bond_between(&self, a: AtomId, b: AtomId) -> Result<Option<BondId>> {
        self.working.bond_between(a, b)
    }

    pub const fn properties(&self) -> &Properties {
        self.working.properties()
    }

    pub const fn atom_properties(&self) -> &PropertyTable {
        self.working.atom_properties()
    }

    pub const fn bond_properties(&self) -> &PropertyTable {
        self.working.bond_properties()
    }

    pub fn atom_property(&self, id: AtomId, key: &PropertyKey) -> Result<Option<PropertyValue>> {
        self.working.atom_property(id, key)
    }

    pub fn bond_property(&self, id: BondId, key: &PropertyKey) -> Result<Option<PropertyValue>> {
        self.working.bond_property(id, key)
    }

    pub fn set_atom_property(
        &mut self,
        id: AtomId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<()> {
        self.working.set_atom_property(id, key, value)
    }

    pub fn set_bond_property(
        &mut self,
        id: BondId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<()> {
        self.working.set_bond_property(id, key, value)
    }

    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>> {
        self.working.insert_property(key, value)
    }

    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.working.remove_property(key)
    }

    pub fn clear_properties(&mut self) {
        self.working.clear_properties()
    }

    pub fn stereo_element(&self, id: StereoElementId) -> Result<&StereoElement> {
        self.working.stereo_element(id)
    }

    pub fn stereo_elements(&self) -> impl Iterator<Item = (StereoElementId, &StereoElement)> {
        self.working.stereo_elements()
    }

    pub fn stereo_element_ids(&self) -> impl Iterator<Item = StereoElementId> + '_ {
        self.working.stereo_element_ids()
    }

    pub fn stereo_group(&self, id: StereoGroupId) -> Result<&StereoGroup> {
        self.working.stereo_group(id)
    }

    pub fn stereo_groups(&self) -> impl Iterator<Item = (StereoGroupId, &StereoGroup)> {
        self.working.stereo_groups()
    }

    pub fn stereo_group_slot_count(&self) -> usize {
        self.working.stereo_group_slot_count()
    }

    pub fn stereo_group_slots(
        &self,
    ) -> impl ExactSizeIterator<Item = (StereoGroupId, Option<&StereoGroup>)> + DoubleEndedIterator + '_
    {
        self.working.stereo_group_slots()
    }

    pub fn perception(&self) -> &Perception {
        self.working.perception()
    }

    pub fn implicit_hydrogens(&self, atom: AtomId) -> Result<Option<u8>> {
        self.working.implicit_hydrogens(atom)
    }

    pub fn atom_is_aromatic(&self, atom: AtomId) -> Result<Option<bool>> {
        self.working.atom_is_aromatic(atom)
    }

    pub fn bond_is_aromatic(&self, bond: BondId) -> Result<Option<bool>> {
        self.working.bond_is_aromatic(bond)
    }

    pub fn cip_descriptor(&self, element: StereoElementId) -> Result<Option<StereoDescriptor>> {
        self.working.cip_descriptor(element)
    }

    pub fn ring_membership(&self) -> Option<&RingMembership> {
        self.working.ring_membership()
    }

    pub fn ring_set(&self) -> Option<&RingSet> {
        self.working.ring_set()
    }

    pub fn clear_perception(&mut self) {
        self.working.clear_perception()
    }

    pub(crate) fn working(&self) -> &Molecule {
        &self.working
    }

    pub(crate) fn working_mut(&mut self) -> &mut Molecule {
        &mut self.working
    }
    /// Returns mutable represented atom state in this private working copy.
    /// Obtaining mutable atom access invalidates perception and owner properties
    /// immediately, even if the guard is forgotten. Use [`Self::replace_atom`]
    /// when an identical replacement should preserve cached state.
    pub fn atom_mut(&mut self, atom: AtomId) -> Result<super::AtomMut<'_>> {
        self.working.atom_mut(atom)
    }

    /// Borrows a bond with checked order editing through [`super::BondMut::set_order`].
    /// Change connectivity with [`Self::set_bond_endpoints`] or [`Self::replace_bond`].
    ///
    /// ```compile_fail
    /// use kekule::core::{AtomId, Bond, BondId, BondOrder, MoleculeEditor};
    /// fn bypass(editor: &mut MoleculeEditor, id: BondId) {
    ///     *editor.bond_mut(id).unwrap() = Bond::new(AtomId::new(0), AtomId::new(99), BondOrder::Single);
    /// }
    /// ```
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

    /// Publishes represented chemistry, clearing draft perception. Install any
    /// reconstructed perception on the resulting [`Molecule::install_perception`].
    ///
    /// ```compile_fail
    /// use kekule::core::{MoleculeEditor, Perception};
    /// let mut editor = MoleculeEditor::new();
    /// editor.install_perception(Perception::default()).unwrap();
    /// ```
    pub fn finish(self) -> std::result::Result<Molecule, MoleculePublicationError> {
        publish_molecule(self.working)
    }

    /// Checks whether a snapshot would publish successfully, leaving this editor
    /// unchanged. This clones the working state to check canonicalization too.
    pub fn validate(&self) -> std::result::Result<(), MoleculePublicationError> {
        publish_molecule(self.working.clone()).map(|_| ())
    }

    /// Publishes, or returns the original editor for repair through
    /// [`MoleculeFinishError::to_editor`]. Keeps a cloned rollback snapshot;
    /// use [`Self::finish`] when recovery is unnecessary.
    pub fn try_finish(self) -> std::result::Result<Molecule, MoleculeFinishError> {
        let snapshot = self.clone();
        publish_molecule(self.working).map_err(|error| MoleculeFinishError {
            error,
            editor: Box::new(snapshot),
        })
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
    if molecule.atom_properties().len() != molecule.graph.atoms.len() {
        return Err(GraphValidationError::AtomPropertySlotCount {
            expected: molecule.graph.atoms.len(),
            actual: molecule.atom_properties().len(),
        });
    }
    if molecule.bond_properties().len() != molecule.graph.bonds.len() {
        return Err(GraphValidationError::BondPropertySlotCount {
            expected: molecule.graph.bonds.len(),
            actual: molecule.bond_properties().len(),
        });
    }
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
    fn publication_checks_property_domains_and_recovery_preserves_malformed_draft() {
        for bonds in [false, true] {
            let mut editor = MoleculeEditor::new();
            editor.add_atom(carbon()).unwrap();
            if bonds {
                editor.working.properties.resize_bonds(1);
            } else {
                editor.working.properties.resize_atoms(0);
            }
            let before = format!("{editor:#?}");
            let error = editor.try_finish().unwrap_err();
            let expected = if bonds {
                GraphValidationError::BondPropertySlotCount {
                    expected: 0,
                    actual: 1,
                }
            } else {
                GraphValidationError::AtomPropertySlotCount {
                    expected: 1,
                    actual: 0,
                }
            };
            assert_eq!(
                error.error(),
                &MoleculePublicationError::InvalidGraph(expected)
            );
            assert_eq!(format!("{:#?}", error.editor()), before);
        }
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
        editor
            .bond_mut(oxo_bond)
            .unwrap()
            .set_order(BondOrder::Double);
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
            editor.bond_mut(bond).unwrap().set_order(BondOrder::Double);
        }
        let draft = format!("{editor:#?}");
        let failure = editor.clone().try_finish().unwrap_err();
        assert_eq!(
            failure.error(),
            &MoleculePublicationError::FormalChargeOutOfRange {
                atom: chlorine,
                charge: 128,
            }
        );
        assert_eq!(format!("{:#?}", failure.editor()), draft);
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
