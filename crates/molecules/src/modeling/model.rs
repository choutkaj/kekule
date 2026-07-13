use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::bio::{AtomSite, AtomSiteId, BioHierarchy, MacroMolecule};
use crate::core::{Atom, AtomId, Bond, BondId, ConformerId, Molecule, Point3, PropMap};
use crate::small::SmallMolecule;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MoleculeInstanceId(u32);

impl MoleculeInstanceId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl fmt::Display for MoleculeInstanceId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "molecule{}", self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InstanceAtomId {
    molecule: MoleculeInstanceId,
    atom: AtomId,
}

impl InstanceAtomId {
    pub const fn new(molecule: MoleculeInstanceId, atom: AtomId) -> Self {
        Self { molecule, atom }
    }

    pub const fn molecule(self) -> MoleculeInstanceId {
        self.molecule
    }

    pub const fn atom(self) -> AtomId {
        self.atom
    }
}

impl fmt::Display for InstanceAtomId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.molecule, self.atom)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InstanceBondId {
    molecule: MoleculeInstanceId,
    bond: BondId,
}

impl InstanceBondId {
    pub const fn new(molecule: MoleculeInstanceId, bond: BondId) -> Self {
        Self { molecule, bond }
    }

    pub const fn molecule(self) -> MoleculeInstanceId {
        self.molecule
    }

    pub const fn bond(self) -> BondId {
        self.bond
    }
}

impl fmt::Display for InstanceBondId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.molecule, self.bond)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ModelAtomIndex(u32);

impl ModelAtomIndex {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum MoleculeRole {
    Polymer,
    Branched,
    NonPolymer,
    Solvent,
    Ion,
    Ligand,
    Cofactor,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MoleculeInstanceMetadata {
    roles: BTreeSet<MoleculeRole>,
    props: PropMap,
}

impl MoleculeInstanceMetadata {
    pub fn roles(&self) -> &BTreeSet<MoleculeRole> {
        &self.roles
    }

    pub fn insert_role(&mut self, role: MoleculeRole) -> bool {
        self.roles.insert(role)
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }
}

#[derive(Debug, Clone, PartialEq)]
enum MoleculeInstancePayload {
    Small(SmallMolecule),
    Macro(MacroMolecule),
}

impl MoleculeInstancePayload {
    fn graph(&self) -> &Molecule {
        match self {
            Self::Small(molecule) => molecule.graph(),
            Self::Macro(molecule) => molecule.graph(),
        }
    }

    fn without_conformers(self) -> Self {
        match self {
            Self::Small(molecule) => Self::Small(molecule.without_conformers()),
            Self::Macro(molecule) => Self::Macro(molecule.without_conformers()),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeInstance {
    id: MoleculeInstanceId,
    payload: MoleculeInstancePayload,
    metadata: MoleculeInstanceMetadata,
}

impl MoleculeInstance {
    pub const fn id(&self) -> MoleculeInstanceId {
        self.id
    }

    pub fn graph(&self) -> &Molecule {
        self.payload.graph()
    }

    pub fn small_molecule(&self) -> Option<&SmallMolecule> {
        match &self.payload {
            MoleculeInstancePayload::Small(molecule) => Some(molecule),
            MoleculeInstancePayload::Macro(_) => None,
        }
    }

    pub fn macro_molecule(&self) -> Option<&MacroMolecule> {
        match &self.payload {
            MoleculeInstancePayload::Macro(molecule) => Some(molecule),
            MoleculeInstancePayload::Small(_) => None,
        }
    }

    pub fn hierarchy(&self) -> Option<&BioHierarchy> {
        self.macro_molecule().map(MacroMolecule::hierarchy)
    }

    pub fn bio_hierarchy(&self) -> Option<InstanceBioHierarchy<'_>> {
        self.hierarchy().map(|hierarchy| InstanceBioHierarchy {
            molecule: self.id,
            hierarchy,
        })
    }

    pub fn roles(&self) -> &BTreeSet<MoleculeRole> {
        self.metadata.roles()
    }

    pub fn has_role(&self, role: MoleculeRole) -> bool {
        self.roles().contains(&role)
    }

    pub fn props(&self) -> &PropMap {
        self.metadata.props()
    }

    pub const fn qualify_atom(&self, atom: AtomId) -> InstanceAtomId {
        InstanceAtomId::new(self.id, atom)
    }

    pub const fn qualify_bond(&self, bond: BondId) -> InstanceBondId {
        InstanceBondId::new(self.id, bond)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct InstanceBioHierarchy<'a> {
    molecule: MoleculeInstanceId,
    hierarchy: &'a BioHierarchy,
}

impl InstanceBioHierarchy<'_> {
    pub const fn molecule(&self) -> MoleculeInstanceId {
        self.molecule
    }

    pub fn hierarchy(&self) -> &BioHierarchy {
        self.hierarchy
    }

    pub fn atom_for_site(&self, site: AtomSiteId) -> Result<InstanceAtomId, ModelError> {
        let site = self
            .hierarchy
            .atom_site(site)
            .map_err(|_| ModelError::InvalidAtomSiteId(site))?;
        Ok(InstanceAtomId::new(self.molecule, site.atom))
    }

    pub fn atom_site_for_atom(&self, atom: InstanceAtomId) -> Option<&AtomSite> {
        (atom.molecule == self.molecule)
            .then(|| self.hierarchy.atom_site_for_atom(atom.atom))
            .flatten()
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ModelTopology {
    molecules: Vec<MoleculeInstance>,
    atom_order: Vec<InstanceAtomId>,
    atom_indexes: BTreeMap<InstanceAtomId, ModelAtomIndex>,
}

impl ModelTopology {
    pub fn molecule(&self, id: MoleculeInstanceId) -> Result<&MoleculeInstance, ModelError> {
        self.molecules
            .get(id.index())
            .ok_or(ModelError::InvalidMoleculeInstanceId(id))
    }

    pub fn molecules(&self) -> impl Iterator<Item = (MoleculeInstanceId, &MoleculeInstance)> {
        self.molecules
            .iter()
            .map(|molecule| (molecule.id, molecule))
    }

    pub fn molecule_count(&self) -> usize {
        self.molecules.len()
    }

    pub fn atom(&self, id: InstanceAtomId) -> Result<&Atom, ModelError> {
        self.molecule(id.molecule)?
            .graph()
            .atom(id.atom)
            .map_err(|_| ModelError::InvalidAtomId(id))
    }

    pub fn bond(&self, id: InstanceBondId) -> Result<&Bond, ModelError> {
        self.molecule(id.molecule)?
            .graph()
            .bond(id.bond)
            .map_err(|_| ModelError::InvalidBondId(id))
    }

    pub fn atoms(&self) -> impl Iterator<Item = (InstanceAtomId, &Atom)> {
        self.molecules.iter().flat_map(|molecule| {
            molecule
                .graph()
                .atoms()
                .map(move |(atom, payload)| (molecule.qualify_atom(atom), payload))
        })
    }

    pub fn bonds(&self) -> impl Iterator<Item = (InstanceBondId, &Bond)> {
        self.molecules.iter().flat_map(|molecule| {
            molecule
                .graph()
                .bonds()
                .map(move |(bond, payload)| (molecule.qualify_bond(bond), payload))
        })
    }

    pub fn atom_ids(&self) -> &[InstanceAtomId] {
        &self.atom_order
    }

    pub fn atom_index(&self, atom: InstanceAtomId) -> Option<ModelAtomIndex> {
        self.atom_indexes.get(&atom).copied()
    }

    pub fn atom_id(&self, index: ModelAtomIndex) -> Option<InstanceAtomId> {
        self.atom_order.get(index.index()).copied()
    }

    pub fn molecule_for_atom(&self, atom: InstanceAtomId) -> Option<&MoleculeInstance> {
        self.atom(atom).ok()?;
        self.molecule(atom.molecule).ok()
    }

    pub fn implicit_hydrogens(&self, atom: InstanceAtomId) -> Result<Option<u8>, ModelError> {
        self.atom(atom)?;
        self.molecule(atom.molecule)
            .expect("validated molecule instance")
            .graph()
            .implicit_hydrogens(atom.atom)
            .map_err(|_| ModelError::InvalidAtomId(atom))
    }

    pub fn atom_is_aromatic(&self, atom: InstanceAtomId) -> Result<Option<bool>, ModelError> {
        self.atom(atom)?;
        self.molecule(atom.molecule)
            .expect("validated molecule instance")
            .graph()
            .atom_is_aromatic(atom.atom)
            .map_err(|_| ModelError::InvalidAtomId(atom))
    }

    pub fn bond_is_aromatic(&self, bond: InstanceBondId) -> Result<Option<bool>, ModelError> {
        self.bond(bond)?;
        self.molecule(bond.molecule)
            .expect("validated molecule instance")
            .graph()
            .bond_is_aromatic(bond.bond)
            .map_err(|_| ModelError::InvalidBondId(bond))
    }
}

#[derive(Debug, PartialEq)]
struct ModelDefinition {
    topology: ModelTopology,
}

/// Opaque identity of a model's immutable topology and molecule-instance definition.
///
/// Clones and coordinate updates preserve the key. Independently built models
/// receive distinct keys even when their topology contents are structurally equal.
#[derive(Clone)]
pub struct ModelDefinitionKey(Arc<ModelDefinition>);

impl fmt::Debug for ModelDefinitionKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("ModelDefinitionKey(..)")
    }
}

impl PartialEq for ModelDefinitionKey {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl Eq for ModelDefinitionKey {}

impl Hash for ModelDefinitionKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        Arc::as_ptr(&self.0).hash(state);
    }
}

#[derive(Clone)]
pub struct MolecularModel {
    definition: ModelDefinitionKey,
    positions: Vec<Point3>,
}

impl fmt::Debug for MolecularModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MolecularModel")
            .field("topology", &self.definition.0.topology)
            .field("positions", &self.positions)
            .finish()
    }
}

impl PartialEq for MolecularModel {
    fn eq(&self, other: &Self) -> bool {
        self.definition.0.as_ref() == other.definition.0.as_ref()
            && self.positions == other.positions
    }
}

impl MolecularModel {
    pub fn builder() -> MolecularModelBuilder {
        MolecularModelBuilder::new()
    }

    pub fn from_small_molecule(
        molecule: &SmallMolecule,
        conformer: ConformerId,
    ) -> Result<Self, ModelBuildError> {
        let mut builder = Self::builder();
        builder.add_small_molecule(molecule, conformer)?;
        builder.build()
    }

    pub fn from_macro_molecule(
        molecule: &MacroMolecule,
        conformer: ConformerId,
    ) -> Result<Self, ModelBuildError> {
        let mut builder = Self::builder();
        builder.add_macro_molecule(molecule, conformer)?;
        builder.build()
    }

    pub fn topology(&self) -> &ModelTopology {
        &self.definition.0.topology
    }

    /// Returns the identity of this model's immutable definition.
    pub fn definition_key(&self) -> &ModelDefinitionKey {
        &self.definition
    }

    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    pub fn positions(&self) -> &[Point3] {
        &self.positions
    }

    pub fn position(&self, atom: InstanceAtomId) -> Result<Point3, PositionError> {
        let index = self
            .topology()
            .atom_index(atom)
            .ok_or(PositionError::InvalidAtomId(atom))?;
        Ok(self.positions[index.index()])
    }

    pub fn position_at(&self, index: ModelAtomIndex) -> Result<Point3, PositionError> {
        self.positions
            .get(index.index())
            .copied()
            .ok_or(PositionError::InvalidAtomIndex(index))
    }

    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Point3,
    ) -> Result<(), PositionError> {
        let index = self
            .topology()
            .atom_index(atom)
            .ok_or(PositionError::InvalidAtomId(atom))?;
        if !point_is_finite(position) {
            return Err(PositionError::NonFinitePosition { atom });
        }
        self.positions[index.index()] = position;
        Ok(())
    }

    pub fn set_positions(&mut self, positions: &[Point3]) -> Result<(), PositionError> {
        if positions.len() != self.positions.len() {
            return Err(PositionError::PositionCountMismatch {
                expected: self.positions.len(),
                actual: positions.len(),
            });
        }
        for (index, point) in positions.iter().copied().enumerate() {
            if !point_is_finite(point) {
                let atom = self.topology().atom_order[index];
                return Err(PositionError::NonFinitePosition { atom });
            }
        }
        self.positions.clone_from_slice(positions);
        Ok(())
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MolecularModelBuilder {
    topology: ModelTopology,
    positions: Vec<Point3>,
}

impl MolecularModelBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_small_molecule(
        &mut self,
        molecule: &SmallMolecule,
        conformer: ConformerId,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        self.add_small_molecule_with_metadata(
            molecule,
            conformer,
            MoleculeInstanceMetadata::default(),
        )
    }

    pub fn add_small_molecule_with_metadata(
        &mut self,
        molecule: &SmallMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        self.add_payload(
            MoleculeInstancePayload::Small(molecule.clone()),
            conformer,
            metadata,
        )
    }

    pub fn add_macro_molecule(
        &mut self,
        molecule: &MacroMolecule,
        conformer: ConformerId,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        self.add_macro_molecule_with_metadata(
            molecule,
            conformer,
            MoleculeInstanceMetadata::default(),
        )
    }

    pub fn add_macro_molecule_with_metadata(
        &mut self,
        molecule: &MacroMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        self.add_payload(
            MoleculeInstancePayload::Macro(molecule.clone()),
            conformer,
            metadata,
        )
    }

    pub fn build(self) -> Result<MolecularModel, ModelBuildError> {
        if self.topology.molecules.is_empty() {
            return Err(ModelBuildError::EmptyModel);
        }
        Ok(MolecularModel {
            definition: ModelDefinitionKey(Arc::new(ModelDefinition {
                topology: self.topology,
            })),
            positions: self.positions,
        })
    }

    fn add_payload(
        &mut self,
        payload: MoleculeInstancePayload,
        conformer_id: ConformerId,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        let mut staged = self.clone();
        let id = staged.add_payload_staged(payload, conformer_id, metadata)?;
        *self = staged;
        Ok(id)
    }

    fn add_payload_staged(
        &mut self,
        payload: MoleculeInstancePayload,
        conformer_id: ConformerId,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        let graph = payload.graph();
        if graph.atom_count() == 0 {
            return Err(ModelBuildError::EmptyMolecule);
        }
        let conformer = graph
            .conformer(conformer_id)
            .map_err(|_| ModelBuildError::InvalidConformerId(conformer_id))?;
        let mut source_positions = BTreeMap::new();
        for (atom, _) in graph.atoms() {
            let point = conformer
                .position(atom)
                .ok_or(ModelBuildError::MissingPosition { atom })?;
            if !point_is_finite(point) {
                return Err(ModelBuildError::NonFinitePosition { atom });
            }
            source_positions.insert(atom, point);
        }

        let id = MoleculeInstanceId::new(self.topology.molecules.len() as u32);
        for (atom, _) in graph.atoms() {
            let qualified = InstanceAtomId::new(id, atom);
            let index = ModelAtomIndex::new(self.positions.len() as u32);
            self.topology.atom_indexes.insert(qualified, index);
            self.topology.atom_order.push(qualified);
            self.positions.push(source_positions[&atom]);
        }
        self.topology.molecules.push(MoleculeInstance {
            id,
            payload: payload.without_conformers(),
            metadata,
        });
        Ok(id)
    }
}

pub(crate) fn point_is_finite(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelError {
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    InvalidAtomId(InstanceAtomId),
    InvalidBondId(InstanceBondId),
    InvalidAtomSiteId(AtomSiteId),
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMoleculeInstanceId(id) => write!(f, "invalid molecule instance: {id}"),
            Self::InvalidAtomId(id) => write!(f, "invalid molecule-instance atom: {id}"),
            Self::InvalidBondId(id) => write!(f, "invalid molecule-instance bond: {id}"),
            Self::InvalidAtomSiteId(id) => write!(f, "invalid atom-site id: {}", id.raw()),
        }
    }
}

impl std::error::Error for ModelError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PositionError {
    InvalidAtomId(InstanceAtomId),
    InvalidAtomIndex(ModelAtomIndex),
    PositionCountMismatch { expected: usize, actual: usize },
    NonFinitePosition { atom: InstanceAtomId },
}

impl fmt::Display for PositionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(atom) => write!(f, "invalid model atom id: {atom}"),
            Self::InvalidAtomIndex(index) => {
                write!(f, "invalid dense model atom index: {}", index.raw())
            }
            Self::PositionCountMismatch { expected, actual } => write!(
                f,
                "model requires {expected} positions, but received {actual}"
            ),
            Self::NonFinitePosition { atom } => {
                write!(f, "model position for atom {atom} is not finite")
            }
        }
    }
}

impl std::error::Error for PositionError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelBuildError {
    EmptyModel,
    EmptyMolecule,
    InvalidConformerId(ConformerId),
    MissingPosition { atom: AtomId },
    NonFinitePosition { atom: AtomId },
}

impl fmt::Display for ModelBuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyModel => write!(f, "a molecular model must contain at least one molecule"),
            Self::EmptyMolecule => write!(f, "a model molecule must contain at least one atom"),
            Self::InvalidConformerId(id) => write!(f, "invalid source conformer id: {id}"),
            Self::MissingPosition { atom } => {
                write!(f, "source conformer has no position for atom {atom}")
            }
            Self::NonFinitePosition { atom } => {
                write!(f, "source conformer position for atom {atom} is not finite")
            }
        }
    }
}

impl std::error::Error for ModelBuildError {}
