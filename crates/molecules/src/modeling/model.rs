use std::collections::BTreeMap;
use std::fmt;

use crate::core::{
    AtomId, AxisStereo, BondId, ConformerId, DoubleBondStereo, Molecule, MoleculeError, Point3,
    PropMap, StereoBondMark, StereoCarrier, StereoElement, StereoElementKind, StereoGroup,
    TetrahedralStereo,
};
use crate::small::SmallMolecule;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
/// Stable identifier for one source-molecule component in a [`MolecularModel`].
pub struct ComponentId(u32);

impl ComponentId {
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

impl fmt::Display for ComponentId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "component{}", self.0)
    }
}

#[derive(Debug, Clone, PartialEq)]
/// Read-only membership and source properties for one model component.
pub struct Component {
    id: ComponentId,
    atoms: Vec<AtomId>,
    props: PropMap,
}

impl Component {
    pub const fn id(&self) -> ComponentId {
        self.id
    }

    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// Source-to-model ID mapping returned when a component is added to a builder.
pub struct ComponentMapping {
    component: ComponentId,
    atoms: BTreeMap<AtomId, AtomId>,
    bonds: BTreeMap<BondId, BondId>,
}

impl ComponentMapping {
    pub const fn component(&self) -> ComponentId {
        self.component
    }

    pub fn atom(&self, source: AtomId) -> Option<AtomId> {
        self.atoms.get(&source).copied()
    }

    pub fn bond(&self, source: BondId) -> Option<BondId> {
        self.bonds.get(&source).copied()
    }

    pub fn atoms(&self) -> &BTreeMap<AtomId, AtomId> {
        &self.atoms
    }

    pub fn bonds(&self) -> &BTreeMap<BondId, BondId> {
        &self.bonds
    }
}

#[derive(Debug, Clone, PartialEq)]
/// A fixed topology with exactly one complete Cartesian coordinate set.
///
/// Positions use angstroms. The topology contains no conformers; coordinates
/// are owned exclusively by this type and may be updated without changing
/// topology or component membership.
pub struct MolecularModel {
    topology: Molecule,
    positions: Vec<Point3>,
    components: Vec<Component>,
}

impl MolecularModel {
    pub fn builder() -> MolecularModelBuilder {
        MolecularModelBuilder::new()
    }

    pub fn from_conformer(
        molecule: &SmallMolecule,
        conformer: ConformerId,
    ) -> Result<Self, ModelBuildError> {
        let mut builder = Self::builder();
        builder.add_component(molecule, conformer)?;
        builder.build()
    }

    pub fn topology(&self) -> &Molecule {
        &self.topology
    }

    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    pub fn positions(&self) -> &[Point3] {
        &self.positions
    }

    pub fn position(&self, atom: AtomId) -> Result<Point3, PositionError> {
        self.topology
            .atom(atom)
            .map_err(|_| PositionError::InvalidAtomId(atom))?;
        Ok(self.positions[atom.index()])
    }

    pub fn set_position(&mut self, atom: AtomId, position: Point3) -> Result<(), PositionError> {
        self.topology
            .atom(atom)
            .map_err(|_| PositionError::InvalidAtomId(atom))?;
        if !point_is_finite(position) {
            return Err(PositionError::NonFinitePosition { atom });
        }
        self.positions[atom.index()] = position;
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
                return Err(PositionError::NonFinitePosition {
                    atom: AtomId::new(index as u32),
                });
            }
        }
        self.positions.clone_from_slice(positions);
        Ok(())
    }

    pub fn component(&self, id: ComponentId) -> Result<&Component, ModelError> {
        self.components
            .get(id.index())
            .ok_or(ModelError::InvalidComponentId(id))
    }

    pub fn components(&self) -> impl Iterator<Item = (ComponentId, &Component)> {
        self.components
            .iter()
            .map(|component| (component.id, component))
    }

    pub fn component_for_atom(&self, atom: AtomId) -> Option<&Component> {
        self.components
            .iter()
            .find(|component| component.atoms.contains(&atom))
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
/// Transactional builder for a fixed-topology [`MolecularModel`].
pub struct MolecularModelBuilder {
    topology: Molecule,
    positions: Vec<Point3>,
    components: Vec<Component>,
}

impl MolecularModelBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_component(
        &mut self,
        molecule: &SmallMolecule,
        conformer_id: ConformerId,
    ) -> Result<ComponentMapping, ModelBuildError> {
        let mut staged = self.clone();
        let mapping = staged.add_component_staged(molecule, conformer_id)?;
        *self = staged;
        Ok(mapping)
    }

    pub fn build(self) -> Result<MolecularModel, ModelBuildError> {
        if self.components.is_empty() {
            return Err(ModelBuildError::EmptyModel);
        }
        Ok(MolecularModel {
            topology: self.topology,
            positions: self.positions,
            components: self.components,
        })
    }

    fn add_component_staged(
        &mut self,
        molecule: &SmallMolecule,
        conformer_id: ConformerId,
    ) -> Result<ComponentMapping, ModelBuildError> {
        let graph = molecule.graph();
        if graph.atom_count() == 0 {
            return Err(ModelBuildError::EmptyComponent);
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

        let component = ComponentId::new(self.components.len() as u32);
        let mut atom_map = BTreeMap::new();
        let mut component_atoms = Vec::with_capacity(graph.atom_count());
        for (source_id, atom) in graph.atoms() {
            let model_id = self.topology.add_atom(atom.clone());
            debug_assert_eq!(model_id.index(), self.positions.len());
            self.positions.push(source_positions[&source_id]);
            atom_map.insert(source_id, model_id);
            component_atoms.push(model_id);
        }

        let mut bond_map = BTreeMap::new();
        for (source_id, bond) in graph.bonds() {
            let (source_a, source_b) = bond.endpoints();
            let a = mapped_atom(&atom_map, source_a)?;
            let b = mapped_atom(&atom_map, source_b)?;
            let model_id = self.topology.add_bond(a, b, bond.order)?;
            {
                let mut copied = self.topology.bond_mut(model_id)?;
                copied.aromatic = bond.aromatic;
                copied.props = bond.props.clone();
            }
            bond_map.insert(source_id, model_id);
        }

        copy_stereo(graph, &mut self.topology, &atom_map, &bond_map)?;

        self.components.push(Component {
            id: component,
            atoms: component_atoms,
            props: graph.props().clone(),
        });
        Ok(ComponentMapping {
            component,
            atoms: atom_map,
            bonds: bond_map,
        })
    }
}

fn copy_stereo(
    source: &Molecule,
    target: &mut Molecule,
    atoms: &BTreeMap<AtomId, AtomId>,
    bonds: &BTreeMap<BondId, BondId>,
) -> Result<(), ModelBuildError> {
    let mut elements = BTreeMap::new();
    for (source_id, element) in source.stereo_elements() {
        let copied = StereoElement {
            kind: remap_stereo_kind(&element.kind, atoms, bonds)?,
            specifiedness: element.specifiedness,
            source: element.source,
            group: None,
            descriptor: None,
        };
        let target_id = target.add_stereo_element(copied)?;
        elements.insert(source_id, target_id);
    }
    for (_, group) in source.stereo_groups() {
        let members = group
            .members
            .iter()
            .map(|member| {
                elements
                    .get(member)
                    .copied()
                    .ok_or(ModelBuildError::InvalidSourceTopology(
                        "stereo group references an unavailable element",
                    ))
            })
            .collect::<Result<Vec<_>, _>>()?;
        target.add_stereo_group(StereoGroup {
            kind: group.kind,
            members,
        })?;
    }
    for mark in source.stereo_bond_marks() {
        target.set_stereo_bond_mark(StereoBondMark {
            bond: mapped_bond(bonds, mark.bond)?,
            kind: mark.kind,
            source: mark.source,
        })?;
    }
    Ok(())
}

fn remap_stereo_kind(
    kind: &StereoElementKind,
    atoms: &BTreeMap<AtomId, AtomId>,
    bonds: &BTreeMap<BondId, BondId>,
) -> Result<StereoElementKind, ModelBuildError> {
    Ok(match kind {
        StereoElementKind::Tetrahedral(stereo) => {
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center: mapped_atom(atoms, stereo.center)?,
                carriers: stereo
                    .carriers
                    .iter()
                    .map(|carrier| remap_carrier(*carrier, atoms))
                    .collect::<Result<_, _>>()?,
                orientation: stereo.orientation,
            })
        }
        StereoElementKind::DoubleBond(stereo) => StereoElementKind::DoubleBond(DoubleBondStereo {
            bond: mapped_bond(bonds, stereo.bond)?,
            left: mapped_atom(atoms, stereo.left)?,
            right: mapped_atom(atoms, stereo.right)?,
            left_carrier: remap_carrier(stereo.left_carrier, atoms)?,
            right_carrier: remap_carrier(stereo.right_carrier, atoms)?,
            orientation: stereo.orientation,
        }),
        StereoElementKind::Axis(stereo) => StereoElementKind::Axis(AxisStereo {
            axis: mapped_bond(bonds, stereo.axis)?,
            carriers: stereo
                .carriers
                .iter()
                .map(|carrier| remap_carrier(*carrier, atoms))
                .collect::<Result<_, _>>()?,
            orientation: stereo.orientation,
        }),
    })
}

fn remap_carrier(
    carrier: StereoCarrier,
    atoms: &BTreeMap<AtomId, AtomId>,
) -> Result<StereoCarrier, ModelBuildError> {
    match carrier {
        StereoCarrier::Atom(atom) => Ok(StereoCarrier::Atom(mapped_atom(atoms, atom)?)),
        StereoCarrier::ImplicitHydrogen => Ok(StereoCarrier::ImplicitHydrogen),
        StereoCarrier::ImplicitLonePair => Ok(StereoCarrier::ImplicitLonePair),
    }
}

fn mapped_atom(
    atoms: &BTreeMap<AtomId, AtomId>,
    source: AtomId,
) -> Result<AtomId, ModelBuildError> {
    atoms
        .get(&source)
        .copied()
        .ok_or(ModelBuildError::InvalidSourceTopology(
            "stereo or bond references an unavailable atom",
        ))
}

fn mapped_bond(
    bonds: &BTreeMap<BondId, BondId>,
    source: BondId,
) -> Result<BondId, ModelBuildError> {
    bonds
        .get(&source)
        .copied()
        .ok_or(ModelBuildError::InvalidSourceTopology(
            "stereo references an unavailable bond",
        ))
}

pub(crate) fn point_is_finite(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelError {
    InvalidComponentId(ComponentId),
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidComponentId(id) => write!(f, "invalid component id: {id}"),
        }
    }
}

impl std::error::Error for ModelError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PositionError {
    InvalidAtomId(AtomId),
    PositionCountMismatch { expected: usize, actual: usize },
    NonFinitePosition { atom: AtomId },
}

impl fmt::Display for PositionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(atom) => write!(f, "invalid model atom id: {atom}"),
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
    EmptyComponent,
    InvalidConformerId(ConformerId),
    MissingPosition { atom: AtomId },
    NonFinitePosition { atom: AtomId },
    InvalidSourceTopology(&'static str),
    Topology(MoleculeError),
}

impl fmt::Display for ModelBuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyModel => write!(f, "a molecular model must contain at least one component"),
            Self::EmptyComponent => write!(f, "a model component must contain at least one atom"),
            Self::InvalidConformerId(id) => write!(f, "invalid source conformer id: {id}"),
            Self::MissingPosition { atom } => {
                write!(f, "source conformer has no position for atom {atom}")
            }
            Self::NonFinitePosition { atom } => {
                write!(f, "source conformer position for atom {atom} is not finite")
            }
            Self::InvalidSourceTopology(message) => {
                write!(f, "invalid source topology: {message}")
            }
            Self::Topology(error) => write!(f, "cannot copy source topology: {error}"),
        }
    }
}

impl std::error::Error for ModelBuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Topology(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MoleculeError> for ModelBuildError {
    fn from(error: MoleculeError) -> Self {
        Self::Topology(error)
    }
}
