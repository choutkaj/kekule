//! Topology-bound coordinate state, models, borrowed views, observations, and
//! finite non-temporal ensembles.

use std::fmt;

use crate::bio::MacroMolecule;
use crate::core::{AtomId, ConformerId, Molecule, PropMap};
use crate::geometry::{PeriodicCell, PeriodicCellError, Point3};
use crate::small::SmallMolecule;
use crate::topology::{
    InstanceAtomId, MoleculeDefinitionId, MoleculeInstanceId, MoleculeInstanceMetadata, Topology,
    TopologyAtomIndex, TopologyBuildError, TopologyBuilder, TopologyIdentity,
};
use crate::units::{Quantity, UnitError, MODEL_LENGTH_UNIT};

/// One complete finite Cartesian array in one topology's dense atom order.
#[derive(Debug, Clone, PartialEq)]
pub struct Positions {
    topology: TopologyIdentity,
    values: Vec<Point3>,
}

impl Positions {
    pub fn new<T>(topology: &Topology, positions: Quantity<T>) -> Result<Self, PositionError>
    where
        T: AsRef<[Point3]>,
    {
        let factor = positions.unit().conversion_factor_to(MODEL_LENGTH_UNIT)?;
        let source = positions.value().as_ref();
        validate_position_count(topology, source.len())?;
        let values = source
            .iter()
            .copied()
            .enumerate()
            .map(|(index, point)| {
                let point = Point3::new(point.x * factor, point.y * factor, point.z * factor);
                if !point.is_finite() {
                    return Err(PositionError::NonFinitePosition {
                        atom: topology.atom_ids()[index],
                    });
                }
                Ok(point)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology: topology.identity(),
            values,
        })
    }

    pub fn zeros(topology: &Topology) -> Self {
        Self {
            topology: topology.identity(),
            values: vec![Point3::origin(); topology.atom_count()],
        }
    }

    pub fn topology_identity(&self) -> &TopologyIdentity {
        &self.topology
    }

    pub fn is_compatible(&self, topology: &Topology) -> bool {
        self.topology == topology.identity()
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    pub fn values(&self) -> Quantity<&[Point3]> {
        Quantity::new(self.values.as_slice(), MODEL_LENGTH_UNIT)
    }

    pub fn position_at(&self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.values
            .get(index.index())
            .copied()
            .map(|point| Quantity::new(point, MODEL_LENGTH_UNIT))
            .ok_or(PositionError::InvalidAtomIndex(index))
    }

    pub fn position(
        &self,
        topology: &Topology,
        atom: InstanceAtomId,
    ) -> Result<Quantity<Point3>, PositionError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(PositionError::InvalidAtomId(atom))?;
        self.position_at(index)
    }

    pub fn set_position(
        &mut self,
        topology: &Topology,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), PositionError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(PositionError::InvalidAtomId(atom))?;
        let point = position.into_unit(MODEL_LENGTH_UNIT)?.into_value();
        if !point.is_finite() {
            return Err(PositionError::NonFinitePosition { atom });
        }
        self.values[index.index()] = point;
        Ok(())
    }

    /// Replaces all positions transactionally while reusing the current
    /// allocation when the capacity permits.
    pub fn set_all<T>(
        &mut self,
        topology: &Topology,
        positions: Quantity<T>,
    ) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        self.ensure_compatible(topology)?;
        let factor = positions.unit().conversion_factor_to(MODEL_LENGTH_UNIT)?;
        let source = positions.value().as_ref();
        validate_position_count(topology, source.len())?;
        for (index, point) in source.iter().copied().enumerate() {
            let converted = Point3::new(point.x * factor, point.y * factor, point.z * factor);
            if !converted.is_finite() {
                return Err(PositionError::NonFinitePosition {
                    atom: topology.atom_ids()[index],
                });
            }
        }
        for (destination, source) in self.values.iter_mut().zip(source.iter().copied()) {
            *destination = Point3::new(source.x * factor, source.y * factor, source.z * factor);
        }
        Ok(())
    }

    pub(crate) fn values_raw(&self) -> &[Point3] {
        &self.values
    }

    fn ensure_compatible(&self, topology: &Topology) -> Result<(), PositionError> {
        if !self.is_compatible(topology) {
            return Err(PositionError::TopologyIdentityMismatch);
        }
        Ok(())
    }
}

fn validate_position_count(topology: &Topology, actual: usize) -> Result<(), PositionError> {
    if actual != topology.atom_count() {
        return Err(PositionError::PositionCountMismatch {
            expected: topology.atom_count(),
            actual,
        });
    }
    if u32::try_from(actual).is_err() && actual != (u32::MAX as usize) + 1 {
        return Err(PositionError::CapacityOverflow);
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PositionError {
    TopologyIdentityMismatch,
    InvalidAtomId(InstanceAtomId),
    InvalidAtomIndex(TopologyAtomIndex),
    PositionCountMismatch { expected: usize, actual: usize },
    NonFinitePosition { atom: InstanceAtomId },
    CapacityOverflow,
    Unit(UnitError),
}

impl fmt::Display for PositionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyIdentityMismatch => {
                formatter.write_str("positions belong to a different topology")
            }
            Self::InvalidAtomId(atom) => write!(formatter, "invalid topology atom: {atom}"),
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid {index}"),
            Self::PositionCountMismatch { expected, actual } => write!(
                formatter,
                "topology requires {expected} positions, but received {actual}"
            ),
            Self::NonFinitePosition { atom } => {
                write!(formatter, "position for atom {atom} is not finite")
            }
            Self::CapacityOverflow => {
                formatter.write_str("position count exceeds fixed-width topology capacity")
            }
            Self::Unit(error) => write!(formatter, "invalid position unit: {error}"),
        }
    }
}

impl std::error::Error for PositionError {}

impl From<UnitError> for PositionError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

/// One complete geometric realization of one exact topology.
#[derive(Debug, Clone, PartialEq)]
pub struct Configuration {
    positions: Positions,
    cell: Option<PeriodicCell>,
}

impl Configuration {
    pub fn new(positions: Positions) -> Self {
        Self {
            positions,
            cell: None,
        }
    }

    pub fn with_cell(positions: Positions, cell: PeriodicCell) -> Self {
        Self {
            positions,
            cell: Some(cell),
        }
    }

    pub fn positions(&self) -> &Positions {
        &self.positions
    }

    pub fn positions_mut(&mut self) -> &mut Positions {
        &mut self.positions
    }

    pub const fn cell(&self) -> Option<&PeriodicCell> {
        self.cell.as_ref()
    }

    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.cell = cell;
    }

    pub fn set_positions(&mut self, positions: Positions) -> Result<(), ConfigurationError> {
        if self.positions.topology != positions.topology {
            return Err(ConfigurationError::TopologyIdentityMismatch);
        }
        self.positions = positions;
        Ok(())
    }

    pub fn view(&self) -> ConfigurationView<'_> {
        ConfigurationView {
            positions: &self.positions,
            cell: self.cell.as_ref(),
        }
    }
}

/// Borrowed configuration state without allocation.
#[derive(Debug, Clone, Copy)]
pub struct ConfigurationView<'a> {
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
}

impl<'a> ConfigurationView<'a> {
    pub fn positions(self) -> &'a Positions {
        self.positions
    }

    pub const fn cell(self) -> Option<&'a PeriodicCell> {
        self.cell
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ConfigurationError {
    TopologyIdentityMismatch,
    InvalidCell(PeriodicCellError),
}

impl fmt::Display for ConfigurationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyIdentityMismatch => {
                formatter.write_str("configuration belongs to a different topology")
            }
            Self::InvalidCell(error) => write!(formatter, "invalid configuration cell: {error}"),
        }
    }
}

impl std::error::Error for ConfigurationError {}

/// Coordinate-model-specific values for one atom observation.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AtomObservation {
    group_pdb: Option<String>,
    source_atom_site_id: Option<String>,
    alternate_location: Option<String>,
    occupancy: Option<f64>,
    occupancy_raw: Option<String>,
    b_factor: Option<f64>,
    b_factor_raw: Option<String>,
    cartesian_raw: [Option<String>; 3],
}

impl AtomObservation {
    pub fn group_pdb(&self) -> Option<&str> {
        self.group_pdb.as_deref()
    }

    pub fn set_group_pdb(&mut self, value: Option<String>) {
        self.group_pdb = value;
    }

    pub fn source_atom_site_id(&self) -> Option<&str> {
        self.source_atom_site_id.as_deref()
    }

    pub fn set_source_atom_site_id(&mut self, value: Option<String>) {
        self.source_atom_site_id = value;
    }

    pub fn alternate_location(&self) -> Option<&str> {
        self.alternate_location.as_deref()
    }

    pub fn set_alternate_location(&mut self, value: Option<String>) {
        self.alternate_location = value;
    }

    pub const fn occupancy(&self) -> Option<f64> {
        self.occupancy
    }

    pub fn set_occupancy(&mut self, value: Option<f64>) -> Result<(), ObservationError> {
        if value.is_some_and(|value| !value.is_finite()) {
            return Err(ObservationError::NonFiniteOccupancy);
        }
        self.occupancy = value;
        Ok(())
    }

    pub fn occupancy_raw(&self) -> Option<&str> {
        self.occupancy_raw.as_deref()
    }

    pub fn set_occupancy_raw(&mut self, value: Option<String>) {
        self.occupancy_raw = value;
    }

    pub const fn b_factor(&self) -> Option<f64> {
        self.b_factor
    }

    pub fn set_b_factor(&mut self, value: Option<f64>) -> Result<(), ObservationError> {
        if value.is_some_and(|value| !value.is_finite()) {
            return Err(ObservationError::NonFiniteBFactor);
        }
        self.b_factor = value;
        Ok(())
    }

    pub fn b_factor_raw(&self) -> Option<&str> {
        self.b_factor_raw.as_deref()
    }

    pub fn set_b_factor_raw(&mut self, value: Option<String>) {
        self.b_factor_raw = value;
    }

    pub fn cartesian_raw(&self) -> &[Option<String>; 3] {
        &self.cartesian_raw
    }

    pub fn set_cartesian_raw(&mut self, value: [Option<String>; 3]) {
        self.cartesian_raw = value;
    }
}

/// Topology-bound observation and provenance state for one configuration.
#[derive(Debug, Clone, PartialEq)]
pub struct StructureObservation {
    topology: TopologyIdentity,
    source_model_id: Option<String>,
    atoms: Vec<AtomObservation>,
    props: PropMap,
}

impl StructureObservation {
    pub fn new(topology: &Topology, atoms: Vec<AtomObservation>) -> Result<Self, ObservationError> {
        if atoms.len() != topology.atom_count() {
            return Err(ObservationError::AtomCountMismatch {
                expected: topology.atom_count(),
                actual: atoms.len(),
            });
        }
        Ok(Self {
            topology: topology.identity(),
            source_model_id: None,
            atoms,
            props: PropMap::new(),
        })
    }

    pub fn empty(topology: &Topology) -> Self {
        Self {
            topology: topology.identity(),
            source_model_id: None,
            atoms: vec![AtomObservation::default(); topology.atom_count()],
            props: PropMap::new(),
        }
    }

    pub fn is_compatible(&self, topology: &Topology) -> bool {
        self.topology == topology.identity()
    }

    pub fn topology_identity(&self) -> &TopologyIdentity {
        &self.topology
    }

    pub fn source_model_id(&self) -> Option<&str> {
        self.source_model_id.as_deref()
    }

    pub fn set_source_model_id(&mut self, value: Option<String>) {
        self.source_model_id = value;
    }

    pub fn atom_at(&self, index: TopologyAtomIndex) -> Option<&AtomObservation> {
        self.atoms.get(index.index())
    }

    pub fn atom(
        &self,
        topology: &Topology,
        atom: InstanceAtomId,
    ) -> Result<&AtomObservation, ObservationError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(ObservationError::InvalidAtomId(atom))?;
        Ok(&self.atoms[index.index()])
    }

    pub fn set_atom(
        &mut self,
        topology: &Topology,
        atom: InstanceAtomId,
        observation: AtomObservation,
    ) -> Result<(), ObservationError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(ObservationError::InvalidAtomId(atom))?;
        self.atoms[index.index()] = observation;
        Ok(())
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    fn ensure_compatible(&self, topology: &Topology) -> Result<(), ObservationError> {
        if !self.is_compatible(topology) {
            return Err(ObservationError::TopologyIdentityMismatch);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ObservationError {
    TopologyIdentityMismatch,
    AtomCountMismatch { expected: usize, actual: usize },
    InvalidAtomId(InstanceAtomId),
    NonFiniteOccupancy,
    NonFiniteBFactor,
}

impl fmt::Display for ObservationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyIdentityMismatch => {
                formatter.write_str("structure observation belongs to a different topology")
            }
            Self::AtomCountMismatch { expected, actual } => write!(
                formatter,
                "structure observation requires {expected} atoms, but received {actual}"
            ),
            Self::InvalidAtomId(atom) => write!(formatter, "invalid observed atom: {atom}"),
            Self::NonFiniteOccupancy => formatter.write_str("structure occupancy must be finite"),
            Self::NonFiniteBFactor => formatter.write_str("structure B-factor must be finite"),
        }
    }
}

impl std::error::Error for ObservationError {}

/// One concrete realization of one immutable topology.
#[derive(Debug, Clone)]
pub struct Model {
    topology: Topology,
    configuration: Configuration,
    observation: Option<StructureObservation>,
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.topology.same_identity(&other.topology)
            && self.configuration == other.configuration
            && self.observation == other.observation
    }
}

impl Model {
    pub fn new(topology: Topology, configuration: Configuration) -> Result<Self, ModelError> {
        if !configuration.positions.is_compatible(&topology) {
            return Err(ModelError::TopologyIdentityMismatch);
        }
        Ok(Self {
            topology,
            configuration,
            observation: None,
        })
    }

    pub fn with_observation(
        topology: Topology,
        configuration: Configuration,
        observation: Option<StructureObservation>,
    ) -> Result<Self, ModelError> {
        let mut model = Self::new(topology, configuration)?;
        model.set_observation(observation)?;
        Ok(model)
    }

    pub fn builder() -> ModelBuilder {
        ModelBuilder::new()
    }

    pub fn from_small_molecule(
        molecule: &SmallMolecule,
        conformer: ConformerId,
    ) -> Result<Self, ModelBuildError> {
        let mut builder = ModelBuilder::new();
        builder.add_small_molecule(molecule, conformer)?;
        builder.build()
    }

    pub fn from_macro_molecule(
        molecule: &MacroMolecule,
        conformer: ConformerId,
    ) -> Result<Self, ModelBuildError> {
        let mut builder = ModelBuilder::new();
        builder.add_macro_molecule(molecule, conformer)?;
        builder.build()
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn configuration(&self) -> &Configuration {
        &self.configuration
    }

    pub fn configuration_mut(&mut self) -> &mut Configuration {
        &mut self.configuration
    }

    pub fn set_configuration(&mut self, configuration: Configuration) -> Result<(), ModelError> {
        if !configuration.positions.is_compatible(&self.topology) {
            return Err(ModelError::TopologyIdentityMismatch);
        }
        self.configuration = configuration;
        Ok(())
    }

    pub fn positions(&self) -> Quantity<&[Point3]> {
        self.configuration.positions.values()
    }

    pub fn position(&self, atom: InstanceAtomId) -> Result<Quantity<Point3>, PositionError> {
        self.configuration.positions.position(&self.topology, atom)
    }

    pub fn position_at(&self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.configuration.positions.position_at(index)
    }

    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), PositionError> {
        self.configuration
            .positions
            .set_position(&self.topology, atom, position)
    }

    pub fn set_positions<T>(&mut self, positions: Quantity<T>) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        self.configuration
            .positions
            .set_all(&self.topology, positions)
    }

    pub const fn cell(&self) -> Option<&PeriodicCell> {
        self.configuration.cell()
    }

    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.configuration.set_cell(cell);
    }

    pub fn observation(&self) -> Option<&StructureObservation> {
        self.observation.as_ref()
    }

    pub fn set_observation(
        &mut self,
        observation: Option<StructureObservation>,
    ) -> Result<(), ModelError> {
        if observation
            .as_ref()
            .is_some_and(|observation| !observation.is_compatible(&self.topology))
        {
            return Err(ModelError::TopologyIdentityMismatch);
        }
        self.observation = observation;
        Ok(())
    }

    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    pub fn view(&self) -> ModelView<'_> {
        ModelView {
            topology: &self.topology,
            configuration: self.configuration.view(),
        }
    }

    /// Copies one instance's current positions to a compatible local conformer.
    pub fn instance_to_conformer(
        &self,
        instance: MoleculeInstanceId,
        target: &mut Molecule,
        conformer: ConformerId,
    ) -> Result<(), InstanceToConformerError> {
        let source = self
            .topology
            .graph_for_instance(instance)
            .map_err(|_| InstanceToConformerError::InvalidMoleculeInstanceId(instance))?;
        let mut updated = target
            .conformer(conformer)
            .map_err(|_| InstanceToConformerError::InvalidConformerId(conformer))?
            .clone();

        for atom in source.atom_ids() {
            if target.atom(atom).is_err() {
                return Err(InstanceToConformerError::MissingTargetAtom(atom));
            }
        }
        for atom in target.atom_ids() {
            if source.atom(atom).is_err() {
                return Err(InstanceToConformerError::UnexpectedTargetAtom(atom));
            }
        }
        for atom in source.atom_ids() {
            let qualified = InstanceAtomId::new(instance, atom);
            updated.set_position(atom, self.position(qualified)?)?;
        }
        *target
            .conformer_mut(conformer)
            .expect("validated conformer remains live") = updated;
        Ok(())
    }

    pub(crate) fn positions_value(&self) -> &[Point3] {
        self.configuration.positions.values_raw()
    }
}

/// Borrowed topology-plus-configuration view for coordinate-dependent kernels.
#[derive(Debug, Clone, Copy)]
pub struct ModelView<'a> {
    topology: &'a Topology,
    configuration: ConfigurationView<'a>,
}

impl<'a> ModelView<'a> {
    pub fn new(
        topology: &'a Topology,
        configuration: ConfigurationView<'a>,
    ) -> Result<Self, ModelError> {
        if !configuration.positions().is_compatible(topology) {
            return Err(ModelError::TopologyIdentityMismatch);
        }
        Ok(Self {
            topology,
            configuration,
        })
    }

    pub const fn topology(self) -> &'a Topology {
        self.topology
    }

    pub const fn configuration(self) -> ConfigurationView<'a> {
        self.configuration
    }

    pub fn positions(self) -> Quantity<&'a [Point3]> {
        self.configuration.positions().values()
    }

    pub fn position(self, atom: InstanceAtomId) -> Result<Quantity<Point3>, PositionError> {
        self.configuration.positions().position(self.topology, atom)
    }

    pub fn position_at(self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.configuration.positions().position_at(index)
    }

    pub const fn cell(self) -> Option<&'a PeriodicCell> {
        self.configuration.cell()
    }

    pub fn atom_count(self) -> usize {
        self.topology.atom_count()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ModelError {
    TopologyIdentityMismatch,
}

impl fmt::Display for ModelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyIdentityMismatch => {
                formatter.write_str("model configuration belongs to a different topology")
            }
        }
    }
}

impl std::error::Error for ModelError {}

/// Convenience builder that assembles topology and one complete configuration.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ModelBuilder {
    topology: TopologyBuilder,
    positions: Vec<Point3>,
}

impl ModelBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn topology_builder(&self) -> &TopologyBuilder {
        &self.topology
    }

    pub fn add_small_molecule_definition(
        &mut self,
        molecule: &SmallMolecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self.topology.add_small_molecule_definition(molecule)?)
    }

    pub fn add_macro_molecule_definition(
        &mut self,
        molecule: &MacroMolecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self.topology.add_macro_molecule_definition(molecule)?)
    }

    pub fn add_instance<T>(
        &mut self,
        definition: MoleculeDefinitionId,
        positions: Quantity<T>,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError>
    where
        T: AsRef<[Point3]>,
    {
        let graph = self.topology.definition(definition)?.graph();
        let staged = stage_instance_positions(graph, positions)?;
        self.positions
            .try_reserve(staged.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let instance = self.topology.add_instance(definition, metadata)?;
        self.positions.extend(staged);
        Ok(instance)
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
        if molecule.graph().atom_count() == 0 {
            return Err(ModelBuildError::Topology(
                TopologyBuildError::EmptyMoleculeDefinition,
            ));
        }
        let staged = stage_conformer_positions(molecule.graph(), conformer)?;
        self.positions
            .try_reserve(staged.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let (_, instance) = self
            .topology
            .add_small_molecule_instance(molecule, metadata)?;
        self.positions.extend(staged);
        Ok(instance)
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
        if molecule.graph().atom_count() == 0 {
            return Err(ModelBuildError::Topology(
                TopologyBuildError::EmptyMoleculeDefinition,
            ));
        }
        molecule.validate().map_err(|error| {
            ModelBuildError::Topology(TopologyBuildError::InvalidMacroMolecule(error))
        })?;
        let staged = stage_conformer_positions(molecule.graph(), conformer)?;
        self.positions
            .try_reserve(staged.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let (_, instance) = self
            .topology
            .add_macro_molecule_instance(molecule, metadata)?;
        self.positions.extend(staged);
        Ok(instance)
    }

    pub fn build(self) -> Result<Model, ModelBuildError> {
        let topology = self.topology.build()?;
        let positions =
            Positions::new(&topology, Quantity::new(self.positions, MODEL_LENGTH_UNIT))?;
        Ok(Model::new(topology, Configuration::new(positions))
            .expect("builder creates exactly topology-bound configuration"))
    }
}

fn stage_conformer_positions(
    graph: &Molecule,
    conformer: ConformerId,
) -> Result<Vec<Point3>, ModelBuildError> {
    let conformer = graph
        .conformer(conformer)
        .map_err(|_| ModelBuildError::InvalidConformerId(conformer))?;
    let mut positions = Vec::new();
    positions
        .try_reserve_exact(graph.atom_count())
        .map_err(|_| ModelBuildError::CapacityOverflow)?;
    for atom in graph.atom_ids() {
        let point = conformer
            .position(atom)
            .ok_or(ModelBuildError::MissingPosition { atom })?
            .into_unit(MODEL_LENGTH_UNIT)?
            .into_value();
        if !point.is_finite() {
            return Err(ModelBuildError::NonFinitePosition { atom });
        }
        positions.push(point);
    }
    Ok(positions)
}

fn stage_instance_positions<T>(
    graph: &Molecule,
    positions: Quantity<T>,
) -> Result<Vec<Point3>, ModelBuildError>
where
    T: AsRef<[Point3]>,
{
    let source = positions.value().as_ref();
    if source.len() != graph.atom_count() {
        return Err(ModelBuildError::InstancePositionCountMismatch {
            expected: graph.atom_count(),
            actual: source.len(),
        });
    }
    let factor = positions.unit().conversion_factor_to(MODEL_LENGTH_UNIT)?;
    source
        .iter()
        .copied()
        .zip(graph.atom_ids())
        .map(|(point, atom)| {
            let point = Point3::new(point.x * factor, point.y * factor, point.z * factor);
            if !point.is_finite() {
                return Err(ModelBuildError::NonFinitePosition { atom });
            }
            Ok(point)
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelBuildError {
    InvalidConformerId(ConformerId),
    MissingPosition { atom: AtomId },
    NonFinitePosition { atom: AtomId },
    InstancePositionCountMismatch { expected: usize, actual: usize },
    CapacityOverflow,
    Topology(TopologyBuildError),
    Position(PositionError),
    Unit(UnitError),
}

impl fmt::Display for ModelBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConformerId(id) => write!(formatter, "invalid source conformer: {id}"),
            Self::MissingPosition { atom } => {
                write!(
                    formatter,
                    "source conformer has no position for atom {atom}"
                )
            }
            Self::NonFinitePosition { atom } => {
                write!(formatter, "source position for atom {atom} is not finite")
            }
            Self::InstancePositionCountMismatch { expected, actual } => write!(
                formatter,
                "definition instance requires {expected} positions, but received {actual}"
            ),
            Self::CapacityOverflow => {
                formatter.write_str("model construction exceeds coordinate capacity")
            }
            Self::Topology(error) => write!(formatter, "cannot build topology: {error}"),
            Self::Position(error) => write!(formatter, "cannot build positions: {error}"),
            Self::Unit(error) => write!(formatter, "invalid source position unit: {error}"),
        }
    }
}

impl std::error::Error for ModelBuildError {}

impl From<TopologyBuildError> for ModelBuildError {
    fn from(error: TopologyBuildError) -> Self {
        Self::Topology(error)
    }
}

impl From<PositionError> for ModelBuildError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<UnitError> for ModelBuildError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum InstanceToConformerError {
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    InvalidConformerId(ConformerId),
    MissingTargetAtom(AtomId),
    UnexpectedTargetAtom(AtomId),
    Position(PositionError),
    Unit(UnitError),
}

impl fmt::Display for InstanceToConformerError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMoleculeInstanceId(id) => {
                write!(formatter, "invalid molecule instance: {id}")
            }
            Self::InvalidConformerId(id) => write!(formatter, "invalid target conformer: {id}"),
            Self::MissingTargetAtom(atom) => {
                write!(formatter, "target molecule is missing source atom {atom}")
            }
            Self::UnexpectedTargetAtom(atom) => {
                write!(formatter, "target molecule contains unexpected atom {atom}")
            }
            Self::Position(error) => write!(formatter, "cannot read model position: {error}"),
            Self::Unit(error) => write!(formatter, "invalid target conformer unit: {error}"),
        }
    }
}

impl std::error::Error for InstanceToConformerError {}

impl From<PositionError> for InstanceToConformerError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<UnitError> for InstanceToConformerError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

/// One finite non-temporal ensemble member.
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleMember {
    configuration: Configuration,
    weight: Option<f64>,
    observation: Option<StructureObservation>,
    props: PropMap,
}

impl EnsembleMember {
    pub fn new(configuration: Configuration) -> Self {
        Self {
            configuration,
            weight: None,
            observation: None,
            props: PropMap::new(),
        }
    }

    pub fn configuration(&self) -> &Configuration {
        &self.configuration
    }

    pub const fn weight(&self) -> Option<f64> {
        self.weight
    }

    pub fn set_weight(&mut self, weight: Option<f64>) -> Result<(), EnsembleError> {
        if weight.is_some_and(|weight| !weight.is_finite() || weight < 0.0) {
            return Err(EnsembleError::InvalidWeight);
        }
        self.weight = weight;
        Ok(())
    }

    pub fn observation(&self) -> Option<&StructureObservation> {
        self.observation.as_ref()
    }

    pub fn set_observation(
        &mut self,
        observation: Option<StructureObservation>,
    ) -> Result<(), EnsembleError> {
        if observation.as_ref().is_some_and(|observation| {
            observation.topology != self.configuration.positions.topology
        }) {
            return Err(EnsembleError::TopologyIdentityMismatch);
        }
        self.observation = observation;
        Ok(())
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    pub fn view<'a>(&'a self, topology: &'a Topology) -> Result<ModelView<'a>, EnsembleError> {
        ModelView::new(topology, self.configuration.view())
            .map_err(|_| EnsembleError::TopologyIdentityMismatch)
    }
}

/// A finite stable-order collection of non-temporal configurations.
#[derive(Debug, Clone)]
pub struct Ensemble {
    topology: Topology,
    members: Vec<EnsembleMember>,
}

impl Ensemble {
    pub fn new(topology: Topology) -> Self {
        Self {
            topology,
            members: Vec::new(),
        }
    }

    pub fn from_members(
        topology: Topology,
        members: impl IntoIterator<Item = EnsembleMember>,
    ) -> Result<Self, EnsembleError> {
        let mut ensemble = Self::new(topology);
        for member in members {
            ensemble.push(member)?;
        }
        Ok(ensemble)
    }

    pub fn from_models(models: &[Model]) -> Result<Self, EnsembleError> {
        let first = models.first().ok_or(EnsembleError::EmptySource)?;
        let mut ensemble = Self::new(first.topology.clone());
        for model in models {
            if !first.topology.same_identity(&model.topology) {
                return Err(EnsembleError::TopologyIdentityMismatch);
            }
            let mut member = EnsembleMember::new(model.configuration.clone());
            member.observation = model.observation.clone();
            ensemble.members.push(member);
        }
        Ok(ensemble)
    }

    pub fn from_small_molecule_conformers(
        molecule: &SmallMolecule,
        conformers: impl IntoIterator<Item = ConformerId>,
    ) -> Result<Self, EnsembleError> {
        let mut topology_builder = TopologyBuilder::new();
        let definition = topology_builder.add_small_molecule_definition(molecule)?;
        topology_builder.add_instance(definition, MoleculeInstanceMetadata::default())?;
        let topology = topology_builder.build()?;
        let mut ensemble = Self::new(topology.clone());
        for conformer in conformers {
            let positions = stage_conformer_positions(molecule.graph(), conformer)
                .map_err(|error| EnsembleError::ModelBuild(Box::new(error)))?;
            let positions = Positions::new(&topology, Quantity::new(positions, MODEL_LENGTH_UNIT))?;
            ensemble.push(EnsembleMember::new(Configuration::new(positions)))?;
        }
        Ok(ensemble)
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn members(&self) -> impl ExactSizeIterator<Item = &EnsembleMember> {
        self.members.iter()
    }

    pub fn member(&self, index: usize) -> Option<&EnsembleMember> {
        self.members.get(index)
    }

    pub fn member_mut(&mut self, index: usize) -> Option<&mut EnsembleMember> {
        self.members.get_mut(index)
    }

    pub fn len(&self) -> usize {
        self.members.len()
    }

    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    pub fn push(&mut self, member: EnsembleMember) -> Result<(), EnsembleError> {
        if !member.configuration.positions.is_compatible(&self.topology)
            || member
                .observation
                .as_ref()
                .is_some_and(|observation| !observation.is_compatible(&self.topology))
        {
            return Err(EnsembleError::TopologyIdentityMismatch);
        }
        self.members.push(member);
        Ok(())
    }

    pub fn views(&self) -> impl ExactSizeIterator<Item = ModelView<'_>> {
        self.members.iter().map(|member| {
            member
                .view(&self.topology)
                .expect("ensemble validates exact topology identity")
        })
    }

    pub fn normalize_weights(&mut self) -> Result<(), EnsembleError> {
        if self.members.is_empty() {
            return Err(EnsembleError::EmptySource);
        }
        let mut total = 0.0;
        for (index, member) in self.members.iter().enumerate() {
            let weight = member
                .weight
                .ok_or(EnsembleError::MissingWeight { member: index })?;
            total += weight;
        }
        if !total.is_finite() || total <= 0.0 {
            return Err(EnsembleError::ZeroTotalWeight);
        }
        for member in &mut self.members {
            member.weight = member.weight.map(|weight| weight / total);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum EnsembleError {
    EmptySource,
    TopologyIdentityMismatch,
    InvalidWeight,
    MissingWeight { member: usize },
    ZeroTotalWeight,
    TopologyBuild(TopologyBuildError),
    ModelBuild(Box<ModelBuildError>),
    Position(PositionError),
}

impl fmt::Display for EnsembleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySource => formatter.write_str("ensemble source is empty"),
            Self::TopologyIdentityMismatch => {
                formatter.write_str("ensemble member belongs to a different topology")
            }
            Self::InvalidWeight => {
                formatter.write_str("ensemble weight must be finite and non-negative")
            }
            Self::MissingWeight { member } => {
                write!(formatter, "ensemble member {member} has no weight")
            }
            Self::ZeroTotalWeight => {
                formatter.write_str("ensemble weights must have a positive finite total")
            }
            Self::TopologyBuild(error) => write!(formatter, "cannot build topology: {error}"),
            Self::ModelBuild(error) => write!(formatter, "cannot build ensemble member: {error}"),
            Self::Position(error) => write!(formatter, "cannot build member positions: {error}"),
        }
    }
}

impl std::error::Error for EnsembleError {}

impl From<TopologyBuildError> for EnsembleError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}

impl From<PositionError> for EnsembleError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Atom, Element, Molecule};
    use crate::geometry::Vector3;
    use crate::units::{ANGSTROM, NANOMETER};

    fn one_atom_topology() -> Topology {
        let mut graph = Molecule::new();
        graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let molecule = SmallMolecule::from_graph(graph);
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&molecule).unwrap();
        builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        builder.build().unwrap()
    }

    fn configuration(topology: &Topology, x: f64) -> Configuration {
        Configuration::new(
            Positions::new(
                topology,
                Quantity::new(vec![Point3::new(x, 0.0, 0.0)], ANGSTROM),
            )
            .unwrap(),
        )
    }

    #[test]
    fn positions_convert_units_reuse_storage_and_update_transactionally() {
        let topology = one_atom_topology();
        let mut positions = Positions::new(
            &topology,
            Quantity::new(vec![Point3::new(0.1, 0.0, 0.0)], NANOMETER),
        )
        .unwrap();
        assert_eq!(positions.values().value()[0], Point3::new(1.0, 0.0, 0.0));
        let pointer = positions.values().value().as_ptr();
        positions
            .set_all(
                &topology,
                Quantity::new(vec![Point3::new(2.0, 0.0, 0.0)], ANGSTROM),
            )
            .unwrap();
        assert_eq!(positions.values().value().as_ptr(), pointer);

        let before = positions.clone();
        assert!(matches!(
            positions.set_all(
                &topology,
                Quantity::new(vec![Point3::new(f64::NAN, 0.0, 0.0)], ANGSTROM)
            ),
            Err(PositionError::NonFinitePosition { .. })
        ));
        assert_eq!(positions, before);
    }

    #[test]
    fn model_requires_exact_topology_and_views_do_not_copy_coordinates() {
        let topology = one_atom_topology();
        let independent = one_atom_topology();
        assert!(topology.structurally_equivalent(&independent));
        assert!(!topology.same_identity(&independent));

        let wrong_configuration = configuration(&independent, 1.0);
        assert_eq!(
            Model::new(topology.clone(), wrong_configuration),
            Err(ModelError::TopologyIdentityMismatch)
        );

        let mut model = Model::new(topology.clone(), configuration(&topology, 1.0)).unwrap();
        let view = model.view();
        assert_eq!(
            view.positions().value().as_ptr(),
            model.positions().value().as_ptr()
        );
        let clone = model.clone();
        assert!(model.topology().same_identity(clone.topology()));
        let identity = model.topology().identity();
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 11.0, 12.0), ANGSTROM),
            [true; 3],
        )
        .unwrap();
        model.set_cell(Some(cell));
        assert_eq!(model.topology().identity(), identity);
        assert_eq!(model.cell(), Some(&cell));
    }

    #[test]
    fn ensembles_validate_identity_observations_order_and_weights() {
        let topology = one_atom_topology();
        let mut first = EnsembleMember::new(configuration(&topology, 1.0));
        first.set_weight(Some(1.0)).unwrap();
        let mut first_observation = StructureObservation::empty(&topology);
        first_observation.set_source_model_id(Some("first".to_owned()));
        first.set_observation(Some(first_observation)).unwrap();

        let mut second = EnsembleMember::new(configuration(&topology, 2.0));
        second.set_weight(Some(3.0)).unwrap();
        let mut invalid_weight = EnsembleMember::new(configuration(&topology, 3.0));
        assert_eq!(
            invalid_weight.set_weight(Some(-1.0)),
            Err(EnsembleError::InvalidWeight)
        );
        let mut ensemble = Ensemble::from_members(topology.clone(), [first, second]).unwrap();
        ensemble.normalize_weights().unwrap();
        assert_eq!(ensemble.member(0).unwrap().weight(), Some(0.25));
        assert_eq!(ensemble.member(1).unwrap().weight(), Some(0.75));
        assert_eq!(
            ensemble
                .views()
                .map(|view| view.positions().value()[0].x)
                .collect::<Vec<_>>(),
            vec![1.0, 2.0]
        );
        assert_eq!(
            ensemble
                .member(0)
                .unwrap()
                .observation()
                .unwrap()
                .source_model_id(),
            Some("first")
        );

        let independent = one_atom_topology();
        assert_eq!(
            ensemble.push(EnsembleMember::new(configuration(&independent, 3.0))),
            Err(EnsembleError::TopologyIdentityMismatch)
        );
    }
}
