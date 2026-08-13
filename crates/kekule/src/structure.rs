//! Topology-bound coordinate state, models, borrowed views, atom data, and
//! finite non-temporal ensembles.

use std::fmt;
use std::sync::Arc;

use crate::bio::MacroMolecule;
use crate::core::{AtomId, ConformerError, ConformerId, Molecule, PropMap};
use crate::geometry::{PeriodicCell, Point3};
use crate::small::SmallMolecule;
use crate::topology::{
    InstanceAtomId, MoleculeDefinitionId, MoleculeInstanceId, MoleculeInstanceMetadata, Topology,
    TopologyAtomIndex, TopologyBuildError, TopologyBuilder, TopologyMapping,
};
use crate::units::{Quantity, UnitError, MODEL_LENGTH_UNIT, SQUARE_ANGSTROM};

/// One complete finite Cartesian array in one topology's dense atom order.
#[derive(Debug, Clone)]
pub struct Positions {
    topology: Arc<Topology>,
    values: Vec<Point3>,
}

impl PartialEq for Positions {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.values == other.values
    }
}

impl Positions {
    pub fn new<T>(topology: &Arc<Topology>, positions: Quantity<T>) -> Result<Self, PositionError>
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
            topology: Arc::clone(topology),
            values,
        })
    }

    pub fn zeros(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            values: vec![Point3::origin(); topology.atom_count()],
        }
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn is_compatible(&self, topology: &Arc<Topology>) -> bool {
        Arc::ptr_eq(&self.topology, topology)
    }

    pub(crate) fn topology_arc(&self) -> &Arc<Topology> {
        &self.topology
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

    /// Copies complete positions through a checked topology lineage mapping.
    pub fn remap_to(
        &self,
        source: &Arc<Topology>,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TopologyRemapError> {
        if !self.is_compatible(source) {
            return Err(TopologyRemapError::SourceTopologyMismatch);
        }
        let values = remap::dense_atom_values(&self.values, source, target, mapping)?;
        Ok(Self {
            topology: Arc::clone(target),
            values,
        })
    }

    /// Copies complete positions through a checked topology lineage mapping
    /// while retaining this array's allocation.
    ///
    /// Validation completes before any destination position changes.
    pub fn copy_remapped_from(
        &mut self,
        source: &Self,
        source_topology: &Arc<Topology>,
        target_topology: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<(), TopologyRemapError> {
        if !source.is_compatible(source_topology) {
            return Err(TopologyRemapError::SourceTopologyMismatch);
        }
        if !self.is_compatible(target_topology) {
            return Err(TopologyRemapError::TargetTopologyMismatch);
        }
        if source.values.len() != source_topology.atom_count() {
            return Err(TopologyRemapError::SourceAtomCountMismatch {
                expected: source_topology.atom_count(),
                actual: source.values.len(),
            });
        }
        remap::validate_complete_atom_mapping(source_topology, target_topology, mapping)?;
        for (source_index, target_index) in mapping.atom_index_pairs() {
            self.values[target_index.index()] = source.values[source_index.index()];
        }
        Ok(())
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
        topology: &Arc<Topology>,
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
        topology: &Arc<Topology>,
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
        topology: &Arc<Topology>,
        positions: Quantity<T>,
    ) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        let factor = self.validate_replacement(topology, &positions)?;
        self.copy_from_validated(positions.value().as_ref(), factor);
        Ok(())
    }

    /// Validates a complete replacement without changing this array.
    ///
    /// This supports external topology-bound containers that must validate
    /// several coupled fields before publishing any of them transactionally.
    pub fn validate_all<T>(
        &self,
        topology: &Arc<Topology>,
        positions: &Quantity<T>,
    ) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        self.validate_replacement(topology, positions).map(drop)
    }

    pub(crate) fn validate_replacement<T>(
        &self,
        topology: &Arc<Topology>,
        positions: &Quantity<T>,
    ) -> Result<f64, PositionError>
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
        Ok(factor)
    }

    pub(crate) fn copy_from_validated(&mut self, source: &[Point3], factor: f64) {
        for (destination, source) in self.values.iter_mut().zip(source.iter().copied()) {
            *destination = Point3::new(source.x * factor, source.y * factor, source.z * factor);
        }
    }

    fn ensure_compatible(&self, topology: &Arc<Topology>) -> Result<(), PositionError> {
        if !self.is_compatible(topology) {
            return Err(PositionError::TopologyMismatch);
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
    crate::core::checked_fixed_id_collection_len(0, actual)
        .map_err(|_| PositionError::CapacityOverflow)?;
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PositionError {
    TopologyMismatch,
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
            Self::TopologyMismatch => {
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

/// Topology-bound model-level per-atom data in dense atom order.
///
/// Each scientific field is an optional column. A field that is absent for all
/// atoms therefore requires no per-atom allocation.
#[derive(Debug, Clone)]
pub struct AtomData {
    topology: Arc<Topology>,
    occupancies: Option<Vec<Option<f64>>>,
    b_factors: Option<Vec<Option<f64>>>,
}

impl PartialEq for AtomData {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.occupancies == other.occupancies
            && self.b_factors == other.b_factors
    }
}

impl AtomData {
    /// Creates atom data with no allocated scientific columns.
    pub fn new(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            occupancies: None,
            b_factors: None,
        }
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn is_compatible(&self, topology: &Arc<Topology>) -> bool {
        Arc::ptr_eq(&self.topology, topology)
    }

    pub(crate) fn topology_arc(&self) -> &Arc<Topology> {
        &self.topology
    }

    /// Returns the atom count of the bound topology.
    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    /// Returns whether every supported scientific column is wholly absent.
    pub fn is_empty(&self) -> bool {
        self.occupancies.is_none() && self.b_factors.is_none()
    }

    pub fn occupancies(&self) -> Option<&[Option<f64>]> {
        self.occupancies.as_deref()
    }

    /// Returns the dense B-factor column in canonical square angstroms.
    pub fn b_factors(&self) -> Option<Quantity<&[Option<f64>]>> {
        self.b_factors
            .as_deref()
            .map(|values| Quantity::new(values, SQUARE_ANGSTROM))
    }

    pub fn occupancy(
        &self,
        topology: &Arc<Topology>,
        atom: InstanceAtomId,
    ) -> Result<Option<f64>, AtomDataError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(AtomDataError::InvalidAtomId(atom))?;
        self.occupancy_at(index)
    }

    pub fn occupancy_at(&self, index: TopologyAtomIndex) -> Result<Option<f64>, AtomDataError> {
        value_at(&self.occupancies, self.atom_count(), index)
    }

    pub fn b_factor(
        &self,
        topology: &Arc<Topology>,
        atom: InstanceAtomId,
    ) -> Result<Option<Quantity<f64>>, AtomDataError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(AtomDataError::InvalidAtomId(atom))?;
        self.b_factor_at(index)
    }

    pub fn b_factor_at(
        &self,
        index: TopologyAtomIndex,
    ) -> Result<Option<Quantity<f64>>, AtomDataError> {
        value_at(&self.b_factors, self.atom_count(), index)
            .map(|value| value.map(|value| Quantity::new(value, SQUARE_ANGSTROM)))
    }

    pub fn set_occupancy(
        &mut self,
        topology: &Arc<Topology>,
        atom: InstanceAtomId,
        value: Option<f64>,
    ) -> Result<(), AtomDataError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(AtomDataError::InvalidAtomId(atom))?;
        self.set_occupancy_at(index, value)
    }

    pub fn set_occupancy_at(
        &mut self,
        index: TopologyAtomIndex,
        value: Option<f64>,
    ) -> Result<(), AtomDataError> {
        validate_index(self.atom_count(), index)?;
        validate_value(value, index, AtomDataField::Occupancy)?;
        let len = self.atom_count();
        set_column_value(&mut self.occupancies, len, index, value);
        Ok(())
    }

    pub fn set_b_factor(
        &mut self,
        topology: &Arc<Topology>,
        atom: InstanceAtomId,
        value: Option<Quantity<f64>>,
    ) -> Result<(), AtomDataError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(AtomDataError::InvalidAtomId(atom))?;
        self.set_b_factor_at(index, value)
    }

    pub fn set_b_factor_at(
        &mut self,
        index: TopologyAtomIndex,
        value: Option<Quantity<f64>>,
    ) -> Result<(), AtomDataError> {
        validate_index(self.atom_count(), index)?;
        let value = value
            .map(|value| value.into_unit(SQUARE_ANGSTROM))
            .transpose()?
            .map(|value| value.into_value());
        validate_value(value, index, AtomDataField::BFactor)?;
        let len = self.atom_count();
        set_column_value(&mut self.b_factors, len, index, value);
        Ok(())
    }

    /// Replaces the complete occupancy column transactionally. An all-absent
    /// column is normalized to no allocation.
    pub fn set_occupancies<T>(&mut self, values: T) -> Result<(), AtomDataError>
    where
        T: AsRef<[Option<f64>]>,
    {
        let values = validate_column(
            &self.topology,
            values.as_ref().to_vec(),
            AtomDataField::Occupancy,
        )?;
        self.occupancies = values;
        Ok(())
    }

    /// Clears the complete occupancy column.
    pub fn clear_occupancies(&mut self) {
        self.occupancies = None;
    }

    /// Replaces the complete B-factor column transactionally. Values are
    /// converted to canonical square angstroms, and an all-absent column is
    /// normalized to no allocation.
    pub fn set_b_factors<T>(&mut self, values: Quantity<T>) -> Result<(), AtomDataError>
    where
        T: AsRef<[Option<f64>]>,
    {
        let factor = values.unit().conversion_factor_to(SQUARE_ANGSTROM)?;
        let values = values
            .value()
            .as_ref()
            .iter()
            .copied()
            .map(|value| value.map(|value| value * factor))
            .collect();
        let values = validate_column(&self.topology, values, AtomDataField::BFactor)?;
        self.b_factors = values;
        Ok(())
    }

    /// Clears the complete B-factor column.
    pub fn clear_b_factors(&mut self) {
        self.b_factors = None;
    }

    /// Remaps every present scientific column through checked topology lineage.
    pub fn remap_to(
        &self,
        source: &Arc<Topology>,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TopologyRemapError> {
        if !self.is_compatible(source) {
            return Err(TopologyRemapError::SourceTopologyMismatch);
        }
        let occupancies = self
            .occupancies
            .as_ref()
            .map(|values| remap::dense_atom_values(values, source, target, mapping))
            .transpose()?;
        let b_factors = self
            .b_factors
            .as_ref()
            .map(|values| remap::dense_atom_values(values, source, target, mapping))
            .transpose()?;
        Ok(Self {
            topology: Arc::clone(target),
            occupancies,
            b_factors,
        })
    }

    fn ensure_compatible(&self, topology: &Arc<Topology>) -> Result<(), AtomDataError> {
        if !self.is_compatible(topology) {
            return Err(AtomDataError::TopologyMismatch);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
enum AtomDataField {
    Occupancy,
    BFactor,
}

fn validate_column(
    topology: &Topology,
    values: Vec<Option<f64>>,
    field: AtomDataField,
) -> Result<Option<Vec<Option<f64>>>, AtomDataError> {
    if values.len() != topology.atom_count() {
        return Err(AtomDataError::AtomCountMismatch {
            expected: topology.atom_count(),
            actual: values.len(),
        });
    }
    for (raw_index, value) in values.iter().copied().enumerate() {
        let index = TopologyAtomIndex::new(raw_index as u32);
        validate_value(value, index, field)?;
    }
    if values.iter().all(Option::is_none) {
        Ok(None)
    } else {
        Ok(Some(values))
    }
}

fn validate_value(
    value: Option<f64>,
    index: TopologyAtomIndex,
    field: AtomDataField,
) -> Result<(), AtomDataError> {
    if value.is_some_and(|value| !value.is_finite()) {
        return Err(match field {
            AtomDataField::Occupancy => AtomDataError::NonFiniteOccupancy { index },
            AtomDataField::BFactor => AtomDataError::NonFiniteBFactor { index },
        });
    }
    Ok(())
}

fn validate_index(len: usize, index: TopologyAtomIndex) -> Result<(), AtomDataError> {
    if index.index() >= len {
        return Err(AtomDataError::InvalidAtomIndex(index));
    }
    Ok(())
}

fn value_at(
    column: &Option<Vec<Option<f64>>>,
    len: usize,
    index: TopologyAtomIndex,
) -> Result<Option<f64>, AtomDataError> {
    validate_index(len, index)?;
    Ok(column.as_ref().and_then(|values| values[index.index()]))
}

fn set_column_value(
    column: &mut Option<Vec<Option<f64>>>,
    len: usize,
    index: TopologyAtomIndex,
    value: Option<f64>,
) {
    if value.is_some() && column.is_none() {
        *column = Some(vec![None; len]);
    }
    let Some(values) = column else {
        return;
    };
    values[index.index()] = value;
    if values.iter().all(Option::is_none) {
        *column = None;
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum AtomDataError {
    TopologyMismatch,
    AtomCountMismatch { expected: usize, actual: usize },
    InvalidAtomId(InstanceAtomId),
    InvalidAtomIndex(TopologyAtomIndex),
    NonFiniteOccupancy { index: TopologyAtomIndex },
    NonFiniteBFactor { index: TopologyAtomIndex },
    Unit(UnitError),
}

impl fmt::Display for AtomDataError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("atom data belongs to a different topology")
            }
            Self::AtomCountMismatch { expected, actual } => write!(
                formatter,
                "atom data requires {expected} values per present column, but received {actual}"
            ),
            Self::InvalidAtomId(atom) => write!(formatter, "invalid topology atom: {atom}"),
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid {index}"),
            Self::NonFiniteOccupancy { index } => {
                write!(formatter, "occupancy at {index} must be finite")
            }
            Self::NonFiniteBFactor { index } => {
                write!(formatter, "B-factor at {index} must be finite")
            }
            Self::Unit(error) => write!(formatter, "invalid B-factor unit: {error}"),
        }
    }
}

impl std::error::Error for AtomDataError {}

impl From<UnitError> for AtomDataError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

/// One concrete realization of one immutable topology.
#[derive(Debug, Clone)]
pub struct Model {
    topology: Arc<Topology>,
    positions: Positions,
    cell: Option<PeriodicCell>,
    atom_data: AtomData,
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.positions == other.positions
            && self.cell == other.cell
            && self.atom_data == other.atom_data
    }
}

impl Model {
    /// Creates a non-periodic model without per-atom scientific data.
    pub fn new(topology: Arc<Topology>, positions: Positions) -> Result<Self, ModelError> {
        if !positions.is_compatible(&topology) {
            return Err(ModelError::TopologyMismatch);
        }
        let atom_data = AtomData::new(&topology);
        Ok(Self {
            topology,
            positions,
            cell: None,
            atom_data,
        })
    }

    /// Creates a model from complete topology-bound state.
    pub fn with_atom_data(
        topology: Arc<Topology>,
        positions: Positions,
        cell: Option<PeriodicCell>,
        atom_data: AtomData,
    ) -> Result<Self, ModelError> {
        if !positions.is_compatible(&topology) || !atom_data.is_compatible(&topology) {
            return Err(ModelError::TopologyMismatch);
        }
        Ok(Self {
            topology,
            positions,
            cell,
            atom_data,
        })
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

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub const fn positions(&self) -> &Positions {
        &self.positions
    }

    pub fn position(&self, atom: InstanceAtomId) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position(&self.topology, atom)
    }

    pub fn position_at(&self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position_at(index)
    }

    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), PositionError> {
        self.positions.set_position(&self.topology, atom, position)
    }

    pub fn set_positions<T>(&mut self, positions: Quantity<T>) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        self.positions.set_all(&self.topology, positions)
    }

    pub const fn cell(&self) -> Option<&PeriodicCell> {
        self.cell.as_ref()
    }

    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.cell = cell;
    }

    pub fn atom_data(&self) -> &AtomData {
        &self.atom_data
    }

    pub fn atom_data_mut(&mut self) -> &mut AtomData {
        &mut self.atom_data
    }

    pub fn set_atom_data(&mut self, atom_data: AtomData) -> Result<(), ModelError> {
        if !atom_data.is_compatible(&self.topology) {
            return Err(ModelError::TopologyMismatch);
        }
        self.atom_data = atom_data;
        Ok(())
    }

    pub fn occupancy(&self, atom: InstanceAtomId) -> Result<Option<f64>, AtomDataError> {
        self.atom_data.occupancy(&self.topology, atom)
    }

    pub fn b_factor(&self, atom: InstanceAtomId) -> Result<Option<Quantity<f64>>, AtomDataError> {
        self.atom_data.b_factor(&self.topology, atom)
    }

    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    pub fn view(&self) -> ModelView<'_> {
        ModelView {
            topology: &self.topology,
            positions: &self.positions,
            cell: self.cell.as_ref(),
            atom_data: &self.atom_data,
        }
    }

    /// Remaps this model to an explicitly related target topology.
    ///
    /// The source model is unchanged. Positions, cell, and atom data are staged
    /// before the target model is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use kekule::core::{Atom, Element, Molecule};
    /// use kekule::geometry::Point3;
    /// use kekule::small::SmallMolecule;
    /// use kekule::structure::{Model, Positions};
    /// use kekule::topology::{
    ///     transform, MoleculeInstanceMetadata, MoleculeRole, TopologyBuilder,
    /// };
    /// use kekule::units::{Quantity, ANGSTROM};
    /// use std::sync::Arc;
    ///
    /// let mut ligand_builder = Molecule::builder();
    /// ligand_builder.add_atom(Atom::new(Element::from_symbol("C").unwrap()))?;
    /// let ligand = SmallMolecule::from_graph(ligand_builder.build()?);
    /// let mut water_builder = Molecule::builder();
    /// water_builder.add_atom(Atom::new(Element::from_symbol("O").unwrap()))?;
    /// let water = SmallMolecule::from_graph(water_builder.build()?);
    /// let mut builder = TopologyBuilder::new();
    /// let ligand_definition = builder.add_small_molecule_definition(&ligand)?;
    /// let water_definition = builder.add_small_molecule_definition(&water)?;
    /// builder.add_instance(ligand_definition, MoleculeInstanceMetadata::default())?;
    /// let mut solvent = MoleculeInstanceMetadata::default();
    /// solvent.insert_role(MoleculeRole::Solvent);
    /// let water_instance = builder.add_instance(water_definition, solvent)?;
    /// let topology = Arc::new(builder.build()?);
    /// let positions = Positions::new(
    ///     &topology,
    ///     Quantity::new(
    ///         vec![Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 0.0, 0.0)],
    ///         ANGSTROM,
    ///     ),
    /// )?;
    /// let model = Model::new(Arc::clone(&topology), positions)?;
    ///
    /// let edit = transform::remove_instances(&topology, [water_instance])?;
    /// let target = edit.shared_topology();
    /// let stripped = model.remap_to(&target, edit.mapping())?;
    /// assert_eq!(stripped.atom_count(), 1);
    /// assert_eq!(stripped.positions().values().value()[0].x, 0.0);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn remap_to(
        &self,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TopologyRemapError> {
        let positions = self.positions.remap_to(&self.topology, target, mapping)?;
        let atom_data = self.atom_data.remap_to(&self.topology, target, mapping)?;
        Ok(Self {
            topology: Arc::clone(target),
            positions,
            cell: self.cell,
            atom_data,
        })
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
}

/// Borrowed topology, positions, cell, and atom data for structural kernels.
#[derive(Debug, Clone, Copy)]
pub struct ModelView<'a> {
    topology: &'a Arc<Topology>,
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
    atom_data: &'a AtomData,
}

impl<'a> ModelView<'a> {
    pub fn new(
        topology: &'a Arc<Topology>,
        positions: &'a Positions,
        cell: Option<&'a PeriodicCell>,
        atom_data: &'a AtomData,
    ) -> Result<Self, ModelError> {
        if !positions.is_compatible(topology) || !atom_data.is_compatible(topology) {
            return Err(ModelError::TopologyMismatch);
        }
        Ok(Self {
            topology,
            positions,
            cell,
            atom_data,
        })
    }

    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub(crate) const fn topology_arc(self) -> &'a Arc<Topology> {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    pub const fn positions(self) -> &'a Positions {
        self.positions
    }

    pub fn position(self, atom: InstanceAtomId) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position(self.topology, atom)
    }

    pub fn position_at(self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position_at(index)
    }

    pub const fn cell(self) -> Option<&'a PeriodicCell> {
        self.cell
    }

    pub const fn atom_data(self) -> &'a AtomData {
        self.atom_data
    }

    pub fn occupancy(self, atom: InstanceAtomId) -> Result<Option<f64>, AtomDataError> {
        self.atom_data.occupancy(self.topology, atom)
    }

    pub fn b_factor(self, atom: InstanceAtomId) -> Result<Option<Quantity<f64>>, AtomDataError> {
        self.atom_data.b_factor(self.topology, atom)
    }

    pub fn atom_count(self) -> usize {
        self.topology.atom_count()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ModelError {
    TopologyMismatch,
}

impl fmt::Display for ModelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("model state belongs to a different topology")
            }
        }
    }
}

impl std::error::Error for ModelError {}

/// Convenience builder that assembles topology and one complete model.
///
/// Macro-molecule insertion validates only the explicitly selected conformer
/// while staging positions; topology insertion separately validates static
/// graph/hierarchy consistency and ignores every unselected conformer.
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

    pub(crate) fn add_small_molecule_definition_unchecked_connectedness(
        &mut self,
        molecule: &SmallMolecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self
            .topology
            .add_small_molecule_definition_unchecked_connectedness(molecule)?)
    }

    pub(crate) fn add_macro_molecule_definition_unchecked_connectedness(
        &mut self,
        molecule: &MacroMolecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self
            .topology
            .add_macro_molecule_definition_unchecked_connectedness(molecule)?)
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

    pub(crate) fn add_small_molecule_with_metadata_unchecked_connectedness(
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
            .add_small_molecule_instance_unchecked_connectedness(molecule, metadata)?;
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

    pub(crate) fn add_macro_molecule_with_metadata_unchecked_connectedness(
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
        let staged = stage_conformer_positions(molecule.graph(), conformer)?;
        self.positions
            .try_reserve(staged.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let (_, instance) = self
            .topology
            .add_macro_molecule_instance_unchecked_connectedness(molecule, metadata)?;
        self.positions.extend(staged);
        Ok(instance)
    }

    pub fn build(self) -> Result<Model, ModelBuildError> {
        let topology = Arc::new(self.topology.build()?);
        let positions =
            Positions::new(&topology, Quantity::new(self.positions, MODEL_LENGTH_UNIT))?;
        Ok(Model::new(topology, positions).expect("builder creates exactly topology-bound state"))
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
    Conformer(ConformerError),
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
            Self::Conformer(error) => write!(formatter, "cannot update target conformer: {error}"),
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

impl From<ConformerError> for InstanceToConformerError {
    fn from(error: ConformerError) -> Self {
        Self::Conformer(error)
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
    positions: Positions,
    cell: Option<PeriodicCell>,
    atom_data: AtomData,
    weight: Option<f64>,
    props: PropMap,
}

impl EnsembleMember {
    pub fn new(positions: Positions) -> Self {
        let atom_data = AtomData::new(positions.topology_arc());
        Self {
            positions,
            cell: None,
            atom_data,
            weight: None,
            props: PropMap::new(),
        }
    }

    pub fn positions(&self) -> &Positions {
        &self.positions
    }

    pub const fn cell(&self) -> Option<&PeriodicCell> {
        self.cell.as_ref()
    }

    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.cell = cell;
    }

    pub fn atom_data(&self) -> &AtomData {
        &self.atom_data
    }

    pub fn atom_data_mut(&mut self) -> &mut AtomData {
        &mut self.atom_data
    }

    pub fn set_atom_data(&mut self, atom_data: AtomData) -> Result<(), EnsembleError> {
        if !Arc::ptr_eq(atom_data.topology_arc(), self.positions.topology_arc()) {
            return Err(EnsembleError::TopologyMismatch);
        }
        self.atom_data = atom_data;
        Ok(())
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

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    pub fn view<'a>(&'a self, topology: &'a Arc<Topology>) -> Result<ModelView<'a>, EnsembleError> {
        ModelView::new(
            topology,
            &self.positions,
            self.cell.as_ref(),
            &self.atom_data,
        )
        .map_err(|_| EnsembleError::TopologyMismatch)
    }
}

/// A finite stable-order collection of non-temporal models.
#[derive(Debug, Clone)]
pub struct Ensemble {
    topology: Arc<Topology>,
    members: Vec<EnsembleMember>,
}

impl Ensemble {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            topology,
            members: Vec::new(),
        }
    }

    pub fn from_members(
        topology: Arc<Topology>,
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
        let mut ensemble = Self::new(Arc::clone(&first.topology));
        for model in models {
            if !Arc::ptr_eq(&first.topology, &model.topology) {
                return Err(EnsembleError::TopologyMismatch);
            }
            let member = EnsembleMember {
                positions: model.positions.clone(),
                cell: model.cell,
                atom_data: model.atom_data.clone(),
                weight: None,
                props: PropMap::new(),
            };
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
        let topology = Arc::new(topology_builder.build()?);
        let mut ensemble = Self::new(Arc::clone(&topology));
        for conformer in conformers {
            let positions = stage_conformer_positions(molecule.graph(), conformer)
                .map_err(|error| EnsembleError::ModelBuild(Box::new(error)))?;
            let positions = Positions::new(&topology, Quantity::new(positions, MODEL_LENGTH_UNIT))?;
            ensemble.push(EnsembleMember::new(positions))?;
        }
        Ok(ensemble)
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
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
        if !member.positions.is_compatible(&self.topology)
            || !member.atom_data.is_compatible(&self.topology)
        {
            return Err(EnsembleError::TopologyMismatch);
        }
        self.members.push(member);
        Ok(())
    }

    pub fn views(&self) -> impl ExactSizeIterator<Item = ModelView<'_>> {
        self.members.iter().map(|member| {
            member
                .view(&self.topology)
                .expect("ensemble validates shared topology allocation")
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

    /// Remaps every member to one exact target topology without renormalizing
    /// weights or changing member order.
    pub fn remap_to(
        &self,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TopologyRemapError> {
        if !mapping.is_source(&self.topology) {
            return Err(TopologyRemapError::MappingSourceMismatch);
        }
        if !mapping.is_target(target) {
            return Err(TopologyRemapError::MappingTargetMismatch);
        }
        let mut members = Vec::with_capacity(self.members.len());
        for (member_index, member) in self.members.iter().enumerate() {
            let remap_member = || -> Result<EnsembleMember, TopologyRemapError> {
                let positions = member.positions.remap_to(&self.topology, target, mapping)?;
                let atom_data = member.atom_data.remap_to(&self.topology, target, mapping)?;
                Ok(EnsembleMember {
                    positions,
                    cell: member.cell,
                    atom_data,
                    weight: member.weight,
                    props: member.props.clone(),
                })
            };
            members.push(remap_member().map_err(|error| TopologyRemapError::Member {
                member: member_index,
                error: Box::new(error),
            })?);
        }
        Ok(Self {
            topology: Arc::clone(target),
            members,
        })
    }
}

/// Focused helpers for external topology-bound state containers.
///
/// These functions validate and apply complete dense atom mappings without
/// exposing mutable topology, structure, or mapping internals. Companion
/// crates such as `kekule-traj` use them to preserve the same shared-allocation
/// and complete-state rules as Kekule's built-in structure containers.
pub mod remap {
    use super::*;

    /// Validates that `mapping` transfers complete per-atom state from `source`
    /// to `target`.
    ///
    /// Added or otherwise unmapped target atoms are rejected because complete
    /// dense state cannot invent values for them.
    pub fn validate_complete_atom_mapping(
        source: &Arc<Topology>,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<(), TopologyRemapError> {
        if !mapping.is_source(source) {
            return Err(TopologyRemapError::MappingSourceMismatch);
        }
        if !mapping.is_target(target) {
            return Err(TopologyRemapError::MappingTargetMismatch);
        }
        if let Some(target_atom) = mapping.added_atoms().first().copied() {
            return Err(TopologyRemapError::AddedAtomsRequireState { target_atom });
        }
        if mapping.atom_index_pairs().len() != target.atom_count() {
            let mapped = mapping
                .atom_pairs()
                .map(|(_, target)| target)
                .collect::<std::collections::BTreeSet<_>>();
            let target_atom = target
                .atom_ids()
                .iter()
                .copied()
                .find(|atom| !mapped.contains(atom))
                .expect("incomplete target mapping has an unmapped target atom");
            return Err(TopologyRemapError::AddedAtomsRequireState { target_atom });
        }
        Ok(())
    }

    /// Remaps one complete dense atom array into the target topology's
    /// authoritative dense order.
    pub fn dense_atom_values<T: Clone>(
        source_values: &[T],
        source: &Arc<Topology>,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Vec<T>, TopologyRemapError> {
        if source_values.len() != source.atom_count() {
            return Err(TopologyRemapError::SourceAtomCountMismatch {
                expected: source.atom_count(),
                actual: source_values.len(),
            });
        }
        validate_complete_atom_mapping(source, target, mapping)?;

        let mut values = std::iter::repeat_with(|| None)
            .take(target.atom_count())
            .collect::<Vec<Option<T>>>();
        for (source_index, target_index) in mapping.atom_index_pairs() {
            let slot = &mut values[target_index.index()];
            if slot.is_some() {
                return Err(TopologyRemapError::DuplicateTargetAssignment { target_index });
            }
            *slot = Some(source_values[source_index.index()].clone());
        }
        values
            .into_iter()
            .enumerate()
            .map(|(index, value)| {
                value.ok_or_else(|| TopologyRemapError::AddedAtomsRequireState {
                    target_atom: target.atom_ids()[index],
                })
            })
            .collect()
    }
}

/// Failure to remap complete topology-bound structure state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologyRemapError {
    /// The source state is not bound to the supplied source topology.
    SourceTopologyMismatch,
    /// The destination state is not bound to the supplied target topology.
    TargetTopologyMismatch,
    /// The mapping is not sourced from the supplied source topology.
    MappingSourceMismatch,
    /// The mapping does not target the supplied target topology.
    MappingTargetMismatch,
    /// A complete source dense array has an invalid length.
    SourceAtomCountMismatch { expected: usize, actual: usize },
    /// A target atom has no source state under this mapping.
    AddedAtomsRequireState { target_atom: InstanceAtomId },
    /// More than one source value was assigned to one target dense index.
    DuplicateTargetAssignment { target_index: TopologyAtomIndex },
    /// One ensemble member could not be remapped.
    Member {
        member: usize,
        error: Box<TopologyRemapError>,
    },
}

impl fmt::Display for TopologyRemapError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SourceTopologyMismatch => {
                formatter.write_str("source state does not belong to the supplied topology")
            }
            Self::TargetTopologyMismatch => {
                formatter.write_str("target state does not belong to the supplied topology")
            }
            Self::MappingSourceMismatch => {
                formatter.write_str("topology mapping does not match the source topology")
            }
            Self::MappingTargetMismatch => {
                formatter.write_str("topology mapping does not match the target topology")
            }
            Self::SourceAtomCountMismatch { expected, actual } => write!(
                formatter,
                "source state requires {expected} atoms, but received {actual}"
            ),
            Self::AddedAtomsRequireState { target_atom } => write!(
                formatter,
                "target atom {target_atom} has no mapped source state"
            ),
            Self::DuplicateTargetAssignment { target_index } => {
                write!(formatter, "target {target_index} received duplicate state")
            }
            Self::Member { member, error } => {
                write!(formatter, "cannot remap ensemble member {member}: {error}")
            }
        }
    }
}

impl std::error::Error for TopologyRemapError {}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum EnsembleError {
    EmptySource,
    TopologyMismatch,
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
            Self::TopologyMismatch => {
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
    use crate::core::{Atom, Conformer, Element, Molecule};
    use crate::geometry::Vector3;
    use crate::units::{ANGSTROM, KELVIN, NANOMETER, SQUARE_ANGSTROM};

    fn one_atom_topology() -> Arc<Topology> {
        let mut graph = Molecule::new();
        graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");
        let molecule = SmallMolecule::from_graph(graph);
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&molecule).unwrap();
        builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        Arc::new(builder.build().unwrap())
    }

    fn positions(topology: &Arc<Topology>, x: f64) -> Positions {
        Positions::new(
            topology,
            Quantity::new(vec![Point3::new(x, 0.0, 0.0)], ANGSTROM),
        )
        .unwrap()
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
        assert!(topology.same_layout(&independent));
        assert!(!Arc::ptr_eq(&topology, &independent));

        let wrong_positions = positions(&independent, 1.0);
        assert_eq!(
            Model::new(Arc::clone(&topology), wrong_positions),
            Err(ModelError::TopologyMismatch)
        );

        let mut model = Model::new(Arc::clone(&topology), positions(&topology, 1.0)).unwrap();
        assert!(model.atom_data().is_empty());
        let view = model.view();
        assert_eq!(
            view.positions().values().value().as_ptr(),
            model.positions().values().value().as_ptr()
        );
        let clone = model.clone();
        assert!(Arc::ptr_eq(
            &model.shared_topology(),
            &clone.shared_topology()
        ));
        let shared = model.shared_topology();
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 11.0, 12.0), ANGSTROM),
            [true; 3],
        )
        .unwrap();
        model.set_cell(Some(cell));
        assert!(Arc::ptr_eq(&model.shared_topology(), &shared));
        assert_eq!(model.cell(), Some(&cell));
    }

    #[test]
    fn ensembles_validate_shared_topology_atom_data_order_and_weights() {
        let topology = one_atom_topology();
        let atom = topology.atom_ids()[0];
        let mut first = EnsembleMember::new(positions(&topology, 1.0));
        first.set_weight(Some(1.0)).unwrap();
        first
            .atom_data_mut()
            .set_occupancy(&topology, atom, Some(0.5))
            .unwrap();

        let mut second = EnsembleMember::new(positions(&topology, 2.0));
        second.set_weight(Some(3.0)).unwrap();
        let mut invalid_weight = EnsembleMember::new(positions(&topology, 3.0));
        assert_eq!(
            invalid_weight.set_weight(Some(-1.0)),
            Err(EnsembleError::InvalidWeight)
        );
        let mut ensemble = Ensemble::from_members(Arc::clone(&topology), [first, second]).unwrap();
        ensemble.normalize_weights().unwrap();
        assert_eq!(ensemble.member(0).unwrap().weight(), Some(0.25));
        assert_eq!(ensemble.member(1).unwrap().weight(), Some(0.75));
        assert_eq!(
            ensemble
                .views()
                .map(|view| view.positions().values().value()[0].x)
                .collect::<Vec<_>>(),
            vec![1.0, 2.0]
        );
        assert_eq!(
            ensemble
                .member(0)
                .unwrap()
                .atom_data()
                .occupancy(&topology, atom)
                .unwrap(),
            Some(0.5)
        );

        let independent = one_atom_topology();
        assert_eq!(
            ensemble.push(EnsembleMember::new(positions(&independent, 3.0))),
            Err(EnsembleError::TopologyMismatch)
        );
    }

    #[test]
    fn atom_data_validates_columns_supports_mutation_and_model_topology_binding() {
        let topology = one_atom_topology();
        let atom = topology.atom_ids()[0];
        let index = topology.atom_index(atom).unwrap();
        let mut data = AtomData::new(&topology);
        assert!(data.is_empty());
        assert_eq!(data.atom_count(), 1);
        assert_eq!(data.occupancy(&topology, atom).unwrap(), None);
        data.set_occupancy(&topology, atom, Some(0.75)).unwrap();
        data.set_b_factor_at(index, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
            .unwrap();
        assert_eq!(data.occupancy_at(index).unwrap(), Some(0.75));
        assert_eq!(
            data.b_factor(&topology, atom).unwrap(),
            Some(Quantity::new(12.5, SQUARE_ANGSTROM))
        );
        data.set_b_factor(
            &topology,
            atom,
            Some(Quantity::new(0.125, NANOMETER.powi(2))),
        )
        .unwrap();
        assert!(data
            .b_factor(&topology, atom)
            .unwrap()
            .unwrap()
            .is_close(&Quantity::new(12.5, SQUARE_ANGSTROM), 1.0e-12, 1.0e-12,)
            .unwrap());
        assert!(matches!(
            data.set_b_factor(&topology, atom, Some(Quantity::new(1.0, KELVIN))),
            Err(AtomDataError::Unit(UnitError::IncompatibleUnits { .. }))
        ));
        assert!(matches!(
            data.set_b_factor(
                &topology,
                atom,
                Some(Quantity::new(f64::INFINITY, SQUARE_ANGSTROM)),
            ),
            Err(AtomDataError::NonFiniteBFactor { .. })
        ));
        assert!(matches!(
            data.set_occupancies(Vec::new()),
            Err(AtomDataError::AtomCountMismatch { .. })
        ));
        data.clear_occupancies();
        data.clear_b_factors();
        assert!(data.is_empty());
        data.set_occupancy(&topology, atom, Some(0.75)).unwrap();
        data.set_b_factor_at(index, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
            .unwrap();

        let independent = one_atom_topology();
        assert_eq!(
            data.occupancy(&independent, atom),
            Err(AtomDataError::TopologyMismatch)
        );
        let mut model = Model::new(Arc::clone(&topology), positions(&topology, 1.0)).unwrap();
        assert_eq!(
            model.set_atom_data(AtomData::new(&independent)),
            Err(ModelError::TopologyMismatch)
        );
        model.set_atom_data(data).unwrap();
        assert_eq!(model.occupancy(atom).unwrap(), Some(0.75));
        assert_eq!(
            model.b_factor(atom).unwrap(),
            Some(Quantity::new(12.5, SQUARE_ANGSTROM))
        );
    }

    #[test]
    fn model_remaps_positions_atom_data_and_cell_together() {
        let source = one_atom_topology();
        let target = one_atom_topology();
        let mapping = TopologyMapping::between_identical_layouts(&source, &target).unwrap();
        let atom = source.atom_ids()[0];
        let mut model = Model::new(Arc::clone(&source), positions(&source, 3.0)).unwrap();
        model
            .atom_data_mut()
            .set_occupancy(&source, atom, Some(0.8))
            .unwrap();
        model
            .atom_data_mut()
            .set_b_factor(&source, atom, Some(Quantity::new(21.0, SQUARE_ANGSTROM)))
            .unwrap();
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(8.0, 9.0, 10.0), ANGSTROM),
            [true; 3],
        )
        .unwrap();
        model.set_cell(Some(cell));

        let remapped = model.remap_to(&target, &mapping).unwrap();
        let target_atom = target.atom_ids()[0];
        assert_eq!(remapped.positions().values().value()[0].x, 3.0);
        assert_eq!(remapped.occupancy(target_atom).unwrap(), Some(0.8));
        assert_eq!(
            remapped.b_factor(target_atom).unwrap(),
            Some(Quantity::new(21.0, SQUARE_ANGSTROM))
        );
        assert_eq!(remapped.cell(), Some(&cell));
        let view = remapped.view();
        assert_eq!(view.position(target_atom).unwrap().value().x, 3.0);
        assert_eq!(view.occupancy(target_atom).unwrap(), Some(0.8));
        assert_eq!(view.atom_data(), remapped.atom_data());
    }

    #[test]
    fn ensemble_from_conformers_preserves_source_order_without_copying_conformers_to_topology() {
        let mut graph = Molecule::new();
        let atom = graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");
        let mut first = Conformer::new(ANGSTROM).unwrap();
        first
            .set_position(atom, Quantity::new(Point3::new(1.0, 0.0, 0.0), ANGSTROM))
            .unwrap();
        let first = graph.add_conformer(first).unwrap();
        let mut second = Conformer::new(ANGSTROM).unwrap();
        second
            .set_position(atom, Quantity::new(Point3::new(2.0, 0.0, 0.0), ANGSTROM))
            .unwrap();
        let second = graph.add_conformer(second).unwrap();
        let molecule = SmallMolecule::from_graph(graph);

        let ensemble =
            Ensemble::from_small_molecule_conformers(&molecule, [second, first]).unwrap();
        assert_eq!(
            ensemble
                .views()
                .map(|view| view.positions().values().value()[0].x)
                .collect::<Vec<_>>(),
            vec![2.0, 1.0]
        );
        assert_eq!(
            ensemble
                .topology()
                .definitions()
                .next()
                .unwrap()
                .1
                .graph()
                .conformers()
                .count(),
            0
        );
    }
}
