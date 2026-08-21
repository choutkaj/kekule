use std::collections::BTreeMap;
use std::fmt;
use std::sync::Arc;

use crate::topology::{
    InstanceAtomId, InstanceBondId, Topology, TopologyAtomIndex, TopologyBondIndex, TopologyMapping,
};
use crate::units::{Quantity, Unit, UnitError, DIMENSIONLESS, SQUARE_ANGSTROM};

use super::{remap, TopologyRemapError};

const MAX_PROPERTY_NAME_LEN: usize = 128;

#[derive(Debug, Clone, PartialEq)]
pub(super) struct ScalarPropertyColumn {
    pub(super) unit: Unit,
    values: Vec<Option<f64>>,
}

impl ScalarPropertyColumn {
    fn quantity(&self) -> Quantity<&[Option<f64>]> {
        Quantity::new(self.values.as_slice(), self.unit)
    }

    fn value(&self, index: usize) -> Option<f64> {
        self.values[index]
    }

    fn stage_value(
        &self,
        index: usize,
        value: Option<Quantity<f64>>,
    ) -> Result<Option<f64>, ScalarPropertyColumnError> {
        stage_property_value(value, self.unit, index)
    }

    /// Replaces one already-validated dense value and reports whether the
    /// column became wholly absent.
    fn replace_value(&mut self, index: usize, value: Option<f64>) -> bool {
        self.values[index] = value;
        self.values.iter().all(Option::is_none)
    }

    fn with_value(unit: Unit, len: usize, index: usize, value: f64) -> Self {
        let mut values = vec![None; len];
        values[index] = Some(value);
        Self { unit, values }
    }

    fn from_values(unit: Unit, values: Vec<Option<f64>>) -> Option<Self> {
        values
            .iter()
            .any(Option::is_some)
            .then_some(Self { unit, values })
    }
}

#[derive(Debug, Clone, PartialEq)]
enum ScalarPropertyColumnError {
    ValueCountMismatch { expected: usize, actual: usize },
    NonFiniteValue { index: usize },
    Unit(UnitError),
}

fn valid_property_name(name: &str) -> bool {
    if name.is_empty() || name.len() > MAX_PROPERTY_NAME_LEN {
        return false;
    }
    let mut bytes = name.bytes();
    let Some(first) = bytes.next() else {
        return false;
    };
    (first.is_ascii_alphabetic() || first == b'_')
        && bytes.all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.'))
}

fn stage_property_column<T>(
    values: Quantity<T>,
    expected: usize,
    stored_unit: Option<Unit>,
) -> Result<Option<ScalarPropertyColumn>, ScalarPropertyColumnError>
where
    T: AsRef<[Option<f64>]>,
{
    let source = values.value().as_ref();
    if source.len() != expected {
        return Err(ScalarPropertyColumnError::ValueCountMismatch {
            expected,
            actual: source.len(),
        });
    }
    let unit = stored_unit.unwrap_or(values.unit());
    let factor = values
        .unit()
        .conversion_factor_to(unit)
        .map_err(ScalarPropertyColumnError::Unit)?;
    let converted = source
        .iter()
        .copied()
        .enumerate()
        .map(|(index, value)| {
            let value = value.map(|value| value * factor);
            if value.is_some_and(|value| !value.is_finite()) {
                return Err(ScalarPropertyColumnError::NonFiniteValue { index });
            }
            Ok(value)
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(ScalarPropertyColumn::from_values(unit, converted))
}

fn converted_property_value(value: Quantity<f64>, stored_unit: Unit) -> Result<f64, UnitError> {
    let value = value.to_unit(stored_unit)?.to_value();
    Ok(value)
}

fn stage_property_value(
    value: Option<Quantity<f64>>,
    stored_unit: Unit,
    index: usize,
) -> Result<Option<f64>, ScalarPropertyColumnError> {
    let converted = value
        .map(|value| converted_property_value(value, stored_unit))
        .transpose()
        .map_err(ScalarPropertyColumnError::Unit)?;
    if converted.is_some_and(|value| !value.is_finite()) {
        return Err(ScalarPropertyColumnError::NonFiniteValue { index });
    }
    Ok(converted)
}

fn replace_named_property_column<T>(
    properties: &mut BTreeMap<String, ScalarPropertyColumn>,
    name: &str,
    expected: usize,
    values: Quantity<T>,
) -> Result<(), ScalarPropertyColumnError>
where
    T: AsRef<[Option<f64>]>,
{
    let stored_unit = properties.get(name).map(|column| column.unit);
    let staged = stage_property_column(values, expected, stored_unit)?;
    match staged {
        Some(column) => {
            properties.insert(name.to_owned(), column);
        }
        None => {
            properties.remove(name);
        }
    }
    Ok(())
}

fn set_named_property_column_value(
    properties: &mut BTreeMap<String, ScalarPropertyColumn>,
    name: &str,
    len: usize,
    index: usize,
    value: Option<Quantity<f64>>,
) -> Result<(), ScalarPropertyColumnError> {
    let remove = if let Some(column) = properties.get_mut(name) {
        let value = column.stage_value(index, value)?;
        column.replace_value(index, value)
    } else {
        let Some(value) = value else {
            return Ok(());
        };
        let unit = value.unit();
        let value = stage_property_value(Some(value), unit, index)?
            .expect("present property value remains present after conversion");
        properties.insert(
            name.to_owned(),
            ScalarPropertyColumn::with_value(unit, len, index, value),
        );
        false
    };
    if remove {
        properties.remove(name);
    }
    Ok(())
}

fn set_optional_property_column_value(
    column: &mut Option<ScalarPropertyColumn>,
    len: usize,
    index: usize,
    value: Option<Quantity<f64>>,
    canonical_unit: Unit,
) -> Result<(), ScalarPropertyColumnError> {
    let staged = match column.as_ref() {
        Some(column) => column.stage_value(index, value)?,
        None => stage_property_value(value, canonical_unit, index)?,
    };
    if let Some(existing) = column.as_mut() {
        let remove = existing.replace_value(index, staged);
        if remove {
            *column = None;
        }
    } else if let Some(value) = staged {
        *column = Some(ScalarPropertyColumn::with_value(
            canonical_unit,
            len,
            index,
            value,
        ));
    }
    Ok(())
}

fn remap_atom_property_column(
    column: &ScalarPropertyColumn,
    source: &Arc<Topology>,
    target: &Arc<Topology>,
    mapping: &TopologyMapping,
) -> Result<Option<ScalarPropertyColumn>, TopologyRemapError> {
    let values = remap::dense_atom_values(&column.values, source, target, mapping)?;
    Ok(ScalarPropertyColumn::from_values(column.unit, values))
}

fn remap_bond_property_column(
    column: &ScalarPropertyColumn,
    source: &Arc<Topology>,
    target: &Arc<Topology>,
    mapping: &TopologyMapping,
) -> Result<Option<ScalarPropertyColumn>, TopologyRemapError> {
    let values = remap::dense_bond_values(&column.values, source, target, mapping)?;
    Ok(ScalarPropertyColumn::from_values(column.unit, values))
}

/// Topology-bound model-level per-atom data in dense atom order.
///
/// Canonical occupancy and B-factor columns retain dedicated APIs and fixed
/// semantic units. Custom scalar properties occupy a separate namespace. Every
/// column is dense and optional, so a field absent for all atoms requires no
/// per-atom allocation.
#[derive(Debug, Clone)]
pub struct AtomData {
    topology: Arc<Topology>,
    pub(super) occupancies: Option<ScalarPropertyColumn>,
    pub(super) b_factors: Option<ScalarPropertyColumn>,
    properties: BTreeMap<String, ScalarPropertyColumn>,
}

impl PartialEq for AtomData {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.occupancies == other.occupancies
            && self.b_factors == other.b_factors
            && self.properties == other.properties
    }
}

impl AtomData {
    /// Creates atom data with no allocated scientific columns.
    pub fn new(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            occupancies: None,
            b_factors: None,
            properties: BTreeMap::new(),
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

    pub(super) fn topology_arc(&self) -> &Arc<Topology> {
        &self.topology
    }

    /// Returns the atom count of the bound topology.
    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    /// Returns whether every supported scientific column is wholly absent.
    pub fn is_empty(&self) -> bool {
        self.occupancies.is_none() && self.b_factors.is_none() && self.properties.is_empty()
    }

    pub fn occupancies(&self) -> Option<&[Option<f64>]> {
        self.occupancies
            .as_ref()
            .map(|column| column.values.as_slice())
    }

    /// Returns the dense B-factor column in canonical square angstroms.
    pub fn b_factors(&self) -> Option<Quantity<&[Option<f64>]>> {
        self.b_factors.as_ref().map(ScalarPropertyColumn::quantity)
    }

    /// Iterates custom scalar properties in stable name order.
    pub fn properties(&self) -> impl ExactSizeIterator<Item = (&str, Quantity<&[Option<f64>]>)> {
        self.properties
            .iter()
            .map(|(name, column)| (name.as_str(), column.quantity()))
    }

    /// Returns a complete custom property column, or `None` when absent.
    pub fn property(&self, name: &str) -> Result<Option<Quantity<&[Option<f64>]>>, AtomDataError> {
        validate_atom_property_name(name)?;
        Ok(self
            .properties
            .get(name)
            .map(ScalarPropertyColumn::quantity))
    }

    /// Replaces a complete custom property transactionally.
    ///
    /// Existing properties retain their stored unit; compatible input values
    /// are converted into it. An all-missing column removes the property.
    pub fn set_property<T>(&mut self, name: &str, values: Quantity<T>) -> Result<(), AtomDataError>
    where
        T: AsRef<[Option<f64>]>,
    {
        validate_atom_property_name(name)?;
        let atom_count = self.atom_count();
        replace_named_property_column(&mut self.properties, name, atom_count, values)
            .map_err(|error| atom_property_column_error(name, error))
    }

    /// Removes a custom property, returning whether it was present.
    pub fn remove_property(&mut self, name: &str) -> Result<bool, AtomDataError> {
        validate_atom_property_name(name)?;
        Ok(self.properties.remove(name).is_some())
    }

    /// Returns one unit-aware custom property value by semantic atom ID.
    pub fn property_value(
        &self,
        topology: &Arc<Topology>,
        name: &str,
        atom: InstanceAtomId,
    ) -> Result<Option<Quantity<f64>>, AtomDataError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(AtomDataError::InvalidAtomId(atom))?;
        self.property_value_at(name, index)
    }

    /// Returns one unit-aware custom property value by dense atom index.
    pub fn property_value_at(
        &self,
        name: &str,
        index: TopologyAtomIndex,
    ) -> Result<Option<Quantity<f64>>, AtomDataError> {
        validate_atom_property_name(name)?;
        validate_index(self.atom_count(), index)?;
        Ok(self.properties.get(name).and_then(|column| {
            column
                .value(index.index())
                .map(|value| Quantity::new(value, column.unit))
        }))
    }

    /// Sets one custom property value by semantic atom ID.
    pub fn set_property_value(
        &mut self,
        topology: &Arc<Topology>,
        name: &str,
        atom: InstanceAtomId,
        value: Option<Quantity<f64>>,
    ) -> Result<(), AtomDataError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(AtomDataError::InvalidAtomId(atom))?;
        self.set_property_value_at(name, index, value)
    }

    /// Sets one custom property value by dense atom index.
    ///
    /// A missing property is created from a present value's unit. Setting
    /// `None` on an absent property is a no-op, and clearing the final present
    /// value removes the property.
    pub fn set_property_value_at(
        &mut self,
        name: &str,
        index: TopologyAtomIndex,
        value: Option<Quantity<f64>>,
    ) -> Result<(), AtomDataError> {
        validate_atom_property_name(name)?;
        validate_index(self.atom_count(), index)?;
        let atom_count = self.atom_count();
        set_named_property_column_value(
            &mut self.properties,
            name,
            atom_count,
            index.index(),
            value,
        )
        .map_err(|error| atom_property_column_error(name, error))
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
        validate_index(self.atom_count(), index)?;
        Ok(self
            .occupancies
            .as_ref()
            .and_then(|column| column.value(index.index())))
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
        validate_index(self.atom_count(), index)?;
        Ok(self.b_factors.as_ref().and_then(|column| {
            column
                .value(index.index())
                .map(|value| Quantity::new(value, column.unit))
        }))
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
        let len = self.atom_count();
        set_optional_property_column_value(
            &mut self.occupancies,
            len,
            index.index(),
            value.map(|value| Quantity::new(value, DIMENSIONLESS)),
            DIMENSIONLESS,
        )
        .map_err(|error| atom_data_column_error(AtomDataField::Occupancy, error))
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
        let len = self.atom_count();
        set_optional_property_column_value(
            &mut self.b_factors,
            len,
            index.index(),
            value,
            SQUARE_ANGSTROM,
        )
        .map_err(|error| atom_data_column_error(AtomDataField::BFactor, error))
    }

    /// Replaces the complete occupancy column transactionally. An all-absent
    /// column is normalized to no allocation.
    pub fn set_occupancies<T>(&mut self, values: T) -> Result<(), AtomDataError>
    where
        T: AsRef<[Option<f64>]>,
    {
        let staged = stage_property_column(
            Quantity::new(values, DIMENSIONLESS),
            self.atom_count(),
            Some(DIMENSIONLESS),
        )
        .map_err(|error| atom_data_column_error(AtomDataField::Occupancy, error))?;
        self.occupancies = staged;
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
        let staged = stage_property_column(values, self.atom_count(), Some(SQUARE_ANGSTROM))
            .map_err(|error| atom_data_column_error(AtomDataField::BFactor, error))?;
        self.b_factors = staged;
        Ok(())
    }

    /// Clears the complete B-factor column.
    pub fn clear_b_factors(&mut self) {
        self.b_factors = None;
    }

    /// Remaps every canonical and custom property column through checked
    /// topology lineage.
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
            .map(|column| remap_atom_property_column(column, source, target, mapping))
            .transpose()?
            .flatten();
        let b_factors = self
            .b_factors
            .as_ref()
            .map(|column| remap_atom_property_column(column, source, target, mapping))
            .transpose()?
            .flatten();
        let mut properties = BTreeMap::new();
        for (name, column) in &self.properties {
            if let Some(column) = remap_atom_property_column(column, source, target, mapping)? {
                properties.insert(name.clone(), column);
            }
        }
        Ok(Self {
            topology: Arc::clone(target),
            occupancies,
            b_factors,
            properties,
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

fn atom_data_column_error(field: AtomDataField, error: ScalarPropertyColumnError) -> AtomDataError {
    match error {
        ScalarPropertyColumnError::ValueCountMismatch { expected, actual } => {
            AtomDataError::AtomCountMismatch { expected, actual }
        }
        ScalarPropertyColumnError::NonFiniteValue { index } => match field {
            AtomDataField::Occupancy => AtomDataError::NonFiniteOccupancy {
                index: TopologyAtomIndex::new(index as u32),
            },
            AtomDataField::BFactor => AtomDataError::NonFiniteBFactor {
                index: TopologyAtomIndex::new(index as u32),
            },
        },
        ScalarPropertyColumnError::Unit(error) => AtomDataError::Unit(error),
    }
}

fn validate_index(len: usize, index: TopologyAtomIndex) -> Result<(), AtomDataError> {
    if index.index() >= len {
        return Err(AtomDataError::InvalidAtomIndex(index));
    }
    Ok(())
}

fn validate_atom_property_name(name: &str) -> Result<(), AtomDataError> {
    if !valid_property_name(name) {
        return Err(AtomDataError::InvalidPropertyName {
            name: name.to_owned(),
        });
    }
    if name.eq_ignore_ascii_case("occupancy") || name.eq_ignore_ascii_case("b_factor") {
        return Err(AtomDataError::ReservedPropertyName {
            name: name.to_owned(),
        });
    }
    Ok(())
}

fn atom_property_column_error(name: &str, error: ScalarPropertyColumnError) -> AtomDataError {
    match error {
        ScalarPropertyColumnError::ValueCountMismatch { expected, actual } => {
            AtomDataError::PropertyValueCountMismatch {
                property: name.to_owned(),
                expected,
                actual,
            }
        }
        ScalarPropertyColumnError::NonFiniteValue { index } => {
            AtomDataError::NonFinitePropertyValue {
                property: name.to_owned(),
                index: TopologyAtomIndex::new(index as u32),
            }
        }
        ScalarPropertyColumnError::Unit(error) => AtomDataError::PropertyUnit {
            property: name.to_owned(),
            error: Box::new(error),
        },
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum AtomDataError {
    TopologyMismatch,
    AtomCountMismatch {
        expected: usize,
        actual: usize,
    },
    InvalidAtomId(InstanceAtomId),
    InvalidAtomIndex(TopologyAtomIndex),
    NonFiniteOccupancy {
        index: TopologyAtomIndex,
    },
    NonFiniteBFactor {
        index: TopologyAtomIndex,
    },
    InvalidPropertyName {
        name: String,
    },
    ReservedPropertyName {
        name: String,
    },
    PropertyValueCountMismatch {
        property: String,
        expected: usize,
        actual: usize,
    },
    NonFinitePropertyValue {
        property: String,
        index: TopologyAtomIndex,
    },
    PropertyUnit {
        property: String,
        error: Box<UnitError>,
    },
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
            Self::InvalidPropertyName { name } => write!(
                formatter,
                "invalid atom property name {name:?}; use a 1-{MAX_PROPERTY_NAME_LEN} character ASCII identifier"
            ),
            Self::ReservedPropertyName { name } => {
                write!(formatter, "atom property name {name:?} is reserved")
            }
            Self::PropertyValueCountMismatch {
                property,
                expected,
                actual,
            } => write!(
                formatter,
                "atom property {property:?} requires {expected} values, but received {actual}"
            ),
            Self::NonFinitePropertyValue { property, index } => write!(
                formatter,
                "atom property {property:?} at {index} must be finite"
            ),
            Self::PropertyUnit { property, error } => {
                write!(formatter, "invalid unit for atom property {property:?}: {error}")
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

/// Topology-bound model-level per-bond annotations in authoritative dense bond
/// order.
///
/// Bond data has no canonical scientific fields. It stores only conservative,
/// user- or analysis-defined unit-aware scalar properties.
#[derive(Debug, Clone)]
pub struct BondData {
    topology: Arc<Topology>,
    properties: BTreeMap<String, ScalarPropertyColumn>,
}

impl PartialEq for BondData {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.properties == other.properties
    }
}

impl BondData {
    /// Creates empty bond data bound to the exact shared topology allocation.
    pub fn new(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            properties: BTreeMap::new(),
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

    pub(super) fn topology_arc(&self) -> &Arc<Topology> {
        &self.topology
    }

    /// Returns the bond count of the bound topology.
    pub fn bond_count(&self) -> usize {
        self.topology.bond_count()
    }

    pub fn is_empty(&self) -> bool {
        self.properties.is_empty()
    }

    /// Iterates custom scalar properties in stable name order.
    pub fn properties(&self) -> impl ExactSizeIterator<Item = (&str, Quantity<&[Option<f64>]>)> {
        self.properties
            .iter()
            .map(|(name, column)| (name.as_str(), column.quantity()))
    }

    /// Returns a complete custom property column, or `None` when absent.
    pub fn property(&self, name: &str) -> Result<Option<Quantity<&[Option<f64>]>>, BondDataError> {
        validate_bond_property_name(name)?;
        Ok(self
            .properties
            .get(name)
            .map(ScalarPropertyColumn::quantity))
    }

    /// Replaces a complete custom property transactionally.
    ///
    /// Existing properties retain their stored unit; compatible input values
    /// are converted into it. An all-missing column removes the property.
    pub fn set_property<T>(&mut self, name: &str, values: Quantity<T>) -> Result<(), BondDataError>
    where
        T: AsRef<[Option<f64>]>,
    {
        validate_bond_property_name(name)?;
        let bond_count = self.bond_count();
        replace_named_property_column(&mut self.properties, name, bond_count, values)
            .map_err(|error| bond_property_column_error(name, error))
    }

    /// Removes a custom property, returning whether it was present.
    pub fn remove_property(&mut self, name: &str) -> Result<bool, BondDataError> {
        validate_bond_property_name(name)?;
        Ok(self.properties.remove(name).is_some())
    }

    /// Returns one unit-aware custom property value by semantic bond ID.
    pub fn property_value(
        &self,
        topology: &Arc<Topology>,
        name: &str,
        bond: InstanceBondId,
    ) -> Result<Option<Quantity<f64>>, BondDataError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .bond_index(bond)
            .ok_or(BondDataError::InvalidBondId(bond))?;
        self.property_value_at(name, index)
    }

    /// Returns one unit-aware custom property value by dense bond index.
    pub fn property_value_at(
        &self,
        name: &str,
        index: TopologyBondIndex,
    ) -> Result<Option<Quantity<f64>>, BondDataError> {
        validate_bond_property_name(name)?;
        validate_bond_index(self.bond_count(), index)?;
        Ok(self.properties.get(name).and_then(|column| {
            column
                .value(index.index())
                .map(|value| Quantity::new(value, column.unit))
        }))
    }

    /// Sets one custom property value by semantic bond ID.
    pub fn set_property_value(
        &mut self,
        topology: &Arc<Topology>,
        name: &str,
        bond: InstanceBondId,
        value: Option<Quantity<f64>>,
    ) -> Result<(), BondDataError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .bond_index(bond)
            .ok_or(BondDataError::InvalidBondId(bond))?;
        self.set_property_value_at(name, index, value)
    }

    /// Sets one custom property value by dense bond index.
    ///
    /// A missing property is created from a present value's unit. Setting
    /// `None` on an absent property is a no-op, and clearing the final present
    /// value removes the property.
    pub fn set_property_value_at(
        &mut self,
        name: &str,
        index: TopologyBondIndex,
        value: Option<Quantity<f64>>,
    ) -> Result<(), BondDataError> {
        validate_bond_property_name(name)?;
        validate_bond_index(self.bond_count(), index)?;
        let bond_count = self.bond_count();
        set_named_property_column_value(
            &mut self.properties,
            name,
            bond_count,
            index.index(),
            value,
        )
        .map_err(|error| bond_property_column_error(name, error))
    }

    /// Remaps every custom property through checked topology lineage.
    pub fn remap_to(
        &self,
        source: &Arc<Topology>,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TopologyRemapError> {
        if !self.is_compatible(source) {
            return Err(TopologyRemapError::SourceTopologyMismatch);
        }
        let mut properties = BTreeMap::new();
        for (name, column) in &self.properties {
            if let Some(column) = remap_bond_property_column(column, source, target, mapping)? {
                properties.insert(name.clone(), column);
            }
        }
        Ok(Self {
            topology: Arc::clone(target),
            properties,
        })
    }

    fn ensure_compatible(&self, topology: &Arc<Topology>) -> Result<(), BondDataError> {
        if !self.is_compatible(topology) {
            return Err(BondDataError::TopologyMismatch);
        }
        Ok(())
    }
}

fn validate_bond_property_name(name: &str) -> Result<(), BondDataError> {
    if !valid_property_name(name) {
        return Err(BondDataError::InvalidPropertyName {
            name: name.to_owned(),
        });
    }
    Ok(())
}

fn validate_bond_index(len: usize, index: TopologyBondIndex) -> Result<(), BondDataError> {
    if index.index() >= len {
        return Err(BondDataError::InvalidBondIndex(index));
    }
    Ok(())
}

fn bond_property_column_error(name: &str, error: ScalarPropertyColumnError) -> BondDataError {
    match error {
        ScalarPropertyColumnError::ValueCountMismatch { expected, actual } => {
            BondDataError::PropertyValueCountMismatch {
                property: name.to_owned(),
                expected,
                actual,
            }
        }
        ScalarPropertyColumnError::NonFiniteValue { index } => {
            BondDataError::NonFinitePropertyValue {
                property: name.to_owned(),
                index: TopologyBondIndex::new(index as u32),
            }
        }
        ScalarPropertyColumnError::Unit(error) => BondDataError::PropertyUnit {
            property: name.to_owned(),
            error: Box::new(error),
        },
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum BondDataError {
    TopologyMismatch,
    InvalidBondId(InstanceBondId),
    InvalidBondIndex(TopologyBondIndex),
    InvalidPropertyName {
        name: String,
    },
    PropertyValueCountMismatch {
        property: String,
        expected: usize,
        actual: usize,
    },
    NonFinitePropertyValue {
        property: String,
        index: TopologyBondIndex,
    },
    PropertyUnit {
        property: String,
        error: Box<UnitError>,
    },
}

impl fmt::Display for BondDataError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("bond data belongs to a different topology")
            }
            Self::InvalidBondId(bond) => write!(formatter, "invalid topology bond: {bond}"),
            Self::InvalidBondIndex(index) => write!(formatter, "invalid {index}"),
            Self::InvalidPropertyName { name } => write!(
                formatter,
                "invalid bond property name {name:?}; use a 1-{MAX_PROPERTY_NAME_LEN} character ASCII identifier"
            ),
            Self::PropertyValueCountMismatch {
                property,
                expected,
                actual,
            } => write!(
                formatter,
                "bond property {property:?} requires {expected} values, but received {actual}"
            ),
            Self::NonFinitePropertyValue { property, index } => write!(
                formatter,
                "bond property {property:?} at {index} must be finite"
            ),
            Self::PropertyUnit { property, error } => {
                write!(formatter, "invalid unit for bond property {property:?}: {error}")
            }
        }
    }
}

impl std::error::Error for BondDataError {}
