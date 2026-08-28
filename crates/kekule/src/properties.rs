//! Unified scalar and columnar annotations for canonical Kekule objects.

use std::collections::BTreeMap;
use std::fmt;
use std::str::FromStr;

use crate::units::{Quantity, Unit, UnitError, DIMENSIONLESS, SQUARE_NANOMETER};

pub const MAX_PROPERTY_KEY_LEN: usize = 128;

/// A validated, deterministic property identifier.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PropertyKey(String);

impl PropertyKey {
    pub fn new(value: impl Into<String>) -> Result<Self, PropertyError> {
        let value = value.into();
        if !valid_property_key(&value) {
            return Err(PropertyError::InvalidKey(value));
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl AsRef<str> for PropertyKey {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

impl fmt::Display for PropertyKey {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl TryFrom<&str> for PropertyKey {
    type Error = PropertyError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Self::new(value)
    }
}

impl TryFrom<String> for PropertyKey {
    type Error = PropertyError;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        Self::new(value)
    }
}

impl FromStr for PropertyKey {
    type Err = PropertyError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        Self::new(value)
    }
}

fn valid_property_key(value: &str) -> bool {
    if value.is_empty() || value.len() > MAX_PROPERTY_KEY_LEN {
        return false;
    }
    let mut bytes = value.bytes();
    let Some(first) = bytes.next() else {
        return false;
    };
    (first.is_ascii_alphabetic() || first == b'_')
        && bytes.all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.'))
}

/// One generic scalar annotation.
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyValue {
    Bool(bool),
    Int(i64),
    Real { value: f64, unit: Unit },
    String(String),
}

impl PropertyValue {
    pub fn real(value: f64, unit: Unit) -> Result<Self, PropertyError> {
        let value = Self::Real { value, unit };
        value.validate()?;
        Ok(value)
    }

    fn validate(&self) -> Result<(), PropertyError> {
        if let Self::Real { value, .. } = self {
            if !value.is_finite() {
                return Err(PropertyError::NonFiniteValue { index: None });
            }
        }
        Ok(())
    }
}

/// One homogeneous optional-valued property column.
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyColumn {
    Bool(Vec<Option<bool>>),
    Int(Vec<Option<i64>>),
    Real {
        unit: Unit,
        values: Vec<Option<f64>>,
    },
    String(Vec<Option<String>>),
}

impl PropertyColumn {
    pub fn len(&self) -> usize {
        match self {
            Self::Bool(values) => values.len(),
            Self::Int(values) => values.len(),
            Self::Real { values, .. } => values.len(),
            Self::String(values) => values.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn is_all_missing(&self) -> bool {
        match self {
            Self::Bool(values) => values.iter().all(Option::is_none),
            Self::Int(values) => values.iter().all(Option::is_none),
            Self::Real { values, .. } => values.iter().all(Option::is_none),
            Self::String(values) => values.iter().all(Option::is_none),
        }
    }

    pub fn value(&self, index: usize) -> Result<Option<PropertyValue>, PropertyError> {
        if index >= self.len() {
            return Err(PropertyError::InvalidIndex {
                len: self.len(),
                index,
            });
        }
        Ok(match self {
            Self::Bool(values) => values[index].map(PropertyValue::Bool),
            Self::Int(values) => values[index].map(PropertyValue::Int),
            Self::Real { unit, values } => {
                values[index].map(|value| PropertyValue::Real { value, unit: *unit })
            }
            Self::String(values) => values[index].clone().map(PropertyValue::String),
        })
    }

    fn validate(&self) -> Result<(), PropertyError> {
        if let Self::Real { values, .. } = self {
            if let Some(index) = values
                .iter()
                .position(|value| value.is_some_and(|value| !value.is_finite()))
            {
                return Err(PropertyError::NonFiniteValue { index: Some(index) });
            }
        }
        Ok(())
    }

    fn converted_to(self, unit: Unit) -> Result<Self, PropertyError> {
        let Self::Real {
            unit: source_unit,
            values,
        } = self
        else {
            return Ok(self);
        };
        let factor = source_unit.conversion_factor_to(unit)?;
        let values = values
            .into_iter()
            .enumerate()
            .map(|(index, value)| {
                let converted = value.map(|value| value * factor);
                if converted.is_some_and(|value| !value.is_finite()) {
                    return Err(PropertyError::NonFiniteValue { index: Some(index) });
                }
                Ok(converted)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self::Real { unit, values })
    }

    fn select_indices(&self, indices: &[usize]) -> Self {
        match self {
            Self::Bool(values) => Self::Bool(indices.iter().map(|index| values[*index]).collect()),
            Self::Int(values) => Self::Int(indices.iter().map(|index| values[*index]).collect()),
            Self::Real { unit, values } => Self::Real {
                unit: *unit,
                values: indices.iter().map(|index| values[*index]).collect(),
            },
            Self::String(values) => {
                Self::String(indices.iter().map(|index| values[*index].clone()).collect())
            }
        }
    }

    fn resize_missing(&mut self, len: usize) {
        match self {
            Self::Bool(values) => values.resize(len, None),
            Self::Int(values) => values.resize(len, None),
            Self::Real { values, .. } => values.resize(len, None),
            Self::String(values) => values.resize(len, None),
        }
    }

    fn clear(&mut self, index: usize) {
        match self {
            Self::Bool(values) => values[index] = None,
            Self::Int(values) => values[index] = None,
            Self::Real { values, .. } => values[index] = None,
            Self::String(values) => values[index] = None,
        }
    }
}

/// Columnar properties for one detached homogeneous entity domain.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyTable {
    len: usize,
    columns: BTreeMap<PropertyKey, PropertyColumn>,
}

impl PropertyTable {
    pub fn new(len: usize) -> Self {
        Self {
            len,
            columns: BTreeMap::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn has_data(&self) -> bool {
        !self.columns.is_empty()
    }

    pub fn row_has_data(&self, index: usize) -> Result<bool, PropertyError> {
        self.validate_index(index)?;
        Ok(self
            .columns
            .values()
            .any(|column| column.value(index).is_ok_and(|value| value.is_some())))
    }

    pub fn get(&self, key: &PropertyKey) -> Option<&PropertyColumn> {
        self.columns.get(key)
    }

    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyKey, &PropertyColumn)> + DoubleEndedIterator {
        self.columns.iter()
    }

    /// Inserts or replaces a complete column transactionally.
    ///
    /// Replacing a real column preserves its existing storage unit and converts
    /// compatible input values into it. An all-missing input removes the key.
    pub fn insert(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, PropertyError> {
        self.insert_validated(key, column)
    }

    fn insert_validated(
        &mut self,
        key: PropertyKey,
        mut column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, PropertyError> {
        if column.len() != self.len {
            return Err(PropertyError::LengthMismatch {
                expected: self.len,
                actual: column.len(),
            });
        }
        column.validate()?;
        if column.is_all_missing() {
            return Ok(self.columns.remove(&key));
        }
        if let Some(PropertyColumn::Real { unit, .. }) = self.columns.get(&key) {
            column = match column {
                PropertyColumn::Real { .. } => column.converted_to(*unit)?,
                _ => return Err(PropertyError::TypeMismatch { key }),
            };
        } else if let Some(existing) = self.columns.get(&key) {
            if std::mem::discriminant(existing) != std::mem::discriminant(&column) {
                return Err(PropertyError::TypeMismatch { key });
            }
        }
        Ok(self.columns.insert(key, column))
    }

    pub fn remove(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.columns.remove(key)
    }

    pub fn value(
        &self,
        key: &PropertyKey,
        index: usize,
    ) -> Result<Option<PropertyValue>, PropertyError> {
        self.validate_index(index)?;
        self.columns
            .get(key)
            .map(|column| column.value(index))
            .transpose()
            .map(Option::flatten)
    }

    /// Sets or clears one cell transactionally, removing an all-missing column.
    pub fn set_value(
        &mut self,
        key: PropertyKey,
        index: usize,
        value: Option<PropertyValue>,
    ) -> Result<(), PropertyError> {
        self.set_value_validated(key, index, value)
    }

    fn set_value_validated(
        &mut self,
        key: PropertyKey,
        index: usize,
        value: Option<PropertyValue>,
    ) -> Result<(), PropertyError> {
        self.validate_index(index)?;
        if let Some(value) = &value {
            value.validate()?;
        }
        let Some(column) = self.columns.get(&key) else {
            let Some(value) = value else {
                return Ok(());
            };
            let mut column = match &value {
                PropertyValue::Bool(_) => PropertyColumn::Bool(vec![None; self.len]),
                PropertyValue::Int(_) => PropertyColumn::Int(vec![None; self.len]),
                PropertyValue::Real { unit, .. } => PropertyColumn::Real {
                    unit: *unit,
                    values: vec![None; self.len],
                },
                PropertyValue::String(_) => PropertyColumn::String(vec![None; self.len]),
            };
            assign_value(&mut column, &key, index, Some(value))?;
            self.columns.insert(key, column);
            return Ok(());
        };
        let mut staged = column.clone();
        assign_value(&mut staged, &key, index, value)?;
        if staged.is_all_missing() {
            self.columns.remove(&key);
        } else {
            self.columns.insert(key, staged);
        }
        Ok(())
    }

    pub fn clear_value(&mut self, key: PropertyKey, index: usize) -> Result<(), PropertyError> {
        self.set_value(key, index, None)
    }

    pub fn select_indices(&self, indices: &[usize]) -> Result<Self, PropertyError> {
        for index in indices {
            self.validate_index(*index)?;
        }
        let columns = self
            .columns
            .iter()
            .filter_map(|(key, column)| {
                let selected = column.select_indices(indices);
                (!selected.is_all_missing()).then(|| (key.clone(), selected))
            })
            .collect();
        Ok(Self {
            len: indices.len(),
            columns,
        })
    }

    pub(crate) fn select_optional_indices(
        &self,
        indices: &[Option<usize>],
    ) -> Result<Self, PropertyError> {
        for index in indices.iter().flatten() {
            self.validate_index(*index)?;
        }
        let mut selected = Self::new(indices.len());
        for (key, column) in &self.columns {
            for (target, source) in indices.iter().enumerate() {
                let value = source
                    .map(|source| column.value(source))
                    .transpose()?
                    .flatten();
                selected.set_value(key.clone(), target, value)?;
            }
        }
        Ok(selected)
    }

    pub(crate) fn resize_missing(&mut self, len: usize) {
        self.len = len;
        for column in self.columns.values_mut() {
            column.resize_missing(len);
        }
    }

    pub(crate) fn clear_index(&mut self, index: usize) {
        if index >= self.len {
            return;
        }
        self.columns.retain(|_, column| {
            column.clear(index);
            !column.is_all_missing()
        });
    }

    fn validate_index(&self, index: usize) -> Result<(), PropertyError> {
        if index >= self.len {
            return Err(PropertyError::InvalidIndex {
                len: self.len,
                index,
            });
        }
        Ok(())
    }
}

fn is_reserved_realization_atom_key(key: &PropertyKey) -> bool {
    matches!(key.as_str(), "occupancy" | "b_factor")
}

fn assign_value(
    column: &mut PropertyColumn,
    key: &PropertyKey,
    index: usize,
    value: Option<PropertyValue>,
) -> Result<(), PropertyError> {
    match (column, value) {
        (PropertyColumn::Bool(values), None) => values[index] = None,
        (PropertyColumn::Bool(values), Some(PropertyValue::Bool(value))) => {
            values[index] = Some(value)
        }
        (PropertyColumn::Int(values), None) => values[index] = None,
        (PropertyColumn::Int(values), Some(PropertyValue::Int(value))) => {
            values[index] = Some(value)
        }
        (PropertyColumn::String(values), None) => values[index] = None,
        (PropertyColumn::String(values), Some(PropertyValue::String(value))) => {
            values[index] = Some(value)
        }
        (PropertyColumn::Real { values, .. }, None) => values[index] = None,
        (
            PropertyColumn::Real { unit, values },
            Some(PropertyValue::Real {
                value,
                unit: source_unit,
            }),
        ) => {
            let value = value * source_unit.conversion_factor_to(*unit)?;
            if !value.is_finite() {
                return Err(PropertyError::NonFiniteValue { index: Some(index) });
            }
            values[index] = Some(value);
        }
        _ => return Err(PropertyError::TypeMismatch { key: key.clone() }),
    }
    Ok(())
}

/// Complete generic property namespace for an owning domain object.
#[derive(Debug, Clone, PartialEq)]
pub struct Properties {
    owner: BTreeMap<PropertyKey, PropertyValue>,
    molecule_instances: PropertyTable,
    atoms: PropertyTable,
    bonds: PropertyTable,
    chains: PropertyTable,
    residues: PropertyTable,
    atom_sites: PropertyTable,
}

impl Default for Properties {
    fn default() -> Self {
        Self::new()
    }
}

impl Properties {
    pub fn new() -> Self {
        Self::with_dimensions(0, 0, 0, 0, 0, 0)
    }

    pub(crate) fn molecule(atom_slots: usize, bond_slots: usize) -> Self {
        Self::with_dimensions(0, atom_slots, bond_slots, 0, 0, 0)
    }

    pub fn realization(atom_count: usize, bond_count: usize) -> Self {
        Self::molecule(atom_count, bond_count)
    }

    fn with_dimensions(
        instance_count: usize,
        atom_count: usize,
        bond_count: usize,
        chain_count: usize,
        residue_count: usize,
        atom_site_count: usize,
    ) -> Self {
        Self {
            owner: BTreeMap::new(),
            molecule_instances: PropertyTable::new(instance_count),
            atoms: PropertyTable::new(atom_count),
            bonds: PropertyTable::new(bond_count),
            chains: PropertyTable::new(chain_count),
            residues: PropertyTable::new(residue_count),
            atom_sites: PropertyTable::new(atom_site_count),
        }
    }

    pub fn get(&self, key: &PropertyKey) -> Option<&PropertyValue> {
        self.owner.get(key)
    }

    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyKey, &PropertyValue)> + DoubleEndedIterator {
        self.owner.iter()
    }

    pub fn insert(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, PropertyError> {
        value.validate()?;
        Ok(self.owner.insert(key, value))
    }

    pub fn remove(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.owner.remove(key)
    }

    pub fn clear_owner(&mut self) {
        self.owner.clear();
    }

    pub fn owner_is_empty(&self) -> bool {
        self.owner.is_empty()
    }

    pub fn is_empty(&self) -> bool {
        self.owner.is_empty()
            && !self.molecule_instances.has_data()
            && !self.atoms.has_data()
            && !self.bonds.has_data()
            && !self.chains.has_data()
            && !self.residues.has_data()
            && !self.atom_sites.has_data()
    }

    pub(crate) const fn molecule_instances(&self) -> &PropertyTable {
        &self.molecule_instances
    }

    pub(crate) fn molecule_instances_mut(&mut self) -> &mut PropertyTable {
        &mut self.molecule_instances
    }

    pub(crate) const fn atoms(&self) -> &PropertyTable {
        &self.atoms
    }

    pub(crate) fn atoms_mut(&mut self) -> &mut PropertyTable {
        &mut self.atoms
    }

    pub(crate) const fn bonds(&self) -> &PropertyTable {
        &self.bonds
    }

    pub(crate) fn bonds_mut(&mut self) -> &mut PropertyTable {
        &mut self.bonds
    }

    pub(crate) const fn chains(&self) -> &PropertyTable {
        &self.chains
    }

    pub(crate) fn chains_mut(&mut self) -> &mut PropertyTable {
        &mut self.chains
    }

    pub(crate) const fn residues(&self) -> &PropertyTable {
        &self.residues
    }

    pub(crate) fn residues_mut(&mut self) -> &mut PropertyTable {
        &mut self.residues
    }

    pub(crate) const fn atom_sites(&self) -> &PropertyTable {
        &self.atom_sites
    }

    pub(crate) fn atom_sites_mut(&mut self) -> &mut PropertyTable {
        &mut self.atom_sites
    }

    /// Reads the dense atom table of a realization owner.
    pub const fn realization_atom_properties(&self) -> &PropertyTable {
        &self.atoms
    }

    /// Reads the dense bond table of a realization owner.
    pub const fn realization_bond_properties(&self) -> &PropertyTable {
        &self.bonds
    }

    /// Sets a non-canonical realization atom property.
    pub fn set_realization_atom_value(
        &mut self,
        key: PropertyKey,
        index: usize,
        value: Option<PropertyValue>,
    ) -> Result<(), PropertyError> {
        reject_reserved_realization_atom_key(&key)?;
        self.atoms.set_value(key, index, value)
    }

    /// Inserts or replaces a non-canonical realization atom column.
    pub fn insert_realization_atom_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, PropertyError> {
        reject_reserved_realization_atom_key(&key)?;
        self.atoms.insert(key, column)
    }

    /// Removes a non-canonical realization atom column.
    pub fn remove_realization_atom_column(
        &mut self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, PropertyError> {
        reject_reserved_realization_atom_key(key)?;
        Ok(self.atoms.remove(key))
    }

    pub fn set_realization_bond_value(
        &mut self,
        key: PropertyKey,
        index: usize,
        value: Option<PropertyValue>,
    ) -> Result<(), PropertyError> {
        self.bonds.set_value(key, index, value)
    }

    pub fn insert_realization_bond_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, PropertyError> {
        self.bonds.insert(key, column)
    }

    pub fn remove_realization_bond_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.bonds.remove(key)
    }

    pub fn occupancy_at(&self, index: usize) -> Result<Option<f64>, PropertyError> {
        match self.atoms.value(&occupancy_key(), index)? {
            None => Ok(None),
            Some(PropertyValue::Real { value, unit }) if unit == DIMENSIONLESS => Ok(Some(value)),
            _ => Err(PropertyError::InvalidCanonicalProperty(occupancy_key())),
        }
    }

    pub fn set_occupancy_at(
        &mut self,
        index: usize,
        value: Option<f64>,
    ) -> Result<(), PropertyError> {
        let value = value
            .map(|value| PropertyValue::real(value, DIMENSIONLESS))
            .transpose()?;
        self.atoms.set_value(occupancy_key(), index, value)
    }

    pub fn b_factor_at(&self, index: usize) -> Result<Option<Quantity<f64>>, PropertyError> {
        match self.atoms.value(&b_factor_key(), index)? {
            None => Ok(None),
            Some(PropertyValue::Real { value, unit }) if unit == SQUARE_NANOMETER => {
                Ok(Some(Quantity::new(value, unit)))
            }
            _ => Err(PropertyError::InvalidCanonicalProperty(b_factor_key())),
        }
    }

    pub fn set_b_factor_at(
        &mut self,
        index: usize,
        value: Option<Quantity<f64>>,
    ) -> Result<(), PropertyError> {
        let value = value
            .map(|value| {
                let value = value.to_unit(SQUARE_NANOMETER)?.to_value();
                PropertyValue::real(value, SQUARE_NANOMETER)
            })
            .transpose()?;
        self.atoms.set_value(b_factor_key(), index, value)
    }

    /// Validates the canonical realization-level atom columns, if present.
    pub fn validate_realization_canonical_properties(&self) -> Result<(), PropertyError> {
        if let Some(column) = self.atoms.get(&occupancy_key()) {
            if !matches!(column, PropertyColumn::Real { unit, .. } if *unit == DIMENSIONLESS) {
                return Err(PropertyError::InvalidCanonicalProperty(occupancy_key()));
            }
        }
        if let Some(column) = self.atoms.get(&b_factor_key()) {
            if !matches!(column, PropertyColumn::Real { unit, .. } if *unit == SQUARE_NANOMETER) {
                return Err(PropertyError::InvalidCanonicalProperty(b_factor_key()));
            }
        }
        Ok(())
    }

    /// Projects dense realization entity columns and drops owner-level values.
    pub fn project_realization(
        &self,
        atom_indices: &[usize],
        bond_indices: &[usize],
    ) -> Result<Self, PropertyError> {
        let properties = Self {
            owner: BTreeMap::new(),
            molecule_instances: PropertyTable::new(0),
            atoms: self.atoms.select_indices(atom_indices)?,
            bonds: self.bonds.select_indices(bond_indices)?,
            chains: PropertyTable::new(0),
            residues: PropertyTable::new(0),
            atom_sites: PropertyTable::new(0),
        };
        properties.validate_realization_canonical_properties()?;
        Ok(properties)
    }

    pub(crate) fn resize_atoms(&mut self, len: usize) {
        self.atoms.resize_missing(len);
    }

    pub(crate) fn resize_bonds(&mut self, len: usize) {
        self.bonds.resize_missing(len);
    }

    pub(crate) fn resize_domains(
        &mut self,
        instance_count: usize,
        atom_count: usize,
        bond_count: usize,
        chain_count: usize,
        residue_count: usize,
        atom_site_count: usize,
    ) {
        self.molecule_instances.resize_missing(instance_count);
        self.atoms.resize_missing(atom_count);
        self.bonds.resize_missing(bond_count);
        self.chains.resize_missing(chain_count);
        self.residues.resize_missing(residue_count);
        self.atom_sites.resize_missing(atom_site_count);
    }

    pub(crate) fn project_topology(
        &self,
        molecule_instances: &[Option<usize>],
        atoms: &[usize],
        bonds: &[usize],
        chains: &[usize],
        residues: &[usize],
        atom_sites: &[usize],
    ) -> Result<Self, PropertyError> {
        Ok(Self {
            owner: BTreeMap::new(),
            molecule_instances: self
                .molecule_instances
                .select_optional_indices(molecule_instances)?,
            atoms: self.atoms.select_indices(atoms)?,
            bonds: self.bonds.select_indices(bonds)?,
            chains: self.chains.select_indices(chains)?,
            residues: self.residues.select_indices(residues)?,
            atom_sites: self.atom_sites.select_indices(atom_sites)?,
        })
    }
}

fn occupancy_key() -> PropertyKey {
    PropertyKey::new("occupancy").expect("canonical property key is valid")
}

fn b_factor_key() -> PropertyKey {
    PropertyKey::new("b_factor").expect("canonical property key is valid")
}

fn reject_reserved_realization_atom_key(key: &PropertyKey) -> Result<(), PropertyError> {
    if is_reserved_realization_atom_key(key) {
        return Err(PropertyError::ReservedKey(key.clone()));
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PropertyError {
    InvalidKey(String),
    LengthMismatch { expected: usize, actual: usize },
    InvalidIndex { len: usize, index: usize },
    TypeMismatch { key: PropertyKey },
    ReservedKey(PropertyKey),
    InvalidCanonicalProperty(PropertyKey),
    NonFiniteValue { index: Option<usize> },
    Unit(UnitError),
}

impl fmt::Display for PropertyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKey(key) => write!(
                formatter,
                "invalid property key {key:?}; use a 1-{MAX_PROPERTY_KEY_LEN} character ASCII identifier"
            ),
            Self::LengthMismatch { expected, actual } => {
                write!(formatter, "property column requires {expected} values, received {actual}")
            }
            Self::InvalidIndex { len, index } => {
                write!(formatter, "property index {index} is outside table length {len}")
            }
            Self::TypeMismatch { key } => write!(formatter, "property {key:?} has a different type"),
            Self::ReservedKey(key) => write!(formatter, "property key {key:?} is reserved for a canonical semantic API"),
            Self::InvalidCanonicalProperty(key) => write!(
                formatter,
                "property {key:?} does not have its canonical realization atom type and unit"
            ),
            Self::NonFiniteValue { index: Some(index) } => {
                write!(formatter, "real property value at index {index} must be finite")
            }
            Self::NonFiniteValue { index: None } => {
                formatter.write_str("real property value must be finite")
            }
            Self::Unit(error) => write!(formatter, "invalid property unit: {error}"),
        }
    }
}

impl std::error::Error for PropertyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

impl From<UnitError> for PropertyError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{
        Quantity, ANGSTROM, DIMENSIONLESS, KELVIN, NANOMETER, SQUARE_ANGSTROM, SQUARE_NANOMETER,
    };

    fn key(value: &str) -> PropertyKey {
        PropertyKey::new(value).unwrap()
    }

    #[test]
    fn property_keys_use_one_conservative_validation_path() {
        for valid in ["x", "_private", "partial_charge", "force-field.v1", "a-2"] {
            assert_eq!(PropertyKey::new(valid).unwrap().as_str(), valid);
            assert_eq!(valid.parse::<PropertyKey>().unwrap().as_str(), valid);
        }
        for invalid in ["", "2bad", "has space", "slash/name", "unicode_µ"] {
            assert!(PropertyKey::new(invalid).is_err(), "accepted {invalid:?}");
        }
        assert!(PropertyKey::new("x".repeat(MAX_PROPERTY_KEY_LEN + 1)).is_err());
    }

    #[test]
    fn scalar_values_cover_the_complete_domain_and_reject_non_finite_reals() {
        let mut properties = Properties::new();
        properties
            .insert(key("bool"), PropertyValue::Bool(true))
            .unwrap();
        properties
            .insert(key("int"), PropertyValue::Int(-7))
            .unwrap();
        properties
            .insert(
                key("real"),
                PropertyValue::Real {
                    value: 2.5,
                    unit: DIMENSIONLESS,
                },
            )
            .unwrap();
        properties
            .insert(key("string"), PropertyValue::String("value".into()))
            .unwrap();
        assert_eq!(properties.iter().count(), 4);
        assert_eq!(
            properties
                .iter()
                .map(|(key, _)| key.as_str())
                .collect::<Vec<_>>(),
            ["bool", "int", "real", "string"]
        );
        assert!(matches!(
            properties.insert(
                key("nan"),
                PropertyValue::Real {
                    value: f64::NAN,
                    unit: DIMENSIONLESS,
                },
            ),
            Err(PropertyError::NonFiniteValue { index: None })
        ));
    }

    #[test]
    fn every_column_type_supports_missing_values_and_deterministic_iteration() {
        let mut table = PropertyTable::new(3);
        table
            .insert(
                key("z_string"),
                PropertyColumn::String(vec![None, Some("x".into()), None]),
            )
            .unwrap();
        table
            .insert(
                key("a_bool"),
                PropertyColumn::Bool(vec![Some(true), None, Some(false)]),
            )
            .unwrap();
        table
            .insert(key("m_int"), PropertyColumn::Int(vec![None, Some(3), None]))
            .unwrap();
        table
            .insert(
                key("r_real"),
                PropertyColumn::Real {
                    unit: DIMENSIONLESS,
                    values: vec![Some(1.0), None, Some(2.0)],
                },
            )
            .unwrap();
        assert_eq!(
            table
                .iter()
                .map(|(key, _)| key.as_str())
                .collect::<Vec<_>>(),
            ["a_bool", "m_int", "r_real", "z_string"]
        );
        assert_eq!(table.value(&key("m_int"), 0).unwrap(), None);
        assert_eq!(
            table.value(&key("m_int"), 1).unwrap(),
            Some(PropertyValue::Int(3))
        );
    }

    #[test]
    fn table_updates_are_transactional_unit_aware_and_normalize_all_missing() {
        let mut table = PropertyTable::new(2);
        let length_error = table.insert(key("bad_len"), PropertyColumn::Bool(vec![Some(true)]));
        assert!(matches!(
            length_error,
            Err(PropertyError::LengthMismatch { .. })
        ));
        assert!(!table.has_data());

        table
            .insert(
                key("distance"),
                PropertyColumn::Real {
                    unit: NANOMETER,
                    values: vec![Some(1.0), None],
                },
            )
            .unwrap();
        table
            .insert(
                key("distance"),
                PropertyColumn::Real {
                    unit: ANGSTROM,
                    values: vec![Some(5.0), Some(10.0)],
                },
            )
            .unwrap();
        let Some(PropertyValue::Real { value, unit }) = table.value(&key("distance"), 0).unwrap()
        else {
            panic!("distance should be a real value");
        };
        assert_eq!(unit, NANOMETER);
        assert!((value - 0.5).abs() < 1.0e-12);

        let before = table.clone();
        assert!(matches!(
            table.insert(
                key("distance"),
                PropertyColumn::Real {
                    unit: KELVIN,
                    values: vec![Some(1.0), Some(2.0)],
                },
            ),
            Err(PropertyError::Unit(_))
        ));
        assert_eq!(table, before);
        assert!(matches!(
            table.insert(
                key("finite"),
                PropertyColumn::Real {
                    unit: DIMENSIONLESS,
                    values: vec![Some(f64::INFINITY), None],
                },
            ),
            Err(PropertyError::NonFiniteValue { index: Some(0) })
        ));

        table.clear_value(key("distance"), 0).unwrap();
        table.clear_value(key("distance"), 1).unwrap();
        assert!(table.get(&key("distance")).is_none());
        table
            .insert(key("missing"), PropertyColumn::Int(vec![None, None]))
            .unwrap();
        assert!(table.get(&key("missing")).is_none());
    }

    #[test]
    fn checked_projection_preserves_columns_and_missing_cells() {
        let mut table = PropertyTable::new(3);
        table
            .insert(
                key("flag"),
                PropertyColumn::Bool(vec![Some(true), None, Some(false)]),
            )
            .unwrap();
        let selected = table.select_indices(&[2, 1]).unwrap();
        assert_eq!(selected.len(), 2);
        assert_eq!(
            selected.value(&key("flag"), 0).unwrap(),
            Some(PropertyValue::Bool(false))
        );
        assert_eq!(selected.value(&key("flag"), 1).unwrap(), None);
        assert!(matches!(
            table.select_indices(&[3]),
            Err(PropertyError::InvalidIndex { .. })
        ));
    }

    #[test]
    fn generic_tables_do_not_reserve_realization_semantic_names() {
        let mut table = PropertyTable::new(1);
        table
            .set_value(key("occupancy"), 0, Some(PropertyValue::Int(7)))
            .unwrap();
        table
            .insert(
                key("b_factor"),
                PropertyColumn::String(vec![Some("generic".into())]),
            )
            .unwrap();
        assert_eq!(
            table.value(&key("occupancy"), 0).unwrap(),
            Some(PropertyValue::Int(7))
        );
        assert_eq!(
            table.value(&key("b_factor"), 0).unwrap(),
            Some(PropertyValue::String("generic".into()))
        );
    }

    #[test]
    fn realization_facade_reserves_and_validates_canonical_atom_properties() {
        let mut properties = Properties::realization(2, 1);
        for reserved in ["occupancy", "b_factor"] {
            let key = key(reserved);
            assert!(matches!(
                properties.set_realization_atom_value(key.clone(), 0, Some(PropertyValue::Int(1))),
                Err(PropertyError::ReservedKey(_))
            ));
            assert!(matches!(
                properties.insert_realization_atom_column(
                    key.clone(),
                    PropertyColumn::Int(vec![Some(1), None])
                ),
                Err(PropertyError::ReservedKey(_))
            ));
            assert!(matches!(
                properties.remove_realization_atom_column(&key),
                Err(PropertyError::ReservedKey(_))
            ));
        }

        properties.set_occupancy_at(0, Some(0.75)).unwrap();
        properties
            .set_b_factor_at(1, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
            .unwrap();
        assert_eq!(properties.occupancy_at(0).unwrap(), Some(0.75));
        let b_factor = properties.b_factor_at(1).unwrap().unwrap();
        assert_eq!(b_factor.unit(), SQUARE_NANOMETER);
        assert!((*b_factor.value() - 0.125).abs() < 1.0e-12);
        assert!(matches!(
            properties.set_occupancy_at(0, Some(f64::NAN)),
            Err(PropertyError::NonFiniteValue { .. })
        ));
        assert!(matches!(
            properties.set_b_factor_at(0, Some(Quantity::new(1.0, KELVIN))),
            Err(PropertyError::Unit(_))
        ));

        let mut malformed = Properties::realization(2, 1);
        malformed
            .atoms_mut()
            .insert(
                key("occupancy"),
                PropertyColumn::Real {
                    unit: KELVIN,
                    values: vec![Some(1.0), None],
                },
            )
            .unwrap();
        assert!(matches!(
            malformed.validate_realization_canonical_properties(),
            Err(PropertyError::InvalidCanonicalProperty(_))
        ));
    }
}
