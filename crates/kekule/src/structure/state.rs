//! Validation and accessors shared by the flat structure owners: structures,
//! geometries, trajectory frames, and structure drafts.
//!
//! Each owner stores its own `positions`, `cell`, and `properties` fields, and
//! structures and geometries also `occupancies` and `b_factors`. The macros
//! below generate the
//! common public accessors from those fields; an owner provides
//! `storage_ref(&self)` and `storage(&mut self)` returning the value that holds
//! them.

use std::fmt;

use crate::math::{PeriodicCell, PeriodicCellError, RigidTransform};
use crate::properties::{GeometryProperties, PropertyError};
use crate::topology::TopologySubset;
use crate::topology::{Topology, TopologyAtomIndex};
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};

use super::{DenseArrayError, Positions};

/// The row of one dense atom index in an owner with `len` atoms.
pub(crate) fn atom_row(len: usize, atom: TopologyAtomIndex) -> Result<usize, GeometryError> {
    let index = atom.index();
    if index >= len {
        return Err(GeometryError::Array(DenseArrayError::IndexOutOfRange {
            index,
            len,
        }));
    }
    Ok(index)
}

/// Checks that structure rows can be bound to `topology` without changing
/// anything: one position per atom, and bond rows either matching or not yet
/// allocated.
pub(crate) fn validate_binding(
    atoms: usize,
    properties: &GeometryProperties,
    topology: &Topology,
) -> Result<(), GeometryError> {
    if atoms != topology.atom_count() {
        return Err(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: topology.atom_count(),
            actual: atoms,
        }));
    }
    properties
        .validate_bond_rows(topology.bond_count())
        .map_err(|_| GeometryError::BondCountMismatch {
            expected: topology.bond_count(),
            actual: properties.bonds().len(),
        })
}

/// Whether structure rows already cover every atom and bond of `topology`.
pub(crate) fn is_bound(atoms: usize, properties: &GeometryProperties, topology: &Topology) -> bool {
    atoms == topology.atom_count() && properties.bonds().len() == topology.bond_count()
}

/// Checks that a replacement keeps the current atom and bond row counts.
pub(crate) fn check_property_rows(
    atoms: usize,
    current: &GeometryProperties,
    replacement: &GeometryProperties,
) -> Result<(), GeometryError> {
    if replacement.atoms().len() != atoms {
        return Err(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: atoms,
            actual: replacement.atoms().len(),
        }));
    }
    if replacement.bonds().len() != current.bonds().len() {
        return Err(GeometryError::BondCountMismatch {
            expected: current.bonds().len(),
            actual: replacement.bonds().len(),
        });
    }
    Ok(())
}

/// Source dense atom and bond rows of a subset, in subset dense order.
pub(crate) fn subset_rows(subset: &TopologySubset) -> (Vec<usize>, Vec<usize>) {
    let correspondence = subset.correspondence();
    let atoms = correspondence
        .source_atom_indices()
        .iter()
        .map(|index| index.index())
        .collect();
    let bonds = correspondence
        .source_bond_indices()
        .iter()
        .map(|index| index.index())
        .collect();
    (atoms, bonds)
}

/// Applies a rigid transform to every position.
pub(crate) fn transform_positions(
    positions: &Positions,
    transform: RigidTransform,
) -> Result<Positions, GeometryError> {
    let points = positions
        .values()
        .value()
        .iter()
        .map(|point| transform.transform_point(*point))
        .collect::<Vec<_>>();
    Ok(Positions::new(Quantity::new(
        points,
        CANONICAL_LENGTH_UNIT,
    ))?)
}

/// Rotates cell vectors; a rigid transform's translation does not move a cell.
pub(crate) fn transform_cell(
    cell: Option<PeriodicCell>,
    transform: RigidTransform,
) -> Result<Option<PeriodicCell>, GeometryError> {
    cell.map(|cell| {
        let vectors = cell
            .vectors()
            .map(|vectors| vectors.map(|vector| transform.transform_vector(vector)));
        Ok(PeriodicCell::new(vectors, cell.periodic_axes())?)
    })
    .transpose()
}

/// Positions, cell, and structure-property readers.
macro_rules! structure_readers {
    () => {
        /// Positions in dense topology atom order.
        pub fn positions(&self) -> &$crate::structure::Positions {
            &self.storage_ref().positions
        }

        /// The periodic cell, or `None` for a non-periodic structure.
        pub fn cell(&self) -> Option<&$crate::math::PeriodicCell> {
            self.storage_ref().cell.as_ref()
        }

        /// Structure annotations: owner values plus one row per dense atom
        /// and bond.
        pub fn properties(&self) -> &$crate::properties::GeometryProperties {
            &self.storage_ref().properties
        }
    };
}
pub(crate) use structure_readers;

/// Occupancy and B-factor readers of structures and geometries.
macro_rules! site_readers {
    () => {
        /// Occupancy of one atom; `None` when it was never recorded or `atom`
        /// is out of range.
        pub fn occupancy(&self, atom: $crate::topology::TopologyAtomIndex) -> Option<f64> {
            self.storage_ref().occupancies.get(atom)
        }

        /// Every occupancy in dense atom order, or `None` when none was
        /// recorded.
        pub fn occupancies(&self) -> Option<&[Option<f64>]> {
            self.storage_ref().occupancies.values()
        }

        /// Isotropic B-factor of one atom; `None` when it was never recorded or
        /// `atom` is out of range.
        pub fn b_factor(
            &self,
            atom: $crate::topology::TopologyAtomIndex,
        ) -> Option<$crate::units::Quantity<f64>> {
            $crate::structure::scalars::b_factor_quantity(self.storage_ref().b_factors.get(atom))
        }

        /// Every B-factor in dense atom order and square nanometres, or `None`
        /// when none was recorded.
        pub fn b_factors(&self) -> Option<$crate::units::Quantity<&[Option<f64>]>> {
            $crate::structure::scalars::b_factor_quantities(&self.storage_ref().b_factors)
        }
    };
}
pub(crate) use site_readers;

/// Dimension-preserving position, cell, and structure-property mutators.
macro_rules! structure_mutators {
    () => {
        /// Replaces every position without changing the atom count. Units,
        /// dimensions, and finite values are checked before mutation.
        pub fn set_positions<T: AsRef<[$crate::math::Point3]>>(
            &mut self,
            positions: $crate::units::Quantity<T>,
        ) -> Result<(), $crate::structure::GeometryError> {
            Ok(self.storage().positions.set_all(positions)?)
        }

        /// Sets one position by dense atom index.
        pub fn set_position_at(
            &mut self,
            atom: $crate::topology::TopologyAtomIndex,
            position: $crate::units::Quantity<$crate::math::Point3>,
        ) -> Result<(), $crate::structure::GeometryError> {
            Ok(self.storage().positions.set(atom.index(), position)?)
        }

        /// Replaces the periodic cell; `None` makes the structure non-periodic.
        /// Positions are not wrapped or rescaled.
        pub fn set_cell(&mut self, cell: Option<$crate::math::PeriodicCell>) {
            self.storage().cell = cell;
        }

        /// Replaces every structure annotation with rows of the same
        /// dimensions, for example ones copied from another structure of the
        /// same topology. Detached payloads accept no bond rows.
        pub fn set_properties(
            &mut self,
            properties: $crate::properties::GeometryProperties,
        ) -> Result<(), $crate::structure::GeometryError> {
            let storage = self.storage();
            $crate::structure::state::check_property_rows(
                storage.positions.len(),
                &storage.properties,
                &properties,
            )?;
            storage.properties = properties;
            Ok(())
        }

        /// Length-preserving mutable access to structure annotations.
        pub fn properties_mut(&mut self) -> $crate::properties::GeometryPropertiesMut<'_> {
            $crate::properties::GeometryPropertiesMut::new(&mut self.storage().properties)
        }
    };
}
pub(crate) use structure_mutators;

/// Occupancy and B-factor mutators of structures and geometries.
macro_rules! site_mutators {
    () => {
        /// Sets or clears one occupancy. Values must be finite.
        pub fn set_occupancy(
            &mut self,
            atom: $crate::topology::TopologyAtomIndex,
            value: Option<f64>,
        ) -> Result<(), $crate::structure::GeometryError> {
            let storage = self.storage();
            let len = storage.positions.len();
            storage.occupancies.set(atom, len, value)
        }

        /// Replaces every occupancy in dense atom order, or clears them all.
        pub fn set_occupancies(
            &mut self,
            values: Option<Vec<Option<f64>>>,
        ) -> Result<(), $crate::structure::GeometryError> {
            let storage = self.storage();
            let len = storage.positions.len();
            storage.occupancies.replace(len, values)
        }

        /// Sets or clears one isotropic B-factor, converting it to square
        /// nanometres.
        pub fn set_b_factor(
            &mut self,
            atom: $crate::topology::TopologyAtomIndex,
            value: Option<$crate::units::Quantity<f64>>,
        ) -> Result<(), $crate::structure::GeometryError> {
            let value = $crate::structure::scalars::canonical_b_factor(value)?;
            let storage = self.storage();
            let len = storage.positions.len();
            storage.b_factors.set(atom, len, value)
        }

        /// Replaces every B-factor in dense atom order, or clears them all.
        pub fn set_b_factors(
            &mut self,
            values: Option<$crate::units::Quantity<Vec<Option<f64>>>>,
        ) -> Result<(), $crate::structure::GeometryError> {
            let values = $crate::structure::scalars::canonical_b_factors(values)?;
            let storage = self.storage();
            let len = storage.positions.len();
            storage.b_factors.replace(len, values)
        }
    };
}
pub(crate) use site_mutators;

/// Invalid geometry state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum GeometryError {
    /// Dense per-atom state is invalid.
    Array(DenseArrayError),
    /// Bound bond property rows have the wrong length.
    BondCountMismatch {
        /// The required bond row count.
        expected: usize,
        /// The supplied bond row count.
        actual: usize,
    },
    /// A frame time is not finite.
    NonFiniteTime,
    /// A state energy is not finite.
    NonFiniteEnergy,
    /// The payload exceeds addressable capacity.
    CapacityOverflow,
    /// A structure property edit failed; see the source.
    Property(Box<PropertyError>),
    /// The periodic cell is invalid; see the source.
    Cell(Box<PeriodicCellError>),
    /// A quantity could not be converted to the required unit; see the source.
    Unit(Box<UnitError>),
}

impl fmt::Display for GeometryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Array(error) => error.fmt(formatter),
            Self::BondCountMismatch { expected, actual } => write!(
                formatter,
                "The bond properties have {actual} rows, but the topology has {expected} bonds."
            ),
            Self::NonFiniteTime => formatter.write_str("The frame time is not finite."),
            Self::NonFiniteEnergy => formatter.write_str("The state energy is not finite."),
            Self::CapacityOverflow => {
                formatter.write_str("The geometry exceeds the addressable capacity.")
            }
            Self::Property(error) => error.fmt(formatter),
            Self::Cell(error) => error.fmt(formatter),
            Self::Unit(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for GeometryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Array(error) => error.source(),
            Self::Property(error) => error.source(),
            Self::Cell(error) => error.source(),
            Self::Unit(error) => error.source(),
            _ => None,
        }
    }
}

impl From<DenseArrayError> for GeometryError {
    fn from(error: DenseArrayError) -> Self {
        Self::Array(error)
    }
}

impl From<PropertyError> for GeometryError {
    fn from(error: PropertyError) -> Self {
        Self::Property(Box::new(error))
    }
}

impl From<PeriodicCellError> for GeometryError {
    fn from(error: PeriodicCellError) -> Self {
        Self::Cell(Box::new(error))
    }
}

impl From<UnitError> for GeometryError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
