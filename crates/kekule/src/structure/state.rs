//! Validation and accessors shared by the flat structure owners: models,
//! ensemble members, trajectory frames, and model drafts.
//!
//! Each owner stores its own `positions`, `cell`, and `properties` fields, and
//! models and ensemble members also `occupancies` and `b_factors`. The macros
//! below generate the
//! common public accessors from those fields; an owner provides
//! `storage_ref(&self)` and `storage(&mut self)` returning the value that holds
//! them.

use std::fmt;

use crate::geometry::{PeriodicCell, PeriodicCellError, RigidTransform};
use crate::properties::{PropertyError, StructureProperties};
use crate::topology::transform::TopologySubset;
use crate::topology::{Topology, TopologyAtomIndex};
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};

use super::{PositionError, Positions};

/// The row of one dense atom index in an owner with `len` atoms.
pub(crate) fn atom_row(len: usize, atom: TopologyAtomIndex) -> Result<usize, StructureError> {
    let index = atom.index();
    if index >= len {
        return Err(StructureError::InvalidAtomIndex { index, len });
    }
    Ok(index)
}

/// Checks that structure rows can be bound to `topology` without changing
/// anything: one position per atom, and bond rows either matching or not yet
/// allocated.
pub(crate) fn validate_binding(
    atoms: usize,
    properties: &StructureProperties,
    topology: &Topology,
) -> Result<(), StructureError> {
    if atoms != topology.atom_count() {
        return Err(StructureError::AtomCountMismatch {
            expected: topology.atom_count(),
            actual: atoms,
        });
    }
    properties
        .validate_bond_rows(topology.bond_count())
        .map_err(|_| StructureError::BondCountMismatch {
            expected: topology.bond_count(),
            actual: properties.bonds().len(),
        })
}

/// Whether structure rows already cover every atom and bond of `topology`.
pub(crate) fn is_bound(
    atoms: usize,
    properties: &StructureProperties,
    topology: &Topology,
) -> bool {
    atoms == topology.atom_count() && properties.bonds().len() == topology.bond_count()
}

/// Checks that a replacement keeps the current atom and bond row counts.
pub(crate) fn check_property_rows(
    atoms: usize,
    current: &StructureProperties,
    replacement: &StructureProperties,
) -> Result<(), StructureError> {
    if replacement.atoms().len() != atoms {
        return Err(StructureError::AtomCountMismatch {
            expected: atoms,
            actual: replacement.atoms().len(),
        });
    }
    if replacement.bonds().len() != current.bonds().len() {
        return Err(StructureError::BondCountMismatch {
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
) -> Result<Positions, StructureError> {
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
) -> Result<Option<PeriodicCell>, StructureError> {
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

        pub fn cell(&self) -> Option<&$crate::geometry::PeriodicCell> {
            self.storage_ref().cell.as_ref()
        }

        /// Structure annotations: owner values plus one row per dense atom
        /// and bond.
        pub fn properties(&self) -> &$crate::properties::StructureProperties {
            &self.storage_ref().properties
        }
    };
}
pub(crate) use structure_readers;

/// Occupancy and B-factor readers of models and ensemble members.
macro_rules! site_readers {
    () => {
        /// Occupancy of one atom; `None` when it was never recorded.
        pub fn occupancy(
            &self,
            atom: $crate::topology::TopologyAtomIndex,
        ) -> Result<Option<f64>, $crate::structure::StructureError> {
            let storage = self.storage_ref();
            storage.occupancies.get(atom, storage.positions.len())
        }

        /// Every occupancy in dense atom order, or `None` when none was
        /// recorded.
        pub fn occupancies(&self) -> Option<&[Option<f64>]> {
            self.storage_ref().occupancies.values()
        }

        /// Isotropic B-factor of one atom; `None` when it was never recorded.
        pub fn b_factor(
            &self,
            atom: $crate::topology::TopologyAtomIndex,
        ) -> Result<Option<$crate::units::Quantity<f64>>, $crate::structure::StructureError> {
            let storage = self.storage_ref();
            Ok($crate::structure::scalars::b_factor_quantity(
                storage.b_factors.get(atom, storage.positions.len())?,
            ))
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
        pub fn set_positions<T: AsRef<[$crate::geometry::Point3]>>(
            &mut self,
            positions: $crate::units::Quantity<T>,
        ) -> Result<(), $crate::structure::StructureError> {
            Ok(self.storage().positions.set_all(positions)?)
        }

        /// Sets one position by dense atom index.
        pub fn set_position_at(
            &mut self,
            atom: $crate::topology::TopologyAtomIndex,
            position: $crate::units::Quantity<$crate::geometry::Point3>,
        ) -> Result<(), $crate::structure::StructureError> {
            Ok(self
                .storage()
                .positions
                .set_position_at(atom.index(), position)?)
        }

        pub fn set_cell(&mut self, cell: Option<$crate::geometry::PeriodicCell>) {
            self.storage().cell = cell;
        }

        /// Replaces every structure annotation with rows of the same
        /// dimensions, for example ones copied from another structure of the
        /// same topology. Detached payloads accept no bond rows.
        pub fn set_properties(
            &mut self,
            properties: $crate::properties::StructureProperties,
        ) -> Result<(), $crate::structure::StructureError> {
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
        pub fn properties_mut(&mut self) -> $crate::properties::StructurePropertiesMut<'_> {
            $crate::properties::StructurePropertiesMut::new(&mut self.storage().properties)
        }
    };
}
pub(crate) use structure_mutators;

/// Occupancy and B-factor mutators of models and ensemble members.
macro_rules! site_mutators {
    () => {
        /// Sets or clears one occupancy. Values must be finite.
        pub fn set_occupancy(
            &mut self,
            atom: $crate::topology::TopologyAtomIndex,
            value: Option<f64>,
        ) -> Result<(), $crate::structure::StructureError> {
            let storage = self.storage();
            let len = storage.positions.len();
            storage.occupancies.set(atom, len, value)
        }

        /// Replaces every occupancy in dense atom order, or clears them all.
        pub fn set_occupancies(
            &mut self,
            values: Option<Vec<Option<f64>>>,
        ) -> Result<(), $crate::structure::StructureError> {
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
        ) -> Result<(), $crate::structure::StructureError> {
            let value = $crate::structure::scalars::canonical_b_factor(value)?;
            let storage = self.storage();
            let len = storage.positions.len();
            storage.b_factors.set(atom, len, value)
        }

        /// Replaces every B-factor in dense atom order, or clears them all.
        pub fn set_b_factors(
            &mut self,
            values: Option<$crate::units::Quantity<Vec<Option<f64>>>>,
        ) -> Result<(), $crate::structure::StructureError> {
            let values = $crate::structure::scalars::canonical_b_factors(values)?;
            let storage = self.storage();
            let len = storage.positions.len();
            storage.b_factors.replace(len, values)
        }
    };
}
pub(crate) use site_mutators;

/// Invalid structure payload state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum StructureError {
    /// Dense atom state has the wrong number of rows.
    AtomCountMismatch {
        expected: usize,
        actual: usize,
    },
    /// Bound bond property rows have the wrong length.
    BondCountMismatch {
        expected: usize,
        actual: usize,
    },
    /// A dense atom index is outside the structure.
    InvalidAtomIndex {
        index: usize,
        len: usize,
    },
    /// A scalar atom value is not finite.
    NonFiniteValue {
        index: usize,
    },
    /// A dense vector value is not finite.
    NonFiniteVector {
        index: usize,
    },
    /// A frame time is not finite.
    NonFiniteTime,
    /// An ensemble weight is not finite and positive.
    InvalidWeight,
    /// The payload exceeds addressable capacity.
    CapacityOverflow,
    Position(PositionError),
    Property(Box<PropertyError>),
    Cell(Box<PeriodicCellError>),
    Unit(Box<UnitError>),
}

impl fmt::Display for StructureError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AtomCountMismatch { expected, actual } => write!(
                formatter,
                "structure requires {expected} atom rows, but received {actual}"
            ),
            Self::BondCountMismatch { expected, actual } => write!(
                formatter,
                "structure requires {expected} bond property rows, but has {actual}"
            ),
            Self::InvalidAtomIndex { index, len } => write!(
                formatter,
                "atom index {index} is outside a structure of {len} atoms"
            ),
            Self::NonFiniteValue { index } => {
                write!(formatter, "atom value at index {index} must be finite")
            }
            Self::NonFiniteVector { index } => {
                write!(formatter, "vector at index {index} must be finite")
            }
            Self::NonFiniteTime => formatter.write_str("frame time must be finite"),
            Self::InvalidWeight => {
                formatter.write_str("ensemble weight must be finite and positive")
            }
            Self::CapacityOverflow => formatter.write_str("structure exceeds addressable capacity"),
            Self::Position(error) => write!(formatter, "invalid positions: {error}"),
            Self::Property(error) => write!(formatter, "invalid structure property: {error}"),
            Self::Cell(error) => write!(formatter, "invalid periodic cell: {error}"),
            Self::Unit(error) => write!(formatter, "invalid quantity unit: {error}"),
        }
    }
}

impl std::error::Error for StructureError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Position(error) => Some(error),
            Self::Property(error) => Some(error.as_ref()),
            Self::Cell(error) => Some(error.as_ref()),
            Self::Unit(error) => Some(error.as_ref()),
            _ => None,
        }
    }
}

impl From<PositionError> for StructureError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<PropertyError> for StructureError {
    fn from(error: PropertyError) -> Self {
        Self::Property(Box::new(error))
    }
}

impl From<PeriodicCellError> for StructureError {
    fn from(error: PeriodicCellError) -> Self {
        Self::Cell(Box::new(error))
    }
}

impl From<UnitError> for StructureError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
