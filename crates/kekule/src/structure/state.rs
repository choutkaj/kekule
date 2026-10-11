//! Validation and accessors shared by the flat realization owners: models,
//! ensemble members, trajectory frames, and model drafts.
//!
//! Each owner stores its own `positions`, `cell`, and `properties` fields, and
//! models and ensemble members also `sites`. The macros below generate the
//! common public accessors from those fields; an owner provides
//! `storage_ref(&self)` and `storage(&mut self)` returning the value that holds
//! them.

use std::fmt;

use crate::geometry::{PeriodicCell, PeriodicCellError, RigidTransform};
use crate::properties::{PropertyError, RealizationProperties};
use crate::topology::{Topology, TopologyAtomIndex};
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};

use super::{PositionError, Positions};

/// The row of one dense atom index in an owner with `len` atoms.
pub(crate) fn atom_row(len: usize, atom: TopologyAtomIndex) -> Result<usize, ConformationError> {
    let index = atom.index();
    if index >= len {
        return Err(ConformationError::InvalidAtomIndex { index, len });
    }
    Ok(index)
}

/// Checks that realization rows can be bound to `topology` without changing
/// anything: one position per atom, and bond rows either matching or not yet
/// allocated.
pub(crate) fn validate_binding(
    atoms: usize,
    properties: &RealizationProperties,
    topology: &Topology,
) -> Result<(), ConformationError> {
    if atoms != topology.atom_count() {
        return Err(ConformationError::AtomCountMismatch {
            expected: topology.atom_count(),
            actual: atoms,
        });
    }
    properties
        .validate_bond_rows(topology.bond_count())
        .map_err(|_| ConformationError::BondCountMismatch {
            expected: topology.bond_count(),
            actual: properties.bonds().len(),
        })
}

/// Whether realization rows already cover every atom and bond of `topology`.
pub(crate) fn is_bound(
    atoms: usize,
    properties: &RealizationProperties,
    topology: &Topology,
) -> bool {
    atoms == topology.atom_count() && properties.bonds().len() == topology.bond_count()
}

/// Checks that a replacement keeps the current atom and bond row counts.
pub(crate) fn check_property_rows(
    atoms: usize,
    current: &RealizationProperties,
    replacement: &RealizationProperties,
) -> Result<(), ConformationError> {
    if replacement.atoms().len() != atoms {
        return Err(ConformationError::AtomCountMismatch {
            expected: atoms,
            actual: replacement.atoms().len(),
        });
    }
    if replacement.bonds().len() != current.bonds().len() {
        return Err(ConformationError::BondCountMismatch {
            expected: current.bonds().len(),
            actual: replacement.bonds().len(),
        });
    }
    Ok(())
}

/// Applies a rigid transform to every position.
pub(crate) fn transform_positions(
    positions: &Positions,
    transform: RigidTransform,
) -> Result<Positions, ConformationError> {
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
) -> Result<Option<PeriodicCell>, ConformationError> {
    cell.map(|cell| {
        let vectors = cell
            .vectors()
            .map(|vectors| vectors.map(|vector| transform.transform_vector(vector)));
        Ok(PeriodicCell::new(vectors, cell.periodic_axes())?)
    })
    .transpose()
}

/// Positions, cell, and realization-property readers.
macro_rules! realization_readers {
    () => {
        /// Positions in dense topology atom order.
        pub fn positions(&self) -> &$crate::structure::Positions {
            &self.storage_ref().positions
        }

        pub fn cell(&self) -> Option<&$crate::geometry::PeriodicCell> {
            self.storage_ref().cell.as_ref()
        }

        /// Realization annotations: owner values plus one row per dense atom
        /// and bond.
        pub fn properties(&self) -> &$crate::properties::RealizationProperties {
            &self.storage_ref().properties
        }
    };
}
pub(crate) use realization_readers;

/// Occupancy and B-factor readers of models and ensemble members.
macro_rules! site_readers {
    () => {
        /// Occupancy of one atom; `None` when it was never recorded.
        pub fn occupancy(
            &self,
            atom: $crate::topology::TopologyAtomIndex,
        ) -> Result<Option<f64>, $crate::structure::ConformationError> {
            let storage = self.storage_ref();
            storage.sites.occupancy(atom, storage.positions.len())
        }

        /// Every occupancy in dense atom order, or `None` when none was
        /// recorded.
        pub fn occupancies(&self) -> Option<&[Option<f64>]> {
            self.storage_ref().sites.occupancies()
        }

        /// Isotropic B-factor of one atom; `None` when it was never recorded.
        pub fn b_factor(
            &self,
            atom: $crate::topology::TopologyAtomIndex,
        ) -> Result<Option<$crate::units::Quantity<f64>>, $crate::structure::ConformationError> {
            let storage = self.storage_ref();
            storage.sites.b_factor(atom, storage.positions.len())
        }

        /// Every B-factor in dense atom order and square nanometres, or `None`
        /// when none was recorded.
        pub fn b_factors(&self) -> Option<$crate::units::Quantity<&[Option<f64>]>> {
            self.storage_ref().sites.b_factors()
        }
    };
}
pub(crate) use site_readers;

/// Dimension-preserving position, cell, and realization-property mutators.
macro_rules! realization_mutators {
    () => {
        /// Replaces every position without changing the atom count. Units,
        /// dimensions, and finite values are checked before mutation.
        pub fn set_positions<T: AsRef<[$crate::geometry::Point3]>>(
            &mut self,
            positions: $crate::units::Quantity<T>,
        ) -> Result<(), $crate::structure::ConformationError> {
            Ok(self.storage().positions.set_all(positions)?)
        }

        /// Sets one position by dense atom index.
        pub fn set_position_at(
            &mut self,
            atom: $crate::topology::TopologyAtomIndex,
            position: $crate::units::Quantity<$crate::geometry::Point3>,
        ) -> Result<(), $crate::structure::ConformationError> {
            Ok(self
                .storage()
                .positions
                .set_position_at(atom.index(), position)?)
        }

        pub fn set_cell(&mut self, cell: Option<$crate::geometry::PeriodicCell>) {
            self.storage().cell = cell;
        }

        /// Replaces every realization annotation with rows of the same
        /// dimensions, for example ones copied from another realization of the
        /// same topology. Detached payloads accept no bond rows.
        pub fn set_properties(
            &mut self,
            properties: $crate::properties::RealizationProperties,
        ) -> Result<(), $crate::structure::ConformationError> {
            let storage = self.storage();
            $crate::structure::state::check_property_rows(
                storage.positions.len(),
                &storage.properties,
                &properties,
            )?;
            storage.properties = properties;
            Ok(())
        }

        /// Length-preserving mutable access to realization annotations.
        pub fn properties_mut(&mut self) -> $crate::properties::RealizationPropertiesMut<'_> {
            $crate::properties::RealizationPropertiesMut::new(&mut self.storage().properties)
        }
    };
}
pub(crate) use realization_mutators;

/// Occupancy and B-factor mutators of models and ensemble members.
macro_rules! site_mutators {
    () => {
        /// Sets or clears one occupancy. Values must be finite.
        pub fn set_occupancy(
            &mut self,
            atom: $crate::topology::TopologyAtomIndex,
            value: Option<f64>,
        ) -> Result<(), $crate::structure::ConformationError> {
            let storage = self.storage();
            let len = storage.positions.len();
            storage.sites.set_occupancy(atom, len, value)
        }

        /// Replaces every occupancy in dense atom order, or clears them all.
        pub fn set_occupancies(
            &mut self,
            values: Option<Vec<Option<f64>>>,
        ) -> Result<(), $crate::structure::ConformationError> {
            let storage = self.storage();
            let len = storage.positions.len();
            storage.sites.set_occupancies(len, values)
        }

        /// Sets or clears one isotropic B-factor, converting it to square
        /// nanometres.
        pub fn set_b_factor(
            &mut self,
            atom: $crate::topology::TopologyAtomIndex,
            value: Option<$crate::units::Quantity<f64>>,
        ) -> Result<(), $crate::structure::ConformationError> {
            let storage = self.storage();
            let len = storage.positions.len();
            storage.sites.set_b_factor(atom, len, value)
        }

        /// Replaces every B-factor in dense atom order, or clears them all.
        pub fn set_b_factors(
            &mut self,
            values: Option<$crate::units::Quantity<Vec<Option<f64>>>>,
        ) -> Result<(), $crate::structure::ConformationError> {
            let storage = self.storage();
            let len = storage.positions.len();
            storage.sites.set_b_factors(len, values)
        }
    };
}
pub(crate) use site_mutators;

/// Invalid realization payload state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ConformationError {
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
    /// A dense atom index is outside the realization.
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

impl fmt::Display for ConformationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AtomCountMismatch { expected, actual } => write!(
                formatter,
                "realization requires {expected} atom rows, but received {actual}"
            ),
            Self::BondCountMismatch { expected, actual } => write!(
                formatter,
                "realization requires {expected} bond property rows, but has {actual}"
            ),
            Self::InvalidAtomIndex { index, len } => write!(
                formatter,
                "atom index {index} is outside a realization of {len} atoms"
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
            Self::CapacityOverflow => {
                formatter.write_str("realization exceeds addressable capacity")
            }
            Self::Position(error) => write!(formatter, "invalid positions: {error}"),
            Self::Property(error) => write!(formatter, "invalid realization property: {error}"),
            Self::Cell(error) => write!(formatter, "invalid periodic cell: {error}"),
            Self::Unit(error) => write!(formatter, "invalid quantity unit: {error}"),
        }
    }
}

impl std::error::Error for ConformationError {
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

impl From<PositionError> for ConformationError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<PropertyError> for ConformationError {
    fn from(error: PropertyError) -> Self {
        Self::Property(Box::new(error))
    }
}

impl From<PeriodicCellError> for ConformationError {
    fn from(error: PeriodicCellError) -> Self {
        Self::Cell(Box::new(error))
    }
}

impl From<UnitError> for ConformationError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
