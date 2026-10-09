use std::fmt;

use crate::geometry::{PeriodicCell, PeriodicCellError, Point3, RigidTransform};
use crate::properties::{PropertyError, RealizationProperties, RealizationPropertiesMut};
use crate::topology::{Topology, TopologyAtomIndex};
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT, SQUARE_NANOMETER};

use super::{PositionError, Positions};

/// One topology-free set of coordinate-dependent atom state.
///
/// A conformation stores positions, an optional periodic cell, optional
/// per-atom occupancies and B-factors, and [`RealizationProperties`]. It is the
/// payload shared by every realization: a [`super::Model`], each ensemble
/// member, and each trajectory frame wrap one conformation and bind it to a
/// topology.
///
/// Atom state is addressed by [`TopologyAtomIndex`] in the dense atom order of
/// the topology the conformation is (or will be) bound to. Every mutation
/// preserves the atom count. Bond property rows are allocated when a detached
/// conformation is first bound to a topology.
///
/// Owners hand out [`ConformationMut`] rather than `&mut Conformation`, so a
/// bound conformation can never be replaced by one with other dimensions.
#[derive(Debug, Clone, PartialEq)]
pub struct Conformation {
    positions: Positions,
    cell: Option<PeriodicCell>,
    occupancies: AtomScalars,
    b_factors: AtomScalars,
    properties: RealizationProperties,
}

impl From<Positions> for Conformation {
    fn from(positions: Positions) -> Self {
        Self::new(positions)
    }
}

impl Conformation {
    /// Creates a non-periodic conformation without occupancies, B-factors, or
    /// properties.
    pub fn new(positions: Positions) -> Self {
        let atoms = positions.len();
        Self {
            positions,
            cell: None,
            occupancies: AtomScalars::default(),
            b_factors: AtomScalars::default(),
            properties: RealizationProperties::new(atoms),
        }
    }

    /// Returns this conformation with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    pub fn positions(&self) -> &Positions {
        &self.positions
    }

    pub fn cell(&self) -> Option<&PeriodicCell> {
        self.cell.as_ref()
    }

    /// Occupancy of one atom; `None` when it was never recorded.
    pub fn occupancy(&self, atom: TopologyAtomIndex) -> Result<Option<f64>, ConformationError> {
        self.occupancies.get(self.atom_row(atom)?)
    }

    /// Every occupancy in dense atom order, or `None` when none was recorded.
    pub fn occupancies(&self) -> Option<&[Option<f64>]> {
        self.occupancies.values()
    }

    /// Isotropic B-factor of one atom; `None` when it was never recorded.
    pub fn b_factor(
        &self,
        atom: TopologyAtomIndex,
    ) -> Result<Option<Quantity<f64>>, ConformationError> {
        Ok(self
            .b_factors
            .get(self.atom_row(atom)?)?
            .map(|value| Quantity::new(value, SQUARE_NANOMETER)))
    }

    /// Every B-factor in dense atom order and square nanometres, or `None` when
    /// none was recorded.
    pub fn b_factors(&self) -> Option<Quantity<&[Option<f64>]>> {
        self.b_factors
            .values()
            .map(|values| Quantity::new(values, SQUARE_NANOMETER))
    }

    pub fn properties(&self) -> &RealizationProperties {
        &self.properties
    }

    fn atom_row(&self, atom: TopologyAtomIndex) -> Result<usize, ConformationError> {
        let index = atom.index();
        if index >= self.positions.len() {
            return Err(ConformationError::InvalidAtomIndex {
                index,
                len: self.positions.len(),
            });
        }
        Ok(index)
    }

    /// Checks that this conformation can be bound to `topology` without
    /// changing anything.
    pub(crate) fn validate_for(&self, topology: &Topology) -> Result<(), ConformationError> {
        let atoms = topology.atom_count();
        if self.positions.len() != atoms {
            return Err(ConformationError::AtomCountMismatch {
                expected: atoms,
                actual: self.positions.len(),
            });
        }
        self.properties
            .validate_bond_rows(topology.bond_count())
            .map_err(|_| ConformationError::BondCountMismatch {
                expected: topology.bond_count(),
                actual: self.properties.bonds().len(),
            })
    }

    /// Whether this conformation already has every row `topology` requires.
    pub(crate) fn is_bound_to(&self, topology: &Topology) -> bool {
        self.positions.len() == topology.atom_count()
            && self.properties.bonds().len() == topology.bond_count()
    }

    /// Allocates bond rows after [`Self::validate_for`] succeeded.
    pub(crate) fn bind(&mut self, topology: &Topology) {
        self.properties.bind_bond_rows(topology.bond_count());
    }

    /// Projects atom and bond rows into a new conformation; owner properties
    /// are dropped because the projection is a new owner.
    pub(crate) fn project(
        &self,
        atoms: &[usize],
        bonds: &[usize],
    ) -> Result<Self, ConformationError> {
        Ok(Self {
            positions: self.positions.select_indices(atoms)?,
            cell: self.cell,
            occupancies: self.occupancies.select(atoms),
            b_factors: self.b_factors.select(atoms),
            properties: self.properties.project(atoms, bonds)?,
        })
    }

    /// Applies a rigid transform to positions and cell vectors.
    pub(crate) fn transformed(&self, transform: RigidTransform) -> Result<Self, ConformationError> {
        let mut result = self.clone();
        result.transform_in_place(transform)?;
        Ok(result)
    }

    pub(crate) fn transform_in_place(
        &mut self,
        transform: RigidTransform,
    ) -> Result<(), ConformationError> {
        let cell = self
            .cell
            .map(|cell| transform_cell(cell, transform))
            .transpose()?;
        let points = self
            .positions
            .values()
            .value()
            .iter()
            .map(|point| transform.transform_point(*point))
            .collect::<Vec<_>>();
        let positions = Positions::new(Quantity::new(points, CANONICAL_LENGTH_UNIT))?;
        self.positions = positions;
        self.cell = cell;
        Ok(())
    }

    /// Appends atoms and bond rows for staged construction. Owner properties
    /// are kept; callers apply their own structural-edit policy.
    pub(crate) fn extend(&mut self, positions: &Positions, added_bonds: usize) {
        self.positions.extend_canonical(positions);
        let atoms = self.positions.len();
        self.occupancies.resize(atoms);
        self.b_factors.resize(atoms);
        let bonds = self.properties.bonds().len() + added_bonds;
        self.properties.resize(atoms, bonds);
    }

    /// Grows private editor slot rows; new atom slots start at the origin and
    /// are positioned by the caller.
    pub(crate) fn resize_slots(&mut self, atoms: usize, bonds: usize) {
        debug_assert!(atoms >= self.positions.len());
        self.positions.resize_canonical(atoms);
        self.occupancies.resize(atoms);
        self.b_factors.resize(atoms);
        self.properties.resize(atoms, bonds);
    }

    /// Clears annotations of editor atom slots that are no longer live.
    pub(crate) fn clear_dead_atom_rows(&mut self, live: &[usize]) {
        let live = live
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>();
        for slot in (0..self.positions.len()).filter(|slot| !live.contains(slot)) {
            self.occupancies.clear(slot);
            self.b_factors.clear(slot);
            self.properties.atoms_raw_mut().clear_index(slot);
        }
    }

    /// Clears annotations of editor bond slots that are no longer live.
    pub(crate) fn clear_dead_bond_rows(&mut self, live: &[usize]) {
        let live = live
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>();
        for slot in (0..self.properties.bonds().len()).filter(|slot| !live.contains(slot)) {
            self.properties.bonds_raw_mut().clear_index(slot);
        }
    }

    /// Copies positions, occupancies, and B-factors of every source atom into
    /// the given destination rows. Property rows are copied by the caller.
    pub(crate) fn copy_atom_state_from(&mut self, source: &Self, rows: &[usize]) {
        debug_assert_eq!(rows.len(), source.atom_count());
        let points = source.positions.values();
        let staged = rows
            .iter()
            .copied()
            .zip(points.value().iter().copied())
            .collect::<Vec<_>>();
        self.positions
            .set_canonical_batch(&staged)
            .expect("source positions are finite and rows are allocated");
        let len = self.positions.len();
        for (scalars, source) in [
            (&mut self.occupancies, &source.occupancies),
            (&mut self.b_factors, &source.b_factors),
        ] {
            if let Some(values) = source.values() {
                for (&row, value) in rows.iter().zip(values) {
                    scalars
                        .set(row, len, *value)
                        .expect("source values are finite");
                }
            }
        }
    }

    pub(crate) fn try_reserve_atoms(&mut self, additional: usize) -> Result<(), ConformationError> {
        self.positions
            .try_reserve(additional)
            .map_err(|_| ConformationError::CapacityOverflow)
    }

    /// Reorders atom rows; `previous[new]` is the old dense index.
    pub(crate) fn reorder_atoms(&mut self, previous: &[usize]) {
        self.positions = self
            .positions
            .select_indices(previous)
            .expect("staged positions follow staged atoms");
        self.occupancies = self.occupancies.select(previous);
        self.b_factors = self.b_factors.select(previous);
        self.properties
            .reorder_atoms(previous)
            .expect("staged realization atom rows follow staged atoms");
    }

    pub(crate) fn positions_mut(&mut self) -> &mut Positions {
        &mut self.positions
    }

    pub(crate) fn properties_storage_mut(&mut self) -> &mut RealizationProperties {
        &mut self.properties
    }

    fn storage(&mut self) -> &mut Self {
        self
    }
}

macro_rules! conformation_mutators {
    () => {
        /// Replaces every position without changing the atom count. Units,
        /// dimensions, and finite values are checked before mutation.
        pub fn set_positions<T: AsRef<[Point3]>>(
            &mut self,
            positions: Quantity<T>,
        ) -> Result<(), ConformationError> {
            Ok(self.storage().positions.set_all(positions)?)
        }

        pub fn set_position(
            &mut self,
            atom: TopologyAtomIndex,
            position: Quantity<Point3>,
        ) -> Result<(), ConformationError> {
            Ok(self
                .storage()
                .positions
                .set_position_at(atom.index(), position)?)
        }

        pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
            self.storage().cell = cell;
        }

        /// Sets or clears one occupancy. Values must be finite.
        pub fn set_occupancy(
            &mut self,
            atom: TopologyAtomIndex,
            value: Option<f64>,
        ) -> Result<(), ConformationError> {
            let storage = self.storage();
            let row = storage.atom_row(atom)?;
            let len = storage.positions.len();
            storage.occupancies.set(row, len, value)
        }

        /// Replaces every occupancy in dense atom order, or clears them all.
        pub fn set_occupancies(
            &mut self,
            values: Option<Vec<Option<f64>>>,
        ) -> Result<(), ConformationError> {
            let storage = self.storage();
            let len = storage.positions.len();
            storage.occupancies.replace(len, values)
        }

        /// Sets or clears one isotropic B-factor, converting it to square nanometres.
        pub fn set_b_factor(
            &mut self,
            atom: TopologyAtomIndex,
            value: Option<Quantity<f64>>,
        ) -> Result<(), ConformationError> {
            let value = value
                .map(|value| value.into_unit(SQUARE_NANOMETER).map(Quantity::into_value))
                .transpose()?;
            let storage = self.storage();
            let row = storage.atom_row(atom)?;
            let len = storage.positions.len();
            storage.b_factors.set(row, len, value)
        }

        /// Replaces every B-factor in dense atom order, or clears them all.
        pub fn set_b_factors(
            &mut self,
            values: Option<Quantity<Vec<Option<f64>>>>,
        ) -> Result<(), ConformationError> {
            let values = values
                .map(|values| {
                    let factor = values.unit().conversion_factor_to(SQUARE_NANOMETER)?;
                    Ok::<_, ConformationError>(
                        values
                            .into_value()
                            .into_iter()
                            .map(|value| value.map(|value| value * factor))
                            .collect::<Vec<_>>(),
                    )
                })
                .transpose()?;
            let storage = self.storage();
            let len = storage.positions.len();
            storage.b_factors.replace(len, values)
        }

        /// Replaces every realization annotation with rows of the same
        /// dimensions, for example ones copied from another conformation of
        /// the same topology. Detached conformations accept no bond rows.
        pub fn set_properties(
            &mut self,
            properties: RealizationProperties,
        ) -> Result<(), ConformationError> {
            let storage = self.storage();
            if properties.atoms().len() != storage.positions.len() {
                return Err(ConformationError::AtomCountMismatch {
                    expected: storage.positions.len(),
                    actual: properties.atoms().len(),
                });
            }
            if properties.bonds().len() != storage.properties.bonds().len() {
                return Err(ConformationError::BondCountMismatch {
                    expected: storage.properties.bonds().len(),
                    actual: properties.bonds().len(),
                });
            }
            storage.properties = properties;
            Ok(())
        }

        /// Length-preserving mutable access to realization annotations.
        pub fn properties_mut(&mut self) -> RealizationPropertiesMut<'_> {
            RealizationPropertiesMut::new(&mut self.storage().properties)
        }
    };
}

impl Conformation {
    conformation_mutators!();
}

/// Length-preserving mutable access to a conformation owned by a realization.
///
/// Reads are available through `Deref`. There is deliberately no `DerefMut`:
/// the owner's dimensions can only change through its own replacement methods.
///
/// ```compile_fail
/// use kekule::{smiles, structure::{Conformation, Model, Positions}};
/// let mut model = Model::new(smiles::to_topology("C").unwrap(), Positions::zeros(1)).unwrap();
/// *model.conformation_mut() = Conformation::new(Positions::zeros(2));
/// ```
#[derive(Debug)]
pub struct ConformationMut<'a> {
    conformation: &'a mut Conformation,
}

impl std::ops::Deref for ConformationMut<'_> {
    type Target = Conformation;

    fn deref(&self) -> &Self::Target {
        self.conformation
    }
}

impl<'a> ConformationMut<'a> {
    pub(crate) fn new(conformation: &'a mut Conformation) -> Self {
        Self { conformation }
    }

    fn storage(&mut self) -> &mut Conformation {
        self.conformation
    }

    conformation_mutators!();
}

fn transform_cell(
    cell: PeriodicCell,
    transform: RigidTransform,
) -> Result<PeriodicCell, ConformationError> {
    let vectors = cell
        .vectors()
        .map(|vectors| vectors.map(|vector| transform.transform_vector(vector)));
    Ok(PeriodicCell::new(vectors, cell.periodic_axes())?)
}

/// Optional dense per-atom scalars, allocated on the first recorded value.
#[derive(Debug, Clone, Default)]
struct AtomScalars(Option<Vec<Option<f64>>>);

impl PartialEq for AtomScalars {
    fn eq(&self, other: &Self) -> bool {
        match (&self.0, &other.0) {
            (Some(left), Some(right)) => left == right,
            (Some(values), None) | (None, Some(values)) => values.iter().all(Option::is_none),
            (None, None) => true,
        }
    }
}

impl AtomScalars {
    fn values(&self) -> Option<&[Option<f64>]> {
        self.0.as_deref()
    }

    fn get(&self, row: usize) -> Result<Option<f64>, ConformationError> {
        Ok(self.0.as_ref().and_then(|values| values[row]))
    }

    fn set(&mut self, row: usize, len: usize, value: Option<f64>) -> Result<(), ConformationError> {
        if value.is_some_and(|value| !value.is_finite()) {
            return Err(ConformationError::NonFiniteValue { index: row });
        }
        match (&mut self.0, value) {
            (Some(values), value) => values[row] = value,
            (None, Some(value)) => {
                let mut values = vec![None; len];
                values[row] = Some(value);
                self.0 = Some(values);
            }
            (None, None) => {}
        }
        Ok(())
    }

    fn replace(
        &mut self,
        len: usize,
        values: Option<Vec<Option<f64>>>,
    ) -> Result<(), ConformationError> {
        if let Some(values) = &values {
            if values.len() != len {
                return Err(ConformationError::AtomCountMismatch {
                    expected: len,
                    actual: values.len(),
                });
            }
            if let Some(index) = values
                .iter()
                .position(|value| value.is_some_and(|value| !value.is_finite()))
            {
                return Err(ConformationError::NonFiniteValue { index });
            }
        }
        self.0 = values.filter(|values| values.iter().any(Option::is_some));
        Ok(())
    }

    fn select(&self, rows: &[usize]) -> Self {
        Self(
            self.0
                .as_ref()
                .map(|values| rows.iter().map(|row| values[*row]).collect()),
        )
    }

    fn resize(&mut self, len: usize) {
        if let Some(values) = &mut self.0 {
            values.resize(len, None);
        }
    }

    fn clear(&mut self, row: usize) {
        if let Some(values) = &mut self.0 {
            values[row] = None;
        }
    }
}

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
    /// A dense atom index is outside the conformation.
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
    /// An ensemble weight is not finite and non-negative.
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
                formatter.write_str("ensemble weight must be finite and non-negative")
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
