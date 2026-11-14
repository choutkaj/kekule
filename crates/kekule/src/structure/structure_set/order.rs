//! What the order of a [`super::StructureSet`]'s geometries means, and each
//! geometry's place in it.

use crate::diagnostics::{Countable, LengthMismatch};
use crate::units::{Quantity, Unit, CANONICAL_TIME_UNIT};

use super::StructureSetError;

/// What the order of a structure set's geometries means, with each geometry's value.
///
/// A set without an order is unordered, as a distribution of conformations
/// is. Order values belong to the order, never to a geometry: a column has a
/// value for every geometry or for none, and a geometry taken out of its set
/// claims no place. Values need not increase, as a scan swept back and forth
/// shows, and selecting geometries keeps each one's value.
///
/// [`super::StructureSet::set_order`] checks that the order has one value per
/// geometry and that every coordinate is finite.
///
/// ```
/// use kekule::structure::{Order, TimeAxis};
/// use kekule::units::{Quantity, DEGREE, PICOSECOND};
/// // A trajectory's frame times, an optimization's steps, a dihedral scan.
/// let time = Order::Time(TimeAxis::new(
///     Some(Quantity::new(vec![0.0, 2.0, 4.0], PICOSECOND)),
///     None,
/// )?);
/// let steps = Order::Step(vec![0, 1, 2]);
/// let scan = Order::Coordinate { unit: DEGREE, values: vec![-60.0, 0.0, 60.0] };
/// assert_ne!(steps, scan);
/// # let _ = time;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Order {
    /// A path in time, as dynamics produces: order and spacing are physical.
    /// The axis may carry each geometry's time and MD step.
    Time(TimeAxis),
    /// Numbered steps: the accepted iterations of a geometry optimization,
    /// or the images of a path. One step per geometry.
    Step(Vec<u64>),
    /// The values of one scanned or reaction coordinate, stored in `unit`,
    /// one per geometry. A coordinate may be dimensionless, as a normalized
    /// path length is.
    Coordinate {
        /// The unit the values are stored in.
        unit: Unit,
        /// One finite value per geometry.
        values: Vec<f64>,
    },
}

/// The values of a time order: each geometry's time and MD step, each column
/// present for every geometry or for none.
///
/// An axis without columns orders geometries in time without stating when they
/// were taken.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TimeAxis {
    /// Finite times in [`CANONICAL_TIME_UNIT`].
    times: Option<Vec<f64>>,
    steps: Option<Vec<u64>>,
}

impl TimeAxis {
    /// An axis with the given columns. Times are converted to the canonical
    /// time unit (ps).
    ///
    /// # Errors
    ///
    /// [`StructureSetError::Unit`] when the times are not times,
    /// [`StructureSetError::NonFiniteOrderValue`] for the first time that is
    /// not finite, and [`StructureSetError::TimeAxisMismatch`] when the columns
    /// have different lengths.
    pub fn new(
        times: Option<Quantity<Vec<f64>>>,
        steps: Option<Vec<u64>>,
    ) -> Result<Self, StructureSetError> {
        let times = times
            .map(|times| {
                let factor = times.unit().conversion_factor_to(CANONICAL_TIME_UNIT)?;
                times
                    .into_value()
                    .into_iter()
                    .enumerate()
                    .map(|(geometry, time)| {
                        let time = time * factor;
                        if time.is_finite() {
                            Ok(time)
                        } else {
                            Err(StructureSetError::NonFiniteOrderValue { geometry })
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .transpose()?;
        if let (Some(times), Some(steps)) = (&times, &steps) {
            if times.len() != steps.len() {
                return Err(StructureSetError::TimeAxisMismatch {
                    times: times.len(),
                    steps: steps.len(),
                });
            }
        }
        Ok(Self { times, steps })
    }

    /// Built from checked columns.
    pub(crate) const fn from_columns(times: Option<Vec<f64>>, steps: Option<Vec<u64>>) -> Self {
        Self { times, steps }
    }

    /// Every geometry's time in [`CANONICAL_TIME_UNIT`], or `None` when the
    /// axis has no times.
    pub fn times(&self) -> Option<Quantity<&[f64]>> {
        self.times
            .as_deref()
            .map(|times| Quantity::new(times, CANONICAL_TIME_UNIT))
    }

    /// Every geometry's MD step, or `None` when the axis has no steps.
    pub fn steps(&self) -> Option<&[u64]> {
        self.steps.as_deref()
    }

    /// One geometry's time, or `None` without times or out of range.
    pub fn time(&self, geometry: usize) -> Option<Quantity<f64>> {
        let time = *self.times.as_ref()?.get(geometry)?;
        Some(Quantity::new(time, CANONICAL_TIME_UNIT))
    }

    /// One geometry's MD step, or `None` without steps or out of range.
    pub fn step(&self, geometry: usize) -> Option<u64> {
        self.steps.as_ref()?.get(geometry).copied()
    }

    fn len(&self) -> Option<usize> {
        self.times
            .as_ref()
            .map(Vec::len)
            .or_else(|| self.steps.as_ref().map(Vec::len))
    }
}

/// A geometry's place in its set's [`Order`]: what a geometry is given when it
/// is inserted into an ordered set, and what its view reports.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Place {
    /// A place on a time order: the geometry's time and MD step, each present
    /// exactly when the axis has that column.
    Time {
        /// The geometry's time.
        time: Option<Quantity<f64>>,
        /// The geometry's MD step.
        step: Option<u64>,
    },
    /// A place on a step order.
    Step(u64),
    /// A place on a coordinate order, converted to the order's unit.
    Coordinate(Quantity<f64>),
}

/// A checked place in an order's own representation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Stored {
    /// A time and step in canonical units, each present when its column is.
    Time(Option<f64>, Option<u64>),
    Step(u64),
    /// In the order's unit.
    Coordinate(f64),
}

impl Order {
    /// Checks that the order has one value per geometry of a set of `len`
    /// and that every coordinate is finite.
    pub(crate) fn validate(&self, len: usize) -> Result<(), StructureSetError> {
        let count = match self {
            Self::Time(axis) => axis.len(),
            Self::Step(steps) => Some(steps.len()),
            Self::Coordinate { values, .. } => {
                if let Some(geometry) = values.iter().position(|value| !value.is_finite()) {
                    return Err(StructureSetError::NonFiniteOrderValue { geometry });
                }
                Some(values.len())
            }
        };
        match count {
            Some(actual) if actual != len => Err(StructureSetError::OrderLengthMismatch(
                LengthMismatch::new(Countable::OrderValue, len, actual),
            )),
            _ => Ok(()),
        }
    }

    /// The place of geometry `index`, which must be in its set.
    pub(crate) fn place(&self, index: usize) -> Place {
        match self {
            Self::Time(axis) => Place::Time {
                time: axis.time(index),
                step: axis.step(index),
            },
            Self::Step(steps) => Place::Step(steps[index]),
            Self::Coordinate { unit, values } => {
                Place::Coordinate(Quantity::new(values[index], *unit))
            }
        }
    }

    /// Checks `place` for geometry `geometry` of this order and converts it.
    pub(crate) fn check(
        &self,
        place: Option<Place>,
        geometry: usize,
    ) -> Result<Stored, StructureSetError> {
        match (self, place) {
            (Self::Time(axis), None) if axis.times.is_none() && axis.steps.is_none() => {
                Ok(Stored::Time(None, None))
            }
            (Self::Time(axis), Some(Place::Time { time, step })) => {
                let time = match (axis.times.is_some(), time) {
                    (true, Some(time)) => {
                        let time = time.value_in(CANONICAL_TIME_UNIT)?;
                        if !time.is_finite() {
                            return Err(StructureSetError::NonFiniteOrderValue { geometry });
                        }
                        Some(time)
                    }
                    (false, None) => None,
                    (true, None) => return Err(StructureSetError::MissingTime { geometry }),
                    (false, Some(_)) => return Err(StructureSetError::UnexpectedTime { geometry }),
                };
                let step = match (axis.steps.is_some(), step) {
                    (true, Some(step)) => Some(step),
                    (false, None) => None,
                    (true, None) => return Err(StructureSetError::MissingStep { geometry }),
                    (false, Some(_)) => return Err(StructureSetError::UnexpectedStep { geometry }),
                };
                Ok(Stored::Time(time, step))
            }
            (Self::Step(_), Some(Place::Step(step))) => Ok(Stored::Step(step)),
            (Self::Coordinate { unit, .. }, Some(Place::Coordinate(coordinate))) => {
                let value = coordinate.value_in(*unit)?;
                if !value.is_finite() {
                    return Err(StructureSetError::NonFiniteOrderValue { geometry });
                }
                Ok(Stored::Coordinate(value))
            }
            _ => Err(StructureSetError::PlaceMismatch { geometry }),
        }
    }

    /// Appends a place [`Self::check`] returned.
    pub(crate) fn push(&mut self, place: Stored) {
        match (self, place) {
            (Self::Time(axis), Stored::Time(time, step)) => {
                if let (Some(times), Some(time)) = (&mut axis.times, time) {
                    times.push(time);
                }
                if let (Some(steps), Some(step)) = (&mut axis.steps, step) {
                    steps.push(step);
                }
            }
            (Self::Step(steps), Stored::Step(step)) => steps.push(step),
            (Self::Coordinate { values, .. }, Stored::Coordinate(value)) => values.push(value),
            _ => unreachable!("checked places match their order"),
        }
    }

    /// Replaces the place of geometry `index` with one [`Self::check`]
    /// returned.
    pub(crate) fn set(&mut self, index: usize, place: Stored) {
        match (self, place) {
            (Self::Time(axis), Stored::Time(time, step)) => {
                if let (Some(times), Some(time)) = (&mut axis.times, time) {
                    times[index] = time;
                }
                if let (Some(steps), Some(step)) = (&mut axis.steps, step) {
                    steps[index] = step;
                }
            }
            (Self::Step(steps), Stored::Step(step)) => steps[index] = step,
            (Self::Coordinate { values, .. }, Stored::Coordinate(value)) => values[index] = value,
            _ => unreachable!("checked places match their order"),
        }
    }

    /// Removes the place of geometry `index`.
    pub(crate) fn remove(&mut self, index: usize) {
        match self {
            Self::Time(axis) => {
                if let Some(times) = &mut axis.times {
                    times.remove(index);
                }
                if let Some(steps) = &mut axis.steps {
                    steps.remove(index);
                }
            }
            Self::Step(steps) => {
                steps.remove(index);
            }
            Self::Coordinate { values, .. } => {
                values.remove(index);
            }
        }
    }

    /// The places of the given geometries, in the given order.
    pub(crate) fn select(&self, indices: &[usize]) -> Self {
        fn pick<T: Copy>(values: &[T], indices: &[usize]) -> Vec<T> {
            indices.iter().map(|index| values[*index]).collect()
        }
        match self {
            Self::Time(axis) => Self::Time(TimeAxis {
                times: axis.times.as_deref().map(|times| pick(times, indices)),
                steps: axis.steps.as_deref().map(|steps| pick(steps, indices)),
            }),
            Self::Step(steps) => Self::Step(pick(steps, indices)),
            Self::Coordinate { unit, values } => Self::Coordinate {
                unit: *unit,
                values: pick(values, indices),
            },
        }
    }
}
