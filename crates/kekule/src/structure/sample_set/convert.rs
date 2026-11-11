//! Lossless conversions between a [`SampleSet`] and the collections it
//! replaces: [`Ensemble`], [`Trajectory`], and [`Series`].

use std::sync::Arc;

use crate::structure::series::ItemPlace;
use crate::structure::{
    Ensemble, EnsembleSample, Series, SeriesAxis, SeriesItem, Trajectory, TrajectoryFrame,
};
use crate::units::Quantity;

use super::{
    ConversionTarget, Distribution, InMemory, Members, Order, Requirement, Sample, SampleSet,
    SampleSetError, SampleStructure, TimeAxis, Unconvertible,
};

impl Sample {
    /// A trajectory frame's positions, cell, properties, velocities, and
    /// forces as a sample, moving them without copying; the time and step
    /// are dropped.
    pub(crate) fn from_frame(frame: TrajectoryFrame) -> Self {
        let mut structure =
            SampleStructure::from_frame_state(frame.positions, frame.cell, frame.properties);
        structure.velocities = frame.velocities;
        structure.forces = frame.forces;
        Self::from_structure(Arc::new(structure))
    }
}

impl SampleSet {
    /// An unordered set of an ensemble's samples with a distribution over
    /// all of them: the ensemble's weights, partition, states, and
    /// temperature. Each sample's structure moves without copying, so a
    /// clone of the ensemble kept elsewhere still shares it.
    ///
    /// [`Self::to_ensemble`] converts back to an equal ensemble.
    pub fn from_ensemble(ensemble: Ensemble) -> Self {
        let Ensemble {
            topology,
            properties,
            samples,
            partition,
            states,
            temperature,
        } = ensemble;
        let (samples, weights) = samples
            .into_iter()
            .map(|sample| (Sample::from_structure(sample.structure), sample.weight))
            .unzip();
        let distribution = Distribution::ensemble(weights, partition, states, temperature);
        Self {
            topology,
            properties,
            storage: InMemory::new(samples),
            order: None,
            distribution: Some(distribution),
        }
    }

    /// A time-ordered set of a trajectory's frames, without a distribution.
    /// Frame times and steps become the order's [`TimeAxis`]; each frame's
    /// positions, cell, properties, velocities, and forces move to its
    /// sample without copying.
    ///
    /// [`Self::to_trajectory`] converts back to an equal trajectory.
    ///
    /// # Errors
    ///
    /// The axis has a time for every frame or for none, and so a step:
    /// [`SampleSetError::MissingTime`] or [`SampleSetError::UnexpectedTime`]
    /// for the first frame whose time differs from the first frame in being
    /// recorded, and [`SampleSetError::MissingStep`] or
    /// [`SampleSetError::UnexpectedStep`] for steps.
    pub fn from_trajectory(trajectory: Trajectory) -> Result<Self, SampleSetError> {
        let (topology, properties, frames) = trajectory.into_parts();
        let axis = TimeAxis::of_frames(&frames)?;
        let samples = frames.into_iter().map(Sample::from_frame).collect();

        Ok(Self {
            topology,
            properties,
            storage: InMemory::new(samples),
            order: Some(Order::Time(axis)),
            distribution: None,
        })
    }

    /// A set of a series' items ordered by their steps ([`Order::Step`]) or
    /// coordinates ([`Order::Coordinate`] in the axis unit), without a
    /// distribution. Each item's positions, cell, and properties move to its
    /// sample without copying.
    ///
    /// [`Self::to_series`] converts back to an equal series.
    pub fn from_series(series: Series) -> Self {
        let (topology, properties, axis, items) = series.into_parts();
        let order = Order::of_items(axis, &items);
        let samples = items
            .into_iter()
            .map(|item| {
                Sample::from_structure(Arc::new(SampleStructure::from_frame_state(
                    item.positions,
                    item.cell,
                    item.properties,
                )))
            })
            .collect();
        Self {
            topology,
            properties,
            storage: InMemory::new(samples),
            order: Some(order),
            distribution: None,
        }
    }

    /// The ensemble of this set's distribution, sharing every sample's
    /// structure. An ensemble always stores weights, so members without
    /// weights weigh `1.0`, and [`Self::from_ensemble`] of the result has
    /// weights where this set had none.
    ///
    /// # Errors
    ///
    /// `Missing(Distribution)` without a distribution, and
    /// [`SampleSetError::NotConvertible`] for an order, samples outside the
    /// distribution, or a sample with velocities or forces.
    pub fn to_ensemble(&self) -> Result<Ensemble, SampleSetError> {
        self.check_ensemble()?;
        Ok(self.clone().into_ensemble_unchecked())
    }

    /// [`Self::to_ensemble`], consuming the set.
    pub fn into_ensemble(self) -> Result<Ensemble, SampleSetError> {
        self.check_ensemble()?;
        Ok(self.into_ensemble_unchecked())
    }

    /// The trajectory of this set's samples, its frames carrying the time
    /// axis' times and steps.
    ///
    /// # Errors
    ///
    /// `Missing(TimeOrder)` without a time order, and
    /// [`SampleSetError::NotConvertible`] for a distribution or a sample
    /// with occupancies or B-factors.
    pub fn to_trajectory(&self) -> Result<Trajectory, SampleSetError> {
        self.check_trajectory()?;
        Ok(self.clone().into_trajectory_unchecked())
    }

    /// [`Self::to_trajectory`], consuming the set and moving every sample's
    /// structure that no other set shares.
    pub fn into_trajectory(self) -> Result<Trajectory, SampleSetError> {
        self.check_trajectory()?;
        Ok(self.into_trajectory_unchecked())
    }

    /// The series of this set's samples, on a step axis for a step order and
    /// a coordinate axis in the order's unit for a coordinate order.
    ///
    /// # Errors
    ///
    /// `Missing(AnyOrder)` without an order, and
    /// [`SampleSetError::NotConvertible`] for a time order, a distribution,
    /// or a sample with occupancies, B-factors, velocities, or forces.
    pub fn to_series(&self) -> Result<Series, SampleSetError> {
        self.check_series()?;
        Ok(self.clone().into_series_unchecked())
    }

    /// [`Self::to_series`], consuming the set and moving every sample's
    /// structure that no other set shares.
    pub fn into_series(self) -> Result<Series, SampleSetError> {
        self.check_series()?;
        Ok(self.into_series_unchecked())
    }

    fn check_ensemble(&self) -> Result<(), SampleSetError> {
        let target = ConversionTarget::Ensemble;
        let distribution = self
            .distribution
            .as_ref()
            .ok_or(Requirement::Distribution)?;
        if self.order.is_some() {
            return Err(not_convertible(target, Unconvertible::Order));
        }
        if distribution.members() != &Members::All {
            return Err(not_convertible(target, Unconvertible::NonMembers));
        }
        self.check_samples(target, false, true)
    }

    fn check_trajectory(&self) -> Result<(), SampleSetError> {
        let target = ConversionTarget::Trajectory;
        if !matches!(self.order, Some(Order::Time(_))) {
            return Err(Requirement::TimeOrder.into());
        }
        if self.distribution.is_some() {
            return Err(not_convertible(target, Unconvertible::Distribution));
        }
        self.check_samples(target, true, false)
    }

    fn check_series(&self) -> Result<(), SampleSetError> {
        let target = ConversionTarget::Series;
        match self.order {
            None => return Err(Requirement::AnyOrder.into()),
            Some(Order::Time(_)) => return Err(not_convertible(target, Unconvertible::Order)),
            Some(_) => {}
        }
        if self.distribution.is_some() {
            return Err(not_convertible(target, Unconvertible::Distribution));
        }
        self.check_samples(target, true, true)
    }

    /// Fails for the first sample state `target` cannot hold: occupancies
    /// and B-factors unless `sites`, velocities and forces unless `dynamics`.
    fn check_samples(
        &self,
        target: ConversionTarget,
        sites: bool,
        dynamics: bool,
    ) -> Result<(), SampleSetError> {
        for (sample, payload) in self.storage.samples.iter().enumerate() {
            let structure = payload.structure();
            let content = if sites && structure.has_occupancies() {
                Unconvertible::Occupancies { sample }
            } else if sites && structure.has_b_factors() {
                Unconvertible::BFactors { sample }
            } else if dynamics && structure.velocities.is_some() {
                Unconvertible::Velocities { sample }
            } else if dynamics && structure.forces.is_some() {
                Unconvertible::Forces { sample }
            } else {
                continue;
            };
            return Err(not_convertible(target, content));
        }
        Ok(())
    }

    fn into_ensemble_unchecked(self) -> Ensemble {
        let Self {
            topology,
            properties,
            storage: InMemory { samples },
            distribution,
            ..
        } = self;
        let Distribution {
            weights,
            partition,
            states,
            temperature,
            ..
        } = distribution.expect("ensembles are checked for a distribution");
        let samples = samples
            .into_iter()
            .enumerate()
            .map(|(index, sample)| EnsembleSample {
                structure: sample.into_structure(),
                weight: weights.as_ref().map_or(1.0, |weights| weights[index]),
            })
            .collect();
        Ensemble {
            topology,
            properties,
            samples,
            partition,
            states,
            temperature,
        }
    }

    fn into_trajectory_unchecked(self) -> Trajectory {
        let Self {
            topology,
            properties,
            storage: InMemory { samples },
            order,
            ..
        } = self;
        let Some(Order::Time(axis)) = order else {
            unreachable!("trajectories are checked for a time order")
        };
        let frames = samples
            .into_iter()
            .enumerate()
            .map(|(index, sample)| {
                let structure = owned(sample);
                TrajectoryFrame {
                    positions: structure.positions,
                    cell: structure.cell,
                    properties: structure.properties,
                    velocities: structure.velocities,
                    forces: structure.forces,
                    time: axis.time(index),
                    step: axis.step(index),
                }
            })
            .collect();
        Trajectory {
            topology,
            properties,
            frames,
        }
    }

    fn into_series_unchecked(self) -> Series {
        let Self {
            topology,
            properties,
            storage: InMemory { samples },
            order,
            ..
        } = self;
        let order = order.expect("series are checked for an order");
        let axis = match &order {
            Order::Step(_) => SeriesAxis::Step,
            Order::Coordinate { unit, .. } => SeriesAxis::Coordinate(*unit),
            Order::Time(_) => unreachable!("series are checked for a step or coordinate order"),
        };
        let items = samples
            .into_iter()
            .enumerate()
            .map(|(index, sample)| {
                let structure = owned(sample);
                let place = match &order {
                    Order::Step(steps) => ItemPlace::Step(steps[index]),
                    Order::Coordinate { unit, values } => {
                        ItemPlace::Coordinate(Quantity::new(values[index], *unit))
                    }
                    Order::Time(_) => {
                        unreachable!("series are checked for a step or coordinate order")
                    }
                };
                SeriesItem {
                    positions: structure.positions,
                    cell: structure.cell,
                    properties: structure.properties,
                    place,
                }
            })
            .collect();
        Series {
            topology,
            properties,
            axis,
            items,
        }
    }
}

impl SampleSet {
    /// Sample `index`, which must exist, copied as a trajectory frame with
    /// its place on a time order; occupancies and B-factors are dropped.
    pub(crate) fn frame(&self, index: usize) -> TrajectoryFrame {
        let structure = self.storage.samples[index].structure();
        let (time, step) = match &self.order {
            Some(Order::Time(axis)) => (axis.time(index), axis.step(index)),
            _ => (None, None),
        };
        TrajectoryFrame {
            positions: structure.positions.clone(),
            cell: structure.cell,
            properties: structure.properties.clone(),
            velocities: structure.velocities.clone(),
            forces: structure.forces.clone(),
            time,
            step,
        }
    }
}

impl Distribution {
    /// The distribution of an ensemble's samples: every one a member, with
    /// the ensemble's weights, partition, states, and temperature.
    fn ensemble(
        weights: Vec<f64>,
        partition: crate::structure::EnsemblePartition,
        states: Vec<crate::structure::EnsembleState>,
        temperature: Option<Quantity<f64>>,
    ) -> Self {
        Self {
            members: Members::All,
            weights: Some(weights),
            partition,
            states,
            temperature,
        }
    }

    /// The distribution [`SampleSet::from_ensemble`] gives `ensemble`,
    /// copying no coordinates.
    pub(crate) fn of_ensemble(ensemble: &Ensemble) -> Self {
        Self::ensemble(
            ensemble
                .samples
                .iter()
                .map(|sample| sample.weight)
                .collect(),
            ensemble.partition.clone(),
            ensemble.states.clone(),
            ensemble.temperature,
        )
    }
}

impl TimeAxis {
    /// The axis of these frames' times and steps.
    ///
    /// # Errors
    ///
    /// The axis has a time for every frame or for none, and so a step:
    /// [`SampleSetError::MissingTime`] or [`SampleSetError::UnexpectedTime`]
    /// for the first frame whose time differs from the first frame in being
    /// recorded, and [`SampleSetError::MissingStep`] or
    /// [`SampleSetError::UnexpectedStep`] for steps.
    fn of_frames(frames: &[TrajectoryFrame]) -> Result<Self, SampleSetError> {
        let timed = frames.first().is_some_and(|frame| frame.time.is_some());
        let stepped = frames.first().is_some_and(|frame| frame.step.is_some());
        for (sample, frame) in frames.iter().enumerate() {
            match (timed, frame.time.is_some()) {
                (true, false) => return Err(SampleSetError::MissingTime { sample }),
                (false, true) => return Err(SampleSetError::UnexpectedTime { sample }),
                _ => {}
            }
            match (stepped, frame.step.is_some()) {
                (true, false) => return Err(SampleSetError::MissingStep { sample }),
                (false, true) => return Err(SampleSetError::UnexpectedStep { sample }),
                _ => {}
            }
        }
        let times = frames
            .iter()
            .map(|frame| frame.time.map(Quantity::into_value))
            .collect::<Option<Vec<_>>>();
        let steps = frames
            .iter()
            .map(|frame| frame.step)
            .collect::<Option<Vec<_>>>();
        Ok(Self::from_columns(
            times.filter(|_| timed),
            steps.filter(|_| stepped),
        ))
    }
}

impl Order {
    /// The time order [`SampleSet::from_trajectory`] gives `trajectory`,
    /// with its errors.
    pub(crate) fn of_trajectory(trajectory: &Trajectory) -> Result<Self, SampleSetError> {
        TimeAxis::of_frames(&trajectory.frames).map(Self::Time)
    }

    /// The step or coordinate order [`SampleSet::from_series`] gives
    /// `series`.
    pub(crate) fn of_series(series: &Series) -> Self {
        Self::of_items(series.axis, &series.items)
    }

    fn of_items(axis: SeriesAxis, items: &[SeriesItem]) -> Self {
        match axis {
            SeriesAxis::Step => Self::Step(
                items
                    .iter()
                    .map(|item| item.step().expect("step items lie on a step axis"))
                    .collect(),
            ),
            SeriesAxis::Coordinate(unit) => Self::Coordinate {
                unit,
                values: items
                    .iter()
                    .map(|item| {
                        item.coordinate()
                            .expect("coordinate items lie on a coordinate axis")
                            .into_value()
                    })
                    .collect(),
            },
        }
    }
}

/// A sample's structure, moved when no other holder shares it.
fn owned(sample: Sample) -> SampleStructure {
    Arc::try_unwrap(sample.into_structure()).unwrap_or_else(|shared| (*shared).clone())
}

const fn not_convertible(target: ConversionTarget, content: Unconvertible) -> SampleSetError {
    SampleSetError::NotConvertible { target, content }
}
