//! Frame state the trajectory writers leave out or give frames themselves.

use crate::diagnostics::{GeometryField, UnrepresentablePolicy, WriteOmission, WriteReport};
use crate::structure::Geometry;
use crate::units::{Quantity, CANONICAL_TIME_UNIT};

use crate::traj::{
    TrajectoryField, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryIoOperation,
};

use super::codec_context;

/// How the XTC and TRR writers, whose formats record a time and a step for
/// every frame, treat a frame without them, such as a frame read from DCD
/// or XYZ.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[non_exhaustive]
pub enum FrameTimingPolicy {
    /// A frame without a time or a step fails to write.
    #[default]
    Require,
    /// Number the frames: frame `n` of the frames written, counted from 0,
    /// is given step `n` when it has no step and time `n × interval` when it
    /// has no time. A frame's own time and step are written as they are.
    Number {
        /// The time between consecutive written frames; a finite,
        /// nonnegative time.
        interval: Quantity<f64>,
    },
}

impl FrameTimingPolicy {
    /// Fails unless a numbering interval is a finite, nonnegative time.
    pub(crate) fn validate(
        self,
        format: TrajectoryFormat,
        source_label: &str,
    ) -> Result<(), TrajectoryIoError> {
        let Self::Number { interval } = self else {
            return Ok(());
        };
        match interval.value_in(CANONICAL_TIME_UNIT) {
            Ok(value) if value.is_finite() && value >= 0.0 => Ok(()),
            _ => Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Open,
                Some(format),
                source_label,
                format!("The {format} frame interval is not a finite, nonnegative time."),
            )),
        }
    }

    /// The step and time written for frame `frame` of the frames written,
    /// whose own are `step` and `time`.
    pub(crate) fn timing(
        self,
        frame: usize,
        step: Option<u64>,
        time: Option<Quantity<f64>>,
    ) -> (Option<u64>, Option<Quantity<f64>>) {
        match self {
            Self::Require => (step, time),
            Self::Number { interval } => (
                step.or(Some(frame as u64)),
                time.or_else(|| {
                    Some(Quantity::new(
                        frame as f64 * *interval.value(),
                        interval.unit(),
                    ))
                }),
            ),
        }
    }
}

/// Checks the frame state a format cannot carry, given as whether a frame
/// has each field and the policy for it, in the order errors name them. The
/// first present field under [`UnrepresentablePolicy::Fail`] fails through
/// `error`; otherwise the returned omissions, zero for absent fields, are for
/// [`record`] once the frame is written.
pub(crate) fn unrepresentable_frame_state<const N: usize>(
    fields: [(bool, TrajectoryField, UnrepresentablePolicy); N],
    error: impl FnOnce(TrajectoryField) -> TrajectoryIoError,
) -> Result<[WriteOmission; N], TrajectoryIoError> {
    if let Some((_, field, _)) = fields
        .iter()
        .find(|(present, _, policy)| *present && *policy == UnrepresentablePolicy::Fail)
    {
        return Err(error(*field));
    }
    Ok(fields.map(|(present, field, _)| omission(field, usize::from(present))))
}

/// Whether `geometry` records an occupancy, and whether it records a
/// B-factor: the site fields no trajectory format stores.
pub(crate) fn site_fields(geometry: &Geometry) -> (bool, bool) {
    (
        geometry
            .occupancies()
            .is_some_and(|values| values.iter().any(Option::is_some)),
        geometry
            .b_factors()
            .is_some_and(|values| values.value().iter().any(Option::is_some)),
    )
}

/// The omission of `field` from `structures` frames.
fn omission(field: TrajectoryField, structures: usize) -> WriteOmission {
    match field {
        TrajectoryField::Occupancies => WriteOmission::Occupancies { structures },
        TrajectoryField::BFactors => WriteOmission::BFactors { structures },
        field => WriteOmission::GeometryState {
            field: geometry_field(field),
            geometries: structures,
        },
    }
}

/// The frame state `field` names.
fn geometry_field(field: TrajectoryField) -> GeometryField {
    match field {
        TrajectoryField::Time => GeometryField::Time,
        TrajectoryField::Step => GeometryField::Step,
        TrajectoryField::Velocities => GeometryField::Velocities,
        TrajectoryField::Forces => GeometryField::Forces,
        TrajectoryField::AtomProperties => GeometryField::AtomProperties,
        TrajectoryField::BondProperties => GeometryField::BondProperties,
        TrajectoryField::Properties => GeometryField::Properties,
        TrajectoryField::OrderStep => GeometryField::OrderStep,
        TrajectoryField::OrderCoordinate => GeometryField::OrderCoordinate,
        other => unreachable!("{other} is not frame state a policy omits"),
    }
}

/// The trajectory field that names the geometry state `field`.
pub(crate) fn trajectory_field(field: GeometryField) -> TrajectoryField {
    match field {
        GeometryField::Weight => TrajectoryField::Weight,
        GeometryField::Energy => TrajectoryField::Energy,
        GeometryField::State => TrajectoryField::State,
        GeometryField::Temperature => TrajectoryField::Temperature,
        GeometryField::StateProperties => TrajectoryField::StateProperties,
        GeometryField::Representative => TrajectoryField::Representative,
        GeometryField::Order => TrajectoryField::Order,
        GeometryField::Membership => TrajectoryField::Membership,
        GeometryField::Time => TrajectoryField::Time,
        GeometryField::Step => TrajectoryField::Step,
        GeometryField::OrderStep => TrajectoryField::OrderStep,
        GeometryField::OrderCoordinate => TrajectoryField::OrderCoordinate,
        GeometryField::Velocities => TrajectoryField::Velocities,
        GeometryField::Forces => TrajectoryField::Forces,
        GeometryField::AtomProperties => TrajectoryField::AtomProperties,
        GeometryField::BondProperties => TrajectoryField::BondProperties,
        GeometryField::Properties => TrajectoryField::Properties,
        GeometryField::SetProperties => TrajectoryField::CollectionProperties,
        GeometryField::Occupancies => TrajectoryField::Occupancies,
        GeometryField::BFactors => TrajectoryField::BFactors,
    }
}

/// Adds `omissions` to `report`, summing counts of a kind it already has.
pub(crate) fn record(report: &mut WriteReport, omissions: impl IntoIterator<Item = WriteOmission>) {
    *report = report
        .omissions()
        .iter()
        .copied()
        .chain(omissions)
        .collect();
}
