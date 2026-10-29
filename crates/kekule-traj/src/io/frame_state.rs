//! Frame state the trajectory writers leave out or give frames themselves.

use kekule::diagnostics::{StructureField, UnrepresentablePolicy, WriteOmission, WriteReport};
use kekule::units::{Quantity, CANONICAL_TIME_UNIT};

use crate::{
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
                format!("{format} frame numbering interval must be a finite, nonnegative time"),
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
    Ok(
        fields.map(|(present, field, _)| WriteOmission::StructureState {
            field: structure_field(field),
            structures: usize::from(present),
        }),
    )
}

/// The frame state `field` names.
fn structure_field(field: TrajectoryField) -> StructureField {
    match field {
        TrajectoryField::Time => StructureField::Time,
        TrajectoryField::Step => StructureField::Step,
        TrajectoryField::Velocities => StructureField::Velocities,
        TrajectoryField::Forces => StructureField::Forces,
        TrajectoryField::AtomProperties => StructureField::AtomProperties,
        TrajectoryField::BondProperties => StructureField::BondProperties,
        TrajectoryField::Properties => StructureField::Properties,
        TrajectoryField::ItemStep => StructureField::ItemStep,
        TrajectoryField::ItemCoordinate => StructureField::ItemCoordinate,
        other => unreachable!("{other} is not frame state a policy omits"),
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
