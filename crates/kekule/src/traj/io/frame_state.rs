//! Frame state the trajectory writers leave out or give frames themselves.

use crate::diagnostics::{GeometryField, WriteOmission, WriteReport};
use crate::io::write_policy::StatePolicies;
use crate::traj::{
    FrameView, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
};
use crate::units::{Quantity, CANONICAL_TIME_UNIT};

use super::codec_context;

/// How the DCD, TRR, and XTC writers, whose formats give every frame a step,
/// treat a frame whose timing does not fit the format, such as a frame read
/// from XYZ or a stride of another trajectory.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum FrameTimingPolicy {
    /// A frame without the step and time its format needs fails to write,
    /// and so does a DCD frame off the header's step sequence.
    #[default]
    Require,
    /// Number the frames: frame `n` of the frames written, counted from 0,
    /// is given what it lacks. TRR and XTC give a frame without a step step
    /// `n`, and one without a time `n` frame intervals
    /// ([`super::trr::TrrWriteOptions::frame_interval`],
    /// [`super::xtc::XtcWriteOptions::frame_interval`]), and write a frame's
    /// own time and step as they are. DCD writes frame `n` at step
    /// `start_step + n × step_interval` of the header's sequence and, under
    /// [`super::dcd::DcdTimePolicy::HeaderDelta`], at that step × `DELTA`; a
    /// frame whose own step or time disagrees has it replaced, and the
    /// writer's report counts the frames whose steps and times were replaced
    /// as [`WriteOmission::GeometryState`] for [`GeometryField::Step`] and
    /// [`GeometryField::Time`].
    Number,
}

impl FrameTimingPolicy {
    /// The step and time a TRR or XTC writer writes for frame `frame` of the
    /// frames written, whose own are `step` and `time`, `interval` apart.
    pub(crate) fn timing(
        self,
        frame: usize,
        interval: Quantity<f64>,
        step: Option<u64>,
        time: Option<Quantity<f64>>,
    ) -> (Option<u64>, Option<Quantity<f64>>) {
        match self {
            Self::Require => (step, time),
            Self::Number => (
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

/// `interval` if it is a finite, nonnegative time, the frame interval of a
/// `format` writer's options.
pub(crate) fn frame_interval(
    format: TrajectoryFormat,
    interval: Quantity<f64>,
) -> Result<Quantity<f64>, TrajectoryIoError> {
    match interval.value_in(CANONICAL_TIME_UNIT) {
        Ok(value) if value.is_finite() && value >= 0.0 => Ok(interval),
        _ => Err(codec_context(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::Open,
            Some(format),
            &format!("{format} writer options"),
            format!("The {format} frame interval is not a finite, nonnegative time."),
        )),
    }
}

/// Checks the state of `frame`, frame `index` of a `format` writer, that the
/// format cannot carry ([`StatePolicies::check`]).
pub(crate) fn check_frame_state(
    format: TrajectoryFormat,
    source_label: &str,
    index: usize,
    frame: FrameView<'_>,
    policies: StatePolicies,
) -> Result<Vec<WriteOmission>, TrajectoryIoError> {
    policies
        .check(frame.as_structure_view())
        .map_err(|field| writer_field_error(format, source_label, index, field))
}

/// The error of a `format` writer for the `field` of frame `frame`, which
/// the format cannot carry.
pub(crate) fn writer_field_error(
    format: TrajectoryFormat,
    source_label: &str,
    frame: usize,
    field: GeometryField,
) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        TrajectoryIoErrorKind::UnsupportedField,
        TrajectoryIoOperation::WriteFrame,
        Some(format),
    )
    .with_source_label(source_label)
    .with_frame(frame)
    .with_field(field)
    .with_detail(format!("{format} cannot store the {field}."))
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
