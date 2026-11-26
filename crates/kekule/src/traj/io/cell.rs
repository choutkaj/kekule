//! The periodic boxes of GROMACS frames, and the cell orientation the
//! GROMACS and DCD writers require.

use crate::traj::trajectory::reason::Reason;

use crate::math::{PeriodicCell, Vector3};
use crate::traj::{
    TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
};
use crate::units::{Quantity, NANOMETER};

use super::codec::frame_error;
use super::codec_context;

/// How an XTC or TRR box of zeros reads, and how the XTC writer writes a
/// frame without a cell.
///
/// GROMACS writes a box of zeros for a system without periodic boundaries:
/// an XTC frame always holds a box, and `mdrun` writes one into every TRR
/// frame it writes. A TRR frame may also hold no box at all, which reads as
/// a frame without a cell under either policy.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ZeroBoxPolicy {
    /// A box of zeros fails to read with
    /// [`TrajectoryIoErrorKind::UnsupportedVariant`], and the XTC writer
    /// refuses a frame without a cell.
    #[default]
    Reject,
    /// A box of zeros reads as a frame without a cell, and the XTC writer
    /// writes a frame without a cell with a box of zeros.
    NoCell,
}

/// The cell of a GROMACS box of finite `vectors`, in nanometres, of frame
/// `frame`: `None` for a box of zeros under [`ZeroBoxPolicy::NoCell`].
pub(crate) fn gromacs_cell(
    vectors: [Vector3; 3],
    policy: ZeroBoxPolicy,
    format: TrajectoryFormat,
    source_label: &str,
    frame: usize,
) -> Result<Option<PeriodicCell>, TrajectoryIoError> {
    let zero = vectors
        .iter()
        .all(|vector| vector.x == 0.0 && vector.y == 0.0 && vector.z == 0.0);
    if zero {
        return match policy {
            ZeroBoxPolicy::NoCell => Ok(None),
            ZeroBoxPolicy::Reject => Err(frame_error(
                format,
                TrajectoryIoErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame,
                Reason::BoxZero,
            )),
        };
    }
    PeriodicCell::new(Quantity::new(vectors, NANOMETER), [true; 3])
        .map(Some)
        .map_err(|error| {
            frame_error(
                format,
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::ReadFrame,
                source_label,
                frame,
                Reason::InvalidPeriodicCell,
            )
            .with_source(error)
        })
}

/// How far from zero the elements above a standard cell's diagonal may lie,
/// relative to its longest vector: the rounding a fit or a unit
/// conversion leaves.
const ORIENTATION_TOLERANCE: f64 = 1e-9;

/// The cell `vectors` in the orientation GROMACS and DCD readers assume:
/// vector a along +x, b in the xy plane with positive y, c with positive z.
/// GROMACS ignores the box elements above its diagonal, and DCD keeps only
/// lengths and angles. Elements above the diagonal within rounding of zero
/// ([`ORIENTATION_TOLERANCE`]) are returned as zeros; others fail.
pub(crate) fn standard_orientation(
    [a, b, c]: [Vector3; 3],
    format: TrajectoryFormat,
    source_label: &str,
) -> Result<[Vector3; 3], TrajectoryIoError> {
    let longest = [a, b, c]
        .iter()
        .map(|vector| vector.norm())
        .fold(0.0, f64::max);
    let tolerance = ORIENTATION_TOLERANCE * longest;
    let tilted = [a.y, a.z, b.z].iter().any(|value| value.abs() > tolerance);
    if tilted || a.x <= 0.0 || b.y <= 0.0 || c.z <= 0.0 {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::WriteFrame,
            Some(format),
            source_label,
            Reason::CellNotStandard,
        ));
    }
    Ok([Vector3::new(a.x, 0.0, 0.0), Vector3::new(b.x, b.y, 0.0), c])
}
