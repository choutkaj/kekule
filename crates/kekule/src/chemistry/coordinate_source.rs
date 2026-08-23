use crate::core::AtomId;
use crate::geometry::Point3;
/// Crate-private coordinate lookup used while normalizing represented stereo.
///
/// Public geometry remains dense `Positions`. Format interpreters may
/// implement this trait for private sparse source-coordinate staging.
pub(crate) trait AtomPositionSource {
    fn position_value(&self, atom: AtomId) -> Option<Point3>;
}
