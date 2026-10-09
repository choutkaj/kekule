//! Atomic Cartesian editing on [`Model`].
//!
//! Setters choose a connected moving fragment; prepared edits resolve that
//! fragment once for repeated absolute targets. Preparation stores no geometry.
//! All operations retain the exact shared topology, properties, and cell. Cells
//! are ignored: image or reconstruct periodic molecules explicitly beforehand.
//! These are geometric operations, not relaxation or chemical rotatability tests.
//!
//! See [`super`] for editing and scan examples.

use std::{
    collections::BTreeSet,
    f64::consts::{PI, TAU},
    fmt,
    sync::Arc,
};

use crate::geometry::{Point3, RigidTransform, Vector3};
use crate::topology::{AtomSelection, InstanceAtomId, Topology};
use crate::units::{Quantity, Unit, UnitError, CANONICAL_ANGLE_UNIT, CANONICAL_LENGTH_UNIT};

use super::{
    measure::{self, MeasurementError},
    Model, ModelView,
};

/// An edit failed before any coordinates were published.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum GeometryEditError {
    InvalidAtomId(InstanceAtomId),
    RepeatedAtom(InstanceAtomId),
    /// Selection or prepared edit belongs to a different shared snapshot.
    TopologyMismatch,
    MissingBond {
        a: InstanceAtomId,
        b: InstanceAtomId,
    },
    /// Removing the traversal bond did not separate moving and fixed references.
    InseparableFragment,
    /// The last reference must be selected and the first must be excluded.
    InvalidMovingSelection,
    Unit(UnitError),
    /// Nonfinite target, nonpositive distance, or angle outside [0, pi].
    InvalidTarget,
    /// A rigid-motion origin, axis, or displacement is nonfinite.
    NonFiniteInput,
    /// A direction, rotation plane, or dihedral is undefined.
    DegenerateGeometry,
    /// Arithmetic overflow or loss of precision prevents achieving the target.
    NumericalFailure,
}

impl fmt::Display for GeometryEditError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(atom) => write!(f, "invalid geometry-edit atom: {atom}"),
            Self::RepeatedAtom(atom) => write!(f, "repeated geometry-edit reference: {atom}"),
            Self::TopologyMismatch => f.write_str("geometry edit belongs to a different topology"),
            Self::MissingBond { a, b } => {
                write!(f, "geometry edit requires a bond between {a} and {b}")
            }
            Self::InseparableFragment => f.write_str("moving fragment reaches a fixed reference"),
            Self::InvalidMovingSelection => {
                f.write_str("select the last reference and exclude the first reference")
            }
            Self::Unit(error) => write!(f, "geometry-edit unit: {error}"),
            Self::InvalidTarget => {
                f.write_str("target must be finite, with distance > 0 or angle in [0, pi]")
            }
            Self::NonFiniteInput => f.write_str("rigid-motion arguments must be finite"),
            Self::DegenerateGeometry => {
                f.write_str("geometry edit has an undefined direction or plane")
            }
            Self::NumericalFailure => {
                f.write_str("geometry edit exceeds numerical range or precision")
            }
        }
    }
}

impl std::error::Error for GeometryEditError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

impl From<UnitError> for GeometryEditError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

impl From<MeasurementError> for GeometryEditError {
    fn from(error: MeasurementError) -> Self {
        match error {
            MeasurementError::InvalidAtomId(atom) => Self::InvalidAtomId(atom),
            MeasurementError::MissingBond { a, b } => Self::MissingBond { a, b },
            MeasurementError::Selection(_) => Self::TopologyMismatch,
            MeasurementError::Unit(error) => Self::Unit(error),
            MeasurementError::InvalidCutoff => Self::InvalidTarget,
            MeasurementError::DegenerateGeometry => Self::DegenerateGeometry,
            MeasurementError::NumericalFailure => Self::NumericalFailure,
        }
    }
}

// Private shared bookkeeping; only the three concrete edit types are public.
#[derive(Debug, Clone)]
struct Prepared<const N: usize> {
    atoms: [InstanceAtomId; N],
    indices: [usize; N],
    moving: AtomSelection,
}

impl<const N: usize> Prepared<N> {
    fn new(
        topology: &Arc<Topology>,
        atoms: [InstanceAtomId; N],
        moving: Option<&AtomSelection>,
    ) -> Result<Self, GeometryEditError> {
        let mut indices = [0; N];
        for (i, &atom) in atoms.iter().enumerate() {
            indices[i] = topology
                .atom_index(atom)
                .ok_or(GeometryEditError::InvalidAtomId(atom))?
                .index();
        }
        for (i, &atom) in atoms.iter().enumerate() {
            if atoms[..i].contains(&atom) {
                return Err(GeometryEditError::RepeatedAtom(atom));
            }
        }
        let moving = if let Some(moving) = moving {
            moving
                .ensure_compatible(topology)
                .map_err(|_| GeometryEditError::TopologyMismatch)?;
            if moving.contains(atoms[0]) || !moving.contains(atoms[N - 1]) {
                return Err(GeometryEditError::InvalidMovingSelection);
            }
            moving.clone()
        } else {
            for pair in atoms.windows(2) {
                if !topology
                    .neighbors(pair[0])
                    .map_err(|_| GeometryEditError::InvalidAtomId(pair[0]))?
                    .any(|a| a == pair[1])
                {
                    return Err(GeometryEditError::MissingBond {
                        a: pair[0],
                        b: pair[1],
                    });
                }
            }
            let (left, right) = if N == 2 {
                (atoms[0], atoms[1])
            } else {
                (atoms[1], atoms[2])
            };
            let mut visited = BTreeSet::from([right]);
            let mut pending = vec![right];
            while let Some(atom) = pending.pop() {
                for neighbor in topology
                    .neighbors(atom)
                    .map_err(|_| GeometryEditError::InvalidAtomId(atom))?
                {
                    if (atom == left && neighbor == right) || (atom == right && neighbor == left) {
                        continue;
                    }
                    if neighbor == left || neighbor == atoms[0] {
                        return Err(GeometryEditError::InseparableFragment);
                    }
                    if visited.insert(neighbor) {
                        pending.push(neighbor);
                    }
                }
            }
            AtomSelection::from_atoms(topology, visited)
                .map_err(|_| GeometryEditError::NumericalFailure)?
        };
        Ok(Self {
            atoms,
            indices,
            moving,
        })
    }

    fn ensure_compatible(&self, model: ModelView<'_>) -> Result<(), GeometryEditError> {
        self.moving
            .ensure_compatible(&model.shared_topology())
            .map_err(|_| GeometryEditError::TopologyMismatch)
    }

    fn points(&self, model: ModelView<'_>) -> [Point3; N] {
        let values = model.positions().values();
        self.indices.map(|i| values.value()[i])
    }

    fn value(points: [Point3; N]) -> Result<f64, GeometryEditError> {
        Ok(match N {
            2 => measure::norm(points[1] - points[0])?,
            3 => measure::angle_points(points[0], points[1], points[2])?,
            4 => measure::dihedral_points(points[0], points[1], points[2], points[3])?,
            _ => unreachable!("only distance, angle, and dihedral edits are constructed"),
        })
    }

    fn unit() -> Unit {
        if N == 2 {
            CANONICAL_LENGTH_UNIT
        } else {
            CANONICAL_ANGLE_UNIT
        }
    }

    fn measure(&self, model: ModelView<'_>) -> Result<Quantity<f64>, GeometryEditError> {
        self.ensure_compatible(model)?;
        Ok(Quantity::new(
            Self::value(self.points(model))?,
            Self::unit(),
        ))
    }

    fn apply(&self, model: &mut Model, target: Quantity<f64>) -> Result<(), GeometryEditError> {
        self.ensure_compatible(model.view())?;
        let mut target = target.into_unit(Self::unit())?.into_value();
        if !target.is_finite()
            || (N == 2 && target <= 0.0)
            || (N == 3 && !(0.0..=PI).contains(&target))
        {
            return Err(GeometryEditError::InvalidTarget);
        }
        if N == 4 {
            target = signed_angle(target);
        }
        let points = self.points(model.view());
        let current = Self::value(points)?;
        let delta = if N == 4 {
            signed_angle(target - current)
        } else {
            target - current
        };
        if delta == 0.0 {
            return Ok(());
        }
        let motion = match N {
            2 => Motion::Translation(measure::normalized(points[1] - points[0])? * delta),
            3 => {
                let u = measure::normalized(points[0] - points[1])?;
                let v = measure::normalized(points[2] - points[1])?;
                Motion::rotation(points[1], measure::normalized(u.cross(v))?, delta)
            }
            4 => Motion::rotation(
                points[1],
                measure::normalized(points[2] - points[1])?,
                delta,
            ),
            _ => unreachable!(),
        };
        // Pivots are kept bit-for-bit, even if included in the moving selection.
        let pivots = &self.indices[1..N - 1];
        let staged = stage(model, &self.moving, pivots, |p| motion.apply(p))?;
        let mut result = points;
        for (point, &index) in result.iter_mut().zip(&self.indices) {
            if let Ok(i) = staged.binary_search_by_key(&index, |&(index, _)| index) {
                *point = staged[i].1;
            }
        }
        // Reject finite but unachievable targets (e.g. displacement below the
        // coordinate resolution at a very large offset), rather than report success.
        let achieved = Self::value(result).map_err(|_| GeometryEditError::NumericalFailure)?;
        let residual = if N == 4 {
            signed_angle(achieved - target).abs()
        } else {
            (achieved - target).abs()
        };
        let tolerance = if N == 2 { 1.0e-10 * target } else { 1.0e-10 };
        if residual > tolerance {
            return Err(GeometryEditError::NumericalFailure);
        }
        publish(model, &staged)
    }
}

macro_rules! prepared_edit {
    ($name:ident, $n:literal, $docs:literal) => {
        #[doc = $docs]
        ///
        /// References and moving atoms bind to one topology layout. Each
        /// application uses current Cartesian coordinates and an absolute target.
        /// Failure is atomic. Properties, represented chemistry, and the cell are
        /// preserved. No periodic imaging or relaxation is performed.
        #[derive(Debug, Clone)]
        pub struct $name(Prepared<$n>);

        impl $name {
            /// Resolves the moving fragment through consecutive bonds. Rings
            /// that connect it to fixed references return `InseparableFragment`.
            /// Bond order does not restrict editing.
            pub fn new(
                topology: &Arc<Topology>,
                atoms: [InstanceAtomId; $n],
            ) -> Result<Self, GeometryEditError> {
                Prepared::new(topology, atoms, None).map(Self)
            }

            /// Moves exactly this selection (except defining pivots, which stay
            /// fixed). Select the last reference and exclude the first. Bonds
            /// are not required; boundary bonds and rings may deliberately deform.
            pub fn with_moving_atoms(
                topology: &Arc<Topology>,
                atoms: [InstanceAtomId; $n],
                moving: &AtomSelection,
            ) -> Result<Self, GeometryEditError> {
                Prepared::new(topology, atoms, Some(moving)).map(Self)
            }

            pub fn atoms(&self) -> [InstanceAtomId; $n] {
                self.0.atoms
            }

            /// Resolved selection, including any selected stationary pivots.
            pub fn moving_atoms(&self) -> &AtomSelection {
                &self.0.moving
            }

            /// Measures the current geometry in canonical units.
            pub fn measure(
                &self,
                model: ModelView<'_>,
            ) -> Result<Quantity<f64>, GeometryEditError> {
                self.0.measure(model)
            }

            /// Sets an absolute target atomically. Accepts compatible units.
            /// Finite staged coordinates must achieve the target within 1e-10
            /// relative error for distances or 1e-10 radians for angles.
            /// Otherwise returns `NumericalFailure` without changing the model.
            pub fn apply(
                &self,
                model: &mut Model,
                target: Quantity<f64>,
            ) -> Result<(), GeometryEditError> {
                self.0.apply(model, target)
            }
        }
    };
}

prepared_edit!(
    DistanceEdit,
    2,
    "An A–B distance edit: translate B's fragment along A→B, keeping A fixed.

Automatic movement cuts A–B. Targets must be finite and positive. A coincident
starting pair has no direction and is rejected."
);
prepared_edit!(
    AngleEdit,
    3,
    "An A–B–C angle edit: rotate C's fragment around B in the current plane.

Automatic movement cuts B–C. A and B stay fixed, including selected B. Targets
lie in [0, pi]. A collinear angle can remain unchanged, but a change needs a
defined rotation plane. Zero-length arms are always rejected."
);
prepared_edit!(
    DihedralEdit,
    4,
    "An A–B–C–D dihedral edit: rotate C's fragment around the B→C axis.

Automatic movement cuts B–C. A, B, and C stay fixed, including selected pivots.
Finite targets wrap modulo 2*pi, following the measurement sign convention.
The shortest rotation is used; an exact half-turn tie chooses +pi. Degenerate
defining triples are rejected even when the requested target is unchanged."
);

/// Maps to (-pi, pi], without adding pi to a potentially huge input.
fn signed_angle(angle: f64) -> f64 {
    let wrapped = angle % TAU;
    if wrapped <= -PI {
        wrapped + TAU
    } else if wrapped > PI {
        wrapped - TAU
    } else {
        wrapped
    }
}

enum Motion {
    Translation(Vector3),
    Rotation {
        origin: Point3,
        axis: Vector3,
        sin: f64,
        cos: f64,
    },
}

impl Motion {
    fn rotation(origin: Point3, axis: Vector3, angle: f64) -> Self {
        let (sin, cos) = signed_angle(angle).sin_cos();
        Self::Rotation {
            origin,
            axis,
            sin,
            cos,
        }
    }

    fn apply(&self, point: Point3) -> Point3 {
        match *self {
            Self::Translation(displacement) => point + displacement,
            Self::Rotation {
                origin,
                axis,
                sin,
                cos,
            } => {
                if sin == 0.0 && cos == 1.0 {
                    return point;
                }
                let v = point - origin;
                origin + v * cos + axis.cross(v) * sin + axis * (axis.dot(v) * (1.0 - cos))
            }
        }
    }
}

fn stage(
    model: &Model,
    selection: &AtomSelection,
    pivots: &[usize],
    transform: impl Fn(Point3) -> Point3,
) -> Result<Vec<(usize, Point3)>, GeometryEditError> {
    selection
        .ensure_compatible(&model.shared_topology())
        .map_err(|_| GeometryEditError::TopologyMismatch)?;
    let values = model.positions().values();
    selection
        .indices()
        .iter()
        .map(|index| index.index())
        .filter(|index| !pivots.contains(index))
        .map(|index| {
            let point = transform(values.value()[index]);
            if !point.is_finite() {
                return Err(GeometryEditError::NumericalFailure);
            }
            Ok((index, point))
        })
        .collect()
}

fn publish(model: &mut Model, staged: &[(usize, Point3)]) -> Result<(), GeometryEditError> {
    model
        .positions
        .set_canonical_batch(staged)
        .map_err(|_| GeometryEditError::NumericalFailure)
}

impl Model {
    /// Sets A–B by translating B's connected fragment. See [`DistanceEdit`].
    pub fn set_distance(
        &mut self,
        a: InstanceAtomId,
        b: InstanceAtomId,
        target: Quantity<f64>,
    ) -> Result<(), GeometryEditError> {
        DistanceEdit::new(&self.shared_topology(), [a, b])?.apply(self, target)
    }

    /// Sets A–B–C by rotating C's connected fragment. See [`AngleEdit`].
    pub fn set_angle(
        &mut self,
        a: InstanceAtomId,
        b: InstanceAtomId,
        c: InstanceAtomId,
        target: Quantity<f64>,
    ) -> Result<(), GeometryEditError> {
        AngleEdit::new(&self.shared_topology(), [a, b, c])?.apply(self, target)
    }

    /// Sets A–B–C–D around B→C. See [`DihedralEdit`] for sign and moving atoms.
    pub fn set_dihedral(
        &mut self,
        a: InstanceAtomId,
        b: InstanceAtomId,
        c: InstanceAtomId,
        d: InstanceAtomId,
        target: Quantity<f64>,
    ) -> Result<(), GeometryEditError> {
        DihedralEdit::new(&self.shared_topology(), [a, b, c, d])?.apply(self, target)
    }

    /// Translates selected Cartesian coordinates atomically, preserving the cell,
    /// topology, and properties. Even empty selections validate units and binding.
    pub fn translate(
        &mut self,
        selection: &AtomSelection,
        displacement: Quantity<Vector3>,
    ) -> Result<(), GeometryEditError> {
        let displacement = displacement.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
        if !displacement.is_finite() {
            return Err(GeometryEditError::NonFiniteInput);
        }
        let staged = stage(self, selection, &[], |point| point + displacement)?;
        publish(self, &staged)
    }

    /// Rotates selected coordinates right-handed around a finite origin and a
    /// finite, nonzero axis (normalized internally). The cell is unchanged even
    /// for an all-atom selection. Empty selections still validate all arguments.
    pub fn rotate(
        &mut self,
        selection: &AtomSelection,
        origin: Quantity<Point3>,
        axis: Vector3,
        angle: Quantity<f64>,
    ) -> Result<(), GeometryEditError> {
        let origin = origin.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
        let angle = angle.into_unit(CANONICAL_ANGLE_UNIT)?.into_value();
        if !origin.is_finite() || !axis.is_finite() {
            return Err(GeometryEditError::NonFiniteInput);
        }
        if !angle.is_finite() {
            return Err(GeometryEditError::InvalidTarget);
        }
        // Scale first so any finite nonzero axis, including subnormal and very
        // large components, can be normalized without norm overflow.
        let scale = axis.x.abs().max(axis.y.abs()).max(axis.z.abs());
        if scale == 0.0 {
            return Err(GeometryEditError::DegenerateGeometry);
        }
        let axis = measure::normalized(axis / scale)?;
        let motion = Motion::rotation(origin, axis, angle);
        let staged = stage(self, selection, &[], |point| motion.apply(point))?;
        publish(self, &staged)
    }

    /// Applies a validated rigid transform atomically to a selection. Translation
    /// is in canonical length units, as in [`RigidTransform`]. The cell, topology,
    /// and properties are preserved, including for an all-atom selection.
    pub fn apply_transform(
        &mut self,
        selection: &AtomSelection,
        transform: &RigidTransform,
    ) -> Result<(), GeometryEditError> {
        let staged = stage(self, selection, &[], |point| {
            transform.transform_point(point)
        })?;
        publish(self, &staged)
    }
}

#[cfg(test)]
mod tests;
