//! Dependency-light three-dimensional geometry shared by structure, analysis,
//! trajectory, and modelling code.
//!
//! [`Point3`], [`Vector3`] and [`Matrix3`] are bare Cartesian values; wrap them
//! in a [`Quantity`] to give them a unit. [`RigidTransform`] is a validated
//! proper rotation plus translation, [`PeriodicCell`] a validated simulation
//! cell, and [`PeriodicGeometry`] its prepared lattice for fractional
//! coordinates and exact nearest-image searches. Cell and lattice outputs use
//! the canonical length unit (nm).
//!
//! # Example
//!
//! ```
//! use kekule::geometry::{PeriodicCell, PeriodicGeometry, Point3, Vector3};
//! use kekule::units::{Quantity, NANOMETER};
//!
//! let a = Point3::new(0.5, 0.0, 0.0);
//! let b = Point3::new(1.75, 0.0, 0.0);
//! assert_eq!(b - a, Vector3::new(1.25, 0.0, 0.0));
//!
//! let cell = PeriodicCell::orthorhombic(
//!     Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
//!     [true; 3],
//! )?;
//! let lattice = PeriodicGeometry::new(cell)?;
//! let (image, shift) = lattice.nearest_image(Quantity::new(b - a, NANOMETER))?;
//! assert_eq!(shift, [1, 0, 0]);
//! assert_eq!(image, Quantity::new(Vector3::new(-0.75, 0.0, 0.0), NANOMETER));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod grid;
mod periodic;
pub(crate) use grid::PointGrid;
pub use periodic::{PeriodicGeometry, PeriodicGeometryError, MAX_LATTICE_IMAGES};

use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use crate::units::{Quantity, ScaleValue, UnitError, CANONICAL_LENGTH_UNIT, CUBIC_NANOMETER};

/// A Cartesian point.
///
/// The coordinates are bare numbers; a [`Quantity<Point3>`](Quantity) carries
/// their length unit.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Point3 {
    /// The x coordinate.
    pub x: f64,
    /// The y coordinate.
    pub y: f64,
    /// The z coordinate.
    pub z: f64,
}

impl Point3 {
    /// The point with coordinates `(x, y, z)`.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The point `(0, 0, 0)`.
    pub const fn origin() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Whether every coordinate is finite (neither NaN nor infinite).
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
}

impl crate::units::sealed::Sealed for Point3 {}

impl ScaleValue for Point3 {
    fn scaled(self, factor: f64) -> Self {
        Self::new(self.x * factor, self.y * factor, self.z * factor)
    }
}

/// A Cartesian vector, such as the displacement between two [`Point3`]s.
///
/// The components are bare numbers; a [`Quantity<Vector3>`](Quantity) carries
/// their unit.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Vector3 {
    /// The x component.
    pub x: f64,
    /// The y component.
    pub y: f64,
    /// The z component.
    pub z: f64,
}

impl Vector3 {
    /// The vector with components `(x, y, z)`.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The zero vector.
    pub const fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Whether every component is finite (neither NaN nor infinite).
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// The dot product `self · other`.
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// The right-handed cross product `self × other`.
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// The squared Euclidean magnitude, `self · self`.
    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean magnitude without intermediate squaring overflow or underflow.
    pub fn norm(self) -> f64 {
        self.x.hypot(self.y).hypot(self.z)
    }
}

impl crate::units::sealed::Sealed for Vector3 {}

impl ScaleValue for Vector3 {
    fn scaled(self, factor: f64) -> Self {
        Self::new(self.x * factor, self.y * factor, self.z * factor)
    }
}

impl Add<Vector3> for Point3 {
    type Output = Self;

    fn add(self, vector: Vector3) -> Self::Output {
        Self::new(self.x + vector.x, self.y + vector.y, self.z + vector.z)
    }
}

impl Sub<Vector3> for Point3 {
    type Output = Self;

    fn sub(self, vector: Vector3) -> Self::Output {
        Self::new(self.x - vector.x, self.y - vector.y, self.z - vector.z)
    }
}

impl Sub for Point3 {
    type Output = Vector3;

    fn sub(self, other: Self) -> Self::Output {
        Vector3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl Add for Vector3 {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl AddAssign for Vector3 {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Sub for Vector3 {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl SubAssign for Vector3 {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Neg for Vector3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl Mul<f64> for Vector3 {
    type Output = Self;

    fn mul(self, scale: f64) -> Self::Output {
        Self::new(self.x * scale, self.y * scale, self.z * scale)
    }
}

impl Mul<Vector3> for f64 {
    type Output = Vector3;

    fn mul(self, vector: Vector3) -> Self::Output {
        vector * self
    }
}

impl Div<f64> for Vector3 {
    type Output = Self;

    fn div(self, scale: f64) -> Self::Output {
        Self::new(self.x / scale, self.y / scale, self.z / scale)
    }
}

/// A three-by-three matrix stored by columns.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix3 {
    columns: [Vector3; 3],
}

impl Matrix3 {
    /// The matrix whose columns are `a`, `b` and `c`, in order.
    pub const fn from_columns(a: Vector3, b: Vector3, c: Vector3) -> Self {
        Self { columns: [a, b, c] }
    }

    /// The identity matrix.
    pub const fn identity() -> Self {
        Self::from_columns(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        )
    }

    /// The three columns, in order.
    pub const fn columns(self) -> [Vector3; 3] {
        self.columns
    }

    /// Whether every entry is finite (neither NaN nor infinite).
    pub fn is_finite(self) -> bool {
        self.columns.into_iter().all(Vector3::is_finite)
    }

    /// The determinant, computed as the triple product `a · (b × c)` of the
    /// columns.
    pub fn determinant(self) -> f64 {
        self.columns[0].dot(self.columns[1].cross(self.columns[2]))
    }

    /// The transposed matrix.
    pub fn transpose(self) -> Self {
        let [a, b, c] = self.columns;
        Self::from_columns(
            Vector3::new(a.x, b.x, c.x),
            Vector3::new(a.y, b.y, c.y),
            Vector3::new(a.z, b.z, c.z),
        )
    }

    /// The matrix-vector product `self * vector`.
    pub fn transform_vector(self, vector: Vector3) -> Vector3 {
        self.columns[0] * vector.x + self.columns[1] * vector.y + self.columns[2] * vector.z
    }

    /// The matrix product `self * other`, which applies `other` first.
    pub fn multiply(self, other: Self) -> Self {
        let [a, b, c] = other.columns;
        Self::from_columns(
            self.transform_vector(a),
            self.transform_vector(b),
            self.transform_vector(c),
        )
    }
}

impl Default for Matrix3 {
    fn default() -> Self {
        Self::identity()
    }
}

impl crate::units::sealed::Sealed for Matrix3 {}

impl ScaleValue for Matrix3 {
    fn scaled(self, factor: f64) -> Self {
        let [a, b, c] = self.columns;
        Self::from_columns(a.scaled(factor), b.scaled(factor), c.scaled(factor))
    }
}

/// A validated periodic simulation cell in Kekule's canonical length unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeriodicCell {
    vectors: Matrix3,
    periodic_axes: [bool; 3],
}

impl PeriodicCell {
    /// Creates an orthorhombic cell from positive axis lengths, with the cell
    /// vectors along x, y and z.
    ///
    /// `periodic_axes[i]` marks whether the cell repeats along vector `i`.
    ///
    /// # Errors
    ///
    /// [`PeriodicCellError::Unit`] when `lengths` is not a length,
    /// [`PeriodicCellError::InvalidOrthorhombicLengths`] when a length is not
    /// finite and strictly positive, and [`PeriodicCellError::NoPeriodicAxes`]
    /// when no axis is periodic.
    pub fn orthorhombic(
        lengths: Quantity<Vector3>,
        periodic_axes: [bool; 3],
    ) -> Result<Self, PeriodicCellError> {
        let lengths = lengths.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
        if !lengths.is_finite() || lengths.x <= 0.0 || lengths.y <= 0.0 || lengths.z <= 0.0 {
            return Err(PeriodicCellError::InvalidOrthorhombicLengths);
        }
        Self::new(
            Quantity::new(
                [
                    Vector3::new(lengths.x, 0.0, 0.0),
                    Vector3::new(0.0, lengths.y, 0.0),
                    Vector3::new(0.0, 0.0, lengths.z),
                ],
                CANONICAL_LENGTH_UNIT,
            ),
            periodic_axes,
        )
    }

    /// Creates an orthorhombic or triclinic cell from three basis vectors.
    ///
    /// `periodic_axes[i]` marks whether the cell repeats along `vectors[i]`.
    ///
    /// # Errors
    ///
    /// [`PeriodicCellError::NoPeriodicAxes`] when no axis is periodic,
    /// [`PeriodicCellError::Unit`] when `vectors` is not a length,
    /// [`PeriodicCellError::NonFiniteVector`] when a component is not finite,
    /// and [`PeriodicCellError::DegenerateVectors`] when the vectors span zero
    /// or non-finite volume.
    pub fn new(
        vectors: Quantity<[Vector3; 3]>,
        periodic_axes: [bool; 3],
    ) -> Result<Self, PeriodicCellError> {
        if !periodic_axes.into_iter().any(|periodic| periodic) {
            return Err(PeriodicCellError::NoPeriodicAxes);
        }
        let [a, b, c] = vectors.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
        let matrix = Matrix3::from_columns(a, b, c);
        if !matrix.is_finite() {
            return Err(PeriodicCellError::NonFiniteVector);
        }
        let volume = matrix.determinant();
        if !volume.is_finite() || volume.abs() <= f64::EPSILON {
            return Err(PeriodicCellError::DegenerateVectors);
        }
        Ok(Self {
            vectors: matrix,
            periodic_axes,
        })
    }

    /// The three cell vectors, in the canonical length unit (nm).
    pub fn vectors(self) -> Quantity<[Vector3; 3]> {
        Quantity::new(self.vectors.columns(), CANONICAL_LENGTH_UNIT)
    }

    /// Whether the cell repeats along each of its three vectors; at least one
    /// is `true`.
    pub const fn periodic_axes(self) -> [bool; 3] {
        self.periodic_axes
    }

    /// The signed volume `a · (b × c)` in nm³; negative for a left-handed
    /// basis.
    pub fn signed_volume(self) -> Quantity<f64> {
        Quantity::new(self.vectors.determinant(), CUBIC_NANOMETER)
    }

    /// Whether the cell vectors lie exactly along x, y and z, in that order.
    pub fn is_orthorhombic(self) -> bool {
        let [a, b, c] = self.vectors.columns();
        a.y == 0.0 && a.z == 0.0 && b.x == 0.0 && b.z == 0.0 && c.x == 0.0 && c.y == 0.0
    }
}

/// A validated rotation and translation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidTransform {
    rotation: Matrix3,
    translation: Vector3,
}

impl RigidTransform {
    /// The transform that leaves every point unchanged.
    pub const fn identity() -> Self {
        Self {
            rotation: Matrix3::identity(),
            translation: Vector3::zero(),
        }
    }

    /// A rotation followed by a translation of any length unit.
    ///
    /// The rotation must be orthonormal with determinant +1, within an
    /// absolute tolerance of 1e-10 per entry.
    ///
    /// # Errors
    ///
    /// [`RigidTransformError::Unit`] when `translation` is not a length,
    /// [`RigidTransformError::NonFinite`] when a component is not finite, and
    /// [`RigidTransformError::NonRigidRotation`] when `rotation` is not a
    /// proper rotation.
    pub fn new(
        rotation: Matrix3,
        translation: Quantity<Vector3>,
    ) -> Result<Self, RigidTransformError> {
        let translation = translation
            .value_in(CANONICAL_LENGTH_UNIT)
            .map_err(|error| RigidTransformError::Unit(Box::new(error)))?;
        Self::from_canonical(rotation, translation)
    }

    /// [`Self::new`] for a translation already in the canonical length unit.
    pub(crate) fn from_canonical(
        rotation: Matrix3,
        translation: Vector3,
    ) -> Result<Self, RigidTransformError> {
        if !rotation.is_finite() || !translation.is_finite() {
            return Err(RigidTransformError::NonFinite);
        }
        let gram = rotation.transpose().multiply(rotation);
        let [a, b, c] = gram.columns();
        let tolerance = 1.0e-10;
        let orthonormal = (a.x - 1.0).abs() <= tolerance
            && a.y.abs() <= tolerance
            && a.z.abs() <= tolerance
            && b.x.abs() <= tolerance
            && (b.y - 1.0).abs() <= tolerance
            && b.z.abs() <= tolerance
            && c.x.abs() <= tolerance
            && c.y.abs() <= tolerance
            && (c.z - 1.0).abs() <= tolerance
            && (rotation.determinant() - 1.0).abs() <= tolerance;
        if !orthonormal {
            return Err(RigidTransformError::NonRigidRotation);
        }
        Ok(Self {
            rotation,
            translation,
        })
    }

    /// The proper rotation, applied before the translation.
    pub const fn rotation(self) -> Matrix3 {
        self.rotation
    }

    /// The translation, in the canonical length unit (nm).
    pub fn translation(self) -> Quantity<Vector3> {
        Quantity::new(self.translation, CANONICAL_LENGTH_UNIT)
    }

    /// Rotates `point` about the origin, then translates it.
    ///
    /// `point` is taken to be in the canonical length unit (nm), the unit of
    /// [`Self::translation`].
    pub fn transform_point(self, point: Point3) -> Point3 {
        let rotated = self
            .rotation
            .transform_vector(Vector3::new(point.x, point.y, point.z));
        Point3::new(
            rotated.x + self.translation.x,
            rotated.y + self.translation.y,
            rotated.z + self.translation.z,
        )
    }

    /// Rotates `vector`; a displacement is unaffected by the translation.
    pub fn transform_vector(self, vector: Vector3) -> Vector3 {
        self.rotation.transform_vector(vector)
    }
}

impl Default for RigidTransform {
    fn default() -> Self {
        Self::identity()
    }
}

/// Cell vectors or axis flags do not describe a valid periodic cell.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum PeriodicCellError {
    /// The vectors could not be converted to the canonical length unit; see the
    /// source.
    Unit(UnitError),
    /// No axis is marked periodic.
    NoPeriodicAxes,
    /// A cell vector component is not finite.
    NonFiniteVector,
    /// The cell vectors span zero or non-finite volume.
    DegenerateVectors,
    /// An orthorhombic length is not finite and strictly positive.
    InvalidOrthorhombicLengths,
}

impl fmt::Display for PeriodicCellError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unit(_) => formatter.write_str("invalid periodic-cell unit"),
            Self::NoPeriodicAxes => formatter.write_str("periodic cell has no periodic axis"),
            Self::NonFiniteVector => formatter.write_str("periodic-cell vectors must be finite"),
            Self::DegenerateVectors => {
                formatter.write_str("periodic-cell vectors must have non-zero finite volume")
            }
            Self::InvalidOrthorhombicLengths => formatter
                .write_str("orthorhombic cell lengths must be finite and strictly positive"),
        }
    }
}

impl std::error::Error for PeriodicCellError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

impl From<UnitError> for PeriodicCellError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

/// A rotation and translation do not form a proper rigid transform.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum RigidTransformError {
    /// The translation is not a length; see the source.
    Unit(Box<UnitError>),
    /// A rotation or translation component is not finite.
    NonFinite,
    /// The rotation matrix is not orthonormal with determinant +1.
    NonRigidRotation,
}

impl fmt::Display for RigidTransformError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unit(_) => formatter.write_str("rigid translation is not a length"),
            Self::NonFinite => formatter.write_str("rigid transform must be finite"),
            Self::NonRigidRotation => {
                formatter.write_str("rotation matrix must be right-handed and orthonormal")
            }
        }
    }
}

impl std::error::Error for RigidTransformError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error.as_ref()),
            Self::NonFinite | Self::NonRigidRotation => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{ANGSTROM, NANOMETER};

    #[test]
    fn point_and_vector_arithmetic_preserves_meaning() {
        let a = Point3::new(1.0, 2.0, 3.0);
        let b = Point3::new(4.0, 6.0, 8.0);
        let displacement = b - a;
        assert_eq!(displacement, Vector3::new(3.0, 4.0, 5.0));
        assert_eq!(a + displacement, b);
        assert_eq!(b - displacement, a);
    }

    #[test]
    fn periodic_cells_validate_units_axes_and_volume() {
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(1.0, 2.0, 3.0), NANOMETER),
            [true; 3],
        )
        .unwrap();
        assert!(cell.is_orthorhombic());
        assert_eq!(
            cell.vectors().value(),
            &[
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, 2.0, 0.0),
                Vector3::new(0.0, 0.0, 3.0)
            ]
        );
        assert_eq!(cell.vectors().unit(), NANOMETER);
        assert_eq!(
            PeriodicCell::new(
                Quantity::new(
                    [
                        Vector3::new(1.0, 0.0, 0.0),
                        Vector3::new(2.0, 0.0, 0.0),
                        Vector3::new(0.0, 0.0, 1.0)
                    ],
                    ANGSTROM
                ),
                [true; 3]
            ),
            Err(PeriodicCellError::DegenerateVectors)
        );
    }

    #[test]
    fn rigid_transform_rejects_scaling() {
        let scaled = Matrix3::from_columns(
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        );
        assert_eq!(
            RigidTransform::new(scaled, Quantity::new(Vector3::zero(), NANOMETER)),
            Err(RigidTransformError::NonRigidRotation)
        );
    }
}
