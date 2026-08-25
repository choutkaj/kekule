//! Dependency-light three-dimensional geometry shared by structure, analysis,
//! trajectory, and modelling code.

use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use crate::units::{Quantity, ScaleValue, UnitError, CANONICAL_LENGTH_UNIT};

/// A Cartesian point.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub const fn origin() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
}

impl ScaleValue for Point3 {
    fn scaled(self, factor: f64) -> Self {
        Self::new(self.x * factor, self.y * factor, self.z * factor)
    }
}

/// A Cartesian vector.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub const fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }

    pub fn norm(self) -> f64 {
        self.norm_squared().sqrt()
    }

    pub(crate) fn add_scaled(&mut self, other: Self, scale: f64) {
        self.x += other.x * scale;
        self.y += other.y * scale;
        self.z += other.z * scale;
    }
}

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
    pub const fn from_columns(a: Vector3, b: Vector3, c: Vector3) -> Self {
        Self { columns: [a, b, c] }
    }

    pub const fn identity() -> Self {
        Self::from_columns(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        )
    }

    pub const fn columns(self) -> [Vector3; 3] {
        self.columns
    }

    pub fn is_finite(self) -> bool {
        self.columns.into_iter().all(Vector3::is_finite)
    }

    pub fn determinant(self) -> f64 {
        self.columns[0].dot(self.columns[1].cross(self.columns[2]))
    }

    pub fn transpose(self) -> Self {
        let [a, b, c] = self.columns;
        Self::from_columns(
            Vector3::new(a.x, b.x, c.x),
            Vector3::new(a.y, b.y, c.y),
            Vector3::new(a.z, b.z, c.z),
        )
    }

    pub fn transform_vector(self, vector: Vector3) -> Vector3 {
        self.columns[0] * vector.x + self.columns[1] * vector.y + self.columns[2] * vector.z
    }

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
    /// Creates an orthorhombic cell from positive axis lengths.
    pub fn orthorhombic(
        lengths: Quantity<Vector3>,
        periodic_axes: [bool; 3],
    ) -> Result<Self, PeriodicCellError> {
        let lengths = lengths.to_unit(CANONICAL_LENGTH_UNIT)?.to_value();
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
    pub fn new(
        vectors: Quantity<[Vector3; 3]>,
        periodic_axes: [bool; 3],
    ) -> Result<Self, PeriodicCellError> {
        if !periodic_axes.into_iter().any(|periodic| periodic) {
            return Err(PeriodicCellError::NoPeriodicAxes);
        }
        let [a, b, c] = vectors.to_unit(CANONICAL_LENGTH_UNIT)?.to_value();
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

    pub fn vectors(self) -> Quantity<[Vector3; 3]> {
        Quantity::new(self.vectors.columns(), CANONICAL_LENGTH_UNIT)
    }

    pub const fn periodic_axes(self) -> [bool; 3] {
        self.periodic_axes
    }

    pub fn signed_volume(self) -> f64 {
        self.vectors.determinant()
    }

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
    pub const fn identity() -> Self {
        Self {
            rotation: Matrix3::identity(),
            translation: Vector3::zero(),
        }
    }

    pub fn new(rotation: Matrix3, translation: Vector3) -> Result<Self, RigidTransformError> {
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

    pub const fn rotation(self) -> Matrix3 {
        self.rotation
    }

    pub const fn translation(self) -> Vector3 {
        self.translation
    }

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

    pub fn transform_vector(self, vector: Vector3) -> Vector3 {
        self.rotation.transform_vector(vector)
    }
}

impl Default for RigidTransform {
    fn default() -> Self {
        Self::identity()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum PeriodicCellError {
    Unit(UnitError),
    NoPeriodicAxes,
    NonFiniteVector,
    DegenerateVectors,
    InvalidOrthorhombicLengths,
}

impl fmt::Display for PeriodicCellError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unit(error) => write!(formatter, "invalid periodic-cell unit: {error}"),
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

impl std::error::Error for PeriodicCellError {}

impl From<UnitError> for PeriodicCellError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RigidTransformError {
    NonFinite,
    NonRigidRotation,
}

impl fmt::Display for RigidTransformError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite => formatter.write_str("rigid transform must be finite"),
            Self::NonRigidRotation => {
                formatter.write_str("rotation matrix must be right-handed and orthonormal")
            }
        }
    }
}

impl std::error::Error for RigidTransformError {}

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
            RigidTransform::new(scaled, Vector3::zero()),
            Err(RigidTransformError::NonRigidRotation)
        );
    }
}
