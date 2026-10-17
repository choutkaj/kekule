//! Exact lattice geometry. Inputs and Cartesian outputs use canonical length units.
use super::{PeriodicCell, Vector3};
use crate::diagnostics::{LimitExceeded, Resource};
use std::fmt;

/// Maximum lattice images an exact nearest-image search examines.
pub const MAX_LATTICE_IMAGES: usize = 1_000_000;

/// A nonfinite, unrepresentable, or excessively expensive periodic calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PeriodicGeometryError {
    /// The calculation exceeds numerical precision or range.
    NumericalFailure,
    /// More than [`MAX_LATTICE_IMAGES`] lattice images would need examination.
    LimitExceeded(LimitExceeded),
}
impl fmt::Display for PeriodicGeometryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NumericalFailure => {
                f.write_str("periodic calculation exceeds numerical precision or range")
            }
            Self::LimitExceeded(limit) => write!(f, "nearest-image search {limit}"),
        }
    }
}

/// The nearest-image search would examine `images` lattice images, more than
/// [`MAX_LATTICE_IMAGES`]; `None` means more than `u64` counts.
fn image_limit(images: Option<u64>) -> PeriodicGeometryError {
    let observed = images
        .and_then(|images| usize::try_from(images).ok())
        .unwrap_or(usize::MAX);
    PeriodicGeometryError::LimitExceeded(LimitExceeded::new(
        Resource::LatticeImages,
        MAX_LATTICE_IMAGES,
        observed.max(MAX_LATTICE_IMAGES + 1),
    ))
}
impl std::error::Error for PeriodicGeometryError {}

/// Prepared reciprocal basis and exact bounded nearest-image search.
///
/// Supports rotated, skewed, left-handed, and partially periodic cells. Ties are
/// deterministic. Cartesian quantities use canonical length units (nm); fractional
/// coordinates and image indices are dimensionless. No approximate image fallback
/// is used when the one-million-candidate bound is exceeded.
#[derive(Debug, Clone)]
pub struct PeriodicGeometry {
    basis: [Vector3; 3],
    reciprocal: [Vector3; 3],
    periodic: [bool; 3],
}

impl PeriodicGeometry {
    pub fn new(cell: PeriodicCell) -> Result<Self, PeriodicGeometryError> {
        let [a, b, c] = cell.vectors().into_value();
        let inverse_volume = 1.0 / a.dot(b.cross(c));
        let reciprocal = [
            b.cross(c) * inverse_volume,
            c.cross(a) * inverse_volume,
            a.cross(b) * inverse_volume,
        ];
        if !reciprocal.iter().all(|v| v.is_finite()) {
            return Err(PeriodicGeometryError::NumericalFailure);
        }
        Ok(Self {
            basis: [a, b, c],
            reciprocal,
            periodic: cell.periodic_axes(),
        })
    }

    pub fn fractional(&self, vector: Vector3) -> Result<[f64; 3], PeriodicGeometryError> {
        let fractional = self.reciprocal.map(|dual| dual.dot(vector));
        if !fractional.iter().all(|v| v.is_finite()) {
            return Err(PeriodicGeometryError::NumericalFailure);
        }
        Ok(fractional)
    }

    pub fn cartesian(&self, fractional: [f64; 3]) -> Vector3 {
        self.basis[0] * fractional[0]
            + self.basis[1] * fractional[1]
            + self.basis[2] * fractional[2]
    }

    pub fn minimum_image(&self, delta: Vector3) -> Result<Vector3, PeriodicGeometryError> {
        self.nearest_image(delta).map(|(vector, _)| vector)
    }

    pub fn nearest_image(
        &self,
        delta: Vector3,
    ) -> Result<(Vector3, [i64; 3]), PeriodicGeometryError> {
        let fractional = self.fractional(delta)?;
        let mut images = [0.0; 3];
        for axis in 0..3 {
            if self.periodic[axis] {
                images[axis] = checked_image(fractional[axis])?;
            }
        }
        let mut best = delta - self.cartesian(images);
        let mut best_images = images.map(|n| n as i64);
        let mut best_squared = best.norm_squared();
        if !best_squared.is_finite() {
            return Err(PeriodicGeometryError::NumericalFailure);
        }
        // In an exactly orthogonal basis the three image choices are independent.
        // Besides avoiding enumeration, this handles very elongated rectangular
        // cells without an unnecessarily large Cartesian search radius.
        if self.basis[0].dot(self.basis[1]) == 0.0
            && self.basis[0].dot(self.basis[2]) == 0.0
            && self.basis[1].dot(self.basis[2]) == 0.0
        {
            return Ok((best, best_images));
        }
        // If r bounds the best Cartesian residual, reciprocal-vector norms bound
        // every fractional residual. Enumerating this finite box therefore includes
        // the true closest lattice image, even for skewed or partially periodic cells.
        let radius = best_squared.sqrt();
        let mut bounds = [(0_i64, 0_i64); 3];
        let mut candidates = 1_u64;
        for axis in 0..3 {
            if self.periodic[axis] {
                let span = radius * self.reciprocal[axis].norm();
                let padding = 64.0 * f64::EPSILON * (1.0 + fractional[axis].abs() + span);
                let low = (fractional[axis] - span - padding).ceil();
                let high = (fractional[axis] + span + padding).floor();
                if !low.is_finite()
                    || !high.is_finite()
                    || low.abs().max(high.abs()) >= 2.0_f64.powi(52)
                {
                    return Err(PeriodicGeometryError::NumericalFailure);
                }
                let (low, high) = (low as i64, high as i64);
                let count = (high - low + 1).max(0) as u64;
                let next = candidates.checked_mul(count);
                candidates = next
                    .filter(|n| *n <= MAX_LATTICE_IMAGES as u64)
                    .ok_or_else(|| image_limit(next))?;
                bounds[axis] = (low, high);
            }
        }
        for a in bounds[0].0..=bounds[0].1 {
            for b in bounds[1].0..=bounds[1].1 {
                for c in bounds[2].0..=bounds[2].1 {
                    let residual = delta - self.cartesian([a as f64, b as f64, c as f64]);
                    let squared = residual.norm_squared();
                    if squared < best_squared {
                        best = residual;
                        best_squared = squared;
                        best_images = [a, b, c];
                    }
                }
            }
        }
        Ok((best, best_images))
    }
}

fn checked_image(value: f64) -> Result<f64, PeriodicGeometryError> {
    if !value.is_finite() || value.abs() >= 2.0_f64.powi(52) {
        return Err(PeriodicGeometryError::NumericalFailure);
    }
    Ok(value.round())
}

impl PeriodicGeometry {
    /// Reciprocal basis lengths, for conservative fractional neighbor bounds.
    pub(crate) fn reciprocal_norms(&self) -> [f64; 3] {
        self.reciprocal.map(|v| v.norm())
    }
}
impl PeriodicGeometry {
    /// Shortest nonzero lattice translation along periodic axes, in canonical units.
    /// Rejects searches requiring more than one million integer candidates.
    pub fn shortest_translation(&self) -> Result<Vector3, PeriodicGeometryError> {
        let mut best = self
            .basis
            .iter()
            .zip(self.periodic)
            .filter(|(_, p)| *p)
            .map(|(v, _)| *v)
            .min_by(|a, b| a.norm().total_cmp(&b.norm()))
            .expect("at least one periodic axis");
        if self.basis[0].dot(self.basis[1]) == 0.0
            && self.basis[0].dot(self.basis[2]) == 0.0
            && self.basis[1].dot(self.basis[2]) == 0.0
        {
            return Ok(best);
        }
        let radius = best.norm();
        let mut bounds = [0_i64; 3];
        let mut count = 1_u64;
        for (i, bound) in bounds.iter_mut().enumerate() {
            if self.periodic[i] {
                let extent =
                    (radius * self.reciprocal[i].norm() * (1.0 + 64.0 * f64::EPSILON)).ceil();
                if !extent.is_finite() || extent > 500_000.0 {
                    // Along this axis alone, 2 * extent + 1 images exceed the
                    // bound; `as` saturates an infinite extent.
                    return Err(image_limit(Some((2.0 * extent + 1.0) as u64)));
                }
                *bound = extent as i64;
                let next = count.checked_mul(2 * (*bound as u64) + 1);
                count = next
                    .filter(|n| *n <= MAX_LATTICE_IMAGES as u64)
                    .ok_or_else(|| image_limit(next))?;
            }
        }
        for a in -bounds[0]..=bounds[0] {
            for b in -bounds[1]..=bounds[1] {
                for c in -bounds[2]..=bounds[2] {
                    if [a, b, c] == [0; 3] {
                        continue;
                    }
                    let candidate = self.cartesian([a as f64, b as f64, c as f64]);
                    if candidate.norm() < best.norm() {
                        best = candidate;
                    }
                }
            }
        }
        if !best.is_finite() || best.norm() == 0.0 {
            return Err(PeriodicGeometryError::NumericalFailure);
        }
        Ok(best)
    }
}
