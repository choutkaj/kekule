//! mmCIF `_cell` parameters and Kekule periodic cells, both in the standard
//! crystallographic frame mmCIF Cartesian coordinates use: a along x, b in
//! the xy plane, c completing a right-handed basis.
use crate::geometry::{PeriodicCell, PeriodicCellError, Vector3};
use crate::units::{Quantity, ANGSTROM};

/// `_cell` lengths in angstroms and angles in degrees, in the order a, b, c,
/// alpha, beta, gamma.
pub(crate) type CellParameters = [f64; 6];

/// Exact for a right angle, so a rectangular cell stays orthorhombic.
fn cos_degrees(angle: f64) -> f64 {
    if angle == 90.0 {
        0.0
    } else {
        angle.to_radians().cos()
    }
}

fn sin_degrees(angle: f64) -> f64 {
    if angle == 90.0 {
        1.0
    } else {
        angle.to_radians().sin()
    }
}

/// The cell, periodic along all three axes, whose basis has these parameters
/// in the standard frame.
pub(crate) fn cell_from_parameters(
    [a, b, c, alpha, beta, gamma]: CellParameters,
) -> Result<PeriodicCell, PeriodicCellError> {
    let (cos_alpha, cos_beta, cos_gamma) =
        (cos_degrees(alpha), cos_degrees(beta), cos_degrees(gamma));
    let sin_gamma = sin_degrees(gamma);
    let cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
    let cz = (1.0 - cos_beta * cos_beta - cy * cy).sqrt();
    PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(a, 0.0, 0.0),
                Vector3::new(b * cos_gamma, b * sin_gamma, 0.0),
                Vector3::new(c * cos_beta, c * cy, c * cz),
            ],
            ANGSTROM,
        ),
        [true; 3],
    )
}

/// Why a cell has no `_cell` parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CellParameterError {
    /// mmCIF cells are periodic along all three axes.
    PartiallyPeriodic,
    /// The basis is not in the standard frame, so lengths and angles alone
    /// would not place it relative to the coordinates.
    NonstandardOrientation,
}

/// The parameters of a fully periodic cell in the standard frame.
pub(crate) fn cell_parameters(cell: &PeriodicCell) -> Result<CellParameters, CellParameterError> {
    if cell.periodic_axes() != [true; 3] {
        return Err(CellParameterError::PartiallyPeriodic);
    }
    let [a, b, c] = cell
        .vectors()
        .into_unit(ANGSTROM)
        .expect("cell vectors are lengths")
        .into_value();
    if a.y != 0.0 || a.z != 0.0 || b.z != 0.0 || a.x <= 0.0 || b.y <= 0.0 || c.z <= 0.0 {
        return Err(CellParameterError::NonstandardOrientation);
    }
    let angle = |u: Vector3, v: Vector3| {
        let dot = u.dot(v);
        if dot == 0.0 {
            90.0
        } else {
            (dot / (u.norm() * v.norm()))
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees()
        }
    };
    Ok([
        a.norm(),
        b.norm(),
        c.norm(),
        angle(b, c),
        angle(a, c),
        angle(a, b),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cell_parameters_round_trip_through_the_standard_frame() {
        for expected in [
            [20.0, 30.0, 40.0, 90.0, 90.0, 90.0],
            [30.0, 30.0, 30.0, 60.0, 60.0, 90.0],
        ] {
            let cell = cell_from_parameters(expected).unwrap();
            assert_eq!(cell.is_orthorhombic(), expected[3..] == [90.0; 3]);
            let actual = cell_parameters(&cell).unwrap();
            // Lengths pass through the canonical nanometre unit.
            for (actual, expected) in actual.into_iter().zip(expected) {
                assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
            }
        }
        let mut partial = cell_from_parameters([20.0, 30.0, 40.0, 90.0, 90.0, 90.0]).unwrap();
        partial = PeriodicCell::new(partial.vectors(), [true, true, false]).unwrap();
        assert_eq!(
            cell_parameters(&partial),
            Err(CellParameterError::PartiallyPeriodic)
        );
        let rotated = PeriodicCell::new(
            Quantity::new(
                [
                    Vector3::new(0.0, 2.0, 0.0),
                    Vector3::new(-2.0, 0.0, 0.0),
                    Vector3::new(0.0, 0.0, 2.0),
                ],
                ANGSTROM,
            ),
            [true; 3],
        )
        .unwrap();
        assert_eq!(
            cell_parameters(&rotated),
            Err(CellParameterError::NonstandardOrientation)
        );
    }
}
