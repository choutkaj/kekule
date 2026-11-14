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
fn cell_parameters_round_trip_through_the_standard_orientation() {
    let parameters = |values: [f64; 6]| CellParameters {
        lengths: Quantity::new([values[0], values[1], values[2]], ANGSTROM),
        angles: Quantity::new([values[3], values[4], values[5]], DEGREE),
    };
    for expected in [
        [20.0, 30.0, 40.0, 90.0, 90.0, 90.0],
        [30.0, 30.0, 30.0, 60.0, 60.0, 90.0],
    ] {
        let cell = PeriodicCell::from_parameters(parameters(expected)).unwrap();
        assert_eq!(cell.is_orthorhombic(), expected[3..] == [90.0; 3]);
        let actual = cell.to_parameters().unwrap();
        assert_eq!(actual.angles.unit(), DEGREE);
        let lengths = actual.lengths.into_value();
        let angles = actual.angles.into_value();
        // Lengths pass through the canonical nanometre unit.
        for (actual, expected) in lengths.into_iter().chain(angles).zip(expected) {
            assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
        }
    }
    let box_ =
        PeriodicCell::from_parameters(parameters([20.0, 30.0, 40.0, 90.0, 90.0, 90.0])).unwrap();
    let partial = PeriodicCell::new(box_.vectors(), [true, true, false]).unwrap();
    assert_eq!(
        partial.to_parameters(),
        Err(PeriodicCellError::PartiallyPeriodic)
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
        rotated.to_parameters(),
        Err(PeriodicCellError::NonstandardOrientation)
    );
    for invalid in [
        [-20.0, 30.0, 40.0, 90.0, 90.0, 90.0],
        [20.0, 0.0, 40.0, 90.0, 90.0, 90.0],
        [20.0, 30.0, f64::NAN, 90.0, 90.0, 90.0],
        [20.0, 30.0, 40.0, 0.0, 90.0, 90.0],
        [20.0, 30.0, 40.0, 90.0, 180.0, 90.0],
        [20.0, 30.0, 40.0, 90.0, 90.0, 270.0],
    ] {
        assert_eq!(
            PeriodicCell::from_parameters(parameters(invalid)),
            Err(PeriodicCellError::InvalidParameters),
            "{invalid:?}"
        );
    }
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
