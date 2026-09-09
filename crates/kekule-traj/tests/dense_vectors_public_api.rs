use kekule::geometry::Vector3;
use kekule::units::{Quantity, CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT};
use kekule_traj::{Forces, FrameBuffer, FrameError, Velocities};

mod support;

struct ChangingSlice {
    calls: std::cell::Cell<usize>,
    valid: [Vector3; 1],
    invalid: [Vector3; 1],
}

impl ChangingSlice {
    fn new() -> Self {
        Self {
            calls: std::cell::Cell::new(0),
            valid: [Vector3::new(1.0, 2.0, 3.0)],
            invalid: [Vector3::new(f64::NAN, 0.0, 0.0)],
        }
    }
}

impl AsRef<[Vector3]> for ChangingSlice {
    fn as_ref(&self) -> &[Vector3] {
        let previous = self.calls.replace(self.calls.get() + 1);
        if previous == 0 {
            &self.valid
        } else {
            &self.invalid
        }
    }
}

#[test]
fn vector_arrays_construct_from_dense_values_without_topology_imports() {
    let values = [Vector3::new(1.0, 2.0, 3.0)];
    let velocities = Velocities::new(Quantity::new(values, CANONICAL_VELOCITY_UNIT)).unwrap();
    let forces = Forces::new(Quantity::new(values, CANONICAL_FORCE_UNIT)).unwrap();
    assert_eq!(velocities.len(), 1);
    assert_eq!(forces.len(), 1);
}

#[test]
fn vector_replacement_validates_and_copies_one_borrowed_slice() {
    let mut velocities = Velocities::zeros(1);
    let mut forces = Forces::zeros(1);
    velocities
        .set_all(Quantity::new(ChangingSlice::new(), CANONICAL_VELOCITY_UNIT))
        .unwrap();
    forces
        .set_all(Quantity::new(ChangingSlice::new(), CANONICAL_FORCE_UNIT))
        .unwrap();
    assert_eq!(velocities.values().value(), &[Vector3::new(1.0, 2.0, 3.0)]);
    assert_eq!(forces.values().value(), &[Vector3::new(1.0, 2.0, 3.0)]);
    let mut buffer = FrameBuffer::new(support::linear_carbon_topology(1));
    buffer
        .set_velocities(Some(Quantity::new(
            ChangingSlice::new(),
            CANONICAL_VELOCITY_UNIT,
        )))
        .unwrap();
    buffer
        .set_forces(Some(Quantity::new(
            ChangingSlice::new(),
            CANONICAL_FORCE_UNIT,
        )))
        .unwrap();
    assert_eq!(
        buffer.frame_view().velocities().unwrap().value(),
        velocities.values().value()
    );
    assert_eq!(
        buffer.frame_view().forces().unwrap().value(),
        forces.values().value()
    );
}

#[test]
fn owned_vector_construction_keeps_allocations_and_validates_converted_values() {
    let values = vec![Vector3::new(1.0, 2.0, 3.0)];
    let pointer = values.as_ptr();
    let velocities = Velocities::from_vec(Quantity::new(values, CANONICAL_VELOCITY_UNIT)).unwrap();
    assert_eq!(velocities.values().value().as_ptr(), pointer);
    let values = velocities.into_values();
    assert_eq!(values.value().as_ptr(), pointer);
    assert_eq!(values.unit(), CANONICAL_VELOCITY_UNIT);
    let values = vec![Vector3::new(1.0, 2.0, 3.0)];
    let pointer = values.as_ptr();
    let source_unit = kekule::units::KILOJOULE_PER_MOLE / kekule::units::ANGSTROM;
    let forces = Forces::from_vec(Quantity::new(values, source_unit)).unwrap();
    assert_eq!(forces.values().value().as_ptr(), pointer);
    assert_eq!(forces.values().value(), &[Vector3::new(10.0, 20.0, 30.0)]);
    let values = forces.into_values();
    assert_eq!(values.value().as_ptr(), pointer);
    assert_eq!(values.unit(), CANONICAL_FORCE_UNIT);
    assert!(matches!(
        Velocities::from_vec(Quantity::new(
            vec![Vector3::new(f64::NAN, 0.0, 0.0)],
            CANONICAL_VELOCITY_UNIT
        )),
        Err(FrameError::NonFiniteVector { index: 0 })
    ));
    assert!(matches!(
        Forces::from_vec(Quantity::new(
            vec![Vector3::new(f64::MAX, 0.0, 0.0)],
            source_unit
        )),
        Err(FrameError::NonFiniteVector { index: 0 })
    ));
    assert!(matches!(
        Forces::from_vec(Quantity::new(
            vec![Vector3::zero()],
            CANONICAL_VELOCITY_UNIT
        )),
        Err(FrameError::Unit(_))
    ));
}

#[test]
fn buffer_vector_clear_methods_need_no_type_annotations_and_preserve_other_state() {
    let mut buffer = FrameBuffer::new(support::linear_carbon_topology(1));
    buffer.set_step(Some(42));
    buffer
        .set_velocities(Some(Quantity::new(
            [Vector3::new(1.0, 2.0, 3.0)],
            CANONICAL_VELOCITY_UNIT,
        )))
        .unwrap();
    buffer
        .set_forces(Some(Quantity::new(
            [Vector3::new(4.0, 5.0, 6.0)],
            CANONICAL_FORCE_UNIT,
        )))
        .unwrap();
    let velocity_pointer = buffer.frame_view().velocities().unwrap().value().as_ptr();
    let force_pointer = buffer.frame_view().forces().unwrap().value().as_ptr();
    buffer.clear_velocities();
    assert!(buffer.frame_view().velocities().is_none());
    assert!(buffer.frame_view().forces().is_some());
    buffer.clear_forces();
    assert!(buffer.frame_view().forces().is_none());
    assert_eq!(buffer.frame_view().step(), Some(42));
    buffer
        .set_velocities(Some(Quantity::new(
            [Vector3::zero()],
            CANONICAL_VELOCITY_UNIT,
        )))
        .unwrap();
    buffer
        .set_forces(Some(Quantity::new([Vector3::zero()], CANONICAL_FORCE_UNIT)))
        .unwrap();
    assert_eq!(
        buffer.frame_view().velocities().unwrap().value().as_ptr(),
        velocity_pointer
    );
    assert_eq!(
        buffer.frame_view().forces().unwrap().value().as_ptr(),
        force_pointer
    );
}
