use kekule::geometry::Vector3;
use kekule::units::{Quantity, CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT};
use kekule_traj::{Forces, Velocities};

#[test]
fn vector_arrays_construct_from_dense_values_without_topology_imports() {
    let values = [Vector3::new(1.0, 2.0, 3.0)];
    let velocities = Velocities::new(Quantity::new(values, CANONICAL_VELOCITY_UNIT)).unwrap();
    let forces = Forces::new(Quantity::new(values, CANONICAL_FORCE_UNIT)).unwrap();
    assert_eq!(velocities.len(), 1);
    assert_eq!(forces.len(), 1);
}
