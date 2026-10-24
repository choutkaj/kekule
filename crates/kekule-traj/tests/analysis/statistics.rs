use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::statistics::{ContactOccupancyAccumulator, RmsfAccumulator};
use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
use kekule::topology::AtomSelection;
use kekule::units::{Quantity, NANOMETER};
use kekule_traj::FrameBuffer;

use crate::support::linear_carbon_topology;

#[test]
fn reused_streaming_buffers_agree_with_loaded_trajectory_statistics() {
    let top = linear_carbon_topology(3);
    let ids = top.atom_ids();
    let pairs = [(ids[0], ids[1]), (ids[1], ids[2])];
    let all = AtomSelection::all(Arc::clone(&top));
    let cutoff = Quantity::new(2.0, NANOMETER);
    let rows = [[0.0, 1.0, 20.0], [2.0, 5.0, 20.0], [4.0, 6.0, 20.0]];
    let t = Trajectory::from_items(
        top.clone(),
        rows.iter().map(|row| {
            let points = row.map(|x| Point3::new(x, 0.0, 0.0));
            TrajectoryFrame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap())
        }),
    )
    .unwrap();

    let mut rmsf = RmsfAccumulator::new(&all).unwrap();
    let mut contacts = ContactOccupancyAccumulator::new(&top, pairs, cutoff).unwrap();
    let mut buffer = FrameBuffer::new(t.shared_topology());
    for (index, frame) in t.iter().enumerate() {
        buffer.copy_from(frame).unwrap();
        rmsf.observe(index, &buffer.as_model_view(), 1.0).unwrap();
        contacts
            .observe(index, &buffer.as_model_view(), 1.0)
            .unwrap();
    }
    assert_eq!(
        rmsf.finish().unwrap().values(),
        t.rmsf(&all).unwrap().values()
    );
    let contacts = contacts.finish().unwrap();
    assert_eq!(contacts.contact_weights(), &[2.0, 0.0]);
    assert_eq!(
        contacts.contact_weights(),
        t.contact_occupancy(pairs, cutoff)
            .unwrap()
            .contact_weights()
    );
}
