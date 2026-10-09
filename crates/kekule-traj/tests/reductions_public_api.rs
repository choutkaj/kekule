use std::sync::Arc;

use kekule::structure::{Trajectory, TrajectoryFrame};
use kekule::{
    geometry::{PeriodicCell, Point3, Vector3},
    structure::Positions,
    topology::{AtomSelection, InstanceAtomId, SelectionError, Topology},
    units::{Quantity, ANGSTROM, NANOMETER, PICOSECOND},
};
use kekule_traj::analysis::{self, ContactOccupancyAccumulator, ReductionError, RmsfAccumulator};
use kekule_traj::{MemoryTrajectoryReader, TrajectoryReader};

mod support;
use support::linear_carbon_topology;

fn trajectory(top: &Arc<Topology>, rows: &[[f64; 3]]) -> Trajectory {
    Trajectory::from_items(
        top.clone(),
        rows.iter().map(|row| {
            TrajectoryFrame::new(
                Positions::new(Quantity::new(
                    row.map(|x| Point3::new(x, 0.0, 0.0)),
                    NANOMETER,
                ))
                .unwrap(),
            )
        }),
    )
    .unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-8, "{a} != {b}");
}

#[test]
fn loaded_and_reused_streaming_buffers_agree_with_independent_statistics() {
    let top = linear_carbon_topology(3);
    let ids = top.atom_ids();
    let pairs = [(ids[0], ids[1]), (ids[1], ids[2])];
    let all = AtomSelection::all(&top);
    // Large common translation exposes cancellation in E[x²] - E[x]².
    let rows = [
        [1.0e9, 1.0e9 + 1.0, 1.0e9 + 20.0],
        [1.0e9 + 2.0, 1.0e9 + 5.0, 1.0e9 + 20.0],
        [1.0e9 + 4.0, 1.0e9 + 6.0, 1.0e9 + 20.0],
    ];
    let t = trajectory(&top, &rows);
    let before = format!("{t:?}");
    let rmsf = analysis::rmsf(&t, &all).unwrap();
    let contacts = analysis::contact_occupancy(&t, pairs, Quantity::new(2.0, NANOMETER)).unwrap();
    assert_eq!(rmsf.frame_count(), 3);
    assert_eq!(rmsf.selection(), &all);
    assert_eq!(contacts.pairs(), &pairs);
    assert!(std::ptr::eq(contacts.topology(), top.as_ref()));
    assert_eq!(contacts.hit_counts(), &[2, 0]);
    close(contacts.contacts().next().unwrap().1, 2.0 / 3.0);
    close(contacts.cutoff().value_in(ANGSTROM).unwrap(), 20.0);
    for (index, (atom, value)) in rmsf.atoms().enumerate() {
        assert_eq!(atom, ids[index]);
        let mut pairwise = 0.0;
        for a in &rows {
            for b in &rows {
                pairwise += (a[index] - b[index]).powi(2);
            }
        }
        close(value.into_value().powi(2), pairwise / (2.0 * 9.0));
    }
    close(rmsf.values().value()[0], (8.0_f64 / 3.0).sqrt());
    close(rmsf.values().value()[1], (14.0_f64 / 3.0).sqrt());
    close(rmsf.values().value()[2], 0.0);
    let subset = analysis::rmsf(
        &t,
        &AtomSelection::from_atoms(&top, [ids[2], ids[0]]).unwrap(),
    )
    .unwrap();
    assert_eq!(
        subset.atoms().map(|(atom, _)| atom).collect::<Vec<_>>(),
        [ids[0], ids[2]]
    );
    assert_eq!(subset.values().value(), &[rmsf.values().value()[0], 0.0]);
    close(
        rmsf.values().value_in(ANGSTROM).unwrap()[0],
        (8.0_f64 / 3.0).sqrt() * 10.0,
    );
    let mut online_rmsf = RmsfAccumulator::new(&all).unwrap();
    let mut online_contacts =
        ContactOccupancyAccumulator::new(&top, pairs, Quantity::new(2.0, NANOMETER)).unwrap();
    let mut reader = MemoryTrajectoryReader::new(&t);
    let mut buffer = reader.frame_buffer();
    let mut index = 0;
    while reader.read_next(&mut buffer).unwrap() {
        online_rmsf.observe(index, &buffer).unwrap();
        online_contacts.observe(index, &buffer).unwrap();
        index += 1;
    }
    assert_eq!(index, 3);
    assert_eq!(online_rmsf.finish().unwrap().values(), rmsf.values());
    assert_eq!(
        online_contacts.finish().unwrap().hit_counts(),
        contacts.hit_counts()
    );
    assert_eq!(
        analysis::contact_occupancy(&t, pairs, Quantity::new(21.0, ANGSTROM))
            .unwrap()
            .hit_counts(),
        &[2, 0]
    );
    assert_eq!(format!("{t:?}"), before);
}

#[test]
fn failed_frames_are_atomic_and_retryable_for_both_accumulators() {
    let top = linear_carbon_topology(3);
    let ids = top.atom_ids();
    let pairs = [(ids[0], ids[1]), (ids[1], ids[2])];
    let good = trajectory(&top, &[[0.0, 1.0, 0.0]]);
    // RMSF can fail after updating the first pending atom; contact calculation
    // can fail after the first pending pair. Neither may publish that prefix.
    let bad = trajectory(&top, &[[0.5, 1.0e308, -1.0e308]]);
    let unrelated = trajectory(&linear_carbon_topology(3), &[[0.0, 1.0, 0.0]]);
    let mut rmsf = RmsfAccumulator::new(&AtomSelection::all(&top)).unwrap();
    let mut contacts =
        ContactOccupancyAccumulator::new(&top, pairs, Quantity::new(2.0, NANOMETER)).unwrap();
    rmsf.observe(10, &good.get(0).unwrap()).unwrap();
    contacts.observe(10, &good.get(0).unwrap()).unwrap();
    assert_eq!(
        rmsf.observe(11, &unrelated.get(0).unwrap()),
        Err(ReductionError::TopologyMismatch { frame: 11 })
    );
    assert_eq!(
        contacts.observe(11, &unrelated.get(0).unwrap()),
        Err(ReductionError::TopologyMismatch { frame: 11 })
    );
    assert!(matches!(
        rmsf.observe(11, &bad.get(0).unwrap()),
        Err(ReductionError::NumericalFailure { frame: 11, .. })
    ));
    assert!(matches!(
        contacts.observe(11, &bad.get(0).unwrap()),
        Err(ReductionError::NumericalFailure { frame: 11, .. })
    ));
    assert_eq!(rmsf.frame_count(), 1);
    assert_eq!(contacts.frame_count(), 1);
    rmsf.observe(11, &good.get(0).unwrap()).unwrap();
    contacts.observe(11, &good.get(0).unwrap()).unwrap();
    let rmsf = rmsf.finish().unwrap();
    let contacts = contacts.finish().unwrap();
    assert_eq!(rmsf.frame_count(), 2);
    assert_eq!(rmsf.values().value(), &[0.0, 0.0, 0.0]);
    assert_eq!(contacts.frame_count(), 2);
    assert_eq!(contacts.hit_counts(), &[2, 2]);
}

#[test]
fn reducers_reject_empty_or_ambiguous_inputs_and_handle_single_frames() {
    let top = linear_carbon_topology(3);
    let ids = top.atom_ids();
    let all = AtomSelection::all(&top);
    let empty = AtomSelection::from_atoms(&top, []).unwrap();
    let cutoff = Quantity::new(1.0, NANOMETER);
    assert!(matches!(
        RmsfAccumulator::new(&empty),
        Err(ReductionError::EmptySelection)
    ));
    assert!(matches!(
        RmsfAccumulator::new(&all).unwrap().finish(),
        Err(ReductionError::NoFrames)
    ));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [], cutoff),
        Err(ReductionError::EmptySelection)
    ));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], ids[1])], cutoff)
            .unwrap()
            .finish(),
        Err(ReductionError::NoFrames)
    ));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], ids[0])], cutoff),
        Err(ReductionError::SelfPair { pair: 0 })
    ));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], ids[1]), (ids[1], ids[0])], cutoff),
        Err(ReductionError::DuplicatePair { pair: 1 })
    ));
    let invalid = InstanceAtomId::new(ids[0].molecule(), kekule::core::AtomId::new(99));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], invalid)], cutoff),
        Err(ReductionError::Selection(SelectionError::InvalidAtomId(_)))
    ));
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            ContactOccupancyAccumulator::new(
                &top,
                [(ids[0], ids[1])],
                Quantity::new(value, NANOMETER)
            ),
            Err(ReductionError::InvalidCutoff)
        ));
    }
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], ids[1])], Quantity::new(1.0, PICOSECOND)),
        Err(ReductionError::Unit(_))
    ));
    let none = Trajectory::new(top.clone());
    assert!(matches!(
        analysis::rmsf(&none, &all),
        Err(ReductionError::NoFrames)
    ));
    assert!(matches!(
        analysis::contact_occupancy(&none, [(ids[0], ids[1])], cutoff),
        Err(ReductionError::NoFrames)
    ));
    let one = trajectory(&top, &[[0.0, 0.0, 1.0]]);
    assert_eq!(
        analysis::rmsf(&one, &all).unwrap().values().value(),
        &[0.0, 0.0, 0.0]
    );
    assert_eq!(
        analysis::contact_occupancy(&one, [(ids[0], ids[1])], Quantity::new(0.0, NANOMETER))
            .unwrap()
            .hit_counts(),
        &[1]
    );
    let foreign = AtomSelection::all(&linear_carbon_topology(3));
    assert!(matches!(
        analysis::rmsf(&one, &foreign),
        Err(ReductionError::Selection(SelectionError::TopologyMismatch))
    ));
}

#[test]
fn stored_coordinate_statistics_do_not_implicitly_image_or_align_periodic_frames() {
    let top = linear_carbon_topology(3);
    let ids = top.atom_ids();
    let mut t = trajectory(&top, &[[0.0, 9.0, 2.0], [1.0, 10.0, 3.0]]);
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    for index in 0..t.len() {
        t.get_mut(index)
            .unwrap()
            .conformation_mut()
            .set_cell(Some(cell));
    }
    // A minimum-image contact would match; stored distance is nine nanometers.
    assert_eq!(
        analysis::contact_occupancy(&t, [(ids[0], ids[1])], Quantity::new(2.0, NANOMETER))
            .unwrap()
            .hit_counts(),
        &[0]
    );
    for value in analysis::rmsf(&t, &AtomSelection::all(&top))
        .unwrap()
        .values()
        .value()
    {
        close(*value, 0.5);
    }
}
