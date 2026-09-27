use kekule::{
    geometry::Point3,
    smiles,
    structure::{
        measure::{self, BondDihedral},
        Model, Positions,
    },
    units::{Quantity, NANOMETER},
};
use kekule_traj::{Trajectory, TrajectoryFrame};

#[test]
fn trajectory_frames_keep_the_same_tied_references_through_degenerate_geometry() {
    let mut molecule = smiles::to_molecules("CC(C)C(C)C").unwrap().remove(0);
    molecule.perceive().unwrap();
    let points = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::origin(),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 1.0),
        Point3::new(1.0, 1.0, 0.0),
    ];
    let model = Model::from_molecule(
        &molecule,
        &Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap();
    let topology = model.shared_topology();
    let axis = model.bond_ids()[2];
    let definition = BondDihedral::new(&topology, axis).unwrap();
    let atoms = definition.atoms().unwrap();
    assert_eq!(
        atoms,
        [
            model.atom_ids()[0],
            model.atom_ids()[1],
            model.atom_ids()[3],
            model.atom_ids()[4]
        ]
    );
    let mut degenerate = points;
    degenerate[0] = Point3::new(0.0, -1.0, 0.0);
    let mut rotated = points;
    rotated[0] = Point3::new(-1.0, 0.0, 0.0);
    let frames = [points, degenerate, rotated].map(|points| {
        TrajectoryFrame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap())
    });
    let trajectory = Trajectory::from_frames(topology, frames).unwrap();
    let angles = trajectory
        .frames()
        .map(|frame| {
            let view = frame.as_model();
            let prepared = definition.measure(view).unwrap();
            assert_eq!(prepared, measure::bond_dihedral(view, axis).unwrap());
            assert_eq!(
                BondDihedral::new(&view.shared_topology(), axis)
                    .unwrap()
                    .atoms(),
                Some(atoms)
            );
            prepared.map(|q| q.into_value())
        })
        .collect::<Vec<_>>();
    assert!((angles[0].unwrap() + std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    assert_eq!(angles[1], None);
    assert!((angles[2].unwrap() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
}
