use kekule::{
    math::Point3,
    smiles,
    structure::{
        measure::{self, BondDihedral},
        Positions, Structure,
    },
    units::{Quantity, NANOMETER},
};

use crate::support::{time_ordered, Frame};

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
    let structure = Structure::from_molecule(
        molecule.clone(),
        &Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap();
    let topology = structure.shared_topology();
    let axis = structure.topology().bond_ids()[2];
    let definition = BondDihedral::new(topology.clone(), axis, Default::default()).unwrap();
    let atoms = definition.atoms().unwrap();
    assert_eq!(
        atoms,
        [
            structure.topology().atom_ids()[0],
            structure.topology().atom_ids()[1],
            structure.topology().atom_ids()[3],
            structure.topology().atom_ids()[4]
        ]
    );
    let mut degenerate = points;
    degenerate[0] = Point3::new(0.0, -1.0, 0.0);
    let mut rotated = points;
    rotated[0] = Point3::new(-1.0, 0.0, 0.0);
    let frames = [points, degenerate, rotated]
        .map(|points| Frame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap()));
    let trajectory = time_ordered(topology, frames);
    let angles = trajectory
        .iter()
        .map(|frame| {
            let view = frame.as_structure_view();
            let prepared = definition.measure(view).unwrap();
            assert_eq!(prepared, measure::bond_dihedral(view, axis).unwrap());
            assert_eq!(
                BondDihedral::new(view.shared_topology(), axis, Default::default())
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
