//! Coordinate files address topology atoms in source atom-row order.

use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::mmcif;
use kekule::structure::Positions;
use kekule::units::{Quantity, ANGSTROM};
use kekule_traj::io::{read_trajectory, write_trajectory};
use kekule_traj::{Trajectory, TrajectoryFrame};

// A covalent link joins the first and third atom-site rows into one molecule,
// so instance-major numbering would move the water between the carbons.
const LINKED_ROWS: &str = r#"
data_linked_rows
loop_
_entity.id
_entity.type
1 non-polymer
2 water
3 non-polymer
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
C 3
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 LIG A 1 . 0.0 0.0 0.0
HETATM 2 O O HOH B 2 . 10.0 0.0 0.0
HETATM 3 C C1 LG2 C 3 . 1.5 0.0 0.0
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
link-1 covale A C1 C C1
"#;

const ROW_ELEMENTS: [&str; 3] = ["C", "O", "C"];
const ROW_X: [f64; 3] = [0.0, 10.0, 1.5];

struct TemporaryPath(PathBuf);

impl Drop for TemporaryPath {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

fn row_x(points: &[Point3]) -> Vec<f64> {
    points
        .iter()
        .map(|point| (point.x * 1.0e4).round() / 1.0e3)
        .collect()
}

#[test]
fn mmcif_topology_reads_external_coordinates_in_atom_site_order() {
    let interpretation =
        mmcif::interpret(&mmcif::parse_str(LINKED_ROWS).unwrap(), Default::default()).unwrap();
    let topology = interpretation.model().shared_topology();
    assert_eq!(topology.instance_count(), 2);
    let instances = topology
        .atom_ids()
        .iter()
        .map(|atom| atom.molecule())
        .collect::<Vec<_>>();
    assert_eq!(instances[0], instances[2]);
    assert_ne!(instances[0], instances[1]);
    assert_eq!(
        topology
            .atoms()
            .map(|(_, atom)| atom.element.symbol())
            .collect::<Vec<_>>(),
        ROW_ELEMENTS
    );
    assert_eq!(
        row_x(interpretation.model().positions().values().value()),
        ROW_X
    );

    // An engine writes coordinates in the structure file's row order. DCD
    // carries no atom identity, so write it through an unrelated topology.
    let external = Arc::new(kekule::smiles::to_topology("C.O.C").unwrap());
    let rows = ROW_X.map(|x| Point3::new(x, 0.0, 0.0)).to_vec();
    let mut frame = TrajectoryFrame::new(Positions::new(Quantity::new(rows, ANGSTROM)).unwrap());
    frame.set_step(Some(0));
    let path = TemporaryPath(
        std::env::temp_dir().join(format!("kekule-atom-order-{}.dcd", std::process::id())),
    );
    write_trajectory(
        &path.0,
        &Trajectory::from_frames(external, [frame]).unwrap(),
    )
    .unwrap();

    let loaded = read_trajectory(&path.0, Arc::clone(&topology)).unwrap();
    let frame = loaded.frame(0).unwrap().as_model();
    assert_eq!(row_x(frame.positions().values().value()), ROW_X);
    for (atom, element) in topology.atom_ids().iter().zip(ROW_ELEMENTS) {
        assert_eq!(frame.atom(*atom).unwrap().element.symbol(), element);
    }
}
