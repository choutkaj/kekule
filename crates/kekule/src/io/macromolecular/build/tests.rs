use super::super::connectivity::FALLBACK_COVALENT_RADIUS_ANGSTROM;
use super::*;
use crate::molecule::BondOrder;

/// The grid counts the unbonded close pairs that measuring every pair
/// counts, for seeded atoms on cell boundaries, exactly at bonding
/// distance, and some bonded; too-large coordinates still fail.
#[test]
fn bond_candidates_match_every_pair() {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut below = move |bound: u64| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) % bound
    };
    let palette = ["C", "N", "O", "S", "Fe", "Xe"]
        .map(|symbol| Element::from_symbol(symbol).expect("element symbol"));
    for spread in [5u64, 9, 30] {
        let mut editor = MoleculeEditor::new();
        let mut sites = Vec::<(AtomId, Element, Point3, usize)>::new();
        for line in 0..400 {
            let element = palette[below(palette.len() as u64) as usize];
            let atom = editor.add_atom(Atom::new(element)).unwrap();
            let point = if line % 6 == 1 {
                // On the 2.1 Å cell edge, or at the largest cutoff.
                let previous = sites[line - 1].2;
                Point3::new(
                    previous.x + COVALENT_BOND_CELL_ANGSTROM,
                    previous.y,
                    previous.z,
                )
            } else if line % 6 == 4 {
                let mut axis = || below(spread * 2) as f64 * 1.05;
                Point3::new(axis(), axis(), axis())
            } else {
                let mut axis = || below(spread * 1000) as f64 * 1e-3;
                Point3::new(axis(), axis(), axis())
            };
            if line > 0 && line % 3 == 0 {
                editor
                    .add_bond(sites[line - 1].0, atom, BondOrder::Single)
                    .unwrap();
            }
            sites.push((atom, element, point, line));
        }
        let graph = editor.graph();
        let mut expected = 0;
        for left in 0..sites.len() {
            for right in left + 1..sites.len() {
                let (left, left_element, left_point, _) = sites[left];
                let (right, right_element, right_point, _) = sites[right];
                let distance_squared = point_distance_squared(left_point, right_point);
                let radius = |element: Element| {
                    element
                        .covalent_radius_angstrom()
                        .unwrap_or(FALLBACK_COVALENT_RADIUS_ANGSTROM)
                };
                let cutoff = (radius(left_element)
                    + radius(right_element)
                    + COVALENT_BOND_TOLERANCE_ANGSTROM)
                    .min(COVALENT_BOND_CELL_ANGSTROM);
                if distance_squared > MIN_COVALENT_BOND_DISTANCE_SQUARED
                    && distance_squared <= cutoff * cutoff
                    && graph.bond_between(left, right).is_none()
                {
                    expected += 1;
                }
            }
        }
        assert!(expected > 0, "{spread}");
        assert_eq!(
            count_bond_candidates(graph, &sites, None).unwrap(),
            expected,
            "{spread}"
        );
        sites[7].2.y = 1.0e300;
        assert!(count_bond_candidates(graph, &sites, None).is_err());
    }
}
