//! Molfile stereo marks and the chiral flag at interpretation.

use kekule::molfile::{
    self, MolfileInterpretIssue, MolfileInterpretation, UnrepresentableStereoReason,
};
use kekule::stereo::StereoElementKind;

/// A V2000 record of `atoms` (symbol, x, y) and `bonds` (from, to, order,
/// stereo code), one-based, with the counts-line chiral flag `chiral`.
fn v2000(atoms: &[(&str, f64, f64)], bonds: &[(usize, usize, u8, u8)], chiral: bool) -> String {
    let mut text = format!(
        "marks\nkekule\n\n{:>3}{:>3}  0  0{:>3}  0            999 V2000\n",
        atoms.len(),
        bonds.len(),
        u8::from(chiral)
    );
    for (symbol, x, y) in atoms {
        text.push_str(&format!(
            "{x:>10.4}{y:>10.4}{:>10.4} {symbol:<3} 0  0  0  0  0  0\n",
            0.0
        ));
    }
    for (from, to, order, stereo) in bonds {
        text.push_str(&format!("{from:>3}{to:>3}{order:>3}{stereo:>3}  0  0  0\n"));
    }
    text.push_str("M  END\n");
    text
}

fn interpret(text: &str) -> MolfileInterpretation {
    molfile::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap_or_else(|error| panic!("{error}\n{text}"))
}

/// Bromochlorofluoromethane with a wedge from carbon to fluorine.
fn wedged_center(chiral: bool) -> String {
    v2000(
        &[
            ("C", 0.0, 0.0),
            ("F", 0.0, 1.0),
            ("Cl", 0.87, -0.5),
            ("Br", -0.87, -0.5),
        ],
        &[(1, 2, 1, 1), (1, 3, 1, 0), (1, 4, 1, 0)],
        chiral,
    )
}

#[test]
fn an_unset_chiral_flag_keeps_wedged_stereo_and_reports_it() {
    let interpretation = interpret(&wedged_center(false));
    let report = &interpretation.reports()[0];
    let molecule = interpretation.molecules().next().unwrap();
    let elements = molecule
        .stereo_elements()
        .filter(|(_, element)| element.is_specified())
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    assert_eq!(elements.len(), 1);
    assert_eq!(
        report.issues(),
        [MolfileInterpretIssue::AbsoluteStereoWithoutChiralFlag { elements }]
    );

    let interpretation = interpret(&wedged_center(true));
    assert!(interpretation.reports()[0].issues().is_empty());
    let molecule = interpretation.molecules().next().unwrap();
    assert_eq!(molecule.stereo_elements().count(), 1);

    // Without stereo, the flag means nothing.
    let unmarked = wedged_center(false).replace("  1  2  1  1", "  1  2  1  0");
    assert!(interpret(&unmarked).reports()[0].issues().is_empty());
}

#[test]
fn wedges_at_atoms_without_stereo_are_dropped_and_reported() {
    let cases = [
        // Ethanol, wedged from the methylene carbon to oxygen.
        (
            vec![("C", 0.0, 0.0), ("C", 1.3, 0.75), ("O", 2.6, 0.0)],
            vec![(1, 2, 1, 0), (2, 3, 1, 1)],
        ),
        // Ethanol, hashed from the methyl carbon.
        (
            vec![("C", 0.0, 0.0), ("C", 1.3, 0.75), ("O", 2.6, 0.0)],
            vec![(1, 2, 1, 6), (2, 3, 1, 0)],
        ),
        // Propene, wedged from the sp2 carbon to the methyl.
        (
            vec![("C", 0.0, 0.0), ("C", 1.3, 0.75), ("C", 2.6, 0.0)],
            vec![(1, 2, 2, 0), (2, 3, 1, 1)],
        ),
        // Trimethylamine, whose nitrogen inverts.
        (
            vec![
                ("N", 0.0, 0.0),
                ("C", 0.0, 1.0),
                ("C", 0.87, -0.5),
                ("C", -0.87, -0.5),
            ],
            vec![(1, 2, 1, 1), (1, 3, 1, 0), (1, 4, 1, 0)],
        ),
    ];
    for (atoms, bonds) in cases {
        let text = v2000(&atoms, &bonds, true);
        let interpretation = interpret(&text);
        let molecule = interpretation.molecules().next().unwrap();
        assert_eq!(molecule.stereo_elements().count(), 0, "{text}");
        let marked = bonds
            .iter()
            .position(|bond| bond.3 != 0)
            .expect("one marked bond");
        let report = &interpretation.reports()[0];
        let mapping = report.bond_mappings()[marked];
        assert_eq!(
            report.issues(),
            [MolfileInterpretIssue::StrayWedge {
                bond: mapping.bond(),
                location: mapping.location(),
            }],
            "{text}"
        );
    }
}

#[test]
fn an_unset_chiral_flag_on_3d_coordinates_is_not_reported() {
    // Bromochlorofluoromethane as a conformation, as RDKit and docking
    // programs write it: the flag unset and no marks.
    let text = "conformation\n  kekule            3D\n\n  5  4  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n    0.6000    0.6000    0.6000 F   0  0  0  0  0  0\n    0.6000   -0.6000   -0.6000 Cl  0  0  0  0  0  0\n   -0.6000    0.6000   -0.6000 Br  0  0  0  0  0  0\n   -0.6000   -0.6000    0.6000 H   0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n";
    let interpretation = interpret(text);
    let molecule = interpretation.molecules().next().unwrap();
    assert_eq!(
        molecule
            .stereo_elements()
            .filter(|(_, element)| element.is_specified())
            .count(),
        1
    );
    assert!(interpretation.reports()[0].issues().is_empty());
}

/// The atoms of a regular hexagon of carbons, 1.4 Å from `center`, the
/// first at `start` degrees and the rest counterclockwise.
fn hexagon(center: (f64, f64), start: f64) -> Vec<(&'static str, f64, f64)> {
    (0..6)
        .map(|k| {
            let angle = (start + 60.0 * f64::from(k)).to_radians();
            (
                "C",
                center.0 + 1.4 * angle.cos(),
                center.1 + 1.4 * angle.sin(),
            )
        })
        .collect()
}

/// A point `distance` beyond `atom` on the line from `center`.
fn outward(
    symbol: &'static str,
    center: (f64, f64),
    atom: (&str, f64, f64),
    distance: f64,
) -> (&'static str, f64, f64) {
    let (dx, dy) = (atom.1 - center.0, atom.2 - center.1);
    let scale = 1.0 + distance / dx.hypot(dy);
    (symbol, center.0 + dx * scale, center.1 + dy * scale)
}

/// The six ring bonds of the hexagon whose first atom is `first`, the bond
/// to the last atom drawn from the first, all of V2000 type `order`.
fn hexagon_bonds(first: usize, order: u8) -> Vec<(usize, usize, u8, u8)> {
    let mut bonds = (0..5)
        .map(|k| (first + k, first + k + 1, order, 0))
        .collect::<Vec<_>>();
    bonds.push((first, first + 5, order, 0));
    bonds
}

/// 2-Chloro-6-fluoro-2'-bromo-6'-iodobiphenyl with aromatic (type 4)
/// rings. Atom 1 ends the axis, atoms 2 and 6 are its ring neighbors.
fn biphenyl(marked: usize, stereo: u8) -> String {
    let (left, right) = ((-1.4, 0.0), (2.9, 0.0));
    let mut atoms = hexagon(left, 0.0);
    atoms.push(outward("Cl", left, atoms[1], 1.0));
    atoms.push(outward("F", left, atoms[5], 1.0));
    let ring = hexagon(right, 180.0);
    atoms.extend(ring.iter().copied());
    atoms.push(outward("Br", right, ring[1], 1.0));
    atoms.push(outward("I", right, ring[5], 1.0));
    let mut bonds = hexagon_bonds(1, 4);
    bonds.extend([(2, 7, 1, 0), (6, 8, 1, 0)]);
    bonds.extend(hexagon_bonds(9, 4));
    bonds.extend([(10, 15, 1, 0), (14, 16, 1, 0), (1, 9, 1, 0)]);
    for bond in &mut bonds {
        if (bond.0, bond.1) == (1, marked) {
            bond.3 = stereo;
        }
    }
    v2000(&atoms, &bonds, true)
}

#[test]
fn axis_wedges_on_aromatic_bonds_do_not_depend_on_the_kekule_structure() {
    // V2000 stereo code 1 is a wedge and 6 a hash.
    let axis = |marked: usize, stereo: u8| {
        let text = biphenyl(marked, stereo);
        let interpretation = interpret(&text);
        assert!(interpretation.reports()[0].issues().is_empty(), "{text}");
        let molecule = interpretation.molecules().next().unwrap();
        let elements = molecule
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>();
        let [element] = elements.as_slice() else {
            panic!("one axis: {elements:?}\n{text}");
        };
        assert!(
            matches!(element.kind, StereoElementKind::Axis(_)) && element.is_specified(),
            "{element:?}"
        );
        element.clone()
    };
    // A wedge to one ring neighbor states what a hash to the other does,
    // whichever of the two ring bonds localization makes double.
    assert_eq!(axis(2, 1), axis(6, 6));
    assert_eq!(axis(2, 6), axis(6, 1));
    assert_ne!(axis(2, 1), axis(2, 6));
}

#[test]
fn wedges_marking_stereo_kekule_cannot_represent_are_reported() {
    let ring = hexagon((1.4, 0.0), 180.0);
    let mut borane = vec![("B", 0.0, 0.0)];
    let mut borane_bonds = Vec::new();
    for (index, angle) in [0.0_f64, 120.0, 240.0].into_iter().enumerate() {
        let (sin, cos) = angle.to_radians().sin_cos();
        let first = borane.len() + 1;
        borane.extend(
            ring.iter()
                .map(|(symbol, x, y)| (*symbol, x * cos - y * sin, x * sin + y * cos)),
        );
        borane_bonds.extend(hexagon_bonds(first, 4));
        // V2000 stereo code 4 is a wavy bond.
        borane_bonds.push((1, first, 1, if index == 0 { 4 } else { 0 }));
    }
    let cases = [
        // A P-stereogenic phosphine bound to rhodium by a coordination
        // (type 9) bond.
        (
            vec![
                ("P", 0.0, 0.0),
                ("C", 0.0, 1.5),
                ("C", 1.3, -0.75),
                ("C", 2.6, 0.0),
                ("F", -1.3, -0.75),
                ("Rh", 0.0, -2.2),
            ],
            vec![
                (1, 2, 1, 1),
                (1, 3, 1, 0),
                (3, 4, 1, 0),
                (1, 5, 1, 0),
                (1, 6, 9, 0),
            ],
            UnrepresentableStereoReason::ZeroOrDativeBond,
        ),
        // A trigonal-bipyramidal phosphorus.
        (
            vec![
                ("P", 0.0, 0.0),
                ("F", 0.0, 1.6),
                ("Cl", 1.5, 0.3),
                ("Br", -1.5, 0.3),
                ("I", 0.9, -1.3),
                ("C", -0.9, -1.3),
            ],
            vec![
                (1, 2, 1, 1),
                (1, 3, 1, 0),
                (1, 4, 1, 0),
                (1, 5, 1, 0),
                (1, 6, 1, 0),
            ],
            UnrepresentableStereoReason::MoreThanFourNeighbors,
        ),
        // Triphenylborane, wavy from boron to one ring: both other
        // boron-carbon bonds could be the axis.
        (
            borane,
            borane_bonds,
            UnrepresentableStereoReason::SeveralCandidateAxes,
        ),
    ];
    for (atoms, bonds, reason) in cases {
        let text = v2000(&atoms, &bonds, true);
        let interpretation = interpret(&text);
        let molecule = interpretation.molecules().next().unwrap();
        assert_eq!(molecule.stereo_elements().count(), 0, "{text}");
        let marked = bonds
            .iter()
            .position(|bond| bond.3 != 0)
            .expect("one marked bond");
        let report = &interpretation.reports()[0];
        let mapping = report.bond_mappings()[marked];
        assert_eq!(
            report.issues(),
            [MolfileInterpretIssue::UnrepresentableStereo {
                atom: report.atom_mappings()[0].atom(),
                bond: mapping.bond(),
                reason,
                location: mapping.location(),
            }],
            "{text}"
        );
    }
}
