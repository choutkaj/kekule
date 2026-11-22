//! Molfile stereo marks and the chiral flag at interpretation.

use kekule::diagnostics::{DiagnosticKind, SourceLocation};
use kekule::molecule::{AtomId, ImplicitHydrogens, Molecule};
use kekule::molfile::{
    self, MolfileInterpretIssue, MolfileInterpretation, MolfileInterpretationReport,
    UnrepresentableStereoReason,
};
use kekule::stereo::{self, StereoCarrier, StereoDescriptor, StereoElementKind};

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

/// The published molecule and report of a one-component record.
fn read_with_report(text: &str) -> (Molecule, MolfileInterpretationReport) {
    let interpretation = interpret(text);
    assert_eq!(interpretation.reports().len(), 1, "{text}");
    let report = interpretation.reports()[0].clone();
    (interpretation.into_molecules().pop().unwrap(), report)
}

fn chain_text(error: &(dyn std::error::Error + 'static)) -> String {
    let mut text = error.to_string();
    let mut source = error.source();
    while let Some(cause) = source {
        text.push_str(": ");
        text.push_str(&cause.to_string());
        source = cause.source();
    }
    text
}

/// A V2000 tetrahedron of `symbol` with V2000 charge code `charge_code`,
/// wedged from the center to fluorine, with the chiral flag set.
fn wedged_tetrahedron(symbol: &str, charge_code: u8, mirror: bool) -> String {
    let y = if mirror { -1.0 } else { 1.0 };
    format!(
        "stereo drawing\nkekule\n\n  4  3  0  0  1  0            999 V2000\n    0.0000    0.0000    0.0000 {symbol:<3} 0  {charge_code}  0  0  0  0\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0\n    0.0000{y:10.4}    0.0000 Br  0  0  0  0  0  0\n  1  2  1  1  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\nM  END\n"
    )
}

/// A carbon at `points[0]` bonded to F, Cl, Br, and I at the other points,
/// wedged to fluorine, with the chiral flag set.
fn wedge_drawing(points: &[[f64; 2]], scale: f64, mirror: f64, v3000: bool) -> String {
    let symbols = ["C", "F", "Cl", "Br", "I"];
    let mut source = format!(
        "drawing\nkekule\n\n{:3}{:3}  0  0  1  0            999 V2000\n",
        points.len(),
        points.len() - 1
    );
    if v3000 {
        source = format!("drawing\nkekule\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {} {} 0 0 1\nM  V30 BEGIN ATOM\n", points.len(), points.len() - 1);
    }
    for (index, [x, y]) in points.iter().enumerate() {
        let (x, y) = (scale * (x + 3.0), scale * (mirror * y - 2.0));
        if v3000 {
            source += &format!("M  V30 {} {} {x} {y} 0 0\n", index + 1, symbols[index]);
        } else {
            source += &format!(
                "{x:10.4}{y:10.4}    0.0000 {:<3} 0  0  0  0  0  0\n",
                symbols[index]
            );
        }
    }
    if v3000 {
        source += "M  V30 END ATOM\nM  V30 BEGIN BOND\n";
    }
    for index in 1..points.len() {
        if v3000 {
            source += &format!(
                "M  V30 {index} 1 1 {}{}\n",
                index + 1,
                if index == 1 { " CFG=1" } else { "" }
            );
        } else {
            source += &format!(
                "  1{:3}  1{:3}  0  0  0\n",
                index + 1,
                usize::from(index == 1)
            );
        }
    }
    if v3000 {
        source += "M  V30 END BOND\nM  V30 END CTAB\n";
    }
    source + "M  END\n"
}

#[test]
fn bond_stereo_needs_enough_source_context() {
    // A wedge at a methyl is dropped with a report; a crossed double bond
    // without carriers fails.
    let input = |order_code: u8, stereo_code: u8| {
        format!(
            "stereo\nkekule\n\n  2  1  0  0  0  0            999 V2000\n   -1.2500    0.0000    0.0000 C   0  0  0  0  0  0\n    1.2500    0.0000    0.0000 C   0  0  0  0  0  0\n  1  2  {order_code}  {stereo_code}  0  0  0\nM  END\n"
        )
    };
    for stereo_code in [1, 4, 6] {
        let interpretation = interpret(&input(1, stereo_code));
        assert!(
            matches!(
                interpretation.reports()[0].issues(),
                [MolfileInterpretIssue::StrayWedge { location, .. }]
                    if *location == SourceLocation::at_line(7)
            ),
            "stereo code {stereo_code}"
        );
    }
    let error = molfile::parse_str(&input(2, 3), Default::default())
        .expect("bond stereo syntax should parse")
        .interpret(Default::default())
        .expect_err("under-specified stereo must not publish a molecule");
    assert_eq!(error.kind(), DiagnosticKind::Chemistry);
    assert_eq!(error.location(), Some(SourceLocation::at_line(7)));
    let text = chain_text(&error);
    assert!(text.contains("on bond 0"), "{text}");
}

#[test]
fn omitted_tetrahedral_hydrogens_come_from_source_valence() {
    for symbol in ["C", "N", "S", "P"] {
        let input = format!(
            "stereo hydrogen\nkekule\n\n  4  3  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 {symbol:<3} 0  0  0  0  0  0\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0\n    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0\n  1  2  1  1  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\nM  END\n"
        );
        let molecule = interpret(&input).into_molecules().pop().unwrap();
        if symbol == "N" {
            // An amine nitrogen inverts, so its wedge marks no stereo center.
            assert_eq!(molecule.stereo_elements().count(), 0);
            continue;
        }
        let expected_carrier = if symbol == "P" {
            StereoCarrier::ImplicitLonePair
        } else {
            assert_eq!(
                molecule.atom(AtomId::new(0)).unwrap().hydrogens(),
                ImplicitHydrogens::Fixed(1)
            );
            StereoCarrier::ImplicitHydrogen
        };
        assert!(molecule.perception().valence().is_none());
        assert_eq!(molecule.stereo_elements().count(), 1);
        assert!(molecule.stereo_elements().any(|(_, element)| {
            matches!(
                &element.kind,
                StereoElementKind::Tetrahedral(stereo)
                    if stereo.carriers.contains(&expected_carrier)
            )
        }));
    }
}

#[test]
fn valid_wedge_geometry_is_scale_and_reflection_consistent() {
    // RDKit 2026.03.3 independently checks these drawing configurations.
    for (points, original) in [
        (
            vec![[0.0, 0.0], [1.0, 0.0], [-1.0, 0.0], [0.0, 1.0]],
            StereoDescriptor::S,
        ),
        (
            vec![[0.0, 0.0], [1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]],
            StereoDescriptor::R,
        ),
        (
            vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0]],
            StereoDescriptor::R,
        ),
        (
            vec![[0.0, 0.0], [0.0, 0.0], [-1.0, 0.0], [0.0, 1.0]],
            StereoDescriptor::S,
        ),
        // Local drawings from PubChem 10524, with distinct test ligands.
        // The three unmarked bonds lie on the same side of the center;
        // their drawing order must not be read as a displaced tetrahedron.
        (
            vec![
                [6.3981, -0.4499],
                [6.3820, -1.4914],
                [7.2641, 0.0501],
                [6.1950, -0.0800],
                [5.5320, 0.0501],
            ],
            StereoDescriptor::S,
        ),
        (
            vec![
                [7.2641, 1.0501],
                [7.6891, 1.7862],
                [6.8422, 1.0410],
                [7.2641, 0.0501],
                [6.3981, 1.5501],
            ],
            StereoDescriptor::S,
        ),
    ] {
        for v3000 in [false, true] {
            for scale in [0.1, 1.0, 10.0] {
                for mirror in [-1.0, 1.0] {
                    let source = wedge_drawing(&points, scale, mirror, v3000);
                    let (mut molecule, report) = read_with_report(&source);
                    assert!(report.issues().is_empty(), "{source}");
                    molecule.perceive().unwrap();
                    let assignment =
                        stereo::cip::assign_cip_descriptors(&mut molecule, Default::default())
                            .unwrap();
                    let expected = if mirror > 0.0 {
                        original
                    } else if original == StereoDescriptor::R {
                        StereoDescriptor::S
                    } else {
                        StereoDescriptor::R
                    };
                    assert_eq!(assignment.assigned().len(), 1, "{source}");
                    assert_eq!(assignment.assigned()[0].descriptor(), expected, "{source}");
                }
            }
        }
    }
}

#[test]
fn redundant_wedges_preserve_consistent_and_unknown_configurations() {
    let source = wedged_tetrahedron("C", 0, false)
        .replace("   -1.0000    0.0000", "   -0.5000    0.8660")
        .replace("    0.0000    1.0000", "   -0.5000   -0.8660");
    let original = interpret(&source).into_molecules().pop().unwrap();
    let redundant = source.replace("  1  3  1  0", "  1  3  1  1");
    let (molecule, report) = read_with_report(&redundant);
    assert!(report.issues().is_empty());
    assert_eq!(
        molecule.stereo_elements().next().unwrap().1.kind,
        original.stereo_elements().next().unwrap().1.kind
    );

    let unknown = source.replace("  1  3  1  0", "  1  3  1  4");
    let molecule = interpret(&unknown).into_molecules().pop().unwrap();
    assert!(molecule
        .stereo_elements()
        .next()
        .unwrap()
        .1
        .is_explicitly_unknown());

    let degenerate = unknown
        .replace("    1.0000    0.0000", "    0.0000    0.0000")
        .replace("   -0.5000    0.8660", "    0.0000    0.0000")
        .replace("   -0.5000   -0.8660", "    0.0000    0.0000");
    let (molecule, report) = read_with_report(&degenerate);
    assert!(report.issues().is_empty());
    assert!(molecule
        .stereo_elements()
        .next()
        .unwrap()
        .1
        .is_explicitly_unknown());

    let conflicting = source.replace("  1  3  1  0", "  1  3  1  6");
    let (molecule, report) = read_with_report(&conflicting);
    assert_eq!(report.issues().len(), 1);
    assert!(molecule.stereo_elements().next().is_none());
}
