//! Molfile stereo marks and the chiral flag at interpretation.

use kekule::molfile::{self, MolfileInterpretIssue, MolfileInterpretation};

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
