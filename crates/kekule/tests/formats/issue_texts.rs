use kekule::diagnostics::SourceLocation;
use kekule::molecule::{AtomId, BondId, Element};
use kekule::molfile::{MolfileInterpretIssue, UnrepresentableStereoReason};
use kekule::perception::valence::{ValenceError, ValenceViolation};
use kekule::stereo::StereoElementId;
use kekule::topology::TopologyAtomIndex;
use kekule::units::{Quantity, ANGSTROM};
use kekule::xyz::XyzInterpretIssue;

const REASONS: [UnrepresentableStereoReason; 3] = [
    UnrepresentableStereoReason::ZeroOrDativeBond,
    UnrepresentableStereoReason::MoreThanFourNeighbors,
    UnrepresentableStereoReason::SeveralCandidateAxes,
];

fn line(line: usize) -> SourceLocation {
    SourceLocation::at_line(line)
}

fn molfile_issues() -> Vec<MolfileInterpretIssue> {
    let mut issues = vec![
        MolfileInterpretIssue::ConflictingAtropisomericWedgeMarks {
            axis: BondId::new(2),
            location: line(14),
            mark_count: 2,
        },
        MolfileInterpretIssue::CoordinateStereoValenceUnsupported {
            error: ValenceError {
                violations: vec![ValenceViolation::ValenceExceeded {
                    atom: AtomId::new(3),
                    explicit_valence: 5,
                    max_allowed: 4,
                }],
            },
        },
        MolfileInterpretIssue::AtomParityContradictsCoordinates {
            center: AtomId::new(1),
            location: line(6),
        },
        MolfileInterpretIssue::ChargeFromBonds {
            atom: AtomId::new(0),
            location: line(5),
            charge: 1,
        },
        MolfileInterpretIssue::ChargeFromBonds {
            atom: AtomId::new(0),
            location: line(5),
            charge: -1,
        },
        MolfileInterpretIssue::ChargesSeparated {
            atom: AtomId::new(0),
            acceptor: AtomId::new(1),
            location: line(5),
        },
        MolfileInterpretIssue::MetalDonor {
            atom: AtomId::new(4),
            location: line(9),
        },
        MolfileInterpretIssue::MetalValenceIncomplete {
            atom: AtomId::new(0),
            location: line(5),
            valence: 3,
        },
        MolfileInterpretIssue::StrayWedge {
            bond: BondId::new(0),
            location: line(12),
        },
    ];
    for mark_count in [1, 3] {
        issues.push(MolfileInterpretIssue::AmbiguousTetrahedralWedgeMarks {
            center: AtomId::new(2),
            location: line(7),
            mark_count,
        });
    }
    for reason in REASONS {
        issues.push(MolfileInterpretIssue::UnrepresentableStereo {
            atom: AtomId::new(0),
            bond: BondId::new(1),
            reason,
            location: line(14),
        });
    }
    for count in [1, 2] {
        issues.push(MolfileInterpretIssue::AbsoluteStereoWithoutChiralFlag {
            elements: (0..count).map(StereoElementId::new).collect(),
        });
    }
    issues
}

fn xyz_issues() -> Vec<XyzInterpretIssue> {
    let atoms = [TopologyAtomIndex::new(1), TopologyAtomIndex::new(4)];
    vec![
        XyzInterpretIssue::BondOrdersAndChargesUnknown,
        XyzInterpretIssue::FallbackCovalentRadius {
            atom: TopologyAtomIndex::new(0),
            element: Element::from_symbol("Og").unwrap(),
        },
        XyzInterpretIssue::CloseContact {
            atoms,
            distance: Quantity::new(0.314, ANGSTROM),
        },
        XyzInterpretIssue::LongContactNotBonded {
            atoms,
            distance: Quantity::new(1.8512, ANGSTROM),
        },
        XyzInterpretIssue::ValenceExceeded {
            atom: TopologyAtomIndex::new(2),
            element: Element::C,
            bonds: 6,
            maximum: 5,
        },
        XyzInterpretIssue::ExtraColumnsOmitted {
            atom: TopologyAtomIndex::new(2),
            location: line(5),
        },
    ]
}

fn assert_sentences(text: &str) {
    assert!(!text.contains("  "), "double space in {text:?}");
    assert!(text.ends_with('.'), "no final period in {text:?}");
    assert!(
        text.starts_with(|first: char| first.is_uppercase()),
        "no capital letter in {text:?}"
    );
}

#[test]
fn issue_texts_are_sentences() {
    let texts = molfile_issues()
        .into_iter()
        .map(|issue| issue.to_string())
        .chain(xyz_issues().into_iter().map(|issue| issue.to_string()))
        .chain(REASONS.iter().map(ToString::to_string));
    for text in texts {
        assert_sentences(&text);
    }
}

#[test]
fn molfile_issues_name_atoms_and_bonds_by_line() {
    let texts = molfile_issues()
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    for expected in [
        "The axis bond at line 14 has 2 conflicting wedge marks. Its configuration is not read.",
        "Atom 3 has valence 5, more than the allowed 4. No stereo is read from the coordinates.",
        "The atom at line 5 takes charge +1 from its bonds.",
        "The atom at line 5 takes charge -1 from its bonds.",
        "The wedge mark at the atom at line 7 does not set one configuration. Its \
         configuration is not read.",
        "The 3 wedge marks at the atom at line 7 do not set one configuration. Its \
         configuration is not read.",
        "The stereo mark on the bond at line 14 is not read. The marked atom has more than \
         four neighbors.",
        "The chiral flag is unset, but 1 drawn configuration is read as absolute.",
        "The chiral flag is unset, but 2 drawn configurations are read as absolute.",
    ] {
        assert!(
            texts.iter().any(|text| text == expected),
            "{expected:?} in {texts:#?}"
        );
    }
}

#[test]
fn xyz_issues_number_atoms_from_one() {
    let texts = xyz_issues()
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    for expected in [
        "Atom 1 (Og) has no tabulated covalent radius. Its bonds use the fallback radius.",
        "Atoms 2 and 5 are only 0.31 Å apart. They are not bonded.",
        "Atoms 2 and 5 are 1.85 Å apart, near bonding distance. They are not bonded.",
        "Atom 3 (C) has 6 single bonds, more than its largest valence of 5. Its distances do \
         not give its connectivity.",
        "The columns after the coordinates of atom 3 at line 5 are not read.",
    ] {
        assert!(
            texts.iter().any(|text| text == expected),
            "{expected:?} in {texts:#?}"
        );
    }
}

#[test]
fn interpreted_xyz_reports_carry_displayable_issues() {
    let interpretation = kekule::xyz::parse_str(
        "3
water
O 0 0 0
H 0.96 0 0
H -0.24 0.93 0
",
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap();
    let issues = interpretation.report().issues();
    assert!(issues.contains(&XyzInterpretIssue::BondOrdersAndChargesUnknown));
    for issue in issues {
        assert_sentences(&issue.to_string());
    }
}
