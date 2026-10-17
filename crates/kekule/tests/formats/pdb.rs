use kekule::diagnostics::{DiagnosticKind, SourceLocation};
use kekule::macromolecular::{
    AltLocPolicy, BondEvidence, CellOmissionReason, ConnectionResolutionReason, EntityKind,
    ModelSelection, SourceResidueId, SourceResiduePosition, StructureInterpretIssue,
    StructureInterpretOptions, StructureInterpretation,
};
use kekule::molecule::BondOrder;
use kekule::pdb::{self, PdbAtomKind, PdbParseOptions};

/// One `ATOM` or `HETATM` record in the 80 fixed columns of the PDB format.
#[derive(Clone, Copy)]
struct Site {
    record: &'static str,
    serial: &'static str,
    /// Columns 13-16, aligned as written.
    name: &'static str,
    alternate_location: char,
    residue: &'static str,
    chain: char,
    sequence: &'static str,
    insertion_code: char,
    position: [f64; 3],
    occupancy: f64,
    temperature_factor: f64,
    segment: &'static str,
    element: &'static str,
    charge: &'static str,
}

impl Site {
    fn line(&self) -> String {
        let [x, y, z] = self.position;
        // Residue names fill columns 18-20, or 18-21 when four long.
        let residue = if self.residue.len() == 4 {
            self.residue.to_owned()
        } else {
            format!("{:>3} ", self.residue)
        };
        format!(
            "{:<6}{:>5} {:<4}{}{residue}{}{:>4}{}   {x:>8.3}{y:>8.3}{z:>8.3}{:>6.2}{:>6.2}      {:<4}{:>2}{:<2}",
            self.record,
            self.serial,
            self.name,
            self.alternate_location,
            self.chain,
            self.sequence,
            self.insertion_code,
            self.occupancy,
            self.temperature_factor,
            self.segment,
            self.element,
            self.charge,
        )
    }

    /// The record without the columns after the coordinates.
    fn coordinates_only(&self) -> String {
        self.line()[..54].to_owned()
    }
}

#[allow(clippy::too_many_arguments)]
fn site(
    record: &'static str,
    serial: &'static str,
    name: &'static str,
    residue: &'static str,
    chain: char,
    sequence: &'static str,
    position: [f64; 3],
    element: &'static str,
) -> Site {
    Site {
        record,
        serial,
        name,
        alternate_location: ' ',
        residue,
        chain,
        sequence,
        insertion_code: ' ',
        position,
        occupancy: 1.0,
        temperature_factor: 20.0,
        segment: "",
        element,
        charge: "",
    }
}

/// A `TER` record ending the chain of the residue it names.
fn ter(serial: &str, residue: &str, chain: char, sequence: &str) -> String {
    format!(
        "{:<6}{serial:>5}{:6}{residue:>3} {chain}{sequence:>4}",
        "TER", ""
    )
}

/// One partner of a `LINK` or `SSBOND` record: residue, chain, sequence.
type Partner<'a> = (&'a str, char, &'a str);

/// A `LINK` record between two named atoms; names are aligned as in
/// columns 13-16 of a coordinate record.
fn link(
    atoms: [&str; 2],
    partners: [Partner<'_>; 2],
    symmetry: [&str; 2],
    distance: f64,
) -> String {
    let [(left, left_chain, left_sequence), (right, right_chain, right_sequence)] = partners;
    format!(
        "{:<12}{:<4} {left:>3} {left_chain}{left_sequence:>4}{:16}{:<4} {right:>3} {right_chain}{right_sequence:>4}{:3}{:>6} {:>6} {distance:>5.2}",
        "LINK", atoms[0], "", atoms[1], "", symmetry[0], symmetry[1],
    )
}

/// An `SSBOND` record between two cysteines.
fn ssbond(serial: &str, partners: [Partner<'_>; 2], distance: f64) -> String {
    let [(left, left_chain, left_sequence), (right, right_chain, right_sequence)] = partners;
    format!(
        "{:<6} {serial:>3} {left:>3} {left_chain} {left_sequence:>4}{:4}{right:>3} {right_chain} {right_sequence:>4}{:24}{:>6} {:>6} {distance:>5.2}",
        "SSBOND", "", "", "1555", "1555",
    )
}

fn entry(lines: &[String]) -> String {
    lines.join("\n") + "\nEND\n"
}

fn interpret(text: &str) -> StructureInterpretation {
    (pdb::parse_str(text, Default::default()).unwrap())
        .interpret(Default::default())
        .unwrap()
}

/// Glycine dipeptide of chain A, in records 1-8.
fn glycine_dipeptide() -> Vec<String> {
    [
        ("1", " N  ", "1", [0.000, 0.000, 0.0], "N"),
        ("2", " CA ", "1", [1.458, 0.000, 0.0], "C"),
        ("3", " C  ", "1", [2.009, 1.420, 0.0], "C"),
        ("4", " O  ", "1", [1.251, 2.390, 0.0], "O"),
        ("5", " N  ", "2", [3.332, 1.536, 0.0], "N"),
        ("6", " CA ", "2", [3.988, 2.836, 0.0], "C"),
        ("7", " C  ", "2", [5.488, 2.736, 0.0], "C"),
        ("8", " O  ", "2", [6.100, 1.680, 0.0], "O"),
    ]
    .into_iter()
    .map(|(serial, name, sequence, position, element)| {
        site(
            "ATOM", serial, name, "GLY", 'A', sequence, position, element,
        )
        .line()
    })
    .collect()
}

/// A cysteine whose `SG` lies at x = 3.8 + `offset`, mirrored in x when
/// `mirror` is -1.
fn cysteine(chain: char, serials: [&'static str; 6], offset: f64, mirror: f64) -> Vec<String> {
    [
        (" N  ", [0.000, 0.000, 0.00], "N"),
        (" CA ", [1.460, 0.000, 0.00], "C"),
        (" C  ", [2.000, 1.420, 0.00], "C"),
        (" O  ", [1.250, 2.390, 0.00], "O"),
        (" CB ", [2.000, -0.750, 1.15], "C"),
        (" SG ", [3.800, -0.800, 1.20], "S"),
    ]
    .into_iter()
    .zip(serials)
    .map(|((name, [x, y, z], element), serial)| {
        let x = 3.8 + offset + mirror * (x - 3.8);
        site("ATOM", serial, name, "CYS", chain, "1", [x, y, z], element).line()
    })
    .collect()
}

/// A water of chain W, in records 1-3.
fn water(offset: f64) -> Vec<String> {
    [
        ("1", " O  ", [0.0, 0.0, 0.0], "O"),
        ("2", " H1 ", [0.96, 0.0, 0.0], "H"),
        ("3", " H2 ", [-0.24, 0.93, 0.0], "H"),
    ]
    .into_iter()
    .map(|(serial, name, [x, y, z], element)| {
        site(
            "HETATM",
            serial,
            name,
            "HOH",
            'W',
            "1",
            [x + offset, y, z],
            element,
        )
        .line()
    })
    .collect()
}

/// Bonds by atom names, order, and evidence, each pair in name order.
fn bonds(result: &StructureInterpretation) -> Vec<(String, String, BondOrder, BondEvidence)> {
    let name = |atom| {
        result
            .report()
            .instances()
            .iter()
            .flat_map(|instance| instance.atoms())
            .find(|candidate| candidate.atom() == atom)
            .unwrap()
            .atom_name()
            .to_owned()
    };
    let mut bonds = result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.bonds())
        .map(|bond| {
            let [a, b] = bond.atoms();
            let molecule = result.topology().instance(a.instance()).unwrap().molecule();
            let id = molecule.bond_between(a.atom(), b.atom()).unwrap();
            let (a, b) = (name(a), name(b));
            let (a, b) = if a <= b { (a, b) } else { (b, a) };
            (
                a,
                b,
                molecule.bond(id).unwrap().order,
                bond.evidence().clone(),
            )
        })
        .collect::<Vec<_>>();
    bonds.sort_by(|left, right| (&left.0, &left.1).cmp(&(&right.0, &right.1)));
    bonds
}

/// Charged atoms by name, in instance order.
fn charges(result: &StructureInterpretation) -> Vec<(String, i8)> {
    result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .filter_map(|atom| {
            let molecule = result
                .topology()
                .instance(atom.atom().instance())
                .unwrap()
                .molecule();
            let charge = molecule.atom(atom.atom().atom()).unwrap().formal_charge();
            (charge != 0).then(|| (atom.atom_name().to_owned(), charge))
        })
        .collect()
}

fn connection(id: Option<&str>, kind: &str) -> BondEvidence {
    BondEvidence::Connection {
        connection_id: id.map(str::to_owned),
        connection_type: kind.to_owned(),
    }
}

fn residue(asym: &str, position: SourceResiduePosition) -> SourceResidueId {
    SourceResidueId::from_parts("1", asym, position, None)
}

#[test]
fn pdb_records_parse_by_fixed_columns() {
    let nitrogen = Site {
        alternate_location: 'A',
        insertion_code: 'A',
        occupancy: 0.6,
        temperature_factor: 15.25,
        segment: "SEGA",
        charge: "1+",
        ..site(
            "ATOM",
            "1",
            " N  ",
            "GLY",
            'A',
            "-1",
            [-1.25, 2.5, 30.125],
            "N",
        )
    };
    let iron = site(
        "HETATM",
        "A0000",
        "FE  ",
        "HEM",
        'B',
        "9999",
        [10.0, -0.5, 1.0],
        "FE",
    );
    let water = site(
        "HETATM",
        "3",
        " O  ",
        "HOH",
        'W',
        "A000",
        [1.0, 2.0, 3.0],
        "O",
    );
    let text = [
        "HEADER    HYDROLASE                               01-JAN-00   1ABC              ",
        "REMARK   2 RESOLUTION.    1.50 ANGSTROMS.",
        &nitrogen.line(),
        &iron.line(),
        &water.coordinates_only(),
        &ter("4", "HEM", 'B', "9999"),
        &link(
            ["FE", " NE2"],
            [("HEM", 'B', "9999"), ("HIS", 'A', "93")],
            ["1555", "2555"],
            2.09,
        ),
        &ssbond("1", [("CYS", 'A', "6"), ("CYS", 'A', "127")], 2.03),
        "CONECT    1A0000    3    3",
        "CONECT    3    1",
        "MASTER        0    0    0    0    0    0    0    0    3    1    2    0",
        // The unit cell of 1AKE.
        "CRYST1   73.200   79.800   85.000  90.00  90.00  90.00 P 21 2 21     8          ",
    ]
    .join("\r\n")
        + "\rEND\r\n";
    let document = pdb::parse_str(&text, Default::default()).unwrap();
    assert_eq!(document.id_code(), Some("1ABC"));
    assert!(document.models().is_empty());

    let atoms = document.atoms();
    assert_eq!(atoms.len(), 3);
    let first = &atoms[0];
    assert_eq!(first.kind(), PdbAtomKind::Atom);
    assert_eq!(first.serial(), Some(1));
    assert_eq!(first.name(), "N");
    assert_eq!(first.name_field(), " N  ");
    assert_eq!(first.alternate_location(), Some('A'));
    assert_eq!(first.residue_name(), "GLY");
    assert_eq!(first.chain_id(), Some('A'));
    assert_eq!(first.residue_sequence(), -1);
    assert_eq!(first.insertion_code(), Some('A'));
    assert_eq!(
        [first.position().x, first.position().y, first.position().z],
        [-1.25, 2.5, 30.125]
    );
    assert_eq!(first.occupancy(), Some(0.6));
    assert_eq!(first.temperature_factor(), Some(15.25));
    assert_eq!(first.segment_id(), Some("SEGA"));
    assert_eq!(first.element(), Some("N"));
    assert_eq!(first.formal_charge(), Some(1));
    assert_eq!(first.model(), None);
    assert_eq!(first.location(), SourceLocation::at_line(3));

    // Hybrid-36 numbers continue past the decimal columns.
    let iron = &atoms[1];
    assert_eq!(iron.kind(), PdbAtomKind::Hetatm);
    assert_eq!(iron.serial(), Some(100_000));
    assert_eq!(iron.name_field(), "FE  ");
    assert_eq!(iron.residue_sequence(), 9999);
    assert_eq!(iron.element(), Some("FE"));
    assert_eq!(iron.formal_charge(), None);
    let water = &atoms[2];
    assert_eq!(water.residue_sequence(), 10_000);
    assert_eq!(water.chain_id(), Some('W'));
    assert_eq!(
        (
            water.occupancy(),
            water.temperature_factor(),
            water.element()
        ),
        (None, None, None)
    );

    let [ter] = document.ter_records() else {
        panic!("one TER record");
    };
    assert_eq!(
        (
            ter.serial(),
            ter.residue_name(),
            ter.chain_id(),
            ter.residue_sequence()
        ),
        (Some(4), Some("HEM"), Some('B'), Some(9999))
    );
    assert_eq!(ter.preceding_atoms(), 3);

    let [link] = document.link_records() else {
        panic!("one LINK record");
    };
    let [iron, histidine] = link.partners();
    assert_eq!(
        (
            iron.atom_name(),
            iron.residue_name(),
            iron.chain_id(),
            iron.residue_sequence(),
            iron.symmetry()
        ),
        (Some("FE"), "HEM", Some('B'), 9999, Some("1555"))
    );
    assert_eq!(
        (
            histidine.atom_name(),
            histidine.residue_name(),
            histidine.chain_id(),
            histidine.residue_sequence(),
            histidine.symmetry()
        ),
        (Some("NE2"), "HIS", Some('A'), 93, Some("2555"))
    );
    assert_eq!(link.distance(), Some(2.09));

    let [disulfide] = document.ssbond_records() else {
        panic!("one SSBOND record");
    };
    assert_eq!(disulfide.serial(), Some(1));
    let [left, right] = disulfide.partners();
    assert_eq!(
        (
            left.atom_name(),
            left.residue_name(),
            left.chain_id(),
            left.residue_sequence()
        ),
        (None, "CYS", Some('A'), 6)
    );
    assert_eq!(right.residue_sequence(), 127);
    assert_eq!(disulfide.distance(), Some(2.03));

    let conect = document.conect_records();
    assert_eq!(conect.len(), 2);
    assert_eq!(conect[0].serial(), 1);
    assert_eq!(conect[0].bonded(), &[100_000, 3, 3]);
    assert_eq!(conect[1].location(), SourceLocation::at_line(10));

    let other = document
        .other_records()
        .iter()
        .map(|record| (record.name(), record.location().line()))
        .collect::<Vec<_>>();
    assert_eq!(other, [("HEADER", 1), ("REMARK", 2), ("MASTER", 11)]);

    let cryst1 = document.cryst1().unwrap();
    assert_eq!(cryst1.lengths(), [Some(73.2), Some(79.8), Some(85.0)]);
    assert_eq!(cryst1.angles(), [Some(90.0); 3]);
    assert_eq!(cryst1.space_group(), Some("P 21 2 21"));
    assert_eq!(cryst1.z(), Some(8));
    assert_eq!(cryst1.location(), SourceLocation::at_line(12));
}

/// A `CRYST1` record in its fixed columns; `None` leaves a field blank.
fn cryst1(lengths: [Option<f64>; 3], space_group: &str) -> String {
    let [a, b, c] =
        lengths.map(|length| length.map_or(String::new(), |length| format!("{length:.3}")));
    format!(
        "CRYST1{a:>9}{b:>9}{c:>9}{:>7}{:>7}{:>7} {space_group:<11}{:>4}",
        "90.00", "90.00", "90.00", "1"
    )
}

#[test]
fn cryst1_is_the_model_cell_only_for_a_whole_unit_cell() {
    let space_group = |name: &str| CellOmissionReason::SpaceGroup {
        name: name.to_owned(),
    };
    for (record, expected, omission) in [
        (Some(cryst1([Some(30.0); 3], "P 1")), Some([30.0; 3]), None),
        // An asymmetric unit's translates are not its crystal.
        (
            Some(cryst1([Some(73.2), Some(79.8), Some(85.0)], "P 21 2 21")),
            None,
            Some(space_group("P 21 2 21")),
        ),
        // The placeholder of structures without a crystal lattice.
        (
            Some(cryst1([Some(1.0); 3], "P 1")),
            None,
            Some(CellOmissionReason::Placeholder),
        ),
        (
            Some(cryst1([Some(30.0), Some(30.0), None], "P 1")),
            None,
            Some(CellOmissionReason::Incomplete),
        ),
        (None, None, None),
    ] {
        let text = entry(
            &record
                .iter()
                .cloned()
                .chain(glycine_dipeptide())
                .collect::<Vec<_>>(),
        );
        let interpreted = interpret(&text);
        let cell = interpreted.model().cell().map(|cell| {
            assert!(cell.is_orthorhombic());
            assert_eq!(cell.periodic_axes(), [true; 3]);
            let [a, b, c] = cell
                .vectors()
                .into_unit(kekule::units::ANGSTROM)
                .unwrap()
                .into_value();
            [a.x, b.y, c.z]
        });
        match (cell, expected) {
            (Some(cell), Some(expected)) => {
                for (cell, expected) in cell.into_iter().zip(expected) {
                    assert!((cell - expected).abs() < 1e-9, "{record:?}");
                }
            }
            (cell, expected) => assert_eq!(cell, expected, "{record:?}"),
        }
        let omissions = interpreted
            .report()
            .issues()
            .iter()
            .filter(|issue| matches!(issue, StructureInterpretIssue::CellOmitted { .. }))
            .cloned()
            .collect::<Vec<_>>();
        let expected = omission
            .map(|reason| StructureInterpretIssue::CellOmitted {
                reason,
                location: SourceLocation::at_line(1),
            })
            .into_iter()
            .collect::<Vec<_>>();
        assert_eq!(omissions, expected, "{record:?}");
    }

    let duplicate = entry(&[
        cryst1([Some(30.0); 3], "P 1"),
        cryst1([Some(30.0); 3], "P 1"),
    ]);
    let error = pdb::parse_str(&duplicate, Default::default()).unwrap_err();
    assert_eq!(
        (error.kind(), error.location(), error.message()),
        (
            DiagnosticKind::Inconsistent,
            Some(SourceLocation::at_line(2)),
            "duplicate CRYST1 record"
        )
    );
}

#[test]
fn ter_records_end_the_polymer_of_their_chain() {
    let ligands = [
        site(
            "HETATM",
            "10",
            " O  ",
            "HOH",
            'A',
            "101",
            [10.0, 10.0, 10.0],
            "O",
        )
        .line(),
        site(
            "HETATM",
            "11",
            "ZN  ",
            "ZN",
            'A',
            "102",
            [-5.0, -5.0, -5.0],
            "ZN",
        )
        .line(),
    ];
    let terminated = [
        glycine_dipeptide(),
        vec![ter("9", "GLY", 'A', "2")],
        ligands.to_vec(),
    ]
    .concat();
    // Without TER, the polymer spans the chain's ATOM records.
    let unterminated = [glycine_dipeptide(), ligands.to_vec()].concat();
    for lines in [terminated, unterminated] {
        let interpreted = interpret(&entry(&lines));
        let instances = interpreted
            .report()
            .instances()
            .iter()
            .map(|instance| (instance.atoms().len(), instance.entity_kinds().to_vec()))
            .collect::<Vec<_>>();
        assert_eq!(
            instances,
            [
                (8, vec![EntityKind::Polymer]),
                (1, vec![EntityKind::Water]),
                (1, vec![EntityKind::NonPolymer]),
            ]
        );
        let report = interpreted.report();
        assert_eq!(
            (
                report.macromolecules(),
                report.small_molecules(),
                report.solvent_molecules()
            ),
            (1, 2, 1)
        );
        assert_eq!(interpreted.topology().bond_count(), 7);
        let peptide = bonds(&interpreted)
            .into_iter()
            .filter(|bond| bond.3 == BondEvidence::PolymerLink)
            .collect::<Vec<_>>();
        assert_eq!(
            peptide,
            [(
                "C".to_owned(),
                "N".to_owned(),
                BondOrder::Single,
                BondEvidence::PolymerLink
            )]
        );
        assert_eq!(charges(&interpreted), [("ZN".to_owned(), 2)]);
        assert_eq!(report.issues(), &[]);
    }
}

#[test]
fn conect_records_list_bonds_that_templates_outrank() {
    let template = |component: &str| BondEvidence::Template {
        component_id: component.to_owned(),
        variant: None,
    };
    let bond = |a: &str, b: &str, order, evidence| (a.to_owned(), b.to_owned(), order, evidence);
    // An untemplated ligand takes the listed bonds; a partner listed twice is
    // a double bond, and a hydrogen the listing omits attaches by distance.
    let ligand = [
        site("HETATM", "1", " C1 ", "LIG", 'L', "1", [0.0, 0.0, 0.0], "C").line(),
        site("HETATM", "2", " C2 ", "LIG", 'L', "1", [1.5, 0.0, 0.0], "C").line(),
        site(
            "HETATM",
            "3",
            " O1 ",
            "LIG",
            'L',
            "1",
            [2.2, 1.05, 0.0],
            "O",
        )
        .line(),
        site(
            "HETATM",
            "4",
            " H1 ",
            "LIG",
            'L',
            "1",
            [-0.55, 0.95, 0.0],
            "H",
        )
        .line(),
        "CONECT    1    2".to_owned(),
        "CONECT    2    1    3    3".to_owned(),
        "CONECT    3    2    2".to_owned(),
    ];
    // A templated residue keeps its template's orders.
    let acetate = [
        site("HETATM", "1", " C  ", "ACT", 'L', "1", [0.0, 0.0, 0.0], "C").line(),
        site(
            "HETATM",
            "2",
            " O  ",
            "ACT",
            'L',
            "1",
            [0.61, 1.06, 0.0],
            "O",
        )
        .line(),
        site(
            "HETATM",
            "3",
            " OXT",
            "ACT",
            'L',
            "1",
            [0.61, -1.06, 0.0],
            "O",
        )
        .line(),
        site(
            "HETATM",
            "4",
            " CH3",
            "ACT",
            'L',
            "1",
            [-1.5, 0.0, 0.0],
            "C",
        )
        .line(),
        "CONECT    1    2    3    4".to_owned(),
        "CONECT    2    1".to_owned(),
        "CONECT    3    1".to_owned(),
        "CONECT    4    1".to_owned(),
    ];
    let cases = [
        (
            &ligand[..],
            vec![
                bond("C1", "C2", BondOrder::Single, connection(None, "CONECT")),
                bond(
                    "C1",
                    "H1",
                    BondOrder::Single,
                    BondEvidence::HydrogenDistance,
                ),
                bond("C2", "O1", BondOrder::Double, connection(None, "CONECT")),
            ],
            vec![],
            vec![StructureInterpretIssue::HydrogenAttachedByDistance {
                residue: residue("L", SourceResiduePosition::Author("1".to_owned())),
                atom_name: "H1".to_owned(),
                heavy_atom_name: "C1".to_owned(),
            }],
            2,
        ),
        (
            &acetate[..],
            vec![
                bond("C", "CH3", BondOrder::Single, template("ACT")),
                bond("C", "O", BondOrder::Double, template("ACT")),
                bond("C", "OXT", BondOrder::Single, template("ACT")),
            ],
            vec![("OXT".to_owned(), -1)],
            vec![],
            3,
        ),
    ];
    for (lines, expected_bonds, expected_charges, expected_issues, applied) in cases {
        let interpreted = interpret(&entry(lines));
        assert_eq!(bonds(&interpreted), expected_bonds);
        assert_eq!(charges(&interpreted), expected_charges);
        assert_eq!(interpreted.report().issues(), expected_issues);
        assert_eq!(interpreted.report().applied_connections(), applied);
        assert_eq!(interpreted.topology().instance_count(), 1);
    }
}

#[test]
fn link_and_ssbond_records_bond_their_partners() {
    let cysteine_a = || cysteine('A', ["1", "2", "3", "4", "5", "6"], 0.0, 1.0);
    let zinc = site(
        "HETATM",
        "7",
        "ZN  ",
        "ZN",
        'A',
        "2",
        [6.1, -0.8, 1.2],
        "ZN",
    )
    .line();
    let zinc_link = |symmetry| {
        link(
            [" SG", "ZN"],
            [("CYS", 'A', "1"), ("ZN", 'A', "2")],
            ["1555", symmetry],
            2.30,
        )
    };
    let cases = [
        // A link to a metal is coordination; a CONECT record repeating it,
        // as wwPDB files write, adds nothing.
        (
            [
                cysteine_a(),
                vec![
                    zinc.clone(),
                    zinc_link("1555"),
                    "CONECT    6    7".to_owned(),
                    "CONECT    7    6".to_owned(),
                ],
            ]
            .concat(),
            Some(("SG", "ZN", BondOrder::Dative, connection(None, "LINK"))),
            1,
            vec![],
        ),
        (
            [
                cysteine_a(),
                cysteine('B', ["7", "8", "9", "10", "11", "12"], 2.04, -1.0),
                vec![ssbond("1", [("CYS", 'A', "1"), ("CYS", 'B', "1")], 2.04)],
            ]
            .concat(),
            Some((
                "SG",
                "SG",
                BondOrder::Single,
                connection(Some("1"), "SSBOND"),
            )),
            1,
            vec![],
        ),
        // Links to symmetry mates are never bonds.
        (
            [cysteine_a(), vec![zinc.clone(), zinc_link("2555")]].concat(),
            None,
            2,
            vec![StructureInterpretIssue::ConnectionUnresolved {
                connection_id: None,
                connection_type: "LINK".to_owned(),
                partner: 2,
                source_line: Some(8),
                reason: ConnectionResolutionReason::UnsupportedSymmetry {
                    symmetry: "2555".to_owned(),
                },
            }],
        ),
    ];
    for (lines, linked, molecules, issues) in cases {
        let interpreted = interpret(&entry(&lines));
        let connections = bonds(&interpreted)
            .into_iter()
            .filter(|bond| matches!(bond.3, BondEvidence::Connection { .. }))
            .collect::<Vec<_>>();
        let expected = linked
            .clone()
            .map(|(a, b, order, evidence)| (a.to_owned(), b.to_owned(), order, evidence))
            .into_iter()
            .collect::<Vec<_>>();
        assert_eq!(connections, expected);
        assert_eq!(interpreted.topology().instance_count(), molecules);
        assert_eq!(
            interpreted.report().applied_connections(),
            usize::from(linked.is_some())
        );
        assert_eq!(interpreted.report().issues(), issues);
    }
}

#[test]
fn coordinate_models_become_ensemble_members() {
    let text = entry(
        &[
            vec!["MODEL        1".to_owned()],
            water(0.0),
            vec!["ENDMDL".to_owned(), "MODEL        2".to_owned()],
            water(5.0),
            vec!["ENDMDL".to_owned()],
        ]
        .concat(),
    );
    let document = pdb::parse_str(&text, Default::default()).unwrap();
    let serials = document
        .models()
        .iter()
        .map(|model| model.serial())
        .collect::<Vec<_>>();
    assert_eq!(serials, [1, 2]);
    let error = document.interpret(Default::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);

    let second = document
        .interpret(
            StructureInterpretOptions::default()
                .with_model_selection(ModelSelection::Select("2".to_owned())),
        )
        .unwrap();
    assert_eq!(second.report().selected_model(), Some("2"));
    assert_eq!(second.report().ignored_coordinate_models(), ["1"]);
    assert_eq!(second.topology().atom_count(), 3);
    assert_eq!(second.topology().bond_count(), 2);

    let ensemble = document.interpret_ensemble(Default::default()).unwrap();
    assert_eq!(ensemble.ensemble().len(), 2);
    let models = ensemble
        .reports()
        .iter()
        .map(|report| report.selected_model())
        .collect::<Vec<_>>();
    assert_eq!(models, [Some("1"), Some("2")]);
}

#[test]
fn recurring_numbers_and_alternate_locations_identify_residues() {
    // Residue numbers wrap around: the third water reuses number 1, and its
    // oxygen has two alternate locations. CONECT names serial 3, the record
    // of label A.
    let alternative = |serial, alternate_location, occupancy, x| {
        Site {
            alternate_location,
            occupancy,
            ..site(
                "HETATM",
                serial,
                " O  ",
                "HOH",
                'W',
                "1",
                [x, 10.0, 0.0],
                "O",
            )
        }
        .line()
    };
    let text = entry(&[
        site("HETATM", "1", " O  ", "HOH", 'W', "1", [0.0, 0.0, 0.0], "O").line(),
        site("HETATM", "2", " O  ", "HOH", 'W', "2", [5.0, 0.0, 0.0], "O").line(),
        alternative("3", 'A', 0.4, 10.0),
        alternative("4", 'B', 0.6, 12.0),
        "CONECT    3    1".to_owned(),
    ]);
    let document = pdb::parse_str(&text, Default::default()).unwrap();
    let repeated = residue(
        "W",
        SourceResiduePosition::RepeatedAuthor {
            sequence: "1".to_owned(),
            repeat: 1,
        },
    );
    let alternatives = document.alternate_locations().unwrap();
    let [alternative] = alternatives.as_slice() else {
        panic!("one residue with alternate locations");
    };
    assert_eq!(alternative.residue(), &repeated);
    assert_eq!(alternative.labels(), ["A", "B"]);

    let cases = [
        (
            AltLocPolicy::HighestOccupancy,
            "B",
            vec![ConnectionResolutionReason::AlternateLocationOmitted {
                alternate_location: "A".to_owned(),
            }],
            0,
        ),
        (AltLocPolicy::SelectLabel("A".to_owned()), "A", vec![], 1),
    ];
    for (policy, label, unresolved, bonds) in cases {
        let interpreted = document
            .interpret(StructureInterpretOptions::default().with_altloc_policy(policy))
            .unwrap();
        let [decision] = interpreted.report().alternate_locations() else {
            panic!("one alternate-location decision");
        };
        assert_eq!(decision.residue(), &repeated);
        assert_eq!(decision.selected_labels(), [label]);
        let mut residues = interpreted
            .report()
            .instances()
            .iter()
            .flat_map(|instance| instance.atoms())
            .map(|atom| (atom.author_sequence_id(), atom.occurrence()))
            .collect::<Vec<_>>();
        residues.sort();
        assert_eq!(
            residues,
            [(Some("1"), None), (Some("1"), Some(1)), (Some("2"), None)]
        );
        let reasons = interpreted
            .report()
            .issues()
            .iter()
            .filter_map(|issue| match issue {
                StructureInterpretIssue::ConnectionUnresolved { reason, .. } => {
                    Some(reason.clone())
                }
                _ => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(reasons, unresolved);
        assert_eq!(interpreted.topology().bond_count(), bonds);
    }
}

#[test]
fn a_ter_record_ends_the_residue_it_follows() {
    // Records after a TER record start a new residue even when they repeat
    // the chain, number, and insertion code of the residue it ends.
    let calpha = |serial, x| site("ATOM", serial, " CA ", "GLY", 'A', "1", [x, 0.0, 0.0], "C");
    let single = vec![
        calpha("1", 0.0).line(),
        ter("2", "GLY", 'A', "1"),
        calpha("3", 10.0).line(),
    ];
    // The second dipeptide, numbered 2-3 in records 10-17, restarts at the
    // number that ends the first.
    let restarted = glycine_dipeptide()
        .into_iter()
        .map(|line| {
            let serial = line[6..11].trim().parse::<u32>().unwrap() + 9;
            let sequence = line[22..26].trim().parse::<i32>().unwrap() + 1;
            let x = line[30..38].trim().parse::<f64>().unwrap() + 20.0;
            format!(
                "{}{serial:>5}{}{sequence:>4}{}{x:>8.3}{}",
                &line[..6],
                &line[11..22],
                &line[26..30],
                &line[38..]
            )
        })
        .collect::<Vec<_>>();
    let segments = [
        glycine_dipeptide(),
        vec![ter("9", "GLY", 'A', "2")],
        restarted,
        vec![ter("18", "GLY", 'A', "3")],
    ]
    .concat();
    // Each instance's residues: author number, occurrence, and source lines.
    type Residues<'a> = Vec<(&'a str, Option<usize>, Vec<usize>)>;
    let cases: [(Vec<String>, EntityKind, Vec<Residues>, usize); 2] = [
        (
            single.clone(),
            EntityKind::NonPolymer,
            vec![vec![("1", None, vec![1])], vec![("1", Some(1), vec![3])]],
            0,
        ),
        (
            segments,
            EntityKind::Polymer,
            vec![
                vec![("1", None, vec![1, 2, 3, 4]), ("2", None, vec![5, 6, 7, 8])],
                vec![
                    ("2", Some(1), vec![10, 11, 12, 13]),
                    ("3", None, vec![14, 15, 16, 17]),
                ],
            ],
            2,
        ),
    ];
    for (lines, kind, expected, peptide_bonds) in cases {
        let interpreted = interpret(&entry(&lines));
        let instances = interpreted
            .report()
            .instances()
            .iter()
            .map(|instance| {
                assert_eq!(instance.entity_kinds(), std::slice::from_ref(&kind));
                let mut residues: Residues = Vec::new();
                for atom in instance.atoms() {
                    let sequence = atom.author_sequence_id().unwrap();
                    match residues.last_mut() {
                        Some((last, occurrence, lines))
                            if *last == sequence && *occurrence == atom.occurrence() =>
                        {
                            lines.push(atom.source_line());
                        }
                        _ => residues.push((sequence, atom.occurrence(), vec![atom.source_line()])),
                    }
                }
                residues
            })
            .collect::<Vec<_>>();
        assert_eq!(instances, expected);
        let hierarchy = interpreted.topology().hierarchy();
        let residues = hierarchy
            .residues()
            .map(|(_, residue)| {
                (
                    residue.name().to_owned(),
                    residue.author_seq_id().map(str::to_owned),
                    residue.atom_sites().len(),
                )
            })
            .collect::<Vec<_>>();
        let sizes = expected
            .iter()
            .flatten()
            .map(|(sequence, _, lines)| {
                ("GLY".to_owned(), Some((*sequence).to_owned()), lines.len())
            })
            .collect::<Vec<_>>();
        assert_eq!(residues, sizes);
        assert_eq!(hierarchy.chains().count(), 1);
        assert_eq!(
            bonds_with_evidence(&interpreted, &BondEvidence::PolymerLink),
            peptide_bonds
        );
        assert_eq!(interpreted.report().issues(), &[]);
    }

    // Every coordinate model numbers its repeats afresh.
    let models = entry(
        &[
            vec!["MODEL        1".to_owned()],
            single.clone(),
            vec!["ENDMDL".to_owned(), "MODEL        2".to_owned()],
            single,
            vec!["ENDMDL".to_owned()],
        ]
        .concat(),
    );
    let ensemble = (pdb::parse_str(&models, Default::default()).unwrap())
        .interpret_ensemble(Default::default())
        .unwrap();
    assert_eq!(ensemble.ensemble().len(), 2);
    let occurrences = ensemble
        .reports()
        .iter()
        .map(|report| {
            report
                .instances()
                .iter()
                .flat_map(|instance| instance.atoms())
                .map(|atom| (atom.source_line(), atom.occurrence()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(
        occurrences,
        [[(2, None), (4, Some(1))], [(7, None), (9, Some(1))]]
    );
}

#[test]
fn atom_names_imply_elements_without_element_columns() {
    let text = entry(&[
        site("ATOM", "1", " CA ", "GLY", 'A', "1", [0.0, 0.0, 0.0], "").coordinates_only(),
        site("HETATM", "2", "CA  ", "CA", 'A', "2", [5.0, 0.0, 0.0], "").coordinates_only(),
        site(
            "HETATM",
            "3",
            "ZN  ",
            "ZN",
            'A',
            "3",
            [10.0, 0.0, 0.0],
            "ZN",
        )
        .line(),
    ]);
    let interpreted = interpret(&text);
    let elements = interpreted
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .map(|atom| {
            let molecule = interpreted
                .topology()
                .instance(atom.atom().instance())
                .unwrap()
                .molecule();
            molecule
                .atom(atom.atom().atom())
                .unwrap()
                .element()
                .symbol()
        })
        .collect::<Vec<_>>();
    assert_eq!(elements, ["C", "Ca", "Zn"]);
    assert_eq!(
        interpreted.report().issues(),
        &[StructureInterpretIssue::ElementsInferredFromAtomNames {
            source_lines: vec![1, 2],
        }]
    );
}

/// Atom names, element symbols, and mass numbers, in instance order.
fn isotopes(result: &StructureInterpretation) -> Vec<(String, &'static str, Option<u16>)> {
    result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .map(|atom| {
            let molecule = result
                .topology()
                .instance(atom.atom().instance())
                .unwrap()
                .molecule();
            let element = molecule.atom(atom.atom().atom()).unwrap();
            (
                atom.atom_name().to_owned(),
                element.element().symbol(),
                element.isotope().map(std::num::NonZero::get),
            )
        })
        .collect()
}

#[test]
fn hydrogen_isotope_symbols_read_as_isotopes_and_write_back() {
    // Heavy water as neutron structures write it, with the CCD `DOD` names.
    let text = entry(
        &[
            ("1", " O  ", [0.0, 0.0, 0.0], "O"),
            ("2", " D1 ", [0.96, 0.0, 0.0], "D"),
            ("3", " D2 ", [-0.24, 0.93, 0.0], "D"),
        ]
        .map(|(serial, name, position, element)| {
            site("HETATM", serial, name, "DOD", 'W', "1", position, element).line()
        }),
    );
    let interpreted = interpret(&text);
    let expected = [
        ("O".to_owned(), "O", None),
        ("D1".to_owned(), "H", Some(2)),
        ("D2".to_owned(), "H", Some(2)),
    ];
    assert_eq!(isotopes(&interpreted), expected);
    let template = BondEvidence::Template {
        component_id: "DOD".to_owned(),
        variant: None,
    };
    assert_eq!(
        bonds(&interpreted),
        [
            (
                "D1".to_owned(),
                "O".to_owned(),
                BondOrder::Single,
                template.clone()
            ),
            ("D2".to_owned(), "O".to_owned(), BondOrder::Single, template),
        ]
    );

    let written =
        kekule::mmcif::write([&interpreted], kekule::mmcif::MmcifWriteOptions::default()).unwrap();
    let type_symbols = written
        .lines()
        .filter(|line| line.starts_with("HETATM"))
        .map(|line| line.split_whitespace().nth(2).unwrap())
        .collect::<Vec<_>>();
    assert_eq!(type_symbols, ["O", "D", "D"]);
    let reread = kekule::mmcif::parse_str(&written, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(isotopes(&reread), expected);
}

#[test]
fn author_chains_holding_ligands_and_water_write_as_separate_mmcif_asyms() {
    // Chain A holds a glycine dipeptide, a sulfate, a zinc that a `LINK`
    // coordinates to the second glycine's carbonyl oxygen, and a water.
    let mut lines = glycine_dipeptide();
    lines.push(ter("9", "GLY", 'A', "2"));
    for (serial, name, position, element) in [
        ("10", " S  ", [10.0, 0.0, 0.0], "S"),
        ("11", " O1 ", [10.849, 0.849, 0.849], "O"),
        ("12", " O2 ", [9.151, -0.849, 0.849], "O"),
        ("13", " O3 ", [9.151, 0.849, -0.849], "O"),
        ("14", " O4 ", [10.849, -0.849, -0.849], "O"),
    ] {
        lines.push(site("HETATM", serial, name, "SO4", 'A', "101", position, element).line());
    }
    lines.push(
        site(
            "HETATM",
            "15",
            "ZN  ",
            "ZN",
            'A',
            "102",
            [7.152, -0.136, 0.0],
            "ZN",
        )
        .line(),
    );
    lines.push(
        site(
            "HETATM",
            "16",
            " O  ",
            "HOH",
            'A',
            "201",
            [20.0, 0.0, 0.0],
            "O",
        )
        .line(),
    );
    lines.push(link(
        [" O  ", "ZN  "],
        [("GLY", 'A', "2"), ("ZN", 'A', "102")],
        ["1555", "1555"],
        2.10,
    ));
    let interpreted = interpret(&entry(&lines));
    let source_bonds = bonds(&interpreted)
        .into_iter()
        .map(|(a, b, order, _)| (a, b, order))
        .collect::<Vec<_>>();
    assert!(source_bonds.contains(&("O".to_owned(), "ZN".to_owned(), BondOrder::Dative)));

    // Residues by component, label asym, author asym, entity, and kind.
    let residues = |result: &StructureInterpretation| {
        let mut atoms = result
            .report()
            .instances()
            .iter()
            .flat_map(|instance| instance.atoms())
            .collect::<Vec<_>>();
        atoms.sort_by_key(|atom| atom.source_line());
        let mut residues = atoms
            .into_iter()
            .map(|atom| {
                (
                    atom.component_id().to_owned(),
                    atom.asym_id().to_owned(),
                    atom.auth_asym_id().map(str::to_owned),
                    atom.entity_id().map(str::to_owned),
                    atom.entity_kind().clone(),
                )
            })
            .collect::<Vec<_>>();
        residues.dedup();
        residues
    };
    let expected = |entities: [&str; 4]| {
        [
            ("GLY", "A", EntityKind::Polymer),
            ("SO4", "A_2", EntityKind::NonPolymer),
            ("ZN", "A_3", EntityKind::NonPolymer),
            ("HOH", "A_4", EntityKind::Water),
        ]
        .into_iter()
        .zip(entities)
        .map(|((component, asym, kind), entity)| {
            (
                component.to_owned(),
                asym.to_owned(),
                Some("A".to_owned()),
                Some(entity.to_owned()),
                kind,
            )
        })
        .collect::<Vec<_>>()
    };
    // Source provenance and canonical classification plan the same asyms;
    // entities without source identifiers are numbered by each plan.
    for (written, entities) in [
        (
            kekule::mmcif::write([&interpreted], kekule::mmcif::MmcifWriteOptions::default()),
            ["K1", "K2", "K3", "K4"],
        ),
        (
            kekule::mmcif::write(
                [interpreted.model()],
                kekule::mmcif::MmcifWriteOptions::default(),
            ),
            ["1", "2", "3", "4"],
        ),
    ] {
        let reread = kekule::mmcif::parse_str(&written.unwrap(), Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap();
        assert_eq!(residues(&reread), expected(entities));
        let reread_bonds = bonds(&reread)
            .into_iter()
            .map(|(a, b, order, _)| (a, b, order))
            .collect::<Vec<_>>();
        assert_eq!(reread_bonds, source_bonds);
    }
}

#[test]
fn pdb_parse_errors_locate_the_record() {
    let oxygen = site("HETATM", "1", " O  ", "HOH", 'W', "1", [0.0, 0.0, 0.0], "O");
    let water = oxygen.line();
    let mut invalid_coordinate = water.clone();
    invalid_coordinate.replace_range(30..38, "   1.2.3");
    let invalid_sequence = Site {
        sequence: "1 2",
        ..oxygen
    }
    .line();
    let invalid_charge = Site {
        charge: "XX",
        ..oxygen
    }
    .line();
    let cases = [
        (
            invalid_coordinate,
            DiagnosticKind::InvalidValue,
            SourceLocation::at_column(1, 31),
        ),
        (
            invalid_sequence,
            DiagnosticKind::InvalidValue,
            SourceLocation::at_column(1, 23),
        ),
        (
            invalid_charge,
            DiagnosticKind::InvalidValue,
            SourceLocation::at_column(1, 79),
        ),
        (
            format!("{water}\nMODEL        1\n{water}\nENDMDL"),
            DiagnosticKind::Syntax,
            SourceLocation::at_line(1),
        ),
        (
            format!("MODEL        1\n{water}\nENDMDL\n{water}"),
            DiagnosticKind::Syntax,
            SourceLocation::at_line(4),
        ),
        (
            format!("MODEL        1\n{water}\nENDMDL\nMODEL        1\nENDMDL"),
            DiagnosticKind::Inconsistent,
            SourceLocation::at_line(4),
        ),
        (
            format!("{water}\nEND\n{water}"),
            DiagnosticKind::Unsupported,
            SourceLocation::at_line(3),
        ),
        (
            water.replace("HOH", "HÖH"),
            DiagnosticKind::Syntax,
            SourceLocation::at_line(1),
        ),
    ];
    for (text, kind, location) in cases {
        let error = pdb::parse_str(&text, Default::default()).unwrap_err();
        assert_eq!(
            (error.kind(), error.location()),
            (kind, Some(location)),
            "{text}"
        );
    }

    let limited = pdb::parse_str(
        &format!("{water}\n{water}"),
        PdbParseOptions::default().with_max_atom_records(1),
    )
    .unwrap_err();
    assert_eq!(
        (limited.kind(), limited.location()),
        (DiagnosticKind::Limit, Some(SourceLocation::at_line(2)))
    );
}

#[test]
fn pdb_interpretation_rejects_contradictory_records() {
    let water = water(0.0).join("\n");
    let unknown = site("HETATM", "1", " X  ", "UNX", 'X', "1", [0.0, 0.0, 0.0], "X").line();
    let cases = [
        (
            format!("{water}\nCONECT    1    1"),
            StructureInterpretOptions::default(),
            DiagnosticKind::InvalidValue,
            Some(4),
        ),
        (
            format!("{water}\nCONECT    1    2    2    2    2"),
            StructureInterpretOptions::default(),
            DiagnosticKind::InvalidValue,
            Some(4),
        ),
        (
            unknown,
            StructureInterpretOptions::default(),
            DiagnosticKind::InvalidValue,
            Some(1),
        ),
        (
            water.clone(),
            StructureInterpretOptions::default().with_entity_metadata(
                kekule::macromolecular::EntityMetadataPolicy::RequireDeclared,
            ),
            DiagnosticKind::Inconsistent,
            None,
        ),
    ];
    for (text, options, kind, line) in cases {
        let document = pdb::parse_str(&text, Default::default()).unwrap();
        let error = document.interpret(options).unwrap_err();
        assert_eq!(
            (
                error.kind(),
                error.location().map(|location| location.line())
            ),
            (kind, line),
            "{text}"
        );
    }
}

/// One atom: the aligned name field, element symbol, and position.
type Atom = (&'static str, &'static str, [f64; 3]);

/// Records of one residue, shifted along x by `offset`.
fn residue_records(
    record: &'static str,
    residue: &'static str,
    chain: char,
    atoms: &[Atom],
    offset: f64,
) -> Vec<String> {
    atoms
        .iter()
        .map(|&(name, element, [x, y, z])| {
            site(
                record,
                "1",
                name,
                residue,
                chain,
                "1",
                [x + offset, y, z],
                element,
            )
            .line()
        })
        .collect()
}

fn inferred_components(result: &StructureInterpretation) -> Vec<String> {
    result
        .report()
        .issues()
        .iter()
        .filter_map(|issue| match issue {
            StructureInterpretIssue::ResidueBondsInferred { component_id, .. } => {
                Some(component_id.clone())
            }
            _ => None,
        })
        .collect()
}

#[test]
fn force_field_libraries_interpret_their_residue_names() {
    use kekule::components::ComponentLibrary;
    use kekule::macromolecular::ConnectivityOptions;
    use std::sync::Arc;

    // Positions are the ideal coordinates of the dictionary's HIS, ASP, and
    // ILE entries, the latter with its atoms under CHARMM names.
    const HISTIDINE: &[Atom] = &[
        (" N  ", "N", [-0.040, -1.210, 0.053]),
        (" CA ", "C", [1.172, -1.709, 0.652]),
        (" C  ", "C", [1.083, -3.207, 0.905]),
        (" O  ", "O", [0.040, -3.770, 1.222]),
        (" CB ", "C", [1.484, -0.975, 1.962]),
        (" CG ", "C", [2.940, -1.060, 2.353]),
        (" ND1", "N", [3.380, -2.075, 3.129]),
        (" CD2", "C", [3.960, -0.251, 2.046]),
        (" CE1", "C", [4.693, -1.908, 3.317]),
        (" NE2", "N", [5.058, -0.801, 2.662]),
        (" OXT", "O", [2.247, -3.882, 0.744]),
        (" H  ", "H", [-0.102, -1.155, -0.950]),
        (" H2 ", "H", [-0.715, -0.741, 0.634]),
        (" HA ", "H", [1.965, -1.558, -0.089]),
        (" HB2", "H", [1.215, 0.087, 1.879]),
        (" HB3", "H", [0.859, -1.368, 2.775]),
        (" HD2", "H", [4.108, 0.647, 1.479]),
        (" HE1", "H", [5.340, -2.550, 3.892]),
        (" HE2", "H", [6.002, -0.428, 2.627]),
        (" HXT", "H", [2.188, -4.848, 0.901]),
    ];
    const ASPARTATE: &[Atom] = &[
        (" N  ", "N", [-0.317, 1.688, 0.066]),
        (" CA ", "C", [-0.470, 0.286, -0.344]),
        (" C  ", "C", [-1.868, -0.180, -0.029]),
        (" O  ", "O", [-2.534, 0.415, 0.786]),
        (" CB ", "C", [0.539, -0.580, 0.413]),
        (" CG ", "C", [1.938, -0.195, 0.004]),
        (" OD1", "O", [2.109, 0.681, -0.810]),
        (" OD2", "O", [2.992, -0.826, 0.543]),
        (" OXT", "O", [-2.374, -1.256, -0.652]),
    ];
    // A zwitterionic isoleucine with CHARMM's N- and C-terminal patches.
    const ISOLEUCINE: &[Atom] = &[
        (" N  ", "N", [-1.944, 0.335, -0.343]),
        (" CA ", "C", [-0.487, 0.519, -0.369]),
        (" C  ", "C", [0.066, -0.032, -1.657]),
        (" OT1", "O", [-0.484, -0.958, -2.203]),
        (" OT2", "O", [1.171, 0.504, -2.197]),
        (" CB ", "C", [0.140, -0.219, 0.814]),
        (" CG1", "C", [-0.421, 0.341, 2.122]),
        (" CG2", "C", [1.658, -0.027, 0.788]),
        (" CD ", "C", [0.206, -0.397, 3.305]),
        (" HA ", "H", [-0.253, 1.582, -0.299]),
        (" HB ", "H", [-0.092, -1.281, 0.744]),
        ("HG11", "H", [-1.502, 0.204, 2.141]),
        ("HG12", "H", [-0.188, 1.403, 2.192]),
        ("HG21", "H", [1.891, 1.034, 0.857]),
        ("HG22", "H", [2.105, -0.554, 1.631]),
        ("HG23", "H", [2.059, -0.427, -0.143]),
        (" HD1", "H", [-0.193, 0.001, 4.237]),
        (" HD2", "H", [-0.026, -1.460, 3.235]),
        (" HD3", "H", [1.287, -0.261, 3.286]),
        (" HT1", "H", [-2.112, -0.656, -0.410]),
        (" HT2", "H", [-2.256, 0.622, 0.572]),
        (" HT3", "H", [-2.300, 0.950, -1.100]),
    ];
    let with = |library: ComponentLibrary| {
        StructureInterpretOptions::default().with_connectivity(
            ConnectivityOptions::default().with_components(Some(Arc::new(library))),
        )
    };

    // AMBER: the NE2-protonated histidine, an aspartate without hydrogens,
    // water, and a sodium ion.
    let amber = entry(
        &[
            residue_records("ATOM", "HIE", 'A', HISTIDINE, 0.0),
            residue_records("ATOM", "ASP", 'B', ASPARTATE, 20.0),
            residue_records(
                "HETATM",
                "WAT",
                'W',
                &[
                    (" O  ", "O", [0.0, 0.0, 0.0]),
                    (" H1 ", "H", [0.96, 0.0, 0.0]),
                    (" H2 ", "H", [-0.24, 0.93, 0.0]),
                ],
                40.0,
            ),
            residue_records(
                "HETATM",
                "Na+",
                'I',
                &[("Na+ ", "NA", [0.0, 0.0, 0.0])],
                60.0,
            ),
        ]
        .concat(),
    );
    let document = pdb::parse_str(&amber, Default::default()).unwrap();
    let templated = document.interpret(with(ComponentLibrary::amber())).unwrap();
    assert_eq!(templated.topology().bond_count(), 30);
    assert_eq!(
        charges(&templated),
        [("OD2".to_owned(), -1), ("Na+".to_owned(), 1)]
    );
    assert_eq!(templated.report().issues(), &[]);
    // The dictionary names neither HIE nor WAT, and its aspartic acid is neutral.
    let untemplated = document.interpret(Default::default()).unwrap();
    assert_eq!(inferred_components(&untemplated), ["HIE", "WAT"]);
    assert_eq!(charges(&untemplated), []);

    // CHARMM: the isoleucine, water, and sodium and chloride ions. Terminal
    // hydrogens HT1-HT3 attach by distance.
    let charmm = entry(
        &[
            residue_records("ATOM", "ILE", 'A', ISOLEUCINE, 0.0),
            residue_records(
                "HETATM",
                "TIP3",
                'W',
                &[
                    (" OH2", "O", [0.0, 0.0, 0.0]),
                    (" H1 ", "H", [0.96, 0.0, 0.0]),
                    (" H2 ", "H", [-0.24, 0.93, 0.0]),
                ],
                20.0,
            ),
            residue_records(
                "HETATM",
                "SOD",
                'I',
                &[("SOD ", "NA", [0.0, 0.0, 0.0])],
                40.0,
            ),
            residue_records(
                "HETATM",
                "CLA",
                'J',
                &[("CLA ", "CL", [0.0, 0.0, 0.0])],
                45.0,
            ),
        ]
        .concat(),
    );
    let document = pdb::parse_str(&charmm, Default::default()).unwrap();
    let templated = document
        .interpret(with(ComponentLibrary::charmm()))
        .unwrap();
    assert_eq!(templated.topology().bond_count(), 23);
    assert_eq!(
        charges(&templated),
        [
            ("N".to_owned(), 1),
            ("OT2".to_owned(), -1),
            ("SOD".to_owned(), 1),
            ("CLA".to_owned(), -1),
        ]
    );
    assert_eq!(inferred_components(&templated), Vec::<String>::new());
    let attached = templated
        .report()
        .issues()
        .iter()
        .filter_map(|issue| match issue {
            StructureInterpretIssue::HydrogenAttachedByDistance {
                atom_name,
                heavy_atom_name,
                ..
            } => Some((atom_name.as_str(), heavy_atom_name.as_str())),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(attached, [("HT1", "N"), ("HT2", "N"), ("HT3", "N")]);
    // Without the library, CHARMM chloride meets the dictionary's chlorophyll a.
    let error = document.interpret(Default::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
}

#[test]
fn force_field_termini_carry_their_charges() {
    use kekule::components::ComponentLibrary;
    use kekule::macromolecular::ConnectivityOptions;
    use std::sync::Arc;

    // Ideal coordinates of the dictionary's PRO and DA entries; the third
    // ammonium hydrogen of proline and the 5'-hydroxyl hydrogen are placed
    // 1.0 A from their heavy atoms.
    let proline = |names: [&'static str; 9]| -> Vec<Atom> {
        let [h2, h3, hb2, hb3, hg2, hg3, hd2, hd3, oxt] = names;
        vec![
            (" N  ", "N", [-0.816, 1.108, 0.254]),
            (" CA ", "C", [0.001, -0.107, 0.509]),
            (" C  ", "C", [1.408, 0.091, 0.005]),
            (" O  ", "O", [1.650, 0.980, -0.777]),
            (" CB ", "C", [-0.703, -1.227, -0.286]),
            (" CG ", "C", [-2.163, -0.753, -0.439]),
            (" CD ", "C", [-2.218, 0.614, 0.276]),
            (oxt, "O", [2.391, -0.721, 0.424]),
            (h2, "H", [-0.707, 1.708, 1.057]),
            (h3, "H", [-0.567, 1.588, -0.599]),
            (" HA ", "H", [0.009, -0.343, 1.573]),
            (hb2, "H", [-0.240, -1.345, -1.266]),
            (hb3, "H", [-0.666, -2.165, 0.267]),
            (hg2, "H", [-2.416, -0.638, -1.493]),
            (hg3, "H", [-2.843, -1.458, 0.040]),
            (hd2, "H", [-2.872, 1.300, -0.263]),
            (hd3, "H", [-2.559, 0.492, 1.304]),
        ]
    };
    let mut amber_proline = proline([
        " H2 ", " H3 ", " HB2", " HB3", " HG2", " HG3", " HD2", " HD3", " OXT",
    ]);
    amber_proline.push((" HXT", "H", [3.293, -0.594, 0.101]));
    let charmm_proline = proline([
        " HN1", " HN2", " HB1", " HB2", " HG1", " HG2", " HD1", " HD2", " OT2",
    ]);
    // A free deoxyadenosine: no phosphate, and hydroxyl hydrogens at 5' and 3'.
    const DEOXYADENOSINE: &[Atom] = &[
        (" O5'", "O", [0.241, -0.771, -3.320]),
        (" C5'", "C", [-0.549, 0.270, -2.744]),
        (" C4'", "C", [-1.239, -0.251, -1.482]),
        (" O4'", "O", [-0.267, -0.564, -0.458]),
        (" C3'", "C", [-2.105, 0.859, -0.835]),
        (" O3'", "O", [-3.409, 0.895, -1.418]),
        (" C2'", "C", [-2.173, 0.398, 0.640]),
        (" C1'", "C", [-0.965, -0.545, 0.797]),
        (" N9 ", "N", [-0.078, -0.047, 1.852]),
        (" C8 ", "C", [0.962, 0.817, 1.689]),
        (" N7 ", "N", [1.535, 1.044, 2.835]),
        (" C5 ", "C", [0.897, 0.346, 3.805]),
        (" C6 ", "C", [1.069, 0.196, 5.191]),
        (" N6 ", "N", [2.079, 0.869, 5.856]),
        (" N1 ", "N", [0.236, -0.603, 5.850]),
        (" C2 ", "C", [-0.729, -1.249, 5.224]),
        (" N3 ", "N", [-0.925, -1.144, 3.927]),
        (" C4 ", "C", [-0.142, -0.368, 3.184]),
        ("HO5'", "H", [0.659, -0.400, -4.113]),
        (" H5'", "H", [-1.302, 0.594, -3.463]),
        ("H5''", "H", [0.092, 1.112, -2.486]),
        (" H4'", "H", [-1.846, -1.126, -1.712]),
        (" H3'", "H", [-1.617, 1.830, -0.918]),
        ("HO3'", "H", [-3.924, 1.538, -0.913]),
        (" H2'", "H", [-3.103, -0.136, 0.831]),
        ("H2''", "H", [-2.079, 1.253, 1.311]),
        (" H1'", "H", [-1.309, -1.549, 1.046]),
        (" H8 ", "H", [1.266, 1.250, 0.748]),
        (" H61", "H", [2.185, 0.761, 6.814]),
        (" H62", "H", [2.683, 1.447, 5.363]),
        (" H2 ", "H", [-1.383, -1.889, 5.798]),
    ];
    let cases = [
        (
            ComponentLibrary::amber(),
            residue_records("ATOM", "DAN", 'D', DEOXYADENOSINE, 0.0),
            vec![],
            vec![],
        ),
        (
            ComponentLibrary::amber(),
            residue_records("ATOM", "NPRO", 'P', &amber_proline, 0.0),
            vec![("N".to_owned(), 1)],
            vec![],
        ),
        // CHARMM's N-terminal proline patch names its hydrogens HN1 and HN2;
        // the C-terminal patch makes a carboxylate.
        (
            ComponentLibrary::charmm(),
            residue_records("ATOM", "PRO", 'P', &charmm_proline, 0.0),
            vec![("N".to_owned(), 1), ("OT2".to_owned(), -1)],
            vec![("HN2".to_owned(), "N".to_owned())],
        ),
    ];
    for (library, lines, expected_charges, expected_attached) in cases {
        let document = pdb::parse_str(&entry(&lines), Default::default()).unwrap();
        let interpreted = document
            .interpret(StructureInterpretOptions::default().with_connectivity(
                ConnectivityOptions::default().with_components(Some(Arc::new(library))),
            ))
            .unwrap();
        let attached = interpreted
            .report()
            .issues()
            .iter()
            .filter_map(|issue| match issue {
                StructureInterpretIssue::HydrogenAttachedByDistance {
                    atom_name,
                    heavy_atom_name,
                    ..
                } => Some((atom_name.clone(), heavy_atom_name.clone())),
                _ => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(charges(&interpreted), expected_charges, "{lines:?}");
        assert_eq!(attached, expected_attached, "{lines:?}");
        assert_eq!(inferred_components(&interpreted), Vec::<String>::new());
        assert_eq!(interpreted.topology().instance_count(), 1);
    }
}

const CRAMBIN: &str = include_str!("../fixtures/pdb/1CRN.pdb");
const ADENYLATE_KINASE: &str = include_str!("../fixtures/pdb/1AKE.pdb");

#[test]
fn pdb_connectivity_matches_an_independent_reference() {
    // Atom, bond, and molecule counts from Biotite 1.2.0
    // (`get_structure(include_bonds=True, altloc="occupancy")`) on the mmCIF
    // files of the same entries, as `mmcif_connectivity_matches_an_independent_reference`
    // asserts. Biotite's PDB reader gives the same counts for 1CRN; for 1AKE it
    // maps `CONECT` serials of the dropped alternate locations of AP5 to the
    // first atom of chain A, adding a MET 1 N self-bond and two bonds 23 A long.
    // Charged atoms count the LYS NZ, ARG NH2, and HIS ND1 sites the dictionary
    // charges (+1), counted from the files. The bundled templates do not cover
    // AP5: its 62 bonds per chain come from CONECT records, none inferred, and
    // the 12 CONECT partners in the omitted alternate location B are reported.
    for (name, text, atoms, bonds, largest, molecules, charged, listed, disulfides, omitted) in [
        ("1CRN", CRAMBIN, 327, 337, vec![327], 1, 2, 0, 3, 0),
        (
            "1AKE",
            ADENYLATE_KINASE,
            3804,
            3484,
            vec![1656, 1656, 57, 57],
            382,
            68,
            124,
            0,
            12,
        ),
    ] {
        let interpreted = interpret(text);
        let topology = interpreted.topology();
        let mut sizes = topology
            .instances()
            .map(|instance| instance.molecule().atom_count())
            .collect::<Vec<_>>();
        sizes.sort_unstable_by(|a, b| b.cmp(a));
        assert_eq!(topology.atom_count(), atoms, "{name}");
        assert_eq!(topology.bond_count(), bonds, "{name}");
        assert_eq!(sizes.len(), molecules, "{name}");
        assert_eq!(&sizes[..largest.len()], largest.as_slice(), "{name}");
        let charges = charges(&interpreted)
            .into_iter()
            .map(|(_, charge)| charge)
            .collect::<Vec<_>>();
        assert_eq!(charges, vec![1; charged], "{name}");

        let report = interpreted.report();
        let count = |wanted: fn(&StructureInterpretIssue) -> bool| {
            report.issues().iter().filter(|issue| wanted(issue)).count()
        };
        let ssbond = (1..=3)
            .map(|serial| {
                bonds_with_evidence(
                    &interpreted,
                    &connection(Some(&serial.to_string()), "SSBOND"),
                )
            })
            .sum::<usize>();
        assert_eq!(
            bonds_with_evidence(&interpreted, &connection(None, "CONECT")),
            listed,
            "{name}"
        );
        assert_eq!(ssbond, disulfides, "{name}");
        assert_eq!(report.applied_connections(), listed + disulfides, "{name}");
        assert_eq!(
            count(|issue| matches!(
                issue,
                StructureInterpretIssue::ConnectionUnresolved {
                    reason: ConnectionResolutionReason::AlternateLocationOmitted { .. },
                    ..
                }
            )),
            omitted,
            "{name}"
        );
        assert_eq!(
            count(|issue| matches!(issue, StructureInterpretIssue::ResidueBondsInferred { .. })),
            0,
            "{name}"
        );

        // The model writes as mmCIF, planned from its source provenance or
        // from canonical classification, and reads back whole.
        let options = kekule::mmcif::MmcifWriteOptions::default;
        for written in [
            kekule::mmcif::write([&interpreted], options()),
            kekule::mmcif::write([interpreted.model()], options()),
        ] {
            let reread = kekule::mmcif::parse_str(&written.unwrap(), Default::default())
                .unwrap()
                .interpret(Default::default())
                .unwrap();
            let topology = reread.topology();
            assert_eq!(
                (
                    topology.atom_count(),
                    topology.bond_count(),
                    topology.instance_count()
                ),
                (atoms, bonds, molecules),
                "{name}"
            );
        }
    }
}

fn bonds_with_evidence(result: &StructureInterpretation, evidence: &BondEvidence) -> usize {
    result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.bonds())
        .filter(|bond| bond.evidence() == evidence)
        .count()
}
