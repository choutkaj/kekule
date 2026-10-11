use kekule::core::BondOrder;
use kekule::diagnostics::{DiagnosticKind, SourceLocation};
use kekule::pdb::{
    self, AltLocPolicy, AltLocResidue, BondEvidence, ConnectionResolutionReason, EntityKind,
    ModelSelection, PdbAtomKind, PdbParseOptions, SourceResidueId, SourceResiduePosition,
    StructureInterpretIssue, StructureInterpretOptions, StructureInterpretation,
};

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
        format!(
            "{:<6}{:>5} {:<4}{}{:>3} {}{:>4}{}   {x:>8.3}{y:>8.3}{z:>8.3}{:>6.2}{:>6.2}      {:<4}{:>2}{:<2}",
            self.record,
            self.serial,
            self.name,
            self.alternate_location,
            self.residue,
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
    pdb::interpret(&pdb::parse_str(text).unwrap()).unwrap()
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
            let molecule = result.topology().molecule(a.molecule()).unwrap().molecule();
            let id = molecule.bond_between(a.atom(), b.atom()).unwrap().unwrap();
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
                .molecule(atom.atom().molecule())
                .unwrap()
                .molecule();
            let charge = molecule.atom(atom.atom().atom()).unwrap().formal_charge;
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
    SourceResidueId {
        model_id: "1".to_owned(),
        asym_id: asym.to_owned(),
        position,
        insertion_code: None,
    }
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
    ]
    .join("\r\n")
        + "\rEND\r\n";
    let document = pdb::parse_str(&text).unwrap();
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
    let document = pdb::parse_str(&text).unwrap();
    let serials = document
        .models()
        .iter()
        .map(|model| model.serial())
        .collect::<Vec<_>>();
    assert_eq!(serials, [1, 2]);
    let error = pdb::interpret(&document).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);

    let second = pdb::interpret_with_options(
        &document,
        StructureInterpretOptions {
            model_selection: ModelSelection::Select("2".to_owned()),
            ..StructureInterpretOptions::default()
        },
    )
    .unwrap();
    assert_eq!(second.report().selected_model(), Some("2"));
    assert_eq!(second.report().ignored_coordinate_models(), ["1"]);
    assert_eq!(second.topology().atom_count(), 3);
    assert_eq!(second.topology().bond_count(), 2);

    let ensemble = pdb::interpret_ensemble(&document).unwrap();
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
    let document = pdb::parse_str(&text).unwrap();
    let repeated = residue(
        "W",
        SourceResiduePosition::RepeatedAuthor {
            sequence: "1".to_owned(),
            repeat: 1,
        },
    );
    assert_eq!(
        document.alternate_locations().unwrap(),
        [AltLocResidue {
            residue: repeated.clone(),
            labels: vec!["A".to_owned(), "B".to_owned()],
        }]
    );

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
        let interpreted = pdb::interpret_with_options(
            &document,
            StructureInterpretOptions {
                altloc_policy: policy,
                ..StructureInterpretOptions::default()
            },
        )
        .unwrap();
        let [decision] = interpreted.report().alternate_locations() else {
            panic!("one alternate-location decision");
        };
        assert_eq!(decision.residue, repeated);
        assert_eq!(decision.selected_labels, [label]);
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
                .molecule(atom.atom().molecule())
                .unwrap()
                .molecule();
            molecule.atom(atom.atom().atom()).unwrap().element.symbol()
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
        let error = pdb::parse_str(&text).unwrap_err();
        assert_eq!((error.kind(), error.location()), (kind, location), "{text}");
    }

    let limited = pdb::parse_str_with_options(
        &format!("{water}\n{water}"),
        PdbParseOptions {
            max_atom_records: 1,
            ..PdbParseOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(
        (limited.kind(), limited.location()),
        (DiagnosticKind::Limit, SourceLocation::at_line(2))
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
            StructureInterpretOptions {
                strict_entity_metadata: true,
                ..StructureInterpretOptions::default()
            },
            DiagnosticKind::Inconsistent,
            None,
        ),
    ];
    for (text, options, kind, line) in cases {
        let document = pdb::parse_str(&text).unwrap();
        let error = pdb::interpret_with_options(&document, options).unwrap_err();
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
            .molecules()
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
