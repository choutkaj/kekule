//! How the PDB and mmCIF readers and writers handle the conventions other
//! programs follow: element columns left blank, CONECT records per model,
//! charges, VMD numbering, molecular-dynamics residue names, and metal bonds
//! within one residue.

use std::sync::Arc;

use kekule::components::ComponentLibrary;
use kekule::diagnostics::DiagnosticKind;
use kekule::diagnostics::{FormatErrorDetail, UnrepresentablePolicy, WriteOmission};
use kekule::macromolecular::{
    ConnectivityOptions, EntityKind, ModelSelection, NumberField, StructureInterpretIssue,
    StructureInterpretOptions, StructureInterpretation,
};
use kekule::mmcif::{MmcifWriteErrorKind, MmcifWriteOptions};
use kekule::molecule::BondOrder;
use kekule::pdb::{PdbWriteErrorKind, PdbWriteOptions};
use kekule::structure::{Positions, Structure};
use kekule::topology::ResidueClass;
use kekule::{mmcif, pdb, smiles};

/// One `ATOM` or `HETATM` record in the 80 fixed columns: `name` fills
/// columns 13-16 as written.
#[allow(clippy::too_many_arguments)]
fn record(
    kind: &str,
    serial: &str,
    name: &str,
    residue: &str,
    chain: char,
    sequence: &str,
    [x, y, z]: [f64; 3],
    element: &str,
    charge: &str,
) -> String {
    let residue = if residue.len() == 4 {
        residue.to_owned()
    } else {
        format!("{residue:>3} ")
    };
    format!(
        "{kind:<6}{serial:>5} {name:<4} {residue}{chain}{sequence:>4}    {x:>8.3}{y:>8.3}{z:>8.3}{:>6.2}{:>6.2}          {element:>2}{charge:<2}",
        1.0, 0.0
    )
    .trim_end()
    .to_owned()
}

fn entry(lines: &[String]) -> String {
    lines.join("\n") + "\nEND\n"
}

fn read(text: &str) -> StructureInterpretation {
    pdb::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

fn read_mmcif(text: &str) -> StructureInterpretation {
    mmcif::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// Element symbols in dense atom order.
fn elements(structure: &Structure) -> Vec<&'static str> {
    structure
        .topology()
        .atoms()
        .map(|atom| atom.atom().element().symbol())
        .collect()
}

/// Bonds by dense atom indices and order.
fn bonds(structure: &Structure) -> Vec<(usize, usize, BondOrder)> {
    let mut bonds = structure
        .topology()
        .bonds()
        .map(|bond| {
            let [left, right] = bond.atoms().map(|atom| atom.index().index());
            (left.min(right), left.max(right), bond.bond().order())
        })
        .collect::<Vec<_>>();
    bonds.sort_unstable_by_key(|&(left, right, _)| (left, right));
    bonds
}

/// Formal charges in dense atom order.
fn charges(structure: &Structure) -> Vec<i8> {
    structure
        .topology()
        .atoms()
        .map(|atom| atom.atom().formal_charge())
        .collect()
}

#[test]
fn blank_element_columns_take_template_elements() {
    // CHARMM ions, and zinc named in the one-letter columns. `CLA` is
    // chlorophyll in the dictionary, so the CHARMM library decides.
    let ions = [
        ("1", "SOD ", "SOD", [0.0, 0.0, 0.0]),
        ("2", "POT ", "POT", [5.0, 0.0, 0.0]),
        ("3", "CLA ", "CLA", [10.0, 0.0, 0.0]),
        ("4", "CAL ", "CAL", [15.0, 0.0, 0.0]),
        ("5", " ZN ", "ZN", [20.0, 0.0, 0.0]),
    ]
    .map(|(serial, name, residue, position)| {
        record(
            "HETATM", serial, name, residue, 'I', serial, position, "", "",
        )
    });
    let charmm = StructureInterpretOptions::default().with_connectivity(
        ConnectivityOptions::default().with_components(Some(Arc::new(ComponentLibrary::charmm()))),
    );
    let read = pdb::parse_str(&entry(&ions), Default::default())
        .unwrap()
        .interpret(charmm)
        .unwrap();
    assert_eq!(elements(read.structure()), ["Na", "K", "Cl", "Ca", "Zn"]);

    // Left-aligned names of residues the dictionary defines.
    let residues = [
        ("1", "N   ", "ILE", "1", [0.000, 0.000, 0.000]),
        ("2", "CA  ", "ILE", "1", [1.458, 0.000, 0.000]),
        ("3", "CD1 ", "ILE", "1", [3.000, 2.000, 0.000]),
        ("4", "HG  ", "SER", "2", [10.000, 0.000, 0.000]),
        ("5", "CE  ", "MET", "3", [20.000, 0.000, 0.000]),
        ("6", "NE2 ", "HIS", "4", [30.000, 0.000, 0.000]),
    ]
    .map(|(serial, name, residue, sequence, position)| {
        record(
            "HETATM", serial, name, residue, 'A', sequence, position, "", "",
        )
    });
    let read = self::read(&entry(&residues));
    assert_eq!(elements(read.structure()), ["N", "C", "C", "H", "C", "N"]);

    // The heme's pyrrole nitrogens are not niobium and neodymium.
    let heme = [
        ("1", "FE  ", [0.0, 0.0, 0.0]),
        ("2", "NB  ", [2.0, 0.0, 0.0]),
        ("3", "ND  ", [-2.0, 0.0, 0.0]),
    ]
    .map(|(serial, name, position)| {
        record("HETATM", serial, name, "HEM", 'A', "1", position, "", "")
    });
    let read = self::read(&entry(&heme));
    assert_eq!(elements(read.structure()), ["Fe", "N", "N"]);
}

#[test]
fn ambiguous_atom_names_without_templates_fail() {
    // `C1` starts in column 13, so the writer left-aligns names, and `CL1`
    // may be carbon or chlorine.
    let text = entry(&[
        record(
            "HETATM",
            "1",
            "C1  ",
            "LIG",
            'A',
            "1",
            [0.0, 0.0, 0.0],
            "",
            "",
        ),
        record(
            "HETATM",
            "2",
            "CL1 ",
            "LIG",
            'A',
            "1",
            [1.75, 0.0, 0.0],
            "",
            "",
        ),
    ]);
    let error = pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::InvalidValue);
    assert!(error.to_string().contains("\"CL1\""), "{error}");

    // Names aligned by the convention read two-letter symbols from
    // column 13.
    let text = entry(&[
        record(
            "HETATM",
            "1",
            " C1 ",
            "LIG",
            'A',
            "1",
            [0.0, 0.0, 0.0],
            "",
            "",
        ),
        record(
            "HETATM",
            "2",
            "CL1 ",
            "LIG",
            'A',
            "1",
            [1.75, 0.0, 0.0],
            "",
            "",
        ),
    ]);
    assert_eq!(elements(read(&text).structure()), ["C", "Cl"]);
}

/// Cisplatin, a component the dictionary bundled with Kekule lacks.
const CISPLATIN: [(&str, &str, [f64; 3], &str); 5] = [
    ("1", "PT1 ", [0.000, 0.000, 0.000], "PT"),
    ("2", " N1 ", [2.050, 0.000, 0.000], "N"),
    ("3", " N2 ", [0.000, 2.050, 0.000], "N"),
    ("4", "CL1 ", [-2.300, 0.000, 0.000], "CL"),
    ("5", "CL2 ", [0.000, -2.300, 0.000], "CL"),
];

#[test]
fn metal_bonds_of_a_residue_without_hydrogens_are_covalent_and_reported() {
    let mut lines = CISPLATIN
        .map(|(serial, name, position, element)| {
            record(
                "HETATM", serial, name, "CPT", 'A', "1", position, element, "",
            )
        })
        .to_vec();
    lines.push("CONECT    1    2    3    4    5".to_owned());
    let read = read(&entry(&lines));
    let structure = read.structure();
    assert_eq!(
        bonds(structure),
        [
            (0, 1, BondOrder::Single),
            (0, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (0, 4, BondOrder::Single),
        ]
    );
    // The chlorides bond platinum, so they carry no hydrogen.
    let mut perceived = structure.clone();
    perceived.perceive().unwrap();
    let hydrogens = perceived
        .topology()
        .atoms()
        .map(|atom| atom.implicit_hydrogens())
        .collect::<Vec<_>>();
    assert_eq!(hydrogens[3..], [Some(0), Some(0)]);
    // Without hydrogens the file cannot tell dative from covalent bonds.
    let unstated = read
        .report()
        .issues()
        .iter()
        .filter(|issue| {
            matches!(
                issue,
                StructureInterpretIssue::MetalBondsUnstated { bonds: 4, .. }
            )
        })
        .count();
    assert_eq!(unstated, 1, "{:?}", read.report().issues());

    // An mmCIF component table without coordination flags agrees.
    let cif = "data_x
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
CPT PT1 N1 SING
CPT PT1 N2 SING
CPT PT1 CL1 SING
CPT PT1 CL2 SING
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
_atom_site.pdbx_PDB_model_num
HETATM 1 PT PT1 CPT A 1 . 0.000 0.000 0.000 1
HETATM 2 N N1 CPT A 1 . 2.050 0.000 0.000 1
HETATM 3 N N2 CPT A 1 . 0.000 2.050 0.000 1
HETATM 4 CL CL1 CPT A 1 . -2.300 0.000 0.000 1
HETATM 5 CL CL2 CPT A 1 . 0.000 -2.300 0.000 1
";
    let read = read_mmcif(cif);
    assert_eq!(bonds(read.structure()), bonds(structure));
}

/// Cisplatin with its ammine hydrogens, without CONECT records: each
/// hydrogen attaches to its nearest nitrogen.
fn cisplatin_with_hydrogens() -> Vec<String> {
    let hydrogens = [
        ("6", "H11 ", [2.390, 0.940, 0.000], "H"),
        ("7", "H12 ", [2.390, -0.470, 0.820], "H"),
        ("8", "H13 ", [2.390, -0.470, -0.820], "H"),
        ("9", "H21 ", [0.940, 2.390, 0.000], "H"),
        ("10", "H22 ", [-0.470, 2.390, 0.820], "H"),
        ("11", "H23 ", [-0.470, 2.390, -0.820], "H"),
    ];
    CISPLATIN
        .iter()
        .chain(&hydrogens)
        .map(|&(serial, name, position, element)| {
            record(
                "HETATM", serial, name, "CPT", 'A', "1", position, element, "",
            )
        })
        .collect()
}

/// The bonds of the first atom, by dense atom indices and order.
fn first_atom_bonds(structure: &Structure) -> Vec<(usize, usize, BondOrder)> {
    bonds(structure)
        .into_iter()
        .filter(|&(left, _, _)| left == 0)
        .collect()
}

fn unstated_metal_bonds(read: &StructureInterpretation) -> Vec<usize> {
    read.report()
        .issues()
        .iter()
        .filter_map(|issue| match issue {
            StructureInterpretIssue::MetalBondsUnstated { bonds, .. } => Some(*bonds),
            _ => None,
        })
        .collect()
}

#[test]
fn explicit_hydrogens_tell_dative_metal_bonds_from_covalent_ones() {
    // CONECT lists only the platinum bonds.
    let mut lines = cisplatin_with_hydrogens();
    lines.push("CONECT    1    2    3    4    5".to_owned());
    let read = read(&entry(&lines));
    let structure = read.structure();
    let metal_bonds = bonds(structure)
        .into_iter()
        .filter(|&(left, _, _)| left == 0)
        .collect::<Vec<_>>();
    assert_eq!(
        metal_bonds,
        [
            (0, 1, BondOrder::Dative),
            (0, 2, BondOrder::Dative),
            (0, 3, BondOrder::Single),
            (0, 4, BondOrder::Single),
        ]
    );
    assert!(charges(structure).iter().all(|&charge| charge == 0));
    let mut perceived = structure.clone();
    perceived.perceive().unwrap();
    let implicit = perceived
        .topology()
        .atoms()
        .take(5)
        .map(|atom| atom.implicit_hydrogens())
        .collect::<Vec<_>>();
    assert_eq!(implicit, [Some(0); 5]);
    assert!(!read
        .report()
        .issues()
        .iter()
        .any(|issue| matches!(issue, StructureInterpretIssue::MetalBondsUnstated { .. })));
}

/// Ethene as residue `ETH`, in two models, each with its own CONECT
/// records listing the double bond.
fn ethene_models(conect: &[&str]) -> String {
    let mut lines = Vec::new();
    for model in 1..=2 {
        lines.push(format!("MODEL     {model:>4}"));
        lines.push(record(
            "HETATM",
            "1",
            " C1 ",
            "ETH",
            'A',
            "1",
            [0.0, 0.0, 0.0],
            "C",
            "",
        ));
        lines.push(record(
            "HETATM",
            "2",
            " C2 ",
            "ETH",
            'A',
            "1",
            [1.33, 0.0, 0.0],
            "C",
            "",
        ));
        lines.extend(conect.iter().map(|line| (*line).to_owned()));
        lines.push("ENDMDL".to_owned());
    }
    entry(&lines)
}

#[test]
fn conect_multiplicity_counts_within_one_model() {
    let first = |text: &str| {
        pdb::parse_str(text, Default::default())
            .unwrap()
            .interpret(
                StructureInterpretOptions::default().with_model_selection(ModelSelection::First),
            )
            .unwrap()
    };
    let double = ethene_models(&["CONECT    1    2    2", "CONECT    2    1    1"]);
    assert_eq!(
        bonds(first(&double).structure()),
        [(0, 1, BondOrder::Double)]
    );
    let single = ethene_models(&["CONECT    1    2", "CONECT    2    1"]);
    assert_eq!(
        bonds(first(&single).structure()),
        [(0, 1, BondOrder::Single)]
    );
}

#[test]
fn charge_columns_read_explicit_zeros_and_unsigned_digits() {
    let charged = |charge: &str| {
        let text = entry(&[record(
            "HETATM",
            "1",
            " N  ",
            "LIG",
            'A',
            "1",
            [0.0, 0.0, 0.0],
            "N",
            charge,
        )]);
        pdb::parse_str(&text, Default::default()).unwrap().atoms()[0].formal_charge()
    };
    for zero in ["0", " 0", "0+", "0-"] {
        assert_eq!(charged(zero), Some(0), "{zero:?}");
    }
    assert_eq!(charged("1"), Some(1));
    assert_eq!(charged("1+"), Some(1));
    assert_eq!(charged(""), None);
}

/// Lysine without hydrogens, every charge stated 0 in mmCIF.
const STATED_ZERO_LYSINE: &str = "data_x
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
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
HETATM 1 N N LYS A 1 . 0.000 0.000 0.000 0 1 A 1
HETATM 2 C CA LYS A 1 . 1.458 0.000 0.000 0 1 A 1
HETATM 3 C C LYS A 1 . 2.009 1.420 0.000 0 1 A 1
HETATM 4 O O LYS A 1 . 1.251 2.390 0.000 0 1 A 1
HETATM 5 C CB LYS A 1 . 1.990 -0.770 -1.210 0 1 A 1
HETATM 6 C CG LYS A 1 . 3.510 -0.770 -1.210 0 1 A 1
HETATM 7 C CD LYS A 1 . 4.040 -1.540 -2.420 0 1 A 1
HETATM 8 C CE LYS A 1 . 5.560 -1.540 -2.420 0 1 A 1
HETATM 9 N NZ LYS A 1 . 6.070 -2.270 -3.590 0 1 A 1
";

#[test]
fn charges_of_templated_residues_round_trip_through_pdb() {
    // A stated 0 in a residue without hydrogens is no evidence: NZ takes
    // its template's +1.
    let structure = read_mmcif(STATED_ZERO_LYSINE).into_structure();
    assert_eq!(charges(&structure)[8], 1);
    let (written, _) = pdb::write(&structure, Default::default()).unwrap();
    let reread = read(&written).into_structure();
    assert_eq!(charges(&reread), charges(&structure));
}

#[test]
fn charges_readers_derive_differently_are_not_written_silently() {
    // A neutral lysine NZ without hydrogens: readers ignore its stated 0
    // and give it the template's +1.
    let source = read_mmcif(STATED_ZERO_LYSINE).into_structure();
    let nz = source.topology().atom_ids()[8];
    let mut editor = source.into_editor();
    let handle = editor.topology().atom_handle(nz).unwrap();
    let neutral = editor
        .topology()
        .atom(handle)
        .unwrap()
        .clone()
        .with_formal_charge(0);
    editor.replace_atom(handle, neutral).unwrap();
    let structure = editor.finish().unwrap();
    assert_eq!(charges(&structure)[8], 0);

    let error = pdb::write(&structure, Default::default()).unwrap_err();
    let Some(FormatErrorDetail::PdbWrite(kind)) = error.detail() else {
        panic!("{error:?}");
    };
    assert_eq!(
        *kind,
        PdbWriteErrorKind::ChargeNotReadBack {
            atom: nz,
            charge: 0,
            read: 1
        }
    );
    let options = PdbWriteOptions::default().with_chemistry(UnrepresentablePolicy::Omit);
    let (_, report) = pdb::write(&structure, options).unwrap();
    assert!(
        report.omissions().contains(&WriteOmission::FormalCharges {
            atoms: 1,
            first: Some(nz)
        }),
        "{report:?}"
    );

    let error = mmcif::write([&structure], Default::default()).unwrap_err();
    let Some(FormatErrorDetail::MmcifWrite(kind)) = error.detail() else {
        panic!("{error:?}");
    };
    assert_eq!(
        *kind,
        MmcifWriteErrorKind::ChargeNotReadBack {
            atom: nz,
            charge: 0,
            read: 1
        }
    );
    let options = MmcifWriteOptions::default().with_chemistry(UnrepresentablePolicy::Omit);
    let (_, report) = mmcif::write([&structure], options).unwrap();
    assert!(
        report.omissions().contains(&WriteOmission::FormalCharges {
            atoms: 1,
            first: Some(nz)
        }),
        "{report:?}"
    );
}

#[test]
fn vmd_hexadecimal_numbers_continue_the_decimal_range() {
    let water = |serial: &str, sequence: &str, x: f64| {
        record(
            "HETATM",
            serial,
            " O  ",
            "HOH",
            'W',
            sequence,
            [x, 0.0, 0.0],
            "O",
            "",
        )
    };
    let text = entry(&[
        water("99998", "9998", 0.0),
        water("99999", "9999", 3.0),
        water("186a0", "2710", 6.0),
        water("186a1", "2711", 9.0),
    ]);
    let document = pdb::parse_str(&text, Default::default()).unwrap();
    let serials = document
        .atoms()
        .iter()
        .map(|atom| atom.serial())
        .collect::<Vec<_>>();
    assert_eq!(
        serials,
        [Some(99_998), Some(99_999), Some(100_000), Some(100_001)]
    );
    let numbers = document
        .atoms()
        .iter()
        .map(|atom| atom.residue_sequence())
        .collect::<Vec<_>>();
    assert_eq!(numbers, [9998, 9999, 10_000, 10_001]);
}

#[test]
fn molecular_dynamics_water_and_ions_are_never_polymer() {
    // GROMACS names water SOL and writes every record as ATOM; the chain
    // ends with TER.
    let mut lines = Vec::new();
    for (index, offset) in [0.0, 4.0, 8.0].into_iter().enumerate() {
        let sequence = (index + 1).to_string();
        for (name, [x, y, z], element) in [
            (" OW ", [0.0, 0.0, 0.0], "O"),
            (" HW1", [0.96, 0.0, 0.0], "H"),
            (" HW2", [-0.24, 0.93, 0.0], "H"),
        ] {
            lines.push(record(
                "ATOM",
                &(lines.len() + 1).to_string(),
                name,
                "SOL",
                'W',
                &sequence,
                [x + offset, y, z],
                element,
                "",
            ));
        }
    }
    lines.push(format!("TER   {:>5}      SOL W   3", lines.len() + 1));
    for (index, (name, residue, element)) in [(" NA ", "NA", "NA"), (" CL ", "CL", "CL")]
        .into_iter()
        .enumerate()
    {
        lines.push(record(
            "ATOM",
            &(lines.len() + 1).to_string(),
            name,
            residue,
            'I',
            &(index + 1).to_string(),
            [20.0 + 5.0 * index as f64, 0.0, 0.0],
            element,
            "",
        ));
    }
    lines.push(format!("TER   {:>5}       CL I   2", lines.len() + 1));
    let read = read(&entry(&lines));
    let kinds = read
        .report()
        .instances()
        .iter()
        .map(|instance| instance.entity_kinds().to_vec())
        .collect::<Vec<_>>();
    assert_eq!(
        kinds,
        [
            vec![EntityKind::Water],
            vec![EntityKind::Water],
            vec![EntityKind::Water],
            vec![EntityKind::NonPolymer],
            vec![EntityKind::NonPolymer],
        ]
    );
    // Classification agrees with reading.
    let classes = read
        .structure()
        .topology()
        .residues()
        .map(|residue| residue.class())
        .collect::<Vec<_>>();
    assert_eq!(
        classes,
        [
            ResidueClass::Water,
            ResidueClass::Water,
            ResidueClass::Water,
            ResidueClass::Ion,
            ResidueClass::Ion,
        ]
    );
}

#[test]
fn mmcif_charges_with_signs_and_coordinates_with_uncertainties_parse() {
    let cif = "data_x
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
_atom_site.pdbx_formal_charge
_atom_site.pdbx_PDB_model_num
HETATM 1 NA NA NA A 1 . 12.345(3) 0.000 0.000 1+ 1
HETATM 2 CL CL CL B 2 . 15.000 0.000(12) 0.000 1- 1
";
    let read = read_mmcif(cif);
    let structure = read.structure();
    assert_eq!(charges(structure), [1, -1]);
    let x = structure
        .positions()
        .get(0)
        .unwrap()
        .value_in(kekule::units::ANGSTROM)
        .unwrap()
        .x;
    assert!((x - 12.345).abs() < 1e-9);
}

/// Glycine residues in chain `chain`, numbered from `first`, the first
/// residue's atoms shifted by `shift`, and the next residue's `N` at
/// `gap` from its `C`.
fn glycines(chain: char, first: usize, serial: usize, shift: [f64; 2], gap: f64) -> Vec<String> {
    let [dx, dy] = shift;
    let second = [gap - 1.328, 0.0];
    [
        (" N  ", 0, [0.000, 0.000], "N"),
        (" CA ", 0, [1.458, 0.000], "C"),
        (" C  ", 0, [2.009, 1.420], "C"),
        (" O  ", 0, [1.251, 2.390], "O"),
        (" N  ", 1, [3.332, 1.536], "N"),
        (" CA ", 1, [3.988, 2.836], "C"),
        (" C  ", 1, [5.488, 2.736], "C"),
        (" O  ", 1, [6.100, 1.680], "O"),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, (name, residue, [x, y], element))| {
        let x = x + dx + if residue == 1 { second[0] } else { 0.0 };
        record(
            "ATOM",
            &(serial + index).to_string(),
            name,
            "GLY",
            chain,
            &(first + residue).to_string(),
            [x, y + dy, 0.0],
            element,
            "",
        )
    })
    .collect()
}

#[test]
fn stretched_polymer_links_are_conect_records() {
    // The peptide bond between the glycines is 3.5 Å long, beyond bonding
    // distance, and a CONECT record states it.
    let mut lines = glycines('A', 1, 1, [0.0, 0.0], 3.5);
    lines.push("CONECT    3    5".to_owned());
    lines.push("CONECT    5    3".to_owned());
    let structure = read(&entry(&lines)).into_structure();
    assert!(bonds(&structure).contains(&(2, 4, BondOrder::Single)));
    let (written, _) = pdb::write(&structure, Default::default()).unwrap();
    let reread = read(&written).into_structure();
    assert_eq!(bonds(&reread), bonds(&structure));
}

/// The atom-site rows of two chains of glycines, both with author chain
/// `A`; the second chain starts at bonding distance from the first.
fn chains_sharing_a_letter() -> String {
    let mut rows = String::new();
    let atoms = [
        ("N", "N", 1, [0.000, 0.000]),
        ("C", "CA", 1, [1.458, 0.000]),
        ("C", "C", 1, [2.009, 1.420]),
        ("O", "O", 1, [1.251, 2.390]),
        ("N", "N", 2, [3.332, 1.536]),
        ("C", "CA", 2, [3.988, 2.836]),
        ("C", "C", 2, [5.488, 2.736]),
        ("O", "O", 2, [6.100, 1.680]),
    ];
    let mut serial = 1;
    for (asym, offset, [dx, dy]) in [("A", 0, [0.0, 0.0]), ("B", 2, [6.818, 2.736])] {
        for (element, name, sequence, [x, y]) in atoms {
            rows.push_str(&format!(
                "ATOM {serial} {element} {name} GLY {asym} 1 {sequence} {:.3} {:.3} 0.000 {} A 1\n",
                x + dx,
                y + dy,
                sequence + offset,
            ));
            serial += 1;
        }
    }
    format!(
        "data_x
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
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
{rows}"
    )
}

#[test]
fn chains_sharing_a_letter_stay_apart() {
    let structure = read_mmcif(&chains_sharing_a_letter()).into_structure();
    assert_eq!(structure.topology().chains().count(), 2);
    let (written, _) = pdb::write(&structure, Default::default()).unwrap();
    assert_eq!(
        written
            .lines()
            .filter(|line| line.starts_with("TER"))
            .count(),
        2
    );
    let reread = read(&written).into_structure();
    assert_eq!(bonds(&reread), bonds(&structure));
}

#[test]
fn pdb_identifiers_must_be_ascii() {
    let cif = |chain: &str, component: &str, code: &str| {
        format!(
            "data_x
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
HETATM 1 O O {component} X 1 . {code} 0.000 0.000 0.000 1 {chain} 1
"
        )
    };
    for text in [
        cif("\u{e9}", "HOH", "."),
        cif("A", "H\u{d6}H", "."),
        cif("A", "HOH", "\u{e9}"),
    ] {
        let structure = read_mmcif(&text).into_structure();
        let error = pdb::write(&structure, Default::default()).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Unsupported, "{text}");
    }
}

#[test]
fn mmcif_states_ligand_bonds_in_component_tables() {
    let structure = read(&entry(&[
        record(
            "HETATM",
            "1",
            " C1 ",
            "ETH",
            'A',
            "1",
            [0.0, 0.0, 0.0],
            "C",
            "",
        ),
        record(
            "HETATM",
            "2",
            " C2 ",
            "ETH",
            'A',
            "1",
            [1.33, 0.0, 0.0],
            "C",
            "",
        ),
        "CONECT    1    2    2".to_owned(),
        "CONECT    2    1    1".to_owned(),
    ]))
    .into_structure();
    let (written, _) = mmcif::write([&structure], Default::default()).unwrap();
    assert!(written.contains("_chem_comp_bond.value_order"), "{written}");
    assert!(written.contains("ETH C1 C2 doub"), "{written}");
    assert!(!written.contains("_struct_conn."), "{written}");
    let reread = read_mmcif(&written).into_structure();
    assert_eq!(bonds(&reread), bonds(&structure));
}

#[test]
fn reports_name_how_atoms_and_numbers_were_read() {
    let text = entry(&[
        record(
            "ATOM",
            "1",
            " CA ",
            "GLY",
            'A',
            "1",
            [0.0, 0.0, 0.0],
            "",
            "",
        ),
        record(
            "HETATM",
            "2",
            " C1 ",
            "LIG",
            'B',
            "2",
            [10.0, 0.0, 0.0],
            "",
            "",
        ),
        record(
            "HETATM",
            "3",
            " N  ",
            "NH4",
            'C',
            "3",
            [20.0, 0.0, 0.0],
            "N",
            "1",
        ),
    ]);
    let issues = read(&text).report().issues().to_vec();
    for issue in [
        StructureInterpretIssue::ElementsFromTemplates {
            source_lines: vec![1],
        },
        StructureInterpretIssue::ElementsInferredFromAtomNames {
            source_lines: vec![2],
        },
        StructureInterpretIssue::UnsignedCharges {
            source_lines: vec![3],
        },
    ] {
        assert!(issues.contains(&issue), "{issue:?} in {issues:?}");
    }

    let water = |serial: &str, sequence: &str, x: f64| {
        record(
            "HETATM",
            serial,
            " O  ",
            "HOH",
            'W',
            sequence,
            [x, 0.0, 0.0],
            "O",
            "",
        )
    };
    let hexadecimal = entry(&[water("99999", "9999", 0.0), water("186a0", "2710", 3.0)]);
    let issues = read(&hexadecimal).report().issues().to_vec();
    for field in [NumberField::AtomSerial, NumberField::ResidueNumber] {
        let issue = StructureInterpretIssue::HexadecimalNumbers {
            field,
            source_line: 2,
        };
        assert!(issues.contains(&issue), "{issue:?} in {issues:?}");
    }

    // Numbers that wrap around to zero, as some writers number past the
    // columns.
    let wrapped = entry(&[water("99999", "9999", 0.0), water("0", "0", 3.0)]);
    let issues = read(&wrapped).report().issues().to_vec();
    for field in [NumberField::AtomSerial, NumberField::ResidueNumber] {
        let issue = StructureInterpretIssue::NumbersRestart {
            field,
            source_line: 2,
        };
        assert!(issues.contains(&issue), "{issue:?} in {issues:?}");
    }
}

#[test]
fn donors_take_the_bond_their_valence_allows() {
    // A phosphonate phosphorus needs the platinum bond for valence five; a
    // phosphine has valence three without it.
    let atom = |serial: &str, name: &str, position: [f64; 3], element: &str| {
        record(
            "HETATM", serial, name, "PPT", 'A', "1", position, element, "",
        )
    };
    let lines = [
        atom("1", "PT1 ", [0.0, 0.0, 0.0], "PT"),
        atom("2", " P1 ", [2.3, 0.0, 0.0], "P"),
        atom("3", " O1 ", [3.0, 1.3, 0.0], "O"),
        atom("4", " O2 ", [3.0, -0.7, 1.2], "O"),
        atom("5", " O3 ", [3.0, -0.7, -1.2], "O"),
        atom("6", " HO2", [3.9, -0.5, 1.3], "H"),
        atom("7", " HO3", [3.9, -0.5, -1.3], "H"),
        atom("8", " P2 ", [0.0, 2.3, 0.0], "P"),
        atom("9", " HP1", [0.0, 3.7, 0.0], "H"),
        atom("10", " HP2", [1.2, 2.8, 0.0], "H"),
        atom("11", " HP3", [-1.2, 2.8, 0.0], "H"),
        atom("12", "CL1 ", [-2.3, 0.0, 0.0], "CL"),
        "CONECT    1    2    8   12".to_owned(),
        "CONECT    2    3    3    4    5".to_owned(),
        "CONECT    4    6".to_owned(),
        "CONECT    5    7".to_owned(),
        "CONECT    8    9   10   11".to_owned(),
    ];
    let read = read(&entry(&lines));
    assert_eq!(
        first_atom_bonds(read.structure()),
        [
            (0, 1, BondOrder::Single),
            (0, 7, BondOrder::Dative),
            (0, 11, BondOrder::Single),
        ]
    );
    assert!(charges(read.structure()).iter().all(|&charge| charge == 0));
    assert!(unstated_metal_bonds(&read).is_empty());
}

#[test]
fn carbene_bonds_stay_covalent_and_are_reported() {
    // A diaminocarbene carbon is short of its valence either way.
    let atom = |serial: &str, name: &str, position: [f64; 3], element: &str| {
        record(
            "HETATM", serial, name, "CBN", 'A', "1", position, element, "",
        )
    };
    let lines = [
        atom("1", "PT1 ", [0.0, 0.0, 0.0], "PT"),
        atom("2", " C1 ", [2.0, 0.0, 0.0], "C"),
        atom("3", " N1 ", [2.7, 1.1, 0.0], "N"),
        atom("4", " N2 ", [2.7, -1.1, 0.0], "N"),
        atom("5", " H11", [3.7, 1.1, 0.0], "H"),
        atom("6", " H12", [2.2, 2.0, 0.0], "H"),
        atom("7", " H21", [3.7, -1.1, 0.0], "H"),
        atom("8", " H22", [2.2, -2.0, 0.0], "H"),
        atom("9", "CL1 ", [-2.3, 0.0, 0.0], "CL"),
        "CONECT    1    2    9".to_owned(),
        "CONECT    2    3    4".to_owned(),
        "CONECT    3    5    6".to_owned(),
        "CONECT    4    7    8".to_owned(),
    ];
    let read = read(&entry(&lines));
    assert_eq!(
        first_atom_bonds(read.structure()),
        [(0, 1, BondOrder::Single), (0, 8, BondOrder::Single)]
    );
    assert_eq!(unstated_metal_bonds(&read), [1]);
}

#[test]
fn link_records_within_one_residue_are_decided_like_conect() {
    // The first ammine's platinum bond is a LINK record.
    let mut lines = cisplatin_with_hydrogens();
    lines.insert(
        0,
        format!("LINK        PT1  CPT A   1{:16} N1  CPT A   1", ""),
    );
    // A bond the residue states makes its listing complete, so the
    // hydrogens are listed too.
    lines.extend([
        "CONECT    1    3    4    5".to_owned(),
        "CONECT    2    6    7    8".to_owned(),
        "CONECT    3    9   10   11".to_owned(),
    ]);
    let read = read(&entry(&lines));
    assert_eq!(
        first_atom_bonds(read.structure()),
        [
            (0, 1, BondOrder::Dative),
            (0, 2, BondOrder::Dative),
            (0, 3, BondOrder::Single),
            (0, 4, BondOrder::Single),
        ]
    );
    assert!(charges(read.structure()).iter().all(|&charge| charge == 0));
}

#[test]
fn pdb_writes_metal_bonds_a_reader_gives_back() {
    let mut lines = cisplatin_with_hydrogens();
    lines.push("CONECT    1    2    3    4    5".to_owned());
    let structure = read(&entry(&lines)).into_structure();
    let (written, _) = pdb::write(&structure, Default::default()).unwrap();
    assert_eq!(bonds(read(&written).structure()), bonds(&structure));

    // A covalent bond to an ammine nitrogen would read back as dative.
    let topology = smiles::to_topology("[H][N]([H])([H])[Pt]Cl").unwrap();
    let structure = Structure::new(topology, Positions::zeros(6)).unwrap();
    let hydrogens = PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let error = pdb::write(&structure, hydrogens.clone()).unwrap_err();
    assert!(
        matches!(
            error.detail(),
            Some(FormatErrorDetail::PdbWrite(
                PdbWriteErrorKind::UnsupportedBondOrder {
                    order: BondOrder::Single,
                    ..
                }
            ))
        ),
        "{error}"
    );
    let (written, report) = pdb::write(
        &structure,
        hydrogens.with_chemistry(UnrepresentablePolicy::Omit),
    )
    .unwrap();
    assert!(report
        .omissions()
        .contains(&WriteOmission::Bonds { bonds: 1 }));
    // The nitrogen lists its hydrogens only.
    assert!(
        written.contains("CONECT    2    1    3    4\n"),
        "{written}"
    );
}

#[test]
fn aromatic_rings_on_metals_take_the_bonds_that_localize() {
    // Cyclopentadienyl iron chloride with the ring given as aromatic bonds:
    // only covalent iron-carbon bonds let the ring localize.
    let mut rows = String::new();
    let mut atoms = vec![
        ("FE1", "FE", [0.0, 0.0, 0.0]),
        ("CL1", "CL", [0.0, 0.0, -2.2]),
    ];
    let mut table = String::from("CPF FE1 CL1 SING\n");
    let ring = ["C1", "C2", "C3", "C4", "C5"];
    let hydrogens = ["H1", "H2", "H3", "H4", "H5"];
    for (index, (carbon, hydrogen)) in ring.iter().zip(hydrogens).enumerate() {
        let angle = std::f64::consts::TAU * index as f64 / 5.0;
        let (sin, cos) = angle.sin_cos();
        atoms.push((carbon, "C", [1.2 * cos, 1.2 * sin, 1.65]));
        atoms.push((hydrogen, "H", [2.2 * cos, 2.2 * sin, 1.65]));
        let next = ring[(index + 1) % 5];
        table.push_str(&format!(
            "CPF FE1 {carbon} SING\nCPF {carbon} {next} AROM\nCPF {carbon} {hydrogen} SING\n"
        ));
    }
    for (index, (name, element, [x, y, z])) in atoms.iter().enumerate() {
        rows.push_str(&format!(
            "HETATM {} {element} {name} CPF A 1 . {x:.3} {y:.3} {z:.3} 1\n",
            index + 1
        ));
    }
    let cif = format!(
        "data_x
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
{table}loop_
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
_atom_site.pdbx_PDB_model_num
{rows}"
    );
    let read = read_mmcif(&cif);
    let iron = first_atom_bonds(read.structure());
    assert_eq!(iron.len(), 6);
    assert!(iron.iter().all(|&(_, _, order)| order == BondOrder::Single));
    assert!(unstated_metal_bonds(&read).is_empty());
}

#[test]
fn stated_covalent_bonds_outrank_unflagged_table_rows() {
    let cif = "data_x
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
PTA PT1 N1 SING
PTA PT1 CL1 SING
PTA N1 H11 SING
PTA N1 H12 SING
PTA N1 H13 SING
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
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_model_num
HETATM 1 PT PT1 PTA A 1 . 0.000 0.000 0.000 1 1
HETATM 2 N N1 PTA A 1 . 2.050 0.000 0.000 1 1
HETATM 3 CL CL1 PTA A 1 . -2.300 0.000 0.000 1 1
HETATM 4 H H11 PTA A 1 . 2.390 0.940 0.000 1 1
HETATM 5 H H12 PTA A 1 . 2.390 -0.470 0.820 1 1
HETATM 6 H H13 PTA A 1 . 2.390 -0.470 -0.820 1 1
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_auth_seq_id
covale1 covale A PTA PT1 1 A PTA N1 1
";
    let read = read_mmcif(cif);
    assert_eq!(
        first_atom_bonds(read.structure()),
        [(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)]
    );
}

#[test]
fn padded_hydrogen_names_read_as_hydrogen() {
    let atom = |serial: &str, name: &str, x: f64| {
        record(
            "HETATM",
            serial,
            name,
            "LIG",
            'A',
            "1",
            [x, 0.0, 0.0],
            "",
            "",
        )
    };
    let read = read(&entry(&[
        atom("1", " C1 ", 0.0),
        atom("2", " O1 ", 1.43),
        atom("3", " HO1", 2.4),
        "CONECT    1    2".to_owned(),
        "CONECT    2    3".to_owned(),
    ]));
    assert_eq!(elements(read.structure()), ["C", "O", "H"]);
    assert!(read.report().issues().iter().any(|issue| matches!(
        issue,
        StructureInterpretIssue::ElementsInferredFromAtomNames { .. }
    )));

    // An atom with its element column shows nothing about alignment.
    let text = entry(&[
        record(
            "HETATM",
            "1",
            "ZN  ",
            " ZN",
            'B',
            "2",
            [9.0, 0.0, 0.0],
            "ZN",
            "",
        ),
        atom("2", " C1 ", 0.0),
        atom("3", " CL1", 1.75),
    ]);
    let error = pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap_err();
    assert!(error.to_string().contains("\"CL1\""), "{error}");
}

#[test]
fn each_model_takes_its_own_conect_records() {
    let mut text = ethene_models(&["CONECT    1    2    2", "CONECT    2    1    1"]);
    // The second model lists the bond once.
    for (twice, once) in [
        ("CONECT    1    2    2", "CONECT    1    2"),
        ("CONECT    2    1    1", "CONECT    2    1"),
    ] {
        let at = text.rfind(twice).unwrap();
        text.replace_range(at..at + twice.len(), once);
    }
    let model = |selection| {
        pdb::parse_str(&text, Default::default())
            .unwrap()
            .interpret(StructureInterpretOptions::default().with_model_selection(selection))
            .unwrap()
            .into_structure()
    };
    assert_eq!(
        bonds(&model(ModelSelection::First)),
        [(0, 1, BondOrder::Double)]
    );
    assert_eq!(
        bonds(&model(ModelSelection::Select("2".to_owned()))),
        [(0, 1, BondOrder::Single)]
    );
}

#[test]
fn padded_names_without_element_columns_are_ambiguous() {
    // GROMACS starts every short name in column 14, so ` CL1` may be carbon
    // or chlorine.
    let text = entry(&[
        record(
            "HETATM",
            "1",
            " C1 ",
            "LIG",
            'A',
            "1",
            [0.0, 0.0, 0.0],
            "",
            "",
        ),
        record(
            "HETATM",
            "2",
            " CL1",
            "LIG",
            'A',
            "1",
            [1.75, 0.0, 0.0],
            "",
            "",
        ),
    ]);
    let error = pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::InvalidValue);
    assert!(error.to_string().contains("\"CL1\""), "{error}");
}

#[test]
fn component_tables_group_identifiers_regardless_of_case() {
    // `LIG` and `lig` are one component to a reader, but bond differently.
    let atom = |serial: &str, name: &str, residue: &str, sequence: &str, x: f64| {
        record(
            "HETATM",
            serial,
            name,
            residue,
            'A',
            sequence,
            [x, 0.0, 0.0],
            "C",
            "",
        )
    };
    let structure = read(&entry(&[
        atom("1", " C1 ", "LIG", "1", 0.0),
        atom("2", " C2 ", "LIG", "1", 1.33),
        atom("3", " C1 ", "lig", "2", 10.0),
        atom("4", " C2 ", "lig", "2", 11.54),
        "CONECT    1    2    2".to_owned(),
        "CONECT    2    1    1".to_owned(),
        "CONECT    3    4".to_owned(),
        "CONECT    4    3".to_owned(),
    ]))
    .into_structure();
    assert_eq!(
        bonds(&structure),
        [(0, 1, BondOrder::Double), (2, 3, BondOrder::Single)]
    );
    let (written, _) = mmcif::write([&structure], Default::default()).unwrap();
    assert!(!written.contains("_chem_comp_bond."), "{written}");
    let reread = read_mmcif(&written).into_structure();
    assert_eq!(bonds(&reread), bonds(&structure));
}

#[test]
fn dative_bonds_without_a_metal_stay_struct_conn_rows() {
    // Ammonia borane: the N-B bond is dative, but neither end is a metal,
    // so a coordination flag would misstate it.
    let cif = "data_x
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
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_model_num
HETATM 1 N N1 NBH A 1 . 0.000 0.000 0.000 1 1
HETATM 2 B B1 NBH A 1 . 1.660 0.000 0.000 1 1
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_auth_seq_id
metalc1 metalc A NBH N1 1 A NBH B1 1
";
    let structure = read_mmcif(cif).into_structure();
    assert_eq!(bonds(&structure), [(0, 1, BondOrder::Dative)]);
    let (written, _) = mmcif::write([&structure], Default::default()).unwrap();
    assert!(!written.contains("_chem_comp_bond."), "{written}");
    assert!(written.contains("metalc"), "{written}");
    let reread = read_mmcif(&written).into_structure();
    assert_eq!(bonds(&reread), bonds(&structure));
}

#[test]
fn aromatic_nitrogens_bond_metals_datively() {
    // Platinum with a chloride and a pyridine, as an mmCIF table without
    // coordination flags states them.
    let cif = "data_x
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
PYP PT1 N1 SING
PYP PT1 CL1 SING
PYP N1 C2 AROM
PYP C2 C3 AROM
PYP C3 C4 AROM
PYP C4 C5 AROM
PYP C5 C6 AROM
PYP C6 N1 AROM
PYP C2 H2 SING
PYP C3 H3 SING
PYP C4 H4 SING
PYP C5 H5 SING
PYP C6 H6 SING
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
_atom_site.pdbx_PDB_model_num
HETATM 1 PT PT1 PYP A 1 . 0.000 0.000 0.000 1
HETATM 2 CL CL1 PYP A 1 . -2.300 0.000 0.000 1
HETATM 3 N N1 PYP A 1 . 2.010 0.000 0.000 1
HETATM 4 C C2 PYP A 1 . 2.705 1.204 0.000 1
HETATM 5 C C3 PYP A 1 . 4.095 1.204 0.000 1
HETATM 6 C C4 PYP A 1 . 4.790 0.000 0.000 1
HETATM 7 C C5 PYP A 1 . 4.095 -1.204 0.000 1
HETATM 8 C C6 PYP A 1 . 2.705 -1.204 0.000 1
HETATM 9 H H2 PYP A 1 . 2.165 2.139 0.000 1
HETATM 10 H H3 PYP A 1 . 4.635 2.139 0.000 1
HETATM 11 H H4 PYP A 1 . 5.870 0.000 0.000 1
HETATM 12 H H5 PYP A 1 . 4.635 -2.139 0.000 1
HETATM 13 H H6 PYP A 1 . 2.165 -2.139 0.000 1
";
    let read = read_mmcif(cif);
    let metal_bonds = bonds(read.structure())
        .into_iter()
        .filter(|&(left, _, _)| left == 0)
        .collect::<Vec<_>>();
    assert_eq!(
        metal_bonds,
        [(0, 1, BondOrder::Single), (0, 2, BondOrder::Dative)]
    );
    assert!(charges(read.structure()).iter().all(|&charge| charge == 0));
}
