//! The bond sources of macromolecular interpretation: connection records,
//! templates, polymer links, and what the report says about each.
use kekule::macromolecular::{
    BondEvidence, SourceResidueId, SourceResiduePosition, StructureInterpretIssue,
    StructureInterpretation,
};
use kekule::molecule::BondOrder;
use kekule::{mmcif, pdb};

/// One `ATOM` or `HETATM` record; `name` is aligned as in columns 13-16.
fn atom(
    record: &str,
    serial: usize,
    name: &str,
    residue: &str,
    sequence: i32,
    [x, y, z]: [f64; 3],
    element: &str,
) -> String {
    format!(
        "{record:<6}{serial:>5} {name:<4} {residue:>3} A{sequence:>4}    {x:>8.3}{y:>8.3}{z:>8.3}  1.00 20.00          {element:>2}"
    )
}

/// A `LINK` record between two atoms of residue `sequence` of chain A.
fn link(atoms: [&str; 2], residue: &str, sequence: i32) -> String {
    format!(
        "{:<12}{:<4} {residue:>3} A{sequence:>4}{:16}{:<4} {residue:>3} A{sequence:>4}{:3}{:>6} {:>6} {:>5.2}",
        "LINK", atoms[0], "", atoms[1], "", "1555", "1555", 1.5,
    )
}

/// `CONECT` records listing each partner of `bonds` from both atoms, each
/// bond `count` times.
fn conect(bonds: &[(usize, usize, usize)]) -> Vec<String> {
    let mut partners = std::collections::BTreeMap::<usize, Vec<usize>>::new();
    for &(left, right, count) in bonds {
        for _ in 0..count {
            partners.entry(left).or_default().push(right);
            partners.entry(right).or_default().push(left);
        }
    }
    partners
        .into_iter()
        .flat_map(|(atom, partners)| {
            partners
                .chunks(4)
                .map(|chunk| {
                    let listed = chunk
                        .iter()
                        .map(|partner| format!("{partner:>5}"))
                        .collect::<String>();
                    format!("CONECT{atom:>5}{listed}")
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

fn read_pdb(lines: &[String]) -> StructureInterpretation {
    let text = lines.join("\n") + "\nEND\n";
    pdb::parse_str(&text, Default::default())
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

/// Bonds by atom names, in name order, with order and evidence, sorted.
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
            let order = molecule.bond(id).unwrap().order();
            (a, b, order, bond.evidence().clone())
        })
        .collect::<Vec<_>>();
    bonds.sort_by(|left, right| (&left.0, &left.1).cmp(&(&right.0, &right.1)));
    bonds
}

fn orders(result: &StructureInterpretation) -> Vec<(String, String, BondOrder)> {
    bonds(result)
        .into_iter()
        .map(|(a, b, order, _)| (a, b, order))
        .collect()
}

/// Nonzero formal charges by atom name, in instance order.
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

/// Report issues other than the entity kinds minimal sources omit.
fn issues(result: &StructureInterpretation) -> Vec<StructureInterpretIssue> {
    result
        .report()
        .issues()
        .iter()
        .filter(|issue| !matches!(issue, StructureInterpretIssue::EntityTypeInferred { .. }))
        .cloned()
        .collect()
}

fn polymer_links(result: &StructureInterpretation) -> Vec<(String, String)> {
    bonds(result)
        .into_iter()
        .filter(|bond| bond.3 == BondEvidence::PolymerLink)
        .map(|(a, b, _, _)| (a, b))
        .collect()
}

fn author(sequence: i32) -> SourceResidueId {
    SourceResidueId::from_parts(
        "1",
        "A",
        SourceResiduePosition::Author(sequence.to_string()),
        None,
    )
}

fn label(sequence: i32) -> SourceResidueId {
    SourceResidueId::from_parts("1", "A", SourceResiduePosition::Label(sequence), None)
}

fn pair(a: &str, b: &str, order: BondOrder) -> (String, String, BondOrder) {
    (a.to_owned(), b.to_owned(), order)
}

/// A benzene ring of six carbons, `C1` to `C6`, without hydrogens.
fn benzene() -> Vec<String> {
    (0..6)
        .map(|index| {
            let angle = f64::from(index) * std::f64::consts::FRAC_PI_3;
            let position = [1.39 * angle.cos(), 1.39 * angle.sin(), 0.0];
            let name = format!(" C{}", index + 1);
            atom("HETATM", index as usize + 1, &name, "LIG", 1, position, "C")
        })
        .collect()
}

#[test]
fn conect_records_without_repeats_state_no_bond_orders() {
    // RCSB files list each bond once from each atom, so a ligand without a
    // template keeps its listed bonds and the report says their orders are
    // unknown.
    let ring = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 6)];
    let once = ring.map(|(left, right)| (left, right, 1));
    let read = read_pdb(&[benzene(), conect(&once)].concat());
    assert!(orders(&read)
        .iter()
        .all(|(_, _, order)| *order == BondOrder::Single));
    assert_eq!(
        issues(&read),
        [StructureInterpretIssue::ResidueBondOrdersUnstated {
            residue: author(1),
            component_id: "LIG".to_owned(),
            bonds: 6,
        }]
    );
    assert_eq!(
        read.report().issues()[0].to_string(),
        "Author residue \"1\" of chain \"A\" in model \"1\" (component \"LIG\") has no \
         template or bond table, and its bond records state no orders. Its 6 bonds between \
         heavy atoms are read as single."
    );

    // Files that list double bonds twice, as Maestro and Kekule write, keep
    // their orders.
    let kekule = ring.map(|(left, right)| {
        (
            left,
            right,
            1 + usize::from(left % 2 == 1 && right == left + 1),
        )
    });
    let read = read_pdb(&[benzene(), conect(&kekule)].concat());
    let double = orders(&read)
        .iter()
        .filter(|(_, _, order)| *order == BondOrder::Double)
        .count();
    assert_eq!(double, 3);
    assert_eq!(issues(&read), []);
}

/// A sulfate ion as RCSB writes it, numbered from 1.
fn sulfate() -> Vec<String> {
    [
        (" S  ", [16.485, 34.924, 5.575], "S"),
        (" O1 ", [17.741, 34.891, 4.826], "O"),
        (" O2 ", [16.787, 35.826, 6.733], "O"),
        (" O3 ", [15.440, 35.551, 4.401], "O"),
        (" O4 ", [16.141, 33.720, 5.930], "O"),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, (name, position, element))| {
        atom("HETATM", index + 1, name, "SO4", 1, position, element)
    })
    .collect()
}

#[test]
fn records_within_a_templated_residue_are_checked_against_its_template() {
    let listed = conect(&[(1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1)]);
    // A LINK within the residue the template agrees with: the template
    // still applies, with its orders and charges.
    let read = read_pdb(
        &[
            sulfate(),
            vec![link([" S", " O1"], "SO4", 1)],
            listed.clone(),
        ]
        .concat(),
    );
    assert_eq!(
        orders(&read),
        [
            pair("O1", "S", BondOrder::Double),
            pair("O2", "S", BondOrder::Double),
            pair("O3", "S", BondOrder::Single),
            pair("O4", "S", BondOrder::Single),
        ]
    );
    assert_eq!(
        charges(&read),
        [("O3".to_owned(), -1), ("O4".to_owned(), -1)]
    );
    assert_eq!(issues(&read), []);

    // A LINK the template lacks is kept and reported.
    let read = read_pdb(&[sulfate(), vec![link([" O1", " O2"], "SO4", 1)], listed].concat());
    assert!(orders(&read).contains(&pair("O1", "O2", BondOrder::Single)));
    assert_eq!(
        issues(&read),
        [StructureInterpretIssue::BondNotInTemplate {
            residue: author(1),
            component_id: "SO4".to_owned(),
            atoms: ["O1".to_owned(), "O2".to_owned()],
        }]
    );
    assert_eq!(
        read.report().issues()[0].to_string(),
        "A connection record bonds atoms \"O1\" and \"O2\" of author residue \"1\" of chain \
         \"A\" in model \"1\" (component \"SO4\"). Its template does not bond them."
    );
}

/// A glycine dipeptide of chain A, entity 1, labelled 1 and 2.
const DIGLYCINE: &[(&str, i32, [f64; 3], &str)] = &[
    ("N", 1, [0.000, 0.000, 0.0], "N"),
    ("CA", 1, [1.458, 0.000, 0.0], "C"),
    ("C", 1, [2.009, 1.420, 0.0], "C"),
    ("O", 1, [1.251, 2.390, 0.0], "O"),
    ("N", 2, [3.332, 1.536, 0.0], "N"),
    ("CA", 2, [3.988, 2.836, 0.0], "C"),
    ("C", 2, [5.488, 2.736, 0.0], "C"),
    ("O", 2, [6.100, 1.680, 0.0], "O"),
];

/// An mmCIF block of polypeptide entity 1 in chain A, with `extra`
/// categories and `atoms` as `(name, component, label_seq_id, position,
/// element)`.
fn polypeptide(
    nstd_linkage: &str,
    atoms: &[(&str, &str, i32, [f64; 3], &str)],
    extra: &str,
) -> String {
    let mut text = format!(
        "data_test\nloop_\n_entity.id\n_entity.type\n1 polymer\nloop_\n_entity_poly.entity_id\n\
         _entity_poly.type\n_entity_poly.nstd_linkage\n1 'polypeptide(L)' {nstd_linkage}\n\
         loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\n{extra}\nloop_\n\
         _atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
         _atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n\
         _atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
    );
    for (index, (name, component, sequence, [x, y, z], element)) in atoms.iter().enumerate() {
        text.push_str(&format!(
            "ATOM {} {element} {name} {component} A 1 {sequence} {x} {y} {z}\n",
            index + 1
        ));
    }
    text
}

fn glycines() -> Vec<(&'static str, &'static str, i32, [f64; 3], &'static str)> {
    DIGLYCINE
        .iter()
        .map(|&(name, sequence, position, element)| (name, "GLY", sequence, position, element))
        .collect()
}

/// A `_struct_conn` loop of `rows`: `id type atom_1 seq_1 atom_2 seq_2
/// order`, all in chain A.
fn struct_conn(rows: &[&str]) -> String {
    let mut text = "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
                    _struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_atom_id\n\
                    _struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr2_label_asym_id\n\
                    _struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_label_seq_id\n\
                    _struct_conn.pdbx_value_order\n"
        .to_owned();
    for row in rows {
        let fields = row.split_whitespace().collect::<Vec<_>>();
        text.push_str(&format!(
            "{} {} A {} {} A {} {} {}\n",
            fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6]
        ));
    }
    text
}

#[test]
fn struct_conn_rows_within_a_residue_keep_its_template() {
    // A row the template agrees with: the template adds the others.
    let read = read_mmcif(&polypeptide(
        "no",
        &glycines()[..4],
        &struct_conn(&["c1 covale N 1 CA 1 ?"]),
    ));
    assert_eq!(
        orders(&read),
        [
            pair("C", "CA", BondOrder::Single),
            pair("C", "O", BondOrder::Double),
            pair("CA", "N", BondOrder::Single),
        ]
    );
    assert_eq!(issues(&read), []);

    // A stated order the template does not give is kept and reported.
    let read = read_mmcif(&polypeptide(
        "no",
        &glycines()[..4],
        &struct_conn(&["c1 covale C 1 O 1 sing"]),
    ));
    assert!(orders(&read).contains(&pair("C", "O", BondOrder::Single)));
    assert_eq!(
        issues(&read),
        [StructureInterpretIssue::BondOrderDiffersFromTemplate {
            residue: label(1),
            component_id: "GLY".to_owned(),
            atoms: ["C".to_owned(), "O".to_owned()],
            order: BondOrder::Single,
            template_order: Some(BondOrder::Double),
        }]
    );
    assert_eq!(
        read.report().issues()[0].to_string(),
        "A connection record states a single bond between atoms \"C\" and \"O\" of residue 1 \
         of chain \"A\" in model \"1\" (component \"GLY\"). Its template has a double bond."
    );
}

#[test]
fn nonstandard_linkage_keeps_the_standard_links_no_record_contradicts() {
    // Head to tail: the record closes the ring, and the standard link
    // between the two residues stays.
    let cyclic = read_mmcif(&polypeptide(
        "yes",
        &glycines(),
        &struct_conn(&["c1 covale C 2 N 1 ?"]),
    ));
    assert_eq!(polymer_links(&cyclic), [("C".to_owned(), "N".to_owned())]);
    assert_eq!(issues(&cyclic), []);

    // A beta-aspartyl link: glycine's N is bonded to aspartate's CG, so the
    // standard link from aspartate's C is not made.
    let atoms = [
        ("N", "ASP", 1, [-0.317, 1.688, 0.066], "N"),
        ("CA", "ASP", 1, [-0.470, 0.286, -0.344], "C"),
        ("C", "ASP", 1, [-1.868, -0.180, -0.029], "C"),
        ("O", "ASP", 1, [-2.534, 0.415, 0.786], "O"),
        ("CB", "ASP", 1, [0.539, -0.580, 0.413], "C"),
        ("CG", "ASP", 1, [1.938, -0.195, 0.004], "C"),
        ("OD1", "ASP", 1, [2.109, 0.681, -0.810], "O"),
        ("OXT", "ASP", 1, [-2.374, -1.256, -0.652], "O"),
        ("N", "GLY", 2, [2.992, -0.826, 0.543], "N"),
        ("CA", "GLY", 2, [4.132, -1.509, 1.126], "C"),
        ("C", "GLY", 2, [5.632, -1.509, 1.126], "C"),
        ("O", "GLY", 2, [6.250, -0.500, 1.126], "O"),
    ];
    for nstd_linkage in ["yes", "no"] {
        let branched = read_mmcif(&polypeptide(
            nstd_linkage,
            &atoms,
            &struct_conn(&["c1 covale CG 1 N 2 sing"]),
        ));
        assert_eq!(polymer_links(&branched), []);
        assert_eq!(
            issues(&branched),
            [StructureInterpretIssue::PolymerLinkSkipped {
                residue: label(1),
                next: label(2),
            }]
        );
    }
}

/// The records of [`DIGLYCINE`] as PDB `ATOM` records, residue 2 numbered
/// `second`.
fn glycine_records(second: i32) -> Vec<String> {
    DIGLYCINE
        .iter()
        .enumerate()
        .map(|(index, &(name, sequence, position, element))| {
            let sequence = if sequence == 1 { 1 } else { second };
            atom(
                "ATOM",
                index + 1,
                &format!(" {name}"),
                "GLY",
                sequence,
                position,
                element,
            )
        })
        .collect()
}

#[test]
fn only_polymer_residues_take_positions_along_a_pdb_chain() {
    // A ligand listed between two residues of a chain is not part of it.
    let mut lines = glycine_records(2);
    lines.splice(
        4..4,
        [
            atom("HETATM", 20, " C1 ", "LIG", 900, [20.0, 0.0, 0.0], "C"),
            atom("HETATM", 21, " O1 ", "LIG", 900, [21.4, 0.0, 0.0], "O"),
        ],
    );
    let read = read_pdb(&lines);
    assert_eq!(polymer_links(&read), [("C".to_owned(), "N".to_owned())]);
    assert_eq!(read.topology().instance_count(), 2);

    // A modified residue in HETATM records is: its template is a peptide's.
    let mut lines = glycine_records(3);
    lines.splice(
        4..4,
        [
            (" N  ", [3.332, 1.536, 0.0], "N"),
            (" CA ", [3.988, 2.836, 0.0], "C"),
            (" C  ", [5.488, 2.736, 0.0], "C"),
            (" O  ", [6.100, 3.792, 0.0], "O"),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, (name, position, element))| {
            atom("HETATM", 20 + index, name, "MSE", 2, position, element)
        }),
    );
    for line in &mut lines[8..] {
        // Residue 3 follows residue 2 along x.
        let x = line[30..38].trim().parse::<f64>().unwrap() + 3.479;
        line.replace_range(30..38, &format!("{x:>8.3}"));
    }
    let read = read_pdb(&lines);
    assert_eq!(polymer_links(&read).len(), 2);
    assert_eq!(read.topology().instance_count(), 1);
}

#[test]
fn absent_linking_atoms_break_the_chain() {
    // Residue 1 without its carbonyl, as a model that leaves it out.
    let mut lines = glycine_records(2);
    lines.drain(2..4);
    let read = read_pdb(&lines);
    assert_eq!(polymer_links(&read), []);
    assert_eq!(
        issues(&read),
        [StructureInterpretIssue::LinkingAtomMissing {
            residue: author(1),
            next: author(2),
        }]
    );
    assert_eq!(
        read.report().issues()[0].to_string(),
        "Author residue \"1\" of chain \"A\" in model \"1\" and author residue \"2\" of chain \
         \"A\" in model \"1\" are not bonded. One of their linking atoms is missing."
    );
}

#[test]
fn charmm_atom_names_under_dictionary_names_read_with_charmm_templates() {
    // A C-terminal isoleucine as CHARMM names it: CD for CD1, and OT1 and
    // OT2 for the carboxyl oxygens.
    let lines = [
        (" N  ", [52.625, 76.235, 68.049], "N"),
        (" CA ", [52.964, 77.620, 67.705], "C"),
        (" C  ", [51.910, 78.234, 66.791], "C"),
        (" OT1", [51.409, 77.508, 65.911], "O"),
        (" CB ", [54.346, 77.727, 66.970], "C"),
        (" CG1", [54.852, 79.179, 66.992], "C"),
        (" CG2", [54.218, 77.237, 65.524], "C"),
        (" CD ", [56.126, 79.382, 66.170], "C"),
        (" OT2", [51.631, 79.444, 66.958], "O"),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, (name, position, element))| {
        atom("ATOM", index + 1, name, "ILE", 1, position, element)
    })
    .collect::<Vec<_>>();
    let read = read_pdb(&lines);
    assert!(orders(&read).contains(&pair("CD", "CG1", BondOrder::Single)));
    assert!(orders(&read).contains(&pair("C", "OT1", BondOrder::Double)));
    assert_eq!(read.topology().instance_count(), 1);
    assert_eq!(
        issues(&read),
        [StructureInterpretIssue::ForceFieldTemplate {
            residue: author(1),
            component_id: "ILE".to_owned(),
            force_field: "CHARMM".to_owned(),
        }]
    );
    assert_eq!(
        read.report().issues()[0].to_string(),
        "Author residue \"1\" of chain \"A\" in model \"1\" (component \"ILE\") fits no \
         dictionary template. It is read with the CHARMM template of that name."
    );
}

#[test]
fn mmcif_writing_fails_on_residues_that_lack_a_template_bond() {
    use kekule::diagnostics::{FormatErrorDetail, UnrepresentablePolicy, WriteOmission};
    use kekule::mmcif::{MmcifWriteErrorKind, MmcifWriteOptions};

    // The block's own table leaves out glycine's C=O bond, which a reader
    // of the written block would add back from the template.
    let table = "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n\
                 _chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nGLY N CA sing\n\
                 GLY CA C sing\n";
    let read = read_mmcif(&polypeptide("no", &glycines()[..4], table));
    assert_eq!(orders(&read).len(), 2);
    let error = mmcif::write([&read], MmcifWriteOptions::default()).unwrap_err();
    assert!(matches!(
        error.detail(),
        Some(FormatErrorDetail::MmcifWrite(
            MmcifWriteErrorKind::DictionaryBondMissing { .. }
        ))
    ));
    let (_, report) = mmcif::write(
        [&read],
        MmcifWriteOptions::default().with_chemistry(UnrepresentablePolicy::Omit),
    )
    .unwrap();
    assert_eq!(
        report.omissions(),
        [WriteOmission::DictionaryBonds { residues: 1 }]
    );
}
