//! The bond sources of macromolecular interpretation: connection records,
//! templates, polymer links, and what the report says about each.
use kekule::macromolecular::{
    BondEvidence, SourceResidueId, SourceResiduePosition, StructureInterpretIssue,
    StructureInterpretation,
};
use kekule::molecule::BondOrder;
use kekule::{mmcif, pdb};
use std::collections::BTreeMap;

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
            neutral_centres: Vec::new(),
        }]
    );
    assert_eq!(
        read.report().issues()[0].to_string(),
        "Author residue \"1\" of chain \"A\" in model \"1\" (component \"LIG\") is bonded \
         by records that state no orders, and no template gives them. Its 6 bonds between \
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
fn struct_conn_rows_without_orders_leave_onium_centres_neutral() {
    // A sulfur bonded to three carbons is a sulfonium only when its bonds
    // are known to be single.
    let atoms = [
        ("S", "XSM", 1, [0.0, 0.0, 0.0], "S"),
        ("C1", "XSM", 1, [1.8, 0.0, 0.0], "C"),
        ("C2", "XSM", 1, [-0.6, 1.7, 0.0], "C"),
        ("C3", "XSM", 1, [-0.6, -0.85, 1.47], "C"),
    ];
    let rows = |order: &str| {
        ["C1", "C2", "C3"]
            .iter()
            .enumerate()
            .map(|(index, carbon)| format!("c{index} covale S 1 {carbon} 1 {order}"))
            .collect::<Vec<_>>()
    };
    let read = |order: &str| {
        let rows = rows(order);
        let rows = rows.iter().map(String::as_str).collect::<Vec<_>>();
        read_mmcif(&polypeptide("no", &atoms, &struct_conn(&rows)))
    };
    let unstated = read("?");
    assert_eq!(charges(&unstated), []);
    let issue = issues(&unstated)
        .into_iter()
        .find(|issue| {
            matches!(
                issue,
                StructureInterpretIssue::ResidueBondOrdersUnstated { .. }
            )
        })
        .expect("the residue's orders are unstated");
    assert_eq!(
        issue,
        StructureInterpretIssue::ResidueBondOrdersUnstated {
            residue: label(1),
            component_id: "XSM".to_owned(),
            bonds: 3,
            neutral_centres: vec![("S".to_owned(), 1)],
        }
    );
    assert_eq!(
        issue.to_string(),
        "Residue 1 of chain \"A\" in model \"1\" (component \"XSM\") is bonded by records \
         that state no orders, and no template gives them. Its 3 bonds between heavy atoms are \
         read as single. Its onium centres took no charge. Atom \"S\" took 1 hydrogen."
    );

    let stated = read("sing");
    assert_eq!(charges(&stated), [("S".to_owned(), 1)]);
    assert!(!issues(&stated).iter().any(|issue| matches!(
        issue,
        StructureInterpretIssue::ResidueBondOrdersUnstated { .. }
    )));
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
    // So is one with an atom named like a linking atom, far from the chain.
    lines[4] = atom("HETATM", 20, " N  ", "LIG", 900, [20.0, 0.0, 0.0], "N");
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

#[test]
fn residues_bonded_to_a_chain_through_linking_atoms_are_part_of_it() {
    use kekule::components::ComponentLibrary;
    use kekule::macromolecular::{ConnectivityOptions, StructureInterpretOptions};
    use std::sync::Arc;

    // A 5'-terminal GTP before the ATOM records of an RNA chain, without
    // TER records. The dictionary's GTP is no polymer residue, but its O3'
    // lies within bonding distance of the next residue's P.
    let lines = [
        ("HETATM", " C3'", "GTP", 1, [0.000, 0.000, 0.0], "C"),
        ("HETATM", " O3'", "GTP", 1, [1.430, 0.000, 0.0], "O"),
        ("ATOM", " P  ", "G", 2, [3.030, 0.000, 0.0], "P"),
        ("ATOM", " OP1", "G", 2, [3.530, 1.400, 0.0], "O"),
        ("ATOM", " OP2", "G", 2, [3.530, -1.400, 0.0], "O"),
        ("ATOM", " O5'", "G", 2, [4.530, 0.500, 0.9], "O"),
    ]
    .into_iter()
    .enumerate()
    .map(
        |(index, (record, name, residue, sequence, position, element))| {
            atom(
                record,
                index + 1,
                name,
                residue,
                sequence,
                position,
                element,
            )
        },
    )
    .collect::<Vec<_>>();
    let read = read_pdb(&lines);
    assert_eq!(polymer_links(&read), [("O3'".to_owned(), "P".to_owned())]);
    assert_eq!(read.topology().instance_count(), 1);

    // An alanine between OpenMM's AMBER caps, ACE and NME, as HETATM
    // records. AMBER's ACE derives from acetate, no polymer residue.
    let lines = [
        ("HETATM", " CH3", "ACE", 1, [-1.500, 0.000, 0.0], "C"),
        ("HETATM", " C  ", "ACE", 1, [0.000, 0.000, 0.0], "C"),
        ("HETATM", " O  ", "ACE", 1, [0.600, 1.040, 0.0], "O"),
        ("ATOM", " N  ", "ALA", 2, [1.330, 0.000, 0.0], "N"),
        ("ATOM", " CA ", "ALA", 2, [2.790, 0.000, 0.0], "C"),
        ("ATOM", " C  ", "ALA", 2, [3.600, 1.200, 0.0], "C"),
        ("ATOM", " O  ", "ALA", 2, [3.000, 2.240, 0.0], "O"),
        ("ATOM", " CB ", "ALA", 2, [3.300, -1.400, 0.0], "C"),
        ("HETATM", " N  ", "NME", 3, [4.930, 1.200, 0.0], "N"),
        ("HETATM", " C  ", "NME", 3, [5.700, 2.400, 0.0], "C"),
    ]
    .into_iter()
    .enumerate()
    .map(
        |(index, (record, name, residue, sequence, position, element))| {
            atom(
                record,
                index + 1,
                name,
                residue,
                sequence,
                position,
                element,
            )
        },
    )
    .collect::<Vec<_>>();
    let text = lines.join("\n") + "\nEND\n";
    let options = StructureInterpretOptions::default().with_connectivity(
        ConnectivityOptions::default().with_components(Some(Arc::new(ComponentLibrary::amber()))),
    );
    let read = pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(options)
        .unwrap();
    assert_eq!(
        polymer_links(&read),
        [
            ("C".to_owned(), "N".to_owned()),
            ("C".to_owned(), "N".to_owned())
        ]
    );
    assert_eq!(read.topology().instance_count(), 1);
}

#[test]
fn residues_without_backbone_atoms_report_no_missing_links() {
    // Ions in ATOM records of one chain, as a GROMACS box writes them, and
    // a CA-only model: neither has linking atoms to miss.
    let ions = (1..=3)
        .map(|index| {
            let x = 5.0 * f64::from(index);
            atom(
                "ATOM",
                index as usize,
                "NA  ",
                "NA",
                index,
                [x, 0.0, 0.0],
                "NA",
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(issues(&read_pdb(&ions)), []);
    let trace = (1..=3)
        .map(|index| {
            let x = 3.8 * f64::from(index);
            atom(
                "ATOM",
                index as usize,
                " CA ",
                "ALA",
                index,
                [x, 0.0, 0.0],
                "C",
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(issues(&read_pdb(&trace)), []);

    // A record bonds residue 1, which lacks its carbonyl, to residue 2.
    let read = read_mmcif(&polypeptide(
        "no",
        &[&glycines()[..2], &glycines()[4..]].concat(),
        &struct_conn(&["c1 covale CA 1 N 2 sing"]),
    ));
    assert!(orders(&read).contains(&pair("CA", "N", BondOrder::Single)));
    assert_eq!(issues(&read), []);
}

#[test]
fn conect_orders_are_decided_per_bonded_group() {
    // An RCSB ligand without repeats and a docked one appended with its
    // double bonds listed twice.
    let ring = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 6)];
    let docked = benzene()
        .into_iter()
        .enumerate()
        .map(|(index, line)| {
            let x = line[30..38].trim().parse::<f64>().unwrap() + 20.0;
            let y = line[38..46].trim().parse::<f64>().unwrap();
            atom(
                "HETATM",
                index + 7,
                &line[12..16],
                "LIG",
                2,
                [x, y, 0.0],
                "C",
            )
        })
        .collect::<Vec<_>>();
    let bonds = ring
        .iter()
        .map(|&(left, right)| (left, right, 1))
        .chain(ring.iter().map(|&(left, right)| {
            (
                left + 6,
                right + 6,
                1 + usize::from(left % 2 == 1 && right == left + 1),
            )
        }))
        .collect::<Vec<_>>();
    let read = read_pdb(&[benzene(), docked, conect(&bonds)].concat());
    let double = orders(&read)
        .iter()
        .filter(|(_, _, order)| *order == BondOrder::Double)
        .count();
    assert_eq!(double, 3);
    assert_eq!(
        issues(&read),
        [StructureInterpretIssue::ResidueBondOrdersUnstated {
            residue: author(1),
            component_id: "LIG".to_owned(),
            bonds: 6,
            neutral_centres: Vec::new(),
        }]
    );
}

/// Glycine `N`, `CA`, `C`, and `O` with an extra `CB`, as a modelled
/// mutation leaves it.
fn glycine_with_cb() -> Vec<String> {
    [
        (" N  ", [0.000, 0.000, 0.0], "N"),
        (" CA ", [1.458, 0.000, 0.0], "C"),
        (" C  ", [2.009, 1.420, 0.0], "C"),
        (" O  ", [1.251, 2.390, 0.0], "O"),
        (" CB ", [1.988, -0.774, -1.205], "C"),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, (name, position, element))| {
        atom("ATOM", index + 1, name, "GLY", 1, position, element)
    })
    .collect()
}

#[test]
fn residues_fitting_no_template_take_the_bonds_records_state() {
    // CONECT records bond a glycine with a CB its template lacks.
    let listed = conect(&[(1, 2, 1), (2, 3, 1), (3, 4, 2), (2, 5, 1)]);
    let read = read_pdb(&[glycine_with_cb(), listed].concat());
    assert_eq!(
        orders(&read),
        [
            pair("C", "CA", BondOrder::Single),
            pair("C", "O", BondOrder::Double),
            pair("CA", "CB", BondOrder::Single),
            pair("CA", "N", BondOrder::Single),
        ]
    );
    assert_eq!(
        issues(&read),
        [StructureInterpretIssue::ResidueFitsNoTemplate {
            residue: author(1),
            component_id: "GLY".to_owned(),
        }]
    );
    assert_eq!(
        read.report().issues()[0].to_string(),
        "Author residue \"1\" of chain \"A\" in model \"1\" (component \"GLY\") fits no \
         template of its name. Its bonds are read from its connection records."
    );
    // Written to PDB, it lists every bond, and reads back the same.
    let (text, _) = pdb::write(&read, Default::default()).unwrap();
    let reread = pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(orders(&reread), orders(&read));
    // Without records, or with records that leave heavy atoms out, it
    // still fails to read.
    for listed in [vec![], conect(&[(2, 5, 1)])] {
        let text = [glycine_with_cb(), listed].concat().join("\n") + "\nEND\n";
        assert!(pdb::parse_str(&text, Default::default())
            .unwrap()
            .interpret(Default::default())
            .is_err());
    }
}

#[test]
fn mmcif_writing_round_trips_residues_fitting_no_template() {
    // The block's own table bonds a glycine with a CB and a sulfur for its
    // oxygen; a written block states their bonds for the reader.
    let atoms = [
        ("N", "GLY", 1, [0.000, 0.000, 0.0], "N"),
        ("CA", "GLY", 1, [1.458, 0.000, 0.0], "C"),
        ("C", "GLY", 1, [2.009, 1.420, 0.0], "C"),
        ("O", "GLY", 1, [1.251, 2.390, 0.0], "O"),
        ("CB", "GLY", 1, [1.988, -0.774, -1.205], "C"),
        ("N", "GLY", 2, [3.332, 1.536, 0.0], "N"),
        ("CA", "GLY", 2, [3.988, 2.836, 0.0], "C"),
        ("C", "GLY", 2, [5.488, 2.736, 0.0], "C"),
        ("O", "GLY", 2, [6.100, 1.680, 0.0], "S"),
    ];
    let table = "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n\
                 _chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nGLY N CA sing\n\
                 GLY CA C sing\nGLY C O doub\nGLY CA CB sing\n";
    let read = read_mmcif(&polypeptide("no", &atoms, table));
    let (text, _) = mmcif::write([&read], Default::default()).unwrap();
    let reread = read_mmcif(&text);
    assert_eq!(orders(&reread), orders(&read));
    assert_eq!(
        issues(&reread),
        [
            StructureInterpretIssue::ResidueFitsNoTemplate {
                residue: label(1),
                component_id: "GLY".to_owned(),
            },
            StructureInterpretIssue::ResidueFitsNoTemplate {
                residue: label(2),
                component_id: "GLY".to_owned(),
            },
        ]
    );
}

/// The bundled histidine variant `variant`.
fn histidine_variant(variant: &str) -> &'static kekule::components::Component {
    kekule::components::ComponentLibrary::bundled()
        .variants("HIS")
        .find(|candidate| candidate.variant() == Some(variant))
        .unwrap()
}

const RING: [&str; 5] = ["CG", "ND1", "CE1", "NE2", "CD2"];

/// The ring bonds of histidine variant `variant` as `_struct_conn` rows,
/// each with `order` or, without one, the variant's.
fn ring_rows(variant: &str, order: Option<&str>) -> Vec<String> {
    use kekule::components::ComponentBondOrder;

    let template = histidine_variant(variant);
    template
        .bonds()
        .iter()
        .enumerate()
        .filter_map(|(index, bond)| {
            let [a, b] = bond.atoms().map(|atom| &template.atoms()[atom]);
            (RING.contains(&a.name()) && RING.contains(&b.name())).then(|| {
                let stated = order.unwrap_or(match bond.order() {
                    ComponentBondOrder::Double => "doub",
                    _ => "sing",
                });
                format!("c{index} covale {} 1 {} 1 {stated}", a.name(), b.name())
            })
        })
        .collect()
}

/// A block of one histidine with the atoms of variant `variant` on a grid,
/// and `rows` as `_struct_conn` rows.
fn histidine(variant: &str, rows: &[String]) -> String {
    let atoms = histidine_variant(variant)
        .atoms()
        .iter()
        .filter(|atom| atom.name() != "HXT" && atom.name() != "OXT")
        .enumerate()
        .map(|(index, atom)| {
            let position = [3.0 * (index % 5) as f64, 3.0 * (index / 5) as f64, 0.0];
            (atom.name(), "HIS", 1, position, atom.element().symbol())
        })
        .collect::<Vec<_>>();
    let rows = rows.iter().map(String::as_str).collect::<Vec<_>>();
    polypeptide("no", &atoms, &struct_conn(&rows))
}

/// Asserts that no atom has more bonds than its element and charge allow.
fn assert_valences(read: &StructureInterpretation) {
    let charge = charges(read).into_iter().collect::<BTreeMap<_, _>>();
    let mut valence = BTreeMap::<String, i32>::new();
    for (a, b, order, _) in bonds(read) {
        let order = match order {
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => 1,
        };
        *valence.entry(a).or_default() += order;
        *valence.entry(b).or_default() += order;
    }
    for (name, valence) in valence {
        let base = match name.chars().next().unwrap() {
            'C' => 4,
            'N' => 3,
            'O' => 2,
            _ => 1,
        };
        let allowed = base + i32::from(charge.get(&name).copied().unwrap_or(0)).abs();
        assert!(valence <= allowed, "{name} has valence {valence}");
    }
}

#[test]
fn stated_orders_are_checked_against_the_variant_a_residue_takes() {
    let ring_orders = |read: &StructureInterpretation| {
        orders(read)
            .into_iter()
            .filter(|(a, b, _)| {
                let ring = ["CG", "ND1", "CE1", "NE2", "CD2"];
                ring.contains(&a.as_str()) && ring.contains(&b.as_str())
            })
            .collect::<Vec<_>>()
    };
    let differs = |read: &StructureInterpretation| {
        issues(read)
            .into_iter()
            .filter(|issue| {
                matches!(
                    issue,
                    StructureInterpretIssue::BondOrderDiffersFromTemplate { .. }
                )
            })
            .count()
    };
    let disagreements = |read: &StructureInterpretation| {
        issues(read)
            .into_iter()
            .filter(|issue| {
                matches!(
                    issue,
                    StructureInterpretIssue::StatedOrdersDisagreeWithHydrogens { .. }
                )
            })
            .collect::<Vec<_>>()
    };
    // Rows that state the ND1-protonated ring of a residue with HD1.
    let read = read_mmcif(&histidine(
        "ND1-protonated",
        &ring_rows("ND1-protonated", None),
    ));
    let doubles = ring_orders(&read)
        .iter()
        .filter(|(_, _, order)| *order == BondOrder::Double)
        .count();
    assert_eq!(doubles, 2);
    assert_eq!(differs(&read), 0);
    assert!(disagreements(&read).is_empty());

    // Rows that state every ring bond single keep them, and report the
    // two the variant makes double. No variant gives them, so the report
    // names the residue.
    let read = read_mmcif(&histidine(
        "ND1-protonated",
        &ring_rows("ND1-protonated", Some("sing")),
    ));
    assert!(ring_orders(&read)
        .iter()
        .all(|(_, _, order)| *order == BondOrder::Single));
    assert_eq!(differs(&read), 2);
    assert_eq!(disagreements(&read).len(), 1);
    assert_valences(&read);

    // Stated orders outrank hydrogens. A double ND1=CE1 bond with HD1 and
    // without HE2 is the NE2-protonated ring with a protonated ND1 and a
    // deprotonated NE2, alone or with the whole ring stated. The report
    // names the residue, since the ND1-protonated variant fits its hydrogens
    // without charges.
    let single = vec!["c1 covale ND1 1 CE1 1 doub".to_owned()];
    for rows in [single, ring_rows("NE2-protonated", None)] {
        let read = read_mmcif(&histidine("ND1-protonated", &rows));
        assert!(ring_orders(&read).contains(&pair("CE1", "ND1", BondOrder::Double)));
        assert!(ring_orders(&read).contains(&pair("CE1", "NE2", BondOrder::Single)));
        assert_eq!(
            charges(&read),
            [("ND1".to_owned(), 1), ("NE2".to_owned(), -1)]
        );
        assert_eq!(differs(&read), 0);
        assert_valences(&read);
        let disagreements = disagreements(&read);
        assert_eq!(disagreements.len(), 1, "{disagreements:?}");
        assert!(disagreements[0].to_string().ends_with(
            "(component \"HIS\") has stated bond orders that disagree with its hydrogens. \
             Its charges are derived from the stated orders and the hydrogens."
        ));
    }
}

#[test]
fn charmm_named_residues_keep_dictionary_charges() {
    // A C-terminal aspartate without hydrogens, as CHARMM names it. CHARMM's
    // own ASP is an anion; the dictionary's, which internal residues take,
    // is neutral.
    let lines = [
        (" N  ", [0.000, 0.000, 0.0], "N"),
        (" CA ", [1.458, 0.000, 0.0], "C"),
        (" C  ", [2.009, 1.420, 0.0], "C"),
        (" OT1", [1.251, 2.390, 0.0], "O"),
        (" OT2", [3.250, 1.600, 0.0], "O"),
        (" CB ", [1.988, -0.774, -1.205], "C"),
        (" CG ", [3.497, -0.774, -1.205], "C"),
        (" OD1", [4.100, 0.250, -1.300], "O"),
        (" OD2", [4.100, -1.850, -1.100], "O"),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, (name, position, element))| {
        atom("ATOM", index + 1, name, "ASP", 1, position, element)
    })
    .collect::<Vec<_>>();
    let read = read_pdb(&lines);
    assert_eq!(charges(&read), []);
    assert_eq!(
        issues(&read),
        [StructureInterpretIssue::ForceFieldTemplate {
            residue: author(1),
            component_id: "ASP".to_owned(),
            force_field: "CHARMM".to_owned(),
        }]
    );
}

#[test]
fn writers_state_standard_links_a_reader_would_not_rebuild() {
    // Glycine 2's N is methylated by a residue of its own, so a reader
    // makes no standard link to it: the peptide bond is stated.
    let atoms = [glycines(), vec![("C1", "XME", 3, [3.900, 0.200, 0.0], "C")]].concat();
    let read = read_mmcif(&polypeptide(
        "no",
        &atoms,
        &struct_conn(&["c1 covale C 1 N 2 sing", "c2 covale N 2 C1 3 sing"]),
    ));
    assert!(orders(&read).contains(&pair("C", "N", BondOrder::Single)));
    assert!(orders(&read).contains(&pair("C1", "N", BondOrder::Single)));

    let (text, _) = mmcif::write([&read], Default::default()).unwrap();
    assert_eq!(orders(&read_mmcif(&text)), orders(&read));

    let (text, _) = pdb::write(&read, Default::default()).unwrap();
    let reread = pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(orders(&reread), orders(&read));
}

#[test]
fn hetatm_residues_with_linking_atoms_stay_in_terminated_chains() {
    // An unknown residue between gaps of a chain a TER record ends: its
    // links are stretched, but it keeps its place, so the report breaks the
    // chain on both sides of it.
    let residue = |record: &str, name: &str, sequence: i32, shift: f64, serial: usize| {
        [
            (" N  ", [0.000, 0.000, 0.0], "N"),
            (" CA ", [1.458, 0.000, 0.0], "C"),
            (" C  ", [2.009, 1.420, 0.0], "C"),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, (atom_name, [x, y, z], element))| {
            atom(
                record,
                serial + index,
                atom_name,
                name,
                sequence,
                [x + shift, y, z],
                element,
            )
        })
        .collect::<Vec<_>>()
    };
    let lines = [
        residue("ATOM", "GLY", 1, 0.0, 1),
        residue("HETATM", "XAA", 2, 10.0, 4),
        residue("ATOM", "GLY", 3, 20.0, 7),
        vec!["TER      10      GLY A   3".to_owned()],
    ]
    .concat();
    let read = read_pdb(&lines);
    let breaks = issues(&read)
        .into_iter()
        .filter(|issue| matches!(issue, StructureInterpretIssue::ChainBreak { .. }))
        .collect::<Vec<_>>();
    assert_eq!(
        breaks,
        [
            StructureInterpretIssue::ChainBreak {
                residue: author(1),
                next: author(2),
            },
            StructureInterpretIssue::ChainBreak {
                residue: author(2),
                next: author(3),
            },
        ]
    );
}
