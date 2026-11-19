use kekule::macromolecular::{
    BondEvidence, ConnectivityOptions, StructureInterpretIssue, StructureInterpretOptions,
    StructureInterpretation,
};
use kekule::molecule::BondOrder;
use kekule::structure::Structure;
use kekule::{mmcif, pdb};

/// Two cysteines of chains A and B joined by a disulfide: the dictionary's
/// ideal `CYS` without `HG`, and a copy turned about the S-S bond's
/// midpoint, so that the bond is 2.04 Å with CB-SG-SG angles of 104°.
const CYSTINE: &[(char, &str, &str, [f64; 3])] = &[
    ('A', "N", "N", [1.585, 0.483, -0.081]),
    ('A', "CA", "C", [0.141, 0.450, 0.186]),
    ('A', "C", "C", [-0.095, 0.006, 1.606]),
    ('A', "O", "O", [0.685, -0.742, 2.143]),
    ('A', "CB", "C", [-0.533, -0.530, -0.774]),
    ('A', "SG", "S", [-0.247, 0.004, -2.484]),
    ('A', "OXT", "O", [-1.174, 0.443, 2.275]),
    ('A', "H", "H", [1.928, -0.454, 0.063]),
    ('A', "H2", "H", [1.693, 0.682, -1.065]),
    ('A', "HA", "H", [-0.277, 1.446, 0.042]),
    ('A', "HB2", "H", [-0.114, -1.526, -0.630]),
    ('A', "HB3", "H", [-1.604, -0.554, -0.575]),
    ('A', "HXT", "H", [-1.326, 0.158, 3.186]),
    ('B', "N", "N", [-0.434, -3.286, -6.013]),
    ('B', "CA", "C", [-1.206, -2.066, -6.281]),
    ('B', "C", "C", [-0.968, -1.623, -7.700]),
    ('B', "O", "O", [0.087, -1.858, -8.237]),
    ('B', "CB", "C", [-0.764, -0.961, -5.320]),
    ('B', "SG", "S", [-1.050, -1.495, -3.610]),
    ('B', "OXT", "O", [-1.929, -0.967, -8.369]),
    ('B', "H", "H", [0.537, -3.052, -6.157]),
    ('B', "H2", "H", [-0.539, -3.486, -5.029]),
    ('B', "HA", "H", [-2.267, -2.270, -6.136]),
    ('B', "HB2", "H", [0.297, -0.758, -5.464]),
    ('B', "HB3", "H", [-1.338, -0.057, -5.519]),
    ('B', "HXT", "H", [-1.776, -0.683, -9.280]),
];

/// The serials of the two `SG` atoms in [`cystine_pdb`] with hydrogens.
const SULFURS: [usize; 2] = [6, 19];

/// The cystine as PDB records of residue name `residue`, with or without its
/// hydrogens, chain B shifted by `shift` Å along x, and `records` before
/// `END`.
fn cystine_pdb(residue: &str, hydrogens: bool, shift: f64, records: &[String]) -> String {
    let mut lines = Vec::new();
    let atoms = CYSTINE
        .iter()
        .filter(|(_, _, element, _)| hydrogens || *element != "H");
    for (serial, &(chain, name, element, [x, y, z])) in atoms.enumerate() {
        let x = if chain == 'B' { x + shift } else { x };
        lines.push(format!(
            "ATOM  {:>5} {:<4} {residue:>3} {chain}   1    {x:>8.3}{y:>8.3}{z:>8.3}  1.00 20.00          {element:>2}",
            serial + 1,
            format!(" {name}"),
        ));
    }
    lines.extend_from_slice(records);
    lines.push("END".to_owned());
    lines.join("\n") + "\n"
}

/// The cystine as an mmCIF block whose `_struct_conn` states the disulfide.
fn cystine_mmcif() -> String {
    let mut text = String::from(
        "data_cystine\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
         _atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
         _atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.auth_seq_id\n\
         _atom_site.auth_asym_id\n_atom_site.pdbx_PDB_model_num\n",
    );
    for (serial, &(chain, name, element, [x, y, z])) in CYSTINE.iter().enumerate() {
        text.push_str(&format!(
            "ATOM {} {element} {name} CYS {chain} 1 . {x} {y} {z} 1 {chain} 1\n",
            serial + 1
        ));
    }
    text.push_str(
        "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n\
         _struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
         _struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_comp_id\n\
         _struct_conn.ptnr2_label_atom_id\ndisulf1 disulf A CYS SG B CYS SG\n",
    );
    text
}

fn ssbond(residue: &str) -> String {
    format!("SSBOND   1 {residue} A    1    {residue} B    1                          1555   1555  2.04  ")
}

fn conect() -> Vec<String> {
    let [a, b] = SULFURS;
    vec![format!("CONECT{a:>5}{b:>5}"), format!("CONECT{b:>5}{a:>5}")]
}

fn read_pdb(text: &str, options: StructureInterpretOptions) -> StructureInterpretation {
    pdb::parse_str(text, Default::default())
        .unwrap()
        .interpret(options)
        .unwrap()
}

fn read_mmcif(text: &str) -> StructureInterpretation {
    mmcif::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// The evidence of each bond between two `SG` atoms.
fn disulfides(result: &StructureInterpretation) -> Vec<BondEvidence> {
    let atoms = result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .collect::<Vec<_>>();
    let name = |atom| {
        atoms
            .iter()
            .find(|candidate| candidate.atom() == atom)
            .unwrap()
            .atom_name()
    };
    result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.bonds())
        .filter(|bond| bond.atoms().map(name) == ["SG", "SG"])
        .map(|bond| bond.evidence().clone())
        .collect()
}

/// The structure perceived, with each charged atom's formal charge.
fn perceived(result: &StructureInterpretation) -> (Structure, Vec<i8>) {
    let mut structure = result.structure().clone();
    structure.perceive().unwrap();
    let charges = structure
        .topology()
        .atoms()
        .map(|atom| atom.formal_charge())
        .filter(|&charge| charge != 0)
        .collect();
    (structure, charges)
}

/// Asserts one molecule of neutral atoms, none with implicit hydrogens.
fn assert_neutral_with_explicit_hydrogens(result: &StructureInterpretation) {
    let (structure, charges) = perceived(result);
    assert_eq!(structure.topology().instance_count(), 1);
    assert_eq!(charges, Vec::<i8>::new());
    assert!(structure
        .topology()
        .atoms()
        .all(|atom| atom.implicit_hydrogens() == Some(0)));
}

fn inferred(result: &StructureInterpretation) -> usize {
    result
        .report()
        .issues()
        .iter()
        .filter(|issue| matches!(issue, StructureInterpretIssue::DisulfideInferred { .. }))
        .count()
}

fn connection(id: Option<&str>, kind: &str) -> BondEvidence {
    BondEvidence::Connection {
        connection_id: id.map(str::to_owned),
        connection_type: kind.to_owned(),
    }
}

#[test]
fn ssbond_conect_and_bonding_distance_give_neutral_cystine() {
    for (records, evidence, inferences) in [
        (vec![ssbond("CYS")], connection(Some("1"), "SSBOND"), 0),
        (conect(), connection(None, "CONECT"), 0),
        // Without records, SG atoms within bonding distance that carry no
        // hydrogen are a disulfide.
        (vec![], BondEvidence::Disulfide, 1),
    ] {
        let result = read_pdb(&cystine_pdb("CYS", true, 0.0, &records), Default::default());
        assert_eq!(disulfides(&result), vec![evidence]);
        assert_eq!(inferred(&result), inferences);
        // The disulfide takes the place of each cysteine's HG: no thiolate.
        assert_neutral_with_explicit_hydrogens(&result);
    }
}

#[test]
fn mmcif_struct_conn_disulfide_gives_neutral_cystine() {
    let result = read_mmcif(&cystine_mmcif());
    assert_eq!(
        disulfides(&result),
        vec![connection(Some("disulf1"), "disulf")]
    );
    assert_eq!(inferred(&result), 0);
    assert_neutral_with_explicit_hydrogens(&result);
}

#[test]
fn cyx_residues_read_with_the_amber_template() {
    for records in [vec![], vec![ssbond("CYX")]] {
        let result = read_pdb(&cystine_pdb("CYX", true, 0.0, &records), Default::default());
        let templates = result
            .report()
            .issues()
            .iter()
            .filter_map(|issue| match issue {
                StructureInterpretIssue::ForceFieldTemplate {
                    component_id,
                    force_field,
                    ..
                } => Some((component_id.as_str(), force_field.as_str())),
                _ => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(templates, [("CYX", "AMBER"), ("CYX", "AMBER")]);
        assert_eq!(disulfides(&result).len(), 1);
        assert_eq!(inferred(&result), usize::from(records.is_empty()));
        assert_neutral_with_explicit_hydrogens(&result);
        // The template gives the carbonyl its double bond.
        let (structure, _) = perceived(&result);
        let doubles = structure
            .topology()
            .bonds()
            .filter(|bond| bond.bond().order() == BondOrder::Double)
            .count();
        assert_eq!(doubles, 2);
    }
}

#[test]
fn force_field_templates_take_hydrogens_they_do_not_name() {
    // N-terminal residues of AMBER files name their amine hydrogens H1 to
    // H3, which the plain CYX template does not.
    let text = cystine_pdb("CYX", true, 0.0, &[]).replace(" H   CYX", " H1  CYX");
    assert_eq!(text.matches(" H1  CYX").count(), 2);
    let result = read_pdb(&text, Default::default());
    let templates = result
        .report()
        .issues()
        .iter()
        .filter(|issue| matches!(issue, StructureInterpretIssue::ForceFieldTemplate { .. }))
        .count();
    assert_eq!(templates, 2);
    assert_eq!(disulfides(&result).len(), 1);
    assert_neutral_with_explicit_hydrogens(&result);
}

#[test]
fn heavy_atom_cysteines_within_bonding_distance_are_bonded() {
    let result = read_pdb(&cystine_pdb("CYS", false, 0.0, &[]), Default::default());
    assert_eq!(disulfides(&result), vec![BondEvidence::Disulfide]);
    assert_eq!(inferred(&result), 1);
    assert_eq!(result.topology().instance_count(), 1);
}

#[test]
fn cysteines_without_a_disulfide_keep_their_thiolates() {
    let apart = read_pdb(&cystine_pdb("CYS", true, 20.0, &[]), Default::default());
    let no_inference = read_pdb(
        &cystine_pdb("CYS", true, 0.0, &[]),
        StructureInterpretOptions::default().with_connectivity(
            ConnectivityOptions::default().with_infer_bonds_from_distance(false),
        ),
    );
    for result in [apart, no_inference] {
        assert_eq!(disulfides(&result), vec![]);
        assert_eq!(inferred(&result), 0);
        let (structure, charges) = perceived(&result);
        assert_eq!(structure.topology().instance_count(), 2);
        assert_eq!(charges, [-1, -1]);
    }
}

#[test]
fn writers_treat_cyx_as_its_template_and_state_the_disulfide() {
    let written = |residue: &str| {
        let read = read_pdb(&cystine_pdb(residue, true, 0.0, &[]), Default::default());
        let (pdb, _) = pdb::write(read.structure(), Default::default()).unwrap();
        let (mmcif, _) = mmcif::write([read.structure()], Default::default()).unwrap();
        (pdb, mmcif)
    };
    let listings = |text: &str| {
        text.lines()
            .filter(|line| line.starts_with("CONECT"))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    let (cyx_pdb, cyx_mmcif) = written("CYX");
    let (cys_pdb, cys_mmcif) = written("CYS");
    // A reader rebuilds the bonds of CYX from its template as those of CYS
    // from the dictionary, so writers state the same bonds for both.
    assert_eq!(listings(&cyx_pdb), listings(&cys_pdb));
    assert_eq!(cyx_mmcif.replace("CYX", "CYS"), cys_mmcif);
    let connections = cyx_mmcif
        .lines()
        .filter(|line| line.starts_with("conn"))
        .collect::<Vec<_>>();
    assert_eq!(connections.len(), 1, "{cyx_mmcif}");
    assert!(connections[0].contains("SG"));
    for reread in [
        read_pdb(&cyx_pdb, Default::default()),
        read_mmcif(&cyx_mmcif),
    ] {
        assert_eq!(disulfides(&reread).len(), 1);
        assert_eq!(inferred(&reread), 0);
        assert_neutral_with_explicit_hydrogens(&reread);
    }
}
