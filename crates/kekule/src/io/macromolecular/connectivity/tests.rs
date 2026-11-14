use super::*;
use crate::topology::MoleculeInstanceId;

const PEPTIDE: &str = r#"
data_peptide
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.nstd_linkage
1 'polypeptide(L)' no
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
A 1 1 GLY
A 1 2 GLY
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
GLY N CA sing
GLY CA C sing
GLY C O doub
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
ATOM 1 N N GLY A 1 1 0.0 0.0 0.0
ATOM 2 C CA GLY A 1 1 1.4 0.0 0.0
ATOM 3 C C GLY A 1 1 2.8 0.0 0.0
ATOM 4 O O GLY A 1 1 3.8 0.0 0.0
ATOM 5 N N GLY A 1 2 4.1 0.0 0.0
ATOM 6 C CA GLY A 1 2 5.5 0.0 0.0
ATOM 7 C C GLY A 1 2 6.9 0.0 0.0
ATOM 8 O O GLY A 1 2 7.9 0.0 0.0
"#;

const BRANCH: &str = r#"
data_branch
loop_
_entity.id
_entity.type
2 branched
loop_
_struct_asym.id
_struct_asym.entity_id
B 2
loop_
_pdbx_entity_branch_list.entity_id
_pdbx_entity_branch_list.comp_id
_pdbx_entity_branch_list.num
_pdbx_entity_branch_list.hetero
2 NAG 1 n
2 GAL 2 n
loop_
_pdbx_branch_scheme.asym_id
_pdbx_branch_scheme.entity_id
_pdbx_branch_scheme.mon_id
_pdbx_branch_scheme.num
_pdbx_branch_scheme.pdb_seq_num
B 2 NAG 1 10
B 2 GAL 2 11
loop_
_pdbx_entity_branch_link.link_id
_pdbx_entity_branch_link.entity_id
_pdbx_entity_branch_link.entity_branch_list_num_1
_pdbx_entity_branch_link.comp_id_1
_pdbx_entity_branch_link.atom_id_1
_pdbx_entity_branch_link.entity_branch_list_num_2
_pdbx_entity_branch_link.comp_id_2
_pdbx_entity_branch_link.atom_id_2
_pdbx_entity_branch_link.value_order
1 2 1 NAG O4 2 GAL C1 sing
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
NAG C4 O4 sing
GAL C1 O5 sing
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C4 NAG B 2 . 10 0.0 0.0 0.0
HETATM 2 O O4 NAG B 2 . 10 1.4 0.0 0.0
HETATM 3 C C1 GAL B 2 . 11 2.8 0.0 0.0
HETATM 4 O O5 GAL B 2 . 11 4.2 0.0 0.0
"#;

#[test]
fn singleton_connectivity_categories_are_syntax_independent() {
    for (input, bonds) in [(PEPTIDE, 7), (BRANCH, 3)] {
        let scalar = crate::test_fixtures::singleton_loops_as_scalars(input);
        let a = crate::mmcif::parse_str(input, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap();
        let b = crate::mmcif::parse_str(&scalar, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap();
        assert_eq!(b.topology().bond_count(), bonds);
        assert_eq!(b.topology().instance_count(), 1);
        assert!(a.topology().same_layout(b.topology()));
        assert_eq!(a.structure().positions(), b.structure().positions());
    }
}

#[test]
fn singleton_component_bond_and_branch_site_are_read_as_scalars() {
    let source = "data_x\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nGLY C O doub\nloop_\n_pdbx_branch_scheme.entity_id\n_pdbx_branch_scheme.asym_id\n_pdbx_branch_scheme.num\n_pdbx_branch_scheme.mon_id\n_pdbx_branch_scheme.auth_seq_num\n2 B 1 NAG 10\n";
    for input in [
        source.to_owned(),
        crate::test_fixtures::singleton_loops_as_scalars(source),
    ] {
        let doc = crate::mmcif::parse_str(&input, Default::default()).unwrap();
        let catalog =
            ConnectivityCatalog::from_block(&doc.blocks()[0], &raw::ConnectivityOptions::default())
                .unwrap();
        let bond = &catalog.component_bonds["GLY"][0];
        assert_eq!((&*bond.atom_1, &*bond.atom_2), ("C", "O"));
        assert!(matches!(
            bond.order,
            SourceBondOrder::Localized(BondOrder::Double)
        ));
        assert_eq!(catalog.branch_sites.len(), 1);
        let site = &catalog.branch_sites[0];
        assert_eq!(
            (
                &*site.entity_id,
                &*site.asym_id,
                site.number,
                &*site.component_id
            ),
            ("2", "B", 1, "NAG")
        );
        assert_eq!(site.author_sequence_id.as_deref(), Some("10"));
    }
}

#[test]
fn component_templates_and_polymer_links_complete_peptide_graph() {
    let document = crate::io::parse_mmcif_str(PEPTIDE, crate::io::MmcifParseOptions::default())
        .expect("parse peptide");
    let interpretation = raw::interpret_mmcif(&document, raw::StructureInterpretOptions::default())
        .expect("interpret peptide");
    let graph = interpretation
        .topology()
        .instance(MoleculeInstanceId::new(0))
        .expect("peptide instance")
        .molecule();
    assert_eq!(graph.bond_count(), 7);
    assert_eq!(
        graph
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Single)
            .count(),
        5
    );
    assert_eq!(
        graph
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .count(),
        2
    );
    assert!(graph.connected_components().len() == 1);
}

#[test]
fn sequence_gap_is_not_bridged_by_a_fake_polymer_bond() {
    let input = PEPTIDE
        .replace("A 1 2 GLY", "A 1 3 GLY")
        .replace("GLY A 1 2 4.1", "GLY A 1 3 4.1")
        .replace("GLY A 1 2 5.5", "GLY A 1 3 5.5")
        .replace("GLY A 1 2 6.9", "GLY A 1 3 6.9")
        .replace("GLY A 1 2 7.9", "GLY A 1 3 7.9");
    let document = crate::io::parse_mmcif_str(&input, crate::io::MmcifParseOptions::default())
        .expect("parse gapped peptide");
    let interpretation = raw::interpret_mmcif(&document, raw::StructureInterpretOptions::default())
        .expect("interpret gapped peptide");
    assert_eq!(interpretation.topology().instance_count(), 2);
    assert!(interpretation
        .topology()
        .instance_records()
        .all(|(instance, _)| interpretation
            .topology()
            .instance(instance)
            .expect("partitioned peptide instance")
            .molecule()
            .connected_components()
            .len()
            == 1));
}

#[test]
fn branch_link_connects_branched_macromolecule() {
    let document = crate::io::parse_mmcif_str(BRANCH, crate::io::MmcifParseOptions::default())
        .expect("parse branch");
    let interpretation = raw::interpret_mmcif(&document, raw::StructureInterpretOptions::default())
        .expect("interpret branch");
    let graph = interpretation
        .topology()
        .instance(MoleculeInstanceId::new(0))
        .expect("branch instance")
        .molecule();
    assert_eq!(graph.bond_count(), 3);
    assert!(graph
        .bonds()
        .all(|(_, bond)| bond.order == BondOrder::Single));
    assert!(graph.connected_components().len() == 1);
}

/// Two copies, B and C, of [`BRANCH`]'s NAG-GAL entity, which a covalent
/// `_struct_conn` from B's `O5` to C's `C4` joins into one molecule
/// group. `rows` are the `_atom_site` rows of [`BRANCH`]'s layout.
fn two_branch_copies(rows: &str) -> String {
    let header = BRANCH
        .split("loop_\n_atom_site.group_PDB")
        .next()
        .unwrap()
        .replace("B 2\n", "B 2\nC 2\n")
        .replace(
            "B 2 GAL 2 11\n",
            "B 2 GAL 2 11\nC 2 NAG 1 20\nC 2 GAL 2 21\n",
        );
    let atom_site = BRANCH.split("loop_\n_atom_site.group_PDB").nth(1).unwrap();
    let columns = atom_site.split("HETATM").next().unwrap();
    format!(
        "{header}loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
             _struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_comp_id\n\
             _struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_label_seq_id\n\
             _struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_comp_id\n\
             _struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_label_seq_id\n\
             covale1 covale B GAL O5 . C NAG C4 .\n\
             loop_\n_atom_site.group_PDB{columns}{rows}"
    )
}

#[test]
fn branch_links_bond_each_copy_within_bonding_distance() {
    let copy_b = "HETATM 1 C C4 NAG B 2 . 10 0.0 0.0 0.0\n\
                      HETATM 2 O O4 NAG B 2 . 10 1.4 0.0 0.0\n";
    let copy_b_gal = "HETATM 3 C C1 GAL B 2 . 11 2.8 0.0 0.0\n\
                          HETATM 4 O O5 GAL B 2 . 11 4.2 0.0 0.0\n";
    let copy_c = "HETATM 5 C C4 NAG C 2 . 20 5.6 0.0 0.0\n\
                      HETATM 6 O O4 NAG C 2 . 20 7.0 0.0 0.0\n";
    let copy_c_gal = "HETATM 7 C C1 GAL C 2 . 21 8.4 0.0 0.0\n\
                          HETATM 8 O O5 GAL C 2 . 21 9.8 0.0 0.0\n";
    // C's galactose 5 Å off its link partner.
    let copy_c_gal_far = "HETATM 7 C C1 GAL C 2 . 21 8.4 5.0 0.0\n\
                              HETATM 8 O O5 GAL C 2 . 21 9.8 5.0 0.0\n";
    let residue = |asym: &str, sequence: &str| SourceResidueId {
        model_id: "1".to_owned(),
        asym_id: asym.to_owned(),
        position: SourceResiduePosition::Author(sequence.to_owned()),
        insertion_code: None,
    };
    for (rows, links, breaks) in [
        // Each copy is linked, not only the first one listed.
        (
            [copy_b, copy_b_gal, copy_c, copy_c_gal].concat(),
            vec![("B", "B"), ("C", "C")],
            vec![],
        ),
        // B's unmodelled NAG leaves its GAL unlinked; C's NAG does not
        // stand in for it.
        (
            [copy_b_gal, copy_c, copy_c_gal].concat(),
            vec![("C", "C")],
            vec![],
        ),
        // A link beyond bonding distance is a reported break.
        (
            [copy_b, copy_b_gal, copy_c, copy_c_gal_far].concat(),
            vec![("B", "B")],
            vec![StructureInterpretIssue::ChainBreak {
                residue: residue("C", "20"),
                next: residue("C", "21"),
            }],
        ),
    ] {
        let text = two_branch_copies(&rows);
        let interpretation = interpret(&text);
        let report = interpretation.report();
        let asym_of = |atom| {
            report
                .instances()
                .iter()
                .flat_map(|instance| instance.atoms())
                .find(|candidate| candidate.atom() == atom)
                .unwrap()
                .asym_id()
                .to_owned()
        };
        let mut actual = report
            .instances()
            .iter()
            .flat_map(|instance| instance.bonds())
            .filter(|bond| *bond.evidence() == BondEvidence::BranchLink)
            .map(|bond| bond.atoms().map(asym_of))
            .collect::<Vec<_>>();
        actual.sort();
        let expected = links
            .iter()
            .map(|(left, right)| [left.to_string(), right.to_string()])
            .collect::<Vec<_>>();
        assert_eq!(actual, expected, "{rows}");
        let actual_breaks = report
            .issues()
            .iter()
            .filter(|issue| matches!(issue, StructureInterpretIssue::ChainBreak { .. }))
            .cloned()
            .collect::<Vec<_>>();
        assert_eq!(actual_breaks, breaks, "{rows}");
    }
}

use crate::components::ComponentLibrary;
use crate::macromolecular::{
    BondEvidence, ConnectivityOptions, SourceResidueId, SourceResiduePosition,
    StructureInterpretIssue, StructureInterpretOptions, StructureInterpretation,
};

/// A minimal block whose `rows` are `type name component asym seq x y z`.
/// Coordinates of standard residues are the dictionary's ideal ones.
fn block(rows: &str, extra: &str) -> String {
    records("HETATM", rows, extra)
}

/// Like [`block`]; `ATOM` records make the chain an inferred polymer.
fn records(group: &str, rows: &str, extra: &str) -> String {
    let mut text = format!(
            "data_test\n{extra}\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        );
    for (index, row) in rows
        .lines()
        .map(str::trim)
        .filter(|row| !row.is_empty())
        .enumerate()
    {
        text.push_str(&format!("{group} {} {row}\n", index + 1));
    }
    text
}

fn interpret(text: &str) -> StructureInterpretation {
    crate::mmcif::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// Report issues other than the entity kinds these minimal blocks omit.
fn issues(result: &StructureInterpretation) -> Vec<StructureInterpretIssue> {
    result
        .report()
        .issues()
        .iter()
        .filter(|issue| !matches!(issue, StructureInterpretIssue::EntityTypeInferred { .. }))
        .cloned()
        .collect()
}

fn line_of(text: &str, needle: &str) -> usize {
    text.lines().position(|line| line.contains(needle)).unwrap() + 1
}

/// Bonds by atom names, order, and evidence.
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
    result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.bonds())
        .map(|bond| {
            let [a, b] = bond.atoms();
            let order = result
                .topology()
                .instance(a.instance())
                .unwrap()
                .molecule()
                .bond_between(a.atom(), b.atom())
                .map(|id| {
                    result
                        .topology()
                        .instance(a.instance())
                        .unwrap()
                        .molecule()
                        .bond(id)
                        .unwrap()
                        .order
                })
                .unwrap();
            (name(a), name(b), order, bond.evidence().clone())
        })
        .collect()
}

fn order_between(result: &StructureInterpretation, a: &str, b: &str) -> BondOrder {
    bonds(result)
        .into_iter()
        .find(|(left, right, _, _)| (left == a && right == b) || (left == b && right == a))
        .unwrap()
        .2
}

/// Nonzero formal charges by atom name.
fn charges(result: &StructureInterpretation) -> Vec<(String, i8)> {
    result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .filter_map(|atom| {
            let charge = result
                .topology()
                .instance(atom.atom().instance())
                .unwrap()
                .molecule()
                .atom(atom.atom().atom())
                .unwrap()
                .formal_charge;
            (charge != 0).then(|| (atom.atom_name().to_owned(), charge))
        })
        .collect()
}

fn residue(asym: &str, sequence: i32) -> SourceResidueId {
    SourceResidueId {
        model_id: "1".to_owned(),
        asym_id: asym.to_owned(),
        position: SourceResiduePosition::Label(sequence),
        insertion_code: None,
    }
}

fn template(component: &str, variant: Option<&str>) -> BondEvidence {
    BondEvidence::Template {
        component_id: component.to_owned(),
        variant: variant.map(str::to_owned),
    }
}

const ASPARTATE: &str = "
        N N ASP A 1 -0.317 1.688 0.066
        C CA ASP A 1 -0.470 0.286 -0.344
        C C ASP A 1 -1.868 -0.180 -0.029
        O O ASP A 1 -2.534 0.415 0.786
        C CB ASP A 1 0.539 -0.580 0.413
        C CG ASP A 1 1.938 -0.195 0.004
        O OD1 ASP A 1 2.109 0.681 -0.810
        O OD2 ASP A 1 2.992 -0.826 0.543
        O OXT ASP A 1 -2.374 -1.256 -0.652
        H H ASP A 1 -0.928 2.289 -0.467
        H H2 ASP A 1 -0.478 1.795 1.056
        H HA ASP A 1 -0.292 0.199 -1.416
        H HB2 ASP A 1 0.419 -0.425 1.485
        H HB3 ASP A 1 0.367 -1.630 0.176
        H HD2 ASP A 1 3.869 -0.545 0.250
        H HXT ASP A 1 -3.275 -1.517 -0.416";

const HISTIDINE: &str = "
        N N HIS A 1 -0.040 -1.210 0.053
        C CA HIS A 1 1.172 -1.709 0.652
        C C HIS A 1 1.083 -3.207 0.905
        O O HIS A 1 0.040 -3.770 1.222
        C CB HIS A 1 1.484 -0.975 1.962
        C CG HIS A 1 2.940 -1.060 2.353
        N ND1 HIS A 1 3.380 -2.075 3.129
        C CD2 HIS A 1 3.960 -0.251 2.046
        C CE1 HIS A 1 4.693 -1.908 3.317
        N NE2 HIS A 1 5.058 -0.801 2.662
        O OXT HIS A 1 2.247 -3.882 0.744
        H H HIS A 1 -0.102 -1.155 -0.950
        H H2 HIS A 1 -0.715 -0.741 0.634
        H HA HIS A 1 1.965 -1.558 -0.089
        H HB2 HIS A 1 1.215 0.087 1.879
        H HB3 HIS A 1 0.859 -1.368 2.775
        H HD1 HIS A 1 2.828 -2.838 3.511
        H HD2 HIS A 1 4.108 0.647 1.479
        H HE1 HIS A 1 5.340 -2.550 3.892
        H HE2 HIS A 1 6.002 -0.428 2.627
        H HXT HIS A 1 2.188 -4.848 0.901";

fn without(rows: &str, atoms: &[&str]) -> String {
    rows.lines()
        .filter(|row| {
            !atoms
                .iter()
                .any(|atom| row.split_whitespace().nth(1) == Some(atom))
        })
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn unstated_coordination_in_a_component_table_follows_template_then_elements() {
    // wwPDB entry files list `_chem_comp_bond` without the coordination
    // flag; the SF4 template flags its Fe-S bonds, so they are dative as
    // a PDB reading of the same structure makes them.
    let table = |flag_column: bool, rows: &str| {
        let mut text = String::from(
                "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n",
            );
        if flag_column {
            text.push_str("_chem_comp_bond.pdbx_metal_coordination_flag\n");
        }
        text.push_str(rows);
        text
    };
    let sf4 = "FE FE1 SF4 A 1 40.988 0.074 22.637
             S S2 SF4 A 1 43.144 0.781 22.394";
    let unstated = interpret(&block(sf4, &table(false, "SF4 FE1 S2 sing\n")));
    assert_eq!(
        bonds(&unstated)
            .iter()
            .map(|(_, _, order, evidence)| (*order, evidence.clone()))
            .collect::<Vec<_>>(),
        [(BondOrder::Dative, BondEvidence::ComponentTable)]
    );
    // A template atom of another element does not lend its flag: two
    // non-metals named like the cluster's atoms keep their single bond.
    let renamed = interpret(&block(
        "C FE1 SF4 A 1 40.988 0.074 22.637
             S S2 SF4 A 1 41.900 0.781 22.394",
        &table(false, "SF4 FE1 S2 sing\n"),
    ));
    assert_eq!(bonds(&renamed)[0].2, BondOrder::Single);
    // A stated flag wins over the template.
    let stated = interpret(&block(sf4, &table(true, "SF4 FE1 S2 sing N\n")));
    assert_eq!(bonds(&stated)[0].2, BondOrder::Single);
    // Without a template, a metal and a non-metal are coordinated, as for
    // PDB `LINK` and `CONECT`; other pairs keep their stated order.
    let unknown = interpret(&block(
        "FE FE XFC A 1 0.000 0.000 0.000
             CL CL1 XFC A 1 2.200 0.000 0.000
             C C1 XFC A 1 -1.900 0.000 0.000
             O O1 XFC A 1 -3.050 0.000 0.000",
        &table(false, "XFC FE CL1 sing\nXFC C1 O1 doub\n"),
    ));
    let orders = bonds(&unknown)
        .into_iter()
        .map(|(a, b, order, _)| (a, b, order))
        .collect::<Vec<_>>();
    assert!(orders.contains(&("FE".to_owned(), "CL1".to_owned(), BondOrder::Dative)));
    assert!(orders.contains(&("C1".to_owned(), "O1".to_owned(), BondOrder::Double)));
}

#[test]
fn templates_supply_bonds_and_metal_coordination_is_dative() {
    // The SF4 template flags every Fe-S bond as metal coordination;
    // its sulfides carry -2 and its iron charges are unknown.
    let cluster = interpret(&block(
        "FE FE1 SF4 A 1 40.988 0.074 22.637
             FE FE2 SF4 A 1 39.644 2.319 22.749
             FE FE3 SF4 A 1 42.036 2.170 23.811
             FE FE4 SF4 A 1 41.755 2.132 21.208
             S S1 SF4 A 1 41.267 3.911 22.551
             S S2 SF4 A 1 43.144 0.781 22.394
             S S3 SF4 A 1 39.806 0.980 20.921
             S S4 SF4 A 1 40.200 1.034 24.551",
        "",
    ));
    assert_eq!(cluster.topology().instance_count(), 1);
    let cluster_bonds = bonds(&cluster);
    assert_eq!(cluster_bonds.len(), 12);
    assert!(cluster_bonds
        .iter()
        .all(|(_, _, order, evidence)| *order == BondOrder::Dative
            && *evidence == template("SF4", None)));
    assert_eq!(
        charges(&cluster),
        [("S1", -2), ("S2", -2), ("S3", -2), ("S4", -2)]
            .map(|(name, charge)| (name.to_owned(), charge))
    );

    // A metal-coordination link between residues is recorded, not
    // bonded: the zinc ion and the cysteine stay separate molecules, and
    // the ion's charge comes from its template.
    let zinc = interpret(&block(
            "N N CYS A 1 1.585 0.483 -0.081
             C CA CYS A 1 0.141 0.450 0.186
             C C CYS A 1 -0.095 0.006 1.606
             O O CYS A 1 0.685 -0.742 2.143
             C CB CYS A 1 -0.533 -0.530 -0.774
             S SG CYS A 1 -0.247 0.004 -2.484
             ZN ZN ZN B . -0.247 0.004 -4.784",
            "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_label_seq_id\nmetalc1 metalc A SG 1 B ZN .\n",
        ));
    assert_eq!(zinc.topology().instance_count(), 2);
    let cysteine = template("CYS", None);
    let as_strings = |rows: &[(&str, &str, BondOrder, BondEvidence)]| {
        rows.iter()
            .map(|(a, b, order, evidence)| {
                ((*a).to_owned(), (*b).to_owned(), *order, evidence.clone())
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(
        bonds(&zinc),
        as_strings(&[
            ("N", "CA", BondOrder::Single, cysteine.clone()),
            ("CA", "C", BondOrder::Single, cysteine.clone()),
            ("CA", "CB", BondOrder::Single, cysteine.clone()),
            ("C", "O", BondOrder::Double, cysteine.clone()),
            ("CB", "SG", BondOrder::Single, cysteine),
        ])
    );
    assert_eq!(charges(&zinc), [("ZN".to_owned(), 2)]);
    let [record] = zinc.report().coordination() else {
        panic!("one coordination record");
    };
    let names = |result: &StructureInterpretation| {
        let [record] = result.report().coordination() else {
            panic!("one coordination record");
        };
        record.atoms().map(|atom| {
            result
                .report()
                .instances()
                .iter()
                .flat_map(|instance| instance.atoms())
                .find(|candidate| candidate.atom() == atom)
                .unwrap()
                .atom_name()
                .to_owned()
        })
    };
    assert_eq!(names(&zinc), ["SG", "ZN"]);
    assert_ne!(record.atoms()[0].instance(), record.atoms()[1].instance());
    assert_eq!(
        (record.connection_id(), record.connection_type()),
        (Some("metalc1"), "metalc")
    );
    assert_eq!(zinc.report().applied_connections(), 1);

    // Written with its report, the coordination is a `metalc` row again;
    // a record naming an atom the structure lacks is left out and reported.
    let written = |report: &crate::macromolecular::StructureInterpretationReport| {
        crate::mmcif::write(
            [crate::mmcif::MmcifBlockSource::structure(zinc.structure())
                .with_reports(std::slice::from_ref(report))],
            Default::default(),
        )
        .unwrap()
    };
    let (text, write_report) = written(zinc.report());
    assert_eq!(text.matches(" metalc ").count(), 1);
    assert!(write_report.omissions().is_empty());
    let reread = interpret(&text);
    assert_eq!(reread.topology().instance_count(), 2);
    assert_eq!(names(&reread), ["SG", "ZN"]);
    let mut dangling = zinc.report().clone();
    let ion = dangling.coordination[0].atoms[1];
    dangling.coordination[0].atoms[1] =
        crate::topology::InstanceAtomId::new(ion.instance(), crate::molecule::AtomId::new(99));
    let (text, write_report) = written(&dangling);
    assert!(!text.contains("metalc"));
    assert_eq!(
        write_report.omissions(),
        [crate::diagnostics::WriteOmission::Coordination { records: 1 }]
    );

    // A structure's own dative bond between residues is written as
    // coordination, which reads back as a record, and is reported so.
    let mut editor = zinc.structure().to_editor();
    let [sulfur, ion] = zinc.report().coordination()[0]
        .atoms()
        .map(|atom| editor.topology().atom_handle(atom).unwrap());
    editor.add_bond(sulfur, ion, BondOrder::Dative).unwrap();
    let bonded = editor.finish().unwrap();
    assert_eq!(bonded.topology().instance_count(), 1);
    let (text, write_report) = crate::mmcif::write([&bonded], Default::default()).unwrap();
    assert_eq!(text.matches(" metalc ").count(), 1);
    assert_eq!(
        write_report.omissions(),
        [crate::diagnostics::WriteOmission::CoordinationBonds { bonds: 1 }]
    );
    let reread = interpret(&text);
    assert_eq!(reread.topology().instance_count(), 2);
    assert_eq!(names(&reread), ["SG", "ZN"]);
}

#[test]
fn model_hydrogens_set_heteroatom_charges_and_select_tautomers() {
    // Without model hydrogens the template's charges apply: the
    // dictionary's histidine is the ND1-protonated imidazolium.
    let heavy = interpret(&block(
        &without(
            HISTIDINE,
            &[
                "H", "H2", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HXT",
            ],
        ),
        "",
    ));
    assert_eq!(charges(&heavy), [("ND1".to_owned(), 1)]);
    assert!(issues(&heavy).is_empty());

    // An extra amine hydrogen named outside the template, or outside the
    // block's own bond table, attaches by distance; with OXT unprotonated
    // the aspartate is a zwitterion.
    let table = "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n\
            _chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n\
            ASP N CA sing\nASP CA C sing\nASP C O doub\nASP C OXT sing\nASP CA CB sing\n\
            ASP CB CG sing\nASP CG OD1 doub\nASP CG OD2 sing\nASP N H sing\nASP N H2 sing\n\
            ASP CA HA sing\nASP CB HB2 sing\nASP CB HB3 sing\nASP OD2 HD2 sing\n";
    for extra in ["", table] {
        let zwitterion = interpret(&block(
            &format!(
                "{}\nH H3 ASP A 1 0.633 1.938 -0.134",
                without(ASPARTATE, &["HXT"])
            ),
            extra,
        ));
        assert_eq!(zwitterion.topology().instance_count(), 1, "{extra}");
        assert_eq!(
            charges(&zwitterion),
            [("N".to_owned(), 1), ("OXT".to_owned(), -1)],
            "{extra}"
        );
        assert_eq!(
            issues(&zwitterion),
            [
                StructureInterpretIssue::HydrogenAttachedByDistance {
                    residue: residue("A", 1),
                    atom_name: "H3".to_owned(),
                    heavy_atom_name: "N".to_owned(),
                },
                StructureInterpretIssue::ProtonationFromHydrogens {
                    residue: residue("A", 1),
                    component_id: "ASP".to_owned(),
                    variant: None,
                    charged_atoms: vec![("N".to_owned(), 1), ("OXT".to_owned(), -1)],
                },
            ],
            "{extra}"
        );
    }

    // A deprotonated side chain carries the charge on its acid oxygen.
    let aspartate = interpret(&block(&without(ASPARTATE, &["HD2"]), ""));
    assert_eq!(charges(&aspartate), [("OD2".to_owned(), -1)]);

    // A ring nitrogen takes or loses one proton whatever its bond orders:
    // the cytosine base (ideal coordinates of the dictionary's DC)
    // protonated at N3, as in i-motif DNA, and the imidazolate histidine.
    let cytosine = "
            N N1 DC A 1 -0.267 -0.584 2.528
            C C2 DC A 1 0.270 0.648 2.563
            O O2 DC A 1 0.052 1.424 1.647
            N N3 DC A 1 1.037 1.035 3.581
            C C4 DC A 1 1.291 0.212 4.589
            N N4 DC A 1 2.085 0.622 5.635
            C C5 DC A 1 0.746 -1.088 4.580
            C C6 DC A 1 -0.035 -1.465 3.541
            H H41 DC A 1 2.461 1.516 5.636
            H H42 DC A 1 2.265 0.019 6.373
            H H5 DC A 1 0.943 -1.771 5.394
            H H6 DC A 1 -0.467 -2.454 3.514
            H H3 DC A 1 1.431 2.041 3.589";
    let imidazolate = without(HISTIDINE, &["HD1", "HE2"]);
    for (rows, component, charged, attached) in [
        (cytosine, "DC", vec![("N3".to_owned(), 1)], Some("N3")),
        (
            imidazolate.as_str(),
            "HIS",
            vec![("ND1".to_owned(), 0), ("NE2".to_owned(), -1)],
            None,
        ),
    ] {
        let result = interpret(&block(rows, ""));
        let nonzero = charged
            .iter()
            .filter(|(_, charge)| *charge != 0)
            .cloned()
            .collect::<Vec<_>>();
        assert_eq!(charges(&result), nonzero, "{component}");
        let mut expected = attached
            .map(
                |heavy| StructureInterpretIssue::HydrogenAttachedByDistance {
                    residue: residue("A", 1),
                    atom_name: "H3".to_owned(),
                    heavy_atom_name: heavy.to_owned(),
                },
            )
            .into_iter()
            .collect::<Vec<_>>();
        expected.push(StructureInterpretIssue::ProtonationFromHydrogens {
            residue: residue("A", 1),
            component_id: component.to_owned(),
            variant: None,
            charged_atoms: charged,
        });
        assert_eq!(issues(&result), expected, "{component}");
    }

    // Missing ring hydrogens select a neutral tautomer, whose ring bonds
    // follow the hydrogen.
    for (removed, variant, nd1_ce1, ce1_ne2) in [
        (
            "HD1",
            "NE2-protonated",
            BondOrder::Double,
            BondOrder::Single,
        ),
        (
            "HE2",
            "ND1-protonated",
            BondOrder::Single,
            BondOrder::Double,
        ),
    ] {
        let tautomer = interpret(&block(&without(HISTIDINE, &[removed]), ""));
        assert_eq!(charges(&tautomer), [], "{variant}");
        assert_eq!(
            issues(&tautomer),
            [StructureInterpretIssue::ProtonationFromHydrogens {
                residue: residue("A", 1),
                component_id: "HIS".to_owned(),
                variant: Some(variant.to_owned()),
                charged_atoms: Vec::new(),
            }]
        );
        assert_eq!(order_between(&tautomer, "ND1", "CE1"), nd1_ce1, "{variant}");
        assert_eq!(order_between(&tautomer, "CE1", "NE2"), ce1_ne2, "{variant}");
    }
}

const UNKNOWN_LIGAND: &str = "
        C C1 XYZ A 1 0.0 0.0 0.0
        C C2 XYZ A 1 1.5 0.0 0.0
        O O1 XYZ A 1 2.0 1.2 0.0";

#[test]
fn residues_without_template_charges_take_the_charges_their_bonds_imply() {
    // wwPDB entry files state component bond orders but no charges. A
    // quaternary nitrogen takes +1 and a tetrahedral boron -1; a carbon
    // with five bonds fits neither and stays neutral.
    let mut table = String::from(
            "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n",
        );
    for (component, center, ligands) in [("XNM", "N", 4), ("XBM", "B", 4), ("XCM", "C", 5)] {
        for ligand in 1..=ligands {
            table.push_str(&format!("{component} {center} C{ligand} sing\n"));
        }
    }
    let result = interpret(&block(
        "N N XNM A 1 0.000 0.000 0.000
             C C1 XNM A 1 0.860 0.860 0.860
             C C2 XNM A 1 -0.860 -0.860 0.860
             C C3 XNM A 1 -0.860 0.860 -0.860
             C C4 XNM A 1 0.860 -0.860 -0.860
             B B XBM A 2 10.000 0.000 0.000
             C C1 XBM A 2 10.900 0.900 0.900
             C C2 XBM A 2 9.100 -0.900 0.900
             C C3 XBM A 2 9.100 0.900 -0.900
             C C4 XBM A 2 10.900 -0.900 -0.900
             C C XCM A 3 20.000 0.000 0.000
             C C1 XCM A 3 21.500 0.000 0.000
             C C2 XCM A 3 18.500 0.000 0.000
             C C3 XCM A 3 20.000 1.500 0.000
             C C4 XCM A 3 20.000 -1.500 0.000
             C C5 XCM A 3 20.000 0.000 1.500",
        &table,
    ));
    assert_eq!(
        charges(&result),
        [("N".to_owned(), 1), ("B".to_owned(), -1)]
    );
    let charged = issues(&result)
        .into_iter()
        .filter(|issue| matches!(issue, StructureInterpretIssue::ChargesFromBonds { .. }))
        .collect::<Vec<_>>();
    assert_eq!(
        charged,
        [
            StructureInterpretIssue::ChargesFromBonds {
                residue: residue("A", 1),
                component_id: "XNM".to_owned(),
                charged_atoms: vec![("N".to_owned(), 1)],
            },
            StructureInterpretIssue::ChargesFromBonds {
                residue: residue("A", 2),
                component_id: "XBM".to_owned(),
                charged_atoms: vec![("B".to_owned(), -1)],
            },
        ]
    );
}

#[test]
fn residues_without_templates_infer_bonds_unless_disabled() {
    let text = block(UNKNOWN_LIGAND, "");
    let inferred = interpret(&text);
    assert_eq!(
        bonds(&inferred),
        [
            ("C1", "C2", BondOrder::Single, BondEvidence::Distance),
            ("C2", "O1", BondOrder::Single, BondEvidence::Distance),
        ]
        .map(|(a, b, order, evidence)| (a.to_owned(), b.to_owned(), order, evidence))
    );
    assert_eq!(
        issues(&inferred),
        [StructureInterpretIssue::ResidueBondsInferred {
            residue: residue("A", 1),
            component_id: "XYZ".to_owned(),
            bonds: 2,
        }]
    );

    let strict = StructureInterpretOptions::default()
        .with_connectivity(ConnectivityOptions::default().with_infer_bonds_from_distance(false));
    let error = crate::mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(strict.clone())
        .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);
    assert_eq!(
        error.location().map(|location| location.line()),
        Some(line_of(&text, " C1 XYZ"))
    );

    // A caller template takes precedence and supplies the bond order.
    let library = ComponentLibrary::from_mmcif(
            &crate::mmcif::parse_str(
                "data_XYZ\n_chem_comp.id XYZ\n_chem_comp.type NON-POLYMER\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\nXYZ C1 C 0\nXYZ C2 C 0\nXYZ O1 O 0\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nXYZ C1 C2 SING\nXYZ C2 O1 DOUB\n",
                Default::default(),
            )
            .unwrap(),
        )
        .unwrap();
    let with_library = strict.clone().with_connectivity(
        strict
            .connectivity()
            .clone()
            .with_components(Some(std::sync::Arc::new(library))),
    );
    let templated = crate::mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(with_library)
        .unwrap();
    assert_eq!(
        bonds(&templated),
        [
            ("C1", "C2", BondOrder::Single, template("XYZ", None)),
            ("C2", "O1", BondOrder::Double, template("XYZ", None)),
        ]
        .map(|(a, b, order, evidence)| (a.to_owned(), b.to_owned(), order, evidence))
    );
    assert!(issues(&templated).is_empty());
}

const DIGLYCINE: &str = "
        N N GLY A 1 0.0 0.0 0.0
        C CA GLY A 1 1.45 0.0 0.0
        C C GLY A 1 2.0 1.4 0.0
        O O GLY A 1 1.3 2.4 0.0
        N N GLY A 2 3.33 1.4 0.0
        C CA GLY A 2 4.0 2.6 0.0
        C C GLY A 2 5.5 2.6 0.0
        O O GLY A 2 6.1 3.6 0.0";

#[test]
fn polymer_links_require_bonding_distance() {
    // Without _entity_poly, both residues' templates identify a peptide.
    let linked = interpret(&records("ATOM", DIGLYCINE, ""));
    assert_eq!(linked.topology().instance_count(), 1);
    assert_eq!(
        bonds(&linked)
            .iter()
            .filter(|(_, _, _, evidence)| *evidence == BondEvidence::PolymerLink)
            .map(|(a, b, _, _)| (a.clone(), b.clone()))
            .collect::<Vec<_>>(),
        [("C".to_owned(), "N".to_owned())]
    );
    assert!(issues(&linked).is_empty());

    // Residue 2 moved 3 A away: adjacent in sequence, but not bonded.
    let shifted = DIGLYCINE
        .lines()
        .map(|row| {
            let fields = row.split_whitespace().collect::<Vec<_>>();
            if fields.get(4) == Some(&"2") {
                let x = fields[5].parse::<f64>().unwrap() + 3.0;
                format!(
                    "{} {} {} {} {} {x} {} {}",
                    fields[0], fields[1], fields[2], fields[3], fields[4], fields[6], fields[7]
                )
            } else {
                row.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let broken = interpret(&records("ATOM", &shifted, ""));
    assert_eq!(broken.topology().instance_count(), 2);
    assert_eq!(broken.topology().bond_count(), 6);
    assert_eq!(
        issues(&broken),
        [StructureInterpretIssue::ChainBreak {
            residue: residue("A", 1),
            next: residue("A", 2),
        }]
    );
}

#[test]
fn untemplated_polymer_residues_link_by_their_neighbours_linkage_or_atom_names() {
    let inferred = |sequence, component: &str| StructureInterpretIssue::ResidueBondsInferred {
        residue: residue("A", sequence),
        component_id: component.to_owned(),
        bonds: 3,
    };
    // Without _entity_poly, a residue without a template defers to its
    // templated neighbour, and two such residues link by atom names.
    for (first, second, expected) in [
        ("GLY", "XGL", vec![inferred(2, "XGL")]),
        ("XGL", "GLY", vec![inferred(1, "XGL")]),
        ("XGL", "YGL", vec![inferred(1, "XGL"), inferred(2, "YGL")]),
    ] {
        let rows = DIGLYCINE
            .replace("GLY A 1", &format!("{first} A 1"))
            .replace("GLY A 2", &format!("{second} A 2"));
        let linked = interpret(&records("ATOM", &rows, ""));
        let polymer_links = bonds(&linked)
            .iter()
            .filter(|(_, _, _, evidence)| *evidence == BondEvidence::PolymerLink)
            .map(|(a, b, _, _)| (a.clone(), b.clone()))
            .collect::<Vec<_>>();
        assert_eq!(polymer_links, [("C".to_owned(), "N".to_owned())]);
        assert_eq!(linked.topology().instance_count(), 1);
        assert_eq!(issues(&linked), expected);
    }
}

#[test]
fn template_mismatches_fail_and_split_residues_are_reported() {
    let unknown_atom = block(
        "N N GLY A 1 0.0 0.0 0.0
             C CA GLY A 1 1.45 0.0 0.0
             C CX GLY A 1 2.0 1.4 0.0",
        "",
    );
    let error = crate::mmcif::parse_str(&unknown_atom, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.location().map(|location| location.line()),
        Some(line_of(&unknown_atom, " CX GLY"))
    );
    assert_eq!(error.message(), "Atom \"CX\" is not in component \"GLY\".");

    // Without its CA, a glycine's amine is cut off from its carbonyl.
    let split = interpret(&block(
        "N N GLY A 1 0.0 0.0 0.0
             C C GLY A 1 2.0 1.4 0.0
             O O GLY A 1 1.3 2.4 0.0",
        "",
    ));
    assert_eq!(split.topology().instance_count(), 2);
    assert_eq!(
        issues(&split),
        [StructureInterpretIssue::ResidueFragmented {
            residue: residue("A", 1),
            fragments: 2,
        }]
    );
}

#[test]
fn stated_internal_bonds_take_precedence_over_templates() {
    // Files that enumerate a residue's bonds, as Kekule's writer does,
    // are taken as stated; the template does not add the others.
    let text = block(
            &DIGLYCINE.lines().take(5).collect::<Vec<_>>().join("\n"),
            "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_label_seq_id\nc1 covale A N 1 A CA 1\n",
        );
    let stated = interpret(&text);
    assert_eq!(
        bonds(&stated),
        [(
            "N".to_owned(),
            "CA".to_owned(),
            BondOrder::Single,
            BondEvidence::Connection {
                connection_id: Some("c1".to_owned()),
                connection_type: "covale".to_owned(),
            },
        )]
    );
}
