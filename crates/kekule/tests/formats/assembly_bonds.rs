//! Bonds between the copies of an assembly: those of the asymmetric unit
//! between generators, links a record states to a symmetry mate, and
//! disulfides between copies.
use kekule::macromolecular::StructureInterpretation;
use kekule::mmcif;
use kekule::stereo::{
    AxisOrientation, AxisStereo, StereoCarrier, StereoElement, StereoElementKind,
};
use kekule::structure::Structure;

/// A cysteine without `HG`, as the dictionary's ideal `CYS` places it.
const CYSTEINE: &[(&str, &str, [f64; 3])] = &[
    ("N", "N", [1.585, 0.483, -0.081]),
    ("CA", "C", [0.141, 0.450, 0.186]),
    ("C", "C", [-0.095, 0.006, 1.606]),
    ("O", "O", [0.685, -0.742, 2.143]),
    ("CB", "C", [-0.533, -0.530, -0.774]),
    ("SG", "S", [-0.247, 0.004, -2.484]),
    ("OXT", "O", [-1.174, 0.443, 2.275]),
];

/// A two-fold axis along z that puts the copy's `SG` 2.04 Å from the
/// original's, with CB-SG-SG angles of 109°.
const TWOFOLD: &str = "2 -1 0 0 0.468 0 -1 0 1.806 0 0 1 0";

/// Chains of [`CYSTEINE`], each as placed or turned by [`TWOFOLD`].
fn block(chains: &[(&str, bool)], connections: &str, generators: &str) -> String {
    let mut text = String::from(
        "data_dimer\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
         _atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
         _atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.auth_seq_id\n\
         _atom_site.auth_asym_id\n_atom_site.pdbx_PDB_model_num\n",
    );
    let mut serial = 0;
    for &(chain, turned) in chains {
        for &(name, element, [x, y, z]) in CYSTEINE {
            serial += 1;
            let (x, y) = if turned {
                (0.468 - x, 1.806 - y)
            } else {
                (x, y)
            };
            text.push_str(&format!(
                "ATOM {serial} {element} {name} CYS {chain} 1 1 {x:.3} {y:.3} {z:.3} 1 {chain} 1\n"
            ));
        }
    }
    text.push_str(
        "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n\
         _struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_seq_id\n\
         _struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_symmetry\n\
         _struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_comp_id\n\
         _struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n\
         _struct_conn.ptnr2_symmetry\n",
    );
    text.push_str(connections);
    text.push_str(&format!(
        "#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n\
         _pdbx_struct_assembly_gen.asym_id_list\n{generators}\n#\nloop_\n_pdbx_struct_oper_list.id\n\
         _pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n\
         _pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n\
         _pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n\
         _pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n\
         _pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n\
         _pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n\
         1 1 0 0 0 0 1 0 0 0 0 1 0\n{TWOFOLD}\n#\n"
    ));
    text
}

fn read(text: &str) -> StructureInterpretation {
    mmcif::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// The `SG`-`SG` bonds of `structure`, by the chain labels of their atoms.
fn disulfides(structure: &Structure) -> Vec<[String; 2]> {
    let topology = structure.topology();
    topology
        .bonds()
        .filter_map(|bond| {
            let [a, b] = bond.atoms().map(|atom| atom.atom_site().unwrap());
            (a.metadata().label_atom_id() == Some("SG")
                && b.metadata().label_atom_id() == Some("SG"))
            .then(|| [a, b].map(|site| site.residue().chain().label_id().to_owned()))
        })
        .collect()
}

#[test]
fn a_disulfide_to_a_symmetry_mate_bonds_the_copies() {
    let text = block(
        &[("A", false)],
        "disulf1 disulf A CYS 1 SG 1_555 A CYS 1 SG 2_555\n",
        "1 '1,2' A",
    );
    let interpretation = read(&text);
    let assembly = &interpretation.report().assemblies()[0];
    let built = assembly.to_structure(interpretation.structure()).unwrap();
    assert_eq!(built.report().created_bonds().len(), 1);
    assert!(built.report().cut_bonds().is_empty());
    let built = built.into_structure();
    assert_eq!(disulfides(&built), [["A".to_owned(), "A-2".to_owned()]]);
    // Neither sulfur keeps a hydrogen.
    let mut perceived = built.clone();
    perceived.perceive().unwrap();
    for atom in perceived.topology().atoms() {
        if atom.atom().element().symbol() == "S" {
            assert_eq!(atom.total_hydrogens(), Some(0));
        }
    }
}

#[test]
fn bonds_between_generators_under_one_operation_stay() {
    // Chains A and B of the asymmetric unit form a disulfide, and the
    // assembly lists them in separate generators under the same operator.
    let a_to_b = "disulf1 disulf A CYS 1 SG 1_555 B CYS 1 SG 1_555\n";
    let interpretation = read(&block(&[("A", false), ("B", true)], a_to_b, "1 1 A\n1 1 B"));
    assert_eq!(disulfides(interpretation.structure()).len(), 1);
    let assembly = &interpretation.report().assemblies()[0];
    assert_eq!(assembly.generators().len(), 2);
    let built = assembly.to_structure(interpretation.structure()).unwrap();
    assert!(built.report().cut_bonds().is_empty());
    assert!(built.report().created_bonds().is_empty());
    assert_eq!(
        disulfides(built.structure()),
        [["A".to_owned(), "B".to_owned()]]
    );
}

#[test]
fn a_bond_to_a_chain_no_copy_holds_is_reported() {
    // Chain B is bonded to A, but the assembly copies A alone.
    let a_to_b = "disulf1 disulf A CYS 1 SG 1_555 B CYS 1 SG 1_555\n";
    let interpretation = read(&block(&[("A", false), ("B", true)], a_to_b, "1 1 A"));
    let assembly = &interpretation.report().assemblies()[0];
    let built = assembly.to_structure(interpretation.structure()).unwrap();
    let cut = built.report().cut_bonds();
    assert_eq!(cut.len(), 1);
    assert_eq!(cut[0].operation(), "1");
    assert!(disulfides(built.structure()).is_empty());
}

/// The disulfides an assembly of chain A and its two-fold copy forms with
/// an ion of `element`, component `component`, on the axis 2.35 A from both
/// sulfurs.
fn disulfides_beside(element: &str, component: &str) -> usize {
    let text = block(
        &[("A", false)],
        &format!(
            "metalc1 metalc A CYS 1 SG 1_555 B {component} . {component} 1_555
"
        ),
        "1 '1,2' A
1 1 B",
    )
    .replace(
        "loop_
_struct_conn.id",
        &format!(
            "HETATM 8 {element} {component} {component} B 2 . 0.234 0.903 -4.600 1 B 1
loop_
_struct_conn.id"
        ),
    );
    let interpretation = read(&text);
    let assembly = &interpretation.report().assemblies()[0];
    let built = assembly.to_structure(interpretation.structure()).unwrap();
    assert_eq!(
        built.report().created_bonds().len(),
        disulfides(built.structure()).len()
    );
    built.report().created_bonds().len()
}

#[test]
fn a_sulfur_beside_zinc_forms_no_disulfide() {
    // The sulfurs coordinate the zinc rather than bond each other.
    assert_eq!(disulfides_beside("Zn", "ZN"), 0);
}

#[test]
fn a_sulfur_beside_sodium_still_forms_a_disulfide() {
    assert_eq!(disulfides_beside("Na", "NA"), 1);
}

#[test]
fn stereo_across_generators_is_reported_lost() {
    // The disulfide's helicity, stated as an axis about its S-S bond, names
    // a carbon of each chain, and the assembly copies the chains apart.
    let a_to_b = "disulf1 disulf A CYS 1 SG 1_555 B CYS 1 SG 1_555\n";
    let interpretation = read(&block(&[("A", false), ("B", true)], a_to_b, "1 1 A\n1 1 B"));
    let structure = interpretation.structure();
    let topology = structure.topology();
    let named = |chain: &str, name: &str| {
        topology
            .atom_sites()
            .find(|site| {
                site.residue().chain().label_id() == chain
                    && site.metadata().label_atom_id() == Some(name)
            })
            .unwrap()
            .atom()
            .id()
    };
    let [sulfur_a, sulfur_b, carbon_a, carbon_b] =
        [("A", "SG"), ("B", "SG"), ("A", "CB"), ("B", "CB")]
            .map(|(chain, name)| named(chain, name));
    let instance = sulfur_a.instance();
    let axis = topology
        .instance(instance)
        .unwrap()
        .molecule()
        .bond_between(sulfur_a.atom(), sulfur_b.atom())
        .unwrap();
    let helicity = StereoElement::new(StereoElementKind::Axis(AxisStereo {
        axis,
        carriers: [
            StereoCarrier::Atom(carbon_a.atom()),
            StereoCarrier::Atom(carbon_b.atom()),
        ],
        orientation: Some(AxisOrientation::Clockwise),
    }));
    let mut editor = structure.to_editor();
    editor
        .replace_source_instance_stereo(instance, &[helicity])
        .unwrap();
    let stated = editor.finish().unwrap();
    let assembly = &interpretation.report().assemblies()[0];
    let built = assembly.to_structure(&stated).unwrap();
    assert!(built.report().cut_bonds().is_empty());
    // Each chain's copy loses it.
    let lost = built.report().lost_stereo();
    assert_eq!(lost.len(), 2);
    let mut axis_atoms = vec![sulfur_a, sulfur_b];
    axis_atoms.sort_unstable();
    for loss in lost {
        let mut atoms = loss.atoms().to_vec();
        atoms.sort_unstable();
        assert_eq!(atoms, axis_atoms);
        assert_eq!(loss.operation(), "1");
    }
}
