//! Charges and hydrogens of mmCIF residues whose component tables give no
//! charges: onium centres, charge pairs, and metals.
use kekule::macromolecular::{
    SourceResidueId, SourceResiduePosition, StructureInterpretIssue, StructureInterpretation,
};
use kekule::mmcif;

/// A minimal block of `HETATM` records whose `rows` are `type name
/// component asym seq x y z`, after the categories in `extra`.
fn block(rows: &str, extra: &str) -> String {
    let mut text = format!(
        "data_test\n{extra}\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
         _atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
         _atom_site.Cartn_y\n_atom_site.Cartn_z\n"
    );
    for (index, row) in rows
        .lines()
        .map(str::trim)
        .filter(|row| !row.is_empty())
        .enumerate()
    {
        text.push_str(&format!("HETATM {} {row}\n", index + 1));
    }
    text
}

fn interpret(text: &str) -> StructureInterpretation {
    mmcif::parse_str(text, Default::default())
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
                .formal_charge();
            (charge != 0).then(|| (atom.atom_name().to_owned(), charge))
        })
        .collect()
}

/// Residue `sequence` of chain A by label.
fn label(sequence: i32) -> SourceResidueId {
    SourceResidueId::from_parts("1", "A", SourceResiduePosition::Label(sequence), None)
}

#[test]
fn onium_centres_without_template_charges_take_a_charge_not_a_hydrogen() {
    // A sulfonium, phosphonium, or iodonium ligand whose component table
    // states no charges is charged +1; it does not gain a hydrogen.
    let mut table = String::from(
            "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n",
        );
    for (component, center, ligands) in [("XSM", "S", 3), ("XPM", "P", 4), ("XIM", "I", 2)] {
        for ligand in 1..=ligands {
            table.push_str(&format!("{component} {center} C{ligand} sing\n"));
        }
    }
    let result = interpret(&block(
        "S S XSM A 1 0.000 0.000 0.000
             C C1 XSM A 1 1.000 1.000 0.800
             C C2 XSM A 1 -1.000 1.000 -0.800
             C C3 XSM A 1 0.000 -1.400 0.000
             P P XPM A 2 10.000 0.000 0.000
             C C1 XPM A 2 10.950 0.950 0.950
             C C2 XPM A 2 9.050 -0.950 0.950
             C C3 XPM A 2 9.050 0.950 -0.950
             C C4 XPM A 2 10.950 -0.950 -0.950
             I I XIM A 3 20.000 0.000 0.000
             C C1 XIM A 3 22.100 0.000 0.000
             C C2 XIM A 3 20.000 2.100 0.000",
        &table,
    ));
    assert_eq!(
        charges(&result),
        [
            ("S".to_owned(), 1),
            ("P".to_owned(), 1),
            ("I".to_owned(), 1)
        ]
    );
    let mut perceived = result.structure().clone();
    perceived.perceive().unwrap();
    for atom in perceived.topology().atoms() {
        if atom.element().atomic_number() != 6 {
            assert_eq!(atom.total_hydrogens(), Some(0), "{}", atom.element());
        }
    }
    let charged = issues(&result)
        .into_iter()
        .filter_map(|issue| match issue {
            StructureInterpretIssue::ChargesFromBonds { charged_atoms, .. } => Some(charged_atoms),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(
        charged,
        [
            vec![("S".to_owned(), 1)],
            vec![("P".to_owned(), 1)],
            vec![("I".to_owned(), 1)],
        ]
    );
}

#[test]
fn onium_oxides_and_boranes_without_template_charges_separate_charges() {
    // DMSO drawn with a single S-O bond, trimethylphosphine-borane,
    // trimethylamine-borane, and trimethylamine N-oxide: the oxygen and
    // boron take -1 with the centre's +1, not hydrogens.
    let mut table = String::from(
            "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n",
        );
    for (component, center, ligands) in [
        ("XSO", "S", ["C1", "C2", "O"].as_slice()),
        ("XPB", "P", &["C1", "C2", "C3", "B"]),
        ("XNB", "N", &["C1", "C2", "C3", "B"]),
        ("XNO", "N", &["C1", "C2", "C3", "O"]),
    ] {
        for ligand in ligands {
            table.push_str(&format!("{component} {center} {ligand} sing\n"));
        }
    }
    let result = interpret(&block(
        "S S XSO A 1 0.000 0.000 0.000
             C C1 XSO A 1 1.800 0.000 0.000
             C C2 XSO A 1 -0.600 1.700 0.000
             O O XSO A 1 -0.500 -0.600 1.300
             P P XPB A 2 10.000 0.000 0.000
             C C1 XPB A 2 11.050 1.050 1.050
             C C2 XPB A 2 8.950 -1.050 1.050
             C C3 XPB A 2 8.950 1.050 -1.050
             B B XPB A 2 11.100 -1.100 -1.100
             N N XNB A 3 20.000 0.000 0.000
             C C1 XNB A 3 20.850 0.850 0.850
             C C2 XNB A 3 19.150 -0.850 0.850
             C C3 XNB A 3 19.150 0.850 -0.850
             B B XNB A 3 20.920 -0.920 -0.920
             N N XNO A 4 30.000 0.000 0.000
             C C1 XNO A 4 30.850 0.850 0.850
             C C2 XNO A 4 29.150 -0.850 0.850
             C C3 XNO A 4 29.150 0.850 -0.850
             O O XNO A 4 30.800 -0.800 -0.800",
        &table,
    ));
    let separated = [("S", "O"), ("P", "B"), ("N", "B"), ("N", "O")]
        .map(|(centre, partner)| vec![(centre.to_owned(), 1), (partner.to_owned(), -1)]);
    assert_eq!(charges(&result), separated.concat());
    let mut perceived = result.structure().clone();
    perceived.perceive().unwrap();
    let hydrogens = perceived
        .topology()
        .atoms()
        .filter(|atom| atom.element().atomic_number() != 6)
        .map(|atom| (atom.element().symbol(), atom.total_hydrogens().unwrap()))
        .collect::<Vec<_>>();
    assert_eq!(
        hydrogens,
        [
            ("S", 0),
            ("O", 0),
            ("P", 0),
            ("B", 3),
            ("N", 0),
            ("B", 3),
            ("N", 0),
            ("O", 0)
        ]
    );
    let charged = issues(&result)
        .into_iter()
        .filter_map(|issue| match issue {
            StructureInterpretIssue::ChargesFromBonds { charged_atoms, .. } => Some(charged_atoms),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(charged, separated);
}

#[test]
fn metals_without_templates_report_valences_their_element_does_not_allow() {
    // Lead is not bundled; the unbonded atom stays neutral, without hydrides.
    // Tin with three carbons and aluminium with two, bonded by distance, are
    // tributylstannane and DIBAL-H without their hydrides.
    let result = interpret(&block(
        "PB PB PB A 1 0.000 0.000 0.000
             SN SN XSN A 2 10.000 0.000 0.000
             C C1 XSN A 2 12.140 0.000 0.000
             C C2 XSN A 2 10.000 2.140 0.000
             C C3 XSN A 2 10.000 0.000 2.140
             AL AL XAL A 3 20.000 0.000 0.000
             C C1 XAL A 3 21.960 0.000 0.000
             C C2 XAL A 3 20.000 1.960 0.000",
        "",
    ));
    assert_eq!(charges(&result), []);
    let mut perceived = result.structure().clone();
    perceived.perceive().unwrap();
    for atom in perceived.topology().atoms() {
        if atom.element().atomic_number() != 6 {
            assert_eq!(atom.total_hydrogens(), Some(0), "{}", atom.element());
        }
    }
    let incomplete = issues(&result)
        .into_iter()
        .filter(|issue| !matches!(issue, StructureInterpretIssue::ResidueBondsInferred { .. }))
        .collect::<Vec<_>>();
    assert_eq!(
        incomplete,
        [("PB", "PB", 1, 0), ("XSN", "SN", 2, 3), ("XAL", "AL", 3, 2)].map(
            |(component, atom, sequence, valence)| {
                StructureInterpretIssue::MetalValenceIncomplete {
                    residue: label(sequence),
                    component_id: component.to_owned(),
                    atom_name: atom.to_owned(),
                    valence,
                }
            }
        )
    );
}
