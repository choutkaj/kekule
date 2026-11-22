//! Charges of residues whose component has no template, from the formal
//! charges of the complete Chemical Component Dictionary.
use kekule::components::DictionaryCharges;
use kekule::macromolecular::{StructureInterpretIssue, StructureInterpretation};
use kekule::{mmcif, pdb};

/// The atoms of TRS (tris, the buffer) with the dictionary's ideal
/// coordinates: name, element, and position.
const TRIS: [(&str, &str, [f64; 3]); 20] = [
    ("C", "C", [0.013, -0.296, 0.014]),
    ("C1", "C", [-0.042, 0.986, 0.848]),
    ("C2", "C", [1.263, -0.276, -0.869]),
    ("C3", "C", [-1.233, -0.384, -0.869]),
    ("N", "N", [0.064, -1.460, 0.908]),
    ("O1", "O", [-0.091, 2.118, -0.022]),
    ("O2", "O", [2.429, -0.311, -0.044]),
    ("O3", "O", [-2.400, -0.403, -0.044]),
    ("H11", "H", [-0.932, 0.972, 1.477]),
    ("H12", "H", [0.846, 1.049, 1.477]),
    ("H21", "H", [1.257, -1.145, -1.527]),
    ("H22", "H", [1.269, 0.633, -1.469]),
    ("H31", "H", [-1.271, 0.480, -1.532]),
    ("H32", "H", [-1.194, -1.297, -1.463]),
    ("HN1", "H", [-0.760, -1.473, 1.490]),
    ("HN2", "H", [0.886, -1.402, 1.490]),
    ("HN3", "H", [0.100, -2.305, 0.358]),
    ("HO1", "H", [-0.125, 2.903, 0.542]),
    ("HO2", "H", [3.192, -0.297, -0.638]),
    ("HO3", "H", [-3.161, -0.459, -0.638]),
];

/// The heavy atoms of BET (betaine) with the dictionary's ideal coordinates.
const BETAINE: [(&str, &str, [f64; 3]); 8] = [
    ("N", "N", [-1.076, -0.083, 0.001]),
    ("CA", "C", [0.112, -0.940, -0.113]),
    ("C", "C", [1.354, -0.091, -0.030]),
    ("O", "O", [1.261, 1.114, 0.004]),
    ("OXT", "O", [2.563, -0.672, 0.006]),
    ("C1", "C", [-1.112, 0.857, -1.127]),
    ("C2", "C", [-1.019, 0.668, 1.262]),
    ("C3", "C", [-2.285, -0.917, -0.018]),
];

/// `HETATM` records of residue `sequence`, without stated charges.
fn hetatm(residue: &str, sequence: usize, atoms: &[(&str, &str, [f64; 3])], offset: f64) -> String {
    atoms
        .iter()
        .enumerate()
        .map(|(index, (name, element, [x, y, z]))| {
            let name = if name.len() < 4 && element.len() == 1 {
                format!(" {name:<3}")
            } else {
                format!("{name:<4}")
            };
            format!(
                "HETATM{:>5} {name} {residue:>3} A{sequence:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00          {element:>2}\n",
                sequence * 100 + index + 1,
                x + offset,
                y,
                z
            )
        })
        .collect()
}

fn read_pdb(text: &str) -> StructureInterpretation {
    pdb::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// The formal charge of each named atom of `component`.
fn charge(result: &StructureInterpretation, component: &str, name: &str) -> i8 {
    let atom = result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .find(|atom| atom.component_id() == component && atom.atom_name() == name)
        .unwrap_or_else(|| panic!("{component} {name}"))
        .atom();
    result
        .topology()
        .instance(atom.instance())
        .unwrap()
        .molecule()
        .atom(atom.atom())
        .unwrap()
        .formal_charge()
}

fn dictionary_issues(result: &StructureInterpretation) -> Vec<StructureInterpretIssue> {
    result
        .report()
        .issues()
        .iter()
        .filter(|issue| {
            matches!(
                issue,
                StructureInterpretIssue::ChargesFromDictionary { .. }
                    | StructureInterpretIssue::DictionaryAtomsUnmatched { .. }
                    | StructureInterpretIssue::ChargesFromBonds { .. }
                    | StructureInterpretIssue::TemplateMissing { .. }
                    | StructureInterpretIssue::TemplateChargeUnknown { .. }
                    | StructureInterpretIssue::ProtonationAssumed { .. }
            )
        })
        .cloned()
        .collect()
}

/// The component, charged atoms, and unheld charges of a
/// `ChargesFromDictionary` issue.
type Applied = (String, Vec<(String, i8)>, Vec<(String, i8)>);

fn applied(result: &StructureInterpretation) -> Vec<Applied> {
    dictionary_issues(result)
        .into_iter()
        .filter_map(|issue| match issue {
            StructureInterpretIssue::ChargesFromDictionary {
                component_id,
                charged_atoms,
                unheld,
                ..
            } => Some((component_id, charged_atoms, unheld)),
            _ => None,
        })
        .collect()
}

fn named(atoms: &[(&str, i8)]) -> Vec<(String, i8)> {
    atoms
        .iter()
        .map(|(name, charge)| ((*name).to_owned(), *charge))
        .collect()
}

#[test]
fn unbundled_charged_components_without_hydrogens_take_their_dictionary_charges() {
    // Tris is the ammonium and platinum(II) is +2. Gold chloride states its
    // gold's charge as unknown, and its chloride's -1 does not fit the bond
    // to gold that distances give. None of them is bundled.
    let heavy = TRIS
        .iter()
        .copied()
        .filter(|(_, element, _)| *element != "H")
        .collect::<Vec<_>>();
    let text = format!(
        "{}{}{}END\n",
        hetatm("TRS", 1, &heavy, 0.0),
        hetatm("PT", 2, &[("PT", "PT", [10.0, 0.0, 0.0])], 0.0),
        hetatm(
            "2ES",
            3,
            &[
                ("AU", "AU", [25.029, 4.847, 33.611]),
                ("CL", "CL", [23.479, 3.320, 34.314])
            ],
            0.0
        ),
    );
    let result = read_pdb(&text);
    assert_eq!(charge(&result, "TRS", "N"), 1);
    assert_eq!(charge(&result, "PT", "PT"), 2);
    assert_eq!(charge(&result, "2ES", "CL"), 0);
    assert_eq!(charge(&result, "2ES", "AU"), 0);
    assert_eq!(
        applied(&result),
        [
            ("TRS".to_owned(), named(&[("N", 1)]), Vec::new()),
            ("PT".to_owned(), named(&[("PT", 2)]), Vec::new()),
            ("2ES".to_owned(), Vec::new(), named(&[("CL", -1)])),
        ]
    );
    let issues = dictionary_issues(&result);
    assert!(issues.iter().any(|issue| matches!(
        issue,
        StructureInterpretIssue::TemplateChargeUnknown { component_id, atoms, .. }
            if component_id == "2ES" && atoms == &["AU"]
    )));
    assert!(!issues
        .iter()
        .any(|issue| matches!(issue, StructureInterpretIssue::TemplateMissing { .. })));
    for issue in &issues {
        assert!(issue.to_string().ends_with('.'), "{issue}");
    }
    let mut perceived = result.structure().clone();
    perceived.perceive().unwrap();
    let nitrogen = perceived
        .topology()
        .atoms()
        .find(|atom| atom.element().atomic_number() == 7)
        .unwrap();
    assert_eq!(nitrogen.formal_charge(), 1);
    assert_eq!(nitrogen.total_hydrogens(), Some(3));
}

#[test]
fn hydrogens_decide_the_charges_of_the_atoms_they_bond() {
    // Tris modelled as the neutral amine keeps it. With hydrogens on its
    // hydroxyls only, its nitrogen has no hydrogen evidence and takes the
    // dictionary's +1.
    let amine = TRIS
        .iter()
        .copied()
        .filter(|(name, _, _)| *name != "HN3")
        .collect::<Vec<_>>();
    let hydroxyls = TRIS
        .iter()
        .copied()
        .filter(|(name, _, _)| !name.starts_with("HN"))
        .collect::<Vec<_>>();
    let text = format!(
        "{}{}END\n",
        hetatm("TRS", 1, &amine, 0.0),
        hetatm("TRS", 2, &hydroxyls, 20.0)
    );
    let result = read_pdb(&text);
    assert_eq!(
        applied(&result),
        [("TRS".to_owned(), named(&[("N", 1)]), Vec::new())]
    );
    let mut perceived = result.structure().clone();
    perceived.perceive().unwrap();
    let nitrogens = perceived
        .topology()
        .atoms()
        .filter(|atom| atom.element().atomic_number() == 7)
        .map(|atom| (atom.formal_charge(), atom.total_hydrogens()))
        .collect::<Vec<_>>();
    assert_eq!(nitrogens, [(0, Some(2)), (1, Some(3))]);
}

#[test]
fn atom_names_the_dictionary_lacks_are_reported_and_take_no_charges() {
    let renamed = TRIS
        .iter()
        .copied()
        .filter(|(_, element, _)| *element != "H")
        .map(|(name, element, position)| (if name == "N" { "N1" } else { name }, element, position))
        .collect::<Vec<_>>();
    let result = read_pdb(&format!("{}END\n", hetatm("TRS", 1, &renamed, 0.0)));
    assert_eq!(charge(&result, "TRS", "N1"), 0);
    let issues = dictionary_issues(&result);
    assert!(matches!(
        &issues[..],
        [StructureInterpretIssue::DictionaryAtomsUnmatched { component_id, atoms, .. }]
            if component_id == "TRS" && atoms == &["N1"]
    ));
    assert_eq!(
        issues[0].to_string(),
        "The Chemical Component Dictionary entry of component \"TRS\" does not name atom \
         \"N1\" of author residue \"1\" of chain \"A\" in model \"1\". It gives the residue no \
         charges."
    );
}

#[test]
fn stated_charges_outrank_the_dictionary() {
    // Writers state every charge, zero included, so a written neutral amine
    // reads back as one.
    let text = hetatm(
        "TRS",
        1,
        &TRIS
            .iter()
            .copied()
            .filter(|(_, element, _)| *element != "H")
            .collect::<Vec<_>>(),
        0.0,
    )
    .lines()
    .map(|line| format!("{line} 0\n"))
    .collect::<String>();
    let result = read_pdb(&format!("{text}END\n"));
    assert_eq!(charge(&result, "TRS", "N"), 0);
    assert!(dictionary_issues(&result).is_empty());
}

#[test]
fn atoms_the_blocks_component_table_omits_take_their_dictionary_charges() {
    // The block's `_chem_comp_atom` states the charge of one tris carbon
    // only; the nitrogen, which it omits, still takes the dictionary's +1.
    let mut text = String::from(
        "data_test\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
         _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\nTRS C C 0\n\
         loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
         _atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    for (index, (name, element, [x, y, z])) in TRIS
        .iter()
        .filter(|(_, element, _)| *element != "H")
        .enumerate()
    {
        text.push_str(&format!(
            "HETATM {} {element} {name} TRS A . {x} {y} {z}\n",
            index + 1
        ));
    }
    let result = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(charge(&result, "TRS", "N"), 1);
    assert_eq!(charge(&result, "TRS", "C"), 0);
    assert_eq!(
        applied(&result),
        vec![("TRS".to_owned(), named(&[("N", 1)]), Vec::new())]
    );
}

/// A block of betaine's heavy atoms, with its component bond table when
/// `orders` holds.
fn betaine_block(orders: bool) -> String {
    let mut text = String::from("data_test\n");
    if orders {
        text.push_str(
            "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n\
             _chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n\
             BET N CA sing\nBET N C1 sing\nBET N C2 sing\nBET N C3 sing\n\
             BET CA C sing\nBET C O doub\nBET C OXT sing\n",
        );
    }
    text.push_str(
        "loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
         _atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    for (index, (name, element, [x, y, z])) in BETAINE.iter().enumerate() {
        text.push_str(&format!(
            "HETATM {} {element} {name} BET A . {x} {y} {z}\n",
            index + 1
        ));
    }
    text
}

#[test]
fn neutral_acid_hydroxyls_of_dictionary_charged_residues_take_their_ph_7_state() {
    // The dictionary states betaine as a cation with a carboxylic acid. With
    // stated orders the acid is deprotonated, as in a templated residue;
    // without them a hydroxyl and a carbonyl look alike, and it stays.
    let read = |orders| {
        mmcif::parse_str(&betaine_block(orders), Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap()
    };
    let stated = read(true);
    assert_eq!(charge(&stated, "BET", "N"), 1);
    assert_eq!(charge(&stated, "BET", "OXT"), -1);
    assert_eq!(charge(&stated, "BET", "O"), 0);
    assert!(dictionary_issues(&stated).iter().any(|issue| matches!(
        issue,
        StructureInterpretIssue::ProtonationAssumed { components }
            if components == &[("BET".to_owned(), 1)]
    )));
    let mut perceived = stated.structure().clone();
    perceived.perceive().unwrap();
    let inferred = read(false);
    assert_eq!(charge(&inferred, "BET", "N"), 1);
    assert_eq!(charge(&inferred, "BET", "OXT"), 0);
    assert!(!dictionary_issues(&inferred)
        .iter()
        .any(|issue| matches!(issue, StructureInterpretIssue::ProtonationAssumed { .. })));
}

#[test]
fn the_dictionary_charge_table_parses_and_states_its_release() {
    let release = DictionaryCharges::release();
    let date = release.as_bytes();
    assert!(
        date.len() == 10
            && date[4] == b'-'
            && date[7] == b'-'
            && release
                .bytes()
                .enumerate()
                .all(|(index, byte)| index == 4 || index == 7 || byte.is_ascii_digit()),
        "{release}"
    );
    let entries = DictionaryCharges::all().collect::<Vec<_>>();
    assert!(entries.len() > 3_000, "{}", entries.len());
    assert!(entries.windows(2).all(|pair| pair[0].id() < pair[1].id()));
    let mut unknown = 0;
    for entry in &entries {
        let atoms = entry.atoms().collect::<Vec<_>>();
        assert!(
            atoms.iter().any(|atom| atom.formal_charge() != Some(0)),
            "{}",
            entry.id()
        );
        let mut names = atoms.iter().map(|atom| atom.name()).collect::<Vec<_>>();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), atoms.len(), "{}", entry.id());
        unknown += usize::from(atoms.iter().any(|atom| atom.formal_charge().is_none()));
        assert_eq!(DictionaryCharges::get(entry.id()), Some(*entry));
    }
    assert!(unknown > 500, "{unknown}");
    let tris = DictionaryCharges::get("trs").unwrap();
    assert_eq!(tris.atom("N").unwrap().formal_charge(), Some(1));
    assert_eq!(tris.atom("N").unwrap().element().atomic_number(), 7);
    assert_eq!(tris.atom("HN1"), None);
    let gold = DictionaryCharges::get("2ES").unwrap().atom("AU").unwrap();
    assert_eq!(
        (gold.element().atomic_number(), gold.formal_charge()),
        (79, None)
    );
}
