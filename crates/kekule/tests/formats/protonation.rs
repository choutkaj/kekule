//! Charges and protonation of templated residues: hydrogens decide where the
//! model has them, the usual pH 7 states where it has none, and every
//! difference from the template or the source is reported.

use kekule::macromolecular::{
    BondEvidence, SourceResidueId, SourceResiduePosition, StructureInterpretIssue as Issue,
    StructureInterpretation,
};
use kekule::mmcif;
use kekule::molecule::BondOrder;

/// One atom-site row.
#[derive(Clone)]
struct Site {
    element: String,
    name: String,
    component: String,
    asym: String,
    sequence: i32,
    point: [f64; 3],
    occupancy: Option<f64>,
    charge: Option<i8>,
}

impl Site {
    fn occupancy(mut self, occupancy: f64) -> Self {
        self.occupancy = Some(occupancy);
        self
    }
}

/// Sets the stated charge of the named site.
fn state_charge(sites: &mut [Site], name: &str, charge: i8) {
    sites
        .iter_mut()
        .find(|site| site.name == name)
        .unwrap()
        .charge = Some(charge);
}

/// The atoms of a bundled template at its ideal coordinates, moved by
/// `offset`, as residue `sequence` of chain `asym`.
fn residue_sites(component: &str, asym: &str, sequence: i32, offset: [f64; 3]) -> Vec<Site> {
    let path = format!(
        "{}/data/components/{component}.cif",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).unwrap();
    let ideal = text
        .lines()
        .filter(|line| line.starts_with("_chem_comp_atom."))
        .position(|line| line.trim_end() == "_chem_comp_atom.pdbx_model_Cartn_x_ideal")
        .unwrap();
    text.lines()
        .filter_map(|line| {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() < ideal + 3 || fields[0] != component {
                return None;
            }
            let coordinate = |index: usize| fields[index].parse::<f64>().ok();
            let point = [
                coordinate(ideal)? + offset[0],
                coordinate(ideal + 1)? + offset[1],
                coordinate(ideal + 2)? + offset[2],
            ];
            Some(Site {
                element: fields[3].to_owned(),
                name: fields[1].trim_matches('"').to_owned(),
                component: component.to_owned(),
                asym: asym.to_owned(),
                sequence,
                point,
                occupancy: None,
                charge: None,
            })
        })
        .collect()
}

fn sites(component: &str) -> Vec<Site> {
    residue_sites(component, "A", 1, [0.0; 3])
}

fn without(sites: Vec<Site>, names: &[&str]) -> Vec<Site> {
    sites
        .into_iter()
        .filter(|site| !names.contains(&site.name.as_str()))
        .collect()
}

fn heavy(sites: Vec<Site>) -> Vec<Site> {
    sites
        .into_iter()
        .filter(|site| site.element != "H")
        .collect()
}

fn site<'a>(sites: &'a [Site], name: &str) -> &'a Site {
    sites.iter().find(|site| site.name == name).unwrap()
}

fn point(sites: &[Site], name: &str) -> [f64; 3] {
    site(sites, name).point
}

/// A hydrogen named `name` on `anchor`, 0.97 Å from it away from `away`.
fn hydrogen(sites: &[Site], name: &str, anchor: &str, away: &str) -> Site {
    let (a, b) = (point(sites, anchor), point(sites, away));
    let direction = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let length = direction
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    let anchor_site = site(sites, anchor);
    Site {
        element: "H".to_owned(),
        name: name.to_owned(),
        point: [0, 1, 2].map(|axis| a[axis] + 0.97 * direction[axis] / length),
        ..anchor_site.clone()
    }
}

fn block(sites: &[Site], extra: &str) -> String {
    let mut text = format!(
        "data_test\n{extra}\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
         _atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
         _atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n\
         _atom_site.pdbx_formal_charge\n"
    );
    for (index, site) in sites.iter().enumerate() {
        let optional = |value: Option<String>| value.unwrap_or_else(|| "?".to_owned());
        text.push_str(&format!(
            "ATOM {} {} \"{}\" {} {} {} {:.3} {:.3} {:.3} {} {}\n",
            index + 1,
            site.element,
            site.name,
            site.component,
            site.asym,
            // A ligand has no sequence position.
            if site.sequence == 0 {
                ".".to_owned()
            } else {
                site.sequence.to_string()
            },
            site.point[0],
            site.point[1],
            site.point[2],
            optional(site.occupancy.map(|value| format!("{value:.2}"))),
            optional(site.charge.map(|value| value.to_string())),
        ));
    }
    text
}

fn read(sites: &[Site], extra: &str) -> StructureInterpretation {
    mmcif::parse_str(&block(sites, extra), Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// Issues other than those these minimal blocks always raise.
fn issues(result: &StructureInterpretation) -> Vec<Issue> {
    result
        .report()
        .issues()
        .iter()
        .filter(|issue| {
            !matches!(
                issue,
                Issue::EntityTypeInferred { .. } | Issue::ConnectivityCandidatesInferred { .. }
            )
        })
        .cloned()
        .collect()
}

fn has_issue(result: &StructureInterpretation, matches: impl Fn(&Issue) -> bool) -> bool {
    issues(result).iter().any(matches)
}

/// Nonzero formal charges as `(chain, atom name, charge)`.
fn charges(result: &StructureInterpretation) -> Vec<(String, String, i8)> {
    let topology = result.structure().topology();
    result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .filter_map(|atom| {
            let charge = topology
                .instance(atom.atom().instance())
                .unwrap()
                .molecule()
                .atom(atom.atom().atom())
                .unwrap()
                .formal_charge();
            (charge != 0).then(|| {
                (
                    atom.asym_id().to_owned(),
                    atom.atom_name().to_owned(),
                    charge,
                )
            })
        })
        .collect()
}

fn charged(rows: &[(&str, &str, i8)]) -> Vec<(String, String, i8)> {
    rows.iter()
        .map(|&(chain, name, charge)| (chain.to_owned(), name.to_owned(), charge))
        .collect()
}

/// The total hydrogens of each named atom of chain `asym` after perception.
fn hydrogens(result: &StructureInterpretation, asym: &str, name: &str) -> usize {
    let mut structure = result.structure().clone();
    structure.perceive().unwrap();
    let atom = result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .find(|atom| atom.asym_id() == asym && atom.atom_name() == name)
        .unwrap()
        .atom();
    structure
        .topology()
        .instance(atom.instance())
        .unwrap()
        .molecule()
        .total_hydrogens(atom.atom())
        .unwrap()
}

/// The order and evidence of the bond between two atoms of chain `asym`.
fn bond(
    result: &StructureInterpretation,
    asym: &str,
    a: &str,
    b: &str,
) -> (BondOrder, BondEvidence) {
    let report = result.report();
    let find = |name: &str| {
        report
            .instances()
            .iter()
            .flat_map(|instance| instance.atoms())
            .find(|atom| atom.asym_id() == asym && atom.atom_name() == name)
            .unwrap()
            .atom()
    };
    let (left, right) = (find(a), find(b));
    let evidence = report
        .instances()
        .iter()
        .flat_map(|instance| instance.bonds())
        .find(|bond| bond.atoms() == [left, right] || bond.atoms() == [right, left])
        .unwrap_or_else(|| panic!("{a} and {b} are bonded"))
        .evidence()
        .clone();
    let molecule = result
        .structure()
        .topology()
        .instance(left.instance())
        .unwrap()
        .molecule();
    let order = molecule
        .bond(molecule.bond_between(left.atom(), right.atom()).unwrap())
        .unwrap()
        .order();
    (order, evidence)
}

fn residue(asym: &str, sequence: i32) -> SourceResidueId {
    SourceResidueId::from_parts("1", asym, SourceResiduePosition::Label(sequence), None)
}

fn names(names: &[&str]) -> Vec<String> {
    names.iter().map(|&name| name.to_owned()).collect()
}

#[test]
fn a_five_prime_hydroxyl_hydrogen_fills_the_absent_phosphorus() {
    // A 5'-terminal nucleotide without its phosphate, as crystal structures
    // of DNA model it: HO5' takes the place of P on O5'.
    let mut rows = without(sites("DA"), &["OP3", "P", "OP1", "OP2", "HOP3", "HOP2"]);
    rows.push(hydrogen(&rows, "HO5'", "O5'", "C5'"));
    let result = read(&rows, "");
    assert_eq!(charges(&result), []);
    assert_eq!(
        bond(&result, "A", "HO5'", "O5'").1,
        BondEvidence::HydrogenName
    );
    assert_eq!(issues(&result), []);
}

#[test]
fn hydrogens_missing_from_nitrogen_and_alcohols_are_reported_not_charged() {
    // A residue whose amine hydrogens and hydroxyl hydrogen are not modelled.
    let rows = without(sites("THR"), &["H", "H2", "HG1"]);
    let result = read(&rows, "");
    assert_eq!(charges(&result), []);
    assert_eq!(hydrogens(&result, "A", "N"), 2);
    assert_eq!(hydrogens(&result, "A", "OG1"), 1);
    assert_eq!(
        issues(&result),
        [Issue::HydrogensMissing {
            residue: residue("A", 1),
            component_id: "THR".to_owned(),
            atoms: names(&["N", "OG1"]),
        }]
    );
    // An acid's missing hydrogen is its deprotonation.
    let acid = read(&without(sites("ASP"), &["HD2"]), "");
    assert_eq!(charges(&acid), charged(&[("A", "OD2", -1)]));
}

#[test]
fn terminal_hydrogens_bond_by_name() {
    for terminal in [["H1", "H2", "H3"], ["HT1", "HT2", "HT3"]] {
        let mut rows = without(sites("GLY"), &["H", "H2", "HXT"]);
        for name in terminal {
            rows.push(hydrogen(&rows, name, "N", "CA"));
        }
        let result = read(&rows, "");
        assert_eq!(
            charges(&result),
            charged(&[("A", "N", 1), ("A", "OXT", -1)]),
            "{terminal:?}"
        );
        assert!(
            !has_issue(&result, |issue| matches!(
                issue,
                Issue::HydrogenAttachedByDistance { .. }
            )),
            "{terminal:?}"
        );
    }
}

#[test]
fn unknown_template_charges_are_reported() {
    for (component, atoms) in [
        ("HEM", vec!["FE"]),
        ("SF4", vec!["FE1", "FE2", "FE3", "FE4"]),
    ] {
        let result = read(&heavy(sites(component)), "");
        assert!(
            has_issue(&result, |issue| matches!(
                issue,
                Issue::TemplateChargeUnknown { component_id, atoms: unknown, .. }
                    if component_id == component && *unknown == names(&atoms)
            )),
            "{component}: {:?}",
            issues(&result)
        );
    }
}

#[test]
fn block_component_charges_apply_without_a_template() {
    // Nitrate, which the bundle lacks, with the block's own component
    // tables: its charges come from `_chem_comp_atom`.
    let tables = "loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
                  _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n\
                  NO3 N N 1\nNO3 O1 O 0\nNO3 O2 O -1\nNO3 O3 O -1\n\
                  loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n\
                  _chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n\
                  NO3 N O1 doub\nNO3 N O2 sing\nNO3 N O3 sing\n";
    let nitrate = [
        ("N", [0.0, 0.0, 0.0]),
        ("O1", [1.25, 0.0, 0.0]),
        ("O2", [-0.625, 1.083, 0.0]),
        ("O3", [-0.625, -1.083, 0.0]),
    ]
    .map(|(name, point)| Site {
        element: name[..1].to_owned(),
        name: name.to_owned(),
        component: "NO3".to_owned(),
        asym: "A".to_owned(),
        sequence: 1,
        point,
        occupancy: None,
        charge: None,
    });
    let result = read(&nitrate, tables);
    assert_eq!(
        charges(&result),
        charged(&[("A", "N", 1), ("A", "O2", -1), ("A", "O3", -1)])
    );
    assert_eq!(hydrogens(&result, "A", "O2"), 0);
    // With bond orders alone, as wwPDB entry files state them, the
    // nitrogen's four bonds give it +1 and one single-bonded oxygen -1
    // beside it. Nothing states the other oxygen's charge, so it is the
    // hydroxyl of nitric acid, and the charges are reported.
    let bonds_only = &tables[tables.find("loop_\n_chem_comp_bond").unwrap()..];
    let result = read(&nitrate, bonds_only);
    assert_eq!(charges(&result), charged(&[("A", "N", 1), ("A", "O2", -1)]));
    assert_eq!(hydrogens(&result, "A", "O3"), 1);
    assert!(
        has_issue(&result, |issue| matches!(
            issue,
            Issue::ChargesFromBonds { .. }
        )),
        "{:?}",
        issues(&result)
    );

    // A single atom without template or stated charge is reported.
    let xenon = Site {
        element: "XE".to_owned(),
        name: "XE".to_owned(),
        component: "XE".to_owned(),
        asym: "B".to_owned(),
        sequence: 1,
        point: [0.0; 3],
        occupancy: None,
        charge: None,
    };
    let result = read(&[xenon], "");
    assert!(
        has_issue(&result, |issue| matches!(
            issue,
            Issue::TemplateMissing { component_id, .. } if component_id == "XE"
        )),
        "{:?}",
        issues(&result)
    );
}

fn metalc(rows: &[(&str, &str, i32, &str, &str)]) -> String {
    let mut text = String::from(
        "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
         _struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_atom_id\n\
         _struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr2_label_asym_id\n\
         _struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_label_seq_id\n",
    );
    for (index, (asym, atom, sequence, metal_asym, metal)) in rows.iter().enumerate() {
        text.push_str(&format!(
            "metalc{} metalc {asym} {atom} {sequence} {metal_asym} {metal} 1\n",
            index + 1
        ));
    }
    text
}

fn zinc(asym: &str, near: [f64; 3]) -> Site {
    Site {
        element: "ZN".to_owned(),
        name: "ZN".to_owned(),
        component: "ZN".to_owned(),
        asym: asym.to_owned(),
        sequence: 1,
        point: [near[0] + 2.3, near[1], near[2]],
        occupancy: None,
        charge: None,
    }
}

#[test]
fn metal_coordination_replaces_a_hydrogen() {
    // Zinc-bound cysteine is a thiolate, and histidine bound through NE2
    // carries its ring hydrogen on ND1.
    let mut rows = heavy(sites("CYS"));
    rows.push(zinc("Z", point(&rows, "SG")));
    let mut histidine = heavy(residue_sites("HIS", "B", 1, [20.0, 0.0, 0.0]));
    histidine.push(zinc("Y", point(&histidine, "NE2")));
    rows.extend(histidine);
    let result = read(
        &rows,
        &metalc(&[("A", "SG", 1, "Z", "ZN"), ("B", "NE2", 1, "Y", "ZN")]),
    );
    let found = charges(&result);
    assert!(
        found.contains(&("A".to_owned(), "SG".to_owned(), -1)),
        "{found:?}"
    );
    assert!(!found
        .iter()
        .any(|(asym, name, _)| asym == "B" && name.starts_with("ND1")));
    assert_eq!(hydrogens(&result, "B", "NE2"), 0);
    assert_eq!(hydrogens(&result, "B", "ND1"), 1);
    assert_eq!(bond(&result, "B", "ND1", "CE1").0, BondOrder::Single);
    assert_eq!(bond(&result, "B", "CE1", "NE2").0, BondOrder::Double);
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::HydrogenReplacedByCoordination { component_id, atoms, .. }
            if component_id == "CYS" && *atoms == names(&["SG"])
    )));

    // PDB states the coordination as LINK records, from which readers give
    // the histidine the same tautomer. Without them it would read as the
    // NE2-protonated tautomer, so PDB cannot carry its bonds.
    let (text, _) = kekule::pdb::write(&result, Default::default()).unwrap();
    let reread = kekule::pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(bond(&reread, "B", "ND1", "CE1").0, BondOrder::Single);
    assert_eq!(bond(&reread, "B", "CE1", "NE2").0, BondOrder::Double);
    assert_eq!(charges(&reread), charges(&result));
    let error = kekule::pdb::write(result.structure(), Default::default()).unwrap_err();
    assert!(matches!(
        error.detail(),
        Some(kekule::diagnostics::FormatErrorDetail::PdbWrite(
            kekule::pdb::PdbWriteErrorKind::DictionaryBondConflict { .. }
        ))
    ));
}

#[test]
fn zero_occupancy_hydrogens_are_not_protonation_evidence() {
    // Both ring hydrogens are placed but not observed: the histidine is
    // not read as the imidazolium cation.
    let rows = sites("HIS")
        .into_iter()
        .map(|site| {
            let ring = site.name == "HD1" || site.name == "HE2";
            site.occupancy(if ring { 0.0 } else { 1.0 })
        })
        .collect::<Vec<_>>();
    let result = read(&rows, "");
    assert!(
        !charges(&result)
            .iter()
            .any(|(_, name, _)| name == "ND1" || name == "NE2"),
        "{:?}",
        charges(&result)
    );
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::ZeroOccupancyHydrogens { atoms, .. } if *atoms == names(&["HD1", "HE2"])
    )));
}

#[test]
fn carboxylic_acids_protonated_on_either_oxygen_are_read() {
    for (component, removed, added, oxygen, other, carbon) in [
        ("ASP", "HD2", "HD1", "OD1", "OD2", "CG"),
        ("GLU", "HE2", "HE1", "OE1", "OE2", "CD"),
    ] {
        let mut rows = without(sites(component), &[removed]);
        rows.push(hydrogen(&rows, added, oxygen, carbon));
        let result = read(&rows, "");
        assert_eq!(charges(&result), [], "{component}");
        assert_eq!(bond(&result, "A", carbon, oxygen).0, BondOrder::Single);
        assert_eq!(bond(&result, "A", carbon, other).0, BondOrder::Double);
    }
}

#[test]
fn residues_without_hydrogens_take_their_ph_7_states() {
    let mut rows = heavy(sites("ASP"));
    rows.extend(heavy(residue_sites("HIS", "B", 1, [20.0, 0.0, 0.0])));
    rows.extend(heavy(residue_sites("LYS", "C", 1, [40.0, 0.0, 0.0])));
    rows.extend(heavy(residue_sites("CYS", "D", 1, [60.0, 0.0, 0.0])));
    let result = read(&rows, "");
    assert_eq!(
        charges(&result),
        charged(&[
            ("A", "N", 1),
            ("A", "OD2", -1),
            ("A", "OXT", -1),
            ("B", "N", 1),
            ("B", "OXT", -1),
            ("C", "N", 1),
            ("C", "NZ", 1),
            ("C", "OXT", -1),
            ("D", "N", 1),
            ("D", "OXT", -1),
        ])
    );
    assert_eq!(hydrogens(&result, "B", "NE2"), 1);
    assert_eq!(hydrogens(&result, "B", "ND1"), 0);
    assert_eq!(hydrogens(&result, "D", "SG"), 1);
    assert_eq!(
        bond(&result, "B", "CE1", "NE2").1,
        BondEvidence::Template {
            component_id: "HIS".to_owned(),
            variant: Some("NE2-protonated".to_owned()),
        }
    );
    let assumed = issues(&result)
        .into_iter()
        .filter(|issue| matches!(issue, Issue::ProtonationAssumed { .. }))
        .collect::<Vec<_>>();
    assert_eq!(
        assumed,
        [Issue::ProtonationAssumed {
            components: ["ASP", "CYS", "HIS", "LYS"]
                .map(|component| (component.to_owned(), 1))
                .to_vec(),
        }]
    );

    // A phosphodiester carries one negative charge.
    let first = without(heavy(sites("DA")), &["OP3", "P", "OP1", "OP2"]);
    let link = point(&first, "O3'");
    let second = heavy(residue_sites("DA", "A", 2, [0.0; 3]));
    let phosphorus = point(&second, "P");
    let shift = [0, 1, 2].map(|axis| link[axis] - phosphorus[axis] + [1.6, 0.0, 0.0][axis]);
    let second = without(residue_sites("DA", "A", 2, shift), &["OP3"]);
    let mut dna = first;
    dna.extend(heavy(second));
    let result = read(&dna, "");
    assert_eq!(charges(&result), charged(&[("A", "OP2", -1)]));
}

#[test]
fn heme_carboxylates_without_their_hydrogens_stay_single_anions() {
    // The dictionary charges heme's O2A and O2D -1 and still lists H2A and
    // H2D on them. Without those hydrogens they stay -1; with them they are
    // neutral acids, which the report names.
    let result = read(&residue_sites("HEM", "A", 0, [0.0; 3]), "");
    let oxygens = charges(&result)
        .into_iter()
        .filter(|(_, name, _)| name.starts_with('O'))
        .collect::<Vec<_>>();
    assert_eq!(oxygens, []);
    result.structure().clone().perceive().unwrap();
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::ProtonationFromHydrogens { component_id, .. } if component_id == "HEM"
    )));
    for rows in [
        without(residue_sites("HEM", "A", 0, [0.0; 3]), &["H2A", "H2D"]),
        heavy(residue_sites("HEM", "A", 0, [0.0; 3])),
    ] {
        let result = read(&rows, "");
        let oxygens = charges(&result)
            .into_iter()
            .filter(|(_, name, _)| name.starts_with('O'))
            .collect::<Vec<_>>();
        assert_eq!(oxygens, charged(&[("A", "O2A", -1), ("A", "O2D", -1)]));
        result.structure().clone().perceive().unwrap();
    }
}

#[test]
fn a_terminal_carboxylate_leaves_histidine_neutral() {
    // Histidine 2 of a sequence of two, without hydrogens: residue 1 is not
    // modelled, so N is an open amide, and OXT ends the chain.
    let rows = heavy(residue_sites("HIS", "A", 2, [0.0; 3]));
    let result = read_sequenced(&rows, &["ALA", "HIS"]);
    assert_eq!(charges(&result), charged(&[("A", "OXT", -1)]));
    assert_eq!(state_at(&result, "A", 2, "NE2"), (0, 1));
    assert_eq!(state_at(&result, "A", 2, "ND1"), (0, 0));
}

#[test]
fn a_stated_cation_after_a_gap_takes_three_hydrogens() {
    // Alanine 3 of three, without hydrogens, states N +1: the gap to
    // residue 2 is no bond, so N is an ammonium with three hydrogens.
    let mut rows = without(heavy(sites("ALA")), &["OXT"]);
    let mut third = heavy(residue_sites("ALA", "A", 3, [10.0, 0.0, 0.0]));
    state_charge(&mut third, "N", 1);
    rows.extend(third);
    let result = read_sequenced(&rows, &["ALA", "ALA", "ALA"]);
    assert_eq!(state_at(&result, "A", 3, "N"), (1, 3));
}

#[test]
fn a_stated_histidine_cation_without_hydrogens_is_the_imidazolium() {
    // Without hydrogens, NE2 +1 fits the doubly protonated histidine by
    // resonance: one hydrogen on each ring nitrogen.
    let mut rows = heavy(sites("HIS"));
    state_charge(&mut rows, "NE2", 1);
    let result = read(&rows, "");
    assert_eq!(state_at(&result, "A", 1, "NE2"), (1, 1));
    assert_eq!(state_at(&result, "A", 1, "ND1"), (0, 1));
}

#[test]
fn missing_template_atoms_are_reported_and_termini_stay_open() {
    // A carboxyl terminus without OXT is not an aldehyde.
    let result = read(&without(heavy(sites("GLY")), &["OXT"]), "");
    assert_eq!(hydrogens(&result, "A", "C"), 0);
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::TerminalAtomMissing { atom_name, anchor, .. } if atom_name == "OXT" && anchor == "C"
    )));
    // A truncated side chain is reported.
    let result = read(&without(heavy(sites("LYS")), &["CE", "NZ"]), "");
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::AtomsMissing { component_id, atoms, .. }
            if component_id == "LYS" && *atoms == names(&["CE", "NZ"])
    )));
    // A 5' phosphate without OP3 has no P-H bond.
    let result = read(&without(heavy(sites("DA")), &["OP3"]), "");
    assert_eq!(hydrogens(&result, "A", "P"), 0);
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::TerminalAtomMissing { atom_name, anchor, .. } if atom_name == "OP3" && anchor == "P"
    )));
}

/// The formal charge and total hydrogens, after perception, of atom `name`
/// of residue `sequence` of chain `asym`.
fn state_at(
    result: &StructureInterpretation,
    asym: &str,
    sequence: i32,
    name: &str,
) -> (i8, usize) {
    let mut structure = result.structure().clone();
    structure.perceive().unwrap();
    let atom = result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .find(|atom| {
            atom.asym_id() == asym
                && atom.label_sequence_id() == Some(sequence)
                && atom.atom_name() == name
        })
        .unwrap_or_else(|| panic!("{asym} {sequence} {name}"))
        .atom();
    let molecule = structure
        .topology()
        .instance(atom.instance())
        .unwrap()
        .molecule();
    (
        molecule.atom(atom.atom()).unwrap().formal_charge(),
        molecule.total_hydrogens(atom.atom()).unwrap(),
    )
}

/// The block with every row in entity `1`.
fn with_entity(text: &str) -> String {
    text.lines()
        .map(|line| {
            if line == "_atom_site.label_asym_id" {
                return format!("{line}\n_atom_site.label_entity_id");
            }
            if !line.starts_with("ATOM ") {
                return line.to_owned();
            }
            let mut fields = line.split(' ').collect::<Vec<_>>();
            fields.insert(6, "1");
            fields.join(" ")
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

fn link_missing(asym: &str, sequence: i32, component: &str, atom: &str) -> Issue {
    Issue::PolymerLinkMissing {
        residue: residue(asym, sequence),
        component_id: component.to_owned(),
        atom_name: atom.to_owned(),
    }
}

/// An `_entity_poly_seq` loop for entity `1`.
fn entity_sequence(monomers: &[&str]) -> String {
    let mut text =
        "loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n"
            .to_owned();
    for (number, monomer) in (1..).zip(monomers) {
        text.push_str(&format!("1 {number} {monomer}\n"));
    }
    text
}

/// The rows as entity `1`, whose sequence is `monomers`.
fn read_sequenced(rows: &[Site], monomers: &[&str]) -> StructureInterpretation {
    let text = with_entity(&block(rows, &entity_sequence(monomers)));
    mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

#[test]
fn a_residue_after_a_gap_is_an_open_amide() {
    // Alanines 1 and 3 of a sequence of three, without hydrogens: only the
    // first is an N terminus.
    let mut rows = without(heavy(sites("ALA")), &["OXT"]);
    rows.extend(heavy(residue_sites("ALA", "A", 3, [10.0, 0.0, 0.0])));
    let unsequenced = read(&rows, "");
    let result = read_sequenced(&rows, &["ALA", "ALA", "ALA"]);
    assert_eq!(state_at(&result, "A", 1, "N"), (1, 3));
    assert_eq!(state_at(&result, "A", 3, "N"), (0, 1));
    assert_eq!(state_at(&result, "A", 1, "C"), (0, 0));
    assert_eq!(state_at(&result, "A", 3, "OXT"), (-1, 0));
    let found = issues(&result);
    assert!(
        found.contains(&link_missing("A", 1, "ALA", "C")),
        "{found:?}"
    );
    assert!(
        found.contains(&link_missing("A", 3, "ALA", "N")),
        "{found:?}"
    );
    assert!(
        !found
            .iter()
            .any(|issue| matches!(issue, Issue::TerminalAtomMissing { .. })),
        "{found:?}"
    );
    // Only a sequence says a polymer lacks residues: without one, residue
    // 3 starts a chain of its own.
    assert_eq!(state_at(&unsequenced, "A", 3, "N"), (1, 3));
    assert!(!has_issue(&unsequenced, |issue| matches!(
        issue,
        Issue::PolymerLinkMissing { .. }
    )));
}

/// Alanines 2 and 3 of an entity of four residues, without hydrogens or
/// `OXT`.
fn inner_alanines() -> StructureInterpretation {
    let mut rows = without(heavy(residue_sites("ALA", "A", 2, [0.0; 3])), &["OXT"]);
    let carbon = point(&rows, "C");
    let next = heavy(residue_sites("ALA", "A", 3, [0.0; 3]));
    let nitrogen = point(&next, "N");
    let shift = [0, 1, 2].map(|axis| carbon[axis] - nitrogen[axis] + [1.33, 0.0, 0.0][axis]);
    rows.extend(without(
        heavy(residue_sites("ALA", "A", 3, shift)),
        &["OXT"],
    ));
    let sequence = "loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n\
                    _entity_poly_seq.mon_id\n1 1 MET\n1 2 ALA\n1 3 ALA\n1 4 GLY\n";
    mmcif::parse_str(&with_entity(&block(&rows, sequence)), Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

#[test]
fn a_sequence_that_continues_past_the_model_opens_its_ends() {
    let result = inner_alanines();
    assert!(!has_issue(&result, |issue| matches!(
        issue,
        Issue::ChainBreak { .. }
    )));
    assert_eq!(state_at(&result, "A", 2, "N"), (0, 1));
    assert_eq!(state_at(&result, "A", 3, "C"), (0, 0));
    let found = issues(&result);
    assert!(
        found.contains(&link_missing("A", 2, "ALA", "N")),
        "{found:?}"
    );
    assert!(
        found.contains(&link_missing("A", 3, "ALA", "C")),
        "{found:?}"
    );
    assert!(
        !found
            .iter()
            .any(|issue| matches!(issue, Issue::TerminalAtomMissing { .. })),
        "{found:?}"
    );

    // mmCIF states the sequence in `_entity_poly_seq`, PDB in `SEQRES`.
    let (text, _) = mmcif::write([&result], Default::default()).unwrap();
    assert!(
        text.contains("1 1 MET\n1 2 ALA\n1 3 ALA\n1 4 GLY\n"),
        "{text}"
    );
    let reread = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(state_at(&reread, "A", 2, "N"), (0, 1));
    assert_eq!(state_at(&reread, "A", 3, "C"), (0, 0));
    let links = reread
        .report()
        .issues()
        .iter()
        .filter(|issue| matches!(issue, Issue::PolymerLinkMissing { .. }))
        .count();
    assert_eq!(links, 2);
    let (text, _) = kekule::pdb::write(&result, Default::default()).unwrap();
    assert!(
        text.starts_with("SEQRES   1 A    4  MET ALA ALA GLY\n"),
        "{text}"
    );
    let reread = kekule::pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let links = reread
        .report()
        .issues()
        .iter()
        .filter(|issue| matches!(issue, Issue::PolymerLinkMissing { .. }))
        .count();
    assert_eq!(links, 2);
    assert_eq!(state_at(&reread, "A", 3, "C"), (0, 0));
    assert_eq!(
        kekule::pdb::write(&reread, Default::default()).unwrap().0,
        text
    );
    // mmCIF written from PDB numbers the chain from 2, after the gap.
    let (text, _) = mmcif::write([&reread], Default::default()).unwrap();
    let reread = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(state_at(&reread, "A", 2, "N"), (0, 1));
}

#[test]
fn ligands_take_the_ph_7_state_of_their_acids_only() {
    // Heavy-atom GTP and a free lysine as ligands: the phosphates lose their
    // acidic hydrogens, and the lysine's amines keep the template's state.
    let mut rows = heavy(residue_sites("GTP", "L", 0, [0.0; 3]));
    rows.extend(heavy(residue_sites("LYS", "M", 0, [30.0, 0.0, 0.0])));
    let result = read(&rows, "");
    let found = charges(&result);
    let gtp = found
        .iter()
        .filter(|(asym, _, _)| asym == "L")
        .map(|(_, _, charge)| i32::from(*charge))
        .sum::<i32>();
    assert_eq!(gtp, -4, "{found:?}");
    assert!(found
        .iter()
        .filter(|(asym, _, _)| asym == "L")
        .all(|(_, name, charge)| *charge == -1 && name.starts_with('O')));
    assert_eq!(
        found
            .iter()
            .filter(|(asym, _, _)| asym == "M")
            .cloned()
            .collect::<Vec<_>>(),
        charged(&[("M", "NZ", 1), ("M", "OXT", -1)])
    );
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::ProtonationAssumed { components }
            if *components == [("GTP".to_owned(), 1), ("LYS".to_owned(), 1)]
    )));
    assert!(issues(&result)
        .iter()
        .find(|issue| matches!(issue, Issue::ProtonationAssumed { .. }))
        .unwrap()
        .to_string()
        .contains("Ligands take nothing else."));
}

#[test]
fn a_cationic_nitrogen_with_its_neutral_hydrogens_is_the_neutral_base() {
    // Lysine with two hydrogens on NZ is neutral, not short of one.
    let lysine = read(&without(sites("LYS"), &["HZ3"]), "");
    assert_eq!(charges(&lysine), []);
    assert_eq!(hydrogens(&lysine, "A", "NZ"), 2);
    assert_eq!(
        issues(&lysine),
        [Issue::ProtonationFromHydrogens {
            residue: residue("A", 1),
            component_id: "LYS".to_owned(),
            variant: None,
            charged_atoms: Vec::new(),
        }]
    );
    // Arginine without HH12 is the neutral base with NH1=CZ.
    let arginine = read(&without(sites("ARG"), &["HH12"]), "");
    assert_eq!(charges(&arginine), []);
    assert_eq!(hydrogens(&arginine, "A", "NH1"), 1);
    assert_eq!(bond(&arginine, "A", "CZ", "NH1").0, BondOrder::Double);
    assert_eq!(bond(&arginine, "A", "CZ", "NH2").0, BondOrder::Single);
    assert!(!has_issue(&arginine, |issue| matches!(
        issue,
        Issue::HydrogensMissing { .. }
    )));
    // Without any of its hydrogens, NZ keeps the template's charge.
    let bare = read(&without(sites("LYS"), &["HZ1", "HZ2", "HZ3"]), "");
    assert_eq!(charges(&bare), charged(&[("A", "NZ", 1)]));
    assert!(has_issue(&bare, |issue| matches!(
        issue,
        Issue::HydrogensMissing { atoms, .. } if *atoms == names(&["NZ"])
    )));
}

#[test]
fn stated_charges_set_the_hydrogens_they_need() {
    // A stated charge of 0 in a residue without hydrogens is no evidence:
    // files without charges state 0 everywhere.
    let mut aspartate = heavy(sites("ASP"));
    state_charge(&mut aspartate, "OD2", 0);
    let result = read(&aspartate, "");
    assert!(charges(&result).contains(&("A".to_owned(), "OD2".to_owned(), -1)));
    assert_eq!(hydrogens(&result, "A", "OD2"), 0);
    assert!(!has_issue(&result, |issue| matches!(
        issue,
        Issue::StatedChargeKept { .. } | Issue::StatedChargeRejected { .. }
    )));
    // A stated phenolate explains its absent hydrogen.
    let mut tyrosine = without(sites("TYR"), &["HH"]);
    state_charge(&mut tyrosine, "OH", -1);
    let result = read(&tyrosine, "");
    assert_eq!(charges(&result), charged(&[("A", "OH", -1)]));
    assert_eq!(
        issues(&result),
        [Issue::StatedChargeKept {
            residue: residue("A", 1),
            component_id: "TYR".to_owned(),
            atom_name: "OH".to_owned(),
            stated: -1,
            derived: 0,
        }]
    );
    // A stated 0 on lysine NZ without hydrogens keeps the ammonium, while
    // a stated nonzero charge stays evidence.
    let mut lysine = heavy(sites("LYS"));
    state_charge(&mut lysine, "NZ", 0);
    let ammonium = read(&lysine, "");
    assert!(charges(&ammonium).contains(&("A".to_owned(), "NZ".to_owned(), 1)));
    assert_eq!(hydrogens(&ammonium, "A", "NZ"), 3);
    let mut glutamate = heavy(sites("GLU"));
    state_charge(&mut glutamate, "OE1", -1);
    let anion = read(&glutamate, "");
    assert!(charges(&anion).contains(&("A".to_owned(), "OE1".to_owned(), -1)));
    assert_eq!(hydrogens(&anion, "A", "OE1"), 0);
    // A stated +1 on a carboxylate oxygen is held as -OH2+ and reported.
    let mut cation = heavy(sites("ASP"));
    state_charge(&mut cation, "OD2", 1);
    let kept = read(&cation, "");
    assert!(charges(&kept).contains(&("A".to_owned(), "OD2".to_owned(), 1)));
    assert_eq!(hydrogens(&kept, "A", "OD2"), 2);
    assert!(has_issue(&kept, |issue| matches!(
        issue,
        Issue::StatedChargeKept { atom_name, stated: 1, derived: -1, .. } if atom_name == "OD2"
    )));
    // No hydrogen count holds -2 on a hydroxyl oxygen.
    let mut serine = heavy(sites("SER"));
    state_charge(&mut serine, "OG", -2);
    let rejected = read(&serine, "");
    assert!(has_issue(&rejected, |issue| matches!(
        issue,
        Issue::StatedChargeRejected { atom_name, stated: -2, kept: 0, .. } if atom_name == "OG"
    )));
    // A stated charge equal to the derived one is no issue of its own.
    let mut aspartate = without(sites("ASP"), &["HD2"]);
    state_charge(&mut aspartate, "OD2", -1);
    let stated = read(&aspartate, "");
    assert_eq!(charges(&stated), charged(&[("A", "OD2", -1)]));
    assert_eq!(
        issues(&stated),
        [Issue::ProtonationFromHydrogens {
            residue: residue("A", 1),
            component_id: "ASP".to_owned(),
            variant: None,
            charged_atoms: vec![("OD2".to_owned(), -1)],
        }]
    );
}

#[test]
fn stated_charges_of_another_resonance_form_move_its_double_bond() {
    let mut rows = heavy(sites("ASP"));
    state_charge(&mut rows, "OD1", -1);
    let mut arginine = heavy(residue_sites("ARG", "B", 1, [20.0, 0.0, 0.0]));
    state_charge(&mut arginine, "NH1", 1);
    rows.extend(arginine);
    let result = read(&rows, "");
    let found = charges(&result);
    for expected in [("A", "OD1", -1), ("B", "NH1", 1)] {
        assert!(found.contains(&charged(&[expected])[0]), "{found:?}");
    }
    for absent in ["OD2", "NH2"] {
        assert!(
            !found.iter().any(|(_, name, _)| name == absent),
            "{found:?}"
        );
    }
    assert_eq!(bond(&result, "A", "CG", "OD1").0, BondOrder::Single);
    assert_eq!(bond(&result, "A", "CG", "OD2").0, BondOrder::Double);
    assert_eq!(bond(&result, "B", "CZ", "NH1").0, BondOrder::Double);
    assert_eq!(bond(&result, "B", "CZ", "NH2").0, BondOrder::Single);
    assert_eq!(hydrogens(&result, "B", "NH1"), 2);
    let resonance = issues(&result)
        .into_iter()
        .filter(|issue| matches!(issue, Issue::StatedChargeResonance { .. }))
        .count();
    assert_eq!(resonance, 2);
    assert!(!has_issue(&result, |issue| matches!(
        issue,
        Issue::StatedChargeRejected { .. }
    )));
    // Writers state the moved bonds' charges, from which readers move them
    // again.
    let (text, _) = mmcif::write([&result], Default::default()).unwrap();
    assert!(!text.contains("CG OD2"), "{text}");
    let reread = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(bond(&reread, "A", "CG", "OD2").0, BondOrder::Double);
}

#[test]
fn a_carboxylate_bound_through_either_oxygen_keeps_its_template() {
    let mut rows = heavy(sites("ASP"));
    rows.push(zinc("Z", point(&rows, "OD2")));
    let result = read(&rows, &metalc(&[("A", "OD2", 1, "Z", "ZN")]));
    assert!(charges(&result).contains(&("A".to_owned(), "OD2".to_owned(), -1)));
    assert_eq!(bond(&result, "A", "CG", "OD1").0, BondOrder::Double);
    assert_eq!(
        bond(&result, "A", "CG", "OD1").1,
        BondEvidence::Template {
            component_id: "ASP".to_owned(),
            variant: None,
        }
    );
    // Without the coordination, the template's bonds are still the ones
    // readers give it.
    kekule::pdb::write(result.structure(), Default::default()).unwrap();
}

#[test]
fn phenols_and_thiols_without_their_hydrogen_stay_neutral() {
    let mut rows = without(sites("TYR"), &["HH"]);
    rows.extend(without(
        residue_sites("CYS", "B", 1, [20.0, 0.0, 0.0]),
        &["HG"],
    ));
    let result = read(&rows, "");
    assert_eq!(charges(&result), []);
    assert_eq!(hydrogens(&result, "A", "OH"), 1);
    assert_eq!(hydrogens(&result, "B", "SG"), 1);
    let missing = issues(&result)
        .into_iter()
        .filter_map(|issue| match issue {
            Issue::HydrogensMissing { atoms, .. } => Some(atoms),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(missing, [names(&["OH"]), names(&["SG"])]);
}

#[test]
fn a_hydrogen_and_a_deuterium_sharing_a_site_are_alternates() {
    // A neutron structure's serine: HG and DG share OG, half each, 0.5 Å
    // apart. The deuterium stays on the tie, and the hydrogen is omitted.
    let mut rows = sites("SER");
    let hg = rows.iter().position(|site| site.name == "HG").unwrap();
    rows[hg].occupancy = Some(0.5);
    let mut deuterium = Site {
        element: "D".to_owned(),
        name: "DG".to_owned(),
        ..rows[hg].clone()
    };
    deuterium.point[0] += 0.5;
    rows.push(deuterium);
    let result = read(&rows, "");
    assert_eq!(charges(&result), []);
    assert_eq!(state_at(&result, "A", 1, "OG"), (0, 1));
    let omitted = issues(&result)
        .into_iter()
        .filter_map(|issue| match issue {
            Issue::AlternateLocationOmitted {
                atom_name, alt_id, ..
            } => Some((atom_name, alt_id)),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(omitted, [("HG".to_owned(), None)]);
    assert!(!has_issue(&result, |issue| matches!(
        issue,
        Issue::ProtonationFromHydrogens { .. }
    )));

    // A larger occupancy stays.
    rows[hg].occupancy = Some(0.6);
    rows.last_mut().unwrap().occupancy = Some(0.4);
    let result = read(&rows, "");
    assert_eq!(state_at(&result, "A", 1, "OG"), (0, 1));
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::AlternateLocationOmitted { atom_name, .. } if atom_name == "DG"
    )));
}

#[test]
fn writers_take_the_component_library_readers_use() {
    // A component that only the caller's library has, without OXT: its C
    // keeps no hydrogen, which only that library's readers fix again.
    let path = format!("{}/data/components/ALA.cif", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path)
        .unwrap()
        .replace("data_ALA", "data_XAL")
        .replace("ALA ", "XAL ");
    let library = std::sync::Arc::new(
        kekule::components::ComponentLibrary::from_mmcif(
            &mmcif::parse_str(&text, Default::default()).unwrap(),
        )
        .unwrap(),
    );
    let mut rows = without(heavy(sites("ALA")), &["OXT"]);
    for row in &mut rows {
        row.component = "XAL".to_owned();
    }
    let options = kekule::macromolecular::StructureInterpretOptions::default().with_connectivity(
        kekule::macromolecular::ConnectivityOptions::default()
            .with_components(Some(library.clone())),
    );
    let result = mmcif::parse_str(&block(&rows, ""), Default::default())
        .unwrap()
        .interpret(options)
        .unwrap();
    assert_eq!(hydrogens(&result, "A", "C"), 0);
    assert!(mmcif::write([&result], Default::default()).is_err());
    assert!(kekule::pdb::write(&result, Default::default()).is_err());
    mmcif::write(
        [&result],
        mmcif::MmcifWriteOptions::default().with_components(Some(library.clone())),
    )
    .unwrap();
    kekule::pdb::write(
        &result,
        kekule::pdb::PdbWriteOptions::default().with_components(Some(library)),
    )
    .unwrap();
}

/// Glycines 2 and 3 of chain A, peptide-bonded, at the template's ideal
/// coordinates and with all their hydrogens, without the leaving atoms of
/// the bond between them.
fn glycine_pair() -> Vec<Site> {
    let mut rows = without(residue_sites("GLY", "A", 2, [0.0; 3]), &["OXT", "HXT"]);
    let carbon = point(&rows, "C");
    let next = residue_sites("GLY", "A", 3, [0.0; 3]);
    let nitrogen = point(&next, "N");
    let shift = [0, 1, 2].map(|axis| carbon[axis] - nitrogen[axis] + [1.33, 0.0, 0.0][axis]);
    rows.extend(without(residue_sites("GLY", "A", 3, shift), &["H2"]));
    rows
}

#[test]
fn explicit_terminal_atoms_end_a_chain_whatever_the_sequence_says() {
    // Glycines 2 and 3 of a sequence of four: residue 2 has the three
    // hydrogens of an N terminus and residue 3 an OXT, so the chain ends
    // at both, though the sequence goes on.
    let mut rows = glycine_pair();
    rows.retain(|site| {
        !(site.sequence == 2 && ["H", "H2"].contains(&site.name.as_str()))
            && !(site.sequence == 3 && site.name == "HXT")
    });
    for name in ["H1", "H2", "H3"] {
        // The first N and CA are residue 2's.
        rows.push(hydrogen(&rows, name, "N", "CA"));
    }
    let result = read_sequenced(&rows, &["MET", "GLY", "GLY", "ALA"]);
    assert_eq!(state_at(&result, "A", 2, "N"), (1, 3));
    assert_eq!(state_at(&result, "A", 3, "OXT"), (-1, 0));
    assert_eq!(state_at(&result, "A", 3, "C"), (0, 0));
    assert!(
        !has_issue(&result, |issue| matches!(
            issue,
            Issue::PolymerLinkMissing { .. }
        )),
        "{:?}",
        issues(&result)
    );
    // Without the terminal hydrogens and OXT, only the sequence speaks.
    let rows = without(heavy(glycine_pair()), &["OXT"]);
    let open = read_sequenced(&rows, &["MET", "GLY", "GLY", "ALA"]);
    assert_eq!(state_at(&open, "A", 2, "N"), (0, 1));
    assert!(has_issue(&open, |issue| *issue == link_missing("A", 2, "GLY", "N")));
}

/// PDB `ATOM` records of the sites in chain `chain`, numbered by their
/// sequence, with a `TER` after each listed sequence number.
fn pdb_text(sites: &[Site], chain: char, ter_after: &[i32]) -> String {
    let mut text = String::new();
    for (index, site) in sites.iter().enumerate() {
        let name = if site.name.len() < 4 {
            format!(" {:<3}", site.name)
        } else {
            site.name.clone()
        };
        text.push_str(&format!(
            "ATOM  {:>5} {name:<4} {:>3} {chain}{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00 20.00          {:>2}\n",
            index + 1,
            site.component,
            site.sequence,
            site.point[0],
            site.point[1],
            site.point[2],
            site.element,
        ));
        let last = sites
            .get(index + 1)
            .is_none_or(|next| next.sequence != site.sequence);
        if last && ter_after.contains(&site.sequence) {
            text.push_str("TER\n");
        }
    }
    text + "END\n"
}

#[test]
fn modelled_residues_elsewhere_make_no_gap() {
    // Two glycines of a chain without `SEQRES`, numbered 1 and 3: each is
    // a whole molecule, with and without hydrogens, in one TER segment or
    // two, in a chain with a letter or without.
    for (chain, with_hydrogens) in [(' ', true), (' ', false), ('A', true), ('A', false)] {
        for ter_after in [&[1, 3][..], &[3]] {
            let mut rows = residue_sites("GLY", "A", 1, [0.0; 3]);
            rows.extend(residue_sites("GLY", "A", 3, [12.0, 0.0, 0.0]));
            rows.retain(|site| with_hydrogens || site.element != "H");
            let text = pdb_text(&rows, chain, ter_after);
            let result = kekule::pdb::parse_str(&text, Default::default())
                .unwrap()
                .interpret(Default::default())
                .unwrap();
            let found = result.report().issues().to_vec();
            assert!(
                !found
                    .iter()
                    .any(|issue| matches!(issue, Issue::PolymerLinkMissing { .. })),
                "{with_hydrogens} {ter_after:?}: {found:?}"
            );
            // Each N is a terminal amine or ammonium, not an open amide.
            let mut structure = result.structure().clone();
            structure.perceive().unwrap();
            let amines = result
                .report()
                .instances()
                .iter()
                .flat_map(|instance| instance.atoms())
                .filter(|atom| atom.atom_name() == "N")
                .map(|atom| {
                    let atom = atom.atom();
                    structure
                        .topology()
                        .instance(atom.instance())
                        .unwrap()
                        .molecule()
                        .total_hydrogens(atom.atom())
                        .unwrap()
                })
                .collect::<Vec<_>>();
            assert_eq!(amines.len(), 2);
            assert!(
                amines.iter().all(|&count| count >= 2),
                "{with_hydrogens} {ter_after:?}: {amines:?}"
            );
            // mmCIF gives each molecule a sequence of its own, and a blank
            // chain a label.
            let (written, _) = mmcif::write([&result], Default::default()).unwrap();
            let reread = mmcif::parse_str(&written, Default::default())
                .unwrap()
                .interpret(Default::default())
                .unwrap();
            assert!(
                !reread
                    .report()
                    .issues()
                    .iter()
                    .any(|issue| matches!(issue, Issue::PolymerLinkMissing { .. })),
                "{written}"
            );
        }
    }
}

/// The bond orders between `centre` and each of `partners` in chain `asym`.
fn orders(
    result: &StructureInterpretation,
    asym: &str,
    centre: &str,
    partners: &[&str],
) -> Vec<BondOrder> {
    partners
        .iter()
        .map(|partner| bond(result, asym, centre, partner).0)
        .collect()
}

#[test]
fn stated_charges_hold_where_valence_allows_them() {
    // Each input reads to a structure that perception accepts.
    let perceives = |result: &StructureInterpretation| {
        let mut structure = result.structure().clone();
        structure
            .perceive()
            .map(|_| ())
            .map_err(|error| error.to_string())
    };
    // Histidine with HD1 and HE2 and +1 stated on NE2: the ring's double
    // bond moves to CE1=NE2.
    let mut histidine = sites("HIS");
    state_charge(&mut histidine, "NE2", 1);
    let result = read(&histidine, "");
    assert_eq!(perceives(&result), Ok(()));
    let found = charges(&result);
    assert!(
        found.contains(&("A".to_owned(), "NE2".to_owned(), 1)),
        "{found:?}"
    );
    assert!(!found.iter().any(|(_, name, _)| name == "ND1"), "{found:?}");
    assert_eq!(bond(&result, "A", "CE1", "NE2").0, BondOrder::Double);
    assert_eq!(hydrogens(&result, "A", "NE2"), 1);
    // Arginine with all its hydrogens and +1 stated on NE: CZ=NE.
    let mut arginine = sites("ARG");
    state_charge(&mut arginine, "NE", 1);
    let result = read(&arginine, "");
    assert_eq!(perceives(&result), Ok(()));
    let found = charges(&result);
    assert!(
        found.contains(&("A".to_owned(), "NE".to_owned(), 1)),
        "{found:?}"
    );
    assert!(
        !found.iter().any(|(_, name, _)| name.starts_with("NH")),
        "{found:?}"
    );
    assert_eq!(
        orders(&result, "A", "CZ", &["NE", "NH1", "NH2"]),
        [BondOrder::Double, BondOrder::Single, BondOrder::Single]
    );
    assert!(!has_issue(&result, |issue| matches!(
        issue,
        Issue::StatedChargeRejected { .. }
    )));
    // Charge-separated DMSO, phosphate, and sulfate: the S=O and P=O
    // bonds become single between the stated charges.
    let mut dmso = sites("DMS");
    state_charge(&mut dmso, "S", 1);
    state_charge(&mut dmso, "O", -1);
    let result = read(&dmso, "");
    assert_eq!(perceives(&result), Ok(()));
    assert_eq!(charges(&result), charged(&[("A", "S", 1), ("A", "O", -1)]));
    assert_eq!(bond(&result, "A", "S", "O").0, BondOrder::Single);
    for (component, centre, charge) in [("PO4", "P", 1), ("SO4", "S", 2)] {
        let mut rows = heavy(sites(component));
        state_charge(&mut rows, centre, charge);
        for oxygen in ["O1", "O2", "O3", "O4"] {
            state_charge(&mut rows, oxygen, -1);
        }
        let result = read(&rows, "");
        assert_eq!(perceives(&result), Ok(()), "{component}");
        assert_eq!(
            orders(&result, "A", centre, &["O1", "O2", "O3", "O4"]),
            [BondOrder::Single; 4],
            "{component}"
        );
        let mut expected = charged(&[
            ("A", "O1", -1),
            ("A", "O2", -1),
            ("A", "O3", -1),
            ("A", "O4", -1),
            ("A", centre, charge),
        ]);
        expected.sort();
        let mut found = charges(&result);
        found.sort();
        assert_eq!(found, expected, "{component}");
        for oxygen in ["O1", "O2", "O3", "O4"] {
            assert_eq!(hydrogens(&result, "A", oxygen), 0, "{component} {oxygen}");
        }
    }
    // With explicit hydrogens, no implicit one is added for a stated
    // charge: +1 on a serine's OG with its one hydrogen is rejected.
    let mut serine = sites("SER");
    state_charge(&mut serine, "OG", 1);
    let result = read(&serine, "");
    assert_eq!(perceives(&result), Ok(()));
    assert_eq!(hydrogens(&result, "A", "OG"), 1);
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::StatedChargeRejected { atom_name, stated: 1, kept: 0, .. } if atom_name == "OG"
    )));
}

#[test]
fn a_hydrogen_and_a_deuterium_at_separate_sites_both_count() {
    // Lysine with DZ3 in place of HZ3: three sites on NZ, so it is the
    // ammonium, not a neutral amine with a shared site.
    let mut rows = sites("LYS");
    let hz3 = rows.iter().position(|site| site.name == "HZ3").unwrap();
    rows[hz3].element = "D".to_owned();
    rows[hz3].name = "DZ3".to_owned();
    let result = read(&rows, "");
    assert_eq!(charges(&result), charged(&[("A", "NZ", 1)]));
    assert_eq!(hydrogens(&result, "A", "NZ"), 3);
}

#[test]
fn a_base_whose_neutral_form_has_no_hydrogen_there_is_neutral_without_it() {
    // Histidine as a component of the caller's library, which has none of
    // the bundle's tautomers: the imidazolium with +1 on ND1. Without HD1,
    // ND1 is the neutral base, which has no hydrogen there.
    let path = format!("{}/data/components/HIS.cif", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path)
        .unwrap()
        .replace("data_HIS", "data_XHS")
        .replace("HIS ", "XHS ");
    let library = std::sync::Arc::new(
        kekule::components::ComponentLibrary::from_mmcif(
            &mmcif::parse_str(&text, Default::default()).unwrap(),
        )
        .unwrap(),
    );
    let mut rows = without(sites("HIS"), &["HD1"]);
    for row in &mut rows {
        row.component = "XHS".to_owned();
    }
    let options = kekule::macromolecular::StructureInterpretOptions::default().with_connectivity(
        kekule::macromolecular::ConnectivityOptions::default().with_components(Some(library)),
    );
    let result = mmcif::parse_str(&block(&rows, ""), Default::default())
        .unwrap()
        .interpret(options)
        .unwrap();
    assert_eq!(charges(&result), []);
    assert_eq!(hydrogens(&result, "A", "ND1"), 0);
    assert_eq!(hydrogens(&result, "A", "NE2"), 1);
    assert!(!has_issue(&result, |issue| matches!(
        issue,
        Issue::HydrogensMissing { .. }
    )));
}

#[test]
fn a_lone_first_residue_takes_its_last_place_before_the_next_segment() {
    // A His-tagged sequence whose modelled residues are a lone histidine
    // numbered 1 and alanines numbered 20 to 21: no one numbering fits
    // both, so the histidine takes the last histidine before the
    // alanines, not the tag's first.
    let sequence = [
        "MET", "HIS", "HIS", "HIS", "HIS", "HIS", "HIS", "GLY", "HIS", "SER", "ALA", "ALA",
    ];
    let mut rows = heavy(residue_sites("HIS", "A", 1, [0.0; 3]));
    rows.extend(heavy(residue_sites("ALA", "A", 20, [15.0, 0.0, 0.0])));
    rows.extend(heavy(residue_sites("ALA", "A", 21, [30.0, 0.0, 0.0])));
    let mut text = String::new();
    for (index, chunk) in sequence.chunks(13).enumerate() {
        text.push_str(&format!(
            "SEQRES {:>3} A {:>4}  {}\n",
            index + 1,
            sequence.len(),
            chunk.join(" ")
        ));
    }
    text.push_str(&pdb_text(&rows, 'A', &[21]));
    let result = kekule::pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let places = result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .filter(|atom| atom.atom_name() == "CA")
        .map(|atom| (atom.component_id().to_owned(), atom.label_sequence_id()))
        .collect::<Vec<_>>();
    assert_eq!(
        places,
        [
            ("HIS".to_owned(), Some(9)),
            ("ALA".to_owned(), Some(11)),
            ("ALA".to_owned(), Some(12)),
        ]
    );
}

/// `SEQRES` records of chain `chain` listing `sequence`.
fn seqres(chain: char, sequence: &[&str]) -> String {
    sequence
        .chunks(13)
        .enumerate()
        .map(|(serial, names)| {
            let names = names
                .iter()
                .map(|name| format!("{name:>3}"))
                .collect::<Vec<_>>()
                .join(" ");
            format!(
                "SEQRES {:>3} {chain} {:>4}  {names}\n",
                serial + 1,
                sequence.len()
            )
        })
        .collect()
}

/// Residues of chain `A` without hydrogens or `OXT`, numbered from
/// `first`, apart from each other.
fn separate_residues(components: &[&str], first: i32) -> Vec<Site> {
    let mut rows = Vec::new();
    for (step, &component) in (0..).zip(components) {
        let offset = [12.0 * f64::from(step), 0.0, 0.0];
        rows.extend(without(
            heavy(residue_sites(component, "A", first + step, offset)),
            &["OXT"],
        ));
    }
    rows
}

/// The formal charge and total hydrogens, after perception, of atom `name`
/// of the PDB residue numbered `number`.
fn pdb_state_at(result: &StructureInterpretation, number: i32, name: &str) -> (i8, usize) {
    let mut structure = result.structure().clone();
    structure.perceive().unwrap();
    let number = number.to_string();
    let atom = result
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .find(|atom| atom.author_sequence_id() == Some(number.as_str()) && atom.atom_name() == name)
        .unwrap_or_else(|| panic!("{number} {name}"))
        .atom();
    let molecule = structure
        .topology()
        .instance(atom.instance())
        .unwrap()
        .molecule();
    (
        molecule.atom(atom.atom()).unwrap().formal_charge(),
        molecule.total_hydrogens(atom.atom()).unwrap(),
    )
}

fn sequence_not_placed(result: &StructureInterpretation) -> bool {
    has_issue(
        result,
        |issue| matches!(issue, Issue::SequenceNotPlaced { chain_id } if chain_id == "A"),
    )
}

#[test]
fn a_sequence_with_a_substituted_residue_places_its_chain() {
    // Residues 2 to 4 of a mutant whose SEQRES lists TRP where the model
    // has PHE: the sequence still places them, so the N of residue 2 is an
    // amide after a gap.
    let rows = separate_residues(&["GLY", "PHE", "LYS"], 2);
    let text = seqres('A', &["SER", "GLY", "TRP", "LYS", "ALA"]) + &pdb_text(&rows, 'A', &[4]);
    let result = read_pdb(&text);
    assert!(!sequence_not_placed(&result));
    assert_eq!(pdb_state_at(&result, 2, "N"), (0, 1));
    assert_eq!(pdb_state_at(&result, 4, "NZ"), (1, 3));
    assert!(has_issue(&result, |issue| matches!(
        issue,
        Issue::PolymerLinkMissing { component_id, atom_name, .. }
            if component_id == "LYS" && atom_name == "C"
    )));

    // A sequence that fits none of them is reported, and the residues keep
    // their record order: the first is the chain's start.
    let text = seqres('A', &["ALA"; 5]) + &pdb_text(&rows, 'A', &[4]);
    let result = read_pdb(&text);
    assert!(sequence_not_placed(&result));
    assert_eq!(pdb_state_at(&result, 2, "N"), (1, 3));
}

#[test]
fn caps_a_sequence_does_not_list_take_the_places_next_to_their_residues() {
    // An acetyl cap before three glycines and an N-methyl cap after them,
    // which SEQRES leaves out.
    let mut rows = vec![
        cap_site("C", "C", "ACE", 0, [-14.0, 0.0, 0.0]),
        cap_site("O", "O", "ACE", 0, [-14.0, 1.2, 0.0]),
        cap_site("CH3", "C", "ACE", 0, [-15.3, -0.7, 0.0]),
    ];
    rows.extend(separate_residues(&["GLY", "GLY", "GLY"], 1));
    rows.push(cap_site("N", "N", "NME", 4, [40.0, 0.0, 0.0]));
    rows.push(cap_site("C", "C", "NME", 4, [41.4, 0.0, 0.0]));
    let text = seqres('A', &["GLY"; 3]) + &pdb_text(&rows, 'A', &[4]);
    let text = text
        .lines()
        .map(|line| {
            if line.contains(" ACE ") || line.contains(" NME ") {
                line.replacen("ATOM  ", "HETATM", 1)
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let result = read_pdb(&text);
    assert!(!sequence_not_placed(&result), "{:?}", issues(&result));
    let places = result
        .structure()
        .topology()
        .chains()
        .flat_map(|chain| chain.residues())
        .map(|residue| (residue.name().to_owned(), residue.label_seq_id()))
        .collect::<Vec<_>>();
    assert_eq!(
        places,
        [("ACE", 1), ("GLY", 2), ("GLY", 3), ("GLY", 4), ("NME", 5)]
            .map(|(name, place)| (name.to_owned(), Some(place)))
    );
}

fn cap_site(name: &str, element: &str, component: &str, sequence: i32, point: [f64; 3]) -> Site {
    Site {
        element: element.to_owned(),
        name: name.to_owned(),
        component: component.to_owned(),
        asym: "A".to_owned(),
        sequence,
        point,
        occupancy: None,
        charge: None,
    }
}

#[test]
fn a_break_inside_a_chain_without_a_sequence_is_a_gap() {
    // Glycines 1 and 5 of chain A without hydrogens or SEQRES, the first
    // without OXT, in one run of records: the break between them is a gap,
    // so the N of the second is an amide.
    let rows = [
        separate_residues(&["GLY"], 1),
        heavy(residue_sites("GLY", "A", 5, [12.0, 0.0, 0.0])),
    ]
    .concat();
    let gap = read_pdb(&pdb_text(&rows, 'A', &[5]));
    assert_eq!(pdb_state_at(&gap, 5, "N"), (0, 1));
    assert!(has_issue(&gap, |issue| matches!(
        issue,
        Issue::PolymerLinkMissing { atom_name, .. } if atom_name == "C"
    )));

    // A TER record between them ends the chain there: each glycine is a
    // molecule of its own, with a free amine.
    let ended = read_pdb(&pdb_text(&rows, 'A', &[1, 5]));
    assert_eq!(pdb_state_at(&ended, 5, "N"), (0, 2));
    assert!(!has_issue(&ended, |issue| matches!(
        issue,
        Issue::PolymerLinkMissing { .. }
    )));

    // So does OXT before the break, whatever the residue after it shows.
    let rows = [
        heavy(residue_sites("GLY", "A", 1, [0.0; 3])),
        heavy(residue_sites("GLY", "A", 5, [12.0, 0.0, 0.0])),
    ]
    .concat();
    let ended = read_pdb(&pdb_text(&rows, 'A', &[5]));
    assert_eq!(pdb_state_at(&ended, 5, "N"), (1, 3));
    assert!(!has_issue(&ended, |issue| matches!(
        issue,
        Issue::PolymerLinkMissing { .. }
    )));
}

#[test]
fn chains_sharing_a_letter_get_no_merged_sequence() {
    // Two chains named A, each a dipeptide whose first N is an open amide:
    // one SEQRES sequence for the letter would join theirs, so none is
    // written. A reader would make the amides ammoniums: writing them
    // fails, or reports them under Omit.
    let sequenced = |asym: &str, offset: f64| {
        let mut rows = separate_residues(&["GLY", "ALA"], 2);
        for site in &mut rows {
            site.asym = asym.to_owned();
            site.point[1] += offset;
        }
        rows
    };
    let rows = [sequenced("A", 0.0), sequenced("B", 30.0)].concat();
    let mut text = with_entity(&block(&rows, &entity_sequence(&["SER", "GLY", "ALA"])));
    text = text.replace(
        "_atom_site.label_asym_id\n",
        "_atom_site.label_asym_id\n_atom_site.auth_asym_id\n",
    );
    text = text
        .lines()
        .map(|line| {
            if !line.starts_with("ATOM ") {
                return line.to_owned();
            }
            let mut fields = line.split(' ').collect::<Vec<_>>();
            fields.insert(6, "A");
            fields.join(" ")
        })
        .collect::<Vec<_>>()
        .join("\n");
    let result = read_mmcif(&text);
    assert_eq!(state_at(&result, "A", 2, "N"), (0, 1));
    assert!(kekule::pdb::write(&result, Default::default()).is_err());
    let omit = kekule::diagnostics::UnrepresentablePolicy::Omit;
    let options = kekule::pdb::PdbWriteOptions::default()
        .with_implicit_hydrogens(omit)
        .with_chemistry(omit);
    let (written, report) = kekule::pdb::write(&result, options).unwrap();
    assert!(!written.contains("SEQRES"), "{written}");
    assert!(
        report
            .omissions()
            .iter()
            .any(|omission| matches!(omission, kekule::diagnostics::WriteOmission::ImplicitHydrogens { atoms } if *atoms > 0)),
        "{:?}",
        report.omissions()
    );
}

const TRNA: &str = include_str!("../fixtures/mmcif/1EHZ.cif");
const CRAMBIN: &str = include_str!("../fixtures/pdb/1CRN.pdb");

fn read_pdb(text: &str) -> StructureInterpretation {
    kekule::pdb::parse_str(text, Default::default())
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

fn links_missing(result: &StructureInterpretation) -> Vec<Issue> {
    result
        .report()
        .issues()
        .iter()
        .filter(|issue| matches!(issue, Issue::PolymerLinkMissing { .. }))
        .cloned()
        .collect()
}

#[test]
fn a_three_prime_gap_round_trips() {
    // 1EHZ without A76, the last of its 76 nucleotides: C75's O3' links to
    // the absent residue, so it keeps no hydrogen.
    let text = TRNA
        .lines()
        .filter(|line| !(line.starts_with("ATOM") && line.contains(" A   A 1 76 ")))
        .collect::<Vec<_>>()
        .join(
            "
",
        );
    let source = read_mmcif(&text);
    let open = [link_missing("A", 75, "C", "O3'")];
    assert_eq!(links_missing(&source), open);
    assert_eq!(state_at(&source, "A", 75, "O3'"), (0, 0));
    // mmCIF states the sequence, from which readers find the gap again.
    let (text, _) = mmcif::write([&source], Default::default()).unwrap();
    assert!(text.contains("_entity_poly_seq.mon_id"), "{text}");
    let reread = read_mmcif(&text);
    assert_eq!(links_missing(&reread), open);
    assert_eq!(state_at(&reread, "A", 75, "O3'"), (0, 0));
    // So does PDB, with `SEQRES`.
    let (text, _) = kekule::pdb::write(&source, Default::default()).unwrap();
    assert!(
        text.starts_with("SEQRES   1 A   76    G   C   G   G"),
        "{text}"
    );
    let reread = read_pdb(&text);
    assert_eq!(links_missing(&reread), open);
    assert_eq!(state_at(&reread, "A", 75, "O3'"), (0, 0));
    assert_eq!(
        kekule::pdb::write(&reread, Default::default()).unwrap().0,
        text
    );
}

#[test]
fn a_chain_break_inside_a_sequenced_chain_round_trips_through_pdb() {
    // Alanines 1 and 2 of a sequence of two, both modelled but 10 Å apart:
    // the sequence opens the C of the first and the N of the second.
    let mut rows = without(heavy(sites("ALA")), &["OXT"]);
    rows.extend(heavy(residue_sites("ALA", "A", 2, [10.0, 0.0, 0.0])));
    let source = read_sequenced(&rows, &["ALA", "ALA"]);
    let open = [
        link_missing("A", 1, "ALA", "C"),
        link_missing("A", 2, "ALA", "N"),
    ];
    assert_eq!(links_missing(&source), open);
    assert_eq!(state_at(&source, "A", 2, "N"), (0, 1));
    // Only `SEQRES` tells a PDB reader that the chain continues.
    let (text, _) = kekule::pdb::write(&source, Default::default()).unwrap();
    assert!(
        text.starts_with(
            "SEQRES   1 A    2  ALA ALA
"
        ),
        "{text}"
    );
    let reread = read_pdb(&text);
    assert_eq!(links_missing(&reread), open);
    assert_eq!(state_at(&reread, "A", 1, "C"), (0, 0));
    assert_eq!(state_at(&reread, "A", 2, "N"), (0, 1));
    assert_eq!(
        kekule::pdb::write(&reread, Default::default()).unwrap().0,
        text
    );
}

#[test]
fn a_sequenced_heterogen_a_reader_cannot_place_follows_the_ter() {
    // Ethylene glycol takes place 3 of the sequence. Its `HETATM` records
    // have no linking atom and no polymer template, so a PDB reader takes
    // it as no polymer residue: the `TER` comes before it.
    let mut rows = without(heavy(sites("ALA")), &["OXT"]);
    let carbon = point(&rows, "C");
    let next = heavy(residue_sites("ALA", "A", 2, [0.0; 3]));
    let nitrogen = point(&next, "N");
    let shift = [0, 1, 2].map(|axis| carbon[axis] - nitrogen[axis] + [1.33, 0.0, 0.0][axis]);
    rows.extend(heavy(residue_sites("ALA", "A", 2, shift)));
    rows.extend(heavy(residue_sites("EDO", "A", 3, [0.0, 10.0, 0.0])));
    let source = read_sequenced(&rows, &["ALA", "ALA", "EDO"]);
    let (text, _) = kekule::pdb::write(&source, Default::default()).unwrap();
    let ter = text
        .lines()
        .position(|line| line.starts_with("TER"))
        .unwrap();
    let glycol = text
        .lines()
        .position(|line| line.contains(" EDO "))
        .unwrap();
    assert!(ter < glycol, "{text}");
    assert!(
        text.starts_with(
            "SEQRES   1 A    3  ALA ALA EDO
"
        ),
        "{text}"
    );
    let reread = read_pdb(&text);
    assert_eq!(
        kekule::pdb::write(&reread, Default::default()).unwrap().0,
        text
    );
}

/// Whether a write report lists charges readers do not give back.
fn omits_charges(report: &kekule::diagnostics::WriteReport) -> bool {
    report.omissions().iter().any(|omission| {
        matches!(
            omission,
            kekule::diagnostics::WriteOmission::FormalCharges { .. }
        )
    })
}

#[test]
fn written_charges_read_back_as_they_are() {
    // Entries as RCSB writes them, without hydrogens.
    for source in [read_pdb(CRAMBIN), read_mmcif(TRNA)] {
        let (text, report) = kekule::pdb::write(&source, Default::default()).unwrap();
        assert!(!omits_charges(&report), "{report:?}");
        let (_, report) = mmcif::write([&source], Default::default()).unwrap();
        assert!(!omits_charges(&report), "{report:?}");
        assert_eq!(sorted_charges(&read_pdb(&text)), sorted_charges(&source));
    }
    // A model with every hydrogen, as Maestro writes it: lysine 1 charged
    // with three hydrogens on NZ, lysine 2 neutral with two.
    let mut rows = residue_sites("LYS", "A", 1, [0.0; 3]);
    rows.extend(without(
        residue_sites("LYS", "B", 1, [10.0, 0.0, 0.0]),
        &["HZ3"],
    ));
    let source = read(&rows, "");
    assert_eq!(state_at(&source, "A", 1, "NZ").0, 1);
    assert_eq!(state_at(&source, "B", 1, "NZ").0, 0);
    let (text, report) = kekule::pdb::write(&source, Default::default()).unwrap();
    assert!(!omits_charges(&report), "{report:?}");
    assert_eq!(sorted_charges(&read_pdb(&text)), sorted_charges(&source));
    let (text, report) = mmcif::write([&source], Default::default()).unwrap();
    assert!(!omits_charges(&report), "{report:?}");
    assert_eq!(sorted_charges(&read_mmcif(&text)), sorted_charges(&source));
}

/// Nonzero formal charges as `(atom name, charge)`, sorted: PDB chains are
/// not mmCIF asyms.
fn sorted_charges(result: &StructureInterpretation) -> Vec<(String, i8)> {
    let mut charges = charges(result)
        .into_iter()
        .map(|(_, name, charge)| (name, charge))
        .collect::<Vec<_>>();
    charges.sort();
    charges
}

#[test]
fn a_metal_bound_amine_keeps_its_hydrogens() {
    // A38's N6 coordinates a magnesium in 1EHZ. An amine with a neutral
    // template does not give its hydrogen to the metal.
    let source = read_mmcif(TRNA);
    assert_eq!(state_at(&source, "A", 38, "N6"), (0, 2));
    assert!(!source.report().issues().iter().any(|issue| matches!(
        issue,
        Issue::HydrogenReplacedByCoordination { atoms, residue: at, .. }
            if atoms.contains(&"N6".to_owned()) && *at == residue("A", 38)
    )));
}

/// Crambin with a C-terminal cap of `cap`'s atoms in place of OXT, before
/// the chain's `TER`: each atom's name, element, and distance from the C.
fn capped_crambin(cap: &str, atoms: &[(&str, &str, f64)]) -> String {
    let carbon = [12.4, 5.982, 9.869];
    let oxt = [12.703, 4.973, 10.746];
    let direction = [0, 1, 2].map(|axis| oxt[axis] - carbon[axis]);
    let length = direction
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    let mut cap_records = String::new();
    for (index, &(name, element, distance)) in atoms.iter().enumerate() {
        let point = [0, 1, 2].map(|axis| carbon[axis] + distance * direction[axis] / length);
        cap_records.push_str(&format!(
            "HETATM{:>5}  {name:<3} {cap} A  47    {:>8.3}{:>8.3}{:>8.3}  1.00  7.86          {element:>2}\n",
            327 + index,
            point[0],
            point[1],
            point[2],
        ));
    }
    let mut text = String::new();
    for line in CRAMBIN.lines() {
        if line.starts_with("ATOM    327  OXT") {
            continue;
        }
        if line.starts_with("TER ") {
            text.push_str(&cap_records);
        }
        text.push_str(line);
        text.push('\n');
    }
    text
}

#[test]
fn a_capped_chain_ends_after_its_cap() {
    for (cap, atoms) in [
        ("NH2", &[("N", "N", 1.33)][..]),
        ("NME", &[("N", "N", 1.33), ("C", "C", 2.75)]),
    ] {
        let source = read_pdb(&capped_crambin(cap, atoms));
        let carbon = source
            .report()
            .instances()
            .iter()
            .flat_map(|instance| instance.atoms())
            .find(|atom| {
                atom.component_id() == "ASN"
                    && atom.atom_name() == "C"
                    && atom.author_sequence_id() == Some("46")
            })
            .unwrap()
            .atom();
        let molecule = source
            .structure()
            .topology()
            .instance(carbon.instance())
            .unwrap()
            .molecule();
        assert!(
            molecule.neighbors(carbon.atom()).unwrap().any(|other| {
                molecule.atom(other).unwrap().element() == kekule::molecule::Element::N
            }),
            "{cap}"
        );
        let (text, _) = kekule::pdb::write(&source, Default::default()).unwrap();
        // `TER` follows the cap, which stays part of the chain.
        let lines = text.lines().collect::<Vec<_>>();
        let ter = lines
            .iter()
            .position(|line| line.starts_with("TER"))
            .unwrap();
        assert!(lines[ter - 1].contains(cap), "{cap}: {}", lines[ter - 1]);
        let reread = read_pdb(&text);
        assert_eq!(
            kekule::pdb::write(&reread, Default::default()).unwrap().0,
            text,
            "{cap}"
        );
    }
}

#[test]
fn a_chain_without_its_first_residues_keeps_its_places() {
    // Crambin without residues 1 to 3: SEQRES places the rest at 4 to 46,
    // in PDB and mmCIF alike.
    let text = CRAMBIN
        .lines()
        .filter(|line| {
            let record = line.get(..6).unwrap_or("");
            let number = line
                .get(22..26)
                .and_then(|text| text.trim().parse::<i32>().ok());
            !(matches!(record, "ATOM  " | "HETATM") && number.is_some_and(|number| number <= 3))
                && !line.starts_with("SSBOND   1")
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let places = |result: &StructureInterpretation| {
        result
            .report()
            .instances()
            .iter()
            .flat_map(|instance| instance.atoms())
            .filter(|atom| atom.atom_name() == "CA")
            .map(|atom| atom.label_sequence_id())
            .collect::<Vec<_>>()
    };
    let expected = (4..=46).map(Some).collect::<Vec<_>>();
    let source = read_pdb(&text);
    assert_eq!(places(&source), expected);
    let (written, _) = kekule::pdb::write(&source, Default::default()).unwrap();
    assert!(
        written.starts_with("SEQRES   1 A   46  THR THR CYS CYS PRO SER"),
        "{written}"
    );
    let reread = read_pdb(&written);
    assert_eq!(places(&reread), expected);
    assert_eq!(
        kekule::pdb::write(&reread, Default::default()).unwrap().0,
        written
    );
    let (written, _) = mmcif::write([&source], Default::default()).unwrap();
    assert_eq!(places(&read_mmcif(&written)), expected);
}
