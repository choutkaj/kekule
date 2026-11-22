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
    text.lines()
        .filter_map(|line| {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() < 18 || fields[0] != component {
                return None;
            }
            let coordinate = |index: usize| fields[index].parse::<f64>().ok();
            let point = [
                coordinate(15)? + offset[0],
                coordinate(16)? + offset[1],
                coordinate(17)? + offset[2],
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
            site.sequence,
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
fn stated_charges_are_checked_against_bonds_and_hydrogens() {
    // A stated +1 on a carboxylate oxygen would make an oxonium ion: the
    // template's charge is kept, as for arginine's NH1, whose stated +1
    // the guanidinium's NH2 already carries.
    let mut rows = heavy(sites("GLU"));
    state_charge(&mut rows, "OE2", 1);
    let mut arginine = heavy(residue_sites("ARG", "B", 1, [20.0, 0.0, 0.0]));
    state_charge(&mut arginine, "NH1", 1);
    rows.extend(arginine);
    let result = read(&rows, "");
    let found = charges(&result);
    assert!(
        !found.contains(&("A".to_owned(), "OE2".to_owned(), 1)),
        "{found:?}"
    );
    assert!(
        !found.contains(&("B".to_owned(), "NH1".to_owned(), 1)),
        "{found:?}"
    );
    assert!(
        found.contains(&("B".to_owned(), "NH2".to_owned(), 1)),
        "{found:?}"
    );
    for (asym, atom) in [("A", "OE2"), ("B", "NH1")] {
        assert!(
            has_issue(&result, |issue| matches!(
                issue,
                Issue::StatedChargeRejected { residue: rejected, atom_name, stated: 1, .. }
                    if *rejected == residue(asym, 1) && atom_name == atom
            )),
            "{asym} {atom}"
        );
    }
    // A stated charge the atom can hold is kept, and a difference reported.
    let mut lysine = heavy(sites("LYS"));
    state_charge(&mut lysine, "NZ", 0);
    let neutral = read(&lysine, "");
    assert!(!charges(&neutral).contains(&("A".to_owned(), "NZ".to_owned(), 1)));
    assert!(has_issue(&neutral, |issue| matches!(
        issue,
        Issue::StatedChargeKept { atom_name, stated: 0, derived: 1, .. } if atom_name == "NZ"
    )));
    // A stated charge equal to the template's is no issue, and no
    // charge is said to come from hydrogens when the source states it.
    let mut aspartate = without(sites("ASP"), &["HD2"]);
    state_charge(&mut aspartate, "OD2", -1);
    let stated = read(&aspartate, "");
    assert_eq!(charges(&stated), charged(&[("A", "OD2", -1)]));
    assert_eq!(issues(&stated), []);
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
    // nitrogen's four bonds give it +1 and its single-bonded oxygens -1,
    // not hydrogens.
    let bonds_only = &tables[tables.find("loop_\n_chem_comp_bond").unwrap()..];
    let result = read(&nitrate, bonds_only);
    assert_eq!(
        charges(&result),
        charged(&[("A", "N", 1), ("A", "O2", -1), ("A", "O3", -1)])
    );
    assert_eq!(hydrogens(&result, "A", "O3"), 0);

    // A single atom without template or stated charge is reported.
    let caesium = Site {
        element: "CS".to_owned(),
        name: "CS".to_owned(),
        component: "CS".to_owned(),
        asym: "B".to_owned(),
        sequence: 1,
        point: [0.0; 3],
        occupancy: None,
        charge: None,
    };
    let result = read(&[caesium], "");
    assert!(
        has_issue(&result, |issue| matches!(
            issue,
            Issue::TemplateMissing { component_id, .. } if component_id == "CS"
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
