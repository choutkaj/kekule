//! Non-covalent interactions of PDB complexes, checked against their
//! geometry: each expected distance and angle below was measured
//! independently from the entry's atom_site coordinates (model 1, first
//! alternate location), and the interactions are the ones the primary
//! publications describe.

use std::sync::Arc;

use kekule::interactions::{
    find_interactions_between, find_interactions_within, Interaction, InteractionCriteria,
    InteractionEnd, InteractionFeatures, InteractionKind, InteractionRole, InteractionSide,
};
use kekule::structure::Structure;
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex};
use kekule::units::{ANGSTROM, DEGREE};

fn load(entry: &str) -> Structure {
    let path = format!(
        "{}/tests/fixtures/mmcif/{entry}.cif",
        env!("CARGO_MANIFEST_DIR")
    );
    kekule::mmcif::parse_str(&std::fs::read_to_string(path).unwrap(), Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_structure()
}

/// `RES<seq>:<atom>` by author residue number and label atom name.
fn label(topology: &Topology, atom: TopologyAtomIndex) -> String {
    let site = topology.atom_at(atom).unwrap().atom_site().unwrap();
    let residue = site.residue();
    format!(
        "{}{}:{}",
        residue.name(),
        residue.author_seq_id().unwrap_or_default(),
        site.metadata().label_atom_id().unwrap_or_default()
    )
}

fn labels(topology: &Topology, end: &InteractionEnd) -> String {
    end.atoms()
        .iter()
        .map(|&atom| label(topology, atom))
        .collect::<Vec<_>>()
        .join(",")
}

type Found = Vec<(InteractionKind, String, String, Interaction)>;

/// The interactions of one residue's atoms with every atom, in a structure
/// whose every molecule perceives.
fn of_ligand(structure: &Structure, ligand: &str) -> Found {
    let topology = structure.shared_topology();
    let features = InteractionFeatures::new(Arc::clone(&topology));
    let unperceived = features
        .unperceived_definitions()
        .iter()
        .map(|&definition| {
            let instance = topology.instances_of(definition).next().unwrap();
            let atom = instance.atoms().next().unwrap();
            atom.residue().unwrap().name().to_owned()
        })
        .collect::<Vec<_>>();
    assert!(unperceived.is_empty(), "{unperceived:?}");
    let ligand = AtomSelection::from_predicate(Arc::clone(&topology), |atom| {
        atom.residue()
            .is_some_and(|residue| residue.name() == ligand)
    });
    let all = AtomSelection::all(Arc::clone(&topology));
    let found = find_interactions_within(
        &features,
        structure,
        &ligand,
        &all,
        &InteractionCriteria::default(),
    )
    .unwrap();
    found
        .interactions()
        .iter()
        .map(|interaction| {
            (
                interaction.kind(),
                labels(&topology, interaction.first()),
                labels(&topology, interaction.second()),
                interaction.clone(),
            )
        })
        .collect()
}

/// The one interaction of `kind` between ends labelled `first` and
/// `second`, at `distance` Å.
fn expect<'a>(
    found: &'a Found,
    kind: InteractionKind,
    first: &str,
    second: &str,
    distance: f64,
) -> &'a Interaction {
    let matches = found
        .iter()
        .filter(|(k, a, b, _)| *k == kind && a == first && b == second)
        .collect::<Vec<_>>();
    assert_eq!(
        matches.len(),
        1,
        "{kind} {first} -- {second} in {:#?}",
        summary(found)
    );
    let interaction = &matches[0].3;
    let measured = interaction.distance().value_in(ANGSTROM).unwrap();
    assert!(
        (measured - distance).abs() < 5e-4,
        "{kind} {first} -- {second}: {measured}"
    );
    interaction
}

fn summary(found: &[(InteractionKind, String, String, Interaction)]) -> Vec<String> {
    found
        .iter()
        .filter(|(kind, ..)| *kind != InteractionKind::HydrophobicContact)
        .map(|(kind, a, b, _)| format!("{kind}: {a} -- {b}"))
        .collect()
}

fn angle(interaction: &Interaction, index: usize) -> f64 {
    interaction
        .angles()
        .nth(index)
        .unwrap()
        .1
        .value_in(DEGREE)
        .unwrap()
}

#[test]
fn kras_g12c_inhibitor_9uoh() {
    // KRas G12C with GDP and the switch-II pocket inhibitor ASP2453
    // (A1L9E): hydrogen bonds to His95 and Lys16, stacking on Tyr96, and a
    // salt bridge of its amine to a sulfate.
    let found = of_ligand(&load("9UOH"), "A1L9E");
    use InteractionKind::*;
    expect(&found, HydrogenBond, "A1L9E203:N15", "HIS95:NE2", 2.8476);
    let lysine = expect(&found, HydrogenBond, "A1L9E203:O43", "LYS16:NZ", 2.6464);
    // CE–NZ…O at the donor, which has no explicit hydrogens.
    assert!((angle(lysine, 0) - 102.017).abs() < 1e-3);
    assert_eq!(lysine.second().role(), InteractionRole::Donor);
    assert!(lysine.second().hydrogen().is_none());
    let tyrosine = expect(
        &found,
        PiStackingTShaped,
        "A1L9E203:C20,A1L9E203:C19,A1L9E203:C18,A1L9E203:C17,A1L9E203:C16,A1L9E203:C21",
        "TYR96:CD1,TYR96:CE1,TYR96:CZ,TYR96:CE2,TYR96:CD2,TYR96:CG",
        5.2868,
    );
    assert!((angle(tyrosine, 0) - 67.569).abs() < 1e-3);
    expect(
        &found,
        SaltBridge,
        "A1L9E203:N37",
        "SO4204:O1,SO4204:O2,SO4204:O3,SO4204:O4,SO4204:S",
        3.8275,
    );
    // Every reported interaction meets its criteria; no halogen bonds (its
    // halogens are fluorines) and no metal near the ligand.
    assert!(found
        .iter()
        .all(|(kind, ..)| !matches!(kind, HalogenBond | MetalCoordination)));
}

#[test]
fn trypsin_benzamidine_3ptb() {
    // Benzamidine's amidinium bridges Asp189 at the bottom of the S1
    // pocket; the calcium site binds Glu70, Asn72, Val75, Glu80, and two
    // waters.
    let structure = load("3PTB");
    let found = of_ligand(&structure, "BEN");
    use InteractionKind::*;
    expect(
        &found,
        SaltBridge,
        "BEN1:N1,BEN1:N2,BEN1:C",
        "ASP189:OD1,ASP189:OD2,ASP189:CG",
        2.8943,
    );
    expect(&found, HydrogenBond, "BEN1:N1", "ASP189:OD1", 2.9188);
    expect(&found, HydrogenBond, "BEN1:N2", "ASP189:OD2", 2.8725);
    expect(&found, HydrogenBond, "BEN1:N1", "GLY219:O", 2.8153);
    expect(&found, HydrogenBond, "BEN1:N2", "SER190:OG", 3.0417);

    let calcium = of_ligand(&structure, "CA");
    let ligands = calcium
        .iter()
        .filter(|(kind, ..)| *kind == MetalCoordination)
        .map(|(_, _, ligand, _)| ligand.split(':').next().unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(ligands.len(), 6, "{ligands:?}");
    for residue in ["GLU70", "ASN72", "VAL75", "GLU80"] {
        assert!(
            ligands.iter().any(|ligand| ligand == residue),
            "{residue} in {ligands:?}"
        );
    }
    expect(&calcium, MetalCoordination, "CA480:CA", "GLU70:OE1", 2.3670);
}

#[test]
fn acetylcholinesterase_tacrine_1acj() {
    // Tacrine stacks between Trp84 and Phe330 in the gorge.
    let found = of_ligand(&load("1ACJ"), "THA");
    let phenylalanine = expect(
        &found,
        InteractionKind::PiStackingParallel,
        "THA999:C1,THA999:C6,THA999:C5,THA999:C4,THA999:C3,THA999:C2",
        "PHE330:CD1,PHE330:CE1,PHE330:CZ,PHE330:CE2,PHE330:CD2,PHE330:CG",
        3.7682,
    );
    assert!((angle(phenylalanine, 0) - 22.737).abs() < 1e-3);
    let offset = phenylalanine.offset().unwrap().value_in(ANGSTROM).unwrap();
    assert!((offset - 0.2602).abs() < 1e-4);
    assert!(found.iter().any(|(kind, _, ring, _)| {
        *kind == InteractionKind::PiStackingParallel && ring.starts_with("TRP84:")
    }));
}

#[test]
fn hydrophobic_contacts_are_reduced_as_plip_reduces_them() {
    // Tacrine stacks between Trp84 and Phe330, and its carbons touch the
    // gorge's apolar side chains.
    let found = of_ligand(&load("1ACJ"), "THA");
    let residue = |label: &str| label.split(':').next().unwrap().to_owned();
    let contacts = found
        .iter()
        .filter(|(kind, ..)| *kind == InteractionKind::HydrophobicContact)
        .map(|(_, ligand, site, _)| (ligand.clone(), site.clone()))
        .collect::<Vec<_>>();
    assert!(contacts.len() > 1, "{contacts:?}");
    // One contact per ligand atom and residue, and per site atom and residue.
    let mut per_ligand_atom = std::collections::BTreeSet::new();
    let mut per_site_atom = std::collections::BTreeSet::new();
    for (ligand, site) in &contacts {
        assert!(
            per_ligand_atom.insert((ligand, residue(site))),
            "{ligand} -- {site}"
        );
        assert!(
            per_site_atom.insert((site, residue(ligand))),
            "{ligand} -- {site}"
        );
    }
    // None between the atoms of two stacked rings.
    let stacked = found
        .iter()
        .filter(|(kind, ..)| {
            matches!(
                kind,
                InteractionKind::PiStackingParallel | InteractionKind::PiStackingTShaped
            )
        })
        .collect::<Vec<_>>();
    assert!(!stacked.is_empty());
    for (_, first, second, _) in stacked {
        let first = first.split(',').collect::<Vec<_>>();
        let second = second.split(',').collect::<Vec<_>>();
        assert!(!contacts.iter().any(|(ligand, site)| {
            first.contains(&ligand.as_str()) && second.contains(&site.as_str())
        }));
    }
}

#[test]
fn hydrophobic_contacts_do_not_depend_on_reporting_pi_stacking() {
    // Tacrine's rings stack on Trp84 and Phe330. A search for hydrophobic
    // contacts alone still drops the contacts between those rings.
    let structure = load("1ACJ");
    let topology = structure.shared_topology();
    let features = InteractionFeatures::new(Arc::clone(&topology));
    let tacrine = AtomSelection::from_predicate(Arc::clone(&topology), |atom| {
        atom.residue()
            .is_some_and(|residue| residue.name() == "THA")
    });
    let all = AtomSelection::all(Arc::clone(&topology));
    let defaults = InteractionCriteria::default();
    let alone = defaults
        .with_hydrogen_bond(defaults.hydrogen_bond().with_enabled(false))
        .with_salt_bridge(defaults.salt_bridge().with_enabled(false))
        .with_pi_stacking(defaults.pi_stacking().with_enabled(false))
        .with_cation_pi(defaults.cation_pi().with_enabled(false))
        .with_halogen_bond(defaults.halogen_bond().with_enabled(false))
        .with_metal_coordination(defaults.metal_coordination().with_enabled(false));
    for (first, second) in [(&tacrine, &all), (&all, &tacrine), (&all, &all)] {
        let search = |criteria: &InteractionCriteria| {
            find_interactions_within(&features, &structure, first, second, criteria).unwrap()
        };
        let reported = search(&defaults);
        assert!(reported.interactions().iter().any(|interaction| matches!(
            interaction.kind(),
            InteractionKind::PiStackingParallel | InteractionKind::PiStackingTShaped
        )));
        let contacts = |found: &kekule::interactions::Interactions| {
            found
                .of_kind(InteractionKind::HydrophobicContact)
                .cloned()
                .collect::<Vec<_>>()
        };
        let found = search(&alone);
        assert!(!contacts(&reported).is_empty());
        assert_eq!(contacts(&found), contacts(&reported));
        assert_eq!(found.interactions(), contacts(&reported).as_slice());
    }
}

#[test]
fn acetylcholinesterase_decamethonium_1acl() {
    // Decamethonium's quaternary ammoniums sit over Trp84 and Trp279. The
    // entry states no charges, so its four-bonded nitrogens take the +1
    // their bonds imply.
    let found = of_ligand(&load("1ACL"), "DME");
    let cation = expect(
        &found,
        InteractionKind::CationPi,
        "DME999:N12",
        "TRP84:CE3,TRP84:CZ3,TRP84:CH2,TRP84:CZ2,TRP84:CE2,TRP84:CD2",
        4.2925,
    );
    let offset = cation.offset().unwrap().value_in(ANGSTROM).unwrap();
    assert!((offset - 0.9205).abs() < 1e-4);
    assert_eq!(cation.first().role(), InteractionRole::Cation);
    assert!(found.iter().any(|(kind, cation, ring, _)| {
        *kind == InteractionKind::CationPi && cation == "DME999:N1" && ring.starts_with("TRP279:")
    }));
}

#[test]
fn egfr_gefitinib_2ito() {
    // Gefitinib's chlorine halogen-bonds the Leu788 carbonyl; its
    // quinazoline accepts from the Met793 hinge.
    let found = of_ligand(&load("2ITO"), "IRE");
    let halogen = expect(
        &found,
        InteractionKind::HalogenBond,
        "IRE2020:CL",
        "LEU788:O",
        3.2850,
    );
    // X…O=C at the acceptor.
    assert!((angle(halogen, 1) - 108.712).abs() < 1e-3);
    expect(
        &found,
        InteractionKind::HydrogenBond,
        "IRE2020:N3",
        "MET793:N",
        2.9164,
    );
}

#[test]
fn rubredoxin_iron_1iro() {
    // Iron bound by four cysteine sulfurs, which the entry also links.
    let found = of_ligand(&load("1IRO"), "FE");
    for (cysteine, distance) in [
        ("CYS6", 2.288),
        ("CYS9", 2.249),
        ("CYS39", 2.283),
        ("CYS42", 2.229),
    ] {
        expect(
            &found,
            InteractionKind::MetalCoordination,
            "FE55:FE",
            &format!("{cysteine}:SG"),
            distance,
        );
    }
}

#[test]
fn split_structures_give_the_same_interactions() {
    // 9UOH's inhibitor and everything else, as one structure and as two.
    let structure = load("9UOH");
    let topology = structure.shared_topology();
    let features = InteractionFeatures::new(Arc::clone(&topology));
    let ligand = AtomSelection::from_predicate(Arc::clone(&topology), |atom| {
        atom.residue()
            .is_some_and(|residue| residue.name() == "A1L9E")
    });
    let rest = ligand.complement();
    let criteria = InteractionCriteria::default();
    let within =
        find_interactions_within(&features, &structure, &ligand, &rest, &criteria).unwrap();

    let ligand_structure = structure.subset(&ligand).unwrap();
    let rest_structure = structure.subset(&rest).unwrap();
    let (ligand_features, rest_features) = (
        InteractionFeatures::new(ligand_structure.shared_topology()),
        InteractionFeatures::new(rest_structure.shared_topology()),
    );
    let (ligand_atoms, rest_atoms) = (
        AtomSelection::all(ligand_structure.shared_topology()),
        AtomSelection::all(rest_structure.shared_topology()),
    );
    let between = find_interactions_between(
        &InteractionSide::new(&ligand_features, &ligand_structure, &ligand_atoms).unwrap(),
        &InteractionSide::new(&rest_features, &rest_structure, &rest_atoms).unwrap(),
        &criteria,
    )
    .unwrap();

    // Subsets keep source order, so a subset's atom i is the i-th selected.
    let source = |selection: &AtomSelection, atoms: &[TopologyAtomIndex]| {
        atoms
            .iter()
            .map(|atom| selection.indices()[atom.index()])
            .collect::<Vec<_>>()
    };
    let key =
        |interaction: &Interaction, first: &[TopologyAtomIndex], second: &[TopologyAtomIndex]| {
            (
                interaction.kind(),
                first.to_vec(),
                second.to_vec(),
                interaction.distance().value_in(ANGSTROM).unwrap(),
                interaction
                    .angles()
                    .map(|(kind, value)| (kind, value.value_in(DEGREE).unwrap()))
                    .collect::<Vec<_>>(),
            )
        };
    let mut expected = within
        .interactions()
        .iter()
        .map(|interaction| {
            key(
                interaction,
                interaction.first().atoms(),
                interaction.second().atoms(),
            )
        })
        .collect::<Vec<_>>();
    let mut actual = between
        .interactions()
        .iter()
        .map(|interaction| {
            key(
                interaction,
                &source(&ligand, interaction.first().atoms()),
                &source(&rest, interaction.second().atoms()),
            )
        })
        .collect::<Vec<_>>();
    let order = |a: &(
        InteractionKind,
        Vec<TopologyAtomIndex>,
        Vec<TopologyAtomIndex>,
        f64,
        Vec<_>,
    ),
                 b: &(
        InteractionKind,
        Vec<TopologyAtomIndex>,
        Vec<TopologyAtomIndex>,
        f64,
        Vec<_>,
    )| { (a.0, &a.1, &a.2).cmp(&(b.0, &b.1, &b.2)) };
    expected.sort_by(order);
    actual.sort_by(order);
    assert!(expected.len() > 10, "{}", expected.len());
    assert_eq!(expected, actual);
}
