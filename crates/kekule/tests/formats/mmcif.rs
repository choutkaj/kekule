use kekule::diagnostics::SourceLocation;
use kekule::macromolecular::{CellOmissionReason, StructureInterpretIssue};
use kekule::mmcif;

const MINIMAL_MMCIF: &str = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
GLY C1 C2 sing
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
C C1 GLY A 1 1 0.0 0.0 0.0
C C2 GLY A 1 1 1.0 0.0 0.0
"#;

/// The cell of an interpreted structure as its diagonal in angstroms, when the
/// cell is orthorhombic, and the report's cell omissions.
fn orthorhombic_cell(text: &str) -> (Option<[f64; 3]>, Vec<StructureInterpretIssue>) {
    let interpreted = mmcif::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let omissions = interpreted
        .report()
        .issues()
        .iter()
        .filter(|issue| matches!(issue, StructureInterpretIssue::CellOmitted { .. }))
        .cloned()
        .collect();
    let cell = interpreted.structure().cell().map(|cell| {
        assert!(cell.is_orthorhombic());
        assert_eq!(cell.periodic_axes(), [true; 3]);
        let [a, b, c] = cell
            .vectors()
            .into_unit(kekule::units::ANGSTROM)
            .unwrap()
            .into_value();
        [a.x, b.y, c.z]
    });
    (cell, omissions)
}

#[test]
fn mmcif_cell_is_the_structure_cell_only_for_a_whole_unit_cell() {
    let cell = |lengths: &str| {
        let [a, b, c]: [&str; 3] = lengths.split(' ').collect::<Vec<_>>().try_into().unwrap();
        format!(
            "_cell.length_a {a}\n_cell.length_b {b}\n_cell.length_c {c}\n\
             _cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 90.0\n"
        )
    };
    // Omissions locate the space group, or else the first stated parameter.
    let omitted = |reason, line| {
        vec![StructureInterpretIssue::CellOmitted {
            reason,
            location: SourceLocation::at_line(line),
        }]
    };
    let space_group = |name: &str| CellOmissionReason::SpaceGroup {
        name: name.to_owned(),
    };
    for (header, expected, omissions) in [
        // No space group declared: the coordinates are taken as a unit cell.
        (cell("20 30 40.0(2)"), Some([20.0, 30.0, 40.0]), vec![]),
        (
            cell("20 30 40") + "_symmetry.space_group_name_H-M 'P 1'\n",
            Some([20.0, 30.0, 40.0]),
            vec![],
        ),
        // An asymmetric unit's translates are not its crystal.
        (
            cell("20 30 40") + "_symmetry.space_group_name_H-M 'P 21 21 21'\n",
            None,
            omitted(space_group("P 21 21 21"), 9),
        ),
        (
            cell("20 30 40") + "_space_group.name_H-M_alt 'P 21 21 21'\n",
            None,
            omitted(space_group("P 21 21 21"), 9),
        ),
        // PDB's placeholder for structures without a crystal lattice.
        (
            cell("1.000 1.000 1.000") + "_symmetry.space_group_name_H-M 'P 1'\n",
            None,
            omitted(CellOmissionReason::Placeholder, 3),
        ),
        (
            cell("20 30 ?"),
            None,
            omitted(CellOmissionReason::Incomplete, 3),
        ),
        // A `_cell` stating no lengths or angles omits nothing.
        ("_cell.entry_id demo\n".to_owned(), None, vec![]),
    ] {
        let text = MINIMAL_MMCIF.replacen("data_demo\n", &format!("data_demo\n{header}"), 1);
        let (actual, actual_omissions) = orthorhombic_cell(&text);
        assert_eq!(actual_omissions, omissions, "{header}");
        match (actual, expected) {
            (Some(actual), Some(expected)) => {
                for (actual, expected) in actual.into_iter().zip(expected) {
                    assert!((actual - expected).abs() < 1e-9, "{header}");
                }
            }
            (actual, expected) => assert_eq!(actual, expected, "{header}"),
        }
    }
    let malformed = MINIMAL_MMCIF.replacen(
        "data_demo\n",
        &format!("data_demo\n{}", cell("20 abc 40")),
        1,
    );
    let error = mmcif::parse_str(&malformed, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap_err();
    assert_eq!(
        error.message(),
        "The _cell.length_b value \"abc\" is not a number."
    );
}

fn with_struct_conn(order: &str) -> String {
    format!(
        r#"{MINIMAL_MMCIF}
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_value_order
duplicate covale A C1 A C2 {order}
"#
    )
}

#[test]
fn mmcif_public_facade_requires_parse_then_interpret() -> Result<(), Box<dyn std::error::Error>> {
    use kekule::macromolecular::{self, EntityKind, StructureInterpretOptions};
    use kekule::mmcif::{
        self, MmcifBlock, MmcifEntityClassifications, MmcifWriteIssue, MmcifWriteOptions,
    };

    let document = mmcif::parse_str(MINIMAL_MMCIF, Default::default())?;
    assert_eq!(
        document,
        mmcif::parse_str(MINIMAL_MMCIF, mmcif::MmcifParseOptions::default())?
    );
    let interpreted = document.interpret(Default::default())?;
    let block: &MmcifBlock = &document.blocks()[0];
    let block_interpreted = block.interpret(Default::default())?;
    let method_default = block.interpret(Default::default())?;
    let method_explicit = block.interpret(StructureInterpretOptions::default())?;

    assert_eq!(document.blocks().len(), 1);
    assert_eq!(interpreted.report().block_name(), "demo");
    assert_eq!(block_interpreted.report(), interpreted.report());
    assert_eq!(method_default.report(), method_explicit.report());
    assert_eq!(
        document.interpret(Default::default())?.report(),
        method_default.report()
    );
    assert!(block_interpreted
        .topology()
        .same_layout(interpreted.topology()));
    assert_eq!(
        block_interpreted.structure().positions(),
        interpreted.structure().positions()
    );
    assert_eq!(interpreted.structure().topology().instance_count(), 1);
    assert_ne!(interpreted.structure().topology().atom_sites().len(), 0);
    assert_eq!(interpreted.structure().positions().len(), 2);
    assert_eq!(interpreted.molecules().count(), 1);
    assert!(interpreted
        .molecules()
        .all(|molecule| molecule.perception() == &kekule::perception::Perception::default()));
    assert!(interpreted
        .topology()
        .same_layout(interpreted.structure().topology()));
    let owned_topology = method_default.clone().into_topology();
    let owned_structure = method_default.clone().into_structure();
    let owned_molecules = method_default.into_molecules();
    assert!(owned_topology.same_layout(owned_structure.topology()));
    assert_eq!(
        owned_topology
            .instances()
            .map(|occurrence| occurrence.molecule())
            .collect::<Vec<_>>(),
        owned_molecules.iter().collect::<Vec<_>>()
    );
    assert!(block
        .interpret(Default::default())?
        .into_topology()
        .same_layout(
            block
                .interpret(Default::default())?
                .into_structure()
                .topology()
        ));
    assert_eq!(
        block.interpret(Default::default())?.into_molecules().len(),
        1
    );
    // Every format namespace pairs `interpret(document)` with
    // `interpret_with_options`, and a single-block document projects like
    // its block.
    let explicit = document.interpret(StructureInterpretOptions::default())?;
    assert_eq!(explicit.report(), interpreted.report());
    for from_document in [
        document.interpret(Default::default())?.into_structure(),
        document
            .interpret(StructureInterpretOptions::default())?
            .into_structure(),
    ] {
        let from_block = block.interpret(Default::default())?.into_structure();
        assert!(from_document.topology().same_layout(from_block.topology()));
        assert_eq!(from_document.positions(), from_block.positions());
    }
    assert_eq!(
        document.interpret(Default::default())?.into_molecules(),
        block.interpret(Default::default())?.into_molecules()
    );
    assert!(document
        .interpret(Default::default())?
        .into_topology()
        .same_layout(
            block
                .interpret(Default::default())?
                .into_structure()
                .topology()
        ));
    let set = document.interpret_structure_set(Default::default())?;
    assert_eq!(
        set.reports(),
        document
            .interpret_structure_set(macromolecular::StructureSetInterpretOptions::default())?
            .reports()
    );
    // The one model is the set's one geometry: two carbons 1 A apart.
    assert_eq!(set.structure_set().len(), 1);
    assert_eq!(set.structure_set().order(), None);
    assert_eq!(set.topology().atom_count(), 2);
    assert_eq!(set.topology().bond_count(), 1);
    let geometry = set.structure_set().get(0).unwrap();
    let positions = geometry.positions();
    let x = |index| {
        positions
            .get(index)
            .unwrap()
            .value_in(kekule::units::ANGSTROM)
            .unwrap()
    };
    assert_eq!(x(0).x, 0.0);
    assert!((x(1).x - 1.0).abs() <= 1e-12);
    assert_eq!((x(1).y, x(1).z), (0.0, 0.0));
    let two_blocks = mmcif::parse_str(
        &format!("{MINIMAL_MMCIF}{}", MINIMAL_MMCIF.replace("demo", "other")),
        Default::default(),
    )?;
    // The same failure, whichever projection follows interpretation.
    assert_eq!(
        format!(
            "{:?}",
            two_blocks
                .interpret(Default::default())
                .map(|interpretation| interpretation.into_structure())
                .unwrap_err()
        ),
        format!(
            "{:?}",
            two_blocks.interpret(Default::default()).unwrap_err()
        )
    );
    assert_eq!(interpreted.report().selected_model(), Some("1"));
    assert_eq!(interpreted.report().instances().len(), 1);
    let provenance = &interpreted.report().instances()[0];
    assert_eq!(provenance.coordinate_model_id(), "1");
    assert_eq!(provenance.atoms().len(), 2);
    assert_eq!(provenance.atoms()[0].atom_name(), "C1");
    assert_eq!(provenance.atoms()[0].label_sequence_id(), None);
    assert_eq!(provenance.atoms()[0].author_sequence_id(), Some("1"));
    assert_eq!(provenance.atoms()[0].insertion_code(), None);
    assert_eq!(provenance.atoms()[0].occurrence(), None);
    assert_eq!(provenance.atoms()[0].selected_alternate_location(), None);
    assert_eq!(
        provenance.atoms()[0].atom().instance(),
        provenance.instance()
    );
    let written = mmcif::write(
        [mmcif::MmcifBlockSource::structure(interpreted.structure())
            .with_reports(std::slice::from_ref(interpreted.report()))],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)?;
    assert!(written.starts_with("data_model\n"));
    assert!(mmcif::parse_str(&written, Default::default()).is_ok());
    let custom_written = mmcif::write(
        [mmcif::MmcifBlockSource::structure(interpreted.structure())
            .with_reports(std::slice::from_ref(interpreted.report()))],
        MmcifWriteOptions::default()
            .with_block_name("custom".to_owned())
            .unwrap(),
    )
    .map(|(text, _)| text)?;
    assert!(custom_written.starts_with("data_custom\n"));
    let unnamed = MmcifWriteOptions::default()
        .with_block_name(String::new())
        .unwrap_err();
    assert_eq!(
        unnamed.detail(),
        Some(&kekule::diagnostics::FormatErrorDetail::MmcifWrite(
            MmcifWriteIssue::InvalidBlockName(String::new())
        ))
    );
    let automatic = mmcif::write([interpreted.structure()], MmcifWriteOptions::default())
        .map(|(text, _)| text)?;
    assert!(automatic.contains("1 non-polymer"));
    let mut classifications = MmcifEntityClassifications::new();
    classifications.insert(provenance.instance(), EntityKind::Polymer)?;
    let written = mmcif::write(
        [mmcif::MmcifBlockSource::structure(interpreted.structure())
            .with_classifications(&classifications)],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)?;
    assert!(written.contains("1 polymer"));
    let structure = interpreted.into_structure();
    assert_eq!(structure.topology().instance_count(), 1);
    assert_eq!(structure.positions().len(), 2);

    Ok(())
}

#[test]
fn duplicate_agreeing_authoritative_bond_evidence_is_idempotent() {
    use kekule::mmcif;

    let document = mmcif::parse_str(&with_struct_conn("sing"), Default::default())
        .expect("duplicate agreeing connectivity parses");
    let interpreted = document
        .interpret(Default::default())
        .expect("duplicate agreeing connectivity is accepted");
    let graph = interpreted
        .structure()
        .topology()
        .definitions()
        .next()
        .expect("one molecule definition")
        .molecule();

    assert_eq!(graph.bond_count(), 1);
}

#[test]
fn duplicate_conflicting_authoritative_bond_evidence_is_rejected() {
    use kekule::diagnostics::FormatError;
    use kekule::mmcif;

    let document = mmcif::parse_str(&with_struct_conn("doub"), Default::default())
        .expect("duplicate conflicting connectivity parses");
    let error: FormatError = document
        .interpret(Default::default())
        .expect_err("conflicting authoritative bond orders must be rejected");

    assert_eq!(error.location(), None);
    assert_eq!(
        error.message(),
        "The bond records give one atom pair two orders: double and single."
    );
}

#[test]
fn mmcif_loop_accessors_reject_out_of_range_indices_in_all_build_profiles() {
    for columns in [1, 2, 3, 5] {
        let mut input = String::from("data_bounds\nloop_\n");
        for column in 0..columns {
            input.push_str(&format!("_x.c{column}\n"));
        }
        for column in 0..columns {
            input.push_str(&format!("v{column} "));
        }
        let document = mmcif::parse_str(&input, Default::default()).unwrap();
        let table = document.blocks()[0].loop_with_tag("_x.c0").unwrap();
        assert_eq!(table.row_count(), 1);
        assert_eq!(table.row(0).unwrap().len(), columns);
        for column in 0..columns {
            let tag = format!("_x.c{column}");
            assert_eq!(table.value(0, &tag).unwrap().text(), format!("v{column}"));
            for row in [
                1,
                usize::MAX / columns,
                (usize::MAX / columns).saturating_add(1),
                usize::MAX,
            ] {
                assert!(
                    table.value(row, &tag).is_none(),
                    "row {row}, columns {columns}"
                );
                assert!(table.row(row).is_none(), "row {row}, columns {columns}");
            }
        }
        assert!(table.value(0, "_x.missing").is_none());
    }
    let document = mmcif::parse_str("data_empty\nloop_\n_x.a\n", Default::default()).unwrap();
    let table = document.blocks()[0].loop_with_tag("_x.a").unwrap();
    assert_eq!(table.row_count(), 0);
    assert!(table.row(0).is_none());
    assert!(table.value(usize::MAX, "_x.a").is_none());
}

const CRAMBIN: &str = include_str!("../fixtures/mmcif/1CRN.cif");
const ADENYLATE_KINASE: &str = include_str!("../fixtures/mmcif/1AKE.cif");

/// The block without its `_chem_comp_bond` loop, as wwPDB archive entries are
/// distributed, so that every residue bond comes from a template or distance.
fn without_component_bonds(text: &str) -> String {
    let lines = text.lines().collect::<Vec<_>>();
    let mut kept = Vec::with_capacity(lines.len());
    let mut index = 0;
    while index < lines.len() {
        if lines[index].trim() == "loop_"
            && lines
                .get(index + 1)
                .is_some_and(|line| line.starts_with("_chem_comp_bond."))
        {
            index += 1;
            while lines[index].starts_with("_chem_comp_bond.") {
                index += 1;
            }
            while index < lines.len() && !lines[index].starts_with(['#', '_']) {
                index += 1;
            }
            continue;
        }
        kept.push(lines[index]);
        index += 1;
    }
    kept.join("\n")
}

#[test]
fn mmcif_connectivity_matches_an_independent_reference() {
    use kekule::macromolecular::{SourceResidueId, SourceResiduePosition, StructureInterpretIssue};

    // Atom, bond, and molecule counts from Biotite 1.2.0
    // (`pdbx.get_structure(include_bonds=True, altloc="occupancy")`), which
    // resolves residue bonds from `_chem_comp_bond` or its own copy of the
    // Chemical Component Dictionary; it gives the same counts with and without
    // the files' bond tables. Charged atoms count the LYS NZ, ARG NH2, and
    // HIS ND1 sites the dictionary charges (+1), counted from the files.
    for (name, text, atoms, bonds, largest, molecules, charged) in [
        ("1CRN", CRAMBIN, 327, 337, vec![327], 1, 2),
        (
            "1AKE",
            ADENYLATE_KINASE,
            3804,
            3484,
            vec![1656, 1656, 57, 57],
            382,
            68,
        ),
    ] {
        for stripped in [false, true] {
            let source = if stripped {
                without_component_bonds(text)
            } else {
                text.to_owned()
            };
            let interpreted = (mmcif::parse_str(&source, Default::default()).unwrap())
                .interpret(Default::default())
                .unwrap();
            let topology = interpreted.topology();
            let mut sizes = topology
                .instances()
                .map(|instance| instance.molecule().atom_count())
                .collect::<Vec<_>>();
            sizes.sort_unstable_by(|a, b| b.cmp(a));
            let charges = topology
                .instances()
                .flat_map(|instance| instance.molecule().atoms())
                .filter(|(_, atom)| atom.formal_charge() != 0)
                .map(|(_, atom)| atom.formal_charge())
                .collect::<Vec<_>>();
            let case = format!("{name} stripped={stripped}");
            assert_eq!(topology.atom_count(), atoms, "{case}");
            assert_eq!(topology.bond_count(), bonds, "{case}");
            assert_eq!(sizes.len(), molecules, "{case}");
            assert_eq!(&sizes[..largest.len()], largest.as_slice(), "{case}");
            assert_eq!(charges, vec![1; charged], "{case}");
            assert!(interpreted.report().instances().iter().all(|instance| {
                instance.bonds().len()
                    == topology
                        .instance(instance.instance())
                        .unwrap()
                        .molecule()
                        .bond_count()
            }));

            // The bundled templates do not cover AP5, so without the file's
            // table its bonds are inferred: as many as the dictionary has.
            let inferred = interpreted
                .report()
                .issues()
                .iter()
                .filter_map(|issue| match issue {
                    StructureInterpretIssue::ResidueBondsInferred {
                        residue,
                        component_id,
                        bonds,
                    } => Some((residue.clone(), component_id.clone(), *bonds)),
                    _ => None,
                })
                .collect::<Vec<_>>();
            let ap5 = |asym: &str| {
                (
                    SourceResidueId::from_parts(
                        "1",
                        asym,
                        SourceResiduePosition::Author("215".to_owned()),
                        None,
                    ),
                    "AP5".to_owned(),
                    62,
                )
            };
            let expected = if stripped && name == "1AKE" {
                vec![ap5("C"), ap5("D")]
            } else {
                Vec::new()
            };
            assert_eq!(inferred, expected, "{case}");
        }
    }
}

/// The atom name and component of an interpreted atom.
fn atom_label(
    report: &kekule::macromolecular::StructureInterpretationReport,
    atom: kekule::topology::InstanceAtomId,
) -> (String, String) {
    let provenance = report
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .find(|candidate| candidate.atom() == atom)
        .expect("every interpreted atom has provenance");
    (
        provenance.component_id().to_owned(),
        provenance.atom_name().to_owned(),
    )
}

/// The coordination records of a report as their partners' components and
/// atom names.
fn coordination_labels(
    report: &kekule::macromolecular::StructureInterpretationReport,
) -> Vec<[(String, String); 2]> {
    report
        .coordination()
        .iter()
        .map(|record| record.atoms().map(|atom| atom_label(report, atom)))
        .collect()
}

#[test]
fn metal_coordination_between_residues_is_recorded_not_bonded() {
    use kekule::descriptors::{molecular_formula, HydrogenCountPolicy};
    use kekule::topology::MoleculeClass;

    // KRas G12C with GDP, Mg2+, and the non-covalent inhibitor ASP2453
    // (A1L9E): the file lists Ser17 OG-Mg and GDP O2B-Mg as `metalc`.
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/mmcif/9UOH.cif");
    let text = std::fs::read_to_string(path).unwrap();
    let interpreted = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let report = interpreted.report();
    let expected = [
        [
            ("SER".to_owned(), "OG".to_owned()),
            ("MG".to_owned(), "MG".to_owned()),
        ],
        [
            ("GDP".to_owned(), "O2B".to_owned()),
            ("MG".to_owned(), "MG".to_owned()),
        ],
    ];
    assert_eq!(coordination_labels(report), expected);
    assert_eq!(
        report
            .coordination()
            .iter()
            .map(|record| (record.connection_id(), record.connection_type()))
            .collect::<Vec<_>>(),
        [(Some("metalc1"), "metalc"), (Some("metalc2"), "metalc")]
    );

    // The protein, GDP, the ion, and the inhibitor are separate molecules.
    let mut structure = interpreted.structure().clone();
    structure.perceive().unwrap();
    let molecule_of = |component: &str| {
        let instances = report
            .instances()
            .iter()
            .filter(|instance| {
                instance
                    .atoms()
                    .iter()
                    .any(|atom| atom.component_id() == component)
            })
            .map(|instance| instance.instance())
            .collect::<Vec<_>>();
        let [instance] = instances[..] else {
            panic!("{component} is in {} molecules", instances.len());
        };
        structure.topology().instance(instance).unwrap()
    };
    let formula = |component: &str| {
        molecular_formula(
            molecule_of(component).molecule(),
            HydrogenCountPolicy::IncludePerceived,
        )
        .unwrap()
        .to_string()
    };
    let protein = molecule_of("SER");
    assert_eq!(protein.class(), MoleculeClass::Protein);
    assert_ne!(protein.id(), molecule_of("GDP").id());
    assert_ne!(protein.id(), molecule_of("MG").id());
    assert_eq!(molecule_of("GDP").class(), MoleculeClass::SmallMolecule);
    assert_eq!(formula("GDP"), "C10H15N5O11P2");
    assert_eq!(molecule_of("MG").class(), MoleculeClass::Ion);
    assert_eq!(formula("MG"), "Mg+2");
    assert_eq!(molecule_of("A1L9E").class(), MoleculeClass::SmallMolecule);
    assert_eq!(formula("A1L9E"), "C40H48F3N7O4");
    assert!(protein
        .molecule()
        .atoms()
        .all(|(_, atom)| atom.element().symbol() != "Mg"));

    // Written with its report, the block states the coordination as
    // `metalc` again, and reading it back records the same partners.
    let (written, write_report) = mmcif::write(
        [mmcif::MmcifBlockSource::from(&interpreted)],
        Default::default(),
    )
    .unwrap();
    assert!(write_report.omissions().is_empty());
    assert_eq!(written.matches(" metalc ").count(), 2);
    let reread = mmcif::parse_str(&written, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(coordination_labels(reread.report()), expected);
    assert_eq!(
        reread.topology().instance_count(),
        interpreted.topology().instance_count()
    );
    assert_eq!(
        reread.topology().bond_count(),
        interpreted.topology().bond_count()
    );

    // The structure alone carries no coordination.
    let (structure_only, _) = mmcif::write(
        [mmcif::MmcifBlockSource::structure(interpreted.structure())],
        Default::default(),
    )
    .unwrap();
    assert!(!structure_only.contains("metalc"));
}

#[test]
fn identical_molecules_share_one_definition_and_edit_alone() {
    use kekule::molecule::Atom;
    use kekule::topology::MoleculeDefinitionId;
    use std::collections::BTreeMap;

    // 9UOH: one protein, GDP, Mg2+, the inhibitor, a sulfate, and waters.
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/mmcif/9UOH.cif");
    let text = std::fs::read_to_string(path).unwrap();
    let interpreted = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let report = interpreted.report();
    let topology = interpreted.topology();
    let component = |instance: kekule::topology::MoleculeInstanceId| {
        report
            .instances()
            .iter()
            .find(|candidate| candidate.instance() == instance)
            .unwrap()
            .atoms()[0]
            .component_id()
            .to_owned()
    };
    let mut counts = BTreeMap::<MoleculeDefinitionId, (String, usize)>::new();
    for instance in topology.instances() {
        counts
            .entry(instance.definition_id())
            .or_insert_with(|| (component(instance.id()), 0))
            .1 += 1;
    }
    let waters = topology
        .instances()
        .filter(|instance| component(instance.id()) == "HOH")
        .count();
    assert!(waters > 50);
    let mut listed = counts.into_values().collect::<Vec<_>>();
    listed.sort();
    assert_eq!(
        listed,
        [
            ("A1L9E", 1),
            ("GDP", 1),
            ("HOH", waters),
            ("MET", 1),
            ("MG", 1),
            ("SO4", 1),
        ]
        .map(|(component, count)| (component.to_owned(), count))
    );
    assert_eq!(topology.definition_count(), 6);
    assert_eq!(topology.instance_count(), waters + 5);

    // Instances of one definition map atom for atom: each water's local
    // atom 0 is its oxygen, named in its own atom site.
    let water_instances = topology
        .instances()
        .filter(|instance| component(instance.id()) == "HOH")
        .map(|instance| instance.id())
        .collect::<Vec<_>>();
    for instance in report
        .instances()
        .iter()
        .filter(|instance| water_instances.contains(&instance.instance()))
    {
        let [oxygen] = instance.atoms() else {
            panic!("a crystallographic water lists its oxygen alone");
        };
        assert_eq!((oxygen.atom().atom().index(), oxygen.atom_name()), (0, "O"));
    }

    // Default perception runs once per definition and serves every water.
    let mut perceived = interpreted.structure().clone();
    perceived.perceive().unwrap();
    let hydrogens = |structure: &kekule::structure::Structure, instance| {
        structure
            .topology()
            .instance(instance)
            .unwrap()
            .molecule()
            .total_hydrogens(kekule::molecule::AtomId::new(0))
            .unwrap()
    };
    assert!(water_instances
        .iter()
        .all(|&water| hydrogens(&perceived, water) == 2));

    // Editing one water copies its definition; the others keep theirs.
    let [edited, untouched, ..] = water_instances[..] else {
        panic!("several waters");
    };
    let mut editor = perceived.to_editor();
    let [oxygen, other_oxygen] = [edited, untouched].map(|instance| {
        editor
            .topology()
            .atom_handle(kekule::topology::InstanceAtomId::new(
                instance,
                kekule::molecule::AtomId::new(0),
            ))
            .unwrap()
    });
    let element = perceived
        .topology()
        .instance(edited)
        .unwrap()
        .molecule()
        .atom(kekule::molecule::AtomId::new(0))
        .unwrap()
        .element();
    editor
        .replace_atom(oxygen, Atom::new(element).with_formal_charge(-1))
        .unwrap();
    let (edited_structure, correspondence) = editor.finish_with_correspondence().unwrap();
    let topology = edited_structure.topology();
    assert_eq!(topology.definition_count(), 7);
    let charge = |instance| {
        topology
            .instance(instance)
            .unwrap()
            .molecule()
            .atom(kekule::molecule::AtomId::new(0))
            .unwrap()
            .formal_charge()
    };
    let published = |handle| correspondence.atom(handle).unwrap().0.instance();
    assert_eq!(charge(published(oxygen)), -1);
    assert_eq!(charge(published(other_oxygen)), 0);
    let untouched_definition = topology
        .instance(published(other_oxygen))
        .unwrap()
        .definition_id();
    assert_eq!(
        topology.instances_of(untouched_definition).count(),
        waters - 1
    );
}
