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

/// The cell of an interpreted model as its diagonal in angstroms, when the
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
    let cell = interpreted.model().cell().map(|cell| {
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
fn mmcif_cell_is_the_model_cell_only_for_a_whole_unit_cell() {
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
    assert_eq!(error.message(), "invalid _cell.length_b `abc`");
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
        block_interpreted.model().positions(),
        interpreted.model().positions()
    );
    assert_eq!(interpreted.model().topology().instance_count(), 1);
    assert!(!interpreted.model().topology().hierarchy().is_empty());
    assert_eq!(interpreted.model().positions().len(), 2);
    assert_eq!(interpreted.molecules().count(), 1);
    assert!(interpreted
        .molecules()
        .all(|molecule| molecule.perception() == &kekule::perception::Perception::default()));
    assert!(interpreted
        .topology()
        .same_layout(interpreted.model().topology()));
    let owned_topology = method_default.clone().into_topology();
    let owned_model = method_default.clone().into_model();
    let owned_molecules = method_default.into_molecules();
    assert!(owned_topology.same_layout(owned_model.topology()));
    assert_eq!(
        owned_topology
            .molecules()
            .map(|occurrence| occurrence.molecule())
            .collect::<Vec<_>>(),
        owned_molecules.iter().collect::<Vec<_>>()
    );
    assert!(block
        .interpret(Default::default())?
        .into_topology()
        .same_layout(block.interpret(Default::default())?.into_model().topology()));
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
        document.interpret(Default::default())?.into_model(),
        document
            .interpret(StructureInterpretOptions::default())?
            .into_model(),
    ] {
        let from_block = block.interpret(Default::default())?.into_model();
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
        .same_layout(block.interpret(Default::default())?.into_model().topology()));
    assert_eq!(
        document.interpret_ensemble(Default::default())?.reports(),
        document
            .interpret_ensemble(macromolecular::StructureEnsembleInterpretOptions::default())?
            .reports()
    );
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
                .map(|interpretation| interpretation.into_model())
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
        provenance.atoms()[0].atom().molecule(),
        provenance.molecule()
    );
    let written = mmcif::write(
        [mmcif::MmcifBlockSource::model(interpreted.model())
            .with_reports(std::slice::from_ref(interpreted.report()))],
        MmcifWriteOptions::default(),
    )?;
    assert!(written.starts_with("data_model\n"));
    assert!(mmcif::parse_str(&written, Default::default()).is_ok());
    let custom_written = mmcif::write(
        [mmcif::MmcifBlockSource::model(interpreted.model())
            .with_reports(std::slice::from_ref(interpreted.report()))],
        MmcifWriteOptions::default()
            .with_block_name("custom".to_owned())
            .unwrap(),
    )?;
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
    let automatic = mmcif::write([interpreted.model()], MmcifWriteOptions::default())?;
    assert!(automatic.contains("1 non-polymer"));
    let mut classifications = MmcifEntityClassifications::new();
    classifications.insert(provenance.molecule(), EntityKind::Polymer)?;
    let written = mmcif::write(
        [mmcif::MmcifBlockSource::model(interpreted.model())
            .with_classifications(&classifications)],
        MmcifWriteOptions::default(),
    )?;
    assert!(written.contains("1 polymer"));
    let model = interpreted.into_model();
    assert_eq!(model.topology().instance_count(), 1);
    assert_eq!(model.positions().len(), 2);

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
        .model()
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
        "conflicting authoritative bond evidence for one atom pair: Double versus Single"
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
                .molecules()
                .map(|instance| instance.molecule().atom_count())
                .collect::<Vec<_>>();
            sizes.sort_unstable_by(|a, b| b.cmp(a));
            let charges = topology
                .molecules()
                .flat_map(|instance| instance.molecule().atoms())
                .filter(|(_, atom)| atom.formal_charge != 0)
                .map(|(_, atom)| atom.formal_charge)
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
                        .molecule(instance.molecule())
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
                    SourceResidueId {
                        model_id: "1".to_owned(),
                        asym_id: asym.to_owned(),
                        position: SourceResiduePosition::Author("215".to_owned()),
                        insertion_code: None,
                    },
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
