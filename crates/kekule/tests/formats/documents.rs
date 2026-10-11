use kekule::core::{Element, Molecule, Perception};
use kekule::diagnostics::{DiagnosticKind, SourceLocation, UnsupportedContent};
use kekule::geometry::Point3;
use kekule::{molfile, sdf, smiles};

const DISCONNECTED_MOLFILE: &str = "salt-like\nkekule\n\n  2  0  0  0  0  0            999 V2000\n    1.2500    2.5000    3.7500 Na  0  3  0  0  0  0  0  0  0  0  0  0\n   -4.0000    5.5000   -6.2500 Cl  0  5  0  0  0  0  0  0  0  0  0  0\nM  END\n";

const INTERLEAVED_COMPONENT_MOLFILE: &str = "interleaved\nkekule\n\n  4  1  0  0  0  0            999 V2000\n   10.0000   11.0000   12.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   20.0000   21.0000   22.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0\n   30.0000   31.0000   32.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   40.0000   41.0000   42.0000 Cl  0  5  0  0  0  0  0  0  0  0  0  0\n  1  3  1  0  0  0  0\nM  END\n";

#[test]
fn molfile_projections_retain_the_interpreted_model_and_reports() {
    let interpretation = molfile::parse_str(INTERLEAVED_COMPONENT_MOLFILE)
        .unwrap()
        .interpret()
        .unwrap();
    let topology = interpretation.model().shared_topology();
    let positions = interpretation.model().positions().values().value().as_ptr();
    let reports = interpretation.reports().to_vec();
    assert!(std::ptr::eq(interpretation.topology(), topology.as_ref()));
    assert!(std::sync::Arc::ptr_eq(
        &interpretation.clone().into_topology(),
        &topology
    ));
    assert!(std::sync::Arc::ptr_eq(
        &interpretation.clone().into_model().shared_topology(),
        &topology
    ));
    let (model, moved_reports, omitted) = interpretation.into_parts();
    assert!(omitted.is_empty());
    assert!(std::sync::Arc::ptr_eq(&model.shared_topology(), &topology));
    assert_eq!(model.positions().values().value().as_ptr(), positions);
    assert_eq!(moved_reports, reports);
    assert_eq!(
        moved_reports
            .iter()
            .map(|report| report
                .atom_mappings()
                .iter()
                .map(|mapping| mapping.location().line())
                .collect::<Vec<_>>())
            .collect::<Vec<_>>(),
        [vec![5, 7], vec![6], vec![8]]
    );
}

#[test]
fn empty_molfile_fails_at_interpretation_before_any_projection() {
    let source = "empty\nkekule\n\n  0  0  0  0  0  0            999 V2000\nM  END\n";
    let document = molfile::parse_str(source).unwrap();
    assert!(document.interpret().is_err());
    assert!(document.to_molecules().is_err());
    assert!(document.to_model().is_err());
    let sdf_document = sdf::parse_str(&format!("{source}$$$$\n")).unwrap();
    let error = sdf_document.records()[0].interpret().unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Chemistry);
    assert_eq!(error.record(), 1);
}

#[test]
fn sdf_model_can_install_perception_without_losing_source_order_or_geometry() {
    let text = format!("{INTERLEAVED_COMPONENT_MOLFILE}$$$$\n");
    let document = sdf::parse_str(&text).unwrap();
    let mut model = document.records()[0].to_model().unwrap();
    let original = model.clone();
    let positions_ptr = model.positions().values().value().as_ptr();
    assert!(model_molecules(&model)
        .iter()
        .all(|m| m.perception() == &Perception::default()));

    model.perceive().unwrap();

    assert!(model.topology().same_layout(original.topology()));
    assert_eq!(model.topology().atom_ids(), original.topology().atom_ids());
    assert_eq!(
        model.topology().hierarchy(),
        original.topology().hierarchy()
    );
    assert_eq!(model.positions(), original.positions());
    assert_eq!(model.positions().values().value().as_ptr(), positions_ptr);
    for (perceived, source) in model_molecules(&model)
        .into_iter()
        .zip(model_molecules(&original))
    {
        let mut expected = source.clone();
        expected.perceive().unwrap();
        assert_eq!(perceived, source);
        assert_eq!(perceived.perception(), expected.perception());
        assert_eq!(source.perception(), &Perception::default());
    }
}

fn element(molecule: &Molecule) -> Element {
    molecule
        .atoms()
        .next()
        .expect("one-atom component")
        .1
        .element
}

fn model_molecules(model: &kekule::structure::Model) -> Vec<&Molecule> {
    model
        .topology()
        .molecules()
        .map(|instance| instance.molecule())
        .collect()
}

fn assert_point_close(actual: Point3, expected: Point3) {
    assert!((actual.x - expected.x).abs() < 1.0e-15);
    assert!((actual.y - expected.y).abs() < 1.0e-15);
    assert!((actual.z - expected.z).abs() < 1.0e-15);
}

fn assert_points_close(actual: &[Point3], expected: &[Point3]) {
    assert_eq!(actual.len(), expected.len());
    for (&actual, &expected) in actual.iter().zip(expected) {
        assert_point_close(actual, expected);
    }
}

#[test]
fn parsed_smiles_document_converts_components_in_source_order() {
    let document = smiles::parse_str("[Na+].[Cl-]").expect("SMILES parses");
    let molecules = document.to_molecules().expect("SMILES interprets");

    assert_eq!(molecules.len(), 2);
    assert_eq!(element(&molecules[0]).symbol(), "Na");
    assert_eq!(element(&molecules[1]).symbol(), "Cl");
    assert!(molecules.iter().all(|molecule| molecule.atom_count() == 1));
    assert!(molecules
        .iter()
        .all(|molecule| molecule.perception() == &Perception::default()));
}

#[test]
fn molfile_document_model_retains_published_component_geometry() {
    let document = molfile::parse_str(DISCONNECTED_MOLFILE).expect("Molfile parses");
    assert_eq!(
        document,
        molfile::parse_str_with_options(
            DISCONNECTED_MOLFILE,
            molfile::MolfileParseOptions::default(),
        )
        .expect("explicit default Molfile parse")
    );
    let interpretation = document.interpret().expect("Molfile report interpretation");
    let expected = [
        Point3::new(0.125, 0.25, 0.375),
        Point3::new(-0.4, 0.55, -0.625),
    ];
    for ((instance, report), expected) in interpretation
        .topology()
        .molecules()
        .zip(interpretation.reports())
        .zip(expected)
    {
        let mapping = report
            .atom_mappings()
            .first()
            .expect("one source atom mapping");
        let position_index = interpretation
            .topology()
            .atom_index(kekule::topology::InstanceAtomId::new(
                instance.id(),
                mapping.atom(),
            ))
            .expect("mapped canonical atom");
        assert_point_close(
            interpretation
                .model()
                .positions()
                .position_at(position_index.index())
                .expect("mapped source coordinate")
                .into_value(),
            expected,
        );
    }
    let molecules = document.to_molecules().expect("Molfile interprets");
    let topology = document.to_topology().expect("Molfile topology builds");
    let model = document.to_model().expect("Molfile model builds");

    assert_eq!(molecules.len(), 2);
    assert_eq!(model.topology().instance_count(), 2);
    assert_eq!(model.atom_count(), 2);
    assert!(topology.same_layout(model.topology()));
    assert_eq!(
        model_molecules(&model),
        molecules.iter().collect::<Vec<_>>()
    );
    assert_points_close(
        model.positions().values().value(),
        &[
            Point3::new(0.125, 0.25, 0.375),
            Point3::new(-0.4, 0.55, -0.625),
        ],
    );
    assert_eq!(element(model_molecules(&model)[0]).symbol(), "Na");
    assert_eq!(element(model_molecules(&model)[1]).symbol(), "Cl");
    assert!(model_molecules(&model)
        .iter()
        .all(|molecule| molecule.perception() == &Perception::default()));
    let hierarchy = model.topology().hierarchy();
    let (_, chain) = hierarchy.chains().next().expect("synthetic chain");
    assert_eq!(chain.label_id(), "A");
    assert_eq!(chain.author_id(), None);
    assert_eq!(chain.residues().len(), 2);
    for (index, residue_id) in chain.residues().iter().enumerate() {
        let residue = hierarchy.residue(*residue_id).unwrap();
        assert_eq!(residue.name(), "UNL");
        assert_eq!(residue.label_seq_id(), Some((index + 1) as i32));
        assert_eq!(
            residue.author_seq_id(),
            Some((index + 1).to_string().as_str())
        );
        assert_eq!(residue.atom_sites().len(), 1);
    }
}

#[test]
fn molfile_model_keeps_interleaved_source_atom_order() {
    let document = molfile::parse_str(INTERLEAVED_COMPONENT_MOLFILE).expect("Molfile parses");
    let molecules = document.to_molecules().expect("Molfile interprets");
    let model = document.to_model().expect("Molfile model builds");

    assert_eq!(molecules.len(), 3);
    assert_eq!(model.topology().instance_count(), 3);
    assert_eq!(
        model_molecules(&model),
        molecules.iter().collect::<Vec<_>>()
    );
    assert_eq!(
        molecules[0]
            .atoms()
            .map(|(_, atom)| atom.element.symbol())
            .collect::<Vec<_>>(),
        ["C", "O"]
    );
    assert_eq!(element(&molecules[1]).symbol(), "Na");
    assert_eq!(element(&molecules[2]).symbol(), "Cl");
    // Dense atom order is the atom-block order even though components interleave,
    // so a coordinate file in this order addresses the same atoms.
    assert_eq!(
        model
            .topology()
            .atoms()
            .map(|atom| atom.element.symbol())
            .collect::<Vec<_>>(),
        ["C", "Na", "O", "Cl"]
    );
    let instances = model
        .topology()
        .atom_ids()
        .iter()
        .map(|atom| atom.molecule().index())
        .collect::<Vec<_>>();
    assert_eq!(instances, [0, 1, 0, 2]);
    assert_points_close(
        model.positions().values().value(),
        &[
            Point3::new(1.0, 1.1, 1.2),
            Point3::new(2.0, 2.1, 2.2),
            Point3::new(3.0, 3.1, 3.2),
            Point3::new(4.0, 4.1, 4.2),
        ],
    );
}

#[test]
fn sdf_parsed_records_remain_independent_conversion_boundaries() {
    let input = format!(
        "{DISCONNECTED_MOLFILE}>  <ROLE>\nions\n\n$$$$\nsecond\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    9.0000    8.0000    7.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n"
    );
    let document = sdf::parse_str(&input).expect("SDF parses");
    assert_eq!(
        document,
        sdf::parse_str_with_options(&input, sdf::SdfParseOptions::default())
            .expect("explicit default SDF parse")
    );
    let records: &[sdf::SdfRecord] = document.records();

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].number(), 1);
    assert_eq!(records[1].number(), 2);
    assert_eq!(records[0].location(), SourceLocation::at_line(1));
    assert_eq!(records[1].location(), SourceLocation::at_line(12));
    assert_eq!(records[0].data_fields()[0].value(), "ions");

    let first_interpretation = records[0].interpret().expect("rich record interpretation");
    assert_eq!(first_interpretation.title(), "salt-like");
    assert_eq!(first_interpretation.data_fields()[0].value(), "ions");
    assert_eq!(first_interpretation.report().record(), 1);
    assert_eq!(first_interpretation.report().molfile_components().len(), 2);
    assert_eq!(first_interpretation.molecules().count(), 2);
    assert!(first_interpretation
        .topology()
        .same_layout(first_interpretation.model().topology()));
    assert_points_close(
        first_interpretation.model().positions().values().value(),
        &[
            Point3::new(0.125, 0.25, 0.375),
            Point3::new(-0.4, 0.55, -0.625),
        ],
    );
    let projected_topology = first_interpretation.clone().into_topology();
    let projected_model = first_interpretation.clone().into_model();
    let projected_molecules = first_interpretation.into_molecules();
    assert!(projected_topology.same_layout(projected_model.topology()));
    assert_eq!(
        projected_topology
            .molecules()
            .map(|occurrence| occurrence.molecule())
            .collect::<Vec<_>>(),
        projected_molecules.iter().collect::<Vec<_>>()
    );
    assert_eq!(projected_topology.hierarchy().chains().count(), 1);
    assert_eq!(projected_topology.hierarchy().residues().count(), 2);

    let first_molecules = records[0].to_molecules().expect("first record interprets");
    let first_model = records[0].to_model().expect("first record model builds");
    let second_molecules = records[1].to_molecules().expect("second record interprets");
    let second_model = records[1].to_model().expect("second record model builds");

    assert_eq!(first_molecules.len(), 2);
    assert_eq!(first_model.topology().instance_count(), 2);
    assert_eq!(
        model_molecules(&first_model),
        first_molecules.iter().collect::<Vec<_>>()
    );
    assert_points_close(
        first_model.positions().values().value(),
        &[
            Point3::new(0.125, 0.25, 0.375),
            Point3::new(-0.4, 0.55, -0.625),
        ],
    );

    assert_eq!(second_molecules.len(), 1);
    assert_eq!(second_model.topology().instance_count(), 1);
    assert_points_close(
        second_model.positions().values().value(),
        &[Point3::new(0.9, 0.8, 0.7)],
    );
    assert!(first_molecules
        .iter()
        .chain(second_molecules.iter())
        .all(|molecule| molecule.perception() == &Perception::default()));
    assert!(model_molecules(&first_model)
        .into_iter()
        .chain(model_molecules(&second_model))
        .all(|molecule| molecule.perception() == &Perception::default()));

    let interpretations = document
        .records()
        .iter()
        .map(sdf::SdfRecord::interpret)
        .collect::<Result<Vec<_>, _>>()
        .expect("every record interprets on its own");
    assert_eq!(interpretations[0].report().record(), 1);
    assert_eq!(interpretations[1].report().record(), 2);
    assert_eq!(
        interpretations[1].report().location(),
        Some(SourceLocation::at_line(12))
    );

    let written = sdf::write(
        &[interpretations[1].clone()],
        sdf::SdfWriteOptions {
            version: sdf::MolfileWriteVersion::V2000,
        },
    )
    .expect("single-component rich record writes");
    let round_trip_document = sdf::parse_str(&written).expect("written record parses");
    let round_trip = round_trip_document.records()[0]
        .interpret()
        .expect("written record interprets");
    assert_eq!(round_trip.title(), "second");
    assert_points_close(
        round_trip.model().positions().values().value(),
        &[Point3::new(0.9, 0.8, 0.7)],
    );
}

#[test]
fn smiles_convenience_preserves_cardinality_and_matches_explicit_pipeline() {
    let ethanol = smiles::to_molecules("CCO").expect("ethanol interprets");
    assert_eq!(ethanol.len(), 1);

    let salt = smiles::to_molecules("[Na+].[Cl-]").expect("salt interprets");
    assert_eq!(salt.len(), 2);

    let document = smiles::parse_str("[Na+].[Cl-]").expect("salt parses");
    let explicit = smiles::interpret(&document)
        .expect("salt interprets explicitly")
        .into_molecules();
    assert_eq!(salt, explicit);
}

#[test]
fn component_interpretation_types_are_nameable_through_format_facades() {
    let smiles_document = smiles::parse_str("CC.O").expect("SMILES parses");
    let smiles_interpretation = smiles_document.interpret().expect("SMILES interprets");
    let smiles_components: &[smiles::SmilesComponentInterpretation] =
        smiles_interpretation.components();
    assert_eq!(smiles_components.len(), 2);
    let count_error: smiles::SmilesComponentCountError = smiles_interpretation
        .molecule()
        .expect_err("singular access rejects multiple components");
    assert_eq!(count_error.actual(), 2);

    let molfile_source = "two components\nkekule\n\n  2  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n    1.0000    0.0000    0.0000 O   0  0  0  0  0  0\nM  END\n";
    let molfile_document = kekule::molfile::parse_str(molfile_source).expect("Molfile parses");
    let molfile_interpretation = molfile_document.interpret().expect("Molfile interprets");
    let reports: &[kekule::molfile::MolfileInterpretationReport] = molfile_interpretation.reports();
    assert_eq!(reports.len(), 2);
    assert_eq!(molfile_interpretation.model().atom_count(), 2);
    assert_eq!(molfile_interpretation.topology().instance_count(), 2);
}

#[test]
fn sdf_interpretation_errors_carry_record_kind_and_input_location() {
    let valid = "valid\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n$$$$\n";
    let unknown = "unknown element\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 Xx  0  0  0  0  0  0\nM  END\n$$$$\n";
    let document = sdf::parse_str(&format!("{valid}{unknown}")).expect("SDF syntax parses");
    let error = document.records()[1]
        .interpret()
        .expect_err("unsupported element fails Molfile interpretation");

    assert_eq!(error.record(), 2);
    assert_eq!(error.kind(), DiagnosticKind::InvalidValue);
    // The atom line of the second record, counted in the SDF input.
    assert_eq!(error.location(), Some(SourceLocation::at_line(12)));
    assert_eq!(error.message(), "unsupported element symbol `Xx`");
    assert_eq!(
        error.to_string(),
        "SDF interpretation error in record 2 at line 12: unsupported element symbol `Xx`"
    );
    let source = std::error::Error::source(&error)
        .and_then(|source| source.downcast_ref::<molfile::MolfileInterpretError>())
        .expect("the Molfile failure is the source");
    assert_eq!(source.location(), error.location());
}

#[test]
fn smiles_default_options_and_topology_projection_are_consistent() {
    let default_document = smiles::parse_str("CCO.[Na+]").expect("default parse");
    let explicit_document =
        smiles::parse_str_with_options("CCO.[Na+]", smiles::SmilesParseOptions::default())
            .expect("explicit default parse");
    assert_eq!(default_document, explicit_document);

    let interpretation = default_document.interpret().expect("document interprets");
    assert_eq!(interpretation.molecules().count(), 2);
    assert!(interpretation
        .molecules()
        .all(|molecule| molecule.perception() == &kekule::core::Perception::default()));

    let topology = interpretation.into_topology().expect("topology builds");
    assert_eq!(topology.instance_count(), 2);
    assert_eq!(topology.atom_count(), 4);
    assert_eq!(
        topology
            .molecules()
            .map(|occurrence| occurrence.molecule().atom_count())
            .collect::<Vec<_>>(),
        vec![3, 1]
    );
    assert!(topology.hierarchy().is_empty());

    let concise = smiles::to_topology("CCO.[Na+]").expect("concise topology builds");
    assert!(concise.same_layout(&topology));
    assert_eq!(
        smiles::to_topology("CCO")
            .expect("connected topology")
            .instance_count(),
        1
    );
}

#[test]
fn many_smiles_components_preserve_source_mappings_and_stereo() {
    // Branches interleave disconnected atoms and ring closures reconnect
    // source fragments. Components cannot be obtained by merely splitting '.'.
    let block = "C1(O.F)CC1.N[C@@H](C)C(=O)O.F/C=C/Cl.c1ccncc1";
    let expected = smiles::parse_str(block).unwrap().interpret().unwrap();
    let source = vec![block; 256].join(".");
    let actual = smiles::parse_str(&source).unwrap().interpret().unwrap();
    assert_eq!(actual.components().len(), 256 * expected.components().len());
    for (index, component) in actual.components().iter().enumerate() {
        let baseline = &expected.components()[index % expected.components().len()];
        let offset = index / expected.components().len() * (block.len() + 1);
        assert_eq!(component.molecule(), baseline.molecule());
        assert_eq!(
            component.source_span(),
            baseline.source_span().start + offset..baseline.source_span().end + offset
        );
        let report = component.report();
        let expected_report = baseline.report();
        assert_eq!(
            report.atom_mappings().len(),
            expected_report.atom_mappings().len()
        );
        assert_eq!(
            report.bond_mappings().len(),
            expected_report.bond_mappings().len()
        );
        for (mapping, expected) in report
            .atom_mappings()
            .iter()
            .zip(expected_report.atom_mappings())
        {
            assert_eq!(mapping.atom(), expected.atom());
            assert_eq!(
                mapping.source_span(),
                expected.source_span().start + offset..expected.source_span().end + offset
            );
        }
        for (mapping, expected) in report
            .bond_mappings()
            .iter()
            .zip(expected_report.bond_mappings())
        {
            assert_eq!(mapping.bond(), expected.bond());
            assert_eq!(mapping.source_offset(), expected.source_offset() + offset);
        }
        assert_eq!(
            report.created_stereo_elements(),
            expected_report.created_stereo_elements()
        );
    }
}

fn v2000_with_properties(properties: &str) -> String {
    format!(
        "properties\nkekule\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n    1.2000    0.0000    0.0000 O   0  0  0  0  0  0\n  1  2  1  0  0  0  0\n{properties}M  END\n"
    )
}

fn v3000_with(atom_options: &str, bond_options: &str, blocks: &str) -> String {
    format!(
        "v3000\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0{atom_options}\nM  V30 2 O 1.2 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2{bond_options}\nM  V30 END BOND\n{blocks}M  V30 END CTAB\nM  END\n"
    )
}

#[test]
fn molfile_unsupported_content_is_rejected_by_default_and_reported_on_request() {
    use molfile::MolfileUnsupportedContent as Content;

    // Record syntax follows the BIOVIA CTfile Formats specification (2020).
    let cases = [
        (
            v2000_with_properties("M  STY  1   1 SUP\n"),
            Content::PropertyRecord("STY".to_owned()),
            8,
        ),
        (
            // An atom alias continues on the next line; both form one record.
            v2000_with_properties("A    1\nHO\n"),
            Content::Record("A".to_owned()),
            8,
        ),
        (
            v2000_with_properties("V    1 note\n"),
            Content::Record("V".to_owned()),
            8,
        ),
        (
            format!("{}trailing\n", v2000_with_properties("")),
            Content::TrailingContent,
            9,
        ),
        (
            v3000_with(" RGROUPS=(1 1)", "", ""),
            Content::AtomOption {
                atom: 0,
                key: "RGROUPS".to_owned(),
            },
            8,
        ),
        (
            v3000_with("", " TOPO=1", ""),
            Content::BondOption {
                bond: 0,
                key: "TOPO".to_owned(),
            },
            12,
        ),
        (
            v3000_with(
                "",
                "",
                "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/HILITE ATOMS=(1 1)\nM  V30 END COLLECTION\n",
            ),
            Content::Collection("MDLV30/HILITE".to_owned()),
            15,
        ),
        (
            v3000_with(
                "",
                "",
                "M  V30 BEGIN SGROUP\nM  V30 1 SUP 0 ATOMS=(1 2) LABEL=CO\nM  V30 END SGROUP\n",
            ),
            Content::Block("SGROUP".to_owned()),
            14,
        ),
        (
            v3000_with("", "", "M  V30 LINKNODE 1 2 2 1 2 1 2\n"),
            Content::Record("LINKNODE".to_owned()),
            14,
        ),
    ];
    let report = molfile::MolfileInterpretOptions {
        unsupported: UnsupportedContent::Report,
    };
    for (input, content, line) in cases {
        let expected = [molfile::MolfileUnsupported::new(
            content.clone(),
            SourceLocation::at_line(line),
        )];
        let document = molfile::parse_str(&input).expect("unsupported content parses");
        assert_eq!(document.unsupported(), expected, "{content}");

        let error = molfile::interpret(&document).expect_err("rejected by default");
        assert_eq!(error.kind(), DiagnosticKind::Unsupported, "{content}");
        assert_eq!(error.location(), Some(SourceLocation::at_line(line)));
        assert_eq!(error.message(), format!("unsupported {content}"));

        let interpretation =
            molfile::interpret_with_options(&document, report).expect("reported on request");
        assert_eq!(interpretation.omitted(), expected, "{content}");
        assert_eq!(interpretation.model().atom_count(), 2, "{content}");
        assert_eq!(interpretation.topology().instance_count(), 1, "{content}");
    }
}

#[test]
fn molfile_documents_expose_decoded_atoms_bonds_and_stereo_groups() {
    use molfile::{
        MolfileAtomParity as Parity, MolfileBondOrder as Order, MolfileBondStereo as Stereo,
        MolfileRadical as Radical,
    };
    type AtomFields = (
        String,
        [f64; 3],
        i32,
        Option<u16>,
        Option<Radical>,
        Option<u8>,
        Option<u8>,
        Option<u32>,
        Option<Parity>,
        usize,
    );
    fn atoms(document: &molfile::MolfileDocument) -> Vec<AtomFields> {
        document
            .atoms()
            .iter()
            .map(|atom| {
                let position = atom.position();
                (
                    atom.symbol().to_owned(),
                    [position.x, position.y, position.z],
                    atom.formal_charge(),
                    atom.isotope(),
                    atom.radical(),
                    atom.hydrogen_count(),
                    atom.valence(),
                    atom.atom_map(),
                    atom.parity(),
                    atom.location().line(),
                )
            })
            .collect()
    }
    fn bonds(
        document: &molfile::MolfileDocument,
    ) -> Vec<([usize; 2], Order, Option<Stereo>, usize)> {
        document
            .bonds()
            .iter()
            .map(|bond| {
                (
                    bond.atoms(),
                    bond.order(),
                    bond.stereo(),
                    bond.location().line(),
                )
            })
            .collect()
    }

    // Codes decoded per the BIOVIA CTfile Formats specification (2020):
    // charge code 3 is +1 and 4 a doublet radical, hydrogen code 2 is H1.
    let v2000 = molfile::parse_str(
        "typed\nkekule\n\n  3  2  0  0  1  0            999 V2000\n    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n    1.5000    0.0000    0.0000 C   0  0  0  2  0  4  0  0  0  7  0  0\n    2.0000    1.0000    0.0000 O   0  4  0  0  0  0  0  0  0  0  0  0\n  1  2  1  1  0  0  0\n  2  3  2  0  0  0  0\nM  ISO  1   1  15\nM  END\n",
    )
    .unwrap();
    assert_eq!(v2000.version(), molfile::MolfileVersion::V2000);
    assert_eq!(v2000.header().title(), "typed");
    assert!(v2000.chiral_flag());
    let expected: [AtomFields; 3] = [
        (
            "N".to_owned(),
            [0.0, 0.0, 0.0],
            1,
            Some(15),
            None,
            None,
            None,
            None,
            None,
            5,
        ),
        (
            "C".to_owned(),
            [1.5, 0.0, 0.0],
            0,
            None,
            None,
            Some(1),
            Some(4),
            Some(7),
            None,
            6,
        ),
        (
            "O".to_owned(),
            [2.0, 1.0, 0.0],
            0,
            None,
            Some(Radical::Doublet),
            None,
            None,
            None,
            None,
            7,
        ),
    ];
    assert_eq!(atoms(&v2000), expected);
    assert_eq!(
        bonds(&v2000),
        [
            ([0, 1], Order::Single, Some(Stereo::Up), 8),
            ([1, 2], Order::Double, None, 9),
        ]
    );
    assert!(v2000.stereo_groups().is_empty());
    assert!(v2000.unsupported().is_empty());

    // A continued logical line is located at its first physical line.
    let v3000 = molfile::parse_str(
        "typed\nkekule\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 1\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 CFG=2 MASS=13\nM  V30 2 N 1.5 0 0 3 CHG=1 HCOUNT=-1 -\nM  V30 VAL=4\nM  V30 3 O -1 0.5 0 0 RAD=2\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 CFG=3\nM  V30 2 1 1 3\nM  V30 END BOND\nM  V30 BEGIN COLLECTION\nM  V30 MDLV30/STERAC2 ATOMS=(1 1)\nM  V30 END COLLECTION\nM  V30 END CTAB\nM  END\n",
    )
    .unwrap();
    assert_eq!(v3000.version(), molfile::MolfileVersion::V3000);
    assert!(v3000.chiral_flag());
    let expected: [AtomFields; 3] = [
        (
            "C".to_owned(),
            [0.0, 0.0, 0.0],
            0,
            Some(13),
            None,
            None,
            None,
            None,
            Some(Parity::Even),
            8,
        ),
        (
            "N".to_owned(),
            [1.5, 0.0, 0.0],
            1,
            None,
            None,
            Some(0),
            Some(4),
            Some(3),
            None,
            9,
        ),
        (
            "O".to_owned(),
            [-1.0, 0.5, 0.0],
            0,
            None,
            Some(Radical::Doublet),
            None,
            None,
            None,
            None,
            11,
        ),
    ];
    assert_eq!(atoms(&v3000), expected);
    assert_eq!(
        bonds(&v3000),
        [
            ([0, 1], Order::Single, Some(Stereo::Down), 14),
            ([0, 2], Order::Single, None, 15),
        ]
    );
    let [group] = v3000.stereo_groups() else {
        panic!("one stereo group");
    };
    assert_eq!(group.kind(), molfile::MolfileStereoGroupKind::And(2));
    assert_eq!(group.atoms(), [0]);
    assert_eq!(group.location(), SourceLocation::at_line(18));
}

fn carbon_record(title: &str) -> String {
    format!(
        "{title}\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n$$$$\n"
    )
}

type StreamItem = Result<(usize, usize, String), (usize, DiagnosticKind, usize)>;

fn stream(records: sdf::SdfRecords<&[u8]>) -> Vec<StreamItem> {
    records
        .map(|record| match record {
            Ok(record) => Ok((
                record.number(),
                record.location().line(),
                record.title().to_owned(),
            )),
            Err(error) => Err((error.record(), error.kind(), error.location().line())),
        })
        .collect()
}

#[test]
fn sdf_record_stream_isolates_malformed_records_and_resumes() {
    let malformed =
        "second\nkekule\n\n  1  0  0  0  0  0            999 V2000\natom record\nM  END\n$$$$\n";
    let mut input = carbon_record("first").into_bytes();
    input.extend_from_slice(malformed.as_bytes());
    input.extend_from_slice(b"th\xffird\n");
    input.extend_from_slice(carbon_record("").as_bytes()[1..].as_ref());
    input.extend_from_slice(carbon_record("fourth").as_bytes());

    // Each record occupies seven lines.
    assert_eq!(
        stream(sdf::records(input.as_slice())),
        [
            Ok((1, 1, "first".to_owned())),
            Err((2, DiagnosticKind::Syntax, 12)),
            Err((3, DiagnosticKind::Io, 15)),
            Ok((4, 22, "fourth".to_owned())),
        ]
    );

    // An oversized record is skipped to its delimiter.
    let long_title = "x".repeat(200);
    let input = format!(
        "{}{}{}",
        carbon_record("first"),
        carbon_record(&long_title),
        carbon_record("third")
    );
    let options = sdf::SdfParseOptions {
        max_record_bytes: 150,
        ..sdf::SdfParseOptions::default()
    };
    assert_eq!(
        stream(sdf::records_with_options(input.as_bytes(), options)),
        [
            Ok((1, 1, "first".to_owned())),
            Err((2, DiagnosticKind::Limit, 8)),
            Ok((3, 15, "third".to_owned())),
        ]
    );

    // Exceeding the record count ends the stream.
    let options = sdf::SdfParseOptions {
        max_records: 1,
        ..sdf::SdfParseOptions::default()
    };
    assert_eq!(
        stream(sdf::records_with_options(input.as_bytes(), options)),
        [
            Ok((1, 1, "first".to_owned())),
            Err((2, DiagnosticKind::Limit, 8)),
        ]
    );

    // The strict document parser fails on the first malformed record.
    let error = sdf::parse_str(&format!("{}{malformed}", carbon_record("first"))).unwrap_err();
    assert_eq!(
        (error.record(), error.kind(), error.location()),
        (2, DiagnosticKind::Syntax, SourceLocation::at_line(12))
    );
}

#[test]
fn sdf_and_molfile_locations_count_lines_of_the_original_input() {
    let unknown = "unknown\nkekule\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n    1.0000    0.0000    0.0000 Xx  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n>  <NAME>\nvalue\n\n$$$$\n";
    let lf = format!("{}{unknown}", carbon_record("first"));
    for input in [lf.clone(), lf.replace('\n', "\r\n"), lf.replace('\n', "\r")] {
        let document = sdf::parse_str(&input).unwrap();
        let record = &document.records()[1];
        assert_eq!(record.location(), SourceLocation::at_line(8));
        assert_eq!(
            record.molfile().atoms()[1].location(),
            SourceLocation::at_line(13)
        );
        assert_eq!(
            record.molfile().bonds()[0].location(),
            SourceLocation::at_line(14)
        );
        assert_eq!(
            record.data_fields()[0].location(),
            Some(SourceLocation::at_line(16))
        );
        let error = record.interpret().unwrap_err();
        assert_eq!(error.location(), Some(SourceLocation::at_line(13)));

        let streamed = sdf::records(input.as_bytes())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(streamed, document.records());
    }

    let molfile = "lone CR\rkekule\r\r  1  0  0  0  0  0            999 V2000\r    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\rM  END\r";
    let interpretation = molfile::parse_str(molfile).unwrap().interpret().unwrap();
    assert_eq!(
        interpretation.reports()[0].atom_mappings()[0].location(),
        SourceLocation::at_line(5)
    );
}

#[test]
fn molfile_and_sdf_readers_match_string_parsing_and_enforce_limits() {
    let molfile = "reader\nkekule\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n";
    assert_eq!(
        molfile::parse_reader(molfile.as_bytes()).unwrap(),
        molfile::parse_str(molfile).unwrap()
    );
    let limited = molfile::MolfileParseOptions {
        max_input_bytes: molfile.len() - 1,
        ..molfile::MolfileParseOptions::default()
    };
    let error = molfile::parse_reader_with_options(molfile.as_bytes(), limited).unwrap_err();
    assert_eq!(
        (error.kind(), error.location()),
        (DiagnosticKind::Limit, SourceLocation::at_line(1))
    );
    let error = molfile::parse_reader(b"\xff\xfe".as_slice()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Io);

    let sdf_text = format!("{molfile}$$$$\n");
    assert_eq!(
        sdf::parse_reader(sdf_text.as_bytes()).unwrap(),
        sdf::parse_str(&sdf_text).unwrap()
    );
    let limited = sdf::SdfParseOptions {
        max_input_bytes: sdf_text.len() - 1,
        ..sdf::SdfParseOptions::default()
    };
    let error = sdf::parse_reader_with_options(sdf_text.as_bytes(), limited).unwrap_err();
    assert_eq!(
        (error.record(), error.kind(), error.location()),
        (1, DiagnosticKind::Limit, SourceLocation::at_line(1))
    );
}
