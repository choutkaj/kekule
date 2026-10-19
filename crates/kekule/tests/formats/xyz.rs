use kekule::diagnostics::{DiagnosticKind, Format, FormatErrorDetail, UnsupportedContent};
use kekule::molecule::{BondOrder, Element, ImplicitHydrogens};
use kekule::topology::TopologyAtomIndex;
use kekule::units::{ANGSTROM, BOHR, KELVIN};
use kekule::xyz::{
    self, XyzEnsembleInterpretOptions, XyzFrameSelection, XyzInterpretIssue, XyzInterpretOptions,
};

const WATER: &str = "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n";

fn interpret(text: &str, options: XyzInterpretOptions) -> xyz::XyzInterpretation {
    xyz::parse_str(text, Default::default())
        .unwrap()
        .interpret(options)
        .unwrap()
}

fn interpret_error(text: &str, options: XyzInterpretOptions) -> kekule::diagnostics::FormatError {
    xyz::parse_str(text, Default::default())
        .unwrap()
        .interpret(options)
        .unwrap_err()
}

fn index(row: u32) -> TopologyAtomIndex {
    TopologyAtomIndex::new(row)
}

#[test]
fn benzene_infers_single_bonds_and_fixes_hydrogens_at_zero() {
    // A regular hexagon of carbons (C-C 1.39 A) with hydrogens 1.08 A out.
    let mut text = String::from("12\nbenzene\n");
    for (radius, element) in [(1.39, "C"), (2.47, "H")] {
        for k in 0..6 {
            let angle = f64::from(k) * std::f64::consts::PI / 3.0;
            text.push_str(&format!(
                "{element} {:.4} {:.4} 0.0\n",
                radius * angle.cos(),
                radius * angle.sin()
            ));
        }
    }
    let interpretation = interpret(&text, Default::default());
    let topology = interpretation.topology();
    assert_eq!(topology.instance_count(), 1);
    assert_eq!(topology.atom_count(), 12);
    assert_eq!(topology.bond_count(), 12);
    let report = interpretation.report();
    assert_eq!(report.inferred_bonds().len(), 12);
    assert!(report
        .inferred_bonds()
        .iter()
        .all(|bond| bond.order() == BondOrder::Single));
    assert!(report.issues().is_empty());
    // Dense atom order is the atom-line order.
    for row in 0..12 {
        let element = topology.atom_at(index(row)).unwrap().element();
        assert_eq!(element, if row < 6 { Element::C } else { Element::H });
    }
    let molecule = topology.definitions().next().unwrap().molecule();
    for atom in molecule.atom_ids() {
        assert_eq!(
            molecule.atom(atom).unwrap().hydrogens(),
            ImplicitHydrogens::Fixed(0)
        );
    }
    let bond = report.inferred_bonds()[0];
    assert_eq!(bond.atoms(), [index(0), index(1)]);
    assert!((bond.distance().value_in(ANGSTROM).unwrap() - 1.39).abs() < 1e-3);
    assert!((bond.cutoff().value_in(ANGSTROM).unwrap() - (0.76 + 0.76 + 0.45)).abs() < 1e-12);
    assert_eq!(report.bond_tolerance().value_in(ANGSTROM).unwrap(), 0.45);
}

#[test]
fn interleaved_molecules_keep_atom_line_order() {
    // Two waters 3 A apart, listed oxygens first.
    let text = "6\n\nO 0 0 0\nO 3 0 0\nH 0.96 0 0\nH 3.96 0 0\nH -0.24 0.93 0\nH 2.76 0.93 0\n";
    let interpretation = interpret(text, Default::default());
    let topology = interpretation.topology();
    assert_eq!(topology.instance_count(), 2);
    assert_eq!(topology.bond_count(), 4);
    let instance_of = |row| topology.atom_id(index(row)).unwrap().instance();
    assert_eq!(instance_of(0), instance_of(2));
    assert_eq!(instance_of(0), instance_of(4));
    assert_eq!(instance_of(1), instance_of(3));
    assert_ne!(instance_of(0), instance_of(1));
    let positions = interpretation.model().positions();
    let x = |row: usize| positions.get(row).unwrap().value_in(ANGSTROM).unwrap().x;
    assert!((x(1) - 3.0).abs() < 1e-12 && (x(3) - 3.96).abs() < 1e-12);
    assert_eq!(topology.residues().count(), 2);
    assert!(topology.residues().all(|residue| residue.name() == "UNL"));
}

#[test]
fn metal_ligand_contacts_are_dative_and_close_contacts_unbonded() {
    let text = "4\n\nFe 0 0 0\nCl 2.2 0 0\nHe 5 5 5\nHe 5.3 5 5\n";
    let interpretation = interpret(text, Default::default());
    let report = interpretation.report();
    assert_eq!(report.inferred_bonds().len(), 1);
    assert_eq!(report.inferred_bonds()[0].order(), BondOrder::Dative);
    assert_eq!(interpretation.topology().instance_count(), 3);
    let [XyzInterpretIssue::CloseContact { atoms, distance }] = report.issues() else {
        panic!("one close contact: {:?}", report.issues())
    };
    assert_eq!(*atoms, [index(2), index(3)]);
    assert!((distance.value_in(ANGSTROM).unwrap() - 0.3).abs() < 1e-9);
}

#[test]
fn elements_without_a_covalent_radius_use_the_reported_fallback() {
    let interpretation = interpret("2\n\nEs 0 0 0\nEs 9 0 0\n", Default::default());
    let report = interpretation.report();
    assert_eq!(
        report
            .fallback_covalent_radius()
            .value_in(ANGSTROM)
            .unwrap(),
        0.77
    );
    let fallback = report
        .issues()
        .iter()
        .filter(|issue| matches!(issue, XyzInterpretIssue::FallbackCovalentRadius { .. }))
        .count();
    assert_eq!(fallback, 2);
}

#[test]
fn element_tokens_accept_any_case_and_atomic_numbers() {
    let interpretation = interpret(
        "3\n\nc 0 0 0\nCL 1.76 0 0\n6 -1.5 0 0\n",
        Default::default(),
    );
    let topology = interpretation.topology();
    let elements = (0..3)
        .map(|row| topology.atom_at(index(row)).unwrap().element())
        .collect::<Vec<_>>();
    assert_eq!(
        elements,
        [Element::C, Element::from_symbol("Cl").unwrap(), Element::C]
    );
    let error = interpret_error("1\n\nXx 0 0 0\n", Default::default());
    assert_eq!(
        (error.format(), error.kind()),
        (Format::Xyz, DiagnosticKind::InvalidValue)
    );
    assert_eq!(error.location().unwrap().line(), 3);
}

#[test]
fn extra_columns_fail_unless_reported_as_omitted() {
    let text = "1\nProperties=species:S:1:pos:R:3:charge:R:1\nHe 0 0 0 0.5\n";
    let error = interpret_error(text, Default::default());
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);
    assert_eq!(error.location().unwrap().line(), 3);
    let options = XyzInterpretOptions::default().with_extra_columns(UnsupportedContent::Report);
    let interpretation = interpret(text, options);
    let [XyzInterpretIssue::ExtraColumnsOmitted { atom, location }] =
        interpretation.report().issues()
    else {
        panic!("one omission")
    };
    assert_eq!((*atom, location.line()), (index(0), 3));
    assert_eq!(
        interpretation.report().comment(),
        "Properties=species:S:1:pos:R:3:charge:R:1"
    );
}

#[test]
fn frame_selection_follows_the_request() {
    let two = format!("{WATER}{}", WATER.replace("water", "second"));
    let error = interpret_error(&two, Default::default());
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    let first = interpret(
        &two,
        XyzInterpretOptions::default().with_frame(XyzFrameSelection::First),
    );
    assert_eq!(first.report().frame(), 0);
    let second = interpret(
        &two,
        XyzInterpretOptions::default().with_frame(XyzFrameSelection::Index(1)),
    );
    assert_eq!(
        (second.report().frame(), second.report().comment()),
        (1, "second")
    );
    let error = interpret_error(
        &two,
        XyzInterpretOptions::default().with_frame(XyzFrameSelection::Index(2)),
    );
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::UnknownRequestedModel {
            model_id: "2".into()
        })
    );
    let error = interpret_error("", Default::default());
    assert_eq!(error.detail(), Some(&FormatErrorDetail::NoCoordinateModels));
}

#[test]
fn frames_of_the_same_connectivity_form_an_ensemble() {
    let stretched = "3\nstretched\nO 0 0 0\nH 1.0 0 0\nH -0.25 0.95 0\n";
    let document = xyz::parse_str(&format!("{WATER}{stretched}"), Default::default()).unwrap();
    let interpretation = document.interpret_ensemble(Default::default()).unwrap();
    let ensemble = interpretation.ensemble();
    assert_eq!(ensemble.len(), 2);
    assert_eq!(interpretation.reports().len(), 2);
    assert_eq!(interpretation.reports()[1].comment(), "stretched");
    let x = ensemble.get(1).unwrap().positions().get(1).unwrap();
    assert!((x.value_in(ANGSTROM).unwrap().x - 1.0).abs() < 1e-12);
    assert_eq!(
        ensemble.get(0).unwrap().weight(),
        ensemble.get(1).unwrap().weight()
    );

    let selected = document
        .interpret_ensemble(
            XyzEnsembleInterpretOptions::default()
                .with_frames(Some(vec![1]))
                .unwrap(),
        )
        .unwrap();
    assert_eq!(selected.reports()[0].frame(), 1);
    assert!(XyzEnsembleInterpretOptions::default()
        .with_frames(Some(Vec::new()))
        .is_err());
    let error = document
        .interpret_ensemble(
            XyzEnsembleInterpretOptions::default()
                .with_frames(Some(vec![0, 0]))
                .unwrap(),
        )
        .unwrap_err();
    assert!(matches!(
        error.detail(),
        Some(FormatErrorDetail::DuplicateRequestedModel { .. })
    ));
}

#[test]
fn frames_that_disagree_fail_the_ensemble() {
    let dissociated = "3\n\nO 0 0 0\nH 0.957 0 0\nH -2.4 9.3 0\n";
    let document = xyz::parse_str(&format!("{WATER}{dissociated}"), Default::default()).unwrap();
    let error = document.interpret_ensemble(Default::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::InconsistentTopology {
            model_id: "1".into()
        })
    );
    let other_atoms = "3\n\nS 0 0 0\nH 0.957 0 0\nH -0.24 0.927 0\n";
    let document = xyz::parse_str(&format!("{WATER}{other_atoms}"), Default::default()).unwrap();
    let error = document.interpret_ensemble(Default::default()).unwrap_err();
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::InconsistentAtomSet {
            model_id: "1".into()
        })
    );
}

#[test]
fn length_units_scale_bond_inference() {
    let bohr = "2\n\nH 0 0 0\nH 1.4 0 0\n";
    let interpretation = interpret(
        bohr,
        XyzInterpretOptions::default()
            .with_length_unit(BOHR)
            .unwrap(),
    );
    assert_eq!(interpretation.report().inferred_bonds().len(), 1);
    let distance = interpretation.report().inferred_bonds()[0].distance();
    assert!((distance.value_in(BOHR).unwrap() - 1.4).abs() < 1e-9);
    let position = interpretation.model().positions().get(1).unwrap();
    assert!((position.value_in(BOHR).unwrap().x - 1.4).abs() < 1e-9);
    // Read as angstrom, the same numbers are too far apart to bond.
    assert!(interpret(bohr, Default::default())
        .report()
        .inferred_bonds()
        .is_empty());
    let error = XyzInterpretOptions::default()
        .with_length_unit(KELVIN)
        .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::InvalidValue);
}

#[test]
fn overlapping_atoms_hit_the_contact_limit() {
    let mut text = String::from("40\n\n");
    for k in 0..40 {
        text.push_str(&format!("C {} 0 0\n", f64::from(k) * 0.01));
    }
    let error = interpret_error(&text, Default::default());
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    assert_eq!(
        error.limit_exceeded().unwrap().limit(),
        xyz::DEFAULT_XYZ_MAX_ATOM_CONTACTS
    );
    let options = XyzInterpretOptions::default()
        .with_max_atom_contacts(64)
        .unwrap();
    assert_eq!(interpret(&text, options).topology().bond_count(), 0);
    assert!(XyzInterpretOptions::default()
        .with_max_atom_contacts(0)
        .is_err());
}

#[test]
fn water_boxes_partition_into_molecules() {
    let mut text = String::new();
    let mut count = 0;
    for i in 0..10 {
        for j in 0..10 {
            for k in 0..10 {
                let (x, y, z) = (f64::from(i) * 3.1, f64::from(j) * 3.1, f64::from(k) * 3.1);
                text.push_str(&format!(
                    "O {x} {y} {z}\nH {} {y} {z}\nH {} {} {z}\n",
                    x + 0.957,
                    x - 0.24,
                    y + 0.927
                ));
                count += 3;
            }
        }
    }
    let text = format!("{count}\nbox\n{text}");
    let interpretation = interpret(&text, Default::default());
    let topology = interpretation.topology();
    assert_eq!(topology.instance_count(), 1000);
    assert_eq!(topology.bond_count(), 2000);
}

#[test]
fn later_frames_keep_element_errors_and_locate_mismatches() {
    // An invalid element in a later frame fails as it would in the first.
    let document = xyz::parse_str(
        &format!("{WATER}3\n\nXx 0 0 0\nH 1 0 0\nH 0 1 0\n"),
        Default::default(),
    )
    .unwrap();
    let error = document.interpret_ensemble(Default::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::InvalidValue);
    assert_eq!(error.location().unwrap().line(), 8);
    // Mismatches name the frame's count line.
    let document = xyz::parse_str(
        &format!("{WATER}3\n\nS 0 0 0\nH 0.957 0 0\nH -0.24 0.927 0\n"),
        Default::default(),
    )
    .unwrap();
    let error = document.interpret_ensemble(Default::default()).unwrap_err();
    assert_eq!(error.location().unwrap().line(), 6);
    let document = xyz::parse_str(
        &format!("{WATER}3\n\nO 0 0 0\nH 0.957 0 0\nH -2.4 9.3 0\n"),
        Default::default(),
    )
    .unwrap();
    let error = document.interpret_ensemble(Default::default()).unwrap_err();
    assert!(matches!(
        error.detail(),
        Some(FormatErrorDetail::InconsistentTopology { .. })
    ));
    assert_eq!(error.location().unwrap().line(), 6);
}

#[test]
fn coordinates_beyond_the_grid_range_fail() {
    let error = interpret_error("2\n\nC 0 0 0\nC 2e15 0 0\n", Default::default());
    assert_eq!(error.kind(), DiagnosticKind::InvalidValue);
    assert_eq!(error.location().unwrap().line(), 4);
    // Far apart but indexable atoms are simply unbonded.
    let interpretation = interpret("2\n\nC 0 0 0\nC 1e14 0 0\n", Default::default());
    assert_eq!(interpretation.topology().instance_count(), 2);
}

#[test]
fn xyz_structures_export_to_pdb_and_mmcif_when_implicit_hydrogens_are_omitted() {
    use kekule::diagnostics::{UnrepresentablePolicy, WriteOmission};
    use kekule::mmcif::{self, MmcifWriteOptions};
    use kekule::pdb::{self, PdbWriteOptions};

    // XYZ lists every atom, so every count is fixed at zero.
    let model = interpret(WATER, Default::default()).into_model();
    let omitted = [WriteOmission::ImplicitHydrogens { atoms: 3 }];

    assert_eq!(
        pdb::write(&model, Default::default()).unwrap_err().kind(),
        DiagnosticKind::Unsupported
    );
    let options = PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let (text, report) = pdb::write(&model, options).unwrap();
    assert_eq!(report.omissions(), omitted);
    let reread = pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(reread.topology().atom_count(), 3);
    assert_eq!(reread.topology().bond_count(), 2);

    assert_eq!(
        mmcif::write([&model], Default::default())
            .unwrap_err()
            .kind(),
        DiagnosticKind::Unsupported
    );
    let options = MmcifWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let (text, report) = mmcif::write([&model], options.clone()).unwrap();
    assert_eq!(report.omissions(), omitted);
    let reread = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(reread.topology().atom_count(), 3);
    // `write_to` returns the same report.
    let mut sink = Vec::new();
    let streamed = mmcif::write_to(&mut sink, [&model], options).unwrap();
    assert_eq!(streamed.omissions(), omitted);
    assert_eq!(String::from_utf8(sink).unwrap(), text);
}

#[test]
fn omitting_implicit_hydrogens_still_rejects_other_chemistry() {
    use kekule::diagnostics::UnrepresentablePolicy;
    use kekule::mmcif::{self, MmcifWriteOptions};
    use kekule::pdb::{self, PdbWriteOptions};
    use kekule::structure::{Model, Positions};

    // Carbon-13 is an isotope neither format can name.
    let topology = kekule::smiles::to_topology("[13CH4]").unwrap();
    let model = Model::new(topology, Positions::zeros(1)).unwrap();
    let pdb_options =
        PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let error = pdb::write(&model, pdb_options).unwrap_err();
    assert!(error.message().contains("isotope"), "{error}");
    let mmcif_options =
        MmcifWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let error = mmcif::write([&model], mmcif_options).unwrap_err();
    assert!(error.message().contains("isotope"), "{error}");
}
