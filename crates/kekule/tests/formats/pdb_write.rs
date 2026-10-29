use std::collections::BTreeSet;
use std::sync::Arc;

use kekule::diagnostics::{
    DiagnosticKind, FormatErrorDetail, StructureField, UnrepresentablePolicy, WriteOmission,
};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::macromolecular::{ModelSelection, StructureInterpretOptions};
use kekule::molecule::{AtomId, BondOrder};
use kekule::pdb::{
    self, PdbAtomField, PdbCellReason, PdbChainIdPolicy, PdbWriteIssue, PdbWriteOptions,
};
use kekule::solvation::{self, BoxShape, PeriodicBoxOptions, SolventOptions};
use kekule::structure::{
    Ensemble, EnsemblePartition, EnsembleSample, EnsembleState, Model, ModelBuilder, Positions,
    Series, SeriesAxis, SeriesItem,
};
use kekule::topology::{InstanceAtomId, Topology};
use kekule::units::{Quantity, ANGSTROM, DEGREE, MOLAR, NANOMETER};
use kekule::{mmcif, smiles};

const CRAMBIN: &str = include_str!("../fixtures/pdb/1CRN.pdb");
const ADENYLATE_KINASE: &str = include_str!("../fixtures/pdb/1AKE.pdb");

fn interpret(text: &str) -> Model {
    pdb::parse_str(text, Default::default())
        .unwrap()
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .unwrap()
        .into_model()
}

fn bonds(topology: &Topology) -> BTreeSet<(usize, usize, String)> {
    topology
        .bonds()
        .map(|bond| {
            let [left, right] = bond.atoms().map(|atom| atom.index().index());
            (
                left.min(right),
                left.max(right),
                format!("{:?}", bond.bond().order),
            )
        })
        .collect()
}

fn issue(error: &kekule::diagnostics::FormatError) -> &PdbWriteIssue {
    match error.detail() {
        Some(FormatErrorDetail::PdbWrite(issue)) => issue,
        other => panic!("not a PDB write issue: {other:?}"),
    }
}

#[test]
fn deposited_entries_round_trip_through_the_writer() {
    for (name, text) in [("1CRN", CRAMBIN), ("1AKE", ADENYLATE_KINASE)] {
        let model = interpret(text);
        let written = pdb::write(&model, Default::default())
            .map(|(text, _)| text)
            .unwrap();
        let reread = interpret(&written);
        let (left, right) = (reread.topology(), model.topology());
        assert_eq!(left.atom_count(), right.atom_count(), "{name}");
        assert_eq!(bonds(left), bonds(right), "{name}");
        assert_eq!(left.instance_count(), right.instance_count(), "{name}");
        for (a, b) in left.residues().zip(right.residues()) {
            assert_eq!(
                (a.name(), a.author_seq_id(), a.insertion_code()),
                (b.name(), b.author_seq_id(), b.insertion_code()),
                "{name}"
            );
            assert_eq!(a.chain().author_id(), b.chain().author_id(), "{name}");
        }
        assert_eq!(reread.cell(), model.cell(), "{name}");
        // PDB -> model -> PDB is a fixed point.
        assert_eq!(
            pdb::write(&reread, Default::default())
                .map(|(text, _)| text)
                .unwrap(),
            written,
            "{name}"
        );
        // Disulfides and ligand bonds are CONECT records; polymer bonds are not.
        assert!(
            written.lines().any(|line| line.starts_with("CONECT")),
            "{name}"
        );
        assert!(
            written.lines().any(|line| line.starts_with("TER")),
            "{name}"
        );
    }
    // 1CRN's P 1 21 1 cell is not the model's cell, so no CRYST1 is written.
    let crambin = pdb::write(&interpret(CRAMBIN), Default::default())
        .map(|(text, _)| text)
        .unwrap();
    assert!(interpret(CRAMBIN).cell().is_none());
    assert!(!crambin.contains("CRYST1"));
    // Three disulfides, each listed from both atoms.
    assert_eq!(
        crambin
            .lines()
            .filter(|line| line.starts_with("CONECT"))
            .count(),
        6
    );
}

/// An entry's text with its chains named by segment, as CHARMM and NAMD
/// name them: chain `A` blank in the chain column and `PROA` in columns
/// 73-76.
fn named_by_segment(text: &str) -> String {
    let mut segmented = String::new();
    for line in text.lines() {
        let mut bytes = line.as_bytes().to_vec();
        if line.starts_with("ATOM") || line.starts_with("HETATM") {
            let chain = bytes[21];
            bytes[21] = b' ';
            bytes[72..76].copy_from_slice(&[b'P', b'R', b'O', chain]);
        } else if line.starts_with("TER") {
            bytes[21] = b' ';
        }
        segmented.push_str(std::str::from_utf8(&bytes).unwrap());
        segmented.push('\n');
    }
    segmented
}

#[test]
fn chains_named_by_segment_round_trip_under_the_segment_policy() {
    let model = interpret(ADENYLATE_KINASE);
    let (named, _) = pdb::write(&model, Default::default()).unwrap();
    let segmented = named_by_segment(&named);
    let reread = interpret(&segmented);
    let chains = reread
        .topology()
        .chains()
        .map(|chain| chain.author_id().unwrap_or(chain.label_id()).to_owned())
        .collect::<BTreeSet<_>>();
    assert_eq!(chains, ["PROA", "PROB"].map(str::to_owned).into());

    // The chains fail by default, and are written back as segments.
    let error = pdb::write(&reread, Default::default()).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::ChainIdTooLong { chain_id, .. } if chain_id == "PROA"
    ));
    assert!(error.message().contains("PdbChainIdPolicy::SegmentId"));
    let options = PdbWriteOptions::default().with_chain_ids(PdbChainIdPolicy::SegmentId);
    let (written, report) = pdb::write(&reread, options).unwrap();
    assert_eq!(written, segmented);
    assert_eq!(report.omissions(), [WriteOmission::ChainIds { chains: 2 }]);
    // Each chain's polymer still ends in its own TER record.
    let ters = |text: &str| text.lines().filter(|line| line.starts_with("TER")).count();
    assert_eq!(ters(&written), 2);
    assert_eq!(ters(&written), ters(&named));
    // One-character chains are written as before under the policy.
    assert_eq!(
        pdb::write(&model, options).unwrap(),
        (named, Default::default())
    );
}

/// Methane solvated in a 1.5 nm padded cube of water at 0.15 M salt.
fn solvated() -> Model {
    let mut model = Model::from_molecule(
        smiles::to_molecules("C").unwrap().remove(0),
        &Positions::zeros(1),
    )
    .unwrap();
    solvation::add_periodic_box(
        &mut model,
        &PeriodicBoxOptions::Padding {
            padding: Quantity::new(1.5, NANOMETER),
            shape: BoxShape::Cube,
        },
    )
    .unwrap();
    solvation::solvate(
        &mut model,
        &SolventOptions::default().with_ionic_strength(Quantity::new(0.15, MOLAR)),
    )
    .unwrap();
    model
}

#[test]
fn solvated_systems_export_their_solvent_chains_through_every_writer() {
    let model = solvated();
    let solvent = model
        .topology()
        .chains()
        .map(|chain| chain.label_id().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(solvent, ["SOL1", "ION1", "ION2"]);
    // Solvent counts no implicit hydrogen, as PDB and mmCIF readers infer.
    for instance in model.topology().instances().skip(1) {
        let molecule = instance.definition().molecule();
        for (atom, _) in molecule.atoms() {
            assert_eq!(molecule.implicit_hydrogens(atom), Some(0));
        }
    }
    // By default PDB refuses the chain identifiers, and nothing else.
    let error = pdb::write(&model, Default::default()).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::ChainIdTooLong { chain_id, .. } if chain_id == "SOL1"
    ));

    // PDB writes the solvent chains as segments, which read back as chains.
    let options = PdbWriteOptions::default().with_chain_ids(PdbChainIdPolicy::SegmentId);
    let (text, report) = pdb::write(&model, options).unwrap();
    assert!(report.omissions().contains(&WriteOmission::ChainIds {
        chains: solvent.len()
    }));
    for line in text.lines().filter(|line| line.starts_with("HETATM")) {
        if &line[17..20] == "HOH" {
            assert_eq!((&line[21..22], &line[72..76]), (" ", "SOL1"), "{line}");
        }
    }
    let reread = interpret(&text);
    assert_eq!(
        reread.topology().atom_count(),
        model.topology().atom_count()
    );
    let chains = reread
        .topology()
        .chains()
        .map(|chain| chain.author_id().unwrap_or(chain.label_id()).to_owned())
        .collect::<BTreeSet<_>>();
    assert!(
        solvent.iter().all(|chain| chains.contains(chain)),
        "{chains:?}"
    );
    assert_eq!(pdb::write(&reread, options).unwrap().0, text);

    // mmCIF names the chains as they are.
    let (cif, _) =
        mmcif::write([mmcif::MmcifBlockSource::model(&model)], Default::default()).unwrap();
    let reread = mmcif::parse_str(&cif, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_model();
    assert_eq!(
        reread.topology().atom_count(),
        model.topology().atom_count()
    );
    let chains = reread
        .topology()
        .chains()
        .map(|chain| chain.author_id().unwrap_or(chain.label_id()).to_owned())
        .collect::<BTreeSet<_>>();
    assert!(
        solvent.iter().all(|chain| chains.contains(chain)),
        "{chains:?}"
    );
}

#[test]
fn chain_ids_longer_than_a_segment_fail_under_the_segment_policy() {
    let mut editor = solvated().to_editor();
    let (solvent, _) = editor.topology().chains().next().unwrap();
    editor
        .set_chain_identifiers(solvent, "SOLVENT", None)
        .unwrap();
    let renamed = editor.finish().unwrap();
    let options = PdbWriteOptions::default().with_chain_ids(PdbChainIdPolicy::SegmentId);
    let error = pdb::write(&renamed, options).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::ChainIdTooLong { chain_id, .. } if chain_id == "SOLVENT"
    ));
}

#[test]
fn bonds_without_a_template_are_conect_listings_by_order() {
    // Bracket atoms fix their implicit hydrogen counts, which PDB cannot
    // carry; here the counts are those a reader infers.
    let topology = smiles::to_topology("CC(=O)[O-].[Na+]").unwrap();
    let points = (0..5)
        .map(|index| Point3::new(f64::from(index) * 1.4, 0.0, 0.0))
        .collect::<Vec<_>>();
    let model = Model::new(
        topology,
        Positions::from_vec(Quantity::new(points, ANGSTROM)).unwrap(),
    )
    .unwrap();
    let omit = PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let (written, report) = pdb::write(&model, omit).unwrap();
    // [O-] and [Na+] fix their counts.
    assert_eq!(
        report.omissions(),
        [WriteOmission::ImplicitHydrogens { atoms: 2 }]
    );
    let lines = written.lines().collect::<Vec<_>>();
    assert_eq!(
        lines[0],
        "HETATM    1  C1  MOL     1       0.000   0.000   0.000                       C"
    );
    assert_eq!(
        lines[3],
        "HETATM    4  O4  MOL     1       4.200   0.000   0.000                       O1-"
    );
    assert_eq!(
        lines[4],
        "HETATM    5 NA1  MOL     2       5.600   0.000   0.000                      NA1+"
    );
    // The double bond is listed twice from each atom.
    assert!(lines.contains(&"CONECT    2    1    3    3    4"));
    assert!(lines.contains(&"CONECT    3    2    2"));
    let reread = interpret(&written);
    assert_eq!(bonds(reread.topology()), bonds(model.topology()));
}

#[test]
fn large_residue_numbers_use_hybrid_36() {
    let text = "\
HETATM    1  O   HOH A9999       0.000   0.000   0.000  1.00  0.00           O
HETATM    2  O   HOH AA000       3.000   0.000   0.000  1.00  0.00           O
END
";
    let model = interpret(text);
    let written = pdb::write(&model, Default::default())
        .map(|(text, _)| text)
        .unwrap();
    assert_eq!(written, text);
    let numbers = model
        .topology()
        .residues()
        .map(|residue| residue.author_seq_id().unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(numbers, ["9999", "10000"]);
}

#[test]
fn identifiers_pdb_cannot_hold_fail_explicitly() {
    let cif = |chain: &str, component: &str| {
        format!(
            "data_x
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
HETATM 1 O O {component} {chain} 1 . 0.000 0.000 0.000 1 {chain} 1
"
        )
    };
    let model = |text: &str| {
        mmcif::parse_str(text, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap()
            .into_model()
    };
    let error = pdb::write(&model(&cif("AB", "HOH")), Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::ChainIdTooLong { chain_id, .. } if chain_id == "AB"
    ));
    let error = pdb::write(&model(&cif("A", "LONG")), Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::ResidueNameTooLong { name, .. } if name == "LONG"
    ));
    // Omitting chemistry does not relax identifiers.
    let options = PdbWriteOptions::default().with_chemistry(UnrepresentablePolicy::Omit);
    assert!(pdb::write(&model(&cif("AB", "HOH")), options)
        .map(|(text, _)| text)
        .is_err());
}

#[test]
fn chemistry_pdb_cannot_carry_fails_unless_omitted() {
    // A bracket atom fixes its implicit hydrogen count.
    let topology = smiles::to_topology("[CH4]").unwrap();
    let model = Model::new(topology, Positions::zeros(1)).unwrap();
    let error = pdb::write(&model, Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::UnsupportedAtomField {
            field: PdbAtomField::Hydrogens,
            ..
        }
    ));
    // The chemistry policy no longer covers implicit hydrogens.
    let chemistry = PdbWriteOptions::default().with_chemistry(UnrepresentablePolicy::Omit);
    assert!(pdb::write(&model, chemistry).is_err());
    let hydrogens = PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let (_, report) = pdb::write(&model, hydrogens).unwrap();
    assert_eq!(
        report.omissions(),
        [WriteOmission::ImplicitHydrogens { atoms: 1 }]
    );

    let topology = smiles::to_topology("C$C").unwrap();
    let model = Model::new(topology, Positions::zeros(2)).unwrap();
    let error = pdb::write(&model, Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::UnsupportedBondOrder {
            order: BondOrder::Quadruple,
            ..
        }
    ));
    // Omitting implicit hydrogens does not omit other chemistry.
    assert!(pdb::write(&model, hydrogens).is_err());
    let (written, report) = pdb::write(&model, chemistry).unwrap();
    assert!(!written.contains("CONECT"));
    assert_eq!(report.omissions(), [WriteOmission::Bonds { bonds: 1 }]);
}

#[test]
fn residues_must_be_consecutive_in_dense_order() {
    let water = smiles::to_molecules("[H]O[H]").unwrap().remove(0);
    let mut builder = ModelBuilder::new();
    let positions = Positions::zeros(3);
    let first = builder.add_molecule(water.clone(), &positions).unwrap();
    let second = builder.add_molecule(water, &positions).unwrap();
    let atom = |instance, index| InstanceAtomId::new(instance, AtomId::new(index));
    builder
        .set_atom_order([
            atom(first, 0),
            atom(second, 0),
            atom(first, 1),
            atom(first, 2),
            atom(second, 1),
            atom(second, 2),
        ])
        .unwrap();
    let model = builder.build().unwrap();
    let omit = PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let error = pdb::write(&model, omit).map(|(text, _)| text).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::NonContiguousResidue { .. }
    ));
}

#[test]
fn ensembles_write_models_and_sample_state_fails_unless_omitted() {
    let topology = smiles::to_topology("C").unwrap();
    let sample = |x: f64| {
        EnsembleSample::new(
            Positions::from_vec(Quantity::new(vec![Point3::new(x, 0.0, 0.0)], ANGSTROM)).unwrap(),
            1.0,
        )
        .unwrap()
    };
    let mut ensemble = Ensemble::new(Arc::clone(&topology));
    ensemble
        .push_state(sample(0.0), EnsembleState::new())
        .unwrap();
    ensemble.push_sample(sample(1.0), 0).unwrap();
    let written = pdb::write(&ensemble, Default::default())
        .map(|(text, _)| text)
        .unwrap();
    assert_eq!(written.matches("MODEL ").count(), 2);
    assert_eq!(written.matches("ENDMDL").count(), 2);
    let reread = pdb::parse_str(&written, Default::default())
        .unwrap()
        .interpret_ensemble(Default::default())
        .unwrap();
    assert_eq!(reread.ensemble().len(), 2);
    assert_eq!(
        pdb::write(reread.ensemble(), Default::default())
            .map(|(text, _)| text)
            .unwrap(),
        written
    );

    // A second state fails, and so does its energy, which comes first.
    ensemble
        .regroup(EnsemblePartition::state_per_sample(2))
        .unwrap();
    let error = pdb::write(&ensemble, Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::UnrepresentableEnsembleSample {
            sample: 2,
            field: StructureField::State
        }
    ));
    ensemble
        .state_mut(1)
        .unwrap()
        .set_energy(Some(Quantity::new(1.0, kekule::units::KILOJOULE_PER_MOLE)))
        .unwrap();
    let error = pdb::write(&ensemble, Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::UnrepresentableEnsembleSample {
            sample: 2,
            field: StructureField::Energy
        }
    ));
    let options = PdbWriteOptions::default().with_ensemble_metadata(UnrepresentablePolicy::Omit);
    let (_, report) = pdb::write(&ensemble, options).unwrap();
    assert_eq!(
        report.omissions(),
        [
            WriteOmission::StructureState {
                field: StructureField::Energy,
                structures: 1
            },
            WriteOmission::StructureState {
                field: StructureField::State,
                structures: 2
            }
        ]
    );
    // A lone sample's weight is relative to samples the entry does not hold.
    let lone = ensemble.get(0).unwrap();
    assert!(
        pdb::write(pdb::PdbSource::structure(&lone), Default::default())
            .map(|(text, _)| text)
            .is_err()
    );
    let (_, report) = pdb::write(pdb::PdbSource::structure(&lone), options).unwrap();
    assert_eq!(
        report.omissions(),
        [WriteOmission::StructureState {
            field: StructureField::Weight,
            structures: 1
        }]
    );
}

#[test]
fn series_write_models_and_their_places_fail_unless_omitted() {
    let topology = smiles::to_topology("C").unwrap();
    let positions = |x: f64| {
        Positions::from_vec(Quantity::new(vec![Point3::new(x, 0.0, 0.0)], ANGSTROM)).unwrap()
    };
    let mut scan = Series::new(Arc::clone(&topology), SeriesAxis::Coordinate(DEGREE));
    for (index, angle) in [0.0, 120.0, 240.0].into_iter().enumerate() {
        let place = Quantity::new(angle, DEGREE);
        scan.push(SeriesItem::at_coordinate(positions(index as f64), place).unwrap())
            .unwrap();
    }
    let error = pdb::write(&scan, Default::default()).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::UnrepresentableSeriesItem {
            item: 1,
            field: StructureField::ItemCoordinate
        }
    ));
    let options = PdbWriteOptions::default().with_series_metadata(UnrepresentablePolicy::Omit);
    let (written, report) = pdb::write(&scan, options).unwrap();
    assert_eq!(
        report.omissions(),
        [WriteOmission::StructureState {
            field: StructureField::ItemCoordinate,
            structures: 3
        }]
    );
    assert_eq!(written.matches("MODEL ").count(), 3);
    assert_eq!(written.matches("ENDMDL").count(), 3);
    // The models read back as an ensemble of the items' positions.
    let reread = pdb::parse_str(&written, Default::default())
        .unwrap()
        .interpret_ensemble(Default::default())
        .unwrap();
    let ensemble = reread.ensemble();
    assert_eq!(ensemble.len(), 3);
    for (sample, item) in ensemble.iter().zip(&scan) {
        assert_eq!(sample.positions(), item.positions());
    }

    // A step history fails on its steps, and a lone item on its own.
    let mut history = Series::new(topology, SeriesAxis::Step);
    history
        .push(SeriesItem::at_step(positions(0.0), 3))
        .unwrap();
    let error = pdb::write(&history, Default::default()).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteIssue::UnrepresentableSeriesItem {
            item: 1,
            field: StructureField::ItemStep
        }
    ));
    let lone = history.get(0).unwrap();
    let error = pdb::write(pdb::PdbSource::structure(&lone), Default::default()).unwrap_err();
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::UnrepresentableState(
            StructureField::ItemStep
        ))
    );
    let (_, report) = pdb::write(pdb::PdbSource::structure(&lone), options).unwrap();
    assert_eq!(
        report.omissions(),
        [WriteOmission::StructureState {
            field: StructureField::ItemStep,
            structures: 1
        }]
    );

    // An empty series is no entry.
    let empty = Series::new(smiles::to_topology("C").unwrap(), SeriesAxis::Step);
    let error = pdb::write(&empty, options).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert!(matches!(issue(&error), PdbWriteIssue::EmptySeries));
}

#[test]
fn cells_crystal_records_cannot_carry_fail_unless_omitted() {
    let topology = smiles::to_topology("C").unwrap();
    let cell = |lengths: Vector3, periodic: [bool; 3]| {
        PeriodicCell::orthorhombic(Quantity::new(lengths, ANGSTROM), periodic).unwrap()
    };
    let model = |cell| {
        Model::new(Arc::clone(&topology), Positions::zeros(1))
            .unwrap()
            .with_cell(cell)
    };
    let written = pdb::write(
        &model(cell(Vector3::new(10.0, 20.0, 30.0), [true; 3])),
        Default::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    assert!(written
        .starts_with("CRYST1   10.000   20.000   30.000  90.00  90.00  90.00 P 1           1\n"));
    for (lengths, periodic, reason) in [
        (
            Vector3::new(10.0, 20.0, 30.0),
            [true, true, false],
            PdbCellReason::PartiallyPeriodic,
        ),
        (
            Vector3::new(1.0, 1.0, 1.0),
            [true; 3],
            PdbCellReason::Placeholder,
        ),
    ] {
        let source = model(cell(lengths, periodic));
        let error = pdb::write(&source, Default::default())
            .map(|(text, _)| text)
            .unwrap_err();
        assert_eq!(
            issue(&error),
            &PdbWriteIssue::UnrepresentableCell { reason },
            "{reason:?}"
        );
        let options = PdbWriteOptions::default().with_cell(UnrepresentablePolicy::Omit);
        let (text, report) = pdb::write(&source, options).unwrap();
        assert!(!text.contains("CRYST1"));
        assert_eq!(report.omissions(), [WriteOmission::Cell { structures: 1 }]);
    }
}

#[test]
fn sink_failures_are_io_errors_after_validation() {
    struct Failing;
    impl std::io::Write for Failing {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("full"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let model = interpret(CRAMBIN);
    let error = pdb::write_to(&mut Failing, &model, Default::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Io);
    let mut sink = Vec::new();
    pdb::write_to(&mut sink, &model, Default::default()).unwrap();
    assert_eq!(
        String::from_utf8(sink).unwrap(),
        pdb::write(&model, Default::default())
            .map(|(text, _)| text)
            .unwrap()
    );
}
