use kekule::diagnostics::{LimitExceeded, Resource};
use kekule::dssp::{
    self, DsspChainBreak, DsspError, DsspLimits, DsspOptions, DsspSecondaryStructure,
};
use kekule::macromolecular::{ModelSelection, StructureInterpretOptions};
use kekule::mmcif;
use kekule::topology::ResidueId;
use kekule::units::{DEGREE, KILOCALORIE_PER_MOLE};

const CRAMBIN_MMCIF: &str = include_str!("../fixtures/mmcif/1CRN.cif");

fn crambin_model() -> kekule::structure::Model {
    let document = mmcif::parse_str(CRAMBIN_MMCIF, Default::default())
        .expect("checked-in RCSB 1CRN fixture parses");
    document
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .expect("checked-in RCSB 1CRN fixture interprets")
        .into_model()
}

#[test]
fn dssp_matches_biopython_mkdssp_4_6_1_for_crambin() {
    let model = crambin_model();
    let result = dssp::assign(&model, DsspOptions::default())
        .expect("crambin has an analyzable protein backbone");
    let residues = result.residues().collect::<Vec<_>>();
    let codes = residues
        .iter()
        .map(|residue| residue.secondary_structure().code())
        .collect::<String>();

    assert_eq!(residues.len(), 46);
    let first_key: ResidueId = residues[0].key();
    assert!(result.get(first_key).is_some());
    assert_eq!(codes, " EE SSHHHHHHHHHHHTTT  HHHHHHHHS EE SSS   GGG  ");
    assert_eq!(residues[0].chain_break(), DsspChainBreak::NewChain);
    assert_eq!(
        residues[1].secondary_structure(),
        DsspSecondaryStructure::ExtendedStrand
    );
    assert!(
        (residues[1]
            .phi()
            .expect("defined phi")
            .value_in(DEGREE)
            .unwrap()
            - -107.8)
            .abs()
            < 0.15
    );
    assert!(
        (residues[1]
            .psi()
            .expect("defined psi")
            .value_in(DEGREE)
            .unwrap()
            - 144.3)
            .abs()
            < 0.15
    );

    let residue_three = residues[2];
    let strongest_acceptor = residue_three.acceptors()[0].expect("reference acceptor");
    let strongest_donor = residue_three.donors()[0].expect("reference donor");
    assert_eq!(
        result
            .get(strongest_acceptor.partner())
            .expect("acceptor partner is in the result")
            .source()
            .label_sequence_id(),
        Some(33)
    );
    assert!(
        (strongest_acceptor
            .energy()
            .value_in(KILOCALORIE_PER_MOLE)
            .unwrap()
            - -2.4)
            .abs()
            < 0.05
    );
    assert_eq!(
        result
            .get(strongest_donor.partner())
            .expect("donor partner is in the result")
            .source()
            .label_sequence_id(),
        Some(33)
    );
    assert!(
        (strongest_donor
            .energy()
            .value_in(KILOCALORIE_PER_MOLE)
            .unwrap()
            - -2.8)
            .abs()
            < 0.05
    );
}

#[test]
fn dssp_is_a_coordinate_snapshot_and_does_not_mutate_the_model() {
    let mut model = crambin_model();
    let before = model.clone();
    let assigned = dssp::assign(&model, DsspOptions::default()).expect("initial assignment");
    assert_eq!(model, before);

    let first_atom = model.topology().atom_ids()[0];
    let position = model
        .position(first_atom)
        .expect("first position")
        .map(|mut point| {
            point.x += 0.25;
            point
        });
    model
        .set_position(first_atom, position)
        .expect("finite coordinate update");

    assert_eq!(assigned.residues().count(), 46);
    assert_ne!(model, before);
}

#[test]
fn dssp_rejects_invalid_options_and_residue_limits() {
    let model = crambin_model();
    assert_eq!(
        DsspOptions::default().with_min_polyproline_stretch(4),
        Err(DsspError::InvalidPolyprolineStretch { value: 4 })
    );

    let options = DsspOptions::default().with_limits(DsspLimits::default().with_max_residues(2));
    assert_eq!(
        dssp::assign(&model, options),
        Err(DsspError::LimitExceeded(LimitExceeded::first_beyond(
            Resource::Residues,
            2
        )))
    );

    let options =
        DsspOptions::default().with_limits(DsspLimits::default().with_max_candidate_pairs(0));
    assert_eq!(
        dssp::assign(&model, options),
        Err(DsspError::LimitExceeded(LimitExceeded::first_beyond(
            Resource::CandidatePairs,
            0
        )))
    );

    let options = DsspOptions::default().with_limits(DsspLimits::default().with_max_ladders(0));
    assert_eq!(
        dssp::assign(&model, options),
        Err(DsspError::LimitExceeded(LimitExceeded::first_beyond(
            Resource::Ladders,
            0
        )))
    );
}

#[test]
fn dssp_rejects_coordinates_outside_its_spatial_index_range_without_panicking() {
    let mut model = crambin_model();
    let first_atom = model.topology().atom_ids()[0];
    let position = model
        .position(first_atom)
        .expect("first position")
        .map(|mut point| {
            point.x = f32::MAX as f64;
            point
        });
    model
        .set_position(first_atom, position)
        .expect("coordinate remains finite in the model");

    assert!(matches!(
        dssp::assign(&model, DsspOptions::default()),
        Err(DsspError::CoordinateOutOfRange { .. })
    ));
}

#[test]
fn dssp_codes_cover_the_complete_dssp4_alphabet() {
    let codes = DsspSecondaryStructure::ALL.map(DsspSecondaryStructure::code);
    assert_eq!(codes, [' ', 'H', 'B', 'E', 'G', 'I', 'P', 'T', 'S']);
    for code in codes {
        assert_eq!(DsspSecondaryStructure::try_from(code).unwrap().code(), code);
    }
    assert_eq!(
        DsspSecondaryStructure::try_from('-').unwrap(),
        DsspSecondaryStructure::Loop
    );
    assert_eq!(
        DsspSecondaryStructure::try_from('X').unwrap_err().code(),
        'X'
    );
}

#[test]
fn dssp_analyzes_canonical_amino_acids_without_label_sequence_numbers() {
    let reference = dssp::assign(&crambin_model(), DsspOptions::default()).unwrap();
    let source = CRAMBIN_MMCIF.replace("_atom_site.label_seq_id", "_atom_site.audit_sequence");
    let document = mmcif::parse_str(&source, Default::default()).unwrap();
    let model = document
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .unwrap()
        .into_model();
    assert!(model
        .topology()
        .residues()
        .all(|residue| residue.label_seq_id().is_none()));
    let assigned = dssp::assign(&model, DsspOptions::default()).unwrap();
    assert_eq!(assigned.statistics(), reference.statistics());
    for (actual, expected) in assigned.residues().zip(reference.residues()) {
        assert_eq!(
            actual.source().author_sequence_id(),
            expected.source().author_sequence_id()
        );
        assert_eq!(actual.secondary_structure(), expected.secondary_structure());
        assert_eq!(actual.phi(), expected.phi());
        assert_eq!(actual.psi(), expected.psi());
        assert_eq!(actual.kappa(), expected.kappa());
    }
}

#[test]
fn dssp_instance_membership_handles_many_reused_solvent_instances() {
    use kekule::geometry::Point3;
    use kekule::molecule::AtomId;
    use kekule::structure::Positions;
    use kekule::topology::{AtomSiteMetadata, InstanceAtomId};
    use kekule::units::{Quantity, ANGSTROM};

    let mut model = crambin_model().into_builder();
    let water = kekule::smiles::to_molecules("O").unwrap().pop().unwrap();
    let definition = model.add_molecule_definition(water.clone()).unwrap();
    let water_position =
        Positions::new(Quantity::new(vec![Point3::new(100., 100., 100.)], ANGSTROM)).unwrap();
    let chain = model.hierarchy_mut().add_chain("water", None).unwrap();
    let mut ignored = Vec::new();
    let count = 10_000;
    for index in 0..count {
        let instance = model.add_instance(definition, &water_position).unwrap();
        if index % 11 == 0 {
            ignored.push(instance);
        } else {
            let residue = model
                .hierarchy_mut()
                .add_residue(chain, "HOH", None, None, None)
                .unwrap();
            model
                .hierarchy_mut()
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, AtomId::new(0)),
                    AtomSiteMetadata::default(),
                )
                .unwrap();
        }
    }
    let model = model.build().unwrap();
    let result = dssp::assign(&model, DsspOptions::default()).unwrap();
    assert_eq!(result.report().ignored_instances(), ignored);
    assert_eq!(result.statistics().analyzed_residues(), 46);
    assert_eq!(
        result.report().non_peptide_residues(),
        count - ignored.len()
    );
}

#[test]
fn dssp_numbering_fallback_is_independent_of_unlabeled_solvent_in_the_chain() {
    use kekule::geometry::Point3;
    use kekule::molecule::AtomId;
    use kekule::structure::{Model, Positions};
    use kekule::topology::{AtomSiteMetadata, InstanceAtomId, TopologyBuilder};
    use kekule::units::{Quantity, ANGSTROM};

    for (labels, authors) in [(true, false), (false, true), (false, false)] {
        let mut builder = TopologyBuilder::new();
        let molecule = kekule::smiles::to_molecules("NCC=O")
            .unwrap()
            .pop()
            .unwrap();
        let definition = builder.add_molecule_definition(molecule.clone()).unwrap();
        let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
        let mut points = Vec::new();
        let mut keys = Vec::new();
        for sequence in [2, 1] {
            let instance = builder.add_instance(definition).unwrap();
            let residue = builder
                .hierarchy_mut()
                .add_residue(
                    chain,
                    "GLY",
                    labels.then_some(sequence),
                    authors.then(|| sequence.to_string()),
                    None,
                )
                .unwrap();
            keys.push(residue);
            for (atom, name) in molecule.atom_ids().zip(["N", "CA", "C", "O"]) {
                builder
                    .hierarchy_mut()
                    .add_atom_site(
                        residue,
                        InstanceAtomId::new(instance, atom),
                        AtomSiteMetadata::default().with_label_atom_id(Some(name.to_owned())),
                    )
                    .unwrap();
            }
            let shift = f64::from(sequence - 1) * 3.;
            points.extend(
                [(0., 0.), (1.4, 0.), (2., 1.2), (3.1, 1.3)]
                    .map(|(x, y)| Point3::new(x + shift, y, 0.)),
            );
        }
        let water = kekule::smiles::to_molecules("O").unwrap().pop().unwrap();
        let instance = builder.add_molecule(water.clone()).unwrap();
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, "HOH", None, None, None)
            .unwrap();
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, AtomId::new(0)),
                AtomSiteMetadata::default(),
            )
            .unwrap();
        points.push(Point3::new(100., 0., 0.));
        let model = Model::new(
            builder.build().unwrap(),
            Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
        )
        .unwrap();
        let result = dssp::assign(&model, DsspOptions::default()).unwrap();
        if labels || authors {
            keys.reverse();
        }
        assert_eq!(
            result
                .residues()
                .map(|residue| residue.key())
                .collect::<Vec<_>>(),
            keys
        );
        assert_eq!(result.report().non_peptide_residues(), 1);
    }
}

const ADENYLATE_KINASE_MMCIF: &str = include_str!("../fixtures/mmcif/1AKE.cif");

fn adenylate_kinase_model() -> kekule::structure::Model {
    mmcif::parse_str(ADENYLATE_KINASE_MMCIF, Default::default())
        .expect("checked-in RCSB 1AKE fixture parses")
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .expect("checked-in RCSB 1AKE fixture interprets")
        .into_model()
}

/// `copies` copies of `model` on a cubic lattice whose spacing is the
/// model's largest extent, so neighboring copies touch: residues of
/// different copies become hydrogen-bond candidates, and copies repeat chain
/// labels.
fn tiled(model: &kekule::structure::Model, copies: usize) -> kekule::structure::Model {
    use kekule::geometry::Point3;
    let positions = model.positions().values();
    let points = positions.value();
    let (low, high) = points.iter().fold(
        (
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
            Point3::new(f64::MIN, f64::MIN, f64::MIN),
        ),
        |(low, high), point| {
            (
                Point3::new(low.x.min(point.x), low.y.min(point.y), low.z.min(point.z)),
                Point3::new(
                    high.x.max(point.x),
                    high.y.max(point.y),
                    high.z.max(point.z),
                ),
            )
        },
    );
    let spacing = (high.x - low.x).max(high.y - low.y).max(high.z - low.z);
    let side = (copies as f64).cbrt().ceil() as usize;
    let mut editor = kekule::structure::ModelEditor::new();
    for copy in 0..copies {
        let [x, y, z] =
            [copy % side, copy / side % side, copy / side / side].map(|step| step as f64 * spacing);
        let mut moved = model.clone();
        moved
            .set_positions(positions.map(|points| {
                points
                    .iter()
                    .map(|point| Point3::new(point.x + x, point.y + y, point.z + z))
                    .collect::<Vec<_>>()
            }))
            .expect("translated positions stay finite");
        editor.append_model(&moved).expect("copies share no cell");
    }
    editor.finish().expect("appended copies form a model")
}

fn with_threads(options: DsspOptions, threads: usize) -> DsspOptions {
    options.with_threads(std::num::NonZeroUsize::new(threads).expect("nonzero thread count"))
}

/// One backbone atom of the analyzed residue at `index`.
fn backbone_atom(
    model: &kekule::structure::Model,
    index: usize,
    atom: &str,
) -> kekule::topology::InstanceAtomId {
    let residue = dssp::assign(model, DsspOptions::default())
        .expect("unmodified model is analyzable")
        .residues()
        .nth(index)
        .expect("residue index within the result")
        .key();
    model
        .topology()
        .residue(residue)
        .expect("analyzed residue is in the hierarchy")
        .atom_sites()
        .find(|site| site.metadata().label_atom_id() == Some(atom))
        .expect("analyzed residues have every backbone atom")
        .atom()
        .id()
}

#[test]
fn dssp_threads_return_the_serial_result() {
    for model in [
        tiled(&crambin_model(), 27),
        tiled(&adenylate_kinase_model(), 8),
    ] {
        let serial = dssp::assign(&model, DsspOptions::default()).expect("tiled model analyzes");
        // Enough residues for every stage to split across the threads.
        assert!(serial.residues().len() > 1000);
        for threads in [2, 3, 5, 8, 16] {
            let parallel = dssp::assign(&model, with_threads(DsspOptions::default(), threads))
                .expect("tiled model analyzes");
            assert_eq!(parallel, serial, "{threads} threads");
        }
    }
}

#[test]
fn dssp_threads_report_the_serial_error() {
    let model = tiled(&adenylate_kinase_model(), 8);
    let serial = dssp::assign(&model, DsspOptions::default()).unwrap();
    let pairs = serial.report().candidate_hydrogen_bond_pairs();
    let ladders = serial.report().generated_ladders();
    let residues = serial.statistics().analyzed_residues();
    for limits in [
        DsspLimits::default().with_max_residues(residues - 1),
        DsspLimits::default().with_max_candidate_pairs(0),
        DsspLimits::default().with_max_candidate_pairs(pairs / 2),
        DsspLimits::default().with_max_candidate_pairs(pairs - 1),
        DsspLimits::default().with_max_ladders(ladders / 2),
        // Bounds that are met exactly succeed.
        DsspLimits::default()
            .with_max_residues(residues)
            .with_max_candidate_pairs(pairs)
            .with_max_ladders(ladders),
    ] {
        let options = DsspOptions::default().with_limits(limits);
        let expected = dssp::assign(&model, options);
        for threads in [2, 7] {
            assert_eq!(
                dssp::assign(&model, with_threads(options, threads)),
                expected
            );
        }
    }

    // A degenerate carbonyl early in the table and a coordinate DSSP cannot
    // index later: extraction fails first, as serially.
    let (oxygen, carbon, far) = (
        backbone_atom(&model, 100, "O"),
        backbone_atom(&model, 100, "C"),
        backbone_atom(&model, 3000, "O"),
    );
    let mut carbonyl = model.clone();
    carbonyl
        .set_position(oxygen, model.position(carbon).unwrap())
        .unwrap();
    let mut broken = carbonyl.clone();
    let position = broken.position(far).unwrap().map(|mut point| {
        point.x = f32::MAX as f64;
        point
    });
    broken.set_position(far, position).unwrap();
    for (model, expected) in [(&carbonyl, "degenerate"), (&broken, "out of range")] {
        let serial = dssp::assign(model, DsspOptions::default());
        match (&serial, expected) {
            (Err(DsspError::DegenerateBackboneGeometry { .. }), "degenerate")
            | (Err(DsspError::CoordinateOutOfRange { .. }), "out of range") => {}
            _ => panic!("{expected}: {serial:?}"),
        }
        for threads in [2, 7] {
            assert_eq!(
                dssp::assign(model, with_threads(DsspOptions::default(), threads)),
                serial
            );
        }
    }
}

#[test]
#[ignore = "timing of DSSP on 100k and 1M atoms; run in release with --ignored --nocapture"]
fn dssp_threads_timing_on_a_hundred_thousand_and_a_million_atoms() {
    let base = adenylate_kinase_model();
    let available = std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get);
    for (atoms, runs) in [(100_000_usize, 30), (1_000_000, 5)] {
        let model = tiled(&base, atoms.div_ceil(base.atom_count()));
        let serial = dssp::assign(&model, DsspOptions::default()).unwrap();
        let mut report = format!(
            "{} atoms, {} residues:",
            model.atom_count(),
            serial.residues().len()
        );
        for threads in [1, 2, 4, 6, 8, available] {
            let options = with_threads(DsspOptions::default(), threads);
            let mut times = (0..runs)
                .map(|_| {
                    let start = std::time::Instant::now();
                    let result = dssp::assign(&model, options).unwrap();
                    let elapsed = start.elapsed();
                    assert_eq!(result, serial, "{threads} threads");
                    elapsed
                })
                .collect::<Vec<_>>();
            times.sort();
            report += &format!(
                " {threads}: {:.1?} (median {:.1?})",
                times[0],
                times[runs / 2]
            );
        }
        eprintln!("{report}");
    }
}

/// Runs work items last first, on the calling thread.
struct Reversed(usize);

impl dssp::DsspExecutor for Reversed {
    fn parallelism(&self) -> std::num::NonZeroUsize {
        std::num::NonZeroUsize::new(self.0).unwrap()
    }

    fn run(&self, count: usize, body: &(dyn Fn(usize) + Sync)) {
        (0..count).rev().for_each(body);
    }
}

/// Runs the odd work items on one spawned thread and the even ones, last
/// first, on another, concurrently.
struct Interleaved;

impl dssp::DsspExecutor for Interleaved {
    fn parallelism(&self) -> std::num::NonZeroUsize {
        std::num::NonZeroUsize::new(5).unwrap()
    }

    fn run(&self, count: usize, body: &(dyn Fn(usize) + Sync)) {
        std::thread::scope(|scope| {
            let odd = scope.spawn(|| (1..count).step_by(2).for_each(body));
            let even = scope.spawn(|| (0..count).step_by(2).rev().for_each(body));
            for worker in [odd, even] {
                if let Err(panic) = worker.join() {
                    std::panic::resume_unwind(panic);
                }
            }
        });
    }
}

#[test]
fn dssp_executors_return_the_serial_result_in_any_work_order() {
    let model = tiled(&adenylate_kinase_model(), 8);
    let serial = dssp::assign(&model, DsspOptions::default()).unwrap();
    let limited = DsspOptions::default().with_limits(
        DsspLimits::default()
            .with_max_candidate_pairs(serial.report().candidate_hydrogen_bond_pairs() / 2),
    );
    let executors: [&dyn dssp::DsspExecutor; 3] = [&Reversed(3), &Reversed(16), &Interleaved];
    for executor in executors {
        assert_eq!(
            dssp::assign_with(&model, DsspOptions::default(), executor).unwrap(),
            serial
        );
        assert_eq!(
            dssp::assign_with(&model, limited, executor),
            dssp::assign(&model, limited)
        );
    }
}

#[test]
fn scoped_threads_pass_a_work_item_panic_to_the_caller() {
    use dssp::{DsspExecutor, ScopedThreads};
    use std::sync::atomic::{AtomicUsize, Ordering};

    for threads in [1, 2, 4, 16] {
        let executor = ScopedThreads::new(std::num::NonZeroUsize::new(threads).unwrap());
        // Item 0 runs on the calling thread; item 5 on a spawned one when
        // there are several.
        for failing in [0, 5] {
            let ran = AtomicUsize::new(0);
            let panic = std::panic::catch_unwind(|| {
                executor.run(8, &|index| {
                    ran.fetch_add(1, Ordering::Relaxed);
                    assert_ne!(index, failing, "work item {index} fails");
                });
            })
            .expect_err("the work item panic reaches the caller");
            let message = panic
                .downcast_ref::<String>()
                .expect("formatted panic message");
            assert!(
                message.contains(&format!("work item {failing} fails")),
                "{message}"
            );
            assert!(ran.into_inner() >= 1);
        }
        let ran = AtomicUsize::new(0);
        executor.run(8, &|_| {
            ran.fetch_add(1, Ordering::Relaxed);
        });
        assert_eq!(ran.into_inner(), 8, "{threads} threads run every item once");
    }
}
