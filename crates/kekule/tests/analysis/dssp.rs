use kekule::diagnostics::{LimitExceeded, Resource};
use kekule::dssp::{
    self, DsspChainBreak, DsspError, DsspLimits, DsspOptions, DsspSecondaryStructure,
};
use kekule::macromolecular::{ModelSelection, StructureInterpretOptions};
use kekule::mmcif;
use kekule::parallel::Run;
use kekule::topology::ResidueId;
use kekule::units::{DEGREE, KILOCALORIE_PER_MOLE};

use crate::support::{Interleaved, Reversed};

const CRAMBIN_MMCIF: &str = include_str!("../fixtures/mmcif/1CRN.cif");

fn crambin_structure() -> kekule::structure::Structure {
    let document = mmcif::parse_str(CRAMBIN_MMCIF, Default::default())
        .expect("checked-in RCSB 1CRN fixture parses");
    document
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .expect("checked-in RCSB 1CRN fixture interprets")
        .into_structure()
}

#[test]
fn dssp_matches_biopython_mkdssp_4_6_1_for_crambin() {
    let structure = crambin_structure();
    let result = dssp::assign(&structure, DsspOptions::default())
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
fn dssp_is_a_coordinate_snapshot_and_does_not_mutate_the_structure() {
    let mut structure = crambin_structure();
    let before = structure.clone();
    let assigned = dssp::assign(&structure, DsspOptions::default()).expect("initial assignment");
    assert_eq!(structure, before);

    let first_atom = structure.topology().atom_ids()[0];
    let position = structure
        .position(first_atom)
        .expect("first position")
        .map(|mut point| {
            point.x += 0.25;
            point
        });
    structure
        .set_position(first_atom, position)
        .expect("finite coordinate update");

    assert_eq!(assigned.residues().count(), 46);
    assert_ne!(structure, before);
}

#[test]
fn dssp_rejects_invalid_options_and_residue_limits() {
    let structure = crambin_structure();
    assert_eq!(
        DsspOptions::default().with_min_polyproline_stretch(4),
        Err(DsspError::InvalidPolyprolineStretch { value: 4 })
    );

    let options = DsspOptions::default().with_limits(DsspLimits::default().with_max_residues(2));
    assert_eq!(
        dssp::assign(&structure, options),
        Err(DsspError::LimitExceeded(LimitExceeded::first_beyond(
            Resource::Residues,
            2
        )))
    );

    let options =
        DsspOptions::default().with_limits(DsspLimits::default().with_max_candidate_pairs(0));
    assert_eq!(
        dssp::assign(&structure, options),
        Err(DsspError::LimitExceeded(LimitExceeded::first_beyond(
            Resource::CandidatePairs,
            0
        )))
    );

    let options = DsspOptions::default().with_limits(DsspLimits::default().with_max_ladders(0));
    assert_eq!(
        dssp::assign(&structure, options),
        Err(DsspError::LimitExceeded(LimitExceeded::first_beyond(
            Resource::Ladders,
            0
        )))
    );
}

#[test]
fn dssp_rejects_coordinates_outside_its_spatial_index_range_without_panicking() {
    let mut structure = crambin_structure();
    let first_atom = structure.topology().atom_ids()[0];
    let position = structure
        .position(first_atom)
        .expect("first position")
        .map(|mut point| {
            point.x = f32::MAX as f64;
            point
        });
    structure
        .set_position(first_atom, position)
        .expect("coordinate remains finite in the structure");

    assert!(matches!(
        dssp::assign(&structure, DsspOptions::default()),
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
    let reference = dssp::assign(&crambin_structure(), DsspOptions::default()).unwrap();
    let source = CRAMBIN_MMCIF.replace("_atom_site.label_seq_id", "_atom_site.audit_sequence");
    let document = mmcif::parse_str(&source, Default::default()).unwrap();
    let structure = document
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .unwrap()
        .into_structure();
    assert!(structure
        .topology()
        .residues()
        .all(|residue| residue.label_seq_id().is_none()));
    let assigned = dssp::assign(&structure, DsspOptions::default()).unwrap();
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
    use kekule::math::Point3;
    use kekule::molecule::AtomId;
    use kekule::structure::Positions;
    use kekule::topology::{AtomSiteMetadata, InstanceAtomId};
    use kekule::units::{Quantity, ANGSTROM};

    let mut structure = crambin_structure().into_builder();
    let water = kekule::smiles::to_molecules("O").unwrap().pop().unwrap();
    let definition = structure.add_molecule_definition(water.clone()).unwrap();
    let water_position =
        Positions::new(Quantity::new(vec![Point3::new(100., 100., 100.)], ANGSTROM)).unwrap();
    let chain = structure.hierarchy_mut().add_chain("water", None).unwrap();
    let mut ignored = Vec::new();
    let count = 10_000;
    for index in 0..count {
        let instance = structure.add_instance(definition, &water_position).unwrap();
        if index % 11 == 0 {
            ignored.push(instance);
        } else {
            let residue = structure
                .hierarchy_mut()
                .add_residue(chain, "HOH", None, None, None)
                .unwrap();
            structure
                .hierarchy_mut()
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, AtomId::new(0)),
                    AtomSiteMetadata::default(),
                )
                .unwrap();
        }
    }
    let structure = structure.build().unwrap();
    let result = dssp::assign(&structure, DsspOptions::default()).unwrap();
    assert_eq!(result.report().ignored_instances(), ignored);
    assert_eq!(result.statistics().analyzed_residues(), 46);
    assert_eq!(
        result.report().non_peptide_residues(),
        count - ignored.len()
    );
}

#[test]
fn dssp_numbering_fallback_is_independent_of_unlabeled_solvent_in_the_chain() {
    use kekule::math::Point3;
    use kekule::molecule::AtomId;
    use kekule::structure::{Positions, Structure};
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
        let structure = Structure::new(
            builder.build().unwrap(),
            Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
        )
        .unwrap();
        let result = dssp::assign(&structure, DsspOptions::default()).unwrap();
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

fn adenylate_kinase_structure() -> kekule::structure::Structure {
    mmcif::parse_str(ADENYLATE_KINASE_MMCIF, Default::default())
        .expect("checked-in RCSB 1AKE fixture parses")
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .expect("checked-in RCSB 1AKE fixture interprets")
        .into_structure()
}

/// `copies` copies of `structure` on a cubic lattice whose spacing is the
/// structure's largest extent, so neighboring copies touch: residues of
/// different copies become hydrogen-bond candidates, and copies repeat chain
/// labels.
fn tiled(structure: &kekule::structure::Structure, copies: usize) -> kekule::structure::Structure {
    use kekule::math::Point3;
    let positions = structure.positions().values();
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
    let mut editor = kekule::structure::StructureEditor::new();
    for copy in 0..copies {
        let [x, y, z] =
            [copy % side, copy / side % side, copy / side / side].map(|step| step as f64 * spacing);
        let mut moved = structure.clone();
        moved
            .set_positions(positions.map(|points| {
                points
                    .iter()
                    .map(|point| Point3::new(point.x + x, point.y + y, point.z + z))
                    .collect::<Vec<_>>()
            }))
            .expect("translated positions stay finite");
        editor
            .append_structure(&moved)
            .expect("copies share no cell");
    }
    editor.finish().expect("appended copies form a structure")
}

/// Assigns `structure` on `threads` scoped threads.
fn assign_on_threads(
    structure: &kekule::structure::Structure,
    options: DsspOptions,
    threads: usize,
) -> Result<kekule::dssp::DsspResult, DsspError> {
    let executor = kekule::parallel::ScopedThreads::new(
        std::num::NonZeroUsize::new(threads).expect("nonzero thread count"),
    );
    dssp::assign_with(structure, options, Run::new().with_executor(&executor))
}

/// One backbone atom of the analyzed residue at `index`.
fn backbone_atom(
    structure: &kekule::structure::Structure,
    index: usize,
    atom: &str,
) -> kekule::topology::InstanceAtomId {
    let residue = dssp::assign(structure, DsspOptions::default())
        .expect("unmodified structure is analyzable")
        .residues()
        .nth(index)
        .expect("residue index within the result")
        .key();
    structure
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
    for structure in [
        tiled(&crambin_structure(), 27),
        tiled(&adenylate_kinase_structure(), 8),
    ] {
        let serial =
            dssp::assign(&structure, DsspOptions::default()).expect("tiled structure analyzes");
        // Enough residues for every stage to split across the threads.
        assert!(serial.residues().len() > 1000);
        for threads in [2, 3, 5, 8, 16] {
            let parallel = assign_on_threads(&structure, DsspOptions::default(), threads)
                .expect("tiled structure analyzes");
            assert_eq!(parallel, serial, "{threads} threads");
        }
    }
}

#[test]
fn dssp_threads_report_the_serial_error() {
    let structure = tiled(&adenylate_kinase_structure(), 8);
    let serial = dssp::assign(&structure, DsspOptions::default()).unwrap();
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
        let expected = dssp::assign(&structure, options);
        for threads in [2, 7] {
            assert_eq!(assign_on_threads(&structure, options, threads), expected);
        }
    }

    // A degenerate carbonyl early in the table and a coordinate DSSP cannot
    // index later: extraction fails first, as serially.
    let (oxygen, carbon, far) = (
        backbone_atom(&structure, 100, "O"),
        backbone_atom(&structure, 100, "C"),
        backbone_atom(&structure, 3000, "O"),
    );
    let mut carbonyl = structure.clone();
    carbonyl
        .set_position(oxygen, structure.position(carbon).unwrap())
        .unwrap();
    let mut broken = carbonyl.clone();
    let position = broken.position(far).unwrap().map(|mut point| {
        point.x = f32::MAX as f64;
        point
    });
    broken.set_position(far, position).unwrap();
    for (structure, expected) in [(&carbonyl, "degenerate"), (&broken, "out of range")] {
        let serial = dssp::assign(structure, DsspOptions::default());
        match (&serial, expected) {
            (Err(DsspError::DegenerateBackboneGeometry { .. }), "degenerate")
            | (Err(DsspError::CoordinateOutOfRange { .. }), "out of range") => {}
            _ => panic!("{expected}: {serial:?}"),
        }
        for threads in [2, 7] {
            assert_eq!(
                assign_on_threads(structure, DsspOptions::default(), threads),
                serial
            );
        }
    }
}

/// The structure with the backbone atom `atom` of analyzed residue `index`
/// moved by `offset` angstroms along x.
fn moved(
    structure: &kekule::structure::Structure,
    index: usize,
    atom: &str,
    offset: f64,
) -> kekule::structure::Structure {
    let atom = backbone_atom(structure, index, atom);
    let mut moved = structure.clone();
    let position = moved.position(atom).unwrap().map(|mut point| {
        point.x += offset / 10.0;
        point
    });
    moved.set_position(atom, position).unwrap();
    moved
}

#[test]
fn dssp_plan_assigns_every_structure_of_its_layout_as_assign_does() {
    use kekule::dssp::DsspPlan;

    let structure = tiled(&adenylate_kinase_structure(), 8);
    let serial = dssp::assign(&structure, DsspOptions::default()).unwrap();
    let mut carbonyl = structure.clone();
    carbonyl
        .set_position(
            backbone_atom(&structure, 100, "O"),
            structure
                .position(backbone_atom(&structure, 100, "C"))
                .unwrap(),
        )
        .unwrap();
    // Structures of the plan's layout: the structure, one whose chain breaks
    // where it did not, one with a degenerate carbonyl, and one with a
    // coordinate DSSP cannot index.
    let frames = [
        structure.clone(),
        moved(&structure, 50, "N", 30.0),
        carbonyl,
        moved(&structure, 3000, "O", f32::MAX as f64),
    ];
    let options = [
        DsspOptions::default(),
        DsspOptions::default()
            .with_min_polyproline_stretch(2)
            .unwrap(),
        DsspOptions::default().with_limits(
            DsspLimits::default()
                .with_max_candidate_pairs(serial.report().candidate_hydrogen_bond_pairs() / 2),
        ),
        DsspOptions::default().with_limits(
            DsspLimits::default().with_max_residues(serial.statistics().analyzed_residues() / 2),
        ),
    ];
    let threads = kekule::parallel::ScopedThreads::new(std::num::NonZeroUsize::new(3).unwrap());
    for options in options {
        let plan = DsspPlan::new(structure.shared_topology(), options);
        assert_eq!(plan.options(), options);
        assert!(plan.topology().shares_layout(structure.topology()));
        for (index, frame) in frames.iter().enumerate() {
            let expected = dssp::assign(frame, options);
            assert_eq!(plan.assign(frame), expected, "frame {index}, {options:?}");
            assert_eq!(
                plan.assign_with(frame, Run::new().with_executor(&Reversed(4))),
                expected,
                "frame {index}, {options:?}"
            );
            assert_eq!(
                plan.assign_with(frame, Run::new().with_executor(&threads)),
                expected,
                "frame {index}, {options:?}, 3 threads"
            );
        }
    }
    assert_ne!(dssp::assign(&frames[1], DsspOptions::default()), Ok(serial));

    let reversed = DsspPlan::new(structure.shared_topology(), DsspOptions::default());
    assert_eq!(
        reversed.assign_with(&structure, Run::new().with_executor(&Reversed(5))),
        dssp::assign(&structure, DsspOptions::default())
    );
}

#[test]
fn dssp_runs_report_their_stages_and_stop_when_cancelled() {
    use kekule::diagnostics::{Cancellation, Cancelled};
    use kekule::dssp::DsspPlan;
    use kekule::parallel::Progress;
    use std::sync::Mutex;

    let structure = adenylate_kinase_structure();
    let serial = dssp::assign(&structure, DsspOptions::default());
    let reports = Mutex::new(Vec::new());
    let progress = |progress: Progress| reports.lock().unwrap().push(progress);
    let run = Run::new()
        .with_executor(&Reversed(3))
        .with_progress(&progress);
    assert_eq!(
        dssp::assign_with(&structure, DsspOptions::default(), run),
        serial
    );
    let stages = std::mem::take(&mut *reports.lock().unwrap());
    assert_eq!(
        stages,
        (0..=5)
            .map(|done| Progress::new(done, 5))
            .collect::<Vec<_>>()
    );

    let cancellation = Cancellation::new(|| true);
    let cancelled = Run::new().with_cancellation(&cancellation);
    let error = dssp::assign_with(&structure, DsspOptions::default(), cancelled).unwrap_err();
    assert!(matches!(error, DsspError::Cancelled(_)), "{error:?}");
    assert_eq!(error.to_string(), "DSSP assignment stopped.");
    assert!(std::error::Error::source(&error).is_some_and(|source| source.is::<Cancelled>()));
    let plan = DsspPlan::new(structure.shared_topology(), DsspOptions::default());
    assert!(matches!(
        plan.assign_with(&structure, cancelled),
        Err(DsspError::Cancelled(_))
    ));
}

#[test]
fn dssp_plan_rejects_structures_of_another_layout() {
    use kekule::dssp::DsspPlan;

    let structure = adenylate_kinase_structure();
    let plan = DsspPlan::new(structure.shared_topology(), DsspOptions::default());
    // An independently read copy has equal contents but not the same layout.
    let copy = adenylate_kinase_structure();
    assert!(copy.topology().same_layout(structure.topology()));
    for other in [copy, crambin_structure()] {
        assert_eq!(plan.assign(&other), Err(DsspError::TopologyMismatch));
    }
    assert_eq!(
        DsspError::TopologyMismatch.to_string(),
        "The structure's topology layout differs from the DSSP plan's."
    );
}

#[test]
#[ignore = "timing of DSSP and DSSP plans on 100k and 1M atoms; run in release with --ignored --nocapture"]
fn dssp_threads_timing_on_a_hundred_thousand_and_a_million_atoms() {
    let base = adenylate_kinase_structure();
    let available = std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get);
    for (atoms, runs) in [(100_000_usize, 30), (1_000_000, 5)] {
        let structure = tiled(&base, atoms.div_ceil(base.atom_count()));
        let serial = dssp::assign(&structure, DsspOptions::default()).unwrap();
        let mut report = format!(
            "{} atoms, {} residues:",
            structure.atom_count(),
            serial.residues().len()
        );
        for threads in [1, 2, 4, 6, 8, available] {
            let options = DsspOptions::default();
            let executor = kekule::parallel::ScopedThreads::new(
                std::num::NonZeroUsize::new(threads).expect("nonzero thread count"),
            );
            let run = Run::new().with_executor(&executor);
            let plan = kekule::dssp::DsspPlan::new(structure.shared_topology(), options);
            for (name, assign) in [
                (
                    "",
                    &(|| dssp::assign_with(&structure, options, run)) as &dyn Fn() -> _,
                ),
                (" plan", &|| plan.assign_with(&structure, run)),
            ] {
                let mut times = (0..runs)
                    .map(|_| {
                        let start = std::time::Instant::now();
                        let result = assign().unwrap();
                        let elapsed = start.elapsed();
                        assert_eq!(result, serial, "{threads} threads{name}");
                        elapsed
                    })
                    .collect::<Vec<_>>();
                times.sort();
                report += &format!(
                    " {threads}{name}: {:.1?} (median {:.1?})",
                    times[0],
                    times[runs / 2]
                );
            }
        }
        eprintln!("{report}");
    }
}

#[test]
fn dssp_executors_return_the_serial_result_in_any_work_order() {
    let structure = tiled(&adenylate_kinase_structure(), 8);
    let serial = dssp::assign(&structure, DsspOptions::default()).unwrap();
    let limited = DsspOptions::default().with_limits(
        DsspLimits::default()
            .with_max_candidate_pairs(serial.report().candidate_hydrogen_bond_pairs() / 2),
    );
    let executors: [&dyn kekule::parallel::Executor; 3] =
        [&Reversed(3), &Reversed(16), &Interleaved];
    for executor in executors {
        let run = Run::new().with_executor(executor);
        assert_eq!(
            dssp::assign_with(&structure, DsspOptions::default(), run).unwrap(),
            serial
        );
        assert_eq!(
            dssp::assign_with(&structure, limited, run),
            dssp::assign(&structure, limited)
        );
    }
}

#[test]
fn scoped_threads_pass_a_work_item_panic_to_the_caller() {
    use kekule::parallel::{Executor, ScopedThreads};
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
