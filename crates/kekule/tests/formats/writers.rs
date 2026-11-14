use kekule::diagnostics::Countable;
use kekule::diagnostics::LengthMismatch;
use std::sync::Arc;
use std::{io, io::Write};

use kekule::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, FormatStage, GeometryField,
    UnrepresentablePolicy, WriteOmission,
};
use kekule::macromolecular::EntityKind;
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::mmcif::{self, MmcifEntityClassifications, MmcifWriteErrorKind, MmcifWriteOptions};
use kekule::molfile::{self, MolfileWriteOptions, MolfileWriteVersion};
use kekule::pdb::{self, PdbWriteOptions};
use kekule::sdf::{self, SdfDataField, SdfRecordInterpretation, SdfWriteOptions};
use kekule::smiles::{self, SmilesWriteMode, SmilesWriteOptions};
use kekule::structure::{
    Distribution, Ensemble, Forces, Geometry, Order, Partition, Positions, Samples, StateData,
    Structure, StructureSet, StructureView, TimeAxis, Velocities,
};
use kekule::topology::{Topology, TopologyBuilder};
use kekule::traj::FrameView;
use kekule::units::{
    Quantity, ANGSTROM, KELVIN, KILOJOULE_PER_MOLE, KILOJOULE_PER_MOLE_PER_NANOMETER, NANOMETER,
    NANOMETER_PER_PICOSECOND, PICOSECOND,
};

fn molecule(source: &str) -> kekule::molecule::Molecule {
    let mut molecules = smiles::to_molecules(source).expect("SMILES interprets");
    assert_eq!(molecules.len(), 1);
    molecules.pop().unwrap()
}

fn structure(source: &str, points: &[[f64; 3]]) -> Structure {
    let topology = smiles::to_topology(source).expect("topology interprets");
    let positions = Positions::new(Quantity::new(
        points
            .iter()
            .map(|point| Point3::new(point[0], point[1], point[2]))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap();
    Structure::new(topology, positions).unwrap()
}

fn classifications(structure: &Structure) -> MmcifEntityClassifications {
    let mut classifications = MmcifEntityClassifications::new();
    for instance in structure
        .topology()
        .instances()
        .map(|molecule| molecule.id())
    {
        classifications
            .insert(instance, EntityKind::NonPolymer)
            .unwrap();
    }
    classifications
}

struct BoundedChunkWriter {
    bytes: Vec<u8>,
    maximum_chunk: usize,
}

impl Write for BoundedChunkWriter {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        if buffer.len() > self.maximum_chunk {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "writer received a buffered document-sized chunk",
            ));
        }
        self.bytes.extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn smiles_writes_molecules_and_repeated_topology_instances_in_authoritative_order() {
    let water = molecule("O");
    let sodium = molecule("[Na+]");
    let mut builder = TopologyBuilder::new();
    let water_definition = builder.add_molecule_definition(water.clone()).unwrap();
    let sodium_definition = builder.add_molecule_definition(sodium.clone()).unwrap();
    builder.add_instance(water_definition).unwrap();
    builder.add_instance(water_definition).unwrap();
    builder.add_instance(sodium_definition).unwrap();
    let topology = builder.build().unwrap();

    assert_eq!(
        smiles::write(&water, SmilesWriteOptions::default()).unwrap(),
        "O"
    );
    assert_eq!(
        smiles::write(&topology, SmilesWriteOptions::default()).unwrap(),
        "O.O.[Na+]"
    );

    let chiral = molecule("F[C@H](Cl)Br");
    assert_eq!(
        smiles::write(
            &chiral,
            SmilesWriteOptions::default().with_mode(SmilesWriteMode::Isomeric)
        )
        .unwrap(),
        smiles::write(&chiral, smiles::SmilesWriteOptions::isomeric()).unwrap()
    );
    assert_eq!(
        smiles::write(
            &molecule("OCC"),
            SmilesWriteOptions::default().with_mode(SmilesWriteMode::Canonical)
        )
        .unwrap(),
        "CCO"
    );
}

#[test]
fn smiles_write_to_emits_each_record_as_one_line() {
    let mut output = Vec::new();
    for source in ["OCC", "[Na+].[Cl-]"] {
        smiles::write_to(
            &mut output,
            &smiles::to_topology(source).unwrap(),
            SmilesWriteOptions::canonical(),
        )
        .unwrap();
    }
    assert_eq!(String::from_utf8(output).unwrap(), "CCO\n[Cl-].[Na+]\n");
}

#[test]
fn molfile_rejects_asserted_stereo_its_coordinates_contradict() {
    // These coordinates are F[C@H](Cl)Br; a 3D reader takes the configuration
    // from them, so writing the opposite assertion would invert it.
    let points = [
        [0.6, 0.6, 0.6],
        [0.0, 0.0, 0.0],
        [0.6, -0.6, -0.6],
        [-0.6, 0.6, -0.6],
    ];
    for version in [MolfileWriteVersion::V2000, MolfileWriteVersion::V3000] {
        molfile::write(
            &structure("F[C@H](Cl)Br", &points),
            MolfileWriteOptions::default().with_version(version),
        )
        .unwrap();
        let error = molfile::write(
            &structure("F[C@@H](Cl)Br", &points),
            MolfileWriteOptions::default().with_version(version),
        )
        .unwrap_err();
        assert_eq!(
            (error.kind(), error.message()),
            (
                DiagnosticKind::Unsupported,
                "The specified tetrahedral stereo at atom 1 contradicts the written coordinates."
            ),
            "{version:?}"
        );
    }
}

#[test]
fn molfile_v3000_continues_records_past_80_columns() {
    // CTfile lines hold at most 80 characters; a longer V3000 record ends
    // each line but its last with `-` and breaks after a space, because
    // readers trim leading whitespace from continuation lines.
    let charged = structure(
        "[13CH3-]",
        &[[-1.2345678901234567, -2.345678901234568, -3.456789012345679]],
    );
    let text = molfile::write(
        &charged,
        MolfileWriteOptions::default().with_version(MolfileWriteVersion::V3000),
    )
    .unwrap()
    .0;
    let stored = charged
        .positions()
        .get(0)
        .unwrap()
        .value_in(ANGSTROM)
        .unwrap();
    let atom_block = text
        .lines()
        .skip_while(|line| *line != "M  V30 BEGIN ATOM")
        .skip(1)
        .take_while(|line| *line != "M  V30 END ATOM")
        .collect::<Vec<_>>();
    // `CHG=-1` would end the first line at column 81, so it starts the second.
    let first = format!("M  V30 1 C {} {} {} 0 -", stored.x, stored.y, stored.z);
    assert_eq!(first.len(), 74);
    assert_eq!(
        atom_block,
        [first.as_str(), "M  V30 CHG=-1 MASS=13 VAL=3"],
        "{text}"
    );
    let reread = molfile::parse_str(&text, Default::default()).unwrap();
    let atom = &reread.atoms()[0];
    assert_eq!(atom.formal_charge(), -1);
    assert_eq!(atom.isotope(), Some(13));
    assert_eq!(atom.position(), stored);
}

#[test]
fn molfile_chiral_flag_states_whether_stereo_is_absolute() {
    // An unset flag makes BIOVIA-convention readers, and RDKit's
    // TranslateChiralFlagToStereoGroups, treat ungrouped centers as one
    // racemic AND group; Kekule reads ungrouped stereo as absolute.
    let layout = [
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, -0.5, 0.0],
        [-1.0, -0.5, 0.0],
    ];
    for (source, version, chiral) in [
        ("CC(N)O", MolfileWriteVersion::V2000, false),
        ("C[C@H](N)O", MolfileWriteVersion::V2000, true),
        ("C[C@H](N)O", MolfileWriteVersion::V3000, true),
        ("C[C@H](N)O |a:1|", MolfileWriteVersion::V3000, true),
        ("C[C@H](N)O |&1:1|", MolfileWriteVersion::V3000, false),
        ("C[C@H](N)O |o1:1|", MolfileWriteVersion::V3000, false),
    ] {
        let text = molfile::write(
            &structure(source, &layout),
            MolfileWriteOptions::default().with_version(version),
        )
        .unwrap()
        .0;
        let document = molfile::parse_str(&text, Default::default()).unwrap();
        assert_eq!(
            document.chiral_flag(),
            chiral,
            "{source} {version:?}\n{text}"
        );
    }
}

#[test]
fn automatic_molfile_version_writes_what_explicit_v3000_writes_when_v2000_cannot() {
    let layout = [
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, -0.5, 0.0],
        [-1.0, -0.5, 0.0],
    ];
    let enhanced = structure("C[C@H](N)O |&1:1|", &layout);
    let sodium = molecule("[Na+]");
    let large = Structure::new(
        Topology::from_molecules(vec![sodium; 1_000]).unwrap(),
        Positions::zeros(1_000),
    )
    .unwrap();
    for structure in [&enhanced, &large] {
        let explicit = molfile::write(
            structure,
            MolfileWriteOptions::default().with_version(MolfileWriteVersion::V3000),
        )
        .unwrap()
        .0;
        assert_eq!(
            molfile::write(structure, MolfileWriteOptions::default())
                .unwrap()
                .0,
            explicit
        );
    }
}

#[test]
fn molfile_structure_flattens_disconnected_topology_and_preserves_coordinates() {
    let structure = structure(
        "CO.O",
        &[[1.25, 2.5, 3.75], [2.0, 2.5, 3.75], [-4.0, 0.5, 8.0]],
    );
    let text = molfile::write(&structure, MolfileWriteOptions::default())
        .unwrap()
        .0;
    assert!(text.contains("V2000"));

    let interpreted = (molfile::parse_str(&text, Default::default()).unwrap())
        .interpret(Default::default())
        .unwrap();
    let round_trip_structure = interpreted.into_structure();
    assert_eq!(round_trip_structure.topology().instance_count(), 2);
    for (index, expected) in [
        Point3::new(1.25, 2.5, 3.75),
        Point3::new(2.0, 2.5, 3.75),
        Point3::new(-4.0, 0.5, 8.0),
    ]
    .into_iter()
    .enumerate()
    {
        let actual = round_trip_structure
            .positions()
            .get(index)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap();
        assert!((actual.x - expected.x).abs() <= 1e-4);
        assert!((actual.y - expected.y).abs() <= 1e-4);
        assert!((actual.z - expected.z).abs() <= 1e-4);
    }
}

#[test]
fn molfile_version_policy_promotes_counts_and_explicit_v2000_fails() {
    let sodium = molecule("[Na+]");
    let molecules = vec![sodium; 1_000];
    let topology = Topology::from_molecules(molecules.clone()).unwrap();
    let structure = Structure::new(topology, Positions::zeros(1_000)).unwrap();

    let error = molfile::write(
        &structure,
        MolfileWriteOptions::default().with_version(MolfileWriteVersion::V2000),
    )
    .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);

    let promoted = molfile::write(&structure, MolfileWriteOptions::default())
        .unwrap()
        .0;
    assert!(promoted.contains("V3000"));
    assert!(promoted.contains("M  V30 COUNTS 1000 0 0 0 0"));
    assert_eq!(
        (molfile::parse_str(&promoted, Default::default()).unwrap())
            .interpret(Default::default())
            .unwrap()
            .into_structure()
            .topology()
            .instance_count(),
        1_000
    );
}

const DRAWING_HEADER: &str = "  kekule            2D";
const CONFORMATION_HEADER: &str = "  kekule            3D";

#[test]
fn molfile_header_dimension_code_states_the_emitted_coordinate_kind() {
    // CTfile header line 2 is `IIPPPPPPPPMMDDYYHHmmdd...`: blank initials,
    // the program name in columns 3-10, a blank timestamp, and the
    // dimensional code in columns 21-22. Readers such as RDKit take any
    // nonzero z as 3D, so the code is 2D exactly when every emitted z is zero.
    use molfile::MolfileSource;
    use MolfileWriteVersion::{V2000, V3000};
    let conformer = structure("CO", &[[0.0, 0.0, 0.0], [1.2, 0.3, -0.4]]);
    let drawing = structure("CO", &[[0.0, 0.0, 0.0], [1.2, 0.3, 0.0]]);
    // An unmarked flat record keeps its plane, so its z is emitted as given.
    let raised_plane = structure("CO", &[[0.0, 0.0, 1.5], [1.2, 0.3, 1.5]]);
    // Both z values print as 0.0000 in V2000 but stay distinct in V3000.
    let sub_precision_relief = structure("CO", &[[0.0, 0.0, 0.0], [1.2, 0.0, 0.00004]]);
    let coordinate_free = molecule("CO");
    for (source, version, expected) in [
        (MolfileSource::from(&conformer), V2000, CONFORMATION_HEADER),
        ((&conformer).into(), V3000, CONFORMATION_HEADER),
        ((&drawing).into(), V2000, DRAWING_HEADER),
        ((&drawing).into(), V3000, DRAWING_HEADER),
        ((&raised_plane).into(), V2000, CONFORMATION_HEADER),
        ((&raised_plane).into(), V3000, CONFORMATION_HEADER),
        ((&sub_precision_relief).into(), V2000, DRAWING_HEADER),
        ((&sub_precision_relief).into(), V3000, CONFORMATION_HEADER),
        ((&coordinate_free).into(), V2000, DRAWING_HEADER),
        ((&coordinate_free).into(), V3000, DRAWING_HEADER),
    ] {
        let text = molfile::write(source, MolfileWriteOptions::default().with_version(version))
            .unwrap()
            .0;
        assert_eq!(text.lines().nth(1), Some(expected), "{version:?}\n{text}");
        assert_eq!(text.lines().nth(2), Some(""));
    }

    let text = sdf::write(
        [
            sdf::SdfRecordSource::structure(&conformer).with_title("conformer"),
            sdf::SdfRecordSource::structure(&drawing).with_title("drawing"),
        ],
        SdfWriteOptions::default(),
    )
    .unwrap()
    .0;
    let headers = text
        .split("$$$$\n")
        .filter(|record| !record.is_empty())
        .map(|record| record.lines().take(2).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    assert_eq!(
        headers,
        [
            ["conformer", CONFORMATION_HEADER],
            ["drawing", DRAWING_HEADER]
        ]
    );
}

#[test]
fn molfile_writes_a_marked_flat_drawing_on_the_z_zero_plane() {
    // RDKit reads any nonzero z as 3D and then ignores bond marks, so a raised
    // drawing whose wedge carries stereo is emitted on z = 0 with XY intact.
    let layout = [
        [0.0, 1.0, 1.5],
        [0.0, 0.0, 1.5],
        [1.0, -0.5, 1.5],
        [-1.0, -0.5, 1.5],
    ];
    let raised_drawing = structure("C[C@H](N)O", &layout);
    let canonical =
        |topology: &Topology| smiles::write(topology, SmilesWriteOptions::canonical()).unwrap();
    assert_eq!(canonical(raised_drawing.topology()), "C[C@H](N)O");
    for version in [MolfileWriteVersion::V2000, MolfileWriteVersion::V3000] {
        let text = molfile::write(
            &raised_drawing,
            MolfileWriteOptions::default().with_version(version),
        )
        .unwrap()
        .0;
        assert_eq!(
            text.lines().nth(1),
            Some(DRAWING_HEADER),
            "{version:?}\n{text}"
        );
        let document = molfile::parse_str(&text, Default::default()).unwrap();
        let z = document
            .atoms()
            .iter()
            .map(|atom| atom.position().z)
            .collect::<Vec<_>>();
        assert_eq!(z, [0.0; 4], "{version:?}\n{text}");
        let marked = document
            .bonds()
            .iter()
            .filter(|bond| bond.stereo().is_some())
            .count();
        assert_eq!(marked, 1, "{version:?}\n{text}");
        let written = document
            .interpret(Default::default())
            .unwrap()
            .into_structure();
        assert_eq!(
            canonical(written.topology()),
            "C[C@H](N)O",
            "{version:?}\n{text}"
        );
    }
}

fn interpret_records(text: &str) -> Vec<SdfRecordInterpretation> {
    sdf::parse_str(text, Default::default())
        .unwrap()
        .records()
        .iter()
        .map(|record| record.interpret(Default::default()).unwrap())
        .collect()
}

#[test]
fn sdf_writes_independent_structures_and_set_geometries_in_order() {
    let first = structure("C", &[[1.0, 0.0, 0.0]]);
    let second = structure("O", &[[2.0, 0.0, 0.0]]);
    let text = sdf::write(&[first.clone(), second.clone()], SdfWriteOptions::default())
        .unwrap()
        .0;
    let interpreted = interpret_records(&text);
    assert_eq!(interpreted.len(), 2);
    assert_eq!(
        interpreted[0]
            .structure()
            .topology()
            .atoms()
            .next()
            .unwrap()
            .element()
            .symbol(),
        "C"
    );
    assert_eq!(
        interpreted[1]
            .structure()
            .topology()
            .atoms()
            .next()
            .unwrap()
            .element()
            .symbol(),
        "O"
    );

    let shared = first.shared_topology();
    let later = Structure::new(
        shared,
        Positions::new(Quantity::new([Point3::new(9.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    let set = StructureSet::from_structures([first, later]).unwrap();
    let set_text = sdf::write(&set, SdfWriteOptions::default()).unwrap().0;
    let records = interpret_records(&set_text);
    assert_eq!(records.len(), 2);
    // An unordered set's geometries carry no data fields.
    assert_eq!(data_fields(&records), [Vec::new(), Vec::new()]);
    assert!(records.iter().all(|record| record.title().is_empty()));
    assert!(
        (records[0]
            .structure()
            .positions()
            .get(0)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            .x
            - 1.0)
            .abs()
            <= 1e-9
    );
    assert!(
        (records[1]
            .structure()
            .positions()
            .get(0)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            .x
            - 9.0)
            .abs()
            <= 1e-9
    );
}

/// Two methane conformers at 1 and 9 Å along x.
fn conformer_set() -> StructureSet {
    let first = structure("C", &[[1.0, 0.0, 0.0]]);
    let later = Structure::new(
        first.shared_topology(),
        Positions::new(Quantity::new([Point3::new(9.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    StructureSet::from_structures([first, later]).unwrap()
}

/// The conformers as one state each, weighted 1:3 (log weights 0 and ln 3),
/// with state energies 0 and 5 kJ/mol.
fn weighted_ensemble() -> Ensemble {
    let states = [0.0, 5.0]
        .map(|energy| {
            StateData::new()
                .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
                .unwrap()
        })
        .to_vec();
    let distribution = Distribution::new(Samples::All, Partition::state_per_sample(2))
        .unwrap()
        .with_log_weights(vec![0.0, 3.0_f64.ln()])
        .unwrap()
        .with_states(states)
        .unwrap();
    Ensemble::new(conformer_set(), distribution).unwrap()
}

fn data_fields(records: &[SdfRecordInterpretation]) -> Vec<Vec<(String, String)>> {
    records
        .iter()
        .map(|record| {
            record
                .data_fields()
                .iter()
                .map(|field| (field.name().to_owned(), field.value().to_owned()))
                .collect::<Vec<_>>()
        })
        .collect()
}

#[test]
fn sdf_never_writes_a_distribution() {
    // The geometries of an ensemble's set write no weight, state, energy, or
    // temperature, whatever the distribution over them.
    let mut ensemble = weighted_ensemble();
    ensemble
        .distribution_mut()
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    let (text, report) = sdf::write(ensemble.set(), SdfWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    assert!(!text.contains("kekule."), "{text}");
    let records = interpret_records(&text);
    assert_eq!(data_fields(&records), [Vec::new(), Vec::new()]);
    let plain = sdf::write(&conformer_set(), SdfWriteOptions::default())
        .unwrap()
        .0;
    assert_eq!(text, plain);
    // So does each geometry alone, and a Molfile of it loses nothing.
    let (single, report) =
        sdf::write([ensemble.set().get(1).unwrap()], SdfWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    assert!(!single.contains("kekule."), "{single}");
    let (_, report) = molfile::write(
        molfile::MolfileSource::structure(&ensemble.set().get(1).unwrap()),
        MolfileWriteOptions::default(),
    )
    .unwrap();
    assert!(report.is_empty());
}

#[test]
fn sdf_reads_the_old_distribution_fields_as_ordinary_data() {
    // Records an earlier writer gave a distribution, a member marker, and a
    // malformed weight.
    let old = [
        ("kekule.weight", "heavy"),
        ("kekule.state", "1"),
        ("kekule.representative", "true"),
        ("kekule.energy.kJ_per_mol", "-2.5"),
        ("kekule.temperature.K", "300"),
    ];
    let outside = [("kekule.member", "false"), ("kekule.weight", "-3")];
    let fields = |pairs: &[(&str, &str)]| {
        pairs
            .iter()
            .map(|(name, value)| SdfDataField::new(*name, *value))
            .collect::<Vec<_>>()
    };
    let set = conformer_set();
    let first = fields(&old);
    let second = fields(&outside);
    let records = [
        sdf::SdfRecordSource::from(set.get(0).unwrap()).with_data_fields(&first),
        sdf::SdfRecordSource::from(set.get(1).unwrap()).with_data_fields(&second),
    ];
    // The writer reserves none of the names: they are a caller's fields.
    let (text, report) = sdf::write(records, SdfWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    let records = interpret_records(&text);
    let owned = |pairs: &[(&str, &str)]| {
        pairs
            .iter()
            .map(|(name, value)| ((*name).to_owned(), (*value).to_owned()))
            .collect::<Vec<_>>()
    };
    assert_eq!(data_fields(&records), [owned(&old), owned(&outside)]);
    // The records form a structure set of both geometries, which has no
    // distribution to state; the fields are neither read nor validated.
    let read = sdf::structure_set(&records).unwrap();
    assert_eq!(read.len(), 2);
    assert_eq!(read.order(), None);
    for (left, right) in read.iter().zip(set.iter()) {
        assert_eq!(left.positions(), right.positions());
    }
    // Writing the set again drops the fields, which belonged to the records.
    let rewritten = sdf::write(&read, SdfWriteOptions::default()).unwrap().0;
    assert!(!rewritten.contains("kekule."), "{rewritten}");
}

/// One-atom methane frames with exactly one piece of frame state each:
/// none, time, MD step, velocities, forces.
struct StatefulFrames {
    topology: Arc<Topology>,
    geometries: Vec<Geometry>,
}

impl StatefulFrames {
    fn new() -> Self {
        let topology = structure("C", &[[0.0; 3]]).shared_topology();
        let mut geometries = vec![Geometry::new(Positions::zeros(1)); 5];
        let zero = || vec![Vector3::new(0.0, 0.0, 0.0)];
        let velocities =
            Velocities::from_vec(Quantity::new(zero(), NANOMETER_PER_PICOSECOND)).unwrap();
        geometries[3].set_velocities(Some(velocities)).unwrap();
        let forces =
            Forces::from_vec(Quantity::new(zero(), KILOJOULE_PER_MOLE_PER_NANOMETER)).unwrap();
        geometries[4].set_forces(Some(forces)).unwrap();
        Self {
            topology,
            geometries,
        }
    }

    /// The frames, the second at 2 ps and the third at MD step 7.
    fn views(&self) -> Vec<FrameView<'_>> {
        self.geometries
            .iter()
            .enumerate()
            .map(|(index, geometry)| {
                let time = (index == 1).then(|| Quantity::new(2.0, PICOSECOND));
                let step = (index == 2).then_some(7);
                FrameView::new(&self.topology, geometry, time, step).unwrap()
            })
            .collect()
    }
}

/// The structure `view` borrows, without the velocities and forces that
/// `to_structure` keeps: what a writer gives once a geometry's state is
/// omitted.
fn bare_structure(view: StructureView<'_>) -> Structure {
    let mut structure = view.to_structure();
    let mut geometry = structure.geometry_mut();
    geometry.take_velocities();
    geometry.take_forces();
    structure
}

/// Writes one structure view as a Molfile, an SDF record and an mmCIF
/// block, each with geometry state set to `policy`.
fn write_structure(
    view: StructureView<'_>,
    policy: UnrepresentablePolicy,
) -> [Result<String, FormatError>; 3] {
    [
        molfile::write(
            molfile::MolfileSource::structure(&view),
            MolfileWriteOptions::default().with_geometry_metadata(policy),
        )
        .map(|(text, _)| text),
        sdf::write(
            [sdf::SdfRecordSource::structure(&view)],
            SdfWriteOptions::from(MolfileWriteOptions::default().with_geometry_metadata(policy)),
        )
        .map(|(text, _)| text),
        mmcif::write(
            [mmcif::MmcifBlockSource::structure(&view)],
            MmcifWriteOptions::default().with_geometry_metadata(policy),
        )
        .map(|(text, _)| text),
    ]
}

#[test]
fn structure_writers_fail_on_geometry_and_frame_state_unless_told_to_omit_it() {
    use GeometryField::{Forces as ForceField, Step, Time, Velocities as VelocityField};
    let frames = StatefulFrames::new();
    let frames = frames.views();
    // Expected detail per writer (Molfile, SDF, mmCIF).
    let rows: [(StructureView<'_>, [Option<GeometryField>; 3]); 5] = [
        (frames[0].as_structure_view(), [None; 3]),
        (frames[1].as_structure_view(), [Some(Time); 3]),
        (frames[2].as_structure_view(), [Some(Step); 3]),
        (frames[3].as_structure_view(), [Some(VelocityField); 3]),
        (frames[4].as_structure_view(), [Some(ForceField); 3]),
    ];
    let formats = [Format::Molfile, Format::Sdf, Format::Mmcif];
    for (row, (view, expected)) in rows.into_iter().enumerate() {
        let bare = bare_structure(view);
        let plain = write_structure(bare.as_structure_view(), UnrepresentablePolicy::Fail);
        let failed = write_structure(view, UnrepresentablePolicy::Fail);
        let omitted = write_structure(view, UnrepresentablePolicy::Omit);
        // A standalone structure keeps the frame's velocities and forces,
        // and the writers hold them to the same policy.
        let structure = view.to_structure();
        let standalone_failed =
            write_structure(structure.as_structure_view(), UnrepresentablePolicy::Fail);
        let standalone_omitted =
            write_structure(structure.as_structure_view(), UnrepresentablePolicy::Omit);
        for (writer, field) in expected.into_iter().enumerate() {
            let format = formats[writer];
            match field {
                Some(field) => {
                    let error = failed[writer].as_ref().unwrap_err();
                    assert_eq!(
                        (error.format(), error.kind(), error.detail()),
                        (
                            format,
                            DiagnosticKind::Unsupported,
                            Some(&FormatErrorDetail::UnrepresentableState(field)),
                        ),
                        "row {row}, {format}"
                    );
                    // Omitting the state writes exactly the bare structure.
                    assert_eq!(
                        omitted[writer].as_ref().unwrap(),
                        plain[writer].as_ref().unwrap(),
                        "row {row}, {format}"
                    );
                    if matches!(field, VelocityField | ForceField) {
                        let error = standalone_failed[writer].as_ref().unwrap_err();
                        assert_eq!(
                            (error.kind(), error.detail()),
                            (
                                DiagnosticKind::Unsupported,
                                Some(&FormatErrorDetail::UnrepresentableState(field)),
                            ),
                            "row {row}, {format}"
                        );
                    } else {
                        assert_eq!(
                            standalone_failed[writer].as_ref().unwrap(),
                            plain[writer].as_ref().unwrap(),
                            "row {row}, {format}"
                        );
                    }
                    assert_eq!(
                        standalone_omitted[writer].as_ref().unwrap(),
                        plain[writer].as_ref().unwrap(),
                        "row {row}, {format}"
                    );
                }
                None => assert_eq!(
                    failed[writer].as_ref().unwrap(),
                    plain[writer].as_ref().unwrap(),
                    "row {row}, {format}"
                ),
            }
        }
    }
    let error = sdf::write(frames.clone(), SdfWriteOptions::default()).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Cannot write the SDF output in record 2. The time of the geometry has no SDF representation."
    );
    let error = molfile::write(
        molfile::MolfileSource::structure(&frames[2].as_structure_view()),
        MolfileWriteOptions::default(),
    )
    .unwrap_err();
    assert_eq!(
        error.to_string(),
        "Cannot write the Molfile output. The MD step of the geometry has no Molfile representation."
    );
}

/// A three-step optimization history and a two-point distance scan of one
/// methane carbon.
fn methane_series() -> (StructureSet, StructureSet) {
    let topology = structure("C", &[[0.0; 3]]).shared_topology();
    let geometries = |count| (0..count).map(|_| Geometry::new(Positions::zeros(1)));
    let mut history = StructureSet::from_geometries(Arc::clone(&topology), geometries(3)).unwrap();
    history.set_order(Some(Order::Step(vec![0, 4, 9]))).unwrap();
    let mut scan = StructureSet::from_geometries(topology, geometries(2)).unwrap();
    scan.set_order(Some(Order::Coordinate {
        unit: ANGSTROM,
        values: vec![1.25, 1.5],
    }))
    .unwrap();
    (history, scan)
}

#[test]
fn sdf_ordered_geometries_carry_their_step_or_coordinate_through_data_fields() {
    let (history, scan) = methane_series();
    let field = |name: &str, value: &str| (name.to_owned(), value.to_owned());
    // A methane carbon at the origin, then `fields`.
    let record = |fields: &[(&str, &str)]| {
        let fields = fields
            .iter()
            .map(|(name, value)| format!(">  <{name}>\n{value}\n\n"))
            .collect::<String>();
        format!(
            "\n  kekule            2D\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    \
             0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n{fields}$$$$\n"
        )
    };
    let text = sdf::write(&history, SdfWriteOptions::default()).unwrap().0;
    assert_eq!(
        text,
        ["0", "4", "9"]
            .map(|step| record(&[(sdf::STEP_FIELD, step)]))
            .concat()
    );
    assert_eq!(
        data_fields(&interpret_records(&text)),
        [
            vec![field(sdf::STEP_FIELD, "0")],
            vec![field(sdf::STEP_FIELD, "4")],
            vec![field(sdf::STEP_FIELD, "9")],
        ]
    );
    let text = sdf::write(&scan, SdfWriteOptions::default()).unwrap().0;
    assert_eq!(
        text,
        ["1.25", "1.5"]
            .map(|coordinate| {
                record(&[
                    (sdf::COORDINATE_FIELD, coordinate),
                    (sdf::COORDINATE_UNIT_FIELD, "A"),
                ])
            })
            .concat()
    );
    assert_eq!(
        data_fields(&interpret_records(&text)),
        [
            vec![
                field(sdf::COORDINATE_FIELD, "1.25"),
                field(sdf::COORDINATE_UNIT_FIELD, "A"),
            ],
            vec![
                field(sdf::COORDINATE_FIELD, "1.5"),
                field(sdf::COORDINATE_UNIT_FIELD, "A"),
            ],
        ]
    );

    // A caller's fields come first, such as a geometry's energy, but cannot
    // give a geometry another place.
    let energy = [SdfDataField::new("energy", "-3.5")];
    let geometry = sdf::SdfRecordSource::from(history.get(1).unwrap()).with_data_fields(&energy);
    let text = sdf::write([geometry], SdfWriteOptions::default())
        .unwrap()
        .0;
    assert_eq!(text, record(&[("energy", "-3.5"), (sdf::STEP_FIELD, "4")]));
    assert_eq!(
        data_fields(&interpret_records(&text)),
        [vec![field("energy", "-3.5"), field(sdf::STEP_FIELD, "4")]]
    );
    for name in [
        sdf::STEP_FIELD,
        sdf::COORDINATE_FIELD,
        sdf::COORDINATE_UNIT_FIELD,
    ] {
        let fields = [SdfDataField::new(name, "1")];
        let geometry =
            sdf::SdfRecordSource::from(history.get(0).unwrap()).with_data_fields(&fields);
        assert_eq!(
            described(sdf::write([geometry], SdfWriteOptions::default()).unwrap_err()),
            (
                DiagnosticKind::Inconsistent,
                Some(FormatErrorDetail::DataField {
                    name: name.to_owned()
                }),
                format!("The data field {name:?} is reserved for the data of the geometry."),
            )
        );
    }
}

#[test]
fn single_structure_writers_fail_on_a_geometry_place_unless_told_to_omit_it() {
    let (history, scan) = methane_series();
    let omit = UnrepresentablePolicy::Omit;
    for (set, field, name) in [
        (&history, GeometryField::OrderStep, "order step"),
        (&scan, GeometryField::OrderCoordinate, "order coordinate"),
    ] {
        let geometry = set.get(1).unwrap();
        let view = geometry.as_structure_view();
        let plain = geometry.to_structure();
        let state = Some(FormatErrorDetail::UnrepresentableState(field));

        let error = molfile::write(
            molfile::MolfileSource::structure(&view),
            MolfileWriteOptions::default(),
        )
        .unwrap_err();
        assert_eq!(
            (error.kind(), error.detail().cloned()),
            (DiagnosticKind::Unsupported, state.clone())
        );
        assert_eq!(
            error.to_string(),
            format!(
                "Cannot write the Molfile output. The {name} of the geometry has no Molfile representation."
            )
        );
        // Omitting the place writes exactly the bare structure.
        assert_eq!(
            molfile::write(
                molfile::MolfileSource::structure(&view),
                MolfileWriteOptions::default().with_geometry_metadata(omit),
            )
            .unwrap()
            .0,
            molfile::write(&plain, MolfileWriteOptions::default())
                .unwrap()
                .0
        );

        let error = mmcif::write(
            [mmcif::MmcifBlockSource::structure(&view)],
            MmcifWriteOptions::default(),
        )
        .unwrap_err();
        assert_eq!(
            (error.kind(), error.detail().cloned()),
            (DiagnosticKind::Unsupported, state)
        );
        assert_eq!(
            error.to_string(),
            format!(
                "Cannot write the mmCIF output. The {name} of the geometry has no mmCIF representation."
            )
        );
        let (text, report) = mmcif::write(
            [mmcif::MmcifBlockSource::structure(&view)],
            MmcifWriteOptions::default().with_geometry_metadata(omit),
        )
        .unwrap();
        assert_eq!(
            report.omissions(),
            [WriteOmission::GeometryState {
                field,
                geometries: 1
            }]
        );
        let (bare, _) = mmcif::write(
            [mmcif::MmcifBlockSource::structure(&plain)],
            MmcifWriteOptions::default(),
        )
        .unwrap();
        assert_eq!(text, bare);
    }
}

/// One methane carbon in a 2 x 3 x 4 nm box periodic along `periodic`.
fn boxed_structure(periodic: [bool; 3]) -> Structure {
    let mut structure = structure("C", &[[1.0, 2.0, 3.0]]);
    structure.geometry_mut().set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(2.0, 3.0, 4.0), NANOMETER),
            periodic,
        )
        .unwrap(),
    ));
    structure
}

#[test]
fn molfile_and_sdf_fail_on_a_periodic_cell_unless_told_to_omit_it() {
    let boxed = boxed_structure([true; 3]);
    let mut open = boxed.clone();
    open.geometry_mut().set_cell(None);
    let message = "The periodic cell of the structure has no Molfile representation.";
    for version in [
        MolfileWriteVersion::Auto,
        MolfileWriteVersion::V2000,
        MolfileWriteVersion::V3000,
    ] {
        let options = MolfileWriteOptions::default().with_version(version);
        let error = molfile::write(&boxed, options).unwrap_err();
        assert_eq!(
            (error.kind(), error.message()),
            (DiagnosticKind::Unsupported, message)
        );
        let omitted = options.with_cell(UnrepresentablePolicy::Omit);
        assert_eq!(
            molfile::write(&boxed, omitted).unwrap().0,
            molfile::write(&open, options).unwrap().0
        );
    }
    let error = sdf::write([&boxed], SdfWriteOptions::default()).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Cannot write the SDF output in record 1. \
         The periodic cell of the structure has no SDF representation."
    );
    let omitted = SdfWriteOptions::from(
        MolfileWriteOptions::default().with_cell(UnrepresentablePolicy::Omit),
    );
    assert_eq!(
        sdf::write([&boxed], omitted).unwrap().0,
        sdf::write([&open], SdfWriteOptions::default()).unwrap().0
    );
}

#[test]
fn mmcif_writes_a_whole_periodic_cell_as_p1_and_reads_it_back() {
    let cell_of = |text: &str| {
        mmcif::parse_str(text, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap()
            .structure()
            .cell()
            .copied()
    };
    let boxed = boxed_structure([true; 3]);
    let text = mmcif::write([&boxed], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    assert!(text.starts_with(
        "data_model\n#\n_cell.length_a 20.000\n_cell.length_b 30.000\n_cell.length_c 40.000\n\
         _cell.angle_alpha 90.000\n_cell.angle_beta 90.000\n_cell.angle_gamma 90.000\n#\n\
         _symmetry.space_group_name_H-M 'P 1'\n#\n"
    ));
    let reread = cell_of(&text).unwrap();
    assert!(reread.is_orthorhombic());
    let lengths = |cell: PeriodicCell| {
        let [a, b, c] = cell.vectors().into_unit(NANOMETER).unwrap().into_value();
        [a.x, b.y, c.z]
    };
    for (actual, expected) in lengths(reread).into_iter().zip([2.0, 3.0, 4.0]) {
        assert!((actual - expected).abs() < 1e-12);
    }

    let slab = boxed_structure([true, true, false]);
    assert_eq!(
        mmcif::write([&slab], MmcifWriteOptions::default())
            .map(|(text, _)| text)
            .unwrap_err()
            .to_string(),
        "Cannot write the mmCIF output. The cell is not periodic along all three axes."
    );
    let omit = MmcifWriteOptions::default().with_cell(UnrepresentablePolicy::Omit);
    assert_eq!(
        cell_of(
            &mmcif::write([&slab], omit.clone())
                .map(|(text, _)| text)
                .unwrap()
        ),
        None
    );

    let mut set = StructureSet::from_structures([boxed.clone(), boxed]).unwrap();
    let text = mmcif::write([&set], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    // One block states the shared cell once, before both models.
    assert_eq!(text.matches("data_").count(), 1, "{text}");
    assert_eq!(text.matches("_cell.length_a 20.000\n").count(), 1, "{text}");
    assert!(
        text.contains("_symmetry.space_group_name_H-M 'P 1'\n"),
        "{text}"
    );
    assert_eq!(
        text.lines()
            .filter(|line| line.starts_with("HETATM"))
            .collect::<Vec<_>>(),
        [
            "HETATM 1 C C1 . MOL M1 1 . . 1.000 2.000 3.000 . . 0 . MOL M1 C1 1",
            "HETATM 2 C C1 . MOL M1 1 . . 1.000 2.000 3.000 . . 0 . MOL M1 C1 2",
        ]
    );
    let members = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap()
        .structure_set()
        .iter()
        .map(|member| member.as_structure_view().cell().is_some())
        .collect::<Vec<_>>();
    assert_eq!(members, [true, true]);
    set.get_mut(1).unwrap().set_cell(None);
    assert_eq!(
        mmcif::write([&set], MmcifWriteOptions::default())
            .map(|(text, _)| text)
            .unwrap_err()
            .to_string(),
        "Cannot write the mmCIF output. The geometries have different cells."
    );
    let pdb_options =
        PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    assert_eq!(
        pdb::write(&set, pdb_options).unwrap_err().to_string(),
        "Cannot write the PDB output. The geometries have different cells."
    );
    let (text, report) = mmcif::write([&set], omit).unwrap();
    assert!(!text.contains("_cell."), "{text}");
    assert!(!text.contains("_symmetry."), "{text}");
    assert_eq!(report.omissions(), [WriteOmission::Cell { structures: 1 }]);
    assert_eq!(
        text.lines()
            .filter(|line| line.starts_with("HETATM"))
            .collect::<Vec<_>>(),
        [
            "HETATM 1 C C1 . MOL M1 1 . . 1.000 2.000 3.000 . . 0 . MOL M1 C1 1",
            "HETATM 2 C C1 . MOL M1 1 . . 1.000 2.000 3.000 . . 0 . MOL M1 C1 2",
        ]
    );
}

#[test]
fn sdf_rejects_titles_and_data_fields_that_would_break_record_framing() {
    let methane = structure("C", &[[0.0, 0.0, 0.0]]);
    let title = |message: &str| {
        (
            DiagnosticKind::Unsupported,
            Some(FormatErrorDetail::Title),
            message.to_owned(),
        )
    };
    let field = |name: &str, message: &str| {
        (
            DiagnosticKind::Unsupported,
            Some(FormatErrorDetail::DataField {
                name: name.to_owned(),
            }),
            message.to_owned(),
        )
    };
    let name = |name: &str| {
        format!(
            "The data field name {name:?} is empty, padded, or holds angle brackets or line breaks."
        )
    };
    for (record_title, fields, expected) in [
        (
            "two\nlines",
            vec![],
            title("The title contains a line break."),
        ),
        (" $$$$ ", vec![], title("The title is \"$$$$\".")),
        (
            "",
            vec![SdfDataField::new("", "value")],
            field("", &name("")),
        ),
        (
            "",
            vec![SdfDataField::new(" ID", "value")],
            field(" ID", &name(" ID")),
        ),
        (
            "",
            vec![SdfDataField::new("A<B", "value")],
            field("A<B", &name("A<B")),
        ),
        (
            "",
            vec![SdfDataField::new("ID", "one\r\ntwo")],
            field(
                "ID",
                "The value of data field \"ID\" contains a carriage return.",
            ),
        ),
        (
            "",
            vec![SdfDataField::new("ID", "one\n\ntwo")],
            field(
                "ID",
                "The value of data field \"ID\" contains a blank line.",
            ),
        ),
        (
            "",
            vec![SdfDataField::new("ID", "one\n$$$$")],
            field(
                "ID",
                "The value of data field \"ID\" contains a \"$$$$\" line.",
            ),
        ),
    ] {
        let record = sdf::SdfRecordSource::structure(&methane)
            .with_title(record_title)
            .with_data_fields(&fields);
        let error = sdf::write([record], SdfWriteOptions::default()).unwrap_err();
        assert_eq!(error.record(), Some(1));
        assert_eq!(described(error), expected);
    }
}

#[test]
fn sdf_explicit_records_preserve_title_and_data_fields() {
    let record = SdfRecordInterpretation::new(
        "named record",
        structure("CO.O", &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
        vec![SdfDataField::new("SOURCE", "round trip")],
    );
    let text = sdf::write(&[record], SdfWriteOptions::default()).unwrap().0;
    let records = interpret_records(&text);
    assert_eq!(records[0].title(), "named record");
    assert_eq!(records[0].data_fields()[0].name(), "SOURCE");
    assert_eq!(records[0].data_fields()[0].value(), "round trip");
    assert_eq!(records[0].structure().topology().instance_count(), 2);
}

#[test]
fn mmcif_blocks_take_their_own_names_and_read_back_under_them() {
    let first = structure("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let second = structure("N", &[[3.0, 0.0, 0.0]]);
    let third = structure("O", &[[4.0, 0.0, 0.0]]);
    let source = |structure, name: Option<&'static str>| {
        let block = mmcif::MmcifBlockSource::structure(structure);
        match name {
            Some(name) => block.with_block_name(name).unwrap(),
            None => block,
        }
    };

    // Named blocks keep their names; an unnamed block keeps its position's
    // number.
    let (text, report) = mmcif::write(
        [
            source(&first, Some("ligand")),
            source(&second, None),
            source(&third, Some("Water.1")),
        ],
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert!(report.is_empty());
    assert_eq!(source(&first, Some("ligand")).block_name(), Some("ligand"));
    assert_eq!(source(&first, None).block_name(), None);
    let document = mmcif::parse_str(&text, Default::default()).unwrap();
    let names = document
        .blocks()
        .iter()
        .map(|block| block.name())
        .collect::<Vec<_>>();
    assert_eq!(names, ["ligand", "model_2", "Water.1"]);
    let atoms = document
        .blocks()
        .iter()
        .map(|block| {
            block
                .interpret(Default::default())
                .unwrap()
                .structure()
                .topology()
                .atom_count()
        })
        .collect::<Vec<_>>();
    assert_eq!(atoms, [2, 1, 1]);
    assert!(document.block("water.1").is_some());

    // A lone named block takes its name rather than the options'.
    let (lone, _) = mmcif::write(
        [source(&first, Some("ligand"))],
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert!(lone.starts_with(
        "data_ligand
"
    ));

    // Names must be valid block names.
    let invalid = mmcif::MmcifBlockSource::structure(&first)
        .with_block_name("two words")
        .unwrap_err();
    assert_eq!(
        invalid.detail(),
        Some(&FormatErrorDetail::MmcifWrite(
            MmcifWriteErrorKind::InvalidBlockName("two words".to_owned())
        ))
    );

    // Names clash case-insensitively, including with generated names, and a
    // clash fails before anything is written.
    for (blocks, clash) in [
        (
            [
                source(&first, Some("Ligand")),
                source(&second, Some("LIGAND")),
            ],
            "LIGAND",
        ),
        (
            [source(&first, None), source(&second, Some("MODEL_1"))],
            "MODEL_1",
        ),
    ] {
        let mut sink = Vec::new();
        let error = mmcif::write_to(&mut sink, blocks, MmcifWriteOptions::default()).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::InvalidValue);
        assert_eq!(
            error.detail(),
            Some(&FormatErrorDetail::MmcifWrite(
                MmcifWriteErrorKind::DuplicateBlockName(clash.to_owned())
            ))
        );
        assert!(sink.is_empty());
    }
}

#[test]
fn mmcif_distinguishes_independent_structures_from_a_structure_set() {
    let first = structure("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let unrelated = structure("N", &[[3.0, 0.0, 0.0]]);
    let independent = mmcif::write(
        [first.clone(), unrelated.clone()]
            .iter()
            .zip(&[classifications(&first), classifications(&unrelated)])
            .map(|(structure, classifications)| {
                mmcif::MmcifBlockSource::structure(structure).with_classifications(classifications)
            }),
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    let independent_document = mmcif::parse_str(&independent, Default::default()).unwrap();
    assert_eq!(independent_document.blocks().len(), 2);
    assert!(independent.starts_with("data_model_1\n"));
    assert!(independent.contains("data_model_2\n"));
    assert!(independent_document
        .interpret_structure_set(Default::default())
        .is_err());

    let shared = first.shared_topology();
    let second = Structure::new(
        Arc::clone(&shared),
        Positions::new(Quantity::new(
            [Point3::new(10.0, 0.0, 0.0), Point3::new(20.0, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    let set = StructureSet::from_structures([first.clone(), second]).unwrap();
    let set_text = mmcif::write(
        [mmcif::MmcifBlockSource::structure_set(&set)
            .with_classifications(&classifications(&first))],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    // One block: both geometries' atoms under their model numbers, and the
    // shared C-O bond once.
    assert_eq!(set_text.matches("data_").count(), 1, "{set_text}");
    assert_eq!(
        set_text
            .lines()
            .filter(|line| line.starts_with("HETATM"))
            .collect::<Vec<_>>(),
        [
            "HETATM 1 C C1 . MOL M1 1 . . 1.000 0.000 0.000 . . 0 . MOL M1 C1 1",
            "HETATM 2 O O2 . MOL M1 1 . . 2.000 0.000 0.000 . . 0 . MOL M1 O2 1",
            "HETATM 3 C C1 . MOL M1 1 . . 10.000 0.000 0.000 . . 0 . MOL M1 C1 2",
            "HETATM 4 O O2 . MOL M1 1 . . 20.000 0.000 0.000 . . 0 . MOL M1 O2 2",
        ]
    );
    assert_eq!(
        set_text
            .lines()
            .filter(|line| line.starts_with("conn"))
            .collect::<Vec<_>>(),
        ["conn1 covale M1 MOL . C1 M1 MOL . O2 . . . . sing"]
    );
    let set_document = mmcif::parse_str(&set_text, Default::default()).unwrap();
    assert_eq!(set_document.blocks().len(), 1);
    let parsed = set_document
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(parsed.reports().len(), 2);
    assert_eq!(parsed.structure_set().len(), 2);
    assert_eq!(parsed.structure_set().order(), None);
    let topology = parsed.structure_set().topology();
    assert_eq!(topology.instance_count(), 1);
    assert_eq!(
        topology
            .atoms()
            .map(|atom| atom.element().symbol())
            .collect::<Vec<_>>(),
        ["C", "O"]
    );
    assert!(
        (parsed
            .structure_set()
            .get(0)
            .unwrap()
            .positions()
            .get(0)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            .x
            - 1.0)
            .abs()
            <= 1e-9
    );
    assert!(
        (parsed
            .structure_set()
            .get(1)
            .unwrap()
            .positions()
            .get(0)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            .x
            - 10.0)
            .abs()
            <= 1e-9
    );

    let round_trip = mmcif::write([&parsed], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    assert_eq!(round_trip, set_text);
    let reparsed = (mmcif::parse_str(&round_trip, Default::default()).unwrap())
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(reparsed.structure_set().len(), 2);
}

#[test]
fn mmcif_structure_is_one_block_with_coordinate_model_one_and_automatic_classification() {
    let structure = structure("C", &[[4.0, 5.0, 6.0]]);
    let automatic = mmcif::write([&structure], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    assert!(automatic.contains("1 non-polymer"));
    let text = mmcif::write(
        [mmcif::MmcifBlockSource::structure(&structure)
            .with_classifications(&classifications(&structure))],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    assert_eq!(
        mmcif::parse_str(&text, Default::default())
            .unwrap()
            .blocks()
            .len(),
        1
    );
    assert!(text.lines().any(|line| line.ends_with(" 1")));
}

#[test]
fn mmcif_export_writes_dense_order_with_reused_atoms_and_unassigned_instances() {
    use kekule::molecule::{Atom, Element, MoleculeEditor};
    use kekule::topology::{AtomSiteMetadata, InstanceAtomId};

    let water = molecule("O");
    let oxygen = water.atom_ids().next().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(water.clone()).unwrap();
    let waters = (0..3)
        .map(|_| builder.add_instance(definition).unwrap())
        .collect::<Vec<_>>();
    let sodium = Atom::new(Element::from_symbol("Na").unwrap()).with_formal_charge(1);
    let mut ion = MoleculeEditor::new();
    ion.add_atom(sodium).unwrap();
    builder.add_molecule(ion.finish().unwrap().clone()).unwrap();
    let hierarchy = builder.hierarchy_mut();
    let chain = hierarchy.add_chain("W", Some("solvent".into())).unwrap();
    // Repeated local AtomIds resolve through their instances; hierarchy order
    // intentionally differs from topology dense order.
    for index in [2, 0, 1] {
        let residue = hierarchy
            .add_residue(chain, "HOH", None, Some((index + 1).to_string()), None)
            .unwrap();
        hierarchy
            .add_atom_site(
                residue,
                InstanceAtomId::new(waters[index], oxygen),
                AtomSiteMetadata::default().with_label_atom_id(Some("O".to_owned())),
            )
            .unwrap();
    }
    let positions = Positions::new(Quantity::new(
        (0..4)
            .map(|i| Point3::new(i as f64 + 1.0, 0.0, 0.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap();
    let mut structure = Structure::new(builder.build().unwrap(), positions).unwrap();
    for index in 0..3 {
        let row = structure
            .topology()
            .atom_index(structure.topology().atom_ids()[index])
            .unwrap();
        structure
            .geometry_mut()
            .set_occupancy(row, Some((index + 1) as f64 / 10.0))
            .unwrap();
    }
    let text = mmcif::write([&structure], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    let document = mmcif::parse_str(&text, Default::default()).unwrap();
    let rows = document.blocks()[0].loop_with_tag("_atom_site.id").unwrap();
    assert_eq!(rows.row_count(), 4);
    // Rows follow dense atom order, not hierarchy order, so coordinate files
    // written in this topology's order stay aligned with the block.
    for (row, index) in [0, 1, 2].into_iter().enumerate() {
        assert_eq!(
            rows.value(row, "_atom_site.label_atom_id").unwrap().text(),
            "O"
        );
        assert_eq!(
            rows.value(row, "_atom_site.auth_seq_id").unwrap().text(),
            (index + 1).to_string()
        );
        assert_eq!(
            rows.value(row, "_atom_site.Cartn_x")
                .unwrap()
                .text()
                .parse::<f64>()
                .unwrap(),
            index as f64 + 1.0
        );
        assert_eq!(
            rows.value(row, "_atom_site.occupancy")
                .unwrap()
                .text()
                .parse::<f64>()
                .unwrap(),
            (index + 1) as f64 / 10.0
        );
    }
    assert_eq!(
        rows.value(3, "_atom_site.type_symbol").unwrap().text(),
        "Na"
    );
    assert_eq!(
        rows.value(3, "_atom_site.Cartn_x")
            .unwrap()
            .text()
            .parse::<f64>()
            .unwrap(),
        4.0
    );
    let reread = document.interpret(Default::default()).unwrap();
    assert_eq!(
        reread
            .structure()
            .positions()
            .values()
            .value()
            .iter()
            .map(|point| (point.x * 10.0).round())
            .collect::<Vec<_>>(),
        [1.0, 2.0, 3.0, 4.0]
    );
}

/// Chain A, a Gly-Gly polymer, and chain B, an ethylene glycol esterified
/// by the second glycine's carboxyl, both in author chain P.
const LINKED_DIPEPTIDE: &str = r#"data_dipeptide
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
covale1 covale A GLY 2 C B EDO . O1
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
ATOM 1 N N GLY A 1 1 0.000 0.000 0.000 1 P
ATOM 2 C CA GLY A 1 1 1.458 0.000 0.000 1 P
ATOM 3 C C GLY A 1 1 2.009 1.420 0.000 1 P
ATOM 4 O O GLY A 1 1 1.251 2.390 0.000 1 P
ATOM 5 N N GLY A 1 2 3.332 1.536 0.000 2 P
ATOM 6 C CA GLY A 1 2 3.988 2.839 0.000 2 P
ATOM 7 C C GLY A 1 2 5.504 2.693 0.000 2 P
ATOM 8 O O GLY A 1 2 6.030 1.580 0.000 2 P
HETATM 9 O O1 EDO B 2 . 6.275 3.789 0.000 101 P
HETATM 10 C C1 EDO B 2 . 7.725 3.789 0.000 101 P
HETATM 11 C C2 EDO B 2 . 8.485 5.105 0.000 101 P
HETATM 12 O O2 EDO B 2 . 9.915 5.105 0.000 101 P
"#;

#[test]
fn mmcif_writes_appended_copies_as_chains_of_their_own_that_keep_author_names() {
    let source = mmcif::parse_str(LINKED_DIPEPTIDE, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .structure()
        .clone();
    // Three bonds in each glycine and in the glycol, the peptide bond, and
    // the ester link.
    assert_eq!(source.topology().bond_count(), 11);
    // Three copies, 20 Å apart, as SeekDelta merges copies of one structure
    // or a crystal's symmetry mates: chains A and B, three times.
    let mut editor = source.to_editor();
    for copy in 1..3 {
        let appended = editor.append_structure(&source).unwrap();
        for &atom in source.topology().atom_ids() {
            let draft = appended.atom(atom).unwrap();
            let point = source.position(atom).unwrap().value_in(ANGSTROM).unwrap();
            editor
                .set_position(
                    draft,
                    Quantity::new(point + Vector3::new(20.0 * copy as f64, 0.0, 0.0), ANGSTROM),
                )
                .unwrap();
        }
    }
    let merged = editor.finish().unwrap();
    let labels = |structure: &Structure| {
        structure
            .topology()
            .chains()
            .map(|chain| {
                (
                    chain.label_id().to_owned(),
                    chain.author_id().map(str::to_owned),
                )
            })
            .collect::<Vec<_>>()
    };
    let author = Some("P".to_owned());
    assert_eq!(
        labels(&merged)
            .into_iter()
            .map(|(label, _)| label)
            .collect::<Vec<_>>(),
        ["A", "B", "A", "B", "A", "B"]
    );

    let text = mmcif::write([&merged], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    let document = mmcif::parse_str(&text, Default::default()).unwrap();
    let block = &document.blocks()[0];
    let asyms = block.loop_with_tag("_struct_asym.id").unwrap();
    assert_eq!(
        (0..asyms.row_count())
            .map(|row| asyms
                .value(row, "_struct_asym.id")
                .unwrap()
                .text()
                .to_owned())
            .collect::<Vec<_>>(),
        ["A", "B", "A-2", "B-2", "A-3", "B-3"]
    );
    // Each copy's ester link names its own copy's asyms.
    let connections = block.loop_with_tag("_struct_conn.id").unwrap();
    assert_eq!(
        (0..connections.row_count())
            .map(|row| {
                [
                    "_struct_conn.ptnr1_label_asym_id",
                    "_struct_conn.ptnr2_label_asym_id",
                ]
                .map(|tag| connections.value(row, tag).unwrap().text().to_owned())
            })
            .collect::<Vec<_>>(),
        [["A", "B"], ["A-2", "B-2"], ["A-3", "B-3"]]
    );

    // The reader makes the written asyms its chains' labels, under the
    // author chain names the copies share.
    let reread = document.interpret(Default::default()).unwrap();
    let reread = reread.structure();
    assert_eq!(
        labels(reread),
        ["A", "B", "A-2", "B-2", "A-3", "B-3"].map(|label| (label.to_owned(), author.clone()))
    );
    assert_eq!(reread.atom_count(), merged.atom_count());
    assert_eq!(
        reread.topology().bond_count(),
        merged.topology().bond_count()
    );
    assert_eq!(reread.topology().bond_count(), 33);
    for (written, read) in merged
        .positions()
        .values()
        .value()
        .iter()
        .zip(reread.positions().values().value().iter())
    {
        assert!((*written - *read).norm() < 1e-3, "{written:?} {read:?}");
    }
}

#[test]
fn mmcif_export_still_rejects_a_partially_assigned_molecule() {
    use kekule::topology::{AtomSiteMetadata, InstanceAtomId};

    let molecule = molecule("CO");
    let atoms = molecule.atom_ids().collect::<Vec<_>>();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(molecule.clone()).unwrap();
    let hierarchy = builder.hierarchy_mut();
    let chain = hierarchy.add_chain("A", None).unwrap();
    let residue = hierarchy
        .add_residue(chain, "UNL", None, None, None)
        .unwrap();
    // Put only the second atom in hierarchy so testing the first atom alone
    // cannot incorrectly choose the hierarchy-free output path.
    hierarchy
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, atoms[1]),
            AtomSiteMetadata::default(),
        )
        .unwrap();
    let structure = Structure::new(builder.build().unwrap(), Positions::zeros(2)).unwrap();
    assert_eq!(
        mmcif_issue(mmcif::write([&structure], MmcifWriteOptions::default()).map(|(text, _)| text)),
        Some(MmcifWriteErrorKind::MissingAtomSite(InstanceAtomId::new(
            instance, atoms[0]
        )))
    );
}

#[test]
fn mmcif_streaming_matches_string_output_for_structures_and_a_structure_set() {
    let first = structure("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let unrelated = structure("N", &[[3.0, 0.0, 0.0]]);
    let structures = [first.clone(), unrelated.clone()];
    let structure_classifications = [classifications(&first), classifications(&unrelated)];
    let expected_structures = mmcif::write(
        structures
            .iter()
            .zip(&structure_classifications)
            .map(|(structure, classifications)| {
                mmcif::MmcifBlockSource::structure(structure).with_classifications(classifications)
            }),
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    let mut streamed_structures = BoundedChunkWriter {
        bytes: Vec::new(),
        maximum_chunk: 64,
    };
    mmcif::write_to(
        &mut streamed_structures,
        structures
            .iter()
            .zip(&structure_classifications)
            .map(|(structure, classifications)| {
                mmcif::MmcifBlockSource::structure(structure).with_classifications(classifications)
            }),
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(streamed_structures.bytes, expected_structures.as_bytes());
    assert_eq!(
        mmcif::parse_str(
            std::str::from_utf8(&streamed_structures.bytes).unwrap(),
            Default::default()
        )
        .unwrap()
        .blocks()
        .len(),
        2
    );

    let second = Structure::new(
        first.shared_topology(),
        Positions::new(Quantity::new(
            [Point3::new(10.0, 0.0, 0.0), Point3::new(20.0, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    let set = StructureSet::from_structures([first.clone(), second]).unwrap();
    let expected_set = mmcif::write(
        [mmcif::MmcifBlockSource::structure_set(&set)
            .with_classifications(&classifications(&first))],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    assert_eq!(expected_set.matches("data_").count(), 1, "{expected_set}");
    let mut streamed_set = BoundedChunkWriter {
        bytes: Vec::new(),
        maximum_chunk: 64,
    };
    mmcif::write_to(
        &mut streamed_set,
        [mmcif::MmcifBlockSource::structure_set(&set)
            .with_classifications(&classifications(&first))],
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(streamed_set.bytes, expected_set.as_bytes());

    let document = mmcif::parse_str(
        std::str::from_utf8(&streamed_set.bytes).unwrap(),
        Default::default(),
    )
    .unwrap();
    assert_eq!(document.blocks().len(), 1);
    assert_eq!(
        document
            .interpret_structure_set(Default::default())
            .unwrap()
            .structure_set()
            .len(),
        2
    );
}

#[test]
fn mmcif_report_count_mismatch_is_context_neutral() {
    let first = structure("C", &[[1.0, 0.0, 0.0]]);
    let second = Structure::new(
        first.shared_topology(),
        Positions::new(Quantity::new([Point3::new(2.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();

    // Each structure block carries its own reports, so a mismatch is reported
    // per block in the same shape as a structure set block's.
    let structure_error = mmcif::write(
        [mmcif::MmcifBlockSource::structure(&first).with_reports(&[])],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap_err();
    assert_eq!(
        structure_error.detail(),
        Some(&FormatErrorDetail::MmcifWrite(
            MmcifWriteErrorKind::ReportCountMismatch(LengthMismatch::new(Countable::Report, 1, 0))
        ))
    );
    assert_eq!(
        structure_error.to_string(),
        "Cannot write the mmCIF output. There are 0 interpretation reports instead of 1."
    );

    let set = StructureSet::from_structures([first, second]).unwrap();
    let set_error = mmcif::write(
        [mmcif::MmcifBlockSource::structure_set(&set).with_reports(&[])],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap_err();
    assert_eq!(
        set_error.detail(),
        Some(&FormatErrorDetail::MmcifWrite(
            MmcifWriteErrorKind::ReportCountMismatch(LengthMismatch::new(Countable::Report, 2, 0))
        ))
    );
    assert_eq!(
        set_error.to_string(),
        "Cannot write the mmCIF output. There are 0 interpretation reports instead of 2."
    );
}

/// The kind, typed detail and message of a format error.
fn described(error: FormatError) -> (DiagnosticKind, Option<FormatErrorDetail>, String) {
    (
        error.kind(),
        error.detail().cloned(),
        error.message().to_owned(),
    )
}

/// The mmCIF writer issue of a failed write.
fn mmcif_issue<T>(result: Result<T, FormatError>) -> Option<MmcifWriteErrorKind> {
    match result.err()?.detail()? {
        FormatErrorDetail::MmcifWrite(issue) => Some(issue.clone()),
        _ => None,
    }
}

struct FailingWriter;

impl Write for FailingWriter {
    fn write(&mut self, _buffer: &[u8]) -> io::Result<usize> {
        Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// The format, stage, kind and sink error kind of a failed write.
fn sink_failure(
    error: FormatError,
) -> (Format, FormatStage, DiagnosticKind, Option<io::ErrorKind>) {
    let sink = std::error::Error::source(&error)
        .and_then(|source| source.downcast_ref::<io::Error>())
        .map(io::Error::kind);
    (error.format(), error.stage(), error.kind(), sink)
}

#[test]
fn streaming_writers_report_structured_io_failures() {
    let structure = structure("C", &[[0.0, 0.0, 0.0]]);
    let expected = |format| {
        (
            format,
            FormatStage::Write,
            DiagnosticKind::Io,
            Some(io::ErrorKind::BrokenPipe),
        )
    };
    assert_eq!(
        sink_failure(
            smiles::write_to(
                &mut FailingWriter,
                &molecule("C"),
                SmilesWriteOptions::default()
            )
            .unwrap_err()
        ),
        expected(Format::Smiles)
    );
    assert_eq!(
        sink_failure(
            molfile::write_to(
                &mut FailingWriter,
                &structure,
                MolfileWriteOptions::default()
            )
            .unwrap_err()
        ),
        expected(Format::Molfile)
    );
    assert_eq!(
        sink_failure(
            sdf::write_to(&mut FailingWriter, [&structure], SdfWriteOptions::default())
                .unwrap_err()
        ),
        expected(Format::Sdf)
    );
    assert_eq!(
        sink_failure(
            mmcif::write_to(
                &mut FailingWriter,
                [mmcif::MmcifBlockSource::structure(&structure)
                    .with_classifications(&classifications(&structure))],
                MmcifWriteOptions::default()
            )
            .unwrap_err()
        ),
        expected(Format::Mmcif)
    );
}

#[test]
fn molfile_titles_and_smiles_names_round_trip_and_reject_unreadable_text() {
    let ethanol = molecule("CCO");
    let positioned = structure("CCO", &[[0.0; 3], [1.5, 0.0, 0.0], [2.0, 1.0, 0.0]]);
    for source in [
        molfile::MolfileSource::molecule(&ethanol),
        molfile::MolfileSource::structure(&positioned),
    ] {
        let text = molfile::write(source.with_title("ethanol 1"), Default::default())
            .unwrap()
            .0;
        let document = molfile::parse_str(&text, Default::default()).unwrap();
        assert_eq!(document.header().title(), "ethanol 1");
        assert_eq!(
            described(molfile::write(source.with_title("a\nb"), Default::default()).unwrap_err()),
            (
                DiagnosticKind::Unsupported,
                Some(FormatErrorDetail::Title),
                "The title contains a line break.".to_owned(),
            )
        );
    }

    let named = smiles::SmilesSource::molecule(&ethanol).with_name("ethanol 1");
    let line = smiles::write(named, SmilesWriteOptions::default()).unwrap();
    assert_eq!(line, "CCO ethanol 1");
    let document = smiles::parse_str(&line, Default::default()).unwrap();
    assert_eq!(document.name(), Some("ethanol 1"));
    let topology = smiles::to_topology("CCO.O").unwrap();
    assert_eq!(
        smiles::write(
            smiles::SmilesSource::topology(&topology).with_name("solution"),
            SmilesWriteOptions::default()
        )
        .unwrap(),
        "CCO.O solution"
    );
    for name in ["", " padded", "two\nlines"] {
        assert_eq!(
            described(
                smiles::write(
                    smiles::SmilesSource::molecule(&ethanol).with_name(name),
                    SmilesWriteOptions::default()
                )
                .unwrap_err()
            ),
            (
                DiagnosticKind::Unsupported,
                Some(FormatErrorDetail::Title),
                "The record name is empty, has surrounding spaces or spans lines.".to_owned(),
            ),
            "{name:?}"
        );
    }
}

/// A write's text and omissions.
fn written(
    result: Result<(String, kekule::diagnostics::WriteReport), FormatError>,
) -> Result<(String, Vec<WriteOmission>), FormatError> {
    result.map(|(text, report)| (text, report.omissions().to_vec()))
}

/// A view, the one piece of state it carries, and the data fields SDF
/// writes for that state.
type ViewCase<'a> = (
    StructureView<'a>,
    Option<GeometryField>,
    Vec<(&'static str, String)>,
);

/// The Molfile, SDF, mmCIF, and PDB writers' outcomes for one view with
/// geometry state set to `policy`.
fn view_outcomes(
    view: StructureView<'_>,
    policy: UnrepresentablePolicy,
) -> [Result<(String, Vec<WriteOmission>), FormatError>; 4] {
    [
        written(molfile::write(
            molfile::MolfileSource::structure(&view),
            MolfileWriteOptions::default().with_geometry_metadata(policy),
        )),
        written(sdf::write(
            [sdf::SdfRecordSource::structure(&view)],
            SdfWriteOptions::from(MolfileWriteOptions::default().with_geometry_metadata(policy)),
        )),
        written(mmcif::write(
            [mmcif::MmcifBlockSource::structure(&view)],
            MmcifWriteOptions::default().with_geometry_metadata(policy),
        )),
        written(pdb::write(
            pdb::PdbSource::structure(&view),
            PdbWriteOptions::default().with_geometry_metadata(policy),
        )),
    ]
}

/// SDF data fields as their text in a record.
fn data_field_text(fields: &[(&str, String)]) -> String {
    fields
        .iter()
        .map(|(name, value)| format!(">  <{name}>\n{value}\n\n"))
        .collect()
}

#[test]
fn geometry_and_frame_views_write_their_state_or_fail_on_it() {
    use GeometryField::{
        Forces as ForceField, OrderCoordinate, OrderStep, Step, Time, Velocities as VelocityField,
    };
    let frames = StatefulFrames::new();
    let frame_views = frames.views();
    let (history, scan) = methane_series();
    // Each view, the one piece of state it carries, and the data fields SDF
    // writes for it: SDF alone carries a step or coordinate place.
    let mut cases: Vec<ViewCase<'_>> = vec![
        (frame_views[0].as_structure_view(), None, Vec::new()),
        (frame_views[1].as_structure_view(), Some(Time), Vec::new()),
        (frame_views[2].as_structure_view(), Some(Step), Vec::new()),
        (
            frame_views[3].as_structure_view(),
            Some(VelocityField),
            Vec::new(),
        ),
        (
            frame_views[4].as_structure_view(),
            Some(ForceField),
            Vec::new(),
        ),
    ];
    for (geometry, step) in history.iter().zip([0, 4, 9]) {
        cases.push((
            geometry.as_structure_view(),
            Some(OrderStep),
            vec![(sdf::STEP_FIELD, step.to_string())],
        ));
    }
    for (geometry, coordinate) in scan.iter().zip(["1.25", "1.5"]) {
        cases.push((
            geometry.as_structure_view(),
            Some(OrderCoordinate),
            vec![
                (sdf::COORDINATE_FIELD, coordinate.to_owned()),
                (sdf::COORDINATE_UNIT_FIELD, "A".to_owned()),
            ],
        ));
    }
    let formats = ["Molfile", "SDF", "mmCIF", "PDB"];
    for (index, (view, field, sdf_fields)) in cases.into_iter().enumerate() {
        // What each writer gives for the bare structure, without state.
        let structure = bare_structure(view);
        let bare = view_outcomes(structure.as_structure_view(), UnrepresentablePolicy::Fail).map(
            |outcome| {
                let (text, omissions) = outcome.unwrap();
                assert!(omissions.is_empty(), "view {index}");
                text
            },
        );
        for policy in [UnrepresentablePolicy::Fail, UnrepresentablePolicy::Omit] {
            let outcomes = view_outcomes(view, policy);
            for (writer, outcome) in outcomes.into_iter().enumerate() {
                let format = formats[writer];
                let context = format!("view {index}, {format}, {policy:?}");
                match field {
                    // SDF writes a place's data fields after the structure.
                    Some(_) if !sdf_fields.is_empty() && format == "SDF" => {
                        let (text, omissions) = outcome.unwrap();
                        let expected = bare[writer]
                            .replace("$$$$\n", &(data_field_text(&sdf_fields) + "$$$$\n"));
                        assert_eq!(text, expected, "{context}");
                        assert!(omissions.is_empty(), "{context}");
                    }
                    None => {
                        let (text, omissions) = outcome.unwrap();
                        assert_eq!(text, bare[writer], "{context}");
                        assert!(omissions.is_empty(), "{context}");
                    }
                    Some(field) if policy == UnrepresentablePolicy::Fail => {
                        let error = outcome.unwrap_err();
                        assert_eq!(error.kind(), DiagnosticKind::Unsupported, "{context}");
                        assert_eq!(
                            error.detail(),
                            Some(&FormatErrorDetail::UnrepresentableState(field)),
                            "{context}"
                        );
                        let record = if format == "SDF" { " in record 1" } else { "" };
                        let verb = if matches!(field, VelocityField | ForceField) {
                            "have"
                        } else {
                            "has"
                        };
                        assert_eq!(
                            error.to_string(),
                            format!(
                                "Cannot write the {format} output{record}. The {field} of the \
                                 geometry {verb} no {format} representation."
                            ),
                            "{context}"
                        );
                    }
                    // Omitting the state writes exactly the bare structure
                    // and reports the omission.
                    Some(field) => {
                        let (text, omissions) = outcome.unwrap();
                        assert_eq!(text, bare[writer], "{context}");
                        assert_eq!(
                            omissions,
                            [WriteOmission::GeometryState {
                                field,
                                geometries: 1
                            }],
                            "{context}"
                        );
                    }
                }
            }
        }
    }
}

/// Two methane geometries, both timed and stepped, the second with velocities.
fn timed_set() -> StructureSet {
    let topology = structure("C", &[[0.0; 3]]).shared_topology();
    let mut geometries = vec![Geometry::new(Positions::zeros(1)); 2];
    let velocities = Velocities::from_vec(Quantity::new(
        vec![Vector3::new(0.0, 0.0, 0.0)],
        NANOMETER_PER_PICOSECOND,
    ))
    .unwrap();
    geometries[1].set_velocities(Some(velocities)).unwrap();
    let mut set = StructureSet::from_geometries(topology, geometries).unwrap();
    set.set_order(Some(Order::Time(
        TimeAxis::new(
            Some(Quantity::new(vec![0.0, 1.0], PICOSECOND)),
            Some(vec![0, 10]),
        )
        .unwrap(),
    )))
    .unwrap();
    set
}

/// Two methane geometries at 1 and 2 Å.
fn plain_structures() -> [Structure; 2] {
    let first = structure("C", &[[1.0, 0.0, 0.0]]);
    let second = Structure::new(
        first.shared_topology(),
        Positions::new(Quantity::new([Point3::new(2.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    [first, second]
}

#[test]
fn structure_sets_write_their_order_or_fail_on_it() {
    let plain = StructureSet::from_structures(plain_structures()).unwrap();
    let (history, scan) = methane_series();
    let timed = timed_set();
    // The writers given each set without its order or frame state.
    let unordered = |set: &StructureSet| {
        StructureSet::from_geometries(
            set.shared_topology(),
            set.iter()
                .map(|geometry| Geometry::new(geometry.positions().clone())),
        )
        .unwrap()
    };
    let atom_rows = |text: &str| {
        text.lines()
            .filter(|line| line.starts_with("HETATM"))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    for policy in [UnrepresentablePolicy::Fail, UnrepresentablePolicy::Omit] {
        let pdb_options = PdbWriteOptions::default()
            .with_implicit_hydrogens(UnrepresentablePolicy::Omit)
            .with_geometry_metadata(policy);
        let mmcif_options = MmcifWriteOptions::default().with_geometry_metadata(policy);
        let sdf_options =
            SdfWriteOptions::from(MolfileWriteOptions::default().with_geometry_metadata(policy));

        // An unordered set writes a model, a block's model number, or a
        // record per geometry, with nothing to omit.
        let (text, omissions) = written(pdb::write(&plain, pdb_options)).unwrap();
        assert_eq!(
            text,
            "MODEL        1\n\
             HETATM    1  C1  MOL     1       1.000   0.000   0.000                       C\n\
             ENDMDL\n\
             MODEL        2\n\
             HETATM    1  C1  MOL     1       2.000   0.000   0.000                       C\n\
             ENDMDL\n\
             END\n"
        );
        assert!(omissions.is_empty());
        let (text, omissions) = written(mmcif::write([&plain], mmcif_options.clone())).unwrap();
        assert_eq!(text.matches("data_").count(), 1, "{text}");
        assert_eq!(
            atom_rows(&text),
            [
                "HETATM 1 C C1 . MOL M1 1 . . 1.000 0.000 0.000 . . 0 . MOL M1 C1 1",
                "HETATM 2 C C1 . MOL M1 1 . . 2.000 0.000 0.000 . . 0 . MOL M1 C1 2",
            ]
        );
        assert!(omissions.is_empty());
        let (text, omissions) = written(sdf::write(&plain, sdf_options)).unwrap();
        let records = interpret_records(&text);
        assert_eq!(data_fields(&records), [Vec::new(), Vec::new()]);
        assert_eq!(
            text,
            sdf::write(&plain_structures(), SdfWriteOptions::default())
                .unwrap()
                .0
        );
        assert!(omissions.is_empty());

        // PDB cannot place the geometries of an ordered set; SDF writes the
        // places as data fields.
        for (set, field, fields) in [
            (
                &history,
                GeometryField::OrderStep,
                ["0", "4", "9"]
                    .map(|step| vec![(sdf::STEP_FIELD, step.to_owned())])
                    .to_vec(),
            ),
            (
                &scan,
                GeometryField::OrderCoordinate,
                ["1.25", "1.5"]
                    .map(|coordinate| {
                        vec![
                            (sdf::COORDINATE_FIELD, coordinate.to_owned()),
                            (sdf::COORDINATE_UNIT_FIELD, "A".to_owned()),
                        ]
                    })
                    .to_vec(),
            ),
        ] {
            let outcome = written(pdb::write(set, pdb_options));
            if policy == UnrepresentablePolicy::Fail {
                let error = outcome.unwrap_err();
                assert_eq!(error.kind(), DiagnosticKind::Unsupported);
                assert_eq!(
                    error.detail(),
                    Some(&FormatErrorDetail::PdbWrite(
                        pdb::PdbWriteErrorKind::UnrepresentableGeometry { geometry: 1, field }
                    ))
                );
                assert_eq!(
                    error.to_string(),
                    format!(
                        "Cannot write the PDB output. The {field} of geometry 1 has no PDB \
                         representation."
                    )
                );
            } else {
                let (text, omissions) = outcome.unwrap();
                assert_eq!(text, pdb::write(&unordered(set), pdb_options).unwrap().0);
                assert_eq!(
                    omissions,
                    [WriteOmission::GeometryState {
                        field,
                        geometries: set.len()
                    }]
                );
            }
            let (text, omissions) = written(sdf::write(set, sdf_options)).unwrap();
            let expected = sdf::write(&unordered(set), SdfWriteOptions::default())
                .unwrap()
                .0
                .split_inclusive("$$$$\n")
                .zip(&fields)
                .map(|(record, fields)| {
                    record.replace("$$$$\n", &(data_field_text(fields) + "$$$$\n"))
                })
                .collect::<String>();
            assert_eq!(text, expected, "{field}");
            assert!(omissions.is_empty());
        }

        // SDF cannot hold the times, steps, and velocities of frames.
        let outcome = written(sdf::write(&timed, sdf_options));
        if policy == UnrepresentablePolicy::Fail {
            let error = outcome.unwrap_err();
            assert_eq!(
                error.detail(),
                Some(&FormatErrorDetail::UnrepresentableState(
                    GeometryField::Time
                ))
            );
            assert_eq!(
                error.to_string(),
                "Cannot write the SDF output in record 1. The time of the geometry has no SDF \
                 representation."
            );
        } else {
            let (text, omissions) = outcome.unwrap();
            assert_eq!(
                text,
                sdf::write(&unordered(&timed), SdfWriteOptions::default())
                    .unwrap()
                    .0
            );
            let state = |field, geometries| WriteOmission::GeometryState { field, geometries };
            assert_eq!(
                omissions,
                [
                    state(GeometryField::Time, 2),
                    state(GeometryField::Step, 2),
                    state(GeometryField::Velocities, 1),
                ]
            );
        }
    }
}

#[test]
fn structure_sets_report_or_refuse_the_state_pdb_and_mmcif_cannot_hold() {
    let set = timed_set();
    let pdb_options =
        PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let error = pdb::write(&set, pdb_options).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);
    assert_eq!(
        error.to_string(),
        "Cannot write the PDB output. The time of geometry 1 has no PDB representation."
    );
    let (text, report) = pdb::write(
        &set,
        pdb_options.with_geometry_metadata(UnrepresentablePolicy::Omit),
    )
    .unwrap();
    let state = |field, geometries| WriteOmission::GeometryState { field, geometries };
    assert!(report.omissions().starts_with(&[
        state(GeometryField::Time, 2),
        state(GeometryField::Step, 2),
        state(GeometryField::Velocities, 1),
    ]));
    // The positions are the geometries' positions.
    let structures = set
        .iter()
        .map(|geometry| bare_structure(geometry.as_structure_view()))
        .collect::<Vec<_>>();
    let bare = StructureSet::from_structures(structures).unwrap();
    assert_eq!(text, pdb::write(&bare, pdb_options).unwrap().0);

    // A structure set without geometries writes nothing.
    let empty = StructureSet::new(bare.shared_topology());
    for error in [
        pdb::write(&empty, pdb_options).unwrap_err(),
        mmcif::write([&empty], MmcifWriteOptions::default()).unwrap_err(),
    ] {
        assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
        assert!(
            error
                .to_string()
                .ends_with("The structure set has no geometries."),
            "{error}"
        );
    }
}

#[test]
fn a_geometry_of_an_ensemble_with_a_place_writes_only_its_place_in_sdf() {
    let mut set = StructureSet::from_structures([structure("C", &[[0.0; 3]])]).unwrap();
    set.set_order(Some(Order::Step(vec![7]))).unwrap();
    let state = StateData::new()
        .with_energy(Quantity::new(-1.5, KILOJOULE_PER_MOLE))
        .unwrap();
    let distribution = Distribution::new(Samples::All, Partition::one_state(1))
        .unwrap()
        .with_log_weights(vec![2.5_f64.ln()])
        .unwrap()
        .with_states(vec![state])
        .unwrap()
        .with_temperature(Quantity::new(300.0, KELVIN))
        .unwrap();
    let ensemble = Ensemble::new(set, distribution).unwrap();
    let view = ensemble.set().get(0).unwrap();
    let text = sdf::write(
        [sdf::SdfRecordSource::structure(&view)],
        SdfWriteOptions::default(),
    )
    .unwrap()
    .0;
    // The place's field alone: the distribution belongs to the ensemble.
    let expected = "
  kekule            2D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <kekule.step>
7

$$$$
";
    assert_eq!(text, expected);
}

#[test]
fn structure_sets_report_or_refuse_the_order_pdb_and_mmcif_lose() {
    use kekule::structure::TimeAxis;
    let plain = StructureSet::from_structures(plain_structures()).unwrap();
    let omit = UnrepresentablePolicy::Omit;
    let pdb_options = PdbWriteOptions::default().with_implicit_hydrogens(omit);
    let (pdb_text, pdb_report) = pdb::write(&plain, pdb_options).unwrap();
    let (mmcif_text, mmcif_report) = mmcif::write([&plain], MmcifWriteOptions::default()).unwrap();
    let state = |field, geometries| WriteOmission::GeometryState { field, geometries };

    // A time order without times or steps.
    let mut bare = plain.clone();
    bare.set_order(Some(Order::Time(TimeAxis::new(None, None).unwrap())))
        .unwrap();

    let cases = [(&bare, GeometryField::Order, 1, 2)];
    for (set, field, geometry, structures) in cases {
        // The order belongs to the whole set.
        let lost = |format| format!("The geometry order has no {format} representation.");
        let error = pdb::write(set, pdb_options).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Unsupported, "{error}");
        assert_eq!(
            error.to_string(),
            format!("Cannot write the PDB output. {}", lost("PDB"))
        );
        let error = mmcif::write([set], MmcifWriteOptions::default()).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!("Cannot write the mmCIF output. {}", lost("mmCIF"))
        );
        assert_eq!(
            error.detail(),
            Some(&FormatErrorDetail::MmcifWrite(
                MmcifWriteErrorKind::UnrepresentableGeometry { geometry, field }
            )),
            "{error}"
        );

        let (text, report) = pdb::write(set, pdb_options.with_geometry_metadata(omit)).unwrap();
        assert_eq!(text, pdb_text, "{field}");
        let lost = state(field, structures);
        assert_eq!(
            report.omissions(),
            [lost]
                .iter()
                .chain(pdb_report.omissions())
                .copied()
                .collect::<Vec<_>>(),
            "{field}"
        );
        let (text, report) = mmcif::write(
            [set],
            MmcifWriteOptions::default().with_geometry_metadata(omit),
        )
        .unwrap();
        assert_eq!(text, mmcif_text, "{field}");
        assert_eq!(
            report.omissions(),
            [lost]
                .iter()
                .chain(mmcif_report.omissions())
                .copied()
                .collect::<Vec<_>>(),
            "{field}"
        );
    }
}

#[test]
fn mmcif_structure_sets_check_the_report_count_before_their_state() {
    let set = conformer_set();
    let reports = Vec::new();
    let error = mmcif::write(
        [mmcif::MmcifBlockSource::structure_set(&set).with_reports(&reports)],
        MmcifWriteOptions::default(),
    )
    .unwrap_err();
    assert!(matches!(
        mmcif_issue::<()>(Err(error)),
        Some(MmcifWriteErrorKind::ReportCountMismatch(payload)) if payload == LengthMismatch::new(Countable::Report, 2, 0)
    ),);
}

#[test]
fn structure_fields_name_the_state_a_format_loses() {
    use GeometryField::*;
    let names = [
        (Order, "geometry order"),
        (Time, "time"),
        (Step, "MD step"),
        (Velocities, "velocities"),
        (Forces, "forces"),
        (AtomProperties, "per-atom geometry properties"),
        (BondProperties, "per-bond geometry properties"),
        (Properties, "geometry properties"),
        (OrderStep, "order step"),
        (OrderCoordinate, "order coordinate"),
        (SetProperties, "structure set properties"),
        (Occupancies, "occupancies"),
        (BFactors, "B factors"),
    ];
    for (field, name) in names {
        assert_eq!(field.to_string(), name, "{field:?}");
    }
}

#[test]
fn mmcif_geometries_share_the_first_geometry_s_entities() {
    // The PDB and mmCIF readings of one entry describe its entity
    // differently, so their reports cannot share one block.
    let from_pdb = pdb::parse_str(include_str!("../fixtures/pdb/1CRN.pdb"), Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let from_mmcif = mmcif::parse_str(
        include_str!("../fixtures/mmcif/1CRN.cif"),
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap();
    let set =
        StructureSet::from_structures([from_pdb.structure().clone(), from_pdb.structure().clone()])
            .unwrap();
    let reports = [from_pdb.report().clone(), from_mmcif.report().clone()];
    let error = mmcif::write(
        [mmcif::MmcifBlockSource::structure_set(&set).with_reports(&reports)],
        MmcifWriteOptions::default(),
    )
    .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::MmcifWrite(
            MmcifWriteErrorKind::IncompatibleGeometry {
                geometry: 2,
                category: kekule::mmcif::MmcifSharedCategory::Entities,
            }
        ))
    );
    assert_eq!(
        error.to_string(),
        "Cannot write the mmCIF output. Geometry 2 has a different entity representation from the first geometry."
    );
}

#[test]
fn atom_field_messages_agree_with_their_subjects() {
    let message = |source: &str, mmcif: bool| {
        let topology = smiles::to_topology(source).unwrap();
        let atoms = topology.atom_count();
        let structure = Structure::new(topology, Positions::zeros(atoms)).unwrap();
        let error = if mmcif {
            mmcif::write([&structure], MmcifWriteOptions::default()).map(|(text, _)| text)
        } else {
            pdb::write(&structure, PdbWriteOptions::default()).map(|(text, _)| text)
        };
        error.unwrap_err().to_string()
    };
    assert_eq!(
        message("[CH4]", true),
        "Cannot write the mmCIF output. The fixed implicit hydrogen count of atom 0 of instance 0 has no mmCIF representation."
    );
    assert_eq!(
        message("[CH3]", true),
        "Cannot write the mmCIF output. The radical electrons of atom 0 of instance 0 have no mmCIF representation."
    );
    assert_eq!(
        message("[13CH4]", true),
        "Cannot write the mmCIF output. The isotope of atom 0 of instance 0 has no mmCIF representation."
    );
    assert_eq!(
        message("[CH3]", false),
        "Cannot write the PDB output. The radical electrons of atom 0 of instance 0 have no PDB representation."
    );
    assert_eq!(
        message("[CH4]", false),
        "Cannot write the PDB output. The fixed implicit hydrogen count of atom 0 of instance 0 has no PDB representation."
    );
}

#[test]
fn molfile_and_sdf_report_what_they_leave_out() {
    use GeometryField::{Forces as ForceField, Step, Time, Velocities as VelocityField};
    let omit = UnrepresentablePolicy::Omit;
    let state = |field, geometries| WriteOmission::GeometryState { field, geometries };
    let plain = structure("C", &[[0.0; 3]]);
    let (_, report) = sdf::write([&plain, &plain], SdfWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    let (_, report) = molfile::write(&plain, MolfileWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    let (_, report) = molfile::write(&molecule("CCO"), MolfileWriteOptions::default()).unwrap();
    assert!(report.is_empty());

    // The cells of the records written without them.
    let boxed = boxed_structure([true; 3]);
    let open = SdfWriteOptions::from(MolfileWriteOptions::default().with_cell(omit));
    let (_, report) = sdf::write([&boxed, &plain, &boxed], open).unwrap();
    assert_eq!(report.omissions(), [WriteOmission::Cell { structures: 2 }]);
    let mut stream = Vec::new();
    let streamed = sdf::write_to(&mut stream, [&boxed, &plain, &boxed], open).unwrap();
    assert_eq!(streamed, report);
    let (_, report) =
        molfile::write(&boxed, MolfileWriteOptions::default().with_cell(omit)).unwrap();
    assert_eq!(report.omissions(), [WriteOmission::Cell { structures: 1 }]);
    let mut stream = Vec::new();
    let streamed = molfile::write_to(
        &mut stream,
        &boxed,
        MolfileWriteOptions::default().with_cell(omit),
    )
    .unwrap();
    assert_eq!(streamed, report);

    // The times, steps, velocities, and forces of frames.
    let frames = StatefulFrames::new();
    let frames = frames.views();
    let metadata =
        SdfWriteOptions::from(MolfileWriteOptions::default().with_geometry_metadata(omit));
    let (_, report) = sdf::write(frames.clone(), metadata).unwrap();
    assert_eq!(
        report.omissions(),
        [
            state(Time, 1),
            state(Step, 1),
            state(VelocityField, 1),
            state(ForceField, 1),
        ]
    );
    // Counted over a time-ordered set's geometries.
    let (_, report) = sdf::write(&timed_set(), metadata).unwrap();
    assert_eq!(
        report.omissions(),
        [state(Time, 2), state(Step, 2), state(VelocityField, 1)]
    );
    let molfile_metadata = MolfileWriteOptions::default().with_geometry_metadata(omit);
    for (frame, field) in frames
        .iter()
        .skip(1)
        .zip([Time, Step, VelocityField, ForceField])
    {
        let (_, report) =
            molfile::write(molfile::MolfileSource::structure(frame), molfile_metadata).unwrap();
        assert_eq!(report.omissions(), [state(field, 1)]);
    }
}

/// `count` methane geometries, geometry `i` at `i` Å along x.
fn methane_set(count: usize) -> StructureSet {
    let first = structure("C", &[[0.0, 0.0, 0.0]]);
    let structures = (0..count).map(|index| {
        Structure::new(
            first.shared_topology(),
            Positions::new(Quantity::new(
                [Point3::new(index as f64, 0.0, 0.0)],
                ANGSTROM,
            ))
            .unwrap(),
        )
        .unwrap()
    });
    StructureSet::from_structures(structures).unwrap()
}

/// A two-atom structure with each kind of state a Molfile and an SDF record
/// cannot carry, one at a time, with the field that names it and the
/// omission that reports it.
fn structures_with_unwritable_state() -> Vec<(Structure, GeometryField, WriteOmission)> {
    use kekule::properties::{PropertyKey, PropertyValue};
    use kekule::topology::{TopologyAtomIndex, TopologyBondIndex};
    use kekule::units::SQUARE_NANOMETER;

    let base = structure("CO", &[[0.0; 3], [1.4, 0.0, 0.0]]);
    let key = || PropertyKey::new("note").unwrap();
    let state = |field| WriteOmission::GeometryState {
        field,
        geometries: 1,
    };
    let mut occupied = base.clone();
    occupied
        .geometry_mut()
        .set_occupancies(Some(vec![Some(0.5), None]))
        .unwrap();
    let mut b_factors = base.clone();
    b_factors
        .geometry_mut()
        .set_b_factors(Some(Quantity::new(vec![None, Some(0.2)], SQUARE_NANOMETER)))
        .unwrap();
    let mut atoms = base.clone();
    atoms
        .geometry_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(
            key(),
            TopologyAtomIndex::new(1),
            Some(PropertyValue::Int(1)),
        )
        .unwrap();
    let mut bonds = base.clone();
    bonds
        .geometry_mut()
        .properties_mut()
        .bonds_mut()
        .set_value(
            key(),
            TopologyBondIndex::new(0),
            Some(PropertyValue::Int(2)),
        )
        .unwrap();
    let mut owner = base.clone();
    owner
        .geometry_mut()
        .properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(3))
        .unwrap();
    vec![
        (
            occupied,
            GeometryField::Occupancies,
            WriteOmission::Occupancies { structures: 1 },
        ),
        (
            b_factors,
            GeometryField::BFactors,
            WriteOmission::BFactors { structures: 1 },
        ),
        (
            atoms,
            GeometryField::AtomProperties,
            state(GeometryField::AtomProperties),
        ),
        (
            bonds,
            GeometryField::BondProperties,
            state(GeometryField::BondProperties),
        ),
        (
            owner,
            GeometryField::Properties,
            state(GeometryField::Properties),
        ),
    ]
}

#[test]
fn molfile_and_sdf_refuse_or_report_sites_and_properties() {
    let omit = UnrepresentablePolicy::Omit;
    let plain = structure("CO", &[[0.0; 3], [1.4, 0.0, 0.0]]);
    let plain_molfile = molfile::write(&plain, MolfileWriteOptions::default())
        .unwrap()
        .0;
    let plain_sdf = sdf::write([&plain], SdfWriteOptions::default()).unwrap().0;
    for (structure, field, omission) in structures_with_unwritable_state() {
        let properties = !matches!(field, GeometryField::Occupancies | GeometryField::BFactors);
        // Each field follows its own policy, and the other one does not
        // omit it.
        let (molfile_omit, molfile_other) = if properties {
            (
                MolfileWriteOptions::default().with_properties(omit),
                MolfileWriteOptions::default().with_geometry_metadata(omit),
            )
        } else {
            (
                MolfileWriteOptions::default().with_geometry_metadata(omit),
                MolfileWriteOptions::default().with_properties(omit),
            )
        };
        let (sdf_omit, sdf_other) = if properties {
            (
                SdfWriteOptions::from(MolfileWriteOptions::default().with_properties(omit)),
                SdfWriteOptions::from(MolfileWriteOptions::default().with_geometry_metadata(omit)),
            )
        } else {
            (
                SdfWriteOptions::from(MolfileWriteOptions::default().with_geometry_metadata(omit)),
                SdfWriteOptions::from(MolfileWriteOptions::default().with_properties(omit)),
            )
        };
        let verb = "have";
        for options in [MolfileWriteOptions::default(), molfile_other] {
            assert_eq!(
                described(molfile::write(&structure, options).unwrap_err()),
                (
                    DiagnosticKind::Unsupported,
                    Some(FormatErrorDetail::UnrepresentableState(field)),
                    format!("The {field} of the structure {verb} no Molfile representation."),
                ),
                "{field:?}"
            );
        }
        for options in [SdfWriteOptions::default(), sdf_other] {
            let error = sdf::write([&plain, &structure], options).unwrap_err();
            assert_eq!(
                error.to_string(),
                format!(
                    "Cannot write the SDF output in record 2. The {field} of the structure {verb} \
                     no SDF representation."
                ),
                "{field:?}"
            );
            assert_eq!(
                error.detail(),
                Some(&FormatErrorDetail::UnrepresentableState(field))
            );
        }
        let (text, report) = molfile::write(&structure, molfile_omit).unwrap();
        assert_eq!(text, plain_molfile, "{field:?}");
        assert_eq!(report.omissions(), [omission], "{field:?}");
        let (text, report) = sdf::write([&structure], sdf_omit).unwrap();
        assert_eq!(text, plain_sdf, "{field:?}");
        assert_eq!(report.omissions(), [omission], "{field:?}");
    }

    // A geometry's sites and properties name the geometry.
    let (structure, _, _) = structures_with_unwritable_state().remove(0);
    let set = StructureSet::from_structures([structure]).unwrap();
    assert_eq!(
        sdf::write(&set, SdfWriteOptions::default())
            .unwrap_err()
            .to_string(),
        "Cannot write the SDF output in record 1. The occupancies of the geometry have no SDF \
         representation."
    );
}

#[test]
fn sdf_and_molfile_refuse_or_report_set_properties() {
    use kekule::properties::{PropertyKey, PropertyValue};
    let omit = UnrepresentablePolicy::Omit;
    let key = || PropertyKey::new("note").unwrap();
    let state = |field, geometries| WriteOmission::GeometryState { field, geometries };

    // The set's own properties, which every record of its geometries loses.
    let mut set = methane_set(2);
    let bare = sdf::write(&set, SdfWriteOptions::default()).unwrap().0;
    set.properties_mut()
        .insert(key(), PropertyValue::Int(1))
        .unwrap();
    let error = sdf::write(&set, SdfWriteOptions::default()).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Cannot write the SDF output in record 1. The structure set properties have no SDF \
         representation."
    );
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::UnrepresentableState(
            GeometryField::SetProperties
        ))
    );
    let (text, report) = sdf::write(
        &set,
        SdfWriteOptions::from(MolfileWriteOptions::default().with_properties(omit)),
    )
    .unwrap();
    assert_eq!(text, bare);
    assert_eq!(report.omissions(), [state(GeometryField::SetProperties, 2)]);
    assert_eq!(
        molfile::write(
            molfile::MolfileSource::structure(&set.get(0).unwrap()),
            MolfileWriteOptions::default()
        )
        .unwrap_err()
        .to_string(),
        "Cannot write the Molfile output. The structure set properties have no Molfile \
         representation."
    );
}

#[test]
fn pdb_and_mmcif_refuse_or_report_structure_and_set_properties() {
    use kekule::properties::{PropertyKey, PropertyValue};
    use kekule::topology::TopologyAtomIndex;
    let omit = UnrepresentablePolicy::Omit;
    let key = || PropertyKey::new("note").unwrap();
    let state = |field, geometries| WriteOmission::GeometryState { field, geometries };

    // A structure's own properties.
    let plain = structure("C", &[[0.0; 3]]);
    let mut noted = plain.clone();
    noted
        .geometry_mut()
        .properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(1))
        .unwrap();
    assert_eq!(
        pdb::write(&noted, PdbWriteOptions::default())
            .unwrap_err()
            .to_string(),
        "Cannot write the PDB output. The geometry properties of the structure have no PDB \
         representation."
    );
    let (text, report) =
        pdb::write(&noted, PdbWriteOptions::default().with_properties(omit)).unwrap();
    assert_eq!(
        text,
        pdb::write(&plain, PdbWriteOptions::default()).unwrap().0
    );
    assert_eq!(report.omissions(), [state(GeometryField::Properties, 1)]);
    let mmcif_default = || MmcifWriteOptions::default();
    assert_eq!(
        mmcif::write([&noted], mmcif_default())
            .unwrap_err()
            .to_string(),
        "Cannot write the mmCIF output. The geometry properties of the structure have no mmCIF \
         representation."
    );
    let (text, report) = mmcif::write([&noted], mmcif_default().with_properties(omit)).unwrap();
    assert_eq!(text, mmcif::write([&plain], mmcif_default()).unwrap().0);
    assert!(report
        .omissions()
        .contains(&state(GeometryField::Properties, 1)));

    // A set's properties and a geometry's.
    let mut set = StructureSet::from_structures([plain.clone(), plain.clone()]).unwrap();
    let bare_pdb = pdb::write(&set, PdbWriteOptions::default()).unwrap().0;
    let bare_mmcif = mmcif::write([&set], mmcif_default()).unwrap().0;
    set.properties_mut()
        .insert(key(), PropertyValue::Int(2))
        .unwrap();
    set.get_mut(1)
        .unwrap()
        .properties_mut()
        .atoms_mut()
        .set_value(
            key(),
            TopologyAtomIndex::new(0),
            Some(PropertyValue::Int(3)),
        )
        .unwrap();
    let error = pdb::write(&set, PdbWriteOptions::default()).unwrap_err();
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::PdbWrite(
            pdb::PdbWriteErrorKind::UnrepresentableGeometry {
                geometry: 1,
                field: GeometryField::SetProperties,
            }
        ))
    );
    assert_eq!(
        error.to_string(),
        "Cannot write the PDB output. The structure set properties have no PDB representation."
    );
    let expected = [
        state(GeometryField::SetProperties, 2),
        state(GeometryField::AtomProperties, 1),
    ];
    let (text, report) =
        pdb::write(&set, PdbWriteOptions::default().with_properties(omit)).unwrap();
    assert_eq!(text, bare_pdb);
    assert_eq!(report.omissions(), expected);
    let (text, report) = mmcif::write([&set], mmcif_default().with_properties(omit)).unwrap();
    assert_eq!(text, bare_mmcif);
    assert!(expected
        .iter()
        .all(|omission| report.omissions().contains(omission)));
    // Without the set's properties, the geometry's name it.
    set.properties_mut().clear();
    assert_eq!(
        mmcif_issue(mmcif::write([&set], mmcif_default())),
        Some(MmcifWriteErrorKind::UnrepresentableGeometry {
            geometry: 2,
            field: GeometryField::AtomProperties,
        })
    );
    // Occupancies are written, as PDB and mmCIF carry them.
    let (occupied, _, _) = structures_with_unwritable_state().remove(0);
    assert!(pdb::write(&occupied, PdbWriteOptions::default()).is_ok());
}
