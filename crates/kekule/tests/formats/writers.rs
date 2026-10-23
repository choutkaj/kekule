use std::sync::Arc;
use std::{io, io::Write};

use kekule::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, FormatStage, StructureField,
    UnrepresentablePolicy,
};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::macromolecular::EntityKind;
use kekule::mmcif::{self, MmcifEntityClassifications, MmcifWriteIssue, MmcifWriteOptions};
use kekule::molfile::{self, MolfileWriteOptions, MolfileWriteVersion};
use kekule::sdf::{self, SdfDataField, SdfRecordInterpretation, SdfWriteOptions};
use kekule::smiles::{self, SmilesWriteMode, SmilesWriteOptions};
use kekule::structure::{
    Ensemble, Forces, Model, ModelView, Positions, Trajectory, TrajectoryFrame, Velocities,
};
use kekule::topology::{Topology, TopologyBuilder};
use kekule::units::{
    Quantity, ANGSTROM, KILOJOULE_PER_MOLE, KILOJOULE_PER_MOLE_PER_NANOMETER, NANOMETER,
    NANOMETER_PER_PICOSECOND, PICOSECOND,
};

fn molecule(source: &str) -> kekule::molecule::Molecule {
    let mut molecules = smiles::to_molecules(source).expect("SMILES interprets");
    assert_eq!(molecules.len(), 1);
    molecules.pop().unwrap()
}

fn model(source: &str, points: &[[f64; 3]]) -> Model {
    let topology = smiles::to_topology(source).expect("topology interprets");
    let positions = Positions::new(Quantity::new(
        points
            .iter()
            .map(|point| Point3::new(point[0], point[1], point[2]))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap();
    Model::new(topology, positions).unwrap()
}

fn classifications(model: &Model) -> MmcifEntityClassifications {
    let mut classifications = MmcifEntityClassifications::new();
    for instance in model.topology().instances().map(|molecule| molecule.id()) {
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
            &model("F[C@H](Cl)Br", &points),
            MolfileWriteOptions::default().with_version(version),
        )
        .unwrap();
        let error = molfile::write(
            &model("F[C@@H](Cl)Br", &points),
            MolfileWriteOptions::default().with_version(version),
        )
        .unwrap_err();
        assert_eq!(
            (error.kind(), error.message()),
            (
                DiagnosticKind::Unsupported,
                "specified tetrahedral stereo at atom 1 contradicts the emitted coordinates, which readers take as the configuration"
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
    let charged = model(
        "[13CH3-]",
        &[[-1.2345678901234567, -2.345678901234568, -3.456789012345679]],
    );
    let text = molfile::write(
        &charged,
        MolfileWriteOptions::default().with_version(MolfileWriteVersion::V3000),
    )
    .unwrap();
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
            &model(source, &layout),
            MolfileWriteOptions::default().with_version(version),
        )
        .unwrap();
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
    let enhanced = model("C[C@H](N)O |&1:1|", &layout);
    let sodium = molecule("[Na+]");
    let large = Model::new(
        Topology::from_molecules(vec![sodium; 1_000]).unwrap(),
        Positions::zeros(1_000),
    )
    .unwrap();
    for model in [&enhanced, &large] {
        let explicit = molfile::write(
            model,
            MolfileWriteOptions::default().with_version(MolfileWriteVersion::V3000),
        )
        .unwrap();
        assert_eq!(
            molfile::write(model, MolfileWriteOptions::default()).unwrap(),
            explicit
        );
    }
}

#[test]
fn molfile_model_flattens_disconnected_topology_and_preserves_coordinates() {
    let model = model(
        "CO.O",
        &[[1.25, 2.5, 3.75], [2.0, 2.5, 3.75], [-4.0, 0.5, 8.0]],
    );
    let text = molfile::write(&model, MolfileWriteOptions::default()).unwrap();
    assert!(text.contains("V2000"));

    let interpreted = (molfile::parse_str(&text, Default::default()).unwrap())
        .interpret(Default::default())
        .unwrap();
    let round_trip_model = interpreted.into_model();
    assert_eq!(round_trip_model.topology().instance_count(), 2);
    for (index, expected) in [
        Point3::new(1.25, 2.5, 3.75),
        Point3::new(2.0, 2.5, 3.75),
        Point3::new(-4.0, 0.5, 8.0),
    ]
    .into_iter()
    .enumerate()
    {
        let actual = round_trip_model
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
    let model = Model::new(topology, Positions::zeros(1_000)).unwrap();

    let error = molfile::write(
        &model,
        MolfileWriteOptions::default().with_version(MolfileWriteVersion::V2000),
    )
    .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);

    let promoted = molfile::write(&model, MolfileWriteOptions::default()).unwrap();
    assert!(promoted.contains("V3000"));
    assert!(promoted.contains("M  V30 COUNTS 1000 0 0 0 0"));
    assert_eq!(
        (molfile::parse_str(&promoted, Default::default()).unwrap())
            .interpret(Default::default())
            .unwrap()
            .into_model()
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
    let conformer = model("CO", &[[0.0, 0.0, 0.0], [1.2, 0.3, -0.4]]);
    let drawing = model("CO", &[[0.0, 0.0, 0.0], [1.2, 0.3, 0.0]]);
    // An unmarked flat record keeps its plane, so its z is emitted as given.
    let raised_plane = model("CO", &[[0.0, 0.0, 1.5], [1.2, 0.3, 1.5]]);
    // Both z values print as 0.0000 in V2000 but stay distinct in V3000.
    let sub_precision_relief = model("CO", &[[0.0, 0.0, 0.0], [1.2, 0.0, 0.00004]]);
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
        let text =
            molfile::write(source, MolfileWriteOptions::default().with_version(version)).unwrap();
        assert_eq!(text.lines().nth(1), Some(expected), "{version:?}\n{text}");
        assert_eq!(text.lines().nth(2), Some(""));
    }

    let text = sdf::write(
        [
            sdf::SdfRecordSource::model(&conformer).with_title("conformer"),
            sdf::SdfRecordSource::model(&drawing).with_title("drawing"),
        ],
        SdfWriteOptions::default(),
    )
    .unwrap();
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
    let raised_drawing = model("C[C@H](N)O", &layout);
    let canonical =
        |topology: &Topology| smiles::write(topology, SmilesWriteOptions::canonical()).unwrap();
    assert_eq!(canonical(raised_drawing.topology()), "C[C@H](N)O");
    for version in [MolfileWriteVersion::V2000, MolfileWriteVersion::V3000] {
        let text = molfile::write(
            &raised_drawing,
            MolfileWriteOptions::default().with_version(version),
        )
        .unwrap();
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
        let written = document.interpret(Default::default()).unwrap().into_model();
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
fn sdf_writes_independent_models_and_ensemble_members_in_order() {
    let first = model("C", &[[1.0, 0.0, 0.0]]);
    let second = model("O", &[[2.0, 0.0, 0.0]]);
    let text = sdf::write(&[first.clone(), second.clone()], SdfWriteOptions::default()).unwrap();
    let interpreted = interpret_records(&text);
    assert_eq!(interpreted.len(), 2);
    assert_eq!(
        interpreted[0]
            .model()
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
            .model()
            .topology()
            .atoms()
            .next()
            .unwrap()
            .element()
            .symbol(),
        "O"
    );

    let shared = first.shared_topology();
    let later = Model::new(
        shared,
        Positions::new(Quantity::new([Point3::new(9.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    let ensemble = Ensemble::from_models([first, later]).unwrap();
    let ensemble_text = sdf::write(&ensemble, SdfWriteOptions::default()).unwrap();
    let records = interpret_records(&ensemble_text);
    assert_eq!(records.len(), 2);
    assert!(
        (records[0]
            .model()
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
            .model()
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

/// A two-member methane ensemble weighted 1:3 with energies 0 and 5 kJ/mol.
fn weighted_ensemble() -> Ensemble {
    let first = model("C", &[[1.0, 0.0, 0.0]]);
    let later = Model::new(
        first.shared_topology(),
        Positions::new(Quantity::new([Point3::new(9.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    let mut ensemble = Ensemble::from_models([first, later]).unwrap();
    for (index, weight, energy) in [(0, 1.0, 0.0), (1, 3.0, 5.0)] {
        let mut member = ensemble.get_mut(index).unwrap();
        member.set_weight(weight).unwrap();
        member
            .set_energy(Some(Quantity::new(energy, KILOJOULE_PER_MOLE)))
            .unwrap();
    }
    ensemble
}

fn member_metadata(ensemble: &Ensemble) -> Vec<(f64, Option<f64>)> {
    ensemble
        .iter()
        .map(|member| {
            (
                member.weight(),
                member
                    .energy()
                    .map(|energy| energy.value_in(KILOJOULE_PER_MOLE).unwrap()),
            )
        })
        .collect()
}

#[test]
fn sdf_ensemble_members_carry_their_weight_and_energy_through_data_fields() {
    let ensemble = weighted_ensemble();
    let text = sdf::write(&ensemble, SdfWriteOptions::default()).unwrap();
    let records = interpret_records(&text);
    let fields = records
        .iter()
        .map(|record| {
            record
                .data_fields()
                .iter()
                .map(|field| (field.name().to_owned(), field.value().to_owned()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let field = |name: &str, value: &str| (name.to_owned(), value.to_owned());
    assert_eq!(
        fields,
        [
            vec![field(sdf::WEIGHT_FIELD, "1"), field(sdf::ENERGY_FIELD, "0")],
            vec![field(sdf::WEIGHT_FIELD, "3"), field(sdf::ENERGY_FIELD, "5")],
        ]
    );
    let reread = sdf::ensemble(&records).unwrap();
    assert_eq!(member_metadata(&reread), member_metadata(&ensemble));
    assert_eq!(reread.probabilities().unwrap(), [0.25, 0.75]);

    // Records without the fields are equally weighted and have no energy.
    let plain = interpret_records(
        &sdf::write(&[model("C", &[[0.0; 3]])], SdfWriteOptions::default()).unwrap(),
    );
    assert_eq!(
        member_metadata(&sdf::ensemble(&plain).unwrap()),
        [(1.0, None)]
    );

    let malformed = text.replacen(">  <kekule.weight>\n3\n", ">  <kekule.weight>\nheavy\n", 1);
    assert_eq!(
        sdf::ensemble(&interpret_records(&malformed))
            .unwrap_err()
            .to_string(),
        "SDF interpretation error in record 2: data field `kekule.weight` is not a number: `heavy`"
    );
    let fields = [SdfDataField::new(sdf::WEIGHT_FIELD, "2")];
    let member = sdf::SdfRecordSource::from(ensemble.get(0).unwrap()).with_data_fields(&fields);
    assert_eq!(
        described(sdf::write([member], SdfWriteOptions::default()).unwrap_err()),
        (
            DiagnosticKind::Inconsistent,
            Some(FormatErrorDetail::DataField {
                name: sdf::WEIGHT_FIELD.to_owned(),
            }),
            "this field carries the ensemble member's own metadata".to_owned(),
        )
    );
}

#[test]
fn mmcif_rejects_ensemble_weights_and_energies_it_cannot_carry() {
    let mut ensemble = weighted_ensemble();
    let error = mmcif::write([&ensemble], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert_eq!(
        error.to_string(),
        "mmCIF write error: mmCIF cannot represent the energy of ensemble member 1"
    );
    for index in 0..2 {
        ensemble.get_mut(index).unwrap().set_energy(None).unwrap();
    }
    let error = mmcif::write([&ensemble], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert_eq!(
        error.to_string(),
        "mmCIF write error: mmCIF cannot represent the weight of ensemble member 2"
    );
    let omitted = mmcif::write(
        [&weighted_ensemble()],
        MmcifWriteOptions::default().with_ensemble_metadata(UnrepresentablePolicy::Omit),
    )
    .map(|(text, _)| text)
    .unwrap();
    let reread = mmcif::parse_str(&omitted, Default::default())
        .unwrap()
        .interpret_ensemble(Default::default())
        .unwrap();
    assert_eq!(
        member_metadata(reread.ensemble()),
        [(1.0, None), (1.0, None)]
    );

    // Equal weights without energies are what interpretation gives back.
    ensemble.get_mut(1).unwrap().set_weight(1.0).unwrap();
    mmcif::write([&ensemble], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
}

/// One-atom methane frames with exactly one piece of frame state each:
/// none, time, step, velocities, forces.
fn stateful_trajectory() -> Trajectory {
    let topology = model("C", &[[0.0; 3]]).shared_topology();
    let mut trajectory = Trajectory::new(topology);
    for _ in 0..5 {
        trajectory
            .push(TrajectoryFrame::new(Positions::zeros(1)))
            .unwrap();
    }
    let zero = || vec![Vector3::new(0.0, 0.0, 0.0)];
    trajectory
        .get_mut(1)
        .unwrap()
        .set_time(Some(Quantity::new(2.0, PICOSECOND)))
        .unwrap();
    trajectory.get_mut(2).unwrap().set_step(Some(7));
    let velocities = Velocities::from_vec(Quantity::new(zero(), NANOMETER_PER_PICOSECOND)).unwrap();
    trajectory
        .get_mut(3)
        .unwrap()
        .set_velocities(Some(velocities))
        .unwrap();
    let forces = Forces::from_vec(Quantity::new(zero(), KILOJOULE_PER_MOLE_PER_NANOMETER)).unwrap();
    trajectory
        .get_mut(4)
        .unwrap()
        .set_forces(Some(forces))
        .unwrap();
    trajectory
}

/// Writes one structure view as a Molfile, an SDF record and an mmCIF model
/// block, each with member and frame state set to `policy`.
fn write_structure(
    view: ModelView<'_>,
    policy: UnrepresentablePolicy,
) -> [Result<String, FormatError>; 3] {
    [
        molfile::write(
            molfile::MolfileSource::model(&view),
            MolfileWriteOptions::default()
                .with_ensemble_metadata(policy)
                .with_frame_metadata(policy),
        ),
        sdf::write(
            [sdf::SdfRecordSource::model(&view)],
            SdfWriteOptions::default().with_frame_metadata(policy),
        ),
        mmcif::write(
            [mmcif::MmcifBlockSource::model(&view)],
            MmcifWriteOptions::default()
                .with_ensemble_metadata(policy)
                .with_frame_metadata(policy),
        )
        .map(|(text, _)| text),
    ]
}

#[test]
fn structure_writers_fail_on_member_and_frame_state_unless_told_to_omit_it() {
    use StructureField::{
        Energy, Forces as ForceField, Step, Time, Velocities as VelocityField, Weight,
    };
    let mut ensemble = weighted_ensemble();
    ensemble.get_mut(1).unwrap().set_energy(None).unwrap();
    let trajectory = stateful_trajectory();
    // Expected detail per writer (Molfile, SDF, mmCIF); SDF writes member
    // weights and energies as data fields.
    let rows: [(ModelView<'_>, [Option<StructureField>; 3]); 7] = [
        (
            ensemble.get(0).unwrap().as_model_view(),
            [Some(Energy), None, Some(Energy)],
        ),
        (
            ensemble.get(1).unwrap().as_model_view(),
            [Some(Weight), None, Some(Weight)],
        ),
        (trajectory.get(0).unwrap().as_model_view(), [None; 3]),
        (trajectory.get(1).unwrap().as_model_view(), [Some(Time); 3]),
        (trajectory.get(2).unwrap().as_model_view(), [Some(Step); 3]),
        (
            trajectory.get(3).unwrap().as_model_view(),
            [Some(VelocityField); 3],
        ),
        (
            trajectory.get(4).unwrap().as_model_view(),
            [Some(ForceField); 3],
        ),
    ];
    let formats = [Format::Molfile, Format::Sdf, Format::Mmcif];
    for (row, (view, expected)) in rows.into_iter().enumerate() {
        let model = view.to_model();
        let plain = write_structure(model.as_model_view(), UnrepresentablePolicy::Fail);
        let failed = write_structure(view, UnrepresentablePolicy::Fail);
        let omitted = write_structure(view, UnrepresentablePolicy::Omit);
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
                    // Omitting the state writes exactly the bare model.
                    assert_eq!(
                        omitted[writer].as_ref().unwrap(),
                        plain[writer].as_ref().unwrap(),
                        "row {row}, {format}"
                    );
                }
                None if format == Format::Sdf && row < 2 => {
                    // Both SDF paths write the member's own data fields.
                    let member = ensemble.get(row).unwrap();
                    let typed = sdf::write([member], SdfWriteOptions::default()).unwrap();
                    assert_eq!(failed[writer].as_ref().unwrap(), &typed);
                    assert!(typed.contains(
                        ">  <kekule.weight>
"
                    ));
                    assert_ne!(&typed, plain[writer].as_ref().unwrap());
                }
                None => assert_eq!(
                    failed[writer].as_ref().unwrap(),
                    plain[writer].as_ref().unwrap(),
                    "row {row}, {format}"
                ),
            }
        }
    }
    let error = sdf::write(&trajectory, SdfWriteOptions::default()).unwrap_err();
    assert_eq!(
        error.to_string(),
        "SDF write error in record 2: SDF cannot represent the trajectory frame time; \
         choose to omit it explicitly"
    );
}

/// One methane carbon in a 2 x 3 x 4 nm box periodic along `periodic`.
fn boxed_model(periodic: [bool; 3]) -> Model {
    let mut model = model("C", &[[1.0, 2.0, 3.0]]);
    model.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(2.0, 3.0, 4.0), NANOMETER),
            periodic,
        )
        .unwrap(),
    ));
    model
}

#[test]
fn molfile_and_sdf_fail_on_a_periodic_cell_unless_told_to_omit_it() {
    let boxed = boxed_model([true; 3]);
    let mut open = boxed.clone();
    open.set_cell(None);
    let message = "a Molfile cannot represent a periodic cell; set the cell option to omit it";
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
            molfile::write(&boxed, omitted).unwrap(),
            molfile::write(&open, options).unwrap()
        );
    }
    let error = sdf::write([&boxed], SdfWriteOptions::default()).unwrap_err();
    assert_eq!(
        error.to_string(),
        format!("SDF write error in record 1: {message}")
    );
    let omitted = SdfWriteOptions::default().with_cell(UnrepresentablePolicy::Omit);
    assert_eq!(
        sdf::write([&boxed], omitted).unwrap(),
        sdf::write([&open], SdfWriteOptions::default()).unwrap()
    );
}

#[test]
fn mmcif_writes_a_whole_periodic_cell_as_p1_and_reads_it_back() {
    let cell_of = |text: &str| {
        mmcif::parse_str(text, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap()
            .model()
            .cell()
            .copied()
    };
    let boxed = boxed_model([true; 3]);
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

    let slab = boxed_model([true, true, false]);
    assert_eq!(
        mmcif::write([&slab], MmcifWriteOptions::default()).map(|(text, _)| text)
            .unwrap_err()
            .to_string(),
        "mmCIF write error: mmCIF cannot represent the periodic cell: an mmCIF cell is periodic along all three axes"
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

    let mut ensemble = Ensemble::from_models([boxed.clone(), boxed]).unwrap();
    let text = mmcif::write([&ensemble], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    let members = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret_ensemble(Default::default())
        .unwrap()
        .ensemble()
        .iter()
        .map(|member| member.as_model_view().cell().is_some())
        .collect::<Vec<_>>();
    assert_eq!(members, [true, true]);
    ensemble.get_mut(1).unwrap().set_cell(None);
    assert_eq!(
        mmcif::write([&ensemble], MmcifWriteOptions::default()).map(|(text, _)| text)
            .unwrap_err()
            .to_string(),
        "mmCIF write error: mmCIF cannot represent the periodic cell: ensemble members have different cells"
    );
    mmcif::write([&ensemble], omit)
        .map(|(text, _)| text)
        .unwrap();
}

#[test]
fn sdf_rejects_titles_and_data_fields_that_would_break_record_framing() {
    let methane = model("C", &[[0.0, 0.0, 0.0]]);
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
    const NAME: &str =
        "SDF data field names must be nonempty, trimmed, and exclude angle brackets or line breaks";
    for (record_title, fields, expected) in [
        (
            "two\nlines",
            vec![],
            title("SDF record titles cannot contain line breaks"),
        ),
        (
            " $$$$ ",
            vec![],
            title("SDF record titles cannot be a record delimiter"),
        ),
        ("", vec![SdfDataField::new("", "value")], field("", NAME)),
        (
            "",
            vec![SdfDataField::new(" ID", "value")],
            field(" ID", NAME),
        ),
        (
            "",
            vec![SdfDataField::new("A<B", "value")],
            field("A<B", NAME),
        ),
        (
            "",
            vec![SdfDataField::new("ID", "one\r\ntwo")],
            field(
                "ID",
                "SDF data field values cannot contain carriage returns",
            ),
        ),
        (
            "",
            vec![SdfDataField::new("ID", "one\n\ntwo")],
            field("ID", "SDF data field values cannot contain blank lines"),
        ),
        (
            "",
            vec![SdfDataField::new("ID", "one\n$$$$")],
            field(
                "ID",
                "SDF data field values cannot contain a record delimiter line",
            ),
        ),
    ] {
        let record = sdf::SdfRecordSource::model(&methane)
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
        model("CO.O", &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
        vec![SdfDataField::new("SOURCE", "round trip")],
    );
    let text = sdf::write(&[record], SdfWriteOptions::default()).unwrap();
    let records = interpret_records(&text);
    assert_eq!(records[0].title(), "named record");
    assert_eq!(records[0].data_fields()[0].name(), "SOURCE");
    assert_eq!(records[0].data_fields()[0].value(), "round trip");
    assert_eq!(records[0].model().topology().instance_count(), 2);
}

#[test]
fn mmcif_blocks_take_their_own_names_and_read_back_under_them() {
    let first = model("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let second = model("N", &[[3.0, 0.0, 0.0]]);
    let third = model("O", &[[4.0, 0.0, 0.0]]);
    let source = |model, name: Option<&'static str>| {
        let block = mmcif::MmcifBlockSource::model(model);
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
                .model()
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
    let invalid = mmcif::MmcifBlockSource::model(&first)
        .with_block_name("two words")
        .unwrap_err();
    assert_eq!(
        invalid.detail(),
        Some(&FormatErrorDetail::MmcifWrite(
            MmcifWriteIssue::InvalidBlockName("two words".to_owned())
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
                MmcifWriteIssue::DuplicateBlockName(clash.to_owned())
            ))
        );
        assert!(sink.is_empty());
    }
}

#[test]
fn mmcif_distinguishes_independent_models_from_an_ensemble() {
    let first = model("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let unrelated = model("N", &[[3.0, 0.0, 0.0]]);
    let independent = mmcif::write(
        [first.clone(), unrelated.clone()]
            .iter()
            .zip(&[classifications(&first), classifications(&unrelated)])
            .map(|(model, classifications)| {
                mmcif::MmcifBlockSource::model(model).with_classifications(classifications)
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
        .interpret_ensemble(Default::default())
        .is_err());

    let shared = first.shared_topology();
    let second = Model::new(
        Arc::clone(&shared),
        Positions::new(Quantity::new(
            [Point3::new(10.0, 0.0, 0.0), Point3::new(20.0, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    let ensemble = Ensemble::from_models([first.clone(), second]).unwrap();
    let ensemble_text = mmcif::write(
        [mmcif::MmcifBlockSource::ensemble(&ensemble)
            .with_classifications(&classifications(&first))],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    let ensemble_document = mmcif::parse_str(&ensemble_text, Default::default()).unwrap();
    assert_eq!(ensemble_document.blocks().len(), 1);
    let parsed = ensemble_document
        .interpret_ensemble(Default::default())
        .unwrap();
    assert_eq!(parsed.ensemble().len(), 2);
    assert!(
        (parsed
            .ensemble()
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
            .ensemble()
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
    let reparsed = (mmcif::parse_str(&round_trip, Default::default()).unwrap())
        .interpret_ensemble(Default::default())
        .unwrap();
    assert_eq!(reparsed.ensemble().len(), 2);
}

#[test]
fn mmcif_model_is_one_block_with_coordinate_model_one_and_automatic_classification() {
    let model = model("C", &[[4.0, 5.0, 6.0]]);
    let automatic = mmcif::write([&model], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    assert!(automatic.contains("1 non-polymer"));
    let text = mmcif::write(
        [mmcif::MmcifBlockSource::model(&model).with_classifications(&classifications(&model))],
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
    let mut model = Model::new(builder.build().unwrap(), positions).unwrap();
    for index in 0..3 {
        let row = model
            .topology()
            .atom_index(model.topology().atom_ids()[index])
            .unwrap();
        model
            .set_occupancy(row, Some((index + 1) as f64 / 10.0))
            .unwrap();
    }
    let text = mmcif::write([&model], MmcifWriteOptions::default())
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
            .model()
            .positions()
            .values()
            .value()
            .iter()
            .map(|point| (point.x * 10.0).round())
            .collect::<Vec<_>>(),
        [1.0, 2.0, 3.0, 4.0]
    );
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
    let model = Model::new(builder.build().unwrap(), Positions::zeros(2)).unwrap();
    assert_eq!(
        mmcif_issue(mmcif::write([&model], MmcifWriteOptions::default()).map(|(text, _)| text)),
        Some(MmcifWriteIssue::MissingAtomSite(InstanceAtomId::new(
            instance, atoms[0]
        )))
    );
}

#[test]
fn mmcif_streaming_matches_string_output_for_models_and_ensemble() {
    let first = model("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let unrelated = model("N", &[[3.0, 0.0, 0.0]]);
    let models = [first.clone(), unrelated.clone()];
    let model_classifications = [classifications(&first), classifications(&unrelated)];
    let expected_models = mmcif::write(
        models
            .iter()
            .zip(&model_classifications)
            .map(|(model, classifications)| {
                mmcif::MmcifBlockSource::model(model).with_classifications(classifications)
            }),
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    let mut streamed_models = BoundedChunkWriter {
        bytes: Vec::new(),
        maximum_chunk: 64,
    };
    mmcif::write_to(
        &mut streamed_models,
        models
            .iter()
            .zip(&model_classifications)
            .map(|(model, classifications)| {
                mmcif::MmcifBlockSource::model(model).with_classifications(classifications)
            }),
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(streamed_models.bytes, expected_models.as_bytes());
    assert_eq!(
        mmcif::parse_str(
            std::str::from_utf8(&streamed_models.bytes).unwrap(),
            Default::default()
        )
        .unwrap()
        .blocks()
        .len(),
        2
    );

    let second = Model::new(
        first.shared_topology(),
        Positions::new(Quantity::new(
            [Point3::new(10.0, 0.0, 0.0), Point3::new(20.0, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    let ensemble = Ensemble::from_models([first.clone(), second]).unwrap();
    let expected_ensemble = mmcif::write(
        [mmcif::MmcifBlockSource::ensemble(&ensemble)
            .with_classifications(&classifications(&first))],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    let mut streamed_ensemble = BoundedChunkWriter {
        bytes: Vec::new(),
        maximum_chunk: 64,
    };
    mmcif::write_to(
        &mut streamed_ensemble,
        [mmcif::MmcifBlockSource::ensemble(&ensemble)
            .with_classifications(&classifications(&first))],
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(streamed_ensemble.bytes, expected_ensemble.as_bytes());

    let document = mmcif::parse_str(
        std::str::from_utf8(&streamed_ensemble.bytes).unwrap(),
        Default::default(),
    )
    .unwrap();
    assert_eq!(document.blocks().len(), 1);
    assert_eq!(
        document
            .interpret_ensemble(Default::default())
            .unwrap()
            .ensemble()
            .len(),
        2
    );
}

#[test]
fn mmcif_report_count_mismatch_is_context_neutral() {
    let first = model("C", &[[1.0, 0.0, 0.0]]);
    let second = Model::new(
        first.shared_topology(),
        Positions::new(Quantity::new([Point3::new(2.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();

    // Each model block carries its own reports, so a mismatch is reported
    // per block in the same shape as an ensemble block's.
    let model_error = mmcif::write(
        [mmcif::MmcifBlockSource::model(&first).with_reports(&[])],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap_err();
    assert_eq!(
        model_error.detail(),
        Some(&FormatErrorDetail::MmcifWrite(
            MmcifWriteIssue::ReportCountMismatch {
                expected: 1,
                actual: 0
            }
        ))
    );
    assert_eq!(
        model_error.to_string(),
        "mmCIF write error: mmCIF writing requires 1 interpretation reports, but received 0"
    );

    let ensemble = Ensemble::from_models([first, second]).unwrap();
    let ensemble_error = mmcif::write(
        [mmcif::MmcifBlockSource::ensemble(&ensemble).with_reports(&[])],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap_err();
    assert_eq!(
        ensemble_error.detail(),
        Some(&FormatErrorDetail::MmcifWrite(
            MmcifWriteIssue::ReportCountMismatch {
                expected: 2,
                actual: 0
            }
        ))
    );
    assert_eq!(
        ensemble_error.to_string(),
        "mmCIF write error: mmCIF writing requires 2 interpretation reports, but received 0"
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
fn mmcif_issue<T>(result: Result<T, FormatError>) -> Option<MmcifWriteIssue> {
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
    let model = model("C", &[[0.0, 0.0, 0.0]]);
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
            molfile::write_to(&mut FailingWriter, &model, MolfileWriteOptions::default())
                .unwrap_err()
        ),
        expected(Format::Molfile)
    );
    assert_eq!(
        sink_failure(
            sdf::write_to(&mut FailingWriter, [&model], SdfWriteOptions::default()).unwrap_err()
        ),
        expected(Format::Sdf)
    );
    assert_eq!(
        sink_failure(
            mmcif::write_to(
                &mut FailingWriter,
                [mmcif::MmcifBlockSource::model(&model)
                    .with_classifications(&classifications(&model))],
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
    let positioned = model("CCO", &[[0.0; 3], [1.5, 0.0, 0.0], [2.0, 1.0, 0.0]]);
    for source in [
        molfile::MolfileSource::molecule(&ethanol),
        molfile::MolfileSource::model(&positioned),
    ] {
        let text = molfile::write(source.with_title("ethanol 1"), Default::default()).unwrap();
        let document = molfile::parse_str(&text, Default::default()).unwrap();
        assert_eq!(document.header().title(), "ethanol 1");
        assert_eq!(
            described(molfile::write(source.with_title("a\nb"), Default::default()).unwrap_err()),
            (
                DiagnosticKind::Unsupported,
                Some(FormatErrorDetail::Title),
                "Molfile titles cannot contain line breaks".to_owned(),
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
                "SMILES record names must be nonempty, trimmed, and on one line".to_owned(),
            ),
            "{name:?}"
        );
    }
}
