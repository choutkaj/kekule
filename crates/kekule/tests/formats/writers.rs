use std::sync::Arc;
use std::{io, io::Write};

use kekule::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, FormatStage, StructureField,
    UnrepresentablePolicy, WriteOmission,
};
use kekule::macromolecular::EntityKind;
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::mmcif::{self, MmcifEntityClassifications, MmcifWriteIssue, MmcifWriteOptions};
use kekule::molfile::{self, MolfileWriteOptions, MolfileWriteVersion};
use kekule::pdb::{self, PdbWriteOptions};
use kekule::sdf::{self, SdfDataField, SdfRecordInterpretation, SdfWriteOptions};
use kekule::smiles::{self, SmilesWriteMode, SmilesWriteOptions};
use kekule::structure::{
    Distribution, Forces, Geometry, Members, Order, Partition, Positions, State, Structure,
    StructureSet, StructureView, TimeAxis, Velocities,
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

fn model(source: &str, points: &[[f64; 3]]) -> Structure {
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

fn classifications(model: &Structure) -> MmcifEntityClassifications {
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
    let charged = model(
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
            &model(source, &layout),
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
    let enhanced = model("C[C@H](N)O |&1:1|", &layout);
    let sodium = molecule("[Na+]");
    let large = Structure::new(
        Topology::from_molecules(vec![sodium; 1_000]).unwrap(),
        Positions::zeros(1_000),
    )
    .unwrap();
    for model in [&enhanced, &large] {
        let explicit = molfile::write(
            model,
            MolfileWriteOptions::default().with_version(MolfileWriteVersion::V3000),
        )
        .unwrap()
        .0;
        assert_eq!(
            molfile::write(model, MolfileWriteOptions::default())
                .unwrap()
                .0,
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
    let text = molfile::write(&model, MolfileWriteOptions::default())
        .unwrap()
        .0;
    assert!(text.contains("V2000"));

    let interpreted = (molfile::parse_str(&text, Default::default()).unwrap())
        .interpret(Default::default())
        .unwrap();
    let round_trip_model = interpreted.into_structure();
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
    let model = Structure::new(topology, Positions::zeros(1_000)).unwrap();

    let error = molfile::write(
        &model,
        MolfileWriteOptions::default().with_version(MolfileWriteVersion::V2000),
    )
    .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);

    let promoted = molfile::write(&model, MolfileWriteOptions::default())
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
    let raised_drawing = model("C[C@H](N)O", &layout);
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
fn sdf_writes_independent_models_and_set_samples_in_order() {
    let first = model("C", &[[1.0, 0.0, 0.0]]);
    let second = model("O", &[[2.0, 0.0, 0.0]]);
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
    let set = distributed([first, later], Partition::one_state(2));
    let set_text = sdf::write(&set, SdfWriteOptions::default()).unwrap().0;
    assert_eq!(crate::fingerprint(&set_text), 0x5a81_ddae_c2c0_94b9);
    let records = interpret_records(&set_text);
    assert_eq!(records.len(), 2);
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

/// `models` as the geometries of one distribution in the states of
/// `partition`, every member weighing 1.
fn distributed(models: impl IntoIterator<Item = Structure>, partition: Partition) -> StructureSet {
    let mut set = StructureSet::from_structures(models).unwrap();
    let members = set.len();
    set.set_distribution(Some(
        Distribution::new(Members::All, partition)
            .unwrap()
            .with_weights(vec![1.0; members])
            .unwrap(),
    ))
    .unwrap();
    set
}

/// Two methane conformers, one state each, weighted 1:3 with state energies
/// 0 and 5 kJ/mol.
fn weighted_set() -> StructureSet {
    let first = model("C", &[[1.0, 0.0, 0.0]]);
    let later = Structure::new(
        first.shared_topology(),
        Positions::new(Quantity::new([Point3::new(9.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    let mut set = distributed([first, later], Partition::state_per_sample(2));
    let mut distribution = set.distribution_mut().unwrap();
    distribution.set_weights(Some(vec![1.0, 3.0])).unwrap();
    for (state, energy) in [(0, 0.0), (1, 5.0)] {
        distribution
            .state_mut(state)
            .unwrap()
            .set_energy(Some(Quantity::new(energy, KILOJOULE_PER_MOLE)))
            .unwrap();
    }
    set
}

/// Four methane geometries weighted 1 to 4 in two states, {1, 3} represented
/// by 3 and {0, 2} represented by 2, the latter at -2.5 kJ/mol, at 300 K.
fn clustered_set() -> StructureSet {
    let first = model("C", &[[0.0, 0.0, 0.0]]);
    let models = (0..4).map(|index| {
        Structure::new(
            first.shared_topology(),
            Positions::new(Quantity::new(
                [Point3::new(f64::from(index), 0.0, 0.0)],
                ANGSTROM,
            ))
            .unwrap(),
        )
        .unwrap()
    });
    let partition = Partition::new(vec![1, 0, 1, 0], vec![3, 2]).unwrap();
    let mut set = distributed(models.collect::<Vec<_>>(), partition);
    let mut distribution = set.distribution_mut().unwrap();
    distribution
        .set_weights(Some(vec![1.0, 2.0, 3.0, 4.0]))
        .unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .set_energy(Some(Quantity::new(-2.5, KILOJOULE_PER_MOLE)))
        .unwrap();
    distribution
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    set
}

/// Every geometry's weight, the partition, every state energy in kJ/mol, and
/// the temperature in kelvin.
type SetMetadata = (
    Vec<Option<f64>>,
    Option<Partition>,
    Vec<Option<f64>>,
    Option<f64>,
);

fn set_metadata(set: &StructureSet) -> SetMetadata {
    let distribution = set.distribution();
    (
        set.iter().map(|geometry| geometry.weight()).collect(),
        distribution.map(|distribution| distribution.partition().clone()),
        set.states()
            .map(|state| {
                state
                    .energy()
                    .map(|energy| energy.value_in(KILOJOULE_PER_MOLE).unwrap())
            })
            .collect(),
        distribution
            .and_then(Distribution::temperature)
            .map(|temperature| temperature.value_in(KELVIN).unwrap()),
    )
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
fn sdf_samples_carry_their_weight_state_and_energy_through_data_fields() {
    let set = weighted_set();
    let text = sdf::write(&set, SdfWriteOptions::default()).unwrap().0;
    assert_eq!(crate::fingerprint(&text), 0x335f_7889_a7ba_d12f);
    let records = interpret_records(&text);
    let field = |name: &str, value: &str| (name.to_owned(), value.to_owned());
    assert_eq!(
        data_fields(&records),
        [
            vec![
                field(sdf::WEIGHT_FIELD, "1"),
                field(sdf::STATE_FIELD, "0"),
                field(sdf::REPRESENTATIVE_FIELD, "true"),
                field(sdf::ENERGY_FIELD, "0"),
            ],
            vec![
                field(sdf::WEIGHT_FIELD, "3"),
                field(sdf::STATE_FIELD, "1"),
                field(sdf::REPRESENTATIVE_FIELD, "true"),
                field(sdf::ENERGY_FIELD, "5"),
            ],
        ]
    );
    let reread = sdf::structure_set(&records).unwrap();
    assert_eq!(set_metadata(&reread), set_metadata(&set));
    assert_eq!(crate::native_fingerprint(&reread), 0x4ae5_cd3f_eab9_4c61);
    assert_eq!(
        reread.distribution().unwrap().probabilities().unwrap(),
        [0.25, 0.75]
    );

    // Records without state fields, as other programs write conformers, are
    // one state each with their energies.
    let foreign = text
        .replace(">  <kekule.state>\n0\n\n", "")
        .replace(">  <kekule.state>\n1\n\n", "")
        .replace(">  <kekule.representative>\ntrue\n\n", "");
    let read = sdf::structure_set(&interpret_records(&foreign)).unwrap();
    assert_eq!(set_metadata(&read), set_metadata(&set));
    assert_eq!(crate::native_fingerprint(&read), 0x4ae5_cd3f_eab9_4c61);
    // Records without any of the fields have no distribution.
    let plain = interpret_records(
        &sdf::write(&[model("C", &[[0.0; 3]])], SdfWriteOptions::default())
            .unwrap()
            .0,
    );
    let read = sdf::structure_set(&plain).unwrap();
    assert_eq!(read.distribution(), None);
    assert_eq!(crate::native_fingerprint(&read), 0x6708_332a_1f3d_942d);

    let malformed = text.replacen(">  <kekule.weight>\n3\n", ">  <kekule.weight>\nheavy\n", 1);
    assert_eq!(
        sdf::structure_set(&interpret_records(&malformed))
            .unwrap_err()
            .to_string(),
        "Cannot interpret the SDF input in record 2. The data field \"kekule.weight\" holds \"heavy\", not a number."
    );
    // An invalid weight names its record alone.
    let negative = text.replacen(">  <kekule.weight>\n3\n", ">  <kekule.weight>\n-3\n", 1);
    let error = sdf::structure_set(&interpret_records(&negative)).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Cannot interpret the SDF input in record 2. The data field \"kekule.weight\" holds an invalid weight."
    );
    assert!(std::error::Error::source(&error).is_none());
    let fields = [SdfDataField::new(sdf::STATE_FIELD, "2")];
    let geometry = sdf::SdfRecordSource::from(set.get(0).unwrap()).with_data_fields(&fields);
    assert_eq!(
        described(sdf::write([geometry], SdfWriteOptions::default()).unwrap_err()),
        (
            DiagnosticKind::Inconsistent,
            Some(FormatErrorDetail::DataField {
                name: sdf::STATE_FIELD.to_owned(),
            }),
            "The data field \"kekule.state\" is reserved for the data of the geometry.".to_owned(),
        )
    );
    // So does a field the geometry does not write, such as a second
    // representative marker on a geometry that represents nothing.
    let clustered = clustered_set();
    let fields = [SdfDataField::new(sdf::REPRESENTATIVE_FIELD, "true")];
    let geometry = sdf::SdfRecordSource::from(clustered.get(0).unwrap()).with_data_fields(&fields);
    assert_eq!(
        described(sdf::write([geometry], SdfWriteOptions::default()).unwrap_err()).1,
        Some(FormatErrorDetail::DataField {
            name: sdf::REPRESENTATIVE_FIELD.to_owned(),
        })
    );
}

#[test]
fn sdf_round_trips_one_state_and_clustered_states() {
    let set = clustered_set();
    let text = sdf::write(&set, SdfWriteOptions::default()).unwrap().0;
    assert_eq!(crate::fingerprint(&text), 0x7c1a_5784_711d_b34e);
    let records = interpret_records(&text);
    let field = |name: &str, value: &str| (name.to_owned(), value.to_owned());
    // Only representatives carry their state's energy; every record carries
    // the temperature.
    assert_eq!(
        data_fields(&records)[2],
        [
            field(sdf::WEIGHT_FIELD, "3"),
            field(sdf::STATE_FIELD, "1"),
            field(sdf::REPRESENTATIVE_FIELD, "true"),
            field(sdf::ENERGY_FIELD, "-2.5"),
            field(sdf::TEMPERATURE_FIELD, "300"),
        ]
    );
    assert_eq!(
        data_fields(&records)[0],
        [
            field(sdf::WEIGHT_FIELD, "1"),
            field(sdf::STATE_FIELD, "1"),
            field(sdf::TEMPERATURE_FIELD, "300"),
        ]
    );
    let reread = sdf::structure_set(&records).unwrap();
    assert_eq!(set_metadata(&reread), set_metadata(&set));
    assert_eq!(crate::native_fingerprint(&reread), 0x1460_5a70_893d_3b93);
    for (left, right) in reread.iter().zip(set.iter()) {
        assert_eq!(left.positions(), right.positions());
    }

    // One state.
    let mut one = clustered_set();
    let mut distribution = one.distribution_mut().unwrap();
    distribution.regroup(Partition::one_state(4)).unwrap();
    distribution.set_temperature(None).unwrap();
    let text_one = sdf::write(&one, SdfWriteOptions::default()).unwrap().0;
    assert_eq!(crate::fingerprint(&text_one), 0x8638_8572_626a_81dc);
    let reread = sdf::structure_set(&interpret_records(&text_one)).unwrap();
    assert_eq!(set_metadata(&reread), set_metadata(&one));
    assert_eq!(crate::native_fingerprint(&reread), 0xca61_de19_5bcd_2aec);
    assert_eq!(reread.distribution().unwrap().state_count(), 1);

    // Records that do not describe a partition are rejected.
    let error = |text: &str| {
        sdf::structure_set(&interpret_records(text))
            .unwrap_err()
            .to_string()
    };
    assert_eq!(
        error(&text.replacen(">  <kekule.representative>\ntrue\n\n", "", 1)),
        "Cannot interpret the SDF input in record 1. The data field \"kekule.representative\" \
         marks no record of state 1."
    );
    assert_eq!(
        error(&text.replacen(">  <kekule.state>\n1\n\n", "", 1)),
        "Cannot interpret the SDF input in record 1. The data field \"kekule.state\" is \
         missing. Other records have it."
    );
    assert_eq!(
        error(&text.replacen(
            ">  <kekule.weight>\n1\n\n",
            ">  <kekule.weight>\n1\n\n>  <kekule.energy.kJ_per_mol>\n4\n\n",
            1
        )),
        "Cannot interpret the SDF input in record 1. The data field \
         \"kekule.energy.kJ_per_mol\" is on a record that does not represent its state."
    );
    assert_eq!(
        error(&text.replacen(
            ">  <kekule.temperature.K>\n300\n",
            ">  <kekule.temperature.K>\n310\n",
            1
        )),
        "Cannot interpret the SDF input in record 2. The data field \"kekule.temperature.K\" \
         differs from that of record 1."
    );
}

#[test]
fn mmcif_rejects_the_weights_states_energies_and_temperature_it_cannot_carry() {
    use StructureField::{Energy, State, Temperature, Weight};
    let mut set = weighted_set();
    let message = |set: &StructureSet| {
        mmcif::write([set], MmcifWriteOptions::default())
            .map(|(text, _)| text)
            .unwrap_err()
            .to_string()
    };
    assert_eq!(
        message(&set),
        "Cannot write the mmCIF output. The state energy of geometry 1 has no mmCIF representation."
    );
    let mut distribution = set.distribution_mut().unwrap();
    for state in 0..2 {
        distribution
            .state_mut(state)
            .unwrap()
            .set_energy(None)
            .unwrap();
    }
    assert_eq!(
        message(&set),
        "Cannot write the mmCIF output. The weight of geometry 2 has no mmCIF representation."
    );
    set.distribution_mut().unwrap().set_weight(1, 1.0).unwrap();
    assert_eq!(
        message(&set),
        "Cannot write the mmCIF output. The state of geometry 2 has no mmCIF representation."
    );
    let mut distribution = set.distribution_mut().unwrap();
    distribution.regroup(Partition::one_state(2)).unwrap();
    distribution
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    assert_eq!(
        message(&set),
        "Cannot write the mmCIF output. The ensemble temperature has no mmCIF representation."
    );
    set.distribution_mut()
        .unwrap()
        .set_temperature(None)
        .unwrap();
    let (omitted, report) = mmcif::write(
        [&clustered_set()],
        MmcifWriteOptions::default().with_set_metadata(UnrepresentablePolicy::Omit),
    )
    .unwrap();
    assert_eq!(
        report.omissions(),
        [Weight, Energy, State, Temperature].map(|field| WriteOmission::StructureState {
            field,
            structures: if field == Energy { 2 } else { 4 },
        })
    );
    assert_eq!(crate::fingerprint(&omitted), 0xa348_4d3d_7813_ac6a);
    let reread = mmcif::parse_str(&omitted, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(crate::native_fingerprint(&reread), 0xab5b_7b59_ab44_b17c);
    assert_eq!(
        set_metadata(reread.structure_set()),
        (
            vec![Some(1.0); 4],
            Some(Partition::one_state(4)),
            vec![None],
            None
        )
    );

    // One state of equal weights without an energy is what interpretation
    // gives back.
    mmcif::write([&set], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
}

/// One-atom methane frames with exactly one piece of frame state each:
/// none, time, MD step, velocities, forces.
struct StatefulFrames {
    topology: Arc<Topology>,
    geometries: Vec<Geometry>,
}

impl StatefulFrames {
    fn new() -> Self {
        let topology = model("C", &[[0.0; 3]]).shared_topology();
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

/// Writes one structure view as a Molfile, an SDF record and an mmCIF model
/// block, each with geometry state set to `policy`.
fn write_structure(
    view: StructureView<'_>,
    policy: UnrepresentablePolicy,
) -> [Result<String, FormatError>; 3] {
    [
        molfile::write(
            molfile::MolfileSource::structure(&view),
            MolfileWriteOptions::default().with_set_metadata(policy),
        )
        .map(|(text, _)| text),
        sdf::write(
            [sdf::SdfRecordSource::structure(&view)],
            SdfWriteOptions::default().with_set_metadata(policy),
        )
        .map(|(text, _)| text),
        mmcif::write(
            [mmcif::MmcifBlockSource::structure(&view)],
            MmcifWriteOptions::default().with_set_metadata(policy),
        )
        .map(|(text, _)| text),
    ]
}

#[test]
fn structure_writers_fail_on_sample_and_frame_state_unless_told_to_omit_it() {
    use StructureField::{
        Energy, Forces as ForceField, Step, Time, Velocities as VelocityField, Weight,
    };
    let mut set = weighted_set();
    set.distribution_mut()
        .unwrap()
        .state_mut(1)
        .unwrap()
        .set_energy(None)
        .unwrap();
    let frames = StatefulFrames::new();
    let frames = frames.views();
    // Expected detail per writer (Molfile, SDF, mmCIF); SDF writes geometry
    // weights, states, and state energies as data fields.
    let rows: [(StructureView<'_>, [Option<StructureField>; 3]); 7] = [
        (
            set.get(0).unwrap().as_structure_view(),
            [Some(Energy), None, Some(Energy)],
        ),
        (
            set.get(1).unwrap().as_structure_view(),
            [Some(Weight), None, Some(Weight)],
        ),
        (frames[0].as_structure_view(), [None; 3]),
        (frames[1].as_structure_view(), [Some(Time); 3]),
        (frames[2].as_structure_view(), [Some(Step); 3]),
        (frames[3].as_structure_view(), [Some(VelocityField); 3]),
        (frames[4].as_structure_view(), [Some(ForceField); 3]),
    ];
    let formats = [Format::Molfile, Format::Sdf, Format::Mmcif];
    for (row, (view, expected)) in rows.into_iter().enumerate() {
        let model = view.to_structure();
        let plain = write_structure(model.as_structure_view(), UnrepresentablePolicy::Fail);
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
                    // Both SDF paths write the geometry's own data fields.
                    let geometry = set.get(row).unwrap();
                    let typed = sdf::write([geometry], SdfWriteOptions::default())
                        .unwrap()
                        .0;
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
    let topology = model("C", &[[0.0; 3]]).shared_topology();
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
fn sdf_ordered_samples_carry_their_step_or_coordinate_through_data_fields() {
    let (history, scan) = methane_series();
    let field = |name: &str, value: &str| (name.to_owned(), value.to_owned());
    let text = sdf::write(&history, SdfWriteOptions::default()).unwrap().0;
    assert_eq!(crate::fingerprint(&text), 0x2543_22fa_c15f_63d1);
    assert_eq!(
        data_fields(&interpret_records(&text)),
        [
            vec![field(sdf::STEP_FIELD, "0")],
            vec![field(sdf::STEP_FIELD, "4")],
            vec![field(sdf::STEP_FIELD, "9")],
        ]
    );
    let text = sdf::write(&scan, SdfWriteOptions::default()).unwrap().0;
    assert_eq!(crate::fingerprint(&text), 0x5065_032e_d48f_d515);
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
    assert_eq!(crate::fingerprint(&text), 0x3b6b_1074_f00f_9ca7);
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
fn single_structure_writers_fail_on_a_sample_place_unless_told_to_omit_it() {
    let (history, scan) = methane_series();
    let omit = UnrepresentablePolicy::Omit;
    for (set, field, name) in [
        (&history, StructureField::OrderStep, "order step"),
        (&scan, StructureField::OrderCoordinate, "order coordinate"),
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
        // Omitting the place writes exactly the bare model.
        assert_eq!(
            molfile::write(
                molfile::MolfileSource::structure(&view),
                MolfileWriteOptions::default().with_set_metadata(omit),
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
            MmcifWriteOptions::default().with_set_metadata(omit),
        )
        .unwrap();
        assert_eq!(
            report.omissions(),
            [WriteOmission::StructureState {
                field,
                structures: 1
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
fn boxed_model(periodic: [bool; 3]) -> Structure {
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
    let message = "A Molfile cannot hold a periodic cell.";
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
        format!("Cannot write the SDF output in record 1. {message}")
    );
    let omitted = SdfWriteOptions::default().with_cell(UnrepresentablePolicy::Omit);
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

    let mut set = distributed([boxed.clone(), boxed], Partition::one_state(2));
    let text = mmcif::write([&set], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    assert_eq!(crate::fingerprint(&text), 0x37db_812f_437b_970d);
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
    let (text, _) = mmcif::write([&set], omit).unwrap();
    assert_eq!(crate::fingerprint(&text), 0xd84b_4cef_cb89_042c);
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
        model("CO.O", &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
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
    let first = model("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let second = model("N", &[[3.0, 0.0, 0.0]]);
    let third = model("O", &[[4.0, 0.0, 0.0]]);
    let source = |model, name: Option<&'static str>| {
        let block = mmcif::MmcifBlockSource::structure(model);
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
fn mmcif_distinguishes_independent_models_from_a_structure_set() {
    let first = model("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let unrelated = model("N", &[[3.0, 0.0, 0.0]]);
    let independent = mmcif::write(
        [first.clone(), unrelated.clone()]
            .iter()
            .zip(&[classifications(&first), classifications(&unrelated)])
            .map(|(model, classifications)| {
                mmcif::MmcifBlockSource::structure(model).with_classifications(classifications)
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
    let set = distributed([first.clone(), second], Partition::one_state(2));
    let set_text = mmcif::write(
        [mmcif::MmcifBlockSource::structure_set(&set)
            .with_classifications(&classifications(&first))],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    assert_eq!(crate::fingerprint(&set_text), 0xd958_9ed3_ad5a_420d);
    let set_document = mmcif::parse_str(&set_text, Default::default()).unwrap();
    assert_eq!(set_document.blocks().len(), 1);
    let parsed = set_document
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(crate::native_fingerprint(&parsed), 0x7b46_6ec6_7f47_b53d);
    assert_eq!(parsed.structure_set().len(), 2);
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
    assert_eq!(crate::fingerprint(&round_trip), 0xd958_9ed3_ad5a_420d);
    let reparsed = (mmcif::parse_str(&round_trip, Default::default()).unwrap())
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(reparsed.structure_set().len(), 2);
}

#[test]
fn mmcif_model_is_one_block_with_coordinate_model_one_and_automatic_classification() {
    let model = model("C", &[[4.0, 5.0, 6.0]]);
    let automatic = mmcif::write([&model], MmcifWriteOptions::default())
        .map(|(text, _)| text)
        .unwrap();
    assert!(automatic.contains("1 non-polymer"));
    let text =
        mmcif::write(
            [mmcif::MmcifBlockSource::structure(&model)
                .with_classifications(&classifications(&model))],
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
    let mut model = Structure::new(builder.build().unwrap(), positions).unwrap();
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
    let labels = |model: &Structure| {
        model
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
    let model = Structure::new(builder.build().unwrap(), Positions::zeros(2)).unwrap();
    assert_eq!(
        mmcif_issue(mmcif::write([&model], MmcifWriteOptions::default()).map(|(text, _)| text)),
        Some(MmcifWriteIssue::MissingAtomSite(InstanceAtomId::new(
            instance, atoms[0]
        )))
    );
}

#[test]
fn mmcif_streaming_matches_string_output_for_models_and_a_structure_set() {
    let first = model("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let unrelated = model("N", &[[3.0, 0.0, 0.0]]);
    let models = [first.clone(), unrelated.clone()];
    let model_classifications = [classifications(&first), classifications(&unrelated)];
    let expected_models = mmcif::write(
        models
            .iter()
            .zip(&model_classifications)
            .map(|(model, classifications)| {
                mmcif::MmcifBlockSource::structure(model).with_classifications(classifications)
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
                mmcif::MmcifBlockSource::structure(model).with_classifications(classifications)
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

    let second = Structure::new(
        first.shared_topology(),
        Positions::new(Quantity::new(
            [Point3::new(10.0, 0.0, 0.0), Point3::new(20.0, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    let set = distributed([first.clone(), second], Partition::one_state(2));
    let expected_set = mmcif::write(
        [mmcif::MmcifBlockSource::structure_set(&set)
            .with_classifications(&classifications(&first))],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    assert_eq!(crate::fingerprint(&expected_set), 0xd958_9ed3_ad5a_420d);
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
    let first = model("C", &[[1.0, 0.0, 0.0]]);
    let second = Structure::new(
        first.shared_topology(),
        Positions::new(Quantity::new([Point3::new(2.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();

    // Each model block carries its own reports, so a mismatch is reported
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
            MmcifWriteIssue::ReportCountMismatch {
                expected: 1,
                actual: 0
            }
        ))
    );
    assert_eq!(
        structure_error.to_string(),
        "Cannot write the mmCIF output. There are 0 interpretation reports for 1 structures."
    );

    let set = distributed([first, second], Partition::one_state(2));
    let set_error = mmcif::write(
        [mmcif::MmcifBlockSource::structure_set(&set).with_reports(&[])],
        MmcifWriteOptions::default(),
    )
    .map(|(text, _)| text)
    .unwrap_err();
    assert_eq!(
        set_error.detail(),
        Some(&FormatErrorDetail::MmcifWrite(
            MmcifWriteIssue::ReportCountMismatch {
                expected: 2,
                actual: 0
            }
        ))
    );
    assert_eq!(
        set_error.to_string(),
        "Cannot write the mmCIF output. There are 0 interpretation reports for 2 structures."
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
                [mmcif::MmcifBlockSource::structure(&model)
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

/// A write's outcome: its output, or its error's kind and detail.
fn outcome<T: std::fmt::Debug>(result: Result<T, FormatError>) -> String {
    match result {
        Ok(value) => format!("Ok({value:?})"),
        Err(error) => format!("Err({:?}, {:?})", error.kind(), error.detail()),
    }
}

/// The fingerprint of every structure writer's outcome for one view with
/// geometry state set to `policy`.
fn view_outcomes(view: StructureView<'_>, policy: UnrepresentablePolicy) -> u64 {
    crate::fingerprint(
        [
            outcome(
                molfile::write(
                    molfile::MolfileSource::structure(&view),
                    MolfileWriteOptions::default().with_set_metadata(policy),
                )
                .map(|(text, _)| text),
            ),
            outcome(
                sdf::write(
                    [sdf::SdfRecordSource::structure(&view)],
                    SdfWriteOptions::default().with_set_metadata(policy),
                )
                .map(|(text, _)| text),
            ),
            outcome(mmcif::write(
                [mmcif::MmcifBlockSource::structure(&view)],
                MmcifWriteOptions::default().with_set_metadata(policy),
            )),
            outcome(pdb::write(
                pdb::PdbSource::structure(&view),
                PdbWriteOptions::default().with_set_metadata(policy),
            )),
        ]
        .join("\n"),
    )
}

#[test]
fn sample_and_frame_views_write_their_state_or_fail_on_it() {
    let mut set = weighted_set();
    set.distribution_mut()
        .unwrap()
        .state_mut(1)
        .unwrap()
        .set_energy(None)
        .unwrap();
    let frames = StatefulFrames::new();
    let (history, scan) = methane_series();
    let views = set
        .iter()
        .map(|geometry| geometry.as_structure_view())
        .chain(
            frames
                .views()
                .into_iter()
                .map(|frame| frame.as_structure_view()),
        )
        .chain(history.iter().map(|geometry| geometry.as_structure_view()))
        .chain(scan.iter().map(|geometry| geometry.as_structure_view()))
        .collect::<Vec<_>>();
    // Molfile, SDF, mmCIF, and PDB outcomes of two weighted geometries, the
    // five stateful frames, and the step and coordinate places.
    let expected = [
        (
            UnrepresentablePolicy::Fail,
            [
                0x3125_0e0b_6f72_259d,
                0x1d1c_ea3c_bea8_caef,
                0xc0c3_c9e3_caca_e18c,
                0x89aa_5014_99a0_f1ed,
                0x186a_5785_5fe1_0661,
                0x16a1_2c31_b4fc_c06d,
                0x611f_a810_d55a_1c81,
                0xfc2c_9016_0697_fed6,
                0x70c5_d730_0069_f862,
                0x2e43_247e_4bed_0415,
                0x0f09_5294_25e6_86a4,
                0x5e14_2e8b_6ca1_318e,
            ],
        ),
        (
            UnrepresentablePolicy::Omit,
            [
                0x518f_49a3_9158_8365,
                0xb2aa_49de_c2e8_003d,
                0xc0c3_c9e3_caca_e18c,
                0x18b2_7334_f2fe_c5e0,
                0x0eb2_c9cc_cfa3_1b58,
                0x956a_216a_f9fc_cd34,
                0xba77_94ba_b8a5_3100,
                0xbdb8_04cd_e3dd_956f,
                0xfd78_26a5_9e7f_6953,
                0xff8f_e988_4e22_efb6,
                0xe8f2_ff23_05a9_cefd,
                0x4397_116d_e774_ff2b,
            ],
        ),
    ];
    for (policy, fingerprints) in expected {
        assert_eq!(views.len(), fingerprints.len());
        for (index, (view, fingerprint)) in views.iter().zip(fingerprints).enumerate() {
            assert_eq!(
                view_outcomes(*view, policy),
                fingerprint,
                "view {index}, {policy:?}"
            );
        }
    }
}

/// Two methane geometries, both timed and stepped, the second with velocities.
fn timed_set() -> StructureSet {
    let topology = model("C", &[[0.0; 3]]).shared_topology();
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
fn plain_models() -> [Structure; 2] {
    let first = model("C", &[[1.0, 0.0, 0.0]]);
    let second = Structure::new(
        first.shared_topology(),
        Positions::new(Quantity::new([Point3::new(2.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    [first, second]
}

#[test]
fn structure_sets_write_their_distribution_and_order_or_fail_on_them() {
    let distributed_sets = [
        (
            distributed(plain_models(), Partition::one_state(2)),
            0xef05_4d2a_760d_29e1,
        ),
        (weighted_set(), 0x543d_8ee2_c93e_9194),
        (clustered_set(), 0x4f1c_d9da_c18c_3e27),
    ];
    let (history, scan) = methane_series();
    let timed = timed_set();
    for (set, expected) in distributed_sets.iter().chain(&[
        (history.clone(), 0x7a22_d2a6_9ef3_764a),
        (scan.clone(), 0xb228_c64b_8f05_5b2b),
        (timed.clone(), 0xea90_f0fe_abe7_d46a),
    ]) {
        assert_eq!(crate::native_fingerprint(set), *expected);
    }
    // PDB, mmCIF, and SDF outcomes of each distributed set; PDB and SDF of
    // each ordered one; SDF of the timed one.
    let expected = [
        (
            UnrepresentablePolicy::Fail,
            [
                0x4eb2_d1a7_b204_40e5,
                0x81d7_5a99_ec97_0c6a,
                0x2b8b_46ef_6646_0b65,
                0x8aba_052b_a262_3bbb,
                0xb666_49e6_a9dd_829d,
                0x33fb_5b41_17c5_491a,
                0x8aba_052b_a262_3bbb,
                0xb666_49e6_a9dd_829d,
                0x5009_da48_133d_9c9b,
                0xf512_9294_4996_18cb,
                0x0d1a_586e_2139_0dc0,
                0x5765_85f7_861d_ffbf,
                0x64ee_b24c_080e_d552,
                0xb8c0_bcdb_13d4_d055,
            ],
        ),
        (
            UnrepresentablePolicy::Omit,
            [
                0x4eb2_d1a7_b204_40e5,
                0x81d7_5a99_ec97_0c6a,
                0x2b8b_46ef_6646_0b65,
                0xb580_6c39_a7a0_1e73,
                0x2a35_a3da_6072_52a4,
                0x33fb_5b41_17c5_491a,
                0xa50e_e374_8865_7588,
                0x3273_d4c8_7d21_0bd0,
                0x5009_da48_133d_9c9b,
                0x28a9_a32a_9522_058f,
                0x0d1a_586e_2139_0dc0,
                0x0e9f_e5ca_171a_3e44,
                0x64ee_b24c_080e_d552,
                0x763c_4a75_cf56_36b4,
            ],
        ),
    ];
    for (policy, fingerprints) in expected {
        let pdb_options = PdbWriteOptions::default()
            .with_implicit_hydrogens(UnrepresentablePolicy::Omit)
            .with_set_metadata(policy);
        let mmcif_options = MmcifWriteOptions::default().with_set_metadata(policy);
        let sdf_options = SdfWriteOptions::default().with_set_metadata(policy);
        let mut outcomes = Vec::new();
        for (set, _) in &distributed_sets {
            outcomes.push(outcome(pdb::write(set, pdb_options)));
            outcomes.push(outcome(mmcif::write([set], mmcif_options.clone())));
            outcomes.push(outcome(sdf::write(set, sdf_options).map(|(text, _)| text)));
        }
        for set in [&history, &scan] {
            outcomes.push(outcome(pdb::write(set, pdb_options)));
            outcomes.push(outcome(sdf::write(set, sdf_options).map(|(text, _)| text)));
        }
        outcomes.push(outcome(
            sdf::write(&timed, sdf_options).map(|(text, _)| text),
        ));
        assert_eq!(outcomes.len(), fingerprints.len());
        for (index, (outcome, fingerprint)) in outcomes.iter().zip(fingerprints).enumerate() {
            assert_eq!(
                crate::fingerprint(outcome),
                fingerprint,
                "case {index}, {policy:?}: {outcome}"
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
        pdb_options.with_set_metadata(UnrepresentablePolicy::Omit),
    )
    .unwrap();
    let state = |field, structures| WriteOmission::StructureState { field, structures };
    assert!(report.omissions().starts_with(&[
        state(StructureField::Time, 2),
        state(StructureField::Step, 2),
        state(StructureField::Velocities, 1),
    ]));
    // The positions are the geometries' positions.
    let models = set
        .iter()
        .map(|geometry| geometry.to_structure())
        .collect::<Vec<_>>();
    let bare = distributed(models, Partition::one_state(2));
    assert_eq!(text, pdb::write(&bare, pdb_options).unwrap().0);

    // A geometry outside the distribution is a model all the same; the file
    // cannot say it is no member.
    let mut subset = bare.clone();
    subset
        .set_distribution(Some(
            Distribution::new(Members::Subset(vec![1]), Partition::one_state(1)).unwrap(),
        ))
        .unwrap();
    let error = mmcif::write([&subset], MmcifWriteOptions::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::MmcifWrite(
            MmcifWriteIssue::UnrepresentableGeometry {
                geometry: 1,
                field: StructureField::Membership,
            }
        ))
    );
    let (text, report) = mmcif::write(
        [&subset],
        MmcifWriteOptions::default().with_set_metadata(UnrepresentablePolicy::Omit),
    )
    .unwrap();
    assert_eq!(report.omissions(), [state(StructureField::Membership, 1)]);
    assert_eq!(
        text,
        mmcif::write([&bare], MmcifWriteOptions::default())
            .unwrap()
            .0
    );

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
fn a_sample_with_a_place_and_a_membership_writes_both_in_sdf() {
    let mut set = StructureSet::from_structures([model("C", &[[0.0; 3]])]).unwrap();
    set.set_order(Some(Order::Step(vec![7]))).unwrap();
    let state = State::new()
        .with_energy(Quantity::new(-1.5, KILOJOULE_PER_MOLE))
        .unwrap();
    let distribution = Distribution::new(Members::All, Partition::one_state(1))
        .unwrap()
        .with_weights(vec![2.5])
        .unwrap()
        .with_states(vec![state])
        .unwrap()
        .with_temperature(Quantity::new(300.0, KELVIN))
        .unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    let view = set.get(0).unwrap();
    let text = sdf::write(
        [sdf::SdfRecordSource::structure(&view)],
        SdfWriteOptions::default(),
    )
    .unwrap()
    .0;
    // The membership's fields, then the place's.
    let expected = "
  kekule            2D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <kekule.weight>
2.5

>  <kekule.state>
0

>  <kekule.representative>
true

>  <kekule.energy.kJ_per_mol>
-1.5

>  <kekule.temperature.K>
300

>  <kekule.step>
7

$$$$
";
    assert_eq!(text, expected);
}

#[test]
fn structure_sets_report_or_refuse_the_order_and_states_pdb_and_mmcif_lose() {
    use kekule::properties::{PropertyKey, PropertyValue};
    use kekule::structure::TimeAxis;
    let plain = distributed(plain_models(), Partition::one_state(2));
    let omit = UnrepresentablePolicy::Omit;
    let pdb_options = PdbWriteOptions::default().with_implicit_hydrogens(omit);
    let (pdb_text, pdb_report) = pdb::write(&plain, pdb_options).unwrap();
    let (mmcif_text, mmcif_report) = mmcif::write([&plain], MmcifWriteOptions::default()).unwrap();
    let state = |field, structures| WriteOmission::StructureState { field, structures };

    // A time order without times or steps.
    let mut bare = plain.clone();
    bare.set_order(Some(Order::Time(TimeAxis::new(None, None).unwrap())))
        .unwrap();
    // A state with properties.
    let mut labelled_state = State::new();
    labelled_state
        .properties_mut()
        .insert(
            PropertyKey::new("label").unwrap(),
            PropertyValue::Bool(true),
        )
        .unwrap();
    let mut labelled = plain.clone();
    labelled
        .set_distribution(Some(
            Distribution::new(Members::All, Partition::one_state(2))
                .unwrap()
                .with_states(vec![labelled_state])
                .unwrap(),
        ))
        .unwrap();
    // A lone state represented by its second member.
    let mut late = plain.clone();
    late.set_distribution(Some(
        Distribution::new(Members::All, Partition::new(vec![0, 0], vec![1]).unwrap()).unwrap(),
    ))
    .unwrap();

    let cases = [
        (&bare, StructureField::Order, 1, 2),
        (&labelled, StructureField::StateProperties, 1, 1),
        (&late, StructureField::Representative, 2, 1),
    ];
    for (set, field, geometry, structures) in cases {
        // The order belongs to the whole set.
        let lost = |format| match field {
            StructureField::Order => format!("The geometry order has no {format} representation."),
            StructureField::StateProperties => {
                format!("The {field} of geometry {geometry} have no {format} representation.")
            }
            _ => format!("The {field} of geometry {geometry} has no {format} representation."),
        };
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
                MmcifWriteIssue::UnrepresentableGeometry { geometry, field }
            )),
            "{error}"
        );

        let (text, report) = pdb::write(set, pdb_options.with_set_metadata(omit)).unwrap();
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
        let (text, report) =
            mmcif::write([set], MmcifWriteOptions::default().with_set_metadata(omit)).unwrap();
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

    // With several states, the lost grouping includes the representatives.
    let mut clustered = plain;
    clustered
        .set_distribution(Some(
            Distribution::new(
                Members::All,
                Partition::new(vec![1, 0], vec![1, 0]).unwrap(),
            )
            .unwrap(),
        ))
        .unwrap();
    let (_, report) = pdb::write(&clustered, pdb_options.with_set_metadata(omit)).unwrap();
    assert!(report
        .omissions()
        .contains(&state(StructureField::State, 2)));
    assert!(!report
        .omissions()
        .iter()
        .any(|omission| *omission == state(StructureField::Representative, 1)));
}

#[test]
fn mmcif_structure_sets_check_the_report_count_before_their_state() {
    // The set's state energies alone would fail to write.
    let set = weighted_set();
    let reports = Vec::new();
    let error = mmcif::write(
        [mmcif::MmcifBlockSource::structure_set(&set).with_reports(&reports)],
        MmcifWriteOptions::default(),
    )
    .unwrap_err();
    assert!(matches!(
        mmcif_issue::<()>(Err(error)),
        Some(MmcifWriteIssue::ReportCountMismatch {
            expected: 2,
            actual: 0
        })
    ),);
}

#[test]
fn structure_fields_name_the_state_a_format_loses() {
    use StructureField::*;
    let names = [
        (Weight, "weight"),
        (Energy, "state energy"),
        (State, "state"),
        (Temperature, "temperature"),
        (StateProperties, "state properties"),
        (Representative, "state representative"),
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
        (Membership, "distribution membership"),
        (SetProperties, "structure set properties"),
        (Occupancies, "occupancies"),
        (BFactors, "B factors"),
    ];
    for (field, name) in names {
        assert_eq!(field.to_string(), name, "{field:?}");
    }
}

#[test]
fn mmcif_samples_share_the_first_sample_s_entities() {
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
    let set = distributed(
        [from_pdb.structure().clone(), from_pdb.structure().clone()],
        Partition::one_state(2),
    );
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
            MmcifWriteIssue::IncompatibleGeometry {
                geometry: 2,
                field: kekule::mmcif::MmcifGeometryField::Entities,
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
        let model = Structure::new(topology, Positions::zeros(atoms)).unwrap();
        let error = if mmcif {
            mmcif::write([&model], MmcifWriteOptions::default()).map(|(text, _)| text)
        } else {
            pdb::write(&model, PdbWriteOptions::default()).map(|(text, _)| text)
        };
        error.unwrap_err().to_string()
    };
    assert_eq!(
        message("[CH4]", true),
        "Cannot write the mmCIF output. The implicit hydrogens of atom 0 of instance 0 have no mmCIF representation."
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
    use StructureField::{Energy, Forces as ForceField, Step, Temperature, Time};
    use StructureField::{Velocities as VelocityField, Weight};
    let omit = UnrepresentablePolicy::Omit;
    let state = |field, structures| WriteOmission::StructureState { field, structures };
    let plain = model("C", &[[0.0; 3]]);
    let (_, report) = sdf::write([&plain, &plain], SdfWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    let (_, report) = molfile::write(&plain, MolfileWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    let (_, report) = molfile::write(&molecule("CCO"), MolfileWriteOptions::default()).unwrap();
    assert!(report.is_empty());

    // The cells of the records written without them.
    let boxed = boxed_model([true; 3]);
    let open = SdfWriteOptions::default().with_cell(omit);
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
    let metadata = SdfWriteOptions::default().with_set_metadata(omit);
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
    let molfile_metadata = MolfileWriteOptions::default().with_set_metadata(omit);
    for (frame, field) in frames
        .iter()
        .skip(1)
        .zip([Time, Step, VelocityField, ForceField])
    {
        let (_, report) =
            molfile::write(molfile::MolfileSource::structure(frame), molfile_metadata).unwrap();
        assert_eq!(report.omissions(), [state(field, 1)]);
    }
    // A Molfile of a geometry cannot carry its membership, which SDF writes.
    let mut weighted = weighted_set();
    weighted
        .distribution_mut()
        .unwrap()
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    let (_, report) = molfile::write(
        molfile::MolfileSource::structure(&weighted.get(0).unwrap()),
        molfile_metadata,
    )
    .unwrap();
    assert_eq!(
        report.omissions(),
        [state(Weight, 1), state(Energy, 1), state(Temperature, 1)]
    );
    let (_, report) = sdf::write(&weighted, SdfWriteOptions::default()).unwrap();
    assert!(report.is_empty());
}

/// `count` methane geometries, geometry `i` at `i` Å along x.
fn methane_set(count: usize) -> StructureSet {
    let first = model("C", &[[0.0, 0.0, 0.0]]);
    let models = (0..count).map(|index| {
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
    StructureSet::from_structures(models).unwrap()
}

/// Four methane geometries with a distribution over geometries 1 and 3 in one
/// state each, weighted, with state energies and a temperature.
fn subset_set() -> StructureSet {
    let mut set = methane_set(4);
    let states = vec![
        State::new()
            .with_energy(Quantity::new(-1.5, KILOJOULE_PER_MOLE))
            .unwrap(),
        State::new()
            .with_energy(Quantity::new(0.1 + 0.2, KILOJOULE_PER_MOLE))
            .unwrap(),
    ];
    set.set_distribution(Some(
        Distribution::new(Members::Subset(vec![1, 3]), Partition::state_per_sample(2))
            .unwrap()
            .with_weights(vec![0.25, 1.0 / 3.0])
            .unwrap()
            .with_states(states)
            .unwrap()
            .with_temperature(Quantity::new(298.15, KELVIN))
            .unwrap(),
    ))
    .unwrap();
    set
}

/// The set the records of `text` form.
fn read_back(text: &str) -> Result<StructureSet, FormatError> {
    sdf::structure_set(&interpret_records(text))
}

/// The data fields of a record outside the distribution.
fn outside_fields() -> Vec<(String, String)> {
    vec![(sdf::MEMBER_FIELD.to_owned(), "false".to_owned())]
}

#[test]
fn sdf_round_trips_a_distribution_over_some_samples() {
    let set = subset_set();
    let (text, report) = sdf::write(&set, SdfWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    let records = interpret_records(&text);
    assert_eq!(records.len(), 4);
    // The non-members' records are marked and carry nothing else.
    let fields = data_fields(&records);
    assert_eq!(fields[0], outside_fields());
    assert_eq!(fields[2], outside_fields());
    assert!(fields[1].iter().any(|(name, _)| name == sdf::WEIGHT_FIELD));
    let read = sdf::structure_set(&records).unwrap();
    assert_eq!(read.len(), 4);
    assert_eq!(read.distribution(), set.distribution());
    let distribution = read.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::Subset(vec![1, 3]));
    assert_eq!(
        distribution.weights().unwrap()[1].to_bits(),
        (1.0_f64 / 3.0).to_bits()
    );
    assert_eq!(distribution.weight(0), None);
    // Writing the set read back gives the same bytes.
    assert_eq!(
        sdf::write(&read, SdfWriteOptions::default()).unwrap().0,
        text
    );

    // A geometry of a set with a distribution owns every distribution field
    // name, whether it is a member or not.
    for (geometry, name) in [(0, sdf::STATE_FIELD), (1, sdf::MEMBER_FIELD)] {
        let fields = [SdfDataField::new(name, "0")];
        let record =
            sdf::SdfRecordSource::from(set.get(geometry).unwrap()).with_data_fields(&fields);
        assert_eq!(
            described(sdf::write([record], SdfWriteOptions::default()).unwrap_err()),
            (
                DiagnosticKind::Inconsistent,
                Some(FormatErrorDetail::DataField {
                    name: name.to_owned()
                }),
                format!("The data field {name:?} is reserved for the data of the geometry."),
            ),
            "{name}"
        );
    }
}

#[test]
fn sdf_takes_only_a_marked_record_as_outside_the_distribution() {
    let error =
        |records: &[SdfRecordInterpretation]| sdf::structure_set(records).unwrap_err().to_string();
    // A record stripped of its fields among the records of a distribution
    // over every geometry is refused, not taken as outside.
    let set = distributed(
        methane_set(3)
            .iter()
            .map(|geometry| geometry.to_structure()),
        Partition::state_per_sample(3),
    );
    let text = sdf::write(&set, SdfWriteOptions::default()).unwrap().0;
    let mut records = interpret_records(&text);
    records[1] = SdfRecordInterpretation::new("", records[1].structure().clone(), Vec::new());
    assert_eq!(
        error(&records),
        "Cannot interpret the SDF input in record 2. The data field \"kekule.state\" is \
         missing. Other records have it."
    );
    // So is a record of another program appended to a written file.
    let mut records = interpret_records(&text);
    records.push(SdfRecordInterpretation::new(
        "",
        records[0].structure().clone(),
        Vec::new(),
    ));
    assert_eq!(
        error(&records),
        "Cannot interpret the SDF input in record 4. The data field \"kekule.state\" is \
         missing. Other records have it."
    );

    // The marker is `false`, on a record with no other distribution field.
    let text = sdf::write(&subset_set(), SdfWriteOptions::default())
        .unwrap()
        .0;
    let marker = ">  <kekule.member>\nfalse\n\n";
    let marked = |replacement: &str| interpret_records(&text.replacen(marker, replacement, 1));
    assert_eq!(
        error(&marked(">  <kekule.member>\ntrue\n\n")),
        "Cannot interpret the SDF input in record 1. The data field \"kekule.member\" holds \
         \"true\", not \"false\"."
    );
    assert_eq!(
        error(&marked(
            ">  <kekule.member>\nfalse\n\n>  <kekule.weight>\n1\n\n"
        )),
        "Cannot interpret the SDF input in record 1. The data field \"kekule.weight\" is on \
         a record outside the distribution."
    );

    // A temperature that differs names the first member's record.
    let temperature = ">  <kekule.temperature.K>\n298.15\n";
    let last = text.rfind(temperature).unwrap();
    let mut changed = text.clone();
    changed.replace_range(
        last..last + temperature.len(),
        ">  <kekule.temperature.K>\n310\n",
    );
    assert_eq!(
        error(&interpret_records(&changed)),
        "Cannot interpret the SDF input in record 4. The data field \"kekule.temperature.K\" \
         differs from that of record 2."
    );
}

#[test]
fn sdf_writes_and_reads_records_outside_the_distribution_alone() {
    let set = subset_set();
    let outside = [set.get(0).unwrap(), set.get(2).unwrap()];
    let (text, report) = sdf::write(outside, SdfWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    let records = interpret_records(&text);
    assert_eq!(data_fields(&records), [outside_fields(), outside_fields()]);
    // Record by record the output reads and writes back the same.
    assert_eq!(
        sdf::write(&records, SdfWriteOptions::default()).unwrap().0,
        text
    );
    // The records form a set whose distribution has no members, which
    // writes back the same.
    let read = sdf::structure_set(&records).unwrap();
    assert_eq!(read.len(), 2);
    assert_eq!(
        read.distribution().unwrap().members(),
        &Members::Subset(Vec::new())
    );
    assert_eq!(
        sdf::write(&read, SdfWriteOptions::default()).unwrap().0,
        text
    );
    // A set with such a distribution round-trips.
    let mut empty = methane_set(2);
    empty
        .set_distribution(Some(
            Distribution::new(Members::Subset(Vec::new()), Partition::state_per_sample(0)).unwrap(),
        ))
        .unwrap();
    let written = sdf::write(&empty, SdfWriteOptions::default()).unwrap().0;
    let reread = read_back(&written).unwrap();
    assert_eq!(reread.distribution(), empty.distribution());
    assert_eq!(
        sdf::write(&reread, SdfWriteOptions::default()).unwrap().0,
        written
    );
    // With a member in the output, the others lie outside a distribution of
    // it.
    let some = [set.get(0).unwrap(), set.get(3).unwrap()];
    let read = read_back(&sdf::write(some, SdfWriteOptions::default()).unwrap().0).unwrap();
    assert_eq!(
        read.distribution().unwrap().members(),
        &Members::Subset(vec![1])
    );
}

#[test]
fn sdf_writes_weights_only_for_a_weighted_distribution() {
    for members in [Members::All, Members::Subset(vec![0, 2])] {
        let mut set = methane_set(3);
        let count = match &members {
            Members::Subset(members) => members.len(),
            _ => 3,
        };
        set.set_distribution(Some(
            Distribution::new(members.clone(), Partition::state_per_sample(count)).unwrap(),
        ))
        .unwrap();
        let (text, report) = sdf::write(&set, SdfWriteOptions::default()).unwrap();
        assert!(report.is_empty());
        assert!(!text.contains(sdf::WEIGHT_FIELD), "{members:?}");
        let read = read_back(&text).unwrap();
        assert_eq!(read.distribution(), set.distribution(), "{members:?}");
        assert_eq!(read.distribution().unwrap().weights(), None);
        assert_eq!(
            sdf::write(&read, SdfWriteOptions::default()).unwrap().0,
            text
        );
    }
}

#[test]
fn sdf_judges_each_record_by_its_own_sample() {
    let set = subset_set();
    // Interpreted records of members, which carry their distribution fields
    // as their own, written with a geometry outside the distribution.
    let members = interpret_records(
        &sdf::write(
            [set.get(1).unwrap(), set.get(3).unwrap()],
            SdfWriteOptions::default(),
        )
        .unwrap()
        .0,
    );
    let mixed: Vec<sdf::SdfRecordSource<'_>> = vec![
        (&members[0]).into(),
        set.get(0).unwrap().into(),
        (&members[1]).into(),
    ];
    let (text, report) = sdf::write(mixed, SdfWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    let read = read_back(&text).unwrap();
    let distribution = read.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::Subset(vec![0, 2]));
    assert_eq!(
        distribution.weights(),
        set.distribution().unwrap().weights()
    );
    assert_eq!(distribution.states(), set.distribution().unwrap().states());

    // A geometry outside the distribution of another set is marked all the
    // same.
    let other = subset_set();
    let records = [
        other.get(0).unwrap(),
        set.get(1).unwrap(),
        set.get(3).unwrap(),
    ];
    let read = read_back(&sdf::write(records, SdfWriteOptions::default()).unwrap().0).unwrap();
    assert_eq!(
        read.distribution().unwrap().members(),
        &Members::Subset(vec![1, 2])
    );
    // Members of two sets each carry their own set's state, so two
    // representatives of state 0 do not read as one distribution.
    let records = [set.get(1).unwrap(), other.get(1).unwrap()];
    assert_eq!(
        read_back(&sdf::write(records, SdfWriteOptions::default()).unwrap().0)
            .unwrap_err()
            .to_string(),
        "Cannot interpret the SDF input in record 2. The data field \"kekule.representative\" \
         marks a second representative of state 0."
    );
}

/// A two-atom model with each kind of state a Molfile and an SDF record
/// cannot carry, one at a time, with the field that names it and the
/// omission that reports it.
fn models_with_unwritable_state() -> Vec<(Structure, StructureField, WriteOmission)> {
    use kekule::properties::{PropertyKey, PropertyValue};
    use kekule::topology::{TopologyAtomIndex, TopologyBondIndex};
    use kekule::units::SQUARE_NANOMETER;

    let base = model("CO", &[[0.0; 3], [1.4, 0.0, 0.0]]);
    let key = || PropertyKey::new("note").unwrap();
    let state = |field| WriteOmission::StructureState {
        field,
        structures: 1,
    };
    let mut occupied = base.clone();
    occupied
        .set_occupancies(Some(vec![Some(0.5), None]))
        .unwrap();
    let mut b_factors = base.clone();
    b_factors
        .set_b_factors(Some(Quantity::new(vec![None, Some(0.2)], SQUARE_NANOMETER)))
        .unwrap();
    let mut atoms = base.clone();
    atoms
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
        .properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(3))
        .unwrap();
    vec![
        (
            occupied,
            StructureField::Occupancies,
            WriteOmission::Occupancies { structures: 1 },
        ),
        (
            b_factors,
            StructureField::BFactors,
            WriteOmission::BFactors { structures: 1 },
        ),
        (
            atoms,
            StructureField::AtomProperties,
            state(StructureField::AtomProperties),
        ),
        (
            bonds,
            StructureField::BondProperties,
            state(StructureField::BondProperties),
        ),
        (
            owner,
            StructureField::Properties,
            state(StructureField::Properties),
        ),
    ]
}

#[test]
fn molfile_and_sdf_refuse_or_report_sites_and_properties() {
    let omit = UnrepresentablePolicy::Omit;
    let plain = model("CO", &[[0.0; 3], [1.4, 0.0, 0.0]]);
    let plain_molfile = molfile::write(&plain, MolfileWriteOptions::default())
        .unwrap()
        .0;
    let plain_sdf = sdf::write([&plain], SdfWriteOptions::default()).unwrap().0;
    for (model, field, omission) in models_with_unwritable_state() {
        let properties = !matches!(
            field,
            StructureField::Occupancies | StructureField::BFactors
        );
        // Each field follows its own policy, and the other one does not
        // omit it.
        let (molfile_omit, molfile_other) = if properties {
            (
                MolfileWriteOptions::default().with_properties(omit),
                MolfileWriteOptions::default().with_set_metadata(omit),
            )
        } else {
            (
                MolfileWriteOptions::default().with_set_metadata(omit),
                MolfileWriteOptions::default().with_properties(omit),
            )
        };
        let (sdf_omit, sdf_other) = if properties {
            (
                SdfWriteOptions::default().with_properties(omit),
                SdfWriteOptions::default().with_set_metadata(omit),
            )
        } else {
            (
                SdfWriteOptions::default().with_set_metadata(omit),
                SdfWriteOptions::default().with_properties(omit),
            )
        };
        let verb = "have";
        for options in [MolfileWriteOptions::default(), molfile_other] {
            assert_eq!(
                described(molfile::write(&model, options).unwrap_err()),
                (
                    DiagnosticKind::Unsupported,
                    Some(FormatErrorDetail::UnrepresentableState(field)),
                    format!("The {field} of the structure {verb} no Molfile representation."),
                ),
                "{field:?}"
            );
        }
        for options in [SdfWriteOptions::default(), sdf_other] {
            let error = sdf::write([&plain, &model], options).unwrap_err();
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
        let (text, report) = molfile::write(&model, molfile_omit).unwrap();
        assert_eq!(text, plain_molfile, "{field:?}");
        assert_eq!(report.omissions(), [omission], "{field:?}");
        let (text, report) = sdf::write([&model], sdf_omit).unwrap();
        assert_eq!(text, plain_sdf, "{field:?}");
        assert_eq!(report.omissions(), [omission], "{field:?}");
    }

    // A geometry's sites and properties name the geometry.
    let (model, _, _) = models_with_unwritable_state().remove(0);
    let set = StructureSet::from_structures([model]).unwrap();
    assert_eq!(
        sdf::write(&set, SdfWriteOptions::default())
            .unwrap_err()
            .to_string(),
        "Cannot write the SDF output in record 1. The occupancies of the geometry have no SDF \
         representation."
    );
}

#[test]
fn sdf_and_molfile_refuse_or_report_set_and_state_properties() {
    use kekule::properties::{PropertyKey, PropertyValue};
    let omit = UnrepresentablePolicy::Omit;
    let key = || PropertyKey::new("note").unwrap();
    let state = |field, structures| WriteOmission::StructureState { field, structures };

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
            StructureField::SetProperties
        ))
    );
    let (text, report) =
        sdf::write(&set, SdfWriteOptions::default().with_properties(omit)).unwrap();
    assert_eq!(text, bare);
    assert_eq!(
        report.omissions(),
        [state(StructureField::SetProperties, 2)]
    );
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

    // A state's properties, on the record of its representative.
    let mut labelled = subset_set();
    let bare = sdf::write(&labelled, SdfWriteOptions::default()).unwrap().0;
    labelled
        .distribution_mut()
        .unwrap()
        .state_mut(1)
        .unwrap()
        .properties_mut()
        .insert(key(), PropertyValue::Int(2))
        .unwrap();
    let error = sdf::write(&labelled, SdfWriteOptions::default()).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Cannot write the SDF output in record 4. The state properties of the geometry have no \
         SDF representation."
    );
    let (text, report) = sdf::write(
        &labelled,
        SdfWriteOptions::default().with_set_metadata(omit),
    )
    .unwrap();
    assert_eq!(text, bare);
    assert_eq!(
        report.omissions(),
        [state(StructureField::StateProperties, 1)]
    );
    // A Molfile loses the membership whole.
    let (_, report) = molfile::write(
        molfile::MolfileSource::structure(&labelled.get(3).unwrap()),
        MolfileWriteOptions::default().with_set_metadata(omit),
    )
    .unwrap();
    assert_eq!(
        report.omissions(),
        [
            state(StructureField::Weight, 1),
            state(StructureField::Energy, 1),
            state(StructureField::Temperature, 1),
            state(StructureField::StateProperties, 1),
        ]
    );
}

#[test]
fn pdb_and_mmcif_refuse_or_report_structure_and_set_properties() {
    use kekule::properties::{PropertyKey, PropertyValue};
    use kekule::topology::TopologyAtomIndex;
    let omit = UnrepresentablePolicy::Omit;
    let key = || PropertyKey::new("note").unwrap();
    let state = |field, structures| WriteOmission::StructureState { field, structures };

    // A model's own properties.
    let plain = model("C", &[[0.0; 3]]);
    let mut noted = plain.clone();
    noted
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
    assert_eq!(report.omissions(), [state(StructureField::Properties, 1)]);
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
        .contains(&state(StructureField::Properties, 1)));

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
            pdb::PdbWriteIssue::UnrepresentableGeometry {
                geometry: 1,
                field: StructureField::SetProperties,
            }
        ))
    );
    assert_eq!(
        error.to_string(),
        "Cannot write the PDB output. The structure set properties have no PDB representation."
    );
    let expected = [
        state(StructureField::SetProperties, 2),
        state(StructureField::AtomProperties, 1),
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
        Some(MmcifWriteIssue::UnrepresentableGeometry {
            geometry: 2,
            field: StructureField::AtomProperties,
        })
    );
    // Occupancies are written, as PDB and mmCIF carry them.
    let (occupied, _, _) = models_with_unwritable_state().remove(0);
    assert!(pdb::write(&occupied, PdbWriteOptions::default()).is_ok());
}
