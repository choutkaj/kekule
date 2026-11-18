//! Native files of every readable format version.
//!
//! `tests/fixtures/native/v<N>/` holds the files the version `N` writer
//! wrote. Files of versions before the current one are frozen, never
//! regenerated. The current writer encodes each decoded file to the bytes of
//! the same file of the current version, so an older version's files read
//! as the current ones do. A builder here gives the object of each file.
//! Stable builders, built by hand, must keep giving the current files;
//! builders that run perception or interpretation serve only to record
//! files. New content in a kind an older version has is a new fixture entry
//! with `since: N`, never an edit of its builder, whose older files are
//! frozen.
//!
//! The current version freezes when a release ships it. Until then, a
//! layout change or a new or edited fixture of the current version rewrites
//! every file of the current version and its fuzz seeds:
//! `cargo test -p kekule --test formats native_versions::record -- --ignored`.
//! A new version records its files the same way.

use std::path::PathBuf;
use std::sync::Arc;

use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{AtomId, BondId, Molecule};
use kekule::native::{
    self, DistributionDescription, EnsembleDescription, NativeKind, NativeObject,
    NativeReadOptions, NativeSource, NativeWriteOptions, StructureSetDescription,
    NATIVE_FORMAT_OLDEST, NATIVE_FORMAT_VERSION,
};
use kekule::pdb;
use kekule::perception::aromaticity::AromaticityModel;
use kekule::perception::rings::RingMembership;
use kekule::perception::valence::{ValenceModel, ValenceValidation};
use kekule::perception::Perception;
use kekule::properties::{PropertyColumn, PropertyKey, PropertyTable, PropertyValue};
use kekule::smiles;
use kekule::structure::{
    Distribution, Ensemble, Forces, Geometry, GeometryIndex, Order, Partition, Positions, Samples,
    Sampling, StateData, Structure, StructureSet, ThermodynamicState, TimeAxis, Velocities,
};
use kekule::topology::{AtomSiteMetadata, InstanceAtomId, Topology, TopologyBuilder};
use kekule::units::{
    Quantity, ANGSTROM, BAR, DALTON, DEGREE, ELEMENTARY_CHARGE, KELVIN, KILOJOULE_PER_MOLE,
    KILOJOULE_PER_MOLE_PER_NANOMETER, NANOMETER, NANOMETER_PER_PICOSECOND, PICOSECOND,
    SQUARE_ANGSTROM,
};

use super::support::fixture;

/// One fixture file: its name, the first version that has its content, the
/// kind, whether it stores perception, and how its object is built.
struct Fixture {
    name: &'static str,
    since: u32,
    kind: NativeKind,
    perception: bool,
    build: Build,
}

/// How a fixture's object is built.
enum Build {
    /// By hand: the builder keeps giving the current file.
    Stable(fn() -> Vec<u8>),
    /// By perception or interpretation, whose output may change: the builder
    /// only records new files.
    Recorded(fn() -> Vec<u8>),
}

const FIXTURES: &[Fixture] = &[
    Fixture {
        name: "molecule",
        since: 7,
        kind: NativeKind::Molecule,
        perception: true,
        build: Build::Stable(|| write(&perceived_molecule(), true)),
    },
    Fixture {
        name: "molecule_plain",
        since: 7,
        kind: NativeKind::Molecule,
        perception: false,
        build: Build::Stable(|| write(&plain_molecule(), false)),
    },
    Fixture {
        name: "molecule_perceived",
        since: 8,
        kind: NativeKind::Molecule,
        perception: true,
        build: Build::Recorded(|| write(&fully_perceived_molecule(), true)),
    },
    Fixture {
        name: "topology",
        since: 7,
        kind: NativeKind::Topology,
        perception: true,
        build: Build::Stable(|| write(&topology(), true)),
    },
    Fixture {
        name: "structure",
        since: 7,
        kind: NativeKind::Structure,
        perception: true,
        build: Build::Stable(|| write(&structure(), true)),
    },
    Fixture {
        name: "interpreted_structure",
        since: 7,
        kind: NativeKind::StructureInterpretation,
        perception: true,
        build: Build::Recorded(|| {
            let interpretation = pdb_entry(false).interpret(Default::default()).unwrap();
            write(&interpretation, true)
        }),
    },
    Fixture {
        name: "interpreted_structure_assemblies",
        since: 8,
        kind: NativeKind::StructureInterpretation,
        perception: true,
        build: Build::Recorded(|| {
            let interpretation = pdb_assemblies_and_alternates()
                .interpret(Default::default())
                .unwrap();
            let report = interpretation.report();
            assert!(!report.assemblies().is_empty());
            assert!(!report.alternate_locations().is_empty());
            write(&interpretation, true)
        }),
    },
    Fixture {
        name: "structure_set",
        since: 7,
        kind: NativeKind::StructureSet,
        perception: true,
        build: Build::Stable(|| write(&structure_set(), true)),
    },
    Fixture {
        name: "structure_set_steps",
        since: 8,
        kind: NativeKind::StructureSet,
        perception: true,
        build: Build::Stable(|| write(&structure_set_steps(), true)),
    },
    Fixture {
        name: "structure_set_columns",
        since: 9,
        kind: NativeKind::StructureSet,
        perception: true,
        build: Build::Stable(|| write(&structure_set_columns(), true)),
    },
    Fixture {
        name: "structure_set_sampling",
        since: 9,
        kind: NativeKind::StructureSet,
        perception: true,
        build: Build::Stable(|| write(&structure_set_sampling(), true)),
    },
    Fixture {
        name: "interpreted_structure_set",
        since: 7,
        kind: NativeKind::StructureSetInterpretation,
        perception: true,
        build: Build::Recorded(|| {
            let interpretation = pdb_entry(true)
                .interpret_structure_set(Default::default())
                .unwrap();
            write(&interpretation, true)
        }),
    },
    Fixture {
        name: "structure_set_description",
        since: 7,
        kind: NativeKind::StructureSet,
        perception: true,
        build: Build::Stable(|| write(&structure_set_description(), true)),
    },
    Fixture {
        name: "structure_set_description_columns",
        since: 9,
        kind: NativeKind::StructureSet,
        perception: true,
        build: Build::Stable(|| write(&structure_set_description_columns(), true)),
    },
    Fixture {
        name: "ensemble",
        since: 7,
        kind: NativeKind::Ensemble,
        perception: true,
        build: Build::Stable(|| write(&ensemble(), true)),
    },
    Fixture {
        name: "ensemble_description",
        since: 7,
        kind: NativeKind::Ensemble,
        perception: true,
        build: Build::Stable(|| write(&ensemble_description(), true)),
    },
    Fixture {
        name: "distribution",
        since: 8,
        kind: NativeKind::Distribution,
        perception: false,
        build: Build::Stable(|| write(&distribution(), false)),
    },
    Fixture {
        name: "distribution_boltzmann",
        since: 9,
        kind: NativeKind::Distribution,
        perception: false,
        build: Build::Stable(|| write(&boltzmann_distribution(), false)),
    },
];

fn write<'a>(source: impl Into<NativeSource<'a>>, perception: bool) -> Vec<u8> {
    native::write(
        source,
        NativeWriteOptions::default().with_perception(perception),
    )
    .unwrap()
}

/// Encodes a decoded object again with the current writer.
fn encode(object: &NativeObject, perception: bool) -> Vec<u8> {
    match object {
        NativeObject::Molecule(value) => write(value.as_ref(), perception),
        NativeObject::Topology(value) => write(value, perception),
        NativeObject::Structure(value) => write(value, perception),
        NativeObject::StructureInterpretation(value) => write(value.as_ref(), perception),
        NativeObject::StructureSet(value) => write(value, perception),
        NativeObject::StructureSetInterpretation(value) => write(value.as_ref(), perception),
        NativeObject::StructureSetDescription(value) => write(value.as_ref(), perception),
        NativeObject::Ensemble(value) => write(value.as_ref(), perception),
        NativeObject::EnsembleDescription(value) => write(value.as_ref(), perception),
        NativeObject::DistributionDescription(value) => write(value.as_ref(), perception),
        _ => unreachable!("every kind is listed"),
    }
}

fn directory(version: u32) -> PathBuf {
    fixture(&format!("native/v{version}"))
}

fn path(version: u32, name: &str) -> PathBuf {
    directory(version).join(format!("{name}.bin"))
}

fn read(version: u32, name: &str) -> Vec<u8> {
    let path = path(version, name);
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// Decodes the file of `fixture` in `version`, checking its header.
fn decode(version: u32, fixture: &Fixture) -> NativeObject {
    let name = format!("v{version}/{}", fixture.name);
    let bytes = read(version, fixture.name);
    assert_eq!(bytes[8..12], version.to_le_bytes(), "{name}");
    let decoded = native::read(bytes.as_slice(), NativeReadOptions::default())
        .unwrap_or_else(|error| panic!("{name}: {error}"));
    assert_eq!(decoded.kind(), fixture.kind, "{name}");
    decoded
}

/// The files of `version`: every fixture whose content the version has.
fn fixtures(version: u32) -> impl Iterator<Item = &'static Fixture> {
    FIXTURES
        .iter()
        .filter(move |fixture| fixture.since <= version)
}

#[test]
fn every_readable_version_has_its_files_and_no_other() {
    for version in NATIVE_FORMAT_OLDEST..=NATIVE_FORMAT_VERSION {
        for fixture in fixtures(version) {
            let path = path(version, fixture.name);
            assert!(path.is_file(), "{} is missing", path.display());
        }
        for entry in std::fs::read_dir(directory(version)).unwrap() {
            let file = entry.unwrap().file_name().into_string().unwrap();
            assert!(
                fixtures(version).any(|fixture| file == format!("{}.bin", fixture.name)),
                "v{version}/{file} is no fixture of version {version}"
            );
        }
    }
}

#[test]
fn files_exist_only_for_readable_versions() {
    for entry in std::fs::read_dir(fixture("native")).unwrap() {
        let name = entry.unwrap().file_name().into_string().unwrap();
        let version = name
            .strip_prefix('v')
            .and_then(|version| version.parse::<u32>().ok())
            .unwrap_or_else(|| panic!("native/{name} is no version directory"));
        assert!(
            (NATIVE_FORMAT_OLDEST..=NATIVE_FORMAT_VERSION).contains(&version),
            "native/{name} is outside the readable versions"
        );
    }
}

#[test]
fn files_of_every_readable_version_read_as_the_current_files() {
    for version in NATIVE_FORMAT_OLDEST..=NATIVE_FORMAT_VERSION {
        for fixture in fixtures(version) {
            let encoded = encode(&decode(version, fixture), fixture.perception);
            // At the current version, a layout change without recorded files
            // fails here.
            assert!(
                encoded == read(NATIVE_FORMAT_VERSION, fixture.name),
                "v{version}/{0} does not encode to v{NATIVE_FORMAT_VERSION}/{0}",
                fixture.name
            );
        }
    }
}

#[test]
fn stable_builders_build_the_current_files() {
    for fixture in fixtures(NATIVE_FORMAT_VERSION) {
        if let Build::Stable(build) = fixture.build {
            assert!(
                build() == read(NATIVE_FORMAT_VERSION, fixture.name),
                "the {} builder no longer builds its file (run native_versions::record)",
                fixture.name
            );
        }
    }
}

#[test]
fn sets_of_versions_before_9_have_no_columns() {
    for version in NATIVE_FORMAT_OLDEST..9 {
        let set = decode(version, named("structure_set"))
            .into_structure_set()
            .unwrap();
        assert_eq!(set.columns().len(), set.len(), "v{version}");
        assert!(!set.columns().has_data(), "v{version}");
        assert!(set.sampling().is_none(), "v{version}");
        let description = decode(version, named("structure_set_description"))
            .into_structure_set_description()
            .unwrap();
        assert_eq!(description.columns().len(), 4, "v{version}");
        assert!(!description.columns().has_data(), "v{version}");
    }
    let columns = decode(9, named("structure_set_columns"))
        .into_structure_set()
        .unwrap();
    assert_eq!(columns.columns(), structure_set_columns().columns());
    let sampled = decode(9, named("structure_set_sampling"))
        .into_structure_set()
        .unwrap();
    assert_eq!(sampled.sampling(), structure_set_sampling().sampling());
}

fn named(name: &str) -> &'static Fixture {
    FIXTURES
        .iter()
        .find(|fixture| fixture.name == name)
        .unwrap()
}

/// Writes every file of the current version, and the same files as seeds
/// of the `native` fuzz target. A file whose content the previous version
/// has is that version's file encoded by the current writer.
#[test]
#[ignore = "records the files of the current format version"]
fn record() {
    let version = NATIVE_FORMAT_VERSION;
    let previous = version - 1;
    let seeds = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fuzz/corpus/native");
    std::fs::create_dir_all(directory(version)).unwrap();
    for fixture in fixtures(version) {
        let bytes = if fixture.since <= previous {
            encode(&decode(previous, fixture), fixture.perception)
        } else {
            match fixture.build {
                Build::Stable(build) | Build::Recorded(build) => build(),
            }
        };
        for path in [
            path(version, fixture.name),
            seeds.join(format!("v{version}_{}.bin", fixture.name)),
        ] {
            std::fs::write(&path, &bytes).unwrap();
        }
    }
}

// The objects the files hold.

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

/// Phenol with valence, ring, and aromaticity perception installed as
/// stated, and properties at every molecule scope.
fn perceived_molecule() -> Molecule {
    let mut molecule = smiles::to_molecules("OC1=CC=CC=C1").unwrap().remove(0);
    let ring = (1..7).map(AtomId::new).collect::<Vec<_>>();
    let ring_bonds = (1..7).map(BondId::new).collect::<Vec<_>>();
    let hydrogens = [(0, 1), (1, 0), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1)]
        .map(|(atom, count)| (AtomId::new(atom), count));
    let mut flags = vec![true; 7];
    flags[0] = false;
    let perception = Perception::builder()
        .with_valence(
            ValenceModel::RdkitLike,
            ValenceValidation::Strict,
            hydrogens.to_vec(),
        )
        .unwrap()
        .with_rings(RingMembership::from_flags(flags.clone(), flags), None)
        .with_aromaticity(AromaticityModel::RdkitLike, ring, ring_bonds)
        .unwrap()
        .build();
    molecule.install_perception(perception).unwrap();
    let mut properties = molecule.properties_mut();
    properties
        .owner_mut()
        .insert(key("mass"), PropertyValue::real(94.11, DALTON).unwrap())
        .unwrap();
    properties
        .owner_mut()
        .insert(key("name"), PropertyValue::String("phenol".into()))
        .unwrap();
    properties
        .atoms_mut()
        .set_value(
            key("charge"),
            AtomId::new(0),
            Some(PropertyValue::real(-0.5, ELEMENTARY_CHARGE).unwrap()),
        )
        .unwrap();
    properties
        .bonds_mut()
        .set_value(
            key("rotatable"),
            BondId::new(0),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    molecule
}

/// Zwitterionic alanine with a stereocenter, an isotope, and an atom map,
/// stored without perception.
fn plain_molecule() -> Molecule {
    let mut molecule = smiles::to_molecules("[13CH3:5][C@H]([NH3+])C(=O)[O-]")
        .unwrap()
        .remove(0);
    molecule
        .properties_mut()
        .owner_mut()
        .insert(key("index"), PropertyValue::Int(-3))
        .unwrap();
    molecule
}

/// A phenylalanine-like zwitterion with two stereo groups and every
/// perception section: valence, rings with their basis, aromaticity,
/// conjugation, resonance, and CIP descriptors.
fn fully_perceived_molecule() -> Molecule {
    let mut molecule = smiles::to_molecules("N[C@@H](Cc1ccccc1)[C@@H](C)C(=O)[O-] |&1:1,o1:9|")
        .unwrap()
        .remove(0);
    molecule.perceive().unwrap();
    kekule::stereo::cip::assign_cip_descriptors(&mut molecule, Default::default()).unwrap();
    assert_eq!(molecule.stereo_groups().count(), 2);
    let perception = molecule.perception();
    assert!(perception.rings().and_then(|rings| rings.basis()).is_some());
    assert!(perception.conjugation().is_some());
    assert!(perception.resonance().is_some());
    assert!(perception.cip().is_some());
    molecule
}

/// Two waters of one definition and an acetate, in an interleaved atom order,
/// with a chain, residues, atom sites, explicit classes, and properties.
fn topology() -> Topology {
    let water = smiles::to_molecules("O").unwrap().remove(0);
    let acetate = smiles::to_molecules("CC(=O)[O-]").unwrap().remove(0);
    let mut builder = TopologyBuilder::new();
    let water_definition = builder.add_molecule_definition(water).unwrap();
    let first = builder.add_instance(water_definition).unwrap();
    let ligand = builder.add_molecule(acetate).unwrap();
    let second = builder.add_instance(water_definition).unwrap();
    builder
        .set_molecule_class(
            kekule::topology::MoleculeDefinitionId::new(1),
            kekule::topology::MoleculeClass::Other,
        )
        .unwrap();
    let atom = |instance, atom| InstanceAtomId::new(instance, AtomId::new(atom));
    let mut order = vec![atom(ligand, 0), atom(second, 0), atom(first, 0)];
    order.extend((1..4).map(|index| atom(ligand, index)));
    builder.set_atom_order(order).unwrap();
    let hierarchy = builder.hierarchy_mut();
    let chain = hierarchy.add_chain("B", Some("A".into())).unwrap();
    let residue = hierarchy
        .add_residue(chain, "ACT", Some(1), Some("101".into()), Some("A".into()))
        .unwrap();
    hierarchy
        .set_residue_class(residue, kekule::topology::ResidueClass::Other)
        .unwrap();
    for index in 0..4 {
        hierarchy
            .add_atom_site(
                residue,
                atom(ligand, index),
                AtomSiteMetadata::default()
                    .with_label_atom_id(Some(format!("C{index}")))
                    .with_type_symbol(Some("C".into())),
            )
            .unwrap();
    }
    let solvent = hierarchy.add_chain("W", None).unwrap();
    for (number, instance) in [(201, second), (200, first)] {
        let residue = hierarchy
            .add_residue(solvent, "HOH", None, Some(number.to_string()), None)
            .unwrap();
        hierarchy
            .add_atom_site(
                residue,
                atom(instance, 0),
                AtomSiteMetadata::default().with_label_atom_id(Some("O".into())),
            )
            .unwrap();
    }
    let mut topology = builder.build().unwrap();
    let mut properties = topology.properties_mut();
    properties
        .owner_mut()
        .insert(key("title"), PropertyValue::String("fixture".into()))
        .unwrap();
    properties
        .molecule_instances_mut()
        .set_value(key("tag"), first, Some(PropertyValue::Int(7)))
        .unwrap();
    properties
        .residues_mut()
        .set_value(
            key("score"),
            residue,
            Some(PropertyValue::real(2.0, ANGSTROM).unwrap()),
        )
        .unwrap();
    topology
}

fn water_topology() -> Arc<Topology> {
    Arc::new(Topology::from_molecules(smiles::to_molecules("O.O").unwrap()).unwrap())
}

fn positions(seed: f64) -> Positions {
    Positions::new(Quantity::new(
        vec![
            Point3::new(seed, -0.0, 1e-300),
            Point3::new(1.0 / 3.0, seed, -2.5),
        ],
        NANOMETER,
    ))
    .unwrap()
}

fn cell() -> PeriodicCell {
    PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(3.0, 0.0, 0.0),
                Vector3::new(0.1, 3.0, 0.0),
                Vector3::new(0.0, 0.0, 3.0),
            ],
            NANOMETER,
        ),
        [true, true, false],
    )
    .unwrap()
}

/// A water geometry with every per-geometry field: a cell, occupancies, B
/// factors, properties, velocities, and forces.
fn rich_geometry(seed: f64) -> Geometry {
    let mut geometry = Geometry::new(positions(seed)).with_cell(cell());
    geometry
        .set_occupancies(Some(vec![Some(0.5), None]))
        .unwrap();
    geometry
        .set_b_factors(Some(Quantity::new(vec![None, Some(-0.0)], SQUARE_ANGSTROM)))
        .unwrap();
    geometry
        .properties_mut()
        .owner_mut()
        .insert(key("replica"), PropertyValue::Int(4))
        .unwrap();
    geometry
        .set_velocities(Some(
            Velocities::from_vec(Quantity::new(
                vec![Vector3::new(1.0, -2.0, 3.0); 2],
                NANOMETER_PER_PICOSECOND,
            ))
            .unwrap(),
        ))
        .unwrap();
    geometry
        .set_forces(Some(
            Forces::from_vec(Quantity::new(
                vec![Vector3::new(-1.0, 0.5, 1e-300); 2],
                KILOJOULE_PER_MOLE_PER_NANOMETER,
            ))
            .unwrap(),
        ))
        .unwrap();
    geometry
}

fn structure() -> Structure {
    Structure::from_parts(water_topology(), rich_geometry(0.25)).unwrap()
}

/// Three water geometries, the first with every field, in a time order with
/// times and steps, and set properties.
fn structure_set() -> StructureSet {
    let geometries = vec![
        rich_geometry(0.5),
        Geometry::new(positions(1.0)),
        Geometry::new(positions(f64::MIN_POSITIVE)),
    ];
    let mut set = StructureSet::from_geometries(water_topology(), geometries).unwrap();
    set.properties_mut()
        .insert(key("method"), PropertyValue::String("REMD".into()))
        .unwrap();
    let times = Quantity::new(vec![0.0, 0.002, 0.004], PICOSECOND);
    set.set_order(Some(Order::Time(
        TimeAxis::new(Some(times), Some(vec![0, 1, u64::MAX])).unwrap(),
    )))
    .unwrap();
    set
}

/// Three water geometries in a step order.
fn structure_set_steps() -> StructureSet {
    let geometries = [0.75, -1.0, 2.0].map(|seed| Geometry::new(positions(seed)));
    let mut set = StructureSet::from_geometries(water_topology(), Vec::from(geometries)).unwrap();
    set.set_order(Some(Order::Step(vec![u64::MAX, 0, 7])))
        .unwrap();
    set
}

/// Four geometries kept elsewhere, in a coordinate order.
fn structure_set_description() -> StructureSetDescription {
    StructureSetDescription::new(water_topology(), 4)
        .with_properties({
            let mut properties = kekule::properties::OwnerProperties::new();
            properties
                .insert(key("source"), PropertyValue::String("scan.xtc".into()))
                .unwrap();
            properties
        })
        .with_order(Some(Order::Coordinate {
            unit: DEGREE,
            values: vec![-0.0, 30.0, 60.0, 90.0],
        }))
        .unwrap()
}

/// A column of every type, each with a missing cell, for three geometries.
fn columns_table() -> PropertyTable<GeometryIndex> {
    let mut set = StructureSet::from_geometries(
        water_topology(),
        [0.0, 1.0, 2.0].map(|seed| Geometry::new(positions(seed))),
    )
    .unwrap();
    let mut columns = set.columns_mut();
    for (name, column) in [
        (
            "bias",
            PropertyColumn::Real {
                unit: KILOJOULE_PER_MOLE,
                values: vec![Some(-0.0), None, Some(1e-300)],
            },
        ),
        (
            "replica",
            PropertyColumn::Int(vec![None, Some(i64::MIN), Some(3)]),
        ),
        (
            "accepted",
            PropertyColumn::Bool(vec![Some(true), None, Some(false)]),
        ),
        (
            "stage",
            PropertyColumn::String(vec![Some("nvt".into()), Some("npt".into()), None]),
        ),
    ] {
        columns.insert(key(name), column).unwrap();
    }
    set.columns().clone()
}

/// Three water geometries in a step order, with columns.
fn structure_set_columns() -> StructureSet {
    let mut set = structure_set_steps();
    let mut columns = set.columns_mut();
    for (key, column) in columns_table().iter() {
        columns.insert(key.clone(), column.clone()).unwrap();
    }
    set
}

/// The set with columns and a sampling record in the isothermal-isobaric
/// ensemble, which names columns without values.
fn structure_set_sampling() -> StructureSet {
    let mut set = structure_set_columns();
    let state = ThermodynamicState::isothermal_isobaric(
        Quantity::new(298.15, KELVIN),
        Quantity::new(1.01325, BAR),
    )
    .unwrap();
    let sampling = Sampling::new(state)
        .with_from_step(7)
        .with_method("GaMD dual boost")
        .with_potential_energy(key("potential_energy"))
        .unwrap()
        .with_bias([key("bias"), key("bias.dihedral")])
        .unwrap();
    set.set_sampling(Some(sampling)).unwrap();
    set
}

/// Three geometries kept elsewhere, in a step order, with columns.
fn structure_set_description_columns() -> StructureSetDescription {
    StructureSetDescription::new(water_topology(), 3)
        .with_order(Some(Order::Step(vec![10, 20, 30])))
        .unwrap()
        .with_columns(columns_table())
        .unwrap()
}

/// Two states with energies and properties, log weights, and a canonical
/// target.
fn clustered_distribution(samples: Samples, sample_count: usize) -> Distribution {
    let assignments = (0..sample_count).map(|sample| sample % 2).collect();
    let mut open = StateData::new()
        .with_energy(Quantity::new(-12.5, KILOJOULE_PER_MOLE))
        .unwrap();
    open.properties_mut()
        .insert(key("cluster"), PropertyValue::String("open".into()))
        .unwrap();
    let log_weights = (0..sample_count)
        .map(|sample| -0.5 * sample as f64)
        .collect();
    Distribution::new(samples, Partition::new(assignments, vec![0, 1]).unwrap())
        .unwrap()
        .with_log_weights(log_weights)
        .unwrap()
        .with_states(vec![open, StateData::new()])
        .unwrap()
        .with_target(ThermodynamicState::canonical(Quantity::new(298.15, KELVIN)).unwrap())
        .unwrap()
}

fn ensemble() -> Ensemble {
    Ensemble::new(structure_set(), clustered_distribution(Samples::All, 3)).unwrap()
}

/// The last three of four geometries kept elsewhere, as two states.
fn ensemble_description() -> EnsembleDescription {
    EnsembleDescription::new(
        structure_set_description(),
        clustered_distribution(Samples::Subset(vec![1, 2, 3]), 3),
    )
    .unwrap()
}

/// Every other geometry of five kept elsewhere, as two states.
fn distribution() -> DistributionDescription {
    DistributionDescription::new(5, clustered_distribution(Samples::Subset(vec![0, 2, 4]), 3))
        .unwrap()
}

/// Four geometries as two states with Boltzmann populations at an
/// isothermal-isobaric target.
fn boltzmann_distribution() -> DistributionDescription {
    let target = ThermodynamicState::isothermal_isobaric(
        Quantity::new(310.0, KELVIN),
        Quantity::new(1.013_25, BAR),
    )
    .unwrap();
    let states = [0.0, 2.5].map(|energy| {
        StateData::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap()
    });
    let geometries = [0.5, 1.0, 1.5, 2.0].map(|seed| Geometry::new(positions(seed)));
    let set = StructureSet::from_geometries(water_topology(), Vec::from(geometries)).unwrap();
    let distribution = Distribution::new(
        Samples::All,
        Partition::new(vec![0, 1, 0, 1], vec![0, 1]).unwrap(),
    )
    .unwrap()
    .with_states(Vec::from(states))
    .unwrap();
    let mut ensemble = Ensemble::new(set, distribution).unwrap();
    ensemble
        .distribution_mut()
        .assign_boltzmann_populations(target)
        .unwrap();
    DistributionDescription::new(4, ensemble.into_parts().1).unwrap()
}

/// A peptide of two glycines and a water, in one model or in two.
fn pdb_entry(two_models: bool) -> pdb::PdbDocument {
    let atoms = [
        ("ATOM", " N  ", "GLY", 1, -1.195, 0.661, "N"),
        ("ATOM", " CA ", "GLY", 1, 0.0, 1.493, "C"),
        ("ATOM", " C  ", "GLY", 1, 1.29, 0.698, "C"),
        ("ATOM", " O  ", "GLY", 1, 1.29, -0.531, "O"),
        ("ATOM", " N  ", "GLY", 2, 2.43, 1.38, "N"),
        ("ATOM", " CA ", "GLY", 2, 3.72, 0.7, "C"),
        ("ATOM", " C  ", "GLY", 2, 4.9, 1.65, "C"),
        ("ATOM", " O  ", "GLY", 2, 4.7, 2.87, "O"),
        ("ATOM", " OXT", "GLY", 2, 6.05, 1.18, "O"),
        ("HETATM", " O  ", "HOH", 101, 0.0, 5.0, "O"),
    ];
    let mut text = String::from("HEADER    TEST\n");
    let models: &[f64] = if two_models { &[0.0, 0.1] } else { &[0.0] };
    for (model, shift) in models.iter().enumerate() {
        if two_models {
            text.push_str(&format!("MODEL     {:>4}\n", model + 1));
        }
        for (serial, (record, name, residue, number, x, y, element)) in atoms.iter().enumerate() {
            text.push_str(&format!(
                "{record:<6}{:>5} {name} {residue} A{number:>4}    {:>8.3}{y:>8.3}{:>8.3}  1.00 10.00          {element:>2}\n",
                serial + 1,
                x + shift,
                0.0,
            ));
        }
        if two_models {
            text.push_str("ENDMDL\n");
        }
    }
    text.push_str("END\n");
    pdb::parse_str(&text, Default::default()).unwrap()
}

/// Two glycines in chains A and B, one with two locations of its alpha
/// carbon, and a biomolecule that copies chain B.
fn pdb_assemblies_and_alternates() -> pdb::PdbDocument {
    let mut text = String::from(
        "\
HEADER    TEST
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350 IN ADDITION APPLY THE FOLLOWING TO CHAINS: B
REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000       10.00000
REMARK 350   BIOMT2   2  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
",
    );
    let atoms = [
        (" N  ", ' ', 'A', -1.195, 0.661, 1.0, "N"),
        (" CA ", 'A', 'A', 0.0, 1.493, 0.6, "C"),
        (" CA ", 'B', 'A', 0.1, 1.593, 0.4, "C"),
        (" C  ", ' ', 'A', 1.29, 0.698, 1.0, "C"),
        (" O  ", ' ', 'A', 1.29, -0.531, 1.0, "O"),
        (" N  ", ' ', 'B', 3.805, 0.661, 1.0, "N"),
        (" CA ", ' ', 'B', 5.0, 1.493, 1.0, "C"),
        (" C  ", ' ', 'B', 6.29, 0.698, 1.0, "C"),
        (" O  ", ' ', 'B', 6.29, -0.531, 1.0, "O"),
    ];
    for (serial, (name, alternate, chain, x, y, occupancy, element)) in atoms.iter().enumerate() {
        text.push_str(&format!(
            "ATOM  {:>5} {name}{alternate}GLY {chain}   1    {x:>8.3}{y:>8.3}{:>8.3}{occupancy:>6.2} 10.00          {element:>2}\n",
            serial + 1,
            0.0,
        ));
    }
    text.push_str("END\n");
    pdb::parse_str(&text, Default::default()).unwrap()
}
