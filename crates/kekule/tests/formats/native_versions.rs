//! Frozen native files of every readable format version.
//!
//! `tests/fixtures/native/v<N>/` holds one file per object kind, written by
//! the version `N` writer and never regenerated. Each decodes to the object
//! its builder here gives, so a builder must keep building the object its
//! files hold. A new format version records its files and fuzz seeds with
//! `cargo test -p kekule --test formats native_versions::record -- --ignored`,
//! which writes only files that are missing.

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
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::smiles;
use kekule::structure::{
    Distribution, Ensemble, Forces, Geometry, Order, Partition, Positions, Samples, StateData,
    Structure, StructureSet, TimeAxis, Velocities,
};
use kekule::topology::{AtomSiteMetadata, InstanceAtomId, Topology, TopologyBuilder};
use kekule::units::{
    Quantity, ANGSTROM, DALTON, DEGREE, ELEMENTARY_CHARGE, KELVIN, KILOJOULE_PER_MOLE,
    KILOJOULE_PER_MOLE_PER_NANOMETER, NANOMETER, NANOMETER_PER_PICOSECOND, PICOSECOND,
    SQUARE_ANGSTROM,
};

use super::support::fixture;

/// One frozen file: its name, the first version that has its kind, the kind,
/// whether it stores perception, and the object it holds.
struct Fixture {
    name: &'static str,
    since: u32,
    kind: NativeKind,
    perception: bool,
    build: fn() -> Vec<u8>,
}

const FIXTURES: &[Fixture] = &[
    Fixture {
        name: "molecule",
        since: 7,
        kind: NativeKind::Molecule,
        perception: true,
        build: || write(&perceived_molecule(), true),
    },
    Fixture {
        name: "molecule_plain",
        since: 7,
        kind: NativeKind::Molecule,
        perception: false,
        build: || write(&plain_molecule(), false),
    },
    Fixture {
        name: "topology",
        since: 7,
        kind: NativeKind::Topology,
        perception: true,
        build: || write(&topology(), true),
    },
    Fixture {
        name: "structure",
        since: 7,
        kind: NativeKind::Structure,
        perception: true,
        build: || write(&structure(), true),
    },
    Fixture {
        name: "interpreted_structure",
        since: 7,
        kind: NativeKind::StructureInterpretation,
        perception: true,
        build: || {
            let interpretation = pdb_entry(false).interpret(Default::default()).unwrap();
            write(&interpretation, true)
        },
    },
    Fixture {
        name: "structure_set",
        since: 7,
        kind: NativeKind::StructureSet,
        perception: true,
        build: || write(&structure_set(), true),
    },
    Fixture {
        name: "interpreted_structure_set",
        since: 7,
        kind: NativeKind::StructureSetInterpretation,
        perception: true,
        build: || {
            let interpretation = pdb_entry(true)
                .interpret_structure_set(Default::default())
                .unwrap();
            write(&interpretation, true)
        },
    },
    Fixture {
        name: "structure_set_description",
        since: 7,
        kind: NativeKind::StructureSet,
        perception: true,
        build: || write(&structure_set_description(), true),
    },
    Fixture {
        name: "ensemble",
        since: 7,
        kind: NativeKind::Ensemble,
        perception: true,
        build: || write(&ensemble(), true),
    },
    Fixture {
        name: "ensemble_description",
        since: 7,
        kind: NativeKind::Ensemble,
        perception: true,
        build: || write(&ensemble_description(), true),
    },
    Fixture {
        name: "distribution",
        since: 8,
        kind: NativeKind::Distribution,
        perception: true,
        build: || write(&distribution(), true),
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

/// The files of `version`: every fixture whose kind the version has.
fn fixtures(version: u32) -> impl Iterator<Item = &'static Fixture> {
    FIXTURES
        .iter()
        .filter(move |fixture| fixture.since <= version)
}

#[test]
fn every_readable_version_has_a_file_of_every_kind() {
    for version in NATIVE_FORMAT_OLDEST..=NATIVE_FORMAT_VERSION {
        for fixture in fixtures(version) {
            let path = path(version, fixture.name);
            assert!(path.is_file(), "{} is missing", path.display());
        }
    }
}

#[test]
fn files_of_every_readable_version_decode_to_their_objects() {
    for version in NATIVE_FORMAT_OLDEST..=NATIVE_FORMAT_VERSION {
        for fixture in fixtures(version) {
            let name = format!("v{version}/{}", fixture.name);
            let bytes = std::fs::read(path(version, fixture.name)).unwrap();
            assert_eq!(bytes[8..12], version.to_le_bytes(), "{name}");
            let decoded = native::read(bytes.as_slice(), NativeReadOptions::default())
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            assert_eq!(decoded.kind(), fixture.kind, "{name}");
            // The current writer encodes the decoded object exactly as it
            // encodes the object built here.
            assert!(
                encode(&decoded, fixture.perception) == (fixture.build)(),
                "{name} does not decode to its object"
            );
        }
    }
}

#[test]
fn files_of_the_current_version_are_what_the_writer_writes() {
    // A layout change without a new version fails here.
    for fixture in fixtures(NATIVE_FORMAT_VERSION) {
        let bytes = std::fs::read(path(NATIVE_FORMAT_VERSION, fixture.name)).unwrap();
        let decoded = native::read(bytes.as_slice(), NativeReadOptions::default()).unwrap();
        assert!(
            encode(&decoded, fixture.perception) == bytes,
            "{} changed; a new layout needs a new version",
            fixture.name
        );
    }
}

/// Writes the current version's missing files, and the same files as seeds
/// of the `native` fuzz target.
#[test]
#[ignore = "records the files of a new format version"]
fn record() {
    let version = NATIVE_FORMAT_VERSION;
    let seeds = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fuzz/corpus/native");
    std::fs::create_dir_all(directory(version)).unwrap();
    for fixture in fixtures(version) {
        let bytes = (fixture.build)();
        for path in [
            path(version, fixture.name),
            seeds.join(format!("v{version}_{}.bin", fixture.name)),
        ] {
            if !path.exists() {
                std::fs::write(&path, &bytes).unwrap();
            }
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

/// Two states with energies and properties, log weights, and a temperature.
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
        .with_temperature(Quantity::new(298.15, KELVIN))
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
