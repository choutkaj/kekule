//! Native serialization: exact round trips, canonical bytes, and rejection of
//! malformed input without panics.

use std::num::NonZeroU32;
use std::sync::Arc;

use kekule::diagnostics::{DiagnosticKind, Format, FormatErrorDetail, Resource};
use kekule::macromolecular::{AltLocPolicy, StructureInterpretOptions};
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::mmcif::MmcifWriteOptions;
use kekule::molecule::{AtomId, Molecule};
use kekule::native::{
    self, NativeKind, NativeObject, NativeReadOptions, NativeSource, NativeWriteOptions,
    StructureSetDescription, NATIVE_FORMAT_VERSION,
};
use kekule::pdb::PdbWriteOptions;
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{
    Distribution, Forces, Geometry, Members, Order, Partition, Place, Positions, State, Structure,
    StructureSet, TimeAxis, Velocities,
};
use kekule::topology::{
    AtomSiteMetadata, InstanceAtomId, MoleculeClass, ResidueClass, Topology, TopologyBuilder,
};
use kekule::units::{
    Quantity, ANGSTROM, DALTON, DEGREE, ELEMENTARY_CHARGE, KILOJOULE_PER_MOLE, NANOMETER,
    PICOSECOND,
};
use kekule::{mmcif, pdb, smiles};

use super::support::{fixture, hierarchy_record};

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

/// Encodes `source`, decodes it, and checks that the decoded object encodes
/// to the same bytes: the encoding is canonical and loses nothing it stores.
fn round_trip<'a>(source: impl Into<NativeSource<'a>>) -> NativeObject {
    let bytes = native::write(source, NativeWriteOptions::default()).unwrap();
    let decoded = native::read(bytes.as_slice(), NativeReadOptions::default()).unwrap();
    let again = match &decoded {
        NativeObject::Molecule(value) => {
            native::write(value.as_ref(), NativeWriteOptions::default())
        }
        NativeObject::Topology(value) => native::write(value, NativeWriteOptions::default()),
        NativeObject::Structure(value) => native::write(value, NativeWriteOptions::default()),
        NativeObject::StructureSet(value) => native::write(value, NativeWriteOptions::default()),
        NativeObject::StructureSetDescription(value) => {
            native::write(value.as_ref(), NativeWriteOptions::default())
        }
        NativeObject::StructureInterpretation(value) => {
            native::write(value.as_ref(), NativeWriteOptions::default())
        }
        NativeObject::SampleSetInterpretation(value) => {
            native::write(value.as_ref(), NativeWriteOptions::default())
        }
        _ => unreachable!("every kind is listed"),
    }
    .unwrap();
    assert!(
        again == bytes,
        "re-encoding a decoded object changed its bytes"
    );
    decoded
}

fn assert_same_molecule(decoded: &Molecule, source: &Molecule) {
    assert_eq!(decoded, source);
    assert_eq!(decoded.perception(), source.perception());
    assert_eq!(decoded.properties(), source.properties());
    for atom in source.atom_ids() {
        let order = |molecule: &Molecule| {
            molecule
                .incident_bonds(atom)
                .unwrap()
                .map(|(id, _)| id)
                .collect::<Vec<_>>()
        };
        assert_eq!(order(decoded), order(source), "adjacency of {atom}");
    }
}

fn assert_same_topology(decoded: &Topology, source: &Topology) {
    assert!(decoded.same_layout(source));
    assert!(
        !decoded.shares_layout(source),
        "decoding publishes a new layout"
    );
    assert_eq!(decoded.atom_ids(), source.atom_ids());
    assert_eq!(decoded.bond_ids(), source.bond_ids());
    assert_eq!(hierarchy_record(decoded), hierarchy_record(source));
    assert_eq!(decoded.properties(), source.properties());
    for (left, right) in decoded.definitions().zip(source.definitions()) {
        assert_same_molecule(left.molecule(), right.molecule());
    }
    for (left, right) in decoded.instances().zip(source.instances()) {
        assert_eq!(left.class_assignment(), right.class_assignment());
    }
    for (left, right) in decoded.residues().zip(source.residues()) {
        assert_eq!(left.class_assignment(), right.class_assignment());
        assert_eq!(left.label_comp_id(), right.label_comp_id());
        assert_eq!(left.author_comp_id(), right.author_comp_id());
    }
}

fn bits(points: &[Point3]) -> Vec<[u64; 3]> {
    points
        .iter()
        .map(|point| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()])
        .collect()
}

fn scalar_bits(values: Option<&[Option<f64>]>) -> Option<Vec<Option<u64>>> {
    values.map(|values| values.iter().map(|value| value.map(f64::to_bits)).collect())
}

fn assert_same_structure(decoded: &Structure, source: &Structure) {
    assert_same_topology(decoded.topology(), source.topology());
    assert_eq!(
        bits(decoded.positions().values().value()),
        bits(source.positions().values().value())
    );
    assert_eq!(decoded.cell(), source.cell());
    assert_eq!(
        scalar_bits(decoded.occupancies()),
        scalar_bits(source.occupancies())
    );
    assert_eq!(
        scalar_bits(decoded.b_factors().map(|values| *values.value())),
        scalar_bits(source.b_factors().map(|values| *values.value()))
    );
    assert_eq!(decoded.properties(), source.properties());
}

fn perceived(input: &str) -> Molecule {
    let mut molecule = smiles::to_molecules(input).unwrap().remove(0);
    molecule.perceive().unwrap();
    molecule
}

#[test]
fn molecules_round_trip_represented_chemistry_perception_and_properties() {
    let inputs = [
        "C[C@H](N)C(=O)O",
        "F/C=C\\Cl",
        "[13CH3][N+](C)(C)C",
        "[CH3:7]O[NH3+]",
        "C[C@H](F)[C@@H](Cl)Br |o1:1,&2:3|",
        "c1ccc2ccccc2c1",
        "c1ccc2[nH]ccc2c1",
        "[O-][N+](=O)c1ccccc1",
        "[Fe+2]",
    ];
    for input in inputs {
        let mut molecule = perceived(input);
        let mut properties = molecule.properties_mut();
        properties
            .owner_mut()
            .insert(key("name"), PropertyValue::String(input.into()))
            .unwrap();
        properties
            .owner_mut()
            .insert(key("mass"), PropertyValue::real(-0.0, DALTON).unwrap())
            .unwrap();
        properties
            .atoms_mut()
            .set_value(
                key("charge"),
                AtomId::new(0),
                Some(PropertyValue::real(0.25, ELEMENTARY_CHARGE).unwrap()),
            )
            .unwrap();
        let decoded = round_trip(&molecule).into_molecule().unwrap();
        assert_same_molecule(&decoded, &molecule);
    }
}

#[test]
fn cip_resonance_and_fixed_hydrogens_survive() {
    let mut molecule = perceived("N[C@@H](C)C(=O)[O-]");
    kekule::stereo::cip::assign_cip_descriptors(&mut molecule, Default::default()).unwrap();
    assert!(molecule.perception().cip().is_some());
    let resonance = molecule.perception().resonance().unwrap();
    assert_eq!(resonance.bonds().len(), 2);
    let decoded = round_trip(&molecule).into_molecule().unwrap();
    assert_same_molecule(&decoded, &molecule);
}

#[test]
fn adjacency_order_is_part_of_the_round_trip() {
    let molecule = smiles::to_molecules("CCC").unwrap().pop().unwrap();
    let atoms = molecule.atom_ids().collect::<Vec<_>>();
    let bond = molecule.bond_ids().next().unwrap();
    let mut editor = molecule.to_editor();
    editor.set_bond_endpoints(bond, atoms[0], atoms[2]).unwrap();
    editor.set_bond_endpoints(bond, atoms[0], atoms[1]).unwrap();
    let rewired = editor.finish().unwrap();
    let decoded = round_trip(&rewired).into_molecule().unwrap();
    assert_same_molecule(&decoded, &rewired);
}

/// A topology with reused definitions, an interleaved dense order, a
/// hierarchy with author and label identifiers, explicit and inferred classes,
/// and annotations at every topology scope.
fn rich_topology() -> Topology {
    let water = smiles::to_molecules("O").unwrap().remove(0);
    let ligand = perceived("OCC(=O)[O-]");
    let mut builder = TopologyBuilder::new();
    let water_definition = builder.add_molecule_definition(water).unwrap();
    let first = builder.add_instance(water_definition).unwrap();
    let ligand_instance = builder.add_molecule(ligand).unwrap();
    let second = builder.add_instance(water_definition).unwrap();
    builder
        .set_molecule_class(
            kekule::topology::MoleculeDefinitionId::new(1),
            MoleculeClass::Other,
        )
        .unwrap();
    let atom = |instance, atom| InstanceAtomId::new(instance, AtomId::new(atom));
    let mut order = vec![atom(ligand_instance, 0), atom(second, 0), atom(first, 0)];
    order.extend((1..5).map(|index| atom(ligand_instance, index)));
    builder.set_atom_order(order).unwrap();
    let hierarchy = builder.hierarchy_mut();
    let chain = hierarchy.add_chain("B", Some("A".into())).unwrap();
    let solvent = hierarchy.add_chain("C", None).unwrap();
    let residue = hierarchy
        .add_residue(chain, "GLC", Some(1), Some("101".into()), Some("A".into()))
        .unwrap();
    hierarchy
        .set_residue_component_ids(residue, Some("GLC".into()), Some("LIG".into()))
        .unwrap();
    hierarchy
        .set_residue_class(residue, ResidueClass::Carbohydrate)
        .unwrap();
    for index in 0..5 {
        hierarchy
            .add_atom_site(
                residue,
                atom(ligand_instance, index),
                AtomSiteMetadata::default()
                    .with_label_atom_id(Some(format!("C{index}")))
                    .with_auth_atom_id(Some(format!("X{index}")))
                    .with_type_symbol(Some("C".into())),
            )
            .unwrap();
    }
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
        .insert(key("title"), PropertyValue::String("rich".into()))
        .unwrap();
    properties
        .molecule_instances_mut()
        .set_value(key("tag"), first, Some(PropertyValue::Int(-3)))
        .unwrap();
    properties
        .atoms_mut()
        .insert(
            key("partial_charge"),
            PropertyColumn::Real {
                unit: ELEMENTARY_CHARGE,
                values: vec![
                    Some(0.1),
                    None,
                    Some(-0.8),
                    None,
                    Some(0.0),
                    None,
                    Some(1.5),
                ],
            },
        )
        .unwrap();
    properties
        .bonds_mut()
        .set_value(
            key("rotatable"),
            kekule::topology::TopologyBondIndex::new(1),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    properties
        .chains_mut()
        .set_value(
            key("note"),
            chain,
            Some(PropertyValue::String("ligand".into())),
        )
        .unwrap();
    properties
        .residues_mut()
        .set_value(
            key("score"),
            residue,
            Some(PropertyValue::real(2.0, ANGSTROM).unwrap()),
        )
        .unwrap();
    properties
        .atom_sites_mut()
        .set_value(
            key("seen"),
            kekule::topology::AtomSiteId::new(2),
            Some(PropertyValue::Bool(false)),
        )
        .unwrap();
    topology
}

#[test]
fn topologies_keep_ids_order_hierarchy_classes_and_properties() {
    let topology = rich_topology();
    let decoded = round_trip(&topology).into_topology().unwrap();
    assert_same_topology(&decoded, &topology);
    let explicit = |topology: &Topology| {
        topology
            .instances()
            .map(|instance| instance.class_assignment().is_explicit())
            .collect::<Vec<_>>()
    };
    assert_eq!(explicit(&decoded), [false, true, false]);
    let unit = decoded
        .properties()
        .residues()
        .get(&key("score"))
        .map(|column| match column {
            PropertyColumn::Real { unit, .. } => unit.symbol(),
            _ => None,
        });
    assert_eq!(unit, Some(Some("A")), "named unit symbols are restored");
}

#[test]
fn macromolecular_models_round_trip_bitwise() {
    for name in ["mmcif/1AKE.cif", "mmcif/1CRN.cif"] {
        let text = std::fs::read_to_string(fixture(name)).unwrap();
        let mut model = mmcif::parse_str(&text, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap()
            .into_structure();
        model.perceive().unwrap();
        model
            .properties_mut()
            .owner_mut()
            .insert(key("source"), PropertyValue::String(name.into()))
            .unwrap();
        let decoded = round_trip(&model).into_structure().unwrap();
        assert_same_structure(&decoded, &model);
    }
}

#[test]
fn interpreted_models_keep_their_reports() {
    let mut unobserved = 0;
    for name in ["mmcif/9UOH.cif", "mmcif/1IRO.cif", "mmcif/1AKE.cif"] {
        let text = std::fs::read_to_string(fixture(name)).unwrap();
        let interpretation = mmcif::parse_str(&text, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap();
        let report = interpretation.report();
        assert!(!report.instances().is_empty(), "{name}");
        assert!(!report.polymer_entities().is_empty(), "{name}");
        let coordinated = !name.contains("1AKE");
        assert_eq!(!report.coordination().is_empty(), coordinated, "{name}");
        unobserved += report
            .polymer_entities()
            .iter()
            .flat_map(|entity| entity.instances())
            .flat_map(|instance| instance.positions())
            .filter(|position| position.residue().is_none())
            .count();

        let decoded = round_trip(&interpretation).into_interpretation().unwrap();
        assert_eq!(decoded.report(), report, "{name}");
        assert_same_structure(decoded.structure(), interpretation.structure());
        // The mmCIF and PDB writers state the same records, metal
        // coordination and entities included, from the decoded report.
        let written = |source| {
            mmcif::write([source], MmcifWriteOptions::default())
                .unwrap()
                .0
        };
        assert_eq!(written(&decoded), written(&interpretation), "{name}");
        // 9UOH has a five-character component, which PDB cannot hold.
        if !name.contains("9UOH") {
            let written = |source| pdb::write(source, PdbWriteOptions::default()).unwrap().0;
            assert_eq!(written(&decoded), written(&interpretation), "{name}");
            assert_eq!(
                written(&decoded).contains("\nLINK  "),
                coordinated,
                "{name}"
            );
        }

        // A model and its report held apart write the same bytes.
        let (model, report) = interpretation.clone().into_parts();
        let apart = native::write(
            NativeSource::interpreted_model(&model, &report),
            NativeWriteOptions::default(),
        )
        .unwrap();
        assert!(apart == native::write(&interpretation, NativeWriteOptions::default()).unwrap());
        // A reader that wants the model alone drops the report.
        let model = native::read(apart.as_slice(), NativeReadOptions::default())
            .unwrap()
            .into_structure()
            .unwrap();
        assert_same_structure(&model, interpretation.structure());
    }
    // The canonical sequences show residues the models lack.
    assert!(unobserved > 0);
}

#[test]
fn interpreted_conformations_keep_one_report_per_sample() {
    let text = std::fs::read_to_string(fixture("mmcif/1AKE.cif")).unwrap();
    let document = mmcif::parse_str(&text, Default::default()).unwrap();
    let selections = ["A", "B"].map(|label| {
        StructureInterpretOptions::default()
            .with_altloc_policy(AltLocPolicy::SelectLabel(label.into()))
    });
    let interpretation = document.interpret_conformations(&selections).unwrap();
    assert_eq!(interpretation.reports().len(), 2);
    assert!(interpretation
        .reports()
        .iter()
        .all(|report| !report.alternate_locations().is_empty()));
    assert_ne!(interpretation.reports()[0], interpretation.reports()[1]);
    // The conformations are one state of equally weighted geometries.
    let distribution = interpretation.structure_set().distribution().unwrap();
    assert_eq!(distribution.state_count(), 1);
    assert_eq!(distribution.weights(), None);
    assert_eq!(
        crate::native_fingerprint(&interpretation),
        0x09ae_12e3_3d1f_f605
    );

    let decoded = round_trip(&interpretation)
        .into_structure_set_interpretation()
        .unwrap();
    assert_eq!(decoded.reports(), interpretation.reports());
    assert_same_structure_set(decoded.structure_set(), interpretation.structure_set());
    let written = |source| {
        mmcif::write([source], MmcifWriteOptions::default())
            .unwrap()
            .0
    };
    assert_eq!(written(&decoded), written(&interpretation));
    assert_eq!(
        crate::fingerprint(written(&interpretation)),
        0x5812_6ef2_c2b1_32c0
    );

    let bytes = native::write(&interpretation, NativeWriteOptions::default()).unwrap();
    let set = native::read(bytes.as_slice(), NativeReadOptions::default())
        .unwrap()
        .into_structure_set()
        .unwrap();
    assert_same_structure_set(&set, interpretation.structure_set());
}

#[test]
fn reports_are_stored_only_when_given() {
    let text = std::fs::read_to_string(fixture("mmcif/1CRN.cif")).unwrap();
    let interpretation = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let plain = native::write(interpretation.structure(), NativeWriteOptions::default()).unwrap();
    let interpreted = native::write(&interpretation, NativeWriteOptions::default()).unwrap();
    // Both are format version 5: a model with its report is kind 7.
    assert_eq!(NATIVE_FORMAT_VERSION, 5);
    for bytes in [&plain, &interpreted] {
        assert_eq!(bytes[8..12], NATIVE_FORMAT_VERSION.to_le_bytes());
    }
    assert_eq!((plain[12], interpreted[12]), (3, 7));
    assert!(interpreted.len() > plain.len());

    let error = native::read(plain.as_slice(), NativeReadOptions::default())
        .unwrap()
        .into_interpretation()
        .unwrap_err();
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::NativeKind {
            expected: NativeKind::StructureInterpretation,
            found: NativeKind::Structure,
        })
    );
    assert_eq!(
        error.to_string(),
        "Cannot interpret the native input. The input holds a model, not an interpreted model."
    );

    // This reader rejects the retired kinds 4, 5, 6, and 8 and the unused
    // kind 11 as unknown.
    for code in [4, 5, 6, 8, 11] {
        let mut unknown = interpreted.clone();
        unknown[12] = code;
        let error = native::read(unknown.as_slice(), NativeReadOptions::default()).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::InvalidValue);
        assert!(
            error
                .to_string()
                .contains(&format!("The object kind {code} is unknown.")),
            "{error}"
        );
    }
}

fn assert_same_structure_set(decoded: &StructureSet, source: &StructureSet) {
    assert_same_topology(decoded.topology(), source.topology());
    assert_eq!(decoded.properties(), source.properties());
    assert_eq!(decoded.len(), source.len());
    for (left, right) in decoded.iter().zip(source.iter()) {
        assert_eq!(left.payload(), right.payload());
        assert_eq!(
            bits(left.payload().positions().values().value()),
            bits(right.payload().positions().values().value())
        );
    }
    assert_eq!(decoded.order(), source.order());
    assert_eq!(decoded.distribution(), source.distribution());
}

/// Three water geometries with every per-geometry field: a cell, occupancies, B
/// factors, properties, velocities, and forces on the first, none on the
/// second, and velocities alone on the third.
fn rich_samples() -> Vec<Geometry> {
    let mut first = Geometry::new(positions(0.5)).with_cell(cell());
    first.set_occupancies(Some(vec![Some(0.5), None])).unwrap();
    first
        .set_b_factors(Some(Quantity::new(
            vec![None, Some(-0.0)],
            kekule::units::SQUARE_ANGSTROM,
        )))
        .unwrap();
    first
        .properties_mut()
        .owner_mut()
        .insert(key("replica"), PropertyValue::Int(4))
        .unwrap();
    let velocities = || {
        Velocities::from_vec(Quantity::new(
            vec![Vector3::new(1.0, -2.0, 3.0); 2],
            kekule::units::NANOMETER_PER_PICOSECOND,
        ))
        .unwrap()
    };
    first.set_velocities(Some(velocities())).unwrap();
    first
        .set_forces(Some(
            Forces::from_vec(Quantity::new(
                vec![Vector3::new(-1.0, 0.5, 1e-300); 2],
                kekule::units::KILOJOULE_PER_MOLE_PER_NANOMETER,
            ))
            .unwrap(),
        ))
        .unwrap();
    let mut third = Geometry::new(positions(f64::MIN_POSITIVE));
    third.set_velocities(Some(velocities())).unwrap();
    vec![first, Geometry::new(positions(1.0)), third]
}

fn rich_set() -> StructureSet {
    let mut set = StructureSet::from_geometries(water_topology(), rich_samples()).unwrap();
    set.properties_mut()
        .insert(key("method"), PropertyValue::String("REMD".into()))
        .unwrap();
    set
}

#[test]
fn structure_sets_keep_every_order_kind_and_sample_field() {
    let times = |values: Vec<f64>| Some(Quantity::new(values, PICOSECOND));
    let orders = [
        None,
        Some(Order::Time(
            TimeAxis::new(times(vec![0.0, 0.002, -0.0]), Some(vec![0, 1, u64::MAX])).unwrap(),
        )),
        Some(Order::Time(
            TimeAxis::new(times(vec![1.0, 2.0, 3.0]), None).unwrap(),
        )),
        Some(Order::Time(
            TimeAxis::new(None, Some(vec![5, 10, 15])).unwrap(),
        )),
        Some(Order::Time(TimeAxis::new(None, None).unwrap())),
        Some(Order::Step(vec![u64::MAX, 0, 7])),
        Some(Order::Coordinate {
            unit: DEGREE,
            values: vec![-0.0, 1e-300, -90.5],
        }),
    ];
    for order in orders {
        let mut set = rich_set();
        set.set_order(order.clone()).unwrap();
        let decoded = round_trip(&set).into_structure_set().unwrap();
        assert_same_structure_set(&decoded, &set);
        assert_eq!(decoded.order(), order.as_ref());
    }
    // A coordinate keeps its sign bit and unit.
    let mut scan = rich_set();
    scan.set_order(Some(Order::Coordinate {
        unit: DEGREE,
        values: vec![-0.0, 1.0, 2.0],
    }))
    .unwrap();
    let decoded = round_trip(&scan).into_structure_set().unwrap();
    let Some(Order::Coordinate { unit, values }) = decoded.order() else {
        panic!("a coordinate order");
    };
    assert_eq!((*unit, values[0].to_bits()), (DEGREE, (-0.0_f64).to_bits()));

    // An empty set keeps its order.
    let mut empty = StructureSet::new(water_topology());
    empty
        .set_order(Some(Order::Coordinate {
            unit: ANGSTROM,
            values: Vec::new(),
        }))
        .unwrap();
    let decoded = round_trip(&empty).into_structure_set().unwrap();
    assert_same_structure_set(&decoded, &empty);
}

#[test]
fn structure_sets_keep_their_distribution() {
    let distributions = [
        // Every geometry, one state, without weights.
        Distribution::new(Members::All, Partition::one_state(3)).unwrap(),
        // Every geometry, weighted, a state each with energies and a temperature.
        Distribution::new(Members::All, Partition::state_per_sample(3))
            .unwrap()
            .with_weights(vec![0.25, 3.0, 1.0e-300])
            .unwrap()
            .with_states(vec![
                State::new()
                    .with_energy(Quantity::new(-12.5, KILOJOULE_PER_MOLE))
                    .unwrap(),
                State::new(),
                State::new()
                    .with_energy(Quantity::new(-0.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
            ])
            .unwrap()
            .with_temperature(Quantity::new(298.15, kekule::units::KELVIN))
            .unwrap(),
        // A subset: geometries 0 and 2, clustered into one state represented by
        // geometry 2, with state properties.
        Distribution::new(
            Members::Subset(vec![0, 2]),
            Partition::new(vec![0, 0], vec![1]).unwrap(),
        )
        .unwrap()
        .with_weights(vec![2.0, 1.0])
        .unwrap()
        .with_states(vec![{
            let mut state = State::new();
            state
                .properties_mut()
                .insert(key("cluster"), PropertyValue::String("open".into()))
                .unwrap();
            state
        }])
        .unwrap(),
    ];
    for distribution in distributions {
        let mut set = rich_set();
        set.set_distribution(Some(distribution)).unwrap();
        set.set_order(Some(Order::Step(vec![1, 2, 3]))).unwrap();
        let decoded = round_trip(&set).into_structure_set().unwrap();
        assert_same_structure_set(&decoded, &set);
    }
}

#[test]
fn interpreted_structure_sets_keep_one_report_per_sample() {
    let text = std::fs::read_to_string(fixture("mmcif/1AKE.cif")).unwrap();
    let document = mmcif::parse_str(&text, Default::default()).unwrap();
    let selections = ["A", "B"].map(|label| {
        StructureInterpretOptions::default()
            .with_altloc_policy(AltLocPolicy::SelectLabel(label.into()))
    });
    let conformations = document.interpret_conformations(&selections).unwrap();
    let bytes = native::write(&conformations, NativeWriteOptions::default()).unwrap();
    assert_eq!(bytes[12], 10);
    let interpretation = native::read(bytes.as_slice(), NativeReadOptions::default())
        .unwrap()
        .into_structure_set_interpretation()
        .unwrap();
    assert_eq!(interpretation.reports(), conformations.reports());
    assert_same_structure_set(
        interpretation.structure_set(),
        conformations.structure_set(),
    );
    let decoded = round_trip(&interpretation)
        .into_structure_set_interpretation()
        .unwrap();
    assert_eq!(decoded.reports(), interpretation.reports());
    assert_same_structure_set(decoded.structure_set(), interpretation.structure_set());

    // A set and its reports held apart write the same bytes; a set needs one
    // report per geometry.
    let (set, reports) = interpretation.clone().into_parts();
    let apart = native::write(
        NativeSource::interpreted_structure_set(&set, &reports),
        NativeWriteOptions::default(),
    )
    .unwrap();
    assert!(apart == bytes);
    let error = native::write(
        NativeSource::interpreted_structure_set(&set, &reports[..1]),
        NativeWriteOptions::default(),
    )
    .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.to_string(),
        "Cannot write the native output. The structure set has 2 geometries and 1 interpretation report. Give one report per geometry."
    );
    // A reader that wants the set alone drops the reports.
    let alone = native::read(bytes.as_slice(), NativeReadOptions::default())
        .unwrap()
        .into_structure_set()
        .unwrap();
    assert_same_structure_set(&alone, &set);
    let error = native::read(
        native::write(&set, NativeWriteOptions::default())
            .unwrap()
            .as_slice(),
        NativeReadOptions::default(),
    )
    .unwrap()
    .into_structure_set_interpretation()
    .unwrap_err();
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::NativeKind {
            expected: NativeKind::SampleSetInterpretation,
            found: NativeKind::StructureSet,
        })
    );
}

#[test]
fn structure_set_descriptions_store_the_number_of_samples() {
    let mut set = rich_set();
    set.set_order(Some(Order::Step(vec![1, 2, 3]))).unwrap();
    set.set_distribution(Some(
        Distribution::new(Members::Subset(vec![1]), Partition::one_state(1)).unwrap(),
    ))
    .unwrap();
    let description = StructureSetDescription::of(&set);
    assert_eq!(description.geometry_count(), 3);
    let decoded = round_trip(&description);
    assert_eq!(decoded.kind(), NativeKind::StructureSet);
    let error = decoded.clone().into_structure_set().unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::NativeExternalGeometries { geometries: 3 })
    );
    assert_eq!(
        error.to_string(),
        "Cannot interpret the native input. The input's structure set stores its 3 geometries elsewhere."
    );
    let decoded = decoded.into_structure_set_description().unwrap();
    assert_same_topology(decoded.topology(), set.topology());
    assert_eq!(decoded.properties(), set.properties());
    assert_eq!(decoded.geometry_count(), 3);
    assert_eq!(decoded.order(), set.order());
    assert_eq!(decoded.distribution(), set.distribution());

    // A full record describes itself too.
    let full = native::read(
        native::write(&set, NativeWriteOptions::default())
            .unwrap()
            .as_slice(),
        NativeReadOptions::default(),
    )
    .unwrap()
    .into_structure_set_description()
    .unwrap();
    assert_eq!(full.geometry_count(), 3);
    assert_eq!(full.order(), set.order());

    // The count is bounded like stored geometries, and the order and
    // distribution are checked against it.
    let bytes = native::write(&description, NativeWriteOptions::default()).unwrap();
    let limited = NativeReadOptions::default().with_max_items(2).unwrap();
    let error = native::read(bytes.as_slice(), limited).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    assert!(StructureSetDescription::new(water_topology(), 2)
        .with_order(Some(Order::Step(vec![1, 2, 3])))
        .is_err());
    for len in 0..bytes.len() {
        assert!(native::read(&bytes[..len], NativeReadOptions::default()).is_err());
    }
}

#[test]
fn distributed_and_ordered_structure_sets_are_one_kind() {
    let clustered = water_set(Partition::new(vec![1, 0, 1, 0], vec![3, 2]).unwrap());
    let mut timed = StructureSet::from_geometries(
        water_topology(),
        (0..2u32).map(|index| Geometry::new(positions(f64::from(index)))),
    )
    .unwrap();
    timed
        .set_order(Some(Order::Time(
            TimeAxis::new(Some(Quantity::new(vec![0.0, 1.0], PICOSECOND)), None).unwrap(),
        )))
        .unwrap();
    let mut scan =
        StructureSet::from_geometries(water_topology(), [Geometry::new(positions(0.0))]).unwrap();
    scan.set_order(Some(Order::Coordinate {
        unit: DEGREE,
        values: vec![-0.0],
    }))
    .unwrap();
    for (set, expected) in [
        (&clustered, 0xe094_8e71_1825_a8c0),
        (&timed, 0xf70e_517b_80a7_996a),
        (&scan, 0xaa6d_2989_3a6c_8868),
    ] {
        let bytes = native::write(set, NativeWriteOptions::default()).unwrap();
        assert_eq!(bytes[12], 9);
        assert_eq!(crate::fingerprint(&bytes), expected);
        assert_same_structure_set(&round_trip(set).into_structure_set().unwrap(), set);
    }
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

/// Four geometries of water: the first with a cell, an occupancy, and an
/// owner property, the last with the smallest positive coordinate.
fn water_samples() -> Vec<Geometry> {
    let mut first = Geometry::new(positions(0.5)).with_cell(cell());
    first.set_occupancies(Some(vec![Some(0.5), None])).unwrap();
    first
        .properties_mut()
        .owner_mut()
        .insert(key("replica"), PropertyValue::Int(4))
        .unwrap();
    vec![
        first,
        Geometry::new(positions(1.0)),
        Geometry::new(positions(2.0)),
        Geometry::new(positions(f64::MIN_POSITIVE)),
    ]
}

/// The water geometries in the states of `partition`, weighted down to the
/// smallest positive weight.
fn water_set(partition: Partition) -> StructureSet {
    let mut set = StructureSet::from_geometries(water_topology(), water_samples()).unwrap();
    set.set_distribution(Some(
        Distribution::new(Members::All, partition)
            .unwrap()
            .with_weights(vec![0.25, 3.0, 0.1, 1.0e-300])
            .unwrap(),
    ))
    .unwrap();
    set
}

/// Asserts that a set round trips and encodes to the bytes `expected`
/// fingerprints, returning the decoded set.
fn round_trip_set(set: &StructureSet, expected: u64) -> StructureSet {
    assert_eq!(crate::native_fingerprint(set), expected);
    let decoded = round_trip(set).into_structure_set().unwrap();
    assert_same_structure_set(&decoded, set);
    for (left, right) in decoded.iter().zip(set.iter()) {
        assert_eq!(
            left.weight().map(f64::to_bits),
            right.weight().map(f64::to_bits)
        );
    }
    decoded
}

#[test]
fn distributions_of_one_state_keep_samples_weights_and_site_values() {
    let mut set = water_set(Partition::one_state(4));
    set.properties_mut()
        .insert(key("method"), PropertyValue::String("REMD".into()))
        .unwrap();
    let decoded = round_trip_set(&set, 0x2e0b_dd44_1d2f_84de);
    assert_eq!(decoded.distribution().unwrap().state_count(), 1);
}

#[test]
fn distributions_of_one_state_per_sample_keep_energies_and_the_temperature() {
    let mut set = water_set(Partition::state_per_sample(4));
    let mut distribution = set.distribution_mut().unwrap();
    for (state, energy) in [(0, -12.5), (2, 0.0)] {
        distribution
            .state_mut(state)
            .unwrap()
            .set_energy(Some(Quantity::new(energy, KILOJOULE_PER_MOLE)))
            .unwrap();
    }
    distribution
        .set_temperature(Some(Quantity::new(298.15, kekule::units::KELVIN)))
        .unwrap();
    let decoded = round_trip_set(&set, 0x6372_8013_2b2e_4f7f);
    assert_eq!(decoded.distribution().unwrap().state_count(), 4);
}

#[test]
fn clustered_distributions_keep_their_partition_and_state_properties() {
    let mut set = water_set(Partition::new(vec![1, 0, 1, 0], vec![3, 2]).unwrap());
    let mut state = State::new()
        .with_energy(Quantity::new(-0.0, KILOJOULE_PER_MOLE))
        .unwrap();
    state
        .properties_mut()
        .insert(key("cluster"), PropertyValue::String("open".into()))
        .unwrap();
    *set.distribution_mut().unwrap().state_mut(1).unwrap() = state;
    let decoded = round_trip_set(&set, 0xd7e9_595d_bac6_e09f);
    assert_eq!(
        decoded
            .state(1)
            .unwrap()
            .energy()
            .unwrap()
            .value()
            .to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(decoded.state(0).unwrap().energy(), None);
}

#[test]
fn time_orders_keep_time_step_velocities_and_forces() {
    let geometries = (0..3u32).map(|index| {
        let mut geometry = Geometry::new(positions(f64::from(index)));
        if index != 1 {
            geometry.set_cell(Some(cell()));
            geometry
                .set_velocities(Some(
                    Velocities::from_vec(Quantity::new(
                        vec![Vector3::new(1.0, 2.0, 3.0); 2],
                        kekule::units::NANOMETER_PER_PICOSECOND,
                    ))
                    .unwrap(),
                ))
                .unwrap();
            geometry
                .set_forces(Some(
                    Forces::from_vec(Quantity::new(
                        vec![Vector3::new(-1.0, 0.5, 0.0); 2],
                        kekule::units::KILOJOULE_PER_MOLE_PER_NANOMETER,
                    ))
                    .unwrap(),
                ))
                .unwrap();
        }
        geometry
    });
    let mut set = StructureSet::from_geometries(water_topology(), geometries).unwrap();
    let times = (0..3u32).map(|index| f64::from(index) * 0.002).collect();
    let steps = (0..3u64).map(|index| u64::MAX - index).collect();
    set.set_order(Some(Order::Time(
        TimeAxis::new(Some(Quantity::new(times, PICOSECOND)), Some(steps)).unwrap(),
    )))
    .unwrap();
    set.get_mut(2)
        .unwrap()
        .properties_mut()
        .atoms_mut()
        .set_value(
            key("flag"),
            kekule::topology::TopologyAtomIndex::new(1),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    round_trip_set(&set, 0x4e84_32d6_85a6_eacb);

    // A time order without times or steps keeps having none.
    let mut untimed = StructureSet::from_geometries(
        water_topology(),
        [Geometry::new(positions(0.0)), Geometry::new(positions(1.0))],
    )
    .unwrap();
    untimed
        .set_order(Some(Order::Time(TimeAxis::new(None, None).unwrap())))
        .unwrap();
    let decoded = round_trip_set(&untimed, 0x1b2a_c843_8dcc_2f53);
    assert!(decoded.iter().all(|geometry| geometry.place()
        == Some(Place::Time {
            time: None,
            step: None
        })));
}

#[test]
fn step_and_coordinate_orders_keep_places_cells_and_values() {
    // An optimization history: steps, energies, a cell, and atom values.
    let geometries = (0..3).map(|index| {
        let mut geometry = Geometry::new(positions(f64::from(index)));
        if index == 1 {
            geometry = geometry.with_cell(cell());
        }
        geometry
            .properties_mut()
            .owner_mut()
            .insert(
                key("energy"),
                PropertyValue::real(-1.5 - f64::from(index), KILOJOULE_PER_MOLE).unwrap(),
            )
            .unwrap();
        geometry
    });
    let mut history = StructureSet::from_geometries(water_topology(), geometries).unwrap();
    history
        .set_order(Some(Order::Step(vec![0, 7, u64::MAX])))
        .unwrap();
    history
        .properties_mut()
        .insert(key("converged"), PropertyValue::Bool(true))
        .unwrap();
    history
        .get_mut(2)
        .unwrap()
        .properties_mut()
        .atoms_mut()
        .set_value(
            key("flag"),
            kekule::topology::TopologyAtomIndex::new(1),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    let decoded = round_trip_set(&history, 0xb76f_4cd8_a1d5_12a7);
    assert_eq!(decoded.get(2).unwrap().place(), Some(Place::Step(u64::MAX)));

    // A relaxed scan in degrees, swept back, with a negative zero.
    let mut scan = StructureSet::from_geometries(
        water_topology(),
        (0..4).map(|index| Geometry::new(positions(f64::from(index)))),
    )
    .unwrap();
    scan.set_order(Some(Order::Coordinate {
        unit: DEGREE,
        values: vec![-0.0, 1e-300, 180.0, -90.5],
    }))
    .unwrap();
    let decoded = round_trip_set(&scan, 0x02a0_cc46_8fab_4ee9);
    let Some(Place::Coordinate(first)) = decoded.get(0).unwrap().place() else {
        panic!("a coordinate order places its geometries by coordinate");
    };
    assert_eq!(first.value().to_bits(), (-0.0_f64).to_bits());
    assert_eq!(first.unit(), DEGREE);

    // An empty set keeps its order's unit.
    let mut empty = StructureSet::new(water_topology());
    empty
        .set_order(Some(Order::Coordinate {
            unit: ANGSTROM,
            values: Vec::new(),
        }))
        .unwrap();
    let decoded = round_trip_set(&empty, 0x93e3_20f6_6457_161b);
    assert!(decoded.is_empty());

    // The geometry count is bounded; every truncation fails cleanly.
    let bytes = native::write(&scan, NativeWriteOptions::default()).unwrap();
    let decoded = native::read(bytes.as_slice(), NativeReadOptions::default()).unwrap();
    assert_eq!(decoded.kind(), NativeKind::StructureSet);
    let limited = NativeReadOptions::default().with_max_items(3).unwrap();
    let error = native::read(bytes.as_slice(), limited).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    for len in 0..bytes.len() {
        assert!(native::read(&bytes[..len], NativeReadOptions::default()).is_err());
    }
}

#[test]
fn models_keep_cells_occupancies_b_factors_and_structure_properties() {
    let topology = rich_topology();
    let points = (0..7)
        .map(|index| Point3::new(f64::from(index), -f64::from(index), 0.125))
        .collect::<Vec<_>>();
    let mut model = Structure::new(
        topology,
        Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
    )
    .unwrap();
    model.set_cell(Some(cell()));
    model
        .set_occupancies(Some(vec![
            Some(1.0),
            None,
            Some(0.5),
            None,
            None,
            None,
            Some(0.0),
        ]))
        .unwrap();
    model
        .set_b_factors(Some(Quantity::new(
            vec![Some(20.0); 7],
            kekule::units::SQUARE_ANGSTROM,
        )))
        .unwrap();
    model
        .properties_mut()
        .bonds_mut()
        .set_value(
            key("order_hint"),
            kekule::topology::TopologyBondIndex::new(0),
            Some(PropertyValue::Int(1)),
        )
        .unwrap();
    let decoded = round_trip(&model).into_structure().unwrap();
    assert_same_structure(&decoded, &model);
}

#[test]
fn perception_is_optional_on_both_sides() {
    let molecule = perceived("c1ccccc1");
    let without = native::write(
        &molecule,
        NativeWriteOptions::default().with_perception(false),
    )
    .unwrap();
    let decoded = native::read(without.as_slice(), NativeReadOptions::default())
        .unwrap()
        .into_molecule()
        .unwrap();
    assert_eq!(decoded, molecule);
    assert_eq!(decoded.perception(), &Default::default());

    let with = native::write(&molecule, NativeWriteOptions::default()).unwrap();
    let discarded = native::read(
        with.as_slice(),
        NativeReadOptions::default().with_perception(false),
    )
    .unwrap()
    .into_molecule()
    .unwrap();
    assert_eq!(discarded.perception(), &Default::default());
}

#[test]
fn objects_follow_each_other_in_one_stream() {
    let molecule = perceived("CCO");
    let topology = water_topology();
    let mut stream = Vec::new();
    native::write_to(&mut stream, &molecule, NativeWriteOptions::default()).unwrap();
    native::write_to(
        &mut stream,
        topology.as_ref(),
        NativeWriteOptions::default(),
    )
    .unwrap();
    let mut reader = stream.as_slice();
    let first = native::read(&mut reader, NativeReadOptions::default()).unwrap();
    let second = native::read(&mut reader, NativeReadOptions::default()).unwrap();
    assert!(reader.is_empty());
    assert_eq!(first.kind(), NativeKind::Molecule);
    let error = second.into_structure().unwrap_err();
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::NativeKind {
            expected: NativeKind::Structure,
            found: NativeKind::Topology,
        })
    );
}

#[test]
fn readers_reject_foreign_magic_other_versions_and_unknown_flags() {
    let bytes = native::write(&perceived("C"), NativeWriteOptions::default()).unwrap();
    let error = native::read(&b"not a native file"[..], NativeReadOptions::default()).unwrap_err();
    assert_eq!(error.format(), Format::Native);
    assert_eq!(error.kind(), DiagnosticKind::Syntax);
    assert_eq!(error.detail(), Some(&FormatErrorDetail::NativeMagic));

    let mut newer = bytes.clone();
    newer[8..12].copy_from_slice(&(NATIVE_FORMAT_VERSION + 1).to_le_bytes());
    let error = native::read(newer.as_slice(), NativeReadOptions::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::NativeVersion {
            found: NATIVE_FORMAT_VERSION + 1,
            supported: NATIVE_FORMAT_VERSION,
        })
    );
    assert_eq!(error.byte_offset(), Some(12));

    // Version 4, whose ensemble, trajectory, and series kinds version 5
    // replaces with structure sets, is rejected with no migration.
    for version in [2u32, 4] {
        let mut older = bytes.clone();
        older[8..12].copy_from_slice(&version.to_le_bytes());
        let error = native::read(older.as_slice(), NativeReadOptions::default()).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Unsupported);
        assert_eq!(
            error.detail(),
            Some(&FormatErrorDetail::NativeVersion {
                found: version,
                supported: NATIVE_FORMAT_VERSION,
            })
        );
        assert!(
            error.to_string().contains(&format!("version {version}")),
            "{error}"
        );
    }

    let mut flagged = bytes;
    flagged[13] |= 0b10;
    let error = native::read(flagged.as_slice(), NativeReadOptions::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);
}

#[test]
fn every_truncation_and_byte_corruption_fails_cleanly() {
    let mut molecule = perceived("C[C@H](N)c1ccccc1");
    molecule
        .properties_mut()
        .atoms_mut()
        .set_value(
            key("label"),
            AtomId::new(1),
            Some(PropertyValue::String("x".into())),
        )
        .unwrap();
    let topology = Topology::from_molecule(molecule).unwrap();
    let bytes = native::write(&topology, NativeWriteOptions::default()).unwrap();
    for len in 0..bytes.len() {
        let error = native::read(&bytes[..len], NativeReadOptions::default()).unwrap_err();
        assert_eq!(error.format(), Format::Native);
    }
    for index in 0..bytes.len() {
        for flip in [0x01, 0x80, 0xFF] {
            let mut corrupted = bytes.clone();
            corrupted[index] ^= flip;
            assert!(
                native::read(corrupted.as_slice(), NativeReadOptions::default()).is_err(),
                "corrupting byte {index} with {flip:#04x} went unnoticed"
            );
        }
    }
}

#[test]
fn limits_bound_declared_counts_before_allocation() {
    let topology = Topology::from_molecules(smiles::to_molecules("CCO.O.O").unwrap()).unwrap();
    let bytes = native::write(&topology, NativeWriteOptions::default()).unwrap();
    let cases = [
        (
            NativeReadOptions::default().with_max_atoms(4).unwrap(),
            Resource::Atoms,
        ),
        (
            NativeReadOptions::default().with_max_bonds(1).unwrap(),
            Resource::Bonds,
        ),
        (
            NativeReadOptions::default()
                .with_max_input_bytes(64)
                .unwrap(),
            Resource::InputBytes,
        ),
    ];
    for (options, resource) in cases {
        let error = native::read(bytes.as_slice(), options).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Limit);
        assert_eq!(
            error.limit_exceeded().map(|limit| limit.resource()),
            Some(resource)
        );
    }

    // A count far beyond the input is rejected without reserving its storage.
    let mut huge = bytes;
    huge[16..20].copy_from_slice(&u32::MAX.to_le_bytes());
    let error = native::read(huge.as_slice(), NativeReadOptions::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    assert!(NativeReadOptions::default().with_max_items(0).is_err());
}

#[test]
fn stored_ids_resolve_against_the_decoded_owner() {
    let topology = rich_topology();
    let stored = topology
        .atom_ids()
        .iter()
        .map(|id| id.to_bits())
        .collect::<Vec<_>>();
    let decoded = round_trip(&topology).into_topology().unwrap();
    for (index, bits) in stored.into_iter().enumerate() {
        let id = InstanceAtomId::from_bits(bits);
        assert_eq!(
            decoded.atom_index(id).map(|index| index.index()),
            Some(index)
        );
    }
    let map = NonZeroU32::new(7);
    let molecule = round_trip(&perceived("[CH3:7]C")).into_molecule().unwrap();
    assert_eq!(molecule.atom(AtomId::new(0)).unwrap().atom_map(), map);
}

/// A solvated system of about one million atoms: one protein-sized molecule
/// and water, every atom placed in the hierarchy.
fn million_atom_model() -> Structure {
    let water = smiles::to_molecules("[H]O[H]").unwrap().remove(0);
    let chain_smiles = "NCC(=O)".repeat(2000) + "O";
    let peptide = smiles::to_molecules(&chain_smiles).unwrap().remove(0);
    let mut builder = TopologyBuilder::new();
    let peptide_instance = builder.add_molecule(peptide).unwrap();
    let water_definition = builder.add_molecule_definition(water).unwrap();
    let peptide_atoms = builder.atom_count();
    let waters = (1_000_000 - peptide_atoms) / 3;
    builder.reserve_instances(waters).unwrap();
    let instances = (0..waters)
        .map(|_| builder.add_instance(water_definition).unwrap())
        .collect::<Vec<_>>();
    let hierarchy = builder.hierarchy_mut();
    let protein = hierarchy.add_chain("A", Some("A".into())).unwrap();
    let residue = hierarchy
        .add_residue(protein, "GLY", Some(1), Some("1".into()), None)
        .unwrap();
    for atom in 0..peptide_atoms {
        let id = InstanceAtomId::new(peptide_instance, AtomId::new(atom as u32));
        hierarchy
            .add_atom_site(residue, id, AtomSiteMetadata::default())
            .unwrap();
    }
    let solvent = hierarchy.add_chain("W", Some("W".into())).unwrap();
    for (number, instance) in instances.into_iter().enumerate() {
        let residue = hierarchy
            .add_residue(solvent, "HOH", None, Some(number.to_string()), None)
            .unwrap();
        for (atom, name) in ["H1", "O", "H2"].into_iter().enumerate() {
            hierarchy
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, AtomId::new(atom as u32)),
                    AtomSiteMetadata::default().with_label_atom_id(Some(name.into())),
                )
                .unwrap();
        }
    }
    let topology = builder.build().unwrap();
    let points = (0..topology.atom_count())
        .map(|index| {
            let index = index as f64;
            Point3::new(index * 0.001, (index * 0.37).sin(), (index * 0.11).cos())
        })
        .collect::<Vec<_>>();
    Structure::new(
        topology,
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

#[test]
#[ignore = "timing of a one-million-atom round trip; run in release with --ignored"]
fn million_atom_model_round_trips_within_a_second() {
    let model = million_atom_model();
    assert!(model.atom_count() >= 999_000);
    let start = std::time::Instant::now();
    let bytes = native::write(&model, NativeWriteOptions::default()).unwrap();
    let encoded = start.elapsed();
    let decoded = native::read(bytes.as_slice(), NativeReadOptions::default())
        .unwrap()
        .into_structure()
        .unwrap();
    let total = start.elapsed();
    eprintln!(
        "{} atoms, {} bytes: encode {:?}, decode {:?}, round trip {:?}",
        model.atom_count(),
        bytes.len(),
        encoded,
        total - encoded,
        total
    );
    assert_same_structure(&decoded, &model);
    if !cfg!(debug_assertions) {
        assert!(
            total < std::time::Duration::from_secs(1),
            "round trip took {total:?}"
        );
    }
}

#[test]
fn empty_hierarchy_nodes_count_against_items_not_atoms() {
    let mut builder = TopologyBuilder::new();
    builder
        .add_molecule(smiles::to_molecules("[Na+]").unwrap().remove(0))
        .unwrap();
    let hierarchy = builder.hierarchy_mut();
    for label in ["A", "B", "C"] {
        let chain = hierarchy.add_chain(label, None).unwrap();
        hierarchy
            .add_residue(chain, "UNK", None, None, None)
            .unwrap();
    }
    let topology = builder.build().unwrap();
    let bytes = native::write(&topology, NativeWriteOptions::default()).unwrap();
    // One atom, three empty chains and residues: the atom limit does not
    // bound hierarchy nodes.
    let options = NativeReadOptions::default().with_max_atoms(1).unwrap();
    let decoded = native::read(bytes.as_slice(), options)
        .unwrap()
        .into_topology()
        .unwrap();
    assert_same_topology(&decoded, &topology);
    let options = NativeReadOptions::default().with_max_items(2).unwrap();
    let error = native::read(bytes.as_slice(), options).unwrap_err();
    assert_eq!(
        error.limit_exceeded().map(|limit| limit.resource()),
        Some(Resource::Residues)
    );
}
