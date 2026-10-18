//! Native serialization: exact round trips, canonical bytes, and rejection of
//! malformed input without panics.

use std::num::NonZeroU32;
use std::sync::Arc;

use kekule::diagnostics::{DiagnosticKind, Format, FormatErrorDetail, Resource};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{AtomId, Molecule};
use kekule::native::{
    self, NativeKind, NativeObject, NativeReadOptions, NativeSource, NativeWriteOptions,
    NATIVE_FORMAT_VERSION,
};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{
    Ensemble, EnsembleMember, Forces, Model, Positions, Trajectory, TrajectoryFrame, Velocities,
};
use kekule::topology::{
    AtomSiteMetadata, InstanceAtomId, MoleculeClass, ResidueClass, Topology, TopologyBuilder,
};
use kekule::units::{
    Quantity, ANGSTROM, DALTON, ELEMENTARY_CHARGE, KILOJOULE_PER_MOLE, NANOMETER, PICOSECOND,
};
use kekule::{mmcif, smiles};

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
        NativeObject::Model(value) => native::write(value, NativeWriteOptions::default()),
        NativeObject::Ensemble(value) => native::write(value, NativeWriteOptions::default()),
        NativeObject::Trajectory(value) => native::write(value, NativeWriteOptions::default()),
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

fn assert_same_model(decoded: &Model, source: &Model) {
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
            .into_model();
        model.perceive().unwrap();
        model
            .properties_mut()
            .owner_mut()
            .insert(key("source"), PropertyValue::String(name.into()))
            .unwrap();
        let decoded = round_trip(&model).into_model().unwrap();
        assert_same_model(&decoded, &model);
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

#[test]
fn ensembles_keep_weights_energies_and_site_values() {
    let topology = water_topology();
    let mut first = EnsembleMember::new(positions(0.5), 0.25)
        .unwrap()
        .with_cell(cell())
        .with_energy(Quantity::new(-12.5, KILOJOULE_PER_MOLE))
        .unwrap();
    first.set_occupancies(Some(vec![Some(0.5), None])).unwrap();
    let second = EnsembleMember::new(positions(f64::MIN_POSITIVE), 3.0).unwrap();
    let mut ensemble = Ensemble::from_items(topology, [first, second]).unwrap();
    ensemble
        .properties_mut()
        .insert(key("method"), PropertyValue::String("REMD".into()))
        .unwrap();
    ensemble
        .get_mut(1)
        .unwrap()
        .properties_mut()
        .owner_mut()
        .insert(key("replica"), PropertyValue::Int(4))
        .unwrap();
    let decoded = round_trip(&ensemble).into_ensemble().unwrap();
    assert_same_topology(decoded.topology(), ensemble.topology());
    assert_eq!(decoded.properties(), ensemble.properties());
    assert_eq!(decoded.len(), 2);
    for (left, right) in decoded.iter().zip(ensemble.iter()) {
        assert_eq!(left.payload(), right.payload());
        assert_eq!(
            left.payload().weight().to_bits(),
            right.payload().weight().to_bits()
        );
    }
}

#[test]
fn trajectories_keep_time_step_velocities_and_forces() {
    let topology = water_topology();
    let mut frames = Vec::new();
    for index in 0..3u32 {
        let mut frame = TrajectoryFrame::new(positions(f64::from(index)));
        if index != 1 {
            frame.set_cell(Some(cell()));
            frame
                .set_time(Some(Quantity::new(f64::from(index) * 0.002, PICOSECOND)))
                .unwrap();
            frame.set_step(Some(u64::MAX - u64::from(index)));
            frame
                .set_velocities(Some(
                    Velocities::from_vec(Quantity::new(
                        vec![Vector3::new(1.0, 2.0, 3.0); 2],
                        kekule::units::NANOMETER_PER_PICOSECOND,
                    ))
                    .unwrap(),
                ))
                .unwrap();
            frame
                .set_forces(Some(
                    Forces::from_vec(Quantity::new(
                        vec![Vector3::new(-1.0, 0.5, 0.0); 2],
                        kekule::units::KILOJOULE_PER_MOLE_PER_NANOMETER,
                    ))
                    .unwrap(),
                ))
                .unwrap();
        }
        frames.push(frame);
    }
    let mut trajectory = Trajectory::from_items(topology, frames).unwrap();
    trajectory
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
    let decoded = round_trip(&trajectory).into_trajectory().unwrap();
    assert_same_topology(decoded.topology(), trajectory.topology());
    for (left, right) in decoded.iter().zip(trajectory.iter()) {
        assert_eq!(left.payload(), right.payload());
    }
}

#[test]
fn models_keep_cells_occupancies_b_factors_and_structure_properties() {
    let topology = rich_topology();
    let points = (0..7)
        .map(|index| Point3::new(f64::from(index), -f64::from(index), 0.125))
        .collect::<Vec<_>>();
    let mut model = Model::new(
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
    let decoded = round_trip(&model).into_model().unwrap();
    assert_same_model(&decoded, &model);
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
    let error = second.into_model().unwrap_err();
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::NativeKind {
            expected: NativeKind::Model,
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
fn million_atom_model() -> Model {
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
    Model::new(
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
        .into_model()
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
    assert_same_model(&decoded, &model);
    if !cfg!(debug_assertions) {
        assert!(
            total < std::time::Duration::from_secs(1),
            "round trip took {total:?}"
        );
    }
}
