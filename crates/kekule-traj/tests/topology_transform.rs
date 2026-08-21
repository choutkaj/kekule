use kekule::bio::SmcraAtomSiteMetadata;
use kekule::core::{Atom, AtomId, BondId, BondOrder, Element, Molecule, PropValue};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::modeling::potential::{HarmonicBondParameter, HarmonicBondPotential, Potential};
use kekule::structure::{AtomData, Ensemble, EnsembleMember, Model, Positions, TopologyRemapError};
use kekule::topology::transform::{
    remove_instances, retain_instances, RemovedSelectionPolicy, SelectionRemapError,
    TopologyTransformError,
};
use kekule::topology::{
    AtomSelection, InstanceAtomId, InstanceBondId, MoleculeDefinitionId, MoleculeInstanceId,
    MoleculeInstanceMetadata, MoleculeRole, Topology, TopologyBuilder, TopologyMapping,
};
use kekule::units::{
    Quantity, ANGSTROM, MODEL_FORCE_CONSTANT_UNIT, MODEL_FORCE_UNIT, MODEL_VELOCITY_UNIT,
    PICOSECOND, SQUARE_ANGSTROM,
};
use kekule_traj::{
    Forces, FrameBuffer, Trajectory, TrajectoryFrame, TrajectoryRemapError, Velocities,
};

struct Fixture {
    topology: Arc<Topology>,
    water_definition: MoleculeDefinitionId,
    ligand_definition: MoleculeDefinitionId,
    macro_definition: MoleculeDefinitionId,
    ion_definition: MoleculeDefinitionId,
    water_first: MoleculeInstanceId,
    ligand: MoleculeInstanceId,
    water_second: MoleculeInstanceId,
    macromolecule: MoleculeInstanceId,
    ion: MoleculeInstanceId,
}

fn element(symbol: &str) -> Element {
    Element::from_symbol(symbol).expect("test element")
}

fn one_atom_small(symbol: &str) -> Molecule {
    let mut graph = kekule::core::MoleculeEditor::new();
    graph
        .add_atom(Atom::new(element(symbol)))
        .expect("test atom capacity");
    graph.finish().expect("connected single atom")
}

fn ligand_with_tombstones() -> Molecule {
    let mut graph = kekule::core::MoleculeEditor::new();
    let carbon = graph
        .add_atom(Atom::new(element("C")))
        .expect("test atom capacity");
    let deleted = graph
        .add_atom(Atom::new(element("H")))
        .expect("test atom capacity");
    let oxygen = graph
        .add_atom(Atom::new(element("O")))
        .expect("test atom capacity");
    graph
        .add_bond(carbon, deleted, BondOrder::Single)
        .expect("test bond capacity");
    graph
        .add_bond(carbon, oxygen, BondOrder::Double)
        .expect("test bond capacity");
    graph.delete_atom(deleted).expect("delete test tombstone");
    graph.finish().expect("connected ligand")
}

fn one_atom_macro() -> Molecule {
    let mut builder = kekule::core::MoleculeEditor::new();
    let atom = builder
        .add_atom(Atom::new(element("C")))
        .expect("test atom capacity");
    let chain = builder
        .hierarchy_mut()
        .add_chain("A", Some("auth-A".to_owned()))
        .expect("test chain");
    let residue = builder
        .hierarchy_mut()
        .add_residue(
            chain,
            "GLY",
            Some(7),
            Some("17".to_owned()),
            Some("B".to_owned()),
        )
        .expect("test residue");
    builder
        .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
        .expect("test atom site");
    builder.finish().expect("valid test macromolecule")
}

fn metadata(role: MoleculeRole, label: &str) -> MoleculeInstanceMetadata {
    let mut metadata = MoleculeInstanceMetadata::default();
    metadata.insert_role(role);
    metadata
        .props_mut()
        .insert("label".to_owned(), PropValue::String(label.to_owned()));
    metadata
}

fn fixture() -> Fixture {
    let water = one_atom_small("O");
    let ligand = ligand_with_tombstones();
    let macromolecule = one_atom_macro();
    let ion = one_atom_small("Na");

    let mut builder = TopologyBuilder::new();
    let water_definition = builder
        .add_molecule_definition(&water)
        .expect("water definition");
    let ligand_definition = builder
        .add_molecule_definition(&ligand)
        .expect("ligand definition");
    let macro_definition = builder
        .add_molecule_definition(&macromolecule)
        .expect("macro definition");
    let ion_definition = builder
        .add_molecule_definition(&ion)
        .expect("ion definition");

    let water_first = builder
        .add_instance(water_definition, metadata(MoleculeRole::Solvent, "water-0"))
        .expect("water instance");
    let ligand = builder
        .add_instance(ligand_definition, metadata(MoleculeRole::Ligand, "ligand"))
        .expect("ligand instance");
    let water_second = builder
        .add_instance(water_definition, metadata(MoleculeRole::Solvent, "water-1"))
        .expect("water instance");
    let macromolecule = builder
        .add_instance(macro_definition, metadata(MoleculeRole::Polymer, "protein"))
        .expect("macro instance");
    let ion = builder
        .add_instance(ion_definition, metadata(MoleculeRole::Ion, "sodium"))
        .expect("ion instance");

    Fixture {
        topology: Arc::new(builder.build().expect("test topology")),
        water_definition,
        ligand_definition,
        macro_definition,
        ion_definition,
        water_first,
        ligand,
        water_second,
        macromolecule,
        ion,
    }
}

fn retained_fixture(fixture: &Fixture) -> kekule::topology::TopologyEditResult {
    retain_instances(
        &fixture.topology,
        [
            fixture.macromolecule,
            fixture.water_second,
            fixture.ligand,
            fixture.water_second,
        ],
    )
    .expect("valid retained subset")
}

#[test]
fn whole_instance_subset_is_deterministic_complete_and_non_mutating() {
    let fixture = fixture();
    let source_clone = Arc::clone(&fixture.topology);
    let edit = retained_fixture(&fixture);
    let target = edit.shared_topology();
    let mapping = edit.mapping();

    assert!(Arc::ptr_eq(&fixture.topology, &source_clone));
    assert_eq!(fixture.topology.definition_count(), 4);
    assert_eq!(fixture.topology.instance_count(), 5);
    assert!(!Arc::ptr_eq(&fixture.topology, &target));
    assert_eq!(target.definition_count(), 3);
    assert_eq!(target.instance_count(), 3);

    assert_eq!(
        mapping.definition_pairs().collect::<Vec<_>>(),
        vec![
            (fixture.water_definition, MoleculeDefinitionId::new(0)),
            (fixture.ligand_definition, MoleculeDefinitionId::new(1)),
            (fixture.macro_definition, MoleculeDefinitionId::new(2)),
        ]
    );
    assert_eq!(
        mapping.instance_pairs().collect::<Vec<_>>(),
        vec![
            (fixture.ligand, MoleculeInstanceId::new(0)),
            (fixture.water_second, MoleculeInstanceId::new(1)),
            (fixture.macromolecule, MoleculeInstanceId::new(2)),
        ]
    );
    assert!(mapping.is_source(&fixture.topology));
    assert!(mapping.is_target(&target));
    assert!(Arc::ptr_eq(&mapping.source_arc(), &fixture.topology));
    assert!(Arc::ptr_eq(&mapping.target_arc(), &target));
    assert!(Arc::ptr_eq(&edit.shared_topology(), &target));
    assert_eq!(mapping.removed_definitions(), &[fixture.ion_definition]);
    assert_eq!(
        mapping.removed_instances(),
        &[fixture.water_first, fixture.ion]
    );
    assert!(mapping.added_definitions().is_empty());
    assert!(mapping.added_instances().is_empty());
    assert!(mapping.added_atoms().is_empty());
    assert!(mapping.added_bonds().is_empty());

    assert_eq!(
        target
            .graph_for_instance(MoleculeInstanceId::new(0))
            .expect("retained ligand")
            .atom_ids()
            .collect::<Vec<_>>(),
        vec![AtomId::new(0), AtomId::new(2)]
    );
    assert_eq!(
        target
            .graph_for_instance(MoleculeInstanceId::new(0))
            .expect("retained ligand")
            .bond_ids()
            .collect::<Vec<_>>(),
        vec![BondId::new(1)]
    );
    assert_eq!(
        target
            .instance(MoleculeInstanceId::new(0))
            .expect("retained ligand")
            .metadata(),
        fixture
            .topology
            .instance(fixture.ligand)
            .expect("source ligand")
            .metadata()
    );
    assert_eq!(
        target
            .definition(MoleculeDefinitionId::new(2))
            .expect("retained macro")
            .hierarchy(),
        fixture
            .topology
            .definition(fixture.macro_definition)
            .expect("source macro")
            .hierarchy()
    );

    assert_eq!(mapping.atom_pairs().len(), target.atom_count());
    assert_eq!(mapping.bond_pairs().len(), target.bond_count());
    assert_eq!(mapping.atom_index_pairs().len(), target.atom_count());
    assert_eq!(mapping.bond_index_pairs().len(), target.bond_count());
    for (source_atom, target_atom) in mapping.atom_pairs() {
        assert_eq!(source_atom.atom(), target_atom.atom());
        assert_eq!(
            mapping.map_atom_index(
                fixture
                    .topology
                    .atom_index(source_atom)
                    .expect("mapped source atom")
            ),
            target.atom_index(target_atom)
        );
    }
    for (source_bond, target_bond) in mapping.bond_pairs() {
        assert_eq!(source_bond.bond(), target_bond.bond());
        assert_eq!(
            mapping.map_bond_index(
                fixture
                    .topology
                    .bond_index(source_bond)
                    .expect("mapped source bond")
            ),
            target.bond_index(target_bond)
        );
    }
}

#[test]
fn subset_normalization_no_ops_and_failures_are_explicit() {
    let fixture = fixture();
    let mut all = fixture
        .topology
        .instances()
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    all.reverse();
    all.push(fixture.ligand);
    let retained = retain_instances(&fixture.topology, all).expect("retain-all no-op");
    assert!(Arc::ptr_eq(&retained.shared_topology(), &fixture.topology));
    assert_eq!(
        retained.mapping().atom_pairs().len(),
        fixture.topology.atom_count()
    );

    let removed = remove_instances(&fixture.topology, []).expect("remove-none no-op");
    assert!(Arc::ptr_eq(&removed.shared_topology(), &fixture.topology));

    assert!(matches!(
        retain_instances(&fixture.topology, []),
        Err(TopologyTransformError::EmptyTargetTopology)
    ));
    assert!(matches!(
        remove_instances(
            &fixture.topology,
            fixture.topology.instances().map(|(id, _)| id),
        ),
        Err(TopologyTransformError::EmptyTargetTopology)
    ));
    let invalid = MoleculeInstanceId::new(99);
    assert!(matches!(
        retain_instances(&fixture.topology, [fixture.ligand, invalid]),
        Err(TopologyTransformError::InvalidSourceInstance(id)) if id == invalid
    ));
    assert!(matches!(
        remove_instances(&fixture.topology, [invalid, invalid]),
        Err(TopologyTransformError::InvalidSourceInstance(id)) if id == invalid
    ));

    let removed = remove_instances(
        &fixture.topology,
        [fixture.water_first, fixture.ion, fixture.water_first],
    )
    .expect("valid remove subset");
    assert_eq!(removed.topology().instance_count(), 3);
    assert_eq!(
        removed.mapping().removed_instances(),
        &[fixture.water_first, fixture.ion]
    );

    let complex = retain_instances(&fixture.topology, [fixture.macromolecule, fixture.ligand])
        .expect("protein-ligand complex");
    assert_eq!(complex.topology().definition_count(), 2);
    assert_eq!(complex.topology().instance_count(), 2);
    assert!(complex
        .topology()
        .instances()
        .all(|(_, instance)| !instance.has_role(MoleculeRole::Solvent)
            && !instance.has_role(MoleculeRole::Ion)));
    assert_eq!(
        complex.mapping().removed_definitions(),
        &[fixture.water_definition, fixture.ion_definition]
    );
}

fn point_values(topology: &Topology, offset: f64) -> Vec<Point3> {
    (0..topology.atom_count())
        .map(|index| Point3::new(offset + index as f64, index as f64 + 0.25, -0.5))
        .collect()
}

fn atom_data(topology: &Arc<Topology>) -> AtomData {
    let mut atom_data = AtomData::new(topology);
    atom_data
        .set_occupancies(
            (0..topology.atom_count())
                .map(|index| Some(0.5 + index as f64 / 100.0))
                .collect::<Vec<_>>(),
        )
        .expect("complete occupancy column");
    atom_data
        .set_b_factors(Quantity::new(
            (0..topology.atom_count())
                .map(|index| Some(10.0 + index as f64))
                .collect::<Vec<_>>(),
            SQUARE_ANGSTROM,
        ))
        .expect("complete B-factor column");
    atom_data
}

fn model(fixture: &Fixture, offset: f64) -> Model {
    let positions = Positions::new(
        &fixture.topology,
        Quantity::new(point_values(&fixture.topology, offset), ANGSTROM),
    )
    .expect("complete test positions");
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(20.0 + offset, 21.0, 22.0), ANGSTROM),
        [true, true, false],
    )
    .expect("valid test cell");
    Model::with_atom_data(
        Arc::clone(&fixture.topology),
        positions,
        Some(cell),
        atom_data(&fixture.topology),
    )
    .expect("valid test model")
}

#[test]
fn model_atom_data_and_selection_remapping_preserve_complete_state() {
    let fixture = fixture();
    let edit = retained_fixture(&fixture);
    let target_topology = edit.shared_topology();
    let source = model(&fixture, 1.0);
    let source_clone = source.clone();
    let target = source
        .remap_to(&target_topology, edit.mapping())
        .expect("model remap");

    assert_eq!(source, source_clone);
    assert!(Arc::ptr_eq(&target.shared_topology(), &target_topology));
    assert_eq!(target.cell(), source.cell());
    let source_atom_data = source.atom_data();
    let target_atom_data = target.atom_data();
    for (source_index, target_index) in edit.mapping().atom_index_pairs() {
        assert_eq!(
            target.position_at(target_index).expect("target position"),
            source.position_at(source_index).expect("source position")
        );
        assert_eq!(
            target_atom_data.occupancy_at(target_index),
            source_atom_data.occupancy_at(source_index)
        );
        assert_eq!(
            target_atom_data.b_factor_at(target_index),
            source_atom_data.b_factor_at(source_index)
        );
    }

    let removed_atom = InstanceAtomId::new(fixture.water_first, AtomId::new(0));
    let retained_atom = InstanceAtomId::new(fixture.ligand, AtomId::new(2));
    let selection = AtomSelection::from_atoms(&fixture.topology, [removed_atom, retained_atom])
        .expect("source selection");
    assert_eq!(
        selection.remap_to(
            &fixture.topology,
            &target_topology,
            edit.mapping(),
            RemovedSelectionPolicy::Error,
        ),
        Err(SelectionRemapError::RemovedSelectedAtom(removed_atom))
    );
    let dropped = selection
        .remap_to(
            &fixture.topology,
            &target_topology,
            edit.mapping(),
            RemovedSelectionPolicy::Drop,
        )
        .expect("explicit selection drop");
    assert_eq!(
        dropped
            .semantic_ids(&target_topology)
            .expect("target selection"),
        vec![edit
            .mapping()
            .map_atom(retained_atom)
            .expect("retained atom")]
    );
    let empty = AtomSelection::from_atoms(&fixture.topology, [removed_atom])
        .expect("removed-only selection")
        .remap_to(
            &fixture.topology,
            &target_topology,
            edit.mapping(),
            RemovedSelectionPolicy::Drop,
        )
        .expect("explicit empty selection");
    assert!(empty.indices().is_empty());
}

fn mapping_with_added_atom() -> (
    Arc<Topology>,
    Arc<Topology>,
    TopologyMapping,
    InstanceAtomId,
) {
    let one = one_atom_small("C");
    let mut source_builder = TopologyBuilder::new();
    let source_definition = source_builder
        .add_molecule_definition(&one)
        .expect("source definition");
    let source_instance = source_builder
        .add_instance(source_definition, MoleculeInstanceMetadata::default())
        .expect("source instance");
    let source_topology = Arc::new(source_builder.build().expect("source topology"));
    let mut target_builder = TopologyBuilder::new();
    let target_definition = target_builder
        .add_molecule_definition(&one)
        .expect("target definition");
    let target_first = target_builder
        .add_instance(target_definition, MoleculeInstanceMetadata::default())
        .expect("target instance");
    let target_added = target_builder
        .add_instance(target_definition, MoleculeInstanceMetadata::default())
        .expect("target added instance");
    let target_topology = Arc::new(target_builder.build().expect("target topology"));
    let mapping = TopologyMapping::from_pairs(
        &source_topology,
        &target_topology,
        [(source_definition, target_definition)],
        [(source_instance, target_first)],
        [(
            InstanceAtomId::new(source_instance, AtomId::new(0)),
            InstanceAtomId::new(target_first, AtomId::new(0)),
        )],
        [],
    )
    .expect("valid deletion/addition lineage");
    (
        source_topology,
        target_topology,
        mapping,
        InstanceAtomId::new(target_added, AtomId::new(0)),
    )
}

#[test]
fn state_remapping_rejects_equal_layout_substitutes_and_unmapped_target_atoms() {
    let fixture = fixture();
    let edit = retained_fixture(&fixture);
    let target_topology = edit.shared_topology();
    let source = model(&fixture, 2.0);
    let independent_source = self::fixture();
    let independent_edit = retained_fixture(&independent_source);
    let independent_target = independent_edit.shared_topology();

    assert!(fixture.topology.same_layout(&independent_source.topology));
    assert!(!Arc::ptr_eq(
        &fixture.topology,
        &independent_source.topology
    ));
    let source_positions = source.positions().clone();
    assert_eq!(
        source_positions.remap_to(
            &independent_source.topology,
            &target_topology,
            edit.mapping(),
        ),
        Err(TopologyRemapError::SourceTopologyMismatch)
    );
    assert!(target_topology.same_layout(&independent_target));
    assert_eq!(
        source_positions.remap_to(&fixture.topology, &independent_target, edit.mapping(),),
        Err(TopologyRemapError::MappingTargetMismatch)
    );

    let mut wrong_destination = Positions::zeros(&independent_target);
    assert_eq!(
        wrong_destination.copy_remapped_from(
            &source_positions,
            &fixture.topology,
            &target_topology,
            edit.mapping(),
        ),
        Err(TopologyRemapError::TargetTopologyMismatch)
    );
    let mut destination = Positions::zeros(&target_topology);
    let allocation = destination.values().value().as_ptr();
    destination
        .copy_remapped_from(
            &source_positions,
            &fixture.topology,
            &target_topology,
            edit.mapping(),
        )
        .expect("valid allocation-reusing position remap");
    assert_eq!(destination.values().value().as_ptr(), allocation);

    let (source_topology, target_topology, mapping, target_added) = mapping_with_added_atom();
    let positions = Positions::new(
        &source_topology,
        Quantity::new(vec![Point3::new(1.0, 2.0, 3.0)], ANGSTROM),
    )
    .expect("source state");
    assert_eq!(
        positions.remap_to(&source_topology, &target_topology, &mapping),
        Err(TopologyRemapError::AddedAtomsRequireState {
            target_atom: target_added,
        })
    );
}

#[test]
fn ensemble_remapping_preserves_member_state_and_reports_member_context() {
    let fixture = fixture();
    let edit = retained_fixture(&fixture);
    let target_topology = edit.shared_topology();
    let first = model(&fixture, 3.0);
    let second = model(&fixture, 30.0);
    let mut first_member = EnsembleMember::new(first.positions().clone());
    first_member.set_cell(first.cell().copied());
    first_member.set_weight(Some(0.25)).expect("valid weight");
    first_member
        .set_atom_data(first.atom_data().clone())
        .expect("valid atom data");
    first_member
        .props_mut()
        .insert("member".to_owned(), PropValue::Int(1));
    let mut second_member = EnsembleMember::new(second.positions().clone());
    second_member.set_cell(second.cell().copied());
    second_member.set_weight(Some(0.75)).expect("valid weight");
    second_member
        .set_atom_data(second.atom_data().clone())
        .expect("valid atom data");
    second_member
        .props_mut()
        .insert("member".to_owned(), PropValue::Int(2));
    let ensemble =
        Ensemble::from_members(Arc::clone(&fixture.topology), [first_member, second_member])
            .expect("valid ensemble");
    let source_weights = ensemble
        .members()
        .map(EnsembleMember::weight)
        .collect::<Vec<_>>();
    let source_props = ensemble
        .members()
        .map(|member| member.props().clone())
        .collect::<Vec<_>>();
    let remapped = ensemble
        .remap_to(&target_topology, edit.mapping())
        .expect("ensemble remap");
    assert!(Arc::ptr_eq(&remapped.shared_topology(), &target_topology));
    assert_eq!(
        remapped
            .members()
            .map(EnsembleMember::weight)
            .collect::<Vec<_>>(),
        source_weights
    );
    assert_eq!(
        remapped
            .members()
            .map(|member| member.props().clone())
            .collect::<Vec<_>>(),
        source_props
    );
    assert_ne!(
        remapped.member(0).expect("first member").cell(),
        remapped.member(1).expect("second member").cell()
    );
    for member in remapped.members() {
        assert!(member.positions().is_compatible(&target_topology));
        assert!(member.atom_data().is_compatible(&target_topology));
    }

    let (added_source, added_target, added_mapping, added_atom) = mapping_with_added_atom();
    let added_positions = Positions::new(
        &added_source,
        Quantity::new(vec![Point3::new(1.0, 0.0, 0.0)], ANGSTROM),
    )
    .expect("added-map source positions");
    let added_ensemble =
        Ensemble::from_members(added_source, [EnsembleMember::new(added_positions)])
            .expect("added-map source ensemble");
    assert!(matches!(
        added_ensemble.remap_to(&added_target, &added_mapping),
        Err(TopologyRemapError::Member { member: 0, error })
            if *error == TopologyRemapError::AddedAtomsRequireState {
                target_atom: added_atom,
            }
    ));
}

fn frame(fixture: &Fixture, offset: f64, step: u64) -> TrajectoryFrame {
    let model = model(fixture, offset);
    let vectors = (0..fixture.topology.atom_count())
        .map(|index| Vector3::new(offset + index as f64, 2.0, 3.0))
        .collect::<Vec<_>>();
    let mut frame = TrajectoryFrame::new(model.positions().clone());
    frame.set_cell(model.cell().copied());
    frame
        .set_velocities(Some(
            Velocities::new(
                &fixture.topology,
                Quantity::new(vectors.clone(), MODEL_VELOCITY_UNIT),
            )
            .expect("velocities"),
        ))
        .expect("compatible velocities");
    frame
        .set_forces(Some(
            Forces::new(&fixture.topology, Quantity::new(vectors, MODEL_FORCE_UNIT))
                .expect("forces"),
        ))
        .expect("compatible forces");
    frame
        .set_time(Some(Quantity::new(offset, PICOSECOND)))
        .expect("finite time");
    frame.set_step(Some(step));
    frame
        .set_atom_data(model.atom_data().clone())
        .expect("compatible atom data");
    frame
        .props_mut()
        .insert("frame".to_owned(), PropValue::Int(step as i64));
    frame
}

#[test]
fn trajectory_and_reusable_buffer_remapping_preserve_every_frame_field() {
    let fixture = fixture();
    let edit = retained_fixture(&fixture);
    let target_topology = edit.shared_topology();
    let positions_only = TrajectoryFrame::new(
        Positions::new(
            &fixture.topology,
            Quantity::new(point_values(&fixture.topology, 6.0), ANGSTROM),
        )
        .expect("positions-only frame"),
    );
    let frames = [
        frame(&fixture, 4.0, 40),
        positions_only,
        frame(&fixture, 8.0, 80),
    ];
    let trajectory = Trajectory::from_frames(Arc::clone(&fixture.topology), frames.clone())
        .expect("valid trajectory");
    let remapped = trajectory
        .remap_to(&target_topology, edit.mapping())
        .expect("trajectory remap");
    assert!(Arc::ptr_eq(&remapped.shared_topology(), &target_topology));
    assert_eq!(remapped.len(), 3);
    for (source_frame, target_frame) in trajectory.frames().zip(remapped.frames()) {
        assert_eq!(target_frame.cell(), source_frame.cell());
        assert_eq!(target_frame.time(), source_frame.time());
        assert_eq!(target_frame.step(), source_frame.step());
        assert_eq!(target_frame.props(), source_frame.props());
        assert_eq!(
            target_frame.atom_data().is_empty(),
            source_frame.atom_data().is_empty()
        );
        for (source_index, target_index) in edit.mapping().atom_index_pairs() {
            assert_eq!(
                target_frame.positions().values().value()[target_index.index()],
                source_frame.positions().values().value()[source_index.index()]
            );
            assert_eq!(
                target_frame
                    .velocities()
                    .map(|values| values.value()[target_index.index()]),
                source_frame
                    .velocities()
                    .map(|values| values.value()[source_index.index()])
            );
            assert_eq!(
                target_frame
                    .forces()
                    .map(|values| values.value()[target_index.index()]),
                source_frame
                    .forces()
                    .map(|values| values.value()[source_index.index()])
            );
        }
    }
    assert!(remapped
        .frame(1)
        .expect("positions-only frame")
        .velocities()
        .is_none());
    assert!(remapped
        .frame(1)
        .expect("positions-only frame")
        .forces()
        .is_none());

    let mut buffer = FrameBuffer::new(Arc::clone(&target_topology));
    buffer
        .copy_remapped_from(
            frames[0]
                .view(&fixture.topology)
                .expect("borrowed source frame"),
            edit.mapping(),
        )
        .expect("buffer remap");
    let buffer_view = buffer.frame_view();
    assert_eq!(buffer_view.time(), frames[0].time());
    assert_eq!(buffer_view.step(), frames[0].step());
    assert_eq!(buffer_view.props(), frames[0].props());
    assert!(Arc::ptr_eq(
        &buffer.model_view().shared_topology(),
        &target_topology
    ));
    let source_bond = InstanceBondId::new(fixture.ligand, BondId::new(1));
    let target_bond = edit
        .mapping()
        .map_bond(source_bond)
        .expect("retained ligand bond");
    let mut potential = HarmonicBondPotential::new(
        &target_topology,
        [HarmonicBondParameter::new(
            target_bond,
            Quantity::new(1.2, ANGSTROM),
            Quantity::new(100.0, MODEL_FORCE_CONSTANT_UNIT),
        )],
    )
    .expect("prepare target potential");
    buffer.set_cell(None);
    assert!(potential
        .evaluate(buffer.frame_view().model_view())
        .expect("evaluate remapped target view")
        .energy()
        .to_value()
        .is_finite());

    let independent = self::fixture();
    let independent_edit = retained_fixture(&independent);
    let mut wrong_buffer = FrameBuffer::new(independent_edit.shared_topology());
    let before_positions = wrong_buffer.positions().values().value().to_vec();
    assert_eq!(
        wrong_buffer.copy_remapped_from(
            frames[0]
                .view(&fixture.topology)
                .expect("borrowed source frame"),
            edit.mapping(),
        ),
        Err(TrajectoryRemapError::IncompatibleDestinationBuffer)
    );
    assert_eq!(
        *wrong_buffer.positions().values().value(),
        before_positions.as_slice()
    );

    let (added_source, added_target, added_mapping, added_atom) = mapping_with_added_atom();
    let added_positions = Positions::new(
        &added_source,
        Quantity::new(vec![Point3::new(2.0, 0.0, 0.0)], ANGSTROM),
    )
    .expect("added-map source positions");
    let added_trajectory =
        Trajectory::from_frames(added_source, [TrajectoryFrame::new(added_positions)])
            .expect("added-map source trajectory");
    assert!(matches!(
        added_trajectory.remap_to(&added_target, &added_mapping),
        Err(TrajectoryRemapError::Frame { frame: 0, error })
            if *error == TrajectoryRemapError::Topology(
                TopologyRemapError::AddedAtomsRequireState {
                    target_atom: added_atom,
                }
            )
    ));
}

#[test]
fn reusable_buffer_remapping_is_transactional_and_clears_stale_state() {
    let fixture = fixture();
    let edit = retained_fixture(&fixture);
    let target_topology = edit.shared_topology();
    let full = frame(&fixture, 4.0, 40);
    let positions_only = TrajectoryFrame::new(
        Positions::new(
            &fixture.topology,
            Quantity::new(point_values(&fixture.topology, 6.0), ANGSTROM),
        )
        .expect("positions-only frame"),
    );
    let mut buffer = FrameBuffer::new(target_topology);
    buffer
        .copy_remapped_from(
            full.view(&fixture.topology).expect("full borrowed frame"),
            edit.mapping(),
        )
        .expect("full buffer remap");
    assert!(buffer.frame_view().cell().is_some());
    assert!(buffer.frame_view().velocities().is_some());
    assert!(buffer.frame_view().forces().is_some());
    assert!(buffer.frame_view().time().is_some());
    assert!(buffer.frame_view().step().is_some());
    assert!(!buffer.frame_view().atom_data().is_empty());
    assert!(!buffer.frame_view().props().is_empty());

    buffer
        .copy_remapped_from(
            positions_only
                .view(&fixture.topology)
                .expect("positions-only borrowed frame"),
            edit.mapping(),
        )
        .expect("positions-only buffer remap");
    let cleared = buffer.frame_view();
    assert_eq!(cleared.cell(), None);
    assert_eq!(cleared.velocities(), None);
    assert_eq!(cleared.forces(), None);
    assert_eq!(cleared.time(), None);
    assert_eq!(cleared.step(), None);
    assert!(cleared.atom_data().is_empty());
    assert!(cleared.props().is_empty());
    for (source_index, target_index) in edit.mapping().atom_index_pairs() {
        assert_eq!(
            cleared.positions().values().value()[target_index.index()],
            positions_only.positions().values().value()[source_index.index()]
        );
    }

    let (source, target, mapping, added_atom) = mapping_with_added_atom();
    let source_frame = TrajectoryFrame::new(
        Positions::new(
            &source,
            Quantity::new(vec![Point3::new(1.0, 2.0, 3.0)], ANGSTROM),
        )
        .expect("source positions"),
    );
    let mut destination = FrameBuffer::new(Arc::clone(&target));
    destination
        .set_positions(Quantity::new(
            vec![Point3::new(11.0, 12.0, 13.0), Point3::new(21.0, 22.0, 23.0)],
            ANGSTROM,
        ))
        .expect("destination positions");
    let destination_cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(31.0, 32.0, 33.0), ANGSTROM),
        [true, false, true],
    )
    .expect("destination cell");
    destination.set_cell(Some(destination_cell));
    destination
        .set_velocities(Some(Quantity::new(
            vec![
                Vector3::new(41.0, 42.0, 43.0),
                Vector3::new(51.0, 52.0, 53.0),
            ],
            MODEL_VELOCITY_UNIT,
        )))
        .expect("destination velocities");
    destination
        .set_forces(Some(Quantity::new(
            vec![
                Vector3::new(61.0, 62.0, 63.0),
                Vector3::new(71.0, 72.0, 73.0),
            ],
            MODEL_FORCE_UNIT,
        )))
        .expect("destination forces");
    destination
        .set_time(Some(Quantity::new(81.0, PICOSECOND)))
        .expect("destination time");
    destination.set_step(Some(91));
    let mut destination_atom_data = AtomData::new(&target);
    destination_atom_data
        .set_b_factor(
            &target,
            target.atom_ids()[0],
            Some(Quantity::new(99.0, SQUARE_ANGSTROM)),
        )
        .expect("destination atom data");
    destination
        .set_atom_data(destination_atom_data)
        .expect("destination atom data");
    destination.props_mut().insert(
        "frame".to_owned(),
        PropValue::String("destination".to_owned()),
    );
    let before = destination.clone();

    assert_eq!(
        destination.copy_remapped_from(
            source_frame.view(&source).expect("borrowed source frame"),
            &mapping,
        ),
        Err(TrajectoryRemapError::Topology(
            TopologyRemapError::AddedAtomsRequireState {
                target_atom: added_atom,
            }
        ))
    );

    let before = before.frame_view();
    let after = destination.frame_view();
    assert!(Arc::ptr_eq(
        &after.shared_topology(),
        &before.shared_topology()
    ));
    assert_eq!(after.positions(), before.positions());
    assert_eq!(after.cell(), before.cell());
    assert_eq!(after.velocities(), before.velocities());
    assert_eq!(after.forces(), before.forces());
    assert_eq!(after.time(), before.time());
    assert_eq!(after.step(), before.step());
    assert_eq!(after.atom_data(), before.atom_data());
    assert_eq!(after.props(), before.props());
}

#[test]
fn solvent_rich_subset_regression_avoids_quadratic_builder_cloning() {
    const WATER_COUNT: usize = 20_000;
    let water = one_atom_small("O");
    let ligand = one_atom_small("C");
    let mut builder = TopologyBuilder::new();
    let water_definition = builder
        .add_molecule_definition(&water)
        .expect("water definition");
    let ligand_definition = builder
        .add_molecule_definition(&ligand)
        .expect("ligand definition");
    builder
        .reserve_instances(WATER_COUNT + 1)
        .expect("synthetic capacity");
    let ligand_instance = builder
        .add_instance(ligand_definition, metadata(MoleculeRole::Ligand, "ligand"))
        .expect("ligand instance");
    let mut retained_waters = Vec::new();
    for index in 0..WATER_COUNT {
        let water_instance = builder
            .add_instance(
                water_definition,
                metadata(MoleculeRole::Solvent, &format!("water-{index}")),
            )
            .expect("water instance");
        if index % 2 == 0 {
            retained_waters.push(water_instance);
        }
    }
    let topology = Arc::new(builder.build().expect("synthetic solvent topology"));
    let requested = std::iter::once(ligand_instance).chain(retained_waters.iter().copied());
    let edit = retain_instances(&topology, requested).expect("large solvent-rich subset");
    assert_eq!(edit.topology().definition_count(), 2);
    assert_eq!(edit.topology().instance_count(), WATER_COUNT / 2 + 1);
    assert_eq!(edit.mapping().instance_pairs().len(), WATER_COUNT / 2 + 1);
    assert_eq!(edit.mapping().atom_pairs().len(), WATER_COUNT / 2 + 1);
}
use std::sync::Arc;
