use molecular::bio::{
    MacroMolecule, SmcraAtomSiteMetadata, SmcraChainId, SmcraHierarchy, SmcraResidueId,
};
use molecular::core::{
    AromaticityModel, AromaticityProvenance, Atom, AtomId, BondOrder, Element, Molecule,
    PerceptionState, PerceptionStateBuildError, PerceptionStateInstallError, PropValue, Ring,
    RingMembership, RingSet, RingWork, StereoCarrier, StereoDescriptor, StereoElement,
    StereoElementId, StereoElementKind, StereoGroup, StereoGroupKind, StereoSource,
    TetrahedralOrientation, TetrahedralStereo, ValenceModel,
};
use molecular::small::SmallMolecule;
use molecular::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};

fn carbon() -> Atom {
    Atom::new(Element::from_symbol("C").expect("carbon is known"))
}

fn export_perception(
    source: &PerceptionState,
) -> Result<PerceptionState, PerceptionStateBuildError> {
    let mut builder = PerceptionState::builder();
    if let Some(valence) = source.valence_state() {
        builder = builder.with_valence(valence.model(), valence.implicit_hydrogens().collect())?;
    }
    if let Some(rings) = source.ring_state() {
        let membership = RingMembership::from_slot_flags(
            rings.membership().atom_slot_flags().to_vec(),
            rings.membership().bond_slot_flags().to_vec(),
        );
        let ring_set = rings
            .ring_set()
            .map(|set| RingSet::from_parts(set.rings().to_vec(), set.work()));
        builder = builder.with_rings(membership, ring_set);
    }
    if let Some(aromaticity) = source.aromaticity_state() {
        builder = builder.with_aromaticity(
            aromaticity.provenance(),
            aromaticity.atoms().collect(),
            aromaticity.bonds().collect(),
        )?;
    }
    builder = builder.with_cip_descriptors(source.cip_descriptors().collect())?;
    Ok(builder.build())
}

fn build_small_topology(molecule: &SmallMolecule) -> Topology {
    let mut builder = TopologyBuilder::new();
    let definition = builder
        .add_small_molecule_definition(molecule)
        .expect("definition");
    builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .expect("instance");
    builder.build().expect("topology")
}

fn tetrahedral_element(center: AtomId, carriers: [AtomId; 4]) -> StereoElement {
    StereoElement::specified(
        StereoElementKind::Tetrahedral(TetrahedralStereo {
            center,
            carriers: carriers.into_iter().map(StereoCarrier::Atom).collect(),
            orientation: TetrahedralOrientation::Clockwise,
        }),
        StereoSource::User,
    )
}

fn stereo_fixture() -> (Molecule, Vec<StereoElementId>) {
    let mut molecule = Molecule::new();
    let center = molecule.add_atom(carbon()).expect("center");
    let carriers = [
        molecule.add_atom(carbon()).expect("carrier"),
        molecule.add_atom(carbon()).expect("carrier"),
        molecule.add_atom(carbon()).expect("carrier"),
        molecule.add_atom(carbon()).expect("carrier"),
    ];
    let elements = (0..4)
        .map(|_| {
            molecule
                .add_stereo_element(tetrahedral_element(center, carriers))
                .expect("stereo element")
        })
        .collect();
    (molecule, elements)
}

#[test]
fn full_installed_perception_round_trips_through_public_api() {
    let mut source =
        SmallMolecule::from_smiles_sanitized("c1ccccc1[C@H](F)Cl").expect("sanitized source");
    molecular::perception::stereo::assign_cip_descriptors(source.graph_mut());
    let perception = source.graph().perception();
    assert_eq!(
        perception.valence_state().and_then(|state| state.model()),
        Some(ValenceModel::RdkitLike)
    );
    assert!(perception
        .ring_state()
        .is_some_and(|state| state.ring_set().is_some()));
    assert_eq!(
        perception
            .aromaticity_state()
            .map(|state| state.provenance()),
        Some(AromaticityProvenance::Perceived(
            AromaticityModel::RdkitLike
        ))
    );
    assert!(perception.cip_descriptors().next().is_some());

    let detached = export_perception(perception).expect("public export");
    assert_eq!(&detached, perception);

    let mut reconstructed_graph = source.graph().clone();
    reconstructed_graph.invalidate_topology();
    reconstructed_graph
        .install_perception_state(detached)
        .expect("public install");
    assert_eq!(
        reconstructed_graph.perception(),
        source.graph().perception()
    );

    let reconstructed = SmallMolecule::from_graph(reconstructed_graph);
    assert!(build_small_topology(&source).same_layout(&build_small_topology(&reconstructed)));
}

#[test]
fn exact_section_presence_and_provenance_round_trip() {
    let mut molecule = Molecule::new();
    let atom = molecule.add_atom(carbon()).expect("atom");

    let model_neutral = PerceptionState::builder()
        .with_valence(None, vec![(atom, 4)])
        .expect("valence")
        .build();
    molecule
        .install_perception_state(model_neutral.clone())
        .expect("install");
    assert!(!molecule.perception().has_valence());
    assert_eq!(
        molecule
            .perception()
            .valence_state()
            .expect("present valence")
            .model(),
        None
    );
    assert_eq!(
        export_perception(molecule.perception()).unwrap(),
        model_neutral
    );

    for provenance in [
        AromaticityProvenance::Imported,
        AromaticityProvenance::Perceived(AromaticityModel::RdkitLike),
    ] {
        let state = PerceptionState::builder()
            .with_aromaticity(provenance, vec![atom], Vec::new())
            .expect("aromaticity")
            .build();
        molecule
            .install_perception_state(state.clone())
            .expect("install");
        assert_eq!(export_perception(molecule.perception()).unwrap(), state);
    }
}

#[test]
fn ring_membership_round_trips_with_and_without_ring_set() {
    let mut molecule = Molecule::new();
    let atoms = [
        molecule.add_atom(carbon()).unwrap(),
        molecule.add_atom(carbon()).unwrap(),
        molecule.add_atom(carbon()).unwrap(),
    ];
    let bonds = [
        molecule
            .add_bond(atoms[0], atoms[1], BondOrder::Single)
            .unwrap(),
        molecule
            .add_bond(atoms[1], atoms[2], BondOrder::Single)
            .unwrap(),
        molecule
            .add_bond(atoms[2], atoms[0], BondOrder::Single)
            .unwrap(),
    ];
    let membership = || RingMembership::from_slot_flags(vec![true; 3], vec![true; 3]);

    let membership_only = PerceptionState::builder()
        .with_rings(membership(), None)
        .build();
    molecule
        .install_perception_state(membership_only.clone())
        .expect("membership-only install");
    assert!(molecule
        .perception()
        .ring_state()
        .is_some_and(|state| state.ring_set().is_none()));
    assert_eq!(
        export_perception(molecule.perception()).unwrap(),
        membership_only
    );

    let work = RingWork {
        atom_count: 3,
        bond_count: 3,
        candidate_cycles: 1,
        total_work: 7,
        ..RingWork::default()
    };
    let with_basis = PerceptionState::builder()
        .with_rings(
            membership(),
            Some(RingSet::from_parts(
                vec![Ring {
                    atoms: atoms.to_vec(),
                    bonds: bonds.to_vec(),
                }],
                work,
            )),
        )
        .build();
    molecule
        .install_perception_state(with_basis.clone())
        .expect("basis install");
    assert_eq!(
        export_perception(molecule.perception()).unwrap(),
        with_basis
    );

    let inconsistent = PerceptionState::builder()
        .with_rings(
            RingMembership::from_slot_flags(vec![false, true, true], vec![true; 3]),
            Some(RingSet::from_parts(
                vec![Ring {
                    atoms: atoms.to_vec(),
                    bonds: bonds.to_vec(),
                }],
                work,
            )),
        )
        .build();
    assert_eq!(
        molecule.install_perception_state(inconsistent),
        Err(PerceptionStateInstallError::InconsistentRingAtomMembership(
            atoms[0]
        ))
    );
    assert_eq!(molecule.perception(), &with_basis);
}

#[test]
fn malformed_perception_is_rejected_transactionally() {
    let mut molecule = Molecule::new();
    let live = molecule.add_atom(carbon()).unwrap();
    let deleted = molecule.add_atom(carbon()).unwrap();
    molecule.delete_atom(deleted).unwrap();
    let previous = PerceptionState::builder()
        .with_valence(Some(ValenceModel::RdkitLike), vec![(live, 4)])
        .unwrap()
        .build();
    molecule
        .install_perception_state(previous.clone())
        .expect("baseline");

    let invalid_atom = PerceptionState::builder()
        .with_valence(None, vec![(deleted, 0)])
        .unwrap()
        .build();
    assert_eq!(
        molecule.install_perception_state(invalid_atom),
        Err(PerceptionStateInstallError::InvalidAtomId(deleted))
    );
    assert_eq!(molecule.perception(), &previous);

    let bad_slots = PerceptionState::builder()
        .with_rings(
            RingMembership::from_slot_flags(vec![false], Vec::new()),
            None,
        )
        .build();
    assert!(matches!(
        molecule.install_perception_state(bad_slots),
        Err(PerceptionStateInstallError::RingAtomSlotCountMismatch { .. })
    ));
    assert_eq!(molecule.perception(), &previous);

    assert!(matches!(
        PerceptionState::builder().with_aromaticity(
            AromaticityProvenance::Imported,
            vec![live, live],
            Vec::new()
        ),
        Err(PerceptionStateBuildError::DuplicateAromaticAtom(id)) if id == live
    ));
}

#[test]
fn malformed_ring_and_stereo_references_are_rejected() {
    let (mut molecule, elements) = stereo_fixture();
    let deleted_bond = molecule
        .add_bond(AtomId::new(0), AtomId::new(1), BondOrder::Single)
        .unwrap();
    molecule.delete_bond(deleted_bond).unwrap();
    let removed = elements[3];
    molecule.remove_stereo_element(removed).unwrap();
    let previous = molecule.perception().clone();

    let invalid_bond = PerceptionState::builder()
        .with_aromaticity(
            AromaticityProvenance::Imported,
            Vec::new(),
            vec![deleted_bond],
        )
        .unwrap()
        .build();
    assert_eq!(
        molecule.install_perception_state(invalid_bond),
        Err(PerceptionStateInstallError::InvalidBondId(deleted_bond))
    );
    assert_eq!(molecule.perception(), &previous);

    let invalid_cip = PerceptionState::builder()
        .with_cip_descriptors(vec![(removed, StereoDescriptor::R)])
        .unwrap()
        .build();
    assert_eq!(
        molecule.install_perception_state(invalid_cip),
        Err(PerceptionStateInstallError::InvalidStereoElementId(removed))
    );
    assert_eq!(molecule.perception(), &previous);

    let membership =
        RingMembership::from_slot_flags(vec![true; molecule.atom_count()], vec![false]);
    let malformed = PerceptionState::builder()
        .with_rings(
            membership,
            Some(RingSet::from_parts(
                vec![Ring {
                    atoms: vec![AtomId::new(0), AtomId::new(1)],
                    bonds: Vec::new(),
                }],
                RingWork {
                    atom_count: molecule.atom_count(),
                    bond_count: 0,
                    ..RingWork::default()
                },
            )),
        )
        .build();
    assert!(matches!(
        molecule.install_perception_state(malformed),
        Err(PerceptionStateInstallError::MalformedRing { .. })
    ));
    assert_eq!(molecule.perception(), &previous);
}

#[test]
fn installed_perception_follows_normal_invalidation_rules() {
    let (mut molecule, elements) = stereo_fixture();
    let state = PerceptionState::builder()
        .with_valence(None, vec![(AtomId::new(0), 0)])
        .unwrap()
        .with_cip_descriptors(vec![(elements[0], StereoDescriptor::R)])
        .unwrap()
        .build();
    molecule
        .install_perception_state(state.clone())
        .expect("install");
    molecule
        .atom_mut(AtomId::new(0))
        .unwrap()
        .props
        .insert("note".into(), PropValue::Bool(true));
    assert_eq!(molecule.perception(), &state);

    let mut replacement = molecule.stereo_element(elements[0]).unwrap().clone();
    replacement.source = StereoSource::Reaction;
    molecule
        .replace_stereo_element(elements[0], replacement)
        .unwrap();
    assert!(molecule.perception().valence_state().is_some());
    assert!(!molecule.perception().has_cip_descriptors());

    molecule.atom_mut(AtomId::new(0)).unwrap().formal_charge = 1;
    assert_eq!(molecule.perception(), &PerceptionState::default());
}

#[test]
fn enriched_smcra_state_round_trips_and_invalid_ids_are_transactional() {
    let mut graph = Molecule::new();
    let atom = graph.add_atom(carbon()).unwrap();
    let mut source = SmcraHierarchy::new();
    let chain = source.add_chain("A", Some("AUTH_A".into())).unwrap();
    let residue = source
        .add_residue(
            chain,
            "display-name",
            Some(7),
            Some("42".into()),
            Some("I".into()),
        )
        .unwrap();
    source
        .set_residue_component_ids(
            residue,
            Some("LABEL_COMPONENT".into()),
            Some("AUTHOR_COMPONENT".into()),
        )
        .unwrap();
    let site = source
        .add_atom_site(
            residue,
            atom,
            SmcraAtomSiteMetadata {
                type_symbol: Some("C".into()),
                label_asym_id: Some("A".into()),
                auth_asym_id: Some("AUTH_A".into()),
                label_atom_id: Some("CA".into()),
                auth_atom_id: Some("C-alpha".into()),
            },
        )
        .unwrap();
    source
        .chain_props_mut(chain)
        .unwrap()
        .insert("chain".into(), PropValue::String("value".into()));
    source
        .residue_props_mut(residue)
        .unwrap()
        .insert("residue".into(), PropValue::Int(7));
    source
        .atom_site_props_mut(site)
        .unwrap()
        .insert("site".into(), PropValue::Float(0.5));

    let mut rebuilt = SmcraHierarchy::new();
    for (id, chain_record) in source.chains() {
        let added = rebuilt
            .add_chain(
                chain_record.label_id(),
                chain_record.author_id().map(str::to_owned),
            )
            .unwrap();
        assert_eq!(added, id);
        rebuilt
            .chain_props_mut(added)
            .unwrap()
            .clone_from(chain_record.props());
    }
    for (id, residue_record) in source.residues() {
        let added = rebuilt
            .add_residue(
                residue_record.chain(),
                residue_record.name(),
                residue_record.label_seq_id(),
                residue_record.author_seq_id().map(str::to_owned),
                residue_record.insertion_code().map(str::to_owned),
            )
            .unwrap();
        assert_eq!(added, id);
        rebuilt
            .set_residue_component_ids(
                added,
                residue_record.label_comp_id().map(str::to_owned),
                residue_record.author_comp_id().map(str::to_owned),
            )
            .unwrap();
        rebuilt
            .residue_props_mut(added)
            .unwrap()
            .clone_from(residue_record.props());
    }
    for (id, site_record) in source.atom_sites() {
        let added = rebuilt
            .add_atom_site(
                site_record.residue(),
                site_record.atom(),
                site_record.metadata().clone(),
            )
            .unwrap();
        assert_eq!(added, id);
        rebuilt
            .atom_site_props_mut(added)
            .unwrap()
            .clone_from(site_record.props());
    }
    assert_eq!(rebuilt, source);
    let before = rebuilt.clone();
    assert!(rebuilt
        .set_residue_component_ids(SmcraResidueId::new(99), None, None)
        .is_err());
    assert!(rebuilt.chain_props_mut(SmcraChainId::new(99)).is_err());
    assert!(rebuilt
        .atom_site_props_mut(molecular::bio::SmcraAtomSiteId::new(99))
        .is_err());
    assert_eq!(rebuilt, before);

    let source_macro = MacroMolecule::try_from_parts(graph.clone(), source).unwrap();
    let rebuilt_macro = MacroMolecule::try_from_parts(graph, rebuilt).unwrap();
    assert_eq!(source_macro.hierarchy(), rebuilt_macro.hierarchy());

    let build = |molecule: &MacroMolecule| {
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_macro_molecule_definition(molecule).unwrap();
        builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        builder.build().unwrap()
    };
    assert!(build(&source_macro).same_layout(&build(&rebuilt_macro)));
}

#[test]
fn stereo_group_tombstone_layout_replays_with_stable_next_id() {
    let (mut source, elements) = stereo_fixture();
    let first = source
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: vec![elements[0]],
        })
        .unwrap();
    let removed = source
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Relative,
            members: vec![elements[1]],
        })
        .unwrap();
    source.remove_stereo_group(removed).unwrap();
    let third = source
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Or,
            members: vec![elements[2]],
        })
        .unwrap();
    let trailing = source.append_stereo_group_tombstone().unwrap();
    assert_eq!(
        (first.raw(), removed.raw(), third.raw(), trailing.raw()),
        (0, 1, 2, 3)
    );

    let (mut rebuilt, rebuilt_elements) = stereo_fixture();
    for (expected, slot) in source.stereo_group_slots() {
        let actual = if let Some(group) = slot {
            rebuilt.add_stereo_group(group.clone()).unwrap()
        } else {
            rebuilt.append_stereo_group_tombstone().unwrap()
        };
        assert_eq!(actual, expected);
    }
    assert_eq!(
        source
            .stereo_group_slots()
            .map(|(id, group)| (id, group.cloned()))
            .collect::<Vec<_>>(),
        rebuilt
            .stereo_group_slots()
            .map(|(id, group)| (id, group.cloned()))
            .collect::<Vec<_>>()
    );
    let source_topology = build_small_topology(&SmallMolecule::from_graph(source.clone()));
    let rebuilt_topology = build_small_topology(&SmallMolecule::from_graph(rebuilt.clone()));
    assert!(source_topology.same_layout(&rebuilt_topology));

    let source_next = source
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::And,
            members: vec![elements[1]],
        })
        .unwrap();
    let rebuilt_next = rebuilt
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::And,
            members: vec![rebuilt_elements[1]],
        })
        .unwrap();
    assert_eq!(source_next, rebuilt_next);
}

#[test]
fn stereo_group_pruning_tombstones_only_empty_groups_and_tombstone_append_preserves_cip() {
    let (mut molecule, elements) = stereo_fixture();
    let single = molecule
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: vec![elements[0]],
        })
        .unwrap();
    let multiple = molecule
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Relative,
            members: vec![elements[1], elements[2]],
        })
        .unwrap();
    let state = PerceptionState::builder()
        .with_cip_descriptors(vec![
            (elements[0], StereoDescriptor::R),
            (elements[1], StereoDescriptor::S),
        ])
        .unwrap()
        .build();
    molecule
        .install_perception_state(state.clone())
        .expect("install CIP");
    let tombstone = molecule.append_stereo_group_tombstone().unwrap();
    assert_eq!(molecule.perception(), &state);
    assert!(molecule
        .stereo_group_slots()
        .nth(tombstone.index())
        .unwrap()
        .1
        .is_none());

    molecule.remove_stereo_element(elements[0]).unwrap();
    assert!(molecule
        .stereo_group_slots()
        .nth(single.index())
        .unwrap()
        .1
        .is_none());
    molecule.remove_stereo_element(elements[1]).unwrap();
    assert_eq!(
        molecule.stereo_group(multiple).unwrap().members,
        vec![elements[2]]
    );
}

#[test]
fn empty_stereo_groups_remain_rejected() {
    let (mut molecule, _) = stereo_fixture();
    assert!(molecule
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: Vec::new(),
        })
        .is_err());
}
