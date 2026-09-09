use std::sync::Arc;

use kekule::core::{
    Atom, AxisStereo, Element, Molecule, MoleculeEditor, MoleculeError, StereoCarrier,
    StereoElement, StereoElementKind, StereoGroup, StereoGroupKind,
};
use kekule::topology::{
    transform, AtomSiteMetadata, InstanceAtomId, InstanceBondId, MoleculeClass,
    MoleculeDefinitionId, ResidueClass, Topology, TopologyBuilder,
};

fn molecule(text: &str) -> Molecule {
    kekule::smiles::to_molecules(text).unwrap().pop().unwrap()
}

#[test]
fn ring_bond_deletion_prunes_invalid_carriers_and_preserves_unrelated_stereo_groups() {
    let mut source = molecule("F[C@]1(Cl)CCNC1C[C@H](Br)I").into_editor();
    let elements = source.stereo_element_ids().collect::<Vec<_>>();
    assert_eq!(elements.len(), 2);
    let group = source
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::And,
            members: elements.clone(),
        })
        .unwrap();
    let source = source.finish().unwrap();
    let StereoElementKind::Tetrahedral(stereo) = &source.stereo_element(elements[0]).unwrap().kind
    else {
        panic!("expected tetrahedral center");
    };
    let carrier = stereo
        .carriers
        .iter()
        .find_map(|carrier| match carrier {
            StereoCarrier::Atom(atom) if source.neighbors(*atom).unwrap().count() > 1 => {
                Some(*atom)
            }
            _ => None,
        })
        .unwrap();
    let removed = source
        .bond_between(stereo.center, carrier)
        .unwrap()
        .unwrap();
    for batch in [false, true] {
        let mut editor = source.edit();
        if batch {
            editor.delete_bonds([removed]).unwrap();
        } else {
            editor.delete_bond(removed).unwrap();
        }
        assert!(editor.stereo_element(elements[0]).is_err());
        assert_eq!(
            editor.stereo_element(elements[1]).unwrap(),
            source.stereo_element(elements[1]).unwrap()
        );
        assert_eq!(
            editor.stereo_group(group).unwrap().members,
            vec![elements[1]]
        );
        editor.validate().unwrap();
        let result = editor.try_finish().unwrap();
        // Enhanced groups have a separate export capability; test carrier export
        // on an ungrouped copy after asserting the preserved group above.
        let mut ungrouped = result.edit();
        ungrouped.remove_stereo_group(group).unwrap();
        let text = kekule::smiles::write_isomeric(&ungrouped.finish().unwrap()).unwrap();
        assert_eq!(molecule(&text).stereo_elements().count(), 1);
    }
    let topology = Arc::new(Topology::from_molecule(&source).unwrap());
    let instance = topology.instances().next().unwrap().0;
    let mut editor = topology.edit();
    let bond = editor
        .bond_handle(InstanceBondId::new(instance, removed))
        .unwrap();
    editor.delete_bond(bond).unwrap();
    let result = editor.finish().unwrap();
    let result = result.molecules().next().unwrap();
    assert_eq!(result.molecule().stereo_elements().count(), 1);
    let mut ungrouped = result.molecule().edit();
    for group in result.molecule().stereo_groups().map(|(id, _)| id) {
        ungrouped.remove_stereo_group(group).unwrap();
    }
    kekule::smiles::write_isomeric(&ungrouped.finish().unwrap()).unwrap();
    assert_eq!(source.stereo_elements().count(), 2);
}

#[test]
fn checked_stereo_boundaries_reject_nonadjacent_carriers_transactionally() {
    let mut axis = molecule("FCCF").into_editor();
    let atoms = axis.atom_ids().collect::<Vec<_>>();
    axis.add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
        axis: axis.bond_between(atoms[1], atoms[2]).unwrap().unwrap(),
        carriers: vec![StereoCarrier::Atom(atoms[0]), StereoCarrier::Atom(atoms[3])],
        orientation: None,
    })))
    .unwrap();
    for source in [
        molecule("F[C@]1(Cl)CCNC1"),
        molecule("F/C=C/F"),
        axis.finish().unwrap(),
    ] {
        let (id, original) = source.stereo_elements().next().unwrap();
        let mut invalid = original.clone();
        match &mut invalid.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                stereo.carriers[0] = StereoCarrier::Atom(stereo.center)
            }
            StereoElementKind::DoubleBond(stereo) => {
                stereo.left_carrier = StereoCarrier::Atom(stereo.right)
            }
            StereoElementKind::Axis(stereo) => {
                stereo.carriers[0] = StereoCarrier::Atom(source.bond(stereo.axis).unwrap().a())
            }
        }
        let mut editor = source.edit();
        let before = format!("{editor:?}");
        assert!(matches!(
            editor.add_stereo_element(invalid.clone()),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        assert!(matches!(
            editor.replace_stereo_element(id, invalid),
            Err(MoleculeError::InvalidStereoReference(_))
        ));
        assert_eq!(format!("{editor:?}"), before);
        assert_eq!(editor.finish().unwrap(), source);
    }
}

fn oxygen_builder() -> (TopologyBuilder, MoleculeDefinitionId) {
    let mut editor = MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder
        .add_molecule_definition_owned(editor.finish().unwrap())
        .unwrap();
    (builder, definition)
}

fn add_instance(
    builder: &mut TopologyBuilder,
    definition: MoleculeDefinitionId,
    component: Option<&str>,
    residue_class: Option<ResidueClass>,
) -> kekule::topology::MoleculeInstanceId {
    let local = builder
        .definition(definition)
        .unwrap()
        .molecule()
        .atom_ids()
        .next()
        .unwrap();
    let instance = builder.add_instance(definition).unwrap();
    if let Some(component) = component {
        let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, component, None, None, None)
            .unwrap();
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, local),
                AtomSiteMetadata::default(),
            )
            .unwrap();
        if let Some(class) = residue_class {
            builder.set_residue_class(residue, class).unwrap();
        }
    }
    instance
}

#[test]
fn resumed_reused_definitions_combine_new_informative_evidence_like_fresh_builders() {
    for explicit in [false, true] {
        let (mut builder, definition) = oxygen_builder();
        add_instance(&mut builder, definition, Some("HOH"), None);
        if explicit {
            builder
                .set_molecule_class(definition, MoleculeClass::Water)
                .unwrap();
        }
        let source = builder.clone().build().unwrap();
        assert_eq!(
            source.definition(definition).unwrap().class(),
            MoleculeClass::Water
        );
        for mut builder in [source.into_builder(), builder] {
            add_instance(
                &mut builder,
                definition,
                Some("UNK"),
                Some(ResidueClass::Ion),
            );
            builder.validate().unwrap();
            let result = builder.try_build().unwrap();
            assert_eq!(
                result.definition(definition).unwrap().class(),
                if explicit {
                    MoleculeClass::Water
                } else {
                    MoleculeClass::Other
                }
            );
            assert_eq!(
                result.residues().map(|r| r.class()).collect::<Vec<_>>(),
                vec![ResidueClass::Water, ResidueClass::Ion]
            );
        }
    }
    let (mut builder, definition) = oxygen_builder();
    add_instance(&mut builder, definition, None, None);
    let mut builder = builder.build().unwrap().into_builder();
    add_instance(&mut builder, definition, Some("HOH"), None);
    assert_eq!(
        builder
            .build()
            .unwrap()
            .definition(definition)
            .unwrap()
            .class(),
        MoleculeClass::Water
    );
}

#[test]
fn uninformative_and_unrelated_appends_preserve_complete_entity_cached_classes() {
    let (mut builder, definition) = oxygen_builder();
    let retained = add_instance(&mut builder, definition, Some("HOH"), None);
    add_instance(
        &mut builder,
        definition,
        Some("UNK"),
        Some(ResidueClass::Ion),
    );
    let source = Arc::new(builder.build().unwrap());
    // The retained class deliberately remembers the complete source definition.
    // Adding an uninformative occurrence must not recompute it from just HOH.
    for component in [None, Some("UNK")] {
        let retained = transform::retain_instances(&source, [retained]).unwrap();
        assert_eq!(
            retained.definition(definition).unwrap().class(),
            MoleculeClass::Other
        );
        let mut builder = Arc::try_unwrap(retained).unwrap().into_builder();
        add_instance(&mut builder, definition, component, None);
        let independent = builder.add_molecule_definition(&molecule("CC")).unwrap();
        builder.add_instance(independent).unwrap();
        let result = builder.build().unwrap();
        assert_eq!(
            result.definition(definition).unwrap().class(),
            MoleculeClass::Other
        );
        assert_eq!(
            result.definition(independent).unwrap().class(),
            MoleculeClass::SmallMolecule
        );
    }
}
