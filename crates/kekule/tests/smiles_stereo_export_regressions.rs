use std::collections::BTreeMap;

use kekule::core::{
    DoubleBondOrientation, Molecule, MoleculeEditor, MoleculeError, StereoCarrier,
    StereoDescriptor, StereoElement, StereoElementKind,
};
use kekule::{smiles, stereo};

fn molecule(source: &str) -> Molecule {
    smiles::to_molecules(source).unwrap().pop().unwrap()
}

#[test]
fn every_conjugated_diene_and_triene_configuration_round_trips() {
    for double_bonds in [2, 3] {
        // Include both equivalent global slash gauges for every E/Z combination.
        for pattern in 0..(1 << (double_bonds + 1)) {
            let direction = |index| {
                if pattern & (1 << index) == 0 {
                    '/'
                } else {
                    '\\'
                }
            };
            let mut source = String::from("F");
            for index in 0..double_bonds {
                source.push(direction(index));
                source.push_str("C=C");
            }
            source.push(direction(double_bonds));
            source.push('F');
            for perceive in [false, true] {
                let mut original = molecule(&source);
                if perceive {
                    original.perceive().unwrap();
                }
                let written = smiles::write_isomeric(&original).unwrap();
                let restored = molecule(&written);
                assert_eq!(
                    original
                        .stereo_elements()
                        .map(|(_, s)| s)
                        .collect::<Vec<_>>(),
                    restored
                        .stereo_elements()
                        .map(|(_, s)| s)
                        .collect::<Vec<_>>(),
                    "{source} -> {written}"
                );
                assert_eq!(restored.stereo_elements().count(), double_bonds);
                assert_eq!(smiles::write_isomeric(&restored).unwrap(), written);
            }
        }
    }
}

fn mapped_polyene(source: &Molecule, reverse: bool) -> Molecule {
    let mut editor = MoleculeEditor::new();
    let mut atom_ids = source.atom_ids().collect::<Vec<_>>();
    let mut bonds = source.bonds().collect::<Vec<_>>();
    if reverse {
        atom_ids.reverse();
        bonds.reverse();
    }
    let atoms = atom_ids
        .into_iter()
        .map(|id| {
            let mut atom = source.atom(id).unwrap().clone();
            atom.atom_map = Some(id.raw() + 1);
            (id, editor.add_atom(atom).unwrap())
        })
        .collect::<BTreeMap<_, _>>();
    let bonds = bonds
        .into_iter()
        .map(|(id, bond)| {
            let (a, b) = if reverse {
                (bond.b(), bond.a())
            } else {
                bond.endpoints()
            };
            (
                id,
                editor.add_bond(atoms[&a], atoms[&b], bond.order).unwrap(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    let carrier = |carrier| match carrier {
        StereoCarrier::Atom(id) => StereoCarrier::Atom(atoms[&id]),
        other => other,
    };
    for (_, element) in source.stereo_elements() {
        let StereoElementKind::DoubleBond(mut stereo) = element.kind.clone() else {
            panic!("expected double-bond assertion");
        };
        stereo.bond = bonds[&stereo.bond];
        stereo.left = atoms[&stereo.left];
        stereo.right = atoms[&stereo.right];
        stereo.left_carrier = carrier(stereo.left_carrier);
        stereo.right_carrier = carrier(stereo.right_carrier);
        editor
            .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(stereo)))
            .unwrap();
    }
    editor.finish().unwrap()
}

fn mapped_cip(molecule: &mut Molecule) -> BTreeMap<(u32, u32), StereoDescriptor> {
    molecule.perceive().unwrap();
    stereo::assign_cip_descriptors(molecule).unwrap();
    molecule
        .stereo_elements()
        .map(|(id, element)| {
            let StereoElementKind::DoubleBond(stereo) = &element.kind else {
                panic!("expected double-bond assertion");
            };
            let left = molecule.atom(stereo.left).unwrap().atom_map.unwrap();
            let right = molecule.atom(stereo.right).unwrap().atom_map.unwrap();
            (
                (left.min(right), left.max(right)),
                molecule.cip_descriptor(id).unwrap().unwrap(),
            )
        })
        .collect()
}

#[test]
fn conjugated_stereo_export_is_independent_of_atom_bond_and_endpoint_numbering() {
    let source = molecule("F/C=C\\C=C/C=C\\F");
    let mut baseline = mapped_polyene(&source, false);
    let expected = mapped_cip(&mut baseline);
    assert_eq!(expected.len(), 3);
    for reverse in [false, true] {
        let mut original = mapped_polyene(&source, reverse);
        assert_eq!(mapped_cip(&mut original), expected);
        let written = smiles::write_isomeric(&original).unwrap();
        let mut restored = molecule(&written);
        assert_eq!(mapped_cip(&mut restored), expected, "{written}");
    }
}

#[test]
fn conjugated_stereo_round_trips_with_implicit_hydrogen_carriers() {
    let mut original = mapped_polyene(&molecule("F/C=C\\C=C/C=C\\F"), false);
    let expected = mapped_cip(&mut original);
    for implicit_carriers in 0..64 {
        let mut editor = original.clone().into_editor();
        for (index, (id, element)) in original.stereo_elements().enumerate() {
            let mut replacement = element.clone();
            let StereoElementKind::DoubleBond(stereo) = &mut replacement.kind else {
                panic!("expected double-bond assertion");
            };
            let left_implicit = implicit_carriers & (1 << (index * 2)) != 0;
            let right_implicit = implicit_carriers & (1 << (index * 2 + 1)) != 0;
            if left_implicit {
                stereo.left_carrier = StereoCarrier::ImplicitHydrogen;
            }
            if right_implicit {
                stereo.right_carrier = StereoCarrier::ImplicitHydrogen;
            }
            // Replacing exactly one substituent by the other one reverses the
            // relative assertion while retaining the same E/Z configuration.
            if left_implicit != right_implicit {
                stereo.orientation = Some(match stereo.orientation.unwrap() {
                    DoubleBondOrientation::Together => DoubleBondOrientation::Opposite,
                    DoubleBondOrientation::Opposite => DoubleBondOrientation::Together,
                });
            }
            editor.replace_stereo_element(id, replacement).unwrap();
        }
        let mut source = editor.finish().unwrap();
        assert_eq!(mapped_cip(&mut source), expected);
        let written = smiles::write_isomeric(&source).unwrap();
        let mut restored = molecule(&written);
        assert_eq!(mapped_cip(&mut restored), expected, "{written}");
    }
}

#[test]
fn contradictory_double_bond_assertions_are_rejected_before_export() {
    let mut original = molecule("F/C=C/F");
    original.perceive().unwrap();
    stereo::assign_cip_descriptors(&mut original).unwrap();
    let (id, element) = original.stereo_elements().next().unwrap();
    let expected = original.cip_descriptor(id).unwrap();
    let expected_export = smiles::write_isomeric(&original).unwrap();
    let mut opposite = element.clone();
    let StereoElementKind::DoubleBond(stereo) = &mut opposite.kind else {
        panic!("expected double-bond assertion");
    };
    stereo.orientation = Some(match stereo.orientation.unwrap() {
        DoubleBondOrientation::Together => DoubleBondOrientation::Opposite,
        DoubleBondOrientation::Opposite => DoubleBondOrientation::Together,
    });
    let mut editor = original.into_editor();
    let before = editor.clone();
    assert!(matches!(
        editor.add_stereo_element(opposite),
        Err(MoleculeError::InvalidStereoReference(_))
    ));
    assert_eq!(editor, before);

    let preserved = editor.finish().unwrap();
    assert_eq!(preserved.stereo_elements().count(), 1);
    assert_eq!(preserved.cip_descriptor(id).unwrap(), expected);
    let written = smiles::write_isomeric(&preserved).unwrap();
    assert_eq!(written, expected_export);
    let mut restored = molecule(&written);
    restored.perceive().unwrap();
    stereo::assign_cip_descriptors(&mut restored).unwrap();
    let restored_id = restored.stereo_element_ids().next().unwrap();
    assert_eq!(restored.cip_descriptor(restored_id).unwrap(), expected);
}
