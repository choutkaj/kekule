use kekule::{core::*, smiles, topology::Topology};

fn molecule(text: &str) -> Molecule {
    let mut mol = smiles::parse_str(text)
        .unwrap()
        .interpret()
        .unwrap()
        .into_molecule()
        .unwrap();
    mol.perceive().unwrap();
    mol
}

fn canonical(text: &str) -> String {
    smiles::write_canonical(&molecule(text)).unwrap()
}

#[test]
fn enhanced_groups_round_trip_and_canonicalize_joint_inversion() {
    for tag in ["&17", "o17"] {
        let a = format!("F[C@H](Cl)[C@@H](Br)I |{tag}:1,3|");
        let inverted = format!("F[C@@H](Cl)[C@H](Br)I |{tag}:3,1|");
        let diastereomer = format!("F[C@H](Cl)[C@H](Br)I |{tag}:1,3|");
        let mol = molecule(&a);
        let output = smiles::write_isomeric(&mol).unwrap();
        let restored = molecule(&output);
        assert_eq!(
            restored.stereo_groups().next().unwrap().1.kind,
            mol.stereo_groups().next().unwrap().1.kind
        );
        assert_eq!(canonical(&a), canonical(&inverted));
        assert_ne!(canonical(&a), canonical(&diastereomer));
        let text = canonical(&a);
        assert_eq!(canonical(&text), text);
    }
    assert_ne!(
        canonical("F[C@H](Cl)Br |o1:1|"),
        canonical("F[C@H](Cl)Br |&1:1|")
    );
    assert_ne!(
        canonical("F[C@H](Cl)Br |a:1|"),
        canonical("F[C@@H](Cl)Br |a:1|")
    );
}

#[test]
fn topology_export_has_one_extension_and_record_global_indices() {
    let first = molecule("F[C@H](Cl)Br |o1:1|");
    let second = molecule("N[C@@H](C)O |&1:1|");
    for mode in [
        smiles::SmilesWriteMode::Isomeric,
        smiles::SmilesWriteMode::Canonical,
    ] {
        let topology = Topology::from_molecules(&[first.clone(), second.clone()]).unwrap();
        let text = smiles::write_topology(&topology, smiles::SmilesWriteOptions { mode }).unwrap();
        assert_eq!(text.matches('|').count(), 2, "{text}");
        let restored = smiles::to_molecules(&text).unwrap();
        assert_eq!(restored.len(), 2);
        let mut kinds = restored
            .iter()
            .map(|m| m.stereo_groups().next().unwrap().1.kind as u8)
            .collect::<Vec<_>>();
        kinds.sort();
        let mut expected = vec![StereoGroupKind::And as u8, StereoGroupKind::Or as u8];
        expected.sort();
        assert_eq!(kinds, expected);
    }
}

#[test]
fn relative_flag_does_not_weaken_ungrouped_absolute_centers() {
    let mut edit = molecule("F[C@H](Cl)[C@@H](Br)I |o1:1|").edit();
    let (id, mut group) = edit
        .stereo_groups()
        .next()
        .map(|(id, g)| (id, g.clone()))
        .unwrap();
    group.kind = StereoGroupKind::Relative;
    edit.replace_stereo_group(id, group).unwrap();
    let mol = edit.finish().unwrap();
    let text = smiles::write_isomeric(&mol).unwrap();
    let restored = molecule(&text);
    assert_eq!(
        restored
            .stereo_groups()
            .filter(|(_, g)| g.kind == StereoGroupKind::Relative)
            .count(),
        1
    );
    assert_eq!(
        restored
            .stereo_groups()
            .filter(|(_, g)| g.kind == StereoGroupKind::Absolute)
            .count(),
        1
    );
    let text = smiles::write_canonical(&mol).unwrap();
    assert_eq!(canonical(&text), text);
}

fn matches(target: &str, query: &str) -> bool {
    use kekule::substructure::{find_substructure_matches_with_options, SubstructureMatchOptions};
    let query = kekule::query::parse_smarts(query).unwrap();
    !find_substructure_matches_with_options(
        &molecule(target),
        &query,
        SubstructureMatchOptions {
            use_enhanced_stereo: true,
            ..Default::default()
        },
    )
    .unwrap()
    .is_empty()
}

#[test]
fn matching_preserves_group_correlation_and_sample_semantics() {
    let base = "F[C@H](Cl)[C@@H](Br)I";
    let inverse = "F[C@@H](Cl)[C@H](Br)I";
    let other = "F[C@H](Cl)[C@H](Br)I";
    for tag in ["o1", "&1"] {
        let target = format!("{base} |{tag}:1,3|");
        assert!(matches(&target, base));
        assert!(matches(&target, inverse));
        assert!(!matches(&target, other));
        assert!(matches(&target, &format!("{inverse} |{tag}:1,3|")));
        assert!(!matches(&target, &format!("{other} |{tag}:1,3|")));
        let independent = format!("{base} |{tag}:1,{tag}2:3|")
            .replace("o12", "o2")
            .replace("&12", "&2");
        assert!(!matches(&target, &independent));
        assert!(matches(&independent, &format!("{other} |{tag}:1,3|")));
        assert!(!matches(base, &format!("{base} |{tag}:1,3|")));
    }
    assert!(matches(
        &format!("{base} |&1:1,3|"),
        &format!("{inverse} |o1:1,3|")
    ));
    assert!(!matches(
        &format!("{base} |o1:1,3|"),
        &format!("{base} |&1:1,3|")
    ));
    assert!(!matches(base, inverse));
}

#[test]
fn cxsmarts_groups_validate_indices_and_syntax_without_panicking() {
    assert_eq!(
        kekule::query::parse_smarts("C ||").unwrap(),
        kekule::query::parse_smarts("C").unwrap()
    );
    for query in [
        "C |",
        "C |||",
        "C |o1:0|",
        "F[C@H](Cl)Br |o1:1,o2:1|",
        "F[C@H](Cl)Br |o+1:1|",
        "F[C@H](Cl)Br |o1:99|",
        "F[C@H](Cl)Br |o1:1,|",
        "F[C@H](Cl)Br |r,r|",
    ] {
        assert!(kekule::query::parse_smarts(query).is_err(), "{query}");
    }
    let query = kekule::query::parse_smarts("F[C@H](Cl)[C@@H](Br)I |a:1,r|").unwrap();
    assert_eq!(query.stereo_groups().len(), 2);
    assert_eq!(query.to_builder().build().unwrap(), query);
}

#[test]
fn topology_edits_reject_separating_correlated_groups() {
    use std::sync::Arc;
    for tag in ["a", "o1", "&1"] {
        let original = molecule(&format!("F[C@H](Cl)CCC[C@H](Br)I |{tag}:1,6|"));
        let source = Arc::new(Topology::from_molecule(&original).unwrap());
        let mut edit = source.edit();
        let instance = source.molecules().next().unwrap().id();
        let local = original
            .bond_between(AtomId::new(3), AtomId::new(4))
            .unwrap()
            .unwrap();
        let cut = edit
            .bond_handle(kekule::topology::InstanceBondId::new(instance, local))
            .unwrap();
        edit.delete_bond(cut).unwrap();
        let result = edit.finish();
        if tag == "a" {
            let restored = result.unwrap();
            assert_eq!(restored.instance_count(), 2);
            assert_eq!(
                restored
                    .molecules()
                    .map(|m| m.molecule().stereo_groups().count())
                    .sum::<usize>(),
                2
            );
        } else {
            assert!(result
                .unwrap_err()
                .to_string()
                .contains("correlated stereo group"));
        }
        assert_eq!(source.molecules().next().unwrap().molecule(), &original);
        let selected = kekule::topology::AtomSelection::from_atoms(
            &source,
            source
                .atom_ids()
                .iter()
                .copied()
                .filter(|a| a.atom().raw() != 4),
        )
        .unwrap();
        // Removing the bridge atom leaves both stereo reference frames intact.
        assert_eq!(source.subset(&selected).is_ok(), tag == "a");
    }
}

#[test]
fn group_aware_symmetry_cleanup_compares_whole_relationships() {
    for (fields, count) in [
        ("o1:1,4", 2),
        ("&1:1,4", 2),
        ("o1:1,o2:4", 2),
        ("o1:1,&1:4", 3),
        ("a:1,o1:4", 3),
    ] {
        let source = molecule(&format!("[C@]([C@H](F)Cl)([C@@H](F)Cl)(Br)I |{fields}|"));
        let mut edit = source.edit();
        let report = kekule::stereo::cleanup_stereo(&mut edit, Default::default()).unwrap();
        assert_eq!(report.removed_elements.len(), 3 - count, "{fields}");
        let cleaned = edit.finish().unwrap();
        assert_eq!(cleaned.stereo_elements().count(), count);
        assert_eq!(
            cleaned.stereo_groups().count(),
            source.stereo_groups().count()
        );
        let text = smiles::write_canonical(&cleaned).unwrap();
        assert_eq!(canonical(&text), text);
    }
}

#[test]
fn hydrogen_transforms_preserve_membership_and_configuration() {
    for tag in ["a", "o1", "&1"] {
        let mut mol = molecule(&format!("F[C@H](Cl)[C@@H](Br)I |{tag}:1,3|"));
        let expected = smiles::write_canonical(&mol).unwrap();
        let groups = mol
            .stereo_groups()
            .map(|(id, g)| (id, g.clone()))
            .collect::<Vec<_>>();
        mol.add_hydrogens().unwrap();
        mol.perceive().unwrap();
        assert_eq!(smiles::write_canonical(&mol).unwrap(), expected);
        assert_eq!(
            mol.stereo_groups()
                .map(|(id, g)| (id, g.clone()))
                .collect::<Vec<_>>(),
            groups
        );
        mol.remove_hydrogens().unwrap();
        mol.perceive().unwrap();
        assert_eq!(smiles::write_canonical(&mol).unwrap(), expected);
    }
}

#[test]
fn topology_matching_does_not_correlate_reused_definitions() {
    use kekule::{query::parse_smarts, substructure::*, topology::TopologyBuilder};
    use std::sync::Arc;
    let mol = molecule("F[C@H](Cl)Br |o1:1|");
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&mol).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    // Each occurrence can choose its own orientation even though the definition
    // and its local group ID are shared.
    let query = parse_smarts("F[C@H](Cl)Br.F[C@@H](Cl)Br").unwrap();
    let matches = find_topology_substructure_matches_complete(
        &topology,
        &query,
        SubstructureMatchOptions {
            use_enhanced_stereo: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!matches.is_empty());
}

fn renumber(source: &Molecule, shift: usize) -> Molecule {
    use std::collections::BTreeMap;
    let mut edit = MoleculeEditor::new();
    let mut order = source.atom_ids().collect::<Vec<_>>();
    let len = order.len();
    order.rotate_left(shift % len);
    if shift % 2 == 1 {
        order.reverse();
    }
    let atoms = order
        .into_iter()
        .map(|a| (a, edit.add_atom(source.atom(a).unwrap().clone()).unwrap()))
        .collect::<BTreeMap<_, _>>();
    for (_, bond) in source.bonds().collect::<Vec<_>>().into_iter().rev() {
        edit.add_bond(atoms[&bond.b()], atoms[&bond.a()], bond.order)
            .unwrap();
    }
    let mut elements = BTreeMap::new();
    for (id, element) in source
        .stereo_elements()
        .collect::<Vec<_>>()
        .into_iter()
        .rev()
    {
        let StereoElementKind::Tetrahedral(mut tetra) = element.kind.clone() else {
            panic!("tetrahedral test graph");
        };
        tetra.center = atoms[&tetra.center];
        for carrier in &mut tetra.carriers {
            if let StereoCarrier::Atom(a) = carrier {
                *a = atoms[a];
            }
        }
        elements.insert(
            id,
            edit.add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(tetra)))
                .unwrap(),
        );
    }
    for (_, group) in source.stereo_groups().collect::<Vec<_>>().into_iter().rev() {
        edit.add_stereo_group(StereoGroup {
            kind: group.kind,
            members: group.members.iter().rev().map(|s| elements[s]).collect(),
        })
        .unwrap();
    }
    let mut mol = edit.finish().unwrap();
    mol.perceive().unwrap();
    mol
}

#[test]
fn canonical_groups_are_invariant_to_atom_bond_element_and_group_numbering() {
    for text in [
        "F[C@H](Cl)[C@@H](Br)I |o1:1,3|",
        "C1[C@H](F)CC[C@@H](Cl)C1 |&9:1,5|",
        "C[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C |o4:1,5,o2:3,7|",
        "C[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C |a:1,&9:3,7,o5:5|",
    ] {
        let mol = molecule(text);
        let expected = smiles::write_canonical(&mol).unwrap();
        for shift in 0..mol.atom_count() {
            assert_eq!(
                smiles::write_canonical(&renumber(&mol, shift)).unwrap(),
                expected,
                "{text}, shift={shift}"
            );
        }
    }
}

#[test]
fn unencodable_group_semantics_fail_explicitly() {
    let source = molecule("F[C@H](Cl)[C@@H](Br)I |o1:1,o2:3|");
    for kind in [StereoGroupKind::Relative, StereoGroupKind::Racemic] {
        let mut edit = source.edit();
        for (id, mut group) in source.stereo_groups().map(|(id, g)| (id, g.clone())) {
            group.kind = kind;
            edit.replace_stereo_group(id, group).unwrap();
        }
        let mol = edit.finish().unwrap();
        assert!(smiles::write_isomeric(&mol).is_err());
        assert!(smiles::write_canonical(&mol).is_err());
    }
    assert!(smiles::to_molecules("F[C@H](Cl)Br.F[C@H](Cl)Br |o1:1,5|").is_err());
}

#[test]
fn matching_agrees_with_independent_rdkit_enhanced_stereo_matrix() {
    // Independently evaluated with RDKit 2026.03.3, useChirality=True and
    // useEnhancedStereo=True. Toy molecules are focused semantic regressions.
    let cases = [
        "F[C@H](Cl)[C@@H](Br)I",
        "F[C@@H](Cl)[C@H](Br)I",
        "F[C@H](Cl)[C@H](Br)I",
        "F[C@H](Cl)[C@@H](Br)I |o1:1,3|",
        "F[C@@H](Cl)[C@H](Br)I |o1:1,3|",
        "F[C@H](Cl)[C@@H](Br)I |&1:1,3|",
        "F[C@@H](Cl)[C@H](Br)I |&1:1,3|",
        "F[C@H](Cl)[C@@H](Br)I |o1:1,o2:3|",
        "F[C@H](Cl)[C@H](Br)I |o1:1,3|",
    ];
    let expected = [
        "100000000",
        "010000000",
        "001000000",
        "110110000",
        "110110000",
        "110111100",
        "110111100",
        "111110011",
        "001000001",
    ];
    for (target, row) in cases.iter().zip(expected) {
        for (query, expected) in cases.iter().zip(row.bytes()) {
            assert_eq!(
                matches(target, query),
                expected == b'1',
                "target={target}, query={query}"
            );
        }
    }
}

#[test]
fn enhanced_matching_respects_boolean_and_partial_carrier_queries() {
    let target = "N[C@](F)(Cl)Br |o1:1|";
    for query in [
        "N[C@@](F)(Cl)Br",
        "N[C!@](F)(Cl)Br",
        "N[C@,C@@](F)(Cl)Br |o1:1|",
        "N[C@] |o1:1|",
    ] {
        assert!(matches(target, query), "{query}");
    }
    assert!(!matches("NC(F)(Cl)Br", "N[C!@](F)(Cl)Br |o1:1|"));
    assert!(!matches(target, "N[C@,C@@](F)(Cl)Br |&1:1|"));
    let query = kekule::query::parse_smarts("N[C@@](F)(Cl)Br").unwrap();
    assert!(
        kekule::substructure::find_substructure_match(&molecule(target), &query)
            .unwrap()
            .is_none()
    );
}

#[test]
fn canonical_projection_merges_absolute_sets_and_shields_relative_groups() {
    let source = molecule("C[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C |o1:1,o2:3,o3:5,o4:7|");
    for relative in [false, true] {
        let mut edit = source.edit();
        for (i, (id, group)) in source.stereo_groups().enumerate() {
            let mut group = group.clone();
            group.kind = if relative && i == 0 {
                StereoGroupKind::Relative
            } else {
                StereoGroupKind::Absolute
            };
            edit.replace_stereo_group(id, group).unwrap();
        }
        let mol = edit.finish().unwrap();
        let expected = smiles::write_canonical(&mol).unwrap();
        assert_eq!(canonical(&expected), expected);
        for shift in 0..mol.atom_count() {
            assert_eq!(
                smiles::write_canonical(&renumber(&mol, shift)).unwrap(),
                expected
            );
        }
    }
    assert_eq!(
        canonical("F[C@H](Cl)[C@@H](Br)I |r|"),
        canonical("F[C@@H](Cl)[C@H](Br)I |r|")
    );
    assert_ne!(
        canonical("F[C@H](Cl)[C@@H](Br)I |r|"),
        canonical("F[C@H](Cl)[C@H](Br)I |r|")
    );
}
