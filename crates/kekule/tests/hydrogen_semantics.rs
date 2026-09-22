use kekule::core::{AtomId, BondOrder, HydrogenDeclaration, Molecule, Perception};
use kekule::descriptors::{molecular_formula, HydrogenCountPolicy};
use kekule::hydrogens::{AddHydrogensOptions, AddedHydrogenOrigin, HydrogenTransformError};
use kekule::topology::{InstanceAtomId, TopologyBuilder};
use kekule::{smiles, stereo};

fn parse(source: &str) -> Molecule {
    smiles::to_molecules(source).unwrap().pop().unwrap()
}

fn carbon(molecule: &Molecule) -> AtomId {
    molecule
        .atoms()
        .find(|(_, a)| a.element.symbol() == "C")
        .unwrap()
        .0
}

#[test]
fn counts_describe_representation_independently_of_how_implicit_h_were_specified() {
    for (source, explicit, implicit) in [
        ("C", 0, 4),
        ("[CH4]", 0, 4),
        ("[H]C([H])([H])[H]", 4, 0),
        ("[H][CH3]", 1, 3),
        ("[2H][CH3]", 1, 3),
    ] {
        let mut molecule = parse(source);
        molecule.perceive().unwrap();
        let c = carbon(&molecule);
        assert_eq!(
            molecule.explicit_hydrogens(c).unwrap(),
            explicit,
            "{source}"
        );
        assert_eq!(
            molecule.implicit_hydrogens(c).unwrap(),
            Some(implicit),
            "{source}"
        );
        assert_eq!(molecule.total_hydrogens(c).unwrap(), Some(4), "{source}");
        let formula = molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived).unwrap();
        assert_eq!(
            formula.count(kekule::core::Element::from_symbol("H").unwrap()),
            4
        );
    }
}

#[test]
fn specified_counts_are_known_without_perception_and_unknown_is_not_zero() {
    for (source, expected) in [
        ("C", None),
        ("[C]", Some(0)),
        ("[CH3]", Some(3)),
        ("[NH4+]", Some(4)),
    ] {
        let molecule = parse(source);
        let atom = molecule.atom_ids().next().unwrap();
        assert!(!molecule.perception().has_valence());
        assert_eq!(
            molecule.implicit_hydrogens(atom).unwrap(),
            expected,
            "{source}"
        );
        assert_eq!(
            molecule.total_hydrogens(atom).unwrap(),
            expected,
            "{source}"
        );
        assert_eq!(molecule.explicit_hydrogens(atom).unwrap(), 0);
    }
    let molecule = parse("C");
    let invalid = AtomId::new(999);
    assert!(molecule.explicit_hydrogens(invalid).is_err());
    assert!(molecule.implicit_hydrogens(invalid).is_err());
    assert!(molecule.total_hydrogens(invalid).is_err());
}

#[test]
fn specified_and_inferred_contributions_are_combined_without_double_counting() {
    let mut editor = parse("C").into_editor();
    let c = editor.atom_ids().next().unwrap();
    editor.atom_mut(c).unwrap().hydrogens = HydrogenDeclaration::Infer { specified: 1 };
    let mut molecule = editor.finish().unwrap();
    assert_eq!(molecule.implicit_hydrogens(c).unwrap(), None);
    molecule.perceive().unwrap();
    assert_eq!(molecule.inferred_hydrogens(c).unwrap(), Some(3));
    assert_eq!(molecule.implicit_hydrogens(c).unwrap(), Some(4));
    let report = molecule
        .add_hydrogens_with_options(AddHydrogensOptions {
            specified_only: true,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(report.added.len(), 1);
    assert_eq!(report.added[0].origin, AddedHydrogenOrigin::Specified);
    molecule.perceive().unwrap();
    assert_eq!(molecule.explicit_hydrogens(c).unwrap(), 1);
    assert_eq!(molecule.implicit_hydrogens(c).unwrap(), Some(3));
    assert_eq!(molecule.total_hydrogens(c).unwrap(), Some(4));
}

#[test]
fn edits_invalidate_inference_and_preserve_specified_counts() {
    for (source, after_edit, after_perception) in [("CC", None, 2), ("[CH3]C", Some(3), 3)] {
        let mut molecule = parse(source);
        molecule.perceive().unwrap();
        let c = carbon(&molecule);
        let bond = molecule.bond_ids().next().unwrap();
        let mut editor = molecule.into_editor();
        editor.bond_mut(bond).unwrap().set_order(BondOrder::Double);
        assert_eq!(editor.implicit_hydrogens(c).unwrap(), after_edit);
        assert_eq!(editor.explicit_hydrogens(c).unwrap(), 0);
        assert_eq!(editor.total_hydrogens(c).unwrap(), after_edit);
        let mut molecule = editor.finish().unwrap();
        if source.starts_with('[') {
            // The edit makes the fixed count overvalent. Strict validation must
            // report it rather than silently removing a specified hydrogen.
            assert!(molecule.perceive().is_err());
            assert_eq!(
                molecule.atom(c).unwrap().hydrogens,
                HydrogenDeclaration::Fixed(3)
            );
        } else {
            molecule.perceive().unwrap();
        }
        assert_eq!(
            molecule.implicit_hydrogens(c).unwrap(),
            Some(after_perception)
        );
    }
}

#[test]
fn topology_counts_follow_each_reused_definition_and_reject_foreign_atoms() {
    let mut molecule = parse("[H][CH3]");
    molecule.perceive().unwrap();
    let c = carbon(&molecule);
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    let first = builder.add_instance(definition).unwrap();
    let second = builder.add_instance(definition).unwrap();
    let topology = builder.build().unwrap();
    for instance in [first, second] {
        let atom = InstanceAtomId::new(instance, c);
        assert_eq!(topology.explicit_hydrogens(atom).unwrap(), 1);
        assert_eq!(topology.implicit_hydrogens(atom).unwrap(), Some(3));
        assert_eq!(topology.total_hydrogens(atom).unwrap(), Some(4));
    }
    let invalid = InstanceAtomId::new(first, AtomId::new(999));
    assert!(topology.explicit_hydrogens(invalid).is_err());
    assert!(topology.implicit_hydrogens(invalid).is_err());
    assert!(topology.total_hydrogens(invalid).is_err());
}

#[test]
fn conversion_preserves_composition_and_reports_explicit_and_implicit_counts() {
    for source in ["C", "[CH4]", "[H][CH3]", "[2H][CH3]", "[H:7][CH3]"] {
        let mut molecule = parse(source);
        molecule.perceive().unwrap();
        let c = carbon(&molecule);
        let formula = molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived).unwrap();
        let original_implicit = molecule.implicit_hydrogens(c).unwrap().unwrap();
        assert_eq!(
            molecule.add_hydrogens().unwrap().added.len(),
            original_implicit
        );
        molecule.perceive().unwrap();
        assert_eq!(molecule.explicit_hydrogens(c).unwrap(), 4);
        assert_eq!(molecule.implicit_hydrogens(c).unwrap(), Some(0));
        let report = molecule.remove_hydrogens().unwrap();
        let retained = usize::from(source.contains("2H") || source.contains(":7"));
        assert_eq!(report.retained.len(), retained);
        assert_eq!(report.adjustments[0].explicit_hydrogens, retained);
        assert_eq!(report.adjustments[0].implicit_hydrogens, 4 - retained);
        molecule.perceive().unwrap();
        assert_eq!(molecule.total_hydrogens(c).unwrap(), Some(4));
        assert_eq!(
            molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived).unwrap(),
            formula
        );
    }
}

#[test]
fn fixed_counts_materialize_without_perception_but_unresolved_counts_fail_atomically() {
    let mut fixed = parse("[CH4]");
    assert_eq!(fixed.add_hydrogens().unwrap().added.len(), 4);
    let mut inferred = parse("C");
    let before = smiles::write_isomeric(&inferred).unwrap();
    assert_eq!(
        inferred.add_hydrogens(),
        Err(HydrogenTransformError::MissingValencePerception)
    );
    assert_eq!(smiles::write_isomeric(&inferred).unwrap(), before);
}

#[test]
fn stereo_and_smiles_round_trips_use_complete_implicit_counts() {
    for source in ["[C]", "[CH3]", "[NH4+]", "c1cc[nH]c1", "F[C@H](Cl)Br"] {
        let mut molecule = parse(source);
        molecule.perceive().unwrap();
        let formula = molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived).unwrap();
        let cip = stereo::assign_cip_descriptors(&mut molecule)
            .unwrap()
            .assigned;
        molecule.add_hydrogens().unwrap();
        molecule.perceive().unwrap();
        molecule.remove_hydrogens().unwrap();
        molecule.perceive().unwrap();
        assert_eq!(
            stereo::assign_cip_descriptors(&mut molecule)
                .unwrap()
                .assigned,
            cip
        );
        for write in [smiles::write_isomeric, smiles::write_canonical] {
            let text = write(&molecule).unwrap();
            let mut restored = parse(&text);
            restored.perceive().unwrap();
            assert_eq!(
                molecular_formula(&restored, HydrogenCountPolicy::IncludePerceived).unwrap(),
                formula,
                "{source} -> {text}"
            );
        }
    }
}

#[test]
fn count_arithmetic_does_not_truncate_large_detached_assignments() {
    let mut editor = parse("C").into_editor();
    let c = editor.atom_ids().next().unwrap();
    editor.atom_mut(c).unwrap().hydrogens = HydrogenDeclaration::Infer { specified: 255 };
    let mut molecule = editor.finish().unwrap();
    let perception = Perception::builder()
        .with_valence(None, vec![(c, 255)])
        .unwrap()
        .build();
    molecule.install_perception(perception).unwrap();
    assert_eq!(molecule.implicit_hydrogens(c).unwrap(), Some(510));
    assert_eq!(molecule.total_hydrogens(c).unwrap(), Some(510));
}

#[test]
fn incomplete_perception_cannot_silently_drop_hydrogens_during_materialization() {
    let mut molecule = parse("CC");
    let atoms = molecule.atom_ids().collect::<Vec<_>>();
    let state = Perception::builder()
        .with_valence(None, vec![(atoms[0], 3)])
        .unwrap()
        .build();
    molecule.install_perception(state).unwrap();
    assert_eq!(molecule.implicit_hydrogens(atoms[0]).unwrap(), Some(3));
    assert_eq!(molecule.implicit_hydrogens(atoms[1]).unwrap(), None);
    let before = molecule.perception().clone();
    assert_eq!(
        molecule.add_hydrogens(),
        Err(HydrogenTransformError::MissingValencePerception)
    );
    assert_eq!(molecule.atom_count(), 2);
    assert_eq!(molecule.perception(), &before);
}
