use kekule::core::Molecule;
use kekule::perception::aromaticity::{
    perceive_aromaticity_with_options, AromaticityError, AromaticityModel, AromaticityOptions,
};
use kekule::perception::rings::RingPerceptionOptions;
use kekule::smiles;
use kekule::stereo::assign_cip_descriptors;

fn molecule(source: &str) -> Molecule {
    smiles::to_molecules(source).unwrap().pop().unwrap()
}

#[test]
fn aromaticity_work_failures_restore_both_missing_and_installed_perception() {
    let mut perceived = molecule("F[C@H](Cl)c1ccccc1-c1ccccc1");
    perceived.perceive().unwrap();
    assign_cip_descriptors(&mut perceived).unwrap();
    for original in [molecule("F[C@H](Cl)c1ccccc1-c1ccccc1"), perceived] {
        let mut failures = 0;
        let mut successes = 0;
        // Exercise failures at successive stages, including after rings have
        // been installed and after aromaticity has replaced the old section.
        for max_total_work in (0..1_000).step_by(7) {
            let mut actual = original.clone();
            let options = AromaticityOptions {
                max_total_work,
                ..Default::default()
            };
            match perceive_aromaticity_with_options(
                &mut actual,
                AromaticityModel::RdkitLike,
                options,
            ) {
                Err(error) => {
                    assert_eq!(
                        error,
                        AromaticityError::ResourceLimit {
                            limit: max_total_work
                        }
                    );
                    assert_eq!(actual.perception(), original.perception());
                    failures += 1;
                }
                Ok(()) => {
                    assert_eq!(
                        actual
                            .atoms()
                            .filter(|(id, _)| actual.atom_is_aromatic(*id).unwrap() == Some(true))
                            .count(),
                        12
                    );
                    assert_eq!(
                        actual
                            .bonds()
                            .filter(|(id, _)| actual.bond_is_aromatic(*id).unwrap() == Some(true))
                            .count(),
                        12
                    );
                    successes += 1;
                }
            }
            assert_eq!(actual, original, "represented chemistry changed");
        }
        assert!(failures > 0 && successes > 0);
    }
}

#[test]
fn installed_rings_skip_enumeration_limits_but_not_aromaticity_limits() {
    let mut actual = molecule("c1ccccc1");
    actual.perceive().unwrap();
    let options = AromaticityOptions {
        ring_options: RingPerceptionOptions {
            max_atoms: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    perceive_aromaticity_with_options(&mut actual, AromaticityModel::RdkitLike, options).unwrap();
    let before = actual.perception().clone();
    assert!(matches!(
        perceive_aromaticity_with_options(
            &mut actual,
            AromaticityModel::RdkitLike,
            AromaticityOptions {
                max_total_work: 0,
                ..options
            }
        ),
        Err(AromaticityError::ResourceLimit { limit: 0 })
    ));
    assert_eq!(actual.perception(), &before);
}

#[test]
fn many_smiles_components_preserve_source_mappings_and_stereo() {
    // Branches interleave disconnected atoms and ring closures reconnect
    // source fragments. Components cannot be obtained by merely splitting '.'.
    let block = "C1(O.F)CC1.N[C@@H](C)C(=O)O.F/C=C/Cl.c1ccncc1";
    let expected = smiles::parse_str(block).unwrap().interpret().unwrap();
    let source = vec![block; 256].join(".");
    let actual = smiles::parse_str(&source).unwrap().interpret().unwrap();
    assert_eq!(actual.components().len(), 256 * expected.components().len());
    for (index, component) in actual.components().iter().enumerate() {
        let baseline = &expected.components()[index % expected.components().len()];
        let offset = index / expected.components().len() * (block.len() + 1);
        assert_eq!(component.molecule(), baseline.molecule());
        assert_eq!(
            component.source_span(),
            baseline.source_span().start + offset..baseline.source_span().end + offset
        );
        let report = component.report();
        let expected_report = baseline.report();
        assert_eq!(
            report.atom_mappings().len(),
            expected_report.atom_mappings().len()
        );
        assert_eq!(
            report.bond_mappings().len(),
            expected_report.bond_mappings().len()
        );
        for (mapping, expected) in report
            .atom_mappings()
            .iter()
            .zip(expected_report.atom_mappings())
        {
            assert_eq!(mapping.atom(), expected.atom());
            assert_eq!(
                mapping.source_span(),
                expected.source_span().start + offset..expected.source_span().end + offset
            );
        }
        for (mapping, expected) in report
            .bond_mappings()
            .iter()
            .zip(expected_report.bond_mappings())
        {
            assert_eq!(mapping.bond(), expected.bond());
            assert_eq!(mapping.source_offset(), expected.source_offset() + offset);
        }
        assert_eq!(
            report.created_stereo_elements(),
            expected_report.created_stereo_elements()
        );
    }
}

#[test]
fn parse_perceive_cip_and_write_are_stable_across_supported_chemistry() {
    for source in [
        "C1CC2CCC1C2",
        "C12C3C4C1C5C2C3C45",
        "c1ccc2ccccc2c1",
        "c1cc[nH]c1",
        "[nH+]1ccccc1",
        "[O-][N+](=O)c1ccccc1",
        "N[C@@H](C)C(=O)O",
        "F/C=C/C=C\\Cl",
        "[2H][C@](F)(Cl)Br",
    ] {
        let mut original = molecule(source);
        original.perceive().unwrap();
        let expected_cip = assign_cip_descriptors(&mut original).unwrap();
        let first = smiles::write_canonical(&original).unwrap();
        let mut restored = molecule(&first);
        restored.perceive().unwrap();
        let restored_cip = assign_cip_descriptors(&mut restored).unwrap();
        let descriptors = |report: &kekule::stereo::CipAssignmentReport| {
            let mut values = report
                .assigned
                .iter()
                .map(|assignment| format!("{:?}", assignment.descriptor))
                .collect::<Vec<_>>();
            values.sort();
            values
        };
        assert_eq!(
            descriptors(&restored_cip),
            descriptors(&expected_cip),
            "{source}"
        );
        assert_eq!(
            restored_cip.skipped.len(),
            expected_cip.skipped.len(),
            "{source}"
        );
        assert_eq!(
            smiles::write_canonical(&restored).unwrap(),
            first,
            "{source}"
        );
        let before = restored.perception().clone();
        restored.perceive().unwrap();
        assign_cip_descriptors(&mut restored).unwrap();
        assert_eq!(restored.perception(), &before, "{source}");
    }
}
