use super::*;
use crate::structure::placement::fragment;

fn candidate(deviation: f64, clearance: f64, crowding: f64) -> Candidate {
    Candidate {
        orientation: 0,
        rotation: Rotation::about(Vector3::new(0.0, 0.0, 1.0), 0.0),
        deviation,
        clearance,
        crowding,
        hydrogens: Default::default(),
    }
}

#[test]
fn the_site_hydrogen_goes_first_then_implicit_then_the_nearest() {
    let toward = |hydrogen: usize| hydrogen as f64;
    // The new bond takes an explicit hydrogen's site.
    assert_eq!(given_up(&[1, 2], 0, 1, Some(1), toward), (Some(1), vec![]));
    // A second unit takes the other explicit hydrogen.
    assert_eq!(given_up(&[1, 2], 0, 2, Some(1), toward), (Some(1), vec![2]));
    // A free site takes an implicit hydrogen, and so does the next unit.
    assert_eq!(given_up::<usize>(&[], 2, 2, None, toward), (None, vec![]));
    // With the implicit one spent, an explicit one follows.
    assert_eq!(given_up(&[7], 1, 2, None, toward), (None, vec![7]));
    // Without free sites or implicit hydrogens, the nearest explicit one.
    assert_eq!(given_up(&[7, 8, 3], 0, 2, None, toward), (Some(8), vec![7]));
}

#[test]
fn the_least_crowded_of_the_aligned_clear_turns_wins() {
    let tolerance = ALIGNMENT_TOLERANCE;
    let candidates = [
        // Best aligned, but in contact.
        candidate(0.0, CLEAR_CONTACT * 0.5, 0.0),
        candidate(0.01, 0.25, 0.5),
        candidate(0.01 + tolerance * 0.5, 0.25, 0.1),
        // Least crowded, but too far off its sites.
        candidate(0.02 + tolerance, 0.3, 0.0),
    ];
    assert_eq!(choose(&candidates), Some(2));
    // Ties go to the earlier.
    let tied = [candidate(0.0, 0.25, 0.1), candidate(0.0, 0.25, 0.1)];
    assert_eq!(choose(&tied), Some(0));
    // None keeps clear: no choice, and the clearest is the fallback.
    let crowded = [
        candidate(0.0, 0.1, 1.0),
        candidate(0.0, 0.15, 1.0),
        candidate(0.0, 0.15, 1.0),
    ];
    assert_eq!(choose(&crowded), None);
    assert_eq!(clearest(&crowded), 1);
}

fn mirror_of(name: &str) -> Option<(Vector3, Vector3)> {
    let fragment = fragment(name).unwrap();
    let molecule = fragment.molecule();
    let bond = molecule.bond(fragment.fusion_bond().unwrap()).unwrap();
    let fused = [bond.a(), bond.b()];
    let dropped = fusion_dropped(molecule, fused);
    let kept = molecule
        .atom_ids()
        .filter(|atom| !dropped.contains(atom))
        .collect::<Vec<_>>();
    let points = molecule
        .atom_ids()
        .map(|atom| (atom, fragment_point(fragment, atom)))
        .collect::<BTreeMap<_, _>>();
    let axis = unit(points[&fused[1]] - points[&fused[0]]).unwrap();
    mirror_plane(molecule, &points, &kept, fused).map(|normal| (normal, axis))
}

#[test]
fn a_flat_ring_has_no_mirror_image_and_a_chair_has_one() {
    assert!(mirror_of("benzene").is_none());
    assert!(mirror_of("pyridine").is_none());
    let (normal, axis) = mirror_of("cyclohexane").unwrap();
    assert!((normal.norm() - 1.0).abs() < 1e-12);
    assert!(normal.dot(axis).abs() < 1e-12);
}

/// (E)-2-butene in the xy plane, its methyls opposite.
fn trans_butene() -> StructureEditor {
    crate::structure::Structure::new(
        crate::smiles::to_topology("C/C=C/C").unwrap(),
        crate::structure::Positions::new(crate::units::Quantity::new(
            vec![
                Point3::new(-0.075, 0.13, 0.0),
                Point3::origin(),
                Point3::new(0.134, 0.0, 0.0),
                Point3::new(0.209, -0.13, 0.0),
            ],
            crate::units::NANOMETER,
        ))
        .unwrap(),
    )
    .unwrap()
    .into_editor()
}

#[test]
fn a_refused_fusion_leaves_the_revision_and_stereo_alone() {
    let mut editor = trans_butene();
    let bond = editor
        .topology
        .bonds()
        .find(|(_, bond)| bond.order == BondOrder::Double)
        .unwrap()
        .0;
    let revision = editor.topology.structural_revision;
    let stereo = |editor: &StructureEditor| {
        editor
            .clone()
            .finish()
            .unwrap()
            .topology()
            .definitions()
            .flat_map(|definition| {
                definition
                    .molecule()
                    .stereo_elements()
                    .map(|(_, element)| element.clone())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let before = stereo(&editor);
    assert_eq!(before.len(), 1);
    for name in ["benzene", "methyl"] {
        assert!(editor.fuse_fragment(bond, fragment(name).unwrap()).is_err());
        assert_eq!(editor.topology.structural_revision, revision, "{name}");
        assert!(editor.topology.pruned_stereo().is_empty(), "{name}");
        assert_eq!(stereo(&editor), before, "{name}");
    }
}
