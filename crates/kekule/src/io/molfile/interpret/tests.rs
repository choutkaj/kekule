use super::*;
use crate::stereo::{AxisOrientation, AxisStereo, StereoCarrier};

#[test]
fn enhanced_group_rejects_an_endpoint_shared_by_two_axes() {
    // This represented graph tests source-member identity, independently
    // of the chemical stability or coordinate perception of either axis.
    let mut editor = MoleculeEditor::new();
    let atoms = (0..5)
        .map(|_| {
            editor
                .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
                .unwrap()
        })
        .collect::<Vec<_>>();
    let bonds = atoms
        .windows(2)
        .map(|pair| {
            editor
                .add_bond(pair[0], pair[1], BondOrder::Single)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let mut axes = Vec::new();
    for (bond, left, right) in [
        (bonds[1], atoms[0], atoms[3]),
        (bonds[2], atoms[1], atoms[4]),
    ] {
        axes.push(
            editor
                .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
                    axis: bond,
                    carriers: [StereoCarrier::Atom(left), StereoCarrier::Atom(right)],
                    orientation: Some(AxisOrientation::Clockwise),
                })))
                .unwrap(),
        );
    }
    let mut molecule = editor.finish().unwrap().into_editor();
    assert_eq!(
        molfile_stereo_group_members_at_atom(&molecule, atoms[2]),
        axes
    );
    assert_eq!(
        molfile_stereo_group_members_at_atom(&molecule, atoms[1]),
        [axes[0]]
    );
    assert_eq!(
        molfile_stereo_group_members_at_atom(&molecule, atoms[3]),
        [axes[1]]
    );
    let atom_map = atoms.iter().copied().map(|atom| (atom, atom)).collect();
    let mut group = MolfileStereoGroup {
        kind: MolfileStereoGroupKind::And(1),
        atoms: vec![atoms[2].index()],
        location: SourceLocation::at_line(73),
    };
    let error =
        install_molfile_stereo_groups(&mut molecule, &atom_map, &[group.clone()]).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Chemistry);
    assert_eq!(error.location(), Some(SourceLocation::at_line(73)));
    assert_eq!(
        error.message(),
        "A V3000 stereo group atom is not one stereocenter or atropisomeric axis."
    );
    assert_eq!(molecule.stereo_groups().count(), 0);

    // The other endpoint of each axis identifies that member uniquely.
    group.atoms = vec![atoms[1].index(), atoms[3].index()];
    install_molfile_stereo_groups(&mut molecule, &atom_map, &[group]).unwrap();
    assert_eq!(molecule.stereo_groups().next().unwrap().1.members, axes);
}
