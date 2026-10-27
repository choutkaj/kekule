//! Builder placement on model drafts: hydrogens with coordinates, ideal
//! geometry, attaching atoms, and names.
use kekule::geometry::{Matrix3, Point3, RigidTransform, Vector3};
use kekule::hydrogens::{AddHydrogensOptions, AddHydrogensScope, RetainedHydrogenReason};
use kekule::molecule::{Atom, BondOrder, Element, ImplicitHydrogens};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::smiles;
use kekule::structure::placement::{
    fragment, fragments, ideal_bond_length, AtomShape, AttachmentSite, FragmentKind, PlacementError,
};
use kekule::structure::{Model, ModelEditor, Positions};
use kekule::topology::edit::{EditAtomId, EditResidueId};
use kekule::units::{Quantity, ANGSTROM, NANOMETER};
use std::collections::{BTreeMap, BTreeSet};

const TETRAHEDRAL: f64 = 109.471_220_634_490_7;

fn element(symbol: &str) -> Element {
    Element::from_symbol(symbol).unwrap()
}

fn point(editor: &ModelEditor, atom: EditAtomId) -> Point3 {
    editor
        .position(atom)
        .unwrap()
        .into_unit(ANGSTROM)
        .unwrap()
        .into_value()
}

fn distance(editor: &ModelEditor, a: EditAtomId, b: EditAtomId) -> f64 {
    (point(editor, a) - point(editor, b)).norm()
}

fn angle(editor: &ModelEditor, a: EditAtomId, center: EditAtomId, b: EditAtomId) -> f64 {
    let u = point(editor, a) - point(editor, center);
    let v = point(editor, b) - point(editor, center);
    (u.dot(v) / (u.norm() * v.norm()))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

fn ideal(a: Element, b: Element, order: BondOrder) -> f64 {
    ideal_bond_length(a, b, order)
        .unwrap()
        .into_unit(ANGSTROM)
        .unwrap()
        .into_value()
}

fn element_of(editor: &ModelEditor, atom: EditAtomId) -> Element {
    editor.topology().atom(atom).unwrap().element()
}

fn neighbors(editor: &ModelEditor, atom: EditAtomId) -> Vec<EditAtomId> {
    editor.topology().neighbors(atom).unwrap().collect()
}

fn hydrogen_count(editor: &ModelEditor, atom: EditAtomId) -> usize {
    neighbors(editor, atom)
        .into_iter()
        .filter(|n| element_of(editor, *n).atomic_number() == 1)
        .count()
}

fn heavy_atoms(editor: &ModelEditor) -> Vec<EditAtomId> {
    editor
        .topology()
        .atoms()
        .filter(|(_, atom)| atom.element().atomic_number() != 1)
        .map(|(id, _)| id)
        .collect()
}

/// Heavy atoms in dense order with their element, charge, and total
/// hydrogens; bonds between them; and the represented stereo count.
type Signature = (
    Vec<(u8, i8, usize)>,
    BTreeSet<(usize, usize, String)>,
    usize,
);

fn signature(model: &Model) -> Signature {
    let mut model = model.clone();
    model.perceive().unwrap();
    let topology = model.topology();
    let mut rank = BTreeMap::new();
    let mut atoms = Vec::new();
    for atom in topology.atoms() {
        if atom.element().atomic_number() != 1 {
            rank.insert(atom.index(), rank.len());
            atoms.push((
                atom.element().atomic_number(),
                atom.formal_charge(),
                atom.explicit_hydrogens() + atom.implicit_hydrogens().unwrap(),
            ));
        }
    }
    let bonds = topology
        .bonds()
        .filter_map(|bond| {
            let [a, b] = bond.atoms();
            let (a, b) = (rank.get(&a.index())?, rank.get(&b.index())?);
            Some((*a.min(b), *a.max(b), format!("{:?}", bond.order)))
        })
        .collect();
    let stereo = topology
        .definitions()
        .map(|definition| definition.molecule().stereo_elements().count())
        .sum();
    (atoms, bonds, stereo)
}

fn shape_angles(shape: AtomShape) -> &'static [f64] {
    match shape {
        AtomShape::Linear => &[180.0],
        AtomShape::TrigonalPlanar => &[120.0],
        AtomShape::Tetrahedral => &[TETRAHEDRAL],
        AtomShape::TrigonalBipyramidal => &[90.0, 120.0, 180.0],
        AtomShape::Octahedral => &[90.0, 180.0],
        _ => &[],
    }
}

/// PDB chemical components' ideal coordinates, with their hydrogens.
fn components() -> Vec<(String, Model)> {
    let document = kekule::sdf::parse_str(
        include_str!("../fixtures/placement/ccd-ideal.sdf"),
        Default::default(),
    )
    .unwrap();
    document
        .records()
        .iter()
        .map(|record| {
            let interpreted = record.interpret(Default::default()).unwrap();
            (record.title().to_owned(), interpreted.into_model())
        })
        .collect()
}

#[test]
fn completed_hydrogens_agree_with_valence_and_collapse_back() {
    for (name, model) in components() {
        let original = signature(&model);
        let mut editor = model.to_editor();
        let heavy = heavy_atoms(&editor);
        let expected = heavy
            .iter()
            .map(|atom| (*atom, hydrogen_count(&editor, *atom)))
            .collect::<BTreeMap<_, _>>();
        let total = expected.values().sum::<usize>();

        let removed = editor.remove_hydrogens(heavy.clone()).unwrap();
        assert_eq!(removed.removed_count(), total, "{name}");
        assert!(removed
            .atoms()
            .iter()
            .all(|atom| atom.retained().is_empty()));
        let stripped = signature(&editor.clone().finish().unwrap());

        let placed = editor
            .add_hydrogens(heavy.clone(), AddHydrogensOptions::default())
            .unwrap();
        assert_eq!(placed.added_count(), total, "{name}");
        for entry in placed.atoms() {
            let atom = entry.atom();
            assert_eq!(entry.hydrogens().len(), expected[&atom], "{name}");
            let shape = editor.atom_geometry(atom).unwrap().shape().unwrap();
            for &hydrogen in entry.hydrogens() {
                let length = ideal(element_of(&editor, atom), element("H"), BondOrder::Single);
                assert!((distance(&editor, atom, hydrogen) - length).abs() < 1e-9);
                // Every angle at the atom involving a new hydrogen is near an
                // ideal angle of its shape; between two new hydrogens or
                // around a single neighbor it is exact.
                let exact = neighbors(&editor, atom).len() - entry.hydrogens().len() <= 1;
                for other in neighbors(&editor, atom) {
                    if other == hydrogen {
                        continue;
                    }
                    let found = angle(&editor, hydrogen, atom, other);
                    let tolerance = if exact { 1e-6 } else { 8.0 };
                    assert!(
                        shape_angles(shape)
                            .iter()
                            .any(|ideal| (found - ideal).abs() < tolerance),
                        "{name}: {found} on {shape:?}"
                    );
                }
            }
        }
        // No new hydrogen clashes with an atom beyond its parent.
        let new = placed
            .atoms()
            .iter()
            .flat_map(|entry| entry.hydrogens().iter().map(|h| (*h, entry.atom())))
            .collect::<Vec<_>>();
        for (hydrogen, parent) in &new {
            for (other, _) in editor.topology().atoms() {
                if other != *hydrogen && other != *parent {
                    assert!(
                        distance(&editor, *hydrogen, other) > 1.5,
                        "{name}: clash at {:.2} Å",
                        distance(&editor, *hydrogen, other)
                    );
                }
            }
        }
        let completed = editor.clone().finish().unwrap();
        assert_eq!(signature(&completed), original, "{name}");
        let mut perceived = completed.clone();
        perceived.perceive().unwrap();

        editor.remove_hydrogens(heavy).unwrap();
        assert_eq!(signature(&editor.finish().unwrap()), stripped, "{name}");
    }
}

/// An isolated atom of `symbol` at the origin, completed with hydrogens.
fn completed(symbol: &str, charge: i8) -> (ModelEditor, EditAtomId, Vec<EditAtomId>) {
    let mut editor = ModelEditor::new();
    let atom = editor
        .place_atom(
            Atom::new(element(symbol)).with_formal_charge(charge),
            Quantity::new(Point3::origin(), NANOMETER),
        )
        .unwrap();
    let report = editor
        .add_hydrogens([atom], AddHydrogensOptions::default())
        .unwrap();
    let hydrogens = report.atoms()[0].hydrogens().to_vec();
    (editor, atom, hydrogens)
}

#[test]
fn selections_with_hydrogens_and_no_changes_are_safe() {
    let (mut editor, carbon, hydrogens) = completed("C", 0);
    let energy = PropertyKey::new("energy").unwrap();
    editor
        .properties_mut()
        .owner_mut()
        .insert(energy.clone(), PropertyValue::Int(7))
        .unwrap();
    // Completing a complete atom, or removing nothing, keeps annotations.
    let none = editor
        .add_hydrogens([carbon], AddHydrogensOptions::default())
        .unwrap();
    assert_eq!(none.added_count(), 0);
    assert_eq!(editor.remove_hydrogens([]).unwrap().removed_count(), 0);
    assert!(editor.properties().owner().get(&energy).is_some());
    // Selecting every atom reports the hydrogens under their parent only.
    let all = editor.topology().atom_ids().collect::<Vec<_>>();
    let report = editor.remove_hydrogens(all).unwrap();
    assert_eq!(report.atoms().len(), 1);
    assert_eq!(report.atoms()[0].atom(), carbon);
    assert_eq!(report.atoms()[0].removed(), hydrogens.as_slice());
    assert!(editor.properties().owner().get(&energy).is_none());
}

#[test]
fn hydrides_take_ideal_lengths_and_angles() {
    for (symbol, charge, count, shape) in [
        ("C", 0, 4, AtomShape::Tetrahedral),
        ("N", 0, 3, AtomShape::Tetrahedral),
        ("O", 0, 2, AtomShape::Tetrahedral),
        ("N", 1, 4, AtomShape::Tetrahedral),
        ("B", 0, 3, AtomShape::TrigonalPlanar),
        ("C", 1, 3, AtomShape::TrigonalPlanar),
        ("S", 0, 2, AtomShape::Tetrahedral),
        ("P", 0, 3, AtomShape::Tetrahedral),
        ("F", 0, 1, AtomShape::Tetrahedral),
    ] {
        let (editor, atom, hydrogens) = completed(symbol, charge);
        assert_eq!(hydrogens.len(), count, "{symbol}{charge}");
        assert_eq!(editor.atom_geometry(atom).unwrap().shape(), Some(shape));
        let length = ideal(element(symbol), element("H"), BondOrder::Single);
        for (i, a) in hydrogens.iter().enumerate() {
            assert!((distance(&editor, atom, *a) - length).abs() < 1e-9);
            for b in &hydrogens[i + 1..] {
                let found = angle(&editor, *a, atom, *b);
                assert!(
                    shape_angles(shape).iter().any(|i| (found - i).abs() < 1e-6),
                    "{symbol}: {found}"
                );
            }
        }
        // C gives methane, in a residue of its own.
        let model = editor.finish().unwrap();
        assert_eq!(model.atom_count(), count + 1);
        let residue = model.topology().residues().next().unwrap();
        assert_eq!(residue.name(), "UNL");
    }
}

/// A chain grown by attaching carbons toward a direction, then completed.
#[test]
fn chains_grow_anti_and_multiple_bonds_set_the_shape() {
    let mut editor = ModelEditor::new();
    let c1 = editor
        .place_atom(
            Atom::new(element("C")),
            Quantity::new(Point3::origin(), NANOMETER),
        )
        .unwrap();
    let toward = Vector3::new(1.0, 0.0, 0.0);
    let c2 = editor
        .attach_atom(
            c1,
            Atom::new(element("C")),
            BondOrder::Single,
            AttachmentSite::Toward(toward),
        )
        .unwrap()
        .atom();
    let c3 = editor
        .attach_atom(
            c2,
            Atom::new(element("C")),
            BondOrder::Single,
            AttachmentSite::Toward(toward),
        )
        .unwrap()
        .atom();
    let c4 = editor
        .attach_atom(
            c3,
            Atom::new(element("C")),
            BondOrder::Single,
            AttachmentSite::Toward(toward),
        )
        .unwrap()
        .atom();
    let cc = ideal(element("C"), element("C"), BondOrder::Single);
    for (a, b) in [(c1, c2), (c2, c3), (c3, c4)] {
        assert!((distance(&editor, a, b) - cc).abs() < 1e-9);
    }
    assert!((angle(&editor, c1, c2, c3) - TETRAHEDRAL).abs() < 1e-6);
    // The toward direction is nearest one site; the butane backbone is anti.
    let dihedral = {
        let [p1, p2, p3, p4] = [c1, c2, c3, c4].map(|a| point(&editor, a));
        let (b1, b2, b3) = (p2 - p1, p3 - p2, p4 - p3);
        let (n1, n2) = (b1.cross(b2), b2.cross(b3));
        let m = n1.cross(b2 / b2.norm());
        m.dot(n2).atan2(n1.dot(n2)).to_degrees()
    };
    assert!((dihedral.abs() - 180.0).abs() < 1e-6, "{dihedral}");
    editor
        .add_hydrogens([c1, c2, c3, c4], AddHydrogensOptions::default())
        .unwrap();
    assert_eq!(editor.topology().atom_count(), 14);

    // Ethylene: a double bond leaves three trigonal sites in one plane.
    let mut editor = ModelEditor::new();
    let a = editor
        .place_atom(
            Atom::new(element("C")),
            Quantity::new(Point3::origin(), NANOMETER),
        )
        .unwrap();
    let b = editor
        .attach_atom(
            a,
            Atom::new(element("C")),
            BondOrder::Double,
            AttachmentSite::Toward(toward),
        )
        .unwrap()
        .atom();
    let report = editor
        .add_hydrogens([a, b], AddHydrogensOptions::default())
        .unwrap();
    assert_eq!(report.added_count(), 4);
    assert!(
        (distance(&editor, a, b) - ideal(element("C"), element("C"), BondOrder::Double)).abs()
            < 1e-9
    );
    for entry in report.atoms() {
        let atom = entry.atom();
        assert_eq!(
            editor.atom_geometry(atom).unwrap().shape(),
            Some(AtomShape::TrigonalPlanar)
        );
        for hydrogen in entry.hydrogens() {
            let other = if atom == a { b } else { a };
            assert!((angle(&editor, *hydrogen, atom, other) - 120.0).abs() < 1e-6);
        }
    }
    let hydrogens = report
        .atoms()
        .iter()
        .flat_map(|entry| entry.hydrogens().iter().copied())
        .collect::<Vec<_>>();
    let normal = (point(&editor, b) - point(&editor, a))
        .cross(point(&editor, hydrogens[0]) - point(&editor, a));
    for h in &hydrogens {
        assert!((point(&editor, *h) - point(&editor, a)).dot(normal).abs() < 1e-9);
    }

    // Acetylene: linear.
    let mut editor = ModelEditor::new();
    let a = editor
        .place_atom(
            Atom::new(element("C")),
            Quantity::new(Point3::origin(), NANOMETER),
        )
        .unwrap();
    let b = editor
        .attach_atom(
            a,
            Atom::new(element("C")),
            BondOrder::Triple,
            AttachmentSite::Toward(toward),
        )
        .unwrap()
        .atom();
    let report = editor
        .add_hydrogens([a, b], AddHydrogensOptions::default())
        .unwrap();
    for entry in report.atoms() {
        let other = if entry.atom() == a { b } else { a };
        assert!((angle(&editor, entry.hydrogens()[0], entry.atom(), other) - 180.0).abs() < 1e-6);
    }
}

/// Positions of a molecule from SMILES with the center at the origin and the
/// given neighbor directions, in atom order after the center.
fn star(smiles_text: &str, directions: &[Vector3]) -> ModelEditor {
    let topology = smiles::to_topology(smiles_text).unwrap();
    let mut points = vec![Point3::origin()];
    for direction in directions {
        points.push(Point3::origin() + *direction * 0.16);
    }
    Model::new(
        topology,
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
    .into_editor()
}

#[test]
fn hypervalent_atoms_take_bipyramidal_and_octahedral_sites() {
    let x = Vector3::new(1.0, 0.0, 0.0);
    let y = Vector3::new(0.0, 1.0, 0.0);
    let z = Vector3::new(0.0, 0.0, 1.0);
    // Five-coordinate sulfur with one implicit hydrogen: octahedral.
    let mut editor = star("S(F)(F)(F)(F)F", &[x, -x, y, -y, z]);
    let sulfur = editor.topology().atom_ids().next().unwrap();
    assert_eq!(
        editor.atom_geometry(sulfur).unwrap().shape(),
        Some(AtomShape::Octahedral)
    );
    let report = editor
        .add_hydrogens([sulfur], AddHydrogensOptions::default())
        .unwrap();
    let hydrogen = report.atoms()[0].hydrogens()[0];
    let direction = point(&editor, hydrogen) - point(&editor, sulfur);
    assert!((direction / direction.norm() - (-z)).norm() < 1e-9);

    // Four-coordinate phosphorus with one implicit hydrogen: the hydrogen
    // takes the free equatorial site of a trigonal bipyramid.
    let (sin, cos) = (2.0 * std::f64::consts::FRAC_PI_3).sin_cos();
    let equatorial = [x, Vector3::new(cos, sin, 0.0)];
    let mut editor = star("P(F)(F)(F)F", &[z, -z, equatorial[0], equatorial[1]]);
    let phosphorus = editor.topology().atom_ids().next().unwrap();
    assert_eq!(
        editor.atom_geometry(phosphorus).unwrap().shape(),
        Some(AtomShape::TrigonalBipyramidal)
    );
    let report = editor
        .add_hydrogens([phosphorus], AddHydrogensOptions::default())
        .unwrap();
    let hydrogen = report.atoms()[0].hydrogens()[0];
    let direction = point(&editor, hydrogen) - point(&editor, phosphorus);
    let expected = Vector3::new(cos, -sin, 0.0);
    assert!((direction / direction.norm() - expected).norm() < 1e-9);
}

#[test]
fn atoms_over_their_valence_are_named() {
    let mut editor = star(
        "C(C)(C)(C)(C)C",
        &[
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, -1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ],
    );
    let atoms = editor.topology().atom_ids().collect::<Vec<_>>();
    let before = editor.topology().atom_count();
    assert_eq!(
        editor.add_hydrogens(atoms.clone(), AddHydrogensOptions::default()),
        Err(PlacementError::InvalidValence { atom: atoms[0] })
    );
    assert_eq!(editor.topology().atom_count(), before);
}

#[test]
fn attaching_replaces_hydrogens_and_keeps_stereo() {
    // L-alanine with explicit hydrogens, completed in a draft.
    let topology = smiles::to_topology("C[C@@H](C(=O)O)N").unwrap();
    let atoms = topology.atom_count();
    let ideal_positions = [
        Point3::new(-0.1, 0.13, 0.0),
        Point3::origin(),
        Point3::new(0.145, -0.05, 0.0),
        Point3::new(0.17, -0.17, 0.0),
        Point3::new(0.25, 0.04, 0.0),
        Point3::new(-0.05, -0.08, 0.12),
    ];
    assert_eq!(atoms, ideal_positions.len());
    let mut editor = Model::new(
        topology,
        Positions::new(Quantity::new(ideal_positions.to_vec(), NANOMETER)).unwrap(),
    )
    .unwrap()
    .into_editor();
    let ids = editor.topology().atom_ids().collect::<Vec<_>>();
    let center = ids[1];
    let report = editor
        .add_hydrogens(ids.clone(), AddHydrogensOptions::default())
        .unwrap();
    let alpha_hydrogen = report.atoms()[1].hydrogens()[0];
    let stereo_before = editor.clone().finish().unwrap();
    // Replace the alpha hydrogen by fluorine; stereo follows it.
    let attached = editor
        .attach_atom(
            center,
            Atom::new(element("F")),
            BondOrder::Single,
            AttachmentSite::Hydrogen(alpha_hydrogen),
        )
        .unwrap();
    assert_eq!(attached.replaced_hydrogen(), Some(alpha_hydrogen));
    assert!(editor.topology().atom(alpha_hydrogen).is_none());
    let direction = point(&editor, attached.atom()) - point(&editor, center);
    assert!((direction.norm() - ideal(element("C"), element("F"), BondOrder::Single)).abs() < 1e-9);
    let model = editor.finish().unwrap();
    let count = |model: &Model| {
        model
            .topology()
            .definitions()
            .map(|definition| definition.molecule().stereo_elements().count())
            .sum::<usize>()
    };
    assert_eq!(count(&stereo_before), 1);
    assert_eq!(count(&model), 1);

    // A heavy atom is not a replaceable hydrogen.
    let mut editor = model.into_editor();
    let ids = editor.topology().atom_ids().collect::<Vec<_>>();
    assert!(matches!(
        editor.attach_atom(
            ids[1],
            Atom::new(element("C")),
            BondOrder::Single,
            AttachmentSite::Hydrogen(ids[0])
        ),
        Err(PlacementError::NotReplaceableHydrogen { .. })
    ));
    // A saturated carbon has no free position.
    assert!(matches!(
        editor.attach_atom(
            ids[1],
            Atom::new(element("C")),
            BondOrder::Single,
            AttachmentSite::Toward(Vector3::new(1.0, 0.0, 0.0))
        ),
        Err(PlacementError::NoFreePosition { .. })
    ));
    assert_eq!(
        editor.attach_atom(
            ids[0],
            Atom::new(element("C")),
            BondOrder::Single,
            AttachmentSite::Toward(Vector3::zero())
        ),
        Err(PlacementError::InvalidDirection)
    );
}

#[test]
fn fixed_counts_are_materialized_and_restored() {
    let topology = smiles::to_topology("[CH3][NH2]").unwrap();
    let mut editor = Model::new(
        topology,
        Positions::new(Quantity::new(
            vec![Point3::origin(), Point3::new(0.147, 0.0, 0.0)],
            NANOMETER,
        ))
        .unwrap(),
    )
    .unwrap()
    .into_editor();
    let ids = editor.topology().atom_ids().collect::<Vec<_>>();
    let none = editor
        .add_hydrogens(
            ids.clone(),
            AddHydrogensOptions::default().with_scope(AddHydrogensScope::FixedOnly),
        )
        .unwrap();
    assert_eq!(none.added_count(), 5);
    for id in &ids {
        assert_eq!(
            editor.topology().atom(*id).unwrap().hydrogens(),
            ImplicitHydrogens::Fixed(0)
        );
    }
    // A hydrogen carrying a structure annotation is kept.
    let annotated = none.atoms()[0].hydrogens()[0];
    editor
        .properties_mut()
        .atoms_mut()
        .set_value(
            PropertyKey::new("charge").unwrap(),
            annotated,
            Some(PropertyValue::Int(1)),
        )
        .unwrap();
    let removed = editor.remove_hydrogens(ids.clone()).unwrap();
    assert_eq!(removed.removed_count(), 4);
    assert_eq!(
        removed.atoms()[0].retained(),
        &[(annotated, RetainedHydrogenReason::AtomProperties)]
    );
    assert_eq!(removed.atoms()[0].implicit_hydrogens(), 2);
    assert_eq!(removed.atoms()[1].implicit_hydrogens(), 2);
    assert_eq!(
        removed.atoms()[1].declaration(),
        ImplicitHydrogens::Fixed(2)
    );
}

/// Every atom name is unique within its residue.
fn assert_unique_names(model: &Model) {
    for residue in model.topology().residues() {
        let mut seen = BTreeSet::new();
        for site in residue.atom_sites() {
            let name = site.metadata().label_atom_id().unwrap().to_owned();
            assert!(
                seen.insert(name.clone()),
                "duplicate {name} in {}",
                residue.name()
            );
        }
    }
}

fn adenylate_kinase() -> Model {
    kekule::mmcif::parse_str(
        include_str!("../fixtures/mmcif/1AKE.cif"),
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap()
    .into_model()
}

/// Copies of `single` side by side, 10 nm apart, until 100,000 atoms.
fn large_model(single: &Model) -> Model {
    let mut large = ModelEditor::new();
    let mut offset = 0.0;
    while large.topology().atom_count() < 100_000 {
        let mut copy = single.clone();
        let shift = RigidTransform::new(
            Matrix3::default(),
            Quantity::new(Vector3::new(offset, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
        copy.apply_rigid_transform(&shift).unwrap();
        large.append_model(&copy).unwrap();
        offset += 10.0;
    }
    large.finish().unwrap()
}

/// The atoms of the first residue named `name`.
fn residue_atoms(editor: &ModelEditor, name: &str) -> (EditResidueId, Vec<EditAtomId>) {
    let topology = editor.topology();
    let residue = topology
        .residues()
        .find(|(_, residue)| residue.name() == name)
        .unwrap()
        .0;
    let atoms = topology
        .atom_sites()
        .filter(|(_, site)| site.residue() == residue)
        .map(|(_, site)| site.atom())
        .collect();
    (residue, atoms)
}

#[test]
fn new_atoms_are_named_uniquely_in_their_residues() {
    let mut editor = adenylate_kinase().into_editor();
    let (residue, ligand) = residue_atoms(&editor, "AP5");
    let report = editor
        .add_hydrogens(ligand.clone(), AddHydrogensOptions::default())
        .unwrap();
    assert!(report.added_count() > 0);
    let names = report
        .atoms()
        .iter()
        .flat_map(|entry| entry.hydrogens().iter().copied())
        .map(|hydrogen| {
            let site = editor
                .topology()
                .atom_site_for_atom(hydrogen)
                .unwrap()
                .unwrap();
            let site = editor.topology().atom_site(site).unwrap();
            assert_eq!(site.residue(), residue);
            site.metadata().label_atom_id().unwrap().to_owned()
        })
        .collect::<Vec<_>>();
    assert_eq!(names[0], "H1");
    assert_eq!(
        editor.next_atom_name(residue, element("H")).unwrap(),
        format!("H{}", names.len() + 1)
    );

    // Attach, remove, re-add and rename in one draft.
    let carbon = ligand
        .iter()
        .copied()
        .find(|atom| {
            element_of(&editor, *atom).atomic_number() == 6 && hydrogen_count(&editor, *atom) > 0
        })
        .unwrap();
    let hydrogen = neighbors(&editor, carbon)
        .into_iter()
        .find(|n| element_of(&editor, *n).atomic_number() == 1)
        .unwrap();
    let methyl = editor
        .attach_atom(
            carbon,
            Atom::new(element("C")),
            BondOrder::Single,
            AttachmentSite::Hydrogen(hydrogen),
        )
        .unwrap()
        .atom();
    editor
        .add_hydrogens([methyl], AddHydrogensOptions::default())
        .unwrap();
    editor.remove_hydrogens(ligand.clone()).unwrap();
    editor
        .add_hydrogens(ligand.clone(), AddHydrogensOptions::default())
        .unwrap();
    let nitrogen = ligand
        .iter()
        .copied()
        .find(|atom| element_of(&editor, *atom).atomic_number() == 6)
        .unwrap();
    let renamed = editor.rename_atom(nitrogen).unwrap().unwrap();
    assert!(renamed.starts_with('C'));
    let (_, protein) = residue_atoms(&editor, "LYS");
    editor
        .add_hydrogens(protein, AddHydrogensOptions::default())
        .unwrap();
    let model = editor.finish().unwrap();
    assert_unique_names(&model);
}

#[test]
fn placing_on_a_ligand_scales_with_the_ligand() {
    let single = adenylate_kinase();
    let large = large_model(&single);
    assert!(large.atom_count() >= 100_000);
    let time = |model: &Model| {
        let mut editor = model.to_editor();
        let (_, ligand) = residue_atoms(&editor, "AP5");
        let start = std::time::Instant::now();
        let report = editor
            .add_hydrogens(ligand.clone(), AddHydrogensOptions::default())
            .unwrap();
        editor.remove_hydrogens(ligand).unwrap();
        (start.elapsed(), report.added_count())
    };
    // Work on the ligand does not grow with the model it sits in.
    // The best of several runs, so load from other processes does not count.
    let best = |model: &Model| (0..5).map(|_| time(model)).min().unwrap();
    let (small, added_small) = best(&single);
    let (big, added_big) = best(&large);
    assert_eq!(added_small, added_big);
    assert!(
        big < small * 10 + std::time::Duration::from_millis(20),
        "{small:?} in one copy, {big:?} in 100k atoms"
    );
}

/// Graph distances, in bonds, from `start` within a draft, up to `limit`.
fn bond_distances(
    editor: &ModelEditor,
    start: EditAtomId,
    limit: usize,
) -> BTreeMap<EditAtomId, usize> {
    let mut distances = BTreeMap::from([(start, 0)]);
    let mut frontier = vec![start];
    for depth in 1..=limit {
        let mut next = Vec::new();
        for atom in frontier {
            for neighbor in neighbors(editor, atom) {
                if let std::collections::btree_map::Entry::Vacant(entry) = distances.entry(neighbor)
                {
                    entry.insert(depth);
                    next.push(neighbor);
                }
            }
        }
        frontier = next;
    }
    distances
}

#[test]
fn every_fragment_places_alone_with_its_ideal_geometry() {
    assert_eq!(fragments().len(), 48);
    for fragment in fragments() {
        let mut editor = ModelEditor::new();
        let target = Point3::new(1.0, 2.0, 3.0);
        let placed = editor
            .place_fragment(fragment, Quantity::new(target, NANOMETER))
            .unwrap();
        let molecule = fragment.molecule();
        assert_eq!(
            placed.atoms().len(),
            molecule.atom_count(),
            "{}",
            fragment.name()
        );
        let handles = placed.atoms().iter().copied().collect::<BTreeMap<_, _>>();
        // Bonds keep the library lengths; the atoms are centered on the point.
        for (_, bond) in molecule.bonds() {
            let library = (fragment.position(bond.a()).unwrap().into_value()
                - fragment.position(bond.b()).unwrap().into_value())
            .norm()
                * 10.0;
            let found = distance(&editor, handles[&bond.a()], handles[&bond.b()]);
            assert!((found - library).abs() < 1e-9, "{}", fragment.name());
        }
        let centroid = placed
            .atoms()
            .iter()
            .fold(Vector3::zero(), |sum, (_, atom)| {
                sum + (point(&editor, *atom) - Point3::origin())
            })
            / placed.atoms().len() as f64;
        assert!((centroid - Vector3::new(10.0, 20.0, 30.0)).norm() < 1e-9);
        let model = editor.finish().unwrap();
        let mut perceived = model.clone();
        perceived.perceive().unwrap();
        assert_unique_names(&model);
        assert_eq!(model.topology().residues().next().unwrap().name(), "UNL");
    }
    // Cyclohexane is a chair: ring torsions alternate near 55 degrees.
    let ring = fragment("cyclohexane").unwrap();
    let carbons = ring
        .molecule()
        .atoms()
        .filter(|(_, atom)| atom.element().atomic_number() == 6)
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    let mut cycle = vec![carbons[0]];
    while cycle.len() < 6 {
        let last = *cycle.last().unwrap();
        let next = ring
            .molecule()
            .neighbors(last)
            .unwrap()
            .find(|n| carbons.contains(n) && !cycle.contains(n))
            .unwrap();
        cycle.push(next);
    }
    let p = |atom| ring.position(atom).unwrap().into_value();
    let mut signs = Vec::new();
    for i in 0..6 {
        let [a, b, c, d] = [0, 1, 2, 3].map(|k| p(cycle[(i + k) % 6]));
        let (b1, b2, b3) = (b - a, c - b, d - c);
        let (n1, n2) = (b1.cross(b2), b2.cross(b3));
        let torsion = n1
            .cross(b2 / b2.norm())
            .dot(n2)
            .atan2(n1.dot(n2))
            .to_degrees();
        assert!((torsion.abs() - 55.0).abs() < 8.0, "{torsion}");
        signs.push(torsion > 0.0);
    }
    assert!(signs.windows(2).all(|pair| pair[0] != pair[1]));
}

#[test]
fn every_fragment_attaches_to_every_other_without_clashes() {
    for host in fragments() {
        for guest in fragments() {
            let mut editor = ModelEditor::new();
            let placed = editor
                .place_fragment(host, Quantity::new(Point3::origin(), NANOMETER))
                .unwrap();
            let handles = placed.atoms().iter().copied().collect::<BTreeMap<_, _>>();
            // The guest takes the place of the host's leaving atom: in place
            // of a leaving hydrogen, or toward where a leaving branch was.
            let parent = handles[&host.attachment_atom()];
            let leaving = handles[&host.leaving_atom()];
            let site = if element_of(&editor, leaving).atomic_number() == 1 {
                AttachmentSite::Hydrogen(leaving)
            } else {
                let toward = point(&editor, leaving) - point(&editor, parent);
                let branch = bond_distances(&editor, leaving, 1)
                    .into_keys()
                    .filter(|atom| *atom != parent);
                editor.delete_atoms(branch.collect::<Vec<_>>()).unwrap();
                AttachmentSite::Toward(toward)
            };
            let host_atoms = editor
                .topology()
                .atom_ids()
                .filter(|atom| *atom != leaving)
                .collect::<Vec<_>>();
            let attached = editor
                .attach_fragment(parent, guest, site)
                .unwrap_or_else(|e| panic!("{} on {}: {e:?}", guest.name(), host.name()));
            let label = format!("{} on {}", guest.name(), host.name());
            let length = ideal(
                element_of(&editor, parent),
                element_of(&editor, attached.attachment()),
                guest.attachment_order(),
            );
            assert!((distance(&editor, parent, attached.attachment()) - length).abs() < 1e-9);
            // Atoms more than two bonds apart across the new bond keep clear;
            // the closest such contacts (hydrogens of rigid rotors, as
            // tert-butyl peri to quinoline) are about 1.6 A.
            for (_, guest_atom) in attached.atoms() {
                let from_guest = bond_distances(&editor, *guest_atom, 2);
                for host_atom in &host_atoms {
                    if from_guest.contains_key(host_atom) {
                        continue;
                    }
                    let found = distance(&editor, *guest_atom, *host_atom);
                    assert!(found > 1.5, "{label}: {found:.2} A");
                }
            }
            let model = editor.finish().unwrap();
            assert_eq!(model.topology().instance_count(), 1, "{label}");
            let mut perceived = model.clone();
            perceived
                .perceive()
                .unwrap_or_else(|e| panic!("{label}: {e:?}"));
            assert_unique_names(&model);
        }
    }
}

#[test]
fn fragments_attach_toward_directions_and_through_other_atoms() {
    let mut editor = ModelEditor::new();
    let carbon = editor
        .place_atom(
            Atom::new(element("C")),
            Quantity::new(Point3::origin(), NANOMETER),
        )
        .unwrap();
    let phenyl = fragment("Phenyl").unwrap();
    assert_eq!(phenyl.kind(), FragmentKind::Group);
    let attached = editor
        .attach_fragment(
            carbon,
            phenyl,
            AttachmentSite::Toward(Vector3::new(0.0, 1.0, 0.0)),
        )
        .unwrap();
    assert_eq!(attached.atoms().len(), 11);
    assert_eq!(attached.replaced_hydrogen(), None);
    let direction = point(&editor, attached.attachment()) - point(&editor, carbon);
    assert!((direction.norm() - ideal(element("C"), element("C"), BondOrder::Single)).abs() < 1e-9);
    editor
        .add_hydrogens([carbon], AddHydrogensOptions::default())
        .unwrap();
    let mut model = editor.finish().unwrap();
    model.perceive().unwrap();
    assert_eq!(model.atom_count(), 15);

    // Pyridine through a carbon next to the nitrogen.
    let pyridine = fragment("pyridine").unwrap();
    let nitrogen = pyridine
        .molecule()
        .atoms()
        .find(|(_, atom)| atom.element().atomic_number() == 7)
        .unwrap()
        .0;
    let ortho = pyridine
        .molecule()
        .neighbors(nitrogen)
        .unwrap()
        .next()
        .unwrap();
    let ortho_pyridine = pyridine.with_attachment(ortho).unwrap();
    assert_eq!(ortho_pyridine.attachment_atom(), ortho);
    assert!(matches!(
        pyridine.with_attachment(nitrogen),
        Err(PlacementError::FragmentAtomWithoutHydrogen { .. })
    ));
    assert!(matches!(
        fragment("adamantane"),
        Err(PlacementError::UnknownFragment { .. })
    ));
}

#[test]
fn attaching_a_fragment_to_a_ligand_scales_with_the_ligand() {
    let single = adenylate_kinase();
    let large = large_model(&single);
    let phenyl = fragment("phenyl").unwrap();
    let time = |model: &Model| {
        let mut editor = model.to_editor();
        let (_, ligand) = residue_atoms(&editor, "AP5");
        editor
            .add_hydrogens(ligand.clone(), AddHydrogensOptions::default())
            .unwrap();
        let carbon = ligand
            .iter()
            .copied()
            .find(|atom| {
                element_of(&editor, *atom).atomic_number() == 6
                    && hydrogen_count(&editor, *atom) > 0
            })
            .unwrap();
        let hydrogens = neighbors(&editor, carbon)
            .into_iter()
            .filter(|n| element_of(&editor, *n).atomic_number() == 1)
            .collect::<Vec<_>>();
        // A first attachment deletes a hydrogen; the second, timed one then
        // runs in a draft with a deleted atom.
        let other = ligand
            .iter()
            .copied()
            .rev()
            .find(|atom| {
                *atom != carbon
                    && element_of(&editor, *atom).atomic_number() == 6
                    && hydrogen_count(&editor, *atom) > 0
            })
            .unwrap();
        let other_hydrogen = neighbors(&editor, other)
            .into_iter()
            .find(|n| element_of(&editor, *n).atomic_number() == 1)
            .unwrap();
        editor
            .attach_fragment(other, phenyl, AttachmentSite::Hydrogen(other_hydrogen))
            .unwrap();
        let start = std::time::Instant::now();
        editor
            .attach_fragment(carbon, phenyl, AttachmentSite::Hydrogen(hydrogens[0]))
            .unwrap();
        start.elapsed()
    };
    // The best of several runs, so load from other processes does not count.
    let best = |model: &Model| (0..5).map(|_| time(model)).min().unwrap();
    let (small, big) = (best(&single), best(&large));
    assert!(
        big < small * 10 + std::time::Duration::from_millis(20),
        "{small:?} in one copy, {big:?} in 100k atoms"
    );
}

#[test]
fn fragments_on_linear_parents_turn_clear_of_their_surroundings() {
    // Phenyl on an acetylene carbon, with an atom an angstrom beside where
    // the ring's first torsion would put one of its carbons.
    let mut editor = ModelEditor::new();
    let placed = editor
        .place_fragment(
            fragment("ethynyl").unwrap(),
            Quantity::new(Point3::origin(), NANOMETER),
        )
        .unwrap();
    let ethynyl = fragment("ethynyl").unwrap();
    let handles = placed.atoms().iter().copied().collect::<BTreeMap<_, _>>();
    let parent = handles[&ethynyl.attachment_atom()];
    let hydrogen = handles[&ethynyl.leaving_atom()];
    assert_eq!(
        editor.atom_geometry(parent).unwrap().shape(),
        Some(AtomShape::Linear)
    );
    let phenyl = fragment("phenyl").unwrap();
    let mut trial = editor.clone();
    let first = trial
        .attach_fragment(parent, phenyl, AttachmentSite::Hydrogen(hydrogen))
        .unwrap();
    let ring_carbon = first
        .atoms()
        .iter()
        .map(|(_, atom)| *atom)
        .find(|atom| *atom != first.attachment() && element_of(&trial, *atom).atomic_number() == 6)
        .unwrap();
    // One angstrom outward from the ring carbon, away from the bond axis.
    let (start, end, carbon) = (
        point(&trial, parent),
        point(&trial, first.attachment()),
        point(&trial, ring_carbon),
    );
    let axis = (end - start) / (end - start).norm();
    let foot = start + axis * (carbon - start).dot(axis);
    let outward = (carbon - foot) / (carbon - foot).norm();
    let obstacle_point = Quantity::new(carbon + outward, ANGSTROM);
    let obstacle = editor
        .place_atom(Atom::new(element("Ar")), obstacle_point)
        .unwrap();
    let attached = editor
        .attach_fragment(parent, phenyl, AttachmentSite::Hydrogen(hydrogen))
        .unwrap();
    let closest = attached
        .atoms()
        .iter()
        .map(|(_, atom)| distance(&editor, *atom, obstacle))
        .fold(f64::INFINITY, f64::min);
    assert!(closest > 2.0, "{closest:.2} A");
}

#[test]
fn ring_fragments_follow_their_numbering() {
    // Heteroatoms at their locants, as 0-based ring positions.
    let heteroatoms: &[(&str, &[(usize, &str)])] = &[
        ("pyridine", &[(0, "N")]),
        ("pyrimidine", &[(0, "N"), (2, "N")]),
        ("pyrazine", &[(0, "N"), (3, "N")]),
        ("pyrrole", &[(0, "N")]),
        ("furan", &[(0, "O")]),
        ("thiophene", &[(0, "S")]),
        ("imidazole", &[(0, "N"), (2, "N")]),
        ("pyrazole", &[(0, "N"), (1, "N")]),
        ("oxazole", &[(0, "O"), (2, "N")]),
        ("thiazole", &[(0, "S"), (2, "N")]),
        ("piperidine", &[(0, "N")]),
        ("piperazine", &[(0, "N"), (3, "N")]),
        ("morpholine", &[(0, "O"), (3, "N")]),
        ("tetrahydrofuran", &[(0, "O")]),
        ("indole", &[(0, "N")]),
        ("quinoline", &[(0, "N")]),
        ("benzimidazole", &[(0, "N"), (2, "N")]),
    ];
    // Attachments at their documented locants, as 0-based ring positions.
    let attachments: &[(&str, usize)] = &[
        ("pyridine", 3),
        ("pyrimidine", 1),
        ("pyrazine", 1),
        ("pyrrole", 0),
        ("furan", 1),
        ("thiophene", 1),
        ("imidazole", 0),
        ("pyrazole", 0),
        ("oxazole", 1),
        ("thiazole", 1),
        ("piperidine", 0),
        ("piperazine", 0),
        ("morpholine", 3),
        ("tetrahydrofuran", 1),
        ("naphthalene", 1),
        ("indole", 2),
        ("quinoline", 3),
        ("benzimidazole", 1),
    ];
    for fragment in fragments() {
        let molecule = fragment.molecule();
        let order = molecule.atom_ids().collect::<Vec<_>>();
        let heavy = |atom| molecule.atom(atom).unwrap().element().atomic_number() != 1;
        let ring = order.iter().take_while(|atom| heavy(**atom)).count();
        // Hydrogens follow the heavy atoms, in their parents' order.
        assert!(
            order[ring..].iter().all(|atom| !heavy(*atom)),
            "{}",
            fragment.name()
        );
        let parents = order[ring..]
            .iter()
            .map(|h| {
                let parent = molecule.neighbors(*h).unwrap().next().unwrap();
                order.iter().position(|atom| *atom == parent).unwrap()
            })
            .collect::<Vec<_>>();
        assert!(
            parents.windows(2).all(|pair| pair[0] <= pair[1]),
            "{}",
            fragment.name()
        );
        if fragment.kind() == FragmentKind::Group {
            if fragment.name() != "phenyl" {
                assert_eq!(order[0], fragment.attachment_atom(), "{}", fragment.name());
            }
            continue;
        }
        // The ring atoms come first, each bonded to the next and the last to
        // the first: around the ring, or along the fused-ring numbering.
        for i in 0..ring {
            let (a, b) = (order[i], order[(i + 1) % ring]);
            assert!(
                molecule.bond_between(a, b).is_some(),
                "{}: ring positions {} and {} are not bonded",
                fragment.name(),
                i + 1,
                (i + 1) % ring + 1
            );
        }
        let symbol = |index: usize| molecule.atom(order[index]).unwrap().element().symbol();
        let expected = heteroatoms
            .iter()
            .find(|(name, _)| *name == fragment.name())
            .map_or(&[][..], |(_, atoms)| *atoms);
        for index in 0..ring {
            let element = expected
                .iter()
                .find(|(position, _)| *position == index)
                .map_or("C", |(_, symbol)| *symbol);
            assert_eq!(
                symbol(index),
                element,
                "{} position {}",
                fragment.name(),
                index + 1
            );
        }
        let attachment = attachments
            .iter()
            .find(|(name, _)| *name == fragment.name())
            .map_or(0, |(_, position)| *position);
        assert_eq!(
            order[attachment],
            fragment.attachment_atom(),
            "{}",
            fragment.name()
        );
    }

    // Placed benzene is named around its ring: C1 and C4 are para.
    let mut editor = ModelEditor::new();
    let placed = editor
        .place_fragment(
            fragment("benzene").unwrap(),
            Quantity::new(Point3::origin(), NANOMETER),
        )
        .unwrap();
    let name = |atom: EditAtomId| {
        let site = editor.topology().atom_site_for_atom(atom).unwrap().unwrap();
        editor
            .topology()
            .atom_site(site)
            .unwrap()
            .metadata()
            .label_atom_id()
            .unwrap()
            .to_owned()
    };
    let carbons = placed
        .atoms()
        .iter()
        .map(|(_, atom)| *atom)
        .filter(|atom| element_of(&editor, *atom).atomic_number() == 6)
        .collect::<Vec<_>>();
    let names = carbons.iter().map(|atom| name(*atom)).collect::<Vec<_>>();
    assert_eq!(names, ["C1", "C2", "C3", "C4", "C5", "C6"]);
    for i in 0..6 {
        assert!(neighbors(&editor, carbons[i]).contains(&carbons[(i + 1) % 6]));
    }
    assert_eq!(bond_distances(&editor, carbons[0], 3)[&carbons[3]], 3);
}
