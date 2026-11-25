//! Builder chemistry on structure drafts: builder valences, settling an
//! atom's hydrogens and angle after an edit, ideal bond lengths, and
//! multiple bonds attached in place of hydrogens.
use kekule::math::{Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, ImplicitHydrogens};
use kekule::structure::placement::{
    builder_charges, builder_hydrogen_count, builder_valences, fragment, ideal_bond_length,
    AttachmentSite, BuilderValenceError, PlacementError, SettledAtomsReport,
};
use kekule::structure::StructureEditor;
use kekule::topology::edit::{EditAtomId, EditBondId};
use kekule::units::{Quantity, ANGSTROM};

const TETRAHEDRAL: f64 = 109.471_220_634_490_7;
/// Lengths agree with the ideal ones within this, in Å.
const LENGTH: f64 = 0.02;
/// Angles agree with the ideal ones within this, in degrees.
const ANGLE: f64 = 3.0;

fn element(symbol: &str) -> Element {
    Element::from_symbol(symbol).unwrap()
}

fn point(editor: &StructureEditor, atom: EditAtomId) -> Point3 {
    editor
        .position(atom)
        .unwrap()
        .into_unit(ANGSTROM)
        .unwrap()
        .into_value()
}

fn distance(editor: &StructureEditor, a: EditAtomId, b: EditAtomId) -> f64 {
    (point(editor, a) - point(editor, b)).norm()
}

fn angle(editor: &StructureEditor, a: EditAtomId, center: EditAtomId, b: EditAtomId) -> f64 {
    let u = point(editor, a) - point(editor, center);
    let v = point(editor, b) - point(editor, center);
    u.cross(v).norm().atan2(u.dot(v)).to_degrees()
}

fn ideal(a: Element, b: Element, order: BondOrder) -> f64 {
    ideal_bond_length(a, b, order)
        .unwrap()
        .into_unit(ANGSTROM)
        .unwrap()
        .into_value()
}

fn atom_of(editor: &StructureEditor, atom: EditAtomId) -> Atom {
    editor.topology().atom(atom).unwrap().clone()
}

fn is_hydrogen(editor: &StructureEditor, atom: EditAtomId) -> bool {
    atom_of(editor, atom).element() == Element::H
}

fn hydrogens(editor: &StructureEditor, atom: EditAtomId) -> Vec<EditAtomId> {
    editor
        .topology()
        .neighbors(atom)
        .unwrap()
        .filter(|neighbor| is_hydrogen(editor, *neighbor))
        .collect()
}

fn heavy_neighbors(editor: &StructureEditor, atom: EditAtomId) -> Vec<EditAtomId> {
    editor
        .topology()
        .neighbors(atom)
        .unwrap()
        .filter(|neighbor| !is_hydrogen(editor, *neighbor))
        .collect()
}

fn bond(editor: &StructureEditor, a: EditAtomId, b: EditAtomId) -> EditBondId {
    editor.topology().bond_between(a, b).unwrap().unwrap()
}

/// A lone atom of `symbol` at the origin, settled with its hydrogens.
fn built(symbol: &str) -> (StructureEditor, EditAtomId) {
    let mut editor = StructureEditor::new();
    let atom = editor
        .place_atom(
            Atom::new(element(symbol)),
            Quantity::new(Point3::origin(), ANGSTROM),
        )
        .unwrap();
    editor.settle_atoms([atom]).unwrap();
    (editor, atom)
}

/// Attaches a single-bonded `symbol` in place of the first hydrogen of
/// `parent` and settles both, as a builder does.
fn grow(editor: &mut StructureEditor, parent: EditAtomId, symbol: &str) -> EditAtomId {
    let hydrogen = hydrogens(editor, parent)[0];
    let atom = editor
        .attach_atom(
            parent,
            Atom::new(element(symbol)),
            BondOrder::Single,
            AttachmentSite::Hydrogen(hydrogen),
        )
        .unwrap()
        .atom();
    editor.settle_atoms([parent, atom]).unwrap();
    atom
}

/// Changes the bond between `a` and `b` to `order` and settles it and both
/// atoms, as a builder does.
fn reorder(editor: &mut StructureEditor, a: EditAtomId, b: EditAtomId, order: BondOrder) {
    let id = bond(editor, a, b);
    editor.set_bond_order(id, order).unwrap();
    editor.settle_bonds([id]).unwrap();
    editor.settle_atoms([a, b]).unwrap();
}

fn charge(editor: &mut StructureEditor, atom: EditAtomId, formal_charge: i8) -> SettledAtomsReport {
    let value = atom_of(editor, atom).with_formal_charge(formal_charge);
    editor.replace_atom(atom, value).unwrap();
    editor.settle_atoms([atom]).unwrap()
}

fn change_element(editor: &mut StructureEditor, atom: EditAtomId, symbol: &str) {
    let value = atom_of(editor, atom).with_element(element(symbol));
    editor.replace_atom(atom, value).unwrap();
    editor.settle_atoms([atom]).unwrap();
}

/// Asserts that `atom` has `count` hydrogens at the ideal length and that
/// every angle between its neighbors is near one of `angles`.
fn assert_settled(editor: &StructureEditor, atom: EditAtomId, count: usize, angles: &[f64]) {
    let own = hydrogens(editor, atom);
    let symbol = atom_of(editor, atom).element().symbol();
    assert_eq!(own.len(), count, "{symbol} hydrogens");
    let length = ideal(
        atom_of(editor, atom).element(),
        Element::H,
        BondOrder::Single,
    );
    for hydrogen in &own {
        assert!(
            (distance(editor, atom, *hydrogen) - length).abs() < LENGTH,
            "{symbol}–H {}",
            distance(editor, atom, *hydrogen)
        );
    }
    let neighbors = editor
        .topology()
        .neighbors(atom)
        .unwrap()
        .collect::<Vec<_>>();
    for (index, a) in neighbors.iter().enumerate() {
        for b in &neighbors[index + 1..] {
            let found = angle(editor, *a, atom, *b);
            assert!(
                angles.iter().any(|ideal| (found - ideal).abs() < ANGLE),
                "{symbol}: {found:.1}° is none of {angles:?}"
            );
        }
    }
    assert_eq!(
        editor.topology().atom(atom).unwrap().hydrogens(),
        ImplicitHydrogens::Inferred
    );
}

fn assert_planar(editor: &StructureEditor, atom: EditAtomId) {
    let neighbors = editor
        .topology()
        .neighbors(atom)
        .unwrap()
        .collect::<Vec<_>>();
    let sum = neighbors
        .iter()
        .enumerate()
        .flat_map(|(index, a)| {
            neighbors[index + 1..]
                .iter()
                .map(move |b| angle(editor, *a, atom, *b))
        })
        .sum::<f64>();
    assert!((sum - 360.0).abs() < ANGLE, "angles sum to {sum}");
}

fn assert_length(editor: &StructureEditor, a: EditAtomId, b: EditAtomId, order: BondOrder) {
    let expected = ideal(
        atom_of(editor, a).element(),
        atom_of(editor, b).element(),
        order,
    );
    let found = distance(editor, a, b);
    assert!(
        (found - expected).abs() < LENGTH,
        "{found} instead of {expected}"
    );
}

#[test]
fn builder_valences_list_charges_and_hydrogen_counts() {
    assert_eq!(builder_charges(Element::C), Some(&[-1, 0, 1][..]));
    assert_eq!(builder_charges(Element::H), Some(&[0][..]));
    assert_eq!(builder_charges(element("Si")), Some(&[0][..]));
    assert_eq!(builder_charges(element("Se")), Some(&[-1, 0, 1][..]));
    assert_eq!(builder_charges(element("Cl")), Some(&[-1, 0][..]));
    assert_eq!(builder_charges(element("Na")), None);
    assert_eq!(builder_valences(Element::P, 0), Some(&[3, 5][..]));
    assert_eq!(builder_valences(Element::B, -1), Some(&[4][..]));
    assert_eq!(builder_valences(Element::F, -1), Some(&[0][..]));
    assert_eq!(builder_valences(Element::C, 2), None);
    let counts = [
        (Element::C, 0, 0, 4),
        (Element::C, 1, 0, 3),
        (Element::C, -1, 0, 3),
        (Element::N, 1, 1, 3),
        (Element::N, -1, 0, 2),
        (Element::O, 1, 0, 3),
        (Element::O, -1, 1, 0),
        (Element::B, 0, 0, 3),
        (Element::B, -1, 0, 4),
        (Element::S, 0, 2, 0),
        (Element::S, 0, 3, 1),
        (Element::S, 0, 5, 1),
        (Element::S, 1, 2, 1),
        (Element::P, 0, 4, 1),
        (Element::P, 1, 0, 4),
        (Element::F, -1, 0, 0),
        (Element::H, 0, 0, 1),
    ];
    for (element, charge, valence, hydrogens) in counts {
        assert_eq!(
            builder_hydrogen_count(element, charge, valence),
            Ok(hydrogens),
            "{} {charge} {valence}",
            element.symbol()
        );
    }
    assert_eq!(
        builder_hydrogen_count(Element::O, 2, 0),
        Err(BuilderValenceError::UnsupportedCharge {
            element: Element::O,
            formal_charge: 2,
            charges: &[-1, 0, 1],
        })
    );
    assert!(matches!(
        builder_hydrogen_count(Element::N, 0, 4),
        Err(BuilderValenceError::ExcessValence { valence: 4, .. })
    ));
}

#[test]
fn errors_name_the_state_and_what_it_takes() {
    let oxygen = builder_hydrogen_count(Element::O, 2, 0).unwrap_err();
    assert_eq!(
        oxygen.to_string(),
        "The atom is O with charge +2. O takes charge -1, 0 or +1."
    );
    let carbon = builder_hydrogen_count(Element::C, 0, 5).unwrap_err();
    assert_eq!(
        carbon.to_string(),
        "The atom has valence 5. C takes valence 4."
    );
    let sulfur = builder_hydrogen_count(Element::S, 1, 4).unwrap_err();
    assert_eq!(
        sulfur.to_string(),
        "The atom has valence 4. S with charge +1 takes valence 3."
    );
    let sulfur = builder_hydrogen_count(Element::S, 0, 7).unwrap_err();
    assert_eq!(
        sulfur.to_string(),
        "The atom has valence 7. S takes valence 2, 4 or 6."
    );
    let iron = builder_hydrogen_count(Element::from_symbol("Fe").unwrap(), 0, 0).unwrap_err();
    assert_eq!(
        iron.to_string(),
        "The atom is Fe. Fe has no builder valences."
    );
}

#[test]
fn methane_charged_and_back_keeps_its_shape() {
    let (mut editor, carbon) = built("C");
    assert_settled(&editor, carbon, 4, &[TETRAHEDRAL]);
    let first = hydrogens(&editor, carbon);

    let report = charge(&mut editor, carbon, 1);
    assert_eq!(report.removed_count(), 1);
    assert_settled(&editor, carbon, 3, &[120.0]);
    assert_planar(&editor, carbon);
    // The hydrogens that stay keep their identity.
    assert!(hydrogens(&editor, carbon).iter().all(|h| first.contains(h)));

    let report = charge(&mut editor, carbon, -1);
    assert_eq!((report.added_count(), report.removed_count()), (0, 0));
    assert!(!report.atoms()[0].moved_hydrogens().is_empty());
    assert_settled(&editor, carbon, 3, &[TETRAHEDRAL]);

    let report = charge(&mut editor, carbon, 0);
    assert_eq!(report.added_count(), 1);
    assert_settled(&editor, carbon, 4, &[TETRAHEDRAL]);

    // Settling again changes nothing.
    let report = editor.settle_atoms([carbon]).unwrap();
    let entry = &report.atoms()[0];
    assert!(entry.added_hydrogens().is_empty() && entry.removed_hydrogens().is_empty());
    assert!(entry.moved_hydrogens().is_empty());
    // A carbon with four explicit hydrogens has no free position, and one
    // that lost a hydrogen gains it back.
    assert!(editor
        .free_directions(carbon, BondOrder::Single)
        .unwrap()
        .is_empty());
    editor.delete_atom(hydrogens(&editor, carbon)[0]).unwrap();
    let report = editor.settle_atoms([carbon]).unwrap();
    assert_eq!(report.added_count(), 1);
    assert_settled(&editor, carbon, 4, &[TETRAHEDRAL]);
    let structure = editor.finish().unwrap();
    assert_eq!(structure.atom_count(), 5);
}

#[test]
fn ring_atoms_report_angles_they_cannot_take() {
    let mut editor = StructureEditor::new();
    let placed = editor
        .place_fragment(
            fragment("benzene").unwrap(),
            Quantity::new(Point3::origin(), ANGSTROM),
        )
        .unwrap();
    let carbon = placed
        .atoms()
        .iter()
        .map(|(_, atom)| *atom)
        .find(|atom| !is_hydrogen(&editor, *atom))
        .unwrap();
    // A phenyl cation's carbon is linear by its domains, which the ring
    // does not allow.
    let report = charge(&mut editor, carbon, 1);
    assert!(report.atoms()[0].angle_in_ring());
    assert!(report.atoms()[0].turned_atoms().is_empty());
    assert_eq!(report.removed_count(), 1);
    let report = charge(&mut editor, carbon, 0);
    assert!(!report.atoms()[0].angle_in_ring());
    assert_settled(&editor, carbon, 1, &[120.0]);
    assert_planar(&editor, carbon);
    // Pyridine's nitrogen takes no hydrogen.
    change_element(&mut editor, carbon, "N");
    assert_settled(&editor, carbon, 0, &[120.0]);
    editor.finish().unwrap();
}

#[test]
fn ions_take_the_shapes_of_their_lone_pairs() {
    let cases = [
        ("B", 0, 3, 120.0),
        ("B", -1, 4, TETRAHEDRAL),
        ("N", 1, 4, TETRAHEDRAL),
        ("N", -1, 2, TETRAHEDRAL),
        ("O", 1, 3, TETRAHEDRAL),
        ("O", -1, 1, 0.0),
        ("F", -1, 0, 0.0),
        ("S", 1, 3, TETRAHEDRAL),
        ("P", 1, 4, TETRAHEDRAL),
    ];
    for (symbol, formal_charge, count, angle) in cases {
        let (mut editor, atom) = built(symbol);
        charge(&mut editor, atom, formal_charge);
        assert_settled(&editor, atom, count, &[angle]);
        editor.finish().unwrap();
    }
    let (editor, boron) = built("B");
    assert_planar(&editor, boron);
}

#[test]
fn oxygen_takes_only_builder_charges() {
    let (mut editor, oxygen) = built("O");
    assert_settled(&editor, oxygen, 2, &[TETRAHEDRAL]);
    charge(&mut editor, oxygen, 1);
    assert_settled(&editor, oxygen, 3, &[TETRAHEDRAL]);
    charge(&mut editor, oxygen, -1);
    assert_settled(&editor, oxygen, 1, &[]);
    charge(&mut editor, oxygen, 0);
    assert_settled(&editor, oxygen, 2, &[TETRAHEDRAL]);

    let before = editor.positions();
    let value = atom_of(&editor, oxygen).with_formal_charge(2);
    editor.replace_atom(oxygen, value).unwrap();
    let error = editor.settle_atoms([oxygen]).unwrap_err();
    assert!(matches!(
        &error,
        PlacementError::BuilderValence {
            atom,
            error: BuilderValenceError::UnsupportedCharge { .. },
        } if *atom == oxygen
    ));
    assert!(error
        .to_string()
        .ends_with("is O with charge +2. O takes charge -1, 0 or +1."));
    // Nothing changed.
    assert_eq!(editor.positions(), before);
    assert_eq!(hydrogens(&editor, oxygen).len(), 2);
}

/// Dimethyl sulfide, built from hydrogen sulfide.
fn dimethyl_sulfide() -> (StructureEditor, EditAtomId, [EditAtomId; 2]) {
    let (mut editor, sulfur) = built("S");
    assert_settled(&editor, sulfur, 2, &[TETRAHEDRAL]);
    let first = grow(&mut editor, sulfur, "C");
    let second = grow(&mut editor, sulfur, "C");
    assert_settled(&editor, sulfur, 0, &[TETRAHEDRAL]);
    (editor, sulfur, [first, second])
}

#[test]
fn sulfide_charged_and_back_is_a_sulfide() {
    let (mut editor, sulfur, carbons) = dimethyl_sulfide();
    for carbon in carbons {
        assert_settled(&editor, carbon, 3, &[TETRAHEDRAL]);
        assert_length(&editor, sulfur, carbon, BondOrder::Single);
    }
    charge(&mut editor, sulfur, 1);
    assert_settled(&editor, sulfur, 1, &[TETRAHEDRAL]);
    let value = atom_of(&editor, sulfur).with_formal_charge(2);
    editor.replace_atom(sulfur, value).unwrap();
    assert!(editor.settle_atoms([sulfur]).is_err());
    charge(&mut editor, sulfur, 0);
    assert_settled(&editor, sulfur, 0, &[TETRAHEDRAL]);
    let value = atom_of(&editor, sulfur).with_formal_charge(-1);
    editor.replace_atom(sulfur, value).unwrap();
    assert!(matches!(
        editor.settle_atoms([sulfur]),
        Err(PlacementError::BuilderValence {
            error: BuilderValenceError::ExcessValence { valence: 2, .. },
            ..
        })
    ));
}

#[test]
fn sulfoxide_and_sulfone_grow_from_a_sulfide() {
    let (mut editor, sulfur, carbons) = dimethyl_sulfide();
    for oxygens in 1..=2 {
        let directions = editor.free_directions(sulfur, BondOrder::Double).unwrap();
        assert_eq!(directions.len(), 1);
        assert!(
            editor
                .free_directions(sulfur, BondOrder::Single)
                .unwrap()
                .len()
                == 1
        );
        let attached = editor
            .attach_atom(
                sulfur,
                Atom::new(Element::O),
                BondOrder::Double,
                AttachmentSite::Toward(directions[0]),
            )
            .unwrap();
        assert!(attached.replaced_hydrogens().is_empty());
        let oxygen = attached.atom();
        editor.settle_atoms([sulfur, oxygen]).unwrap();
        assert_settled(&editor, sulfur, 0, &[TETRAHEDRAL]);
        assert_settled(&editor, oxygen, 0, &[]);
        assert_length(&editor, sulfur, oxygen, BondOrder::Double);
        assert_eq!(heavy_neighbors(&editor, sulfur).len(), 2 + oxygens);
        for carbon in carbons {
            assert_settled(&editor, carbon, 3, &[TETRAHEDRAL]);
        }
    }
    // A sulfone has no room for another oxygen.
    assert!(editor
        .free_directions(sulfur, BondOrder::Double)
        .unwrap()
        .is_empty());
    assert!(matches!(
        editor.attach_atom(
            sulfur,
            Atom::new(Element::O),
            BondOrder::Double,
            AttachmentSite::Toward(Vector3::new(1.0, 0.0, 0.0)),
        ),
        Err(PlacementError::NoFreePosition { atom }) if atom == sulfur
    ));
    editor.finish().unwrap();
}

#[test]
fn a_double_bond_in_place_of_a_hydrogen_takes_two() {
    let (mut editor, carbon) = built("C");
    let hydrogen = hydrogens(&editor, carbon)[0];
    let attached = editor
        .attach_atom(
            carbon,
            Atom::new(Element::O),
            BondOrder::Double,
            AttachmentSite::Hydrogen(hydrogen),
        )
        .unwrap();
    assert_eq!(attached.replaced_hydrogens().len(), 2);
    assert_eq!(attached.replaced_hydrogens()[0], hydrogen);
    let oxygen = attached.atom();
    editor.settle_atoms([carbon, oxygen]).unwrap();
    assert_settled(&editor, carbon, 2, &[120.0]);
    assert_planar(&editor, carbon);
    assert_settled(&editor, oxygen, 0, &[]);
    assert_length(&editor, carbon, oxygen, BondOrder::Double);

    // A nitrile from methane: three hydrogens go.
    let (mut editor, carbon) = built("C");
    let hydrogen = hydrogens(&editor, carbon)[0];
    let attached = editor
        .attach_atom(
            carbon,
            Atom::new(Element::N),
            BondOrder::Triple,
            AttachmentSite::Hydrogen(hydrogen),
        )
        .unwrap();
    assert_eq!(attached.replaced_hydrogens().len(), 3);
    let nitrogen = attached.atom();
    editor.settle_atoms([carbon, nitrogen]).unwrap();
    assert_settled(&editor, carbon, 1, &[180.0]);
    assert_settled(&editor, nitrogen, 0, &[]);
    editor.finish().unwrap();
}

#[test]
fn a_multiple_bond_past_the_valence_has_no_position() {
    let (mut editor, carbon) = built("C");
    let oxygen = grow(&mut editor, carbon, "O");
    let hydrogen = hydrogens(&editor, oxygen)[0];
    let before = editor.positions();
    let error = editor
        .attach_atom(
            oxygen,
            Atom::new(Element::C),
            BondOrder::Double,
            AttachmentSite::Hydrogen(hydrogen),
        )
        .unwrap_err();
    assert!(matches!(error, PlacementError::NoFreePosition { atom } if atom == oxygen));
    assert_eq!(editor.positions(), before);
}

#[test]
fn pyramids_keep_their_face_and_lone_atoms_their_frame() {
    let (mut editor, nitrogen) = built("N");
    let first = grow(&mut editor, nitrogen, "C");
    let second = grow(&mut editor, nitrogen, "C");
    let hydrogen = hydrogens(&editor, nitrogen)[0];
    // Mirror the hydrogen through the plane of the nitrogen and carbons.
    let center = point(&editor, nitrogen);
    let normal = (point(&editor, first) - center).cross(point(&editor, second) - center);
    let normal = normal / normal.norm();
    let offset = point(&editor, hydrogen) - center;
    let mirrored = center + (offset - normal * (2.0 * offset.dot(normal)));
    editor
        .set_position(hydrogen, Quantity::new(mirrored, ANGSTROM))
        .unwrap();
    let report = editor.settle_atoms([nitrogen]).unwrap();
    assert!(report.atoms()[0].moved_hydrogens().is_empty());
    assert!((point(&editor, hydrogen) - mirrored).norm() < 1e-6);

    let (mut editor, oxygen) = built("O");
    let hydrogen = hydrogens(&editor, oxygen)[0];
    let direction = point(&editor, hydrogen) - point(&editor, oxygen);
    let direction = direction / direction.norm();
    charge(&mut editor, oxygen, 1);
    let after = point(&editor, hydrogen) - point(&editor, oxygen);
    let after = after / after.norm();
    assert!((after - direction).norm() < 1e-9);
}

#[test]
fn ethane_to_ethyne_and_back() {
    let (mut editor, first) = built("C");
    let second = grow(&mut editor, first, "C");
    let steps = [
        (BondOrder::Single, 3, TETRAHEDRAL),
        (BondOrder::Double, 2, 120.0),
        (BondOrder::Triple, 1, 180.0),
        (BondOrder::Double, 2, 120.0),
        (BondOrder::Single, 3, TETRAHEDRAL),
    ];
    for (order, count, angle) in steps {
        reorder(&mut editor, first, second, order);
        assert_length(&editor, first, second, order);
        for carbon in [first, second] {
            assert_settled(&editor, carbon, count, &[angle]);
        }
        if order == BondOrder::Double {
            assert_planar(&editor, first);
            assert_planar(&editor, second);
            // The four hydrogens lie in one plane with the carbons.
            let axis = point(&editor, second) - point(&editor, first);
            let on_first = point(&editor, hydrogens(&editor, first)[0]) - point(&editor, first);
            let normal = axis.cross(on_first);
            for hydrogen in hydrogens(&editor, second) {
                let offset = point(&editor, hydrogen) - point(&editor, first);
                assert!(normal.dot(offset).abs() / normal.norm() < 0.01);
            }
        }
    }
    let structure = editor.finish().unwrap();
    assert_eq!(structure.atom_count(), 8);
}

#[test]
fn propyne_bends_and_straightens_at_its_middle_carbon() {
    let (mut editor, methyl) = built("C");
    let middle = grow(&mut editor, methyl, "C");
    let end = grow(&mut editor, middle, "C");
    let steps = [
        (BondOrder::Triple, 0, 1, 180.0),
        (BondOrder::Double, 1, 2, 120.0),
        (BondOrder::Single, 2, 3, TETRAHEDRAL),
        (BondOrder::Triple, 0, 1, 180.0),
        (BondOrder::Single, 2, 3, TETRAHEDRAL),
    ];
    for (order, middle_count, end_count, middle_angle) in steps {
        reorder(&mut editor, middle, end, order);
        assert_length(&editor, middle, end, order);
        assert_length(&editor, methyl, middle, BondOrder::Single);
        assert_settled(&editor, middle, middle_count, &[middle_angle]);
        assert!((angle(&editor, methyl, middle, end) - middle_angle).abs() < ANGLE);
        let end_angle = match order {
            BondOrder::Triple => 180.0,
            BondOrder::Double => 120.0,
            _ => TETRAHEDRAL,
        };
        assert_settled(&editor, end, end_count, &[end_angle]);
        assert_settled(&editor, methyl, 3, &[TETRAHEDRAL]);
    }
    // An internal triple bond straightens both of its atoms.
    let further = grow(&mut editor, end, "C");
    reorder(&mut editor, middle, end, BondOrder::Triple);
    assert!((angle(&editor, methyl, middle, end) - 180.0).abs() < ANGLE);
    assert!((angle(&editor, middle, end, further) - 180.0).abs() < ANGLE);
    assert_settled(&editor, further, 3, &[TETRAHEDRAL]);
    editor.finish().unwrap();
}

#[test]
fn element_changes_give_each_element_its_hydrogens() {
    let (mut editor, atom) = built("C");
    let steps = [
        ("N", 3, TETRAHEDRAL),
        ("O", 2, TETRAHEDRAL),
        ("F", 1, 0.0),
        ("S", 2, TETRAHEDRAL),
        ("P", 3, TETRAHEDRAL),
        ("B", 3, 120.0),
        ("Si", 4, TETRAHEDRAL),
        ("C", 4, TETRAHEDRAL),
    ];
    for (symbol, count, angle) in steps {
        change_element(&mut editor, atom, symbol);
        assert_settled(&editor, atom, count, &[angle]);
    }
    let structure = editor.finish().unwrap();
    assert_eq!(structure.atom_count(), 5);
}

#[test]
fn element_changes_in_a_carbonyl_settle_its_bond() {
    let (mut editor, carbon) = built("C");
    let oxygen = grow(&mut editor, carbon, "O");
    reorder(&mut editor, carbon, oxygen, BondOrder::Double);
    assert_settled(&editor, carbon, 2, &[120.0]);
    assert_settled(&editor, oxygen, 0, &[]);
    for symbol in ["S", "N"] {
        let value = atom_of(&editor, oxygen).with_element(element(symbol));
        editor.replace_atom(oxygen, value).unwrap();
        editor
            .settle_bonds([bond(&editor, carbon, oxygen)])
            .unwrap();
        editor.settle_atoms([carbon, oxygen]).unwrap();
        assert_length(&editor, carbon, oxygen, BondOrder::Double);
    }
    // An imine: the nitrogen takes one hydrogen, bent.
    assert_settled(&editor, oxygen, 1, &[120.0]);
    editor.finish().unwrap();
}

#[test]
fn a_cap_hydrogen_takes_the_angle_of_its_atom() {
    let (mut editor, sulfur, [first, second]) = dimethyl_sulfide();
    let mut methyl = hydrogens(&editor, second);
    methyl.push(second);
    editor.delete_atoms(methyl).unwrap();
    let report = editor.settle_atoms([sulfur]).unwrap();
    assert_eq!(report.added_count(), 1);
    assert_settled(&editor, sulfur, 1, &[TETRAHEDRAL]);
    assert_settled(&editor, first, 3, &[TETRAHEDRAL]);
    editor.finish().unwrap();
}

#[test]
fn opening_a_ring_settles_both_atoms() {
    let mut editor = StructureEditor::new();
    let placed = editor
        .place_fragment(
            fragment("benzene").unwrap(),
            Quantity::new(Point3::origin(), ANGSTROM),
        )
        .unwrap();
    let carbons = placed
        .atoms()
        .iter()
        .map(|(_, atom)| *atom)
        .filter(|atom| !is_hydrogen(&editor, *atom))
        .collect::<Vec<_>>();
    let (a, b) = carbons
        .iter()
        .flat_map(|a| carbons.iter().map(move |b| (*a, *b)))
        .find(|(a, b)| {
            editor
                .topology()
                .bond_between(*a, *b)
                .unwrap()
                .and_then(|id| editor.topology().bond(id))
                .is_some_and(|bond| bond.order == BondOrder::Single)
        })
        .unwrap();

    // A ring bond keeps its length and says so.
    let ring_bond = bond(&editor, a, b);
    editor.set_bond_order(ring_bond, BondOrder::Double).unwrap();
    let report = editor.settle_bonds([ring_bond]).unwrap();
    assert!(report.bonds()[0].in_ring());
    assert!(report.bonds()[0].moved_atoms().is_empty());
    editor.set_bond_order(ring_bond, BondOrder::Single).unwrap();

    editor.delete_bond(ring_bond).unwrap();
    editor.settle_atoms([a, b]).unwrap();
    for atom in [a, b] {
        assert_settled(&editor, atom, 2, &[120.0]);
        assert_planar(&editor, atom);
    }
    editor.finish().unwrap();
}

#[test]
fn bonds_settle_by_moving_their_smaller_side() {
    let (mut editor, first) = built("C");
    let second = grow(&mut editor, first, "C");
    let third = grow(&mut editor, second, "C");
    let before = point(&editor, first);
    let id = bond(&editor, first, second);
    editor.set_bond_order(id, BondOrder::Double).unwrap();
    let report = editor.settle_bonds([id]).unwrap();
    let settled = &report.bonds()[0];
    assert_eq!(settled.bond(), id);
    assert!(!settled.in_ring());
    // The methyl end and its hydrogens move; the ethyl side stays.
    assert_eq!(settled.moved_atoms().len(), 4);
    assert!(settled.moved_atoms().contains(&first));
    assert!(!settled.moved_atoms().contains(&third));
    assert!((point(&editor, first) - before).norm() > 0.1);
    assert_length(&editor, first, second, BondOrder::Double);
    // An ideal bond does not move.
    let report = editor.settle_bonds([id]).unwrap();
    assert!(report.bonds()[0].moved_atoms().is_empty());
}

#[test]
fn settling_rejects_states_outside_the_builder_and_changes_nothing() {
    let (mut editor, carbon) = built("C");
    let other = grow(&mut editor, carbon, "C");
    // Five units of bond valence on a carbon.
    let extra = editor
        .add_atom(
            Atom::new(Element::C),
            Quantity::new(Point3::new(0.0, 0.0, 3.0), ANGSTROM),
        )
        .unwrap();
    editor
        .add_bond(carbon, extra, BondOrder::Quadruple)
        .unwrap();
    let before = editor.positions();
    let count = editor.topology().atom_count();
    let error = editor.settle_atoms([other, carbon]).unwrap_err();
    assert!(matches!(
        &error,
        PlacementError::BuilderValence {
            atom,
            error: BuilderValenceError::ExcessValence { valence: 5, .. },
        } if *atom == carbon
    ));
    assert!(error
        .to_string()
        .ends_with("has valence 5. C takes valence 4."));
    assert_eq!(editor.positions(), before);
    assert_eq!(editor.topology().atom_count(), count);
}

#[test]
fn elements_outside_the_builder_follow_the_valence_model() {
    let (editor, germanium) = built("Ge");
    assert_eq!(hydrogens(&editor, germanium).len(), 4);
    let (editor, sodium) = built("Na");
    assert!(hydrogens(&editor, sodium).is_empty());
}

/// Attaches a single-bonded `symbol` to `parent` at its first free
/// position and settles both, as a builder does.
fn grow_free(editor: &mut StructureEditor, parent: EditAtomId, symbol: &str) -> EditAtomId {
    let direction = editor.free_directions(parent, BondOrder::Single).unwrap()[0];
    let atom = editor
        .attach_atom(
            parent,
            Atom::new(element(symbol)),
            BondOrder::Single,
            AttachmentSite::Toward(direction),
        )
        .unwrap()
        .atom();
    editor.settle_atoms([parent, atom]).unwrap();
    atom
}

/// The angles between the neighbors of `atom`, in degrees, ascending.
fn neighbor_angles(editor: &StructureEditor, atom: EditAtomId) -> Vec<f64> {
    let neighbors = editor
        .topology()
        .neighbors(atom)
        .unwrap()
        .collect::<Vec<_>>();
    let mut angles = neighbors
        .iter()
        .enumerate()
        .flat_map(|(index, a)| {
            neighbors[index + 1..]
                .iter()
                .map(move |b| angle(editor, *a, atom, *b))
        })
        .collect::<Vec<_>>();
    angles.sort_by(f64::total_cmp);
    angles
}

/// Asserts that the angles at `atom` are, within [`ANGLE`], `count` each of
/// 90°, 120°, and 180°.
fn assert_hypervalent(editor: &StructureEditor, atom: EditAtomId, counts: [usize; 3]) {
    let angles = neighbor_angles(editor, atom);
    let found = [90.0, 120.0, 180.0].map(|ideal| {
        angles
            .iter()
            .filter(|angle| (**angle - ideal).abs() < ANGLE)
            .count()
    });
    assert_eq!(found, counts, "{angles:?}");
}

/// Methanesulfenic acid, CH3–S–OH, built from hydrogen sulfide.
fn methanesulfenic_acid() -> (StructureEditor, EditAtomId) {
    let (mut editor, sulfur) = built("S");
    grow(&mut editor, sulfur, "C");
    grow(&mut editor, sulfur, "O");
    assert_settled(&editor, sulfur, 0, &[TETRAHEDRAL]);
    (editor, sulfur)
}

#[test]
fn a_third_neighbor_makes_a_seesaw_with_its_lone_pair_equatorial() {
    let (mut editor, sulfur) = methanesulfenic_acid();
    let chlorine = grow_free(&mut editor, sulfur, "Cl");
    // Sulfur(IV) takes a hydrogen: three heavy neighbors, a hydrogen, and an
    // equatorial lone pair, so one axial pair, one equatorial pair, and four
    // right angles.
    assert_settled(&editor, sulfur, 1, &[90.0, 120.0, 180.0]);
    assert_hypervalent(&editor, sulfur, [4, 1, 1]);
    assert!(neighbor_angles(&editor, sulfur)[0] > 85.0);
    assert_length(&editor, sulfur, chlorine, BondOrder::Single);
    editor.finish().unwrap();
}

#[test]
fn a_fourth_neighbor_makes_phosphorus_a_trigonal_bipyramid() {
    let (mut editor, phosphorus) = built("P");
    let carbons = [(); 3].map(|()| grow(&mut editor, phosphorus, "C"));
    assert_settled(&editor, phosphorus, 0, &[TETRAHEDRAL]);
    grow_free(&mut editor, phosphorus, "Cl");
    assert_settled(&editor, phosphorus, 1, &[90.0, 120.0, 180.0]);
    assert_hypervalent(&editor, phosphorus, [6, 3, 1]);
    for carbon in carbons {
        assert_settled(&editor, carbon, 3, &[TETRAHEDRAL]);
    }
    editor.finish().unwrap();
}

#[test]
fn sulfur_four_to_six_settles_octahedral() {
    let (mut editor, sulfur) = methanesulfenic_acid();
    grow_free(&mut editor, sulfur, "Cl");
    // In place of the hydrogen, a fourth heavy neighbor: sulfur(IV) without
    // hydrogens.
    grow(&mut editor, sulfur, "F");
    assert_settled(&editor, sulfur, 0, &[90.0, 120.0, 180.0]);
    assert_hypervalent(&editor, sulfur, [4, 1, 1]);
    // A fifth makes sulfur(VI), which takes a hydrogen.
    let report = {
        let direction = editor.free_directions(sulfur, BondOrder::Single).unwrap()[0];
        let fluorine = editor
            .attach_atom(
                sulfur,
                Atom::new(Element::F),
                BondOrder::Single,
                AttachmentSite::Toward(direction),
            )
            .unwrap()
            .atom();
        editor.settle_atoms([sulfur, fluorine]).unwrap()
    };
    assert_eq!(report.added_count(), 1);
    assert!(!report.atoms()[0].turned_atoms().is_empty());
    assert_settled(&editor, sulfur, 1, &[90.0, 180.0]);
    assert_hypervalent(&editor, sulfur, [12, 0, 3]);
    editor.finish().unwrap();
}

#[test]
fn a_ring_holds_a_hypervalent_angle_and_says_so() {
    let mut editor = StructureEditor::new();
    let placed = editor
        .place_fragment(
            fragment("cyclopentane").unwrap(),
            Quantity::new(Point3::origin(), ANGSTROM),
        )
        .unwrap();
    let atom = placed
        .atoms()
        .iter()
        .map(|(_, atom)| *atom)
        .find(|atom| !is_hydrogen(&editor, *atom))
        .unwrap();
    change_element(&mut editor, atom, "S");
    let ring = heavy_neighbors(&editor, atom);
    let positions = ring
        .iter()
        .map(|neighbor| point(&editor, *neighbor))
        .collect::<Vec<_>>();
    let direction = editor.free_directions(atom, BondOrder::Single).unwrap()[0];
    let chlorine = editor
        .attach_atom(
            atom,
            Atom::new(element("Cl")),
            BondOrder::Single,
            AttachmentSite::Toward(direction),
        )
        .unwrap()
        .atom();
    let report = editor.settle_atoms([atom, chlorine]).unwrap();
    // The ring's angle at sulfur, near the pentagon's 108°, is far from the
    // seesaw's equatorial 120° and axial 90°; the ring stays, and the
    // chlorine turns.
    let settled = &report.atoms()[0];
    assert!(settled.angle_in_ring());
    assert_eq!(settled.turned_atoms(), [chlorine]);
    assert_eq!(settled.added_hydrogens().len(), 1);
    for (neighbor, before) in ring.iter().zip(&positions) {
        assert!((point(&editor, *neighbor) - *before).norm() < 1e-9);
    }
    editor.finish().unwrap();
}

/// Seeded xorshift.
struct Random(u64);

impl Random {
    fn below(&mut self, bound: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound as u64) as usize
    }
}

/// One random builder edit of `editor` and its settling; `false` when the
/// edit was refused or leaves a state outside the builder, which the caller
/// undoes.
fn random_edit(editor: &mut StructureEditor, random: &mut Random) -> bool {
    const ELEMENTS: [&str; 8] = ["C", "N", "O", "S", "P", "Se", "Cl", "F"];
    let heavy = editor
        .topology()
        .atom_ids()
        .filter(|atom| !is_hydrogen(editor, *atom))
        .collect::<Vec<_>>();
    // Half the edits are of the first atom, so that it grows hypervalent,
    // and half of all edits attach an atom.
    let atom = heavy[random.below(2 * heavy.len()).saturating_sub(heavy.len())];
    let settled = match random.below(6) {
        0..=2 => {
            let order = [BondOrder::Single, BondOrder::Single, BondOrder::Double][random.below(3)];
            let own = hydrogens(editor, atom);
            let site = if !own.is_empty() && random.below(2) == 0 {
                AttachmentSite::Hydrogen(own[random.below(own.len())])
            } else {
                match editor.free_directions(atom, order) {
                    Ok(directions) if !directions.is_empty() => {
                        AttachmentSite::Toward(directions[random.below(directions.len())])
                    }
                    _ => return false,
                }
            };
            let symbol = ELEMENTS[random.below(ELEMENTS.len())];
            let Ok(attached) = editor.attach_atom(atom, Atom::new(element(symbol)), order, site)
            else {
                return false;
            };
            editor.settle_atoms([atom, attached.atom()])
        }
        3 => {
            let symbol = ELEMENTS[random.below(ELEMENTS.len())];
            let value = atom_of(editor, atom).with_element(element(symbol));
            editor.replace_atom(atom, value).unwrap();
            let bonds = heavy_neighbors(editor, atom)
                .into_iter()
                .map(|neighbor| bond(editor, atom, neighbor))
                .collect::<Vec<_>>();
            editor.settle_bonds(bonds).unwrap();
            editor.settle_atoms([atom])
        }
        4 => {
            let value = atom_of(editor, atom);
            let charges = builder_charges(value.element()).unwrap();
            let charge = charges[random.below(charges.len())];
            editor
                .replace_atom(atom, value.with_formal_charge(charge))
                .unwrap();
            editor.settle_atoms([atom])
        }
        _ => {
            let neighbors = heavy_neighbors(editor, atom);
            if neighbors.is_empty() {
                return false;
            }
            let other = neighbors[random.below(neighbors.len())];
            let order = [BondOrder::Single, BondOrder::Double, BondOrder::Triple][random.below(3)];
            let id = bond(editor, atom, other);
            editor.set_bond_order(id, order).unwrap();
            editor.settle_bonds([id]).unwrap();
            editor.settle_atoms([atom, other])
        }
    };
    match settled {
        Ok(_) => true,
        Err(
            PlacementError::BuilderValence { .. }
            | PlacementError::InvalidValence { .. }
            | PlacementError::UnsupportedGeometry { .. },
        ) => false,
        Err(error) => panic!("settling failed: {error}"),
    }
}

#[test]
fn random_edits_keep_every_angle_of_a_tree_open() {
    let mut random = Random(0x2545_f491_4f6c_dd1d);
    for sequence in 0..300 {
        let symbol = ["S", "P", "Se", "N", "C"][sequence % 5];
        let (mut editor, _) = built(symbol);
        for step in 0..12 {
            let before = editor.clone();
            if !random_edit(&mut editor, &mut random) {
                editor = before;
                continue;
            }
            for atom in editor.topology().atom_ids() {
                let angles = neighbor_angles(&editor, atom);
                assert!(
                    angles.first().is_none_or(|smallest| *smallest >= 85.0),
                    "sequence {sequence}, step {step}: {} has {angles:?}",
                    atom_of(&editor, atom).element().symbol()
                );
            }
        }
        editor.finish().unwrap();
    }
}
