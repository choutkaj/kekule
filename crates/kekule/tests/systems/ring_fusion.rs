//! Fusing ring fragments onto bonds of structure drafts.
use kekule::math::{Point3, Vector3};
use kekule::molecule::{BondOrder, ImplicitHydrogens, Molecule};
use kekule::smiles::{self, SmilesWriteOptions};
use kekule::structure::placement::{
    fragment, fragments, ideal_bond_length, FragmentKind, PlacementError,
};
use kekule::structure::{Positions, Structure, StructureEditor};
use kekule::topology::edit::{EditAtomId, EditBondId, PrunedStereo};
use kekule::units::{Quantity, ANGSTROM, NANOMETER};
use std::collections::{BTreeMap, BTreeSet};

/// A draft holding the library fragment `name` alone, with its hydrogens.
fn placed(name: &str) -> StructureEditor {
    let mut editor = StructureEditor::new();
    editor
        .place_fragment(
            fragment(name).unwrap(),
            Quantity::new(Point3::origin(), NANOMETER),
        )
        .unwrap();
    editor
}

fn is_heavy(editor: &StructureEditor, atom: EditAtomId) -> bool {
    editor
        .topology()
        .atom(atom)
        .unwrap()
        .element()
        .atomic_number()
        != 1
}

/// The draft's bonds between two heavy atoms, in handle order.
fn heavy_bonds(editor: &StructureEditor) -> Vec<EditBondId> {
    editor
        .topology()
        .bonds()
        .filter(|(_, bond)| is_heavy(editor, bond.a()) && is_heavy(editor, bond.b()))
        .map(|(id, _)| id)
        .collect()
}

/// The first bond between two heavy atoms.
fn heavy_bond(editor: &StructureEditor) -> EditBondId {
    heavy_bonds(editor)[0]
}

fn point(editor: &StructureEditor, atom: EditAtomId) -> Point3 {
    editor
        .position(atom)
        .unwrap()
        .into_unit(ANGSTROM)
        .unwrap()
        .into_value()
}

/// The torsion of four points, in degrees.
fn torsion([a, b, c, d]: [Point3; 4]) -> f64 {
    let (b1, b2, b3) = (b - a, c - b, d - c);
    let (n1, n2) = (b1.cross(b2), b2.cross(b3));
    n1.cross(b2 / b2.norm())
        .dot(n2)
        .atan2(n1.dot(n2))
        .to_degrees()
}

/// The position of the atom at dense `index`, in Ă….
fn angstroms(structure: &Structure, index: usize) -> Point3 {
    structure
        .positions()
        .get(index)
        .unwrap()
        .into_unit(ANGSTROM)
        .unwrap()
        .into_value()
}

/// The canonical SMILES of the one molecule of `structure`, with its
/// hydrogens implicit.
fn canonical(structure: &Structure) -> String {
    assert_eq!(structure.topology().instance_count(), 1);
    let mut molecule = perceived_molecule(structure);
    molecule.remove_hydrogens().unwrap();
    molecule.perceive().unwrap();
    smiles::write(&molecule, SmilesWriteOptions::canonical()).unwrap()
}

fn canonical_of(text: &str) -> String {
    let mut molecule = smiles::to_molecules(text).unwrap().remove(0);
    molecule.perceive().unwrap();
    smiles::write(&molecule, SmilesWriteOptions::canonical()).unwrap()
}

fn perceived_molecule(structure: &Structure) -> Molecule {
    let mut structure = structure.clone();
    structure.perceive().unwrap();
    let topology = structure.topology();
    let molecule = topology.definitions().next().unwrap().molecule().clone();
    molecule
}

/// The sizes of the perceived rings and the number of aromatic atoms.
fn rings_and_aromatic_atoms(structure: &Structure) -> (Vec<usize>, usize) {
    let molecule = perceived_molecule(structure);
    let mut sizes = molecule
        .perception()
        .rings()
        .unwrap()
        .ring_set()
        .unwrap()
        .rings()
        .iter()
        .map(|ring| ring.atoms.len())
        .collect::<Vec<_>>();
    sizes.sort_unstable();
    let aromaticity = molecule.perception().aromaticity().unwrap();
    let aromatic = molecule
        .atom_ids()
        .filter(|atom| aromaticity.atom_is_aromatic(*atom))
        .count();
    (sizes, aromatic)
}

/// Bond lengths lie within 85 % to 112 % of an ideal single bond's, and atoms
/// three or more bonds apart keep 2.3 Ă… from each other between heavy atoms and
/// 1.6 Ă… otherwise; unrelaxed bay hydrogens, as in phenanthrene, come to 1.75 Ă….
fn assert_sensible_geometry(structure: &Structure, label: &str) {
    let topology = structure.topology();
    let points = (0..structure.atom_count())
        .map(|index| angstroms(structure, index))
        .collect::<Vec<_>>();
    let mut neighbors = BTreeMap::<usize, BTreeSet<usize>>::new();
    for bond in topology.bonds() {
        let [a, b] = bond.atoms().map(|atom| atom.index().index());
        let length = (points[a] - points[b]).norm();
        let [first, second] = bond.atoms().map(|atom| atom.element());
        let single = ideal_bond_length(first, second, BondOrder::Single)
            .unwrap()
            .into_unit(ANGSTROM)
            .unwrap()
            .into_value();
        let ratio = length / single;
        assert!(
            (0.85..1.12).contains(&ratio),
            "{label}: bond {a}-{b} is {length:.3} A"
        );
        neighbors.entry(a).or_default().insert(b);
        neighbors.entry(b).or_default().insert(a);
    }
    let heavy = topology
        .atoms()
        .map(|atom| atom.element().atomic_number() != 1)
        .collect::<Vec<_>>();
    let near = |atom: usize| {
        let mut seen = BTreeSet::from([atom]);
        for _ in 0..2 {
            let next = seen
                .iter()
                .flat_map(|atom| neighbors.get(atom).into_iter().flatten().copied())
                .collect::<Vec<_>>();
            seen.extend(next);
        }
        seen
    };
    for a in 0..points.len() {
        let close = near(a);
        for b in a + 1..points.len() {
            if !close.contains(&b) {
                let found = (points[a] - points[b]).norm();
                let limit = if heavy[a] && heavy[b] { 2.3 } else { 1.6 };
                assert!(
                    found > limit,
                    "{label}: atoms {a} and {b} are {found:.2} A apart"
                );
            }
        }
    }
}

/// The largest distance of a heavy atom from the plane through the first
/// three heavy atoms, in Ă….
fn heavy_atoms_off_plane(structure: &Structure) -> f64 {
    let points = structure
        .topology()
        .atoms()
        .filter(|atom| atom.element().atomic_number() != 1)
        .map(|atom| angstroms(structure, atom.index().index()))
        .collect::<Vec<_>>();
    let normal = (points[1] - points[0]).cross(points[2] - points[0]);
    let normal = normal / normal.norm();
    points
        .iter()
        .map(|point| (*point - points[0]).dot(normal).abs())
        .fold(0.0, f64::max)
}

/// A draft's atoms (element, hydrogen declaration, position, name), bonds
/// (dense atoms and order), and represented stereo once published, with the
/// stereo the draft pruned.
type Snapshot = (
    Vec<(u8, String, Point3, Option<String>)>,
    Vec<(usize, usize, String)>,
    Vec<String>,
    Vec<PrunedStereo>,
);

fn snapshot(editor: StructureEditor) -> Snapshot {
    let pruned = editor.topology().pruned_stereo().to_vec();
    let structure = editor.finish().unwrap();
    let topology = structure.topology();
    let atoms = topology
        .atoms()
        .map(|atom| {
            (
                atom.element().atomic_number(),
                format!("{:?}", atom.atom().hydrogens()),
                angstroms(&structure, atom.index().index()),
                atom.atom_site()
                    .map(|site| site.metadata().label_atom_id().unwrap().to_owned()),
            )
        })
        .collect();
    let bonds = topology
        .bonds()
        .map(|bond| {
            let [a, b] = bond.atoms().map(|atom| atom.index().index());
            (a, b, format!("{:?}", bond.order()))
        })
        .collect();
    let stereo = topology
        .definitions()
        .flat_map(|definition| {
            definition
                .molecule()
                .stereo_elements()
                .map(|(_, element)| format!("{element:?}"))
                .collect::<Vec<_>>()
        })
        .collect();
    (atoms, bonds, stereo, pruned)
}

/// A draft of the molecule `text` with its atoms at `points`, in nm, in
/// SMILES order.
fn drafted(text: &str, points: &[Point3]) -> StructureEditor {
    let topology = smiles::to_topology(text).unwrap();
    assert_eq!(topology.atom_count(), points.len());
    Structure::new(
        topology,
        Positions::new(Quantity::new(points.to_vec(), NANOMETER)).unwrap(),
    )
    .unwrap()
    .into_editor()
}

/// The point bonded to `c` at `length` with the angle `b`â€“`c`â€“new and the
/// torsion `a`â€“`b`â€“`c`â€“new, in degrees.
fn bonded_point(a: Point3, b: Point3, c: Point3, length: f64, angle: f64, torsion: f64) -> Point3 {
    let bc = (c - b) / (c - b).norm();
    let normal = (b - a).cross(bc);
    let normal = normal / normal.norm();
    let across = normal.cross(bc);
    let (angle, torsion) = (angle.to_radians(), torsion.to_radians());
    c + bc * (-length * angle.cos())
        + across * (length * angle.sin() * torsion.cos())
        + normal * (length * angle.sin() * torsion.sin())
}

/// The ethylene of 2-butene in the xy plane, its methyls on one side
/// (`cis`) or opposite, in SMILES atom order.
fn butene(cis: bool) -> StructureEditor {
    let (text, y) = if cis {
        (r"C/C=C\C", 0.13)
    } else {
        ("C/C=C/C", -0.13)
    };
    drafted(
        text,
        &[
            Point3::new(-0.075, 0.13, 0.0),
            Point3::origin(),
            Point3::new(0.134, 0.0, 0.0),
            Point3::new(0.209, y, 0.0),
        ],
    )
}

/// The draft's double bond.
fn double_bond(editor: &StructureEditor) -> EditBondId {
    editor
        .topology()
        .bonds()
        .find(|(_, bond)| bond.order == BondOrder::Double)
        .unwrap()
        .0
}

/// The bond between the atoms named `first` and `second`.
fn named_bond(editor: &StructureEditor, first: &str, second: &str) -> EditBondId {
    let topology = editor.topology();
    let name = |atom: EditAtomId| {
        let site = topology.atom_site_for_atom(atom).unwrap().unwrap();
        topology
            .atom_site(site)
            .unwrap()
            .metadata()
            .label_atom_id()
            .unwrap()
            .to_owned()
    };
    topology
        .bonds()
        .find(|(_, bond)| {
            let names = [name(bond.a()), name(bond.b())];
            names == [first, second] || names == [second, first]
        })
        .unwrap()
        .0
}

#[test]
fn benzene_fused_on_any_benzene_bond_gives_naphthalene() {
    let benzene = fragment("benzene").unwrap();
    let source = placed("benzene");
    let bonds = heavy_bonds(&source);
    assert_eq!(bonds.len(), 6);
    let orders = bonds
        .iter()
        .map(|bond| format!("{:?}", source.topology().bond(*bond).unwrap().order))
        .collect::<BTreeSet<_>>();
    // Both KekulĂ© orders of the target bond are covered.
    assert_eq!(orders.len(), 2);
    for bond in bonds {
        let mut editor = source.clone();
        let fused = editor.fuse_fragment(bond, benzene).unwrap();
        // Four carbons and their four hydrogens are new; each target carbon
        // gives up its hydrogen.
        assert_eq!(fused.atoms().len(), 8);
        assert_eq!(fused.removed_hydrogens().len(), 2);
        assert_eq!(fused.bond(), bond);
        let target = editor.topology().bond(bond).unwrap();
        let merged = fused.fused_atoms().map(|(_, atom)| atom);
        assert_eq!(
            BTreeSet::from(merged),
            BTreeSet::from([target.a(), target.b()])
        );
        for hydrogen in fused.removed_hydrogens() {
            assert!(editor.topology().atom(*hydrogen).is_none());
        }
        let structure = editor.finish().unwrap();
        assert_eq!(structure.atom_count(), 18);
        assert_eq!(canonical(&structure), canonical_of("c1ccc2ccccc2c1"));
        assert_eq!(rings_and_aromatic_atoms(&structure), (vec![6, 6], 10));
        assert_sensible_geometry(&structure, "naphthalene");
        // In the plane of the first ring, on the side away from it.
        assert!(heavy_atoms_off_plane(&structure) < 0.05);
    }
}

#[test]
fn cyclopentadiene_fused_on_any_benzene_bond_gives_indene() {
    let cyclopentadiene = fragment("cyclopentadiene").unwrap();
    // The fusion bond is the C1=C2 double bond, next to the CH2.
    let fusion = cyclopentadiene
        .molecule()
        .bond(cyclopentadiene.fusion_bond().unwrap())
        .unwrap();
    assert_eq!(fusion.order(), BondOrder::Double);
    let source = placed("benzene");
    for bond in heavy_bonds(&source) {
        let mut editor = source.clone();
        let fused = editor.fuse_fragment(bond, cyclopentadiene).unwrap();
        // Three carbons and their four hydrogens are new; each target carbon
        // gives up its hydrogen.
        assert_eq!(fused.atoms().len(), 7);
        assert_eq!(fused.removed_hydrogens().len(), 2);
        let structure = editor.finish().unwrap();
        assert_eq!(structure.atom_count(), 17);
        assert_eq!(canonical(&structure), canonical_of("C1C=Cc2ccccc21"));
        assert_eq!(rings_and_aromatic_atoms(&structure), (vec![5, 6], 6));
        // C9H8: one saturated CH2, six CH, and the two fused carbons bare.
        let molecule = perceived_molecule(&structure);
        let hydrogens = |atom| {
            molecule
                .neighbors(atom)
                .unwrap()
                .filter(|n| molecule.atom(*n).unwrap().element().atomic_number() == 1)
                .count()
        };
        let carbons = molecule
            .atoms()
            .filter(|(_, atom)| atom.element().atomic_number() == 6)
            .map(|(id, _)| id)
            .collect::<Vec<_>>();
        let mut counts = carbons.iter().map(|c| hydrogens(*c)).collect::<Vec<_>>();
        counts.sort_unstable();
        assert_eq!(counts, [0, 0, 1, 1, 1, 1, 1, 1, 2]);
        let methylene = carbons.iter().find(|c| hydrogens(**c) == 2).unwrap();
        assert!(molecule
            .incident_bonds(*methylene)
            .unwrap()
            .all(|(_, bond)| bond.order() == BondOrder::Single));
        assert_sensible_geometry(&structure, "indene");
        assert!(heavy_atoms_off_plane(&structure) < 0.05);
    }
}

#[test]
fn cyclohexane_fused_on_cyclohexane_gives_trans_decalin() {
    let mut editor = placed("cyclohexane");
    let bond = heavy_bond(&editor);
    let fused = editor
        .fuse_fragment(bond, fragment("cyclohexane").unwrap())
        .unwrap();
    // Four carbons with eight hydrogens; each target carbon keeps one.
    assert_eq!(fused.atoms().len(), 12);
    assert_eq!(fused.removed_hydrogens().len(), 2);
    // Both rings take equatorial bonds of each other, so the bridgehead
    // hydrogens are anti.
    let [first, second] = fused.fused_atoms().map(|(_, atom)| atom);
    let hydrogen_on = |atom: EditAtomId| {
        editor
            .topology()
            .neighbors(atom)
            .unwrap()
            .find(|n| !is_heavy(&editor, *n))
            .unwrap()
    };
    let torsion = torsion(
        [hydrogen_on(first), first, second, hydrogen_on(second)].map(|atom| point(&editor, atom)),
    );
    assert!(torsion.abs() > 150.0, "{torsion:.1}");
    let structure = editor.finish().unwrap();
    assert_eq!(structure.atom_count(), 28);
    assert_eq!(canonical(&structure), canonical_of("C1CCC2CCCCC2C1"));
    assert_eq!(rings_and_aromatic_atoms(&structure), (vec![6, 6], 0));
    assert_sensible_geometry(&structure, "decalin");
}

#[test]
fn five_membered_rings_fuse_on_benzene() {
    let cases = [
        ("cyclopentane", "C1Cc2ccccc2C1", vec![5, 6], 6),
        ("cyclopentadiene", "C1C=Cc2ccccc21", vec![5, 6], 6),
        ("pyrrole", "c1ccc2[nH]ccc2c1", vec![5, 6], 9),
        ("thiophene", "c1ccc2sccc2c1", vec![5, 6], 9),
    ];
    for (name, expected, rings, aromatic) in cases {
        let mut editor = placed("benzene");
        let bond = heavy_bond(&editor);
        editor.fuse_fragment(bond, fragment(name).unwrap()).unwrap();
        let structure = editor.finish().unwrap();
        assert_eq!(canonical(&structure), canonical_of(expected), "{name}");
        assert_eq!(
            rings_and_aromatic_atoms(&structure),
            (rings, aromatic),
            "{name}"
        );
        assert_sensible_geometry(&structure, name);
    }
}

#[test]
fn every_ring_fragment_fuses_on_benzene() {
    let expected = [
        ("benzene", "c1ccc2ccccc2c1"),
        ("pyridine", "c1ccc2ncccc2c1"),
        ("pyrimidine", "c1ccc2ncncc2c1"),
        ("pyrazine", "c1ccc2nccnc2c1"),
        ("pyrrole", "c1ccc2[nH]ccc2c1"),
        ("furan", "c1ccc2occc2c1"),
        ("thiophene", "c1ccc2sccc2c1"),
        ("imidazole", "c1ccc2[nH]cnc2c1"),
        ("pyrazole", "c1ccc2[nH]ncc2c1"),
        ("oxazole", "c1ccc2ocnc2c1"),
        ("thiazole", "c1ccc2scnc2c1"),
        ("cyclopropane", "C1c2ccccc21"),
        ("cyclobutane", "C1Cc2ccccc21"),
        ("cyclopentane", "C1Cc2ccccc2C1"),
        ("cyclohexane", "C1CCc2ccccc2C1"),
        ("cyclopentadiene", "C1C=Cc2ccccc21"),
        ("piperidine", "C1CNc2ccccc2C1"),
        ("piperazine", "C1CNc2ccccc2N1"),
        ("morpholine", "C1COc2ccccc2N1"),
        ("tetrahydrofuran", "C1COc2ccccc21"),
        ("naphthalene", "c1ccc2c(c1)ccc1ccccc12"),
        ("indole", "c1ccc2c(c1)[nH]c1ccccc12"),
        ("quinoline", "c1ccc2cc3ccccc3nc2c1"),
        ("benzimidazole", "c1ccc2c(c1)ccc1[nH]cnc12"),
    ];
    let expected = expected.into_iter().collect::<BTreeMap<_, _>>();
    let rings = fragments()
        .iter()
        .filter(|fragment| fragment.kind() == FragmentKind::Ring)
        .collect::<Vec<_>>();
    assert_eq!(rings.len(), expected.len());
    for ring in rings {
        let source = placed("benzene");
        for bond in heavy_bonds(&source) {
            let mut editor = source.clone();
            let label = ring.name();
            editor
                .fuse_fragment(bond, ring)
                .unwrap_or_else(|error| panic!("{label}: {error}"));
            let structure = editor.finish().unwrap();
            assert_eq!(
                canonical(&structure),
                canonical_of(expected[label]),
                "{label}"
            );
            assert_sensible_geometry(&structure, label);
        }
    }
}

#[test]
fn benzene_fused_on_cyclohexane_gives_tetralin() {
    let mut editor = placed("cyclohexane");
    let bond = heavy_bond(&editor);
    let fused = editor
        .fuse_fragment(bond, fragment("benzene").unwrap())
        .unwrap();
    // Each target carbon gives up both its hydrogens to the aromatic ring.
    assert_eq!(fused.removed_hydrogens().len(), 4);
    let structure = editor.finish().unwrap();
    assert_eq!(canonical(&structure), canonical_of("c1ccc2c(c1)CCCC2"));
    assert_eq!(rings_and_aromatic_atoms(&structure), (vec![6, 6], 6));
    assert_sensible_geometry(&structure, "tetralin");
}

#[test]
fn implicit_hydrogens_follow_the_fusion() {
    for fixed in [false, true] {
        // Benzene with implicit hydrogens only, inferred or fixed.
        let mut editor = placed("benzene");
        let heavy = editor
            .topology()
            .atom_ids()
            .filter(|atom| is_heavy(&editor, *atom))
            .collect::<Vec<_>>();
        editor.remove_hydrogens(heavy.clone()).unwrap();
        assert_eq!(editor.topology().atom_count(), 6);
        if fixed {
            for atom in &heavy {
                let carbon = editor
                    .topology()
                    .atom(*atom)
                    .unwrap()
                    .clone()
                    .with_hydrogens(ImplicitHydrogens::Fixed(1));
                editor.replace_atom(*atom, carbon).unwrap();
            }
        }
        let bond = heavy_bond(&editor);
        let fused = editor
            .fuse_fragment(bond, fragment("benzene").unwrap())
            .unwrap();
        assert!(fused.removed_hydrogens().is_empty());
        if fixed {
            for (_, atom) in fused.fused_atoms() {
                assert_eq!(
                    editor.topology().atom(atom).unwrap().hydrogens(),
                    ImplicitHydrogens::Fixed(0)
                );
            }
        }
        let structure = editor.finish().unwrap();
        // Fixed counts write as bracket atoms, so compare inferred ones.
        if !fixed {
            assert_eq!(canonical(&structure), canonical_of("c1ccc2ccccc2c1"));
        }
        assert_eq!(rings_and_aromatic_atoms(&structure), (vec![6, 6], 10));
        let molecule = perceived_molecule(&structure);
        let hydrogens = molecule
            .atoms()
            .map(|(id, atom)| {
                usize::from(atom.element().atomic_number() == 1)
                    + molecule.implicit_hydrogens(id).unwrap()
            })
            .sum::<usize>();
        assert_eq!(hydrogens, 8, "fixed: {fixed}");
    }
}

#[test]
fn fused_atoms_join_the_target_residue_under_free_names() {
    let mut editor = placed("benzene");
    let bond = heavy_bond(&editor);
    let fused = editor
        .fuse_fragment(bond, fragment("pyridine").unwrap())
        .unwrap();
    let topology = editor.topology();
    let residue_of = |atom: EditAtomId| {
        let site = topology.atom_site_for_atom(atom).unwrap().unwrap();
        topology.atom_site(site).unwrap().residue()
    };
    let residue = residue_of(editor.topology().bond(bond).unwrap().a());
    for (_, atom) in fused.atoms() {
        assert_eq!(residue_of(*atom), residue);
    }
    let structure = editor.finish().unwrap();
    let residue = structure.topology().residues().next().unwrap();
    let names = residue
        .atom_sites()
        .map(|site| site.metadata().label_atom_id().unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(names.len(), structure.atom_count());
    let unique = names.iter().collect::<BTreeSet<_>>();
    assert_eq!(unique.len(), names.len(), "{names:?}");
    for name in ["N1", "C7", "C8", "C9"] {
        assert!(names.iter().any(|found| found == name), "{names:?}");
    }
}

#[test]
fn fusion_refusals_change_nothing() {
    let refused = |editor: &StructureEditor, bond: EditBondId, name: &str| {
        let mut trial = editor.clone();
        let error = trial
            .fuse_fragment(bond, fragment(name).unwrap())
            .unwrap_err();
        assert_eq!(snapshot(trial), snapshot(editor.clone()));
        error
    };

    let benzene = placed("benzene");
    let error = refused(&benzene, heavy_bond(&benzene), "methyl");
    assert!(matches!(error, PlacementError::NotRingFragment { .. }));
    assert_eq!(
        error.to_string(),
        "Fragment \"methyl\" is a group, not a ring."
    );

    let to_hydrogen = benzene
        .topology()
        .bonds()
        .find(|(_, bond)| !is_heavy(&benzene, bond.b()) || !is_heavy(&benzene, bond.a()))
        .unwrap()
        .0;
    let error = refused(&benzene, to_hydrogen, "benzene");
    assert_eq!(
        error,
        PlacementError::FusionAtHydrogen { bond: to_hydrogen }
    );
    assert_eq!(
        error.to_string(),
        format!("Cannot fuse a ring onto {to_hydrogen}. It ends at a hydrogen.")
    );

    let acetylene = placed("ethynyl");
    let triple = heavy_bond(&acetylene);
    assert_eq!(
        acetylene.topology().bond(triple).unwrap().order,
        BondOrder::Triple
    );
    let error = refused(&acetylene, triple, "benzene");
    assert_eq!(
        error.to_string(),
        format!("Cannot fuse a ring onto {triple}. It is a triple bond.")
    );

    // The bond naphthalene's rings share.
    let naphthalene = placed("naphthalene");
    let shared = heavy_bonds(&naphthalene)
        .into_iter()
        .find(|bond| {
            let bond = naphthalene.topology().bond(*bond).unwrap();
            [bond.a(), bond.b()].iter().all(|atom| {
                naphthalene
                    .topology()
                    .neighbors(*atom)
                    .unwrap()
                    .filter(|n| is_heavy(&naphthalene, *n))
                    .count()
                    == 3
            })
        })
        .unwrap();
    let error = refused(&naphthalene, shared, "benzene");
    assert_eq!(error, PlacementError::SharedFusionBond { bond: shared });
    assert_eq!(
        error.to_string(),
        format!("Cannot fuse a ring onto {shared}. It is already shared by two rings.")
    );

    // (E)-2-butene: the ring would close across a trans double bond.
    let trans = butene(false);
    let double = double_bond(&trans);
    let error = refused(&trans, double, "benzene");
    assert_eq!(error, PlacementError::TransFusionBond { bond: double });
    assert_eq!(
        error.to_string(),
        format!("Cannot fuse a ring onto {double}. Its stereo puts the ring bonds trans.")
    );

    // Toluene's ring bond at the methyl carbon: the ring carbon there carries
    // no hydrogen.
    let mut toluene = placed("benzene");
    let ring_bond = heavy_bond(&toluene);
    let ipso = toluene.topology().bond(ring_bond).unwrap().a();
    let hydrogen = toluene
        .topology()
        .neighbors(ipso)
        .unwrap()
        .find(|n| !is_heavy(&toluene, *n))
        .unwrap();
    toluene
        .attach_fragment(
            ipso,
            fragment("methyl").unwrap(),
            kekule::structure::placement::AttachmentSite::Hydrogen(hydrogen),
        )
        .unwrap();
    let error = refused(&toluene, ring_bond, "benzene");
    assert_eq!(
        error,
        PlacementError::FusionAtomWithoutHydrogen {
            bond: ring_bond,
            atom: ipso
        }
    );
    assert_eq!(
        error.to_string(),
        format!(
            "Cannot fuse a ring onto {ring_bond}. E{} carries no hydrogen.",
            &ipso.to_string()[1..]
        )
    );
}

#[test]
fn a_cis_double_bond_fuses_and_its_stereo_is_dropped() {
    // (Z)-2-butene with benzene gives o-xylene; the ring fixes the
    // configuration, so the double bond's stereo goes.
    let mut editor = butene(true);
    let bond = double_bond(&editor);
    assert!(editor.topology().pruned_stereo().is_empty());
    editor
        .fuse_fragment(bond, fragment("benzene").unwrap())
        .unwrap();
    assert!(matches!(
        editor.topology().pruned_stereo(),
        [PrunedStereo::DoubleBond { .. }]
    ));
    let structure = editor.finish().unwrap();
    let stereo = structure
        .topology()
        .definitions()
        .map(|definition| definition.molecule().stereo_elements().count())
        .sum::<usize>();
    assert_eq!(stereo, 0);
    assert_eq!(canonical(&structure), canonical_of("Cc1ccccc1C"));
    assert_sensible_geometry(&structure, "o-xylene");
    assert!(heavy_atoms_off_plane(&structure) < 0.05);
}

#[test]
fn cyclopropane_on_cyclohexane_aligns_with_the_hydrogens_it_replaces() {
    let mut editor = placed("cyclohexane");
    let bond = heavy_bond(&editor);
    let target = editor.topology().bond(bond).unwrap();
    let ends = [target.a(), target.b()];
    let [a, b] = ends.map(|atom| point(&editor, atom));
    let axis = (b - a) / (b - a).norm();
    let perpendicular = |v: Vector3| v - axis * v.dot(axis);
    let hydrogens = ends.map(|end| {
        editor
            .topology()
            .neighbors(end)
            .unwrap()
            .filter(|n| !is_heavy(&editor, *n))
            .map(|hydrogen| point(&editor, hydrogen))
            .collect::<Vec<_>>()
    });
    let reference = perpendicular(hydrogens[0][0] - a);
    let azimuth = |v: Vector3| {
        let v = perpendicular(v);
        axis.dot(reference.cross(v)).atan2(reference.dot(v))
    };
    let fused = editor
        .fuse_fragment(bond, fragment("cyclopropane").unwrap())
        .unwrap();
    // The aligned turns lay the ring bond along a hydrogen's place or halfway
    // between one of each bond atom; a 15Â° step would miss them.
    let sites = [a, b].map(|from| {
        hydrogens[usize::from(from != a)]
            .iter()
            .map(|hydrogen| azimuth(*hydrogen - from))
            .collect::<Vec<_>>()
    });
    let mut aligned = sites.iter().flatten().copied().collect::<Vec<_>>();
    for first in &sites[0] {
        for second in &sites[1] {
            aligned.push((first.sin() + second.sin()).atan2(first.cos() + second.cos()));
        }
    }
    let apex = fused
        .atoms()
        .iter()
        .map(|(_, atom)| *atom)
        .find(|atom| is_heavy(&editor, *atom))
        .unwrap();
    let found = azimuth(point(&editor, apex) - a);
    let miss = aligned
        .iter()
        .map(|angle| {
            let difference = found - angle;
            difference.sin().atan2(difference.cos()).abs().to_degrees()
        })
        .fold(f64::INFINITY, f64::min);
    assert!(miss < 0.5, "{miss:.1} degrees from an aligned turn");
    let structure = editor.finish().unwrap();
    assert_eq!(canonical(&structure), canonical_of("C1CCC2CC2C1"));
    assert_eq!(rings_and_aromatic_atoms(&structure), (vec![3, 6], 0));
    assert_sensible_geometry(&structure, "norcarane");
}

#[test]
fn pyrrole_on_cyclohexane_turns_the_bond_double() {
    for fixed in [false, true] {
        let mut editor = placed("cyclohexane");
        if fixed {
            let heavy = editor
                .topology()
                .atom_ids()
                .filter(|atom| is_heavy(&editor, *atom))
                .collect::<Vec<_>>();
            editor.remove_hydrogens(heavy.clone()).unwrap();
            for atom in &heavy {
                let carbon = editor
                    .topology()
                    .atom(*atom)
                    .unwrap()
                    .clone()
                    .with_hydrogens(ImplicitHydrogens::Fixed(2));
                editor.replace_atom(*atom, carbon).unwrap();
            }
        }
        let bond = heavy_bond(&editor);
        let fused = editor
            .fuse_fragment(bond, fragment("pyrrole").unwrap())
            .unwrap();
        // Each bond atom gives up both hydrogens to the aromatic ring.
        let removed = if fixed { 0 } else { 4 };
        assert_eq!(fused.removed_hydrogens().len(), removed);
        assert_eq!(
            editor.topology().bond(bond).unwrap().order,
            BondOrder::Double
        );
        if fixed {
            for (_, atom) in fused.fused_atoms() {
                assert_eq!(
                    editor.topology().atom(atom).unwrap().hydrogens(),
                    ImplicitHydrogens::Fixed(0)
                );
            }
        }
        let structure = editor.finish().unwrap();
        // Fixed counts write as bracket atoms, so compare inferred ones.
        if !fixed {
            assert_eq!(canonical(&structure), canonical_of("c1cc2c([nH]1)CCCC2"));
        }
        assert_eq!(rings_and_aromatic_atoms(&structure), (vec![5, 6], 5));
        let molecule = perceived_molecule(&structure);
        let hydrogens = molecule
            .atoms()
            .map(|(id, atom)| {
                usize::from(atom.element().atomic_number() == 1)
                    + molecule.implicit_hydrogens(id).unwrap()
            })
            .sum::<usize>();
        assert_eq!(hydrogens, 11, "fixed: {fixed}");
        assert_sensible_geometry(&structure, "tetrahydroindole");
    }
}

#[test]
fn a_ring_atom_takes_the_hydrogen_place_of_a_stereocenter() {
    // F[C@H](Cl)C(C)C, tetrahedral at both carbons, with the hydrogens of
    // the bond's atoms gauche.
    let reference = Point3::new(0.0, 0.1, 0.0);
    let center = Point3::origin();
    let next = Point3::new(0.153, 0.0, 0.0);
    let tetrahedral = 109.471_220_634;
    let hydrogen = bonded_point(reference, next, center, 0.109, tetrahedral, 120.0);
    let points = [
        bonded_point(reference, next, center, 0.135, tetrahedral, 0.0),
        center,
        bonded_point(reference, next, center, 0.177, tetrahedral, 240.0),
        next,
        bonded_point(hydrogen, center, next, 0.153, tetrahedral, 180.0),
        bonded_point(hydrogen, center, next, 0.153, tetrahedral, 300.0),
    ];
    let expected = canonical_of("F[C@]1(Cl)C(C)(C)CCCC1");
    assert_ne!(expected, canonical_of("F[C@@]1(Cl)C(C)(C)CCCC1"));
    for explicit in [false, true] {
        let mut editor = drafted("F[C@H](Cl)C(C)C", &points);
        let ids = editor.topology().atom_ids().collect::<Vec<_>>();
        if explicit {
            editor
                .add_hydrogens(ids.clone(), Default::default())
                .unwrap();
        }
        let bond = editor
            .topology()
            .bonds()
            .find(|(_, bond)| [bond.a(), bond.b()] == [ids[1], ids[3]])
            .unwrap()
            .0;
        let fused = editor
            .fuse_fragment(bond, fragment("cyclohexane").unwrap())
            .unwrap();
        let removed = if explicit { 2 } else { 0 };
        assert_eq!(fused.removed_hydrogens().len(), removed);
        // The ring bond lies along the hydrogen it replaces.
        let ring = editor
            .topology()
            .neighbors(ids[1])
            .unwrap()
            .find(|atom| fused.atoms().iter().any(|(_, new)| new == atom))
            .unwrap();
        let along = point(&editor, ring) - point(&editor, ids[1]);
        let toward = hydrogen - center;
        let angle = (along.dot(toward) / along.norm() / toward.norm())
            .acos()
            .to_degrees();
        assert!(angle < 15.0, "explicit: {explicit}, {angle:.1} degrees");
        let structure = editor.finish().unwrap();
        assert_eq!(canonical(&structure), expected, "explicit: {explicit}");
    }
}

#[test]
fn benzene_on_the_nh_bond_of_imidazole_gives_imidazopyridine() {
    let mut editor = placed("imidazole");
    let bond = named_bond(&editor, "N1", "C2");
    let fused = editor
        .fuse_fragment(bond, fragment("benzene").unwrap())
        .unwrap();
    // The nitrogen and the carbon each give up their hydrogen.
    assert_eq!(fused.removed_hydrogens().len(), 2);
    let structure = editor.finish().unwrap();
    assert_eq!(canonical(&structure), canonical_of("c1ccn2ccnc2c1"));
    assert_eq!(rings_and_aromatic_atoms(&structure), (vec![5, 6], 9));
    assert_sensible_geometry(&structure, "imidazo[1,2-a]pyridine");
    assert!(heavy_atoms_off_plane(&structure) < 0.05);
}
