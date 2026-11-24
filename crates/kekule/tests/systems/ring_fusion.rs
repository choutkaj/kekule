//! Fusing ring fragments onto bonds of structure drafts.
use kekule::math::Point3;
use kekule::molecule::{BondOrder, ImplicitHydrogens, Molecule};
use kekule::smiles::{self, SmilesWriteOptions};
use kekule::structure::placement::{
    fragment, fragments, ideal_bond_length, FragmentKind, PlacementError,
};
use kekule::structure::{Structure, StructureEditor};
use kekule::topology::edit::{EditAtomId, EditBondId};
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

/// The position of the atom at dense `index`, in Å.
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
/// three or more bonds apart keep 2.3 Å from each other between heavy atoms and
/// 1.6 Å otherwise; unrelaxed bay hydrogens, as in phenanthrene, come to 1.75 Å.
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
/// three heavy atoms, in Å.
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

/// A draft's atoms (element, hydrogen declaration, position, name) and
/// bonds (dense atoms and order) once published.
type Snapshot = (
    Vec<(u8, String, Point3, String)>,
    Vec<(usize, usize, String)>,
);

fn snapshot(editor: StructureEditor) -> Snapshot {
    let structure = editor.finish().unwrap();
    let topology = structure.topology();
    let atoms = topology
        .atoms()
        .map(|atom| {
            let site = atom.atom_site().unwrap();
            (
                atom.element().atomic_number(),
                format!("{:?}", atom.atom().hydrogens()),
                angstroms(&structure, atom.index().index()),
                site.metadata().label_atom_id().unwrap().to_owned(),
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
    (atoms, bonds)
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
    // Both Kekulé orders of the target bond are covered.
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
