//! The mmCIF writer's chemistry policy.

use kekule::diagnostics::{
    AtomField, FormatError, FormatErrorDetail, UnrepresentablePolicy, WriteOmission,
};
use kekule::math::Point3;
use kekule::mmcif::{self, MmcifWriteErrorKind, MmcifWriteOptions};
use kekule::molecule::{Atom, AtomRadical, BondOrder, Element, MoleculeEditor};
use kekule::sdf;
use kekule::stereo::coordinates::{infer_coordinate_stereo, CoordinateStereoOptions};
use kekule::structure::{Positions, Structure, StructureSet};
use kekule::topology::Topology;
use kekule::units::{Quantity, ANGSTROM};

/// A 3D glucoside with five stereocenters and the chiral flag set.
const GLUCOSIDE: &str = "
  kekule            3D

 32 33  0  0  1  0            999 V2000
  -19.6811    9.0976   -1.3097 C   0  0  0  0  0  0  0  0  0  0  0  0
  -18.1470    9.3470   -1.1788 C   0  0  0  0  0  0  0  0  0  0  0  0
  -17.3557    8.0491   -1.4963 C   0  0  0  0  0  0  0  0  0  0  0  0
  -17.9117    6.8503   -0.6932 C   0  0  0  0  0  0  0  0  0  0  0  0
  -19.4213    6.7153   -0.9860 C   0  0  0  0  0  0  0  0  0  0  0  0
  -20.0318    5.5069   -0.2506 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.7874    8.0704   -1.7876 N   0  0  0  0  0  0  0  0  0  0  0  0
  -20.6302   10.5227   -0.6788 S   0  0  0  0  0  0  0  0  0  0  0  0
  -20.0677    7.9264   -0.5810 O   0  0  0  0  0  0  0  0  0  0  0  0
  -21.3398    5.2784   -0.7128 O   0  0  0  0  0  0  0  0  0  0  0  0
  -17.6819    6.9561    0.7028 O   0  0  0  0  0  0  0  0  0  0  0  0
  -17.7312   10.4012   -2.0300 O   0  0  0  0  0  0  0  0  0  0  0  0
  -15.9069    8.1962   -1.2925 N   0  0  0  0  0  0  0  0  0  0  0  0
  -15.0649    7.9011   -2.2872 N   0  0  0  0  0  0  0  0  0  0  0  0
  -13.8629    8.4275   -0.5134 C   0  0  0  0  0  0  0  0  0  0  0  0
  -15.2059    8.5375   -0.1613 C   0  0  0  0  0  0  0  0  0  0  0  0
  -19.9625    8.9765   -2.3582 H   0  0  0  0  0  0  0  0  0  0  0  0
  -17.9283    9.6594   -0.1559 H   0  0  0  0  0  0  0  0  0  0  0  0
  -17.5244    7.8281   -2.5539 H   0  0  0  0  0  0  0  0  0  0  0  0
  -17.4037    5.9392   -1.0109 H   0  0  0  0  0  0  0  0  0  0  0  0
  -19.5784    6.5669   -2.0566 H   0  0  0  0  0  0  0  0  0  0  0  0
  -20.0459    5.6564    0.8295 H   0  0  0  0  0  0  0  0  0  0  0  0
  -19.4469    4.6074   -0.4456 H   0  0  0  0  0  0  0  0  0  0  0  0
  -21.9273    5.9081   -0.2714 H   0  0  0  0  0  0  0  0  0  0  0  0
  -16.8294    6.5543    0.9071 H   0  0  0  0  0  0  0  0  0  0  0  0
  -16.8420   10.6641   -1.7820 H   0  0  0  0  0  0  0  0  0  0  0  0
  -15.7056    8.7879    0.7629 H   0  0  0  0  0  0  0  0  0  0  0  0
  -22.3820   10.1958   -0.8474 C   0  0  0  0  0  0  0  0  0  0  0  0
  -22.9489   11.0471   -0.4705 H   0  0  0  0  0  0  0  0  0  0  0  0
  -22.6228   10.0361   -1.8984 H   0  0  0  0  0  0  0  0  0  0  0  0
  -22.6425    9.3048   -0.2761 H   0  0  0  0  0  0  0  0  0  0  0  0
  -12.9600    8.5799    0.0592 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  2  3  1  1  0  0  0
  3  4  1  1  0  0  0
  4  5  1  1  0  0  0
  5  6  1  1  0  0  0
  8  1  1  0  0  0  0
  9  1  1  0  0  0  0
  9  5  1  0  0  0  0
 10  6  1  0  0  0  0
 11  4  1  0  0  0  0
 12  2  1  0  0  0  0
 13  3  1  0  0  0  0
 14  7  2  0  0  0  0
 14 13  1  0  0  0  0
 15  7  1  0  0  0  0
 16 13  1  0  0  0  0
 16 15  2  0  0  0  0
  1 17  1  0  0  0  0
  2 18  1  0  0  0  0
  3 19  1  0  0  0  0
  4 20  1  0  0  0  0
  5 21  1  0  0  0  0
 22  6  1  0  0  0  0
 23  6  1  0  0  0  0
 24 10  1  0  0  0  0
 25 11  1  0  0  0  0
 26 12  1  0  0  0  0
 27 16  1  0  0  0  0
 28  8  1  0  0  0  0
 29 28  1  0  0  0  0
 30 28  1  0  0  0  0
 31 28  1  0  0  0  0
 32 15  1  0  0  0  0
M  END
$$$$
";

fn issue<T: std::fmt::Debug>(result: Result<T, FormatError>) -> MmcifWriteErrorKind {
    match result.unwrap_err().detail() {
        Some(FormatErrorDetail::MmcifWrite(issue)) => issue.clone(),
        other => panic!("not an mmCIF write issue: {other:?}"),
    }
}

fn glucoside() -> Structure {
    sdf::parse_str(GLUCOSIDE, Default::default())
        .unwrap()
        .records()[0]
        .interpret(Default::default())
        .unwrap()
        .structure()
        .clone()
}

fn read_back(text: &str) -> Structure {
    mmcif::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_structure()
}

fn omit_chemistry() -> MmcifWriteOptions {
    MmcifWriteOptions::default().with_chemistry(UnrepresentablePolicy::Omit)
}

#[test]
fn chemistry_fails_by_default() {
    assert_eq!(
        MmcifWriteOptions::default().chemistry(),
        UnrepresentablePolicy::Fail
    );
    let structure = glucoside();
    let error = mmcif::write([&structure], MmcifWriteOptions::default()).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Cannot write the mmCIF output. The stereochemistry of instance 0 has no mmCIF representation."
    );
    assert!(matches!(
        issue(mmcif::write([&structure], MmcifWriteOptions::default())),
        MmcifWriteErrorKind::UnsupportedStereo(_)
    ));
}

#[test]
fn omitted_stereo_reads_back_from_the_coordinates() {
    let structure = glucoside();
    let original = structure
        .topology()
        .instances()
        .next()
        .unwrap()
        .definition()
        .molecule()
        .clone();
    assert_eq!(original.stereo_elements().count(), 5);
    let options = omit_chemistry().with_coordinate_precision(4).unwrap();
    let (text, report) = mmcif::write([&structure], options).unwrap();
    assert_eq!(report.omissions(), [WriteOmission::Stereo { molecules: 1 }]);
    // Models share one topology, so a block reports its stereo once.
    let set = StructureSet::from_structures([structure.clone(), structure.clone()]).unwrap();
    let (_, report) = mmcif::write([&set], omit_chemistry()).unwrap();
    assert_eq!(report.omissions(), [WriteOmission::Stereo { molecules: 1 }]);

    let read = read_back(&text);
    assert_eq!(read.topology().atom_count(), 32);
    assert_eq!(read.topology().bond_count(), 33);
    for (written, read) in structure
        .positions()
        .values()
        .value()
        .iter()
        .zip(read.positions().values().value().iter())
    {
        assert!((*written - *read).norm() < 1e-9, "{written:?} {read:?}");
    }
    // The reader gives no stereo; the coordinates give back every element.
    let molecule = read
        .topology()
        .instances()
        .next()
        .unwrap()
        .definition()
        .molecule()
        .clone();
    assert_eq!(molecule.stereo_elements().count(), 0);
    let inferred = infer_coordinate_stereo(
        &molecule,
        read.positions(),
        CoordinateStereoOptions::default(),
    )
    .unwrap();
    assert_eq!(inferred.elements().len(), 5);
    for (_, element) in original.stereo_elements() {
        assert!(inferred.elements().contains(element), "{element:?}");
    }
}

fn atom_structure(atom: Atom) -> Structure {
    let mut editor = MoleculeEditor::new();
    editor.add_atom(atom).unwrap();
    let topology = Topology::from_molecule(editor.finish().unwrap()).unwrap();
    Structure::new(topology, Positions::zeros(1)).unwrap()
}

#[test]
fn radicals_and_isotopes_fail_unless_omitted() {
    let carbon = || Atom::new(Element::from_symbol("C").unwrap());
    let radical = atom_structure(carbon().with_radical(AtomRadical::new(1, Some(2))));
    let isotope = atom_structure(carbon().with_isotope(std::num::NonZeroU16::new(13)));
    for (structure, field, omission) in [
        (
            &radical,
            AtomField::Radical,
            WriteOmission::Radicals { atoms: 1 },
        ),
        (
            &isotope,
            AtomField::Isotope,
            WriteOmission::Isotopes { atoms: 1 },
        ),
    ] {
        assert!(matches!(
            issue(mmcif::write([structure], MmcifWriteOptions::default())),
            MmcifWriteErrorKind::UnsupportedAtomField { field: found, .. } if found == field
        ));
        // Omitting implicit hydrogens does not omit other chemistry.
        let hydrogens =
            MmcifWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
        assert!(mmcif::write([structure], hydrogens).is_err());
        let (text, report) = mmcif::write([structure], omit_chemistry()).unwrap();
        assert_eq!(report.omissions(), [omission]);
        // The atom reads back as plain carbon.
        let read = read_back(&text);
        let atom = read.topology().atoms().next().unwrap().atom().clone();
        assert_eq!(atom.element().symbol(), "C");
        assert_eq!((atom.isotope(), atom.radical()), (None, None));
    }
    // Deuterium has an mmCIF symbol, so the policy does not touch it.
    let deuterium = atom_structure(
        Atom::new(Element::from_symbol("H").unwrap()).with_isotope(std::num::NonZeroU16::new(2)),
    );
    let (text, report) = mmcif::write([&deuterium], MmcifWriteOptions::default()).unwrap();
    assert!(report.is_empty());
    assert!(text.contains(" D "), "{text}");
}

#[test]
fn atom_maps_and_zero_order_bonds_fail_unless_omitted() {
    let mapped = atom_structure(
        Atom::new(Element::from_symbol("O").unwrap()).with_atom_map(std::num::NonZeroU32::new(1)),
    );
    assert!(matches!(
        issue(mmcif::write([&mapped], MmcifWriteOptions::default())),
        MmcifWriteErrorKind::UnsupportedAtomField {
            field: AtomField::AtomMap,
            ..
        }
    ));
    let (_, report) = mmcif::write([&mapped], omit_chemistry()).unwrap();
    assert_eq!(report.omissions(), [WriteOmission::AtomMaps { atoms: 1 }]);

    let mut editor = MoleculeEditor::new();
    let ions = [("Na", 1), ("Cl", -1)].map(|(symbol, charge)| {
        let atom = Atom::new(Element::from_symbol(symbol).unwrap()).with_formal_charge(charge);
        editor.add_atom(atom).unwrap()
    });
    editor.add_bond(ions[0], ions[1], BondOrder::Zero).unwrap();
    let topology = Topology::from_molecule(editor.finish().unwrap()).unwrap();
    let positions = [Point3::new(0.0, 0.0, 0.0), Point3::new(8.0, 0.0, 0.0)];
    let positions = Positions::new(Quantity::new(positions.to_vec(), ANGSTROM)).unwrap();
    let pair = Structure::new(topology, positions).unwrap();
    assert!(matches!(
        issue(mmcif::write([&pair], MmcifWriteOptions::default())),
        MmcifWriteErrorKind::UnsupportedBondOrder {
            order: BondOrder::Zero,
            ..
        }
    ));
    let (text, report) = mmcif::write([&pair], omit_chemistry()).unwrap();
    assert_eq!(report.omissions(), [WriteOmission::Bonds { bonds: 1 }]);
    assert_eq!(read_back(&text).topology().bond_count(), 0);
}
