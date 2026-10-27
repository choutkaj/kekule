//! Failures name every atom they concern, and the charge model's element
//! domain is readable before parameterizing.

use kekule::molecule::{AtomId, Element, Molecule};
use kekule::smiles;
use kekule::topology::TopologyBuilder;
use kekule_openff::{ChargeMethod, ErrorKind, ForceField, NaglModel};
use std::sync::Arc;

// Only water's bond and angle have parameters, so other molecules are left
// unassigned.
const WATER: &str = r#"<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
<Bonds version="0.4"><Bond smirks="[#8:1]-[#1:2]" id="water-bond" length="0.1*nanometer" k="100*kilojoule_per_mole/nanometer**2"/></Bonds>
<Angles version="0.3"><Angle smirks="[*:1]~[*:2]~[*:3]" id="any-angle" angle="90*degree" k="100*kilojoule_per_mole/radian**2"/></Angles>
<ProperTorsions version="0.4"/>
<vdW version="0.4"><Atom smirks="[*:1]" id="any-vdw" sigma="0.3*nanometer" epsilon="0.1*kilojoule_per_mole"/></vdW>
<Electrostatics version="0.4" scale14="0.8333333333"/>
</SMIRNOFF>"#;

fn molecule(text: &str) -> Molecule {
    smiles::to_molecules(text).unwrap().remove(0)
}

fn ids(indices: &[u32]) -> Vec<AtomId> {
    indices.iter().map(|&i| AtomId::new(i)).collect()
}

#[test]
fn elements_outside_the_charge_model_name_every_such_atom() {
    let ash = NaglModel::ash().unwrap();
    let magnesium = Element::from_symbol("Mg").unwrap();
    // H0 Mg1 O2 Mg3 H4
    let error = ash
        .assign_charges(&molecule("[H][Mg]O[Mg][H]"))
        .unwrap_err();
    assert_eq!(error.kind(), ErrorKind::UnsupportedMolecule);
    assert_eq!(error.atoms(), ids(&[1, 3]));
    assert_eq!(
        error.to_string(),
        "element Mg is outside the model domain: atoms 1, 3"
    );
    assert!(!ash.supports_element(magnesium));

    let elements = ash.elements();
    for symbol in ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"] {
        let element = Element::from_symbol(symbol).unwrap();
        assert!(elements.contains(&element), "{symbol}");
        assert!(ash.supports_element(element));
    }
    assert!(!elements.contains(&magnesium));
    assert!(elements
        .windows(2)
        .all(|w| w[0].atomic_number() < w[1].atomic_number()));
}

#[test]
fn system_failures_name_atoms_of_the_failing_definition() {
    let mut builder = TopologyBuilder::new();
    let water = builder
        .add_molecule_definition(molecule("[H]O[H]"))
        .unwrap();
    let magnesium = builder
        .add_molecule_definition(molecule("[H][Mg][H]"))
        .unwrap();
    builder.add_instance(water).unwrap();
    builder.add_instance(magnesium).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let ash = NaglModel::ash().unwrap();
    let error = ForceField::rosemary()
        .unwrap()
        .parameterize(topology, &ash)
        .unwrap_err();
    assert_eq!(error.definition(), Some(magnesium));
    assert!(error.atoms().contains(&AtomId::new(1)), "{error}");
}

#[test]
fn unassigned_interactions_name_all_their_atoms() {
    let ff = ForceField::from_offxml(WATER).unwrap();
    // H0 C1 H2 H3 O4 H5: three C-H bonds and C-O lack parameters; O-H has one.
    let error = ff
        .label_molecule(&molecule("[H]C([H])([H])O[H]"))
        .unwrap_err();
    assert_eq!(error.kind(), ErrorKind::Unparameterized);
    assert_eq!(error.atoms(), ids(&[0, 1, 2, 3, 4]));
    assert_eq!(
        error.detail(),
        "4 unassigned Bonds interactions, the first of atoms 0, 1"
    );
}

#[test]
fn implicit_hydrogens_and_radicals_name_their_atoms() {
    let ff = ForceField::rosemary().unwrap();
    // C0 O1 C2: both carbons have implicit hydrogens.
    let error = ff.label_molecule(&molecule("COC")).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::UnsupportedMolecule);
    assert_eq!(error.atoms(), ids(&[0, 2]));
    assert_eq!(
        error.detail(),
        "atoms 0, 2 have implicit hydrogens; expand hydrogens before parameterization"
    );
    // H0 C1 H2 H3 C4 H5 H6: C4 is a radical.
    let error = ff
        .label_molecule(&molecule("[H]C([H])([H])[C]([H])[H]"))
        .unwrap_err();
    assert_eq!(error.kind(), ErrorKind::UnsupportedMolecule);
    assert_eq!(error.atoms(), ids(&[4]), "{error}");
}

#[test]
fn failures_without_particular_atoms_name_none() {
    let error = ForceField::from_offxml("<SMIRNOFF").unwrap_err();
    assert!(error.atoms().is_empty());
    let mut water = molecule("O");
    water.perceive().unwrap();
    water.add_hydrogens(Default::default()).unwrap();
    let error = ForceField::from_offxml(WATER)
        .unwrap()
        .parameterize_molecule(water, ChargeMethod::LibraryOnly)
        .unwrap_err();
    assert_eq!(error.kind(), ErrorKind::Charges);
    assert!(error.atoms().is_empty());
}
