//! Molecular formulas and masses.
use kekule::descriptors::{
    average_mass, molecular_formula, monoisotopic_mass, MolecularDescriptorError,
};
use kekule::diagnostics::UnknownHydrogenCount;
use kekule::molecule::{
    Atom, AtomId, BondOrder, Element, ImplicitHydrogens, Molecule, MoleculeEditor,
};
use kekule::smiles;
use kekule::units::Quantity;

fn element(symbol: &str) -> Element {
    Element::from_symbol(symbol).expect("test element")
}

fn mass_value(value: Quantity<f64>) -> f64 {
    *value.value()
}

fn read_smiles(source: &str) -> Result<Molecule, Box<dyn std::error::Error>> {
    let mut molecules = smiles::parse_str(source, Default::default())?
        .interpret(Default::default())?
        .into_molecules();
    assert_eq!(molecules.len(), 1, "{source}");
    Ok(molecules.pop().unwrap())
}

#[test]
fn inferred_hydrogens_are_counted_or_the_formula_fails() {
    // An inferred count is never omitted: before perception the formula and
    // masses fail rather than describe methane as "C".
    let mut molecule = read_smiles("C").expect("methane parses");
    let unknown = Err(MolecularDescriptorError::UnknownHydrogenCount(
        UnknownHydrogenCount::new(AtomId::new(0)),
    ));
    assert_eq!(molecular_formula(&molecule), unknown);
    assert!(matches!(
        average_mass(&molecule),
        Err(MolecularDescriptorError::UnknownHydrogenCount(_))
    ));
    assert!(matches!(
        monoisotopic_mass(&molecule),
        Err(MolecularDescriptorError::UnknownHydrogenCount(_))
    ));

    molecule.perceive().expect("methane perceives");
    assert_eq!(
        molecular_formula(&molecule)
            .expect("perceived formula")
            .to_string(),
        "CH4"
    );
}

#[test]
fn fixed_counts_need_no_perception_and_inferred_ones_do() {
    let fixed = Atom::new(element("C")).with_hydrogens(ImplicitHydrogens::Fixed(3));
    let mut graph = MoleculeEditor::new();
    let fixed = graph.add_atom(fixed).expect("fixed carbon");
    let inferred = graph.add_atom(Atom::new(element("C"))).expect("carbon");
    graph
        .add_bond(fixed, inferred, BondOrder::Single)
        .expect("bond");
    let mut molecule = graph.finish().expect("ethane");

    assert_eq!(
        molecular_formula(&molecule),
        Err(MolecularDescriptorError::UnknownHydrogenCount(
            UnknownHydrogenCount::new(inferred)
        ))
    );
    molecule.perceive().expect("ethane perceives");
    assert_eq!(
        molecular_formula(&molecule)
            .expect("complete formula")
            .to_string(),
        "C2H6"
    );
}

#[test]
fn no_carbon_formula_and_charge_are_reported_per_connected_entity() {
    let document = smiles::parse_str("[NH4+].[Cl-]", Default::default()).expect("salt parses");
    let components = document
        .interpret(Default::default())
        .expect("salt interprets")
        .into_molecules();
    let formulas = components
        .iter()
        .map(|component| {
            molecular_formula(component)
                .expect("component formula")
                .to_string()
        })
        .collect::<Vec<_>>();
    assert_eq!(formulas, vec!["H4N+", "Cl-"]);

    let ammonium = read_smiles("[NH4+]").expect("ammonium parses");
    assert_eq!(
        molecular_formula(&ammonium).expect("formula").to_string(),
        "H4N+"
    );

    let chloride = read_smiles("[Cl-]").expect("chloride parses");
    assert_eq!(
        molecular_formula(&chloride).expect("formula").to_string(),
        "Cl-"
    );

    let iron = read_smiles("[Fe+2]").expect("iron ion parses");
    assert_eq!(
        molecular_formula(&iron).expect("formula").to_string(),
        "Fe+2"
    );
}

#[test]
fn average_and_monoisotopic_mass_use_pinned_tables() {
    let mut water = read_smiles("O").expect("water parses");
    water.perceive().expect("water perceives");
    let average = mass_value(average_mass(&water).expect("average mass"));
    let monoisotopic = mass_value(monoisotopic_mass(&water).expect("monoisotopic mass"));
    assert!((average - 18.015).abs() < 1.0e-12);
    assert!((monoisotopic - 18.010_564_683_056).abs() < 1.0e-12);

    let carbon_13 = read_smiles("[13C]").expect("isotope parses");
    assert!(
        (mass_value(monoisotopic_mass(&carbon_13).expect("isotope mass")) - 13.003_354_835_34)
            .abs()
            < 1.0e-12
    );
}
