use super::*;
use crate::molecule::{Atom, AtomRadical, BondOrder, ImplicitHydrogens};

fn element(symbol: &str) -> Element {
    Element::from_symbol(symbol).expect("test element")
}

fn mass_value(value: Quantity<f64>) -> f64 {
    *value.value()
}

#[test]
fn stored_and_perceived_hydrogens_are_explicit_policies() {
    let mut molecule = crate::test_fixtures::read_smiles("C").expect("methane parses");
    assert_eq!(
        molecular_formula(&molecule, HydrogenCountPolicy::StoredOnly)
            .expect("stored formula")
            .to_string(),
        "C"
    );
    assert_eq!(
        molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived),
        Err(MolecularDescriptorError::UnknownHydrogenCount(
            UnknownHydrogenCount::new(AtomId::new(0))
        ))
    );

    molecule.perceive().expect("methane perceives");
    assert_eq!(
        molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived)
            .expect("perceived formula")
            .to_string(),
        "CH4"
    );
}

#[test]
fn stored_only_counts_fixed_hydrogens_and_omits_inferred_ones() {
    let mut fixed = Atom::new(element("C"));
    fixed.hydrogens = ImplicitHydrogens::Fixed(3);
    let mut graph = crate::molecule::MoleculeEditor::new();
    let fixed = graph.add_atom(fixed).expect("fixed carbon");
    let inferred = graph.add_atom(Atom::new(element("C"))).expect("carbon");
    graph
        .add_bond(fixed, inferred, crate::molecule::BondOrder::Single)
        .expect("bond");
    let mut molecule = graph.finish().expect("ethane");

    assert_eq!(
        molecular_formula(&molecule, HydrogenCountPolicy::StoredOnly)
            .expect("stored formula")
            .to_string(),
        "C2H3"
    );
    molecule.perceive().expect("ethane perceives");
    assert_eq!(
        molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived)
            .expect("complete formula")
            .to_string(),
        "C2H6"
    );
}

#[test]
fn hill_order_isotope_order_and_accessors_are_stable() {
    let mut graph = crate::molecule::MoleculeEditor::new();
    let mut carbon_13_atom = Atom::new(element("C"));
    carbon_13_atom.isotope = std::num::NonZeroU16::new(13);
    let mut carbon_12_atom = carbon_13_atom.clone();
    carbon_12_atom.isotope = std::num::NonZeroU16::new(12);
    let carbon_13 = graph
        .add_atom(carbon_13_atom)
        .expect("atom identifier capacity");
    let oxygen = graph
        .add_atom(Atom::new(element("O")))
        .expect("atom identifier capacity");
    let carbon = graph
        .add_atom(Atom::new(element("C")))
        .expect("atom identifier capacity");
    let carbon_12 = graph
        .add_atom(carbon_12_atom)
        .expect("atom identifier capacity");
    graph
        .add_bond(carbon_13, oxygen, BondOrder::Single)
        .expect("first connecting bond");
    graph
        .add_bond(oxygen, carbon, BondOrder::Single)
        .expect("second connecting bond");
    graph
        .add_bond(carbon, carbon_12, BondOrder::Single)
        .expect("third connecting bond");
    let graph = graph.into_unchecked();
    let formula = molecular_formula(&graph, HydrogenCountPolicy::StoredOnly).expect("formula");

    assert_eq!(formula.to_string(), "C[12C][13C]O");
    assert_eq!(formula.count(element("C")), 3);
    assert_eq!(formula.unlabeled_count(element("C")), 1);
    assert_eq!(formula.isotope_count(element("C"), 13), 1);
    assert_eq!(
        formula.terms().collect::<Vec<_>>(),
        vec![
            (element("C"), None, 1),
            (element("C"), Some(12), 1),
            (element("C"), Some(13), 1),
            (element("O"), None, 1),
        ]
    );
}

#[test]
fn no_carbon_formula_and_charge_are_reported_per_connected_entity() {
    let document =
        crate::smiles::parse_str("[NH4+].[Cl-]", Default::default()).expect("salt parses");
    let components = document
        .interpret(Default::default())
        .expect("salt interprets")
        .into_molecules();
    let formulas = components
        .iter()
        .map(|component| {
            molecular_formula(component, HydrogenCountPolicy::StoredOnly)
                .expect("component formula")
                .to_string()
        })
        .collect::<Vec<_>>();
    assert_eq!(formulas, vec!["H4N+", "Cl-"]);

    let ammonium = crate::test_fixtures::read_smiles("[NH4+]").expect("ammonium parses");
    assert_eq!(
        molecular_formula(&ammonium, HydrogenCountPolicy::StoredOnly)
            .expect("formula")
            .to_string(),
        "H4N+"
    );

    let chloride = crate::test_fixtures::read_smiles("[Cl-]").expect("chloride parses");
    assert_eq!(
        molecular_formula(&chloride, HydrogenCountPolicy::StoredOnly)
            .expect("formula")
            .to_string(),
        "Cl-"
    );

    let iron = crate::test_fixtures::read_smiles("[Fe+2]").expect("iron ion parses");
    assert_eq!(
        molecular_formula(&iron, HydrogenCountPolicy::StoredOnly)
            .expect("formula")
            .to_string(),
        "Fe+2"
    );
}

#[test]
fn atomic_masses_are_neutral_and_do_not_require_a_molecule() {
    let carbon = element("C");
    assert_eq!(mass_value(atomic_mass(carbon, Some(12)).unwrap()), 12.0);
    assert!((mass_value(atomic_mass(carbon, Some(13)).unwrap()) - 13.003_354_835).abs() < 1e-7);
    assert!(atomic_mass(carbon, Some(999)).is_none());
}

#[test]
fn average_and_monoisotopic_mass_use_pinned_tables() {
    let mut water = crate::test_fixtures::read_smiles("O").expect("water parses");
    water.perceive().expect("water perceives");
    let average = mass_value(
        average_mass(&water, HydrogenCountPolicy::IncludePerceived).expect("average mass"),
    );
    let monoisotopic = mass_value(
        monoisotopic_mass(&water, HydrogenCountPolicy::IncludePerceived)
            .expect("monoisotopic mass"),
    );
    assert!((average - 18.015).abs() < 1.0e-12);
    assert!((monoisotopic - 18.010_564_683_056).abs() < 1.0e-12);

    let carbon_13 = crate::test_fixtures::read_smiles("[13C]").expect("isotope parses");
    assert!(
        (mass_value(
            monoisotopic_mass(&carbon_13, HydrogenCountPolicy::StoredOnly).expect("isotope mass")
        ) - 13.003_354_835_34)
            .abs()
            < 1.0e-12
    );
}

#[test]
fn formal_charge_applies_electron_mass_correction() {
    let sodium = crate::test_fixtures::read_smiles("[Na+]").expect("sodium parses");
    let chloride = crate::test_fixtures::read_smiles("[Cl-]").expect("chloride parses");
    let sodium_mass =
        mass_value(average_mass(&sodium, HydrogenCountPolicy::StoredOnly).expect("sodium mass"));
    let chloride_mass = mass_value(
        average_mass(&chloride, HydrogenCountPolicy::StoredOnly).expect("chloride mass"),
    );
    assert!((sodium_mass - (22.990 - data::ELECTRON_MASS_DA)).abs() < 1.0e-12);
    assert!((chloride_mass - (35.45 + data::ELECTRON_MASS_DA)).abs() < 1.0e-12);
}

#[test]
fn radicals_need_no_mass_correction_beyond_charge() {
    let mut graph = crate::molecule::MoleculeEditor::new();
    let mut carbon = Atom::new(element("C"));
    carbon.radical = AtomRadical::new(1, Some(2));
    graph.add_atom(carbon).expect("atom identifier capacity");
    let molecule = graph;
    let molecule = molecule.into_unchecked();
    assert_eq!(
        molecular_formula(&molecule, HydrogenCountPolicy::StoredOnly)
            .expect("formula")
            .to_string(),
        "C"
    );
    assert_eq!(
        mass_value(
            monoisotopic_mass(&molecule, HydrogenCountPolicy::StoredOnly).expect("radical mass")
        ),
        12.0
    );
}

#[test]
fn unavailable_weights_and_isotopes_are_structured_errors() {
    let technetium = crate::test_fixtures::read_smiles("[Tc]").expect("technetium parses");
    assert!(matches!(
        average_mass(&technetium, HydrogenCountPolicy::StoredOnly),
        Err(MolecularDescriptorError::MissingStandardAtomicWeight {
            element: actual_element,
            ..
        }) if actual_element == element("Tc")
    ));
    assert!(matches!(
        monoisotopic_mass(&technetium, HydrogenCountPolicy::StoredOnly),
        Err(MolecularDescriptorError::MissingNaturalIsotope {
            element: actual_element,
            ..
        }) if actual_element == element("Tc")
    ));

    let technetium_99 =
        crate::test_fixtures::read_smiles("[99Tc]").expect("technetium isotope parses");
    let average = mass_value(
        average_mass(&technetium_99, HydrogenCountPolicy::StoredOnly)
            .expect("labeled isotope has an AME mass"),
    );
    let monoisotopic = mass_value(
        monoisotopic_mass(&technetium_99, HydrogenCountPolicy::StoredOnly)
            .expect("labeled isotope has an AME mass"),
    );
    assert_eq!(average, monoisotopic);
    assert!((average - 98.906_249_681).abs() < 1.0e-12);

    let mut graph = crate::molecule::MoleculeEditor::new();
    let mut impossible = Atom::new(element("C"));
    impossible.isotope = std::num::NonZeroU16::new(999);
    graph
        .add_atom(impossible)
        .expect("atom identifier capacity");
    let graph = graph.into_unchecked();
    assert!(matches!(
        monoisotopic_mass(&graph, HydrogenCountPolicy::StoredOnly),
        Err(MolecularDescriptorError::UnknownIsotope {
            mass_number: 999,
            ..
        })
    ));
}

#[test]
fn empty_molecule_has_empty_formula_and_zero_masses() {
    let molecule = crate::molecule::MoleculeEditor::new();
    let molecule = molecule.into_unchecked();
    let formula = molecular_formula(&molecule, HydrogenCountPolicy::StoredOnly).expect("formula");
    assert!(formula.is_empty());
    assert_eq!(formula.to_string(), "");
    assert_eq!(
        mass_value(average_mass(&molecule, HydrogenCountPolicy::StoredOnly).expect("mass")),
        0.0
    );
    assert_eq!(
        mass_value(monoisotopic_mass(&molecule, HydrogenCountPolicy::StoredOnly).expect("mass")),
        0.0
    );
}

#[test]
fn formula_count_overflow_is_checked() {
    let mut counts = BTreeMap::from([((element("H"), None), u64::MAX)]);
    assert_eq!(
        add_count(&mut counts, element("H"), None, 1),
        Err(MolecularDescriptorError::CountOverflow {
            element: element("H"),
            mass_number: None,
        })
    );
}
