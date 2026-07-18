//! Read-only small-molecule formula and mass descriptors.
//!
//! Descriptor calculation never mutates, sanitizes, or perceives the input.
//! Callers select explicitly whether installed implicit-hydrogen state is part
//! of the calculation.

mod data;

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

use crate::core::{AtomId, Element};
use crate::small::SmallMolecule;
use crate::units::{Quantity, DALTON};

use data::{exact_isotope_mass, most_abundant_isotope, standard_atomic_weight};

/// Selects which non-atom hydrogen counts contribute to a descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HydrogenCountPolicy {
    /// Count live hydrogen atoms and stored explicit-hydrogen declarations.
    StoredOnly,
    /// Also count installed implicit-hydrogen perception state.
    IncludePerceived,
}

/// A stable, structured molecular formula.
///
/// Terms are stored in Hill order and keep unlabeled and isotope-labeled atoms
/// separate. The aggregate formal charge is the sum over all live atoms.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolecularFormula {
    terms: Vec<FormulaTerm>,
    formal_charge: i64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct FormulaTerm {
    element: Element,
    isotope: Option<u16>,
    count: u64,
}

impl MolecularFormula {
    /// Returns the total count of an element across unlabeled and labeled terms.
    pub fn count(&self, element: Element) -> u64 {
        self.terms
            .iter()
            .filter(|term| term.element == element)
            .map(|term| term.count)
            .sum()
    }

    /// Returns the count of an unlabeled element.
    pub fn unlabeled_count(&self, element: Element) -> u64 {
        self.terms
            .iter()
            .find(|term| term.element == element && term.isotope.is_none())
            .map_or(0, |term| term.count)
    }

    /// Returns the count of one explicitly labeled isotope.
    pub fn isotope_count(&self, element: Element, mass_number: u16) -> u64 {
        self.terms
            .iter()
            .find(|term| term.element == element && term.isotope == Some(mass_number))
            .map_or(0, |term| term.count)
    }

    /// Iterates `(element, isotope, count)` terms in stable Hill order.
    pub fn terms(&self) -> impl ExactSizeIterator<Item = (Element, Option<u16>, u64)> + '_ {
        self.terms
            .iter()
            .map(|term| (term.element, term.isotope, term.count))
    }

    /// Returns the aggregate asserted formal charge.
    pub const fn formal_charge(&self) -> i64 {
        self.formal_charge
    }

    /// Returns whether the formula contains no atoms or hydrogen declarations.
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

impl fmt::Display for MolecularFormula {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        for term in &self.terms {
            match term.isotope {
                Some(mass_number) => write!(formatter, "[{mass_number}{}]", term.element.symbol())?,
                None => formatter.write_str(term.element.symbol())?,
            }
            if term.count != 1 {
                write!(formatter, "{}", term.count)?;
            }
        }
        match self.formal_charge.cmp(&0) {
            Ordering::Greater => {
                formatter.write_str("+")?;
                if self.formal_charge != 1 {
                    write!(formatter, "{}", self.formal_charge)?;
                }
            }
            Ordering::Less => {
                formatter.write_str("-")?;
                let magnitude = self.formal_charge.unsigned_abs();
                if magnitude != 1 {
                    write!(formatter, "{magnitude}")?;
                }
            }
            Ordering::Equal => {}
        }
        Ok(())
    }
}

/// A structured failure from molecular formula or mass calculation.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum MolecularDescriptorError {
    /// Installed implicit-hydrogen state is absent for a live atom that permits it.
    MissingImplicitHydrogens { atom: AtomId },
    /// No AME 2020 atomic mass is available for an asserted isotope label.
    UnknownIsotope {
        atom: AtomId,
        element: Element,
        mass_number: u16,
    },
    /// CIAAW does not define an abridged standard atomic weight for an element.
    MissingStandardAtomicWeight { atom: AtomId, element: Element },
    /// CIAAW does not define a naturally occurring isotope for an element.
    MissingNaturalIsotope { atom: AtomId, element: Element },
    /// Formula count accumulation exceeded `u64`.
    CountOverflow {
        element: Element,
        mass_number: Option<u16>,
    },
    /// The calculated mass exceeded the finite range of `f64`.
    MassOverflow,
}

impl fmt::Display for MolecularDescriptorError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingImplicitHydrogens { atom } => write!(
                formatter,
                "implicit-hydrogen state is not installed for atom {atom}"
            ),
            Self::UnknownIsotope {
                atom,
                element,
                mass_number,
            } => write!(
                formatter,
                "AME 2020 has no atomic mass for [{mass_number}{element}] at atom {atom}"
            ),
            Self::MissingStandardAtomicWeight { atom, element } => write!(
                formatter,
                "CIAAW 2024 has no abridged standard atomic weight for {element} at atom {atom}"
            ),
            Self::MissingNaturalIsotope { atom, element } => write!(
                formatter,
                "CIAAW 2024 has no naturally occurring isotope for {element} at atom {atom}"
            ),
            Self::CountOverflow {
                element,
                mass_number,
            } => match mass_number {
                Some(mass_number) => write!(
                    formatter,
                    "formula count overflowed for [{mass_number}{element}]"
                ),
                None => write!(formatter, "formula count overflowed for {element}"),
            },
            Self::MassOverflow => formatter.write_str("molecular mass exceeded finite f64 range"),
        }
    }
}

impl Error for MolecularDescriptorError {}

/// Constructs a structured molecular formula under an explicit hydrogen policy.
pub fn molecular_formula(
    molecule: &SmallMolecule,
    hydrogen_policy: HydrogenCountPolicy,
) -> Result<MolecularFormula, MolecularDescriptorError> {
    let mut counts = BTreeMap::new();
    visit_constituents(molecule, hydrogen_policy, |_, element, isotope, count| {
        add_count(&mut counts, element, isotope, count)
    })?;

    let contains_carbon = counts.keys().any(|(element, _)| element.symbol() == "C");
    let mut terms = counts
        .into_iter()
        .map(|((element, isotope), count)| FormulaTerm {
            element,
            isotope,
            count,
        })
        .collect::<Vec<_>>();
    terms.sort_by(|left, right| hill_term_cmp(*left, *right, contains_carbon));

    Ok(MolecularFormula {
        terms,
        formal_charge: molecule.graph().formal_charge(),
    })
}

/// Calculates average molecular mass in daltons under an explicit hydrogen policy.
pub fn average_mass(
    molecule: &SmallMolecule,
    hydrogen_policy: HydrogenCountPolicy,
) -> Result<Quantity<f64>, MolecularDescriptorError> {
    molecular_mass(molecule, hydrogen_policy, MassKind::Average)
}

/// Calculates monoisotopic molecular mass in daltons under an explicit hydrogen policy.
pub fn monoisotopic_mass(
    molecule: &SmallMolecule,
    hydrogen_policy: HydrogenCountPolicy,
) -> Result<Quantity<f64>, MolecularDescriptorError> {
    molecular_mass(molecule, hydrogen_policy, MassKind::Monoisotopic)
}

fn add_count(
    counts: &mut BTreeMap<(Element, Option<u16>), u64>,
    element: Element,
    isotope: Option<u16>,
    count: u64,
) -> Result<(), MolecularDescriptorError> {
    let entry = counts.entry((element, isotope)).or_default();
    *entry = entry
        .checked_add(count)
        .ok_or(MolecularDescriptorError::CountOverflow {
            element,
            mass_number: isotope,
        })?;
    Ok(())
}

fn hill_term_cmp(left: FormulaTerm, right: FormulaTerm, contains_carbon: bool) -> Ordering {
    hill_element_key(left.element, contains_carbon)
        .cmp(&hill_element_key(right.element, contains_carbon))
        .then_with(|| left.isotope.cmp(&right.isotope))
}

fn hill_element_key(element: Element, contains_carbon: bool) -> (u8, &'static str) {
    let symbol = element.symbol();
    if contains_carbon {
        match symbol {
            "C" => (0, ""),
            "H" => (1, ""),
            _ => (2, symbol),
        }
    } else {
        (0, symbol)
    }
}

fn visit_constituents(
    molecule: &SmallMolecule,
    hydrogen_policy: HydrogenCountPolicy,
    mut visit: impl FnMut(AtomId, Element, Option<u16>, u64) -> Result<(), MolecularDescriptorError>,
) -> Result<(), MolecularDescriptorError> {
    let graph = molecule.graph();
    let hydrogen = Element::from_atomic_number(1).expect("hydrogen is a supported element");
    for (atom_id, atom) in graph.atoms() {
        visit(atom_id, atom.element, atom.isotope, 1)?;
        if atom.explicit_hydrogens != 0 {
            visit(atom_id, hydrogen, None, u64::from(atom.explicit_hydrogens))?;
        }
        if hydrogen_policy == HydrogenCountPolicy::IncludePerceived && !atom.no_implicit_hydrogens {
            let implicit = graph
                .implicit_hydrogens(atom_id)
                .expect("a live atom identifier remains valid during read-only traversal")
                .ok_or(MolecularDescriptorError::MissingImplicitHydrogens { atom: atom_id })?;
            if implicit != 0 {
                visit(atom_id, hydrogen, None, u64::from(implicit))?;
            }
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
enum MassKind {
    Average,
    Monoisotopic,
}

fn molecular_mass(
    molecule: &SmallMolecule,
    hydrogen_policy: HydrogenCountPolicy,
    kind: MassKind,
) -> Result<Quantity<f64>, MolecularDescriptorError> {
    let mut mass = 0.0;
    visit_constituents(
        molecule,
        hydrogen_policy,
        |atom, element, isotope, count| {
            let constituent_mass = match (kind, isotope) {
                (_, Some(mass_number)) => exact_isotope_mass(element, mass_number).ok_or(
                    MolecularDescriptorError::UnknownIsotope {
                        atom,
                        element,
                        mass_number,
                    },
                )?,
                (MassKind::Average, None) => standard_atomic_weight(element).ok_or(
                    MolecularDescriptorError::MissingStandardAtomicWeight { atom, element },
                )?,
                (MassKind::Monoisotopic, None) => {
                    let mass_number = most_abundant_isotope(element)
                        .ok_or(MolecularDescriptorError::MissingNaturalIsotope { atom, element })?;
                    exact_isotope_mass(element, mass_number).ok_or(
                        MolecularDescriptorError::UnknownIsotope {
                            atom,
                            element,
                            mass_number,
                        },
                    )?
                }
            };
            mass += constituent_mass * count as f64;
            Ok(())
        },
    )?;
    mass -= molecule.graph().formal_charge() as f64 * data::ELECTRON_MASS_DA;
    if !mass.is_finite() {
        return Err(MolecularDescriptorError::MassOverflow);
    }
    Ok(Quantity::new(mass, DALTON))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Atom, AtomRadical, Molecule};
    use crate::perception;

    fn element(symbol: &str) -> Element {
        Element::from_symbol(symbol).expect("test element")
    }

    fn mass_value(value: Quantity<f64>) -> f64 {
        *value.value()
    }

    #[test]
    fn stored_and_perceived_hydrogens_are_explicit_policies() {
        let mut molecule = SmallMolecule::from_smiles("C").expect("methane parses");
        assert_eq!(
            molecular_formula(&molecule, HydrogenCountPolicy::StoredOnly)
                .expect("stored formula")
                .to_string(),
            "C"
        );
        assert_eq!(
            molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived),
            Err(MolecularDescriptorError::MissingImplicitHydrogens {
                atom: AtomId::new(0)
            })
        );

        perception::sanitize(&mut molecule).expect("methane sanitizes");
        assert_eq!(
            molecular_formula(&molecule, HydrogenCountPolicy::IncludePerceived)
                .expect("perceived formula")
                .to_string(),
            "CH4"
        );
    }

    #[test]
    fn hill_order_isotope_order_and_accessors_are_stable() {
        let mut graph = Molecule::new();
        let mut carbon_13 = Atom::new(element("C"));
        carbon_13.isotope = Some(13);
        graph.add_atom(carbon_13.clone());
        graph.add_atom(Atom::new(element("O")));
        graph.add_atom(Atom::new(element("C")));
        let mut carbon_12 = carbon_13;
        carbon_12.isotope = Some(12);
        graph.add_atom(carbon_12);
        let formula = molecular_formula(
            &SmallMolecule::from_graph(graph),
            HydrogenCountPolicy::StoredOnly,
        )
        .expect("formula");

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
    fn no_carbon_formula_and_disconnected_charge_use_one_entity() {
        let salt = SmallMolecule::from_smiles("[NH4+].[Cl-]").expect("salt parses");
        let formula = molecular_formula(&salt, HydrogenCountPolicy::StoredOnly).expect("formula");
        assert_eq!(formula.to_string(), "ClH4N");
        assert_eq!(formula.formal_charge(), 0);

        let ammonium = SmallMolecule::from_smiles("[NH4+]").expect("ammonium parses");
        assert_eq!(
            molecular_formula(&ammonium, HydrogenCountPolicy::StoredOnly)
                .expect("formula")
                .to_string(),
            "H4N+"
        );

        let chloride = SmallMolecule::from_smiles("[Cl-]").expect("chloride parses");
        assert_eq!(
            molecular_formula(&chloride, HydrogenCountPolicy::StoredOnly)
                .expect("formula")
                .to_string(),
            "Cl-"
        );

        let iron = SmallMolecule::from_smiles("[Fe+2]").expect("iron ion parses");
        assert_eq!(
            molecular_formula(&iron, HydrogenCountPolicy::StoredOnly)
                .expect("formula")
                .to_string(),
            "Fe+2"
        );
    }

    #[test]
    fn average_and_monoisotopic_mass_use_pinned_tables() {
        let mut water = SmallMolecule::from_smiles("O").expect("water parses");
        perception::sanitize(&mut water).expect("water sanitizes");
        let average = mass_value(
            average_mass(&water, HydrogenCountPolicy::IncludePerceived).expect("average mass"),
        );
        let monoisotopic = mass_value(
            monoisotopic_mass(&water, HydrogenCountPolicy::IncludePerceived)
                .expect("monoisotopic mass"),
        );
        assert!((average - 18.015).abs() < 1.0e-12);
        assert!((monoisotopic - 18.010_564_683_056).abs() < 1.0e-12);

        let carbon_13 = SmallMolecule::from_smiles("[13C]").expect("isotope parses");
        assert!(
            (mass_value(
                monoisotopic_mass(&carbon_13, HydrogenCountPolicy::StoredOnly)
                    .expect("isotope mass")
            ) - 13.003_354_835_34)
                .abs()
                < 1.0e-12
        );
    }

    #[test]
    fn formal_charge_applies_electron_mass_correction() {
        let sodium = SmallMolecule::from_smiles("[Na+]").expect("sodium parses");
        let chloride = SmallMolecule::from_smiles("[Cl-]").expect("chloride parses");
        let sodium_mass = mass_value(
            average_mass(&sodium, HydrogenCountPolicy::StoredOnly).expect("sodium mass"),
        );
        let chloride_mass = mass_value(
            average_mass(&chloride, HydrogenCountPolicy::StoredOnly).expect("chloride mass"),
        );
        assert!((sodium_mass - (22.990 - data::ELECTRON_MASS_DA)).abs() < 1.0e-12);
        assert!((chloride_mass - (35.45 + data::ELECTRON_MASS_DA)).abs() < 1.0e-12);
    }

    #[test]
    fn radicals_need_no_mass_correction_beyond_charge() {
        let mut graph = Molecule::new();
        let mut carbon = Atom::new(element("C"));
        carbon.radical = Some(AtomRadical::Doublet);
        graph.add_atom(carbon);
        let molecule = SmallMolecule::from_graph(graph);
        assert_eq!(
            molecular_formula(&molecule, HydrogenCountPolicy::StoredOnly)
                .expect("formula")
                .to_string(),
            "C"
        );
        assert_eq!(
            mass_value(
                monoisotopic_mass(&molecule, HydrogenCountPolicy::StoredOnly)
                    .expect("radical mass")
            ),
            12.0
        );
    }

    #[test]
    fn unavailable_weights_and_isotopes_are_structured_errors() {
        let technetium = SmallMolecule::from_smiles("[Tc]").expect("technetium parses");
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
            SmallMolecule::from_smiles("[99Tc]").expect("technetium isotope parses");
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

        let mut graph = Molecule::new();
        let mut impossible = Atom::new(element("C"));
        impossible.isotope = Some(999);
        graph.add_atom(impossible);
        assert!(matches!(
            monoisotopic_mass(
                &SmallMolecule::from_graph(graph),
                HydrogenCountPolicy::StoredOnly
            ),
            Err(MolecularDescriptorError::UnknownIsotope {
                mass_number: 999,
                ..
            })
        ));
    }

    #[test]
    fn empty_molecule_has_empty_formula_and_zero_masses() {
        let molecule = SmallMolecule::new();
        let formula =
            molecular_formula(&molecule, HydrogenCountPolicy::StoredOnly).expect("formula");
        assert!(formula.is_empty());
        assert_eq!(formula.to_string(), "");
        assert_eq!(
            mass_value(average_mass(&molecule, HydrogenCountPolicy::StoredOnly).expect("mass")),
            0.0
        );
        assert_eq!(
            mass_value(
                monoisotopic_mass(&molecule, HydrogenCountPolicy::StoredOnly).expect("mass")
            ),
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
}
