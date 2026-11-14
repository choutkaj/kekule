//! Read-only small-molecule descriptors: formula and masses, atomic masses,
//! canonical atom ranking, and rotatable bonds.
//!
//! Descriptor calculation never mutates or perceives the input; a descriptor
//! that needs perception it cannot reuse computes it on a temporary copy.
//! Callers select explicitly whether installed inferred-hydrogen state is part
//! of the calculation.
//!
//! Average masses use CIAAW 2024 abridged standard atomic weights for unlabeled
//! atoms. Monoisotopic masses use AME 2020 atomic masses for the most abundant
//! naturally occurring isotope selected from CIAAW 2024 compositions. Explicit
//! isotope labels use their AME 2020 mass in either calculation. Both mass
//! calculations subtract the aggregate formal charge times the CODATA 2022
//! electron mass. Missing standard weights or natural isotopes are reported as
//! structured errors; no mass-number substitute is used.
//!
//! # Example
//!
//! ```
//! use kekule::descriptors::{molecular_formula, HydrogenCountPolicy};
//!
//! let mut ethanol = kekule::smiles::to_molecules("CCO")?.remove(0);
//! // Stored counts omit the hydrogens SMILES leaves to inference.
//! let stored = molecular_formula(&ethanol, HydrogenCountPolicy::StoredOnly)?;
//! assert_eq!(stored.to_string(), "C2O");
//!
//! ethanol.perceive()?;
//! let formula = molecular_formula(&ethanol, HydrogenCountPolicy::IncludePerceived)?;
//! assert_eq!(formula.to_string(), "C2H6O");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod data;

pub use crate::algorithms::{
    canonical_atom_ranking, rotatable_bonds, CanonicalAtomRanking, RotatableBondOptions,
    RotatableBondSet,
};

/// Neutral atomic mass used for inertia and vibrational analysis. Unlabeled
/// elements use standard atomic weights; labeled atoms use exact isotope masses.
/// No implicit hydrogens or formal-charge electron corrections are included.
///
/// Returns `None` when the element has no standard atomic weight (unlabeled)
/// or the isotope has no AME 2020 mass (labeled).
pub fn atomic_mass(element: Element, isotope: Option<u16>) -> Option<Quantity<f64>> {
    let mass = match isotope {
        Some(number) => data::exact_isotope_mass(element, number),
        None => data::standard_atomic_weight(element),
    }?;
    Some(Quantity::new(mass, DALTON))
}

/// The mass number a legacy CTfile mass difference counts from: the most
/// abundant natural isotope, as RDKit reads it. `None` for an element without
/// a naturally abundant isotope.
pub(crate) fn mass_difference_reference(element: Element) -> Option<u16> {
    data::most_abundant_isotope(element)
}

use crate::diagnostics::UnknownHydrogenCount;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

use crate::molecule::{AtomId, Element, Molecule};
use crate::units::{Quantity, DALTON};

use data::{exact_isotope_mass, most_abundant_isotope, standard_atomic_weight};

/// Selects which non-atom hydrogen counts contribute to a descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum HydrogenCountPolicy {
    /// Count explicit graph hydrogens and fixed implicit counts only.
    /// Inferred counts are omitted, so this can describe a partial composition.
    StoredOnly,
    /// Count explicit hydrogens and the complete implicit hydrogen count.
    /// Requires perception for atoms whose declarations enable inference.
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
    UnknownHydrogenCount(UnknownHydrogenCount),
    /// No AME 2020 atomic mass is available for an asserted isotope label.
    UnknownIsotope {
        /// The labelled atom.
        atom: AtomId,
        /// The atom's element.
        element: Element,
        /// The asserted mass number.
        mass_number: u16,
    },
    /// CIAAW does not define an abridged standard atomic weight for an element.
    MissingStandardAtomicWeight {
        /// The atom of that element.
        atom: AtomId,
        /// The element without a standard atomic weight.
        element: Element,
    },
    /// CIAAW does not define a naturally occurring isotope for an element.
    MissingNaturalIsotope {
        /// The atom of that element.
        atom: AtomId,
        /// The element without a naturally occurring isotope.
        element: Element,
    },
    /// Formula count accumulation exceeded `u64`.
    CountOverflow {
        /// The element whose count overflowed.
        element: Element,
        /// The isotope label counted, or `None` for unlabelled atoms.
        mass_number: Option<u16>,
    },
    /// The calculated mass exceeded the finite range of `f64`.
    MassOverflow,
}

impl fmt::Display for MolecularDescriptorError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownHydrogenCount(error) => error.fmt(formatter),
            Self::UnknownIsotope {
                atom,
                element,
                mass_number,
            } => write!(
                formatter,
                "AME 2020 has no atomic mass for {mass_number}{element}, the isotope of {atom}."
            ),
            Self::MissingStandardAtomicWeight { atom, element } => write!(
                formatter,
                "CIAAW 2024 has no standard atomic weight for {element}, the element of {atom}."
            ),
            Self::MissingNaturalIsotope { atom, element } => write!(
                formatter,
                "CIAAW 2024 lists no naturally occurring isotope of {element}, the element of {atom}."
            ),
            Self::CountOverflow {
                element,
                mass_number,
            } => match mass_number {
                Some(mass_number) => write!(
                    formatter,
                    "The formula count of {mass_number}{element} is too large."
                ),
                None => write!(formatter, "The formula count of {element} is too large."),
            },
            Self::MassOverflow => formatter.write_str("The molecular mass is too large."),
        }
    }
}

impl Error for MolecularDescriptorError {}

/// Constructs a structured molecular formula under an explicit hydrogen policy.
///
/// # Errors
///
/// [`MolecularDescriptorError::UnknownHydrogenCount`] when
/// [`HydrogenCountPolicy::IncludePerceived`] meets an atom without installed
/// implicit-hydrogen state, and [`MolecularDescriptorError::CountOverflow`]
/// when a count exceeds `u64`.
pub fn molecular_formula(
    molecule: &Molecule,
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
        formal_charge: molecule.formal_charge(),
    })
}

/// Calculates average molecular mass in daltons under an explicit hydrogen policy.
///
/// # Errors
///
/// [`MolecularDescriptorError::UnknownHydrogenCount`] when
/// [`HydrogenCountPolicy::IncludePerceived`] meets an atom without installed
/// implicit-hydrogen state, [`MolecularDescriptorError::UnknownIsotope`] for an
/// isotope label without an AME 2020 mass,
/// [`MolecularDescriptorError::MissingStandardAtomicWeight`] for an
/// unlabeled element without a standard atomic weight, and
/// [`MolecularDescriptorError::MassOverflow`] when the mass is not finite.
pub fn average_mass(
    molecule: &Molecule,
    hydrogen_policy: HydrogenCountPolicy,
) -> Result<Quantity<f64>, MolecularDescriptorError> {
    molecular_mass(molecule, hydrogen_policy, MassKind::Average)
}

/// Calculates monoisotopic molecular mass in daltons under an explicit hydrogen policy.
///
/// # Errors
///
/// [`MolecularDescriptorError::UnknownHydrogenCount`] when
/// [`HydrogenCountPolicy::IncludePerceived`] meets an atom without installed
/// implicit-hydrogen state, [`MolecularDescriptorError::UnknownIsotope`] for an
/// isotope label without an AME 2020 mass,
/// [`MolecularDescriptorError::MissingNaturalIsotope`] for an
/// unlabeled element without a naturally occurring isotope, and
/// [`MolecularDescriptorError::MassOverflow`] when the mass is not finite.
pub fn monoisotopic_mass(
    molecule: &Molecule,
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
    molecule: &Molecule,
    hydrogen_policy: HydrogenCountPolicy,
    mut visit: impl FnMut(AtomId, Element, Option<u16>, u64) -> Result<(), MolecularDescriptorError>,
) -> Result<(), MolecularDescriptorError> {
    let graph = molecule;
    let hydrogen = Element::from_atomic_number(1).expect("hydrogen is a supported element");
    for (atom_id, atom) in graph.atoms() {
        visit(
            atom_id,
            atom.element,
            atom.isotope.map(std::num::NonZero::get),
            1,
        )?;
        let implicit = match hydrogen_policy {
            HydrogenCountPolicy::StoredOnly => usize::from(atom.hydrogens.represented_count()),
            HydrogenCountPolicy::IncludePerceived => graph.implicit_hydrogens(atom_id).ok_or(
                MolecularDescriptorError::UnknownHydrogenCount(UnknownHydrogenCount::new(atom_id)),
            )?,
        };
        if implicit != 0 {
            // At most two u8 contributions, widened by the public count API.
            visit(atom_id, hydrogen, None, implicit as u64)?;
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
    molecule: &Molecule,
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
    mass -= molecule.formal_charge() as f64 * data::ELECTRON_MASS_DA;
    if !mass.is_finite() {
        return Err(MolecularDescriptorError::MassOverflow);
    }
    Ok(Quantity::new(mass, DALTON))
}

#[cfg(test)]
mod tests;
