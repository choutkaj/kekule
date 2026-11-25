//! The valence states a builder makes: per element and formal charge, the
//! valences an atom takes, from which its hydrogen count follows.
use crate::molecule::Element;
use std::fmt;

/// The formal charges of each builder element, ascending, by atomic number.
const CHARGES: &[(u8, &[i8])] = &[
    (1, &[0]),
    (5, &[-1, 0]),
    (6, &[-1, 0, 1]),
    (7, &[-1, 0, 1]),
    (8, &[-1, 0, 1]),
    (9, &[-1, 0]),
    (14, &[0]),
    (15, &[0, 1]),
    (16, &[-1, 0, 1, 2]),
    (17, &[-1, 0]),
    (34, &[-1, 0, 1, 2]),
    (35, &[-1, 0]),
    (53, &[-1, 0]),
];

/// The valences of each builder state, ascending: the octet's for second-row
/// elements, the expanded ones of phosphorus, sulfur, and selenium in their
/// neutral states, four for the sulfur and selenium dications of a
/// charge-separated sulfone, and none for a halide.
const VALENCES: &[(u8, i8, &[u8])] = &[
    (1, 0, &[1]),
    (5, -1, &[4]),
    (5, 0, &[3]),
    (6, -1, &[3]),
    (6, 0, &[4]),
    (6, 1, &[3]),
    (7, -1, &[2]),
    (7, 0, &[3]),
    (7, 1, &[4]),
    (8, -1, &[1]),
    (8, 0, &[2]),
    (8, 1, &[3]),
    (9, -1, &[0]),
    (9, 0, &[1]),
    (14, 0, &[4]),
    (15, 0, &[3, 5]),
    (15, 1, &[4]),
    (16, -1, &[1]),
    (16, 0, &[2, 4, 6]),
    (16, 1, &[3]),
    (16, 2, &[4]),
    (17, -1, &[0]),
    (17, 0, &[1]),
    (34, -1, &[1]),
    (34, 0, &[2, 4, 6]),
    (34, 1, &[3]),
    (34, 2, &[4]),
    (35, -1, &[0]),
    (35, 0, &[1]),
    (53, -1, &[0]),
    (53, 0, &[1]),
];

/// The formal charges a builder gives atoms of `element`, ascending, or
/// `None` for an element outside the builder's table: H 0; B −1, 0; C, N,
/// and O −1, 0, +1; S and Se −1, 0, +1, +2; F, Cl, Br, and I −1, 0; Si 0;
/// P 0, +1. A host steps an
/// atom's charge through these and names them when it refuses another.
///
/// ```
/// use kekule::{molecule::Element, structure::placement::builder_charges};
/// assert_eq!(builder_charges(Element::O), Some(&[-1, 0, 1][..]));
/// assert_eq!(builder_charges(Element::from_symbol("Fe").unwrap()), None);
/// ```
pub fn builder_charges(element: Element) -> Option<&'static [i8]> {
    CHARGES
        .iter()
        .find(|(number, _)| *number == element.atomic_number())
        .map(|(_, charges)| *charges)
}

/// The valences, ascending, of an atom of `element` with `formal_charge` in
/// a builder, or `None` when that is not a builder state: C 4, C⁺ and C⁻ 3,
/// N 3, N⁺ 4, O 2, O⁺ 3, B 3, B⁻ 4, P 3 and 5, S 2, 4, and 6, S²⁺ 4, a
/// halide 0.
///
/// ```
/// use kekule::{molecule::Element, structure::placement::builder_valences};
/// assert_eq!(builder_valences(Element::S, 0), Some(&[2, 4, 6][..]));
/// assert_eq!(builder_valences(Element::N, 1), Some(&[4][..]));
/// assert_eq!(builder_valences(Element::O, 2), None);
/// ```
pub fn builder_valences(element: Element, formal_charge: i8) -> Option<&'static [u8]> {
    VALENCES
        .iter()
        .find(|(number, charge, _)| *number == element.atomic_number() && *charge == formal_charge)
        .map(|(_, _, valences)| *valences)
}

/// The hydrogens an atom of `element` with `formal_charge` takes in a
/// builder when `valence` is the bond orders to its other neighbors plus its
/// radical electrons: its lowest [builder valence](builder_valences) at or
/// above `valence`, less `valence`. The atom's own hydrogens do not count,
/// so a sulfur with two carbons takes none at charge 0 and one at +1.
///
/// ```
/// use kekule::{molecule::Element, structure::placement::builder_hydrogen_count};
/// assert_eq!(builder_hydrogen_count(Element::C, 0, 1)?, 3);
/// assert_eq!(builder_hydrogen_count(Element::S, 0, 3)?, 1);
/// assert_eq!(builder_hydrogen_count(Element::N, 1, 0)?, 4);
/// assert!(builder_hydrogen_count(Element::C, 0, 5).is_err());
/// # Ok::<(), kekule::structure::placement::BuilderValenceError>(())
/// ```
///
/// # Errors
///
/// [`BuilderValenceError::UnknownElement`] outside the builder's table,
/// [`BuilderValenceError::UnsupportedCharge`] for a charge it does not
/// give the element, and [`BuilderValenceError::ExcessValence`] when
/// `valence` passes every valence of the state.
pub fn builder_hydrogen_count(
    element: Element,
    formal_charge: i8,
    valence: usize,
) -> Result<usize, BuilderValenceError> {
    let charges =
        builder_charges(element).ok_or(BuilderValenceError::UnknownElement { element })?;
    let valences =
        builder_valences(element, formal_charge).ok_or(BuilderValenceError::UnsupportedCharge {
            element,
            formal_charge,
            charges,
        })?;
    valences
        .iter()
        .map(|allowed| usize::from(*allowed))
        .find(|allowed| *allowed >= valence)
        .map(|allowed| allowed - valence)
        .ok_or(BuilderValenceError::ExcessValence {
            element,
            formal_charge,
            valence,
            valences,
        })
}

/// The largest builder valence of `element` with `formal_charge`, or `None`
/// outside the builder's states.
pub(crate) fn builder_max_valence(element: Element, formal_charge: i8) -> Option<usize> {
    builder_valences(element, formal_charge)
        .and_then(|valences| valences.last())
        .map(|valence| usize::from(*valence))
}

/// An atom's element, charge, or valence is outside the builder's states
/// ([`builder_hydrogen_count`]).
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum BuilderValenceError {
    /// The element is not in the builder's table ([`builder_charges`]).
    UnknownElement {
        /// The element.
        element: Element,
    },
    /// The builder does not give the element this charge.
    UnsupportedCharge {
        /// The element.
        element: Element,
        /// The charge asked for.
        formal_charge: i8,
        /// The charges the builder gives the element, ascending.
        charges: &'static [i8],
    },
    /// The valence passes every valence of the state; the atom has too many
    /// bonds for its element and charge.
    ExcessValence {
        /// The element.
        element: Element,
        /// The charge.
        formal_charge: i8,
        /// The bond orders to the atom's other neighbors plus its radical
        /// electrons.
        valence: usize,
        /// The valences of the state, ascending.
        valences: &'static [u8],
    },
}

impl BuilderValenceError {
    /// Writes the error as sentences about `subject`, a capitalized noun
    /// phrase for the atom, such as `Edit atom 3`.
    pub(crate) fn write_about(
        &self,
        formatter: &mut fmt::Formatter<'_>,
        subject: impl fmt::Display,
    ) -> fmt::Result {
        match self {
            Self::UnknownElement { element } => write!(
                formatter,
                "{subject} is {}. {} has no builder valences.",
                element.symbol(),
                element.symbol()
            ),
            Self::UnsupportedCharge {
                element,
                formal_charge,
                charges,
            } => write!(
                formatter,
                "{subject} is {} with charge {}. {} takes charge {}.",
                element.symbol(),
                Signed(*formal_charge),
                element.symbol(),
                alternatives(&charges.iter().map(|c| Signed(*c)).collect::<Vec<_>>())
            ),
            Self::ExcessValence {
                element,
                formal_charge,
                valence,
                valences,
            } => {
                write!(
                    formatter,
                    "{subject} has valence {valence}. {}",
                    element.symbol()
                )?;
                if *formal_charge != 0 {
                    write!(formatter, " with charge {}", Signed(*formal_charge))?;
                }
                write!(formatter, " takes valence {}.", alternatives(valences))
            }
        }
    }
}

impl fmt::Display for BuilderValenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_about(formatter, "The atom")
    }
}

impl std::error::Error for BuilderValenceError {}

/// A charge with its sign, `0` without one.
struct Signed(i8);

impl fmt::Display for Signed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            0 => formatter.write_str("0"),
            charge => write!(formatter, "{charge:+}"),
        }
    }
}

/// Items joined as alternatives, such as `-1, 0 or +1`.
fn alternatives<T: fmt::Display>(items: &[T]) -> String {
    match items {
        [] => String::new(),
        [only] => only.to_string(),
        [rest @ .., last] => {
            let rest = rest.iter().map(ToString::to_string).collect::<Vec<_>>();
            format!("{} or {last}", rest.join(", "))
        }
    }
}

#[cfg(test)]
mod tests;
