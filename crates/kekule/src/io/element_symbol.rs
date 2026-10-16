//! Element symbols as macromolecular formats write them.
//!
//! The Chemical Component Dictionary, mmCIF `_atom_site.type_symbol`, and the
//! PDB element column write symbols in any case and name the hydrogen isotopes
//! deuterium and tritium `D` and `T`. One reader and one writer keep the three
//! formats in agreement.

use crate::molecule::Element;

const HYDROGEN: u8 = 1;

/// Reads a macromolecular element symbol in any case. `D` and `T` read as
/// hydrogen with mass number 2 and 3.
pub(crate) fn parse_structure_element_symbol(symbol: &str) -> Option<(Element, Option<u16>)> {
    let hydrogen = || Element::from_atomic_number(HYDROGEN);
    match symbol {
        "D" | "d" => Some((hydrogen()?, Some(2))),
        "T" | "t" => Some((hydrogen()?, Some(3))),
        _ => {
            let mut chars = symbol.chars();
            let first = chars.next()?;
            let mut canonical = first.to_ascii_uppercase().to_string();
            canonical.extend(chars.map(|char| char.to_ascii_lowercase()));
            Element::from_symbol(&canonical).map(|element| (element, None))
        }
    }
}

/// The symbol a macromolecular format writes for an atom, or `None` when the
/// format cannot name its isotope: only deuterium and tritium have symbols.
pub(crate) fn structure_element_symbol(
    element: Element,
    isotope: Option<u16>,
) -> Option<&'static str> {
    match (element.atomic_number(), isotope) {
        (_, None) => Some(element.symbol()),
        (HYDROGEN, Some(2)) => Some("D"),
        (HYDROGEN, Some(3)) => Some("T"),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn element(symbol: &str) -> Element {
        Element::from_symbol(symbol).expect("element")
    }

    #[test]
    fn reads_symbols_in_any_case_and_hydrogen_isotopes() {
        let cases = [
            ("C", Some((element("C"), None))),
            ("FE", Some((element("Fe"), None))),
            ("zn", Some((element("Zn"), None))),
            ("H", Some((element("H"), None))),
            ("D", Some((element("H"), Some(2)))),
            ("d", Some((element("H"), Some(2)))),
            ("T", Some((element("H"), Some(3)))),
            ("Ds", Some((element("Ds"), None))),
            ("TI", Some((element("Ti"), None))),
            ("X", None),
            ("", None),
        ];
        for (symbol, expected) in cases {
            assert_eq!(parse_structure_element_symbol(symbol), expected, "{symbol}");
        }
    }

    #[test]
    fn writes_deuterium_and_tritium_and_no_other_isotope() {
        let cases = [
            (element("C"), None, Some("C")),
            (element("H"), None, Some("H")),
            (element("H"), Some(1), None),
            (element("H"), Some(2), Some("D")),
            (element("H"), Some(3), Some("T")),
            (element("C"), Some(13), None),
        ];
        for (element, isotope, expected) in cases {
            assert_eq!(structure_element_symbol(element, isotope), expected);
        }
    }
}
