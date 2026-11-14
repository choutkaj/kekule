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
