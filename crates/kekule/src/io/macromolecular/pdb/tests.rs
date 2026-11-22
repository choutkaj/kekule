use super::*;

fn element(symbol: &str) -> NameElement {
    NameElement::Element(element_symbol(symbol).unwrap())
}

#[test]
fn atom_name_alignment_implies_the_element() {
    for (name, residue, symbol) in [
        (" CA ", "ALA", "C"),
        ("CA  ", "CA", "Ca"),
        ("FE  ", "HEM", "Fe"),
        (" N  ", "ALA", "N"),
        ("1HB ", "ALA", "H"),
        ("HG21", "THR", "H"),
        ("HG  ", "HG", "Hg"),
        (" HG ", "CYS", "H"),
        ("ZN  ", "ZN", "Zn"),
        ("CL1 ", "LIG", "Cl"),
        ("OXT ", "ALA", "O"),
        ("SE  ", "MSE", "Se"),
    ] {
        assert_eq!(
            element_from_atom_name(name, residue, true),
            element(symbol),
            "{name:?} in {residue}"
        );
    }
    assert_eq!(
        element_from_atom_name(" 1  ", "LIG", true),
        NameElement::None
    );
}

#[test]
fn names_that_read_as_two_elements_are_ambiguous() {
    let ambiguous = |one: &str, two: &str| {
        NameElement::Ambiguous([element_symbol(one).unwrap(), element_symbol(two).unwrap()])
    };
    // Without the alignment convention, column 13 may start either symbol.
    assert_eq!(
        element_from_atom_name("CL1 ", "LIG", false),
        ambiguous("C", "Cl")
    );
    assert_eq!(
        element_from_atom_name("FE  ", "HEM", false),
        ambiguous("F", "Fe")
    );
    // Four-character names fill the columns either way.
    assert_eq!(
        element_from_atom_name("CA12", "LIG", true),
        ambiguous("C", "Ca")
    );
    // Writers that pad every name into column 14 leave ` CL1` either.
    assert_eq!(
        element_from_atom_name(" CL1", "LIG", false),
        ambiguous("C", "Cl")
    );
    assert_eq!(element_from_atom_name(" C1 ", "LIG", false), element("C"));
    // A name led by `H` is a hydrogen, as four-character names are.
    for name in [" HO1", " HG1", " HE2", "HG1 "] {
        assert_eq!(element_from_atom_name(name, "LIG", false), element("H"));
    }
    assert_eq!(element_from_atom_name("1HG1", "LIG", false), element("H"));
}
