use super::*;

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
            element_from_atom_name(name, residue).map(Element::symbol),
            Some(symbol),
            "{name:?} in {residue}"
        );
    }
    assert_eq!(element_from_atom_name(" 1  ", "LIG"), None);
}
