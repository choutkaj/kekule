use super::numbered;

#[test]
fn only_symbol_and_digits_are_numbered() {
    assert_eq!(numbered("C21", "C"), Some(21));
    assert_eq!(numbered("c7", "C"), Some(7));
    assert_eq!(numbered("CL1", "C"), None);
    assert_eq!(numbered("CA", "C"), None);
    assert_eq!(numbered("C", "C"), None);
    assert_eq!(numbered("CL12", "CL"), Some(12));
    assert_eq!(numbered("H", "HE"), None);
}

#[test]
fn fragments_placed_alone_take_their_locants_and_distinct_names() {
    use super::super::{fragment, fragments, FragmentKind};
    use super::fragment_names;
    let names = |name: &str| fragment_names(fragment(name).unwrap()).join(" ");
    assert_eq!(names("benzene"), "C1 C2 C3 C4 C5 C6 H1 H2 H3 H4 H5 H6");
    assert_eq!(names("pyridine"), "N1 C2 C3 C4 C5 C6 H2 H3 H4 H5 H6");
    assert_eq!(
        names("naphthalene"),
        "C1 C2 C3 C4 C4A C5 C6 C7 C8 C8A H1 H2 H3 H4 H5 H6 H7 H8"
    );
    assert_eq!(
        names("indole"),
        "N1 C2 C3 C3A C4 C5 C6 C7 C7A H1 H2 H3 H4 H5 H6 H7"
    );
    assert_eq!(
        names("piperidine"),
        "N1 C2 C3 C4 C5 C6 H1 H21 H22 H31 H32 H41 H42 H51 H52 H61 H62"
    );
    assert_eq!(names("methylamino"), "N1 C1 HN11 HN12 H11 H12 H13");
    for fragment in fragments() {
        let molecule = fragment.molecule();
        let placed = fragment_names(fragment);
        let distinct = placed.iter().collect::<std::collections::BTreeSet<_>>();
        assert_eq!(distinct.len(), placed.len(), "{}", fragment.name());
        assert!(placed.iter().all(|name| name.len() <= 4), "{placed:?}");
        for (atom, name) in molecule.atom_ids().zip(&placed) {
            let element = molecule.atom(atom).unwrap().element();
            match fragment.locant(atom) {
                Some(locant) => assert_eq!(
                    *name,
                    format!("{}{}", element.symbol(), locant).to_ascii_uppercase()
                ),
                // Every heavy atom of a ring fragment is a ring atom.
                None => assert!(
                    fragment.kind() == FragmentKind::Group || element.atomic_number() == 1,
                    "{} {name}",
                    fragment.name()
                ),
            }
        }
    }
}

#[test]
fn new_atoms_never_take_a_name_the_template_places_elsewhere() {
    use super::ResidueNames;
    use crate::components::ComponentLibrary;
    use crate::molecule::Element;
    let library = ComponentLibrary::bundled();
    let templates = Some((
        library.get("LEU").unwrap(),
        library.variants("LEU").collect(),
    ));
    let mut names = ResidueNames::new(templates, ["N", "CA", "C", "O", "CB"].map(str::to_owned));
    let hydrogen = Element::from_atomic_number(1).unwrap();
    let mut allocate = |parent| names.allocate(hydrogen, Some(parent));
    // The template's own names, then the N terminus's conventional one.
    assert_eq!(
        ["N", "N", "N", "CA", "CB", "CB"].map(&mut allocate),
        ["H", "H2", "H1", "HA", "HB2", "HB3"]
    );
    // Elsewhere, numbers that no template name or terminal convention uses.
    assert_eq!(allocate("CB"), "H4");
    assert_eq!(allocate("N"), "H3");
}
