use super::*;
use crate::io::ComponentLibrary;

#[test]
fn table_is_sorted_and_unique() {
    let entries = entries();
    assert!(entries.len() > 2_000);
    assert!(entries.windows(2).all(|pair| pair[0].id < pair[1].id));
}

#[test]
fn modified_residues_carry_their_parent_codes() {
    let threonine = dictionary_code("tpo").expect("phosphothreonine");
    assert_eq!(threonine.linking, ComponentLinking::Peptide);
    assert_eq!(threonine.code, Some("T"));
    assert_eq!(threonine.parents().collect::<Vec<_>>(), ["THR"]);
    let pseudouridine = dictionary_code("PSU").expect("pseudouridine");
    assert_eq!(pseudouridine.linking, ComponentLinking::Nucleotide);
    assert_eq!(pseudouridine.code, Some("U"));
    assert_eq!(dictionary_code("HOH"), None);
}

#[test]
fn bundled_polymer_templates_agree_with_the_table() {
    for component in ComponentLibrary::bundled().components() {
        let linking = component.linking();
        if !matches!(
            linking,
            ComponentLinking::Peptide | ComponentLinking::Nucleotide
        ) {
            continue;
        }
        let entry = dictionary_code(component.id()).expect("bundled polymer component");
        assert_eq!(entry.linking, linking, "{}", component.id());
        assert_eq!(
            entry.code,
            component.one_letter_code(),
            "{}",
            component.id()
        );
        assert!(
            entry
                .parents()
                .eq(component.parent_ids().iter().map(String::as_str)),
            "{}",
            component.id()
        );
    }
}
